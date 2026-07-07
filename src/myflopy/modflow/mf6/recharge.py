"""Recharge helpers that translate GIS/tabular sources into MF6 recharge inputs."""

from __future__ import annotations
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Real
from typing import TYPE_CHECKING, Any

import flopy.utils.binaryfile
import geopandas as gpd
import numpy as np
import pandas as pd
import pickle
import shapely as shp
from pathlib import Path

from myflopy.advanced import rch_spec
from myflopy.modflow.mf6.boundary_support import (
    build_cell_id,
    filter_inactive_cells,
    merge_stress_period_data,
    normalize_grid_type,
)
from myflopy.modflow.mf6.boundaries import Boundaries
from myflopy.modflow.utils.prism_ppt import PrismPrecipScaling
from myflopy.specs import ModelContext, PackageSpec

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

idxx = pd.IndexSlice
inches_to_feet = 1 / 12
CellId = tuple[int, int]


@dataclass(frozen=True, slots=True)
class RCHBuilder:
    """Prepare a list-based RCH package from explicit recharge configuration.

    ``recharge`` accepts the same practical shapes used in model setup:

    - a scalar or time-series name applied to every selected cell and period
    - a sequence with one value per selected cell, repeated every period
    - a mapping of ``cellid -> value``, repeated every period
    - a mapping of ``period -> scalar | sequence | cell mapping``

    Cell ids are DISV-style ``(layer, cell)`` tuples. Integer cell ids are also
    accepted and are placed on ``layer``. Use tuple keys when a mapping might
    otherwise look like period keys, such as ``{0: ...}``.
    """

    context: ModelContext
    nper: int
    recharge: Any
    cells: str | Sequence[int | CellId] = "top_active"
    layer: int = 0
    name_by_cell: Mapping[int | CellId, str] | None = None
    boundnames: bool = False
    name: str = "rch"
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate ``nper >= 1`` and that a string ``cells`` selector is a known keyword."""

        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        if isinstance(self.cells, str) and self.cells not in {"top_active", "all_active", "surface_only"}:
            raise ValueError("cells must be 'top_active', 'surface_only', 'all_active', or explicit cells.")

    def with_updates(self, **updates: Any) -> "RCHBuilder":
        """Return a changed builder without modifying the original."""

        return replace(self, **updates)

    @property
    def domain(self):
        """Return the active-domain array carried by the model context."""

        return None if self.context.domain is None else np.asarray(self.context.domain)

    @property
    def selected_cells(self) -> tuple[CellId, ...]:
        """Return recharge target cells as DISV-style ``(layer, cell)`` IDs."""

        if not isinstance(self.cells, str):
            return tuple(self._normalize_cell(cell) for cell in self.cells)

        domain = self.domain
        if domain is None:
            grid = self.context.grid
            if grid is None or not hasattr(grid, "ncpl"):
                raise ValueError("context.domain or context.grid.ncpl is required for generated RCH cells.")
            if self.cells == "all_active":
                return tuple((self.layer, cell) for cell in range(int(grid.ncpl)))
            return tuple((self.layer, cell) for cell in range(int(grid.ncpl)))

        if domain.ndim == 1:
            active = np.flatnonzero(domain > 0)
            return tuple((self.layer, int(cell)) for cell in active)

        selected: list[CellId] = []
        for cell in range(domain.shape[1]):
            active_layers = np.flatnonzero(domain[:, cell] > 0)
            if len(active_layers) == 0:
                continue
            if self.cells == "all_active":
                selected.extend((int(layer), int(cell)) for layer in active_layers)
            else:
                selected.append((int(active_layers[0]), int(cell)))
        return tuple(selected)

    def _normalize_cell(self, cell: int | CellId) -> CellId:
        """Coerce a cell spec to a ``(layer, cell)`` pair (a bare int lands on ``self.layer``)."""

        if isinstance(cell, Real):
            return (int(self.layer), int(cell))
        if not isinstance(cell, Sequence) or isinstance(cell, str) or len(cell) != 2:
            raise ValueError("Explicit RCH cells must be integer cells or (layer, cell) pairs.")
        layer, node = cell
        return (int(layer), int(node))

    def _cell_lookup_keys(self, cellid: CellId) -> tuple[Any, ...]:
        """The keys to try when looking up a cell in a user mapping: the pair, then the bare cell."""

        return (cellid, cellid[1])

    def _name_for_cell(self, cellid: CellId) -> str:
        """The boundname for a cell -- from ``name_by_cell`` if present, else ``<name>_<layer>_<cell>``."""

        if not self.name_by_cell:
            return f"{self.name}_{cellid[0]}_{cellid[1]}"
        for key in self._cell_lookup_keys(cellid):
            if key in self.name_by_cell:
                return str(self.name_by_cell[key])
        return f"{self.name}_{cellid[0]}_{cellid[1]}"

    @staticmethod
    def _period_keys(value: Mapping[Any, Any], nper: int) -> bool:
        """Whether a mapping is keyed by stress-period indices (all int keys in ``range(nper)``)."""

        return bool(value) and all(isinstance(key, int) and key in range(nper) for key in value)

    def _value_for_cell(self, value: Any, cellid: CellId, selected_cells: tuple[CellId, ...]) -> Any:
        """Resolve one cell's recharge from a cell mapping, a scalar/field, or a per-cell sequence."""

        if isinstance(value, Mapping):
            for key in self._cell_lookup_keys(cellid):
                if key in value:
                    return value[key]
            raise ValueError(f"recharge is missing cell {cellid}.")
        if isinstance(value, Real) or isinstance(value, str):
            return value
        values = list(value)
        if len(values) != len(selected_cells):
            raise ValueError("recharge sequence must contain one value per selected cell.")
        return values[selected_cells.index(cellid)]

    def _value_for_period_cell(self, period: int, cellid: CellId, selected_cells: tuple[CellId, ...]) -> Any:
        """Resolve one cell's recharge for a given period (unwrapping a per-period mapping first)."""

        value = self.recharge
        if isinstance(value, Mapping) and self._period_keys(value, self.nper):
            if period not in value:
                raise ValueError(f"recharge is missing period {period}.")
            value = value[period]
        return self._value_for_cell(value, cellid, selected_cells)

    @property
    def stress_period_data(self) -> dict[int, list[list[Any]]]:
        """Return normalized MF6 RCH stress-period data."""

        cells = self.selected_cells
        data: dict[int, list[list[Any]]] = {period: [] for period in range(self.nper)}
        for period in range(self.nper):
            for cellid in cells:
                value = self._value_for_period_cell(period, cellid, cells)
                if isinstance(value, Real) and float(value) < 0.0:
                    raise ValueError("recharge must be nonnegative.")
                row = [cellid, value]
                if self.boundnames:
                    row.append(self._name_for_cell(cellid))
                data[period].append(row)
        return data

    def validate(self) -> None:
        """Raise clear errors for invalid recharge inputs."""

        cells = self.selected_cells
        if not cells:
            raise ValueError("RCHBuilder selected no active recharge cells.")
        domain = self.domain
        if domain is not None and domain.ndim == 1:
            for _layer, cell in cells:
                if cell >= domain.shape[0] or domain[cell] <= 0:
                    raise ValueError(f"RCH cell {cell} is inactive or outside the domain.")
        if domain is not None and domain.ndim > 1:
            for layer, cell in cells:
                if layer >= domain.shape[0] or cell >= domain.shape[1] or domain[layer, cell] <= 0:
                    raise ValueError(f"RCH cell {(layer, cell)} is inactive or outside the domain.")
        self.stress_period_data

    def build(self) -> PackageSpec:
        """Validate stored configuration and return its package spec."""

        self.validate()
        return rch_spec(
            self.stress_period_data,
            name=self.name,
            boundnames=self.boundnames,
            **dict(self.options),
        ).with_metadata(
            builder="RCHBuilder",
            cells=list(self.selected_cells),
            recharge_form=type(self.recharge).__name__,
        )


class RechargeFromShp(Boundaries):
    """Build recharge stress-period data from polygon GIS inputs."""

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            shp_gpkg: Path = None,
            uid: str = None,
            rch_fields: list | slice = None,
            rch_fields_to_pers: list = None,
            xlsx_rch: Path = None,
            background_rch: int | float = 0.0,
            apply_background_rch: bool = True,
            rch_in_vol: bool = False,
            multiplier=1,
            grid_type: str = 'disv',
            limit_to_k33: bool = True,
            limit_to_k33_by: int | float = 1,
            verbose: bool = False,
            **kwargs

    ):
        """Parameters
        ----------
        model
            Model to which the recharge applies.
        vor
            Grid helper used to intersect polygons with model cells.
        shp_gpkg
            Polygon geometry source for recharge zones.
        uid
            Unique-id field in the geometry attributes.
        rch_fields, rch_fields_to_pers
            Attribute fields containing recharge values and their per-period mapping.
        xlsx_rch
            Optional spreadsheet overriding recharge attribute values.
        background_rch
            Background recharge value used when requested.
        apply_background_rch
            Whether unspecified periods should receive ``background_rch``.
        rch_in_vol
            Whether source values are volumes that should be converted to rates.
        multiplier
            Scalar applied to the source recharge values.
        grid_type
            MODFLOW grid type, usually ``disv`` or ``disu``.
        limit_to_k33, limit_to_k33_by
            Optional vertical-conductivity-based cap on recharge values.
        verbose
            Whether to emit extra progress/details while building data.
        kwargs
            Additional arguments passed through to :class:`Boundaries`.
        """
        super().__init__(model, vor, shp_gpkg, uid, **kwargs)
        self.bound_type = 'rch'
        self.fields = rch_fields
        self.rch_fields_to_pers = rch_fields_to_pers
        self.background_rch = background_rch
        self.apply_background_rch = apply_background_rch
        self.xlsx_rch = xlsx_rch
        self._cell_ids = None
        self._rch_fields = None
        self.uid = uid
        self._fields_to_pers = None
        self._recharges = None
        self.rch_in_vol = rch_in_vol
        self.multiplier = multiplier
        self.limit_to_k33 = limit_to_k33
        self.limit_to_k33_by = limit_to_k33_by
        self.grid_type = normalize_grid_type(grid_type)
        self.verbose = verbose

    @property
    def cell_ids(self):
        """gets cell ids for each recharge polygon in the shapefile as a dict. Keys are
        unique ids (uids)"""
        if self._cell_ids is None:
            cell_ids = self.intersections_no_duplicates
            cell_ids = cell_ids['no_dup'].to_dict()
            self._cell_ids = cell_ids
        return self._cell_ids

    @property
    def rch_fields(self):
        """gets a DataFrame of just the recharge data fields from the shapefile or from excel file if proivded.
        Used to build the recharge dict for input into a flopy modflow model"""
        if self._rch_fields is None:
            if self.xlsx_rch:
                """use excel if it exists, otherwise get from shapefile"""
                rch_fields = pd.read_excel(self.xlsx_rch).set_index(self.uid)
                assert len(rch_fields) == len(
                    self.gdf), 'number of rows in excel file and number of shapefile polys must be the equal'
            else:
                rch_fields = self.gdf.loc[:, self.fields]
            if self.rch_in_vol:  # if recharge is in volumes, divide by the voronoi area of each polygon
                rch_fields = rch_fields.apply(lambda row: row / self.gdf.area)
            rch_fields = rch_fields * self.multiplier
            self._rch_fields = rch_fields
        return self._rch_fields

    @property
    def fields_to_pers(self):
        """creates a list of indices and values that correspond to the columns/fields of the rch_fields DataFrame.
        In the case that the length of the fields is not long enough, -1 or -2 is added with correspond to either,
        use the last index given for the remaining stress periods or set the remaining stress periods to a background
        recharge, defined by setting the background recharge class attribute and setting apply_background_rch to True."""
        if self._fields_to_pers is None:
            fields_to_pers = self.rch_fields_to_pers.copy()
            fields_to_pers = [] if fields_to_pers is None else fields_to_pers
            for per in range(self.nper):
                if len(fields_to_pers) <= per:
                    if self.apply_background_rch:
                        fields_to_pers.append(-2)
                    else:
                        fields_to_pers.append(-1)
            self._fields_to_pers = fields_to_pers
        return self._fields_to_pers

    @property
    def recharges(self):
        """gets a dict of recharge values for each recharge area (each uid) for each stress period. Pass to
        get_rch() to generate a recharge dict to pass to the flopy recharge class."""
        nper = self.nper
        uids = self.gdf.index.to_list()
        scaled_rch_fields = self.rch_fields.mul(self.shp_to_vor_poly_scale, axis=0)
        rch_fields_dict = scaled_rch_fields.to_dict()
        rch_fields_cols = self.rch_fields.columns

        if self._recharges is None:
            recharges = {}
            for uid in uids:
                recharges[uid] = []
            for per in range(nper):
                if self.fields_to_pers[per] == -1:
                    for uid in uids:
                        #  if -1 then apply the last field indicated in the provided rch_fields_to_pers agrument
                        field_for_per = rch_fields_cols[self.rch_fields_to_pers[-1]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
                elif self.fields_to_pers[per] == -2:
                    for uid in uids:
                        recharges[uid].append(self.background_rch)
                else:
                    for uid in uids:
                        field_for_per = rch_fields_cols[self.fields_to_pers[per]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
            self._recharges = recharges
        return self._recharges

    def from_polygons(
            self,
            cell_ids: dict = None,
            recharges: dict = None,
            background_rch: int | float = None,
            register_regions: bool = False,
            region_name_prefix: str | None = None,
            combined_region_name: str | None = None,
            region_tags: list[str] | None = None,
            overwrite_regions: bool = False,
    ) -> dict:
        """Build recharge stress-period data from polygon features on the builder."""
        rch_dict = {}
        grid_type = normalize_grid_type(self.grid_type)
        cell_ids = self.cell_ids if cell_ids is None else cell_ids
        recharges = self.recharges if recharges is None else recharges
        background_rch = self.background_rch if background_rch is None else background_rch
        k33 = self.model.gwf.npf.k33.data[0] if self.limit_to_k33 else None
        nper = self.nper
        assert nper == len(list(recharges.values())[0]), 'Number of periods and length of recharge values must match'
        region_cellids_by_name: dict[str, list] = {}
        region_layers_by_name: dict[str, int] = {}
        region_geometries_by_name: dict[str, object] = {}
        region_metadata_by_name: dict[str, dict] = {}

        for per in range(nper):
            cell_list = []
            all_rch_cells = set()

            for name, cell_nums in cell_ids.items():
                active_cells = filter_inactive_cells(cell_nums, self.inactive_cells)
                recharge = recharges[name][per]
                if register_regions:
                    region_layers_by_name[name] = 0
                    if self.gdf is not None and name in self.gdf.index:
                        region_geometries_by_name[name] = self.gdf.loc[name, "geometry"]
                    region_metadata_by_name.setdefault(name, {"background_rch": background_rch})
                    region_cellids_by_name.setdefault(name, [])
                for cell in active_cells:
                    all_rch_cells.add(cell)
                    cell_id = build_cell_id(cell, grid_type=grid_type, layer=0)
                    if self.limit_to_k33 and k33[cell] < recharge:
                        limited_recharge = k33[cell] * self.limit_to_k33_by
                        if self.verbose:
                            print(
                                f'cell {cell_id} has k33 {k33[cell]}, which is less than given recharge {recharge}.'
                                f' Changing recharge to {k33[cell] * self.limit_to_k33_by}'
                            )
                        cell_list.append([cell_id, limited_recharge])
                    else:
                        cell_list.append([cell_id, recharge])
                    if register_regions:
                        region_cellids_by_name[name].append(cell_id)

            if background_rch is not None:
                for cell in filter_inactive_cells(range(self.vor.ncpl), self.inactive_cells):
                    if cell not in all_rch_cells:
                        cell_id = build_cell_id(cell, grid_type=grid_type, layer=0)
                        cell_list.append([cell_id, background_rch])

            rch_dict[per] = cell_list

        if register_regions and region_cellids_by_name:
            self._register_boundary_groups(
                cellids_by_name=region_cellids_by_name,
                layers_by_name=region_layers_by_name,
                geometries_by_name=region_geometries_by_name,
                region_name_prefix=region_name_prefix or self.bound_type,
                combined_region_name=combined_region_name,
                tags=region_tags,
                metadata_by_name=region_metadata_by_name,
                overwrite=overwrite_regions,
            )

        return rch_dict

    def from_vector(self, **kwargs) -> dict:
        """Alias for :meth:`from_polygons` for shapefile/geopackage workflows."""

        return self.from_polygons(**kwargs)

    def get_rch(self, **kwargs) -> dict:
        """Backward-compatible alias for :meth:`from_polygons`."""

        return self.from_polygons(**kwargs)

    @staticmethod
    def add_to_rch_dict(rch_dict: dict, rch_to_add: dict, replace: bool = True) -> dict:
        """Merge ``rch_to_add`` into a stress-period RCH dict (``replace`` overwrites collisions)."""

        return merge_stress_period_data(rch_dict, rch_to_add, replace=replace)


class RechargeFromPrism(Boundaries):

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            prism_raster: Path = None,
            weather_station_location: Path = None,
            weather_station_precip: list | pd.Series = None,
            et_dict: dict = None,
            period_months: list = None,
    ):
        """Build recharge from PRISM-scaled station precipitation minus monthly ET.

        Scales a weather-station precip series to each Voronoi cell by a PRISM
        raster ratio, then subtracts per-month ET; validates the precip length,
        ET months/lengths, and period-months on assignment.
        """

        super().__init__(model, vor)
        self.bound_type = 'rch'
        self.prism_raster = prism_raster
        self.weather_station_location = weather_station_location

        self._prism_scaling = None
        self._weather_station_precip = None
        self._et_dict = None
        self._scaled_precip = None
        self._period_months = None

        self.weather_station_precip = weather_station_precip
        self.et_dict = et_dict
        self.period_months = period_months

    @property
    def prism_scaling(self):
        """returns a Pandas Series where the values are the scaling factors to
        use for precipitation for each voronoi cell, listed by voronoi cell number index"""
        if self._prism_scaling is None:
            self._prism_scaling = PrismPrecipScaling(
                self.vor,
                self.prism_raster,
                self.weather_station_location).scaling
        return self._prism_scaling

    @property
    def weather_station_precip(self):
        """The per-stress-period weather-station precipitation series."""

        return self._weather_station_precip

    @weather_station_precip.setter
    def weather_station_precip(self, value):
        """Set the station precip series, asserting one value per stress period."""

        assert len(value) == self.model.nper, \
            'length of weather station precip must equal number of stress periods'
        self._weather_station_precip = value

    @property
    def et_dict(self):
        """dictionary of evapotranspiration values for each month for each cell
        in the model. Keys are month integers between 1 and 12"""
        return self._et_dict

    @et_dict.setter
    def et_dict(self, value):
        """Set the month->per-cell ET dict, asserting month keys (1-12) and one value per cell."""

        assert isinstance(value, dict), 'et_dict must be a dict'
        assert all(month in range(1, 13) for month in value.keys()), \
            'keys of et_dict must be month integers between 1 and 12'
        assert all(len(ets) == self.model.modelgrid.ncpl for ets in value.values()), \
            'length of et_dict values must equal number of model cells'
        self._et_dict = value

    @property
    def scaled_precip(self):
        """scales the weather station precipitation for each voronoi cell
        prism scaling factor"""
        if self._scaled_precip is None:
            scaled_precip = []
            for per in range(self.model.nper):
                scaled_precip.append(self.prism_scaling * self.weather_station_precip[per])
            self._scaled_precip = scaled_precip
        return self._scaled_precip

    @property
    def period_months(self):
        """a list of the months (int), one int for each stress period"""
        return self._period_months

    @period_months.setter
    def period_months(self, value):
        """Set the per-period month list, asserting one month (1-12) per stress period."""

        assert len(value) == self.model.nper, \
            'length of period_months must equal number of stress periods'
        assert all(month in range(1, 13) for month in value), \
            'values in period_months must be integers between 1 and 12'
        self._period_months = value

    @property
    def rch_dict(self):
        """returns a dictionary of recharge values for each stress period"""
        rch_dict = {}
        for per in range(self.model.nper):
            rch_dict[per] = [
                self.scaled_precip[per] - self.et_dict[self.period_months[per]]
            ]
        return rch_dict


# Preferred alias for the vector-driven recharge builder API.
RCHFromVector = RechargeFromShp

__all__ = ["RCHBuilder", "RCHFromVector", "RechargeFromPrism", "RechargeFromShp"]

