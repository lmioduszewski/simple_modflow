"""GeoPackage-first construction of reusable MODFLOW model inputs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np

from myflopy.advanced import (
    chd_spec,
    drn_spec,
    evt_spec,
    ghb_spec,
    rch_spec,
    riv_spec,
    wel_spec,
)
from myflopy.specs import ModelContext, PackageSpec

SurfaceReference = str


def _period_value(row, field: str | list[str] | tuple[str, ...], period: int):
    """Read a constant or per-period value from one GeoPackage feature."""

    if isinstance(field, str):
        return row[field]
    if not field:
        raise ValueError("At least one value field is required.")
    return row[field[min(period, len(field) - 1)]]


@dataclass(frozen=True, slots=True)
class CellSurfaceOffset:
    """Resolve a value from a cell surface plus a scalar or row-field offset.

    ``offset`` and ``minimum`` may be numeric constants or GeoPackage field
    names. A positive offset raises the surface-derived value; a negative offset
    lowers it.
    """

    reference: SurfaceReference
    offset: int | float | str = 0.0
    minimum: int | float | str | None = None

    def __post_init__(self) -> None:
        """Lowercase and validate ``reference`` (one of model_top/cell_top/cell_bottom)."""

        reference = self.reference.lower()
        valid = {"model_top", "cell_top", "cell_bottom"}
        if reference not in valid:
            raise ValueError(f"reference must be one of: {', '.join(sorted(valid))}.")
        object.__setattr__(self, "reference", reference)

    @staticmethod
    def _row_or_value(row, value: int | float | str | None) -> float | None:
        """Resolve a numeric constant or a GeoPackage field name to a float (``None`` -> ``None``)."""

        if value is None:
            return None
        if isinstance(value, str):
            return float(row[value])
        return float(value)

    def required_fields(self) -> set[str]:
        """Return GeoPackage fields required to resolve this value."""

        fields = set()
        if isinstance(self.offset, str):
            fields.add(self.offset)
        if isinstance(self.minimum, str):
            fields.add(self.minimum)
        return fields

    def resolve(self, source: GeoPackageSource, row, *, layer: int, cell: int) -> float:
        """Return the resolved elevation/head for one mapped cell."""

        surface = source.surface_value(self.reference, layer=layer, cell=cell)
        value = surface + self._row_or_value(row, self.offset)
        minimum = self._row_or_value(row, self.minimum)
        if minimum is not None:
            value = max(value, minimum)
        return float(value)


RowValue = str | int | float | list[str] | tuple[str, ...] | CellSurfaceOffset


def _required_value_fields(value: RowValue) -> set[str]:
    """Return GeoPackage fields needed by one value specification."""

    if isinstance(value, CellSurfaceOffset):
        return value.required_fields()
    if isinstance(value, str):
        return {value}
    return set(value) if isinstance(value, (list, tuple)) else set()


def _metadata_value(value: RowValue):
    """Return a manifest-friendly representation of one value specification."""

    if isinstance(value, CellSurfaceOffset):
        return {
            "type": "CellSurfaceOffset",
            "reference": value.reference,
            "offset": value.offset,
            "minimum": value.minimum,
        }
    return value


@dataclass(slots=True)
class GeoPackageSource:
    """Map one GeoPackage layer onto a model grid and build reusable inputs.

    Layer values are assumed to be one-based by default because that is the
    usual convention in external GIS inputs. Stress-period values are assumed
    to be zero-based.
    """

    path: Path | str
    context: ModelContext
    nper: int
    layer: str | None = None
    name_field: str | None = "name"
    layer_field: str | None = "layer"
    period_field: str | None = None
    layer_base: int = 1
    period_base: int = 0
    _gdf: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        """Coerce ``path`` and validate the file exists, a grid is present, and ``nper >= 1``."""

        self.path = Path(self.path)
        if not self.path.exists():
            raise FileNotFoundError(f"GeoPackage not found: {self.path}")
        if self.context.grid is None:
            raise ValueError("ModelContext.grid is required for GeoPackage mapping.")
        if self.nper < 1:
            raise ValueError("nper must be at least 1.")

    @property
    def grid(self):
        """Grid helper carried by the model context."""

        return self.context.grid

    @property
    def gdf(self) -> gpd.GeoDataFrame:
        """Read and align the GeoPackage layer with the model grid CRS."""

        if self._gdf is None:
            kwargs = {} if self.layer is None else {"layer": self.layer}
            gdf = gpd.read_file(self.path, **kwargs)
            grid_crs = getattr(self.grid, "crs", None)
            if grid_crs is not None and gdf.crs is not None and gdf.crs != grid_crs:
                gdf = gdf.to_crs(grid_crs)
            self._gdf = gdf
        return self._gdf

    def _active(self, layer: int, cell: int) -> bool:
        """Whether ``(layer, cell)`` is active in the context's idomain (True if none set)."""

        domain = self.context.domain
        if domain is None:
            return True
        values = np.asarray(domain)
        if values.ndim == 1:
            return bool(values[cell])
        if layer >= values.shape[0]:
            raise ValueError(
                f"GeoPackage layer {layer} is outside domain with {values.shape[0]} layers."
            )
        return bool(values[layer, cell])

    def _layer(self, row) -> int:
        """The zero-based grid layer for a feature (from ``layer_field`` minus ``layer_base``)."""

        if self.layer_field is None:
            return 0
        layer = int(row[self.layer_field]) - self.layer_base
        if layer < 0:
            raise ValueError(f"GeoPackage layer resolves to negative index: {layer}")
        return layer

    def _periods(self, row) -> range | tuple[int]:
        """The stress periods a feature applies to: all periods, or its single ``period_field``."""

        if self.period_field is None:
            return range(self.nper)
        period = int(row[self.period_field]) - self.period_base
        if period not in range(self.nper):
            raise ValueError(f"GeoPackage period is outside the simulation: {period}")
        return (period,)

    def _cells(self, geometry, *, layer: int, edges_only: bool = False) -> list[int]:
        """The active grid cells a feature's ``geometry`` intersects (optionally edge cells only)."""

        cells = self.grid.gdf_vorPolys.index[
            self.grid.gdf_vorPolys.intersects(geometry)
        ].tolist()
        if edges_only:
            edge_cells = set(self.grid.get_grid_edge())
            cells = [cell for cell in cells if cell in edge_cells]
        return [int(cell) for cell in cells if self._active(layer, int(cell))]

    def _surfaces(self):
        """The surface source for :class:`CellSurfaceOffset`: the context's, else the grid's top/botm."""

        surfaces = self.context.surfaces
        if surfaces is None:
            surfaces = getattr(self.grid, "gdf_topbtm", None)
        if surfaces is None:
            raise ValueError(
                "ModelContext.surfaces or grid.gdf_topbtm is required for CellSurfaceOffset."
            )
        return surfaces

    @staticmethod
    def _cell_sequence(value: Any, cell: int) -> float:
        """One cell's value from a scalar (broadcast) or a per-cell array."""

        array = np.asarray(value, dtype=float)
        if array.ndim == 0:
            return float(array)
        return float(array.reshape(-1)[cell])

    def surface_value(self, reference: SurfaceReference, *, layer: int, cell: int) -> float:
        """Return a model/cell surface value for a mapped boundary cell."""

        surfaces = self._surfaces()
        if isinstance(surfaces, dict):
            if reference in {"model_top", "cell_top"} and layer == 0:
                return self._cell_sequence(surfaces["top"], cell)
            if reference == "model_top":
                return self._cell_sequence(surfaces["top"], cell)
            bottom = surfaces.get("bottom", surfaces.get("botm"))
            if bottom is None:
                raise ValueError("Surface dictionary requires 'bottom' or 'botm'.")
            bottom_array = np.asarray(bottom, dtype=float)
            if reference == "cell_top":
                if bottom_array.ndim == 1:
                    raise ValueError("cell_top for layers below 0 requires multilayer bottom surfaces.")
                return self._cell_sequence(bottom_array[layer - 1], cell)
            if bottom_array.ndim == 1:
                if layer != 0:
                    raise ValueError("cell_bottom for layers below 0 requires multilayer bottom surfaces.")
                return self._cell_sequence(bottom_array, cell)
            return self._cell_sequence(bottom_array[layer], cell)

        columns = getattr(surfaces, "columns", ())
        if reference == "model_top":
            candidates = (0, "top", "model_top")
        elif reference == "cell_top":
            candidates = ((0, "top", "model_top") if layer == 0 else (layer, f"layer_{layer}_top"))
        else:
            candidates = (layer + 1, "bottom" if layer == 0 else f"layer_{layer}_bottom")
        for column in candidates:
            if column in columns:
                return float(surfaces.loc[cell, column])
        raise ValueError(
            f"Could not resolve {reference!r} for layer {layer}; available surface columns are {list(columns)!r}."
        )

    def _value(self, row, value: RowValue, *, period: int, layer: int, cell: int):
        """Resolve one value spec for a cell: surface offset, numeric constant, or (per-period) field."""

        if isinstance(value, CellSurfaceOffset):
            return value.resolve(self, row, layer=layer, cell=cell)
        if isinstance(value, (int, float)):
            return value
        return _period_value(row, value, period)

    def _boundary_data(
        self,
        *fields: RowValue,
        edges_only: bool = False,
        boundnames: bool = False,
    ) -> dict[int, list[list[Any]]]:
        """Build MF6 stress-period data by mapping every feature to its cells and value fields.

        Validates the required columns exist, then for each feature emits one
        ``[(layer, cell), *values(, boundname)]`` record per applicable period and
        intersected cell, keyed by stress period.
        """

        required = {
            field_name
            for field in fields
            for field_name in _required_value_fields(field)
        }
        if self.layer_field is not None:
            required.add(self.layer_field)
        if self.period_field is not None:
            required.add(self.period_field)
        if boundnames and self.name_field is not None:
            required.add(self.name_field)
        missing = sorted(required - set(self.gdf.columns))
        if missing:
            raise ValueError(
                f"GeoPackage '{self.path.name}' is missing fields: {', '.join(missing)}"
            )

        data = {period: [] for period in range(self.nper)}
        for index, row in self.gdf.iterrows():
            layer = self._layer(row)
            cells = self._cells(row.geometry, layer=layer, edges_only=edges_only)
            name = (
                str(row[self.name_field])
                if boundnames and self.name_field is not None
                else str(index)
            )
            for period in self._periods(row):
                for cell in cells:
                    values = [
                        self._value(row, field, period=period, layer=layer, cell=cell)
                        for field in fields
                    ]
                    record = [(layer, cell), *values]
                    if boundnames:
                        record.append(name)
                    data[period].append(record)
        return data

    def _metadata(self, method: str, **fields: Any) -> dict[str, Any]:
        """Return serializable provenance for a generated model input."""

        return {
            "source_type": "geopackage",
            "source_path": str(self.path),
            "source_layer": self.layer,
            "builder": method,
            "fields": {name: _metadata_value(value) for name, value in fields.items()},
            "name_field": self.name_field,
            "layer_field": self.layer_field,
            "period_field": self.period_field,
            "layer_base": self.layer_base,
            "period_base": self.period_base,
        }

    def chd(
        self,
        *,
        head: RowValue = "head",
        name: str = "chd",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a CHD package spec from GeoPackage features."""

        return chd_spec(
            self._boundary_data(head, edges_only=edges_only, boundnames=boundnames),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(**self._metadata("chd", head=head))

    def ghb(
        self,
        *,
        head: RowValue = "head",
        conductance: RowValue = "conductance",
        name: str = "ghb",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a GHB package spec from GeoPackage features."""

        return ghb_spec(
            self._boundary_data(
                head,
                conductance,
                edges_only=edges_only,
                boundnames=boundnames,
            ),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(**self._metadata("ghb", head=head, conductance=conductance))

    def drn(
        self,
        *,
        elevation: RowValue = "elevation",
        conductance: RowValue = "conductance",
        name: str = "drn",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a DRN package spec from GeoPackage features."""

        return drn_spec(
            self._boundary_data(
                elevation,
                conductance,
                edges_only=edges_only,
                boundnames=boundnames,
            ),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(
            **self._metadata("drn", elevation=elevation, conductance=conductance)
        )

    def riv(
        self,
        *,
        stage: RowValue = "stage",
        conductance: RowValue = "conductance",
        rbot: RowValue = "rbot",
        name: str = "riv",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a RIV package spec from GeoPackage features."""

        return riv_spec(
            self._boundary_data(
                stage,
                conductance,
                rbot,
                edges_only=edges_only,
                boundnames=boundnames,
            ),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(
            **self._metadata("riv", stage=stage, conductance=conductance, rbot=rbot)
        )

    def wel(
        self,
        *,
        rate: RowValue = "rate",
        name: str = "wel",
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a WEL package spec from GeoPackage features."""

        return wel_spec(
            self._boundary_data(rate, boundnames=boundnames),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(**self._metadata("wel", rate=rate))

    def rch(
        self,
        *,
        recharge: RowValue = "recharge",
        name: str = "rch",
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a list-based RCH package spec from GeoPackage features."""

        return rch_spec(
            self._boundary_data(recharge, boundnames=boundnames),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(**self._metadata("rch", recharge=recharge))

    def evt(
        self,
        *,
        surface: RowValue = "surface",
        rate: RowValue = "rate",
        depth: RowValue = "depth",
        name: str = "evt",
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a list-based EVT package spec from GeoPackage features.

        Single-segment ET only: feature mapping emits fixed
        ``(cellid, surface, rate, depth)`` records, which is exactly the
        ``nseg=1`` record shape. Segmented ET needs ``nseg - 1`` extra
        ``pxdp``/``petm`` values per record that no feature mapping supplies,
        so it is rejected here rather than failing later inside FloPy.
        """

        if int(options.get("nseg", 1)) != 1:
            raise ValueError(
                "GeoPackage-driven EVT supports nseg=1 only: mapped features "
                "yield (cellid, surface, rate, depth) records with no pxdp/petm "
                "values. For segmented ET, assemble the records yourself and "
                "use mf.evt(...) / mf.evt.flopy(...) with nseg=."
            )

        return evt_spec(
            self._boundary_data(surface, rate, depth, boundnames=boundnames),
            name=name,
            boundnames=boundnames,
            **options,
        ).with_metadata(
            **self._metadata("evt", surface=surface, rate=rate, depth=depth)
        )

    def k_array(
        self,
        *,
        value: str = "k",
        nlay: int,
        defaults: float | list[float],
    ) -> np.ndarray:
        """Return an ``(nlay, ncpl)`` K array from GeoPackage polygons."""

        fallback = [defaults] * nlay if np.isscalar(defaults) else list(defaults)
        if len(fallback) != nlay:
            raise ValueError("defaults must be a scalar or contain one value per layer.")
        result = np.array(
            [[float(fallback[layer])] * self.grid.ncpl for layer in range(nlay)],
            dtype=float,
        )
        for _, row in self.gdf.iterrows():
            layer = self._layer(row)
            if layer >= nlay:
                raise ValueError(f"GeoPackage layer {layer} is outside nlay={nlay}.")
            for cell in self._cells(row.geometry, layer=layer):
                result[layer, cell] = float(row[value])
        return result


__all__ = ["CellSurfaceOffset", "GeoPackageSource"]
