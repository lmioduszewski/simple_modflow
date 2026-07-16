"""General-head-boundary helpers built from polygons, lines, and tabular inputs."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.boundaries import Boundaries
from myflopy.modflow.mf6.boundary_support import (
    build_cell_id,
    merge_stress_period_data,
    normalize_grid_type,
)

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

idxx = pd.IndexSlice
inches_to_feet = 1 / 12


class GHB(Boundaries):
    """Build MF6 GHB stress-period data from GIS features and per-period values."""

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            shp_gpkg: Path = None,
            line_shp: Path = None,
            uid: str = 'name',
            crs: int = 2927,
            fields: list | slice = None,
            fields_to_pers: list = None,
            xlsx: Path = None,
            idomain: list[int] | pd.Series = None,
            idomain_path: Path = None,
            verbose: bool = False,
    ):
        """Parameters
        ----------
        model
            Model to which the boundary applies.
        vor
            Grid helper used to intersect geometry with model cells.
        shp_gpkg
            Polygon geometry source for the GHB features.
        line_shp
            Optional line geometry source for line-based GHB segments.
        uid
            Unique-id field in the geometry attributes.
        crs
            EPSG code for the geometry inputs.
        fields, fields_to_pers
            Attribute fields and period mapping for time-varying values.
        xlsx
            Optional spreadsheet overriding geometry-attribute values.
        idomain, idomain_path
            Optional active-domain definition used to filter inactive cells.
        verbose
            Whether to emit extra progress/details while building data.
        """
        super().__init__(model, vor, shp_gpkg, uid, crs, idomain=idomain, idomain_path=idomain_path)
        self.bound_type = 'ghb'
        self._field_names = fields
        self.xlsx = xlsx
        self._cell_ids = None
        self._fields = None
        self.uid = uid
        self._fields_to_pers = fields_to_pers
        self._recharges = None
        self.line_shp = line_shp
        self._line_ghb = None
        self.line_fields = {'cond': ['cond_strt', 'cond_end'], 'elev': ['elev_strt', 'elev_end'], 'layer': 'layer'}
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
    def fields(self):
        """gets a DataFrame of just the data fields from the shapefile or from excel file if provided.
        Used to build the dict for input into a flopy modflow model"""
        if self._fields is None:
            if self.xlsx:
                """use excel if it exists, otherwise get from shapefile"""
                fields = pd.read_excel(self.xlsx).set_index(self.uid)
                assert len(fields) == len(
                    self.gdf), 'number of rows in excel file and number of shapefile polys must be the equal'
            else:
                fields = self.gdf.loc[:, self._field_names]
            self._fields = fields
        return self._fields

    @property
    def fields_to_pers(self):
        """creates a list of indices and values that correspond to the columns/fields of the fields DataFrame.
        In the case that the length of the fields is not long enough, -1 is added which tells the class to fill the
        remaining stress periods with the last field value given"""
        if self._fields_to_pers is None:
            self._fields_to_pers = []
        fields_to_pers = list(self._fields_to_pers)
        for per in range(self.nper):
            if len(fields_to_pers) <= per:
                fields_to_pers.append(-1)
        self._fields_to_pers = fields_to_pers
        return self._fields_to_pers

    @property
    def data(self):
        """gets a dict of values for each boundary area (each uid) for each stress period. Pass to
        self.get() to generate a dict to pass to the flopy modflow boundary class."""
        nper = self.nper
        uids = self.gdf.index.to_list()
        scaled_rch_fields = self.fields.mul(self.shp_to_vor_poly_scale, axis=0)
        rch_fields_dict = scaled_rch_fields.to_dict()
        rch_fields_cols = self.fields.columns

        if self._recharges is None:
            recharges = {}
            for uid in uids:
                recharges[uid] = []
            for per in range(nper):
                if self.fields_to_pers[per] == -1:
                    for uid in uids:
                        #  if -1 then apply the last field indicated in the provided fields_to_pers argument
                        field_for_per = rch_fields_cols[self.fields_to_pers[-1]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
                else:
                    for uid in uids:
                        field_for_per = rch_fields_cols[self.fields_to_pers[per]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
            self._recharges = recharges
        return self._recharges

    @property
    def line_ghb(self):
        """A line-driven GHB :class:`Boundaries` built from ``line_shp`` (cached)."""

        if self._line_ghb is None:
            assert isinstance(self.line_shp, Path), 'No line shapefile provided. Set the line_shp attribute'
            line_ghb = Boundaries(self.model, self.vor, self.line_shp, self.uid, self.crs, bound_type='ghb')
            self._line_ghb = line_ghb
        return self._line_ghb

    def add_line_ghb(self, existing_ghb_dict: dict = None):
        """adds a line based ghb boundary to an existing ghb dict"""
        elev_strt = self.line_ghb.gdf[self.line_fields['elev'][0]]
        elev_end = self.line_ghb.gdf[self.line_fields['elev'][1]]
        elev_delta = elev_end - elev_strt
        cond_strt = self.line_ghb.gdf[self.line_fields['cond'][0]]
        cond_end = self.line_ghb.gdf[self.line_fields['cond'][1]]
        cond_delta = cond_end - cond_strt
        ghb_line_lists = []
        for idx, _delta in enumerate(elev_delta):
            line_len = self.line_ghb.gdf.geometry[idx].length
            sorted_cells = self.line_ghb.sorted_cells_along_line(idx)
            line_layer = self.line_ghb.gdf['layer'].iloc[idx] - 1
            for cell in sorted_cells.index:
                len_ratio = sorted_cells.loc[cell, 'distance_along_line'] / line_len
                elev_delta_ratio = elev_delta.iloc[idx] * len_ratio
                this_elev = elev_strt.iloc[idx] + elev_delta_ratio
                cond_delta_ratio = cond_delta.iloc[idx] * len_ratio
                this_cond = cond_strt.iloc[idx] + cond_delta_ratio
                this_cell = (line_layer, cell)
                this_ghb_list = [this_cell, this_elev, this_cond]
                ghb_line_lists.append(this_ghb_list)
        # TODO this only works assumiung constant ghb boundary. Adds the ghb data to stress period 1. Need to make more generic
        to_add = {0: ghb_line_lists}
        updated_ghb_dict = self.add_to_dict(existing_ghb_dict, to_add)
        return updated_ghb_dict

    def from_polygons(
            self,
            grid_type: str = 'disv',
            fields: dict = None,
            elev_reference: dict = None,
            reference_offset: float | int = 0,
            edges_only: bool = False,
            register_regions: bool = False,
            region_name_prefix: str | None = None,
            combined_region_name: str | None = None,
            region_tags: list[str] | None = None,
            overwrite_regions: bool = False,
    ) -> dict:
        """Build GHB stress-period data from polygon features on the builder."""
        elev_reference = {} if elev_reference is None else elev_reference
        grid_type = normalize_grid_type(grid_type)
        nper = self.nper if self.nper is not None else 1
        if fields is None:
            fields = {
                'name': 'name',
                'elevation': 'elev',
                'height_over_btm': 'height',
                'conductance': 'cond',
                'layer': 'layer',
                'min_elev': 'min_elev',
            }
        # get bottoms of model layers from voronoi grid
        vor = self._require_vor()
        lyr_botms = vor.gdf_topbtm.drop('geometry', axis='columns').iloc[:, 1:]
        region_cellids_by_name: dict[str, list] = {}
        region_layers_by_name: dict[str, int] = {}
        region_geometries_by_name: dict[str, object] = {}
        region_metadata_by_name: dict[str, dict] = {}
        boundary_dict = {}

        for per in range(nper):  # FOR EACH STRESS PERIOD
            cell_list = []

            for name, row, active_cells in self.iter_polygon_boundary_features(
                name_field=fields["name"],
                edges_only=edges_only,  # default False = every cell in the polygon (as
                                        # the legacy get_ghb_from_shp did); True keeps
                                        # only the polygon cells on the model GRID edge
                                        # (get_grid_edge), for edge-boundary polygons.
            ):  # FOR EACH GHB BOUDARY
                # conductance + layer are required; the others are optional and
                # may be absent from the source attributes (e.g. an elevation-only
                # GHB shapefile), so read them defensively -> None when missing.
                conductance = row[fields['conductance']]
                layer = row[fields['layer']]
                boundary_height = row.get(fields['height_over_btm'])
                min_elev = row.get(fields['min_elev'])
                elev = row.get(fields['elevation'])
                min_elev = min_elev if min_elev is not None else 0
                # adjust layer number for zero-based indexing
                layer_idx = layer - 1
                region_layers_by_name[name] = layer_idx
                region_geometries_by_name[name] = row['geometry']
                region_metadata_by_name[name] = {
                    "reference_offset": reference_offset,
                    "uses_elev_reference": name in elev_reference,
                }
                region_cellids_by_name.setdefault(name, [])

                for cell in active_cells:  # FOR EACH CELL IN THIS GHB BOUNDARY
                    if elev is not None:
                        if name in elev_reference.keys():
                            try:
                                e = elev_reference[name][per]
                            except:
                                boundary_elev = elev + reference_offset
                                e = None
                            if isinstance(e, float | int):
                                if self.verbose:
                                    print(f'using {name} as an elevation reference for cell {cell}')
                                boundary_elev = e + reference_offset
                            else:
                                boundary_elev = elev + reference_offset
                        else:
                            boundary_elev = elev
                    elif boundary_height is not None:
                        boundary_elev = lyr_botms.iloc[cell, layer_idx] + boundary_height
                    else:
                        boundary_elev = lyr_botms.iloc[cell, layer_idx]
                    if boundary_elev < min_elev:
                        boundary_elev = min_elev  # adjusts drn elev to minimum allowed if specified
                    cell_id = build_cell_id(cell, grid_type=grid_type, layer=layer_idx)
                    cell_list.append([cell_id, boundary_elev, conductance])
                    region_cellids_by_name[name].append(cell_id)

            boundary_dict[per] = cell_list

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

        return boundary_dict

    def from_vector(self, **kwargs) -> dict:
        """Alias for :meth:`from_polygons` for shapefile/geopackage workflows."""

        return self.from_polygons(**kwargs)

    def get_from_poly(self, **kwargs) -> dict:
        """Backward-compatible alias for :meth:`from_polygons`."""

        return self.from_polygons(**kwargs)

    def add_to_dict(self, existing_dict: dict = None, dict_to_add: dict = None):
        """Merge a validated GHB stress-period dict into ``existing_dict`` (replacing collisions)."""

        assert self._verify_boundary_dict_structure(dict_to_add), 'dict to add not valid or missing'
        if existing_dict is None:
            return dict_to_add
        return merge_stress_period_data(existing_dict, dict_to_add, replace=True)

    def _verify_boundary_dict_structure(self, dict_to_check: dict = None):
        """Validate a boundary dict: every period, layer, and cell id is within the model grid."""

        if dict_to_check is None:
            return dict_to_check
        for per, period_list in dict_to_check.items():
            assert per in range(self.model.nper), f'stress period {per} not in model'
            cell_dict = {cell_list[0]: cell_list[1:] for cell_list in period_list}
            for k in cell_dict.keys():
                assert k[0] in range(self.model.modelgrid.nlay), f'layer {k[0]} not in model grid'
                assert k[1] in range(self.model.modelgrid.ncpl), f'cell {k} not a valid model cell id'
            # TODO check cell_dict values based on type of boundary
        return True


# Preferred alias for the vector-driven general-head-boundary builder API.
GHBFromVector = GHB

