"""Drain-boundary helpers built on polygon/line GIS inputs and Voronoi cells."""

from __future__ import annotations
from typing import TYPE_CHECKING

import flopy.utils.binaryfile
import shapely as shp
from pathlib import Path
import geopandas as gpd
import pandas as pd
import pickle
from simple_modflow.modflow.mf6.boundary_support import (
    build_cell_id,
    coerce_values_by_cell,
    filter_inactive_cells,
    normalize_grid_type,
)
from simple_modflow.modflow.mf6.boundaries import Boundaries

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

idxx = pd.IndexSlice
inches_to_feet = 1 / 12


class DRN(Boundaries):
    """Build drain stress-period data from polygons or explicit cell selections."""

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            shp_gpkg: Path = None,
            uid: str = None,
            crs: int = 2927,
            grid_type: str = 'disv',
            idomain: list[int] | pd.Series = None,
            idomain_path: Path = None,
    ):
        """Parameters
        ----------
        model
            Model to which the boundary applies.
        vor
            Grid helper; defaults to the model grid if omitted.
        shp_gpkg
            Polygon/geometry file describing the drain features.
        uid
            Unique-id field in the geometry attributes.
        crs
            EPSG code for the boundary geometry.
        grid_type
            MODFLOW grid type, usually ``disv`` or ``disu``.
        idomain, idomain_path
            Optional active-domain definition used to filter inactive cells.
        """
        super().__init__(model, vor, shp_gpkg, uid, crs, idomain=idomain, idomain_path=idomain_path)
        self.bound_type = 'drn'
        self.grid_type = normalize_grid_type(grid_type)

    def get_drn_stress_period_data(
            self,
            cells: list,
            bottom_addition: float | int = 0,
            conductance: float | int | list = 100,
            disMf: str = 'disv',
            bottoms: dict = None,
            layer: int = None,
            region_name: str | None = None,
            region_tags: list[str] | None = None,
            region_metadata: dict | None = None,
            overwrite_region: bool = False,
    ) -> list:
        """Returns a list of lists. Each nested list corresponds to the DRN package
        boundary data for a particular voronoi cell in the grid, which includes cell
        ID, elevation of drain, and conductance. Can be passed to the flopy DRN package.

        Args:
            cells (list): list of cell IDs in this drain
            bottoms (dict): dict of bottoms, keys are the cell indices, values are the bottom elevations. If no bottom addition is provided, bottoms is just the elevation of the boundary.
            bottom_addition (float): height above the bottom of cell for the drain. This is added to the bottom of cell elevation derived from the Voronoi grid object.
            conductance (float | int | list): conductance for the cell of this DRN, can also provide a list
            of conductance values that is of equal length as the list of cells. Function will create a
            dict of cell:conductance key:value pairs.
            disMf (str, optional): Either 'disu' or 'disv' works, and refers to the MODFLOW6 discretization package being used. Defaults to 'disu'.
            layer: layer to apply drains

        Returns:
            list: List of lists that contain the data for this drain and can be passed to the flopy DRN package
        """
        if self.vor is None:
            raise ValueError("No voronoi grid defined")

        grid_type = normalize_grid_type(disMf)
        layer_idx = 0 if layer is None else int(layer)
        conductance_by_cell = coerce_values_by_cell(cells, conductance, name="conductance")
        active_cells = filter_inactive_cells(cells, self.inactive_cells)
        drn_values = []

        # generate list of lists for all DRN boundary cells
        for cell in active_cells:
            cell_id = build_cell_id(cell, grid_type=grid_type, layer=layer_idx)
            if bottoms:
                drain_elevation = bottoms[cell] + bottom_addition
            elif self.vor.gdf_topbtm is not None:
                try:
                    drain_elevation = self.vor.gdf_topbtm.loc[cell, layer_idx + 1] + bottom_addition
                except Exception:
                    drain_elevation = bottom_addition
            else:
                drain_elevation = bottom_addition
            drn_values.append([cell_id, drain_elevation, conductance_by_cell[cell]])

        if region_name is not None and drn_values:
            self._register_region(
                region_name,
                cellids=[row[0] for row in drn_values],
                layer=layer_idx,
                tags=region_tags,
                metadata=region_metadata,
                overwrite=overwrite_region,
            )
        return drn_values

    def from_polygons(
            self,
            grid_type: str = 'disv',
            fields: dict = None,
            edges_only: bool = False,
            top_drain: bool = False,
            register_regions: bool = False,
            region_name_prefix: str | None = None,
            combined_region_name: str | None = None,
            region_tags: list[str] | None = None,
            overwrite_regions: bool = False,
            # top_minus = 0
    ) -> dict:
        """Build drain stress-period data from polygon features on the builder.

        Parameters
        ----------
        grid_type
            MODFLOW grid type, usually ``"disv"`` or ``"disu"``.
        fields
            Mapping of logical field names to geometry attribute names.
        edges_only
            When ``True``, only edge-cell intersections are used.
        top_drain
            When ``True``, interpret the height field relative to model top
            instead of cell bottom.
        register_regions, region_name_prefix, combined_region_name, region_tags,
        overwrite_regions
            Optional model-region registration settings.
        """
        vor = self._require_vor()
        nper = self.nper if self.nper is not None else 1
        grid_type = normalize_grid_type(grid_type)
        if fields is None:
            fields = {
                'name': 'name',
                'height_over_btm': 'height',
                'conductance': 'cond',
                'layer': 'layer',
                'min_elev': 'min_elev'
            }
        # get bottoms of model layers from voronoi grid
        lyr_botms = vor.gdf_topbtm.drop('geometry', axis='columns').iloc[:, 1:]
        if top_drain:
            model_top = vor.gdf_topbtm.drop('geometry', axis='columns').iloc[:, 0]
        region_cellids_by_name: dict[str, list] = {}
        region_layers_by_name: dict[str, int] = {}
        region_geometries_by_name: dict[str, object] = {}
        region_metadata_by_name: dict[str, dict] = {}
        drn_dict = {}
        for per in range(nper):
            cell_list = []
            for name, row, active_cells in self.iter_polygon_boundary_features(
                name_field=fields["name"],
                edges_only=edges_only,
            ):
                boundary_height = row[fields['height_over_btm']]
                conductance = row[fields['conductance']]
                layer = row[fields['layer']]
                min_elev = row[fields['min_elev']]
                min_elev = min_elev if min_elev is not None else 0
                # adjust layer number for zero-based indexing
                layer_idx = layer - 1
                region_layers_by_name[name] = layer_idx
                region_geometries_by_name[name] = row["geometry"]
                region_metadata_by_name[name] = {"edges_only": edges_only, "top_drain": top_drain}
                region_cellids_by_name.setdefault(name, [])
                for cell in active_cells:
                    if top_drain:
                        boundary_elev = model_top.iloc[cell] + boundary_height
                    else:
                        boundary_elev = lyr_botms.iloc[cell, layer_idx] + boundary_height
                    if boundary_elev < min_elev:
                        boundary_elev = min_elev  # adjusts drn elev to minimum allowed if specified
                    cell_id = build_cell_id(cell, grid_type=grid_type, layer=layer_idx)
                    cell_list.append([cell_id, boundary_elev, conductance])
                    region_cellids_by_name[name].append(cell_id)
            drn_dict[per] = cell_list

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
        return drn_dict

    def from_vector(self, **kwargs) -> dict:
        """Alias for :meth:`from_polygons` for shapefile/geopackage workflows."""

        return self.from_polygons(**kwargs)

    def get_drn_from_poly(self, **kwargs) -> dict:
        """Backward-compatible alias for :meth:`from_polygons`."""

        return self.from_polygons(**kwargs)

    @staticmethod
    def update_drn_dict(drn_dict: dict, update_dict: dict, update_existing_only: bool = True):
        assert all(key in drn_dict.keys() for key in update_dict.keys()), 'update_dict keys must be in drn_dict keys'
        for key, updater in update_dict.items():
            updated = pd.DataFrame(drn_dict[key]).set_index(0)
            updated_idx = list(updated.index)
            updater = pd.DataFrame(updater).set_index(0)
            updater_idx = list(updater.index)
            updater_idx = [idx for idx in updater_idx if idx in updated_idx]
            updated.loc[updater_idx] = updater.loc[updater_idx]  # update the Dataframe with new elevs and conductances
            if update_existing_only is False:
                new_idx = [idx for idx in list(updater.index) if idx not in updated_idx]
                updated = pd.concat([updated, updater.loc[new_idx]])
            updated = updated.reset_index()
            updated = updated.to_numpy().tolist()  # recreate list then update the dict
            drn_dict[key] = updated
        return drn_dict


# Preferred alias for the vector-driven drain builder API.
DRNFromVector = DRN
