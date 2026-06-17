"""Constant-head builders that translate vector features into MF6 CHD inputs."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.boundaries import Boundaries
from myflopy.modflow.mf6.boundary_support import build_cell_id, normalize_grid_type

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase


class CHDFromVector(Boundaries):
    """Build MF6 CHD stress-period data from polygon or geopackage features."""

    def __init__(
        self,
        model: "SimulationBase" = None,
        vor: "Vor" = None,
        shp_gpkg: Path = None,
        uid: str = "name",
        crs: int = 2927,
        idomain: list[int] | pd.Series = None,
        idomain_path: Path = None,
    ):
        """Parameters
        ----------
        model
            Model to which the CHD package will be attached.
        vor
            Voronoi grid helper used to intersect vector features with cells.
        shp_gpkg
            Polygon or geopackage path describing CHD features.
        uid
            Unique-id field in the feature attributes.
        crs
            EPSG code for the feature geometry.
        idomain, idomain_path
            Optional active-domain definition used to filter inactive cells.
        """

        super().__init__(model, vor, shp_gpkg, uid, crs, idomain=idomain, idomain_path=idomain_path)
        self.bound_type = "chd"

    def from_polygons(
        self,
        grid_type: str = "disv",
        fields: dict | None = None,
        head_reference: dict | None = None,
        reference_offset: float | int = 0,
        edges_only: bool = False,
        register_regions: bool = False,
        region_name_prefix: str | None = None,
        combined_region_name: str | None = None,
        region_tags: list[str] | None = None,
        overwrite_regions: bool = False,
    ) -> dict:
        """Build CHD stress-period data from polygon features.

        Parameters
        ----------
        grid_type
            MODFLOW grid type, usually ``"disv"`` or ``"disu"``.
        fields
            Mapping of logical field names to geometry attribute names. Defaults
            to ``{"name": "name", "elevation": "elev", "layer": "layer"}``.
        head_reference
            Optional per-feature or per-period override dictionary. When
            provided, values from ``head_reference[name]`` replace the static
            elevation field for that feature. A scalar applies to every period;
            a list applies per stress period.
        reference_offset
            Offset added to any value drawn from ``head_reference``.
        edges_only
            When ``True``, only grid-edge intersections are used.
        register_regions, region_name_prefix, combined_region_name, region_tags,
        overwrite_regions
            Optional model-region registration settings.

        Returns
        -------
        dict
            MF6 CHD stress-period data keyed by stress period.
        """

        grid_type = normalize_grid_type(grid_type)
        nper = self.nper if self.nper is not None else 1
        head_reference = {} if head_reference is None else head_reference
        if fields is None:
            fields = {"name": "name", "elevation": "elev", "layer": "layer"}

        region_cellids_by_name: dict[str, list] = {}
        region_layers_by_name: dict[str, int] = {}
        region_geometries_by_name: dict[str, object] = {}
        region_metadata_by_name: dict[str, dict] = {}
        chd_dict = {}

        for per in range(nper):
            cell_list = []
            for name, row, active_cells in self.iter_polygon_boundary_features(
                name_field=fields["name"],
                edges_only=edges_only,
            ):
                layer_idx = int(row[fields["layer"]]) - 1
                head = row[fields["elevation"]]
                if name in head_reference:
                    override = head_reference[name]
                    if isinstance(override, (int, float)):
                        head = override + reference_offset
                    else:
                        try:
                            head = override[per] + reference_offset
                        except Exception:
                            head = head + reference_offset
                region_layers_by_name[name] = layer_idx
                region_geometries_by_name[name] = row["geometry"]
                region_metadata_by_name[name] = {
                    "edges_only": edges_only,
                    "uses_head_reference": name in head_reference,
                    "reference_offset": reference_offset,
                }
                region_cellids_by_name.setdefault(name, [])

                for cell in active_cells:
                    cell_id = build_cell_id(cell, grid_type=grid_type, layer=layer_idx)
                    cell_list.append([cell_id, head])
                    region_cellids_by_name[name].append(cell_id)
            chd_dict[per] = cell_list

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

        return chd_dict

    def from_vector(self, **kwargs) -> dict:
        """Alias for :meth:`from_polygons` for shapefile/geopackage workflows."""

        return self.from_polygons(**kwargs)

    def get_from_poly(self, **kwargs) -> dict:
        """Backward-compatible alias for :meth:`from_polygons`."""

        return self.from_polygons(**kwargs)
