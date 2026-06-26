"""Voronoi grid wrapper and convenience helpers for MODFLOW workflows.

``VoronoiGridPlus`` extends FloPy's Voronoi grid representation with cached
GeoDataFrame views, DISU connectivity helpers, plotting shortcuts, raster-based
surface helpers, and selection utilities used throughout ``myflopy``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid

from myflopy.modflow.mf6.grid.connectivity import build_disu_connectivity
from myflopy.modflow.mf6.grid.geometry import (
    adjust_cells_by_id as geometry_adjust_cells_by_id,
    adjust_top_btm_overlaps as geometry_adjust_top_btm_overlaps,
    calculate_distance as geometry_calculate_distance,
    find_adjacent_cells as geometry_find_adjacent_cells,
    find_adjacent_polygons as geometry_find_adjacent_polygons,
    generate_grid_around_point as geometry_generate_grid_around_point,
    generate_grid_coordinates as geometry_generate_grid_coordinates,
    generate_grid_polygons as geometry_generate_grid_polygons,
    get_centroids as geometry_get_centroids,
    get_gdf_latlon as geometry_get_gdf_latlon,
    get_gdf_vor_polys as geometry_get_gdf_vor_polys,
    get_grid_centroid as geometry_get_grid_centroid,
    get_latlon as geometry_get_latlon,
    get_overlapping_area as geometry_get_overlapping_area,
    get_vor_idx_from_geometry as geometry_get_vor_idx_from_geometry,
    get_vor_idx_from_geometry_idx as geometry_get_vor_idx_from_geometry_idx,
    get_voronoi_polygons as geometry_get_voronoi_polygons,
    reconcile_surfaces as geometry_reconcile_surfaces,
    set_k_vor as geometry_set_k_vor,
    shared_face_length as geometry_shared_face_length,
    voronoi_refine_by_point as geometry_voronoi_refine_by_point,
)
from myflopy.modflow.mf6.grid.helpers import get_griddata_from_disu
from myflopy.modflow.mf6.grid.plotting import (
    build_choropleth as plotting_build_choropleth,
    build_grid_section as plotting_build_grid_section,
    get_dash_selector as plotting_get_dash_selector,
    map_nodes as plotting_map_nodes,
    mapit as plotting_mapit,
    plot2d as plotting_plot2d,
    plot3d as plotting_plot3d,
    plottri as plotting_plottri,
    show as plotting_show,
    show_overlapping_geometry as plotting_show_overlapping_geometry,
    show_selected_cells as plotting_show_selected_cells,
)
from myflopy.modflow.mf6.grid.selection import (
    get_grid_edge_cells,
    get_model_boundary_polygons,
    get_vor_cells_as_dict as selection_get_vor_cells_as_dict,
    get_vor_cells_as_series as selection_get_vor_cells_as_series,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_cell_areas as surface_get_cell_areas,
    get_domain as surface_get_domain,
    get_gdf_topbtm_multilyr as surface_get_gdf_topbtm_multilyr,
    get_normal_from_strike_and_dip as surface_get_normal_from_strike_and_dip,
    get_origin_xy as surface_get_origin_xy,
    get_raster_from_strike_dip as surface_get_raster_from_strike_dip,
    get_raster_vals_at_centroids as surface_get_raster_vals_at_centroids,
)

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


class VoronoiGridPlus(VoronoiGrid):
    """An unstructured Voronoi (DISV) grid with MODFLOW-modeling conveniences.

    Extends FloPy's ``VoronoiGrid`` with the geometry and helpers myflopy needs to
    build and post-process unstructured models: DISV grid properties
    (:meth:`get_disv_gridprops` -> ``ncpl``/``vertices``/``cell2d``), DISU
    connectivity (``iac``/``ja``/``cl12``/``hwva``), cell-centroid and adjacency
    lookups, CRS handling, and the per-cell ``gdf_vorPolys`` / ``gdf_topbtm``
    GeoDataFrames the surface-aware builders read.

    Build one from a :class:`~myflopy.modflow.mf6.grid.triangle.TriangleGrid`
    triangulation, or declaratively from GeoPackages via
    ``mf.GridSpec.voronoi(...).resolve(ws)``. Put the result in a
    :class:`~myflopy.specs.ModelContext` so the package-first GIS helpers can map
    features onto its cells, and pass ``get_disv_gridprops()`` to ``mf.disv``.

    Parameters
    ----------
    tri
        A built ``Triangle`` triangulation defining the Voronoi tessellation.
    crs
        Coordinate reference system (default ``"EPSG:2927"``); the single source of
        truth for aligning vector/raster inputs.
    rasters
        Optional raster(s) sampled for cell elevations/properties.
    """

    def __init__(
        self,
        tri: Triangle = None,
        crs: str = 'EPSG:2927',
        rasters: list | Path = None,
        name: str = 'voronoi_grid',
        qhull_options: str = None,
        idomain: list = None,
        idomain_path: Path = None,
        verts: np.ndarray = None,
        iverts: list[list] = None,
        xcyc: np.ndarray = None,
        **kwargs,
    ):
        """Create a Voronoi view from Triangle output or DISU-style arrays.

        Parameters
        ----------
        tri
            Triangle mesh used to generate the Voronoi dual.
        crs
            Coordinate reference system string for geometry outputs.
        rasters
            Optional list of raster paths used later for top/bottom sampling.
        name
            Friendly grid name for diagnostics and exports.
        qhull_options
            Optional QHull options passed through to FloPy's Voronoi builder.
        idomain, idomain_path
            Optional active-cell selections used by edge-detection helpers.
        verts, iverts, xcyc
            Explicit vertex-grid arrays used when reconstructing from DISU data
            instead of from a Triangle mesh.
        """
        print("VoronoiGrid initializing.")
        if tri:
            # FloPy's tri2vor helper currently emits harmless RuntimeWarnings from
            # point_in_polygon when domain edges are horizontal. Suppress only that
            # known warning pattern so realistic Triangle-built grids stay quiet.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="invalid value encountered in divide",
                    category=RuntimeWarning,
                    module=r"flopy\.utils\.geometry",
                )
                warnings.filterwarnings(
                    "ignore",
                    message="divide by zero encountered in divide",
                    category=RuntimeWarning,
                    module=r"flopy\.utils\.geometry",
                )
                super().__init__(tri, qhull_options=qhull_options, **kwargs)
        elif all(v is not None for v in [verts, iverts, xcyc]):
            self.verts = verts
            self.iverts = iverts
            self.points = xcyc
            self.ncpl = len(iverts)
            self.nverts = verts.shape[0]
        else:
            raise ValueError('Either tri or (verts, iverts, xcyc) must be provided.')

        self.tri = tri
        self.crs_latlon = "EPSG:4326"
        self._centroids = None
        self._gdf_latlon = None
        self._nlay = None
        self._latlon = None
        self.rasters = rasters
        self.crs = crs
        self.name = name
        self.x_coords_by_node = []
        self.y_coords_by_node = []
        self._gdf_vorPolys = None
        self._iac = None
        self._nja = None
        self._idomain_path = None
        self._idomain = None
        self._gdf_topbtm = None
        self._area_list = None
        self._adjacent_cells_idx = None
        self._ja, self._cl12, self._hwva = None, None, None
        self.centroids_x, self.centroids_y = self.get_centroids()
        self.idomain_path = idomain_path
        self.config = {'scrollZoom': True}
        self.scatt_layout = {'height': 1000, 'width': 1000, 'dragmode': 'pan'}
        self.grid_centroid = self.get_grid_centroid()
        print('Voronoi grid initialized.')

    # Attributes dropped from pickles: the Triangle builder plus caches that
    # recompute from the stored mesh geometry. Keeping them would bloat the
    # pickle and tie it more tightly to dependency versions.
    _PICKLE_VOLATILE = (
        "tri",
        "_centroids",
        "_gdf_latlon",
        "_latlon",
        "_gdf_vorPolys",
        "_iac",
        "_nja",
        "_ja",
        "_cl12",
        "_hwva",
        "_area_list",
        "_adjacent_cells_idx",
    )

    def __getstate__(self):
        """Return a lean picklable state, dropping ``tri`` and recomputable caches."""

        return {
            key: value
            for key, value in self.__dict__.items()
            if key not in self._PICKLE_VOLATILE
        }

    def __setstate__(self, state):
        """Restore state, leaving dropped caches empty for lazy recomputation."""

        self.__dict__.update(state)
        for key in self._PICKLE_VOLATILE:
            self.__dict__.setdefault(key, None)

    @property
    def vor_list(self):
        """Return the Voronoi polygons as a plain geometry list."""
        return self.gdf_vorPolys.geometry.to_list()

    @property
    def cell_list(self):
        """Return the cell indices for the Voronoi grid."""
        return self.gdf_vorPolys.index.to_list()

    @property
    def area_list(self):
        """Return cached polygon areas for all Voronoi cells."""
        if self._area_list is None:
            self._area_list = [cell.area for cell in self.vor_list]
        return self._area_list

    @property
    def idomain_path(self):
        """Path used by some edge helpers to read external idomain data."""
        return self._idomain_path

    @idomain_path.setter
    def idomain_path(self, value):
        """Set the optional idomain path used by edge-detection helpers."""
        if value is not None:
            assert isinstance(value, Path), 'idomain path must be a Path instance'
        self._idomain_path = value

    @property
    def idomain(self):
        """List of active Voronoi cell indices, if one has been provided."""
        return self._idomain

    @idomain.setter
    def idomain(self, value):
        """Set the active-cell list used by some edge-detection helpers."""
        if value is not None:
            assert isinstance(value, list), 'idomain must be a list'
            assert all(idx in self.cell_list for idx in value), 'idomain indices must be valid Voronoi grid cells'
        self._idomain = value

    @property
    def nlay(self):
        """Infer the number of layers from ``gdf_topbtm``."""
        if self._nlay is None:
            self._nlay = len(self.gdf_topbtm.drop('geometry', axis=1).columns) - 1
        return self._nlay

    @property
    def gdf_topbtm(self):
        """Return cached top/bottom surfaces, building them from rasters if needed."""
        if self._gdf_topbtm is None and self.rasters is not None:
            self._gdf_topbtm = self.get_gdf_topbtm_multilyr(rasters=self.rasters)
        return self._gdf_topbtm

    @gdf_topbtm.setter
    def gdf_topbtm(self, value):
        """Set the cached top/bottom surface GeoDataFrame."""
        self._gdf_topbtm = value

    def get_disu_connectivity(self, *, use_representative_point: bool = False, tol: float = 0.0, validate: bool = True):
        """Compute DISU-style connectivity arrays for the current Voronoi grid."""
        iac, ja, cl12, hwva, nja = build_disu_connectivity(
            self.gdf_vorPolys,
            self.adjacent_cells_idx,
            use_representative_point=use_representative_point,
            tol=tol,
            validate=validate,
        )
        self._iac = iac
        self._ja = ja
        self._cl12 = cl12
        self._hwva = hwva
        self._nja = nja
        return iac, ja, cl12, hwva, nja

    @property
    def iac(self):
        """Return the cached DISU ``iac`` array, computing it if needed."""
        if self._iac is None:
            self.get_disu_connectivity()
        return self._iac

    @property
    def ja(self):
        """Return the cached DISU ``ja`` array, computing it if needed."""
        if self._ja is None:
            self.get_disu_connectivity()
        return self._ja

    @property
    def cl12(self):
        """Return the cached half-cell connection lengths."""
        if self._cl12 is None:
            self.get_disu_connectivity()
        return self._cl12

    @property
    def hwva(self):
        """Return the cached shared-face widths/areas."""
        if self._hwva is None:
            self.get_disu_connectivity()
        return self._hwva

    @property
    def nja(self):
        """Return the cached number of ``ja`` entries."""
        if self._nja is None:
            self.get_disu_connectivity()
        return self._nja

    @property
    def gdf_vorPolys(self):
        """Return the Voronoi polygons as a GeoDataFrame."""
        if self._gdf_vorPolys is None:
            self._gdf_vorPolys = self.get_gdf_vorPolys(crs=self.crs)
        return self._gdf_vorPolys

    @gdf_vorPolys.setter
    def gdf_vorPolys(self, value):
        """Set the cached Voronoi polygon GeoDataFrame."""
        self._gdf_vorPolys = value

    @property
    def adjacent_cells_idx(self):
        """Return cached adjacency lists for each Voronoi cell."""
        if self._adjacent_cells_idx is None:
            self._adjacent_cells_idx = self.find_adjacent_polygons(self.gdf_vorPolys)
        return self._adjacent_cells_idx

    def get_disu_props(self):
        """Populate and cache the standard DISU connectivity properties."""
        self.iac
        self.ja
        self.nja
        self.cl12
        self.hwva

    def get_voronoi_polygons(self):
        """Return the Voronoi polygons as shapely geometry objects."""
        return geometry_get_voronoi_polygons(self)

    def mapit(self, crs='EPSG:2927'):
        """Open a quick interactive map view of the Voronoi polygons."""
        return plotting_mapit(self, crs=crs)

    def choropleth(
        self,
        model: SimulationBase = None,
        kstpkper: tuple = None,
        per: int = None,
        layer: int = 0,
        type: str = 'hds',
        custom_hover: dict = None,
        custom_zs: list = None,
        zmin: float | int = None,
        zmax: float | int = None,
        zoom: int = 13,
        show_layer_elevs: bool = False,
        show_mounding: bool = False,
        hover_heads: bool = True,
        hover_ks: bool = False,
        locs: Path = None,
    ):
        """Build the standard Voronoi choropleth wrapper used by the package."""
        return plotting_build_choropleth(
            self,
            model=model,
            kstpkper=kstpkper,
            per=per,
            layer=layer,
            type=type,
            custom_hover=custom_hover,
            custom_zs=custom_zs,
            zmin=zmin,
            zmax=zmax,
            zoom=zoom,
            show_layer_elevs=show_layer_elevs,
            show_mounding=show_mounding,
            hover_heads=hover_heads,
            hover_ks=hover_ks,
            locs=locs,
        )

    @property
    def dash_selector(self):
        """Return a Dash-compatible selector for choropleth exploration."""
        return plotting_get_dash_selector(self)

    def show(self):
        """Display the default choropleth view."""
        return plotting_show(self)

    def map_nodes(self):
        """Display cell polygons with node IDs in an interactive map."""
        return plotting_map_nodes(self)

    def get_vor_cells_as_series(
        self,
        overlapping_geometry: shp.Polygon | shp.Point | gpd.GeoSeries | Path = None,
        predicate: str = 'intersects',
        return_dict: bool = False,
        name_field: str = 'ExploName',
    ) -> pd.Series | dict:
        """Map overlapping geometry to intersected Voronoi cell IDs."""
        return selection_get_vor_cells_as_series(
            self.gdf_vorPolys,
            overlapping_geometry=overlapping_geometry,
            predicate=predicate,
            return_dict=return_dict,
            name_field=name_field,
        )

    def get_vor_cells_as_dict(
        self,
        locs: Path,
        crs: str = None,
        predicate: str = 'intersects',
        loc_name_field: str = None,
        return_gdf: bool = False,
    ) -> dict:
        """Return intersected Voronoi cell IDs keyed by feature name."""
        return selection_get_vor_cells_as_dict(
            self,
            locs=locs,
            crs=crs,
            predicate=predicate,
            loc_name_field=loc_name_field,
            return_gdf=return_gdf,
        )

    def show_selected_cells(self, cell_list: list = None, **kwargs):
        """Show a quick map of a selected list of cells."""
        return plotting_show_selected_cells(self, cell_list=cell_list, **kwargs)

    def show_overlapping_geometry(self, shp_gpkg):
        """Show cells that overlap the provided geometry."""
        return plotting_show_overlapping_geometry(self, shp_gpkg)

    def get_model_boundary_polygons(self) -> dict:
        """Return polygon geometries for the model boundary components."""
        return get_model_boundary_polygons(self.gdf_vorPolys)

    def get_grid_edge(self, idomain: list = None, idomain_path: Path = None, include_interiors: bool = True) -> list:
        """Return cell IDs that lie on the active-model edge."""
        return get_grid_edge_cells(
            self,
            idomain=idomain,
            idomain_path=idomain_path,
            include_interiors=include_interiors,
        )

    def plot3d(self, z=None):
        """Plot the grid as a simple 3D line or centroid view."""
        return plotting_plot3d(self, z=z)

    def plot2d(self):
        """Plot the Voronoi edges in 2D using Plotly."""
        return plotting_plot2d(self)

    def plottri(self):
        """Plot the underlying Triangle mesh used to build the Voronoi grid."""
        return plotting_plottri(self)

    def generate_grid_coordinates(self, grid_spacing: int) -> list:
        """Generate regularly spaced coordinates over the grid extent."""
        return geometry_generate_grid_coordinates(self, grid_spacing)

    @property
    def centroids(self):
        """Return cached x/y centroid coordinate arrays."""
        if self._centroids is None:
            self._centroids = self.get_centroids()
        return self._centroids

    def get_centroids(self) -> tuple[list[float], list[float]]:
        """Compute x and y centroid coordinates for all Voronoi cells."""
        return geometry_get_centroids(self)

    def find_adjacent_polygons(self, gdf: gpd.GeoDataFrame) -> list:
        """Return adjacency lists for the polygons in ``gdf``."""
        return geometry_find_adjacent_polygons(gdf)

    def find_adjacent_cells(self, cell_id):
        """Return neighboring cell IDs for one Voronoi cell."""
        return geometry_find_adjacent_cells(self, cell_id)

    def calculate_distance(self, gdf, poly_idx1, poly_idx2):
        """Calculate the distance between two polygons in a GeoDataFrame."""
        return geometry_calculate_distance(gdf, poly_idx1, poly_idx2)

    def shared_face_length(self, poly1, poly2):
        """Return the length of the shared face between two polygons."""
        return geometry_shared_face_length(poly1, poly2)

    def get_gdf_vorPolys(self, crs=None):
        """Build a GeoDataFrame of Voronoi polygons."""
        return geometry_get_gdf_vor_polys(self, crs=crs)

    def get_gdf_topbtm(self, rasters: list, labels: list = None):
        """Backward-compatible alias for ``get_gdf_topbtm_multilyr``."""
        return self.get_gdf_topbtm_multilyr(rasters=rasters, labels=labels)

    def get_raster_vals_at_centroids(
        self,
        raster_files: Path | list[Path | int] = None,
        labels: str | list[str] = None,
    ) -> gpd.GeoDataFrame | None:
        """Sample raster values at Voronoi cell centroids."""
        return surface_get_raster_vals_at_centroids(self, raster_files=raster_files, labels=labels)

    def get_gdf_topbtm_multilyr(self, rasters: list, labels: list = None):
        """Build top/bottom surfaces for one or more layers from rasters."""
        return surface_get_gdf_topbtm_multilyr(self, rasters=rasters, labels=labels)

    def get_cell_areas(self):
        """Return per-cell polygon areas."""
        return surface_get_cell_areas(self)

    def get_origin_xy(self):
        """Return an origin point suitable for plotting or raster work."""
        return surface_get_origin_xy(self)

    @property
    def gdf_latlon(self):
        """Return the Voronoi polygons reprojected to latitude/longitude."""
        return geometry_get_gdf_latlon(self)

    @property
    def latlon(self):
        """Return centroid coordinates in latitude/longitude."""
        return geometry_get_latlon(self)

    def get_grid_centroid(self):
        """Return the centroid of the overall grid extent."""
        return geometry_get_grid_centroid(self)

    def get_overlapping_area(self, shp_gpkg=None, cell_list=None):
        """Calculate overlap area between cells and external geometry."""
        return geometry_get_overlapping_area(self, shp_gpkg=shp_gpkg, cell_list=cell_list)

    def get_vor_idx_from_geometry(
        self,
        shp_to_query: shp = None,
        gdf_to_query: gpd.GeoDataFrame = None,
        crs: str = "EPSG:2927",
        crs_latlon: str = "EPSG:4326",
        name_col: str = None,
        predicate: str = "intersects",
    ) -> dict:
        """Map external geometry to intersecting Voronoi cell IDs."""
        return geometry_get_vor_idx_from_geometry(
            self,
            shp_to_query=shp_to_query,
            gdf_to_query=gdf_to_query,
            crs=crs,
            crs_latlon=crs_latlon,
            name_col=name_col,
            predicate=predicate,
        )

    def get_vor_idx_from_geometry_idx(
        self,
        gdf_to_query: gpd.GeoDataFrame = None,
        idx: int = 0,
        predicate: str = "intersects",
    ) -> list:
        """Return intersecting cell IDs for one row of a GeoDataFrame."""
        return geometry_get_vor_idx_from_geometry_idx(
            self,
            gdf_to_query=gdf_to_query,
            idx=idx,
            predicate=predicate,
        )

    def set_k_vor(self, k_dict: dict = None, k_default=100) -> list:
        """Build a per-cell hydraulic-conductivity list from region mappings."""
        return geometry_set_k_vor(self, k_dict=k_dict, k_default=k_default)

    def get_raster_from_strike_dip(
        self,
        strike: int,
        dip: int,
        known_point: tuple,
        pixel_size: int = 1,
        output_filename: Path = Path.cwd().joinpath('raster.tif'),
    ):
        """Create a raster surface from strike/dip assumptions."""
        return surface_get_raster_from_strike_dip(
            self,
            strike=strike,
            dip=dip,
            known_point=known_point,
            pixel_size=pixel_size,
            output_filename=output_filename,
        )

    def to_shapefile(self, filepath: str | Path = 'vor_shp.shp'):
        """Write the Voronoi polygons to a shapefile or other GeoPandas target."""
        return self.gdf_vorPolys.to_file(filepath)

    def ugrid2d(self):
        """Return this grid's topology as an :class:`xugrid.Ugrid2d` (no data).

        Builds the UGRID 2-D unstructured mesh -- nodes (cell vertices), faces
        (cells), and the face-node connectivity -- from :meth:`get_disv_gridprops`,
        carrying the grid :attr:`crs`. This is the shared topology builder behind
        :meth:`to_xugrid` and ``HeadsPlus.to_xugrid`` / ``model.to_xugrid``; call
        it directly when you want to attach your own ``xarray`` data on the face
        dimension (``grid.face_dimension``).

        Returns
        -------
        xugrid.Ugrid2d
            The mesh topology, with ``n_face == ncpl`` and ``n_node == nvert``.

        Raises
        ------
        ImportError
            If the optional ``xugrid`` package is not installed.
        """

        try:
            import xugrid as xu
        except ImportError as err:  # pragma: no cover - optional dependency
            raise ImportError(
                "ugrid2d() requires the optional 'xugrid' package. Install it "
                "with `pip install xugrid xarray`."
            ) from err

        gp = self.get_disv_gridprops()
        nvert = int(gp["nvert"])
        ncpl = int(gp["ncpl"])

        # Node coordinates, placed at their vertex id so connectivity lines up.
        node_x = np.full(nvert, np.nan, dtype=float)
        node_y = np.full(nvert, np.nan, dtype=float)
        for iv, x, y in gp["vertices"]:
            node_x[int(iv)] = float(x)
            node_y[int(iv)] = float(y)

        # Ragged face -> node connectivity, padded to a rectangular array.
        fill_value = -1
        rows = [
            [int(v) for v in cell[4:4 + int(cell[3])]]
            for cell in gp["cell2d"]
        ]
        max_nodes = max(len(row) for row in rows)
        face_node_connectivity = np.full((ncpl, max_nodes), fill_value, dtype=np.int64)
        for i, row in enumerate(rows):
            face_node_connectivity[i, : len(row)] = row

        return xu.Ugrid2d(
            node_x,
            node_y,
            fill_value,
            face_node_connectivity,
            name="mesh2d",
            is_projected=True,
            crs=self.crs,
        )

    def to_xugrid(self, data=None, *, name: str = "data", layer_dim: str = "layer"):
        """Export this DISV grid (and optional per-cell data) as an xugrid object.

        Builds a UGRID 2-D unstructured-mesh representation of the Voronoi grid --
        nodes (cell vertices), faces (cells), and the face-node connectivity --
        wrapped in an `xugrid <https://deltares.github.io/xugrid/>`_ object. xugrid
        is xarray extended for unstructured grids, so the result plugs directly
        into xarray-style analysis, unstructured plotting/cross-sections, and a
        UGRID-NetCDF export that QGIS, ParaView, and other tools can open -- which
        makes DISV results far easier to share than raw MODFLOW binaries.

        The grid topology is taken from :meth:`get_disv_gridprops` and the grid's
        :attr:`crs` is carried onto the mesh.

        Parameters
        ----------
        data
            Optional per-cell values to attach on the face dimension. One of:

            - ``None`` (default) -- return just the mesh topology;
            - a 1-D array of length ``ncpl`` -- one field, e.g. a head layer
              (``model.hds.array(layer=0)``) or a K array;
            - a 2-D array shaped ``(nlay, ncpl)`` -- a layered field, given the
              extra ``layer_dim`` dimension;
            - a ``dict`` of ``{variable_name: array}`` -- several fields at once
              (each 1-D or 2-D as above).
        name
            Variable name used when ``data`` is a single array (ignored for a
            ``dict``). Defaults to ``"data"``.
        layer_dim
            Dimension name given to the first axis of 2-D ``(nlay, ncpl)`` inputs.

        Returns
        -------
        xugrid.UgridDataArray | xugrid.UgridDataset
            A ``UgridDataArray`` when ``data`` is a single array, otherwise a
            ``UgridDataset`` (also for ``data=None`` -- topology only).

        Raises
        ------
        ImportError
            If the optional ``xugrid`` / ``xarray`` packages are not installed.
        ValueError
            If an input array's cell axis does not match ``ncpl``.

        Examples
        --------
        >>> uda = vor.to_xugrid(model.hds.array(layer=0), name="head")
        >>> uda.ugrid.plot()                      # unstructured choropleth
        >>> uda.ugrid.to_netcdf("heads.nc")       # UGRID NetCDF for QGIS/ParaView

        >>> # several fields, including a layered one
        >>> uds = vor.to_xugrid({"head": heads_2d, "k": k_2d})   # (nlay, ncpl)
        >>> uds.ugrid.to_netcdf("model.nc")
        """

        try:
            import xarray as xr
            import xugrid as xu
        except ImportError as err:  # pragma: no cover - optional dependency
            raise ImportError(
                "to_xugrid() requires the optional 'xugrid' and 'xarray' "
                "packages. Install them with `pip install xugrid xarray`."
            ) from err

        grid = self.ugrid2d()
        ncpl = int(self.ncpl)
        face_dim = grid.face_dimension

        def _face_dataarray(values, varname):
            array = np.asarray(values, dtype=float)
            if array.ndim == 1:
                if array.shape[0] != ncpl:
                    raise ValueError(
                        f"to_xugrid(): '{varname}' has length {array.shape[0]}, "
                        f"expected ncpl={ncpl}."
                    )
                dims = (face_dim,)
            elif array.ndim == 2:
                if array.shape[-1] != ncpl:
                    raise ValueError(
                        f"to_xugrid(): '{varname}' last axis is {array.shape[-1]}, "
                        f"expected ncpl={ncpl}."
                    )
                dims = (layer_dim, face_dim)
            else:
                raise ValueError(
                    f"to_xugrid(): '{varname}' must be 1-D (ncpl,) or 2-D "
                    f"(nlay, ncpl); got {array.ndim} dimensions."
                )
            return xr.DataArray(array, dims=dims, name=varname)

        if data is None:
            return xu.UgridDataset(grids=[grid])
        if isinstance(data, dict):
            dataset = xr.Dataset(
                {key: _face_dataarray(values, key) for key, values in data.items()}
            )
            return xu.UgridDataset(dataset, grids=[grid])
        return xu.UgridDataArray(_face_dataarray(data, name), grid)

    def get_domain(self):
        """Return the overall polygonal domain of the Voronoi grid."""
        return surface_get_domain(self)

    @staticmethod
    def get_normal_from_strike_and_dip(strike: int, dip: int) -> np.array:
        """Return a unit normal vector from strike and dip."""
        return surface_get_normal_from_strike_and_dip(strike, dip)

    @staticmethod
    def generate_grid_around_point(center_point: shp.Point, spacing: float, size: int, crs: str) -> gpd.GeoSeries:
        """Generate a simple regular grid around a center point."""
        return geometry_generate_grid_around_point(center_point, spacing, size, crs)

    @staticmethod
    def generate_grid_polygons(center_point, spacing, size, gap):
        """Generate regular polygons around a center point."""
        return geometry_generate_grid_polygons(center_point, spacing, size, gap)

    @staticmethod
    def voronoi_refine_by_point(point: shp.Point, spacing: int, tri: Triangle) -> gpd.GeoDataFrame:
        """Generate local refinement geometry around a point."""
        return geometry_voronoi_refine_by_point(point, spacing, tri)

    def reconcile_surfaces(self, df: pd.DataFrame = None, min_sep=0.1, trigger_sep=1, which='bottom'):
        """Reconcile top/bottom surfaces so adjacent layers do not overlap."""
        return geometry_reconcile_surfaces(
            self,
            df=df,
            min_sep=min_sep,
            trigger_sep=trigger_sep,
            which=which,
        )

    def adjust_cells_by_id(
        self,
        cell_ids: list,
        adjustment: int | float,
        df: pd.DataFrame = None,
        layer: int = 0,
        reconcile: bool = True,
    ):
        """Adjust elevations for selected cells by a fixed amount."""
        return geometry_adjust_cells_by_id(
            self,
            cell_ids=cell_ids,
            adjustment=adjustment,
            df=df,
            layer=layer,
            reconcile=reconcile,
        )

    def adjust_top_btm_overlaps(
        self,
        elev_df=None,
        shp: Path = None,
        buffer=1,
        layer_bottom_name=1,
        min_sep=None,
    ) -> pd.DataFrame:
        """Adjust overlapping top/bottom surfaces within selected geometry."""
        return geometry_adjust_top_btm_overlaps(
            self,
            elev_df=elev_df,
            shp_path=shp,
            buffer=buffer,
            layer_bottom_name=layer_bottom_name,
            min_sep=min_sep,
        )

    def cross_section(self, line: shp.LineString | Path):
        """Build a cross-section helper for the provided line."""
        return plotting_build_grid_section(self, line=line)

    @classmethod
    def vor_from_disu(
        cls,
        disu_path: Path,
        crs: str = 'EPSG:2927',
        rasters: list | Path = None,
        name: str = 'voronoi_grid',
        qhull_options: str = None,
        idomain: list = None,
        idomain_path: Path = None,
    ):
        """Reconstruct a Voronoi-style grid view from a DISU input file."""
        verts, iverts, xcyc = get_griddata_from_disu(disu_path)
        return cls(
            verts=verts,
            iverts=iverts,
            xcyc=xcyc,
            crs=crs,
            name=name,
            rasters=rasters,
            idomain=idomain,
            idomain_path=idomain_path,
            qhull_options=qhull_options,
        )
