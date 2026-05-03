from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase

import numpy as np
from scipy.spatial import Voronoi
import pandas as pd
import geopandas as gpd
from flopy.utils.voronoi import VoronoiGrid, tri2vor
from flopy.utils.triangle import Triangle
from flopy.utils.geospatial_utils import GeoSpatialUtil
import shapely as shp
from shapely.geometry import Polygon, MultiLineString, Point, LineString
from pathlib import Path
from simple_modflow.modflow.mf6.grid.helpers import (
    densify_poly,
    flatten,
    get_griddata_from_disu,
    rows_truncate_at_first_missing,
    signed_area,
)
from simple_modflow.modflow.mf6.grid.connectivity import build_disu_connectivity
from simple_modflow.modflow.mf6.grid.selection import (
    get_grid_edge_cells,
    get_model_boundary_polygons,
    get_vor_cells_as_dict as selection_get_vor_cells_as_dict,
    get_vor_cells_as_series as selection_get_vor_cells_as_series,
)
from simple_modflow.modflow.mf6.grid.geometry import (
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
from simple_modflow.modflow.mf6.grid.plotting import (
    GridSection,
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
from simple_modflow.modflow.mf6.grid.surfaces import (
    get_cell_areas as surface_get_cell_areas,
    get_domain as surface_get_domain,
    get_gdf_topbtm_multilyr as surface_get_gdf_topbtm_multilyr,
    get_normal_from_strike_and_dip as surface_get_normal_from_strike_and_dip,
    get_origin_xy as surface_get_origin_xy,
    get_raster_from_strike_dip as surface_get_raster_from_strike_dip,
    get_raster_vals_at_centroids as surface_get_raster_vals_at_centroids,
)
from simple_modflow.modflow.mf6.grid.triangle import TriangleGrid

class VoronoiGridPlus(VoronoiGrid):
    def __init__(self,
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

                 **kwargs
                 ):
        """
        Initializes a Voronoi grid object with all necessary parameters and configurations.

        The constructor sets up a Voronoi grid using the input triangulation and optional
        arguments such as the coordinate reference system, raster inputs, and idomain. It
        also initializes many attributes related to grid geometry, assigns centroids for
        grid cells, and prepares configurations for visual representations.

        :param tri: The triangulation instance used to create the Voronoi grid. Optional.
        :param crs: A string representing the coordinate reference system. Default
            is 'EPSG:2927'.
        :param rasters: A list or Path instance pointing to raster file(s) that can
            be associated with the Voronoi grid. Optional.
        :param name: Name of the Voronoi grid instance. Default is 'voronoi_grid'.
        :param qhull_options: Options for scipy's qhull algorithm, used for creating
            the Voronoi grid. Optional.
        :param idomain: A list used to define active and inactive areas in the grid.
            Optional. Should be a list of cell indices that are inactive
        :param idomain_path: A Path to an external input defining the idomain. Path
                should be to a geometry (ex. shapefile) that overlaps inactive cells. Optional.
        :param verts: A numpy array of coordinates for each vertex in the grid. Optional.
        :param iverts: A list of lists containing the indices of each vertex for each voronoi cell in the grid.
        :param xcyc: A numpy array of coordinates for the centroids of each cell in the grid. Optional.
        :param kwargs: Additional keyword arguments that may be passed to the parent
            class or for extended functionality.
        """

        print("VoronoiGrid initializing.")
        if tri:
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

        self.config = {
            'scrollZoom': True,
        }
        self.scatt_layout = {
            'height': 1000,
            'width': 1000,
            'dragmode': 'pan'
        }

        self.grid_centroid = self.get_grid_centroid()

        print('Voronoi grid initialized.')

    @property
    def vor_list(self):
        """list of voronoi cell geometries"""
        return self.gdf_vorPolys.geometry.to_list()

    @property
    def cell_list(self):
        """list of cell indices in the voronoi grid"""
        return self.gdf_vorPolys.index.to_list()

    @property
    def area_list(self):
        """list of cell areas in the voronoi grid"""
        if self._area_list is None:
            self._area_list = [cell.area for cell in self.vor_list]
        return self._area_list

    @property
    def idomain_path(self):
        return self._idomain_path

    @idomain_path.setter
    def idomain_path(self, value):
        if value is not None:
            assert isinstance(value, Path), 'idomain path must be a Path instance'
        self._idomain_path = value

    @property
    def idomain(self):
        return self._idomain

    @idomain.setter
    def idomain(self, value):
        if value is not None:
            assert isinstance(value, list), 'idomain must be a list'
            assert all(idx in self.cell_list for idx in value), 'idomain indices must be valid Voronoi grid cells'
        self._idomain = value

    @property
    def nlay(self):
        if self._nlay is None:
            self._nlay = len(self.gdf_topbtm.drop('geometry', axis=1).columns) - 1
        return self._nlay

    @property
    def gdf_topbtm(self):
        if self._gdf_topbtm is None and self.rasters is not None:
            self._gdf_topbtm = self.get_gdf_topbtm_multilyr(rasters=self.rasters)
        return self._gdf_topbtm

    @gdf_topbtm.setter
    def gdf_topbtm(self, value):
        self._gdf_topbtm = value

    """def get_disu_connectivity(self):

        df = self.gdf_vorPolys
        geoms = df.geometry
        centroids = df.geometry.centroid
        coords = np.array([(pt.x, pt.y) for pt in centroids])

        iac, ja, cl12, hwva = [], [], [], []

        print('getting connectivity properties (iac, ja, cl12, hwva, nja)')
        for i, neighbors in enumerate(self.adjacent_cells_idx):
            sorted_neighbors = sorted(neighbors)

            # iac is number of connections + 1 (self)
            iac.append(len(sorted_neighbors) + 1)

            # ja: self index followed by neighbors
            ja.extend([i] + sorted_neighbors)

            # cl12: 0 for self, then Euclidean distance to each neighbor
            dists = np.linalg.norm(coords[sorted_neighbors] - coords[i], axis=1)
            cl12.extend([0] + dists.tolist())

            # hwva: 0 for self, then shared edge length
            poly1 = geoms[i]
            prep_poly1 = prep(poly1)
            face_lengths = [
                poly1.intersection(geoms[j]).length if prep_poly1.intersects(geoms[j]) else 0
                for j in sorted_neighbors
            ]
            hwva.extend([0] + face_lengths)

        nja = sum(iac)
        self._iac = iac
        self._ja = ja
        self._cl12 = cl12
        self._hwva = hwva
        self._nja = nja

        return iac, ja, cl12, hwva, nja"""

    def get_disu_connectivity(self, *, use_representative_point: bool = False, tol: float = 0.0, validate: bool = True):
        """
        Build DISU connectivity vectors (iac, ja, cl12, hwva, nja) from planar polygons.

        Definitions (MF6-style):
          - iac[i] : number of entries for row i in ja/cl12/hwva (self + neighbors)
          - ja     : concatenation of each row's [i, neighbors...]
          - cl12   : 0 for self; for neighbors = |(c_j - c_i) · n_hat| (normal distance across shared edge)
          - hwva   : 0 for self; for neighbors = shared edge length
          - nja    : len(ja) == sum(iac)

        Args:
          use_representative_point : if True, use polygon.representative_point() for centers (robust for concave cells);
                                     otherwise use geometry.centroid
          tol : minimum face length to keep a connection (e.g., 1e-9 for numeric noise)
          validate : run consistency checks (recommended in development)

        Returns:
          (iac, ja, cl12, hwva, nja) as NumPy arrays (int64, int64, float64, float64, int)
        """
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
        if self._iac is None:
            self.get_disu_connectivity()
        return self._iac

    @property
    def ja(self):
        if self._ja is None:
            self.get_disu_connectivity()
        return self._ja

    @property
    def cl12(self):
        if self._cl12 is None:
            self.get_disu_connectivity()
        return self._cl12

    @property
    def hwva(self):
        if self._hwva is None:
            self.get_disu_connectivity()
        return self._hwva

    @property
    def nja(self):
        if self._nja is None:
            self.get_disu_connectivity()
        return self._nja

    @property
    def gdf_vorPolys(self):
        if self._gdf_vorPolys is None:
            self._gdf_vorPolys = self.get_gdf_vorPolys(crs=self.crs)
        return self._gdf_vorPolys

    @gdf_vorPolys.setter
    def gdf_vorPolys(self, value):
        self._gdf_vorPolys = value

    @property
    def adjacent_cells_idx(self):
        if self._adjacent_cells_idx is None:
            self._adjacent_cells_idx = self.find_adjacent_polygons(self.gdf_vorPolys)
        return self._adjacent_cells_idx

    def get_disu_props(self):
        self.iac
        self.ja
        self.nja
        self.cl12
        self.hwva

    def get_voronoi_polygons(self):
        """get polygons for each Voronoi cell, returns
        a list of Shapely polygons objects

        Args:
            verts (list): list of vertices
            iverts (list): list of lists of indices, each corresponding to a voronoi region

        Returns:
            list: list of Shapely polygon objects
        """
        return geometry_get_voronoi_polygons(self)

    def mapit(self, crs='EPSG:2927'):
        """maps voronoi grid based on x,y coords
        and a defined crs projection

        Args:
            crs (str, optional): string of coordinate reference system. Defaults to 'EPSG:2927'.

        Returns:
            geodataframe: returns a geodataframe explore method
        """
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
        return plotting_get_dash_selector(self)

    def show(self):
        return plotting_show(self)

    def map_nodes(self):
        return plotting_map_nodes(self)

    def get_vor_cells_as_series(
            self,
            overlapping_geometry: shp.Polygon | shp.Point | gpd.GeoSeries | Path = None,
            predicate: str = 'intersects',
            return_dict: bool = False,
            name_field: str = 'ExploName'
    ) -> pd.Series | dict:
        """
        Identify Voronoi cells matching a spatial predicate against provided geometries.

        :param overlapping_geometry: Shapely geometry, GeoSeries, GeoDataFrame, or a file path.
        :param predicate: Spatial query type (e.g., 'intersects', 'contains').
        :param return_dict: If True, returns a dictionary keyed by geometry name_field or index.
        :param name_field: Field name for geometry names when using GeoDataFrame or file input.
        :return: Series of intersecting cell indices or a dictionary mapping names to cell indices.
        """
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
            return_gdf: bool = False
    ) -> dict:
        """
        Provide a shapefile (locs) and get back the voronoi cells that contains them, by default.
        But you can change the predicate to search by something else, like 'intersects'.
        :param locs: shapefile of features to check against the vornoi grid
        :param crs: crs of method output. defaults to EPSG:2927
        :param predicate: None, “contains”, “contains_properly”, “covered_by”, “covers”,
        “crosses”, “intersects”, “intersects”, “touches”, “within”.
        :param loc_name_field: field in the loc shapefile that will be the dict key
        :return: dict of voronoi cells indices (values) that contain each location in locs (keys)
        """
        return selection_get_vor_cells_as_dict(
            self,
            locs=locs,
            crs=crs,
            predicate=predicate,
            loc_name_field=loc_name_field,
            return_gdf=return_gdf,
        )

    def show_selected_cells(
            self,
            cell_list: list = None,
            **kwargs

    ):
        """Method to show selected cells of the voronoi grid.
        Just provide a list of cell indices."""

        return plotting_show_selected_cells(self, cell_list=cell_list, **kwargs)

    def show_overlapping_geometry(self, shp_gpkg):
        """convenience method to show overlapping geometries just by providing a shapefile or geopackage"""
        return plotting_show_overlapping_geometry(self, shp_gpkg)

    def get_model_boundary_polygons(self) -> dict:
        """Returns a dict of the polygons that form the model domain boundary.
        The keys of the dict are voronoi cell indices"""
        return get_model_boundary_polygons(self.gdf_vorPolys)

    def get_grid_edge(self, idomain: list = None, idomain_path: Path = None, include_interiors: bool = True) -> list:
        """
        get edge cells for the voronoi grid. If a shapefile or geopackage of the idomain
        is provided, the returned edge cells are adjusted for the inactive cells
        :param include_interiors: if True, will include interior holes when returning grid edge cells
        :param idomain: provide list of cells NOT in the domain (inactive), instead of a geometry path - idomain_path
        :param idomain_path: provide to remove idomain cells from the returned grid edge cells
        :return:  list of grid edge cells
        """

        return get_grid_edge_cells(
            self,
            idomain=idomain,
            idomain_path=idomain_path,
            include_interiors=include_interiors,
        )

        """def is_edge(cell, idx):
            num_ja_cells = len(self.adjacent_cells_idx[idx])
            num_cell_faces = len(cell.geometry.exterior.coords) - 1
            edge_bool = True if num_cell_faces > num_ja_cells else False
            return edge_bool

        df = self.gdf_vorPolys
        edges = list(df[df.apply(lambda x: is_edge(x, x.name), axis=1)].index)

        if idomain_path:  # if idomain, determine new edge cells after removing idomain cells
            idomain = read_shp_gpkg(idomain_path)
            icells = self.get_vor_cells_as_series(idomain.geometry).to_list()
        elif idomain:
            icells = idomain
        if idomain_path or idomain:
            all_cells = self.gdf_vorPolys.copy()
            not_icells = [cell for cell in all_cells.index if cell not in icells]
            new_cells = all_cells.loc[not_icells, :]
            new_exterior = new_cells.union_all().exterior
            new_edge_cells = self.get_vor_cells_as_series(new_exterior)
            if include_interiors:
                interiors = new_cells.union_all().interiors
                interior_cells = []
                for interior in interiors:
                    interior_cells.append(self.get_vor_cells_as_series(interior))
                interior_cells.append(new_edge_cells)
                new_edge_cells = pd.concat(interior_cells)
            new_edge_cells = [cell for cell in new_edge_cells if cell not in icells]
            return new_edge_cells

        return edges"""

    def plot3d(self, z=None):
        return plotting_plot3d(self, z=z)

    def plot2d(self):
        return plotting_plot2d(self)

    def plottri(self):
        return plotting_plottri(self)

    def generate_grid_coordinates(self, grid_spacing: int) -> list:
        """Generates a grid of evenly spaced x- and y- coordinates
        based on an unstructured grid passed as 'self'. Returns 
        a list of x-coords and a list of y-coords. Grid will be
        rectangular, regardless of the input grid.  

        Args:
            grid_spacing (int): spacing between points on the generated grid

        Returns:
            list: two lists - one of x-coords and one of y-coords
        """
        return geometry_generate_grid_coordinates(self, grid_spacing)

    @property
    def centroids(self):
        """
        The `centroids` are internally calculated by invoking the `get_centroids()`
        method if they are not already computed.

        :return: The computed centroids of the voronoi grid cells.
        :rtype: Same type as returned by `get_centroids()` method.
        """
        if self._centroids is None:
            self._centroids = self.get_centroids()
        return self._centroids

    def get_centroids(self) -> tuple[list[float], list[float]]:
        """
        Returns coordinates of centroids of each Voronoi polygon.

        Returns:
            tuple: (list of x-coords, list of y-coords)
        """
        return geometry_get_centroids(self)

    def find_adjacent_polygons(self, gdf: gpd.GeoDataFrame) -> list:
        """Iterates though a GeoDataFrame of Voronoi polygons
        and finds the adjacent Voronoi polygon for each polygon. Returns
        a list of lists with the indices of all adjacent polygons for each
        polygon

        Args:
            gdf (gpd.GeoDataFrame): GeoDataFrame of Voronoi polygons

        Returns:
            list: list of lists, each containing the indices of adjacent
            polygons for each polygon
        """
        return geometry_find_adjacent_polygons(gdf)

    def find_adjacent_cells(self, cell_id):
        return geometry_find_adjacent_cells(self, cell_id)

    def calculate_distance(self, gdf, poly_idx1, poly_idx2):
        """Calculates the distance between the centroid of one polygon
        and the shared face of an adjacent polygon

        Args:
            gdf (GeoDataFrame): GeoDataFrame of Voronoi Grid
            poly_idx1 (int): index of polygon with centroid to measure from
            poly_idx2 (int): index of adjacent polygon

        Returns:
            float: returns distance from the centroid to the shared face
        """
        return geometry_calculate_distance(gdf, poly_idx1, poly_idx2)

    def shared_face_length(self, poly1, poly2):
        """Calculates the length of the shared face between
        two adjacent polygons in a Voronoi grid

        Args:
            poly1 (Polygon): Polygon object from GeoDataFrame
            poly2 (Polygon): Polygon object from GeoDataFrame

        Returns:
            float: length of shared face
        """
        return geometry_shared_face_length(poly1, poly2)

    def get_gdf_vorPolys(self, crs=None):
        return geometry_get_gdf_vor_polys(self, crs=crs)

    def get_gdf_topbtm(self, rasters: list, labels: list = None):
        """written for a one layer model with a top and bottom. More general funtion needed.
        Use get_gdf_topbtm_multilyr"""

        return self.get_gdf_topbtm_multilyr(rasters=rasters, labels=labels)

    def get_raster_vals_at_centroids(
            self,
            raster_files: Path | list[Path | int] = None,
            labels: str | list[str] = None
    ) -> gpd.GeoDataFrame | None:

        """Get centroid elevations of the Voronoi grid for each elevation raster or constant value.

        :param raster_files: Path or list of Paths to elevation raster(s) or constant values
        :param labels: str or list of labels for each raster/value
        :return: GeoDataFrame with elevation values at centroids"""

        return surface_get_raster_vals_at_centroids(self, raster_files=raster_files, labels=labels)

    """def get_raster_vals_at_centroids(
            self,
            raster_files: Path | list[Path | int | float] = None,
            labels: str | list[str] = None,
            max_search_radius: int = 256  # pixels; increase if your gaps are wider
    ) -> gpd.GeoDataFrame | None:
        
        if raster_files is None:
            print('No raster files provided')
            return None

        if isinstance(raster_files, (int, float, Path, str)):  # allow single arg
            raster_files = [raster_files]

        if labels is None:
            labels = list(range(len(raster_files)))
        if isinstance(labels, str):
            labels = [labels]

        if len(labels) != len(raster_files):
            raise ValueError("Labels length must match raster_files length")

        xs, ys = np.asarray(self.centroids[0]), np.asarray(self.centroids[1])
        centroid_vals = self.gdf_vorPolys.copy()

        # Separate numeric constants vs raster paths
        numeric_vals = {labels[i]: float(val)
                        for i, val in enumerate(raster_files)
                        if isinstance(val, (int, float))}
        raster_items = [(labels[i], Path(path))
                        for i, path in enumerate(raster_files)
                        if isinstance(path, (str, Path))]

        def nearest_valid_pixel_value(band, mask_valid, r0, c0, max_r):
    
            H, W = band.shape
            if mask_valid[r0, c0]:
                return float(band[r0, c0])

            # Expand search window
            for rad in range(1, max_r + 1):
                rmin = max(0, r0 - rad);
                rmax = min(H - 1, r0 + rad)
                cmin = max(0, c0 - rad);
                cmax = min(W - 1, c0 + rad)

                # window of candidate valid pixels
                win_mask = mask_valid[rmin:rmax + 1, cmin:cmax + 1]
                if not win_mask.any():
                    continue

                # coordinates of valid pixels within window
                rr, cc = np.nonzero(win_mask)
                rr = rr + rmin
                cc = cc + cmin

                # choose nearest by Euclidean distance
                d2 = (rr - r0) * (rr - r0) + (cc - c0) * (cc - c0)
                k = np.argmin(d2)
                return float(band[rr[k], cc[k]])

            # none found
            return np.nan

        # Process rasters
        for label, raster_path in raster_items:
            print(f'reading raster file {raster_path}')
            with rasterio.open(raster_path) as src:
                # --- 1) Reproject centroids if CRS differs ---
                xs_in, ys_in = xs, ys
                try:
                    crs_grid = getattr(self.gdf_vorPolys, "crs", None)
                    if crs_grid is not None and crs_grid != src.crs:
                        from pyproj import Transformer
                        tfm = Transformer.from_crs(crs_grid, src.crs, always_xy=True)
                        xs_in, ys_in = tfm.transform(xs_in, ys_in)
                except Exception:
                    # If transform fails, fall back to original coords (better to proceed)
                    pass

                # --- 2) Identify OOB centroids BEFORE clamping ---
                left, bottom, right, top = src.bounds
                eps = 1e-9
                oob = (xs_in < left) | (xs_in > right) | (ys_in < bottom) | (ys_in > top)

                # Clamp coords to bounds for indexing
                xs_c = np.clip(xs_in, left, right - eps)
                ys_c = np.clip(ys_in, bottom, top - eps)

                # --- 3) Map to pixel indices (ensure integers) ---
                rows, cols = rasterio.transform.rowcol(src.transform, xs_c, ys_c, op=np.floor)
                rows = np.clip(np.asarray(rows, dtype=np.int64), 0, src.height - 1)
                cols = np.clip(np.asarray(cols, dtype=np.int64), 0, src.width - 1)

                band = src.read(1)
                nodata = src.nodata
                vals = band[rows, cols].astype(float)

                # --- 4) Build validity masks ---
                # Explicit nodata mask (covers numeric nodata and NaN nodata)
                if nodata is None:
                    is_nodata = np.zeros_like(vals, dtype=bool)
                    band_valid_base = ~np.isnan(band) if np.issubdtype(band.dtype, np.floating) else np.ones_like(band,
                                                                                                                  dtype=bool)
                else:
                    if np.isnan(nodata):
                        is_nodata = np.isnan(vals)
                        band_valid_base = ~np.isnan(band)
                    else:
                        is_nodata = (vals == nodata)
                        band_valid_base = (band != nodata)

                # --- 5) Only treat ZERO as invalid for OOB points (edge behavior) ---
                zero_is_bad_for = oob  # boolean array, True only for OOB centroids
                is_zero = (vals == 0)
                need_zero_fix = is_zero & zero_is_bad_for

                # Anything nodata anywhere also needs fixing
                need_nodata_fix = is_nodata

                need_search = need_zero_fix | need_nodata_fix

                if np.any(need_search):
                    # For the search, define "valid" pixels:
                    # - Always exclude nodata
                    # - Exclude zeros ONLY when searching for an OOB point
                    # To avoid branching per-pixel, build two masks and pick per case.
                    band_valid_nozero = band_valid_base & (band != 0)
                    band_valid_allowzero = band_valid_base  # zeros allowed

                    r_bad = rows[need_search]
                    c_bad = cols[need_search]
                    repaired = np.empty_like(r_bad, dtype=float)

                    for i in range(r_bad.size):
                        # Choose which validity mask to use for this pixel
                        use_nozero = zero_is_bad_for[need_search][i]
                        mask_valid = band_valid_nozero if use_nozero else band_valid_allowzero

                        repaired[i] = nearest_valid_pixel_value(
                            band, mask_valid, int(r_bad[i]), int(c_bad[i]), max_search_radius
                        )

                    vals[need_search] = repaired

                # Normalize any remaining explicit nodata to NaN (numeric nodata only)
                if (nodata is not None) and (not np.isnan(nodata)):
                    vals = np.where(vals == nodata, np.nan, vals)

                centroid_vals[label] = vals

        # Assign static numeric values
        for label, val in numeric_vals.items():
            centroid_vals[label] = val

        # Ensure ordered columns: geometry first, then labels
        ordered_cols = ['geometry'] + labels
        centroid_vals = centroid_vals[ordered_cols]

        # Set geometry to centroids (points) for visualization/joins
        centroid_vals.geometry = centroid_vals.centroid

        return centroid_vals"""

    def get_gdf_topbtm_multilyr(self, rasters: list, labels: list = None):
        """
        Get a GeoDataFrame of top and bottom elevations for each model layer.

        :param rasters: list of rasters or static values ordered from top to bottom.
        :return: GeoDataFrame with layer elevations.
        """
        return surface_get_gdf_topbtm_multilyr(self, rasters=rasters, labels=labels)

    """def get_gdf_topbtm_multilyr(self, rasters: list):
        """
    """get a GeoDataFrame that has the elevations for the top of the model and the bottom of evey model layer for
    every voronoi cell in the model grid. For this to work as intended, the list of rasters need to be provided
    in order of top to bottom.
    :param rasters: this assumes that the list of raster is ordered from top to bottom
    :return: GeoDataFrame of top of model and bottom of every layer. The DataFrame columns are labeled starting at
    zero for the top of model, and all subsequent numbers correspond to elevations at the bottom of that numbered
    layer, for example label number '1' refers to the bottom of layer 1."""
    """
    labels = list(range(0, len(rasters)))
    gdf_topbtm = self.get_raster_vals_at_centroids(
        raster_files=rasters,
        labels=labels
    )
    #  TODO check difference between layers

    return gdf_topbtm

def get_raster_vals_at_centroids(
        self,
        raster_files: Path | list[Path | int],
        labels: str | list[str]
) -> gpd.GeoDataFrame:
    """
    """Get centroid elevations of the Voronoi grid for each elevation raster.
    Vectorized and optimized version.

    :param raster_files: Path or list of Paths to elevation raster(s)
    :param labels: str or list of labels for each raster
    :return: GeoDataFrame with elevation values at centroids"""
    """
    if isinstance(raster_files, Path):
        raster_files = [raster_files]
    if isinstance(labels, str):
        labels = [labels]

    if len(labels) != len(raster_files):
        raise ValueError("Labels length must match elevation files length")

    # separate provided elevations by single values or Paths for processing below
    is_int = [(i, obj) for i, obj in enumerate(raster_files) if isinstance(obj, int | float)]
    non_ints = [(i, obj) for i, obj in enumerate(raster_files) if isinstance(obj, Path)]
    if len(is_int) > 0:
        int_idx, int_layers = zip(*is_int)
    else:
        int_layers = []
    if len(non_ints) > 0:
        raster_idx, raster_layers = zip(*non_ints)
    else:
        raster_layers = []

    # Open raster files and store arrays and metadata
    srcs = [rasterio.open(path) for path in raster_layers]
    arrays = [src.read(1) for src in srcs]
    transforms = [src.transform for src in srcs]
    bounds = [src.bounds for src in srcs]

    # Get centroid coordinates as NumPy arrays
    centroids = self.gdf_vorPolys.geometry.centroid
    xs, ys = centroids.x.to_numpy(), centroids.y.to_numpy()

    centroid_vals = self.gdf_vorPolys.copy()

    # Create result GeoDataFrame and populate elevation values for RASTER layers
    if len(raster_layers) > 0:

        raster_labels = [labels[i] for i in raster_idx]
        for i, (label, arr, tfm, bds) in enumerate(zip(raster_labels, arrays, transforms, bounds)):
            # Convert x,y to raster row/col indices (float by default)
            rows, cols = map(np.array, rasterio.transform.rowcol(tfm, xs, ys, op=np.floor))

            # Identify valid points that fall within the raster bounds
            valid = (
                    (rows >= 0) & (cols >= 0) &
                    (rows < arr.shape[0]) & (cols < arr.shape[1])
            )

            # Initialize elevation values with NaN, then assign only valid entries
            values = np.full(xs.shape, np.nan)
            valid_rows = rows[valid].astype(int)
            valid_cols = cols[valid].astype(int)
            values[valid] = arr[valid_rows, valid_cols]

            # Add to the GeoDataFrame under the given label
            centroid_vals[label] = values.astype(float)

    # Close all raster files
    for src in srcs:
        src.close()

    # for layers where only int or float elevations provided, make all centroids that elev
    if len(int_layers) > 0:
        int_labels = [labels[i] for i in int_idx]
        for i, (label, val) in enumerate(zip(int_labels, int_layers)):
            centroid_vals[label] = float(val)

    # make sure layer columns are in correct order
    labels = ['geometry'] + labels
    centroid_vals = centroid_vals[labels]
    centroid_vals.geometry = centroid_vals.centroid

    return centroid_vals"""
    """def get_centroid_elevations(
            self,
            elevations_files: Path | list[Path],
            labels: str | list[str]
    ) -> gpd.GeoDataFrame:
    
        if isinstance(elevations_files, Path):
            elevations_files = [elevations_files]
        if isinstance(labels, str):
            labels = [labels]

        n_rasters = len(elevations_files)
        if len(labels) != n_rasters:
            raise ValueError("Length of labels must match number of elevation files")

        # Open rasters and read metadata
        srcs = [rasterio.open(path) for path in elevations_files]
        elevations = [src.read(1) for src in srcs]
        shapes = [arr.shape for arr in elevations]  # (rows, cols)

        def get_elevations(x, y):
            elevs = []
            for i in range(n_rasters):
                try:
                    row, col = srcs[i].index(x, y)
                except ValueError:
                    # x, y is outside raster bounds
                    elevs.append(np.nan)
                    continue

                if 0 <= row < shapes[i][0] and 0 <= col < shapes[i][1]:
                    elev = elevations[i][row, col]
                else:
                    elev = np.nan  # out of bounds
                elevs.append(elev)
            return tuple(elevs)

        # Copy centroid geometry
        centroids_gdf = self.gdf_vorPolys.copy()
        centroids_gdf.geometry = centroids_gdf.centroid

        # Apply elevation sampling per raster
        elev_data = centroids_gdf.geometry.apply(lambda pt: get_elevations(pt.x, pt.y))

        for i, label in enumerate(labels):
            centroids_gdf[label] = elev_data.apply(lambda t: t[i])

        for src in srcs:
            src.close()

        return centroids_gdf"""

    def get_cell_areas(self):
        return surface_get_cell_areas(self)

    def get_origin_xy(self):
        """Get the x,y coordinates for the origin of the 
        model grid.

        Returns:
            tuple: return a tuple of the form - (x,y) 
        """

        return surface_get_origin_xy(self)

    @property
    def gdf_latlon(self):
        return geometry_get_gdf_latlon(self)

    @property
    def latlon(self):
        return geometry_get_latlon(self)

    def get_grid_centroid(self):
        """Gets the Shapley representation of the centroid of the defined voronoi grid

        Returns:
            shp: returns Shapely representation of the grid centroid
        """
        return geometry_get_grid_centroid(self)

    def get_overlapping_area(self, shp_gpkg=None, cell_list=None):
        """
        simple method to get the voronoi cell area of an overlapping geometry. Can provide
        a shapefile or geopackage of the geometry. Function will determine what voronoi cells
        the geometry overlaps and then calculate the area of those voronoi cells. Alternatively,
        if you know the cell ids, you can provide a list of cell indices.
        :param shp_gpkg: shapefile or geopackage of geometry to check
        :param cell_list: list of cell ids. If both shp_gpkg and cell_list are provided, shp_gpkg takes precedence
        :return:
        """
        return geometry_get_overlapping_area(self, shp_gpkg=shp_gpkg, cell_list=cell_list)

    def get_vor_idx_from_geometry(
            self,
            shp_to_query: shp = None,
            gdf_to_query: gpd.GeoDataFrame = None,
            crs: str = "EPSG:2927",
            crs_latlon: str = "EPSG:4326",
            name_col: str = None,
            predicate: str = "intersects"
    ) -> dict:
        """Method to do a spatial query to determine the voronoi cells that intersect
        the given geometries. Geometries should be a shapefile of points or polygons. This
        function returns a dictionary with keys consisting of shapefile indices or names in
        given name_col field and values consisting of the indices of the intersecting
        voronoi cells for each key. The predicate corresponds to the sindex.query function
        in the GeoPandas package.
        
        Args:
            shp_to_query: shapefile of points or polygons to query
            gdf_to_query: GeoDataFrame with geometry to query; this will take precedence
            over shp_to_query if it is passed to the function
            crs (string): coordinate reference system for the shapefile
            crs_latlon (string): crs for lat/lon, EPSG:4326. Shouldn't need to change this
            name_col (string): string corresponding to the field name in the shapefile to be used
            as for the keys in the return dict.
            predicate (string): method for query, defaults to intersects, but can use any predicate allowed
            by GeoPandas sindex.query
        """

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
            predicate: str = "intersects"
    ) -> list:
        """Method to do a spatial query to determine the voronoi cells that intersect
        the single geometry in a GeoDataFrame. Use get_vor_idx_from_geometry for multiple
        geometries. This function returns a list consisting indices of the intersecting
        voronoi cells for the given geometry. The predicate corresponds to the sindex.query function
        in the GeoPandas package.
        
        Args:
            gdf_to_query: GeoDataFrame with geometry to query
            idx: GeoDataFrame index with geometry to query
            predicate (string): method for query, defaults to intersects, but can use any predicate allowed
            by GeoPandas sindex.query
        """
        return geometry_get_vor_idx_from_geometry_idx(
            self,
            gdf_to_query=gdf_to_query,
            idx=idx,
            predicate=predicate,
        )

    def set_k_vor(
            self,
            k_dict: dict = None,
            k_default=100
    ) -> list:
        """Method to set the hydraulic conductivity for all voronoi cells
        in the model grid. k_dict should be generated using the get_vor_idx_from_geometry
        method. The k_default is assigned to all cells prior to k_dict, in case some
        cells are not in k_dict."""
        return geometry_set_k_vor(self, k_dict=k_dict, k_default=k_default)

    """def get_raster_from_strike_dip(
            self,
            strike: int,
            dip: int,
            known_point: tuple,
            pixel_size: int = 1,
            output_filename: Path = Path.cwd().joinpath('raster.tif')
    ):
        """
    """

        hull_bounds = shp.MultiPolygon(self.gdf_vorPolys.geometry.to_list()).convex_hull.bounds
        # Unpack the known point and hull bounds
        known_x, known_y, known_elevation = known_point
        min_x, min_y, max_x, max_y = hull_bounds
        # Calculate the dimensions of the raster
        width = int((max_x - min_x) / pixel_size)
        height = int((max_y - min_y) / pixel_size)
        dip_rad = np.radians(90 - dip)

        # Calculate the normal vector to the plane
        normal = self.get_normal_from_strike_and_dip(strike, dip)
        # Create an affine transform for the raster
        transform = from_origin(min_x, max_y, pixel_size, pixel_size)
        # Initialize the raster array
        # elevation_data = np.zeros((height, width), dtype=rasterio.float32)

        # Calculate the elevation values using vectorized operations
        print('getting elevations for raster')

        # Create a meshgrid of x and y coordinates
        x_coords = min_x + np.arange(width) * pixel_size
        y_coords = max_y - np.arange(height) * pixel_size
        x_grid, y_grid = np.meshgrid(x_coords, y_coords)

        # Calculate the position vectors of all pixels relative to the known point
        point_vectors = np.stack([x_grid - known_x, y_grid - known_y, np.zeros_like(x_grid)], axis=-1)

        # Calculate the dot product for all points
        norm_normal = np.linalg.norm(normal)
        distance_along_normal = np.dot(point_vectors, normal) / norm_normal

        # Calculate the elevation for all points
        elevation_data = known_elevation - distance_along_normal * np.cos(dip_rad)
        print('got raster from strike and dip')

        # Write the raster to a file
        with rasterio.open(
                output_filename,
                'w',
                driver='GTiff',
                height=height,
                width=width,
                count=1,
                dtype=rasterio.float32,
                crs=self.crs,
                transform=transform,
        ) as dst:
            dst.write(elevation_data, 1)

        centroids_gdf = self.get_centroid_elevations([output_filename], ['elev'])

        return centroids_gdf"""

    def get_raster_from_strike_dip(
            self,
            strike: int,
            dip: int,
            known_point: tuple,
            pixel_size: int = 1,
            output_filename: Path = Path.cwd().joinpath('raster.tif')
    ):
        """
        Generate a raster file representing elevations of a sloping plane.

        Parameters:
        - strike: The strike of the plane in degrees, measured from north.
        - dip: The dip of the plane in degrees, measured from the horizontal.
        - known_point: A tuple (x, y, elevation) for a known point on the plane.
        - pixel_size: The size of each pixel in spatial units.
        - output_filename: The filename for the output raster.
        """

        return surface_get_raster_from_strike_dip(
            self,
            strike=strike,
            dip=dip,
            known_point=known_point,
            pixel_size=pixel_size,
            output_filename=output_filename,
        )

    def to_shapefile(self, filepath: str | Path = 'vor_shp.shp'):
        return self.gdf_vorPolys.to_file(filepath)

    def get_domain(self):
        return surface_get_domain(self)

    @staticmethod
    def get_normal_from_strike_and_dip(strike: int, dip: int) -> np.array:
        """
        Returns the normal vector of a plane based on strike and dip
        :param strike: strike in degrees
        :param dip: dip in degrees
        :return: normal vector
        """
        return surface_get_normal_from_strike_and_dip(strike, dip)

    @staticmethod
    def generate_grid_around_point(center_point: shp.Point, spacing: float, size: int, crs: str) -> gpd.GeoSeries:
        """Generate a GeoPandas GeoSeries of points in a grid pattern.
        
        Args:
            center_point (shapely.geometry.Point): The center point of the grid.
            spacing (float): The spacing between grid points.
            size (int): The number of points in one dimension of the grid.
            crs (str): Set the crs of the generated points

        Returns:
            geopandas.GeoSeries: The GeoSeries of points.
        """
        return geometry_generate_grid_around_point(center_point, spacing, size, crs)

    @staticmethod
    def generate_grid_polygons(center_point, spacing, size, gap):
        """Generate a MultiPolygon object with polygons in a grid pattern with a set spacing.
        
        Args:
            center_point (shapely.geometry.Point): The center point of the grid.
            spacing (float): The spacing between grid points.
            size (int): The number of points in one dimension of the grid.
            gap (float): The gap between polygons.

        Returns:
            shapely.geometry.MultiPolygon: The MultiPolygon of grid squares.
        """
        return geometry_generate_grid_polygons(center_point, spacing, size, gap)

    @staticmethod
    def voronoi_refine_by_point(point: shp.Point, spacing: int, tri: Triangle) -> gpd.GeoDataFrame:
        return geometry_voronoi_refine_by_point(point, spacing, tri)

    def reconcile_surfaces(self, df: pd.DataFrame = None, min_sep=0.1, trigger_sep=1, which='bottom'):
        """
        helper to iterate through surface elevations and check for layers that are above the overlying
        layer, then adjust so they don't overlap.
        :param which: 'bottom' or 'top'. If bottom, will adjust bottom layer to maintain min_sep. Same with top.
        :param trigger_sep: trigger separation, if separation is less than this the layers will be adjusted
        :param df: dataframe of surface elevations at each voronoi cell, column names are the surface names
        :param min_sep: surfaces that are too high will be reduced below the overlying surface by this minimum separation
        :return: new dataframe with adjusted surface elevations
        """
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
            reconcile: bool = True
    ):
        """
        method to adjust the elevations of provided cells by given adjustment amount
        :param cell_ids: list of cell ids to adjust
        :param adjustment: amount to adjust the cells by
        :param df: DataFrame of layer elevations. Index are cell ids, columns are layers
        :param layer: which layer to adjust, 0 = top of 1st layer, 1 = 1st layer botom, 2 = 2nd layer bottom, etc.
        :param reconcile: boolean value to indicate if the adjusted surfaces df should be sent to self.reconciled_surfaces
        :return: df
        """
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
            shp: Path = None, buffer=1,
            layer_bottom_name=1,
            min_sep=None
    ) -> pd.DataFrame:
        """
        method to easily adjust the bottoms of certain voronoi cells so that the bottoms are not higher than any of
        the tops of the adjacent cells. In order to have a continually overlapping layer of cells horiontally. This
        is mostly an issue for steeply sloping surfaces with thinner layer thicknesses.
        :param elev_df: DataFrame with cell elevations to adjust. Defaults to the return df of self.reconcile_surfaces()
        :param shp: Path of shapefile that identifies cells to be adjusted
        :param buffer: how much extra below the lowest adjacent cell top to lower each cell, defaults to 1
        :param layer_bottom_name: name of the layer we are adjusting, defaults to layer 1
        :return: returns a new dataframe of reconciled surfaces
        """
        return geometry_adjust_top_btm_overlaps(
            self,
            elev_df=elev_df,
            shp_path=shp,
            buffer=buffer,
            layer_bottom_name=layer_bottom_name,
            min_sep=min_sep,
        )

    def cross_section(self, line: shp.LineString | Path):
        """return a GridSetion object from provided line. can plot with the .show method"""
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

        verts, iverts, xcyc = get_griddata_from_disu(disu_path)
        vor = cls(verts=verts, iverts=iverts, xcyc=xcyc, crs=crs, name=name,
                  rasters=rasters, idomain=idomain, idomain_path=idomain_path,
                  qhull_options=qhull_options)

        return vor
