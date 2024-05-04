import numpy as np
import plotly.graph_objects as go
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point
import pandas as pd
import geopandas as gpd
from flopy.utils.voronoi import VoronoiGrid, tri2vor
from flopy.utils.triangle import Triangle as Triangle
import rasterio
from rasterio.transform import from_origin
import shapely as shp
from shapely.geometry import LineString
from shapely.geometry import Polygon, MultiLineString, Point
import json
from pathlib import Path
from . import mf2Dplots
from figs import create_hover


def flatten(l):
    return [item for sublist in l for item in sublist]

class TriangleGrid(Triangle):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def add_circle(
            self,
            radius: float = 100.0,
            center_coords: tuple = (0, 0),
            polygon_to_add=None,
            point_to_add=None,
            point_region_size_max=1,
            return_only: bool = False,
            radians_step: int = 0.1
    ):
        """Create a circular grid. Center coords must be a single (x,y) tuple."""
        theta = np.arange(0.0, 2 * np.pi, radians_step)
        x = radius * np.cos(theta) + center_coords[0]
        y = radius * np.sin(theta) + center_coords[1]
        circle_poly = [(x, y) for x, y in zip(x, y)]
        if not return_only:
            self.add_polygon(circle_poly)

        if polygon_to_add:
            polygon_to_add = Polygon(shell=polygon_to_add)
            self.add_polygon(polygon=polygon_to_add)
        if point_to_add:
            self.add_region(point=point_to_add, maximum_area=point_region_size_max)

        return circle_poly

    def add_rectangle(
            self,
            x_dist=100,
            y_dist=100,
            origin=(0, 0),
            return_only = False
    ):
        x_min, y_min = origin[0], origin[1]
        x_max, y_max = x_min + x_dist, y_min + y_dist
        polygon_coords = ((x_min, y_min), (x_min, y_max),
                          (x_max, y_max), (x_max, y_min))
        polygon = shp.Polygon(polygon_coords)
        if not return_only:
            self.add_polygon(polygon)

        return polygon

    def add_regions(self, points, attributes=None, maximum_areas=None):
        """
        Add points that will become regions with maximum areas, if
        specified.

        Parameters
        ----------
        points : tuples
            (x, y)...

        attributes : integer or float
            integer values assigned to output elements

        maximum_areas : float
            maximum area of elements in regions

        Returns
        -------
        None

        """
        if attributes is None:
            attributes = [0 for _ in points]
        if maximum_areas is None:
            maximum_areas = [None for _ in points]

        regions = zip(points, attributes, maximum_areas)
        for region in regions:
            self._regions.append(region)
        return


class VoronoiGridPlus(VoronoiGrid):
    def __init__(self,
                 tri: Triangle,
                 crs: str = 'EPSG:2927',
                 rasters: list | Path = None,
                 name: str = 'voronoi_grid'
                 ):

        super().__init__(tri)
        self.vor = VoronoiGrid(tri=tri)
        print("VoronoiGrid initialized.")
        self.verts = self.vor.verts
        self.iverts = self.vor.iverts
        self.tri = tri
        self.crs_latlon = "EPSG:4326"
        self.rasters = rasters
        self.crs = crs
        self.name = name
        self.points = self.tri.verts
        self.x_coords_by_node = []
        self.y_coords_by_node = []
        self._gdf_vorPolys = None
        self._iac = None
        self._nja = None
        self._gdf_topbtm = None
        self._adjacent_cells_idx = None
        self._ja, self._cl12, self._hwva = None, None, None
        self.centroids_x, self.centroids_y = self.get_centroids()

        print('Getting SciPy voronoi grid')
        self.vor_scipy = Voronoi(self.tri.verts)
        self.scipy_points = self.vor_scipy.points
        self.scipy_ridge_points = self.vor_scipy.ridge_points
        self.scipy_ridge_vertices = self.vor_scipy.ridge_vertices
        self.scipy_verts = self.vor_scipy.vertices
        self.scipy_regions = self.vor_scipy.regions
        print('Got SciPy voronoi grid')

        self.config = {
            'scrollZoom': True,
        }
        self.scatt_layout = {
            'height': 1000,
            'width': 1000,
            'dragmode': 'pan'
        }

        self.x_vor_regions = []
        self.y_vor_regions = []
        self.x_vor = self.verts[:, 0]
        self.y_vor = self.verts[:, 1]
        for region in self.iverts:
            x_verts_in_region = [self.verts[region[i]][0] for i in range(len(region))] + [None]
            y_verts_in_region = [self.verts[region[i]][1] for i in range(len(region))] + [None]
            self.x_vor_regions.append(x_verts_in_region)
            self.y_vor_regions.append(y_verts_in_region)

        # x-coords and y-coords for voronoi regions, separated by None values
        self.x_vor_regions = flatten(self.x_vor_regions)
        self.y_vor_regions = flatten(self.y_vor_regions)

        # i,j,k for each triangle in the triangulation (argument = tri)
        self.i = [tri_idx[0] for tri_idx in self.iverts]
        self.j = [tri_idx[1] for tri_idx in self.iverts]
        self.k = [tri_idx[2] for tri_idx in self.iverts]

        self.gdf_latlon = None
        self.latslons = None
        if self.crs is not None:
            print('getting lats and lons')
            self.get_latslons()  # generate json of grid and save to self
            print('got lats and lons')
        self.grid_centroid = self.get_grid_centroid()

        self.vor_list = self.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]

    @property
    def gdf_topbtm(self):
        if self._gdf_topbtm is None and self.rasters is not None:
            self._gdf_topbtm = self.get_gdf_topbtm(rasters=self.rasters)
        return self._gdf_topbtm

    @gdf_topbtm.setter
    def gdf_topbtm(self, value):
        self._gdf_topbtm = value

    @property
    def ja(self):
        if self._ja is None:
            ja_cl12_hwva = self.get_ja_cl12_hwva()
            self._ja = ja_cl12_hwva[0]
            self._cl12 = ja_cl12_hwva[1]
            self._hwva = ja_cl12_hwva[2]
        return self._ja

    @property
    def cl12(self):
        if self._cl12 is None:
            ja_cl12_hwva = self.get_ja_cl12_hwva()
            self._ja = ja_cl12_hwva[0]
            self._cl12 = ja_cl12_hwva[1]
            self._hwva = ja_cl12_hwva[2]
        return self._cl12

    @property
    def hwva(self):
        if self._hwva is None:
            ja_cl12_hwva = self.get_ja_cl12_hwva()
            self._ja = ja_cl12_hwva[0]
            self._cl12 = ja_cl12_hwva[1]
            self._hwva = ja_cl12_hwva[2]
        return self._hwva

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

    @property
    def iac(self):
        if self._iac is None:
            self._iac = self.get_iac()
        return self._iac

    @property
    def nja(self):
        if self._nja is None:
            self._nja = self.get_nja()
        return self._nja

    def get_voronoi_polygons(self):
        """get polygons for each Voronoi cell, returns
        a list of Shapely polygons objects

        Args:
            verts (list): list of vertices
            iverts (list): list of lists of indices, each corresponding to a voronoi region

        Returns:
            list: list of Shapely polygon objects
        """
        polygons = []
        for region in self.iverts:
            vertices = [self.verts[j] for j in region]
            polygon = Polygon(vertices)
            polygons.append(polygon)
        return polygons

    def get_iac(self):

        # Initialize list of adjacent cell counts
        adjacent_counts = [0] * len(self.iverts)

        # Loop through each ridge and increment adjacent cell counts
        for ridge in self.scipy_ridge_points:
            if -1 not in ridge:
                i, j = ridge
                adjacent_counts[i] += 1
                adjacent_counts[j] += 1

        return adjacent_counts

    def mapit(self, crs='EPSG:2927'):
        """maps voronoi grid based on x,y coords
        and a defined crs projection

        Args:
            crs (str, optional): string of coordinate reference system. Defaults to 'EPSG:2927'.

        Returns:
            geodataframe: returns a geodataframe explore method
        """
        poly = self.get_voronoi_polygons()

        return gpd.GeoDataFrame(geometry=poly, crs=crs).explore()

    def plot_choropleth(self, zmin=None, zmax=None, zoom=18):
        """Plot choropleth of Voronoi grid.

            Args:

            """
        if zmax is None:
            zmax = len(self.latslons['features'])
        if zmin is None:
            zmin = 0

        fig_mbox = mf2Dplots.ChoroplethPlot(vor=self, zoom=zoom)
        custom_data, hover_template = create_hover(
            {
                'Cell No.': self.cell_list,
                'Area': self.area_list,
                'x': self.x_list,
                'y': self.y_list
            }
        )
        fig_mbox.add_choroplethmapbox(
            geojson=self.latslons,
            featureidkey="id",
            locations=self.gdf_latlon.cell.to_list(),
            z=self.gdf_latlon.cell.to_list(),
            customdata=custom_data,
            colorscale="earth",
            zmax=zmax,
            zmin=zmin,
            hovertemplate=hover_template
        )

        return fig_mbox

    def map_nodes(self):

        latlonselect = self.gdf_vorPolys.to_crs('EPSG:4326')
        latlonselect['cellidx'] = latlonselect.index.astype(str)
        geojsonselect = json.loads(latlonselect['geometry'].to_json())
        centroid_grid = self.grid_centroid
        fig_sel = go.Figure(go.Choroplethmapbox(
            geojson=geojsonselect,
            locations=latlonselect['cellidx'].to_list(),
            featureidkey='id',
            z=latlonselect['cellidx'].to_list(),
            colorscale='earth'
        ))

        fig_sel.update_layout(
            margin={"r": 0, "t": 20, "l": 0, "b": 0},
            mapbox_style="carto-positron",
            mapbox_zoom=15,
            mapbox_center={"lat": centroid_grid.y, "lon": centroid_grid.x},
        )

        return fig_sel

    def get_overlapping_grid_cells(
            self,
            overlapping_geometry: shp.Polygon | shp.Point | gpd.GeoSeries = None,
            predicate: str = 'covered_by'
    ):
        """
        Compares given geometries to the voronoi polygon geometries and returns
        cells that math the given predicate (as defined by GeoPandas spatial query).
        :param overlapping_geometry: Geometry to compare to voronoi gird, can be a
        shapely polygon, point, or a GeoSeries
        :param predicate: options defined by GeoPandas spatial index query
        :return: Pandas Series of intersecting cells
        """
        if isinstance(overlapping_geometry, (shp.Polygon, shp.Point)):
            print('ok')

        elif isinstance(overlapping_geometry, gpd.GeoSeries):
            geocount = overlapping_geometry.count()

            if isinstance(overlapping_geometry.array, gpd.array.GeometryArray) and geocount != 0:
                print(f'You gave me {geocount} geometries to check')
            else:
                print('The GeoSeries is empty!')
                return

        else:
            print('Wrong data types')
            return

        print('well done')

        geometries = gpd.GeoSeries(overlapping_geometry)
        intersecting_cells = geometries.array.sindex.query(
            self.gdf_vorPolys["geometry"],
            predicate=predicate,
        )[0]

        return pd.Series(intersecting_cells)

    def show_selected_cells(
            self,
            cell_list: list = None
    ):
        """Method to show selected cells of the voronoi grid.
        Just provide a list of cell indices."""

        choro = self.plot_choropleth()
        choro.data[0].selectedpoints = (tuple(cell_list))

        return go.Figure(choro)

    def get_model_boundary_polygons(self) -> dict:
        """Returns a dict of the polygons that form the model domain boundary.
        The keys of the dict are voronoi cell indices"""

        grid = shp.MultiPolygon(self.gdf_vorPolys.geometry.to_list())
        polygons = self.gdf_vorPolys.geometry.to_list()
        convex_hull = grid.convex_hull
        convex_hull_boundary = convex_hull.boundary

        # This will hold the polygons that are on the boundary of the convex hull
        boundary_polygons = []
        boundary_polygons_idx = []

        for idx, polygon in enumerate(polygons):
            # Check if the polygon boundary intersects with the convex hull boundary
            if convex_hull_boundary.dwithin(polygon, 0.01):
                boundary_polygons.append(polygon)
                boundary_polygons_idx.append(idx)
        boundary_polygon_dict = dict(zip(boundary_polygons_idx, boundary_polygons))
        return boundary_polygon_dict

    def plot3d(self, z=None) -> go.Figure:
        if z == None:
            z = [0 for x in self.x_vor_regions]
        else:
            z = z
        fig3d = go.Figure(
            data=go.Scatter3d(
                x=self.x_vor_regions,
                y=self.y_vor_regions,
                z=z,
                opacity=0.5, mode='lines',
                line_color='black',
                # marker_size=2,
            ),
            layout={
                'height': 1000,
            },
        )
        fig3d.add_trace(
            go.Mesh3d(
                x=self.x_vor,
                y=self.y_vor,
                z=[0 for z in range(len(self.x_vor))],
                i=self.i,
                j=self.j,
                k=self.k,
                colorscale='Viridis',
                intensity=self.x_vor

            )
        )
        fig3d.update_layout(title='Voronoi Diagram',
                            scene=dict(
                                xaxis=dict(title='X'),
                                yaxis=dict(title='Y'),
                                zaxis=dict(title='')
                            )
                            )
        return fig3d

    def plot2d(self) -> go.Figure:

        fig2d = go.Figure(layout=self.scatt_layout)
        for cell in range(len(self.x_coords_by_node)):
            fig2d.add_scattergl(
                x=self.x_coords_by_node[cell],
                y=self.y_coords_by_node[cell],
                opacity=1,
                mode='lines',
                line_color='black',
                line_width=1,
                # fill='toself',
                # fillcolor='gray',
                # hoveron='points+fills'
            )

        self.fig2d = fig2d

        return fig2d

    def plottri(self) -> go.Figure:
        """plot the triangulated mesh generated by the Triangle program

        Returns:
            plotly fig: plotly figure of triangulated mesh
        """
        df_tricells = pd.DataFrame(self.tri.get_cell2d(),
                                   columns=[
                                       'triidx',
                                       'x',
                                       'y',
                                       'numverts',
                                       'idx1',
                                       'idx2',
                                       'idx3'
                                   ])
        df_triverts = pd.DataFrame(self.tri.verts, columns=['x', 'y'])
        ilist = df_tricells['idx1'].to_list()
        jlist = df_tricells['idx2'].to_list()
        klist = df_tricells['idx3'].to_list()
        ijkzip = zip(ilist, jlist, klist)
        xtriverts = df_triverts['x'].to_list()
        ytriverts = df_triverts['y'].to_list()
        xtricells = []
        ytricells = []
        for v in list(ijkzip):
            i, j, k = v[0], v[1], v[2]
            xtricur = [xtriverts[i]] + [xtriverts[j]] + [xtriverts[k]] + [None]
            xtricells = xtricells + xtricur
            ytricur = [ytriverts[i]] + [ytriverts[j]] + [ytriverts[k]] + [None]
            ytricells = ytricells + ytricur
        trifig2d = go.Figure(go.Scattergl(x=xtricells,
                                          y=ytricells,
                                          mode='lines',
                                          line_color='black',
                                          line_width=1, ),
                             layout=self.scatt_layout,
                             )
        return trifig2d.show(config=self.config)

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
        x_min = self.x_vor.min()
        x_max = self.x_vor.max()
        y_min = self.y_vor.min()
        y_max = self.y_vor.max()
        x_coords = []
        y_coords = []
        for x in range(int(x_min), int(x_max) + grid_spacing, grid_spacing):
            for y in range(int(y_min), int(y_max) + grid_spacing, grid_spacing):
                x_coords.append(x)
                y_coords.append(y)
        return x_coords, y_coords

    def get_centroids(self) -> list:
        """Returns coordinates of centroids of each Voronoi
        polygon.

        Returns:
            list: lists of x-coords and y-coords
        """
        vor_xverts_per_cells = []
        vor_yverts_per_cells = []
        vor_centroids_x = []
        vor_centroids_y = []

        for cell in self.vor.iverts:
            thiscellx = []
            thiscelly = []
            for idx in cell:
                thiscellx += [self.vor.verts[idx][0]]
                thiscelly += [self.vor.verts[idx][1]]
                thiscell_centroidx = np.array(thiscellx).mean()
                thiscell_centroidy = np.array(thiscelly).mean()
            vor_xverts_per_cells.append(thiscellx)
            vor_yverts_per_cells.append(thiscelly)
            vor_centroids_x.append(thiscell_centroidx)
            vor_centroids_y.append(thiscell_centroidy)
        return vor_centroids_x, vor_centroids_y

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

        # Create a list to hold the adjacent polygon indices for each polygon
        adjacent = [[] for _ in range(len(gdf))]

        # Create a spatial index to speed up the search for adjacent polygons
        gdf_sindex = gdf.sindex

        # Iterate over each polygon in the GeoDataFrame
        for i, poly in gdf.geometry.items():
            # Get the indices of all polygons that intersect the bounding box of the current polygon
            candidates = list(gdf_sindex.intersection(poly.bounds))

            # Iterate over the candidate polygons and check for adjacency
            for j in candidates:
                # Skip the current polygon
                if i == j:
                    continue

                # Check if the current polygon shares an edge with the candidate polygon
                shared_edge = poly.intersection(gdf.iloc[j].geometry)
                if isinstance(shared_edge, LineString):
                    # Add the index of the adjacent polygon to the list for the current polygon
                    adjacent[i].append(j)

        return adjacent

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
        # Get the Polygon objects for the two specified polygons
        poly1 = gdf.iloc[poly_idx1].geometry
        poly2 = gdf.iloc[poly_idx2].geometry

        # Calculate the shared edge between the two polygons
        shared_edge = poly1.intersection(poly2)

        # Calculate the centroid of the first polygon
        centroid = poly1.centroid

        # Calculate the distance between the centroid and the shared edge
        distance = centroid.distance(shared_edge)

        return distance

    def shared_face_length(self, poly1, poly2):
        """Calculates the length of the shared face between
        two adjacent polygons in a Voronoi grid

        Args:
            poly1 (Polygon): Polygon object from GeoDataFrame
            poly2 (Polygon): Polygon object from GeoDataFrame

        Returns:
            float: length of shared face
        """
        # Calculate the intersection of the two polygons
        intersection = poly1.intersection(poly2)
        # If the intersection is not a LineString, return 0
        if not isinstance(intersection, LineString):
            return 0

        # Calculate the length of the intersection
        length = intersection.length
        return length

    def get_gdf_vorPolys(self, crs=None):
        vertices_by_cells = []
        xvertices_by_cells = []
        yvertices_by_cells = []
        polys_ListbyShapely = []

        """find lists of cell vertices's x,y-coords, 
        x-coords, and y-coords arranged by cell"""
        for cell in self.iverts:
            thiscell = []
            thiscellx = []
            thiscelly = []
            for vidx in cell:
                thisvertx = self.verts[vidx][0].tolist()
                thisverty = self.verts[vidx][1].tolist()
                thisvert = self.verts[vidx].tolist()
                thiscellx += [thisvertx]
                thiscelly += [thisverty]
                thiscell += [thisvert]
            xvertices_by_cells += [thiscellx]
            yvertices_by_cells += [thiscelly]
            vertices_by_cells += [thiscell]

        """find list of Shapely Polygon objects for all cells"""
        for cell in vertices_by_cells:
            thispoly = shp.Polygon(cell)
            polys_ListbyShapely += [thispoly]
        """instantiate GeoDataFrame of Voronoi polygons"""
        self.gdf_vorPolys = gpd.GeoDataFrame(geometry=polys_ListbyShapely, crs=crs)
        self.gdf_vorPolys["cell"] = self.gdf_vorPolys.index.astype(str)
        """set x and y list attributes"""
        self.x_coords_by_node = xvertices_by_cells
        self.y_coords_by_node = yvertices_by_cells

        return self.gdf_vorPolys

    def get_gdf_topbtm(self, rasters: list, labels: list | tuple = ('bottom', 'top')):

        gdf_topbtm = self.get_centroid_elevations(
            elevations_files=rasters,
            labels=labels
        )
        top_botm_diff = gdf_topbtm['top'] - gdf_topbtm['bottom']
        for i in range(len(gdf_topbtm)):
            if top_botm_diff[i] < 0:
                print(f'{i} is bad. Setting bottom to 1 foot below top')
                gdf_topbtm.loc[i, 'bottom'] = gdf_topbtm.loc[i, 'top'] - 1
            gdf_topbtm.loc[i, 'diff'] = top_botm_diff[i]
        return gdf_topbtm

    def get_centroid_elevations(self, elevations_files: list, labels: list):
        n_rasters = len(elevations_files)

        # Open the elevations raster files
        srcs = []
        for i in range(n_rasters):
            srcs.append(rasterio.open(elevations_files[i]))
        # Get the elevations for each raster file as a numpy array
        elevations = []
        for i in range(n_rasters):
            elevations.append(srcs[i].read(1))

        # Create a function to get the elevation at a point for each raster file
        def get_elevations(x, y):
            elevs = []
            for i in range(n_rasters):
                row, col = srcs[i].index(x, y)
                elev = elevations[i][row-1, col-1]  # subtract 1 so rows and cols start at zero, else Python error
                elevs.append(elev)
            return tuple(elevs)

        # Get the centroids of the polygons as a geodataframe
        centroids_gdf = self.gdf_vorPolys.copy()
        centroids_gdf.geometry = centroids_gdf.centroid
        # Get the elevation at each centroid for each raster file
        for i in range(n_rasters):
            label = labels[i]
            centroids_gdf[label] = centroids_gdf.apply(
                lambda row: (get_elevations(
                    row.geometry.x, row.geometry.y)[i]
                ), axis=1
            )
        # Close the raster files
        for i in range(n_rasters):
            srcs[i].close()

        return centroids_gdf

    def get_cell_areas(self):
        print('getting cell areas')
        ### Find cell areas for DISU Package
        poly_area_list = []
        for poly in self.gdf_vorPolys['geometry']:
            poly_area_list += [poly.area]
        return poly_area_list

    def get_iac(self):
        ### Find Parameter iac for DISU Package
        print('getting iac for disu')
        num_adjacent_by_cell = [None for cell in self.adjacent_cells_idx]
        for i, cell in enumerate(self.adjacent_cells_idx):
            num_adjacent_by_cell[i] = len(cell) + 1
        return num_adjacent_by_cell

    def get_nja(self):
        print('getting nja for disu')
        ### Find Parameter nja for DISU Package ###
        sumiac = 0
        for i, num in enumerate(self.iac):
            sumiac += self.iac[i]
        nja_for_disu = sumiac
        return nja_for_disu

    def get_ja_cl12_hwva(self):
        print('getting ja, cl12, and hwva for disu')
        ### Find Parameter ja, cl12, and hwva for DISU Package ###
        ja_for_disu = []
        cl12_for_disu = []
        hwva_for_disu = []
        for i, cell in enumerate(self.adjacent_cells_idx):
            """get a sorted list of the indices of adjacent cells to this cell"""
            sortedcell = sorted(self.adjacent_cells_idx[i])
            """get list of indices starting with the reference cell and then 
            appending the indices of all cells adjacent to that reference cell"""
            thisja = [i] + sortedcell
            """get a list of distances between the centroid of the reference cell
            and the shared face of each adjacent cell"""
            thiscl12 = [0] + [self.calculate_distance(
                gdf=self.gdf_vorPolys,
                poly_idx1=i,
                poly_idx2=j,
            ) for j in sortedcell
            ]
            """get a list of the shared face lengths for each index in the 
            thisja list"""
            thishwva = [0] + [self.shared_face_length(
                poly1=self.gdf_vorPolys.loc[i].geometry,
                poly2=self.gdf_vorPolys.loc[j].geometry,
            ) for j in sortedcell
            ]
            ja_for_disu += thisja
            cl12_for_disu += thiscl12
            hwva_for_disu += thishwva
        return ja_for_disu, cl12_for_disu, hwva_for_disu

    def get_origin_xy(self):
        """Get the x,y coordinates for the origin of the 
        model grid.

        Returns:
            tuple: return a tuple of the form - (x,y) 
        """

        df_verts = pd.DataFrame(self.verts)
        xmin = df_verts.iloc[:, 0].min()
        ymin = df_verts[df_verts.iloc[:, 0] == xmin].iloc[:, 1].min()
        origin_xy = (xmin, ymin)
        self.origin_xy = origin_xy

        return origin_xy

    def get_latslons(self):
        """Generate a JSON representation of the Voronoi grid Polygons
        """
        gdf_latlon = self.gdf_vorPolys.to_crs(self.crs_latlon)
        latslons = json.loads(gdf_latlon["geometry"].to_json())

        self.gdf_latlon = gdf_latlon
        self.latslons = latslons

        return latslons

    def get_grid_centroid(self):
        """Gets the Shapley representation of the centroid of the defined voronoi grid

        Returns:
            shp: returns Shapely representation of the grid centroid
        """
        grid_centroid = shp.MultiPolygon(self.gdf_latlon['geometry'].to_list()).centroid
        return grid_centroid

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
            gdf_to_query: GeoDataFrame with geometry to query; this will take precedence over shp_to_query if it is passed to the function
            crs (string): coordinate reference system for the shapefile
            crs_latlon (string): crs for lat/lon, EPSG:4326. Shouldn't need to change this
            name_col (string): string corresponding to the field name in the shapefile to be used
            as for the keys in the return dict.
            predicate (string): method for query, defaults to intersects, but can use any predicate allowed
            by GeoPandas sindex.query
        """

        if shp_to_query:
            gdf_query = gpd.read_file(shp_to_query).to_crs(crs)
        if gdf_to_query:
            gdf_query = gdf_to_query

        """Create dictionary with Vornoi cell indices that contain/intersect each location"""
        VorIdx_dict = {}
        for idx in gdf_query.index:
            thiswell = (
                gpd.GeoSeries(gdf_query.to_crs(crs_latlon).iloc[idx]["geometry"])
                .sindex.query(
                    self.gdf_vorPolys["geometry"].to_crs(crs_latlon), predicate=predicate
                )[0]
                .tolist()
            )
            if name_col is not None:
                VorIdx_dict[gdf_query.iloc[idx][name_col]] = thiswell
            elif name_col is None:
                VorIdx_dict[idx] = thiswell

        return VorIdx_dict

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
        this_idx = (
            self.gdf_vorPolys["geometry"]).sindex.query(
            gdf_to_query['geometry'].iloc[idx], predicate=predicate, sort=True).tolist()
        return this_idx

    def set_k_vor(
            self,
            k_dict: dict = None,
            k_default=100
    ) -> list:
        """Method to set the hydraulic conductivity for all voronoi cells
        in the model grid. k_dict should be generated using the get_vor_idx_from_geometry
        method. The k_default is assigned to all cells prior to k_dict, in case some
        cells are not in k_dict."""

        """set initial default Kh"""
        self.gdf_vorPolys["Kh"] = k_default
        """set K values based on a dictionary of k values imported from shapefile"""
        for key in k_dict.keys():
            self.gdf_vorPolys.loc[self.gdf_vorPolys.index.isin(k_dict[key]), ["Kh"]] = key
        """make list of K values for each Voronoi cell"""
        k_vorcell_list = self.gdf_vorPolys["Kh"].to_list()

        return k_vorcell_list

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
        - hull_bounds: A tuple (min_x, min_y, max_x, max_y) representing the bounds of the area.
        - pixel_size: The size of each pixel in the same units as the hull_bounds.
        - output_filename: The filename for the output raster.
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
        elevation_data = np.zeros((height, width), dtype=rasterio.float32)

        # Calculate the elevation values
        print('getting elevations for raster')
        for row in range(height):
            for col in range(width):
                # Calculate the position of the current pixel
                x = min_x + col * pixel_size
                y = max_y - row * pixel_size

                # Project the point onto the normal vector
                point_vector = np.array([x - known_x, y - known_y, 0])
                distance_along_normal = np.dot(point_vector, normal) / np.linalg.norm(normal)
                elevation = known_elevation - distance_along_normal * np.cos(dip_rad)

                # Set the elevation value in the raster array
                elevation_data[row, col] = elevation
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
                crs='EPSG:2927',
                transform=transform,
        ) as dst:
            dst.write(elevation_data, 1)

        centroids_gdf = self.get_centroid_elevations([output_filename], ['elev'])

        return centroids_gdf

    @staticmethod
    def get_normal_from_strike_and_dip(strike: int, dip: int) -> np.array:
        """
        Returns the normal vector of a plane based on strike and dip
        :param strike: strike in degrees
        :param dip: dip in degrees
        :return: normal vector
        """
        strike_rad = np.radians(strike)
        dip_rad = np.radians(90 - dip)  # Convert dip into inclination from the vertical
        normal = np.array([
            np.sin(dip_rad) * np.sin(strike_rad),
            np.sin(dip_rad) * np.cos(strike_rad),
            np.cos(dip_rad)
        ])
        return normal

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

        # Create grid of points
        minx, miny, maxx, maxy = (center_point.x - size / 2 * spacing, center_point.y - size / 2 * spacing,
                                  center_point.x + size / 2 * spacing, center_point.y + size / 2 * spacing)

        x_coords = list(range(int(minx), int(maxx) + 1, spacing))
        y_coords = list(range(int(miny), int(maxy) + 1, spacing))

        points = [shp.Point(x, y) for x in x_coords for y in y_coords]

        # Convert to GeoSeries
        geoseries = gpd.GeoSeries(points).set_crs(crs)

        return geoseries

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

        # Create grid of points
        minx, miny, maxx, maxy = (center_point.x - size / 2 * spacing, center_point.y - size / 2 * spacing,
                                  center_point.x + size / 2 * spacing, center_point.y + size / 2 * spacing)

        x_coords = list(range(int(minx), int(maxx) + 1, spacing))
        y_coords = list(range(int(miny), int(maxy) + 1, spacing))

        polygons = []
        for x in x_coords[:-1]:  # We exclude the last coordinate because we're creating boxes "between" the points
            for y in y_coords[:-1]:
                # Create a box (polygon) for each pair of coordinates
                polygons.append(shp.box(x + gap / 2, y + gap / 2, x + spacing - gap / 2, y + spacing - gap / 2))

        # Convert to MultiPolygon
        multipolygon = shp.MultiPolygon(polygons)

        return multipolygon

    @staticmethod
    def voronoi_refine_by_point(point: shp.Point, spacing: int, tri: Triangle) -> gpd.GeoDataFrame:

        polypoints = [(point.x, point.y),
                      (point.x + spacing, point.y),
                      (point.x + spacing, point.y - spacing),
                      (point.x, point.y - spacing)]
        poly_main = shp.Polygon(polypoints)
        transform_dist = spacing * 2

        polyE = shp.transform(poly_main, lambda x: x + [transform_dist, 0])
        polyW = shp.transform(poly_main, lambda x: x - [transform_dist, 0])
        polyN = shp.transform(poly_main, lambda x: x + [0, transform_dist])
        polyS = shp.transform(poly_main, lambda x: x - [0, transform_dist])
        polyNE = shp.transform(poly_main, lambda x: x + [transform_dist, transform_dist])
        polyNW = shp.transform(poly_main, lambda x: x + [-transform_dist, transform_dist])
        polySE = shp.transform(poly_main, lambda x: x + [transform_dist, -transform_dist])
        polySW = shp.transform(poly_main, lambda x: x + [-transform_dist, -transform_dist])
        allPolys = [poly_main, polyE, polyW, polyN, polyS, polyNE, polyNW, polySE, polySW]

        gdf_allPolys = gpd.GeoDataFrame(geometry=allPolys)
        for poly in range(len(gdf_allPolys)):
            tri.add_polygon(gdf_allPolys.loc[poly, 'geometry'])

        return gdf_allPolys
