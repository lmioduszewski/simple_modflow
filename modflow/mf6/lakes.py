import geopandas as gpd
from shapely.geometry import Polygon
import rasterio
import numpy as np
from rasterio.mask import mask
from skimage.measure import find_contours
from rasterio.features import geometry_mask
from pathlib import Path
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
from flopy.mf6.modflow.mfutllaktab import ModflowUtllaktab

class LakeAreaVolumeRelationship:
    """Defines the lake elevation, surface area, volume relationship for a given lake. The lake table
    will be generated and associated with the provided GWF model provided. The relationship
    calculations are based on a provided DEM, lake footprint shapefile, and a list of elevations. The lake
    number should be defined (default is zero), and needs to be the same lake number referenced in the other
    lake package classes during setup."""
    def __init__(
            self,
            dem_path: Path,
            lake_shapefile: Path,
            elevations: list = None,
            nlake: int = 0,
            model: SimulationBase = None
    ):
        self.dem_path = dem_path
        self.lake_shapefile = lake_shapefile
        self.load_dem()
        self.load_lake_footprint()
        self._lake_table = None
        self.elevations = elevations
        self.nlake = nlake,
        self.model = model

        if self.model:
            self.add_to_model()

    @property
    def lake_table(self):
        if self._lake_table is None:
            if self.elevations:
                self._lake_table = self.create_table()
            else:
                raise ValueError('No elevations were provided')
        return self._lake_table

    def load_dem(self):
        with rasterio.open(self.dem_path) as src:
            self.dem = src.read(1)
            self.transform = src.transform
            self.bounds = src.bounds

    def load_lake_footprint(self):
        lake_gdf = gpd.read_file(self.lake_shapefile)
        self.lake_polygon = lake_gdf.geometry.iloc[0]
        self.clip_dem_to_lake()

    def clip_dem_to_lake(self):
        with rasterio.open(self.dem_path) as src:
            out_image, out_transform = mask(src, [self.lake_polygon], crop=True)
            self.dem = out_image[0]
            self.transform = out_transform
        self.lake_area = self.lake_polygon.area

    def get_polygon_for_elevation(self, elevation):
        contours = find_contours(self.dem, level=elevation)
        polygons = []
        for contour in contours:
            coords = [self.transform * (c[1], c[0]) for c in contour]
            polygon = Polygon(coords)
            if polygon.is_valid and polygon.area < self.lake_area * 0.9:  # Exclude entire footprint
                polygons.append(polygon)
        return gpd.GeoSeries(polygons).unary_union

    def analyze_lake(self):
        elevations = self.elevations
        results = []
        for elevation in elevations:
            polygon = self.get_polygon_for_elevation(elevation)
            area = polygon.area
            volume = self.calculate_volume(polygon, elevation)
            results.append({
                'elevation': elevation,
                'area': area,
                'volume': volume,
                'geometry': polygon
            })
        return results

    def create_table(self):
        results = self.analyze_lake()
        lake_table = []
        for result in results:
            row = [result['elevation'], result['volume'], result['area']]
            lake_table.append(row)
        return lake_table

    def calculate_volume(self, polygon, elevation):
        # Create a mask for the entire polygon area
        mask = geometry_mask([polygon], transform=self.transform, invert=True, out_shape=self.dem.shape)
        # Apply the mask to the DEM to isolate the lake area
        lake_dem = np.where(mask, self.dem, np.nan)
        # Calculate the depth of water above the DEM
        depth = elevation - lake_dem
        # Calculate the cell size (area) from the affine transformation
        cell_area = self.transform.a * -self.transform.e  # a is the pixel width, e is the pixel height
        # Calculate the volume of water above the DEM
        volume = np.nansum(depth[depth > 0]) * cell_area
        return volume

    def add_to_model(self):

        laktab = ModflowUtllaktab(
            model=self.model.gwf,
            nrow=len(self.lake_table),
            ncol=3,
            table=self.lake_table,
            filename='lake_table.lak',
            pname='lak'
        )


class LakeConnectionData:

    def __init__(self, vor: Vor, lake_shp: Path, bed_leakance = 1, nlake: int = 0):

        self.vor = vor
        self.lake_shp = lake_shp
        self.lake_vor_cells = vor.get_vor_cells_as_series(self.lake_shp)
        self._connection_data = None
        self.lake_no = nlake
        self.bed_leakance = bed_leakance

    @property
    def connection_data(self):
        if self._connection_data is None:
            self._connection_data = self.get_connection_data()
        return self._connection_data

    def get_connection_data(self, use_reconciled_surfaces=True):
        """gets the connection data for the MODFLOW 6 LAK package based on the provided vornoi grid
        and a shapefile of the lake footprint"""

        idx_conn = 0
        vor = self.vor
        lake_no = self.lake_no
        bed_leakance = self.bed_leakance
        connection_data = []
        elev_df = vor.reconcile_surfaces() if use_reconciled_surfaces else vor.gdf_topbtm

        for cell_id in self.lake_vor_cells:
            start_index = sum(vor.iac[:cell_id])
            #  set vertical connection for cell
            vert_conn = [lake_no, idx_conn, (0, cell_id), 'VERTICAL', bed_leakance, 0, 0, 0, 0]
            connection_data.append(vert_conn)
            idx_conn += 1

            adjacent_cells = vor.find_adjacent_cells(cell_id)
            ja_cell_top_elevs = elev_df[0].loc[adjacent_cells]
            cell_top_elev = elev_df[0].loc[cell_id]
            for idx, ja_cell in enumerate(adjacent_cells):
                """iterate through all adjacent cells and if the adjacent cell is lower in elevation
                then get parameters for a horizontal lake connection"""
                if ja_cell_top_elevs[ja_cell] < cell_top_elev:
                    top = cell_top_elev
                    botm = ja_cell_top_elevs[ja_cell]
                    conn_len = vor.cl12[start_index + idx + 1]  # add 1 to skip the cell itself
                    conn_width = vor.hwva[start_index + idx + 1]
                    #  set horizontal connection for cell
                    horz_conn = [lake_no, idx_conn, (0, cell_id), 'HORIZONTAL', bed_leakance, botm, top, conn_len,
                                 conn_width]
                    connection_data.append(horz_conn)
                    idx_conn += 1
        return connection_data

class LakePackageData:

    def __init__(
            self,
            nlakes: int,
            starting_stage: list,
            connectiondata: list,
            boundnames: list

    ):
        self._packagedata = None
        self.nlakes = nlakes
        self.starting_stage = starting_stage
        self.connectiondata = connectiondata
        self.boundnames = boundnames

    @property
    def packagedata(self):
        if self._packagedata is None:
            self._packagedata = self.get_packagedata()
        return self._packagedata

    def get_packagedata(self):

        nlakes = self.nlakes
        starting_stage = self.starting_stage
        connectiondata = self.connectiondata
        boundnames = self.boundnames
        packagedata = []

        for lake in range(nlakes):
            lak_starting_stage = starting_stage[lake]
            lakeconn = connectiondata[lake]
            lak_packagedata = [lake, lak_starting_stage, len(lakeconn)]
            packagedata.append(lak_packagedata)
        return packagedata