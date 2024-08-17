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
            lake_num: int = 0,
            model: SimulationBase = None,
            shapefile_buffer=0
    ):
        self.dem_path = dem_path
        self.lake_shapefile = lake_shapefile
        self.shapefile_buffer = shapefile_buffer
        self.load_dem()
        self.load_lake_footprint()
        self._lake_table = None
        self.elevations = elevations
        self.lake_num = lake_num
        self.model = model
        self.filename = f'lake_table.lak'
        self.table_block = [self.lake_num, self.filename]

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
        self.lake_polygon = lake_gdf.geometry.iloc[0].buffer(self.shapefile_buffer)
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
        return gpd.GeoSeries(polygons).union_all()

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
            filename=self.filename,
            pname='lak'
        )


class LakeConnectionData:

    def __init__(
            self, vor: Vor,
            lake_shp: Path,
            bed_leakance = 1,
            lake_num: int = 0,
            only_layer = None,
    ):

        self.vor = vor
        self.lake_shp = lake_shp
        self.lake_vor_cells = vor.get_vor_cells_as_series(self.lake_shp)
        self._connection_data = None
        self.lake_num = lake_num
        self.bed_leakance = bed_leakance
        self.only_layer = only_layer

    @property
    def connection_data(self):
        if self._connection_data is None:
            self._connection_data = self.get_connection_data()
        return self._connection_data

    def get_connection_data(self, use_reconciled_surfaces=True, alt_surface_df=None, min_sep=None, only_vertical=False):
        """gets the connection data for the MODFLOW 6 LAK package based on the provided vornoi grid
        and a shapefile of the lake footprint"""

        only_layer = self.only_layer
        idx_conn = 0  # starting index for numbering connections
        vor = self.vor
        lake_num = self.lake_num
        bed_leakance = self.bed_leakance
        connection_data = []
        min_sep = 0.1 if min_sep is None else min_sep
        if alt_surface_df is not None:
            elev_df = alt_surface_df
        else:
            elev_df = vor.reconcile_surfaces(min_sep=min_sep) if use_reconciled_surfaces else vor.gdf_topbtm
        adjusted_cells = {}

        for cell_id in self.lake_vor_cells:
            start_index = sum(vor.iac[:cell_id])
            #  set vertical connection for cell
            vert_conn = [lake_num, idx_conn, (0, cell_id), 'VERTICAL', bed_leakance, 0, 0, 0, 0]
            connection_data.append(vert_conn)
            idx_conn += 1

            if only_vertical:
                continue

            adjacent_cells = vor.find_adjacent_cells(cell_id)
            ja_cell_tops = elev_df[0].loc[adjacent_cells]
            cell_top = elev_df[0].loc[cell_id]

            for idx, ja_cell in enumerate(adjacent_cells):
                ja_cell_top = elev_df.loc[ja_cell, 0]
                conn_len = vor.cl12[start_index + idx + 1]  # add 1 to skip the cell itself
                conn_width = vor.hwva[start_index + idx + 1]
                #  get list of booleans that indicate what layers for this cell are higher than adjacent cell top
                cell_comp = (elev_df.loc[cell_id] > ja_cell_top)
                for layer, boolean in enumerate(cell_comp):
                    if only_layer is not None:
                        if layer != only_layer:
                            continue
                    if boolean:  # if cell's top of this layer higher than adjacent cell top
                        top = elev_df.loc[cell_id, layer]
                        # if cell's bottom of this layer is not higher than adjacent cell top
                        if not cell_comp[layer + 1]:
                            botm = ja_cell_top
                        else:  # but if this cell's bottom is higher than the adjacent cell top
                            botm = elev_df.loc[cell_id, layer + 1]
                    else:  # if this cell's top is not higher than adjacent cell top then no more connection data needed
                        break
                    horz_conn = [lake_num, idx_conn, (layer, cell_id), 'HORIZONTAL', bed_leakance, botm, top, conn_len,
                                 conn_width]
                    connection_data.append(horz_conn)
                    idx_conn += 1
        return connection_data

class LakePackageData:

    def __init__(
            self,
            nlakes: int,
            starting_stage: list,
            connectiondata: list = None,
            boundnames: list = None

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


class LakePeriodData:

    def __init__(
            self,
            model: SimulationBase = None,
            lake_ids = [0],
            lake_stages = None,
            rainfall_rates = None,
            evaporation_rates = None,
            withdrawals = None,
            nper = 1
    ):
        self.model = model
        self.lake_ids = lake_ids
        self.lake_stages = lake_stages
        self.rainfall_rates = rainfall_rates
        self.evaporation_rates = evaporation_rates
        self.withdrawals = withdrawals
        self.nper = nper if self.model is None else self.model.nper
        self._period_data = None

    @property
    def perioddata(self):
        self._period_data = self.get_lak_period_data(
            nper=self.nper,
            lake_ids=self.lake_ids,
            lake_stages=self.lake_stages,
            rainfall_rates=self.rainfall_rates,
            evaporation_rates=self.evaporation_rates,
            withdrawals=self.withdrawals
        )
        return self._period_data


    def get_lak_period_data(
            self,
            nper: int = None,
            lake_ids: list = [0],
            lake_stages: list = None,
            rainfall_rates: list = None,
            evaporation_rates: list = None,
            withdrawals: list = None
    ) -> list:
        """
        Generate period data for the MODFLOW 6 LAK package.

        Parameters:
        num_stress_periods (int): Number of stress periods.
        lake_ids (list of int): List of lake IDs.
        lake_stages (list of float, optional): List of lake stages for each stress period.
        rainfall_rates (list of float, optional): List of rainfall rates for each stress period.
        evaporation_rates (list of float, optional): List of evaporation rates for each stress period.
        withdrawals (list of float, optional): List of withdrawal rates for each stress period.

        Returns:
        dict: Dictionary of period data for the LAK package. Each key is a stress period number, and the value is a list of lists.
        """
        period_data = {}
        nper = self.nper if nper is None else nper

        for period in range(nper):
            period_data[period] = []

            for lake_id in lake_ids:
                lake_number = lake_id  # lake_id is 0-based

                laksetting = [lake_number]
                if lake_stages is not None and period < len(lake_stages):
                    print(lake_stages[period])
                    print(len(lake_stages))
                    laksetting.extend(['stage', lake_stages[period]])
                if rainfall_rates is not None:
                    laksetting.extend(['rainfall', rainfall_rates[period]])
                if evaporation_rates is not None:
                    laksetting.extend(['evaporation', evaporation_rates[period]])
                if withdrawals is not None:
                    laksetting.extend(['withdrawal', withdrawals[period]])

                period_data[period].append(laksetting)

        return period_data


