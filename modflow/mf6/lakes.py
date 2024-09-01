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
from simple_modflow import read_gpkg, read_shp_gpkg
import pandas as pd


class LakeAreaVolumeRelationship:
    """Defines the lake elevation, surface area, volume relationship for a given lake. The lake table
    will be generated and associated with the provided GWF model provided. The relationship
    calculations are based on a provided DEM, lake footprint shapefile, and a list of elevations. The lake
    number should be defined (default is zero), and needs to be the same lake number referenced in the other
    lake package classes during setup."""

    def __init__(
            self,
            dem_path: Path = None,
            lake_shapefile: Path = None,
            elevations: list = None,
            lake_num: int = 0,
            model: SimulationBase = None,
            shapefile_buffer=0,
            filename=None,
            area: int | float = None,
            top: int | float = None,
            btm: int | float = None,
            storage_coeff: int | float = 1,
            get_elevation_every: int | float = 1
    ):
        self.dem_path = dem_path
        self.lake_shapefile = lake_shapefile
        self.shapefile_buffer = shapefile_buffer
        if dem_path and lake_shapefile:
            self.load_dem()
            self.load_lake_footprint()
        self._lake_table = None
        self.elevations = elevations
        self.lake_num = lake_num
        self.model = model
        self.filename = f'lake_table.lak' if filename is None else filename
        self.table_block = [self.lake_num, self.filename]
        self.area = area
        self.top = top
        self.btm = btm
        self.storage_coeff = storage_coeff
        self.get_elevation_every = get_elevation_every

        if all([self.model, self.dem_path, self.lake_shapefile]):
            self.add_to_model()
            if any([self.top, self.btm, self.area]):
                print(f'Ignoring top, btm, and area arugments, since a dem path and shapefile were provided')
        elif all([self.top, self.btm, self.area]):
            self.rectangular_facility()
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

    def rectangular_facility(
            self,
            area: int | float = None,
            top: int | float = None,
            btm: int | float = None,
            storage_coeff: int | float = None,
            get_elevation_every: int = None,
            model: SimulationBase = None
    ):
        """gets a table for a simple rectangular-shaped buried infiltration facility"""
        model = model if model is not None else self.model
        area = area if area is not None else self.area
        top = top if top is not None else self.top
        btm = btm if btm is not None else self.btm
        storage_coeff = storage_coeff if storage_coeff is not None else self.storage_coeff
        get_elevation_every = get_elevation_every if get_elevation_every is not None else self.get_elevation_every

        tot_vol = area * (top - btm) * storage_coeff
        elev = btm + get_elevation_every
        table = []
        while elev <= top:
            this_vol = area * (elev - btm) * storage_coeff
            table.append([elev, this_vol, area])
            elev = elev + get_elevation_every
        self._lake_table = table
        return table

    def add_to_model(self):

        laktab = ModflowUtllaktab(
            model=self.model.gwf,
            nrow=len(self.lake_table),
            ncol=3,
            table=self.lake_table,
            filename=self.filename,
            pname=f'laktab_{self.lake_num}'
        )


class LakeConnectionData:

    def __init__(
            self,
            vor: Vor,
            paths: Path | list = None,
            bed_leakance: list | int = 1,
            only_layer=None,
            horizontal_connections: dict = None,
            use_reconciled_surfaces: bool = True,
            alt_surface_df=None,
            min_sep=0.1,
            only_vertical: bool = False,
            verbose: bool = False
    ):
        """
        Provide one or more Path objects that are shapefiles. Use property connection_data as an argument in the
        LakePackageData class.
        :param vor: VoronoiGridPlus voronoi grid object
        :param paths: one or more Path objects, a geopackage, a shapefile, or a list of shapefiles.
                      Not a list of geopackages.
        :param bed_leakance: an integer defining the lakebed leakance, or a list of integers if multiple lakes.
                             If only one value is provided for multiple lakes, it will be used for all lakes.
        :param only_layer: if you only want the lake to be connected to one particular model layer. This is a zero index.
        :param horizontal_connections: dict of lake numbers that have horizontal connections and top and bottom of those lakes.

                for example {0: [380, 370]} for a one lake model with a top of 380 and bottom of 370.
        """

        self.vor = vor
        self._lakes = None
        self._num_lakes = None
        self._bed_leakance = None
        self._connection_data = None
        self._lakes_vor_cells = None
        self._horizontal_connections = None
        self.use_reconciled_surfaces = use_reconciled_surfaces
        self.alt_surface_df = alt_surface_df
        self.min_sep = min_sep
        self.only_vertical = only_vertical

        if paths:
            self.lakes = paths
        if bed_leakance:
            self.bed_leakance = bed_leakance
        if horizontal_connections:
            self.horizontal_connections = horizontal_connections

        self.only_layer = only_layer

    @property
    def horizontal_connections(self) -> dict:
        return self._horizontal_connections

    @horizontal_connections.setter
    def horizontal_connections(self, val):
        assert isinstance(val, dict), 'horizontal_connections must be a dict'
        for key, value in val.items():
            assert key in list(range(self.num_lakes)), 'dict keys must be integers less than num lakes'
            assert isinstance(value, list), 'dict values must be lists'
            assert len(value) == 2, 'dict value lists must be length 2, top and bottom of lake horizontal connections'
            assert all(isinstance(i, int | float) for i in value), 'value list items must be integers or floats'
        self._horizontal_connections = val

    @property
    def num_lakes(self):
        self._num_lakes = len(self.lakes) if self.lakes is not None else 0
        return self._num_lakes

    @property
    def lakes(self) -> gpd.GeoDataFrame:
        """A GeoDataFrame containing the lake footprints."""
        return self._lakes

    @lakes.setter
    def lakes(self, val):

        if isinstance(val, Path):
            if val.suffix == 'shp':
                try:
                    self._lakes = gpd.read_file(val)
                except:
                    raise ValueError(f'cannot read shapefile {val}')
            elif val.suffix == 'gpkg':
                try:
                    self._lakes = read_gpkg(val)
                except:
                    raise ValueError(f'cannot read geopackage {val}')
            else:
                raise ValueError('path must be a .shp or .gpkg')

        elif isinstance(val, list):
            are_paths = all(isinstance(i, Path) for i in val)
            assert are_paths, 'all items in list must be Path'
            all_are_shp = all(path.suffix == '.shp' or path.suffix == '.gpkg' for path in val)
            assert all_are_shp, 'all lakes must have .shp or .gpkg extension'
            try:
                lakes = [read_shp_gpkg(lake_shp) for lake_shp in val]
                geoms = [lake.geometry for lake in lakes]
                self._lakes = gpd.GeoDataFrame.from_records(geoms).set_geometry(0)
            except:
                raise ValueError(f'cannot read lakes list {val}')

    @property
    def bed_leakance(self) -> list:
        """lakebed leakance, one value for each lake"""
        return self._bed_leakance

    @bed_leakance.setter
    def bed_leakance(self, val) -> list:
        if isinstance(val, int):
            if self.num_lakes > 1:
                print(f'applying bed leakance {val} to all {self.num_lakes} lakes')
            self._bed_leakance = [val] * self.num_lakes
        elif isinstance(val, list):
            assert len(val) == self.num_lakes, f'list must have {self.num_lakes} leakance values, one for each lake'
            assert all(isinstance(num, int | float) for num in val), 'all lakes must have integer leakance values'
            self._bed_leakance = val

    @property
    def lakes_vor_cells(self) -> dict:
        """Dict of lake cells for each lake. Keys are lake numbers and values are lake cells"""
        if self.num_lakes == 0:
            return None
        lakes_vor_cells = self.vor.get_vor_cells_as_series(self.lakes, return_dict=True)
        self._lakes_vor_cells = lakes_vor_cells
        return self._lakes_vor_cells

    @property
    def connection_data(self):
        if self._connection_data is None:
            self._connection_data = self.get_connection_data()
        return self._connection_data

    def custom_horizontal_connections(
            self,
            cell_id,
            lak_idx_conn,
            only_layer,
            bed_leakance,
            connection_data,
            elev_df,
            start_index,
            lake_num,
            adjacent_cells,
            top: int | float = None,
            botm: int | float = None

    ) -> tuple:
        """
        If you want to force certain horizontal connections of the 'lake'. Useful if using this class to define a
        rectangular lake to simulate an infiltration vault, trench, or other rectangular lake-like structure, such
        as buried facilities with storage rock. You can define a top and bottom elevation for all lake connections to
        create a box-like lake which may contain some sort of media.
        :param cell_id: cell id
        :param lak_idx_conn:
        :param only_layer:
        :param bed_leakance:
        :param connection_data:
        :param elev_df:
        :param start_index:
        :param lake_num:
        :param adjacent_cells:
        :param top:
        :param botm:
        :return:
        """

        vor = self.vor
        if not top or not botm:
            top, botm = self.horizontal_connections[lake_num][0], self.horizontal_connections[lake_num][1]

        for idx, ja_cell in enumerate(adjacent_cells):
            conn_len = vor.cl12[start_index + idx + 1]  # add 1 to skip the cell itself
            conn_width = vor.hwva[start_index + idx + 1]
            n_layers = len(vor.gdf_topbtm.drop('geometry', axis=1).columns) - 1

            for layer in range(n_layers):
                if only_layer is not None:
                    if layer != only_layer:
                        continue

                horiz_conn = [lake_num, lak_idx_conn, (layer, cell_id), 'HORIZONTAL',
                              bed_leakance[lake_num], botm, top, conn_len, conn_width]
                connection_data.append(horiz_conn)
                lak_idx_conn += 1

        return cell_id, lak_idx_conn, only_layer, bed_leakance, connection_data, elev_df, start_index, lake_num

    def get_connection_data(self):

        """The meat...gets the connection data for the MODFLOW 6 LAK package based on the provided voronoi grid
        and a shapefile of the lake footprint"""

        use_reconciled_surfaces = self.use_reconciled_surfaces
        alt_surface_df = self.alt_surface_df
        min_sep = self.min_sep
        only_vertical = self.only_vertical
        only_layer = self.only_layer
        vor = self.vor
        bed_leakance = self.bed_leakance

        connection_data = []

        if alt_surface_df is not None:
            elev_df = alt_surface_df
        else:
            elev_df = vor.reconcile_surfaces(self.min_sep) \
                if use_reconciled_surfaces else vor.gdf_topbtm

        for lake_num, lake_cells in self.lakes_vor_cells.items():
            print(f'Getting connection data for lake {lake_num}')
            if self.horizontal_connections:
                if lake_num in self.horizontal_connections.keys():
                    print(f'Using custom horizontal connections for lake {lake_num}')
            lak_idx_conn = 0  # starting index for numbering each lake's connections
            for cell_id in lake_cells:
                start_index = sum(vor.iac[:cell_id])  # index to start when looking up conn_len and conn_width
                #  set vertical connection for cell
                vert_conn = [lake_num, lak_idx_conn, (0, cell_id), 'VERTICAL', bed_leakance[lake_num], 0, 0, 0, 0]
                connection_data.append(vert_conn)
                lak_idx_conn += 1
                adjacent_cells = vor.find_adjacent_cells(cell_id)

                if only_vertical:
                    continue
                if self.horizontal_connections:
                    if lake_num in self.horizontal_connections.keys():
                        (cell_id, lak_idx_conn, only_layer, bed_leakance,
                         connection_data, elev_df, start_index, lake_num) = self.custom_horizontal_connections(
                            cell_id, lak_idx_conn, only_layer, bed_leakance,
                            connection_data, elev_df, start_index, lake_num, adjacent_cells)
                        continue

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
                        else:
                            # if this cell's top is not higher than adjacent cell top
                            # then no more connection data needed
                            break
                        horz_conn = [lake_num, lak_idx_conn, (layer, cell_id), 'HORIZONTAL',
                                     bed_leakance[lake_num], botm, top, conn_len, conn_width]
                        connection_data.append(horz_conn)
                        lak_idx_conn += 1

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

        # get the number of connections for each lake
        length_conns = pd.DataFrame(connectiondata[0]).loc[:, 0].value_counts()

        print(f'Getting package data for {nlakes} lakes')
        for lake in range(nlakes):
            lak_starting_stage = starting_stage[lake]
            lak_packagedata = [lake, lak_starting_stage, length_conns[lake]]
            packagedata.append(lak_packagedata)

        return packagedata


class LakePeriodData:

    def __init__(
            self,
            model: SimulationBase = None,
            lake_ids=[0],
            lake_stages=None,
            rainfall_rates=None,
            evaporation_rates=None,
            withdrawals=None,
            inflow=None,
            status=None,
            nper=1
    ):
        self.model = model
        self.lake_ids = lake_ids
        self.lake_stages = lake_stages
        self.rainfall_rates = rainfall_rates
        self.evaporation_rates = evaporation_rates
        self.withdrawals = withdrawals
        self.inflow = inflow
        self.status = status
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
            withdrawals=self.withdrawals,
            inflow=self.inflow,
            status=self.status
        )
        return self._period_data

    def get_lak_period_data(
            self,
            nper: int = None,
            lake_ids: list = [0],
            lake_stages: list = None,
            rainfall_rates: list = None,
            evaporation_rates: list = None,
            withdrawals: list = None,
            inflow: list = None,
            status: list = None,
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
        inflow (list of float, optional): List of inflow rates for each stress period.
        status (list of float, optional): List of lake status for each stress period.

        Returns:
        dict: Dictionary of period data for the LAK package. Each key is a stress period number, and the value is a list of lists.
        """
        period_data = {}
        nper = self.nper if nper is None else nper

        for period in range(nper):
            period_data[period] = []

            if len(lake_ids) == 1:

                for lake_id in lake_ids:
                    lake_number = lake_id  # lake_id is 0-based

                    laksetting = [lake_number]
                    if lake_stages is not None and period < len(lake_stages):
                        laksetting.extend(['stage', lake_stages[period]])
                    if rainfall_rates is not None:
                        laksetting.extend(['rainfall', rainfall_rates[period]])
                    if evaporation_rates is not None:
                        laksetting.extend(['evaporation', evaporation_rates[period]])
                    if withdrawals is not None:
                        laksetting.extend(['withdrawal', withdrawals[period]])
                    if inflow is not None:
                        laksetting.extend(['inflow', inflow[period]])
                    if status is not None:
                        laksetting.extend(['status', status[period]])

                    if len(laksetting) > 1:
                        period_data[period].append(laksetting)

            elif len(lake_ids) > 1:

                for lake_number in lake_ids:
                    if lake_stages is not None and lake_number not in lake_stages.keys():
                        continue
                    laksetting = [lake_number]
                    if lake_stages is not None and period < len(lake_stages[lake_number]):
                        laksetting.extend(['stage', lake_stages[lake_number][period]])
                    if rainfall_rates is not None:
                        laksetting.extend(['rainfall', rainfall_rates[lake_number][period]])
                    if evaporation_rates is not None:
                        laksetting.extend(['evaporation', evaporation_rates[lake_number][period]])
                    if withdrawals is not None:
                        laksetting.extend(['withdrawal', withdrawals[lake_number][period]])
                    if inflow is not None:
                        laksetting.extend(['inflow', inflow[lake_number][period]])
                    if status is not None:
                        laksetting.extend(['status', status[lake_number][period]])

                    if len(laksetting) > 1:
                        period_data[period].append(laksetting)

        return period_data
