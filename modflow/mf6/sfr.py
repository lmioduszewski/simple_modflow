import flopy
import pandas as pd
import geopandas as gpd
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
from pathlib import Path
import shapely as shp
import numpy as np

def flatten(l):
    return [item for sublist in l for item in sublist]

class SFR:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            stream_path: Path = None
    ):
        print('initing sfr')
        self.model = model
        self.vor = vor
        self.stream_path = stream_path
        self.stream_gdf = gpd.read_file(stream_path)
        self.stream_geom = self.stream_gdf.geometry
        self.stream_points = {}
        self.vor_intersect = {}
        self.vor_intersect_lists = []
        self._reach_lens = None
        self.sfr = None

        for i, geom in enumerate(self.stream_geom):
            print(f'\rGetting geometry and intersect lists {i}', flush=True)
            x, y = geom.coords.xy
            self.stream_points[i] = [shp.Point(point) for point in list(zip(x,y))]
            self.vor_intersect[i] = vor.get_vor_cells_as_series(geom).to_list()

        print(self.vor_intersect[0])
        for i in range(len(self.vor_intersect)):
            print(f'\radding lists together {i}', flush=True)
            self.vor_intersect_lists.extend(self.vor_intersect[i])
        print('sorting')
        self.vor_intersect_flat_list = sorted(self.vor_intersect_lists)
        self.nreaches = len(self.vor_intersect_flat_list)  # Number of stream reaches
        self.stream_cells = None
        self.sfr_reach_data = None
        self._sfr_connection_data = None
        self.sfr_period_data = None

        print('getting reach data')
        self.get_reach_data()
        # self.add_sfr()


    @property
    def reach_lens(self):
        if self.reach_lens is None:
            reach_lens = self.get_reach_lens()
            self._reach_lens = reach_lens
        return self._reach_lens


    def sfr_connection_data(self, idx: int = 0, reverse: bool = False):
        if self._sfr_connection_data is None:
            stream_cells, _ = self.get_sorted_cells_along_stream(idx, reverse)
            self._sfr_connection_data = self.get_connection_data(stream_cells)
        return self._sfr_connection_data

    def get_gradient(self):
        pass

    def get_reach_data(self):

        # Define the stream network data
        nreaches = self.nreaches
        sfr_cells = [(0, i) for i in self.vor_intersect_flat_list]  # cells where the stream reaches are located in lyr 1

        # Define SFR package data
        sfr_reach_data = np.zeros(nreaches, dtype=[
            ('rno', int),  # Reach number
            ('cellid', object),  # Cell ID
            ('rlen', float),  # Reach length in feet
            ('rwid', float),  # Reach width in feet
            ('rgrd', float),  # Reach gradient (dimensionless)
            ('rtp', float),  # Reach top elevation in feet
            ('rbth', float),  # Reach bottom thickness in feet
            ('rhk', float),  # Reach hydraulic conductivity in feet/day
            ('man', float),  # Manning's roughness coefficient
            ('ncon', int),  # Number of connections
            ('ustrf', float),  # Upstream fraction
            ('ndv', int)  # Number of downstream diversions
        ])

        # Populate the reach data
        for i, cell in enumerate(sfr_cells):
            sfr_reach_data['rno'][i] = i + 1
            sfr_reach_data['cellid'][i] = cell
            sfr_reach_data['rlen'][i] = self.reach_lens[i]  # Length of each reach in feet
            sfr_reach_data['rwid'][i] = 10  # Width of each reach in feet
            sfr_reach_data['rgrd'][i] = 0.0015  # Gradient of each reach (dimensionless)
            sfr_reach_data['rtp'][i] = self.vor.gdf_topbtm[0][self.vor_intersect_flat_list[idx]]  # Top elevation of each reach in feet
            sfr_reach_data['rbth'][i] = 1  # Thickness of the streambed in feet
            sfr_reach_data['rhk'][i] = 5  # Hydraulic conductivity of the streambed in feet/day
            sfr_reach_data['man'][i] = 0.03  # Manning's roughness coefficient (example value)
            sfr_reach_data['ncon'][i] = 1  # Number of connections (example value)
            sfr_reach_data['ustrf'][i] = 1.0  # Upstream fraction (example value)
            sfr_reach_data['ndv'][i] = 0  # Number of downstream diversions (example value)
            print(f'\r{i}, {cell}', flush=True)

        self.sfr_period_data = {0: [(0, 'inflow', 5.0)]}  # Inflow of 5 cubic feet per day at the first reach
        self.sfr_reach_data = sfr_reach_data

    def add_sfr(self):

        # Create the SFR package
        self.sfr = flopy.mf6.ModflowGwfsfr(
            self.model.gwf,
            save_flows=True,
            print_input=True,
            print_flows=True,
            pname='sfr',
            nreaches=self.nreaches,
            packagedata=self.sfr_reach_data,
            connectiondata=self._sfr_connection_data,
            perioddata=self.sfr_period_data
        )
        return self.sfr

    def get_reach_lens(self):
        vor_idxs = self.vor_intersect_flat_list
        reach_lens = []
        for idx in vor_idxs:
            reach_len = self.stream_geom.unary_union.intersection(vor.gdf_vorPolys.loc[idx]).geometry.length
            reach_lens.append(reach_len)
        return reach_lens

    def get_sorted_cells_along_stream(self, idx=0, reverse=False):
        """
        method to get a list of voronoi model cells along the length of a strewam in order from
        the beginning to end of the stream linestring.
        :param idx: index of the stream, if there was more than one provided in the shapefile
        :param reverse: set to True if the list needs to be reversed because the stream line
        start and end need to be reversed
        :return: list of stream voronoi cells, dataframe of intersecting cells
        """
        stream = self.stream_geom[idx]
        intersecting_cells = self.vor.gdf_vorPolys[self.vor.gdf_vorPolys.intersects(stream)].copy()
        # Calculate the centroid of each intersecting cell
        intersecting_cells['centroid'] = intersecting_cells.centroid

        # Calculate the distance along the stream line for each cell centroid
        intersecting_cells['distance_along_stream'] = intersecting_cells.centroid.apply(
            lambda point: stream.project(point)
        )
        # Sort the cells based on the distance along the stream line
        sorted_cells = intersecting_cells.sort_values('distance_along_stream')
        sorted_cell_ids = sorted_cells.index.tolist()
        if reverse:
            sorted_cell_ids = sorted_cell_ids[::-1]

        return sorted_cell_ids, sorted_cells

    def get_smoothed_reach_elevs(self, idx=0, reverse=False):

        # get cells along the stream
        cells = self.get_sorted_cells_along_stream(idx=idx, reverse=reverse)[0]
        #  get a copy of a dataframe with the top of model elevations for each of the stream cells
        stream_cells = self.vor.gdf_topbtm.loc[cells, :][0].copy()
        stream_cells = stream_cells.to_dict()
        stream_cell_keys = list(stream_cells.keys())
        for idx, (cell, elev) in enumerate(stream_cells.items()):
            next_cell_idx = idx + 1
            try:
                next_cell = stream_cell_keys[next_cell_idx]
            except:
                continue
            next_cell_elev = stream_cells[next_cell]
            if next_cell_elev > elev:
                stream_cells[next_cell] = elev
        return stream_cells

    def get_connection_data(self, cell_ids):
        """
        Computes the SFR connection data for the MODFLOW 6 SFR package.

        Parameters:
        - cell_ids: List of cell IDs ordered from the beginning to the end of the stream.

        Returns:
        - connection_data: List of connection data in the format required by the SFR package.
        """
        nreaches = len(cell_ids)
        connection_data = []

        for i in range(nreaches):
            if i == 0:
                # First reach has no upstream connections
                connection_data.append([i + 1, -(i + 2)])  # Connection to the downstream reach
            elif i == nreaches - 1:
                # Last reach has no downstream connections
                connection_data.append([i + 1, 0])  # No further connections
            else:
                # Middle reaches have both upstream and downstream connections
                connection_data.append([i + 1, -(i + 2)])

        return connection_data