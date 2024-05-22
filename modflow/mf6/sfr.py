import flopy
import pandas as pd
import geopandas as gpd
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
from pathlib import Path
import shapely as shp

def flatten(l):
    return [item for sublist in l for item in sublist]

class SFR:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            stream_path: Path = None
    ):
        self.model = model
        self.vor = vor
        self.stream_path = stream_path
        self.stream_gdf = gpd.read_file(stream_path)
        self.stream_geom = self.stream_gdf.geometry
        self.stream_points = {}
        self.vor_intersect = {}
        self.vor_intersect_lists = []
        self._reach_lens = None

        for i, geom in enumerate(self.stream_geom):
            x, y = geom.coords.xy
            self.stream_points[i] = [shp.Point(point) for point in list(zip(x,y))]
            self.vor_intersect[i] = vor.get_vor_cells_as_series(geom)

        self.vor_intersect_lists = [
            self.vor_intersect_lists + self.vor_intersect[i].to_list()
            for i in range(len(self.vor_intersect))
        ]
        self.vor_intersect_flat_list = sorted(flatten(self.vor_intersect_lists))
        self.nreaches = len(self.vor_intersect_flat_list)  # Number of stream reaches
        self.sfr_reach_data = None
        self.sfr_connection_data = None
        self.sfr_period_data = None

        self.get_reach_data()
        self.add_sfr()


    @property
    def reach_lens(self):
        if self.reach_lens is None:
            reach_lens = self.get_reach_lens()
            self._reach_lens = reach_lens
        return self._reach_lens
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

        # Populate the reach data (replace with your actual data)
        for i, cell in enumerate(sfr_cells):
            ja_cells = self.vor.ja
            sfr_reach_data['rno'][i] = i + 1
            sfr_reach_data['cellid'][i] = cell
            sfr_reach_data['rlen'][i] = self.reach_lens[i]  # Length of each reach in feet
            sfr_reach_data['rwid'][i] = 10  # Width of each reach in feet
            sfr_reach_data['rgrd'][i] = 0.001  # Gradient of each reach (dimensionless)
            sfr_reach_data['rtp'][i] = self.vor.gdf_topbtm[0][self.vor_intersect_flat_list[idx]]  # Top elevation of each reach in feet
            sfr_reach_data['rbth'][i] = 1  # Thickness of the streambed in feet
            sfr_reach_data['rhk'][i] = 5  # Hydraulic conductivity of the streambed in feet/day
            sfr_reach_data['man'][i] = 0.03  # Manning's roughness coefficient (example value)
            sfr_reach_data['ncon'][i] = 1  # Number of connections (example value)
            sfr_reach_data['ustrf'][i] = 1.0  # Upstream fraction (example value)
            sfr_reach_data['ndv'][i] = 0  # Number of downstream diversions (example value)

        # Define connection data for each reach (example connections)
        sfr_connection_data = [
            [1, -2],  # Reach 1 connects downstream to reach 2
            [2, -3],  # Reach 2 connects downstream to reach 3
            [3, -4],  # Reach 3 connects downstream to reach 4
            [4, -5],  # Reach 4 connects downstream to reach 5
            [5, 0]  # Reach 5 has no further downstream connections (outflow)
        ]
        self.sfr_period_data = {0: [(0, 'inflow', 100.0)]}  # Inflow of 100 cubic feet per day at the first reach
        self.sfr_reach_data = sfr_reach_data
        self.sfr_connection_data = sfr_connection_data

    def add_sfr(self):

        # Create the SFR package
        sfr = flopy.mf6.ModflowGwfsfr(
            self.model.gwf,
            save_flows=True,
            print_input=True,
            print_flows=True,
            pname='sfr',
            nreaches=self.nreaches,
            packagedata=self.sfr_reach_data,
            connectiondata=self.sfr_connection_data,
            perioddata=self.sfr_period_data
        )

    def get_reach_lens(self):
        vor_idxs = self.vor_intersect_flat_list
        reach_lens = []
        for idx in vor_idxs:
            reach_len = self.stream_geom.unary_union.intersection(vor.gdf_vorPolys.loc[idx]).geometry.length
            reach_lens.append(reach_len)
        return reach_lens

    def get_sorted_cells_along_stream(self, idx):

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

        return sorted_cell_ids