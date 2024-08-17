import flopy
import pandas as pd
import geopandas as gpd
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
from pathlib import Path
import shapely as shp
import numpy as np
import pickle

"""def flatten(l):
    return [item for sublist in l for item in sublist]"""

class SFR:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            stream_path: Path = None,
            reverse_stream: bool = False,
            stream_idx = 0,
            inflows: dict | int | float = None,
            widths: list | int | float = None,
            gradients: list | float = None,
            mannings: list | float = None,
            streambed_k: list | int | float = None,
            streambed_thickness: list | int | float = None,
            add_sfr = True
    ):
        """

        :param model: modflow model file, should be class SimulationBase
        :param vor: voronoi grid file, should be class Vor
        :param stream_path:
        :param reverse_stream:
        :param stream_idx: arbitary index for the stream
        :param inflows: dict where keys are all stress periods and each value is a list of tuples with len 2. Each tuple = (reach id, inflow)
        :param widths: list of reach widths for the stream. List length must be equal to the number of reaches. Or may provide single value (int or float) for all reaches
        :param gradients: list of reach gradients. List length must be equal to the number of reaches or may provide single float to apply to all reaches
        :param mannings: list of reach manning's coefficients, or a float to apply to all reaches
        :param streambed_k:
        :param streambed_thickness:
        :param add_sfr:
        """
        print('initing sfr')
        self.model = model
        self.vor = vor
        self.stream_idx = stream_idx
        self.stream_path = stream_path
        self.reverse = reverse_stream
        self._stream_cells = None
        self.sfr_reach_data = None
        self._sfr_connection_data = None
        self._sfr_period_data = None
        self.stream_gdf = gpd.read_file(stream_path)
        self.stream_geom = self.stream_gdf.geometry.union_all()
        self.stream_points = {}
        self._reach_lens = None
        self.sfr = None

        self.inflows = inflows
        self.widths = widths
        self.gradients = gradients
        self.mannings = mannings
        self.streambed_k = streambed_k
        self.streambed_thickness = streambed_thickness

        x, y = self.stream_geom.coords.xy
        self.stream_points = [shp.Point(point) for point in list(zip(x,y))]
        self.nrch = len(self.stream_cells)  # Number of stream reaches

        if add_sfr:
            print('getting reach data')
            self.get_reach_data()
            self.add_sfr()


    @property
    def reach_lens(self):
        if self._reach_lens is None:
            reach_lens = self.get_reach_lens()
            self._reach_lens = reach_lens
        return self._reach_lens

    @property
    def stream_cells(self):
        if self._stream_cells is None:
            stream_cells = self.get_sorted_cells_along_stream(self.stream_idx, self.reverse)[0]
            self._stream_cells = stream_cells
        return self._stream_cells

    @property
    def num_sfr_cells(self):
        return len(self.stream_cells)

    @property
    def sfr_connection_data(self):
        if self._sfr_connection_data is None:
            self._sfr_connection_data = self.get_connection_data(self.stream_cells)
        return self._sfr_connection_data

    @property
    def inflows(self):
        return self._inflows

    @inflows.setter
    def inflows(self, val):
        """setter that adds 'inflow' in the middle of each tuple in the provided val."""
        if isinstance(val, dict):
            for per in val.keys():
                tupls = []
                for tup in val[per]:
                    tupl = (tup[0], 'inflow', tup[1])
                    tupls.append(tupl)
                val[per] = tupls
        elif isinstance(val, int | float):
            print(f'applying {val} as starting inflow to the stream in all stress periods')
            starting_sfr_cell = self.stream_cells[0]
            val = {per: [(starting_sfr_cell, 'inflow', val)] for per in range(self.model.nper)}
        else:
            raise ValueError('must provide at least one inflow to the stream. The starting inflow')
        self._inflows = val

    @property
    def period_data(self):
        periodd = {}
        for per in range(self.model.nper):
            periodd[per] = []
        if self.inflows is not None:
            for per in periodd.keys():
                for rch_inflow in self.inflows[per]:
                    periodd[per].append(rch_inflow)
        self._sfr_period_data = periodd
        return self._sfr_period_data

    @property
    def widths(self):
        return self._widths

    @widths.setter
    def widths(self, val):
        if isinstance(val, list):
            assert len(val) == self.num_sfr_cells, 'length of width list must be equal to number of reaches'
            self._widths = val
        elif isinstance(val, int | float):
            print(f'applying {val} as width to all reaches')
            wid = [val for reach in range(self.num_sfr_cells)]
            self._widths = wid
        else:
            raise ValueError(f'width provided invalid. Must be list, int, or float, not type: {type(val)}')

    @property
    def gradients(self):
        return self._gradients

    @gradients.setter
    def gradients(self, val):
        if isinstance(val, list):
            assert len(val) == self.num_sfr_cells, 'length of gradients list must be equal to number of reaches'
            self._gradients = val
        elif isinstance(val, float):
            print(f'applying {val} as gradient to all reaches')
            grad = [val for reach in range(self.num_sfr_cells)]
            self._gradients = grad
        else:
            raise ValueError(f'gradients provided invalid. Must be list of floats or float, not type: {type(val)}')

    @property
    def mannings(self):
        return self._mannings

    @mannings.setter
    def mannings(self, val):
        if isinstance(val, list):
            assert len(val) == self.num_sfr_cells, 'length of mannings coefficients list must be equal to number of reaches'
            self._mannings = val
        elif isinstance(val, float):
            print(f'applying {val} as mannings coefficient to all reaches')
            man = [val for reach in range(self.num_sfr_cells)]
            self._mannings = man
        else:
            raise ValueError(f'mannings coeff provided invalid. Must be list of floats or float, not type: {type(val)}')

    def get_gradient(self):
        raise NotImplementedError

    def get_reach_data(self, elev_add=0.1, streambed_k = 1):

        # Define the stream network data
        nreaches = self.nrch
        sfr_cells = [(0, i) for i in self.stream_cells]  # cells where the stream reaches are located in lyr 1
        top_elevs = self.get_smoothed_reach_elevs()

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
            if i == 0:
                nconn = 1
            elif i == len(sfr_cells) - 1:
                nconn = 1
            else:
                nconn = 2
            sfr_reach_data['rno'][i] = i
            sfr_reach_data['cellid'][i] = cell
            sfr_reach_data['rlen'][i] = self.reach_lens[i]  # Length of each reach in feet
            sfr_reach_data['rwid'][i] = self.widths[i]  # Width of each reach in feet
            sfr_reach_data['rgrd'][i] = self.gradients[i]  # Gradient of each reach (dimensionless)
            sfr_reach_data['rtp'][i] = top_elevs[self.stream_cells[i]] + elev_add  # Top elevation of each reach in feet
            sfr_reach_data['rbth'][i] = 1  # Thickness of the streambed in feet
            sfr_reach_data['rhk'][i] = streambed_k  # Hydraulic conductivity of the streambed in feet/day
            sfr_reach_data['man'][i] = self.mannings[i]  # Manning's roughness coefficient (example value)
            sfr_reach_data['ncon'][i] = nconn  # Number of connections (example value)
            sfr_reach_data['ustrf'][i] = 1.0  # Upstream fraction (example value)
            sfr_reach_data['ndv'][i] = 0  # Number of downstream diversions (example value)
        self.sfr_reach_data = sfr_reach_data.tolist()

    def add_sfr(self):

        # Create the SFR package
        self.sfr = flopy.mf6.ModflowGwfsfr(
            self.model.gwf,
            save_flows=True,
            print_input=True,
            print_flows=True,
            pname='sfr',
            nreaches=self.nrch,
            packagedata=self.sfr_reach_data,
            connectiondata=self.sfr_connection_data,
            perioddata=self.period_data,
            maximum_picard_iterations=1,
            maximum_iterations=1000,
            maximum_depth_change=0.01,
            budget_filerecord='sfr_budget.sfr',
            stage_filerecord='sfr_stage.sfr',
            length_conversion=3.28081,  # since we are using feet instead of meters
            time_conversion=86_400  # since we are using days instead of seconds

        )
        return self.sfr

    def get_reach_lens(self):
        vor_idxs = sorted(self.stream_cells)
        reach_lens = []
        for idx in vor_idxs:
            reach_len = self.stream_geom.intersection(self.vor.gdf_vorPolys.loc[idx]).geometry.length
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
        stream = self.stream_geom
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

    def get_smoothed_reach_elevs(self, idx=0):

        # get cells along the stream
        cells = self.get_sorted_cells_along_stream(idx=idx, reverse=self.reverse)[0]
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
            connections = [i]  # Start with the reach number (0-based index)

            if i > 0:
                # Add upstream connection (positive index)
                connections.append(i - 1)
            if i < nreaches - 1:
                # Add downstream connection (negative index)
                connections.append(-(i + 1))

            connection_data.append(connections)

        return connection_data

    def show_stream(self):
        cells = self.stream_cells
        self.vor.show_selected_cells(cells)



if __name__ == "__main__":
    vor_path = Path(r'C:\Users\lukem\Python\MODFLOW\LakePointe\new_vor_lakepointe.vor')
    with open(vor_path, 'rb') as file:
        vor: Vor = pickle.load(file)
    with open(Path(r"C:\Users\lukem\Python\MODFLOW\LakePointe\LakePointe.model"), 'rb') as file:
        model = pickle.load(file)
    input_path = Path(r'C:\Users\lukem\Python\MODFLOW\LakePointe\inputs')
    stream = input_path / r'shp\rivers and streams\jenkins.shp'
    sfr = SFR(vor=vor, stream_path=stream, model=model)