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
            stream_paths: list[Path] = None,
            reverse_streams: list[bool] = None,
            inflows: dict | int | float | list = None,
            widths: list | int | float = None,
            gradients: list | float = None,
            mannings: list | float = None,
            streambed_k: list | int | float = None,
            streambed_thickness: list | int | float = None,
            upstream_fraction: list | int | float = None,
            stream_end_conn: dict = None,
            diversions: dict = None,
            mover: bool = False,
            add_sfr=True
    ):
        """

        :param model: modflow model file, should be class SimulationBase
        :param vor: voronoi grid file, should be class Vor
        :param stream_paths:
        :param reverse_streams:
        :param stream_idx: arbitary index for the stream
        :param inflows: dict where keys are all stress periods and each value is a list of tuples with len 2. Each tuple = (reach id, inflow)
        :param widths: list of reach widths for the stream. List length must be equal to the number of reaches. Or may provide single value (int or float) for all reaches
        :param gradients: list of reach gradients. List length must be equal to the number of reaches or may provide single float to apply to all reaches
        :param mannings: list of reach manning's coefficients, or a float to apply to all reaches
        :param streambed_k:
        :param streambed_thickness:
        :param upstream_fraction:
        :param stream_end_conn: dict of sfr connections. Keys are stream indexes with an end point connection, values
        are tuples of length three. First tuple values are the stream indexes that connect to that end point. Second
        tuple values indicate which end of the key stream connects to the value stream, +1 for the start of the stream,
        -1 for the end of the stream. Thirds tuple values indicates whether the stream end upstream or downstream of
        the connecting stream, +1 for downstream and -1 for upstream, e.g. {0: (1, -1, 1)} means the 'end' end of
        stream 0 connects to and is downstream of stream 1.
        :param diversions: dict of sfr diversions. Keys are start reaches, values are target reaches, e.g., {3: 10} diverts from reach 3 to reach 10)
        :param mover: boolean value to indicate that this SFR package can be used with the water mover (MVR) package
        :param add_sfr:
        """
        print('initing sfr')

        self.valid_packagedata_names = ['rlen', 'rwid', 'rgrd', 'rtp', 'rbth', 'rhk',
                                        'man', 'ncon', 'ustrf', 'ndv', 'aux', 'boundname']
        self.model = model
        self.vor = model.vor if vor is None else vor
        self.stream_paths = stream_paths
        self.reverse_streams = reverse_streams if reverse_streams else [False] * len(stream_paths)
        self._stream_cells = None
        self._stream_polys = None
        self._sfr_connection_data = None
        self._sfr_period_data = None
        self._reach_lens = None
        self.sfr = None
        self._stream_geoms = None
        self._ndv = None
        self._rhk = None
        self._rbth = None
        self._ustrf = None
        self._stream_conn = None
        self._mapped_connections = None
        self._mapped_diversions = None
        self.mover = mover
        self._rno_to_cell_dict = None

        self.inflows = inflows
        self.widths = widths
        self.gradients = gradients
        self.mannings = mannings
        self.rhk = streambed_k
        self.rbth = streambed_thickness
        self.ustrf = upstream_fraction
        self.stream_endpoint_connections = stream_end_conn
        self.diversions = diversions

        if add_sfr:
            print('Processing streams')
            self.process_streams()
            print('Adding SFR package')
            self.add_sfr()

    @property
    def stream_geoms(self):
        if self._stream_geoms is None:
            geoms = [gpd.read_file(path).geometry.union_all() for path in self.stream_paths]
            self._stream_geoms = geoms
        return self._stream_geoms

    @property
    def reach_lens(self):
        if self._reach_lens is None:
            reach_lens = self.get_reach_lens()
            self._reach_lens = reach_lens
        return self._reach_lens

    @property
    def stream_cells(self):
        """list of stream cell indices for each stream"""
        if self._stream_cells is None:
            stream_cells = self.get_sorted_cells_along_stream(self.reverse_streams)
            stream_cells = [cells.index.to_list() for cells in stream_cells]
            self._stream_cells = stream_cells
        return self._stream_cells

    @property
    def stream_polys(self):
        """list of stream voronoi polygons for each stream"""
        if self._stream_polys is None:
            stream_polys = self.get_sorted_cells_along_stream(self.reverse_streams)
            stream_polys = [gdf.geometry for gdf in stream_polys]
            self._stream_polys = stream_polys
        return self._stream_polys

    @property
    def rno_to_cell_dict(self):
        """Returns a dict of reach numbers (keys) and their corresponding stream cell indices (values)"""
        if self._rno_to_cell_dict is None:
            rno_dict = {}
            rno = 0
            for s in sfr.stream_cells:
                for cell in s:
                    rno_dict[rno] = cell
                    rno += 1
            self._rno_to_cell_dict = rno_dict
        return self._rno_to_cell_dict

    @property
    def cell_to_rno_dict(self):
        """Gets dict of cell indices (keys) and corresponding reach numbers (values)"""
        return {v: k for k, v in self.rno_to_cell_dict.items()}

    @property
    def num_streams(self):
        return len(self.stream_cells)

    @property
    def stream_indices(self):
        return list(range(self.num_streams))

    @property
    def num_reach_cells_per_stream(self):
        num_cells_per_stream = []
        for one_stream in self.stream_cells:
            num_cells_per_stream.append(len(one_stream))
        return num_cells_per_stream

    @property
    def total_nreaches(self):
        return sum(self.num_reach_cells_per_stream)

    @property
    def stream_endpoint_connections(self):
        return self._stream_conn

    @stream_endpoint_connections.setter
    def stream_endpoint_connections(self, d: dict):
        if d is not None:
            assert isinstance(d, dict), 'stream_endpoint_connections must be a dict'
            assert all(key in self.stream_indices for key in list(d.keys())), f'dict keys must be stream indices'
            assert all(isinstance(value, tuple) for value in list(d.values())), 'dict values must be tuples'
            assert all(len(tup) == 3 for tup in list(d.values())), 'tuple values must be length 3'
            assert all(tup[0] in self.stream_indices for tup in d.values()), f'first tuple value must be a stream index'
            assert all(tup[1] in [1, -1] for tup in d.values()), f'second tuple value must be 1 or -1'
            assert all(tup[2] in [1, -1] for tup in d.values()), f'third tuple value must be 1 or -1'

        self._stream_conn = d

    @property
    def mapped_connections(self):
        """return a dict of mapped stream connections based on the stream geometries given
        and the stream_endpoint_connections"""
        if self._mapped_connections is None:
            mapped_connections = {}
            mapped_diversions = {}
            for stream_idx, stream_conn in self.stream_endpoint_connections.items():
                which_end = 0 if stream_conn[1] == 1 else -1
                end_poly = self.stream_polys[stream_idx].iloc[which_end]
                end_cell = self.stream_cells[stream_idx][which_end]
                # find distances from end poly to all connecting stream polys, and then find minimum distance
                distances = self.stream_polys[stream_conn[0]].apply(lambda x: x.distance(end_poly))
                closest_cell = distances.sort_values().index[0]
                upstream_reach = end_cell if stream_conn[2] == -1 else closest_cell
                downstream_reach = end_cell if stream_conn[2] == 1 else closest_cell
                mapped_connections[upstream_reach] = downstream_reach
                if downstream_reach == end_cell:
                    mapped_diversions[upstream_reach] = downstream_reach
            self._mapped_connections = mapped_connections
            self._mapped_diversions = mapped_diversions
        return self._mapped_connections

    @property
    def mapped_diversions(self):
        if self._mapped_diversions is None:
            _ = self.mapped_connections  # get mapped connections property, which also sets mapped diversions
        return self._mapped_diversions

    @property
    def connectiondata(self):
        if self._sfr_connection_data is None:
            self._sfr_connection_data = self.get_connection_data()
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
            assert len(val) >= self.model.nper
        elif isinstance(val, int | float):
            print(f'applying {val} as starting inflow to the first stream (index of 0) in all stress periods')
            val = {per: [(0, 'inflow', val)] for per in range(self.model.nper)}
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

    def package_data_validator(
            self,
            name: str = None,
            data: list[list[int | float] | int | float] = None
    ):
        """checks given data for the sfr packagedata. Makes data a list of length num_streams, with
        nested lists of length stream_cells"""

        name = self.packagedata_name_validator(name)
        if isinstance(data, int | float):
            data = [data for _ in range(self.num_streams)]
        if isinstance(data, list):
            assert len(data) == self.num_streams, f'length of {name} list must be equal to number of streams'
            for stream_idx in range(self.num_streams):
                if isinstance(data[stream_idx], int | float):
                    print(f'Stream {stream_idx}: applying {data[stream_idx]} for {name} to all reaches')
                    data[stream_idx] = [data[stream_idx] for _ in range(self.num_reach_cells_per_stream[stream_idx])]
                assert self.num_reach_cells_per_stream[stream_idx] == len(data[stream_idx]), \
                    (f'length of {name} data for stream {stream_idx} must be equal to number of stream cells:'
                     f' {self.num_reach_cells_per_stream[stream_idx]}')
            return data
        else:
            raise ValueError(f'{name} data provided invalid. Must be list, not type: {type(data)}')

    def packagedata_name_validator(self, name: str = None):
        assert name in self.valid_packagedata_names, \
            f'{name} not a valid name, can be one of {self.valid_packagedata_names}'
        return name

    @property
    def widths(self):
        return self._widths

    @widths.setter
    def widths(self, val):
        widths = self.package_data_validator(name='rwid', data=val)
        self._widths = widths

    @property
    def gradients(self):
        return self._gradients

    @gradients.setter
    def gradients(self, val):
        grad = self.package_data_validator(name='rgrd', data=val)
        self._gradients = grad

    @property
    def mannings(self):
        return self._mannings

    @mannings.setter
    def mannings(self, val):
        man = self.package_data_validator(name='man', data=val)
        self._mannings = man

    @property
    def ndv(self):
        return self._ndv

    @ndv.setter
    def ndv(self, val):
        ndv = self.package_data_validator(name='ndv', data=val)
        self._ndv = ndv

    @property
    def rhk(self):
        return self._rhk

    @rhk.setter
    def rhk(self, val):
        rhk = self.package_data_validator(name='rhk', data=val)
        self._rhk = rhk

    @property
    def rbth(self):
        return self._rbth

    @rbth.setter
    def rbth(self, val):
        rbth = self.package_data_validator(name='rbth', data=val)
        self._rbth = rbth

    @property
    def ustrf(self):
        return self._ustrf

    @ustrf.setter
    def ustrf(self, val):
        usrtf = self.package_data_validator(name='ustrf', data=val)
        self._ustrf = usrtf

    def get_gradient(self):
        """get reach gradients from the reach elevations and distances"""
        raise NotImplementedError

    def get_reach_data(self):

        # Define SFR package data
        sfr_reach_data = np.zeros(self.total_nreaches, dtype=[
            ('rno', int),  # Reach number
            ('cellid', object),  # Cell ID
            ('rlen', float),  # Reach length
            ('rwid', float),  # Reach width
            ('rgrd', float),  # Reach gradient (dimensionless)
            ('rtp', float),  # Reach top elevation
            ('rbth', float),  # Reach bottom thickness
            ('rhk', float),  # Reach hydraulic conductivity
            ('man', float),  # Manning's roughness coefficient
            ('ncon', int),  # Number of connections
            ('ustrf', float),  # Upstream fraction
            ('ndv', int)  # Number of downstream diversions
        ])

        rno = 0  # the starting reach number
        # Populate the reach data
        for stream_idx, stream_cells in enumerate(self.stream_cells):
            for cell_idx, cell in enumerate(stream_cells):
                i = 1
                if i == 0:
                    nconn = 1
                elif i == self.num_reach_cells_per_stream[stream_idx] - 1:
                    nconn = 1
                else:
                    nconn = 2
                sfr_reach_data['rno'][rno] = rno
                sfr_reach_data['cellid'][rno] = (0, cell)
                sfr_reach_data['rlen'][rno] = self.reach_lens[stream_idx][cell_idx]
                sfr_reach_data['rwid'][rno] = self.widths[stream_idx][cell_idx]
                sfr_reach_data['rgrd'][rno] = self.gradients[stream_idx][cell_idx]
                sfr_reach_data['rtp'][rno] = self.get_smoothed_reach_elevs()[stream_idx][cell]
                sfr_reach_data['rbth'][rno] = self.rbth[stream_idx][cell_idx]
                sfr_reach_data['rhk'][rno] = self.rhk[stream_idx][cell_idx]
                sfr_reach_data['man'][rno] = self.mannings[stream_idx][cell_idx]
                sfr_reach_data['ncon'][rno] = nconn
                sfr_reach_data['ustrf'][rno] = self.ustrf[stream_idx][cell_idx]
                sfr_reach_data['ndv'][rno] = 0

                rno += 1
        assert int(rno) == int(
            self.total_nreaches), f'something is wrong, total number of reaches in reach data {rno} is incorrect, total_nreaches  is {self.total_nreaches}'

        return sfr_reach_data.tolist()

    @property
    def packagedata(self):
        return self.get_reach_data()

    def add_sfr(self):
        # Create the SFR package
        self.sfr = flopy.mf6.ModflowGwfsfr(
            self.model.gwf,
            save_flows=True,
            print_input=True,
            print_flows=True,
            pname='sfr',
            nreaches=self.total_nreaches,
            packagedata=self.packagedata,
            connectiondata=self.connectiondata,
            perioddata=self.period_data,
            maximum_picard_iterations=1,
            maximum_iterations=1000,
            maximum_depth_change=0.01,
            budget_filerecord='sfr_budget.sfr',
            stage_filerecord='sfr_stage.sfr',
            length_conversion=3.28081,  # since we are using feet instead of meters
            time_conversion=86_400,  # since we are using days instead of seconds
            mover=self.mover
        )
        return self.sfr

    def get_reach_lens(self):
        # vor_idxs = [sorted(cells) for cells in self.stream_cells]
        all_reach_lens = []
        for stream_num, geom in enumerate(self.stream_geoms):
            reach_lens = []
            reach_vor_polys = geom.intersection(self.vor.gdf_vorPolys.geometry)
            reach_vor_polys = reach_vor_polys.apply(lambda x: np.nan if x.length == 0 else x).dropna()
            for vor_idx in self.stream_cells[stream_num]:
                reach_len = reach_vor_polys[vor_idx].length
                reach_lens.append(reach_len)
            all_reach_lens.append(reach_lens)
        return all_reach_lens

    def get_sorted_cells_along_stream(self, reverse=None):
        """
        method to get a list of voronoi model cells along the length of a stream in order from
        the beginning to end of the stream linestring.
        :param idx: index of the stream, if there was more than one provided in the shapefile
        :param reverse: set to True if the list needs to be reversed because the stream line
        start and end need to be reversed
        :return: list of stream voronoi cells, dataframe of intersecting cells
        """

        reverse = self.reverse_streams if reverse is None else reverse
        sorted_cells = []
        for idx, geom in enumerate(self.stream_geoms):
            intersecting_cells = self.vor.gdf_vorPolys[self.vor.gdf_vorPolys.intersects(geom)].copy()
            # Calculate the centroid of each intersecting cell
            intersecting_cells['centroid'] = intersecting_cells.centroid

            # Calculate the distance along the stream line for each cell centroid
            intersecting_cells['distance_along_stream'] = intersecting_cells.centroid.apply(
                lambda point: geom.project(point)
            )
            # Sort the cells based on the distance along the stream line
            sort = intersecting_cells.sort_values('distance_along_stream')
            if reverse[idx] is True:
                sort = sort.iloc[::-1]
            sorted_cells.append(sort)

        return sorted_cells

    def get_smoothed_reach_elevs(self):
        """adjust so that downstream reaches are the same elevation or lower elevation than upstream reaches"""
        all_streams_smoothed = []
        for stream_cell_list in self.stream_cells:
            #  get a copy of a dataframe with the top of model elevations for each of the stream cells
            stream_cells = self.vor.gdf_topbtm.loc[stream_cell_list, :][0].copy()
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
            all_streams_smoothed.append(stream_cells)
        return all_streams_smoothed

    def get_connection_data(self):
        """
        Computes the SFR connection data for the MODFLOW 6 SFR package.

        Returns:
        - connection_data: List of connection data in the format required by the SFR package.
        """

        connection_data = []
        rno = 0  # starting reach number
        for stream_idx, stream_cells in enumerate(self.stream_cells):
            for cell_idx, cell in enumerate(stream_cells):
                connections = [rno]  # Start with the reach number
                if cell_idx > 0:  # if not the first reach
                    connections.append(rno - 1)  # add upstream connection (positive index)
                if cell_idx < self.num_reach_cells_per_stream[stream_idx] - 1:  # if not the last reach in this stream
                    connections.append(-(rno + 1))  # add downstream connection (negative index)
                if cell in self.mapped_connections.keys():  # reach has mapped downstream connection to another stream
                    connections.append(-self.cell_to_rno_dict[cell])  # add downstream connection (negative index)
                if cell in self.mapped_connections.values():  # reach has mapped upstream connection to another stream
                    connections.append(self.cell_to_rno_dict[cell])  # add upstream connection (positive index)

                connection_data.append(connections)
                rno += 1

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
    sfr = SFR(vor=vor, stream_paths=stream, model=model)
