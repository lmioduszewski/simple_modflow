"""Streamflow-routing helpers for building MF6 SFR packages from geometries."""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

import flopy
import geopandas as gpd
from pathlib import Path
import numpy as np
import pandas as pd
import pickle
from simple_modflow.modflow.mf6.simulation.packages import _maybe_create_package_artifact
from simple_modflow.modflow.mf6.surface_water_validation import validate_sfr_configuration


class SFR:
    """Build SFR reach, connection, and period data from stream geometries."""

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            stream_paths: list[Path] = None,
            reverse_streams: list[bool] = None,
            inflows: dict | int | float | list = None,
            diversion_perioddata: dict = None,
            widths: list | int | float = 10,
            gradients: list | float = 0.001,
            mannings: list | float = 0.03,
            streambed_k: list | int | float = 1,
            streambed_thickness: list | int | float = 1,
            stream_end_conn: dict = None,
            div_prioritization='FRACTION',
            mover: bool = False,
            add_sfr=True,
            validate: bool = True,
            nper=None,
            register_regions: bool = False,
            region_name_prefix: str | None = None,
            combined_region_name: str | None = None,
            region_tags: list[str] | None = None,
            overwrite_regions: bool = False,
            artifact_id: str | None = None,
            artifact_catalog=None,
            artifact_description: str | None = None,
            artifact_tags: list[str] | None = None,
            artifact_metadata: dict | None = None,
            artifact_overwrite: bool = False,
    ):
        """Parameters
        ----------
        model
            Target model that will receive the SFR package.
        vor
            Grid helper used to map stream geometries to cells.
        stream_paths
            Path(s) to line geometries describing stream centerlines.
        reverse_streams
            Optional per-stream flags indicating whether stream ordering should be reversed.
        inflows
            Optional stress-period inflow data by reach.
        diversion_perioddata
            Optional explicit diversion period records.
        widths, gradients, mannings, streambed_k, streambed_thickness
            Per-reach or scalar hydraulic properties.
        stream_end_conn
            Optional rules describing how stream endpoints connect/divert.
        div_prioritization
            Diversion prioritization mode such as ``FRACTION``.
        mover
            Whether the resulting SFR package should be MVR-compatible.
        add_sfr
            Whether to build and attach the SFR package immediately.
        validate
            Whether to validate the derived SFR reach and connection inputs
            before attaching the package.
        nper
            Optional explicit number of stress periods.
        register_regions, region_name_prefix, combined_region_name, region_tags,
        overwrite_regions
            Optional region-registration behavior for stream footprints/reaches.
        artifact_id, artifact_catalog, artifact_description, artifact_tags,
        artifact_metadata, artifact_overwrite
            Optional reusable-artifact capture settings.
        """
        self.valid_packagedata_names = ['rlen', 'rwid', 'rgrd', 'rtp', 'rbth', 'rhk',
                                        'man', 'ncon', 'ustrf', 'ndv', 'aux', 'boundname']
        self.model = model
        self.vor = model.vor if vor is None else vor
        self.nper = nper if model is None else model.nper
        self.stream_paths = stream_paths
        self.reverse_streams = reverse_streams if reverse_streams else [False] * len(stream_paths)
        self._stream_cells = None
        self._stream_reaches = None
        self._stream_polys = None
        self._ncon = None
        self._connectiondata = None
        self._div_prioritization = None
        self._diversion_perioddata = None
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
        self._diversions = None
        self.register_regions = register_regions
        self.region_name_prefix = region_name_prefix
        self.combined_region_name = combined_region_name
        self.region_tags = [] if region_tags is None else list(region_tags)
        self.overwrite_regions = overwrite_regions
        self.artifact_id = artifact_id
        self.artifact_catalog = artifact_catalog
        self.artifact_description = artifact_description
        self.artifact_tags = artifact_tags
        self.artifact_metadata = artifact_metadata
        self.artifact_overwrite = artifact_overwrite
        self.package_artifact = None
        self.validation_report = None

        self.inflows = inflows
        self.widths = widths
        self.gradients = gradients
        self.mannings = mannings
        self.rhk = streambed_k
        self.rbth = streambed_thickness
        self.stream_endpoint_connections = stream_end_conn
        self.div_prioritization = div_prioritization
        self.diversion_perioddata = diversion_perioddata

        if add_sfr:
            self.add_sfr(validate=validate)
        elif self.register_regions and self.model is not None:
            self.register_model_regions()

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
    def stream_reaches(self):
        """return a nested list of stream reach indices, one list per stream"""
        if self._stream_reaches is None:
            stream_reaches = []
            for riv in self.stream_cells:
                cell_to_rno = pd.DataFrame.from_dict(self.cell_to_rno_dict, orient='index')
                riv_rno = cell_to_rno.loc[riv][0].to_list()
                riv_rno = [i + 1 for i in riv_rno]
                stream_reaches.append(riv_rno)
            self._stream_reaches = stream_reaches
        return self._stream_reaches

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
            for s in self.stream_cells:
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
            if self.stream_endpoint_connections is not None:
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
    def div_prioritization(self):
        return self._div_prioritization

    @div_prioritization.setter
    def div_prioritization(self, prior: str):
        allowed_proirs = ['FRACTION', 'EXCESS', 'THRESHOLD', 'UPTO']
        assert prior in allowed_proirs, f'Diversion prioritization must be one of {allowed_proirs}'
        self._div_prioritization = prior

    @property
    def diversions(self):
        """gets diversions data for direct input into the SFR flopy package"""
        if self._diversions is None:
            if self.mapped_diversions is None:
                return None
            diversions = []
            for ustr_cell, dstr_cell in self.mapped_diversions.items():
                ustr_reach = self.cell_to_rno_dict[ustr_cell]
                dstr_reach = self.cell_to_rno_dict[dstr_cell]
                idiv = 0  # assumes only one diversion per reach allowed
                # TODO update to allow multiple
                diversions.append([ustr_reach, idiv, dstr_reach, self.div_prioritization])
                # TODO allow more than one type of diversion in package
            self._diversions = diversions
        return self._diversions

    @property
    def connectiondata(self):
        if self._connectiondata is None:
            self.get_connection_data()
        return self._connectiondata

    @property
    def inflows(self):
        return self._inflows

    @inflows.setter
    def inflows(self, val):
        """setter that adds 'inflow' in the middle of each tuple in the provided val."""
        val = self.perioddata_validator('inflow', val)
        self._inflows = val

    @property
    def diversion_perioddata(self):
        return self._diversion_perioddata

    @diversion_perioddata.setter
    def diversion_perioddata(self, val):
        val = self.perioddata_validator('diversion', val)
        self._diversion_perioddata = val

    @property
    def perioddata(self):
        """builds perioddata for input to flopy sfr class"""
        perioddata = {}
        for per in range(self.nper):
            perioddata[per] = []
        if self.inflows is not None:
            for per in perioddata.keys():
                for setting in self.inflows[per]:
                    perioddata[per].append(setting)
        if self.diversion_perioddata is not None:
            for per in perioddata.keys():
                for setting in self.diversion_perioddata[per]:
                    perioddata[per].append(setting)
        self._sfr_period_data = perioddata
        return self._sfr_period_data

    def perioddata_validator(self, name: str, data: dict):

        allowed_settings = ['STATUS', 'MANNING', 'STAGE', 'INFLOW', 'RAINFALL', 'EVAPORATION',
                            'RUNOFF', 'DIVERSION', 'UPSTREAM_FRACTION', 'AUXILIARY']
        if data is None:
            return data
        assert isinstance(data, dict), f'data must be dict type'
        assert name.upper() in allowed_settings, f'Perioddata {name} must be one of {allowed_settings}'
        assert all(per in data.keys() for per in list(range(self.nper))), 'data keys must be valid stress periods'
        for per, settings in data.items():
            assert isinstance(settings, list), f'setting for {name} must be a list of tuples or lists'
            new_settings = []
            for setting in settings:
                assert setting[0] in list(range(self.total_nreaches)), 'first item in dict value must be a valid reach'
                new_data = (setting[0], name.upper(), *setting[1:])
                new_settings.append(new_data)
            data[per] = new_settings
        assert len(data) >= self.nper
        return data

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
        """Sets number of diversions for each reach. Sets to 0 unless the reach is upstream of a diversion,
                then ndv is 1"""
        if self.mapped_diversions is None:
            return None
        if self._ndv is None:
            ndv = []
            for s in self.stream_cells:
                s_ndv = []
                for cell in s:
                    if self.mapped_diversions is not None:
                        if cell in self.mapped_diversions.keys():
                            s_ndv.append(1)  # set to 1 if the cell is upstream of a diversion
                        else:
                            s_ndv.append(0)
                    else:
                        s_ndv.append(0)
                ndv.append(s_ndv)
            self._ndv = ndv
        return self._ndv

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
        """Sets upstream flow fraction for each reach. Sets to 1.0 unless the reach is downstream of a diversion,
        then ustrf is 0.0, and the upstream flow is set in the diversions in perioddata"""
        if self.mapped_diversions is None:
            return None
        if self._ustrf is None:
            ustrf = []
            for s in self.stream_cells:
                s_ustrf = []
                for cell in s:
                    if self.mapped_connections is not None:
                        if cell in self.mapped_diversions.values():
                            s_ustrf.append(0.0)  # set to 0.0 if the cell is downstream of a diversion
                        else:
                            s_ustrf.append(1.0)
                    else:
                        s_ustrf.append(1.0)
                ustrf.append(s_ustrf)
            self._ustrf = ustrf
        return self._ustrf

    @property
    def ncon(self):
        if self._ncon is None:
            self.get_connection_data()
        return self._ncon

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
                sfr_reach_data['rno'][rno] = rno
                sfr_reach_data['cellid'][rno] = (0, cell)
                sfr_reach_data['rlen'][rno] = self.reach_lens[stream_idx][cell_idx]
                sfr_reach_data['rwid'][rno] = self.widths[stream_idx][cell_idx]
                sfr_reach_data['rgrd'][rno] = self.gradients[stream_idx][cell_idx]
                sfr_reach_data['rtp'][rno] = self.get_smoothed_reach_elevs()[stream_idx][cell]
                sfr_reach_data['rbth'][rno] = self.rbth[stream_idx][cell_idx]
                sfr_reach_data['rhk'][rno] = self.rhk[stream_idx][cell_idx]
                sfr_reach_data['man'][rno] = self.mannings[stream_idx][cell_idx]
                sfr_reach_data['ncon'][rno] = self.ncon[stream_idx][cell_idx]
                sfr_reach_data['ustrf'][rno] = self.ustrf[stream_idx][cell_idx] if self.ustrf is not None else 1
                if self.ndv is not None:
                    sfr_reach_data['ndv'][rno] = self.ndv[stream_idx][cell_idx]

                rno += 1
        assert int(rno) == int(self.total_nreaches), \
            (f'something is wrong, total number of reaches in reach data {rno} '
             f'is incorrect, total_nreaches  is {self.total_nreaches}')

        return sfr_reach_data.tolist()

    @property
    def packagedata(self):
        return self.get_reach_data()

    def add_sfr(self, validate: bool = True):
        """Build and optionally validate the MF6 SFR package before attaching it."""

        self.validation_report = validate_sfr_configuration(self)
        if validate:
            self.validation_report.raise_for_errors("SFR validation failed.")
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
            perioddata=self.perioddata,
            maximum_picard_iterations=1,
            maximum_iterations=1000,
            maximum_depth_change=0.01,
            budget_filerecord='sfr_budget.sfr',
            stage_filerecord='sfr_stage.sfr',
            length_conversion=3.28081,  # since we are using feet instead of meters
            time_conversion=86_400,  # since we are using days instead of seconds
            mover=self.mover,
            diversions=self.diversions
        )
        self.model._sfr_input = self
        if self.register_regions and self.model is not None:
            self.register_model_regions()
        self.package_artifact = _maybe_create_package_artifact(
            self.model,
            "sfr",
            artifact_id=self.artifact_id,
            artifact_catalog=self.artifact_catalog,
            artifact_description=self.artifact_description,
            artifact_tags=self.artifact_tags,
            artifact_metadata=self.artifact_metadata,
            artifact_overwrite=self.artifact_overwrite,
        )
        return self.sfr

    def register_model_regions(self):
        if self.model is None:
            return {}

        registered = {}
        combined_cellids = []
        for stream_idx, stream_cells in enumerate(self.stream_cells):
            if not stream_cells:
                continue
            base_name = None
            if self.stream_paths is not None and stream_idx < len(self.stream_paths):
                base_name = Path(self.stream_paths[stream_idx]).stem
            if base_name is None:
                base_name = f"stream_{stream_idx}"
            region_name = (
                f"{self.region_name_prefix}_{base_name}"
                if self.region_name_prefix is not None
                else base_name
            )
            cellids = [(0, cell) for cell in stream_cells]
            geometry = self.stream_geoms[stream_idx]
            metadata = {"stream_index": stream_idx, "num_reaches": len(stream_cells)}
            registered[base_name] = self.model.add_region_from_cells(
                region_name,
                cellids=cellids,
                category="boundary",
                package="sfr",
                tags=self.region_tags or ["sfr"],
                geometry=geometry,
                metadata=metadata,
                overwrite=self.overwrite_regions,
            )
            combined_cellids.extend(cellids)

        if self.combined_region_name is not None and combined_cellids:
            registered["__combined__"] = self.model.add_region_from_cells(
                self.combined_region_name,
                cellids=combined_cellids,
                category="boundary",
                package="sfr",
                tags=self.region_tags or ["sfr"],
                metadata={"num_streams": self.num_streams},
                overwrite=self.overwrite_regions,
            )

        return registered

    def get_reach_lens(self):
        # vor_idxs = [sorted(cells) for cells in self.stream_cells]
        all_reach_lens = []
        for stream_num, geom in enumerate(self.stream_geoms):
            reach_lens = []
            reach_vor_polys = geom.intersection(self.vor.gdf_vorPolys.geometry)
            reach_vor_polys = reach_vor_polys.apply(lambda x: np.nan if x.length == 0 else x).dropna()
            for vor_idx in self.stream_cells[stream_num]:
                if vor_idx in reach_vor_polys.index:
                    reach_len = reach_vor_polys.loc[vor_idx].length
                else:
                    # Some branched or kinked stream lines can still map a valid
                    # stream cell whose first-pass aligned intersection collapses
                    # to zero length. Recompute directly against that cell so the
                    # reach-length builder remains stable for compact synthetic
                    # workflows and real diverted networks.
                    reach_len = self.vor.gdf_vorPolys.geometry.loc[vor_idx].intersection(geom).length
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
            if not intersecting_cells.empty:
                lengths = intersecting_cells.geometry.intersection(geom).length
                intersecting_cells = intersecting_cells.loc[lengths > 0].copy()
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
        ncon = []  # list to store number of connections per reach for packagedata
        rno = 0  # starting reach number
        if self.mapped_connections is not None:
            inverse_mapped_connections = {v: k for k, v in self.mapped_connections.items()}
        else:
            inverse_mapped_connections = {}
        for stream_idx, stream_cells in enumerate(self.stream_cells):
            stream_ncon = []
            for cell_idx, cell in enumerate(stream_cells):
                connections = [rno]  # Start with the reach number
                if cell_idx > 0:  # if not the first reach
                    connections.append(rno - 1)  # add upstream connection (positive index)
                if cell_idx < self.num_reach_cells_per_stream[stream_idx] - 1:  # if not the last reach in this stream
                    connections.append(-(rno + 1))  # add downstream connection (negative index)
                if self.mapped_connections is not None:
                    if cell in self.mapped_connections.keys():  # reach has mapped downstream connection to another stream
                        downstream_reach = self.cell_to_rno_dict[self.mapped_connections[cell]]
                        connections.append(-downstream_reach)  # add downstream connection (negative index)
                if cell in inverse_mapped_connections.keys():  # reach has mapped upstream connection to another stream
                    upstream_reach = self.cell_to_rno_dict[inverse_mapped_connections[cell]]
                    connections.append(upstream_reach)  # add upstream connection (positive index)

                stream_ncon.append(len(connections) - 1)
                connection_data.append(connections)
                rno += 1
            ncon.append(stream_ncon)
        self._ncon = ncon
        self._connectiondata = connection_data

        return connection_data

    def show_stream(self):
        all_cells = []
        for one_stream in self.stream_cells:
            all_cells = all_cells + one_stream
        self.vor.show_selected_cells(all_cells)


if __name__ == "__main__":
    vor_path = Path(r'C:\Users\lukem\Python\MODFLOW\LakePointe\new_vor_lakepointe.vor')
    with open(vor_path, 'rb') as file:
        vor: Vor = pickle.load(file)
    with open(Path(r"C:\Users\lukem\Python\MODFLOW\LakePointe\LakePointe.model"), 'rb') as file:
        model = pickle.load(file)
    input_path = Path(r'C:\Users\lukem\Python\MODFLOW\LakePointe\inputs')
    stream = input_path / r'shp\rivers and streams\jenkins.shp'
    sfr = SFR(vor=vor, stream_paths=stream, model=model)
