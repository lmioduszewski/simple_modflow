from __future__ import annotations
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from .voronoiplus import VoronoiGridPlus as Vor
    from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
    from simple_modflow.modflow.mf6.sfr import SFR

import flopy
from pathlib import Path
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as Hp
from simple_modflow.modflow.utils.datatypes.surface_data import ModelSurface
from simple_modflow.modflow.utils.datatypes.choros import Choro
from simple_modflow.modflow.utils.datatypes.xsections import XSection
from shapely import LineString
import pickle
from simple_modflow.modflow.mf6.sfr import SFR
from simple_modflow.modflow.utils.inputs import Inputs
from simple_modflow.modflow.utils.outputs import LakOutputData, SFROutputData
from simple_modflow.modflow.mf6.budget import Budget
from simple_modflow.modflow.utils.datatypes.modelgrid import create_custom_modelgrid
import numpy as np


class SimulationBase:

    def __init__(
            self,
            name: str = 'mf6_model',
            mf_folder_path: Path = Path().home().joinpath('mf6'),
            nper: int = 1,
            vor: Vor = None,
            per_dates: list | pd.DatetimeIndex = None
    ):
        self.name = name
        self.vor = vor
        self.nlay = None
        self.nper = nper
        self.num_steps = None
        self.per_len = None
        self.model_output_folder_path = mf_folder_path.joinpath(f'{name}')
        self._master_celld = {}
        self._hds = None
        self._bud = None
        self._sfr_input = None
        self._per_dates = None
        self.per_dates = per_dates
        self._node_to_lni = None
        self._lni_to_node = None
        self._kstpkper = None

        self.sim = flopy.mf6.MFSimulation(
            sim_name=self.name,
            exe_name="mf6",
            version="mf6",
            sim_ws=self.model_output_folder_path,
        )
        self.gwf = flopy.mf6.ModflowGwf(
            self.sim,
            modelname=self.name,
            model_nam_file=f"{self.name}.nam",
            print_flows=False,
            newtonoptions='under_relaxation',
            save_flows=True
        )
        self.ims = flopy.mf6.modflow.mfims.ModflowIms(
            self.sim,
            print_option='ALL',
            pname="ims",
            complexity="COMPLEX",
            under_relaxation="DBD",
            under_relaxation_theta=0.72,
            under_relaxation_kappa=0.1,
            # under_relaxation_gamma=0.2,
            under_relaxation_momentum=0.001,
            backtracking_number=20,
            backtracking_tolerance=1.1,
            backtracking_reduction_factor=0.3,
            backtracking_residual_limit=75,
            outer_maximum=200,
            inner_maximum=200,
            outer_dvclose=0.3,
            inner_dvclose=0.2,
            # rcloserecord=[300_000, 'strict'],
            relaxation_factor=0.97,
            linear_acceleration='BICGSTAB',
        )

    @property
    def per_dates(self):
        return self._per_dates

    @per_dates.setter
    def per_dates(self, per_dates):
        if per_dates is not None:
            if isinstance(per_dates, list):
                try:
                    per_dates = pd.to_datetime(per_dates)
                except ValueError:
                    print("Can't convert provided period dates to Pandas DateTime")
            else:
                assert isinstance(per_dates, pd.DatetimeIndex), 'Cannot recognize valid dates in provide period dates'
        self._per_dates = per_dates

    @property
    def modelgrid(self) -> flopy.discretization.vertexgrid.VertexGrid:
        modelgrid: flopy.discretization.vertexgrid.VertexGrid = self.gwf.modelgrid
        return modelgrid

    @property
    def node_to_lni(self) -> dict:
        """
        create a dict where keys are model nodes and values are corresponding layer node indices,
        i.e. layer specific index of node
        :return: dict
        """
        if self._node_to_lni is None:

            node_to_lni = {}
            for node in range(self.modelgrid.nnodes):
                node_to_lni[node] = self.modelgrid.get_lni([node])[0]
            self._node_to_lni = node_to_lni

        return self._node_to_lni

    @property
    def lni_to_node(self) -> dict:
        """
        create a dict where keys are layer node indices and values are corresponding model nodes
        :return: dict
        """
        if self._lni_to_node is None:

            lni_to_node = {}
            for node in range(self.modelgrid.nnodes):
                lni_to_node[self.modelgrid.get_lni([node])[0]] = node
            self._lni_to_node = lni_to_node

        return self._lni_to_node

    @property
    def cellids(self):
        return list(self.lni_to_node.keys())

    @property
    def sfr_input(self) -> SFR:
        if self._sfr_input is not None:
            assert isinstance(self._sfr_input, SFR), 'sfr input not SFR class'
            sfr_input: SFR = self._sfr_input
            return sfr_input
        else:
            return None

    @property
    def hds(self):
        self._hds = Hp(model=self, vor=self.vor)
        return self._hds

    @property
    def all_heads(self):
        return self.hds.all_heads

    @property
    def surf(self):
        return ModelSurface(model=self)

    def choro(
            self,
            kstpkper: tuple = None,
            per: int = None,
            layer: int = 0,
            choro_type: str = 'hds',
            custom_hover: dict = None,
            custom_zs: list = None,
            zmin: float | int = None,
            zmax: float | int = None,
            zoom: int = 13,
            show_layer_elevs: bool = True,
            show_mounding: bool = False,
            hover_heads: bool = True,
            hover_ks: bool = False,
            locs=None,
            rch_scale=None,
            bgs=False,
            **kwargs
    ) -> Choro:
        """
        kwargs can be any allowable keyword arguments from the Choro class
        Class defining the basic choropleth plots generated from a modflow model.

        :param rch_scale:
        :param kstpkper:
        :param layer:
        :param choro_type:
        :param custom_hover:
        :param custom_zs:
        :param zmin:
        :param zmax:
        :param zoom:
        :param show_layer_elevs:
        :param show_mounding:
        :param hover_heads:
        :param hover_ks:
        :param locs:
        :param kwargs:
        :return:
        """
        return Choro(
            model=self,
            kstpkper=kstpkper,
            per=per,
            layer=layer,
            choro_type=choro_type,
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
            rch_scale=rch_scale,
            bgs=bgs,
            **kwargs
        )

    def xsect(
            self,
            per: int = None,
            kstpkper: tuple = None,
            layer: int = 0,
            cells: int | list[int] = None,
            x_or_y: str = None,
            spacing: int = 10,
            num_points: int = 100,
            extrapolate_beyond_section_ends: bool = False,
            interpolate: bool = False,
            use_rbf: bool = False
    ):
        """
        Use to plot a cross-section of heads through a model. Can be used to create an animation
        of head changes for all stress periods. The cross-section line can be defined by providing
        one cell (the 'cells' parameter) or as two ends by providing two cells to the 'cells'
        parameter. If just one cell is given, 'x_or_y' parameter defines whether the cross-section
        is vertical (along 'y' axis) or horizontal (along 'x' axis).

        Examples:

            Show an animated cross-section of all stress periods:

                XSection(model, cells=[1653, 651, 1241]).ani.show() ...OR...
                XSection(model, cells=69, layer=2, x_or_y='y').ani.show()

            Show just a cross-section of one stress period, no animation:

                XSection(model, cells=[1653, 651, 1241], kstpkper=(9, 50)).show()

        :param model: model (SimulationBase object) instance
        :param per: stress period number, O-based index; will take presedence over kstpkper if provided
        :param kstpkper: defaults to the first model stress period if not provided
        :param layer: defaults to 0
        :param cells: defines cross-section location. Can provide any number of cells
        :param x_or_y: only used if one cell is given, defines whether
        the cross-section is vertical (along 'y' axis) or horizontal (along 'x' axis).
        :param spacing: x distance between points on the plot
        :param num_points: number of points in the cross-section plot
        :param extrapolate_beyond_section_ends: not implemented
        :param surf_type: can be hds (default) or lyr (for model layers)

        """
        return XSection(
            model=self,
            per=per,
            kstpkper=kstpkper,
            layer=layer,
            cells=cells,
            x_or_y=x_or_y,
            spacing=spacing,
            num_points=num_points,
            extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
            interpolate=interpolate,
            use_rbf=use_rbf
        )

    @property
    def inputs(self):
        return Inputs(self)

    @property
    def lak(self):
        return LakOutputData(self)

    @property
    def sfr(self):
        return SFROutputData(self)

    def bud(self, package: str = None):
        if package is None:
            return Budget(self)
        else:
            return Budget(self, package)

    @property
    def master_celld(self):
        return self._master_celld

    @property
    def kstpkper(self):
        if self._kstpkper is None:
            kstpkper = self.hds.kstpkper
            self._kstpkper = kstpkper
        return self._kstpkper

    def run_simulation(self):
        # Write the datasets
        self.sim.write_simulation()
        # Save the model object to a .model file
        model_file_path = self.model_output_folder_path / f'{self.name}.model'
        with open(model_file_path, 'wb') as file:
            pickle.dump(self, file)

        # Run the simulation
        success, buff = self.sim.run_simulation(silent=False, report=True)
        print("\nSuccess is: ", success)

    """def plot_hds(self, kstpkper, zoom=13, plot_mounding=False, layer=0, zmin=None, zmax=None):
        layer_nums = self.vor.gdf_topbtm.columns[2:].to_list()
        hover = {"": ["" for cell in range(self.vor.ncpl)]}
        hover.update({
            f'Top of Model': self.vor.gdf_topbtm.loc[:, 0].to_list()
        })
        hover.update({
            f'Layer {layer} Bottom': self.vor.gdf_topbtm.loc[:, layer_num].to_list() for layer_num in layer_nums
        })

        self.hds.plot_choropleth(
            kstpkper=kstpkper,
            zoom=zoom,
            plot_mounding=plot_mounding,
            custom_hover=hover,
            layer=layer,
            zmax=zmax,
            zmin=zmin
        )"""


class OutputControl:

    def __init__(
            self,
            model: SimulationBase,
            save_record=(("HEAD", "LAST"), ("BUDGET", "LAST")),
            print_record=None
    ):
        head_file = f"{model.name}.hds"
        budget_file = f"{model.name}.cbc"
        self.oc = flopy.mf6.modflow.ModflowGwfoc(
            model.gwf,
            pname="oc",
            filename=f"{model.name}.oc",
            saverecord=save_record,
            head_filerecord=head_file,
            budget_filerecord=budget_file,
            printrecord=print_record,
        )


class DisuGrid:

    def __init__(
            self,
            vor: Vor,
            model: SimulationBase,
            top=None,
            bottom=None
    ):
        self.nlay = 1
        if top is None:
            try:
                top = vor.gdf_topbtm["top"].to_list()
            except ValueError:
                print('no top in voronoi grid. cannot find top')
        if bottom is None:
            try:
                bottom = vor.gdf_topbtm["bottom"].to_list()
            except ValueError:
                print('no bottom in voronoi grid. cannot find bottom')
        grid_props = vor.get_disv_gridprops()
        self.disu = flopy.mf6.ModflowGwfdisu(
            model.gwf,
            vertices=grid_props["vertices"],
            cell2d=grid_props["cell2d"],
            length_units="FEET",
            top=top,
            bot=bottom,
            filename=f"{model.name}.disu",
            nvert=len(grid_props["vertices"]),
            nodes=len(grid_props["cell2d"]),
            nja=vor.nja,
            iac=vor.iac,
            ja=vor.ja,
            area=vor.get_cell_areas(),
            ihc=1,
            cl12=vor.cl12,
            hwva=vor.hwva,
            idomain=[1 for i in range(vor.ncpl)],
        )


class DisvGrid:

    def __init__(
            self,
            vor: Vor = None,
            model: SimulationBase = None,
            top=None,
            bottom=None,
            nlay=1,
            idomain=None
    ):
        self.nlay = nlay
        model.nlay = nlay
        vor = model.vor if vor is None else vor
        grid_props = vor.get_disv_gridprops()
        self.disv = flopy.mf6.ModflowGwfdisv(
            model.gwf,
            length_units="FEET",
            nlay=nlay,
            ncpl=grid_props['ncpl'],
            nvert=len(grid_props["vertices"]),
            vertices=grid_props['vertices'],
            cell2d=grid_props['cell2d'],
            pname='disv',
            filename=f'{model.name}.disv',
            top=top,
            botm=bottom,
            idomain=idomain
        )


class InitialConditions:

    def __init__(
            self,
            model: SimulationBase,
            vor: Vor,
            botm_cells: list = None,
            initial_sat_thickness: float = 0.5,
            nlay=1,
            strt=None
    ):
        if botm_cells is None:
            botm_cells = [0 for cell in range(vor.ncpl * nlay)]
        if strt is None:
            strt = [cell_elev + initial_sat_thickness for cell_elev in botm_cells]
        self.ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(
            model.gwf,
            pname="ic",
            strt=strt,
            filename=f"{model.name}.ic",
        )


class TemporalDiscretization:

    def __init__(
            self,
            model: SimulationBase,
            time_units: str = 'DAYS',
            per_len: int = 1,
            period_data: list = None,
            num_steps=10,
            multiplier=1.1
    ):
        nper = model.nper
        if period_data is None:
            period_data = [[per_len, num_steps, multiplier] for per in range(nper)]
        model.num_steps = num_steps
        model.per_len = per_len
        self.tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
            model.sim,
            pname="tdis",
            time_units=time_units,
            nper=nper,
            perioddata=period_data,
            filename=f"{model.name}.tdis"
        )


class KFlow:

    def __init__(
            self,
            model: SimulationBase,
            k: list = None,
            k33_vert=None,
            perched: bool = False

    ):
        self.npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            model.gwf,
            pname="npf",
            icelltype=1,
            k=k,
            perched=perched,
            k33=k33_vert,
            save_flows=True,
            filename=f"{model.name}.npf",
            # perched=True,
        )


class Storage:

    def __init__(
            self,
            model: SimulationBase,
            specific_storage: float = 0.0001,
            specific_yield: float = 0.2,
            sto_steady: dict = None,
            sto_transient: dict = None,

    ):
        if sto_steady is None and sto_transient is None:
            sto_steady = {0: True}
        if sto_transient is None:
            sto_transient = {1: True}
        self.sto = flopy.mf6.ModflowGwfsto(
            model.gwf,
            pname="sto",
            filename=f"{model.name}.sto",
            save_flows=True,
            iconvert=1,
            ss=specific_storage,
            sy=specific_yield,
            steady_state=sto_steady,
            transient=sto_transient,
        )


class Recharge:

    def __init__(
            self,
            model: SimulationBase,
            vor: Vor = None,
            rch_dict: dict = None
    ):
        vor = model.vor if vor is None else vor
        print(vor.ncpl)
        self.rch = flopy.mf6.ModflowGwfrch(
            model.gwf,
            pname="rch",
            print_input=False,
            print_flows=False,
            save_flows=True,
            maxbound=len(vor.iverts),
            stress_period_data=rch_dict,
            filename=f"{model.name}.rch"
        )


class Drains:

    def __init__(
            self,
            model: SimulationBase,
            stress_period_data: list
    ):

        self.drn = flopy.mf6.ModflowGwfdrn(
            model=model.gwf,
            pname="drn",
            filename=f"{model.name}.drn",
            save_flows=True,
            print_flows=False,
            print_input=False,
            stress_period_data=stress_period_data,
        )


class GHB:

    def __init__(
            self,
            model: SimulationBase,
            stress_period_data
    ):
        self.ghb = flopy.mf6.ModflowGwfghb(
            model=model.gwf,
            print_input=False,
            print_flows=False,
            save_flows=True,
            filename=f"{model.name}.ghb",
            pname='ghb',
            stress_period_data=stress_period_data
        )


class CHD:

    def __init__(
            self,
            model: SimulationBase,
            stress_period_data
    ):
        self.chd = flopy.mf6.ModflowGwfchd(
            model=model.gwf,
            print_input=False,
            print_flows=False,
            save_flows=True,
            filename=f"{model.name}.chd",
            pname='chd',
            stress_period_data=stress_period_data
        )


class LAK:

    def __init__(
            self,
            model: SimulationBase,
            nlakes: int = 1,
            noutlets: int = 0,
            ntables: int = 1,
            packagedata=None,
            connectiondata=None,
            tables=None,
            outlets=None,
            perioddata=None,
            print_input=False,
            print_flows=False,
            print_stage=True,
            mover=False
    ):
        self.lak = flopy.mf6.ModflowGwflak(
            model=model.gwf,
            print_input=print_input,
            print_flows=print_flows,
            print_stage=print_stage,
            save_flows=True,
            stage_filerecord=f'{model.name}_stage.lak',
            budget_filerecord=f'{model.name}_budget.lak',
            budgetcsv_filerecord=f'{model.name}_lake_budget.csv',
            package_convergence_filerecord=f'{model.name}_lake_convergence.csv',
            mover=mover,
            surfdep=0,
            time_conversion=86_400.0,  #  assumes model time units are DAYS
            length_conversion=3.28081,  #  assumes model length units are FEET
            nlakes=nlakes,
            noutlets=noutlets,
            ntables=ntables,
            packagedata=packagedata,
            connectiondata=connectiondata,
            tables=tables,
            outlets=outlets,
            perioddata=perioddata,
            filename=f'{model.name}.lak',
            pname='lak',
            maximum_iterations=10000,
            maximum_stage_change=0.001,
        )


class UZF:

    def __init__(
            self,
            model: SimulationBase,
            packagedata=None,
            perioddata=None,
            print_input=False,
            print_flows=True,
            save_flows=True,
            mover=False,
            simulate_et=False,
            linear_gwet=False,
            square_gwet=False,
            simulate_gwseep=False,
            unsat_etwc=False,
            unsat_etae=False,
            nuzfcells=None,
            ntrailwaves=7,
            nwavesets=40,
    ):

        if nuzfcells is None:
            # sets num uzf cells to num model active cells
            nuzfcells = int(np.bincount(model.modelgrid.idomain[0])[1])

        self.uzf = flopy.mf6.ModflowGwfuzf(
            model=model.gwf,
            print_input=print_input,
            print_flows=print_flows,
            save_flows=save_flows,
            budget_filerecord=f'{model.name}_budget.uzf',
            budgetcsv_filerecord=f'{model.name}_uzf_budget.csv',
            package_convergence_filerecord=f'{model.name}_uzf_package_convergence.csv',
            mover=mover,
            simulate_et=simulate_et,
            linear_gwet=linear_gwet,
            square_gwet=square_gwet,
            simulate_gwseep=simulate_gwseep,
            unsat_etwc=unsat_etwc,
            unsat_etae=unsat_etae,
            nuzfcells=nuzfcells,
            ntrailwaves=ntrailwaves,
            nwavesets=nwavesets,
            packagedata=packagedata,
            perioddata=perioddata,
            filename=f'{model.name}.uzf',
            pname='uzf',
        )
