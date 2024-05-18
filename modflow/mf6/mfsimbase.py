import flopy
from pathlib import Path
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor


class SimulationBase:

    def __init__(
            self,
            name: str = 'mf6_model',
            mf_folder_path: Path = Path().home().joinpath('mf6'),
    ):
        self.name = name
        self.nlay = None
        self.model_output_folder_path = mf_folder_path.joinpath(f'{name}')
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
            newtonoptions="under_relaxation",
            print_flows=True,
            save_flows=True
        )
        self.ims = flopy.mf6.modflow.mfims.ModflowIms(
            self.sim,
            pname="ims",
            complexity="COMPLEX",
            under_relaxation="DBD",
            under_relaxation_theta=0.72,
            under_relaxation_kappa=0.24,
            under_relaxation_momentum=0.0001,
            backtracking_number=20,
            backtracking_tolerance=100,
            backtracking_reduction_factor=0.3,
            backtracking_residual_limit=100,
        )

    def run_simulation(
            self,
    ):
        # Write the datasets
        self.sim.write_simulation()
        # Run the simulation
        success, buff = self.sim.run_simulation()
        print("\nSuccess is: ", success)


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
            vor: Vor,
            model: SimulationBase,
            top=None,
            bottom=None,
            nlay=1,
    ):
        self.nlay = nlay
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
            botm=bottom
        )


class InitialConditions:

    def __init__(
            self,
            model: SimulationBase,
            vor: Vor,
            botm_cells: list = None,
            initial_sat_thickness: float = 0.5,
            nlay = 1,
            strt = None
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
            nper: int = 1,
            per_len: int = 1,
            period_data: list = None
    ):
        if period_data is None:
            period_data = [[per_len, 20, 1.2] for per in range(nper)]
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

    ):
        self.npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            model.gwf,
            pname="npf",
            icelltype=1,
            k=k,
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
        if sto_steady is None:
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
            vor: Vor,
            rch_dict: dict = None
    ):
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
            packagedata = None,
            connectiondata = None,
            tables = None,
            outlets = None,
            perioddata = None
    ):

        self.lak = flopy.mf6.ModflowGwflak(
            model=model.gwf,
            print_input=False,
            print_flows=False,
            print_stage=False,
            save_flows=True,
            stage_filerecord=f'{model.name}_stage.lak',
            budget_filerecord=f'{model.name}_budget.lak',
            budgetcsv_filerecord=f'{model.name}_lake_budget.csv',
            package_convergence_filerecord=f'{model.name}_lake_convergence.csv',
            mover=True,
            surfdep=0,
            time_conversion=86_400.0, #  assumes model time units are DAYS
            length_conversion=3.28081, #  assumes model length units are FEET
            nlakes=nlakes,
            noutlets=noutlets,
            ntables=ntables,
            packagedata=packagedata,
            connectiondata=connectiondata,
            tables=tables,
            outlets=outlets,
            perioddata=perioddata,
            filename=f'{model.name}_lak',
            pname='lak'
        )