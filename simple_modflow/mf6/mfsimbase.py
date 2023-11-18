import flopy
from pathlib import Path
from voronoiplus import VoronoiGridPlus as Vor


class SimulationBase:

    def __init__(
            self,
            name: str = 'mf6_model',
            mf_folder_path: Path = Path().home().joinpath('mf6')
    ):
        self.name = name
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


class InitialConditions:

    def __init__(
            self,
            model: SimulationBase,
            vor: Vor,
            botm_cells: list = None,
            initial_sat_thickness: float = 0.5
    ):
        if botm_cells is None:
            botm_cells = [0 for cell in range(vor.ncpl)]
        self.ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(
            model.gwf,
            pname="ic",
            strt=[cell_elev + initial_sat_thickness for cell_elev in botm_cells],
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
            k_list: list = None,

    ):
        self.npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            model.gwf,
            pname="npf",
            icelltype=1,
            k=k_list,
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
            rch_dict: dict
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
