from __future__ import annotations

from typing import TYPE_CHECKING

import flopy
import numpy as np

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor


class OutputControl:

    def __init__(
            self,
            model: "SimulationBase",
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


class InitialConditions:

    def __init__(
            self,
            model: "SimulationBase",
            vor: "Vor",
            botm_cells: list = None,
            initial_sat_thickness: float = 0.5,
            nlay=1,
            strt=None
    ):
        if botm_cells is None:
            botm_cells = [0 for _ in range(vor.ncpl * nlay)]
        if strt is None:
            strt = [cell_elev + initial_sat_thickness for cell_elev in botm_cells]
        self.ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(
            model.gwf,
            pname="ic",
            strt=strt,
            filename=f"{model.name}.ic",
        )


class KFlow:

    def __init__(
            self,
            model: "SimulationBase",
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
            save_saturation=True,
            save_specific_discharge=True,
            filename=f"{model.name}.npf",
        )


class Storage:

    def __init__(
            self,
            model: "SimulationBase",
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
            model: "SimulationBase",
            vor: "Vor" = None,
            rch_dict: dict = None,
            auxiliary: list[str] = None
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
            filename=f"{model.name}.rch",
            auxiliary=auxiliary
        )


class Drains:

    def __init__(
            self,
            model: "SimulationBase",
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
            model: "SimulationBase",
            stress_period_data,
            auxiliary=None
    ):
        self.ghb = flopy.mf6.ModflowGwfghb(
            model=model.gwf,
            print_input=False,
            print_flows=False,
            save_flows=True,
            filename=f"{model.name}.ghb",
            pname='ghb',
            stress_period_data=stress_period_data,
            auxiliary=auxiliary
        )


class CHD:

    def __init__(
            self,
            model: "SimulationBase",
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
            model: "SimulationBase",
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
            time_conversion=86_400.0,
            length_conversion=3.28081,
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
            maximum_iterations=100,
            maximum_stage_change=1e-5,
        )


class UZF:

    def __init__(
            self,
            model: "SimulationBase",
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
