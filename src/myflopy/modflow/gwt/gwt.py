from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

from pathlib import Path

import flopy
import numpy as np
import pandas as pd

from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg


class GWT:

    def __init__(
            self,
            gwf_model: SimulationBase,
            simulation=None,
            name: str = None,
            ic: np.array | pd.DataFrame = None,
            save_record: tuple = (("BUDGET", "LAST"), ("CONCENTRATION", "LAST")),
    ):
        self.gwf_model = gwf_model
        self.name = name
        simulation = self.gwf_model.sim if simulation is None else simulation
        self.gwt = flopy.mf6.ModflowGwt(
            modelname=self.name,
            model_nam_file=f"{self.name}.gwt.nam",
            # version='mf6',
            # exe_name='mf6',
            simulation=simulation,
            save_flows=True,
            print_flows=False
        )
        self.ic = flopy.mf6.ModflowGwtic(
            model=self.gwt,
            filename=f"{self.name}.gwt.ic",
            strt=ic
        )
        self.oc = flopy.mf6.ModflowGwtoc(
            model=self.gwt,
            filename=f"{self.name}.gwt.oc",
            concentration_filerecord=f"{self.name}.gwt.unc",
            saverecord=save_record,
            printrecord=[("BUDGET", "ALL"), ("CONCENTRATION", "LAST")],
            budget_filerecord=f"{name}.gwt.cbc",
        )


class ADV:

    def __init__(self, gwt_model: GWT, adv_scheme: str = 'upstream'):
        self.gwt_model = gwt_model
        self.adv = flopy.mf6.ModflowGwtadv(
            model=self.gwt_model.gwt,
            filename=f"{self.gwt_model.name}.gwt.adv",
            scheme=adv_scheme,
        )


class DSP:

    def __init__(
            self,
            gwt_model: GWT,
            diffc=None,
            alh=None,
            alv=None,
            ath1=None,
            ath2=None,
            atv=None
    ):
        self.gwt_model = gwt_model

        self.dsp = flopy.mf6.ModflowGwtdsp(
            model=self.gwt_model.gwt,
            filename=f"{self.gwt_model.name}.gwt.dsp",
            xt3d_off=False,
            diffc=diffc,
            alh=alh,
            alv=alv,
            ath1=ath1,
            ath2=ath2,
            atv=atv,

        )

class CNC:

    def __init__(
            self,
            gwt_model: GWT,
            stress_period_data=None,
    ):
        self.gwt_model = gwt_model

        self.cnc = flopy.mf6.ModflowGwtcnc(
            model=self.gwt_model.gwt,
            filename=f"{self.gwt_model.name}.gwt.cnc",
            stress_period_data=stress_period_data,
        )

class FMI:

    def __init__(
            self,
            gwt_model: GWT,
            packagedata=None
    ):
        self.gwt_model = gwt_model

        if packagedata is None:
            gwf = self.gwt_model.gwf_model

        self.fmi = flopy.mf6.ModflowGwtfmi(
            model=self.gwt_model.gwt,
            flow_imbalance_correction=True,
            filename=f"{self.gwt_model.name}.gwt.fmi",
            packagedata=packagedata,
        )


if __name__ == '__main__':
    import pickle

    import flopy
    with open(Path(r"C:\Users\lukem\mf6\ssb_Phs2_nPrch\ssb_Phs2_nPrch.model"), 'rb') as f:
        gwf_model: SimulationBase = pickle.load(f)
    itus_phs2_path = Path(r"C:\Users\lukem\mf6\SSB data\itus_phs2_no_reserve_UICs.gpkg")
    itus_phs2 = read_shp_gpkg(itus_phs2_path)
    uic_cells = gwf_model.vor.get_vor_cells_as_series(itus_phs2).apply(lambda x: x[0]).to_list()
    stress_period_data = {0: [[(1, cell), 1.0] for cell in uic_cells]}

    sim = flopy.mf6.MFSimulation(
        sim_name='ssb_gwt',
        sim_ws=gwf_model.model_output_folder_path / 'ssb_gwt',
    )
    tdis = flopy.mf6.ModflowTdis(
        simulation=sim,
        filename=f"{sim.name}.tdis",
        time_units="DAYS",
        nper=gwf_model.nper,
        perioddata=gwf_model.sim.tdis.perioddata.get_data().tolist()
    )
    ims = flopy.mf6.ModflowIms(
        simulation=sim,
        filename=f"{sim.name}.ims",
        complexity="COMPLEX",
        # Nonlinear controls
        under_relaxation="DBD",
        under_relaxation_theta=0.7,  # 0.7–0.8 typical
        under_relaxation_kappa=0.20,
        under_relaxation_momentum=0.0,  # disable momentum for stability
        backtracking_number=50,
        backtracking_tolerance=1.1,
        backtracking_reduction_factor=0.2,
        backtracking_residual_limit=1000,

        # Iteration budgets
        outer_maximum=300,
        inner_maximum=300,

        # Convergence criteria
        outer_dvclose=1e-3,
        inner_dvclose=1e-3,
        # rcloserecord=[1e-6, "strict"],

        # Linear solver + conditioning
        linear_acceleration="BICGSTAB",
        scaling_method="L2NORM",
        reordering_method="RCM",
        preconditioner_levels=2,
        preconditioner_drop_tolerance=0.0,

        relaxation_factor=0.97,
    )
    transport = GWT(
        gwf_model=gwf_model,
        simulation=sim,
        name='ssb_gwt',
        ic=0.5
    )
    disv = flopy.mf6.ModflowGwtdisv(
        model=transport.gwt,
        filename=f"{transport.name}.gwt.disv",
        length_units="FEET",
        nlay=7,
        ncpl=gwf_model.vor.ncpl,
        nvert=gwf_model.gwf.disv.nvert.data,
        vertices=gwf_model.gwf.disv.vertices.get_data(),
        top=gwf_model.gwf.disv.top.data,
        botm=gwf_model.gwf.disv.botm.data,
        idomain=gwf_model.gwf.disv.idomain.data,
        cell2d=gwf_model.gwf.disv.cell2d.get_data()
    )
    adv = ADV(transport)
    dsp = DSP(transport, alh=31.25, alv=31.25, ath1=3.125, atv=0.3125)
    cnc = CNC(transport, stress_period_data=stress_period_data)
    fmi = FMI(transport,
              packagedata=[['GWFBUDGET', gwf_model.gwf.output.budget().filename],
                           ['GWFHEAD', gwf_model.gwf.output.head().filename]])
    mst = flopy.mf6.ModflowGwtmst(
        transport.gwt,
        porosity=0.25,
        filename=f"{transport.name}.gwt.mst",
    )
    ssm = flopy.mf6.ModflowGwtssm(
        transport.gwt,
        filename=f"{transport.name}.gwt.ssm",


    )
    sim.write_simulation()
    sim.run_simulation()
