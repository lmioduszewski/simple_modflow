"""An SFR reach connected to NO cell (MF6 ``NONE``) is a row with a missing cell.

``split_cellid_columns`` reads MF6's ``NONE`` -- FloPy's string ``'none'`` -- as a
missing ``layer`` and ``cell`` since ledger 419, and every SFR table then cast both
with ``astype(int)``, so the reach table, the input table, the stage table and the
long profile raised ``IntCastingNaNError`` and every SFR view with them. The columns
are nullable integers now, and the reach keeps its row.

The model is a real MF6 run on a 2 x 3 x 4 DIS: three reaches, the third on no cell.
"""

from __future__ import annotations

from pathlib import Path

import flopy
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6 import package_budget, package_tables
from myflopy.modflow.mf6.package_explorer_utils import CellNumbering, split_cellid_columns


@pytest.fixture(scope="module")
def run(tmp_path_factory) -> Path:
    workspace = tmp_path_factory.mktemp("sfr_none")
    sim = flopy.mf6.MFSimulation(sim_name="sim", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim, complexity="COMPLEX")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True, newtonoptions="NEWTON")
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=2, nrow=3, ncol=4, delr=10.0, delc=10.0, top=10.0, botm=[0.0, -10.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=[1, 0])
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[((0, r, 0), 5.0) for r in range(3)] + [((0, r, 3), 1.0) for r in range(3)]
    )
    cells = [(0, 1, 1), (0, 1, 2), "none"]
    packagedata = [
        [i, cell, 10.0, 1.0, 0.001, 4.0 - i * 0.5, 1.0, 1.0, 0.03, 1 if i in (0, 2) else 2, 1.0, 0]
        for i, cell in enumerate(cells)
    ]
    flopy.mf6.ModflowGwfsfr(
        gwf, nreaches=3, packagedata=packagedata, connectiondata=[[0, -1], [1, 0, -2], [2, 1]],
        perioddata={0: [[0, "inflow", 10.0]]}, stage_filerecord="m.sfr.stg",
        budget_filerecord="m.sfr.cbc", save_flows=True,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    success, report = sim.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])
    assert "  3  none  " in (workspace / "m.sfr").read_text()
    return workspace


def test_the_reach_tables_keep_a_reach_on_no_cell(run):
    model = mf.load_mf6_run(str(run))
    reaches = package_tables.build_sfr_reach_table(model).set_index("reach")
    assert list(reaches.index) == [0, 1, 2]
    assert reaches.loc[0, ["layer", "cell"]].tolist() == [0, 5]  # (0, 1, 1) on a 4-column DIS
    assert pd.isna(reaches.loc[2, "layer"]) and pd.isna(reaches.loc[2, "cell"])
    assert str(reaches["cell"].dtype) == "Int64"

    inputs = package_tables.build_sfr_input_table(model)
    assert sorted(inputs["reach"].unique()) == [0, 1, 2]
    assert inputs.loc[inputs["reach"] == 2, "cell"].isna().all()

    stage = package_budget.build_sfr_stage_result_table(model)
    assert sorted(stage["reach"].unique()) == [0, 1, 2]
    profile = package_budget.build_sfr_long_profile_table(model, per=0)
    assert profile.set_index("reach").loc[1, "cell"] == 6


def test_the_views_read_the_model(run):
    model = mf.load_mf6_run(str(run))
    assert len(model.sfr.results.q.get()) == 4  # the two reaches ON a cell, two periods
    assert sorted(model.sfr.results.profile.get()["reach"]) == [0, 1, 2]
    assert model.sfr.inputs.rtp.get()["cell"].isna().sum() == 2  # reach 2, both periods


def test_a_table_with_every_cell_keeps_plain_integers():
    numbering = CellNumbering("structured", nlay=1, ncpl=4, ncol=2)
    frame = split_cellid_columns(pd.DataFrame({"cellid": [(0, 0, 1), (0, 1, 1)]}), numbering)
    assert frame["cell"].tolist() == [1, 3] and frame["cell"].dtype == "int64"
    missing = split_cellid_columns(pd.DataFrame({"cellid": [(0, 1, 1), "none", None]}), numbering)
    assert missing["cell"].dtype == "Int64" and missing["cell"].isna().tolist() == [False, True, True]
