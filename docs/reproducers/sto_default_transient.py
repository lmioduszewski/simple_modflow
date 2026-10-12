"""Which stress periods MF6 runs steady state, and what ``model.calendar`` reads (ledger 529).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/sto_default_transient.py

Needs `mf6` on PATH (measured with 6.8.0.dev0, FloPy 3.10). About 2 s.

A 1 x 3 unconfined FloPy DIS model (one CHD cell, recharge, three 10-day periods),
run with four STO declarations. MF6's storage flow per period, from its budget file,
is the measure: a steady-state period carries none.

    no STO                          storage: no records         -> steady, steady, steady
    STO, TRANSIENT from period 2    storage: >0, >0, >0         -> transient x3
    STO, TRANSIENT from period 1    storage: identical to above -> transient x3
    STO, STEADY 1 then TRANSIENT 2  storage: 0, >0, >0          -> steady, transient x2

So a period BEFORE the first STO block runs TRANSIENT -- MF6's listing names only
the periods with a block ("STRESS PERIOD 2 IS TRANSIENT"), so the listing alone does
not show it; the storage flows do, equal to the bit to the declaration that names
period 1. ``period_calendar.steady_periods`` reads the four the same way, from the
run's FILES (``mf.load_mf6_run``).
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile

import flopy
import numpy as np

import myflopy as mf


def _run(workspace: pathlib.Path, sto: dict | None) -> pathlib.Path:
    sim = flopy.mf6.MFSimulation(sim_ws=str(workspace), exe_name="mf6", verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(10.0, 1, 1.0)] * 3, time_units="days",
                          start_date_time="2016-09-30")
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nrow=1, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, 0, 0), 6.0)])
    flopy.mf6.ModflowGwfrcha(gwf, recharge=1e-3)
    if sto is not None:
        flopy.mf6.ModflowGwfsto(gwf, sy=0.2, ss=1e-5, iconvert=1, **sto)
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="m.cbc", head_filerecord="m.hds",
                           saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success, workspace
    return workspace


def _storage(workspace: pathlib.Path) -> list[float] | str:
    budget = flopy.utils.CellBudgetFile(workspace / "m.cbc", precision="double")
    names = [name.decode().strip() for name in budget.get_unique_record_names()]
    if "STO-SY" not in names:
        return "no records"
    return [
        float(sum(np.abs(np.asarray(budget.get_data(text=text, kstpkper=key)[0], dtype=float)).sum()
                  for text in ("STO-SY", "STO-SS")))
        for key in budget.get_kstpkper()
    ]


if __name__ == "__main__":
    root = pathlib.Path(tempfile.mkdtemp(prefix="sto_default_transient_"))
    try:
        cases = {
            "no STO": None,
            "TRANSIENT from period 2": {"transient": {1: True}},
            "TRANSIENT from period 1": {"transient": {0: True}},
            "STEADY 1, TRANSIENT 2": {"steady_state": {0: True}, "transient": {1: True}},
        }
        for index, (label, sto) in enumerate(cases.items()):
            workspace = _run(root / f"case{index}", sto)
            listing = (workspace / "m.lst").read_text().splitlines()
            said = [line.strip() for line in listing
                    if line.strip().startswith("STRESS PERIOD") and " IS " in line]
            flags = mf.load_mf6_run(workspace).calendar.steady
            print(f"{label:<26} storage={_storage(workspace)}  listing={said}  calendar.steady={flags}")
    finally:
        shutil.rmtree(root, ignore_errors=True)
