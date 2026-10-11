"""What MF6 runs when a split simulation is run before it is written.

Ledger 398. FloPy's ``Mf6Splitter`` builds the split simulation with the SOURCE
simulation's ``sim_path``; only ``ParallelSplitRun.write`` moves it to the split
workspace. ``ParallelModelWorkflow.run(write=False)`` and
``ParallelSplitRun.run_serial/run(write=False)`` before ``write()`` used to call
``run_simulation`` on it anyway. This runs FloPy's own ``run_simulation`` on such an
unwritten split -- bypassing myflopy's refusal -- and reports what it did.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/split_run_before_write.py

THE SETUP. A 1-layer 3x6 DIS GWF model between two CHDs, written and run in
``<tmp>/source``; split into two 3x3 halves with ``split_model(write=False)``.

MEASURED (FloPy 3.10, MF6 6.8.0.dev0, 2026-10-11):
  the unwritten split's sim_path  == the source directory
  run_simulation()                -> success True, "Normal termination"
  files rewritten in the source   gwf.cbc, gwf.dis.grb, gwf.hds, gwf.lst, mfsim.lst
  models in the run's mfsim.lst   ``gwf`` only -- the SOURCE model; no gwf_0/gwf_1
  the split workspace             never created
So the old ``ParallelModelWorkflow.run(write=False)`` returned normally, having re-run
the source model over its own outputs and run no part of the split. Today myflopy
refuses: ``ValueError`` from the workflow, ``RuntimeError`` from the split run, printed
at the end.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace

import flopy
import numpy as np

from myflopy.modflow.mf6.parallel import ParallelModelWorkflow


def _source(workspace: Path) -> SimpleNamespace:
    sim = flopy.mf6.MFSimulation(sim_name="src", sim_ws=workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    ims = flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=6, delr=1.0, delc=1.0, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, 0, 0), 6.0), ((0, 2, 5), 4.0)])
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success
    return SimpleNamespace(sim=sim, gwf=gwf)


def _state(directory: Path) -> dict[str, tuple[int, bytes]]:
    return {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in directory.iterdir() if p.is_file()}


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = _source(root / "source")
        before = _state(root / "source")
        mask = np.array([[0, 0, 0, 1, 1, 1]] * 3)
        split = ParallelModelWorkflow(source).split_model(workspace=root / "split", mask=mask, write=False)
        print("unwritten split sim_path:", split.simulation.sim_path)
        print("is the source directory: ", Path(split.simulation.sim_path).resolve() == (root / "source").resolve())

        # FloPy's own run, as the old write=False path made it
        success, report = split.simulation.run_simulation(silent=True, report=True)
        after = _state(root / "source")
        print("run_simulation success:  ", success, "|", report[-1].strip() if report else "")
        print("rewritten in the source: ", sorted(k for k in after if before.get(k) != after[k]))
        listing = (root / "source" / "mfsim.lst").read_text(errors="replace")
        names = sorted({n for n in ("gwf_0", "gwf_1") if n.upper() in listing.upper()})
        print("split models in mfsim.lst:", names or "none (the source model ran)")
        print("split workspace exists:  ", (root / "split").exists())

        for call in (
            lambda: ParallelModelWorkflow(source).run(workspace=root / "w", mask=mask, write=False),
            lambda: split.run_serial(write=False, silent=True),
        ):
            try:
                call()
            except (ValueError, RuntimeError) as error:
                print(f"today: {type(error).__name__}: {str(error)[:110]}...")


if __name__ == "__main__":
    main()
