"""Two records on one cell against the ONE record a USG export writes for them, through MF6.

Ledger 493 (amending ledger 467). Run it::

    python docs/reproducers/usg_export_shared_cell_elevation.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given (a
temporary directory, removed at the end, by default). A few seconds.

WHAT IT MEASURES
  One 10 x 10 ft cell whose head a GHB of conductance 1e9 holds at ``h``; on it a DRN,
  GHB or RIV given as two records, then as the one record the export writes for them
  (``usg/export._combine_records``: conductances add, every elevation is the
  CONDUCTANCE-weighted mean), then as the old plain mean. Each line prints MF6's flow
  through the package (MF6's sign: a drain's
  outflow is negative).

  * DRN ``elev 7.5 / cond 100`` + ``7.0 / 50``: at h = 8.0 both drain and the weighted
    record is exact (two records -100.0, weighted -100.0, mean -112.5); at h = 7.2 the
    head is BETWEEN the elevations and neither combination is (two records -10.0, both
    0.0) -- the case the export warns about.
  * GHB ``9.0 / 100`` + ``6.0 / 50``: linear in the head, so the weighted record is
    exact at every head (150.0 / 150.0 / 75.0 at h = 7.0).
  * RIV ``stage 9.0, cond 100, rbot 7.5`` + ``8.0, 50, 7.0``: exact above both rbots
    (h = 8.5: 25.0 / 25.0 / 0.0) and below both (h = 6.0: 200.0 / 200.0 / 187.5), not
    between (h = 7.2: 190.0 / 200.0 / 187.5).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import flopy
import numpy as np
import pandas as pd

from myflopy.modflow.usg.export import _combine_records

COLUMNS = {"DRN": ["elev", "cond"], "GHB": ["bhead", "cond"], "RIV": ["stage", "cond", "rbot"]}


def flow(workspace: Path, package: str, records, hold: float) -> float:
    sim = flopy.mf6.MFSimulation(sim_name="s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE", inner_dvclose=1e-10, outer_dvclose=1e-10)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=1, delr=10.0, delc=10.0, top=20.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=hold)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data={0: [[(0, 0, 0), hold, 1e9]]}, pname="hold")
    getattr(flopy.mf6, f"ModflowGwf{package.lower()}")(
        gwf, stress_period_data={0: [[(0, 0, 0), *r] for r in records]}, pname="bc",
        filename=f"m.bc.{package.lower()}",
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="m.cbc", saverecord=[("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok, workspace
    return float(flopy.utils.CellBudgetFile(str(workspace / "m.cbc")).get_data(text=package)[-1]["q"].sum())


def one_record(package: str, records, rule: str):
    columns = COLUMNS[package]
    if rule == "mean":
        values = np.mean(np.asarray(records, float), axis=0)
        values[columns.index("cond")] = sum(r[columns.index("cond")] for r in records)
        return [list(values)]
    frame = pd.DataFrame([[0, 0, 0, *r] for r in records], columns=["period", "layer", "cell", *columns])
    combined, _inexact, _unweighted = _combine_records(frame, package, columns)
    return [list(combined.loc[0, columns])]


def main(root: Path) -> None:
    cases = [
        ("DRN", [(7.5, 100.0), (7.0, 50.0)], (8.0, 7.2)),
        ("GHB", [(9.0, 100.0), (6.0, 50.0)], (7.0, 10.0)),
        ("RIV", [(9.0, 100.0, 7.5), (8.0, 50.0, 7.0)], (8.5, 6.0, 7.2)),
    ]
    for package, records, heads in cases:
        for hold in heads:
            flows = {
                label: flow(root / f"{package}_{hold}_{label}", package, recs, hold)
                for label, recs in (
                    ("two records", records),
                    ("exported (weighted)", one_record(package, records, "weighted")),
                    ("plain mean", one_record(package, records, "mean")),
                )
            }
            print(f"{package} h = {hold}: " + "; ".join(f"{k} {v:.3f}" for k, v in flows.items()))


if __name__ == "__main__":
    given = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    root = given or Path(tempfile.mkdtemp())
    try:
        main(root)
    finally:
        if given is None:
            shutil.rmtree(root, ignore_errors=True)
