"""What MODFLOW 6 does with a boundary elevation below its own cell bottom.

The measurements behind ledger 291 (and a second reading of ledger 290's RIV-pair
rule), on a one-cell bench small enough to check by eye. Run it::

    python docs/reproducers/riv_pair_and_icelltype.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given.

**The bench.** One layer, DISV, a 2x2 plan of 100 ft squares. ``top = 60.0``,
``botm = 48.0``. Cell 0 is a CHD at 55.0 -- the driver -- and cell 3 carries the
package under test. Steady state, ``k = 10``.

**What it shows.**

1. ``icelltype`` decides everything. At 0 -- what ``mf.npf(k=...)`` writes, because it
   names no icelltype at all -- MODFLOW 6 ACCEPTS a below-bottom DRN, GHB, CHD or RIV,
   terminates normally, and over-discharges: the drain at 45.0 pulls -545.45 against
   -272.73 for the same drain at 50.0, **2.00x**. At 1 every one of them is REFUSED by
   name, and NEWTON does not change it. The MODFLOW-USG importer carries ``icelltype``
   straight from the LPF ``LAYTYP``, so a converted model lands on either side of this
   split -- which is why its ``validate()`` reads the severity from the model rather
   than asserting it.

2. **A river is a pair.** MODFLOW 6 compares ``rbot`` against the cell bottom and
   enforces ``stage >= rbot`` separately. With a declared stage 40 / rbot 35 -- BOTH
   below a cell bottomed at 48 -- as declared q = -818.18 (refused at icelltype 1);
   stage lifted alone -381.82 (**still refused at icelltype 1**: ``rbot`` is untouched,
   and ``rbot`` is what MODFLOW 6 looks at); both clamped to the floor -381.82 with the
   declared 5 ft of depth collapsed to zero; the stage CARRIED up by the bed's delta
   -109.09, and it runs at both. Four readings of one river, 7.5x apart end to end.

3. **But carrying the stage is right only when the stage is itself illegal.** Stage 55
   with rbot 40 in the same cell is an ordinary river -- water surface well inside the
   cell, bed below it. Nothing about the stage is wrong, so a check that reads the first
   value in the record reports nothing, and MODFLOW 6 refuses it at icelltype 1 all the
   same. Repair it by carrying the stage up with the bed and the surface lands at 63 --
   **three feet above land surface** -- turning q = 0.00 into **+436.37**; with a deeper
   bed (stage 50 / rbot 20) the sign reverses, -272.73 to +1254.55. Lifting the bed
   alone reproduces the declared answer exactly and runs at both icelltypes. Ledger 290
   measured only case 2 and generalised; ledger 291 measured this one.

One detail worth keeping: a DISV ``cell2d`` vertex ring must be CLOCKWISE. Counter-
clockwise gets "Calculated CELL2D area less than zero", which is not what it sounds like.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import flopy
import numpy as np

TOP, BOT = 60.0, 48.0
CELL = 3  # the cell carrying the package under test; cell 0 drives the model

VERTICES = [[j * 3 + i, i * 100.0, j * 100.0] for j in range(3) for i in range(3)]
CELL2D = []
for cell, (row, col) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
    lower_left = row * 3 + col
    ring = [lower_left, lower_left + 3, lower_left + 4, lower_left + 1]  # clockwise
    CELL2D.append([cell, (col + 0.5) * 100.0, (row + 0.5) * 100.0, 4, *ring])

PACKAGES = {
    "drn": flopy.mf6.ModflowGwfdrn,
    "ghb": flopy.mf6.ModflowGwfghb,
    "riv": flopy.mf6.ModflowGwfriv,
    "chd": flopy.mf6.ModflowGwfchd,
}


def run(workspace: Path, tag: str, *, icelltype: int, package: str, values, newton=False):
    """Build and run one bench case; return its verdict."""

    ws = workspace / tag
    if ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True)

    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="COMPLEX", outer_maximum=200, inner_maximum=200)
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="t", save_flows=True,
        newtonoptions="NEWTON UNDER_RELAXATION" if newton else None,
    )
    flopy.mf6.ModflowGwfdisv(
        gwf, nlay=1, ncpl=4, nvert=len(VERTICES), vertices=VERTICES, cell2d=CELL2D,
        top=TOP, botm=[BOT],
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=10.0, icelltype=icelltype, save_flows=True)
    flopy.mf6.ModflowGwfic(gwf, strt=55.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0), 55.0]]}, pname="chd")
    PACKAGES[package](
        gwf, stress_period_data={0: [[(0, CELL), *values]]}, pname=f"{package}x", save_flows=True
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord="t.hds", budget_filerecord="t.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    ok, buff = sim.run_simulation(silent=True, report=True)

    if not ok:
        listing = (ws / "t.lst").read_text() if (ws / "t.lst").is_file() else ""
        blame = [
            line.strip()
            for line in ("\n".join(buff) + "\n" + listing).splitlines()
            if "LESS THAN" in line.upper()
        ]
        print(f"{tag:<34} REFUSED  {blame[0] if blame else '(see listing)'}")
        return None

    budget = flopy.utils.CellBudgetFile(str(ws / "t.cbc"))
    flux = 0.0
    for name in budget.get_unique_record_names():
        if name.decode().strip().upper() != package.upper():
            continue
        # Sum only the cell under test. The driver is a CHD, so in the `chd` case TWO
        # packages write under the text "CHD" and `get_data(text=...)[0]` returns the
        # DRIVER's block -- reporting +1200 (its inflow) for a boundary that is in fact
        # pulling -1200 out. Node numbers are 1-based; the cell under test is 3.
        for block in budget.get_data(text=name):
            flux += float(np.sum([row[2] for row in block if int(row[0]) == CELL + 1]))
    print(f"{tag:<34} RUNS     q = {flux:>11.4f}")
    return flux


def main(workspace: Path) -> None:
    print(f"one DISV cell, top {TOP}, bottom {BOT}\n")
    for icelltype in (0, 1):
        print(f"--- icelltype = {icelltype} " + "-" * 40)
        inside = run(workspace, f"drn_inside_ict{icelltype}", icelltype=icelltype,
                     package="drn", values=[50.0, 100.0])
        below = run(workspace, f"drn_below_ict{icelltype}", icelltype=icelltype,
                    package="drn", values=[45.0, 100.0])
        if inside and below:
            print(f"{'':<34}          ^ {below / inside:.2f}x the drain inside its cell")
        run(workspace, f"ghb_below_ict{icelltype}", icelltype=icelltype,
            package="ghb", values=[45.0, 100.0])
        run(workspace, f"chd_below_ict{icelltype}", icelltype=icelltype,
            package="chd", values=[45.0])

        # One river, four readings. Declared stage 40 / rbot 35 -- BOTH below the cell
        # bottom. The lifts put the floored number at 48.000048, = 48 + the clearance.
        run(workspace, f"riv_both_below_as_declared_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[40.0, 100.0, 35.0])
        run(workspace, f"riv_both_below_stage_only_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[48.000048, 100.0, 35.0])
        run(workspace, f"riv_both_below_clamped_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[48.000048, 100.0, 48.000048])
        run(workspace, f"riv_both_below_carried_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[53.000048, 100.0, 48.000048])

        # The ORDINARY river: surface inside the cell, bed below it. A stage-only check
        # is blind to it, and carrying the stage up with the bed wrecks it.
        run(workspace, f"riv_ordinary_as_declared_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[55.0, 100.0, 40.0])
        run(workspace, f"riv_ordinary_stage_carried_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[63.000048, 100.0, 48.000048])
        run(workspace, f"riv_ordinary_bed_only_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[55.0, 100.0, 48.000048])
        # And the inversion, which MODFLOW 6 refuses whatever icelltype is.
        run(workspace, f"riv_stage_under_bed_ict{icelltype}", icelltype=icelltype,
            package="riv", values=[50.0, 100.0, 55.0])
        print()

    print("--- icelltype = 1, NEWTON " + "-" * 34)
    run(workspace, "drn_below_newton", icelltype=1, package="drn",
        values=[45.0, 100.0], newton=True)
    run(workspace, "riv_stage_lifted_newton", icelltype=1, package="riv",
        values=[48.000048, 100.0, 35.0], newton=True)


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "riv_pair_and_icelltype_ws"))
