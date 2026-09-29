"""A flow target is compared in MF6's own sign, on both routes (ledger 350).

``SfrFlowTargets`` reads a reach's streamflow on one of two routes: MF6's own
``DOWNSTREAM-FLOW`` observation when the SFR package's OBS declares the target, else
the SFR budget file (outbound ``FLOW-JA-FACE`` + ``TO-MVR`` + ``EXT-OUTFLOW``).
``DrnFlowTargets`` likewise reads MF6's ``DRN`` observations or the cell-budget file.
This script measures which sign each route hands on, what a gauged flow entered either
way scores, and what a PEST build writes.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/flow_target_sign.py

It uses the public API only, so it measures whichever tree is on ``PYTHONPATH`` --
the "before" column below is this script run against the parent commit's ``src``
(``git archive 80af79d src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``).
Needs ``mf6``; section C needs pyEMU.

THE SETUP. The toy of ``tests/test_flow_targets_period_end.py``: a 6 x 6 irregular
Voronoi grid of ~100 ft cells, one confined layer, CHD at 102 ft west and east, heads
filling from 90 ft; a four-reach stream (reach 3 the outlet), a three-cell DRN zone and
a two-cell lake. Three 10-day periods of three steps, ``SAVE ALL``. SFR
``DOWNSTREAM-FLOW`` is declared for every reach under ``r0`` .. ``r3`` and DRN for the
zone under ``seep``; the same features under ``g0`` .. ``g3`` / ``zone`` read the
budget files. Two runs:

  ``plain``  as above (the scoping toy of the finding, F3).
  ``mover``  plus an MVR: reach 1 passes 30% of its flow to the lake (``FACTOR 0.3``)
             and the outlet up to 20 (``UPTO 20``) -- the quantity then has all three
             budget pieces (``tests/test_flow_target_sign.py``'s toy).

  A. Every reach on both routes, at each period end: the two numbers and their largest
     relative difference; the budget pieces of the mover reaches as MF6 wrote them.
  B. ``compare()`` given the model's own flow as a gauged ``Q``, entered ``+Q`` and
     ``-Q``, on both routes, SFR (reach 0 of ``plain``) and DRN; phi = sum of squared
     residuals.
  C. A PEST build (noptmax 0) of the SFR observation-route target entered ``+Q``, and of
     every route entered ``-Q``: phi of the control file's ``obsval`` against the base
     run's table as PEST's instruction file reads it.

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, 2026-09-29)                before           after
  A  plain, reach 0: observation route, periods 0/1/2    -165.733, -81.129,   the same
                                                         -218.466
     plain, reach 0: budget route, periods 0/1/2         +165.733, +81.129,   -165.733, -81.129,
                                                         +218.466             -218.466
     max |obs - budget| / |obs|, 4 reaches, each run     2 (plain, mover)     0 (plain, mover)
     mover, period 0 end: reach 1 FLOW-JA-FACE / TO-MVR  -34.836 / -14.926    (MF6's files)
     mover, period 0 end: outlet EXT-OUTFLOW / TO-MVR    -83.269 / -20.000    (MF6's files)
  B  SFR +Q, observation route: residual / phi           -2Q / 327,106.2      REFUSED at
     SFR +Q, budget route: residual                      0                    construction,
     DRN +Q, either route: residual / phi                -2Q / 7,916,474.4    naming -Q
     SFR -Q, observation / budget route: residual        0 / +2Q (phi         0 / 0
                                                         327,106.2)
     DRN -Q, either route: residual                      0                    0
  C  build, SFR observation route +Q: phi                327,106.2            REFUSED
     build, -Q on each route: phi                        SFR budget route     0.0 on each
                                                         327,106.2, the
                                                         other three 0.0
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile
import warnings

import numpy as np
import pandas as pd
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]
LAST = len(REACHES) - 1
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
OBSERVED = {f"r{reach}": reach for reach in range(len(REACHES))}
GAUGED = {f"g{reach}": reach for reach in range(len(REACHES))}


def toy(root: pathlib.Path, *, mover: bool):
    """The toy in the module docstring, run; returns the live run model."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, LAST) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < LAST else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud", mover=mover,
        observations=mf.SfrFlowTargets(locations=OBSERVED).to_flopy_obs(
            csv_name="sfr_flow_targets.csv"),
    )
    drain = mf.drn(
        stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
        observations=mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).to_flopy_obs(
            None, csv_name="drn_flow_targets.csv"),
    )
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True, mover=mover,
        stage_filerecord="m.lak.stg",
    )
    packages = [
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream, drain, lake,
    ]
    if mover:
        packages.append(mf.mvr(nper=1, moves=(
            mf.Move(mf.MoverConnection("sfr", 1), mf.MoverConnection("lak", 0), "FACTOR", 0.3),
            mf.Move(mf.MoverConnection("sfr", LAST), mf.MoverConnection("lak", 0), "UPTO", 20.0),
        )))
    packages.append(mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
                          saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]))
    simulation = mf.simulation(mf.gwf("m", grid=grid, packages=packages), name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def show(values) -> str:
    return ", ".join(f"{float(v):+.3f}" for v in values)


def attempt(make):
    """``make()``, or the first line of what refused it."""

    try:
        return make()
    except ValueError as error:
        return f"REFUSED: {str(error)[:160]}..."


def section_a(runs) -> None:
    print("A. every reach on both routes, at each period end")
    for label, model in runs.items():
        observed = mf.SfrFlowTargets(locations=OBSERVED).simulated_series(model)
        budget = mf.SfrFlowTargets(locations=GAUGED).simulated_series(model)
        worst = max(
            float(np.max(np.abs(observed[f"r{r}"].to_numpy() - budget[f"g{r}"].to_numpy())
                         / np.abs(observed[f"r{r}"].to_numpy())))
            for r in range(len(REACHES)))
        print(f"   {label}: reach 0 observation route  {show(observed['r0'])}")
        print(f"   {label}: reach 0 budget route       {show(budget['g0'])}")
        print(f"   {label}: max |obs - budget| / |obs| over 4 reaches: {worst:.3g}")
    reader = CellBudgetFile(pathlib.Path(runs["mover"].workspace) / "m.sfr.bud", precision="double")
    for text in ("FLOW-JA-FACE", "TO-MVR", "EXT-OUTFLOW"):
        data = reader.get_data(text=text, kstpkper=(NSTP - 1, 0))[0]
        rows = [(int(n) - 1, round(float(q), 3)) for n, q in zip(data["node"], data["q"], strict=True) if q < 0]
        print(f"   mover, period 0 end, {text:>12} (reach, q < 0): {rows}")
    reader.close()


def section_b(model) -> None:
    print("\nB. compare() given the model's own flow as a gauged Q (plain run)")
    q = -mf.SfrFlowTargets(locations={"r0": 0}).simulated_series(model)["r0"].to_numpy()
    d = -mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(model)["seep"].to_numpy()
    print(f"   Q (reach 0) = {show(q)};  Q (zone) = {show(d)}")
    for family, make, name, values in (
        ("SFR observation", mf.SfrFlowTargets, "r0", q), ("SFR budget", mf.SfrFlowTargets, "g0", q),
        ("DRN observation", mf.DrnFlowTargets, "seep", d), ("DRN budget", mf.DrnFlowTargets, "zone", d),
    ):
        where = 0 if make is mf.SfrFlowTargets else DRN_CELLS
        for sign in (+1.0, -1.0):
            targets = attempt(lambda: make(locations={name: where}, values=pd.DataFrame(
                {"time": range(NPER), name: sign * values})))
            if isinstance(targets, str):
                print(f"   {family:>15}, entered {'+Q' if sign > 0 else '-Q'}: {targets}")
                continue
            residual = targets.compare(model).sort_values("per")["residual"].to_numpy()
            print(f"   {family:>15}, entered {'+Q' if sign > 0 else '-Q'}: residual {show(residual)}"
                  f"  phi {float(np.sum(residual ** 2)):,.1f}")


def section_c(model, root: pathlib.Path) -> None:
    import pyemu

    print("\nC. what a PEST build hands PEST (plain run, noptmax 0)")
    q = -mf.SfrFlowTargets(locations={"r0": 0}).simulated_series(model)["r0"].to_numpy()
    d = -mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(model)["seep"].to_numpy()
    cases = {
        "SFR observation route, +Q": [("up", mf.SfrFlowTargets, "r0", 0, q)],
        "every route, -Q": [("sfr_obs", mf.SfrFlowTargets, "r0", 0, -q),
                            ("sfr_bud", mf.SfrFlowTargets, "g0", 0, -q),
                            ("drn_obs", mf.DrnFlowTargets, "seep", DRN_CELLS, -d),
                            ("drn_bud", mf.DrnFlowTargets, "zone", DRN_CELLS, -d)],
    }
    for label, sets in cases.items():
        targets = [(prefix, attempt(lambda make=make, name=name, where=where, values=values: make(
            locations={name: where}, values=pd.DataFrame({"time": range(NPER), name: values}))))
            for prefix, make, name, where, values in sets]
        refused = [t for _, t in targets if isinstance(t, str)]
        if refused:
            print(f"   {label}: {refused[0]}")
            continue
        cal = model.pest(label[:3].lower(), workspace=root / label.replace(" ", "_").replace(",", ""),
                         start_datetime="2024-01-01")
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
        for prefix, t in targets:
            cal.observe(t, prefix=prefix)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pst = cal.build("q.pst", noptmax=0)
        template = pathlib.Path(cal.template_workspace)
        obs = pst.observation_data
        weighted = obs.loc[obs["weight"] > 0]
        for prefix, t in targets:
            kind = "sfr_flow" if isinstance(t, mf.SfrFlowTargets) else "drn_flow"
            table = template / f"{prefix}_simulated_{kind}.csv"
            read = pyemu.pst_utils.InstructionFile(f"{table}.ins").read_output_file(str(table))["obsval"]
            rows = weighted.index.intersection(read.index)
            phi = float(np.sum((read.loc[rows].astype(float) - weighted.loc[rows, "obsval"]) ** 2))
            print(f"   {label}: {prefix}: {len(rows)} weighted rows, phi {phi:,.1f}")


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="flow_target_sign_"))
    try:
        runs = {"plain": toy(root / "plain", mover=False), "mover": toy(root / "mover", mover=True)}
        section_a(runs)
        section_b(runs["plain"])
        section_c(runs["plain"], root / "pest")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
