"""What PEST is handed matches the targets: a PEST build, read back through PEST++ itself.

Ledger 340. The build rebuilt each target's observation NAME from the spelling of
its key and skipped any name the control file did not hold; it left every simulated
row no target named at pyEMU's defaults (weight 1, the BASE run's own value); it wrote
the base head table through a pivot, whose columns pandas sorts, while every forward
run wrote mapping order; and its lake/SFR/DRN forward run reloaded the simulation's
FIRST GWF model. This script measures each through what PEST++ actually reads.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_handed_targets.py

It uses only the public API, so it measures whichever tree is on ``PYTHONPATH`` --
the "before" column below is this script run against the parent commit's ``src``
(``git archive b46366b src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``).
Needs ``mf6`` and ``pestpp-glm``.

THE SETUP. ``tests/test_pest_handed_targets.py``'s model: a 6 x 6 irregular Voronoi
grid of ~100 ft cells, one confined layer, CHD at 102 ft west and east, heads filling
from 90 ft, a four-reach stream, a three-cell DRN zone and a two-cell lake; three
10-day periods of three steps, ``SAVE HEAD LAST`` + ``SAVE BUDGET LAST``, no MF6 OBS,
so every family reads its FILE route, keyed by stress period. One parameter, ``k``.

  A. Observed: heads W2, W1, W10 (cells 1, 2, 34 -- names NOT in sorted order) keyed
     ``0.0, 1.0`` (floats); SFR flow ``a`` (reach 0) keyed ``1.0``. Forecast: head P
     (cell 23) at period ``2.0``. ``pestpp-glm`` NOPTMAX=0 runs the forward run once;
     its ``.rei`` says what PEST read for each observation, which is compared with the
     column the observation is NAMED after in the forward run's own CSV.
  B. One refusal or crash per defect: a target at period 5 of 3; two readings of one
     well in one period; a wide table with a blank; a value that is not a number; a
     forecast given no values.
  C. Two GWF models, calibrating the SECOND against a DRN zone: the base table and a
     forward run's.

MEASURED (mf6 6.8.0.dev0, pestpp-glm, pyEMU 1.4.0, 2026-09-27)   before          after
  A  target rows carrying their value in the pst (of 7)             0              7
     observations at weight > 0 (nnz_obs, of nobs 15)              15              7
     forecast rows listed / zero-weighted (P has 3 rows)          0 / 0          1 / 3
     head observations PEST read from ANOTHER well's column      9 of 9         0 of 9
  B  target at period 5 of 3                          built (nnz_obs 3)      REFUSED
     two readings (99.0, 95.0) of A in period 0       built, obsval 95.0      REFUSED
     a blank in a wide table                  CRASH, pyEMU: "NaNs in       built, weight 0
                                              obs_data" (and it writes
                                              pst.obs_data.nans.csv into
                                              the CALLER's cwd)
     "12,5" as a head                           read as a blank (NaN)       REFUSED
     forecast given no values                   CRASH, KeyError: 'time'    built, 3 rows listed
  C  forward-run DRN flow, calibrating model n     m's: -16.619, -17.705   n's: -804.167,
                                                   (base: n's -804.167,       -853.364, equal
                                                   -853.364)                  to the base

"9 of 9" is every head: the base header was the sorted pivot ``per,W1,W10,W2`` and the
forward run wrote ``per,W2,W1,W10``, so PEST's instruction file, which reads by
position and names what it reads after the base header, read W2's head as W1's, W1's
as W10's and W10's as W2's in every realization -- ``pestpp-glm`` exits 0 either way.
Why each is so, and what was decided: ledger 340.
"""

from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
import pyemu

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
PESTPP = shutil.which("pestpp-glm") or str(pathlib.Path.home() / ".local/bin/pestpp-glm")


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    disv = mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                   vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0])
    return grid, disv, edge


def _run(root, simulation):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run


def _toy(root):
    grid, disv, edge = _grid()
    last = len(REACHES) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud")
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 10.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg")
    model = mf.gwf("m", grid=grid, packages=[
        disv, mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}), stream,
        mf.drn(stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    return _run(root, simulation).model()


def _heads(names_cells, values):
    return mf.HeadTargets(
        locations=pd.DataFrame({"name": [n for n, _ in names_cells],
                                "cell": [c for _, c in names_cells], "layer": 0}),
        values=values)


def _build(model, workspace, observe=(), forecast=()):
    cal = model.pest("handed", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in observe:
        cal.observe(targets)
    for targets in forecast:
        cal.forecast(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("handed.pst", noptmax=0)
    return cal, pst


def _observation_at(template, table, period, column):
    """The observation PEST reads from one cell of a table: write a probe, read it back."""

    frame = pd.read_csv(template / table)
    probe = frame.astype({column: float})
    probe.loc[frame.iloc[:, 0].astype(float) == period, column] = 123456.789
    probe.to_csv(template / "probe.csv", index=False)
    read = pyemu.pst_utils.InstructionFile(str(template / f"{table}.ins")).read_output_file(
        str(template / "probe.csv"))["obsval"]
    return read.index[np.isclose(read.to_numpy(dtype=float), 123456.789)][0]


def section_a(model, root):
    wells = _heads([("W2", 1), ("W1", 2), ("W10", 34)], pd.DataFrame(
        {"time": [0.0, 1.0], "W2": [101.0, 101.5], "W1": [100.0, 100.5],
         "W10": [99.0, 99.5]}))
    flow = mf.SfrFlowTargets(locations={"a": 0}, values=pd.DataFrame({"time": [1.0], "a": [7.0]}))
    predicted = _heads([("P", 23)], pd.DataFrame({"time": [2.0], "P": [0.0]}))
    cal, pst = _build(model, root / "a", observe=(wells, flow), forecast=(predicted,))
    template = pathlib.Path(cal.template_workspace)
    obs = pst.observation_data

    wanted = [("hds_simulated_heads.csv", per, well, value)
              for well, values in (("W2", (101.0, 101.5)), ("W1", (100.0, 100.5)),
                                   ("W10", (99.0, 99.5)))
              for per, value in enumerate(values)]
    wanted.append(("sfr_flow_simulated_sfr_flow.csv", 1, "a", 7.0))
    carried = sum(obs.loc[_observation_at(template, t, p, c), "obsval"] == v
                  for t, p, c, v in wanted)
    forecast_rows = [_observation_at(template, "fore1_simulated_heads.csv", p, "P")
                     for p in range(NPER)]
    listed = [name for name in pst.pestpp_options.get("forecasts", "").split(",") if name]
    print(f"A  target rows carrying their value: {carried} of {len(wanted)}")
    print(f"   nnz_obs {pst.nnz_obs} of nobs {pst.nobs}")
    print(f"   forecast rows listed {len(listed)}, zero-weighted "
          f"{int((obs.loc[forecast_rows, 'weight'] == 0).sum())} of {len(forecast_rows)}")

    done = subprocess.run([PESTPP, "handed.pst"], cwd=template, capture_output=True,
                          text=True, timeout=600)
    rei = pyemu.pst_utils.read_resfile(str(next(template.glob("handed*.rei"))))
    written = pd.read_csv(template / "hds_simulated_heads.csv")
    heads = rei[rei["name"].str.startswith("oname:hds_")]
    wrong = []
    for name, modelled in zip(heads["name"], heads["modelled"], strict=True):
        label, per = re.search(r"usecol:(.+)_per:(\d+)$", name).groups()
        column = next(c for c in written.columns if c.lower() == label)
        if not np.isclose(modelled, written.loc[int(per), column], rtol=1e-9):
            wrong.append(name)
    print(f"   pestpp-glm exit {done.returncode}; forward-run header "
          f"{','.join(written.columns)}; PEST read another well's column for "
          f"{len(wrong)} of {len(heads)} head observations")


def _attempt(label, action):
    try:
        outcome = action()
    except Exception as error:  # noqa: BLE001 - reports whatever each tree does, crash included
        text = str(error).splitlines()[0]
        print(f"B  {label}: {type(error).__name__}: {text[:150]}")
        return
    print(f"B  {label}: {outcome}")


def section_b(model, root):
    late = _heads([("W", 22)], pd.DataFrame({"time": [5], "W": [50.0]}))
    _attempt("target at period 5 of 3",
             lambda: f"built, nnz_obs {_build(model, root / 'b1', observe=(late,))[1].nnz_obs}")

    def twice():
        heads = mf.HeadTargets(
            locations=pd.DataFrame({"name": ["A"], "cell": [22], "layer": [0]}),
            values=pd.DataFrame({"time": [0, 0], "name": ["A", "A"], "head": [99.0, 95.0]}))
        cal, pst = _build(model, root / "b2", observe=(heads,))
        name = _observation_at(pathlib.Path(cal.template_workspace), "hds_simulated_heads.csv", 0, "A")
        return f"built, obsval {pst.observation_data.loc[name, 'obsval']}"

    _attempt("two readings (99.0, 95.0) of A in period 0", twice)

    def blank():
        gappy = _heads([("W1", 22), ("W2", 23)], pd.DataFrame(
            {"time": [0, 1], "W1": [50.0, None], "W2": [52.0, 53.0]}))
        cal, pst = _build(model, root / "b3", observe=(gappy,))
        name = _observation_at(pathlib.Path(cal.template_workspace), "hds_simulated_heads.csv", 1, "W1")
        return f"built, the blank row's weight {pst.observation_data.loc[name, 'weight']}"

    _attempt("a blank in a wide table", blank)

    def typo():
        targets = _heads([("W1", 22)], pd.DataFrame({"time": [0, 1], "W1": ["50", "12,5"]}))
        return f"constructed, the value read as {targets.get()['head_target'].tolist()}"

    _attempt('"12,5" as a head', typo)

    def valueless():
        cal, pst = _build(model, root / "b5",
                          forecast=(mf.SfrFlowTargets(locations={"out": len(REACHES) - 1}),))
        listed = [n for n in pst.pestpp_options.get("forecasts", "").split(",") if n]
        return f"built, {len(listed)} forecast rows listed"

    _attempt("forecast given no values", valueless)


def section_c(root):
    grid, disv, edge = _grid()

    def flow_model(name, cond):
        return mf.gwf(name, grid=grid, packages=[
            disv, mf.npf(k=10.0, icelltype=0),
            mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=100.0),
            mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
            mf.drn(stress_period_data={0: [[(0, c), 96.0, cond] for c in DRN_CELLS]}),
            mf.oc(head_filerecord=f"{name}.hds", budget_filerecord=f"{name}.cbc",
                  saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])])

    run = _run(root / "c", mf.simulation(
        flow_model("m", 1.0), flow_model("n", 50.0), name="s",
        packages=[mf.tdis(nper=2, perioddata=[(PERLEN, 1, 1.0)] * 2), mf.ims()]))
    targets = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                                values=pd.DataFrame({"time": [0, 1], "seep": [-500.0, -500.0]}))
    cal = run.model("n").pest("two", workspace=root / "c" / "pest", start_datetime="2024-01-01")
    cal.parameterize("drn.cond", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 1000.0))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("two.pst", noptmax=0)
    table = pathlib.Path(cal.template_workspace) / "drn_flow_simulated_drn_flow.csv"
    base = pd.read_csv(table)["seep"].tolist()
    subprocess.run([sys.executable, "forward_run.py"], cwd=table.parent, check=True,
                   capture_output=True, text=True)
    forward = pd.read_csv(table)["seep"].tolist()
    first = targets.simulated_series(run.model("m"))["seep"].tolist()
    print(f"C  base (model n) {[round(v, 3) for v in base]}; forward run "
          f"{[round(v, 3) for v in forward]}; model m reads {[round(v, 3) for v in first]}")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as scratch:
        root = pathlib.Path(scratch)
        toy = _toy(root / "toy")
        section_a(toy, root)
        section_b(toy, root)
        section_c(root)
