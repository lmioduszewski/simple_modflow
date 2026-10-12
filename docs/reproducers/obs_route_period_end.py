"""The observation route is one row per stress period, at its end (ledger 346).

MF6 writes a continuous observation at the end of EVERY time step, labelled by the
simulated time ``totim`` and nothing else. ``LakeStageTargets``, ``SfrStageTargets``,
``SfrFlowTargets`` and ``DrnFlowTargets`` handed every row on, keyed by that time,
where every other route is one row per stress period at its end, keyed by the
period (ledger 336). This script measures what that did and what reads now.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/obs_route_period_end.py

It uses the public API (and ``period_end_rows`` where the tree has it), so it measures
whichever tree is on ``PYTHONPATH`` -- the "before" column below is this script run
against the parent commit's ``src`` (``git archive 0c38bf5 src | tar -x -C <dir>``, then
``PYTHONPATH=<dir>/src``). Needs ``mf6`` and ``pestpp-glm``.

THE SETUP. ``tests/test_obs_route_period_end.py``'s toy: a 6 x 6 irregular Voronoi grid
of ~100 ft cells, one confined layer, CHD at 102 ft west and east, heads filling from
90 ft; a four-reach stream, a three-cell DRN zone and a two-cell lake holding water
over its 110 ft bed (115 ft, drained by leakage), so every observed quantity moves
within a period. OBS declares LAK ``STAGE`` ``pond``, SFR ``STAGE`` ``gage`` (reach 1),
SFR ``DOWNSTREAM-FLOW`` ``up`` (reach 0) and DRN ``seep``, so those names take the
observation route. Three 10-day periods, ``SAVE ALL``.

  A. Three steps per period: each family's observation route -- its rows, and
     ``compare()`` keyed by period (0, 1, 2) and by MF6's time (10, 20, 30).
  B. A PEST build of all four keyed by period.
  C. ATS, one declared step per period: the rows MF6 writes at K = 10, 1 and 2, and
     what PEST's instruction file, built on the K = 10 base run, reads from a K = 2
     realization's forward-run CSV (a capture template, so both trees build it).
  D. Six ONE-HOUR periods, read live: MF6's period-end time against the cumulative
     length the model holds, the cumulative rule, and what the route reads.
  E. ``pestpp-glm`` NOPTMAX=0, end to end: a K = 20 "truth" run supplies each family's
     value and a head at each period's end; the K = 10 base is calibrated against them
     keyed by period; phi from PEST++'s own ``.rei`` at the base (multiplier 1) and at
     the truth (multiplier 2).
  F. How pandas reads MF6's observation text: the 27 DRN cell values of A against the
     cell-budget file.
  G. ``canonical_feature_signals``' seepage peak on the canonical model (testing
     profile, two steps per period, as the notebooks run it).

MEASURED (mf6 6.8.0.dev0, pestpp-glm, pyEMU 1.4.0, 2026-09-29)     before           after
  A  rows per family (9 steps, 3 periods)                        9                3
     compare() keyed by period 0, 1, 2 (every family)        NaN, NaN, NaN    each period's end
                                                                              (DRN -710.391,
                                                                              -855.018, -865.013)
     compare() keyed by time 10, 20, 30                       those period     REFUSED, naming
                                                              ends             "10 is the end of
                                                                              period 0 ... key it 0"
  B  build keyed by period, all four families                 REFUSED (340's   built: 12 rows,
                                                              interim)         12 weighted
  C  rows MF6 writes at K = 10 / 1 / 2                          12 / 13 / 14     12 / 13 / 14
     base CSV rows (K = 10) / realization CSV rows (K = 2)      12 / 14          3 / 3
     .ins reading the K = 2 realization's period ends        0.0, -274.426,   0.0, -381.261,
     (its true ends: 0.0, -381.261, -581.729)                -462.621         -581.729
  D  |MF6's end - cumulative held length|, periods 0..5      3.3e-9 .. 2.0e-8 (the same files)
     the cumulative rule at END_RTOL 1e-8                    refuses 3, 4, 5  --
     route rows                                               18 (every step)  6, every period
  E  phi at base K = 10 / at truth K = 20                    refused at build  18,754.7416 /
                                                              (period keys)     0.0000
  F  DRN cell values pandas' default parser reads exactly    23 of 27         (round_trip: 27)
     the route's one-cell zones equal to the budget file     6 of 9           9 of 9
  G  seepage_peak north_springs / south_springs              1332.889 /       1327.946 /
                                                              894.203          891.231
                                                              (max, 12 steps)  (max, 6 period ends)

The "before" of B and E is ledger 340's interim: period keys on this route are refused
(time keys 10, 20, 30 build there). F1's scoping measured the tree before 340
(``b46366b``), on a toy whose lake was dry: the same period-keyed targets built with
every row at weight 1 and the base run's value, phi 0.3456 at the base K and 178,835.03
at the truth. Why each is so, and what was decided: ledger 346.
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
import pyemu
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.period_end import END_RTOL, period_lengths

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

# PEST++ starts each forward run in the template, where a relative PYTHONPATH would
# import another myflopy: made absolute here, before anything changes directory
# (ledger 425).
make_pythonpath_absolute()

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
HEAD_CELL = 22
PESTPP = shutil.which("pestpp-glm") or str(pathlib.Path.home() / ".local/bin/pestpp-glm")


def _toy(root: pathlib.Path, *, ats=False, k=10.0, perlen=PERLEN, nper=NPER):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    last = len(REACHES) - 1
    stream_obs = mf.SfrFlowTargets(locations={"up": 0}).to_flopy_obs(
        csv_name="sfr_flow_targets.csv")
    stream_obs.update(mf.SfrStageTargets(locations={"gage": 1}).to_flopy_obs(
        csv_name="sfr_stage_targets.csv"))
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud", observations=stream_obs)
    drain = mf.drn(
        stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
        observations=mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).to_flopy_obs(
            None, csv_name="drn_flow_targets.csv"))
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 115.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.001, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg",
        observations=mf.LakeStageTargets(locations={"pond": 0}).to_flopy_obs(
            csv_name="lake_stage_targets.csv"))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=k, icelltype=0), mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0), mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream, drain, lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    timing = {"ats": True, "ats_dt0": 1.0, "ats_dtmin": 0.01, "ats_dtmax": 4.0} if ats else {}
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(perlen, 1 if ats else NSTP, 1.0)] * nper, **timing)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _families():
    return {
        "lake stage": (mf.LakeStageTargets, {"pond": 0}, "pond", "sim_stage"),
        "sfr stage": (mf.SfrStageTargets, {"gage": 1}, "gage", "sim_stage"),
        "sfr flow": (mf.SfrFlowTargets, {"up": 0}, "up", "sim_flow"),
        "drn flow": (mf.DrnFlowTargets, {"seep": DRN_CELLS}, "seep", "sim_flow"),
    }


def _try(call) -> str:
    try:
        return call()
    except ValueError as error:
        return f"REFUSED -- {str(error)[:260]}"


def _round(values) -> list:
    return [round(float(v), 3) for v in values]


def _period_end_values(model, csv: str, columns: list[str], perlen=PERLEN, nper=NPER):
    raw = pd.read_csv(pathlib.Path(model.workspace) / csv, float_precision="round_trip")
    rows = [int(np.flatnonzero(np.isclose(raw["time"], perlen * (per + 1), rtol=1e-6))[0])
            for per in range(nper)]
    return raw.iloc[rows].loc[:, columns].sum(axis=1).to_numpy()


def _placeholder(make) -> float:
    """A target value that only has to be accepted: 1.0, or -1.0 for a FLOW, which is
    entered in MF6's own sign since ledger 350 (a value above zero is refused)."""

    return -1.0 if make in (mf.SfrFlowTargets, mf.DrnFlowTargets) else 1.0


def section_a(toy):
    print("A. three steps per period: the observation route")
    for family, (make, locations, name, simulated) in _families().items():
        series = make(locations=locations).simulated_series(toy)
        print(f"   {family}: {len(series)} rows, time {_round(series['time'])}")
        for label, keys in (("period", [0, 1, 2]), ("time", [10, 20, 30])):
            targets = make(locations=locations,
                           values=pd.DataFrame({"time": keys, name: _placeholder(make)}))
            got = _try(lambda targets=targets: str(_round(targets.compare(toy)[simulated])))
            print(f"      compare() keyed by {label} {keys}: {got}")


def _build(model, workspace, observe):
    cal = model.pest("obsper", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for prefix, targets in observe.items():
        cal.observe(targets, prefix=prefix)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("obsper.pst", noptmax=0)
    return cal, pst


def section_b(toy, root):
    print("\nB. a PEST build of all four, keyed by period")
    observe = {
        family.replace(" ", "_"): make(locations=locations, values=pd.DataFrame(
            {"time": [0, 1, 2], name: _placeholder(make)}))
        for family, (make, locations, name, _sim) in _families().items()
    }

    def build():
        cal, pst = _build(toy, root / "b", observe)
        obs = pst.observation_data
        rows = {prefix: len(pd.read_csv(path))
                for prefix in observe
                for path in pathlib.Path(cal.template_workspace).glob(f"{prefix}_simulated_*.csv")}
        return f"built: {len(obs)} observations, {int((obs.weight > 0).sum())} weighted; rows {rows}"

    print(f"   {_try(build)}")


def section_c(root):
    print("\nC. ATS, one declared step per period")
    runs = {k: _toy(root / f"ats_k{k:g}", ats=True, k=k) for k in (10.0, 1.0, 2.0)}
    for k, model in runs.items():
        raw = pd.read_csv(pathlib.Path(model.workspace) / "drn_flow_targets.csv")
        print(f"   K = {k:g}: MF6 wrote {len(raw)} rows, times {_round(raw['time'])}")
    template = mf.DrnFlowTargets(locations={"seep": DRN_CELLS})  # a capture template
    cal, _pst = _build(runs[10.0], root / "c", {"drn": template})
    workspace = pathlib.Path(cal.template_workspace)
    table = workspace / "drn_simulated_drn_flow.csv"
    base_rows = len(pd.read_csv(table))
    realization = root / "c" / "realization.csv"
    from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

    write_named_series_targets(
        [{"kind": "drn_flow", "locations_file": workspace / "drn_target_locations.csv",
          "output_csv": realization}],
        sim_ws=runs[2.0].workspace, model_name="m")
    read = pyemu.pst_utils.InstructionFile(str(workspace / f"{table.name}.ins")).read_output_file(
        str(realization))["obsval"]
    truth = _period_end_values(runs[2.0], "drn_flow_targets.csv",
                               [f"SEEP_C{c}" for c in DRN_CELLS])
    print(f"   base CSV rows {base_rows}; K = 2 realization CSV rows "
          f"{len(pd.read_csv(realization))}")
    # The observations the base run named for the three period ends: t = 10, 20, 30 on
    # a route keyed by MF6's time, periods 0, 1, 2 on one keyed by period.
    keys = {float(name.split("time:")[1]): name for name in read.index}
    ends = [keys[k] for k in ((10.0, 20.0, 30.0) if 30.0 in keys else (0.0, 1.0, 2.0))]
    print(f"   PEST reads from it, under the base run's period-end names: "
          f"{[(name.split('usecol:')[1], round(float(read[name]), 3)) for name in ends]}")
    print(f"   the K = 2 run's period ends: {_round(truth)}")


def section_d(root):
    print("\nD. six one-hour periods, read live")
    hourly = _toy(root / "hourly", perlen=1.0 / 24.0, nper=6)
    held = np.cumsum(period_lengths(hourly))
    raw = pd.read_csv(pathlib.Path(hourly.workspace) / "drn_flow_targets.csv",
                      float_precision="round_trip")
    ends = raw["time"].to_numpy()[NSTP - 1::NSTP]
    gap = np.abs(ends - held)
    tol = END_RTOL * np.maximum(held, 1.0)
    for per in range(6):
        print(f"   period {per}: MF6 ends at {float(ends[per])!r}, the model holds "
              f"{float(held[per])!r}: "
              f"{gap[per]:.2e} apart -- the cumulative rule "
              f"{'reads' if gap[per] <= tol[per] else 'REFUSES'} it")
    series = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(hourly)
    print(f"   the route: {len(series)} rows, time {_round(series['time'])}")


def section_e(root):
    print("\nE. pestpp-glm NOPTMAX=0: phi at the base and at the truth")
    truth = _toy(root / "truth", k=20.0)
    measured = {
        "lake_stage": ("pond", _period_end_values(truth, "lake_stage_targets.csv", ["POND"])),
        "sfr_stage": ("gage", _period_end_values(truth, "sfr_stage_targets.csv", ["GAGE"])),
        "sfr_flow": ("up", _period_end_values(truth, "sfr_flow_targets.csv", ["UP"])),
        "drn_flow": ("seep", _period_end_values(truth, "drn_flow_targets.csv",
                                                [f"SEEP_C{c}" for c in DRN_CELLS])),
    }
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [HEAD_CELL], "layer": [0]}),
        values=pd.DataFrame({"time": [0, 1, 2], "name": "h", "head": 0.0}))
    true_heads = heads.compare(truth).sort_values("per")["sim_head"].to_numpy()
    base = _toy(root / "base", k=10.0)
    observe = {
        prefix: make(locations=locations,
                     values=pd.DataFrame({"time": [0, 1, 2], name: measured[prefix][1]}))
        for prefix, (make, locations, name, _sim) in zip(
            measured, _families().values(), strict=True)
    }
    observe["hds"] = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [HEAD_CELL], "layer": [0]}),
        values=pd.DataFrame({"time": [0, 1, 2], "name": "h", "head": true_heads}))

    def run():
        cal, _pst = _build(base, root / "e", observe)
        template = pathlib.Path(cal.template_workspace)
        phis = []
        for multiplier in (1.0, 2.0):
            work = root / f"e_glm_{multiplier:g}"
            shutil.copytree(template, work)
            pst = pyemu.Pst(str(work / "obsper.pst"))
            pst.parameter_data["parval1"] = multiplier
            pst.write(str(work / "obsper.pst"), version=2)
            done = subprocess.run([PESTPP, "obsper.pst"], cwd=work, capture_output=True,
                                  text=True, timeout=600)
            rei = pyemu.pst_utils.read_resfile(str(work / "obsper.rei"))
            weighted = rei.loc[rei["weight"] > 0]
            phis.append(float(((weighted["weight"] * weighted["residual"]) ** 2).sum()))
            assert done.returncode == 0, done.stdout[-2000:]
        return (f"phi at base K = 10: {phis[0]:.4f}; at truth K = 20: {phis[1]:.4f} "
                f"({len(weighted)} weighted observations)")

    print(f"   {_try(run)}")


def section_f(toy):
    print("\nF. how MF6's observation text is read")
    path = pathlib.Path(toy.workspace) / "drn_flow_targets.csv"
    cbc = CellBudgetFile(pathlib.Path(toy.workspace) / "m.cbc", precision="double")
    for label, frame in (("pandas default", pd.read_csv(path)),
                         ("round_trip", pd.read_csv(path, float_precision="round_trip"))):
        equal = 0
        for row, (kstp, kper) in enumerate(cbc.get_kstpkper()):
            drains = cbc.get_data(text="DRN", kstpkper=(kstp, kper))[0]
            equal += sum(int(drains["q"][drains["node"] == c + 1][0] == frame.loc[row, f"SEEP_C{c}"])
                         for c in DRN_CELLS)
        print(f"   {label}: {equal} of {len(frame) * len(DRN_CELLS)} cell values equal the budget file")
    cbc.close()
    cbc = CellBudgetFile(pathlib.Path(toy.workspace) / "m.cbc", precision="double")
    equal = 0
    for cell in DRN_CELLS:
        values = mf.DrnFlowTargets(locations={"seep": [cell]}).simulated_series(toy)["seep"]
        values = values.to_numpy()
        if len(values) == NPER * NSTP:  # every step: the period ends are every third
            values = values[NSTP - 1::NSTP]
        for per in range(NPER):
            drains = cbc.get_data(text="DRN", kstpkper=(NSTP - 1, per))[0]
            equal += int(drains["q"][drains["node"] == cell + 1][0] == values[per])
    cbc.close()
    print(f"   the route's one-cell zones at the period ends: {equal} of "
          f"{NPER * len(DRN_CELLS)} equal the budget file")


def section_g(root):
    print("\nG. canonical_feature_signals' seepage peak")
    from myflopy.modflow.mf6.canonical import canonical_feature_signals
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_transient_model

    config = CanonicalModelConfig(nrow=21, ncol=21, nlay=4, nper=6, cell_size=100.0,
                                  steps_per_period=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = build_transient_model(root / "canonical", config=config)
        success, _ = model.run_simulation()
        assert success
        series = model.targets.drn_flow.targets.simulated_series(model)
        signals = canonical_feature_signals(model)
    print(f"   DRN series rows {len(series)}; seepage_peak "
          f"{ {name: round(value, 3) for name, value in signals['seepage_peak'].items()} }")


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="obs_route_period_end_"))
    try:
        toy = _toy(root / "toy")
        section_a(toy)
        section_b(toy, root)
        section_c(root)
        section_d(root)
        section_e(root)
        section_f(toy)
        section_g(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
