"""A period's flow is its END-of-period flow: the flow targets added every saved step.

Ledger 336 (the user's rule: "a period's flow should be end of period value").
``SfrFlowTargets`` and ``DrnFlowTargets`` read their budget-file route by ADDING every
time step a stress period saved; ledger 331 measured 3.0x at ``nstp=3``. They read
the record that ENDS each period now, through the one rule every target family uses
(``myflopy.modflow.mf6.period_end``), and a period whose end the file does not hold
is refused.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/flow_targets_period_end.py

THE SETUP. The model ``tests/test_flow_targets_period_end.py`` runs: a 6 x 6 irregular
Voronoi grid of ~100 ft cells, one confined layer, CHD at 102 ft on the west and east
edges, heads filling from 90 ft so every flow moves within a period; a four-reach
stream (reach 3 is the outlet, with no downstream reach), a three-cell DRN zone and a
two-cell lake; three 10-day periods of THREE steps each. SFR ``DOWNSTREAM-FLOW`` and
per-cell DRN observations are written for ``up``/``down``/``seep``; any other name
reads the budget file. Then the canonical model (``CanonicalModelConfig`` testing
profile, 21 x 21 x 4, six 30-day periods) at THREE steps per period, ``SAVE ALL``.

The "before" of each flow is the old rule -- a period's saved steps ADDED -- applied
by hand to the budget file read with FloPy, so this still shows it after the fix.

WHAT IT MEASURES
  1. ``SAVE ALL``: each flow per saved step, the old rule, the new route, and MF6's
     own observation at each period end. The old SFR route could not run here at
     all: it asked the budget for ``TO-MVR``, which a stream with no mover never
     writes, and FloPy raises.
  2. ``SAVE HEAD ALL`` + ``SAVE BUDGET LAST``: 9 head records, 3 budget records; the
     flows read the budget file's own clock and match 1.
  3. ``SAVE BUDGET FIRST``: what the old rule read (step 1 of 3, silently) and the
     refusal now; heads beside it (``SAVE HEAD LAST``) still read.
  4. ``SAVE HEAD FIRST``: the same decision for heads, SFR stage and LAK stage.
  5. ATS, one declared step per period: the steps MF6 took, and the record read --
     found by ``pertim == perlen``, which ``kstp == nstp - 1`` would have refused.
  6. The observation route, recorded not changed: its SFR flow is MF6's NEGATIVE
     ``DOWNSTREAM-FLOW``; a DRN target keyed by PERIOD matches nothing in
     ``compare()``, and a PEST build keeps the simulated value as its ``obsval``.
  7. The canonical model: the outlet reach and a one-cell DRN zone, before and
     after, and MF6's ``DOWNSTREAM-FLOW`` at each period end.
  8. What the branch's review found, on the toy: (a) ONE-HOUR periods read live --
     the length the model holds against the one MF6 ended each period at, which
     FloPy wrote at 8 decimals, and the tolerance that refused it; (b) a DRN zone
     given as TWO ROWS on the observation route, which summed its last row only;
     (c) ``model.hds.get()`` under ``SAVE HEAD FIRST``, which read step 1 of 3 as
     each period's head; (d) why the group comparisons were LEFT: their ``per=``
     is the saved-record index (``_coerce_kstpkper``, ledger 301), so a period-end
     filter under ``_reduce_to_period_end`` alone, patched in here, breaks them.
"""

from __future__ import annotations

import pathlib
import tempfile
import warnings

import numpy as np
import pandas as pd
from flopy.utils import CellBudgetFile, HeadFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_transient_model
from myflopy.modflow.mf6.period_end import at_period_end, period_ends, period_lengths, saved_steps
from myflopy.project.group import spatial as group_spatial
from myflopy.project.group.core import ModelGroup

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]


def _toy(root: pathlib.Path, saverecord, *, ats: bool = False, perlen: float = PERLEN):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    last = len(REACHES) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]], 2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
        observations=mf.SfrFlowTargets(locations={"up": 0, "down": last}).to_flopy_obs(
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
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg",
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0), mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream, drain, lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc", saverecord=saverecord),
    ])
    timing = {"ats": True, "ats_dt0": 1.0, "ats_dtmin": 0.01, "ats_dtmax": 4.0} if ats else {}
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(perlen, 1 if ats else NSTP, 1.0)] * NPER, **timing)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def _flows(workspace, *, sfr_file: str, cbc_file: str, reaches, drn_cells):
    """``{(kstp, kper): {reach: flow, "zone": q}}`` from the budget files, zero-based."""

    sfr = CellBudgetFile(pathlib.Path(workspace) / sfr_file, precision="double")
    cbc = CellBudgetFile(pathlib.Path(workspace) / cbc_file, precision="double")
    texts = {name.decode().strip() for name in sfr.get_unique_record_names()}
    flows = {}
    for kstp, kper in sfr.get_kstpkper():
        key = (int(kstp), int(kper))
        routed = sfr.get_data(text="FLOW-JA-FACE", kstpkper=key)[0]
        row = {}
        for reach in reaches:
            here = routed[routed["node"] == reach + 1]
            total = -here["q"][here["q"] < 0].sum()
            old = total
            for text in ("TO-MVR", "EXT-OUTFLOW"):
                if text in texts:
                    term = sfr.get_data(text=text, kstpkper=key)[0]
                    part = -term["q"][term["node"] == reach + 1].sum()
                    total += part
                    old += part if text == "TO-MVR" else 0.0
            row[reach] = (float(total), float(old))
        drains = cbc.get_data(text="DRN", kstpkper=key)[0]
        row["zone"] = float(drains["q"][np.isin(drains["node"], np.asarray(drn_cells) + 1)].sum())
        flows[key] = row
    sfr.close()
    cbc.close()
    return flows


def _by_period(flows, pick):
    """``pick(row)`` per saved step, grouped by period in step order."""

    out = {}
    for (_kstp, kper), row in sorted(flows.items(), key=lambda item: (item[0][1], item[0][0])):
        out.setdefault(kper, []).append(pick(row))
    return out


def _show(values) -> str:
    return "[" + ", ".join(f"{v:.2f}" for v in values) + "]"


def _try(call) -> str:
    try:
        result = call()
    except ValueError as error:
        return f"REFUSED -- {error}"
    return _show(result)


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="flow_targets_period_end_"))
    last = len(REACHES) - 1

    print("1. SAVE ALL, nstp=3: a period's steps added, before; the step that ends it, after")
    toy = _toy(root / "all", [("HEAD", "ALL"), ("BUDGET", "ALL")])
    flows = _flows(toy.workspace, sfr_file="m.sfr.bud", cbc_file="m.cbc",
                   reaches=(0, last), drn_cells=DRN_CELLS)
    try:
        toy.outputs.sfr.bud.get("TO-MVR")
    except ValueError as error:
        print(f"   the old SFR route asked for TO-MVR on a stream with no mover: {error}")
    sfr = mf.SfrFlowTargets(locations={"a": 0, "b": last}).simulated_series(toy)
    drn = mf.DrnFlowTargets(locations={"zone": DRN_CELLS}).simulated_series(toy)
    observed_sfr = mf.SfrFlowTargets(locations={"up": 0, "down": last}).simulated_series(toy)
    observed_drn = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(toy)
    for label, pick, series, observed in (
        ("reach 0", lambda row: row[0][0], sfr["a"], observed_sfr["up"]),
        (f"reach {last} (outlet)", lambda row: row[last][0], sfr["b"], observed_sfr["down"]),
        ("DRN zone (3 cells)", lambda row: row["zone"], drn["zone"], observed_drn["seep"]),
    ):
        steps = _by_period(flows, pick)
        print(f"   {label}: steps by period {[_show(v) for v in steps.values()]}")
        if label.startswith("reach 0"):
            old = _by_period(flows, lambda row: row[0][1])
            print(f"      before (the old rule, had it run: steps added) "
                  f"{_show([sum(v) for v in old.values()])}")
        elif label.startswith("reach"):
            print("      before (the old rule, had it run): no row -- an outlet sends "
                  "nothing through FLOW-JA-FACE, only EXT-OUTFLOW")
        else:
            print(f"      before (sum of steps) {_show([sum(v) for v in steps.values()])}")
        print(f"      after                 {_show(series)}")
        print(f"      MF6's observation at each period end {_show(observed.iloc[NSTP - 1::NSTP])}")

    print("\n2. SAVE HEAD ALL + SAVE BUDGET LAST: each file's own clock")
    budget_last = _toy(root / "budget_last", [("HEAD", "ALL"), ("BUDGET", "LAST")])
    heads = HeadFile(pathlib.Path(budget_last.workspace) / "m.hds", precision="double")
    print(f"   head records {len(heads.get_kstpkper())}, budget records "
          f"{len(CellBudgetFile(pathlib.Path(budget_last.workspace) / 'm.cbc').get_kstpkper())}")
    heads.close()
    print(f"   DRN zone {_show(mf.DrnFlowTargets(locations={'zone': DRN_CELLS}).simulated_series(budget_last)['zone'])}"
          f" (1's after: {_show(drn['zone'])})")

    print("\n3. SAVE BUDGET FIRST: no period's end in the budget file")
    first = _toy(root / "budget_first", [("HEAD", "LAST"), ("BUDGET", "FIRST")])
    flows = _flows(first.workspace, sfr_file="m.sfr.bud", cbc_file="m.cbc",
                   reaches=(0, last), drn_cells=DRN_CELLS)
    print(f"   what the old rule read (step 1 of 3): DRN zone "
          f"{_show([row['zone'] for _, row in sorted(flows.items(), key=lambda i: i[0][1])])}")
    print(f"   now: {_try(lambda: mf.DrnFlowTargets(locations={'zone': DRN_CELLS}).simulated_series(first)['zone'])}")
    heads = mf.HeadTargets(locations=pd.DataFrame({"name": ["h"], "cell": [22], "layer": [0]}),
                           values=pd.DataFrame({"time": range(NPER), "name": "h", "head": 1.0}))
    print(f"   heads beside it (SAVE HEAD LAST): {_show(heads.simulated_heads(first)['h'].astype(float))}")

    print("\n4. SAVE HEAD FIRST: the same decision for heads and stages")
    head_first = _toy(root / "head_first", [("HEAD", "FIRST"), ("BUDGET", "LAST")])
    print(f"   heads:     {_try(lambda: heads.simulated_heads(head_first)['h'])[:160]}")
    print(f"   SFR stage: {_try(lambda: mf.SfrStageTargets(locations={'a': 0}).simulated_series(head_first)['a'])[:160]}")
    print(f"   LAK stage: {_try(lambda: mf.LakeStageTargets(locations={'p': 0}).simulated_series(head_first)['p'])[:160]}")
    print(f"   DRN zone beside them (SAVE BUDGET LAST): "
          f"{_show(mf.DrnFlowTargets(locations={'zone': DRN_CELLS}).simulated_series(head_first)['zone'])}")

    print("\n5. ATS, one declared step per period")
    ats = _toy(root / "ats", [("HEAD", "ALL"), ("BUDGET", "ALL")], ats=True)
    reader = CellBudgetFile(pathlib.Path(ats.workspace) / "m.cbc", precision="double")
    saved = saved_steps(reader)
    reader.close()
    for per in range(NPER):
        mine = [(kstp + 1, round(pertim, 6)) for kstp, kper, pertim in saved if kper == per]
        print(f"   period {per}: (step, pertim) saved {mine}")
    flows = _flows(ats.workspace, sfr_file="m.sfr.bud", cbc_file="m.cbc",
                   reaches=(0, last), drn_cells=DRN_CELLS)
    ends = {kper: row["zone"] for (kstp, kper), row in sorted(flows.items())}
    print(f"   DRN zone read {_show(mf.DrnFlowTargets(locations={'zone': DRN_CELLS}).simulated_series(ats)['zone'])}"
          f"; at each period's last step {_show(ends.values())}")

    print("\n6. the observation route (recorded, not changed)")
    print(f"   SFR 'up' observation at period ends {_show(observed_sfr['up'].iloc[NSTP - 1::NSTP])} "
          f"against the budget route {_show(sfr['a'])}")
    for label, times in (("by period", [0, 1, 2]), ("by time", [10, 20, 30])):
        compare = mf.DrnFlowTargets(
            locations={"seep": DRN_CELLS},
            values=pd.DataFrame({"time": times, "seep": -800.0})).compare(toy)
        print(f"   DRN target keyed {label} {times}: compare() sim_flow {list(compare['sim_flow'].round(2))}")
    targets = mf.DrnFlowTargets(locations={"seep": DRN_CELLS},
                                values=pd.DataFrame({"time": [0, 1, 2], "seep": -800.0}))
    cal = toy.pest("q", workspace=root / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        observations = cal.build("q.pst", noptmax=0).observation_data
    seep = observations[observations.obsnme.str.contains("usecol:seep")]
    print("   PEST build of that by-period target (obsnme -> obsval, weight):")
    for name, value, weight in zip(seep.obsnme, seep.obsval, seep.weight, strict=True):
        print(f"     {name.split('usecol:')[1]:<22} {value:12.4f}  {weight:g}")

    print("\n7. the canonical model, three steps per period, SAVE ALL")
    workspace = root / "canonical"
    config = CanonicalModelConfig(nrow=21, ncol=21, nlay=4, nper=6, cell_size=100.0,
                                  steps_per_period=3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = build_transient_model(workspace, config=config)
        success, _ = model.run_simulation()
    assert success
    run = mf.load_mf6_run(model.model_output_folder_path, verbosity_level=0)
    outlet, spring = 54, 317
    flows = _flows(run.workspace, sfr_file=f"{run.name}_sfr_budget.sfr",
                   cbc_file=f"{run.name}.cbc", reaches=(outlet,), drn_cells=[spring])
    old = _by_period(flows, lambda row: row[outlet][1])
    zone = _by_period(flows, lambda row: row["zone"])
    sfr = mf.SfrFlowTargets(locations={"outlet": outlet}).simulated_series(run)
    drn = mf.DrnFlowTargets(locations={"spring": [spring]}).simulated_series(run)
    observed = pd.read_csv(pathlib.Path(run.workspace) / "master_sfr_flow.csv")
    print(f"   outlet reach {outlet}: steps of period 0 {_show(_by_period(flows, lambda r: r[outlet][0])[0])}")
    print(f"      before {_show([sum(v) for v in old.values()])}")
    print(f"      after  {_show(sfr['outlet'])}")
    print(f"      MF6 DOWNSTREAM-FLOW at period ends {_show(observed['SFR_OUTFLOW'].iloc[2::3])}")
    print(f"   one-cell DRN zone, cell {spring}: steps of period 0 {_show(zone[0])}")
    print(f"      before {_show([sum(v) for v in zone.values()])}")
    print(f"      after  {_show(drn['spring'])}")

    print("\n8. what the review found")
    hour = 1.0 / 24.0
    hourly = _toy(root / "hourly", [("HEAD", "ALL"), ("BUDGET", "ALL")], perlen=hour)
    held = period_lengths(hourly)[0]
    reader = HeadFile(pathlib.Path(hourly.workspace) / "m.hds", precision="double")
    ended = saved_steps(reader)[NSTP - 1][2]
    reader.close()
    tdis = (pathlib.Path(hourly.workspace) / "s.tdis").read_text()
    written = tdis.split("BEGIN perioddata")[1].split()[0]
    print(f"   (a) one-hour periods, live: the model holds {held!r}; FloPy wrote {written}; "
          f"MF6 ended period 0 at {ended!r} -- {abs(ended - held):.2e} apart")
    print(f"       the old tolerance, 1e-9 x max(1, perlen) = {1e-9 * max(abs(held), 1.0):.0e}: "
          f"{'refused' if abs(ended - held) > 1e-9 * max(abs(held), 1.0) else 'read'} "
          f"every period; now: heads {_try(lambda: heads.simulated_heads(hourly)['h'].astype(float))}")
    rows = pd.DataFrame({"name": ["seep", "seep"], "cells": [[20], [21, 27]]})
    observed = pd.read_csv(pathlib.Path(toy.workspace) / "drn_flow_targets.csv")
    ends = observed.iloc[NSTP - 1::NSTP]
    print("   (b) a zone given as rows [20] and [21, 27], observation route at each period end:")
    print(f"       the last row's cells (what it summed) "
          f"{_show(ends['SEEP_C21'] + ends['SEEP_C27'])}; all three "
          f"{_show(ends['SEEP_C20'] + ends['SEEP_C21'] + ends['SEEP_C27'])}")
    series = mf.DrnFlowTargets(locations=rows).simulated_series(toy)
    print(f"       now {_show(series['seep'].iloc[NSTP - 1::NSTP])}; the budget route "
          f"{_show(mf.DrnFlowTargets(locations=rows.assign(name='zone')).simulated_series(toy)['zone'])}")
    first_step = HeadFile(pathlib.Path(head_first.workspace) / "m.hds", precision="double")
    step_one = [float(first_step.get_data(kstpkper=(0, per)).ravel()[22]) for per in range(NPER)]
    first_step.close()
    full = HeadFile(pathlib.Path(toy.workspace) / "m.hds", precision="double")
    period_end = [float(full.get_data(kstpkper=(NSTP - 1, per)).ravel()[22]) for per in range(NPER)]
    full.close()
    print(f"   (c) model.hds.get(cells=[22]) under SAVE HEAD FIRST read step 1 of 3 as each "
          f"period: {_show(step_one)}, where the periods end at {_show(period_end)} (SAVE ALL run)")
    print(f"       now: {_try(lambda: head_first.hds.get(cells=[22], layer=0)['head'])[:170]}")
    group = ModelGroup({"ref": str(toy.workspace), "var": str(toy.workspace)})
    shipped = group_spatial._reduce_to_period_end

    def period_end_alone(frame):
        pieces = []
        for name, rows in frame.groupby("model", sort=False):
            member = group.models[name]
            ends = period_ends(member, reader=lambda member=member: member.hds,
                               file="the head file", saves="HEAD", fallback=rows["kstpkper"])
            pieces.append(at_period_end(rows, ends))
        return shipped(pd.concat(pieces) if pieces else frame)

    print("   (d) group.hds.compare(per=, cells=[22]) on the SAVE ALL run grouped with itself:")
    for label, reduce in (("as shipped", shipped), ("period-end filter alone", period_end_alone)):
        group_spatial._reduce_to_period_end = reduce
        try:
            for per in range(NPER):
                try:
                    compare = group.hds.compare(per=per, cells=[22], layer=0)
                    got = f"kstpkper {list(compare['kstpkper'])} labelled per {list(compare['per'])}"
                except KeyError as error:
                    got = f"KeyError {error}"
                print(f"       {label}: per={per} -> {got}")
        finally:
            group_spatial._reduce_to_period_end = shipped


if __name__ == "__main__":
    main()
