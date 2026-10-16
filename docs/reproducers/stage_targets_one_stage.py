"""A stage target added a lake's stage once per connected cell, and a reach's once per step.

Ledger 331. ``LakeStageTargets``' stage-file route read ``packages.lak.results.stage``
-- one row per CONNECTED CELL, each a copy of the lake's one stage -- and
``observations._shared._aggregate_named_series`` SUMMED them. Measured by the lead on
the user's run ``sw_mounding_flow_pitk_9232026_QpogK15_exc`` (an excavated pond with
612 vertical and 1,202 horizontal connections on 763 distinct cells): period 1's stage
526.003 ft came back as 401,339.97 ft, 763 times itself, in the series PEST's base run
and every forward run wrote. ``SfrStageTargets`` shared the helper, and its table is one
row per reach per saved TIME STEP, so a period saving three steps read three stages
added together.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/stage_targets_one_stage.py [RUN_WS]

``RUN_WS`` (optional) is a run workspace holding a LAK stage file -- the user's pond, or
a copy of it; nothing there is written. Without it only the toys below run.

THE SETUP. The toy runs ``tests/test_stage_targets_one_stage.py`` uses: a 6 x 6
irregular Voronoi grid of ~100 ft cells, one confined layer, CHD on the west and east
edges.

  lake      six 10-day periods, one step each. Lake 0 has THREE vertical connections
            and a 900 ft3/d withdrawal (dry in records 1-3, refilled by rain in 4);
            lake 1 has TWO and is held CONSTANT at 108 ft. LAK STAGE observations
            are written for the names ``drying``/``held``.
  stream    four reaches, three periods of THREE time steps each, ``SAVE HEAD ALL``,
            heads filling from 90 ft toward 102 so the stage moves within a period.
  steps     one lake of THREE vertical connections and a 100 ft3/d withdrawal,
            three periods of THREE time steps each, ``SAVE HEAD ALL``, no LAK
            observations -- the stage-file route's CLOCK, not its sum.

The "before" of each series is the old rule -- the rows of one (per, id) ADDED,
``sum(min_count=1)`` -- applied by hand to the same table, so this still shows it
after the fix. A name the observation CSV does not hold takes the stage-file route.

WHAT IT MEASURES
  1. Rows per (record, lake) in the stage table, and the stage-file route before and
     after, against the stage file and against MF6's own STAGE observation.
  2. The same for SFR: rows per (period, reach), and the route before and after,
     against the stage at each period's LAST step and the observation there.
  3. The LAK stage-file route's clock. It was keyed by SAVED RECORD (ledger 331's
     limit, open under 301): at three saves per period, targets keyed by period 0, 1,
     2 read records 0, 1, 2 -- all inside period 0 -- in ``compare()`` and in a PEST
     build, where ``HeadTargets`` on the same run reads each period's last step. Since
     ledger 348 it is keyed by period and reads records 2, 5, 8, the period ends;
     ``docs/reproducers/lak_stage_clock.py`` has the before and after.
  4. With ``RUN_WS``: the same numbers as 1 on that run.
"""

from __future__ import annotations

import atexit
import pathlib
import shutil
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

NPER = 6


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    disv = dict(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], botm=[0.0])
    return grid, disv, edge


def _run(root: pathlib.Path, simulation):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def _lake_model(root: pathlib.Path):
    grid, disv, edge = _grid()
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    drying, held = interior[:3], interior[-2:]
    connections = [[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                   for i, c in enumerate(drying)]
    connections += [[1, i, (0, c), "VERTICAL", 0.01, 0.0, 0.0, 0.0, 0.0]
                    for i, c in enumerate(held)]
    keep = [[1, "STATUS", "CONSTANT"], [1, "STAGE", 108.0]]
    lake = mf.lak.flopy(
        nlakes=2,
        packagedata=[[0, 100.5, 3, "drying"], [1, 108.0, 2, "held"]],
        connectiondata=connections,
        perioddata={0: [[0, "WITHDRAWAL", 900.0], *keep], 4: [[0, "RAINFALL", 0.5], *keep],
                    5: [[0, "RAINFALL", 0.0], *keep]},
        boundnames=True, stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
        observations=mf.LakeStageTargets(locations={"drying": 0, "held": 1}).to_flopy_obs(
            csv_name="lake_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(top=100.0, **disv), mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=50.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 50.0] for c in edge]}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 1, 1.0)] * NPER)]))


def _sfr_model(root: pathlib.Path):
    grid, disv, edge = _grid()
    row = [13, 14, 15, 16]
    last = len(row) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(row)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(row))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]], 2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
        observations=mf.SfrStageTargets(locations={"up": 0, "down": last}).to_flopy_obs(
            csv_name="sfr_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(top=110.0, **disv), mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}), stream,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(10.0, 3, 1.0)] * 3)]))


def _steps_model(root: pathlib.Path):
    grid, disv, edge = _grid()
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    connections = [[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                   for i, c in enumerate(interior[:3])]
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, 3, "pond"]], connectiondata=connections,
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]}, boundnames=True,
        stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(top=100.0, **disv), mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=50.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 50.0] for c in edge]}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(10.0, 3, 1.0)] * 3)])), interior[0]


def _show(values) -> str:
    return "[" + ", ".join("nan" if np.isnan(v) else f"{v:.10g}" for v in values) + "]"


def _old_rule(table: pd.DataFrame, id_column: str) -> pd.DataFrame:
    """The reduction before ledger 331: every row of one (per, id) added."""

    return table.groupby(["per", id_column])["stage"].sum(min_count=1).unstack(id_column)


def _lake(model, names: dict[str, int], *, label: str) -> None:
    table = model.packages.lak.results.stage.get()
    rows = table.groupby(["per", "lake"]).size().unstack("lake")
    stage = model.outputs.lak.stage.get()
    print(f"   {label}: rows per (record, lake): "
          + ", ".join(f"lake {lake} {sorted(set(rows[lake]))}" for lake in rows.columns))
    old = _old_rule(table, "lake")
    series = mf.LakeStageTargets(locations=names).simulated_series(model)
    for name, lake in names.items():
        print(f"   lake {lake} stage file   {_show(stage[:8, lake])}")
        print(f"   {'':>6} before       {_show(old[lake].to_numpy()[:8])}")
        print(f"   {'':>6} after        {_show(series[name].to_numpy()[:8])}")


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="stage_targets_one_stage_"))
    atexit.register(shutil.rmtree, root, ignore_errors=True)  # /tmp is a quota'd tmpfs (ledger 460)

    print("1. a lake: its stage once per connected cell, before; once, after")
    lake = _lake_model(root / "lake")
    _lake(lake, {"a": 0, "b": 1}, label="toy")
    observed = mf.LakeStageTargets(locations={"drying": 0, "held": 1}).simulated_series(lake)
    print(f"   MF6's STAGE observation, lake 0 {_show(observed['drying'])} "
          f"(the bottom, 100, where the file says dry)")
    print(f"   {'':>24} lake 1 {_show(observed['held'])}")
    print(f"   time labels: observation {list(observed['time'])}, "
          f"stage file {list(mf.LakeStageTargets(locations={'a': 0}).simulated_series(lake)['time'])}")

    print("\n2. a reach: every saved step of a period added, before; the last step, after")
    stream = _sfr_model(root / "sfr")
    table = stream.packages.sfr.results.stage.get()
    print(f"   rows per (period, reach): {sorted(set(table.groupby(['per', 'reach']).size()))}")
    table = table.assign(kstp=table["kstpkper"].map(lambda value: int(value[0])))
    old = _old_rule(table, "reach")
    series = mf.SfrStageTargets(locations={"a": 0, "b": 3}).simulated_series(stream)
    observed = mf.SfrStageTargets(locations={"up": 0, "down": 3}).simulated_series(stream)
    for name, obs_name, reach in (("a", "up", 0), ("b", "down", 3)):
        steps = table[table["reach"] == reach].pivot_table(index="per", columns="kstp",
                                                           values="stage")
        print(f"   reach {reach} steps 0/1/2 by period {steps.round(4).to_numpy().tolist()}")
        print(f"   {'':>7} before     {_show(old[reach])}")
        print(f"   {'':>7} after      {_show(series[name])}")
        # One row per period at its end since ledger 346 (every step, [2::3], before).
        print(f"   {'':>7} observation at each period end "
              f"{_show(observed[obs_name].to_numpy())}")
    print(f"   time labels: observation {observed['time'].round(3).tolist()}, "
          f"stage file {series['time'].tolist()}")

    print("\n3. a lake at three saves per period: the stage-file route's clock (ledger 348)")
    steps, cell = _steps_model(root / "steps")
    stage = steps.outputs.lak.stage.get()[:, 0]
    series = mf.LakeStageTargets(locations={"a": 0}).simulated_series(steps)
    print(f"   stage file: {len(stage)} records; route time {list(series['time'])}")
    print(f"   period ends (records 2, 5, 8) {_show(stage[2::3])}")
    targets = mf.LakeStageTargets(
        locations={"a": 0}, values=pd.DataFrame({"time": [0, 1, 2], "a": [111.0, 222.0, 333.0]}))
    compare = targets.compare(steps).sort_values("per")
    print(f"   compare() keyed by period 0/1/2 {_show(compare['sim_stage'])} "
          f"(records 2, 5, 8 since ledger 348; records 0, 1, 2 -- steps 1-3 of period 0 "
          f"-- before)")
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [cell], "layer": [0]}),
        values=pd.DataFrame({"time": [0, 1, 2], "name": "h", "head": 50.0}))
    table = steps.all_heads.reset_index()
    table = table[(table["layer"] == 0) & (table["cell"] == cell)]
    head_at = {(int(kper), int(kstp)): float(elev)
               for (kstp, kper), elev in zip(table["kstpkper"], table["elev"], strict=True)}
    print(f"   HeadTargets beside it, cell {cell}: {_show(heads.simulated_heads(steps)['h'])}; "
          f"head at each period's step 3 {_show([head_at[(p, 2)] for p in range(3)])}")
    cal = steps.pest("lk", workspace=root / "steps_pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        observations = cal.build("lk.pst").observation_data
    lake = observations[observations.obsnme.str.contains("usecol:a_")]
    print("   PEST build (obsnme -> obsval, weight):")
    for name, value, weight in zip(lake.obsnme, lake.obsval, lake.weight, strict=True):
        print(f"     {name.split('usecol:')[1]:<10} {value:12.6f}  {weight:g}")

    if len(sys.argv) > 1:
        print(f"\n4. {sys.argv[1]}")
        model = mf.load_mf6_run(sys.argv[1], verbosity_level=0)
        nlakes = int(model.fl.lak.nlakes.data)
        _lake(model, {f"lake{lake}": lake for lake in range(nlakes)}, label="run")


if __name__ == "__main__":
    main()
