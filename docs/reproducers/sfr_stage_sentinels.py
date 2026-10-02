"""Every SFR stage view drew MF6's -1e30 (dry) and +1e30 (inactive) as a stage.

Ledger 353, LAK's ledger 329 for SFR. The user's call: "the plots for SFR should just
not plot a stage if the reach is dry" -- and leave the SFR TARGETS as they are.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/sfr_stage_sentinels.py

THE SETUP. The same toy model ``tests/test_sfr_stage_sentinels.py`` runs: a 6 x 6
irregular Voronoi grid of ~100 ft cells, one confined layer (top 110, bottom 0, K 10),
CHD at 90 ft on the west and east edges, six 10-day transient periods of one step,
and five reaches whose beds (top 100 ft falling 0.1 ft a reach, 1 ft thick, K 0.001)
sit above the water table, so a reach with no inflow has no water:

  reaches 0-1-2  one stream, 50 ft3/d into reach 0; the inflow stops in periods 1-3,
                 so all three are DRY there (depth 0), and resumes in period 4.
  reach 3        alone, 30 ft3/d; STATUS INACTIVE in period 2.
  reach 4        alone, 100 ft3/d throughout: a control.

A second run ("wet") never stops the inflow and keeps reach 3 active, for the group
diff. SFR STAGE observations are declared under r0 .. r4. The "before" of each view is
its OLD rule applied by hand to the raw file (``outputs.sfr.stage.get(masked=False)``),
so this still shows it after the fix.

WHAT IT MEASURES
  1. What MF6 writes, in its stage file and in its STAGE observation.
  2. Each view, before and after: the series, a map's colour range, the long profile,
     the exchange hover, ``canonical_sfr_signals``, and the group diff summary.
  3. The targets, which do NOT change: the stage-file route still reads -1e30, and
     the observation route reads the streambed top, at a reach with no stage.
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import canonical_sfr_signals, irregular_voronoi_grid

NPER = 6
CELLS = (13, 14, 15, 20, 27)
BED_TOP = tuple(100.0 - 0.1 * reach for reach in range(5))
OBSERVED = {f"r{reach}": reach for reach in range(5)}


def _run(root: pathlib.Path, *, drying: bool):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    ncon = (1, 2, 1, 0, 0)
    stream = mf.sfr.flopy(
        packagedata=[[reach, (0, cell), 100.0, 5.0, 0.001, BED_TOP[reach], 1.0, 0.001,
                      0.03, ncon[reach], 1.0, 0] for reach, cell in enumerate(CELLS)],
        connectiondata=[[0, -1], [1, 0, -2], [2, 1], [3], [4]],
        perioddata={
            0: [[0, "INFLOW", 50.0], [3, "INFLOW", 30.0], [4, "INFLOW", 100.0]],
            1: [[0, "INFLOW", 0.0 if drying else 50.0]],
            2: [[3, "STATUS", "INACTIVE" if drying else "ACTIVE"]],
            3: [[3, "STATUS", "ACTIVE"]],
            4: [[0, "INFLOW", 50.0]],
        },
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
        observations=mf.SfrStageTargets(locations=OBSERVED).to_flopy_obs(
            csv_name="sfr_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0), mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        stream,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 1, 1.0)] * NPER)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def _show(values) -> str:
    return "[" + ", ".join("nan" if np.isnan(v) else f"{v:.4g}" for v in values) + "]"


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="sfr_stage_sentinels_"))
    try:
        _measure(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def _measure(root: pathlib.Path) -> None:
    dry = _run(root / "dry", drying=True)
    wet = _run(root / "wet", drying=False)
    raw = dry.outputs.sfr.stage.get(masked=False).to_numpy(dtype=float)  # (reach, record)

    print("1. what MF6 writes (rows = reaches 0..4, columns = periods 0..5)")
    for reach, row in enumerate(raw):
        print(f"   reach {reach}: {_show(row)}")
    obs = pd.read_csv(pathlib.Path(dry.workspace) / "sfr_stage_targets.csv")
    print(f"   STAGE observation, reach 0 (dry in 1-3): {_show(obs['R0'])}")
    print(f"   STAGE observation, reach 3 (inactive in 2): {_show(obs['R3'])}")

    print("\n2. each view, before (old rule by hand) and after")
    stage = dry.sfr.results.stage
    series = {t.name: np.asarray(t.y, float) for t in stage.plot().data}
    print(f"   series, reach 0   before {_show(raw[0])}")
    print(f"   {'':>17} after {_show(series['Reach 0'])}")
    for per in (1, 2):
        # Unconnected cells take the default fill_value, 0.0, before and after.
        values = np.asarray(stage.map(per=per).linear_zs, float)
        old = np.append(raw[:, per], 0.0)
        print(f"   map(per={per}) colour range before {old.min():.4g} .. {old.max():.4g}; "
              f"after {np.nanmin(values):.4g} .. {np.nanmax(values):.4g}")
    for per in (1, 2):
        trace = next(t for t in dry.sfr.results.profile.plot(per=per).data if t.name == "Stage")
        print(f"   long profile stage, per {per}: before {_show(raw[:, per])}, "
              f"after {_show(np.asarray(trace.y, float))}")
    trace = dry.sfr.results.q.map(per=2).fig.data[0]
    labels = [piece.split("%{customdata[")[0] for piece in trace.hovertemplate.split("<br>")]
    slot = next(i for i, label in enumerate(labels) if ">stage<" in label)
    shown = [trace.customdata[cell][slot] for cell in CELLS]
    print(f"   q.map(per=2) hover stage, reaches 0..4: before "
          f"{[f'{v:.4g} ft' for v in raw[:, 2]]}, after {shown}")
    signals = canonical_sfr_signals(dry, per=2)
    old_wet = np.isfinite(raw[:, 2]) & (raw[:, 2] > -1.0e20)  # the old rule on the raw file
    old_depth = (raw[:, 2] - np.array(BED_TOP))[old_wet]
    print(f"   canonical_sfr_signals(per=2) before: wet {int(old_wet.sum())}, "
          f"maximum_depth {old_depth.max():.4g}")
    print(f"   {'':>27} after: wet {signals['wet_reach_count']}, "
          f"maximum_depth {signals['maximum_depth']:.4g}")

    diff = mf.ModelGroup({"wet": wet, "dry": dry}, reference="wet").diff()
    raw_wet = wet.outputs.sfr.stage.get(masked=False).to_numpy(dtype=float)
    old = (raw - raw_wet).ravel()  # the old summary read the raw stage table
    print(f"   group diff summary before: n {old.size}, max |d| {np.abs(old).max():.3g}, "
          f"mean |d| {np.abs(old).mean():.3g}")
    summary = diff.packages.sfr.results.stage.summary().iloc[0]
    print(f"   {'':>19} after: n {summary['n']}, n_one_sided {summary['n_one_sided']}, "
          f"max |d| {summary['max_abs_diff']:.3g}, within tolerance "
          f"{summary['within_tolerance']}")

    print("\n3. the targets, unchanged (the user's call)")
    file_route = mf.SfrStageTargets(locations={"a": 0, "b": 3}).simulated_series(dry)
    print(f"   stage-file route, reach 0: {_show(file_route['a'])}")
    print(f"   stage-file route, reach 3: {_show(file_route['b'])}")
    obs_route = mf.SfrStageTargets(locations={"r0": 0, "r3": 3}).simulated_series(dry)
    print(f"   observation route, reach 0: {_show(obs_route['r0'])} (bed top {BED_TOP[0]})")
    print(f"   observation route, reach 3: {_show(obs_route['r3'])} (bed top {BED_TOP[3]:.4g})")


if __name__ == "__main__":
    main()
