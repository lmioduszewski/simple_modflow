"""Every LAK stage view drew MF6's -1e30 (dry) and +1e30 (inactive) as a stage.

Ledger 329. Reported against a real run whose excavated pond went dry on 7 of 38 saved
stages: ``packages.lak.results.stage.plot()`` drew those days at -1e30.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/lak_stage_sentinels.py

THE SETUP. The same toy model ``tests/test_lak_stage_sentinels.py`` runs: a 6 x 6
irregular Voronoi grid of ~100 ft cells, one confined layer (top 100, bottom 0, K 10),
CHD at 50 ft on the west and east edges, six 10-day transient periods, and three
one-cell VERTICAL lakes, so each lake's bottom is its cell's top, 100 ft:

  lake 0  starts 0.5 ft deep, drawn down by a 300 ft3/d withdrawal; dry in periods
          1-3, refilled by 0.5 ft/d of rain in period 4.
  lake 1  STATUS INACTIVE in period 2 only.
  lake 2  CONSTANT at 108 ft: a control.

A second run ("wet") has no withdrawal and lake 1 always active, for the group diff.
The "before" of each view is its OLD rule applied here by hand to the raw file
(``outputs.lak.stage.get(masked=False)``), so this still shows it after the fix.

WHAT IT MEASURES
  1. What MF6 writes, in its stage file and in its STAGE observation.
  2. Each view, before and after: the series, a map's colour range, the stage change,
     the observation route of ``LakeStageTargets``, and the group diff summary.
  3. The group series, which sums, and pandas sums an all-NaN group to 0.0.
  4. The same sum in a map payload (``stage.map(agg="sum")``), and the stage-file route
     of ``LakeStageTargets`` (target names the observation CSV does not hold), which
     dropped the masked rows -- the series PEST reads.
"""

from __future__ import annotations

import pathlib
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

NPER = 6


def _run(root: pathlib.Path, *, withdrawal: float, switch: bool):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    cells = (interior[0], interior[5], interior[-1])
    held = [[2, "STATUS", "CONSTANT"], [2, "STAGE", 108.0]]
    targets = mf.LakeStageTargets(locations={"drying": 0, "switched": 1, "wet": 2})
    lake = mf.lak.flopy(
        nlakes=3,
        packagedata=[[0, 100.5, 1, "drying"], [1, 104.0, 1, "switched"], [2, 108.0, 1, "wet"]],
        connectiondata=[[i, 0, (0, c), "VERTICAL", 0.5 if i == 0 else 0.01, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(cells)],
        perioddata={
            0: [[0, "WITHDRAWAL", withdrawal], [1, "STATUS", "ACTIVE"], *held],
            2: [[1, "STATUS", "INACTIVE" if switch else "ACTIVE"], *held],
            3: [[1, "STATUS", "ACTIVE"], *held],
            4: [[0, "RAINFALL", 0.5], *held],
            5: [[0, "RAINFALL", 0.0], *held],
        },
        boundnames=True, stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
        observations=targets.to_flopy_obs(csv_name="lake_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=50.0), mf.chd(stress_period_data={0: [[(0, int(c)), 50.0] for c in edge]}),
        lake,
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
    root = pathlib.Path(tempfile.mkdtemp(prefix="lak_stage_sentinels_"))
    dry = _run(root / "dry", withdrawal=300.0, switch=True)
    wet = _run(root / "wet", withdrawal=0.0, switch=False)
    raw = dry.outputs.lak.stage.get(masked=False)

    print("1. what MF6 writes (rows = saved records, columns = lakes 0, 1, 2)")
    for per, row in enumerate(raw):
        print(f"   per {per}: {_show(row)}")
    obs = pd.read_csv(pathlib.Path(dry.workspace) / "lake_stage_targets.csv")
    print(f"   STAGE observation, lake 0 (dry in 1-3): {_show(obs['DRYING'])}")
    print(f"   STAGE observation, lake 1 (inactive in 2): {_show(obs['SWITCHED'])}")

    print("\n2. each view, before (old rule by hand) and after")
    series = {t.name: np.asarray(t.y, float) for t in dry.packages.lak.results.stage.plot().data}
    print(f"   series, lake 0   before {_show(raw[:, 0])}")
    print(f"   {'':>16} after {_show(series['Lake 0'])}")
    for per in (1, 2):
        # Unconnected cells take the default fill_value, 0.0, before and after.
        values = np.asarray(dry.packages.lak.results.stage.map(per=per).linear_zs, float)
        old = np.append(raw[per], 0.0)
        print(f"   map(per={per}) colour range before {old.min():.4g} .. {old.max():.4g}; "
              f"after {np.nanmin(values):.4g} .. {np.nanmax(values):.4g}")
    change = dry.packages.lak.results.stage_change.get(lake=0)
    old_change = np.diff(raw[:, 0])
    print(f"   stage change, lake 0  before {_show(old_change)} (rows kept: 5)")
    print(f"   {'':>20} after {_show(change['stage_change'])} (rows kept: {len(change)})")
    masked_rows = pd.DataFrame({"per0": range(NPER), "stage0": dry.outputs.lak.stage.get()[:, 0]})
    masked_rows["stage1"] = masked_rows["stage0"].shift(-1)
    masked_rows["per1"] = masked_rows["per0"].shift(-1)
    print(f"   {'':>9} masked but dropna(stage1) kept: "
          f"{len(masked_rows.dropna(subset=['per1', 'stage1']))} of 5 transitions")
    targets = mf.LakeStageTargets(locations={"drying": 0, "switched": 1, "wet": 2})
    after = targets.simulated_series(dry)["switched"]
    print(f"   LakeStageTargets, lake 1 (obs route)  before {_show(obs['SWITCHED'].mask(obs['SWITCHED'].abs() > 1e20, 0.0))}")
    print(f"   {'':>32} after {_show(after)}")

    diff = mf.ModelGroup({"wet": wet, "dry": dry}, reference="wet").diff()
    table = diff.packages.lak.results.stage.get()
    raw_wet = wet.outputs.lak.stage.get(masked=False)
    old = (raw - raw_wet).ravel()  # the old summary read the raw stage table
    print(f"   group diff summary before: n {old.size}, max |d| {np.abs(old).max():.3g}, "
          f"mean |d| {np.abs(old).mean():.3g}")
    summary = diff.packages.lak.results.stage.summary().iloc[0]
    print(f"   {'':>19} after: n {summary['n']}, n_one_sided {summary['n_one_sided']}, "
          f"max |d| {summary['max_abs_diff']:.3g}, mean |d| {summary['mean_abs_diff']:.3g}, "
          f"within tolerance {summary['within_tolerance']} ({len(table)} rows)")

    print("\n3. the group series sums; an all-NaN group sums to 0.0")
    frame = diff.group.packages.lak.results.stage.get(model_name="dry", lake=0)
    print(f"   pandas sum of the masked lake-0 series: "
          f"{_show(frame.groupby('per')['stage'].agg('sum'))}")
    lines = {t.name: np.asarray(t.y, float)
             for t in diff.group.packages.lak.results.stage.plot().data}
    print(f"   group plot after: {_show(lines['dry / Lake 0'])}")

    print("\n4. a map that sums, and the targets' stage-file route")
    stage = dry.packages.lak.results.stage
    cell = int(stage.get().query("lake == 0")["cell"].iloc[0])  # the drying lake's cell
    table = stage.get(per=1)
    old = table.groupby("cell")["stage"].agg("sum")  # the payload's old aggregation
    print(f"   map(per=1, agg='sum') dry lake's cell before {old.loc[cell]:.4g}, "
          f"colour range {old.min():.4g} .. {old.max():.4g}")
    values = np.asarray(stage.map(per=1, agg="sum", fill_value=float("nan")).linear_zs, float)
    print(f"   {'':>22} after {values[cell]:.4g}, "
          f"colour range {np.nanmin(values):.4g} .. {np.nanmax(values):.4g}")
    for locations in ({"foo": 0}, {"foo": 0, "bar": 1}):
        masked = stage.get()
        kept = masked[masked["lake"].isin(locations.values())].dropna(subset=["stage"])
        series = mf.LakeStageTargets(locations=locations).simulated_series(dry)
        print(f"   stage-file route {sorted(locations)}: rows before "
              f"{kept['per'].nunique()} of {NPER} (dropna on the stage), after {len(series)}")


if __name__ == "__main__":
    main()
