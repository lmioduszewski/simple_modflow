"""A DRY lake is calibrated against its BOTTOM, on both routes (ledger 351).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/dry_lake_bottom.py

Needs `mf6` and `pestpp-ies` on PATH. Measured with MF6 6.8.0.dev0+0418efe4, FloPy
3.10.0, pyEMU 1.4.0 and PESTPP-IES 5.2.16. About 80 s.

MF6's two outputs disagree at a dry lake. The stage file writes -1e30 once the stage is
at or below the lake's bottom (`gwf-lak.f90`, `lak_ot_dv`); the STAGE observation
writes the stage the solver solved (`lak_bd_obs`), which it never lets below the bottom
(`xnewpak` is clamped at `lakebot`) -- so a dry lake's observation IS its bottom. MF6
writes that bottom nowhere else. `LakeStageTargets` read the observation on one route
and the masked stage file on the other, so a dry period was the bottom on one and `NaN`
on the other -- and 1e30 in the file PEST read.

PART 1 -- `package_tables.lake_bottoms` (MF6's `lakebot` rule, `gwf-lak.f90`
`lak_read_lake_connections`: a VERTICAL connection's cell top, a HORIZONTAL one's BELEV
or its cell's bottom where BELEV == TELEV, the lowest over the lake's connections)
against MF6's STAGE observation at every record the stage file wrote -1e30, on three
runs of six 10-day periods (one step each), every lake with a STAGE observation:
  vertical    6 x 6 irregular Voronoi, one confined layer (top 100, one cell
              99.123456789123), CHD 50. `one`: one VERTICAL connection. `two`: two,
              under tops 100 and 99.123456789123. Both drawn down by a withdrawal and
              refilled by rain in period 4. `off`: STATUS INACTIVE throughout.
  horizontal  the same grid in two confined layers (100/60/40), CHD 30. `hneq`: one
              HORIZONTAL connection into a layer-1 cell, BELEV 45 < TELEV 58; `heq`:
              BELEV = TELEV = 55. No area, so no water: dry in every record.
  excavated   12 x 12 grid of 20 ft cells, a 100/80/40 stack, a pit of radius 30 ft
              excavated to 82.3456789012345 (`LayerStack.excavate`) under `mf.lak`
              (ledger 325's rule), NEWTON, CHD 60 in the lower layer.
  implicit    the excavated run with `mf.lak(implicit=True)`, MF6's LAK IMPLICIT
              (ledger 326): the stage is a matrix unknown, held at the bottom by the
              SAME clamp -- `lak_cq` re-runs `lak_solve(update=.false.)`, which raises
              `xnewpak` to `lakebot` in the solution vector the implicit stage lives in
              (`gwf-lak.f90`, lak_cq, and lak_solve's `xnewpak < lakebot` guards), before
              either output is written (read in the source). `lak_nur` acts only
              under NEWTON UNDER_RELAXATION, and this run's name file reads NEWTON.
Each on the model LOADED from its files and on the one still in memory from its build.

PART 2 -- the two routes and what compare() and PEST are handed, BEFORE (the old rule
applied by hand: the masked table, `NaN` for a dry lake, 1e30 in the PEST file) and
AFTER, on the vertical run: a gauged lake that never ran dry,
[100.45, 100.30, 100.20, 100.10, 100.35, 100.40], against lake `one`.

PART 3 -- what PESTPP-IES does with a dry period written as 1e30 against the bottom.
Not MF6: a one-parameter "lake", stage_t = 99.5 + p * c_t with c = (1.0, 0.8, 0.6) and
bottom 100 (a stage at or below it is dry, written as 1e30 or as 100.0); truth p = 1
gives observations 100.5, 100.3, 100.1 at weight 10. 40 realizations drawn from
N(0.7, 0.3) clipped to [0.05, 2.0] (seed 7), `ies_no_noise`, NOPTMAX 4.

MEASURED
  PART 1  every dry record: lake_bottoms on the LOADED model == MF6's observation, to
          the bit -- vertical one 100.0 (records 1-3), two 99.12345679 (record 3; the
          LOWER top, as FloPy WROTE it), horizontal hneq 45.0 and heq 40.0 (the cell's
          bottom, not the 55 written; records 0-5), excavated and implicit 82.3456789
          (records 1-3).
          The model still in memory reads its unrounded elevations: two 99.123456789123
          (-8.8e-10 from MF6's), excavated 82.3456789012345 (+1.2e-9).
  PART 2  routes AFTER: equal at every record, dry and wet (max |diff| 0.0). BEFORE:
          the stage-file route NaN at records 1-3. compare()/stats() BEFORE: n 6 (rmse
          0.6960) on the observation route, n 3 (rmse 0.9603) on the stage-file route;
          AFTER: n 6 and rmse 0.6960 on both. The PEST file on the stage-file route:
          1e30 at records 1-3 BEFORE, 100.0 AFTER; only lake `off` is 1e30.
  PART 3  bottom: phi mean 22.66 -> 2.295 -> 0.1426 -> 0.1059 -> 0.08456, posterior p
          0.9975 +/- 0.0207, the .rec's lambda 1 -> 1 -> 0.075 -> 0.0056 -> 0.00042.
          1e30: phi mean 1.725e62 -> 7.5e60 -> 1.111 -> 1.111 -> 1.111, posterior p
          0.9630 +/- 0.0655, lambda 1e61 -> 1e61 -> 7.5e59 -> 5.6e58 -> 4.2e57: set by
          the 1e60 phi of the first iteration, it never came back down once every
          realization was wet, so iterations 3-4 moved nothing. It does not fail; it
          stops calibrating.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import shutil
import subprocess
import sys
import tempfile
import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.package_tables import lake_bottoms

NPER = 6
GAUGED = [100.45, 100.30, 100.20, 100.10, 100.35, 100.40]


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    return grid, props, [int(c) for c in edge], interior


def _observations(names):
    return mf.LakeStageTargets(locations=names).to_flopy_obs(csv_name="lake_stage_targets.csv")


def _run(root, model, packages=()):
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 1, 1.0)] * NPER), *packages])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run


def vertical(root):
    names = {"one": 0, "two": 1, "off": 2}
    grid, props, edge, interior = _grid()
    c0, c1, c2, c3 = interior[0], interior[5], interior[6], interior[10]
    top = np.full(grid.ncpl, 100.0)
    top[c2] = 99.123456789123
    lake = mf.lak.flopy(
        nlakes=3,
        packagedata=[[0, 100.5, 1, "one"], [1, 100.5, 2, "two"], [2, 104.0, 1, "off"]],
        connectiondata=[
            [0, 0, (0, c0), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 0, (0, c1), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 1, (0, c2), "VERTICAL", 0.05, 0.0, 0.0, 0.0, 0.0],
            [2, 0, (0, c3), "VERTICAL", 0.01, 0.0, 0.0, 0.0, 0.0],
        ],
        perioddata={
            0: [[0, "WITHDRAWAL", 300.0], [1, "WITHDRAWAL", 300.0], [2, "STATUS", "INACTIVE"]],
            4: [[0, "RAINFALL", 0.5], [1, "RAINFALL", 0.5]],
            5: [[0, "RAINFALL", 0.0], [1, "RAINFALL", 0.0]],
        },
        boundnames=True, stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
        observations=_observations(names),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=top, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=50.0), mf.chd(stress_period_data={0: [[(0, c), 50.0] for c in edge]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, model), names


def horizontal(root):
    names = {"hneq": 0, "heq": 1}
    grid, props, edge, interior = _grid()
    lake = mf.lak.flopy(
        nlakes=2,
        packagedata=[[0, 57.0, 1, "hneq"], [1, 57.0, 1, "heq"]],
        connectiondata=[
            [0, 0, (1, interior[1]), "HORIZONTAL", 0.5, 45.0, 58.0, 50.0, 100.0],
            [1, 0, (1, interior[9]), "HORIZONTAL", 0.5, 55.0, 55.0, 50.0, 100.0],
        ],
        perioddata={0: [[0, "WITHDRAWAL", 3000.0], [1, "WITHDRAWAL", 3000.0]]},
        boundnames=True, stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
        observations=_observations(names),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=[60.0, 40.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=30.0),
        mf.chd(stress_period_data={0: [[(k, c), 30.0] for c in edge for k in (0, 1)]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, model), names


def excavated(root, *, implicit=False):
    names = {"pit": 0}
    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs,
                            geometry=[Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 64)])
    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(80.0))
             .add("lower", bottom=mf.Flat(40.0)))
    result = stack.excavate(bottom=mf.Flat(82.3456789012345), within=pond, name="pit").build(grid)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    lake = mf.lak(context=context, nper=NPER, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=0.5, unexcavated="raise", implicit=implicit,
                  withdrawals=[5000.0, 5000.0, 5000.0, 0.0, 0.0, 0.0],
                  rainfall=[0.0, 0.0, 0.0, 0.0, 2.0, 0.0],
                  stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
                  observations=_observations(names))
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, transient={0: True}),
        mf.ic(strt=60.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 60.0] for c in edge if result.idomain[1, c] > 0]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, model, [mf.ims(models=["m"], complexity="COMPLEX",
                                     linear_acceleration="BICGSTAB")]), names


def part_one(runs):
    print("PART 1 -- lake_bottoms against MF6's STAGE observation at every dry record")
    for label, (run, names) in runs.items():
        raw = run.model().outputs.lak.stage.get(masked=False)
        observed = pd.read_csv(run.workspace / "lake_stage_targets.csv", float_precision="round_trip")
        loaded = lake_bottoms(mf.load_mf6_run(run.workspace, model_name="m", verbosity_level=0))
        built = lake_bottoms(run.model())
        for name, lake in names.items():
            dry = np.flatnonzero(raw[:, lake] == -1e30)
            if not len(dry):
                state = "inactive throughout" if (raw[:, lake] == 1e30).all() else "never dry"
                print(f"  {label:<10} {name:<5} {state}")
                continue
            mf6 = observed[name.upper()].to_numpy()[dry]
            print(f"  {label:<10} {name:<5} dry records {dry.tolist()}: MF6 observation "
                  f"{sorted(set(mf6.tolist()))}; lake_bottoms loaded {float(loaded[lake])!r} "
                  f"(equal: {bool((mf6 == loaded[lake]).all())}); in memory "
                  f"{float(built[lake])!r} ({built[lake] - mf6[0]:+.2g})")


def _old_stage_file_route(model, names):
    """The rule before ledger 351, by hand: the MASKED stage file, NaN for a dry lake."""

    masked = model.outputs.lak.stage.get()
    return pd.DataFrame({"time": range(NPER), **{name: masked[:, lake] for name, lake in names.items()}})


def _stats(simulated, gauged):
    residual = np.asarray(simulated, dtype=float) - np.asarray(gauged)
    kept = residual[~np.isnan(residual)]
    return len(kept), float(np.sqrt(np.mean(kept ** 2)))


def part_two(run, names):
    print("\nPART 2 -- the two routes, compare()/stats() and the PEST file (vertical run)")
    model = mf.load_mf6_run(run.workspace, model_name="m", verbosity_level=0)
    file_names = {f"f_{name}": lake for name, lake in names.items()}
    observed = mf.LakeStageTargets(locations=names).simulated_series(model)
    from_file = mf.LakeStageTargets(locations=file_names).simulated_series(model)
    before = _old_stage_file_route(model, names)
    for name in ("one", "two"):
        after_diff = np.nanmax(np.abs(observed[name].to_numpy() - from_file[f"f_{name}"].to_numpy()))
        print(f"  {name}: observation route {np.round(observed[name].to_numpy(), 4).tolist()}")
        print(f"  {'':<5}stage file BEFORE  {np.round(before[name].to_numpy(), 4).tolist()}")
        print(f"  {'':<5}stage file AFTER   {np.round(from_file[f'f_{name}'].to_numpy(), 4).tolist()}"
              f"  (max |obs - file| {after_diff})")
    print(f"  off: observation route all NaN {observed['off'].isna().all()}, "
          f"stage file all NaN {from_file['f_off'].isna().all()}")

    for route, name in (("observation", "one"), ("stage file", "f_one")):
        targets = mf.LakeStageTargets(locations={name: 0},
                                      values=pd.DataFrame({"time": range(NPER), name: GAUGED}))
        stats = targets.stats(model).iloc[0]
        print(f"  stats AFTER, {route:<11} n {int(stats['n'])} rmse {stats['rmse']:.4f}")
    n, rmse = _stats(before["one"], GAUGED)
    print(f"  stats BEFORE, stage file  n {n} rmse {rmse:.4f}")
    written_before = before.fillna(mf.LakeStageTargets._PEST_NO_STAGE)
    written_after = from_file.fillna(mf.LakeStageTargets._PEST_NO_STAGE)
    print(f"  PEST file, stage-file route, lake one: BEFORE {written_before['one'].tolist()}")
    print(f"  {'':<38}AFTER  {written_after['f_one'].round(4).tolist()}")
    print(f"  1e30 AFTER only in: {[c for c in written_after.columns if (written_after[c] == 1e30).any()]}")


IES_MODEL = '''
import numpy as np
p = float(open("par.dat").read().split()[0])
stage = 99.5 + p * np.array([1.0, 0.8, 0.6])
rows = ["time,lake"]
for t, v in enumerate(stage):
    if v <= 100.0:
        v = {dry}
    rows.append(f"{{t}},{{float(v)!r}}")
open("out.csv", "w").write("\\n".join(rows) + "\\n")
'''


def _ies(root, dry_value, *, reals=40, noptmax=4):
    import pyemu

    ws = root / f"ies_{dry_value}"
    ws.mkdir(parents=True)
    (ws / "model.py").write_text(IES_MODEL.format(dry=dry_value))
    (ws / "run.sh").write_text(f'#!/bin/sh\nexec "{sys.executable}" model.py\n')
    (ws / "par.dat.tpl").write_text("ptf ~\n~ p      ~\n")
    (ws / "par.dat").write_text("0.7\n")
    subprocess.run([sys.executable, "model.py"], cwd=ws, check=True)
    with contextlib.redirect_stdout(io.StringIO()):
        pyemu.pst_utils.csv_to_ins_file(str(ws / "out.csv"), ins_filename=str(ws / "out.csv.ins"),
                                        only_cols=["lake"], prefix="lk")
        pst = pyemu.Pst.from_io_files(str(ws / "par.dat.tpl"), str(ws / "par.dat"),
                                      str(ws / "out.csv.ins"), str(ws / "out.csv"), pst_path=".")
    par = pst.parameter_data
    par.loc[:, ["parval1", "parlbnd", "parubnd"]] = [0.7, 0.05, 2.0]
    par.loc[:, "partrans"] = "none"
    obs = pst.observation_data
    obs.loc[:, "obsval"] = 99.5 + np.array([1.0, 0.8, 0.6])
    obs.loc[:, "weight"] = 10.0
    pst.model_command = ["sh run.sh"]
    pst.control_data.noptmax = noptmax
    draws = np.clip(np.random.default_rng(7).normal(0.7, 0.3, reals), 0.05, 2.0)
    pd.DataFrame({"p": draws}, index=[str(i) for i in range(reals)]).to_csv(ws / "prior.csv")
    pst.pestpp_options.update({"ies_par_en": "prior.csv", "ies_num_reals": reals,
                               "ies_no_noise": True})
    with contextlib.redirect_stdout(io.StringIO()):
        pst.write(str(ws / "t.pst"), version=2)
    subprocess.run(["pestpp-ies", "t.pst"], cwd=ws, check=True, capture_output=True, text=True)
    phi = pd.read_csv(ws / "t.phi.actual.csv")["mean"].tolist()
    posterior = pd.read_csv(ws / f"t.{noptmax}.par.csv", index_col=0)["p"]
    lambdas = [line.split(":")[-1].strip() for line in (ws / "t.rec").read_text().splitlines()
               if "current lambda" in line]
    return phi, posterior, lambdas


def part_three(root):
    print("\nPART 3 -- PESTPP-IES, a dry period written as the bottom or as 1e30")
    for label, value in (("bottom", "100.0"), ("1e30", "1e30")):
        phi, posterior, lambdas = _ies(root, value)
        print(f"  {label:<6} phi mean by iteration {[float(f'{v:.4g}') for v in phi]}")
        print(f"  {'':<6} posterior p {posterior.mean():.4f} +/- {posterior.std():.4f}; "
              f"lambda {lambdas}")


def main():
    root = pathlib.Path(tempfile.mkdtemp(prefix="dry_lake_bottom_"))
    try:
        runs = {"vertical": vertical(root / "vertical"), "horizontal": horizontal(root / "horizontal"),
                "excavated": excavated(root / "excavated"),
                "implicit": excavated(root / "implicit", implicit=True)}
        part_one(runs)
        part_two(*runs["vertical"])
        part_three(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
