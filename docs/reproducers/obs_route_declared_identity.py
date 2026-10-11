"""The observation route reads the package's DECLARED observations, by identity.

Ledger 344. ``LakeStageTargets``, ``SfrStageTargets``, ``SfrFlowTargets`` and
``DrnFlowTargets`` chose their observation route, and read it, by the SPELLING of a
target's name: the first ``*.csv`` anywhere in the model workspace, in sorted order,
whose header held every target's label. The observation's TYPE and its lake, reach
or cell were never consulted, and the file need not belong to the package -- or to
this run. PEST's build and every forward run call the same ``simulated_series``.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/obs_route_declared_identity.py

It uses the public API (plus the forward run's own writer, which a PEST template
calls), so it measures whichever tree is on ``PYTHONPATH`` -- the "before" column
below is this script run against the parent commit's ``src`` (``git archive f0fe9c7
src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``). Needs ``mf6``.

THE SETUP. The 6 x 6 irregular Voronoi grid of ``tests/test_flow_targets_period_end.py``
(``ncpl`` 36), one confined layer, CHD at 102 ft on the west and east edges, heads
filling from 90 ft, three 10-day periods of one step.

``declared``  Lake 0 ``drying`` (cells 8, 9; starts at 111 ft on a 110 ft bed, a
              900 withdrawal empties it from period 1) and lake 1 ``held`` (cells 25,
              26; CONSTANT 118 ft), boundnames on. A four-reach stream on cells
              13-16 (reach 3 the outlet) and drains on cells 20, 21, 27. Declared:
              LAK  lake_stage_targets.csv  drying STAGE 1, held STAGE 2,
                                           bn_held STAGE held, pond VOLUME 1
              SFR  sfr_stage_targets.csv   up STAGE 1, gauge STAGE 4, down STAGE 4,
                                           and two 42-character names sharing 40
                                           characters (STAGE 1, STAGE 4)
                   sfr_flow_targets.csv    gauge DOWNSTREAM-FLOW 4, up
                                           DOWNSTREAM-FLOW 1, xchg SFR 4
              DRN  drn_flow_targets.csv    seep_c20 / _c21 / _c27 DRN (0, cell)
              GWF  heads.csv               held HEAD (0, 10) -- sorts before the
                                           LAK file (a ``PostBuildHook``)
``reused``    ``drying`` alone, run as ``Project.run("r")`` with its STAGE
              observation, then ``Project.run("r", overwrite=True)`` WITHOUT it and
              with the withdrawal cut to 300: the directory keeps the first CSV.
``two``       Two LAK packages ``lak_a`` (cells 8, 9) and ``lak_b`` (cells 25, 26),
              one lake each, each declaring ``pond STAGE 1`` (a.csv, b.csv); two DRN
              packages ``drn_a``/``drn_b`` on cell 20, each declaring ``s_c20 DRN
              (0, 20)`` (da.csv, db.csv).

  A. ``{"drying": 0, "held": 1}`` (correct); ``{"drying": 1}`` (another lake's id).
  B. ``{"drying": 0, "a": 1}``: one name declared, one not.
  C. ``{"held": 1}``: a GWF HEAD observation ``held`` sorts first.
  D. ``{"pond": 0}``: ``pond`` is a VOLUME; ``{"bn_held": 1}`` and ``{"bn_held": 0}``:
     a boundname id.
  E. ``to_flopy_obs`` of ``{"Pond A": 0, "pond-a": 1}``, and with ``kind="VOLUME"``.
  F. SFR stage ``{"gauge": 3}`` (a DOWNSTREAM-FLOW under the same name sorts first),
     ``{"up": 3}`` (declared for reach 0), the two 42-character names; SFR flow
     ``{"xchg": 3}`` (an ``SFR`` exchange observation) and ``{"gauge": 3, "up": 0}``;
     the DRN zone ``seep``.
  G. A copy of ``declared``'s workspace holding ``archive_lake.csv`` -- the LAK file
     with every stage 20 ft lower: ``{"drying": 0, "held": 1}``, and the PEST
     forward-run writer on that copy.
  H. ``reused``: ``{"drying": 0}``, against the second run's own stage file.
  I. ``two``: ``{"pond": 0}``, then with ``package="lak_b"``; ``{"x": 0}`` (no label
     declared); DRN ``{"s": [20]}``, then ``package="drn_b"``.
  J. A copy of ``declared`` whose SFR stage file holds MF6's no-data value, 3e30,
     for ``down`` (reach 3) in period 1: ``{"down": 3}``. (MF6 6.8 writes an
     INACTIVE reach's STAGE as its bed top -- 94.7 ft for a reach with STATUS
     INACTIVE, measured, the stage file writing 1e30 -- not 3e30, so this copy is
     edited, not run.)
  K. (ledger 345) A copy of ``declared`` whose LAK OBS writes ``FILEOUT
     lake_stage_targets.bsv BINARY``, RUN by MF6 there: ``{"drying": 0, "held": 1}``.
  L. (ledger 345) ``two``'s LAKs renamed ``lak`` and ``lak_2`` -- the first NAMED
     exactly its type, as a default-named package is: ``{"pond": 0}``, ``{"x": 0}``
     (nothing declared), ``{"pond": 0}`` with ``package="lak_2"``.

MEASURED (mf6 6.8.0.dev0, 2026-09-29; "before" is f0fe9c7, K and L's "after" ledger 345's tree)
                                             before                   after
  A  {drying: 0, held: 1}                    t 10/20/30: drying       the same, to the bit
                                             110.239 / 110.0 / 110.0
                                             (dry: the bed), held 118
     {drying: 1}                             lake 0's 110.239 / 110.0 ValueError: 'drying'
                                             / 110.0, under lake 1    (lake 1) is declared
                                                                      for lake 0, `drying
                                                                      STAGE 1`
  B  {drying: 0, a: 1}                       BOTH from the stage      ValueError naming
                                             file, t 0/1/2: drying    'drying' (declared)
                                             110.239 / NaN / NaN,     and 'a' (not)
                                             a 118 -- drying moved
                                             off its own route
  C  {held: 1}                               the HEAD, 100.239 /      118 x3, the lake
                                             101.880 / 102.211
  D  {pond: 0}                               the VOLUME, 4861.62 /    ValueError: declared
                                             0 / 0, as a stage        as VOLUME, not STAGE
     {bn_held: 1}                            118 x3                   118 x3
     {bn_held: 0}                            118 x3, as lake 0        ValueError: declared
                                                                      for lake 1
  E  {Pond A, pond-a}.to_flopy_obs           Pond_A, pond_a written   ValueError naming both
     kind="VOLUME"                           ('drying', 'VOLUME', 1)  ValueError
  F  stage {gauge: 3}                        -89.669 / -155.999 /     100.786 / 101.986 /
                                             -104.024: the            102.536, reach 3's
                                             DOWNSTREAM-FLOW          stage (= its stage file)
     stage {up: 3}                           -145.345 / ...: reach    ValueError: declared
                                             0's DOWNSTREAM-FLOW      for reach 0
     the two 42-character names              the stage FILE, t 0/1/2  the observation route,
                                             (a cut label is found    t 10/20/30, each its
                                             in no header)            own reach's values
     flow {xchg: 3}                          34.53 / -23.346 /        ValueError: declared
                                             117.849: the EXCHANGE    as SFR, not
                                                                      DOWNSTREAM-FLOW
     flow {gauge: 3, up: 0}; DRN {seep}      the same on both trees
  G  {drying, held} beside archive_lake.csv  the ARCHIVE: 90.239,     the live file:
                                             98 (20 ft low)           110.239, 118
     forward run on the copy                 the ARCHIVE              the live file
  H  {drying: 0}, second run                 the FIRST run's CSV,     the second run's stage
                                             t 10/20/30: 110.239 /    file, t 0/1/2: 110.365 /
                                             110.0 / 110.0            110.093 / NaN
  I  {pond: 0}, two LAK                      a.csv (lak_a), by sort   ValueError naming
                                             order: 111.693 / ...     'lak_a', 'lak_b'
     package="lak_b"                         TypeError: no such       b.csv: 112.549 / ...
                                             parameter
     {x: 0}, two LAK                         AttributeError: 'list'   ValueError naming both
                                             object has no attribute  and package=
                                             'nlakes'
     DRN {s: [20]}, two DRN                  da.csv, by sort order:   ValueError naming
                                             -205.856 / ...           'drn_a', 'drn_b'
     package="drn_b"                         TypeError                db.csv: -155.856 / ...
  J  {down: 3}, 3e30 in period 1             0.0 ft                   NaN
  K  {drying, held}, BINARY (mf6 rc 0,       the stage FILE, t 0/1/2  t 10/20/30: 110.239 /
     the file opens b'cont double')          (not a *.csv); 344's     110.0 / 110.0, 118 --
                                             commit (09a61a4):        the CSV run's values
                                             UnicodeDecodeError       to the bit
  L  {pond: 0}, LAK 'lak' + 'lak_2'          a.csv ('lak') by sort    ValueError naming
                                             order: 111.693 / ...     'lak', 'lak_2'
     {x: 0}, LAK 'lak' + 'lak_2'             'lak''s stage file, in   ValueError naming both
                                             silence (FloPy answers   and package= (344 and
                                             get_package("lak") with  after: found by TYPE)
                                             the package NAMED lak)
     package="lak_2"                         TypeError                b.csv: 112.549 / ...
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
DRYING, HELD = [8, 9], [25, 26]
LONG_UP = "a" * 39 + "_up"     # 42 characters; MF6 keeps 40, the same for both
LONG_DOWN = "a" * 39 + "_dn"


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    return grid, props, [int(c) for c in edge]


def _base(grid, props, edge):
    return [
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
    ]


def _oc():
    return mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
                 saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])


def _head_obs(model, packages, context):
    """A GWF HEAD observation named ``held``, in ``heads.csv`` (sorts before the LAK file)."""

    import flopy

    return flopy.mf6.ModflowUtlobs(model, pname="head_obs", filename="m.head.obs",
                                   continuous={"heads.csv": [("held", "HEAD", (0, 10))]})


def _lake(*, name="lak", lakes=((DRYING, 111.0, "drying", 900.0), (HELD, 118.0, "held", None)),
          observations=None):
    connections, packagedata, perioddata = [], [], []
    for number, (cells, strt, boundname, withdrawal) in enumerate(lakes):
        connections += [[number, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(cells)]
        packagedata.append([number, strt, len(cells), boundname])
        if withdrawal is None:
            perioddata += [[number, "STATUS", "CONSTANT"], [number, "STAGE", strt]]
        else:
            perioddata.append([number, "WITHDRAWAL", withdrawal])
    options = {"observations": observations} if observations else {}
    return mf.lak.flopy(
        name=name, filename=f"m.{name}", nlakes=len(lakes), packagedata=packagedata,
        connectiondata=connections, perioddata={0: perioddata}, boundnames=True,
        stage_filerecord=f"m.{name}.stg", budget_filerecord=f"m.{name}.bud", **options)


def _stream(observations):
    last = len(REACHES) - 1
    packagedata = [
        [i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
         1 if i in (0, last) else 2, 1.0, 0]
        for i, c in enumerate(REACHES)
    ]
    connectiondata = [
        [i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
        for i in range(len(REACHES))
    ]
    return mf.sfr.flopy(
        packagedata=packagedata, connectiondata=connectiondata,
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]], 2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud", observations=observations)


def _drain(name="drn", cells=DRN_CELLS, observations=None, elevation=96.0):
    options = {"observations": observations} if observations else {}
    return mf.drn.flopy(name=name, filename=f"m.{name}",
                        stress_period_data={0: [[(0, c), elevation, 50.0] for c in cells]}, **options)


def _run(root, packages, *, overwrite=False, project=None, hooks=()):
    grid, props, edge = _grid()
    model = mf.gwf("m", grid=grid, packages=[*_base(grid, props, edge), *packages], hooks=hooks)
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=3, perioddata=[(10.0, 1, 1.0)] * 3)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        project = project or mf.Project(root / "p", name="p")
        run = project.run("r", simulation, overwrite=overwrite)
    assert run.success, "\n".join(run.report[-20:])
    return run.model(), project


def declared(root):
    lak_obs = {"lake_stage_targets.csv": [
        ("drying", "STAGE", 1), ("held", "STAGE", 2), ("bn_held", "STAGE", "held"), ("pond", "VOLUME", 1)]}
    sfr_obs = {
        "sfr_stage_targets.csv": [("up", "STAGE", 1), ("gauge", "STAGE", 4), ("down", "STAGE", 4),
                                  (LONG_UP, "STAGE", 1), (LONG_DOWN, "STAGE", 4)],
        "sfr_flow_targets.csv": [("gauge", "DOWNSTREAM-FLOW", 4), ("up", "DOWNSTREAM-FLOW", 1),
                                 ("xchg", "SFR", 4)],
    }
    drn_obs = {"drn_flow_targets.csv": [(f"seep_c{c}", "DRN", (0, c)) for c in DRN_CELLS]}
    model, _ = _run(root, [_lake(observations=lak_obs), _stream(sfr_obs),
                           _drain(observations=drn_obs), _oc()],
                    hooks=[mf.PostBuildHook("head_obs", _head_obs)])
    return model


def _show(frame):
    if isinstance(frame, pd.DataFrame):
        return "; ".join(
            f"{column}={np.round(frame[column].to_numpy(dtype=float), 3).tolist()}" for column in frame.columns)
    return repr(frame)


def _attempt(label, what, call):
    try:
        result = call()
    except Exception as exc:  # noqa: BLE001 - a reproducer reports whatever each tree raises
        result = f"{type(exc).__name__}: {str(exc)[:400]}"
    else:
        result = _show(result)
    print(f"{label}  {what}\n     {result}")


def main():
    warnings.simplefilter("ignore", FutureWarning)
    root = pathlib.Path(tempfile.mkdtemp(prefix="obs_route_identity_"))
    try:
        model = declared(root / "declared")
        L, SS, SF, D = mf.LakeStageTargets, mf.SfrStageTargets, mf.SfrFlowTargets, mf.DrnFlowTargets
        print("stage file lake 0/1:", np.round(model.outputs.lak.stage.get(), 3).tolist())
        _attempt("A", "{drying: 0, held: 1}", lambda: L(locations={"drying": 0, "held": 1}).simulated_series(model))
        _attempt("A", "{drying: 1}", lambda: L(locations={"drying": 1}).simulated_series(model))
        _attempt("B", "{drying: 0, a: 1}", lambda: L(locations={"drying": 0, "a": 1}).simulated_series(model))
        _attempt("C", "{held: 1}", lambda: L(locations={"held": 1}).simulated_series(model))
        _attempt("D", "{pond: 0}", lambda: L(locations={"pond": 0}).simulated_series(model))
        _attempt("D", "{bn_held: 1}", lambda: L(locations={"bn_held": 1}).simulated_series(model))
        _attempt("D", "{bn_held: 0}", lambda: L(locations={"bn_held": 0}).simulated_series(model))
        _attempt("E", "{Pond A: 0, pond-a: 1}.to_flopy_obs()",
                 lambda: L(locations={"Pond A": 0, "pond-a": 1}).to_flopy_obs())
        _attempt("E", "{drying: 0}.to_flopy_obs(kind='VOLUME')",
                 lambda: L(locations={"drying": 0}).to_flopy_obs(kind="VOLUME"))
        print("   SFR stage file reach 3:",
              np.round(model.sfr.results.stage.get().query("reach == 3")["stage"].to_numpy(), 3).tolist())
        _attempt("F", "stage {gauge: 3}", lambda: SS(locations={"gauge": 3}).simulated_series(model))
        _attempt("F", "stage {up: 3}", lambda: SS(locations={"up": 3}).simulated_series(model))
        _attempt("F", "stage, two 42-character names",
                 lambda: SS(locations={LONG_UP: 0, LONG_DOWN: 3}).simulated_series(model))
        _attempt("F", "flow {xchg: 3}", lambda: SF(locations={"xchg": 3}).simulated_series(model))
        _attempt("F", "flow {gauge: 3, up: 0}", lambda: SF(locations={"gauge": 3, "up": 0}).simulated_series(model))
        _attempt("F", "drn {seep: cells}", lambda: D(locations={"seep": DRN_CELLS}).simulated_series(model))

        # G: an archived copy beside the live file.
        copy = root / "copy"
        shutil.copytree(model.workspace, copy)
        live = pd.read_csv(copy / "lake_stage_targets.csv")
        live[["DRYING", "HELD"]].sub(20.0).assign(time=live["time"]).loc[:, ["time", "DRYING", "HELD"]].to_csv(
            copy / "archive_lake.csv", index=False)
        reloaded = mf.load_mf6_run(copy, verbosity_level=0)
        _attempt("G", "{drying: 0, held: 1} beside archive_lake.csv",
                 lambda: L(locations={"drying": 0, "held": 1}).simulated_series(reloaded))

        def forward():
            from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets

            pd.DataFrame({"name": ["drying", "held"], "lake": [0, 1]}).to_csv(copy / "loc.csv", index=False)
            write_named_series_targets("lake_stage", copy / "loc.csv", copy / "fw.csv", sim_ws=copy)
            return pd.read_csv(copy / "fw.csv")

        _attempt("G", "forward-run writer on the copy", forward)

        # H: a reused run directory.
        lak_obs = {"lake_stage_targets.csv": [("drying", "STAGE", 1)]}
        _first, project = _run(root / "reused", [_lake(lakes=((DRYING, 111.0, "drying", 900.0),),
                                                       observations=lak_obs), _oc()])
        second, _ = _run(root / "reused", [_lake(lakes=((DRYING, 111.0, "drying", 300.0),)), _oc()],
                         overwrite=True, project=project)
        print("   second run's stage file:", np.round(second.outputs.lak.stage.get()[:, 0], 3).tolist())
        _attempt("H", "{drying: 0}, second run", lambda: L(locations={"drying": 0}).simulated_series(second))

        # I: two packages of one type.
        two, _ = _run(root / "two", [
            _lake(name="lak_a", lakes=((DRYING, 114.0, "pa", 100.0),), observations={"a.csv": [("pond", "STAGE", 1)]}),
            _lake(name="lak_b", lakes=((HELD, 116.0, "pb", 100.0),), observations={"b.csv": [("pond", "STAGE", 1)]}),
            _drain("drn_a", [20], {"da.csv": [("s_c20", "DRN", (0, 20))]}, elevation=96.0),
            _drain("drn_b", [20], {"db.csv": [("s_c20", "DRN", (0, 20))]}, elevation=97.0),
            _oc(),
        ])
        for name in ("a.csv", "b.csv", "da.csv", "db.csv"):
            print(f"   {name}:", pd.read_csv(pathlib.Path(two.workspace) / name).iloc[:, 1].round(3).tolist())
        _attempt("I", "{pond: 0}", lambda: L(locations={"pond": 0}).simulated_series(two))
        _attempt("I", "{pond: 0}, package='lak_b'", lambda: L(locations={"pond": 0}, package="lak_b").simulated_series(two))
        _attempt("I", "{x: 0}", lambda: L(locations={"x": 0}).simulated_series(two))
        _attempt("I", "DRN {s: [20]}", lambda: D(locations={"s": [20]}).simulated_series(two))
        _attempt("I", "DRN {s: [20]}, package='drn_b'",
                 lambda: D(locations={"s": [20]}, package="drn_b").simulated_series(two))

        # J: MF6's no-data value in a declared SFR stage file.
        nodata = root / "nodata"
        shutil.copytree(model.workspace, nodata)
        text = (nodata / "sfr_stage_targets.csv").read_text().splitlines()
        cells = text[2].split(",")
        cells[3] = "3.0000000000000000E+30"  # period 1, `down` (reach 3)
        text[2] = ",".join(cells)
        (nodata / "sfr_stage_targets.csv").write_text("\n".join(text) + "\n")
        _attempt("J", "{down: 3} with 3e30 in period 1",
                 lambda: SS(locations={"down": 3}).simulated_series(mf.load_mf6_run(nodata, verbosity_level=0)))

        # K: the LAK observations written BINARY, by MF6 (ledger 345).
        binary = root / "binary"
        shutil.copytree(model.workspace, binary)
        obs = (binary / "m.lak.obs").read_text()
        (binary / "m.lak.obs").write_text(
            obs.replace("BEGIN continuous  FILEOUT  lake_stage_targets.csv",
                        "BEGIN continuous  FILEOUT  lake_stage_targets.bsv  BINARY")
            .replace("END continuous  FILEOUT  lake_stage_targets.csv", "END continuous"))
        (binary / "lake_stage_targets.csv").unlink()
        ran = subprocess.run(["mf6"], cwd=binary, capture_output=True, text=True)
        print(f"   mf6 on the BINARY declaration: return code {ran.returncode}, file opens "
              f"{(binary / 'lake_stage_targets.bsv').read_bytes()[:11]!r}")
        _attempt("K", "{drying: 0, held: 1}, BINARY",
                 lambda: L(locations={"drying": 0, "held": 1}).simulated_series(
                     mf.load_mf6_run(binary, verbosity_level=0)))

        # L: two packages, one NAMED like its type (ledger 345).
        named, _ = _run(root / "named", [
            _lake(name="lak", lakes=((DRYING, 114.0, "pa", 100.0),), observations={"a.csv": [("pond", "STAGE", 1)]}),
            _lake(name="lak_2", lakes=((HELD, 116.0, "pb", 100.0),), observations={"b.csv": [("pond", "STAGE", 1)]}),
            _oc(),
        ])
        for name in ("a.csv", "b.csv"):
            print(f"   {name}:", pd.read_csv(pathlib.Path(named.workspace) / name).iloc[:, 1].round(3).tolist())
        print("   m.lak.stg (package 'lak'):",
              np.round(mf.load_mf6_run(named.workspace, verbosity_level=0).outputs.lak(package="lak").stage.get()[:, 0], 3).tolist())
        _attempt("L", "{pond: 0}, LAK 'lak' and 'lak_2'", lambda: L(locations={"pond": 0}).simulated_series(named))
        _attempt("L", "{x: 0}, LAK 'lak' and 'lak_2'", lambda: L(locations={"x": 0}).simulated_series(named))
        _attempt("L", "{pond: 0}, package='lak_2'",
                 lambda: L(locations={"pond": 0}, package="lak_2").simulated_series(named))
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
