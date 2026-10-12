"""The observation route reads the package's DECLARED observations, by identity (ledger 344).

``LakeStageTargets``, ``SfrStageTargets``, ``SfrFlowTargets`` and ``DrnFlowTargets``
chose their observation route, and read it, by the spelling of a target's name:
the first ``*.csv`` anywhere in the workspace, in sorted order, whose header held
every label. A lake target read another lake's stage, a GWF head or a LAK volume;
an SFR stage target read a streamflow; a stale copy left by an earlier run, or a
backup, froze every PEST realization. They now read what the owning package's OBS
DECLARES -- a label counts only as the family's quantity, for the target's own lake,
reach or cell -- only from the files it names, and a set declared in part, or a
model with several packages of the type and no ``package=``, is refused by name.

Real MF6 runs of the model in ``docs/reproducers/obs_route_declared_identity.py``:
the 6 x 6 irregular Voronoi grid (``ncpl`` 36), one confined layer, CHD at 102 ft on
the west and east edges, three 10-day periods of one step.

``declared``  lake 0 ``drying`` (cells 8, 9; dries from period 1, its bed 110 ft)
              and lake 1 ``held`` (cells 25, 26; CONSTANT 118 ft); a four-reach
              stream on cells 13-16; drains on cells 20, 21, 27. Declared: LAK
              ``drying STAGE 1``, ``held STAGE 2``, ``bn_held STAGE held``,
              ``pond VOLUME 1``; SFR stage ``up STAGE 1``, ``gauge STAGE 4``,
              ``down STAGE 4`` and two 42-character names sharing 40 characters; SFR
              flow ``gauge DOWNSTREAM-FLOW 4``, ``up DOWNSTREAM-FLOW 1``, ``xchg SFR
              4``; DRN ``seep_c<cell> DRN (0, cell)``; and a GWF ``held HEAD``
              in ``heads.csv``, which sorts before the LAK file.
``reused``    ``drying`` alone, run with its STAGE observation and then again in the
              same directory WITHOUT it (``overwrite=True``).
``two``       two LAK packages ``lak_a``/``lak_b`` and two DRN packages
              ``drn_a``/``drn_b``, each declaring the same label.

What each file holds is read with pandas directly, not through the code under test.
"""

from __future__ import annotations

import shutil
import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
from myflopy.modflow.mf6.pest.observations import prepare_lake_stage_observations
from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec

REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
DRYING, HELD = [8, 9], [25, 26]
LONG_UP = "a" * 39 + "_up"     # 42 characters: MF6 keeps 40, so both head "AAA..._"
LONG_DOWN = "a" * 39 + "_dn"
#: The observation route's rows: one per stress period, at its end (ledger 346). The
#: three 10-day periods end at MF6's t = 10, 20 and 30, one time step each.
PERIODS = [0, 1, 2]


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    return grid, props, [int(c) for c in edge]


def _lake(*, name="lak", lakes, observations=None):
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


def _drain(name, cells, observations=None, elevation=96.0):
    options = {"observations": observations} if observations else {}
    return mf.drn.flopy(name=name, filename=f"m.{name}",
                        stress_period_data={0: [[(0, c), elevation, 50.0] for c in cells]}, **options)


def _head_obs(model, packages, context):
    """A GWF HEAD observation named ``held``, in ``heads.csv`` (sorts before the LAK file)."""

    import flopy

    return flopy.mf6.ModflowUtlobs(model, pname="head_obs", filename="m.head.obs",
                                   continuous={"heads.csv": [("held", "HEAD", (0, 10))]})


def _run(root, packages, *, hooks=(), overwrite=False, project=None):
    grid, props, edge = _grid()
    base = [
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ]
    model = mf.gwf("m", grid=grid, packages=[*base, *packages], hooks=hooks)
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(10.0, 1, 1.0)] * 3)])
    with warnings.catch_warnings():
        # A built grid object (and a hook defined in a test module) cannot be frozen
        # into the run's snapshot, and these runs are never reopened from it.
        warnings.simplefilter("ignore", UserWarning)
        project = project or mf.Project(root / "p", name="p")
        run = project.run("r", simulation, overwrite=overwrite)
    assert run.success, "\n".join(run.report[-20:])
    return run.model(), project


@pytest.fixture(scope="module")
def declared(tmp_path_factory):
    """The ``declared`` model of the module docstring, run."""

    lake = _lake(
        lakes=((DRYING, 111.0, "drying", 900.0), (HELD, 118.0, "held", None)),
        observations={"lake_stage_targets.csv": [
            ("drying", "STAGE", 1), ("held", "STAGE", 2), ("bn_held", "STAGE", "held"),
            ("pond", "VOLUME", 1)]},
    )
    stream = _stream({
        "sfr_stage_targets.csv": [("up", "STAGE", 1), ("gauge", "STAGE", 4), ("down", "STAGE", 4),
                                  (LONG_UP, "STAGE", 1), (LONG_DOWN, "STAGE", 4)],
        "sfr_flow_targets.csv": [("gauge", "DOWNSTREAM-FLOW", 4), ("up", "DOWNSTREAM-FLOW", 1),
                                 ("xchg", "SFR", 4)],
    })
    drain = _drain("drn", DRN_CELLS, {"drn_flow_targets.csv": [
        (f"seep_c{c}", "DRN", (0, c)) for c in DRN_CELLS]})
    model, _ = _run(tmp_path_factory.mktemp("declared"), [lake, stream, drain],
                    hooks=[mf.PostBuildHook("head_obs", _head_obs)])
    return model


@pytest.fixture(scope="module")
def two(tmp_path_factory):
    """Two LAK and two DRN packages, each declaring the same label."""

    model, _ = _run(tmp_path_factory.mktemp("two"), [
        _lake(name="lak_a", lakes=((DRYING, 114.0, "pa", 100.0),),
              observations={"a.csv": [("pond", "STAGE", 1)]}),
        _lake(name="lak_b", lakes=((HELD, 116.0, "pb", 100.0),),
              observations={"b.csv": [("pond", "STAGE", 1)]}),
        _drain("drn_a", [20], {"da.csv": [("s_c20", "DRN", (0, 20))]}, elevation=96.0),
        _drain("drn_b", [20], {"db.csv": [("s_c20", "DRN", (0, 20))]}, elevation=97.0),
    ])
    return model


def _csv(model, name):
    """A declared observation file as MF6 wrote it, read with pandas to the last digit.

    ``round_trip``, as the reader does since ledger 346: MF6 writes every digit, and
    pandas' default parser reads some values one unit in the last place off.
    """

    return pd.read_csv(f"{model.workspace}/{name}", float_precision="round_trip")


# -- what a correct declaration reads ------------------------------------------


@pytest.mark.slow
def test_a_correctly_declared_set_reads_its_own_columns(declared):
    """Nothing that read correctly moves: each family's columns equal MF6's file."""

    lakes = _csv(declared, "lake_stage_targets.csv")
    series = mf.LakeStageTargets(locations={"drying": 0, "held": 1}).simulated_series(declared)
    assert list(series.columns) == ["time", "drying", "held"]
    np.testing.assert_array_equal(series["time"], PERIODS)
    np.testing.assert_array_equal(series["drying"], lakes["DRYING"])
    np.testing.assert_array_equal(series["held"], lakes["HELD"])
    np.testing.assert_array_equal(series["held"], 118.0)

    flows = _csv(declared, "sfr_flow_targets.csv")
    series = mf.SfrFlowTargets(locations={"gauge": 3, "up": 0}).simulated_series(declared)
    np.testing.assert_array_equal(series["gauge"], flows["GAUGE"])
    np.testing.assert_array_equal(series["up"], flows["UP"])

    drains = _csv(declared, "drn_flow_targets.csv")
    series = mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).simulated_series(declared)
    np.testing.assert_allclose(series["seep"], drains.iloc[:, 1:].sum(axis=1), rtol=1e-14)


@pytest.mark.slow
def test_a_head_observation_sorting_first_is_not_read_as_a_lake(declared):
    """``heads.csv`` holds a GWF HEAD named ``held`` and sorts before the LAK file:
    the name scan read the aquifer's head (about 101 ft) as the lake held at 118."""

    heads = _csv(declared, "heads.csv")
    assert not np.isclose(heads["HELD"], 118.0).any()
    series = mf.LakeStageTargets(locations={"held": 1}).simulated_series(declared)
    np.testing.assert_array_equal(series["held"], 118.0)


@pytest.mark.slow
def test_a_boundname_id_resolves_to_its_lake(declared):
    """``bn_held STAGE held`` names lake 1 by its boundname."""

    series = mf.LakeStageTargets(locations={"bn_held": 1}).simulated_series(declared)
    np.testing.assert_array_equal(series["bn_held"], 118.0)
    with pytest.raises(ValueError, match=r"'bn_held' \(lake 0.*declared for lake 1.*bn_held STAGE held"):
        mf.LakeStageTargets(locations={"bn_held": 0}).simulated_series(declared)


@pytest.mark.slow
def test_sfr_stage_reads_its_stage_not_a_same_named_flow(declared):
    """``gauge`` is declared twice -- a STAGE and a DOWNSTREAM-FLOW -- and the flow
    file sorts first: the stage target read the (negative) streamflow as a stage."""

    stage = declared.sfr.results.stage.get().query("reach == 3").sort_values("per")["stage"]
    series = mf.SfrStageTargets(locations={"gauge": 3}).simulated_series(declared)
    np.testing.assert_allclose(series["gauge"], stage, rtol=1e-12)
    assert (series["gauge"] > 90.0).all()


@pytest.mark.slow
def test_names_mf6_cuts_to_one_header_read_their_own_columns(declared):
    """Two 42-character names MF6 writes under ONE 40-character header: each is read
    from its own column, by position (the name scan found neither label in the file
    and fell to the stage file -- ledger 343's found-not-fixed, for this route)."""

    header = (declared.workspace / "sfr_stage_targets.csv").read_text().splitlines()[0].split(",")
    assert header[4] == header[5] == "A" * 39 + "_"
    stage = declared.sfr.results.stage.get()
    series = mf.SfrStageTargets(locations={LONG_UP: 0, LONG_DOWN: 3}).simulated_series(declared)
    np.testing.assert_array_equal(series["time"], PERIODS)
    for name, reach in ((LONG_UP, 0), (LONG_DOWN, 3)):
        np.testing.assert_allclose(
            series[name], stage.loc[stage["reach"] == reach].sort_values("per")["stage"], rtol=1e-12)


# -- what is refused -------------------------------------------------------------


@pytest.mark.slow
def test_a_target_naming_another_lakes_observation_raises(declared):
    """``drying`` is declared for lake 0: a target claiming lake 1 read lake 0's stage."""

    for locations in ({"drying": 1}, {"drying": 1, "held": 0}):
        targets = mf.LakeStageTargets(locations=locations, values={"time": [0], "drying": [111.0]})
        for read in (targets.simulated_series, targets.compare):
            with pytest.raises(ValueError, match=(
                r"'drying' \(lake 1, label DRYING\) is declared for lake 0: "
                r"`drying STAGE 1` in lake_stage_targets.csv")):
                read(declared)


@pytest.mark.slow
def test_a_label_declared_as_another_quantity_raises(declared):
    """``pond`` is a LAK VOLUME and ``xchg`` an SFR exchange: each was read as the
    target's quantity (a volume of 4861.62 as a stage)."""

    with pytest.raises(ValueError, match=r"'pond'.*another quantity, not STAGE: `pond VOLUME 1`"):
        mf.LakeStageTargets(locations={"pond": 0}).simulated_series(declared)
    with pytest.raises(ValueError, match=r"'xchg'.*another quantity, not DOWNSTREAM-FLOW: `xchg SFR 4`"):
        mf.SfrFlowTargets(locations={"xchg": 3}).simulated_series(declared)
    with pytest.raises(ValueError, match=r"'up' \(reach 3, label UP\) is declared for reach 0"):
        mf.SfrStageTargets(locations={"up": 3}).simulated_series(declared)


@pytest.mark.slow
def test_a_set_declared_in_part_raises_naming_both_groups(declared):
    """One declared name moved the WHOLE set to the file route, on another clock."""

    with pytest.raises(ValueError, match=r"1 of 2 targets are declared.*\('drying'\) while 'a' is not"):
        mf.LakeStageTargets(locations={"drying": 0, "a": 1}).simulated_series(declared)
    zones = mf.DrnFlowTargets(locations={"seep": DRN_CELLS, "other": [22]})
    with pytest.raises(ValueError, match=r"\('seep'\) while 'other' is not \(labels OTHER_C22\)"):
        zones.simulated_series(declared)
    # A set with NO declared name still reads the file, as before.
    series = mf.LakeStageTargets(locations={"a": 0, "b": 1}).simulated_series(declared)
    np.testing.assert_array_equal(series["time"], [0, 1, 2])


def test_to_flopy_obs_refuses_two_names_one_label_and_another_kind():
    """MF6 accepted ``Pond_A`` and ``pond_a`` and wrote the header ``POND_A,POND_A``;
    a ``kind`` other than the family's quantity was written under the targets' names."""

    with pytest.raises(ValueError, match=r"POND_A from 'Pond A', 'pond-a'"):
        mf.LakeStageTargets(locations={"Pond A": 0, "pond-a": 1}).to_flopy_obs()
    with pytest.raises(ValueError, match=r"GAGE_A from 'Gage A', 'gage-a'"):
        mf.SfrStageTargets(locations={"Gage A": 0, "gage-a": 1}).to_flopy_obs()
    with pytest.raises(ValueError, match=r"GAGE_A from 'Gage A', 'gage-a'"):
        mf.SfrFlowTargets(locations={"Gage A": 0, "gage-a": 1}).to_flopy_obs()
    for targets, wrong in (
        (mf.LakeStageTargets(locations={"a": 0}), "VOLUME"),
        (mf.SfrStageTargets(locations={"a": 0}), "DOWNSTREAM-FLOW"),
        (mf.SfrFlowTargets(locations={"a": 0}), "SFR"),
    ):
        with pytest.raises(ValueError, match=rf"kind='{wrong}'"):
            targets.to_flopy_obs(kind=wrong)
    with pytest.raises(ValueError, match=r"kind='TO-MVR'"):
        mf.DrnFlowTargets(locations={"a": [0]}).to_flopy_obs(None, kind="TO-MVR")
    # The family's own quantity, in any case, still writes.
    assert mf.LakeStageTargets(locations={"a": 0}).to_flopy_obs(kind="stage") == {
        "lake_stage_targets.csv": [("a", "STAGE", 1)]}


# -- only the declared files are read -----------------------------------------------


@pytest.mark.slow
def test_a_stale_copy_beside_the_declared_file_is_not_read(declared, tmp_path):
    """A copy of the LAK file 20 ft low, ``archive_lake.csv``, sorts first: the name
    scan read it, and a PEST forward run on such a template wrote it every realization."""

    copy = tmp_path / "ws"
    shutil.copytree(declared.workspace, copy)
    live = pd.read_csv(copy / "lake_stage_targets.csv", float_precision="round_trip")
    live.assign(DRYING=live["DRYING"] - 20.0, HELD=live["HELD"] - 20.0).loc[
        :, ["time", "DRYING", "HELD"]].to_csv(copy / "archive_lake.csv", index=False)
    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    series = mf.LakeStageTargets(locations={"drying": 0, "held": 1}).simulated_series(reloaded)
    np.testing.assert_array_equal(series["held"], 118.0)
    np.testing.assert_array_equal(series["drying"], live["DRYING"])

    pd.DataFrame({"name": ["drying", "held"], "lake": [0, 1]}).to_csv(copy / "loc.csv", index=False)
    write_named_series_targets(
        [{"kind": "lake_stage", "locations_file": copy / "loc.csv",
          "output_csv": copy / "fw.csv"}],
        sim_ws=copy)
    written = pd.read_csv(copy / "fw.csv", float_precision="round_trip")
    np.testing.assert_array_equal(written["held"], 118.0)


@pytest.mark.slow
def test_a_declared_file_that_is_missing_raises(declared, tmp_path):
    """The declaration says where MF6 writes; a file that is not there is refused."""

    copy = tmp_path / "ws"
    shutil.copytree(declared.workspace, copy)
    (copy / "lake_stage_targets.csv").unlink()
    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    with pytest.raises(FileNotFoundError, match=(
            r"`drying STAGE 1` in lake_stage_targets.csv \(package 'lak'\) declares its "
            r"output in lake_stage_targets.csv, which is not in")):
        mf.LakeStageTargets(locations={"drying": 0}).simulated_series(reloaded)


@pytest.mark.slow
def test_a_declared_file_that_does_not_line_up_raises(declared, tmp_path):
    """Each record is read from its column BY POSITION, so a file whose columns are not
    the ones its OBS declares -- one column short, or two swapped -- is refused rather
    than read from the wrong column."""

    copy = tmp_path / "ws"
    shutil.copytree(declared.workspace, copy)
    frame = pd.read_csv(copy / "lake_stage_targets.csv")
    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    frame.drop(columns="POND").to_csv(copy / "lake_stage_targets.csv", index=False)
    with pytest.raises(ValueError, match=r"holds 3 observation column\(s\) where its OBS declares 4"):
        mf.LakeStageTargets(locations={"drying": 0}).simulated_series(reloaded)
    frame.loc[:, ["time", "HELD", "DRYING", "BN_HELD", "POND"]].to_csv(
        copy / "lake_stage_targets.csv", index=False)
    with pytest.raises(ValueError, match=r"column 1 of .* is headed 'HELD' where `drying STAGE 1`"):
        mf.LakeStageTargets(locations={"drying": 0}).simulated_series(reloaded)


@pytest.mark.slow
def test_a_reused_run_directory_reads_its_own_stage_file(tmp_path):
    """``Project.run(..., overwrite=True)`` reuses the directory, and MF6's CSV from the
    first run -- which declared the observation -- stays. The second run declares
    none: the name scan read the FIRST run's lake."""

    lake = dict(lakes=((DRYING, 111.0, "drying", 900.0),))
    first, project = _run(tmp_path, [_lake(**lake, observations={
        "lake_stage_targets.csv": [("drying", "STAGE", 1)]})])
    second, _ = _run(tmp_path, [_lake(lakes=((DRYING, 111.0, "drying", 300.0),))],
                     overwrite=True, project=project)
    stale = pd.read_csv(second.workspace / "lake_stage_targets.csv")
    raw = second.outputs.lak.stage.get(masked=False)[:, 0]
    # A dry record is the lake's bottom -- the 110 ft top of both its cells -- as MF6's
    # own STAGE observation would report it (ledger 351).
    own = np.where(raw == -1e30, 110.0, raw)
    assert (raw == -1e30).any()
    assert not np.allclose(stale["DRYING"].to_numpy()[:2], own[:2])

    series = mf.LakeStageTargets(locations={"drying": 0}).simulated_series(second)
    np.testing.assert_array_equal(series["time"], [0, 1, 2])
    np.testing.assert_array_equal(series["drying"], own)


@pytest.mark.slow
def test_an_sfr_stage_no_data_value_reads_nan(declared, tmp_path):
    """MF6's no-data value, 3e30, read as a stage of 0.0 ft on this route; a lake's
    read NaN since ledger 329. (MF6 6.8 writes an INACTIVE reach's stage as its bed
    top, measured, so the value is put in a copy of MF6's own file.)"""

    copy = tmp_path / "ws"
    shutil.copytree(declared.workspace, copy)
    lines = (copy / "sfr_stage_targets.csv").read_text().splitlines()
    assert lines[0].split(",")[3] == "DOWN"
    cells = lines[2].split(",")
    cells[3] = "3.0000000000000000E+30"
    lines[2] = ",".join(cells)
    (copy / "sfr_stage_targets.csv").write_text("\n".join(lines) + "\n")
    series = mf.SfrStageTargets(locations={"down": 3}).simulated_series(
        mf.load_mf6_run(copy, verbosity_level=0))
    assert np.isnan(series["down"].iloc[1])
    assert np.isfinite(series["down"].iloc[[0, 2]]).all()


# -- several packages of one type ----------------------------------------------------


@pytest.mark.slow
def test_several_packages_of_a_type_need_the_targets_package(two):
    """Each LAK declares ``pond STAGE 1`` -- its own lake 0 -- and each DRN ``s_c20``:
    the name scan read whichever file sorted first (``a.csv``, ``da.csv``)."""

    b, db = _csv(two, "b.csv"), _csv(two, "db.csv")
    assert not np.allclose(_csv(two, "a.csv")["POND"], b["POND"])
    with pytest.raises(ValueError, match=r"2 LAK packages \('lak_a', 'lak_b'\).*package=<its name>"):
        mf.LakeStageTargets(locations={"pond": 0}).simulated_series(two)
    series = mf.LakeStageTargets(locations={"pond": 0}, package="LAK_B").simulated_series(two)
    np.testing.assert_array_equal(series["pond"], b["POND"])
    with pytest.raises(ValueError, match=r"package='lak_c' names no LAK package.*'lak_a', 'lak_b'"):
        mf.LakeStageTargets(locations={"pond": 0}, package="lak_c").simulated_series(two)
    # No label declared: the LAK stage file belongs to ONE package, and it failed on
    # FloPy's list of two (AttributeError: 'list' object has no attribute 'nlakes').
    with pytest.raises(ValueError, match=r"2 LAK packages.*none of the targets is declared"):
        mf.LakeStageTargets(locations={"x": 0}).simulated_series(two)

    with pytest.raises(ValueError, match=r"2 DRN packages \('drn_a', 'drn_b'\)"):
        mf.DrnFlowTargets(locations={"s": [20]}).simulated_series(two)
    series = mf.DrnFlowTargets(locations={"s": [20]}, package="drn_b").simulated_series(two)
    np.testing.assert_array_equal(series["s"], db["S_C20"])
    # A zone no DRN declares reads the budget file, which sums EVERY package's drains
    # in its cells -- and, named, reads that package's alone (ledger 534:
    # it was refused, the budget table carrying no package to tell them apart by).
    both = mf.DrnFlowTargets(locations={"z": [20]}).simulated_series(two)
    np.testing.assert_allclose(
        both["z"].to_numpy(), (_csv(two, "da.csv")["S_C20"] + db["S_C20"]).to_numpy(), rtol=1e-10)
    alone = mf.DrnFlowTargets(locations={"z": [20]}, package="drn_a").simulated_series(two)
    np.testing.assert_allclose(alone["z"].to_numpy(), _csv(two, "da.csv")["S_C20"].to_numpy(),
                               rtol=1e-10)


@pytest.mark.slow
def test_the_targets_package_reaches_attach_and_the_pest_forward_run(two, tmp_path):
    """``package=`` names the package on every path: the OBS a target set attaches,
    the PEST build's forward-run call, and the forward run's own read."""

    targets = mf.LakeStageTargets(locations={"pond": 0}, package="lak_b")
    project = SimpleNamespace(model=two, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_lake_stage_observations(
        project, LakeStageObservationSpec(targets=targets, prefix="lk"))
    config = prepared["named_series_forward_run_config"]
    assert config["package"] == "lak_b"
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        [{"kind": "lake_stage", "locations_file": tmp_path / config["locations_file"],
          "output_csv": forward, "package": config["package"]}],
        sim_ws=two.workspace)
    assert forward.read_text() == (tmp_path / config["output_csv"]).read_text()
    np.testing.assert_array_equal(
        pd.read_csv(forward, float_precision="round_trip")["pond"], _csv(two, "b.csv")["POND"])

    fresh = mf.load_mf6_run(two.workspace, verbosity_level=0)
    attached = mf.LakeStageTargets(locations={"other": 0}, package="lak_b").attach_flopy_obs(
        fresh, filename="other.obs")
    assert attached.parent_file.package_name == "lak_b"


# -- the review of ledger 344 (ledger 345) -------------------------------------------


@pytest.fixture(scope="module")
def named_like_the_type(tmp_path_factory):
    """``two`` with the first LAK and DRN named EXACTLY their type (``lak``, ``drn``) --
    what a second ``mf.lak``/``mf.drn`` beside a default-named first one gives."""

    model, _ = _run(tmp_path_factory.mktemp("named_like_the_type"), [
        _lake(name="lak", lakes=((DRYING, 114.0, "pa", 100.0),),
              observations={"a.csv": [("pond", "STAGE", 1)]}),
        _lake(name="lak_2", lakes=((HELD, 116.0, "pb", 100.0),),
              observations={"b.csv": [("pond", "STAGE", 1)]}),
        _drain("drn", [20], {"da.csv": [("s_c20", "DRN", (0, 20))]}, elevation=96.0),
        _drain("drn_2", [20], {"db.csv": [("s_c20", "DRN", (0, 20))]}, elevation=97.0),
    ])
    return model


@pytest.mark.slow
def test_a_package_named_like_its_type_does_not_hide_the_others(named_like_the_type):
    """FloPy's ``get_package("lak")`` answers the package NAMED ``lak`` before the type,
    so a lookup by that spelling saw ONE LAK here: ``pond`` read ``lak``'s lake in
    silence, and a set with nothing declared read ``lak``'s stage file as the model's
    (ledger 344 found the packages by type; nothing tested it until ledger 345)."""

    model = named_like_the_type
    a, b = _csv(model, "a.csv"), _csv(model, "b.csv")
    assert not np.allclose(a["POND"], b["POND"])
    with pytest.raises(ValueError, match=r"2 LAK packages \('lak', 'lak_2'\).*package=<its name>"):
        mf.LakeStageTargets(locations={"pond": 0}).simulated_series(model)
    with pytest.raises(ValueError, match=r"2 LAK packages \('lak', 'lak_2'\); none of the targets"):
        mf.LakeStageTargets(locations={"x": 0}).simulated_series(model)
    series = mf.LakeStageTargets(locations={"pond": 0}, package="lak_2").simulated_series(model)
    np.testing.assert_array_equal(series["pond"], b["POND"])
    series = mf.LakeStageTargets(locations={"pond": 0}, package="LAK").simulated_series(model)
    np.testing.assert_array_equal(series["pond"], a["POND"])

    with pytest.raises(ValueError, match=r"2 DRN packages \('drn', 'drn_2'\)"):
        mf.DrnFlowTargets(locations={"s": [20]}).simulated_series(model)
    series = mf.DrnFlowTargets(locations={"s": [20]}, package="drn_2").simulated_series(model)
    np.testing.assert_array_equal(series["s"], _csv(model, "db.csv")["S_C20"])


@pytest.mark.slow
def test_a_binary_observation_file_is_read_as_mf6_wrote_it(declared, tmp_path):
    """``FILEOUT <file> BINARY`` is MF6's other output form: the reader opened it as
    text and died on ``UnicodeDecodeError``. MF6 runs the declaration here, and the
    binary file reads the CSV's values to the bit, by position and header."""

    import subprocess

    copy = tmp_path / "ws"
    shutil.copytree(declared.workspace, copy)
    text = (copy / "m.lak.obs").read_text()
    assert "FILEOUT  lake_stage_targets.csv\n" in text
    (copy / "m.lak.obs").write_text(
        text.replace("BEGIN continuous  FILEOUT  lake_stage_targets.csv",
                     "BEGIN continuous  FILEOUT  lake_stage_targets.bsv  BINARY")
        .replace("END continuous  FILEOUT  lake_stage_targets.csv", "END continuous"))
    (copy / "lake_stage_targets.csv").unlink()
    run = subprocess.run(["mf6"], cwd=copy, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout[-2000:]
    assert (copy / "lake_stage_targets.bsv").read_bytes()[:11] == b"cont double"

    csv = _csv(declared, "lake_stage_targets.csv")
    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    series = mf.LakeStageTargets(locations={"drying": 0, "held": 1}).simulated_series(reloaded)
    np.testing.assert_array_equal(series["time"], PERIODS)
    np.testing.assert_array_equal(series["drying"], csv["DRYING"])
    np.testing.assert_array_equal(series["held"], 118.0)
    with pytest.raises(ValueError, match=r"another quantity, not STAGE: `pond VOLUME 1` in lake_stage_targets.bsv"):
        mf.LakeStageTargets(locations={"pond": 0}).simulated_series(reloaded)

    whole = (copy / "lake_stage_targets.bsv").read_bytes()
    (copy / "lake_stage_targets.bsv").write_bytes(whole[:-8])
    with pytest.raises(ValueError, match=r"not a whole number of rows of 5 \(time and 4 records\)"):
        mf.LakeStageTargets(locations={"drying": 0}).simulated_series(reloaded)


@pytest.mark.slow
def test_the_pest_forward_run_reads_the_package_the_targets_name(two, tmp_path):
    """The forward run PEST executes must be handed ``package=``. Dropping it from the
    call the build writes into ``forward_run.py`` failed no test (ledger 345), and every
    realization of a model with two LAK (or DRN) packages would then refuse the set
    the base run read. The generated ``forward_run.py`` runs here, as PEST++ runs it."""

    import subprocess
    import sys

    pytest.importorskip("pyemu")
    lakes = mf.LakeStageTargets(locations={"pond": 0}, package="lak_b",
                                values={"time": [0], "pond": [112.0]})
    seeps = mf.DrnFlowTargets(locations={"s": [20]}, package="drn_b",
                              values=pd.DataFrame({"time": [0], "s": [-150.0]}))
    cal = two.pest("pkg", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(lakes)
    cal.observe(seeps)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("pkg.pst", noptmax=0)
    template = cal.template_workspace
    forward = (template / "forward_run.py").read_text()
    assert "'package': 'lak_b'" in forward and "'package': 'drn_b'" in forward
    tables = {name: (template / name).read_text() for name in (
        "lak_stage_simulated_lake_stage.csv", "drn_flow_simulated_drn_flow.csv")}

    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)
    for name, base in tables.items():
        assert (template / name).read_text() == base, name
    np.testing.assert_array_equal(
        pd.read_csv(template / "lak_stage_simulated_lake_stage.csv", float_precision="round_trip")["pond"],
        _csv(two, "b.csv")["POND"])
    np.testing.assert_array_equal(
        pd.read_csv(template / "drn_flow_simulated_drn_flow.csv", float_precision="round_trip")["s"],
        _csv(two, "db.csv")["S_C20"])
