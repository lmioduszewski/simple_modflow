"""A stage target reads ONE stage per feature per saved record (ledger 331).

``LakeStageTargets``' stage-file route read ``packages.lak.results.stage`` -- one
row per CONNECTED CELL, each a copy of the lake's one stage -- and
``_aggregate_named_series`` ADDED the copies: the user's excavated pond, 763
connected cells, read 401,339.97 ft for a stage of 526.00, and that was the series
PEST's base run and every forward run wrote. ``SfrStageTargets`` shared the helper,
and its table is one row per reach per saved TIME STEP, so a period saving three
steps read three stages added together. A stage is an elevation, and an elevation
is never added (ledger 297's rule, ``package_registry.agg_for_kind``).

Three small real MF6 runs:

``lake_runs``  a 6 x 6 irregular Voronoi grid, one confined layer, CHD at 50 ft on
               the west and east edges, six 10-day periods, one step each. Lake
               ``0`` has THREE vertical connections and a withdrawal that drains
               it: dry in records 1-3, refilled by rain in record 4. Lake ``1`` has
               TWO and is held CONSTANT at 108 ft -- the ledger-329 case that read
               216.0. LAK STAGE observations are attached for the names
               ``drying``/``held``, so those names take the observation route and
               any other name the stage-file route.
``sfr_run``    the same grid with a four-reach stream, THREE time steps per period
               saved (``SAVE HEAD ALL``), heads filling from 90 ft toward 102 so the
               stage moves within each period; SFR STAGE observations for
               ``up``/``down``.
``lake_steps_run``  the same grid with ONE lake of three vertical connections and a
               steady withdrawal, three periods of THREE time steps each saved
               (``SAVE HEAD ALL``); no LAK observations, so every name takes the
               stage-file route. It pins a limit, not a fix: that route is keyed by
               saved record, not by period (ledger 331, open under 301).

What each file holds is asserted before it is used, not assumed.
"""

from __future__ import annotations

import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations._shared import _aggregate_named_series
from myflopy.modflow.mf6.package_registry import FieldSpec, agg_for_kind
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
from myflopy.modflow.mf6.pest.observations import (
    prepare_lake_stage_observations,
    prepare_sfr_stage_observations,
)
from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec, SfrStageObservationSpec

NPER = 6
DRY_RECORDS = (1, 2, 3)  # lake 0
SFR_NPER, SFR_NSTP = 3, 3
LAKE_STEPS_NPER, LAKE_STEPS_NSTP = 3, 3


def _grid():
    """The 6 x 6 grid every run uses, its DISV props, and its west/east edge cells."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    return grid, props, edge


def _disv(grid, props, *, top: float):
    """One confined layer from ``top`` down to 0."""

    return mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                   vertices=props["vertices"], cell2d=props["cell2d"], top=top, botm=[0.0])


def _run(root, simulation):
    """Run ``simulation`` in a fresh project and return the run model."""

    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


@pytest.fixture(scope="module")
def lake_runs(tmp_path_factory):
    """Two lakes of three and two connected cells; see the module docstring."""

    grid, props, edge = _grid()
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    drying, held = interior[:3], interior[-2:]
    connections = [[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                   for i, c in enumerate(drying)]
    connections += [[1, i, (0, c), "VERTICAL", 0.01, 0.0, 0.0, 0.0, 0.0]
                    for i, c in enumerate(held)]
    keep = [[1, "STATUS", "CONSTANT"], [1, "STAGE", 108.0]]
    observed = mf.LakeStageTargets(locations={"drying": 0, "held": 1})
    lake = mf.lak.flopy(
        nlakes=2,
        packagedata=[[0, 100.5, len(drying), "drying"], [1, 108.0, len(held), "held"]],
        connectiondata=connections,
        perioddata={0: [[0, "WITHDRAWAL", 900.0], *keep],
                    4: [[0, "RAINFALL", 0.5], *keep],
                    5: [[0, "RAINFALL", 0.0], *keep]},
        boundnames=True,
        stage_filerecord="m.lak.stg",
        budget_filerecord="m.lak.bud",
        observations=observed.to_flopy_obs(csv_name="lake_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props, top=100.0),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=50.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 50.0] for c in edge]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 1, 1.0)] * NPER)])
    return _run(tmp_path_factory.mktemp("lake_one_stage"), simulation)


@pytest.fixture(scope="module")
def lake_steps_run(tmp_path_factory):
    """One lake of three connected cells saving three steps per period; see the module docstring."""

    grid, props, edge = _grid()
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    connections = [[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                   for i, c in enumerate(interior[:3])]
    lake = mf.lak.flopy(
        nlakes=1,
        packagedata=[[0, 104.0, len(connections), "pond"]],
        connectiondata=connections,
        perioddata={0: [[0, "WITHDRAWAL", 100.0]]},
        boundnames=True,
        stage_filerecord="m.lak.stg",
        budget_filerecord="m.lak.bud",
    )
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props, top=100.0),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=50.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 50.0] for c in edge]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=LAKE_STEPS_NPER, perioddata=[(10.0, LAKE_STEPS_NSTP, 1.0)] * LAKE_STEPS_NPER)])
    return _run(tmp_path_factory.mktemp("lake_steps"), simulation)


@pytest.fixture(scope="module")
def sfr_run(tmp_path_factory):
    """A four-reach stream saving three time steps per period; see the module docstring."""

    grid, props, edge = _grid()
    row = [13, 14, 15, 16]  # one interior row, west to east
    last = len(row) - 1
    packagedata = [
        [i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
         1 if i in (0, last) else 2, 1.0, 0]
        for i, c in enumerate(row)
    ]
    connectiondata = [
        [i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
        for i in range(len(row))
    ]
    observed = mf.SfrStageTargets(locations={"up": 0, "down": last})
    stream = mf.sfr.flopy(
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                    2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg",
        budget_filerecord="m.sfr.bud",
        observations=observed.to_flopy_obs(csv_name="sfr_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props, top=110.0),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=SFR_NPER, perioddata=[(10.0, SFR_NSTP, 1.0)] * SFR_NPER)])
    return _run(tmp_path_factory.mktemp("sfr_one_stage"), simulation)


def _lake_stage_file(model):
    """The masked LAK stage file, asserting the fixture is what the docstring says."""

    raw = model.outputs.lak.stage.get(masked=False)
    assert raw.shape == (NPER, 2)
    assert np.array_equal(np.flatnonzero(raw[:, 0] == -1e30), DRY_RECORDS)
    np.testing.assert_array_equal(raw[:, 1], 108.0)
    rows = model.lak.results.stage.get().groupby(["per", "lake"]).size()
    assert set(rows.xs(0, level="lake")) == {3} and set(rows.xs(1, level="lake")) == {2}
    return model.outputs.lak.stage.get()


def _sfr_stage_file(model):
    """``{(per, reach): [stage at kstp 0, 1, 2]}``, asserting three steps saved per period."""

    table = model.sfr.results.stage.get()
    assert set(table.groupby(["per", "reach"]).size()) == {SFR_NSTP}
    table = table.assign(kstp=table["kstpkper"].map(lambda value: int(value[0])))
    return {key: group.sort_values("kstp")["stage"].to_numpy()
            for key, group in table.groupby(["per", "reach"])}


# -- the lake ------------------------------------------------------------------


@pytest.mark.slow
def test_a_lake_reads_its_stage_once_not_once_per_connected_cell(lake_runs):
    """The stage-file route added a lake's copies: 3x lake 0's stage, and 216.0 for a
    lake held at 108.0 with two connections."""

    stage = _lake_stage_file(lake_runs)
    series = mf.LakeStageTargets(locations={"a": 0, "b": 1}).simulated_series(lake_runs)

    assert list(series["time"]) == list(range(NPER))
    np.testing.assert_array_equal(series["a"], stage[:, 0])
    np.testing.assert_array_equal(series["b"], 108.0)


@pytest.mark.slow
def test_the_two_lake_routes_agree_wherever_the_lake_holds_water(lake_runs):
    """The observation route always read the lake itself; the stage-file route now
    reads the same number. They still differ where MF6's two outputs do -- a DRY
    lake is its bottom in the observation, ``NaN`` in the file (ledger 329) -- and
    in how ``time`` is labelled (MF6's time against the record index)."""

    stage = _lake_stage_file(lake_runs)
    observed = mf.LakeStageTargets(locations={"drying": 0, "held": 1}).simulated_series(lake_runs)
    # Any name the observation CSV does not hold takes the stage-file route; the
    # route is chosen by NAME, so a trailing "_" (dropped from a label) would not do.
    from_file = mf.LakeStageTargets(locations={"a": 0, "b": 1}).simulated_series(lake_runs)

    np.testing.assert_array_equal(observed["time"], 10.0 * np.arange(1, NPER + 1))
    np.testing.assert_array_equal(from_file["time"], np.arange(NPER))
    wet = np.isfinite(stage[:, 0])
    np.testing.assert_allclose(observed["drying"][wet], from_file["a"][wet], rtol=1e-14)
    np.testing.assert_array_equal(observed["drying"][~wet], 100.0)
    assert from_file["a"][~wet].isna().all()
    np.testing.assert_array_equal(observed["held"], from_file["b"])


@pytest.mark.slow
def test_a_lake_compare_reads_the_same_one_stage(lake_runs):
    """``compare()``'s stage-file fallback went through the same sum."""

    stage = _lake_stage_file(lake_runs)
    targets = mf.LakeStageTargets(
        locations={"a": 0, "b": 1},
        values=pd.DataFrame({"time": range(NPER), "a": 100.5, "b": 108.0}))
    compare = targets.compare(lake_runs)

    by_lake = {name: frame.sort_values("per") for name, frame in compare.groupby("name")}
    np.testing.assert_array_equal(by_lake["a"]["sim_stage"], stage[:, 0])
    np.testing.assert_array_equal(by_lake["b"]["sim_stage"], 108.0)
    np.testing.assert_array_equal(by_lake["b"]["residual"], 0.0)


@pytest.mark.parametrize("names", [("drying", "held"), ("a", "b")],
                         ids=["observation-route", "stage-file-route"])
@pytest.mark.slow
def test_the_files_pest_reads_carry_each_lake_once(lake_runs, tmp_path, names):
    """What PEST calibrated against: the base run's CSV and every forward run's."""

    stage = _lake_stage_file(lake_runs)
    targets = mf.LakeStageTargets(
        locations=dict(zip(names, (0, 1), strict=True)),
        values=pd.DataFrame({"time": range(NPER), names[0]: 100.5, names[1]: 108.0}))
    project = SimpleNamespace(model=lake_runs, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_lake_stage_observations(
        project, LakeStageObservationSpec(targets=targets, prefix="lk"))
    config = prepared["named_series_forward_run_config"]
    base = tmp_path / config["output_csv"]
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        "lake_stage", tmp_path / config["locations_file"], forward, sim_ws=lake_runs.workspace)

    assert forward.read_text() == base.read_text()
    written = pd.read_csv(base)
    wet = np.isfinite(stage[:, 0])
    np.testing.assert_allclose(written[names[0]][wet], stage[wet, 0], rtol=1e-14)
    np.testing.assert_array_equal(written[names[1]], 108.0)
    # A dry lake: the bottom MF6 observes on the observation route, PEST's no-value
    # fill (ledger 329) on the stage-file route.
    expected_dry = 100.0 if names[0] == "drying" else 1e30
    np.testing.assert_array_equal(written[names[0]][~wet], expected_dry)


@pytest.mark.slow
def test_a_lake_target_keyed_by_period_reads_the_record_clock_at_several_saves(
    lake_steps_run,
):
    """A LIMIT, pinned so the docs cannot outlive it (ledger 331, open under 301).

    The LAK stage table labels ``per`` by SAVED-RECORD index, so at three saves per
    period the stage-file route has nine rows labelled 0..8, and a target keyed by
    stress period 0, 1, 2 reads records 0, 1, 2 -- all inside period 0 -- where
    ``HeadTargets`` on the same run reads each period's LAST step. Nothing warns.
    When ``build_lak_stage_result_table`` carries the stage file's own ``kstpkper``
    this fails: update ledger 331's time-labels bullet, the ``LakeStageTargets``
    docstrings and ``docs/package_api_reference.md`` with it.
    """

    records = LAKE_STEPS_NPER * LAKE_STEPS_NSTP
    stage = lake_steps_run.outputs.lak.stage.get()
    assert stage.shape == (records, 1)
    # The stage moves at every step, so which step is read is visible.
    assert len(np.unique(stage[:, 0])) == records
    rows = lake_steps_run.lak.results.stage.get().groupby(["per", "lake"]).size()
    assert set(rows) == {3}

    series = mf.LakeStageTargets(locations={"a": 0}).simulated_series(lake_steps_run)
    assert list(series["time"]) == list(range(records))
    np.testing.assert_array_equal(series["a"], stage[:, 0])

    targets = mf.LakeStageTargets(
        locations={"a": 0},
        values=pd.DataFrame({"time": range(LAKE_STEPS_NPER), "a": 100.0}))
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        compare = targets.compare(lake_steps_run).sort_values("per")
    in_period_0 = stage[:LAKE_STEPS_NPER, 0]
    period_end = stage[LAKE_STEPS_NSTP - 1::LAKE_STEPS_NSTP, 0]
    np.testing.assert_array_equal(compare["sim_stage"], in_period_0)
    assert not np.isclose(compare["sim_stage"].to_numpy()[1:], period_end[1:]).any()

    # The head target beside it, keyed the same way, reads period end.
    cell = int(lake_steps_run.lak.results.stage.get()["cell"].iloc[0])
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["h"], "cell": [cell], "layer": [0]}),
        values=pd.DataFrame({"time": range(LAKE_STEPS_NPER), "name": "h", "head": 50.0}))
    simulated = heads.simulated_heads(lake_steps_run)
    table = lake_steps_run.all_heads.reset_index()
    table = table[(table["layer"] == 0) & (table["cell"] == cell)]
    head_at = {(int(kper), int(kstp)): float(elev)
               for (kstp, kper), elev in zip(table["kstpkper"], table["elev"], strict=True)}
    assert len(head_at) == records
    period_end_head = [head_at[(per, LAKE_STEPS_NSTP - 1)] for per in range(LAKE_STEPS_NPER)]
    assert list(simulated["per"]) == list(range(LAKE_STEPS_NPER))
    np.testing.assert_array_equal(simulated["h"], period_end_head)


# -- the stream ----------------------------------------------------------------


@pytest.mark.slow
def test_a_reach_reads_its_period_end_stage_not_the_sum_of_every_step(sfr_run):
    """Three saved steps per period were added: a reach read about three times its
    stage. It reads the LAST step now, the period-end rule ``HeadTargets`` uses."""

    stage = _sfr_stage_file(sfr_run)
    series = mf.SfrStageTargets(locations={"a": 0, "b": 3}).simulated_series(sfr_run)

    assert list(series["time"]) == list(range(SFR_NPER))
    for column, reach in (("a", 0), ("b", 3)):
        expected = [stage[(per, reach)][-1] for per in range(SFR_NPER)]
        np.testing.assert_array_equal(series[column], expected)
    # The stage moves within a period, so the choice of step is exercised.
    assert any(not np.isclose(stage[key][0], stage[key][-1]) for key in stage)


@pytest.mark.slow
def test_the_two_sfr_routes_agree_at_each_period_end(sfr_run):
    """The observation route writes every time step at MF6's time; its period-end
    rows are the stage-file route's rows."""

    observed = mf.SfrStageTargets(locations={"up": 0, "down": 3}).simulated_series(sfr_run)
    from_file = mf.SfrStageTargets(locations={"a": 0, "b": 3}).simulated_series(sfr_run)

    assert len(observed) == SFR_NPER * SFR_NSTP
    period_end = observed.iloc[SFR_NSTP - 1::SFR_NSTP].reset_index(drop=True)
    np.testing.assert_allclose(period_end["time"], 10.0 * np.arange(1, SFR_NPER + 1))
    assert len(from_file) == SFR_NPER
    np.testing.assert_allclose(period_end["up"], from_file["a"], rtol=1e-14)
    np.testing.assert_allclose(period_end["down"], from_file["b"], rtol=1e-14)


@pytest.mark.slow
def test_an_sfr_compare_and_the_files_pest_reads_carry_the_period_end_stage(sfr_run, tmp_path):
    """``compare()``'s fallback, the base run's CSV and a forward run's, together."""

    stage = _sfr_stage_file(sfr_run)
    expected = [stage[(per, 0)][-1] for per in range(SFR_NPER)]
    targets = mf.SfrStageTargets(
        locations={"a": 0},
        values=pd.DataFrame({"time": range(SFR_NPER), "a": 100.0}))
    compare = targets.compare(sfr_run).sort_values("per")
    np.testing.assert_array_equal(compare["sim_stage"], expected)

    project = SimpleNamespace(model=sfr_run, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_sfr_stage_observations(
        project, SfrStageObservationSpec(targets=targets, prefix="st"))
    config = prepared["named_series_forward_run_config"]
    base = tmp_path / config["output_csv"]
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        "sfr_stage", tmp_path / config["locations_file"], forward, sim_ws=sfr_run.workspace)

    assert forward.read_text() == base.read_text()
    np.testing.assert_allclose(pd.read_csv(base)["a"], expected, rtol=1e-14)


# -- the reduction itself --------------------------------------------------------


def _stub_lake_model(table: pd.DataFrame):
    """A model whose only output is this LAK stage table: the stage-file route."""

    stage = SimpleNamespace(get=lambda: table.copy())
    return SimpleNamespace(lak=SimpleNamespace(results=SimpleNamespace(stage=stage)))


def test_copies_that_disagree_are_refused_naming_the_lake_the_record_and_the_spread():
    """MF6's table cannot produce this; averaging it would invent a stage."""

    table = pd.DataFrame({"per": [0, 0, 1, 1], "lake": [0, 0, 0, 0], "cell": [4, 9, 4, 9],
                          "stage": [101.0, 101.0, 101.0, 102.5]})
    targets = mf.LakeStageTargets(locations={"pond": 0})
    with pytest.raises(ValueError, match=r"lake 0 at per=1: 2 rows span 101\.\.102\.5 \(spread 1\.5\)"):
        targets.simulated_series(_stub_lake_model(table))
    with pytest.raises(ValueError, match=r"lake 0 at per=1"):
        targets.compare(_stub_lake_model(table))


def test_an_sfr_stage_is_an_elevation_in_the_series_and_in_compare():
    """MF6's SFR table has one row per reach per record, so a sum and the one value
    agree there and only a table with copies tells the kinds apart: a reach's two
    rows in one record are refused, not added."""

    table = pd.DataFrame({"kstpkper": [(0, 0), (0, 0)], "per": [0, 0], "reach": [2, 2],
                          "stage": [95.0, 96.0]})
    stage = SimpleNamespace(get=lambda: table.copy())
    model = SimpleNamespace(sfr=SimpleNamespace(results=SimpleNamespace(stage=stage)))
    targets = mf.SfrStageTargets(locations={"gage": 2})
    with pytest.raises(ValueError, match=r"reach 2 at per=0: 2 rows span 95\.\.96"):
        targets.simulated_series(model)
    with pytest.raises(ValueError, match=r"reach 2 at per=0"):
        targets.compare(model)


def test_a_record_where_only_some_copies_hold_a_stage_is_refused():
    table = pd.DataFrame({"per": [0, 0], "lake": [0, 0], "cell": [4, 9],
                          "stage": [101.0, np.nan]})
    with pytest.raises(ValueError, match=r"lake 0 at per=0: 1 of 2 rows hold a stage"):
        mf.LakeStageTargets(locations={"pond": 0}).simulated_series(_stub_lake_model(table))


def test_copies_within_rounding_and_a_record_with_no_stage_pass():
    table = pd.DataFrame({"per": [0, 0, 1, 1], "lake": [0, 0, 0, 0], "cell": [4, 9, 4, 9],
                          "stage": [526.0025, 526.0025 + 1e-10, np.nan, np.nan]})
    series = mf.LakeStageTargets(locations={"pond": 0}).simulated_series(_stub_lake_model(table))
    assert series["pond"].iloc[0] == 526.0025 and np.isnan(series["pond"].iloc[1])


def test_only_an_extensive_value_adds_and_the_rule_is_the_registrys():
    """One rule, ``agg_for_kind``, serves ``FieldSpec.agg`` and the targets."""

    table = pd.DataFrame({"per": [0, 0, 0, 1, 1], "zone": [0, 0, 1, 0, 0],
                          "q": [1.5, 2.0, np.nan, np.nan, np.nan]})
    summed = _aggregate_named_series(table, id_column="zone", value_column="q", kind="extensive")
    assert summed["q"].iloc[0] == 3.5
    assert np.isnan(summed.loc[summed["zone"] == 1, "q"]).all()  # all-NaN is not 0.0
    assert np.isnan(summed["q"].iloc[-1])
    for kind in ("extensive", "elevation", "intensive"):
        assert FieldSpec("x", kind=kind).agg == agg_for_kind(kind)
    assert agg_for_kind("extensive") == "sum" and agg_for_kind("elevation") == "mean"
    with pytest.raises(ValueError, match="kind 'stage' is not one of"):
        _aggregate_named_series(table, id_column="zone", value_column="q", kind="stage")


def test_a_table_saving_several_steps_keeps_each_periods_last():
    """The SFR stage table's shape: one row per reach per saved (kstp, kper)."""

    table = pd.DataFrame({
        "kstpkper": [(0, 0), (1, 0), (2, 0), (0, 1)],
        "per": [0, 0, 0, 1], "reach": [0, 0, 0, 0], "stage": [10.0, 11.0, 12.0, 20.0]})
    reduced = _aggregate_named_series(table, id_column="reach", value_column="stage",
                                      kind="elevation")
    assert reduced["stage"].tolist() == [12.0, 20.0]
