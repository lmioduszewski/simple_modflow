"""A lake's stage file is read on its own clock, and a period is its last saved step (ledger 348).

``build_lak_stage_result_table`` labelled ``per`` by the INDEX of the saved stage
record (ledger 301's stage clock), so at three saves per period the nine records of
a three-period run read as periods 0..8: a ``LakeStageTargets`` keyed by period
read steps 1-3 of period 0 and a PEST build wrote nine rows for three periods; the
stage map, series, stage change, exchange hover and the group tier sat on the same
wrong clock. It carries the stage file's own ``kstpkper`` now and ``per`` is the
stress period. The views read period ``N`` at its LAST saved step (ledger 337's
rule) through ``period_end.at_last_saved_step(s)``, and so do the pictures that had
been ADDING a period's saved budget steps (the LAK/SFR/cell-budget maps and every
series) and the tables that did the same (the SFR long profile, the LAK exchange
budget, the combined surface-water table, the group diff maps).

Five small real MF6 runs on one 6 x 6 irregular Voronoi grid (one confined layer,
CHD at 50 ft on the west and east edges), ten-day periods:

``steps``    one lake of three vertical connections and a 100 ft3/d withdrawal,
             three periods of THREE steps, ``SAVE HEAD ALL`` + ``SAVE BUDGET ALL``.
``uneven``   the same lake, periods of 1, 3 and 2 steps, everything saved -- the
             record index differs from the period at periods 1 and 2.
``asym``     the same lake at ``nstp=3``, ``SAVE HEAD LAST`` + ``SAVE BUDGET ALL``:
             three stage records against nine budget records.
``stream``   a four-reach stream, three periods of three steps, everything saved
             (the SFR stage and exchange views share the fixes).
``ats``      the lake of ``steps`` under ATS, one declared step per period and
             everything saved: MF6 takes four steps a period, decided at run time.

What each file holds is read from the file and asserted before it is used.
"""

from __future__ import annotations

import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.package_budget import (
    build_lak_stage_result_table,
    build_sfr_stage_result_table,
    build_surface_water_exchange_cell_table,
)
from myflopy.modflow.mf6.period_end import saved_steps
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
from myflopy.modflow.mf6.pest.observations import prepare_lake_stage_observations
from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec
from myflopy.modflow.utils.outputs import LakStageOutput

NPER = 3
PERLEN = 10.0
SAVE_ALL = (("HEAD", "ALL"), ("BUDGET", "ALL"))
REACH_CELLS = [13, 14, 15, 16]


def _grid():
    """The 6 x 6 grid every run uses, its DISV props, and its west/east edge cells."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    return grid, props, edge


def _disv(grid, props, *, top: float):
    return mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                   vertices=props["vertices"], cell2d=props["cell2d"], top=top, botm=[0.0])


def _run(root, simulation):
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and these
        # runs are never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _lake_run(root, *, nstp, saverecord, ats=False):
    grid, props, edge = _grid()
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    lake = mf.lak.flopy(
        nlakes=1,
        packagedata=[[0, 104.0, 3, "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(interior[:3])],
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
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc", saverecord=list(saverecord)),
    ])
    # ATS: MF6 picks the steps at run time (1, 2, 4 and 3 days here), so a period's
    # last saved step is found by time, never by its declared step count.
    timing = {"ats": True, "ats_dt0": 1.0, "ats_dtmin": 0.01, "ats_dtmax": 4.0} if ats else {}
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=len(nstp), perioddata=[(PERLEN, int(n), 1.0) for n in nstp], **timing)])
    return _run(root, simulation)


def _stream_run(root):
    grid, props, edge = _grid()
    last = len(REACH_CELLS) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACH_CELLS)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACH_CELLS))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]], 2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg",
        budget_filerecord="m.sfr.bud",
    )
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props, top=110.0),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc", saverecord=list(SAVE_ALL)),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, 3, 1.0)] * NPER)])
    return _run(root, simulation)


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """The five runs of the module docstring, by name."""

    root = tmp_path_factory.mktemp("lak_stage_clock")
    return {
        "steps": _lake_run(root / "steps", nstp=(3, 3, 3), saverecord=SAVE_ALL),
        "uneven": _lake_run(root / "uneven", nstp=(1, 3, 2), saverecord=SAVE_ALL),
        "asym": _lake_run(root / "asym", nstp=(3, 3, 3),
                          saverecord=(("HEAD", "LAST"), ("BUDGET", "ALL"))),
        "stream": _stream_run(root / "stream"),
        "ats": _lake_run(root / "ats", nstp=(1, 1, 1), saverecord=SAVE_ALL, ats=True),
    }


# -- what the files hold, read from the files -----------------------------------


def _stage_file(model):
    """``[(kstp, kper)]`` and the stage per record, from the LAK stage file itself."""

    reader = model.fl.lak.output.stage()
    clock = [(kstp, kper) for kstp, kper, _ in saved_steps(reader)]
    values = np.asarray(reader.get_alldata(), dtype=float).reshape(len(clock), -1)[:, 0]
    return clock, values


def _stage_at(model):
    """``{(kstp, kper): stage}`` of the one lake, asserting it moves at every step."""

    clock, values = _stage_file(model)
    assert len(np.unique(values)) == len(values), "the fixture cannot tell steps apart"
    return dict(zip(clock, values, strict=True))


def _last_steps(clock):
    """Each period's last saved ``(kstp, kper)``, worked out by hand."""

    last = {}
    for kstp, kper in clock:
        last[kper] = max(last.get(kper, -1), kstp)
    return {kper: (kstp, kper) for kper, kstp in sorted(last.items())}


def _lak_q(model, kstpkper):
    """``{cell: (q, flow_area)}`` of the LAK-GWF exchange at one saved step."""

    table = model.lak.results.q.get()
    rows = table[table["kstpkper"].map(lambda value: tuple(map(int, value))) == kstpkper]
    return {int(c): (float(q), float(a))
            for c, q, a in zip(rows["cell"], rows["q_lake"], rows["flow_area"], strict=True)}


def _z(picture):
    return np.asarray(picture.fig.data[0].z, dtype=float)


# -- (a) the table ---------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("name", ["steps", "uneven", "asym"])
def test_the_lak_stage_table_is_labelled_by_the_stage_files_own_clock(runs, name):
    """``kstpkper`` from the stage file's record headers, ``per`` its stress period;
    it was the saved record's INDEX. The asymmetric run is asserted asymmetric first:
    its LAK budget holds nine steps where the stage file holds three, so a table
    labelled from the budget clock (or ``per`` from a record index) cannot pass."""

    model = runs[name]
    clock, values = _stage_file(model)
    if name == "asym":
        assert len(saved_steps(model.fl.lak.output.budget())) == 9
        assert clock == [(2, 0), (2, 1), (2, 2)]
    table = model.lak.results.stage.get()
    assert list(table.columns) == [
        "model", "package", "kstpkper", "per", "lake", "layer", "cell", "stage"]
    assert list(table.columns[:4]) == list(build_sfr_stage_result_table(runs["stream"]).columns[:4])
    one = table.drop_duplicates(subset=["kstpkper", "lake"])
    assert list(one["kstpkper"]) == clock == model.outputs.lak.stage.kstpkper()
    assert list(one["per"]) == [kper for _, kper in clock]
    np.testing.assert_array_equal(one["stage"], values)
    assert set(table.groupby("kstpkper").size()) == {3}  # three connected cells
    if name != "asym":
        assert list(one["per"]) != list(range(len(clock)))  # not the record index


@pytest.mark.slow
def test_a_stage_file_whose_headers_do_not_line_up_is_refused(runs, monkeypatch):
    """A reader whose records do not line up RAISES (ledger 301's rule) rather
    than zipping the records against a shorter clock -- on the table and on the
    group's raw stage listing, which paired the records with the HEAD file's clock
    and cut them to its length."""

    model = runs["steps"]
    real = LakStageOutput.kstpkper
    monkeypatch.setattr(LakStageOutput, "kstpkper", lambda self: real(self)[:-1])
    refusal = r"holds 9 record\(s\) of 1 lake\(s\) but names 8 saved"
    with pytest.raises(ValueError, match=refusal):
        build_lak_stage_result_table(model)
    with pytest.raises(ValueError, match=refusal):
        mf.ModelGroup({"a": model}, reference="a").outputs.lak.stage()


# -- (b) the targets -------------------------------------------------------------


@pytest.mark.slow
def test_uneven_steps_read_each_periods_end_not_the_first_records(runs):
    """At 1, 3 and 2 steps the period ends are records 0, 3 and 5, where the record
    index read records 0, 1, 2 -- periods 1 and 2 at the wrong step."""

    model = runs["uneven"]
    at = _stage_at(model)
    clock, values = _stage_file(model)
    ends = _last_steps(clock)
    assert [clock.index(ends[per]) for per in range(NPER)] == [0, 3, 5]
    expected = [at[ends[per]] for per in range(NPER)]
    assert not np.isclose(values[:NPER], expected)[1:].any()

    series = mf.LakeStageTargets(locations={"pond": 0}).simulated_series(model)
    assert list(series["time"]) == list(range(NPER))
    np.testing.assert_array_equal(series["pond"], expected)
    targets = mf.LakeStageTargets(
        locations={"pond": 0}, values=pd.DataFrame({"time": range(NPER), "pond": 100.0}))
    np.testing.assert_array_equal(targets.compare(model).sort_values("per")["sim_stage"], expected)


@pytest.mark.slow
def test_a_pest_build_holds_one_lake_row_per_period_carrying_the_targets(runs, tmp_path):
    """The build wrote one observation per saved RECORD -- nine for three periods --
    with the targets on records 0-2 (steps inside period 0) and the base run's own
    stage on the rest. It writes one per period now, the targets on them, and every
    forward run writes the file the base run did."""

    model = runs["steps"]
    at = _stage_at(model)
    ends = _last_steps(list(at))
    observed = [111.0, 112.0, 113.0]
    targets = mf.LakeStageTargets(
        locations={"a": 0}, values=pd.DataFrame({"time": range(NPER), "a": observed}))
    cal = model.pest("lk", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obs = cal.build("lk.pst", noptmax=0).observation_data
    lake = obs[obs.obsnme.str.contains("usecol:a")]
    assert len(lake) == NPER
    np.testing.assert_array_equal(np.sort(lake["obsval"].to_numpy()), observed)
    assert (lake["weight"] > 0).all()

    project = SimpleNamespace(model=model, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_lake_stage_observations(
        project, LakeStageObservationSpec(targets=mf.LakeStageTargets(locations={"a": 0}),
                                          prefix="lk"))
    config = prepared["named_series_forward_run_config"]
    base = tmp_path / config["output_csv"]
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        "lake_stage", tmp_path / config["locations_file"], forward, sim_ws=model.workspace)
    assert forward.read_text() == base.read_text()
    written = pd.read_csv(base)
    assert list(written["time"]) == list(range(NPER))
    np.testing.assert_allclose(written["a"], [at[ends[per]] for per in range(NPER)], rtol=1e-14)


# -- (c) the stage views -----------------------------------------------------------


@pytest.mark.slow
def test_the_lak_stage_views_draw_each_periods_last_saved_step(runs):
    """map(per=1) read record 1 (step 2 of period 0), plot() nine "periods", and
    stage_change eight transitions between steps. ``get(per=1)`` stays a FILTER."""

    model = runs["steps"]
    at = _stage_at(model)
    ends = _last_steps(list(at))
    stage = model.lak.results.stage
    cells = sorted(set(stage.get()["cell"]))

    period_one = stage.get(per=1)
    assert sorted(set(period_one["kstpkper"])) == [(0, 1), (1, 1), (2, 1)]
    z = _z(stage.map(per=1))
    np.testing.assert_array_equal(z[cells], at[ends[1]])

    (line,) = stage.plot().data
    assert list(line.x) == list(range(NPER))
    np.testing.assert_array_equal(line.y, [at[ends[per]] for per in range(NPER)])

    change = model.lak.results.stage_change.get()
    assert list(zip(change["per0"], change["per1"], strict=True)) == [(0, 1), (1, 2)]
    np.testing.assert_array_equal(change["stage0"], [at[ends[0]], at[ends[1]]])
    np.testing.assert_array_equal(change["stage1"], [at[ends[1]], at[ends[2]]])


@pytest.mark.slow
@pytest.mark.parametrize("name", ["steps", "asym"])
def test_the_lak_exchange_views_read_one_step_and_the_stage_on_its_own_clock(runs, name):
    """q.map(per=1) ADDED the period's three flows (and three copies of each area),
    and its hover joined record 1 of the stage file -- a step of period 0. Each file
    is read at ITS OWN last saved step of period 1: under ``SAVE HEAD LAST`` +
    ``SAVE BUDGET ALL`` the budget's is (2, 1) of nine and the stage file's (2, 1) of
    three. The LAK exchange budget and the combined surface-water table added the
    steps too."""

    model = runs[name]
    at = _stage_at(model)
    stage_end = at[_last_steps(list(at))[1]]
    budget_clock = [(k, p) for k, p, _ in saved_steps(model.fl.lak.output.budget())]
    assert len(budget_clock) == 9
    end = _last_steps(budget_clock)[1]
    flows = _lak_q(model, end)
    first = _lak_q(model, (0, 1))
    assert all(not np.isclose(flows[c][0], first[c][0]) for c in flows)

    picture = model.lak.results.q.map(per=1)
    z = _z(picture)
    hover = picture._custom_hover
    for cell, (q, area) in flows.items():
        assert z[cell] == pytest.approx(q / area, rel=1e-12)
        assert hover["q"][cell] == pytest.approx(q, rel=1e-12)
        assert hover["flow_area"][cell] == pytest.approx(area, rel=1e-12)
        assert hover["stage"][cell] == stage_end

    (line,) = model.lak.results.q.plot().data
    totals = [sum(q for q, _ in _lak_q(model, _last_steps(budget_clock)[p]).values())
              for p in range(NPER)]
    np.testing.assert_allclose(line.y, totals, rtol=1e-12)

    budget = model.lak.results.q.budget.get(per=1)
    assert budget["record_count"].sum() == len(flows)
    assert budget["q_lake"].sum() == pytest.approx(sum(q for q, _ in flows.values()), rel=1e-12)

    combined = build_surface_water_exchange_cell_table(model, per=1, include="lak")
    by_cell = combined.set_index("cell")["exchange_intensity"]
    for cell, (q, area) in flows.items():
        assert by_cell[cell] == pytest.approx(q / area, rel=1e-12)


@pytest.mark.slow
def test_a_cell_budget_noun_draws_one_step_per_period(runs):
    """Every cell-budget map and series added a period's saved steps: CHD's too."""

    model = runs["steps"]
    table = model.chd.results.q.get()
    last = table[table["kstpkper"].map(lambda value: tuple(map(int, value))) == (2, 1)]
    first = table[table["kstpkper"].map(lambda value: tuple(map(int, value))) == (0, 1)]
    assert not np.isclose(last["q_gwf"].sum(), first["q_gwf"].sum())
    z = _z(model.chd.results.q.map(per=1))
    np.testing.assert_allclose(z[last["cell"].to_numpy()], last["q_gwf"].to_numpy(), rtol=1e-12)

    (line,) = model.chd.results.q.plot().data
    ends = table[table["kstpkper"].map(lambda value: int(value[0])) == 2]
    np.testing.assert_allclose(line.y, ends.groupby("per")["q_gwf"].sum().to_numpy(), rtol=1e-12)


@pytest.mark.slow
def test_the_sfr_stage_and_exchange_views_read_one_step(runs):
    """The SFR noun was on a period clock already, and still combined the steps:
    map(per=1) painted the FIRST step (``agg="first"``), plot() averaged them, and
    the exchange map and the long profile ADDED them."""

    model = runs["stream"]
    reader = model.fl.sfr.output.stage()
    clock = [tuple(int(v) for v in key) for key in reader.get_kstpkper()]
    data = np.asarray(reader.get_alldata(), dtype=float).reshape(len(clock), -1)
    stage_at = {key: data[i] for i, key in enumerate(clock)}
    assert not np.isclose(stage_at[(0, 1)][0], stage_at[(2, 1)][0])

    z = _z(model.sfr.results.stage.map(per=1))
    np.testing.assert_array_equal(z[REACH_CELLS], stage_at[(2, 1)])
    lines = {trace.name: np.asarray(trace.y, dtype=float) for trace in
             model.sfr.results.stage.plot().data}
    np.testing.assert_array_equal(lines["Reach 0"], [stage_at[(2, per)][0] for per in range(NPER)])

    q = model.sfr.results.q.get()
    at_end = q[q["kstpkper"].map(lambda value: tuple(map(int, value))) == (2, 1)]
    q_end = dict(zip(at_end["reach"].astype(int), at_end["q_gwf"], strict=True))
    profile = model.sfr.results.profile.get(per=1).set_index("reach")
    np.testing.assert_allclose(profile["q_gwf"].to_numpy(),
                               [q_end[r] for r in profile.index], rtol=1e-12)
    np.testing.assert_array_equal(profile["stage"].to_numpy(), stage_at[(2, 1)][profile.index])
    z = _z(model.sfr.results.q.map(per=1))
    rlen = dict(zip(at_end["reach"].astype(int), at_end["rlen"], strict=True))
    np.testing.assert_allclose(z[REACH_CELLS], [q_end[r] / rlen[r] for r in range(4)], rtol=1e-12)


# -- the review of ledger 348: the tables and profiles the choke points missed ------


@pytest.mark.slow
def test_the_sfr_field_profiles_draw_one_step_per_reach(runs):
    """``sfr.results.stage.profile`` and ``sfr.results.q.profile`` held every step
    period N saved -- three values per reach at ``nstp=3`` -- and their line
    zig-zagged through them, beside a long profile reading the last step."""

    model = runs["stream"]
    reader = model.fl.sfr.output.stage()
    clock = [tuple(int(v) for v in key) for key in reader.get_kstpkper()]
    data = np.asarray(reader.get_alldata(), dtype=float).reshape(len(clock), -1)
    stage_end = data[clock.index((2, 1))]
    q = model.sfr.results.q.get()
    at_end = q[q["kstpkper"].map(lambda value: tuple(map(int, value))) == (2, 1)]
    q_end = dict(zip(at_end["reach"].astype(int), at_end["q_gwf"], strict=True))
    assert len(model.sfr.results.stage.get(per=1)) == 3 * len(REACH_CELLS)  # a filter

    stage = model.sfr.results.stage.profile.get(per=1)
    assert list(stage["reach"]) == [0, 1, 2, 3]
    np.testing.assert_array_equal(stage["stage"], stage_end)
    (line,) = model.sfr.results.stage.profile(per=1).plot().data
    np.testing.assert_array_equal(line.y, stage_end)

    exchange = model.sfr.results.q.profile.get(per=1)
    assert list(exchange["reach"]) == [0, 1, 2, 3]
    np.testing.assert_allclose(exchange["q_gwf"], [q_end[r] for r in range(4)], rtol=1e-12)


@pytest.mark.slow
def test_the_budget_tables_by_period_read_the_number_the_map_draws(runs):
    """``q.wide()`` and the LAK/SFR ``budget`` namespaces' ``summary``/``wide``
    pivot or group by PERIOD, and ADDED every step the period saved -- a CHD
    ``wide(per=1)`` read -2152.678 beside a ``map(per=1)`` of -452.712 once ledger
    348 moved the map. Grouping by ``kstpkper`` keeps every step."""

    model = runs["steps"]
    table = model.chd.results.q.get()
    last = table[table["kstpkper"].map(lambda value: tuple(map(int, value))) == (2, 1)]
    wide = model.chd.results.q.wide(per=1).set_index("cell")["per_1"]
    z = _z(model.chd.results.q.map(per=1))
    for cell, value in zip(last["cell"].astype(int), last["q_gwf"], strict=True):
        assert wide[cell] == pytest.approx(value, rel=1e-12)
        assert z[cell] == pytest.approx(value, rel=1e-12)
    every = model.chd.results.q.wide(per=1, index=("kstpkper", "layer", "cell"))
    assert len(every) == 3 * len(last)

    flows = _lak_q(model, (2, 1))
    q_end = sum(q for q, _ in flows.values())
    area_end = sum(area for _, area in flows.values())
    summary = model.lak.budget.gwf.summary().set_index("per")
    assert summary.loc[1, "q"] == pytest.approx(q_end, rel=1e-12)
    assert summary.loc[1, "flow_area"] == pytest.approx(area_end, rel=1e-12)
    assert summary.loc[1, "record_count"] == len(flows)
    assert summary.loc[1, "q_per_area"] == pytest.approx(q_end / area_end, rel=1e-12)
    budget_view = model.lak.results.q.budget.get(per=1)
    assert summary.loc[1, "q"] == pytest.approx(budget_view["q_lake"].sum(), rel=1e-12)
    wide = model.lak.budget.wide(term="GWF").set_index("per")
    assert wide.loc[1, "GWF"] == pytest.approx(q_end, rel=1e-12)
    per_step = model.lak.budget.gwf.summary(by=["kstpkper", "lake"])
    assert len(per_step) == 9

    stream = runs["stream"]
    rows = stream.sfr.budget.get(term="GWF")
    step = rows[rows["kstpkper"].map(lambda value: tuple(map(int, value))) == (2, 1)]
    expected = dict(zip(step["reach"].astype(int), step["q"], strict=True))
    sfr = stream.sfr.budget.gwf.summary(per=1).set_index("reach")
    for reach, value in expected.items():
        assert sfr.loc[reach, "q"] == pytest.approx(value, rel=1e-12)
    sfr_wide = stream.sfr.budget.wide(term="GWF", per=1).set_index("reach")
    for reach, value in expected.items():
        assert sfr_wide.loc[reach, "GWF"] == pytest.approx(value, rel=1e-12)


@pytest.mark.slow
def test_an_ats_lake_reads_each_periods_end_by_time(runs):
    """Under ATS the steps a period takes are decided at run time, so its end is the
    record whose ``pertim`` is the period's length -- read by the targets (checked
    against the TDIS) and by every view (the last saved step), not by the declared
    ``nstp=1``."""

    model = runs["ats"]
    saved = saved_steps(model.fl.lak.output.stage())
    assert len(saved) > NPER  # ATS took several steps a period
    at = _stage_at(model)
    ends = {kper: (kstp, kper) for kstp, kper, pertim in saved if pertim == PERLEN}
    assert sorted(ends) == list(range(NPER)) and all(kstp > 0 for kstp, _ in ends.values())
    expected = [at[ends[per]] for per in range(NPER)]

    series = mf.LakeStageTargets(locations={"pond": 0}).simulated_series(model)
    assert list(series["time"]) == list(range(NPER))
    np.testing.assert_array_equal(series["pond"], expected)
    stage = model.lak.results.stage
    cells = sorted(set(stage.get()["cell"]))
    np.testing.assert_array_equal(_z(stage.map(per=1))[cells], expected[1])
    # The series takes the mean of the lake's three identical cells: one ulp off.
    np.testing.assert_allclose(stage.plot().data[0].y, expected, rtol=1e-14)
    change = model.lak.results.stage_change.get()
    np.testing.assert_array_equal(change["stage1"], expected[1:])


# -- the group and diff tiers --------------------------------------------------------


@pytest.mark.slow
def test_the_group_stage_tables_read_each_members_own_period_end(runs):
    """The group LAK stage table was keyed by saved record; its SFR twin decoded a
    reach-by-time frame as time-by-reach, so ``per`` held the reach and ``reach`` the
    record. Each member now reads its own last saved step of each period, so a
    three-step member and a 1/3/2-step member align on the period."""

    steps, uneven = runs["steps"], runs["uneven"]
    group = mf.ModelGroup({"steps": steps, "uneven": uneven}, reference="steps")
    stage = group.packages.lak.results.stage
    expected = {}
    for name, model in (("steps", steps), ("uneven", uneven)):
        at = _stage_at(model)
        ends = _last_steps(list(at))
        expected[name] = [at[ends[per]] for per in range(NPER)]
        frame = stage.get(model_name=name).sort_values("per")
        assert list(frame["per"]) == list(range(NPER))
        np.testing.assert_array_equal(frame["stage"], expected[name])
        np.testing.assert_array_equal(stage.get(model_name=name, per=1)["stage"],
                                      [expected[name][1]])
    compare = stage.compare().sort_values("per")
    np.testing.assert_allclose(compare["stage_diff"],
                               np.subtract(expected["uneven"], expected["steps"]), rtol=1e-12)
    summary = group.diff().packages.lak.results.stage.summary().iloc[0]
    assert summary["n"] == NPER and summary["argmax_per"] in range(NPER)

    outputs = group.outputs.lak.stage()
    for name, model in (("steps", steps), ("uneven", uneven)):
        clock, values = _stage_file(model)
        mine = outputs[outputs["model"] == name]
        assert list(mine["kstpkper"]) == clock
        np.testing.assert_array_equal(mine["stage"], values)

    stream = runs["stream"]
    reader = stream.fl.sfr.output.stage()
    clock = [tuple(int(v) for v in key) for key in reader.get_kstpkper()]
    data = np.asarray(reader.get_alldata(), dtype=float).reshape(len(clock), -1)
    sfr = mf.ModelGroup({"a": stream, "b": stream}, reference="a").packages.sfr.results.stage
    table = sfr.get(model_name="a")
    assert sorted(set(table["reach"])) == [0, 1, 2, 3]
    assert sorted(set(table["per"])) == list(range(NPER))
    np.testing.assert_array_equal(sfr.get(model_name="a", per=1, reach=3)["stage"],
                                  [data[clock.index((2, 1)), 3]])


@pytest.mark.slow
def test_a_group_exchange_diff_map_reads_each_members_last_step(runs):
    """The group LAK diff map summed each member's saved steps before differencing.

    ``asym`` is the reference, not ``steps``: a PEST build rewrites its model's input
    files in place (external arrays), which changes the grid signature this map
    checks, and ``steps`` is the run this module builds PEST on."""

    reference_run, uneven = runs["asym"], runs["uneven"]
    group = mf.ModelGroup({"asym": reference_run, "uneven": uneven}, reference="asym")
    mine = _lak_q(uneven, (2, 1))
    reference = _lak_q(reference_run, (2, 1))
    assert not np.isclose(sum(q for q, _ in _lak_q(uneven, (0, 1)).values()),
                          sum(q for q, _ in mine.values()))
    z = _z(group.packages.lak.results.q.compare_map(model_name="uneven", per=1))
    for cell in mine:
        expected = mine[cell][0] / mine[cell][1] - reference[cell][0] / reference[cell][1]
        assert z[cell] == pytest.approx(expected, rel=1e-9)


@pytest.mark.slow
def test_region_heads_per_reads_the_last_saved_step(runs):
    """``region_heads(per=N)`` kept every saved step of period N -- three heads per
    cell at nstp=3 -- beside ``kstpkper=``, which names one time."""

    model = runs["steps"]
    model.add_region("west", cellids=[0, 1, 2], layer=0, overwrite=True)
    heads = model.region_heads("west", per=1)
    assert set(heads["kstpkper"].map(lambda value: tuple(map(int, value)))) == {(2, 1)}
    assert sorted(heads["cell"]) == [0, 1, 2]
    same = model.region_heads("west", kstpkper=(2, 1))
    np.testing.assert_array_equal(heads["elev"].to_numpy(), same["elev"].to_numpy())
    with pytest.raises(ValueError, match="Stress period 7 is unavailable"):
        model.region_heads("west", per=7)
