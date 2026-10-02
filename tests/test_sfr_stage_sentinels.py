"""A reach with no stage is a gap in every view of the SFR stage file (ledger 353).

MF6 writes a reach's stage as a number even when it has none: -1e30 for a DRY reach
(its depth is zero) and +1e30 for an INACTIVE one (``gwf-sfr.f90``'s ``sfr_ot_dv``)
-- the two numbers ledger 329 masked for LAK. Nothing masked them for SFR, so the
stage series drew a dry reach 1e30 ft down, every map, mosaic and animation coloured
from -1e30 to +1e30, the long profile and the stage profile plotted the sentinel,
the exchange hover read "-1e+30 ft", ``canonical_sfr_signals`` counted an inactive
reach WET at a depth of 1e30, and a group diff averaged the sentinel. The user's
call: "the plots for SFR should just not plot a stage if the reach is dry". The rule
is the head file's, ``headsplus._is_sentinel`` (``|value| >= 1e29``), applied once
where the stage file is read -- ``model.outputs.sfr.stage.get()`` -- with
``masked=False`` for the raw file, exactly as LAK's.

**A reach TARGET is NOT masked** (the user: leave SFR targets as they are):
``SfrStageTargets``' stage-file route reads the file raw, so a dry reach still reads
-1e30 there and in the files PEST reads. The last three tests pin that.

Every MF6-backed test reads ``sfr_stage_runs``: two runs of one small model -- a 6 x
6 irregular Voronoi grid of ~100 ft cells, one confined layer (top 110, bottom 0, K
10), CHD at 90 ft on the west and east edges, six 10-day periods of one step, and
five reaches whose beds (top 100 ft, 1 ft thick) sit above the water table, so a
reach with no inflow has no water:

  reaches 0-1-2  one stream, 50 ft3/d into reach 0; in ``dry`` the inflow stops in
                 periods 1-3, so all three are DRY there, and resumes in period 4.
  reach 3        alone, 30 ft3/d; STATUS INACTIVE in period 2 of ``dry``.
  reach 4        alone, 100 ft3/d throughout: a control.

``wet`` is the same model with the inflow never stopped and reach 3 always active.
SFR STAGE observations are declared under ``r0`` .. ``r4``: those names take the
observation route, any other name the stage file. What the stage file holds is
asserted here, not assumed.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import canonical_sfr_signals, irregular_voronoi_grid
from myflopy.modflow.mf6.package_budget import build_sfr_stage_result_table

NPER = 6
DRY_PERIODS = (1, 2, 3)  # reaches 0-2
DRYING = (0, 1, 2)
INACTIVE_PERIOD = 2  # reach 3
INACTIVE = 3
NREACH = 5
#: The model cell of each reach -- one interior cell each, none on a CHD edge.
CELLS = (13, 14, 15, 20, 27)
#: Each reach's streambed top (``rtp``): 100.0 falling 0.1 ft a reach.
BED_TOP = tuple(100.0 - 0.1 * reach for reach in range(NREACH))
OBSERVED = {f"r{reach}": reach for reach in range(NREACH)}


def _run(root, *, drying: bool):
    """One small transient model with the five reaches of the module docstring."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    assert not set(CELLS) & set(edge.tolist()), "a reach sits on a CHD cell"
    ncon = (1, 2, 1, 0, 0)
    stream = mf.sfr.flopy(
        # reach, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, man, ncon, ustrf, ndv
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
        stage_filerecord="m.sfr.stg",
        budget_filerecord="m.sfr.bud",
        observations=mf.SfrStageTargets(locations=OBSERVED).to_flopy_obs(
            csv_name="sfr_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        stream,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 1, 1.0)] * NPER)])
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


@pytest.fixture(scope="module")
def sfr_stage_runs(tmp_path_factory):
    """``{"dry": model, "wet": model}`` -- see the module docstring."""

    return {
        "dry": _run(tmp_path_factory.mktemp("sfr_dry"), drying=True),
        "wet": _run(tmp_path_factory.mktemp("sfr_wet"), drying=False),
    }


def _raw(model) -> np.ndarray:
    """The stage file as MF6 wrote it, ``(reach, record)``, asserting the fixture's claims."""

    raw = model.outputs.sfr.stage.get(masked=False).to_numpy(dtype=float)
    assert raw.shape == (NREACH, NPER)
    dry = np.zeros(raw.shape, dtype=bool)
    dry[np.ix_(DRYING, DRY_PERIODS)] = True
    assert np.array_equal(raw == -1e30, dry), "reaches 0-2 must be DRY in periods 1-3 only"
    inactive = np.zeros(raw.shape, dtype=bool)
    inactive[INACTIVE, INACTIVE_PERIOD] = True
    assert np.array_equal(raw == 1e30, inactive), "reach 3 must be INACTIVE in period 2 only"
    return raw


def _no_stage(raw) -> np.ndarray:
    """``raw`` with both sentinels ``NaN`` -- what every view should draw."""

    return np.where(np.abs(raw) >= 1e29, np.nan, raw)


def _strings(customdata) -> list[str]:
    """Every entry of a trace's ``customdata``, as text."""

    return [str(value) for row in customdata for value in np.atleast_1d(row)]


# -- the views: a gap --------------------------------------------------------------


def test_the_stage_file_is_nan_exactly_where_mf6_wrote_no_stage(sfr_stage_runs):
    """``masked=True`` (default) is ``NaN`` at both sentinels and the file elsewhere;
    ``masked=False`` is the file, the only way left to tell dry from inactive."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    stage = model.outputs.sfr.stage.get()

    assert list(stage.columns) == [(0, per) for per in range(NPER)]
    np.testing.assert_array_equal(stage.to_numpy(dtype=float), _no_stage(raw))
    unmasked = build_sfr_stage_result_table(model, masked=False)
    np.testing.assert_array_equal(
        unmasked.pivot(index="reach", columns="per", values="stage").to_numpy(), raw)


def test_the_stage_table_keeps_every_row_and_blanks_the_dry_ones(sfr_stage_runs):
    """A dry or inactive period is a ROW with no stage, never a missing row."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    frame = model.sfr.results.stage.get()

    assert len(frame) == NPER * NREACH
    for (per, reach), value in frame.set_index(["per", "reach"])["stage"].items():
        if abs(raw[reach, per]) >= 1e29:
            assert np.isnan(value), (per, reach, value)
        else:
            assert value == raw[reach, per]


@pytest.mark.parametrize("per", [1, 2])
def test_a_stage_map_is_blank_where_a_reach_has_no_stage(sfr_stage_runs, per):
    """The drawn values, the trace, the hover and the static norm all leave the dry
    (and, in period 2, the inactive) reaches blank. The trace ran from -1e30."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    real = raw[np.abs(raw) < 1e29]
    no_stage = [CELLS[reach] for reach in range(NREACH) if abs(raw[reach, per]) >= 1e29]

    choro = model.sfr.results.stage.map(per=per, fill_value=float("nan"))
    values = np.asarray(choro.linear_zs, dtype=float)
    assert len(no_stage) == (4 if per == INACTIVE_PERIOD else 3)
    assert np.isnan(values[no_stage]).all()
    assert values[CELLS[4]] == raw[4, per]
    drawn = np.asarray(choro.fig.data[0].z, dtype=float)
    assert np.nanmin(drawn) >= real.min() and np.nanmax(drawn) <= real.max()
    assert not [s for s in _strings(choro.fig.data[0].customdata) if "e+30" in s]

    static = model.sfr.results.stage.map(per=per, backend="mpl")
    norms = [c.norm for ax in static.axes for c in ax.collections if c.get_array() is not None]
    assert norms, "the static map drew no coloured cells"
    for norm in norms:  # fill_value 0.0 is the default backdrop, a real 0
        assert 0.0 <= norm.vmin and norm.vmax <= real.max(), (norm.vmin, norm.vmax)


@pytest.mark.parametrize("agg", ["first", "mean", "min", "max", "last", "sum"])
def test_a_stage_map_is_blank_there_under_every_documented_aggregation(sfr_stage_runs, agg):
    """``agg="sum"`` would paint a masked reach 0.0 ft (pandas sums an all-NaN group
    to 0.0); the map payload's ``_aggregate_leaving_gaps`` keeps it blank."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    real = raw[np.abs(raw) < 1e29]
    choro = model.sfr.results.stage.map(per=1, agg=agg, fill_value=float("nan"))
    values = np.asarray(choro.linear_zs, dtype=float)

    assert np.isnan(values[[CELLS[reach] for reach in DRYING]]).all()
    assert real.min() <= np.nanmin(values) and np.nanmax(values) <= real.max()


def test_mosaic_and_animate_colour_real_stages_only(sfr_stage_runs):
    """The pooled colour axis of both composers ran from -1e30 to +1e30."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    top = raw[np.abs(raw) < 1e29].max()
    stage = model.sfr.results.stage
    for figure in (stage.mosaic(per=INACTIVE_PERIOD), stage.animate()):
        axis = figure.layout.coloraxis
        assert 0.0 <= axis.cmin and axis.cmax <= top, (axis.cmin, axis.cmax)


@pytest.mark.parametrize("backend", ["plotly", "mpl"])
def test_a_stage_series_leaves_a_gap_on_a_day_with_no_stage(sfr_stage_runs, backend):
    """The user's words: the plots should not draw a stage where the reach is dry."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    figure = model.sfr.results.stage.plot(backend=backend)
    if backend == "plotly":
        lines = {trace.name: np.asarray(trace.y, dtype=float) for trace in figure.data}
    else:
        lines = {line.get_label(): np.asarray(line.get_ydata(), dtype=float)
                 for line in figure.axes[0].get_lines()}
        low, high = figure.axes[0].get_ylim()
        assert 90.0 < low and high < 110.0, (low, high)
    for reach in range(NREACH):
        np.testing.assert_array_equal(lines[f"Reach {reach}"], _no_stage(raw)[reach])


def test_the_long_profile_and_the_stage_profile_leave_a_dry_reach_blank(sfr_stage_runs):
    """``results.profile`` and ``results.stage.profile`` drew the sentinel as a stage
    point: the long profile's y-axis ran down to -1e30 (and up to +1e30 at the
    inactive reach in period 2)."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    for per in DRY_PERIODS:
        expected = _no_stage(raw)[:, per]
        profile = model.sfr.results.profile.get(per=per).sort_values("reach")
        np.testing.assert_array_equal(profile["stage"].to_numpy(dtype=float), expected)
        figure = model.sfr.results.profile.plot(per=per)
        stage_trace = next(trace for trace in figure.data if trace.name == "Stage")
        np.testing.assert_array_equal(np.asarray(stage_trace.y, dtype=float), expected)

        field = model.sfr.results.stage.profile
        np.testing.assert_array_equal(
            field.get(per=per)["stage"].to_numpy(dtype=float), expected)
        (line,) = field.plot(per=per).data
        np.testing.assert_array_equal(np.asarray(line.y, dtype=float), expected)


@pytest.mark.parametrize("per", [1, 2])
def test_the_exchange_map_hover_has_no_stage_at_a_reach_without_one(sfr_stage_runs, per):
    """``results.q.map`` joins the reach stage into its hover; it read "-1e+30 ft"
    for a dry reach and "1e+30 ft" for an inactive one."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    trace = model.sfr.results.q.map(per=per).fig.data[0]
    assert not [s for s in _strings(trace.customdata) if "e+30" in s]
    labels = [piece.split("%{customdata[")[0] for piece in trace.hovertemplate.split("<br>")]
    stage_slot = next(i for i, label in enumerate(labels) if ">stage<" in label)
    for reach, cell in enumerate(CELLS):
        shown = trace.customdata[cell][stage_slot]
        if abs(raw[reach, per]) >= 1e29:
            assert shown == "", (reach, shown)
        else:
            assert shown.startswith(f"{raw[reach, per]:.4g}"), (reach, shown)


def test_an_inactive_reach_is_not_counted_wet_by_the_stream_signals(sfr_stage_runs):
    """``canonical_sfr_signals`` called a reach wet when its stage was above -1e20, so
    an INACTIVE reach's +1e30 was wet at a depth of 1e30 (its ``maximum_depth``)."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    signals = canonical_sfr_signals(model, per=INACTIVE_PERIOD)

    assert (signals["wet_reach_count"], signals["dry_reach_count"]) == (1, 4)
    assert signals["maximum_depth"] == pytest.approx(raw[4, INACTIVE_PERIOD] - BED_TOP[4])


def test_the_group_stage_views_leave_the_same_gaps(sfr_stage_runs):
    """The group reads the same stage file; its series engine SUMS by default, and
    pandas sums an all-NaN group to 0.0."""

    group = mf.ModelGroup({"wet": sfr_stage_runs["wet"], "dry": sfr_stage_runs["dry"]},
                          reference="wet")
    raw = _raw(sfr_stage_runs["dry"])
    stage = group.packages.sfr.results.stage
    frame = stage.get(model_name="dry")
    assert len(frame) == NPER * NREACH and frame["stage"].isna().sum() == 10

    lines = {trace.name: np.asarray(trace.y, dtype=float) for trace in stage.plot().data}
    for reach in range(NREACH):
        np.testing.assert_array_equal(lines[f"dry / Reach {reach}"], _no_stage(raw)[reach])
    values = np.asarray(stage.map("dry", per=INACTIVE_PERIOD).linear_zs, dtype=float)
    assert np.nanmin(values) >= 0.0 and np.nanmax(values) <= raw[np.abs(raw) < 1e29].max()


def test_a_stage_diff_counts_a_one_sided_dry_reach_instead_of_averaging_it(sfr_stage_runs):
    """A reach dry in one model and wet in the other is COUNTED (``n_one_sided``) and
    fails the tolerance; two dry reaches agree. The sentinel was averaged into the
    statistics: a max |Δstage| of 1e30."""

    dry, wet = sfr_stage_runs["dry"], sfr_stage_runs["wet"]
    _raw(dry)
    diff = mf.ModelGroup({"wet": wet, "dry": dry}, reference="wet").diff()
    summary = diff.packages.sfr.results.stage.summary().iloc[0]
    both = diff.packages.sfr.results.stage.get().dropna(subset=["stage", "reference_stage"])

    assert summary["n_one_sided"] == 10 and summary["n"] == len(both) == 20
    assert summary["max_abs_diff"] == pytest.approx(both["stage_diff"].abs().max())
    assert summary["max_abs_diff"] < 1.0
    assert not summary["within_tolerance"]

    same = mf.ModelGroup({"a": dry, "b": dry}, reference="a").diff()
    agreed = same.packages.sfr.results.stage.summary().iloc[0]
    assert agreed["n_one_sided"] == 0 and agreed["n"] == 20 and agreed["within_tolerance"]


# -- the targets: unchanged (the user's call) ---------------------------------------


def test_a_reach_target_reads_the_raw_stage_file_as_before(sfr_stage_runs):
    """The stage-file route (names no OBS declares) reads the file RAW: a dry reach
    is -1e30 and an inactive one +1e30, in the series and in ``compare()``, exactly as
    before the views were masked -- the user kept SFR targets as they were."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    locations = {"a": 0, "b": INACTIVE, "c": 4}
    series = mf.SfrStageTargets(locations=locations).simulated_series(model)
    assert list(series["time"]) == list(range(NPER))
    for name, reach in locations.items():
        np.testing.assert_array_equal(series[name], raw[reach])

    compare = mf.SfrStageTargets(
        locations={"a": 0}, values=pd.DataFrame({"time": range(NPER), "a": 101.0}),
    ).compare(model).sort_values("per")
    np.testing.assert_array_equal(compare["sim_stage"], raw[0])
    assert (compare["sim_stage"] == -1e30).sum() == len(DRY_PERIODS)


def test_the_files_pest_reads_still_hold_the_raw_stage_at_a_dry_reach(sfr_stage_runs, tmp_path):
    """The base run's CSV and a forward run's agree row for row and carry MF6's -1e30
    at the dry reach and +1e30 at the inactive one -- what they carried before."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
    from myflopy.modflow.mf6.pest.observations import prepare_sfr_stage_observations
    from myflopy.modflow.mf6.pest.specs import SfrStageObservationSpec

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    targets = mf.SfrStageTargets(
        locations={"a": 0, "b": INACTIVE},
        values=pd.DataFrame({"time": range(NPER), "a": 101.0, "b": 101.0}))
    project = SimpleNamespace(model=model, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_sfr_stage_observations(
        project, SfrStageObservationSpec(targets=targets, prefix="st"))
    config = prepared["named_series_forward_run_config"]
    base = tmp_path / config["output_csv"]
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        "sfr_stage", tmp_path / config["locations_file"], forward, sim_ws=model.workspace)

    assert forward.read_text() == base.read_text()
    written = pd.read_csv(base, float_precision="round_trip")
    np.testing.assert_array_equal(written["a"], raw[0])
    np.testing.assert_array_equal(written["b"], raw[INACTIVE])


def test_the_observation_route_reads_the_streambed_top_where_the_file_has_no_stage(
        sfr_stage_runs):
    """MF6's SFR STAGE observation writes the reach's streambed top for a DRY reach
    and -- measured, not read in the source -- for an INACTIVE one, where the stage
    file writes -1e30 and +1e30: the two routes of one target differ there, as they
    did before (no rule was decided for a reach)."""

    model = sfr_stage_runs["dry"]
    raw = _raw(model)
    series = mf.SfrStageTargets(locations=OBSERVED).simulated_series(model)

    for reach in DRYING:
        np.testing.assert_allclose(
            series[f"r{reach}"].iloc[list(DRY_PERIODS)], BED_TOP[reach], rtol=0, atol=1e-9)
    assert series[f"r{INACTIVE}"].iloc[INACTIVE_PERIOD] == pytest.approx(BED_TOP[INACTIVE])
    wet = np.abs(raw) < 1e29
    np.testing.assert_allclose(
        series[[f"r{reach}" for reach in range(NREACH)]].to_numpy().T[wet], raw[wet],
        rtol=1e-9)
