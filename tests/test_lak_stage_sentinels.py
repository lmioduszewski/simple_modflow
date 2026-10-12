"""A lake with no stage is a gap in every view of the LAK stage file (ledger 329).

MF6 writes a lake's stage as a number even when it has none: -1e30 for a DRY lake
(its stage at or below the lake bottom) and +1e30 for an INACTIVE one
(``gwf-lak.f90``'s ``lak_ot_dv``). Nothing masked either, so the stage series drew
a dry lake 1e30 ft down, every map, mosaic and animation coloured from -1e30 to
+1e30, the exchange hover read "-1e+30 ft", a stage change into a dry period read
-1e30, and a group diff averaged the sentinel into a mean |Δstage| of 2.2e29. The
rule is the head file's, ``headsplus._is_sentinel`` (``|value| >= 1e29``), applied
once where the stage file is read -- ``model.outputs.lak.stage.get()`` -- as
``hds.array()`` applies it, with ``masked=False`` for the raw file.

Every MF6-backed test reads ``lak_stage_runs``: two runs of one small model. In
``dry`` lake 0 is drawn down by a withdrawal until it goes dry in periods 1-3 and
refills in period 4, lake 1 is set INACTIVE in period 2, and lake 2 is held at a
CONSTANT 108 ft. ``wet`` is the same model with no withdrawal and lake 1 always
active. What the stage file holds is asserted here, not assumed.

The stage VIEWS keep both as a gap. A lake TARGET does not (ledger 351, the user's
rule): a DRY lake is calibrated against its bottom -- 100 ft here, MF6's own STAGE
observation of it -- on both routes, and only an INACTIVE lake has no stage there
(``tests/test_dry_lake_calibration_stage.py`` covers every connection type).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations._shared import _clean_observation_output_values
from myflopy.modflow.mf6.package_plotting import (
    SpatialView,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
)

DRY_PERIODS = (1, 2, 3)  # lake 0, dry
INACTIVE_PERIOD = 2  # lake 1, STATUS INACTIVE
NPER = 6
#: Every lake's bottom by MF6's rule: one VERTICAL connection under a 100 ft top.
BOTTOM = 100.0


def _run(root, *, withdrawal: float, switch: bool):
    """One small transient model with three one-cell lakes; returns ``(model, cells)``."""

    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    cells = {"dry": interior[0], "switched": interior[5], "wet": interior[-1]}
    held = [[2, "STATUS", "CONSTANT"], [2, "STAGE", 108.0]]
    targets = mf.LakeStageTargets(locations={"drying": 0, "switched": 1, "wet": 2})
    lake = mf.lak.flopy(
        nlakes=3,
        packagedata=[[0, 100.5, 1, "drying"], [1, 104.0, 1, "switched"], [2, 108.0, 1, "wet"]],
        # VERTICAL, so each lake's bottom is its cell's top, 100 ft.
        connectiondata=[
            [0, 0, (0, cells["dry"]), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 0, (0, cells["switched"]), "VERTICAL", 0.01, 0.0, 0.0, 0.0, 0.0],
            [2, 0, (0, cells["wet"]), "VERTICAL", 0.01, 0.0, 0.0, 0.0, 0.0],
        ],
        perioddata={
            0: [[0, "WITHDRAWAL", withdrawal], [1, "STATUS", "ACTIVE"], *held],
            2: [[1, "STATUS", "INACTIVE" if switch else "ACTIVE"], *held],
            3: [[1, "STATUS", "ACTIVE"], *held],
            4: [[0, "RAINFALL", 0.5], *held],
            5: [[0, "RAINFALL", 0.0], *held],
        },
        boundnames=True,
        stage_filerecord="m.lak.stg",
        budget_filerecord="m.lak.bud",
        observations=targets.to_flopy_obs(csv_name="lake_stage_targets.csv"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=[0.0]),
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
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model(), cells


@pytest.fixture(scope="module")
def lak_stage_runs(tmp_path_factory):
    """``{"dry": model, "wet": model, "cells": {...}}`` -- see the module docstring."""

    dry, cells = _run(tmp_path_factory.mktemp("lak_dry"), withdrawal=300.0, switch=True)
    wet, _ = _run(tmp_path_factory.mktemp("lak_wet"), withdrawal=0.0, switch=False)
    return {"dry": dry, "wet": wet, "cells": cells}


def _raw(model):
    """The stage file as MF6 wrote it, asserting it holds what the fixture claims."""

    raw = model.outputs.lak.stage.get(masked=False)
    expected = np.zeros(raw.shape, dtype=bool)
    expected[list(DRY_PERIODS), 0] = True
    assert np.array_equal(raw == -1e30, expected), "lake 0 must be DRY in periods 1-3 only"
    inactive = np.zeros(raw.shape, dtype=bool)
    inactive[INACTIVE_PERIOD, 1] = True
    assert np.array_equal(raw == 1e30, inactive), "lake 1 must be INACTIVE in period 2 only"
    return raw


def _strings(customdata) -> list[str]:
    """Every entry of a trace's ``customdata``, as text."""

    return [str(value) for row in customdata for value in np.atleast_1d(row)]


def test_the_stage_array_is_nan_exactly_where_mf6_wrote_no_stage(lak_stage_runs):
    """``masked=True`` (default) is ``NaN`` at both sentinels and the file elsewhere;
    ``masked=False`` is the file, the only way left to tell dry from inactive."""

    model = lak_stage_runs["dry"]
    raw = _raw(model)
    stage = model.outputs.lak.stage.get()

    assert stage.shape == (NPER, 3)
    assert np.array_equal(np.isnan(stage), np.abs(raw) >= 1e29)
    np.testing.assert_array_equal(stage[~np.isnan(stage)], raw[np.abs(raw) < 1e29])


def test_the_stage_table_keeps_every_row_and_blanks_the_dry_ones(lak_stage_runs):
    """A dry or inactive period is a ROW with no stage, never a missing row."""

    model = lak_stage_runs["dry"]
    raw = _raw(model)
    frame = model.lak.results.stage.get()

    assert len(frame) == NPER * 3
    for (per, lake), value in frame.set_index(["per", "lake"])["stage"].items():
        if abs(raw[per, lake]) >= 1e29:
            assert np.isnan(value), (per, lake, value)
        else:
            assert value == raw[per, lake]


@pytest.mark.parametrize("per", [1, 2])
def test_a_stage_map_is_blank_where_a_lake_has_no_stage(lak_stage_runs, per):
    """The drawn values, the trace, the hover and the static norm all leave the dry
    (and, in period 2, the inactive) lake blank. The trace ran from -1e30 to +1e30."""

    model, cells = lak_stage_runs["dry"], lak_stage_runs["cells"]
    raw = _raw(model)
    real = raw[np.abs(raw) < 1e29]
    no_stage = [cells["dry"]] + ([cells["switched"]] if per == INACTIVE_PERIOD else [])

    choro = model.lak.results.stage.map(per=per, fill_value=float("nan"))
    values = np.asarray(choro.linear_zs, dtype=float)
    assert np.isnan(values[no_stage]).all()
    assert values[cells["wet"]] == 108.0
    drawn = np.asarray(choro.fig.data[0].z, dtype=float)
    assert np.nanmin(drawn) >= real.min() and np.nanmax(drawn) <= real.max()
    text = _strings(choro.fig.data[0].customdata)
    assert not [s for s in text if "e+30" in s], "the hover still prints a sentinel"

    static = model.lak.results.stage.map(per=per, backend="mpl")
    norms = [c.norm for ax in static.axes for c in ax.collections if c.get_array() is not None]
    assert norms, "the static map drew no coloured cells"
    for norm in norms:  # fill_value 0.0 is the default backdrop, a real 0
        assert 0.0 <= norm.vmin and norm.vmax <= real.max(), (norm.vmin, norm.vmax)


@pytest.mark.parametrize("agg", ["first", "mean", "min", "max", "last", "sum"])
def test_a_stage_map_is_blank_there_under_every_documented_aggregation(lak_stage_runs, agg):
    """``agg="sum"`` painted the dry lake's cell at 0.0 ft and ran the colour range
    from 0 -- pandas sums an all-NaN group to 0.0, which the map payload took."""

    model, cells = lak_stage_runs["dry"], lak_stage_runs["cells"]
    raw = _raw(model)
    real = raw[np.abs(raw) < 1e29]
    choro = model.lak.results.stage.map(per=1, agg=agg, fill_value=float("nan"))
    values = np.asarray(choro.linear_zs, dtype=float)

    assert np.isnan(values[cells["dry"]]), values[cells["dry"]]
    assert real.min() <= np.nanmin(values) and np.nanmax(values) <= real.max()


def test_a_map_payload_leaves_a_cell_with_no_value_blank():
    """A cell PRESENT in the table whose rows hold no value is NaN under every
    aggregation but a counting one; an absent cell still takes ``fill_value``."""

    frame = pd.DataFrame({"package": "lak", "cell": [0, 0, 1, 1],
                          "stage": [5.0, 5.0, np.nan, np.nan]})
    for agg in ("sum", "mean", "first", "max"):
        values, _ = build_cell_input_map_payload(
            frame, ncpl=3, value_column="stage", per=0, layer=0, fill_value=-7.0, agg=agg)
        assert np.isnan(values[1]) and values[2] == -7.0, (agg, values)
    values, _ = build_cell_input_map_payload(
        frame, ncpl=3, value_column="stage", per=0, layer=0, agg="count")
    assert values == [2.0, 0.0, 0.0]

    compare = frame.assign(stage_diff=frame["stage"] - 5.0, model="b", reference_model="a")
    diffs, _, _ = build_group_input_compare_map_payload(
        compare, ncpl=3, value_column="stage", diff_column="stage_diff", model_name="b",
        reference_model="a", per=0, layer=0, agg="sum")
    assert diffs[0] == 0.0 and np.isnan(diffs[1]), diffs  # undefined is not "no difference"


def test_mosaic_and_animate_colour_real_stages_only(lak_stage_runs):
    """The pooled colour axis of both composers ran from -1e30 to +1e30."""

    model = lak_stage_runs["dry"]
    raw = _raw(model)
    top = raw[np.abs(raw) < 1e29].max()
    stage = model.lak.results.stage
    for figure in (stage.mosaic(per=INACTIVE_PERIOD), stage.animate()):
        axis = figure.layout.coloraxis
        assert 0.0 <= axis.cmin and axis.cmax <= top, (axis.cmin, axis.cmax)


@pytest.mark.parametrize("backend", ["plotly", "mpl"])
def test_a_stage_series_leaves_a_gap_on_a_day_with_no_stage(lak_stage_runs, backend):
    """The user's report: ``lak.results.stage.plot()`` drew the dry days at -1e30."""

    model = lak_stage_runs["dry"]
    raw = _raw(model)
    figure = model.lak.results.stage.plot(backend=backend)
    if backend == "plotly":
        lines = {trace.name: np.asarray(trace.y, dtype=float) for trace in figure.data}
    else:
        lines = {line.get_label(): np.asarray(line.get_ydata(), dtype=float)
                 for line in figure.axes[0].get_lines()}
        low, high = figure.axes[0].get_ylim()
        assert 90.0 < low and high < 110.0, (low, high)
    for lake in range(3):
        expected = np.where(np.abs(raw[:, lake]) >= 1e29, np.nan, raw[:, lake])
        np.testing.assert_array_equal(lines[f"Lake {lake}"], expected)


def test_a_stage_change_into_or_out_of_a_dry_period_is_kept_and_undefined(lak_stage_runs):
    """Every transition stays a row; one touching a period with no stage has no
    change. Dropping on ``stage1`` deleted the transitions INTO a dry period while
    keeping the way back out, so a lake that went dry looked as if it never did."""

    model = lak_stage_runs["dry"]
    raw = _raw(model)
    frame = model.lak.results.stage_change.get()

    assert len(frame) == 3 * (NPER - 1)
    for row in frame.itertuples():
        before, after = raw[row.per0, row.lake], raw[row.per1, row.lake]
        if max(abs(before), abs(after)) >= 1e29:
            assert np.isnan(row.stage_change), row
        else:
            assert row.stage_change == pytest.approx(after - before)
    figure = model.lak.results.stage_change.plot()
    assert all(np.nanmax(np.abs(np.asarray(t.y, dtype=float))) < 10.0 for t in figure.data)


def test_the_exchange_map_hover_has_no_stage_at_a_dry_lake(lak_stage_runs):
    """``results.q.map`` joins the lake stage into its hover; it read "-1e+30 ft"."""

    model, cells = lak_stage_runs["dry"], lak_stage_runs["cells"]
    _raw(model)
    choro = model.lak.results.q.map(per=1)
    trace = choro.fig.data[0]
    assert not [s for s in _strings(trace.customdata) if "e+30" in s]
    labels = [piece.split("%{customdata[")[0] for piece in trace.hovertemplate.split("<br>")]
    stage_slot = next(i for i, label in enumerate(labels) if ">stage<" in label)
    assert trace.customdata[cells["dry"]][stage_slot] == ""
    assert trace.customdata[cells["switched"]][stage_slot].startswith("103.")


def test_lake_stage_targets_read_an_inactive_lake_as_no_stage(lak_stage_runs):
    """The observation CSV route read MF6's 3e30 for an inactive lake as a stage of
    0.0 ft. A DRY lake there is MF6's own observed stage -- its bottom, 100 ft --
    which the stage-file route now reads too, where the file writes -1e30 (ledger
    351)."""

    model = lak_stage_runs["dry"]
    _raw(model)
    targets = mf.LakeStageTargets(locations={"drying": 0, "switched": 1, "wet": 2})
    series = targets.simulated_series(model)

    assert np.isnan(series["switched"].iloc[INACTIVE_PERIOD])
    assert series["switched"].drop(index=INACTIVE_PERIOD).notna().all()
    np.testing.assert_array_equal(series["drying"].iloc[list(DRY_PERIODS)], 100.0)
    compare = targets.compare(model)
    assert compare["sim_stage"].abs().max() < 1e3
    assert np.isnan(compare.loc[compare["name"] == "switched", "sim_stage"].iloc[INACTIVE_PERIOD])


def test_lake_stage_targets_keep_every_record_on_the_stage_file_route(lak_stage_runs):
    """The stage-file route DROPPED a record with no stage, so the wide series was
    3 rows for 6 records -- and 5 for two lakes, since period 2 has none for either.
    Names the observation CSV does not hold, so the route is the stage file's.

    A DRY record is the lake's bottom there now, and only an INACTIVE one is ``NaN``
    (ledger 351): the file's -1e30 and +1e30 are told apart, where the masked table
    this route read merged them."""

    model = lak_stage_runs["dry"]
    raw = _raw(model)
    for locations in ({"foo": 0}, {"foo": 0, "bar": 1}):
        series = mf.LakeStageTargets(locations=locations).simulated_series(model)
        assert list(series["time"]) == list(range(NPER))
        for name, lake in locations.items():
            expected = np.where(raw[:, lake] == -1e30, BOTTOM, raw[:, lake])
            expected = np.where(expected == 1e30, np.nan, expected)
            np.testing.assert_array_equal(series[name], expected)
    compare = mf.LakeStageTargets(locations={"foo": 0}).compare(model)
    assert len(compare) == NPER and compare["sim_stage"].notna().all()
    np.testing.assert_array_equal(
        compare.sort_values("time")["sim_stage"].iloc[list(DRY_PERIODS)], BOTTOM)


def test_a_lake_inactive_in_every_record_keeps_its_column_and_its_rows(monkeypatch):
    """The pivot's default drops a column that is all NaN and a row that is all NaN:
    a lake with no stage for the whole base run would lose its column, so a later
    realization in which it holds water would write one column MORE than the
    instruction file reads. Only an INACTIVE lake has no stage for a target (ledger
    351; a dry one reads its bottom). A double handing over the RAW stage table, in
    MF6's +1e30; it declares no observation, so the stage-file route."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.observations import lake as lake_module

    table = pd.DataFrame({"per": [0, 0, 1, 1, 2, 2], "lake": [0, 1] * 3, "cell": [4, 9] * 3,
                          "stage": [1e30, 101.0, 1e30, 1e30, 1e30, 102.0]})

    def build(model, *, masked=True):
        assert masked is False
        return table.copy()

    monkeypatch.setattr(lake_module, "build_lak_stage_result_table", build)
    series = mf.LakeStageTargets(locations={"always_off": 0, "wet": 1}).simulated_series(
        SimpleNamespace())

    assert list(series.columns) == ["time", "always_off", "wet"]
    assert list(series["time"]) == [0, 1, 2]
    assert series["always_off"].isna().all()
    np.testing.assert_array_equal(series["wet"], [101.0, np.nan, 102.0])


@pytest.mark.parametrize("names", [("drying", "switched"), ("foo", "bar")],
                         ids=["observation-route", "stage-file-route"])
def test_the_files_pest_reads_agree_row_for_row_and_have_no_empty_field(
        lak_stage_runs, tmp_path, names):
    """PEST's instruction file is built from the BASE run, one line per row. The
    forward run wrote fewer rows (EOF: the realization fails) and empty fields for a
    lake with no stage (PEST++ refuses those, and ``nan``). Both files now carry
    every record, and 1e30 where the series says NaN -- an INACTIVE lake only: a DRY
    one is written as its bottom on both routes (ledger 351)."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
    from myflopy.modflow.mf6.pest.observations import prepare_lake_stage_observations
    from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec

    model = lak_stage_runs["dry"]
    _raw(model)
    # A template (no values): the FILES are what this pins, not which rows a target
    # names. A template used to die in the build on a missing ``time`` column (ledger
    # 329); it registers now.
    targets = mf.LakeStageTargets(locations=dict(zip(names, (0, 1), strict=True)))
    series = targets.simulated_series(model)
    project = SimpleNamespace(model=model, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_lake_stage_observations(
        project, LakeStageObservationSpec(targets=targets, prefix="lk"))
    base = tmp_path / prepared["named_series_forward_run_config"]["output_csv"]
    forward = tmp_path / "forward.csv"
    config = prepared["named_series_forward_run_config"]
    write_named_series_targets(
        [{"kind": "lake_stage", "locations_file": tmp_path / config["locations_file"],
          "output_csv": forward}],
        sim_ws=model.workspace)

    assert forward.read_text() == base.read_text()
    lines = base.read_text().splitlines()
    assert len(lines) == 1 + NPER
    assert not [line for line in lines if ",," in line or line.endswith(",")]
    written = pd.read_csv(base)
    no_stage = series[list(names)].isna().to_numpy()
    assert no_stage.any()
    np.testing.assert_array_equal(written[list(names)].to_numpy()[no_stage], 1e30)
    np.testing.assert_allclose(written[list(names)].to_numpy()[~no_stage],
                               series[list(names)].to_numpy()[~no_stage], rtol=1e-14)
    np.testing.assert_array_equal(written[names[0]].iloc[list(DRY_PERIODS)], BOTTOM)
    assert written[names[1]].iloc[INACTIVE_PERIOD] == 1e30
    assert (written.to_numpy() == 1e30).sum() == 1


def test_an_observation_flow_keeps_zero_and_a_stage_asks_for_nan():
    """``no_value`` is what MF6's 3e30 becomes: 0.0 for a flow, as it always was."""

    values = pd.Series([1.5, 3.0e30, -2.0])
    np.testing.assert_array_equal(_clean_observation_output_values(values), [1.5, 0.0, -2.0])
    cleaned = _clean_observation_output_values(values, no_value=np.nan)
    assert cleaned.iloc[0] == 1.5 and np.isnan(cleaned.iloc[1]) and cleaned.iloc[2] == -2.0


def test_the_group_stage_views_leave_the_same_gaps(lak_stage_runs):
    """The group reads the same stage file. Its series engine SUMS by default, and
    pandas sums an all-NaN group to 0.0, so once masked it drew a dry lake at 0 ft."""

    group = mf.ModelGroup({"wet": lak_stage_runs["wet"], "dry": lak_stage_runs["dry"]},
                          reference="wet")
    raw = _raw(lak_stage_runs["dry"])
    stage = group.packages.lak.results.stage
    frame = stage.get(model_name="dry")
    assert frame["stage"].isna().sum() == 4

    lines = {trace.name: np.asarray(trace.y, dtype=float) for trace in stage.plot().data}
    for lake in range(3):
        expected = np.where(np.abs(raw[:, lake]) >= 1e29, np.nan, raw[:, lake])
        np.testing.assert_array_equal(lines[f"dry / Lake {lake}"], expected)
    outputs = group.outputs.lak.stage()
    assert outputs["stage"].abs().max() < 1e3 and outputs["stage"].isna().sum() == 4
    values = np.asarray(stage.map("dry", per=INACTIVE_PERIOD).linear_zs, dtype=float)
    assert np.nanmin(values) >= 0.0 and np.nanmax(values) <= 108.0


def test_a_stage_diff_counts_a_one_sided_dry_lake_instead_of_averaging_it(lak_stage_runs):
    """The diff summary averaged the sentinel: max |Δstage| 1e30, mean 2.2e29. A dry
    lake against a wet one is a difference no number measures, so it is COUNTED
    (``n_one_sided``) and fails the tolerance; two dry lakes agree."""

    dry, wet = lak_stage_runs["dry"], lak_stage_runs["wet"]
    _raw(dry)
    diff = mf.ModelGroup({"wet": wet, "dry": dry}, reference="wet").diff()
    summary = diff.packages.lak.results.stage.summary().iloc[0]
    both = diff.packages.lak.results.stage.get().dropna(subset=["stage", "reference_stage"])

    assert summary["n_one_sided"] == 4 and summary["n"] == len(both) == 14
    assert summary["max_abs_diff"] == pytest.approx(both["stage_diff"].abs().max())
    assert summary["max_abs_diff"] < 1.0
    assert not summary["within_tolerance"]

    same = mf.ModelGroup({"a": dry, "b": dry}, reference="a").diff()
    agreed = same.packages.lak.results.stage.summary().iloc[0]
    assert agreed["n_one_sided"] == 0 and agreed["n"] == 14 and agreed["within_tolerance"]


def test_a_one_sided_row_fails_the_tolerance_even_when_every_other_row_agrees():
    """A lake dry in one model only is not "within tolerance" because the rows that
    CAN be compared agree. Built on a compare table, so nothing else can differ."""

    from types import SimpleNamespace

    from myflopy.project.model_results_diff import StageResultDiff

    table = pd.DataFrame({
        "model": "b", "reference_model": "a", "per": [0, 1, 2], "lake": 0,
        "stage": [101.0, np.nan, np.nan], "reference_stage": [101.0, 100.5, np.nan],
    })
    table["stage_diff"] = table["stage"] - table["reference_stage"]
    accessor = SimpleNamespace(compare=lambda **_: table.copy())
    group = SimpleNamespace(models={"a": None, "b": None}, reference="a")
    row = StageResultDiff(SimpleNamespace(group=group), accessor, entity="lake").summary().iloc[0]

    assert (row["n"], row["n_one_sided"]) == (1, 1)
    assert row["max_abs_diff"] == 0.0 and not row["within_tolerance"]
    table.loc[1, "reference_stage"] = np.nan  # both dry: they agree
    row = StageResultDiff(SimpleNamespace(group=group), accessor, entity="lake").summary().iloc[0]
    assert (row["n"], row["n_one_sided"]) == (1, 0) and row["within_tolerance"]
    table.loc[0, "stage"] = np.nan  # nothing left to compare on either side
    table.loc[0, "reference_stage"] = np.nan
    row = StageResultDiff(SimpleNamespace(group=group), accessor, entity="lake").summary().iloc[0]
    assert row["n"] == 0 and np.isnan(row["max_abs_diff"]) and row["argmax_lake"] == -1


def test_a_series_point_with_no_value_stays_a_gap_under_any_aggregation():
    """``sum`` of an all-NaN group is 0.0 in pandas; a counting aggregation's 0 stays."""

    frame = pd.DataFrame({"lake": [0, 0, 0, 0], "per": [0, 0, 1, 1],
                          "stage": [5.0, 5.0, np.nan, np.nan]})
    for agg in ("sum", "mean", "first", "max"):
        (_, _, y), = SpatialView._series_lines(SpatialView(), frame, ["lake"], "stage", agg)
        assert not np.isnan(y[0]) and np.isnan(y[1]), (agg, y)
    (_, _, y), = SpatialView._series_lines(SpatialView(), frame, ["lake"], "stage", "count")
    np.testing.assert_array_equal(y, [2.0, 0.0])
