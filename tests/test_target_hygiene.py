"""Observation-side hygiene before new target families land (ledgers 499 .. 506).

Every test runs on the real MF6 toy of ``tests/test_pest_handed_targets.py`` -- a
6 x 6 Voronoi grid, one layer, three 10-day periods of three steps, a four-reach
stream, a three-cell DRN zone, a two-cell lake named ``pond``, no MF6 OBS package --
so every family reads its FILE route, keyed by zero-based stress period.

What each test pins, in the order the defects were found:

* a date, or a duration, in a target's key was read as a period numbered in
  nanoseconds (ledger 499);
* ``compare()``/``stats()`` on every FILE route answered ``NaN`` for a key past the
  run's periods and dropped it from ``stats()`` (ledger 500);
* a blank, unreadable or fractional well layer became layer 0 (ledger 501),
  and a name repeated in one layer two full-weight observations (ledger 502);
* ``cal.observe(model.targets.heads)`` raised, and a point target's ``compare()``
  returned ``x_x``/``y_x``;
* a lake or reach id was ``int()``-ed and never checked, and a boundname refused
  (ledger 503); a lake or reach given no value lost its column in every PEST
  realization;
* a per-reading ``weight`` column was dropped, and SFR sets had no weight (ledger
  504);
* the PEST head writer read ``<model>.hds`` whatever the OC declared (ledger
  505), and a cell the base run left empty built, whatever its weight (ledger
  506); a DRN reading's weight beside its zone's split into `weight_x`/`weight_y`.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.specs import (
    HeadTargetObservationSpec,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
)
from tests.test_pest_handed_targets import (
    NPER,
    NSTP,
    PERLEN,
    REACHES,
    _build,
    _disv,
    _grid_and_edges,
    _head_targets,
    _observation_at,
    _run,
    run_toy,
)

pyemu = pytest.importorskip("pyemu")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The toy model of ``test_pest_handed_targets``, run."""

    return run_toy(tmp_path_factory.mktemp("target_hygiene"))


# -- ledger 499: a date is not a period ---------------------------------

_DATED = r"reads `time` as a zero-based stress period NUMBER, and {} (is|are) a date or a duration"


def test_a_date_key_is_refused_by_compare_stats_to_wide_and_the_build(toy, tmp_path):
    """pandas reads a datetime as nanoseconds since 1970 -- finite and whole -- so the
    target was keyed by period 1706659200000000000: compare() NaN, stats() n = 0."""

    dated = pd.DataFrame({"time": pd.to_datetime(["2024-01-31", "2024-02-20"]), "W": [95.0, 96.0]})
    heads = _head_targets([("W", 22)], dated)
    for call in (lambda: heads.compare(toy), lambda: heads.stats(toy), heads.to_wide):
        with pytest.raises(ValueError, match=_DATED.format("2024-01-31, 2024-02-20")):
            call()
    with pytest.raises(ValueError, match=_DATED.format("2024-01-31, 2024-02-20")):
        _build(toy, tmp_path / "pest", observe=(heads,))

    # A Series indexed by date reaches the same table as Timestamp OBJECTS.
    series = pd.Series([99.0], index=pd.DatetimeIndex(["2024-01-31"], name="time"), name="pond")
    pond = mf.LakeStageTargets(locations={"pond": 0}, values=series.reset_index().astype(object))
    with pytest.raises(ValueError, match=_DATED.format("2024-01-31")):
        pond.compare(toy)
    # A duration too, which pandas reads as nanoseconds the same way.
    flow = mf.SfrFlowTargets(locations={"a": 0},
                             values=pd.DataFrame({"time": pd.to_timedelta([10], unit="D"), "a": [-7.0]}))
    with pytest.raises(ValueError, match=_DATED.format("10 days")):
        flow.compare(toy)


def test_a_negative_key_is_refused(toy):
    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0, -1], "W": [95.0, 96.0]}))
    with pytest.raises(ValueError, match=r"and -1 is negative: the first stress period is 0"):
        heads.stats(toy)


# -- ledger 500: every file route refuses a period the run lacks --------

_OFF = (
    r"{caller}: a target is keyed by zero-based stress period and read at its end, and "
    r"the periods {route} holds are 0 \.\. 2; 1 target row\(s\) name none of them: "
    r"'{name}' at time 7"
)


def test_every_file_route_refuses_a_key_past_the_runs_periods(toy):
    """Keys 0 and 7 on a three-period run: compare() gave period 7 a NaN and stats()
    counted one row, with no message -- while the PEST build refused the same key."""

    def values(name, value):
        return pd.DataFrame({"time": [0, 7], name: [value, value]})

    cases = (
        (_head_targets([("W", 22)], values("W", 95.0)), "HeadTargets", "the head file", "W"),
        (mf.LakeStageTargets(locations={"pond": 0}, values=values("pond", 99.0)),
         "LakeStageTargets", "the LAK stage file", "pond"),
        (mf.SfrStageTargets(locations={"g": 1}, values=values("g", 95.0)),
         "SfrStageTargets", "the SFR stage file", "g"),
        (mf.SfrFlowTargets(locations={"a": 0}, values=values("a", -7.0)),
         "SfrFlowTargets", "the SFR budget file", "a"),
    )
    for targets, owner, route, name in cases:
        refusal = _OFF.format(caller=rf"{owner}\.compare\(\.\.\.\)", route=route, name=name)
        with pytest.raises(ValueError, match=refusal):
            targets.compare(toy)
        with pytest.raises(ValueError, match=refusal):
            targets.stats(toy)
    # The periods the run holds still compare.
    assert _head_targets([("W", 22)], pd.DataFrame({"time": [0, 2], "W": [95.0, 96.0]})).stats(toy)["n"].item() == 2


# -- ledgers 501 / 502: a well's layer and name ----------------


def test_a_present_layer_column_must_hold_a_layer_number_on_every_row(toy):
    """['2,3', None, 'L3', 0.7] all became layer 0, in silence."""

    locations = pd.DataFrame({"name": ["A", "B", "C", "D", "E"], "cell": [1, 2, 3, 4, 5],
                              "layer": ["2,3", None, "L3", 0.7, -1]})
    values = pd.DataFrame({"time": [0], "A": [1.0]})
    with pytest.raises(ValueError, match=(
            r"the layer column 'layer' must give every target a zero-based model layer.*"
            r"5 row\(s\) do not: 'A' \('2,3'\); 'B' \(blank\); 'C' \('L3'\); 'D' \(0\.7\); "
            r"'E' \(-1\)")):
        mf.HeadTargets(locations=locations, values=values)
    with pytest.raises(ValueError, match=r"'D' \(0\.7\)"):
        mf.HeadTargets.from_cells(cells=[4], names=["D"], layers=0.7, values=[[1.0]])
    # 1.0 is layer 1; no column at all is layer 0, as documented.
    whole = mf.HeadTargets(locations=pd.DataFrame({"name": ["A"], "cell": [1], "layer": [1.0]}),
                           values=values)
    assert whole.locations_gdf["layer"].tolist() == [1]
    absent = mf.HeadTargets(locations=pd.DataFrame({"name": ["A"], "cell": [1]}), values=values)
    assert absent.locations_gdf["layer"].tolist() == [0]


def test_a_name_repeated_in_one_layer_is_refused_and_in_two_layers_is_two_screens(toy, tmp_path):
    """W1 twice in layer 0 became W1_lay0 and W1_lay0_c5 -- two full-weight observations
    of one series -- and compare() cross-joined it: stats() n = 5 for 2 measurements."""

    values = pd.DataFrame({"time": [0], "W1": [95.0], "W2": [97.0]})
    for cells in ([5, 5, 22], [5, 15, 22]):
        locations = pd.DataFrame({"name": ["W1", "w1 ", "W2"], "cell": cells, "layer": 0})
        with pytest.raises(ValueError, match=(
                rf"a target name may appear once per layer, and 'W1' in layer 0 \(cell "
                rf"{cells[0]}, cell {cells[1]}\) is given more than once")):
            mf.HeadTargets(locations=locations, values=values)
    # Two layers of one well are two observations, as ledger 234 made them.
    two = mf.HeadTargets(locations=pd.DataFrame({"name": ["W1", "W1", "W2"], "cell": [5, 5, 22],
                                                 "layer": [0, 0, 0]}).assign(layer=[0, 1, 0]),
                         values=values)
    assert sorted(two.locations_gdf["layer"]) == [0, 0, 1]


# -- small ones: bound sets observed, compare()'s x/y ----------------------------


def test_a_target_set_read_back_from_the_model_is_observed_and_forecast(toy, tmp_path):
    """cal.observe(model.targets.heads) raised TypeError: the bound helper was not unwrapped."""

    toy.targets.heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0], "W": [95.0]}))
    try:
        cal = toy.pest("bound", workspace=tmp_path / "pest", start_datetime="2024-01-01")
        spec = cal.observe(toy.targets.heads)
        assert spec.targets is toy.targets.heads.targets
        assert cal.forecast(toy.targets.heads).targets is toy.targets.heads.targets
    finally:
        del toy.targets["heads"]


def test_compare_of_point_targets_returns_x_and_y_in_the_grids_frame(toy):
    """A point-located set returned x_x/y_x (its own CRS) and x_y/y_y, and no `x`."""

    centres = toy.vor.gdf_vorPolys.geometry.centroid
    points = gpd.GeoDataFrame({"name": ["A", "B"], "layer": [0, 0]},
                              geometry=[centres.iloc[5], centres.iloc[22]],
                              crs=toy.vor.gdf_vorPolys.crs)
    frame = mf.HeadTargets(locations=points, values=pd.DataFrame({"time": [0], "A": [95.0], "B": [97.0]})).compare(toy)
    assert {"x", "y"} <= set(frame.columns)
    assert not {"x_x", "y_x", "x_y", "y_y"} & set(frame.columns)
    assert frame.set_index("name")["x"].round(6).to_dict() == {
        "A": round(centres.iloc[5].x, 6), "B": round(centres.iloc[22].x, 6)}


# -- ledger 503: lake and reach ids checked; boundnames accepted --------


def test_a_lake_or_reach_id_must_be_a_whole_number_from_0_or_a_boundname():
    """int() truncated 0.7 to lake 0 and 2.6 to reach 2, in silence."""

    for family, feature in ((mf.LakeStageTargets, "lake"), (mf.SfrStageTargets, "reach"),
                            (mf.SfrFlowTargets, "reach")):
        for bad, why in ((0.7, "is not a whole number"), (-1, "is negative"), (None, "is blank"),
                         (" ", "is blank")):
            with pytest.raises(ValueError, match=rf"target 'x' names {feature} .* which {why}"):
                family(locations={"x": bad})
    # A whole number, as a float or as text, is the number; other text is a boundname.
    assert mf.LakeStageTargets(locations={"a": 1.0, "b": "2", "c": "Keevie"}).get()["lake"].tolist() == [1, 2, "Keevie"]


def test_a_lake_or_reach_the_package_does_not_have_is_refused_naming_it(toy):
    """Lake 1 on a one-lake model and reach 9 on four reaches read an empty frame, and
    with values a NaN that stats() dropped -- only the PEST build refused them. MF6
    numbers from ONE, so typing 1 for the first lake is the likely mistake."""

    pond = mf.LakeStageTargets(locations={"p": 0, "past": 1},
                               values=pd.DataFrame({"time": [0], "p": [99.0], "past": [99.0]}))
    for call in (pond.compare, pond.stats, pond.simulated_series):
        with pytest.raises(ValueError, match=(
                r"LakeStageTargets: 1 target\(s\) name a lake the LAK package 'lak' does not "
                r"have \(its lakes are 0 \.\. 0, numbered from 0\): 'past': lake 1 -- one past "
                r"the last: MF6's files number lakes from 1, so the last is 0 here")):
            call(toy)
    gage = mf.SfrFlowTargets(locations={"far": 9})
    with pytest.raises(ValueError, match=r"'far': reach 9\. Such a target read nothing"):
        gage.simulated_series(toy)


def test_a_lake_given_by_its_boundname_reads_that_lake(toy):
    by_number = mf.LakeStageTargets(locations={"p": 0}, values=pd.DataFrame({"time": [0, 1], "p": [99.0, 98.0]}))
    by_name = mf.LakeStageTargets(locations={"p": "POND"}, values=pd.DataFrame({"time": [0, 1], "p": [99.0, 98.0]}))
    pd.testing.assert_frame_equal(by_name.compare(toy), by_number.compare(toy))
    pd.testing.assert_frame_equal(by_name.simulated_series(toy), by_number.simulated_series(toy))
    # Written for MF6 as the boundname, which MF6 resolves in the LAK packagedata.
    assert by_name.to_flopy_obs() == {"lake_stage_targets.csv": [("p", "STAGE", "POND")]}
    with pytest.raises(ValueError, match=r"'q': boundname 'nowhere' names no lake of it"):
        mf.LakeStageTargets(locations={"q": "nowhere"}).simulated_series(toy)


def test_a_reach_boundname_naming_one_reach_reads_it_and_one_naming_two_is_refused(tmp_path):
    grid, props, edge = _grid_and_edges()
    last = len(REACHES) - 1
    names = ["up", "up", "mid", "down"]
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0, names[i]] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]]}, boundnames=True,
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
    )
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props), mf.npf(k=10.0, icelltype=0), mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}), stream,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    run = _run(tmp_path, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)]))
    m = run.model()
    for family in (mf.SfrFlowTargets, mf.SfrStageTargets):
        pd.testing.assert_frame_equal(family(locations={"g": "Mid"}).simulated_series(m),
                                      family(locations={"g": 2}).simulated_series(m))
        with pytest.raises(ValueError, match=r"'g': boundname 'up' names 2 reaches \(0, 1\), and a target reads ONE"):
            family(locations={"g": "up"}).compare(m)


def test_a_lake_or_reach_given_no_value_keeps_its_column_in_every_realization(toy, tmp_path):
    """flows-3: the forward run rebuilt the targets from a snapshot of get() -- the VALUE
    rows -- so a location with no value lost its column in every realization: base
    'time,pond,pond_unmeasured', forward 'time,pond', and pyEMU's instruction read
    failed. The build succeeded."""

    pond = mf.LakeStageTargets(locations={"pond": "pond", "pond_unmeasured": 0},
                               values=pd.DataFrame({"time": [0], "pond": [99.0]}))
    stage = mf.SfrStageTargets(locations={"g": 3, "g_up": 1}, values=pd.DataFrame({"time": [0], "g": [95.0]}))
    flow = mf.SfrFlowTargets(locations={"gage": 3, "gage_up": 1},
                             values=pd.DataFrame({"time": [0], "gage": [-40.0]}))
    cal, _pst = _build(toy, tmp_path / "pest", observe=(pond, stage, flow))
    template = cal.template_workspace
    tables = ("lak_stage_simulated_lake_stage.csv", "sfr_stage_simulated_sfr_stage.csv",
              "sfr_flow_simulated_sfr_flow.csv")
    base = {table: (template / table).read_text().splitlines()[0] for table in tables}
    assert base == {tables[0]: "time,pond,pond_unmeasured", tables[1]: "time,g,g_up",
                    tables[2]: "time,gage,gage_up"}
    for table in tables:
        (template / table).unlink()
    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)
    for table in tables:
        assert (template / table).read_text().splitlines()[0] == base[table]
        pyemu.pst_utils.InstructionFile(str(template / f"{table}.ins")).read_output_file(
            str(template / table))


# -- ledger 504: a reading's own weight ---------------------------------


def test_a_readings_weight_beats_its_locations_which_beats_the_sets_which_beats_one(toy, tmp_path):
    """A `weight` column in a LONG values table was dropped without a word; SFR sets had
    no weight at all; a lake set's weight overwrote every row."""

    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["A", "B"], "cell": [22, 5], "layer": 0, "weight": [2.0, None]}),
        values=pd.DataFrame({"time": [0, 1, 0, 1], "name": ["A", "A", "B", "B"],
                             "head": [95.0, 96.0, 97.0, 98.0], "weight": [0.1, None, 5.0, None]}))
    assert heads.get().set_index(["name", "time"])["weight"].to_dict() == {
        ("A", 0): 0.1, ("A", 1): 2.0, ("B", 0): 5.0, ("B", 1): pytest.approx(float("nan"), nan_ok=True)}
    pond = mf.LakeStageTargets(locations={"pond": 0}, values=pd.DataFrame(
        {"time": [0, 1], "name": ["pond", "pond"], "stage": [99.0, 98.0], "WEIGHT": [0.5, ""]}))
    flow = mf.SfrFlowTargets(locations={"gage": 3}, values=pd.DataFrame(
        {"time": [0, 1], "name": ["gage", "gage"], "flow": [-40.0, -41.0], "weight": [0.01, None]}))
    stage = mf.SfrStageTargets(locations={"g": 1}, values=pd.DataFrame(
        {"time": [0, 1], "stage": [95.0, 95.5], "weight": [3.0, None]}))
    seep = mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["zone"], "cells": [[20, 21, 27]], "weight": [0.2]}),
        values=pd.DataFrame({"time": [0, 1], "name": ["zone", "zone"], "flow": [-3.0, -4.0],
                             "weight": [7.0, None]}))
    cal, pst = _build(toy, tmp_path / "pest", observe=(
        heads, LakeStageObservationSpec(targets=pond, weight=4.0),
        SfrFlowObservationSpec(targets=flow, weight=0.03), stage, seep))
    template = cal.template_workspace
    obs = pst.observation_data
    for table, period, column, weight in (
        ("hds_simulated_heads.csv", 0, "A", 0.1),   # the reading's
        ("hds_simulated_heads.csv", 1, "A", 2.0),   # the location's
        ("hds_simulated_heads.csv", 0, "B", 5.0),
        ("hds_simulated_heads.csv", 1, "B", 1.0),   # neither
        ("lak_stage_simulated_lake_stage.csv", 0, "pond", 0.5),
        ("lak_stage_simulated_lake_stage.csv", 1, "pond", 4.0),   # the set's
        ("sfr_flow_simulated_sfr_flow.csv", 0, "gage", 0.01),
        ("sfr_flow_simulated_sfr_flow.csv", 1, "gage", 0.03),
        ("sfr_stage_simulated_sfr_stage.csv", 0, "g", 3.0),
        ("sfr_stage_simulated_sfr_stage.csv", 1, "g", 1.0),
        ("drn_flow_simulated_drn_flow.csv", 0, "zone", 7.0),
        ("drn_flow_simulated_drn_flow.csv", 1, "zone", 0.2),   # the zone's
    ):
        name = _observation_at(template, table, period, column)
        assert obs.loc[name, "weight"] == pytest.approx(weight), (table, period, column)


def test_a_weight_that_is_not_a_number_from_0_or_sits_in_a_wide_table_is_refused():
    long = pd.DataFrame({"time": [0, 1], "name": ["A", "A"], "head": [1.0, 2.0], "weight": [-1, "x"]})
    locations = pd.DataFrame({"name": ["A"], "cell": [1]})
    with pytest.raises(ValueError, match=r"a reading's weight \(column 'weight'\) must be blank.*"
                                         r"2 are not: 'A' at time 0 weighs -1; 'A' at time 1 weighs 'x'"):
        mf.HeadTargets(locations=locations, values=long)
    with pytest.raises(ValueError, match=r"a WIDE table \(one column per target\) has no room for a "
                                         r"per-reading weight"):
        mf.LakeStageTargets(locations={"pond": 0},
                            values=pd.DataFrame({"time": [0], "pond": [99.0], "weight": [2.0]}))


# -- ledger 505: the head file the OC declares --------------------------


def _renamed_head_run(root, *, chd, head_file):
    """A one-layer run on the toy's grid whose OC writes its heads to ``head_file``."""

    grid, props, edge = _grid_and_edges()
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props), mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.chd(stress_period_data={0: [[(0, int(c)), chd] for c in edge]}),
        mf.oc(head_filerecord=head_file, budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    return _run(root, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])).model()


def test_the_head_writers_read_the_file_the_oc_declares_and_never_a_stale_one(tmp_path):
    """The writers opened `<model>.hds`. With `HEAD FILEOUT renamed.hds` and an old
    `m.hds` left in the run directory, the build took the STALE file's heads as the
    base run and every realization then failed (the template holds no m.hds)."""

    stale = _renamed_head_run(tmp_path / "old", chd=50.0, head_file="m.hds")
    model = _renamed_head_run(tmp_path / "new", chd=102.0, head_file="renamed.hds")
    workspace = Path(model.workspace)
    shutil.copy(Path(stale.workspace) / "m.hds", workspace / "m.hds")
    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0, 1, 2], "W": [99.0, 99.5, 100.0]}))
    expected = heads.simulated_heads(_renamed_head_run(tmp_path / "check", chd=102.0, head_file="m.hds"))

    cal, _pst = _build(model, tmp_path / "pest", observe=(heads,))
    table = cal.template_workspace / "hds_simulated_heads.csv"
    base = pd.read_csv(table)
    np.testing.assert_allclose(base["W"], expected["W"].astype(float), rtol=1e-10)
    table.unlink()
    subprocess.run([sys.executable, "forward_run.py"], cwd=cal.template_workspace, check=True,
                   capture_output=True, text=True)
    pd.testing.assert_frame_equal(pd.read_csv(table), base)

    (workspace / "renamed.hds").unlink()
    with pytest.raises(ValueError, match=(
            r"model 'm' writes its head to 'renamed\.hds' -- what its OC declares \(HEAD "
            r"FILEOUT\).*holds no such file")):
        _build(model, tmp_path / "pest2", observe=(heads,))


# -- ledger 506: every cell PEST reads needs a base value --------------


def test_a_cell_the_base_run_left_empty_is_refused_whatever_its_weight(toy, tmp_path):
    """The build accepted an empty base cell, and PEST++ fails a run on reading an
    empty field at weight 0 exactly as at weight 1
    (docs/reproducers/pest_no_value_fields.py). The first check refused it only under
    a positive weight (review of ledger 506)."""

    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0, 1], "W": [99.0, 99.5]}))
    table = pd.DataFrame({"per": [0, 1, 2], "W": [100.0, np.nan, 100.0]})
    spec = HeadTargetObservationSpec(targets=heads, simulated_values=table)
    with pytest.raises(ValueError, match=(
            r"HeadTargets \(PEST build, prefix 'hds'\): the base run's table "
            r"hds_simulated_heads\.csv holds no number in 1 cell\(s\) PEST's instruction "
            r"file reads: 'W' at 1 \(oname:hds_otype:lst_usecol:w_per:1, weight 1\)")):
        _build(toy, tmp_path / "pest", observe=(spec,))
    # At weight 0 -- the location's, and an unmeasured period -- it is refused too.
    unweighted = heads.locations_gdf.assign(weight=0.0)
    quiet = mf.HeadTargets(locations=unweighted, values=pd.DataFrame({"time": [0, 1], "W": [99.0, 99.5]}))
    for targets, gap in ((quiet, 1), (heads, 2)):
        holed = table.assign(W=[100.0, 100.0, 100.0])
        holed.loc[gap, "W"] = np.nan
        with pytest.raises(ValueError, match=rf"'W' at {gap} \(oname:hds_otype:lst_usecol:w_per:{gap}, weight 0\)"):
            _build(toy, tmp_path / f"pest{gap}",
                   observe=(HeadTargetObservationSpec(targets=targets, simulated_values=holed),))


def test_a_drn_zone_whose_drains_are_off_in_an_unmeasured_period_is_refused(tmp_path):
    """Measured on a real run: a zone whose drains are off in periods 1-2, measured in
    period 0 only, left those cells of the base table empty; the build succeeded and
    every PESTPP-IES realization failed on `error tokenizing output line ('1,0.0,')`."""

    grid, props, edge = _grid_and_edges()
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props), mf.npf(k=10.0, icelltype=0), mf.ic(strt=100.0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        mf.drn(stress_period_data={0: [[(0, c), 96.0, 50.0] for c in (20, 21, 27, 30)],
                                   1: [[(0, 30), 96.0, 50.0]]}),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    run = _run(tmp_path, mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])).model()
    seep = mf.DrnFlowTargets(locations=pd.DataFrame({"name": ["zone", "b"], "cells": [[20, 21, 27], [30]]}),
                             values=pd.DataFrame({"time": [0], "zone": [-3.0]}))
    with pytest.raises(ValueError, match=(
            r"DrnFlowTargets \(PEST build, prefix 'drn_flow'\): the base run's table "
            r"drn_flow_simulated_drn_flow\.csv holds no number in 2 cell\(s\).*"
            r"'zone' at 1 \(oname:drn_flow_otype:lst_usecol:zone_time:1, weight 0\); "
            r"'zone' at 2 ")):
        _build(run, tmp_path / "pest", observe=(seep,))


def test_a_drn_readings_weight_and_its_zones_are_one_weight_column(toy):
    """get() and compare() returned weight_x/weight_y once a values table carried a
    weight, and residuals() raised KeyError naming `weight`."""

    seep = mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["zone"], "cells": [[20, 21, 27]], "weight": [0.2]}),
        values=pd.DataFrame({"time": [0, 1], "name": ["zone", "zone"], "flow": [-3.0, -4.0],
                             "weight": [7.0, None]}))
    for frame in (seep.get(), seep.get(toy), seep.compare(toy), seep.residuals(toy)):
        assert not {"weight_x", "weight_y"} & set(frame.columns)
        assert frame["weight"].tolist() == [7.0, 0.2]
