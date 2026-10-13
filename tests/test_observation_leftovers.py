"""The observation-side leftovers of the PEST-preparation lanes (ledgers NEW-L2OBS-*).

Most tests run on the real MF6 toy of ``tests/test_pest_handed_targets.py`` (a 6 x 6
Voronoi grid, one layer, three 10-day periods of three steps, a four-reach stream, a
three-cell DRN zone and a two-cell lake). What each pins:

* a boolean time key or layer was read as period/layer 1 (``pd.to_numeric(True)``);
* an OC that is readable but declares no ``HEAD FILEOUT`` made every target reader
  open ``<model>.hds`` -- which MF6 never writes then, so any such file is another
  run's (ledger 576);
* a repeated name with an unreadable cell raised a bare ``int()`` error, and a
  fractional cell was truncated;
* a table ``get()`` returned crashed in pandas' ``melt`` (ledger 577);
* the IES review's ``conflict()``/``obs_residuals()`` carried no date (ledger 578);
* a screened well reached the PEST snapshot, the review and the heads calibration
  plot as layer ``NaN``/0 (ledger 580);
* two head sets copied the forward-run writer twice (ledger 583);
* ``LakeStageTargets.from_series``/``from_records`` dropped ``package=``;
* a model with no run behind it was read from a value table (ledger 581).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from tests.test_pest_handed_targets import (
    NPER,
    NSTP,
    PERLEN,
    _build,
    _disv,
    _grid_and_edges,
    _head_targets,
    _run,
    run_toy,
)

pyemu = pytest.importorskip("pyemu")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The toy model of ``test_pest_handed_targets``, run."""

    return run_toy(tmp_path_factory.mktemp("observation_leftovers"))


# -- a boolean is not a period or a layer ---------------------------------


def test_a_boolean_time_key_is_refused_not_read_as_period_one(toy):
    """``pd.to_numeric(True)`` is 1: a flag column keyed the target by period 1."""

    for key in ([True], pd.Series([True], dtype=object), pd.Series([True], dtype="boolean")):
        heads = _head_targets([("W", 22)], pd.DataFrame({"time": key, "W": [95.0]}))
        with pytest.raises(ValueError, match=r"reads `time` as a zero-based stress period "
                                             r"NUMBER, and True is True/False"):
            heads.compare(toy)
    lake = mf.LakeStageTargets(locations={"pond": 0},
                               values=pd.DataFrame({"time": [False, True], "pond": [99.0, 99.5]}))
    with pytest.raises(ValueError, match=r"and False, True are True/False"):
        lake.compare(toy)


def test_a_boolean_layer_is_refused_not_read_as_layer_one():
    locations = pd.DataFrame({"name": ["A", "B"], "cell": [1, 2], "layer": [True, False]})
    with pytest.raises(ValueError, match=(
            r"the layer column 'layer' must give every target a zero-based model layer, and "
            r"2 row\(s\) hold True/False instead: 'A' \(True\); 'B' \(False\)")):
        mf.HeadTargets(locations=locations, values=pd.DataFrame({"time": [0], "A": [1.0]}))
    mixed = locations.astype({"layer": object}).assign(layer=[np.bool_(True), 2])
    with pytest.raises(ValueError, match=r"1 row\(s\) hold True/False instead: 'A' \(True\)"):
        mf.HeadTargets(locations=mixed, values=pd.DataFrame({"time": [0], "A": [1.0]}))


# -- ledger 576: an OC that declares no head file -----------------


@pytest.fixture(scope="module")
def headless(tmp_path_factory, toy):
    """The toy's flow model run with an OC saving the budget only, plus a STALE ``m.hds``.

    MF6 refuses ``SAVE HEAD`` without ``HEAD FILEOUT``, so a run whose OC declares
    none saved no head; the ``m.hds`` copied in is the toy's, the file a reused
    run directory would still hold.
    """

    grid, props, edge = _grid_and_edges()
    model = mf.gwf("m", grid=grid, packages=[
        _disv(grid, props), mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        mf.oc(budget_filerecord="m.cbc", saverecord=[("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    run = _run(tmp_path_factory.mktemp("headless"), simulation)
    shutil.copy(Path(toy.workspace) / "m.hds", Path(run.workspace) / "m.hds")
    return run.model()


_NO_HEAD_FILE = (
    r"model 'm''s output control 'm.oc' declares no HEAD FILEOUT, so the run saved no head "
    r"-- MF6 refuses SAVE HEAD without one -- and a m.hds in the run directory is not this "
    r"run's \(ledger 576\)"
)


def test_targets_on_a_run_whose_oc_declares_no_head_file_are_refused_not_read_from_a_stale_hds(
    headless, tmp_path
):
    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0, 1], "W": [95.0, 96.0]}))
    for call in (lambda: heads.compare(headless), lambda: heads.simulated_heads(headless),
                 lambda: heads.stats(headless)):
        with pytest.raises(ValueError, match="HeadTargets: " + _NO_HEAD_FILE):
            call()
    with pytest.raises(ValueError, match=r"HeadTargets \(PEST build, prefix 'hds'\): "
                                         + _NO_HEAD_FILE):
        _build(headless, tmp_path / "pest", observe=(heads,))


def test_a_transport_oc_declaring_no_concentration_file_is_refused(tmp_path):
    """The rule reads the OC file, so a GWT model's is checked the same way."""

    from myflopy.modflow.mf6.output_files import _declared_field_output

    (tmp_path / "t.nam").write_text("BEGIN packages\n  OC6  t.oc  oc\nEND packages\n")
    (tmp_path / "t.oc").write_text("BEGIN options\n  BUDGET FILEOUT t.cbc\nEND options\n")
    (tmp_path / "t.ucn").write_bytes(b"stale")
    with pytest.raises(ValueError, match=(
            r"ConcTargets: model 't''s output control 't.oc' declares no CONCENTRATION FILEOUT"
            r".*mf.oc\(concentration_filerecord=\"t.ucn\"")):
        _declared_field_output(tmp_path, "t", "concentration", caller="ConcTargets")
    (tmp_path / "t.oc").write_text("BEGIN options\n  CONCENTRATION FILEOUT other.ucn\nEND options\n")
    assert _declared_field_output(tmp_path, "t", "concentration", caller="x") == tmp_path / "other.ucn"
    # An OC that cannot be read (no name file) still reads <model>.<suffix>.
    assert _declared_field_output(tmp_path, "u", "concentration", caller="x") == tmp_path / "u.ucn"


# -- a cell id is checked before the repeated-name message prints it -------


def test_a_repeated_name_with_an_unreadable_cell_names_the_cell_not_a_bare_int_error():
    values = pd.DataFrame({"time": [0], "W": [1.0]})
    for cells, shown in ((["x", "y"], r"'W' \('x'\); 'w' \('y'\)"), ([1.5, 2], r"'W' \(1\.5\)"),
                         ([True, 3], r"'W' \(True\)")):
        locations = pd.DataFrame({"name": ["W", "w"], "cell": cells, "layer": [0, 0]})
        with pytest.raises(ValueError, match=r"Head target cell/node values must be numeric: each "
                                             r"row names its cell by a whole number, and \d row\(s\) "
                                             rf"do not: {shown}"):
            mf.HeadTargets(locations=locations, values=values)
    # A whole number written as text or a float is still a cell.
    read = mf.HeadTargets(locations=pd.DataFrame({"name": ["A", "B"], "cell": ["3", 4.0]}),
                          values=values.rename(columns={"W": "A"}))
    assert read.locations_gdf["cell"].tolist() == [3, 4]


# -- ledger 577: a table get() returned reads back -------------------


def test_a_long_table_whose_value_column_is_head_target_reads_back(toy):
    """``get()``'s own table crashed in pandas' melt (``value_name (head_target) cannot
    match an element in the DataFrame columns``)."""

    locations = pd.DataFrame({"name": ["W", "V"], "cell": [22, 15]})
    heads = mf.HeadTargets(locations=locations,
                           values=pd.DataFrame({"time": [0, 1], "W": [95.0, 96.0], "V": [97.0, None]}))
    again = mf.HeadTargets(locations=locations, values=heads.get())
    columns = ["time", "name", "head_target", "cell", "layer"]
    pd.testing.assert_frame_equal(again.get()[columns], heads.get()[columns])
    pd.testing.assert_frame_equal(again.compare(toy), heads.compare(toy), check_like=True)
    # One series under that name reads as one series.
    single = mf.HeadTargets(locations=locations.iloc[:1],
                            values=pd.DataFrame({"time": [0], "head_target": [95.0]}))
    assert single.get()[["name", "head_target"]].values.tolist() == [["W", 95.0]]
    # The lake family's own value column, ``stage``, too.
    pond = mf.LakeStageTargets(locations={"pond": 0},
                               values=pd.DataFrame({"time": [0, 1], "pond": [99.0, 99.5]}))
    assert mf.LakeStageTargets(locations={"pond": 0}, values=pond.get()).get().equals(pond.get())


def test_a_wide_table_with_a_head_target_column_beside_targets_is_refused_naming_it():
    with pytest.raises(ValueError, match=(
            r"HeadTargets values: a WIDE table \(one column per target\) cannot carry a column "
            r"named 'head_target'")):
        mf.HeadTargets(locations=pd.DataFrame({"name": ["W"], "cell": [3]}),
                       values=pd.DataFrame({"time": [0], "head_target": [1.0], "W": [2.0]}))


# -- ledger 344: from_series / from_records name the LAK package ------------------


def test_lake_from_series_and_from_records_carry_package(toy):
    """Neither forwarded ``package=``, so a set they built could not name one of
    several LAK packages (ledger 344, 352)."""

    name = "LAK"  # the toy's one LAK package, named in another case than FloPy holds it
    built = (
        mf.LakeStageTargets.from_series(lake="pond", lake_id=0, values=[99.0], package=name),
        mf.LakeStageTargets.from_series(lake="pond", lake_id=0, package=name,
                                        values=pd.Series([99.0], index=[0])),
        mf.LakeStageTargets.from_records([{"name": "pond", "lake": 0, "time": 0, "stage": 99.0}],
                                         package=name),
    )
    for targets in built:
        assert targets.package == name
        assert targets.compare(toy)["sim_stage"].notna().all()
    wrong = mf.LakeStageTargets.from_series(lake="pond", lake_id=0, values=[99.0], package="nope")
    with pytest.raises(ValueError, match=r"package='nope' names no LAK package"):
        wrong.compare(toy)


# -- ledger 581: one head reader --------------------------------------


def test_a_model_with_no_run_behind_it_is_refused_not_read_from_a_value_table():
    """A double holding an ``all_heads`` table was read through a second derivation
    kept for doubles alone; the targets read a run's head file, or refuse."""

    from types import SimpleNamespace

    from tests.test_mf6_pest import _two_cell_vor_clockwise

    index = pd.MultiIndex.from_tuples([((0, 0), 0, 0), ((0, 0), 0, 1)],
                                      names=["kstpkper", "layer", "cell"])
    double = SimpleNamespace(vor=_two_cell_vor_clockwise(),
                             all_heads=pd.DataFrame({"elev": [10.0, 9.5]}, index=index))
    heads = mf.HeadTargets(locations={"A": 0}, values={"time": [0], "A": [10.0]})
    with pytest.raises(ValueError, match=(
            r"HeadTargets reads its simulated values from the head file a RUN model's output "
            r"control declares, .* and this model has no workspace \(ledger 581\)")):
        heads.compare(double)
