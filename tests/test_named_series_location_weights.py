"""Lake and SFR targets weigh by their locations, and a name is given once (ledger 609).

``LakeStageTargets``, ``SfrStageTargets`` and ``SfrFlowTargets`` read a locations
table as ``name -> id`` and dropped every other column, so a ``weight`` column -- the
spelling heads, concentrations and DRN/RIV zones take -- was lost in silence, and
those observations reached the control file at 1.0 (a lake known to 0.02 ft,
weight 50, then counted 1/2,500 of what was meant in phi). A name given twice was
resolved in silence: the LAST row of a lake table or of any ``Series``, the FIRST of
an SFR table, and ``{1: 0, "1": 1}`` read lake 1. No family checked a location's
weight: -3 reached the pst as -3.0, ``'x'`` on a zone as 1.0. And three refusals
named spellings that fail -- ``cal.observe(..., weight=)``, which does not exist, and
"the locations' weight column", the one that was dropped.

The build tests build a real control file on the toy run of
``tests/test_pest_handed_targets.py`` (no PEST++ run) and find which observation
holds which cell the way PEST does, through pyEMU's own instruction file.
"""

from __future__ import annotations

import inspect
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.observations.hydrographs import read_readings
from myflopy.modflow.mf6.pest import PestProject
from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec, SfrFlowObservationSpec
from tests.test_pest_handed_targets import _observation_at, run_toy

pyemu = pytest.importorskip("pyemu")

SRC = Path(__file__).resolve().parents[1] / "src" / "myflopy"


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The toy of ``tests/test_pest_handed_targets.py``, run: lake ``pond`` (lake 0), reaches 0-3."""

    return run_toy(tmp_path_factory.mktemp("location_weights"))


def _build(model, workspace: Path, *, observe=()):
    """A calibration of ``k`` against ``observe``, built with the default control settings."""

    cal = model.pest("weights", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in observe:
        cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("weights.pst")
    return cal, pst


def _weights(cal, pst, table, column, periods):
    return [
        float(pst.observation_data.loc[_observation_at(cal.template_workspace, table, period, column), "weight"])
        for period in periods
    ]


def _pond():
    return mf.LakeStageTargets(
        locations=pd.DataFrame({"name": ["pond", "pond2"], "lake": [0, 0], "weight": [50.0, None]}),
        values=pd.DataFrame({
            "time": [0, 1, 2, 0, 1], "name": ["pond", "pond", "pond", "pond2", "pond2"],
            "stage": [99.0, 98.0, 97.0, 99.0, 98.0], "weight": [0.5, None, None, None, None],
        }),
    )


def _stage():
    return mf.SfrStageTargets(
        locations=pd.DataFrame({"name": ["g"], "reach": [1], "weight": [20.0]}),
        values=pd.DataFrame({"time": [0, 1], "g": [95.0, 95.5]}),
    )


def _flow():
    return mf.SfrFlowTargets(
        locations=[{"name": "gage", "reach": 3, "weight": 0.01}],
        values=pd.DataFrame({"time": [0, 1], "gage": [-40.0, -41.0]}),
    )


def test_a_lake_or_reach_locations_weight_reaches_the_pst_below_a_readings_and_above_the_sets(toy, tmp_path):
    """Ledger 504's order, now for lakes and reaches: reading, location, set, 1.0.

    Before, the location's weight was dropped: ``pond`` read ``[0.5, 4, 4]``, ``g``
    ``[1, 1]`` and ``gage`` the set's ``[0.03, 0.03]``.
    """

    cal, pst = _build(toy, tmp_path / "pest", observe=(
        LakeStageObservationSpec(targets=_pond(), weight=4.0), _stage(),
        SfrFlowObservationSpec(targets=_flow(), weight=0.03),
    ))
    lake = "lak_stage_simulated_lake_stage.csv"
    assert _weights(cal, pst, lake, "pond", (0, 1, 2)) == [0.5, 50.0, 50.0]
    assert _weights(cal, pst, lake, "pond2", (0, 1)) == [4.0, 4.0]   # no location weight: the set's
    assert _weights(cal, pst, "sfr_stage_simulated_sfr_stage.csv", "g", (0, 1)) == [20.0, 20.0]
    assert _weights(cal, pst, "sfr_flow_simulated_sfr_flow.csv", "gage", (0, 1)) == pytest.approx([0.01, 0.01])


def _no_suffixed_weight(frame):
    assert not {"weight_x", "weight_y", "_reading_weight"} & set(frame.columns), list(frame.columns)


def test_get_compare_and_residuals_carry_one_weight_column(toy):
    pond = _pond()
    for frame in (pond.get(), pond.compare(toy), pond.residuals(toy)):
        _no_suffixed_weight(frame)
        weights = frame.set_index(["name", "time"])["weight"]
        assert weights[("pond", 0)] == 0.5          # the reading's
        assert weights[("pond", 1)] == 50.0         # the location's
        assert np.isnan(weights[("pond2", 0)])      # neither: the set's, else 1.0, in the build
    for targets, value in ((_stage(), 20.0), (_flow(), 0.01)):
        for frame in (targets.get(), targets.compare(toy), targets.residuals(toy)):
            _no_suffixed_weight(frame)
            assert frame["weight"].tolist() == [value, value]
    # A set given no values: its definitions, and every simulated row, carry the location's.
    for targets, value in (
        (mf.LakeStageTargets(locations=pd.DataFrame({"name": ["pond"], "lake": [0], "weight": [50.0]})), 50.0),
        (mf.SfrStageTargets(locations=pd.DataFrame({"name": ["g"], "reach": [1], "weight": [20.0]})), 20.0),
        (mf.SfrFlowTargets(locations=[{"name": "gage", "reach": 3, "weight": 0.01}]), 0.01),
    ):
        assert targets.get()["weight"].tolist() == [value]
        for frame in (targets.compare(toy), targets.residuals(toy)):
            _no_suffixed_weight(frame)
            assert len(frame) == 3 and frame["weight"].eq(value).all()


@pytest.mark.parametrize("family, id_column", [
    (mf.LakeStageTargets, "lake"), (mf.SfrStageTargets, "reach"), (mf.SfrFlowTargets, "reach"),
])
def test_weight_column_names_the_column_and_none_reads_none(family, id_column):
    table = pd.DataFrame({"name": ["a"], id_column: [0], "w": [7.0], "weight": [3.0]})
    assert family(locations=table, weight_column="w").get()["weight"].tolist() == [7.0]
    assert family(locations=table).get()["weight"].tolist() == [3.0]
    assert np.isnan(family(locations=table, weight_column=None).get()["weight"].iloc[0])
    assert np.isnan(family(locations={"a": 0}).get()["weight"].iloc[0])
    assert np.isnan(family(locations=pd.Series({"a": 0})).get()["weight"].iloc[0])


@pytest.mark.parametrize("family, id_column", [
    (mf.LakeStageTargets, "lake"), (mf.SfrStageTargets, "reach"), (mf.SfrFlowTargets, "reach"),
])
@pytest.mark.parametrize("shape", ["frame", "records", "series", "dict", "identical"])
def test_a_name_given_twice_is_refused_naming_both_rows(family, id_column, shape):
    """Last-wins for a lake table and any Series, first-wins for an SFR table, in silence."""

    rows = [{"name": "a", id_column: 0}, {"name": "a", id_column: 1}]
    locations = {
        "frame": pd.DataFrame(rows),
        "records": rows,
        "series": pd.Series([0, 1], index=["a", "a"]),
        "dict": {1: 0, "1": 1},
        "identical": pd.DataFrame([rows[0], rows[0]]),
    }[shape]
    name = "'1'" if shape == "dict" else "'a'"
    second = 0 if shape == "identical" else 1
    with pytest.raises(ValueError, match=(
        rf"{family.__name__} locations: a target name may be given once, and {name} "
        rf"\(row 0: {id_column} 0, row 1: {id_column} {second}\) is given more than once"
    )):
        family(locations=locations)


def test_a_dataframes_rows_are_named_by_its_own_labels():
    with pytest.raises(ValueError, match=r"'a' \(row 'x': lake 0, row 'z': lake 2\)"):
        mf.LakeStageTargets(locations=pd.DataFrame(
            {"name": ["a", "b", "a"], "lake": [0, 1, 2]}, index=["x", "y", "z"]))


def test_from_records_refuses_a_name_on_two_lakes_and_keeps_each_records_weight():
    records = [
        {"name": "a", "lake": 0, "time": 0, "stage": 1.0},
        {"name": "a", "lake": 1, "time": 1, "stage": 2.0},
    ]
    with pytest.raises(ValueError, match=r"'a' \(row 0: lake 0, row 1: lake 1\) is given more than once"):
        mf.LakeStageTargets.from_records(records)
    kept = mf.LakeStageTargets.from_records([
        {"name": "a", "lake": 0, "time": 0, "stage": 1.0, "weight": 50.0},
        {"name": "a", "lake": 0, "time": 1, "stage": 2.0, "weight": None},
    ])
    frame = kept.get()
    assert frame["lake"].tolist() == [0, 0]
    assert frame["weight"].iloc[0] == 50.0 and np.isnan(frame["weight"].iloc[1])
    with pytest.raises(ValueError, match=r"a reading's weight \(column 'weight'\) must be blank"):
        mf.LakeStageTargets.from_records([{"name": "a", "lake": 0, "time": 0, "stage": 1.0, "weight": -2}])


_LOCATION_WEIGHT = (
    r"a location's weight \(column 'weight'\) must be blank -- none of its own -- or a finite "
    r"number from 0, and 2 are not: 'a' weighs -1; 'b' weighs 'x'"
)


@pytest.mark.parametrize("family", [
    mf.HeadTargets, mf.ConcTargets, mf.LakeStageTargets, mf.SfrStageTargets, mf.SfrFlowTargets,
    mf.DrnFlowTargets, mf.RivFlowTargets,
])
def test_a_location_weight_that_is_not_a_number_from_0_is_refused_on_every_family(family):
    """A head location weight of -3 reached the pst as -3.0, and a zone's 'x' as 1.0."""

    columns = {
        mf.HeadTargets: {"cell": [1, 2]}, mf.ConcTargets: {"cell": [1, 2]},
        mf.LakeStageTargets: {"lake": [0, 0]}, mf.SfrStageTargets: {"reach": [0, 1]},
        mf.SfrFlowTargets: {"reach": [0, 1]},
        mf.DrnFlowTargets: {"cells": [[1], [2]]}, mf.RivFlowTargets: {"cells": [[1], [2]]},
    }[family]
    locations = pd.DataFrame({"name": ["a", "b"], **columns, "weight": [-1, "x"]})
    values = pd.DataFrame({"time": [0], "a": [-1.0], "b": [-1.0]})
    with pytest.raises(ValueError, match=rf"{family.__name__} locations: {_LOCATION_WEIGHT}"):
        family(locations=locations, values=values)
    # A blank and a number from 0 pass, read as floats.
    locations["weight"] = ["", "2.5"]
    got = family(locations=locations, values=values).get().sort_values("name")
    assert got["weight"].fillna(-1).tolist() == [-1.0, 2.5]


def test_from_cells_weights_pass_the_same_rule():
    with pytest.raises(ValueError, match=r"HeadTargets locations: a location's weight .* 'OBS_001' weighs -1.0"):
        mf.HeadTargets.from_cells(cells=[1], values=[[95.0]], weights=-1.0)


def test_no_refusal_names_a_weight_spelling_that_fails():
    """``cal.observe(..., weight=)`` raised TypeError, and the WIDE refusal named the dropped column."""

    assert "weight" not in inspect.signature(PestProject.observe).parameters
    offenders = [
        str(path.relative_to(SRC)) for path in SRC.rglob("*.py")
        if re.search(r"observe\(\.\.\., weight=", path.read_text())
    ]
    assert offenders == []
    wide = pd.DataFrame({"time": [0], "pond": [99.0], "weight": [2.0]})
    with pytest.raises(ValueError, match=r"in the weight column of its locations table \(weight_column=, default 'weight'\)"):
        mf.LakeStageTargets(locations={"pond": 0}, values=wide)
    with pytest.raises(ValueError, match=r"through the weight column of its locations table "
                                         r"\(weight_column=, default 'weight'\)"):
        read_readings(pd.DataFrame({"date": ["2024-01-01"], "pond": [1.0], "weight": [2.0]}),
                      value_column="stage", owner="LakeStageTargets")
    # Doing what the WIDE refusal says carries the weight.
    followed = mf.LakeStageTargets(
        locations=pd.DataFrame({"name": ["pond"], "lake": [0], "weight": [2.0]}), values=wide.drop(columns="weight"))
    assert followed.get()["weight"].tolist() == [2.0]
