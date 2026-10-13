"""A target set can be a one-sided BOUND, and the review scores it as PEST++ does.

Ledger 564. ``inequality="at_most"`` / ``"at_least"`` on a target set
(heads, concentration, lake and reach stage, reach flow, DRN and RIV zones) puts each
MEASURED observation in PESTPP-IES's inequality group -- ``l_`` / ``g_`` before the
group it would otherwise be in -- which scores 0 wherever the bound holds. Ledger
358's refusal of a USER label or prefix starting ``l_``/``less``/``g_``/``greater``
stays: ``inequality=`` is the one spelling, and the refusal names it.

The review (ies-6): ``phi_contributions`` reads PESTPP-IES's own ``phi.group.csv`` and
falls back to scoring through ONE helper using pyEMU's tag rule; ``conflict()`` is
one-sided for a bound; ``obs_residuals`` scores a bound 0 where it holds; the panels
draw a bound as a triangle. Before, both review tools scored every group two-sided,
so a bound satisfied by every realization read as the dominant misfit and a 100%
prior-data conflict.

The PEST tests build on ``test_pest_handed_targets.py``'s toy (6 x 6 Voronoi, three
10-day periods, heads ~100-102 ft) and run PESTPP-IES's prior on it.
"""

from __future__ import annotations

import shutil

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.ies import IesResults
from tests.test_pest_handed_targets import DRN_CELLS, _build, run_toy

pyemu = pytest.importorskip("pyemu")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return run_toy(tmp_path_factory.mktemp("bound_targets"))


def _heads(values, *, groups=(None, None), inequality=None, names=("W1", "W2"), cells=(1, 2)):
    locations = pd.DataFrame({"name": list(names), "cell": list(cells), "layer": 0, "group": list(groups)})
    return mf.HeadTargets(locations=locations, values=pd.DataFrame(values), inequality=inequality)


# -- the spelling --------------------------------------------------------------------


@pytest.mark.parametrize("value", ["less", "AT_MOST", "upper", 1, True])
@pytest.mark.parametrize("family", ["head", "conc", "lake", "sfr_stage", "sfr_flow", "drn", "riv"])
def test_an_inequality_that_is_not_one_of_the_two_is_refused_naming_them(family, value):
    kwargs = {
        "head": (mf.HeadTargets, {"locations": {"a": 1}, "values": {"time": [0], "a": [1.0]}}),
        "conc": (mf.ConcTargets, {"locations": {"a": 1}, "values": {"time": [0], "a": [1.0]}}),
        "lake": (mf.LakeStageTargets, {"locations": {"a": 0}}),
        "sfr_stage": (mf.SfrStageTargets, {"locations": {"a": 0}}),
        "sfr_flow": (mf.SfrFlowTargets, {"locations": {"a": 0}}),
        "drn": (mf.DrnFlowTargets, {"locations": pd.DataFrame({"name": ["z"], "cells": [[1]]})}),
        "riv": (mf.RivFlowTargets, {"locations": pd.DataFrame({"name": ["z"], "cells": [[1]]})}),
    }
    cls, base = kwargs[family]
    assert cls(**base, inequality="at_most").inequality == "at_most"
    assert cls(**base, inequality="at_least").inequality == "at_least"
    with pytest.raises(ValueError, match=r"is not a bound.*'at_most'.*'at_least'"):
        cls(**base, inequality=value)


def test_every_constructor_forwards_inequality():
    cells = mf.HeadTargets.from_cells(cells=[1], values=[[1.0]], inequality="at_least")
    records = mf.HeadTargets.from_records(
        [{"name": "a", "cell": 1, "time": 0, "head": 1.0}], inequality="at_most")
    series = mf.LakeStageTargets.from_series(lake="p", lake_id=0, values=[1.0], inequality="at_least")
    lake_records = mf.LakeStageTargets.from_records(
        [{"name": "p", "lake": 0, "time": 0, "stage": 1.0}], inequality="at_most")
    conc = mf.ConcTargets.from_cells(cells=[1], values=[[1.0]], inequality="at_most")
    assert [t.inequality for t in (cells, records, series, lake_records, conc)] == [
        "at_least", "at_most", "at_least", "at_most", "at_most"]


# -- compare(): PEST++'s residual, and the raw one beside it --------------------------


def test_compare_scores_a_satisfied_bound_zero_and_a_violated_one_its_misfit(toy):
    """The toy's heads are ~100-102 ft: at most 1000 holds, at most 1 is violated."""

    values = {"time": [0, 1], "W1": [1000.0, 1.0], "W2": [1000.0, None]}
    plain = _heads(values).compare(toy).set_index(["name", "per"])
    bound = _heads(values, inequality="at_most").compare(toy).set_index(["name", "per"])
    assert "raw_residual" not in plain.columns
    assert bound["raw_residual"].equals(plain["residual"])
    assert bound.loc[("W1", 0), "residual"] == 0.0
    assert bound.loc[("W2", 0), "residual"] == 0.0
    assert bound.loc[("W1", 1), "residual"] == pytest.approx(plain.loc[("W1", 1), "residual"])
    assert bound.loc[("W1", 1), "residual"] > 90.0
    assert np.isnan(bound.loc[("W2", 1), "residual"])  # blank stays blank
    stats = _heads(values, inequality="at_most").stats(toy).iloc[0]
    assert stats["n"] == 3
    assert stats["rmse"] == pytest.approx(np.sqrt(plain.loc[("W1", 1), "residual"] ** 2 / 3))
    assert "raw_residual" in _heads(values, inequality="at_most").residuals(toy).columns
    # at least 1000 is violated everywhere, at least 1 holds.
    above = _heads(values, inequality="at_least").compare(toy).set_index(["name", "per"])
    assert above.loc[("W1", 0), "residual"] < -890.0
    assert above.loc[("W1", 1), "residual"] == 0.0


def test_a_flow_bound_is_on_the_number_as_written_in_mf6s_sign(toy):
    """'At least 1e9 leaves the reach' is -1e9 at_most: the reach carries far less, so
    it is violated; 'at most 1e9 leaves' is -1e9 at_least and holds (ledger 350)."""

    flows = {"time": [0.0], "a": [-1.0e9]}
    leaves_at_least = mf.SfrFlowTargets(locations={"a": 0}, values=flows, inequality="at_most")
    leaves_at_most = mf.SfrFlowTargets(locations={"a": 0}, values=flows, inequality="at_least")
    violated = leaves_at_least.compare(toy).iloc[0]
    assert violated["residual"] == pytest.approx(violated["raw_residual"]) and violated["residual"] > 0
    assert leaves_at_most.compare(toy).iloc[0]["residual"] == 0.0
    with pytest.raises(ValueError, match="Enter a gauged streamflow Q as -Q"):
        mf.SfrFlowTargets(locations={"a": 0}, values={"time": [0.0], "a": [2.0]}, inequality="at_most")
    zone = mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["zone"], "cells": [DRN_CELLS]}),
        values=pd.DataFrame({"time": [2.0], "zone": [0.0]}), inequality="at_most",
    )
    scored = zone.compare(toy).iloc[0]
    # The drains discharge (negative), so "at most 0" holds.
    assert scored["raw_residual"] < 0 and scored["residual"] == 0.0


# -- the build ---------------------------------------------------------------------


def _groups(pst) -> dict[str, list[str]]:
    obs = pst.observation_data
    return {group: sorted(obs.index[obs["obgnme"] == group]) for group in sorted(set(obs["obgnme"]))}


def test_a_bound_sets_measured_observations_are_in_pestpp_ies_inequality_groups(toy, tmp_path):
    cap = _heads({"time": [0, 1], "W1": [1000.0, 1.0], "W2": [1000.0, None]},
                 groups=("North", None), inequality="at_most")
    floor = _heads({"time": [2], "W3": [50.0]}, names=("W3",), cells=(34,), groups=(None,),
                   inequality="at_least")
    stage = mf.LakeStageTargets(locations={"pond": 0}, values={"time": [1.0], "pond": [200.0]},
                                inequality="at_most")
    cal = toy.pest("bounds", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(cap)
    cal.observe(floor, prefix="floor")
    cal.observe(stage)
    pst = cal.build("bounds.pst", noptmax=0)
    groups = _groups(pst)
    usecol = "oname:hds_otype:lst_usecol:{}"
    assert groups["l_north"] == [usecol.format("w1_per:0"), usecol.format("w1_per:1")]
    assert groups["l_hds"] == [usecol.format("w2_per:0")]
    assert groups["g_floor"] == ["oname:floor_otype:lst_usecol:w3_per:2"]
    assert groups["l_lak_stage"] == ["oname:lak_stage_otype:lst_usecol:pond_time:1"]
    # Unmeasured rows are in the set's `_unmeasured`, never in a bound group.
    obs = pst.observation_data
    for name in ("hds_unmeasured", "floor_unmeasured", "lak_stage_unmeasured"):
        assert (obs.loc[obs["obgnme"] == name, "weight"] == 0).all()
    # pyEMU's own tag rule reads exactly these as the inequalities.
    assert sorted(pst.less_than_obs_constraints) == sorted(
        groups["l_north"] + groups["l_hds"] + groups["l_lak_stage"])
    assert list(pst.greater_than_obs_constraints) == groups["g_floor"]
    # The build records it for the review.
    sets = {entry["prefix"]: entry for entry in cal._prepared_observations}
    assert sets["hds"]["metadata"]["inequality"] == "at_most"
    assert sets["floor"]["metadata"]["inequality"] == "at_least"
    snapshot = pd.read_csv(cal.template_workspace / "hds_target_values.csv")
    assert set(snapshot["inequality"]) == {"at_most"}


@pytest.mark.parametrize("label", ["less_x", "l_upper", "Greater Seattle"])
def test_a_user_label_that_reads_as_an_inequality_is_still_refused_naming_inequality(toy, tmp_path, label):
    heads = _heads({"time": [0], "W1": [90.0], "W2": [91.0]}, groups=(label, None))
    with pytest.raises(ValueError, match=r"INEQUALITY.*inequality='at_most' or 'at_least'"):
        _build(toy, tmp_path / "pest", observe=(heads,))


def test_a_prefix_that_reads_as_an_inequality_is_still_refused_naming_inequality(toy, tmp_path):
    cal = toy.pest("ineq", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(_heads({"time": [0], "W1": [90.0], "W2": [91.0]}), prefix="l_wells")
    with pytest.raises(ValueError, match=r"INEQUALITY.*inequality='at_most' or 'at_least'"):
        cal.build("ineq.pst", noptmax=0)


def test_a_forecast_of_a_bound_set_is_refused(toy, tmp_path):
    cal = toy.pest("fore", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    with pytest.raises(ValueError, match=r"forecast\(\).*BOUND \(inequality='at_least'\)"):
        cal.forecast(_heads({"time": [0], "W1": [90.0], "W2": [91.0]}, inequality="at_least"))
    assert cal._forecast_specs == []


# -- PESTPP-IES and the review --------------------------------------------------------


@pytest.fixture(scope="module")
def bound_run(toy, tmp_path_factory):
    """Two heads bounds through PESTPP-IES's prior: ``cap`` (at most 1000 ft, held by
    every realization) and ``low`` (at most 1 ft, violated by every one), beside an
    ordinary set."""

    workspace = tmp_path_factory.mktemp("bound_ies")
    cal = toy.pest("bounded", workspace=workspace / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(_heads({"time": [0, 1], "W1": [1000.0, 1000.0], "W2": [1000.0, None]},
                       inequality="at_most"), prefix="cap")
    cal.observe(_heads({"time": [2], "W3": [1.0]}, names=("W3",), cells=(34,), groups=(None,),
                       inequality="at_most"), prefix="low")
    cal.observe(_heads({"time": [0], "W4": [95.0]}, names=("W4",), cells=(23,), groups=(None,)))
    cal.build("bounded.pst", noptmax=0)
    return cal, cal.prior(reals=3)


def _by_hand(cal, prefix: str, column: str, per: int, target: float) -> float:
    simulated = pd.read_csv(cal.template_workspace / f"{prefix}_simulated_heads.csv").set_index("per")
    return float((simulated.loc[per, column] - target) ** 2)


def test_pestpp_ies_scores_a_held_bound_zero_and_a_violated_one_its_full_misfit(bound_run):
    """Read from PESTPP-IES's own file: the measurement this lane rests on."""

    cal, ies = bound_run
    table = pd.read_csv(ies.workspace / "bounded.phi.group.csv")
    base = table.loc[table["par_realization"].astype(str) == "base"].iloc[0]
    assert base["l_cap"] == 0.0
    assert base["l_low"] == pytest.approx(_by_hand(cal, "low", "W3", 2, 1.0), rel=1e-4)
    assert base["hds"] == pytest.approx(_by_hand(cal, "hds", "W4", 0, 95.0), rel=1e-4)
    # The review reads that file ...
    contributions = ies.phi_contributions()
    assert contributions["l_cap"] == 0.0
    assert contributions["l_low"] == pytest.approx(base["l_low"])
    assert list(contributions.index) == ["l_low", "hds", "l_cap"]


def test_phi_contributions_scores_a_bound_one_sided_without_the_file(bound_run, tmp_path):
    """A run whose phi-by-group file is gone is scored through the review's one helper,
    PEST++'s rule: the held bound 0 -- it read 2,423,185, two-sided -- and the
    others what PESTPP-IES wrote, to its six figures."""

    _cal, ies = bound_run
    copy = tmp_path / "copy"
    shutil.copytree(ies.workspace, copy)
    (copy / "bounded.phi.group.csv").unlink()
    scored = IesResults(copy, case_name="bounded").phi_contributions()
    written = ies.phi_contributions()
    assert scored["l_cap"] == 0.0
    assert scored["l_low"] == pytest.approx(written["l_low"], rel=1e-5)
    assert scored["hds"] == pytest.approx(written["hds"], rel=1e-5)


def test_phi_contributions_is_what_pestpp_ies_wrote(bound_run, tmp_path):
    """The file is the authority when it covers every weighted group: a number written
    there is the number reported, at the posterior iteration and the realization asked."""

    _cal, ies = bound_run
    copy = tmp_path / "copy"
    shutil.copytree(ies.workspace, copy)
    path = copy / "bounded.phi.group.csv"
    table = pd.read_csv(path)
    table.loc[table["par_realization"].astype(str) == "base", "l_low"] = 12345.0
    table.to_csv(path, index=False)
    review = IesResults(copy, case_name="bounded")
    assert review.phi_contributions()["l_low"] == 12345.0
    other = str(table.loc[table["par_realization"].astype(str) != "base", "par_realization"].iloc[0])
    assert review.phi_contributions(realization=other)["l_low"] != 12345.0
    # A group the file does not hold (the control file edited since): scored here.
    table.drop(columns=["hds"]).to_csv(path, index=False)
    assert IesResults(copy, case_name="bounded").phi_contributions()["l_low"] != 12345.0


def test_conflict_is_one_sided_for_a_bound(bound_run):
    _cal, ies = bound_run
    table = ies.conflict()
    cap = table[table["obgnme"] == "l_cap"]
    low = table[table["obgnme"] == "l_low"]
    assert len(cap) == 3 and set(cap["bound"]) == {"at_most"} and not cap["in_conflict"].any()
    # Every prior head is above 1 ft: the whole band on the wrong side.
    assert low["in_conflict"].all()
    assert table.loc[table["obgnme"] == "hds", "bound"].isna().all()


def test_obs_residuals_scores_a_held_bound_zero(bound_run):
    _cal, ies = bound_run
    frame = ies.obs_residuals().set_index(["prefix", "location"])
    assert frame.loc[("cap", "w1"), "residual"] == 0.0
    assert frame.loc[("cap", "w1"), "raw_residual"] < -890.0
    assert frame.loc[("low", "w3"), "residual"] == pytest.approx(frame.loc[("low", "w3"), "raw_residual"])
    assert frame.loc[("hds", "w4"), "residual"] == pytest.approx(frame.loc[("hds", "w4"), "raw_residual"])


def test_the_review_draws_a_bound_as_a_triangle_and_a_measurement_as_a_circle(bound_run):
    _cal, ies = bound_run
    for figure in (ies.plot_vs_obs(), ies.plot_prior_vs_obs(show_conflict=False)):
        markers = {(trace.name, trace.marker.symbol) for trace in figure.data if trace.mode == "markers"}
        assert markers == {("bound: at most", "triangle-down"), ("measured", "circle")}
    static = ies.plot_vs_obs(backend="matplotlib")
    drawn = {line.get_marker() for ax in static.axes for line in ax.lines if line.get_linestyle() == "None"}
    assert drawn == {"v", "o"}
