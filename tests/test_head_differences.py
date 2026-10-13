"""Difference observations give PEST each hydrograph's SHAPE (ledger 565).

``HeadTargets.differences(kind=...)`` returns a ``HeadDifferences`` set (``ConcTargets``
a ``ConcDifferences``): ``kind="mean"`` -- each observation label minus its own mean
over ITS measured periods, at the source weight -- or ``kind="previous"`` -- each
measured period minus the previous one of that label, at the two readings' combined
weight (``w / sqrt(2)`` for two equal weights). It answers ``compare()``/``stats()``
off-PEST and is registered with ``cal.observe(...)`` under the prefix ``hdiff``
(``cdiff``), one PEST group per set. The forward run builds the differences from the
ONE period-end head reader (ledger 557), so a screened well is one label and its
difference is of the T-weighted head (ledger 556); rows with no defined difference are
written 0.0 at weight 0 (PEST++ fails on an empty field, ledger 351).

The model-backed tests use ``test_pest_handed_targets.py``'s toy (three 10-day periods)
and ``test_screened_head_targets.py``'s screened model; PESTPP-IES runs the prior.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.specs import (
    ConcDifferenceObservationSpec,
    HeadDifferenceObservationSpec,
)
from myflopy.modflow.mf6.target_values import target_differences
from tests.test_pest_handed_targets import run_toy
from tests.test_screened_head_targets import WELL_CELL, K, _by_hand, _forward_run, _write_run
from tests.test_screened_head_targets import _heads as _screen_heads

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return run_toy(tmp_path_factory.mktemp("head_differences"))


def _simulated(model) -> pd.DataFrame:
    """The toy's own period-end heads at W1 (cell 1) and W2 (cell 2)."""

    base = mf.HeadTargets(locations={"W1": 1, "W2": 2}, values={"time": [0], "W1": [1.0], "W2": [1.0]})
    return base.simulated_heads(model).set_index("per")


def _targets(values, *, weights=None, groups=None, **kwargs):
    locations = pd.DataFrame({"name": ["W1", "W2"], "cell": [1, 2], "layer": 0})
    if weights is not None:
        locations["weight"] = weights
    if groups is not None:
        locations["group"] = groups
    return mf.HeadTargets(locations=locations, values=pd.DataFrame(values), **kwargs)


# -- the rule, by hand -------------------------------------------------------------------


def test_the_one_rule_by_hand():
    table = pd.DataFrame({"per": [0, 1, 2, 3], "a": [1.0, 2.0, 4.0, 8.0], "b": [5.0, 5.0, 9.0, 1.0]})
    periods = pd.DataFrame({"obsname": ["a", "a", "a", "b", "b"], "per": [0, 2, 3, 1, 2]})
    mean = target_differences(table, periods, kind="mean")
    assert mean["a"].tolist()[0] == pytest.approx(1.0 - 13 / 3)
    assert np.isnan(mean.loc[1, "a"])
    assert mean["b"].tolist()[1:3] == [-2.0, 2.0]
    previous = target_differences(table, periods, kind="previous")
    assert np.isnan(previous.loc[0, "a"]) and np.isnan(previous.loc[1, "a"])
    assert previous["a"].tolist()[2:] == [3.0, 4.0]
    assert np.isnan(previous.loc[1, "b"]) and previous.loc[2, "b"] == 4.0
    assert list(previous.columns) == ["per", "a", "b"]



def test_a_measured_row_with_no_number_stays_empty_for_pest(tmp_path, monkeypatch):
    """Only a row with NO difference is written 0.0 (ledger 565): a model
    that wrote NaN at a measured period left an empty field, so PEST++ fails that run --
    filling it 0.0 scored a broken realization as a flat hydrograph."""

    from myflopy.modflow.mf6 import target_values
    from myflopy.modflow.mf6.pest.forward_run import _write_target_difference_csv

    table = pd.DataFrame({"per": [0, 1, 2], "W1": [10.0, np.nan, 14.0], "W2": [5.0, 6.0, 9.0]})
    monkeypatch.setattr(target_values, "read_target_values", lambda *args, **kwargs: table.copy())
    pd.DataFrame({"obsname": ["W1", "W1", "W1", "W2", "W2"], "per": [0, 1, 2, 0, 2]}).to_csv(
        tmp_path / "periods.csv", index=False)
    pd.DataFrame({"obsname": ["W1", "W2"], "name": ["W1", "W2"], "layer": [0, 0], "cell": [1, 2]}).to_csv(
        tmp_path / "map.csv", index=False)
    for kind in ("mean", "previous"):
        _write_target_difference_csv("m", tmp_path / "map.csv", tmp_path / "out.csv",
                                     tmp_path / "periods.csv", kind, sim_ws=tmp_path)
        written = pd.read_csv(tmp_path / "out.csv").set_index("per")
        # W1's NaN reaches every difference it is part of, and stays empty.
        assert written["W1"].isna().sum() == (3 if kind == "mean" else 2)
        # W2 is unmeasured in period 1 (and has no previous in period 0): 0.0 there.
        assert written.loc[1, "W2"] == 0.0
        assert written.loc[2, "W2"] == pytest.approx(2.0 if kind == "mean" else 4.0)
    assert target_differences(table, pd.read_csv(tmp_path / "periods.csv"), kind="previous",
                              undefined=0.0)["W1"].tolist()[0] == 0.0

def test_mean_differences_ignore_a_constant_offset_where_absolute_heads_do_not(toy):
    """The feature's point: a datum error moves every reading of a well alike."""

    sim = _simulated(toy)
    values = {"time": [0, 1, 2], "W1": (sim["W1"] + 7.0).tolist(), "W2": (sim["W2"] - 3.0).tolist()}
    targets = _targets(values)
    absolute = targets.compare(toy)
    assert absolute.set_index("name")["residual"].round(9).to_dict() == {"W1": -7.0, "W2": 3.0}
    diffs = targets.differences()
    assert isinstance(diffs, mf.HeadDifferences) and diffs.kind == "mean"
    frame = diffs.compare(toy)
    assert len(frame) == 6
    assert np.abs(frame["residual"]).max() < 1e-9
    expected = sim["W1"] + 7.0 - (sim["W1"] + 7.0).mean()
    assert frame.loc[frame["name"] == "W1", "target_difference"].to_numpy() == pytest.approx(expected.to_numpy())
    assert diffs.stats(toy).iloc[0]["n"] == 6


def test_previous_differences_and_their_weights_by_hand(toy):
    values = {"time": [0, 1, 2], "W1": [10.0, 12.0, 11.0], "W2": [5.0, None, 9.0]}
    weights = [2.0, 1.0]
    frame = _targets(values, weights=weights).differences("previous").compare(toy).set_index(["name", "per"])
    sim = _simulated(toy)
    assert np.isnan(frame.loc[("W1", 0), "target_difference"])
    assert frame.loc[("W1", 1), "target_difference"] == 2.0
    assert frame.loc[("W1", 2), "target_difference"] == -1.0
    assert frame.loc[("W2", 2), "target_difference"] == 4.0  # 9 - 5: the previous MEASURED
    assert frame.loc[("W2", 2), "sim_difference"] == pytest.approx(sim.loc[2, "W2"] - sim.loc[0, "W2"])
    assert frame.loc[("W1", 1), "weight"] == pytest.approx(2.0 / np.sqrt(2.0))
    assert frame.loc[("W2", 2), "weight"] == pytest.approx(1.0 / np.sqrt(2.0))
    # Two readings of unequal weight: 1 / sqrt(1/w1^2 + 1/w2^2).
    reading = pd.DataFrame({"name": ["W1"] * 3, "time": [0, 1, 2], "head": [10.0, 12.0, 11.0],
                            "weight": [1.0, 2.0, 2.0]})
    mixed = mf.HeadTargets(locations={"W1": 1}, values=reading).differences("previous").get()
    assert mixed.set_index("per").loc[1, "weight"] == pytest.approx(1.0 / np.sqrt(1.0 + 0.25))
    mean = _targets(values, weights=weights).differences().get().set_index(["name", "per"])
    assert mean.loc[("W1", 1), "weight"] == 2.0
    override = _targets(values, weights=weights).differences("previous", weight=0.5).get()
    assert set(override["weight"].dropna()) == {0.5}


def test_a_label_with_one_reading_is_left_out_and_named_and_none_left_is_refused():
    values = {"time": [0, 1], "W1": [10.0, 12.0], "W2": [5.0, None]}
    with pytest.warns(UserWarning, match=r"fewer than 2 measured periods.*left out of the differences: W2"):
        diffs = _targets(values).differences()
    assert diffs.names == ["W1"] and diffs.left_out == ["W2"]
    assert set(diffs.get()["name"]) == {"W1"}
    with pytest.raises(ValueError, match="no label has 2 measured periods"):
        _targets({"time": [0], "W1": [1.0], "W2": [2.0]}).differences()


def test_the_refusals():
    values = {"time": [0, 1], "W1": [10.0, 12.0], "W2": [5.0, 6.0]}
    with pytest.raises(ValueError, match=r"kind='diff'.*'mean'.*'previous'"):
        _targets(values).differences("diff")
    with pytest.raises(ValueError, match=r"a BOUND.*has no mean"):
        _targets(values, inequality="at_most").differences()
    with pytest.raises(ValueError, match=r"weight=-1"):
        _targets(values).differences(weight=-1)


def test_to_wide_and_summary():
    values = {"time": [0, 1, 2], "W1": [10.0, 12.0, 14.0], "W2": [5.0, None, 9.0]}
    diffs = _targets(values).differences()
    wide = diffs.to_wide()
    assert list(wide.columns) == ["time", "W1", "W2"]
    assert wide["W1"].tolist() == [-2.0, 0.0, 2.0]
    assert np.isnan(wide.loc[1, "W2"])
    summary = diffs.summary().iloc[0]
    assert summary["kind"] == "mean" and summary["n_locations"] == 2 and summary["n_rows"] == 5


# -- PEST ---------------------------------------------------------------------------------


def _build(model, workspace, *observed, noptmax=0):
    cal = model.pest("diffs", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in observed:
        cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("diffs.pst", noptmax=noptmax)
    return cal, pst


def test_the_build_hands_pest_the_differences_in_one_group(toy, tmp_path):
    values = {"time": [0, 1, 2], "W1": [10.0, 12.0, 14.0], "W2": [5.0, None, 9.0]}
    diffs = _targets(values, weights=[2.0, None], groups=["north", "south"]).differences()
    cal, pst = _build(toy, tmp_path / "pest", diffs)
    obs = pst.observation_data
    measured = obs.loc[obs["weight"] > 0]
    usecol = "oname:hdiff_otype:lst_usecol:{}_per:{}"
    assert measured["obsval"].to_dict() == pytest.approx({
        usecol.format("w1", 0): -2.0, usecol.format("w1", 1): 0.0, usecol.format("w1", 2): 2.0,
        usecol.format("w2", 0): -2.0, usecol.format("w2", 2): 2.0,
    })
    assert set(measured["obgnme"]) == {"hdiff"}  # one group per set, labels not read
    assert measured.loc[usecol.format("w1", 1), "weight"] == 2.0
    assert measured.loc[usecol.format("w2", 0), "weight"] == 1.0
    unmeasured = obs.loc[usecol.format("w2", 1)]
    assert unmeasured["weight"] == 0 and unmeasured["obgnme"] == "hdiff_unmeasured"
    # The base table is the forward run's: the simulated differences, 0.0 where none.
    base = pd.read_csv(cal.template_workspace / "hdiff_simulated_differences.csv").set_index("per")
    sim = _simulated(toy)
    assert base["W1"].to_numpy() == pytest.approx((sim["W1"] - sim["W1"].mean()).to_numpy())
    assert base.loc[1, "W2"] == 0.0
    assert base.loc[0, "W2"] == pytest.approx(sim.loc[0, "W2"] - sim.loc[[0, 2], "W2"].mean())
    periods = pd.read_csv(cal.template_workspace / "hdiff_difference_periods.csv")
    assert periods.to_dict("list") == {"obsname": ["W1", "W1", "W1", "W2", "W2"], "per": [0, 1, 2, 0, 2]}
    entry = next(item["metadata"] for item in cal._prepared_observations if item["prefix"] == "hdiff")
    assert entry["kind"] == "head_differences" and entry["difference_kind"] == "mean"
    assert entry["geometry"] == "points"



def test_a_reading_of_weight_0_moves_no_difference(toy, tmp_path):
    """A silenced outlier is no measurement of a difference (ledger 566): it
    moved the MEAN, so the target of every weighted reading of its well -- 10 and 12 beside
    an outlier of 50 at weight 0 read -14 and -12, where they are -1 and +1."""

    readings = pd.DataFrame({"name": ["W1"] * 3 + ["W2"] * 2, "time": [0, 1, 2, 0, 2],
                             "head": [10.0, 12.0, 50.0, 5.0, 9.0],
                             "weight": [1.0, 1.0, 0.0, 1.0, 1.0]})
    targets = mf.HeadTargets(locations={"W1": 1, "W2": 2}, values=readings)
    mean = targets.differences().get().set_index(["name", "per"])
    assert mean.loc[("W1", 0), "target_difference"] == -1.0
    assert mean.loc[("W1", 1), "target_difference"] == 1.0
    assert np.isnan(mean.loc[("W1", 2), "target_difference"])
    # "previous" differences ACROSS it, at the two counted readings' weight.
    middle = readings.assign(head=[10.0, 50.0, 11.0, 5.0, 9.0], weight=[1.0, 0.0, 1.0, 1.0, 1.0])
    previous = mf.HeadTargets(locations={"W1": 1, "W2": 2}, values=middle).differences("previous")
    step = previous.get().set_index(["name", "per"])
    assert np.isnan(step.loc[("W1", 1), "target_difference"])
    assert step.loc[("W1", 2), "target_difference"] == 1.0
    assert step.loc[("W1", 2), "weight"] == pytest.approx(1.0 / np.sqrt(2.0))
    # A well whose second reading is silenced has one reading left: it is left out.
    with pytest.warns(UserWarning, match=r"weight 0 does not count.*left out of the differences: W1"):
        lone = mf.HeadTargets(locations={"W1": 1, "W2": 2},
                              values=readings.assign(weight=[1.0, 0.0, 0.0, 1.0, 1.0])).differences()
    assert lone.names == ["W2"]
    # The build hands PEST the same targets, and its forward run the same periods.
    cal, pst = _build(toy, tmp_path / "pest", targets.differences())
    obs = pst.observation_data
    usecol = "oname:hdiff_otype:lst_usecol:w1_per:{}"
    assert obs.loc[usecol.format(0), "obsval"] == pytest.approx(-1.0)
    assert obs.loc[usecol.format(1), "obsval"] == pytest.approx(1.0)
    assert obs.loc[usecol.format(2), ["weight", "obgnme"]].tolist() == [0.0, "hdiff_unmeasured"]
    periods = pd.read_csv(cal.template_workspace / "hdiff_difference_periods.csv")
    assert periods.to_dict("list") == {"obsname": ["W1", "W1", "W2", "W2"], "per": [0, 1, 0, 2]}
    base = pd.read_csv(cal.template_workspace / "hdiff_simulated_differences.csv").set_index("per")
    sim = _simulated(toy)
    assert base.loc[0, "W1"] == pytest.approx(sim.loc[0, "W1"] - sim.loc[[0, 1], "W1"].mean())
    assert base.loc[2, "W1"] == 0.0


@pytest.fixture(scope="module")
def offset_run(toy, tmp_path_factory):
    """Targets = the model's own heads + a per-well constant, absolute and differenced."""

    sim = _simulated(toy)
    values = {"time": [0, 1, 2], "W1": (sim["W1"] + 7.0).tolist(), "W2": (sim["W2"] - 3.0).tolist()}
    targets = _targets(values)
    cal, _pst = _build(toy, tmp_path_factory.mktemp("offset") / "pest", targets, targets.differences())
    return cal, cal.prior(reals=3)


def test_pestpp_ies_scores_an_offset_model_zero_on_differences_and_not_on_heads(offset_run):
    _cal, ies = offset_run
    table = pd.read_csv(ies.workspace / "diffs.phi.group.csv")
    base = table.loc[table["par_realization"].astype(str) == "base"].iloc[0]
    assert base["hdiff"] < 1e-6
    assert base["hds"] == pytest.approx(3 * 49.0 + 3 * 9.0, rel=1e-4)
    # The other realizations move K, and the hydrographs' shapes with it.
    assert (table.loc[table["par_realization"].astype(str) != "base", "hdiff"] > 0).any()
    titles = {annotation.text for annotation in ies.plot_vs_obs().layout.annotations}
    assert {"hdiff: W1 [hdiff]", "hdiff: W2 [hdiff]"} <= titles


def test_a_difference_set_of_a_screened_dated_head_set_builds_and_runs(tmp_path):
    """Where three batches meet: dated readings (ledger 543), a screened well read by
    the one period-end reader (ledgers 556-557), differenced in the forward run with
    the REALIZATION's K."""

    workspace = _write_run(tmp_path / "dated", nper=3)
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0)
    sim.tdis.start_date_time.set_data("2020-01-01")
    sim.tdis.write()
    model = mf.load_mf6_run(workspace)
    locations = pd.DataFrame({"name": ["W1"], "cell": [WELL_CELL], "screen_top": [70.0],
                              "screen_bottom": [20.0]})
    readings = pd.DataFrame({"date": ["2020-01-02", "2020-01-03", "2020-01-04"], "W1": [85.0, 86.0, 84.0]})
    targets = mf.HeadTargets.from_hydrographs(locations, readings, model=model)
    diffs = targets.differences("previous")
    cal = model.pest("hy", workspace=tmp_path / "hy", start_datetime="2020-01-01")
    cal.observe(diffs)
    cal.parameterize("k", layers=[2], bounds=(0.1, 20.0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("hy.pst", noptmax=0)
    measured = pst.observation_data.loc[pst.observation_data["weight"] > 0]
    assert measured["obsval"].tolist() == [1.0, -2.0]
    assert measured["weight"].tolist() == pytest.approx([1 / np.sqrt(2)] * 2)
    entry = next(item["metadata"] for item in cal._prepared_observations)
    assert entry["screened"] and len(entry["period_end_dates"]) == 3
    heads = _screen_heads(workspace)
    t_head = [_by_hand(heads[per], K, 70.0, 20.0) for per in range(3)]
    base = pd.read_csv(cal.template_workspace / "hdiff_simulated_differences.csv")
    assert base["W1"].tolist() == pytest.approx([0.0, t_head[1] - t_head[0], t_head[2] - t_head[1]])
    group = next(g for g in pst.parameter_data["pargp"].unique() if g.startswith("k"))
    worker = _forward_run(cal, pst, tmp_path / "worker", {group: 10.0})
    realized = (K[0], K[1], 10.0 * K[2])
    heads = _screen_heads(worker)
    t_head = [_by_hand(heads[per], realized, 70.0, 20.0) for per in range(3)]
    written = pd.read_csv(worker / "hdiff_simulated_differences.csv")
    assert written["W1"].tolist() == pytest.approx(
        [0.0, t_head[1] - t_head[0], t_head[2] - t_head[1]], rel=1e-10)


def test_a_concentration_set_differences_under_cdiff():
    conc = mf.ConcTargets(locations={"a": 1}, values={"time": [0, 1], "a": [1.0, 3.0]})
    diffs = conc.differences()
    assert isinstance(diffs, mf.ConcDifferences)
    assert diffs.get()["target_difference"].tolist() == [-1.0, 1.0]
    from myflopy.modflow.mf6.pest.project import PestProject

    spec = PestProject._coerce_observation_spec(None, diffs, prefix=None)
    assert isinstance(spec, ConcDifferenceObservationSpec) and spec.prefix == "cdiff"
    head = PestProject._coerce_observation_spec(None, _targets(
        {"time": [0, 1], "W1": [1.0, 2.0], "W2": [1.0, 2.0]}).differences(), prefix=None)
    assert isinstance(head, HeadDifferenceObservationSpec) and head.prefix == "hdiff"


def test_a_concentration_difference_set_builds_from_the_transport_models_file(tmp_path):
    """``cdiff`` reads the GWT model's concentration file through the one reader and
    writes, per label, the simulated concentrations minus their measured-period mean."""

    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
    from myflopy.modflow.mf6.canonical_transport import build_canonical_transport_calibration_demo

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r"mf\.lak: lake\(s\) NOT excavated")
        demo = build_canonical_transport_calibration_demo(
            tmp_path / "model", config=CanonicalModelConfig.testing()
        )
    diffs = demo.conc_targets.differences()
    cal = demo.model.pest("conc_diffs", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.05, 2.0), physical=(0.01, 300.0))
    cal.observe(diffs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obs = cal.build("conc_diffs.pst", noptmax=0).observation_data
    measured = obs.loc[obs["weight"] > 0]
    assert set(measured["obgnme"]) == {"cdiff"}
    expected = diffs.get().dropna(subset=["target_difference"])
    assert sorted(measured["obsval"].round(9)) == sorted(expected["target_difference"].round(9))
    base = pd.read_csv(cal.template_workspace / "cdiff_simulated_differences.csv").set_index("per")
    simulated = demo.conc_targets.simulated_conc(demo.transport.transport_view()).set_index("per")
    periods = pd.read_csv(cal.template_workspace / "cdiff_difference_periods.csv")
    label = str(periods["obsname"].iloc[0])
    pers = periods.loc[periods["obsname"] == label, "per"].tolist()
    assert base.loc[pers, label].to_numpy() == pytest.approx(
        (simulated.loc[pers, label] - simulated.loc[pers, label].mean()).to_numpy())


def test_a_screened_difference_set_has_its_k_restored_per_layer_with_no_k_parameter(tmp_path):
    """The build re-stores a whole-grid NPF K one file per layer whenever a set reads
    screens (ledger 556) -- a DIFFERENCE set of a screened source included, which
    reads the same T-weighted heads; with only a CHD parameter nothing else would."""

    workspace = _write_run(tmp_path / "flat", nper=2)
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0)
    sim.get_model("m").npf.k.set_data(5.0, layer=None)
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    model = mf.load_mf6_run(workspace)
    locations = pd.DataFrame({"name": ["W1"], "cell": [WELL_CELL], "screen_top": [70.0],
                              "screen_bottom": [20.0]})
    diffs = mf.HeadTargets(locations=locations, values={"time": [0, 1], "W1": [85.0, 86.0]}).differences()
    cal = model.pest("flat", workspace=tmp_path / "flat_pest", start_datetime="2020-01-01")
    cal.observe(diffs)
    cal.parameterize("chd", style="constant")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("flat.pst", noptmax=0)
    heads = _screen_heads(workspace)
    t_head = [_by_hand(heads[per], (5.0, 5.0, 5.0), 70.0, 20.0) for per in range(2)]
    base = pd.read_csv(cal.template_workspace / "hdiff_simulated_differences.csv")
    mean = np.mean(t_head)
    assert base["W1"].tolist() == pytest.approx([t_head[0] - mean, t_head[1] - mean])
