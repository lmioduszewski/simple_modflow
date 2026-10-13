"""The recalibration rehearsal: a Ten Trails-shaped truth calibrated back (ledger 589).

``myflopy.modflow.mf6.canonical_rehearsal`` builds a start model and a truth that
differs from it by known multipliers, declared through the package-first API;
``examples/mf6/notebooks/canonical_08_recalibration_rehearsal.ipynb`` calibrates the
one back to the other with the API the Ten Trails recalibration will use. Pinned here:

- the case has Ten Trails' shape (two RIV packages ``riv``/``crisp``, an array RCH,
  layered STO, NEWTON, one IMPLICIT LAK of two lakes), and the truth differs from the
  start by EXACTLY its multipliers, package by package (no mf6 needed);
- the synthetic readings bin back to the truth exactly -- average-year medians and the
  history mode's period ends -- while the history mode's distractors would move a median;
- the TRUTH ROUND TRIP (the notebook up to it): the template at the truth's multipliers
  reproduces the noise-free data of every family; the start does not. It is the test that
  catches a misrouted mapping: a ``package="crisp"`` multiplier sent to ``riv``, the wet
  and dry months swapped, average-year targets binned onto the wrong cycle, or the two
  lakes' bed leakances swapped each fail it;
- the whole notebook, run once: the user's pass criteria (Q9) on its
  ``rehearsal_report`` -- the round trips, the equilibrium, phi and crisp's log gap
  asserted, the two posterior-coverage criteria held as a non-strict xfail until the user
  decides (ledger 589: not met at any ensemble size measured; ledger
  590: the ensemble has not converged).
"""

from __future__ import annotations

import shutil
import warnings

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal
from tests.test_examples_run import execute_notebook

NOTEBOOK = "canonical_08_recalibration_rehearsal"
NEEDS_MF6 = pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
NEEDS_PEST = pytest.mark.skipif(
    shutil.which("mf6") is None or shutil.which("pestpp-ies") is None,
    reason="needs mf6 and pestpp-ies on PATH",
)


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return rehearsal.build_rehearsal(tmp_path_factory.mktemp("rehearsal"))


def _built(case, simulation):
    """``simulation``'s flow model, built in memory (no MF6 run)."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = case.project.prepare_run(f"shape_{simulation}", simulation, overwrite=True)
    return run.model(case.model_name)


# -- the case ---------------------------------------------------------------------------


def test_the_start_has_ten_trails_shape(case):
    fl = _built(case, case.start_simulation).fl
    rivers = fl.riv
    assert sorted(river.package_name for river in rivers) == ["crisp", "riv"]
    assert type(fl.rch).__name__ == "ModflowGwfrcha", "an ARRAY RCH (READASARRAYS)"
    assert type(fl.lak).__name__ == "ModflowGwflakImplicit"
    assert len(fl.lak.packagedata.get_data()) == 2, "one LAK package, two lakes"
    sy = fl.sto.sy.array
    assert sy.shape == (5, 441) and len(np.unique(sy.max(axis=1))) == 5, "Sy per layer"
    assert fl.riv[0].model_or_sim.name_file.newtonoptions.get_data()[0][0] == "NEWTON"
    assert case.year.nper == 1 + 12 * 3
    assert case.wells["screen_top"].notna().sum() == 2, "two wells screened across layers 1-2"


def test_the_truth_differs_from_the_start_by_exactly_its_multipliers(case):
    model = _built(case, case.start_simulation)
    start, truth = model.fl, _built(case, case.truth_simulation).fl
    # The cells the pilot-point zones label (by centroid) are the cells the truth's
    # mf.Polygons multipliers took (by majority, with default=1.0): a zone is the union
    # of its cells, so the two rules cannot differ.
    labels = np.asarray(case.zone_spec(model).values).reshape(-1)[:441]
    ratio = truth.npf.k.array / start.npf.k.array
    expected = np.ones_like(ratio)
    for label, factor in rehearsal.TRUTH_MULTIPLIERS["k"].items():
        expected[rehearsal.K_LAYER, labels == label] = factor
    np.testing.assert_allclose(ratio, expected, rtol=1e-12)
    assert (labels > 0).sum() > 40, "both zones hold cells"
    np.testing.assert_allclose(truth.sto.sy.array / start.sto.sy.array, 0.6, rtol=1e-12)
    for package, factor in rehearsal.TRUTH_MULTIPLIERS["riv.cond"].items():
        pick = {river.package_name: river for river in start.riv}[package]
        moved = {river.package_name: river for river in truth.riv}[package]
        np.testing.assert_allclose(moved.stress_period_data.get_data(0)["cond"],
                                   factor * pick.stress_period_data.get_data(0)["cond"])
    months = case.year.periods
    for row in months.itertuples():
        mine = truth.rch.recharge.get_data(row.per) / start.rch.recharge.get_data(row.per)
        if row.steady:
            factor = 1.0  # the steady period keeps the start's mean: no month group reaches it
        else:
            factor = rehearsal.TRUTH_MULTIPLIERS["recharge"][
                "wet" if int(row.month) in rehearsal.WET_MONTHS else "dry"]
        np.testing.assert_allclose(mine, factor, rtol=1e-12)
    lakes = {row["boundname"]: row["ifno"] for row in start.lak.packagedata.get_data()}
    before = start.lak.connectiondata.get_data()
    after = truth.lak.connectiondata.get_data()
    for name, number in lakes.items():
        rows = before["ifno"] == number
        np.testing.assert_allclose(after["bedleak"][rows].astype(float)
                                   / before["bedleak"][rows].astype(float),
                                   rehearsal.TRUTH_MULTIPLIERS["lak.bedleak"][name])


def test_two_rivers_are_never_picked_in_silence(case):
    cal = _built(case, case.start_simulation).pest("refusal")
    with pytest.raises(mf.SeveralPackagesError, match="package='crisp'"):
        cal.parameterize("riv.cond")


@pytest.mark.parametrize(("kwargs", "match"), [
    ({"cycles": 1}, "at least 2 cycles"),
    ({"nlay": 2}, "3 to 5 layers"),
    ({"nlay": 6}, "3 to 5 layers"),
])
def test_build_rehearsal_refuses(tmp_path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        rehearsal.build_rehearsal(tmp_path, **kwargs)


@pytest.mark.parametrize(("kwargs", "match"), [
    ({"mode": "average year"}, "average_year"),
    ({"readings_per_month": 2}, "ODD count"),
    ({"mode": "history", "years": 3}, "neither can act"),
    ({"noise_sd": -0.1}, ">= 0"),
])
def test_synthetic_hydrographs_refuses_before_running_anything(case, kwargs, match, monkeypatch):
    # The case is module-scoped, so whether its truth already ran depends on test order:
    # make running it fail instead of reading case._runs.
    monkeypatch.setattr(rehearsal, "_truth_heads", lambda *_: pytest.fail("a refusal ran a model"))
    with pytest.raises(ValueError, match=match):
        rehearsal.synthetic_hydrographs(case, **kwargs)


# -- the data ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def truth_heads(case):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return rehearsal._truth_heads(case)


@pytest.mark.slow
@NEEDS_MF6
def test_average_year_readings_bin_to_exactly_the_truth(case, truth_heads):
    path = rehearsal.synthetic_hydrographs(case, seed=7)
    # round_trip: pandas' default parser reads some 17-digit numbers one ulp off.
    readings = pd.read_csv(path, parse_dates=["date"], float_precision="round_trip").set_index("date")
    months = case.year.months_of(-1)
    for name in readings.columns:
        column = readings[name].dropna()
        by_year_month = column.groupby([column.index.year, column.index.month])
        assert (by_year_month.size() == 3).all(), "three readings a month, every year read"
        for (_, month), median in by_year_month.median().items():
            assert median == truth_heads.loc[months[month], name], (name, month)
    targets = mf.HeadTargets.from_hydrographs(case.wells, path, model=case.truth_model(),
                                              mode="average_year")
    table = targets.to_wide().set_index("time")
    for month, period in months.items():
        for name in truth_heads.columns:
            value = table.loc[period, name]
            if month in rehearsal.GAPS.get(name, ()):
                assert np.isnan(value), (name, month)
            else:  # from_hydrographs reads the CSV with pandas' default parser: an ulp
                assert value == pytest.approx(truth_heads.loc[period, name], rel=1e-15, abs=0)


@pytest.mark.slow
@NEEDS_MF6
def test_history_readings_bin_to_the_period_ends_and_the_distractors_would_not(case, truth_heads):
    path = rehearsal.synthetic_hydrographs(case, mode="history", seed=7)
    readings = pd.read_csv(path, parse_dates=["date"]).set_index("date")
    ends = set(pd.to_datetime(case.year.periods["end"]))
    off = readings.index[~readings.index.isin(list(ends))]
    assert len(off) == case.year.nper - 1, "one distractor in each transient period"
    model = case.truth_model()
    nearest = mf.HeadTargets.from_hydrographs(case.wells, path, model=model).to_wide().set_index("time")
    median = mf.HeadTargets.from_hydrographs(case.wells, path, model=model, rule="median")
    median = median.to_wide().set_index("time")
    gaps = {name: set(months) for name, months in rehearsal.GAPS.items()}
    for row in case.year.periods.itertuples():
        for name in truth_heads.columns:
            truth = truth_heads.loc[row.per, name]
            if not row.steady and int(row.month) in gaps.get(name, ()):
                assert np.isnan(nearest.loc[row.per, name]), (name, row.per)
                continue
            assert nearest.loc[row.per, name] == pytest.approx(truth, rel=1e-15, abs=0)
            if not row.steady:  # the median of the period's two readings is moved
                assert median.loc[row.per, name] == pytest.approx(
                    truth + rehearsal.DISTRACTOR_OFFSET / 2, rel=1e-12)


# -- the round trips --------------------------------------------------------------------


@pytest.mark.slow
@NEEDS_PEST
def test_truth_round_trip_reproduces_the_data(tmp_path, monkeypatch):
    """The notebook up to its round trip: the truth's multipliers reproduce every family."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ns = execute_notebook(NOTEBOOK, tmp_path, monkeypatch, until="# THE TRUTH ROUND TRIP")
    truth, start = ns["truth_trip"], ns["start_trip"]
    assert set(truth["group"]) == {"hds", "hdiff", "lak_stage", "riv_flow", "l_rim"}
    worst = truth.groupby("group")["ratio"].max()
    assert (worst <= 1.0).all(), worst
    # ... and not because nothing moves: the start misses every family but the bound.
    missed = start.groupby("group")["ratio"].max()
    assert (missed.drop("l_rim") > 100).all(), missed
    assert ns["equilibrium"].converged


@pytest.fixture(scope="module")
def rehearsed(tmp_path_factory):
    """The whole notebook, run once at its own knobs; its ``rehearsal_report``."""

    with pytest.MonkeyPatch.context() as patch, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        namespace = execute_notebook(NOTEBOOK, tmp_path_factory.mktemp("canonical_08"), patch)
    return namespace["rehearsal_report"]


#: The user's criteria (Q9) a run of the notebook's size meets on every seed measured
#: (``docs/reproducers/rehearsal_recovery_budget.py``).
MET = ("truth_round_trip", "write_back", "equilibrium", "phi_ratio", "crisp_gap_closed")
#: The two it does not: the scalar truths inside their posterior 5-95% bands, and the
#: forecast's truth inside the posterior's range. Measured at 16x2 to 48x3 over three
#: seeds, one to six of the six scalar truths fall outside and the forecast range misses
#: on some seeds -- an ensemble that has not converged at this budget: the realization
#: IES ends on fits the noisy data 1.6-3x worse than the truth's multipliers do (ledger
#: 589, ledger 590), pending the user's decision.
UNDECIDED = ("scalars_inside", "forecast_covered")


@pytest.mark.slow
@NEEDS_PEST
def test_rehearsal_notebook_recovers_the_truth(rehearsed):
    assert set(rehearsed["passes"]) == set(MET) | set(UNDECIDED)
    failed = sorted(name for name in MET if not rehearsed["passes"][name])
    assert failed == [], {name: rehearsed.get(name) for name in failed}
    assert rehearsed["start_misfit"] > 100, "the start is far from the data"


@pytest.mark.slow
@NEEDS_PEST
@pytest.mark.xfail(strict=False, reason=(
    "an unconverged ensemble at 16x2: the posterior 5-95% bands miss 5-6 of 6 scalar truths "
    "and the forecast range misses on some seeds -- ledger 590, the user's "
    "decision"))
def test_rehearsal_posterior_covers_the_truth(rehearsed):
    failed = sorted(name for name in UNDECIDED if not rehearsed["passes"][name])
    assert failed == [], {"scalars": rehearsed["scalars"],
                          "forecast": (rehearsed["forecast_truth"], rehearsed["forecast_posterior"])}
