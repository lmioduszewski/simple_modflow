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
- the whole notebook, run once at its CI budget (16 realizations x 2 iterations, one
  lambda): the user's per-run pass criteria (Q9) on its ``rehearsal_report`` -- the
  round trips, the equilibrium, phi and crisp's log gap. Coverage is not asked of it:
  at that budget the ensemble has not converged (ledgers 589, 590);
- the recovery's ONE scoring rule (ledger 593): a scalar's band and its ``z``
  come from the normal fitted to the posterior REALIZATIONS' log10 values (never the
  ``base`` row), so ``inside`` is ``|z| <= Z_BAND``; and the pooled verdict over several
  runs, whose thresholds an exact posterior fails at most 5% of the time -- checked by
  Monte Carlo of exact posteriors scored as the rehearsal scores them, each ``z`` against
  a FINITE ensemble's own mean and sd (ledger 594: a standard normal's
  thresholds failed 4.6%, not the 3.1% they claimed);
- OPT-IN, the pooled coverage of CONVERGED runs: with
  ``SIMPLE_MODFLOW_REHEARSAL_COVERAGE=1`` the notebook runs at 64 realizations x 3
  iterations with PESTPP-IES's own lambda testing on six seeds (half an hour at 12
  workers) and ``rehearsal.pooled_recovery`` must pass --
  the user's criterion in place of "every scalar truth inside its band in every run",
  which an exact posterior meets about half the time (0.9^6). Measured with the
  notebook's declaration (levels and their differences both at the reading's weight):
  32 of 36 inside, largest |z| 3.73 against a bound of 3.88, forecast 6 of 6.
"""

from __future__ import annotations

import json
import os
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal
from tests.test_examples_run import execute_notebook

NOTEBOOK = "canonical_08_recalibration_rehearsal"
ROOT = Path(__file__).resolve().parents[1]
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
    ({"head_noise": 0.0}, "head_noise=0.0"),
    ({"head_noise": float("nan")}, "head_noise=nan"),
    ({"head_noise": True}, "head_noise=True"),
])
def test_build_rehearsal_refuses(tmp_path, kwargs, match):
    with pytest.raises(ValueError, match=match):
        rehearsal.build_rehearsal(tmp_path, **kwargs)


def test_the_head_weight_follows_the_head_noise(case, tmp_path):
    """Each well's weight is 1 / the noise drawn on its readings (ledger 606).

    It was a literal 1 / 0.05: editing the notebook's ``NOISE["head"]`` moved the data's
    noise and left every head weight at 20.
    """

    assert (case.wells["weight"] == 1.0 / 0.05).all(), "the default is 0.05 ft"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        noisier = rehearsal.build_rehearsal(tmp_path, head_noise=0.2)
    assert (noisier.wells["weight"] == 5.0).all()
    # The notebook hands its one knob to both: the data's noise and the weights.
    source = "\n".join(
        "".join(cell["source"]) for cell in json.loads(
            (ROOT / "examples" / "mf6" / "notebooks" / f"{NOTEBOOK}.ipynb").read_text(
                encoding="utf-8"))["cells"])
    assert 'build_rehearsal(WORKSPACE, cycles=CYCLES, head_noise=NOISE["head"])' in source
    assert 'synthetic_hydrographs(case, noise_sd=NOISE["head"]' in source


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


#: The user's per-run criteria (Q9), which a run of the notebook's CI size meets on every
#: seed measured (``docs/reproducers/rehearsal_recovery_budget.py``). The scalar truths
#: inside their bands and the forecast covered are judged POOLED over several converged
#: runs instead (ledger 593): the opt-in test below.
MET = ("truth_round_trip", "write_back", "equilibrium", "phi_ratio", "crisp_gap_closed")


@pytest.mark.slow
@NEEDS_PEST
def test_rehearsal_notebook_recovers_the_truth(rehearsed):
    assert set(rehearsed["passes"]) == set(MET)
    failed = sorted(name for name in MET if not rehearsed["passes"][name])
    assert failed == [], {name: rehearsed.get(name) for name in failed}
    assert rehearsed["start_misfit"] > 100, "the start is far from the data"
    # What the pooled verdict reads off each run is in the report, by the one rule.
    for group, row in rehearsed["scalars"].items():
        assert rehearsed["scalars_inside"][group] == (abs(row["z"]) <= rehearsal.Z_BAND), group
    forecast = rehearsed["forecast"]
    assert forecast["members"] >= 10 and set(forecast) >= {"truth", "covered"}
    rehearsal.pooled_recovery([rehearsed])  # a report it can read


# -- the scoring rule and the pooled verdict (ledger 593) ----------------------


class _Ensemble:
    def __init__(self, frame):
        self._df = frame


def _scored(monkeypatch, posterior, prior, truth):
    """``recovery`` over hand-made ensembles of one scalar ``p`` (group ``g``)."""

    monkeypatch.setattr(rehearsal, "truth_parameter_values", lambda case, cal: {"p": truth})
    cal = SimpleNamespace(pst=SimpleNamespace(
        parameter_data=pd.DataFrame({"parnme": ["p"], "pargp": ["g"]})))
    ies = SimpleNamespace(posterior_parameters=_Ensemble(posterior),
                          prior_parameters=_Ensemble(prior))
    return rehearsal.recovery(None, cal, ies).loc["g"]


def test_recovery_scores_the_realizations_by_one_rule(monkeypatch):
    rng = np.random.default_rng(3)
    logs = rng.normal(0.0, 0.1, 63)
    realizations = pd.DataFrame({"P": 10 ** logs}, index=[str(i) for i in range(63)])
    # A base row far from the realizations: it is not a draw and must move nothing.
    posterior = pd.concat([realizations, pd.DataFrame({"P": [100.0]}, index=["base"])])
    prior = pd.concat([realizations * 3.0, pd.DataFrame({"P": [1.0]}, index=["base"])])
    mean, sd = logs.mean(), logs.std(ddof=1)
    for z in (-1.7, -1.6, 0.4, 1.62, 1.66, 4.0):
        row = _scored(monkeypatch, posterior, prior, 10 ** (mean + z * sd))
        assert row["z"] == pytest.approx(z, abs=1e-9)
        assert row["inside"] == (abs(z) <= rehearsal.Z_BAND), z
        assert row["posterior_p05"] == pytest.approx(10 ** (mean - rehearsal.Z_BAND * sd))
        assert row["posterior_p95"] == pytest.approx(10 ** (mean + rehearsal.Z_BAND * sd))
        assert row["posterior_p05"] <= row["truth"] <= row["posterior_p95"] or not row["inside"]
    assert row["posterior_mean"] == pytest.approx(10 ** mean)
    assert row["base"] == 100.0 and row["members"] == 63, "the base row is no member"
    assert row["prior_p95"] == pytest.approx(3.0 * row["posterior_p95"])
    assert not _scored(monkeypatch, prior, posterior, 1.0)["narrowed"], "a wider posterior"


def test_recovery_z_of_a_collapsed_ensemble(monkeypatch):
    flat = pd.DataFrame({"p": [2.0] * 8}, index=range(8))
    assert _scored(monkeypatch, flat, flat * 1.5, 2.0)["z"] == 0.0
    row = _scored(monkeypatch, flat, flat * 1.5, 1.0)
    assert row["z"] == -np.inf and not row["inside"]


def test_forecast_recovery_reads_the_realizations_range():
    columns = ["oname:fore_otype:lst_usecol:f01_per:36"]
    posterior = pd.DataFrame({columns[0]: [10.0, 10.5, 11.0, 30.0]}, index=["0", "1", "2", "base"])
    prior = pd.DataFrame({columns[0]: [8.0, 12.0, 9.0, 10.0]}, index=["0", "1", "2", "base"])
    ies = SimpleNamespace(posterior=_Ensemble(posterior), prior=_Ensemble(prior),
                          forecast=lambda name: SimpleNamespace(name=columns[0]))
    inside = rehearsal.forecast_recovery(ies, pd.Series({columns[0]: 10.9}))
    assert inside == {"name": columns[0], "truth": 10.9, "posterior_low": 10.0,
                      "posterior_high": 11.0, "prior_low": 8.0, "prior_high": 12.0,
                      "members": 3, "covered": True}
    # The base realization does not widen the range.
    assert not rehearsal.forecast_recovery(ies, pd.Series({columns[0]: 20.0}))["covered"]


def _report(zs, covered=True, members=63):
    return {"scalars": {f"p{i}": {"z": float(z), "members": members} for i, z in enumerate(zs)},
            "forecast": {"covered": covered, "members": members}}


def _exact_zs(rng, runs, scalars, members):
    """The ``z`` of exact posteriors as :func:`recovery` scores them: ``members``
    realizations and one truth, all drawn from one posterior, the truth measured against
    the realizations' own mean and sd. Shape ``(runs, scalars)``."""

    draws = rng.normal(size=(runs, scalars, members + 1))
    ensemble, truth = draws[..., :members], draws[..., members]
    return (truth - ensemble.mean(axis=-1)) / ensemble.std(axis=-1, ddof=1)


def test_an_exact_posteriors_z_is_wider_than_a_normal():
    """Ledger 594: measured against ESTIMATED moments, an exact posterior's
    ``z`` is ``sqrt(1 + 1/n)`` x Student t(n - 1) -- its band holds the truth 0.892 of the
    time at 63 realizations, not 0.9, and its 1-in-3600 tail is 3.88, not a normal's 3.63.
    The pooled thresholds stand on these two numbers."""

    rng = np.random.default_rng(7)
    for members in (15, 63):
        z = np.abs(_exact_zs(rng, 4000, 50, members)).ravel()  # 200,000 values
        assert (z <= rehearsal.Z_BAND).mean() == pytest.approx(
            rehearsal.band_coverage(members), abs=0.003), members
        for edge in (2.5, 3.0):
            assert (z >= edge).mean() == pytest.approx(rehearsal.z_tail(edge, members),
                                                       rel=0.1), (members, edge)
    assert rehearsal.band_coverage(63) == pytest.approx(0.8922, abs=1e-4)
    assert rehearsal.band_coverage(15) == pytest.approx(0.8665, abs=1e-4)
    assert rehearsal.z_tail(float("nan"), 63) == 0.0
    with pytest.raises(ValueError, match="at least 2"):
        rehearsal.z_tail(1.0, 1)


def test_pooled_thresholds_and_their_false_fail_rates():
    verdict = rehearsal.pooled_recovery([_report([0.0] * 6)] * 6)
    # 36 truths, each inside with its 63-member band's coverage 0.892: fewer than 27 inside has
    # probability 0.0038, fewer than 28 0.0123 (over the 0.01 budget); all 36 |z| inside 3.884
    # with probability 0.99 (a normal null would say 3.634, ledger 594); 5 of 6
    # forecasts at 2/64 a miss.
    assert (verdict.inside_min, verdict.forecasts_min) == (27, 5)
    assert verdict.z_max == pytest.approx(3.884, abs=1e-3)
    assert verdict.false_fail["coverage"] == pytest.approx(0.00381, abs=1e-5)
    assert verdict.false_fail["forecast"] == pytest.approx(0.01347, abs=1e-5)
    assert verdict.false_fail["total"] == pytest.approx(0.0273, abs=1e-4)
    assert verdict.passed and verdict.fraction == 1.0
    one = rehearsal.pooled_recovery([_report([0.0] * 6)])
    assert (one.inside_min, one.forecasts_min) == (3, 0), "one run: 3 of 6, any forecast"
    # Fewer realizations, a wider z: at 15 the bound is 4.96.
    few = rehearsal.pooled_recovery([_report([0.0] * 6, members=15)] * 6)
    assert few.z_max == pytest.approx(4.96, abs=0.01) and few.inside_min == 26


def test_pooled_verdict_fails_each_criterion_by_name():
    good = [_report([0.0] * 6)] * 6
    nine_out = [_report([2.0] * 6), _report([2.0] * 3 + [0.0] * 3)] + good[2:]
    verdict = rehearsal.pooled_recovery(nine_out)
    assert verdict.inside == 27 and verdict.passes["coverage"]
    worse = rehearsal.pooled_recovery([_report([2.0] * 6), _report([2.0] * 4 + [0.0] * 2)]
                                      + good[2:])
    assert worse.inside == 26 and not worse.passes["coverage"] and not worse.passed
    # The as-written notebook's Sy at +3.73 (ledger 594) is inside the bound.
    near = rehearsal.pooled_recovery([_report([0.0] * 5 + [3.73])] + good[1:])
    assert near.passed and near.largest_z == pytest.approx(3.73)
    far = rehearsal.pooled_recovery([_report([0.0] * 5 + [-3.95])] + good[1:])
    assert far.passes["coverage"] and not far.passes["largest_z"]
    assert far.largest_z_at == (0, "p5") and far.largest_z == pytest.approx(3.95)
    # Judged against its OWN run's realizations: 4.5 at 15 is inside that run's 4.96,
    # 3.95 at 63 is outside its 3.88 -- the second is the one nearest failing.
    mixed = rehearsal.pooled_recovery([_report([0.0] * 5 + [4.5], members=15),
                                       _report([0.0] * 5 + [-3.95])] + good[2:])
    assert mixed.largest_z_at == (1, "p5") and not mixed.passes["largest_z"]
    assert rehearsal.pooled_recovery([_report([0.0] * 5 + [4.5], members=15)]
                                     + good[1:]).passes["largest_z"]
    missed = rehearsal.pooled_recovery([_report([0.0] * 6, covered=False)] * 2 + good[2:])
    assert missed.forecasts_covered == 4 and not missed.passes["forecast"]
    assert rehearsal.pooled_recovery([_report([float("nan")] + [0.0] * 5)] + good[1:]).inside == 35


def test_an_exact_posterior_fails_the_pooled_verdict_rarely():
    """Monte Carlo of exact, converged posteriors scored as the rehearsal scores them --
    six runs of six independent scalars, each ``z`` against 63 realizations' own mean and
    sd, a forecast missing 2/64 of the time -- fail the verdict at its stated rate, 0.027
    (the reproducer's ``correlated_false_fail`` measures 0.034 with the rehearsal's
    correlations, ledger 594). A standard normal's thresholds failed 0.046."""

    rng = np.random.default_rng(11)
    trials = 3000
    zs = _exact_zs(rng, trials * 6, 6, 63).reshape(trials, 6, 6)
    covered = rng.random((trials, 6)) >= 2 / 64
    failed = sum(
        not rehearsal.pooled_recovery([_report(zs[t, r], covered=bool(covered[t, r]))
                                       for r in range(6)]).passed
        for t in range(trials))
    stated = rehearsal.pooled_recovery([_report([0.0] * 6)] * 6).false_fail["total"]
    assert failed / trials == pytest.approx(stated, abs=0.008)
    assert failed / trials < 0.05


@pytest.mark.parametrize(("reports", "error", "match"), [
    ([], ValueError, "at least one"),
    ([{"scalars": {"p": {"z": 0.0, "members": 9}}}], KeyError, "'forecast'"),
    ([{"scalars": {}, "forecast": {"covered": True, "members": 9}}], ValueError, "no scalar"),
    ([{"scalars": {"p": {"z": 0.0}}, "forecast": {"covered": True, "members": 9}}], KeyError,
     "'members'"),
    ([{"scalars": {"p": {"z": 0.0, "members": 1}}, "forecast": {"covered": True, "members": 1}}],
     ValueError, "at least 2"),
])
def test_pooled_recovery_refuses(reports, error, match):
    with pytest.raises(error, match=match):
        rehearsal.pooled_recovery(reports)


#: The opt-in pooled coverage check (ledger 593): set to ``1`` to run it.
COVERAGE_ENV = "SIMPLE_MODFLOW_REHEARSAL_COVERAGE"
#: A converging budget (ledger 591: phi has nearly stopped falling after the third
#: iteration; 64 realizations covered the forecast on every seed) and the seeds pooled.
COVERAGE_BUDGET = "64x3+lambda"
COVERAGE_SEEDS = (1, 2, 3, 4, 5, 6)
#: PESTPP-IES workers a run: about five minutes a run at 12.
COVERAGE_WORKERS = min(12, os.cpu_count() or 1)


@pytest.mark.slow
@NEEDS_PEST
@pytest.mark.skipif(os.environ.get(COVERAGE_ENV) != "1",
                    reason=f"opt-in: half an hour of PESTPP-IES; set {COVERAGE_ENV}=1")
def test_converged_rehearsals_cover_the_truth_pooled(tmp_path, monkeypatch):
    import IPython.display

    budget = _budget_reproducer()
    reals, iterations, flags = budget.parse_config(COVERAGE_BUDGET)
    here = tmp_path / "examples" / "mf6" / "notebooks"
    here.mkdir(parents=True)
    monkeypatch.chdir(here)  # what the notebook writes beside itself lands under tmp_path
    monkeypatch.syspath_prepend(str(budget.NOTEBOOKS))
    monkeypatch.setattr(IPython.display, "display", IPython.display.display)  # run() swaps it
    reports = []
    for seed in COVERAGE_SEEDS:
        workspace = tmp_path / f"seed_{seed}"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = budget.run(reals, iterations, seed, workspace, flags=flags,
                                workers=COVERAGE_WORKERS)
        failed = sorted(name for name in MET if not report["passes"][name])
        assert failed == [], (seed, {name: report.get(name) for name in failed})
        reports.append(report)
        shutil.rmtree(workspace)  # hundreds of MB a run on a quota'd /tmp
    verdict = rehearsal.pooled_recovery(reports)
    assert verdict.passed, verdict.summary().to_string()


# -- the budget reproducer ----------------------------------------------------------------


def _budget_reproducer():
    """``docs/reproducers/rehearsal_recovery_budget.py``, imported by path."""

    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "docs" / "reproducers" / "rehearsal_recovery_budget.py"
    spec = importlib.util.spec_from_file_location("_rehearsal_recovery_budget", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_budget_reproducers_flags_still_edit_the_notebook():
    """Each ``+flag`` of the converging-budget measurement (ledger 591)
    edits text the notebook still holds, once: an edit of the notebook that moves it
    fails here, not hours into a measurement."""

    budget = _budget_reproducer()
    plain = budget._cells(())
    for flag, edits in budget.FLAGS.items():
        edited = budget._cells((flag,))
        changed = [index for index, (a, b) in enumerate(zip(plain, edited, strict=True)) if a != b]
        assert len(changed) == 1, flag
        for old, _new in edits:
            assert old in plain[changed[0]] and old not in edited[changed[0]], flag
        compile(edited[changed[0]], f"<{flag}>", "exec")
    assert budget.parse_config("64x6+lambda+nodiff") == (64, 6, ("lambda", "nodiff"))
    with pytest.raises(SystemExit, match="unknown flag"):
        budget.parse_config("64x6+lamda")
