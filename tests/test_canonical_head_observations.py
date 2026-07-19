"""The canonical model's synthetic head observations.

Reported 2026-07-18: ``model.targets.heads.calibration_plot()`` drew empty axes
because every head target carried ``head: np.nan`` -- the model is synthetic, so
there were no measured values to compare against, and ``stats()``/``residuals()``
were vacuous for the same reason.

The targets now sample the **regional water table** (the conceptual surface that
also seeds initial conditions), so the calibration tier is exercised end to end.
Contracts pinned here:

1. Measured values exist and pair with simulated ones, so compare/stats/plot are
   all non-empty.
2. They are NOT copies of the simulated heads -- residuals must be real, or the
   calibration plot is a tautology that would pass while teaching nothing.
3. Pumping wells are observed once, pre-pumping. A regional survey value inside
   an actively drawn-down well is not a measurement, and carrying it forward
   would plant an outlier that dominates every residual statistic.
4. The values are deterministic across runs (the notebooks must not drift).
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def compare(canonical_run):
    """The canonical model's observed-vs-simulated head table."""

    return canonical_run.targets.heads.compare()


def test_head_targets_carry_measured_values(compare):
    """The defect: ``head_target`` used to be NaN for every row."""

    assert compare["head_target"].notna().any(), (
        "head targets carry no measured values -- calibration_plot(), stats() "
        "and residuals() are all vacuous when this regresses"
    )


def test_observed_and_simulated_actually_pair(compare):
    """Points survive the dropna inside the cross plot."""

    paired = compare.dropna(subset=["head_target", "sim_head"])
    assert len(paired) >= 8, f"only {len(paired)} paired points"


def test_stats_are_computed_not_nan(canonical_run):
    stats = canonical_run.targets.heads.stats()
    assert not stats.empty
    row = stats.iloc[0]
    assert int(row["n"]) >= 8
    for column in ("mean_error", "mae", "rmse"):
        assert np.isfinite(float(row[column])), column


def test_calibration_plot_has_points_and_no_empty_warning(canonical_run):
    """The reported symptom, inverted: a populated plot with no diagnosis box."""

    import warnings

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        figure = canonical_run.targets.heads.calibration_plot()

    empty_warnings = [
        str(record.message)
        for record in caught
        if issubclass(record.category, UserWarning)
        and "no points to plot" in str(record.message)
    ]
    assert not empty_warnings, empty_warnings

    model_trace = figure.data[0]
    assert model_trace.x is not None and len(model_trace.x) >= 8


def test_observations_are_not_copies_of_the_simulated_heads(compare):
    """Residuals must be real -- otherwise the plot is a tautology.

    Sampling the regional surface (rather than the model's own output) is the
    whole point: the simulated heads depart from the conceptual water table
    exactly where the stresses bite, and that departure is the lesson.
    """

    paired = compare.dropna(subset=["head_target", "sim_head"])
    residual = (paired["sim_head"] - paired["head_target"]).abs()
    assert residual.max() > 0.5, "observations look copied from simulated heads"
    # ...but not so far off that the demo looks broken
    assert residual.max() < 25.0, f"residuals implausibly large: {residual.max()}"


def test_production_wells_carry_no_measured_value(compare):
    """A survey does not report a head measured inside a pumping well.

    They stay registered as observation LOCATIONS -- their simulated series and
    head-change signals are still used -- but the measured side is absent.
    Including them planted +11 ft and +23 ft outliers that dominated every
    residual statistic.
    """

    for name in ("shallow_pumping", "deep_pumping"):
        rows = compare[compare["name"] == name]
        assert not rows.empty, f"{name} must remain an observation location"
        assert rows["head_target"].isna().all(), f"{name} must carry no measured value"
        assert rows["sim_head"].notna().any(), f"{name} must still simulate"


def test_monitoring_wells_are_observed_every_period(compare, canonical_run):
    """The regional network carries the full series."""

    monitors = [name for name in compare["name"].unique() if str(name).startswith("monitor_")]
    assert len(monitors) >= 6, f"thin monitoring network: {monitors}"
    for name in [*monitors, "regional_center", "pond_mound"]:
        rows = compare[compare["name"] == name].dropna(subset=["head_target"])
        assert len(rows) == canonical_run.nper, name


# ---------------------------------------------------------------------------
# the user requirement: the fast tour must show a WELL-calibrated model
# ---------------------------------------------------------------------------


def test_the_model_reads_as_well_calibrated(canonical_run):
    """RMSE under 5% of the observed head range, and essentially unbiased.

    This is the contract behind the fast tour's calibration scatter (ledger
    entry 25). 5% of range is the usual "good fit" benchmark; the mean-error
    bound is what stops a surface that is merely *precise* while sitting
    systematically high or low, which is what the first attempt did (RMSE 6.6 ft
    = 20% of range, with every well biased the same direction).
    """

    stats = canonical_run.targets.heads.stats().iloc[0]
    paired = canonical_run.targets.heads.compare().dropna(
        subset=["head_target", "sim_head"]
    )
    simulated = paired["sim_head"].astype(float)
    head_range = float(simulated.max() - simulated.min())

    rmse = float(stats["rmse"])
    assert head_range > 20.0, "head range too small for the ratio to mean anything"
    assert rmse / head_range < 0.05, (
        f"RMSE {rmse:.2f} ft is {100 * rmse / head_range:.1f}% of the "
        f"{head_range:.1f} ft head range; the fast tour must show a "
        "well-calibrated model (ledger entry 25)"
    )
    assert abs(float(stats["mean_error"])) < 0.75, (
        f"mean error {float(stats['mean_error']):+.2f} ft -- the observation "
        "surface has drifted into a systematic bias"
    )


def test_observations_are_deterministic(canonical_config):
    """No global RNG: rebuilding must reproduce the values exactly."""

    from myflopy.modflow.mf6.canonical_example import _synthetic_head_observations

    xn = np.linspace(0.0, 1.0, 400)
    relief = np.zeros(400)
    locations = {"a": (0, 10), "b": (2, 250)}

    first = _synthetic_head_observations(canonical_config, locations, xn, relief)
    second = _synthetic_head_observations(canonical_config, locations, xn, relief)
    assert first.equals(second)


def test_the_observed_surface_declines_down_valley_and_with_depth():
    """The fitted surface must keep the physical shape it stands in for."""

    from myflopy.modflow.mf6.canonical_example import _observed_water_table

    xn = np.linspace(0.05, 0.95, 25)
    surface = _observed_water_table(xn, np.zeros_like(xn), 0)
    assert surface[0] > surface[-1] + 30.0, "no down-valley decline"
    assert np.all(np.diff(surface) < 0.0), "the water table must fall monotonically"

    deeper = _observed_water_table(xn, np.zeros_like(xn), 3)
    assert np.all(deeper < surface), "head must decline with depth"


def test_monitoring_wells_avoid_every_stress_cell(canonical_run):
    """A monitor reading seepage or drawdown is not a regional baseline."""

    locations = canonical_run.targets.heads.get()
    monitors = locations[locations["name"].astype(str).str.startswith("monitor_")]
    assert not monitors.empty

    model = canonical_run
    stressed = set()
    for record in model.gwf.lak.connectiondata.get_data():
        stressed.add(int(record[2][-1]))
    for record in model.gwf.sfr.packagedata.get_data():
        stressed.add(int(record[1][-1]))
    for package in ("riv", "wel", "drn"):
        data = getattr(model.gwf, package).stress_period_data.get_data(0)
        for record in data if data is not None else []:
            stressed.add(int(record[0][-1]))

    overlap = sorted(set(monitors["cell"].astype(int)) & stressed)
    assert not overlap, f"monitoring wells sit on stress cells: {overlap}"
