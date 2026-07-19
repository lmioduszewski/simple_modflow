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


def test_pumping_wells_are_observed_only_before_pumping(compare):
    """A regional value inside a drawn-down well is not a measurement."""

    for name in ("shallow_pumping", "deep_pumping"):
        rows = compare[compare["name"] == name]
        observed = rows.dropna(subset=["head_target"])
        assert set(observed["time"]) == {0}, (
            f"{name} should be observed only at period 0, got {sorted(set(observed['time']))}"
        )


def test_unpumped_wells_are_observed_every_period(compare, canonical_run):
    """The monitoring wells carry the full series."""

    for name in ("regional_center", "pond_mound"):
        rows = compare[compare["name"] == name].dropna(subset=["head_target"])
        assert len(rows) == canonical_run.nper


def test_observations_are_deterministic(canonical_config):
    """No global RNG: rebuilding must reproduce the values exactly."""

    from myflopy.modflow.mf6.canonical_example import _synthetic_head_observations

    regional = np.linspace(150.0, 90.0, 400)
    locations = {"a": (0, 10), "b": (2, 250)}

    first = _synthetic_head_observations(canonical_config, locations, regional)
    second = _synthetic_head_observations(canonical_config, locations, regional)
    assert first.equals(second)


def test_layer_decline_matches_initial_conditions():
    """The sampled value follows the regional table down the column."""

    from myflopy.modflow.mf6.canonical_example import _synthetic_head_observations

    class _Config:
        nper = 1

    regional = np.full(300, 120.0)
    # Same well name in both, so the per-well survey offset cancels and only
    # the layer term is under test.
    shallow = _synthetic_head_observations(_Config(), {"w": (0, 5)}, regional)
    deep = _synthetic_head_observations(_Config(), {"w": (3, 5)}, regional)

    drop = float(shallow["head"].iloc[0]) - float(deep["head"].iloc[0])
    assert drop == pytest.approx(3 * 0.6, abs=1e-9)
