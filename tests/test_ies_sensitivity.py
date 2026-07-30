"""Ensemble-based sensitivity on ``IesResults`` (plan §5.8 item 4, Tier A).

Tier A deliberately answers a *different question* from composite scaled
sensitivity, and does it without a jacobian — which matters because a completed
PESTPP-IES run does not produce one, so every jacobian-based measure (pyEMU's
``Schur``/``ErrVar``) needs a PESTPP-GLM run mode this library does not launch.

Everything here reads the prior and posterior ensembles the run already wrote.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest

from myflopy.modflow.mf6.pest.ies import IesResults

REALS = 40


def _stub_sensitivity_results(*, transform="none", n_reals=REALS):
    """Two parameters with known spreads and a forecast driven by one of them.

    ``informed`` has three quarters of its prior spread removed; ``idle`` has
    none. The forecast is a copy of ``informed``'s posterior draw, so its
    correlation is 1 by construction and ``idle``'s is whatever the sample noise
    of ``n_reals`` draws happens to be.
    """

    rng = np.random.default_rng(20260730)
    index = [f"r{i}" for i in range(n_reals)]
    unit = rng.normal(size=n_reals)
    unit = (unit - unit.mean()) / unit.std(ddof=1)
    other = rng.normal(size=n_reals)
    other = (other - other.mean()) / other.std(ddof=1)

    prior = pd.DataFrame({"informed": unit, "idle": other}, index=index)
    posterior = pd.DataFrame({"informed": unit * 0.25, "idle": other}, index=index)
    if transform == "log":
        prior = 10.0 ** prior
        posterior = 10.0 ** posterior

    parameter_data = pd.DataFrame(
        {
            "parnme": ["informed", "idle"],
            "pargp": ["ga", "gb"],
            "partrans": [transform, transform],
        }
    ).set_index("parnme")

    observations = pd.DataFrame(
        {"pred": posterior["informed"].to_numpy()}, index=index
    )

    class _Stub(IesResults):
        def __init__(self):
            self.pst = SimpleNamespace(parameter_data=parameter_data)

        @property
        def prior_parameters(self):
            return SimpleNamespace(_df=prior)

        @property
        def posterior_parameters(self):
            return SimpleNamespace(_df=posterior)

        @property
        def posterior(self):
            return SimpleNamespace(_df=observations)

        @property
        def forecast_names(self):
            return ["pred"]

    return _Stub()


def test_learned_is_the_share_of_prior_spread_the_data_removed():
    """`1 - posterior_sd / prior_sd`, per parameter group. A group the data did
    not constrain scores 0, not "small" -- so the number is readable as "three
    quarters of this group's uncertainty came from the data", and a group at
    zero is visibly untouched rather than merely low."""

    table = _stub_sensitivity_results().sensitivity()

    assert list(table.index) == ["ga", "gb"], "most informed group first"
    assert table.loc["ga", "learned"] == pytest.approx(0.75)
    assert table.loc["gb", "learned"] == pytest.approx(0.0, abs=1e-12)
    assert table.loc["ga", "n_reals"] == REALS


def test_log_parameters_are_scored_in_log_space():
    """A log-transformed multiplier's spread is a RATIO. Measuring it linearly
    would make a group that moved from x10 to x2 look differently informed
    depending on where its mean sat."""

    linear = _stub_sensitivity_results(transform="none").sensitivity()
    logged = _stub_sensitivity_results(transform="log").sensitivity()

    assert logged.loc["ga", "learned"] == pytest.approx(
        linear.loc["ga", "learned"], abs=1e-9
    )
    assert logged.loc["gb", "learned"] == pytest.approx(0.0, abs=1e-9)


def test_a_forecast_ranks_the_parameters_that_drive_it():
    """The 'what should I measure better?' half. The forecast here IS the
    informed parameter, so its correlation is 1 by construction; the other
    parameter's is sampling noise."""

    results = _stub_sensitivity_results()
    table = results.sensitivity(forecast="pred")

    assert table.index[0] == "ga"
    assert table.loc["ga", "forecast_corr"] == pytest.approx(1.0, abs=1e-9)
    assert abs(table.loc["gb", "forecast_corr"]) < 0.5

    with pytest.raises(ValueError, match="No forecast named"):
        results.sensitivity(forecast="nope")


def _stub_with_an_opposed_group(n_reals=REALS):
    """One group holding two parameters that push a forecast opposite ways.

    Needed because a one-parameter-per-group fixture cannot tell a signed mean
    from an absolute one -- the first version of this test passed with either.
    """

    rng = np.random.default_rng(4242)
    index = [f"r{i}" for i in range(n_reals)]
    unit = rng.normal(size=n_reals)
    unit = (unit - unit.mean()) / unit.std(ddof=1)

    # up and down are perfectly anti-correlated, and both track the forecast.
    posterior = pd.DataFrame({"up": unit, "down": -unit}, index=index)
    prior = posterior * 2.0
    parameter_data = pd.DataFrame(
        {"parnme": ["up", "down"], "pargp": ["pair", "pair"],
         "partrans": ["none", "none"]}
    ).set_index("parnme")
    observations = pd.DataFrame({"pred": unit}, index=index)

    class _Stub(IesResults):
        def __init__(self):
            self.pst = SimpleNamespace(parameter_data=parameter_data)

        @property
        def prior_parameters(self):
            return SimpleNamespace(_df=prior)

        @property
        def posterior_parameters(self):
            return SimpleNamespace(_df=posterior)

        @property
        def posterior(self):
            return SimpleNamespace(_df=observations)

        @property
        def forecast_names(self):
            return ["pred"]

    return _Stub()


def test_group_correlation_uses_absolute_values():
    """A group whose members push a forecast in OPPOSITE directions is still
    influential; averaging signed correlations cancels it to zero and hides the
    group that matters most."""

    results = _stub_with_an_opposed_group()

    per_parameter = results.sensitivity(forecast="pred", by_group=False)
    assert per_parameter.loc["up", "forecast_corr"] == pytest.approx(1.0, abs=1e-9)
    assert per_parameter.loc["down", "forecast_corr"] == pytest.approx(-1.0, abs=1e-9)
    # Signed values would average to 0.0 here; magnitudes keep the group at 1.0.
    assert sum(per_parameter["forecast_corr"]) == pytest.approx(0.0, abs=1e-9)

    grouped = results.sensitivity(forecast="pred")
    assert grouped.loc["pair", "forecast_corr"] == pytest.approx(1.0, abs=1e-9)


def test_the_noise_floor_is_drawn_on_the_figure():
    """`1/sqrt(n_reals)` is the level below which a correlation is
    indistinguishable from zero. Documenting it in a docstring is not enough --
    the reader is looking at bars, so the threshold belongs on the chart."""

    results = _stub_sensitivity_results()
    figure = results.plot_sensitivity(forecast="pred")

    assert isinstance(figure, go.Figure)
    lines = [
        shape for shape in figure.layout.shapes
        if getattr(shape, "line", None) is not None
    ]
    assert lines, "no noise-floor line drawn"
    assert lines[0].x0 == pytest.approx(1.0 / np.sqrt(REALS), rel=1e-6)


def test_fixed_and_tied_parameters_are_excluded():
    """They do not vary, so their 'spread reduction' is 0/0 -- a NaN that would
    sort into the table as if it were a real, uninformed group."""

    results = _stub_sensitivity_results()
    results.pst.parameter_data.loc["idle", "partrans"] = "fixed"

    table = results.sensitivity()
    assert list(table.index) == ["ga"]


def test_plot_sensitivity_supports_both_backends():
    results = _stub_sensitivity_results()
    assert isinstance(results.plot_sensitivity(), go.Figure)

    static = results.plot_sensitivity(backend="matplotlib")
    assert static.axes, "matplotlib backend returned no axes"
