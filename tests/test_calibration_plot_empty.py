"""An empty calibration plot must say why, not render blank styled axes.

Reported 2026-07-18: ``model.targets.heads.calibration_plot()`` on the canonical
model drew a fully formatted cross plot with no points. The cause was upstream of
the plot -- the canonical head targets carry ``head: np.nan``, so the
``dropna(subset=[target, simulated])`` inside ``from_obs_vs_sim`` emptied the
frame. The figure had no way to say so, and the caller was left guessing.

These pin the diagnosis: which side is missing, in the caller's own column
names, both as a ``UserWarning`` (catchable in scripts) and as an on-figure
annotation (unmissable in a notebook). Healthy data must stay silent -- a plot
that cried wolf would be worse than the blank one.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.calcs.calibration import CalibrationPlot


def _plot(frame, **kwargs):
    """Build the cross plot, returning ``(figure, user-warning messages)``."""

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        figure = CalibrationPlot.from_obs_vs_sim(frame, **kwargs)
    messages = [
        str(record.message)
        for record in caught
        if issubclass(record.category, UserWarning)
    ]
    return figure, messages


def _annotations(figure):
    return [note.text for note in figure.layout.annotations]


def _point_counts(figure):
    return [0 if trace.x is None else len(trace.x) for trace in figure.data]


# ---------------------------------------------------------------------------
# the four ways a cross plot can come up empty
# ---------------------------------------------------------------------------


def test_all_nan_observations_name_the_target_column():
    """The reported case: targets exist but carry no measured values."""

    figure, messages = _plot(
        pd.DataFrame({"head_target": [np.nan] * 5, "sim_head": [1.0, 2, 3, 4, 5]})
    )

    assert len(messages) == 1
    assert "'head_target' is entirely NaN across 5 rows" in messages[0]
    assert "no measured values" in messages[0]
    assert any("No paired observed/simulated values" in text for text in _annotations(figure))
    assert _point_counts(figure) == [0]


def test_all_nan_simulated_points_at_the_model():
    """The opposite failure reads as a model problem, not a data problem."""

    _, messages = _plot(
        pd.DataFrame({"head_target": [1.0, 2, 3], "sim_head": [np.nan] * 3})
    )

    assert "'sim_head' is entirely NaN across 3 rows" in messages[0]
    assert "was the model run?" in messages[0]


def test_non_overlapping_rows_are_named_as_a_join_problem():
    """Both columns have values, but never on the same row."""

    _, messages = _plot(
        pd.DataFrame({"head_target": [1.0, np.nan], "sim_head": [np.nan, 2.0]})
    )

    assert "no row has BOTH" in messages[0]
    assert "1 observed and 1 simulated" in messages[0]
    assert "did not line up" in messages[0]


def test_a_completely_empty_compare_table_is_named_as_such():
    _, messages = _plot(pd.DataFrame({"head_target": [], "sim_head": []}))

    assert "the compare table is empty" in messages[0]


def test_both_columns_nan_is_reported_together():
    _, messages = _plot(
        pd.DataFrame({"head_target": [np.nan] * 3, "sim_head": [np.nan] * 3})
    )

    assert "both 'head_target' and 'sim_head' are entirely NaN" in messages[0]


# ---------------------------------------------------------------------------
# the diagnosis must not fire on healthy data
# ---------------------------------------------------------------------------


def test_healthy_data_plots_silently():
    """No warning, no annotation, and the points are actually there."""

    figure, messages = _plot(
        pd.DataFrame({"head_target": [1.0, 2, 3], "sim_head": [1.1, 2.1, 2.9]})
    )

    assert messages == []
    assert _annotations(figure) == []
    assert _point_counts(figure)[0] == 3


def test_partially_missing_data_plots_the_rows_that_pair():
    """Some NaN is normal; only a total absence of pairs is a problem."""

    figure, messages = _plot(
        pd.DataFrame(
            {"head_target": [1.0, np.nan, 3.0], "sim_head": [1.1, 2.1, 2.9]}
        )
    )

    assert messages == []
    assert _point_counts(figure)[0] == 2


# ---------------------------------------------------------------------------
# both backends
# ---------------------------------------------------------------------------


def test_the_matplotlib_backend_annotates_too():
    """The static backend must not be the one that stays silent."""

    pytest.importorskip("seaborn")
    frame = pd.DataFrame({"head_target": [np.nan] * 4, "sim_head": [1.0, 2, 3, 4]})

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        figure = CalibrationPlot.from_obs_vs_sim(frame, backend="matplotlib")

    assert any(issubclass(record.category, UserWarning) for record in caught)
    texts = [text.get_text() for axis in figure.axes for text in axis.texts]
    assert any("No paired observed/simulated values" in text for text in texts)


def test_custom_column_names_appear_in_the_diagnosis():
    """The message must use the caller's column names, not hardcoded ones."""

    _, messages = _plot(
        pd.DataFrame({"observed_q": [np.nan] * 3, "modelled_q": [1.0, 2, 3]}),
        target_column="observed_q",
        simulated_column="modelled_q",
    )

    assert "'observed_q' is entirely NaN" in messages[0]
