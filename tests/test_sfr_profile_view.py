"""The SFR long-profile view -- the first package view on the house shape.

``model.packages.sfr.results.profile`` is a view object rather than a pair of
loose ``long_profile`` / ``plot_long_profile`` methods (see
``docs/view_layer_conventions.md``). Contracts pinned here:

1. The view's shape: ``get`` / ``plot`` / ``summary``, and ``__call__``
   rebinding the stress period so the two spellings agree.
2. ``signed_exchange`` colors gaining reaches BLUE and losing reaches RED,
   taken from the same scale the SFR maps use -- the convention documented in
   ``docs/mf6io_reference.md`` and plan §1, and the one the old hand-rolled
   notebook cell got wrong (it used green for gaining).
3. The deprecated ``long_profile`` / ``plot_long_profile`` spellings still
   return exactly what they returned before, and stay hidden per D12.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.mf6.package_plotting import (
    _blue_white_red_diverging_colorscale,
)
from myflopy.modflow.mf6.package_surface_water import (
    SfrProfileView,
    SfrResultsNamespace,
)

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def profile(canonical_run):
    """The canonical model's SFR profile view."""

    return canonical_run.packages.sfr.results.profile


# ---------------------------------------------------------------------------
# 1. the view's shape
# ---------------------------------------------------------------------------


def test_profile_is_a_view_not_a_frame(profile):
    """``results.profile`` is the view; the frame comes from ``.get()``."""

    assert isinstance(profile, SfrProfileView)
    assert isinstance(profile.get(), pd.DataFrame)


def test_profile_frame_carries_the_merged_fields(profile):
    """The table merges geometry, streambed, stage and exchange by reach."""

    frame = profile.get()
    assert not frame.empty
    for column in ("reach", "cell", "distance_mid", "streambed_top", "stage", "q"):
        assert column in frame.columns, column
    # one row per reach -- the merge must not fan out
    assert frame["reach"].is_unique


def test_calling_the_view_rebinds_the_period(canonical_run):
    """``profile(per=n).get()`` equals ``profile.get(per=n)``."""

    results = canonical_run.packages.sfr.results
    last = canonical_run.nper - 1
    rebound = results.profile(per=last)

    assert rebound.per == last
    assert results.profile.per == 0, "rebinding must not mutate the original view"
    pd.testing.assert_frame_equal(rebound.get(), results.profile.get(per=last))


def test_plot_honors_both_period_spellings(canonical_run):
    """``profile.plot(per=n)`` and ``profile(per=n).plot()`` agree."""

    results = canonical_run.packages.sfr.results
    last = canonical_run.nper - 1
    assert (
        results.profile.plot(per=last).layout.title.text
        == results.profile(per=last).plot().layout.title.text
        == f"SFR long profile (per={last})"
    )


def test_summary_digests_the_profile_fields(profile):
    """``.summary()`` covers the fields worth digesting, and nothing empty."""

    summary = profile.summary()
    assert not summary.empty


def test_plot_returns_the_house_figure_class(profile):
    """The figure must be myflopy's ``Fig`` -- scrollZoom, pan, house template.

    This is the defect that started this work: the notebook hand-rolled
    matplotlib and silently gave up all of it.
    """

    from myflopy import viz as figs

    assert isinstance(profile.plot(), figs.Fig)


# ---------------------------------------------------------------------------
# 2. the signed-exchange convention
# ---------------------------------------------------------------------------


def _exchange_trace(fig):
    """Return the figure's exchange trace."""

    matches = [trace for trace in fig.data if "Exchange q" in str(trace.name)]
    assert len(matches) == 1, [trace.name for trace in fig.data]
    return matches[0]


def test_signed_exchange_colors_gaining_blue_and_losing_red(profile):
    """Blue where the reach gains (q<0), red where it loses (q>0).

    SFR keeps MF6's RAW sign: the cell record is flow FROM the reach TO the
    aquifer (the "gwf" frame), so NEGATIVE q means the reach gains. myflopy does
    not normalize it (``docs/mf6io_reference.md``); the frame is declared, and
    the colours orient to it.
    """

    frame = profile.get()
    q = pd.to_numeric(frame["q"], errors="coerce").to_numpy(float)
    assert (q < 0).any() and (q > 0).any(), (
        "the canonical stream must both gain and lose for this test to bite"
    )

    gaining, losing = _blue_white_red_diverging_colorscale()[0][1], (
        _blue_white_red_diverging_colorscale()[-1][1]
    )
    colors = np.asarray(_exchange_trace(profile.plot()).marker.color)
    np.testing.assert_array_equal(colors, np.where(q < 0.0, gaining, losing))


def test_signed_bar_colors_track_the_shared_map_scale(profile):
    """The bars read their colors off the SFR map's scale, so they cannot drift."""

    scale = _blue_white_red_diverging_colorscale()
    colors = set(np.asarray(_exchange_trace(profile.plot()).marker.color).tolist())
    assert colors <= {scale[0][1], scale[-1][1]}


def test_signed_exchange_draws_bars_and_unsigned_draws_a_line(profile):
    """The two modes differ in trace type, not just color."""

    assert _exchange_trace(profile.plot(signed_exchange=True)).type == "bar"
    assert _exchange_trace(profile.plot(signed_exchange=False)).type == "scattergl"


def test_exchange_rides_the_secondary_axis(profile):
    """Exchange is a flux against an elevation axis -- it needs its own scale."""

    assert _exchange_trace(profile.plot()).yaxis == "y2"
    assert profile.plot().layout.yaxis2.overlaying == "y"


def test_include_flags_drop_their_traces(profile):
    """Each ``include_*`` flag removes exactly its own trace."""

    names = lambda fig: {str(trace.name) for trace in fig.data}  # noqa: E731

    assert not any("Exchange" in n for n in names(profile.plot(include_exchange=False)))
    assert not any("Stage" == n for n in names(profile.plot(include_stage=False)))
    assert not any(
        "Streambed" in n for n in names(profile.plot(include_streambed=False))
    )


# ---------------------------------------------------------------------------
# 3. the deprecated spellings (D12)
# ---------------------------------------------------------------------------


def test_long_profile_alias_warns_and_returns_the_frame(canonical_run):
    """The old name still returns a DataFrame, identical to ``profile.get()``."""

    results = canonical_run.packages.sfr.results
    with pytest.warns(DeprecationWarning, match="long_profile is deprecated"):
        frame = results.long_profile(per=0)
    pd.testing.assert_frame_equal(frame, results.profile.get(per=0))


def test_plot_long_profile_alias_keeps_the_unsigned_line(canonical_run):
    """The old plot name must keep drawing the line it always drew.

    New callers get signed bars by default; the alias must not change shape
    under existing callers.
    """

    results = canonical_run.packages.sfr.results
    with pytest.warns(DeprecationWarning, match="plot_long_profile is deprecated"):
        fig = results.plot_long_profile(per=0)
    assert _exchange_trace(fig).type == "scattergl"


def test_deprecated_spellings_stay_out_of_completion():
    """D12: warned aliases resolve only via ``__getattr__``."""

    for name in ("long_profile", "plot_long_profile"):
        assert name not in dir(SfrResultsNamespace)
        assert not hasattr(SfrResultsNamespace, name)


def test_retired_spellings_are_invisible_to_ide_completion():
    """D12, the IDE-facing half: nothing an IDE reads may echo a retired name.

    myflopy ships ``py.typed`` and no stubs, so editors read this source
    directly. That makes four surfaces IDE-visible: ``dir()``, ``__all__``,
    every attribute in the class dict (**including** underscore members, which
    editors do offer), and every docstring shown in a tooltip. A retired name
    surviving in any of them re-creates exactly the confusion D12 removes.
    """

    import myflopy.modflow.mf6.package_surface_water as module

    for retired in ("long_profile", "plot_long_profile"):
        assert retired not in module.__all__
        for cls in (SfrResultsNamespace, SfrProfileView):
            assert retired not in dir(cls)

            echoing = [name for name in vars(cls) if retired in name]
            assert not echoing, f"{cls.__name__} attribute names echo {retired!r}: {echoing}"

            tooltips = {"<class>": cls.__doc__ or ""}
            for name in dir(cls):
                member = getattr(cls, name, None)
                tooltips[name] = getattr(member, "__doc__", "") or ""
            leaking = [name for name, doc in tooltips.items() if retired in doc]
            assert not leaking, f"{cls.__name__} docstrings mention {retired!r}: {leaking}"


def test_unknown_attributes_still_raise_attribute_error(canonical_run):
    """The alias ``__getattr__`` must not swallow genuine typos."""

    with pytest.raises(AttributeError, match="no attribute 'no_such_thing'"):
        canonical_run.packages.sfr.results.no_such_thing
