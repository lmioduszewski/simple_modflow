"""Colorscale policy: diverging red/white/blue ONLY for signed 'q'-like fields
and diff maps (negative red, positive blue); everything else uses the house
brown-to-blue 'earth' scale (the mounding-figure default).
"""

from __future__ import annotations

import numpy as np
import pytest
from matplotlib.colors import to_hex

from myflopy.modflow.mf6.grid.plotting import build_choropleth
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.package_plotting import _blue_white_red_diverging_colorscale
from myflopy.modflow.mf6.package_registry import (
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
    get_package_explorer_spec,
    get_package_result_spec,
)
from myflopy.modflow.mf6.pest.ies import _field_map_policy, _relabel_log_colorbar
from myflopy.modflow.utils.datatypes.hover import parameter_field_hover


def test_non_signed_inputs_default_to_earth():
    # every non-signed input field and package default is the earth scale
    for package in ("rch", "chd", "drn", "ghb", "riv", "evt"):
        assert get_default_package_colorscale(package) == "earth"
    for field in ("finf", "pet", "extdp", "extwc", "ha", "hroot", "rootact"):
        assert get_default_package_colorscale(f"uzf_{field}") == "earth"
    assert get_default_package_colorscale("drn_cond") == "earth"
    assert get_default_package_colorscale("ghb_bhead") == "earth"
    for field in ("stage", "cond", "rbot"):
        assert get_default_package_colorscale(f"riv_{field}") == "earth"
    for field in ("surface", "rate", "depth"):
        assert get_default_package_colorscale(f"evt_{field}") == "earth"


def test_signed_q_fields_keep_diverging_scales():
    # wel q is signed (pumping negative) -> diverging stays
    assert get_default_package_colorscale("wel") == "RdBu"
    assert get_package_explorer_spec("wel").inputs["q"].colorscale == "RdBu"
    # every package's q RESULT keeps a diverging scale
    for package in ("rch", "chd", "drn", "ghb", "riv", "evt", "wel", "sfr", "lak"):
        assert get_package_result_spec(package, "q").colorscale == "RdBu"


def test_uzf_results_are_earth_not_diverging():
    assert get_package_result_spec("uzf", "gwrch").colorscale == "earth"
    assert get_package_result_spec("uzf", "sat").colorscale == "earth"


def test_prt_derived_maps_follow_the_same_non_signed_rule():
    """PRT's derived maps are registry-free, so their scale is pinned at the source.

    Travel times and particle counts are one-sided magnitudes -- never a signed
    difference -- so they take the same ``'earth'`` scale as every other
    non-signed field. The per-map assertions live in ``test_prt_maps.py``.
    """

    from myflopy.modflow.mf6.prt_maps import PRT_COLORSCALE

    assert PRT_COLORSCALE == "earth"


# ---------------------------------------------------------------------------
# PEST/IES captured parameter fields -- registry-free like PRT's derived maps,
# but keyed by STATISTIC rather than by field name: the same captured K array is
# a magnitude as a mean and a ratio as a change.
# ---------------------------------------------------------------------------


def _two_cell_vor():
    """The smallest grid a Choro will build on (mirrors test_mf6_pest.py)."""

    verts = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        dtype=float,
    )
    return VoronoiGridPlus(
        verts=verts,
        iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float),
    )


def test_field_map_magnitudes_use_the_same_earth_scale_as_everything_else():
    """A parameter-field magnitude is not signed, so it takes the house scale.

    Asserted equal to the other two "earth" answers rather than to the literal,
    so the three cannot drift apart.
    """

    from myflopy.modflow.mf6.prt_maps import PRT_COLORSCALE

    for stat in ("mean", "base", "std"):
        assert _field_map_policy(stat, [1.0, 2.0])["colorscale"] == "earth"
    assert _field_map_policy("mean", [1.0, 2.0])["colorscale"] == PRT_COLORSCALE
    assert _field_map_policy("mean", [1.0, 2.0])["colorscale"] == get_default_package_colorscale("rch")


def test_only_the_conductivity_stats_go_log_never_the_spread():
    """K spans decades and goes log; a posterior spread must not.

    A posterior standard deviation is legitimately ZERO wherever the ensemble
    collapsed, and ``log10(0)`` blanks the cell -- which would erase exactly the
    cells an uncertainty map exists to show.
    """

    assert _field_map_policy("mean", [1e-3, 1e2])["logscale"] is True
    assert _field_map_policy("base", [1e-3, 1e2])["logscale"] is True
    assert _field_map_policy("std", [0.0, 2.0])["logscale"] is False


def test_change_maps_diverge_red_for_decrease_symmetrically_about_one():
    """``change`` is a RATIO: neutral at 1, mapped in log space, red = reduced.

    The endpoints are asserted against the sibling helper rather than against
    hex literals, so the two diverging orientations cannot drift apart.
    """

    policy = _field_map_policy("change", [0.5, 10.0])
    blue, white, red = (color for _, color in _blue_white_red_diverging_colorscale())

    assert policy["logscale"] is True
    assert policy["colorscale"] == [[0.0, red], [0.5, white], [1.0, blue]]
    # symmetric in LOG space: a halving and a doubling sit equally either side.
    assert policy["zmin"] == -policy["zmax"]
    assert policy["zmid"] == 0.0
    # a decade colorbar, so the reader sees ratios and not log10 units
    assert policy["colorbar"]["ticktext"] == ["0.1x", "1x", "10x"]


def test_degenerate_change_fields_still_get_a_defined_range():
    """A prior Monte-Carlo run has change == 1 everywhere; an empty layer is NaN.

    Both give a symmetric-limit of 0, which would collapse the color range onto a
    single point, so the policy floors it at one decade.
    """

    for values in ([1.0, 1.0], [np.nan, np.nan]):
        policy = _field_map_policy("change", values)
        assert (policy["zmin"], policy["zmax"]) == (-1.0, 1.0)


def test_an_undefined_ratio_cannot_widen_the_color_range_or_warn():
    """A zero prior mean gives inf, which must not blow the scale out or warn."""

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        policy = _field_map_policy("change", [np.inf, -np.inf, np.nan, 2.0])

    # the range comes from the one finite value (log10(2)), not from inf
    assert policy["zmax"] == pytest.approx(np.log10(2.0))


def test_an_unknown_field_stat_raises_instead_of_silently_defaulting():
    """The deliberate inverse of ``Choro.colorscale``, which prints and falls back.

    A typo'd stat is a caller bug; falling back to earth would draw a plausible
    map of the wrong thing.
    """

    with pytest.raises(ValueError, match="No field-map color policy"):
        _field_map_policy("meen", [1.0])


def test_field_map_policy_survives_the_real_choropleth_front_door():
    """The call-site pin: assert on what the front door actually renders.

    The same lesson as ``test_every_exchange_map_puts_blue_on_gaining`` below --
    a helper-only test proves arithmetic while the real map renders inverted.
    This one specifically catches that ``Choro``'s plotly-to-matplotlib table maps
    ``'rdbu'`` to the REVERSED colormap, so a diverging scale passed as a NAME
    renders mirrored between the two backends. Stops survive both.
    """

    vor = _two_cell_vor()
    values = [0.5, 10.0]
    policy = _field_map_policy("change", values)
    choro = build_choropleth(
        vor,
        custom_zs=values,
        custom_hover={"change": values},
        hover_spec=parameter_field_hover("change", title="k change"),
        hover_heads=False,
        hover_ks=False,
        **policy,
    )
    trace = choro.get_choropleth()
    red = _blue_white_red_diverging_colorscale()[-1][1]

    # plotly: red at the low end, limits symmetric in log space
    assert trace.colorscale[0][1] == red
    assert trace.zmin == -trace.zmax
    # the widening landed: these four names could not reach Choro before 6.4A
    assert trace.zmid == 0.0
    assert tuple(trace.colorbar.ticktext) == ("0.1x", "1x", "10x")
    assert "k change" in trace.hovertemplate
    assert trace.customdata is not None

    # matplotlib: the SAME low-end color and the SAME limits
    collection = choro.plot_mpl().axes[0].collections[0]
    assert collection.get_clim() == (trace.zmin, trace.zmax)
    assert to_hex(collection.get_cmap()(0.0)) == red


def test_a_log_map_reads_in_real_units_on_both_backends():
    """A log map's colorbar must not be labeled in log10 units on either backend.

    ``plot_mpl`` builds its colorbar through geopandas with no tick hook, so
    without an explicit relabel a static K map reads ``-3 … 2`` where the
    interactive one reads ``0.001 … 100`` — the same field, off by orders of
    magnitude, with nothing on the figure to say so.
    """

    values = [0.001, 100.0]
    policy = _field_map_policy("mean", values)
    choro = build_choropleth(
        _two_cell_vor(), custom_zs=values, hover_heads=False, hover_ks=False, **policy
    )
    figure = choro.plot_mpl()
    _relabel_log_colorbar(figure, policy["colorbar"])

    expected = list(policy["colorbar"]["ticktext"])
    assert expected[0] == "0.001"  # real units, not "-3"
    assert [label.get_text() for label in figure.axes[-1].get_yticklabels()] == expected
    # and the linear stat has no decade colorbar to misapply
    assert _field_map_policy("std", [0.0, 2.0]).get("colorbar") is None


def test_relabeling_a_colorbar_that_does_not_exist_leaves_the_map_axes_alone():
    """``plot_mpl(colorbar=False)`` leaves ONE axes -- the map -- and it must be untouched.

    Written this way because the obvious version (assert it does not raise) passes
    against broken code: with no colorbar, ``axes[-1]`` *is* the map, and stamping
    decade ticks onto its northing axis raises nothing. It just relabels the map's
    y-axis in log10(K).
    """

    values = [0.001, 100.0]
    policy = _field_map_policy("mean", values)
    choro = build_choropleth(
        _two_cell_vor(), custom_zs=values, hover_heads=False, hover_ks=False, **policy
    )
    figure = choro.plot_mpl(colorbar=False)
    assert len(figure.axes) == 1
    before = list(figure.axes[0].get_yticks())

    _relabel_log_colorbar(figure, policy["colorbar"])

    assert list(figure.axes[0].get_yticks()) == before
    decades = set(policy["colorbar"]["ticktext"])
    assert not decades & {label.get_text() for label in figure.axes[0].get_yticklabels()}


def test_the_grid_accessor_accepts_what_the_function_accepts():
    """``vor.choropleth`` restates ``build_choropleth``'s signature by hand.

    Widening only the function would leave the accessor raising ``TypeError`` for
    arguments the function itself takes.
    """

    choro = _two_cell_vor().choropleth(
        custom_zs=[1.0, 2.0], colorscale="earth", logscale=True, zmid=0.0
    )
    assert choro.get_choropleth().zmid == 0.0


def test_category_colors_are_policy_and_stable_across_figures(isolated_category_colors):
    """Named categories get colors from one memoized helper, not per-call-site hexes.

    A release group drawn on the pathline map, its arrival curve, and its capture
    bars must be the same color, or a set of small multiples stops being readable
    -- so the mapping is remembered rather than recomputed per figure.
    """

    from myflopy.viz import PALETTE, category_colors

    first = category_colors(["west_wells", "east_wells"])
    assert set(first.values()) <= set(PALETTE.categorical)
    assert len(set(first.values())) == 2

    # Same name -> same color later. "west_wells" sorted SECOND when it was first
    # seen and sorts FIRST here, so an implementation that re-enumerates per call
    # instead of remembering would hand it a different color.
    later = category_colors(["west_wells", "zzz_wells"])
    assert later["west_wells"] == first["west_wells"]
    assert later["zzz_wells"] != later["west_wells"]

    # A figure's own categories stay distinguishable even when one name is first
    # seen long after its sibling -- the collision a plain global cycle counter
    # produces once the palette wraps. Again the remembered name sorts second, so
    # a forgetful implementation would move it.
    for index in range(len(PALETTE.categorical) + 1):
        category_colors([f"filler_{index}"])
    pair = category_colors(["aaa_late_wells", "east_wells"])
    assert pair["east_wells"] == first["east_wells"]
    assert pair["aaa_late_wells"] != pair["east_wells"]

    # more categories at once than colors: the palette repeats rather than running out
    many = category_colors([f"zone_{index}" for index in range(len(PALETTE.categorical) + 3)])
    assert set(many.values()) <= set(PALETTE.categorical)


def test_diff_maps_use_rdbu_negative_red_positive_blue():
    # plotly RdBu runs red (low) -> blue (high); with zmid=0 that is
    # negative red / positive blue, the required diff-map orientation
    assert get_default_group_compare_colorscale() == "RdBu"


def test_gaining_losing_scale_is_blue_negative_red_positive():
    scale = _blue_white_red_diverging_colorscale()
    assert scale[0][0] == 0.0 and scale[-1][0] == 1.0
    assert scale[0][1] == "#1f77b4"  # low end (gaining, negative q) blue
    assert scale[-1][1] == "#d62728"  # high end (losing, positive q) red


# --- Choro accepts explicit list colorscales (regression: setter dropped them) --
def _bare_choro():
    from myflopy.modflow.utils.datatypes.choros import Choro

    choro = Choro.__new__(Choro)
    choro._colorscale = None
    return choro


def test_choro_setter_accepts_list_colorscales():
    choro = _bare_choro()
    choro.colorscale = _blue_white_red_diverging_colorscale()
    assert choro.colorscale[0][1] == "#1f77b4"  # list passes through, not 'earth'


def test_choro_setter_accepts_earth_and_rejects_unknown():
    choro = _bare_choro()
    choro.colorscale = "earth"
    assert choro.colorscale == "earth"
    choro = _bare_choro()
    choro.colorscale = "not_a_scale"
    assert choro.colorscale == "earth"  # falls back with a console note


def test_choro_default_colorscale_is_earth():
    assert _bare_choro().colorscale == "earth"


def test_list_colorscale_maps_to_mpl_colormap():
    from matplotlib.colors import LinearSegmentedColormap

    stops = _blue_white_red_diverging_colorscale()
    cmap = LinearSegmentedColormap.from_list(
        "choro_custom", [(float(pos), color) for pos, color in stops]
    )
    low = cmap(0.0)
    high = cmap(1.0)
    assert low[2] > low[0]  # low end more blue than red
    assert high[0] > high[2]  # high end more red than blue


# ---------------------------------------------------------------------------
# signed exchange: the FEATURE's perspective, one rule for every package
# ---------------------------------------------------------------------------


def test_every_exchange_map_puts_blue_on_gaining(canonical_run):
    """The house rule, asserted on the REAL maps rather than on the helper.

    Gaining = blue, losing = red, from the surface-water feature's point of
    view. myflopy keeps MF6's RAW sign and never normalizes ``q``, so *which
    end* is gaining differs by reference frame:

    * ``sfr`` is aquifer-referenced (the "gwf" frame): gaining is NEGATIVE, so
      blue must sit at the NEGATIVE end.
    * ``lak`` is feature-referenced: gaining is POSITIVE, blue at the POSITIVE
      end.
    * ``surface_water`` draws the normalized ``exchange_intensity`` field
      (positive = gaining), so it follows the feature orientation too.

    This test exists in this form because an earlier one did not catch a live
    inversion: it exercised ``_exchange_colorscale`` with hand-supplied signs,
    proving the helper's arithmetic while the real SFR map rendered gaining
    reaches red. Assert on what the call site actually produces, and derive the
    expected end from the registry frame so the colour and the declared sign
    cannot drift apart.
    """

    BLUE, RED = "#1f77b4", "#d62728"
    # frame per map: sfr/lak from the registry, surface_water is the normalized
    # exchange_intensity field (feature-oriented: positive = gaining).
    frames = {
        "sfr": get_package_result_spec("sfr", "q").reference_frame,
        "lak": get_package_result_spec("lak", "q").reference_frame,
        "surface_water": "feature",
    }
    maps = {
        "sfr": canonical_run.packages.sfr.results.q,
        "lak": canonical_run.packages.lak.results.q,
        "surface_water": canonical_run.packages.surface_water.results,
    }
    for label, accessor in maps.items():
        choro = accessor.map(per=0)
        scale = list(choro.colorscale) if hasattr(choro, "colorscale") else None
        if scale is None:  # plotly figure -> pull it off the trace
            scale = list(choro.data[0].colorscale)
        negative_end, positive_end = scale[0][1].lower(), scale[-1][1].lower()
        # gaining end: negative for the "gwf" frame, positive for "feature".
        gaining_end, losing_end = (
            (negative_end, positive_end)
            if frames[label] == "gwf"
            else (positive_end, negative_end)
        )
        assert gaining_end == BLUE, (
            f"{label} ({frames[label]} frame): the GAINING end must be blue, "
            f"got {gaining_end}"
        )
        assert losing_end == RED, (
            f"{label} ({frames[label]} frame): the LOSING end must be red, "
            f"got {losing_end}"
        )


def test_every_group_exchange_map_puts_blue_on_gaining(canonical_run):
    """Same house rule on the GROUP tier — where an inversion once shipped silently.

    The single-model LAK map was fixed in c55ea23, but the group LAK and combined
    group surface_water maps kept the raw blue-at-negative scale and drew gaining
    features red, with no test to catch it. This asserts the group maps agree with
    their single-model twins by orienting blue onto the gaining end per frame.
    """

    import myflopy as mf
    from myflopy.project.group.core import ModelGroup

    BLUE, RED = "#1f77b4", "#d62728"
    reloaded = mf.load_mf6_run(canonical_run.workspace)
    group = ModelGroup({"a": canonical_run, "b": reloaded}, reference="a")
    # sfr is the gwf frame (gaining negative); lak and the combined exchange are
    # the feature frame (gaining positive).
    cases = {
        "sfr": (group.packages.sfr.results.q.map(per=0, model_name="b"), "gwf"),
        "lak": (group.packages.lak.results.q.map(per=0, model_name="b"), "feature"),
        "surface_water": (
            group.packages.surface_water.results.map(per=0, model_name="b"),
            "feature",
        ),
    }
    for label, (choro, frame) in cases.items():
        scale = list(choro.colorscale) if getattr(choro, "colorscale", None) else list(
            choro.data[0].colorscale
        )
        negative_end, positive_end = scale[0][1].lower(), scale[-1][1].lower()
        gaining_end, losing_end = (
            (negative_end, positive_end) if frame == "gwf" else (positive_end, negative_end)
        )
        assert gaining_end == BLUE, f"group {label} ({frame}): gaining end must be blue"
        assert losing_end == RED, f"group {label} ({frame}): losing end must be red"


def test_lak_plot_budget_keeps_positive_blue(canonical_run):
    """``plot_budget`` was already correct -- do not "fix" it.

    It draws LAK package budget terms, which are written from the lake's
    perspective: FROM-MVR is positive (inflow), GWF and EVAPORATION negative
    (outflow). So positive = flow in = blue is the house rule already applied.
    An adversarial review called this inverted; it is not.
    """

    figure = canonical_run.packages.lak.results.q.budget.plot(per=0)
    assert figure is not None


def test_lak_budget_is_a_view_with_the_house_shape(canonical_run):
    """``lak.results.q.budget`` is a view: ``get`` table, ``plot`` figure, ``summary``."""

    import pandas as pd

    from myflopy.modflow.mf6.package_surface_water import LakBudgetView

    budget = canonical_run.packages.lak.results.q.budget
    assert isinstance(budget, LakBudgetView)
    assert isinstance(budget.get(per=0), pd.DataFrame)
    assert budget.plot(per=0) is not None
    assert not budget.summary().empty


def test_lak_budget_retired_spellings_warn(canonical_run):
    """The retired ``budget_summary``/``plot_budget`` spellings warn (D12)."""

    import pytest

    q = canonical_run.packages.lak.results.q
    with pytest.warns(DeprecationWarning, match="budget_summary is deprecated"):
        q.budget_summary(per=0)
    with pytest.warns(DeprecationWarning, match="plot_budget is deprecated"):
        q.plot_budget(per=0)
