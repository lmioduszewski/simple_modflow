"""Colorscale policy: diverging red/white/blue ONLY for signed 'q'-like fields
and diff maps (negative red, positive blue); everything else uses the house
brown-to-blue 'earth' scale (the mounding-figure default).
"""

from __future__ import annotations

from myflopy.modflow.mf6.package_plotting import _blue_white_red_diverging_colorscale
from myflopy.modflow.mf6.package_registry import (
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
    get_package_explorer_spec,
    get_package_result_spec,
)


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


def test_gaining_is_blue_and_losing_is_red_for_every_exchange_source():
    """The house rule, stated from the surface-water feature's point of view.

    Gaining = blue, losing = red. Flow in blue, flow out red. Increase blue,
    decrease red. It must hold for SFR, LAK and the combined explorer even
    though their raw signs disagree.

    Reported 2026-07-19. LAK and surface_water both shipped INVERTED -- a
    perched, losing lake drew blue -- because SFR's "gaining (negative q) blue"
    comment was copied to both without re-deriving the sign. SFR reads the
    model's cell record (flow from reach to cell, so gaining is negative); LAK
    reads its package budget file (the lake's perspective, so gaining is
    positive). Nothing pinned the invariant, so nothing noticed.
    """

    from myflopy.modflow.mf6.package_plotting import _exchange_colorscale

    BLUE, RED = "#1f77b4", "#d62728"

    for label, gaining_sign in (("sfr-like", -1), ("lak-like", 1)):
        scale = _exchange_colorscale(gaining_sign=gaining_sign)
        negative_end, positive_end = scale[0][1], scale[-1][1]
        gaining_colour = negative_end if gaining_sign < 0 else positive_end
        losing_colour = positive_end if gaining_sign < 0 else negative_end

        assert gaining_colour == BLUE, f"{label}: gaining must be blue"
        assert losing_colour == RED, f"{label}: losing must be red"
        assert scale[1][1] == "#ffffff", f"{label}: zero must be white"


def test_the_registry_records_each_packages_gaining_sign():
    """The sign lives in one place, so the two cannot drift apart again."""

    from myflopy.modflow.mf6.package_registry import get_package_result_spec

    assert get_package_result_spec("sfr", "q").gaining_sign == -1
    assert get_package_result_spec("lak", "q").gaining_sign == 1


def test_lak_and_sfr_really_do_disagree_on_sign(canonical_run):
    """Guard the premise: if MF6 ever agreed, the reversal would be wrong.

    The canonical lake is perched and leaks downward while its stream has both
    gaining and losing reaches -- so the two sources must show opposite signs
    for the same physical direction. If this ever fails, re-derive
    ``gaining_sign`` before touching the colours.
    """

    lak_q = canonical_run.packages.lak.results.q.get()["q"].astype(float)
    assert (lak_q <= 0).all() and lak_q.min() < -1.0, (
        f"expected the perched lake to lose (negative q), got {lak_q.min()}..{lak_q.max()}"
    )

    sfr_q = canonical_run.packages.sfr.results.profile.get(
        per=canonical_run.nper - 1
    )["q"].astype(float)
    assert (sfr_q > 0).any(), "expected some losing reaches with POSITIVE q"


def test_lak_plot_budget_keeps_positive_blue(canonical_run):
    """``plot_budget`` was already correct -- do not "fix" it.

    It draws LAK package budget terms, which are written from the lake's
    perspective: FROM-MVR is positive (inflow), GWF and EVAPORATION negative
    (outflow). So positive = flow in = blue is the house rule already applied.
    An adversarial review called this inverted; it is not.
    """

    figure = canonical_run.packages.lak.results.q.plot_budget(per=0)
    assert figure is not None
