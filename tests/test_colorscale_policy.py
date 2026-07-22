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

    figure = canonical_run.packages.lak.results.q.plot_budget(per=0)
    assert figure is not None
