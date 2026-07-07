"""Unit tests for the sectioned hover spec engine (``utils.datatypes.hover``).

These construct :class:`HoverContext` directly (no model build) and assert on the
generated ``customdata`` / ``hovertemplate`` / ``hoverlabel`` -- the testable
contract behind the visual hover. Rendering fidelity in a browser is verified
separately by the canonical map tests.
"""

from __future__ import annotations

import numpy as np
import pytest

from myflopy.modflow.utils.datatypes.hover import (
    Fields,
    HoverContext,
    HoverSpec,
    HoverStyle,
    LayerTable,
    cell_input_hover,
    compare_hover,
    format_number,
    head_hover,
    lak_hover,
    result_hover,
    sfr_hover,
    surface_water_hover,
)


def _render(spec, ctx):
    customdata, template, hoverlabel = spec.render(ctx)
    return customdata, template, hoverlabel


# --- formatting -----------------------------------------------------------------
def test_format_number_precision_minus_and_units():
    assert format_number(16.108, precision=4) == "16.11"
    assert format_number(-38.1, precision=4).startswith("−")  # real minus
    assert format_number(0.0021, precision=4) == "0.0021"
    assert format_number(9.0, precision=4, unit="ft") == "9 ft"
    assert format_number(np.nan) == ""
    assert format_number(None) == ""
    assert format_number("VERTICAL") == "VERTICAL"


# --- head hover: layer modes ----------------------------------------------------
def _head_ctx(nlay=5, ncpl=2, active=2):
    heads = [[10.0 + layer + cell for cell in range(ncpl)] for layer in range(nlay)]
    botm = [[5.0 - layer for _ in range(ncpl)] for layer in range(nlay)]
    return HoverContext(
        ncpl=ncpl,
        active_layer=active,
        layer_fields={"head": heads},
        top=[21.3 for _ in range(ncpl)],
        botm=botm,
        cells=list(range(ncpl)),
        period=3,
        date="Jun 2026",
        area=[9812.0 for _ in range(ncpl)],
    )


def test_head_active_only_shows_primary_not_other_layers():
    ctx = _head_ctx()
    spec = head_hover(layers="active")
    customdata, template, hoverlabel = _render(spec, ctx)
    assert "Head" in template and "cell " in template
    assert "Layer 3" in template  # active layer annotation (0-based 2 -> "Layer 3")
    assert "L1" not in template and "L5" not in template  # no strip, no table
    # header id + primary value are customdata columns (2 columns, one row per cell)
    assert len(customdata) == ctx.ncpl
    assert "Period 3" in template and "Jun 2026" in template


def test_head_active_plus_strip_lists_other_layers_inline():
    ctx = _head_ctx()
    _, template, _ = _render(head_hover(layers="active+strip"), ctx)
    assert "L1" in template and "L2" in template and "L4" in template and "L5" in template
    assert "L3" not in template  # active layer omitted from the strip
    # strip is one muted line, not a stacked table
    assert template.count("<br>") <= 4


def test_head_all_renders_layer_table_with_active_row_accent():
    ctx = _head_ctx()
    spec = head_hover(layers="all")
    _, template, _ = _render(spec, ctx)
    for layer in range(1, 6):
        assert f"L{layer}" in template
    assert spec.style.accent in template  # active row is accent-colored


def test_head_all_with_surfaces_merges_and_marks_dry():
    # make the active layer dry in cell 0: head < that layer's bottom
    ctx = _head_ctx(active=2)
    ctx.layer_fields["head"][4][0] = -100.0  # L5 cell0 far below its bottom (5-4=1.0)
    spec = head_hover(layers="all", surfaces=True)
    customdata, template, _ = _render(spec, ctx)
    assert "head" in template and "bot" in template  # merged column headers
    assert "Top" in template  # surfaces add the Top row
    # per-cell values (incl. the dry dagger) live in customdata, not the template
    flat = "".join(str(v) for row in customdata for v in row)
    assert "†" in flat  # dry marker present somewhere


def test_head_surfaces_only_without_all_layers():
    ctx = _head_ctx()
    spec = head_hover(layers="active", surfaces=True)
    _, template, _ = _render(spec, ctx)
    # surfaces requested -> a table with bottoms appears even in active mode
    assert "bot" in template and "Top" in template


# --- package hovers -------------------------------------------------------------
def test_lak_hover_primary_and_feature_block():
    ctx = HoverContext(
        ncpl=2,
        payload={
            "q_per_area": [-0.42, -0.10],
            "q": [-38.1, -9.0],
            "flow_area": [91.0, 40.0],
            "stage": [12.31, 12.4],
            "claktype": ["VERTICAL", "HORIZONTAL"],
        },
        cells=[901, 902],
        period=3,
    )
    spec = lak_hover()
    customdata, template, _ = _render(spec, ctx)
    assert "Lake exchange" in template
    assert "q / area" in template  # relabeled primary
    assert "stage" in template and "type" in template  # feature block + relabel
    flat = "".join(str(v) for row in customdata for v in row)
    assert "ft/d" in flat  # primary unit is baked into the formatted value
    assert len(customdata) == 2


def test_sfr_hover_has_reach_fields():
    ctx = HoverContext(
        ncpl=1,
        payload={
            "q_per_length": [-1.6],
            "q": [-48.0],
            "rlen": [30.0],
            "stage": [15.21],
            "depth": [0.83],
        },
        cells=[771],
        period=3,
    )
    _, template, _ = _render(sfr_hover(), ctx)
    assert "Stream exchange" in template
    assert "stage" in template and "depth" in template and "rlen" in template


def test_cell_input_hover_uses_blue_accent_and_extra_fields():
    ctx = HoverContext(
        ncpl=1,
        payload={"elev": [14.2], "cond": [250.0]},
        cells=[512],
        period=0,
    )
    spec = cell_input_hover("elev", extra_fields=("cond",))
    _, template, _ = _render(spec, ctx)
    assert "#185FA5" in template  # input (blue) accent
    assert "cond" in template


# --- overrides + footer ---------------------------------------------------------
def test_with_fields_appends_inline_block():
    ctx = HoverContext(ncpl=1, payload={"head": [9.0], "stage": [8.5]}, cells=[1])
    spec = HoverSpec(primary="head", title="Head").with_fields("stage")
    _, template, _ = _render(spec, ctx)
    assert "stage" in template


def test_with_style_overrides_accent_and_reaches_hoverlabel():
    ctx = HoverContext(ncpl=1, payload={"head": [9.0]}, cells=[1])
    spec = HoverSpec(primary="head", title="Head").with_style(accent="#123456", font_size=13)
    _, template, hoverlabel = _render(spec, ctx)
    assert "#123456" in template
    assert hoverlabel["font"]["size"] == 13
    assert hoverlabel["bgcolor"] and hoverlabel["bordercolor"]


def test_footer_area_and_period_render():
    ctx = HoverContext(ncpl=1, payload={"head": [9.0]}, cells=[1], period=2, area=[9812.0])
    spec = HoverSpec(primary="head", title="Head", footer=("period", "area"))
    customdata, template, _ = _render(spec, ctx)
    assert "Period 2" in template
    assert "ft²" in template  # area unit in footer


def test_customdata_row_count_matches_cells_and_column_length_validated():
    ctx = HoverContext(ncpl=3, payload={"head": [1.0, 2.0, 3.0]}, cells=[0, 1, 2])
    customdata, _, _ = _render(HoverSpec(primary="head", title="Head"), ctx)
    assert len(customdata) == 3
    assert all(isinstance(row, list) for row in customdata)


def test_mismatched_column_length_raises():
    ctx = HoverContext(ncpl=3, payload={"head": [1.0, 2.0]}, cells=[0, 1, 2])
    with pytest.raises(ValueError):
        HoverSpec(primary="head", title="Head").render(ctx)


def test_hoverlabel_carries_house_style_defaults():
    hoverlabel = HoverStyle().to_hoverlabel()
    assert "Calibri" in hoverlabel["font"]["family"]
    assert hoverlabel["align"] == "left"


# --- compare / result / surface-water factories ----------------------------------
def test_compare_hover_renders_diff_primary_and_model_context():
    ctx = HoverContext(
        ncpl=2,
        payload={
            "elev": [18.4, 17.9],
            "reference_elev": [18.1, 18.0],
            "diff": [0.3, -0.1],
            "Model": ["F9b", "F9b"],
            "Reference Model": ["F9 new", "F9 new"],
        },
        cells=[0, 1],
        period=3,
    )
    spec = compare_hover("elev", "diff", title="Δ head vs reference", units={"diff": "ft"})
    customdata, template, hoverlabel = spec.render(ctx)
    assert "Δ head vs reference" in template
    assert "model" in template and "reference" in template  # value/reference labels
    assert "vs" in template  # model-context block
    assert "#534AB7" in template  # purple delta accent
    assert hoverlabel["font"]["family"].startswith("Calibri")


def test_result_hover_and_surface_water_hover_defaults():
    ctx = HoverContext(
        ncpl=1,
        payload={
            "q": [-12.4],
            "surface_water_exchange": [-1.1],
            "sfr_exchange": [-0.8],
            "lak_exchange": [-0.3],
            "source": ["sfr"],
        },
        cells=[88],
        period=3,
    )
    _, q_template, _ = result_hover("q", title="GHB q", units={"q": "ft³/d"}).render(ctx)
    assert "GHB q" in q_template

    _, sw_template, _ = surface_water_hover().render(ctx)
    assert "Surface-water exchange" in sw_template
    assert "sfr" in sw_template and "lak" in sw_template  # breakdown labels


# --- stage passthrough in LAK/SFR q payload builders ------------------------------
def test_sfr_q_payload_passes_joined_stage_through_to_hover():
    import pandas as pd
    from myflopy.modflow.mf6.package_plotting import build_sfr_q_map_payload

    frame = pd.DataFrame(
        {
            "cell": [0, 0, 1],
            "q": [-10.0, -5.0, 2.0],
            "rlen": [10.0, 5.0, 20.0],
            "reach": [1, 2, 3],
            "package": ["sfr"] * 3,
            "stage": [15.2, 15.0, 14.1],  # joined per-feature stage
        }
    )
    values, hover = build_sfr_q_map_payload(frame, ncpl=3, per=0, layer=0)
    assert "stage" in hover
    assert hover["stage"][0] == pytest.approx(15.1)  # mean of the two reaches in cell 0
    assert hover["stage"][1] == pytest.approx(14.1)
    assert np.isnan(hover["stage"][2])  # cell without stream keeps NaN (blank hover)


def test_lak_q_payload_passes_stage_and_skips_when_absent():
    import pandas as pd
    from myflopy.modflow.mf6.package_plotting import build_lak_q_map_payload

    with_stage = pd.DataFrame(
        {
            "cell": [0],
            "q": [-38.1],
            "flow_area": [91.0],
            "lake": [2],
            "stage": [12.31],
        }
    )
    _, hover = build_lak_q_map_payload(with_stage, ncpl=1, per=0, layer=0)
    assert hover["stage"][0] == pytest.approx(12.31)

    without = with_stage.drop(columns=["stage"])
    _, hover = build_lak_q_map_payload(without, ncpl=1, per=0, layer=0)
    assert "stage" not in hover  # column absent -> no hover entry, no error


def test_join_feature_stage_maps_mean_stage_per_cell():
    import pandas as pd
    from myflopy.modflow.mf6.package_surface_water import _join_feature_stage

    frame = pd.DataFrame({"cell": [5, 9], "q": [-1.0, 2.0]})
    stage_table = pd.DataFrame(
        {
            "per": [0, 0, 0, 1],
            "cell": [5, 5, 9, 5],
            "stage": [10.0, 12.0, 8.0, 99.0],  # per=1 row must be excluded
        }
    )
    joined = _join_feature_stage(frame, stage_table, per=0)
    assert joined["stage"].tolist() == pytest.approx([11.0, 8.0])
    # already-joined frames and empty stage tables are left alone
    assert _join_feature_stage(joined, stage_table, per=0) is joined
    assert _join_feature_stage(frame, stage_table.iloc[0:0], per=0) is frame


# --- call-site sugar resolution (Choro merges hover_* kwargs into the spec) ------
def _bare_choro(**attrs):
    from myflopy.modflow.utils.datatypes.choros import Choro

    choro = Choro.__new__(Choro)  # bypass heavy __init__
    choro._hover_override = None
    choro.hover_spec = None
    choro._hover_layers = None
    choro._hover_surfaces = None
    choro._hover_fields = None
    for key, value in attrs.items():
        setattr(choro, key, value)
    return choro


def test_hover_sugar_layers_and_surfaces_override_base_spec():
    choro = _bare_choro(
        hover_spec=head_hover(layers="active"),
        _hover_layers="all",
        _hover_surfaces=True,
    )
    spec = choro._resolved_hover_spec()
    assert spec.layers == "all" and spec.surfaces is True


def test_hover_full_override_wins_over_base_spec():
    choro = _bare_choro(hover_spec=head_hover(), _hover_override=lak_hover())
    assert choro._resolved_hover_spec().primary == "q_per_area"


def test_hover_fields_appends_a_block():
    base = lak_hover()
    choro = _bare_choro(hover_spec=base, _hover_fields=["stage", "depth"])
    resolved = choro._resolved_hover_spec()
    assert len(resolved.blocks) == len(base.blocks) + 1


def test_no_sugar_returns_base_spec_unchanged():
    base = head_hover(layers="active+strip")
    assert _bare_choro(hover_spec=base)._resolved_hover_spec() is base
    assert _bare_choro()._resolved_hover_spec() is None  # no spec at all
