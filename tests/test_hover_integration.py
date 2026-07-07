"""End-to-end wiring of the hover spec into real choropleth maps.

Runs on the canonical model (slow) to confirm the sectioned hover reaches the
Plotly trace with real per-cell data, correct customdata width, and the house
hoverlabel -- the part the pure-unit tests can't cover.
"""

from __future__ import annotations

import pytest


@pytest.mark.slow
def test_head_map_uses_sectioned_hover_with_surfaces(canonical_run):
    model = canonical_run

    choro = model.hds.map(hover_layers="all", hover_surfaces=True)
    trace = choro.get_choropleth()

    template = trace.hovertemplate
    assert "Head" in template  # header title
    assert "L1" in template and "L2" in template  # per-layer table rows
    assert "bot" in template and "Top" in template  # surfaces merged in
    assert trace.hoverlabel.bgcolor  # house-style hoverlabel applied
    assert "Calibri" in trace.hoverlabel.font.family

    ncpl = int(model.vor.ncpl)
    assert len(trace.customdata) == ncpl  # one row per cell
    assert all(len(row) == len(trace.customdata[0]) for row in trace.customdata)


@pytest.mark.slow
def test_head_map_active_strip_is_default_and_compact(canonical_run):
    model = canonical_run
    choro = model.hds.map()  # defaults to active+strip
    trace = choro.get_choropleth()
    template = trace.hovertemplate
    assert "Head" in template
    # active+strip lists other layers on one muted line but not a full stacked table
    assert "Layer" in template  # active-layer annotation
    assert "bot" not in template  # no surfaces unless asked


@pytest.mark.slow
def test_head_map_custom_hover_spec_overrides(canonical_run):
    from myflopy.modflow.utils.datatypes.hover import HoverSpec

    model = canonical_run
    spec = HoverSpec(primary="head", title="Water level", layers="active").with_style(
        accent="#123456"
    )
    choro = model.hds.map(hover=spec)
    trace = choro.get_choropleth()
    assert "Water level" in trace.hovertemplate
    assert "#123456" in trace.hovertemplate


@pytest.mark.slow
def test_lak_q_map_uses_lake_hover_spec(canonical_run):
    model = canonical_run
    choro = model.packages.lak.results.map(field="q")
    trace = choro.get_choropleth()
    assert "Lake exchange" in trace.hovertemplate  # lak default spec title
    assert "q / area" in trace.hovertemplate  # relabeled primary
    assert trace.hoverlabel.bgcolor  # house-style hoverlabel applied
    assert len(trace.customdata) == int(model.vor.ncpl)


@pytest.mark.slow
def test_package_map_hover_spec_is_overridable(canonical_run):
    from myflopy.modflow.utils.datatypes.hover import lak_hover

    model = canonical_run
    choro = model.packages.lak.results.map(
        field="q", hover_spec=lak_hover().with_style(accent="#abcdef")
    )
    assert "#abcdef" in choro.get_choropleth().hovertemplate


@pytest.mark.slow
def test_group_lak_q_mosaic_panels_carry_sectioned_hover(canonical_run):
    # the user's original example: gp.packages.lak.results.mosaic(field='q')
    import myflopy as mf
    from myflopy.project.model_group import ModelGroup

    model = canonical_run
    reloaded = mf.load_mf6_run(model.workspace)
    group = ModelGroup({"a": model, "b": reloaded}, reference="a")

    fig = group.packages.lak.results.mosaic(field="q")
    map_traces = [t for t in fig.data if getattr(t, "type", "") == "choroplethmap"]
    assert map_traces  # one map per member
    assert all("Lake exchange" in t.hovertemplate for t in map_traces)


@pytest.mark.slow
def test_lak_map_hover_includes_joined_stage(canonical_run):
    model = canonical_run
    trace = model.packages.lak.results.map(field="q").get_choropleth()
    assert "stage" in trace.hovertemplate  # feature block field now lights up
    flat = "".join(str(v) for row in trace.customdata for v in row)
    assert "ft" in flat  # stage values formatted with units in customdata


@pytest.mark.slow
def test_heads_diff_map_uses_purple_diff_spec(canonical_run):
    import myflopy as mf
    from myflopy.project.model_group import ModelGroup

    model = canonical_run
    reloaded = mf.load_mf6_run(model.workspace)
    group = ModelGroup({"a": model, "b": reloaded}, reference="a")

    # public route: the ONE diff() verb (compare_map is internal plumbing)
    trace = group.diff().hds.map("b").get_choropleth()
    assert "Δ head" in trace.hovertemplate
    assert "#534AB7" in trace.hovertemplate  # delta accent
    assert "vs" in trace.hovertemplate  # model-context block


@pytest.mark.slow
def test_group_package_result_and_compare_maps_carry_specs(canonical_run):
    import myflopy as mf
    from myflopy.project.model_group import ModelGroup

    model = canonical_run
    reloaded = mf.load_mf6_run(model.workspace)
    group = ModelGroup({"a": model, "b": reloaded}, reference="a")

    member = group.packages.ghb.results.map("b").get_choropleth()
    assert "GHB q" in member.hovertemplate  # teal result spec

    # the ONE public diff verb: diff() surface maps carry the purple compare spec
    diff_map = group.diff().packages.ghb.results.q.map("b").get_choropleth()
    assert "Δ GHB q" in diff_map.hovertemplate


@pytest.mark.slow
def test_stage_and_input_field_maps_carry_specs(canonical_run):
    model = canonical_run
    stage = model.packages.lak.results.stage.map().get_choropleth()
    assert "LAK stage" in stage.hovertemplate

    ghb_input = model.packages.ghb.inputs.map().get_choropleth()
    assert "#185FA5" in ghb_input.hovertemplate  # blue input accent


@pytest.mark.slow
def test_mosaic_hover_sugar_threads_through_without_leaking(canonical_run):
    # hover_* sugar must reach the spec on the mosaic/package path, not the trace
    model = canonical_run

    # heads mosaic: hover_layers builds the full per-layer table in each panel
    heads = model.hds.mosaic(kind="map", by="layer", hover_layers="all")
    heads_traces = [t for t in heads.data if getattr(t, "type", "") == "choroplethmap"]
    assert heads_traces and all("L1" in t.hovertemplate for t in heads_traces)

    # package mosaic: hover_fields appends without erroring; hover_layers is a
    # harmless no-op (q has no per-layer profile)
    lak = model.packages.lak.results.mosaic(
        field="q", hover_fields=["claktype"], hover_layers="all"
    )
    lak_traces = [t for t in lak.data if getattr(t, "type", "") == "choroplethmap"]
    assert lak_traces and all("Lake exchange" in t.hovertemplate for t in lak_traces)
