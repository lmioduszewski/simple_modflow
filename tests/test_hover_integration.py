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


# ---------------------------------------------------------------------------
# `show_layer_elevs` on a map that carries no per-layer field (ledger 258).
#
# Reported 2026-09-14: `packages.npf.k.mosaic(logscale=True, features=[...],
# show_layer_elevs=True)` drew everything asked for EXCEPT the elevations. Two
# defects, one behind the other. `Choro._build_hover_context` fetched `top`/
# `botm` inside its HEADS branch, so every input and static-array map -- each of
# which passes `hover_heads=False` -- reached the renderer with no elevations at
# all. And `HoverSpec._effective_blocks` answered a surfaces-only request with
# `LayerTable(fields=("__surface__",))`, a key nothing has ever populated, so
# `LayerTable.build` returned at its first line. The layer frame is the GRID's;
# it has nothing to do with the mapped field.
# ---------------------------------------------------------------------------


def _elevation_rows(choro) -> list[str]:
    """The `Top` / `L<n>` rows the rendered hover template carries."""

    import re

    template = choro.get_choropleth().hovertemplate
    text = re.sub(r"<[^>]+>", "", template.replace("<br>", "\n"))
    return [line for line in text.splitlines() if re.match(r"^\s*(Top|L\d)", line.strip())]


@pytest.mark.slow
@pytest.mark.parametrize("flag", ["show_layer_elevs", "hover_surfaces"])
def test_a_static_array_map_can_show_the_grids_layer_elevations(canonical_run, flag):
    """`npf.k` carries one value per cell, no per-layer field -- and the
    elevations are still the grid's, so both spellings draw them."""

    nlay = int(canonical_run.vor.gdf_topbtm.shape[1] - 2)
    rows = _elevation_rows(canonical_run.packages.npf.k.map(**{flag: True}))
    assert len(rows) == nlay + 1, rows          # a Top row plus one per layer
    assert rows[0].strip().startswith("Top")
    assert not _elevation_rows(canonical_run.packages.npf.k.map()), "on by default"


@pytest.mark.slow
def test_a_record_noun_map_can_show_them_too(canonical_run):
    """The same holds for a package's period table reduced to one value per cell."""

    assert _elevation_rows(canonical_run.packages.drn.inputs.cond.map(show_layer_elevs=True))
    assert _elevation_rows(canonical_run.packages.drn.results.q.map(show_layer_elevs=True))


@pytest.mark.slow
def test_a_head_map_still_merges_them_into_its_per_layer_table(canonical_run):
    """The regression guard: where there IS a per-layer field, the elevations
    stay a COLUMN of that table rather than becoming a table of their own."""

    trace = canonical_run.hds.map(show_layer_elevs=True).get_choropleth()
    assert "Top" in trace.hovertemplate and "bot" in trace.hovertemplate
    # The head column and the bottom column, side by side on each table row --
    # TWO values, where a surfaces-only table has one. (The head hover also
    # carries an inline "active+strip" summary row, which is not one of these.)
    rows = _elevation_rows(canonical_run.hds.map(show_layer_elevs=True))
    paired = [r for r in rows if r.strip().startswith("L") and r.count("%{customdata[") == 2]
    assert len(paired) >= 2, rows
    assert [r for r in rows if r.strip().startswith("Top")]


@pytest.mark.slow
def test_the_elevations_are_not_fetched_unless_a_map_asks(canonical_run):
    """`gdf_topbtm` is lazy and samples rasters, so a map that does not ask must
    not pay for it. Asserted on the CONTEXT, which is where the cost is."""

    asked = canonical_run.packages.npf.k.map(show_layer_elevs=True)._build_hover_context()
    plain = canonical_run.packages.npf.k.map()._build_hover_context()
    assert asked.botm is not None and asked.top is not None
    assert plain.botm is None and plain.top is None


@pytest.mark.slow
def test_every_mosaic_panel_carries_the_elevations_and_its_own_active_layer(canonical_run):
    """The reported call. Each panel is a different layer, so each bolds its own
    row -- which is only meaningful because the table reaches every panel."""

    mosaic = canonical_run.packages.npf.k.mosaic(logscale=True, show_layer_elevs=True)
    panels = [t for t in mosaic.data if getattr(t, "z", None) is not None]
    assert len(panels) > 1
    for trace in panels:
        assert "Top" in trace.hovertemplate, "a panel lost its layer elevations"


# ---------------------------------------------------------------------------
# A static array map lists the package's OTHER arrays (ledger 259).
#
# Asked for 2026-09-14: "when I plot k or k33, I also want to show the other in
# the hover ... not bold like k, but still there for easy reference". The same
# rule the list packages already follow through `_sibling_input_fields` -- a GHB
# map answers "what head, against what conductance" -- reaching the arrays,
# where the payload carried one field and the registry names none.
# ---------------------------------------------------------------------------


def _hover_labels(choro) -> list[str]:
    """The label of each row in the rendered hover, in order."""

    import re

    template = choro.get_choropleth().hovertemplate
    return re.findall(r'<span style="color:#[0-9a-f]{6}">([^<]+?)\s*</span>', template)


def _hover_row(template: str, name: str) -> str:
    """The rendered row whose LABEL is ``name``.

    Keyed on the muted label span, not on ``>name<``: the hover's title is the
    mapped field's own name, so the loose match returns the header instead of
    the value row and any assertion about styling then reads the wrong line.
    """

    import re

    rows = [
        line for line in template.split("<br>")
        if re.search(r'<span style="color:#[0-9a-f]{6}">' + re.escape(name) + r'\s*</span>', line)
    ]
    assert rows, f"no hover row labelled {name!r}"
    return rows[0]


@pytest.mark.slow
@pytest.mark.parametrize(("field", "sibling"), [("k", "k33"), ("k33", "k")])
def test_a_conductivity_map_lists_its_counterpart(canonical_run, field, sibling):
    """Plotting either of the pair shows the other, with the mapped one primary."""

    import numpy as np

    npf = canonical_run.packages.npf
    choro = getattr(npf, field).map(layer=1)
    labels = _hover_labels(choro)
    assert field in labels and sibling in labels, labels
    assert labels.index(field) < labels.index(sibling), "the mapped field comes first"

    # the sibling's column carries the sibling's OWN values, at this layer
    trace = choro.get_choropleth()
    column = _hover_row(trace.hovertemplate, sibling)
    index = int(column.split("customdata[")[1].split("]")[0])
    shown = np.array([float(row[index]) for row in trace.customdata])
    truth = getattr(npf, sibling).get(layer=1)[sibling].to_numpy(dtype=float)
    assert np.allclose(shown, truth, rtol=1e-3)


@pytest.mark.slow
def test_the_mapped_field_is_the_bold_one_and_the_sibling_is_not(canonical_run):
    """"Not bold like k, but still there for easy reference"."""

    template = canonical_run.packages.npf.k.map().get_choropleth().hovertemplate
    primary, sibling = (_hover_row(template, name) for name in ("k", "k33"))
    assert "<b>" in primary and "font-size:15px" in primary
    assert "<b>" not in sibling and "font-size:15px" not in sibling


@pytest.mark.slow
def test_a_declared_sibling_with_no_values_is_not_listed(canonical_run):
    """FloPy answers `npf.k22` with a live object on every model and its `.array`
    is None on the many that never set it, which reached the hover as a row of
    blanks under a real field name."""

    npf = canonical_run.packages.npf
    assert "k22" in npf.k.siblings, "k22 is declared, so the skip is measured not vacuous"
    assert getattr(canonical_run.package("npf"), "k22", None) is not None
    assert canonical_run.package("npf").k22.array is None
    assert "k22" not in _hover_labels(npf.k.map())


@pytest.mark.slow
def test_storage_arrays_pair_the_same_way_and_a_lone_array_is_unchanged(canonical_run):
    """One rule, not an npf special case -- and a package with one array is untouched."""

    assert "sy" in _hover_labels(canonical_run.packages.sto.ss.map())
    assert _hover_labels(canonical_run.packages.ic.strt.map()) == ["cell", "strt", "layer"]
