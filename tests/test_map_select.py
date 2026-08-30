"""Highlighting cells on a map, at every scope (ledger 160).

`select=` was grid-only and drew its highlight by DIMMING everything else --
Plotly's `selectedpoints`, whose styling on a choropleth exposes only `opacity`.
On a bare node-id grid that is the right picture. On a field it is not: dimming
the unselected cells to 20% costs a measured 6.9x of readable contrast across
the cells you did not select, while the colorbar still advertises the full range.

So the default became a dissolved-boundary outline, which leaves the field at
full opacity, survives the composers, and cannot be clobbered by the user's own
box/lasso gesture. `select_style="dim"` brings the old picture back verbatim.
"""

from __future__ import annotations

import numpy as np
import pytest

from myflopy import viz


def _types(picture):
    """The trace types on an assembled figure."""

    return [trace.type for trace in picture.fig.data]


# --- reach: the request that motivated the work -------------------------------


def test_a_model_result_map_takes_select(canonical_run):
    """`model.plot.map(select=...)` -- the call that used to die inside Plotly.

    It raised `ValueError: Invalid property ... 'select'` because the name fell
    through the `**trace_kwargs` tail into the trace, and only when `.fig` was
    touched, several frames from the mistake.
    """

    cells = [0, 1, 2, 3]
    picture = canonical_run.plot.map(layer=0, select=cells)
    assert _types(picture) == ["choroplethmap", "scattermap"]


def test_the_grid_scope_still_takes_select(canonical_run):
    """`vor.plot.map(select=...)` keeps working; only its picture changed."""

    picture = canonical_run.vor.plot.map(select=[0, 1, 2])
    assert _types(picture) == ["choroplethmap", "scattermap"]


def test_the_free_verb_takes_select(canonical_run):
    """The free function is the one the bound verbs mirror."""

    from myflopy import plot

    assert _types(plot.map(canonical_run.vor, select=[0, 1])) == [
        "choroplethmap",
        "scattermap",
    ]


def test_every_scope_offers_the_same_names():
    """The narrower scope may offer fewer parameters, never other ones.

    `select` was the ONLY name the grid scope had that the free verb lacked,
    which is why `model.plot.map(select=...)` failed while `vor.plot.map` worked.
    """

    import inspect

    from myflopy import plot
    from myflopy.plot import GridPlots, ModelPlots

    def names(func):
        return {
            n
            for n, p in inspect.signature(func).parameters.items()
            if p.kind not in (p.VAR_KEYWORD, p.VAR_POSITIONAL) and n not in ("self", "source")
        }

    free = names(plot.map)
    for label, namespace in (("model", ModelPlots), ("grid", GridPlots)):
        extra = names(namespace.map) - free
        assert not extra, f"{label}.plot.map invented parameters the free verb lacks: {extra}"
    for name in ("select", "select_style", "select_color"):
        assert name in free and name in names(ModelPlots.map) and name in names(GridPlots.map)


# --- the three styles ---------------------------------------------------------


def test_outline_leaves_the_field_at_full_opacity(canonical_run):
    """The default must not touch the cells it did not select."""

    picture = canonical_run.vor.plot.map(select=[0, 1, 2])
    cells = picture.fig.data[0]
    assert cells.selectedpoints is None
    assert cells.unselected is None or cells.unselected.marker.opacity is None


def test_dim_reproduces_the_old_picture(canonical_run):
    """`select_style="dim"` is the pre-160 behaviour, stated explicitly."""

    picture = canonical_run.vor.plot.map(select=[0, 1, 2], select_style="dim")
    cells = picture.fig.data[0]
    assert tuple(cells.selectedpoints) == (0, 1, 2)
    assert cells.unselected.marker.opacity == viz.HIGHLIGHT_DIM_OPACITY
    assert _types(picture) == ["choroplethmap"]


def test_both_draws_each(canonical_run):
    picture = canonical_run.vor.plot.map(select=[0, 1, 2], select_style="both")
    assert _types(picture) == ["choroplethmap", "scattermap"]
    assert tuple(picture.fig.data[0].selectedpoints) == (0, 1, 2)


def test_an_unknown_style_says_so(canonical_run):
    with pytest.raises(ValueError, match="select_style"):
        canonical_run.vor.plot.map(select=[0], select_style="glow")


def test_the_outline_cannot_be_lassoed_away(canonical_run):
    """A lines-only scatter is not selectable, so a user gesture cannot clobber it.

    Plotly treats a scatter-like trace as selectable only when it has markers or
    text; `selectedpoints` -- which the browser's box/lasso writes -- is the same
    property `select_style="dim"` uses, so the dim IS clobberable and the outline
    is not.
    """

    outline = canonical_run.vor.plot.map(select=[0, 1, 2]).fig.data[1]
    assert outline.mode == "lines"
    assert outline.marker is None or outline.marker.size is None


def test_the_highlight_colour_is_overridable(canonical_run):
    """The palette red collides with a red-blue diverging colorscale."""

    default = canonical_run.vor.plot.map(select=[0]).fig.data[1]
    assert default.line.color == viz.PALETTE.highlight
    chosen = canonical_run.vor.plot.map(select=[0], select_color="#0072B2").fig.data[1]
    assert chosen.line.color == "#0072B2"


# --- composition: what a post-render highlight silently lost ------------------


def test_the_highlight_survives_overlay_traces(canonical_run):
    """`mosaic` and `animate` rebuild from `overlay_traces()`, never from `.fig`.

    The old implementation patched the finished figure, so a highlighted map
    placed in a mosaic came back with no highlight and no warning.
    """

    picture = canonical_run.vor.plot.map(select=[0, 1, 2])
    assert [t.type for t in picture.overlay_traces()] == ["scattermap"]


def test_the_dim_survives_the_cell_trace_rebuild(canonical_run):
    """`get_choropleth()` is the rebuild path; the dim must be set there."""

    picture = canonical_run.vor.plot.map(select=[0, 1, 2], select_style="dim")
    rebuilt = picture.get_choropleth()
    assert tuple(rebuilt.selectedpoints) == (0, 1, 2)


# --- the resolver, and the bugs it must not inherit ---------------------------


def test_a_boolean_mask_selects_the_true_cells(canonical_run):
    mask = np.zeros(canonical_run.vor.ncpl, dtype=bool)
    mask[[2, 5, 7]] = True
    picture = canonical_run.vor.plot.map(select=mask)
    assert picture._select_cells == [2, 5, 7]


def test_cells_are_sorted_and_deduplicated(canonical_run):
    picture = canonical_run.vor.plot.map(select=[5, 2, 5, 2, 9])
    assert picture._select_cells == [2, 5, 9]


def test_an_empty_selection_is_a_no_op_not_a_blanked_map(canonical_run):
    """`select=[]` used to dim every cell.

    An empty tuple serialises to `[]`, which is truthy in JavaScript, so the
    selection branch fired with nothing selected and the whole map rendered at
    20% -- an intersection that matched nothing produced a blank-looking map.
    """

    picture = canonical_run.vor.plot.map(select=[])
    assert picture._select_cells == []
    assert _types(picture) == ["choroplethmap"]
    assert picture.fig.data[0].selectedpoints is None


def test_an_out_of_range_cell_is_refused(canonical_run):
    """Out-of-range indices were accepted silently and highlighted nothing."""

    with pytest.raises(ValueError, match="outside the grid"):
        canonical_run.vor.plot.map(select=[0, 10**9])


def test_a_bad_selection_raises_at_the_call_not_at_the_figure(canonical_run):
    """Resolution happens in `__init__`, so the traceback points at the call.

    Lazy resolution surfaced the error from a notebook's display hook, cells
    away from the line that was actually wrong.
    """

    with pytest.raises(ValueError):
        canonical_run.vor.plot.map(select=[-1])


def test_a_region_name_resolves_through_the_model(canonical_run):
    """`select="all_streams"` -- the registry the package builders populate."""

    listing = canonical_run.list_regions()
    names = list(listing["name"]) if hasattr(listing, "columns") else list(listing)
    assert names, "the canonical model should register regions"
    name = next((r for r in names if canonical_run.get_region_cells(r)), None)
    if name is None:
        pytest.skip("every registered region resolved to zero cells")
    picture = canonical_run.plot.map(layer=0, select=name)
    assert picture._select_cells == sorted(set(canonical_run.get_region_cells(name)))
    # the legend is labelled with the region, not "selection"
    assert picture.fig.data[1].name == name


# --- the adjacent bug this pass also closed -----------------------------------


def test_an_imported_usg_model_is_drawable(canonical_run):
    """`plot.map(usg)` died with `AttributeError: no attribute 'gdf_vorPolys'`.

    A `UsgModel` carries its grid as `.grid` and a `BuiltModel` as
    `.context.grid`; neither answers `.vor`, so both were returned AS IF they
    were grids and failed several frames from the call. Reachable the moment
    `mf.read_usg` shipped.
    """

    from types import SimpleNamespace

    from myflopy import plot

    vor = canonical_run.vor
    assert plot._grid_of(SimpleNamespace(grid=vor)) == (vor, False)
    assert plot._grid_of(SimpleNamespace(context=SimpleNamespace(grid=vor))) == (vor, False)
    # and the dispatch stays duck-typed for anything that answers the protocol
    assert plot._grid_of(vor) == (vor, False)


def test_the_matplotlib_backend_draws_the_highlight_too(canonical_run):
    """Every grammar leaf offers `backend="mpl"`; without parity it silently drops.

    The two backends of one map object used to disagree about what a highlight
    is -- which is why `canonical_02` had to drop out of the interactive map into
    `plot_mpl(outline_regions=...)` to get cell outlines at all.
    """

    import matplotlib

    matplotlib.use("Agg")
    picture = canonical_run.vor.plot.map(select=[0, 1, 2])
    figure = picture.plot_mpl()
    # cells + the highlight boundary
    assert len(figure.axes[0].collections) == 2


def test_outline_regions_and_select_resolve_the_same_way(canonical_run):
    """`plot_mpl(outline_regions=...)` delegates to the shared resolver."""

    import matplotlib

    matplotlib.use("Agg")
    listing = canonical_run.list_regions()
    names = list(listing["name"]) if hasattr(listing, "columns") else list(listing)
    name = next((r for r in names if canonical_run.get_region_cells(r)), None)
    if name is None:
        pytest.skip("every registered region resolved to zero cells")

    from myflopy.modflow.utils.datatypes.choros import _resolve_cells

    cells, label = _resolve_cells(name, vor=canonical_run.vor, model=canonical_run)
    assert label == name
    assert cells == sorted(set(canonical_run.get_region_cells(name)))


# --- layer elevations: the model has them, so the grid need not ---------------


def test_layer_elevation_hover_does_not_need_a_grid_frame(canonical_run):
    """`show_layer_elevs=True` used to crash on any grid without `gdf_topbtm`.

    It opened by reading `self.vor.gdf_topbtm.columns` into `layer_nums`, which
    was then unconditionally overwritten and never read -- so a grid built from
    a `.gsf` died with `'NoneType' object has no attribute 'columns'` for a value
    that was thrown away. The model branch takes its elevations from the MODEL,
    which necessarily has them or it could not have run.
    """

    saved = canonical_run.vor.gdf_topbtm
    try:
        canonical_run.vor.gdf_topbtm = None
        picture = canonical_run.plot.map(layer=0, show_layer_elevs=True)
        hover = picture.hover_dict
        assert "Top of Model" in hover
        assert "Layer 1 Bottom" in hover
    finally:
        canonical_run.vor.gdf_topbtm = saved


# --- mounding: in the hover, and clear about its datum ------------------------


def _hover_lines(picture, needle):
    template = picture.fig.data[0].hovertemplate or ""
    return [line for line in template.replace("<br>", "\n").split("\n") if needle in line]


def test_mounding_reaches_the_hover(canonical_run):
    """`show_mounding=True` drew the mounding but never showed the number.

    `hover_dict` carried it, but a model map renders through `hover_spec`, which
    reads the context payload and never looks at `hover_dict` -- and a spec only
    renders fields it names.
    """

    picture = canonical_run.plot.map(layer=0, show_mounding=True)
    assert _hover_lines(picture, "Mounding")


def test_the_hover_says_what_the_mounding_is_measured_from(canonical_run):
    """"Layer 1 Mounding" did not say above WHAT, and above ground it was wrong.

    `layer=-1` measures from the model top but sets `self.layer = 0`, so the old
    label claimed layer 1's bottom as the datum when the top was used.
    """

    assert canonical_run.plot.map(layer=0, show_mounding=True)._mounding_label == (
        "Mounding above layer 1 bottom"
    )
    assert canonical_run.plot.map(layer=1, show_mounding=True)._mounding_label == (
        "Mounding above layer 2 bottom"
    )
    assert canonical_run.plot.map(layer=-1, show_mounding=True)._mounding_label == (
        "Mounding above ground surface"
    )


def test_the_mounding_datum_does_not_depend_on_what_you_ask_for_first(canonical_run):
    """`layer=-1` is resolved at construction, not on first use.

    Reading the label before the numbers used to report layer 1's bottom, because
    the above-ground flag was set as a side effect of computing the series.
    """

    picture = canonical_run.plot.map(layer=-1, show_mounding=True)
    label_first = picture._mounding_label
    picture._mounding_series()
    assert picture._mounding_label == label_first == "Mounding above ground surface"


def test_mounding_is_never_negative(canonical_run):
    """Head below the datum is no mounding, not negative mounding."""

    series = canonical_run.plot.map(layer=0, show_mounding=True)._mounding_series()
    assert (series >= 0).all()


def test_a_map_without_mounding_does_not_mention_it(canonical_run):
    assert not _hover_lines(canonical_run.plot.map(layer=0), "Mounding")


# --- layer elevations reach the sectioned hover too ---------------------------


def test_show_layer_elevs_controls_the_sectioned_hover(canonical_run):
    """`show_layer_elevs` fed only `hover_dict`; a model map renders from the spec.

    `HoverSpec._effective_blocks` decides from the SPEC's own `surfaces` flag,
    and `head_hover()` defaults it to False -- so asking for layer elevations on
    a model map produced a hover without them, and asking to turn them OFF did
    nothing either.
    """

    assert canonical_run.plot.map(show_layer_elevs=True)._resolved_hover_spec().surfaces
    assert not canonical_run.plot.map(show_layer_elevs=False)._resolved_hover_spec().surfaces


def test_hover_surfaces_stays_the_explicit_override(canonical_run):
    """`hover_surfaces=` is the precise knob and must beat the coarse one."""

    spec = canonical_run.plot.map(
        show_layer_elevs=True, hover_surfaces=False
    )._resolved_hover_spec()
    assert spec.surfaces is False


def test_the_default_map_stays_compact(canonical_run):
    """Only an EXPLICIT request opens the surfaces table.

    `show_layer_elevs=None` resolves to True merely because the grid HAS a layer
    frame -- that means "we could", not "you asked". Coupling the resolved value
    to the hover turned every head map into the full stacked table, which
    `test_head_map_active_strip_is_default_and_compact` forbids by name.
    """

    assert not canonical_run.plot.map()._resolved_hover_spec().surfaces
    assert canonical_run.vor.gdf_topbtm is not None, "and not because the frame is missing"


def test_the_elevation_rows_need_the_grids_layer_frame(canonical_run):
    """Without `gdf_topbtm` the layer table renders heads only -- no Top, no bot.

    This is what an imported USG model looked like before `read_usg` published
    the elevations: a layer table that plainly had no elevations in it.
    """

    saved = canonical_run.vor.gdf_topbtm
    try:
        canonical_run.vor.gdf_topbtm = None
        context = canonical_run.plot.map(show_layer_elevs=True)._build_hover_context()
        assert context.botm is None
    finally:
        canonical_run.vor.gdf_topbtm = saved

    context = canonical_run.plot.map(show_layer_elevs=True)._build_hover_context()
    assert context.botm is not None
    assert len(context.botm) == canonical_run.vor.nlay
