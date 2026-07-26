"""myflopy plotting front door -- import every figure from here.

One backend per kind, one place to change a default or a color:

    from myflopy import viz

    fig  = viz.Fig()              # interactive Plotly (pan + scroll-zoom, house template)
    fig  = viz.subplots(2, 1)     # Plotly subplot grid, same backend
    f, ax    = viz.mpl_axes()     # static matplotlib (seaborn-whitegrid house style)
    f, axes  = viz.mpl_axes(2, 2) # matplotlib grid

Backends
--------
- **Plotly** is :class:`figs.Fig` -- the same wrapper the Choropleth maps and SFR
  profiles use: pan-to-drag, scroll-to-zoom (inline too), no logo, house template.
  Build every Plotly figure with :data:`Fig` / :func:`subplots`.
- **Matplotlib** is the shared seaborn-whitegrid helper :func:`mpl_axes` (single
  axes or a grid). For publication / report-quality matplotlib, use the figs
  report theme re-exported here (:func:`report_axes`, :data:`REPORT`, :class:`Theme`).

Per-plot-type colors live in :class:`PALETTE` so a palette change is one edit in
one findable place; pass them explicitly (e.g. ``marker_color=PALETTE.prior``).
Custom per-plot themes are fine -- keep them next to the plot, sourced from here.

Deliberate exceptions (kept on raw ``plotly.graph_objects`` by design, *not*
``Fig``): **3-D scenes** (layer/surface ``surface_trace`` plots in ``layers``,
``surfaces``), **mapbox maps** (the node-id / cell debug plots in
``grid/plotting``; the model-data choropleth already runs through ``Fig`` via the
``Choro`` class), and **animation re-wraps** (``interactive_plotting`` rebuilds a
figure from existing data + frames). The 2-D house template (paper-anchored
border, x/y axis styling) does not belong on those, so they stay raw.
"""

from __future__ import annotations

from collections import Counter

# Re-export the figs primitives the project uses, so `myflopy.viz` is a superset
# drop-in for `figs`: a module can `from myflopy import viz as f` (or
# `from myflopy.viz import Fig, create_hover`) and never import figs directly.
# External-first: the author's machine exercises the live figs project; installed
# environments (and CI) fall back to the vendored snapshot in myflopy._vendor.
try:
    from figs import Fig, Subplot, Template, create_hover
    from figs.mpl import REPORT, Theme, get_mplfig, plot_cross_section
except ImportError:  # vendored fallback for installed environments
    from myflopy._vendor.figs import Fig, Subplot, Template, create_hover
    from myflopy._vendor.figs.mpl import REPORT, Theme, get_mplfig, plot_cross_section
from plotly.subplots import make_subplots as _make_subplots

__all__ = [
    "Fig",
    "Subplot",
    "Template",
    "create_hover",
    "subplots",
    "mosaic",
    "shared_map_view",
    "mpl_axes",
    "report_axes",
    "Theme",
    "REPORT",
    "PALETTE",
    "category_colors",
    "plot_cross_section",
]


def subplots(rows: int = 1, cols: int = 1, **kwargs) -> Fig:
    """Return a Plotly subplot grid on the shared :class:`figs.Fig` backend.

    Thin wrapper over :func:`plotly.subplots.make_subplots` that yields a
    ``figs.Fig`` (so the subplot figure carries scroll-zoom + the house template).
    ``dragmode='pan'`` is set explicitly because the subplot layout copy would
    otherwise drop the template default.
    """

    fig = Fig(subplot=_make_subplots(rows=rows, cols=cols, **kwargs))
    fig.update_layout(dragmode="pan")
    return fig


def shared_map_view(panels, *, bounds=None):
    """Return one ``{style, center, zoom}`` view fitting all geo panels, or ``None``.

    Unions the WGS84 bounds of every panel that exposes ``map_view()`` (real
    ``Choro`` maps) so a set of small multiples -- or the frames of a map
    animation -- all start framed to the same site. ``bounds`` overrides the
    computed union. Returns ``None`` when no panel carries geometry.
    """

    import numpy as np

    viewable = [panel for panel in panels if hasattr(panel, "map_view")]
    if not viewable:
        return None
    if bounds is None:
        extents = [
            panel.latlon_bounds
            for panel in viewable
            if getattr(panel, "latlon_bounds", None) is not None
        ]
        if extents:
            arr = np.asarray(extents, dtype=float)
            bounds = (
                float(arr[:, 0].min()),
                float(arr[:, 1].min()),
                float(arr[:, 2].max()),
                float(arr[:, 3].max()),
            )
    view = viewable[0].map_view(bounds=bounds)
    return view or None


def _fit_map_views(fig, map_cells) -> list[str]:
    """Frame every map subplot to one shared extent; return the framed ids.

    Plotly subplots default every map cell to the world view; a composed map
    only knows where to zoom once we set its ``map``/``map2``/... subplot
    layout. All map panels are framed to the same (unioned) extent so the small
    multiples line up "at the site". ``map_cells`` is ``[(subplot_id, panel),
    ...]``; the returned subplot ids drive optional live view-syncing.
    """

    view = shared_map_view([panel for _, panel in map_cells])
    if not view:
        return []
    fig.update_layout(**{sid: view for sid, _ in map_cells})
    return [sid for sid, _ in map_cells]


def _map_sync_post_script(subplot_ids) -> str:
    """JS that keeps several map subplots panned/zoomed together after render.

    Listens for ``plotly_relayout`` on any map subplot and mirrors the new
    center/zoom onto the others, with a re-entrancy lock so the mirrored update
    does not feed back. ``{plot_id}`` is substituted with the plot div id by
    Plotly's ``post_script`` machinery.
    """

    import json

    maps = json.dumps([str(sid) for sid in subplot_ids])
    return (
        "(function(){"
        "var gd=document.getElementById('{plot_id}');"
        "if(!gd){return;}"
        "var MAPS=" + maps + ";"
        "if(MAPS.length<2){return;}"
        "var lock=false;"
        "gd.on('plotly_relayout',function(ev){"
        "if(lock){return;}"
        "var src=null,center=null,zoom=null,i,m;"
        "for(i=0;i<MAPS.length;i++){m=MAPS[i];"
        "if(ev[m+'.center']!==undefined){center=ev[m+'.center'];src=m;}"
        "if(ev[m+'.zoom']!==undefined){zoom=ev[m+'.zoom'];src=m;}"
        "}"
        "if(src===null){return;}"
        "var L=gd.layout[src]||{};"
        "if(center===null){center=L.center;}"
        "if(zoom===null){zoom=L.zoom;}"
        "var upd={},j;"
        "for(j=0;j<MAPS.length;j++){if(MAPS[j]===src){continue;}"
        "if(center!==undefined&&center!==null){upd[MAPS[j]+'.center']=center;}"
        "if(zoom!==undefined&&zoom!==null){upd[MAPS[j]+'.zoom']=zoom;}"
        "}"
        "if(Object.keys(upd).length===0){return;}"
        "lock=true;"
        "Plotly.relayout(gd,upd).then(function(){lock=false;})"
        ".catch(function(){lock=false;});"
        "});"
        "})();"
    )


def mosaic(
    panels,
    *,
    ncols: int = 3,
    title: str | None = None,
    diff: bool = False,
    sync_views: bool = True,
):
    """Compose arbitrary panel objects into one Plotly grid.

    The free-form composer of the unified view grammar: pass any mix of
    ``Choro`` choropleth maps and Plotly figures (``viz.Fig`` timeseries,
    cross-sections, ...) and get one figure back. The leaf
    ``.mosaic(by=...)`` verbs are sugar over this same composition.

    Parameters
    ----------
    panels
        A list of panels, or of ``(label, panel)`` pairs. A panel is either a
        ``Choro`` (anything exposing ``get_choropleth()``) or a Plotly figure
        whose traces are copied into its grid cell. Map panels contribute their
        cell trace **and** their overlays (contours, location markers,
        pathlines) via ``overlay_traces()``. Unlabeled panels take their figure
        title, else ``Panel <n>``.
    ncols
        Grid width; rows grow as needed.
    title
        Overall figure title.
    diff
        When ``True``, the shared map color scale is centered at zero
        (diverging), as used by the diff surfaces.
    sync_views
        Every map panel always *starts* framed to the same shared extent (the
        union of the panels' grid bounds) so the small multiples line up at the
        site. When ``sync_views`` is ``True`` (default), the panels are also
        wired to pan and zoom **together** live -- dragging or zooming one map
        moves the others (via a ``plotly_relayout`` handler injected at
        ``show()`` / ``write_html()`` time; inline notebook display shows the
        shared start view but not the live linking). Set ``False`` to let each
        map pan/zoom independently after the shared start. Panels without
        geometry (raw Plotly figures) are unaffected.

    Examples
    --------
    >>> viz.mosaic([
    ...     group.hds.map("F9b"),                       # a choropleth
    ...     group.packages.lak.results.stage.plot(),     # a timeseries
    ... ], ncols=2)
    """

    import numpy as np
    import plotly.graph_objects as go

    normalized = []
    for index, item in enumerate(panels):
        if isinstance(item, (tuple, list)) and len(item) == 2:
            label, panel = item
        else:
            label, panel = None, item
        if label is None:
            layout_title = getattr(getattr(panel, "layout", None), "title", None)
            label = getattr(layout_title, "text", None) or f"Panel {index + 1}"
        normalized.append((str(label), panel))
    if not normalized:
        raise ValueError("mosaic requires at least one panel.")

    kinds, cell_traces = [], []
    for _label, panel in normalized:
        if hasattr(panel, "get_choropleth"):
            kinds.append("map")
            # The cell trace plus everything drawn over it (contours, location
            # markers, pathlines). Copying only the cells used to silently drop
            # every overlay, so a mosaic of contoured maps lost its contours.
            overlays = panel.overlay_traces() if hasattr(panel, "overlay_traces") else []
            cell_traces.append([panel.get_choropleth(), *overlays])
        elif isinstance(panel, go.Figure):
            kinds.append("xy")
            cell_traces.append(list(panel.data))
        else:
            raise TypeError(
                f"Cannot compose a panel of type {type(panel).__name__}; pass "
                "Choro maps or Plotly figures."
            )

    ncols = min(int(ncols), len(normalized)) or 1
    nrows = -(-len(normalized) // ncols)
    specs = [[{"type": "xy"} for _ in range(ncols)] for _ in range(nrows)]
    for index, kind in enumerate(kinds):
        if kind == "map":
            specs[index // ncols][index % ncols] = {"type": "map"}
    fig = subplots(
        nrows,
        ncols,
        specs=specs,
        subplot_titles=[label for label, _ in normalized],
    )

    map_values = []
    map_colorscale = None
    map_cells = []  # (subplot_id, panel) for each map cell, in add order
    for index, (kind, traces) in enumerate(zip(kinds, cell_traces, strict=False)):
        row, col = index // ncols + 1, index % ncols + 1
        subplot_id = None
        for trace in traces:
            # Only the cell trace carries the shared color scale; overlay traces
            # (Scattermap lines/markers) have no `z` and no top-level coloraxis.
            is_field = kind == "map" and getattr(trace, "z", None) is not None
            if is_field:
                trace.coloraxis = "coloraxis"
                if map_colorscale is None:
                    map_colorscale = trace.colorscale
                values = np.asarray(trace.z, dtype=float)
                if np.isfinite(values).any():
                    map_values.append(values[np.isfinite(values)])
            fig.add_trace(trace, row=row, col=col)
            if kind == "map":
                subplot_id = getattr(fig.data[-1], "subplot", None)
        if kind == "map" and subplot_id:
            map_cells.append((subplot_id, normalized[index][1]))

    if map_cells:
        synced_ids = _fit_map_views(fig, map_cells)
        if sync_views and len(synced_ids) >= 2 and hasattr(fig, "add_post_script"):
            fig.add_post_script(_map_sync_post_script(synced_ids))

    if "map" in kinds:
        coloraxis = {"colorscale": map_colorscale or "Earth"}
        if map_values:
            finite = np.concatenate(map_values)
            if diff:
                absmax = float(np.nanmax(np.abs(finite))) or 1.0
                coloraxis.update(cmin=-absmax, cmax=absmax, cmid=0.0, cauto=False)
            else:
                coloraxis.update(
                    cmin=float(np.nanmin(finite)),
                    cmax=float(np.nanmax(finite)),
                    cauto=False,
                )
        fig.update_layout(coloraxis=coloraxis)
    fig.update_layout(title=title, uirevision="lock")
    return fig


def mpl_axes(nrows: int = 1, ncols: int = 1, *, figsize=None, **kwargs):
    """Create the shared static-matplotlib figure/axes (seaborn whitegrid).

    The everyday matplotlib backend for the project: a seaborn ``whitegrid``
    figure with no global side effects. Returns ``(fig, ax)`` for a single axes
    or ``(fig, axes)`` for a grid, matching :func:`matplotlib.pyplot.subplots`.
    For report-quality output use :func:`report_axes` instead.
    """

    import matplotlib.pyplot as plt
    import seaborn as sns

    with sns.axes_style("whitegrid"):
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, **kwargs)
    return fig, axes


def report_axes(**kwargs):
    """Publication / report-quality matplotlib ``(fig, ax)`` via the figs theme.

    Delegates to :func:`figs.mpl.get_mplfig` (house fonts, 300 dpi, styled
    gridlines). Use for figures destined for a document rather than the screen.
    """

    return get_mplfig(**kwargs)


class PALETTE:
    """Shared plot colors -- the one findable place to change the project palette.

    Plotly (rgba/rgb strings) and the matplotlib equivalents live together so a
    given semantic color (prior, posterior, measured, truth) stays consistent
    across both backends.
    """

    # Plotly
    prior = "rgba(150,150,150,0.45)"
    posterior = "rgba(31,119,180,0.55)"
    # The opaque form of `posterior`, for the summary line drawn ON TOP of a
    # translucent ensemble. Exactly `mpl_posterior` in rgb, so the two backends
    # agree; callers must not reconstruct it by string-editing `posterior`'s alpha.
    posterior_solid = "rgb(31,119,180)"
    noise = "rgba(214,39,40,0.45)"
    measured = "rgb(214,39,40)"
    truth = "rgb(214,39,40)"
    conflict = "darkorange"
    # One faint line per realization, drawn many times over. Deliberately not
    # `prior`: these are every realization of any iteration, not the prior series.
    ensemble = "rgba(80,80,80,0.35)"

    # Qualitative sequence for *named categories* -- release groups, zones,
    # scenarios: things with no order and no midpoint, where a colorscale would
    # imply a ranking that is not there. Colorblind-safe: Okabe-Ito, minus the
    # yellow that vanishes on a light basemap and with its black softened to
    # #4D4D4D. Hex is backend-neutral, so matplotlib reads the same tuple; go
    # through :func:`category_colors` rather than indexing it, so a category
    # keeps its color across figures.
    categorical = (
        "#0072B2",  # blue
        "#D55E00",  # vermillion
        "#009E73",  # green
        "#CC79A7",  # reddish purple
        "#56B4E9",  # sky blue
        "#E69F00",  # orange
        "#4D4D4D",  # dark grey
    )

    # Matplotlib equivalents
    mpl_prior = "0.6"
    mpl_posterior = "#1f77b4"
    mpl_measured = "crimson"
    mpl_truth = "crimson"
    mpl_conflict = "darkorange"
    mpl_ensemble = "0.5"
    mpl_categorical = categorical


_CATEGORY_COLORS: dict[str, str] = {}


def category_colors(names, *, memoize: bool = True) -> dict[str, str]:
    """Map category names to stable colors from :attr:`PALETTE.categorical`.

    The assignment is **memoized for the life of the process**, so a release group
    drawn blue on the pathline map is blue again on its arrival curve and its
    capture bars -- the property that makes a set of small multiples readable, and
    the reason call sites should not index the palette themselves.

    A name new to this call takes the least-used color that none of the *other*
    names in the same call already hold, so the categories of one figure stay
    distinguishable even after many unrelated names have been registered. Two
    names first seen in **separate** calls can still collide once more than
    ``len(PALETTE.categorical)`` names exist -- a 7-color palette cannot promise
    otherwise (compromise ledger 65).

    ``memoize=False`` colors this call only, leaving the shared memo untouched.
    Use it for labels that are *not* a recurring category -- per-particle ids,
    row keys -- which would otherwise fill the memo and shift the colors every
    later figure gets.
    """

    registry = _CATEGORY_COLORS if memoize else dict(_CATEGORY_COLORS)
    requested = {str(name) for name in names}
    taken = {registry[name] for name in requested & registry.keys()}
    for name in sorted(requested - registry.keys()):
        free = [color for color in PALETTE.categorical if color not in taken]
        usage = Counter(registry.values())
        color = min(
            free or PALETTE.categorical,
            key=lambda candidate: (usage[candidate], PALETTE.categorical.index(candidate)),
        )
        registry[name] = color
        taken.add(color)
    return {str(name): registry[str(name)] for name in names}
