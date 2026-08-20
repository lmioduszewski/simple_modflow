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
``surfaces``) and **animation re-wraps** (``interactive_plotting`` rebuilds a
figure from existing data + frames). The 2-D house template (paper-anchored
border, x/y axis styling) does not belong on those, so they stay raw.

The mapbox-map exception is GONE as of plan 8.4b: the node-id and cell debug
plots that claimed it (``map_nodes``, ``plot2d``, ``plot3d``) were deleted, and
every map now runs through ``Fig`` via ``Choro`` or ``GridMesh``.

Not every picture is Plotly. :class:`MplPicture` answers the same verbs over a
Matplotlib Axes, for drawings -- filled geologic cross-sections -- that have no
Plotly equivalent to defer to.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

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
import plotly.graph_objects as go
from plotly.subplots import make_subplots as _make_subplots

__all__ = [
    "Fig",
    "Picture",
    "Subplot",
    "Template",
    "create_hover",
    "subplots",
    "mosaic",
    "shared_map_view",
    "MplPicture",
    "VtkScene",
    "FrameAnimation",
    "mpl_axes",
    "report_axes",
    "Theme",
    "REPORT",
    "PALETTE",
    "category_colors",
    "plot_cross_section",
]


class Picture:
    """The one contract every drawable myflopy object answers.

    A picture is anything you can look at: a map, a cross-section, a surface, a
    time series. Whatever kind it is, and whichever verb produced it, it answers
    the same four things::

        picture            # renders itself in Jupyter
        picture.fig        # the underlying figure, to modify before display
        picture.show()     # display it explicitly
        picture.save(path) # write it out (.html, or .png/.svg/.pdf)

    Before plan 8.1 there were four classes with four different names for their
    figure (``Choro.choropleth``, ``XSection.fig``, ``GridSection.figure``, and
    ``InterpolatedSurface`` with none at all) and three incompatible meanings for
    ``.plot()``: return the figure, show the figure, or open a browser window and
    return None. Callers learned each class separately, and ``model.cor().plot()
    .show()`` -- three calls to see one map -- was the cost.

    **Plotly subclasses supply exactly one thing: a ``fig`` property.** It must
    be idempotent -- repeated access returns the same assembled figure, never one
    that has accumulated its traces twice.

    Most pictures here are Plotly, so the four methods below are written in terms
    of ``fig``. That is a DEFAULT, not the contract: a picture whose native
    renderer is something else answers the same four verbs by overriding them --
    see :class:`MplPicture`. What callers are promised is the verbs, not the
    figure object behind them.
    """

    @property
    def fig(self) -> Fig:
        """The assembled figure. Subclasses must implement this."""

        raise NotImplementedError(
            f"{type(self).__name__} is a Picture but does not define `fig`."
        )

    def show(self, *args, **kwargs):
        """Display the picture."""

        return self.fig.show(*args, **kwargs)

    def html(self, path, *, include_plotlyjs: str = "cdn", **kwargs):
        """Write a standalone HTML file and return its path.

        ``include_plotlyjs="cdn"`` keeps the file small; pass ``True`` to inline
        plotly.js for a file that works with no network.
        """

        from pathlib import Path as _Path

        path = _Path(path)
        self.fig.write_html(str(path), include_plotlyjs=include_plotlyjs, **kwargs)
        return path

    def save(self, path, **kwargs):
        """Write the picture out, choosing the format from the suffix.

        ``.html`` goes through :meth:`html`; raster and vector formats go through
        Plotly's static export, which needs ``kaleido`` -- an optional dependency
        this package does not pin, so the error names it rather than surfacing
        Plotly's own.
        """

        from pathlib import Path as _Path

        path = _Path(path)
        if path.suffix.lower() in {".html", ".htm"}:
            return self.html(path, **kwargs)

        try:
            import kaleido  # noqa: F401  - presence check only
        except ImportError as error:
            # Raised inline rather than through `myflopy._optional.require`:
            # `viz` is deliberately externals-only (Layer 0 in the import map),
            # and importing any myflopy module here -- at module level or
            # deferred -- would either push every L0 leaf up a layer or trip the
            # deferred-import ratchet, which only moves down.
            raise ImportError(
                "kaleido is required for saving a figure as a static image "
                "(install it, or save to .html instead, which needs nothing "
                "extra)."
            ) from error
        self.fig.write_image(str(path), **kwargs)
        return path

    def _repr_mimebundle_(self, *args, **kwargs):
        """Render in Jupyter without an explicit call.

        Delegates to the figure's own mimebundle, which is how Plotly renders --
        so a picture displays exactly as its figure would, honouring the same
        renderer settings.
        """

        return self.fig._repr_mimebundle_(*args, **kwargs)


class MplPicture(Picture):
    """A :class:`Picture` whose native renderer is Matplotlib, not Plotly.

    Some drawings genuinely are Matplotlib: a filled, layer-coloured geologic
    cross-section is built by FloPy's ``PlotCrossSection``, and there is no
    Plotly equivalent to defer to. Rather than exempt those from the picture
    grammar -- or fake a ``fig`` that is not a :class:`Fig` -- this answers the
    same four verbs over an Axes.

    Subclasses implement :meth:`draw`. Everything else follows::

        picture                 # renders inline
        picture.axes            # the Matplotlib Axes, to adjust before display
        picture.show()
        picture.save("s.png")   # .png/.pdf/.svg via savefig; .html embeds a PNG

    ``fig`` deliberately RAISES here. It is documented package-wide as "the
    Plotly figure", and returning an ``mpl.Figure`` from it would break every
    caller that reasonably expects ``.add_trace``/``.update_layout``. The error
    names the alternative instead of pretending.
    """

    #: Cached Axes from the first :meth:`draw`, so the picture is idempotent the
    #: way `Picture` requires -- repeated access must not redraw.
    _axes = None

    def draw(self, ax=None, **kwargs):
        """Render into ``ax`` (or a new one) and return the Axes."""

        raise NotImplementedError(
            f"{type(self).__name__} is an MplPicture but does not define `draw`."
        )

    @property
    def axes(self):
        """The rendered Matplotlib Axes (drawn once, then cached)."""

        if self._axes is None:
            self._axes = self.draw()
        return self._axes

    @property
    def fig(self) -> Fig:
        """Not available: this picture is Matplotlib, not Plotly."""

        raise TypeError(
            f"{type(self).__name__} is drawn with Matplotlib, so it has no Plotly "
            "`fig`. Use `.axes` (or `.plot_mpl(ax=...)`) to adjust it, `.show()` "
            "to display it, and `.save(path)` to write it out."
        )

    def plot_mpl(self, ax=None, **kwargs):
        """Draw into a specific Axes -- the mosaic/panel entry point.

        Named to match :meth:`Choro.plot_mpl` and ``GridSection.plot_mpl``, which
        are the Matplotlib BACKEND of a Plotly picture. Here it is the only
        renderer, but the spelling is the same so callers need not care which.
        """

        return self.draw(ax=ax, **kwargs)

    def show(self, *args, **kwargs):
        """Display the picture."""

        import matplotlib.pyplot as plt

        self.axes  # ensure it is drawn
        return plt.show(*args, **kwargs)

    def save(self, path, *, dpi: int = 150, **kwargs):
        """Write the picture out; ``.html`` embeds the PNG in a minimal page."""

        from pathlib import Path as _Path

        path = _Path(path)
        if path.suffix.lower() in {".html", ".htm"}:
            return self.html(path, dpi=dpi, **kwargs)
        self.axes.figure.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
        return path

    def html(self, path, *, dpi: int = 150, **kwargs):
        """Write a standalone HTML page with the rendered figure inlined as a PNG.

        Self-contained and needs no network -- there is no plotly.js to fetch,
        because there is no Plotly figure. The trade against an interactive
        export is deliberate: this is the artifact you email someone.
        """

        import base64
        import io
        from pathlib import Path as _Path

        path = _Path(path)
        buffer = io.BytesIO()
        self.axes.figure.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight", **kwargs)
        encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
        title = getattr(self, "title", None) or type(self).__name__
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "<!doctype html><meta charset='utf-8'>"
            f"<title>{title}</title>"
            "<body style='margin:0;display:flex;justify-content:center'>"
            f"<img alt='{title}' style='max-width:100%' src='data:image/png;base64,{encoded}'>"
            "</body>",
            encoding="utf-8",
        )
        return path

    def _repr_mimebundle_(self, *args, **kwargs):
        """Render inline in Jupyter as a PNG."""

        import base64
        import io

        buffer = io.BytesIO()
        self.axes.figure.savefig(buffer, format="png", dpi=150, bbox_inches="tight")
        return {"image/png": base64.b64encode(buffer.getvalue()).decode("ascii")}


def _normalize_frames(frames):
    """``[picture | (label, picture), ...]`` -> ``[(label, picture), ...]``.

    The same normalization :func:`myflopy.viz.mosaic` does, so the two
    combinators accept the same shapes.
    """

    normalized = []
    for index, item in enumerate(frames):
        if isinstance(item, (tuple, list)) and len(item) == 2:
            label, picture = item
        else:
            label, picture = None, item
        if label is None:
            layout_title = getattr(getattr(picture, "layout", None), "title", None)
            label = getattr(layout_title, "text", None) or f"Frame {index + 1}"
        normalized.append((str(label), picture))
    if not normalized:
        raise ValueError("animate requires at least one frame.")
    return normalized



def _build_frame_figure(frames, *, title=None):
    """One plotly figure whose frames flip between ``[(label, picture), ...]``.

    Promoted from ``SpatialView._plotly_animation``, which had this exact shape
    but was private and reachable only from the view grammar.

    Choropleths contribute ``get_choropleth()`` -- the cell trace PLUS whatever
    is drawn over it (contours, location markers, pathlines). Reading only the
    cell trace is how an overlay silently disappears between the static picture
    and its animation, which is the bug ``viz.mosaic`` carries a comment about.
    """

    labels = [label for label, _ in frames]
    pictures = [picture for _, picture in frames]

    traces = []
    for picture in pictures:
        if hasattr(picture, "get_choropleth"):
            cells = [picture.get_choropleth()]
            overlays = list(picture.overlay_traces()) if hasattr(picture, "overlay_traces") else []
            traces.append(cells + overlays)
        else:
            traces.append(list(picture.fig.data))

    widths = {len(group) for group in traces}
    if len(widths) > 1:
        raise ValueError(
            "backend='plotly' needs every frame to have the same trace "
            f"structure; got frames with {sorted(widths)} traces. These frames "
            'are not the same picture with different data -- use backend="png".'
        )

    figure = Fig(data=traces[0])
    figure.frames = [
        go.Frame(data=group, name=name)
        for name, group in zip(labels, traces, strict=False)
    ]
    play = {"frame": {"duration": 600, "redraw": True}, "fromcurrent": True}
    pause = {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}
    layout = {
        "title": title,
        "uirevision": "lock",
        "updatemenus": [
            {
                "type": "buttons",
                "showactive": False,
                "buttons": [
                    {"label": "Play", "method": "animate", "args": [None, play]},
                    {"label": "Pause", "method": "animate", "args": [[None], pause]},
                ],
            }
        ],
        "sliders": [
            {
                "active": 0,
                "steps": [
                    {
                        "method": "animate",
                        "args": [[name], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                        "label": name,
                    }
                    for name in labels
                ],
            }
        ],
    }
    # A map flipped across frames still needs its view fitted to the data, or it
    # renders zoomed out to the world -- the same fix mosaic's subplots needed.
    if all(hasattr(picture, "get_choropleth") for picture in pictures):
        shared = shared_map_view(pictures)
        if shared:
            layout["map"] = shared
    figure.update_layout(**layout)
    return figure



_DEFAULT_PLOTLY_CONFIG = {
    "scrollZoom": True,
    "displaylogo": False,
}


def _plotly_config(fig, config: dict[str, Any] | None = None) -> dict[str, Any]:
    """Merge the standard figs Plotly config with optional export overrides."""

    merged = dict(_DEFAULT_PLOTLY_CONFIG)
    figure_config = getattr(fig, "_config", None)
    if figure_config:
        merged.update(figure_config)
    if config:
        merged.update(config)
    return merged


def _write_plotly_choropleth_restyle_html(
    fig,
    output_path: str | Path,
    *,
    include_plotlyjs: bool | str = True,
    interval_ms: int = 250,
    config: dict[str, Any] | None = None,
) -> Path:
    """Export a choropleth animation using repeatable trace-only restyles."""

    import json

    import plotly.io as pio
    from plotly.utils import PlotlyJSONEncoder

    if not fig.frames:
        raise ValueError("A choropleth restyle export requires at least one frame.")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame_payload = [
        {
            "name": str(frame.name),
            "z": frame.data[0].z,
            "customdata": frame.data[0].customdata,
        }
        for frame in fig.frames
    ]
    payload_json = json.dumps(frame_payload, cls=PlotlyJSONEncoder)

    export_fig = go.Figure(data=fig.data, layout=fig.layout)
    export_fig.frames = ()
    export_fig.update_layout(
        updatemenus=[
            {
                "type": "buttons",
                "buttons": [
                    {"label": "Play", "method": "skip", "args": []},
                    {"label": "Pause", "method": "skip", "args": []},
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "currentvalue": {"prefix": "Frame:"},
                "steps": [
                    {"label": frame["name"], "method": "skip", "args": []}
                    for frame in frame_payload
                ],
            }
        ],
    )
    post_script = f"""
const smGraph = document.getElementById("{{plot_id}}");
const smFrames = {payload_json};
let smFrameIndex = 0;
let smTimer = null;
let smPlayToken = 0;

function smStop() {{
  smPlayToken += 1;
  if (smTimer !== null) {{
    clearTimeout(smTimer);
    smTimer = null;
  }}
}}

function smApplyFrame(index) {{
  smFrameIndex = index;
  const frame = smFrames[index];
  const z = frame.z.slice();
  const customdata = frame.customdata.map(
    (row) => Array.isArray(row) ? row.slice() : row
  );
  return Plotly.restyle(
    smGraph,
    {{z: [z], customdata: [customdata]}},
    [0]
  ).then(() => Plotly.relayout(smGraph, {{"sliders[0].active": index}}));
}}

function smPlayFrame(index, token) {{
  if (token !== smPlayToken) return;
  smApplyFrame(index).then(() => {{
    if (token !== smPlayToken || index >= smFrames.length - 1) return;
    smTimer = setTimeout(() => smPlayFrame(index + 1, token), {int(interval_ms)});
  }});
}}

smGraph.on("plotly_sliderchange", (event) => {{
  smStop();
  const index = smFrames.findIndex((frame) => frame.name === String(event.step.label));
  if (index >= 0) smApplyFrame(index);
}});

smGraph.on("plotly_buttonclicked", (event) => {{
  if (event.button.label === "Pause") {{
    smStop();
    return;
  }}
  if (event.button.label !== "Play") return;
  smStop();
  const token = smPlayToken;
  smPlayFrame(0, token);
}});
smGraph.dataset.simpleModflowRestyleReady = "true";
"""
    pio.write_html(
        export_fig,
        file=output_path,
        include_plotlyjs=include_plotlyjs,
        post_script=post_script,
        config=_plotly_config(fig, config),
        auto_play=False,
        auto_open=False,
    )
    return output_path



class FrameAnimation(Picture):
    """Frames flipped in one interactive Plotly figure (plan 8.6a).

    The fast, live form: every frame is a plotly trace in a single figure with
    play/pause and a slider. Requires the frames to share a trace structure --
    they are the same picture with different data -- so it suits a field over
    stress periods, not a mixed bag.

    ``.html(path)`` does NOT write this figure verbatim. A choropleth animation
    re-embeds the whole geojson in every frame, which is why a 20-frame map of a
    2,000-cell grid weighs ~15 MB; the exported page instead ships the geometry
    once and swaps only the values, for ~1.6 MB. Same picture, one tenth the
    file -- and no colorbar flash, which the frame-based page has.
    """

    def __init__(self, frames, *, title=None):
        """Hold ``[(label, picture), ...]`` for a single flipped figure."""

        self.frames = _normalize_frames(frames)
        self.title = title
        self._fig = None

    def __repr__(self):
        """Name the frame count, which is what you check when one is missing."""

        return f"FrameAnimation({len(self.frames)} frames)"

    @property
    def labels(self) -> list[str]:
        """The per-frame slider labels, in order."""

        return [label for label, _ in self.frames]

    @property
    def fig(self):
        """The assembled play/slider figure (built once, then reused)."""

        if self._fig is None:
            self._fig = _build_frame_figure(self.frames, title=self.title)
        return self._fig

    def html(self, path, *, include_plotlyjs: bool | str = True, **kwargs):
        """Write a standalone page and return its path.

        For a CHOROPLETH animation this deliberately does not write ``self.fig``
        verbatim. A plotly frames figure re-embeds the whole cell geometry in
        EVERY frame, and its native playback repaints the colorbar each step,
        which flashes and makes a second play unreliable. The exported page
        ships the geometry once and restyles only the values.

        The saving is therefore per-frame geometry, so it grows with
        frames x cells: measured at 6 frames on a 441-cell grid it is ~17%
        (5.7 MB vs 6.9 MB, mostly inlined plotly.js either way), and at 20
        frames on a 2,000-cell grid roughly an order of magnitude. Anything
        that is not a choropleth falls through to plotly's own writer.

        ``include_plotlyjs`` defaults to ``True`` (inlined) rather than
        ``Picture``'s ``"cdn"``: this is the artifact you send someone, so it
        should open with no network.
        """

        from pathlib import Path as _Path

        figure = self.fig
        path = _Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if figure.frames and getattr(figure.frames[0].data[0], "z", None) is not None:
            return _write_plotly_choropleth_restyle_html(
                figure, path, include_plotlyjs=include_plotlyjs, **kwargs
            )
        figure.write_html(str(path), include_plotlyjs=include_plotlyjs, **kwargs)
        return path



class VtkScene(Picture):
    """A :class:`Picture` over an interactive 3-D PyVista scene (plan 8.5b).

    The third renderer, after Plotly and Matplotlib. A layered grid VOLUME and a
    bundle of pathline TUBES are neither a Plotly figure nor a Matplotlib Axes,
    but they answer the same questions -- show me, save this, write me a file I
    can send -- so they answer the same verbs::

        scene                   # renders inline
        scene.scene             # the PyVista Plotter, to adjust before display
        scene.show()
        scene.html("grid.html") # self-contained vtk.js page, ~1 MB, no network
        scene.save("grid.png")  # a screenshot

    Wraps an already-built ``Plotter``: the *building* is the caller's job and is
    where ``pyvista`` gets imported (through ``myflopy._optional.require``, which
    names the ``viz3d`` extra). This module stays externals-only, so the one
    optional-dependency error it raises itself -- for ``trame``, needed only by
    HTML export -- is written inline, exactly as :meth:`Picture.save` does for
    ``kaleido``.
    """

    def __init__(self, plotter, *, meshes=(), title: str | None = None):
        """Wrap a configured PyVista ``plotter`` (and the meshes in it)."""

        self.scene = plotter
        self.meshes = tuple(meshes)
        self.title = title

    @property
    def fig(self) -> Fig:
        """Not available: this picture is PyVista, not Plotly."""

        raise TypeError(
            f"{type(self).__name__} is a 3-D PyVista scene, so it has no Plotly "
            "`fig`. Use `.scene` (the Plotter) to adjust it, `.show()` to display "
            "it, and `.save(path)` / `.html(path)` to write it out."
        )

    def _export_html(self, filename):
        """``Plotter.export_html``, with the trame hint attached on failure."""

        try:
            return self.scene.export_html(filename)
        except ImportError as error:
            # Inline rather than via `myflopy._optional.require`: `viz` is
            # deliberately externals-only (Layer 0), and importing a myflopy
            # module here would push every L0 leaf up a layer.
            raise ImportError(
                "Writing a 3-D scene to HTML needs PyVista's trame extras. "
                "Install them with `pip install 'myflopy[viz3d]'` (or "
                "`pip install trame trame-vtk trame-vuetify`)."
            ) from error

    def html(self, path, **kwargs):
        """Write a self-contained interactive vtk.js page and return its path."""

        from pathlib import Path as _Path

        path = _Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self._export_html(str(path))
        return path

    def save(self, path, **kwargs):
        """Write the scene out; ``.html`` is interactive, image formats are stills."""

        from pathlib import Path as _Path

        path = _Path(path)
        if path.suffix.lower() in {".html", ".htm"}:
            return self.html(path, **kwargs)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.scene.screenshot(str(path), **kwargs)
        return path

    def show(self, *args, **kwargs):
        """Display the scene in its own window (or inline, per PyVista config)."""

        return self.scene.show(*args, **kwargs)

    def _repr_mimebundle_(self, *args, **kwargs):
        """Render inline in Jupyter as a self-contained vtk.js page.

        Exported to a string rather than a file, so nothing is written to disk
        just by looking at a scene -- which is what the old ``vtk_3d`` did,
        dropping an HTML file into the working directory on every call.
        """

        return {"text/html": self._export_html(None).getvalue()}


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
    colorbar=None,
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
    colorbar
        Colorbar settings for the shared map color axis -- a dict of Plotly
        ``colorbar`` properties, or a **callable** ``(cmin, cmax) -> dict``.
        Panels are pooled onto one ``coloraxis``, which discards each panel's
        own ``colorbar``; without this a log-scaled mosaic silently reads in
        log10 units. The callable form exists because the useful labels depend
        on the pooled limits, and those are only known here -- e.g.
        ``colorbar=lambda lo, hi: {"tickvals": ..., "ticktext": ...}``.
        Ignored when no panel is a map.

    Examples
    --------
    >>> viz.mosaic([
    ...     group.hds.map("F9b"),                       # a choropleth
    ...     group.packages.lak.results.stage.plot(),     # a timeseries
    ... ], ncols=2)
    """

    import numpy as np

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
        if colorbar is not None:
            # Resolved AFTER the limits: a log mosaic's real-unit ticks depend on
            # the pooled cmin/cmax, which nothing outside this function knows.
            resolved = (
                colorbar(coloraxis.get("cmin"), coloraxis.get("cmax"))
                if callable(colorbar) else colorbar
            )
            if resolved:
                coloraxis["colorbar"] = dict(resolved)
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
