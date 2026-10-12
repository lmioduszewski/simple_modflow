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

import re
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
    "HIGHLIGHT_DIM_OPACITY",
    "HIGHLIGHT_WIDTH",
    "PALETTE",
    "category_colors",
    "log_decade_colorbar",
    "mpl_color",
    "mpl_colormap_for",
    "plotly_colorscale_stops",
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
        """Display the picture through its figure's own ``show``.

        For a Plotly picture that is the house ``figs.Fig.show``, whose renderer
        defaults to ``"browser"``: the call opens the figure in a web browser
        tab, not inline. To see it inline in Jupyter, leave the picture as the
        last expression of a cell, or pass ``renderer="notebook"``.

        Parameters
        ----------
        *args
            Positional arguments to ``figs.Fig.show``, whose first two are
            ``renderer`` and ``config``. A third raises ``TypeError``. Use
            keywords instead.
        **kwargs
            Keywords to ``figs.Fig.show`` and then ``plotly.io.show``:
            ``renderer`` (a Plotly renderer name, default ``"browser"``),
            ``config`` (a Plotly config dict, default the house config with
            scroll-zoom on), ``validate``, and the options Plotly passes to the
            renderer (``auto_play``, ``post_script``, ``animation_opts``, ...).

        Returns
        -------
        None
        """

        return self.fig.show(*args, **kwargs)

    def html(self, path, *, include_plotlyjs: str = "cdn", **kwargs):
        """Write a standalone HTML file and return its path.

        ``include_plotlyjs="cdn"`` keeps the file small; pass ``True`` to inline
        plotly.js for a file that works with no network.

        Parameters
        ----------
        path : str or os.PathLike
            The file to write. An existing file is overwritten, and missing
            directories are created, as every picture's ``html`` does (this one
            raised ``FileNotFoundError`` until ledger 481).
        include_plotlyjs : bool or str, default "cdn"
            How the page gets plotly.js, as Plotly's ``write_html`` reads it.
            ``"cdn"`` links the CDN copy, so the file is small but needs a
            network to open. ``True`` inlines it (about 4.8 MB with Plotly
            6.9), so the file opens offline. ``False`` leaves it out, ``"directory"`` links a
            ``plotly.min.js`` next to the file, and ``"require"`` loads it
            through require.js. A string ending in ``.js`` is used as the
            script's URL.
        **kwargs
            Keywords to ``figs.Fig.write_html`` and then Plotly's
            ``write_html``: ``config`` (default the figure's own house config),
            ``auto_play``, ``include_mathjax``, ``post_script`` (default the
            figure's own, when it has one), ``full_html``, ``animation_opts``,
            ``validate``, ``default_width``, ``default_height``, ``auto_open``
            and ``div_id``. Any other keyword raises ``TypeError``.

        Returns
        -------
        pathlib.Path
            ``path``.
        """

        from pathlib import Path as _Path

        path = _Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fig.write_html(str(path), include_plotlyjs=include_plotlyjs, **kwargs)
        return path

    def save(self, path, **kwargs):
        """Write the picture out, choosing the format from the suffix.

        ``.html`` goes through :meth:`html`; raster and vector formats go through
        Plotly's static export, which needs ``kaleido`` -- an optional dependency
        this package does not pin, so the error names it rather than surfacing
        Plotly's own.

        Parameters
        ----------
        path : str or os.PathLike
            The file to write. Its suffix chooses the route, case-insensitively.
            ``.html`` or ``.htm`` goes to this picture's own :meth:`html`. Any
            other suffix goes to Plotly's ``write_image``. It reads the format
            from the suffix (``.png``, ``.jpg``/``.jpeg``, ``.webp``, ``.svg``,
            ``.pdf``) unless ``format=`` is given, refuses a format it does not
            know, and raises ``ValueError`` for a path with no suffix and no
            ``format=``. Without kaleido installed, the image route raises
            ``ImportError`` naming it before anything is written. The image
            route never creates a missing directory; the HTML route does.
        **kwargs
            For ``.html``: the keywords of this picture's :meth:`html`
            (``include_plotlyjs`` and Plotly's ``write_html`` keywords). For an
            image: Plotly's ``write_image`` keywords, which are ``format``,
            ``scale``, ``width`` and ``height`` (pixels), ``validate`` and
            ``engine``. On a map (a picture with ``export_map_view``, such as
            a ``Choro``), ``width``/``height`` also set the size the view is
            refitted for: given, else the figure's own, else Plotly's 700x500.
            The screen view is put back after the write. A view the caller
            moved on ``.fig`` is written as it stands.

        Returns
        -------
        pathlib.Path
            ``path``.
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
        figure = self.fig
        export_view = getattr(self, "export_map_view", None)
        view = None
        if export_view is not None:
            # Kaleido writes at its own default size (700x500) unless told, and a
            # view fitted for the screen cropped the top and bottom of every map
            # saved to PNG -- 78.6% of the grid's height (ledger 307). Refit for
            # the size kaleido WILL write -- an explicit width/height wins over the
            # figure's, as plotly's own `to_image` does (ledger 311) -- then put the
            # screen view back.
            import plotly.io as pio

            size = (
                kwargs.get("width") or figure.layout.width or pio.defaults.default_width,
                kwargs.get("height") or figure.layout.height or pio.defaults.default_height,
            )
            view = export_view(size)
        if not view:
            figure.write_image(str(path), **kwargs)
            return path
        # PLAIN copies: `layout.map.center` is a live object that `update_layout`
        # mutates in place, so holding it "restored" the export view (ledger 311).
        screen = {}
        for key in view:
            value = figure.layout.map[key]
            screen[key] = value.to_plotly_json() if hasattr(value, "to_plotly_json") else value
        figure.update_layout(map=view)
        try:
            figure.write_image(str(path), **kwargs)
        finally:
            figure.update_layout(map=screen)
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
        """Render into ``ax`` (or a new one) and return the Axes.

        Abstract: this base raises ``NotImplementedError``, and each subclass
        supplies the drawing. :attr:`axes` calls it once with no arguments and
        caches the result. :meth:`plot_mpl` calls it again on every call.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The Axes to draw into, such as one panel of a mosaic. The default
            ``None`` makes a new figure and Axes.
        **kwargs
            The subclass's drawing options. ``LayerSection`` merges them over
            the options it was built with and passes them to
            ``LayerBuildResult._draw_cross_section`` (``color_by``, ``cmap``,
            ``show_grid``, ``legend``, ``title``). ``LayerThicknessMap``
            accepts none, and any keyword raises ``TypeError``.

        Returns
        -------
        matplotlib.axes.Axes
            The Axes drawn into.
        """

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

        It draws AGAIN on every call and leaves the cached :attr:`axes` alone.
        It returns the Axes, where ``Choro.plot_mpl`` returns a Figure;
        ``as_mpl_figure`` turns either into a Figure.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The Axes to draw into. The default ``None`` makes a new figure.
        **kwargs
            Passed to :meth:`draw`, whose entry says what each subclass accepts.

        Returns
        -------
        matplotlib.axes.Axes
            The Axes drawn into.
        """

        return self.draw(ax=ax, **kwargs)

    def show(self, *args, **kwargs):
        """Display the picture, drawing it first if it has not been drawn.

        This is ``matplotlib.pyplot.show``, which shows EVERY open pyplot
        figure, this one included. On a non-interactive backend such as Agg
        it shows nothing, and Matplotlib issues a ``UserWarning`` that the
        figure cannot be shown.

        Parameters
        ----------
        *args
            Passed to ``matplotlib.pyplot.show``, whose backends take no
            positional argument, so any raises ``TypeError``.
        **kwargs
            Passed to ``matplotlib.pyplot.show``: ``block`` (bool or None,
            whether to wait until the windows are closed).

        Returns
        -------
        None
        """

        import matplotlib.pyplot as plt

        self.axes  # ensure it is drawn
        return plt.show(*args, **kwargs)

    def save(self, path, *, dpi: int = 150, **kwargs):
        """Write the picture out; ``.html`` embeds the PNG in a minimal page.

        Parameters
        ----------
        path : str or os.PathLike
            The file to write. Its suffix chooses the route, case-insensitively.
            ``.html`` or ``.htm`` goes to :meth:`html`, which creates a missing
            directory. Any other suffix goes to Matplotlib's ``savefig``, which
            reads the format from the suffix (``.png``, ``.pdf``, ``.svg``,
            ``.jpg``, ...) and does not create the directory. An existing file
            is overwritten.
        dpi : int, default 150
            Resolution in dots per inch, for a raster format and for the PNG
            embedded in an ``.html`` page. A vector format uses it only for any
            raster it embeds.
        **kwargs
            Keywords to ``matplotlib.figure.Figure.savefig`` (``transparent``,
            ``facecolor``, ``pad_inches``, ``metadata``, ``format``, ...). This
            method fixes ``bbox_inches="tight"``, so passing ``bbox_inches``
            raises ``TypeError``, and so does ``format`` on the ``.html`` route.

        Returns
        -------
        pathlib.Path
            ``path``.
        """

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

        Parameters
        ----------
        path : str or os.PathLike
            The page to write, whatever its suffix. Missing directories are
            created and an existing file is overwritten. The page title is the
            picture's ``title`` attribute, else its class name.
        dpi : int, default 150
            Resolution of the embedded PNG, in dots per inch.
        **kwargs
            Keywords to ``matplotlib.figure.Figure.savefig`` for the embedded
            PNG (``transparent``, ``facecolor``, ``pad_inches``, ...). This
            method fixes ``format="png"`` and ``bbox_inches="tight"``, so
            passing either raises ``TypeError``.

        Returns
        -------
        pathlib.Path
            ``path``.
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
        shared = shared_map_view(
            pictures, size=map_canvas(figure.layout, colorbar=True, controls=True)
        )
        if shared:
            layout["map"] = shared
        _pool_frame_color_axis(figure, pictures)
    figure.update_layout(**layout)
    return figure


def _pool_frame_color_axis(figure, pictures) -> None:
    """One shared colour axis for a choropleth animation (mirrors ``Choro.ani``).

    Without it every frame autoscales its own range: a colour in one frame
    means something else in the next, the bar jumps, and the restyle HTML
    export -- which restyles only ``z`` -- keeps frame 0's decade labels over a
    range that has moved, so a log frame in another decade band shows none.
    Explicit ``zmin``/``zmax`` on the first frame's cell trace (already in trace
    space) win over the pooled range. Log-scaled frames get decade labels over
    the pooled range, under any ``colorbar`` keys the first picture's caller
    gave; mixing log and linear frames raises, as ``mosaic`` does.
    """

    import numpy as np

    cells = [figure.data[0], *(frame.data[0] for frame in figure.frames)]
    if any(getattr(trace, "z", None) is None for trace in cells):
        return
    log_flags = {bool(getattr(picture, "logscale", False)) for picture in pictures}
    if len(log_flags) > 1:
        raise ValueError(
            "animate cannot pool log-scaled and linear maps onto one color axis; "
            "give every frame the same logscale=."
        )
    parts = [np.asarray(trace.z, dtype=float) for trace in cells[1:]]
    parts = [values[np.isfinite(values)] for values in parts]
    finite = np.concatenate(parts) if parts else np.asarray([])
    first = cells[1]
    cmin = first.zmin if first.zmin is not None else (float(finite.min()) if finite.size else None)
    cmax = first.zmax if first.zmax is not None else (float(finite.max()) if finite.size else None)
    colorbar = first.colorbar.to_plotly_json() if first.colorbar is not None else {}
    if log_flags == {True} and cmin is not None and cmax is not None:
        caller = getattr(pictures[0], "kwargs", {}).get("colorbar")
        caller = {} if caller is None else (
            caller if isinstance(caller, dict) else caller.to_plotly_json()
        )
        thousands = getattr(pictures[0], "number_style", "sig") == "thousands"
        colorbar = {**colorbar, **log_decade_colorbar(cmin, cmax, thousands=thousands), **caller}
    coloraxis = {"colorscale": first.colorscale, "colorbar": colorbar}
    if cmin is not None and cmax is not None:
        coloraxis.update(cmin=cmin, cmax=cmax, cauto=False)
    for trace in cells:
        trace.coloraxis = "coloraxis"
        trace.colorscale = None
        trace.zmin = None
        trace.zmax = None
    figure.update_layout(coloraxis=coloraxis)



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



#: The trace types :meth:`FrameAnimation.html`'s restyle page is written for.
_CHOROPLETH_TRACE_TYPES = ("choropleth", "choroplethmap", "choroplethmapbox")


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
        frames on a 2,000-cell grid roughly an order of magnitude.

        The restyle page is for a frame whose first trace is a CHOROPLETH;
        anything else -- 2-D lines and markers, a 3-D ``surface``, a
        ``scatter3d``, a heatmap -- goes to plotly's own writer. Until ledger
        481 the test was "the first trace carries a ``z``", which sent a
        surface or heatmap to the restyle page, whose script calls
        ``customdata.map`` on a trace with none: Play and the slider failed in
        the browser while frame 0 still drew.

        ``include_plotlyjs`` defaults to ``True`` (inlined) rather than
        ``Picture``'s ``"cdn"``: this is the artifact you send someone, so it
        should open with no network.

        Parameters
        ----------
        path : str or os.PathLike
            The page to write. Missing directories are created and an existing
            file is overwritten.
        include_plotlyjs : bool or str, default True
            How the page gets plotly.js, as Plotly's ``write_html`` reads it.
            ``True`` inlines it so the page opens offline, and ``"cdn"`` links
            it so the page is smaller but needs a network. ``False``,
            ``"directory"``, ``"require"`` and a ``.js`` URL are as in
            :meth:`Picture.html`.
        **kwargs
            On the restyle page (a choropleth animation), two keywords only,
            and any other raises ``TypeError``:

            - ``interval_ms`` (int, default 250): milliseconds between frames
              during playback.
            - ``config`` (dict): Plotly config keys merged over the figure's
              house config (scroll-zoom on, no logo).

            Otherwise, Plotly's ``write_html`` keywords via ``figs.Fig.write_html``,
            as in :meth:`Picture.html`.

        Returns
        -------
        pathlib.Path
            ``path``.
        """

        from pathlib import Path as _Path

        figure = self.fig
        path = _Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if figure.frames and getattr(figure.frames[0].data[0], "type", None) in _CHOROPLETH_TRACE_TYPES:
            return _write_plotly_choropleth_restyle_html(
                figure, path, include_plotlyjs=include_plotlyjs, **kwargs
            )
        figure.write_html(str(path), include_plotlyjs=include_plotlyjs, **kwargs)
        return path



#: The still formats PyVista's ``Plotter.screenshot`` writes (lower case).
_SCREENSHOT_SUFFIXES = frozenset({".png", ".jpeg", ".jpg", ".bmp", ".tif", ".tiff"})


class VtkScene(Picture):
    """A :class:`Picture` over an interactive 3-D PyVista scene (plan 8.5b).

    The third renderer, after Plotly and Matplotlib. A layered grid VOLUME and a
    bundle of pathline TUBES are neither a Plotly figure nor a Matplotlib Axes,
    but they answer the same questions -- show me, save this, write me a file I
    can send -- so they answer the same verbs::

        scene.show()            # IN JUPYTER, USE THIS -- see below
        scene.scene             # the PyVista Plotter, to adjust before display
        scene.html("grid.html") # self-contained vtk.js page, ~1 MB, no network
        scene.save("grid.png")  # a screenshot
        scene                   # a static inline export; read the caveat first

    **This is the one Picture where the bare expression is the wrong habit.**
    `.show()` gives PyVista's trame widget: interactive, server-backed, and
    nothing kept in the notebook file. The `_repr_mimebundle_` below instead
    emits a self-contained 1-1.5 MB HTML document, which JupyterLab commonly
    discards -- the iopub data-rate limit is 1 MB/s by default, and embedded
    vtk.js scripts can be sanitized -- and when it does render it writes those
    megabytes into the `.ipynb` for every scene. It stays because it is the right
    behaviour outside a live kernel, and because `.html(path)` shares its export
    path; it is not the display route to reach for.

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

    #: PyVista emits its vtk.js bundle as `<script type="module">`, which makes the
    #: written page unopenable from disk. See `_make_file_url_safe`.
    _MODULE_SCRIPT = '<script type="module">'

    @staticmethod
    def _make_file_url_safe(html: str) -> str:
        """Turn the exported page into one that also works from ``file://``.

        PyVista writes the vtk.js bundle as ``<script type="module">``, and a
        module script is not executed when the page is opened from disk. The
        bundle ends with ``window.OfflineLocalView = {...}``; a following CLASSIC
        script calls ``OfflineLocalView.load(...)``. So over ``file://`` the
        viewer shell runs, the loader never does, and vtk.js falls back to its
        "Drop File / Explore Scene" placeholder -- with
        ``ReferenceError: OfflineLocalView is not defined`` in the console.
        Served over HTTP the same file is fine, which is what makes this look
        like a browser bug rather than a packaging one.

        Dropping the marker is safe for this bundle and checked before use: it
        contains no ``import``/``export``/``import.meta`` and no top-level
        ``await``, so it is valid as a classic script. Ordering survives too --
        a classic script runs at parse time, i.e. BEFORE the consumer below it,
        where the module was merely deferred and happened to win a
        ``setTimeout(..., 0)`` race.
        """

        return html.replace(VtkScene._MODULE_SCRIPT, "<script>", 1)

    def html(self, path):
        """Write a self-contained interactive vtk.js page and return its path.

        The page opens from disk as well as over HTTP -- see
        :meth:`_make_file_url_safe` for why that needs help.

        Parameters
        ----------
        path : str or os.PathLike
            The page to write, whatever its suffix. Missing directories are
            created and an existing file is overwritten. Needs PyVista's trame
            extras (``ImportError`` naming them otherwise). It takes no other
            argument: PyVista's ``Plotter.export_html`` takes only the file name,
            and the keywords this accepted and ignored are refused (ledger
            481).

        Returns
        -------
        pathlib.Path
            ``path``.
        """

        from pathlib import Path as _Path

        path = _Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self._export_html(str(path))
        path.write_text(
            self._make_file_url_safe(path.read_text(encoding="utf-8")),
            encoding="utf-8",
        )
        return path

    def save(self, path, **kwargs):
        """Write the scene out; ``.html`` is interactive, image formats are stills.

        Parameters
        ----------
        path : str or os.PathLike
            The file to write. Missing directories are created on both routes.
            ``.html`` or ``.htm`` (any case) goes to :meth:`html`. Any other
            suffix is a screenshot by PyVista's ``Plotter.screenshot``:
            ``.png``, ``.jpeg``/``.jpg``, ``.bmp`` or ``.tif``/``.tiff``, any case
            (PyVista refused ``.PNG``). Another suffix, or none, raises
            ``ValueError`` -- a path with none was written with ``.png`` appended,
            so the path returned named no file (ledger 481).
        **kwargs
            For a screenshot, ``Plotter.screenshot``'s keywords:
            ``transparent_background`` (bool), ``window_size`` (``(width,
            height)`` in pixels) and ``scale`` (int, multiplies the window
            size). Passing ``return_img`` changes only what PyVista computes;
            the path is returned either way. For ``.html`` none is taken
            (``TypeError``).

        Returns
        -------
        pathlib.Path
            ``path``.
        """

        import os as _os
        from pathlib import Path as _Path

        path = _Path(path)
        suffix = path.suffix.lower()
        if suffix in {".html", ".htm"}:
            if kwargs:
                raise TypeError(
                    f"An HTML scene takes no options; {sorted(kwargs)} are screenshot "
                    "keywords (ledger 481)."
                )
            return self.html(path)
        if suffix not in _SCREENSHOT_SUFFIXES:
            raise ValueError(
                f"{path.name!r} names no image format a screenshot writes; end it in "
                f"one of {sorted(_SCREENSHOT_SUFFIXES)} (any case), or .html."
            )
        path.parent.mkdir(parents=True, exist_ok=True)
        # PyVista reads the suffix case-sensitively; write the lower-case name and
        # move it onto the one asked for.
        written = path.with_suffix(suffix)
        self.scene.screenshot(str(written), **kwargs)
        if written != path:
            _os.replace(written, path)
        return path

    def show(self, *args, **kwargs):
        """Display the scene in its own window (or inline, per PyVista config).

        Parameters
        ----------
        *args
            Positional arguments to PyVista's ``Plotter.show``, in its order:
            ``title``, ``window_size``, ``interactive``, ``auto_close``, ...
            Keywords are clearer.
        **kwargs
            Keywords to ``Plotter.show``. ``jupyter_backend`` picks the
            notebook viewer (``"trame"``, ``"html"``, ``"static"``, ...; default
            PyVista's global theme). The others include ``window_size``
            (``(width, height)`` in pixels), ``cpos`` (a camera position),
            ``screenshot`` (a path to save a still on close) and
            ``return_viewer``.

        Returns
        -------
        object
            Whatever ``Plotter.show`` returns: None, a camera position, an
            image array or a notebook viewer widget, depending on the
            arguments and the environment.
        """

        return self.scene.show(*args, **kwargs)

    def _repr_mimebundle_(self, *args, **kwargs):
        """Render inline in Jupyter as a self-contained vtk.js page.

        Exported to a string rather than a file, so nothing is written to disk
        just by looking at a scene -- which is what the old ``vtk_3d`` did,
        dropping an HTML file into the working directory on every call.
        """

        return {"text/html": self._make_file_url_safe(self._export_html(None).getvalue())}


def _is_one_extent(value) -> bool:
    """Whether ``value`` is ONE ``extent=`` by its own shape (duck-typed; see below)."""

    import numbers

    import numpy as np

    if value is None or isinstance(value, (str, bytes)):
        return True
    # A geometry, a GeoSeries/GeoDataFrame, or an area object.
    if any(hasattr(value, name) for name in ("geom_type", "total_bounds", "to_geoseries")):
        return True
    try:
        flat = np.asarray(value, dtype=object).reshape(-1).tolist()
    except (TypeError, ValueError):
        return False
    if not flat or not all(
        isinstance(item, numbers.Real) and not isinstance(item, bool) for item in flat
    ):
        return False
    shape = np.asarray(value, dtype=float).shape
    # Four numbers, or a run of (x, y) points -- NOT a run of 4-number boxes.
    return (len(shape) == 1 and shape[0] == 4) or (len(shape) == 2 and shape[1] == 2)


def is_extent_list(value) -> bool:
    """Whether ``value`` is a LIST of ``extent=`` values -- one per panel -- rather than one.

    Told apart by what the ITEMS are, never by length alone: a flat run of four
    numbers is always one box, and a run of ``(x, y)`` pairs is one set of
    points. A list whose items are each a box, a geometry, an area or a name is
    one per panel, so a 4-panel mosaic given four names frames each on its own.
    Lives here, duck-typed, because :func:`mosaic` needs it and this module
    imports nothing from myflopy (ledger 309).
    """

    if _is_one_extent(value) or not isinstance(value, (list, tuple)) or not value:
        return False
    return all(_is_one_extent(item) for item in value)


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


#: What plotly.js draws when nothing sets a size -- its own ``layout.width`` /
#: ``layout.height`` defaults. A map is fitted for THIS canvas unless the figure
#: names its own, because the real one is unknowable from Python (an inline
#: notebook figure is as wide as its cell) and fitting a SMALLER canvas than the
#: real one only adds margin, where fitting a larger one crops. The old fit
#: assumed 1000x700 and cropped a single map to 88.6% of the grid's height at
#: 1084x584, and every mosaic panel to ~10% of the grid (ledger 307).
PLOTLY_DEFAULT_SIZE = (700, 450)
#: plotly.js's default ``layout.margin`` -- what applies when neither the figure
#: nor its template names one.
PLOTLY_DEFAULT_MARGIN = {"l": 80, "r": 80, "t": 100, "b": 80}
#: Pixel height of one row of MAP panels in :func:`mosaic`.
MOSAIC_ROW_HEIGHT = 320
#: The right margin plotly GROWS to at draw time to fit a colorbar (measured 81 px
#: in plotly 6.9 on a standalone map declaring ``r=0``; 100 leaves room for wider
#: tick labels). A declared margin is not the drawn one, and fitting the declared
#: canvas cropped a wide window's edges (ledger 311).
COLORBAR_MARGIN = 100
#: The side margins plotly grows to for an animation's slider and play buttons
#: (measured 80 px each, against a declared 50/0).
CONTROLS_MARGIN = 80


def map_canvas(
    layout,
    subplot: str | None = None,
    *,
    size=None,
    colorbar: bool = False,
    controls: bool = False,
) -> tuple[float, float]:
    """``(width, height)`` in pixels of the map area plotly will draw ``layout`` into.

    ``size`` when given -- an EXPORT size, which plotly lets win over the
    figure's own (``to_image(width=)`` beats ``layout.width``) -- else the
    figure's own ``width``/``height``, else :data:`PLOTLY_DEFAULT_SIZE`; less the
    resolved margins (the figure's, else its template's, else plotly's), grown to
    what plotly adds at draw time for a ``colorbar`` and for animation
    ``controls``; times ``subplot``'s domain when one is named.
    """

    if size is not None:
        width, height = size
    else:
        width = layout.width or PLOTLY_DEFAULT_SIZE[0]
        height = layout.height or PLOTLY_DEFAULT_SIZE[1]

    template = getattr(getattr(layout, "template", None), "layout", None)
    margins = (layout.margin, getattr(template, "margin", None))

    def _margin(side):
        for margin in margins:
            value = getattr(margin, side, None) if margin is not None else None
            if value is not None:
                return float(value)
        return float(PLOTLY_DEFAULT_MARGIN[side])

    left, right = _margin("l"), _margin("r")
    if colorbar:
        right = max(right, COLORBAR_MARGIN)
    if controls:
        left, right = max(left, CONTROLS_MARGIN), max(right, CONTROLS_MARGIN)
    plot_w = float(width) - left - right
    plot_h = float(height) - _margin("t") - _margin("b")
    if subplot is not None:
        domain = layout[subplot].domain
        x0, x1 = domain.x if domain.x is not None else (0.0, 1.0)
        y0, y1 = domain.y if domain.y is not None else (0.0, 1.0)
        plot_w *= float(x1) - float(x0)
        plot_h *= float(y1) - float(y0)
    return max(plot_w, 1.0), max(plot_h, 1.0)


def shared_map_view(panels, *, bounds=None, size=None):
    """Return one ``{style, center, zoom}`` view fitting all geo panels, or ``None``.

    Unions the WGS84 bounds of every panel that exposes ``map_view()`` (real
    ``Choro`` maps) so a set of small multiples -- or the frames of a map
    animation -- all start framed to the same site. ``bounds`` overrides the
    computed union. ``size`` is the ``(width, height)`` map canvas the view is
    fitted for (see :func:`map_canvas`). Returns ``None`` when no panel carries
    geometry.
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
    view = viewable[0].map_view(bounds=bounds, size=size)
    return view or None


def _fit_map_views(fig, map_cells, boxes=None) -> tuple[list[str], bool]:
    """Frame every map subplot; return ``(framed ids, whether they DIFFER)``.

    Plotly subplots default every map cell to the world view; a composed map
    only knows where to zoom once we set its ``map``/``map2``/... subplot
    layout. ``map_cells`` is ``[(subplot_id, panel), ...]``; ``boxes`` one WGS84
    window per cell, ``None`` for "the panel's own" (its ``extent``, else its
    grid). Panels that all share one window -- the usual case, including no
    extent at all -- are framed to the UNION, so the small multiples line up "at
    the site"; panels framed to different windows (a list of extents, ledger
    309) each get their own view, and cannot be synced.
    """

    # Every map cell of one grid has the same domain size, so the first one's
    # canvas is every one's. The map_view used to be fitted for a 1000x700
    # single map, which showed ~10% of the grid in each panel (ledger 307).
    size = map_canvas(fig.layout, map_cells[0][0], colorbar=True)
    panels = [panel for _, panel in map_cells]
    boxes = list(boxes) if boxes is not None else [None] * len(panels)
    windows = [
        box if box is not None else getattr(panel, "latlon_bounds", None)
        for box, panel in zip(boxes, panels, strict=True)
    ]
    own = [
        box is not None or getattr(panel, "extent", None) is not None
        for box, panel in zip(boxes, panels, strict=True)
    ]
    distinct = any(own) and len({tuple(w) for w in windows if w is not None}) > 1
    if not distinct:
        shared = next((w for w, o in zip(windows, own, strict=True) if o), None)
        view = shared_map_view(panels, bounds=shared, size=size)
        if not view:
            return [], False
        fig.update_layout(**{sid: view for sid, _ in map_cells})
        return [sid for sid, _ in map_cells], False
    for (sid, panel), window in zip(map_cells, windows, strict=True):
        if hasattr(panel, "map_view") and window is not None:
            fig.update_layout(**{sid: panel.map_view(bounds=window, size=size)})
    return [sid for sid, _ in map_cells], True


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


_RGB_PATTERN = re.compile(r"rgba?\(([^)]*)\)")


def mpl_color(color: str):
    """A matplotlib-acceptable color for ``color``, Plotly ``rgb(...)`` included.

    Matplotlib refuses Plotly's ``rgb(214,39,40)`` spelling outright, which is
    how a highlight color that worked on one backend raised on the other. Hex and
    named colors pass straight through as an ``(r, g, b[, a])`` tuple; anything
    matplotlib cannot parse raises ``ValueError``.
    """

    from matplotlib.colors import to_rgba

    match = _RGB_PATTERN.fullmatch(str(color).strip())
    if match:
        channels = [float(c) for c in match.group(1).split(",")]
        rgb = tuple(c / 255.0 for c in channels[:3])
        alpha = channels[3] if len(channels) > 3 else 1.0
        return (*rgb, alpha)
    return to_rgba(color)


def plotly_colorscale_stops(colorscale):
    """``[(position, color), ...]`` for a NAMED Plotly colorscale, or ``None``.

    Any of Plotly's named scales, case-insensitively, with an ``_r`` suffix
    reversing it (``"viridis_r"``). ``None`` for a name Plotly does not know,
    so a caller can decide between warning and raising.
    """

    import plotly.colors as pc
    from plotly.exceptions import PlotlyError

    try:
        stops = pc.get_colorscale(str(colorscale))
    except PlotlyError:
        return None
    return [(float(pos), color) for pos, color in stops]


def mpl_colormap_for(colorscale):
    """The matplotlib colormap that renders exactly as a Plotly ``colorscale``.

    Shared so that anything drawn *over* a static choropleth -- residual
    markers, category points -- colors on the same ramp as the cells beneath
    it. Building a second colormap at the call site is how a point and the
    cell under it end up different colors for the same value.

    Both a name and explicit ``[[position, color], ...]`` stops become a
    ``LinearSegmentedColormap`` built from Plotly's OWN stops for that scale, so
    the two backends of one map agree by construction. This replaced a
    hand-written name table that mapped ``'rdbu'`` to matplotlib's REVERSED
    colormap and sent every other name to ``gist_earth`` (ledger 69/70): a
    named diverging scale rendered mirrored between backends, and
    ``colorscale="Turbo"`` was silently brown on the static picture. An
    unknown name falls back to ``'earth'``, as the ``Choro`` setter does.
    """

    from matplotlib.colors import LinearSegmentedColormap

    if isinstance(colorscale, (list, tuple)):
        stops = [(float(pos), color) for pos, color in colorscale]
        name = "choro_custom"
    else:
        stops = plotly_colorscale_stops(colorscale)
        name = f"plotly_{colorscale}"
        if stops is None:
            stops = plotly_colorscale_stops("earth")
            name = "plotly_earth"
    return LinearSegmentedColormap.from_list(
        name, [(pos, mpl_color(color)) for pos, color in stops]
    )


def log_decade_colorbar(
    low: float, high: float, *, unit: str = "", thousands: bool = False
) -> dict:
    """Colorbar ticks for a log10-scaled trace, labelled in REAL units.

    Plotly has no log color axis: a ``logscale`` map colours ``log10(value)``, so
    the bar would read ``-3 ... 2`` for a field that runs 0.001 to 100 -- the
    same numbers the hover shows, off by orders of magnitude. ``low``/``high``
    are the trace's own limits, in log10 space; the returned ``tickvals`` sit at
    whole decades inside them and the ``ticktext`` says what each decade IS.

    Fewer than three whole decades in range gives a bar with one or two ticks,
    which reads as broken, so that case falls back to five evenly spaced labels
    (still in real units). ``unit`` is appended to every label -- ``"x"`` for a
    ratio map. The one helper behind ``Choro``, ``viz.mosaic`` and the PEST field
    maps, so the three cannot label the same field three ways.
    """

    import numpy as np

    low, high = float(low), float(high)
    if not np.isfinite([low, high]).all():
        return {}
    decades = [
        decade for decade in range(int(np.floor(low)), int(np.ceil(high)) + 1)
        if low <= decade <= high
    ]
    ticks = [float(d) for d in decades] if len(decades) >= 3 else [
        float(t) for t in np.linspace(low, high, 5)
    ]
    return {
        "tickvals": ticks,
        "ticktext": [_decade_label(tick, unit, thousands=thousands) for tick in ticks],
    }


def _decade_label(tick: float, unit: str, *, thousands: bool = False) -> str:
    """``10 ** tick`` written the way a reader expects: ``1000000``, not ``1e+06``.

    A whole decade is written out in full -- ``0.00001`` and ``1000000`` alike,
    since ``:g`` reverts to exponent notation past six digits either way; a
    non-integer tick (the evenly spaced fallback) keeps three significant
    figures.
    """

    tick = float(tick)
    value = 10.0 ** tick
    if tick.is_integer():
        if tick >= 0:
            # `thousands`: 100,000 rather than 100000 (a transmissivity map, ledger 313).
            text = f"{value:,.0f}" if thousands else f"{value:.0f}"
        else:
            text = f"{value:.{int(-tick)}f}"
    elif thousands and value >= 1.0:
        text = f"{value:,.0f}"  # the evenly spaced fallback, whole: 1,673 not 1.67e+03
    else:
        text = f"{value:.3g}"
    return f"{text}{unit}"


def _log_colorbar_for_mosaic(low, high, *, thousands: bool = False) -> dict:
    """``mosaic``'s default colorbar callback for a grid of log-scaled panels.

    The pooled limits are only known inside ``mosaic``, hence the callback form;
    they are ``None`` when no panel has finite data, and then there is nothing
    to label.
    """

    if low is None or high is None:
        return {}
    return log_decade_colorbar(low, high, thousands=thousands)


def mosaic(
    panels,
    *,
    ncols: int = 3,
    title: str | None = None,
    diff: bool = False,
    sync_views: bool | None = None,
    colorbar=None,
    extent=None,
):
    """Compose arbitrary panel objects into one Plotly grid.

    The free-form composer of the unified view grammar: pass any mix of
    ``Choro`` choropleth maps and Plotly figures (``viz.Fig`` timeseries,
    cross-sections, ...) and get one figure back. The leaf
    ``.mosaic(by=...)`` verbs are sugar over this same composition.

    Parameters
    ----------
    panels : iterable of Choro, plotly Figure or (label, panel) pairs
        Required. The panels in reading order (filled row by row), each either a
        panel or a two-item ``(label, panel)`` tuple or list. A panel is either a
        ``Choro`` (anything exposing ``get_choropleth()``) or a
        ``plotly.graph_objects.Figure`` (``viz.Fig`` included) whose traces are
        copied into its grid cell; anything else raises ``TypeError``, and an
        empty iterable raises ``ValueError``. Map panels contribute their
        cell trace **and** their overlays (contours, location markers,
        pathlines) via ``overlay_traces()``. A label is drawn through ``str()``;
        an unlabeled panel (or a ``None`` label) takes its figure title, else
        ``Panel <n>``.
    ncols : int, default 3
        Grid width in panels, capped at the number of panels (``0`` means one
        column); rows grow as needed. A grid holding maps is
        :data:`MOSAIC_ROW_HEIGHT` pixels per row (plotly's default 450 px is the
        whole figure), and each map subplot is framed for its own canvas.
    title : str or None, default None
        Overall figure title; ``None`` draws none. Each panel's own title is its
        subplot label (see ``panels``), not this.
    diff : bool, default False
        When ``True``, the shared map color scale is centered at zero
        (diverging), as used by the diff surfaces: a side no panel's
        ``zmin``/``zmax`` names mirrors the pooled magnitude. Has no effect when
        no panel is a map.
    sync_views : bool or None, default None
        Map panels framed to ONE window (the union of their grids, or one shared
        ``extent``) start there together and, by default (``None``), are wired to
        pan and zoom **together** live -- dragging or zooming one map moves the
        others (via a ``plotly_relayout`` handler injected at ``show()`` /
        ``write_html()`` time; inline notebook display shows the shared start view
        but not the live linking). Panels framed to DIFFERENT windows -- a list of
        extents -- are not synced by default, since the first pan would snap every
        panel to one window; ``True`` then raises rather than doing that. Set
        ``False`` to let each map pan/zoom independently. Panels without geometry
        (raw Plotly figures) are unaffected.
    colorbar : dict, callable or None, default None
        Colorbar settings for the shared map color axis -- a dict of Plotly
        ``colorbar`` properties, or a **callable** ``(cmin, cmax) -> dict``.
        Panels are pooled onto one ``coloraxis``, which discards each panel's
        own ``colorbar``; the axis runs over the pooled data, except that a
        ``zmin``/``zmax`` a panel was given wins on its side (the widest across
        panels), so ``mosaic(zmin=100, zmax=120)`` holds the scale it names. When every map panel is log-scaled and this is not
        given, the pooled bar is relabelled in real units at whole decades
        automatically (:func:`myflopy.viz.log_decade_colorbar`); mixing
        log-scaled and linear map panels raises, since one axis cannot mean
        both. The callable form exists because other useful labels also depend
        on the pooled limits, and those are only known here -- e.g.
        ``colorbar=lambda lo, hi: {"tickvals": ..., "ticktext": ...}``.
        Ignored when no panel is a map. ``None`` (the default) leaves the bar
        as Plotly draws it, apart from the automatic log-decade labels above and
        a ``","`` thousands format when every map panel writes its numbers as
        whole thousands (transmissivity).
    extent : tuple, geometry, GeoDataFrame, str, list or None, default None
        The window every map panel opens on (see :func:`myflopy.plot.map`), or a
        LIST of them, one per panel in order -- ``None`` in the list for "this
        panel's own". A panel built with its own ``map(extent=...)`` keeps it
        when this is ``None``. Sets the VIEW only. An extent is anything
        :func:`myflopy.plot.map` takes as ``extent=`` (four numbers
        ``(xmin, ymin, xmax, ymax)``, points, a geometry, a GeoDataFrame or an
        area name). Ignored entirely when no panel is a map; otherwise a list
        must hold exactly one entry per panel (``ValueError`` otherwise), a
        non-map panel's entry must be ``None`` (``TypeError``), and a map panel
        drawn with ``fit_bounds=False`` refuses a window (``ValueError``).

    Examples
    --------
    >>> viz.mosaic([
    ...     group.hds.map("F9b"),                       # a choropleth
    ...     group.packages.lak.results.stage.plot(),     # a timeseries
    ... ], ncols=2)

    Returns
    -------
    Fig
        A single Plotly figure holding every panel as a subplot. Not a
        :class:`Picture` -- it is already the assembled figure, so use it
        directly (``.show()``, ``.write_html(...)``).

    See Also
    --------
    myflopy.plot.animate : the same panels in sequence rather than side by side.
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

    # A log-scaled panel carries log10 values, a linear one real values; pooled
    # onto ONE color axis they would share a scale that means two different
    # things, with a colorbar that can only be labelled for one of them.
    log_flags = {
        bool(getattr(panel, "logscale", False))
        for (_label, panel), kind in zip(normalized, kinds, strict=False)
        if kind == "map"
    }
    if len(log_flags) > 1:
        raise ValueError(
            "mosaic cannot pool log-scaled and linear maps onto one color axis; "
            "give every map panel the same logscale=."
        )
    # Panels that all write their numbers as whole thousands (transmissivity,
    # ledger 313) keep that on the pooled bar, which discards their own colorbars.
    thousands = {
        getattr(panel, "number_style", "sig")
        for (_label, panel), kind in zip(normalized, kinds, strict=False)
        if kind == "map"
    } == {"thousands"}
    if colorbar is None and log_flags == {True}:
        # The panels' own decade ticks are discarded with their colorbars when
        # they are pooled, so they are rebuilt here from the pooled limits.
        colorbar = (
            (lambda low, high: _log_colorbar_for_mosaic(low, high, thousands=True))
            if thousands else _log_colorbar_for_mosaic
        )
    elif colorbar is None and thousands:
        colorbar = {"tickformat": ","}

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
    if "map" in kinds:
        # Plotly's default 450 px is the WHOLE figure: a 2x2 map mosaic got ~115 px
        # per row once plotly's own margins were taken, too small to read. A
        # height per row, so framing each panel has a canvas worth fitting.
        margin = PLOTLY_DEFAULT_MARGIN
        fig.update_layout(height=MOSAIC_ROW_HEIGHT * nrows + margin["t"] + margin["b"])

    map_values = []
    map_limits = ([], [])  # caller-set trace zmin / zmax, already in trace space
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
                if getattr(trace, "zmin", None) is not None:
                    map_limits[0].append(float(trace.zmin))
                if getattr(trace, "zmax", None) is not None:
                    map_limits[1].append(float(trace.zmax))
            fig.add_trace(trace, row=row, col=col)
            if kind == "map":
                subplot_id = getattr(fig.data[-1], "subplot", None)
        if kind == "map" and subplot_id:
            map_cells.append((subplot_id, normalized[index][1]))

    if map_cells:
        # One window per MAP cell, from the composer's extent= (a list is one per
        # PANEL, map or not); None means "the panel's own".
        if extent is None:
            per_panel = [None] * len(normalized)
        elif is_extent_list(extent):
            per_panel = list(extent)
            if len(per_panel) != len(normalized):
                raise ValueError(
                    f"extent= lists {len(per_panel)} windows for {len(normalized)} "
                    "panels; give one per panel (None keeps a panel's own)."
                )
        else:
            per_panel = [extent] * len(normalized)
        for index, kind in enumerate(kinds):
            if kind != "map" and per_panel[index] is not None and is_extent_list(extent):
                raise TypeError(
                    f"panel {index + 1} is not a map, so it has no window to frame; "
                    "give it None in extent=."
                )
        map_panels = [panel for kind, (_, panel) in zip(kinds, normalized, strict=True) if kind == "map"]
        map_values_in = [value for kind, value in zip(kinds, per_panel, strict=True) if kind == "map"]
        for number, (panel, value) in enumerate(zip(map_panels, map_values_in, strict=True), 1):
            if value is not None and not getattr(panel, "fit_bounds", True):
                raise ValueError(
                    f"map panel {number} was drawn with fit_bounds=False, and extent= is "
                    "a window to FIT -- the refusal a single map makes. Draw that "
                    "panel with the default fit_bounds, or give it None in extent=."
                )
        boxes = [
            None if value is None else panel.resolve_extent(value).wgs84
            for panel, value in zip(map_panels, map_values_in, strict=True)
        ]
        synced_ids, distinct = _fit_map_views(fig, map_cells, boxes)
        if distinct and sync_views:
            raise ValueError(
                "sync_views=True with panels framed to DIFFERENT extents: the first "
                "pan would snap every panel to one window. Drop sync_views, or give "
                "every panel the same extent."
            )
        sync = (not distinct) if sync_views is None else bool(sync_views)
        if sync and len(synced_ids) >= 2 and hasattr(fig, "add_post_script"):
            fig.add_post_script(_map_sync_post_script(synced_ids))

    if "map" in kinds:
        coloraxis = {"colorscale": map_colorscale or "Earth"}
        # A panel's own zmin/zmax WIN over the pooled range, each side on its
        # own (widest across panels). They were pooled over the data regardless,
        # so `hds.mosaic(zmin=100, zmax=120)` drew 79.04..153.76 (ledger 306).
        lo = min(map_limits[0]) if map_limits[0] else None
        hi = max(map_limits[1]) if map_limits[1] else None
        if map_values or (lo is not None and hi is not None):
            finite = np.concatenate(map_values) if map_values else np.asarray([], dtype=float)
            if diff:
                # A missing side mirrors the pooled magnitude, as the static
                # composer does -- filling it with the raw data max dropped the
                # zero centre, so an all-zero diff drew in the top colour.
                absmax = (float(np.nanmax(np.abs(finite))) if finite.size else 1.0) or 1.0
                coloraxis.update(
                    cmin=lo if lo is not None else -absmax,
                    cmax=hi if hi is not None else absmax,
                    cauto=False,
                )
                if lo is None and hi is None:
                    coloraxis["cmid"] = 0.0
            else:
                coloraxis.update(
                    cmin=lo if lo is not None else float(np.nanmin(finite)),
                    cmax=hi if hi is not None else float(np.nanmax(finite)),
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


#: Line width of a `map(select=...)` highlight outline. One constant rather than a
#: `select_width` parameter: a contour SET has variable density and earns
#: `contour_width`, a highlight boundary has one job.
HIGHLIGHT_WIDTH = 2.0

#: Opacity applied to unselected cells by `select_style="dim"`. This is Plotly's
#: own implicit `DESELECTDIM`; setting it explicitly is what makes the number
#: findable and tunable in one place, at zero visual cost.
HIGHLIGHT_DIM_OPACITY = 0.2


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

    # The cell-highlight outline drawn by `map(select=...)`. Red because it is the
    # most separable colour in this palette against the default `earth` colorscale
    # (min dE76 48.6, vs 33.3 for the vermillion and 20.8 for white) -- but it
    # collides head-on with a red-blue diverging scale, which is why `select_color`
    # exists on the verb rather than this being the only answer.
    highlight = "rgb(214,39,40)"

    # Matplotlib equivalents
    mpl_prior = "0.6"
    mpl_posterior = "#1f77b4"
    mpl_measured = "crimson"
    mpl_truth = "crimson"
    mpl_conflict = "darkorange"
    mpl_ensemble = "0.5"
    mpl_highlight = "crimson"
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
