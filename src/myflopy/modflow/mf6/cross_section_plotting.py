"""Matplotlib cross-section plotting for MODFLOW model results.

This module is intentionally model-aware. It focuses on report-style and
exploration figures that show discretization, layers, and simulated heads,
while optionally borrowing a matplotlib theme from :mod:`figs`.
"""

from __future__ import annotations

from collections.abc import Sequence
from contextlib import nullcontext
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import shapely as shp
from flopy.plot import PlotCrossSection
from matplotlib.colors import ListedColormap

from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg
from myflopy.viz import Fig, Picture, category_colors, mpl_axes

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from myflopy.modflow.mf6.simulation.base import SimulationBase


@dataclass(frozen=True)
class ModelCrossSectionStyle:
    """Styling controls for report-oriented MODFLOW cross sections."""

    figsize: tuple[float, float] = (20, 12)
    xlabel: str = "Distance (ft)"
    ylabel: str = "Elevation (ft)"
    title: str | None = None
    title_fontsize: float = 24
    label_fontsize: float = 16
    tick_labelsize: float = 12
    grid_linewidth: float = 0.2
    grid_color: str = "k"
    layer_alpha: float = 0.85
    head_color: str = "blue"
    head_linewidth: float = 1.5
    head_label: str = "Simulated Head"
    legend_loc: str = "upper right"
    #: None means "pick from the background". A fixed colour is inherited from
    #: `text.color`, which a notebook theme sets -- white text on the white
    #: figure `mpl_axes` builds is a legend nobody can read.
    legend_labelcolor: str | None = None
    legend_frameon: bool = False
    legend_fontsize: float = 16
    use_figs_theme: bool = False
    layer_colors: tuple[str, ...] = field(default_factory=tuple)


def _theme_context(style: ModelCrossSectionStyle):
    """A matplotlib theme context for the figs house style (a no-op if disabled/unavailable)."""

    if not style.use_figs_theme:
        return nullcontext()
    try:
        from myflopy.viz import REPORT
    except ImportError:
        # The house theme is optional styling; a model cross-section must still
        # draw without it.
        return nullcontext()
    return REPORT.context()


def normalize_line(line) -> dict:
    """Normalize a section line to a ``{"line": coords}`` dict.

    Accepts a dict, a ``LineString``, a sequence of ``(x, y)`` pairs, or a path
    to a vector file to read one from -- the same set every ``section(line=...)``
    entry point takes, so the filled branch and the profile branch cannot
    disagree about what a line is.
    """

    if isinstance(line, dict):
        return line

    if isinstance(line, (str, Path)):
        geometry = read_shp_gpkg(Path(line)).union_all()
        return normalize_line(geometry)

    if isinstance(line, shp.MultiLineString):
        merged = shp.line_merge(line)
        return normalize_line(
            merged if isinstance(merged, shp.LineString) else max(merged.geoms, key=lambda g: g.length)
        )

    if isinstance(line, shp.LineString):
        return {"line": list(line.coords)}

    if isinstance(line, Sequence) and not isinstance(line, (str, bytes)):
        return {"line": list(line)}

    raise TypeError(f"Unsupported line specification: {type(line)}")


def _resolve_layer_array(grid) -> np.ndarray:
    """A ``(nlay, ncpl)`` layer-index array for coloring the section by layer."""

    nlay = grid.nlay
    if hasattr(grid, "ncpl"):
        ncpl = grid.ncpl
        return np.arange(nlay)[:, None] * np.ones((1, ncpl), dtype=int)
    raise ValueError("layered cross sections currently expect a grid with an 'ncpl' attribute.")


def _masked_to_layers(modelgrid, values, chosen):
    """``(layer_id, values)`` as a filled section paints them: NaN outside ``chosen``.

    Masked, not dropped: FloPy's ``plot_array`` skips NaN cells, which keeps every
    other cell's geometry exactly where it was. Both renderers of a filled section
    paint from these arrays, so they cannot disagree about which cells are drawn.
    """

    nlay = modelgrid.nlay
    layer_id = _resolve_layer_array(modelgrid)
    if chosen is not None:
        keep = np.zeros(nlay, dtype=bool)
        keep[chosen] = True
        layer_id = np.where(keep[:, None], layer_id.astype(float), np.nan)
        if values is not None:
            values = np.where(keep[:, None], np.asarray(values, dtype=float), np.nan)
    return layer_id, values


def _resolve_head_data(
    model: SimulationBase,
    *,
    kstpkper: tuple | None,
    head_data,
):
    """Resolve the head array to overlay: explicit ``head_data``, else read at ``kstpkper``, else ``None``."""

    if head_data is not None:
        return np.asarray(head_data)
    if kstpkper is None:
        return None
    return np.asarray(model.gwf.output.head().get_data(kstpkper=kstpkper))


def _resolve_head_surface(head_data, head_layer: int):
    """The single head surface for the water-table line: the array itself, or its ``head_layer`` slice."""

    if head_data is None:
        return None
    arr = np.asarray(head_data)
    if arr.ndim == 1:
        return arr
    if arr.ndim >= 2:
        return arr[head_layer]
    return arr


def _resolve_layer_colors(nlay: int, layer_colors: Sequence[str] | None, style: ModelCrossSectionStyle) -> list[str]:
    """Exactly ``nlay`` layer fill colors from the explicit list, the style, or defaults (last repeated)."""

    if layer_colors is not None:
        colors = list(layer_colors)
    elif style.layer_colors:
        colors = list(style.layer_colors)
    else:
        defaults = ["#fff6cc", "#d4ac6e", "#c7d7b5", "#b9cbe1", "#d6c0e8"]
        colors = defaults[:]
    if len(colors) < nlay:
        colors.extend([colors[-1]] * (nlay - len(colors)))
    return colors[:nlay]


def _resolve_layer_labels(nlay: int, layer_labels: Sequence[str] | None) -> list[str]:
    """Exactly ``nlay`` layer labels from the explicit list, padded with ``Layer N`` defaults."""

    if layer_labels is None:
        return [f"Layer {i + 1}" for i in range(nlay)]
    labels = list(layer_labels)
    if len(labels) < nlay:
        labels.extend([f"Layer {i + 1}" for i in range(len(labels), nlay)])
    return labels[:nlay]


def _build_legend_handles(
    *,
    layer_colors: Sequence[str],
    layer_labels: Sequence[str],
    surfaces: Sequence[tuple[str, str]] = (),
    head_linewidth: float = 1.5,
) -> list:
    """Legend handles: a colour patch per drawn layer, then one line per surface."""

    handles = [
        mpatches.Patch(color=color, label=label)
        for color, label in zip(layer_colors, layer_labels, strict=False)
    ]
    handles += [
        mlines.Line2D([], [], color=colour, linewidth=head_linewidth, label=label)
        for label, colour in surfaces
    ]
    return handles


def readable_on(background) -> str:
    """``"black"`` or ``"white"``, whichever can be read on ``background``.

    Reported as "the legend text is unreadable because it is white on a white
    background". Matplotlib takes label colour from the `text.color` rcParam,
    which a dark notebook theme sets globally -- while `mpl_axes` builds a
    seaborn *whitegrid* figure, so the text went white on white. Inheriting an
    ambient setting is the same class of defect as `PALETTE.highlight` being a
    Plotly colour matplotlib rejects: state from elsewhere decides the picture.

    Derived rather than forced to black, so a genuinely dark figure does not get
    black-on-black in exchange.
    """

    from matplotlib.colors import to_rgb

    try:
        red, green, blue = to_rgb(background)
    except (ValueError, TypeError):                        # pragma: no cover
        return "black"
    # Rec. 601 luma -- the standard weighting for perceived brightness.
    return "black" if (0.299 * red + 0.587 * green + 0.114 * blue) > 0.5 else "white"


def _resolve_drawn_layers(layers, nlay: int) -> list[int] | None:
    """The zero-based layers to draw, or None for all of them.

    Validated against ``nlay`` here rather than left to produce an empty picture:
    a section that silently draws nothing is the same defect as one that draws
    the wrong thing.
    """

    if layers is None:
        return None
    chosen = [layers] if isinstance(layers, (int, np.integer)) else list(layers)
    chosen = [int(k) for k in chosen]
    bad = [k for k in chosen if not 0 <= k < nlay]
    if bad:
        raise ValueError(
            f"layers={bad} are outside this model's {nlay} layers (0..{nlay - 1})."
        )
    if not chosen:
        raise ValueError("layers= selected nothing; omit it to draw them all.")
    return sorted(set(chosen))


def _resolve_head_surfaces(head_surfaces, head_surface, style) -> list:
    """Normalize the surface arguments to ``[(label, array), ...]``."""

    if head_surfaces:
        return [(str(label), np.asarray(a, dtype=float)) for label, a in head_surfaces]
    if head_surface is not None:
        return [(style.head_label, np.asarray(head_surface, dtype=float))]
    return []


def _surface_colours(labels: Sequence[str], style) -> list[str]:
    """One colour per drawn surface, keyed by LABEL.

    A single surface keeps ``style.head_color`` exactly, so every figure drawn
    before several were possible is unchanged. Several go through
    :func:`~myflopy.viz.category_colors`, which memoizes by name -- so a given
    layer's water level is the same colour on every section it appears in, and
    the palette is already hex (matplotlib reads it; Plotly's ``rgb(...)`` form
    would not).
    """

    labels = list(labels)
    if len(labels) <= 1:
        return [style.head_color]
    mapping = category_colors(labels)
    return [mapping[label] for label in labels]


def _drawn_elevation_bounds(axes):
    """``(ymin, ymax)`` of what is actually ON the axes, padded, or None.

    Measured from the drawn polygons rather than from the grid's own
    ``top``/``botm``, which was the first cut and was wrong in a way that looked
    almost right: those are whole-GRID statistics, so the limits came from cells
    the section line never crosses. Measured on the canonical model,
    ``layers=[0]`` gave an axis of 66.3..164.2 for content spanning 72.0..146.2 --
    a quarter of the height empty, which reads as "the axis still thinks the
    other layers are there".
    """

    values = [
        path.vertices[:, 1]
        for collection in axes.collections
        for path in collection.get_paths()
        if len(path.vertices)
    ]
    # The water surfaces are `Line2D`, not collections. Leaving them out cropped
    # them off the picture -- `layers=[3]` with the default `head_layers=0` drew
    # a legend entry for a line above the top of the axis.
    values += [
        np.asarray(line.get_ydata(), dtype=float)
        for line in axes.lines
        if len(line.get_ydata())
    ]
    return _padded_elevation_bounds(values)


def _padded_elevation_bounds(values):
    """``(ymin, ymax)`` over ``values`` -- arrays of drawn elevations -- padded 5%.

    The arithmetic behind :func:`_drawn_elevation_bounds`, separate so the Plotly
    renderer frames a cropped section from the same numbers rather than from a
    second rule of its own.
    """

    if not values:                                         # pragma: no cover
        return None
    stacked = np.concatenate(values)
    finite = stacked[np.isfinite(stacked)]
    if finite.size == 0:                                   # pragma: no cover
        return None
    lower, upper = float(finite.min()), float(finite.max())
    if upper <= lower:                                     # pragma: no cover
        return None
    pad = 0.05 * (upper - lower)
    return (lower - pad, upper + pad)


def layer_labels_from_model(model, nlay: int) -> list[str] | None:
    """Layer names carried by the model's own build context, if it has any.

    `ModelSpec.build` stores the `ModelContext` on the model, so a model declared
    with `ModelContext(surfaces=stack.build(vor))` knows what its layers are
    CALLED -- "sand", "clay" -- rather than only how many there are. A model
    built the imperative way, or one whose `surfaces` is a plain frame, carries
    no names and gets the `Layer N` default.
    """

    names = getattr(getattr(getattr(model, "myflopy_context", None), "surfaces", None),
                    "names", None)
    if names is None:
        return None
    names = list(names)
    return names if len(names) == nlay else None


#: Friendly placements, mapped onto what matplotlib's `loc` actually accepts.
#: `"bottom"` is the obvious word and is not one of its values; `"outside ..."`
#: is not a placement it has at all, and is the one that matters on a section --
#: a layer legend sits on top of the geology everywhere inside the axes.
_LEGEND_ALIASES = {
    "auto": "best",
    "top": "upper center",
    "bottom": "lower center",
    "left": "center left",
    "right": "center right",
    "topleft": "upper left",
    "topright": "upper right",
    "bottomleft": "lower left",
    "bottomright": "lower right",
}

#: `(loc, bbox_to_anchor)` for placements that sit OUTSIDE the axes.
_LEGEND_OUTSIDE = {
    "right": ("center left", (1.02, 0.5)),
    "left": ("center right", (-0.02, 0.5)),
    "bottom": ("upper center", (0.5, -0.08)),
    "top": ("lower center", (0.5, 1.02)),
}


def legend_placement(legend):
    """``(show, loc, bbox_to_anchor)`` from a friendly ``legend=`` value.

    Accepts ``True``/``False``/``None``, ``"auto"``, a plain side (``"bottom"``,
    ``"left"``, ``"right"``, ``"top"``), a corner (``"topright"`` or matplotlib's
    own ``"upper right"``), or ``"outside <side>"`` to put it beside the axes.
    Anything matplotlib's ``loc`` understands passes straight through, so this
    narrows nothing.

    Raises on an unknown word rather than falling back to ``"best"``: a legend
    that silently ignored where you told it to go is the defect this whole
    section's worth of work has been about.
    """

    if legend is None or legend is False:
        return False, None, None
    if legend is True:
        return True, None, None

    text = str(legend).strip().lower()
    if text.startswith("outside"):
        side = text.replace("outside", "").strip().replace("_", "") or "right"
        if side not in _LEGEND_OUTSIDE:
            raise ValueError(
                f"outside legend must name {sorted(_LEGEND_OUTSIDE)}, not {side!r}."
            )
        return (True, *_LEGEND_OUTSIDE[side])

    key = text.replace("_", "").replace(" ", "")
    if key in _LEGEND_ALIASES:
        return True, _LEGEND_ALIASES[key], None
    known = {
        "best", "upper right", "upper left", "lower left", "lower right",
        "right", "center left", "center right", "lower center", "upper center",
        "center",
    }
    if text in known:
        return True, text, None
    raise ValueError(
        f"legend={legend!r} is not a placement. Use True/False, 'auto', a side "
        f"({', '.join(sorted(_LEGEND_ALIASES))}), 'outside <side>', or one of "
        f"matplotlib's own: {', '.join(sorted(known))}."
    )


def mirror_horizontally(axes) -> None:
    """Reflect the section about its own centre, so it reads from the other end.

    The filled section cannot be reversed by reversing its LINE: FloPy's
    `PlotCrossSection` normalizes direction internally, and a forward and a
    reversed line produce byte-identical geometry (measured -- same 31 nodes,
    same `xcenters`, `direction` 'x' either way). The profile branch DOES honour
    line order, because `XSection` samples the line itself.

    Done by inverting the VIEW and relabelling the axis, not by rewriting the
    drawn geometry. Two attempts at the latter failed on matplotlib internals --
    `PatchCollection.get_paths()` hands back paths derived from its patches, so
    writing to `path.vertices` is silently discarded, and its `set_paths`
    expects patches rather than paths. Inverting the limits needs neither, works
    for every artist including the ones added later, and leaves the data in one
    consistent frame.

    The tick labels are rewritten to measure from the NEW start, so the
    numbering still counts up left to right: a feature that was `d` from the old
    start is `L - d` from this one, and saying otherwise would make the two
    pictures incomparable.

    Mirrors the CURRENT x view and does not refit it. The refit belongs to
    whichever caller needs one: the filled branch mirrors before its limits are
    set, so it autoscales first (a `PatchCollection` does not update the data
    limits on its own). The grid's geometry section arrives already framed, and
    autoscaling it here added matplotlib's 5% margins -- widening the picture as
    a side effect of reversing it, which is not what the argument says it does.
    """

    from matplotlib.ticker import FuncFormatter

    low, high = axes.get_xlim()
    total = low + high
    axes.set_xlim(high, low)
    axes.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _pos: f"{total - value:g}")
    )


def _cross_section(flopy_model, modelgrid, line_spec):
    """FloPy's ``PlotCrossSection``, with a readable error when the line misses.

    FloPy raises a BARE ``StopIteration`` here -- its own guard reads
    ``next(iter(self.xypts.values()))`` to build the message, which is exactly
    the call that fails when no cell was intersected. The traceback names
    `crosssection.py` and says nothing at all, so a mistyped line (``zip(xy[0],
    xy[0])`` pairs X with X) looks like a library fault.
    """

    try:
        return PlotCrossSection(
            model=flopy_model, modelgrid=modelgrid, line=line_spec
        )
    except StopIteration as error:
        points = np.asarray(line_spec.get("line", []), dtype=float)
        extent = modelgrid.extent
        where = (
            f"the line spans x {points[:, 0].min():.6g}..{points[:, 0].max():.6g}, "
            f"y {points[:, 1].min():.6g}..{points[:, 1].max():.6g}"
            if points.ndim == 2 and len(points) else "the line has no points"
        )
        raise ValueError(
            f"the section line crosses no cells, so there is nothing to draw. "
            f"The grid spans x {extent[0]:.6g}..{extent[1]:.6g}, "
            f"y {extent[2]:.6g}..{extent[3]:.6g}, and {where}. "
            f"A line built as `zip(xy[0], xy[0])` pairs X with X and will miss "
            f"the grid; `zip(xy[0], xy[1])` is the one you want."
        ) from error


def plot_layered_cross_section(
    modelgrid,
    line,
    *,
    ax: Axes | None = None,
    style: ModelCrossSectionStyle | None = None,
    layer_colors: Sequence[str] | None = None,
    layer_labels: Sequence[str] | None = None,
    head_surface=None,
    head_surfaces=None,
    layers=None,
    values=None,
    values_cmap: str = "viridis",
    values_label: str | None = None,
    colorbar: bool = True,
    flopy_model=None,
    ylim: tuple[float, float] | None = None,
    xlim: tuple[float, float] | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    show_grid: bool = True,
    show_layers: bool = True,
    show_head: bool = True,
    show_legend: bool = True,
    legend_loc: str | None = None,
    legend_anchor=None,
    reverse: bool = False,
):
    """Render a layer-colored cross section from any flopy modelgrid + line.

    This is the grid-level rendering core shared by
    :func:`plot_model_cross_section` (model-aware) and
    ``myflopy.layers.LayerBuildResult.cross_section`` (pre-model, off raw
    ``top``/``botm``/``idomain`` arrays). It colors cells by layer index, draws
    the grid, an optional head/water surface line, and a layer legend.

    Parameters
    ----------
    modelgrid
        Any flopy modelgrid exposing ``nlay`` and ``ncpl`` (e.g. a ``VertexGrid``).
    line
        Section line accepted by ``flopy.plot.PlotCrossSection`` -- a dict, a
        shapely ``LineString``, or a sequence of ``(x, y)`` pairs.
    head_surface
        Optional 1-D array drawn as a line over the section (e.g. simulated heads).
    head_surfaces
        Optional ``[(label, array), ...]`` drawn as several labelled surfaces --
        one water level per layer, say. Takes precedence over ``head_surface``.
    layers
        Zero-based layers to DRAW. With none, every layer. Cells outside the
        selection are masked out and the vertical extent is cropped to what
        remains, because leaving the axis at full height puts the layers you
        asked for in a thin band with empty space above and below.
    values
        Optional ``(nlay, ncpl)`` array to colour the CELLS by -- heads,
        concentration, K, a zone id. Replaces the layer colouring rather than
        adding to it: a cell has one fill, and drawing both would mean the legend
        and the colorbar describe the same patch differently.
    values_cmap
        Matplotlib colormap for ``values``.
    values_label
        Colorbar label for ``values``.
    colorbar
        Draw the colorbar for ``values``. Ignored without ``values``.
    flopy_model
        Optional flopy model passed through to ``PlotCrossSection`` (only needed
        for model-aware head rendering).

    Returns ``(fig, ax)``.
    """
    style = ModelCrossSectionStyle() if style is None else style
    line_spec = normalize_line(line)
    nlay = modelgrid.nlay
    chosen = _resolve_drawn_layers(layers, nlay)
    layer_id, values = _masked_to_layers(modelgrid, values, chosen)
    resolved_layer_colors = _resolve_layer_colors(nlay, layer_colors, style)
    resolved_layer_labels = _resolve_layer_labels(nlay, layer_labels)

    with _theme_context(style):
        if ax is None:
            fig, ax = mpl_axes(figsize=style.figsize)
        else:
            fig = ax.figure

        xsect = _cross_section(flopy_model, modelgrid, line_spec)

        # `plot_grid` draws EVERY cell's edges and takes no layer filter, so with
        # a subset it outlined the layers just masked out of the fill -- a bottom
        # layer excluded from `layers=` still appeared, as outlines. Cropping the
        # y-axis cannot fix that: a layer's elevation range overlaps its
        # neighbours' (measured, the canonical bottom layer spans -1.3..59.8
        # inside a retained band of 28.5..166.0). So with a subset the edges come
        # from the FILL collection itself, which is already masked.
        if show_grid and chosen is None:
            xsect.plot_grid(ax=ax, linewidths=style.grid_linewidth, color=style.grid_color)
        edges = (
            {"edgecolor": style.grid_color, "linewidth": style.grid_linewidth}
            if show_grid and chosen is not None else {}
        )

        if values is not None:
            # A continuous field REPLACES the layer colouring: one fill per cell,
            # so the legend and the colorbar cannot describe the same patch two
            # different ways. The layer legend is suppressed for the same reason.
            painted = xsect.plot_array(
                np.asarray(values, dtype=float),
                ax=ax,
                cmap=values_cmap,
                alpha=style.layer_alpha,
                **edges,
            )
            if colorbar:
                bar = fig.colorbar(painted, ax=ax, fraction=0.025, pad=0.02)
                if values_label:
                    bar.set_label(values_label, size=style.label_fontsize)
        elif show_layers:
            xsect.plot_array(
                layer_id,
                ax=ax,
                cmap=ListedColormap(resolved_layer_colors),
                vmin=0,
                vmax=nlay - 1,
                alpha=style.layer_alpha,
                **edges,
            )

        drawn_surfaces = _resolve_head_surfaces(head_surfaces, head_surface, style)
        if show_head:
            surface_colours = _surface_colours([n for n, _ in drawn_surfaces], style)
            for (_label, surface), colour in zip(
                drawn_surfaces, surface_colours, strict=False
            ):
                xsect.plot_surface(
                    surface, ax=ax, color=colour, lw=style.head_linewidth,
                )

        if reverse:
            # Before the limits are set, so the crop below measures the mirrored
            # drawing rather than the original -- which is also why the refit
            # is here: at this point the axes carry only collections, and those
            # leave the data limits alone.
            ax.relim()
            ax.autoscale_view(scalex=True, scaley=False)
            mirror_horizontally(ax)

        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)
        elif chosen is not None:
            # Only for a SUBSET. Drawing every layer keeps FloPy's own framing,
            # which `plot_model_cross_section` and every figure built on it
            # already use -- tightening that unasked would move them all.
            bounds = _drawn_elevation_bounds(ax)
            if bounds is not None:
                ax.set_ylim(*bounds)

        ax.set_xlabel(xlabel or style.xlabel, size=style.label_fontsize)
        ax.set_ylabel(ylabel or style.ylabel, size=style.label_fontsize)
        ax.tick_params(axis="both", which="major", labelsize=style.tick_labelsize)

        if title or style.title:
            ax.set_title(title or style.title, fontsize=style.title_fontsize)

        if show_legend and values is None:
            shown = chosen if chosen is not None else range(nlay)
            handles = _build_legend_handles(
                layer_colors=[resolved_layer_colors[k] for k in shown],
                layer_labels=[resolved_layer_labels[k] for k in shown],
                surfaces=list(zip(
                    [label for label, _ in drawn_surfaces],
                    _surface_colours([n for n, _ in drawn_surfaces], style),
                    strict=False,
                )) if show_head else [],
                head_linewidth=style.head_linewidth,
            )
            ax.legend(
                handles=handles,
                loc=legend_loc or style.legend_loc,
                frameon=style.legend_frameon,
                fontsize=style.legend_fontsize,
                labelcolor=(
                    style.legend_labelcolor
                    or readable_on(fig.get_facecolor())
                ),
                **({"bbox_to_anchor": legend_anchor} if legend_anchor else {}),
            )
            if legend_anchor:
                # An outside legend needs room made for it, or it is drawn off
                # the canvas edge and saved figures clip it.
                fig.tight_layout()

    return fig, ax


#: What `section(fill=...)` accepts as a named fill.
FILL_KINDS = ("layer", "results")


def section_values(model, fill, *, per=None, kstpkper=None):
    """The ``(nlay, ncpl)`` array a ``fill=`` asks for, or None to colour by layer.

    ``"results"`` resolves through the model's own field reader rather than
    ``.hds``, so it is heads on GWF, concentration on GWT and temperature on GWE
    -- the rule ledger 99 established for cross-sections generally. Anything
    array-like is taken as given.
    """

    if isinstance(fill, str):
        if fill == "layer":
            return None
        if fill not in FILL_KINDS:
            raise ValueError(
                f"fill must be one of {FILL_KINDS} or a per-cell array, not {fill!r}."
            )
        reader = getattr(model, "_field_reader", None)
        if reader is None:
            raise ValueError(
                "fill='results' needs a model with results; a bare grid has "
                "none. Use fill='layer' for the geometry, or pass an array."
            )
        nlay = model.gwf.modelgrid.nlay
        return np.asarray(
            [reader.array(layer=k, per=per, kstpkper=kstpkper) for k in range(nlay)],
            dtype=float,
        )

    values = np.asarray(fill, dtype=float)
    if values.ndim == 1:
        raise ValueError(
            f"fill= needs one value per cell PER LAYER -- a (nlay, ncpl) array. "
            f"Got a flat {values.shape} array, which would colour every layer "
            f"the same. Stack the layers, or pass fill='results'."
        )
    return values


def section_time(model, per, kstpkper):
    """The ``(kstp, kper)`` a filled section reads, from ``per``/``kstpkper``/the end."""

    if kstpkper is not None:
        return tuple(kstpkper)
    times = list(model.kstpkper)
    if per is None:
        return times[-1]
    matching = [t for t in times if int(t[1]) == int(per)]
    if not matching:
        raise ValueError(
            f"per={per} has no saved output; this model wrote periods "
            f"{sorted({int(t[1]) for t in times})}."
        )
    return max(matching)


def filled_section(
    modelgrid,
    line,
    *,
    model=None,
    fill="layer",
    layers=None,
    head_layers=0,
    layer_labels=None,
    show_grid: bool = True,
    legend=True,
    reverse: bool = False,
    per=None,
    kstpkper=None,
    cmap: str = "viridis",
    label: str | None = None,
    title: str | None = None,
    backend: str = "mpl",
):
    """The shared renderer behind ``section(fill=...)`` at every scope.

    Lives here, at layer 1, so both `myflopy.plot.section` (layer 7) and
    `GridPlots.section` (layer 3) reach it DOWNWARD. The obvious alternative --
    the grid scope calling up into `myflopy.plot` -- is an upward import the
    layering ratchet rightly refuses, and it would have been the only one in the
    package.

    ``model`` is optional and duck-typed: with one, ``fill="results"`` can read a
    field and ``fill="layer"`` draws the simulated head as a water surface over
    the geology; without one (a bare grid) neither is available and
    ``fill="results"`` raises rather than quietly degrading to ``"layer"``.

    Every argument but ``backend`` is :class:`FilledSection`'s, documented
    there. ``backend="mpl"`` returns a Matplotlib ``Figure``; ``"plotly"``
    returns the :class:`FilledSection` picture itself. Both come from one
    picture, so the two renderers cannot be handed different arguments.
    """

    if backend not in ("mpl", "plotly"):
        raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")
    picture = FilledSection(
        modelgrid, line, model=model, fill=fill, layers=layers,
        head_layers=head_layers, layer_labels=layer_labels, show_grid=show_grid,
        legend=legend, reverse=reverse, per=per, kstpkper=kstpkper, cmap=cmap,
        label=label, title=title,
    )
    return picture.plot_mpl() if backend == "mpl" else picture


#: Where a Plotly legend sits for each placement matplotlib's ``loc`` names:
#: ``(x, y, xanchor, yanchor)`` in the plot area's own coordinates, which is
#: what matplotlib's axes coordinates are. ``"best"`` has no Plotly equivalent
#: and is left to Plotly's own default.
_PLOTLY_LEGEND_INSIDE = {
    "upper right": (1.0, 1.0, "right", "top"),
    "upper left": (0.0, 1.0, "left", "top"),
    "lower left": (0.0, 0.0, "left", "bottom"),
    "lower right": (1.0, 0.0, "right", "bottom"),
    "right": (1.0, 0.5, "right", "middle"),
    "center left": (0.0, 0.5, "left", "middle"),
    "center right": (1.0, 0.5, "right", "middle"),
    "lower center": (0.5, 0.0, "center", "bottom"),
    "upper center": (0.5, 1.0, "center", "top"),
    "center": (0.5, 0.5, "center", "middle"),
}

#: The four ``"outside <side>"`` placements. Farther out than matplotlib's
#: anchors on the left and bottom, because Plotly's axis titles sit in the
#: margin where matplotlib's anchor would put the legend on top of them.
_PLOTLY_LEGEND_OUTSIDE = {
    "right": {"x": 1.02, "y": 0.5, "xanchor": "left", "yanchor": "middle"},
    "left": {"x": -0.12, "y": 0.5, "xanchor": "right", "yanchor": "middle"},
    "bottom": {"x": 0.5, "y": -0.18, "xanchor": "center", "yanchor": "top",
               "orientation": "h"},
    "top": {"x": 0.5, "y": 1.02, "xanchor": "center", "yanchor": "bottom",
            "orientation": "h"},
}


def _plotly_legend(loc, anchor, style) -> dict | None:
    """The ``layout.legend`` for a placement :func:`legend_placement` resolved.

    Driven by the SAME resolved ``(loc, anchor)`` the Matplotlib renderer
    receives, so one ``legend=`` word means one placement on both backends and an
    unknown word has already been refused. ``None`` means Plotly's own default.
    """

    if anchor is not None:
        x, y = anchor
        side = "right" if x > 1 else "left" if x < 0 else "bottom" if y < 0 else "top"
        return dict(_PLOTLY_LEGEND_OUTSIDE[side])
    loc = loc or style.legend_loc
    if loc == "best":
        return None
    x, y, xanchor, yanchor = _PLOTLY_LEGEND_INSIDE[loc]
    # Inside the plot the legend sits on the geology it describes. Matplotlib's
    # is frameless; an opaque Plotly box would hide the cells under it.
    return {"x": x, "y": y, "xanchor": xanchor, "yanchor": yanchor,
            "bgcolor": "rgba(255,255,255,0.6)"}


def _plotly_rgba(rgba) -> str:
    """A Matplotlib RGBA tuple as Plotly's ``rgba(r,g,b,a)``."""

    red, green, blue, alpha = (float(c) for c in rgba)
    return (
        f"rgba({round(red * 255)},{round(green * 255)},{round(blue * 255)},"
        f"{alpha:.6g})"
    )


def _points_to_px(points: float) -> float:
    """A Matplotlib line width in points as a Plotly width in pixels (96 dpi).

    Floored at half a pixel: the house hairline is 0.2 pt, which browsers render
    as nothing at all.
    """

    return max(0.5, float(points) * 96.0 / 72.0)


class FilledSection(Picture):
    """The cells a section line cuts, coloured by layer or by a field.

    The picture behind ``section(fill=...)`` at every scope, on BOTH renderers:
    ``.fig`` is the interactive Plotly figure, and ``.plot_mpl()`` draws the same
    section with Matplotlib through :func:`plot_layered_cross_section`. The two
    read ONE FloPy ``PlotCrossSection`` -- its ``polygons`` for the cells and its
    ``projpts`` for the water levels -- and ONE resolved set of arrays, colours
    and labels, so they draw the same cells in the same colours with the same
    water levels. Plotly re-derives none of the geometry; a second derivation is
    what would let the two drift.

    Everything is resolved when the picture is built, so a bad argument raises
    at the call rather than at the first ``.fig``: the line (one that crosses no
    cell), ``fill``, ``layers``, ``head_layers``, ``legend`` and ``cmap``.

    What Plotly adds is hover: every drawn cell names its layer, its cell number
    and its top and bottom along the line, plus its value under a field fill.

    Parameters
    ----------
    modelgrid : flopy.discretization.VertexGrid
        The grid to section. Needs ``nlay`` and ``ncpl``.
    line : LineString, dict, sequence of (x, y), or path
        The section line -- anything :func:`normalize_line` accepts: a shapely
        ``LineString``, a flopy ``{"line": [...]}`` dict, ``[(x, y), ...]``, or
        a vector file to read one from.
    model : SimulationBase, optional
        Whose results ``fill="results"`` paints and whose water levels are drawn
        over ``fill="layer"``. Without one -- a bare grid -- there are no
        results: ``"layer"`` draws the geometry alone and ``"results"`` raises.
    fill : {'layer', 'results'} or array-like, default 'layer'
        What colours the cells. ``"layer"`` colours each cell by its layer and
        draws the water levels chosen by ``head_layers`` on top. ``"results"``
        colours them by the model's own dependent variable -- heads on GWF,
        concentration on GWT, temperature on GWE -- with a colorbar. An
        ``(nlay, ncpl)`` array colours them by anything else (K, a zone id, a
        residual); a flat ``(ncpl,)`` array raises rather than painting every
        layer the same. A cell whose value is NaN (a dry or inactive cell under
        ``"results"``) is left unfilled.
    layers : int or sequence of int, optional
        Which layers' cells to draw, zero-based; every layer by default. The rest
        are masked out and the vertical extent is cropped to what is drawn. A
        layer outside the model raises.
    head_layers : int, sequence of int, or None, default 0
        Whose water levels to draw over ``fill="layer"``: one line per layer,
        labelled when there are several. ``None`` draws none. Ignored under
        ``"results"`` or an array, which already paint the field onto the cells.
    layer_labels : sequence of str, optional
        Legend names for the layers, one per layer in model order. Defaults to
        the names the model's build context carries, else ``Layer 1..N``.
    show_grid : bool, default True
        Draw the cell edges.
    legend : bool or str, default True
        Where the layer legend goes: ``False`` drops it, ``True`` puts it at the
        upper right, ``"auto"`` leaves it to the renderer, a side or corner
        (``"bottom"``, ``"topleft"``, matplotlib's ``"upper right"``) places it
        inside, and ``"outside right"``/``"outside bottom"`` (also ``left`` and
        ``top``) place it beside the plot. An unknown word raises.
    reverse : bool, default False
        Draw the line from its other end, measuring distance from the new start.
    per : int, optional
        Stress period (zero-based) of the results painted or drawn. Defaults to
        the last saved time.
    kstpkper : tuple of (int, int), optional
        Exact ``(time step, period)``; takes precedence over ``per``.
    cmap : str, default 'viridis'
        A MATPLOTLIB colormap name for ``"results"`` or an array fill. Both
        renderers colour from its lookup table, so a name Matplotlib does not
        know raises here. Unused by ``fill="layer"``, which has its own palette.
    label : str, optional
        Colorbar title for ``"results"`` or an array fill.
    title : str, optional
        The figure title.
    """

    #: Class attributes, so an instance built with `object.__new__` (a test
    #: double) still answers the idempotency contract, as `XSection` does.
    _assembled = False
    _fig = None

    def __init__(
        self,
        modelgrid,
        line,
        *,
        model=None,
        fill="layer",
        layers=None,
        head_layers=0,
        layer_labels=None,
        show_grid: bool = True,
        legend=True,
        reverse: bool = False,
        per=None,
        kstpkper=None,
        cmap: str = "viridis",
        label: str | None = None,
        title: str | None = None,
    ):
        from matplotlib import colormaps

        self._grid = modelgrid
        self._line = normalize_line(line)
        self._model = model
        self._style = ModelCrossSectionStyle()
        self._show_legend, self._legend_loc, self._legend_anchor = legend_placement(legend)
        self._values = section_values(model, fill, per=per, kstpkper=kstpkper)
        names = layer_labels
        if names is None and model is not None:
            names = layer_labels_from_model(model, modelgrid.nlay)
        self._layer_labels = names
        # The water levels are named from the SAME list as the layers, so
        # `layer_labels=["sand", ...]` reaches "Head, sand" too -- it stopped at
        # the layer entries before, and the lines read "Head, Layer 1".
        self._surfaces = section_head_surfaces(
            model, head_layers, per=per, kstpkper=kstpkper,
            names=_resolve_layer_labels(modelgrid.nlay, names) if names is not None else None,
        ) if self._values is None else []
        self._layers = layers
        self._chosen = _resolve_drawn_layers(layers, modelgrid.nlay)
        if cmap not in colormaps:
            raise ValueError(
                f"cmap={cmap!r} is not a Matplotlib colormap; both renderers of a "
                f"filled section colour from Matplotlib's lookup table. Try "
                f"'viridis', 'magma', 'RdBu' or 'terrain'."
            )
        self._cmap = cmap
        self._label = label
        self._title = title
        self._show_grid = bool(show_grid)
        self._reverse = bool(reverse)
        self._flopy_model = getattr(model, "gwf", None)
        # Built now: it is what refuses a line that crosses no cell, and a
        # picture that raised on first `.fig` instead of at the call would be
        # the lazy-failure defect `Picture` exists to rule out.
        self._xsect = _cross_section(self._flopy_model, modelgrid, self._line)

    def plot_mpl(self, ax=None):
        """Draw this section with Matplotlib and return the ``Figure``.

        The renderer ``section(fill=..., backend="mpl")`` uses; a fresh figure
        per call, like ``Choro.plot_mpl``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Draw into this Axes rather than a new figure.

        Returns
        -------
        matplotlib.figure.Figure
        """

        figure, _axes = plot_layered_cross_section(
            self._grid,
            self._line,
            ax=ax,
            values=self._values,
            values_cmap=self._cmap,
            values_label=self._label,
            layers=self._layers,
            layer_labels=self._layer_labels,
            show_grid=self._show_grid,
            reverse=self._reverse,
            show_legend=self._show_legend,
            legend_loc=self._legend_loc,
            legend_anchor=self._legend_anchor,
            head_surfaces=self._surfaces,
            flopy_model=self._flopy_model,
            title=self._title,
            show_head=bool(self._surfaces),
        )
        return figure

    @property
    def fig(self) -> Fig:
        """The Plotly figure, assembled once and then reused."""

        if not self._assembled:
            self._fig = self._assemble()
            self._assembled = True
        return self._fig

    def _assemble(self) -> Fig:
        """Build the Plotly figure from the SAME inputs `plot_mpl` draws.

        Mirrors :func:`plot_layered_cross_section` step for step: the grid edges
        under the fill when every layer is drawn (FloPy's ``plot_grid`` outlines
        EVERY cut cell, filled or not), the edges on the fill itself when a
        subset is (so excluded layers are not outlined), the water levels from
        FloPy's own ``projpts`` exactly as ``plot_surface`` walks them, the
        extent FloPy frames, and ``reverse`` as distance from the other end.
        """

        from matplotlib import colormaps
        from matplotlib.colors import Normalize, to_rgba

        style, grid, xsect = self._style, self._grid, self._xsect
        nlay, ncpl = int(grid.nlay), int(grid.ncpl)
        xmin, xmax, ymin, ymax = (float(v) for v in xsect.extent)
        mirror = xmin + xmax

        def along(x):
            """Distance along the line as DRAWN -- from the far end when reversed."""

            x = np.asarray(x, dtype=float)
            return mirror - x if self._reverse else x

        layer_id, values = _masked_to_layers(grid, self._values, self._chosen)
        painted = np.ravel(np.asarray(values if values is not None else layer_id, dtype=float))
        polygons = [
            (int(node), np.asarray(polygon.get_xy(), dtype=float))
            for node, cell_polygons in sorted(xsect.polygons.items())
            for polygon in cell_polygons
        ]
        filled = [
            (node, vertices) for node, vertices in polygons
            if node < painted.size and np.isfinite(painted[node])
        ]
        labels = _resolve_layer_labels(nlay, self._layer_labels)
        edge_colour = _plotly_rgba(to_rgba(style.grid_color))
        edge_width = _points_to_px(style.grid_linewidth)
        subset = self._chosen is not None

        def outlines(pieces):
            """``(x, y)`` with each polygon closed and separated by ``None``."""

            xs, ys = [], []
            for _node, vertices in pieces:
                xs += [*along(vertices[:, 0]), None]
                ys += [*vertices[:, 1], None]
            return xs, ys

        figure = Fig()
        edges = (
            {"color": edge_colour, "width": edge_width}
            if self._show_grid and subset else {"width": 0}
        )

        if values is not None:
            drawn = np.array([painted[node] for node, _ in filled], dtype=float)
            colormap = colormaps[self._cmap]
            norm = Normalize()
            norm.autoscale_None(drawn)
            groups: dict[str, list] = {}
            for (node, vertices), value in zip(filled, drawn, strict=True):
                colour = _plotly_rgba(to_rgba(colormap(norm(value)), style.layer_alpha))
                groups.setdefault(colour, []).append((node, vertices))
            for colour, pieces in groups.items():
                xs, ys = outlines(pieces)
                figure.add_scatter(
                    x=xs, y=ys, mode="lines", fill="toself", fillcolor=colour,
                    line=edges, hoverinfo="skip", showlegend=False, meta="cells",
                )
            stops = [
                [index / (colormap.N - 1), _plotly_rgba(colormap(index))]
                for index in range(colormap.N)
            ]
            figure.add_scatter(
                x=[None, None], y=[None, None], mode="markers", hoverinfo="skip",
                showlegend=False, meta="colorbar",
                marker={
                    "color": [norm.vmin, norm.vmax], "cmin": norm.vmin,
                    "cmax": norm.vmax, "colorscale": stops, "showscale": True,
                    "colorbar": {"title": {"text": self._label or ""}},
                },
            )
        else:
            palette = ListedColormap(_resolve_layer_colors(nlay, None, style))
            norm = Normalize(vmin=0, vmax=nlay - 1)
            shown = self._chosen if subset else range(nlay)
            for layer in shown:
                pieces = [(n, v) for n, v in filled if n // ncpl == layer]
                xs, ys = outlines(pieces) if pieces else ([None], [None])
                figure.add_scatter(
                    x=xs, y=ys, mode="lines", fill="toself",
                    fillcolor=_plotly_rgba(to_rgba(palette(norm(layer)), style.layer_alpha)),
                    line=edges, hoverinfo="skip", name=labels[layer],
                    legendgroup=f"layer {layer}", showlegend=self._show_legend,
                    meta="cells",
                )

        if self._show_grid and not subset:
            # Every cut cell, filled or not -- FloPy's `plot_grid` outlines them
            # all, so a dry cell under `fill="results"` keeps its edges. Drawn
            # OVER the fill: under an 85%-opaque fill a hairline does not show.
            xs, ys = outlines(polygons)
            figure.add_scatter(
                x=xs, y=ys, mode="lines", line={"color": edge_colour, "width": edge_width},
                hoverinfo="skip", showlegend=False, name="cell edges", meta="grid",
            )

        surface_colours = _surface_colours([name for name, _ in self._surfaces], style)
        extents = {
            node: (float(np.min(np.asarray(verts)[:, 0])), float(np.max(np.asarray(verts)[:, 0])))
            for node, verts in sorted(xsect.projpts.items())
        }
        drawn_levels = []
        for (name, surface), colour in zip(self._surfaces, surface_colours, strict=False):
            level = np.ravel(np.asarray(surface, dtype=float))
            xs, ys = [], []
            for node, (left, right) in extents.items():
                if node >= level.size or np.isnan(level[node]):  # FloPy's own test
                    continue
                xs += [*along([left, right]), None]
                ys += [level[node], level[node], None]
                drawn_levels.append(level[node])
            figure.add_scatter(
                x=xs, y=ys, mode="lines",
                line={"color": colour, "width": _points_to_px(style.head_linewidth)},
                name=name, showlegend=self._show_legend, meta="surface",
                hovertemplate=f"{name}<br>%{{y:.2f}}<extra></extra>",
            )

        if filled:
            rows = []
            for node, vertices in filled:
                corners = vertices[:-1] if len(vertices) > 1 and np.allclose(
                    vertices[0], vertices[-1]) else vertices
                rows.append((
                    float(np.mean(along(corners[:, 0]))), float(np.mean(corners[:, 1])),
                    labels[node // ncpl], node % ncpl,
                    float(np.max(corners[:, 1])), float(np.min(corners[:, 1])),
                    float(painted[node]),
                ))
            valued = values is not None
            name = self._label or "value"
            figure.add_scatter(
                x=[r[0] for r in rows], y=[r[1] for r in rows], mode="markers",
                marker={"size": 10, "opacity": 0}, showlegend=False, meta="hover",
                customdata=[r[2:] for r in rows],
                hovertemplate=(
                    "%{customdata[0]}, cell %{customdata[1]}<br>"
                    "top %{customdata[2]:.2f}, bottom %{customdata[3]:.2f}"
                    + (f"<br>{name} %{{customdata[4]:.4g}}" if valued else "")
                    + "<extra></extra>"
                ),
            )

        if subset:
            bounds = _padded_elevation_bounds(
                [vertices[:, 1] for _node, vertices in filled]
                + ([np.asarray(drawn_levels, dtype=float)] if drawn_levels else [])
            )
            y_range = list(bounds) if bounds is not None else [ymin, ymax]
        else:
            y_range = [ymin, ymax]
        legend = _plotly_legend(self._legend_loc, self._legend_anchor, style)
        figure.update_layout(
            title=self._title or style.title,
            xaxis={"title": {"text": style.xlabel}, "range": [xmin, xmax]},
            yaxis={"title": {"text": style.ylabel}, "range": y_range},
            showlegend=self._show_legend and values is None,
            hovermode="closest",
            **({"legend": legend} if legend else {}),
        )
        return figure


def section_head_surfaces(model, head_layers, *, per=None, kstpkper=None, names=None) -> list:
    """``[(label, array)]`` -- one water surface per requested layer.

    ``head_layers=None`` draws none, an int or a list draws those layers. The
    default is layer 0, which is what a single water table means and what every
    figure drawn before 2026-09-02 got.

    Labelled per layer rather than "Simulated Head" once, because several
    unlabelled lines on one section are unreadable -- and the label carries the
    layer NAME: ``names`` (one per layer, the caller's ``layer_labels``) when
    given, else the names the model's build context knows, else ``Layer N``.
    """

    if model is None or head_layers is None:
        return []
    reader = getattr(model, "_field_reader", None)
    if reader is None:
        return []
    nlay = model.gwf.modelgrid.nlay
    chosen = _resolve_drawn_layers(head_layers, nlay)
    if chosen is None:                                     # pragma: no cover
        chosen = list(range(nlay))
    names = list(names) if names is not None else (
        layer_labels_from_model(model, nlay) or [f"Layer {k + 1}" for k in range(nlay)]
    )
    return [
        (f"Head, {names[k]}" if len(chosen) > 1 else "Simulated Head",
         reader.array(layer=k, per=per, kstpkper=kstpkper))
        for k in chosen
    ]


def plot_model_cross_section(
    model: SimulationBase,
    line,
    *,
    kstpkper: tuple | None = None,
    head_data=None,
    head_layer: int = 0,
    ax: Axes | None = None,
    style: ModelCrossSectionStyle | None = None,
    ylim: tuple[float, float] | None = None,
    xlim: tuple[float, float] | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    layer_colors: Sequence[str] | None = None,
    layer_labels: Sequence[str] | None = None,
    show_grid: bool = True,
    show_layers: bool = True,
    show_head: bool = True,
    show_legend: bool = True,
):
    """
    Build a report-oriented matplotlib cross section for a MODFLOW model.

    Resolves heads from the model and delegates the rendering to
    :func:`plot_layered_cross_section`.

    Parameters
    ----------
    model
        Parent simulation/model wrapper containing ``gwf``.
    line
        Cross-section line specification accepted by ``flopy.plot.PlotCrossSection``.
        This may be a dict, a shapely ``LineString``, or a sequence of ``(x, y)``
        pairs.
    kstpkper
        Stress-period / time-step tuple used to read heads from the model output
        when ``head_data`` is not provided.
    head_data
        Optional explicit head array. When provided it takes precedence over
        ``kstpkper``.
    head_layer
        Zero-based layer index used when selecting the head surface from a
        multi-layer head array.
    ax
        Existing matplotlib axes to draw into. When omitted a new figure is
        created.
    style
        Optional style preset for report plotting.
    """
    head_array = _resolve_head_data(model, kstpkper=kstpkper, head_data=head_data)
    head_surface = _resolve_head_surface(head_array, head_layer=head_layer)
    return plot_layered_cross_section(
        model.gwf.modelgrid,
        line,
        ax=ax,
        style=style,
        layer_colors=layer_colors,
        layer_labels=layer_labels,
        head_surface=head_surface,
        flopy_model=model.gwf,
        ylim=ylim,
        xlim=xlim,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        show_grid=show_grid,
        show_layers=show_layers,
        show_head=show_head,
        show_legend=show_legend,
    )
