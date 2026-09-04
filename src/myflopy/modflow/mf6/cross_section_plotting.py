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
from myflopy.viz import category_colors, mpl_axes

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


def _normalize_line(line) -> dict:
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
        return _normalize_line(geometry)

    if isinstance(line, shp.MultiLineString):
        merged = shp.line_merge(line)
        return _normalize_line(
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
    line_spec = _normalize_line(line)
    nlay = modelgrid.nlay
    layer_id = _resolve_layer_array(modelgrid)
    chosen = _resolve_drawn_layers(layers, nlay)
    if chosen is not None:
        # Masked, not dropped: FloPy's `plot_array` skips NaN cells, which keeps
        # every other cell's geometry exactly where it was.
        keep = np.zeros(nlay, dtype=bool)
        keep[chosen] = True
        layer_id = np.where(keep[:, None], layer_id.astype(float), np.nan)
        if values is not None:
            values = np.where(keep[:, None], np.asarray(values, dtype=float), np.nan)
    resolved_layer_colors = _resolve_layer_colors(nlay, layer_colors, style)
    resolved_layer_labels = _resolve_layer_labels(nlay, layer_labels)

    with _theme_context(style):
        if ax is None:
            fig, ax = mpl_axes(figsize=style.figsize)
        else:
            fig = ax.figure

        xsect = PlotCrossSection(model=flopy_model, modelgrid=modelgrid, line=line_spec)

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
    per=None,
    kstpkper=None,
    cmap: str = "viridis",
    label: str | None = None,
    title: str | None = None,
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
    """

    show, loc, anchor = legend_placement(legend)
    values = section_values(model, fill, per=per, kstpkper=kstpkper)
    surfaces = section_head_surfaces(
        model, head_layers, per=per, kstpkper=kstpkper
    ) if values is None else []
    names = layer_labels
    if names is None and model is not None:
        names = layer_labels_from_model(model, modelgrid.nlay)
    figure, _axes = plot_layered_cross_section(
        modelgrid,
        _normalize_line(line),
        values=values,
        values_cmap=cmap,
        values_label=label,
        layers=layers,
        layer_labels=names,
        show_grid=show_grid,
        show_legend=show,
        legend_loc=loc,
        legend_anchor=anchor,
        head_surfaces=surfaces,
        flopy_model=getattr(model, "gwf", None),
        title=title,
        show_head=bool(surfaces),
    )
    return figure


def section_head_surfaces(model, head_layers, *, per=None, kstpkper=None) -> list:
    """``[(label, array)]`` -- one water surface per requested layer.

    ``head_layers=None`` draws none, an int or a list draws those layers. The
    default is layer 0, which is what a single water table means and what every
    figure drawn before 2026-09-02 got.

    Labelled per layer rather than "Simulated Head" once, because several
    unlabelled lines on one section are unreadable -- and the label carries the
    layer NAME when the model's build context knows it.
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
    names = layer_labels_from_model(model, nlay) or [f"Layer {k + 1}" for k in range(nlay)]
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
