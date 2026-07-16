"""Matplotlib cross-section plotting for MODFLOW model results.

This module is intentionally model-aware. It focuses on report-style and
exploration figures that show discretization, layers, and simulated heads,
while optionally borrowing a matplotlib theme from :mod:`figs`.
"""

from __future__ import annotations
from myflopy.viz import mpl_axes

from collections.abc import Sequence
from contextlib import nullcontext
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import shapely as shp
from flopy.plot import PlotCrossSection
from matplotlib.colors import ListedColormap

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

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
    except Exception:
        return nullcontext()
    return REPORT.context()


def _normalize_line(line) -> dict:
    """Normalize a section line (dict, ``LineString``, or coord sequence) to a ``{"line": coords}`` dict."""

    if isinstance(line, dict):
        return line

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
    model: "SimulationBase",
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
    head_label: str,
    head_color: str,
    head_linewidth: float,
    include_head: bool,
) -> list:
    """Legend handles for the section: a color patch per layer, plus the head line when included."""

    handles = [
        mpatches.Patch(color=color, label=label)
        for color, label in zip(layer_colors, layer_labels)
    ]
    if include_head:
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=head_color,
                linewidth=head_linewidth,
                label=head_label,
            )
        )
    return handles


def plot_layered_cross_section(
    modelgrid,
    line,
    *,
    ax: "Axes" | None = None,
    style: ModelCrossSectionStyle | None = None,
    layer_colors: Sequence[str] | None = None,
    layer_labels: Sequence[str] | None = None,
    head_surface=None,
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
    flopy_model
        Optional flopy model passed through to ``PlotCrossSection`` (only needed
        for model-aware head rendering).

    Returns ``(fig, ax)``.
    """
    style = ModelCrossSectionStyle() if style is None else style
    line_spec = _normalize_line(line)
    nlay = modelgrid.nlay
    layer_id = _resolve_layer_array(modelgrid)
    resolved_layer_colors = _resolve_layer_colors(nlay, layer_colors, style)
    resolved_layer_labels = _resolve_layer_labels(nlay, layer_labels)

    with _theme_context(style):
        if ax is None:
            fig, ax = mpl_axes(figsize=style.figsize)
        else:
            fig = ax.figure

        xsect = PlotCrossSection(model=flopy_model, modelgrid=modelgrid, line=line_spec)

        if show_grid:
            xsect.plot_grid(ax=ax, linewidths=style.grid_linewidth, color=style.grid_color)

        if show_layers:
            xsect.plot_array(
                layer_id,
                ax=ax,
                cmap=ListedColormap(resolved_layer_colors),
                vmin=0,
                vmax=nlay - 1,
                alpha=style.layer_alpha,
            )

        if show_head and head_surface is not None:
            xsect.plot_surface(
                head_surface,
                ax=ax,
                color=style.head_color,
                lw=style.head_linewidth,
            )

        if xlim is not None:
            ax.set_xlim(*xlim)
        if ylim is not None:
            ax.set_ylim(*ylim)

        ax.set_xlabel(xlabel or style.xlabel, size=style.label_fontsize)
        ax.set_ylabel(ylabel or style.ylabel, size=style.label_fontsize)
        ax.tick_params(axis="both", which="major", labelsize=style.tick_labelsize)

        if title or style.title:
            ax.set_title(title or style.title, fontsize=style.title_fontsize)

        if show_legend:
            handles = _build_legend_handles(
                layer_colors=resolved_layer_colors,
                layer_labels=resolved_layer_labels,
                head_label=style.head_label,
                head_color=style.head_color,
                head_linewidth=style.head_linewidth,
                include_head=show_head and head_surface is not None,
            )
            ax.legend(
                handles=handles,
                loc=style.legend_loc,
                frameon=style.legend_frameon,
                fontsize=style.legend_fontsize,
            )

    return fig, ax


def plot_model_cross_section(
    model: "SimulationBase",
    line,
    *,
    kstpkper: tuple | None = None,
    head_data=None,
    head_layer: int = 0,
    ax: "Axes" | None = None,
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
