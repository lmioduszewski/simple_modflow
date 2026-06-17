"""Matplotlib cross-section plotting for MODFLOW model results.

This module is intentionally model-aware. It focuses on report-style and
exploration figures that show discretization, layers, and simulated heads,
while optionally borrowing a matplotlib theme from :mod:`figs`.
"""

from __future__ import annotations

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
    if not style.use_figs_theme:
        return nullcontext()
    try:
        from figs.mpl import REPORT
    except Exception:
        return nullcontext()
    return REPORT.context()


def _normalize_line(line) -> dict:
    if isinstance(line, dict):
        return line

    if isinstance(line, shp.LineString):
        return {"line": list(line.coords)}

    if isinstance(line, Sequence) and not isinstance(line, (str, bytes)):
        return {"line": list(line)}

    raise TypeError(f"Unsupported line specification: {type(line)}")


def _resolve_layer_array(model) -> np.ndarray:
    grid = model.gwf.modelgrid
    nlay = grid.nlay
    if hasattr(grid, "ncpl"):
        ncpl = grid.ncpl
        return np.arange(nlay)[:, None] * np.ones((1, ncpl), dtype=int)
    raise ValueError("plot_model_cross_section currently expects a grid with an 'ncpl' attribute.")


def _resolve_head_data(
    model: "SimulationBase",
    *,
    kstpkper: tuple | None,
    head_data,
):
    if head_data is not None:
        return np.asarray(head_data)
    if kstpkper is None:
        return None
    return np.asarray(model.gwf.output.head().get_data(kstpkper=kstpkper))


def _resolve_head_surface(head_data, head_layer: int):
    if head_data is None:
        return None
    arr = np.asarray(head_data)
    if arr.ndim == 1:
        return arr
    if arr.ndim >= 2:
        return arr[head_layer]
    return arr


def _resolve_layer_colors(nlay: int, layer_colors: Sequence[str] | None, style: ModelCrossSectionStyle) -> list[str]:
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
    style = ModelCrossSectionStyle() if style is None else style
    line_spec = _normalize_line(line)
    head_array = _resolve_head_data(model, kstpkper=kstpkper, head_data=head_data)
    head_surface = _resolve_head_surface(head_array, head_layer=head_layer)
    layer_id = _resolve_layer_array(model)
    nlay = model.gwf.modelgrid.nlay
    resolved_layer_colors = _resolve_layer_colors(nlay, layer_colors, style)
    resolved_layer_labels = _resolve_layer_labels(nlay, layer_labels)

    with _theme_context(style):
        if ax is None:
            fig, ax = plt.subplots(figsize=style.figsize)
        else:
            fig = ax.figure

        xsect = PlotCrossSection(
            model=model.gwf,
            modelgrid=model.gwf.modelgrid,
            line=line_spec,
        )

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
