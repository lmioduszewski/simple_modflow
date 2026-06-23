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
``surfaces``, ``mf3dplots``), **mapbox maps** (the node-id / cell debug plots in
``grid/plotting``; the model-data choropleth already runs through ``Fig`` via the
``Choro`` class), and **animation re-wraps** (``interactive_plotting`` rebuilds a
figure from existing data + frames). The 2-D house template (paper-anchored
border, x/y axis styling) does not belong on those, so they stay raw.
"""

from __future__ import annotations

# Re-export the figs primitives the project uses, so `myflopy.viz` is a superset
# drop-in for `figs`: a module can `from myflopy import viz as f` (or
# `from myflopy.viz import Fig, create_hover`) and never import figs directly.
from figs import Fig, Subplot, Template, create_hover
from figs.mpl import REPORT, Theme, get_mplfig
from plotly.subplots import make_subplots as _make_subplots

__all__ = [
    "Fig",
    "Subplot",
    "Template",
    "create_hover",
    "subplots",
    "mpl_axes",
    "report_axes",
    "Theme",
    "REPORT",
    "PALETTE",
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
    noise = "rgba(214,39,40,0.45)"
    measured = "rgb(214,39,40)"
    truth = "rgb(214,39,40)"
    conflict = "darkorange"

    # Matplotlib equivalents
    mpl_prior = "0.6"
    mpl_posterior = "#1f77b4"
    mpl_measured = "crimson"
    mpl_truth = "crimson"
    mpl_conflict = "darkorange"
