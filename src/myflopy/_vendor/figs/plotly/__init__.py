"""Trimmed vendored copy: only the core Fig/Subplot/Template surface.

figs' real ``plotly/__init__`` also exposes aq-test/layout/preset/scaling
helpers whose imports require bokeh/cairosvg/reportlab/svglib/svgpathtools;
myflopy does not use them, so the vendored copy imports only ``core``.
"""

from .core import Fig, PlotlyExpressProxy, Subplot, Template, snsfig

__all__ = ["Fig", "PlotlyExpressProxy", "Subplot", "Template", "snsfig"]
