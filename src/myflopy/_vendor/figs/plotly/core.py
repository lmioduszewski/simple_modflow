"""
Core Plotly figure classes for `figs`.

This module re-exports the established `Fig`, `Subplot`, and `Template`
objects from the legacy implementation module under the preferred
`figs.plotly` namespace.
"""

from .._fig import Fig, PlotlyExpressProxy, Subplot, Template, snsfig

__all__ = ["Fig", "PlotlyExpressProxy", "Subplot", "Template", "snsfig"]
