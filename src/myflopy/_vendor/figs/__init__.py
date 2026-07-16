"""Vendored snapshot of figs (see myflopy/_vendor/README.md).

Trimmed: the aq-test/scaling/cross-section stacks are omitted.
"""


from typing import Any

from ._datatypes import ExcelData, ExcelDateData
from .plotly.core import Fig, Subplot, Template

__version__ = "0.1.0"

__all__ = [
    "ExcelData",
    "ExcelDateData",
    "Fig",
    "Subplot",
    "Template",
    "add_to_hover_dict",
    "create_hover",
    "plot",
]


def create_hover(name_dict: dict[str, list[Any]] | None = None):
    """Create Plotly hover metadata and a hovertemplate string."""
    return Template.create_hover(name_dict)


def add_to_hover_dict(existing_hover_dict: dict[str, Any], dict_to_add: dict[str, Any]):
    """Merge additional hover fields into an existing hover-data dictionary."""
    return Template.add_to_hover_dict(existing_hover_dict, dict_to_add)


def plot(data, index_is_x=False):
    """Tries to plot whatever data is passed in."""
    import pandas as pd

    fig = Fig()
    if isinstance(data, pd.Series):
        fig.add_scattergl(
            x=data.index,
            y=data,
            name=data.name,
        )
    elif isinstance(data, pd.DataFrame):
        if index_is_x:
            data = data.reset_index(drop=False)
        for i, col in enumerate(data.columns):
            if i == 0:
                continue
            fig.add_scattergl(
                x=data.iloc[:, 0],
                y=data.iloc[:, i],
                name=col,
            )
    return fig.show()
