"""Plotting helpers for :mod:`myflopy.modflow.mf6.budget`."""

from __future__ import annotations

from itertools import cycle
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors

from myflopy.modflow.mf6.budget_tables import budget_obs_df, coerce_plot_times
from myflopy.viz import Fig

if TYPE_CHECKING:
    from myflopy.modflow.mf6.budget import Budget
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def plot_budget_obs(
    budget: Budget,
    *,
    model: SimulationBase,
    shp_gpkg: Path,
    q: str = "q",
    plot_fig: bool = True,
    return_fig: bool = False,
    name_field: str = "name",
    times=None,
    multiplier: float = 24 * 60 * 60,
    y_range: list[float] | None = None,
):
    """Plot or return observation-area budget totals for one package.

    Parameters
    ----------
    budget
        Parent budget accessor providing the selected groundwater package.
    model
        Parent model supplying grid and stress-period metadata.
    shp_gpkg
        Polygon dataset defining the observation areas.
    q
        Budget value column to aggregate.
    plot_fig
        If ``True``, display the figure before returning.
    return_fig
        If ``True``, return the constructed figure object.
    name_field
        Attribute field used to label each observation polygon.
    times
        Optional datetime labels for the x-axis.
    multiplier
        Unit-conversion divisor applied after summing flows.
    y_range
        Y-axis range used when plotting.
    """

    times = coerce_plot_times(times)
    frame = budget_obs_df(
        budget,
        model=model,
        shp_gpkg=shp_gpkg,
        q=q,
        name_field=name_field,
        multiplier=multiplier,
    )

    if y_range is None:
        y_range = [0, 4]
    if plot_fig or return_fig:
        fig = Fig()
        for col in frame.columns:
            x = frame.index if times is None else times
            fig.add_scattergl(x=x, y=frame[col], name=col)
        fig.update_yaxes(range=y_range)
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
    return frame


def multimodel_plot_budget_obs(
    models: list[SimulationBase],
    *,
    model_package: str = "drn",
    shp_gpkg: Path,
    q: str = "q",
    name_field: str = "name",
    times: pd.DatetimeIndex | None = None,
    plot_fig: bool = True,
    return_frame: bool = False,
):
    """Plot observation-area budget totals for several models at once.

    Parameters
    ----------
    models
        Models whose observation-area budget totals should be compared.
    model_package
        Groundwater package to summarize, such as ``"drn"`` or ``"rch"``.
    shp_gpkg
        Polygon dataset defining the observation areas.
    q
        Budget value column to aggregate.
    name_field
        Attribute field used to label each observation polygon.
    times
        Optional datetime labels for the x-axis.
    plot_fig
        If ``True``, display the figure before returning it.
    return_frame
        If ``True``, return the dictionary of dataframe data, not the figure
    """

    times = coerce_plot_times(times)
    obs_dict = {}
    color_cycle = cycle(colors)
    fig = Fig()

    for model in models:
        obs_dict[model.name] = model.bud(model_package).plot_budget_obs(
            shp_gpkg=shp_gpkg,
            q=q,
            name_field=name_field,
            plot_fig=False,
        )

    for model_name, frame in obs_dict.items():
        color = next(color_cycle)
        x = frame.index if times is None else times
        for col in frame.columns:
            fig.add_scattergl(
                x=x,
                y=frame[col],
                name=f"{col}-{model_name}",
                line=dict(color=color),
            )

    if plot_fig:
        fig.show()
    if return_frame:
        return obs_dict
    return fig
