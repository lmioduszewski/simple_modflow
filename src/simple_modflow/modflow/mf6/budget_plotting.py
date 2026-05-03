"""Plotting helpers for :mod:`simple_modflow.modflow.mf6.budget`."""

from __future__ import annotations

from itertools import cycle
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from figs import Fig
from pandas import IndexSlice as idxx
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors

from simple_modflow.modflow.mf6.budget_tables import budget_obs_df, coerce_plot_times

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.budget import Budget
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase


def plot_budget_obs(
    budget: "Budget",
    *,
    model: "SimulationBase",
    shp_gpkg: Path,
    q: str = "q",
    plot_fig: bool = True,
    return_fig: bool = False,
    name_field: str = "name",
    times=None,
    multiplier: float = 24 * 60 * 60,
    y_range: list[float] = [0, 4],
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
    models: list["SimulationBase"],
    *,
    model_package: str = "drn",
    shp_gpkg: Path,
    q: str = "q",
    name_field: str = "name",
    times: pd.DatetimeIndex | None = None,
    plot_fig: bool = True,
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
    return fig


def plot_drn_choropleth(model: "SimulationBase", *, per: int = 0, zmax=None):
    """Plot DRN flows for one zero-based stress-period index as a choropleth.

    Parameters
    ----------
    model
        Parent model whose DRN budget should be visualized.
    per
        Zero-based stress-period index into ``model.kstpkper``.
    zmax
        Optional explicit upper bound for the colorscale.
    """

    kstpkper = model.kstpkper[per]
    drn_df = model.bud("drn").df
    drn_flows = drn_df.loc[idxx[:, kstpkper], :].q.droplevel(1) * -1
    if zmax is None:
        zmax = drn_flows.max()
    node = drn_flows.reset_index().drop_duplicates("node")
    node.set_index("node", inplace=True)
    full_idx = range(model.vor.ncpl)
    drn_flows = node.reindex(full_idx, fill_value=0)
    model.choro(per=per, custom_zs=drn_flows.q.to_list(), zmin=0, zmax=zmax).plot()
