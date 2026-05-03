"""Plotting helpers for :mod:`simple_modflow.modflow.mf6.headsplus`.

The functions here intentionally focus on presentation concerns so the
``HeadsPlus`` reader can remain centered on data access and shaping.
"""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import pandas as pd

import figs

from . import mf2Dplots

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.headsplus import HeadsPlus
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

idxx = pd.IndexSlice
crs_latlon = "EPSG:4326"


def multimodel_plot_heads(models: list["SimulationBase"], locs: int | list[int] | Path, **kwargs):
    """Combine one heads plot from each model into a single figure."""

    figures = []
    for model in models:
        fig = model.hds.plot_heads(locs=locs, plot_fig=False, return_fig=True, **kwargs)
        figures.append(fig)

    model_order = {model.name: index for index, model in enumerate(models)}
    grouped: "OrderedDict[str, list[tuple[int, object]]]" = OrderedDict()

    for model, fig in zip(models, figures):
        for trace in fig.data:
            base_name = trace.name or "trace"
            trace_json = trace.to_plotly_json()
            trace_json["name"] = f"{base_name}-{model.name}"
            trace_json["legendgroup"] = base_name
            grouped.setdefault(base_name, []).append((model_order[model.name], trace_json))

    out = figs.Fig()
    out.update_layout(figures[0].layout)
    for base_name, items in grouped.items():
        for _, trace_json in sorted(items, key=lambda item: item[0]):
            out.add_trace(figs.Fig(data=[trace_json]).data[0])
    return out.show()


def _plot_x_values(
    heads: "HeadsPlus",
    *,
    times: pd.DatetimeIndex | None,
    show_times: bool,
    show_dates: bool,
) -> list | pd.DatetimeIndex:
    """Resolve x-axis values for heads time-series plots."""

    if times is not None:
        if not isinstance(times, pd.DatetimeIndex):
            times = pd.DatetimeIndex(times)
        return times
    if show_times:
        periods = heads.all_heads.index.get_level_values("kstpkper").unique().to_list()
        return [round(heads.model.times[per]) for per in periods]
    if show_dates:
        return heads.model.per_dates if heads.model.per_dates is not None else list(range(len(heads.kstpkper)))
    return list(range(len(heads.kstpkper)))


def plot_heads(
    heads: "HeadsPlus",
    locs: Path | int | list,
    *,
    crs: str | None = None,
    layer: int = 0,
    loc_name_field: str = "ExploName",
    plot_fig: bool = True,
    return_fig: bool = False,
    show_dates: bool = False,
    show_times: bool = False,
    start_period: int = 0,
    loc_names: list | None = None,
    times: pd.DatetimeIndex | None = None,
):
    """Plot heads at specific observation locations or cell ids."""

    crs = heads.crs if crs is None else crs
    if locs is None:
        return None

    fig = figs.Fig()
    data = heads.all_heads

    if start_period > 0:
        per_tuples = data.index.get_level_values("kstpkper")
        mask = [kper >= start_period for (_, kper) in per_tuples]
        data = data[mask]

    if isinstance(locs, Path):
        obs_dict = heads.vor.get_vor_cells_as_dict(
            locs=locs,
            crs=crs,
            predicate="contains",
            loc_name_field=loc_name_field,
        )
        obs_dict = {key: value[0] for key, value in obs_dict.items() if len(value[0]) > 0}
        obs_df = pd.DataFrame.from_dict(obs_dict).transpose()
        obs_locs = obs_df.index
    elif isinstance(locs, int):
        obs_locs = [locs]
        loc_names = loc_names if loc_names is not None else [locs]
        assert len(obs_locs) == len(loc_names), "loc_names must be the same length as locs"
    elif isinstance(locs, list):
        obs_locs = locs
        loc_names = loc_names if loc_names is not None else locs
        assert len(obs_locs) == len(loc_names), "loc_names must be the same length as locs"
    else:
        raise ValueError("locs must be a file path, integer, or list of integers")

    for index, obs_loc in enumerate(obs_locs):
        if isinstance(locs, Path):
            obs_heads = data.loc[idxx[:, layer, obs_df.loc[obs_loc]], "elev"]
        else:
            obs_heads = data.loc[idxx[:, layer, obs_loc], "elev"]

        xs = _plot_x_values(heads, times=times, show_times=show_times, show_dates=show_dates)
        fig.add_scattergl(
            x=xs,
            y=obs_heads,
            name=loc_names[index] if loc_names is not None else obs_loc,
        )

    if plot_fig:
        fig.show()
    if return_fig:
        return fig
    return None


def choropleth(
    heads: "HeadsPlus",
    *,
    kstpkper: tuple = (0, 0),
    plot_mounding: bool = False,
    zmin=None,
    zmax=None,
    zoom=13,
    custom_hover: dict | None = None,
    bottom=None,
    bottom_array=None,
    all_layers: bool = False,
    layer: int = 1,
    obs: Path | None = None,
    obs_name: str = "ExploName",
):
    """Instantiate a choropleth figure for one stress period and layer."""

    stp_to_plot = kstpkper[0]
    per_to_plot = kstpkper[1]
    kstpkper_key = f"sp{per_to_plot}ts{stp_to_plot}"
    choro_dict = {}
    choro_heads = {}
    vor = heads.vor

    if bottom:
        bottom_elev = bottom
    elif plot_mounding:
        if bottom_array is not None:
            bottom_elev = bottom_array
        else:
            bottom_elev = vor.gdf_topbtm.loc[:, layer].to_numpy()
    else:
        bottom_elev = 0

    if all_layers is True:
        pass
    choro_heads[kstpkper_key] = heads.all_heads.loc[idxx[(stp_to_plot, per_to_plot), layer], :]
    choro_dict[kstpkper_key] = choro_heads[kstpkper_key]["elev"] - bottom_elev

    if zmax is None:
        zmax = choro_dict[kstpkper_key].max()
    if zmin is None:
        zmin = choro_dict[kstpkper_key].min()

    fig_mbox = mf2Dplots.ChoroplethPlot(vor=vor, zoom=zoom)

    hover_dict = {
        "Cell No.": heads.cell_list,
        "Area": heads.area_list,
        "x": heads.x_list,
        "y": heads.y_list,
    }
    for lyr in range(heads.nlay):
        lyr_heads = heads.all_heads.loc[idxx[(stp_to_plot, per_to_plot), lyr], "elev"].to_list()
        hover_dict[f"Layer {lyr + 1} Heads"] = lyr_heads
    if plot_mounding:
        hover_dict["Mounding"] = choro_dict[kstpkper_key].to_list()
    if custom_hover:
        for name, data in custom_hover.items():
            hover_dict[str(name)] = data

    custom_data, hover_template = figs.create_hover(hover_dict)
    geojson = getattr(vor, "latslons", None)
    if geojson is None:
        geojson = vor.gdf_latlon.__geo_interface__

    fig_mbox.add_choroplethmap(
        geojson=geojson,
        featureidkey="id",
        locations=vor.gdf_latlon.index.to_list(),
        z=choro_dict[kstpkper_key],
        hovertemplate=hover_template,
        customdata=custom_data,
        colorscale="earth",
        zmax=zmax,
        zmin=zmin,
    )
    if obs:
        obs_gdf = gpd.read_file(obs).to_crs(crs_latlon)
        fig_mbox.add_scattermap(
            lat=obs_gdf.geometry.y,
            lon=obs_gdf.geometry.x,
            text=obs_gdf[obs_name],
            hoverinfo="text",
            marker_color="red",
        )

    return fig_mbox


def plot_choropleth(heads: "HeadsPlus", *args, **kwargs):
    """Plot a choropleth immediately and return ``None``."""

    fig = choropleth(heads, *args, **kwargs)
    fig.show()
