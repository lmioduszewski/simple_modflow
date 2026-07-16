from __future__ import annotations

"""
Plotting helpers for source-agnostic Matplotlib cross sections.

This module focuses on rendering and reusable style presets. Data loading and
normalization live in `figs.mpl.cross_section_data`, which keeps the plotting
surface lighter and easier to extend.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import shapely as shp

from .cross_section_data import (
    line_to_profile_frame,
    load_cross_section_data,
    read_section_line,
)
from .mpl import get_mplfig, set_axis_scale

DEFAULT_XLABEL = "Distance"
DEFAULT_YLABEL = "Elevation"


@dataclass(slots=True)
class CrossSectionStyle:
    """
    Reusable style preset for `plot_cross_section(...)`.

    This dataclass bundles the most common presentation and print-scale options
    for cross-section plots so they can be reused across multiple figures.
    Explicit keyword arguments passed to `plot_cross_section(...)` always take
    precedence over values stored on the style object.

    Parameters
    ----------
    title : str, optional
        Default title for the figure.
    xlabel : str, default "Distance"
        Default x-axis label.
    ylabel : str, default "Elevation"
        Default y-axis label.
    show_legend : bool, default True
        Whether grouped plots should show a legend by default.
    x_units_per_inch : float, optional
        Default horizontal engineering scale in data units per inch.
    y_units_per_inch : float, optional
        Default vertical engineering scale in data units per inch.
    figure_kwargs : dict, optional
        Extra keyword arguments forwarded to `get_mplfig(...)`, such as
        `figsize`, `font_size`, or `sns_theme`.
    line_kwargs : dict, optional
        Extra keyword arguments forwarded to `Axes.plot(...)`, such as
        `color`, `linewidth`, or `linestyle`.
    scale_kwargs : dict, optional
        Extra keyword arguments forwarded to `set_axis_scale(...)`, such as
        `x_anchor`, `y_anchor`, `x_value`, `y_value`, or `tick_inches`.

    Notes
    -----
    A common pattern is to define one `CrossSectionStyle` for a report or map
    series and reuse it across many figures, overriding only the title or scale
    on individual calls when needed.
    """

    title: str | None = None
    xlabel: str = DEFAULT_XLABEL
    ylabel: str = DEFAULT_YLABEL
    show_legend: bool = True
    x_units_per_inch: float | None = None
    y_units_per_inch: float | None = None
    figure_kwargs: dict[str, Any] = field(default_factory=dict)
    line_kwargs: dict[str, Any] = field(default_factory=dict)
    scale_kwargs: dict[str, Any] = field(default_factory=dict)

    def merged_figure_kwargs(self, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
        merged = dict(self.figure_kwargs)
        if overrides:
            merged.update(overrides)
        return merged

    def merged_line_kwargs(self, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
        merged = dict(self.line_kwargs)
        if overrides:
            merged.update(overrides)
        return merged

    def merged_scale_kwargs(self, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
        merged = dict(self.scale_kwargs)
        if overrides:
            merged.update(overrides)
        return merged


def plot_cross_section(
    data: pd.DataFrame | Path | str | shp.LineString | shp.MultiLineString | gpd.GeoDataFrame,
    x: str = "x",
    y: str = "y",
    series_col: str | None = None,
    style: CrossSectionStyle | None = None,
    show_legend: bool | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    x_units_per_inch: float | None = None,
    y_units_per_inch: float | None = None,
    layer: str | int | None = None,
    simplify_tolerance: float | None = None,
    sheet_name: str | int = 0,
    read_kwargs: dict[str, Any] | None = None,
    scale_kwargs: dict | None = None,
    line_kwargs: dict | None = None,
    **figure_kwargs,
):
    """
    Plot a source-agnostic cross section and optionally apply an exact print scale.

    Parameters
    ----------
    data : DataFrame, Path, str, LineString, MultiLineString, or GeoDataFrame
        Cross-section source. This can be:

        - a `pandas.DataFrame`
        - a GIS path such as `.gpkg`, `.shp`, or `.geojson`
        - an Excel path such as `.xlsx`
        - a CSV or TSV path
        - a Shapely `LineString` or `MultiLineString`
        - a `geopandas.GeoDataFrame`

        Non-dataframe sources are normalized internally by
        `load_cross_section_data(...)`.
    x : str, default "x"
        Column name to use for horizontal values. For tabular inputs, this
        defaults to `x`. If `x` and `y` are left at their defaults and the
        input does not contain literal `x`/`y` columns, the first two columns
        are assumed to represent x and y.
    y : str, default "y"
        Column name to use for vertical values. For geometry-derived sources,
        this will usually be `y` or `elevation`. For tabular data, the same
        first-two-columns fallback applies when the default names are used.
    series_col : str, optional
        Optional grouping column. When provided, one line is drawn for each
        distinct value in `series_col`.
    style : CrossSectionStyle, optional
        Reusable style preset. Explicit keyword arguments passed directly to
        this function override the style values.
    show_legend : bool, optional
        Whether grouped plots should show a legend. If omitted, the value from
        `style.show_legend` is used.
    title : str, optional
        Figure title. Overrides `style.title`.
    xlabel : str, optional
        X-axis label. Overrides `style.xlabel`.
    ylabel : str, optional
        Y-axis label. Overrides `style.ylabel`.
    x_units_per_inch : float, optional
        Horizontal engineering scale in data units per inch. When both x and y
        scale values are provided, `set_axis_scale(...)` is applied after
        plotting.
    y_units_per_inch : float, optional
        Vertical engineering scale in data units per inch.
    layer : str or int, optional
        GIS layer name or index used when `data` points to a GIS dataset.
    simplify_tolerance : float, optional
        Optional simplification tolerance applied to GIS-derived geometry before
        plotting.
    sheet_name : str or int, default 0
        Excel sheet name or index used when `data` points to an Excel file.
    read_kwargs : dict, optional
        Additional keyword arguments passed to pandas when reading tabular file
        sources such as Excel, CSV, or TSV.
    scale_kwargs : dict, optional
        Extra keyword arguments passed to `set_axis_scale(...)`.
    line_kwargs : dict, optional
        Extra keyword arguments passed to Matplotlib `Axes.plot(...)`.
    **figure_kwargs
        Additional keyword arguments forwarded to `get_mplfig(...)`.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        The created figure and axes.

    Notes
    -----
    This function is intentionally lightweight and source-agnostic. It is meant
    for exact-scale profile plotting and general cross-section presentation, not
    for MODFLOW-specific results exploration. Model-aware cross sections that
    show discretization, arrays, or simulated heads should stay in
    `simple_modflow`.

    Examples
    --------
    Plot a dataframe that already has `x` and `y` columns:

    ```python
    fig, ax = plot_cross_section(profile_df)
    ```

    Plot a dataframe with explicit column names:

    ```python
    fig, ax = plot_cross_section(profile_df, x="distance", y="elevation")
    ```

    Plot a GeoPackage line directly:

    ```python
    fig, ax = plot_cross_section(
        "section.gpkg",
        layer=1,
        simplify_tolerance=0.5,
        x_units_per_inch=2000,
        y_units_per_inch=100,
    )
    ```

    Plot an Excel file directly, using the first two columns as x and y:

    ```python
    fig, ax = plot_cross_section("section_profile.xlsx")
    ```

    Reuse a style preset:

    ```python
    style = CrossSectionStyle(
        title="Lake Sawyer",
        x_units_per_inch=2000,
        y_units_per_inch=100,
        scale_kwargs={"y_anchor": "max", "y_value": 750},
    )
    fig, ax = plot_cross_section("section.gpkg", style=style, layer=1)
    ```
    """
    style = style or CrossSectionStyle()
    data = load_cross_section_data(
        data,
        x=x,
        y=y,
        layer=layer,
        simplify_tolerance=simplify_tolerance,
        sheet_name=sheet_name,
        read_kwargs=read_kwargs,
    )
    resolved_title = title if title is not None else style.title
    resolved_xlabel = xlabel if xlabel is not None else style.xlabel
    resolved_ylabel = ylabel if ylabel is not None else style.ylabel
    resolved_show_legend = show_legend if show_legend is not None else style.show_legend
    resolved_x_units_per_inch = (
        x_units_per_inch if x_units_per_inch is not None else style.x_units_per_inch
    )
    resolved_y_units_per_inch = (
        y_units_per_inch if y_units_per_inch is not None else style.y_units_per_inch
    )

    fig, ax = get_mplfig(
        title=resolved_title,
        xlabel=resolved_xlabel,
        ylabel=resolved_ylabel,
        **style.merged_figure_kwargs(figure_kwargs),
    )
    plot_kwargs = style.merged_line_kwargs(line_kwargs)

    if series_col is None:
        ax.plot(data[x], data[y], **plot_kwargs)
    else:
        for series_name, series_df in data.groupby(series_col, sort=False):
            ax.plot(
                series_df[x],
                series_df[y],
                label=series_name,
                **plot_kwargs,
            )
        if resolved_show_legend:
            ax.legend(title="")

    if resolved_x_units_per_inch is not None and resolved_y_units_per_inch is not None:
        set_axis_scale(
            ax,
            x_units_per_inch=resolved_x_units_per_inch,
            y_units_per_inch=resolved_y_units_per_inch,
            **style.merged_scale_kwargs(scale_kwargs),
        )

    return fig, ax


__all__ = [
    "CrossSectionStyle",
    "line_to_profile_frame",
    "load_cross_section_data",
    "plot_cross_section",
    "read_section_line",
]
