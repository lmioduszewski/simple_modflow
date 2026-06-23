"""Plot payload builders and plotting mixins for package explorers."""

from __future__ import annotations
from myflopy.viz import mpl_axes

from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from flopy.plot import PlotMapView
from plotly.subplots import make_subplots

if TYPE_CHECKING:
    pass
from myflopy.modflow.mf6.package_explorer_utils import (
    _aggregate_hover_strings,
)


def _symmetric_color_limit(values: Iterable[float]) -> float:
    """Return the symmetric absolute color bound for a numeric sequence."""

    array = np.asarray(list(values), dtype=float)
    if array.size == 0:
        return 0.0
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return 0.0
    return float(np.max(np.abs(finite)))


def _blue_white_red_diverging_colorscale() -> list[list[object]]:
    """Return a blue-white-red diverging colorscale for signed maps."""

    return [
        [0.0, "#1f77b4"],
        [0.5, "#ffffff"],
        [1.0, "#d62728"],
    ]


def _as_layer_cell_property(values, *, nlay: int, ncpl: int, label: str) -> np.ndarray:
    """Return static package arrays as a consistent ``(nlay, ncpl)`` array."""

    arr = np.asarray(values)
    arr = np.squeeze(arr)
    if arr.ndim == 0:
        return np.full((int(nlay), int(ncpl)), arr.item())
    if arr.ndim == 1:
        if arr.size == int(ncpl):
            return np.tile(arr.reshape(1, int(ncpl)), (int(nlay), 1))
        if arr.size == int(nlay) * int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
    if arr.ndim == 2:
        if arr.shape == (int(nlay), int(ncpl)):
            return arr
        if arr.shape == (int(ncpl), int(nlay)):
            return arr.T
        if arr.size == int(nlay) * int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
    if arr.ndim >= 3:
        if arr.shape[0] == int(nlay) and int(np.prod(arr.shape[1:])) == int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
        if arr.size == int(nlay) * int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
    raise ValueError(
        f"{label} must resolve to one value per layer/cell; got shape={arr.shape}."
    )


def build_sfr_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
) -> tuple[list[float], dict[str, list]]:
    """Convert SFR exchange rows into a cell choropleth normalized by reach length.

    The mapped value is ``sum(q) / sum(rlen)`` within each cell, which is more
    meaningful than raw exchange magnitude alone when reach lengths vary across
    cells.
    """

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name="q_per_length")
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            "q_per_length": values.tolist(),
            "q": [0.0 for _ in full_index],
            "rlen": [0.0 for _ in full_index],
        }
        return values.tolist(), hover

    selected = frame.copy()
    selected["q"] = pd.to_numeric(selected["q"], errors="coerce")
    selected["rlen"] = pd.to_numeric(selected["rlen"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    q_sum = grouped["q"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    rlen_sum = grouped["rlen"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    values = pd.Series(float(fill_value), index=full_index, dtype=float)
    valid = rlen_sum > 0.0
    values.loc[valid] = (q_sum.loc[valid] / rlen_sum.loc[valid]).astype(float)
    values = values * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        "q_per_length": values.tolist(),
        "q": q_sum.fillna(0.0).astype(float).tolist(),
        "rlen": rlen_sum.fillna(0.0).astype(float).tolist(),
    }
    if "reach" in selected.columns:
        reach_strings = selected.assign(
            reach=selected["reach"].astype("string")
        ).groupby("cell", dropna=False)["reach"]
        hover["reach"] = (
            reach_strings.agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    return values.tolist(), hover


def build_lak_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
) -> tuple[list[float], dict[str, list]]:
    """Convert LAK exchange rows into a cell choropleth normalized by area.

    The mapped value is ``sum(q) / sum(flow_area)`` within each cell, giving a
    lake-groundwater exchange intensity with units of length per time.
    """

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name="q_per_area")
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            "q_per_area": values.tolist(),
            "q": [0.0 for _ in full_index],
            "flow_area": [0.0 for _ in full_index],
        }
        return values.tolist(), hover

    selected = frame.copy()
    selected["q"] = pd.to_numeric(selected["q"], errors="coerce")
    selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    q_sum = grouped["q"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    area_sum = (
        grouped["flow_area"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    )
    values = pd.Series(float(fill_value), index=full_index, dtype=float)
    valid = area_sum > 0.0
    values.loc[valid] = (q_sum.loc[valid] / area_sum.loc[valid]).astype(float)
    values = values * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        "q_per_area": values.tolist(),
        "q": q_sum.fillna(0.0).astype(float).tolist(),
        "flow_area": area_sum.fillna(0.0).astype(float).tolist(),
    }
    if "lake" in selected.columns:
        lake_strings = selected.assign(lake=selected["lake"].astype("string")).groupby(
            "cell", dropna=False
        )["lake"]
        hover["lake"] = (
            lake_strings.agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    if "claktype" in selected.columns:
        hover["claktype"] = (
            grouped["claktype"]
            .agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    return values.tolist(), hover


def build_surface_water_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
) -> tuple[list[float], dict[str, list]]:
    """Convert combined SFR/LAK exchange rows into one shared cell map payload."""

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(
            float(fill_value), index=full_index, name="exchange_intensity"
        )
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "surface_water_exchange": values.tolist(),
            "sfr_exchange": [0.0 for _ in full_index],
            "lak_exchange": [0.0 for _ in full_index],
        }
        return values.tolist(), hover

    selected = frame.copy()
    selected["exchange_intensity"] = pd.to_numeric(
        selected["exchange_intensity"], errors="coerce"
    )
    grouped = selected.groupby("cell", dropna=False)
    values = (
        grouped["exchange_intensity"]
        .sum(min_count=1)
        .reindex(full_index, fill_value=np.nan)
    )
    values = values.fillna(float(fill_value)).astype(float) * float(multiplier)

    sfr_component = selected.loc[selected["source"] == "sfr"].groupby(
        "cell", dropna=False
    )["exchange_intensity"].sum(min_count=1).reindex(full_index, fill_value=0.0).astype(
        float
    ) * float(multiplier)
    lak_component = selected.loc[selected["source"] == "lak"].groupby(
        "cell", dropna=False
    )["exchange_intensity"].sum(min_count=1).reindex(full_index, fill_value=0.0).astype(
        float
    ) * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        "surface_water_exchange": values.tolist(),
        "sfr_exchange": sfr_component.tolist(),
        "lak_exchange": lak_component.tolist(),
    }
    if "source" in selected.columns:
        hover["source"] = (
            grouped["source"]
            .agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    return values.tolist(), hover


def build_cell_input_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    value_column: str,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    agg: str = "sum",
) -> tuple[list[float], dict[str, list]]:
    """Convert a normalized cell-input table into choropleth values and hover.

    Parameters
    ----------
    frame
        Normalized table from :func:`build_cell_package_input_table` or a
        compatible filtered table.
    ncpl
        Number of cells per layer in the target grid.
    value_column
        Numeric column to render as choropleth values.
    per, layer
        Metadata added to the hover table.
    multiplier
        Optional scalar multiplier applied to the mapped values.
    fill_value
        Value used for cells not present in the input data for the selected
        period/layer.
    agg
        Aggregation name passed to ``groupby().agg(...)`` for duplicate cells.
    """

    if value_column not in frame.columns:
        raise KeyError(f"Value column {value_column!r} was not found.")

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name=value_column)
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            value_column: values.astype(float).tolist(),
        }
        return values.astype(float).tolist(), hover

    selected = frame.copy()
    grouped = selected.groupby("cell", dropna=False)
    values = grouped[value_column].agg(agg).reindex(
        full_index, fill_value=fill_value
    ).astype(float) * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        value_column: values.tolist(),
    }

    excluded = {"model", "package", "per", "layer", "cell", value_column}
    for column in selected.columns:
        if column in excluded:
            continue
        column_group = grouped[column]
        if pd.api.types.is_numeric_dtype(selected[column]):
            hover[column] = (
                column_group.agg("first")
                .reindex(full_index, fill_value=np.nan)
                .tolist()
            )
        else:
            hover[column] = (
                column_group.agg(_aggregate_hover_strings)
                .reindex(full_index, fill_value="")
                .tolist()
            )
    return values.tolist(), hover


def build_group_input_compare_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    value_column: str,
    diff_column: str,
    model_name: str,
    reference_model: str,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    agg: str = "sum",
) -> tuple[list[float], dict[str, list], float]:
    """Convert grouped comparison rows into diff-map values and hover data.

    Parameters
    ----------
    frame
        Filtered comparison table produced by a grouped ``compare()`` helper.
    ncpl
        Number of cells per layer in the reference grid.
    value_column
        Raw model value column, such as ``"recharge"`` or ``"finf"``.
    diff_column
        Difference column to map, such as ``"recharge_diff"``.
    model_name, reference_model
        Added to hover metadata.
    per, layer
        Added to hover metadata.
    multiplier
        Optional multiplier applied to the mapped difference values.
    fill_value
        Value used where a cell has no comparison row.
    agg
        Aggregation used when more than one row maps to the same cell.

    Returns
    -------
    tuple[list[float], dict[str, list], float]
        Choropleth z-values, hover payload, and the symmetric absolute color
        range suggested for diverging diff maps.
    """

    if diff_column not in frame.columns:
        raise KeyError(f"Difference column {diff_column!r} was not found.")
    if value_column not in frame.columns:
        raise KeyError(f"Value column {value_column!r} was not found.")

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name=diff_column)
        hover = {
            "Model": [model_name for _ in full_index],
            "Reference Model": [reference_model for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            diff_column: values.astype(float).tolist(),
        }
        return values.astype(float).tolist(), hover, float(abs(fill_value))

    grouped = frame.groupby("cell", dropna=False)
    diff_values = grouped[diff_column].agg(agg).reindex(
        full_index, fill_value=fill_value
    ).astype(float) * float(multiplier)
    current_values = (
        grouped[value_column].agg("first").reindex(full_index, fill_value=np.nan)
    )
    reference_column = f"reference_{value_column}"
    reference_values = (
        grouped[reference_column].agg("first").reindex(full_index, fill_value=np.nan)
        if reference_column in frame.columns
        else pd.Series(np.nan, index=full_index)
    )

    hover: dict[str, list] = {
        "Model": [model_name for _ in full_index],
        "Reference Model": [reference_model for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        value_column: current_values.tolist(),
        reference_column: reference_values.tolist(),
        diff_column: diff_values.tolist(),
    }

    excluded = {
        "model",
        "reference_model",
        "package",
        "per",
        "layer",
        "cell",
        value_column,
        reference_column,
        diff_column,
    }
    for column in frame.columns:
        if column in excluded:
            continue
        column_group = grouped[column]
        if pd.api.types.is_numeric_dtype(frame[column]):
            hover[column] = (
                column_group.agg("first")
                .reindex(full_index, fill_value=np.nan)
                .tolist()
            )
        else:
            hover[column] = (
                column_group.agg(_aggregate_hover_strings)
                .reindex(full_index, fill_value="")
                .tolist()
            )

    absmax = (
        float(np.nanmax(np.abs(diff_values.to_numpy(dtype=float))))
        if len(diff_values)
        else 0.0
    )
    return diff_values.tolist(), hover, absmax


class MappedFieldVisualizationMixin:
    """Shared spatial views for cell-mapped package input and result fields."""

    @property
    def _mapped_value_name(self) -> str:
        for attribute in ("field_name", "value_name"):
            value = getattr(self, attribute, None)
            if value is not None:
                return str(value)
        return "stage"

    def _mapped_periods(self) -> list[int]:
        frame = self.get()
        if "per" not in frame.columns or frame.empty:
            return [0]
        return sorted({int(value) for value in frame["per"].dropna().tolist()})

    def _mapped_layers(self, layers=None) -> list[int]:
        if layers is None:
            return list(range(int(self.model.gwf.modelgrid.nlay)))
        if isinstance(layers, (int, np.integer)):
            return [int(layers)]
        resolved = [int(layer) for layer in layers]
        if not resolved:
            raise ValueError("layers must contain at least one layer.")
        return resolved

    def _mapped_values(self, *, per: int, layer: int, **map_kwargs) -> np.ndarray:
        return np.asarray(self.map(per=per, layer=layer, **map_kwargs).zs, dtype=float)

    def plot(
        self,
        *,
        per: int = 0,
        layers=None,
        ncols: int = 3,
        figsize: tuple[float, float] | None = None,
        cmap: str = "viridis",
        vmin: float | None = None,
        vmax: float | None = None,
        show_grid: bool = True,
        show_colorbar: bool = True,
        title: str | None = None,
        **map_kwargs,
    ):
        """Plot one layer or a selected-layer Matplotlib mosaic."""

        resolved_layers = self._mapped_layers(layers)
        ncols = min(int(ncols), len(resolved_layers))
        nrows = int(np.ceil(len(resolved_layers) / ncols))
        fig, axes = mpl_axes(
            nrows,
            ncols,
            figsize=figsize or (5.0 * ncols, 4.0 * nrows),
            squeeze=False,
        )
        arrays = [
            self._mapped_values(per=per, layer=layer, **map_kwargs)
            for layer in resolved_layers
        ]
        finite_parts = [
            values[np.isfinite(values)]
            for values in arrays
            if np.isfinite(values).any()
        ]
        finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
        if finite.size:
            vmin = float(np.nanmin(finite)) if vmin is None else vmin
            vmax = float(np.nanmax(finite)) if vmax is None else vmax
        image = None
        for ax, layer, values in zip(axes.flat, resolved_layers, arrays, strict=False):
            view = PlotMapView(
                model=self.model.gwf,
                modelgrid=self.model.gwf.modelgrid,
                layer=layer,
                ax=ax,
            )
            image = view.plot_array(values, cmap=cmap, vmin=vmin, vmax=vmax)
            if show_grid:
                view.plot_grid(color="#3c4652", linewidth=0.2)
            ax.set_title(f"Layer {layer + 1}")
            ax.set_aspect("equal")
        for ax in axes.flat[len(resolved_layers) :]:
            ax.set_visible(False)
        if show_colorbar and image is not None:
            fig.colorbar(
                image,
                ax=list(axes.flat[: len(resolved_layers)]),
                shrink=0.75,
                label=self._mapped_value_name,
            )
        fig.suptitle(title or f"{self._mapped_value_name} | stress period {per}")
        return fig

    def plotly_mosaic(
        self,
        *,
        per: int = 0,
        layers=None,
        ncols: int = 3,
        title: str | None = None,
        **map_kwargs,
    ):
        """Return a selected-layer Plotly choropleth mosaic."""

        resolved_layers = self._mapped_layers(layers)
        ncols = min(int(ncols), len(resolved_layers))
        nrows = int(np.ceil(len(resolved_layers) / ncols))
        fig = make_subplots(
            rows=nrows,
            cols=ncols,
            specs=[[{"type": "map"} for _ in range(ncols)] for _ in range(nrows)],
            subplot_titles=[f"Layer {layer + 1}" for layer in resolved_layers],
        )
        traces = []
        for index, layer in enumerate(resolved_layers):
            trace = self.map(per=per, layer=layer, **map_kwargs).get_choropleth()
            trace.coloraxis = "coloraxis"
            traces.append(trace)
            fig.add_trace(trace, row=(index // ncols) + 1, col=(index % ncols) + 1)
        finite_parts = [
            np.asarray(trace.z, dtype=float)[
                np.isfinite(np.asarray(trace.z, dtype=float))
            ]
            for trace in traces
            if np.isfinite(np.asarray(trace.z, dtype=float)).any()
        ]
        finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
        coloraxis = {"colorscale": traces[0].colorscale if traces else "Viridis"}
        if finite.size:
            coloraxis.update(
                cmin=float(np.nanmin(finite)),
                cmax=float(np.nanmax(finite)),
                cauto=False,
            )
        fig.update_layout(
            title=title or f"{self._mapped_value_name} | stress period {per}",
            coloraxis=coloraxis,
            uirevision="lock",
        )
        return fig

    def slider_html(
        self,
        output_path: str | Path,
        *,
        periods=None,
        layers=None,
        ncols: int = 3,
        dpi: int = 160,
        title: str | None = None,
        **map_kwargs,
    ):
        """Export selected layers through stress periods as standalone Matplotlib HTML."""

        from myflopy.modflow.mf6.interactive_plotting import (
            export_matplotlib_slider_html,
        )

        periods = (
            self._mapped_periods()
            if periods is None
            else [int(period) for period in periods]
        )
        labels = [f"Stress period {period}" for period in periods]

        def render(period, index):
            return self.plot(
                per=period,
                layers=layers,
                ncols=ncols,
                title=labels[index],
                **map_kwargs,
            )

        return export_matplotlib_slider_html(
            render,
            periods,
            output_path,
            labels=labels,
            title=title or self._mapped_value_name,
            dpi=dpi,
        )

    def plotly_animation(
        self,
        *,
        periods=None,
        layers=None,
        ncols: int = 3,
        output_path: str | Path | None = None,
        title: str | None = None,
        **map_kwargs,
    ):
        """Return or export a Plotly stress-period animation for selected layers."""

        from myflopy.modflow.mf6.interactive_plotting import _plotly_config
        from myflopy.modflow.utils.animations import Animation
        import plotly.io as pio

        periods = (
            self._mapped_periods()
            if periods is None
            else [int(period) for period in periods]
        )
        figures = [
            self.plotly_mosaic(
                per=period, layers=layers, ncols=ncols, title=title, **map_kwargs
            )
            for period in periods
        ]
        fig = go.Figure(data=figures[0].data, layout=figures[0].layout)
        fig.frames = [
            go.Frame(data=frame.data, name=str(period))
            for period, frame in zip(periods, figures, strict=False)
        ]
        animation = Animation(self.model, periods=periods, redraw=True)
        fig.update_layout(
            updatemenus=animation.updatemenus,
            sliders=animation.sliders,
            uirevision="lock",
        )
        if output_path is not None:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            pio.write_html(
                fig,
                file=output_path,
                include_plotlyjs=True,
                config=_plotly_config(fig),
                auto_play=False,
                auto_open=False,
            )
        return fig


__all__ = [
    "_symmetric_color_limit",
    "_blue_white_red_diverging_colorscale",
    "_as_layer_cell_property",
    "build_sfr_q_map_payload",
    "build_lak_q_map_payload",
    "build_surface_water_q_map_payload",
    "build_cell_input_map_payload",
    "build_group_input_compare_map_payload",
    "MappedFieldVisualizationMixin",
]
