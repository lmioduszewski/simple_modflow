"""Plot payload builders and plotting mixins for package explorers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from myflopy.viz import Fig, mpl_axes, shared_map_view, subplots
from myflopy.viz import mosaic as viz_mosaic

if TYPE_CHECKING:
    pass
from myflopy.modflow.mf6.package_explorer_utils import (
    _aggregate_hover_strings,
)


def _apply_backend(choro, backend: str = "plotly"):
    """Return a ``Choro`` as an interactive Plotly figure or a static mpl one.

    Lets a single-model atomic ``map(..., backend=)`` honor the same
    ``backend="plotly"|"mpl"`` switch the :class:`SpatialView` grammar uses,
    without the map having to inherit the engine.
    """

    if str(backend).lower() in ("mpl", "matplotlib", "static"):
        return choro.plot_mpl()
    if str(backend).lower() in ("plotly", "interactive"):
        return choro
    raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")


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
    """Return a blue-white-red diverging colorscale for signed maps.

    Plotly maps position 0.0 to ``zmin`` and 1.0 to ``zmax``, so this is
    blue at the NEGATIVE end and red at the POSITIVE end. Use
    :func:`_exchange_colorscale` for surface-water exchange rather than
    reaching for this directly -- which sign means "gaining" is a per-package
    fact, not a constant.
    """

    return [
        [0.0, "#1f77b4"],
        [0.5, "#ffffff"],
        [1.0, "#d62728"],
    ]


def _red_white_blue_diverging_colorscale() -> list[list[object]]:
    """Return the diff-map orientation as STOPS: red at ``zmin``, blue at ``zmax``.

    The same endpoint colors as :func:`_blue_white_red_diverging_colorscale`,
    reversed -- the "negative red, positive blue" rule that
    ``get_default_group_compare_colorscale() == 'RdBu'`` declares for difference
    maps. It exists as stops rather than as that NAME because
    ``_PLOTLY_TO_MPL_CMAP['rdbu']`` is ``'RdBu_r'`` (choros.py:36), the REVERSED
    colormap, so a named diverging scale renders mirrored between
    ``Choro.plot()`` and ``Choro.plot_mpl()``. Stops survive both backends: the
    colorscale setter passes lists through, and ``plot_mpl`` rebuilds them with
    ``LinearSegmentedColormap.from_list``.

    Derived from its sibling rather than restating the hexes, so the two
    orientations cannot drift apart.
    """

    colors = [color for _, color in reversed(_blue_white_red_diverging_colorscale())]
    return [[position, color] for position, color in zip((0.0, 0.5, 1.0), colors, strict=True)]


def _exchange_colorscale(frame: str) -> list[list[object]]:
    """Return the exchange scale oriented for ``frame``: gaining blue, losing red.

    ``frame`` is a ``ResultSpec.reference_frame`` value and MUST be supplied --
    which sign means "gaining" is a property of the source file, not a constant,
    because myflopy keeps MF6's raw sign and never normalizes it:

    * ``"gwf"`` -- an aquifer-referenced cell (.cbc) record: MF6 writes flow FROM
      the feature TO the cell, so a GAINING feature is NEGATIVE. Blue sits at the
      negative end -> ``[blue, white, red]``. (SFR and the list BCs.)
    * ``"feature"`` -- a feature-referenced package budget file (LAK's ``GWF``
      record) OR myflopy's normalized ``exchange_intensity``: a gaining feature
      is POSITIVE. Blue sits at the positive end -> ``[red, white, blue]``.

    Plotly maps position 0.0 to ``zmin`` and 1.0 to ``zmax``. Passing a stale
    frame is what inverted the SFR map on 2026-07-19; callers now derive it from
    the registry, which describes what MF6 wrote and so cannot drift from the
    data.
    """

    blue, white, red = "#1f77b4", "#ffffff", "#d62728"
    if frame == "gwf":
        return [[0.0, blue], [0.5, white], [1.0, red]]
    return [[0.0, red], [0.5, white], [1.0, blue]]

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


#: per-feature result columns passed through to LAK/SFR q-map hovers when a
#: caller has joined them onto the exchange frame (see ``join_*_stage``).
_FEATURE_HOVER_COLUMNS = ("stage", "depth")


def _append_feature_hover_columns(hover, selected, full_index) -> None:
    """Add joined per-feature columns (stage, depth) to a q-map hover payload.

    Stage repeats on every cell a feature touches, so a per-cell mean is the
    faithful aggregate (matches ``StageResultsExplorer``'s series default).
    """

    for column in _FEATURE_HOVER_COLUMNS:
        if column not in selected.columns or column in hover:
            continue
        numeric = pd.to_numeric(selected[column], errors="coerce")
        hover[column] = (
            numeric.groupby(selected["cell"], dropna=False)
            .mean()
            .reindex(full_index, fill_value=np.nan)
            .tolist()
        )


def build_sfr_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    value_column: str = "q_gwf",
) -> tuple[list[float], dict[str, list]]:
    """Convert SFR exchange rows into a cell choropleth normalized by reach length.

    The mapped value is ``sum(q) / sum(rlen)`` within each cell, which is more
    meaningful than raw exchange magnitude alone when reach lengths vary across
    cells. ``value_column`` is the frame's exchange column (``q_gwf``); it keeps
    MF6's raw aquifer-referenced sign, so gaining reaches are negative.
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
    selected[value_column] = pd.to_numeric(selected[value_column], errors="coerce")
    selected["rlen"] = pd.to_numeric(selected["rlen"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    q_sum = grouped[value_column].sum(min_count=1).reindex(full_index, fill_value=np.nan)
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
    _append_feature_hover_columns(hover, selected, full_index)
    return values.tolist(), hover


def build_lak_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    value_column: str = "q_lake",
) -> tuple[list[float], dict[str, list]]:
    """Convert LAK exchange rows into a cell choropleth normalized by area.

    The mapped value is ``sum(q) / sum(flow_area)`` within each cell, giving a
    lake-groundwater exchange intensity with units of length per time.
    ``value_column`` is the frame's exchange column (``q_lake``); it keeps MF6's
    raw feature-referenced sign, so a gaining lake is positive.
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
    selected[value_column] = pd.to_numeric(selected[value_column], errors="coerce")
    selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    q_sum = grouped[value_column].sum(min_count=1).reindex(full_index, fill_value=np.nan)
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
    _append_feature_hover_columns(hover, selected, full_index)
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


def _xy_panels_bounds(panels):
    """Global finite x/y bounds across ``[(label, [(name, x, y)])]`` panels."""

    xs_parts, ys_parts = [], []
    for _, lines in panels:
        for _, x_values, y_values in lines:
            x_array = np.asarray(x_values, dtype=float)
            y_array = np.asarray(y_values, dtype=float)
            xs_parts.append(x_array[np.isfinite(x_array)])
            ys_parts.append(y_array[np.isfinite(y_array)])
    xs_all = np.concatenate(xs_parts) if xs_parts else np.asarray([])
    ys_all = np.concatenate(ys_parts) if ys_parts else np.asarray([])
    x_bounds = (float(xs_all.min()), float(xs_all.max())) if xs_all.size else None
    y_bounds = (float(ys_all.min()), float(ys_all.max())) if ys_all.size else None
    return x_bounds, y_bounds


def _xy_mosaic_plotly(panels, *, ncols, title, xaxis_title, yaxis_title, markers):
    """Shared-scale grid of xy panels (series or section profiles)."""

    if not panels:
        raise ValueError("mosaic requires at least one panel.")
    labels = [label for label, _ in panels]
    ncols = min(int(ncols), len(panels)) or 1
    nrows = int(np.ceil(len(panels) / ncols))
    fig = subplots(
        nrows,
        ncols,
        shared_xaxes="all",
        shared_yaxes="all",
        subplot_titles=labels,
    )
    mode = "lines+markers" if markers else "lines"
    seen = set()
    for index, (_, lines) in enumerate(panels):
        row, col = index // ncols + 1, index % ncols + 1
        for name, x_values, y_values in lines:
            fig.add_scatter(
                x=x_values,
                y=y_values,
                mode=mode,
                name=str(name),
                legendgroup=str(name),
                showlegend=str(name) not in seen,
                row=row,
                col=col,
            )
            seen.add(str(name))
    fig.update_xaxes(title_text=xaxis_title, row=nrows)
    fig.update_yaxes(title_text=yaxis_title, col=1)
    fig.update_layout(title=title, uirevision="lock")
    return fig


def _xy_mosaic_mpl(panels, *, ncols, title, xlabel, ylabel, markers):
    """Matplotlib grid of xy panels with shared axes and one legend."""

    if not panels:
        raise ValueError("mosaic requires at least one panel.")
    labels = [label for label, _ in panels]
    ncols = min(int(ncols), len(panels)) or 1
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = mpl_axes(
        nrows,
        ncols,
        figsize=(5.0 * ncols, 3.8 * nrows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    marker = "o" if markers else None
    handles: dict[str, object] = {}
    for axis, label, (_, lines) in zip(axes.flat, labels, panels, strict=False):
        for name, x_values, y_values in lines:
            (handle,) = axis.plot(
                x_values, y_values, marker=marker, linewidth=1.8, label=str(name)
            )
            handles.setdefault(str(name), handle)
        axis.set_title(str(label))
    for axis in axes.flat[len(panels):]:
        axis.set_visible(False)
    for axis in axes[-1, :]:
        axis.set_xlabel(xlabel)
    for axis in axes[:, 0]:
        axis.set_ylabel(ylabel)
    if handles:
        fig.legend(handles.values(), handles.keys(), loc="upper right")
    fig.suptitle(title)
    return fig


def _xy_animation_plotly(frames, *, title, xaxis_title, yaxis_title, markers):
    """Plotly play/slider animation over ``[(label, [(name, x, y)])]`` frames."""

    if not frames:
        raise ValueError("animate requires at least one frame.")
    x_bounds, y_bounds = _xy_panels_bounds(frames)
    mode = "lines+markers" if markers else "lines"

    def _traces(lines):
        """Build one ``go.Scatter`` per ``(name, x, y)`` line for a single frame."""

        return [
            go.Scatter(x=x_values, y=y_values, mode=mode, name=str(name))
            for name, x_values, y_values in lines
        ]

    names = [str(label) for label, _ in frames]
    fig = Fig(data=_traces(frames[0][1]))
    fig.frames = [
        go.Frame(data=_traces(lines), name=str(label)) for label, lines in frames
    ]
    play = {"frame": {"duration": 600, "redraw": True}, "fromcurrent": True}
    pause = {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}
    fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        uirevision="lock",
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "buttons": [
                    {"label": "Play", "method": "animate", "args": [None, play]},
                    {"label": "Pause", "method": "animate", "args": [[None], pause]},
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "steps": [
                    {
                        "method": "animate",
                        "args": [[name], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                        "label": name,
                    }
                    for name in names
                ],
            }
        ],
    )
    # freeze axes so the animation does not rescale frame to frame
    if x_bounds is not None:
        span = (x_bounds[1] - x_bounds[0]) or 1.0
        fig.update_xaxes(range=[x_bounds[0] - 0.05 * span, x_bounds[1] + 0.05 * span])
    if y_bounds is not None:
        span = (y_bounds[1] - y_bounds[0]) or 1.0
        fig.update_yaxes(range=[y_bounds[0] - 0.05 * span, y_bounds[1] + 0.05 * span])
    return fig


def _xy_animation_mpl(frames, *, title, xlabel, ylabel, markers):
    """Matplotlib ``FuncAnimation`` over ``[(label, [(name, x, y)])]`` frames."""

    from matplotlib.animation import FuncAnimation

    if not frames:
        raise ValueError("animate requires at least one frame.")
    labels = [str(label) for label, _ in frames]
    x_bounds, y_bounds = _xy_panels_bounds(frames)
    marker = "o" if markers else None
    fig, axis = mpl_axes(1, 1, figsize=(8.0, 5.0))

    def _draw(index):
        """Redraw the axes for animation frame ``index`` with fixed shared limits."""

        axis.clear()
        for name, x_values, y_values in frames[index][1]:
            axis.plot(x_values, y_values, marker=marker, linewidth=1.8, label=str(name))
        axis.set_title(f"{title} - {labels[index]}")
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        if frames[index][1]:
            axis.legend(loc="best")
        if x_bounds is not None:
            span = (x_bounds[1] - x_bounds[0]) or 1.0
            axis.set_xlim(x_bounds[0] - 0.05 * span, x_bounds[1] + 0.05 * span)
        if y_bounds is not None:
            span = (y_bounds[1] - y_bounds[0]) or 1.0
            axis.set_ylim(y_bounds[0] - 0.05 * span, y_bounds[1] + 0.05 * span)

    _draw(0)
    return FuncAnimation(fig, _draw, frames=len(frames), interval=600, blit=False)


class SpatialView:
    """Unified view grammar for cell-mapped package fields.

    Panel verbs :meth:`map` (plan view) and :meth:`plot` (series view) plus
    composers :meth:`mosaic` and :meth:`animate`, exposed identically on
    single-model, model-group, and diff surfaces (inputs and results). A host
    supplies the atomic styled panel via :meth:`_spatial_map` (and optionally
    the series hooks) plus the dimensions it has; this mixin composes single
    panels, faceted mosaics (by ``layer`` or ``model``), and animations (over
    ``period`` or ``model``) for both backends.

    ``backend="plotly"`` (default) returns an interactive figure/``Choro``;
    ``backend="mpl"`` returns a matplotlib ``Figure`` (or, from :meth:`animate`,
    a ``FuncAnimation``). The mapped *content* -- raw value vs difference -- is
    decided by the node: a ``diff.*`` node maps deltas, an inputs/results node
    maps raw values.
    """

    _PLOTLY_BACKENDS = ("plotly", "interactive")
    _MPL_BACKENDS = ("mpl", "matplotlib", "static")

    # -- hooks the host provides ------------------------------------------
    def _spatial_map(self, *, per, layer, model=None, **kwargs):
        """Return one styled Plotly ``Choro`` for a (per, layer, model).

        Default: delegate to the host's atomic ``map(per=, layer=, ...)`` builder
        -- the single-model surfaces, whose ``map`` already returns a ``Choro``,
        so they gain :meth:`mosaic`/:meth:`animate` for free. Group/diff surfaces
        -- which need a ``model`` selector plus their own delta/colorscale policy
        -- override this hook and let :meth:`map` come from this mixin.
        """

        return self.map(per=int(per), layer=int(layer), **kwargs)

    def _spatial_reference_model(self):
        """Model whose grid/layers frame the view (single-model default)."""

        return self.model

    def _spatial_layers(self) -> list[int]:
        """Layers to offer for faceting: those present in the data, else every grid layer."""

        # Prefer the layers actually present in this field's data (many BC
        # packages live in a single layer); fall back to every grid layer.
        frame = self.get()
        columns = getattr(frame, "columns", [])
        if "layer" in columns and not frame.empty:
            present = sorted({int(value) for value in frame["layer"].dropna().tolist()})
            if present:
                return present
        return list(range(int(self._spatial_reference_model().gwf.modelgrid.nlay)))

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the data (``[0]`` when there is no period axis)."""

        frame = self.get()
        columns = getattr(frame, "columns", [])
        if "per" not in columns or frame.empty:
            return [0]
        return sorted({int(value) for value in frame["per"].dropna().tolist()})

    def _spatial_models(self) -> list[str] | None:
        """Model names for the model axis; ``None`` for a single-model surface."""

        return None  # single-model surfaces have no model axis

    def _spatial_default_facet(self) -> str:
        """Default mosaic facet: ``"model"`` for a group, else ``"layer"``."""

        return "model" if self._spatial_models() else "layer"

    def _spatial_value_label(self) -> str:
        """Human label for the mapped quantity (field/value/package name, else "value")."""

        for attribute in ("field_name", "value_name", "package_name"):
            value = getattr(self, attribute, None)
            if value is not None:
                return str(value)
        return "value"

    def _spatial_is_diff(self) -> bool:
        """Whether panels show a signed difference (drives diverging, zero-centered scales)."""

        return False

    # -- series hooks (the plot verb) ---------------------------------------
    #: entity id columns that get one line each in ``plot`` when present
    _SERIES_ENTITY_COLUMNS = ("lake", "reach", "well")

    def _series_table(self) -> pd.DataFrame:
        """Long table backing :meth:`plot`; defaults to the node's ``get()``."""

        get = getattr(self, "get", None)
        if get is None:
            raise NotImplementedError(
                "plot() requires a tabular get(); this node does not provide one."
            )
        return get()

    def _series_value_column(self, frame: pd.DataFrame) -> str:
        """Column of :meth:`_series_table` that :meth:`plot` draws."""

        value = getattr(self, "value_name", None)
        if value is not None and value in frame.columns:
            return str(value)
        label = self._spatial_value_label()
        if label in frame.columns:
            return label
        raise KeyError(
            "Could not infer the value column for plot(); available columns: "
            f"{list(frame.columns)}."
        )

    def _series_default_agg(self) -> str:
        """How :meth:`plot` collapses cells within a line (sum for fluxes)."""

        return "sum"

    # -- dimension resolution ---------------------------------------------
    def _resolve_layers(self, layer) -> list[int]:
        """Normalize the ``layer`` selector to a list: ``None`` -> all present, int -> one, iterable -> many."""

        if layer is None:
            return self._spatial_layers()
        if isinstance(layer, (int, np.integer)):
            return [int(layer)]
        resolved = [int(value) for value in layer]
        if not resolved:
            raise ValueError("layer must select at least one layer.")
        return resolved

    def _single_layer(self, layer) -> int:
        """Reduce the ``layer`` selector to one layer (0 by default; first of a set)."""

        if layer is None:
            return 0
        if isinstance(layer, (int, np.integer)):
            return int(layer)
        return self._resolve_layers(layer)[0]

    def _resolve_periods(self, per) -> list[int]:
        """Normalize the ``per`` selector to a list: ``None`` -> all present, int -> one, iterable -> many."""

        if per is None:
            return self._spatial_periods()
        if isinstance(per, (int, np.integer)):
            return [int(per)]
        resolved = [int(value) for value in per]
        if not resolved:
            raise ValueError("per must select at least one period.")
        return resolved

    def _resolve_models(self, model) -> list:
        """Normalize the ``model`` selector against this surface's model axis.

        Single-model surfaces accept only ``None``/``"*"`` (returning ``[None]``);
        group surfaces map ``None``/``"*"`` to all models and validate any explicit
        names against the group.
        """

        models = self._spatial_models()
        if not models:
            if model not in (None, "*"):
                raise ValueError(
                    "This surface has a single model; 'model' is not applicable."
                )
            return [None]
        if model in (None, "*"):
            return list(models)
        names = [str(model)] if isinstance(model, str) else [str(m) for m in model]
        missing = [name for name in names if name not in models]
        if missing:
            raise KeyError(f"Models not in this group: {missing!r}")
        return names

    @classmethod
    def _normalize_backend(cls, backend: str) -> str:
        """Canonicalize a backend alias to ``"plotly"`` or ``"mpl"`` (raises if unknown)."""

        value = str(backend).lower()
        if value in cls._PLOTLY_BACKENDS:
            return "plotly"
        if value in cls._MPL_BACKENDS:
            return "mpl"
        raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")

    def map(
        self,
        model=None,
        *,
        per: int = 0,
        layer: int = 0,
        backend: str = "plotly",
        **map_kwargs,
    ):
        """One choropleth panel.

        ``backend="plotly"`` (default) returns an interactive ``Choro``;
        ``backend="mpl"`` returns a static matplotlib ``Figure``. On a
        group/diff surface ``model`` (first positional) selects which model to
        draw (defaults to the reference).
        """

        choro = self._spatial_map(per=int(per), layer=int(layer), model=model, **map_kwargs)
        if self._normalize_backend(backend) == "plotly":
            return choro
        return choro.plot_mpl()

    @staticmethod
    def _series_line_label(keys, key) -> str:
        """Human label for one plotted line (model / entity / layer / cell)."""

        parts = []
        for column, value in zip(keys, key, strict=False):
            if pd.isna(value):
                continue
            if column == "model":
                parts.append(str(value))
            elif column == "layer":
                parts.append(f"L{int(value)}")
            elif column == "cell":
                parts.append(f"C{int(value)}")
            else:  # entity columns: lake, reach, well
                parts.append(f"{column.capitalize()} {int(value)}")
        return " / ".join(parts) if parts else "All"

    def _series_selection(self, *, model=None, per=None, layer=None, cells=None):
        """Filter the series table; return ``(frame, value_column, line_keys)``."""

        selected_models = self._resolve_models(model)
        frame = self._series_table()
        if frame is None or "per" not in getattr(frame, "columns", []):
            raise ValueError(
                "plot() requires period-indexed rows; this node has no 'per' data."
            )
        frame = frame.copy()

        # column-based filters (a filter is skipped when its column is absent)
        if model is not None and "model" in frame.columns:
            frame = frame[frame["model"].isin([str(name) for name in selected_models])]
        if per is not None:
            per_values = (
                [int(per)] if isinstance(per, (int, np.integer)) else [int(v) for v in per]
            )
            frame = frame[frame["per"].isin(per_values)]
        if layer is not None and "layer" in frame.columns:
            layer_values = (
                [int(layer)]
                if isinstance(layer, (int, np.integer))
                else [int(v) for v in layer]
            )
            frame = frame[frame["layer"].isin(layer_values)]
        if cells is not None and "cell" in frame.columns:
            cell_values = (
                [int(cells)]
                if isinstance(cells, (int, np.integer))
                else [int(v) for v in cells]
            )
            frame = frame[frame["cell"].isin(cell_values)]

        value_column = self._series_value_column(frame)
        frame[value_column] = pd.to_numeric(frame[value_column], errors="coerce")

        # line keys: model axis, entity ids, multi-layer split, explicit cells.
        # A single-model surface may still carry a ``model`` column in its
        # table; only treat it as a line axis when the surface has one.
        keys = ["model"] if self._spatial_models() and "model" in frame.columns else []
        keys += [c for c in self._SERIES_ENTITY_COLUMNS if c in frame.columns]
        if "layer" in frame.columns and frame["layer"].nunique() > 1:
            keys.append("layer")
        if cells is not None and "cell" in frame.columns:
            keys.append("cell")
        return frame, value_column, keys

    def _series_lines(self, frame, keys, value_column, aggregation):
        """Aggregate a filtered series frame into ``[(label, x, y)]`` lines."""

        grouped = frame.groupby([*keys, "per"], dropna=False, as_index=False)[
            value_column
        ].agg(aggregation)
        lines = []
        if keys:
            for key, sub in grouped.groupby(keys, dropna=False):
                key = key if isinstance(key, tuple) else (key,)
                sub = sub.sort_values("per")
                lines.append(
                    (
                        self._series_line_label(keys, key),
                        sub["per"].astype(int).to_numpy(),
                        sub[value_column].astype(float).to_numpy(),
                    )
                )
        elif not grouped.empty:
            sub = grouped.sort_values("per")
            lines.append(
                (
                    self._spatial_value_label(),
                    sub["per"].astype(int).to_numpy(),
                    sub[value_column].astype(float).to_numpy(),
                )
            )
        return lines

    def plot(
        self,
        model=None,
        *,
        per=None,
        layer=None,
        cells=None,
        agg: str | None = None,
        backend: str = "plotly",
        title: str | None = None,
    ):
        """Series panel: this field by stress period.

        One line per model (groups/diffs), per entity (lake/reach) when the
        package has one, and per cell when ``cells`` selects specific cells;
        otherwise cells are aggregated (``agg``, default sum for fluxes).
        ``backend="plotly"`` (default) returns an interactive ``viz.Fig``;
        ``backend="mpl"`` a matplotlib ``Figure``. Stress periods are zero-based.
        """

        backend_kind = self._normalize_backend(backend)
        frame, value_column, keys = self._series_selection(
            model=model, per=per, layer=layer, cells=cells
        )
        aggregation = self._series_default_agg() if agg is None else str(agg)
        lines = self._series_lines(frame, keys, value_column, aggregation)
        heading = title or (
            f"{self._spatial_value_label()} by stress period"
            + (" (model - reference)" if self._spatial_is_diff() else "")
        )
        if backend_kind == "plotly":
            fig = Fig()
            for label, xs_values, ys_values in lines:
                fig.add_scatter(
                    x=xs_values, y=ys_values, mode="lines+markers", name=label
                )
            fig.update_layout(
                title=heading, xaxis_title="Stress Period", yaxis_title=value_column
            )
            return fig
        fig, ax = mpl_axes(figsize=(8, 4))
        for label, xs_values, ys_values in lines:
            ax.plot(xs_values, ys_values, marker="o", linewidth=2.0, label=label)
        ax.set_title(heading)
        ax.set_xlabel("Stress Period")
        ax.set_ylabel(value_column)
        if lines:
            ax.legend()
        fig.tight_layout()
        return fig

    # -- cross-section hooks (the section verb + kind="section" composers) --
    def _sections(self, model=None, **kwargs):
        """Return ``{name: XSection}`` for this node (heads leaves override)."""

        raise NotImplementedError(
            "This node has no cross-section view; section/kind='section' is available on "
            "heads leaves (model.hds, group.hds, diff heads)."
        )

    def section(
        self,
        model=None,
        *,
        line=None,
        cells=None,
        per: int | None = None,
        layer=0,
        backend: str = "plotly",
        title: str | None = None,
        **kwargs,
    ):
        """Cross-section panel: this field vs distance along a section line.

        ``line`` is a shapely ``LineString`` / ``(x, y)`` pairs / flopy dict;
        alternatively ``cells`` traces the section through cell centroids. On a
        group surface ``model=None`` overlays every member; a diff overlays the
        reference with each compared model. ``backend="plotly"`` returns a
        ``viz.Fig``, ``backend="mpl"`` a matplotlib figure.
        """

        from myflopy.modflow.utils.datatypes.xsections import render_xsections

        sections = self._sections(
            model=model, line=line, cells=cells, per=per, layer=layer, **kwargs
        )
        return render_xsections(sections, backend=backend, title=title)

    @staticmethod
    def _section_lines(sections) -> list:
        """Flatten sections to ``[(series_label, x, y)]`` at their current period."""

        from myflopy.modflow.utils.datatypes.xsections import combined_section_frame

        data = combined_section_frame(sections)
        return [
            (
                str(name),
                sub["distance"].to_numpy(dtype=float),
                sub["elevation"].to_numpy(dtype=float),
            )
            for name, sub in data.groupby("series", sort=False)
        ]

    @staticmethod
    def _period_end_kstpkper(model, period: int):
        """Return the period-end ``(kstp, kper)`` saved for ``period`` (or None)."""

        candidates = [
            tuple(int(v) for v in key)
            for key in model.kstpkper
            if int(key[1]) == int(period)
        ]
        return max(candidates) if candidates else None

    def mosaic(
        self,
        *,
        kind: str = "map",
        by: str | None = None,
        per=None,
        layer=None,
        model=None,
        cells=None,
        agg: str | None = None,
        backend: str = "plotly",
        ncols: int = 3,
        title: str | None = None,
        sync_views: bool = True,
        **kwargs,
    ):
        """Faceted grid of panels sharing one scale.

        ``kind`` picks the panel type: ``"map"`` (choropleths, default),
        ``"plot"`` (series panels), or ``"section"`` (cross sections, heads leaves).
        ``by`` picks the facet axis -- maps: ``"layer"``/``"model"``; series:
        ``"model"``/``"layer"``/an entity (``"lake"``, ``"reach"``); sections:
        ``"model"``/``"period"``. Defaults to the surface's natural axis.
        For ``kind="map"``, ``sync_views`` (default ``True``) frames every panel
        to one shared extent so the maps line up; set it ``False`` to fit each
        map to its own data. It is ignored for non-map kinds.
        """

        kind_key = str(kind).lower()
        if kind_key == "map":
            axis = (by or self._spatial_default_facet()).lower()
            fixed_per = self._single_period(per)
            panels = self._facet_panels(
                axis, per=fixed_per, layer=layer, model=model, **kwargs
            )
            heading = title or f"{self._spatial_value_label()} by {axis}"
            if self._normalize_backend(backend) == "plotly":
                return self._plotly_mosaic(
                    panels, ncols=int(ncols), title=heading, sync_views=sync_views
                )
            return self._mpl_mosaic(panels, ncols=int(ncols), title=heading)
        if kind_key == "plot":
            return self._plot_mosaic(
                by=by,
                per=per,
                layer=layer,
                model=model,
                cells=cells,
                agg=agg,
                backend=backend,
                ncols=ncols,
                title=title,
            )
        if kind_key == "section":
            return self._xs_mosaic(
                by=by,
                per=per,
                layer=layer,
                model=model,
                backend=backend,
                ncols=ncols,
                title=title,
                **kwargs,
            )
        raise ValueError(f"kind must be 'map', 'plot', or 'section', got {kind!r}.")

    def _plot_mosaic(
        self, *, by, per, layer, model, cells, agg, backend, ncols, title
    ):
        """Grid of series panels faceted by model, layer, or an entity."""

        frame, value_column, keys = self._series_selection(
            model=model, per=per, layer=layer, cells=cells
        )
        default_axis = (
            "model"
            if (self._spatial_models() and "model" in frame.columns)
            else next(
                (c for c in self._SERIES_ENTITY_COLUMNS if c in frame.columns),
                "layer",
            )
        )
        axis = str(by or default_axis).lower()
        if axis == "period":
            raise ValueError(
                "Series panels already have period on the x-axis; facet by "
                "'model', 'layer', or an entity (e.g. 'lake') instead."
            )
        if axis not in frame.columns:
            raise ValueError(
                f"Cannot facet series by {axis!r}; available facet columns: "
                f"{[c for c in ('model', 'layer', *self._SERIES_ENTITY_COLUMNS) if c in frame.columns]}."
            )
        if axis == "model":
            in_frame = set(frame["model"])
            facet_values = [n for n in self._resolve_models(model) if n in in_frame]
        else:
            facet_values = sorted(frame[axis].dropna().unique().tolist())
        aggregation = self._series_default_agg() if agg is None else str(agg)
        panel_keys = [key for key in keys if key != axis]
        panels = [
            (
                self._series_line_label([axis], (value,)),
                self._series_lines(
                    frame[frame[axis] == value], panel_keys, value_column, aggregation
                ),
            )
            for value in facet_values
        ]
        heading = title or (
            f"{self._spatial_value_label()} by stress period per {axis}"
            + (" (model - reference)" if self._spatial_is_diff() else "")
        )
        if self._normalize_backend(backend) == "plotly":
            return _xy_mosaic_plotly(
                panels,
                ncols=int(ncols),
                title=heading,
                xaxis_title="Stress Period",
                yaxis_title=value_column,
                markers=True,
            )
        return _xy_mosaic_mpl(
            panels,
            ncols=int(ncols),
            title=heading,
            xlabel="Stress Period",
            ylabel=value_column,
            markers=True,
        )

    def _xs_mosaic(
        self, *, by, per, layer, model, backend, ncols, title, **kwargs
    ):
        """Grid of cross-section panels faceted by model or period."""

        axis = str(by or ("model" if self._spatial_models() else "period")).lower()
        base_layer = 0 if layer is None else layer
        panels = []
        if axis == "model":
            for name in self._resolve_models(model):
                sections = self._sections(
                    model=name, per=self._single_period(per), layer=base_layer, **kwargs
                )
                panels.append((str(name), self._section_lines(sections)))
        elif axis == "period":
            sections = self._sections(model=model, layer=base_layer, **kwargs)
            for period in self._resolve_periods(per):
                for section in sections.values():
                    key = self._period_end_kstpkper(section.model, period)
                    if key is not None:
                        section.kstpkper = key
                panels.append((f"Period {period}", self._section_lines(sections)))
        else:
            raise ValueError(f"by must be 'model' or 'period' for kind='section', got {axis!r}.")
        heading = title or f"{self._spatial_value_label()} cross sections by {axis}"
        if self._normalize_backend(backend) == "plotly":
            return _xy_mosaic_plotly(
                panels,
                ncols=int(ncols),
                title=heading,
                xaxis_title="Distance",
                yaxis_title="Elevation",
                markers=False,
            )
        return _xy_mosaic_mpl(
            panels,
            ncols=int(ncols),
            title=heading,
            xlabel="Distance",
            ylabel="Elevation",
            markers=False,
        )

    def _facet_panels(self, axis, *, per, layer, model, **map_kwargs):
        """Build ``[(label, Choro), ...]`` mosaic panels faceted by ``"layer"`` or ``"model"``.

        A layer facet fixes the model (reference) and varies the layer; a model
        facet fixes the layer and varies the model. Both hold ``per`` constant.
        """

        if axis == "layer":
            reference = None if self._spatial_models() is None else self._resolve_models(model)[0]
            return [
                (
                    f"Layer {value + 1}",
                    self._spatial_map(per=int(per), layer=value, model=reference, **map_kwargs),
                )
                for value in self._resolve_layers(layer)
            ]
        if axis == "model":
            base_layer = self._single_layer(layer)
            return [
                (
                    str(name),
                    self._spatial_map(per=int(per), layer=base_layer, model=name, **map_kwargs),
                )
                for name in self._resolve_models(model)
            ]
        raise ValueError(f"by must be 'layer' or 'model', got {axis!r}.")

    def _plotly_mosaic(self, panels, *, ncols, title, sync_views=True):
        """Compose ``Choro`` panels into one synchronized Plotly small-multiples figure."""

        return viz_mosaic(
            list(panels),
            ncols=int(ncols),
            title=title,
            diff=self._spatial_is_diff(),
            sync_views=sync_views,
        )

    def _mpl_mosaic(self, panels, *, ncols, title):
        """Render ``Choro`` panels as a matplotlib grid sharing one color scale + colorbar.

        Diff surfaces get a symmetric ``RdBu`` scale centered at zero; other
        surfaces get a ``viridis`` scale spanning the panels' finite value range.
        """

        from matplotlib import cm
        from matplotlib import colors as mcolors

        labels = [label for label, _ in panels]
        choros = [choro for _, choro in panels]
        if not choros:
            raise ValueError("mosaic requires at least one panel.")
        ncols = min(int(ncols), len(choros)) or 1
        nrows = int(np.ceil(len(choros) / ncols))
        fig, axes = mpl_axes(nrows, ncols, figsize=(5.0 * ncols, 4.5 * nrows), squeeze=False)
        arrays = [np.asarray(choro.zs, dtype=float) for choro in choros]
        finite_parts = [values[np.isfinite(values)] for values in arrays if np.isfinite(values).any()]
        finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
        if self._spatial_is_diff():
            absmax = (float(np.nanmax(np.abs(finite))) if finite.size else 1.0) or 1.0
            vmin, vmax, cmap = -absmax, absmax, "RdBu"
        else:
            vmin = float(np.nanmin(finite)) if finite.size else 0.0
            vmax = float(np.nanmax(finite)) if finite.size else 1.0
            cmap = "viridis"
        for axis, label, choro, values in zip(axes.flat, labels, choros, arrays, strict=False):
            gdf = choro.vor.gdf_vorPolys.copy()
            gdf["_spatial_value"] = values
            gdf.plot(
                column="_spatial_value",
                ax=axis,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                linewidth=0.2,
                edgecolor="#666666",
            )
            axis.set_title(str(label))
            axis.set_axis_off()
            axis.set_aspect("equal")
        for axis in axes.flat[len(choros):]:
            axis.set_visible(False)
        mappable = cm.ScalarMappable(norm=mcolors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        mappable.set_array([])
        fig.colorbar(
            mappable,
            ax=list(axes.flat[: len(choros)]),
            shrink=0.85,
            label=self._spatial_value_label(),
        )
        fig.suptitle(title)
        return fig

    def _single_period(self, per) -> int:
        """Reduce the ``per`` selector to one period (0 by default; first of a set)."""

        if per is None:
            return 0
        if isinstance(per, (int, np.integer)):
            return int(per)
        return self._resolve_periods(per)[0]

    def _frame_panels(self, axis, *, per, layer, model, **map_kwargs):
        """Build ``[(label, Choro), ...]`` animation frames over ``"period"`` or ``"model"``.

        A period animation fixes the model and layer and steps through periods; a
        model animation fixes the period and layer and steps through models.
        """

        base_layer = self._single_layer(layer)
        if axis == "period":
            reference = self._resolve_models(model)[0]
            return [
                (
                    f"Period {value}",
                    self._spatial_map(per=value, layer=base_layer, model=reference, **map_kwargs),
                )
                for value in self._resolve_periods(per)
            ]
        if axis == "model":
            fixed_per = self._single_period(per)
            return [
                (
                    str(name),
                    self._spatial_map(per=fixed_per, layer=base_layer, model=name, **map_kwargs),
                )
                for name in self._resolve_models(model)
            ]
        raise ValueError(f"over must be 'period' or 'model', got {axis!r}.")

    def animate(
        self,
        *,
        kind: str = "map",
        over: str | None = None,
        per=None,
        layer: int = 0,
        model=None,
        cells=None,
        agg: str | None = None,
        backend: str = "plotly",
        title: str | None = None,
        **kwargs,
    ):
        """Animate a panel across a dimension.

        ``kind`` picks the panel type: ``"map"`` (default), ``"section"`` (cross
        sections, heads leaves), or ``"plot"`` (series). ``over="period"``
        sweeps stress periods (maps and sections); ``over="model"`` sweeps the
        group's models (any kind). ``backend="plotly"`` returns an interactive
        house ``viz.Fig`` (pan-drag, scroll-zoom, house template) with
        play/slider; ``backend="mpl"`` a matplotlib ``FuncAnimation``. Both are
        objects for the caller to display or save -- nothing is written to
        disk.
        """

        kind_key = str(kind).lower()
        if kind_key == "map":
            axis = str(over or "period").lower()
            frames = self._frame_panels(axis, per=per, layer=layer, model=model, **kwargs)
            heading = title or f"{self._spatial_value_label()} over {axis}"
            if self._normalize_backend(backend) == "plotly":
                return self._plotly_animation(frames, title=heading)
            return self._mpl_animation(frames, title=heading)
        if kind_key == "plot":
            return self._plot_animation(
                over=over,
                per=per,
                layer=layer,
                model=model,
                cells=cells,
                agg=agg,
                backend=backend,
                title=title,
            )
        if kind_key == "section":
            return self._xs_animation(
                over=over,
                per=per,
                layer=layer,
                model=model,
                backend=backend,
                title=title,
                **kwargs,
            )
        raise ValueError(f"kind must be 'map', 'plot', or 'section', got {kind!r}.")

    def _plot_animation(
        self, *, over, per, layer, model, cells, agg, backend, title
    ):
        """Flip series panels over the group's models."""

        axis = str(over or "model").lower()
        if axis != "model":
            raise ValueError(
                "Series panels already sweep periods on the x-axis; use "
                "over='model', or animate kind='map'/'section' over periods."
            )
        frame, value_column, keys = self._series_selection(
            model=model, per=per, layer=layer, cells=cells
        )
        if not (self._spatial_models() and "model" in frame.columns):
            raise ValueError("over='model' requires a model axis (groups/diffs).")
        aggregation = self._series_default_agg() if agg is None else str(agg)
        panel_keys = [key for key in keys if key != "model"]
        frames = [
            (
                str(name),
                self._series_lines(
                    frame[frame["model"] == name], panel_keys, value_column, aggregation
                ),
            )
            for name in self._resolve_models(model)
        ]
        heading = title or (
            f"{self._spatial_value_label()} by stress period over model"
            + (" (model - reference)" if self._spatial_is_diff() else "")
        )
        if self._normalize_backend(backend) == "plotly":
            return _xy_animation_plotly(
                frames,
                title=heading,
                xaxis_title="Stress Period",
                yaxis_title=value_column,
                markers=True,
            )
        return _xy_animation_mpl(
            frames,
            title=heading,
            xlabel="Stress Period",
            ylabel=value_column,
            markers=True,
        )

    def _xs_animation(self, *, over, per, layer, model, backend, title, **kwargs):
        """Animate cross sections over stress periods or the group's models."""

        axis = str(over or "period").lower()
        base_layer = 0 if layer is None else layer
        frames = []
        if axis == "period":
            sections = self._sections(model=model, layer=base_layer, **kwargs)
            for period in self._resolve_periods(per):
                for section in sections.values():
                    key = self._period_end_kstpkper(section.model, period)
                    if key is not None:
                        section.kstpkper = key
                frames.append((f"Period {period}", self._section_lines(sections)))
        elif axis == "model":
            for name in self._resolve_models(model):
                sections = self._sections(
                    model=name, per=self._single_period(per), layer=base_layer, **kwargs
                )
                frames.append((str(name), self._section_lines(sections)))
        else:
            raise ValueError(
                f"over must be 'period' or 'model' for kind='section', got {axis!r}."
            )
        heading = title or f"{self._spatial_value_label()} cross section over {axis}"
        if self._normalize_backend(backend) == "plotly":
            return _xy_animation_plotly(
                frames,
                title=heading,
                xaxis_title="Distance",
                yaxis_title="Elevation",
                markers=False,
            )
        return _xy_animation_mpl(
            frames,
            title=heading,
            xlabel="Distance",
            ylabel="Elevation",
            markers=False,
        )

    def _plotly_animation(self, frames, *, title):
        """Build a Plotly play/slider map animation over ``[(label, Choro), ...]`` frames.

        All frames share one data-fitted map view so the map does not reset to a
        world view between frames.
        """

        if not frames:
            raise ValueError("animate requires at least one frame.")
        names = [str(label) for label, _ in frames]
        traces = [choro.get_choropleth() for _, choro in frames]
        fig = Fig(data=[traces[0]])
        fig.frames = [
            go.Frame(data=[trace], name=name)
            for name, trace in zip(names, traces, strict=False)
        ]
        # A single map flipped across frames still needs its view fitted to the
        # data -- otherwise it renders zoomed out to the world, like the mosaic
        # subplots did. All frames share the site, so one shared view suffices.
        map_view = shared_map_view([choro for _, choro in frames])
        play = {"frame": {"duration": 600, "redraw": True}, "fromcurrent": True}
        pause = {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}
        fig.update_layout(
            title=title,
            uirevision="lock",
            map=map_view or {},
            updatemenus=[
                {
                    "type": "buttons",
                    "showactive": False,
                    "buttons": [
                        {"label": "Play", "method": "animate", "args": [None, play]},
                        {"label": "Pause", "method": "animate", "args": [[None], pause]},
                    ],
                }
            ],
            sliders=[
                {
                    "active": 0,
                    "steps": [
                        {
                            "method": "animate",
                            "args": [[name], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                            "label": name,
                        }
                        for name in names
                    ],
                }
            ],
        )
        return fig

    def _mpl_animation(self, frames, *, title):
        """Build a matplotlib ``FuncAnimation`` over ``[(label, Choro), ...]`` map frames.

        Uses one fixed color scale across all frames (symmetric ``RdBu`` for diff
        surfaces, ``viridis`` otherwise).
        """

        from matplotlib.animation import FuncAnimation

        if not frames:
            raise ValueError("animate requires at least one frame.")
        labels = [str(label) for label, _ in frames]
        choros = [choro for _, choro in frames]
        arrays = [np.asarray(choro.zs, dtype=float) for choro in choros]
        finite_parts = [values[np.isfinite(values)] for values in arrays if np.isfinite(values).any()]
        finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
        if self._spatial_is_diff():
            absmax = (float(np.nanmax(np.abs(finite))) if finite.size else 1.0) or 1.0
            vmin, vmax, cmap = -absmax, absmax, "RdBu"
        else:
            vmin = float(np.nanmin(finite)) if finite.size else 0.0
            vmax = float(np.nanmax(finite)) if finite.size else 1.0
            cmap = "viridis"
        fig, axis = mpl_axes(1, 1, figsize=(7.0, 6.0))

        def _draw(index):
            """Redraw the map for animation frame ``index`` under the shared color scale."""

            axis.clear()
            gdf = choros[index].vor.gdf_vorPolys.copy()
            gdf["_spatial_value"] = arrays[index]
            gdf.plot(
                column="_spatial_value",
                ax=axis,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                linewidth=0.2,
                edgecolor="#666666",
            )
            axis.set_title(f"{title} - {labels[index]}")
            axis.set_axis_off()
            axis.set_aspect("equal")

        _draw(0)
        return FuncAnimation(fig, _draw, frames=len(frames), interval=600, blit=False)


class LeafFieldSugar:
    """``field=`` keyword and ``.field(name)`` nodes on an aggregator *leaf*.

    For hosts that are themselves a :class:`SpatialView` leaf (their verbs draw
    a default field) but carry several numeric input fields -- e.g. a grouped
    GHB input accessor with ``bhead`` and ``cond``. ``field=None`` keeps the
    leaf's own verb; ``field="cond"`` (or the ``inputs.cond`` attribute)
    dispatches to a field-pinned node from :meth:`_field_node`. Namespace-style
    hosts without their own verbs use :class:`FieldMappable` instead.

    Place this mixin *before* the SpatialView base in the MRO.
    """

    def _field_names(self) -> list[str]:
        """Return the selectable field names (override per host)."""

        raise NotImplementedError

    def field_names(self) -> list[str]:
        """List the mappable field names available on this node."""

        return list(self._field_names())

    def _field_node(self, name: str):
        """Return this leaf pinned to one field (override per host)."""

        raise NotImplementedError

    def field(self, name: str):
        """Return the field-pinned node for ``name`` (validated)."""

        key = str(name).lower()
        available = self._field_names()
        if available and key not in available:
            raise ValueError(
                f"Unknown field {name!r} for this package; choose from {available}."
            )
        return self._field_node(key)

    #: attributes _field_names itself may read -- never treat these as fields
    #: (prevents __getattr__ recursion when they are genuinely missing)
    _FIELD_SUGAR_GUARD = ("package_name", "group", "model", "field_name")

    def __getattr__(self, name: str):
        """Resolve ``node.<field>`` to that field's pinned node (else raise ``AttributeError``).

        Only names that are actual field names (and not private or guarded
        attributes) are dispatched, so a genuinely missing attribute still errors
        cleanly instead of recursing.
        """

        if (
            not name.startswith("_")
            and name not in self._FIELD_SUGAR_GUARD
            and name in self._field_names()
        ):
            return self._field_node(name)
        raise AttributeError(
            f"{type(self).__name__!s} has no attribute or input field {name!r}"
        )

    def map(self, *args, field=None, **kwargs):
        """Map this leaf's default field, or the one named by ``field=``."""

        if field is not None:
            return self.field(field).map(*args, **kwargs)
        return super().map(*args, **kwargs)

    def plot(self, *args, field=None, **kwargs):
        """Series-plot this leaf's default field, or the one named by ``field=``."""

        if field is not None:
            return self.field(field).plot(*args, **kwargs)
        return super().plot(*args, **kwargs)

    def mosaic(self, *args, field=None, **kwargs):
        """Mosaic this leaf's default field, or the one named by ``field=``."""

        if field is not None:
            return self.field(field).mosaic(*args, **kwargs)
        return super().mosaic(*args, **kwargs)

    def animate(self, *args, field=None, **kwargs):
        """Animate this leaf's default field, or the one named by ``field=``."""

        if field is not None:
            return self.field(field).animate(*args, **kwargs)
        return super().animate(*args, **kwargs)


class DiffSpatialView(SpatialView):
    """:class:`SpatialView` for a diff node: panels are per compared-model deltas.

    The mapped content is a difference (model - reference), so the shared color
    scale is diverging and centered at zero, and faceting defaults to one panel
    per compared model. The host supplies ``_spatial_map(per, layer, model)``
    that returns the delta ``Choro`` for one compared model plus the compared
    model names via ``_spatial_models``.
    """

    def _spatial_is_diff(self) -> bool:
        """Diff nodes always map signed deltas (diverging, zero-centered scale)."""

        return True

    def _spatial_default_facet(self) -> str:
        """Diff mosaics default to one panel per compared model."""

        return "model"


class FieldMappable:
    """Namespace-level ``map`` / ``mosaic`` / ``animate`` with a ``field=`` selector.

    A results/inputs namespace groups several mappable *fields* -- e.g. a lake's
    ``q`` (exchange), ``stage``, ``stage_change``; UZF's ``gwrch`` / ``sat``.
    Each field is a first-class accessor (``ns.stage``) with a :class:`SpatialView`
    surface. This mixin adds the same verbs at the namespace level with a
    ``field=`` selector that simply dispatches to the named accessor -- so
    ``ns.map(field="stage")`` is exactly ``ns.stage.map()``. ``field=None`` uses
    the package's default field; :meth:`field_names` lists the choices.
    """

    _default_field: str = "q"

    def _field_names(self) -> list[str]:
        """Return the selectable field names (override per namespace)."""

        raise NotImplementedError

    def field_names(self) -> list[str]:
        """List the mappable field names available on this namespace."""

        return list(self._field_names())

    def _field_accessor(self, field=None):
        """Return the accessor for ``field`` (or the default field), validating the name."""

        name = field if field is not None else self._default_field
        available = self._field_names()
        if available and name not in available:
            raise ValueError(
                f"Unknown field {field!r} for this package; choose from {available}."
            )
        return getattr(self, name)

    def map(self, *args, field=None, **kwargs):
        """Map one result/input field (``field=`` selects it; default otherwise).

        Extra positional/keyword args pass through to the field's ``map`` -- e.g.
        a group namespace forwards a positional model name: ``results.map("F9b")``.
        """

        return self._field_accessor(field).map(*args, **kwargs)

    def plot(self, *args, field=None, **kwargs):
        """Series plot of one field by stress period (``field=`` selects it)."""

        return self._field_accessor(field).plot(*args, **kwargs)

    def section(self, *args, field=None, **kwargs):
        """Cross-section of one field along a line (``field=`` selects it)."""

        return self._field_accessor(field).section(*args, **kwargs)

    def mosaic(self, *args, field=None, **kwargs):
        """Shared-scale mosaic of one field over layers/models."""

        return self._field_accessor(field).mosaic(*args, **kwargs)

    def animate(self, *args, field=None, **kwargs):
        """Animate one field over periods/models."""

        return self._field_accessor(field).animate(*args, **kwargs)


__all__ = [
    "_symmetric_color_limit",
    "_blue_white_red_diverging_colorscale",
    "_red_white_blue_diverging_colorscale",
    "_as_layer_cell_property",
    "build_sfr_q_map_payload",
    "build_lak_q_map_payload",
    "build_surface_water_q_map_payload",
    "build_cell_input_map_payload",
    "build_group_input_compare_map_payload",
    "SpatialView",
    "DiffSpatialView",
    "FieldMappable",
    "LeafFieldSugar",
]
