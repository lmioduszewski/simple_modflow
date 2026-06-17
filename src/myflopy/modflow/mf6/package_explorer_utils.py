"""Cross-cutting helpers for package explorer modules."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def _get_package_explorer_cache(model: "SimulationBase") -> dict[tuple, pd.DataFrame]:
    """Return a per-model cache for expensive normalized explorer tables."""

    cache = getattr(model, "_package_explorer_cache", None)
    if cache is None:
        cache = {}
        setattr(model, "_package_explorer_cache", cache)
    return cache


def _extract_structured_column(data, column_name: str):
    """Return one column from structured MF6 record data with a safe fallback."""

    dtype = getattr(data, "dtype", None)
    names = getattr(dtype, "names", None)
    if names and column_name in names:
        return np.asarray(data[column_name])
    frame = pd.DataFrame(data).copy()
    if column_name not in frame.columns:
        raise KeyError(column_name)
    return frame[column_name].to_numpy()


def _coerce_numeric_like_columns(
    frame: pd.DataFrame, *, exclude: Iterable[str] = ()
) -> pd.DataFrame:
    """Convert object/string columns to numeric when every non-null value is numeric-like."""

    excluded = {str(column) for column in exclude}
    coerced = frame.copy()
    for column in coerced.columns:
        if str(column) in excluded:
            continue
        series = coerced[column]
        if not (
            pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series)
        ):
            continue
        non_null = series.notna()
        if not non_null.any():
            continue
        converted = pd.to_numeric(series, errors="coerce")
        if converted.loc[non_null].notna().all():
            coerced[column] = converted
    return coerced


def split_cellid_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Expand a FloPy ``cellid`` column into zero-based ``layer`` and ``cell``.

    Parameters
    ----------
    frame
        Input DataFrame that may contain a FloPy-style ``cellid`` object
        column. The function accepts structured tuples such as ``(layer, cell)``
        as well as simpler scalar cell identifiers.

    Returns
    -------
    pandas.DataFrame
        A copy of ``frame`` where ``cellid`` has been replaced with zero-based
        ``layer`` and ``cell`` columns when possible.
    """

    if "cellid" not in frame.columns:
        return frame.copy()

    expanded = frame.copy()
    layers: list[int | None] = []
    cells: list[int | None] = []
    for cellid in expanded["cellid"].tolist():
        if isinstance(cellid, tuple) and len(cellid) >= 2:
            layers.append(int(cellid[0]))
            cells.append(int(cellid[1]))
        elif cellid is None:
            layers.append(None)
            cells.append(None)
        else:
            layers.append(0)
            cells.append(int(cellid))
    expanded["layer"] = layers
    expanded["cell"] = cells
    return expanded.drop(columns=["cellid"])


def _normalize_term_filter(term: str | Iterable[str] | None) -> list[str] | None:
    """Normalize an optional package-budget term filter to uppercase strings."""

    if term is None:
        return None
    if isinstance(term, str):
        return [str(term).strip().upper()]
    return [str(value).strip().upper() for value in term]


def _normalize_connection_type_filter(
    connection_type: str | Iterable[str] | None,
) -> list[str] | None:
    """Normalize optional connection-type filters to uppercase labels."""

    if connection_type is None:
        return None
    if isinstance(connection_type, str):
        return [connection_type.upper()]
    normalized: list[str] = []
    for value in connection_type:
        normalized.append(str(value).upper())
    return normalized


def _normalize_surface_water_include(include: str | Iterable[str] | None) -> list[str]:
    """Normalize selected surface-water package names."""

    if include is None:
        normalized = ["sfr", "lak"]
    elif isinstance(include, str):
        normalized = [include.lower()]
    else:
        normalized = [str(value).lower() for value in include]
    allowed = {"sfr", "lak"}
    invalid = sorted(set(normalized) - allowed)
    if invalid:
        raise ValueError(
            f"Unsupported surface-water packages: {invalid!r}. Allowed values are 'sfr' and 'lak'."
        )
    ordered: list[str] = []
    for package_name in ("sfr", "lak"):
        if package_name in normalized:
            ordered.append(package_name)
    return ordered


def _infer_default_value_column(
    frame: pd.DataFrame, *, fallback: str | None = None
) -> str:
    """Infer the primary numeric value column for a normalized input table."""

    if fallback is not None and fallback in frame.columns:
        return fallback

    excluded = {"model", "package", "per", "layer", "cell", "ifno"}
    numeric = [
        column
        for column in frame.columns
        if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
    ]
    if not numeric:
        raise ValueError(
            "Could not infer a numeric value column for this package table."
        )
    return numeric[0]


def _normalize_iterable_filter(values: int | Iterable[int] | None) -> list[int] | None:
    """Normalize optional scalar-or-iterable selectors to integer lists."""

    if values is None:
        return None
    if isinstance(values, Iterable) and not isinstance(values, (str, bytes)):
        return [int(value) for value in values]
    return [int(values)]


def _filter_normalized_table(
    frame: pd.DataFrame,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
) -> pd.DataFrame:
    """Apply standard ``per/layer/cell`` filters to a normalized input table."""

    selected = frame.copy()
    if per is not None and "per" in selected.columns:
        selected = selected[selected["per"] == int(per)]
    layer_values = _normalize_iterable_filter(layer)
    if layer_values is not None and "layer" in selected.columns:
        selected = selected[selected["layer"].isin(layer_values)]
    cell_values = _normalize_iterable_filter(cells)
    if cell_values is not None and "cell" in selected.columns:
        selected = selected[selected["cell"].isin(cell_values)]
    return selected.reset_index(drop=True)


def _aggregate_hover_strings(series: pd.Series) -> str:
    """Aggregate a cell's non-numeric values into a concise hover string."""

    values = [str(value) for value in series.dropna().tolist()]
    if not values:
        return ""
    unique_values = list(dict.fromkeys(values))
    return ", ".join(unique_values)


def _default_show_layer_elevs(model) -> bool:
    """Return whether choropleths should show layer elevations by default."""

    return getattr(model.vor, "gdf_topbtm", None) is not None


__all__ = [
    "_get_package_explorer_cache",
    "_extract_structured_column",
    "_coerce_numeric_like_columns",
    "split_cellid_columns",
    "_normalize_term_filter",
    "_normalize_connection_type_filter",
    "_normalize_surface_water_include",
    "_infer_default_value_column",
    "_normalize_iterable_filter",
    "_filter_normalized_table",
    "_aggregate_hover_strings",
    "_default_show_layer_elevs",
]
