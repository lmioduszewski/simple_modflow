"""GIS loading and bounds helpers for ``myflopy`` calibration specs."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.pest.specs import VectorParameterSource


def load_vector_parameter_source(source: VectorParameterSource) -> gpd.GeoDataFrame:
    """Load and validate a vector parameter source."""

    gdf = gpd.read_file(source.path)
    required = [source.value_column]
    if source.feature_id_column is not None:
        required.append(source.feature_id_column)
    missing = [column for column in required if column not in gdf.columns]
    if missing:
        raise ValueError(
            f"Vector parameter source {Path(source.path)!s} is missing required columns: "
            + ", ".join(missing)
        )
    if gdf.geometry.is_empty.any():
        raise ValueError(f"Vector parameter source {Path(source.path)!s} contains empty geometries.")
    return gdf


def derive_bounds(
    base_values: pd.Series,
    bounds: tuple[float, float] | None,
    bounds_mode: str,
    *,
    lower_bound_column: pd.Series | None = None,
    upper_bound_column: pd.Series | None = None,
) -> tuple[pd.Series, pd.Series]:
    """Derive per-parameter lower and upper bounds.

    Parameters
    ----------
    base_values
        Base parameter values.
    bounds
        Global bounds tuple used in ``absolute`` and ``multiplier`` modes.
    bounds_mode
        One of ``"absolute"``, ``"multiplier"``, ``"from_columns"``, or
        ``"multiplier_from_columns"``.
    lower_bound_column, upper_bound_column
        Optional per-row bound columns used by the column-based modes.
    """

    base = pd.to_numeric(base_values, errors="coerce")
    mode = str(bounds_mode).strip().lower()

    if mode == "absolute":
        if bounds is None:
            raise ValueError("Absolute bounds mode requires a bounds tuple.")
        lower = pd.Series(float(bounds[0]), index=base.index)
        upper = pd.Series(float(bounds[1]), index=base.index)
    elif mode == "multiplier":
        if bounds is None:
            raise ValueError("Multiplier bounds mode requires a bounds tuple.")
        lower = base * float(bounds[0])
        upper = base * float(bounds[1])
    elif mode == "from_columns":
        if lower_bound_column is None or upper_bound_column is None:
            raise ValueError("from_columns mode requires lower_bound_column and upper_bound_column.")
        lower = pd.to_numeric(lower_bound_column, errors="coerce")
        upper = pd.to_numeric(upper_bound_column, errors="coerce")
    elif mode == "multiplier_from_columns":
        if lower_bound_column is None or upper_bound_column is None:
            raise ValueError(
                "multiplier_from_columns mode requires lower_bound_column and upper_bound_column."
            )
        lower = base * pd.to_numeric(lower_bound_column, errors="coerce")
        upper = base * pd.to_numeric(upper_bound_column, errors="coerce")
    else:
        raise ValueError(f"Unsupported bounds_mode {bounds_mode!r}.")

    invalid = lower.isna() | upper.isna()
    if invalid.any():
        raise ValueError("Derived bounds contain missing values.")
    swapped = lower > upper
    if swapped.any():
        raise ValueError("Derived bounds contain lower values greater than upper values.")
    return lower.astype(float), upper.astype(float)
