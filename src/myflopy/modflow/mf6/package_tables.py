"""Normalized package input and static table builders."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_explorer_utils import (
    _coerce_numeric_like_columns,
    _extract_structured_column,
    _get_package_explorer_cache,
    _normalize_iterable_filter,
    split_cellid_columns,
)


def build_cell_package_input_table(
    model: "SimulationBase",
    package_name: str,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
) -> pd.DataFrame:
    """Return a normalized table for a cell-based MF6 stress-period package.

    Parameters
    ----------
    model
        Live or file-backed model object.
    package_name
        MF6 package name such as ``"rch"``, ``"chd"``, ``"drn"``, or ``"ghb"``.

    Returns
    -------
    pandas.DataFrame
        Tidy DataFrame with at least ``package``, ``per``, ``layer``, and
        ``cell`` columns plus the original package data fields.
    """

    package = model.package(package_name)
    layer_values = _normalize_iterable_filter(layer)
    cell_values = _normalize_iterable_filter(cells)
    if per is None and layer_values is None and cell_values is None:
        cache_key = ("cell_package_input_table", str(package_name).lower())
        cached = _get_package_explorer_cache(model).get(cache_key)
        if cached is not None:
            return cached.copy()
    frames: list[pd.DataFrame] = []
    if per is not None:
        period_items = [(int(per), package.stress_period_data.get_data(key=int(per)))]
    else:
        period_items = package.stress_period_data.get_data().items()
    for period, data in period_items:
        frame = pd.DataFrame(data).copy()
        if frame.empty:
            continue
        frame = split_cellid_columns(frame)
        if layer_values is not None and "layer" in frame.columns:
            frame = frame[frame["layer"].isin(layer_values)]
        if cell_values is not None and "cell" in frame.columns:
            frame = frame[frame["cell"].isin(cell_values)]
        if frame.empty:
            continue
        frame["package"] = str(package_name).lower()
        frame["per"] = int(period)
        frame["model"] = model.name
        frames.append(frame)

    if not frames:
        return pd.DataFrame(columns=["model", "package", "per", "layer", "cell"])

    combined = pd.concat(frames, ignore_index=True)
    for column in ("layer", "cell", "per"):
        if column in combined.columns:
            combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(combined, exclude={"model", "package"})
    if per is None and layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[
            ("cell_package_input_table", str(package_name).lower())
        ] = combined.copy()
    return combined


def build_uzf_field_input_table(
    model: "SimulationBase",
    field_name: str,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
) -> pd.DataFrame:
    """Return a normalized table for one UZF perioddata field.

    Parameters
    ----------
    model
        Live or file-backed model object.
    field_name
        Perioddata field name such as ``"finf"``.

    Returns
    -------
    pandas.DataFrame
        Tidy DataFrame with ``per``, ``ifno``, zero-based ``layer``/``cell``,
        and the requested field.
    """

    package = model.package("uzf")
    layer_values = _normalize_iterable_filter(layer)
    cell_values = _normalize_iterable_filter(cells)
    if per is None and layer_values is None and cell_values is None:
        cache_key = ("uzf_field_input_table", str(field_name).lower())
        cached = _get_package_explorer_cache(model).get(cache_key)
        if cached is not None:
            return cached.copy()
    ifno_to_cellid = model.outputs.uzf.ifno_to_cellid
    layer_lookup = ifno_to_cellid["layer"].to_numpy()
    cell_lookup = ifno_to_cellid["cellid"].to_numpy()
    frames: list[pd.DataFrame] = []
    if per is not None:
        period_items = [(int(per), package.perioddata.get_data(key=int(per)))]
    else:
        period_items = package.perioddata.get_data().items()
    for period, data in period_items:
        if data is None or len(data) == 0:
            continue
        try:
            ifno_raw = _extract_structured_column(data, "ifno")
            value_raw = _extract_structured_column(data, field_name)
        except KeyError as exc:
            raise KeyError(f"UZF field {field_name!r} was not found.") from exc

        ifno_series = pd.to_numeric(pd.Series(ifno_raw), errors="coerce")
        valid_mask = ifno_series.notna().to_numpy()
        if not valid_mask.any():
            continue
        ifno_values = ifno_series.loc[valid_mask].astype(int).to_numpy()
        layer_values_for_rows = np.take(layer_lookup, ifno_values)
        cell_values_for_rows = np.take(cell_lookup, ifno_values)
        field_values = np.asarray(value_raw)[valid_mask]

        row_mask = np.ones(len(ifno_values), dtype=bool)
        if layer_values is not None:
            row_mask &= np.isin(layer_values_for_rows, layer_values)
        if cell_values is not None:
            row_mask &= np.isin(cell_values_for_rows, cell_values)
        if not row_mask.any():
            continue

        ifno_values = ifno_values[row_mask]
        layer_values_for_rows = layer_values_for_rows[row_mask]
        cell_values_for_rows = cell_values_for_rows[row_mask]
        field_values = field_values[row_mask]
        frames.append(
            pd.DataFrame(
                {
                    "model": model.name,
                    "package": "uzf",
                    "per": int(period),
                    "ifno": ifno_values,
                    "layer": layer_values_for_rows,
                    "cell": cell_values_for_rows,
                    field_name: field_values,
                }
            )
        )

    if not frames:
        return pd.DataFrame(
            columns=["model", "package", "per", "ifno", "layer", "cell", field_name]
        )

    combined = pd.concat(frames, ignore_index=True)
    for column in ("per", "ifno", "layer", "cell"):
        combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(combined, exclude={"model", "package"})
    if per is None and layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[
            ("uzf_field_input_table", str(field_name).lower())
        ] = combined.copy()
    return combined


def build_uzf_field_input_wide_table(
    model: "SimulationBase",
    field_name: str,
    *,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
) -> pd.DataFrame:
    """Return a compact wide UZF table with one row per UZF record/cell."""

    cache_key = ("uzf_field_input_wide_table", str(field_name).lower())
    layer_values = _normalize_iterable_filter(layer)
    cell_values = _normalize_iterable_filter(cells)
    if layer_values is None and cell_values is None:
        cached = _get_package_explorer_cache(model).get(cache_key)
        if cached is not None:
            return cached.copy()

    package = model.package("uzf")
    ifno_to_cellid = model.outputs.uzf.ifno_to_cellid.copy()
    ifno_to_cellid.index.name = "ifno"
    base = ifno_to_cellid.reset_index().rename(columns={"cellid": "cell"})
    if layer_values is not None:
        base = base[base["layer"].isin(layer_values)]
    if cell_values is not None:
        base = base[base["cell"].isin(cell_values)]
    if base.empty:
        return pd.DataFrame(columns=["ifno", "layer", "cell"])

    wide = base.set_index("ifno").copy()
    period_items = package.perioddata.get_data().items()
    for period, data in period_items:
        if data is None or len(data) == 0:
            continue
        try:
            ifno_raw = _extract_structured_column(data, "ifno")
            value_raw = _extract_structured_column(data, field_name)
        except KeyError as exc:
            raise KeyError(f"UZF field {field_name!r} was not found.") from exc
        ifno_series = pd.to_numeric(pd.Series(ifno_raw), errors="coerce")
        value_series = pd.Series(value_raw)
        valid_mask = ifno_series.notna().to_numpy()
        if not valid_mask.any():
            continue
        period_series = pd.Series(
            value_series.loc[valid_mask].to_numpy(),
            index=ifno_series.loc[valid_mask].astype(int).to_numpy(),
            name=f"per_{int(period)}",
        )
        wide = wide.join(period_series, how="left")

    wide = wide.reset_index()
    wide = _coerce_numeric_like_columns(wide)
    if layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[cache_key] = wide.copy()
    return wide


def build_sfr_reach_table(model: "SimulationBase") -> pd.DataFrame:
    """Return per-reach geometry metadata used by SFR tables and profiles."""

    packagedata = pd.DataFrame(model.sfr.packagedata.get_data()).copy()
    packagedata = split_cellid_columns(packagedata)
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "reach"})
    packagedata["reach"] = packagedata["reach"].astype(int)
    packagedata["layer"] = packagedata["layer"].astype(int)
    packagedata["cell"] = packagedata["cell"].astype(int)
    packagedata["rlen"] = pd.to_numeric(packagedata["rlen"], errors="coerce")
    packagedata = packagedata.sort_values("reach").reset_index(drop=True)
    packagedata["distance_start"] = packagedata["rlen"].fillna(
        0.0
    ).cumsum() - packagedata["rlen"].fillna(0.0)
    packagedata["distance_end"] = packagedata["distance_start"] + packagedata[
        "rlen"
    ].fillna(0.0)
    packagedata["distance_mid"] = packagedata["distance_start"] + (
        packagedata["rlen"].fillna(0.0) / 2.0
    )
    return packagedata[
        [
            "reach",
            "layer",
            "cell",
            "rlen",
            "distance_start",
            "distance_mid",
            "distance_end",
        ]
    ].drop_duplicates()


def build_lak_connection_table(model: "SimulationBase") -> pd.DataFrame:
    """Return a normalized table of LAK connection geometry by cell.

    Notes
    -----
    The returned table keeps one row per LAK connection record and computes a
    ``connection_area`` field that can be mapped directly as a cell
    choropleth:

    - vertical connections use the plan-view Voronoi cell area
    - horizontal/embedded horizontal connections use ``connwidth * (telev - belev)``

    The values therefore represent the physical interface area available for
    exchange between the lake and groundwater.
    """

    connectiondata = pd.DataFrame(model.lak.connectiondata.get_data()).copy()
    if connectiondata.empty:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "lake",
                "iconn",
                "layer",
                "cell",
                "claktype",
                "belev",
                "telev",
                "connlen",
                "connwidth",
                "connection_area",
            ]
        )

    connectiondata = split_cellid_columns(connectiondata)
    if "ifno" in connectiondata.columns:
        connectiondata = connectiondata.rename(columns={"ifno": "lake"})
    connectiondata["lake"] = pd.to_numeric(
        connectiondata["lake"], errors="coerce"
    ).astype("Int64")
    if "iconn" in connectiondata.columns:
        connectiondata["iconn"] = pd.to_numeric(
            connectiondata["iconn"], errors="coerce"
        ).astype("Int64")
    connectiondata["layer"] = pd.to_numeric(
        connectiondata["layer"], errors="coerce"
    ).astype(int)
    connectiondata["cell"] = pd.to_numeric(
        connectiondata["cell"], errors="coerce"
    ).astype(int)
    connectiondata["package"] = "lak"
    connectiondata["model"] = model.name
    if "claktype" in connectiondata.columns:
        connectiondata["claktype"] = (
            connectiondata["claktype"].astype("string").str.upper()
        )
    else:
        connectiondata["claktype"] = "UNKNOWN"

    cell_area_lookup = {
        int(cell): float(area)
        for cell, area in enumerate(
            np.asarray(model.vor.area_list, dtype=float).reshape(-1)
        )
    }
    top = pd.to_numeric(connectiondata.get("telev"), errors="coerce")
    bottom = pd.to_numeric(connectiondata.get("belev"), errors="coerce")
    width = pd.to_numeric(connectiondata.get("connwidth"), errors="coerce")
    vertical_mask = connectiondata["claktype"].astype(str).str.upper().eq("VERTICAL")
    thickness = (top - bottom).clip(lower=0.0)
    horizontal_area = width * thickness
    vertical_area = connectiondata["cell"].map(cell_area_lookup).astype(float)
    connectiondata["connection_area"] = np.where(
        vertical_mask, vertical_area, horizontal_area
    )

    ordered = [
        "model",
        "package",
        "lake",
        "iconn",
        "layer",
        "cell",
        "claktype",
        "belev",
        "telev",
        "connlen",
        "connwidth",
        "connection_area",
    ]
    remaining = [column for column in connectiondata.columns if column not in ordered]
    return connectiondata[ordered + remaining]


def build_sfr_input_table(model: "SimulationBase") -> pd.DataFrame:
    """Return normalized static and stress-period SFR inputs mapped to cells."""

    base = pd.DataFrame(model.sfr.packagedata.get_data()).copy()
    base = split_cellid_columns(base)
    if "ifno" in base.columns:
        base = base.rename(columns={"ifno": "reach"})
    base["reach"] = pd.to_numeric(base["reach"], errors="coerce").astype(int)
    base["layer"] = pd.to_numeric(base["layer"], errors="coerce").astype(int)
    base["cell"] = pd.to_numeric(base["cell"], errors="coerce").astype(int)
    settings = _numeric_period_settings(model.sfr, id_column="ifno", nper=model.nper)
    frames = []
    for period in range(int(model.nper)):
        frame = base.copy()
        for reach, values in settings[period].items():
            mask = frame["reach"] == int(reach)
            for name, value in values.items():
                frame.loc[mask, name] = value
        frame["per"] = period
        frame["package"] = "sfr"
        frame["model"] = model.name
        frames.append(frame)
    return _coerce_numeric_like_columns(
        pd.concat(frames, ignore_index=True), exclude={"model", "package"}
    )


def build_lak_input_table(model: "SimulationBase") -> pd.DataFrame:
    """Return normalized static and stress-period LAK inputs mapped to connection cells."""

    connections = build_lak_connection_table(model)
    packagedata = pd.DataFrame(model.lak.packagedata.get_data()).copy()
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "lake"})
    packagedata["lake"] = pd.to_numeric(packagedata["lake"], errors="coerce").astype(
        int
    )
    base = connections.merge(
        packagedata, on="lake", how="left", suffixes=("", "_package")
    )
    settings = _numeric_period_settings(model.lak, id_column="number", nper=model.nper)
    frames = []
    for period in range(int(model.nper)):
        frame = base.copy()
        for lake, values in settings[period].items():
            mask = frame["lake"] == int(lake)
            for name, value in values.items():
                frame.loc[mask, name] = value
        frame["per"] = period
        frame["package"] = "lak"
        frame["model"] = model.name
        frames.append(frame)
    return _coerce_numeric_like_columns(
        pd.concat(frames, ignore_index=True), exclude={"model", "package"}
    )


def summarize_input_table(
    frame: pd.DataFrame, *, label: str, value_columns: list[str]
) -> pd.DataFrame:
    """Build a compact one-row summary for a normalized input table."""

    if frame.empty:
        return pd.DataFrame(
            [
                {
                    "label": label,
                    "records": 0,
                    "periods": 0,
                    "layers": 0,
                    "cells": 0,
                    "value_columns": value_columns,
                }
            ]
        )

    return pd.DataFrame(
        [
            {
                "label": label,
                "records": int(len(frame)),
                "periods": int(frame["per"].nunique()) if "per" in frame.columns else 0,
                "layers": int(frame["layer"].nunique())
                if "layer" in frame.columns
                else 0,
                "cells": int(frame["cell"].nunique()) if "cell" in frame.columns else 0,
                "value_columns": value_columns,
            }
        ]
    )


def _numeric_period_settings(
    package, *, id_column: str, nper: int
) -> dict[int, dict[int, dict[str, float]]]:
    """Return carried-forward numeric package settings keyed by period and feature id."""

    current: dict[int, dict[str, float]] = {}
    by_period: dict[int, dict[int, dict[str, float]]] = {}
    perioddata = package.perioddata.get_data()
    for period in range(int(nper)):
        records = perioddata.get(period)
        if records is not None:
            frame = pd.DataFrame(records).copy()
            if not frame.empty:
                identifier = (
                    id_column if id_column in frame.columns else frame.columns[0]
                )
                setting_column = next(
                    (
                        column
                        for column in frame.columns
                        if str(column).lower().endswith("setting")
                    ),
                    None,
                )
                data_column = next(
                    (
                        column
                        for column in frame.columns
                        if str(column).lower().endswith("setting_data")
                    ),
                    None,
                )
                if setting_column is not None and data_column is not None:
                    for row in frame.itertuples(index=False):
                        feature_id = int(getattr(row, identifier))
                        setting = (
                            str(getattr(row, setting_column))
                            .strip()
                            .lower()
                            .replace("-", "_")
                        )
                        raw_value = getattr(row, data_column)
                        numeric = pd.to_numeric(
                            pd.Series([raw_value]), errors="coerce"
                        ).iloc[0]
                        if pd.notna(numeric):
                            current.setdefault(feature_id, {})[setting] = float(numeric)
        by_period[period] = {
            feature_id: dict(settings) for feature_id, settings in current.items()
        }
    return by_period


__all__ = [
    "build_cell_package_input_table",
    "build_uzf_field_input_table",
    "build_uzf_field_input_wide_table",
    "build_sfr_reach_table",
    "build_lak_connection_table",
    "build_sfr_input_table",
    "build_lak_input_table",
    "summarize_input_table",
    "_numeric_period_settings",
]
