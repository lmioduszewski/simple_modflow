"""Preferred package-input exploration API for MF6 models and loaded runs.

This module provides the normalized table and choropleth helpers behind the
``model.packages...`` surface. The goal is to make package inputs easier to
inspect consistently across live :class:`SimulationBase` models and lazy
file-backed :class:`LoadedMf6Run` objects without duplicating package-specific
logic in many places.

The first slice focuses on:

- cell-based stress-period packages: ``RCH``, ``CHD``, ``DRN``, ``GHB``
- registry-backed stress-period packages such as ``WEL`` and ``UZF``
- static layer/cell arrays such as ``IC``, ``NPF``, and ``STO``

Each explorer follows the same shape where practical:

- ``get()`` returns a normalized tidy DataFrame
- ``summary()`` returns a compact one-row summary
- ``map()`` returns the standard choropleth wrapper used elsewhere in the
  package

The first results slice adds:

- budget-backed cell results for ``RCH``, ``CHD``, ``DRN``, ``GHB``
- budget-backed groundwater recharge results for ``UZF`` via ``UZF-GWRCH``
- stage-result explorers for ``LAK`` and ``SFR``
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import figs
import plotly.graph_objects as go
from flopy.plot import PlotMapView
from plotly.subplots import make_subplots

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


@dataclass(frozen=True)
class FieldSpec:
    """Registry metadata for one package input field."""

    name: str
    label: str | None = None
    colorscale: str | None = None
    fill_value: float = 0.0
    agg: str = "sum"


@dataclass(frozen=True)
class ResultSpec:
    """Registry metadata for one package result term."""

    name: str
    budget_text: str
    value_name: str = "q"
    label: str | None = None
    colorscale: str | None = None
    diverging: bool = True


@dataclass(frozen=True)
class PackageExplorerSpec:
    """Registry metadata for one MF6 package explorer."""

    name: str
    kind: str = "cell_stress"
    default_input: str | None = None
    colorscale: str | None = None
    inputs: dict[str, FieldSpec] = field(default_factory=dict)
    results: dict[str, ResultSpec] = field(default_factory=dict)


_PACKAGE_EXPLORER_SPECS: dict[str, PackageExplorerSpec] = {
    "rch": PackageExplorerSpec(
        name="rch",
        default_input="recharge",
        colorscale="Viridis",
        inputs={
            "recharge": FieldSpec("recharge", label="Recharge", colorscale="Viridis"),
        },
        results={
            "q": ResultSpec("q", budget_text="RCH", value_name="q", colorscale="RdBu"),
        },
    ),
    "chd": PackageExplorerSpec(
        name="chd",
        default_input="head",
        colorscale="Blues",
        inputs={
            "head": FieldSpec("head", label="Constant head", colorscale="Blues"),
        },
        results={
            "q": ResultSpec("q", budget_text="CHD", value_name="q", colorscale="RdBu"),
        },
    ),
    "drn": PackageExplorerSpec(
        name="drn",
        default_input="elev",
        colorscale="YlOrRd",
        inputs={
            "elev": FieldSpec("elev", label="Drain elevation", colorscale="YlOrRd"),
            "cond": FieldSpec("cond", label="Drain conductance", colorscale="Viridis"),
        },
        results={
            "q": ResultSpec("q", budget_text="DRN", value_name="q", colorscale="RdBu"),
        },
    ),
    "ghb": PackageExplorerSpec(
        name="ghb",
        default_input="bhead",
        colorscale="Portland",
        inputs={
            "bhead": FieldSpec("bhead", label="Boundary head", colorscale="Portland"),
            "cond": FieldSpec("cond", label="Boundary conductance", colorscale="Viridis"),
        },
        results={
            "q": ResultSpec("q", budget_text="GHB", value_name="q", colorscale="RdBu"),
        },
    ),
    "wel": PackageExplorerSpec(
        name="wel",
        default_input="q",
        colorscale="RdBu",
        inputs={
            "q": FieldSpec("q", label="Well flow", colorscale="RdBu"),
        },
        results={
            "q": ResultSpec("q", budget_text="WEL", value_name="q", colorscale="RdBu"),
        },
    ),
    "uzf": PackageExplorerSpec(
        name="uzf",
        kind="uzf",
        inputs={
            "finf": FieldSpec("finf", label="UZF infiltration", colorscale="Viridis"),
            "pet": FieldSpec("pet", label="Potential evapotranspiration", colorscale="YlOrRd"),
            "extdp": FieldSpec("extdp", label="ET extinction depth", colorscale="Blues"),
            "extwc": FieldSpec("extwc", label="ET extinction water content", colorscale="Viridis"),
            "ha": FieldSpec("ha", label="Surface depression storage depth", colorscale="Blues"),
            "hroot": FieldSpec("hroot", label="Root zone thickness", colorscale="Blues"),
            "rootact": FieldSpec("rootact", label="Root activity", colorscale="Viridis"),
        },
        results={
            "gwrch": ResultSpec(
                "gwrch",
                budget_text="UZF-GWRCH",
                value_name="gwrch",
                label="UZF groundwater recharge",
                colorscale="Viridis",
                diverging=False,
            ),
            "sat": ResultSpec(
                "sat",
                budget_text="DATA-SAT",
                value_name="sat",
                label="UZF saturation",
                colorscale="Viridis",
                diverging=False,
            ),
        },
    ),
    "sfr": PackageExplorerSpec(
        name="sfr",
        kind="surface_water",
        results={
            "q": ResultSpec("q", budget_text="SFR", value_name="q", colorscale="RdBu"),
        },
    ),
    "lak": PackageExplorerSpec(
        name="lak",
        kind="surface_water",
        results={
            "q": ResultSpec("q", budget_text="GWF", value_name="q", colorscale="RdBu"),
        },
    ),
}

_GROUP_COMPARE_DEFAULT_COLORSCALE = "RdBu"


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


def _coerce_numeric_like_columns(frame: pd.DataFrame, *, exclude: Iterable[str] = ()) -> pd.DataFrame:
    """Convert object/string columns to numeric when every non-null value is numeric-like."""

    excluded = {str(column) for column in exclude}
    coerced = frame.copy()
    for column in coerced.columns:
        if str(column) in excluded:
            continue
        series = coerced[column]
        if not (pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series)):
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
        _get_package_explorer_cache(model)[("cell_package_input_table", str(package_name).lower())] = combined.copy()
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
        return pd.DataFrame(columns=["model", "package", "per", "ifno", "layer", "cell", field_name])

    combined = pd.concat(frames, ignore_index=True)
    for column in ("per", "ifno", "layer", "cell"):
        combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(combined, exclude={"model", "package"})
    if per is None and layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[("uzf_field_input_table", str(field_name).lower())] = combined.copy()
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


def _normalize_budget_nodes(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize budget node columns to zero-based indexing."""

    normalized = frame.copy()
    for column in ("node", "node2"):
        if column not in normalized.columns:
            continue
        numeric = pd.to_numeric(normalized[column], errors="coerce")
        mask = numeric.notna()
        if not mask.any():
            continue
        ints = numeric.loc[mask].astype(int)
        normalized[column] = numeric
        if int(ints.min()) >= 1:
            normalized.loc[mask, column] = (ints - 1).astype(float)
        else:
            normalized.loc[mask, column] = ints.astype(float)
    return normalized


def build_budget_result_table(
    model: "SimulationBase",
    *,
    budget_text: str,
    package_name: str,
    value_name: str = "q",
) -> pd.DataFrame:
    """Return a normalized cell-result table from one budget record type.

    Parameters
    ----------
    model
        Live or file-backed model object.
    budget_text
        Binary budget record name such as ``"RCH"`` or ``"UZF-GWRCH"``.
    package_name
        Logical package name reported in the normalized table.
    value_name
        Public column name to use for the normalized budget quantity.

    Returns
    -------
    pandas.DataFrame
        Tidy DataFrame with ``kstpkper``, zero-based ``per``, ``layer``,
        ``cell``, and the normalized result field.
    """

    reader = model._get_budget_reader()
    frames: list[pd.DataFrame] = []
    for kstpkper in model._get_budget_kstpkper():
        records = reader.get_data(text=budget_text, kstpkper=kstpkper)
        if not records:
            continue
        frame = pd.DataFrame.from_records(records[0]).copy()
        if frame.empty or "node" not in frame.columns:
            continue
        frame = _normalize_budget_nodes(frame)
        nodes = frame["node"].astype(int).tolist()
        layers: list[int] = []
        cells: list[int] = []
        for node in nodes:
            layer, cell = model.node_to_lni[int(node)]
            layers.append(int(layer))
            cells.append(int(cell))
        frame["layer"] = layers
        frame["cell"] = cells
        frame["model"] = model.name
        frame["package"] = str(package_name).lower()
        frame["kstpkper"] = [tuple(int(value) for value in kstpkper)] * len(frame)
        frame["per"] = int(kstpkper[1])
        if value_name != "q" and "q" in frame.columns:
            frame[value_name] = pd.to_numeric(frame["q"], errors="coerce")
        frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=["model", "package", "kstpkper", "per", "layer", "cell", value_name]
        )

    combined = pd.concat(frames, ignore_index=True)
    for column in ("per", "layer", "cell"):
        combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(
        combined,
        exclude={"model", "package", "kstpkper"},
    )
    return combined


def build_sfr_stage_result_table(model: "SimulationBase") -> pd.DataFrame:
    """Return a normalized table of SFR stage results by reach and period."""

    stage_frame = model.outputs.sfr.stage.get().copy()
    stage_frame.index.name = "reach"
    stage_long = stage_frame.reset_index().melt(
        id_vars="reach", var_name="kstpkper", value_name="stage"
    )
    stage_long["reach"] = stage_long["reach"].astype(int)
    stage_long["kstpkper"] = stage_long["kstpkper"].apply(tuple)
    stage_long["per"] = stage_long["kstpkper"].apply(lambda value: int(value[1]))

    reach_table = build_sfr_reach_table(model)
    stage_long = stage_long.merge(reach_table, on="reach", how="left")
    stage_long["model"] = model.name
    stage_long["package"] = "sfr"
    stage_long["layer"] = stage_long["layer"].astype(int)
    stage_long["cell"] = stage_long["cell"].astype(int)
    ordered = [
        "model",
        "package",
        "kstpkper",
        "per",
        "reach",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        "stage",
    ]
    return stage_long[ordered]


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
    packagedata["distance_start"] = packagedata["rlen"].fillna(0.0).cumsum() - packagedata["rlen"].fillna(0.0)
    packagedata["distance_end"] = packagedata["distance_start"] + packagedata["rlen"].fillna(0.0)
    packagedata["distance_mid"] = packagedata["distance_start"] + (packagedata["rlen"].fillna(0.0) / 2.0)
    return packagedata[
        ["reach", "layer", "cell", "rlen", "distance_start", "distance_mid", "distance_end"]
    ].drop_duplicates()


def build_sfr_long_profile_table(model: "SimulationBase", *, per: int = 0) -> pd.DataFrame:
    """Return one merged SFR long-profile table for a selected stress period.

    Parameters
    ----------
    model
        Live or file-backed model object.
    per
        Zero-based stress period to extract.

    Returns
    -------
    pandas.DataFrame
        Reach-ordered table combining geometry metadata, packagedata fields,
        stage results, and stream-groundwater exchange results. The table is
        intended as the canonical input to higher-level SFR profile plots.
    """

    packagedata = pd.DataFrame(model.sfr.packagedata.get_data()).copy()
    packagedata = split_cellid_columns(packagedata)
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "reach"})
    packagedata["reach"] = packagedata["reach"].astype(int)
    packagedata["layer"] = packagedata["layer"].astype(int)
    packagedata["cell"] = packagedata["cell"].astype(int)
    for column in ("rlen", "rwid", "rgrd", "rtp", "rbth", "rhk", "man"):
        if column in packagedata.columns:
            packagedata[column] = pd.to_numeric(packagedata[column], errors="coerce")
    packagedata["streambed_top"] = pd.to_numeric(packagedata.get("rtp"), errors="coerce")
    packagedata["streambed_bottom"] = (
        pd.to_numeric(packagedata.get("rtp"), errors="coerce")
        - pd.to_numeric(packagedata.get("rbth"), errors="coerce")
    )

    reach_table = build_sfr_reach_table(model)
    packagedata = packagedata.merge(
        reach_table,
        on=["reach", "layer", "cell", "rlen"],
        how="left",
    )

    stage = build_sfr_stage_result_table(model)
    stage = (
        stage.loc[stage["per"] == int(per)]
        .sort_values("kstpkper")
        .groupby("reach", as_index=False)
        .tail(1)
        .loc[:, ["reach", "stage"]]
        .drop_duplicates(subset=["reach"])
        .reset_index(drop=True)
    )

    q_frame = build_sfr_budget_result_table(model, budget_text="SFR", value_name="q")
    if not q_frame.empty:
        q_frame = (
            q_frame.loc[q_frame["per"] == int(per), ["reach", "q"]]
            .dropna(subset=["reach"])
            .groupby("reach", as_index=False)["q"]
            .sum()
        )
        q_frame["reach"] = q_frame["reach"].astype(int)
    else:
        q_frame = pd.DataFrame(columns=["reach", "q"])

    profile = packagedata.merge(stage, on="reach", how="left").merge(q_frame, on="reach", how="left")
    profile["q_per_length"] = np.where(
        pd.to_numeric(profile["rlen"], errors="coerce") > 0.0,
        pd.to_numeric(profile["q"], errors="coerce") / pd.to_numeric(profile["rlen"], errors="coerce"),
        np.nan,
    )
    profile["model"] = model.name
    profile["package"] = "sfr"
    profile["per"] = int(per)
    ordered = [
        "model",
        "package",
        "per",
        "reach",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        "rwid",
        "rgrd",
        "rtp",
        "rbth",
        "rhk",
        "man",
        "streambed_top",
        "streambed_bottom",
        "stage",
        "q",
        "q_per_length",
    ]
    remaining = [column for column in profile.columns if column not in ordered]
    return profile[ordered + remaining].sort_values(["reach", "cell"]).reset_index(drop=True)


def build_sfr_budget_result_table(
    model: "SimulationBase",
    *,
    budget_text: str = "SFR",
    value_name: str = "q",
) -> pd.DataFrame:
    """Return a normalized SFR exchange-result table with reach ids attached."""

    frame = build_budget_result_table(
        model,
        budget_text=budget_text,
        package_name="sfr",
        value_name=value_name,
    )
    if frame.empty:
        return pd.DataFrame(
            columns=["model", "package", "kstpkper", "per", "reach", "layer", "cell", value_name]
        )

    reach_map = build_sfr_reach_table(model)
    frame = frame.merge(reach_map, on=["layer", "cell"], how="left", suffixes=("", "_pkg"))
    frame["reach"] = frame["reach"].astype("Int64")
    if value_name == "q":
        q_series = pd.to_numeric(frame["q"], errors="coerce")
        rlen_series = pd.to_numeric(frame["rlen"], errors="coerce")
        frame["q_per_length"] = np.where(rlen_series > 0.0, q_series / rlen_series, np.nan)
    ordered = [
        "model",
        "package",
        "kstpkper",
        "per",
        "reach",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        value_name,
    ]
    if value_name == "q":
        ordered.append("q_per_length")
    remaining = [column for column in frame.columns if column not in ordered]
    return frame[ordered + remaining]


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
        "Package": grouped["package"].agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size().reindex(full_index, fill_value=0).astype(int).tolist(),
        "q_per_length": values.tolist(),
        "q": q_sum.fillna(0.0).astype(float).tolist(),
        "rlen": rlen_sum.fillna(0.0).astype(float).tolist(),
    }
    if "reach" in selected.columns:
        reach_strings = selected.assign(reach=selected["reach"].astype("string")).groupby("cell", dropna=False)["reach"]
        hover["reach"] = reach_strings.agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
    return values.tolist(), hover


def build_lak_stage_result_table(model: "SimulationBase") -> pd.DataFrame:
    """Return a normalized table of LAK stage results mapped to connected cells."""

    stage_array = np.asarray(model.outputs.lak.stage.get(), dtype=float)
    nlakes = int(model.lak.nlakes.data)
    periods = pd.DataFrame(stage_array.reshape(-1, nlakes))
    periods["per"] = periods.index.astype(int)
    stage_long = periods.melt(id_vars="per", var_name="lake", value_name="stage")
    stage_long["lake"] = stage_long["lake"].astype(int)

    connectiondata = pd.DataFrame(model.lak.connectiondata.get_data()).copy()
    connectiondata = split_cellid_columns(connectiondata)
    if "ifno" in connectiondata.columns:
        connectiondata = connectiondata.rename(columns={"ifno": "lake"})
    cell_map = connectiondata.loc[:, ["lake", "layer", "cell"]].drop_duplicates().copy()
    stage_long = stage_long.merge(cell_map, on="lake", how="left")
    stage_long["model"] = model.name
    stage_long["package"] = "lak"
    stage_long["layer"] = stage_long["layer"].astype(int)
    stage_long["cell"] = stage_long["cell"].astype(int)
    return stage_long[["model", "package", "per", "lake", "layer", "cell", "stage"]]


def build_lak_budget_result_table(
    model: "SimulationBase",
    *,
    budget_text: str = "GWF",
    value_name: str = "q",
) -> pd.DataFrame:
    """Return a normalized LAK exchange-result table with connection geometry.

    Notes
    -----
    The returned table carries one row per LAK budget record together with the
    matching LAK connection metadata. Per the MF6 I/O specification, the LAK
    package-output ``GWF`` budget rows do not include the input-file ``iconn``
    identifier. They only carry the lake id, the connected GWF cell/node id, the
    simulated flow ``q``, and the auxiliary ``FLOW-AREA`` value. The helper
    therefore reconstructs connection-specific metadata by matching repeated
    LAK budget rows back to ``lak.connectiondata`` in the written within-cell
    connection order for each stress period. When ``value_name == "q"``, the
    helper also computes:

    - ``flow_area``: the physical lake-groundwater exchange area
    - ``q_per_area``: ``q / flow_area`` with units of length per time

    For vertical lake connections, ``flow_area`` is the plan-view Voronoi cell
    area. For horizontal lake connections, it is ``connwidth * (telev - belev)``
    after clipping negative thicknesses to zero.
    """

    del budget_text
    frame = model.outputs.lak.bud.get("GWF").copy()
    if frame.empty:
        return pd.DataFrame(
            columns=["model", "package", "kstpkper", "per", "lake", "layer", "cell", value_name]
        )
    frame = _normalize_budget_nodes(frame)
    frame["lake"] = pd.to_numeric(frame["node"], errors="coerce").astype("Int64")
    node2 = pd.to_numeric(frame["node2"], errors="coerce")
    layers: list[int] = []
    cells: list[int] = []
    for node in node2.astype(int).tolist():
        layer_id, cell_id = model.node_to_lni[int(node)]
        layers.append(int(layer_id))
        cells.append(int(cell_id))
    frame["layer"] = layers
    frame["cell"] = cells
    frame["model"] = model.name
    frame["package"] = "lak"
    if "kstpkper" in frame.columns:
        frame["kstpkper"] = frame["kstpkper"].apply(lambda values: tuple(int(value) for value in values))
        frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
    else:
        frame["kstpkper"] = [(0, 0)] * len(frame)
        frame["per"] = 0
    if value_name != "q" and "q" in frame.columns:
        frame[value_name] = pd.to_numeric(frame["q"], errors="coerce")

    connectiondata = pd.DataFrame(model.lak.connectiondata.get_data()).copy()
    connectiondata = split_cellid_columns(connectiondata)
    if "ifno" in connectiondata.columns:
        connectiondata = connectiondata.rename(columns={"ifno": "lake"})
    if "iconn" in connectiondata.columns:
        connectiondata["iconn"] = pd.to_numeric(connectiondata["iconn"], errors="coerce").astype("Int64")
    connectiondata["lake"] = pd.to_numeric(connectiondata["lake"], errors="coerce").astype("Int64")
    connectiondata["layer"] = pd.to_numeric(connectiondata["layer"], errors="coerce").astype(int)
    connectiondata["cell"] = pd.to_numeric(connectiondata["cell"], errors="coerce").astype(int)
    if "claktype" in connectiondata.columns:
        connectiondata["claktype"] = connectiondata["claktype"].astype("string").str.upper()
    raw_iconn_column = next((column for column in frame.columns if str(column).lower() == "iconn"), None)
    use_raw_iconn = False
    if raw_iconn_column is not None and "iconn" in connectiondata.columns:
        if raw_iconn_column != "iconn":
            frame = frame.rename(columns={raw_iconn_column: "iconn"})
        frame["iconn"] = pd.to_numeric(frame["iconn"], errors="coerce").astype("Int64")
        grouped_iconn = frame.groupby(["lake", "layer", "cell"], dropna=False)["iconn"]
        ambiguous_groups = (grouped_iconn.size() > 1) & (grouped_iconn.nunique(dropna=False) <= 1)
        use_raw_iconn = not bool(ambiguous_groups.any())

    if use_raw_iconn:
        frame = frame.merge(
            connectiondata,
            on=["lake", "iconn"],
            how="left",
            suffixes=("", "_conn"),
        )
        for column in ("layer", "cell"):
            connection_column = f"{column}_conn"
            if connection_column in frame.columns:
                frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("Int64").combine_first(
                    pd.to_numeric(frame[connection_column], errors="coerce").astype("Int64")
                )
                frame = frame.drop(columns=[connection_column])
        if "lake_conn" in frame.columns:
            frame = frame.drop(columns=["lake_conn"])
    else:
        connectiondata["_connection_row_order"] = connectiondata.groupby(["lake", "layer", "cell"]).cumcount()
        frame_group_keys = ["lake", "layer", "cell"]
        if "kstpkper" in frame.columns:
            # LAK package-output budgets repeat the same lake-cell connection rows
            # every time step, so the within-cell row order must reset each period.
            frame_group_keys = ["kstpkper", *frame_group_keys]
        frame["_connection_row_order"] = frame.groupby(frame_group_keys).cumcount()
        frame = frame.merge(
            connectiondata,
            on=["lake", "layer", "cell", "_connection_row_order"],
            how="left",
            suffixes=("", "_conn"),
        )
        frame = frame.drop(columns=["_connection_row_order"], errors="ignore")

    fallback_columns = [
        column
        for column in ("iconn", "claktype", "bedleak", "belev", "telev", "connlen", "connwidth")
        if column in connectiondata.columns
    ]
    if fallback_columns:
        fallback_source = connectiondata.loc[:, ["lake", "layer", "cell", *fallback_columns]].copy()

        def _first_non_null(series: pd.Series):
            values = series.dropna()
            if values.empty:
                return np.nan
            return values.iloc[0]

        fallback_agg: dict[str, object] = {}
        for column in fallback_columns:
            if column == "claktype":
                fallback_agg[column] = _aggregate_hover_strings
            else:
                fallback_agg[column] = _first_non_null
        fallback = (
            fallback_source.groupby(["lake", "layer", "cell"], dropna=False, as_index=False)
            .agg(fallback_agg)
            .rename(columns={column: f"{column}_fallback" for column in fallback_columns})
        )
        frame = frame.merge(fallback, on=["lake", "layer", "cell"], how="left")
        for column in fallback_columns:
            fallback_column = f"{column}_fallback"
            if column in frame.columns and fallback_column in frame.columns:
                missing_mask = frame[column].isna()
                frame.loc[missing_mask, column] = frame.loc[missing_mask, fallback_column]
                frame = frame.drop(columns=[fallback_column])

    if value_name == "q":
        flow_area = pd.to_numeric(frame.get("FLOW-AREA"), errors="coerce")
        q_series = pd.to_numeric(frame["q"], errors="coerce")
        frame["flow_area"] = flow_area
        frame["q_per_area"] = np.where(flow_area > 0.0, q_series / flow_area, np.nan)

    ordered = [
        "model",
        "package",
        "kstpkper",
        "per",
        "lake",
        "layer",
        "cell",
        "claktype",
        "belev",
        "telev",
        "connlen",
        "connwidth",
        value_name,
    ]
    if value_name == "q":
        ordered.extend(["flow_area", "q_per_area"])
    remaining = [column for column in frame.columns if column not in ordered]
    return frame[ordered + remaining]


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
    area_sum = grouped["flow_area"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    values = pd.Series(float(fill_value), index=full_index, dtype=float)
    valid = area_sum > 0.0
    values.loc[valid] = (q_sum.loc[valid] / area_sum.loc[valid]).astype(float)
    values = values * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"].agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size().reindex(full_index, fill_value=0).astype(int).tolist(),
        "q_per_area": values.tolist(),
        "q": q_sum.fillna(0.0).astype(float).tolist(),
        "flow_area": area_sum.fillna(0.0).astype(float).tolist(),
    }
    if "lake" in selected.columns:
        lake_strings = selected.assign(lake=selected["lake"].astype("string")).groupby("cell", dropna=False)["lake"]
        hover["lake"] = lake_strings.agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
    if "claktype" in selected.columns:
        hover["claktype"] = (
            grouped["claktype"].agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
        )
    return values.tolist(), hover


def _normalize_term_filter(term: str | Iterable[str] | None) -> list[str] | None:
    """Normalize an optional package-budget term filter to uppercase strings."""

    if term is None:
        return None
    if isinstance(term, str):
        return [str(term).strip().upper()]
    return [str(value).strip().upper() for value in term]


def build_lak_budget_term_table(
    model: "SimulationBase",
    *,
    term: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Return one canonical long dataframe of LAK package-output budget terms.

    Parameters
    ----------
    model
        Live or file-backed model object exposing ``model.outputs.lak.bud``.
    term
        Optional single term or iterable of terms to include. Terms follow the
        MF6 package-output identifiers such as ``"GWF"``, ``"STORAGE"``, or
        ``"TO-MVR"``.

    Returns
    -------
    pandas.DataFrame
        One normalized long table containing a ``term`` column and common
        metadata columns such as ``model``, ``package``, ``kstpkper``, ``per``,
        ``lake``, and ``q``. Term-specific columns, including ``layer``,
        ``cell``, ``iconn``, ``FLOW-AREA``, ``flow_area``, ``q_per_area``, and
        ``VOLUME``, are preserved when applicable.
    """

    requested_terms = _normalize_term_filter(term)
    available_terms = [str(value).strip().upper() for value in model.outputs.lak.bud.types]
    selected_terms = available_terms if requested_terms is None else [value for value in available_terms if value in requested_terms]

    frames: list[pd.DataFrame] = []
    for current_term in selected_terms:
        if current_term == "GWF":
            frame = build_lak_budget_result_table(model, budget_text=current_term, value_name="q").copy()
            if frame.empty:
                continue
            frame["term"] = current_term
            frame["node"] = pd.to_numeric(frame.get("lake"), errors="coerce")
            frames.append(frame)
            continue

        raw = model.outputs.lak.bud.get(current_term)
        if not isinstance(raw, pd.DataFrame):
            continue
        frame = raw.copy()
        if frame.empty:
            continue
        frame = _normalize_budget_nodes(frame)
        frame["model"] = model.name
        frame["package"] = "lak"
        frame["term"] = current_term
        if "kstpkper" in frame.columns:
            frame["kstpkper"] = frame["kstpkper"].apply(lambda values: tuple(int(value) for value in values))
            frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        else:
            frame["kstpkper"] = [(0, 0)] * len(frame)
            frame["per"] = 0
        if "node" in frame.columns:
            frame["lake"] = pd.to_numeric(frame["node"], errors="coerce").astype("Int64")
        if current_term == "FLOW-JA-FACE" and "node2" in frame.columns:
            frame["lake_to"] = pd.to_numeric(frame["node2"], errors="coerce").astype("Int64")
        frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "term",
                "kstpkper",
                "per",
                "lake",
                "q",
            ]
        )

    combined = pd.concat(frames, ignore_index=True, sort=False)
    ordered = [
        "model",
        "package",
        "term",
        "kstpkper",
        "per",
        "lake",
        "lake_to",
        "layer",
        "cell",
        "iconn",
        "claktype",
        "q",
        "FLOW-AREA",
        "flow_area",
        "q_per_area",
        "VOLUME",
        "node",
        "node2",
    ]
    remaining = [column for column in combined.columns if column not in ordered]
    return combined[[column for column in ordered if column in combined.columns] + remaining]


def build_sfr_budget_term_table(
    model: "SimulationBase",
    *,
    term: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Return one canonical long dataframe of SFR package-output budget terms."""

    requested_terms = _normalize_term_filter(term)
    available_terms = [str(value).strip().upper() for value in model.outputs.sfr.bud.types]
    selected_terms = available_terms if requested_terms is None else [value for value in available_terms if value in requested_terms]
    reach_table = build_sfr_reach_table(model)

    frames: list[pd.DataFrame] = []
    for current_term in selected_terms:
        raw = model.outputs.sfr.bud.get(current_term)
        if not isinstance(raw, pd.DataFrame):
            continue
        frame = raw.copy()
        if frame.empty:
            continue
        frame = _normalize_budget_nodes(frame)
        frame["model"] = model.name
        frame["package"] = "sfr"
        frame["term"] = current_term
        if "kstpkper" in frame.columns:
            frame["kstpkper"] = frame["kstpkper"].apply(lambda values: tuple(int(value) for value in values))
            frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        else:
            frame["kstpkper"] = [(0, 0)] * len(frame)
            frame["per"] = 0
        if "node" in frame.columns:
            frame["reach"] = pd.to_numeric(frame["node"], errors="coerce").astype("Int64")
        if "reach" in frame.columns:
            frame = frame.merge(reach_table, on="reach", how="left", suffixes=("", "_pkg"))
        if current_term == "GWF" and "node2" in frame.columns:
            node2_numeric = pd.to_numeric(frame["node2"], errors="coerce")
            layers: list[int | None] = []
            cells: list[int | None] = []
            for node_value in node2_numeric.tolist():
                if pd.isna(node_value):
                    layers.append(None)
                    cells.append(None)
                    continue
                layer_id, cell_id = model.node_to_lni[int(node_value)]
                layers.append(int(layer_id))
                cells.append(int(cell_id))
            frame["layer"] = layers
            frame["cell"] = cells
            q_series = pd.to_numeric(frame.get("q"), errors="coerce")
            rlen_series = pd.to_numeric(frame.get("rlen"), errors="coerce")
            frame["q_per_length"] = np.where(rlen_series > 0.0, q_series / rlen_series, np.nan)
        elif current_term == "FLOW-JA-FACE" and "node2" in frame.columns:
            frame["reach_to"] = pd.to_numeric(frame["node2"], errors="coerce").astype("Int64")
        else:
            if "q" in frame.columns and "rlen" in frame.columns:
                q_series = pd.to_numeric(frame.get("q"), errors="coerce")
                rlen_series = pd.to_numeric(frame.get("rlen"), errors="coerce")
                frame["q_per_length"] = np.where(rlen_series > 0.0, q_series / rlen_series, np.nan)
        frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "term",
                "kstpkper",
                "per",
                "reach",
                "q",
            ]
        )

    combined = pd.concat(frames, ignore_index=True, sort=False)
    ordered = [
        "model",
        "package",
        "term",
        "kstpkper",
        "per",
        "reach",
        "reach_to",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        "q",
        "FLOW-AREA",
        "q_per_length",
        "VOLUME",
        "node",
        "node2",
    ]
    remaining = [column for column in combined.columns if column not in ordered]
    return combined[[column for column in ordered if column in combined.columns] + remaining]


def build_surface_water_exchange_cell_table(
    model: "SimulationBase",
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    include: str | Iterable[str] | None = None,
    lak_connection_type: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Return a normalized cell table combining SFR and LAK exchange intensity.

    Notes
    -----
    The returned ``exchange_intensity`` field uses one unified physical sign
    convention across packages:

    - positive = groundwater gaining into the surface-water feature
    - negative = surface-water losing to groundwater

    LAK already follows that convention in its raw ``q`` sign, so
    ``exchange_intensity = q_per_area``. SFR reports ``GWF`` exchange from the
    stream to groundwater, so its sign is inverted here and
    ``exchange_intensity = -q_per_length``.
    """

    include_packages = _normalize_surface_water_include(include)
    frames: list[pd.DataFrame] = []

    if "sfr" in include_packages:
        sfr_frame = build_sfr_budget_result_table(model, budget_text="SFR", value_name="q")
        sfr_frame = _filter_normalized_table(sfr_frame, per=per, layer=layer, cells=None)
        if not sfr_frame.empty:
            grouped = (
                sfr_frame.groupby(["model", "package", "per", "layer", "cell"], as_index=False)
                .agg(q=("q", "sum"), rlen=("rlen", "sum"))
            )
            grouped["q_per_length"] = np.where(
                pd.to_numeric(grouped["rlen"], errors="coerce") > 0.0,
                pd.to_numeric(grouped["q"], errors="coerce") / pd.to_numeric(grouped["rlen"], errors="coerce"),
                np.nan,
            )
            grouped["exchange_intensity"] = -pd.to_numeric(grouped["q_per_length"], errors="coerce")
            grouped["source"] = "sfr"
            frames.append(grouped)

    if "lak" in include_packages:
        lak_frame = build_lak_budget_result_table(model, budget_text="GWF", value_name="q")
        lak_frame = _filter_normalized_table(lak_frame, per=per, layer=layer, cells=None)
        connection_types = _normalize_connection_type_filter(lak_connection_type)
        if connection_types is not None and "claktype" in lak_frame.columns:
            lak_frame = lak_frame.loc[
                lak_frame["claktype"].astype("string").str.upper().isin(connection_types)
            ].copy()
        if not lak_frame.empty:
            grouped = (
                lak_frame.groupby(["model", "package", "per", "layer", "cell"], as_index=False)
                .agg(q=("q", "sum"), flow_area=("flow_area", "sum"))
            )
            grouped["q_per_area"] = np.where(
                pd.to_numeric(grouped["flow_area"], errors="coerce") > 0.0,
                pd.to_numeric(grouped["q"], errors="coerce") / pd.to_numeric(grouped["flow_area"], errors="coerce"),
                np.nan,
            )
            grouped["exchange_intensity"] = pd.to_numeric(grouped["q_per_area"], errors="coerce")
            grouped["source"] = "lak"
            frames.append(grouped)

    if not frames:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "per",
                "layer",
                "cell",
                "source",
                "exchange_intensity",
            ]
        )

    combined = pd.concat(frames, ignore_index=True, sort=False)
    for column in ("per", "layer", "cell"):
        if column in combined.columns:
            combined[column] = pd.to_numeric(combined[column], errors="coerce").astype(int)
    return combined.reset_index(drop=True)


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
        values = pd.Series(float(fill_value), index=full_index, name="exchange_intensity")
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
    selected["exchange_intensity"] = pd.to_numeric(selected["exchange_intensity"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    values = grouped["exchange_intensity"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    values = values.fillna(float(fill_value)).astype(float) * float(multiplier)

    sfr_component = (
        selected.loc[selected["source"] == "sfr"]
        .groupby("cell", dropna=False)["exchange_intensity"]
        .sum(min_count=1)
        .reindex(full_index, fill_value=0.0)
        .astype(float)
        * float(multiplier)
    )
    lak_component = (
        selected.loc[selected["source"] == "lak"]
        .groupby("cell", dropna=False)["exchange_intensity"]
        .sum(min_count=1)
        .reindex(full_index, fill_value=0.0)
        .astype(float)
        * float(multiplier)
    )

    hover: dict[str, list] = {
        "Package": grouped["package"].agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size().reindex(full_index, fill_value=0).astype(int).tolist(),
        "surface_water_exchange": values.tolist(),
        "sfr_exchange": sfr_component.tolist(),
        "lak_exchange": lak_component.tolist(),
    }
    if "source" in selected.columns:
        hover["source"] = grouped["source"].agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
    return values.tolist(), hover


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
    connectiondata["lake"] = pd.to_numeric(connectiondata["lake"], errors="coerce").astype("Int64")
    if "iconn" in connectiondata.columns:
        connectiondata["iconn"] = pd.to_numeric(connectiondata["iconn"], errors="coerce").astype("Int64")
    connectiondata["layer"] = pd.to_numeric(connectiondata["layer"], errors="coerce").astype(int)
    connectiondata["cell"] = pd.to_numeric(connectiondata["cell"], errors="coerce").astype(int)
    connectiondata["package"] = "lak"
    connectiondata["model"] = model.name
    if "claktype" in connectiondata.columns:
        connectiondata["claktype"] = connectiondata["claktype"].astype("string").str.upper()
    else:
        connectiondata["claktype"] = "UNKNOWN"

    cell_area_lookup = {
        int(cell): float(area)
        for cell, area in enumerate(np.asarray(model.vor.area_list, dtype=float).reshape(-1))
    }
    top = pd.to_numeric(connectiondata.get("telev"), errors="coerce")
    bottom = pd.to_numeric(connectiondata.get("belev"), errors="coerce")
    width = pd.to_numeric(connectiondata.get("connwidth"), errors="coerce")
    vertical_mask = connectiondata["claktype"].astype(str).str.upper().eq("VERTICAL")
    thickness = (top - bottom).clip(lower=0.0)
    horizontal_area = width * thickness
    vertical_area = connectiondata["cell"].map(cell_area_lookup).astype(float)
    connectiondata["connection_area"] = np.where(vertical_mask, vertical_area, horizontal_area)

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


def _numeric_period_settings(package, *, id_column: str, nper: int) -> dict[int, dict[int, dict[str, float]]]:
    """Return carried-forward numeric package settings keyed by period and feature id."""

    current: dict[int, dict[str, float]] = {}
    by_period: dict[int, dict[int, dict[str, float]]] = {}
    perioddata = package.perioddata.get_data()
    for period in range(int(nper)):
        records = perioddata.get(period)
        if records is not None:
            frame = pd.DataFrame(records).copy()
            if not frame.empty:
                identifier = id_column if id_column in frame.columns else frame.columns[0]
                setting_column = next(
                    (column for column in frame.columns if str(column).lower().endswith("setting")),
                    None,
                )
                data_column = next(
                    (column for column in frame.columns if str(column).lower().endswith("setting_data")),
                    None,
                )
                if setting_column is not None and data_column is not None:
                    for row in frame.itertuples(index=False):
                        feature_id = int(getattr(row, identifier))
                        setting = str(getattr(row, setting_column)).strip().lower().replace("-", "_")
                        raw_value = getattr(row, data_column)
                        numeric = pd.to_numeric(pd.Series([raw_value]), errors="coerce").iloc[0]
                        if pd.notna(numeric):
                            current.setdefault(feature_id, {})[setting] = float(numeric)
        by_period[period] = {
            feature_id: dict(settings)
            for feature_id, settings in current.items()
        }
    return by_period


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
    return _coerce_numeric_like_columns(pd.concat(frames, ignore_index=True), exclude={"model", "package"})


def build_lak_input_table(model: "SimulationBase") -> pd.DataFrame:
    """Return normalized static and stress-period LAK inputs mapped to connection cells."""

    connections = build_lak_connection_table(model)
    packagedata = pd.DataFrame(model.lak.packagedata.get_data()).copy()
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "lake"})
    packagedata["lake"] = pd.to_numeric(packagedata["lake"], errors="coerce").astype(int)
    base = connections.merge(packagedata, on="lake", how="left", suffixes=("", "_package"))
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
    return _coerce_numeric_like_columns(pd.concat(frames, ignore_index=True), exclude={"model", "package"})


def build_lak_stage_change_table(model: "SimulationBase") -> pd.DataFrame:
    """Return per-transition lake-stage changes for one model.

    The resulting table carries one row per lake and stress-period transition,
    with ``stage_change = stage_1 - stage_0``.
    """

    stage_table = build_lak_stage_result_table(model)
    if stage_table.empty:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "lake",
                "per0",
                "per1",
                "stage0",
                "stage1",
                "stage_change",
            ]
        )

    base = (
        stage_table.sort_values(["lake", "per", "cell"])
        .drop_duplicates(["lake", "per"])
        .loc[:, ["model", "package", "lake", "per", "stage"]]
        .reset_index(drop=True)
    )
    shifted = base.copy()
    shifted["per0"] = shifted["per"]
    shifted["stage0"] = shifted["stage"]
    shifted["per1"] = shifted.groupby("lake")["per0"].shift(-1)
    shifted["stage1"] = shifted.groupby("lake")["stage0"].shift(-1)
    shifted = shifted.dropna(subset=["per1", "stage1"]).copy()
    shifted["per0"] = shifted["per0"].astype(int)
    shifted["per1"] = shifted["per1"].astype(int)
    shifted["stage_change"] = shifted["stage1"].astype(float) - shifted["stage0"].astype(float)
    return shifted[["model", "package", "lake", "per0", "per1", "stage0", "stage1", "stage_change"]]


def summarize_input_table(frame: pd.DataFrame, *, label: str, value_columns: list[str]) -> pd.DataFrame:
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
                "layers": int(frame["layer"].nunique()) if "layer" in frame.columns else 0,
                "cells": int(frame["cell"].nunique()) if "cell" in frame.columns else 0,
                "value_columns": value_columns,
            }
        ]
    )


def get_default_package_value_column(package_name: str) -> str | None:
    """Return the preferred primary numeric input column for one package."""

    spec = get_package_explorer_spec(package_name)
    return spec.default_input if spec is not None else None


def get_default_package_colorscale(package_name: str) -> str | None:
    """Return the preferred choropleth colorscale for one package or field."""

    normalized = str(package_name).lower()
    if "_" in normalized:
        package, field_name = normalized.split("_", 1)
        field_spec = get_package_input_field_spec(package, field_name)
        if field_spec is not None and field_spec.colorscale is not None:
            return field_spec.colorscale

    spec = get_package_explorer_spec(normalized)
    return spec.colorscale if spec is not None else None


def get_default_budget_term(package_name: str) -> tuple[str, str] | None:
    """Return the preferred budget text and public value name for a package."""

    normalized = str(package_name).lower()
    if "_" in normalized:
        package, result_name = normalized.split("_", 1)
        result_spec = get_package_result_spec(package, result_name)
    else:
        result_spec = get_package_result_spec(normalized, "q")
    if result_spec is None:
        return None
    return result_spec.budget_text, result_spec.value_name


def get_package_explorer_spec(package_name: str) -> PackageExplorerSpec | None:
    """Return registry metadata for one package, if it is known."""

    return _PACKAGE_EXPLORER_SPECS.get(str(package_name).lower())


def get_package_input_field_spec(package_name: str, field_name: str) -> FieldSpec | None:
    """Return registry metadata for one input field, if it is known."""

    spec = get_package_explorer_spec(package_name)
    if spec is None:
        return None
    return spec.inputs.get(str(field_name).lower())


def get_package_result_spec(package_name: str, result_name: str) -> ResultSpec | None:
    """Return registry metadata for one result field, if it is known."""

    spec = get_package_explorer_spec(package_name)
    if spec is None:
        return None
    return spec.results.get(str(result_name).lower())


def get_default_group_compare_colorscale() -> str:
    """Return the default diverging colorscale for grouped difference maps."""

    return _GROUP_COMPARE_DEFAULT_COLORSCALE


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


def _normalize_connection_type_filter(connection_type: str | Iterable[str] | None) -> list[str] | None:
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
        raise ValueError(f"Unsupported surface-water packages: {invalid!r}. Allowed values are 'sfr' and 'lak'.")
    ordered: list[str] = []
    for package_name in ("sfr", "lak"):
        if package_name in normalized:
            ordered.append(package_name)
    return ordered


def _infer_default_value_column(frame: pd.DataFrame, *, fallback: str | None = None) -> str:
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
        raise ValueError("Could not infer a numeric value column for this package table.")
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
    raise ValueError(f"{label} must resolve to one value per layer/cell; got shape={arr.shape}.")


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
    values = grouped[value_column].agg(agg).reindex(full_index, fill_value=fill_value).astype(float) * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"].agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size().reindex(full_index, fill_value=0).astype(int).tolist(),
        value_column: values.tolist(),
    }

    excluded = {"model", "package", "per", "layer", "cell", value_column}
    for column in selected.columns:
        if column in excluded:
            continue
        column_group = grouped[column]
        if pd.api.types.is_numeric_dtype(selected[column]):
            hover[column] = column_group.agg("first").reindex(full_index, fill_value=np.nan).tolist()
        else:
            hover[column] = column_group.agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()
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
    diff_values = grouped[diff_column].agg(agg).reindex(full_index, fill_value=fill_value).astype(float) * float(multiplier)
    current_values = grouped[value_column].agg("first").reindex(full_index, fill_value=np.nan)
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
        "Record Count": grouped.size().reindex(full_index, fill_value=0).astype(int).tolist(),
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
            hover[column] = column_group.agg("first").reindex(full_index, fill_value=np.nan).tolist()
        else:
            hover[column] = column_group.agg(_aggregate_hover_strings).reindex(full_index, fill_value="").tolist()

    absmax = float(np.nanmax(np.abs(diff_values.to_numpy(dtype=float)))) if len(diff_values) else 0.0
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
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=figsize or (5.0 * ncols, 4.0 * nrows),
            squeeze=False,
        )
        arrays = [self._mapped_values(per=per, layer=layer, **map_kwargs) for layer in resolved_layers]
        finite_parts = [values[np.isfinite(values)] for values in arrays if np.isfinite(values).any()]
        finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
        if finite.size:
            vmin = float(np.nanmin(finite)) if vmin is None else vmin
            vmax = float(np.nanmax(finite)) if vmax is None else vmax
        image = None
        for ax, layer, values in zip(axes.flat, resolved_layers, arrays, strict=False):
            view = PlotMapView(model=self.model.gwf, modelgrid=self.model.gwf.modelgrid, layer=layer, ax=ax)
            image = view.plot_array(values, cmap=cmap, vmin=vmin, vmax=vmax)
            if show_grid:
                view.plot_grid(color="#3c4652", linewidth=0.2)
            ax.set_title(f"Layer {layer + 1}")
            ax.set_aspect("equal")
        for ax in axes.flat[len(resolved_layers):]:
            ax.set_visible(False)
        if show_colorbar and image is not None:
            fig.colorbar(image, ax=list(axes.flat[:len(resolved_layers)]), shrink=0.75, label=self._mapped_value_name)
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
            np.asarray(trace.z, dtype=float)[np.isfinite(np.asarray(trace.z, dtype=float))]
            for trace in traces
            if np.isfinite(np.asarray(trace.z, dtype=float)).any()
        ]
        finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
        coloraxis = {"colorscale": traces[0].colorscale if traces else "Viridis"}
        if finite.size:
            coloraxis.update(cmin=float(np.nanmin(finite)), cmax=float(np.nanmax(finite)), cauto=False)
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

        from myflopy.modflow.mf6.interactive_plotting import export_matplotlib_slider_html

        periods = self._mapped_periods() if periods is None else [int(period) for period in periods]
        labels = [f"Stress period {period}" for period in periods]

        def render(period, index):
            return self.plot(per=period, layers=layers, ncols=ncols, title=labels[index], **map_kwargs)

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

        periods = self._mapped_periods() if periods is None else [int(period) for period in periods]
        figures = [
            self.plotly_mosaic(per=period, layers=layers, ncols=ncols, title=title, **map_kwargs)
            for period in periods
        ]
        fig = go.Figure(data=figures[0].data, layout=figures[0].layout)
        fig.frames = [
            go.Frame(data=frame.data, name=str(period))
            for period, frame in zip(periods, figures, strict=False)
        ]
        animation = Animation(self.model, periods=periods, redraw=True)
        fig.update_layout(updatemenus=animation.updatemenus, sliders=animation.sliders, uirevision="lock")
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


class CellPackageInputsExplorer:
    """Normalized input explorer for one cell-based MF6 stress-period package."""

    def __init__(self, model: "SimulationBase", package_name: str):
        self.model = model
        self.package_name = str(package_name).lower()

    def __getattr__(self, field_name: str) -> "CellPackageInputFieldExplorer":
        """Return a field-specific explorer for registry-backed input fields."""

        field_spec = get_package_input_field_spec(self.package_name, field_name)
        if field_spec is None:
            raise AttributeError(f"{type(self).__name__!s} has no input field {field_name!r}")
        return CellPackageInputFieldExplorer(self, field_spec)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a normalized input table for this package.

        Parameters
        ----------
        per
            Optional zero-based stress period.
        layer
            Optional zero-based layer or layers to keep.
        cells
            Optional zero-based cell ids to keep.
        """

        frame = build_cell_package_input_table(
            self.model,
            self.package_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available inputs for this package."""

        frame = self.get()
        value_columns = [
            column
            for column in frame.columns
            if column not in {"model", "package", "per", "layer", "cell"}
        ]
        return summarize_input_table(frame, label=f"{self.package_name}.inputs", value_columns=value_columns)

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a package-input choropleth using the preferred shared style.

        Parameters
        ----------
        per, layer
            Zero-based stress period and layer to map.
        value_column
            Numeric field to display. If omitted, a package-specific default is
            used when available.
        multiplier
            Optional multiplier applied to mapped values before plotting.
        fill_value
            Value used for cells with no package record in the selected
            period/layer.
        agg
            Aggregation used when a package contains multiple rows for one cell
            in the selected period/layer.
        colorscale
            Optional Plotly colorscale name. When omitted, a package-specific
            default is used.
        kwargs
            Forwarded to :meth:`SimulationBase.cor`.
        """

        selected = self.get(per=per, layer=layer)
        chosen_value_column = _infer_default_value_column(
            selected if not selected.empty else self.get(),
            fallback=get_default_package_value_column(self.package_name),
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=value_column or chosen_value_column,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(self.package_name),
            **kwargs,
        )

    @property
    def default(self) -> "CellPackageInputFieldExplorer":
        """Return the registry-defined preferred input field."""

        return getattr(self, get_default_package_value_column(self.package_name))

    def plot(self, **kwargs):
        return self.default.plot(**kwargs)

    def plotly_mosaic(self, **kwargs):
        return self.default.plotly_mosaic(**kwargs)

    def slider_html(self, *args, **kwargs):
        return self.default.slider_html(*args, **kwargs)

    def plotly_animation(self, **kwargs):
        return self.default.plotly_animation(**kwargs)


class CellPackageInputFieldExplorer(MappedFieldVisualizationMixin):
    """Field-specific view over a cell package's normalized input table."""

    def __init__(self, inputs: CellPackageInputsExplorer, field_spec: FieldSpec):
        self.inputs = inputs
        self.field_spec = field_spec
        self.field_name = field_spec.name

    @property
    def model(self):
        """Return the underlying model."""

        return self.inputs.model

    @property
    def package_name(self) -> str:
        """Return the underlying package name."""

        return self.inputs.package_name

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return rows for this input field, preserving package metadata."""

        frame = self.inputs.get(per=per, layer=layer, cells=cells)
        required_columns = ["model", "package", "per", "layer", "cell", self.field_name]
        if frame.empty:
            return pd.DataFrame(columns=required_columns)
        missing = [column for column in required_columns if column not in frame.columns]
        if missing:
            raise KeyError(f"Input field {self.field_name!r} is missing required columns: {missing}")
        return frame.loc[:, required_columns].copy()

    def summary(self) -> pd.DataFrame:
        """Return a compact summary for this input field."""

        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float | None = None,
        agg: str | None = None,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a choropleth for this specific input field."""

        return self.inputs.map(
            per=per,
            layer=layer,
            value_column=self.field_name,
            multiplier=multiplier,
            fill_value=self.field_spec.fill_value if fill_value is None else fill_value,
            agg=self.field_spec.agg if agg is None else agg,
            colorscale=colorscale or self.field_spec.colorscale,
            **kwargs,
        )


class UzfFieldInputsExplorer(MappedFieldVisualizationMixin):
    """Normalized explorer for one UZF perioddata field."""

    def __init__(self, model: "SimulationBase", field_name: str):
        self.model = model
        self.field_name = str(field_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized UZF field table for the selected rows."""

        frame = build_uzf_field_input_table(
            self.model,
            self.field_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of this UZF field's available input data."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"uzf.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def wide(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a wide DataFrame with one row per UZF record and one column per period."""

        return build_uzf_field_input_wide_table(
            self.model,
            self.field_name,
            layer=layer,
            cells=cells,
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a choropleth for the selected UZF field."""

        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(f"uzf_{self.field_name}") or "Viridis",
            **kwargs,
        )


class UzfInputsNamespace:
    """Namespace for normalized UZF input explorers."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported UZF perioddata fields."""

        spec = get_package_explorer_spec("uzf")
        if spec is None:
            return pd.DataFrame(columns=["field", "label", "colorscale", "fill_value", "agg"])
        rows = [
            {
                "field": field_spec.name,
                "label": field_spec.label,
                "colorscale": field_spec.colorscale,
                "fill_value": field_spec.fill_value,
                "agg": field_spec.agg,
            }
            for field_spec in spec.inputs.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported UZF input field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(self._field(field_name).summary())
        if not frames:
            return pd.DataFrame(columns=["label", "records", "periods", "layers", "cells", "value_columns"])
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(lambda values: values[0] if values else None)
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def _field(self, field_name: str) -> UzfFieldInputsExplorer:
        """Return one registry-backed UZF perioddata field explorer."""

        field_spec = get_package_input_field_spec("uzf", field_name)
        if field_spec is None:
            raise AttributeError(f"{type(self).__name__!s} has no UZF input field {field_name!r}")
        return UzfFieldInputsExplorer(self.model, field_spec.name)

    def __getattr__(self, field_name: str) -> UzfFieldInputsExplorer:
        """Return a registry-backed UZF perioddata field explorer."""

        return self._field(field_name)

    @property
    def default(self) -> UzfFieldInputsExplorer:
        """Return the preferred UZF input field."""

        return self.finf

    def map(self, **kwargs):
        return self.default.map(**kwargs)

    def plot(self, **kwargs):
        return self.default.plot(**kwargs)

    def plotly_mosaic(self, **kwargs):
        return self.default.plotly_mosaic(**kwargs)

    def slider_html(self, *args, **kwargs):
        return self.default.slider_html(*args, **kwargs)

    def plotly_animation(self, **kwargs):
        return self.default.plotly_animation(**kwargs)

    @property
    def finf(self) -> UzfFieldInputsExplorer:
        """Return the preferred infiltration-rate explorer."""

        return self._field("finf")

    @property
    def pet(self) -> UzfFieldInputsExplorer:
        """Return the potential evapotranspiration explorer."""

        return self._field("pet")

    @property
    def extdp(self) -> UzfFieldInputsExplorer:
        """Return the ET extinction-depth explorer."""

        return self._field("extdp")

    @property
    def extwc(self) -> UzfFieldInputsExplorer:
        """Return the ET extinction-water-content explorer."""

        return self._field("extwc")

    @property
    def ha(self) -> UzfFieldInputsExplorer:
        """Return the surface-depression-storage-depth explorer."""

        return self._field("ha")

    @property
    def hroot(self) -> UzfFieldInputsExplorer:
        """Return the root-zone-thickness explorer."""

        return self._field("hroot")

    @property
    def rootact(self) -> UzfFieldInputsExplorer:
        """Return the root-activity explorer."""

        return self._field("rootact")


class CellBudgetResultsExplorer(MappedFieldVisualizationMixin):
    """Normalized explorer for one cell-based package result term."""

    def __init__(self, model: "SimulationBase", package_name: str, budget_text: str, value_name: str):
        self.model = model
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized result table for this budget term."""

        frame = build_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            package_name=self.package_name,
            value_name=self.value_name,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the available result rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"{self.package_name}.results.{self.value_name}",
            value_columns=[self.value_name],
        )

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.DataFrame:
        """Pivot this result term to one column per stress period."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        if frame.empty:
            return pd.DataFrame(columns=[*index])
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Wide result index columns were not found: {missing_index}")

        wide = frame.pivot_table(
            index=list(index),
            columns="per",
            values=value_column,
            aggfunc=agg,
        )
        wide.columns = [f"per_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def long(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.Series:
        """Return a long result series indexed by ``kstpkper/layer/cell``."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        index_columns = ["kstpkper", "layer", "cell"]
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays(
                [[] for _ in index_columns],
                names=index_columns,
            )
            return pd.Series([], index=empty_index, dtype=float, name=value_column)
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [column for column in index_columns if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Long result index columns were not found: {missing_index}")

        series = (
            frame.groupby(index_columns, dropna=False)[value_column]
            .agg(agg)
            .sort_index()
        )
        series.name = value_column
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`."""

        return self.long(**kwargs)

    def plot_timeseries(
        self,
        *,
        cells: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        agg: str = "sum",
        ax=None,
        return_fig: bool = True,
    ):
        """Plot this cell result by stress period for selected cells."""

        selected_cells = [int(cells)] if isinstance(cells, (int, np.integer)) else cells
        frame = self.get(layer=layer, cells=selected_cells)
        result_spec = get_package_result_spec(self.package_name, self.value_name)
        display_label = (
            result_spec.label
            if result_spec is not None and result_spec.label is not None
            else f"{self.package_name.upper()} {self.value_name}"
        )
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if frame.empty:
            ax.set_title(f"{display_label} by stress period")
            ax.set_xlabel("Stress Period")
            ax.set_ylabel(self.value_name)
            if return_fig:
                return fig
            return None

        grouped_keys = [column for column in ("layer", "cell") if column in frame.columns]
        for key, group in frame.groupby(grouped_keys, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            key_map = dict(zip(grouped_keys, key, strict=False))
            series = group.groupby("per", as_index=False)[self.value_name].agg(agg).sort_values("per")
            layer_label = f"L{int(key_map['layer'])} " if "layer" in key_map and pd.notna(key_map["layer"]) else ""
            cell_label = f"C{int(key_map['cell'])}" if "cell" in key_map and pd.notna(key_map["cell"]) else "All cells"
            ax.plot(
                series["per"].astype(int).to_numpy(),
                series[self.value_name].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"{layer_label}{cell_label}",
            )

        ax.set_title(f"{display_label} by stress period")
        ax.set_xlabel("Stress Period")
        ax.set_ylabel(self.value_name)
        ax.legend()
        fig.tight_layout()
        if return_fig:
            return fig
        return None

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a choropleth for this cell-based result field."""

        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        if self.value_name == "q":
            absmax = _symmetric_color_limit(values)
            kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
            kwargs.setdefault("zmax", absmax if absmax > 0 else None)
            kwargs.setdefault("zmid", 0.0)
        result_spec = get_package_result_spec(self.package_name, self.value_name)
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or (result_spec.colorscale if result_spec is not None else None)
                or ("RdBu" if self.value_name == "q" else None)
                or get_default_package_colorscale(self.package_name)
                or "Viridis"
            ),
            **kwargs,
        )


class SfrBudgetResultsExplorer(CellBudgetResultsExplorer):
    """SFR-specific result explorer with reach-profile helpers."""

    def __init__(self, model: "SimulationBase", *, budget_text: str = "SFR", value_name: str = "q"):
        super().__init__(model, "sfr", budget_text, value_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized SFR budget-result table for selected rows."""

        frame = build_sfr_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def profile(self, *, per: int = 0) -> pd.DataFrame:
        """Return one reach-ordered profile table for the selected period."""

        frame = self.get(per=per)
        if frame.empty:
            return frame
        return frame.sort_values(["reach", "cell"]).reset_index(drop=True)

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build an SFR exchange choropleth normalized by total reach length.

        Notes
        -----
        The mapped value is ``sum(q) / sum(rlen)`` within each cell. This makes
        SFR exchange maps less sensitive to cells that only appear larger
        because they contain longer stream reaches. The ``agg`` argument is
        accepted for API compatibility but is not used because the
        normalization is computed explicitly from total exchange and total
        length per cell.

        MF6 reports the SFR ``GWF`` budget term as flow from the stream reach
        to the groundwater cell. Positive values therefore indicate losing
        reaches, while negative values indicate gaining reaches. The default
        diverging colorscale is defined explicitly so gaining reaches plot blue
        and losing reaches plot red.
        """

        del agg
        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_sfr_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or _blue_white_red_diverging_colorscale(),
            **kwargs,
        )

    def plot_profile(
        self,
        *,
        per: int = 0,
        x: str = "distance",
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot one SFR exchange profile by reach or cumulative stream distance.

        Parameters
        ----------
        per
            Zero-based stress period to plot.
        x
            Either ``"distance"`` for cumulative stream distance or
            ``"reach"`` for raw reach number.
        plot_fig
            If ``True``, call ``show()`` on the created figure.
        return_fig
            If ``True``, return the created figure.
        """

        frame = self.profile(per=per)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Profile x-axis column {x_column!r} was not found.")
        fig.add_scattergl(
            x=frame[x_column],
            y=frame[self.value_name],
            mode="lines+markers",
            name=f"SFR {self.value_name} per {per}",
            customdata=np.column_stack([frame["reach"], frame["cell"]]),
            hovertemplate=(
                "reach=%{customdata[0]}<br>"
                "cell=%{customdata[1]}<br>"
                f"{self.value_name}=%{{y}}<extra></extra>"
            ),
        )
        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title=self.value_name,
            title=f"SFR {self.value_name} profile (per={per})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None


class LakBudgetResultsExplorer(CellBudgetResultsExplorer):
    """LAK-specific result explorer with area-normalized exchange maps."""

    def __init__(self, model: "SimulationBase", *, budget_text: str = "GWF", value_name: str = "q"):
        super().__init__(model, "lak", budget_text, value_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK budget-result table for selected rows."""

        frame = build_lak_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
        )
        frame = _filter_normalized_table(frame, per=per, layer=layer, cells=cells)
        connection_types = _normalize_connection_type_filter(connection_type)
        if connection_types is not None and "claktype" in frame.columns:
            frame = frame.loc[frame["claktype"].astype("string").str.upper().isin(connection_types)].copy()
        return frame.reset_index(drop=True)

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        connection_type: str | Iterable[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a LAK exchange choropleth normalized by flow-surface area.

        Notes
        -----
        The mapped value is ``sum(q) / sum(flow_area)`` within each cell. This
        yields a signed exchange intensity in length-per-time units instead of
        raw volumetric exchange, which would otherwise scale with lake
        connection area.
        """

        del agg
        selected = self.get(per=per, layer=layer, connection_type=connection_type)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_lak_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "RdBu",
            **kwargs,
        )

    def budget_summary(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Summarize lake-groundwater exchange by period, lake, and connection type.

        Parameters
        ----------
        per
            Optional zero-based stress period filter. When omitted, all periods
            are summarized.

        Returns
        -------
        pandas.DataFrame
            Summary table with signed volumetric exchange ``q`` and
            area-normalized exchange ``q_per_area`` grouped by period, lake,
            and connection type.
        """

        frame = self.get(per=per, connection_type=connection_type)
        if frame.empty:
            return pd.DataFrame(
                columns=[
                    "per",
                    "lake",
                    "claktype",
                    "record_count",
                    "q",
                    "flow_area",
                    "q_per_area",
                ]
            )
        summary = (
            frame.groupby(["per", "lake", "claktype"], dropna=False, as_index=False)
            .agg(
                record_count=("cell", "size"),
                q=("q", "sum"),
                flow_area=("flow_area", "sum"),
            )
            .sort_values(["per", "lake", "claktype"])
            .reset_index(drop=True)
        )
        summary["q_per_area"] = np.where(
            pd.to_numeric(summary["flow_area"], errors="coerce") > 0.0,
            pd.to_numeric(summary["q"], errors="coerce") / pd.to_numeric(summary["flow_area"], errors="coerce"),
            np.nan,
        )
        return summary

    def plot_budget(
        self,
        *,
        per: int = 0,
        connection_type: str | Iterable[str] | None = None,
        value: str = "q_per_area",
        ax=None,
        return_fig: bool = True,
    ):
        """Plot a compact LAK budget summary by connection type for one period.

        Parameters
        ----------
        per
            Zero-based stress period to plot.
        value
            Summary field to visualize. Supported values are ``"q"``,
            ``"flow_area"``, and ``"q_per_area"``.
        ax
            Optional Matplotlib axes object to draw onto.
        return_fig
            If ``True``, return the created figure.
        """

        summary = self.budget_summary(per=per, connection_type=connection_type)
        if value not in {"q", "flow_area", "q_per_area"}:
            raise ValueError("value must be one of: 'q', 'flow_area', 'q_per_area'")
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if summary.empty:
            ax.set_title(f"LAK {value} summary (per={per})")
            ax.set_xlabel("Connection Type")
            ax.set_ylabel(value)
            if return_fig:
                return fig
            return None

        summary = summary.copy()
        labels = summary.apply(
            lambda row: f"Lake {int(row['lake'])}\n{row['claktype']}",
            axis=1,
        )
        colors = []
        if value == "flow_area":
            colors = ["#4c78a8" for _ in range(len(summary))]
        else:
            for current in pd.to_numeric(summary[value], errors="coerce").fillna(0.0):
                colors.append("#1f77b4" if current >= 0.0 else "#d62728")
        ax.bar(labels, summary[value].astype(float), color=colors)
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
        ax.set_title(f"LAK {value} summary (per={per})")
        ax.set_xlabel("Lake / Connection Type")
        ax.set_ylabel(value)
        ax.tick_params(axis="x", rotation=0)
        fig.tight_layout()
        if return_fig:
            return fig
        return None


class StageResultsExplorer(MappedFieldVisualizationMixin):
    """Normalized explorer for cell-mapped stage results such as LAK and SFR."""

    def __init__(self, model: "SimulationBase", package_name: str, builder):
        self.model = model
        self.package_name = str(package_name).lower()
        self._builder = builder

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized stage table for the selected rows."""

        frame = self._builder(self.model)
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available stage results."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"{self.package_name}.results.stage",
            value_columns=["stage"],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "first",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a stage choropleth mapped to cells."""

        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column="stage",
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "Blues",
            **kwargs,
        )


class LakStageResultsExplorer(StageResultsExplorer):
    """LAK stage explorer with simple period-oriented plotting helpers."""

    def __init__(self, model: "SimulationBase"):
        super().__init__(model, "lak", build_lak_stage_result_table)

    def plot_timeseries(
        self,
        *,
        lake: int | None = None,
        ax=None,
        return_fig: bool = True,
    ):
        """Plot stage by stress period for one lake or all lakes.

        Parameters
        ----------
        lake
            Optional zero-based lake id. When omitted, all lakes are plotted.
        ax
            Optional Matplotlib axes object to draw onto.
        return_fig
            If ``True``, return the created figure.
        """

        frame = self.get()
        if lake is not None:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        frame = frame.sort_values(["lake", "per", "cell"]).drop_duplicates(["lake", "per"])
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if frame.empty:
            ax.set_title("LAK stage by stress period")
            ax.set_xlabel("Stress Period")
            ax.set_ylabel("Stage")
            if return_fig:
                return fig
            return None

        for lake_id, group in frame.groupby("lake", dropna=False):
            ax.plot(
                group["per"].astype(int).to_numpy(),
                group["stage"].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"Lake {int(lake_id)}",
            )
        ax.set_title("LAK stage by stress period")
        ax.set_xlabel("Stress Period")
        ax.set_ylabel("Stage")
        ax.legend()
        fig.tight_layout()
        if return_fig:
            return fig
        return None


class LakStageChangeExplorer:
    """Explorer for lake-stage changes between stress periods."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def get(
        self,
        *,
        lake: int | None = None,
        per0: int | None = None,
        per1: int | None = None,
    ) -> pd.DataFrame:
        """Return stage-change rows for one lake and/or one period transition."""

        frame = build_lak_stage_change_table(self.model)
        if lake is not None:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        if per0 is not None:
            frame = frame.loc[frame["per0"] == int(per0)].copy()
        if per1 is not None:
            frame = frame.loc[frame["per1"] == int(per1)].copy()
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available lake-stage transitions."""

        frame = self.get()
        if frame.empty:
            return pd.DataFrame(
                [
                    {
                        "label": "lak.results.stage_change",
                        "records": 0,
                        "lakes": 0,
                        "transitions": 0,
                    }
                ]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.results.stage_change",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "transitions": int(frame[["per0", "per1"]].drop_duplicates().shape[0]),
                }
            ]
        )

    def plot_timeseries(
        self,
        *,
        lake: int | None = None,
        ax=None,
        return_fig: bool = True,
    ):
        """Plot stage changes by stress-period transition."""

        frame = self.get(lake=lake)
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if frame.empty:
            ax.set_title("LAK stage change by transition")
            ax.set_xlabel("Stress-Period Transition")
            ax.set_ylabel("Stage Change")
            if return_fig:
                return fig
            return None

        for lake_id, group in frame.groupby("lake", dropna=False):
            labels = [f"{int(start)}->{int(end)}" for start, end in zip(group["per0"], group["per1"], strict=False)]
            ax.plot(
                labels,
                group["stage_change"].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"Lake {int(lake_id)}",
            )
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
        ax.set_title("LAK stage change by transition")
        ax.set_xlabel("Stress-Period Transition")
        ax.set_ylabel("Stage Change")
        ax.legend()
        fig.tight_layout()
        if return_fig:
            return fig
        return None


class LakConnectionsExplorer:
    """Explorer for LAK connection geometry and exchange interface area."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def get(
        self,
        *,
        lake: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK connection table."""

        frame = build_lak_connection_table(self.model)
        frame = _filter_normalized_table(frame, per=None, layer=layer, cells=cells)
        if lake is not None and "lake" in frame.columns:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of LAK connection geometry."""

        frame = self.get()
        if frame.empty:
            return pd.DataFrame(
                [
                    {
                        "label": "lak.connections",
                        "records": 0,
                        "lakes": 0,
                        "layers": 0,
                        "cells": 0,
                        "total_connection_area": 0.0,
                    }
                ]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.connections",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "layers": int(frame["layer"].nunique()),
                    "cells": int(frame["cell"].nunique()),
                    "total_connection_area": float(pd.to_numeric(frame["connection_area"], errors="coerce").sum()),
                }
            ]
        )

    def map(
        self,
        *,
        lake: int | None = None,
        layer: int = 0,
        value_column: str = "connection_area",
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a choropleth of LAK connection geometry by cell.

        Parameters
        ----------
        lake
            Optional zero-based lake id filter.
        layer
            Zero-based model layer to render.
        value_column
            Connection field to map. Common choices are ``"connection_area"``
            and ``"connwidth"``.
        agg
            Aggregation passed through to the generic cell-input map builder.
        multiplier
            Optional scalar multiplier applied to the mapped values.
        fill_value
            Fill value for cells without lake connections.
        colorscale
            Optional choropleth colorscale override.
        """

        selected = self.get(lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        return self.model.cor(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "Blues",
            **kwargs,
        )


class SfrStageResultsExplorer(StageResultsExplorer):
    """SFR stage explorer with reach-profile helpers."""

    def __init__(self, model: "SimulationBase"):
        super().__init__(model, "sfr", build_sfr_stage_result_table)

    def profile(self, *, per: int = 0) -> pd.DataFrame:
        """Return one reach-ordered stage profile for the selected period."""

        frame = self.get(per=per)
        if frame.empty:
            return frame
        return frame.sort_values(["reach", "cell"]).reset_index(drop=True)

    def plot_profile(
        self,
        *,
        per: int = 0,
        x: str = "distance",
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot one SFR stage profile by reach or cumulative stream distance."""

        frame = self.profile(per=per)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Profile x-axis column {x_column!r} was not found.")
        fig.add_scattergl(
            x=frame[x_column],
            y=frame["stage"],
            mode="lines+markers",
            name=f"SFR stage per {per}",
            customdata=np.column_stack([frame["reach"], frame["cell"]]),
            hovertemplate=(
                "reach=%{customdata[0]}<br>"
                "cell=%{customdata[1]}<br>"
                "stage=%{y}<extra></extra>"
            ),
        )
        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title="stage",
            title=f"SFR stage profile (per={per})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None


class CellPackageResultsNamespace:
    """Namespace for cell-based package results represented by one budget term."""

    def __init__(self, model: "SimulationBase", package_name: str):
        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported package result fields."""

        spec = get_package_explorer_spec(self.package_name)
        if spec is None:
            return pd.DataFrame(columns=["field", "budget_text", "value_name", "label", "colorscale", "diverging"])
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported package result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(getattr(self, field_name).summary())
        if not frames:
            return pd.DataFrame(columns=["label", "records", "periods", "layers", "cells", "value_columns"])
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(lambda values: values[0] if values else None)
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def __getattr__(self, result_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed named result explorer."""

        result_spec = get_package_result_spec(self.package_name, result_name)
        if result_spec is None:
            raise AttributeError(f"{type(self).__name__!s} has no result {result_name!r}")
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
        )

    @property
    def q(self) -> CellBudgetResultsExplorer:
        """Return the primary package-exchange result explorer."""

        result_spec = get_package_result_spec(self.package_name, "q")
        if result_spec is None:
            return CellBudgetResultsExplorer(self.model, self.package_name, self.package_name.upper(), "q")
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
        )


class UzfResultsNamespace:
    """Namespace for UZF result explorers."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported UZF result fields."""

        spec = get_package_explorer_spec("uzf")
        if spec is None:
            return pd.DataFrame(columns=["field", "budget_text", "value_name", "label", "colorscale", "diverging"])
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported UZF result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(self._field(field_name).summary())
        if not frames:
            return pd.DataFrame(columns=["label", "records", "periods", "layers", "cells", "value_columns"])
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(lambda values: values[0] if values else None)
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def _field(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return one registry-backed UZF result explorer."""

        result_spec = get_package_result_spec("uzf", field_name)
        if result_spec is None:
            raise AttributeError(f"{type(self).__name__!s} has no UZF result field {field_name!r}")
        return CellBudgetResultsExplorer(self.model, "uzf", result_spec.budget_text, result_spec.value_name)

    def __getattr__(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed UZF result explorer."""

        return self._field(field_name)

    @property
    def gwrch(self) -> CellBudgetResultsExplorer:
        """Return groundwater recharge from the UZF package."""

        return self._field("gwrch")

    @property
    def sat(self) -> CellBudgetResultsExplorer:
        """Return normalized unsaturated-zone saturation results."""

        return self._field("sat")


class LakResultsNamespace:
    """Namespace for LAK result explorers."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def stage(self) -> LakStageResultsExplorer:
        """Return the lake stage explorer mapped to connected cells."""

        return LakStageResultsExplorer(self.model)

    @property
    def stage_change(self) -> LakStageChangeExplorer:
        """Return the lake-stage change explorer."""

        return LakStageChangeExplorer(self.model)

    @property
    def q(self) -> LakBudgetResultsExplorer:
        """Return the lake-groundwater exchange result explorer.

        Exchange maps are normalized to lake connection area, so
        ``map(...)`` renders ``sum(q) / sum(flow_area)`` by cell.
        """

        budget_text, value_name = get_default_budget_term("lak") or ("GWF", "q")
        return LakBudgetResultsExplorer(self.model, budget_text=budget_text, value_name=value_name)


class LakBudgetNamespace:
    """Namespace for all MF6-defined LAK package-output budget terms."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the LAK package-output budget term names available for the model."""

        return [str(value).strip().upper() for value in self.model.outputs.lak.bud.types]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of LAK package-output budget terms.

        Parameters
        ----------
        term
            Optional LAK budget term filter such as ``"GWF"`` or
            ``["GWF", "STORAGE"]``.
        per
            Optional zero-based stress period filter.
        lakes
            Optional iterable of zero-based lake ids to keep.
        """

        frame = build_lak_budget_term_table(self.model, term=term)
        if per is not None:
            if isinstance(per, Iterable) and not isinstance(per, (str, bytes)):
                periods = {int(value) for value in per}
                frame = frame.loc[frame["per"].isin(periods)].copy()
            else:
                frame = frame.loc[frame["per"] == int(per)].copy()
        if lakes is not None and "lake" in frame.columns:
            if isinstance(lakes, Iterable) and not isinstance(lakes, (str, bytes)):
                lake_ids = {int(value) for value in lakes}
            else:
                lake_ids = {int(lakes)}
            frame = frame.loc[pd.to_numeric(frame["lake"], errors="coerce").isin(lake_ids)].copy()
        return frame.reset_index(drop=True)

    def summary(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
    ) -> pd.DataFrame:
        """Summarize LAK budget terms by selected grouping columns."""

        frame = self.get(term=term, per=per, lakes=lakes)
        group_columns = list(by) if by is not None else ["per", "lake", "term"]
        if frame.empty:
            columns = [*group_columns, "record_count", "q"]
            return pd.DataFrame(columns=columns)

        agg_map: dict[str, tuple[str, str]] = {
            "record_count": ("q", "size"),
            "q": ("q", "sum"),
        }
        for column in ("FLOW-AREA", "flow_area", "VOLUME"):
            if column in frame.columns:
                agg_map[column] = (column, "sum")
        summary = (
            frame.groupby(group_columns, dropna=False, as_index=False)
            .agg(**agg_map)
            .sort_values(group_columns)
            .reset_index(drop=True)
        )
        if "flow_area" in summary.columns:
            summary["q_per_area"] = np.where(
                pd.to_numeric(summary["flow_area"], errors="coerce") > 0.0,
                pd.to_numeric(summary["q"], errors="coerce") / pd.to_numeric(summary["flow_area"], errors="coerce"),
                np.nan,
            )
        return summary

    def wide(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
    ) -> pd.DataFrame:
        """Pivot the long LAK budget table to one wide table by term."""

        frame = self.get(term=term, per=per, lakes=lakes)
        if frame.empty:
            return pd.DataFrame()
        if values not in frame.columns:
            raise KeyError(f"LAK budget column {values!r} was not found.")
        wide = (
            frame.pivot_table(
                index=list(index),
                columns="term",
                values=values,
                aggfunc="sum",
            )
            .sort_index()
        )
        if isinstance(wide.columns, pd.Index):
            wide.columns.name = None
        return wide.reset_index()

    @property
    def gwf(self) -> "PackageBudgetTermExplorer":
        """Lake-groundwater exchange term helper."""

        return PackageBudgetTermExplorer(self, term="GWF", label="lak.budget.gwf")

    @property
    def storage(self) -> "PackageBudgetTermExplorer":
        """Lake storage term helper."""

        return PackageBudgetTermExplorer(self, term="STORAGE", label="lak.budget.storage")

    @property
    def runoff(self) -> "PackageBudgetTermExplorer":
        """Lake runoff term helper."""

        return PackageBudgetTermExplorer(self, term="RUNOFF", label="lak.budget.runoff")

    @property
    def rainfall(self) -> "PackageBudgetTermExplorer":
        """Lake rainfall term helper."""

        return PackageBudgetTermExplorer(self, term="RAINFALL", label="lak.budget.rainfall")

    @property
    def evaporation(self) -> "PackageBudgetTermExplorer":
        """Lake evaporation term helper."""

        return PackageBudgetTermExplorer(self, term="EVAPORATION", label="lak.budget.evaporation")

    @property
    def withdrawal(self) -> "PackageBudgetTermExplorer":
        """Lake withdrawal term helper."""

        return PackageBudgetTermExplorer(self, term="WITHDRAWAL", label="lak.budget.withdrawal")

    @property
    def constant(self) -> "PackageBudgetTermExplorer":
        """Lake constant-stage balancing flow term helper."""

        return PackageBudgetTermExplorer(self, term="CONSTANT", label="lak.budget.constant")

    @property
    def ext_inflow(self) -> "PackageBudgetTermExplorer":
        """External inflow term helper."""

        return PackageBudgetTermExplorer(self, term="EXT-INFLOW", label="lak.budget.ext_inflow")

    @property
    def ext_outflow(self) -> "PackageBudgetTermExplorer":
        """External outflow term helper."""

        return PackageBudgetTermExplorer(self, term="EXT-OUTFLOW", label="lak.budget.ext_outflow")

    @property
    def from_mvr(self) -> "PackageBudgetTermExplorer":
        """Mover inflow term helper."""

        return PackageBudgetTermExplorer(self, term="FROM-MVR", label="lak.budget.from_mvr")

    @property
    def to_mvr(self) -> "PackageBudgetTermExplorer":
        """Mover outflow term helper."""

        return PackageBudgetTermExplorer(self, term="TO-MVR", label="lak.budget.to_mvr")

    @property
    def flow_ja_face(self) -> "PackageBudgetTermExplorer":
        """Lake-to-lake outlet/routing connection term helper."""

        return PackageBudgetTermExplorer(self, term="FLOW-JA-FACE", label="lak.budget.flow_ja_face")

    @property
    def auxiliary(self) -> "PackageBudgetTermExplorer":
        """Auxiliary term helper."""

        return PackageBudgetTermExplorer(self, term="AUXILIARY", label="lak.budget.auxiliary")

    @property
    def mvr(self) -> "PackageBudgetTermExplorer":
        """Combined mover-related LAK budget term helper."""

        return PackageBudgetTermExplorer(self, term=["FROM-MVR", "TO-MVR"], label="lak.budget.mvr")

    @property
    def lake_fluxes(self) -> "PackageBudgetTermExplorer":
        """Combined lake-level flux term helper excluding connection-level GWF rows."""

        return PackageBudgetTermExplorer(
            self,
            term=[
                "EXT-INFLOW",
                "RUNOFF",
                "RAINFALL",
                "EVAPORATION",
                "WITHDRAWAL",
                "STORAGE",
                "CONSTANT",
                "EXT-OUTFLOW",
                "FROM-MVR",
                "TO-MVR",
            ],
            label="lak.budget.lake_fluxes",
        )


class SfrBudgetNamespace:
    """Namespace for all MF6-defined SFR package-output budget terms."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the SFR package-output budget term names available for the model."""

        return [str(value).strip().upper() for value in self.model.outputs.sfr.bud.types]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of SFR package-output budget terms."""

        frame = build_sfr_budget_term_table(self.model, term=term)
        if per is not None:
            if isinstance(per, Iterable) and not isinstance(per, (str, bytes)):
                periods = {int(value) for value in per}
                frame = frame.loc[frame["per"].isin(periods)].copy()
            else:
                frame = frame.loc[frame["per"] == int(per)].copy()
        if reaches is not None and "reach" in frame.columns:
            if isinstance(reaches, Iterable) and not isinstance(reaches, (str, bytes)):
                reach_ids = {int(value) for value in reaches}
            else:
                reach_ids = {int(reaches)}
            frame = frame.loc[pd.to_numeric(frame["reach"], errors="coerce").isin(reach_ids)].copy()
        return frame.reset_index(drop=True)

    def summary(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
    ) -> pd.DataFrame:
        """Summarize SFR budget terms by selected grouping columns."""

        frame = self.get(term=term, per=per, reaches=reaches)
        group_columns = list(by) if by is not None else ["per", "reach", "term"]
        if frame.empty:
            columns = [*group_columns, "record_count", "q"]
            return pd.DataFrame(columns=columns)

        agg_map: dict[str, tuple[str, str]] = {
            "record_count": ("q", "size"),
            "q": ("q", "sum"),
        }
        if "rlen" in frame.columns and "reach" in group_columns:
            agg_map["rlen"] = ("rlen", "first")
        for column in ("FLOW-AREA", "VOLUME"):
            if column in frame.columns:
                agg_map[column] = (column, "sum")
        summary = (
            frame.groupby(group_columns, dropna=False, as_index=False)
            .agg(**agg_map)
            .sort_values(group_columns)
            .reset_index(drop=True)
        )
        if "rlen" in summary.columns:
            summary["q_per_length"] = np.where(
                pd.to_numeric(summary["rlen"], errors="coerce") > 0.0,
                pd.to_numeric(summary["q"], errors="coerce") / pd.to_numeric(summary["rlen"], errors="coerce"),
                np.nan,
            )
        return summary

    def wide(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "reach"),
        values: str = "q",
    ) -> pd.DataFrame:
        """Pivot the long SFR budget table to one wide table by term."""

        frame = self.get(term=term, per=per, reaches=reaches)
        if frame.empty:
            return pd.DataFrame()
        if values not in frame.columns:
            raise KeyError(f"SFR budget column {values!r} was not found.")
        wide = (
            frame.pivot_table(
                index=list(index),
                columns="term",
                values=values,
                aggfunc="sum",
            )
            .sort_index()
        )
        if isinstance(wide.columns, pd.Index):
            wide.columns.name = None
        return wide.reset_index()

    @property
    def gwf(self) -> "PackageBudgetTermExplorer":
        """Stream-groundwater exchange term helper."""

        return PackageBudgetTermExplorer(self, term="GWF", label="sfr.budget.gwf")

    @property
    def flow_ja_face(self) -> "PackageBudgetTermExplorer":
        """Reach-to-reach routing connection term helper."""

        return PackageBudgetTermExplorer(self, term="FLOW-JA-FACE", label="sfr.budget.flow_ja_face")

    @property
    def ext_inflow(self) -> "PackageBudgetTermExplorer":
        """External inflow term helper."""

        return PackageBudgetTermExplorer(self, term="EXT-INFLOW", label="sfr.budget.ext_inflow")

    @property
    def runoff(self) -> "PackageBudgetTermExplorer":
        """Runoff term helper."""

        return PackageBudgetTermExplorer(self, term="RUNOFF", label="sfr.budget.runoff")

    @property
    def rain(self) -> "PackageBudgetTermExplorer":
        """Rainfall term helper."""

        return PackageBudgetTermExplorer(self, term="RAIN", label="sfr.budget.rain")

    @property
    def evaporation(self) -> "PackageBudgetTermExplorer":
        """Evaporation term helper."""

        return PackageBudgetTermExplorer(self, term="EVAPORATION", label="sfr.budget.evaporation")

    @property
    def ext_outflow(self) -> "PackageBudgetTermExplorer":
        """External outflow term helper."""

        return PackageBudgetTermExplorer(self, term="EXT-OUTFLOW", label="sfr.budget.ext_outflow")

    @property
    def storage(self) -> "PackageBudgetTermExplorer":
        """Storage term helper."""

        return PackageBudgetTermExplorer(self, term="STORAGE", label="sfr.budget.storage")

    @property
    def from_mvr(self) -> "PackageBudgetTermExplorer":
        """Mover inflow term helper."""

        return PackageBudgetTermExplorer(self, term="FROM-MVR", label="sfr.budget.from_mvr")

    @property
    def to_mvr(self) -> "PackageBudgetTermExplorer":
        """Mover outflow term helper."""

        return PackageBudgetTermExplorer(self, term="TO-MVR", label="sfr.budget.to_mvr")

    @property
    def auxiliary(self) -> "PackageBudgetTermExplorer":
        """Auxiliary term helper."""

        return PackageBudgetTermExplorer(self, term="AUXILIARY", label="sfr.budget.auxiliary")

    @property
    def mvr(self) -> "PackageBudgetTermExplorer":
        """Combined mover-related SFR budget term helper."""

        return PackageBudgetTermExplorer(self, term=["FROM-MVR", "TO-MVR"], label="sfr.budget.mvr")

    @property
    def stream_fluxes(self) -> "PackageBudgetTermExplorer":
        """Combined reach-level flux term helper excluding GWF and routing rows."""

        return PackageBudgetTermExplorer(
            self,
            term=[
                "EXT-INFLOW",
                "RUNOFF",
                "RAIN",
                "EVAPORATION",
                "EXT-OUTFLOW",
                "STORAGE",
                "FROM-MVR",
                "TO-MVR",
            ],
            label="sfr.budget.stream_fluxes",
        )


class PackageBudgetTermExplorer:
    """Filtered helper for one package budget term or a small term family."""

    def __init__(
        self,
        namespace,
        *,
        term: str | Iterable[str],
        label: str,
    ):
        self._namespace = namespace
        self.term = term
        self.label = label

    @property
    def types(self) -> list[str]:
        """Return the normalized MF6 LAK term names covered by this helper."""

        return _normalize_term_filter(self.term) or []

    def get(
        self,
        *,
        per: int | Iterable[int] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Return the filtered package budget-term dataframe."""

        return self._namespace.get(term=self.term, per=per, **filters)

    def summary(
        self,
        *,
        per: int | Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Summarize the filtered package budget terms."""

        return self._namespace.summary(term=self.term, per=per, by=by, **filters)

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
        **filters,
    ) -> pd.DataFrame:
        """Pivot the filtered package budget terms to a wide dataframe."""

        return self._namespace.wide(term=self.term, per=per, index=index, values=values, **filters)


class SfrResultsNamespace:
    """Namespace for SFR result explorers."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def stage(self) -> StageResultsExplorer:
        """Return the stream stage explorer mapped to reach cells."""

        return SfrStageResultsExplorer(self.model)

    @property
    def q(self) -> SfrBudgetResultsExplorer:
        """Return the stream-groundwater exchange result explorer."""

        budget_text, value_name = get_default_budget_term("sfr") or ("SFR", "q")
        return SfrBudgetResultsExplorer(self.model, budget_text=budget_text, value_name=value_name)

    def long_profile(self, *, per: int = 0) -> pd.DataFrame:
        """Return a merged SFR long-profile table for one stress period.

        The returned table aligns stage, stream-groundwater exchange, and
        packagedata-derived geometry fields by reach. It is useful for custom
        analysis as well as for the higher-level ``plot_long_profile`` helper.
        """

        return build_sfr_long_profile_table(self.model, per=per)

    def plot_long_profile(
        self,
        *,
        per: int = 0,
        x: str = "distance",
        include_stage: bool = True,
        include_streambed: bool = True,
        include_exchange: bool = True,
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot a richer SFR long profile with common hydrologic overlays.

        Parameters
        ----------
        per
            Zero-based stress period to plot.
        x
            Either ``"distance"`` for cumulative stream distance or
            ``"reach"`` for raw reach number.
        include_stage
            Whether to show simulated stream stage.
        include_streambed
            Whether to show streambed top and bottom elevations.
        include_exchange
            Whether to show stream-groundwater exchange on a secondary axis.
        plot_fig
            If ``True``, call ``show()`` on the created figure.
        return_fig
            If ``True``, return the created figure.
        """

        frame = self.long_profile(per=per)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Long-profile x-axis column {x_column!r} was not found.")

        customdata = np.column_stack([frame["reach"], frame["cell"]])
        if include_streambed and "streambed_top" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["streambed_top"],
                mode="lines",
                name="Streambed Top",
                line={"color": "#8c564b", "width": 2},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "streambed_top=%{y}<extra></extra>"
                ),
            )
        if include_streambed and "streambed_bottom" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["streambed_bottom"],
                mode="lines",
                name="Streambed Bottom",
                line={"color": "#c49c94", "width": 2, "dash": "dash"},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "streambed_bottom=%{y}<extra></extra>"
                ),
            )
        if include_stage and "stage" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["stage"],
                mode="lines+markers",
                name="Stage",
                line={"color": "#1f77b4", "width": 3},
                marker={"size": 7},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "stage=%{y}<extra></extra>"
                ),
            )
        if include_exchange and "q" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["q"],
                mode="lines+markers",
                name="Exchange q",
                line={"color": "#d62728", "width": 2},
                marker={"size": 6, "symbol": "diamond"},
                yaxis="y2",
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "q=%{y}<extra></extra>"
                ),
            )
            fig.update_layout(
                yaxis2={
                    "title": "Exchange q",
                    "overlaying": "y",
                    "side": "right",
                    "showgrid": False,
                    "zeroline": True,
                }
            )

        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title="Elevation / Stage",
            title=f"SFR long profile (per={per})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None


class SurfaceWaterExchangeResultsExplorer:
    """Combined SFR/LAK exchange explorer with one shared physical sign scale."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        include: str | Iterable[str] | None = None,
        lak_connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return combined SFR/LAK exchange rows in a unified L/T convention."""

        return build_surface_water_exchange_cell_table(
            self.model,
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of combined surface-water exchange rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label="surface_water.results.q",
            value_columns=["exchange_intensity"],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        include: str | Iterable[str] | None = None,
        lak_connection_type: str | Iterable[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build one combined SFR/LAK exchange map with a shared L/T scale.

        Notes
        -----
        The mapped value uses a unified physical sign convention across SFR and
        LAK:

        - positive = groundwater gaining into the surface-water feature
        - negative = surface-water losing to groundwater
        """

        selected = self.get(
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_surface_water_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "RdBu",
            **kwargs,
        )


class SurfaceWaterResultsNamespace:
    """Namespace for combined surface-water result explorers."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def q(self) -> SurfaceWaterExchangeResultsExplorer:
        """Return one shared SFR/LAK exchange explorer."""

        return SurfaceWaterExchangeResultsExplorer(self.model)


class StaticArrayFieldExplorer(MappedFieldVisualizationMixin):
    """Explorer for static layer/cell arrays such as IC, NPF, and STO fields."""

    def __init__(
        self,
        model: "SimulationBase",
        package_name: str,
        field_name: str,
        *,
        label: str | None = None,
        colorscale: str = "Viridis",
    ):
        self.model = model
        self.package_name = str(package_name).lower()
        self.field_name = str(field_name)
        self.label = label or f"{self.package_name}.{self.field_name}"
        self.colorscale = colorscale

    def _array(self) -> np.ndarray:
        package = self.model.package(self.package_name)
        data = getattr(package, self.field_name)
        values = getattr(data, "array", None)
        if values is None:
            values = getattr(data, "data", data)
        nlay = int(getattr(self.model.gwf.modelgrid, "nlay", getattr(self.model, "nlay", 1)))
        ncpl = int(self.model.vor.ncpl)
        return _as_layer_cell_property(values, nlay=nlay, ncpl=ncpl, label=self.label)

    def get(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return this static array field as a normalized layer/cell table."""

        arr = self._array()
        rows = []
        for layer_index in range(arr.shape[0]):
            for cell in range(arr.shape[1]):
                rows.append(
                    {
                        "model": self.model.name,
                        "package": self.package_name,
                        "field": self.field_name,
                        "layer": layer_index,
                        "cell": cell,
                        self.field_name: arr[layer_index, cell],
                    }
                )
        frame = pd.DataFrame(rows)
        return _filter_normalized_table(frame, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary for this array field."""

        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.arrays.{self.field_name}",
            value_columns=[self.field_name],
        )

    def long(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.Series:
        """Return this field as a series indexed by ``layer/cell``."""

        frame = self.get(layer=layer, cells=cells)
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays([[], []], names=["layer", "cell"])
            return pd.Series([], index=empty_index, dtype=float, name=self.field_name)
        series = frame.set_index(["layer", "cell"])[self.field_name].sort_index()
        series.name = self.field_name
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`."""

        return self.long(**kwargs)

    def wide(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Pivot this field to one row per cell and one column per layer."""

        frame = self.get(layer=layer, cells=cells)
        if frame.empty:
            return pd.DataFrame(columns=["cell"])
        wide = frame.pivot_table(index="cell", columns="layer", values=self.field_name, aggfunc="first")
        wide.columns = [f"layer_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a choropleth for one layer of this static array field."""

        del per
        selected = self.get(layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=None,
            layer=layer,
            agg="first",
        )
        return self.model.cor(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or self.colorscale,
            **kwargs,
        )


class StaticArrayPackageExplorer:
    """Namespace for static array fields in one package."""

    def __init__(self, model: "SimulationBase", package_name: str, fields: Mapping[str, dict[str, str]]):
        self.model = model
        self.package_name = str(package_name).lower()
        self._fields = dict(fields)

    def _available_field_items(self) -> list[tuple[str, dict[str, str]]]:
        package = self.model.package(self.package_name)
        available = []
        for field_name, metadata in self._fields.items():
            data = getattr(package, field_name, None)
            if data is None:
                continue
            available.append((field_name, metadata))
        return available

    @property
    def fields(self) -> pd.DataFrame:
        """Return supported static array fields available on this package."""

        return pd.DataFrame(
            [
                {
                    "field": field_name,
                    "label": metadata.get("label"),
                    "colorscale": metadata.get("colorscale", "Viridis"),
                }
                for field_name, metadata in self._available_field_items()
            ]
        ).reindex(columns=["field", "label", "colorscale"])

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported array field."""

        frames = [self._field(field_name).summary() for field_name in self.fields["field"].tolist()]
        if not frames:
            return pd.DataFrame(columns=["label", "records", "periods", "layers", "cells", "value_columns"])
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(lambda values: values[0] if values else None)
        return summary.merge(self.fields.rename(columns={"field": "field_name"}), on="field_name", how="left")

    def _field(self, field_name: str) -> StaticArrayFieldExplorer:
        metadata = self._fields.get(str(field_name))
        if metadata is None:
            raise AttributeError(f"{type(self).__name__!s} has no array field {field_name!r}")
        package = self.model.package(self.package_name)
        if getattr(package, str(field_name), None) is None:
            raise AttributeError(f"Package {self.package_name!r} has no available array field {field_name!r}")
        return StaticArrayFieldExplorer(
            self.model,
            self.package_name,
            str(field_name),
            label=metadata.get("label"),
            colorscale=metadata.get("colorscale", "Viridis"),
        )

    def __getattr__(self, field_name: str) -> StaticArrayFieldExplorer:
        """Return a supported static array field explorer."""

        return self._field(field_name)


class PackageExplorer:
    """Namespace for one package's preferred exploration helpers."""

    def __init__(self, model: "SimulationBase", package_name: str):
        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def inputs(self) -> CellPackageInputsExplorer:
        """Return normalized input helpers for this package."""

        return CellPackageInputsExplorer(self.model, self.package_name)

    @property
    def results(self) -> CellPackageResultsNamespace:
        """Return normalized result helpers for this package."""

        return CellPackageResultsNamespace(self.model, self.package_name)


class ModelPackages:
    """Preferred package exploration namespace for one model/run.

    Examples
    --------
    ``model.packages.rch.inputs.get()``
        Normalized recharge input table.
    ``model.packages.rch.inputs.map(per=0)``
        Recharge choropleth using the shared choropleth styling.
    ``model.packages.uzf.inputs.finf.summary()``
        Compact summary of UZF infiltration inputs.
    """

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def __getattr__(self, package_name: str) -> PackageExplorer:
        """Return a registry-backed generic package explorer."""

        spec = get_package_explorer_spec(package_name)
        if spec is None or spec.kind != "cell_stress":
            raise AttributeError(f"{type(self).__name__!s} has no package {package_name!r}")
        return PackageExplorer(self.model, spec.name)

    @property
    def rch(self) -> PackageExplorer:
        """Recharge package exploration helpers."""

        return PackageExplorer(self.model, "rch")

    @property
    def chd(self) -> PackageExplorer:
        """Constant-head package exploration helpers."""

        return PackageExplorer(self.model, "chd")

    @property
    def drn(self) -> PackageExplorer:
        """Drain package exploration helpers."""

        return PackageExplorer(self.model, "drn")

    @property
    def ghb(self) -> PackageExplorer:
        """General-head boundary package exploration helpers."""

        return PackageExplorer(self.model, "ghb")

    @property
    def wel(self) -> PackageExplorer:
        """Well package exploration helpers."""

        return PackageExplorer(self.model, "wel")

    @property
    def uzf(self) -> "UzfPackageExplorer":
        """UZF package exploration helpers."""

        return UzfPackageExplorer(self.model)

    @property
    def ic(self) -> StaticArrayPackageExplorer:
        """Initial conditions array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "ic",
            {
                "strt": {"label": "Starting head", "colorscale": "Viridis"},
            },
        )

    @property
    def npf(self) -> StaticArrayPackageExplorer:
        """Node property flow array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "npf",
            {
                "k": {"label": "Horizontal hydraulic conductivity", "colorscale": "Viridis"},
                "k22": {"label": "Horizontal hydraulic conductivity K22", "colorscale": "Viridis"},
                "k33": {"label": "Vertical hydraulic conductivity", "colorscale": "Viridis"},
            },
        )

    @property
    def sto(self) -> StaticArrayPackageExplorer:
        """Storage package array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "sto",
            {
                "ss": {"label": "Specific storage", "colorscale": "Viridis"},
                "sy": {"label": "Specific yield", "colorscale": "Viridis"},
            },
        )

    @property
    def lak(self) -> "LakPackageExplorer":
        """LAK package exploration helpers."""

        return LakPackageExplorer(self.model)

    @property
    def sfr(self) -> "SfrPackageExplorer":
        """SFR package exploration helpers."""

        return SfrPackageExplorer(self.model)

    @property
    def surface_water(self) -> "SurfaceWaterPackageExplorer":
        """Combined SFR/LAK exploration helpers."""

        return SurfaceWaterPackageExplorer(self.model)


class UzfPackageExplorer:
    """Top-level UZF package explorer namespace."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def inputs(self) -> UzfInputsNamespace:
        """Return the UZF input exploration namespace."""

        return UzfInputsNamespace(self.model)

    @property
    def results(self) -> UzfResultsNamespace:
        """Return the UZF result exploration namespace."""

        return UzfResultsNamespace(self.model)


class SurfaceWaterInputFieldExplorer(MappedFieldVisualizationMixin):
    """One cell-mapped LAK or SFR input field."""

    def __init__(self, inputs: "SurfaceWaterInputsNamespace", field_name: str):
        self.inputs = inputs
        self.model = inputs.model
        self.package_name = inputs.package_name
        self.field_name = str(field_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return rows containing this mapped input field."""

        frame = self.inputs.get(per=per, layer=layer, cells=cells)
        metadata = [
            column
            for column in ("model", "package", "per", "lake", "reach", "layer", "cell")
            if column in frame.columns
        ]
        return frame.loc[:, [*metadata, self.field_name]].copy()

    def summary(self) -> pd.DataFrame:
        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str | None = None,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a Plotly choropleth for this LAK or SFR input field."""

        selected = self.get(per=per, layer=layer)
        if agg is None:
            agg = "sum" if self.field_name in {"connection_area", "rlen", "inflow", "runoff"} else "first"
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "Viridis",
            **kwargs,
        )


class SurfaceWaterInputsNamespace:
    """Consistent input exploration namespace for LAK and SFR."""

    def __init__(self, model: "SimulationBase", package_name: str):
        self.model = model
        self.package_name = str(package_name).lower()

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return normalized package inputs mapped to groundwater cells."""

        builder = build_lak_input_table if self.package_name == "lak" else build_sfr_input_table
        return _filter_normalized_table(builder(self.model), per=per, layer=layer, cells=cells)

    @property
    def fields(self) -> pd.DataFrame:
        """Return numeric input fields that can be mapped."""

        frame = self.get()
        excluded = {
            "model", "package", "per", "lake", "reach", "layer", "cell",
            "ifno", "iconn", "ncon", "ndv", "nlakeconn",
        }
        fields = [
            column
            for column in frame.columns
            if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
        ]
        return pd.DataFrame({"field": fields})

    def summary(self) -> pd.DataFrame:
        frames = [getattr(self, field).summary() for field in self.fields["field"].tolist()]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def __getattr__(self, field_name: str) -> SurfaceWaterInputFieldExplorer:
        if field_name not in set(self.fields["field"].tolist()):
            raise AttributeError(f"{type(self).__name__!s} has no input field {field_name!r}")
        return SurfaceWaterInputFieldExplorer(self, field_name)

    @property
    def default(self) -> SurfaceWaterInputFieldExplorer:
        """Return the preferred package input field."""

        preferred = "connection_area" if self.package_name == "lak" else "rhk"
        return getattr(self, preferred)

    def map(self, **kwargs):
        return self.default.map(**kwargs)

    def plot(self, **kwargs):
        return self.default.plot(**kwargs)

    def plotly_mosaic(self, **kwargs):
        return self.default.plotly_mosaic(**kwargs)

    def slider_html(self, *args, **kwargs):
        return self.default.slider_html(*args, **kwargs)

    def plotly_animation(self, **kwargs):
        return self.default.plotly_animation(**kwargs)


class LakPackageExplorer:
    """Top-level LAK package explorer namespace."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped LAK starting-stage, connection, and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "lak")

    @property
    def connections(self) -> LakConnectionsExplorer:
        """Return LAK connection-geometry exploration helpers."""

        return LakConnectionsExplorer(self.model)

    @property
    def budget(self) -> LakBudgetNamespace:
        """Return LAK package-output budget helpers for all MF6-defined terms."""

        return LakBudgetNamespace(self.model)

    @property
    def results(self) -> LakResultsNamespace:
        """Return the LAK result exploration namespace."""

        return LakResultsNamespace(self.model)


class SfrPackageExplorer:
    """Top-level SFR package explorer namespace."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped SFR reach-hydraulic and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "sfr")

    @property
    def budget(self) -> SfrBudgetNamespace:
        """Return SFR package-output budget helpers for all MF6-defined terms."""

        return SfrBudgetNamespace(self.model)

    @property
    def results(self) -> SfrResultsNamespace:
        """Return the SFR result exploration namespace."""

        return SfrResultsNamespace(self.model)


class SurfaceWaterPackageExplorer:
    """Top-level combined surface-water explorer namespace."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def results(self) -> SurfaceWaterResultsNamespace:
        """Return combined SFR/LAK result explorers."""

        return SurfaceWaterResultsNamespace(self.model)
