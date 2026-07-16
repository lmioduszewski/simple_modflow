"""Budget and output-derived table builders for package explorers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_explorer_utils import (
    _aggregate_hover_strings,
    _coerce_numeric_like_columns,
    _filter_normalized_table,
    _normalize_connection_type_filter,
    _normalize_surface_water_include,
    _normalize_term_filter,
    split_cellid_columns,
)
from myflopy.modflow.mf6.package_tables import (
    build_sfr_reach_table,
)


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
    model: SimulationBase,
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


def build_sfr_stage_result_table(model: SimulationBase) -> pd.DataFrame:
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


def build_sfr_long_profile_table(
    model: SimulationBase, *, per: int = 0
) -> pd.DataFrame:
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
    packagedata["streambed_top"] = pd.to_numeric(
        packagedata.get("rtp"), errors="coerce"
    )
    packagedata["streambed_bottom"] = pd.to_numeric(
        packagedata.get("rtp"), errors="coerce"
    ) - pd.to_numeric(packagedata.get("rbth"), errors="coerce")

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

    profile = packagedata.merge(stage, on="reach", how="left").merge(
        q_frame, on="reach", how="left"
    )
    profile["q_per_length"] = np.where(
        pd.to_numeric(profile["rlen"], errors="coerce") > 0.0,
        pd.to_numeric(profile["q"], errors="coerce")
        / pd.to_numeric(profile["rlen"], errors="coerce"),
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
    return (
        profile[ordered + remaining]
        .sort_values(["reach", "cell"])
        .reset_index(drop=True)
    )


def build_sfr_budget_result_table(
    model: SimulationBase,
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
            columns=[
                "model",
                "package",
                "kstpkper",
                "per",
                "reach",
                "layer",
                "cell",
                value_name,
            ]
        )

    reach_map = build_sfr_reach_table(model)
    frame = frame.merge(
        reach_map, on=["layer", "cell"], how="left", suffixes=("", "_pkg")
    )
    frame["reach"] = frame["reach"].astype("Int64")
    if value_name == "q":
        q_series = pd.to_numeric(frame["q"], errors="coerce")
        rlen_series = pd.to_numeric(frame["rlen"], errors="coerce")
        frame["q_per_length"] = np.where(
            rlen_series > 0.0, q_series / rlen_series, np.nan
        )
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


def build_lak_stage_result_table(model: SimulationBase) -> pd.DataFrame:
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
    model: SimulationBase,
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
            columns=[
                "model",
                "package",
                "kstpkper",
                "per",
                "lake",
                "layer",
                "cell",
                value_name,
            ]
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
        frame["kstpkper"] = frame["kstpkper"].apply(
            lambda values: tuple(int(value) for value in values)
        )
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
        connectiondata["iconn"] = pd.to_numeric(
            connectiondata["iconn"], errors="coerce"
        ).astype("Int64")
    connectiondata["lake"] = pd.to_numeric(
        connectiondata["lake"], errors="coerce"
    ).astype("Int64")
    connectiondata["layer"] = pd.to_numeric(
        connectiondata["layer"], errors="coerce"
    ).astype(int)
    connectiondata["cell"] = pd.to_numeric(
        connectiondata["cell"], errors="coerce"
    ).astype(int)
    if "claktype" in connectiondata.columns:
        connectiondata["claktype"] = (
            connectiondata["claktype"].astype("string").str.upper()
        )
    raw_iconn_column = next(
        (column for column in frame.columns if str(column).lower() == "iconn"), None
    )
    use_raw_iconn = False
    if raw_iconn_column is not None and "iconn" in connectiondata.columns:
        if raw_iconn_column != "iconn":
            frame = frame.rename(columns={raw_iconn_column: "iconn"})
        frame["iconn"] = pd.to_numeric(frame["iconn"], errors="coerce").astype("Int64")
        grouped_iconn = frame.groupby(["lake", "layer", "cell"], dropna=False)["iconn"]
        ambiguous_groups = (grouped_iconn.size() > 1) & (
            grouped_iconn.nunique(dropna=False) <= 1
        )
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
                frame[column] = (
                    pd.to_numeric(frame[column], errors="coerce")
                    .astype("Int64")
                    .combine_first(
                        pd.to_numeric(frame[connection_column], errors="coerce").astype(
                            "Int64"
                        )
                    )
                )
                frame = frame.drop(columns=[connection_column])
        if "lake_conn" in frame.columns:
            frame = frame.drop(columns=["lake_conn"])
    else:
        connectiondata["_connection_row_order"] = connectiondata.groupby(
            ["lake", "layer", "cell"]
        ).cumcount()
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
        for column in (
            "iconn",
            "claktype",
            "bedleak",
            "belev",
            "telev",
            "connlen",
            "connwidth",
        )
        if column in connectiondata.columns
    ]
    if fallback_columns:
        fallback_source = connectiondata.loc[
            :, ["lake", "layer", "cell", *fallback_columns]
        ].copy()

        def _first_non_null(series: pd.Series):
            """The first non-null value of a series (NaN if all are null)."""

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
            fallback_source.groupby(
                ["lake", "layer", "cell"], dropna=False, as_index=False
            )
            .agg(fallback_agg)
            .rename(
                columns={column: f"{column}_fallback" for column in fallback_columns}
            )
        )
        frame = frame.merge(fallback, on=["lake", "layer", "cell"], how="left")
        for column in fallback_columns:
            fallback_column = f"{column}_fallback"
            if column in frame.columns and fallback_column in frame.columns:
                missing_mask = frame[column].isna()
                frame.loc[missing_mask, column] = frame.loc[
                    missing_mask, fallback_column
                ]
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


def build_lak_budget_term_table(
    model: SimulationBase,
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
    available_terms = [
        str(value).strip().upper() for value in model.outputs.lak.bud.types
    ]
    selected_terms = (
        available_terms
        if requested_terms is None
        else [value for value in available_terms if value in requested_terms]
    )

    frames: list[pd.DataFrame] = []
    for current_term in selected_terms:
        if current_term == "GWF":
            frame = build_lak_budget_result_table(
                model, budget_text=current_term, value_name="q"
            ).copy()
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
            frame["kstpkper"] = frame["kstpkper"].apply(
                lambda values: tuple(int(value) for value in values)
            )
            frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        else:
            frame["kstpkper"] = [(0, 0)] * len(frame)
            frame["per"] = 0
        if "node" in frame.columns:
            frame["lake"] = pd.to_numeric(frame["node"], errors="coerce").astype(
                "Int64"
            )
        if current_term == "FLOW-JA-FACE" and "node2" in frame.columns:
            frame["lake_to"] = pd.to_numeric(frame["node2"], errors="coerce").astype(
                "Int64"
            )
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
    return combined[
        [column for column in ordered if column in combined.columns] + remaining
    ]


def build_sfr_budget_term_table(
    model: SimulationBase,
    *,
    term: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Return one canonical long dataframe of SFR package-output budget terms."""

    requested_terms = _normalize_term_filter(term)
    available_terms = [
        str(value).strip().upper() for value in model.outputs.sfr.bud.types
    ]
    selected_terms = (
        available_terms
        if requested_terms is None
        else [value for value in available_terms if value in requested_terms]
    )
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
            frame["kstpkper"] = frame["kstpkper"].apply(
                lambda values: tuple(int(value) for value in values)
            )
            frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        else:
            frame["kstpkper"] = [(0, 0)] * len(frame)
            frame["per"] = 0
        if "node" in frame.columns:
            frame["reach"] = pd.to_numeric(frame["node"], errors="coerce").astype(
                "Int64"
            )
        if "reach" in frame.columns:
            frame = frame.merge(
                reach_table, on="reach", how="left", suffixes=("", "_pkg")
            )
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
            frame["q_per_length"] = np.where(
                rlen_series > 0.0, q_series / rlen_series, np.nan
            )
        elif current_term == "FLOW-JA-FACE" and "node2" in frame.columns:
            frame["reach_to"] = pd.to_numeric(frame["node2"], errors="coerce").astype(
                "Int64"
            )
        else:
            if "q" in frame.columns and "rlen" in frame.columns:
                q_series = pd.to_numeric(frame.get("q"), errors="coerce")
                rlen_series = pd.to_numeric(frame.get("rlen"), errors="coerce")
                frame["q_per_length"] = np.where(
                    rlen_series > 0.0, q_series / rlen_series, np.nan
                )
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
    return combined[
        [column for column in ordered if column in combined.columns] + remaining
    ]


def build_surface_water_exchange_cell_table(
    model: SimulationBase,
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
        sfr_frame = build_sfr_budget_result_table(
            model, budget_text="SFR", value_name="q"
        )
        sfr_frame = _filter_normalized_table(
            sfr_frame, per=per, layer=layer, cells=None
        )
        if not sfr_frame.empty:
            grouped = sfr_frame.groupby(
                ["model", "package", "per", "layer", "cell"], as_index=False
            ).agg(q=("q", "sum"), rlen=("rlen", "sum"))
            grouped["q_per_length"] = np.where(
                pd.to_numeric(grouped["rlen"], errors="coerce") > 0.0,
                pd.to_numeric(grouped["q"], errors="coerce")
                / pd.to_numeric(grouped["rlen"], errors="coerce"),
                np.nan,
            )
            grouped["exchange_intensity"] = -pd.to_numeric(
                grouped["q_per_length"], errors="coerce"
            )
            grouped["source"] = "sfr"
            frames.append(grouped)

    if "lak" in include_packages:
        lak_frame = build_lak_budget_result_table(
            model, budget_text="GWF", value_name="q"
        )
        lak_frame = _filter_normalized_table(
            lak_frame, per=per, layer=layer, cells=None
        )
        connection_types = _normalize_connection_type_filter(lak_connection_type)
        if connection_types is not None and "claktype" in lak_frame.columns:
            lak_frame = lak_frame.loc[
                lak_frame["claktype"]
                .astype("string")
                .str.upper()
                .isin(connection_types)
            ].copy()
        if not lak_frame.empty:
            grouped = lak_frame.groupby(
                ["model", "package", "per", "layer", "cell"], as_index=False
            ).agg(q=("q", "sum"), flow_area=("flow_area", "sum"))
            grouped["q_per_area"] = np.where(
                pd.to_numeric(grouped["flow_area"], errors="coerce") > 0.0,
                pd.to_numeric(grouped["q"], errors="coerce")
                / pd.to_numeric(grouped["flow_area"], errors="coerce"),
                np.nan,
            )
            grouped["exchange_intensity"] = pd.to_numeric(
                grouped["q_per_area"], errors="coerce"
            )
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
            combined[column] = pd.to_numeric(combined[column], errors="coerce").astype(
                int
            )
    return combined.reset_index(drop=True)


def build_lak_stage_change_table(model: SimulationBase) -> pd.DataFrame:
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
    shifted["stage_change"] = shifted["stage1"].astype(float) - shifted[
        "stage0"
    ].astype(float)
    return shifted[
        ["model", "package", "lake", "per0", "per1", "stage0", "stage1", "stage_change"]
    ]


__all__ = [
    "_normalize_budget_nodes",
    "build_budget_result_table",
    "build_sfr_stage_result_table",
    "build_sfr_long_profile_table",
    "build_sfr_budget_result_table",
    "build_lak_stage_result_table",
    "build_lak_budget_result_table",
    "build_lak_budget_term_table",
    "build_sfr_budget_term_table",
    "build_surface_water_exchange_cell_table",
    "build_lak_stage_change_table",
]
