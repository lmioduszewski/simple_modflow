"""Observation builders for the first ``simple_modflow`` PEST slice."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from simple_modflow.modflow.mf6.observations import _normalize_row_labels
from simple_modflow.modflow.mf6.pest.specs import (
    DrnFlowObservationSpec,
    HeadTargetObservationSpec,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
)


def _load_values_table(values, *, time_column: str, value_column: str) -> pd.DataFrame:
    """Normalize a long or wide table of target values."""

    if values is None:
        return pd.DataFrame(columns=["time", "name", value_column])
    if isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(f"Target values must include {time_column!r}.")
    if "name" in frame.columns and value_column in frame.columns:
        long = frame.rename(columns={time_column: "time", value_column: "value"}).loc[
            :, ["time", "name", "value"]
        ]
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name="value",
        ).rename(columns={time_column: "time"})
    long["name"] = long["name"].astype(str)
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    long["row_label"] = _normalize_row_labels(long["time"])
    long["col_label"] = long["name"].astype(str).str.strip().str.lower()
    return long


def _observation_name(prefix: str, column_label: str, row_label: str) -> str:
    """Reproduce pyEMU's list-style observation naming convention."""

    return f"oname:{prefix.lower()}_otype:lst_usecol:{column_label}_{row_label}"


def _build_index_row_labels(values: pd.Series, index_name: str) -> pd.Series:
    """Build pyEMU-style row labels for a list-style output index column."""

    normalized_index_name = str(index_name).strip().lower()
    return values.astype(str).str.strip().str.lower().map(lambda value: f"{normalized_index_name}:{value}")


def _assign_target_values(pst, target_frame: pd.DataFrame, *, prefix: str, weight_column: str | None = None):
    """Assign target values and weights to a built ``pyemu.Pst`` object."""

    if target_frame.empty:
        return
    obs = pst.observation_data
    for row in target_frame.itertuples(index=False):
        obsnme = _observation_name(prefix, row.col_label, row.row_label)
        if obsnme not in obs.index:
            continue
        value = getattr(row, "value", getattr(row, "head_target", np.nan))
        obs.loc[obsnme, "obsval"] = float(value)
        if weight_column is not None and hasattr(row, weight_column):
            weight = getattr(row, weight_column)
            if pd.notna(weight):
                obs.loc[obsnme, "weight"] = float(weight)


def _write_frame(frame: pd.DataFrame, path: Path):
    """Write a DataFrame to CSV with a predictable location."""

    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    return path


def _write_locations_snapshot(gdf, path: Path):
    """Write a normalized observation-location snapshot to a vector file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _prepare_named_series_observations(
    project,
    *,
    targets,
    prefix: str,
    value_column: str,
    simulated_wide: pd.DataFrame,
    metadata_kind: str,
    locations_snapshot: pd.DataFrame | gpd.GeoDataFrame | None = None,
):
    """Create one pyEMU list-style observation source from a named-series target set."""

    simulated_path = project.template_workspace / f"{prefix}_simulated_{metadata_kind}.csv"
    _write_frame(simulated_wide, simulated_path)
    use_cols = [column for column in simulated_wide.columns if column != simulated_wide.columns[0]]
    project.pf.add_observations(
        simulated_path.name,
        insfile=f"{simulated_path.name}.ins",
        index_cols=simulated_wide.columns[0],
        use_cols=use_cols,
        prefix=prefix,
    )

    target_frame = targets.to_long().copy()
    if not target_frame.empty:
        target_frame["col_label"] = target_frame["name"].astype(str).str.strip().str.lower()
        target_frame["row_label"] = _build_index_row_labels(target_frame["time"], simulated_wide.columns[0])
        values_path = project.template_workspace / f"{prefix}_target_values.csv"
        values_snapshot = target_frame.loc[:, ["time", "name", value_column]].rename(columns={value_column: "value"})
        _write_frame(values_snapshot.rename(columns={"value": value_column}), values_path)
    else:
        values_path = None

    metadata = {
        "kind": metadata_kind,
        "prefix": prefix,
        "time_column": "time",
        "value_column": value_column,
        "n_locations": int(len(use_cols)),
        "n_rows": int(len(target_frame)),
    }
    if values_path is not None:
        metadata["values_file"] = values_path.name
    if locations_snapshot is not None:
        if isinstance(locations_snapshot, gpd.GeoDataFrame):
            locations_path = project.template_workspace / f"{prefix}_target_locations.gpkg"
            _write_locations_snapshot(locations_snapshot, locations_path)
        else:
            locations_path = project.template_workspace / f"{prefix}_target_locations.csv"
            _write_frame(pd.DataFrame(locations_snapshot), locations_path)
        metadata["locations_file"] = locations_path.name

    return {
        "prefix": prefix,
        "target_frame": target_frame.rename(columns={value_column: "value"}),
        "named_series_forward_run_config": {
            "kind": metadata_kind,
            "locations_file": locations_path.name if locations_snapshot is not None else None,
            "output_csv": simulated_path.name,
        },
        "metadata": metadata,
    }


def prepare_head_target_observations(project, spec: HeadTargetObservationSpec):
    """Create simulated-head CSVs and register pyEMU observations."""

    if spec.simulated_values is None:
        simulated = spec.targets.simulated_heads(project.model)
    elif isinstance(spec.simulated_values, pd.DataFrame):
        simulated = spec.simulated_values.copy()
    else:
        simulated = pd.read_csv(spec.simulated_values)

    simulated_path = project.template_workspace / f"{spec.prefix}_simulated_heads.csv"
    _write_frame(simulated, simulated_path)
    use_cols = [column for column in simulated.columns if column != simulated.columns[0]]
    project.pf.add_observations(
        simulated_path.name,
        insfile=f"{simulated_path.name}.ins",
        index_cols=simulated.columns[0],
        use_cols=use_cols,
        prefix=spec.prefix,
    )
    mapping = spec.targets.match_to_model(project.model).loc[:, ["name", "layer", "cell"]]
    mapping_path = project.template_workspace / f"{spec.prefix}_head_target_map.csv"
    mapping.to_csv(mapping_path, index=False)
    locations_frame = spec.targets.locations_gdf.copy()
    if isinstance(locations_frame, gpd.GeoDataFrame):
        locations_snapshot = locations_frame.loc[
            :,
            ["name", "layer", "group", "weight", locations_frame.geometry.name],
        ].copy()
    else:
        snapshot = spec.targets.match_to_model(project.model).loc[:, ["name", "layer", "group", "weight", "x", "y"]].copy()
        geometry = gpd.points_from_xy(snapshot["x"], snapshot["y"], crs=getattr(project.model.vor, "crs", None))
        locations_snapshot = gpd.GeoDataFrame(
            snapshot.drop(columns=["x", "y"]),
            geometry=geometry,
            crs=getattr(project.model.vor, "crs", None),
        )
    locations_path = project.template_workspace / f"{spec.prefix}_target_locations.gpkg"
    _write_locations_snapshot(locations_snapshot, locations_path)
    values_snapshot = spec.targets.to_long().loc[:, ["time", "name", "head_target"]].rename(
        columns={"head_target": "head"}
    )
    values_path = project.template_workspace / f"{spec.prefix}_target_values.csv"
    _write_frame(values_snapshot, values_path)
    target_frame = spec.targets.get().copy()
    target_frame["col_label"] = target_frame["name"].astype(str).str.strip().str.lower()
    target_frame["row_label"] = _build_index_row_labels(
        target_frame["time"],
        simulated.columns[0],
    )
    return {
        "prefix": spec.prefix,
        "target_frame": target_frame,
        "forward_run_config": {
            "mapping_csv": mapping_path.name,
            "output_csv": simulated_path.name,
        },
        "metadata": {
            "kind": "head_targets",
            "prefix": spec.prefix,
            "locations_file": locations_path.name,
            "values_file": values_path.name,
            "name_column": "name",
            "layer_column": "layer",
            "group_column": "group",
            "weight_column": "weight",
            "time_column": "time",
            "value_column": "head",
            "n_locations": int(locations_snapshot["name"].nunique()),
            "n_rows": int(len(values_snapshot)),
        },
    }


def prepare_lake_stage_observations(project, spec: LakeStageObservationSpec):
    """Create simulated lake-stage CSVs and register pyEMU observations."""

    targets = spec.targets
    if targets is None:
        raise ValueError(
            "LakeStageObservationSpec now expects targets=LakeStageTargets(...) "
            "as the canonical workflow."
        )

    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.get().loc[:, ["name", "lake"]].drop_duplicates()
    prepared = _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="stage",
        simulated_wide=pivot,
        metadata_kind="lake_stage",
        locations_snapshot=locations_snapshot,
    )
    if spec.weight is not None and not prepared["target_frame"].empty:
        prepared["target_frame"]["weight"] = float(spec.weight)
    return prepared


def prepare_sfr_stage_observations(project, spec: SfrStageObservationSpec):
    """Create simulated SFR-stage CSVs and register pyEMU observations."""

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.get().loc[:, ["name", "reach"]].drop_duplicates()
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="stage_target",
        simulated_wide=pivot,
        metadata_kind="sfr_stage",
        locations_snapshot=locations_snapshot,
    )


def prepare_sfr_flow_observations(project, spec: SfrFlowObservationSpec):
    """Create simulated SFR-flow CSVs and register pyEMU observations."""

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.get().loc[:, ["name", "reach"]].drop_duplicates()
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="flow_target",
        simulated_wide=pivot,
        metadata_kind="sfr_flow",
        locations_snapshot=locations_snapshot,
    )


def prepare_drn_flow_observations(project, spec: DrnFlowObservationSpec):
    """Create simulated DRN-zone seepage CSVs and register pyEMU observations."""

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.zone_definitions(project.model)
    if isinstance(locations_snapshot, gpd.GeoDataFrame):
        snapshot = locations_snapshot.loc[:, ["name", "group", "weight", "cells", locations_snapshot.geometry.name]].copy()
    else:
        snapshot = locations_snapshot.loc[:, [column for column in ["name", "group", "weight", "cells"] if column in locations_snapshot.columns]].copy()
    if "cells" in snapshot.columns:
        snapshot["cells"] = snapshot["cells"].apply(lambda values: ",".join(str(value) for value in values))
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="flow_target",
        simulated_wide=pivot,
        metadata_kind="drn_flow",
        locations_snapshot=snapshot,
    )


def finalize_observations(project, prepared: list[dict]):
    """Assign observation values and weights after ``pst`` has been built."""

    for item in prepared:
        target_frame = item.get("target_frame")
        if target_frame is None or target_frame.empty:
            continue
        if "weight" in target_frame.columns:
            _assign_target_values(project.pst, target_frame, prefix=item["prefix"], weight_column="weight")
            continue
        _assign_target_values(project.pst, target_frame, prefix=item["prefix"])
        weight = item.get("weight")
        if weight is None:
            continue
        names = [
            _observation_name(item["prefix"], row.col_label, row.row_label)
            for row in target_frame.itertuples(index=False)
        ]
        existing = [name for name in names if name in project.pst.observation_data.index]
        project.pst.observation_data.loc[existing, "weight"] = float(weight)
