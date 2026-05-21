"""Observation builders for the first ``simple_modflow`` PEST slice."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from simple_modflow.modflow.mf6.observations import _normalize_row_labels
from simple_modflow.modflow.mf6.pest.specs import HeadTargetObservationSpec, LakeStageObservationSpec


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

    stage = project.model.packages.lak.results.stage.get().copy()
    if stage.empty:
        raise ValueError("No LAK stage results are available for the model.")
    if "per" in stage.columns:
        stage = stage.rename(columns={"per": spec.time_column})
    name_column = "lake_name" if "lake_name" in stage.columns else "name"
    if name_column not in stage.columns:
        stage[name_column] = stage["lake"].astype(str)
    pivot = stage.pivot_table(
        index=spec.time_column,
        columns=name_column,
        values="stage",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    simulated_path = project.template_workspace / f"{spec.prefix}_simulated_lake_stage.csv"
    _write_frame(pivot, simulated_path)
    use_cols = [column for column in pivot.columns if column != pivot.columns[0]]
    project.pf.add_observations(
        simulated_path.name,
        insfile=f"{simulated_path.name}.ins",
        index_cols=pivot.columns[0],
        use_cols=use_cols,
        prefix=spec.prefix,
    )
    target_frame = _load_values_table(
        spec.values,
        time_column=spec.time_column,
        value_column=spec.value_column,
    )
    metadata = {
        "kind": "lake_stage",
        "prefix": spec.prefix,
        "time_column": "time",
        "value_column": "stage",
        "weight": spec.weight,
    }
    if not target_frame.empty:
        target_frame["row_label"] = _build_index_row_labels(
            target_frame["time"],
            pivot.columns[0],
        )
        values_path = project.template_workspace / f"{spec.prefix}_target_values.csv"
        values_snapshot = target_frame.loc[:, ["time", "name", "value"]].rename(columns={"value": "stage"})
        _write_frame(values_snapshot, values_path)
        metadata["values_file"] = values_path.name
        metadata["n_rows"] = int(len(values_snapshot))
    return {"prefix": spec.prefix, "target_frame": target_frame, "weight": spec.weight, "metadata": metadata}


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
