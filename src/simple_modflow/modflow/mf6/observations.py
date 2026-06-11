"""Reusable observation and target data helpers for MF6 workflows.

This module intentionally sits outside the calibration-specific ``pest``
package. Observation targets such as groundwater heads are useful for several
tasks in ``simple_modflow``:

- PEST/pyEMU observation setup
- calibration residual review
- plotting and exploratory statistics
- comparing simulated heads against measured targets without a PEST workflow

The first concrete target object implemented here is :class:`HeadTargets`,
which combines a vector file of observation locations with a CSV/DataFrame of
target heads through time or stress period.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _copy_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a shallow copy with a clean integer index."""

    return frame.copy().reset_index(drop=True)


def _normalize_identifier_series(values: pd.Series) -> pd.Series:
    """Normalize identifiers for robust joins and observation naming.

    Parameters
    ----------
    values
        Series of location or column identifiers.

    Returns
    -------
    pandas.Series
        Lower-case, stripped string identifiers.
    """

    return values.astype(str).str.strip().str.lower()


def _normalize_row_labels(values: pd.Series) -> pd.Series:
    """Normalize time/per labels using pyEMU-compatible row naming rules.

    pyEMU lower-cases row labels, trims whitespace, and appends a suffix such as
    ``_2`` when the same row label appears multiple times. Reusing that rule
    lets :class:`HeadTargets` assign target values to the exact observation
    names constructed by ``pyemu.utils.PstFrom.add_observations(...)``.
    """

    labels = values.astype(str).str.strip().str.lower().tolist()
    visit: dict[str, int] = {}
    normalized: list[str] = []
    for label in labels:
        if label in visit:
            visit[label] += 1
            normalized.append(f"{label}_{visit[label]}")
        else:
            visit[label] = 1
            normalized.append(label)
    return pd.Series(normalized, index=values.index, name="row_label")


def _observation_label(value: str) -> str:
    """Return a compact, MF6-friendly observation label fragment."""

    normalized = "".join(
        character if str(character).isalnum() else "_"
        for character in str(value).strip()
    )
    normalized = normalized.strip("_")
    return normalized or "obs"


def _workspace_csv_frames(model) -> list[tuple[Path, pd.DataFrame]]:
    """Return readable CSV tables from the model workspace.

    This intentionally stays shallow: the canonical observation workflow writes
    MF6 observation CSVs directly into the active model workspace, so those
    files are the first place target objects should look when reconstructing
    simulated series after a run or reload.
    """

    workspace = getattr(model, "workspace", None)
    if workspace is None:
        workspace = getattr(model, "model_output_folder_path", None)
    if workspace is None:
        return []
    workspace = Path(workspace)
    if not workspace.exists():
        return []

    frames: list[tuple[Path, pd.DataFrame]] = []
    for path in sorted(workspace.glob("*.csv")):
        lower_name = path.name.lower()
        if "simulated_" in lower_name or "_target_" in lower_name:
            continue
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if frame.empty:
            continue
        frames.append((path, frame))
    return frames


def _clean_observation_output_values(values: pd.Series) -> pd.Series:
    """Normalize MF6 observation CSV values into ordinary numeric series."""

    numeric = pd.to_numeric(values, errors="coerce")
    # MF6 observation outputs often use huge sentinels such as 3e30 for
    # inactive/dry values. For review and calibration aggregation, treat those
    # as zero-flow / missing contributions rather than literal magnitudes.
    numeric = numeric.mask(numeric.abs() > 1.0e20, 0.0)
    return numeric.astype(float)


def _find_named_observation_output(
    model,
    *,
    expected_columns: dict[str, str],
) -> pd.DataFrame | None:
    """Find one workspace CSV containing all requested observation columns.

    Parameters
    ----------
    expected_columns
        Mapping from public target names to the raw MF6 observation column names
        expected in one CSV.
    """

    if not expected_columns:
        return None

    expected_lower = {str(column).strip().lower() for column in expected_columns.values()}
    for _path, frame in _workspace_csv_frames(model):
        lower_map = {str(column).strip().lower(): column for column in frame.columns}
        if not expected_lower.issubset(lower_map):
            continue
        time_column = frame.columns[0]
        result = pd.DataFrame({"time": pd.to_numeric(frame[time_column], errors="coerce")})
        for name, raw_column in expected_columns.items():
            actual = lower_map[str(raw_column).strip().lower()]
            result[str(name)] = _clean_observation_output_values(frame[actual])
        return result
    return None


def _ensure_unique_observation_names(
    frame: pd.DataFrame,
    *,
    base_column: str = "name",
    layer_column: str = "layer",
    cell_column: str = "cell",
) -> pd.Series:
    """Build deterministic unique observation names for FloPy OBS records."""

    counts = frame[base_column].value_counts(dropna=False)
    names: list[str] = []
    used: dict[str, int] = {}
    for row in frame.itertuples(index=False):
        base_name = _observation_label(getattr(row, base_column))
        if counts.get(getattr(row, base_column), 0) > 1:
            base_name = f"{base_name}_lay{int(getattr(row, layer_column))}"
        candidate = base_name
        if candidate in used:
            candidate = f"{candidate}_c{int(getattr(row, cell_column))}"
        if candidate in used:
            used[candidate] += 1
            candidate = f"{candidate}_{used[candidate]}"
        else:
            used[candidate] = 1
        names.append(candidate)
    return pd.Series(names, index=frame.index, name="obsname")


def _load_locations(
    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple,
    *,
    name_column: str,
    layer_column: str | None,
    group_column: str | None,
    weight_column: str | None,
) -> pd.DataFrame:
    """Load and validate a head-target locations layer."""

    if isinstance(locations, gpd.GeoDataFrame):
        frame = locations.copy()
    elif isinstance(locations, pd.Series):
        index_name = locations.index.name or name_column
        frame = locations.rename("cell").rename_axis(index_name).reset_index()
        if index_name != name_column:
            frame = frame.rename(columns={index_name: name_column})
    elif isinstance(locations, dict):
        lower_keys = {str(key).strip().lower() for key in locations.keys()}
        if {"name", "cell"} & lower_keys or {"name", "node"} & lower_keys:
            frame = pd.DataFrame(locations)
        else:
            frame = pd.DataFrame({name_column: list(locations.keys()), "cell": list(locations.values())})
    elif isinstance(locations, (list, tuple)):
        frame = pd.DataFrame(locations)
    elif isinstance(locations, pd.DataFrame):
        frame = locations.copy()
    else:
        frame = gpd.read_file(locations)

    required = [name_column]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(
            "Head target locations are missing required columns: "
            + ", ".join(missing)
        )

    keep_columns = [name_column]
    cell_column = None
    for candidate in ("cell", "node"):
        if candidate in frame.columns:
            cell_column = candidate
            keep_columns.append(candidate)
            break
    for column in (layer_column, group_column, weight_column):
        if column is not None and column in frame.columns:
            keep_columns.append(column)
    if isinstance(frame, gpd.GeoDataFrame):
        if frame.geometry.is_empty.any():
            raise ValueError("Head target locations contain empty geometries.")
        keep_columns.append(frame.geometry.name)
    elif cell_column is None:
        raise ValueError(
            "Non-geospatial head target locations must include a 'cell' or 'node' column."
        )
    frame = frame.loc[:, list(dict.fromkeys(keep_columns))].copy()

    frame = frame.rename(columns={name_column: "name"})
    if cell_column is not None:
        frame = frame.rename(columns={cell_column: "cell"})
    if layer_column and layer_column in frame.columns:
        frame = frame.rename(columns={layer_column: "layer"})
    if group_column and group_column in frame.columns:
        frame = frame.rename(columns={group_column: "group"})
    if weight_column and weight_column in frame.columns:
        frame = frame.rename(columns={weight_column: "weight"})

    if "layer" not in frame.columns:
        frame["layer"] = 0
    else:
        frame["layer"] = pd.to_numeric(frame["layer"], errors="coerce").fillna(0).astype(int)

    if "group" not in frame.columns:
        frame["group"] = pd.NA
    if "weight" not in frame.columns:
        frame["weight"] = np.nan

    frame["name"] = frame["name"].astype(str)
    frame["name_key"] = _normalize_identifier_series(frame["name"])
    if "cell" in frame.columns:
        frame["cell"] = pd.to_numeric(frame["cell"], errors="coerce")
        if frame["cell"].isna().any():
            raise ValueError("Head target cell/node values must be numeric.")
        frame["cell"] = frame["cell"].astype(int)
    if isinstance(frame, gpd.GeoDataFrame):
        frame["x"] = frame.geometry.centroid.x
        frame["y"] = frame.geometry.centroid.y
    else:
        frame["x"] = pd.to_numeric(frame.get("x", np.nan), errors="coerce")
        frame["y"] = pd.to_numeric(frame.get("y", np.nan), errors="coerce")
    return frame


def _normalize_values_frame(
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple,
    *,
    time_column: str,
    value_column: str,
    times=None,
) -> pd.DataFrame:
    """Normalize long or wide target tables to a consistent long format."""

    if isinstance(values, pd.Series):
        index_name = values.index.name or time_column
        frame = values.rename(value_column).rename_axis(index_name).reset_index()
        if index_name != time_column:
            frame = frame.rename(columns={index_name: time_column})
    elif isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            frame.insert(0, time_column, list(times))
    elif isinstance(values, (list, tuple)):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            if frame.shape[1] == 1:
                frame.columns = [value_column]
            if frame.shape[0] != len(list(times)):
                raise ValueError(
                    f"Provided times length {len(list(times))} does not match values length {frame.shape[0]}."
                )
            frame.insert(0, time_column, list(times))
    elif isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(
            f"Head target values must include the time/per column {time_column!r}."
        )

    lower_columns = {str(column).strip().lower(): column for column in frame.columns}
    name_column = lower_columns.get("name")
    explicit_value_column = lower_columns.get(value_column.strip().lower())

    if name_column is not None and explicit_value_column is not None:
        long = frame.rename(
            columns={
                time_column: "time",
                name_column: "name",
                explicit_value_column: "head_target",
            }
        ).loc[:, ["time", "name", "head_target"]]
    elif name_column is None and explicit_value_column is not None and len(frame.columns) == 2:
        long = frame.rename(columns={time_column: "time", explicit_value_column: "head_target"}).loc[
            :, ["time", "head_target"]
        ]
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name="head_target",
        ).rename(columns={time_column: "time"})

    if "name" in long.columns:
        long["name"] = long["name"].astype(str)
        long["name_key"] = _normalize_identifier_series(long["name"])
    long["row_label"] = _normalize_row_labels(long["time"])
    long["head_target"] = pd.to_numeric(long["head_target"], errors="coerce")
    return _copy_frame(long)


def _simulated_heads_by_period(model: Any) -> pd.DataFrame:
    """Build one simulated-head record per period/layer/cell.

    The current implementation keeps the last available time step within each
    stress period because that is typically the most useful calibration target
    when the target data are defined per stress period.
    """

    frame = model.all_heads.reset_index().copy()
    frame["kstp"] = frame["kstpkper"].apply(lambda item: int(item[0]))
    frame["per"] = frame["kstpkper"].apply(lambda item: int(item[1]))
    frame = frame.sort_values(["per", "kstp"]).drop_duplicates(
        subset=["per", "layer", "cell"],
        keep="last",
    )
    frame = frame.rename(columns={"elev": "sim_head"})
    return frame.loc[:, ["per", "layer", "cell", "sim_head"]].reset_index(drop=True)


def _default_obs_csv_name(filename: str | Path, default_name: str) -> str:
    """Derive a default MF6 observation CSV name from an OBS control filename."""

    if filename is None:
        return default_name
    stem = Path(filename).stem
    return f"{stem}.csv"


def _coerce_label_list(values, *, count: int, default_prefix: str) -> list[str]:
    """Normalize optional names/groups to a length-matched list."""

    if values is None:
        return [f"{default_prefix}_{idx + 1:03d}" for idx in range(count)]
    if isinstance(values, str):
        if count != 1:
            raise ValueError(f"{default_prefix} label scalar only valid for one target.")
        return [values]
    items = list(values)
    if len(items) != count:
        raise ValueError(f"Expected {count} {default_prefix} labels, received {len(items)}.")
    return [str(item) for item in items]


def _coerce_scalar_or_list(values, *, count: int, default=None):
    """Broadcast scalars or validate sequence lengths."""

    if values is None:
        return [default] * count
    if isinstance(values, str) or not isinstance(values, (list, tuple, np.ndarray, pd.Series)):
        return [values] * count
    items = list(values)
    if len(items) != count:
        raise ValueError(f"Expected {count} values, received {len(items)}.")
    return items


def _build_values_from_cells(
    values,
    *,
    names: list[str],
    times,
    time_column: str,
    value_column: str,
) -> pd.DataFrame:
    """Normalize cell-based target values into the standard long-format input."""

    if isinstance(values, pd.Series):
        if len(names) != 1:
            raise ValueError("Series input for from_cells(...) only supports one target name.")
        frame = values.rename(names[0]).reset_index()
        frame.columns = [time_column, names[0]]
        return frame
    if isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns:
            if times is None:
                raise ValueError("times is required when values dict does not include the time column.")
            frame.insert(0, time_column, list(times))
        return frame
    if isinstance(values, pd.DataFrame):
        frame = values.copy()
        if time_column not in frame.columns:
            if times is None:
                raise ValueError("times is required when values DataFrame does not include the time column.")
            frame.insert(0, time_column, list(times))
        unnamed = [column for column in frame.columns if column != time_column]
        if unnamed == list(range(len(unnamed))):
            rename = {column: names[idx] for idx, column in enumerate(unnamed)}
            frame = frame.rename(columns=rename)
        return frame

    array = np.asarray(values, dtype=object)
    if array.ndim == 1:
        if len(names) != 1:
            raise ValueError("One-dimensional values for from_cells(...) require exactly one target.")
        if times is None:
            times = list(range(len(array)))
        return pd.DataFrame({time_column: list(times), names[0]: array.tolist()})
    if array.ndim != 2:
        raise ValueError("from_cells(...) values must be 1D or 2D.")
    if array.shape[1] != len(names):
        raise ValueError(
            f"from_cells(...) expected {len(names)} value columns, received {array.shape[1]}."
        )
    if times is None:
        times = list(range(array.shape[0]))
    return pd.DataFrame(array.tolist(), columns=names).assign(**{time_column: list(times)})[
        [time_column, *names]
    ]


def _normalize_lake_stage_values(
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None,
    *,
    time_column: str,
    value_column: str,
    times=None,
) -> pd.DataFrame:
    """Normalize lake-stage target values to a long format."""

    if values is None:
        return pd.DataFrame(columns=["time", "name", "stage"])
    if isinstance(values, pd.Series):
        index_name = values.index.name or time_column
        frame = values.rename(value_column).rename_axis(index_name).reset_index()
        if index_name != time_column:
            frame = frame.rename(columns={index_name: time_column})
    elif isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            frame.insert(0, time_column, list(times))
    elif isinstance(values, (list, tuple)):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            if frame.shape[1] == 1:
                frame.columns = [value_column]
            if frame.shape[0] != len(list(times)):
                raise ValueError(
                    f"Provided times length {len(list(times))} does not match values length {frame.shape[0]}."
                )
            frame.insert(0, time_column, list(times))
    elif isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(f"Lake-stage target values must include {time_column!r}.")
    lower_columns = {str(column).strip().lower(): column for column in frame.columns}
    name_column = lower_columns.get("name")
    explicit_value_column = lower_columns.get(value_column.strip().lower())
    if name_column is not None and explicit_value_column is not None:
        long = frame.rename(
            columns={
                time_column: "time",
                name_column: "name",
                explicit_value_column: "stage",
            }
        ).loc[:, ["time", "name", "stage"]]
    elif name_column is None and explicit_value_column is not None and len(frame.columns) == 2:
        long = frame.rename(columns={time_column: "time", explicit_value_column: "stage"}).loc[
            :, ["time", "stage"]
        ]
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name="stage",
        ).rename(columns={time_column: "time"})
    if "name" in long.columns:
        long["name"] = long["name"].astype(str)
    long["stage"] = pd.to_numeric(long["stage"], errors="coerce")
    return _copy_frame(long)


def _normalize_named_integer_locations(
    locations,
    *,
    id_column: str,
    target_name: str,
) -> pd.DataFrame:
    """Normalize simple ``name -> integer id`` location definitions."""

    source = {} if locations is None else locations
    if isinstance(source, pd.Series):
        source = source.to_dict()
    elif isinstance(source, pd.DataFrame):
        required = {"name", id_column}
        if required - set(source.columns):
            missing = sorted(required - set(source.columns))
            raise ValueError(f"{target_name} locations DataFrame is missing columns: {', '.join(missing)}")
        source = source.loc[:, ["name", id_column]].drop_duplicates(subset=["name"]).set_index("name")[id_column].to_dict()
    elif isinstance(source, (list, tuple)):
        source_frame = pd.DataFrame(source)
        required = {"name", id_column}
        if required - set(source_frame.columns):
            missing = sorted(required - set(source_frame.columns))
            raise ValueError(f"{target_name} location records are missing columns: {', '.join(missing)}")
        source = (
            source_frame.loc[:, ["name", id_column]]
            .drop_duplicates(subset=["name"])
            .set_index("name")[id_column]
            .to_dict()
        )

    normalized = {
        str(name): int(value)
        for name, value in dict(source).items()
    }
    return pd.DataFrame([{"name": name, id_column: value} for name, value in normalized.items()])


def _normalize_named_series_values(
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None,
    *,
    time_column: str,
    value_column: str,
    output_column: str,
    times=None,
) -> pd.DataFrame:
    """Normalize long or wide named-series targets to one long table."""

    if values is None:
        return pd.DataFrame(columns=["time", "name", output_column])
    if isinstance(values, pd.Series):
        index_name = values.index.name or time_column
        frame = values.rename(value_column).rename_axis(index_name).reset_index()
        if index_name != time_column:
            frame = frame.rename(columns={index_name: time_column})
    elif isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            frame.insert(0, time_column, list(times))
    elif isinstance(values, (list, tuple)):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            if frame.shape[1] == 1:
                frame.columns = [value_column]
            if frame.shape[0] != len(list(times)):
                raise ValueError(
                    f"Provided times length {len(list(times))} does not match values length {frame.shape[0]}."
                )
            frame.insert(0, time_column, list(times))
    elif isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(f"Target values must include {time_column!r}.")

    lower_columns = {str(column).strip().lower(): column for column in frame.columns}
    name_column = lower_columns.get("name")
    explicit_value_column = lower_columns.get(value_column.strip().lower())

    if name_column is not None and explicit_value_column is not None:
        long = frame.rename(
            columns={
                time_column: "time",
                name_column: "name",
                explicit_value_column: output_column,
            }
        ).loc[:, ["time", "name", output_column]]
    elif name_column is None and explicit_value_column is not None and len(frame.columns) == 2:
        long = frame.rename(columns={time_column: "time", explicit_value_column: output_column}).loc[
            :, ["time", output_column]
        ]
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name=output_column,
        ).rename(columns={time_column: "time"})
    if "name" in long.columns:
        long["name"] = long["name"].astype(str)
    long[output_column] = pd.to_numeric(long[output_column], errors="coerce")
    return _copy_frame(long)


def _numeric_periods(frame: pd.DataFrame, *, caller: str) -> pd.Series:
    """Return zero-based stress periods parsed from a target/result frame."""

    if "per" in frame.columns:
        per = pd.to_numeric(frame["per"], errors="coerce")
    elif "time" in frame.columns:
        per = pd.to_numeric(frame["time"], errors="coerce")
    else:
        raise ValueError(f"{caller} requires a 'per' or 'time' column.")
    if per.isna().any():
        raise ValueError(f"{caller} currently expects numeric zero-based stress periods.")
    return per.astype(int)


def _aggregate_named_series(
    frame: pd.DataFrame,
    *,
    id_column: str,
    value_column: str,
) -> pd.DataFrame:
    """Aggregate model result rows to one value per period and integer id."""

    data = frame.copy()
    data["per"] = _numeric_periods(data, caller="_aggregate_named_series(...)")
    data[id_column] = pd.to_numeric(data[id_column], errors="coerce")
    data[value_column] = pd.to_numeric(data[value_column], errors="coerce")
    data = data.dropna(subset=[id_column, value_column]).copy()
    data[id_column] = data[id_column].astype(int)
    return (
        data.groupby(["per", id_column], dropna=False, as_index=False)[value_column]
        .sum()
        .sort_values(["per", id_column])
        .reset_index(drop=True)
    )


def _default_compare_stats(frame: pd.DataFrame, *, residual_column: str = "residual") -> pd.DataFrame:
    """Return basic residual statistics for a compare-style frame."""

    data = frame.dropna(subset=[residual_column]).copy()
    if data.empty:
        return pd.DataFrame([{"n": 0, "mean_error": np.nan, "mae": np.nan, "rmse": np.nan}])
    residual = data[residual_column].to_numpy(dtype=float)
    return pd.DataFrame(
        [
            {
                "n": int(len(data)),
                "mean_error": float(np.mean(residual)),
                "mae": float(np.mean(np.abs(residual))),
                "rmse": float(np.sqrt(np.mean(residual**2))),
            }
        ]
    )


@dataclass
class HeadTargets:
    """Observation locations plus target heads through time or stress period.

    Parameters
    ----------
    locations
        Path to a vector dataset or an in-memory ``GeoDataFrame`` describing
        observation locations.
    values
        CSV path or DataFrame containing target heads. Both long and wide forms
        are supported.
    name_column
        Name of the observation identifier column in ``locations``.
    layer_column
        Optional layer column in ``locations``. Defaults to layer 0 when absent.
    group_column
        Optional observation-group column in ``locations``.
    weight_column
        Optional observation-weight column in ``locations``.
    time_column
        Column in ``values`` representing stress period, time stamp, or another
        row identifier shared with simulated results.
    value_column
        Value column name used when ``values`` is supplied in long form.
    """

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple
    name_column: str = "name"
    layer_column: str | None = "layer"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "head"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Load and normalize locations and target tables."""

        self._locations = _load_locations(
            self.locations,
            name_column=self.name_column,
            layer_column=self.layer_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
        )
        self._values = _normalize_values_frame(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            times=self.times,
        )
        if "name" not in self._values.columns:
            if len(self._locations["name_key"].unique()) != 1:
                raise ValueError(
                    "HeadTargets values without explicit names require exactly one target location."
                )
            only_name = str(self._locations["name"].iloc[0])
            self._values["name"] = only_name
            self._values["name_key"] = _normalize_identifier_series(self._values["name"])
        missing = sorted(
            set(self._values["name_key"].unique()) - set(self._locations["name_key"].unique())
        )
        if missing:
            raise ValueError(
                "Head target values reference locations not found in the locations layer: "
                + ", ".join(missing[:10])
            )

    @classmethod
    def from_cells(
        cls,
        *,
        cells,
        values,
        times=None,
        names=None,
        layers=0,
        groups=None,
        weights=None,
        time_column: str = "time",
        value_column: str = "head",
    ):
        """Create head targets directly from explicit model cell ids."""

        cell_list = [int(value) for value in list(cells)]
        count = len(cell_list)
        if count == 0:
            raise ValueError("from_cells(...) requires at least one cell.")
        names_list = _coerce_label_list(names, count=count, default_prefix="OBS")
        layer_list = [int(value) for value in _coerce_scalar_or_list(layers, count=count, default=0)]
        group_list = _coerce_scalar_or_list(groups, count=count, default=pd.NA)
        weight_list = _coerce_scalar_or_list(weights, count=count, default=np.nan)
        locations = pd.DataFrame(
            {
                "name": names_list,
                "layer": layer_list,
                "cell": cell_list,
                "group": group_list,
                "weight": weight_list,
            }
        )
        values_frame = _build_values_from_cells(
            values,
            names=names_list,
            times=times,
            time_column=time_column,
            value_column=value_column,
        )
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
        )

    @classmethod
    def from_records(
        cls,
        records,
        *,
        time_column: str = "time",
        value_column: str = "head",
    ):
        """Create head targets from long-format record dictionaries/dataframes."""

        frame = pd.DataFrame(records).copy()
        if frame.empty:
            raise ValueError("from_records(...) requires at least one record.")
        if "name" not in frame.columns:
            raise ValueError("from_records(...) requires a 'name' column.")
        if time_column not in frame.columns:
            raise ValueError(f"from_records(...) requires the time column {time_column!r}.")
        if value_column not in frame.columns:
            raise ValueError(f"from_records(...) requires the value column {value_column!r}.")
        if "cell" not in frame.columns and "geometry" not in frame.columns and "node" not in frame.columns:
            raise ValueError("from_records(...) requires 'cell', 'node', or 'geometry' columns.")
        locations_columns = [
            column
            for column in ["name", "layer", "cell", "node", "group", "weight", "geometry", "x", "y"]
            if column in frame.columns
        ]
        locations = frame.loc[:, locations_columns].drop_duplicates(subset=["name"]).copy()
        if "geometry" in locations.columns:
            locations = gpd.GeoDataFrame(locations, geometry="geometry")
        values_frame = frame.loc[:, [time_column, "name", value_column]].copy()
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
        )

    @property
    def locations_gdf(self) -> gpd.GeoDataFrame:
        """Normalized observation locations table."""

        return self._locations.copy()

    def get(self) -> pd.DataFrame:
        """Return the merged long-format target table."""

        if isinstance(self._locations, gpd.GeoDataFrame):
            location_frame = self._locations.drop(columns=self._locations.geometry.name)
        else:
            location_frame = self._locations.copy()
        merged = self._values.merge(
            location_frame,
            on="name_key",
            how="left",
            suffixes=("", "_loc"),
        )
        if "name_loc" in merged.columns:
            merged = merged.drop(columns=["name"]).rename(columns={"name_loc": "name"})
        return _copy_frame(merged)

    def to_long(self) -> pd.DataFrame:
        """Return the long-format target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per location."""

        frame = self.get().pivot_table(
            index="time",
            columns="name",
            values="head_target",
            aggfunc="first",
        ).reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the target dataset."""

        frame = self.get()
        summary = {
            "n_locations": int(frame["name_key"].nunique()),
            "n_rows": int(len(frame)),
            "n_groups": int(frame["group"].dropna().nunique()) if "group" in frame else 0,
            "n_weighted": int(frame["weight"].notna().sum()) if "weight" in frame else 0,
        }
        return pd.DataFrame([summary])

    def plot_locations(self, *, ax=None, **kwargs):
        """Plot target locations and return the matplotlib axis."""

        if not isinstance(self._locations, gpd.GeoDataFrame):
            raise ValueError("plot_locations() requires geospatial target locations.")
        if ax is None:
            _, ax = plt.subplots()
        self._locations.plot(ax=ax, **kwargs)
        ax.set_title("Head Target Locations")
        return ax

    def match_to_model(self, model) -> pd.DataFrame:
        """Map target locations to zero-based model cell ids.

        Parameters
        ----------
        model
            ``SimulationBase`` or another model-like object exposing
            ``vor.gdf_vorPolys``.

        Returns
        -------
        pandas.DataFrame
            One row per target location with mapped ``cell`` and ``layer``.
        """

        if "cell" in self._locations.columns:
            matches = self._locations.copy()
            if matches["x"].isna().any() or matches["y"].isna().any():
                vor = model.vor.gdf_vorPolys.reset_index().rename(columns={"index": "cell"})
                centroids = pd.DataFrame(
                    {
                        "cell": vor["cell"].astype(int),
                        "x": vor.geometry.centroid.x,
                        "y": vor.geometry.centroid.y,
                    }
                )
                matches = matches.drop(columns=["x", "y"], errors="ignore").merge(centroids, on="cell", how="left")
            return _copy_frame(
                matches.loc[:, ["name", "name_key", "layer", "group", "weight", "cell", "x", "y"]]
            )

        vor = model.vor.gdf_vorPolys.reset_index().rename(columns={"index": "cell"})
        joined = gpd.sjoin(
            self._locations,
            vor.loc[:, ["cell", vor.geometry.name]],
            how="left",
            predicate="intersects",
        )
        if joined["cell"].isna().any():
            unmatched = joined["cell"].isna()
            nearest = gpd.sjoin_nearest(
                joined.loc[unmatched, self._locations.columns],
                vor.loc[:, ["cell", vor.geometry.name]],
                how="left",
            )
            joined.loc[unmatched, "cell"] = nearest["cell"].to_numpy()
        joined["cell"] = joined["cell"].astype(int)
        return _copy_frame(
            joined.drop(columns=["index_right"]).loc[
                :, ["name", "name_key", "layer", "group", "weight", "cell", "x", "y"]
            ]
        )

    def simulated_heads(self, model) -> pd.DataFrame:
        """Return simulated heads for the target locations by stress period.

        The result is a wide table suitable for pyEMU list-style observation
        setup: one row per stress period and one column per observation name.
        """

        matches = self.match_to_model(model)
        sim = _simulated_heads_by_period(model)
        merged = matches.merge(sim, on=["layer", "cell"], how="left")
        wide = merged.pivot_table(
            index="per",
            columns="name",
            values="sim_head",
            aggfunc="first",
        ).reset_index()
        wide.columns.name = None
        return wide

    def compare(self, model) -> pd.DataFrame:
        """Compare target heads against simulated heads by stress period.

        The first implementation assumes the target table's ``time_column`` can
        be interpreted directly as a zero-based stress period integer. This
        matches the initial Cumberland calibration use case where targets are
        defined per stress period.
        """

        targets = self.get().copy()
        per = pd.to_numeric(targets["time"], errors="coerce")
        if per.isna().any():
            raise ValueError(
                "HeadTargets.compare(...) currently expects the time column to "
                "contain zero-based stress period integers."
            )
        targets["per"] = per.astype(int)
        sim = _simulated_heads_by_period(model)
        if "cell" in targets.columns and targets["cell"].notna().all():
            merged = targets.merge(sim, on=["per", "layer", "cell"], how="left")
            if (targets["x"].isna().any() or targets["y"].isna().any()) and hasattr(model, "vor"):
                matches = self.match_to_model(model).loc[:, ["name", "name_key", "layer", "cell", "x", "y"]]
                merged = merged.drop(columns=["x", "y"], errors="ignore").merge(
                    matches,
                    on=["name", "name_key", "layer", "cell"],
                    how="left",
                )
        else:
            matches = self.match_to_model(model).loc[:, ["name", "name_key", "layer", "cell", "x", "y"]]
            merged = (
                targets.merge(matches, on=["name", "name_key", "layer"], how="left")
                .merge(sim, on=["per", "layer", "cell"], how="left")
            )
        merged["residual"] = merged["sim_head"] - merged["head_target"]
        merged["abs_residual"] = merged["residual"].abs()
        return _copy_frame(merged)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table."""

        frame = self.compare(model)
        return frame.loc[
            :,
            [
                "name",
                "group",
                "time",
                "per",
                "layer",
                "cell",
                "head_target",
                "sim_head",
                "residual",
                "abs_residual",
                "weight",
            ],
        ]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the target set."""

        frame = self.compare(model).dropna(subset=["residual"])
        if frame.empty:
            return pd.DataFrame(
                [{"n": 0, "mean_error": np.nan, "mae": np.nan, "rmse": np.nan}]
            )
        residual = frame["residual"].to_numpy(dtype=float)
        stats = {
            "n": int(len(frame)),
            "mean_error": float(np.mean(residual)),
            "mae": float(np.mean(np.abs(residual))),
            "rmse": float(np.sqrt(np.mean(residual**2))),
        }
        return pd.DataFrame([stats])

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        """Return a FloPy ``continuous`` dict for MF6 head observations."""

        matches = self.match_to_model(model).copy()
        matches["obsname"] = _ensure_unique_observation_names(matches)
        records = [
            (str(row.obsname), str(kind).upper(), (int(row.layer), int(row.cell)))
            for row in matches.itertuples(index=False)
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str = "HEAD",
        package=None,
    ):
        """Attach an MF6 utility observation package for these head targets."""

        import flopy

        owner = model.gwf if package is None else package
        csv_name = _default_obs_csv_name(filename, "head_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(model, csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

    def calibration_plot(self, model, *, type: str = "calibration"):
        """Build a calibration plot directly from these targets and one model."""

        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        return CalibrationPlot.from_targets(self, model=model, type=type)


@dataclass
class LakeStageTargets:
    """Named lake-stage observation definitions independent of PEST setup.

    Parameters
    ----------
    series
        Mapping of user-facing lake series names to zero-based MF6 lake ids.
    """

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "stage"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        normalized: dict[str, int] = {}
        source = self.locations if self.locations is not None else {}
        if isinstance(source, pd.Series):
            source = source.to_dict()
        elif isinstance(source, pd.DataFrame):
            if {"name", "lake"} - set(source.columns):
                raise ValueError("LakeStageTargets locations DataFrame requires 'name' and 'lake' columns.")
            source = source.set_index("name")["lake"].to_dict()
        elif isinstance(source, (list, tuple)):
            source_frame = pd.DataFrame(source)
            if {"name", "lake"} - set(source_frame.columns):
                raise ValueError("LakeStageTargets locations records require 'name' and 'lake' keys.")
            source = source_frame.set_index("name")["lake"].to_dict()
        for name, lake_id in dict(source).items():
            normalized[str(name)] = int(lake_id)
        self._series = normalized
        self._values = _normalize_lake_stage_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            times=self.times,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._series) != 1:
                raise ValueError(
                    "LakeStageTargets values without explicit names require exactly one named lake target."
                )
            only_name = next(iter(self._series))
            self._values["name"] = only_name

    @classmethod
    def from_series(
        cls,
        *,
        lake: str,
        lake_id: int,
        values,
        times=None,
        time_column: str = "time",
        value_column: str = "stage",
    ):
        """Create one named lake-stage target from a Series/list-like input."""

        if isinstance(values, pd.Series):
            series = values.copy()
            series.index.name = time_column
            frame = pd.DataFrame(
                {
                    time_column: series.index.to_list(),
                    "name": [lake] * len(series),
                    value_column: series.to_list(),
                }
            )
            return cls(locations={lake: lake_id}, values=frame, time_column=time_column, value_column=value_column)
        if times is None:
            times = list(range(len(values)))
        frame = pd.DataFrame({time_column: list(times), lake: list(values)})
        return cls(locations={lake: lake_id}, values=frame, time_column=time_column, value_column=value_column)

    @classmethod
    def from_records(
        cls,
        records,
        *,
        time_column: str = "time",
        value_column: str = "stage",
    ):
        """Create lake-stage targets from long-format records."""

        frame = pd.DataFrame(records).copy()
        if frame.empty:
            raise ValueError("from_records(...) requires at least one record.")
        if {"name", "lake", time_column, value_column} - set(frame.columns):
            missing = {"name", "lake", time_column, value_column} - set(frame.columns)
            raise ValueError("from_records(...) is missing required columns: " + ", ".join(sorted(missing)))
        locations = (
            frame.loc[:, ["name", "lake"]]
            .drop_duplicates(subset=["name"])
            .set_index("name")["lake"]
            .to_dict()
        )
        values_frame = frame.loc[:, [time_column, "name", value_column]].copy()
        return cls(locations=locations, values=values_frame, time_column=time_column, value_column=value_column)

    def get(self) -> pd.DataFrame:
        """Return the named lake-stage definitions or merged values."""

        definitions = pd.DataFrame(
            [{"name": name, "lake": lake_id} for name, lake_id in self._series.items()]
        )
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions, on="name", how="left"))

    def to_long(self) -> pd.DataFrame:
        """Return the long-format lake-stage target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per lake series."""

        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._series.keys()])
        frame = self.get().pivot_table(
            index="time",
            columns="name",
            values="stage",
            aggfunc="first",
        ).reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the lake-stage target set."""

        return pd.DataFrame(
            [{"n_lakes": int(len(self._series)), "n_rows": int(len(self._values))}]
        )

    def _obs_column_map(self) -> dict[str, str]:
        return {
            str(name): _observation_label(name)
            for name in self._series.keys()
        }

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-stage table indexed by stress period."""

        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            return observed
        frame = model.packages.lak.results.stage.get().copy()
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._series.keys()])
        frame = _aggregate_named_series(frame, id_column="lake", value_column="stage")
        definitions = pd.DataFrame([{"name": name, "lake": lake_id} for name, lake_id in self._series.items()])
        merged = definitions.merge(frame, on="lake", how="left")
        wide = merged.pivot_table(index="per", columns="name", values="stage", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target lake stages against simulated lake stages."""

        definitions = pd.DataFrame([{"name": name, "lake": lake_id} for name, lake_id in self._series.items()])
        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = definitions.merge(simulated, on="name", how="left")
        else:
            simulated = model.packages.lak.results.stage.get().copy()
            simulated = _aggregate_named_series(simulated, id_column="lake", value_column="stage")
            simulated = definitions.merge(simulated, on="lake", how="left").rename(columns={"stage": "sim_stage"})
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["stage_target"] = np.nan
        else:
            targets = self.get().copy()
            targets["per"] = _numeric_periods(targets, caller="LakeStageTargets.compare(...)")
            if "stage" in targets.columns and "stage_target" not in targets.columns:
                targets = targets.rename(columns={"stage": "stage_target"})
            compare = targets.merge(simulated, on=["name", "lake", "per"], how="left")
        compare["residual"] = compare["sim_stage"] - compare["stage_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table for lake-stage targets."""

        frame = self.compare(model)
        return frame.loc[:, ["name", "lake", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the lake-stage target set."""

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 lake-stage observations."""

        records = [
            (_observation_label(name), str(kind).upper(), int(lake_id) + 1)
            for name, lake_id in self._series.items()
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "lak_obs",
        filename: str = "lake_stage_targets.obs",
        csv_name: str | None = None,
        kind: str = "STAGE",
        package=None,
    ):
        """Attach an MF6 utility observation package for lake-stage targets."""

        import flopy

        owner = model.gwf.lak if package is None else package
        csv_name = _default_obs_csv_name(filename, "lake_stage_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

    def calibration_plot(self, model, *, type: str = "calibration"):
        """Build a calibration plot directly from these targets and one model."""

        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        compare = self.compare(model)
        return CalibrationPlot.from_compare(
            compare,
            type="obs_vs_sim" if type == "calibration" else type,
            target_column="stage_target",
            simulated_column="sim_stage",
            title="Observed vs simulated lake stage" if type != "heads" else "Lake stage",
            yaxis_title="Simulated stage",
        )


@dataclass
class SfrStageTargets:
    """Named SFR stage targets keyed by zero-based reach number."""

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "stage"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        self._locations = _normalize_named_integer_locations(self.locations, id_column="reach", target_name="SfrStageTargets")
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="stage_target",
            times=self.times,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations) != 1:
                raise ValueError("SfrStageTargets values without explicit names require exactly one target reach.")
            self._values["name"] = str(self._locations["name"].iloc[0])

    def get(self) -> pd.DataFrame:
        definitions = self._locations.copy()
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions, on="name", how="left"))

    def to_long(self) -> pd.DataFrame:
        return self.get()

    def to_wide(self) -> pd.DataFrame:
        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = self.get().pivot_table(index="time", columns="name", values="stage_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame([{"n_reaches": int(len(self._locations)), "n_rows": int(len(self._values))}])

    def _obs_column_map(self) -> dict[str, str]:
        return {
            str(name): _observation_label(name)
            for name in self._locations["name"].astype(str).tolist()
        }

    def simulated_series(self, model) -> pd.DataFrame:
        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            return observed
        frame = model.packages.sfr.results.stage.get().copy()
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = _aggregate_named_series(frame, id_column="reach", value_column="stage")
        merged = self._locations.merge(frame, on="reach", how="left")
        wide = merged.pivot_table(index="per", columns="name", values="stage", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = self._locations.merge(simulated, on="name", how="left")
        else:
            simulated = _aggregate_named_series(model.packages.sfr.results.stage.get().copy(), id_column="reach", value_column="stage")
            simulated = self._locations.merge(simulated, on="reach", how="left").rename(columns={"stage": "sim_stage"})
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["stage_target"] = np.nan
        else:
            targets = self.get().copy()
            targets["per"] = _numeric_periods(targets, caller="SfrStageTargets.compare(...)")
            compare = targets.merge(simulated, on=["name", "reach", "per"], how="left")
        compare["residual"] = compare["sim_stage"] - compare["stage_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        frame = self.compare(model)
        return frame.loc[:, ["name", "reach", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "sfr_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, int]]]:
        records = [
            (_observation_label(name), str(kind).upper(), int(reach) + 1)
            for name, reach in self._locations[["name", "reach"]].itertuples(index=False)
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "sfr_stage_obs",
        filename: str = "sfr_stage_targets.obs",
        csv_name: str | None = None,
        kind: str = "STAGE",
        package=None,
    ):
        import flopy

        owner = model.gwf.sfr if package is None else package
        csv_name = _default_obs_csv_name(filename, "sfr_stage_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )


@dataclass
class SfrFlowTargets:
    """Named SFR flow targets keyed by zero-based reach number."""

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "flow"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        self._locations = _normalize_named_integer_locations(self.locations, id_column="reach", target_name="SfrFlowTargets")
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="flow_target",
            times=self.times,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations) != 1:
                raise ValueError("SfrFlowTargets values without explicit names require exactly one target reach.")
            self._values["name"] = str(self._locations["name"].iloc[0])

    def get(self) -> pd.DataFrame:
        definitions = self._locations.copy()
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions, on="name", how="left"))

    def to_long(self) -> pd.DataFrame:
        return self.get()

    def to_wide(self) -> pd.DataFrame:
        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = self.get().pivot_table(index="time", columns="name", values="flow_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        return pd.DataFrame([{"n_reaches": int(len(self._locations)), "n_rows": int(len(self._values))}])

    def _obs_column_map(self) -> dict[str, str]:
        return {
            str(name): _observation_label(name)
            for name in self._locations["name"].astype(str).tolist()
        }

    def _package_flow_table(self, model) -> pd.DataFrame:
        if not hasattr(model, "outputs"):
            try:
                flow = model.packages.sfr.results.q.get().copy()
            except Exception:
                return pd.DataFrame(columns=["per", "reach", "sim_flow"])
            if flow.empty:
                return pd.DataFrame(columns=["per", "reach", "sim_flow"])
            flow["per"] = pd.to_numeric(flow["per"], errors="coerce").astype(int)
            flow["reach"] = pd.to_numeric(flow["reach"], errors="coerce").astype(int)
            flow["sim_flow"] = _clean_observation_output_values(flow["q"]).abs()
            return flow.loc[:, ["per", "reach", "sim_flow"]]

        flow = model.outputs.sfr.bud.get("FLOW-JA-FACE")
        if not isinstance(flow, pd.DataFrame) or flow.empty:
            try:
                fallback = model.packages.sfr.results.q.get().copy()
            except Exception:
                return pd.DataFrame(columns=["per", "reach", "sim_flow"])
            if fallback.empty:
                return pd.DataFrame(columns=["per", "reach", "sim_flow"])
            fallback["per"] = pd.to_numeric(fallback["per"], errors="coerce").astype(int)
            fallback["reach"] = pd.to_numeric(fallback["reach"], errors="coerce").astype(int)
            fallback["sim_flow"] = _clean_observation_output_values(fallback["q"]).abs()
            return fallback.loc[:, ["per", "reach", "sim_flow"]]

        frame = flow.copy()
        frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        frame["reach"] = pd.to_numeric(frame["node"], errors="coerce")
        frame["reach_to"] = pd.to_numeric(frame["node2"], errors="coerce")
        if frame["reach"].dropna().min() >= 1:
            frame["reach"] = frame["reach"] - 1
        if frame["reach_to"].dropna().min() >= 1:
            frame["reach_to"] = frame["reach_to"] - 1
        frame["q"] = _clean_observation_output_values(frame["q"])

        # For in-channel through-flow, use the outbound FLOW-JA-FACE rows
        # leaving each reach. Those are the negative rows in the SFR package
        # budget. Taking absolute values keeps the reported flow intuitive.
        outbound = frame.loc[frame["q"] <= 0.0, ["per", "reach", "q"]].copy()
        outbound["sim_flow"] = outbound["q"].abs()
        grouped = outbound.groupby(["per", "reach"], as_index=False)["sim_flow"].sum()

        to_mvr = model.outputs.sfr.bud.get("TO-MVR")
        if isinstance(to_mvr, pd.DataFrame) and not to_mvr.empty:
            mvr = to_mvr.copy()
            mvr["per"] = mvr["kstpkper"].apply(lambda values: int(values[1]))
            mvr["reach"] = pd.to_numeric(mvr["node"], errors="coerce")
            if mvr["reach"].dropna().min() >= 1:
                mvr["reach"] = mvr["reach"] - 1
            mvr["sim_flow"] = _clean_observation_output_values(mvr["q"]).abs()
            grouped = (
                pd.concat([grouped, mvr.loc[:, ["per", "reach", "sim_flow"]]], ignore_index=True)
                .groupby(["per", "reach"], as_index=False)["sim_flow"]
                .sum()
            )

        grouped["reach"] = grouped["reach"].astype(int)
        grouped["per"] = grouped["per"].astype(int)
        return grouped

    def simulated_series(self, model) -> pd.DataFrame:
        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            return observed

        frame = self._package_flow_table(model)
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].astype(str).tolist()])
        merged = self._locations.merge(frame, on="reach", how="left")
        value_column = "sim_flow" if "sim_flow" in merged.columns else "q"
        wide = merged.pivot_table(index="per", columns="name", values=value_column, aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_flow")
            simulated["per"] = simulated["time"]
            simulated = self._locations.merge(simulated, on="name", how="left")
        else:
            simulated = self._package_flow_table(model)
            simulated = self._locations.merge(simulated, on="reach", how="left")
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["flow_target"] = np.nan
        else:
            targets = self.get().copy()
            targets["per"] = _numeric_periods(targets, caller="SfrFlowTargets.compare(...)")
            compare = targets.merge(simulated, on=["name", "reach", "per"], how="left")
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        frame = self.compare(model)
        return frame.loc[:, ["name", "reach", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "sfr_flow_targets.csv",
        kind: str = "DOWNSTREAM-FLOW",
    ) -> dict[str, list[tuple[str, str, int]]]:
        records = [
            (_observation_label(name), str(kind).upper(), int(reach) + 1)
            for name, reach in self._locations[["name", "reach"]].itertuples(index=False)
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "sfr_flow_obs",
        filename: str = "sfr_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "DOWNSTREAM-FLOW",
        package=None,
    ):
        import flopy

        owner = model.gwf.sfr if package is None else package
        csv_name = _default_obs_csv_name(filename, "sfr_flow_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )


def _normalize_drn_zone_locations(
    locations,
    *,
    name_column: str = "name",
    group_column: str | None = "group",
    weight_column: str | None = "weight",
) -> pd.DataFrame:
    """Normalize drain-zone definitions from cells or polygons."""

    if isinstance(locations, pd.Series):
        frame = locations.rename("cells").rename_axis(name_column).reset_index()
        if frame.columns[0] != name_column:
            frame = frame.rename(columns={frame.columns[0]: name_column})
    elif isinstance(locations, dict):
        lower_keys = {str(key).strip().lower() for key in locations.keys()}
        if "name" in lower_keys or "cells" in lower_keys or "geometry" in lower_keys:
            frame = pd.DataFrame(locations)
        else:
            frame = pd.DataFrame({name_column: list(locations.keys()), "cells": list(locations.values())})
    elif isinstance(locations, (list, tuple)):
        frame = pd.DataFrame(locations)
    elif isinstance(locations, (pd.DataFrame, gpd.GeoDataFrame)):
        frame = locations.copy()
    else:
        frame = gpd.read_file(locations)

    if name_column not in frame.columns:
        raise ValueError("DrnFlowTargets locations must include a name column.")
    keep = [name_column]
    for column in ("cells", group_column, weight_column, "layer"):
        if column is not None and column in frame.columns:
            keep.append(column)
    if isinstance(frame, gpd.GeoDataFrame):
        keep.append(frame.geometry.name)
    frame = frame.loc[:, list(dict.fromkeys(keep))].copy()
    frame = frame.rename(columns={name_column: "name"})
    if group_column and group_column in frame.columns:
        frame = frame.rename(columns={group_column: "group"})
    if weight_column and weight_column in frame.columns:
        frame = frame.rename(columns={weight_column: "weight"})
    if "group" not in frame.columns:
        frame["group"] = pd.NA
    if "weight" not in frame.columns:
        frame["weight"] = np.nan
    if "layer" not in frame.columns:
        frame["layer"] = 0
    if "cells" in frame.columns:
        def _coerce_cells(values):
            if isinstance(values, str):
                return [int(value) for value in values.split(",") if str(value).strip() != ""]
            if isinstance(values, (list, tuple, np.ndarray, pd.Series)):
                return [int(value) for value in values]
            if pd.isna(values):
                return []
            return [int(values)]

        frame["cells"] = frame["cells"].apply(_coerce_cells)
    return frame


def _resolve_drn_zone_cells(locations: pd.DataFrame, model) -> pd.DataFrame:
    """Attach explicit cell lists to each DRN zone, intersecting polygons when needed."""

    frame = locations.copy()
    if "cells" in frame.columns and frame["cells"].notna().all():
        return frame
    if not isinstance(frame, gpd.GeoDataFrame):
        raise ValueError("DrnFlowTargets polygon zones require a GeoDataFrame or GIS layer input.")
    vor = model.vor.gdf_vorPolys.reset_index().rename(columns={"index": "cell"})
    joined = gpd.sjoin(frame, vor.loc[:, ["cell", vor.geometry.name]], how="left", predicate="intersects")
    cells = (
        joined.groupby("name", dropna=False)["cell"]
        .apply(lambda series: sorted({int(value) for value in series.dropna().tolist()}))
        .rename("cells")
        .reset_index()
    )
    merged = frame.drop(columns=["cells"], errors="ignore").merge(cells, on="name", how="left")
    merged["cells"] = merged["cells"].apply(lambda value: [] if not isinstance(value, list) else value)
    return merged


@dataclass
class DrnFlowTargets:
    """Named seepage-zone targets aggregated from DRN budget cells."""

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    name_column: str = "name"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "flow"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        self._locations = _normalize_drn_zone_locations(
            self.locations,
            name_column=self.name_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
        )
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="flow_target",
            times=self.times,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations["name"].unique()) != 1:
                raise ValueError("DrnFlowTargets values without explicit names require exactly one zone.")
            self._values["name"] = str(self._locations["name"].iloc[0])

    @property
    def locations_gdf(self):
        return self._locations.copy()

    def zone_definitions(self, model) -> pd.DataFrame:
        return _resolve_drn_zone_cells(self._locations, model)

    def get(self, model=None) -> pd.DataFrame:
        definitions = self._locations.copy()
        if model is not None:
            definitions = self.zone_definitions(model)
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions.drop(columns="geometry", errors="ignore"), on="name", how="left"))

    def to_long(self, model=None) -> pd.DataFrame:
        return self.get(model=model)

    def to_wide(self) -> pd.DataFrame:
        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].astype(str).tolist()])
        frame = self._values.pivot_table(index="time", columns="name", values="flow_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self, model=None) -> pd.DataFrame:
        definitions = self.zone_definitions(model) if model is not None else self._locations.copy()
        zone_count = int(definitions["name"].nunique())
        cell_count = int(sum(len(value) for value in definitions.get("cells", pd.Series(dtype=object))))
        return pd.DataFrame([{"n_zones": zone_count, "n_cells": cell_count, "n_rows": int(len(self._values))}])

    def _drn_budget_frame(self, model) -> pd.DataFrame:
        budget = model.bud("drn").df.reset_index().copy()
        if "per" not in budget.columns:
            if "kstpkper" not in budget.columns:
                raise ValueError("DRN budget dataframe requires 'per' or 'kstpkper' columns.")
            budget["per"] = budget["kstpkper"].apply(lambda item: int(item[1]))
        if "node" not in budget.columns:
            raise ValueError("DRN budget dataframe requires a zero-based 'node' column.")
        budget["node"] = pd.to_numeric(budget["node"], errors="coerce").astype(int)
        budget["q"] = pd.to_numeric(budget["q"], errors="coerce")
        return budget

    def _obs_column_map(self, model) -> tuple[pd.DataFrame, dict[str, list[str]]]:
        zones = self.zone_definitions(model)
        zone_columns: dict[str, list[str]] = {}
        for row in zones.itertuples(index=False):
            zone_columns[str(row.name)] = [
                _observation_label(f"{row.name}_c{int(cell)}")
                for cell in getattr(row, "cells", [])
            ]
        return zones, zone_columns

    def _simulated_series_from_obs(self, model) -> pd.DataFrame | None:
        zones, zone_columns = self._obs_column_map(model)
        expected = {
            f"{zone_name}__{index}": column
            for zone_name, columns in zone_columns.items()
            for index, column in enumerate(columns)
        }
        observed = _find_named_observation_output(model, expected_columns=expected)
        if observed is None:
            return None

        result = pd.DataFrame({"time": observed["time"]})
        for zone_name, columns in zone_columns.items():
            if not columns:
                result[zone_name] = 0.0
                continue
            zone_frame = pd.DataFrame(
                {
                    column: _clean_observation_output_values(
                        observed[f"{zone_name}__{index}"]
                    )
                    for index, column in enumerate(columns)
                }
            )
            result[zone_name] = zone_frame.sum(axis=1)
        return result

    def _compare_from_obs(self, model) -> pd.DataFrame | None:
        wide = self._simulated_series_from_obs(model)
        if wide is None:
            return None

        zones = self.zone_definitions(model)
        simulated = wide.melt(id_vars=["time"], var_name="name", value_name="sim_flow")
        simulated["per"] = simulated["time"]
        simulated = zones.merge(simulated, on="name", how="left")
        if self._values.empty:
            compare = simulated.copy()
            compare["flow_target"] = np.nan
        else:
            targets = self._values.copy()
            targets["per"] = _numeric_periods(targets, caller="DrnFlowTargets.compare(...)")
            compare = targets.merge(simulated, on=["name", "per"], how="left")
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def simulated_series(self, model) -> pd.DataFrame:
        observed = self._simulated_series_from_obs(model)
        if observed is not None:
            return observed

        zones = self.zone_definitions(model)
        budget = self._drn_budget_frame(model)
        rows = []
        for row in zones.itertuples(index=False):
            zone_cells = {int(value) for value in getattr(row, "cells", [])}
            if not zone_cells:
                continue
            subset = budget.loc[budget["node"].isin(zone_cells)].copy()
            grouped = subset.groupby("per", as_index=False)["q"].sum()
            grouped["name"] = str(row.name)
            rows.append(grouped)
        if not rows:
            return pd.DataFrame(columns=["time", *zones["name"].astype(str).tolist()])
        frame = pd.concat(rows, ignore_index=True)
        wide = frame.pivot_table(index="per", columns="name", values="q", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        observed = self._compare_from_obs(model)
        if observed is not None:
            return observed

        zones = self.zone_definitions(model)
        budget = self._drn_budget_frame(model)
        rows = []
        for row in zones.itertuples(index=False):
            zone_cells = {int(value) for value in getattr(row, "cells", [])}
            subset = budget.loc[budget["node"].isin(zone_cells)].copy()
            grouped = subset.groupby("per", as_index=False)["q"].sum()
            grouped["name"] = str(row.name)
            grouped["group"] = getattr(row, "group", pd.NA)
            grouped["weight"] = getattr(row, "weight", np.nan)
            grouped["cells"] = [list(zone_cells)] * len(grouped)
            rows.append(grouped)
        if rows:
            simulated = pd.concat(rows, ignore_index=True).rename(columns={"q": "sim_flow"})
        else:
            simulated = pd.DataFrame(columns=["per", "name", "sim_flow", "group", "weight", "cells"])
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["flow_target"] = np.nan
        else:
            targets = self._values.copy()
            targets["per"] = _numeric_periods(targets, caller="DrnFlowTargets.compare(...)")
            compare = targets.merge(simulated, on=["name", "per"], how="left")
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        frame = self.compare(model)
        columns = ["name", "group", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual", "weight"]
        if "cells" in frame.columns:
            columns.append("cells")
        return frame.loc[:, columns]

    def stats(self, model) -> pd.DataFrame:
        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "drn_flow_targets.csv",
        kind: str = "DRN",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        zones = self.zone_definitions(model)
        records: list[tuple[str, str, tuple[int, int]]] = []
        for row in zones.itertuples(index=False):
            for cell in getattr(row, "cells", []):
                obsname = _observation_label(f"{row.name}_c{int(cell)}")
                records.append((obsname, str(kind).upper(), (int(getattr(row, "layer", 0)), int(cell))))
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "drn_flow_obs",
        filename: str = "drn_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "DRN",
        package=None,
    ):
        import flopy

        owner = model.gwf.drn if package is None else package
        csv_name = _default_obs_csv_name(filename, "drn_flow_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(model, csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

class BoundHeadTargets:
    """Model-bound head-target helper returned from ``model.targets.heads``."""

    def __init__(self, model, targets: HeadTargets):
        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        return self.targets.summary()

    def plot_locations(self, *, ax=None, **kwargs):
        return self.targets.plot_locations(ax=ax, **kwargs)

    def match_to_model(self, model=None) -> pd.DataFrame:
        return self.targets.match_to_model(self.model if model is None else model)

    def simulated_heads(self, model=None) -> pd.DataFrame:
        return self.targets.simulated_heads(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ):
        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str = "HEAD",
        package=None,
    ):
        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        """Plotting helpers for this model-bound target set."""

        if self._plot is None:
            self._plot = BoundHeadTargetPlots(self)
        return self._plot


class BoundHeadTargetPlots:
    """Plotting facade for model-bound head targets."""

    def __init__(self, bound_targets: BoundHeadTargets):
        self.bound_targets = bound_targets

    @property
    def model(self):
        return self.bound_targets.model

    @property
    def targets(self):
        return self.bound_targets.targets

    def locations(self, *, ax=None, **kwargs):
        """Plot target locations on a map axis."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def calibration(self, *, type: str = "calibration"):
        """Return a calibration/heads plot for the bound targets."""

        return self.bound_targets.calibration_plot(type=type)

    def obs_vs_sim(self, *, baseline=None, ax=None):
        """Return a target-vs-simulated cross plot for one or two models."""

        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(current, baseline_compare=baseline_frame)

    def timeseries(self, name: str | None = None, *, baseline=None, ax=None):
        """Return a time-series plot for one target location."""

        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
        )

    def residuals_by_period(self, *, baseline=None):
        """Return a by-period residual summary plot."""

        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(
            current,
            baseline_compare=baseline_frame,
        )


class BoundLakeStageTargets:
    """Model-bound lake-stage helper returned from ``model.targets.lake_stage``."""

    def __init__(self, model, targets: LakeStageTargets):
        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        return self.targets.summary()

    def simulated_series(self, model=None) -> pd.DataFrame:
        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ):
        return self.targets.to_flopy_obs(csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "lak_obs",
        filename: str = "lake_stage_targets.obs",
        csv_name: str | None = None,
        kind: str = "STAGE",
        package=None,
    ):
        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="stage_target",
                simulated_column="sim_stage",
                title="Observed vs simulated lake stage",
                yaxis_title="Stage",
            )
        return self._plot


class BoundSfrStageTargets:
    """Model-bound SFR stage helper returned from ``model.targets.sfr_stage``."""

    def __init__(self, model, targets: SfrStageTargets):
        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        return self.targets.summary()

    def simulated_series(self, model=None) -> pd.DataFrame:
        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "sfr_stage_targets.csv", kind: str = "STAGE"):
        return self.targets.to_flopy_obs(csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "sfr_stage_obs",
        filename: str = "sfr_stage_targets.obs",
        csv_name: str | None = None,
        kind: str = "STAGE",
        package=None,
    ):
        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    @property
    def plot(self):
        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="stage_target",
                simulated_column="sim_stage",
                title="Observed vs simulated SFR stage",
                yaxis_title="Stage",
            )
        return self._plot


class BoundSfrFlowTargets:
    """Model-bound SFR flow helper returned from ``model.targets.sfr_flow``."""

    def __init__(self, model, targets: SfrFlowTargets):
        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        return self.targets.summary()

    def simulated_series(self, model=None) -> pd.DataFrame:
        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "sfr_flow_targets.csv", kind: str = "DOWNSTREAM-FLOW"):
        return self.targets.to_flopy_obs(csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "sfr_flow_obs",
        filename: str = "sfr_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "DOWNSTREAM-FLOW",
        package=None,
    ):
        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    @property
    def plot(self):
        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="flow_target",
                simulated_column="sim_flow",
                title="Observed vs simulated SFR flow",
                yaxis_title="Flow",
            )
        return self._plot


class BoundDrnFlowTargets:
    """Model-bound DRN seepage-zone helper returned from ``model.targets.drn_flow``."""

    def __init__(self, model, targets: DrnFlowTargets):
        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        return self.targets.get(model=self.model)

    def to_long(self) -> pd.DataFrame:
        return self.targets.to_long(model=self.model)

    def to_wide(self) -> pd.DataFrame:
        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        return self.targets.summary(model=self.model)

    def simulated_series(self, model=None) -> pd.DataFrame:
        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "drn_flow_targets.csv", kind: str = "DRN"):
        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "drn_flow_obs",
        filename: str = "drn_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "DRN",
        package=None,
    ):
        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    @property
    def plot(self):
        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="flow_target",
                simulated_column="sim_flow",
                title="Observed vs simulated DRN seepage",
                yaxis_title="Flow",
            )
        return self._plot


class BoundNamedSeriesTargetPlots:
    """Shared plotting facade for model-bound stage/flow target sets."""

    def __init__(
        self,
        bound_targets,
        *,
        target_column: str,
        simulated_column: str,
        title: str,
        yaxis_title: str,
    ):
        self.bound_targets = bound_targets
        self.target_column = target_column
        self.simulated_column = simulated_column
        self.title = title
        self.yaxis_title = yaxis_title

    def obs_vs_sim(self, *, baseline=None):
        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(
            current,
            baseline_compare=baseline_frame,
            target_column=self.target_column,
            simulated_column=self.simulated_column,
            title=self.title,
            yaxis_title=self.yaxis_title,
        )

    def timeseries(self, name: str | None = None, *, baseline=None):
        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
            target_column=self.target_column,
            simulated_column=self.simulated_column,
            yaxis_title=self.yaxis_title,
            title=None,
        )

    def residuals_by_period(self, *, baseline=None):
        from simple_modflow.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(current, baseline_compare=baseline_frame)


class TargetRegistry:
    """Lightweight model-bound registry for reusable calibration targets."""

    _INTERNAL_NAMES = {"model", "_targets"}

    def __init__(self, model):
        object.__setattr__(self, "model", model)
        object.__setattr__(self, "_targets", {})

    def _bind(self, target):
        if isinstance(target, HeadTargets):
            return BoundHeadTargets(self.model, target)
        if isinstance(target, LakeStageTargets):
            return BoundLakeStageTargets(self.model, target)
        if isinstance(target, SfrStageTargets):
            return BoundSfrStageTargets(self.model, target)
        if isinstance(target, SfrFlowTargets):
            return BoundSfrFlowTargets(self.model, target)
        if isinstance(target, DrnFlowTargets):
            return BoundDrnFlowTargets(self.model, target)
        return target

    def _coerce_target(self, value):
        if isinstance(
            value,
            (
                BoundHeadTargets,
                BoundLakeStageTargets,
                BoundSfrStageTargets,
                BoundSfrFlowTargets,
                BoundDrnFlowTargets,
            ),
        ):
            return value.targets
        if isinstance(value, (HeadTargets, LakeStageTargets, SfrStageTargets, SfrFlowTargets, DrnFlowTargets)):
            return value
        raise TypeError(
            "Model targets currently support HeadTargets, LakeStageTargets, "
            "SfrStageTargets, SfrFlowTargets, and DrnFlowTargets."
        )

    def keys(self) -> list[str]:
        return sorted(self._targets)

    def summary(self) -> pd.DataFrame:
        rows = [
            {"name": name, "type": type(target).__name__}
            for name, target in sorted(self._targets.items())
        ]
        return pd.DataFrame(rows)

    def _named_target(self, name: str):
        if name not in self._targets:
            raise AttributeError(f"{type(self).__name__!r} has no target set {name!r}")
        return self[name]

    @property
    def heads(self) -> BoundHeadTargets:
        """Return model-bound head targets with IDE-visible completion."""

        return self._named_target("heads")

    @heads.setter
    def heads(self, value: HeadTargets | BoundHeadTargets):
        self["heads"] = value

    @property
    def lake_stage(self) -> BoundLakeStageTargets:
        """Return model-bound lake-stage targets."""

        return self._named_target("lake_stage")

    @lake_stage.setter
    def lake_stage(self, value: LakeStageTargets | BoundLakeStageTargets):
        self["lake_stage"] = value

    @property
    def sfr_stage(self) -> BoundSfrStageTargets:
        """Return model-bound SFR-stage targets."""

        return self._named_target("sfr_stage")

    @sfr_stage.setter
    def sfr_stage(self, value: SfrStageTargets | BoundSfrStageTargets):
        self["sfr_stage"] = value

    @property
    def sfr_flow(self) -> BoundSfrFlowTargets:
        """Return model-bound SFR-flow targets."""

        return self._named_target("sfr_flow")

    @sfr_flow.setter
    def sfr_flow(self, value: SfrFlowTargets | BoundSfrFlowTargets):
        self["sfr_flow"] = value

    @property
    def drn_flow(self) -> BoundDrnFlowTargets:
        """Return model-bound DRN seepage-flow targets."""

        return self._named_target("drn_flow")

    @drn_flow.setter
    def drn_flow(self, value: DrnFlowTargets | BoundDrnFlowTargets):
        self["drn_flow"] = value

    def __getitem__(self, key: str):
        return self._bind(self._targets[key])

    def __setitem__(self, key: str, value):
        self._targets[str(key)] = self._coerce_target(value)

    def __contains__(self, key: str) -> bool:
        return str(key) in self._targets

    def __delitem__(self, key: str):
        del self._targets[str(key)]

    def __getattr__(self, name: str):
        targets = object.__getattribute__(self, "__dict__").get("_targets", {})
        if name in targets:
            return self._bind(targets[name])
        raise AttributeError(f"{type(self).__name__!r} has no target set {name!r}")

    def __setattr__(self, name: str, value):
        if name in self._INTERNAL_NAMES:
            object.__setattr__(self, name, value)
            return
        targets = object.__getattribute__(self, "__dict__").get("_targets")
        if targets is None:
            object.__setattr__(self, name, value)
            return
        targets[name] = self._coerce_target(value)

    def __dir__(self):
        targets = object.__getattribute__(self, "__dict__").get("_targets", {})
        return sorted(set(super().__dir__()) | set(targets))
