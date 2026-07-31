"""Shared normalization/series helpers for observation targets."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy._logging import get_logger

logger = get_logger(__name__)


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
        except (OSError, UnicodeDecodeError, pd.errors.EmptyDataError,
                pd.errors.ParserError):
            # This scans a workspace for ANY csv that might be observations, so
            # non-observation csvs are expected: unreadable, empty, wrong
            # encoding, or not comma-shaped at all. Skipping them is the job.
            logger.debug("skipping %s: not a readable observation csv", path)
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


def _simulated_field_by_period(
    model: Any,
    *,
    table_attribute: str = "all_heads",
    store_column: str = "elev",
) -> pd.DataFrame:
    """Build one simulated record per period/layer/cell for a dependent variable.

    Keeps the last available time step within each stress period, because that is
    typically the most useful calibration target when the target data are defined
    per stress period.

    ``table_attribute``/``store_column`` select the field: heads
    (``all_heads``/``elev``) or concentration (``all_conc``/``conc``). The output
    column is ``sim_head`` for BOTH -- it is a frame-internal name meaning "the
    simulated value", and keeping it shared is what lets the calibration-plot and
    PEST machinery in ``calcs/calibration.py`` serve concentration unchanged
    (ledger 102).
    """

    frame = getattr(model, table_attribute).reset_index().copy()
    frame["kstp"] = frame["kstpkper"].apply(lambda item: int(item[0]))
    frame["per"] = frame["kstpkper"].apply(lambda item: int(item[1]))
    frame = frame.sort_values(["per", "kstp"]).drop_duplicates(
        subset=["per", "layer", "cell"],
        keep="last",
    )
    frame = frame.rename(columns={store_column: "sim_head"})
    return frame.loc[:, ["per", "layer", "cell", "sim_head"]].reset_index(drop=True)


def _simulated_heads_by_period(model: Any) -> pd.DataFrame:
    """Simulated heads per period/layer/cell (the heads instance of the above)."""

    return _simulated_field_by_period(model)


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


def _normalize_compare_time(frame: pd.DataFrame) -> pd.DataFrame:
    """Return one canonical ``time`` column after target/result merges."""

    result = frame.copy()
    if "time" not in result.columns:
        candidates = [
            result[column]
            for column in ("time_x", "time_y")
            if column in result.columns
        ]
        if candidates:
            time = candidates[0].copy()
            for candidate in candidates[1:]:
                time = time.combine_first(candidate)
            result["time"] = time
    return result.drop(columns=["time_x", "time_y"], errors="ignore")


