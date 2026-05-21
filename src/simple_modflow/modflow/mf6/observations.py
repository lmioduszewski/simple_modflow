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

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the lake-stage target set."""

        return pd.DataFrame(
            [{"n_lakes": int(len(self._series)), "n_rows": int(len(self._values))}]
        )

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, tuple[int]]]]:
        """Return a FloPy ``continuous`` dict for MF6 lake-stage observations."""

        records = [
            (_observation_label(name), str(kind).upper(), (int(lake_id),))
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
        """Plot target heads against simulated heads for one or two models."""

        current = self.targets.compare(self.model)
        if ax is None:
            _, ax = plt.subplots()
        ax.scatter(current["head_target"], current["sim_head"], alpha=0.8, label="Model")
        if baseline is not None:
            baseline_frame = self.targets.compare(baseline)
            ax.scatter(
                baseline_frame["head_target"],
                baseline_frame["sim_head"],
                alpha=0.8,
                label="Baseline",
            )
            values = pd.concat(
                [
                    current["head_target"],
                    current["sim_head"],
                    baseline_frame["sim_head"],
                ],
                axis=0,
            ).dropna()
        else:
            values = pd.concat([current["head_target"], current["sim_head"]], axis=0).dropna()
        if not values.empty:
            lower = float(values.min())
            upper = float(values.max())
            ax.plot([lower, upper], [lower, upper], linestyle="--", color="black", linewidth=1)
        ax.set_title("Observed vs simulated heads")
        ax.set_xlabel("Observed head")
        ax.set_ylabel("Simulated head")
        ax.grid(True, alpha=0.3)
        ax.legend()
        return ax

    def timeseries(self, name: str | None = None, *, baseline=None, ax=None):
        """Plot target and simulated heads through time for one observation."""

        current = self.targets.compare(self.model)
        names = current["name"].astype(str)
        if name is None:
            unique_names = sorted(names.unique().tolist())
            if len(unique_names) != 1:
                raise ValueError(
                    "plot.timeseries(...) requires name= when more than one target location is present."
                )
            name = unique_names[0]
        mask = names.str.lower() == str(name).strip().lower()
        current = current.loc[mask].copy()
        if current.empty:
            raise ValueError(f"No target rows found for observation name {name!r}.")
        sort_values = pd.to_numeric(current["time"], errors="coerce")
        if sort_values.notna().all():
            current = current.assign(_sort=sort_values).sort_values("_sort").drop(columns="_sort")
        else:
            current = current.sort_values("time")
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(current["time"], current["head_target"], marker="o", label="Target")
        ax.plot(current["time"], current["sim_head"], marker="o", label="Model")
        if baseline is not None:
            baseline_frame = self.targets.compare(baseline)
            baseline_frame = baseline_frame.loc[
                baseline_frame["name"].astype(str).str.lower() == str(name).strip().lower()
            ].copy()
            baseline_sort = pd.to_numeric(baseline_frame["time"], errors="coerce")
            if baseline_sort.notna().all():
                baseline_frame = baseline_frame.assign(_sort=baseline_sort).sort_values("_sort").drop(columns="_sort")
            else:
                baseline_frame = baseline_frame.sort_values("time")
            ax.plot(baseline_frame["time"], baseline_frame["sim_head"], marker="o", label="Baseline")
        ax.set_title(str(current["name"].iloc[0]))
        ax.set_xlabel("Time / period")
        ax.set_ylabel("Head")
        ax.grid(True, alpha=0.3)
        ax.legend()
        return ax


class BoundLakeStageTargets:
    """Model-bound lake-stage helper returned from ``model.targets.lake_stage``."""

    def __init__(self, model, targets: LakeStageTargets):
        self.model = model
        self.targets = targets

    def get(self) -> pd.DataFrame:
        return self.targets.get()

    def summary(self) -> pd.DataFrame:
        return self.targets.summary()

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
        return target

    def _coerce_target(self, value):
        if isinstance(value, (BoundHeadTargets, BoundLakeStageTargets)):
            return value.targets
        if isinstance(value, (HeadTargets, LakeStageTargets)):
            return value
        raise TypeError(
            "Model targets currently support HeadTargets and LakeStageTargets."
        )

    def keys(self) -> list[str]:
        return sorted(self._targets)

    def summary(self) -> pd.DataFrame:
        rows = [
            {"name": name, "type": type(target).__name__}
            for name, target in sorted(self._targets.items())
        ]
        return pd.DataFrame(rows)

    def __getitem__(self, key: str):
        return self._bind(self._targets[key])

    def __setitem__(self, key: str, value):
        self._targets[str(key)] = self._coerce_target(value)

    def __contains__(self, key: str) -> bool:
        return str(key) in self._targets

    def __delitem__(self, key: str):
        del self._targets[str(key)]

    def __getattr__(self, name: str):
        if name in self._targets:
            return self._bind(self._targets[name])
        raise AttributeError(f"{type(self).__name__!r} has no target set {name!r}")

    def __setattr__(self, name: str, value):
        if name in self._INTERNAL_NAMES:
            object.__setattr__(self, name, value)
            return
        self._targets[name] = self._coerce_target(value)

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(self._targets))
