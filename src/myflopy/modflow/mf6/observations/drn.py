"""DRN flow observation targets (`DrnFlowTargets` + bound view)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _clean_observation_output_values,
    _copy_frame,
    _default_compare_stats,
    _default_obs_csv_name,
    _find_named_observation_output,
    _normalize_compare_time,
    _normalize_named_series_values,
    _numeric_periods,
    _observation_label,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots


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
            """Coerce one zone's cell spec (CSV string / sequence / scalar / NaN) to a list of ints."""

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
    """Observed seepage/discharge targets aggregated from DRN budget cells.

    Sums the simulated DRN flow over a *group* of drain cells to form one flow
    series per named seepage zone, then aligns it with observed discharge -- the
    pattern for calibrating spring/seep or drain-discharge measurements where the
    field datum is a lumped flux over many cells. Unlike the SFR/LAK targets,
    ``locations`` maps each ``name`` to a set of drain cells (via geometry or an
    explicit cell list), and the comparison aggregates the DRN cell-by-cell
    budget. Pass an instance to ``cal.observe(...)`` / ``cal.forecast(...)``, or
    call ``.compare(model)`` / ``.stats(model)`` for standalone review.

    Parameters
    ----------
    locations
        Seepage-zone definitions mapping ``name`` to its member drain cells --
        a GeoDataFrame of points/polygons to spatially match, or a table with a
        ``name`` column plus ``cell``/``node`` ids.
    values
        Observed discharge table (long or wide). Optional (capture template).
    name_column, group_column, weight_column
        Column names in ``locations`` for the zone id, observation group, and
        observation weight.
    time_column, value_column
        Column names in ``values`` for the period/time index and the flow value.
    times
        Optional explicit period/time labels when ``values`` carries none.

    Examples
    --------
    >>> DrnFlowTargets(locations=seep_zones_gdf,
    ...                values=observed_seepage_df)       # -> cal.observe(...)

    See Also
    --------
    SfrFlowTargets, HeadTargets, LakeStageTargets
    """

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    name_column: str = "name"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "flow"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Normalize the zone definitions and observed discharge into internal tables."""

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
        """A copy of the normalized zone-definition table (geometry preserved)."""

        return self._locations.copy()

    def zone_definitions(self, model) -> pd.DataFrame:
        """Resolve each zone's member cells against ``model`` (intersecting polygons if needed)."""

        return _resolve_drn_zone_cells(self._locations, model)

    def get(self, model=None) -> pd.DataFrame:
        """Return the zone definitions (cells resolved when ``model`` is given), merged with values."""

        definitions = self._locations.copy()
        if model is not None:
            definitions = self.zone_definitions(model)
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions.drop(columns="geometry", errors="ignore"), on="name", how="left"))

    def to_long(self, model=None) -> pd.DataFrame:
        """Return the long-format DRN-flow target table."""

        return self.get(model=model)

    def to_wide(self) -> pd.DataFrame:
        """Return the observed discharge in wide format with one column per zone."""

        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].astype(str).tolist()])
        frame = self._values.pivot_table(index="time", columns="name", values="flow_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self, model=None) -> pd.DataFrame:
        """Return a compact summary (zone count, resolved cell count, row count)."""

        definitions = self.zone_definitions(model) if model is not None else self._locations.copy()
        zone_count = int(definitions["name"].nunique())
        cell_count = int(sum(len(value) for value in definitions.get("cells", pd.Series(dtype=object))))
        return pd.DataFrame([{"n_zones": zone_count, "n_cells": cell_count, "n_rows": int(len(self._values))}])

    def _drn_budget_frame(self, model) -> pd.DataFrame:
        """The DRN cell-by-cell budget as a tidy per/node/q frame (adds ``per`` if only kstpkper)."""

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
        """Return the resolved zones and, per zone, its per-cell MF6 observation labels."""

        zones = self.zone_definitions(model)
        zone_columns: dict[str, list[str]] = {}
        for row in zones.itertuples(index=False):
            zone_columns[str(row.name)] = [
                _observation_label(f"{row.name}_c{int(cell)}")
                for cell in getattr(row, "cells", [])
            ]
        return zones, zone_columns

    def _simulated_series_from_obs(self, model) -> pd.DataFrame | None:
        """Sum per-cell observation-CSV outputs into one wide series per zone, or ``None`` if absent."""

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
        """Build the target-vs-simulated compare table from observation output, or ``None``."""

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
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-discharge table per zone, summed over each zone's DRN cells."""

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
        """Compare observed discharge against simulated per-zone DRN flow (with residuals)."""

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
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table (target/sim/residual, with group/weight) for DRN zones."""

        frame = self.compare(model)
        columns = ["name", "group", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual", "weight"]
        if "cells" in frame.columns:
            columns.append("cells")
        return frame.loc[:, columns]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the DRN-flow target set."""

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "drn_flow_targets.csv",
        kind: str = "DRN",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        """Return a FloPy ``continuous`` dict with one MF6 DRN observation per member cell."""

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
        """Attach an MF6 utility observation package to DRN for these seepage-zone targets."""

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

class BoundDrnFlowTargets:
    """Model-bound DRN seepage-zone helper returned from ``model.targets.drn_flow``."""

    def __init__(self, model, targets: DrnFlowTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The zone definitions with member cells resolved against the bound model."""

        return self.targets.get(model=self.model)

    def to_long(self) -> pd.DataFrame:
        """The long-format DRN-flow target table (resolved against the bound model)."""

        return self.targets.to_long(model=self.model)

    def to_wide(self) -> pd.DataFrame:
        """The observed discharge in wide format (one column per zone)."""

        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        """A compact summary of the seepage-zone target set (resolved against the bound model)."""

        return self.targets.summary(model=self.model)

    def simulated_series(self, model=None) -> pd.DataFrame:
        """The wide simulated-discharge series per zone (defaults to the bound model)."""

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated DRN discharge (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the seepage-zone targets (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "drn_flow_targets.csv", kind: str = "DRN"):
        """A FloPy ``continuous`` observation dict for the DRN seepage-zone targets."""

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
        """Attach an MF6 DRN observation package to the bound model for these zones."""

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
        """Plotting helpers (obs-vs-sim, time series) for these DRN seepage targets."""

        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="flow_target",
                simulated_column="sim_flow",
                title="Observed vs simulated DRN seepage",
                yaxis_title="Flow",
            )
        return self._plot


