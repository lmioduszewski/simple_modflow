"""Head observation targets (`HeadTargets` + bound views/plots)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _build_values_from_cells,
    _coerce_label_list,
    _coerce_scalar_or_list,
    _copy_frame,
    _default_obs_csv_name,
    _ensure_unique_observation_names,
    _load_locations,
    _normalize_identifier_series,
    _normalize_values_frame,
    _simulated_heads_by_period,
)
from myflopy.viz import mpl_axes


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
            _, ax = mpl_axes()
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

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        return CalibrationPlot.from_targets(self, model=model, type=type)


class BoundHeadTargets:
    """Model-bound head-target helper returned from ``model.targets.heads``."""

    def __init__(self, model, targets: HeadTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The target definitions (delegates to the wrapped :class:`HeadTargets`)."""

        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        """The long-format target table (delegates to the wrapped targets)."""

        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        """The wide-format target table (delegates to the wrapped targets)."""

        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        """A compact summary of the target set (delegates to the wrapped targets)."""

        return self.targets.summary()

    def plot_locations(self, *, ax=None, **kwargs):
        """Plot the target locations on a map axis (delegates to the wrapped targets)."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def match_to_model(self, model=None) -> pd.DataFrame:
        """Match each target location to its grid cell (defaults to the bound model)."""

        return self.targets.match_to_model(self.model if model is None else model)

    def simulated_heads(self, model=None) -> pd.DataFrame:
        """Simulated heads at the target locations (defaults to the bound model)."""

        return self.targets.simulated_heads(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated heads (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the head targets (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ):
        """A FloPy ``continuous`` observation dict for the head targets on the bound model."""

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
        """Attach an MF6 head observation package to the bound model."""

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the head targets on the bound model."""

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
        """Wrap a :class:`BoundHeadTargets` to expose its plotting helpers."""

        self.bound_targets = bound_targets

    @property
    def model(self):
        """The model the underlying targets are bound to."""

        return self.bound_targets.model

    @property
    def targets(self):
        """The underlying :class:`HeadTargets` being plotted."""

        return self.bound_targets.targets

    def locations(self, *, ax=None, **kwargs):
        """Plot target locations on a map axis."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def calibration(self, *, type: str = "calibration"):
        """Return a calibration/heads plot for the bound targets."""

        return self.bound_targets.calibration_plot(type=type)

    def obs_vs_sim(self, *, baseline=None, ax=None, backend: str = "plotly"):
        """Return a target-vs-simulated cross plot for one or two models.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(current, baseline_compare=baseline_frame, backend=backend)

    def timeseries(self, name: str | None = None, *, baseline=None, ax=None, backend: str = "plotly"):
        """Return a time-series plot for one target location.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
            backend=backend,
        )

    def residuals_by_period(self, *, baseline=None, backend: str = "plotly"):
        """Return a by-period residual summary plot.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(
            current,
            baseline_compare=baseline_frame,
            backend=backend,
        )


