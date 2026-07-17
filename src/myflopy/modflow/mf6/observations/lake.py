"""Lake stage observation targets (`LakeStageTargets` + bound view)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _aggregate_named_series,
    _copy_frame,
    _default_compare_stats,
    _default_obs_csv_name,
    _find_named_observation_output,
    _normalize_compare_time,
    _normalize_lake_stage_values,
    _numeric_periods,
    _observation_label,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots


@dataclass
class LakeStageTargets:
    """Observed lake-stage time series tied to named lakes, for review or PEST.

    The lake-stage analogue of :class:`HeadTargets`: it pairs named lakes
    (``locations`` -- a ``name -> lake`` mapping, zero-based MF6 lake ids) with
    observed stage ``values`` through time, then aligns them with the model's
    simulated LAK stage (read from the LAK stage output, or a workspace
    observation CSV when present). Pass an instance directly to
    ``cal.observe(...)`` for calibration, or call ``.compare(model)`` /
    ``.stats(model)`` to review residuals without a PEST workflow. Omitting
    ``values`` yields an observation *template* useful for capturing simulated
    series only.

    Parameters
    ----------
    locations
        ``name -> lake`` definitions: a dict, ``Series``, list of records, or a
        DataFrame with ``name``/``lake`` columns. Lake ids are zero-based.
    values
        Observed lake-stage table (long or wide). Optional (template mode).
    time_column, value_column
        Column names in ``values`` for the period/time index and the stage value.
    times
        Optional explicit period/time labels when ``values`` carries none.

    Examples
    --------
    >>> LakeStageTargets(locations={"valley_lake": 0},
    ...                  values=observed_stage_df)        # -> cal.observe(...)
    >>> LakeStageTargets(locations={"valley_lake": 0}).compare(model)

    See Also
    --------
    HeadTargets, SfrStageTargets, SfrFlowTargets, DrnFlowTargets
    """

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "stage"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Normalize the ``name -> lake`` locations and observed values into internal tables."""

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
        """Map each target name to its MF6 observation label (for matching output CSVs)."""

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
        compare = _normalize_compare_time(compare)
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

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        compare = self.compare(model)
        return CalibrationPlot.from_compare(
            compare,
            type="obs_vs_sim" if type == "calibration" else type,
            target_column="stage_target",
            simulated_column="sim_stage",
            title="Observed vs simulated lake stage" if type != "heads" else "Lake stage",
            yaxis_title="Simulated stage",
        )


class BoundLakeStageTargets:
    """Model-bound lake-stage helper returned from ``model.targets.lake_stage``."""

    def __init__(self, model, targets: LakeStageTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The target definitions (delegates to the wrapped targets)."""

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

    def simulated_series(self, model=None) -> pd.DataFrame:
        """The wide simulated series (defaults to the bound model)."""

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated values (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ):
        """A FloPy ``continuous`` observation dict for the lake-stage targets."""

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
        """Attach an MF6 lake-stage observation package to the bound model."""

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the lake-stage targets on the bound model."""

        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        """Plotting helpers (obs-vs-sim, time series) for these lake-stage targets."""

        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="stage_target",
                simulated_column="sim_stage",
                title="Observed vs simulated lake stage",
                yaxis_title="Stage",
            )
        return self._plot


