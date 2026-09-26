"""SFR stage/flow observation targets (+ bound views)."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _aggregate_named_series,
    _clean_observation_output_values,
    _copy_frame,
    _default_compare_stats,
    _default_obs_csv_name,
    _find_named_observation_output,
    _normalize_compare_time,
    _normalize_named_integer_locations,
    _normalize_named_series_values,
    _numeric_periods,
    _observation_label,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots


def _exchange_column(frame: pd.DataFrame) -> str:
    """Return the SFR exchange column, tolerating the frame name (``q_gwf``).

    ``sfr.results.q.get()`` names the exchange column by its reference frame
    (``q_gwf``); older/mock frames may still use the bare ``q``. This flow
    observation takes ``.abs()`` of it, so the frame's sign is irrelevant here.
    """

    return "q_gwf" if "q_gwf" in frame.columns else "q"


@dataclass
class SfrStageTargets:
    """Observed SFR stage time series tied to named reaches, for review or PEST.

    The SFR-stage analogue of :class:`HeadTargets`: it pairs named stream reaches
    (``locations`` -- a ``name -> reach`` mapping, zero-based reach numbers) with
    observed stage ``values`` through time, then aligns them with the model's
    simulated SFR stage. Pass an instance straight to ``cal.observe(...)`` for
    calibration, or use ``.compare(model)`` / ``.stats(model)`` to review
    residuals without PEST. Omitting ``values`` yields an observation *template*.

    Parameters
    ----------
    locations
        ``name -> reach`` definitions: a dict, ``Series``, list of records, or a
        DataFrame with ``name``/``reach`` columns. Reach numbers are zero-based.
    values
        Observed stage table (long or wide). Optional (template mode).
    time_column, value_column
        Column names in ``values`` for the period/time index and the stage value.
    times
        Optional explicit period/time labels when ``values`` carries none.

    Examples
    --------
    >>> SfrStageTargets(locations={"gage_main": 12},
    ...                 values=observed_stage_df)        # -> cal.observe(...)
    >>> SfrStageTargets(locations={"gage_main": 12}).compare(model)

    See Also
    --------
    SfrFlowTargets, HeadTargets, LakeStageTargets
    """

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "stage"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Normalize the ``name -> reach`` locations and observed values into internal tables."""

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
        """Return the named reach definitions, or the observed values merged with them."""

        definitions = self._locations.copy()
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions, on="name", how="left"))

    def to_long(self) -> pd.DataFrame:
        """Return the long-format SFR-stage target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per reach series."""

        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = self.get().pivot_table(index="time", columns="name", values="stage_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary (reach count, row count) of the target set."""

        return pd.DataFrame([{"n_reaches": int(len(self._locations)), "n_rows": int(len(self._values))}])

    def _obs_column_map(self) -> dict[str, str]:
        """Map each target name to its MF6 observation label (for matching output CSVs)."""

        return {
            str(name): _observation_label(name)
            for name in self._locations["name"].astype(str).tolist()
        }

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-stage table: a ``time`` column and one column per target.

        Read from the model's SFR stage observations when one workspace CSV holds
        every target's column (MF6's observation labels, as :meth:`to_flopy_obs`
        writes them), else from ``packages.sfr.results.stage``. The two routes
        label ``time`` differently: the observation route carries MF6's simulated
        time, one row per time step; the stage-file route the zero-based stress
        period.

        **On the stage-file route a period holds ONE stage per reach: its LAST
        saved time step** (ledger 331), the period-end rule ``HeadTargets`` uses.
        The route used to ADD every time step a period saved, so under ``SAVE HEAD
        ALL`` with three steps per period a reach at 141.43 ft read 424.29 --
        a stage is not added (ledger 297). One saved step per period (``SAVE HEAD
        LAST``, or one time step) reads as before. A dry reach is still MF6's
        -1e30 on this route: SFR's stage is not masked yet (ledger 329).

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns. The observation route needs its workspace;
            the stage-file route needs an SFR ``stage_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per target name. On the stage-file
            route a target naming a reach the model does not have gets no column.
            Empty (those columns, no rows) when the stage file holds no record.

        Raises
        ------
        ValueError
            On the stage-file route, when the rows of one reach in one saved
            record are not copies of one stage -- naming the reach, the period
            and the spread. MF6's stage table cannot produce that; it guards a
            table that is not the one this reads.

        Examples
        --------
        >>> targets = mf.SfrStageTargets(locations={"gage": 12})
        >>> series = targets.simulated_series(model)   # one row per stress period
        """

        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            return observed
        frame = model.packages.sfr.results.stage.get().copy()
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = _aggregate_named_series(
            frame, id_column="reach", value_column="stage", kind="elevation"
        )
        merged = self._locations.merge(frame, on="reach", how="left")
        wide = merged.pivot_table(index="per", columns="name", values="stage", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target reach stages against simulated SFR stages (with residuals).

        ``sim_stage`` comes from the same two routes as :meth:`simulated_series`:
        on the stage-file route, one stage per reach per stress period, its last
        saved time step (ledger 331). Targets are matched on ``per``, read from
        their ``time`` as a zero-based integer -- compared with MF6's simulated
        time on the observation route and with the stress period on the
        stage-file route.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row when the targets hold no
            values): ``name``, ``reach``, ``time``, ``per``, ``stage_target``,
            ``sim_stage``, ``residual`` (simulated minus target) and
            ``abs_residual``.

        Raises
        ------
        ValueError
            On the stage-file route, when a reach's copies in one record
            disagree (as :meth:`simulated_series`), or when a target's ``time``
            is not a number.
        """

        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = self._locations.merge(simulated, on="name", how="left")
        else:
            simulated = _aggregate_named_series(
                model.packages.sfr.results.stage.get().copy(),
                id_column="reach",
                value_column="stage",
                kind="elevation",
            )
            simulated = self._locations.merge(simulated, on="reach", how="left").rename(columns={"stage": "sim_stage"})
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["stage_target"] = np.nan
        else:
            targets = self.get().copy()
            targets["per"] = _numeric_periods(targets, caller="SfrStageTargets.compare(...)")
            compare = targets.merge(simulated, on=["name", "reach", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_stage"] - compare["stage_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table (target/sim/residual) for SFR-stage targets."""

        frame = self.compare(model)
        return frame.loc[:, ["name", "reach", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the SFR-stage target set."""

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "sfr_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 SFR-stage observations (one-based reaches)."""

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
        """Attach an MF6 utility observation package to SFR for these stage targets."""

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
    """Observed SFR flow time series tied to named reaches, for review or PEST.

    The streamflow analogue of :class:`SfrStageTargets`: it pairs named stream
    reaches (``locations`` -- a ``name -> reach`` mapping, zero-based reach
    numbers, typically gage reaches) with observed flow ``values`` through time,
    then aligns them with simulated SFR flow. Use it for gaged-discharge
    calibration targets. Pass an instance to ``cal.observe(...)`` /
    ``cal.forecast(...)``, or call ``.compare(model)`` / ``.stats(model)`` for
    standalone residual review. Omitting ``values`` yields a capture template.

    Parameters
    ----------
    locations
        ``name -> reach`` definitions: a dict, ``Series``, list of records, or a
        DataFrame with ``name``/``reach`` columns. Reach numbers are zero-based.
    values
        Observed flow table (long or wide). Optional (template mode).
    time_column, value_column
        Column names in ``values`` for the period/time index and the flow value.
    times
        Optional explicit period/time labels when ``values`` carries none.

    Examples
    --------
    >>> SfrFlowTargets(locations={"outlet_gage": 27},
    ...                values=observed_q_df)             # -> cal.observe(...)

    See Also
    --------
    SfrStageTargets, DrnFlowTargets, HeadTargets
    """

    locations: dict[str, int] | pd.Series | pd.DataFrame | list | tuple | None = None
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    time_column: str = "time"
    value_column: str = "flow"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Normalize the ``name -> reach`` locations and observed flows into internal tables."""

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
        """Return the named reach definitions, or the observed flows merged with them."""

        definitions = self._locations.copy()
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions, on="name", how="left"))

    def to_long(self) -> pd.DataFrame:
        """Return the long-format SFR-flow target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per reach series."""

        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = self.get().pivot_table(index="time", columns="name", values="flow_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary (reach count, row count) of the target set."""

        return pd.DataFrame([{"n_reaches": int(len(self._locations)), "n_rows": int(len(self._values))}])

    def _obs_column_map(self) -> dict[str, str]:
        """Map each target name to its MF6 observation label (for matching output CSVs)."""

        return {
            str(name): _observation_label(name)
            for name in self._locations["name"].astype(str).tolist()
        }

    def _package_flow_table(self, model) -> pd.DataFrame:
        """Build a per-period, per-reach simulated-flow table from the SFR budget output.

        Uses the FLOW-JA-FACE in-channel *routing* through-flow (outbound rows, plus
        TO-MVR) and returns ``per``/``reach``/``sim_flow``. If routing flow is
        unavailable on a real model it returns an empty frame with a warning — it does
        **not** substitute ``sfr.results.q`` (the stream-aquifer *exchange*, a different
        physical quantity, ~1e2 leakage), because relabeling that as ``sim_flow`` would
        silently feed PEST the wrong quantity under the streamflow column name.
        """

        if not hasattr(model, "outputs"):
            # Mock-only accommodation: a real model is a ``SimulationBase`` and always
            # exposes ``.outputs`` (a property), so this branch never runs in
            # production. Test doubles without ``.outputs`` feed flow-shaped data
            # through ``results.q``; keep reading it so those doubles resolve.
            try:
                flow = model.packages.sfr.results.q.get().copy()
            except (AttributeError, TypeError):
                # A test double that exposes neither `.outputs` NOR a usable
                # `results.q` -- an empty frame is the documented answer for
                # "this model has no streamflow to observe".
                return pd.DataFrame(columns=["per", "reach", "sim_flow"])
            if flow.empty:
                return pd.DataFrame(columns=["per", "reach", "sim_flow"])
            flow["per"] = pd.to_numeric(flow["per"], errors="coerce").astype(int)
            flow["reach"] = pd.to_numeric(flow["reach"], errors="coerce").astype(int)
            flow["sim_flow"] = _clean_observation_output_values(
                flow[_exchange_column(flow)]
            ).abs()
            return flow.loc[:, ["per", "reach", "sim_flow"]]

        flow = model.outputs.sfr.bud.get("FLOW-JA-FACE")
        if not isinstance(flow, pd.DataFrame) or flow.empty:
            # In-channel routing flow is unavailable -- the SFR budget file holds
            # no FLOW-JA-FACE rows. Do NOT substitute ``sfr.results.q`` here: that
            # is the stream-aquifer exchange, and feeding it to PEST as
            # ``sim_flow`` would be a silent wrong-quantity calibration. Surface
            # the gap instead.
            #
            # The ``isinstance`` half is now unreachable and kept only as a belt:
            # until ledger 301 ``SFRBudget.get`` warned and handed back a raw
            # record list whenever the SFR budget was saved on a different cadence
            # than heads, which is the case this guard was written for. That
            # comparison is gone -- the records are labelled from the SFR budget
            # file's own clock -- so an empty frame is the only route left.
            warnings.warn(
                "SFR in-channel routing flow (FLOW-JA-FACE) is unavailable, so no "
                "simulated SFR flow can be produced. Returning no simulated values "
                "rather than substituting the stream-aquifer exchange "
                "(sfr.results.q), which is a different physical quantity.",
                UserWarning,
                stacklevel=2,
            )
            return pd.DataFrame(columns=["per", "reach", "sim_flow"])

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
        """Return one wide simulated-flow table indexed by stress period."""

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
        """Compare target reach flows against simulated SFR flows (with residuals)."""

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
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table (target/sim/residual) for SFR-flow targets."""

        frame = self.compare(model)
        return frame.loc[:, ["name", "reach", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the SFR-flow target set."""

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "sfr_flow_targets.csv",
        kind: str = "DOWNSTREAM-FLOW",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 SFR-flow observations (one-based reaches)."""

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
        """Attach an MF6 utility observation package to SFR for these flow targets."""

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


class BoundSfrStageTargets:
    """Model-bound SFR stage helper returned from ``model.targets.sfr_stage``."""

    def __init__(self, model, targets: SfrStageTargets):
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

    def to_flopy_obs(self, *, csv_name: str = "sfr_stage_targets.csv", kind: str = "STAGE"):
        """A FloPy ``continuous`` observation dict for the SFR-stage targets."""

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
        """Attach an MF6 SFR-stage observation package to the bound model."""

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
        """Plotting helpers (obs-vs-sim, time series) for these SFR-stage targets."""

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

    def to_flopy_obs(self, *, csv_name: str = "sfr_flow_targets.csv", kind: str = "DOWNSTREAM-FLOW"):
        """A FloPy ``continuous`` observation dict for the SFR-flow targets."""

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
        """Attach an MF6 SFR-flow observation package to the bound model."""

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
        """Plotting helpers (obs-vs-sim, time series) for these SFR-flow targets."""

        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="flow_target",
                simulated_column="sim_flow",
                title="Observed vs simulated SFR flow",
                yaxis_title="Flow",
            )
        return self._plot


