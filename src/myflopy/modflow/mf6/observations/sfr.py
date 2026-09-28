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
from myflopy.modflow.mf6.period_end import period_ends


def _exchange_column(frame: pd.DataFrame) -> str:
    """Return the SFR exchange column, tolerating the frame name (``q_gwf``).

    ``sfr.results.q.get()`` names the exchange column by its reference frame
    (``q_gwf``); older/mock frames may still use the bare ``q``. This flow
    observation takes ``.abs()`` of it, so the frame's sign is irrelevant here.
    """

    return "q_gwf" if "q_gwf" in frame.columns else "q"


def _stage_at_period_end(model, frame: pd.DataFrame) -> pd.DataFrame:
    """One stage per reach per stress period: the one that ENDS it (ledgers 331, 336).

    The SFR stage table holds one row per reach per saved step, labelled from the
    stage file's own clock -- which MF6 writes when it saves HEAD, not BUDGET
    (ledger 301). The period-end record is read off that file and checked against
    the TDIS, so a period whose end was not saved is refused rather than read at
    an earlier step.

    Parameters
    ----------
    model : SimulationBase
        The run model the table was read from.
    frame : pandas.DataFrame
        ``packages.sfr.results.stage.get()``: ``kstpkper``, ``per``, ``reach``,
        ``stage``.

    Returns
    -------
    pandas.DataFrame
        ``per``, ``reach``, ``stage``: one row per reach per period.

    Raises
    ------
    ValueError
        When a period the stage file saved lacks its final step, or a reach's
        rows in one record disagree.
    """

    ends = period_ends(
        model,
        reader=lambda: model.sfr.output.stage(),
        file="the SFR stage file",
        saves="HEAD",
        fallback=frame["kstpkper"] if "kstpkper" in frame.columns else (),
    )
    return _aggregate_named_series(
        frame, id_column="reach", value_column="stage", kind="elevation", ends=ends
    )


#: SFR budget terms that carry water OUT of a reach downstream -- to the next
#: reach, to a mover, or out of the network -- and so make up its streamflow.
_OUTFLOW_TERMS = ("TO-MVR", "EXT-OUTFLOW")


@dataclass
class SfrStageTargets:
    """Observed SFR stage time series tied to named reaches, for review or PEST.

    The SFR-stage analogue of :class:`HeadTargets`: it pairs named stream reaches
    (``locations`` -- a ``name -> reach`` mapping, zero-based reach numbers) with
    observed stage ``values`` through time, then aligns them with the model's
    simulated SFR stage. Pass an instance straight to ``cal.observe(...)`` for
    calibration, or use ``.compare(model)`` / ``.stats(model)`` to review
    residuals without PEST. Omitting ``values`` yields an observation *template*.

    A target keyed by stress period reads the stage at that period's END (ledger
    336): on the stage file, the record that ends the period, and a period whose
    final step output control did not save is refused rather than read early.

    Parameters
    ----------
    locations : dict, pandas.Series, pandas.DataFrame, list, tuple or None, optional
        The named reaches: a ``{name: reach}`` dict, a ``Series`` indexed by
        name, or a ``DataFrame`` / list of record dicts with ``name`` and
        ``reach`` columns (``ValueError`` naming whichever is missing; a
        repeated name keeps its FIRST row). ``reach`` is MF6's reach number,
        ZERO-based (:meth:`to_flopy_obs` writes it one-based); each is passed
        through ``int`` and each name through ``str``. Default ``None``: no
        targets.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed stage through time, in the model's length units. A path is
        read with ``pandas.read_csv``. LONG form has ``name`` and
        ``value_column`` columns (headers matched ignoring case); WIDE form has
        one column per target name beside ``time_column``; a table of only
        ``time_column`` and ``value_column``, or a ``Series`` indexed by time,
        is one series and needs exactly one target reach (``ValueError``
        otherwise). Names are NOT checked against ``locations`` and are joined
        to them exactly (case counts): a row naming no target gets no ``reach``
        and no simulated stage in :meth:`compare`.
        Default ``None`` makes a template that only captures simulated series.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` reads it as a zero-based stress
        period on the stage-file route.
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed stage, matched
        ignoring case. Default ``"stage"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.

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

        **On the stage-file route a period holds ONE stage per reach: the stage at
        the period's END** (ledgers 331, 336) -- the record that ends the period on
        the stage file's own clock (MF6 writes it when it saves HEAD), the one rule
        head and flow targets read too. The route used to ADD every time step a
        period saved, so under ``SAVE HEAD ALL`` with three steps per period a
        reach at 141.43 ft read 424.29 -- a stage is not added (ledger 297). One
        saved step per period that ends it (``SAVE HEAD LAST``, or one time step)
        reads as before; a period whose final step was not saved (``FIRST``,
        ``STEPS`` short of it) is refused, because its end-of-period stage does not
        exist in the file. A dry reach is still MF6's -1e30 on this route: SFR's
        stage is not masked yet (ledger 329).

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
            On the stage-file route: when a period the stage file saved lacks
            its final step -- naming the period, the steps saved and the ``SAVE
            HEAD LAST`` setting that fixes it; or when the rows of one reach in
            one saved record are not copies of one stage -- naming the reach, the
            period and the spread. MF6's stage table cannot produce the second;
            it guards a table that is not the one this reads.

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
        frame = _stage_at_period_end(model, frame)
        merged = self._locations.merge(frame, on="reach", how="left")
        wide = merged.pivot_table(index="per", columns="name", values="stage", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target reach stages against simulated SFR stages (with residuals).

        ``sim_stage`` comes from the same two routes as :meth:`simulated_series`:
        on the stage-file route, one stage per reach per stress period, the stage
        at the period's END (ledgers 331, 336). Targets are matched on ``per``, read from
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
            On the stage-file route, when a period the stage file saved lacks
            its final step or a reach's copies in one record disagree (as
            :meth:`simulated_series`), or when a target's ``time`` is not a
            number.
        """

        observed = _find_named_observation_output(model, expected_columns=self._obs_column_map())
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = self._locations.merge(simulated, on="name", how="left")
        else:
            simulated = _stage_at_period_end(model, model.packages.sfr.results.stage.get().copy())
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
        """Return a compact residual table (target/sim/residual) for SFR-stage targets.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``reach``, ``time``,
            ``per``, ``stage_target``, ``sim_stage``, ``residual`` and
            ``abs_residual``.
        """

        frame = self.compare(model)
        return frame.loc[:, ["name", "reach", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the SFR-stage target set.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual), ``mean_error``, ``mae`` and
            ``rmse`` of simulated minus target, in model length units. ``n`` is 0
            and the rest ``NaN`` when no row has a residual.
        """

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "sfr_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 SFR-stage observations (one-based reaches).

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"sfr_stage_targets.csv"``. :meth:`simulated_series` reads
            observations back from a workspace CSV holding every target's
            column, and skips any file whose name contains ``simulated_`` or
            ``_target_``.
        kind : str, optional
            The MF6 SFR observation type written on every record, upper-cased.
            Default ``"STAGE"``. The observation route of :meth:`simulated_series`
            and :meth:`compare` reads whatever that CSV holds as the stage.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, kind, reach + 1), ...]}``, one record per
            target: its name with every character that is not a letter or digit
            made ``_`` and leading and trailing ``_`` stripped (``"obs"`` if
            nothing is left), and its reach ONE-based. The shape
            ``ModflowUtlobs(continuous=)`` and an SFR package's
            ``observations=`` take.
        """

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
        """Attach an MF6 utility observation package to SFR for these stage targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase
            The model whose FloPy SFR package, ``model.gwf.sfr``, receives the
            OBS package unless ``package`` is given (a model with no SFR raises
            ``AttributeError`` there). The records themselves need no model.
        pname : str, optional
            FloPy package name. Default ``"sfr_stage_obs"``.
        filename : str, optional
            The OBS input file name FloPy writes, in the owning model's
            directory. Default ``"sfr_stage_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"sfr_stage_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 SFR observation type, as in :meth:`to_flopy_obs`. Default
            ``"STAGE"``.
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            ``model.gwf.sfr``. Default ``None``.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

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

    **A period's flow is its END-of-period flow** (ledger 336): on the SFR budget
    file, the flow in the time step that ends the period -- never the sum of the
    steps saved -- and a period whose final step output control did not save is
    refused rather than read early. A reach's flow is what it sends downstream,
    to the next reach, to a mover or out of the network, as MF6's own
    ``DOWNSTREAM-FLOW`` observation measures it.

    Parameters
    ----------
    locations : dict, pandas.Series, pandas.DataFrame, list, tuple or None, optional
        The named (gage) reaches: a ``{name: reach}`` dict, a ``Series`` indexed
        by name, or a ``DataFrame`` / list of record dicts with ``name`` and
        ``reach`` columns (``ValueError`` naming whichever is missing; a
        repeated name keeps its FIRST row). ``reach`` is MF6's reach number,
        ZERO-based (:meth:`to_flopy_obs` writes it one-based); each is passed
        through ``int`` and each name through ``str``. Default ``None``: no
        targets.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed streamflow through time, in model length cubed per model time
        unit. :meth:`compare` subtracts it from the simulated flow, which is a
        positive MAGNITUDE on the budget-file route and carries MF6's sign on
        the observation route. A path is read with ``pandas.read_csv``. LONG
        form has ``name`` and ``value_column`` columns (headers matched
        ignoring case); WIDE form has one column per target name beside
        ``time_column``; a table of only ``time_column`` and ``value_column``,
        or a ``Series`` indexed by time, is one series and needs exactly one
        target reach (``ValueError`` otherwise). Names are NOT checked against
        ``locations`` and are joined to them exactly (case counts): a row
        naming no target gets no ``reach`` and no simulated flow. Default
        ``None`` makes a capture template.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` reads it as a zero-based stress
        period on the budget-file route.
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed flow, matched
        ignoring case. Default ``"flow"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.

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
        """Each reach's streamflow at each stress period's END, from the SFR budget file.

        A reach's streamflow is the water it sends downstream in one time step:
        its outbound FLOW-JA-FACE rows (to the next reaches), plus ``TO-MVR`` and
        ``EXT-OUTFLOW`` (out of the network) where the budget file holds them --
        the quantity MF6's own ``DOWNSTREAM-FLOW`` observation reports (measured
        equal on the canonical model and on a stream with no mover), as a
        magnitude. Those rows ADD, within one time step (``kind="extensive"``).

        **A period's value is its END-of-period value** (ledger 336): the record
        read is the one that ends the period on the SFR budget file's own clock
        (``SAVE BUDGET``, not ``SAVE HEAD`` -- ledger 301), through the one rule
        every target reads (``_aggregate_named_series``). It used to ADD every
        step a period saved: the canonical model at ``nstp=3`` under ``SAVE
        BUDGET ALL`` read its outlet reach at 113,453.29 where each step held
        37,817.76. A period whose final step the file does not hold is refused.

        If routing flow is unavailable on a real model it returns an empty frame
        with a warning -- it does **not** substitute ``sfr.results.q`` (the
        stream-aquifer *exchange*, a different physical quantity, ~1e2 leakage),
        because relabeling that as ``sim_flow`` would silently feed PEST the wrong
        quantity under the streamflow column name.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR ``budget_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``per``, ``reach`` (zero-based) and ``sim_flow``, one row per reach
            that sends water downstream, per stress period.

        Raises
        ------
        ValueError
            When a period the SFR budget file saved lacks its final step --
            naming the period, the steps saved and ``SAVE BUDGET LAST``.
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

        budget = model.outputs.sfr.bud
        flow = budget.get("FLOW-JA-FACE")
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

        ends = period_ends(
            model,
            reader=lambda: budget.budget,
            file="the SFR budget file",
            saves="BUDGET",
            fallback=flow["kstpkper"],
        )
        frame = flow.copy()
        frame["q"] = _clean_observation_output_values(frame["q"])
        # For in-channel through-flow, use the outbound FLOW-JA-FACE rows
        # leaving each reach. Those are the negative rows in the SFR package
        # budget. Taking absolute values keeps the reported flow intuitive.
        pieces = [frame.loc[frame["q"] <= 0.0, ["kstpkper", "node", "q"]]]
        # Only the terms the file HOLDS: FloPy raises on a record text the budget
        # file does not have, and a stream with no mover writes no TO-MVR -- so
        # asking for it unconditionally failed every such model on this route.
        held = {str(term).strip().upper() for term in budget.types}
        for term in _OUTFLOW_TERMS:
            if term not in held:
                continue
            extra = budget.get(term)
            if isinstance(extra, pd.DataFrame) and not extra.empty:
                extra = extra.loc[:, ["kstpkper", "node", "q"]].copy()
                extra["q"] = _clean_observation_output_values(extra["q"])
                pieces.append(extra)
        rows = pd.concat(pieces, ignore_index=True)
        rows["reach"] = pd.to_numeric(rows["node"], errors="coerce")
        if rows["reach"].dropna().min() >= 1:
            rows["reach"] = rows["reach"] - 1
        rows["sim_flow"] = rows["q"].abs()
        rows["per"] = rows["kstpkper"].map(lambda value: int(value[1]))
        return _aggregate_named_series(
            rows, id_column="reach", value_column="sim_flow", kind="extensive", ends=ends
        )

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-flow table: a ``time`` column and one column per target.

        Read from the model's SFR flow observations when one workspace CSV holds
        every target's column (MF6's observation labels, as :meth:`to_flopy_obs`
        writes them), else from the SFR budget file. The two routes label
        ``time`` differently: the observation route carries MF6's simulated time,
        one row per time step, and MF6's own sign (``DOWNSTREAM-FLOW`` is
        negative); the budget-file route the zero-based stress period, and the
        flow's magnitude.

        **On the budget-file route a period holds ONE flow per reach: the flow
        at the period's END** (ledger 336) -- the record that ends the period on
        the SFR budget file's own clock, the same rule head and stage targets
        read. The route used to ADD every time step a period saved (3.0x at
        ``nstp=3`` under ``SAVE BUDGET ALL``). One saved step per period that ends
        it (``SAVE BUDGET LAST``, or one time step) reads as before; a period
        whose final step was not saved (``FIRST``, ``STEPS`` short of it) is
        refused, because its end-of-period flow does not exist in the file.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns. The observation route needs its workspace;
            the budget-file route an SFR ``budget_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per target name. On the budget-file
            route a target naming a reach that sends no water downstream gets no
            column.

        Raises
        ------
        ValueError
            On the budget-file route, when a period the SFR budget file saved
            lacks its final step -- naming the period, the steps saved and the
            ``SAVE BUDGET LAST`` setting that fixes it.

        Examples
        --------
        >>> targets = mf.SfrFlowTargets(locations={"gage": 27})
        >>> series = targets.simulated_series(model)   # one row per stress period
        """

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
        """Compare target reach flows against simulated SFR flows (with residuals).

        ``sim_flow`` comes from the same two routes as :meth:`simulated_series`:
        on the budget-file route, one flow per reach per stress period, the flow
        at the period's END (ledger 336). Targets are matched on ``per``, read
        from their ``time`` as a zero-based integer -- compared with the stress
        period on the budget-file route, and with MF6's simulated time on the
        observation route, where a target keyed by period matches no row.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row when the targets hold no
            values): ``name``, ``reach``, ``time``, ``per``, ``flow_target``,
            ``sim_flow``, ``residual`` (simulated minus target) and
            ``abs_residual``.

        Raises
        ------
        ValueError
            On the budget-file route, when a period the SFR budget file saved
            lacks its final step (as :meth:`simulated_series`), or when a
            target's ``time`` is not a number.
        """

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
        """Return a compact residual table (target/sim/residual) for SFR-flow targets.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``reach``, ``time``,
            ``per``, ``flow_target``, ``sim_flow``, ``residual`` and
            ``abs_residual``.
        """

        frame = self.compare(model)
        return frame.loc[:, ["name", "reach", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the SFR-flow target set.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual), ``mean_error``, ``mae`` and
            ``rmse`` of simulated minus target, in model flow units (length
            cubed per time). ``n`` is 0 and the rest ``NaN`` when no row has a
            residual.
        """

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "sfr_flow_targets.csv",
        kind: str = "DOWNSTREAM-FLOW",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 SFR-flow observations (one-based reaches).

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"sfr_flow_targets.csv"``. :meth:`simulated_series` reads
            observations back from a workspace CSV holding every target's
            column, and skips any file whose name contains ``simulated_`` or
            ``_target_``.
        kind : str, optional
            The MF6 SFR observation type written on every record, upper-cased.
            Default ``"DOWNSTREAM-FLOW"``, which MF6 reports NEGATIVE; the
            observation route of :meth:`simulated_series` and :meth:`compare`
            reads whatever that CSV holds, sign included, as the flow.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, kind, reach + 1), ...]}``, one record per
            target: its name with every character that is not a letter or digit
            made ``_`` and leading and trailing ``_`` stripped (``"obs"`` if
            nothing is left), and its reach ONE-based. The shape
            ``ModflowUtlobs(continuous=)`` and an SFR package's
            ``observations=`` take.
        """

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
        """Attach an MF6 utility observation package to SFR for these flow targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase
            The model whose FloPy SFR package, ``model.gwf.sfr``, receives the
            OBS package unless ``package`` is given (a model with no SFR raises
            ``AttributeError`` there). The records themselves need no model.
        pname : str, optional
            FloPy package name. Default ``"sfr_flow_obs"``.
        filename : str, optional
            The OBS input file name FloPy writes, in the owning model's
            directory. Default ``"sfr_flow_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"sfr_flow_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 SFR observation type, as in :meth:`to_flopy_obs`. Default
            ``"DOWNSTREAM-FLOW"``.
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            ``model.gwf.sfr``. Default ``None``.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

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


