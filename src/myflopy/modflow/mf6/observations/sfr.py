"""SFR stage/flow observation targets (+ bound views)."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _aggregate_named_series,
    _check_inequality,
    _clean_observation_output_values,
    _compare_without_values,
    _copy_frame,
    _declared_route,
    _default_compare_stats,
    _default_obs_csv_name,
    _measured_periods,
    _normalize_compare_time,
    _normalize_named_integer_locations,
    _normalize_named_series_values,
    _observation_label,
    _ObservedQuantity,
    _owning_package,
    _read_declared_observations,
    _reading_weight_first,
    _readings,
    _refuse_flows_leaving_the_frame,
    _refuse_keys_off_the_periods,
    _refuse_observation_kind,
    _refuse_reserved_target_names,
    _refuse_shared_labels,
    _resolve_feature_ids,
    _score_bounds,
    _stage_at_period_end,
    _wide_by_period,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.package_budget import build_sfr_stage_result_table
from myflopy.modflow.mf6.period_end import period_ends

#: What a reach-stage target reads on the observation route: an SFR package's
#: STAGE observations (ledger 344).
_SFR_STAGE = _ObservedQuantity(
    owner="SfrStageTargets",
    package_type="sfr",
    obstype="stage",
    feature="reach",
    file_route="the SFR stage file",
)

#: What a streamflow target reads there: an SFR package's DOWNSTREAM-FLOW
#: observations (ledger 344).
_SFR_FLOW = _ObservedQuantity(
    owner="SfrFlowTargets",
    package_type="sfr",
    obstype="downstream-flow",
    feature="reach",
    file_route="the SFR budget file",
)


def _calibration_stage_table(model) -> pd.DataFrame:
    """The SFR stage table a reach TARGET reads: the stage file RAW, a dry reach -1e30.

    The stage VIEWS read the file masked -- a dry (-1e30) or inactive (+1e30) reach is
    a gap (ledger 353) -- and a target does not: the user kept SFR targets as they
    were, so this reads ``build_sfr_stage_result_table(model, masked=False)`` and a
    dry reach reaches ``compare()`` and the files PEST reads as the -1e30 MF6 wrote,
    exactly as before the views were masked. (A lake target reads its file raw too,
    and gives a dry lake its bottom -- ledger 351; no rule was decided for a reach.)

    Parameters
    ----------
    model : SimulationBase
        A run model with one SFR package declaring a ``stage_filerecord``.

    Returns
    -------
    pandas.DataFrame
        :func:`~myflopy.modflow.mf6.package_budget.build_sfr_stage_result_table`'s
        columns, ``stage`` as MF6 wrote it, model length units.

    Raises
    ------
    AttributeError
        As the table builder: no SFR package, or no stage file.
    SeveralPackagesError
        As the table builder: the model carries several SFR packages.
    """

    return build_sfr_stage_result_table(model, masked=False)


def _sfr_stage_at_period_end(model, frame: pd.DataFrame) -> pd.DataFrame:
    """The SFR stage file read at each period's END, by the rule both stage families share.

    :func:`~myflopy.modflow.mf6.observations._shared._stage_at_period_end` on
    ``model.fl.sfr.output.stage()`` -- ``per``, ``reach``, ``stage``, one row per
    reach per period.

    Parameters
    ----------
    model : SimulationBase
        The run model the table was read from.
    frame : pandas.DataFrame
        :func:`_calibration_stage_table`'s rows (the stage file, raw).

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

    return _stage_at_period_end(
        model, frame, reader=lambda: model.fl.sfr.output.stage(),
        file="the SFR stage file", id_column="reach",
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
    final step output control did not save is refused rather than read early; on
    the observation route, MF6's observation at the time step that ends it (ledger
    346).

    **The observation route reads what the SFR package DECLARES** (ledger 344):
    a target reads MF6's STAGE observation only when its label is declared in
    the SFR package's own OBS as a STAGE of the target's own reach, from the file
    that OBS names; declared as another quantity (a ``DOWNSTREAM-FLOW`` under
    the same name is another observation, and is left alone), for another reach
    or twice, it raises, and so does a set declared only in part. It used to
    read the first workspace CSV whose header held the labels: a flow column
    read as a stage (-93.87 ft), reach 0 for a target naming reach 3. MF6's
    no-data value (3e30) reads as ``NaN`` there, as a lake's does.

    **A reach with no water surface is NOT a gap here** (ledger 353, the user's
    call): the stage VIEWS read the stage file masked, a dry (-1e30) or inactive
    (+1e30) reach left blank, while a target reads it RAW, so on the stage-file
    route a dry reach is -1e30 in :meth:`compare` and in the files PEST reads, as
    it always was. MF6's STAGE observation reports the reach's streambed top for
    a dry reach and, measured, for an inactive one, so the two routes differ there.

    **In a PEST build, its measured observations are in the set's prefix group**
    (``sfr_stage``, or the prefix ``observe(..., prefix=)`` names; ledger 358): a
    reach target carries no group label, so to see two sets of reaches apart in
    PEST++'s phi by group, register them as two sets under two prefixes. Its
    unmeasured periods are in ``<prefix>_unmeasured``.

    Parameters
    ----------
    locations : dict, pandas.Series, pandas.DataFrame, list, tuple or None, optional
        The named reaches: a ``{name: reach}`` dict, a ``Series`` indexed by
        name, or a ``DataFrame`` / list of record dicts with ``name`` and
        ``reach`` columns (``ValueError`` naming whichever is missing).
        ``reach`` is MF6's reach number,
        ZERO-based (:meth:`to_flopy_obs` writes it one-based). A whole number (``1.0``, or text reading as one, ``"3"``)
        is the number; other text is the reach's BOUNDNAME, found (any case) in
        the SFR packagedata when the targets are read and refused unless it
        names exactly one reach. A blank, a fraction or a negative raises
        ``ValueError`` here, and a number the package does not have raises when
        the targets are read (ledger 503 -- ``int()`` truncated
        ``0.7`` to 0 and a missing reach read nothing, in silence; MF6's files
        number reaches from ONE, so one past the last is named so). Each name is
        passed through ``str``. A name may be given ONCE: one given twice --
        ``1`` and ``"1"`` are one name -- raises ``ValueError`` naming each row
        and its reach (it kept the FIRST row of a table and the LAST of a
        ``Series``, in silence; ledger 609). Several names may share
        a reach. A table may carry the ``weight_column``, each reach target's
        PEST weight. Default ``None``: no targets.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed stage through time, in the model's length units. A path is
        read with ``pandas.read_csv``. LONG form has ``name`` and
        ``value_column`` columns (headers matched ignoring case); WIDE form has
        one column per target name beside ``time_column``; a table of only
        ``time_column`` and ``value_column``, or a ``Series`` indexed by time,
        is one series and needs exactly one target reach (``ValueError``
        otherwise). Names are NOT checked against ``locations`` and are joined
        to them exactly (case counts): a row naming no target gets no ``reach``
        and no simulated stage in :meth:`compare`, and the PEST build refuses
        it. A blank value (missing, or empty text) is unmeasured; any other
        entry that is not a finite number raises ``ValueError`` naming the target and
        its column (ledger 340).
        A LONG or one-series table may carry a ``weight`` column (matched
        ignoring case): each READING's own PEST weight (ledger 504). It
        beats the location's weight (``weight_column``), which beats the set's
        (the observation spec's ``weight``), which beats 1.0; a blank means "none of its own", and
        anything but a blank or a finite number from 0 raises ``ValueError``.
        PESTPP-IES draws a reading's noise with standard deviation 1/weight. A
        WIDE table cannot carry one (its ``weight`` column would be a target of
        that name) and raises, naming the LONG form. Any other extra column of
        a LONG table is not read.
        Default ``None`` makes a template that only captures simulated series:
        observed, its rows register at weight 0; forecast, every row is a
        forecast.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` and the PEST build read it as a
        zero-based stress period, by number (``1.0`` is period 1; a fraction,
        a negative key, a date or a duration raises -- ledger 499), on
        both routes (the observation route since ledger 346 -- MF6's simulated
        time was its key); a key naming no period the route's table holds is
        refused on both (ledger 500).
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed stage, matched
        ignoring case. Default ``"stage"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.
    package : str or None, optional
        The NAME of the SFR package the reaches belong to (FloPy's
        ``package_name``, any case), needed when the model carries several SFR
        packages, whose reach numbers each start at 0: only its observations are
        read, and :meth:`attach_flopy_obs` attaches to it. Default ``None``: the
        model's one SFR package. With several and no ``package``, reading raises
        ``ValueError`` naming them; so does a ``package`` the model does not
        have, and a named package that declares none of the targets (the stage
        file cannot be chosen by package).
    inequality : str or None, optional
        Whether every value of the set is a one-sided BOUND (ledger
        564): ``"at_most"`` -- the simulated stage must not exceed
        it -- or ``"at_least"`` -- must not fall below it. Default ``None``: an
        ordinary target; anything else raises ``ValueError``. In a PEST build each
        measured observation is in PESTPP-IES's ``l_``/``g_`` group
        (``l_sfr_stage``), which scores 0 wherever the bound holds; :meth:`compare`
        scores it so too (``residual`` 0 there, ``raw_residual`` simulated minus
        target). ``cal.forecast(...)`` refuses a bound set.
    weight_column : str or None, optional
        Column of a ``DataFrame`` or list-of-records ``locations`` holding each
        reach target's PEST weight, matched exactly and carried as ``weight``.
        Default ``"weight"``; ``None``, or a table without the default column,
        gives no location weights, and a ``{name: reach}`` dict or a ``Series``
        carries none. Any OTHER name the table lacks, or given with such a dict
        or ``Series``, raises ``ValueError`` naming the columns the table has
        (ledger 609). A blank is "none of its own"; anything but a blank or a finite
        number from 0 raises ``ValueError`` naming the target. It is the weight
        of each of the target's readings that holds a value and has no weight of
        its own in ``values`` -- a reading's own weight beats it, and it beats the
        set's (``SfrStageObservationSpec``'s ``weight``), which beats 1.0
        (ledger 504) -- in :meth:`get`, :meth:`compare` and the PEST build; a
        period with no target value carries weight 0 there (ledger 340).
        PESTPP-IES draws a reading's noise with standard deviation 1/weight, so
        a stage known to 0.02 ft is weight 50. Until ledger 609 the column was dropped in silence
        and the build gave 1.0.

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
    package: str | None = None
    inequality: str | None = None
    weight_column: str | None = "weight"

    #: What the family reads, for resolving and checking its reaches (ledger 503).
    _QUANTITY: ClassVar[_ObservedQuantity] = _SFR_STAGE

    def __post_init__(self):
        """Normalize the ``name -> reach`` locations and observed values into internal tables."""

        self.inequality = _check_inequality(self.inequality, owner=type(self).__name__)
        # The ONE reading of lake or reach locations (ledger 609): a
        # name given twice is refused (a table kept its FIRST row and a Series
        # its LAST, in silence), and a table's weight column is each reach
        # target's weight (it was dropped).
        self._locations = _normalize_named_integer_locations(
            self.locations, id_column="reach", target_name=type(self).__name__,
            weight_column=self.weight_column,
        )
        _refuse_reserved_target_names(
            self._locations["name"], owner=type(self).__name__,
            columns=("stage_target", "sim_stage"),
        )
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="stage_target",
            times=self.times,
            owner=type(self).__name__,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations) != 1:
                raise ValueError("SfrStageTargets values without explicit names require exactly one target reach.")
            self._values["name"] = str(self._locations["name"].iloc[0])

    def get(self) -> pd.DataFrame:
        """Return the named reach definitions, or the observed values merged with them.

        Returns
        -------
        pandas.DataFrame
            With no values, ``name``, ``reach`` and ``weight`` (the location's,
            ``NaN`` for none); with values, each values row joined on ``name``
            with ONE ``weight`` -- the reading's own, else its location's (ledger
            504) -- the weight the PEST build gives the observation, ``NaN``
            meaning the set's weight, else 1.0.
        """

        definitions = self._locations.copy()
        if self._values.empty:
            return definitions
        return _copy_frame(
            _reading_weight_first(_readings(self._values).merge(definitions, on="name", how="left"))
        )

    def to_long(self) -> pd.DataFrame:
        """Return the long-format SFR-stage target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per reach series.

        One row per stress period a ``time`` names -- ``time`` is the zero-based
        period as ``int``, read as :meth:`compare` and the PEST build read it, so ``"1"``
        and ``1`` are one row (they were two, ledger 422) -- and one column per target
        name holding its observed stage, in model length units; a set given no values gives ``time`` and one
        empty column per target. Two rows HOLDING a value for one target in one
        stress period are REFUSED -- names compared stripped and lower-cased --
        by the check :meth:`compare` and the PEST build read (ledger 340): this
        used to keep the first in silence (``pivot_table(aggfunc="first")``), so
        the wide table and the build disagreed about the same set (ledger
        371). A blank row is unmeasured and never collides; a target whose
        every row is blank keeps its column, all blank, as :meth:`compare` keeps
        its rows (ledger 422).

        Returns
        -------
        pandas.DataFrame
            ``time`` (``int64``, the period) plus one float column per target
            name, sorted.

        Raises
        ------
        ValueError
            When a ``time`` label is not a zero-based stress period -- not a
            number, or one with a fraction (``1.0`` is period 1) -- or two rows
            hold a value for one target in one period, naming the target, the
            period and the count.
        """

        if self._values.empty:
            return _wide_by_period(
                self._values, value_column="stage_target", caller="SfrStageTargets.to_wide()",
                names=self._locations["name"],
            )
        _measured_periods(self._values, value_column="stage_target", caller="SfrStageTargets.to_wide()")
        # After the check a (target, time) cell holds at most ONE value, so
        # "first" only prefers it over blank rows.
        return _wide_by_period(self.get(), value_column="stage_target", caller="SfrStageTargets.to_wide()")

    def summary(self) -> pd.DataFrame:
        """Return a compact summary (reach count, row count) of the target set."""

        return pd.DataFrame([{"n_reaches": int(len(self._locations)), "n_rows": int(len(self._values))}])

    def _locations_on(self, model) -> pd.DataFrame:
        """``name``/``reach`` with every reach a checked ZERO-based number in ``model`` (ledger 503)."""

        reaches = _resolve_feature_ids(
            model, self._QUANTITY, dict(self._locations[["name", "reach"]].itertuples(index=False)),
            package=self.package,
        )
        return pd.DataFrame(
            [{"name": name, "reach": reach} for name, reach in reaches.items()],
            columns=["name", "reach"],
        )

    def _observed(self, model):
        """The observation route's wide table, or ``None`` when the set reads the stage file.

        A reach's no-data value (3e30) is ``NaN`` here, as a lake's is: a stage
        of 0.0 ft is never the reading of "no stage" (ledger 344).
        """

        chosen = _declared_route(
            model,
            _SFR_STAGE,
            {
                str(name): (str(name), _observation_label(name), int(reach))
                for name, reach in self._locations_on(model)[["name", "reach"]].itertuples(index=False)
            },
            package=self.package,
        )
        if chosen is None:
            return None
        return _read_declared_observations(model, chosen, owner=_SFR_STAGE.owner, no_value=np.nan)

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-stage table: a ``time`` column and one column per target.

        Read from the SFR package's own STAGE observations when EVERY target's
        label (its name as :meth:`to_flopy_obs` writes it) is declared there, for
        the target's own reach -- only the files that OBS names are read, each
        record from its own column, MF6's no-data value (3e30) as ``NaN`` --
        else, when NO target is declared, from the SFR stage file read RAW
        (ledger 344; see the class notes for what is refused). Both routes label
        ``time`` by the zero-based stress period, one row per period at its END:
        the observation route takes MF6's observation at the time step that ends
        each period (ledger 346 -- it handed on every step, keyed by MF6's
        simulated time, so a target keyed by period matched nothing there), and
        refuses a file that does not end every period the TDIS declares.

        **On the stage-file route a period holds ONE stage per reach: the stage at
        the period's END** (ledgers 331, 336) -- the record that ends the period on
        the stage file's own clock (MF6 writes it when it saves HEAD), the one rule
        head and flow targets read too. The route used to ADD every time step a
        period saved, so under ``SAVE HEAD ALL`` with three steps per period a
        reach at 141.43 ft read 424.29 -- a stage is not added (ledger 297). One
        saved step per period that ends it (``SAVE HEAD LAST``, or one time step)
        reads as before; a period whose final step was not saved (``FIRST``,
        ``STEPS`` short of it) is refused, because its end-of-period stage does not
        exist in the file. **A dry reach is MF6's -1e30 on this route, and an
        inactive one +1e30**: the stage VIEWS leave both a gap (ledger 353), and a
        target reads the file raw because the user kept SFR targets as they were.
        The observation route reads MF6's STAGE observation there instead, which
        is the reach's streambed top for a dry reach and -- measured -- for an
        inactive one too, so the two routes differ at such a reach.

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns. The observation route reads the SFR
            package's declared OBS and the files it names in the model's
            workspace; the stage-file route needs an SFR ``stage_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per target name. Empty (those
            columns, no rows) when the stage file holds no record. A target
            naming a reach the SFR package does not have raises (it got no
            column on the stage-file route, ledger 503).

        Raises
        ------
        FileNotFoundError
            On the observation route, when a file the OBS declares is not in the
            workspace.
        ValueError
            When a target's reach is not one the SFR package has, or its
            boundname names none or several (ledger 503); when a
            target's label is declared in the SFR's OBS as another
            quantity only, for another reach or twice as a STAGE; when some
            targets of the set are declared and others are not; when the model
            carries several SFR packages and ``package`` does not name one, or
            names one that declares none of the targets; or when a declared
            file's columns do not line up with its records (ledger 344), or it
            does not end every stress period (a run that stopped, a truncated
            copy -- ledger 346). On the
            stage-file route: when a period the stage file saved lacks
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

        observed = self._observed(model)
        if observed is not None:
            return observed
        frame = _calibration_stage_table(model)
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].tolist()])
        frame = _sfr_stage_at_period_end(model, frame)
        merged = self._locations_on(model).merge(frame, on="reach", how="left")
        wide = merged.pivot_table(index="per", columns="name", values="stage", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target reach stages against simulated SFR stages (with residuals).

        ``sim_stage`` comes from the same two routes as :meth:`simulated_series`:
        on the stage-file route, one stage per reach per stress period, the stage
        at the period's END (ledgers 331, 336). Targets are matched on ``per``, read from
        their ``time`` as a zero-based integer -- the stress period, on both routes
        (ledger 346). On the observation route a key naming no period of the model
        is refused, naming the period each simulated-time key ends: MF6's time
        (``10, 20, 30``) was that route's key until ledger 346, and such a set
        would otherwise read ``NaN``. The stage-file route refuses a key naming a
        period its file does not hold too (ledger 500): it read ``NaN``, and
        :meth:`stats` left the row out in silence.

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
            ``weight`` (as :meth:`get`'s; for a set given no values it is the
            location's, though the build gives every such row weight 0),
            ``sim_stage``, ``residual`` (simulated minus target) and
            ``abs_residual``. For a BOUND set (``inequality``) ``residual`` is PEST++'s, 0 where the
            bound holds, and ``raw_residual`` after it is simulated minus target
            (ledger 564). With no values ``stage_target`` and the residuals are
            ``NaN``, and ``time`` is the period, equal to ``per``, on both
            routes: the table had no ``per`` -- it was renamed ``time`` -- so
            :meth:`residuals` raised ``KeyError`` (found in ledger 342).

        Raises
        ------
        FileNotFoundError
            As :meth:`simulated_series`: a declared observation file is missing.
        ValueError
            As :meth:`simulated_series`: a reach the package does not have,
            a target's label declared as something
            else, a set declared in part, several SFR packages and no
            ``package``, or a file that does not line up or does not end every
            period; on either route, when a target's key is no stress period the
            route's table holds (ledgers 346, 500); on the stage-file
            route, when a period the stage file saved lacks its final step or a
            reach's copies in one record disagree; or when a target's ``time`` is
            not a whole number, or two rows hold a value for one target in one
            period (the period rule PEST's build reads too, ledger 340).
        """

        locations = self._locations_on(model)
        # Each target's location weight, joined by name so the simulated side
        # stays id-only (ledger 609).
        weights = self._locations.loc[:, ["name", "weight"]]
        observed = self._observed(model)
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = locations.merge(simulated, on="name", how="left")
        else:
            simulated = _sfr_stage_at_period_end(model, _calibration_stage_table(model))
            simulated = locations.merge(simulated, on="reach", how="left").rename(columns={"stage": "sim_stage"})
        if self._values.empty:
            compare = _compare_without_values(
                simulated.merge(weights, on="name", how="left"), target_column="stage_target"
            )
        else:
            # The values joined to the CHECKED reaches, a boundname resolved --
            # what ``get()`` shows as declared -- with ONE weight: the reading's
            # own, else its location's (ledger 504).
            targets = _copy_frame(_reading_weight_first(_readings(self._values).merge(
                locations.merge(weights, on="name"), on="name", how="left"
            )))
            targets["per"] = _measured_periods(
                targets, value_column="stage_target", caller="SfrStageTargets.compare(...)"
            )
            if observed is not None:
                _refuse_keys_off_the_periods(
                    targets, observed["time"], model=model, caller="SfrStageTargets.compare(...)"
                )
            else:
                # Every file route refuses a key past the run's periods too
                # (ledger 500): it read NaN, and stats() left it out.
                _refuse_keys_off_the_periods(
                    targets, simulated["per"].dropna().unique(), model=model,
                    caller="SfrStageTargets.compare(...)", route=_SFR_STAGE.file_route,
                )
            compare = targets.merge(simulated, on=["name", "reach", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_stage"] - compare["stage_target"]
        compare["abs_residual"] = compare["residual"].abs()
        # A bound scores 0 where it holds (ledger 564).
        return _copy_frame(_score_bounds(compare, self.inequality))

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
            ``per``, ``stage_target``, ``sim_stage``, ``residual``,
            ``abs_residual`` and ``weight`` (as :meth:`compare`'s) -- and ``raw_residual`` for a bound set (ledger
            564).
            A set given no ``values`` -- a capture template, a forecast set --
            gets one row per simulated row, its targets and residuals ``NaN``:
            it raised ``KeyError`` naming ``per`` on both routes until ledger
            342's fix for ``DrnFlowTargets`` reached this family.
        """

        frame = self.compare(model)
        columns = [
            "name", "reach", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual", "weight",
        ]
        return frame.loc[:, columns + (["raw_residual"] if "raw_residual" in frame.columns else [])]

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
            Default ``"sfr_stage_targets.csv"``. Any name: :meth:`simulated_series`
            reads the file the SFR package's OBS declares, whatever it is called
            (ledger 344).
        kind : str, optional
            The MF6 SFR observation type written on every record. Must be
            ``"STAGE"`` (any case), the default -- the only type the observation
            route reads back as a stage -- and any other raises ``ValueError``
            (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "STAGE", reach + 1), ...]}``, one record per
            target: its name with every character that is not a letter or digit
            made ``_`` and leading and trailing ``_`` stripped (``"obs"`` if
            nothing is left), and its reach ONE-based. The shape
            ``ModflowUtlobs(continuous=)`` and an SFR package's
            ``observations=`` take.

        Raises
        ------
        ValueError
            When ``kind`` is not ``"STAGE"``, or two target names would be
            written as one label once MF6 upper-cases it (``"Gage A"`` and
            ``"gage-a"``).
        """

        kind = _refuse_observation_kind(kind, _SFR_STAGE)
        _refuse_shared_labels(self._locations["name"].astype(str), _SFR_STAGE)
        # A number is written from one, as MF6 numbers reaches; a boundname as it
        # is, which MF6 resolves in the SFR's packagedata.
        records = [
            (_observation_label(name), kind, reach if isinstance(reach, str) else int(reach) + 1)
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
            The model whose FloPy SFR package receives the OBS package unless
            ``package`` is given: the one the targets name (their ``package``
            field, ``ValueError`` when the model has no SFR package of that
            name), else the model's only SFR package, found by type (several
            raise ``ValueError`` naming them and ``package=``, none
            ``AttributeError``; ledger 352). The records themselves need no model.
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
            The MF6 SFR observation type, as in :meth:`to_flopy_obs`: ``"STAGE"``,
            the default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            the one ``model`` resolves (see ``model``). Default ``None``. Reading
            the observations back needs the targets' own ``package`` field to
            name it when the model carries several SFR packages.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.

        Raises
        ------
        ValueError
            As :meth:`to_flopy_obs`; or when the targets name a package the model
            does not have.
        """

        import flopy

        if package is not None:
            owner = package
        else:
            owner = _owning_package(model, _SFR_STAGE, self.package)
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
    refused rather than read early; on the observation route, MF6's observation at
    the time step that ends it (ledger 346). A reach's flow is what it sends downstream,
    to the next reach, to a mover or out of the network, as MF6's own
    ``DOWNSTREAM-FLOW`` observation measures it.

    **A flow is in MF6's own sign, on both routes** (ledger 350, the user's
    rule): MF6 writes water leaving the reach as NEGATIVE, so ``sim_flow`` is zero
    or negative, and a gauged streamflow ``Q`` is entered as ``-Q``. A target
    value above zero can never be matched and is REFUSED when the targets are
    built, naming ``-Q``; zero (a dry reach) is accepted. The budget-file route
    reported the magnitude until then while the observation route read MF6's
    negative number, so one reach read with opposite signs under two names, and
    a positive ``Q`` matched exactly scored a residual of ``-2Q`` -- a misfit only
    a smaller simulated flow lowers, which pulls a calibration toward drying the
    reach.

    **The observation route reads what the SFR package DECLARES** (ledger 344):
    a target reads MF6's ``DOWNSTREAM-FLOW`` observation only when its label is
    declared in the SFR package's own OBS as a ``DOWNSTREAM-FLOW`` of the
    target's own reach, from the file that OBS names; declared as another
    quantity only (``SFR``, the stream-aquifer exchange), for another reach or
    twice, it raises, and so does a set declared only in part. It used to read
    the first workspace CSV whose header held the labels.

    **In a PEST build, its measured observations are in the set's prefix group**
    (``sfr_flow``, or the prefix ``observe(..., prefix=)`` names; ledger 358): a
    reach target carries no group label, so to see two sets of reaches apart in
    PEST++'s phi by group, register them as two sets under two prefixes. Its
    unmeasured periods are in ``<prefix>_unmeasured``.

    Parameters
    ----------
    locations : dict, pandas.Series, pandas.DataFrame, list, tuple or None, optional
        The named (gage) reaches: a ``{name: reach}`` dict, a ``Series`` indexed
        by name, or a ``DataFrame`` / list of record dicts with ``name`` and
        ``reach`` columns (``ValueError`` naming whichever is missing).
        ``reach`` is MF6's reach number,
        ZERO-based (:meth:`to_flopy_obs` writes it one-based). A whole number (``1.0``, or text reading as one, ``"3"``)
        is the number; other text is the reach's BOUNDNAME, found (any case) in
        the SFR packagedata when the targets are read and refused unless it
        names exactly one reach. A blank, a fraction or a negative raises
        ``ValueError`` here, and a number the package does not have raises when
        the targets are read (ledger 503 -- ``int()`` truncated
        ``0.7`` to 0 and a missing reach read nothing, in silence; MF6's files
        number reaches from ONE, so one past the last is named so). Each name is
        passed through ``str``. A name may be given ONCE: one given twice --
        ``1`` and ``"1"`` are one name -- raises ``ValueError`` naming each row
        and its reach (it kept the FIRST row of a table and the LAST of a
        ``Series``, in silence; ledger 609). Several names may share
        a reach. A table may carry the ``weight_column``, each reach target's
        PEST weight. Default ``None``: no targets.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed streamflow through time, in model length cubed per model time
        unit, in MF6's own sign: water leaving the reach is NEGATIVE, so enter a
        gauged streamflow ``Q`` as ``-Q``. :meth:`compare` subtracts it, as given,
        from the simulated flow, which is MF6's ``DOWNSTREAM-FLOW`` on both
        routes (ledger 350). A value above zero is refused (``ValueError``
        naming the target, the time and ``-Q``); zero is accepted. A path is
        read with ``pandas.read_csv``. LONG
        form has ``name`` and ``value_column`` columns (headers matched
        ignoring case); WIDE form has one column per target name beside
        ``time_column``; a table of only ``time_column`` and ``value_column``,
        or a ``Series`` indexed by time, is one series and needs exactly one
        target reach (``ValueError`` otherwise). Names are NOT checked against
        ``locations`` and are joined to them exactly (case counts): a row
        naming no target gets no ``reach`` and no simulated flow, and the PEST
        build refuses it. A blank value (missing, or empty text) is unmeasured; any other
        entry that is not a finite number raises ``ValueError`` naming the target and
        its column (ledger 340).
        A LONG or one-series table may carry a ``weight`` column (matched
        ignoring case): each READING's own PEST weight (ledger 504). It
        beats the location's weight (``weight_column``), which beats the set's
        (the observation spec's ``weight``), which beats 1.0; a blank means "none of its own", and
        anything but a blank or a finite number from 0 raises ``ValueError``.
        PESTPP-IES draws a reading's noise with standard deviation 1/weight. A
        WIDE table cannot carry one (its ``weight`` column would be a target of
        that name) and raises, naming the LONG form. Any other extra column of
        a LONG table is not read. Default ``None`` makes a
        capture template: observed, its rows register at weight 0; forecast,
        every row is a forecast.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` and the PEST build read it as a
        zero-based stress period, by number (``1.0`` is period 1; a fraction,
        a negative key, a date or a duration raises -- ledger 499), on
        both routes (the observation route since ledger 346 -- MF6's simulated
        time was its key); a key naming no period the route's table holds is
        refused on both (ledger 500).
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed flow, matched
        ignoring case. Default ``"flow"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.
    package : str or None, optional
        The NAME of the SFR package the reaches belong to (FloPy's
        ``package_name``, any case), needed when the model carries several SFR
        packages, whose reach numbers each start at 0: only its observations are
        read, and :meth:`attach_flopy_obs` attaches to it. Default ``None``: the
        model's one SFR package. With several and no ``package``, reading raises
        ``ValueError`` naming them; so does a ``package`` the model does not
        have, and a named package that declares none of the targets (the budget
        file cannot be chosen by package).
    inequality : str or None, optional
        Whether every value of the set is a one-sided BOUND (ledger
        564), on the number AS WRITTEN -- MF6's sign, a reach's
        streamflow zero or negative (ledger 350). ``"at_most"``: the simulated
        flow must not exceed the value; ``"at_least"``: must not fall below it.
        So "at least 2 cfs leaves the reach" is ``-2.0`` (in the model's units)
        with ``inequality="at_most"`` -- simulated at most -2, i.e. at least 2
        leaving -- and "at most 5 leaves" is ``-5.0`` with ``"at_least"``.
        Default ``None``: an ordinary target; anything else raises ``ValueError``.
        A value above zero is refused either way. In a PEST build each measured
        observation is in PESTPP-IES's ``l_``/``g_`` group (``l_sfr_flow``), which
        scores 0 wherever the bound holds; :meth:`compare` scores it so too.
        ``cal.forecast(...)`` refuses a bound set.
    weight_column : str or None, optional
        Column of a ``DataFrame`` or list-of-records ``locations`` holding each
        reach target's PEST weight, matched exactly and carried as ``weight``.
        Default ``"weight"``; ``None``, or a table without the default column,
        gives no location weights, and a ``{name: reach}`` dict or a ``Series``
        carries none. Any OTHER name the table lacks, or given with such a dict
        or ``Series``, raises ``ValueError`` naming the columns the table has
        (ledger 609). A blank is "none of its own"; anything but a blank or a finite
        number from 0 raises ``ValueError`` naming the target. It is the weight
        of each of the target's readings that holds a value and has no weight of
        its own in ``values`` -- a reading's own weight beats it, and it beats the
        set's (``SfrFlowObservationSpec``'s ``weight``), which beats 1.0
        (ledger 504) -- in :meth:`get`, :meth:`compare` and the PEST build; a
        period with no target value carries weight 0 there (ledger 340).
        PESTPP-IES draws a reading's noise with standard deviation 1/weight, so
        a flow known to 10 ft3/d is weight 0.1, in 1/(model length3/time).
        Until ledger 609 the column was dropped in silence and the
        build gave 1.0.

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
    package: str | None = None
    inequality: str | None = None
    weight_column: str | None = "weight"

    #: What the family reads, for resolving and checking its reaches (ledger 503).
    _QUANTITY: ClassVar[_ObservedQuantity] = _SFR_FLOW

    def __post_init__(self):
        """Normalize the ``name -> reach`` locations and observed flows into internal tables."""

        self.inequality = _check_inequality(self.inequality, owner=type(self).__name__)
        # The ONE reading of lake or reach locations (ledger 609): a
        # name given twice is refused (a table kept its FIRST row and a Series
        # its LAST, in silence), and a table's weight column is each reach
        # target's weight (it was dropped).
        self._locations = _normalize_named_integer_locations(
            self.locations, id_column="reach", target_name=type(self).__name__,
            weight_column=self.weight_column,
        )
        _refuse_reserved_target_names(
            self._locations["name"], owner=type(self).__name__,
            columns=("flow_target", "sim_flow"),
        )
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="flow_target",
            times=self.times,
            owner=type(self).__name__,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations) != 1:
                raise ValueError("SfrFlowTargets values without explicit names require exactly one target reach.")
            self._values["name"] = str(self._locations["name"].iloc[0])
        _refuse_flows_leaving_the_frame(
            self._values,
            "flow_target",
            owner=type(self).__name__,
            quantity="MF6's DOWNSTREAM-FLOW",
            leaving="the reach",
            measured="a gauged streamflow",
        )

    def get(self) -> pd.DataFrame:
        """Return the named reach definitions, or the observed flows merged with them.

        Returns
        -------
        pandas.DataFrame
            With no values, ``name``, ``reach`` and ``weight`` (the location's,
            ``NaN`` for none); with values, each values row joined on ``name``
            with ONE ``weight`` -- the reading's own, else its location's (ledger
            504) -- the weight the PEST build gives the observation, ``NaN``
            meaning the set's weight, else 1.0.
        """

        definitions = self._locations.copy()
        if self._values.empty:
            return definitions
        return _copy_frame(
            _reading_weight_first(_readings(self._values).merge(definitions, on="name", how="left"))
        )

    def to_long(self) -> pd.DataFrame:
        """Return the long-format SFR-flow target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per reach series.

        One row per stress period a ``time`` names -- ``time`` is the zero-based
        period as ``int``, read as :meth:`compare` and the PEST build read it, so ``"1"``
        and ``1`` are one row (they were two, ledger 422) -- and one column per target
        name holding its observed flow, in model flow units and MF6's sign (ledger 350); a set given no values gives ``time`` and one
        empty column per target. Two rows HOLDING a value for one target in one
        stress period are REFUSED -- names compared stripped and lower-cased --
        by the check :meth:`compare` and the PEST build read (ledger 340): this
        used to keep the first in silence (``pivot_table(aggfunc="first")``), so
        the wide table and the build disagreed about the same set (ledger
        371). A blank row is unmeasured and never collides; a target whose
        every row is blank keeps its column, all blank, as :meth:`compare` keeps
        its rows (ledger 422).

        Returns
        -------
        pandas.DataFrame
            ``time`` (``int64``, the period) plus one float column per target
            name, sorted.

        Raises
        ------
        ValueError
            When a ``time`` label is not a zero-based stress period -- not a
            number, or one with a fraction (``1.0`` is period 1) -- or two rows
            hold a value for one target in one period, naming the target, the
            period and the count.
        """

        if self._values.empty:
            return _wide_by_period(
                self._values, value_column="flow_target", caller="SfrFlowTargets.to_wide()",
                names=self._locations["name"],
            )
        _measured_periods(self._values, value_column="flow_target", caller="SfrFlowTargets.to_wide()")
        # After the check a (target, time) cell holds at most ONE value, so
        # "first" only prefers it over blank rows.
        return _wide_by_period(self.get(), value_column="flow_target", caller="SfrFlowTargets.to_wide()")

    def summary(self) -> pd.DataFrame:
        """Return a compact summary (reach count, row count) of the target set."""

        return pd.DataFrame([{"n_reaches": int(len(self._locations)), "n_rows": int(len(self._values))}])

    def _locations_on(self, model) -> pd.DataFrame:
        """``name``/``reach`` with every reach a checked ZERO-based number in ``model`` (ledger 503)."""

        reaches = _resolve_feature_ids(
            model, self._QUANTITY, dict(self._locations[["name", "reach"]].itertuples(index=False)),
            package=self.package,
        )
        return pd.DataFrame(
            [{"name": name, "reach": reach} for name, reach in reaches.items()],
            columns=["name", "reach"],
        )

    def _observed(self, model):
        """The observation route's wide table, or ``None`` when the set reads the budget file."""

        chosen = _declared_route(
            model,
            _SFR_FLOW,
            {
                str(name): (str(name), _observation_label(name), int(reach))
                for name, reach in self._locations_on(model)[["name", "reach"]].itertuples(index=False)
            },
            package=self.package,
        )
        if chosen is None:
            return None
        return _read_declared_observations(model, chosen, owner=_SFR_FLOW.owner)

    def _package_flow_table(self, model) -> pd.DataFrame:
        """Each reach's streamflow at each stress period's END, from the SFR budget file, in MF6's sign.

        A reach's streamflow is the water it sends downstream in one time step:
        its outbound FLOW-JA-FACE rows (to the next reaches), plus ``TO-MVR`` and
        ``EXT-OUTFLOW`` (out of the network) where the budget file holds them --
        the quantity MF6's own ``DOWNSTREAM-FLOW`` observation reports. Those rows
        ADD, within one time step (``kind="extensive"``), each taken as a
        magnitude so the sum does not depend on the sign the file gives a term.

        **The flow is reported in MF6's own sign** (ledger 350, the user's rule):
        ``sim_flow`` is MINUS that sum -- water leaving the reach is negative, as
        MF6 writes ``DOWNSTREAM-FLOW`` -- so it equals the observation route's
        number (measured equal to the last digit at every period end on a plain
        reach, an outlet, and a reach and an outlet each passing water to a
        mover; ``docs/reproducers/flow_target_sign.py``). This
        route reported the magnitude, POSITIVE, while the observation route
        beside it read MF6's negative number, so the same reach read with
        opposite signs under two names, and a gauged flow matched exactly on the
        observation route scored a residual of ``-2Q``.

        **A period's value is its END-of-period value** (ledger 336): the record
        read is the one that ends the period on the SFR budget file's own clock
        (``SAVE BUDGET``, not ``SAVE HEAD`` -- ledger 301), through the one rule
        every target reads (``_aggregate_named_series``). It used to ADD every
        step a period saved: the canonical model at ``nstp=3`` under ``SAVE
        BUDGET ALL`` read its outlet reach at 113,453.29 where each step held
        37,817.76. A period whose final step the file does not hold is refused.

        If routing flow is unavailable it returns an empty frame with a warning
        -- it does **not** substitute ``sfr.results.q`` (the stream-aquifer
        *exchange*, a different physical quantity), because relabeling that as
        ``sim_flow`` would silently feed PEST the wrong quantity under the
        streamflow column name. A model without ``outputs`` -- only a test double
        -- is not read at all: a branch for such doubles read ``sfr.results.q``,
        that exchange, as the streamflow (ledger 350).

        Parameters
        ----------
        model : SimulationBase
            A run model with an SFR ``budget_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``per``, ``reach`` (zero-based) and ``sim_flow`` -- MF6's sign, zero
            or negative -- one row per reach the budget file holds an outbound
            record for, per stress period.

        Raises
        ------
        ValueError
            When a period the SFR budget file saved lacks its final step --
            naming the period, the steps saved and ``SAVE BUDGET LAST``.
        """

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
        # leaving each reach: the negative rows in the SFR package budget.
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
        # Each piece as a magnitude, so the sum is the reach's DOWNSTREAM-FLOW
        # whatever sign the file gives a term, then that sum in MF6's own sign.
        rows["sim_flow"] = rows["q"].abs()
        rows["per"] = rows["kstpkper"].map(lambda value: int(value[1]))
        table = _aggregate_named_series(
            rows, id_column="reach", value_column="sim_flow", kind="extensive", ends=ends
        )
        # ``0.0 - x``, not ``-x``: a dry reach reads 0.0, as MF6 writes it, not -0.0.
        table["sim_flow"] = 0.0 - table["sim_flow"]
        return table

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-flow table: a ``time`` column and one column per target.

        Read from the SFR package's own ``DOWNSTREAM-FLOW`` observations when
        EVERY target's label (its name as :meth:`to_flopy_obs` writes it) is
        declared there, for the target's own reach -- only the files that OBS
        names are read, each record from its own column -- else, when NO target
        is declared, from the SFR budget file (ledger 344; see the class notes
        for what is refused). Both routes label ``time`` by the zero-based stress
        period, one row per period at its END -- the observation route takes
        MF6's observation at the time step that ends each period (ledger 346; it
        handed on every step, keyed by MF6's simulated time) and refuses a file
        that does not end every period. Both carry MF6's own sign --
        ``DOWNSTREAM-FLOW``, negative when the reach sends water on -- and read
        the same number for one reach (ledger 350; the budget-file route
        reported the magnitude, so the two read opposite signs).

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
            ``run.model()`` returns. The observation route reads the SFR
            package's declared OBS and the files it names in the model's
            workspace; the budget-file route needs an SFR ``budget_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per target name, each a reach's
            streamflow in MF6's sign: zero or NEGATIVE, water leaving the reach
            (ledger 350). On the budget-file route a target naming a reach the
            SFR budget file holds no outbound record for gets no column.

        Raises
        ------
        FileNotFoundError
            On the observation route, when a file the OBS declares is not in the
            workspace.
        ValueError
            When a target's reach is not one the SFR package has, or its
            boundname names none or several (ledger 503); when a
            target's label is declared in the SFR's OBS as another
            quantity only, for another reach or twice as a ``DOWNSTREAM-FLOW``;
            when some targets of the set are declared and others are not; when
            the model carries several SFR packages and ``package`` does not name
            one, or names one that declares none of the targets; or when a
            declared file's columns do not line up with its records (ledger
            344), or it does not end every stress period (a run that stopped, a
            truncated copy -- ledger 346). On the budget-file route, when a
            period the SFR budget file saved
            lacks its final step -- naming the period, the steps saved and the
            ``SAVE BUDGET LAST`` setting that fixes it.

        Examples
        --------
        >>> targets = mf.SfrFlowTargets(locations={"gage": 27})
        >>> series = targets.simulated_series(model)   # one row per stress period
        """

        observed = self._observed(model)
        if observed is not None:
            return observed

        frame = self._package_flow_table(model)
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].astype(str).tolist()])
        merged = self._locations_on(model).merge(frame, on="reach", how="left")
        value_column = "sim_flow" if "sim_flow" in merged.columns else "q"
        wide = merged.pivot_table(index="per", columns="name", values=value_column, aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target reach flows against simulated SFR flows (with residuals).

        ``sim_flow`` comes from the same two routes as :meth:`simulated_series`:
        on the budget-file route, one flow per reach per stress period, the flow
        at the period's END (ledger 336). Targets are matched on ``per``, read
        from their ``time`` as a zero-based integer -- the stress period, on both
        routes (ledger 346; the observation route was keyed by MF6's simulated
        time, where a target keyed by period matched no row). On the observation
        route a key naming no period of the model is refused, naming the period
        each simulated-time key ends; on the budget-file route a key naming a
        period the file does not hold is refused too (ledger 500 -- it read
        ``NaN``, and :meth:`stats` left the row out in silence).

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
            ``weight`` (as :meth:`get`'s; for a set given no values it is the
            location's, though the build gives every such row weight 0),
            ``sim_flow``, ``residual`` (simulated minus target) and
            ``abs_residual``. ``flow_target`` and ``sim_flow`` are both in MF6's
            sign -- zero or negative, water leaving the reach (ledger 350) -- so a
            residual above zero means the model sends LESS water downstream than
            was gauged. For a BOUND set (``inequality``) ``residual`` is PEST++'s, 0 where the
            bound holds, and ``raw_residual`` after it is simulated minus target
            (ledger 564). With no values ``flow_target`` and the residuals are
            ``NaN``, and ``time`` is the period, equal to ``per``, on both
            routes: the table had no ``per`` -- it was renamed ``time`` -- so
            :meth:`residuals` raised ``KeyError`` (found in ledger 342).

        Raises
        ------
        FileNotFoundError
            As :meth:`simulated_series`: a declared observation file is missing.
        ValueError
            As :meth:`simulated_series`: a reach the package does not have,
            a target's label declared as something
            else, a set declared in part, several SFR packages and no
            ``package``, or a file that does not line up or does not end every
            period; on either route, when a target's key is no stress period the
            route's table holds (ledgers 346, 500); on the budget-file
            route, when a period the SFR budget file saved lacks its final step;
            or when a target's ``time`` is not a whole number, or two rows hold a value
            for one target in one period (the period rule PEST's build reads too,
            ledger 340).
        """

        locations = self._locations_on(model)
        # Each target's location weight, joined by name so the simulated side
        # stays id-only (ledger 609).
        weights = self._locations.loc[:, ["name", "weight"]]
        observed = self._observed(model)
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_flow")
            simulated["per"] = simulated["time"]
            simulated = locations.merge(simulated, on="name", how="left")
        else:
            simulated = self._package_flow_table(model)
            simulated = locations.merge(simulated, on="reach", how="left")
        if self._values.empty:
            compare = _compare_without_values(
                simulated.merge(weights, on="name", how="left"), target_column="flow_target"
            )
        else:
            # The values joined to the CHECKED reaches, a boundname resolved --
            # what ``get()`` shows as declared -- with ONE weight: the reading's
            # own, else its location's (ledger 504).
            targets = _copy_frame(_reading_weight_first(_readings(self._values).merge(
                locations.merge(weights, on="name"), on="name", how="left"
            )))
            targets["per"] = _measured_periods(
                targets, value_column="flow_target", caller="SfrFlowTargets.compare(...)"
            )
            if observed is not None:
                _refuse_keys_off_the_periods(
                    targets, observed["time"], model=model, caller="SfrFlowTargets.compare(...)"
                )
            else:
                # Every file route refuses a key past the run's periods too
                # (ledger 500): it read NaN, and stats() left it out.
                _refuse_keys_off_the_periods(
                    targets, simulated["per"].dropna().unique(), model=model,
                    caller="SfrFlowTargets.compare(...)", route=_SFR_FLOW.file_route,
                )
            compare = targets.merge(simulated, on=["name", "reach", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        # A bound scores 0 where it holds (ledger 564).
        return _copy_frame(_score_bounds(compare, self.inequality))

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
            ``per``, ``flow_target``, ``sim_flow``, ``residual``,
            ``abs_residual`` and ``weight`` (as :meth:`compare`'s) -- and ``raw_residual`` for a bound set (ledger
            564).
            A set given no ``values`` -- a capture template, a forecast set --
            gets one row per simulated row, its targets and residuals ``NaN``:
            it raised ``KeyError`` naming ``per`` on both routes until ledger
            342's fix for ``DrnFlowTargets`` reached this family.
        """

        frame = self.compare(model)
        columns = [
            "name", "reach", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual", "weight",
        ]
        return frame.loc[:, columns + (["raw_residual"] if "raw_residual" in frame.columns else [])]

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
            Default ``"sfr_flow_targets.csv"``. Any name: :meth:`simulated_series`
            reads the file the SFR package's OBS declares, whatever it is called
            (ledger 344).
        kind : str, optional
            The MF6 SFR observation type written on every record. Must be
            ``"DOWNSTREAM-FLOW"`` (any case), the default -- which MF6 reports
            NEGATIVE, the sign both routes read (ledger 350) -- and any other
            raises ``ValueError``: ``kind="SFR"`` wrote the
            stream-aquifer EXCHANGE under the targets' names, which that route
            read as streamflow (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "DOWNSTREAM-FLOW", reach + 1), ...]}``, one
            record per target: its name with every character that is not a
            letter or digit made ``_`` and leading and trailing ``_`` stripped
            (``"obs"`` if nothing is left), and its reach ONE-based. The shape
            ``ModflowUtlobs(continuous=)`` and an SFR package's
            ``observations=`` take.

        Raises
        ------
        ValueError
            When ``kind`` is not ``"DOWNSTREAM-FLOW"``, or two target names would
            be written as one label once MF6 upper-cases it.
        """

        kind = _refuse_observation_kind(kind, _SFR_FLOW)
        _refuse_shared_labels(self._locations["name"].astype(str), _SFR_FLOW)
        # A number is written from one, as MF6 numbers reaches; a boundname as it
        # is, which MF6 resolves in the SFR's packagedata.
        records = [
            (_observation_label(name), kind, reach if isinstance(reach, str) else int(reach) + 1)
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
            The model whose FloPy SFR package receives the OBS package unless
            ``package`` is given: the one the targets name (their ``package``
            field, ``ValueError`` when the model has no SFR package of that
            name), else the model's only SFR package, found by type (several
            raise ``ValueError`` naming them and ``package=``, none
            ``AttributeError``; ledger 352). The records themselves need no model.
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
            The MF6 SFR observation type, as in :meth:`to_flopy_obs`:
            ``"DOWNSTREAM-FLOW"``, the default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            the one ``model`` resolves (see ``model``). Default ``None``. Reading
            the observations back needs the targets' own ``package`` field to
            name it when the model carries several SFR packages.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.

        Raises
        ------
        ValueError
            As :meth:`to_flopy_obs`; or when the targets name a package the model
            does not have.
        """

        import flopy

        if package is not None:
            owner = package
        else:
            owner = _owning_package(model, _SFR_FLOW, self.package)
        csv_name = _default_obs_csv_name(filename, "sfr_flow_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )


class BoundSfrStageTargets:
    """Model-bound SFR stage helper returned from ``model.targets.sfr_stage``.

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from.
    """

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
        """The wide simulated series (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model with an SFR package to read. Default ``None`` uses the
            bound model; another run reads that run's files for the same
            targets. Read by :meth:`SfrStageTargets.simulated_series`'s two
            routes -- the SFR's declared STAGE observations when every target is
            declared there, else the SFR stage file read RAW (a dry reach -1e30,
            an inactive one +1e30) -- one row per stress period at its END, whose
            ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            ``time`` (the zero-based stress period) plus one stage column per
            target, in model length units.
        """

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated values (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model with an SFR package to compare against. Default
            ``None`` uses the bound model; another run compares the same targets
            against it -- what the plots' ``baseline=`` does. Read as
            :meth:`SfrStageTargets.compare` reads it (each target's ``time`` a
            zero-based stress period), whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row, ``stage_target`` and the
            residuals ``NaN``, when the set holds no values): ``name``,
            ``reach``, ``time``, ``per``, ``stage_target``, ``weight`` (as
            :meth:`SfrStageTargets.compare`'s), ``sim_stage``,
            ``residual`` (simulated minus target) and ``abs_residual``.
        """

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``reach``, ``time``,
            ``per``, ``stage_target``, ``sim_stage``, ``residual``,
            ``abs_residual`` and ``weight`` (as :meth:`compare`'s) -- and ``raw_residual`` for a bound set (ledger
            564).
        """

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            One row: ``n``, ``mean_error``, ``mae`` and ``rmse`` of simulated
            minus target stage, in model length units
            (:meth:`SfrStageTargets.stats`).
        """

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "sfr_stage_targets.csv", kind: str = "STAGE"):
        """A FloPy ``continuous`` observation dict for the SFR-stage targets.

        The records need no model -- each names its reach by number -- so the
        bound model takes no part.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"sfr_stage_targets.csv"``. Any name: the observation route
            reads whichever file the SFR package's OBS declares (ledger 344).
        kind : str, optional
            The MF6 SFR observation type written on every record. Must be
            ``"STAGE"`` (any case), the default; any other raises ``ValueError``
            (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "STAGE", reach + 1), ...]}``, the reach
            ONE-based, as :meth:`SfrStageTargets.to_flopy_obs` returns it (whose
            ``Raises`` apply).
        """

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
        """Attach an MF6 SFR-stage observation package to the bound model.

        Added to the bound model's FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        pname : str, optional
            FloPy package name of the OBS package. Default ``"sfr_stage_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory; its stem names the default ``csv_name``. Default
            ``"sfr_stage_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"sfr_stage_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 SFR observation type, as in :meth:`to_flopy_obs`:
            ``"STAGE"``, the default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to. Default
            ``None`` takes the bound model's SFR package: the one the targets
            name (their ``package`` field; ``ValueError`` when the model has no
            SFR package of that name), else its only SFR package, found by type
            (several raise ``ValueError`` naming them, none ``AttributeError``;
            ledger 352). Reading the observations back needs the targets' own
            ``package`` field to name it when the model carries several.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

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
    """Model-bound SFR flow helper returned from ``model.targets.sfr_flow``.

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from. Flows are in MF6's
    sign -- zero or negative, water the reach sends on (ledger 350).
    """

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
        """The wide simulated series (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model with an SFR package to read. Default ``None`` uses the
            bound model; another run reads that run's files for the same
            targets. Read by :meth:`SfrFlowTargets.simulated_series`'s two
            routes -- the SFR's declared ``DOWNSTREAM-FLOW`` observations when
            every target is declared there, else the SFR budget file -- one row
            per stress period at its END, whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            ``time`` (the zero-based stress period) plus one flow column per
            target, in model length cubed per time and MF6's sign (zero or
            negative).
        """

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated values (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model with an SFR package to compare against. Default
            ``None`` uses the bound model; another run compares the same targets
            against it -- what the plots' ``baseline=`` does. Read as
            :meth:`SfrFlowTargets.compare` reads it (each target's ``time`` a
            zero-based stress period), whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row, ``flow_target`` and the
            residuals ``NaN``, when the set holds no values): ``name``,
            ``reach``, ``time``, ``per``, ``flow_target``, ``weight`` (as
            :meth:`SfrFlowTargets.compare`'s), ``sim_flow``,
            ``residual`` (simulated minus target, above zero when the model
            sends LESS water downstream than was gauged) and ``abs_residual``.
        """

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``reach``, ``time``,
            ``per``, ``flow_target``, ``sim_flow``, ``residual``,
            ``abs_residual`` and ``weight`` (as :meth:`compare`'s) -- and ``raw_residual`` for a bound set (ledger
            564).
        """

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            One row: ``n``, ``mean_error``, ``mae`` and ``rmse`` of simulated
            minus target flow, in model length cubed per time
            (:meth:`SfrFlowTargets.stats`).
        """

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "sfr_flow_targets.csv", kind: str = "DOWNSTREAM-FLOW"):
        """A FloPy ``continuous`` observation dict for the SFR-flow targets.

        The records need no model -- each names its reach by number -- so the
        bound model takes no part.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"sfr_flow_targets.csv"``. Any name: the observation route
            reads whichever file the SFR package's OBS declares (ledger 344).
        kind : str, optional
            The MF6 SFR observation type written on every record. Must be
            ``"DOWNSTREAM-FLOW"`` (any case), the default, which MF6 reports
            negative; any other raises ``ValueError`` (``"SFR"`` would declare
            the stream-aquifer exchange under the targets' names, which the
            observation route refuses as another quantity when read back;
            ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "DOWNSTREAM-FLOW", reach + 1), ...]}``, the
            reach ONE-based, as :meth:`SfrFlowTargets.to_flopy_obs` returns it
            (whose ``Raises`` apply).
        """

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
        """Attach an MF6 SFR-flow observation package to the bound model.

        Added to the bound model's FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        pname : str, optional
            FloPy package name of the OBS package. Default ``"sfr_flow_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory; its stem names the default ``csv_name``. Default
            ``"sfr_flow_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"sfr_flow_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 SFR observation type, as in :meth:`to_flopy_obs`:
            ``"DOWNSTREAM-FLOW"``, the default, and nothing else
            (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to. Default
            ``None`` takes the bound model's SFR package: the one the targets
            name (their ``package`` field; ``ValueError`` when the model has no
            SFR package of that name), else its only SFR package, found by type
            (several raise ``ValueError`` naming them, none ``AttributeError``;
            ledger 352). Reading the observations back needs the targets' own
            ``package`` field to name it when the model carries several.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

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
