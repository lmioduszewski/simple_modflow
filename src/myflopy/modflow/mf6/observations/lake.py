"""Lake stage observation targets (`LakeStageTargets` + bound view)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.headsplus import _is_sentinel
from myflopy.modflow.mf6.observations._shared import (
    _compare_without_values,
    _copy_frame,
    _declared_route,
    _default_compare_stats,
    _default_obs_csv_name,
    _feature_id,
    _measured_periods,
    _normalize_compare_time,
    _normalize_lake_stage_values,
    _observation_label,
    _ObservedQuantity,
    _owning_package,
    _read_declared_observations,
    _refuse_keys_off_the_periods,
    _refuse_observation_kind,
    _refuse_reserved_target_names,
    _refuse_shared_labels,
    _resolve_feature_ids,
    _stage_at_period_end,
    _wide_by_period,
)
from myflopy.modflow.mf6.observations.hydrographs import (
    bin_readings,
    model_clock,
    read_readings,
    refuse_another_clock,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.package_budget import build_lak_stage_result_table
from myflopy.modflow.mf6.package_tables import lake_bottoms

#: What a lake target reads on the observation route: a LAK package's STAGE
#: observations (ledger 344).
_LAKE_STAGE = _ObservedQuantity(
    owner="LakeStageTargets",
    package_type="lak",
    obstype="stage",
    feature="lake",
    file_route="the LAK stage file",
)


def _calibration_stage_table(model) -> pd.DataFrame:
    """The LAK stage table a lake TARGET reads: a DRY lake at MF6's bottom, an INACTIVE one ``NaN``.

    The user's rule (ledger 351): a dry lake is CALIBRATED against its bottom, on
    both routes. MF6's STAGE observation already reports it there -- the solver
    clamps a lake's stage at its bottom, and the observation writes that stage --
    while the stage file writes -1e30, which the stage VIEWS keep as a gap (ledger
    329). So this reads the file RAW, where the two sentinels still differ: -1e30
    (dry) becomes the lake's bottom by MF6's own rule
    (:func:`~myflopy.modflow.mf6.package_tables.lake_bottoms`), and +1e30 (MF6's
    ``STATUS INACTIVE``) stays ``NaN``, which the PEST files write as 1e30.

    Parameters
    ----------
    model : SimulationBase
        A run model with one LAK package declaring a ``stage_filerecord``.

    Returns
    -------
    pandas.DataFrame
        :func:`~myflopy.modflow.mf6.package_budget.build_lak_stage_result_table`'s
        columns, ``stage`` in model length units: the file's stage, the lake's
        bottom where it wrote -1e30, ``NaN`` where it wrote +1e30.

    Raises
    ------
    AttributeError
        As the table builder: no LAK package, or no stage file. As
        ``lake_bottoms``, only when a lake is dry: several LAK packages, or no
        discretization package.
    ValueError
        As the table builder, when the file's records do not line up; as
        ``lake_bottoms``, when a DRY lake connects ``EMBEDDEDH``/``EMBEDDEDV``,
        whose bottom comes from a lake table this does not read.
    """

    frame = build_lak_stage_result_table(model, masked=False)
    stage = frame["stage"].to_numpy(dtype=float)
    no_stage = _is_sentinel(stage)
    dry = no_stage & (stage < 0.0)
    frame["stage"] = np.where(no_stage, np.nan, stage)
    if dry.any():
        lakes = frame.loc[dry, "lake"]
        bottoms = lake_bottoms(model, lakes=sorted(set(lakes.astype(int))))
        frame.loc[dry, "stage"] = lakes.astype(int).map(bottoms).to_numpy(dtype=float)
    return frame


def _lake_stage_at_period_end(model, frame: pd.DataFrame) -> pd.DataFrame:
    """The LAK stage file read at each period's END, by the rule both stage families share.

    :func:`~myflopy.modflow.mf6.observations._shared._stage_at_period_end` on
    ``model.fl.lak.output.stage()``, which MF6 writes when it saves HEAD. The LAK
    table was keyed by SAVED RECORD until ledger 348 (ledger 301's stage clock), so
    this only refused a schedule saving no period's end and kept every record: at
    ``nstp=3`` under ``SAVE HEAD ALL`` a target keyed by period 0, 1, 2 read steps
    1-3 of period 0.

    Parameters
    ----------
    model : SimulationBase
        A run model with a LAK ``stage_filerecord``.
    frame : pandas.DataFrame
        :func:`_calibration_stage_table`'s rows (the stage-file table with a dry
        lake at its bottom).

    Returns
    -------
    pandas.DataFrame
        ``per``, ``lake``, ``stage``: one row per lake per period, ``NaN`` where the
        lake is inactive.

    Raises
    ------
    ValueError
        When a period the LAK stage file saved lacks its final step -- naming the
        period, the steps saved and ``SAVE HEAD LAST`` -- or a lake's copies in one
        record disagree.
    """

    return _stage_at_period_end(
        model, frame, reader=lambda: model.fl.lak.output.stage(),
        file="the LAK stage file", id_column="lake",
    )


@dataclass
class LakeStageTargets:
    """Observed lake-stage time series tied to named lakes, for review or PEST.

    The lake-stage analogue of :class:`HeadTargets`: it pairs named lakes
    (``locations`` -- a ``name -> lake`` mapping, zero-based MF6 lake ids) with
    observed stage ``values`` through time, then aligns them with the model's
    simulated LAK stage (read from the LAK package's own STAGE observations
    when every target is declared there, else from the LAK stage file). Pass an
    instance directly to
    ``cal.observe(...)`` for calibration, or call ``.compare(model)`` /
    ``.stats(model)`` to review residuals without a PEST workflow. Omitting
    ``values`` yields an observation *template* useful for capturing simulated
    series only.

    A target keyed by stress period reads the stage at that period's END on both
    routes, as heads, reach stages and flows do (ledger 336): on the stage file,
    the record that ends the period on the file's own clock (ledger 348 -- that
    route was keyed by SAVED RECORD until then, so at several saves per period a
    period key read a step inside an EARLIER period), and a period whose final
    step output control did not save is refused rather than read early; on the
    observation route, MF6's observation at the time step that ends it (ledger
    346).

    **A DRY lake reads its BOTTOM, on both routes** (ledger 351, the user's
    rule): the stage MF6's own STAGE observation reports for it, since the
    solver never lets a lake's stage below its bottom. The stage file writes
    -1e30 there instead, so its route puts the lake's bottom in its place, by
    MF6's rule (:func:`~myflopy.modflow.mf6.package_tables.lake_bottoms`), and
    the two routes read one number record for record. An INACTIVE lake (MF6's
    ``STATUS INACTIVE``) has no stage on either route: ``NaN`` here, 1e30 in
    the files PEST reads. The stage VIEWS (``model.lak.results.stage`` and its
    pictures) still leave a dry day as a gap (ledger 329).

    **The observation route reads what the LAK package DECLARES** (ledger 344):
    a target reads MF6's STAGE observation only when its label (its name as
    :meth:`to_flopy_obs` writes it) is declared in the LAK package's own OBS as
    a STAGE of the target's own lake, and only the file that OBS names is read.
    A label declared as another quantity, for another lake or twice raises,
    naming the target and the record; so does a set of which some targets are
    declared and some are not. It used to read the first CSV anywhere in the
    workspace, in sorted order, whose header held the labels -- a GWF head, a
    LAK volume, another lake's stage, or a stale copy a previous run left there.

    **In a PEST build, its measured observations are in the set's prefix group**
    (``lak_stage``, or the prefix ``observe(..., prefix=)`` names; ledger 358 --
    ``stage`` until ledger 372): a lake
    target carries no group label, so to see two sets of lakes apart in PEST++'s phi
    by group, register them as two sets under two prefixes. Its unmeasured periods
    are in ``<prefix>_unmeasured``.

    Parameters
    ----------
    locations : dict, pandas.Series, pandas.DataFrame, list, tuple or None, optional
        The named lakes: a ``{name: lake}`` dict, a ``Series`` indexed by name,
        or a ``DataFrame`` / list of record dicts with ``name`` and ``lake``
        columns (``ValueError`` when either is missing). ``lake`` is MF6's lake
        number, ZERO-based (:meth:`to_flopy_obs` writes it one-based). A whole number (``1.0``, or text reading as one, ``"3"``)
        is the number; other text is the lake's BOUNDNAME, found (any case) in
        the LAK packagedata when the targets are read and refused unless it
        names exactly one lake. A blank, a fraction or a negative raises
        ``ValueError`` here, and a number the package does not have raises when
        the targets are read (ledger 503 -- ``int()`` truncated
        ``0.7`` to 0 and a missing lake read nothing, in silence; MF6's files
        number lakes from ONE, so one past the last is named so). Each name is
        passed through ``str``.
        Several names may share a lake. Default ``None``: no targets.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed lake stage through time, in the model's length units. A path is
        read with ``pandas.read_csv``. LONG form has ``name`` and
        ``value_column`` columns (headers matched ignoring case); WIDE form has
        one column per target name beside ``time_column``; a table of only
        ``time_column`` and ``value_column``, or a ``Series`` indexed by time,
        is one series and needs exactly one named lake (``ValueError``
        otherwise). Names are NOT checked against ``locations`` and are joined
        to them exactly (case counts): a row naming no target gets no ``lake``
        and no simulated stage in :meth:`compare`, and the PEST build refuses
        it. A blank value (missing, or empty text) is unmeasured; any other
        entry that is not a finite number raises ``ValueError`` naming the target and
        its column (ledger 340).
        A LONG or one-series table may carry a ``weight`` column (matched
        ignoring case): each READING's own PEST weight (ledger 504). It
        beats the set's weight (a ``LakeStageObservationSpec``'s ``weight``), which beats 1.0; a blank means "none of its own", and
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
        a negative key, a date or a duration raises -- ledger 499; dated
        readings are binned onto periods by :meth:`from_hydrographs`) --
        see :meth:`compare` for what it is matched against on each route.
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed stage, matched
        ignoring case. Default ``"stage"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.
    package : str or None, optional
        The NAME of the LAK package the lakes belong to -- FloPy's
        ``package_name``, as ``model.package_names`` lists it, any case --
        needed when the model carries several LAK packages, whose lake numbers
        each start at 0: only that package's STAGE observations are read, and
        :meth:`attach_flopy_obs` attaches to it. Default ``None``: the model's
        one LAK package. With several and no ``package``, reading raises
        ``ValueError`` naming them; so does a ``package`` the model does not
        have, and a named package that declares none of the targets (the stage
        file cannot be chosen by package).

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
    package: str | None = None

    #: What an INACTIVE lake -- the one lake with no stage for calibration -- is
    #: WRITTEN as in the CSV PEST reads (ledger 329), by both writers -- the base
    #: run and ``write_named_series_targets``, which reads it off the instance
    #: because pyEMU copies that function's source and not its module. The series
    #: itself says ``NaN``; PEST++ fails the run on an empty field and on ``nan``
    #: (``docs/reproducers/pest_no_value_fields.py``), so a lake inactive in a
    #: period nobody even weights would fail every realization. MF6's own no-value
    #: magnitude (the stage file's inactive lake, the head file's inactive cell). A
    #: DRY lake is NOT this: it reads its bottom (ledger 351), because 1e30 there
    #: stalled PESTPP-IES (lambda stuck near 1e57, phi 1.11 against 0.085 with the
    #: bottom; ``docs/reproducers/dry_lake_bottom.py``).
    _PEST_NO_STAGE: ClassVar[float] = 1.0e30
    #: The :class:`~myflopy.modflow.mf6.observations.hydrographs.Binning` of a
    #: :meth:`from_hydrographs` set; ``None`` otherwise.
    _binned: ClassVar[object] = None

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
            normalized[str(name)] = _feature_id(
                lake_id, name=name, owner=type(self).__name__, feature="lake"
            )
        _refuse_reserved_target_names(
            normalized, owner=type(self).__name__, columns=("stage", "sim_stage")
        )
        self._series = normalized
        self._values = _normalize_lake_stage_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            times=self.times,
            owner=type(self).__name__,
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
        """Create one named lake-stage target from a Series/list-like input.

        Parameters
        ----------
        lake : str
            The target's name.
        lake_id : int
            MF6's lake number for it, ZERO-based, or its boundname -- read as
            the constructor's ``locations`` reads a lake.
        values : pandas.Series or sequence of float
            Observed stages, in the model's length units. A ``Series`` carries
            its time labels as its index; any other sized sequence is one stage
            per entry of ``times``.
        times : sequence or None, optional
            Time labels for a non-``Series`` ``values``, one per value (another
            length raises ``ValueError`` from pandas). Default ``None``:
            ``0, 1, 2, ...``. Ignored for a ``Series``.
        time_column : str, optional
            Name of the time column in the table built here, stored on the
            result. Default ``"time"``.
        value_column : str, optional
            Name of the stage column in the table built from a ``Series``,
            stored on the result. Default ``"stage"``.

        Returns
        -------
        LakeStageTargets
            A target set holding the one lake.
        """

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
        """Create lake-stage targets from long-format records.

        Parameters
        ----------
        records : list of dict or pandas.DataFrame
            One row per target per time, as ``pandas.DataFrame(records)`` reads
            it, with columns ``name``, ``lake`` (MF6's lake number, ZERO-based),
            ``time_column`` and ``value_column``, all matched exactly.
            ``ValueError`` when there are no records or a column is missing. A
            name's lake is taken from its FIRST row.
        time_column : str, optional
            Column holding each record's time label. Default ``"time"``.
        value_column : str, optional
            Column holding each record's observed stage. Default ``"stage"``.

        Returns
        -------
        LakeStageTargets
            The targets, one per distinct ``name``.
        """

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

    @classmethod
    def from_hydrographs(
        cls,
        locations,
        readings,
        *,
        model,
        package: str | None = None,
        mode: str = "history",
        rule: str | None = None,
        date_column: str = "date",
        value_column: str = "stage",
        max_offset=None,
        between=None,
        min_years: int | None = None,
    ):
        """Lake-stage targets from DATED readings, binned onto ``model``'s stress periods.

        The lake twin of :meth:`HeadTargets.from_hydrographs
        <myflopy.modflow.mf6.observations.heads.HeadTargets.from_hydrographs>`,
        binning by the same rule (:mod:`~myflopy.modflow.mf6.observations.hydrographs`,
        ledger 543) against ``model``'s calendar. The result is an
        ordinary period-keyed set carrying :attr:`binning` and the CLOCK STAMP
        :attr:`calendar`, which :meth:`compare`, :meth:`simulated_series` and the
        PEST build check against the model they meet.

        Parameters
        ----------
        locations : dict, pandas.Series, pandas.DataFrame, list or tuple
            The named lakes, read as the constructor's ``locations``: ``{name:
            lake}`` with MF6's lake number ZERO-based, or the lake's boundname.
        readings : str, pathlib.Path, pandas.DataFrame or pandas.Series
            The dated stage readings, in the model's length units. A path is read
            with ``pandas.read_csv``. WIDE: ``date_column`` plus one column per
            lake target, headed by its name; LONG: ``date_column``, ``name`` and
            ``value_column``; ONE SERIES (or a ``Series`` indexed by date) for a
            single lake. Dates may be the index. A blank value is no reading;
            anything else that is not a finite number raises ``ValueError``. A
            name that is not a target name (matched EXACTLY, as lake values are)
            raises ``ValueError`` naming it; two readings of one lake at one
            instant collapse when equal and raise when not; a ``weight`` column
            raises (weight the set through ``cal.observe(..., weight=)``).
        model : SimulationBase
            The model the targets will be compared on, whose calendar the
            readings are binned against (its TDIS must declare
            ``START_DATE_TIME`` and ``TIME_UNITS``).
        package : str or None, default None
            The LAK package's name when the model carries several (the
            constructor's ``package``).
        mode : {"history", "average_year"}, default "history"
            As :meth:`HeadTargets.from_hydrographs
            <myflopy.modflow.mf6.observations.heads.HeadTargets.from_hydrographs>`:
            ``"history"`` bins the record onto the model's own dates (readings
            outside its span are left out and warned); ``"average_year"`` needs
            an ``mf.AverageYear`` model and puts each calendar month's target on
            its period of the last cycle.
        rule : str or None, default None
            History: ``"nearest_end"`` (default), ``"median"`` or ``"mean"``;
            average year: ``"median"`` (default; the median over years of each
            year's monthly median) or ``"mean"`` -- as for heads.
        date_column : str, default "date"
            The readings' date column, matched exactly; each entry a date with a
            four-digit year and no time zone. Text ``d/m/y`` or ``m/d/y`` is read
            MONTH first, and a column mixing day-first and month-first rows raises
            ``ValueError`` (ledger 545).
        value_column : str, default "stage"
            The stage column of a LONG or one-series table, matched ignoring case;
            the result's ``value_column``.
        max_offset : str, pandas.Timedelta, datetime.timedelta, numpy.timedelta64 or None, default None
            History with ``rule="nearest_end"`` only: a reading farther than this
            from its period's end is not used. Elsewhere it raises.
        between : tuple of two dates or None, default None
            Average year only: the readings used, ``(first, last)`` inclusive,
            either ``None`` for open; a date is its midnight. In history mode it
            raises.
        min_years : int or None, default None
            Average year only: a month read in fewer calendar years has no
            target; ``None`` means 1. In history mode it raises.

        Returns
        -------
        LakeStageTargets
            Keyed by zero-based stress period, carrying :attr:`binning` and
            :attr:`calendar`.

        Raises
        ------
        ValueError
            As listed under each parameter, and as the constructor raises.

        Warns
        -----
        UserWarning
            As :meth:`HeadTargets.from_hydrographs
            <myflopy.modflow.mf6.observations.heads.HeadTargets.from_hydrographs>`.
        """

        owner = cls.__name__
        known = cls(locations=locations, package=package)
        calendar, months = model_clock(model, mode=mode, owner=owner)
        long = read_readings(readings, date_column=date_column, value_column=value_column, owner=owner)
        if "name" not in long.columns:
            if len(known._series) != 1:
                raise ValueError(
                    f"{owner}.from_hydrographs(...): the readings name no lake (one series), "
                    f"and there are {len(known._series)} lake targets; give a WIDE table, one "
                    "column per target, or a LONG one with a 'name' column."
                )
            long["name"] = next(iter(known._series))
        unknown = sorted(set(long["name"].astype(str)) - set(known._series))
        if unknown:
            raise ValueError(
                f"{owner}.from_hydrographs(...): the readings name {len(unknown)} lake "
                "target(s) not among the locations (matched exactly, case counts): "
                + ", ".join(unknown[:10]) + (" ..." if len(unknown) > 10 else "") + "."
            )
        binning = bin_readings(
            long, calendar, mode=mode, rule=rule, max_offset=max_offset, between=between,
            min_years=min_years, months=months, owner=owner,
        )
        targets = cls(
            locations=locations,
            values=binning.values.rename(columns={"value": value_column}),
            time_column="time",
            value_column=value_column,
            package=package,
        )
        targets._binned = binning
        return targets

    @property
    def binning(self) -> pd.DataFrame | None:
        """How each target was made from dated readings, or ``None`` for a set keyed by period.

        Set by :meth:`from_hydrographs`; the columns are
        :attr:`HeadTargets.binning
        <myflopy.modflow.mf6.observations.heads.HeadTargets.binning>`'s.
        """

        binned = self._binned
        return None if binned is None else binned.report.copy()

    @property
    def calendar(self):
        """The :class:`~myflopy.modflow.mf6.period_calendar.PeriodCalendar` the targets were binned against, or ``None``.

        The CLOCK STAMP of a :meth:`from_hydrographs` set (ledger 543).
        """

        binned = self._binned
        return None if binned is None else binned.calendar

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
        """Return the target table in wide format with one column per lake series.

        One row per stress period a ``time`` names -- ``time`` is the zero-based
        period as ``int``, read as :meth:`compare` and the PEST build read it, so ``"1"``
        and ``1`` are one row (they were two, ledger 422) -- and one column per target
        name holding its observed stage, in model length units; a set given no
        values gives ``time`` and one empty column per target. Two rows HOLDING
        a value for one target in one stress period are REFUSED, by the check
        :meth:`compare` and the PEST build read (ledger 340): this used to keep
        the first in silence (``pivot_table(aggfunc="first")``), so the wide
        table and the build disagreed about the same set (ledger 371). A
        blank row is unmeasured and never collides; a target whose every row is
        blank keeps its column, all blank, as :meth:`compare` keeps its rows
        (ledger 422).

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
                self._values, value_column="stage", caller="LakeStageTargets.to_wide()",
                names=self._series.keys(),
            )
        _measured_periods(self._values, value_column="stage", caller="LakeStageTargets.to_wide()")
        # After the check a (target, period) cell holds at most ONE value, so
        # "first" only prefers it over blank rows.
        return _wide_by_period(self.get(), value_column="stage", caller="LakeStageTargets.to_wide()")

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the lake-stage target set."""

        return pd.DataFrame(
            [{"n_lakes": int(len(self._series)), "n_rows": int(len(self._values))}]
        )

    def _lakes_on(self, model) -> dict[str, int]:
        """Each target's ZERO-based lake number in ``model``'s LAK package, checked (ledger 503)."""

        return _resolve_feature_ids(model, _LAKE_STAGE, self._series, package=self.package)

    def _definitions_on(self, model) -> pd.DataFrame:
        """``name``/``lake`` with every lake a checked zero-based number in ``model``."""

        return pd.DataFrame(
            [{"name": name, "lake": lake} for name, lake in self._lakes_on(model).items()],
            columns=["name", "lake"],
        )

    def _observation_route(self, model):
        """The LAK STAGE observation each target reads, or ``None`` for the stage file (ledger 344)."""

        return _declared_route(
            model,
            _LAKE_STAGE,
            {name: (name, _observation_label(name), lake) for name, lake in self._lakes_on(model).items()},
            package=self.package,
        )

    def _observed(self, model):
        """The observation route's wide table, or ``None`` when the set reads the stage file."""

        chosen = self._observation_route(model)
        if chosen is None:
            return None
        return _read_declared_observations(model, chosen, owner=_LAKE_STAGE.owner, no_value=np.nan)

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-stage table: a ``time`` column and one column per target.

        Read from the LAK package's own STAGE observations when EVERY target's
        label (its name as :meth:`to_flopy_obs` writes it) is declared there, for
        the target's own lake -- only the files that OBS names are read, each
        record from its own column -- else, when NO target is declared, from
        the LAK stage file (ledger 344; see the class notes for what is
        refused). Both routes label ``time`` by the zero-based STRESS PERIOD, one
        row per period at its END: the observation route takes MF6's observation
        at the time step that ends each period (ledger 346; it carried every step,
        keyed by MF6's simulated time), the stage-file route the record that ends
        it on the stage file's own clock (ledger 348).

        **The stage-file route was keyed by SAVED RECORD until ledger 348**
        (ledgers 301, 331): its table labelled ``per`` by the record's index, so at
        ``SAVE HEAD ALL`` with ``nstp=3`` a target keyed by period 0, 1, 2 read
        steps 1-3 of period 0 while ``SfrStageTargets`` and ``HeadTargets`` beside
        it read each period's end, and a PEST build wrote nine rows for three
        periods. It reads the period-end record now, by the rule they read
        (``period_end.period_end_steps`` on the stage file's own clock); one saved
        step per period that ends it (``SAVE HEAD LAST``, or one time step) reads
        as before. A period whose final step was not saved (``SAVE HEAD FIRST``,
        ``STEPS`` short of it) is refused, because its end-of-period stage does not
        exist in the file.

        **A DRY lake reads its BOTTOM on both routes** (ledger 351, the user's
        rule), which is MF6's own number for it. The observation reports the stage
        MF6 solved, and the solver clamps a lake's stage at its bottom, so a dry
        lake's observation IS the bottom (100 ft on the reproducer's one-cell lake,
        measured); the stage file writes -1e30 there instead, and this route puts
        the lake's bottom in its place, computed by MF6's rule
        (:func:`~myflopy.modflow.mf6.package_tables.lake_bottoms`: a VERTICAL
        connection's cell top, a HORIZONTAL one's ``belev``, the lowest over the
        lake's connections). The two routes then read the same number at every
        record, dry or wet (measured equal to the bit on vertical, horizontal and
        excavated lakes loaded from their files, under LAK ``IMPLICIT`` too; a
        model still in memory from its build carries its elevations unrounded
        where FloPy wrote them at 8 decimals, so there a dry lake can differ from
        MF6's bottom by up to 5e-9).
        It was ``NaN`` on this route until then, and 1e30 in the file PEST read,
        which stalled PESTPP-IES (``docs/reproducers/dry_lake_bottom.py``). The
        stage VIEWS keep a dry day as a gap (ledger 329).

        **An INACTIVE lake is a row holding ``NaN``, never a missing row**
        (ledger 329), on both routes: the observation file's 3e30 had been read as
        a stage of 0.0 ft, and the stage file's +1e30 as itself.

        The stage-file route keeps one row per period and one column per
        target, even a target inactive in every period, so the shape does not
        depend on which lakes are active. That is what PEST needs: pyEMU builds the instruction
        file from the BASE run's copy of this table, one line advance per row and
        each value found by counting commas past a skipped ``time``, so a
        realization one row short fails at end of file and one with an extra
        column reads the wrong column without error (both measured,
        ``docs/reproducers/pest_no_value_fields.py``). The file PEST reads carries
        1e30 in place of ``NaN``, because PEST++ refuses an empty field and ``nan``
        alike (see ``write_named_series_targets``).

        **A lake reads its ONE stage on both routes** (ledger 331). The stage
        table repeats a lake's stage on every cell the lake connects to, and the
        stage-file route used to ADD those copies: the user's pond, 763 connected
        cells, read 401,339.97 ft for a stage of 526.00. It takes the one number
        per record now, as a stage is not added (ledger 297), and refuses
        copies that disagree rather than averaging them.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns. The observation route reads the LAK
            package's declared OBS and the files it names in the model's
            workspace; the stage-file route needs a LAK ``stage_filerecord``.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per target name, in model length
            units: the lake's bottom where it is dry, ``NaN`` where it is
            inactive. Empty (those columns, no rows) when the stage file holds
            no record. A target naming a lake the LAK package does not have
            raises (it got no column on the stage-file route, ledger
            503).

        Raises
        ------
        AttributeError
            When no target is declared and there is no stage file to fall back
            on: the model has no LAK package ("Package 'lak' is not attached to
            model ...") or its LAK declares no ``stage_filerecord``.
        FileNotFoundError
            On the observation route, when a file the OBS declares is not in the
            workspace.
        ValueError
            When a target's lake is not one the LAK package has, or its
            boundname names none or several (ledger 503); when a
            target's label is declared in the LAK's OBS as another
            quantity, for another lake or more than once; when some targets of
            the set are declared and others are not; when the model carries
            several LAK packages and ``package`` does not name one, or names one
            that declares none of the targets; or when a declared file's columns
            do not line up with its records (ledger 344). On the observation
            route, when the file does not END every stress period the TDIS
            declares -- a run that stopped, a truncated copy -- naming the period
            (ledger 346). On the stage-file
            route: when a period the stage file saved lacks its
            final step -- naming the period, the steps saved and the ``SAVE HEAD
            LAST`` setting that fixes it (ledger 336); or when the rows of one
            lake in one saved record are not copies of one stage -- naming the
            lake, the record and the spread. MF6's stage table cannot produce the
            second; it guards a table that is not the one this reads. On the
            stage-file route, when a DRY lake connects ``EMBEDDEDH`` or
            ``EMBEDDEDV``: MF6 takes its bottom from its lake table, which is not
            read (declare a STAGE observation for it instead).

        Examples
        --------
        >>> targets = mf.LakeStageTargets(locations={"pond": 0})
        >>> series = targets.simulated_series(model)   # one row per stress period
        >>> series["pond"].isna()                       # True where the pond is inactive
        """

        refuse_another_clock(self, model, caller="LakeStageTargets.simulated_series(...)")
        observed = self._observed(model)
        if observed is not None:
            return observed
        frame = _calibration_stage_table(model)
        if frame.empty:
            return pd.DataFrame(columns=["time", *self._series.keys()])
        frame = _lake_stage_at_period_end(model, frame)
        # Every lake is one the package has (ledger 503): one it did not
        # have matched no row here and lost its column in silence.
        definitions = self._definitions_on(model)
        merged = definitions.merge(frame, on="lake", how="left").dropna(subset=["per"])
        merged["per"] = merged["per"].astype(int)
        # ``dropna=False`` keeps a record in which every target is dry, and a
        # target dry in EVERY record -- the default drops both, and the second
        # would give a later realization one column more than PEST's base run.
        wide = merged.pivot_table(
            index="per", columns="name", values="stage", aggfunc="first", dropna=False
        ).reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare target lake stages against simulated lake stages.

        ``sim_stage`` comes from the same two routes as :meth:`simulated_series`,
        and is the lake's ONE stage on both: the stage-file route used to add the
        copies the stage table holds for each connected cell (ledger 331). A DRY
        lake's is its BOTTOM on both routes (ledger 351), so a dry period is
        compared -- ``residual`` is the bottom minus the target -- and counted in
        :meth:`stats` alike on both; the stage-file route read ``NaN`` there and
        dropped it (measured: ``n`` 3 against 6 for one run). ``sim_stage`` is
        ``NaN`` only for an INACTIVE lake, and so are ``residual`` and
        ``abs_residual``. Targets are matched on ``per``, read from their ``time`` as a
        zero-based integer -- the stress period, on both routes (ledgers 346, 348),
        at its END. On the observation route a key naming no period of the model is
        refused, naming the period each simulated-time key ends: that route was
        keyed by MF6's time (``10, 20, 30``) until ledger 346, and such a set would
        otherwise read ``NaN``. The stage-file route refuses a key naming a period
        its file does not hold too (ledger 500) -- it read ``NaN``, and
        :meth:`stats` left the row out in silence. That route was keyed by
        saved-record index until ledger 348: at ``nstp=3`` periods 0, 1, 2 read steps
        1, 2, 3 of period 0 (measured), so a set keyed by record index for that route
        -- as this docstring advised -- now matches its first ``nper`` keys to PERIODS
        and refuses the rest; re-key it by period. A schedule that saves no period's
        end is refused (ledger 336).

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row when the targets hold no
            values): ``name``, ``lake``, ``time``, ``per``, ``stage_target``,
            ``sim_stage``, ``residual`` (simulated minus target) and
            ``abs_residual``. With no values ``stage_target`` and the residuals are
            ``NaN``, and ``time`` is the period, equal to ``per``, on both
            routes: the table had no ``per`` -- it was renamed ``time`` -- so
            :meth:`residuals` raised ``KeyError`` (found in ledger 342).

        Raises
        ------
        AttributeError
            As :meth:`simulated_series`: no target is declared and the model has
            no LAK stage file.
        FileNotFoundError
            As :meth:`simulated_series`: a declared observation file is missing.
        ValueError
            As :meth:`simulated_series`: a lake the package does not have,
            a target's label declared as something
            else, a set declared in part, several LAK packages and no
            ``package``, a file that does not line up or one short of a period's
            end; on either route, when a target's key is no stress period the
            route's table holds (ledgers 346, 500); on the stage-file
            route, when a period the stage file saved lacks its final step or a
            lake's copies in one record disagree, or a DRY lake connects
            ``EMBEDDEDH``/``EMBEDDEDV``; or when a target's ``time`` is not a
            whole number, or two rows hold a value for one target in one period
            (the period rule PEST's build reads too, ledger 340).
        """

        refuse_another_clock(self, model, caller="LakeStageTargets.compare(...)")
        definitions = self._definitions_on(model)
        observed = self._observed(model)
        if observed is not None:
            simulated = observed.melt(id_vars=["time"], var_name="name", value_name="sim_stage")
            simulated["per"] = simulated["time"]
            simulated = definitions.merge(simulated, on="name", how="left")
        else:
            simulated = _lake_stage_at_period_end(model, _calibration_stage_table(model))
            simulated = definitions.merge(simulated, on="lake", how="left").rename(columns={"stage": "sim_stage"})
        if self._values.empty:
            compare = _compare_without_values(simulated, target_column="stage_target")
        else:
            # The values joined to the CHECKED lakes, a boundname resolved -- what
            # ``get()`` shows as declared.
            targets = _copy_frame(self._values.merge(definitions, on="name", how="left"))
            targets["per"] = _measured_periods(
                targets, value_column="stage", caller="LakeStageTargets.compare(...)"
            )
            if observed is not None:
                _refuse_keys_off_the_periods(
                    targets, observed["time"], model=model, caller="LakeStageTargets.compare(...)"
                )
            else:
                # Every file route refuses a key past the run's periods too
                # (ledger 500): it read NaN, and stats() left it out.
                _refuse_keys_off_the_periods(
                    targets, simulated["per"].dropna().unique(), model=model,
                    caller="LakeStageTargets.compare(...)", route=_LAKE_STAGE.file_route,
                )
            if "stage" in targets.columns and "stage_target" not in targets.columns:
                targets = targets.rename(columns={"stage": "stage_target"})
            compare = targets.merge(simulated, on=["name", "lake", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_stage"] - compare["stage_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table for lake-stage targets.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``lake``, ``time``,
            ``per``, ``stage_target``, ``sim_stage``, ``residual`` and
            ``abs_residual``.
            A set given no ``values`` -- a capture template, a forecast set --
            gets one row per simulated row, its targets and residuals ``NaN``:
            it raised ``KeyError`` naming ``per`` on both routes until ledger
            342's fix for ``DrnFlowTargets`` reached this family.
        """

        frame = self.compare(model)
        return frame.loc[:, ["name", "lake", "time", "per", "stage_target", "sim_stage", "residual", "abs_residual"]]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the lake-stage target set.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual -- an INACTIVE lake has none; a
            dry one is compared at its bottom), ``mean_error``, ``mae`` and ``rmse`` of simulated minus
            target stage, in model length units. ``n`` is 0 and the rest
            ``NaN`` when no row has a residual.
        """

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ) -> dict[str, list[tuple[str, str, int]]]:
        """Return a FloPy ``continuous`` dict for MF6 lake-stage observations.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"lake_stage_targets.csv"``. Any name: :meth:`simulated_series`
            reads the file the LAK package's OBS declares, whatever it is called
            (ledger 344).
        kind : str, optional
            The MF6 LAK observation type written on every record. Must be
            ``"STAGE"`` (any case), the default -- the only type the observation
            route reads back -- and any other raises ``ValueError``: it wrote,
            say, a ``VOLUME`` under the targets' names, which that route then
            read as the stage (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "STAGE", lake + 1), ...]}``, one record per
            target: its name with every character that is not a letter or digit
            made ``_`` and leading and trailing ``_`` stripped (``"obs"`` if
            nothing is left), and its lake number ONE-based. The shape
            ``ModflowUtlobs(continuous=)`` and a LAK package's ``observations=``
            take.

        Raises
        ------
        ValueError
            When ``kind`` is not ``"STAGE"``, or two target names would be
            written as one label once MF6 upper-cases it (``"Pond A"`` and
            ``"pond-a"``: MF6 wrote the header ``POND_A,POND_A`` and both targets
            read the first lake).
        """

        kind = _refuse_observation_kind(kind, _LAKE_STAGE)
        _refuse_shared_labels(self._series.keys(), _LAKE_STAGE)
        # A number is written from one, as MF6 numbers lakes; a boundname as it is,
        # which MF6 resolves in the LAK's packagedata.
        records = [
            (_observation_label(name), kind, lake_id if isinstance(lake_id, str) else int(lake_id) + 1)
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
        """Attach an MF6 utility observation package for lake-stage targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase
            The model whose FloPy LAK package receives the OBS package unless
            ``package`` is given: the one the targets name (their ``package``
            field, ``ValueError`` when the model has no LAK package of that
            name), else the model's only LAK package, found by type (several
            raise ``ValueError`` naming them and ``package=``, none
            ``AttributeError``; ledger 352). The records themselves need no model.
        pname : str, optional
            FloPy package name. Default ``"lak_obs"``.
        filename : str, optional
            The OBS input file name FloPy writes, in the owning model's
            directory. Default ``"lake_stage_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"lake_stage_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 LAK observation type, as in :meth:`to_flopy_obs`: ``"STAGE"``,
            the default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            the one ``model`` resolves (see ``model``). Default ``None``. Reading
            the observations back needs the targets' own ``package`` field to
            name it when the model carries several LAK packages.

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
            owner = _owning_package(model, _LAKE_STAGE, self.package)
        csv_name = _default_obs_csv_name(filename, "lake_stage_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

    def calibration_plot(self, model, *, type: str = "calibration"):
        """Build a calibration plot directly from these targets and one model.

        Parameters
        ----------
        model : SimulationBase
            A run model with a LAK package; the plot is drawn from
            :meth:`compare`, whose ``Raises`` apply.
        type : str, optional
            Which plot:

            - ``"calibration"`` (default) or ``"obs_vs_sim"``: target stage
              against simulated stage, one point per row, titled "Observed vs
              simulated lake stage", without fit statistics;
            - ``"residuals_by_period"``: mean absolute residual per ``time``.

            ``"heads"`` RAISES ``ValueError`` (the lake comparison has no
            ``layer`` column, which that plot requires), and any other string
            draws the cross plot WITH fit statistics.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

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
    """Model-bound lake-stage helper returned from ``model.targets.lake_stage``.

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from.
    """

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
        """The wide simulated series (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model with a LAK package to read. Default ``None`` uses the
            bound model; another run reads that run's files for the same
            targets. Read by :meth:`LakeStageTargets.simulated_series`'s two
            routes -- the LAK's declared STAGE observations when every target is
            declared there, else the LAK stage file -- one row per stress period
            at its END, whose ``Raises`` apply.

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
            The run model with a LAK package to compare against. Default
            ``None`` uses the bound model; another run compares the same targets
            against it -- what the plots' ``baseline=`` does. Read as
            :meth:`LakeStageTargets.compare` reads it (each target's ``time`` a
            zero-based stress period, a DRY lake at its bottom, ledger 351),
            whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row, ``stage_target`` and the
            residuals ``NaN``, when the set holds no values): ``name``, ``lake``,
            ``time``, ``per``, ``stage_target``, ``sim_stage``, ``residual``
            (simulated minus target) and ``abs_residual``.
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
            :meth:`compare`'s rows, reduced to ``name``, ``lake``, ``time``,
            ``per``, ``stage_target``, ``sim_stage``, ``residual`` and
            ``abs_residual``.
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
            (:meth:`LakeStageTargets.stats`).
        """

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "lake_stage_targets.csv",
        kind: str = "STAGE",
    ):
        """A FloPy ``continuous`` observation dict for the lake-stage targets.

        The records need no model -- each names its lake by number -- so the
        bound model takes no part.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"lake_stage_targets.csv"``. Any name: the observation route
            reads whichever file the LAK package's OBS declares (ledger 344).
        kind : str, optional
            The MF6 LAK observation type written on every record. Must be
            ``"STAGE"`` (any case), the default; any other raises ``ValueError``:
            declared under the targets' names, the observation route would refuse
            it as another quantity when read back (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "STAGE", lake + 1), ...]}``, the lake number
            ONE-based, as :meth:`LakeStageTargets.to_flopy_obs` returns it (whose
            ``Raises`` apply).
        """

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
        """Attach an MF6 lake-stage observation package to the bound model.

        Added to the bound model's FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        pname : str, optional
            FloPy package name of the OBS package. Default ``"lak_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory; its stem names the default ``csv_name``. Default
            ``"lake_stage_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"lake_stage_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 LAK observation type, as in :meth:`to_flopy_obs`:
            ``"STAGE"``, the default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to. Default
            ``None`` takes the bound model's LAK package: the one the targets
            name (their ``package`` field; ``ValueError`` when the model has no
            LAK package of that name), else its only LAK package, found by type
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

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the lake-stage targets on the bound model.

        Parameters
        ----------
        type : str, optional
            Which plot, drawn from :meth:`compare` on the bound model (whose
            ``Raises`` apply):

            - ``"calibration"`` (default) or ``"obs_vs_sim"``: target stage
              against simulated stage, one point per row, titled "Observed vs
              simulated lake stage", without fit statistics;
            - ``"residuals_by_period"``: mean absolute residual per ``time``.

            ``"heads"`` RAISES ``ValueError`` (the lake comparison has no
            ``layer`` column, which that plot requires), and any other string
            draws the cross plot WITH fit statistics.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

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
