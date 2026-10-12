"""Dated readings binned onto stress periods: the one rule behind ``from_hydrographs``.

A target set is keyed by ZERO-BASED STRESS PERIOD (ledger 340), and a date is
refused as a key (ledger 499). A hydrograph table is keyed by DATE. This module
is the one place a date becomes a period, shared by
:meth:`HeadTargets.from_hydrographs <myflopy.modflow.mf6.observations.heads.HeadTargets.from_hydrographs>`
(and ``ConcTargets``, which inherits it) and
:meth:`LakeStageTargets.from_hydrographs <myflopy.modflow.mf6.observations.lake.LakeStageTargets.from_hydrographs>`,
so a well and a lake bin by the same rule (ledger 543).

The dates of a period are the model's :class:`~myflopy.modflow.mf6.period_calendar.PeriodCalendar`
(ledger 529) -- its TDIS ``START_DATE_TIME``, ``TIME_UNITS`` and period lengths.
Two modes, the user's answer to Q4:

* ``mode="history"`` -- the model runs the dated record. ``rule="nearest_end"``
  (default) takes, for each period, the reading NEAREST its END, because the
  simulated value a target is compared with is the end-of-period value (ledger
  336); ``"median"``/``"mean"`` reduce the readings inside the period.
* ``mode="average_year"`` -- the model is an ``mf.AverageYear`` cycle (ledger 530)
  and the readings come from many years. Each calendar month's target is the
  MEDIAN over years of each year's monthly median (``rule="mean"``: the mean of
  yearly means), placed on that month's period of the LAST cycle, the one targets
  are compared on.

Binning is EAGER against one model, and the result carries that model's
calendar as a CLOCK STAMP: :func:`refuse_another_clock` refuses the set wherever
it meets a model whose periods are other dates (another cycle count, another
start), because a period number would then name another date in silence.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package
from myflopy.modflow.mf6.observations._shared import ROW_WEIGHT_COLUMN, _long_values
from myflopy.modflow.mf6.period_calendar import (
    PeriodCalendar,
    format_start,
    parse_start,
    period_calendar,
)

logger = get_logger(__name__)

#: The binning modes, and each one's rules (the first is the default).
RULES = {
    "history": ("nearest_end", "median", "mean"),
    "average_year": ("median", "mean"),
}

#: The report's columns, in order (:attr:`Binning.report`).
REPORT_COLUMNS = (
    "name", "per", "month", "period_start", "period_end", "n_readings", "n_years",
    "first_reading", "last_reading", "reading", "offset_days", "value", "measured",
)


@dataclass(frozen=True)
class Binning:
    """How a set of dated readings became period-keyed target values.

    Returned by :func:`bin_readings`; a ``from_hydrographs`` target set carries
    its ``report`` as ``.binning`` and its ``calendar`` as ``.calendar``.

    Attributes
    ----------
    values : pandas.DataFrame
        The targets: ``time`` (the zero-based stress period), ``name`` and
        ``value``, one row per (name, period) that has a target -- a period with no
        reading, and an average-year month under ``min_years``, have none.
    report : pandas.DataFrame
        One row per (name, period) the readings reached -- see :func:`bin_readings`.
    calendar : PeriodCalendar
        The model's calendar the readings were binned against: the CLOCK STAMP.
    mode : str
        ``"history"`` or ``"average_year"``.
    rule : str
        The rule applied: ``"nearest_end"``, ``"median"`` or ``"mean"``.
    options : dict
        ``max_offset`` (a ``pandas.Timedelta`` or ``None``), ``between`` (two
        timestamps or ``None`` each, or ``None``) and ``min_years`` (an ``int``, or
        ``None`` in history mode), as applied.
    """

    values: pd.DataFrame
    report: pd.DataFrame
    calendar: PeriodCalendar
    mode: str
    rule: str
    options: dict

    def record(self) -> dict:
        """What the PEST build writes to its metadata: mode, rule and options as text."""

        between = self.options.get("between")
        offset = self.options.get("max_offset")
        return {
            "mode": self.mode,
            "rule": self.rule,
            "max_offset": None if offset is None else str(offset),
            "between": None if between is None else [
                None if stamp is None else stamp.isoformat() for stamp in between
            ],
            "min_years": self.options.get("min_years"),
        }


def _warn(message: str) -> None:
    """Say ``message`` through ``warnings.warn`` AND the log (ledger 288: a log line alone reaches nobody)."""

    logger.warning(message)
    warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())


def _listing(items, limit: int = 8) -> str:
    items = list(items)
    shown = ", ".join(str(item) for item in items[:limit])
    return shown + (f", and {len(items) - limit} more" if len(items) > limit else "")


def read_readings(
    readings: Any,
    *,
    date_column: str = "date",
    value_column: str,
    owner: str,
) -> pd.DataFrame:
    """A hydrograph table, WIDE or LONG, as one long table of dated readings.

    Parameters
    ----------
    readings : str or pathlib.Path or pandas.DataFrame or pandas.Series
        The readings. A path is read with ``pandas.read_csv``. A table is:

        - WIDE: ``date_column`` plus one column per location, headed by its name;
        - LONG: ``date_column``, ``name`` and ``value_column`` (both headers matched
          ignoring case), one row per reading;
        - ONE SERIES: ``date_column`` and ``value_column`` only, for a single
          location (the caller supplies its name).

        A table whose dates are its INDEX (a ``DatetimeIndex``, or an index named
        ``date_column``) and that has no ``date_column`` column is read with the
        index as the dates; a ``Series`` is one location's readings indexed by date.
        A blank value is no reading; any other entry that is not a finite number
        raises ``ValueError`` naming it (ledger 340's rule).
    date_column : str, default "date"
        The column holding each row's date, matched exactly. Every entry must be a
        date pandas parses with a four-digit year (``"2018-01-15"``,
        ``"2018-01-15 08:30"``, a ``Timestamp``), with no time zone -- MF6's
        ``START_DATE_TIME`` is a wall-clock string (ledger 529). Text written
        ``d/m/y`` or ``m/d/y`` is read MONTH first, as pandas reads it; a column
        written day first in some rows (``"13/02/2018"``) and read month first in
        others (``"01/02/2018"``) raises ``ValueError`` -- parse such a column
        yourself (``pandas.to_datetime(..., dayfirst=True)``) and pass the table
        (ledger 545).
    value_column : str
        The value column of a LONG or one-series table, matched ignoring case;
        unused by a WIDE table.
    owner : str
        The target class, named in every refusal (``"HeadTargets"``).

    Returns
    -------
    pandas.DataFrame
        ``name`` (text, absent for one series), ``date`` (timezone-naive
        ``Timestamp``) and ``value`` (float), one row per non-blank reading. Two
        readings of one location at one instant holding the SAME value are one row.

    Raises
    ------
    ValueError
        When the table has no ``date_column`` and no date index; when a date does
        not parse, has no four-digit year or carries a time zone (naming up to five
        rows); when text dates mix day-first and month-first rows (naming one of
        each); when a value is neither blank nor a finite number; when two readings
        of one location at one instant hold DIFFERENT values (naming them); or when
        the table carries a ``weight`` column -- several readings become one target,
        so a reading's own weight has nothing to attach to.
    TypeError
        When ``readings`` is none of the accepted types.
    """

    caller = f"{owner}.from_hydrographs(...)"
    if isinstance(readings, pd.Series):
        name = readings.name
        frame = readings.rename(value_column).rename_axis(date_column).reset_index()
        if name is not None and not isinstance(name, (int, np.integer)):
            frame = frame.rename(columns={value_column: str(name)})
    elif isinstance(readings, pd.DataFrame):
        frame = readings.copy()
    elif isinstance(readings, (str, bytes)) or hasattr(readings, "__fspath__"):
        frame = pd.read_csv(readings)
    else:
        raise TypeError(
            f"{caller} reads `readings` as a CSV path, a DataFrame or a Series; got "
            f"{type(readings).__name__}."
        )
    if date_column not in frame.columns:
        index = frame.index
        if isinstance(index, pd.DatetimeIndex) or index.name == date_column:
            frame = frame.rename_axis(date_column).reset_index()
        else:
            raise ValueError(
                f"{caller}: the readings have no {date_column!r} column (columns: "
                f"{_listing(frame.columns)}); name the column holding each row's date "
                "with date_column=."
            )
    lower = {str(column).strip().lower(): column for column in frame.columns}
    weight = lower.get(ROW_WEIGHT_COLUMN)
    if weight is not None and weight not in (date_column,) and str(weight).lower() != value_column.lower():
        raise ValueError(
            f"{caller}: the readings carry a {weight!r} column, and binning makes ONE "
            "target from several readings, so a reading's own weight has no "
            "observation to go to. Weight a location through its locations table, or "
            "the whole set through its observation spec (cal.observe(..., weight=) for "
            "a lake set)."
        )

    dates = frame[date_column]
    parsed: dict[Any, pd.Timestamp] = {}
    bad: list[str] = []
    for position, raw in enumerate(dates.tolist()):
        key = raw if _hashable(raw) else repr(raw)
        if key in parsed:
            continue
        try:
            parsed[key] = parse_start(raw, what=f"date {raw!r}")
        except ValueError as error:
            if len(bad) < 5:
                bad.append(f"row {position}: {error}")
            parsed[key] = pd.NaT
    failed = [
        position for position, raw in enumerate(dates.tolist())
        if parsed[raw if _hashable(raw) else repr(raw)] is pd.NaT
    ]
    if failed:
        raise ValueError(
            f"{caller}: {len(failed)} row(s) of {date_column!r} are not a date it can place "
            "on the model's calendar -- " + "; ".join(bad)
            + (f"; and {len(failed) - len(bad)} more" if len(failed) > len(bad) else "")
            + ". Give each row a date with its year and no time zone, e.g. '2018-01-15'."
        )
    _refuse_mixed_day_order(dates.tolist(), parsed, column=date_column, caller=caller)
    frame[date_column] = [parsed[raw if _hashable(raw) else repr(raw)] for raw in dates.tolist()]
    long = _long_values(
        frame, time_column=date_column, value_column=value_column,
        output_column="value", owner=owner,
    ).rename(columns={"time": "date"})
    long = long.loc[long["value"].notna()].copy()
    long["date"] = pd.DatetimeIndex(long["date"])
    keys = ["name", "date"] if "name" in long.columns else ["date"]
    spread = long.groupby(keys, sort=True)["value"].agg(["min", "max"])
    clash = spread.loc[spread["min"] != spread["max"]]
    if not clash.empty:
        listing = [
            (f"{index[0]!s} at {index[1].isoformat()}" if len(keys) == 2 else index.isoformat())
            + f" ({row['min']:g} and {row['max']:g})"
            for index, row in clash.head(5).iterrows()
        ]
        raise ValueError(
            f"{caller}: {len(clash)} reading(s) are given twice at one instant with "
            f"different values: {'; '.join(listing)}"
            + (f"; and {len(clash) - 5} more" if len(clash) > 5 else "")
            + ". Which one is right is a question about the data; keep one."
        )
    return long.drop_duplicates(subset=keys).reset_index(drop=True)


#: A text date written day-month-year or month-day-year with a four-digit year:
#: ``"13/02/2018"``, ``"2.13.2018"``, ``"01-02-2018 08:30"``.
_NUMERIC_DMY = re.compile(r"^\s*(\d{1,2})[/.\-](\d{1,2})[/.\-]\d{4}(?!\d)")


def _refuse_mixed_day_order(raws: list, parsed: Mapping, *, column: str, caller: str) -> None:
    """Refuse a text date column that is day-first in some rows and read month-first in others.

    Each entry is parsed on its own (:func:`parse_start`), and pandas reads an
    ambiguous ``"01/02/2018"`` MONTH first while ``"13/02/2018"`` can only be day
    first -- so a column written day first came back with every day 1-12 in the
    wrong month and every day 13-31 right, in silence (ledger 545).
    A column whose text dates admit ONE order -- month first throughout, day first
    throughout with no ambiguous row, or ISO ``YYYY-MM-DD`` -- passes; a column of
    only ambiguous rows reads month first, as pandas does, and cannot be told apart.
    """

    day_first: tuple[int, str] | None = None
    month_first: tuple[int, str, pd.Timestamp] | None = None
    for position, raw in enumerate(raws):
        if not isinstance(raw, str):
            continue
        match = _NUMERIC_DMY.match(raw)
        if match is None:
            continue
        first, second = int(match.group(1)), int(match.group(2))
        if first > 12 >= second:
            day_first = day_first or (position, raw)
        elif first <= 12 and first != second:
            month_first = month_first or (position, raw, parsed[raw])
        if day_first is not None and month_first is not None:
            break
    if day_first is None or month_first is None:
        return
    raise ValueError(
        f"{caller}: {column!r} is written DAY first in some rows (row {day_first[0]}: "
        f"{day_first[1]!r}) and is read MONTH first in others (row {month_first[0]}: "
        f"{month_first[1]!r} -> {month_first[2].date().isoformat()}); a date such as "
        "'01/02/2018' reads either way, and binning it to the wrong month would be silent. "
        f"Parse the column yourself and pass the table -- readings[{column!r}] = "
        f"pandas.to_datetime(readings[{column!r}], dayfirst=True) -- or write the dates "
        "as YYYY-MM-DD."
    )


def _hashable(value: Any) -> bool:
    try:
        hash(value)
    except TypeError:
        return False
    return True


def _parse_offset(max_offset: Any, *, caller: str) -> pd.Timedelta | None:
    if max_offset is None:
        return None
    if not isinstance(max_offset, np.timedelta64) and isinstance(
        max_offset, (bool, np.bool_, int, float, np.number)
    ):
        # numpy's timedelta64 subclasses its signed integers, and it HAS a unit.
        raise ValueError(
            f"{caller}: max_offset={max_offset!r} has no unit (pandas would read a bare "
            "number as nanoseconds); give a duration such as '10D' or '36h'."
        )
    try:
        offset = pd.Timedelta(max_offset)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{caller}: max_offset={max_offset!r} is not a duration ({error}); give one "
            "such as '10D' or '36h'."
        ) from None
    if offset is pd.NaT or offset <= pd.Timedelta(0):
        raise ValueError(f"{caller}: max_offset={max_offset!r} must be a positive duration.")
    return offset


def _parse_between(between: Any, *, caller: str) -> tuple[pd.Timestamp | None, pd.Timestamp | None] | None:
    if between is None:
        return None
    if isinstance(between, (str, bytes)) or not isinstance(between, (tuple, list)) or len(between) != 2:
        raise ValueError(
            f"{caller}: between= takes two dates, (first, last), either of them None for "
            f"open; got {between!r}."
        )
    first, last = (
        None if value is None else parse_start(value, what=f"{caller} between=")
        for value in between
    )
    if first is not None and last is not None and first > last:
        raise ValueError(
            f"{caller}: between=({format_start(first)}, {format_start(last)}) ends before it starts."
        )
    return first, last


def _resolve_options(mode, rule, max_offset, between, min_years, *, caller: str):
    """``(rule, options)``, refusing a mode, a rule or an option the mode does not read."""

    if mode not in RULES:
        raise ValueError(f"{caller}: mode={mode!r}; it is 'history' or 'average_year'.")
    if rule is None:
        rule = RULES[mode][0]
    if rule not in RULES[mode]:
        if mode == "average_year" and rule == "nearest_end":
            raise ValueError(
                f"{caller}: rule='nearest_end' picks one reading at a period's END, and an "
                "average-year target is a statistic over many years' readings of its month; "
                "use rule='median' (the default) or 'mean'."
            )
        raise ValueError(
            f"{caller}: rule={rule!r}; mode={mode!r} takes "
            + ", ".join(repr(name) for name in RULES[mode]) + "."
        )
    if mode == "history":
        given = [name for name, value in (("between", between), ("min_years", min_years)) if value is not None]
        if given:
            raise ValueError(
                f"{caller}: {', '.join(given)} {'is' if len(given) == 1 else 'are'} read in "
                "mode='average_year' only; in history mode the model's own span decides "
                "which readings count."
            )
        if max_offset is not None and rule != "nearest_end":
            raise ValueError(
                f"{caller}: max_offset limits how far the reading nearest a period's end "
                f"may be from it, so it acts under rule='nearest_end' only, not {rule!r}."
            )
        return rule, {"max_offset": _parse_offset(max_offset, caller=caller), "between": None,
                      "min_years": None}
    if max_offset is not None:
        raise ValueError(
            f"{caller}: max_offset is read in mode='history' with rule='nearest_end' only; "
            "an average-year target is a statistic over each month's readings."
        )
    if min_years is None:
        min_years = 1
    if isinstance(min_years, (bool, np.bool_)) or not isinstance(min_years, (int, np.integer)) or min_years < 1:
        raise ValueError(f"{caller}: min_years={min_years!r} must be a whole number from 1.")
    return rule, {"max_offset": None, "between": _parse_between(between, caller=caller),
                  "min_years": int(min_years)}


def bin_readings(
    readings: pd.DataFrame,
    calendar: PeriodCalendar,
    *,
    mode: str = "history",
    rule: str | None = None,
    max_offset: Any = None,
    between: Any = None,
    min_years: int | None = None,
    months: Mapping[int, int] | None = None,
    owner: str = "HeadTargets",
) -> Binning:
    """Bin dated readings onto a model's stress periods, by one stated rule.

    Parameters
    ----------
    readings : pandas.DataFrame
        Long readings as :func:`read_readings` returns them: ``name``, ``date``
        (timezone-naive) and ``value`` (finite), one row per reading.
    calendar : PeriodCalendar
        The model's calendar (``model.calendar``): a period's start, end and
        midpoint come from it, and it becomes the result's clock stamp.
    mode : {"history", "average_year"}, default "history"
        ``"history"``: the model runs the dated record, and a period's target comes
        from the readings at its dates. ``"average_year"``: the model is an
        average-year cycle; a calendar month's target comes from that month's
        readings in every year, and lands on the month's period of the LAST cycle
        (``months``).
    rule : str or None, default None
        History: ``"nearest_end"`` (the default) -- every reading in the model's
        span, ``(start, last end]``, goes to its NEAREST period end (a reading
        midway between two goes to the end of the period holding it), and each
        period keeps the reading nearest its end (a tie goes to the EARLIER
        reading, the one inside the period); ``"median"``/``"mean"`` -- the
        statistic of the readings inside the period, ``(start, end]``.
        Average year: ``"median"`` (the default) -- each year's monthly median,
        then the median of those over the years; ``"mean"`` -- the same with means.
        ``"nearest_end"`` in average-year mode raises ``ValueError`` naming
        ``rule='median'``.
    max_offset : str or pandas.Timedelta or datetime.timedelta or numpy.timedelta64 or None, default None
        History with ``rule="nearest_end"`` only: a reading farther than this from
        its period's end is not used, so a period with no reading that near has no
        target. ``"10D"``, ``"36h"``, a ``numpy.timedelta64``; a bare number
        (nanoseconds to pandas) or a duration that is not positive raises.
        ``None``: no limit.
    between : tuple of two dates or None, default None
        Average year only: ``(first, last)``, the readings used, both ends
        included; either may be ``None`` for open. Dates parse as ``readings``'
        do, and a date is its midnight: ``last="2005-12-31"`` leaves out a reading
        at 08:00 that day (give ``"2005-12-31 23:59:59"`` or the next day's date
        to keep it). ``None``: every reading.
    min_years : int or None, default None
        Average year only: a month read in fewer distinct calendar years than this
        has no target (it stays in the report, ``measured`` False). ``None``
        means 1 in average-year mode; a whole number from 1.
    months : mapping of int to int or None, default None
        Average year: ``{calendar month: zero-based period}`` of the cycle the
        targets land on -- ``AverageYear.months_of(-1)``; required in that mode
        and must name the twelve months. Ignored in history mode.
    owner : str, default "HeadTargets"
        The target class named in every message.

    Returns
    -------
    Binning
        ``values`` (``time``, ``name``, ``value``), ``report`` and the stamp.
        ``report`` has one row per (name, period) the readings reached:
        ``name``, ``per``, ``month`` (the period's midpoint month),
        ``period_start``, ``period_end``, ``n_readings`` (readings binned to it),
        ``n_years`` (their distinct calendar years), ``first_reading`` and
        ``last_reading``, ``reading`` (the one used, ``nearest_end`` only, else
        ``NaT``), ``offset_days`` (that reading minus the period end, in days;
        negative before it), ``value`` and ``measured``.

    Raises
    ------
    ValueError
        For a mode, rule or option the mode does not read (each named); in
        history mode when no reading falls in the model's span (printing both
        spans), or none survives ``max_offset``; in average-year mode when
        ``months`` is missing or does not name the twelve months, or no reading
        falls in ``between``.

    Warns
    -----
    UserWarning
        History: when readings fall outside the model's span (counted per name)
        and when a name has no reading left -- through ``warnings.warn`` and the
        log (ledger 288). Average year: when a name has no reading in ``between``.
    """

    caller = f"{owner}.from_hydrographs(...)"
    rule, options = _resolve_options(mode, rule, max_offset, between, min_years, caller=caller)
    frame = readings.loc[:, ["name", "date", "value"]].copy()
    frame["date"] = pd.DatetimeIndex(frame["date"])
    names = list(dict.fromkeys(frame["name"]))
    starts, ends = calendar.starts, calendar.ends
    if mode == "history":
        report = _bin_history(frame, calendar, rule=rule, max_offset=options["max_offset"], caller=caller)
    else:
        if months is None or sorted(months) != list(range(1, 13)):
            raise ValueError(
                f"{caller}: mode='average_year' needs the twelve months' periods of the last "
                f"cycle; got {months!r}."
            )
        report = _bin_average_year(
            frame, months, rule=rule, between=options["between"],
            min_years=options["min_years"], caller=caller,
        )
    used = set(report["name"]) if not report.empty else set()
    silent = [name for name in names if name not in used]
    if silent:
        _warn(
            f"{caller}: {len(silent)} of {len(names)} location(s) have no reading the "
            f"binning could use and get no target: {_listing(silent)}."
        )
    per = report["per"].to_numpy(dtype=int)
    report["month"] = calendar.months[per]
    report["period_start"] = starts[per]
    report["period_end"] = ends[per]
    report = report.loc[:, list(REPORT_COLUMNS)].sort_values(["name", "per"], kind="stable")
    report = report.reset_index(drop=True)
    measured = report.loc[report["measured"]]
    if measured.empty:
        raise ValueError(
            f"{caller}: no month was read in at least min_years={options['min_years']} "
            "calendar year(s), so no target is left; lower min_years or widen between=."
        )
    values = pd.DataFrame(
        {"time": measured["per"].astype(int), "name": measured["name"].astype(str),
         "value": measured["value"].astype(float)}
    ).reset_index(drop=True)
    logger.info(
        "%s: %d reading(s) of %d location(s) binned (%s, %s) into %d target(s).",
        caller, len(frame), len(names), mode, rule, len(values),
    )
    return Binning(values=values, report=report, calendar=calendar, mode=mode, rule=rule, options=options)


def _bin_history(frame: pd.DataFrame, calendar: PeriodCalendar, *, rule: str, max_offset, caller: str) -> pd.DataFrame:
    """History: each period's reading nearest its end, or its readings' statistic."""

    start, ends = calendar.start, calendar.ends
    last = ends[-1]
    inside = (frame["date"] > start) & (frame["date"] <= last)
    if not inside.any():
        raise ValueError(
            f"{caller}: no reading falls in the model's span, {start.isoformat()} to "
            f"{last.isoformat()} (the readings run {frame['date'].min().isoformat()} to "
            f"{frame['date'].max().isoformat()}). Bin them against the model that runs "
            "those dates -- or, for readings from many years against an average-year "
            "model, use mode='average_year'."
        )
    outside = frame.loc[~inside]
    if not outside.empty:
        counts = outside.groupby("name", sort=True).size()
        _warn(
            f"{caller}: {len(outside)} reading(s) fall outside the model's span "
            f"({start.isoformat()} to {last.isoformat()}) and are left out: "
            + _listing(f"{name} ({count})" for name, count in counts.items()) + "."
        )
    frame = frame.loc[inside].copy()
    ends_ns = ends.asi8
    stamps = frame["date"].to_numpy(dtype="datetime64[ns]").astype("int64")
    after = np.searchsorted(ends_ns, stamps, side="left")  # first end at or after the reading
    if rule == "nearest_end":
        before = np.maximum(after - 1, 0)
        to_after = ends_ns[after] - stamps
        to_before = np.where(after > 0, stamps - ends_ns[before], np.iinfo(np.int64).max)
        # Midway between two ends goes to the end of the period holding it.
        per = np.where(to_after <= to_before, after, before)
        frame["per"] = per
        frame["offset"] = stamps - ends_ns[per]
        if max_offset is not None:
            far = np.abs(frame["offset"].to_numpy()) > max_offset.value
            if far.any():
                logger.info("%s: %d reading(s) farther than %s from their period end are not used.",
                            caller, int(far.sum()), max_offset)
            frame = frame.loc[~far]
            if frame.empty:
                raise ValueError(
                    f"{caller}: no reading lies within max_offset={max_offset} of a period end."
                )
        # Nearest first; a tie goes to the EARLIER reading, the one inside the period.
        frame["distance"] = np.abs(frame["offset"].to_numpy())
        ordered = frame.sort_values(["name", "per", "distance", "date"], kind="stable")
        groups = ordered.groupby(["name", "per"], sort=False)
        chosen = groups.head(1).set_index(["name", "per"])
        summary = groups.agg(
            n_readings=("value", "size"), n_years=("date", lambda d: d.dt.year.nunique()),
            first_reading=("date", "min"), last_reading=("date", "max"),
        )
        summary["reading"] = chosen["date"]
        summary["offset_days"] = chosen["offset"] / 86_400e9
        summary["value"] = chosen["value"]
    else:
        frame["per"] = after
        groups = frame.groupby(["name", "per"], sort=False)
        summary = groups.agg(
            n_readings=("value", "size"), n_years=("date", lambda d: d.dt.year.nunique()),
            first_reading=("date", "min"), last_reading=("date", "max"),
            value=("value", rule),
        )
        summary["reading"] = pd.NaT
        summary["offset_days"] = np.nan
    summary["measured"] = True
    return summary.reset_index()


def _bin_average_year(frame, months, *, rule: str, between, min_years: int, caller: str) -> pd.DataFrame:
    """Average year: per calendar month, the statistic over years of each year's monthly statistic."""

    if between is not None:
        first, last = between
        keep = pd.Series(True, index=frame.index)
        if first is not None:
            keep &= frame["date"] >= first
        if last is not None:
            keep &= frame["date"] <= last
        if not keep.any():
            raise ValueError(
                f"{caller}: no reading falls in between=({_shown_date(first)}, {_shown_date(last)}); "
                f"the readings run {frame['date'].min().isoformat()} to "
                f"{frame['date'].max().isoformat()}."
            )
        dropped = int((~keep).sum())
        if dropped:
            logger.info("%s: %d reading(s) outside between= are not used.", caller, dropped)
        frame = frame.loc[keep]
    frame = frame.assign(year=frame["date"].dt.year, month=frame["date"].dt.month)
    yearly = frame.groupby(["name", "month", "year"], sort=True)["value"].agg(rule)
    monthly = yearly.groupby(level=["name", "month"]).agg([rule, "size"])
    spans = frame.groupby(["name", "month"]).agg(
        n_readings=("value", "size"), first_reading=("date", "min"), last_reading=("date", "max"),
    )
    summary = spans.join(monthly.rename(columns={rule: "value", "size": "n_years"}))
    summary["measured"] = summary["n_years"] >= int(min_years)
    summary.loc[~summary["measured"], "value"] = np.nan
    summary = summary.reset_index()
    summary["per"] = summary["month"].map(lambda month: int(months[int(month)]))
    summary["reading"] = pd.NaT
    summary["offset_days"] = np.nan
    return summary.drop(columns=["month"])


def _shown_date(stamp) -> str:
    return "open" if stamp is None else format_start(stamp)


def model_clock(model: Any, *, mode: str, owner: str) -> tuple[PeriodCalendar, dict[int, int] | None]:
    """The calendar to bin against and, in average-year mode, the last cycle's months.

    Parameters
    ----------
    model : SimulationBase
        The model the targets will be compared on -- ``run.model(name)`` or
        ``mf.load_mf6_run(ws)``.
    mode : str
        ``"history"`` or ``"average_year"``.
    owner : str
        The target class named in every message.

    Returns
    -------
    tuple
        ``(model.calendar, months)``: ``months`` is
        ``mf.AverageYear.from_model(model).months_of(-1)`` in average-year mode and
        ``None`` in history mode.

    Raises
    ------
    ValueError
        When the model's TDIS cannot date its periods (no ``START_DATE_TIME``, a
        ``TIME_UNITS`` MF6 reads as undefined -- naming ``mf.tdis(...)``); in
        average-year mode when the model is not an ``mf.AverageYear`` cycle, naming
        the period that breaks the pattern and ``mf.AverageYear``.
    """

    caller = f"{owner}.from_hydrographs(...)"
    try:
        calendar = period_calendar(model)
    except ValueError as error:
        raise ValueError(f"{caller} dates the stress periods from the model's TDIS: {error}") from None
    if mode != "average_year":
        return calendar, None
    # cycles.py imports the observation layer, so it cannot be imported at the top of
    # this module; AverageYear.from_model is the ONE reading of a run's cycle.
    from myflopy.cycles import AverageYear

    try:
        year = AverageYear.from_model(model)
    except ValueError as error:
        raise ValueError(
            f"{caller}: mode='average_year' places each calendar month's target on that "
            f"month's period of the model's LAST cycle, and this model is not an average-year "
            f"cycle ({error}). Build its time setup with mf.AverageYear(...) -- year.tdis() "
            "and the year's mf.sto flags -- or bin the record with mode='history'."
        ) from None
    return calendar, year.months_of(-1)


def clock_difference(stamp: PeriodCalendar, on_model: PeriodCalendar) -> str | None:
    """What makes ``on_model``'s periods other dates than ``stamp``'s, in words, or ``None``.

    Parameters
    ----------
    stamp : PeriodCalendar
        The calendar a target set was binned against.
    on_model : PeriodCalendar
        The calendar of the model the set meets.

    Returns
    -------
    str or None
        The first difference -- the period count, the start, the unit or a
        period's length -- or ``None`` when every period has the same dates
        (:meth:`PeriodCalendar.same_clock`; the steady flags are not compared).
    """

    if on_model.nper != stamp.nper:
        return f"this model runs {on_model.nper} stress period(s) where the targets were binned on {stamp.nper}"
    if on_model.start != stamp.start:
        return (
            f"this model starts {on_model.start.isoformat()} where the targets' model started "
            f"{stamp.start.isoformat()}"
        )
    if on_model.time_units != stamp.time_units:
        return f"this model counts time in {on_model.time_units} where the targets' counted {stamp.time_units}"
    if not on_model.same_clock(stamp):
        gaps = np.abs(np.asarray(on_model.perlen) - np.asarray(stamp.perlen))
        per = int(np.argmax(gaps))
        return (
            f"this model's period {per} is {on_model.perlen[per]:g} {on_model.time_units} long "
            f"where the targets' was {stamp.perlen[per]:g}"
        )
    return None


def refuse_another_clock(targets: Any, model: Any, *, caller: str) -> None:
    """Refuse a ``from_hydrographs`` target set meeting a model on another clock.

    A set built from dated readings is keyed by period NUMBER, binned against one
    model's calendar. On a model whose periods are other dates -- an average-year
    run of 5 cycles where the set was binned on 3, a TDIS starting another day --
    the same number names another date, and the comparison would be silently about
    the wrong time. A set not built from dates (no stamp) passes.

    Parameters
    ----------
    targets : object
        A target set; the :class:`Binning` a ``from_hydrographs`` set carries
        (``_binned``) holds the stamp, and a set without one passes.
    model : SimulationBase
        The model it meets.
    caller : str
        What is reading it, named in the error (``"HeadTargets.compare(...)"``).

    Raises
    ------
    ValueError
        When the model's calendar differs from the stamp (naming the difference
        and ``from_hydrographs(..., model=<this model>)``), or when the model's
        TDIS can no longer date its periods.
    """

    binned = getattr(targets, "_binned", None)
    if binned is None:
        return
    stamp = binned.calendar
    owner = type(targets).__name__
    try:
        on_model = period_calendar(model)
    except ValueError as error:
        raise ValueError(
            f"{caller}: these targets were binned from dated readings against a model's "
            f"calendar, and this model's periods have no dates ({error})."
        ) from None
    difference = clock_difference(stamp, on_model)
    if difference is not None:
        raise ValueError(
            f"{caller}: these targets were binned from dated readings onto another model's "
            f"stress periods -- {difference} -- so a period number here names another date. "
            f"Bin the readings again against this model: {owner}.from_hydrographs(..., "
            "model=<this model>)."
        )


__all__ = [
    "Binning",
    "REPORT_COLUMNS",
    "RULES",
    "bin_readings",
    "clock_difference",
    "model_clock",
    "read_readings",
    "refuse_another_clock",
]
