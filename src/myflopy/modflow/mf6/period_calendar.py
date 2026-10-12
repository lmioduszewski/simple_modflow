"""When each stress period starts and ends: the one reading of a model's clock (ledger 529).

MF6 keeps no dates. Its TDIS holds a ``START_DATE_TIME`` string it stores and never
uses, a ``TIME_UNITS`` word and one length per stress period; a period's dates are
those three put together, and every place that needs one -- PEST's start date, a
period's calendar month, a hydrograph reading binned to a period -- must put them
together the same way. :class:`PeriodCalendar` is that one computation and
:func:`period_calendar` the one reading of a model's TDIS into it;
:func:`tdis_start_date_time` is the one reader of ``START_DATE_TIME``.

**A period's calendar month is the month of its MIDPOINT.** A TDIS that starts one
day before a month -- the Ten Trails model's ``2017-09-30`` -- ends every period at
00:00 on the LAST day of its month, so its START date names the month before;
the midpoint survives that offset, and on calendar-aligned periods
(``mf.AverageYear``) start, midpoint and end name the same month.

**The calendar also carries each period's steady-state flag** -- read from the flow
model's STO by MF6's rule: a GWF model with no STO is steady in every period, and
a STO period block switches the state until the next one. **A period before the
first STO block is TRANSIENT**, not steady (measured on MF6 6.8: a STO declaring
``TRANSIENT`` from period 2 alone ran period 1 with the same storage flows as one
declaring it from period 1; ``docs/reproducers/sto_default_transient.py``).

**How two calendars differ is one rule too** (ledger 546):
:meth:`PeriodCalendar.clock_difference` names the first difference -- the period
count, the start, the unit, a period's length, and (when asked) a period's steady
flag -- and :meth:`PeriodCalendar.same_clock`, a ``from_hydrographs`` set's clock
stamp and ``AverageYear.equilibrium`` all ask it; each says the answer in its own
words. **And the one reading of a run's average-year cycle is here**,
:func:`monthly_cycle` (``mf.AverageYear.from_model`` is it), so the observation layer
reads a cycle without importing ``myflopy.cycles``, which sits above it.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy._time_units import mf6_time_conversion, mf6_time_word
from myflopy.modflow.mf6.hydraulics import _simulation_tdis, declared_time_unit
from myflopy.modflow.mf6.package_explorer_utils import written_period_blocks
from myflopy.modflow.mf6.package_registry import model_packages_of_type
from myflopy.modflow.mf6.period_end import END_RTOL, period_lengths, period_records

logger = get_logger(__name__)

#: Every date is rounded to this. FloPy writes a ``perlen`` at 8 decimals, so an
#: hour read back in days is ``0.04166667`` -- 3600.0003 s -- and a date carrying
#: that 0.3 ms is not the date the modeller declared.
_DATE_RESOLUTION = "ms"

#: Days in each calendar month of a NON-leap year, January = 1. An average year has
#: no 29 February: the user's own 72-period TDIS runs 28-day Februaries throughout.
NONLEAP_DAYS = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}
_MONTH_NAMES = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


def _is_int(value: Any) -> bool:
    return isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))


def _month_name(month: int) -> str:
    return _MONTH_NAMES[int(month) - 1]


def parse_start(value: Any, *, what: str) -> pd.Timestamp:
    """``value`` as a timezone-naive :class:`pandas.Timestamp`, or ``ValueError`` naming ``what``.

    Parameters
    ----------
    value : str or datetime.date or datetime.datetime or pandas.Timestamp or numpy.datetime64
        A date, or a date and time: ``"2017-10-01"``, ``"2017-09-30T12:00:00"``,
        ``datetime.date(2017, 10, 1)``.
    what : str
        What ``value`` is, for the refusal (``"mf.AverageYear(start=)"``).

    Returns
    -------
    pandas.Timestamp
        The instant, timezone-naive.

    Raises
    ------
    ValueError
        When ``value`` is not a date pandas parses -- ``None``, a number, a bool
        and text with no four-digit year such as ``"October"`` (which pandas would
        read as year 1) included -- or carries a time zone (MF6's
        ``START_DATE_TIME`` is a local wall-clock string; a zone would be dropped
        on the way into the TDIS file).
    """

    if value is None or isinstance(value, (bool, np.bool_, int, float, np.number)):
        raise ValueError(f"{what} takes a date such as '2017-10-01', got {value!r}.")
    if isinstance(value, str) and not re.search(r"\d{4}", value):
        # pandas reads "October" as 0001-10-01: a month with no year is not a date.
        raise ValueError(f"{what} takes a date with its year, such as '2017-10-01'; {value!r} is not one.")
    try:
        stamp = pd.Timestamp(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"{what} takes a date such as '2017-10-01'; {value!r} is not one ({error})."
        ) from None
    if stamp is pd.NaT:
        raise ValueError(f"{what} takes a date such as '2017-10-01', got {value!r}.")
    if stamp.tzinfo is not None:
        raise ValueError(
            f"{what} = {value!r} carries a time zone; give the local date and time "
            f"without one ({stamp.tz_localize(None).isoformat()!r}): MF6's START_DATE_TIME "
            "is a wall-clock string."
        )
    return stamp


def format_start(stamp: pd.Timestamp) -> str:
    """A TDIS ``START_DATE_TIME`` string: ``2017-09-30`` at midnight, else ISO date and time.

    Parameters
    ----------
    stamp : pandas.Timestamp
        A timezone-naive instant.

    Returns
    -------
    str
        ``"YYYY-MM-DD"`` when ``stamp`` is a midnight -- the form the Ten Trails TDIS
        carries -- and ``"YYYY-MM-DDTHH:MM:SS"`` (fractional seconds when present)
        otherwise.
    """

    if stamp == stamp.normalize():
        return stamp.strftime("%Y-%m-%d")
    return stamp.isoformat()


@dataclass(frozen=True)
class CalendarDifference:
    """The first way one calendar dates its periods otherwise than another.

    Returned by :meth:`PeriodCalendar.clock_difference`, which decides WHAT differs;
    every caller says it in its own words (a target set's clock stamp, an average
    year's run).

    Attributes
    ----------
    what : str
        What differs, the first of these in this order: ``"nper"`` (the period
        count), ``"start"`` (the instant period 0 begins), ``"time_units"`` (the unit
        of the lengths), ``"perlen"`` (a period's length, beyond FloPy's 8-decimal
        write) or ``"steady"`` (a period's steady-state flag; only when asked).
    mine : int or pandas.Timestamp or str or float or bool
        The value of the calendar ``clock_difference`` was called on: the period
        count, the start, the unit word, period ``per``'s length (in that unit) or
        period ``per``'s flag (``True`` = steady state).
    theirs : int or pandas.Timestamp or str or float or bool
        The other calendar's value of the same thing.
    per : int or None
        The period, zero-based, for ``"perlen"`` (the one whose lengths differ MOST)
        and ``"steady"`` (the first whose flags differ); ``None`` otherwise.
    """

    what: str
    mine: Any
    theirs: Any
    per: int | None = None


@dataclass(frozen=True)
class PeriodCalendar:
    """The dates of every stress period: a start instant, the period lengths and their unit.

    Built from a model by :func:`period_calendar` (``model.calendar``), or directly
    from the three TDIS values -- which is what a reader of a run's FILES does. Every
    date is ``start`` plus the cumulative lengths, converted at MF6's own seconds per
    unit (a ``years`` period is 365.25 days), rounded to the millisecond.

    Parameters
    ----------
    start : str or datetime.date or datetime.datetime or pandas.Timestamp
        The instant period 0 begins: the TDIS ``START_DATE_TIME``. Parsed by
        pandas (``"2017-09-30"``, ``"2017-09-30T12:00:00"``); a value that is not a
        date, or one carrying a time zone, raises ``ValueError``. Stored as a
        timezone-naive :class:`pandas.Timestamp`.
    perlen : sequence of float
        Each stress period's length, zero-based, in ``time_units`` -- the TDIS
        ``PERLEN`` column. Must hold at least one period, each finite and ``>= 0``
        (``ValueError`` otherwise). Stored as a tuple of float.
    time_units : str, default "days"
        The unit of ``perlen``: an MF6 ``TIME_UNITS`` word (``"seconds"``,
        ``"minutes"``, ``"hours"``, ``"days"``, ``"years"``) or any spelling
        :func:`myflopy._time_units.mf6_time_word` takes (``"d"``, ``"day"``), stored
        as MF6's word. ``None``, ``""`` and ``"unknown"`` -- MF6's undefined -- raise
        ``ValueError``: lengths in no unit cannot become dates.
    steady : sequence of bool or None, default None
        Each period's steady-state flag, as the flow model's STO sets it (see the
        module notes for MF6's rule), or ``None`` when not known -- a transport
        model's calendar, or one built from TDIS values alone. A sequence must have
        one entry per period (``ValueError`` otherwise). Stored as a tuple of bool.

    Attributes
    ----------
    nper : int
        The number of stress periods.
    starts, ends, midpoints : pandas.DatetimeIndex
        Each period's start, end and midpoint instants.
    months : numpy.ndarray
        Each period's calendar month, 1-12: its MIDPOINT's (module notes).

    Examples
    --------
    >>> cal = mf.PeriodCalendar("2017-09-30", [31.0, 30.0, 31.0], "days")
    >>> cal.table()[["per", "end", "month"]]
       per        end  month
    0    0 2017-10-31     10
    1    1 2017-11-30     11
    2    2 2017-12-31     12
    """

    start: pd.Timestamp
    perlen: tuple[float, ...]
    time_units: str = "days"
    steady: tuple[bool, ...] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "start", parse_start(self.start, what="PeriodCalendar(start=)"))
        if isinstance(self.perlen, (str, bytes)) or not isinstance(self.perlen, (Sequence, np.ndarray)):
            raise ValueError(
                f"PeriodCalendar(perlen=) takes one length per stress period, got {self.perlen!r}."
            )
        perlen = tuple(float(value) for value in self.perlen)
        if not perlen:
            raise ValueError("PeriodCalendar(perlen=) holds no stress period.")
        bad = [(per, value) for per, value in enumerate(perlen) if not np.isfinite(value) or value < 0]
        if bad:
            per, value = bad[0]
            raise ValueError(
                f"PeriodCalendar(perlen=): period {per} has length {value!r}; every length "
                "must be a finite number >= 0."
            )
        object.__setattr__(self, "perlen", perlen)
        word = mf6_time_word(self.time_units)
        if word is None:
            raise ValueError(
                f"PeriodCalendar(time_units={self.time_units!r}) names no unit, so the period "
                "lengths cannot become dates; declare one, e.g. mf.tdis(time_units='days', ...)."
            )
        object.__setattr__(self, "time_units", word)
        if self.steady is not None:
            steady = tuple(bool(flag) for flag in self.steady)
            if len(steady) != len(perlen):
                raise ValueError(
                    f"PeriodCalendar(steady=) holds {len(steady)} flags for {len(perlen)} periods."
                )
            object.__setattr__(self, "steady", steady)

    @property
    def nper(self) -> int:
        """The number of stress periods."""

        return len(self.perlen)

    def _bounds(self) -> pd.DatetimeIndex:
        """``nper + 1`` instants: period 0's start, then every period's end."""

        seconds = np.concatenate([[0.0], np.cumsum(self.perlen)]) * mf6_time_conversion(self.time_units)
        offsets = pd.to_timedelta(seconds, unit="s").round(_DATE_RESOLUTION)
        return pd.DatetimeIndex(self.start + offsets)

    @property
    def starts(self) -> pd.DatetimeIndex:
        """Each period's start instant, zero-based."""

        return self._bounds()[:-1]

    @property
    def ends(self) -> pd.DatetimeIndex:
        """Each period's end instant, zero-based: the start of the next."""

        return self._bounds()[1:]

    @property
    def midpoints(self) -> pd.DatetimeIndex:
        """Each period's midpoint instant, zero-based, rounded to the millisecond."""

        bounds = self._bounds()
        return pd.DatetimeIndex(bounds[:-1] + (bounds[1:] - bounds[:-1]) / 2).round(_DATE_RESOLUTION)

    @property
    def months(self) -> np.ndarray:
        """Each period's calendar month, 1-12: the month of its MIDPOINT."""

        return np.asarray(self.midpoints.month, dtype=int)

    def table(self) -> pd.DataFrame:
        """One row per stress period: its dates, month, length and steady-state flag.

        Returns
        -------
        pandas.DataFrame
            ``per`` (zero-based), ``start``, ``end``, ``mid`` (timestamps),
            ``month`` (1-12, the midpoint's), ``perlen`` (in :attr:`time_units`) and
            ``steady`` (bool, or ``<NA>`` when the calendar carries no flags).
        """

        frame = pd.DataFrame(
            {
                "per": np.arange(self.nper),
                "start": self.starts,
                "end": self.ends,
                "mid": self.midpoints,
                "month": self.months,
                "perlen": np.asarray(self.perlen, dtype=float),
            }
        )
        frame["steady"] = pd.array(
            [pd.NA] * self.nper if self.steady is None else list(self.steady), dtype="boolean"
        )
        return frame

    def same_clock(self, other: PeriodCalendar) -> bool:
        """Whether ``other`` has the same start, unit and period lengths.

        The lengths agree to :data:`~myflopy.modflow.mf6.period_end.END_RTOL` of each
        period (floored at one time unit): FloPy writes a ``perlen`` at 8 decimals,
        so a calendar read back from a run's files and the one declared in memory
        differ by up to 5e-9 of a period. The steady-state flags are not compared --
        they are the flow model's, not the clock's.

        Parameters
        ----------
        other : PeriodCalendar
            The calendar to compare.

        Returns
        -------
        bool
            ``True`` when both calendars date every period the same.

        Raises
        ------
        TypeError
            When ``other`` is not a :class:`PeriodCalendar`.
        """

        if not isinstance(other, PeriodCalendar):
            raise TypeError(f"same_clock() compares two PeriodCalendars, got {type(other).__name__}.")
        return self.clock_difference(other) is None

    def clock_difference(self, other: PeriodCalendar, *, steady: bool = False) -> CalendarDifference | None:
        """The first way ``other`` dates its periods otherwise than this calendar, or ``None``.

        The ONE comparison of two calendars (ledger 546): :meth:`same_clock` is
        ``clock_difference(other) is None``, a ``from_hydrographs`` target set checks
        its clock stamp through it, and ``AverageYear.equilibrium`` checks a run
        through it with ``steady=True``. Each caller words the answer itself. Checked
        in this order, and the first that differs is returned: the period count, the
        start, the unit, the period lengths -- to
        :data:`~myflopy.modflow.mf6.period_end.END_RTOL` of each period, floored at
        one time unit, as FloPy writes a ``perlen`` at 8 decimals -- and, with
        ``steady=True``, the steady-state flags.

        Parameters
        ----------
        other : PeriodCalendar
            The calendar compared against; its values are the difference's ``theirs``,
            this calendar's its ``mine``.
        steady : bool, default False
            Whether the steady-state flags count. ``False`` compares the CLOCK only --
            what dates a period, which is all a period number needs to name the same
            dates. ``True`` also compares the flags, after every date, when BOTH
            calendars carry them: a transport model's are ``None`` (its periods are as
            steady as the flow model it reads), and then nothing is compared. A
            non-bool raises ``TypeError``.

        Returns
        -------
        CalendarDifference or None
            What differs first (:class:`CalendarDifference`: ``what``, ``mine``,
            ``theirs`` and the period ``per`` for a length or a flag); ``None`` when
            every compared value agrees.

        Raises
        ------
        TypeError
            When ``other`` is not a :class:`PeriodCalendar`, or ``steady`` is not a
            bool.
        """

        if not isinstance(other, PeriodCalendar):
            raise TypeError(
                f"clock_difference() compares two PeriodCalendars, got {type(other).__name__}."
            )
        if not isinstance(steady, (bool, np.bool_)):
            raise TypeError(f"clock_difference(steady=) takes True or False, got {steady!r}.")
        if self.nper != other.nper:
            return CalendarDifference("nper", self.nper, other.nper)
        if self.start != other.start:
            return CalendarDifference("start", self.start, other.start)
        if self.time_units != other.time_units:
            return CalendarDifference("time_units", self.time_units, other.time_units)
        mine = np.asarray(self.perlen)
        theirs = np.asarray(other.perlen)
        gaps = np.abs(mine - theirs)
        if np.any(gaps > END_RTOL * np.maximum(1.0, np.abs(mine))):
            per = int(np.argmax(gaps))
            return CalendarDifference("perlen", self.perlen[per], other.perlen[per], per)
        if steady and self.steady is not None and other.steady is not None and self.steady != other.steady:
            per = next(p for p, (a, b) in enumerate(zip(self.steady, other.steady, strict=True)) if a != b)
            return CalendarDifference("steady", self.steady[per], other.steady[per], per)
        return None


def tdis_start_date_time(model: Any) -> str | None:
    """The TDIS ``START_DATE_TIME`` of the simulation ``model`` belongs to, or ``None``.

    The ONE reader of it: :func:`period_calendar` dates the periods from it and
    ``model.pest()`` takes PEST's ``start_datetime`` from it when none is given.

    Parameters
    ----------
    model : object
        A myflopy model (its ``sim`` or ``simulation``) or anything whose ``gwf`` is a
        FloPy model (its ``gwf.simulation``); the first of those three that has a
        ``tdis`` is read.

    Returns
    -------
    str or None
        The string as the TDIS holds it (``"2017-09-30"``); ``None`` when there is
        no TDIS or it declares no start, an empty string included.
    """

    tdis = _simulation_tdis(model)
    held = getattr(tdis, "start_date_time", None) if tdis is not None else None
    if held is None:
        return None
    getter = getattr(held, "get_data", None)
    value = getter() if callable(getter) else held
    if value is None:
        return None
    return str(value).strip() or None


def steady_periods(model: Any, nper: int) -> tuple[bool, ...] | None:
    """Each stress period's steady-state flag by MF6's STO rule, or ``None`` for a non-flow model.

    A GWF model with no STO package is steady in every period. With one, a period
    block's ``STEADY-STATE`` or ``TRANSIENT`` holds until the next block, and a
    period BEFORE the first block is TRANSIENT (module notes; measured). Read through
    :func:`~myflopy.modflow.mf6.package_explorer_utils.written_period_blocks`, so
    asking never adds a block (ledger 449).

    Parameters
    ----------
    model : SimulationBase or object
        A myflopy model; its ``model_type`` (``"gwf6"`` when absent) decides whether it
        has storage flags of its own, and its STO is found by TYPE
        (:func:`~myflopy.modflow.mf6.package_registry.model_packages_of_type`).
    nper : int
        The number of stress periods.

    Returns
    -------
    tuple of bool or None
        ``True`` for a steady-state period, zero-based; ``None`` for a GWT, GWE or PRT
        model, whose periods are as steady as the flow model they read.

    Raises
    ------
    ValueError
        When the model carries more than one STO package (MF6 allows one).
    """

    if str(getattr(model, "model_type", "gwf6") or "gwf6").lower() != "gwf6":
        return None
    packages = model_packages_of_type(model, "sto")
    if not packages:
        return (True,) * int(nper)
    if len(packages) > 1:
        raise ValueError(
            f"The model carries {len(packages)} STO packages "
            f"({', '.join(str(p.package_name) for p in packages)}); MF6 reads one."
        )
    sto = packages[0]
    steady_blocks = written_period_blocks(getattr(sto, "steady_state", None))
    transient_blocks = written_period_blocks(getattr(sto, "transient", None))
    state = False
    flags = []
    for per in range(int(nper)):
        if steady_blocks.get(per):
            state = True
        elif transient_blocks.get(per):
            state = False
        flags.append(state)
    return tuple(flags)


def period_calendar(model: Any) -> PeriodCalendar:
    """The :class:`PeriodCalendar` of ``model``'s simulation: its TDIS dates and STO flags.

    ``model.calendar`` is this. The TDIS of a live model and of one read back from
    its files (``mf.load_mf6_run``) give the same calendar
    (:meth:`PeriodCalendar.same_clock`).

    Parameters
    ----------
    model : SimulationBase or object
        A myflopy model whose simulation carries a TDIS (``model.sim.tdis``).

    Returns
    -------
    PeriodCalendar
        Its start from the TDIS ``START_DATE_TIME``, its lengths and unit from the
        period block and ``TIME_UNITS``, and its steady flags from the model's STO
        (:func:`steady_periods`; ``None`` for a non-flow model).

    Raises
    ------
    ValueError
        When the TDIS declares no ``START_DATE_TIME`` or a ``TIME_UNITS`` MF6 reads as
        undefined -- naming the ``mf.tdis(start_date_time=..., time_units='days')``
        that declares them -- when it has no period block, or when the start is not
        a date (:func:`parse_start`).

    Examples
    --------
    >>> model.calendar.table().head()
    >>> mf.period_calendar(mf.load_mf6_run("runs/baseline")).ends[-1]
    """

    start = tdis_start_date_time(model)
    if start is None:
        raise ValueError(
            "The model's TDIS declares no START_DATE_TIME, so its periods have no dates; "
            "declare one: mf.tdis(start_date_time='2017-10-01', ...) (or "
            "mf.AverageYear(...).tdis(), which writes it)."
        )
    perlen = period_lengths(model)
    if perlen is None:
        raise ValueError("The model carries no TDIS period block, so its periods have no lengths.")
    units = declared_time_unit(model)
    if units is None:
        raise ValueError(
            "The model's TDIS declares no TIME_UNITS MF6 reads, so its period lengths "
            "cannot become dates; declare one: mf.tdis(time_units='days', ...)."
        )
    stamp = parse_start(start, what="The TDIS START_DATE_TIME")
    return PeriodCalendar(stamp, perlen, units, steady_periods(model, len(perlen)))


@dataclass(frozen=True)
class MonthlyCycle:
    """The period layout of an average water year: an optional steady period, then cycles of twelve months.

    What ``mf.AverageYear`` declares and :func:`monthly_cycle` reads back from a run,
    in the one place the import layers let both reach: ``AverageYear`` keeps one and
    answers ``months_of`` and its midpoint check through it, and the observation
    layer (``from_hydrographs(mode="average_year")``) reads a run's last cycle
    through :func:`monthly_cycle` without importing ``myflopy.cycles``. Nothing here
    is validated: ``AverageYear`` checks what it is given, and :func:`monthly_cycle`
    builds one only from a run that IS such a cycle.

    Parameters
    ----------
    cycles : int
        How many times the twelve months repeat after the steady period, ``>= 1``.
    start : pandas.Timestamp
        The first day of the first cycle, at midnight; its month is the cycle's
        first month.
    steady : bool
        Whether a steady-state period, period 0, comes first.
    steady_length : float or None
        The steady period's length in days; ``None`` without one.
    nstp : tuple of int
        Twelve time-step counts by CALENDAR month, January first.
    tsmult : float
        The time-step multiplier of every monthly period.

    Attributes
    ----------
    months : tuple of int
        The twelve calendar months in cycle order, ``start``'s month first.
    nper : int
        The number of stress periods, ``int(steady) + 12 * cycles``.
    """

    cycles: int
    start: pd.Timestamp
    steady: bool
    steady_length: float | None
    nstp: tuple[int, ...]
    tsmult: float

    @property
    def months(self) -> tuple[int, ...]:
        """The twelve calendar months in cycle order, ``start``'s month first."""

        return tuple((self.start.month - 1 + k) % 12 + 1 for k in range(12))

    @property
    def nper(self) -> int:
        """The number of stress periods, ``int(steady) + 12 * cycles``."""

        return int(self.steady) + 12 * self.cycles

    def period_months(self) -> list[int | None]:
        """Each stress period's calendar month, zero-based; ``None`` for the steady period.

        Returns
        -------
        list of int or None
            ``nper`` entries: ``None`` first when there is a steady period, then
            :attr:`months` once per cycle.
        """

        return ([None] if self.steady else []) + list(self.months) * self.cycles

    def months_of(self, cycle: int = -1) -> dict[int, int]:
        """``{calendar month: period}`` for one cycle, in cycle order.

        Parameters
        ----------
        cycle : int, default -1
            The cycle, zero-based; negative counts from the last (``-1`` is the last
            cycle, the one targets are compared on). Out of range raises
            ``IndexError``; a bool or non-integer ``TypeError``.

        Returns
        -------
        dict of int to int
            Each calendar month (January = 1) to its zero-based stress period.
        """

        if not _is_int(cycle):
            raise TypeError(f"months_of(cycle=) takes an int, got {cycle!r}.")
        index = int(cycle) + self.cycles if int(cycle) < 0 else int(cycle)
        if not 0 <= index < self.cycles:
            raise IndexError(f"months_of(cycle={cycle!r}): this year runs {self.cycles} cycle(s).")
        first = int(self.steady) + 12 * index
        return {month: first + offset for offset, month in enumerate(self.months)}

    def check_midpoints(self, calendar: PeriodCalendar) -> None:
        """Refuse a cycle whose leap-day drift moves a period's midpoint out of its month.

        The months run NON-leap lengths, so each leap year starts the later periods
        one day early; after about 14 leap years (some 56 cycles) a period's MIDPOINT
        -- which names its month on a :class:`PeriodCalendar` -- leaves the month
        the cycle gives it.

        Parameters
        ----------
        calendar : PeriodCalendar
            The calendar of these periods (``nper`` of them): ``AverageYear.calendar``
            for a declared year, the run's own for one read back.

        Raises
        ------
        ValueError
            Naming the first period whose midpoint month is not its cycle month.
        """

        first = int(self.steady)
        assigned = np.asarray(self.period_months()[first:])
        found = calendar.months[first:]
        off = np.flatnonzero(assigned != found)
        if off.size:
            per = int(off[0]) + first
            raise ValueError(
                f"AverageYear: period {per} is {_month_name(assigned[off[0]])} in the cycle but "
                f"its midpoint falls in {_month_name(found[off[0]])} -- the 28-day Februaries "
                "have drifted the periods past half a month over this many leap years; run "
                "fewer cycles, or start in a later year."
            )


def monthly_cycle(model: Any) -> MonthlyCycle:
    """The average-year cycle a run was built on, read back from its own TDIS and STO.

    The ONE reading of a run's cycle: ``mf.AverageYear.from_model`` is this, and
    ``from_hydrographs(mode="average_year")`` places its targets on
    ``monthly_cycle(model).months_of(-1)``. ``mf.load_mf6_run(ws)`` and the PEST build
    read files only (ledger 356), so the cycle is inferred, not recorded. The run must
    be days, with ``START_DATE_TIME`` declared; its periods 0 or 1 steady-state period,
    then transient periods of non-leap calendar month lengths starting on the first of
    a month, in whole twelve-month cycles; each month's step count the same in every
    cycle and one multiplier throughout. The steady period's own steps and multiplier
    are not read (``AverageYear`` always writes 1 and 1.0).

    Parameters
    ----------
    model : SimulationBase
        A myflopy GWF model -- ``run.model(name)`` or ``mf.load_mf6_run(ws)``.

    Returns
    -------
    MonthlyCycle
        The layout ``mf.AverageYear`` would declare for this run.

    Raises
    ------
    ValueError
        When the run is not such a cycle, naming the first period that breaks the
        pattern and its length; when its TDIS has no ``START_DATE_TIME``, no
        ``TIME_UNITS`` or units other than days; when it is not a flow model (the STO
        flags are the flow model's); or when a period's midpoint has drifted out of its
        month (:meth:`MonthlyCycle.check_midpoints`).
    """

    calendar = period_calendar(model)
    if calendar.time_units != "days":
        raise ValueError(
            f"The run's TDIS counts time in {calendar.time_units}; an AverageYear is in days."
        )
    if calendar.steady is None:
        raise ValueError(
            "AverageYear.from_model reads the steady-state flags of a FLOW model's STO; "
            "pass the GWF model."
        )
    steady = bool(calendar.steady[0])
    first = int(steady)
    late = [per for per in range(first, calendar.nper) if calendar.steady[per]]
    if late:
        raise ValueError(
            f"Period {late[0]} is steady state; an AverageYear has at most one steady "
            "period, period 0, and every period after it transient."
        )
    count = calendar.nper - first
    if count < 12 or count % 12:
        raise ValueError(
            f"The run has {count} transient period(s), not a whole number of twelve-month "
            "cycles."
        )
    begin = calendar.starts[first]
    if begin != begin.normalize().replace(day=1):
        raise ValueError(
            f"Period {first}, the first month, starts at {begin.isoformat()}, not at midnight "
            "on the first of a month."
        )
    records = period_records(model) or []
    months = tuple((begin.month - 1 + k) % 12 + 1 for k in range(12))
    steps: dict[int, int] = {}
    multipliers = set()
    for offset in range(count):
        per = first + offset
        month = months[offset % 12]
        length = calendar.perlen[per]
        expected = NONLEAP_DAYS[month]
        if abs(length - expected) > END_RTOL * max(1.0, expected):
            raise ValueError(
                f"Period {per} is {length:g} days long where an AverageYear's "
                f"{_month_name(month)} is {expected}; the run is not a calendar cycle "
                f"from {format_start(begin)}."
            )
        _, nstp, tsmult = records[per] if per < len(records) else (None, None, None)
        if nstp is None or tsmult is None:
            raise ValueError(f"The run's TDIS holds no nstp/tsmult for period {per}.")
        if steps.setdefault(month, nstp) != nstp:
            raise ValueError(
                f"Period {per} ({_month_name(month)}) takes {nstp} steps where an earlier "
                f"{_month_name(month)} took {steps[month]}; an AverageYear repeats each "
                "month's steps."
            )
        multipliers.add(round(tsmult, 6))
    if len(multipliers) > 1:
        raise ValueError(
            f"The run's monthly periods use several time-step multipliers ({sorted(multipliers)}); "
            "an AverageYear uses one."
        )
    cycle = MonthlyCycle(
        cycles=count // 12,
        start=begin,
        steady=steady,
        steady_length=float(calendar.perlen[0]) if steady else None,
        nstp=tuple(int(steps[month]) for month in range(1, 13)),
        tsmult=float(records[first][2]),
    )
    cycle.check_midpoints(calendar)
    return cycle


__all__ = [
    "NONLEAP_DAYS",
    "CalendarDifference",
    "MonthlyCycle",
    "PeriodCalendar",
    "format_start",
    "monthly_cycle",
    "parse_start",
    "period_calendar",
    "steady_periods",
    "tdis_start_date_time",
]
