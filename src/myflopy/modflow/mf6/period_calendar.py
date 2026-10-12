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
from myflopy.modflow.mf6.period_end import END_RTOL, period_lengths

logger = get_logger(__name__)

#: Every date is rounded to this. FloPy writes a ``perlen`` at 8 decimals, so an
#: hour read back in days is ``0.04166667`` -- 3600.0003 s -- and a date carrying
#: that 0.3 ms is not the date the modeller declared.
_DATE_RESOLUTION = "ms"


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
        if self.start != other.start or self.time_units != other.time_units or self.nper != other.nper:
            return False
        mine = np.asarray(self.perlen)
        theirs = np.asarray(other.perlen)
        return bool(np.all(np.abs(mine - theirs) <= END_RTOL * np.maximum(1.0, np.abs(mine))))


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


__all__ = [
    "PeriodCalendar",
    "format_start",
    "parse_start",
    "period_calendar",
    "steady_periods",
    "tdis_start_date_time",
]
