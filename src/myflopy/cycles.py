"""An average water year: one steady period, then repeats of the twelve calendar months.

A model calibrated to an AVERAGE year runs a leading steady-state period under the
year's mean forcing, then the same twelve months ``cycles`` times until the heads
repeat from one year to the next; targets are compared on the LAST year.
:class:`AverageYear` states that period-to-month map ONCE, and every period-keyed
input reads it from there (ledger 530):

* the TDIS -- ``year.tdis()`` -- one period per calendar month, its length the month's
  days in a NON-leap year, starting on the first of ``start``'s month;
* the STO flags -- ``mf.sto(steady_state=year.steady_state, transient=year.transient)``;
* every per-period forcing -- ``year.by_month({1: jan, ..., 12: dec})`` gives each period
  its month's value, and the steady period the month-length-weighted annual mean
  (the user's rule, Q3);
* where the run stands -- ``year.equilibrium(model)`` reports how far the last two
  cycles differ, period end against period end.

A run READ BACK from its files recovers its cycle from its own TDIS and STO
(:meth:`AverageYear.from_model`); nothing is recorded beside the run, because the
PEST build and ``mf.load_mf6_run`` read the files and nothing else (ledger 356).

Dates go through the one calendar, :class:`~myflopy.modflow.mf6.period_calendar.PeriodCalendar`:
a period's month is its MIDPOINT's there, and the cycle a year declares is checked
to agree with it.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from numbers import Real
from typing import Any

import numpy as np
import pandas as pd

from myflopy import viz
from myflopy._logging import get_logger
from myflopy.modflow.mf6.headsplus import _is_sentinel
from myflopy.modflow.mf6.observations._shared import _stage_at_period_end
from myflopy.modflow.mf6.package_budget import build_lak_stage_result_table
from myflopy.modflow.mf6.package_plotting import normalize_backend
from myflopy.modflow.mf6.package_registry import model_packages_of_type
from myflopy.modflow.mf6.period_calendar import (
    PeriodCalendar,
    format_start,
    parse_start,
    period_calendar,
)
from myflopy.modflow.mf6.period_end import END_RTOL, period_ends, period_records
from myflopy.package_api import tdis as _tdis

logger = get_logger(__name__)

#: Days in each calendar month of a NON-leap year, January = 1. An average year has
#: no 29 February: the user's own 72-period TDIS runs 28-day Februaries throughout.
NONLEAP_DAYS = {1: 31, 2: 28, 3: 31, 4: 30, 5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31, 11: 30, 12: 31}
_MONTH_NAMES = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

#: What ``AverageYear.tdis()`` sets itself, and the argument that sets it.
_TDIS_OWNED = {
    "nper": "cycles= and steady=",
    "perioddata": "the calendar month lengths with nstp= and tsmult=",
    "time_units": "AverageYear itself (always 'days': a month length is a number of days)",
    "start_date_time": "start= and steady_length=",
}


def _is_int(value: Any) -> bool:
    return isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))


def _is_number(value: Any) -> bool:
    return isinstance(value, (Real, np.number)) and not isinstance(value, (bool, np.bool_))


def _month_name(month: int) -> str:
    return _MONTH_NAMES[int(month) - 1]


class AverageYear:
    """One steady-state period, then ``cycles`` repeats of the twelve calendar months.

    The one statement of a cyclic, average-year time setup (module notes): build the
    TDIS with :meth:`tdis`, the STO flags from :attr:`steady_state` /
    :attr:`transient`, every per-period forcing with :meth:`by_month`, and check the
    run with :meth:`equilibrium`. ``nper`` is ``int(steady) + 12 * cycles``.

    Parameters
    ----------
    cycles : int
        How many times the twelve months repeat after the steady period; ``>= 1``
        (``ValueError``; a bool or non-integer raises ``TypeError``). The run time
        is roughly proportional to it, and :meth:`equilibrium` -- which needs two --
        says whether the last two agree.
    start : str or datetime.date or datetime.datetime or pandas.Timestamp
        The first day of the first cycle, at midnight: ``"2016-10-01"``. Its month
        is the cycle's first month -- October gives a water year. Anything else
        (another day of the month, a time of day) raises ``ValueError`` naming the
        first of that month, because the periods end on month ends only from a
        month's first day (until the cycles pass a 29 February: each leap year the
        28-day Februaries end the later periods one day early; see Raises). Parsed as :func:`~myflopy.modflow.mf6.period_calendar.parse_start`
        parses it (a time zone raises).
    steady : bool, default True
        Whether a steady-state period comes first. ``True`` makes period 0 steady
        state, run under :meth:`by_month`'s annual-mean forcing, and every cycle
        transient after it; ``False`` makes every period transient (the run then
        starts from the initial heads). A non-bool raises ``TypeError``.
    steady_length : float, optional
        The steady period's length in days, ``> 0`` and finite (``ValueError``).
        Default ``None``: 1.0 day when ``steady`` is ``True``. It only sets the dates
        -- a steady period's length changes no head -- and the TDIS
        ``START_DATE_TIME`` is ``start`` minus it, so the first month begins on
        ``start``. Given with ``steady=False`` it raises ``ValueError``: there is
        no steady period for it to set.
    nstp : int or mapping of int to int
        Time steps per month, ``>= 1``: one int for every month, or
        ``{calendar month: nstp}`` naming all twelve months (January = 1; a missing
        month or a key outside 1-12 raises ``ValueError``). REQUIRED -- it is the
        run-time lever and changes the answer (a 72-period run took 3 h 53 min at
        10 steps), so no default suits both accuracy and a calibration's run count.
        The steady period always takes 1 step with multiplier 1.0.
    tsmult : float, default 1.0
        The time-step multiplier of every monthly period, ``> 0`` and finite
        (``ValueError``).

    Attributes
    ----------
    nper : int
        The number of stress periods, ``int(steady) + 12 * cycles``.
    months : tuple of int
        The twelve calendar months in cycle order, ``start``'s month first.
    periods : pandas.DataFrame
        One row per stress period: ``per``, ``cycle`` (zero-based, ``<NA>`` for the
        steady period), ``month`` (1-12, ``<NA>`` for it), ``start``, ``end``,
        ``perlen`` (days), ``nstp``, ``tsmult`` and ``steady``.
    perioddata : list of tuple
        ``(perlen, nstp, tsmult)`` per period, as :meth:`tdis` writes it.
    start_date_time : str
        The TDIS ``START_DATE_TIME``: ``start`` minus the steady period.
    calendar : PeriodCalendar
        The dates of every period, as ``model.calendar`` reads them back.
    steady_state, transient : dict or None
        The ``mf.sto`` flags (see each property).

    Raises
    ------
    ValueError
        When a period's MIDPOINT would fall in another month than the one the cycle
        assigns it: the calendar runs 28-day Februaries, so each leap year the
        periods start one day earlier than the true month, which crosses a midpoint
        only after about 14 leap years (some 56 cycles).

    Examples
    --------
    >>> year = mf.AverageYear(3, start="2016-10-01", nstp=10, tsmult=1.1)
    >>> sim = mf.simulation(
    ...     mf.gwf("flow", packages=[
    ...         ..., mf.sto(sy=0.15, steady_state=year.steady_state, transient=year.transient),
    ...         mf.rch.array(recharge=year.by_month(monthly_recharge)),
    ...     ]),
    ...     packages=[year.tdis(ats=True), mf.ims(complexity="complex")],
    ... )
    >>> year.equilibrium(run.model("flow"), tolerance=0.01).summary()
    """

    def __init__(
        self,
        cycles: int,
        *,
        start: Any,
        steady: bool = True,
        steady_length: float | None = None,
        nstp: int | Mapping[int, int],
        tsmult: float = 1.0,
    ) -> None:
        if not _is_int(cycles):
            raise TypeError(f"AverageYear(cycles=) takes a whole number of years, got {cycles!r}.")
        if int(cycles) < 1:
            raise ValueError(f"AverageYear(cycles={cycles!r}): at least one cycle is needed.")
        stamp = parse_start(start, what="AverageYear(start=)")
        first = stamp.normalize().replace(day=1)
        if stamp != first:
            raise ValueError(
                f"AverageYear(start={start!r}) is not the first day of a month at midnight; "
                f"a cycle starts on one so every period ends on a month end: "
                f"start={format_start(first)!r}."
            )
        if not isinstance(steady, (bool, np.bool_)):
            raise TypeError(f"AverageYear(steady=) takes True or False, got {steady!r}.")
        steady = bool(steady)
        if steady_length is not None and not steady:
            raise ValueError(
                f"AverageYear(steady_length={steady_length!r}) with steady=False: there is no "
                "steady period for it to set; drop it, or pass steady=True."
            )
        if steady:
            length = 1.0 if steady_length is None else steady_length
            if not _is_number(length) or not np.isfinite(length) or float(length) <= 0:
                raise ValueError(
                    f"AverageYear(steady_length={steady_length!r}) must be a finite number of days > 0."
                )
            steady_length = float(length)
        self._cycles = int(cycles)
        self._start = stamp
        self._steady = steady
        self._steady_length = steady_length
        self._nstp = self._check_nstp(nstp)
        if not _is_number(tsmult) or not np.isfinite(tsmult) or float(tsmult) <= 0:
            raise ValueError(f"AverageYear(tsmult={tsmult!r}) must be a finite number > 0.")
        self._tsmult = float(tsmult)
        self._months = tuple((stamp.month - 1 + k) % 12 + 1 for k in range(12))
        self._check_midpoints()

    # -- construction -----------------------------------------------------------

    @staticmethod
    def _check_nstp(nstp: Any) -> tuple[int, ...]:
        """``nstp`` as twelve step counts by calendar month (index 0 = January)."""

        if _is_int(nstp):
            if int(nstp) < 1:
                raise ValueError(f"AverageYear(nstp={nstp!r}) must be >= 1.")
            return (int(nstp),) * 12
        if isinstance(nstp, Mapping):
            stray = [key for key in nstp if not _is_int(key) or not 1 <= int(key) <= 12]
            if stray:
                raise ValueError(
                    f"AverageYear(nstp=) keys are calendar months 1-12 (January = 1); "
                    f"{stray!r} is not one."
                )
            missing = [month for month in range(1, 13) if month not in {int(key) for key in nstp}]
            if missing:
                raise ValueError(
                    "AverageYear(nstp=) names every calendar month, or is one int for all; "
                    f"missing: {', '.join(_month_name(m) for m in missing)}."
                )
            counts = {int(key): value for key, value in nstp.items()}
            bad = [m for m, value in counts.items() if not _is_int(value) or int(value) < 1]
            if bad:
                raise ValueError(
                    f"AverageYear(nstp=): {_month_name(bad[0])} has {counts[bad[0]]!r}; "
                    "every step count must be an int >= 1."
                )
            return tuple(int(counts[month]) for month in range(1, 13))
        raise TypeError(
            f"AverageYear(nstp=) takes an int or {{calendar month: int}}, got {nstp!r}."
        )

    def _check_midpoints(self) -> None:
        """Refuse a cycle whose leap-day drift moves a midpoint out of its month."""

        calendar = self.calendar
        first = int(self._steady)
        assigned = np.asarray(self._period_months()[first:])
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

    # -- what it declares -------------------------------------------------------

    @property
    def cycles(self) -> int:
        """How many times the twelve months repeat."""

        return self._cycles

    @property
    def start(self) -> pd.Timestamp:
        """The first day of the first cycle."""

        return self._start

    @property
    def steady(self) -> bool:
        """Whether a steady-state period comes first."""

        return self._steady

    @property
    def steady_length(self) -> float | None:
        """The steady period's length in days, ``None`` without one."""

        return self._steady_length

    @property
    def nstp(self) -> int | dict[int, int]:
        """Time steps per month: one int, or ``{calendar month: nstp}`` when they differ."""

        if len(set(self._nstp)) == 1:
            return self._nstp[0]
        return {month: self._nstp[month - 1] for month in range(1, 13)}

    @property
    def tsmult(self) -> float:
        """The time-step multiplier of every monthly period."""

        return self._tsmult

    @property
    def months(self) -> tuple[int, ...]:
        """The twelve calendar months in cycle order, ``start``'s month first."""

        return self._months

    @property
    def nper(self) -> int:
        """The number of stress periods, ``int(steady) + 12 * cycles``."""

        return int(self._steady) + 12 * self._cycles

    def _period_months(self) -> list[int | None]:
        """Each period's calendar month, ``None`` for the steady period."""

        return ([None] if self._steady else []) + list(self._months) * self._cycles

    @property
    def perioddata(self) -> list[tuple[float, int, float]]:
        """``(perlen, nstp, tsmult)`` per stress period, as :meth:`tdis` writes them."""

        rows = [(self._steady_length, 1, 1.0)] if self._steady else []
        for month in self._period_months()[int(self._steady):]:
            rows.append((float(NONLEAP_DAYS[month]), self._nstp[month - 1], self._tsmult))
        return rows

    @property
    def start_date_time(self) -> str:
        """The TDIS ``START_DATE_TIME``: ``start`` minus the steady period's length."""

        begin = self._start
        if self._steady:
            begin = begin - pd.to_timedelta(self._steady_length, unit="D")
        return format_start(begin)

    @property
    def calendar(self) -> PeriodCalendar:
        """The dates of every period -- what ``model.calendar`` reads back from the run."""

        return PeriodCalendar(
            self.start_date_time,
            [row[0] for row in self.perioddata],
            "days",
            steady=[self._steady] + [False] * (12 * self._cycles) if self._steady else [False] * self.nper,
        )

    @property
    def periods(self) -> pd.DataFrame:
        """One row per stress period: ``per``, ``cycle``, ``month``, dates, length, steps, flag.

        ``cycle`` and ``month`` are ``<NA>`` for the steady period; ``start``/``end``
        come from :attr:`calendar`.
        """

        table = self.calendar.table()
        first = int(self._steady)
        cycle = [pd.NA] * first + [index // 12 for index in range(12 * self._cycles)]
        rows = self.perioddata
        return pd.DataFrame(
            {
                "per": table["per"],
                "cycle": pd.array(cycle, dtype="Int64"),
                "month": pd.array(self._period_months(), dtype="Int64"),
                "start": table["start"],
                "end": table["end"],
                "perlen": [row[0] for row in rows],
                "nstp": [row[1] for row in rows],
                "tsmult": [row[2] for row in rows],
                "steady": table["steady"].astype(bool),
            }
        )

    @property
    def steady_state(self) -> dict[int, bool] | None:
        """``mf.sto(steady_state=)``: ``{0: True}``, or ``None`` with no steady period.

        ``None`` rather than ``{}``: FloPy writes an EMPTY dict as empty period
        blocks, which drops the transient flags MF6 needed (CLAUDE.md, USG import).
        """

        return {0: True} if self._steady else None

    @property
    def transient(self) -> dict[int, bool]:
        """``mf.sto(transient=)``: ``{1: True}`` after a steady period, else ``{0: True}``.

        One block: MF6 keeps a period's STO state until the next block, so every
        monthly period after it is transient.
        """

        return {int(self._steady): True}

    def tdis(self, **options: Any) -> Any:
        """The ``mf.tdis`` spec of this cycle: its periods, ``days`` and its start date.

        Parameters
        ----------
        **options
            Every other ``mf.tdis`` argument -- ``ats=True`` and the ``ats_*``
            settings, ``filename``, ``pname``, ``name``, extra FloPy options -- is
            passed through. ``nper``, ``perioddata``, ``time_units`` and
            ``start_date_time`` are this object's and raise ``ValueError`` naming
            the argument that sets them.

        Returns
        -------
        PackageSpec
            ``mf.tdis(nper=self.nper, perioddata=self.perioddata, time_units="days",
            start_date_time=self.start_date_time, **options)``.
        """

        owned = sorted(set(options) & set(_TDIS_OWNED))
        if owned:
            raise ValueError(
                "AverageYear.tdis() sets "
                + "; ".join(f"{key} (from {_TDIS_OWNED[key]})" for key in owned)
                + " itself; change the AverageYear instead."
            )
        return _tdis(
            nper=self.nper,
            perioddata=self.perioddata,
            time_units="days",
            start_date_time=self.start_date_time,
            **options,
        )

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
        index = int(cycle) + self._cycles if int(cycle) < 0 else int(cycle)
        if not 0 <= index < self._cycles:
            raise IndexError(f"months_of(cycle={cycle!r}): this year runs {self._cycles} cycle(s).")
        first = int(self._steady) + 12 * index
        return {month: first + offset for offset, month in enumerate(self._months)}

    # -- forcing ----------------------------------------------------------------

    def by_month(self, values: Mapping[int, Any] | pd.Series, *, steady: Any = None) -> dict[int, Any]:
        """``{period: value}`` for every period: each month's value, and the steady period's.

        Pass the result wherever a helper takes a value per period:
        ``mf.rch.array(recharge=year.by_month(...))``, ``mf.lak(..., rainfall=year.by_month(...))``,
        ``mf.chd(stress_period_data=year.by_month(...))``. Every period gets a block of
        its own, so nothing rides on MF6's rule of repeating an earlier one.

        Parameters
        ----------
        values : mapping or pandas.Series
            Keyed by CALENDAR month, January = 1, whatever month ``start`` is; all
            twelve must be named (a missing month raises ``ValueError`` naming it,
            a key outside 1-12 or a non-integer key too). A list, tuple or array
            raises ``TypeError``: twelve values in a row are ambiguous between
            calendar order and cycle order when the cycle starts in October. Each
            value is anything the receiving slot takes -- a number, an ``(ncpl,)``
            array, a ``{lake_id: value}`` dict for ``mf.lak``, a records list --
            and the SAME object is placed in every period of its month (an array is
            not copied). A ``None`` value raises ``ValueError``.
        steady : object, optional
            The steady period's value. Default ``None``: the MONTH-LENGTH-WEIGHTED
            ANNUAL MEAN of the twelve values (days of a non-leap year, /365) when
            every value is numeric -- all numbers, all arrays of one shape, or all
            dicts of numbers with the same keys, averaged key by key -- the steady
            state of the mean forcing, from which a periodic run settles soonest (the
            user's rule, Q3). Any other values (records, time-series names, recipes,
            mixed kinds) need it given (``ValueError``). Given when the year has no
            steady period, it raises ``ValueError``.

        Returns
        -------
        dict of int to object
            ``{period: value}`` for every zero-based period, the steady period first.
        """

        if isinstance(values, pd.Series):
            items = list(values.items())
        elif isinstance(values, Mapping):
            items = list(values.items())
        else:
            raise TypeError(
                f"by_month() takes {{calendar month: value}} (January = 1), got "
                f"{type(values).__name__}: a sequence of twelve is ambiguous between calendar "
                "order and cycle order. Key it by month: dict(zip(range(1, 13), values))."
            )
        stray = [key for key, _ in items if not _is_int(key) or not 1 <= int(key) <= 12]
        if stray:
            raise ValueError(
                f"by_month() keys are calendar months 1-12 (January = 1); {stray[:5]!r} is not."
            )
        monthly = {int(key): value for key, value in items}
        if len(monthly) != len(items):
            raise ValueError("by_month() names a calendar month twice.")
        missing = [month for month in range(1, 13) if month not in monthly]
        if missing:
            raise ValueError(
                "by_month() needs a value for every calendar month; missing: "
                + ", ".join(_month_name(m) for m in missing) + "."
            )
        empty = [month for month, value in monthly.items() if value is None]
        if empty:
            raise ValueError(
                f"by_month(): {', '.join(_month_name(m) for m in empty)} hold None; give every "
                "month a value (MF6 would repeat the previous block there)."
            )
        result: dict[int, Any] = {}
        if self._steady:
            result[0] = _annual_mean(monthly) if steady is None else steady
        elif steady is not None:
            raise ValueError(
                "by_month(steady=) given, but this AverageYear has no steady period (steady=False)."
            )
        for per, month in enumerate(self._period_months()[int(self._steady):], start=int(self._steady)):
            result[per] = monthly[month]
        return result

    # -- reading a run ----------------------------------------------------------

    @classmethod
    def from_model(cls, model: Any) -> AverageYear:
        """The cycle a run was built on, read back from its own TDIS and STO.

        How a run recovers its year: ``mf.load_mf6_run(ws)`` and the PEST build read
        files only (ledger 356), so the cycle is inferred, not recorded. The run must
        be days, with ``START_DATE_TIME`` declared; its periods 0 or 1 steady-state
        period, then transient periods of non-leap calendar month lengths starting on
        the first of a month, in whole twelve-month cycles; each month's step count
        the same in every cycle and one multiplier throughout. The steady period's own
        steps and multiplier are not read (this object always writes 1 and 1.0).

        Parameters
        ----------
        model : SimulationBase
            A myflopy GWF model -- ``run.model(name)`` or ``mf.load_mf6_run(ws)``.

        Returns
        -------
        AverageYear
            Equal (``==``) to the one that declared the run.

        Raises
        ------
        ValueError
            When the run is not such a cycle, naming the first period that breaks the
            pattern and its length; when its TDIS has no ``START_DATE_TIME``, no
            ``TIME_UNITS`` or units other than days; or when it is not a flow model
            (the STO flags are the flow model's).
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
        return cls(
            count // 12,
            start=begin,
            steady=steady,
            steady_length=calendar.perlen[0] if steady else None,
            nstp={month: steps[month] for month in range(1, 13)},
            tsmult=float(records[first][2]),
        )

    def equilibrium(
        self, model: Any, *, targets: Any = None, tolerance: float | None = None
    ) -> CycleEquilibrium:
        """How far the last two cycles of a run differ, period end against period end.

        Parameters
        ----------
        model : SimulationBase
            The RUN of this year: ``run.model(name)`` or ``mf.load_mf6_run(ws)`` --
            the flow model, or, for ``ConcTargets``, the transport model holding
            their concentration file. Its calendar must be this year's: same
            number of periods, start, unit and lengths, and, for a flow model, the
            same STO steady flags (``ValueError`` naming the first difference).
        targets : target set, sequence of target sets, or None, optional
            Where to compare. Default ``None``: EVERYWHERE -- every cell's head
            (any discretization; a flow model only, ``ValueError`` otherwise) and
            every lake's stage (each LAK package saving a stage file). A
            ``HeadTargets``/``ConcTargets`` (``simulated_heads``) or a
            ``LakeStageTargets``/``SfrStageTargets``/``SfrFlowTargets``/``DrnFlowTargets``
            (``simulated_series``), or a list of them: each set's own simulated values
            at period end, read as that set reads them.
        tolerance : float or None, optional
            The largest change between the last two cycles that counts as settled, in
            each source's own units (length for heads and stages). Default ``None``:
            numbers only, no verdict -- what counts as equilibrium is the modeller's
            call, relative to the target noise. Given, ``summary()`` gains ``passes``
            and :attr:`CycleEquilibrium.converged` answers. Must be finite and ``>= 0``.

        Returns
        -------
        CycleEquilibrium
            ``get()`` / ``summary()`` / ``plot()``.

        Raises
        ------
        ValueError
            When the year has fewer than two cycles (there is nothing to compare: run
            at least two), when the model's calendar is not this year's, or when
            ``tolerance`` is negative or not finite, or when ``targets=None`` on a
            model that is not a flow model or saves no head file.
        """

        return CycleEquilibrium(self, model, targets=targets, tolerance=tolerance)

    # -- identity ---------------------------------------------------------------

    def _key(self) -> tuple:
        return (self._cycles, self._start, self._steady, self._steady_length, self._nstp, self._tsmult)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AverageYear):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self) -> int:
        return hash(self._key())

    def __repr__(self) -> str:
        nstp = self.nstp
        steady = (
            f", steady_length={self._steady_length:g}" if self._steady else ", steady=False"
        )
        return (
            f"AverageYear({self._cycles}, start={format_start(self._start)!r}{steady}, "
            f"nstp={nstp!r}, tsmult={self._tsmult:g})  # nper={self.nper}, "
            f"{_month_name(self._months[0])}-{_month_name(self._months[-1])}"
        )


def _value_kind(value: Any) -> tuple[str, Any] | None:
    """``("number"|"array"|"dict", detail)`` for a value the annual mean can average, else None."""

    if _is_number(value):
        return ("number", None)
    if isinstance(value, np.ndarray) and value.dtype.kind in "iuf":
        return ("array", value.shape)
    if isinstance(value, Mapping) and value and all(_is_number(v) for v in value.values()):
        return ("dict", frozenset(value))
    return None


def _annual_mean(monthly: Mapping[int, Any]) -> Any:
    """The month-length-weighted mean of twelve numeric values (numbers, arrays, or dicts)."""

    kinds = {month: _value_kind(value) for month, value in monthly.items()}
    if any(kind is None for kind in kinds.values()) or len(set(kinds.values())) != 1:
        found = sorted({type(value).__name__ for value in monthly.values()})
        raise ValueError(
            "by_month() averages the twelve values for the steady period only when they are "
            "all numbers, all arrays of one shape, or all {key: number} dicts with the same "
            f"keys; these are {', '.join(found)}. Give the steady period's value: "
            "by_month(..., steady=<value>)."
        )
    weights = {month: NONLEAP_DAYS[month] / 365.0 for month in monthly}
    kind = next(iter(kinds.values()))[0]
    if kind == "number":
        return float(sum(weights[m] * float(v) for m, v in monthly.items()))
    if kind == "array":
        return sum(weights[m] * np.asarray(v, dtype=float) for m, v in monthly.items())
    keys = list(next(iter(monthly.values())))
    return {key: float(sum(weights[m] * float(v[key]) for m, v in monthly.items())) for key in keys}


class CycleEquilibrium:
    """How far a cyclic run's last two cycles differ, read at every period end.

    Made by :meth:`AverageYear.equilibrium`. Each SOURCE -- the model's heads, each
    LAK package's lake stages, or each target set given -- is read at the END of
    every period of the last two cycles (ledger 336: a period's value is its
    end-of-period value; the head file through
    :func:`~myflopy.modflow.mf6.period_end.period_ends`, a target set through its own
    ``simulated_heads``/``simulated_series``), and each month's value in the last
    cycle is compared with the same month's in the one before. A cell, lake or target
    with no value in either (inactive, dry) is left out of that month; one with a
    value in ONE of the two only is unsettled, and counted (``n_one_sided``).

    Parameters
    ----------
    year : AverageYear
        The cycle the run was built on; at least two cycles (``ValueError``).
    model : SimulationBase
        The run model; its :attr:`~SimulationBase.calendar` must be ``year``'s
        (same number of periods, start, unit and lengths, and -- for a flow model,
        whose calendar carries them -- the same steady flags), or ``ValueError``
        naming the first difference.
    targets : target set, sequence of target sets, or None, optional
        As :meth:`AverageYear.equilibrium` documents; default ``None`` reads every
        cell's head and every lake's stage. A set answering neither
        ``simulated_heads`` nor ``simulated_series`` raises ``TypeError``.
    tolerance : float or None, optional
        As :meth:`AverageYear.equilibrium` documents. ``None`` (default) gives
        numbers without a verdict.
    """

    def __init__(self, year: AverageYear, model: Any, *, targets: Any = None, tolerance: float | None = None):
        if not isinstance(year, AverageYear):
            raise TypeError(f"CycleEquilibrium takes an AverageYear, got {type(year).__name__}.")
        if year.cycles < 2:
            raise ValueError(
                f"{year!r} runs one cycle, and equilibrium compares the last two: declare "
                "AverageYear(cycles=2) or more and run it."
            )
        if tolerance is not None and (not _is_number(tolerance) or not np.isfinite(tolerance) or tolerance < 0):
            raise ValueError(f"equilibrium(tolerance={tolerance!r}) must be a finite number >= 0, or None.")
        on_model = period_calendar(model)
        declared = year.calendar
        differs = _calendar_difference(on_model, declared)
        if differs:
            raise ValueError(
                f"The model's periods are not {year!r}'s: {differs}. "
                "AverageYear.from_model(model) reads the model's own cycle."
            )
        self.year = year
        self.model = model
        self.tolerance = None if tolerance is None else float(tolerance)
        self._sources = _sources(model, targets)
        self._cache: dict[tuple[str, int], np.ndarray] = {}

    # -- the engine -------------------------------------------------------------

    def _values(self, source: _Source, per: int) -> np.ndarray:
        key = (source.label, per)
        if key not in self._cache:
            self._cache[key] = np.asarray(source.read(per), dtype=float)
        return self._cache[key]

    def _pair(self, previous: int, last: int) -> pd.DataFrame:
        """The long table of one cycle pair: one row per source, location and month."""

        before = self.year.months_of(previous)
        after = self.year.months_of(last)
        frames = []
        for source in self._sources:
            for month in self.year.months:
                old = self._values(source, before[month])
                new = self._values(source, after[month])
                # A location with a value in ONE cycle only -- a cell or lake dry or
                # inactive at one period end and not the other -- is the plainest
                # sign of a run still moving: kept, with no `change`, and counted by
                # summary() (ledger 531).
                keep = np.isfinite(old) | np.isfinite(new)
                frame = source.locations.loc[keep].copy()
                frame.insert(0, "source", source.label)
                frame["month"] = month
                frame["per_previous"] = before[month]
                frame["per_last"] = after[month]
                frame["previous"] = old[keep]
                frame["last"] = new[keep]
                frame["change"] = np.where(
                    np.isfinite(old[keep]) & np.isfinite(new[keep]), new[keep] - old[keep], np.nan
                )
                frames.append(frame)
        columns = ["source", "name", "layer", "cell", "month", "per_previous", "per_last",
                   "previous", "last", "change"]
        if not frames:
            return pd.DataFrame(columns=columns)
        table = pd.concat(frames, ignore_index=True)[columns]
        table["layer"] = table["layer"].astype("Int64")
        table["cell"] = table["cell"].astype("Int64")
        return table

    # -- the verbs --------------------------------------------------------------

    def get(self) -> pd.DataFrame:
        """Every location's change between the last two cycles, month by month.

        Returns
        -------
        pandas.DataFrame
            One row per source, location and calendar month: ``source`` (``"heads"``,
            ``"lake stage (<package>)"`` or the target set's label), ``name`` (a lake
            as ``"<package>:<lake>"``, zero-based; a target's label; ``None`` for a
            cell), ``layer``/``cell`` (zero-based; ``<NA>`` off the head grid),
            ``month`` (1-12), ``per_previous``/``per_last`` (the zero-based periods
            compared), ``previous``/``last`` (their period-end values) and
            ``change`` (``last - previous``). A location with no value in EITHER
            period (inactive, or dry or without stage at both ends) is not listed
            for that month; one with a value in ONE of them only -- dry at one
            period end and not at the other -- is listed with ``NaN`` on its empty
            side and in ``change``, and :meth:`summary` counts it as unsettled.
        """

        return self._pair(self.year.cycles - 2, self.year.cycles - 1)

    def summary(self) -> pd.DataFrame:
        """One row per source: the largest and mean change, where, and the annual-mean drift.

        Returns
        -------
        pandas.DataFrame
            ``source``, ``n`` (locations compared in at least one month),
            ``max_abs_change`` (the largest ``|change|`` over every month and
            location), ``month`` and ``at`` (where it is: ``"layer L, cell C"``, a
            lake or a target label), ``mean_abs_change`` (over every month and
            location), ``annual_mean_drift`` (the largest ``|change|`` of a
            location's month-length-weighted mean over the months it was compared,
            how far the year's average level still moves) and ``n_one_sided``
            (location-months with a value in one of the two cycles only -- a cell
            or lake dry at one period end and not the other; ledger
            531). With a ``tolerance``, ``passes``: ``max_abs_change <=
            tolerance`` with ``n_one_sided == 0`` and something compared.
        """

        table = self.get()
        rows = []
        for source in self._sources:
            rows_of = table[table["source"] == source.label]
            part = rows_of[rows_of["change"].notna()]
            row: dict[str, Any] = {"source": source.label}
            if part.empty:
                row.update(n=0, max_abs_change=np.nan, month=pd.NA, at=None,
                           mean_abs_change=np.nan, annual_mean_drift=np.nan)
            else:
                size = part["change"].abs()
                worst = part.loc[size.idxmax()]
                weights = part["month"].map(NONLEAP_DAYS).astype(float)
                keys = ["name", "layer", "cell"]
                grouped = (
                    part.assign(_w=weights, _wc=weights * part["change"])
                    .groupby(keys, dropna=False)[["_w", "_wc"]].sum()
                )
                drift = (grouped["_wc"] / grouped["_w"]).abs()
                row.update(
                    n=int(len(grouped)),
                    max_abs_change=float(size.max()),
                    month=int(worst["month"]),
                    at=_where(worst),
                    mean_abs_change=float(size.mean()),
                    annual_mean_drift=float(drift.max()),
                )
            row["n_one_sided"] = int(len(rows_of) - len(part))
            if self.tolerance is not None:
                row["passes"] = bool(
                    row["n"] and row["max_abs_change"] <= self.tolerance and not row["n_one_sided"]
                )
            rows.append(row)
        return pd.DataFrame(rows)

    @property
    def converged(self) -> bool:
        """Whether every source's largest change is within ``tolerance``, none one-sided.

        Raises ``ValueError`` without a ``tolerance``: the verdict is the
        modeller's, made by naming one.
        """

        if self.tolerance is None:
            raise ValueError(
                "No tolerance was given, so there is no verdict; "
                "year.equilibrium(model, tolerance=<length>) gives one."
            )
        return bool(self.summary()["passes"].all())

    def plot(self, *, backend: str = "plotly") -> Any:
        """The largest change per calendar month, for the last two cycle pairs.

        One line per source and cycle pair -- the last pair solid, the pair before it
        dashed when the year runs three cycles or more -- so a line falling between
        the pairs says the run is still settling. The tolerance, when given, is a
        dashed horizontal line.

        Parameters
        ----------
        backend : str, default "plotly"
            ``"plotly"`` returns an interactive ``viz.Fig``; ``"mpl"`` (or
            ``"matplotlib"``) a Matplotlib ``Figure``. Any other value raises
            ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            Max ``|change|`` against the calendar month, in cycle order.
        """

        renderer = normalize_backend(backend)
        cycles = self.year.cycles
        pairs = [(cycles - 2, cycles - 1)]
        if cycles >= 3:
            pairs.insert(0, (cycles - 3, cycles - 2))
        labels = [_month_name(m) for m in self.year.months]
        colors = viz.category_colors([source.label for source in self._sources])
        lines = []
        for previous, last in pairs:
            table = self._pair(previous, last)
            for source in self._sources:
                part = table[table["source"] == source.label]
                worst = part.assign(size=part["change"].abs()).groupby("month")["size"].max()
                y = [float(worst.get(month, np.nan)) for month in self.year.months]
                solid = last == cycles - 1
                lines.append((source.label, f"cycle {last} vs {previous}", y, solid))
        title = "Change between cycles at period end"
        if renderer == "mpl":
            fig, ax = viz.mpl_axes(figsize=(8, 4))
            for label, pair, y, solid in lines:
                ax.plot(labels, y, color=colors[label], ls="-" if solid else "--",
                        marker="o", label=f"{label}: {pair}")
            if self.tolerance is not None:
                ax.axhline(self.tolerance, color=viz.PALETTE.mpl_measured, ls=":", label="tolerance")
            ax.set_xlabel("month")
            ax.set_ylabel("max |change|")
            ax.set_title(title)
            ax.legend()
            return fig
        fig = viz.Fig()
        for label, pair, y, solid in lines:
            fig.add_scatter(
                x=labels, y=y, mode="lines+markers", name=f"{label}: {pair}",
                line=dict(color=colors[label], dash="solid" if solid else "dash"),
            )
        if self.tolerance is not None:
            fig.add_hline(y=self.tolerance, line=dict(color=viz.PALETTE.measured, dash="dot"),
                          annotation_text="tolerance")
        fig.update_layout(title=title, xaxis_title="month", yaxis_title="max |change|", dragmode="pan")
        return fig

    def __repr__(self) -> str:
        verdict = "" if self.tolerance is None else f", tolerance={self.tolerance:g}"
        return (
            f"CycleEquilibrium({self.year!r}{verdict}; sources: "
            f"{', '.join(source.label for source in self._sources)})"
        )


def _calendar_difference(on_model: PeriodCalendar, declared: PeriodCalendar) -> str | None:
    """What makes a model's calendar another year's, in words, or ``None`` when it is this one.

    The steady flags are compared only when the model carries them: a transport
    model's periods are as steady as the flow model it reads (``steady`` is ``None``).
    """

    if on_model.nper != declared.nper:
        return f"it runs {on_model.nper} period(s) where the year runs {declared.nper}"
    if on_model.start != declared.start:
        return (
            f"it starts {format_start(on_model.start)} where the year starts "
            f"{format_start(declared.start)}"
        )
    if on_model.time_units != declared.time_units:
        return f"it counts time in {on_model.time_units} where the year counts days"
    if not on_model.same_clock(declared):
        lengths = np.abs(np.asarray(on_model.perlen) - np.asarray(declared.perlen))
        per = int(np.argmax(lengths))
        return (
            f"its period {per} is {on_model.perlen[per]:g} {on_model.time_units} long where "
            f"the year's is {declared.perlen[per]:g} days"
        )
    if on_model.steady is not None and on_model.steady != declared.steady:
        per = next(p for p, (a, b) in enumerate(zip(on_model.steady, declared.steady, strict=True)) if a != b)
        state = "steady state" if on_model.steady[per] else "transient"
        return f"its STO runs period {per} {state}, where the year's does not"
    return None


def _where(row: pd.Series) -> str:
    if row["name"] is not None and not (isinstance(row["name"], float) and np.isnan(row["name"])):
        return str(row["name"])
    return f"layer {int(row['layer'])}, cell {int(row['cell'])}"


class _Source:
    """One thing compared: a label, its locations and a reader of one period's values."""

    def __init__(self, label: str, locations: pd.DataFrame, read):
        self.label = label
        self.locations = locations.reset_index(drop=True)
        self.read = read


def _sources(model: Any, targets: Any) -> list[_Source]:
    if targets is None:
        return [_head_source(model), *_lake_sources(model)]
    sets = list(targets) if isinstance(targets, (list, tuple)) else [targets]
    if not sets:
        raise ValueError("equilibrium(targets=) is empty; pass None to compare everywhere.")
    sources = []
    seen: dict[str, int] = {}
    for target_set in sets:
        name = type(target_set).__name__
        seen[name] = seen.get(name, 0) + 1
        label = name if seen[name] == 1 else f"{name}#{seen[name]}"
        sources.append(_target_source(model, target_set, label))
    return sources


def _head_source(model: Any) -> _Source:
    """Every cell's head at each period end, through the one period-end rule.

    The head file is FloPy's reader of the file the model's OC names
    (``model.gwf.output.head()``), so a DIS, DISV or DISU model is read alike;
    ``model.hds`` is DISV-only (ledger 531).
    """

    if str(getattr(model, "model_type", "gwf6") or "gwf6").lower() != "gwf6":
        raise ValueError(
            f"equilibrium(targets=None) compares a FLOW model's heads, and this is a "
            f"{getattr(model, 'model_type', '?')} model: pass the GWF model, or pass this "
            "model's targets (ConcTargets) as targets=."
        )
    output = getattr(getattr(model, "gwf", None), "output", None)
    heads = output.head() if output is not None else None
    if heads is None:
        raise ValueError(
            "The model saves no head file; equilibrium reads every period of the last two "
            "cycles: mf.oc(head_filerecord=..., saverecord=[('HEAD', 'LAST')])."
        )
    ends = period_ends(model, reader=lambda: heads, file="the head file", saves="HEAD")
    if not ends:
        raise ValueError(
            "The head file holds no record; equilibrium reads every period of the last two "
            "cycles: save them with mf.oc(saverecord=[('HEAD', 'LAST')])."
        )
    first = np.asarray(heads.get_data(kstpkper=next(iter(ends.values()))), dtype=float)
    nlay = int(first.shape[0]) if first.ndim > 1 else 1
    ncpl = first.size // nlay
    layer, cell = np.divmod(np.arange(nlay * ncpl), ncpl)
    locations = pd.DataFrame({"name": None, "layer": layer, "cell": cell})

    def read(per: int) -> np.ndarray:
        if per not in ends:
            raise ValueError(
                f"The head file holds no record of period {per}; equilibrium reads every period "
                "of the last two cycles: save them with mf.oc(saverecord=[('HEAD', 'LAST')])."
            )
        values = np.asarray(heads.get_data(kstpkper=ends[per]), dtype=float).ravel()
        values[_is_sentinel(values)] = np.nan
        return values

    return _Source("heads", locations, read)


def _lake_sources(model: Any) -> list[_Source]:
    """Each LAK package's lakes at each period end; a package saving no stage file is warned about."""

    sources = []
    for package in model_packages_of_type(model, "lak"):
        name = str(package.package_name)
        record = getattr(package, "stage_filerecord", None)
        held = record.get_data() if record is not None and hasattr(record, "get_data") else None
        if held is None or len(held) == 0:
            message = (
                f"LAK package {name!r} saves no stage file, so its lakes are not in the "
                "equilibrium report; declare stage_filerecord on it, or pass its lakes as "
                "LakeStageTargets."
            )
            logger.warning(message)
            warnings.warn(message, stacklevel=4)
            continue
        frame = build_lak_stage_result_table(model, package=name)
        stages = _stage_at_period_end(
            model, frame, reader=lambda package=package: package.output.stage(),
            file=f"the stage file of LAK package {name!r}", id_column="lake",
        )
        table = stages.pivot(index="per", columns="lake", values="stage")
        lakes = list(table.columns)
        locations = pd.DataFrame(
            {"name": [f"{name}:{int(lake)}" for lake in lakes], "layer": None, "cell": None}
        )

        def read(per: int, table=table, name=name) -> np.ndarray:
            if per not in table.index:
                raise ValueError(
                    f"The stage file of LAK package {name!r} holds no record of period {per}."
                )
            return table.loc[per].to_numpy(dtype=float)

        sources.append(_Source(f"lake stage ({name})", locations, read))
    return sources


def _target_source(model: Any, target_set: Any, label: str) -> _Source:
    """A target set's own simulated values at each period end, read as the set reads them."""

    if hasattr(target_set, "simulated_heads"):
        table = target_set.simulated_heads(model)
        key = "per"
    elif hasattr(target_set, "simulated_series"):
        table = target_set.simulated_series(model)
        key = "time"
    else:
        raise TypeError(
            f"equilibrium(targets=) takes target sets (HeadTargets, LakeStageTargets, ...); "
            f"{type(target_set).__name__} answers neither simulated_heads nor simulated_series."
        )
    table = table.set_index(key)
    table.index = [int(round(float(value))) for value in table.index]
    names = [str(column) for column in table.columns]
    locations = pd.DataFrame({"name": names, "layer": None, "cell": None})

    def read(per: int) -> np.ndarray:
        if per not in table.index:
            raise ValueError(f"{label} has no simulated value at period {per}.")
        values = table.loc[per].to_numpy(dtype=float)
        values[_is_sentinel(values)] = np.nan
        return values

    return _Source(label, locations, read)


__all__ = ["NONLEAP_DAYS", "AverageYear", "CycleEquilibrium"]
