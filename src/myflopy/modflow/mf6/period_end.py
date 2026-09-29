"""Which saved record IS a stress period's end: the one rule every target reads (ledger 336).

A target keyed by stress period reads ONE value per period, and that value is the
period's END -- the user's rule, for a flow as for a head or a stage. A binary
output file can hold several records per period (``SAVE ... ALL`` at ``nstp > 1``),
one (``LAST``), or ones that stop short of the end (``FIRST``, ``STEPS 1 3`` at
``nstp=5``, ``FREQUENCY``); which one ends the period is read from the file's own
clock, never from another file's (ledger 301: a budget file and a head file are
written on independent MF6 schedules).

A record ends its period when its ``pertim`` -- the time into the period, which MF6
writes in every record header -- is the period's TDIS ``perlen``. NOT when its
``kstp`` is the period's TDIS ``nstp``: under ATS the step count is decided at run
time (measured: two periods declared ``nstp=1`` took 7 and 6 steps), and MF6 still
writes the final step's ``pertim`` as ``perlen`` exactly (measured on ``tsmult`` 1.2
and 0.9 periods and on those two ATS periods).

A period whose end the file does not hold is REFUSED, naming the period, the steps
output control did save and the setting that fixes it: the value asked for does
not exist, and an earlier step is a different number.

An MF6 CONTINUOUS OBSERVATION file carries no ``kstp``/``kper`` at all -- one row per
time step, labelled by simulated time -- so :func:`period_end_rows` works out which
row ends each period from those times and the TDIS, and hands the result to the same
check (ledger 346).

Deliberately light -- numpy, pandas and myflopy's logger at import -- because PEST's
forward run imports it on every model run (``pest/forward_run.py``), where importing
the observation package would cost ~2.4 s a run for nothing. :func:`period_lengths`
imports FloPy's exception type when it is CALLED, and the forward run never calls it:
the build hands it each period's length.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from myflopy._logging import get_logger

logger = get_logger(__name__)

#: How close a record's ``pertim`` must be to its period's ``perlen`` to END that
#: period, as a fraction of ``perlen`` floored at one time unit. The gap it absorbs
#: is FloPy's WRITE rounding, not a reader's: MF6 writes the final step's
#: ``pertim`` as the ``perlen`` it READ from the TDIS file, and FloPy writes that
#: file at 8 decimals (``{:15.8f}`` from 0.001 to 1e5, ``{:.8E}`` outside), while a
#: caller compares the ``perlen`` FloPy was GIVEN. That is at most 5e-9 x
#: max(1, ``perlen``) -- measured: 1/24 d written ``0.04166667``, 3.3e-9 off;
#: 100000.00049 written ``1.00000000E+05``, 4.9e-9 relative -- and this is twice it.
#: At 1e-9 a live model with an hourly, 1/3- or 4/3-unit period refused every
#: period while its file held each one's end (ledger 336). The step before the end
#: is short of it by a whole time step, so a record reads as the end wrongly only
#: when that last step is shorter than 1e-8 of the period (``tsmult=0.5`` over 27+
#: steps); a FloPy ``float_precision`` set below its default of 8 is not covered.
END_RTOL = 1e-8


def _end_tolerance(length: float) -> float:
    """How far a time into a period may sit from the period's length and still END it.

    :data:`END_RTOL` of the length, floored at one time unit -- the one tolerance
    both :func:`period_end_steps` (does this record end its period?) and
    :func:`period_end_rows` (is this row past its period's end?) apply.
    """

    return END_RTOL * max(abs(float(length)), 1.0)


def saved_steps(reader) -> list[tuple[int, int, float]]:
    """Every time a binary MF6 output file saved, on that file's own clock.

    Parameters
    ----------
    reader : flopy.utils.HeadFile or flopy.utils.CellBudgetFile
        An open reader of the file -- a head, concentration or stage file, the
        model's cell-budget file, or an advanced package's budget file. Its
        ``recordarray`` is read, not its data.

    Returns
    -------
    list of tuple of (int, int, float)
        ``(kstp, kper, pertim)`` for each distinct saved time, zero-based (the
        file is one-based) and in file order. ``pertim`` is the time into the
        stress period at the end of that step, as MF6 wrote it in the record
        header.

    Examples
    --------
    >>> saved_steps(model.fl.sfr.output.budget())
    [(0, 0, 3.333...), (1, 0, 6.666...), (2, 0, 10.0), (0, 1, 3.333...), ...]
    """

    records = reader.recordarray
    seen: dict[tuple[int, int], float] = {}
    for kstp, kper, pertim in zip(
        records["kstp"], records["kper"], records["pertim"], strict=True
    ):
        seen.setdefault((int(kstp) - 1, int(kper) - 1), float(pertim))
    return [(kstp, kper, pertim) for (kstp, kper), pertim in seen.items()]


def last_saved_steps(saved: Iterable[Sequence]) -> dict[int, tuple[int, int]]:
    """Each stress period's LAST saved step -- the one rule for it (ledger 337).

    What ``per=N`` means wherever a picture or a table takes one: the last step
    the file saved in stress period ``N``, zero-based. Every reader of that rule
    calls this function (:func:`last_saved_step` for one period,
    :func:`period_end_steps` before it checks that the step ENDS its period, the
    map's ``per``, both sections, ``surface``, the section composers, a group's
    ``per`` and its period-end reduction, ``HeadsPlus.array``, and through
    :func:`at_last_saved_step`/:func:`at_last_saved_steps` every results map, series,
    stage change and group stage table, ledger 348); there were five
    copies of it, and a sixth reading ``per`` as a SAVED-RECORD index, so
    ``per=1`` drew step (2, 1) on one picture and (1, 0) on the next.

    Parameters
    ----------
    saved : iterable of tuple
        The file's saved times, zero-based: ``(kstp, kper)`` pairs -- a reader's
        or a model's ``kstpkper``, or a table's ``kstpkper`` column -- or
        ``(kstp, kper, pertim)`` as :func:`saved_steps` returns them (``pertim``
        is not read here).

    Returns
    -------
    dict of int to tuple of (int, int)
        ``{kper: (kstp, kper)}``, sorted by period: the highest-numbered step
        saved in each period that saved anything. A period that saved nothing is
        absent.

    Examples
    --------
    >>> last_saved_steps([(0, 0), (2, 0), (1, 0), (0, 1)])
    {0: (2, 0), 1: (0, 1)}
    """

    last: dict[int, int] = {}
    for item in saved:
        kstp, kper = int(item[0]), int(item[1])
        if kstp > last.get(kper, -1):
            last[kper] = kstp
    return {kper: (kstp, kper) for kper, kstp in sorted(last.items())}


def last_saved_step(saved: Iterable[Sequence], per: int) -> tuple[int, int]:
    """The last step saved in stress period ``per``, refusing a period with none.

    :func:`last_saved_steps` for one period: what ``per=`` selects on every
    picture and table that takes it (ledger 337, the user's rule -- ``per=0`` is
    the last saved step of the FIRST stress period).

    Parameters
    ----------
    saved : iterable of tuple
        The file's saved times, zero-based, as :func:`last_saved_steps` takes them.
    per : int
        The stress period, zero-based. Not a position in ``saved``: at ``nstp=3``
        with every step saved, ``per=1`` is ``(2, 1)``, the sixth saved time.

    Returns
    -------
    tuple of (int, int)
        ``(kstp, kper)``, zero-based.

    Raises
    ------
    ValueError
        When the file saved nothing in ``per`` -- a period past the end, a
        negative one, or one output control skipped -- naming the periods it has.

    Examples
    --------
    >>> last_saved_step([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)], per=1)
    (2, 1)
    """

    steps = last_saved_steps(saved)
    try:
        return steps[int(per)]
    except KeyError:
        raise ValueError(
            f"Stress period {per} is unavailable. Available periods: {sorted(steps)}"
        ) from None


def period_end_steps(
    saved: Iterable[Sequence],
    *,
    perlen: Sequence[float] | None = None,
    file: str = "the output file",
    saves: str = "BUDGET",
    fix: str | None = None,
) -> dict[int, tuple[int, int]]:
    """The saved record at each stress period's END, refusing a period whose end is missing.

    Each period that saved anything is represented by its LAST saved step. With
    ``perlen`` that step must also END the period -- its ``pertim`` equal to the
    period's length, within :data:`END_RTOL` -- or the period is refused: output
    control saved only steps before the end (``FIRST``, ``STEPS``, ``FREQUENCY``),
    so the period-end value does not exist in the file, and an earlier step is a
    different number (ledger 336). A period that saved NOTHING is simply absent:
    no value is read for it, so none can be wrong.

    Parameters
    ----------
    saved : iterable of tuple
        The file's saved times, zero-based: ``(kstp, kper, pertim)`` as
        :func:`saved_steps` returns them, or bare ``(kstp, kper)`` pairs when no
        ``perlen`` is given (a table's ``kstpkper`` column, with nothing to check
        it against).
    perlen : sequence of float, optional
        Each stress period's length, from the simulation's TDIS, in the model's
        time units. ``None`` (the default) skips the check and takes each
        period's last saved step as given: the route for a table with no file
        behind it, never for a run model.
    file : str, default "the output file"
        What the file is, for the refusal (``"the SFR budget file"``).
    saves : str, default "BUDGET"
        The output-control word that schedules this file -- ``"BUDGET"`` for a
        budget file, ``"HEAD"`` for a head or stage file (LAK and SFR write their
        stage when MF6 saves HEAD), ``"CONCENTRATION"`` or ``"TEMPERATURE"`` for a
        concentration or temperature file -- for the fix the refusal names.
    fix : str or None, optional
        The refusal's last sentence: what to do about a missing end. Default
        ``None`` names ``SAVE <saves> LAST`` in the model's output control, the fix
        for a binary output file. :func:`period_end_rows` passes its own, because
        output control does not schedule a continuous observation file -- MF6
        writes one at every time step -- and ``saves`` is then not read.

    Returns
    -------
    dict of int to tuple of (int, int)
        ``{kper: (kstp, kper)}``, zero-based, sorted by period: the record that
        ends each period the file saved.

    Raises
    ------
    ValueError
        With ``perlen``: when a period's last saved step does not end it --
        naming up to five such periods, the steps saved (one-based, as output
        control counts them), where the last of them ends against the period's
        length, and the ``SAVE <saves> LAST`` setting that fixes it; when a
        saved time carries no ``pertim``; or when the file names a period the
        TDIS does not have.

    Examples
    --------
    >>> from myflopy.modflow.mf6.period_end import period_end_steps, saved_steps
    >>> period_end_steps(saved_steps(reader), perlen=[10.0, 10.0, 10.0],
    ...                  file="the SFR budget file", saves="BUDGET")
    {0: (2, 0), 1: (2, 1), 2: (2, 2)}
    >>> period_end_steps([(0, 0), (2, 0), (1, 0)])       # nothing to check against
    {0: (2, 0)}
    """

    saved = [tuple(item) for item in saved]
    ends = last_saved_steps(saved)
    if perlen is None:
        return ends
    steps: dict[int, dict[int, float | None]] = {}
    for item in saved:
        kstp, kper = int(item[0]), int(item[1])
        pertim = float(item[2]) if len(item) > 2 else None
        steps.setdefault(kper, {})[kstp] = pertim

    lengths = [float(value) for value in perlen]
    short: list[str] = []
    for kper, (kstp, _) in ends.items():
        if kper >= len(lengths):
            raise ValueError(
                f"{file} holds a record in stress period {kper} (zero-based), but the "
                f"simulation's TDIS declares {len(lengths)} periods."
            )
        pertim = steps[kper][kstp]
        if pertim is None:
            raise ValueError(
                f"Cannot tell which record of {file} ends stress period {kper}: its "
                "saved times carry no pertim. Pass them as saved_steps(reader) "
                "returns them."
            )
        length = lengths[kper]
        if abs(pertim - length) <= _end_tolerance(length):
            continue
        saved_here = sorted(steps[kper])
        noun = "step" if len(saved_here) == 1 else "steps"
        listing = ", ".join(str(step + 1) for step in saved_here)
        ended, long = _told_apart(pertim, length)
        short.append(
            f"period {kper} saved {noun} {listing}, the last ending {ended} into "
            f"a period {long} long"
        )
    if short:
        more = len(short) - 5
        listing = "; ".join(short[:5]) + (f"; and {more} more" if more > 0 else "")
        if fix is None:
            fix = (
                f"Save each period's last step: SAVE {saves} LAST (or ALL) in the model's "
                f"output control -- mf.oc(saverecord=[..., ('{saves}', 'LAST')])."
            )
        raise ValueError(
            f"A stress period's END is missing from {file}: {listing} (periods "
            "zero-based, steps one-based as output control counts them). A period's "
            "value is its end-of-period value (ledger 336), which this file does not "
            f"hold, and an earlier step is a different number. {fix}"
        )
    return ends


#: What a continuous-observation file short of a period's end means, for the refusal.
_OBSERVATION_FIX = (
    "MF6 writes a continuous observation at the end of EVERY time step, whatever output "
    "control saves, so a file short of a period's end is a run that stopped before it "
    "(read the model's listing file) or a file cut after MF6 wrote it. Run the model "
    "to its end."
)


def period_end_rows(
    totim: Sequence[float],
    perlen: Sequence[float],
    *,
    file: str = "the observation file",
) -> dict[int, int]:
    """Which row of an MF6 continuous-observation file ENDS each stress period (ledger 346).

    MF6 writes a continuous observation at the end of every time step, one row per
    step, labelled by the simulated time ``totim`` and carrying no step or period
    number. So the rows are assigned to periods here, from those times and the
    TDIS, and each period's last row is checked for being its END by the one rule
    every target reads (:func:`period_end_steps`, :data:`END_RTOL`): the observation
    route is keyed by zero-based stress period, one row per period at its end, like
    every other route.

    **The time into a period is taken by DIFFERENCING MF6's own clock**, never by
    comparing ``totim`` with the cumulative ``perlen``: MF6 accumulates the lengths
    it READ, which FloPy wrote at 8 decimals, so on a live model with hourly
    periods ``totim`` runs ahead of the cumulative length the model holds by
    3.3e-9 more every period and outgrows :data:`END_RTOL` from the fourth period
    (``docs/reproducers/obs_route_period_end.py``). Each period's time is measured
    from the row that ENDS the period before it -- MF6's own ``totim`` there -- so
    the rounding never accumulates past one period's.

    Parameters
    ----------
    totim : sequence of float
        The file's ``time`` column: MF6's simulated time at the end of each time
        step, in file order, in the model's time units. MF6 writes a TEXT file's
        time at 13 significant digits (a binary file's in full; OBS has no
        single-precision output), and that rounding stays inside :data:`END_RTOL`
        for any realistic model -- ten years of hourly periods, 87,600 ends, read
        each within 4.4e-9 of its length -- but not for a period shorter than one
        time unit past ``totim`` ~1e4, where a row can be put in the wrong period
        and the file refused (ledger 347).
    perlen : sequence of float
        Each stress period's length from the simulation's TDIS
        (:func:`period_lengths`), one per period, in the same time units.
    file : str, default "the observation file"
        What the file is, for the refusals (``"the observation file drn.csv
        (package 'drn')"``).

    Returns
    -------
    dict of int to int
        ``{kper: row}``: for EVERY stress period, zero-based and in order, the
        position (from 0) of the row that ends it.

    Raises
    ------
    ValueError
        When a time is not a number or the times decrease (not MF6's file); when a
        row lies past the end of the last period (a file from another simulation);
        when a period's last row does not end it -- the run stopped inside it, or the
        file was cut -- naming the period, the steps it holds and where the last of
        them ends (:func:`period_end_steps`); or when a period holds no row at all,
        naming it and where the file stops.

    Examples
    --------
    >>> period_end_rows([10 / 3, 20 / 3, 10.0, 40 / 3, 50 / 3, 20.0], [10.0, 10.0])
    {0: 2, 1: 5}
    """

    times = np.asarray(pd.to_numeric(pd.Series(list(totim), dtype=object), errors="coerce"),
                       dtype=float)
    if not np.isfinite(times).all():
        bad = int(np.flatnonzero(~np.isfinite(times))[0])
        raise ValueError(
            f"{file} holds a time that is not a number at row {bad + 1}: MF6 labels each "
            "row of a continuous observation file with its simulated time."
        )
    if times.size > 1 and (np.diff(times) < 0).any():
        bad = int(np.flatnonzero(np.diff(times) < 0)[0]) + 1
        raise ValueError(
            f"{file}'s times go back at row {bad + 1} ({float(times[bad - 1])!r} then "
            f"{float(times[bad])!r}): MF6 writes one row per time step in time order, so this "
            "is not the file MF6 wrote."
        )
    lengths = [float(value) for value in perlen]
    saved: list[tuple[int, int, float]] = []
    position: dict[tuple[int, int], int] = {}
    row = 0
    start = 0.0
    for kper, length in enumerate(lengths):
        kstp = 0
        while row < times.size and times[row] - start <= length + _end_tolerance(length):
            pertim = float(times[row] - start)
            saved.append((kstp, kper, pertim))
            position[(kstp, kper)] = row
            kstp += 1
            row += 1
            if abs(pertim - length) <= _end_tolerance(length):
                # The row that ends a period is its last: the next one, even at the
                # same time (a period of length 0 after it), is the next period's.
                break
        if kstp and abs(saved[-1][2] - length) <= _end_tolerance(length):
            # The period ends on MF6's own clock here; the next one is measured from it.
            start = float(times[row - 1])
        else:
            # No end to measure from: its nominal end, so the periods after it are
            # still told apart (the check below refuses this one).
            start += length
    if row < times.size:
        raise ValueError(
            f"{file} holds {times.size - row} row(s) past the end of the simulation's last "
            f"stress period (t = {float(times[row])!r}, where the {len(lengths)} periods of its "
            f"TDIS end at {sum(lengths)!r}): it was not written by this simulation."
        )
    ends = period_end_steps(saved, perlen=lengths, file=file, fix=_OBSERVATION_FIX)
    absent = [kper for kper in range(len(lengths)) if kper not in ends]
    if absent:
        where = (
            f"its last row is t = {float(times[-1])!r}, in period {saved[-1][1]}"
            if saved else "it holds no row at all"
        )
        raise ValueError(
            f"{file} holds no row in stress period(s) "
            f"{', '.join(str(kper) for kper in absent[:10])}"
            f"{' ...' if len(absent) > 10 else ''} (zero-based) of the {len(lengths)} its "
            f"TDIS declares: {where}. {_OBSERVATION_FIX}"
        )
    return {kper: position[step] for kper, step in ends.items()}


def _told_apart(first: float, second: float) -> tuple[str, str]:
    """Two numbers printed with the fewest digits (six at least) that tell them apart.

    A refusal that printed both at six figures read "ending 0.0416667 into a period
    0.0416667 long" for two numbers 3e-9 apart; this prints as many as it takes.

    Parameters
    ----------
    first : float
        The first number (a record's ``pertim``).
    second : float
        The second number (its period's ``perlen``).

    Returns
    -------
    tuple of (str, str)
        Both numbers in ``g`` format at the same number of significant digits, the
        fewest from 6 to 17 at which the two strings differ; their ``repr`` when
        even 17 digits print them alike.

    Examples
    --------
    >>> _told_apart(9.9999998, 10.0)
    ('9.9999998', '10')
    """

    for digits in range(6, 18):
        a, b = f"{first:.{digits}g}", f"{second:.{digits}g}"
        if a != b:
            return a, b
    return repr(first), repr(second)


def period_lengths(model: Any) -> list[float] | None:
    """Each stress period's length from ``model``'s TDIS, or ``None`` when it has none.

    Every run model carries a TDIS (MF6 will not run without one) -- a live one
    (``run.model()``) and one read back from its files (``mf.load_mf6_run``) alike --
    so ``None`` means a duck-typed double, a table with no simulation behind it, and
    the period-end rule then has nothing to check a file against. The lengths are the
    ones FloPy holds, which for a live model differ from what MF6 read by FloPy's
    8-decimal write rounding; :data:`END_RTOL` absorbs that.

    Parameters
    ----------
    model : SimulationBase or object
        A run model (``model.sim.tdis``), or a double without one.

    Returns
    -------
    list of float or None
        ``perlen`` per stress period, zero-based, in the model's time units.
    """

    from flopy.mf6.mfbase import MFDataException

    tdis = getattr(getattr(model, "sim", None), "tdis", None)
    if tdis is None:
        return None
    try:
        data = tdis.perioddata.get_data()
    except (AttributeError, MFDataException):
        # A TDIS object with no period block behind it: nothing declares how long
        # a period is, so there is nothing to check a file's times against.
        logger.debug("TDIS carries no perioddata; the period-end rule goes unchecked")
        return None
    if data is None:
        return None
    return [float(value) for value in data["perlen"]]


def period_ends(
    model: Any,
    *,
    reader: Callable[[], Any],
    file: str,
    saves: str,
    fallback: Iterable[Sequence] = (),
) -> dict[int, tuple[int, int]]:
    """The record that ends each stress period in one of ``model``'s output files.

    :func:`period_end_steps` read on THAT file's own clock (ledger 301), checked
    against the model's TDIS (:func:`period_lengths`): a period whose end the file
    does not hold is refused.

    Parameters
    ----------
    model : SimulationBase or object
        The run model the file belongs to.
    reader : callable
        Takes no arguments and returns the file's open FloPy reader
        (``lambda: model.fl.sfr.output.budget()``). Called only when the model has a
        TDIS -- every run model does -- so a double without the file is never asked
        for it.
    file : str
        What the file is, for the refusal (``"the SFR budget file"``).
    saves : str
        The output-control word scheduling it (``"BUDGET"``, ``"HEAD"``,
        ``"CONCENTRATION"``, ``"TEMPERATURE"``).
    fallback : iterable of tuple, optional
        The table's own ``(kstp, kper)`` labels, used UNCHECKED when the model has
        no TDIS (a double). Empty by default.

    Returns
    -------
    dict of int to tuple of (int, int)
        ``{kper: (kstp, kper)}``, zero-based.

    Raises
    ------
    ValueError
        When a period the file saved lacks its final step (see
        :func:`period_end_steps`).

    Examples
    --------
    >>> period_ends(model, reader=model._get_budget_reader,
    ...             file="the model's cell-budget file", saves="BUDGET")
    {0: (2, 0), 1: (2, 1), 2: (2, 2)}
    """

    perlen = period_lengths(model)
    if perlen is None:
        return period_end_steps(fallback)
    return period_end_steps(saved_steps(reader()), perlen=perlen, file=file, saves=saves)


def at_period_end(
    frame: pd.DataFrame,
    ends: Mapping[int, tuple[int, int]],
    *,
    column: str = "kstpkper",
) -> pd.DataFrame:
    """The rows of ``frame`` saved at a period's end.

    Parameters
    ----------
    frame : pandas.DataFrame
        A table labelled by its own file's clock: one ``(kstp, kper)`` tuple per
        row, zero-based, in ``column``.
    ends : mapping of int to tuple of (int, int)
        The records that end each period, as :func:`period_end_steps` returns
        them.
    column : str, default "kstpkper"
        The column holding each row's ``(kstp, kper)``.

    Returns
    -------
    pandas.DataFrame
        The rows whose ``(kstp, kper)`` is one of ``ends``' values, in their
        original order and with their original index.

    Examples
    --------
    >>> table = model.outputs.sfr.bud.get("FLOW-JA-FACE")
    >>> at_period_end(table, {0: (2, 0), 1: (2, 1)})
    """

    keep = set(ends.values())
    mask = np.fromiter(
        ((int(value[0]), int(value[1])) in keep for value in frame[column]),
        dtype=bool,
        count=len(frame),
    )
    return frame.loc[mask]


def at_last_saved_step(
    frame: pd.DataFrame,
    per: int,
    *,
    column: str = "kstpkper",
) -> pd.DataFrame:
    """The rows of ``frame`` saved at stress period ``per``'s LAST saved step (ledger 337).

    What ``per=N`` draws on a picture -- a results map, the exchange hover's stage --
    read on the TABLE's own clock, the file it was read from (ledger 301): a LAK or SFR
    stage table follows ``SAVE HEAD``, a budget table ``SAVE BUDGET``. The steps a
    period saved are never combined: at ``nstp=3`` under ``SAVE BUDGET ALL`` a map
    that took every row of the period ADDED three flows (ledger 348). It is
    :func:`last_saved_step` and :func:`at_period_end` together, so the rule stays
    the one function.

    Parameters
    ----------
    frame : pandas.DataFrame
        A table labelled by its own file's clock: one ``(kstp, kper)`` tuple per
        row, zero-based, in ``column``. It may hold every period or only
        ``per``'s rows (a ``get(per=N)``); the step is found among the rows given.
    per : int
        The stress period, zero-based.
    column : str, default "kstpkper"
        The column holding each row's ``(kstp, kper)``.

    Returns
    -------
    pandas.DataFrame
        The rows saved at that one step, in their original order and with their
        original index.

    Raises
    ------
    ValueError
        When no row of ``frame`` was saved in ``per`` (an empty frame included),
        naming the periods it has -- :func:`last_saved_step`'s refusal.

    Examples
    --------
    >>> table = model.lak.results.q.get(per=1)         # three steps at nstp=3
    >>> at_last_saved_step(table, 1)["kstpkper"].unique()
    [(2, 1)]
    """

    step = last_saved_step(frame[column], per)
    return at_period_end(frame, {int(per): step}, column=column)


def at_last_saved_steps(
    frame: pd.DataFrame,
    *,
    column: str = "kstpkper",
    by: str | None = None,
) -> pd.DataFrame:
    """The rows of ``frame`` saved at each stress period's LAST saved step (ledger 337).

    The reduction every per-period view of a results table reads before it combines
    rows -- a series by stress period, a stage change between periods, a group
    table: a period is represented by the last step its file saved, and its other
    steps are dropped rather than added or averaged (ledger 348). Deliberately the
    UNCHECKED rule, as every picture's is: a run saved with ``SAVE HEAD FIRST`` still
    draws; the targets read :func:`period_end_steps`, which refuses it.

    Parameters
    ----------
    frame : pandas.DataFrame
        A table labelled by its own file's clock: one ``(kstp, kper)`` tuple per
        row, zero-based, in ``column``.
    column : str, default "kstpkper"
        The column holding each row's ``(kstp, kper)``.
    by : str or None, optional
        A column whose groups each keep THEIR OWN last steps -- ``"model"`` for a
        group table, whose members may save different steps (``nstp=3`` beside
        ``nstp=2``). Default ``None``: the whole frame is one clock. A column the
        frame does not have raises ``KeyError``.

    Returns
    -------
    pandas.DataFrame
        The rows saved at a period's last saved step (of their own ``by`` group),
        in their original order and with their original index. A row with no
        record at that step -- a cell present only earlier in the period -- is
        dropped, never read at the earlier step. An empty frame comes back as it is.

    Examples
    --------
    >>> table = model.lak.results.stage.get()          # nine records at nstp=3
    >>> sorted(set(at_last_saved_steps(table)["kstpkper"]))
    [(2, 0), (2, 1), (2, 2)]
    """

    if frame.empty:
        return frame
    keys = [(int(value[0]), int(value[1])) for value in frame[column]]
    if by is None:
        keep = set(last_saved_steps(keys).values())
        mask = [key in keep for key in keys]
    else:
        groups = frame[by].to_numpy()
        owned: dict[object, list[tuple[int, int]]] = {}
        for key, group in zip(keys, groups, strict=True):
            owned.setdefault(group, []).append(key)
        ends = {group: set(last_saved_steps(saved).values()) for group, saved in owned.items()}
        mask = [key in ends[group] for key, group in zip(keys, groups, strict=True)]
    return frame.loc[np.asarray(mask, dtype=bool)]


__all__ = [
    "END_RTOL",
    "at_last_saved_step",
    "at_last_saved_steps",
    "at_period_end",
    "last_saved_step",
    "last_saved_steps",
    "period_end_rows",
    "period_end_steps",
    "period_ends",
    "period_lengths",
    "saved_steps",
]
