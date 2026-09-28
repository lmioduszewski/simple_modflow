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
    ``per`` and its period-end reduction, ``HeadsPlus.array``); there were five
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
        if abs(pertim - length) <= END_RTOL * max(abs(length), 1.0):
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
        raise ValueError(
            f"A stress period's END is missing from {file}: {listing} (periods "
            "zero-based, steps one-based as output control counts them). A period's "
            "value is its end-of-period value (ledger 336), which this file does not "
            "hold, and an earlier step is a different number. Save each period's last "
            f"step: SAVE {saves} LAST (or ALL) in the model's output control -- "
            f"mf.oc(saverecord=[..., ('{saves}', 'LAST')])."
        )
    return ends


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


__all__ = [
    "END_RTOL",
    "at_period_end",
    "last_saved_step",
    "last_saved_steps",
    "period_end_steps",
    "period_ends",
    "period_lengths",
    "saved_steps",
]
