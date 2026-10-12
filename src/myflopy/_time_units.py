"""A model's time unit: which spellings are accepted, and what MF6 reads.

The ONE place the time-unit tables live (ledger 456), as
:mod:`myflopy._length_units` is for length. A time unit enters a model in two
places, and both are checked where they enter:

* a TDIS ``time_units`` option -- ``mf.tdis(time_units=)``, an edit of a TDIS spec
  (``with_options``/``packages.fork``/``packages.edit``), a PRT simulation's own
  TDIS and the legacy ``TemporalDiscretization`` -- which MF6 reads only as one of :data:`MF6_TIME_UNITS`, in any case. Any
  other word -- ``d``, ``day``, ``undefined`` -- STOPS the run ("Unrecognized input
  value for TIME_UNITS option"), where an unknown ``LENGTH_UNITS`` word is read as
  undefined (``docs/reproducers/tdis_time_units.py``, ledger 442). FloPy writes any
  word it is given, so it is refused at the declaration naming MF6's word;
* :class:`myflopy.specs.ModelContext` ``time_units`` -- the unit SFR/LAK's
  ``TIME_CONVERSION`` and the build warnings' flux labels read -- which takes any
  spelling in :data:`_TIME_SPELLINGS`, as the context's ``length_units`` takes any
  spelling a ``LayerStack`` does.

Every label abbreviates through :data:`MF6_TIME_UNITS`, so a TDIS read back from a
run and a context declared in a short spelling label the same unit the same way.
"""

from __future__ import annotations

from typing import Any

#: MF6's TDIS ``TIME_UNITS`` words, each with its abbreviation. MF6 reads these six in
#: any case and STOPS on any other word ("Unrecognized input value for TIME_UNITS
#: option") -- unlike ``LENGTH_UNITS``, which it reads as undefined (ledger 404).
#: ``unknown`` is MF6's own "undefined" and has no abbreviation
#: (``docs/reproducers/tdis_time_units.py``, ledger 442).
MF6_TIME_UNITS = {
    "unknown": None,
    "seconds": "s",
    "minutes": "min",
    "hours": "h",
    "days": "d",
    "years": "yr",
}

#: Every spelling of a time unit the library takes, mapped to MF6's word. A
#: ``ModelContext`` takes any of these; a TDIS takes MF6's words only, and a short
#: spelling's refusal names its word.
_TIME_SPELLINGS = {
    "seconds": "seconds", "second": "seconds", "sec": "seconds", "s": "seconds",
    "minutes": "minutes", "minute": "minutes", "min": "minutes",
    "hours": "hours", "hour": "hours", "hr": "hours", "h": "hours",
    "days": "days", "day": "days", "d": "days",
    "years": "years", "year": "years", "yr": "years",
}

#: MF6's ``TIME_CONVERSION`` (SFR, LAK outlets) per word: seconds -> model units. A
#: year is 365.25 days.
_SECONDS_PER = {
    "seconds": 1.0,
    "minutes": 60.0,
    "hours": 3600.0,
    "days": 86400.0,
    "years": 31557600.0,
}


def _spelling(units: Any) -> str | None:
    """The table key of a time unit: lower case, surrounding spaces stripped."""

    return units.strip().lower() if isinstance(units, str) else None


def _declares_no_unit(units: Any) -> bool:
    """``True`` for ``None``, ``""`` and MF6's ``"unknown"`` (any case, spaces ignored)."""

    return units is None or _spelling(units) in ("", "unknown")


def check_tdis_time_units(value: Any, *, helper: str) -> Any:
    """Refuse a TDIS ``time_units``, AS DECLARED, that MF6 does not read.

    The declaration rule of ``mf.tdis``, of an edit of a TDIS spec and of
    ``PRTProject(time_units=)`` (ledger 456): MF6's own word --
    ``unknown``, ``seconds``, ``minutes``, ``hours``, ``days`` or ``years``, in any
    case -- or ``None``, which writes no ``TIME_UNITS``. A short spelling is refused
    rather than mapped, naming the word to write, so one simulation has one spelling
    and one digest -- the rule ledger 404 set for a discretization's
    ``length_units``.

    Parameters
    ----------
    value : str or None
        The ``time_units`` option as passed. ``None`` is returned as is.
    helper : str
        What the caller wrote, for the message: ``"mf.tdis"``,
        ``"tdis.with_options"``, ``"PRTProject"`` or the legacy
        ``"TemporalDiscretization"``.

    Returns
    -------
    str or None
        ``value`` unchanged.

    Raises
    ------
    ValueError
        ``value`` is not one of MF6's six words in any case (a non-``str``
        included). MF6 stops a run on such a word; the message names MF6's words,
        and MF6's word for a known spelling (``time_units='days'`` for ``'d'``).
    """

    if value is None or (isinstance(value, str) and value.lower() in MF6_TIME_UNITS):
        return value
    hint = ""
    if _spelling(value) in _TIME_SPELLINGS:
        hint = f" Write time_units={_TIME_SPELLINGS[_spelling(value)]!r}."
    raise ValueError(
        f"{helper}(time_units={value!r}) is not a word MF6 reads: it takes "
        f"{list(MF6_TIME_UNITS)} (any case), and STOPS the run on any other word "
        f"('Unrecognized input value for TIME_UNITS option').{hint}"
    )


def mf6_time_word(units: Any) -> str | None:
    """MF6's word for a time unit in any spelling the library takes, or None for none.

    Parameters
    ----------
    units : str or None
        A time unit in any spelling of :data:`_TIME_SPELLINGS` (``"d"``,
        ``"day"``, ``"days"``, ...), any case, surrounding spaces ignored; or
        ``None`` / ``""`` / ``"unknown"`` for no unit.

    Returns
    -------
    str or None
        ``"seconds"``, ``"minutes"``, ``"hours"``, ``"days"`` or ``"years"``;
        ``None`` for no unit.

    Raises
    ------
    ValueError
        Any other value, a non-``str`` included, naming the spellings.
    """

    if _declares_no_unit(units):
        return None
    word = _TIME_SPELLINGS.get(_spelling(units))
    if word is None:
        raise ValueError(
            f"Unknown time unit {units!r}; known: {sorted(_TIME_SPELLINGS)} (any case), "
            "or None / '' / 'unknown' for no unit."
        )
    return word


def check_model_time_units(value: Any, *, owner: str = "ModelContext") -> Any:
    """Refuse a ``ModelContext(time_units=)`` or ``LayerStack(time_units=)`` nothing reads.

    The context's unit is read by ``mf.sfr``/``mf.lak`` for ``TIME_CONVERSION`` and
    by the build warnings' flux labels. ``mf.sfr``/``mf.lak`` refused an unknown
    name only on the branch that derived a conversion, so ``"weeks"`` was carried
    through every other model (ledger 456), as ``length_units`` was until
    ledger 433.

    Parameters
    ----------
    value : str or None
        Any spelling :func:`mf6_time_word` takes, or ``None`` / ``""`` /
        ``"unknown"`` for no declared unit (``TIME_CONVERSION`` 1.0).
    owner : str, optional
        The class the unit was given to, named in the message:
        ``"ModelContext"`` (the default) or ``"LayerStack"``.

    Returns
    -------
    str or None
        ``value`` unchanged.

    Raises
    ------
    ValueError
        Any other value, a non-``str`` included, naming the spellings.
    """

    try:
        mf6_time_word(value)
    except ValueError:
        raise ValueError(
            f"{owner}(time_units={value!r}) is not a time unit; known: "
            f"{sorted(_TIME_SPELLINGS)} (any case), or None / '' / 'unknown' for no unit."
        ) from None
    return value


def mf6_time_conversion(units: str | None) -> float:
    """Return MF6 ``TIME_CONVERSION`` (seconds -> model units) for a unit name.

    Parameters
    ----------
    units : str or None
        The model's time unit, in any spelling :func:`mf6_time_word` takes.
        ``None``, ``""`` or ``"unknown"`` means no declared unit.

    Returns
    -------
    float
        86400.0 for days, 3600.0 for hours, 60.0 for minutes, 1.0 for seconds,
        31557600.0 for years; 1.0 with no declared unit.

    Raises
    ------
    ValueError
        An unrecognized unit name.
    """

    word = mf6_time_word(units)
    return 1.0 if word is None else _SECONDS_PER[word]


def time_abbreviation(units: Any) -> str | None:
    """The label abbreviation of a time unit in any accepted spelling, or None.

    ``"d"`` for ``"days"``/``"day"``/``"d"``, and so on through
    :data:`MF6_TIME_UNITS`; ``None`` for no unit AND for a word that is not a time
    unit -- a label never raises.

    Parameters
    ----------
    units : str or None
        A time unit as a ``ModelContext`` or a TDIS holds it.

    Returns
    -------
    str or None
        ``"s"``, ``"min"``, ``"h"``, ``"d"`` or ``"yr"``, or ``None``.
    """

    word = _TIME_SPELLINGS.get(_spelling(units))
    return None if word is None else MF6_TIME_UNITS[word]
