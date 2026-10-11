"""A model's length unit: which spellings are accepted, and what MF6 reads.

The ONE place the length-unit tables live (ledger 432). They sat in
``surfaces.py`` (import layer 2), out of reach of ``specs.py`` (layer 1), so the two
places a length unit enters a model -- a DIS/DISV/DISU ``length_units`` option and
:class:`myflopy.specs.ModelContext` ``length_units`` -- could not ask them, and
``specs.mf6_length_conversion`` kept a third table of its own that disagreed with
both (``"metres"`` converted a surface and raised in ``mf.sfr``).

MF6 reads a DIS/DISV/DISU ``LENGTH_UNITS`` word other than ``FEET``/``METERS``/
``CENTIMETERS``/``UNKNOWN`` (any case) as unit 0, undefined, without a message, and
uses the unit for nothing but LABELS: ``lenuni`` is read by the listing line and the
NetCDF exporters alone (``Dis.f90``/``Disv.f90``/``Disu.f90``; measured: the
canonical model under ``FT`` and under ``feet`` writes 21 byte-identical output
files, ``docs/reproducers/mf6_length_units.py``). So a wrong word changes no number
MF6 computes; what it changes is every label read back from the model's own
declaration.
"""

from __future__ import annotations

from typing import Any

# Elevations are converted by a single scale factor (same datum, different unit), so
# only a per-unit "value in feet" table is needed.
_LENGTH_TO_FEET = {
    "feet": 1.0, "ft": 1.0, "foot": 1.0, "us-ft": 1.0, "usft": 1.0,
    "meters": 3.280839895013123, "meter": 3.280839895013123,
    "metres": 3.280839895013123, "m": 3.280839895013123,
    "centimeters": 0.03280839895013123, "centimeter": 0.03280839895013123,
    "cm": 0.03280839895013123,
}
# Model length unit -> MF6 LENGTH_UNITS token, for EVERY spelling `_LENGTH_TO_FEET`
# converts (ledger 404). It held only the three long names, and anything else was
# written upper-cased -- `"ft"` as `FT` -- which MF6 reads as length unit 0, undefined,
# under Normal termination (docs/reproducers/mf6_length_units.py). The two tables are
# held together by a test: a spelling the conversion accepts must have a token.
_MF6_LENGTH_UNITS = {
    "feet": "FEET", "ft": "FEET", "foot": "FEET", "us-ft": "FEET", "usft": "FEET",
    "meters": "METERS", "meter": "METERS", "metres": "METERS", "m": "METERS",
    "centimeters": "CENTIMETERS", "centimeter": "CENTIMETERS", "cm": "CENTIMETERS",
}
#: The words MF6's DIS/DISV/DISU ``LENGTH_UNITS`` reads (any case); every other word,
#: ``FT`` included, it reads as 0 -- undefined -- and says nothing.
_MF6_LENGTH_TOKENS = frozenset({"FEET", "METERS", "CENTIMETERS", "UNKNOWN"})

#: MF6's ``LENGTH_CONVERSION`` (SFR, LAK outlets) per token: metres -> model units.
_MF6_LENGTH_CONVERSION = {"FEET": 3.28081, "METERS": 1.0, "CENTIMETERS": 100.0}


def _is_unknown(units: Any) -> bool:
    """``True`` for MF6's ``UNKNOWN`` in any case and surrounding spaces."""

    return isinstance(units, str) and units.strip().upper() == "UNKNOWN"


def _declares_no_unit(units: Any) -> bool:
    """``True`` for the three ways a context or a call says "no length unit".

    ``None``, ``""`` and MF6's ``"unknown"`` (any case, surrounding spaces ignored)
    -- what :func:`mf6_length_conversion` reads as 1.0.
    """

    return units is None or (isinstance(units, str) and units.strip().lower() in ("", "unknown"))


def _spelling(units: Any) -> str | None:
    """The table key of a length unit: lower case, surrounding spaces stripped."""

    return units.strip().lower() if isinstance(units, str) else None


def _mf6_length_units(units: str) -> str:
    """The MF6 ``LENGTH_UNITS`` word for a model length unit, in any spelling the stack takes.

    The ONE mapping behind ``LayerSurfaces.to_disv`` and ``LayerBuildResult.to_disv``
    (so ``LayerStack.to_disv`` and ``mf.disv(stack)`` too), behind what the
    DIS/DISV/DISU builders write (:func:`mf6_length_option`) and behind the
    ``ModelContext(length_units=)`` check: every spelling :data:`_LENGTH_TO_FEET`
    converts maps to ``FEET``/``METERS``/``CENTIMETERS``.

    Parameters
    ----------
    units : str
        A length unit as a ``LayerStack(length_units=)`` or ``to_disv(length_units=)``
        takes it -- ``"feet"``, ``"ft"``, ``"foot"``, ``"us-ft"``, ``"usft"``,
        ``"meters"``, ``"meter"``, ``"metres"``, ``"m"``, ``"centimeters"``,
        ``"centimeter"``, ``"cm"`` -- in any case.

    Returns
    -------
    str
        ``"FEET"``, ``"METERS"`` or ``"CENTIMETERS"``.

    Raises
    ------
    ValueError
        ``units`` is not one of those spellings (a non-``str`` included), naming them.
        It used to be written upper-cased, and MF6 reads any word it does not know as
        length unit 0 without a message. MF6's own ``"unknown"`` is refused too, by a
        message naming ``length_units=None`` -- the stack's spelling of "no unit",
        which writes no ``LENGTH_UNITS`` and which MF6 reads as undefined just the same.
    """

    token = _MF6_LENGTH_UNITS.get(_spelling(units))
    if token is None and _is_unknown(units):
        # MF6 KNOWS this word, so the generic message below would be false for it.
        # A stack has no factor for it; "no unit" is spelled `None`, which writes no
        # LENGTH_UNITS -- what MF6 reads as UNKNOWN (ledger 404).
        raise ValueError(
            f"length_units={units!r} is not a unit a stack can convert into. To declare "
            "no length unit, pass length_units=None: nothing is converted and no "
            "LENGTH_UNITS is written, which MF6 reads as undefined (unit 0) exactly as "
            f"it reads UNKNOWN. Known units: {sorted(_MF6_LENGTH_UNITS)}."
        )
    if token is None:
        raise ValueError(
            f"Unknown length unit {units!r}; known: {sorted(_MF6_LENGTH_UNITS)}. MF6 "
            "reads a LENGTH_UNITS word it does not know as 0 (undefined) and says "
            "nothing, so it is refused here rather than written."
        )
    return token


def _check_mf6_length_token(value: Any, *, helper: str) -> Any:
    """Refuse a DIS/DISV/DISU ``length_units`` option, AS DECLARED, that is not MF6's word.

    The declaration rule (ledger 404, extended to every edit route by ledger
    432): ``mf.dis``/``mf.disv``/``mf.disu`` and ``with_options``/
    ``packages.fork``/``packages.edit`` on a discretization take MF6's own word --
    ``feet``, ``meters``, ``centimeters`` or ``unknown`` in any case -- and refuse
    anything else naming the word to write. An abbreviation is refused rather than
    mapped at the declaration so one model has one spelling, and one digest; what
    an EXISTING spec holds is mapped at build (:func:`mf6_length_option`).

    Parameters
    ----------
    value : str or None
        The ``length_units`` option as passed. ``None`` writes nothing and is
        returned as is.
    helper : str
        What the caller wrote, for the message: ``"mf.disv"``, or
        ``"disv.with_options"``.

    Returns
    -------
    str or None
        ``value`` unchanged.

    Raises
    ------
    ValueError
        ``value`` is not one of MF6's four words; an abbreviation's message names its
        word (``length_units='feet'`` for ``'ft'``).
    """

    if value is None or (isinstance(value, str) and value.upper() in _MF6_LENGTH_TOKENS):
        return value
    hint = ""
    if _spelling(value) in _MF6_LENGTH_UNITS:
        hint = f" Write length_units={_MF6_LENGTH_UNITS[_spelling(value)].lower()!r}."
    raise ValueError(
        f"{helper}(length_units={value!r}) is not a word MF6 reads: it takes "
        f"{sorted(token.lower() for token in _MF6_LENGTH_TOKENS)} (any case), and reads "
        f"every other word as length unit 0, undefined, without a message.{hint}"
    )


def mf6_length_option(value: Any) -> Any:
    """What a DIS/DISV/DISU builder WRITES for the ``length_units`` its spec holds.

    The build-time half of the rule (ledger 432). A spec built through a
    declaration route already holds MF6's word; one that was not -- a frozen run
    snapshot, a saved project or a spec edited before the declaration check reached
    ``with_options``/``packages.fork``/``packages.edit`` -- may hold an abbreviation
    such as ``"ft"``, which FloPy wrote as ``FT`` and MF6 read as unit 0. It is
    MAPPED here rather than refused, so such a run still rebuilds, now declaring
    the unit it meant. MF6 computes nothing from the unit (the module docstring),
    so a rebuilt run's heads and budgets are byte-identical; its listing and every
    label read from the model's declaration now name the unit.

    Parameters
    ----------
    value : str or None
        The spec's ``length_units`` option. ``None`` writes nothing.

    Returns
    -------
    str or None
        ``None`` for ``None``; ``value`` unchanged when it is already one of MF6's
        four words in any case (so a build writes exactly what was declared); else
        the MF6 word of an accepted spelling (``"ft"`` -> ``"FEET"``).

    Raises
    ------
    ValueError
        A word that is neither MF6's nor a spelling :data:`_LENGTH_TO_FEET` takes
        (a non-``str`` included), naming the spellings.
    """

    if value is None or (isinstance(value, str) and value.upper() in _MF6_LENGTH_TOKENS):
        return value
    if _spelling(value) in _MF6_LENGTH_UNITS:
        return _MF6_LENGTH_UNITS[_spelling(value)]
    raise ValueError(
        f"length_units={value!r} is not a length unit: MF6 reads "
        f"{sorted(token.lower() for token in _MF6_LENGTH_TOKENS)} (any case) and these "
        f"spellings map to them: {sorted(_MF6_LENGTH_UNITS)}. Any other word is read "
        "as length unit 0, undefined, without a message."
    )


def check_model_length_units(value: Any) -> Any:
    """Refuse a ``ModelContext(length_units=)`` no consumer of it can read.

    The context's unit is read by ``mf.sfr``/``mf.lak`` for ``LENGTH_CONVERSION``,
    by ``mf.lak`` to convert a lakebed surface's ``units=``, and for labels. Each
    consumer refused an unknown name only on the branch that used it, and a
    ``Surface`` with no ``units=`` never asked -- so ``"inches"`` was carried
    through a whole model (ledger 433).

    Parameters
    ----------
    value : str or None
        ``None``, ``""`` or MF6's ``"unknown"`` -- no declared unit:
        ``LENGTH_CONVERSION`` is 1.0, and a lakebed surface declaring ``units=``
        is refused -- or any spelling :data:`_LENGTH_TO_FEET` takes. Case and
        surrounding spaces are ignored.

    Returns
    -------
    str or None
        ``value`` unchanged.

    Raises
    ------
    ValueError
        Any other value, a non-``str`` included, naming the spellings.
    """

    if _declares_no_unit(value) or _spelling(value) in _MF6_LENGTH_UNITS:
        return value
    raise ValueError(
        f"ModelContext(length_units={value!r}) is not a length unit; known: "
        f"{sorted(_MF6_LENGTH_UNITS)} (any case), or None / '' / 'unknown' for no unit."
    )


def mf6_length_conversion(units: str | None) -> float:
    """Return MF6 ``LENGTH_CONVERSION`` (meters -> model units) for a unit name.

    Parameters
    ----------
    units : str or None
        The model's length unit, in any spelling :data:`_LENGTH_TO_FEET` takes and
        any case, surrounding spaces ignored. ``None``, ``""`` or ``"unknown"``
        means no declared unit.

    Returns
    -------
    float
        3.28081 for feet, 1.0 for metres, 100.0 for centimetres; 1.0 with no
        declared unit.

    Raises
    ------
    ValueError
        An unrecognized unit name.
    """

    if _declares_no_unit(units):
        return 1.0
    token = _MF6_LENGTH_UNITS.get(_spelling(units))
    if token is None:
        raise ValueError(
            f"Unknown length unit {units!r}; known: {sorted(_MF6_LENGTH_UNITS)}."
        )
    return _MF6_LENGTH_CONVERSION[token]


def _length_factor(from_units: str, to_units: str) -> float:
    """Multiplier converting an elevation from ``from_units`` to ``to_units``."""

    try:
        return _LENGTH_TO_FEET[from_units.strip().lower()] / _LENGTH_TO_FEET[to_units.strip().lower()]
    except KeyError as exc:
        raise ValueError(
            f"Unknown length unit {exc.args[0]!r}; known: {sorted(_LENGTH_TO_FEET)}."
        ) from None
