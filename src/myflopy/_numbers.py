"""The one check for a length, step or size that must be a finite number above zero.

Three copies of it sat in the grid modules (ledger 457), and they disagreed:
``_lattice_axis`` cast first, so ``"2.5"`` and ``True`` passed as a spacing of 2.5
and 1.0, where ``generate_grid_coordinates`` and the strike/dip raster refused both.
The same cast let a lattice ``size`` of ``"2"`` or ``True`` through, and the box
``gap`` held a fourth copy; both read this rule with ``allow_zero=True``.

Also the one rule for MODFLOW's NO-VALUE numbers (:func:`is_no_value`, :func:`is_dry`).
It lived in ``modflow/mf6/headsplus.py`` as ``_is_sentinel`` -- which re-exports it
still, for the eleven modules importing it from there -- until a grid module needed it
too: ``modflow/mf6/grid/remap.py`` carries a solved field between meshes, the heads
reader imports that module, and importing the rule back from the reader would be a
cycle. Down here, at the bottom of the import layering, every module can reach it.
"""

from __future__ import annotations

import numbers

import numpy as np


def positive_number(
    value, name: str, *, units: str | None = None, allow_zero: bool = False
) -> float:
    """``value`` as a ``float``, or a ``ValueError`` naming ``name``.

    Parameters
    ----------
    value : int or float
        A real number -- a Python or numpy integer or float -- finite and above
        zero (zero or above with ``allow_zero=True``). A ``bool`` (a Python
        ``int`` subclass), a string, ``NaN``, an infinity and a negative are
        refused, and so is zero unless ``allow_zero``.
    name : str
        The argument's name, for the message (``"spacing"``).
    units : str, optional
        What the number is measured in, for the message (``"model length
        units"``). Default ``None``: not named.
    allow_zero : bool, optional
        ``True`` accepts ``0`` as well -- a lattice ``size`` or a box ``gap``, which
        may be nothing but not less. The message then says "zero or more" rather
        than "above zero". Default ``False``.

    Returns
    -------
    float
        ``float(value)``.

    Raises
    ------
    ValueError
        ``value`` is anything else, naming ``name`` and what was given.
    """

    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, numbers.Real)
        or not np.isfinite(value)
        or value < 0
        or (value == 0 and not allow_zero)
    ):
        of = f" of {units}" if units else ""
        bound = ", zero or more," if allow_zero else " above zero,"
        raise ValueError(f"{name} must be a finite number{of}{bound} got {value!r}.")
    return float(value)


#: MODFLOW's no-value threshold: a number at or beyond it in magnitude is a marker,
#: not a value. MF6 writes +1e30 for an INACTIVE cell and -1e30 for a DRY one, and no
#: real head, concentration or temperature comes near either.
NO_VALUE = 1.0e29


def is_no_value(values) -> np.ndarray:
    """True where ``values`` is one of MODFLOW's no-value markers: ``|value| >= 1e29``.

    MF6 writes a cell with no value as a number, not NaN: +1e30 for an inactive
    cell and -1e30 for a DRY one. No real head, concentration or temperature
    comes near either, so the test is the magnitude, and it has to take BOTH
    signs. A head surface cut at ``> 10_000`` instead, which kept every dry cell
    as a crater to -1e30 and erased every real head above 10,000 ft (ledger 315).

    Parameters
    ----------
    values : array_like
        Numbers as a MODFLOW output file holds them; anything ``numpy.asarray``
        turns into floats.

    Returns
    -------
    numpy.ndarray
        A boolean array of the same shape, True at each marker.
    """

    return np.abs(np.asarray(values, dtype=float)) >= NO_VALUE


def is_dry(values) -> np.ndarray:
    """True where ``values`` is MF6's DRY marker: ``value <= -1e29``.

    The negative half of :func:`is_no_value`. MF6 writes ``-1e30`` (HDRY) for a
    cell that went dry under the standard formulation and ``+1e30`` for an
    inactive one; a reader that treats the two differently -- a dry cell's water
    table is in the layer below it, an inactive cell has nothing -- tells them
    apart here. Under NEWTON MF6 writes a dry cell's computed head, below its
    bottom, as an ordinary number, so this is never true there.

    Parameters
    ----------
    values : array_like
        Numbers as a MODFLOW output file holds them.

    Returns
    -------
    numpy.ndarray
        A boolean array of the same shape, True at each dry marker.
    """

    return np.asarray(values, dtype=float) <= -NO_VALUE
