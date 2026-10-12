"""The one check for a length, step or size that must be a finite number above zero.

Three copies of it sat in the grid modules (ledger 457), and they disagreed:
``_lattice_axis`` cast first, so ``"2.5"`` and ``True`` passed as a spacing of 2.5
and 1.0, where ``generate_grid_coordinates`` and the strike/dip raster refused both.
The same cast let a lattice ``size`` of ``"2"`` or ``True`` through, and the box
``gap`` held a fourth copy; both read this rule with ``allow_zero=True``.
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
