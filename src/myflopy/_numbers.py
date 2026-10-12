"""The one check for a length, step or size that must be a finite number above zero.

Three copies of it sat in the grid modules (ledger 457), and they disagreed:
``_lattice_axis`` cast first, so ``"2.5"`` and ``True`` passed as a spacing of 2.5
and 1.0, where ``generate_grid_coordinates`` and the strike/dip raster refused both.
"""

from __future__ import annotations

import numbers

import numpy as np


def positive_number(value, name: str, *, units: str | None = None) -> float:
    """``value`` as a ``float``, or a ``ValueError`` naming ``name``.

    Parameters
    ----------
    value : int or float
        A real number -- a Python or numpy integer or float -- finite and above
        zero. A ``bool`` (a Python ``int`` subclass), a string, ``NaN``, an
        infinity, zero and a negative are refused.
    name : str
        The argument's name, for the message (``"spacing"``).
    units : str, optional
        What the number is measured in, for the message (``"model length
        units"``). Default ``None``: not named.

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
        or value <= 0
    ):
        of = f" of {units}" if units else ""
        raise ValueError(f"{name} must be a finite number{of} above zero, got {value!r}.")
    return float(value)
