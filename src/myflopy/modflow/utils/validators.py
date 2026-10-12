from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase



def valid_package_names_3ltr(model: SimulationBase) -> list | None:
    """
    Returns a list of valid package names for a modflow model.
    :param model: model object
    :return: list
    """
    # Duck-typed. `SimulationBase` is imported for annotations only -- it sits in a
    # higher layer -- so the `isinstance` this used to be raised NameError on every
    # call, taking `package_name_validator` and the legacy datalists with it.
    types = getattr(model, "package_types", None)
    if types is None:
        return None
    # TYPES, not names: every caller compares lowercase three-letter types, and a
    # model whose CHD is named `chd_0` has no package NAMED `chd` (ledger 228).
    return [str(name).lower() for name in types]


def package_name_validator(model: SimulationBase, package: str):
    """
    Checks that a package name is valid for a given modflow model.
    :param model: modflow model object
    :param package: str of package name, 3 letters (ex. 'drn', 'ghb', 'lak')
    :return: returns package name if valid, otherwise raises error
    """
    valid_names = valid_package_names_3ltr(model)
    assert package in valid_names, \
        f'Package not valid, must be one of {valid_names}'
    return package


def valid_list_of_cell_ints(model: SimulationBase = None, cells: list = None):
    """Assert ``cells`` is a list of ints, each a valid grid cell when ``model`` is given."""

    assert all(isinstance(loc, int) for loc in cells), 'cells must be a list of ints'
    if model:
        valid_cells = model.vor.gdf_vorPolys.index.to_list()
        assert all(loc in valid_cells for loc in cells), 'Not valid cells'

    return True


def whole_number(value, *, what: str) -> int:
    """``value`` as a Python ``int`` when it IS a whole number; anything else raises.

    The one rule for an index a caller types -- a cell number, a layer (ledger
    478). ``int()`` was the rule before, and it TRUNCATES a float (``5.7`` read
    as cell 5, ``layer=1.7`` as layer 1) and reads a numeric string (``"5"``), both in
    silence. A Python or numpy integer is taken; so is an integral float (``5.0``, what
    a pandas column holding a NaN makes of its integers). A ``bool``, a fractional or
    non-finite float, a string and anything else raise ``TypeError`` naming ``what``.

    Parameters
    ----------
    value : int, numpy.integer or float
        The index as given.
    what : str
        What the index is, for the message: ``"cell"``, ``"layer"``.

    Returns
    -------
    int
        The index.

    Raises
    ------
    TypeError
        When ``value`` is not a whole number.

    Examples
    --------
    >>> whole_number(5.0, what="cell")
    5
    >>> whole_number(5.7, what="cell")
    Traceback (most recent call last):
    ...
    TypeError: A cell is a whole number, not 5.7 (float).
    """

    import math
    from numbers import Integral, Real

    if isinstance(value, Integral) and not isinstance(value, bool):
        return int(value)
    if (
        isinstance(value, Real) and not isinstance(value, bool)
        and math.isfinite(float(value)) and float(value).is_integer()
    ):
        return int(value)
    raise TypeError(f"A {what} is a whole number, not {value!r} ({type(value).__name__}).")
