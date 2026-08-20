"""Import helpers for optional third-party dependencies.

Some myflopy features depend on packages that are not core requirements: PEST
calibration needs :mod:`pyemu`, the unstructured xarray/NetCDF interchange
needs :mod:`xugrid`, and interactive 3-D scenes need :mod:`pyvista`/:mod:`trame`. These are declared as install *extras* in
``pyproject.toml`` and imported lazily, so the base package imports cleanly
without them. This module centralizes the "import it or explain how to install
it" pattern behind :func:`require`.
"""

from __future__ import annotations

import importlib
from types import ModuleType

# Maps an importable module name to the pip extra that provides it.
_EXTRA_FOR_MODULE = {
    "pyemu": "pest",
    "xugrid": "xugrid",
    "xarray": "xugrid",
    "pyvista": "viz3d",
    "trame": "viz3d",
    "dash": "viz",
}


def require(module: str, *, feature: str | None = None) -> ModuleType:
    """Import and return ``module``, or raise a clear install hint if missing.

    Parameters
    ----------
    module
        The importable module name (e.g. ``"pyemu"``).
    feature
        Optional short description of the feature that needs it, woven into the
        error message (e.g. ``"PEST calibration workflows"``).

    Returns
    -------
    module
        The imported module object.

    Raises
    ------
    ModuleNotFoundError
        If the module is not installed, with a ``pip install "myflopy[extra]"``
        hint when the module maps to a known extra.
    """

    try:
        return importlib.import_module(module)
    except ModuleNotFoundError as exc:  # pragma: no cover - environment guard
        extra = _EXTRA_FOR_MODULE.get(module)
        hint = f' Install it with:  pip install "myflopy[{extra}]"' if extra else ""
        what = f" for {feature}" if feature else ""
        raise ModuleNotFoundError(
            f"The optional dependency {module!r} is required{what}.{hint}"
        ) from exc
