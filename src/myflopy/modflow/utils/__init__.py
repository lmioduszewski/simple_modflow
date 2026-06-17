"""Convenience exports for utility helpers."""

from __future__ import annotations

from importlib import import_module


_EXPORTS = {
    "geotiff_to_contours": ("myflopy.modflow.utils.raster", "geotiff_to_contours"),
    "get_iheads": ("myflopy.modflow.utils.iheads", "get_iheads"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    return sorted(list(globals().keys()) + __all__)
