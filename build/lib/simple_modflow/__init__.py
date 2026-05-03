"""Public package surface for simple_modflow."""

from __future__ import annotations

from importlib import import_module

try:
    from importlib.metadata import version

    __version__ = version("simple_modflow")
except Exception:
    __version__ = "0.1.0"


_EXPORTS = {
    "SimpleModel": ("simple_modflow.modflow.mf6.simplemodel", "SimpleModel"),
    "SimulationBase": ("simple_modflow.modflow.mf6.mfsimbase", "SimulationBase"),
    "TriangleGrid": ("simple_modflow.modflow.mf6.grid.triangle", "TriangleGrid"),
    "VoronoiGridPlus": ("simple_modflow.modflow.mf6.voronoiplus", "VoronoiGridPlus"),
    "modflow": ("simple_modflow", "modflow"),
    "read_gpkg": ("simple_modflow.modflow.utils.datatypes.readers", "read_gpkg"),
    "read_shp_gpkg": ("simple_modflow.modflow.utils.datatypes.readers", "read_shp_gpkg"),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    if name == "modflow":
        return import_module("simple_modflow.modflow")

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    return sorted(list(globals().keys()) + __all__)
