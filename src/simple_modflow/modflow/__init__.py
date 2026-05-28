"""Convenience exports for the modflow package namespace."""

from __future__ import annotations

from importlib import import_module


_SUBMODULES = {
    "calcs": "simple_modflow.modflow.calcs",
    "gwt": "simple_modflow.modflow.gwt",
    "mf6": "simple_modflow.modflow.mf6",
    "mp3du": "simple_modflow.modflow.mp3du",
    "utils": "simple_modflow.modflow.utils",
}

_EXPORTS = {
    "get_iheads": ("simple_modflow.modflow.utils.iheads", "get_iheads"),
    "geotiff_to_contours": ("simple_modflow.modflow.utils.raster", "geotiff_to_contours"),
    "ParticleTrackingInput": ("simple_modflow.modflow.mp3du", "ParticleTrackingInput"),
    "prepare_particle_tracking": ("simple_modflow.modflow.mp3du", "prepare_particle_tracking"),
    "run_particle_tracking": ("simple_modflow.modflow.mp3du", "run_particle_tracking"),
}

__all__ = sorted([*_SUBMODULES.keys(), *_EXPORTS.keys()])


def __getattr__(name: str):
    if name in _SUBMODULES:
        return import_module(_SUBMODULES[name])

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    return sorted(list(globals().keys()) + __all__)
