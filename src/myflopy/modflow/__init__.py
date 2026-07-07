"""Convenience exports for the modflow package namespace."""

from __future__ import annotations

from importlib import import_module


_SUBMODULES = {
    "calcs": "myflopy.modflow.calcs",
    "gwt": "myflopy.modflow.gwt",
    "mf6": "myflopy.modflow.mf6",
    "mp3du": "myflopy.modflow.mp3du",
    "utils": "myflopy.modflow.utils",
}

_EXPORTS = {
    "get_iheads": ("myflopy.modflow.utils.iheads", "get_iheads"),
    "geotiff_to_contours": ("myflopy.modflow.utils.raster", "geotiff_to_contours"),
    "ParticleTrackingInput": ("myflopy.modflow.mp3du", "ParticleTrackingInput"),
    "prepare_particle_tracking": ("myflopy.modflow.mp3du", "prepare_particle_tracking"),
    "run_particle_tracking": ("myflopy.modflow.mp3du", "run_particle_tracking"),
}

__all__ = sorted([*_SUBMODULES.keys(), *_EXPORTS.keys()])


def __getattr__(name: str):
    """Lazily import a submodule or public export on first attribute access."""

    if name in _SUBMODULES:
        return import_module(_SUBMODULES[name])

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    """Advertise the lazily-exported names for tab-completion."""

    return sorted(list(globals().keys()) + __all__)
