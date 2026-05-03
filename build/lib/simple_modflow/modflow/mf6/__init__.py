"""Public MF6-facing exports for simple_modflow."""

from __future__ import annotations

from importlib import import_module


_SUBMODULES = {
    "grid": "simple_modflow.modflow.mf6.grid",
    "simulation": "simple_modflow.modflow.mf6.simulation",
}

_EXPORTS = {
    "Boundaries": ("simple_modflow.modflow.mf6.boundaries", "Boundaries"),
    "CHD": ("simple_modflow.modflow.mf6.simulation.packages", "CHD"),
    "DRN": ("simple_modflow.modflow.mf6.drn", "DRN"),
    "DisuGrid": ("simple_modflow.modflow.mf6.simulation.discretization", "DisuGrid"),
    "DisvGrid": ("simple_modflow.modflow.mf6.simulation.discretization", "DisvGrid"),
    "Drains": ("simple_modflow.modflow.mf6.simulation.packages", "Drains"),
    "GHB": ("simple_modflow.modflow.mf6.ghb", "GHB"),
    "InitialConditions": ("simple_modflow.modflow.mf6.simulation.packages", "InitialConditions"),
    "KFlow": ("simple_modflow.modflow.mf6.simulation.packages", "KFlow"),
    "LAK": ("simple_modflow.modflow.mf6.simulation.packages", "LAK"),
    "LakeAreaVolumeRelationship": ("simple_modflow.modflow.mf6.lakes", "LakeAreaVolumeRelationship"),
    "OutputControl": ("simple_modflow.modflow.mf6.simulation.packages", "OutputControl"),
    "Recharge": ("simple_modflow.modflow.mf6.simulation.packages", "Recharge"),
    "RechargeFromShp": ("simple_modflow.modflow.mf6.recharge", "RechargeFromShp"),
    "SFR": ("simple_modflow.modflow.mf6.sfr", "SFR"),
    "SimpleModel": ("simple_modflow.modflow.mf6.simplemodel", "SimpleModel"),
    "SimulationBase": ("simple_modflow.modflow.mf6.mfsimbase", "SimulationBase"),
    "Storage": ("simple_modflow.modflow.mf6.simulation.packages", "Storage"),
    "TemporalDiscretization": ("simple_modflow.modflow.mf6.simulation.discretization", "TemporalDiscretization"),
    "TriangleGrid": ("simple_modflow.modflow.mf6.grid.triangle", "TriangleGrid"),
    "UZF": ("simple_modflow.modflow.mf6.simulation.packages", "UZF"),
    "UZFPackageData": ("simple_modflow.modflow.mf6.uzf", "UZFPackageData"),
    "VoronoiGridPlus": ("simple_modflow.modflow.mf6.voronoiplus", "VoronoiGridPlus"),
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
