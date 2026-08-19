"""Simulation helpers extracted from the MF6 base model module."""

from .accessors import (
    get_all_heads,
    get_budget,
    get_budget_cumulative,
    get_budget_incremental,
    get_hds,
    get_inputs,
    get_kstpkper,
    get_lak_output,
    get_sfr_output,
    get_uzf_output,
)
from .base import SimulationBase
from .discretization import DisuGrid, DisvGrid, TemporalDiscretization
from .indexing import (
    build_idomain,
    build_model_times,
    build_ncpl_arr,
    build_node_to_lni,
    build_offsets,
    coerce_per_dates,
)
from .packages import (
    CHD,
    GHB,
    UZF,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
    Wells,
)
from .regions import ModelRegion, RegionGroup, RegionRegistry
from .runtime import run_simulation

__all__ = [
    "CHD",
    "DisuGrid",
    "DisvGrid",
    "Drains",
    "GHB",
    "InitialConditions",
    "KFlow",
    "ModelRegion",
    "RegionGroup",
    "OutputControl",
    "Recharge",
    "Wells",
    "RegionRegistry",
    "SimulationBase",
    "Storage",
    "TemporalDiscretization",
    "UZF",
    "build_idomain",
    "build_model_times",
    "build_ncpl_arr",
    "build_node_to_lni",
    "build_offsets",
    "coerce_per_dates",
    "get_all_heads",
    "get_budget",
    "get_budget_cumulative",
    "get_budget_incremental",
    "get_hds",
    "get_inputs",
    "get_kstpkper",
    "get_lak_output",
    "get_sfr_output",
    "get_uzf_output",
    "run_simulation",
]
