"""Reusable observation and target data helpers for MF6 workflows.

This package intentionally sits outside the calibration-specific ``pest``
package. Observation targets such as groundwater heads are useful for several
tasks in ``myflopy``:

- PEST/pyEMU observation setup
- calibration residual review
- plotting and exploratory statistics
- comparing simulated heads against measured targets without a PEST workflow

Split from the former single-module ``observations.py`` (implementation plan
4.1); this ``__init__`` re-exports the full public surface so both import
paths keep working: ``from myflopy.modflow.mf6.observations import
HeadTargets`` and ``from myflopy.modflow.mf6.observations.heads import
HeadTargets``.
"""

from __future__ import annotations

from myflopy.modflow.mf6.observations.conc import (
    ConcDifferences,
    ConcTargets,
)
from myflopy.modflow.mf6.observations.drn import BoundDrnFlowTargets, DrnFlowTargets
from myflopy.modflow.mf6.observations.heads import (
    BoundHeadTargetPlots,
    BoundHeadTargets,
    HeadDifferences,
    HeadTargets,
)
from myflopy.modflow.mf6.observations.lake import BoundLakeStageTargets, LakeStageTargets
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.observations.registry import TargetRegistry
from myflopy.modflow.mf6.observations.riv import BoundRivFlowTargets, RivFlowTargets
from myflopy.modflow.mf6.observations.sfr import (
    BoundSfrFlowTargets,
    BoundSfrStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)

__all__ = [
    "BoundDrnFlowTargets",
    "BoundHeadTargetPlots",
    "BoundHeadTargets",
    "BoundLakeStageTargets",
    "BoundNamedSeriesTargetPlots",
    "BoundRivFlowTargets",
    "BoundSfrFlowTargets",
    "BoundSfrStageTargets",
    "ConcDifferences",
    "ConcTargets",
    "DrnFlowTargets",
    "HeadDifferences",
    "HeadTargets",
    "LakeStageTargets",
    "RivFlowTargets",
    "SfrFlowTargets",
    "SfrStageTargets",
    "TargetRegistry",
]
