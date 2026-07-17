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

# Redundant alias: private but part of the package-root contract — the pest
# package reuses pyEMU's row-labeling rule through it (plan 4.1 note).
from myflopy.modflow.mf6.observations._shared import (
    _normalize_row_labels as _normalize_row_labels,
)
from myflopy.modflow.mf6.observations.drn import BoundDrnFlowTargets, DrnFlowTargets
from myflopy.modflow.mf6.observations.heads import (
    BoundHeadTargetPlots,
    BoundHeadTargets,
    HeadTargets,
)
from myflopy.modflow.mf6.observations.lake import BoundLakeStageTargets, LakeStageTargets
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.observations.registry import TargetRegistry
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
    "BoundSfrFlowTargets",
    "BoundSfrStageTargets",
    "DrnFlowTargets",
    "HeadTargets",
    "LakeStageTargets",
    "SfrFlowTargets",
    "SfrStageTargets",
    "TargetRegistry",
]
