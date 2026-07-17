"""Group-oriented helpers for comparing multiple models with one lazy API.

Split from the former single-module ``project/model_group.py`` (implementation
plan 4.2). ``myflopy.project.model_group`` remains as a facade over this
package, so both import paths keep working.
"""

from __future__ import annotations

from myflopy.project.group.budget import (
    GroupBudget,
)
from myflopy.project.group.core import (
    ModelGroup,
)
from myflopy.project.group.inputs import (
    GroupPackageInputField,
    GroupPackageInputs,
)
from myflopy.project.group.lak import (
    GroupLakBudgetResults,
    GroupLakConnections,
    GroupLakOutputs,
    GroupLakPackageAccessor,
    GroupLakResultsNamespace,
    GroupLakStageResults,
)
from myflopy.project.group.packages import (
    GroupOutputs,
    GroupPackageAccessor,
    GroupPackages,
    GroupResultsOnlyPackageAccessor,
)
from myflopy.project.group.results import (
    GroupCellPackageResults,
    GroupCellPackageResultsNamespace,
)
from myflopy.project.group.sfr import (
    GroupSfrBudgetResults,
    GroupSfrResultsNamespace,
    GroupSfrStageResults,
)
from myflopy.project.group.spatial import (
    GroupHeads,
)
from myflopy.project.group.surface_water import (
    GroupSurfaceWaterExchangeResults,
    GroupSurfaceWaterResultsNamespace,
)
from myflopy.project.group.uzf import (
    GroupUzfFieldAccessor,
    GroupUzfInputs,
    GroupUzfPackageAccessor,
    GroupUzfResultsNamespace,
)

__all__ = [
    "GroupBudget",
    "GroupCellPackageResults",
    "GroupCellPackageResultsNamespace",
    "GroupHeads",
    "GroupLakBudgetResults",
    "GroupLakConnections",
    "GroupLakOutputs",
    "GroupLakPackageAccessor",
    "GroupLakResultsNamespace",
    "GroupLakStageResults",
    "GroupOutputs",
    "GroupPackageAccessor",
    "GroupPackageInputField",
    "GroupPackageInputs",
    "GroupPackages",
    "GroupResultsOnlyPackageAccessor",
    "GroupSfrBudgetResults",
    "GroupSfrResultsNamespace",
    "GroupSfrStageResults",
    "GroupSurfaceWaterExchangeResults",
    "GroupSurfaceWaterResultsNamespace",
    "GroupUzfFieldAccessor",
    "GroupUzfInputs",
    "GroupUzfPackageAccessor",
    "GroupUzfResultsNamespace",
    "ModelGroup",
]
