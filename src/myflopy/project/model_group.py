"""Facade over :mod:`myflopy.project.group` (split there by plan 4.2).

Every name the old single-module ``model_group`` defined stays importable
from this path — including the private helpers and views that tests and
downstream code reach for. New code should import from
``myflopy.project.group`` (or its submodules) directly.
"""

from __future__ import annotations

from myflopy.project.group._shared import (
    TResultsNamespace,
    _coerce_kstpkper,
    _coerce_models,
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _grid_signature_for_model,
    _load_group_model,
    _normalize_iterable_filter,
    _reduce_to_period_end,
    _resolve_group_compare_target,
    _stable_compare_keys,
)
from myflopy.project.group.budget import (
    GroupBudget,
)
from myflopy.project.group.conc import (
    GroupConc,
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
    _GroupFieldView,
    _GroupSpatialView,
)
from myflopy.project.group.surface_water import (
    GroupSurfaceWaterExchangeResults,
    GroupSurfaceWaterResultsNamespace,
)
from myflopy.project.group.temp import (
    GroupTemp,
)
from myflopy.project.group.uzf import (
    GroupUzfFieldAccessor,
    GroupUzfInputs,
    GroupUzfPackageAccessor,
    GroupUzfResultsNamespace,
)

__all__ = ['GroupBudget', 'GroupCellPackageResults', 'GroupCellPackageResultsNamespace', 'GroupConc', 'GroupHeads', 'GroupLakBudgetResults', 'GroupLakConnections', 'GroupLakOutputs', 'GroupLakPackageAccessor', 'GroupLakResultsNamespace', 'GroupLakStageResults', 'GroupOutputs', 'GroupPackageAccessor', 'GroupPackageInputField', 'GroupPackageInputs', 'GroupPackages', 'GroupResultsOnlyPackageAccessor', 'GroupSfrBudgetResults', 'GroupSfrResultsNamespace', 'GroupSfrStageResults', 'GroupSurfaceWaterExchangeResults', 'GroupSurfaceWaterResultsNamespace', 'GroupTemp', 'GroupUzfFieldAccessor', 'GroupUzfInputs', 'GroupUzfPackageAccessor', 'GroupUzfResultsNamespace', 'ModelGroup', 'TResultsNamespace', '_GroupFieldView', '_GroupSpatialView', '_coerce_kstpkper', '_coerce_models', '_ensure_group_map_compatible', '_filter_group_input_table', '_grid_signature_for_model', '_load_group_model', '_normalize_iterable_filter', '_reduce_to_period_end', '_resolve_group_compare_target', '_stable_compare_keys']
