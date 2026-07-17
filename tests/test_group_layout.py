"""The project/group/ package layout (implementation plan 4.2).

The historical facade (``myflopy.project.model_group``), the package root,
and the real submodules must all resolve to the same objects.
"""

from __future__ import annotations

import importlib

import myflopy.project.group as group_pkg
import myflopy.project.model_group as facade

EXPECTED_HOMES = {
    "GroupHeads": "spatial",
    "GroupBudget": "budget",
    "GroupPackageInputField": "inputs",
    "GroupPackageInputs": "inputs",
    "GroupCellPackageResults": "results",
    "GroupCellPackageResultsNamespace": "results",
    "GroupSfrBudgetResults": "sfr",
    "GroupSfrStageResults": "sfr",
    "GroupSfrResultsNamespace": "sfr",
    "GroupLakBudgetResults": "lak",
    "GroupLakOutputs": "lak",
    "GroupLakStageResults": "lak",
    "GroupLakConnections": "lak",
    "GroupLakResultsNamespace": "lak",
    "GroupLakPackageAccessor": "lak",
    "GroupUzfFieldAccessor": "uzf",
    "GroupUzfInputs": "uzf",
    "GroupUzfPackageAccessor": "uzf",
    "GroupUzfResultsNamespace": "uzf",
    "GroupOutputs": "packages",
    "GroupPackageAccessor": "packages",
    "GroupResultsOnlyPackageAccessor": "packages",
    "GroupPackages": "packages",
    "GroupSurfaceWaterExchangeResults": "surface_water",
    "GroupSurfaceWaterResultsNamespace": "surface_water",
    "ModelGroup": "core",
}


def test_package_root_exports_all_public_classes():
    assert sorted(group_pkg.__all__) == sorted(EXPECTED_HOMES)


def test_facade_package_and_submodules_are_the_same_objects():
    for name, stem in EXPECTED_HOMES.items():
        submodule = importlib.import_module(f"myflopy.project.group.{stem}")
        real = getattr(submodule, name)
        assert getattr(group_pkg, name) is real, name
        assert getattr(facade, name) is real, name
        assert real.__module__ == submodule.__name__, name


def test_facade_keeps_the_private_surface_alive():
    # Tests and downstream code reach for these through the old module path.
    from myflopy.project.model_group import (  # noqa: F401
        TResultsNamespace,
        _coerce_models,
        _GroupSpatialView,
    )


def test_model_group_deprecated_shortcuts_survive_the_move():
    # The Phase 3 registry names the facade path; the split must not change it.
    import myflopy
    from myflopy._deprecation import registered_deprecations

    registered = set(registered_deprecations())
    declared = set(myflopy.__compatibility__)
    for pkg in ("rch", "chd", "drn", "ghb", "wel", "uzf"):
        old = f"myflopy.project.model_group.ModelGroup.{pkg}"
        assert old in registered and old in declared, old
