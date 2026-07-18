"""Single<->group inspection symmetry (Phase 4).

The group side mirrors the single-model side: ``group.packages.<pkg>.inputs`` /
``.results`` matches ``model.packages.<pkg>.inputs`` / ``.results``. The old flat
``group.<pkg>`` shortcuts (which have no single-model equivalent) are deprecated
but still work, returning the very same accessor as the mirrored path.
"""

from __future__ import annotations

import warnings

import pandas as pd
import pytest

from myflopy.project.model_config import ModelConfig
from myflopy.project.model_group import GroupPackageInputs, ModelGroup


class _FakeModel:
    def __init__(self, name, package_names=("GHB",)):
        self.name = name
        self.package_names = list(package_names)
        self.config = ModelConfig(pd.DataFrame(columns=["section", "setting", "value"]))


@pytest.fixture
def group():
    return ModelGroup(
        {"ref": _FakeModel("ref"), "alt": _FakeModel("alt")}, reference="ref"
    )


def test_group_packages_expose_inputs_and_results(group):
    # riv/evt are included: they are cell-stress list BCs like the rest, so the
    # group accessors must exist. They deliberately have NO deprecated flat
    # shortcut (see the parametrized test below) -- new packages ship only the
    # mirrored group.packages.<pkg>.inputs path.
    for pkg in ("rch", "chd", "drn", "ghb", "riv", "wel", "evt"):
        accessor = getattr(group.packages, pkg)
        assert hasattr(accessor, "inputs"), pkg
        assert hasattr(accessor, "results"), pkg
        # each accessor is wired to its OWN package (a copy-paste of another
        # package's private attr would otherwise pass the hasattr checks)
        assert accessor.inputs.package_name == pkg
    assert hasattr(group.packages.uzf, "inputs")
    assert hasattr(group.packages.uzf, "results")


def test_diff_tier_covers_exactly_the_group_cell_bc_accessors(group):
    """The hardcoded diff package lists must track the group's list-BC accessors.

    Regression guard (2026-07-18): adding ``mf.riv``/``mf.evt`` added group
    accessors but not the diff lists, so ``group.diff().packages.riv`` raised
    AttributeError while ``group.packages.riv`` worked -- and the comment above
    ``_DIFF_PACKAGES`` claimed the two matched. Neither list is registry-driven,
    so pin the invariant instead of trusting the comment.
    """

    from myflopy.project.model_diff import _DIFF_PACKAGES
    from myflopy.project.model_results_diff import _CELL_BUDGET_PACKAGES

    group_bcs = {
        value.package_name
        for value in vars(group).values()
        if isinstance(value, GroupPackageInputs)
    }
    assert group_bcs, "expected the group to expose cell-stress BC accessors"
    assert group_bcs == set(_DIFF_PACKAGES), (
        "every group list-BC accessor must be diffable: "
        f"group-only={sorted(group_bcs - set(_DIFF_PACKAGES))}, "
        f"diff-only={sorted(set(_DIFF_PACKAGES) - group_bcs)}"
    )
    # the results side carries the same BCs (plus the advanced sfr/lak terms)
    assert group_bcs <= set(_CELL_BUDGET_PACKAGES), (
        "cell-budget diff is missing: "
        f"{sorted(group_bcs - set(_CELL_BUDGET_PACKAGES))}"
    )


def test_inputs_accessor_exposes_the_standard_verbs(group):
    inputs = group.packages.ghb.inputs
    for verb in ("get", "summary", "map"):
        assert callable(getattr(inputs, verb)), verb


def test_mirror_path_emits_no_deprecation_warning(group):
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        accessor = group.packages.rch.inputs
    assert isinstance(accessor, GroupPackageInputs)


@pytest.mark.parametrize("pkg", ["rch", "chd", "drn", "ghb", "wel", "uzf"])
def test_flat_shortcut_deprecated_but_returns_the_mirror_object(group, pkg):
    with pytest.warns(DeprecationWarning):
        flat = getattr(group, pkg)
    # The deprecated flat path and the mirrored path are the same accessor.
    assert flat is getattr(group.packages, pkg).inputs


def test_group_inputs_summary_counts_cells_and_periods(monkeypatch):
    # Patch where the builder is USED: group/inputs.py (plan 4.2 split).
    from myflopy.project.group import inputs as group_inputs

    def stub(model, package_name, *, per=None, layer=None, cells=None):
        data = {
            "ref": [(0, 0, 10), (0, 0, 11), (1, 0, 10)],
            "alt": [(0, 0, 10)],
        }[model.name]
        return pd.DataFrame(data, columns=["per", "layer", "cell"]).assign(
            package="ghb", model=model.name
        )

    monkeypatch.setattr(group_inputs, "build_cell_package_input_table", stub)
    grp = ModelGroup({"ref": _FakeModel("ref"), "alt": _FakeModel("alt")}, reference="ref")
    summary = grp.packages.ghb.inputs.summary().set_index("model")
    assert summary.loc["ref", "cells"] == 2  # cells 10, 11
    assert summary.loc["ref", "periods"] == 2  # periods 0, 1
    assert summary.loc["alt", "cells"] == 1


def test_diff_engine_does_not_trip_the_deprecated_flat_path(group):
    # ModelDiff builds its input accessors directly; exercising it must not emit
    # the flat-shortcut DeprecationWarning.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        diff = group.diff()
        diff.summary()
        diff.report()
