"""Grouped package/outputs namespaces (`GroupPackages`, `GroupOutputs`)."""

from __future__ import annotations

from typing import TYPE_CHECKING, Generic

from myflopy.modflow.mf6.package_explorer import get_default_budget_term
from myflopy.project.group._shared import TResultsNamespace
from myflopy.project.group.lak import (
    GroupLakBudgetResults,
    GroupLakOutputs,
    GroupLakPackageAccessor,
    GroupLakResultsNamespace,
)
from myflopy.project.group.results import GroupCellPackageResults, GroupCellPackageResultsNamespace
from myflopy.project.group.sfr import GroupSfrBudgetResults, GroupSfrResultsNamespace
from myflopy.project.group.surface_water import GroupSurfaceWaterResultsNamespace
from myflopy.project.group.uzf import GroupUzfPackageAccessor

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupOutputs:
    """Namespace for grouped package-specific output accessors."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped package-output namespace to ``group``."""

        self.group = group

    @property
    def lak(self) -> GroupLakOutputs:
        """Return grouped LAK output helpers."""

        return GroupLakOutputs(self.group)


class GroupPackageAccessor:
    """Namespace for one grouped package's preferred exploration helpers."""

    def __init__(self, accessor):
        """Wrap a grouped inputs ``accessor`` and expose ``.inputs`` / ``.results``."""

        self.inputs = accessor

    @property
    def results(self):
        """Return grouped result helpers for this package."""

        budget_text, value_name = get_default_budget_term(self.inputs.package_name) or (
            self.inputs.package_name.upper(),
            "q",
        )
        return GroupCellPackageResultsNamespace(
            GroupCellPackageResults(
                self.inputs.group,
                self.inputs.package_name,
                budget_text=budget_text,
                value_name=value_name,
            )
        )


class GroupResultsOnlyPackageAccessor(Generic[TResultsNamespace]):
    """Namespace for grouped packages that currently expose results only."""

    def __init__(self, results_namespace: TResultsNamespace):
        """Wrap a grouped results namespace for a package exposing results only."""

        self._results_namespace = results_namespace

    @property
    def results(self) -> TResultsNamespace:
        """Return grouped result helpers for this package."""

        return self._results_namespace


class GroupPackages:
    """Preferred grouped package exploration namespace.

    This mirrors the single-model ``model.packages`` surface where practical
    while reusing the existing grouped ``get()/compare()`` accessors.
    """

    def __init__(self, group: ModelGroup):
        """Bind the preferred grouped ``packages`` namespace to ``group``."""

        self.group = group

    @property
    def rch(self) -> GroupPackageAccessor:
        """Grouped recharge input helpers."""

        return GroupPackageAccessor(self.group._rch)

    @property
    def chd(self) -> GroupPackageAccessor:
        """Grouped constant-head input helpers."""

        return GroupPackageAccessor(self.group._chd)

    @property
    def drn(self) -> GroupPackageAccessor:
        """Grouped drain input helpers."""

        return GroupPackageAccessor(self.group._drn)

    @property
    def ghb(self) -> GroupPackageAccessor:
        """Grouped general-head-boundary input helpers."""

        return GroupPackageAccessor(self.group._ghb)

    @property
    def wel(self) -> GroupPackageAccessor:
        """Grouped well package input helpers."""

        return GroupPackageAccessor(self.group._wel)

    @property
    def uzf(self) -> GroupUzfPackageAccessor:
        """Grouped UZF input helpers."""

        return GroupUzfPackageAccessor(self.group._uzf)

    @property
    def sfr(self) -> GroupResultsOnlyPackageAccessor[GroupSfrResultsNamespace]:
        """Grouped SFR result helpers."""

        budget_text, value_name = get_default_budget_term("sfr") or ("SFR", "q")
        return GroupResultsOnlyPackageAccessor(
            GroupSfrResultsNamespace(
                GroupSfrBudgetResults(
                    self.group,
                    budget_text=budget_text,
                    value_name=value_name,
                )
            )
        )

    @property
    def lak(self) -> GroupLakPackageAccessor:
        """Grouped LAK geometry and result helpers."""

        budget_text, value_name = get_default_budget_term("lak") or ("GWF", "q")
        results_namespace = GroupLakResultsNamespace(
            GroupLakBudgetResults(
                self.group,
                budget_text=budget_text,
                value_name=value_name,
            )
        )
        return GroupLakPackageAccessor(self.group, results_namespace)

    @property
    def surface_water(self) -> GroupResultsOnlyPackageAccessor[GroupSurfaceWaterResultsNamespace]:
        """Grouped combined surface-water result helpers."""

        return GroupResultsOnlyPackageAccessor(GroupSurfaceWaterResultsNamespace(self.group))


