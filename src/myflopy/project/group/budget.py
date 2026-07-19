"""Grouped model-budget comparison (`GroupBudget`)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupBudget:
    """Budget accessor for :class:`ModelGroup`."""

    def __init__(self, group: ModelGroup, package: str | None = None):
        """Bind the group budget accessor, optionally pinned to one ``package``."""

        self.group = group
        self.package = None if package is None else str(package).lower()

    # NOTE: there is deliberately no node re-normalization here. ``model.bud()``
    # already returns zero-based nodes (``budget_tables._zero_base_budget_frame``).
    # This class used to re-apply an ``if node.min() >= 1: subtract 1`` heuristic
    # on top of that, which is only valid on RAW MF6 records; on already-normalized
    # frames it silently shifted every package that does not happen to touch cell 0
    # one cell low. Pinned by tests/test_budget_node_basing.py.

    def get(
        self,
        *,
        package: str | None = None,
    ) -> pd.DataFrame:
        """Return aligned budget data for all models in the group.

        Parameters
        ----------
        package
            Optional package filter such as ``"rch"`` or ``"drn"``. If omitted,
            the package passed to :meth:`ModelGroup.bud` is used.
        """

        package_name = self.package if package is None else str(package).lower()
        if package_name is None:
            raise ValueError("A budget package name is required, for example group.bud('rch').get().")

        frames = []
        for model_name, model in self.group.models.items():
            budget = model.bud(package_name)
            frame = budget.df.reset_index().copy()
            renamed = {}
            if "level_0" in frame.columns:
                renamed["level_0"] = "node"
            if "level_1" in frame.columns:
                renamed["level_1"] = "kstpkper"
            if renamed:
                frame = frame.rename(columns=renamed)
            frame["model"] = model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)

    def compare(
        self,
        *,
        package: str | None = None,
        value_column: str = "q",
    ) -> pd.DataFrame:
        """Compare budget values for all models against the reference model.

        Parameters
        ----------
        package
            Optional package override.
        value_column
            Budget value column to compare. Defaults to ``"q"``.
        """

        data = self.get(package=package)
        if data.empty:
            return pd.DataFrame()

        if value_column not in data.columns:
            raise KeyError(f"Budget value column {value_column!r} was not found.")

        reference = self.group.reference
        key_columns = [column for column in data.columns if column not in {"model", value_column}]
        ref = (
            data[data["model"] == reference]
            .rename(columns={value_column: f"reference_{value_column}"})
            .drop(columns=["model"])
        )
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["diff"] = comp[value_column].astype(float) - comp[f"reference_{value_column}"].astype(float)
        ordered = ["model", "reference_model", *key_columns, value_column, f"reference_{value_column}", "diff"]
        return comp[ordered]


