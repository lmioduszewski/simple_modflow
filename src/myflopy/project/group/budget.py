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

    @staticmethod
    def _normalize_nodes(frame: pd.DataFrame) -> pd.DataFrame:
        """Normalize budget node columns to zero-based indexing when needed."""

        normalized = frame.copy()
        for column in ("node", "node2"):
            if column in normalized.columns:
                numeric = pd.to_numeric(normalized[column], errors="coerce")
                mask = numeric.notna()
                if mask.any():
                    ints = numeric.loc[mask].astype(int)
                    if int(ints.min()) >= 1:
                        normalized[column] = numeric
                        normalized.loc[mask, column] = (ints - 1).astype(float)
                    else:
                        normalized[column] = numeric
                        normalized.loc[mask, column] = ints.astype(float)
        return normalized

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
            frame = self._normalize_nodes(frame)
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


