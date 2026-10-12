"""Grouped model-budget comparison (`GroupBudget`)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from myflopy.project.group._shared import _reduce_to_period_end, _stable_compare_keys

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

        One row per record entry of every member's cell-budget file, as each
        member's ``model.bud(package).df`` holds it, with a ``model`` column naming
        the member.

        Parameters
        ----------
        package : str or None, optional
            The budget record to read, such as ``"rch"`` or ``"drn"``, lower-cased
            and handed to every member's ``model.bud(package)`` -- matched there,
            case-insensitively, as a SUBSTRING of the record names in that
            member's budget file (see
            :meth:`~myflopy.modflow.mf6.simulation.base.SimulationBase.bud`). Default
            ``None`` uses the package given to :meth:`ModelGroup.bud`; with
            neither, ``ValueError``. A member whose file holds no matching record
            raises ``ValueError`` from its own ``bud()``, which ends the whole call.

        Returns
        -------
        pandas.DataFrame
            ``node`` (the model-wide node, zero-based: ``layer * ncpl + cell`` on
            a DISV), ``kstpkper`` (zero-based ``(kstp, kper)`` on that member's
            own budget-file clock), the record's own columns (``node2``, ``q`` in
            model length^3/time with MF6's sign, and any auxiliary column) and
            ``model``.
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

        Each member is first reduced to its own LAST saved step of every period
        (``_reduce_to_period_end``, ledgers 337 and 348), so members of different
        ``nstp`` align by stress period, and a record is paired with the
        reference's on its integer and text identity -- ``node``, ``node2``,
        ``per`` and any non-float column -- by an inner join, so a record with no
        partner is dropped. Until ledger 485 the join matched ``kstpkper``
        EXACTLY and every float column, so with one step against two it paired
        the reference's period end with the other member's FIRST step, and an
        auxiliary that moved with the solution (a wetted area) dropped the row.

        Parameters
        ----------
        package : str or None, optional
            The budget record to compare, as :meth:`get` takes it (a
            case-insensitive substring of the member files' record names).
            Default ``None`` uses the package given to :meth:`ModelGroup.bud`;
            with neither, ``ValueError``.
        value_column : str, default "q"
            The column of :meth:`get` to difference -- ``"q"``, the flow in model
            length^3/time with MF6's sign, or a FLOAT auxiliary column. A name
            :meth:`get` does not have raises ``KeyError``. The other float
            columns are neither keys nor carried.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``reference_model``, ``kstpkper`` (the compared member's
            own step), the join keys, ``value_column``,
            ``reference_<value_column>`` and ``diff`` (member minus reference).
            Empty, with those columns, when no row pairs.
        """

        data = self.get(package=package)
        if data.empty:
            return pd.DataFrame()

        if value_column not in data.columns:
            raise KeyError(f"Budget value column {value_column!r} was not found.")

        data = _reduce_to_period_end(data)
        reference = self.group.reference
        key_columns = _stable_compare_keys(data, [value_column])
        ref = (
            data[data["model"] == reference][[*key_columns, value_column]]
            .rename(columns={value_column: f"reference_{value_column}"})
        )
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["diff"] = comp[value_column].astype(float) - comp[f"reference_{value_column}"].astype(float)
        ordered = [
            "model", "reference_model", "kstpkper", *key_columns,
            value_column, f"reference_{value_column}", "diff",
        ]
        return comp[ordered]
