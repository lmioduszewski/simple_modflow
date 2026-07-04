"""Unified multi-model difference engine -- the ``diff`` verb.

:class:`ModelDiff` is the single entry point for "what is different between these
models?" It is built from a :class:`~myflopy.project.model_group.ModelGroup`
(reference-star topology: the group's reference vs each other model) and reports
tiers of difference. Phase 1 ships two:

* **structural** -- which packages, and which cells within a package, exist in one
  model but not another. This is a true set difference: unlike a value-only
  compare (which inner-joins on shared cells) it surfaces cells that appear in
  only one model instead of silently dropping them.
* **value** -- numeric field differences on the cells shared by both models
  (this reuses the group's existing aligned-value ``compare``).

Two doors reach the same engine, and both return a :class:`ModelDiff`:

* ``model.diff(other)`` / ``model.diff([m2, m3])`` -- quick and ad hoc; ``self``
  is the reference. Accepts model objects or workspace paths.
* ``ModelGroup(...).diff()`` -- durable and named, with explicit reference control.

All stress-period and layer arguments are zero-based, matching the rest of the
myflopy API.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.package_tables import build_cell_package_input_table
from myflopy.project.model_group import GroupPackageInputs

if TYPE_CHECKING:  # pragma: no cover - typing only
    from myflopy.project.model_group import ModelGroup

# Cell-based stress-period BC packages the Phase-1 diff understands. These match
# the packages the ModelGroup exposes as GroupPackageInputs accessors.
_DIFF_PACKAGES: tuple[str, ...] = ("rch", "chd", "drn", "ghb", "wel")

_SUMMARY_COLUMNS = [
    "package",
    "model",
    "present_in_reference",
    "present_in_model",
    "cells_only_in_reference",
    "cells_only_in_model",
    "cells_shared",
    "value_cells_changed",
    "identical",
]


def _package_present(model, package_name: str) -> bool:
    """Return whether ``model`` has a package of the given type attached."""

    try:
        names = [str(name).lower() for name in model.package_names]
    except Exception:  # pragma: no cover - defensive
        return False
    pkg = package_name.lower()
    return any(name == pkg or name.startswith(pkg) for name in names)


def _package_keys(model, package_name: str, *, per=None, layer=None):
    """Return ``(present, keys)`` for a BC package on ``model``.

    ``keys`` is a set of zero-based ``(per, layer, cell)`` tuples identifying
    every cell the package applies to (optionally filtered by ``per``/``layer``).
    A package that is absent, empty, or unreadable yields an empty set; the
    ``present`` flag distinguishes "attached but no cells" from "not attached".
    """

    try:
        table = build_cell_package_input_table(model, package_name, per=per, layer=layer)
    except Exception:
        return _package_present(model, package_name), set()
    if table is None or table.empty:
        return _package_present(model, package_name), set()
    n = len(table)
    per_col = table["per"].to_numpy(dtype=int) if "per" in table.columns else np.zeros(n, dtype=int)
    layer_col = table["layer"].to_numpy(dtype=int) if "layer" in table.columns else np.zeros(n, dtype=int)
    cell_col = table["cell"].to_numpy(dtype=int)
    keys = set(zip(per_col.tolist(), layer_col.tolist(), cell_col.tolist()))
    return True, keys


class PackageDiff:
    """Difference view for one cell-based BC package across the group.

    Reached via ``model_diff.packages.<package>`` (e.g. ``.ghb``). Exposes the
    structural tier (:meth:`cells`), the value tier (:meth:`values`, :meth:`map`),
    and per-model counts (:meth:`summary`), always relative to the group's
    reference model.
    """

    def __init__(self, diff: "ModelDiff", package_name: str):
        self._diff = diff
        self.group = diff.group
        self.package_name = str(package_name).lower()
        self._inputs = GroupPackageInputs(self.group, self.package_name)

    def _targets(self, model_name=None) -> list[str]:
        """Return the non-reference model names to diff (or just ``model_name``)."""

        if model_name is not None:
            name = str(model_name)
            if name not in self.group.models:
                raise KeyError(f"Model {name!r} is not in the group.")
            if name == self.group.reference:
                raise ValueError("The reference model cannot be diffed against itself.")
            return [name]
        return [name for name in self.group.models if name != self.group.reference]

    # -- value tier (reuses the group's aligned-value compare) ----------------
    def values(self, **kwargs) -> pd.DataFrame:
        """Aligned value table: each numeric field as ``x`` / ``reference_x`` /
        ``x_diff`` on the cells shared with the reference model.

        Accepts the same ``model_name`` / ``per`` / ``layer`` / ``cells``
        filters as the group input accessor.
        """

        return self._inputs.compare(**kwargs)

    def map(self, **kwargs):
        """Diverging choropleth of a field's difference vs the reference model."""

        return self._inputs.compare_map(**kwargs)

    # -- structural tier (the new set-difference engine) ----------------------
    def cells(
        self,
        *,
        model_name=None,
        per=None,
        layer=None,
        include_shared: bool = False,
    ) -> pd.DataFrame:
        """Structural set difference of applied cells vs the reference model.

        Returns a tidy frame with columns ``model, per, layer, cell,
        membership`` where ``membership`` is ``only_in_reference`` (the cell is
        in the reference but not this model), ``only_in_model`` (in this model
        but not the reference), or -- when ``include_shared=True`` --
        ``shared``. This is the difference a value-only compare cannot show,
        because its inner join drops cells present in only one model.
        """

        reference = self.group.reference
        _, ref_keys = _package_keys(
            self.group.models[reference], self.package_name, per=per, layer=layer
        )
        rows: list[tuple] = []
        for name in self._targets(model_name):
            _, model_keys = _package_keys(
                self.group.models[name], self.package_name, per=per, layer=layer
            )
            for key in sorted(ref_keys - model_keys):
                rows.append((name, *key, "only_in_reference"))
            for key in sorted(model_keys - ref_keys):
                rows.append((name, *key, "only_in_model"))
            if include_shared:
                for key in sorted(ref_keys & model_keys):
                    rows.append((name, *key, "shared"))
        return pd.DataFrame(rows, columns=["model", "per", "layer", "cell", "membership"])

    def summary(self, *, model_name=None) -> pd.DataFrame:
        """Per-model structural + value difference counts for this package."""

        reference = self.group.reference
        ref_present, ref_keys = _package_keys(self.group.models[reference], self.package_name)
        rows = []
        for name in self._targets(model_name):
            model_present, model_keys = _package_keys(self.group.models[name], self.package_name)
            only_ref = len(ref_keys - model_keys)
            only_model = len(model_keys - ref_keys)
            shared = len(ref_keys & model_keys)
            changed = self._value_cells_changed(name) if shared else 0
            identical = (
                ref_present == model_present
                and only_ref == 0
                and only_model == 0
                and changed == 0
            )
            rows.append(
                {
                    "package": self.package_name,
                    "model": name,
                    "present_in_reference": ref_present,
                    "present_in_model": model_present,
                    "cells_only_in_reference": only_ref,
                    "cells_only_in_model": only_model,
                    "cells_shared": shared,
                    "value_cells_changed": changed,
                    "identical": identical,
                }
            )
        return pd.DataFrame(rows, columns=_SUMMARY_COLUMNS)

    def _value_cells_changed(self, model_name: str) -> int:
        """Count shared cells whose numeric fields differ from the reference."""

        try:
            comparison = self._inputs.compare(model_name=model_name)
        except Exception:
            return 0
        if comparison.empty:
            return 0
        diff_columns = [column for column in comparison.columns if column.endswith("_diff")]
        if not diff_columns:
            return 0
        nonzero = (comparison[diff_columns].fillna(0.0) != 0.0).any(axis=1)
        return int(nonzero.sum())


class _PackageDiffNamespace:
    """Attribute access ``diff.packages.<package>`` -> :class:`PackageDiff`."""

    def __init__(self, diff: "ModelDiff"):
        self._diff = diff

    def __getattr__(self, name: str) -> PackageDiff:
        return self._diff.package(name)

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(self._diff.package_names))


class _FocusedModelDiff:
    """A :class:`ModelDiff` narrowed to a single non-reference model."""

    def __init__(self, diff: "ModelDiff", model_name: str):
        self._diff = diff
        self.model_name = str(model_name)
        if self.model_name not in diff.group.models:
            raise KeyError(f"Model {self.model_name!r} is not in the group.")
        if self.model_name == diff.group.reference:
            raise ValueError("The reference model cannot be diffed against itself.")

    def summary(self) -> pd.DataFrame:
        return self._diff.summary(model_name=self.model_name)

    def report(self) -> str:
        return self._diff._render_report(model_name=self.model_name)

    def cells(self, package: str, **kwargs) -> pd.DataFrame:
        return self._diff.package(package).cells(model_name=self.model_name, **kwargs)

    def values(self, package: str, **kwargs) -> pd.DataFrame:
        return self._diff.package(package).values(model_name=self.model_name, **kwargs)


class ModelDiff:
    """Reference-star difference across a :class:`ModelGroup` (the ``diff`` verb)."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

    @property
    def reference(self) -> str:
        """The baseline model every other model is compared against."""

        return self.group.reference

    @property
    def model_names(self) -> list[str]:
        """The non-reference model names being compared to the reference."""

        return [name for name in self.group.models if name != self.group.reference]

    @property
    def package_names(self) -> list[str]:
        """BC packages present in at least one model in the group."""

        return [
            pkg
            for pkg in _DIFF_PACKAGES
            if any(_package_present(model, pkg) for model in self.group.models.values())
        ]

    @property
    def packages(self) -> _PackageDiffNamespace:
        """Namespace for per-package diffs: ``diff.packages.ghb.cells()`` etc."""

        return _PackageDiffNamespace(self)

    def package(self, name: str) -> PackageDiff:
        """Return the :class:`PackageDiff` for one BC package."""

        pkg = str(name).lower()
        if pkg not in _DIFF_PACKAGES:
            raise AttributeError(
                f"ModelDiff does not diff package {name!r}; supported: {', '.join(_DIFF_PACKAGES)}."
            )
        return PackageDiff(self, pkg)

    def model(self, name: str) -> _FocusedModelDiff:
        """Focus the diff on a single non-reference model."""

        return _FocusedModelDiff(self, name)

    def summary(self, *, model_name=None) -> pd.DataFrame:
        """Package x model matrix of structural + value difference counts.

        One row per (package, non-reference model). ``identical`` is ``True``
        when a package matches the reference in presence, cell membership, and
        all field values.
        """

        frames = [
            self.package(pkg).summary(model_name=model_name) for pkg in self.package_names
        ]
        frames = [frame for frame in frames if not frame.empty]
        if not frames:
            return pd.DataFrame(columns=_SUMMARY_COLUMNS)
        return pd.concat(frames, ignore_index=True)

    def report(self) -> str:
        """Return a readable Markdown 'faithful-copy' report."""

        return self._render_report()

    def _render_report(self, *, model_name=None) -> str:
        reference = self.group.reference
        summary = self.summary(model_name=model_name)
        targets = [model_name] if model_name is not None else self.model_names
        lines = [f"# Model diff -- reference: `{reference}`", ""]
        if summary.empty:
            lines.append("_No comparable BC packages found in the group._")
            return "\n".join(lines)
        for name in targets:
            block = summary[summary["model"] == name]
            if not block.empty and bool(block["identical"].all()):
                lines.append(f"## `{name}` -- identical to reference")
                lines.append("")
                continue
            lines.append(f"## `{name}` -- differs from reference")
            lines.append("")
            lines.append(
                "| package | present (ref/model) | only-ref cells | only-model cells "
                "| shared | value cells changed |"
            )
            lines.append("|---|---|---|---|---|---|")
            for _, row in block.iterrows():
                lines.append(
                    f"| {row['package']} "
                    f"| {row['present_in_reference']}/{row['present_in_model']} "
                    f"| {row['cells_only_in_reference']} "
                    f"| {row['cells_only_in_model']} "
                    f"| {row['cells_shared']} "
                    f"| {row['value_cells_changed']} |"
                )
            lines.append("")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"ModelDiff(reference={self.group.reference!r}, "
            f"models={self.model_names!r})"
        )
