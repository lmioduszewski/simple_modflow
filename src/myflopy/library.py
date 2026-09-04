"""Typed, browsable libraries for a project's reusable specs (design step 3).

A :class:`~myflopy.workspace.Project` holds three registries -- packages, grids
and simulations -- keyed by free-form library keys like ``"npf/base"`` or
``"rch_pond"``. They were plain ``dict``s, so an editor could offer ``.keys()``
and nothing else, and there was no way to see what a project contained without
knowing the keys already.

These are ``dict`` **subclasses**: everything a dict does still works
(subscript, ``in``, ``len``, iteration, ``.get``, ``.items``, ``.pop``,
assignment), and the entry type is fixed at the class statement so
``project.packages["npf/base"].with_options(...)`` still resolves in PyCharm.
On top of that each answers the same three questions -- what is defined
(:meth:`~SpecLibrary.names`), what is each one (:meth:`~SpecLibrary.summary`),
and give me that one (subscript, with a did-you-mean on a miss).

Library keys are **not** attribute-completable and never will be: they contain
``/``, so they are not identifiers. They are answered at runtime instead --
``names()``, ``summary()``, the repr, and the ``KeyError`` message.
"""

from __future__ import annotations

import difflib
from collections.abc import Iterable, Mapping
from html import escape
from typing import TYPE_CHECKING, Any, ClassVar, NoReturn, TypeVar

from myflopy.specs import (
    GridRef,
    GridSpec,
    PackageRef,
    PackageSpec,
    SimulationSpec,
    _builder_label,
    _package_entry_label,
    _summarize_value,
)

if TYPE_CHECKING:
    import pandas as pd

_T = TypeVar("_T")

__all__ = [
    "GridLibrary",
    "PackageLibrary",
    "SimulationLibrary",
    "SpecLibrary",
]


class SpecLibrary(dict[str, _T]):
    """A named library of specs: a ``dict``, plus the verbs every library answers.

    Subclasses declare ``_kind``/``_columns`` and override :meth:`_row`. Nothing
    else -- four hand-written near-identical classes is the maintenance cost
    ``CLAUDE.md`` §4.7.6 measured and rejected.

    ``_project_name`` is a **class-level default rather than a constructor
    argument** on purpose: ``dict.copy()`` and ``dict.__or__`` return a plain
    ``dict``, so a library carrying required instance state would silently
    degrade to a bare dict through either. :meth:`copy` and :meth:`__or__` are
    overridden to preserve the type, and the name is attached by the project
    after construction; absent, messages simply omit it.
    """

    _kind: ClassVar[str] = "entry"
    _columns: ClassVar[tuple[str, ...]] = ("key",)
    _project_name: str | None = None

    # -- the single override point -------------------------------------------

    def _row(self, key: str, entry: _T) -> dict[str, Any]:
        """Return one :meth:`summary` row for ``entry``. Subclasses override this."""

        return {"key": key}

    # -- discover -------------------------------------------------------------

    def names(self) -> list[str]:
        """Return the defined keys, sorted.

        The spelling for "what is in here?" -- ``sorted(library)`` also works,
        this just says what it is for.
        """

        return sorted(self)

    def summary(self) -> pd.DataFrame:
        """Return one row per entry: what each key holds.

        Columns are ``_columns``, in that order, and are present even when the
        library is empty -- so ``.empty`` is the test, never a missing column.
        Reads nothing from disk and builds nothing.
        """

        import pandas as pd

        rows = [self._row(key, self[key]) for key in self.names()]
        return pd.DataFrame(rows, columns=list(self._columns))

    # -- inspect --------------------------------------------------------------

    def __missing__(self, key: str) -> NoReturn:
        """Raise a ``KeyError`` naming the near miss and everything defined.

        ``dict.__getitem__`` calls this on a miss, so the message lands on the
        expression people already type. ``.get(key, default)`` is untouched --
        ``dict.get`` never calls ``__missing__``.
        """

        where = (
            f"Project '{self._project_name}'"
            if self._project_name
            else f"This {self._kind} library"
        )
        lines = [f"{where} has no {self._kind} '{key}'."]
        close = difflib.get_close_matches(key, self.names(), n=1, cutoff=0.6)
        if close:
            lines.append(f"Did you mean: '{close[0]}'?")
        if self:
            lines.append(
                f"{len(self)} {self._kind}s: " + ", ".join(self.names())
            )
        else:
            lines.append(f"No {self._kind}s are defined.")
        raise KeyError(" ".join(lines))

    def __repr__(self) -> str:
        """Compact one-liner: how many, and the first few keys."""

        if not self:
            return f"{type(self).__name__}(empty)"
        shown = self.names()[:4]
        more = "" if len(self) <= 4 else f", +{len(self) - 4}"
        listed = ", ".join(repr(name) for name in shown)
        return f"{type(self).__name__}({len(self)} {self._kind}s: {listed}{more})"

    def _repr_html_(self) -> str:
        """Return the :meth:`summary` table, for notebooks."""

        if not self:
            return (
                f"<div><em>No {escape(self._kind)}s defined.</em></div>"
            )
        return (
            f"<div><strong>{escape(type(self).__name__)}</strong> "
            f"({len(self)} {escape(self._kind)}s)"
            f"{self.summary().to_html(index=False)}</div>"
        )

    # -- dict behaviour that must not silently degrade ------------------------

    def copy(self) -> SpecLibrary[_T]:
        """Return a copy that is still a library. ``dict.copy`` returns a plain dict."""

        clone = type(self)(self)
        clone._project_name = self._project_name
        return clone

    def __or__(self, other: Mapping[str, _T]) -> SpecLibrary[_T]:
        """Merge, preserving the library type. ``dict.__or__`` returns a plain dict."""

        merged = self.copy()
        merged.update(other)
        return merged

    def __ror__(self, other: Mapping[str, _T]) -> SpecLibrary[_T]:
        """Reflected merge, preserving the library type."""

        merged = type(self)(other)
        merged._project_name = self._project_name
        merged.update(self)
        return merged


class PackageLibrary(SpecLibrary[PackageSpec]):
    """The project's reusable package specs, keyed by library key."""

    _kind: ClassVar[str] = "package"
    _columns: ClassVar[tuple[str, ...]] = (
        "key",
        "package",
        "builder",
        "enabled",
        "options",
        "inputs",
        "requires",
    )

    def _row(self, key: str, entry: PackageSpec) -> dict[str, Any]:
        """Summarize one package: what it builds and what it sets."""

        return {
            "key": key,
            "package": entry.name,
            "builder": _builder_label(entry.builder),
            "enabled": entry.enabled,
            "options": ", ".join(sorted(entry.options)) or "—",
            "inputs": ", ".join(sorted(entry.inputs)) or "—",
            "requires": ", ".join(entry.requires) or "—",
        }

    def persistence(self) -> pd.DataFrame:
        """Report how each package would be persisted by ``Project.save()``.

        A deliberate method rather than a :meth:`summary` column: deciding it
        means attempting ``to_dict()``, which is **119 ms** on a 100 000-record
        list BC (measured), and Jupyter calls ``_repr_html_`` on every bare
        ``project``. A package carrying arrays short-circuits immediately; one
        that is fully JSON-able does not, so the cost lands exactly where it is
        least affordable.

        ``pickle`` is a cache, not an archive -- it is the column worth watching.
        """

        import pandas as pd

        rows = []
        for key in self.names():
            try:
                self[key].to_dict()
            except ValueError as error:
                rows.append({"key": key, "persists": "pickle", "reason": str(error)})
            else:
                rows.append({"key": key, "persists": "json", "reason": "—"})
        return pd.DataFrame(rows, columns=["key", "persists", "reason"])


class GridLibrary(SpecLibrary["GridSpec | GridRef"]):
    """The project's reusable grids -- recipes, or already-built meshes."""

    _kind: ClassVar[str] = "grid"
    _columns: ClassVar[tuple[str, ...]] = (
        "key",
        "grid",
        "type",
        "method",
        "crs",
        "built",
        "cells",
    )

    def _row(self, key: str, entry: GridSpec | GridRef) -> dict[str, Any]:
        """Summarize one grid. Never resolves a recipe -- a summary that meshes
        a grid is a summary nobody dares print, so ``cells`` is blank until the
        grid is actually built."""

        if isinstance(entry, GridRef):
            return {
                "key": key,
                "grid": entry.key,
                "type": "ref",
                "method": "—",
                "crs": "—",
                "built": False,
                "cells": "—",
            }
        built = entry.method == "object" and entry.obj is not None
        cells = "—"
        if built:
            ncpl = getattr(entry.obj, "ncpl", None)
            cells = "—" if ncpl is None else ncpl
        return {
            "key": key,
            "grid": entry.name,
            "type": entry.grid_type,
            "method": entry.method,
            "crs": str(entry.crs) if entry.crs else "—",
            "built": built,
            "cells": cells,
        }


class SimulationLibrary(SpecLibrary[SimulationSpec]):
    """The project's simulation specs, keyed by simulation name."""

    _kind: ClassVar[str] = "simulation"
    _columns: ClassVar[tuple[str, ...]] = (
        "name",
        "models",
        "packages",
        "exchanges",
        "derived_from",
    )

    def _row(self, key: str, entry: SimulationSpec) -> dict[str, Any]:
        """Summarize one simulation: its models, what they carry, its lineage."""

        packages = [
            _package_entry_label(package)
            for model in entry.models
            for package in model.packages
        ]
        return {
            "name": key,
            "models": ", ".join(entry.model_names) or "—",
            "packages": len(packages),
            "exchanges": len(entry.exchanges),
            "derived_from": entry.derived_from or "—",
        }
