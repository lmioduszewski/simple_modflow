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
from collections.abc import Mapping
from html import escape
from typing import TYPE_CHECKING, Any, ClassVar, NoReturn, TypeVar

from myflopy.specs import (
    GridRef,
    GridSpec,
    PackageSpec,
    SimulationSpec,
    _builder_label,
    _package_entry_label,
    fork_package,
    package_digest,
    ref,
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

    def _where(self) -> str:
        """How to name this library in a message -- with its project when known."""

        return (
            f"Project '{self._project_name}'"
            if self._project_name
            else f"This {self._kind} library"
        )

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

        lines = [f"{self._where()} has no {self._kind} '{key}'."]
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
        # Copy ALL instance state, not just the name: a subclass that grows
        # another binding (RunLibrary's package library) would otherwise come
        # back detached and answer as if it had nothing to compare against.
        clone.__dict__.update(self.__dict__)
        return clone

    def __or__(self, other: Mapping[str, _T]) -> SpecLibrary[_T]:
        """Merge, preserving the library type. ``dict.__or__`` returns a plain dict."""

        merged = self.copy()
        merged.update(other)
        return merged

    def __ror__(self, other: Mapping[str, _T]) -> SpecLibrary[_T]:
        """Reflected merge, preserving the library type."""

        merged = type(self)(other)
        merged.__dict__.update(self.__dict__)
        merged.update(self)
        return merged


class PackageLibrary(SpecLibrary[PackageSpec]):
    """The project's reusable package specs, keyed by library key."""

    _kind: ClassVar[str] = "package"
    _digest_cache: dict[str, str | None] | None = None
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

    def digests(self) -> dict[str, str | None]:
        """Fingerprint every package, caching until the library changes.

        Digesting is cheap for array-valued options (a 72x9405 recharge array is
        ~5 ms) and expensive for record lists, because the cost is pickling
        Python objects rather than bytes: the USG-derived 72-period, 9,090-record
        EVT measures **1.3 s**. ``RunLibrary.summary`` fingerprints the whole
        library on every call and Jupyter calls that on every bare
        ``project.runs``, so without this the browse surface stalls on exactly
        the projects worth browsing.

        The cache is dropped whenever the library is mutated. It is keyed on the
        library, not on spec identity, so replacing an entry through any route --
        ``edit``, ``fork``, plain assignment, ``pop`` -- invalidates it. Mutating
        a spec's ``options`` dict IN PLACE bypasses it, as it also bypasses
        ``_check_options``; that is not a supported way to change a package.
        """

        cached = getattr(self, "_digest_cache", None)
        if cached is None:
            cached = {key: package_digest(entry) for key, entry in self.items()}
            self._digest_cache = cached
        return cached

    def _invalidate(self) -> None:
        """Drop the digest cache. Called from every mutating dict operation."""

        self._digest_cache = None

    def __setitem__(self, key: str, value: PackageSpec) -> None:
        """Register a package, dropping any cached fingerprints."""

        self._invalidate()
        super().__setitem__(key, value)

    def __delitem__(self, key: str) -> None:
        """Remove a package, dropping any cached fingerprints."""

        self._invalidate()
        super().__delitem__(key)

    def pop(self, *args: Any) -> Any:
        """Remove and return a package, dropping any cached fingerprints."""

        self._invalidate()
        return super().pop(*args)

    def popitem(self) -> tuple[str, PackageSpec]:
        """Remove and return one package, dropping any cached fingerprints."""

        self._invalidate()
        return super().popitem()

    def update(self, *args: Any, **kwargs: Any) -> None:
        """Merge packages in, dropping any cached fingerprints."""

        self._invalidate()
        super().update(*args, **kwargs)

    def setdefault(self, key: str, default: Any = None) -> Any:
        """Fetch or insert a package, dropping any cached fingerprints."""

        self._invalidate()
        return super().setdefault(key, default)

    def clear(self) -> None:
        """Empty the library, dropping any cached fingerprints."""

        self._invalidate()
        super().clear()

    def __ior__(self, other: Mapping[str, PackageSpec]) -> PackageLibrary:
        """In-place merge, dropping any cached fingerprints.

        ``|=`` does NOT route through ``__or__`` or ``update`` on a dict
        subclass, so it slipped past every other invalidation hook and left a
        stale fingerprint behind -- which reports a changed package as
        unchanged, the one thing drift must never do.
        """

        self._invalidate()
        super().__ior__(other)
        return self

    def fork(
        self,
        source: str,
        key: str,
        *,
        name: str | None = None,
        deep: bool = True,
        note: str | None = None,
        overwrite: bool = False,
        **overrides: Any,
    ) -> PackageSpec:
        """Register a copy of ``source`` under ``key``, with ``overrides`` applied.

        The copy is BUILT through the ordinary construction path, so a typo in
        ``overrides`` raises here exactly as it would in ``mf.npf(...)`` -- see
        :func:`~myflopy.specs.fork_package`. ``deep=True`` also unaliases the
        options the fork did not touch, so editing the fork cannot write through
        into the original's arrays.

        Mutates the library and returns the new spec. Refuses an existing ``key``
        unless ``overwrite=True``: silently replacing a package that other
        simulations reference is the accident this verb exists to avoid.
        """

        if key in self and not overwrite:
            raise KeyError(
                f"{self._where()} already has a package '{key}'. "
                "Pass overwrite=True to replace it, or pick another key."
            )
        forked = fork_package(
            self[source], name=name, deep=deep, note=note, source=source, **overrides
        )
        self[key] = forked
        return forked

    def edit(
        self,
        key: str,
        *,
        name: str | None = None,
        deep: bool = True,
        note: str | None = None,
        **overrides: Any,
    ) -> PackageSpec:
        """Replace the package at ``key`` with an edited copy of itself.

        The same operation as :meth:`fork` with a different destination, spelled
        apart on purpose: :meth:`fork` leaves the original alone, this changes
        what an existing key MEANS. Every run that already used ``key`` recorded
        the key, not its contents, so this is the verb that can invalidate
        history -- ``project.runs.drift()`` is what makes that visible.

        Editing is never refused. It is your library, and revising ``npf/base``
        is the normal case; the point is that the consequence stops being silent.
        """

        edited = fork_package(
            self[key],
            name=name,
            deep=deep,
            note=note,
            source=key,
            stamp="edited",
            **overrides,
        )
        self[key] = edited
        return edited

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

    def fork(
        self,
        source: str,
        name: str,
        *,
        swap: Mapping[str, Any] | None = None,
        register: bool = True,
        overwrite: bool = False,
    ) -> SimulationSpec:
        """Register a variant of ``source`` under ``name``, swapping packages.

        Wraps :meth:`~myflopy.specs.SimulationSpec.derive`, which already records
        ``derived_from`` and appends to ``lineage``; each swap adds its own
        ``replace_package`` entry, so the fork carries an auditable trail.

        ``swap`` maps an existing package label to the library key that should
        take its slot, in either of two forms::

            swap={"rch/series": "rch_pond"}            # every model that has it
            swap={"flow": {"rch/series": "rch_pond"}}  # one named model

        The flat form applies to **every model carrying that label**, which is
        what you want for a single-model simulation and for genuinely shared
        packages. A coupled GWF/GWT simulation has ``ic``, ``oc`` and ``dis`` on
        both models, so name the model when you mean one of them. A label that
        matches nothing raises -- a typo that silently changed no packages would
        hand back a "variant" identical to its parent.

        ``register=False`` returns the spec without adding it. Note that a spec
        whose ``derived_from`` names a simulation NOT in the library makes
        :meth:`Project.save` refuse it, so an unregistered parent is a save-time
        error rather than a fork-time one.
        """

        if name in self and not overwrite:
            raise KeyError(
                f"{self._where()} already has a simulation '{name}'. "
                "Pass overwrite=True to replace it, or pick another name."
            )
        # Resolved against the SOURCE, so a bad label names the simulation the
        # caller asked about rather than the derived name they have not seen yet.
        targets = _swap_targets(self[source], swap or {})
        spec = self[source].derive(name)
        for model_name, old_label, new_key in targets:
            # Target the FULL LABEL, not the entry's slot hint: a slot is only
            # the key's first segment for a ref, so wel/north and wel/south both
            # answer to "wel" and the swap landed on whichever came first.
            # slot= keeps replace_package's must-exist assertion: a swap that
            # misses would ADD, and MF6 sums two RCH packages under a clean run.
            spec = spec.with_packages(
                ref(new_key), on=model_name, slot=old_label
            )
        if register:
            self[name] = spec
        return spec

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


def _swap_targets(
    spec: SimulationSpec, swap: Mapping[str, Any]
) -> list[tuple[str, str, str]]:
    """Resolve a ``swap`` mapping to ``(model, old label, new key)`` triples.

    Accepts the flat form (``{"rch/series": "rch_pond"}``, every model carrying
    the label) and the model-qualified form
    (``{"flow": {"rch/series": "rch_pond"}}``, exactly that model); a mapping
    VALUE marks the qualified form, so the two cannot be confused.

    Every label is validated in BOTH forms. A label matching nothing raises: a
    typo that silently swapped nothing would hand back a "variant" identical to
    its parent, which is the failure this check exists for -- and the qualified
    form is the one the docs tell you to use for coupled models, so it needs the
    guarantee at least as much.

    Two entries targeting the same model and label also raise rather than letting
    dict order decide which wins.
    """

    model_names = list(spec.model_names)
    simulation_labels = {
        _package_entry_label(entry) for entry in spec.packages
    }
    targets: list[tuple[str, str, str]] = []
    claimed: dict[tuple[str, str], str] = {}

    def labels_of(model_name: str) -> list[str]:
        return [
            _package_entry_label(entry)
            for entry in spec.model(model_name).packages
        ]

    def inventory() -> str:
        carried = sorted(
            {label for model in model_names for label in labels_of(model)}
            | simulation_labels
        )
        return ", ".join(carried)

    def claim(model_name: str, old_label: str, new_key: str) -> None:
        if not isinstance(new_key, str):
            raise TypeError(
                f"swap[{old_label!r}] must be a library key (a string), not "
                f"{type(new_key).__name__}. Pass the key the package is "
                "registered under."
            )
        previous = claimed.get((model_name, old_label))
        if previous is not None and previous != new_key:
            raise ValueError(
                f"swap targets '{old_label}' on model '{model_name}' twice, "
                f"with '{previous}' and '{new_key}'. Name it once -- otherwise "
                "which one applies depends on the order the mapping was written."
            )
        claimed[(model_name, old_label)] = new_key
        targets.append((model_name, old_label, new_key))

    def refuse(old_label: str, where: str) -> None:
        if old_label in simulation_labels:
            raise KeyError(
                f"'{old_label}' is a simulation-level package of "
                f"'{spec.name}' (tdis/ims and friends), which swap= cannot "
                "replace yet. Build the variant with mf.simulation(...) instead."
            )
        raise KeyError(
            f"{where} has no package '{old_label}' to swap. "
            f"It carries: {inventory()}."
        )

    # Qualified entries first, so a mixed mapping resolves the same either way.
    for outer, value in swap.items():
        if not isinstance(value, Mapping):
            continue
        if outer not in model_names:
            raise KeyError(
                f"Simulation '{spec.name}' has no model '{outer}'. "
                f"Models: {', '.join(sorted(model_names))}."
            )
        for old_label, new_key in value.items():
            if old_label not in labels_of(outer):
                refuse(old_label, f"Model '{outer}' of simulation '{spec.name}'")
            claim(outer, old_label, new_key)

    for outer, value in swap.items():
        if isinstance(value, Mapping):
            continue
        matched = [model for model in model_names if outer in labels_of(model)]
        if not matched:
            refuse(outer, f"Simulation '{spec.name}'")
        for model_name in matched:
            claim(model_name, outer, value)

    return targets
