"""Composable specifications for building MODFLOW 6 simulations.

The classes in this module are intentionally small. A spec stores a named
builder and its keyword arguments, while the builder remains ordinary Python.
This keeps model assembly explicit and makes package alternatives easy to
replace without introducing a framework around FloPy.
"""

from __future__ import annotations

import pickle
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field, replace
from enum import Enum
from html import escape
from importlib import import_module
from pathlib import Path
from typing import Any, TypeAlias, TypeVar, cast

import flopy

from myflopy.sources import (
    DataSourceSpec,
    LiteralSource,
    source_from_dict,
)

Builder = Callable[..., Any]
Hook = Callable[[Any, dict[str, Any], "ModelContext"], Any]
Mf6Simulation: TypeAlias = flopy.mf6.MFSimulation
GwfModel: TypeAlias = flopy.mf6.ModflowGwf
GwtModel: TypeAlias = flopy.mf6.ModflowGwt
GweModel: TypeAlias = flopy.mf6.ModflowGwe
PrtModel: TypeAlias = flopy.mf6.ModflowPrt
Mf6Model: TypeAlias = GwfModel | GwtModel | GweModel | PrtModel
_ModelT = TypeVar("_ModelT", bound=Mf6Model)


@dataclass(frozen=True, slots=True)
class PackageRef:
    """A by-key placeholder for a package defined once in a project's library.

    Lets a model's package list point at a shared :class:`PackageSpec` registered
    on a :class:`~myflopy.workspace.Project` (via ``project.add_package(key,
    spec)``) instead of inlining it, so several models can reuse the same package
    definition. The reference is resolved to the real spec when the run is built.
    Construct one with the :func:`ref` helper rather than directly. It is
    serializable (``to_dict``/``from_dict``) so a recipe carrying refs round-trips.

    Attributes
    ----------
    key
        The project package-library key this reference resolves to.

    See Also
    --------
    ref : The preferred constructor.
    GridRef : The grid-library equivalent.
    """

    key: str

    def __post_init__(self) -> None:
        """Strip the key and reject an empty one (on this frozen dataclass)."""

        key = self.key.strip()
        if not key:
            raise ValueError("PackageRef key cannot be empty.")
        object.__setattr__(self, "key", key)

    @property
    def name(self) -> str:
        """Return the ref's spec-time identity: the full project key.

        The ref does not claim a package type -- the package it resolves to keeps
        its own ``PackageSpec.name``. The key prefix is used only as a
        slot-targeting hint by ``replace_package`` (see :func:`_package_slot`).
        """

        return self.key

    @property
    def enabled(self) -> bool:
        """Package references are active until resolved or removed."""

        return True

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready package reference."""

        return {"kind": "PackageRef", "key": self.key}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PackageRef:
        """Recreate a package reference from :meth:`to_dict` output."""

        return cls(data["key"])


def ref(key: str) -> PackageRef:
    """Return a deferred reference to a project-level package spec by ``key``.

    Use this in a model's package list to reuse a package that is defined once
    in a :class:`~myflopy.workspace.Project` package library, for example
    ``mf.gwf("flow", packages=[mf.ref("npf/base"), ...])``. The reference is
    resolved against the project's package library when the run is built.
    """

    return PackageRef(key)


@dataclass(frozen=True, slots=True)
class GridRef:
    """A by-key placeholder for a grid recipe/built grid in a project's library.

    The grid-library counterpart to :class:`PackageRef`: lets a model use a grid
    registered once on a :class:`~myflopy.workspace.Project` (via
    ``project.add_grid(key, ...)``) -- either an unbuilt :class:`GridSpec` recipe
    or an already-built grid -- instead of carrying its own. Resolved against the
    project's grid library when the run is built. Construct one with the
    :func:`grid_ref` helper. Serializable so recipes round-trip.

    Attributes
    ----------
    key
        The project grid-library key this reference resolves to.

    See Also
    --------
    grid_ref : The preferred constructor.
    PackageRef : The package-library equivalent.
    """

    key: str

    def __post_init__(self) -> None:
        """Strip the key and reject an empty one (on this frozen dataclass)."""

        key = self.key.strip()
        if not key:
            raise ValueError("GridRef key cannot be empty.")
        object.__setattr__(self, "key", key)

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready grid reference."""

        return {"kind": "GridRef", "key": self.key}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> GridRef:
        """Recreate a grid reference from :meth:`to_dict` output."""

        return cls(data["key"])


def grid_ref(key: str) -> GridRef:
    """Return a deferred reference to a project-level grid by ``key``.

    Use this as a model's ``grid`` to reuse a grid defined once in a
    :class:`~myflopy.workspace.Project` grid library, for example
    ``mf.gwf("flow", grid=mf.grid_ref("base"))``. The reference is resolved
    against the project's grid library when the run is built.
    """

    return GridRef(key)


class ModelType(str, Enum):
    """The MODFLOW 6 model kinds a :class:`ModelSpec` can build.

    A string enum naming the four supported model types; each maps to its FloPy
    constructor (``GWF`` -> ``ModflowGwf``, and likewise GWT/GWE/PRT). It is a
    ``str`` subclass, so ``ModelType.GWF == "gwf"`` and the member can be used
    anywhere the lowercase string is expected. The package-first helpers
    (``mf.gwf``/``mf.gwt``/``mf.gwe``/``mf.prt``) set this for you.

    Members
    -------
    GWF
        Groundwater flow.
    GWT
        Groundwater solute transport.
    GWE
        Groundwater energy/heat transport.
    PRT
        Particle tracking.
    """

    GWF = "gwf"
    GWT = "gwt"
    GWE = "gwe"
    PRT = "prt"


_DEFAULT_MODEL_BUILDERS: dict[ModelType, Builder] = {
    ModelType.GWF: flopy.mf6.ModflowGwf,
    ModelType.GWT: flopy.mf6.ModflowGwt,
    ModelType.GWE: flopy.mf6.ModflowGwe,
    ModelType.PRT: flopy.mf6.ModflowPrt,
}


def _require_unique_names(items: Iterable[Any], *, item_type: str) -> None:
    """Raise when a collection contains duplicate ``name`` values."""

    names = [item.name for item in items]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate {item_type} names: {', '.join(duplicates)}")


def _builder_label(builder: Builder | None) -> str:
    """Return a short readable label for a builder callable."""

    if builder is None:
        return "default"
    module = getattr(builder, "__module__", "")
    name = getattr(builder, "__name__", builder.__class__.__name__)
    if module.startswith("flopy."):
        return f"flopy.{name}"
    return name


def _callable_ref(builder: Builder) -> str:
    """Serialize a builder callable to a ``"module:qualname"`` reference.

    Raises if the callable is not importable (e.g. a local/lambda), since such a
    reference could not be resolved back later.
    """

    module = getattr(builder, "__module__", None)
    qualname = getattr(builder, "__qualname__", None)
    if not module or not qualname or "<locals>" in qualname:
        raise ValueError(
            f"Builder {builder!r} cannot be serialized. Use an importable builder."
        )
    return f"{module}:{qualname}"


def _resolve_callable(reference: str) -> Builder:
    """Import and return the callable named by a ``"module:qualname"`` reference."""

    module_name, _, qualname = reference.partition(":")
    if not module_name or not qualname:
        raise ValueError(f"Invalid builder reference: {reference!r}")
    value: Any = import_module(module_name)
    for part in qualname.split("."):
        value = getattr(value, part)
    return cast(Builder, value)


def _json_value(value: Any) -> Any:
    """Recursively convert a spec field to a JSON-serializable value.

    Sources serialize via ``to_dict``, paths to POSIX strings, tuples/lists/dicts
    element-wise; primitives pass through and anything else raises.
    """

    if isinstance(value, (DataSourceSpec, LiteralSource)):
        return value.to_dict()
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise ValueError(
        f"Value of type {type(value).__module__}.{type(value).__name__} "
        "cannot be serialized in a spec."
    )


def _spec_value(value: Any) -> Any:
    """Recursively rehydrate a deserialized value, rebuilding tagged source specs.

    The inverse of :func:`_json_value`: any ``{"kind": "...Source", ...}`` dict
    becomes a source spec; lists/dicts are walked; everything else passes through.
    """

    if isinstance(value, dict) and "kind" in value:
        kind = value["kind"]
        if kind.endswith("Source"):
            return source_from_dict(value)
    if isinstance(value, list):
        return [_spec_value(item) for item in value]
    if isinstance(value, dict):
        return {key: _spec_value(item) for key, item in value.items()}
    return value


def _preview(values: Iterable[Any], *, limit: int = 4) -> str:
    """Return a short comma-separated preview of an iterable."""

    items = list(values)
    shown = ", ".join(repr(value) for value in items[:limit])
    if len(items) > limit:
        shown = f"{shown}, ..."
    return shown


def _shape_label(value: Any) -> str | None:
    """Return a compact shape label for array/table-like values."""

    shape = getattr(value, "shape", None)
    if shape is None:
        return None
    class_name = value.__class__.__name__
    return f"{class_name}(shape={tuple(shape)!r})"


def _summarize_mapping(value: dict[Any, Any]) -> str:
    """Return a compact summary for mappings, including stress-period data."""

    keys = list(value.keys())
    if not keys:
        return "dict(len=0)"
    if all(isinstance(key, int) for key in keys):
        record_counts = [
            len(records)
            for records in value.values()
            if isinstance(records, (list, tuple))
        ]
        if record_counts:
            return (
                "period_data("
                f"periods={len(keys)}, "
                f"records={sum(record_counts)}, "
                f"max_per_period={max(record_counts)}"
                ")"
            )
    return f"dict(len={len(keys)}, keys=[{_preview(keys)}])"


def _summarize_value(value: Any) -> str:
    """Return a compact, notebook-friendly value summary."""

    shape = _shape_label(value)
    if shape is not None:
        return shape
    if isinstance(value, dict):
        return _summarize_mapping(value)
    if isinstance(value, (list, tuple)):
        if not value:
            return f"{type(value).__name__}(len=0)"
        return f"{type(value).__name__}(len={len(value)}, preview=[{_preview(value)}])"
    if isinstance(value, str):
        return repr(value if len(value) <= 80 else f"{value[:77]}...")
    text = repr(value)
    return text if len(text) <= 80 else f"{text[:77]}..."


def _html_table(rows: Iterable[tuple[str, str]], *, empty: str = "None") -> str:
    """Return a small HTML table for notebook rich display."""

    rows = list(rows)
    if not rows:
        return f"<p>{escape(empty)}</p>"
    body = "\n".join(
        f"<tr><td><code>{escape(key)}</code></td><td>{escape(value)}</td></tr>"
        for key, value in rows
    )
    return (
        "<table>"
        "<thead><tr><th>Name</th><th>Summary</th></tr></thead>"
        f"<tbody>{body}</tbody>"
        "</table>"
    )


def _require_model_type(model: Any, model_type: type[_ModelT], label: str) -> _ModelT:
    """Return ``model`` as ``model_type`` or explain the mismatch clearly."""

    if not isinstance(model, model_type):
        raise TypeError(
            f"Built model is not a {label}. "
            f"Found {model.__class__.__module__}.{model.__class__.__name__}."
        )
    return cast(_ModelT, model)


@dataclass(frozen=True, slots=True)
class PackageSpec:
    """Reusable instructions for building one package on a model or simulation.

    Builders receive the target as their first positional argument followed by
    the values in ``options``. A builder can be a FloPy package class or a
    small project-specific function.
    """

    name: str
    builder: Builder
    options: dict[str, Any] = field(default_factory=dict)
    enabled: bool = True
    requires: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)
    inputs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze ``requires`` to a tuple for hashable, immutable dependency ordering."""

        object.__setattr__(self, "requires", tuple(self.requires))

    def build(self, target: Any) -> Any | None:
        """Build the package on ``target``, or return ``None`` when disabled."""

        if not self.enabled:
            return None
        return self.builder(target, **self.inputs, **self.options)

    def with_options(self, **overrides: Any) -> PackageSpec:
        """Return a copy with selected options added or replaced."""

        return replace(self, options={**self.options, **overrides})

    def disabled(self) -> PackageSpec:
        """Return a disabled copy of this package specification."""

        return replace(self, enabled=False)

    def with_metadata(self, **updates: Any) -> PackageSpec:
        """Return a copy with selected provenance metadata added or replaced."""

        return replace(self, metadata={**self.metadata, **updates})

    def with_input(self, name: str, source: Any) -> PackageSpec:
        """Return a copy with one named package input source added."""

        return replace(self, inputs={**self.inputs, name: source})

    def with_inputs(self, **sources: Any) -> PackageSpec:
        """Return a copy with named package input sources added or replaced."""

        return replace(self, inputs={**self.inputs, **sources})

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready representation of this package recipe."""

        return {
            "kind": "PackageSpec",
            "name": self.name,
            "builder": _callable_ref(self.builder),
            "inputs": _json_value(self.inputs),
            "options": _json_value(self.options),
            "enabled": self.enabled,
            "requires": list(self.requires),
            "metadata": _json_value(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PackageSpec:
        """Recreate a package specification from :meth:`to_dict` output."""

        return cls(
            name=data["name"],
            builder=_resolve_callable(data["builder"]),
            inputs=_spec_value(dict(data.get("inputs", {}))),
            options=_spec_value(dict(data.get("options", {}))),
            enabled=bool(data.get("enabled", True)),
            requires=tuple(data.get("requires", ())),
            metadata=_spec_value(dict(data.get("metadata", {}))),
        )

    @property
    def option_summary(self) -> dict[str, str]:
        """Return a compact summary of package options for display."""

        return {key: _summarize_value(value) for key, value in self.options.items()}

    def __repr__(self) -> str:
        """Compact representation: name, builder, status, inputs, and option summary."""

        status = "enabled" if self.enabled else "disabled"
        pieces = [
            f"name={self.name!r}",
            f"builder={_builder_label(self.builder)!r}",
            f"status={status!r}",
            f"inputs={tuple(self.inputs)!r}",
            f"options={self.option_summary!r}",
        ]
        if self.requires:
            pieces.append(f"requires={self.requires!r}")
        if self.metadata:
            pieces.append(f"metadata_keys={tuple(self.metadata)!r}")
        return f"PackageSpec({', '.join(pieces)})"

    def _repr_html_(self) -> str:
        """Return a compact rich representation for notebooks."""

        metadata = ""
        if self.metadata:
            metadata = (
                "<p><strong>Metadata:</strong> "
                f"{escape(', '.join(self.metadata.keys()))}</p>"
            )
        requires = ""
        if self.requires:
            requires = (
                f"<p><strong>Requires:</strong> {escape(', '.join(self.requires))}</p>"
            )
        return (
            "<div>"
            f"<h4>PackageSpec: <code>{escape(self.name)}</code></h4>"
            f"<p><strong>Builder:</strong> {escape(_builder_label(self.builder))} "
            f"<strong>Status:</strong> {'enabled' if self.enabled else 'disabled'}</p>"
            f"{requires}"
            f"{metadata}"
            f"{_html_table(self.option_summary.items(), empty='No options')}"
            "</div>"
        )


PackageEntry: TypeAlias = PackageSpec | PackageRef


def _package_entry(value: PackageEntry | str) -> PackageEntry:
    """Normalize a package entry: a bare string becomes a :class:`PackageRef`."""

    if isinstance(value, str):
        return PackageRef(value)
    return value


def _package_entry_from_dict(data: dict[str, Any]) -> PackageEntry:
    """Rebuild a package entry from a payload, dispatching on its ``kind`` tag."""

    kind = data.get("kind")
    if kind == "PackageRef":
        return PackageRef.from_dict(data)
    if kind == "PackageSpec":
        return PackageSpec.from_dict(data)
    raise ValueError(f"Unknown package entry kind: {kind!r}")


def _package_entry_label(package: PackageEntry) -> str:
    """The display identity of an entry: a ref's key, else a concrete package's name."""

    if isinstance(package, PackageRef):
        return package.key
    return package.name


def _package_slot(package: PackageEntry) -> str:
    """Return the MF6 package slot an entry targets.

    A concrete package targets its own name. A reference targets the first
    segment of its key as a *hint* (``mf.ref("npf/high")`` -> ``"npf"``), so a
    ref can replace a same-typed package without being resolved first. When the
    key prefix is not the package type (a semantic key like ``"k/calibrated"``),
    pass an explicit ``name`` to :meth:`SimulationSpec.replace_package`.
    """

    if isinstance(package, PackageRef):
        return package.key.split("/", 1)[0]
    return package.name


def _package_entry_matches(package: PackageEntry, name: str) -> bool:
    """Match an entry by its full identity (name / key) or by its slot hint.

    So ``model.package("npf/base")`` finds that exact ref, and
    ``model.package("npf")`` finds whatever occupies the ``npf`` slot.
    """

    return package.name == name or _package_slot(package) == name


def _package_entry_is_enabled(package: PackageEntry) -> bool:
    """Whether an entry is enabled -- concrete packages carry the flag; refs are always on."""

    return package.enabled if isinstance(package, PackageSpec) else True


def _validate_concrete_packages(packages: Iterable[PackageSpec]) -> None:
    """Validate a resolved package list: unique names, satisfied ``requires``, correct order.

    Raises if names collide, an enabled package requires a missing/disabled one,
    or a required dependency is declared after the package that needs it.
    """

    packages = tuple(packages)
    _require_unique_names(packages, item_type="package")
    available = {package.name for package in packages if package.enabled}
    package_positions = {
        package.name: index for index, package in enumerate(packages) if package.enabled
    }
    for package in packages:
        missing = sorted(set(package.requires) - available)
        if package.enabled and missing:
            raise ValueError(
                f"Package '{package.name}' requires missing packages: {', '.join(missing)}"
            )
        if package.enabled:
            later = sorted(
                name
                for name in package.requires
                if package_positions[name] > package_positions[package.name]
            )
            if later:
                raise ValueError(
                    f"Package '{package.name}' must be declared after required packages: "
                    f"{', '.join(later)}"
                )


@dataclass(frozen=True, slots=True)
class GridSpec:
    """A deferred recipe for the grid a model is built on.

    Rather than a concrete grid object, a ``GridSpec`` describes *how* to build one
    -- from GeoPackage boundary/refinement/breakline layers (``GridSpec.voronoi``),
    a Python builder script (``GridSpec.python``), or a serialized dict
    (``GridSpec.from_dict``). It is resolved into a real grid at run time (or
    eagerly with ``.resolve(workspace)``).

    Use it two ways:

    - **Deferred** -- ``mf.gwf(...).with_grid(GridSpec.voronoi(...))``; the
      :class:`~myflopy.workspace.Project` builds the grid into
      ``run.workspace/_grid/<model>`` and populates ``model.context.grid``. Composes
      with disv + simple BCs.
    - **Eager** -- ``vor = GridSpec.voronoi(...).resolve(ws)``; then put ``vor`` in a
      ``ModelContext`` so the GIS package helpers (``mf.uzf``/``mf.sfr``/``.gpkg``)
      can resolve cells immediately.

    Examples
    --------
    >>> grid = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec("in.gpkg", layer="boundary"),
    ...                            refinement=mf.GeoPackageSourceSpec("in.gpkg", layer="refine"),
    ...                            boundary_max_area=200.0)
    >>> vor = grid.resolve("runs/_grid")     # eager: a concrete VoronoiGridPlus
    """

    name: str
    grid_type: str
    method: str
    script: Path | str | None = None
    function: str | None = None
    source: DataSourceSpec | None = None
    boundary: DataSourceSpec | None = None
    refinement: DataSourceSpec | None = None
    breaklines: tuple[DataSourceSpec, ...] = ()
    points: tuple[DataSourceSpec, ...] = ()
    inputs: tuple[Any, ...] = ()
    crs: str | None = None
    engine: str | None = None
    options: dict[str, Any] = field(default_factory=dict)
    triangle_options: dict[str, Any] = field(default_factory=dict)
    mesh_options: dict[str, Any] = field(default_factory=dict)
    voronoi_options: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    obj: Any = field(default=None, repr=False)
    persist: str = "pickle"
    pickle_path: str | None = None

    def __post_init__(self) -> None:
        """Freeze the source sequences (breaklines/points/inputs) to tuples for immutability."""

        object.__setattr__(self, "breaklines", tuple(self.breaklines))
        object.__setattr__(self, "points", tuple(self.points))
        object.__setattr__(self, "inputs", tuple(self.inputs))

    @classmethod
    def from_object(
        cls,
        obj: Any,
        *,
        name: str = "grid",
        grid_type: str = "disv",
        persist: str = "pickle",
    ) -> GridSpec:
        """Wrap an already-built grid object (such as a ``VoronoiGridPlus``).

        Build the grid however you like (typically with ``TriangleGrid`` /
        ``VoronoiGridPlus``), look at it, then register it. The object is held
        in memory; persisting it to disk is handled by the project, not here.
        """

        return cls(name=name, grid_type=grid_type, method="object", obj=obj, persist=persist)

    @classmethod
    def from_pickle(
        cls,
        path: Path | str,
        *,
        name: str = "grid",
        grid_type: str = "disv",
    ) -> GridSpec:
        """Reference a previously pickled grid on disk.

        The pickle is loaded lazily by :meth:`resolve`. Relative paths are
        resolved against the ``project_root`` supplied at resolve time.
        """

        return cls(name=name, grid_type=grid_type, method="pickle", pickle_path=str(path))

    @classmethod
    def python(
        cls,
        script: Path | str,
        *,
        function: str,
        name: str = "grid",
        grid_type: str = "disv",
        inputs: Iterable[Any] = (),
        crs: str | None = None,
        options: dict[str, Any] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> GridSpec:
        """Return a grid spec backed by a project-local Python builder script.

        The escape hatch for grids you build with your own code: the spec records
        the ``script`` path and the ``function`` name to call; :meth:`resolve`
        imports the script and invokes it to produce the grid object.

        Parameters
        ----------
        script : Path or str
            Path to a Python file containing the grid-builder function.
        function : str
            Name of the callable in ``script`` that returns the grid (required).
        name : str, default "grid"
            Grid name.
        grid_type : str, default "disv"
            Discretization type the builder produces.
        inputs : Iterable, optional
            Input files the builder needs (validated to exist at resolve time).
        crs : str, optional
            CRS metadata for the produced grid.
        options : dict, optional
            Keyword options forwarded to the builder function.
        metadata : dict, optional
            Free-form provenance metadata.

        Returns
        -------
        GridSpec

        Examples
        --------
        >>> mf.GridSpec.python("build_grid.py", function="make_grid",
        ...                    inputs=[mf.ShapeSource("domain.shp")])
        """

        if not function:
            raise ValueError("GridSpec.python requires a function name.")
        return cls(
            name=name,
            grid_type=grid_type,
            method="python",
            script=script,
            function=function,
            inputs=tuple(inputs),
            crs=crs,
            options={} if options is None else dict(options),
            metadata={} if metadata is None else dict(metadata),
        )

    @classmethod
    def from_geopackage(
        cls,
        path: Path | str,
        *,
        name: str = "grid",
        layer: str | None = None,
        id_column: str | None = None,
        grid_type: str = "disv",
        crs: str | None = None,
        **options: Any,
    ) -> GridSpec:
        """Load an existing grid from a GeoPackage. **Not implemented yet.**

        Planned API. :meth:`resolve` currently wires only generated Voronoi
        grids (:meth:`voronoi`) and Python grid builders (:meth:`python`), so a
        spec created here would raise at resolve time -- far from this call.
        This constructor therefore fails fast. To use an existing grid today,
        wrap a built grid object with :meth:`from_object`, or load a pickled
        grid with :meth:`from_pickle`.
        """

        raise NotImplementedError(
            "GridSpec.from_geopackage(...) is not implemented yet: myflopy "
            "resolves only generated Voronoi grids (GridSpec.voronoi) and Python "
            "grid builders (GridSpec.python). To use an existing grid today, wrap "
            "a built grid object with GridSpec.from_object(grid), or load one with "
            "GridSpec.from_pickle(path)."
        )

    @classmethod
    def structured(
        cls,
        *,
        name: str = "grid",
        nlay: int,
        nrow: int,
        ncol: int,
        delr: float,
        delc: float,
        top: Any,
        botm: Any,
        crs: str | None = None,
        **options: Any,
    ) -> GridSpec:
        """Build a simple structured DIS grid. **Not implemented yet.**

        Planned API. myflopy is Voronoi/DISV-first: :meth:`resolve` wires only
        generated Voronoi grids (:meth:`voronoi`) and Python grid builders
        (:meth:`python`), and the downstream package/visualization stack assumes
        an unstructured grid view. A spec created here would raise at resolve
        time, so this constructor fails fast instead. To use a structured grid
        today, build it yourself (e.g. via FloPy) and wrap it with
        :meth:`from_object`, or generate the geometry in a :meth:`python` builder.
        """

        raise NotImplementedError(
            "GridSpec.structured(...) is not implemented yet: myflopy resolves "
            "only generated Voronoi grids (GridSpec.voronoi) and Python grid "
            "builders (GridSpec.python). To use a structured grid today, build it "
            "yourself and wrap it with GridSpec.from_object(grid), or construct "
            "the geometry in a GridSpec.python(...) builder."
        )

    @classmethod
    def voronoi(
        cls,
        *,
        name: str = "grid",
        boundary: DataSourceSpec,
        refinement: DataSourceSpec | None = None,
        breaklines: Iterable[DataSourceSpec] = (),
        points: Iterable[DataSourceSpec] = (),
        crs: str | None = None,
        engine: str = "triangle_voronoi_plus",
        options: dict[str, Any] | None = None,
        triangle_options: dict[str, Any] | None = None,
        mesh_options: dict[str, Any] | None = None,
        voronoi_options: dict[str, Any] | None = None,
        **engine_options: Any,
    ) -> GridSpec:
        """Return a durable spec for a generated DISV/Voronoi grid.

        Records *where* the grid geometry comes from (boundary + refinement /
        breakline / point sources) and *how* to mesh it, so the grid can be built
        later with :meth:`resolve` (or by a :class:`Project` at run time). The
        engine builds a :class:`TriangleGrid` and wraps it in a
        :class:`VoronoiGridPlus`.

        Parameters
        ----------
        boundary : DataSourceSpec
            The model-domain polygon source (required).
        name : str, default "grid"
            Grid name.
        refinement : DataSourceSpec, optional
            Polygon source whose features refine the mesh (smaller cells inside).
        breaklines : Iterable[DataSourceSpec], optional
            Line sources buffered into refinement regions (e.g. streams, faults).
        points : Iterable[DataSourceSpec], optional
            Point sources pinned as fixed mesh vertices.
        crs : str, optional
            Target CRS for the grid (sources are reprojected to it).
        engine : str, default "triangle_voronoi_plus"
            Meshing engine identifier.
        options, triangle_options, mesh_options, voronoi_options : dict, optional
            Advanced per-stage option dicts forwarded to the Triangle / mesh /
            Voronoi stages.
        **engine_options
            Convenience options routed to the right stage by name (e.g.
            ``min_angle``, ``maximum_area``, ``profile``, ``idomain``).

        Returns
        -------
        GridSpec
            A durable grid recipe; call ``.resolve(workspace)`` for an eager grid,
            or pass it to ``mf.gwf(...).with_grid(...)`` / ``project.add_grid(...)``.

        Examples
        --------
        >>> spec = mf.GridSpec.voronoi(boundary=mf.ShapeSource("domain.shp"),
        ...                            refinement=mf.ShapeSource("wellfield.shp"),
        ...                            maximum_area=1.0e5, min_angle=30)
        >>> vor = spec.resolve(workspace="runs/_grid")
        """

        merged_options = {**({} if options is None else options)}
        merged_mesh_options = {**({} if mesh_options is None else mesh_options)}
        merged_triangle_options = {
            **({} if triangle_options is None else triangle_options)
        }
        merged_voronoi_options = {
            **({} if voronoi_options is None else voronoi_options)
        }
        for key, value in engine_options.items():
            if key == "min_angle":
                merged_triangle_options["angle"] = value
            elif key in {
                "angle",
                "exe_name",
                "maximum_area",
                "nodes",
                "additional_args",
                "region_point_tolerance",
            }:
                merged_triangle_options[key] = value
            elif key in {
                "cleanup",
                "damping",
                "max_optimization_points",
                "min_feature_area",
                "min_move",
                "optimization_iterations",
                "optimize",
                "profile",
                "protect_sources",
                "protected_labels",
                "resample_domain_boundary",
                "resample_region_sources",
                "simplify_tolerance",
                "snap_tolerance",
                "target_segment_length",
                "verbose",
            }:
                merged_mesh_options[key] = value
            elif key in {"idomain", "idomain_path", "name", "qhull_options", "rasters"}:
                merged_voronoi_options[key] = value
            else:
                merged_options[key] = value
        return cls(
            name=name,
            grid_type="disv",
            method="voronoi",
            boundary=boundary,
            refinement=refinement,
            breaklines=tuple(breaklines),
            points=tuple(points),
            crs=crs,
            engine=engine,
            options=merged_options,
            triangle_options=merged_triangle_options,
            mesh_options=merged_mesh_options,
            voronoi_options=merged_voronoi_options,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready representation of this grid recipe."""

        if self.method == "object":
            raise ValueError(
                "A GridSpec wrapping a built grid object cannot be serialized "
                "directly. Persist it with project.save(), or use "
                "GridSpec.python(...) for a reproducible recipe."
            )
        payload: dict[str, Any] = {
            "kind": "GridSpec",
            "name": self.name,
            "grid_type": self.grid_type,
            "method": self.method,
            "options": _json_value(self.options),
        }
        if self.pickle_path is not None:
            payload["pickle_path"] = self.pickle_path
        if self.script is not None:
            payload["script"] = _json_value(self.script)
        if self.function is not None:
            payload["function"] = self.function
        if self.inputs:
            payload["inputs"] = _json_value(self.inputs)
        if self.triangle_options:
            payload["triangle_options"] = _json_value(self.triangle_options)
        if self.mesh_options:
            payload["mesh_options"] = _json_value(self.mesh_options)
        if self.voronoi_options:
            payload["voronoi_options"] = _json_value(self.voronoi_options)
        for key in ("source", "boundary", "refinement"):
            value = getattr(self, key)
            if value is not None:
                payload[key] = value.to_dict()
        if self.breaklines:
            payload["breaklines"] = [source.to_dict() for source in self.breaklines]
        if self.points:
            payload["points"] = [source.to_dict() for source in self.points]
        if self.crs is not None:
            payload["crs"] = self.crs
        if self.engine is not None:
            payload["engine"] = self.engine
        if self.metadata:
            payload["metadata"] = _json_value(self.metadata)
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> GridSpec:
        """Recreate a grid specification from :meth:`to_dict` output."""

        return cls(
            name=data["name"],
            grid_type=data["grid_type"],
            method=data["method"],
            pickle_path=data.get("pickle_path"),
            script=data.get("script"),
            function=data.get("function"),
            source=(
                source_from_dict(data["source"])
                if data.get("source") is not None
                else None
            ),
            boundary=(
                source_from_dict(data["boundary"])
                if data.get("boundary") is not None
                else None
            ),
            refinement=(
                source_from_dict(data["refinement"])
                if data.get("refinement") is not None
                else None
            ),
            breaklines=tuple(
                source_from_dict(item) for item in data.get("breaklines", ())
            ),
            points=tuple(source_from_dict(item) for item in data.get("points", ())),
            inputs=tuple(_spec_value(list(data.get("inputs", ())))),
            crs=data.get("crs"),
            engine=data.get("engine"),
            options=_spec_value(dict(data.get("options", {}))),
            triangle_options=_spec_value(dict(data.get("triangle_options", {}))),
            mesh_options=_spec_value(dict(data.get("mesh_options", {}))),
            voronoi_options=_spec_value(dict(data.get("voronoi_options", {}))),
            metadata=_spec_value(dict(data.get("metadata", {}))),
        )

    def resolve(
        self,
        *,
        project_root: Path | str | None = None,
        workspace: Path | str | None = None,
        build: bool = True,
        return_triangle: bool = False,
    ) -> Any:
        """Resolve this grid recipe into the configured grid implementation.

        Generated Voronoi specs currently resolve through the existing
        ``TriangleGrid`` plus ``VoronoiGridPlus`` workflow. Set ``build=False``
        to prepare and inspect the Triangle setup without running Triangle.
        A spec created with :meth:`from_object` returns its held grid directly,
        and one created with :meth:`from_pickle` unpickles it from disk.
        """

        if self.method == "object":
            if self.obj is None:
                raise ValueError("This GridSpec has no materialized grid object.")
            return self.obj

        if self.method == "pickle":
            if self.pickle_path is None:
                raise ValueError("This GridSpec has no pickle path to load.")
            path = Path(self.pickle_path)
            if not path.is_absolute() and project_root is not None:
                path = Path(project_root) / path
            with path.open("rb") as handle:
                return pickle.load(handle)

        from myflopy.grid_spec_resolver import resolve_grid_spec

        return resolve_grid_spec(
            self,
            project_root=project_root,
            workspace=workspace,
            build=build,
            return_triangle=return_triangle,
        )


def _grid_entry(value: Any) -> GridSpec | GridRef | None:
    """Normalize a model/library grid value into a spec or reference.

    Accepts a :class:`GridSpec`, a :class:`GridRef`, ``None``, or an
    already-built grid object (such as a ``VoronoiGridPlus``), wrapping a bare
    object with :meth:`GridSpec.from_object`.
    """

    if value is None or isinstance(value, (GridSpec, GridRef)):
        return value
    return GridSpec.from_object(value)


def _grid_from_dict(data: dict[str, Any] | None) -> GridSpec | GridRef | None:
    """Rebuild a model grid from a serialized GridSpec or GridRef payload."""

    if data is None:
        return None
    if data.get("kind") == "GridRef":
        return GridRef.from_dict(data)
    return GridSpec.from_dict(data)


# MF6 Manning/weir conversion factors. LENGTH_CONVERSION scales a length in METERS to
# the model's length units; TIME_CONVERSION scales a time in SECONDS to the model's time
# units. Used by SFR (and LAK outlets) so the streamflow equations run in model units.
_MF6_LENGTH_TO_MODEL = {
    "meters": 1.0, "meter": 1.0, "m": 1.0,
    "feet": 3.28081, "foot": 3.28081, "ft": 3.28081,
    "centimeters": 100.0, "centimeter": 100.0, "cm": 100.0,
}
_MF6_TIME_TO_MODEL = {
    "seconds": 1.0, "second": 1.0, "sec": 1.0, "s": 1.0,
    "minutes": 60.0, "minute": 60.0, "min": 60.0,
    "hours": 3600.0, "hour": 3600.0, "hr": 3600.0, "h": 3600.0,
    "days": 86400.0, "day": 86400.0, "d": 86400.0,
    "years": 31557600.0, "year": 31557600.0, "yr": 31557600.0,
}


def mf6_length_conversion(units: str | None) -> float:
    """Return MF6 ``LENGTH_CONVERSION`` (meters -> model units) for a unit name.

    ``feet`` -> 3.28081, ``meters`` -> 1.0, ``centimeters`` -> 100.0. ``None`` or
    ``"unknown"`` -> 1.0. Raises ``ValueError`` for an unrecognized unit.
    """
    if units is None:
        return 1.0
    key = str(units).strip().lower()
    if key in ("", "unknown"):
        return 1.0
    try:
        return _MF6_LENGTH_TO_MODEL[key]
    except KeyError as error:
        raise ValueError(
            f"Unknown length unit {units!r}; expected feet, meters, or centimeters."
        ) from error


def mf6_time_conversion(units: str | None) -> float:
    """Return MF6 ``TIME_CONVERSION`` (seconds -> model units) for a unit name.

    ``days`` -> 86400, ``hours`` -> 3600, ``minutes`` -> 60, ``seconds`` -> 1.0,
    ``years`` -> 31557600. ``None`` or ``"unknown"`` -> 1.0. Raises for unrecognized.
    """
    if units is None:
        return 1.0
    key = str(units).strip().lower()
    if key in ("", "unknown"):
        return 1.0
    try:
        return _MF6_TIME_TO_MODEL[key]
    except KeyError as error:
        raise ValueError(
            f"Unknown time unit {units!r}; expected seconds, minutes, hours, days, or years."
        ) from error


@dataclass(frozen=True, slots=True)
class ModelContext:
    """The geometry a model is built against -- the bridge to the GIS builders.

    Context rides on the **model** (``mf.gwf(name, context=ctx, ...)``), not on the
    project, and is kept separate from FloPy package options on purpose: it gives
    the data-aware package builders and post-build hooks access to the grid,
    surfaces, dates and active domain. The built model exposes it as
    ``model.myflopy_context``.

    The GIS package helpers (``mf.uzf``/``mf.sfr``/``mf.lak`` and the ``.gpkg``
    boundary forms) take ``context=`` so they can map features/cells onto the grid.
    Because they resolve cells eagerly, build the grid first and put it here.

    Attributes
    ----------
    grid
        The grid object (e.g. a ``VoronoiGridPlus``) features are mapped onto.
    surfaces
        Per-cell layer elevations (``vor.gdf_topbtm`` or equivalent) read by
        surface-aware builders (SFR reach tops, LAK lake-cell layering).
    domain
        The active-domain (idomain) array used to pick/validate cells.
    dates
        Stress-period datetimes, when time-aware builders need them.
    metadata
        Free-form project metadata passed through to hooks.

    Examples
    --------
    >>> vor = mf.GridSpec.voronoi(...).resolve("runs/_grid")     # or a VoronoiGridPlus
    >>> layers = stack.build(attach=True)
    >>> ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)
    >>> flow = mf.gwf("flow", context=ctx, packages=[...])
    """

    grid: Any = None
    surfaces: Any = None
    dates: Any = None
    domain: Any = None
    length_units: str = "feet"
    time_units: str = "days"
    metadata: dict[str, Any] = field(default_factory=dict)

    def with_metadata(self, **updates: Any) -> ModelContext:
        """Return a copy with selected metadata added or replaced."""

        return replace(self, metadata={**self.metadata, **updates})


@dataclass(frozen=True, slots=True)
class SpecBuildContext:
    """Filesystem + library context threaded through a spec build.

    Carries the information a :class:`SimulationSpec`/:class:`ModelSpec` needs while
    it materializes into FloPy objects: where to write (project root, simulation
    workspace, grid workspace) and the resolved project libraries used to look up
    :class:`PackageRef`/:class:`GridRef` placeholders. The
    :class:`~myflopy.workspace.Run` lifecycle constructs and passes this for you;
    you rarely build one directly. ``build_grids`` toggles whether deferred grid
    recipes are materialized during the build.

    Attributes
    ----------
    project_root, simulation_workspace, grid_workspace
        Output locations (stored as :class:`~pathlib.Path` when set). When
        ``grid_workspace`` is unset, per-model grids land under
        ``<simulation_workspace>/_grid/<model>``.
    package_specs, grid_specs
        Resolved project libraries keyed by name, used to resolve refs.
    build_grids
        Whether to build deferred :class:`GridSpec` grids during this build.
    """

    project_root: Path | str | None = None
    simulation_workspace: Path | str | None = None
    grid_workspace: Path | str | None = None
    package_specs: dict[str, PackageSpec] = field(default_factory=dict)
    grid_specs: dict[str, GridSpec] = field(default_factory=dict)
    build_grids: bool = True

    def __post_init__(self) -> None:
        """Coerce the workspace fields to ``Path`` and copy the spec-library dicts."""

        for name in ("project_root", "simulation_workspace", "grid_workspace"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, Path(value))
        object.__setattr__(self, "package_specs", dict(self.package_specs))
        object.__setattr__(self, "grid_specs", dict(self.grid_specs))

    def with_simulation_workspace(
        self,
        workspace: Path | str | None,
    ) -> SpecBuildContext:
        """Return a copy with a simulation workspace when one is known."""

        if workspace is None or self.simulation_workspace is not None:
            return self
        return replace(self, simulation_workspace=Path(workspace))

    def model_grid_workspace(self, model_name: str) -> Path | None:
        """Return the generated grid workspace for one model, when available."""

        if self.grid_workspace is not None:
            return self.grid_workspace / model_name
        if self.simulation_workspace is not None:
            return self.simulation_workspace / "_grid" / model_name
        return None

    def package_spec(self, key: str) -> PackageSpec:
        """Return a resolved project-level package spec by key."""

        try:
            return self.package_specs[key]
        except KeyError as error:
            raise KeyError(
                f"Package reference '{key}' is unresolved in the build context."
            ) from error

    def grid_spec(self, key: str) -> GridSpec:
        """Return a resolved project-level grid spec by key."""

        try:
            return self.grid_specs[key]
        except KeyError as error:
            raise KeyError(
                f"Grid reference '{key}' is unresolved in the build context."
            ) from error


@dataclass(frozen=True, slots=True)
class PostBuildHook:
    """A named callback that runs once a model and all its packages are built.

    An extension point on :class:`ModelSpec`: after the FloPy model and every
    package have been created, each registered hook's ``callback`` is invoked with
    ``(model, packages, context)`` and its return value is stored under ``name`` in
    the :class:`BuiltModel`'s ``hook_results``. Use a hook to attach extra FloPy
    objects (e.g. an OBS utility package), tweak the assembled model, or compute a
    derived artifact -- without subclassing the builder.

    Attributes
    ----------
    name
        Key under which the hook's result is recorded in ``hook_results``.
    callback
        Callable ``(model, packages, context) -> Any`` run after the build.
    """

    name: str
    callback: Hook

    def run(self, model: Any, packages: dict[str, Any], context: ModelContext) -> Any:
        """Run the callback with the completed FloPy model and its context."""

        return self.callback(model, packages, context)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation of this hook.

        The ``callback`` must be an importable (module-level) function;
        lambdas and closures raise via :func:`_callable_ref`.
        """

        return {
            "name": self.name,
            "callback": _callable_ref(self.callback),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PostBuildHook:
        """Recreate a hook from :meth:`to_dict` output."""

        return cls(name=data["name"], callback=_resolve_callable(data["callback"]))


@dataclass(frozen=True, slots=True)
class BuiltModel:
    """The result of building one :class:`ModelSpec`: the FloPy model + its parts.

    Returned (inside a :class:`BuiltSimulation`) when a model spec is built. Bundles
    the live FloPy ``model`` with the package objects created from its spec
    (``packages``, keyed by package name), the :class:`ModelContext` it carried
    (geometry/dates), any ``hook_results`` from :class:`PostBuildHook`\\ s, and a
    back-reference to the parent ``simulation``. Convenience accessors such as
    ``.gwf`` / ``.gwt`` / ``.sim`` return the model under its concrete FloPy type.

    Attributes
    ----------
    model
        The built FloPy model object.
    packages
        Built package objects keyed by package name.
    context
        The :class:`ModelContext` (grid/domain/surfaces/dates) for this model.
    hook_results
        Values returned by post-build hooks, keyed by hook name.
    simulation
        The parent FloPy simulation, when known.
    """

    model: Mf6Model
    packages: dict[str, Any]
    context: ModelContext = field(default_factory=ModelContext)
    hook_results: dict[str, Any] = field(default_factory=dict)
    simulation: Mf6Simulation | None = None

    @property
    def sim(self) -> Mf6Simulation:
        """Return the parent FloPy simulation for this built model."""

        if self.simulation is not None:
            return self.simulation
        model_simulation = getattr(self.model, "sim", None)
        if isinstance(model_simulation, Mf6Simulation):
            return model_simulation
        raise AttributeError("This built model does not have a parent MFSimulation.")

    @property
    def gwf(self) -> GwfModel:
        """Return the model as a typed FloPy GWF model."""

        return self.as_gwf()

    @property
    def gwt(self) -> GwtModel:
        """Return the model as a typed FloPy GWT model."""

        return self.as_gwt()

    @property
    def gwe(self) -> GweModel:
        """Return the model as a typed FloPy GWE model."""

        return self.as_gwe()

    @property
    def prt(self) -> PrtModel:
        """Return the model as a typed FloPy PRT model."""

        return self.as_prt()

    def as_gwf(self) -> GwfModel:
        """Return the model as ``flopy.mf6.ModflowGwf``."""

        return _require_model_type(self.model, GwfModel, "GWF model")

    def as_gwt(self) -> GwtModel:
        """Return the model as ``flopy.mf6.ModflowGwt``."""

        return _require_model_type(self.model, GwtModel, "GWT model")

    def as_gwe(self) -> GweModel:
        """Return the model as ``flopy.mf6.ModflowGwe``."""

        return _require_model_type(self.model, GweModel, "GWE model")

    def as_prt(self) -> PrtModel:
        """Return the model as ``flopy.mf6.ModflowPrt``."""

        return _require_model_type(self.model, PrtModel, "PRT model")

    def package(self, name: str) -> Any:
        """Return one built package by name."""

        try:
            return self.packages[name]
        except KeyError as error:
            raise KeyError(f"Built model has no package named '{name}'.") from error


@dataclass(frozen=True, slots=True)
class ModelSpec:
    """Definition of one GWF, GWT, GWE, or PRT model.

    Package names are unique within the model. Calling :meth:`with_package`
    replaces an existing package with the same name, which makes scenario
    variants explicit and inexpensive to create.
    """

    name: str
    model_type: ModelType | str
    packages: tuple[PackageEntry, ...] = ()
    options: dict[str, Any] = field(default_factory=dict)
    builder: Builder | None = None
    context: ModelContext = field(default_factory=ModelContext)
    hooks: tuple[PostBuildHook, ...] = ()
    grid: GridSpec | GridRef | Any = None

    def __post_init__(self) -> None:
        """Normalize and validate the model spec's fields on construction.

        Coerces ``model_type`` to the enum, normalizes package/hook/grid entries,
        freezes them to tuples, and requires unique package + hook names.
        """

        object.__setattr__(self, "model_type", ModelType(self.model_type))
        object.__setattr__(
            self,
            "packages",
            tuple(_package_entry(package) for package in self.packages),
        )
        object.__setattr__(self, "hooks", tuple(self.hooks))
        object.__setattr__(self, "grid", _grid_entry(self.grid))
        _require_unique_names(self.packages, item_type="package")
        _require_unique_names(self.hooks, item_type="post-build hook")
        concrete_packages = [
            package for package in self.packages if isinstance(package, PackageSpec)
        ]
        if len(concrete_packages) == len(self.packages):
            _validate_concrete_packages(concrete_packages)

    def with_package(
        self, package: PackageEntry | str, *, slot: str | None = None
    ) -> ModelSpec:
        """Return a copy with ``package`` added or replaced in its slot.

        The slot defaults to the entry's own slot (a concrete package's name, or
        a reference's key prefix). Pass ``slot`` to target a slot explicitly.
        """

        package = _package_entry(package)
        target = slot if slot is not None else _package_slot(package)
        packages = list(self.packages)
        for index, existing in enumerate(packages):
            if _package_slot(existing) == target:
                packages[index] = package
                break
        else:
            packages.append(package)
        return replace(self, packages=tuple(packages))

    def with_packages(self, *packages: PackageEntry | str) -> ModelSpec:
        """Return a copy with each supplied package added or replaced by name."""

        spec = self
        for package in packages:
            spec = spec.with_package(package)
        return spec

    def package(self, name: str) -> PackageEntry:
        """Return one package specification by name."""

        for package in self.packages:
            if _package_entry_matches(package, name):
                return package
        raise KeyError(f"Model '{self.name}' has no package named '{name}'.")

    def without_package(self, name: str) -> ModelSpec:
        """Return a copy without the named package."""

        return replace(
            self,
            packages=tuple(
                package
                for package in self.packages
                if not _package_entry_matches(package, name)
            ),
        )

    def resolved_packages(
        self,
        build_context: SpecBuildContext | None = None,
    ) -> tuple[PackageSpec, ...]:
        """Return concrete package specs, resolving project package refs."""

        packages: list[PackageSpec] = []
        for package in self.packages:
            if isinstance(package, PackageSpec):
                packages.append(package)
                continue
            if build_context is None:
                raise ValueError(
                    f"Package reference '{package.key}' cannot be resolved without "
                    "a SpecBuildContext package library."
                )
            packages.append(build_context.package_spec(package.key))
        _validate_concrete_packages(packages)
        return tuple(packages)

    def with_options(self, **overrides: Any) -> ModelSpec:
        """Return a copy with selected model options added or replaced."""

        return replace(self, options={**self.options, **overrides})

    def with_context(self, context: ModelContext) -> ModelSpec:
        """Return a copy carrying the supplied model context."""

        return replace(self, context=context)

    def with_grid(self, grid: GridSpec | GridRef | Any) -> ModelSpec:
        """Return a copy with a grid recipe, reference, or built object attached."""

        return replace(self, grid=grid)

    def with_hook(self, hook: PostBuildHook) -> ModelSpec:
        """Return a copy with a post-build hook added or replaced by name."""

        hooks = list(self.hooks)
        for index, existing in enumerate(hooks):
            if existing.name == hook.name:
                hooks[index] = hook
                break
        else:
            hooks.append(hook)
        return replace(self, hooks=tuple(hooks))

    def build(
        self,
        simulation: Any,
        *,
        build_context: SpecBuildContext | None = None,
    ) -> BuiltModel:
        """Build the model and all enabled packages in declaration order."""

        context = self.context
        if self.grid is not None:
            build_context = (
                SpecBuildContext() if build_context is None else build_context
            )
            grid_spec = self.grid
            if isinstance(grid_spec, GridRef):
                grid_spec = build_context.grid_spec(grid_spec.key)
            grid = grid_spec.resolve(
                project_root=build_context.project_root,
                workspace=build_context.model_grid_workspace(self.name),
                build=build_context.build_grids,
            )
            context = replace(context, grid=grid)

        builder = self.builder or _DEFAULT_MODEL_BUILDERS[self.model_type]
        model = builder(simulation, modelname=self.name, **self.options)
        model.myflopy_context = context
        packages = self.resolved_packages(build_context)
        built_packages = {
            package.name: built
            for package in packages
            if (built := package.build(model)) is not None
        }
        hook_results = {
            hook.name: hook.run(model, built_packages, context) for hook in self.hooks
        }
        return BuiltModel(
            model=model,
            packages=built_packages,
            context=context,
            hook_results=hook_results,
            simulation=simulation,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready representation of this model recipe.

        Post-build hooks are serialized by importable reference; a hook whose
        callback is a lambda or closure raises via :func:`_callable_ref`.
        """

        payload: dict[str, Any] = {
            "kind": "ModelSpec",
            "name": self.name,
            "model_type": self.model_type.value,
            "options": _json_value(self.options),
            "packages": [package.to_dict() for package in self.packages],
            "context": {
                "metadata": _json_value(self.context.metadata),
            },
        }
        if self.builder is not None:
            payload["builder"] = _callable_ref(self.builder)
        if self.grid is not None:
            payload["grid"] = self.grid.to_dict()
        if self.hooks:
            payload["hooks"] = [hook.to_dict() for hook in self.hooks]
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ModelSpec:
        """Recreate a model specification from :meth:`to_dict` output."""

        context_data = data.get("context", {})
        return cls(
            name=data["name"],
            model_type=data["model_type"],
            packages=tuple(
                _package_entry_from_dict(package)
                for package in data.get("packages", ())
            ),
            options=_spec_value(dict(data.get("options", {}))),
            builder=(
                _resolve_callable(data["builder"])
                if data.get("builder") is not None
                else None
            ),
            context=ModelContext(
                metadata=_spec_value(dict(context_data.get("metadata", {}))),
            ),
            hooks=tuple(
                PostBuildHook.from_dict(hook) for hook in data.get("hooks", ())
            ),
            grid=_grid_from_dict(data.get("grid")),
        )

    @property
    def package_names(self) -> tuple[str, ...]:
        """Return enabled package names in declaration order."""

        return tuple(
            _package_entry_label(package)
            for package in self.packages
            if _package_entry_is_enabled(package)
        )

    def __repr__(self) -> str:
        """Compact representation: name, model type, package names, options, and hooks."""

        return (
            "ModelSpec("
            f"name={self.name!r}, "
            f"model_type={self.model_type.value!r}, "
            f"packages={self.package_names!r}, "
            f"options={_summarize_mapping(self.options)!r}, "
            f"hooks={tuple(hook.name for hook in self.hooks)!r}"
            ")"
        )

    def _repr_html_(self) -> str:
        """Return a compact rich representation for notebooks."""

        package_rows = [
            (
                _package_entry_label(package),
                (
                    "project package reference"
                    if isinstance(package, PackageRef)
                    else f"{_builder_label(package.builder)}; "
                    f"{len(package.options)} options"
                    + ("; disabled" if not package.enabled else "")
                    + (
                        f"; requires {', '.join(package.requires)}"
                        if package.requires
                        else ""
                    )
                ),
            )
            for package in self.packages
        ]
        hooks = ", ".join(hook.name for hook in self.hooks) or "None"
        return (
            "<div>"
            f"<h4>ModelSpec: <code>{escape(self.name)}</code></h4>"
            f"<p><strong>Type:</strong> {escape(self.model_type.value)} "
            f"<strong>Packages:</strong> {len(self.packages)} "
            f"<strong>Hooks:</strong> {escape(hooks)}</p>"
            f"{_html_table(package_rows, empty='No packages')}"
            "</div>"
        )


@dataclass(frozen=True, slots=True)
class ExchangeSpec:
    """Definition of an exchange between two or more models.

    The exchange builder receives the simulation, a tuple containing the built
    model objects named by ``models``, and then the configured keyword options.
    """

    name: str
    builder: Builder
    models: tuple[str, ...]
    options: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze ``models`` to a tuple and require at least two participants."""

        object.__setattr__(self, "models", tuple(self.models))
        if len(self.models) < 2:
            raise ValueError("An exchange must reference at least two models.")

    def build(self, simulation: Any, built_models: dict[str, BuiltModel]) -> Any:
        """Build the exchange after resolving its model names."""

        missing = [name for name in self.models if name not in built_models]
        if missing:
            raise KeyError(
                f"Exchange '{self.name}' references unknown models: {', '.join(missing)}"
            )
        models = tuple(built_models[name].model for name in self.models)
        return self.builder(simulation, models, **self.options)

    def with_options(self, **overrides: Any) -> ExchangeSpec:
        """Return a copy with selected exchange options added or replaced."""

        return replace(self, options={**self.options, **overrides})

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation of this exchange.

        The ``builder`` must be an importable (module-level) callable;
        lambdas and closures raise via :func:`_callable_ref`.
        """

        return {
            "name": self.name,
            "builder": _callable_ref(self.builder),
            "models": list(self.models),
            "options": _json_value(self.options),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ExchangeSpec:
        """Recreate an exchange from :meth:`to_dict` output."""

        return cls(
            name=data["name"],
            builder=_resolve_callable(data["builder"]),
            models=tuple(data.get("models", ())),
            options=_spec_value(dict(data.get("options", {}))),
        )

    def __repr__(self) -> str:
        """Compact representation: name, builder, participating models, and options."""

        return (
            "ExchangeSpec("
            f"name={self.name!r}, "
            f"builder={_builder_label(self.builder)!r}, "
            f"models={self.models!r}, "
            f"options={_summarize_mapping(self.options)!r}"
            ")"
        )


@dataclass(frozen=True, slots=True)
class BuiltSimulation:
    """The result of building a :class:`SimulationSpec`: the FloPy sim + its parts.

    The top-level product of ``SimulationSpec.build(...)`` (and what a
    :class:`~myflopy.workspace.Run` exposes as ``run.built``). Bundles the live
    FloPy ``simulation`` with its built models (``models``, each a
    :class:`BuiltModel`, keyed by model name), the simulation-level packages
    (``packages`` -- TDIS/IMS), and the inter-model ``exchanges``. Use ``.sim`` to
    reach the FloPy ``MFSimulation`` for writing/running, or ``run.model(name)`` to
    reach a specific built model.

    Attributes
    ----------
    simulation
        The built FloPy ``MFSimulation``.
    models
        Built models keyed by model name.
    packages
        Simulation-level packages (TDIS, IMS) keyed by name.
    exchanges
        Inter-model exchange objects keyed by name.
    """

    simulation: Mf6Simulation
    models: dict[str, BuiltModel]
    packages: dict[str, Any]
    exchanges: dict[str, Any]

    @property
    def sim(self) -> Mf6Simulation:
        """Return the built FloPy ``MFSimulation``."""

        return self.simulation

    def model(self, name: str) -> Mf6Model:
        """Return one built FloPy model by name."""

        return self.built_model(name).model

    def built_model(self, name: str) -> BuiltModel:
        """Return one ``BuiltModel`` wrapper by name."""

        try:
            return self.models[name]
        except KeyError as error:
            raise KeyError(f"Built simulation has no model named '{name}'.") from error

    def gwf(self, name: str) -> GwfModel:
        """Return a named model as ``flopy.mf6.ModflowGwf``."""

        return self.built_model(name).as_gwf()

    def gwt(self, name: str) -> GwtModel:
        """Return a named model as ``flopy.mf6.ModflowGwt``."""

        return self.built_model(name).as_gwt()

    def gwe(self, name: str) -> GweModel:
        """Return a named model as ``flopy.mf6.ModflowGwe``."""

        return self.built_model(name).as_gwe()

    def prt(self, name: str) -> PrtModel:
        """Return a named model as ``flopy.mf6.ModflowPrt``."""

        return self.built_model(name).as_prt()

    def package(self, name: str) -> Any:
        """Return one built simulation-level package by name."""

        try:
            return self.packages[name]
        except KeyError as error:
            raise KeyError(
                f"Built simulation has no simulation-level package named '{name}'."
            ) from error


@dataclass(frozen=True, slots=True)
class SimulationSpec:
    """One complete MF6 simulation: timing, solver(s), and the model(s).

    The top-level declarative object. It holds the simulation-wide packages
    (``mf.tdis`` and one or more ``mf.ims`` solvers) plus the ``ModelSpec`` models
    (built with ``mf.gwf``/``mf.gwt``/``mf.gwe``/``mf.prt``) and any inter-model
    ``exchanges``. ``build_flopy(workspace)`` materializes it into a FloPy
    simulation; or register it on a :class:`~myflopy.workspace.Project` with
    ``project.add_simulation(sim)`` and run it via ``project.prepare_run(...)``.

    Parameters
    ----------
    name
        Simulation name (and default run name).
    models
        The model specs (``mf.gwf(...)`` etc.).
    packages
        Simulation-wide packages -- ``mf.tdis(...)`` and ``mf.ims(...)`` solver(s).
    exchanges
        Inter-model exchanges (e.g. ``mf.build_gwf_gwt_exchange("flow", "transport")``).
    workspace, executable, run_name
        Optional run location / MF6 executable / run name overrides.

    Examples
    --------
    >>> sim = mf.SimulationSpec("baseline", models=(flow,),
    ...     packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
    ...               mf.ims(models=("valley",), complexity="MODERATE")])
    >>> built = sim.build_flopy("runs/baseline")     # -> FloPy simulation
    """

    name: str
    models: tuple[ModelSpec, ...]
    packages: tuple[PackageEntry, ...] = ()
    exchanges: tuple[ExchangeSpec, ...] = ()
    options: dict[str, Any] = field(default_factory=dict)
    builder: Builder = flopy.mf6.MFSimulation
    workspace: Path | str | None = None
    executable: str = "mf6"
    run_name: str | None = None
    derived_from: str | None = None
    lineage: tuple[dict[str, Any], ...] = ()

    def __post_init__(self) -> None:
        """Normalize and validate the simulation spec's fields on construction.

        Freezes models/packages/exchanges/lineage to tuples (normalizing package
        entries) and requires unique model and simulation-package names.
        """

        object.__setattr__(self, "models", tuple(self.models))
        object.__setattr__(
            self,
            "packages",
            tuple(_package_entry(package) for package in self.packages),
        )
        object.__setattr__(self, "exchanges", tuple(self.exchanges))
        object.__setattr__(self, "lineage", tuple(self.lineage))
        _require_unique_names(self.models, item_type="model")
        _require_unique_names(self.packages, item_type="simulation package")
        _require_unique_names(self.exchanges, item_type="exchange")
        concrete_packages = [
            package for package in self.packages if isinstance(package, PackageSpec)
        ]
        if len(concrete_packages) == len(self.packages):
            _validate_concrete_packages(concrete_packages)

    def with_model(self, model: ModelSpec) -> SimulationSpec:
        """Return a copy with ``model`` added or replaced by name."""

        models = list(self.models)
        for index, existing in enumerate(models):
            if existing.name == model.name:
                models[index] = model
                break
        else:
            models.append(model)
        return replace(self, models=tuple(models))

    def with_models(self, *models: ModelSpec) -> SimulationSpec:
        """Return a copy with each supplied model added or replaced by name."""

        spec = self
        for model in models:
            spec = spec.with_model(model)
        return spec

    def model(self, name: str) -> ModelSpec:
        """Return one model specification by name."""

        for model in self.models:
            if model.name == name:
                return model
        raise KeyError(f"Simulation '{self.name}' has no model named '{name}'.")

    def with_package(self, package: PackageEntry | str) -> SimulationSpec:
        """Return a copy with a simulation-level package added or replaced."""

        package = _package_entry(package)
        packages = list(self.packages)
        for index, existing in enumerate(packages):
            if existing.name == package.name:
                packages[index] = package
                break
        else:
            packages.append(package)
        return replace(self, packages=tuple(packages))

    def with_packages(self, *packages: PackageEntry | str) -> SimulationSpec:
        """Return a copy with each simulation-level package added or replaced."""

        spec = self
        for package in packages:
            spec = spec.with_package(package)
        return spec

    def with_exchange(self, exchange: ExchangeSpec) -> SimulationSpec:
        """Return a copy with ``exchange`` added or replaced by name."""

        exchanges = list(self.exchanges)
        for index, existing in enumerate(exchanges):
            if existing.name == exchange.name:
                exchanges[index] = exchange
                break
        else:
            exchanges.append(exchange)
        return replace(self, exchanges=tuple(exchanges))

    def with_workspace(self, workspace: Path | str) -> SimulationSpec:
        """Return a copy with a default run workspace."""

        return replace(self, workspace=workspace)

    def derive(self, name: str) -> SimulationSpec:
        """Return a named version of this fully resolved simulation recipe."""

        return replace(
            self,
            name=name,
            run_name=name,
            derived_from=self.name,
            lineage=(
                *self.lineage,
                {"operation": "derive", "from": self.name, "to": name},
            ),
        )

    def replace_package(
        self, model_name: str, package: PackageEntry | str, *, name: str | None = None
    ) -> SimulationSpec:
        """Return a copy with a model package replaced in its slot.

        The slot defaults to the replacement's own slot (a concrete package's
        name, or a reference's key prefix). Pass ``name`` to target a slot
        explicitly -- needed when a semantically-keyed ref such as
        ``mf.ref("k/calibrated")`` replaces the ``npf`` slot.

        The target slot must already exist; use :meth:`add_package` to add.
        """

        package = _package_entry(package)
        model = self.model(model_name)
        target = name if name is not None else _package_slot(package)
        if not any(_package_slot(existing) == target for existing in model.packages):
            raise KeyError(
                f"Model '{model_name}' has no '{target}' package to replace. "
                "Pass name=<slot> to target a slot explicitly, or use add_package."
            )
        return replace(
            self.with_model(model.with_package(package, slot=target)),
            lineage=(
                *self.lineage,
                {
                    "operation": "replace_package",
                    "model": model_name,
                    "package": _package_entry_label(package),
                },
            ),
        )

    def replace_grid(
        self, model_name: str, grid: GridSpec | GridRef | Any
    ) -> SimulationSpec:
        """Return a copy with one model's grid replaced.

        ``grid`` may be a :class:`GridSpec`, a :func:`grid_ref` reference into a
        project grid library, or an already-built grid object.
        """

        model = self.model(model_name)
        if isinstance(grid, GridRef):
            grid_label = grid.key
        elif isinstance(grid, GridSpec):
            grid_label = grid.name
        else:
            grid_label = "object"
        return replace(
            self.with_model(model.with_grid(grid)),
            lineage=(
                *self.lineage,
                {"operation": "replace_grid", "model": model_name, "grid": grid_label},
            ),
        )

    def add_package(
        self, model_name: str, package: PackageEntry | str, *, name: str | None = None
    ) -> SimulationSpec:
        """Return a copy with a package added or replaced on one model."""

        package = _package_entry(package)
        model = self.model(model_name)
        target = name if name is not None else _package_slot(package)
        existed = any(_package_slot(existing) == target for existing in model.packages)
        operation = "replace_package" if existed else "add_package"
        return replace(
            self.with_model(model.with_package(package, slot=target)),
            lineage=(
                *self.lineage,
                {
                    "operation": operation,
                    "model": model_name,
                    "package": _package_entry_label(package),
                },
            ),
        )

    def resolved_packages(
        self,
        build_context: SpecBuildContext | None = None,
    ) -> tuple[PackageSpec, ...]:
        """Return concrete simulation-level package specs."""

        packages: list[PackageSpec] = []
        for package in self.packages:
            if isinstance(package, PackageSpec):
                packages.append(package)
                continue
            if build_context is None:
                raise ValueError(
                    f"Package reference '{package.key}' cannot be resolved without "
                    "a SpecBuildContext package library."
                )
            packages.append(build_context.package_spec(package.key))
        _validate_concrete_packages(packages)
        return tuple(packages)

    def build(
        self,
        workspace: Path | str | None = None,
        *,
        run_name: str | None = None,
        build_context: SpecBuildContext | None = None,
    ):
        """Build this specification as a workflow ``Run``.

        Use :meth:`build_flopy` when you need only the raw in-memory FloPy
        objects without a run lifecycle wrapper.
        """

        from myflopy.workspace import Run

        run_workspace = workspace if workspace is not None else self.workspace
        if run_workspace is None:
            run_workspace = Path(self.run_name or self.name)
        run = Run(
            name=run_name or self.run_name or self.name,
            workspace=run_workspace,
            spec=self,
            executable=self.executable,
            build_context=build_context,
        )
        run.build()
        return run

    def build_flopy(
        self,
        workspace: Path | str | None = None,
        *,
        build_context: SpecBuildContext | None = None,
    ) -> BuiltSimulation:
        """Build and return the raw FloPy simulation objects."""

        options = dict(self.options)
        run_workspace = workspace if workspace is not None else self.workspace
        if run_workspace is not None:
            options["sim_ws"] = str(run_workspace)
        simulation = self.builder(sim_name=self.name, **options)
        build_context = (build_context or SpecBuildContext()).with_simulation_workspace(
            run_workspace
        )

        packages = self.resolved_packages(build_context)
        # MF6 build order: timing (TDIS) before models, then exchanges, then
        # solvers (IMS). Solver packages register their models, so the models
        # must already exist -- otherwise multiple solutions collapse into one,
        # which MF6 rejects for coupled (e.g. GWF-GWT) simulations.
        timing = [p for p in packages if p.builder is flopy.mf6.ModflowTdis]
        solvers = [p for p in packages if p.builder is not flopy.mf6.ModflowTdis]

        simulation_packages: dict[str, Any] = {}
        for package in timing:
            if (built := package.build(simulation)) is not None:
                simulation_packages[package.name] = built
        built_models = {
            model.name: model.build(simulation, build_context=build_context)
            for model in self.models
        }
        built_exchanges = {
            exchange.name: exchange.build(simulation, built_models)
            for exchange in self.exchanges
        }
        for package in solvers:
            if (built := package.build(simulation)) is not None:
                simulation_packages[package.name] = built
        return BuiltSimulation(
            simulation=simulation,
            models=built_models,
            packages=simulation_packages,
            exchanges=built_exchanges,
        )

    @property
    def model_names(self) -> tuple[str, ...]:
        """Return model names in declaration order."""

        return tuple(model.name for model in self.models)

    @property
    def package_names(self) -> tuple[str, ...]:
        """Return simulation-level package names in declaration order."""

        return tuple(
            _package_entry_label(package)
            for package in self.packages
            if _package_entry_is_enabled(package)
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready representation of this simulation recipe.

        Inter-model exchanges are serialized by importable builder reference; an
        exchange whose builder is a lambda or closure raises via
        :func:`_callable_ref`.
        """

        payload: dict[str, Any] = {
            "kind": "SimulationSpec",
            "name": self.name,
            "models": [model.to_dict() for model in self.models],
            "packages": [package.to_dict() for package in self.packages],
            "exchanges": [exchange.to_dict() for exchange in self.exchanges],
            "options": _json_value(self.options),
            "builder": _callable_ref(self.builder),
            "executable": self.executable,
            "lineage": _json_value(list(self.lineage)),
        }
        if self.workspace is not None:
            payload["workspace"] = str(Path(self.workspace).as_posix())
        if self.run_name is not None:
            payload["run_name"] = self.run_name
        if self.derived_from is not None:
            payload["derived_from"] = self.derived_from
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SimulationSpec:
        """Recreate a simulation specification from :meth:`to_dict` output."""

        return cls(
            name=data["name"],
            models=tuple(
                ModelSpec.from_dict(model) for model in data.get("models", ())
            ),
            packages=tuple(
                _package_entry_from_dict(package)
                for package in data.get("packages", ())
            ),
            exchanges=tuple(
                ExchangeSpec.from_dict(exchange)
                for exchange in data.get("exchanges", ())
            ),
            options=_spec_value(dict(data.get("options", {}))),
            builder=(
                _resolve_callable(data["builder"])
                if data.get("builder") is not None
                else flopy.mf6.MFSimulation
            ),
            workspace=data.get("workspace"),
            executable=data.get("executable", "mf6"),
            run_name=data.get("run_name"),
            derived_from=data.get("derived_from"),
            lineage=tuple(_spec_value(list(data.get("lineage", ())))),
        )

    def __repr__(self) -> str:
        """Compact representation: name, model names, package names, and exchange names."""

        return (
            "SimulationSpec("
            f"name={self.name!r}, "
            f"models={self.model_names!r}, "
            f"packages={self.package_names!r}, "
            f"exchanges={tuple(exchange.name for exchange in self.exchanges)!r}"
            ")"
        )

    def _repr_html_(self) -> str:
        """Return a compact rich representation for notebooks."""

        model_rows = [
            (
                model.name,
                f"{model.model_type.value}; packages={len(model.packages)}",
            )
            for model in self.models
        ]
        package_rows = [
            (
                package.name,
                (
                    "project package reference"
                    if isinstance(package, PackageRef)
                    else f"{_builder_label(package.builder)}; {len(package.options)} options"
                ),
            )
            for package in self.packages
            if _package_entry_is_enabled(package)
        ]
        exchange_rows = [
            (
                exchange.name,
                f"{_builder_label(exchange.builder)}; models={', '.join(exchange.models)}",
            )
            for exchange in self.exchanges
        ]
        return (
            "<div>"
            f"<h4>SimulationSpec: <code>{escape(self.name)}</code></h4>"
            f"<p><strong>Models:</strong> {len(self.models)} "
            f"<strong>Simulation Packages:</strong> {len(self.packages)} "
            f"<strong>Exchanges:</strong> {len(self.exchanges)}</p>"
            "<h5>Models</h5>"
            f"{_html_table(model_rows, empty='No models')}"
            "<h5>Simulation Packages</h5>"
            f"{_html_table(package_rows, empty='No simulation packages')}"
            "<h5>Exchanges</h5>"
            f"{_html_table(exchange_rows, empty='No exchanges')}"
            "</div>"
        )
