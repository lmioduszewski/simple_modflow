"""Composable specifications for building MODFLOW 6 simulations.

The classes in this module are intentionally small. A spec stores a named
builder and its keyword arguments, while the builder remains ordinary Python.
This keeps model assembly explicit and makes package alternatives easy to
replace without introducing a framework around FloPy.
"""

from __future__ import annotations

import copy
import difflib
import functools
import hashlib
import inspect
import os
import pickle
import sys
import warnings
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from enum import Enum
from html import escape
from importlib import import_module
from numbers import Integral
from pathlib import Path, PurePath
from types import ModuleType
from typing import Any, TypeAlias, TypeVar, cast

import flopy

from myflopy._deprecation import deprecated_instance_getattr
from myflopy._length_units import (
    _check_mf6_length_token,
    check_model_length_units,
)
from myflopy._logging import get_logger
from myflopy._time_units import check_model_time_units, check_tdis_time_units
from myflopy.areas import _same_crs

# Re-exported so no caller changed when these moved out (ledger 264). Each name is
# marked unused-import-exempt individually: ruff cannot tell a re-export from a stray
# import, and the alternative -- a per-file ignore, as `myflopy/__init__.py` has --
# would switch the check off for all 3,500 lines of this module rather than for these
# twelve names. (Spelling the directive out in prose here would itself be read as one.)
from myflopy.grid_identity import (
    CELLWISE_OPTIONS,
    DROPPED_CONTEXT_FIELDS,  # noqa: F401
    GENERATED_SOLVER,
    GRID_PIN,  # noqa: F401
    GRID_PIN_OPTIONS,  # noqa: F401
    INLINED_FROM_REF,
    MOVE_RECIPE,
    MOVED_FROM,
    MOVED_VIA_REF,
    PIN_VOIDED,
    REMAPPED_FROM,
    SOURCE_SHA,  # noqa: F401
    _as_int,  # noqa: F401
    _build_refusal,
    _built_grid,
    _canonical_ring,  # noqa: F401
    _cellids_in,  # noqa: F401
    _context_metadata,  # noqa: F401
    _discretization_complaint_families,
    _discretization_complaints,  # noqa: F401
    _mesh_residual,  # noqa: F401
    _mesh_scale,  # noqa: F401
    active_cell_warnings,
    cellid_form_complaints,
    cellwise_options,
    coupled_mesh_complaints,
    disv_pin,  # noqa: F401
    dropped_context_fields,  # noqa: F401
    grid_pin,  # noqa: F401
    implicit_lake_models,
    implicit_solver_complaints,
    lak_build_warnings,
    move_report,
    moved_package_complaints,
    package_move,  # noqa: F401 -- re-exported: the public spelling is myflopy.specs
    package_moves,
    source_fingerprint,  # noqa: F401
    unwritable_faces,
    void_pin,
)
from myflopy.sources import (
    ArrayProduct,
    DataSourceSpec,
    Heads,
    LiteralSource,
    Polygons,
    RasterValues,
    is_source_payload,
    source_from_dict,
)
from myflopy.specs_io import simulation_from_yaml, simulation_to_yaml

logger = get_logger(__name__)

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

    Parameters
    ----------
    key : str
        The project package-library key this reference resolves to, as given to
        ``project.add_package(key, spec)``. Leading and trailing whitespace is
        stripped; an empty or all-whitespace key raises ``ValueError``. It is
        not checked against any library here -- an unknown key fails when the
        run is built and the reference is resolved.
    moved_from : str, optional
        The pin of the mesh a ``with_grid`` moved this reference AWAY from,
        set by that move when the library package it names could not be
        verified on the new mesh -- per-cell values, or cellids it typed rather
        than resolved. The build copies it onto the package it resolves and
        refuses it, as it refuses an inline package the move stamped. Not
        something to set by hand: replace the entry with a fresh ``mf.ref(key)``
        (or name the key in ``grid_free=``) to clear it. Default ``None``.

    See Also
    --------
    ref : The preferred constructor.
    GridRef : The grid-library equivalent.
    """

    key: str
    moved_from: str | None = None

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

    @property
    def metadata(self) -> dict[str, Any]:
        """The move stamp in the shape a :class:`PackageSpec` carries it, else empty.

        So every reader of a stamp -- :func:`~myflopy.grid_identity.package_move`'s
        report, the build's refusal -- reads a stamped ref as it reads a stamped
        package, without the library the ref names.
        """

        return {} if self.moved_from is None else {MOVED_FROM: self.moved_from}

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready package reference.

        ``moved_from`` is written only when set, so an unmoved reference
        serializes exactly as it did before the field existed.
        """

        payload: dict[str, Any] = {"kind": "PackageRef", "key": self.key}
        if self.moved_from is not None:
            payload["moved_from"] = self.moved_from
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PackageRef:
        """Recreate a package reference from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A mapping carrying a ``"key"`` entry, as :meth:`to_dict` writes
            (``{"kind": "PackageRef", "key": ...}``), and optionally
            ``"moved_from"``. A missing ``"key"`` raises ``KeyError``; other
            entries are ignored.
        """

        return cls(data["key"], moved_from=data.get("moved_from"))


def ref(key: str) -> PackageRef:
    """Return a deferred reference to a project-level package spec by ``key``.

    Use this in a model's package list to reuse a package that is defined once
    in a :class:`~myflopy.workspace.Project` package library, for example
    ``mf.gwf("flow", packages=[mf.ref("npf/base"), ...])``. The reference is
    resolved against the project's package library when the run is built.

    Parameters
    ----------
    key : str
        The package-library key, as registered with ``project.add_package(key,
        spec)``. Whitespace is stripped; an empty key raises ``ValueError``.
        The key is not looked up here, so an unregistered one fails only when
        the run is built.

    Returns
    -------
    PackageRef
        A reference whose ``name`` is ``key``.
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

    Parameters
    ----------
    key : str
        The project grid-library key this reference resolves to, as given to
        ``project.add_grid(key, ...)``. Leading and trailing whitespace is
        stripped; an empty or all-whitespace key raises ``ValueError``. It is
        not checked against any library here.

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
        """Recreate a grid reference from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A mapping carrying a ``"key"`` entry, as :meth:`to_dict` writes
            (``{"kind": "GridRef", "key": ...}``). Only ``"key"`` is read; a
            missing one raises ``KeyError``.
        """

        return cls(data["key"])

    def resolve(
        self,
        *,
        project_root: Path | str | None = None,
        workspace: Path | str | None = None,
        build: bool = True,
        return_triangle: bool = False,
    ) -> Any:
        """Always raises: a reference has nothing to resolve without its project.

        Present so that ``project.grids[key]`` -- typed ``GridSpec | GridRef``
        -- offers ``.resolve`` on both members, and so the failure names the fix
        instead of arriving as a bare ``AttributeError`` on a union an editor
        already flagged. A ``GridRef`` points INTO a project's grid library, so
        it cannot build anything on its own; look the key up, or let
        ``prepare_run`` resolve it.

        Parameters
        ----------
        project_root : Path or str, optional
            Accepted only so the signature matches :meth:`GridSpec.resolve`;
            ignored. Default ``None``.
        workspace : Path or str, optional
            Accepted to match :meth:`GridSpec.resolve`; ignored. Default ``None``.
        build : bool, optional
            Accepted to match :meth:`GridSpec.resolve`; ignored. Default ``True``.
        return_triangle : bool, optional
            Accepted to match :meth:`GridSpec.resolve`; ignored. Default
            ``False``.

        Raises
        ------
        TypeError
            Always, whatever the arguments, naming ``project.grids[key]`` and
            ``project.prepare_run(...)`` as the two ways to resolve it.
        """

        raise TypeError(
            f"GridRef({self.key!r}) is a by-key reference, not a recipe -- it "
            f"has no sources to build from. Resolve the entry it points at "
            f"(project.grids[{self.key!r}] holds the GridSpec), or let "
            f"project.prepare_run(...) resolve it into the model."
        )


def grid_ref(key: str) -> GridRef:
    """Return a deferred reference to a project-level grid by ``key``.

    Use this as a model's ``grid`` to reuse a grid defined once in a
    :class:`~myflopy.workspace.Project` grid library, for example
    ``mf.gwf("flow", grid=mf.grid_ref("base"))``. The reference is resolved
    against the project's grid library when the run is built.

    Parameters
    ----------
    key : str
        The grid-library key, as registered with ``project.add_grid(key, ...)``.
        Whitespace is stripped; an empty key raises ``ValueError``. The key is
        not looked up here, so an unregistered one fails only when the run is
        built.

    Returns
    -------
    GridRef
        A by-key reference, serialized on the model by ``ModelSpec.to_dict``.
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
    """Return a short readable label for a builder callable.

    A ``functools.partial`` is unwrapped to the FloPy class it binds -- the seven
    list BCs all build through ``partial(_build_named, ModflowGwfchd)`` and so
    every one of them otherwise labelled itself ``'partial'``, in every repr and
    summary.
    """

    if builder is None:
        return "default"
    if isinstance(builder, functools.partial):
        bound = next(
            (arg for arg in builder.args if isinstance(arg, type)),
            None,
        )
        return _builder_label(bound) if bound is not None else _builder_label(builder.func)
    module = getattr(builder, "__module__", "")
    name = getattr(builder, "__name__", builder.__class__.__name__)
    if module.startswith("flopy."):
        return f"flopy.{name}"
    return name


def _callable_ref(builder: Builder) -> Any:
    """Serialize a builder callable to a JSON-ready reference.

    A plain importable callable becomes a ``"module:qualname"`` string. A
    :func:`functools.partial` over an importable function becomes a dict
    ``{"partial": <func ref>, "args": [...], "keywords": {...}}`` -- this is the
    shape the list-BC ``*_spec`` factories use (``partial(_build_named,
    ModflowGwfchd)``), so CHD/GHB/DRN/RIV/WEL/RCH/EVT specs round-trip. Each bound
    argument is itself a reference (callables/classes, e.g. the FloPy package
    class) or a plain JSON value. Raises ``ValueError`` if any piece is not
    importable -- a lambda (module-level too), a local function, a name that no
    longer imports anything callable, or a method bound to an instance (or reached
    through a subclass) whose binding the reference would drop -- since such a
    reference could not be resolved later, and a saved project holding one would
    fail to load as a whole. A plain function or class is written by NAME: one
    rebound since (a notebook cell re-run) is written as what the name holds now.
    """

    if isinstance(builder, functools.partial):
        return {
            "partial": _callable_ref(builder.func),
            "args": [_bound_arg_ref(arg) for arg in builder.args],
            "keywords": {key: _bound_arg_ref(val) for key, val in builder.keywords.items()},
        }
    module = getattr(builder, "__module__", None)
    qualname = getattr(builder, "__qualname__", None)
    # `<` marks a name no import can reach: `<lambda>` (a MODULE-LEVEL lambda's
    # qualname is just that, with no `<locals>`), `<locals>`, `<genexpr>`. Checking
    # `<locals>` alone let `f = lambda ...` at module level through as
    # "mod:<lambda>", so `Project.save()` wrote a project `Project.load()` then died
    # on with an AttributeError.
    if not module or not qualname or "<" in qualname:
        raise ValueError(
            f"Builder {builder!r} cannot be serialized: {qualname or 'it'} has no "
            "importable name (a lambda or a function defined inside another one). "
            "Define it with `def` at module level."
        )
    resolved = _resolve_reference(module, qualname)
    if resolved is None or not callable(resolved):
        raise ValueError(
            f"Builder {builder!r} cannot be serialized: '{module}:{qualname}' "
            "imports nothing callable now (deleted, or rebound to a value), so a "
            "saved project could not load it back."
        )
    bound_to = getattr(builder, "__self__", None)
    if (
        bound_to is not None
        and not isinstance(bound_to, ModuleType)
        and resolved is not builder
        and resolved != builder
    ):
        # A METHOD names its function, not what it is bound to: an instance method
        # resolves to the plain function and a classmethod reached through a
        # subclass to the base class's binding, so the reference would DROP the
        # instance or the subclass. `==` admits a classmethod read off its own class,
        # a fresh bound object on every attribute read.
        raise ValueError(
            f"Builder {builder!r} cannot be serialized: '{module}:{qualname}' names "
            f"the function, not what it is bound to ({bound_to!r}), so the reference "
            "would load a different callable. Use a module-level function, a class, "
            "or a classmethod read off its own class."
        )
    # A plain function or class is a NAME, not an object (ledger 387): what
    # the name holds now is what a reload gets, exactly as for a module edited
    # between save and load -- `package_digest` likewise never sees a builder's
    # body. So a notebook cell re-run, which rebinds `hook` to a fresh but
    # identical function, still saves; refusing it dropped the simulation from the
    # save with nothing but a log line.
    return f"{module}:{qualname}"


def _resolve_reference(module_name: str, qualname: str) -> Any:
    """Return what ``module:qualname`` names now, or None when nothing is there."""

    module = sys.modules.get(module_name)
    if module is None:
        try:
            module = import_module(module_name)
        except ImportError:
            logger.debug(
                "Module %r is not importable; the builder reference is refused.",
                module_name,
            )
            return None
    value: Any = module
    for part in qualname.split("."):
        value = getattr(value, part, None)
        if value is None:
            return None
    return value


def _bound_arg_ref(value: Any) -> Any:
    """Serialize one ``functools.partial`` bound argument.

    A callable/class becomes a tagged reference (``{"$callable": <ref>}``) so it
    survives as an import path; anything else goes through :func:`_json_value`.
    """

    if callable(value):
        return {"$callable": _callable_ref(value)}
    return _json_value(value)


def _resolve_callable(reference: Any) -> Builder:
    """Import and return the callable named by a reference from :func:`_callable_ref`.

    Accepts a ``"module:qualname"`` string, or a ``{"partial": ...}`` dict that
    rebuilds a :func:`functools.partial`.
    """

    if isinstance(reference, dict):
        func = _resolve_callable(reference["partial"])
        args = tuple(_revive_bound_arg(arg) for arg in reference.get("args", ()))
        keywords = {key: _revive_bound_arg(val) for key, val in reference.get("keywords", {}).items()}
        return cast(Builder, functools.partial(func, *args, **keywords))
    module_name, _, qualname = reference.partition(":")
    if not module_name or not qualname:
        raise ValueError(f"Invalid builder reference: {reference!r}")
    value: Any = import_module(module_name)
    for part in qualname.split("."):
        value = getattr(value, part)
    return cast(Builder, value)


def _revive_bound_arg(value: Any) -> Any:
    """Inverse of :func:`_bound_arg_ref`."""

    if isinstance(value, dict) and set(value) == {"$callable"}:
        return _resolve_callable(value["$callable"])
    # A bound argument is builder input, so it gets the same int-key restoration
    # as ``options``/``inputs`` -- see :func:`_spec_options`.
    return _spec_value(value, int_keys=True)


def _json_value(value: Any) -> Any:
    """Recursively convert a spec field to a JSON-serializable value.

    Sources serialize via ``to_dict``, paths to POSIX strings, tuples/lists/dicts
    element-wise, numpy scalars to the Python type they wrap; primitives pass
    through and anything else raises.

    The output is a TRUE Python primitive, not merely something ``json.dumps``
    accepts -- ``to_yaml`` goes through PyYAML, which dispatches on the exact type
    and refuses a ``numpy.float64`` that JSON takes happily.
    """

    if (token := getattr(type(value), "_SPEC_TOKEN", None)) is not None:
        # A `$`-token, not a `kind` tag: `Surface` has a `kind` field of its own and
        # the source dispatcher keys on `kind`, so tagging one as the other made a
        # surface come back a plain dict (ledger 269). Duck-typed on the class
        # attribute rather than imported -- `specs` is import layer 1 and `layers`
        # is layer 6.
        return {token: _json_value(value.to_dict())}
    if isinstance(value, ArrayProduct):
        # Its own token, and NOT a `_SPEC_TOKEN`: its factors are SOURCES, so every
        # walker below treats them as sources (role "source", `external` honoured, the
        # package helpers stamp their hashes) rather than as a token payload whose
        # `path` keys would be repointed as a surface's (ledger 541).
        return {ArrayProduct.TOKEN: _json_value(value.to_dict())}
    if isinstance(value, (DataSourceSpec, LiteralSource)):
        # Recursed, not handed straight back. A source's own `to_dict` writes its
        # fields raw -- `metadata`, `table`, `crs`, `band`, `value` -- so this was the
        # one branch above the normalization that could still return numpy, and
        # `PackageSpec.inputs` is exactly where source specs live. Measured: a
        # `np.int64` inside one made `Project.validate()` return `[]` and
        # `Project.save()` die inside `json.dumps`, having already written the
        # project and simulation files, so the next `Project.load` raised
        # `FileNotFoundError` on the package it never got to (ledger 265).
        return _json_value(value.to_dict())
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, tuple):
        # Tagged, not flattened to a list. JSON has no tuple, and the difference
        # is load-bearing: FloPy re-tuples cellids when it builds a recarray, but
        # ``ModflowUtlobs`` writes each record field with ``repr``, so an
        # untagged cellid comes back a list and is written ``west chd [0, 0, 0]``
        # instead of ``west chd 1 1 1``. MF6 parses that, terminates normally,
        # and returns DNODATA for every observation. Unlike an integer key there
        # is no invertible read-side test -- a list that WAS a tuple is
        # indistinguishable from one that never was -- so the tag is the fix.
        return {"$tuple": [_json_value(item) for item in value]}
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if getattr(value, "ndim", None) == 0 and hasattr(value, "item"):
        # A numpy SCALAR, which is a number and always was. `.item()` is the scalar
        # protocol, so this covers int64, float32 and bool_ alike, and a 0-d array
        # with them.
        #
        # This must sit ABOVE the primitive check, not below it, and the first
        # version of it did not (ledger 263, corrected in 265). `np.float64`
        # SUBCLASSES `float`, so it satisfied `isinstance(value, float)` and was
        # returned unconverted, while `np.int64` does not subclass `int` and fell
        # through to the raise. Converting only the half that raised left the other
        # half in the output: `json.dumps` accepts a `np.float64`, so `to_dict` and
        # the frozen snapshot looked fine, but PyYAML dispatches on the EXACT type
        # and refuses it -- so `SimulationSpec.to_yaml()` raised `RepresenterError`
        # for every model carrying voronoi gridprops, on any grid, which is why the
        # real-grid fixture could not have found it either.
        unwrapped = value.item()
        if type(unwrapped) is not type(value):
            return _json_value(unwrapped)
        # `np.longdouble.item()` hands back a `np.longdouble`, so recursing on it
        # never terminates -- a `RecursionError` in place of this function's own
        # message naming the type. Fall through to that message instead.
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise ValueError(
        f"Value of type {type(value).__module__}.{type(value).__name__} "
        "cannot be serialized in a spec."
    )


def _is_int_key(key: Any) -> bool:
    """Whether ``key`` is exactly what :func:`_json_value` writes for an ``int`` key.

    The test is a canonical round trip -- ``str(int(key)) == key`` -- rather than
    "looks like digits", which makes the promotion in :func:`_spec_value`
    provably invertible and collision-free. ``"0"`` and ``"-1"`` qualify;
    ``"007"``, ``"1_0"``, ``"1.0"``, ``" 1"`` and ``"head.csv"`` do not, so a
    restored dict can never merge two distinct keys onto one.
    """

    if not isinstance(key, str):
        return False
    try:
        return str(int(key)) == key
    except ValueError:
        return False


def _spec_value(value: Any, *, int_keys: bool = False) -> Any:
    """Recursively rehydrate a deserialized value, rebuilding tagged source specs.

    The inverse of :func:`_json_value`: any ``{"kind": "...Source", ...}`` dict
    becomes a source spec; lists/dicts are walked; everything else passes through.

    ``int_keys`` restores the ``int`` dict keys that :func:`_json_value` had to
    stringify to reach JSON. MF6 period-keyed options (``stress_period_data``,
    ``saverecord``, ``steady_state``, ...) are keyed by stress period, and FloPy
    silently writes an EMPTY period block for a string key -- so without this a
    saved list BC comes back with no rows under a normal termination. Only whole
    dicts whose keys are *all* canonical integers are promoted; see
    :func:`_is_int_key`. Off by default because free-form fields (``metadata``,
    ``lineage``) must keep the keys the user gave them.
    """

    if isinstance(value, dict) and set(value) == {"$tuple"}:
        return tuple(_spec_value(item, int_keys=int_keys) for item in value["$tuple"])
    if isinstance(value, dict) and len(value) == 1:
        # Imported lazily: `myflopy.surfaces` is import layer 2 and `myflopy.layers`
        # is layer 6, so a module-level import here would be an upward edge.
        #
        # Every payload is `_untag_tuples`-ed first. `_json_value` tags a tuple
        # wherever it meets one, INSIDE a token's payload included, but this branch
        # handed the payload to `from_dict` raw -- so the tag was never removed. A
        # `split=(0.3, 0.7)` came back as `{"$tuple": [0.3, 0.7]}`, and the saved
        # stack died building with `could not convert string to float: '$tuple'`
        # while its digest no longer matched the live declaration's.
        if "$surface" in value:
            from myflopy.surfaces import Surface

            return Surface.from_dict(_untag_tuples(value["$surface"]))
        if "$layer_stack" in value:
            from myflopy.layers import LayerStack

            return LayerStack.from_dict(_untag_tuples(value["$layer_stack"]))
        if "$spread" in value:
            from myflopy.geopackage import Spread

            return Spread.from_dict(_untag_tuples(value["$spread"]))
        if "$cell_surface_offset" in value:
            from myflopy.geopackage import CellSurfaceOffset

            return CellSurfaceOffset.from_dict(_untag_tuples(value["$cell_surface_offset"]))
        if ArrayProduct.TOKEN in value:
            # Each factor through THIS function, so a source factor is revived exactly
            # as a bare one is -- `from_dict` keeps an already-live recipe as it is.
            payload = value[ArrayProduct.TOKEN]
            return ArrayProduct.from_dict(
                {"factors": [_spec_value(item) for item in payload["factors"]]}
            )
    if isinstance(value, dict) and "kind" in value:
        # The registry decides, not a name heuristic: ``DataSourceSpec``'s kind
        # ends in "Spec", so an ``endswith("Source")`` test handed the public
        # base class back as a plain dict. The endswith arm is kept so an
        # UNKNOWN "...Source" kind still raises rather than degrading silently.
        if is_source_payload(value) or str(value["kind"]).endswith("Source"):
            if int_keys:
                # A source wraps its payload, so the promotion has to reach
                # inside it -- ``LiteralSource({0: True})`` is period-keyed and
                # otherwise came back string-keyed through JSON while YAML kept
                # the ints, leaving the two formats disagreeing.
                value = {
                    key: _spec_value(item, int_keys=True)
                    for key, item in value.items()
                }
            return source_from_dict(value)
    if isinstance(value, list):
        return [_spec_value(item, int_keys=int_keys) for item in value]
    if isinstance(value, dict):
        restored = {
            key: _spec_value(item, int_keys=int_keys) for key, item in value.items()
        }
        if int_keys and all(map(_is_int_key, restored)):
            return {int(key): item for key, item in restored.items()}
        return restored
    return value


def _untag_tuples(payload: Any) -> Any:
    """Turn every ``{"$tuple": [...]}`` inside a token's payload back into a tuple.

    ONLY that. A token's payload is the object's OWN format -- ``LayerStack.to_dict``
    nests each surface as a plain ``{"kind": "flat", ...}`` dict and its ``from_dict``
    expects exactly that back -- so the rest of :func:`_spec_value` must not run on
    it: no source revival keyed on ``kind``, and no integer-key promotion, which is
    for period-keyed builder options and would rewrite a ``grass_kwargs`` key the
    user typed. Removing the one tag :func:`_json_value` adds is the exact inverse
    of what that function did to the payload.
    """

    if isinstance(payload, dict):
        if set(payload) == {"$tuple"}:
            return tuple(_untag_tuples(item) for item in payload["$tuple"])
        return {key: _untag_tuples(item) for key, item in payload.items()}
    if isinstance(payload, list):
        return [_untag_tuples(item) for item in payload]
    return payload


def _spec_options(data: Any) -> dict[str, Any]:
    """Rehydrate a builder-kwargs mapping, restoring ``int`` keys inside its values.

    The mapping's own keys are FloPy parameter names and are left alone; it is
    the dicts *nested* in a value (``stress_period_data={0: [...]}``) that carry
    stress-period keys. Used for :class:`PackageSpec`'s ``options``/``inputs``,
    the only spec fields that reach a package builder able to take period data,
    and for :class:`ExchangeSpec`'s options, because ``ModflowGwfgwf.perioddata``
    carries an exchange-level MVR. Grid, model and simulation options have no
    period dimension -- neither constructor takes a dict-valued parameter.
    """

    return {key: _spec_value(item, int_keys=True) for key, item in dict(data).items()}


#: Exact leaf types that can be neither a path, a source, nor a token-bearing object.
#: Membership is tested by EXACT type, not `isinstance`: a subclass may carry a
#: `_SPEC_TOKEN` or be a `PurePath`, and must fall through to the branches below.
_CANONICAL_LEAF_TYPES = frozenset((float, int, bool, type(None), bytes, complex))


def _canonical(value: Any) -> Any:
    """Normalize a value for fingerprinting, so representation is not content.

    Dict key ORDER is not content -- ``Project.save()``/``load()`` round-trips a
    package's options through JSON and hands them back in a different order, so
    without this every reloaded package fingerprints as edited. Sets are ordered
    by ``repr`` for the same reason one step further: their pickle varies BETWEEN
    PROCESSES, which would make a digest unstable in exactly the way its contract
    promises it is not.

    Lists and tuples keep their order (it is meaningful in a record list) but
    collapse to one type, since MF6 cannot tell a list of records from a tuple of
    them. Arrays and unknown objects pass through; four kinds of value do NOT,
    because the pickle would otherwise tell them apart by representation alone:

    - a NumPy scalar becomes its Python value -- ``np.float64`` subclasses float,
      so JSON keeps it and hands back a plain ``float``, which pickles differently;
    - a string is interned. A pickle records object IDENTITY: one object used
      twice is written once and back-referenced, two equal objects are written
      twice. A boundname a GeoPackage feature shares across its records, or OC's
      two ``"ALL"``, is one object before ``Project.save`` and several after, and
      every such package read as edited (ledger 222);
    - a ``PurePath`` becomes its posix string, which is what ``_json_value`` writes;
    - a token-bearing object (``Surface``, ``LayerStack``) or a source is
      fingerprinted by its ``to_dict()`` -- its DECLARATION.

    The last two exist because the digest's job is to survive the round trip
    ``_json_value`` defines, so the two functions have to agree about what a value
    IS. They disagreed for anything ``_json_value`` rewrote and this passed through,
    and ledger 222's defect recurred one level deeper inside such an object.

    Fingerprinting the declaration has a consequence worth knowing: whatever
    ``to_dict`` does not write is INVISIBLE. That cuts the right way for a cached
    preview grid (``LayerStack._draft``, which ``.plot`` populates -- previewing a
    stack used to mark its package edited) and the wrong way if a field ever stops
    being written, which is why the digest tests assert discrimination on every
    conditionally-written field, not just round-trip stability.
    """

    # Normalize and then FALL THROUGH, rather than returning. Returning here meant
    # `np.str_` -- which IS a `numpy.generic` -- took this branch and never reached
    # the interning below, so the very fix this function documents did not apply to
    # the strings a GeoPackage actually produces: `np.unique` on a boundname column
    # gives `np.str_`, and a package using one reported EDITED after a plain
    # save/load, through `project.runs.drift()`. Testing `ndim`/`item` rather than
    # `numpy.generic` also covers a 0-d ndarray, which `_json_value` normalizes and
    # this did not -- so the two disagreed across a save and moved the digest on
    # their own (ledger 265).
    if getattr(value, "ndim", None) == 0 and hasattr(value, "item"):
        value = value.item()
    if type(value) is str:
        return sys.intern(value)
    if type(value) in _CANONICAL_LEAF_TYPES:
        # Pure speed, no behaviour: an EXACT-type hit here can be neither a path, a
        # source nor a token-bearing object, so it skips every branch below. Tested by
        # exact type and not `isinstance`, because a subclass may be either. Measured
        # digesting the user's real 12-package library: 153.9 ms with it, 221.7 ms
        # without, against 182.9 ms before this commit -- so it more than pays for the
        # recursion the commit adds.
        return value

    # These three mirror `_json_value`'s FIRST three branches, in its order and ABOVE
    # the containers -- which is the whole point of the commit. The digest has to
    # survive the round trip `_json_value` defines, so the two must agree about what a
    # value IS. They disagreed for anything `_json_value` rewrites and `_canonical`
    # passed through, and a `LayerStack` in `inputs` is exactly that: the SAME stack
    # digested one way live and another after a reload, because `'constant_thickness'`
    # is one interned string object in the live graph and several after JSON, so pickle
    # back-references it once and writes it repeatedly (ledger 222, one level deeper --
    # inside an object this function never entered). A ONE-LAYER stack does not show
    # it. Order matters for a token-bearing Mapping, which the dict branch would
    # otherwise claim first, fingerprinting a body `_json_value` throws away.
    if (token := getattr(type(value), "_SPEC_TOKEN", None)) is not None:
        # Duck-typed on the class attribute rather than imported: `specs` is import
        # layer 1 and `layers` is layer 6. The tag is kept because `Surface.to_dict()`
        # is a bare `{'kind': ..., 'value': ...}` and `LayerStack.to_dict()` carries no
        # kind at all, so an untagged payload is indistinguishable from an ordinary
        # dict. It does NOT make the result unforgeable -- a literal tuple of the same
        # shape still collides, exactly as it does at `_json_value`.
        # `digest_payload()` where the class defines one, exactly as the source branch
        # below does. A `Surface` substitutes its file's content hash for the path and
        # drops the GRASS cache; a `LayerStack` reaches its surfaces' payloads, which
        # `to_dict` has already flattened into plain dicts. Without this the LAYER STACK
        # carrier kept a path-keyed digest while the `.gpkg` carrier beside it did not
        # -- so bundling a layered model moved its disv digest on a byte-identical
        # archive, which is the churn ledger 293 exists to remove (found in review).
        payload = getattr(value, "digest_payload", value.to_dict)()
        return (sys.intern(token), _canonical(payload))
    if isinstance(value, ArrayProduct):
        # Its factors as each one digests -- a raster or polygon map by its CONTENT
        # (ledger 293) -- in order: a product is not commutative in what it SAYS, and
        # a reordered declaration reading as edited is the safe direction.
        return (
            sys.intern(ArrayProduct.TOKEN),
            tuple(_canonical(factor) for factor in value.factors),
        )
    if isinstance(value, (DataSourceSpec, LiteralSource)):
        # A CONSTANT tag, not the class name: a source's own `to_dict` already carries
        # the `kind` its reader dispatches on, so renaming the Python class must not
        # move a digest.
        #
        # `digest_payload()`, NOT `to_dict()`. They differ for exactly one source --
        # a `BoundarySource`, which substitutes its file's content hash for its path,
        # so that bundling a project into an archive does not report every lazy
        # package as edited and an edited GeoPackage is not reported as unchanged
        # (ledger 293). Every other source returns `to_dict()` unchanged, so no
        # existing digest moves; `tests/test_boundary_source.py` pins both halves.
        return (sys.intern("$source"), _canonical(value.digest_payload()))
    if isinstance(value, PurePath):
        return sys.intern(value.as_posix())

    if isinstance(value, dict):
        return tuple(
            sorted(
                ((_canonical(key), _canonical(item)) for key, item in value.items()),
                key=lambda pair: repr(pair[0]),
            )
        )
    if isinstance(value, (set, frozenset)):
        return tuple(sorted((_canonical(item) for item in value), key=repr))
    if isinstance(value, (list, tuple)):
        return tuple(_canonical(item) for item in value)
    return value


#: Every ``Surface`` field that can name a file the stack READS. ``zone`` is typed
#: ``Any`` and takes EITHER a shapely geometry or a path to a vector file
#: (``Surface.where``/``within``, taught in the cheatsheet), which is why it was missed
#: on the first pass -- and missing it shipped an archive that reported success and
#: died `DataSourceError` on a file nobody had carried. The value-type test below is
#: what makes including it safe: a geometry is not a ``str``, so it is left alone.
#: Matched by KEY inside a
#: serialized payload rather than by attribute, because ``specs`` is import layer 1 and
#: ``surfaces``/``layers`` are 2 and 6 -- the same reason ``_canonical`` duck-types on
#: ``_SPEC_TOKEN`` instead of importing the classes. ``Spread`` and ``CellSurfaceOffset``,
#: the other two token-bearing payloads that reach ``inputs``, declare no field of any of
#: these names, so the key match cannot claim one of theirs.
SURFACE_INPUT_KEYS = frozenset(
    {"path", "contours", "region_raster", "region_vector", "zone"}
)

#: ``Surface.out`` -- the GRASS interpolation CACHE, which the stack WRITES rather than
#: reads. Separate because the two want opposite treatment when a file is absent: a
#: missing input is a fault to report, a missing cache is just one that has not been
#: produced yet.
SURFACE_WRITE_KEYS = frozenset({"out"})


#: The one payload key whose own KEYS are field names, so a naive walk reads its
#: VALUES -- content hashes -- as if they were the paths they stand in for. Found by
#: test: the bundler reported `proj/06464ba5355f3e8e` as a missing source.
SHA_KEY = "source_shas"


def _path_text(value: Any) -> str:
    """A path as it is WRITTEN in a payload: posix text, never a `Path` object."""

    return value.as_posix() if isinstance(value, PurePath) else str(value)


def _repoint_payload(payload: Any, transform, notes: list[tuple[str, str]]) -> Any:
    """Recurse through a serialized token payload, transforming the paths it names."""

    if isinstance(payload, dict):
        rebuilt = {}
        for key, item in payload.items():
            if key == SHA_KEY:
                rebuilt[key] = item
            elif isinstance(item, str) and key in SURFACE_INPUT_KEYS:
                stored = transform(PurePath(item), role="surface", notes=notes)
                rebuilt[key] = item if stored is None else _path_text(stored)
            elif isinstance(item, str) and key in SURFACE_WRITE_KEYS:
                stored = transform(PurePath(item), role="write_target", notes=notes)
                rebuilt[key] = item if stored is None else _path_text(stored)
            else:
                rebuilt[key] = _repoint_payload(item, transform, notes)
        return rebuilt
    if isinstance(payload, list):
        return [_repoint_payload(item, transform, notes) for item in payload]
    return payload


def repoint_source_paths(
    value: Any,
    transform,
    *,
    honor_external: bool = False,
    notes: list[tuple[str, str]] | None = None,
) -> Any:
    """Return ``value`` with every file path it names passed through ``transform``.

    ONE walker for the two callers that must agree about what a path IS -- the
    archive bundler (``project.save(bundle_sources=True)``) and the build-time anchor
    in :meth:`ModelSpec.build`. They were written apart in the first draft and
    disagreed immediately: bundling rewrote a layer stack's raster to a
    project-relative path that the anchor, which only knew about ``DataSourceSpec``,
    then left for ``rasterio`` to resolve against the current working directory. Half
    a feature is worse than neither, because the save reports success.

    **Two carriers, and the one that is easy to forget fails FIRST** (ledger 293). A
    lazy ``.gpkg`` recipe is a :class:`~myflopy.sources.DataSourceSpec` in ``inputs``;
    a lazy ``mf.disv`` holds a gridless ``LayerStack``, which is not a dataclass, names
    no source spec, and per ledger 270 is the ONLY shape a layered model saves in. A
    model carrying both dies on the raster before the boundary condition is reached,
    so bundling only the recipe ships an archive whose headline is false.

    Parameters
    ----------
    value : DataSourceSpec, ArrayProduct, Surface, LayerStack, Spread, CellSurfaceOffset, dict, list, tuple, or any other object
        What to walk. A :class:`~myflopy.sources.DataSourceSpec` has its ``path``
        transformed, and each further path field its class names in
        ``EXTRA_PATH_FIELDS`` (``Heads.grid_path``) -- each walked as a source. An :class:`~myflopy.sources.ArrayProduct` has each factor walked
        as a source, ``honor_external`` included, and is rebuilt when one moved. A
        token-bearing object (one with a ``_SPEC_TOKEN``) is written
        down with its own ``to_dict()``, every ``path``/``contours``/
        ``region_raster``/``region_vector``/``zone`` string and every ``out`` string
        in that payload is transformed, and the object is rebuilt from the result
        with its class's ``from_dict`` -- the ``{"$tuple": [...]}`` tag
        ``_json_value`` leaves in a payload taken off first, as a reload does. A
        ``dict``, ``list`` or ``tuple`` is walked and rebuilt as the same type.
        Anything else is returned as it is.
    transform : callable
        ``transform(path, *, role, notes)`` returns the path to store (a ``str`` or
        ``Path``), or ``None`` to leave the path alone. ``path`` is the source's own
        ``Path`` (or an ``EXTRA_PATH_FIELDS`` one) for a ``DataSourceSpec`` and a
        ``PurePath`` inside a token payload.
        ``role`` is ``"source"`` (a ``DataSourceSpec``), ``"surface"`` (a file a
        surface READS) or ``"write_target"`` (``Surface.out``, the GRASS cache it
        WRITES); see :data:`SURFACE_WRITE_KEYS` for why the last is apart. ``notes``
        is the list below, for the transform to append its own.
    honor_external : bool, default False
        ``True`` leaves a ``DataSourceSpec`` flagged ``external`` untouched -- that
        flag means "never anchor me to the project root", which is exactly the
        anchor's job and exactly what bundling makes false, so the two callers pass
        opposite values (the anchor ``True``, the bundler ``False``). A source the
        bundler does move has the flag CLEARED, since a copy inside the project is
        no longer external. It does not reach a token payload, which has no such
        flag.
    notes : list of (str, str), optional
        Appended to IN PLACE with ``(type name, message)`` for each token-bearing
        object left unchanged because it could not be written down or rebuilt. A
        ``Surface`` holding per-cell arrays cannot be serialized at all (ledger 193),
        so a stack containing one is returned UNCHANGED with a note rather than
        raised over: such a model is skipped by the save anyway, and a bundler that
        raised would cost the caller every other simulation in the project. Omitted,
        a private list is used and the notes are discarded.

    Returns
    -------
    object
        ``value`` with the paths transformed. Nothing is changed in place: a moved
        source or rebuilt token object is a NEW object, and when ``transform``
        changed nothing the SAME object comes back (for a ``DataSourceSpec`` and a
        token-bearing object; a container is always rebuilt).

    Raises
    ------
    Exception
        Whatever ``transform`` raises, unchanged. A ``from_dict`` refusing the
        repointed payload does NOT raise: it becomes a note and ``value`` comes back
        unmoved.

    Examples
    --------
    >>> import myflopy as mf
    >>> from myflopy.specs import repoint_source_paths
    >>> stack = mf.LayerStack(top=mf.Surface.raster("gis/ground.tif")).add(
    ...     "sand", thickness=10.0)
    >>> moved = repoint_source_paths(
    ...     stack, lambda path, *, role, notes: "inputs/" + path.name)
    >>> moved._top.path.as_posix(), stack._top.path.as_posix()
    ('inputs/ground.tif', 'gis/ground.tif')
    >>> repoint_source_paths(stack, lambda path, *, role, notes: None) is stack
    True
    """

    notes = [] if notes is None else notes
    if isinstance(value, DataSourceSpec):
        if honor_external and value.external:
            return value
        # `path`, and any further file a source names (`EXTRA_PATH_FIELDS`: an
        # `mf.Heads` recipe reads a head file AND the run's binary grid file), each a
        # source -- so the bundler copies, the anchor joins and `source_paths` lists all.
        changes = {}
        for name in ("path", *getattr(type(value), "EXTRA_PATH_FIELDS", ())):
            current = getattr(value, name)
            if current is None:
                continue
            stored = transform(Path(current), role="source", notes=notes)
            if stored is not None and str(stored) != str(current):
                changes[name] = stored
        if not changes:
            return value
        # `external` means "never anchor to the project root". Only the bundler gets
        # here with it set, and only by copying the file inside -- which is what makes
        # the flag false. The anchor returned above without touching it.
        return replace(value, **changes, external=False)
    if isinstance(value, ArrayProduct):
        # Its factors are sources, walked as sources: the bundler copies a factor's
        # raster, the anchor joins a relative one to the project root and honours
        # `external`, and `source_paths` lists each -- with no wiring of their own.
        factors = tuple(
            repoint_source_paths(
                factor, transform, honor_external=honor_external, notes=notes
            )
            for factor in value.factors
        )
        if all(new is old for new, old in zip(factors, value.factors, strict=True)):
            return value
        return ArrayProduct(factors=factors)
    if getattr(type(value), "_SPEC_TOKEN", None) is not None:
        try:
            payload = value.to_dict()
        except (ValueError, TypeError) as error:
            notes.append((type(value).__name__, str(error)))
            return value
        repointed = _repoint_payload(payload, transform, notes)
        if repointed == payload:
            return value
        try:
            # Untagged for the reason `_walk_surfaces` gives: this is another funnel from
            # a payload back to an object, and the tag is `_json_value`'s, not the class's.
            rebuilt = type(value).from_dict(_untag_tuples(repointed))
        except (ValueError, TypeError, KeyError) as error:
            notes.append((type(value).__name__, f"could not be rebuilt: {error}"))
            return value
        return _carry_grid(value, rebuilt)
    if isinstance(value, dict):
        return {
            key: repoint_source_paths(
                item, transform, honor_external=honor_external, notes=notes
            )
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        rebuilt = [
            repoint_source_paths(
                item, transform, honor_external=honor_external, notes=notes
            )
            for item in value
        ]
        return type(value)(rebuilt) if isinstance(value, list) else tuple(rebuilt)
    return value


def source_paths(value: Any) -> list[str]:
    """Every file path ``value`` names, in declaration order, without changing it.

    Reuses :func:`repoint_source_paths` with a transform that RECORDS and returns
    ``None`` -- "leave this path alone" -- so the walk is the same one the bundler and
    the build-time anchor take. One walker means these three can never disagree about
    what counts as a source, which they did within a day of being written apart.

    ``write_target`` paths are excluded: a GRASS cache is a file this package WRITES,
    and what it contains is decided by the inputs, which are listed.
    """

    found: list[str] = []

    def _collect(path, *, role, notes):
        """Record the path; never change it."""

        if role != "write_target":
            found.append(str(path))
        return None

    repoint_source_paths(value, _collect, notes=[])
    return list(dict.fromkeys(found))


def stamp_source_hashes(value: Any) -> Any:
    """Return a COPY of ``value`` in which each readable file it names carries a content hash.

    Called where a surface becomes part of a LAZY package -- ``mf.disv`` on a gridless
    stack, a lazy array-form ``mf.evt``. That placement is the whole point: it mirrors
    ``_lazy_bc_spec``, which hashes its GeoPackage at declaration, while leaving a bare
    ``Surface.raster(path)`` doing no I/O at all. Hashing in ``Surface`` itself would
    put a full read of a multi-hundred-megabyte DEM behind the cheapest call in the
    layer API, for a fingerprint most surfaces never need.

    Parameters
    ----------
    value : Surface, LayerStack, Spread, CellSurfaceOffset, dict, list, tuple, or any other object
        What to stamp. A token-bearing object (one with a ``_SPEC_TOKEN``: the four
        types named) is rebuilt from its own ``to_dict()`` payload with a
        ``source_shas`` entry added to every surface that names a readable file, and
        that path RESOLVED to an absolute one (the build anchors a relative path to
        the project root, the declaration read it against the working directory --
        ledger 293). A ``dict``, ``list`` or ``tuple`` is walked and rebuilt as the
        same type. Anything else is returned as it is.

    Returns
    -------
    object
        A token-bearing object comes back as a NEW object whether or not anything was
        stamped, so the lazy package holding it owns it: a later edit to the caller's
        ``LayerStack`` cannot reach an already-declared ``mf.disv`` (ledger 270,
        amended 2026-09-25). With nothing stamped the copy EQUALS the original and
        digests as it does -- a tuple inside the payload (``mf.Spread(("c0", "c1"))``,
        a stack's ``split=(0.3, 0.7)``) comes back a tuple, not the ``{"$tuple": ...}``
        tag the payload carries. Two exceptions come back
        as the SAME object: one whose ``to_dict`` refuses (a surface holding per-cell
        arrays, ``Surface.from_array``/``from_points``), which has no payload to
        rebuild from, and anything that is not token-bearing. A path that cannot be read
        is left unhashed and keeps digesting by path -- the same fallback
        :meth:`~myflopy.sources.BoundarySource.digest_payload` makes, for the same
        reason: two unreadable sources must not digest alike.

    Raises
    ------
    ValueError, TypeError or KeyError
        From the object's own ``from_dict``, if it refuses the payload its ``to_dict``
        wrote -- a round-trip defect in that class, which used to be swallowed into an
        uncopied, unstamped return.

    Examples
    --------
    >>> import myflopy as mf
    >>> from myflopy.specs import stamp_source_hashes
    >>> stack = mf.LayerStack(top=mf.Surface.flat(100.0)).add("sand", thickness=10.0)
    >>> held = stamp_source_hashes(stack)
    >>> held is stack, held.to_dict() == stack.to_dict()
    (False, True)
    >>> _ = stack.add("clay", thickness=5.0)         # the caller edits its stack...
    >>> [layer.name for layer in held._layers]        # ...and the copy is untouched
    ['sand']
    """

    return _walk_surfaces(value)


def _stamp_payload(payload: Any) -> Any:
    """Add ``source_shas`` to every Surface payload inside ``payload``."""

    if isinstance(payload, dict):
        rebuilt = {
            key: item if key == SHA_KEY else _stamp_payload(item)
            for key, item in payload.items()
        }
        shas = {}
        for key, item in payload.items():
            if isinstance(item, str) and key in SURFACE_INPUT_KEYS:
                sha = source_fingerprint(item)
                if sha is not None:
                    shas[key] = sha
                    # RESOLVED at the same moment, and for the same reason
                    # `_lazy_bc_spec` resolves its GeoPackage: a relative path here
                    # means one thing to the declaration (the cwd) and another to
                    # `_anchored_packages` (the project root), and the build silently
                    # takes the second. Only a path we just READ is rewritten, so a
                    # surface pointing at a file that does not exist yet is untouched.
                    rebuilt[key] = PurePath(item).as_posix() if PurePath(item).is_absolute() \
                        else Path(item).resolve().as_posix()
        if shas:
            rebuilt["source_shas"] = {**payload.get("source_shas", {}), **shas}
        return rebuilt
    if isinstance(payload, list):
        return [_stamp_payload(item) for item in payload]
    return payload


def _walk_surfaces(value: Any) -> Any:
    """Rebuild every token-bearing object in ``value`` with its file hashes stamped."""

    if getattr(type(value), "_SPEC_TOKEN", None) is not None:
        try:
            payload = value.to_dict()
        except (ValueError, TypeError) as error:
            # A surface holding per-cell arrays (`Surface.from_array`/`from_points`)
            # refuses `to_dict`, so there is no payload to hash or to rebuild from. It
            # is handed back as it is; a caller that must own a copy makes one itself,
            # as `mf.disv` does.
            logger.debug(
                "%s cannot be written down (%s), so its file paths are not hashed and "
                "it is returned uncopied",
                type(value).__name__,
                error,
            )
            return value
        # ALWAYS rebuilt, stamped or not: the lazy package receiving the result has to
        # own it. This returned the CALLER'S object whenever nothing was stamped -- an
        # all-inline stack -- so a later `stack.add(...)` reached an already-declared
        # `mf.disv`, a library copy of it and its cached digest, while a stack naming a
        # file was rebuilt and did not (ledger 270, amended). A `from_dict` refusing the
        # payload its own `to_dict` wrote is a defect in that class, so it now RAISES
        # instead of silently handing back the caller's unstamped object.
        #
        # `_untag_tuples` first, as `_spec_value` does: a payload can carry the
        # `{"$tuple": [...]}` tag (`Spread.to_dict` runs `_json_value` on itself), and
        # without it this second payload-to-object funnel handed back
        # `Spread(field={'$tuple': ['a', 'b']})` -- unequal to the original while
        # digesting identically.
        rebuilt = type(value).from_dict(_untag_tuples(_stamp_payload(payload)))
        return _carry_grid(value, rebuilt)
    if isinstance(value, dict):
        return {key: _walk_surfaces(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        rebuilt = [_walk_surfaces(item) for item in value]
        return type(value)(rebuilt) if isinstance(value, list) else tuple(rebuilt)
    return value


def _carry_grid(original: Any, rebuilt: Any) -> Any:
    """Re-attach a bound grid the ``to_dict``/``from_dict`` round trip dropped.

    A ``LayerStack``'s payload is its DECLARATION and deliberately excludes ``vor``
    (ledger 269), so rebuilding one hands back a GRIDLESS stack -- and the failure is
    `LayerStack.build() needs a grid. This stack was declared without one`, about a
    stack that was declared with one. Not reachable through a lazy ``mf.disv``, which
    requires a gridless stack to begin with, but these walkers are generic and shared.
    """

    if getattr(original, "vor", None) is not None and getattr(rebuilt, "vor", None) is None:
        rebuilt.vor = original.vor
    return rebuilt


#: Options the RECIPE supplies, so an eager spec's copies must be dropped when it is
#: re-lazified or the lazy builder is handed each of them twice --
#: `GeoPackageSource.ghb() got multiple values for keyword argument 'boundnames'`.
#: Everything else in `options` is a real MF6 option the caller set and is carried over.
_RECIPE_OWNED_OPTIONS = frozenset({"stress_period_data", "boundnames"})


def _mesh_parts(grid: Any, *, verb: str) -> Any:
    """A built grid's cell polygons, MF6 cell centres and vertex rings, from whatever was passed.

    A built grid, a `GridSpec` (resolved), or a `GridRef` are all things a caller
    reasonably has in hand at a move; only the first carries its mesh already. A
    `GridSpec` is meshed in a TEMPORARY directory that is gone on return, so what is
    wanted -- the polygons, the centres (``points``) and the rings (``iverts``), which a
    solved field is interpolated on (ledger 600) -- is copied out first; the recipe
    itself is what the model keeps and the build meshes again in the run's own `_grid/`
    -- `GridSpec.resolve` refuses to pick a directory nobody named (ledger 599). ``verb``
    names the caller in the refusal.
    """

    if all(hasattr(grid, name) for name in ("gdf_vorPolys", "points", "iverts")):
        return grid
    if isinstance(grid, GridSpec):
        import tempfile
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory(prefix="myflopy-grid-") as scratch:
            built = grid.resolve(workspace=scratch)
            return SimpleNamespace(
                gdf_vorPolys=built.gdf_vorPolys,
                points=built.points,
                iverts=[list(ring) for ring in built.iverts],
                crs=getattr(built, "crs", None),
                ncpl=built.ncpl,
            )
    if getattr(grid, "gdf_vorPolys", None) is not None:
        return grid                  # polygons alone: enough for an area overlay
    raise TypeError(
        f"{verb} needs a built grid or a GridSpec, not {type(grid).__name__}. A "
        "`mf.grid_ref(key)` names a library entry and is resolved at build, so pass "
        "`project.grid(key)` instead."
    )


def _polygons_of(grid: Any, *, verb: str) -> Any:
    """The cell polygons of whatever a caller passed as a grid (:func:`_mesh_parts`)."""

    return _mesh_parts(grid, verb=verb).gdf_vorPolys


def _mesh_in_hand(model: Any) -> Any:
    """The BUILT mesh a model holds at this call, or ``None``: the one rule for "the old grid".

    ``grid=<built grid>`` -- stored as :meth:`GridSpec.from_object` -- is what the build
    resolves the model onto, REPLACING ``context.grid``, so it comes first;
    ``context=ModelContext(grid=...)`` is the other spelling of the same thing. Reading
    the context alone made a per-cell array on a model declared ``grid=`` read
    ``unchanged`` across a move -- it BUILT, on the old mesh's values -- where the same
    array declared through the context read ``unverifiable`` (ledger 599). A DEFERRED
    grid (a recipe, a pickle path, a ``grid_ref``) has no mesh until the build, so the
    context's grid is what is left, as before; with neither there is nothing to measure
    a payload against (ledger 295's recorded limit).

    The callers that measure WITHOUT resolving anything -- both ``coverage_on``\\ s and
    ``Run.resolution`` -- read it here; :func:`_declared_mesh`, behind ``with_grid``,
    ``move_report`` and ``remap=True``, reads it for every grid already in hand and
    resolves only a deferred one (ledger 597), so the two cannot disagree about a built
    mesh (ledger 303).
    """

    grid = getattr(model, "grid", None)
    if isinstance(grid, GridSpec) and grid.method == "object" and grid.obj is not None:
        return grid.obj
    return getattr(getattr(model, "context", None), "grid", None)


def _relazified(package: Any) -> Any:
    """An eagerly declared `.gpkg` package, rebuilt as its lazy self.

    An eager `.gpkg` resolves at declaration and holds finished records, so a grid move
    classifies it ``stale`` -- and that is the LARGEST stale class, because declaring
    with a ``context=`` is what every pre-B2 document teaches. It needs no overlay to
    move: it recorded the recipe it was built from (``metadata['move_recipe']``), so it
    is re-declared against the new mesh EXACTLY, which is what an overlay could only
    approximate (ledger 296).

    Returns the package unchanged when there is no recipe, when the builder cannot be
    resolved, or when the package already carries one in ``inputs`` -- a lazily
    declared package rebuilds itself and must not be rebuilt twice.
    """

    metadata = getattr(package, "metadata", None) or {}
    recipe = metadata.get(MOVE_RECIPE)
    if recipe is None or getattr(package, "inputs", None):
        return package
    try:
        # Deferred because `package_api` is import layer 5 and this is layer 1. The
        # lazy BUILDER is the one thing that cannot be derived from data alone: the
        # recipe names its package, but binding the FloPy class is what keeps the
        # declaration-time option check byte-identical to the eager form (ledger 279).
        from myflopy.package_api import _lazy_factory

        builder = _lazy_factory(recipe.package)
    except (ImportError, KeyError, AttributeError) as error:
        logger.debug(
            "package %r records a move recipe but its lazy builder could not be "
            "resolved (%s); leaving it to be classified as it stands",
            getattr(package, "name", "?"), error,
        )
        return package
    return replace(
        package,
        builder=builder,
        options={k: v for k, v in (package.options or {}).items()
                 if k not in _RECIPE_OWNED_OPTIONS},
        inputs={"stress_period_data": recipe},
        # The pin and the coverage counts describe the OLD resolution. Keeping them
        # would leave a package claiming records it no longer holds -- and the pin
        # branch of `package_move` reads `GRID_PIN_OPTIONS <= inputs`, so a stale pin
        # beside a fresh recipe is exactly the collision ledger 275 voided.
        metadata={k: v for k, v in metadata.items()
                  if k not in (GRID_PIN, GRID_PIN_OPTIONS, "coverage")},
    )


def _ref_resolved(entry: PackageRef, package: PackageSpec) -> PackageSpec:
    """The library ``package`` an ``mf.ref`` names, carrying that ref's move stamp.

    The ONE rule for reading a ref through its library, shared by the build
    (:meth:`ModelSpec.resolved_packages`), the move (``with_grid``/``move_report``) and
    ``Run.resolution``, so a stamped ref cannot be refused by one and waved through by
    another. The library package itself is never stamped: it is shared with every other
    model and run that names it.
    """

    if entry.moved_from is None:
        return package
    return package.with_metadata(
        **{MOVED_FROM: entry.moved_from, MOVED_VIA_REF: entry.key}
    )


def _checked_library(library: Any, verb: str) -> Any:
    """``library=`` as given, refused unless it can answer what a ref or grid ref names."""

    if library is None:
        return None
    if not (
        callable(getattr(library, "package", None))
        and callable(getattr(library, "grid", None))
    ):
        raise TypeError(
            f"{verb}(library=...) takes the Project whose library the model's `mf.ref` "
            f"and `mf.grid_ref` entries name, not {type(library).__name__}."
        )
    return library


def _move_entries(
    model: Any, library: Any, *, verb: str, required: bool = True
) -> list[tuple[Any, Any]]:
    """Every package entry of ``model``, beside what classifying it for a move must read.

    An inline package is its own answer. An ``mf.ref`` carries no options, inputs or
    metadata, so classifying the REF answered ``unchanged`` -- "names no cells" -- for
    every one of them, and a library ``mf.ic(strt=<array>)`` moved to a mesh of the same
    size BUILT on the old mesh's values: ledger 281(a)'s defect, on the spec side, which
    has no library. So it is resolved through ``library=`` instead (user decision,
    2026-10-15). A ref an earlier move already stamped needs no library -- its verdict
    was decided by that move -- and an unstamped one with none is REFUSED naming
    ``library=`` rather than guessed.

    ``required=False`` passes refs through as themselves instead of refusing: what
    ``with_grid`` asks for when the model holds NO mesh, so there is no pin a verdict
    could stamp -- ``mf.gwf(...).with_grid(recipe)`` gives a model its FIRST grid, and
    refusing that over a ref would refuse a call that moves nothing.
    """

    unread = [
        entry.key
        for entry in model.packages
        if isinstance(entry, PackageRef) and entry.moved_from is None
    ]
    if unread and library is None and required:
        # All of them at once: naming the first sent the caller back once per ref.
        named = ", ".join(f"mf.ref({key!r})" for key in unread)
        raise ValueError(
            f"Model {model.name!r}: {verb} cannot classify {named} -- a ref names a "
            "package in a project library, and this spec does not hold the library, so "
            "the package it names would read `unchanged` whatever it holds. Pass "
            "library=project, the Project these keys are registered in."
        )
    return [
        (entry, _ref_resolved(entry, library.package(entry.key)))
        if isinstance(entry, PackageRef) and library is not None
        else (entry, entry)
        for entry in model.packages
    ]


def _grid_free_names(pairs: list[tuple[Any, Any]], grid_free: tuple[str, ...]) -> tuple[str, ...]:
    """``grid_free=`` as the names the classifier sees: a ref's KEY becomes its package's name.

    The caller names a ref by the library key it typed (``"ic_regrid"``); the classifier
    is handed the resolved package, which carries its own name (``"ic"``).
    """

    return tuple(
        dict.fromkeys(
            name
            for entry, resolved in pairs
            if _package_entry_label(entry) in grid_free
            and (name := getattr(resolved, "name", None)) is not None
        )
    )


def _declared_mesh(model: Any, library: Any, *, verb: str) -> Any:
    """The built mesh ``model``'s packages were declared on, resolved if deferred; or ``None``.

    ``grid=`` first, because it is what the build resolves the model onto, REPLACING
    ``context.grid``; then ``context.grid``. A DEFERRED grid -- a recipe, a pickle, an
    ``mf.grid_ref`` -- used to count as "no mesh", so nothing could measure a per-cell
    array or pin a hand-typed cellid against it, and both moved as ``unchanged``: the
    limit ledgers 271 and 295 recorded, measured as an inline ``mf.ic(strt=<array>)``
    building on the old mesh's values. It is resolved here instead, as the build would
    resolve it (user decision, 2026-10-15): a :meth:`GridSpec.from_object` costs nothing,
    a recipe is meshed once in a temporary directory removed on return -- a move must not
    write a mesh into the project -- anchored to ``library.root`` when one is given. A
    ``grid_ref`` needs ``library=`` and is refused naming it without one.

    ``None`` when the model holds no grid at all, or one whose route does not resolve
    (``NotImplementedError``: ``GridSpec.structured``/``from_geopackage``) -- a model
    that cannot be built either, so there was never a mesh to measure against.
    """

    grid = getattr(model, "grid", None)
    if grid is None or (isinstance(grid, GridSpec) and grid.method == "object"):
        # A mesh already IN HAND -- the one rule `coverage_on` and `Run.resolution` read
        # too; only a deferred grid needs resolving below.
        return _mesh_in_hand(model)
    if isinstance(grid, GridRef):
        if library is None:
            raise ValueError(
                f"Model {model.name!r}: {verb} measures each package against the mesh it "
                f"was declared on, and that mesh is mf.grid_ref({grid.key!r}) -- a project "
                "library entry this spec does not hold. Pass library=project so it can be "
                "resolved."
            )
        grid = _grid_entry(library.grid(grid.key))
        if not isinstance(grid, GridSpec):
            # The build's own rule (`_library_grid_spec`): a library entry must be a
            # resolvable GridSpec, not another reference.
            raise TypeError(
                f"Project grid library entry {model.grid.key!r} must be a GridSpec, "
                f"got {type(grid).__name__}"
            )
    if grid.method == "object":
        return grid.obj
    import tempfile

    with tempfile.TemporaryDirectory(prefix="myflopy-declared-grid-") as scratch:
        try:
            return grid.resolve(
                project_root=getattr(library, "root", None), workspace=scratch
            )
        except NotImplementedError as error:
            logger.debug(
                "model %r: its grid %r does not resolve (%s), so %s has no mesh to "
                "measure its packages against",
                getattr(model, "name", "?"), grid.name, error, verb,
            )
            return None


def _move_stamped(package: Any, pin: str) -> Any:
    """``package`` stamped as moved away from mesh ``pin`` -- unless a move already did.

    A second move must not overwrite the first stamp: the values still mean the mesh
    they were DECLARED on, and the refusal names that one.
    """

    if MOVED_FROM in (package.metadata or {}):
        return package
    return package.with_metadata(**{MOVED_FROM: pin})


def _remap_engine():
    """The overlay remap module, imported on first use (it pulls in geopandas).

    One deferred import serving both callers, so the import ratchet does not move.
    """

    from myflopy.modflow.mf6.grid import remap

    return remap


def _remap_choices(
    remap: bool, rule: str | None, uncovered: str | None
) -> tuple[str | None, str | None]:
    """Check ``with_grid``'s ``rule=``/``uncovered=`` AT THE CALL; return them resolved.

    Both act only under ``remap=True``, so given without it they are REFUSED, naming
    ``remap=True`` -- a ``rule=`` the move then ignored read as a rule applied. Under
    ``remap=True`` an unknown value raises here, before any work, instead of only once
    some package turned out to hold an array to carry: a typo was accepted whenever
    there was nothing to remap.
    """

    if not remap:
        given = [
            f"{label}={value!r}"
            for label, value in (("rule", rule), ("uncovered", uncovered))
            if value is not None
        ]
        if given:
            raise ValueError(
                f"{' and '.join(given)} only act{'s' if len(given) == 1 else ''} with "
                "remap=True, which carries per-cell arrays onto the new mesh by area "
                "overlay. Pass remap=True, or drop "
                + ("it." if len(given) == 1 else "them.")
            )
        return None, None
    engine = _remap_engine()
    rule = engine.DEFAULT_RULE if rule is None else rule
    uncovered = engine.DEFAULT_UNCOVERED if uncovered is None else uncovered
    engine.check_remap_choices(rule, uncovered)
    return rule, uncovered


def _remapped_package(
    package, *, old_polys, new_polys, ncpl, domain, rule, uncovered, old_grid=None,
    new_grid=None,
):
    """Carry one package's per-cell arrays onto the new mesh, or leave it alone.

    Returns ``(package, remapped)``. A package is a candidate only when
    :func:`~myflopy.grid_identity.package_move` would call it ``unverifiable`` FOR THE
    ARRAY REASON -- it carries ``CELLWISE_OPTIONS`` that no recipe supplies. Everything
    else is left exactly as it was: a pinned package is ``stale`` and re-declared, a
    hand-typed cellid list cannot be remapped at all (there are no per-cell source values
    to weight), and a lazy one rebuilds itself and is more faithful than any overlay.

    A refusal from the engine is allowed to propagate. The alternative -- catching it and
    leaving the package stamped ``unverifiable`` -- would make ``remap=True`` silently do
    nothing for an integer flag or a recharge array, which is ledger 295's failure over
    again: the caller asked for a remap and would be told the move was fine.
    """

    metadata = getattr(package, "metadata", None) or {}
    # Two routes, because only SOME packages are stamped. `mf.rch.array` records
    # `CELLWISE_OPTIONS` at declaration; `npf`/`ic`/`sto` take no `context=` and so have
    # nothing to measure against then -- `package_move` derives them from the OLD grid's
    # ncpl instead (ledger 295), and so must this. Reading only the stamp made
    # `remap=True` silently do NOTHING for the very packages the engine exists for:
    # measured, `npf` and `ic` still read `unverifiable` and not one warning fired.
    stamped = tuple(metadata.get(CELLWISE_OPTIONS) or ())
    derived = cellwise_options(package.options, ncpl) if ncpl else ()
    supplied = set(package.inputs or {})
    cellwise = tuple(
        name for name in dict.fromkeys(stamped + derived) if name not in supplied
    )
    if not cellwise:
        return package, False
    if PIN_VOIDED in metadata or metadata.get(GRID_PIN) is not None:
        # A voided package's values were never verified on the OLD mesh either, so an
        # overlay would carry an unverified guess forward and stamp it certified; a
        # pinned one is `stale` and has a better route (re-declaration). Ledger 275.
        return package, False
    if MOVED_FROM in metadata:
        # An EARLIER move stamped it, so its values mean the mesh before THAT one, not
        # the mesh in hand, and weighting them by this mesh's polygons carries the wrong
        # cells. Measured without this: the report read `remapped` while the surviving
        # stamp still refused the build -- the disagreement `move_report` exists to
        # prevent. Unreachable while a moved model held no mesh; reading `grid=` as the
        # declared mesh (2026-10-15) made it live. It stays stamped, and says so.
        return package, False
    moved_options, report = _remap_engine().remap_arrays(
        dict(package.options),
        old_polys=old_polys,
        new_polys=new_polys,
        cellwise=cellwise,
        package=package.name,
        rule=rule,
        uncovered=uncovered,
        domain=domain,
        old_grid=old_grid,
        new_grid=new_grid,
    )
    stamp = {
        "was": grid_pin_of_polygons(old_polys, ncpl),
        "filled": report["filled"] + report.get("strt_filled", 0),
        "filled_area": report["filled_area"],
        "rule": rule,
    }
    # `replace`, not `with_options`: that verb VOIDS the pin of a certified package and
    # `package_move` then reports "its records were edited after the library resolved
    # them", which is a lie about a correct remap (measured). Nothing was edited.
    remapped = replace(package, options=moved_options)
    return remapped.with_metadata(**{REMAPPED_FROM: stamp}), True


def grid_pin_of_polygons(polys, ncpl) -> str:
    """A short label for the mesh a remap took its values FROM.

    Not `grid_pin`: that needs the grid OBJECT's connectivity and a remap is handed cell
    polygons. The label only has to identify the source mesh in a report, so it is the
    cell count and a hash of the bounding box -- enough to tell two meshes apart in a
    sentence, and explicitly not a certificate.
    """

    import hashlib

    bounds = tuple(round(float(one), 6) for one in polys.total_bounds)
    digest = hashlib.sha256(repr((int(ncpl or len(polys)), bounds)).encode()).hexdigest()
    return digest[:16]


def _anchored_packages(
    packages: tuple[PackageSpec, ...],
    build_context: SpecBuildContext | None,
) -> tuple[PackageSpec, ...]:
    """Resolve every relative source path in ``packages`` against the project root.

    ``project.save(bundle_sources=True)`` stores paths RELATIVE to the project root --
    that is what makes a bundled project portable. Nothing resolved them on the package
    side, so a bundled project rebuilt from any directory but its own root opened
    ``inputs/boundaries.gpkg`` against the current working directory and died. Grids
    never had the bug: ``grid_spec_resolver._source_path`` has anchored them all along,
    and this is the same rule one tier over.

    **After ``resolved_packages``, never before.** A ``PackageRef`` is not a
    ``PackageSpec`` and carries no ``inputs``, so anchoring the DECLARED list silently
    skips every package reached through ``mf.ref`` -- the documented library route, and
    ledger 280's exact miss in a different function.

    Anchoring here rather than at load time is deliberate and measured: rewriting the
    stored spec to an absolute path MOVED the package digest, and the next plain
    ``save()`` then wrote that absolute path back, quietly un-bundling the archive. The
    library keeps the relative path; only this build's copy is absolute.
    """

    root = getattr(build_context, "project_root", None) if build_context else None
    if root is None:
        return packages

    def _anchor(path, *, role, notes):
        """Absolutize one relative path; leave an absolute one alone."""

        return None if PurePath(path).is_absolute() else Path(root) / path

    anchored = []
    for package in packages:
        if not package.inputs:
            anchored.append(package)
            continue
        inputs = repoint_source_paths(package.inputs, _anchor, honor_external=True)
        anchored.append(package if inputs == package.inputs else replace(package, inputs=inputs))
    return tuple(anchored)


def package_digest(package: PackageSpec) -> str | None:
    """Return a content fingerprint of what ``package`` BUILDS.

    Fingerprints the FIELDS -- name, builder reference, options, inputs, enabled,
    requires -- rather than the spec object. That distinction is load-bearing: a
    ``functools.partial`` builder defines no ``__eq__``, so hashing the object
    graph gives a different answer after a JSON round trip even when every field
    compares equal, which would make every reloaded package look edited.

    Measured stable across processes, ``PYTHONHASHSEED``, ``to_dict``/
    ``from_dict``, ``to_yaml``/``from_yaml`` and ``Project.save``/``load``; it
    distinguishes ``k=1.0`` from ``k=50.0`` and one cell changed by ``1e-9`` in a
    9,405-cell array.

    Cost has the same shape as forking: it is pickling PYTHON OBJECTS, not bytes.
    A scalar package is ~0.04 ms and a 72x9,405 array package ~5 ms, but the
    USG-derived 72-period, 9,090-record list EVT measures **1.3 s**. That is why
    :meth:`~myflopy.library.PackageLibrary.digests` caches until the library
    changes -- ``RunLibrary.summary`` fingerprints the whole library per call,
    and Jupyter calls it on every bare ``project.runs``.

    ``metadata`` is deliberately EXCLUDED: it carries provenance (``forked_from``,
    ``edited``) that says where a spec came from, not what it builds. Including it
    would make stamping a note look like a change to the model.

    Returns ``None`` when the package cannot be fingerprinted at all -- a builder
    that is not importable, or an option that does not pickle. Such a package is
    still perfectly buildable, so this reports rather than raises.

    What it does NOT see: the BODY of an importable builder (only its reference),
    and the CONTENTS of a file a source points at (only its path). A digest says
    the recipe is unchanged, not that the world is. For a ``Surface``, a
    ``LayerStack`` or a source it sees the serialized DECLARATION, so an attribute
    ``to_dict`` drops -- ``LayerStack._draft``, ``vor`` -- is invisible; two stacks
    bound to different grids therefore digest alike.
    """

    try:
        # EVERY field goes through `_canonical`, including the three scalar ones.
        # `options` and `inputs` always did; `name`, `enabled` and `requires` did not,
        # and `to_dict`/`from_dict` normalize them (`bool(...)`, and the JSON round
        # trip for the strings). So a `np.bool_` `enabled`, a `np.str_` `name` or a
        # `np.str_` in `requires` digested one way live and another after a
        # `Project.save()` + `load()`, and `runs.drift()` reported a package nobody
        # touched as EDITED -- permanently, since the reloaded value never goes back.
        # `np.bool_` is what any numpy predicate returns (`arr.any()`, `mask[i]`, a
        # comparison), so this is reachable without anyone meaning to. Ledger 265.
        payload = (
            _canonical(package.name),
            _callable_ref(package.builder),
            pickle.dumps(_canonical(package.options), protocol=5),
            pickle.dumps(_canonical(package.inputs), protocol=5),
            _canonical(package.enabled),
            _canonical(tuple(sorted(package.requires))),
        )
    except Exception as error:  # noqa: BLE001
        # BROAD ON PURPOSE, and it became necessary in this commit. `_canonical` now
        # calls `to_dict()` on any token-bearing or source object, which is arbitrary
        # code -- so the failure set stopped being the four bounded ones below
        # (`ValueError`, `TypeError`, `AttributeError`, `pickle.PicklingError`) and
        # became open. Measured with a token class whose `to_dict` raises `KeyError`:
        # `project.packages.digests()` and `project.runs.summary()` both DIED, where
        # before they reported a digest. That breaks this function's documented
        # contract twice over -- it "reports rather than raises", and `workspace.py`
        # requires that a bookkeeping column never make `prepare_run` fail. Jupyter
        # calls `RunLibrary.summary` on every bare `project.runs`, so the blast radius
        # of a raise is the user's project browser, not a test.
        #
        # A builder that is not importable (a closure, a lambda) or an option
        # that will not pickle (an open file handle) cannot be fingerprinted.
        # That is NOT a reason to refuse the package: it was buildable and
        # runnable before drift existed, and a bookkeeping column must never
        # gate whether a model can be prepared. Callers report it as
        # "unfingerprinted" rather than guessing at "unchanged".
        logger.debug(
            "Package '%s' cannot be fingerprinted (%s); drift will report it "
            "as unfingerprinted.",
            package.name,
            error,
        )
        return None
    return hashlib.sha256(pickle.dumps(payload, protocol=5)).hexdigest()[:16]


def _renamed_options(
    options: dict[str, Any], old: str, new: str, supplied: Iterable[str] = ()
) -> dict[str, Any]:
    """Carry a package rename into the option values that embed the old name.

    A package name is not just ``PackageSpec.name``: the helpers bake it into
    ``options`` too (``advanced._named_options`` sets ``pname=name`` and
    ``filename="{model_name}.{name}"``). A fork that renamed only the dataclass
    field left those stale, and MF6 then keeps ONE package of two under a FloPy
    "Replacing existing package" warning -- measured, two CHD specs on one model
    collapsed to a single ``flow.chd`` holding only the second.

    Only those two are rewritten, and only when they still hold the value the
    helper wrote. Anything else naming the package RAISES rather than being
    guessed at, because the conventions genuinely conflict and ``PackageSpec``
    carries no package-type field to dispatch on: uzf/lak/sfr write
    ``{model_name}_{name}_{suffix}`` (``package_registry``) while mvr writes
    ``{model_name}.{name}.bud`` (``advanced``). A ``str.replace`` across them
    would corrupt ``{model_name}_sfr_stage.sfr`` when renaming a package called
    ``sfr``.

    NOT handled, and out of reach from one spec: a sibling's ``requires=`` tuple,
    MVR's ``options["packages"]`` and its ``perioddata`` source/destination
    strings, and ``MoverConnection.package``. Renaming a package another one
    names is the caller's to fix.
    """

    renamed = dict(options)
    if renamed.get("pname") == old:
        renamed["pname"] = new
    if renamed.get("filename") == f"{{model_name}}.{old}":
        renamed["filename"] = f"{{model_name}}.{new}"
    # A containment test cannot decide these on its own: the corrected value
    # "{model_name}_sfr_north_budget.sfr" still contains "sfr" twice -- once
    # because the new name embeds the old, once because ".sfr" is the file
    # suffix. So the question asked is whether the CALLER took responsibility for
    # the key in this call, which is answerable exactly.
    handled = set(supplied)
    stale = sorted(
        key
        for key, value in renamed.items()
        if key.endswith("_filerecord")
        and key not in handled
        and isinstance(value, str)
        and old in value
    )
    if stale:
        raise ValueError(
            f"Cannot rename package '{old}' to '{new}': "
            f"{', '.join(stale)} still names '{old}', and the output-file "
            "conventions differ per package type (uzf/lak/sfr write "
            "'{model_name}_{name}_<suffix>', mvr writes '{model_name}.{name}.bud'), "
            "so rewriting them here would be a guess. Pass the corrected "
            f"{'value' if len(stale) == 1 else 'values'} explicitly, or fork "
            "without name= and rename at the source."
        )
    return renamed


def fork_package(
    package: PackageSpec,
    /,
    *,
    name: str | None = None,
    deep: bool = True,
    note: str | None = None,
    source: str | None = None,
    stamp: str = "forked_from",
    **overrides: Any,
) -> PackageSpec:
    """Return a copy of ``package`` with ``overrides`` applied, BUILT not patched.

    The obvious implementation is silently wrong. :class:`PackageSpec` is
    ``frozen=True, slots=True`` and validates in ``__post_init__`` ->
    :func:`_check_options`, but ``copy.deepcopy`` reconstructs a slots dataclass
    through ``__reduce_ex__``, so ``__post_init__`` never runs. Deep-copying the
    spec and then writing into ``copy.options`` therefore accepts a typo silently
    and surfaces it at MF6 build time. Measured:

    - ``deepcopy`` + in-place option write -> ``kk_typo`` accepted
    - ``replace`` over deep-copied option VALUES -> ``TypeError: ... does not
      accept 'kk_typo'. Did you mean 'k'?``

    So the values are copied and handed to :func:`dataclasses.replace`, which is
    the ordinary construction path -- every check guarding a hand-written
    ``mf.npf(...)`` guards a forked one. An edited discretization ``length_units`` must be MF6's
    own word, as at ``mf.disv`` (ledger 432), and an edited TDIS ``time_units`` one of
    MF6's six, as at ``mf.tdis`` (ledger 456).

    An override REPLACES a recipe in ``inputs`` that fills the same option -- an
    ``mf.Polygons``/``mf.RasterValues`` array, a product or per-layer list of them, any
    other source spec -- which is dropped and named in the stamp's
    ``replaced_inputs`` (:func:`replaced_sources`, ledger 561). A recipe
    the builder itself requires (a lazy ``.gpkg`` declaration's ``stress_period_data``)
    is not a source for an option and is still refused.

    ``deep=True`` (the default) also unaliases options the fork never touched,
    which is where the real hazard lives: :meth:`PackageSpec.with_options`
    rebuilds the options dict but not its values, so a "copy" can write through
    into the parent's arrays. Measured: cost scales with the number of PYTHON OBJECTS copied, not with cell count:
    a 9,405-cell ``npf`` forks in 0.1 ms and a 72x9,405 recharge ARRAY in 1.7 ms
    (arrays copy at memcpy speed), while the USG-derived 72-period, 9,090-record
    list EVT takes **808 ms** for its 654,480 records. ``deep=False`` is 0.1 ms
    on that same package and opts back into sharing.
    """

    _check_unit_overrides(
        package.builder,
        overrides,
        route="packages.edit" if stamp == "edited" else "packages.fork",
    )
    options = copy.deepcopy(package.options) if deep else dict(package.options)
    # Overrides merge FIRST, so a caller who supplies a corrected filerecord
    # satisfies the rename check below rather than being refused for a value they
    # just replaced.
    options.update(overrides)
    if name is not None and name != package.name:
        options = _renamed_options(options, package.name, name, overrides)
    inputs = copy.deepcopy(package.inputs) if deep else dict(package.inputs)
    # An override REPLACES the recipe that would fill the same option (ledger
    # 561): left beside it, the pair built "supplied twice".
    inputs, replaced = replaced_sources(package.builder, inputs, overrides)
    metadata = copy.deepcopy(package.metadata) if deep else dict(package.metadata)
    # An override REPLACES a resolved option wholesale (measured: 11 records -> 1), so
    # the library's pin stops being true of what is there now. Keyed on the override
    # NAMES, not on a diff of values -- so a `name=` rename, which rewrites
    # `pname`/`filename` through `_renamed_options`, does not void (ledger 275).
    metadata = void_pin(metadata, touched=set(overrides))
    provenance: dict[str, Any] = {
        "at": datetime.now(timezone.utc).isoformat(timespec="seconds")
    }
    if source is not None:
        provenance["key"] = source
    if note is not None:
        provenance["note"] = note
    if replaced:
        provenance["replaced_inputs"] = replaced
    metadata[stamp] = provenance
    return replace(
        package,
        name=package.name if name is None else name,
        options=options,
        inputs=inputs,
        metadata=metadata,
    )


def _fills_an_option(value: Any) -> bool:
    """Whether an ``inputs`` entry is a SOURCE for the option of its name.

    A source spec (``DataSourceSpec``/``LiteralSource`` -- ``mf.Polygons`` and
    ``mf.RasterValues`` among them), an ``mf.RasterValues * factor`` product, or a
    per-layer list / per-period dict holding an array recipe -- the shapes the array
    helpers declare. A lazy declaration's plain builder keywords are not.
    """

    if isinstance(value, (DataSourceSpec, LiteralSource)) or _is_array_recipe(value):
        return True
    if isinstance(value, (list, tuple)):
        return any(_is_array_recipe(item) for item in value)
    if isinstance(value, dict):
        return any(_is_array_recipe(item) for item in value.values())
    return False


def replaced_sources(
    builder, inputs: Mapping[str, Any], overrides: Mapping[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """``inputs`` without the recipes ``overrides`` replace, and their names.

    The ONE rule for an edit that names an option a recipe fills (ledger
    561): the override says what the caller wants, so the recipe goes.
    Left beside it, the pair was accepted at the edit and refused at BUILD as
    "supplied twice" (``fork_package(mf.npf(k=mf.Polygons(...)), k=array)``,
    measured) -- or, for a product or a per-layer list, refused at the edit as "part
    of this package's recipe". A recipe the builder REQUIRES as a keyword (a lazy
    ``.gpkg`` declaration's ``stress_period_data``) is the builder's input, not a
    source for an option, and stays; :func:`_check_input_collisions` refuses that
    override as before.

    Parameters
    ----------
    builder : callable
        The package's builder (``PackageSpec.builder``), whose required keyword-only
        parameters are never replaced.
    inputs : mapping
        The package's ``inputs``; not modified.
    overrides : mapping
        The option values the edit writes; only their NAMES are read.

    Returns
    -------
    tuple of (dict, list of str)
        A new ``inputs`` dict without the replaced entries, and their names, sorted.
    """

    resolves = set(_required_keyword_only(builder))
    replaced = sorted(
        key for key in set(inputs) & set(overrides)
        if key not in resolves and _fills_an_option(inputs[key])
    )
    return {key: value for key, value in inputs.items() if key not in replaced}, replaced


def substitute_values(
    package: PackageSpec,
    values: Mapping[str, Any],
    *,
    provenance: Mapping[str, Any] | None = None,
    stamp: str = "calibrated_from",
    metadata: Mapping[str, Any] | None = None,
) -> PackageSpec:
    """``package`` with option ``values`` replaced, its grid pin and lake rules KEPT.

    The write-back's substitution (ledger 560). :func:`fork_package`
    voids a grid pin an override speaks for, because a caller's records are no longer
    the ones the library resolved; here the CALLER has already checked that the
    values name the same cells (a calibrated column substituted record by record
    after each cellid was compared, or an array of the same mesh), so the pin is
    still true and is kept, and so are ``mf.lak``'s ``LAKE_RULES``/``LAKE_CONNECTIONS``
    -- voiding them made a correct value edit read as "its records were edited"
    (ledger 300's trap). A recipe in ``inputs`` filling a replaced option is dropped
    (:func:`replaced_sources`), and so is an eager ``.gpkg`` package's
    ``move_recipe``: a grid move would re-declare the package from its file and lose
    the substituted values without a word, so it reads ``stale`` instead (ledger
    296). Built through :func:`dataclasses.replace`, so every construction check runs.

    Parameters
    ----------
    package : PackageSpec
        The package to copy; not modified.
    values : mapping of str to object
        Option name to its new value, each replacing that option's whole value. The
        caller asserts each names the same cells as the value it replaces; nothing
        here can check that.
    provenance : mapping or None, optional
        Recorded under ``metadata[stamp]``, with ``replaced_inputs`` added when a
        recipe was dropped. ``None`` (default) records only ``replaced_inputs``, and
        nothing when no recipe was dropped.
    stamp : str, optional
        The metadata key the provenance goes under. Default ``"calibrated_from"``.
    metadata : mapping or None, optional
        Other metadata entries to set (``source_path``, ``source_sha``), merged over
        the package's. Metadata is not digested, so none of this moves the package's
        digest. Default ``None``.

    Returns
    -------
    PackageSpec
        The copy.

    Raises
    ------
    TypeError
        When a value's name is not an option the builder accepts, or is one a recipe
        the builder requires already supplies.
    """

    inputs, replaced = replaced_sources(package.builder, dict(package.inputs), values)
    # A move recipe re-declares the package from its source on another mesh, which
    # would rebuild the DECLARED values over the substituted ones in silence; without
    # it the pinned package reads `stale` and is re-declared by hand (ledger 296).
    stamped = {
        key: value for key, value in package.metadata.items() if key != MOVE_RECIPE
    }
    stamped.update(dict(metadata or {}))
    record = dict(provenance or {})
    if replaced:
        record["replaced_inputs"] = replaced
    if record:
        stamped[stamp] = record
    return replace(
        package,
        options={**package.options, **dict(values)},
        inputs=inputs,
        metadata=stamped,
    )


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
    """Return a compact summary for mappings, including stress-period data.

    A period key is any integer -- a ``numpy.int64`` read out of an array too, which
    ``isinstance(key, int)`` does not see (ledger 265's class) -- and never a ``bool``.
    """

    keys = list(value.keys())
    if not keys:
        return "dict(len=0)"
    if all(isinstance(key, Integral) and not isinstance(key, bool) for key in keys):
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


def _accepted_options(builder) -> set[str] | None:
    """Every keyword ``builder`` accepts, or ``None`` when that cannot be known.

    A ``PackageSpec``'s builder is one of three things: a FloPy package class
    (inspect it), a wrapper declaring its target via ``mf6_targets`` (inspect
    those, plus the wrapper's own named parameters), or a ``functools.partial``
    over one of those (unwrap and recurse).

    **A FloPy package class's own ``**kwargs`` tail is not permission.**
    ``MFPackage.__init__`` raises ``FlopyException: Extraneous kwargs`` for
    anything it does not consume, so the NAMED parameters are the accepted set.
    Treating the tail as "accepts anything" is what made the first version of
    this check pass everything through.

    A dispatching wrapper such as ``build_oc`` may construct any of four classes
    depending on a model kind not known until build time, so the accepted set is
    their UNION -- deliberately permissive. This catches a name that belongs to
    no MODFLOW 6 package of that type at all, not the wrong kind of one.
    """

    import functools

    if isinstance(builder, functools.partial):
        accepted: set[str] = set()
        for value in builder.args:
            nested = _accepted_options(value)
            if nested:
                accepted |= nested
        for value in builder.keywords.values():
            nested = _accepted_options(value)
            if nested:
                accepted |= nested
        return accepted or None

    targets = getattr(builder, "mf6_targets", None)
    if targets is None:
        if not inspect.isclass(builder):
            # A plain function with no declared target: nothing to check against.
            return None
        targets = (builder,)

    accepted = set()
    for target in targets:
        try:
            parameters = inspect.signature(target).parameters
        except (TypeError, ValueError):
            # Unintrospectable callable; without a signature there is no check.
            return None
        accepted |= {
            name
            for name, parameter in parameters.items()
            if parameter.kind not in (parameter.VAR_KEYWORD, parameter.VAR_POSITIONAL)
        }
    if not inspect.isclass(builder):
        try:
            accepted |= set(inspect.signature(builder).parameters)
        except (TypeError, ValueError):
            pass
    return accepted - {"self", "model", "simulation", "parent", "loading_package", "options"}


def _builder_label_for_error(builder) -> str:
    """A name a reader recognises: the FloPy class, not ``functools.partial``."""

    import functools

    if isinstance(builder, functools.partial):
        for value in (*builder.args, *builder.keywords.values()):
            if inspect.isclass(value):
                return value.__name__
        return _builder_label_for_error(builder.func)
    targets = getattr(builder, "mf6_targets", None)
    if targets:
        return " / ".join(target.__name__ for target in targets)
    return getattr(builder, "__name__", type(builder).__name__)


def _check_unit_overrides(builder, overrides: dict[str, Any], *, route: str) -> None:
    """Refuse an edited ``length_units`` or ``time_units`` that is not MF6's own word.

    The declaration rule of ``mf.dis``/``mf.disv``/``mf.disu`` (ledger 404) for the
    routes that edit an existing spec -- ``with_options``, ``packages.fork``,
    ``packages.edit`` -- which built a ``FT`` MF6 read as unit 0 (ledger 432); and
    ``mf.tdis``'s for a TDIS, whose ``"d"`` MF6 stops on (ledger 456).
    Keyed on the builder ACCEPTING the option, which only the DIS/DISV/DISU classes
    (``length_units``) and the TDIS (``time_units``) do, and on the OVERRIDE only:
    a length unit the spec already held is the builder's to map
    (``_length_units.mf6_length_option``), so a spec loaded from a saved project can
    still be forked.
    """

    accepted = _accepted_options(builder) or ()
    if "length_units" in overrides and "length_units" in accepted:
        _check_mf6_length_token(overrides["length_units"], helper=route)
    if "time_units" in overrides and "time_units" in accepted:
        check_tdis_time_units(overrides["time_units"], helper=route)


def _check_options(name: str, builder, options: dict[str, Any]) -> None:
    """Raise if ``options`` names something the builder cannot accept.

    Every ``mf.*`` helper ends in a ``**kwargs`` tail so that any FloPy option
    stays reachable without myflopy having to mirror hundreds of signatures.
    The cost, unchecked, is that a MISSPELLED or MISPLACED option is swallowed
    silently and surfaces much later from inside FloPy -- ``model_nam_file=``
    passed to ``mf.ims`` produced ``FlopyException: Extraneous kwargs`` at
    ``prepare_run``, five frames down and twenty lines from the mistake.

    So the tail stays open, and the names are checked here instead: at the call,
    against the target's real signature, with the closest valid name suggested.
    """

    accepted = _accepted_options(builder)
    if accepted is None:
        return
    unknown = sorted(set(options) - accepted)
    if not unknown:
        return

    import difflib

    label = _builder_label_for_error(builder)
    problems = []
    for key in unknown:
        close = difflib.get_close_matches(key, sorted(accepted), n=1, cutoff=0.7)
        problems.append(f"{key!r}" + (f" (did you mean {close[0]!r}?)" if close else ""))
    raise TypeError(
        f"{name}: {label} does not accept " + ", ".join(problems) + ". "
        "Every myflopy package helper forwards unknown keywords to its MODFLOW 6 "
        "package, so a name that package does not have is almost always a typo or "
        "an argument meant for a different one."
    )


def _is_array_recipe(value: Any) -> bool:
    """Whether ``value`` is ONE array recipe: ``mf.RasterValues``, ``mf.Polygons``, a product,
    or ``mf.Heads``."""

    return isinstance(value, (RasterValues, Polygons, ArrayProduct, Heads))


def _refuse_recipe_options(name: str, options: dict[str, Any]) -> None:
    """Refuse an array recipe sitting in a package's OPTIONS (ledger 542).

    A recipe (``mf.RasterValues``, ``mf.Polygons``, a product of them) is resolved at
    BUILD only from ``inputs``, and only by the builder the ``mf.npf``/``mf.ic``/
    ``mf.sto``/``mf.rch.array``/``mf.evt.array`` helpers wrap around FloPy's. In
    ``options`` it reaches FloPy unresolved, which dies five frames down with an
    ``MFDataException`` naming neither the recipe nor the edit. Every edit route --
    ``packages.fork``, ``packages.edit``, ``with_options``, a hand-built spec -- writes
    ``options``, so they all put a recipe there; measured: ``packages.fork("npf",
    "npf2", k=mf.RasterValues(p))`` constructed fine and failed at build.

    Shallow on purpose: the value itself, the items of a list/tuple, and the values of a
    dict -- the three shapes a recipe is ever declared in. Records are never walked: a
    654,480-record EVT would cost a scan on every construction.
    """

    for key, value in options.items():
        if isinstance(value, (list, tuple)):
            items = value
        elif isinstance(value, dict):
            items = list(value.values())
        else:
            items = (value,)
        if isinstance(value, Heads) or any(isinstance(item, Heads) for item in items):
            raise TypeError(
                f"{name}: option {key!r} holds mf.Heads(...), which is resolved at BUILD from "
                "a package's recipe inputs -- never from its options -- and fills one slot: "
                "a GWF model's starting heads, mf.ic(strt=mf.Heads(...)), declared through "
                "mf.ic and given alone."
            )
        if _is_array_recipe(value) or any(_is_array_recipe(item) for item in items):
            raise TypeError(
                f"{name}: option {key!r} holds an array recipe (mf.RasterValues, "
                "mf.Polygons or a product of them), which is resolved at BUILD from a "
                "package's recipe inputs -- never from its options, where FloPy would be "
                "handed it unresolved. Declare the package again through its helper "
                f"(mf.npf/mf.ic/mf.sto({key}=...), mf.rch.array/mf.evt.array) and "
                "register that, e.g. project.add_package(key, spec, overwrite=True); "
                "packages.fork, packages.edit and with_options cannot route a recipe."
            )


def _check_input_collisions(
    name: str, builder, inputs: dict[str, Any], options: dict[str, Any]
) -> None:
    """Raise if one keyword is set twice -- once as a recipe input, once as an option.

    ``PackageSpec.build`` splats both into the builder
    (``self.builder(target, **self.inputs, **self.options)``), so a name in both is
    always fatal. Unchecked it surfaces at BUILD as
    ``_build_lazy_disv() got multiple values for keyword argument 'layer_stack'`` --
    which names the private builder rather than the call the user wrote, and arrives
    a whole `prepare_run` away from the mistake.

    It cannot be caught by :func:`_check_options`, whose question is whether a name
    is ACCEPTABLE. ``layer_stack`` is perfectly acceptable to the lazy DISV builder;
    the defect is that the recipe already supplied it.

    **A SOURCE spec is exempt, UNLESS the builder itself requires that keyword.**
    ``inputs`` carries two different things. A lazy recipe's entries are builder keyword
    arguments, and a duplicate there is a mistake. A ``DataSourceSpec``/``LiteralSource``
    entry is an UNRESOLVED value for the option of the same name -- ``inputs["k"]`` says
    where ``options["k"]`` comes from -- so there the two names matching is the point,
    and a future resolver will consume one into the other.

    The refinement is what keeps the exemption honest once a builder RESOLVES the
    recipe itself. ``_build_lazy_list_bc`` takes ``stress_period_data`` as a required
    keyword-only parameter, so both edit routes -- ``with_options`` and
    ``fork_package``/``packages.edit`` -- constructed fine and died at BUILD with
    ``functools.partial(...) got multiple values for keyword argument
    'stress_period_data'``, naming a partial object. That is precisely the failure this
    function exists to prevent, let through by the exemption it grants. A builder that
    requires the keyword will be handed BOTH, so the collision is real however the value
    is spelled.
    """

    resolves = set(_required_keyword_only(builder))
    collisions = sorted(
        key
        for key in set(inputs) & set(options)
        if key in resolves
        or not isinstance(inputs[key], (DataSourceSpec, LiteralSource))
    )
    if not collisions:
        return
    raise TypeError(
        f"{name}: " + ", ".join(repr(key) for key in collisions) + " "
        + ("is" if len(collisions) == 1 else "are")
        + " part of this package's recipe and cannot also be passed as an option -- "
        "the recipe already supplies "
        + ("it" if len(collisions) == 1 else "them")
        + ". Change the declaration that built the recipe instead."
    )


@functools.cache
def _required_keyword_only(builder) -> tuple[str, ...]:
    """Keyword-only parameters of ``builder`` that have no default.

    Cached because it runs on EVERY ``PackageSpec`` construction and
    ``inspect.signature`` is not cheap -- uncached it measured `mf.npf` at 44.2 ->
    82.5 us/call, for a check that can fire on one builder in the whole library.
    """

    try:
        parameters = inspect.signature(builder).parameters
    except (TypeError, ValueError):  # C-implemented or unintrospectable
        return ()
    return tuple(
        name
        for name, parameter in parameters.items()
        if parameter.kind is parameter.KEYWORD_ONLY
        and parameter.default is parameter.empty
    )


@functools.cache
def _wants_spec(builder) -> bool:
    """Whether ``builder`` declares a keyword-only ``_spec`` parameter.

    A builder is called as ``builder(target, **inputs, **options)`` and so learns
    nothing about the spec it came from -- not its name, not its metadata. That is fine
    for a FloPy class, which is told everything through options, and wrong for a builder
    that RESOLVES a recipe: the package it produces has to be named, and the only honest
    source of that name is the ``PackageSpec``.

    Declaring the name in the recipe instead was measured worse than useless.
    ``fork_package`` renames through ``_renamed_options``, which rewrites ``options`` and
    never touches ``inputs``, so a frozen name survives a rename unchanged -- and two
    renamed forks then collapse into ONE MF6 package under a FloPy "Replacing existing
    package" warning, which is the defect ``_renamed_options`` exists to prevent.

    Opt-in by signature rather than by a flag, so a builder that does not ask is called
    exactly as before. Cached because it runs on every build.
    """

    try:
        parameters = inspect.signature(builder).parameters
    except (TypeError, ValueError):  # C-implemented or unintrospectable
        return False
    parameter = parameters.get("_spec")
    return parameter is not None and parameter.kind is parameter.KEYWORD_ONLY


def _check_recipe_complete(name: str, builder, inputs: dict, options: dict) -> None:
    """Raise when an edit has left a recipe that can no longer build.

    Clearing ``inputs`` on a lazy package removes the very thing that resolves it.
    Measured before this, the result was worse than an error: ``package_move`` called
    it ``"unchanged"`` -- the verdict that means MOVES UNTOUCHED -- so a spec that
    could not build on ANY mesh sailed through a grid move, and the failure arrived
    later as ``_build_lazy_disv() missing 2 required keyword-only arguments``.

    Refused at construction rather than given a fifth verdict, because a label cannot
    refuse anything and the state has no legitimate producer. That follows the house
    pattern -- ledger 208 refuses a duplicate registration outright instead of
    reporting it.

    Keyword-only is the right narrowness: a positional parameter is filled by the
    build target or bound by a ``functools.partial``. Measured, it refuses nothing that
    exists today -- 108 ``flopy.mf6.Modflow*`` classes have zero required keyword-only
    parameters, and ``_build_lazy_disv`` is the only builder in the library with any.
    """

    required = _required_keyword_only(builder)
    if not required:
        return
    # `_spec` is filled by `build` itself, not by the recipe, so it is never missing in
    # the sense this check means. Without the exclusion a builder that asks for its own
    # identity is refused at construction -- ledger 275's guard firing on ledger 278's
    # mechanism.
    missing = sorted(set(required) - set(inputs) - set(options) - {"_spec"})
    if not missing:
        return
    # `builder.__name__`, NOT `_builder_label_for_error`: since ledger 272 the lazy
    # builder declares `mf6_targets`, so that helper answers with the four FloPy DISV
    # classes -- none of which requires these names. Naming the wrong thing in an
    # error is how a reader loses an afternoon.
    label = getattr(builder, "__name__", type(builder).__name__)
    raise TypeError(
        f"{name}: this package's recipe no longer supplies "
        + ", ".join(repr(key) for key in missing)
        + f", which {label} requires. Clearing `inputs` removes the recipe that "
        "resolves the package, and what is left cannot build on ANY mesh. Re-declare "
        "it rather than editing `inputs`."
    )


@dataclass(frozen=True, slots=True)
class PackageSpec:
    """Reusable instructions for building one package on a model or simulation.

    Builders receive the target as their first positional argument followed by
    the values in ``options``. A builder can be a FloPy package class or a
    small project-specific function.

    The ``mf.*`` helpers (``mf.npf``, ``mf.chd.gpkg``, ...) construct these for
    you; build one by hand only for a builder no helper covers.

    Parameters
    ----------
    name : str
        The package's identity inside its model (or simulation): the key
        :meth:`ModelSpec.package`, ``with_packages(slot=...)`` and
        ``without_packages`` match on, the key of ``BuiltModel.packages``, and
        the name ``requires`` entries refer to. It must be unique among the
        packages of one model. It is NOT passed to the builder, so it is not the
        FloPy package name: the helpers put that in ``options`` (``pname``), and
        a bare spec's package is named by FloPy (``chd_0``) unless ``options``
        says otherwise. A builder that declares a keyword-only ``_spec``
        parameter receives the whole spec and can read it.
    builder : callable
        Called at build time as ``builder(target, **inputs, **options)``, where
        ``target`` is the FloPy model (or ``MFSimulation`` for a
        simulation-level package); a builder declaring a keyword-only ``_spec``
        parameter is also passed ``_spec=self``. Usually a FloPy package class
        (``flopy.mf6.ModflowGwfnpf``), a ``functools.partial`` over one, or a
        module-level function. :meth:`to_dict` stores it as an import path, so
        it must be importable by that path to round-trip: a lambda (at module
        level too) or a function defined inside another function builds but
        raises ``ValueError`` in :meth:`to_dict`, so ``Project.save`` reports the
        simulation as skipped instead of writing one it could not load.
    options : dict, optional
        Keyword arguments for the builder. Checked at construction against what
        the builder accepts -- the named parameters of a FloPy class, of the
        classes a wrapper lists in ``mf6_targets`` (plus the wrapper's own), or
        of the class a ``functools.partial`` binds -- and an unknown name raises
        ``TypeError``, suggesting the closest valid name when one is close. A
        plain function with no ``mf6_targets``, or a partial binding no class,
        is not checked. An array recipe (``mf.RasterValues``, ``mf.Polygons`` or
        a product of them) as a value, a list item or a dict value raises
        ``TypeError``: it is resolved only from ``inputs``, by the builder the
        array helpers wrap (ledger 542). Default: empty.
    enabled : bool, optional
        ``False`` keeps the package declared but skipped: :meth:`build` returns
        ``None``, and an enabled sibling that ``requires`` it is refused with
        ``ValueError`` as requiring a missing package, wherever ``requires`` is
        checked. Default ``True``.
    requires : sequence of str, optional
        ``name`` values of packages in the same model that must be present,
        enabled, and declared EARLIER in the package list; each violation raises
        ``ValueError`` when the packages are validated -- at ``ModelSpec`` /
        ``SimulationSpec`` construction when every entry is a concrete spec,
        else when the refs are resolved at build. Checked only for an enabled
        package. Stored as a tuple. Default: empty.
    metadata : dict, optional
        Provenance and bookkeeping written by the helpers and the grid machinery
        -- the grid pin (``"grid_pin"``) and the options it certifies, move and
        remap stamps, ``forked_from``/``edited`` notes, source hashes. It is not
        passed to the builder as arguments (a builder declaring ``_spec`` can
        read it off the spec), and :func:`package_digest` excludes it, so
        stamping a note does not read as a change to the model. Default: empty.
    inputs : dict, optional
        The package's unresolved RECIPE, splatted into the builder beside
        ``options``. Either keyword arguments for a recipe-resolving builder
        (the lazy ``mf.disv``/``.gpkg``/``rch``/``evt`` forms), or a
        source spec (``DataSourceSpec``/``LiteralSource``, such as
        ``mf.Polygons``) standing for the option of the same name. A key in
        both ``inputs`` and ``options`` raises ``TypeError`` -- at construction,
        unless the ``inputs`` value is a source spec the builder does not itself
        require, and then at :meth:`build`. Leaving out a keyword-only argument
        the builder requires also raises ``TypeError`` at construction.
        Default: empty.
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
        _check_options(self.name, self.builder, self.options)
        _refuse_recipe_options(self.name, self.options)
        _check_input_collisions(self.name, self.builder, self.inputs, self.options)
        _check_recipe_complete(self.name, self.builder, self.inputs, self.options)

        object.__setattr__(self, "requires", tuple(self.requires))

    def build(self, target: Any) -> Any | None:
        """Build the package on ``target``, or return ``None`` when disabled.

        Parameters
        ----------
        target : flopy.mf6.MFModel or flopy.mf6.MFSimulation
            What the package is built on, passed to the builder as its first
            positional argument: the FloPy model for a model package, the
            simulation for ``tdis``/``ims``/``ems``.

        Returns
        -------
        object or None
            Whatever the builder returns (normally the FloPy package), or
            ``None`` when ``enabled`` is ``False``.

        Raises
        ------
        TypeError
            A key present in both ``inputs`` and ``options``.
        """

        if not self.enabled:
            return None
        doubled = sorted(set(self.inputs) & set(self.options))
        if doubled:
            # `_check_input_collisions` lets a SOURCE in `inputs` shadow the option of
            # the same name through, on the premise that a resolver would consume one
            # into the other. The resolver arrived (ledger 299) and it REPLACES the
            # option rather than merging with it, so this state is unbuildable -- and
            # what it produced was FloPy's own `got multiple values for keyword
            # argument`, naming neither the package nor the recipe. Measured: the
            # SUPPORTED shape got the worse error and the unsupported one got the good
            # message, which is the guard exactly inverted.
            raise TypeError(
                f"{self.name}: "
                + ", ".join(repr(key) for key in doubled)
                + (" is" if len(doubled) == 1 else " are")
                + " supplied twice -- once by this package's recipe and once as an "
                "option. A recipe REPLACES the option it fills; it cannot sit beside "
                "it. Drop the option, or drop the recipe."
            )
        if _wants_spec(self.builder):
            # A recipe-resolving builder needs the identity of the spec it came from --
            # its NAME above all, since the package it produces has to be named and the
            # recipe deliberately does not carry one (ledger 277).
            return self.builder(target, _spec=self, **self.inputs, **self.options)
        return self.builder(target, **self.inputs, **self.options)

    def with_options(self, **overrides: Any) -> PackageSpec:
        """Return a copy with selected options added or replaced.

        Voids a grid pin the overrides speak for, like every other edit route
        (ledger 275). Measured safe for the five internal callers: every one of them
        operates on an unpinned spec (`metadata == {}` for `to_disv`, `mf.disv(result)`
        and `mf.ims()`), so nothing in the library voids its own pin here.

        Parameters
        ----------
        **overrides
            Option names and values merged over ``options`` key by key (an
            override replaces the whole value of its key; others are kept). The
            copy is re-validated like a new spec, so a name the builder does not
            accept raises ``TypeError``, and so does one a recipe in ``inputs``
            already supplies as a keyword the builder requires (a lazy ``.gpkg``
            declaration's ``stress_period_data``). A recipe that FILLS the option
            -- an ``mf.Polygons``/``mf.RasterValues`` array, a product or per-layer
            list of them -- is replaced by the override and dropped from
            ``inputs`` (ledger 561).
            Overriding an option the grid pin certifies
            -- any option, for a pin that names none -- replaces the pin with a
            ``grid_pin_voided`` record naming the edit. On a DIS/DISV/DISU, a
            ``length_units`` must be MF6's own word -- ``"feet"``, ``"meters"``,
            ``"centimeters"`` or ``"unknown"``, any case -- as at ``mf.disv``:
            another (``"ft"``, which MF6 reads as undefined) raises ``ValueError``
            naming the word to write (ledger 432). On a TDIS, a ``time_units``
            must be one of MF6's six words, as at ``mf.tdis`` -- MF6 stops the run on
            any other (ledger 456). An array recipe (``mf.RasterValues``,
            ``mf.Polygons``, a product) raises ``TypeError``: it resolves only from
            ``inputs``, so declare it through the package's helper (ledger
            542).
        """

        _check_unit_overrides(
            self.builder, overrides, route=f"{self.name}.with_options"
        )
        inputs, _replaced = replaced_sources(self.builder, self.inputs, overrides)
        return replace(
            self,
            options={**self.options, **overrides},
            inputs=inputs,
            metadata=void_pin(self.metadata, touched=set(overrides)),
        )

    def with_metadata(self, **updates: Any) -> PackageSpec:
        """Return a copy with selected provenance metadata added or replaced.

        Parameters
        ----------
        **updates
            Keys and values merged over ``metadata`` (an update replaces its
            key's value; other keys are kept). The keys are not checked and no
            pin is voided. Metadata is not passed to the builder as arguments.
        """

        return replace(self, metadata={**self.metadata, **updates})

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready representation of this package recipe."""

        return {
            "kind": "PackageSpec",
            "name": self.name,
            "builder": _callable_ref(self.builder),
            "inputs": _json_value(self.inputs),
            "options": _json_value(self.options),
            # Through the converter like every other field: `np.bool_` is NOT a
            # `bool` subclass, so a hand-built spec carrying one wrote a value
            # `json.dumps` refuses -- and numpy 2 calls its own type `bool`, so
            # the error read "Object of type bool is not JSON serializable".
            "enabled": _json_value(self.enabled),
            "requires": list(self.requires),
            "metadata": _json_value(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PackageSpec:
        """Recreate a package specification from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A payload as :meth:`to_dict` writes it. ``"name"`` and ``"builder"``
            (an import-path reference) are required -- a missing one raises
            ``KeyError``, a malformed reference ``ValueError``, and one that no
            longer imports ``ImportError``/``AttributeError``.
            ``"inputs"``, ``"options"``, ``"enabled"``, ``"requires"`` and
            ``"metadata"`` default to empty/``True``. The rebuilt spec is
            validated like a new one.
        """

        return cls(
            name=data["name"],
            builder=_resolve_callable(data["builder"]),
            inputs=_spec_options(data.get("inputs", {})),
            options=_spec_options(data.get("options", {})),
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


def _flatten_packages(items: Iterable[Any]) -> list[Any]:
    """Flatten one level of nesting, so a verb takes one thing or a list of them.

    ``with_packages(npf)``, ``with_packages(npf, ic)`` and
    ``with_packages([npf, ic])`` all mean the same thing. Strings are NOT
    iterated apart -- a package name is one item, not a list of characters.
    """

    flat: list[Any] = []
    for item in items:
        if isinstance(item, str) or not isinstance(item, Iterable):
            flat.append(item)
        else:
            flat.extend(item)
    return flat



def _package_scope(builder) -> str | None:
    """``"model"``, ``"simulation"``, or ``None`` when it cannot be told.

    Read from the BUILDER's first constructor parameter, never from the package
    name. Measured across every ``flopy.mf6`` package class, the first parameter
    partitions them cleanly -- ``model`` (78), ``simulation`` (15),
    ``parent_model_or_package`` (8, genuinely either: MVR/MVT/GNC),
    ``parent_package`` (7, sub-packages) -- while a name-prefix heuristic is
    WRONG on 17 of them: all six exchange classes (``gwfgwt``, ``gwfgwf``, ...)
    are simulation-scoped but carry a model-kind prefix, and eleven ``utl*``
    classes (``utllaktab``, ``utlsfrtab``, ``utlats``, ...) are model- or
    sub-package-scoped but read as simulation.

    Unwraps ``functools.partial`` for the same reason :func:`_accepted_options`
    does, and it matters more here: the seven list BCs all build through
    ``partial(_build_named, ModflowGwfchd)`` and measure as unknown without it --
    so a check meant to keep ``chd`` off the simulation would have waved through
    exactly the packages it was written for while catching ``npf``.

    ``None`` for anything undecidable, including a ``PackageRef`` (no builder at
    all): a ref cannot be classified before it resolves, and refusing one would
    break the library indirection.
    """

    import functools

    if builder is None:
        return None
    if isinstance(builder, functools.partial):
        for value in (*builder.args, *builder.keywords.values()):
            scope = _package_scope(value)
            if scope is not None:
                return scope
        return None

    targets = getattr(builder, "mf6_targets", None)
    if targets is None:
        if not inspect.isclass(builder):
            return None
        targets = (builder,)

    scopes = set()
    for target in targets:
        try:
            parameters = list(inspect.signature(target.__init__).parameters)
        except (ValueError, TypeError):  # C-implemented or unintrospectable
            continue
        if len(parameters) > 1:
            scopes.add(parameters[1])
    if not scopes:
        return None
    if scopes == {"model"}:
        return "model"
    if "simulation" in scopes:
        return "simulation"
    return None


def _refuse_model_packages(
    packages: Iterable[PackageSpec], simulation_name: str
) -> None:
    """Refuse a MODEL package declared at the simulation tier.

    A simulation carries timing and solvers -- ``tdis``, ``ims``, ``ems`` -- and a
    ``chd`` handed to it was accepted at declaration and died much later inside
    FloPy as ``Invalid package type "chd". Unable to register package.``, which
    names the type but not the mistake. This raises at the call site and names
    the fix.

    Only entries whose scope is decidably ``"model"`` are refused; anything
    undecidable (a ref, MVR/GNC, a plain callable) passes, because a check that
    guesses is worse here than one that abstains.
    """

    stray = [
        package for package in packages if _package_scope(package.builder) == "model"
    ]
    if not stray:
        return
    names = ", ".join(package.name for package in stray)
    raise ValueError(
        f"Simulation '{simulation_name}' was given model package(s): {names}. "
        "A simulation carries timing and solvers (mf.tdis, mf.ims, mf.ems); "
        "model packages belong on the model -- pass them to mf.gwf(packages=[...]) "
        "or use model.with_packages(...)."
    )


def _package_entry_indices(
    packages: Sequence[PackageEntry], target: str
) -> list[int]:
    """Every entry ``target`` names: by FULL LABEL first, then by slot hint.

    The label pass is what makes ``<type>/<variant>`` keys addressable, and it
    can match at most one entry because names are unique within a model. The slot
    pass is a HINT -- only the key's first segment for a ref -- so it can match
    several, and that ambiguity belongs to the caller rather than to list order.

    ONE resolver, because three had drifted apart. On a model carrying
    ``mf.ref("wel/north")`` and ``mf.ref("wel/south")``, ``package("wel")``
    returned the north one, ``with_package`` replaced the north one, and
    ``without_package("wel")`` deleted BOTH -- three answers to one question, and
    nobody chose any of them. The label-first pass was added after a swap deleted
    the wrong well package and MF6 reported normal termination on the wrong
    model; that fix reached this lookup and neither of the other two.
    """

    labelled = [
        index
        for index, existing in enumerate(packages)
        if _package_entry_label(existing) == target
    ]
    if labelled:
        return labelled
    return [
        index
        for index, existing in enumerate(packages)
        if _package_slot(existing) == target
    ]


def _package_entry_index(
    packages: Sequence[PackageEntry], target: str, *, where: str | None = None
) -> int | None:
    """The single entry ``target`` names, or ``None``; raises when it names several.

    Resolving an ambiguous name by taking the first match is what let a swap land
    on the wrong package. ``where`` names the owner in the error so the caller can
    see which model it is being asked about.
    """

    matches = _package_entry_indices(packages, target)
    if len(matches) > 1:
        labels = ", ".join(_package_entry_label(packages[i]) for i in matches)
        owner = f"{where} has" if where else "There are"
        raise KeyError(
            f"{owner} {len(matches)} packages matching '{target}': {labels}. "
            "Name one of them exactly -- matching the slot alone would act on "
            "whichever happened to be declared first."
        )
    return matches[0] if matches else None


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
    eagerly with ``.resolve(workspace=...)``).

    Use it two ways:

    - **Deferred** -- ``mf.gwf(...).with_grid(GridSpec.voronoi(...))``; the
      :class:`~myflopy.workspace.Project` builds the grid into
      ``run.workspace/_grid/<model>`` and populates the built model's
      ``myflopy_context.grid``. Composes with disv + simple BCs.
    - **Eager** -- ``vor = GridSpec.voronoi(...).resolve(workspace=ws)``; then put
      ``vor`` in a ``ModelContext`` so the GIS package helpers
      (``mf.uzf``/``mf.sfr``/``.gpkg``) can resolve cells immediately.

    Build one with a constructor -- :meth:`voronoi`, :meth:`python`,
    :meth:`from_object`, :meth:`from_pickle` -- rather than directly; each sets
    ``method`` and the fields its route reads. Fields another route reads are
    stored but ignored.

    Parameters
    ----------
    name : str
        The grid's name: the label a library, ``with_grid`` lineage and error
        messages use, the default ``name`` of the built ``VoronoiGridPlus``, and
        the directory a :meth:`python` builder's default workspace is named after
        (``<project_root>/_grid/<name>``).
    grid_type : str
        A label for the discretization the grid feeds (the constructors default
        to ``"disv"``, and :meth:`voronoi` always writes it). Shown in the grid
        library's summary and kept by
        ``project.save()``; it does not choose how the grid is resolved.
    method : str
        The resolution route: ``"voronoi"``, ``"python"``, ``"object"`` or
        ``"pickle"``. Any other value raises ``NotImplementedError`` at
        :meth:`resolve`.
    script : Path or str, optional
        ``"python"`` route: the builder script, anchored to ``project_root`` at
        resolve when relative. Default ``None``.
    function : str, optional
        ``"python"`` route: the name of the callable in ``script`` to call.
        Default ``None``.
    source : DataSourceSpec, optional
        Serialized and restored by :meth:`to_dict`/:meth:`from_dict`, but set
        by no constructor and read by no resolution route. Default ``None``.
    boundary : DataSourceSpec, optional
        ``"voronoi"`` route: the model-domain polygon source; resolving a
        Voronoi spec without one raises ``ValueError``. Default ``None``.
    refinement : DataSourceSpec, optional
        ``"voronoi"`` route: one polygon source of refinement regions. See
        :meth:`voronoi`. Default ``None``.
    breaklines : sequence of DataSourceSpec, optional
        ``"voronoi"`` route: line sources buffered into refinement corridors.
        Stored as a tuple. Default: empty.
    points : sequence of DataSourceSpec, optional
        ``"voronoi"`` route: point sources pinned as mesh vertices. Stored as a
        tuple. Default: empty.
    inputs : sequence, optional
        ``"python"`` route: files the builder needs. Each ``DataSourceSpec``
        (unless ``external=True``) or path is checked to exist at resolve,
        anchored to ``project_root``, and a missing one raises
        ``FileNotFoundError``. Stored as a tuple. Default: empty.
    crs : str or int, optional
        ``"voronoi"`` route: the CRS every source is reprojected to and the
        grid is built in. When ``None`` no source is reprojected, and the grid
        takes ``voronoi_options["crs"]``, else ``options["crs"]``, else
        ``"EPSG:2927"``; two of the three naming different CRSs raise
        ``TypeError`` at resolve (ledger 491). The ``"python"`` route only
        stores it, for the builder to read as ``spec.crs``. Default ``None``.
    engine : str, optional
        ``"voronoi"`` route: the meshing engine. Only
        ``"triangle_voronoi_plus"`` (:meth:`voronoi`'s default) resolves;
        anything else, ``None`` included, raises ``NotImplementedError`` at
        resolve. Default ``None``.
    options : dict, optional
        ``"voronoi"`` route: the named meshing options read by key -- the
        ``boundary_*``/``refinement_*``/``breakline_*`` settings of
        :meth:`voronoi`, ``crs``, ``model_ws``/``workspace`` -- and a fallback
        for the stage keys below when the stage dict does not set them. On this
        route a key it never reads, or two spellings of one option (``max_area``
        beside ``boundary_max_area``), raises ``TypeError`` naming it -- at
        :meth:`voronoi` and again at :meth:`resolve`, so a spec set directly or
        loaded from a saved project is checked too (ledger 417). The
        ``"python"`` route stores it for the builder to read as ``spec.options``,
        unchecked. The removed ``protect_sources``/``protected_labels`` raise
        ``TypeError`` naming the key (ledger 445). Default: empty.
    triangle_options : dict, optional
        ``"voronoi"`` route: keyword arguments for ``TriangleGrid(...)``
        (``angle``, ``maximum_area`` -- Triangle's global cap on every triangle,
        not the background size ``options["boundary_max_area"]``/``max_area`` --
        ``nodes``, ``exe_name``, ``additional_args``), plus
        ``region_point_tolerance`` and
        ``model_ws``/``workspace`` (Triangle's working directory), which are
        taken out first. An unknown key raises ``TypeError`` (from FloPy's
        ``Triangle``) at resolve. Default: empty.
    mesh_options : dict, optional
        ``"voronoi"`` route: options for ``TriangleGrid.build_mesh`` --
        ``profile`` (default ``"balanced"``), ``optimize``, ``verbose`` and the
        other keys :meth:`voronoi` routes here. Any other key raises
        ``TypeError`` naming it and the accepted keys, at :meth:`voronoi` and at
        :meth:`resolve`; a key since removed (``max_optimization_points``, ledger
        410; ``protect_sources``/``protected_labels``, ledger 445) says so.
        Default: empty.
    voronoi_options : dict, optional
        ``"voronoi"`` route: keyword arguments for ``VoronoiGridPlus(tri,
        ...)`` -- ``idomain_path``, ``rasters``, ``qhull_options``,
        ``tessellation`` (an ``"idomain"`` key raises ``TypeError`` at resolve:
        a grid takes no active domain, ledger 388). ``crs`` and ``name`` default to the spec's, and
        ``idomain_path``/``rasters`` are anchored to ``project_root``. Default:
        empty.
    metadata : dict, optional
        Free-form provenance, serialized with the spec; a :meth:`python`
        builder can read it as ``spec.metadata``. No resolution route reads it.
        Default: empty.
    obj : object, optional
        ``"object"`` route: the already-built grid :meth:`resolve` returns
        unchanged. Not serializable -- :meth:`to_dict` raises for an
        ``"object"`` spec; ``project.save()`` pickles ``obj`` for a
        grid-library entry and skips a simulation holding one inline. Default
        ``None``.
    persist : str, optional
        Stored and read by nothing: ``project.save()`` pickles an ``"object"``
        grid whatever it says, and :meth:`to_dict` does not write it. Default
        ``"pickle"``.
    pickle_path : str, optional
        ``"pickle"`` route: the pickled grid :meth:`resolve` loads, anchored to
        ``project_root`` when relative. Default ``None``.

    Examples
    --------
    >>> grid = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec("in.gpkg", layer="boundary"),
    ...                            refinement=mf.GeoPackageSourceSpec("in.gpkg", layer="refine"),
    ...                            boundary_max_area=200.0)
    >>> vor = grid.resolve(workspace="runs/_grid")   # eager: a concrete VoronoiGridPlus
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
    crs: str | int | None = None
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

        Parameters
        ----------
        obj : object
            The built grid. :meth:`resolve` returns it unchanged (no copy), and
            ``resolve`` raises ``ValueError`` when it is ``None``. Not checked
            here. Not serializable: :meth:`to_dict` raises; ``project.save()``
            pickles it for a grid-library entry and skips a simulation holding
            one inline.
        name : str, optional
            The spec's name (library and lineage label). Default ``"grid"``.
        grid_type : str, optional
            A discretization label kept on the spec; it does not change what
            :meth:`resolve` returns. Default ``"disv"``.
        persist : str, optional
            Stored on the spec and read by nothing -- ``project.save()`` pickles
            the object whatever this says. Default ``"pickle"``.

        Returns
        -------
        GridSpec
            A spec with ``method="object"``.
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

        Parameters
        ----------
        path : Path or str
            The pickle file, stored as a string. Not checked here; each
            :meth:`resolve` opens and unpickles it (relative to ``project_root``
            when given, else to the current working directory). Only unpickle a
            file you trust.
        name : str, optional
            The spec's name (library and lineage label). Default ``"grid"``.
        grid_type : str, optional
            A discretization label kept on the spec. Default ``"disv"``.

        Returns
        -------
        GridSpec
            A spec with ``method="pickle"``.
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
        crs: str | int | None = None,
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
            Path to a Python file containing the grid-builder function. Anchored
            to ``project_root`` at resolve when relative (the current working
            directory when none is given); a missing file raises
            ``FileNotFoundError`` there.
        function : str
            Name of the callable in ``script`` that returns the grid (required;
            an empty name raises ``ValueError`` here). It is called as
            ``function(project_root=<Path>, workspace=<Path>, spec=<this
            GridSpec>)`` and must not return ``None``.
        name : str, default "grid"
            Grid name; also names the builder's default workspace,
            ``<project_root>/_grid/<name>``.
        grid_type : str, default "disv"
            A label for the discretization the builder produces. Not checked
            against what it returns.
        inputs : Iterable, optional
            Input files the builder needs: ``DataSourceSpec`` objects or paths.
            Each one not marked ``external=True`` is checked to exist at resolve
            (anchored to ``project_root``), raising ``FileNotFoundError``. The
            builder reads them itself, as ``spec.inputs``.
        crs : str or int, optional
            Stored as ``spec.crs`` for the builder to read. Nothing applies it to
            the grid the builder returns.
        options : dict, optional
            Stored as ``spec.options``. NOT passed to the builder as keyword
            arguments -- the builder receives the whole spec and reads them from
            it.
        metadata : dict, optional
            Free-form provenance metadata, stored as ``spec.metadata`` and
            serialized with the spec.

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
        crs: str | int | None = None,
        **options: Any,
    ) -> GridSpec:
        """Load an existing grid from a GeoPackage. **Not implemented yet.**

        Planned API. :meth:`resolve` currently wires only generated Voronoi
        grids (:meth:`voronoi`) and Python grid builders (:meth:`python`), so a
        spec created here would raise at resolve time -- far from this call.
        This constructor therefore fails fast. To use an existing grid today,
        wrap a built grid object with :meth:`from_object`, or load a pickled
        grid with :meth:`from_pickle`.

        Parameters
        ----------
        path : Path or str
            The GeoPackage holding the grid's cell polygons. Unused: the call
            raises before reading anything.
        name : str, optional
            The planned grid name. Unused. Default ``"grid"``.
        layer : str, optional
            The planned GeoPackage layer holding the cells. Unused. Default
            ``None``.
        id_column : str, optional
            The planned column carrying each cell's number. Unused. Default
            ``None``.
        grid_type : str, optional
            The planned discretization label. Unused. Default ``"disv"``.
        crs : str or int, optional
            The planned CRS. Unused. Default ``None``.
        **options
            Accepted and discarded; the call raises before any is read.

        Raises
        ------
        NotImplementedError
            Always.
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
        crs: str | int | None = None,
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

        Parameters
        ----------
        name : str, optional
            The planned grid name. Unused. Default ``"grid"``.
        nlay : int
            The planned number of layers. Required by the signature, unused.
        nrow : int
            The planned number of rows. Required by the signature, unused.
        ncol : int
            The planned number of columns. Required by the signature, unused.
        delr : float
            The planned column width, in model length units. Required by the
            signature, unused.
        delc : float
            The planned row height, in model length units. Required by the
            signature, unused.
        top : float or array_like
            The planned model top elevation(s). Required by the signature,
            unused.
        botm : array_like
            The planned layer bottom elevations. Required by the signature,
            unused.
        crs : str or int, optional
            The planned CRS. Unused. Default ``None``.
        **options
            Accepted and discarded; the call raises before any is read.

        Raises
        ------
        NotImplementedError
            Always.
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
        crs: str | int | None = None,
        boundary_max_area: float | None = None,
        boundary_label: str | None = None,
        boundary_buffer: float | None = None,
        boundary_simplify_tolerance: float | None = None,
        boundary_densify_dist: float | None = None,
        refinement_max_area: float | None = None,
        refinement_priority: int | None = None,
        refinement_buffer: float | None = None,
        refinement_simplify_tolerance: float | None = None,
        refinement_densify_dist: float | None = None,
        breakline_max_area: float | None = None,
        breakline_priority: int | None = None,
        breakline_buffer: float | None = None,
        breakline_negative_buffer_after_clipping: float | None = None,
        breakline_simplify_tolerance: float | None = None,
        breakline_densify_dist: float | None = None,
        region_point_tolerance: float | None = None,
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
            Grid name -- the spec's, and the built grid's. An ``options`` or
            ``voronoi_options`` ``"name"`` key that differs raises ``TypeError``
            naming both, here and at :meth:`resolve`: the grid it built was named
            after the key while the spec kept this one (ledger 491). The
            same name given there too passes.
        refinement : DataSourceSpec, optional
            **One** polygon source. Every polygon in it becomes its own region
            with its own target cell area -- so "several refinement areas" means
            several features in one layer, not several sources. Passing a list
            here does not work.

            Size each one with ``fields={"area": "<column>"}`` on a
            :class:`~myflopy.sources.GeoPackageSourceSpec`; features with no
            value fall back to ``refinement_max_area``. One of the two is
            required. ``"label"`` and ``"priority"`` are mappable the same way.

            Polygons only, unless you set ``refinement_buffer``. A line in this
            source raises *"A LineString region has no area"*, and the buffer
            that would fix it applies to **every** feature in the source --
            measured, a 400x400 polygon in a buffered layer came out 67% larger.
            Send lines through ``breaklines`` instead.
        breaklines : Iterable[DataSourceSpec], optional
            Line sources buffered into refinement regions -- streams, faults,
            any corridor you want meshed finely along its length. Each line
            becomes one region, buffered by ``breakline_buffer`` (default 10,
            i.e. a corridor 20 units wide) and sized by its ``"area"`` field or
            ``breakline_max_area``.

            Unlike ``refinement`` this takes a LIST, which makes it the channel
            for several sources. Polygon sources are legal here too: with
            ``breakline_buffer=0`` a polygon passes through unbuffered, so N
            polygon geopackages can be one ``refinement`` plus the rest as
            breaklines. Mind the priority default below if you do that.
        points : Iterable[DataSourceSpec], optional
            Point sources pinned as fixed mesh vertices -- a well that must sit
            on a node, a gauge, a boundary corner. These are LOCATIONS, not
            sizes: a point contributes no cell-area target and no region. A
            non-point geometry contributes its ``representative_point()``.
        crs : str or int, optional
            Target CRS for the grid: every source is reprojected to it, and the
            built grid is labelled with it. Without it no source is reprojected;
            ``options["crs"]`` (or ``voronoi_options["crs"]``) then only labels
            the grid, and with none the grid is built as ``"EPSG:2927"`` (see
            :meth:`resolve`, Notes). Two of the three naming DIFFERENT CRSs raise
            ``TypeError`` naming both, here and at :meth:`resolve` -- with
            ``voronoi_options["crs"]`` the grid held ``crs=`` coordinates labelled
            as the other CRS, in silence (ledger 491); one CRS spelled twice
            (``2927`` and ``"EPSG:2927"``) passes.
        engine : str, default "triangle_voronoi_plus"
            Meshing engine identifier. It is the only one :meth:`resolve`
            supports; any other value raises ``NotImplementedError`` there.
        options : dict, optional
            Named meshing options as a dict -- the same keys as the named
            arguments under Other Parameters (and ``crs``,
            ``model_ws``/``workspace``), plus any stage key below, read when the
            stage dict does not set it. The named arguments are merged over it.
            A key nothing reads raises ``TypeError`` naming it, the nearest
            match and every accepted key -- here and again at :meth:`resolve`,
            where a saved spec arrives (ledger 417) -- and so do two
            spellings of one option: ``max_area``/``default_cell_area`` beside
            ``boundary_max_area``, a ``default_*_area`` beside its
            ``*_max_area``, ``line_buffer`` beside ``breakline_buffer``,
            ``angle`` beside ``min_angle``, ``workspace`` beside ``model_ws``.
            So do those spellings split between ``options`` and
            ``triangle_options`` -- ``options["min_angle"]`` beside
            ``triangle_options["angle"]``, ``options["model_ws"]`` beside
            ``triangle_options["workspace"]`` -- where the stage dict's won in
            silence (ledger 458).
            ``"idomain"`` raises ``TypeError`` (see ``**engine_options``), and so do
            the removed ``"protect_sources"``/``"protected_labels"`` (ledger 445).
        triangle_options : dict, optional
            Keyword arguments for ``TriangleGrid(...)`` -- ``angle``,
            ``maximum_area`` (the global cap, not ``max_area``), ``nodes``,
            ``exe_name``, ``additional_args``, plus
            ``model_ws``/``workspace`` for Triangle's working directory. A key
            also given through ``**engine_options`` takes that value; another
            SPELLING of one of these keys given there or in ``options``
            (``min_angle`` for ``angle``, ``model_ws`` for ``workspace``) raises
            ``TypeError``, here and at :meth:`resolve` (ledger 458). Not
            otherwise validated here; an unknown key raises ``TypeError`` at
            resolve.
        mesh_options : dict, optional
            Options for the mesh-build step (``profile``, ``optimize``,
            ``verbose``, ...). A key also given through ``**engine_options``
            takes that value. Any other key raises ``TypeError`` naming it and
            the accepted keys, here and at :meth:`resolve` (ledger 417); the removed
            ``protect_sources``/``protected_labels`` say so (ledger 445).
        voronoi_options : dict, optional
            Keyword arguments for ``VoronoiGridPlus(tri, ...)`` --
            ``idomain_path``, ``rasters``, ``qhull_options``, ``tessellation``,
            ``name``, ``crs``. A key also given through ``**engine_options`` takes
            that value. An ``"idomain"`` key raises ``TypeError``, as below; so do
            a ``"crs"`` that differs from ``crs=``/``options["crs"]`` and a
            ``"name"`` that differs from ``name=`` (ledger 491).
        **engine_options
            Stage options routed by name -- ``min_angle``, ``maximum_area``,
            ``profile``, ``idomain_path``, ``optimize`` and friends.
            ``min_angle`` is Triangle's ``angle``: giving both raises
            ``TypeError`` (the later keyword won, in silence), and so does
            ``min_angle=`` beside a ``triangle_options`` or ``options`` ``angle``
            (or ``angle=`` beside an ``options`` ``min_angle``) -- one replaced
            the other in silence (ledger 458).
            **``maximum_area`` is not ``max_area``**: ``maximum_area`` is
            Triangle's GLOBAL cap on every triangle's area, applied ON TOP of
            the boundary's and each refinement's own size, so the smaller limit
            wins (ledger 408); ``max_area`` is an alias of ``boundary_max_area``,
            the BACKGROUND cell size where no refinement applies. **An unknown
            keyword raises**, with the nearest match suggested; before that,
            a misspelled option was accepted, never read, and silently meshed at
            the default. ``tessellation="flopy"`` keeps FloPy's cells unrepaired
            (the default ``"dual"`` repairs them -- ledgers 316, 317; see
            :class:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus`).
            ``idomain=`` -- here, or as a key of ``options`` or
            ``voronoi_options`` -- raises ``TypeError`` naming
            ``mf.ModelContext(domain=...)``: it was routed to the grid and dropped
            there in silence, since a grid holds no active domain (ledger
            388).

        Other Parameters
        ----------------
        boundary_max_area : float, optional
            Target cell area for the background region, i.e. everywhere no
            refinement applies. In CRS units squared. Aliases
            ``default_cell_area`` / ``max_area`` reach the same option, and
            giving two of the three raises ``TypeError`` (only the first was
            read). Not ``maximum_area``, Triangle's global cap on every triangle,
            refinements included (see ``**engine_options``). One
            larger than the domain can hold -- no point of it more than
            ``0.5 * sqrt(max_area)`` from the boundary, e.g. ``>= H**2`` for a
            rectangle H wide -- gives a handful of cells and warns.
        boundary_label : str, default "domain"
            Name for the background region in diagnostics.
        boundary_buffer : float, default 0
            Buffer applied to the domain polygon before meshing.
        refinement_max_area : float, optional
            Target cell area for every ``refinement`` feature that does not
            carry its own via ``fields={"area": ...}``. One of the two is
            required, or resolution raises naming both. Alias
            ``default_refinement_area``; giving both raises ``TypeError``.
        refinement_priority : int, default 0
            Who claims the overlap where refinement regions cross others.
            Higher wins; a per-feature ``"priority"`` field beats this.
        refinement_buffer : float, default 0
            Buffer applied to **every** feature in the refinement source -- not
            only lines that need width. A polygon layer buffered to admit one
            line grows too; send lines through ``breaklines`` instead.
        breakline_max_area : float, optional
            Target cell area inside every breakline corridor. Alias
            ``default_breakline_area``; giving both raises ``TypeError``.
        breakline_priority : int, default 1
            Note this differs from ``refinement_priority`` (0), so a corridor
            takes the overlap from a refinement polygon by default.
        breakline_buffer : float, default 10
            Half-width of the corridor a line is buffered into, in CRS units --
            so the default corridor is 20 wide. ``0`` leaves a polygon source
            unbuffered, which is how ``breaklines`` doubles as the multi-source
            polygon channel. Alias ``line_buffer``; giving both raises
            ``TypeError``.
        breakline_negative_buffer_after_clipping : float, default 0
            Post-clip shrink; must be <= 0.
        boundary_simplify_tolerance, refinement_simplify_tolerance, breakline_simplify_tolerance : float, optional
            Douglas-Peucker tolerance applied per channel before meshing.
            Breaklines default to 10; the others to none.
        boundary_densify_dist, refinement_densify_dist, breakline_densify_dist : float, optional
            Maximum spacing to densify region outlines to, per channel.
        region_point_tolerance : float, optional
            Minimum spacing between the interior marker points Triangle needs,
            one per region. Auto-scaled from the extent when omitted.

        Returns
        -------
        GridSpec
            A durable grid recipe; call ``.resolve(workspace=...)`` for an eager grid,
            or pass it to ``mf.gwf(...).with_grid(...)`` / ``project.add_grid(...)``.

        Examples
        --------
        >>> spec = mf.GridSpec.voronoi(boundary=mf.ShapeSource("domain.shp"),
        ...                            refinement=mf.ShapeSource("wellfield.shp"),
        ...                            boundary_max_area=1.0e5,
        ...                            refinement_max_area=2.5e3, min_angle=30)
        >>> vor = spec.resolve(workspace="runs/_grid")

        Polygons and lines from one GeoPackage, each channel to its own layer,
        with a per-feature cell size on the polygons:

        >>> spec = mf.GridSpec.voronoi(
        ...     boundary=mf.ShapeSource("domain.gpkg", crs=2927),
        ...     refinement=mf.GeoPackageSourceSpec(
        ...         "refine.gpkg", layer="areas",
        ...         fields={"area": "max_area", "label": "name"}),
        ...     breaklines=[mf.GeoPackageSourceSpec("refine.gpkg", layer="creeks")],
        ...     breakline_max_area=800.0, breakline_buffer=40,
        ...     crs=2927, boundary_max_area=20_000.0,
        ... )

        Notes
        -----
        **Overlaps are resolved before meshing, by priority then size.** Regions
        sort by ``(-priority, area)``, and each claims only what higher-priority
        regions left it. So at equal priority a small polygon inside a big one
        claims first and the big one gets the donut -- nesting works. But the
        defaults are NOT equal: breaklines are priority 1 and refinement is 0,
        so a breakline corridor takes the overlap from a refinement polygon.
        Usually right for a stream; wrong if you used ``breaklines`` merely to
        pass a second polygon file, so set ``priority`` explicitly there.

        A region left with no unique area raises *"does not have any unique
        interior area left after resolving overlaps"* -- that is what a
        duplicated polygon gets you.

        **Inspect before meshing** with ``spec.resolve(build=False)``: each
        entry of ``tri._prepared_regions`` carries the ``label``, ``max_area``,
        ``source`` (``refinement`` / ``line`` / ``domain``) and both
        ``geometry`` and ``claim_geometry``, so the difference between the two
        is exactly what overlap resolution took away.
        """

        for label, given in (
            ("", engine_options),
            ("options=", options),
            ("voronoi_options=", voronoi_options),
        ):
            if given and "idomain" in given:
                where = f"{label}{{'idomain': ...}}" if label else "idomain="
                raise TypeError(
                    f"GridSpec.voronoi() takes no {where}: a grid holds no active "
                    "domain, and the option was routed to VoronoiGridPlus and dropped "
                    "there in silence (ledger 388). Give the model its active "
                    "domain with mf.ModelContext(domain=...) on mf.gwf(..., "
                    "context=...); to exclude cells from the grid's edge, use "
                    "idomain_path= (a vector file covering the inactive cells)."
                )
        boundary = _grid_source(boundary, "boundary")
        if refinement is not None:
            refinement = _grid_source(refinement, "refinement")
        breaklines = _grid_sources(breaklines, "breaklines")
        points = _grid_sources(points, "points")

        merged_options = {**({} if options is None else options)}
        # `None` means "not given" so the resolver's own default still applies --
        # inserting a mirrored default here would shadow the alias fallbacks
        # (`boundary_max_area` present-as-None would beat `default_cell_area`).
        for key, value in (
            ("boundary_max_area", boundary_max_area),
            ("boundary_label", boundary_label),
            ("boundary_buffer", boundary_buffer),
            ("boundary_simplify_tolerance", boundary_simplify_tolerance),
            ("boundary_densify_dist", boundary_densify_dist),
            ("refinement_max_area", refinement_max_area),
            ("refinement_priority", refinement_priority),
            ("refinement_buffer", refinement_buffer),
            ("refinement_simplify_tolerance", refinement_simplify_tolerance),
            ("refinement_densify_dist", refinement_densify_dist),
            ("breakline_max_area", breakline_max_area),
            ("breakline_priority", breakline_priority),
            ("breakline_buffer", breakline_buffer),
            ("breakline_negative_buffer_after_clipping",
             breakline_negative_buffer_after_clipping),
            ("breakline_simplify_tolerance", breakline_simplify_tolerance),
            ("breakline_densify_dist", breakline_densify_dist),
            ("region_point_tolerance", region_point_tolerance),
        ):
            if value is not None:
                merged_options[key] = value
        merged_mesh_options = {**({} if mesh_options is None else mesh_options)}
        merged_triangle_options = {
            **({} if triangle_options is None else triangle_options)
        }
        merged_voronoi_options = {
            **({} if voronoi_options is None else voronoi_options)
        }
        _check_grid_option_places(
            [("keyword", engine_options), ("options", merged_options),
             ("triangle_options", merged_triangle_options)],
            where="GridSpec.voronoi()",
        )
        if "min_angle" in engine_options and "angle" in engine_options:
            # Both write Triangle's `angle`; the later keyword won in silence.
            raise TypeError(
                "GridSpec.voronoi() got min_angle= and angle=, which are spellings of "
                "ONE option (Triangle's minimum angle); only the later keyword was "
                "ever read. Give one."
            )
        for key, value in engine_options.items():
            if key == "min_angle":
                merged_triangle_options["angle"] = value
            elif key in _GRID_TRIANGLE_OPTION_KEYS:
                merged_triangle_options[key] = value
            elif key in _GRID_MESH_OPTION_KEYS:
                merged_mesh_options[key] = value
            elif key in _GRID_VORONOI_ROUTED_KEYS:
                merged_voronoi_options[key] = value
            elif key in _VORONOI_OPTION_KEYS:
                merged_options[key] = value
            elif key in _GRID_REMOVED_OPTION_KEYS:
                raise TypeError(
                    f"GridSpec.voronoi() got {key}=, which "
                    f"{_GRID_REMOVED_OPTION_KEYS[key]}; delete it."
                )
            else:
                near = difflib.get_close_matches(
                    key, sorted(_VORONOI_OPTION_KEYS), n=1, cutoff=0.6
                )
                hint = f" Did you mean {near[0]!r}?" if near else ""
                raise TypeError(
                    f"GridSpec.voronoi() got an unknown option {key!r}.{hint} An "
                    f"unrecognised option used to be accepted and never read, so a "
                    f"typo silently meshed at the default."
                )
        _check_grid_option_keys(
            merged_options, merged_mesh_options, where="GridSpec.voronoi()"
        )
        _check_grid_identity(
            name=name, crs=crs, options=merged_options,
            voronoi_options=merged_voronoi_options, where="GridSpec.voronoi()",
        )
        return cls(
            name=name,
            grid_type="disv",
            method="voronoi",
            boundary=boundary,
            refinement=refinement,
            breaklines=breaklines,
            points=points,
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
            # list() first: this field is a tuple whose ELEMENT list is the
            # serialized shape, and ``from_dict`` rebuilds the tuple itself.
            # Passing the tuple in would tag the field as ``{"$tuple": [...]}``,
            # which the reader then walks as a mapping (``lineage`` already
            # does the same, :meth:`SimulationSpec.to_dict`).
            payload["inputs"] = _json_value(list(self.inputs))
        if self.triangle_options:
            payload["triangle_options"] = _json_value(self.triangle_options)
        if self.mesh_options:
            payload["mesh_options"] = _json_value(self.mesh_options)
        if self.voronoi_options:
            payload["voronoi_options"] = _json_value(self.voronoi_options)
        for key in ("source", "boundary", "refinement"):
            value = getattr(self, key)
            if value is not None:
                payload[key] = _json_value(value)
        if self.breaklines:
            payload["breaklines"] = [_json_value(source) for source in self.breaklines]
        if self.points:
            payload["points"] = [_json_value(source) for source in self.points]
        if self.crs is not None:
            # `crs` documents itself as taking `"EPSG:2927"` OR `2927`, and pandas
            # hands back `numpy.int64` for an integer column -- so a CRS read from a
            # config table is the ordinary producer of the value that broke the save.
            payload["crs"] = _json_value(self.crs)
        if self.engine is not None:
            payload["engine"] = self.engine
        if self.metadata:
            payload["metadata"] = _json_value(self.metadata)
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> GridSpec:
        """Recreate a grid specification from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A payload as :meth:`to_dict` writes it. ``"name"``, ``"grid_type"``
            and ``"method"`` are required (a missing one raises ``KeyError``);
            every other field is optional and defaults as on the constructor.
            Source entries are rebuilt into their ``DataSourceSpec`` classes.
            An ``"object"`` spec has no payload to read, since :meth:`to_dict`
            refuses to write one.
        """

        return cls(
            name=data["name"],
            grid_type=data["grid_type"],
            method=data["method"],
            pickle_path=data.get("pickle_path"),
            script=data.get("script"),
            function=data.get("function"),
            source=(
                _grid_data_source(data["source"])
                if data.get("source") is not None
                else None
            ),
            boundary=(
                _grid_data_source(data["boundary"])
                if data.get("boundary") is not None
                else None
            ),
            refinement=(
                _grid_data_source(data["refinement"])
                if data.get("refinement") is not None
                else None
            ),
            breaklines=tuple(
                _grid_data_source(item) for item in data.get("breaklines", ())
            ),
            points=tuple(_grid_data_source(item) for item in data.get("points", ())),
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
        """Build the grid this recipe describes, and return it.

        A :class:`GridSpec` is a durable *recipe*: it records where the geometry
        comes from and how to mesh it, but holds no grid. ``resolve`` is what
        turns it into an object you can hand to ``ModelContext(grid=...)``.
        Calling it yourself is the **eager** path; a :class:`Project` calls it
        for you at run-build time, which is why deferred specs exist.

        What comes back depends on how the spec was made. A
        :meth:`voronoi` spec builds a :class:`~myflopy.modflow.mf6.grid.triangle.TriangleGrid`
        from its boundary / refinement / breakline / point sources, meshes it,
        and wraps it in a ``VoronoiGridPlus``. :meth:`from_object` hands back the
        grid it already holds and :meth:`from_pickle` unpickles one, both
        ignoring every argument below except ``project_root`` (which anchors the
        pickle path). :meth:`python` imports the builder script and calls it.

        Parameters
        ----------
        project_root : Path or str, optional
            The directory that **relative paths in the recipe are anchored to**
            -- every source path (``boundary``, ``refinement``, ``breaklines``,
            ``points``), the ``idomain_path`` and ``rasters`` options, a
            :meth:`from_pickle` path, and a :meth:`python` spec's script and
            declared inputs. Absolute paths are used unchanged, and a source
            marked ``external=True`` is never anchored -- that flag exists for
            data that legitimately lives outside the project.

            Pass ``project_root=project.root`` so a spec written with paths
            relative to the project stays portable. Omitting it leaves relative
            source paths to resolve against the **current working directory**,
            which is fine for a script that never changes directory and a
            silent source of "file does not exist" the moment one does. (A
            :meth:`python` spec differs: it falls back to the CWD explicitly and
            passes the result to the builder as ``project_root=``.)
        workspace : Path or str, optional
            Where **Triangle's working files** are written -- the generated
            ``_triangle.0.poly`` / ``.node`` it is given, and the
            ``_triangle.1.node`` / ``.ele`` / ``.edge`` / ``.neigh`` / ``.poly``
            it produces. The directory is created if it does not exist. Reading
            ``_triangle.0.poly`` is the honest way to check what Triangle was
            actually asked for: its region section lists each refinement point
            and the maximum area attached to it.

            Resolution order, first hit wins: this argument, then
            ``triangle_options['model_ws'|'workspace']``, then
            ``options['model_ws'|'workspace']``. With none of them, meshing
            raises ``ValueError`` naming ``workspace=`` rather than writing into
            the current working directory; ``build=False`` writes nothing and
            is accepted, its ``TriangleGrid`` carrying ``<cwd>/_triangle`` for a
            mesh you build from it yourself. Give each grid its own directory --
            the file names are fixed, so two grids sharing a workspace overwrite
            each other's meshes.

            For a :meth:`python` spec this is instead the directory handed to
            the builder as ``workspace=`` (and created), defaulting to
            ``<project_root>/_grid/<name>``; with neither it nor
            ``project_root`` it raises ``ValueError``.
        build : bool, default True
            Whether to actually mesh. ``True`` runs Triangle and returns the
            finished ``VoronoiGridPlus``.

            ``False`` returns the **prepared but unmeshed** ``TriangleGrid``
            instead -- sources read, regions resolved against one another,
            interior points chosen, :meth:`~TriangleGrid.prepare` called, but
            Triangle never run. That is the cheap way to inspect what the
            recipe compiled to before paying for a mesh: ``tri.region_specs``
            is what you declared and ``tri._prepared_regions`` is what survived
            overlap resolution, each carrying the ``max_area``, ``priority`` and
            interior ``point`` Triangle will be given. Voronoi specs only; a
            :meth:`python` spec ignores it.
        return_triangle : bool, default False
            Return ``(vor, tri)`` rather than just ``vor``, keeping the
            ``TriangleGrid`` that produced the mesh -- useful for region
            diagnostics after the fact, since the ``VoronoiGridPlus`` does not
            carry them. Has no effect when ``build=False`` (the ``TriangleGrid``
            is already what you get), and raises for a :meth:`python` spec.

        Returns
        -------
        VoronoiGridPlus
            For a :meth:`voronoi` spec with ``build=True``. Tagged with
            ``grid.myflopy_grid_spec = spec``, so a built grid remembers the
            recipe it came from.
        tuple of (VoronoiGridPlus, TriangleGrid)
            The same, with ``return_triangle=True``.
        TriangleGrid
            For a :meth:`voronoi` spec with ``build=False`` -- prepared, unmeshed.
        Any
            For :meth:`from_object` (the held grid), :meth:`from_pickle` (the
            unpickled grid), or :meth:`python` (whatever the builder returned).

        Raises
        ------
        ValueError
            An ``object`` spec holding no grid, a ``pickle`` spec with no path,
            ``return_triangle=True`` on a :meth:`python` spec, or no directory
            named for the files resolving writes (see ``workspace``).
        NotImplementedError
            A method this does not wire. ``GridSpec.structured`` and
            ``GridSpec.from_geopackage`` construct but do not resolve -- wrap an
            already-built grid with :meth:`from_object` instead. Also raised for
            an ``engine`` other than ``'triangle_voronoi_plus'``.
        FileNotFoundError
            A :meth:`python` spec whose builder script does not exist.

        Notes
        -----
        **The CRS has a hardcoded fallback.** If neither the spec nor its
        options name one, the returned grid is built as ``"EPSG:2927"``
        (Washington State Plane South, feet). Pass ``crs=`` to
        :meth:`voronoi` rather than relying on it.

        Resolving is not cached -- each call re-reads the sources and re-runs
        Triangle, overwriting the files in ``workspace``. Hold the returned grid
        in a variable rather than calling ``resolve()`` twice.

        Examples
        --------
        >>> spec = mf.GridSpec.voronoi(
        ...     name="tentrails1",
        ...     boundary=mf.ShapeSource("domain.gpkg", crs=2927),
        ...     crs=2927,
        ...     boundary_max_area=20_000.0,
        ... )
        >>> vor = spec.resolve(project_root=proj.root,
        ...                    workspace=proj.root / "grids" / "base")

        Check the refinement before meshing:

        >>> tri = spec.resolve(project_root=proj.root, build=False)
        >>> [(r.label, r.max_area) for r in tri._prepared_regions]
        [('halfmoon', 400.0), ('domain_size', 20000.0)]

        Keep the ``TriangleGrid`` alongside the built grid:

        >>> vor, tri = spec.resolve(project_root=proj.root,
        ...                         workspace=ws, return_triangle=True)
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


#: Every key :mod:`myflopy.grid_spec_resolver` reads out of ``GridSpec.options``.
#: ``GridSpec.voronoi`` rejects anything else, so a misspelled option raises at the
#: call instead of silently falling back to a default. Kept in lockstep with the
#: resolver by ``test_gridspec_fail_fast.py``, which fails in BOTH directions.
_VORONOI_OPTION_KEYS: frozenset[str] = frozenset({
    "boundary_max_area", "default_cell_area", "max_area", "boundary_label",
    "boundary_buffer", "boundary_simplify_tolerance", "boundary_densify_dist",
    "refinement_max_area", "default_refinement_area", "refinement_priority",
    "refinement_buffer", "refinement_simplify_tolerance", "refinement_densify_dist",
    "breakline_max_area", "default_breakline_area", "breakline_priority",
    "breakline_buffer", "line_buffer", "breakline_negative_buffer_after_clipping",
    "breakline_simplify_tolerance", "breakline_densify_dist",
    "region_point_tolerance", "model_ws", "workspace", "crs", "profile",
})


#: The ``TriangleGrid(...)`` keys the resolver also takes out of ``GridSpec.options``
#: (``min_angle`` arrives as ``angle``), and the ``build_mesh`` keys a ``GridSpec``
#: carries in ``mesh_options`` -- or in ``options`` as a fallback. ONE copy each:
#: ``GridSpec.voronoi`` routes ``**engine_options`` by them and the resolver merges
#: by them (ledger 417).
_GRID_TRIANGLE_OPTION_KEYS: frozenset[str] = frozenset({
    "additional_args", "angle", "exe_name", "maximum_area", "nodes",
    "region_point_tolerance",
})
_GRID_MESH_OPTION_KEYS: frozenset[str] = frozenset({
    "cleanup", "damping", "min_feature_area", "min_move", "optimization_iterations",
    "optimize", "profile",
    "resample_domain_boundary", "resample_region_sources", "simplify_tolerance",
    "snap_tolerance", "target_segment_length", "verbose",
})
#: The ``VoronoiGridPlus(...)`` keys the resolver takes out of ``GridSpec.options``.
#: ``idomain`` is routed so a saved spec carrying it reaches the grid's own refusal
#: (ledger 388); ``GridSpec.voronoi`` refuses it at the call.
_GRID_VORONOI_ROUTED_KEYS: frozenset[str] = frozenset({
    "idomain", "idomain_path", "name", "qhull_options", "rasters", "tessellation",
})
#: Every key the voronoi route reads out of ``GridSpec.options``.
_GRID_OPTIONS_KEYS: frozenset[str] = (
    _VORONOI_OPTION_KEYS
    | {"min_angle"}
    | _GRID_TRIANGLE_OPTION_KEYS
    | _GRID_MESH_OPTION_KEYS
    | _GRID_VORONOI_ROUTED_KEYS
)
#: Option keys that USED to exist, each with why it is gone. A saved spec carrying
#: one is refused naming it rather than dropped (ledger 417).
_GRID_REMOVED_OPTION_KEYS: dict[str, str] = {
    "max_optimization_points": (
        "was removed (ledger 410): it was only echoed into the optimizer's report and "
        "never applied, so deleting it changes no mesh"
    ),
    "protect_sources": (
        "was removed (ledger 445): which vertices mesh optimization holds fixed is no "
        "longer configurable -- every refinement ring's, and those near a 'line' band's "
        "ring, as the old default held them"
    ),
    "protected_labels": (
        "was removed (ledger 445): which vertices mesh optimization holds fixed is no "
        "longer configurable -- every refinement ring's, and those near a 'line' band's "
        "ring, as the old default held them"
    ),
}
#: Spellings of ONE option, each group in the order the resolver reads them: it
#: takes the first present and ignored the rest, so giving two is refused (ledger
#: 417). ``max_area`` is the BACKGROUND cell size, not ``maximum_area`` --
#: Triangle's global cap. ``min_angle`` is read as Triangle's ``angle`` and wins over
#: an ``angle`` key in the same dict; ``model_ws`` wins over ``workspace``.
_GRID_OPTION_ALIASES: tuple[tuple[str, ...], ...] = (
    ("boundary_max_area", "default_cell_area", "max_area"),
    ("refinement_max_area", "default_refinement_area"),
    ("breakline_max_area", "default_breakline_area"),
    ("breakline_buffer", "line_buffer"),
    ("min_angle", "angle"),
    ("model_ws", "workspace"),
)


def _check_grid_option_keys(
    options: dict[str, Any], mesh_options: dict[str, Any], *, where: str
) -> None:
    """Refuse a key the voronoi route would never read, or two spellings of one option.

    ``where`` names the caller in the message: ``GridSpec.voronoi()`` at the call, or
    the spec's name at resolve -- where a SAVED spec arrives, so a key removed since,
    or one stored before keys were checked, is refused naming it rather than dropped
    in silence as it was (ledger 417).
    """

    for label, given, accepted in (
        ("options", options, _GRID_OPTIONS_KEYS),
        ("mesh_options", mesh_options, _GRID_MESH_OPTION_KEYS),
    ):
        unknown = sorted(str(key) for key in given if key not in accepted)
        if not unknown:
            continue
        # A removed key is explained even when another unknown key sorts first.
        removed = [name for name in unknown if name in _GRID_REMOVED_OPTION_KEYS]
        key = removed[0] if removed else unknown[0]
        if key in _GRID_REMOVED_OPTION_KEYS:
            reason = (
                f" {key!r} {_GRID_REMOVED_OPTION_KEYS[key]}; delete it from the "
                f"spec's {label}."
            )
        else:
            near = difflib.get_close_matches(key, sorted(accepted), n=1, cutoff=0.6)
            reason = f" Did you mean {near[0]!r}?" if near else ""
        raise TypeError(
            f"{where}: {label} holds {', '.join(map(repr, unknown))}, which the voronoi "
            f"route never reads, so it would be dropped in silence.{reason} "
            f"Accepted {label} keys: {', '.join(sorted(accepted))}."
        )
    for group in _GRID_OPTION_ALIASES:
        given = [name for name in group if name in options]
        if len(given) > 1:
            raise TypeError(
                f"{where}: options give {' and '.join(map(repr, given))}, which are "
                f"spellings of ONE option ({group[0]!r}); only {given[0]!r} was ever "
                "read. Give one."
            )


#: Spellings of one option that live in DIFFERENT places -- ``min_angle`` (a keyword or
#: an ``options`` key) is Triangle's ``angle`` (a ``triangle_options`` key), and
#: ``model_ws``/``workspace`` are read from ``triangle_options`` before ``options``.
_GRID_CROSS_PLACE_SPELLINGS: tuple[tuple[str, ...], ...] = (
    ("min_angle", "angle"),
    ("model_ws", "workspace"),
)


def _check_grid_option_places(places: list[tuple[str, dict[str, Any]]], *, where: str) -> None:
    """Refuse one option given under two spellings in two places (ledger 458).

    ``places`` is ``[(label, dict), ...]`` -- the keywords, ``options`` and
    ``triangle_options`` at :meth:`GridSpec.voronoi`, and the spec's ``options`` and
    ``triangle_options`` at resolve. ``GridSpec.voronoi(min_angle=30,
    triangle_options={"angle": 20})`` meshed at 30 and a spec holding
    ``options={"min_angle": 30}`` beside ``triangle_options={"angle": 20}`` at 20,
    each with no word; ``triangle_options`` holding both ``model_ws`` and
    ``workspace`` read the first. The SAME key in two places is not refused: a
    keyword is documented to win over the dict it lands in, and a stage dict over
    ``options``. Two spellings inside ``options``, or as two keywords, are refused
    elsewhere with their own message (:func:`_check_grid_option_keys`).
    """

    for group in _GRID_CROSS_PLACE_SPELLINGS:
        found = [
            (label, key) for label, given in places for key in group if key in given
        ]
        labels = {label for label, _ in found}
        if len({key for _, key in found}) < 2:
            continue
        if len(labels) == 1 and labels <= {"options", "keyword"}:
            # Refused by `_check_grid_option_keys` / `GridSpec.voronoi`, with their
            # own message naming the spelling that was read.
            continue
        named = " and ".join(
            f"{key}=" if label == "keyword" else f"{label}[{key!r}]"
            for label, key in found
        )
        raise TypeError(
            f"{where}: {named} are spellings of ONE option ({group[0]!r}); one "
            "silently replaced the other. Give one."
        )


def _check_grid_identity(
    *,
    name: str,
    crs: str | int | None,
    options: dict[str, Any],
    voronoi_options: dict[str, Any],
    where: str,
) -> None:
    """Refuse two different CRSs, or two different names, for one Voronoi grid.

    A grid's CRS can be spelled three ways -- ``crs=``, ``options["crs"]`` and
    ``voronoi_options["crs"]`` -- and so can its name (``name=``,
    ``options["name"]``, ``voronoi_options["name"]``). The resolver reprojects
    every source to ``crs=`` but LABELLED the grid from ``voronoi_options["crs"]``
    when it was there, so ``crs="EPSG:2927", voronoi_options={"crs":
    "EPSG:32610"}`` built a grid in 2927 coordinates labelled 32610, and
    ``options["crs"]`` beside ``crs=`` was ignored; ``voronoi_options["name"]``
    renamed the built grid away from the spec's ``name=``. Two DIFFERENT values
    are refused naming both (ledger 491, the two-spellings rule of ledger
    417); one value given twice -- ``2927`` and ``"EPSG:2927"`` included -- passes.
    Called by :meth:`GridSpec.voronoi` and again at resolve, where a spec built
    directly or loaded from a saved project arrives.
    """

    crs_given = [
        (label, value)
        for label, value in (
            ("crs", crs),
            ("options['crs']", options.get("crs")),
            ("voronoi_options['crs']", voronoi_options.get("crs")),
        )
        if value is not None
    ]
    for label, value in crs_given[1:]:
        first_label, first = crs_given[0]
        if not _same_crs(first, value):
            raise TypeError(
                f"{where}: {first_label}={first!r} and {label}={value!r} name two CRSs "
                "for one grid. Every source is reprojected to crs= while "
                "voronoi_options['crs'] (else options['crs']) labels the built grid, "
                "so the grid would hold coordinates in one CRS labelled as another "
                "(or one of them was ignored). Give the grid's CRS once, as crs=."
            )
    names_given = [
        (label, value)
        for label, value in (
            ("options['name']", options.get("name")),
            ("voronoi_options['name']", voronoi_options.get("name")),
        )
        if value is not None
    ]
    for label, value in names_given:
        if str(value) != str(name):
            raise TypeError(
                f"{where}: name={name!r} and {label}={value!r} name one grid twice; "
                f"the spec would be {name!r} and the grid it builds {value!r}. Give "
                "the name once, as name=."
            )


def _grid_source(value: Any, argument: str) -> DataSourceSpec:
    """One grid source, or a ``TypeError`` naming what to do instead.

    Every rejected form here used to reach the resolver and die as
    ``AttributeError: 'X' object has no attribute 'path'``, which names neither
    the argument at fault nor the fix.
    """

    if isinstance(value, DataSourceSpec):
        return value
    if isinstance(value, (list, tuple, set, frozenset)):
        raise TypeError(
            f"{argument}= takes ONE source, not a {type(value).__name__}. Every "
            f"FEATURE in that one source becomes its own region, so several "
            f"refinement areas belong in one layer. For several FILES: merge "
            f"them into one layer carrying an 'area' column, or pass the extras "
            f"through breaklines=[...] with breakline_buffer=0 (which leaves a "
            f"polygon unbuffered)."
        )
    extra = ""
    if type(value).__name__ in {"GeoDataFrame", "DataFrame", "GeoSeries"}:
        extra = (
            f" A {type(value).__name__} is not file-backed, and a GridSpec is a "
            f"serializable recipe -- it has to survive to_dict()/from_yaml(). "
            f"Write it out first, e.g. "
            f"gdf.to_file(path, layer='areas', driver='GPKG'), then point a "
            f"GeoPackageSourceSpec at it."
        )
    raise TypeError(
        f"{argument}= takes a file-backed source spec (ShapeSource, "
        f"GeoPackageSourceSpec or TableSource), got {type(value).__name__}.{extra}"
    )


def _grid_sources(value: Any, argument: str) -> tuple[DataSourceSpec, ...]:
    """A tuple of grid sources from an iterable argument (``breaklines``/``points``)."""

    if isinstance(value, DataSourceSpec):
        raise TypeError(
            f"{argument}= takes a SEQUENCE of sources, not one source -- write "
            f"{argument}=[source]. (Unlike refinement=, this one is a list.)"
        )
    try:
        items = list(value)
    except TypeError:
        raise TypeError(
            f"{argument}= takes a sequence of source specs, got "
            f"{type(value).__name__}."
        ) from None
    return tuple(_grid_source(item, argument) for item in items)


def _grid_data_source(data: dict[str, Any]) -> DataSourceSpec:
    """Rebuild a grid data source, rejecting non-file-backed kinds."""

    source = source_from_dict(data)
    if not isinstance(source, DataSourceSpec):
        raise ValueError(
            f"Grid sources must be file-backed data sources, got {type(source).__name__}"
        )
    return source


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
    Given one, they resolve cells eagerly, so build the grid first and put it here.
    ``mf.sfr``/``mf.lak`` always need it; the ``.gpkg`` forms, ``mf.rch``/``mf.evt``
    and ``mf.uzf`` can instead be declared without ``context=`` and resolve at build
    against the model's own context.

    Parameters
    ----------
    grid : VoronoiGridPlus or grid object, optional
        The built grid features are mapped onto -- normally a ``VoronoiGridPlus``
        (from ``GridSpec.voronoi(...).resolve(...)`` or built directly). Read by
        the eager helpers, the build-time grid checks, :meth:`ModelSpec.coverage_on`
        and ``with_grid(remap=True)``. A model declared with ``grid=`` (a
        :class:`GridSpec`/:class:`GridRef`) has this REPLACED by the resolved grid
        at build. Not serialized: :meth:`ModelSpec.to_dict` keeps only its
        connectivity fingerprint (``"grid_pin"`` in ``metadata``). Default
        ``None``.
    surfaces : pandas.DataFrame or LayerBuildResult, optional
        Per-cell layer elevations read by the surface-aware builders (SFR reach
        tops, LAK connections, ``CellSurfaceOffset``): a table shaped like
        ``vor.gdf_topbtm`` (geometry first, then model top, then each layer
        bottom), or a ``LayerBuildResult``, which is turned into that table. When
        ``None`` the grid's live ``gdf_topbtm`` is used; when given it wins over
        the grid's table, and a disagreement between the two is WARNED about.
        Not serialized. Default ``None``.
    dates : object, optional
        Stress-period datetimes, carried for post-build hooks and user code.
        Nothing in the library reads it, and neither :meth:`ModelSpec.to_dict`
        nor a run's frozen snapshot keeps it. Default ``None``.
    domain : array_like of int, optional
        The IDOMAIN: shaped ``(nlay, ncpl)``, or a flat ``(ncpl,)`` array meaning
        the same activity in every layer. Cells ``> 0`` are active. The helpers
        pick and validate cells against it, and the build refuses a domain whose
        cell count (or, when 2-D, layer count) disagrees with the built DISV.
        ``None`` leaves activity to the DISV / grid. Not serialized. Default
        ``None``.
    length_units : str or None, optional
        The model's length unit, in any spelling a ``LayerStack`` takes
        (``"feet"``/``"foot"``/``"ft"``/``"us-ft"``/``"usft"``,
        ``"meters"``/``"meter"``/``"metres"``/``"m"``,
        ``"centimeters"``/``"centimeter"``/``"cm"``, any case and surrounding
        spaces ignored), or ``None`` / ``""`` / ``"unknown"`` for no declared
        unit. Anything else -- a non-``str`` included -- raises ``ValueError``
        HERE, at construction (ledger 433). Read by
        ``mf.sfr``/``mf.lak`` to derive MF6's ``LENGTH_CONVERSION`` when none is
        given (3.28081 for feet, 1.0 for metres, 100.0 for centimetres, 1.0 with
        no unit), by ``mf.lak`` to convert a lakebed ``Surface`` declaring
        ``units=`` (with no unit declared here such a surface is REFUSED), and
        for flux labels in build warnings. It does NOT set the DISV's
        ``length_units``. Serialized. Default ``"feet"``.
    time_units : str or None, optional
        The model's time unit: ``"seconds"``/``"second"``/``"sec"``/``"s"``,
        ``"minutes"``/``"minute"``/``"min"``, ``"hours"``/``"hour"``/``"hr"``/
        ``"h"``, ``"days"``/``"day"``/``"d"`` or ``"years"``/``"year"``/``"yr"``,
        any case and surrounding spaces ignored, or ``None`` / ``""`` /
        ``"unknown"`` for no declared unit. Anything else -- a non-``str``
        included -- raises ``ValueError`` HERE, at construction (ledger
        456). Read by ``mf.sfr``/``mf.lak`` to derive MF6's
        ``TIME_CONVERSION`` when none is given (86400.0 for days, 3600.0 for
        hours, 60.0 for minutes, 1.0 for seconds, 31557600.0 for years, 1.0 with
        no unit) and for flux labels in build warnings. It does NOT set TDIS's
        ``time_units``. Serialized. Default ``"days"``.
    metadata : dict, optional
        Free-form metadata, serialized with the model and visible to post-build
        hooks through the context. :meth:`ModelSpec.to_dict` adds the grid's
        ``"grid_pin"`` to the copy it writes. Default: empty.

    Examples
    --------
    >>> vor = mf.GridSpec.voronoi(...).resolve(workspace="runs/_grid")  # or a VoronoiGridPlus
    >>> layers = stack.build(attach=True)
    >>> ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)
    >>> flow = mf.gwf("flow", context=ctx, packages=[...])
    """

    grid: Any = None
    surfaces: Any = None
    dates: Any = None
    domain: Any = None
    length_units: str | None = "feet"
    time_units: str | None = "days"
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Refuse a ``length_units`` (ledger 433) or ``time_units`` (ledger 456) no consumer can read."""

        check_model_length_units(self.length_units)
        check_model_time_units(self.time_units)

#: The ``attrs`` key and value a lazy ``mf.disv`` stamps on the elevation table it
#: publishes onto the model's context: that table IS the built discretization's, so a
#: grid table disagreeing with it is the one out of date (``attach=False`` leaves the
#: grid's table alone on purpose), and ``_layer_elevations`` does not warn about it.
_PUBLISHED_BY_KEY, _PUBLISHED_BY_DISV = "published_by", "disv"


def _layer_elevations(context: Any) -> Any:
    """The ONE per-cell layer-elevation table every surface-aware builder reads.

    Three builders asked three ways. SFR read ``grid.gdf_topbtm`` directly and
    never looked at ``context.surfaces``; LAK and ``CellSurfaceOffset`` read
    ``context.surfaces`` first. ``ModelContext(surfaces=vor.gdf_topbtm)`` captures
    the frame OBJECT, so re-running a layer-build cell in a notebook -- which
    re-publishes ``vor.gdf_topbtm`` -- left two packages in one model on different
    elevations. Measured: SFR reach top 70.0 against ``CellSurfaceOffset`` 60.0,
    truth 30.0, no exception, normal termination. The only route found to wrong
    numbers in an MF6 input file.

    ``context.surfaces`` wins when given, because it is what the caller said and
    there is no way to tell a deliberate pin from a stale snapshot. What changes
    is that every builder now agrees on it, and a DISAGREEMENT with the grid's
    live table is announced rather than written into the model -- except for a table
    a lazy ``mf.disv`` published (``attrs[_PUBLISHED_BY_KEY]``): that one is the built
    discretization's own, so with ``attach=False`` it is the GRID's table that is out of
    date, and the warning's advice -- rebuild the context -- would put the old
    elevations back (ledger 323).

    A ``LayerBuildResult`` is accepted and turned into the table -- it is a
    documented ``surfaces`` value, and LAK used to die on it with
    ``'LayerBuildResult' object has no attribute 'drop'``.
    """

    surfaces = context.surfaces
    live = getattr(context.grid, "gdf_topbtm", None)
    if surfaces is None:
        return live
    if not hasattr(surfaces, "columns") and hasattr(surfaces, "botm"):
        surfaces = _elevation_table_from_layers(surfaces)
    published = getattr(surfaces, "attrs", None)
    published = isinstance(published, dict) and published.get(_PUBLISHED_BY_KEY) == _PUBLISHED_BY_DISV
    if published:
        if live is not None and live is not surfaces and _elevations_disagree(surfaces, live):
            logger.debug(
                "the grid's gdf_topbtm disagrees with the table a lazy DISV published "
                "(attach=False keeps the grid's); the DISV's is used, and is the model's"
            )
    elif live is not None and live is not surfaces and _elevations_disagree(surfaces, live):
        # Through `warnings.warn` as well as the logger: `import myflopy` installs a
        # NullHandler and root has none, so the log line alone reached NOBODY -- 0 bytes
        # on stderr (ledger 288; the excavation design's D11, ledger 323).
        message = (
            "ModelContext.surfaces and the grid's gdf_topbtm disagree -- the context "
            "was most likely built before the layers were re-published (a notebook "
            "cell re-run, or a stack excavated after the context was made). Every "
            "surface-aware package in this model uses context.surfaces; rebuild the "
            "ModelContext if the grid's table is the one you meant."
        )
        logger.warning(message)
        warnings.warn(message, UserWarning, stacklevel=3)
    return surfaces


def _layer_elevation_source(context: Any) -> tuple[Any, dict[str, Any]]:
    """The elevation table :func:`_layer_elevations` returns AND the excavation records
    that ride with it, from the SAME place, so the two can never disagree.

    ``mf.lak`` chooses each lake's connection rule from an excavation record (the
    excavation design's F6, ledger 325), and a record read from somewhere other than
    the elevations it describes is the one way the rule and the geometry could part.
    So this follows ``_layer_elevations``' precedence exactly:

    * ``context.surfaces`` given as a ``LayerBuildResult`` -- its ``excavations``;
    * given as a table -- the ``attrs["excavations"]`` a lazy DISV (or a thawed run)
      published on it, or, when the table IS the grid's own live ``gdf_topbtm``
      object, the grid's ``excavations``; any other bare table carries none;
    * not given -- the grid's ``gdf_topbtm`` and the grid's ``excavations``, which
      ``LayerBuildResult.attach_to_grid`` sets beside it (and the ``gdf_topbtm``
      setter clears, so another writer leaves no stale record).

    Parameters
    ----------
    context : ModelContext
        The model's context (anything with ``surfaces`` and ``grid``).

    Returns
    -------
    tuple of (pandas.DataFrame or None, dict of str to ExcavationRecord)
        The table exactly as :func:`_layer_elevations` returns it (its warnings
        included), and a copy of the records riding with it -- ``{}`` when none do.
    """

    table = _layer_elevations(context)
    surfaces = getattr(context, "surfaces", None)
    grid = getattr(context, "grid", None)
    grid_records = getattr(grid, "excavations", None) or {}
    if surfaces is None:
        return table, dict(grid_records) if table is not None else {}
    if not hasattr(surfaces, "columns") and hasattr(surfaces, "botm"):
        return table, dict(getattr(surfaces, "excavations", None) or {})
    attrs = getattr(surfaces, "attrs", None)
    published = attrs.get("excavations") if isinstance(attrs, dict) else None
    if published:
        return table, dict(published)
    if surfaces is getattr(grid, "gdf_topbtm", None):
        return table, dict(grid_records)
    return table, {}


def _elevation_table_from_layers(layers: Any) -> Any:
    """``{0: top, 1..nlay: bottoms}`` from a ``LayerBuildResult``, no geometry.

    The same column contract ``LayerBuildResult.attach_to_grid`` publishes, built
    without publishing it: resolving a builder's input must not mutate the grid.
    """

    import numpy as np
    import pandas as pd

    columns = {0: np.asarray(layers.top, dtype=float).reshape(-1)}
    for index, bottom in enumerate(np.asarray(layers.botm, dtype=float)):
        columns[index + 1] = bottom.reshape(-1)
    return pd.DataFrame(columns)


def _elevations_disagree(first: Any, second: Any) -> bool:
    """Whether two elevation tables carry different numbers, geometry ignored."""

    import numpy as np

    try:
        a = first.select_dtypes("number").to_numpy(dtype=float)
        b = second.select_dtypes("number").to_numpy(dtype=float)
    except (AttributeError, TypeError, ValueError):  # not both tables -- cannot compare
        return False
    return a.shape != b.shape or not np.allclose(a, b, equal_nan=True)


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

    Parameters
    ----------
    project_root : Path or str, optional
        The directory relative paths are anchored to: a model's deferred
        :class:`GridSpec` sources (passed as ``GridSpec.resolve(project_root=)``)
        and the relative source paths inside every package's ``inputs`` recipe.
        ``None`` leaves relative paths to the current working directory. Stored
        as a ``Path``. Default ``None``.
    simulation_workspace : Path or str, optional
        The simulation's run directory. When ``grid_workspace`` is unset,
        deferred grids are built under ``<simulation_workspace>/_grid/<model>``
        -- and with this unset too, under the FloPy simulation's own
        ``sim_path`` the same way. :meth:`SimulationSpec.build_flopy` fills it
        from its workspace when unset. Stored as a ``Path``. Default ``None``.
    grid_workspace : Path or str, optional
        The directory deferred grids are built under, one subdirectory per
        model (``<grid_workspace>/<model>``); takes precedence over
        ``simulation_workspace``. Stored as a ``Path``. Default ``None``.
    package_specs : dict of str to PackageSpec, optional
        The project package library, keyed by library key: what each
        :class:`PackageRef` resolves to. A ref whose key is absent raises
        ``KeyError`` at build. Copied. Default: empty.
    grid_specs : dict of str to GridSpec, optional
        The project grid library, keyed by library key: what each
        :class:`GridRef` resolves to. A ref whose key is absent raises
        ``KeyError`` at build. Copied. Default: empty.
    build_grids : bool, optional
        Passed to ``GridSpec.resolve(build=...)`` for each model's deferred grid:
        ``False`` resolves a Voronoi spec to its prepared, UNMESHED
        ``TriangleGrid`` rather than a built grid. Default ``True``.
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
        """Return a copy with a simulation workspace when one is known.

        Parameters
        ----------
        workspace : Path or str or None
            The simulation workspace to record. Ignored -- ``self`` is returned
            unchanged -- when it is ``None`` or when this context already has a
            ``simulation_workspace``, so an explicit one is never overwritten.
        """

        if workspace is None or self.simulation_workspace is not None:
            return self
        return replace(self, simulation_workspace=Path(workspace))

    def model_grid_workspace(self, model_name: str) -> Path | None:
        """Return the generated grid workspace for one model, when available.

        Parameters
        ----------
        model_name : str
            The model whose deferred grid is being built.

        Returns
        -------
        Path or None
            ``<grid_workspace>/<model_name>`` when ``grid_workspace`` is set,
            else ``<simulation_workspace>/_grid/<model_name>``, else ``None``.
            Nothing is created on disk here.
        """

        if self.grid_workspace is not None:
            return Path(self.grid_workspace) / model_name
        if self.simulation_workspace is not None:
            return Path(self.simulation_workspace) / "_grid" / model_name
        return None

    def package_spec(self, key: str) -> PackageSpec:
        """Return a resolved project-level package spec by key.

        Parameters
        ----------
        key : str
            A package-library key, as a :class:`PackageRef` carries it. Matched
            exactly against ``package_specs``; an absent key raises
            ``KeyError``.
        """

        try:
            return self.package_specs[key]
        except KeyError as error:
            raise KeyError(
                f"Package reference '{key}' is unresolved in the build context."
            ) from error

    def grid_spec(self, key: str) -> GridSpec:
        """Return a resolved project-level grid spec by key.

        Parameters
        ----------
        key : str
            A grid-library key, as a :class:`GridRef` carries it. Matched exactly
            against ``grid_specs``; an absent key raises ``KeyError``.
        """

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

    Parameters
    ----------
    name : str
        Key under which the hook's result is recorded in
        ``BuiltModel.hook_results``. Must be unique among one model's hooks
        (``ModelSpec`` raises ``ValueError`` on a duplicate).
    callback : callable
        Called as ``callback(model, packages, context)`` after every package is
        built and the build-time grid checks have passed; hooks run in
        declaration order. Its return value is stored under ``name``. Must be an
        importable module-level function to round-trip: a lambda (at module
        level too) or one defined inside another function runs but raises
        ``ValueError`` in :meth:`to_dict`.
    """

    name: str
    callback: Hook

    def run(self, model: Any, packages: dict[str, Any], context: ModelContext) -> Any:
        """Run the callback with the completed FloPy model and its context.

        Parameters
        ----------
        model : flopy.mf6.MFModel
            The built FloPy model.
        packages : dict of str to object
            The built FloPy packages keyed by ``PackageSpec.name`` (disabled
            packages absent).
        context : ModelContext
            The model's context as the build left it -- after a deferred grid
            was resolved into it and a lazy DISV published its domain and
            surfaces.

        Returns
        -------
        object
            Whatever ``callback`` returns.
        """

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
        """Recreate a hook from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            ``{"name": ..., "callback": "<module>:<qualname>"}`` as
            :meth:`to_dict` writes it. Both keys are required (``KeyError``);
            the callback is imported, so one that no longer imports raises
            ``ImportError``/``AttributeError``.
        """

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
    model : flopy.mf6.ModflowGwf or ModflowGwt or ModflowGwe or ModflowPrt
        The built FloPy model object. It also carries the context as
        ``model.myflopy_context``.
    packages : dict of str to object
        Built package objects keyed by ``PackageSpec.name`` -- for an ``mf.ref``
        entry, the resolved spec's name, not the library key. A disabled
        package is absent.
    context : ModelContext
        The :class:`ModelContext` (grid/domain/surfaces/dates) for this model as
        the build left it: a deferred grid resolved into it, and a lazy DISV's
        domain and surfaces published on it.
    hook_results : dict of str to object
        Values returned by post-build hooks, keyed by hook name.
    simulation : flopy.mf6.MFSimulation or None
        The parent FloPy simulation, when known; :attr:`sim` falls back to the
        model's own ``sim``.
    resolved : dict of str to PackageSpec
        What a LAZY package actually built, as an eager ``PackageSpec``, keyed
        by package name -- and also by the ``mf.ref`` key it was declared
        under. Filled only by builders that resolve at build time -- the lazy
        ``mf.disv``, ``.gpkg``, ``rch``/``evt`` (list and array) and ``uzf``
        forms and ``mf.Polygons`` arrays; a run freezes these specs rather than
        the recipe.
    """

    model: Mf6Model
    packages: dict[str, Any]
    context: ModelContext = field(default_factory=ModelContext)
    hook_results: dict[str, Any] = field(default_factory=dict)
    simulation: Mf6Simulation | None = None
    #: What a LAZY package actually built, as an eager ``PackageSpec``, keyed by
    #: package name. Populated only by builders that resolve something at build time
    #: (today: the lazy DISV). It is what lets a run FREEZE the arrays rather than the
    #: recipe -- freezing the declared spec would store a `LayerStack` whose rasters
    #: the reopened run then has to go back and read, which is the very dependence the
    #: snapshot exists to remove (ledger 251, 270).
    resolved: dict[str, Any] = field(default_factory=dict)

    @property
    def sim(self) -> Mf6Simulation:
        """Return the parent FloPy simulation for this built model."""

        if self.simulation is not None:
            return self.simulation
        model_simulation = getattr(self.model, "sim", None)
        if isinstance(model_simulation, Mf6Simulation):  # type: ignore[misc]  # flopy is untyped
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
        """Return one built package by name.

        Parameters
        ----------
        name : str
            A key of :attr:`packages` -- the ``PackageSpec.name``, matched
            exactly (no slot-prefix matching, unlike :meth:`ModelSpec.package`).
            An absent name raises ``KeyError``.
        """

        try:
            return self.packages[name]
        except KeyError as error:
            raise KeyError(f"Built model has no package named '{name}'.") from error


























@dataclass(frozen=True, slots=True)
class ModelSpec:
    """Definition of one GWF, GWT, GWE, or PRT model.

    Package names are unique within the model. Calling :meth:`with_packages`
    replaces an existing package with the same name, which makes scenario
    variants explicit and inexpensive to create. Build one with ``mf.gwf`` /
    ``mf.gwt`` / ``mf.gwe`` / ``mf.prt`` rather than directly.

    Parameters
    ----------
    name : str
        The model name: passed to the FloPy model as ``modelname=``, and the
        name :meth:`SimulationSpec.model`, ``on=`` and an exchange's ``models``
        refer to. Must be unique within a simulation.
    model_type : ModelType or str
        ``"gwf"``, ``"gwt"``, ``"gwe"`` or ``"prt"`` (lowercase; any other value
        raises ``ValueError``). Coerced to :class:`ModelType`; picks the FloPy
        model class when ``builder`` is ``None``.
    packages : sequence of PackageSpec, PackageRef or str, optional
        The model's packages, built in this order. A bare ``str`` becomes a
        :class:`PackageRef` to that library key. Names must be unique
        (``ValueError``); when every entry is a concrete :class:`PackageSpec`,
        each one's ``requires`` is checked here too. Stored as a tuple.
        Default: empty.
    options : dict, optional
        Keyword arguments for the model constructor, called as
        ``builder(simulation, modelname=name, **options)`` -- for example
        ``newtonoptions``. Not checked here; FloPy refuses an unknown one at
        build (``FlopyException: Extraneous kwargs``). Default: empty.
    builder : callable, optional
        The model constructor. ``None`` uses the FloPy class for
        ``model_type`` (``flopy.mf6.ModflowGwf`` etc.); a custom one must accept
        ``(simulation, modelname=..., **options)`` and return the FloPy model,
        and must be importable for :meth:`to_dict`. Default ``None``.
    context : ModelContext, optional
        The geometry the model is built against (grid, domain, surfaces).
        Default: an empty :class:`ModelContext`.
    hooks : sequence of PostBuildHook, optional
        Callbacks run after the model and its packages are built, in order.
        Names must be unique (``ValueError``). Stored as a tuple. Default:
        empty.
    grid : GridSpec, GridRef or grid object, optional
        A DEFERRED grid, resolved at every :meth:`build` (a :class:`GridRef`
        through the build context's grid library) and put into the built
        model's context as ``grid``, replacing ``context.grid``. A bare built
        grid is wrapped with :meth:`GridSpec.from_object`. ``None`` builds on
        ``context`` as given. Default ``None``.
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
        # NOTE: after this point ``model_type`` is always the enum; use
        # ``type_enum`` where a typed value is needed.
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

    def with_packages(
        self, *packages: PackageEntry | str | Iterable[Any], slot: str | None = None
    ) -> ModelSpec:
        """Return a copy with each package added or replaced in its slot.

        One verb, one or many, varargs or a list::

            model.with_packages(mf.npf(k=50.0))
            model.with_packages(npf, ic, oc)
            model.with_packages([npf, ic, oc])

        The slot defaults to each entry's own (a concrete package's name, or a
        reference's key prefix). ``slot=`` targets one explicitly -- needed when
        a semantically-keyed ref such as ``mf.ref("k/calibrated")`` fills the
        ``npf`` slot, since its key prefix says ``k``.

        **With ``slot=`` the slot must already exist**, and that assertion is
        load-bearing rather than pedantic: without it a swap that misses silently
        ADDS, and MF6 accepts two RCH packages on one model and sums them --
        measured at 0.8000 against an intended 0.4000, under normal termination.
        Without ``slot=`` this adds or replaces freely, which is what declaring a
        package means.

        Parameters
        ----------
        *packages : PackageSpec, PackageRef, str, or iterables of them
            The entries to add or replace, applied in order. A list or tuple is
            flattened one level; a ``str`` is one :class:`PackageRef` key, never
            split into characters. Each one replaces the entry its slot matches
            -- by full name first, then by a ref's key prefix -- or is appended
            when nothing matches. A slot matching several entries raises
            ``KeyError``.
        slot : str, optional
            The existing entry to replace, named by its full name or slot. Only
            valid with exactly ONE package (``ValueError`` otherwise), and it
            must match an existing entry (``KeyError`` otherwise). Default
            ``None``: each entry's own slot.

        Returns
        -------
        ModelSpec
            The copy; this spec is unchanged.
        """

        entries = [_package_entry(entry) for entry in _flatten_packages(packages)]
        if slot is not None and len(entries) != 1:
            raise ValueError(
                f"slot='{slot}' targets ONE package; {len(entries)} were given. "
                "Call with_packages once per targeted slot."
            )
        spec = self
        for entry in entries:
            target = slot if slot is not None else _package_slot(entry)
            current = list(spec.packages)
            index = _package_entry_index(
                current, target, where=f"Model '{spec.name}'"
            )
            if index is None:
                if slot is not None:
                    raise KeyError(
                        f"Model '{spec.name}' has no '{target}' package to replace. "
                        f"It carries: {', '.join(spec.package_names) or 'nothing'}. "
                        "Drop slot= to add it instead."
                    )
                current.append(entry)
            else:
                current[index] = entry
            spec = replace(spec, packages=tuple(current))
        return spec

    def package(self, name: str) -> PackageEntry:
        """Return one package specification by name.

        Resolves through the shared lookup, so an ambiguous name raises here
        rather than quietly handing back whichever was declared first.

        Parameters
        ----------
        name : str
            A package's full name (a ref's full key), else a slot: a concrete
            package's name or the first ``/`` segment of a ref's key. Disabled
            packages are found too. No match, or several on the slot, raises
            ``KeyError``.

        Returns
        -------
        PackageSpec or PackageRef
            The entry as declared -- a ref is not resolved.
        """

        index = _package_entry_index(
            list(self.packages), name, where=f"Model '{self.name}'"
        )
        if index is None:
            raise KeyError(
                f"Model '{self.name}' has no package named '{name}'. "
                f"It carries: {', '.join(self.package_names) or 'nothing'}."
            )
        return self.packages[index]

    def without_packages(self, *names: str | Iterable[str]) -> ModelSpec:
        """Return a copy without the named packages. One or many.

        RAISES when a name matches nothing. It used to be a silent no-op that
        still returned a new object, so a typo -- or a package already removed --
        read exactly like success, and the only way to know was to compare
        ``package_names`` before and after.

        Removes exactly ONE entry per name. Matching the slot alone used to
        delete every entry sharing it, so ``"wel"`` on a model carrying
        ``wel/north`` and ``wel/south`` took both; now that names two packages
        and raises.

        Parameters
        ----------
        *names : str or iterables of str
            The packages to remove, each matched as :meth:`package` matches it.
            Lists are flattened one level. A name matching nothing, or several
            entries by slot, raises ``KeyError``.
        """

        spec = self
        for name in _flatten_packages(names):
            current = list(spec.packages)
            index = _package_entry_index(current, name, where=f"Model '{spec.name}'")
            if index is None:
                raise KeyError(
                    f"Model '{spec.name}' has no package '{name}' to remove. "
                    f"It carries: {', '.join(spec.package_names) or 'nothing'}."
                )
            del current[index]
            spec = replace(spec, packages=tuple(current))
        return spec

    def resolved_packages(
        self,
        build_context: SpecBuildContext | None = None,
    ) -> tuple[PackageSpec, ...]:
        """Return concrete package specs, resolving project package refs.

        Parameters
        ----------
        build_context : SpecBuildContext, optional
            Supplies the package library each :class:`PackageRef` is looked up
            in (``package_specs``). Required only when the model holds a ref:
            without one a ref raises ``ValueError``, and a key missing from the
            library raises ``KeyError``. Default ``None``.

        Returns
        -------
        tuple of PackageSpec
            In declaration order, disabled ones included, after checking that
            names are unique and every ``requires`` is present, enabled and
            declared earlier (``ValueError`` otherwise).
        """

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
            # A ref `with_grid` stamped hands its stamp to the package it resolves, so
            # the build refuses it as it refuses a stamped inline package.
            packages.append(_ref_resolved(package, build_context.package_spec(package.key)))
        _validate_concrete_packages(packages)
        return tuple(packages)

    def with_options(self, **overrides: Any) -> ModelSpec:
        """Return a copy with selected model options added or replaced.

        Parameters
        ----------
        **overrides
            Model-constructor keyword arguments merged over ``options`` key by
            key (``newtonoptions=...``). Not checked here; FloPy refuses an
            unknown one at build.
        """

        return replace(self, options={**self.options, **overrides})

    def with_context(self, context: ModelContext) -> ModelSpec:
        """Return a copy carrying the supplied model context.

        Parameters
        ----------
        context : ModelContext
            Replaces the whole context -- grid, domain, surfaces, units and
            metadata together. Packages already declared eagerly against the old
            context keep what they resolved; nothing is re-resolved or checked
            until build.
        """

        return replace(self, context=context)

    def with_grid(
        self,
        grid: GridSpec | GridRef | Any,
        *,
        grid_free: tuple[str, ...] | list[str] = (),
        remap: bool = False,
        rule: str | None = None,
        uncovered: str | None = None,
        library: Any = None,
    ) -> ModelSpec:
        """Return a copy on ``grid``, with everything derived from the old one dropped.

        ``context.grid``/``domain``/``surfaces`` all describe the OLD mesh, and keeping
        them is how a grid swap used to build clean on cell numbers that meant something
        else (ledger 249). They are dropped -- the same three
        :data:`~myflopy.grid_identity.DROPPED_CONTEXT_FIELDS` a save cannot carry.

        Packages are NOT re-resolved. Each is classified instead, and a package that
        cannot be verified on the new mesh is STAMPED so the build refuses it:

        - a lazy package (a recipe in ``inputs``) rebuilds itself -- nothing to do;
        - a PINNED package is refused by the existing pin check, naming its old mesh;
        - a package naming cells with NO pin -- hand-written cellids, a ``.flopy``
          form -- can be told apart from a correct one by nothing at all, since every
          number is still in range on a larger mesh. It is stamped and refused.

        ``grid_free=`` names packages the caller asserts are mesh-independent; their
        stamp is cleared rather than applied. Read :meth:`move_report` first.

        **An ``mf.ref`` is classified as the library package it NAMES**, which needs
        ``library=project``: a ref carries no options, so classified as itself it read
        ``unchanged`` and a library ``mf.ic(strt=<array>)`` built on the old mesh's values.
        A ref that cannot be verified is stamped ON THE REF -- the library package is
        shared and stays untouched -- and the build refuses it. One the move must CHANGE
        (arrays ``remap=True`` carries, an eager ``.gpkg`` re-declared) becomes an INLINE
        copy in place of the ref, since the changed package is no longer the library's.
        A model holding an unstamped ref and no ``library=`` raises naming it.

        **The mesh measured against is the one the packages were DECLARED on** --
        ``grid=``, else ``context.grid`` -- and a DEFERRED one (a recipe, a pickle, an
        ``mf.grid_ref``) is resolved here to measure it, as the build would: a
        :meth:`GridSpec.from_object` costs nothing, a recipe is meshed once in a
        temporary directory. Before that a deferred-grid model had nothing to measure,
        and its per-cell arrays and typed cellids moved as ``unchanged``.

        ``remap=True`` CARRIES a package's per-cell arrays onto the new mesh by area
        overlay instead of stamping them (ledger 300), turning ``unverifiable`` into
        ``remapped``. Opt-in, not the default, for a measured reason: where a package has
        a SOURCE, re-declaring it on the new mesh is more faithful than interpolating it
        -- 7.1% worst-case boundary-flux error against this engine's 12.5% at best -- so
        the remap is for an array with no recipe, such as a calibrated K field. It serves
        ``npf`` ``k``/``k22``/``k33`` and the intensive arrays ``ic.strt``/``sto.ss``/
        ``sto.sy``; an integer flag and anything extensive are refused by name.

        ``rule=`` is the conductivity rule and exists because the right average for K
        depends on the FLOW DIRECTION -- a property of the solution, not of the field --
        so nothing declared per-field can be right everywhere. Measured in MF6 on one 2:1
        coarsening: with the contrast ACROSS the flow the arithmetic mean gave 264x the
        fine-mesh boundary flux; rotated 90 degrees it was exact and the harmonic mean was
        wrong by about 250x. ``"geometric"`` is the default because it is never best and
        never catastrophic, and the remap WARNS naming the cells where the arithmetic and
        harmonic means diverge -- so the choice is surfaced where it matters rather than
        needing to be known in advance.

        ``uncovered=`` says what a new cell with no old cell under it means:
        ``"nearest"`` fills it from the nearest old cell and says so in detail,
        ``"raise"`` refuses.

        Parameters
        ----------
        grid : GridSpec, GridRef or grid object
            The new mesh, stored as the model's ``grid`` (a bare built grid is
            wrapped with :meth:`GridSpec.from_object`) and resolved at build.
            ``remap=True`` needs real cell polygons, so there it must be a built
            grid or a :class:`GridSpec` (resolved on the spot in a temporary
            directory, without a project root, and meshed again at build); a
            :class:`GridRef` raises ``TypeError``.
        grid_free : sequence of str, optional
            Package names (a ref's library key) the caller asserts are
            mesh-independent: they read ``unchanged`` and are not stamped.
            Default: empty.
        remap : bool, optional
            Carry per-cell arrays (``npf`` ``k``/``k22``/``k33``, ``ic.strt``,
            ``sto.ss``/``sto.sy``) from the mesh the model was declared on --
            ``grid=``, else ``context.grid``, resolved if deferred -- onto
            ``grid`` by area overlay; a pinned, pin-voided or already-stamped
            package is left alone. A ref'd package it carries becomes an inline
            copy (its ``move_report`` row says so). Raises ``ValueError`` when
            the model holds no grid at all, and -- once a package has an array
            to carry -- when the two meshes are in different frames (CRS, or no
            overlap), or the array is an integer flag or one with no measured
            rule. Default ``False``.
        rule : str or None, optional
            The averaging rule for the conductivities (``k``/``k22``/``k33``)
            under ``remap=True``: ``"geometric"``, ``"mean"``, ``"harmonic"`` or
            ``"majority"``; ``ss``/``sy`` always take the area-weighted mean, and
            ``strt`` -- a head field -- is carried linear between cell centres, as
            ``model.hds.remap`` carries one (ledger 600). Any other value raises
            ``ValueError`` at the call, whether or
            not a package has an array to carry. Given without ``remap=True`` it
            raises ``ValueError`` naming ``remap=True``. Default ``None``:
            ``"geometric"`` under ``remap=True``.
        uncovered : str or None, optional
            For ``remap=True``, what a new cell no old cell reaches gets:
            ``"nearest"`` (the nearest old cell's value, reported) or
            ``"raise"`` (``ValueError``). Checked at the call like ``rule``, and
            likewise refused without ``remap=True``. Default ``None``:
            ``"nearest"`` under ``remap=True``.
        library : Project, optional
            The project whose libraries the model's ``mf.ref`` and
            ``mf.grid_ref`` entries name: each ref is classified as the package
            it names, a ``grid_ref`` grid is resolved to measure against, and
            the project root anchors a grid recipe's relative paths. Required
            when the model's grid is a ``grid_ref``, and when the model is on a
            mesh and holds a ref no earlier move stamped -- ``ValueError``
            naming it otherwise. A model on no mesh at all (the first
            ``with_grid`` of ``mf.gwf(...)``) needs none: nothing could be
            stamped. Anything without ``package``/``grid`` methods raises
            ``TypeError``. Nothing in the library is changed. Default ``None``.

        Returns
        -------
        ModelSpec
            The copy, with ``context.grid``, ``domain`` and ``surfaces`` set to
            ``None`` and the packages re-lazified, remapped or stamped as
            described above (a stamped ref carries ``moved_from``).

        Raises
        ------
        ValueError
            For a ref or ``grid_ref`` with no ``library=``, and as ``remap``,
            ``rule`` and ``uncovered`` describe.
        TypeError
            For a ``library=`` that is not a project, or a ``grid`` that
            ``remap=True`` cannot take polygons from.
        KeyError
            For a ref or ``grid_ref`` key the library does not hold.
        """

        rule, uncovered = _remap_choices(remap, rule, uncovered)
        library = _checked_library(library, "with_grid")
        free = tuple(grid_free)
        # The OLD mesh, still in hand: `with_grid` returns a copy and does not mutate,
        # so the grid these packages were declared on is present at this call -- and a
        # deferred one is resolved now, as the build would. It is what lets
        # `package_move` see a per-cell array that nothing stamped (ledger 295).
        old = _declared_mesh(self, library, verb="with_grid")
        declared_on = getattr(old, "ncpl", None)
        # Each entry beside what it RESOLVES to: a ref is classified as the library
        # package it names, never as itself (refused without `library=`). A model on
        # NO mesh has no pin to stamp against, so a ref's verdict could change nothing
        # and the library is not demanded -- this is also the spelling that gives a
        # model its first grid.
        pairs = _move_entries(self, library, verb="with_grid", required=old is not None)
        # An EAGERLY declared `.gpkg` is re-lazified from the recipe it recorded, so it
        # moves EXACTLY rather than being refused (ledger 296). Done first, because it
        # turns a `stale` package into a `rebuilt` one and the stamping below must see
        # the result, not the original.
        relazified = tuple(_relazified(package) for _entry, package in pairs)
        if remap:
            relazified = self._remap_onto(
                relazified, grid, old_grid=old, rule=rule, uncovered=uncovered
            )
        # `package_moves`, not `package_move`: a MOVER's staleness lives in its siblings,
        # and one package is not enough to see it (ledger 303). All three callers share
        # this function so the pre-build report, the stamp and `Run.resolution` cannot
        # disagree about the same mover.
        verdicts = package_moves(
            relazified, grid_free=_grid_free_names(pairs, free), ncpl=declared_on
        )
        pin = grid_pin(old)
        moved_entries = []
        for (entry, resolved), package in zip(pairs, relazified, strict=True):
            # `grid_free` is honoured inside `package_move`, which answers
            # "unchanged" for a named package -- so testing it again here was dead
            # code that read as load-bearing. Mutation found it.
            unverifiable = (
                pin is not None
                and verdicts.get(getattr(package, "name", None), ("", ""))[0]
                == "unverifiable"
                # A VOIDED package is already refused by its own complaint, which names
                # the mesh the edit invalidated. Stamping it too made `with_grid` list
                # the same package twice with two different explanations -- measured on
                # a real ncpl 36 -> 100 move. One defect, one sentence.
                and PIN_VOIDED not in (package.metadata or {})
            )
            if isinstance(entry, PackageRef) and package is resolved:
                # The library package is shared, so the stamp rides on the REF and the
                # build hands it across (`_ref_resolved`). An earlier move's stamp is
                # kept: the values still mean the mesh they were declared on.
                if unverifiable and entry.moved_from is None:
                    entry = replace(entry, moved_from=pin)
                moved_entries.append(entry)
                continue
            if isinstance(entry, PackageRef):
                # The move CHANGED what the ref names -- arrays remapped, an eager
                # `.gpkg` re-declared -- and the changed copy is not the library's
                # package, so it goes INLINE in place of the ref (user decision,
                # 2026-10-15), saying where it came from. Not "via a ref" any more.
                metadata = {
                    key: value for key, value in (package.metadata or {}).items()
                    if key != MOVED_VIA_REF
                }
                package = replace(
                    package, metadata={**metadata, INLINED_FROM_REF: entry.key}
                )
            moved_entries.append(_move_stamped(package, pin) if unverifiable else package)
        return replace(
            self,
            grid=grid,
            packages=tuple(moved_entries),
            context=replace(self.context, grid=None, domain=None, surfaces=None),
        )

    def _remap_onto(self, packages, grid, *, old_grid, rule: str, uncovered: str):
        """Carry every remappable package's arrays onto ``grid``; see ``with_grid``.

        ``old_grid`` is the mesh the arrays were declared on, read by ``with_grid`` at
        the call and not later, because ``with_grid`` does not mutate: that grid is
        present at this call and absent from the result. That is the whole reason the
        remap is eager rather than a recipe resolved at build.
        """

        if old_grid is None:
            raise ValueError(
                f"{self.name}: remap=True carries this model's per-cell arrays from the "
                "mesh they were declared on onto the new one, and this model holds no "
                "grid to carry them from -- neither `grid=` nor "
                "`context=ModelContext(grid=...)`, or one whose route does not resolve. "
                "Move it without `remap=`, or declare the grid it was built on."
            )
        old_mesh = _mesh_parts(old_grid, verb="with_grid(remap=True)")
        new_mesh = _mesh_parts(grid, verb="with_grid(remap=True)")
        old_polys, new_polys = old_mesh.gdf_vorPolys, new_mesh.gdf_vorPolys
        carried = []
        for package in packages:
            if isinstance(package, PackageSpec):
                package, _ = _remapped_package(
                    package,
                    old_polys=old_polys,
                    new_polys=new_polys,
                    old_grid=old_mesh,
                    new_grid=new_mesh,
                    ncpl=getattr(old_grid, "ncpl", None),
                    domain=getattr(self.context, "domain", None),
                    rule=rule,
                    uncovered=uncovered,
                )
            carried.append(package)
        return tuple(carried)

    def coverage_on(self, grid: Any) -> Any:
        """How much of ``grid`` this model's mesh has NO data under, in model units.

        Read it BEFORE `with_grid`, because this is the question `move_report` cannot
        answer: that one says what happens to each PACKAGE, and this says what the new
        mesh would have no data for at all. It is pure geometry and knows nothing about
        packages, so it is the same number whatever is being moved.

        Reported as AREA in the model's own units rather than as a coverage fraction
        (user decision, ledger 298): a fraction needs a tolerance before it means
        anything, and on a perfectly nesting refinement `coverage < 1` flags a third of
        the cells as short by 4e-15 of themselves.

        The model's mesh is its ``grid=<built grid>``, else ``context.grid``. Raises
        rather than guessing when it holds neither: a DEFERRED-grid model has nothing to
        overlay, and answering "nothing is covered" would be indistinguishable from a
        real total miss.

        Parameters
        ----------
        grid : grid object or GridSpec
            The candidate new mesh: a built grid carrying ``gdf_vorPolys``, or a
            :class:`GridSpec`, which is resolved here in a temporary directory
            removed on return (no project root, so relative sources resolve
            against the current working directory). A :class:`GridRef` raises
            ``TypeError`` -- pass ``project.grid(key)``.

        Returns
        -------
        CoverageReport
            ``total_area``, ``uncovered_area`` and per-cell ``uncovered`` areas
            in the model's units squared, plus ``disjoint`` for two meshes that
            share no area at all.

        Raises
        ------
        ValueError
            When the model holds no built grid -- neither ``grid=<built grid>``
            nor ``context.grid``.
        """

        # DEFERRED, and measured worth it: hoisting this works, but it pushes `specs`
        # from import layer 1 to 2 and 54 modules with it -- and it would make a bare
        # `import myflopy.specs` pull geopandas, which almost everything imports and
        # almost nothing needs geometry from.
        from myflopy.modflow.mf6.grid.geometry import cell_coverage

        old = _mesh_in_hand(self)
        if old is None:
            raise ValueError(
                f"Model {self.name!r} carries no grid to measure coverage against -- "
                "neither a built `grid=` nor `context=ModelContext(grid=...)`. Declare "
                "one, or resolve its GridSpec first -- a deferred grid has no mesh "
                "until build."
            )
        return cell_coverage(old.gdf_vorPolys, _polygons_of(grid, verb="coverage_on"))

    def move_report(
        self,
        *,
        grid_free: tuple[str, ...] | list[str] = (),
        library: Any = None,
    ) -> list[dict[str, str]]:
        """What a grid move would do to each package, readable BEFORE the build.

        The build is what refuses, so this is how you find out what it will refuse
        without waiting for `prepare_run` to raise. It classifies exactly as
        :meth:`with_grid` does -- an ``mf.ref`` as the library package it names, against
        the mesh the packages were declared on, a deferred one resolved -- so the two
        cannot disagree; which is why it needs ``library=`` in the same cases.

        Parameters
        ----------
        grid_free : sequence of str, optional
            Package names (a ref's library key) to report as mesh-independent,
            as :meth:`with_grid` would treat them. Default: empty.
        library : Project, optional
            The project whose libraries the model's ``mf.ref`` and
            ``mf.grid_ref`` entries name, as for :meth:`with_grid`. Required
            when the model holds a ref no earlier move stamped, or a
            ``grid_ref`` grid -- ``ValueError`` naming it otherwise. Default
            ``None``.

        Returns
        -------
        list of dict
            One ``{"package", "how", "why"}`` row per package -- a ref's row
            under its library key -- where ``how`` is ``"rebuilt"``, ``"stale"``,
            ``"unverifiable"``, ``"remapped"`` or ``"unchanged"``.
        """

        library = _checked_library(library, "move_report")
        pairs = _move_entries(self, library, verb="move_report")
        old = _declared_mesh(self, library, verb="move_report")
        return move_report(
            [(_package_entry_label(entry), resolved) for entry, resolved in pairs],
            grid_free=_grid_free_names(pairs, tuple(grid_free)),
            ncpl=getattr(old, "ncpl", None),
        )

    def build(
        self,
        simulation: Any,
        *,
        build_context: SpecBuildContext | None = None,
    ) -> BuiltModel:
        """Build the model and all enabled packages in declaration order.

        Resolves a deferred ``grid`` first, then constructs the FloPy model,
        builds each enabled package, checks the result against the grid it was
        built on, and runs the post-build hooks.

        Parameters
        ----------
        simulation : flopy.mf6.MFSimulation
            The FloPy simulation the model is added to, passed to the model
            constructor as its first argument.
        build_context : SpecBuildContext, optional
            The project libraries and locations: needed to resolve an
            ``mf.ref``/``mf.grid_ref`` entry, to anchor relative source paths to
            the project root, and to place a deferred grid's working files --
            when it names no directory for them, the FloPy simulation's own
            ``<sim_path>/_grid/<model>``. ``None`` uses an empty one, so a ref
            raises. Default ``None``.

        Returns
        -------
        BuiltModel

        Raises
        ------
        ValueError
            When a package or the context disagrees with the grid it was built
            against, all listed at once and grouped by what each is about -- the
            mesh (DISV against the context grid, a grid pin, a move stamp), the
            context's domain shape, a cellid no cell answers to (negative, or past
            the grid), a cell the discretization has removed, a lake connection it
            contradicts -- each group under its own heading and advice, read on a
            DIS, DISV or DISU through its own cellid form (ledger 469); for a
            package whose cellids are not that form (a DISV pair on a DIS), where
            FloPy stopped with ``IndexError`` -- unless the package's own builder
            refuses it first, as ``mf.hfb``'s does; and for a package ref with no ``build_context`` on a model without
            a deferred ``grid``.
        KeyError
            For a ref whose key is not in the build context's library -- which
            includes every ref when ``build_context`` is ``None`` and the model
            has a deferred ``grid``, since an empty context is used then.
        """

        context = self.context
        if self.grid is not None:
            build_context = (
                SpecBuildContext() if build_context is None else build_context
            )
            grid_spec = self.grid
            if isinstance(grid_spec, GridRef):
                grid_spec = build_context.grid_spec(grid_spec.key)
            grid_workspace = build_context.model_grid_workspace(self.name)
            if grid_workspace is None and getattr(simulation, "sim_path", None) is not None:
                # No run directory was named, so the simulation's OWN -- FloPy's
                # `sim_ws`, the cwd by its default -- under the documented
                # `<simulation_workspace>/_grid/<model>`. `GridSpec.resolve` refuses to
                # choose a directory itself (ledger 599); it used to choose
                # `<cwd>/_triangle`, shared by every model.
                grid_workspace = Path(simulation.sim_path) / "_grid" / self.name
            grid = grid_spec.resolve(
                project_root=build_context.project_root,
                workspace=grid_workspace,
                build=build_context.build_grids,
            )
            context = replace(context, grid=grid)

        builder = self.builder or _DEFAULT_MODEL_BUILDERS[self.type_enum]
        model = builder(simulation, modelname=self.name, **self.options)
        model.myflopy_context = context
        packages = _anchored_packages(self.resolved_packages(build_context), build_context)
        built_packages = {}
        for package in packages:
            # A cellid of another form than the built discretization's -- a DISV pair on
            # a DIS, a DIS triple on a DISV -- is refused by name (ledger 469).
            # FloPy stops on a SHORTER one with `IndexError: tuple index out of range`,
            # which is translated here; a LONGER one builds, FloPy writing only its
            # leading parts, and the complaints below refuse it before anything is
            # written. Not asked before the build: a builder with its own refusal --
            # `mf.hfb`'s, naming `pairs=` in the grid's form -- speaks first.
            try:
                built = package.build(model)
            except IndexError as error:
                wrong = cellid_form_complaints(model, [package])
                if not wrong:
                    raise
                built_grid = _built_grid(model)
                raise ValueError(
                    _build_refusal(
                        self.name, {"cellid": wrong}, built_grid.name if built_grid else "DISV"
                    )
                ) from error
            if built is not None:
                built_packages[package.name] = built
        # RE-READ the context: a lazy DISV resolves the layer stack against the grid
        # and publishes the domain and surfaces it produced. `context` above is the
        # DECLARED one, and passing it on from here condemns a correct lazy build --
        # measured, the guard raised "context domain is shaped (2, 28) against
        # ncpl=92" whether or not the DISV had just published a correct (2, 92) one
        # (ledger 270).
        context = getattr(model, "myflopy_context", context)
        resolved_specs = dict(getattr(model, "_myflopy_resolved", {}) or {})
        # ALSO key a resolved spec by the `mf.ref` KEY it was declared under, because
        # that is the name the frozen snapshot will look it up by. `freeze._for` reads
        # `getattr(package, "name")` off the DECLARED entry, and a `PackageRef`'s name
        # is its library key ("chd/west"), while a builder only ever knows the resolved
        # spec's own name ("chd"). The mismatch is COMPLETELY silent -- the snapshot
        # just keeps the recipe, and the run dies on reopen looking for a file that has
        # moved. `mf.ref` is the documented library path, so this is the common case
        # rather than an edge (ledger 280).
        for entry, spec in zip(self.packages, packages, strict=True):
            if isinstance(entry, PackageRef) and spec.name in resolved_specs:
                resolved_specs[entry.key] = resolved_specs[spec.name]
        # Before the hooks, so that a stale grid is reported by name rather than as
        # whatever a hook happens to raise first. A hook attaches FloPy objects
        # directly to the model rather than as specs, so cellids it adds are outside
        # what this can read -- the grid pin is what will cover those.
        live = [p for p in packages if p.name in built_packages]
        # `built_packages=`: pair each spec with the FloPy object it BUILT, not by name -- a bare
        # PackageSpec with no pname is named `chd_0` by FloPy and went unchecked.
        complaints = _discretization_complaint_families(
            model, context, live, built_packages=built_packages
        )
        complaints["mesh"] += moved_package_complaints(live)
        # One paragraph per kind of complaint, each saying what it is about: a negative
        # cellid or a removed cell is not "resolved on another mesh" (ledger 427).
        built_grid = _built_grid(model)
        refusal = _build_refusal(self.name, complaints, built_grid.name if built_grid else "DISV")
        if refusal is not None:
            raise ValueError(refusal)
        # A face MF6 cannot compute once the file is written: zero length, or both cell
        # centres on it -- NaN heads under normal termination, or a division by zero (ledger
        # 316). Its own refusal, since no package is at fault.
        disv = model.get_package("disv") if hasattr(model, "get_package") else None
        flat = unwritable_faces(disv, getattr(model, "simulation_data", None)) if disv is not None else []
        if flat:
            listing = "; ".join(f"cells {a} and {b} at ({x:.10g}, {y:.10g}): {why}" for a, b, x, y, why in flat[:5])
            if len(flat) > 5:
                listing += f"; and {len(flat) - 5} more"
            raise ValueError(
                f"Model {self.name!r}'s DISV would be written with {len(flat)} face(s) MF6 cannot "
                f"compute -- {listing}. FloPy writes a coordinate above 1e5 to nine significant "
                "digits (0.01 ft at x = 1.2e6), and MF6 returns NaN for every head across a "
                "zero-length face under normal termination, or stops on one whose cell centres "
                "both lie on it. The grid has a feature narrower than the file resolves there: a "
                "boundary or region vertex a hair from another, a hairline notch or crack in the "
                "boundary. Remove or simplify it (profiles 'balanced' and 'clean' simplify the "
                "boundary by 2 ft), or build the model in local coordinates. A grid rebuilt from a "
                "snapshot or a pickle is the grid that ran: re-resolve it."
            )
        # The documented opt-outs the guard does not refuse: an array-form flux put in
        # one layer everywhere (`irch=None`) where that layer's cell is removed. Legal,
        # and MF6 drops the column's flux without a word, so it is SAID (ledger 322) --
        # through `warnings.warn` as well as the logger, which alone reaches nobody
        # (ledger 288).
        for message in active_cell_warnings(model, context, live, built_packages=built_packages):
            message = f"Model {self.name!r}: {message}"
            logger.warning(message)
            warnings.warn(message, UserWarning, stacklevel=2)
        # An excavated lake's build-time warnings: no NEWTON, walls into icelltype 0,
        # recharge/ET switched off under a wet lake, a transparent bed over uneven
        # cells (ledger 325). Each has a legitimate use, so none refuses.
        for message in lak_build_warnings(model, live):
            message = f"Model {self.name!r}: {message}"
            logger.warning(message)
            warnings.warn(message, UserWarning, stacklevel=2)
        hook_results = {
            hook.name: hook.run(model, built_packages, context) for hook in self.hooks
        }
        return BuiltModel(
            model=model,
            packages=built_packages,
            context=context,
            hook_results=hook_results,
            simulation=simulation,
            resolved=resolved_specs,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a durable JSON-ready representation of this model recipe.

        Post-build hooks are serialized by importable reference; a hook whose
        callback is a lambda or closure raises via :func:`_callable_ref`.
        """

        payload: dict[str, Any] = {
            "kind": "ModelSpec",
            "name": self.name,
            "model_type": self.type_enum.value,
            "options": _json_value(self.options),
            "packages": [package.to_dict() for package in self.packages],
            "context": {
                "metadata": _json_value(_context_metadata(self.context)),
                # Two plain strings: the refusal to serialize grids and domains
                # (ledger 193 -- 227 and this comment both cited 180, which ledger
                # 194's renumbering had already moved) never covered them. They went
                # with the arrays, and a metres/seconds model came back as feet/days
                # (ledger 227).
                "length_units": self.context.length_units,
                "time_units": self.context.time_units,
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
        """Recreate a model specification from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A payload as :meth:`to_dict` writes it. ``"name"`` and
            ``"model_type"`` are required (``KeyError``). Packages, options,
            builder, hooks and grid are rebuilt when present. The context comes
            back with only its ``metadata`` (including the recorded
            ``"grid_pin"``) and units -- grid, domain and surfaces were never
            written -- and a payload without units reads as feet/days.
        """

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
                # A file saved before 227 has neither; it meant the defaults.
                length_units=context_data.get("length_units", "feet"),
                time_units=context_data.get("time_units", "days"),
            ),
            hooks=tuple(
                PostBuildHook.from_dict(hook) for hook in data.get("hooks", ())
            ),
            grid=_grid_from_dict(data.get("grid")),
        )

    @property
    def type_enum(self) -> ModelType:
        """``model_type`` as the enum (coerced in ``__post_init__``)."""

        return ModelType(self.model_type)

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
            f"model_type={self.type_enum.value!r}, "
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
            f"<p><strong>Type:</strong> {escape(self.type_enum.value)} "
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
    Build one with ``mf.exchange(a, b)`` (or its ``.gwfgwt``/``.gwfgwe``/
    ``.gwfprt``/``.gwfgwf`` forms) rather than directly.

    Parameters
    ----------
    name : str
        The exchange's name, unique among a simulation's exchanges
        (``ValueError`` otherwise); the key of ``BuiltSimulation.exchanges`` and
        what :meth:`SimulationSpec.with_exchange` replaces by.
    builder : callable
        Called at build as ``builder(simulation, models, **options)``, where
        ``models`` is the tuple of built FloPy models; normally one of the
        ``myflopy.builders.build_*_exchange`` functions. The simulation checks
        that a GWF-GWT/GWE/PRT exchange's models share a mesh by this builder's
        ``__name__``. Must be importable for :meth:`to_dict`.
    models : sequence of str
        Names of the coupled models, in the order the builder expects (for the
        ``build_*`` functions, the flow model first). At least two, or
        ``ValueError``; a name no built model carries raises ``KeyError`` at
        build. Stored as a tuple.
    options : dict, optional
        Keyword arguments for the builder, which the ``build_*`` functions pass
        on to the FloPy exchange class. Not checked here. Default: empty.
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
        """Build the exchange after resolving its model names.

        Parameters
        ----------
        simulation : flopy.mf6.MFSimulation
            The FloPy simulation the exchange is added to, passed to the builder
            first.
        built_models : dict of str to BuiltModel
            The simulation's built models keyed by model name; each name in
            ``models`` is looked up here (``KeyError`` naming the missing ones)
            and its FloPy ``model`` passed on.

        Returns
        -------
        object
            Whatever the builder returns -- normally the FloPy exchange package.
        """

        missing = [name for name in self.models if name not in built_models]
        if missing:
            raise KeyError(
                f"Exchange '{self.name}' references unknown models: {', '.join(missing)}"
            )
        models = tuple(built_models[name].model for name in self.models)
        return self.builder(simulation, models, **self.options)

    def with_options(self, **overrides: Any) -> ExchangeSpec:
        """Return a copy with selected exchange options added or replaced.

        Parameters
        ----------
        **overrides
            Builder keyword arguments merged over ``options`` key by key. Not
            checked here.
        """

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
        """Recreate an exchange from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A payload as :meth:`to_dict` writes it. ``"name"`` and ``"builder"``
            are required (``KeyError``); ``"models"`` and ``"options"`` default
            to empty, and fewer than two models raises ``ValueError``.
        """

        return cls(
            name=data["name"],
            builder=_resolve_callable(data["builder"]),
            models=tuple(data.get("models", ())),
            options=_spec_options(data.get("options", {})),
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
    simulation : flopy.mf6.MFSimulation
        The built FloPy ``MFSimulation``.
    models : dict of str to BuiltModel
        Built models keyed by model name.
    packages : dict of str to object
        Simulation-level FloPy packages (TDIS, IMS/EMS) keyed by their
        ``PackageSpec.name``; a disabled package is absent.
    exchanges : dict of str to object
        Inter-model exchange objects keyed by ``ExchangeSpec.name``.
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
        """Return one built FloPy model by name.

        Parameters
        ----------
        name : str
            A model name, a key of :attr:`models`; an absent one raises
            ``KeyError``.
        """

        return self.built_model(name).model

    def built_model(self, name: str) -> BuiltModel:
        """Return one ``BuiltModel`` wrapper by name.

        Parameters
        ----------
        name : str
            A model name, a key of :attr:`models`; an absent one raises
            ``KeyError``.
        """

        try:
            return self.models[name]
        except KeyError as error:
            raise KeyError(f"Built simulation has no model named '{name}'.") from error

    def gwf(self, name: str) -> GwfModel:
        """Return a named model as ``flopy.mf6.ModflowGwf``.

        Parameters
        ----------
        name : str
            A model name, a key of :attr:`models`; an absent one raises
            ``KeyError``, and a model of another type raises ``TypeError``.
        """

        return self.built_model(name).as_gwf()

    def gwt(self, name: str) -> GwtModel:
        """Return a named model as ``flopy.mf6.ModflowGwt``.

        Parameters
        ----------
        name : str
            A model name, a key of :attr:`models`; an absent one raises
            ``KeyError``, and a model of another type raises ``TypeError``.
        """

        return self.built_model(name).as_gwt()

    def gwe(self, name: str) -> GweModel:
        """Return a named model as ``flopy.mf6.ModflowGwe``.

        Parameters
        ----------
        name : str
            A model name, a key of :attr:`models`; an absent one raises
            ``KeyError``, and a model of another type raises ``TypeError``.
        """

        return self.built_model(name).as_gwe()

    def prt(self, name: str) -> PrtModel:
        """Return a named model as ``flopy.mf6.ModflowPrt``.

        Parameters
        ----------
        name : str
            A model name, a key of :attr:`models`; an absent one raises
            ``KeyError``, and a model of another type raises ``TypeError``.
        """

        return self.built_model(name).as_prt()

    def package(self, name: str) -> Any:
        """Return one built simulation-level package by name.

        Parameters
        ----------
        name : str
            The package's ``PackageSpec.name``, a key of :attr:`packages`,
            matched exactly; an absent one raises ``KeyError``.
        """

        try:
            return self.packages[name]
        except KeyError as error:
            raise KeyError(
                f"Built simulation has no simulation-level package named '{name}'."
            ) from error


def _size_ats(simulation, records) -> None:
    """Set ``MAXATS`` on the ATS package FloPy creates from ``ats_perioddata``.

    FloPy does not size the ``MAXATS`` dimension from the record list, so the
    written file declares ``MAXATS 1`` however many periods were supplied and MF6
    reads only the first. Measured on a 4-period stack: 4 records in, ``MAXATS 1``
    out. The legacy ``TemporalDiscretization`` path already corrected this; the
    spec path did not, so ``ats_perioddata=`` was quietly broken here.
    """

    if not records:
        return
    from myflopy.modflow.mf6.simulation.discretization import _set_maxats

    _set_maxats(simulation, len(records))


def _bicgstab_for_implicit_lakes(solvers, simulation) -> list:
    """BICGSTAB for a solver ``mf.simulation`` generated, once its model shows an implicit lake.

    ``mf.simulation`` picks BICGSTAB for a model whose lake is declared ``implicit=True``
    inline, but it cannot see through an ``mf.ref`` -- the library is not reachable
    there -- so a library lake got the SIMPLE solver, CG by default, and the build then
    refused a solver the caller never declared (ledger 326; measured with
    ``project.add_package('lak/pit', ...)`` + ``mf.ref('lak/pit')``). Here the models
    are built and the refs resolved, so a solver carrying ``GENERATED_SOLVER`` that
    names no ``linear_acceleration`` is given BICGSTAB when a model it solves holds a
    lake writing ``IMPLICIT``. A solver the caller declared is left alone, so a CG one
    is still refused by name.
    """

    implicit = set(implicit_lake_models(simulation))
    if not implicit:
        return list(solvers)
    chosen = []
    for package in solvers:
        options = package.options or {}
        solved = implicit & {str(name) for name in options.get("models") or ()}
        if (
            (package.metadata or {}).get(GENERATED_SOLVER)
            and "linear_acceleration" not in options
            and solved
        ):
            logger.info(
                "solver %r: BICGSTAB, because model(s) %s hold an implicit lake (ledger 326)",
                package.name, ", ".join(sorted(solved)),
            )
            package = replace(package, options={**options, "linear_acceleration": "BICGSTAB"})
        chosen.append(package)
    return chosen


def executable_name(value: Any, *, where: str) -> str:
    """The MF6 executable as a ``str``: ONE rule wherever one is stored (ledger 447).

    A ``str`` is kept as its text and an ``os.PathLike`` is read through
    ``os.fspath``, so ``Path(shutil.which("mf6"))`` is stored as its text. A run
    manifest and :meth:`SimulationSpec.to_dict` write JSON, which refuses a
    ``Path`` -- stored as-is it failed only later, writing the manifest. Anything
    else (``None`` included, and a ``bytes`` path) is refused, naming ``where``.
    The result is always a plain ``str``: a ``numpy.str_`` (a ``str`` subclass, an
    element of a numpy string array) passes ``isinstance``, and PyYAML, which
    dispatches on the exact type, refused it in :meth:`SimulationSpec.to_yaml`.

    Parameters
    ----------
    value : str or os.PathLike
        The executable as the caller gave it. A ``str`` (or ``str`` subclass such
        as ``numpy.str_``) is returned as a plain ``str`` of the same text; an
        ``os.PathLike`` whose ``os.fspath`` is a ``str`` returns that text.
        Anything else raises ``TypeError``.
    where : str
        What is storing it (``"Run 'baseline'"``), named in the refusal.

    Returns
    -------
    str
        The executable's name or path as text.
    """

    if isinstance(value, str):
        return str(value)
    if isinstance(value, os.PathLike):
        name = os.fspath(value)
        if isinstance(name, str):
            return str(name)
    raise TypeError(
        f"executable= ({where}) takes the MF6 executable's name or path as a str "
        f"or os.PathLike (e.g. 'mf6' or Path('/opt/mf6/bin/mf6')), not {value!r}."
    )


@dataclass(frozen=True, slots=True)
class SimulationSpec:
    """One complete MF6 simulation: timing, solver(s), and the model(s).

    The top-level declarative object. It holds the simulation-wide packages
    (``mf.tdis`` and one or more ``mf.ims`` solvers) plus the ``ModelSpec`` models
    (built with ``mf.gwf``/``mf.gwt``/``mf.gwe``/``mf.prt``) and any inter-model
    ``exchanges``. ``build_flopy(workspace)`` materializes it into a FloPy
    simulation; or register it on a :class:`~myflopy.workspace.Project` with
    ``project.add_simulation(sim)`` and run it via ``project.prepare_run(...)``.

    ``mf.simulation(*models, packages=[...])`` builds one and fills in a default
    TDIS and solver; construct ``SimulationSpec`` directly only to control every
    field.

    Parameters
    ----------
    name : str
        Simulation name: passed to the builder as ``sim_name=``, and the default
        run name and run directory (relative to the current working directory)
        when neither ``run_name`` nor a workspace is given.
    models : sequence of ModelSpec
        The model specs (``mf.gwf(...)`` etc.), built in this order. Names must
        be unique (``ValueError``). Stored as a tuple.
    packages : sequence of PackageSpec, PackageRef or str, optional
        Simulation-wide packages -- ``mf.tdis(...)`` and ``mf.ims(...)``/
        ``mf.ems(...)`` solver(s); a ``str`` becomes a :class:`PackageRef`. A
        MODEL package (``chd``, ``npf``, ...) raises ``ValueError``, as do
        duplicate names. At build TDIS is built first, then the models and
        exchanges, then every other package. Stored as a tuple. Default: empty.
    exchanges : sequence of ExchangeSpec, optional
        Inter-model exchanges (e.g. ``mf.exchange(flow, transport)``). Names
        must be unique. Nothing is wired automatically: coupled models without
        one build and run uncoupled. Default: empty.
    options : dict, optional
        Keyword arguments for the simulation constructor, called as
        ``builder(sim_name=name, **options)``. ``sim_ws`` is overwritten by the
        workspace whenever one is known. Default: empty.
    builder : callable, optional
        The simulation constructor. Default ``flopy.mf6.MFSimulation``; must be
        importable for :meth:`to_dict`.
    workspace : Path or str, optional
        The default run directory for :meth:`build` and :meth:`build_flopy`
        when they are given none (passed to FloPy as ``sim_ws``). Default
        ``None``.
    executable : str or os.PathLike, optional
        The MF6 executable a :meth:`build` run is given (set as the FloPy
        simulation's ``exe_name`` there). :meth:`build_flopy` does not use it.
        Default ``"mf6"``.
        A path is stored as its text (``os.fspath``); anything but a ``str`` or
        ``os.PathLike`` raises ``TypeError``.
    run_name : str, optional
        The run name :meth:`build` uses when not given one, and the default run
        directory name; :meth:`derive` sets it. ``None`` falls back to
        ``name``. Default ``None``.
    derived_from : str, optional
        Provenance: the name of the simulation this one was derived from, set
        by :meth:`derive` (and ``with_packages(name=)``). Serialized; nothing
        builds from it. Default ``None``.
    lineage : sequence of dict, optional
        Provenance: one record per recorded edit (``{"operation": "derive" |
        "with_packages" | "with_grid", ...}``), appended by the verbs that
        record themselves. Serialized; nothing builds from it. Stored as a
        tuple. Default: empty.

    Examples
    --------
    >>> sim = mf.SimulationSpec("baseline", models=(flow,),
    ...     packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
    ...               mf.ims(models=("valley",), complexity="MODERATE")])
    >>> built = sim.build_flopy("runs/baseline")     # -> BuiltSimulation (FloPy objects)
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
            "executable",
            executable_name(self.executable, where=f"SimulationSpec {self.name!r}"),
        )
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
        _refuse_model_packages(concrete_packages, self.name)
        if len(concrete_packages) == len(self.packages):
            _validate_concrete_packages(concrete_packages)

    def with_model(self, model: ModelSpec) -> SimulationSpec:
        """Return a copy with ``model`` added or replaced by name.

        Parameters
        ----------
        model : ModelSpec
            Replaces the model of the same ``name`` in place (keeping its
            position), or is appended when no model has that name.
        """

        models = list(self.models)
        for index, existing in enumerate(models):
            if existing.name == model.name:
                models[index] = model
                break
        else:
            models.append(model)
        return replace(self, models=tuple(models))

    def with_models(self, *models: ModelSpec) -> SimulationSpec:
        """Return a copy with each supplied model added or replaced by name.

        Parameters
        ----------
        *models : ModelSpec
            Each one applied in order through :meth:`with_model`. Not flattened:
            pass models as separate arguments (``*my_list``), not a list.
        """

        spec = self
        for model in models:
            spec = spec.with_model(model)
        return spec

    def model(self, name: str) -> ModelSpec:
        """Return one model specification by name.

        Parameters
        ----------
        name : str
            The ``ModelSpec.name``, matched exactly; an absent one raises
            ``KeyError``.
        """

        for model in self.models:
            if model.name == name:
                return model
        raise KeyError(f"Simulation '{self.name}' has no model named '{name}'.")

    def with_packages(
        self,
        *packages: PackageEntry | str | Iterable[Any],
        on: str | None = None,
        slot: str | None = None,
        name: str | None = None,
    ) -> SimulationSpec:
        """Return a copy with packages added or replaced. One verb, both tiers.

        Without ``on=`` these are SIMULATION packages -- timing and solvers. A
        model package raises; it belongs on a model.

        ``on="flow"`` targets that model's packages instead, which is what
        ``add_package`` and ``replace_package`` used to do between them. They
        differed only in the absent-slot case, and their shared assertion now
        rides on ``slot=``::

            sim.with_packages(mf.tdis(nper=12, perioddata=spd))
            sim.with_packages(mf.ref("npf/high"), on="flow", slot="npf/base")
            sim.with_packages(chd, drn, on="flow")

        ``name=`` returns the copy under a new simulation name, routed through
        :meth:`derive` so ``derived_from`` and ``lineage`` are recorded -- one
        step from a simulation to a named variant of it.

        Parameters
        ----------
        *packages : PackageSpec, PackageRef, str, or iterables of them
            The entries to add or replace; lists are flattened one level and a
            ``str`` is a :class:`PackageRef` key. With ``on=`` they go to
            :meth:`ModelSpec.with_packages` on that model. Without it each
            replaces the simulation package of the same full ``name`` (a ref's
            key) or is appended, and a model package raises ``ValueError``.
        on : str, optional
            A model name: edit that model's packages instead of the
            simulation's, and record a ``with_packages`` lineage entry. An
            unknown model raises ``KeyError``. Default ``None``.
        slot : str, optional
            With ``on=`` only (``ValueError`` without it): the existing model
            package to replace, as :meth:`ModelSpec.with_packages` takes it.
            Default ``None``.
        name : str, optional
            Return the copy through :meth:`derive` under this simulation name
            (also the new ``run_name``). Default ``None`` keeps the name.
        """

        entries = [_package_entry(entry) for entry in _flatten_packages(packages)]
        if on is not None:
            model = self.model(on)
            spec = self.with_model(model.with_packages(*entries, slot=slot))
            labels = [_package_entry_label(entry) for entry in entries]
            spec = replace(
                spec,
                lineage=(
                    *spec.lineage,
                    {
                        "operation": "with_packages",
                        "model": on,
                        "packages": labels,
                    },
                ),
            )
            return spec.derive(name) if name is not None else spec

        if slot is not None:
            raise ValueError(
                "slot= targets a package inside a MODEL; pass on=<model name> "
                "with it. A simulation's own packages are matched by full name."
            )
        current = list(self.packages)
        for entry in entries:
            for index, existing in enumerate(current):
                if existing.name == entry.name:
                    current[index] = entry
                    break
            else:
                current.append(entry)
        spec = replace(self, packages=tuple(current))
        return spec.derive(name) if name is not None else spec

    def without_packages(
        self, *names: str | Iterable[str], on: str | None = None
    ) -> SimulationSpec:
        """Return a copy without the named packages, at either tier.

        ``on="flow"`` removes from that model. Without it, removes the
        simulation's own -- which had NO verb at all before, so dropping a tdis
        or a solver meant rebuilding the ``packages`` tuple by hand.

        Raises when a name matches nothing, at both tiers.

        Parameters
        ----------
        *names : str or iterables of str
            The packages to remove; lists are flattened one level. With ``on=``
            each is matched as :meth:`ModelSpec.without_packages` matches it;
            without, against the simulation packages' full names. No match
            raises ``KeyError``.
        on : str, optional
            A model name whose packages to remove from; an unknown model raises
            ``KeyError``. Default ``None``: the simulation's own packages.
        """

        wanted = list(_flatten_packages(names))
        if on is not None:
            model = self.model(on)
            return self.with_model(model.without_packages(*wanted))

        current = list(self.packages)
        for name in wanted:
            for index, existing in enumerate(current):
                if existing.name == name:
                    del current[index]
                    break
            else:
                carried = ", ".join(entry.name for entry in current) or "nothing"
                raise KeyError(
                    f"Simulation '{self.name}' has no package '{name}' to remove. "
                    f"It carries: {carried}. Pass on=<model name> for a model's "
                    "packages."
                )
        return replace(self, packages=tuple(current))

    def with_exchange(self, exchange: ExchangeSpec) -> SimulationSpec:
        """Return a copy with ``exchange`` added or replaced by name.

        Parameters
        ----------
        exchange : ExchangeSpec
            Replaces the exchange of the same ``name`` in place, or is appended.
            Its model names are not checked until build.
        """

        exchanges = list(self.exchanges)
        for index, existing in enumerate(exchanges):
            if existing.name == exchange.name:
                exchanges[index] = exchange
                break
        else:
            exchanges.append(exchange)
        return replace(self, exchanges=tuple(exchanges))

    def with_workspace(self, workspace: Path | str) -> SimulationSpec:
        """Return a copy with a default run workspace.

        Parameters
        ----------
        workspace : Path or str
            Stored as ``workspace``: the directory :meth:`build` and
            :meth:`build_flopy` use when called without one. Not created or
            checked here.
        """

        return replace(self, workspace=workspace)

    def derive(self, name: str) -> SimulationSpec:
        """Return a named version of this fully resolved simulation recipe.

        Parameters
        ----------
        name : str
            The variant's simulation name, also set as its ``run_name``; the
            current name becomes ``derived_from`` and a ``derive`` record is
            appended to ``lineage``. Not checked for clashes with anything.
        """

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

    def with_grid(
        self,
        grid: GridSpec | GridRef | Any,
        *,
        on: str | None = None,
        name: str | None = None,
        grid_free: tuple[str, ...] | list[str] = (),
        remap: bool = False,
        rule: str | None = None,
        uncovered: str | None = None,
        library: Any = None,
    ) -> SimulationSpec:
        """Return a copy moved to ``grid``, with everything derived from the old one gone.

        ``grid`` may be a :class:`GridSpec`, a :func:`grid_ref` reference into a
        project grid library, or an already-built grid object.

        ``on`` names one model; the default moves **every model that has a grid**,
        because coupled models must share a mesh -- MF6 requires it for
        ``gwfgwt``/``gwfgwe``/``gwfprt`` and cannot detect a same-size violation
        (ledger 268), so moving one and stranding its sibling is the failure that check
        exists to catch. ``name`` derives a named variant in the same step, as
        ``with_packages(name=)`` does.

        ``grid_free`` names packages the caller asserts are mesh-independent, so they
        are not stamped. Everything else that names cells without a pin IS stamped and
        the build refuses it -- see :meth:`ModelSpec.with_grid`.

        A model built from a project names its packages by ``mf.ref`` and often its grid
        by ``mf.grid_ref``; classifying those needs the library they name, so pass
        ``library=project`` -- a moved model holding either raises naming it otherwise.

        **This reports and stamps; the BUILD refuses.** Read :meth:`move_report` to see
        what will happen before it does.

        Parameters
        ----------
        grid : GridSpec, GridRef or grid object
            The new mesh, handed to :meth:`ModelSpec.with_grid` for each moved
            model (a bare built grid is wrapped with :meth:`GridSpec.from_object`
            there). ``remap=True`` needs a built grid or a :class:`GridSpec`.
        on : str, optional
            One model to move; an unknown name raises ``KeyError``. Default
            ``None``: every model that has a ``grid`` or a ``context.grid`` --
            a model with neither is left alone, so a simulation of gridless
            models moves nothing (and says so only in ``lineage``).
        name : str, optional
            Derive a named variant in the same step, through :meth:`derive` as
            ``with_packages(name=)`` does: ``name`` becomes the simulation name
            AND its ``run_name``, the current name its ``derived_from``, and a
            ``derive`` record follows the ``with_grid`` one in ``lineage``.
            Default ``None`` keeps the name.
        grid_free : sequence of str, optional
            Package names (a ref's library key) asserted mesh-independent in
            every moved model; recorded in the ``with_grid`` lineage entry.
            Default: empty.
        remap : bool, optional
            Carry per-cell arrays onto the new mesh by area overlay, per
            :meth:`ModelSpec.with_grid`. Default ``False``.
        rule : str or None, optional
            The averaging rule for the conductivities (``k``/``k22``/``k33``)
            under ``remap=True``: ``"geometric"``, ``"mean"``, ``"harmonic"`` or
            ``"majority"``; ``ss``/``sy`` always take the area-weighted mean, and
            ``strt`` -- a head field -- is carried linear between cell centres, as
            ``model.hds.remap`` carries one (ledger 600). Any other value raises
            ``ValueError`` at the call, even when no
            model moves; given without ``remap=True`` it raises ``ValueError``
            naming ``remap=True``. Default ``None``: ``"geometric"`` under
            ``remap=True``.
        uncovered : str or None, optional
            For ``remap=True``: ``"nearest"`` or ``"raise"`` for a new cell no old
            cell reaches, per :meth:`ModelSpec.with_grid`. Checked and refused
            like ``rule``. Default ``None``: ``"nearest"`` under ``remap=True``.
        library : Project, optional
            The project whose libraries the moved models' ``mf.ref`` and
            ``mf.grid_ref`` entries name, handed to :meth:`ModelSpec.with_grid`
            for each: refs are classified as the packages they name, a
            ``grid_ref`` grid is resolved to measure against. Required when a
            moved model holds a ref no earlier move stamped, or a ``grid_ref``
            grid (``ValueError`` naming it otherwise). Nothing in the library is
            changed. Default ``None``.

        Returns
        -------
        SimulationSpec
            The moved copy, with a ``with_grid`` record appended to ``lineage``
            (and a ``derive`` record after it when ``name`` is given).
        """

        if isinstance(grid, GridRef):
            grid_label = grid.key
        elif isinstance(grid, GridSpec):
            grid_label = grid.name
        else:
            grid_label = "object"

        rule, uncovered = _remap_choices(remap, rule, uncovered)
        free = tuple(grid_free)
        if on is not None:
            moved = self.with_model(
                self.model(on).with_grid(
                    grid, grid_free=free, remap=remap, rule=rule, uncovered=uncovered,
                    library=library,
                )
            )
            targets = [on]
        else:
            targets = [
                model.name
                for model in self.models
                if model.grid is not None or getattr(model.context, "grid", None) is not None
            ]
            moved = self
            for model_name in targets:
                moved = moved.with_model(
                    moved.model(model_name).with_grid(
                        grid, grid_free=free, remap=remap, rule=rule,
                        uncovered=uncovered, library=library,
                    )
                )
        moved = replace(
            moved,
            lineage=(
                *self.lineage,
                {
                    "operation": "with_grid",
                    "models": targets,
                    "grid": grid_label,
                    "grid_free": list(free),
                },
            ),
        )
        # `derive`, as `with_packages(name=)`: a bare rename kept the SOURCE's
        # `run_name`, and `build()` defaults its workspace to `run_name`, so the
        # variant wrote over the source variant's run directory.
        return moved.derive(name) if name is not None else moved

    def move_report(
        self,
        *,
        on: str | None = None,
        grid_free: tuple[str, ...] | list[str] = (),
        library: Any = None,
    ) -> list[dict[str, str]]:
        """One row per package per model: what a grid move does to it.

        Readable on the SPEC, before anything is built, because the build is what
        refuses and `prepare_run` is a late place to find out.

        Parameters
        ----------
        on : str, optional
            Report one model only; an unknown name raises ``KeyError``. Default
            ``None``: every model.
        grid_free : sequence of str, optional
            Package names to report as mesh-independent, as :meth:`with_grid`
            would treat them. Default: empty.
        library : Project, optional
            The project whose libraries the models' ``mf.ref`` and
            ``mf.grid_ref`` entries name, as for :meth:`with_grid` -- required in
            the same cases, ``ValueError`` naming it otherwise. Default ``None``.

        Returns
        -------
        list of dict
            One ``{"model", "package", "how", "why"}`` row per package per
            model; see :meth:`ModelSpec.move_report`.
        """

        models = [self.model(on)] if on is not None else list(self.models)
        return [
            {"model": model.name, **row}
            for model in models
            for row in model.move_report(grid_free=grid_free, library=library)
        ]

    def coverage_on(self, grid: Any, *, on: str | None = None) -> dict[str, Any]:
        """Per model, how much of ``grid`` has no data under it -- in model units.

        One entry per model, because a coupled simulation's models share a mesh
        (ledger 268) but a simulation can hold models that do not. Read it before
        `with_grid`; `move_report` says what happens to each PACKAGE, this says what
        the new mesh would have nothing for.

        Parameters
        ----------
        grid : grid object or GridSpec
            The candidate new mesh, as :meth:`ModelSpec.coverage_on` takes it.
        on : str, optional
            Measure one model only; an unknown name raises ``KeyError``. Default
            ``None``: every model.

        Returns
        -------
        dict of str to CoverageReport
            Keyed by model name. A model holding no built grid (neither
            ``grid=<built grid>`` nor ``context.grid``) is left out rather than
            raising (unlike :meth:`ModelSpec.coverage_on`), so a deferred-grid
            model -- even one named by ``on`` -- yields no entry.
        """

        models = [self.model(on)] if on is not None else list(self.models)
        return {
            model.name: model.coverage_on(grid)
            for model in models
            if _mesh_in_hand(model) is not None
        }

    def _replace_grid_compat(
        self, model_name: str, grid: GridSpec | GridRef | Any
    ) -> SimulationSpec:
        """The old `replace_grid`, reached only through the warned alias below.

        It swapped the grid and re-mapped NOTHING -- every package kept cellids
        resolved on the old mesh, and nothing looked (ledger 249). `with_grid` drops
        the derived context and stamps what cannot be verified.
        """

        return self.with_grid(grid, on=model_name)

    __getattr__ = deprecated_instance_getattr(
        {
            "replace_grid": (
                "_replace_grid_compat",
                "SimulationSpec.with_grid(grid, on=model_name)",
                "0.1.0",
            )
        },
        "myflopy.specs.SimulationSpec",
    )

    def resolved_packages(
        self,
        build_context: SpecBuildContext | None = None,
    ) -> tuple[PackageSpec, ...]:
        """Return concrete simulation-level package specs.

        Parameters
        ----------
        build_context : SpecBuildContext, optional
            Supplies the package library each simulation-level
            :class:`PackageRef` is looked up in. Required only when one is
            present: without it a ref raises ``ValueError``, and a key missing
            from the library raises ``KeyError``. Default ``None``.

        Returns
        -------
        tuple of PackageSpec
            In declaration order, after the same name/``requires`` validation as
            :meth:`ModelSpec.resolved_packages`.
        """

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

        The run builds in memory: its directory is created and its manifest and
        frozen snapshot are written there, but no MF6 input file is written and
        MF6 is not run -- ``run.execute()`` does both.

        Parameters
        ----------
        workspace : Path or str, optional
            The run directory. ``None`` falls back to ``self.workspace``, then
            to ``Path(self.run_name or self.name)`` relative to the current
            working directory -- the ``run_name`` argument does not move it.
            Default ``None``.
        run_name : str, optional
            The run's name. ``None`` falls back to ``self.run_name``, then
            ``self.name``. Default ``None``.
        build_context : SpecBuildContext, optional
            Project libraries and locations for resolving refs and anchoring
            relative sources; see :meth:`build_flopy`. Default ``None``.

        Returns
        -------
        myflopy.workspace.Run
            The built run (``status`` ``"built"``), with ``run.built`` holding
            the :class:`BuiltSimulation`.
        """

        from myflopy.workspace import Run

        run_workspace = workspace if workspace is not None else self.workspace
        if run_workspace is None:
            run_workspace = Path(self.run_name or self.name)
        run = Run(
            name=run_name or self.run_name or self.name,
            workspace=Path(run_workspace),
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
        """Build and return the raw FloPy simulation objects.

        Builds TDIS, then each model (resolving deferred grids and refs), checks
        that coupled models share a mesh, builds the exchanges, then the solvers.
        No MF6 input file is written; resolving a deferred grid does write its
        working files under the grid workspace (see ``workspace``).

        Parameters
        ----------
        workspace : Path or str, optional
            The simulation directory, passed to the builder as ``sim_ws`` and
            recorded as the build context's ``simulation_workspace`` (where
            deferred grids are built, ``<workspace>/_grid/<model>``). ``None``
            falls back to ``self.workspace``; with neither, ``options`` decides
            ``sim_ws`` (FloPy's default is the current working directory) and
            the grids go under it the same way. Default ``None``.
        build_context : SpecBuildContext, optional
            Project libraries (for ``mf.ref``/``mf.grid_ref``) and the project
            root (for relative source paths). ``None`` uses an empty one, so a
            ref raises. A ``simulation_workspace`` already set on it is kept
            over ``workspace``. Default ``None``.

        Returns
        -------
        BuiltSimulation

        Raises
        ------
        ValueError
            A model that disagrees with its grid, coupled models on different
            meshes, or an implicit lake with a solver that cannot solve it.
        """

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
            _size_ats(simulation, package.options.get("ats_perioddata"))
        built_models = {
            model.name: model.build(simulation, build_context=build_context)
            for model in self.models
        }
        # Before the exchanges are built, not after: MF6 would accept this and run it.
        if mismatched := coupled_mesh_complaints(self.exchanges, built_models):
            raise ValueError(
                f"Simulation {self.name!r} couples models that are not on the same mesh:\n"
                + "\n".join(f"  - {complaint}" for complaint in mismatched)
                + "\nRe-declare the coupled models against one grid, or move them together "
                "with `with_grid(..., on=None)`."
            )
        built_exchanges = {
            exchange.name: exchange.build(simulation, built_models)
            for exchange in self.exchanges
        }
        solvers = _bicgstab_for_implicit_lakes(solvers, simulation)
        for package in solvers:
            if (built := package.build(simulation)) is not None:
                simulation_packages[package.name] = built
        # A lake solved IMPLICITLY makes the matrix asymmetric and MF6 stops on CG --
        # refused here, where the solver is known, rather than by MF6 after the files
        # are written (ledger 326).
        if refused := implicit_solver_complaints(simulation):
            raise ValueError(
                f"Simulation {self.name!r} cannot solve its implicit lakes:\n"
                + "\n".join(f"  - {complaint}" for complaint in refused)
            )
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
        """Recreate a simulation specification from :meth:`to_dict` output.

        Parameters
        ----------
        data : dict
            A payload as :meth:`to_dict` writes it. ``"name"`` is required
            (``KeyError``); models, packages, exchanges, options and lineage
            default to empty, ``"builder"`` to ``flopy.mf6.MFSimulation`` and
            ``"executable"`` to ``"mf6"``. Each model is rebuilt with
            :meth:`ModelSpec.from_dict`, so its context comes back without its
            grid, domain and surfaces (a deferred ``grid`` recipe is restored).
        """

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

    def to_yaml(self, path: str | Path | None = None) -> str:
        """Serialize this simulation to a YAML string (plan §5.6).

        Round-trips through :meth:`to_dict`, so the same rules apply: builders
        must be importable and every value JSON-representable (a computed numpy
        array baked into a package's options is not -- that raises in
        :meth:`to_dict`). If ``path`` is given the YAML is also written there.
        See :mod:`myflopy.specs_io`. TOML is not supported yet (see that module).

        Parameters
        ----------
        path : str or Path, optional
            A file to write the YAML to as well (overwritten if it exists; the
            parent directory must exist). ``None`` only returns the text.
            Default ``None``.

        Returns
        -------
        str
            The YAML text, whether or not it was written.
        """

        text = simulation_to_yaml(self)
        if path is not None:
            Path(path).write_text(text)
        return text

    @classmethod
    def from_yaml(cls, source: str | Path) -> SimulationSpec:
        """Load a simulation from YAML text or a ``.yaml`` file path (plan §5.6).

        A :class:`~pathlib.Path` (or a ``str`` naming an existing file) is read
        from disk; any other ``str`` is parsed as YAML text. The inverse of
        :meth:`to_yaml`.

        Parameters
        ----------
        source : str or Path
            A ``Path`` is always read as a file. A ``str`` is read as a file
            when it names an existing one, otherwise parsed as YAML text -- so a
            mistyped file name is parsed as YAML and fails there. The document
            must be a mapping (``ValueError`` otherwise) in the
            :meth:`to_dict` shape.
        """

        return simulation_from_yaml(source)

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
                f"{model.type_enum.value}; packages={len(model.packages)}",
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
