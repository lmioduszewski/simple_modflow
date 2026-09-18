"""Project and run lifecycle objects for spec-driven MODFLOW workflows."""

from __future__ import annotations

import json
import pickle
import re
import warnings
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from html import escape
from importlib import import_module
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import flopy

from myflopy._logging import get_logger
from myflopy.library import (
    GridLibrary,
    PackageLibrary,
    SimulationLibrary,
    SpecLibrary,
)
from myflopy.run_record import ListingOutcome, pest_tree_root, read_listing_outcome

if TYPE_CHECKING:
    import pandas as pd
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.project.bundle import bundle_grid_sources
from myflopy.project.freeze import (
    discard_snapshot,
    freeze_simulation,
    frozen_summary,
    thaw_simulation,
)
from myflopy.project.run_model import load_mf6_run, patch_simulation_plot
from myflopy.specs import (
    SOURCE_SHA,
    BuiltSimulation,
    GridRef,
    GridSpec,
    ModelSpec,
    PackageRef,
    PackageSpec,
    SimulationSpec,
    SpecBuildContext,
    _grid_entry,
    _json_value,
    _package_entry_label,
    dropped_context_fields,
    package_digest,
    package_move,
    source_fingerprint,
)

logger = get_logger(__name__)

RUN_MANIFEST_NAME = "run.json"
# A run's own copy of the library packages it references, beside its `_grid/`.
RUN_LIBRARY_DIR = "_library"
PROJECT_MANIFEST_NAME = "project.json"
PROJECT_SPEC_NAME = "project_spec.json"


class ModelView(SimulationBase):
    """Preferred myflopy view over one live built GWF model.

    The view keeps the raw FloPy objects available as ``sim`` and ``gwf`` while
    exposing the familiar myflopy helpers such as ``cor()``, ``hds``,
    ``packages``, ``outputs``, and budget accessors.
    """

    def __init__(self, *, run: Run, model_name: str):
        """Build a model view over one model of a completed ``run``."""

        self._initialize_from_built_run(run, model_name)


def load_run(
    workspace: Path | str,
    *,
    executable: str = "mf6",
    load: bool = True,
) -> Run:
    """Reopen a previously prepared/executed :class:`Run` from its workspace.

    Reads the run's JSON manifest and (with ``load=True``) the MF6 input files,
    returning a :class:`Run` you can inspect, re-execute, or read results from --
    e.g. ``run.model("flow")`` for the model view, or ``run.pest_runs`` for any
    calibrations done on it. Use this to come back to a run in a later session.

    Parameters
    ----------
    workspace
        The run directory (``<project>/runs/<name>``).
    executable
        MF6 executable name/path (default ``"mf6"``).
    load
        Eagerly load the FloPy simulation (default ``True``); ``False`` defers it.

    Examples
    --------
    >>> run = mf.load_run("project/runs/baseline")
    >>> heads = run.model("valley").hds.array()
    """

    return Run.load(workspace, executable=executable, load=load)


def _refuse_overwrite(
    library: Any, key: str, verb: str, overwrite: bool, project_name: str
) -> None:
    """Refuse to replace an existing library entry unless asked to.

    ``add_package`` used to replace silently, while ``packages.fork`` refused the
    same collision -- so whether re-using a key destroyed the old entry or raised
    depended on which verb you happened to reach for. The registration verbs were
    the unguarded ones, which is backwards: a fork is a deliberate new thing,
    while re-registering a key is usually a second script run or a copied line.

    ``overwrite=True`` is the way to mean it. ``packages.edit(key, ...)`` is the
    way to change what a key MEANS, and is deliberately never refused.
    """

    if overwrite or key not in library:
        return
    kind = getattr(library, "_kind", "entry")
    raise KeyError(
        f"Project '{project_name}' already has a {kind} '{key}'. Pass "
        f"overwrite=True to replace it, or pick another key. To change what "
        f"'{key}' MEANS, use {kind}s.edit(...) -- that records the edit and "
        f"project.runs.drift() will show which runs it invalidated."
    )


def _refuse_wrong_library(value: Any, expected: type, verb: str) -> None:
    """Refuse a spec handed to the wrong project library.

    ``project.add_package("flow_hi", mf.gwf(...))`` was accepted in full:
    ``validate()`` returned clean, ``save()`` reported success while silently
    writing no simulations, and the NEXT session's ``Project.load`` raised
    ``KeyError: 'builder'`` -- the whole project unopenable, from one wrong verb
    and with every signal in between saying it had worked.

    Only a spec that clearly belongs to a different library is refused. A built
    mesh handed to ``add_grid`` is none of these types and passes, which is the
    point: ``add_grid`` legitimately takes a ``GridSpec``, a ``GridRef``, or an
    already-built grid object.
    """

    elsewhere = {
        ModelSpec: (
            "a model",
            "wrap it with mf.simulation(model) and pass that to "
            "project.add_simulation()",
        ),
        SimulationSpec: ("a simulation", "use project.add_simulation(sim)"),
        PackageSpec: ("a package", "use project.add_package(key, spec)"),
        GridSpec: ("a grid", "use project.add_grid(key, grid)"),
    }
    for kind, (label, fix) in elsewhere.items():
        if kind is not expected and isinstance(value, kind):
            raise TypeError(
                f"project.{verb}() was given {label} ({type(value).__name__}) -- {fix}."
            )


def _library_grid_spec(key: str, value) -> GridSpec:
    """Normalize a Project grid-library entry, requiring a resolvable GridSpec."""

    entry = _grid_entry(value)
    if not isinstance(entry, GridSpec):
        raise TypeError(
            f"Project grid library entry {key!r} must be a GridSpec, "
            f"got {type(entry).__name__}"
        )
    return entry


def _utc_now() -> str:
    """Return a current UTC timestamp suitable for a manifest."""

    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write a readable JSON manifest, creating its parent directory."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _read_json(path: Path) -> dict[str, Any]:
    """Read a JSON document from ``path``."""

    return json.loads(path.read_text(encoding="utf-8"))


def _slug(value: str) -> str:
    """Return a filesystem-safe slug for a project, simulation, or key part."""

    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._-")
    if not slug:
        raise ValueError(
            "Project, simulation, and package names must contain at least one "
            "filesystem-safe character."
        )
    return slug


def _package_key_parts(key: str) -> tuple[str, ...]:
    """Split a package key such as ``"npf/base"`` into safe path segments."""

    parts = tuple(_slug(part) for part in str(key).split("/") if part)
    if not parts:
        raise ValueError("Package spec keys must contain at least one safe segment.")
    return parts


def _colliding_keys(
    keys: Iterable[str], path_for: Callable[[str], Path]
) -> list[tuple[str, str, Path]]:
    """Find library keys that map to the same file, in key order.

    ``_slug`` is lossy -- it rewrites every unsafe character to ``_`` -- so
    ``"npf/high k"`` and ``"npf/high_k"`` are distinct keys that name one file.
    Left undetected, ``save()`` writes both and the second silently wins; the
    reload then hands the survivor's spec back under BOTH keys. Grids collide
    more readily than packages: their whole key is slugged flat, so
    ``"base/coarse"`` and ``"base_coarse"`` collide too.
    """

    seen: dict[Path, str] = {}
    collisions: list[tuple[str, str, Path]] = []
    for key in sorted(keys):
        path = path_for(key)
        first = seen.setdefault(path, key)
        if first != key:
            collisions.append((first, key, path))
    return collisions


def _artifact_versions() -> dict[str, str]:
    """Return library versions recorded beside a pickled grid or package.

    Built grids and array-bearing packages are pickled as regenerable caches,
    not an archival format. These versions let :meth:`Project.load` warn when a
    pickle was written under a different stack that may no longer unpickle
    cleanly.
    """

    versions: dict[str, str] = {}
    for name in ("flopy", "geopandas", "shapely", "numpy", "pandas"):
        try:
            module = import_module(name)
        except ImportError:
            # An optional stack member is not installed; it simply does not
            # appear in the recorded versions.
            continue
        versions[name] = getattr(module, "__version__", "unknown")
    return versions


def _materialize_models(spec: SimulationSpec) -> SimulationSpec:
    """Give each model its own subdirectory when a simulation has multiple models.

    A coupled simulation (for example GWF with GWT, GWE, or PRT) writes each
    model's input files into ``<workspace>/<model name>/`` so it is clear which
    files belong to which model. Single-model simulations stay flat, and a model
    that already declares a ``model_rel_path`` is left untouched.
    """

    if len(spec.models) <= 1:
        return spec
    models = []
    changed = False
    for model in spec.models:
        if model.options.get("model_rel_path") in (None, ".", ""):
            models.append(model.with_options(model_rel_path=model.name))
            changed = True
        else:
            models.append(model)
    return spec.with_models(*models) if changed else spec


@dataclass(frozen=True, slots=True)
class ProjectLayout:
    """Directory layout for a durable project on disk.

    This is internal infrastructure, but it can be supplied to :class:`Project`
    to customize where specs, package definitions, and inputs are written.
    """

    root: Path
    specs_dir_name: str = "specs"
    simulations_dir_name: str = "simulations"
    inputs_dir_name: str = "inputs"
    packages_dir_name: str = "packages"
    grids_dir_name: str = "grids"

    def __post_init__(self) -> None:
        """Coerce ``root`` to a :class:`~pathlib.Path` (on this frozen dataclass)."""

        object.__setattr__(self, "root", Path(self.root))

    @classmethod
    def for_project(cls, name: str, root: Path | str | None = None) -> ProjectLayout:
        """Return the default layout for a project name (defaults under ~/mf6)."""

        project_root = Path.home() / "mf6" / _slug(name) if root is None else Path(root)
        return cls(project_root)

    @property
    def manifest_path(self) -> Path:
        """Path to the project manifest file under the project root."""

        return self.root / PROJECT_MANIFEST_NAME

    @property
    def specs_dir(self) -> Path:
        """Directory holding all durable spec files."""

        return self.root / self.specs_dir_name

    @property
    def project_spec_path(self) -> Path:
        """Path to the serialized project spec file."""

        return self.specs_dir / PROJECT_SPEC_NAME

    @property
    def simulation_specs_dir(self) -> Path:
        """Directory holding one JSON file per durable simulation spec."""

        return self.specs_dir / self.simulations_dir_name

    @property
    def package_specs_dir(self) -> Path:
        """Directory holding durable project package specs."""

        return self.specs_dir / self.packages_dir_name

    @property
    def inputs_dir(self) -> Path:
        """Directory holding project-owned input data files."""

        return self.root / self.inputs_dir_name

    def simulation_spec_path(self, name: str) -> Path:
        """Return the durable spec path for one simulation name."""

        return self.simulation_specs_dir / f"{_slug(name)}.json"

    def package_spec_path(self, key: str) -> Path:
        """Return the durable spec path for one project package key."""

        *parents, filename = _package_key_parts(key)
        return self.package_specs_dir.joinpath(*parents, f"{filename}.json")

    def package_pickle_path(self, key: str) -> Path:
        """Return the pickle path for one array-bearing project package key."""

        *parents, filename = _package_key_parts(key)
        return self.package_specs_dir.joinpath(*parents, f"{filename}.pkl")

    def package_sidecar_path(self, key: str) -> Path:
        """Return the version-sidecar path for one pickled project package key."""

        *parents, filename = _package_key_parts(key)
        return self.package_specs_dir.joinpath(*parents, f"{filename}.versions.json")

    @property
    def grid_specs_dir(self) -> Path:
        """Directory holding durable project grid specs."""

        return self.specs_dir / self.grids_dir_name

    def grid_spec_path(self, key: str) -> Path:
        """Return the durable spec path for one project grid key."""

        return self.grid_specs_dir / f"{_slug(key)}.json"

    def grid_pickle_path(self, key: str) -> Path:
        """Return the pickle path for one project grid key."""

        return self.grid_specs_dir / f"{_slug(key)}.pkl"

    def grid_sidecar_path(self, key: str) -> Path:
        """Return the version-sidecar path for one pickled project grid key."""

        return self.grid_specs_dir / f"{_slug(key)}.versions.json"

    def ensure(self) -> None:
        """Create the durable spec directories."""

        for path in (
            self.root,
            self.specs_dir,
            self.simulation_specs_dir,
            self.package_specs_dir,
            self.inputs_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> dict[str, Any]:
        """Serialize the root and directory names for the project manifest."""

        return {
            "root": str(self.root.as_posix()),
            "specs_dir": self.specs_dir_name,
            "simulations_dir": self.simulations_dir_name,
            "inputs_dir": self.inputs_dir_name,
            "packages_dir": self.packages_dir_name,
            "grids_dir": self.grids_dir_name,
        }

    @classmethod
    def from_dict(cls, root: Path | str, data: dict[str, Any]) -> ProjectLayout:
        """Rebuild a layout from saved directory names, rooted at ``root``.

        The current ``root`` is used (the project may have moved); only the
        directory names are restored.
        """

        return cls(
            Path(root),
            specs_dir_name=data.get("specs_dir", "specs"),
            simulations_dir_name=data.get("simulations_dir", "simulations"),
            inputs_dir_name=data.get("inputs_dir", "inputs"),
            packages_dir_name=data.get("packages_dir", "packages"),
            grids_dir_name=data.get("grids_dir", "grids"),
        )


def _spec_summary(spec: SimulationSpec) -> dict[str, Any]:
    """Return stable, readable provenance for a simulation specification."""

    def package_metadata(package: PackageSpec | PackageRef) -> dict[str, Any]:
        """The provenance metadata of a concrete package (empty for a bare ref)."""

        # Through the converter: a run manifest is written with `json.dumps` and no
        # `default=`, and neither metadata dict goes through `PackageSpec.to_dict`
        # on this path. `spec.with_metadata(calibrated=np.True_)` or a context
        # metadata value read from a table serialized fine through `to_dict` and
        # then made `prepare_run` raise `TypeError: Object of type bool is not JSON
        # serializable` -- numpy 2's name for its own type -- so the run could not
        # be prepared at all (ledger 265).
        return _json_value(package.metadata) if isinstance(package, PackageSpec) else {}

    return {
        "name": spec.name,
        "models": [
            {
                "name": model.name,
                "type": model.type_enum.value,
                # Use the entry label so a reference keeps its library key
                # (e.g. "npf/high_k"), preserving which variant was used.
                "packages": [
                    _package_entry_label(package)
                    for package in model.packages
                    if package.enabled
                ],
                "package_metadata": {
                    _package_entry_label(package): metadata
                    for package in model.packages
                    if package.enabled and (metadata := package_metadata(package))
                },
                "context_metadata": _json_value(model.context.metadata),
                "post_build_hooks": [hook.name for hook in model.hooks],
            }
            for model in spec.models
        ],
        "simulation_packages": [
            _package_entry_label(package) for package in spec.packages if package.enabled
        ],
        "exchanges": [
            {"name": exchange.name, "models": list(exchange.models)}
            for exchange in spec.exchanges
        ],
    }


def _package_ref_keys(spec: SimulationSpec) -> list[str]:
    """The library keys ``spec`` references by ``mf.ref``, each once, in order."""

    entries = [
        *spec.packages,
        *(entry for model in spec.models for entry in model.packages),
    ]
    return list(
        dict.fromkeys(entry.key for entry in entries if isinstance(entry, PackageRef))
    )


def _persist_package(
    package: PackageSpec,
    *,
    spec_path: Path,
    pickle_path: Path,
    sidecar_path: Path,
    root: Path,
) -> bool:
    """Write one package as readable JSON, pickling it when it carries arrays.

    Shared by the project library and a run's own copy of it, so the two hold
    packages in one format. The pickle is named in the JSON document relative to
    ``root``. Returns whether it fell back to pickling.
    """

    try:
        doc = package.to_dict()
    except ValueError:
        doc = None
    if doc is not None:
        try:
            _write_json(spec_path, doc)
            return False
        except TypeError:
            # The fallback below exists for a package that cannot be written as
            # JSON, and it could not be reached for the commonest way that happens:
            # `to_dict` returns, and `json.dumps` raises `TypeError` a line later,
            # outside the guard. So the exception escaped `Project.save()` and broke
            # its documented "persist what it can rather than refusing the whole
            # project" contract, leaving a half-written project directory behind
            # (ledger 265). `_json_value` closes the numpy routes into this; the
            # fallback is for whatever is left.
            logger.debug(
                "Package '%s' has a field JSON cannot write; pickling it instead",
                package.name, exc_info=True,
            )
            spec_path.unlink(missing_ok=True)
    pickle_path.parent.mkdir(parents=True, exist_ok=True)
    with pickle_path.open("wb") as handle:
        pickle.dump(package, handle)
    _write_json(sidecar_path, _artifact_versions())
    _write_json(
        spec_path,
        {
            "kind": "PackageSpec",
            "name": package.name,
            "method": "pickle",
            "pickle_path": pickle_path.relative_to(root).as_posix(),
        },
    )
    return True


def _warn_artifact_versions(
    kind: str, key: str, sidecar_path: Path, *, stacklevel: int
) -> None:
    """Warn if a pickled artifact was written under a different library stack."""

    if not sidecar_path.exists():
        return
    stored = _read_json(sidecar_path)
    current = _artifact_versions()
    mismatched = {
        name: (stored.get(name), current.get(name))
        for name in current
        if name in stored and stored[name] != current[name]
    }
    if mismatched:
        warnings.warn(
            f"Pickled {kind} '{key}' was written under different library "
            f"versions {mismatched}. If it fails to load, rebuild it.",
            stacklevel=stacklevel,
        )


def _revive_package(
    key: str,
    doc: dict[str, Any],
    *,
    root: Path,
    sidecar_path: Path,
    stacklevel: int,
) -> PackageSpec:
    """Reconstruct one package from the document :func:`_persist_package` wrote."""

    if doc.get("method") == "pickle":
        _warn_artifact_versions("package", key, sidecar_path, stacklevel=stacklevel)
        with (root / doc["pickle_path"]).open("rb") as handle:
            return pickle.load(handle)
    return PackageSpec.from_dict(doc)


def _write_package_snapshot(
    workspace: Path,
    packages: dict[str, PackageSpec],
    *,
    run_name: str,
    existing: dict[str, str] | None = None,
) -> dict[str, str]:
    """Copy library packages into a run's own folder; return the manifest's index.

    The index maps each key to its document, relative to the run folder, so a
    reopened run finds a copy by KEY and never by guessing a file name: ``_slug``
    is lossy, and ``"npf/high k"`` and ``"npf/high_k"`` would share one. A clash
    gets ``~2``, ``~3``, which no slug can produce. ``existing`` is an index
    already on disk, whose files are left alone.

    A package that is neither JSON-encodable nor picklable is left out with a
    WARNING: preparing a run must not fail over a copy. A key with no copy
    resolves from the library when the run is rebuilt, which :meth:`Run.build`
    reports if the library has changed.
    """

    taken = {relative.removesuffix(".json") for relative in (existing or {}).values()}
    index: dict[str, str] = {}
    for key, package in packages.items():
        natural = Path(RUN_LIBRARY_DIR, "packages", *_package_key_parts(key))
        stem, clash = natural, 1
        while stem.as_posix() in taken:
            clash += 1
            stem = natural.with_name(f"{natural.name}~{clash}")
        folder, name = workspace / stem.parent, stem.name
        try:
            _persist_package(
                package,
                spec_path=folder / f"{name}.json",
                pickle_path=folder / f"{name}.pkl",
                sidecar_path=folder / f"{name}.versions.json",
                root=workspace,
            )
        except (pickle.PicklingError, TypeError, AttributeError) as error:
            (folder / f"{name}.pkl").unlink(missing_ok=True)
            logger.warning(
                "Run '%s': library package '%s' could not be copied into the run "
                "(%s). Rebuilding the run in a later session resolves it from the "
                "library as it is then.",
                run_name,
                key,
                error,
            )
            continue
        taken.add(stem.as_posix())
        index[key] = f"{stem.as_posix()}.json"
    return index


def _read_package_snapshot(
    workspace: Path, index: dict[str, str], *, run_name: str
) -> dict[str, PackageSpec]:
    """Read a run's own copies of its library packages back, by key.

    A missing or unreadable copy RAISES rather than falling back to the library:
    the library may hold a different package by now, and building that silently
    is the defect these copies exist to prevent.
    """

    packages: dict[str, PackageSpec] = {}
    for key, relative in index.items():
        path = workspace / relative
        if not path.exists():
            raise FileNotFoundError(
                f"Run '{run_name}' was prepared with its own copy of library "
                f"package '{key}', and {path} is missing. It is not rebuilt from "
                "the library instead, which may now hold a different package: "
                "restore the file, or prepare a new run."
            )
        stem = path.name.removesuffix(".json")
        try:
            packages[key] = _revive_package(
                key,
                _read_json(path),
                root=workspace,
                sidecar_path=path.with_name(f"{stem}.versions.json"),
                stacklevel=6,
            )
        except (
            ValueError,
            TypeError,
            KeyError,
            AttributeError,
            ImportError,
            EOFError,
            pickle.UnpicklingError,
        ) as error:
            raise ValueError(
                f"Run '{run_name}': its copy of library package '{key}' ({path}) "
                f"could not be read back ({error}). It is not rebuilt from the "
                "library instead, which may now hold a different package."
            ) from error
    return packages


def _spec_from_recipe(
    recipe: dict[str, Any] | None, *, name: str
) -> SimulationSpec | None:
    """Rebuild a run's :class:`SimulationSpec` from its manifest recipe.

    Returns ``None`` when the manifest has no recipe (written before recipes
    existed, or by a run whose spec embedded raw arrays) or when the recipe
    cannot be revived -- a builder whose module moved, a source kind this version
    no longer knows. A run that cannot be rebuilt is the status quo, so this
    degrades to summary-only rather than making the run unopenable.
    """

    if not recipe:
        return None
    try:
        return SimulationSpec.from_dict(recipe)
    except (ValueError, TypeError, KeyError, AttributeError, ImportError) as error:
        logger.warning(
            "Run '%s': stored spec could not be revived (%s); the manifest's recipe "
            "is unusable. The run still opens, and rebuilds from its frozen snapshot "
            "if it has one.",
            name,
            error,
        )
        return None


_MISSING_KEY = object()
"""Distinguishes "the library has no such key" from "the key digests to None"."""


def _drift_state(ran_with: str | None, now: Any) -> str:
    """Classify one recorded package against what its key holds now.

    Four states, and the two Nones are why there are four rather than three:
    a key the library no longer has is ``removed``, while a key it still has but
    cannot fingerprint (a closure builder, an unpicklable option) is
    ``unfingerprinted``. Calling either one ``unchanged`` would be the same class
    of lie as ``status="written"`` reading as success.
    """

    if now is _MISSING_KEY:
        return "removed"
    if ran_with is None or now is None:
        return "unfingerprinted"
    return "unchanged" if now == ran_with else "edited"


class RunLibrary(SpecLibrary["Run"]):
    """Every run under a project -- with a manifest, without one, or broken.

    Unlike the three spec libraries this is **discovered from disk**, not
    registered: a run is created by ``prepare_run``, never added. It is keyed the
    way each run is addressable -- by directory name for a project-owned run, and
    by path relative to the project root for a bare MF6 workspace, so
    ``runs/from_usg`` and ``scratch/from_usg`` stay distinct.

    Manifest-less workspaces are included deliberately. Measured under
    ``examples/``: 39 MF6 workspaces against 14 ``run.json`` files, and most
    notebooks call ``sim.run_simulation`` directly. A ledger that silently omitted
    two thirds of real runs would be worse than none.
    """

    _kind: ClassVar[str] = "run"
    _columns: ClassVar[tuple[str, ...]] = (
        "key",
        "status",
        "outcome",
        "stale",
        "drift",
        "elapsed",
        "simulation",
        "models",
        "packages",
        "updated",
        "path",
    )
    # The library to compare recorded digests against. A class default rather
    # than constructor state, for the reason in SpecLibrary: dict.copy() and
    # dict.__or__ drop instance attributes.
    _packages: PackageLibrary | None = None
    _digest_memo: dict[str, str | None] | None = None
    _calibrations: dict[str, int] | None = None

    @classmethod
    def discover(
        cls, project: Project, *, load: bool = False
    ) -> RunLibrary:
        """Scan a project for runs. Never raises on a bad manifest."""

        library = cls()
        library._project_name = project.name
        library._packages = project.packages
        root = project.root
        seen: set[Path] = set()

        manifests = sorted(project.runs_dir.rglob(RUN_MANIFEST_NAME))
        for manifest_path in manifests:
            workspace = manifest_path.parent
            seen.add(workspace)
            key = _relative_key(workspace, project.runs_dir)
            try:
                run = Run.reopen(workspace, load=load)
            except (OSError, ValueError, KeyError) as error:
                # An unreadable manifest is a ROW, not an exception: one truncated
                # file must not hide every other run. `_write_json` is a
                # non-atomic write_text called several times per run, so this is
                # a real state, not a hypothetical one.
                logger.debug(
                    "Run '%s': manifest unreadable (%s); listed as unreadable.",
                    key,
                    error,
                )
                run = Run(name=workspace.name, workspace=workspace, status="unreadable")
            library[key] = run

        calibrations: dict[str, int] = {}
        for name_file in sorted(root.rglob("mfsim.nam")):
            workspace = name_file.parent
            if workspace in seen:
                continue
            # A calibration's template, masters and workers are each a full copy
            # of the model workspace. They are counted and reported, never listed
            # as runs -- one calibration would otherwise add a dozen rows to a
            # ledger whose job is telling you what you ran.
            calibration = pest_tree_root(workspace, root)
            if calibration is not None:
                key = _relative_key(calibration, root)
                calibrations[key] = calibrations.get(key, 0) + 1
                continue
            # Key by directory name when it sits under runs/, exactly as a
            # manifest-bearing run is keyed -- whether a run happens to have a
            # manifest must not change how it is addressed. Anything elsewhere in
            # the project keys by its path, so runs/from_usg and scratch/from_usg
            # stay distinct.
            base = project.runs_dir if workspace.is_relative_to(project.runs_dir) else root
            key = _relative_key(workspace, base)
            # Timestamps from the FILES, not from now. `created_at` defaults to
            # the moment the Run is constructed, so a scanned run reported the
            # scan's own wall clock in the `updated` column -- on the real
            # project two discovered runs differed by 55 microseconds, which is
            # the loop, not the runs.
            stamp = _workspace_timestamp(name_file)
            run = Run(
                name=workspace.name,
                workspace=workspace,
                status="discovered",
                created_at=stamp,
                updated_at=stamp,
            )
            if load:
                run.open()
            library[key] = run
        library._calibrations = calibrations
        return library

    @property
    def calibrations(self) -> pd.DataFrame:
        """The PEST calibration trees under this project, which are NOT runs.

        A calibration's template, IES/prior masters and every pyEMU worker are
        each a full copy of the model workspace, so each carries an
        ``mfsim.nam``. Listing them as runs floods the ledger -- one calibration
        with 8 workers adds ~11 rows -- and conflates a tool's scratch space with
        work you did. They are excluded and counted here instead, and the run
        repr says how many, so nothing is silently absent.

        Open a calibration with ``model.pest_runs`` / ``mf.find_pest_runs``,
        which return handles that understand phi, ensembles and forecasts.
        """

        import pandas as pd

        rows = [
            {"calibration": key, "workspaces": count, "open_with": "model.pest_runs"}
            for key, count in sorted((self._calibrations or {}).items())
        ]
        return pd.DataFrame(rows, columns=["calibration", "workspaces", "open_with"])

    def __repr__(self) -> str:
        """Compact one-liner, naming excluded calibrations when there are any."""

        base = super().__repr__()
        count = len(self._calibrations or {})
        if not count:
            return base
        noun = "calibration" if count == 1 else "calibrations"
        return f"{base[:-1]}; {count} {noun} excluded -- see .calibrations)"

    def _repr_html_(self) -> str:
        """The run table, with the excluded calibrations named beneath it."""

        html = super()._repr_html_()
        if not self._calibrations:
            return html
        listed = ", ".join(
            f"<code>{escape(key)}</code> ({count})"
            for key, count in sorted(self._calibrations.items())
        )
        return (
            f"{html}<p><em>Excluded as PEST calibrations, not runs:</em> {listed}. "
            "Open them with <code>model.pest_runs</code>.</p>"
        )

    def drift(self, key: str | None = None) -> pd.DataFrame:
        """Report which library packages a run used have changed since it ran.

        The manifest records package KEYS, so ``project.packages.edit("npf/base",
        ...)`` leaves every earlier run still naming ``npf/base`` -- true of the
        key, false of the model it actually used. Each run also records a
        fingerprint of what those keys MEANT at the time; this compares them.

        Five states, and the four that are not ``unchanged`` all mean "do not
        assume": ``unchanged`` (reproducible from the library as it stands),
        ``edited`` (the key still exists and means something else), ``removed``
        (the key is gone from the library), ``unfingerprinted`` (the package is
        still perfectly buildable but cannot be fingerprinted -- a closure
        builder, an option that does not pickle), and ``unknown`` for a run
        recorded before digests existed. A missing fingerprint is not evidence
        of sameness, and reporting any of the last three as ``unchanged`` would
        be the same class of lie as ``status="written"`` reading as success.

        A run keeps its own copy of these packages (``<run>/_library/``) and
        rebuilds from it, so ``edited`` does not mean the run can no longer be
        rebuilt as it ran: it means preparing it again from the library would
        give a different model. Runs prepared before the copies existed are the
        exception, and warn when they are rebuilt from a changed library.
        """

        import pandas as pd

        if self._packages is None:
            raise RuntimeError(
                "This run library is not bound to a project, so there is no "
                "package library to compare against. Use project.runs."
            )
        if key is None:
            keys = self.names()
        else:
            self[key]  # raise the did-you-mean KeyError for an unknown run
            keys = [key]
        # Fingerprint the library ONCE per call. Digesting per (run, package)
        # made summary() O(runs x packages) and re-hashed the same arrays for
        # every run -- measured ~1 ms to ~3 s on a realistic project. The memo is
        # per-call and discarded, so an edit is picked up on the next call rather
        # than needing invalidation.
        now_by_key = self._library_digest_memo()
        source_by_key = self._library_source_memo()
        rows: list[dict[str, Any]] = []
        for run_key in keys:
            recorded = self[run_key].package_digests
            if not recorded:
                rows.append(
                    {
                        "run": run_key,
                        "key": "—",
                        "ran_with": "—",
                        "library_now": "—",
                        "state": "unknown" if recorded is None else "unchanged",
                        "source": "—",
                    }
                )
                continue
            for package_key, ran_with in sorted(recorded.items()):
                now = now_by_key.get(package_key, _MISSING_KEY)
                state = _drift_state(ran_with, now)
                rows.append(
                    {
                        "run": run_key,
                        "key": package_key,
                        "ran_with": ran_with or "—",
                        "library_now": (now if now is not _MISSING_KEY else None) or "—",
                        "state": state,
                        "source": source_by_key.get(package_key, "—"),
                    }
                )
        return pd.DataFrame(
            rows, columns=["run", "key", "ran_with", "library_now", "state", "source"]
        )

    def _library_source_memo(self) -> dict[str, str]:
        """Per package key, whether the FILE it was built from still says the same.

        A separate column rather than a fifth ``state``, because it answers a different
        question. ``state`` is about the LIBRARY -- does this key still mean what the
        run used? ``source`` is about the world outside it -- is the file that key was
        built from still the one it was built from? A package can be ``unchanged`` and
        ``edited`` at once, and collapsing them would make one of the two invisible.

        Four values. ``ok`` the recorded hash matches the file now; ``edited`` it
        differs; ``missing`` the path is gone; ``unfingerprinted`` the package records a
        source but no hash -- every package declared before ledger 276, and any whose
        file could not be read at declaration. ``—`` means the package names no source
        at all, which is most of them.
        """

        answers: dict[str, str] = {}
        for key, spec in (self._packages or {}).items():
            metadata = getattr(spec, "metadata", None) or {}
            path = metadata.get("source_path")
            if path is None:
                continue
            recorded = metadata.get(SOURCE_SHA)
            if recorded is None:
                answers[key] = "unfingerprinted"
                continue
            now = source_fingerprint(path)
            answers[key] = "missing" if now is None else ("ok" if now == recorded else "edited")
        return answers

    def _library_digest_memo(self) -> dict[str, str | None]:
        """The current library fingerprints, computed once and cached per call."""

        memo = getattr(self, "_digest_memo", None)
        if memo is None:
            memo = {} if self._packages is None else self._packages.digests()
        return memo

    def _drift_label(self, entry: Run, now_by_key: dict[str, str | None]) -> str:
        """The one-cell drift summary for :meth:`summary`."""

        recorded = entry.package_digests
        if recorded is None:
            return "unknown"
        if self._packages is None:
            # Detached by dict.copy(): reporting "—" here would be a clean bill
            # of health for a run whose packages may have changed underneath it.
            return "unbound"
        if not recorded:
            return "—"
        moved = [
            package_key
            for package_key, ran_with in sorted(recorded.items())
            if _drift_state(ran_with, now_by_key.get(package_key, _MISSING_KEY))
            != "unchanged"
        ]
        return ", ".join(moved) if moved else "—"

    def summary(self) -> pd.DataFrame:
        """One row per run, with the library fingerprinted once for the whole table."""

        import pandas as pd

        memo = {} if self._packages is None else self._packages.digests()
        self._digest_memo = memo
        try:
            rows = [self._row(key, self[key]) for key in self.names()]
        finally:
            self._digest_memo = None
        return pd.DataFrame(rows, columns=list(self._columns))

    def _row(self, key: str, entry: Run) -> dict[str, Any]:
        """Summarize one run: what myflopy did, and what MF6 did."""

        outcome = entry.outcome
        block = entry.spec_block or {}
        models = [model.get("name", "?") for model in block.get("models", ())]
        packages = sum(len(model.get("packages", ())) for model in block.get("models", ()))
        return {
            "key": key,
            "status": entry.status,
            "outcome": outcome.outcome,
            "stale": outcome.stale,
            "drift": self._drift_label(entry, self._library_digest_memo()),
            "elapsed": outcome.elapsed or "—",
            "simulation": block.get("name", "—"),
            "models": ", ".join(models) or "—",
            "packages": packages or "—",
            "updated": entry.updated_at,
            "path": str(entry.workspace),
        }


def _workspace_timestamp(name_file: Path) -> str:
    """When a manifest-less run's files were last written, ISO-8601 UTC."""

    try:
        newest = max(
            (item.stat().st_mtime for item in name_file.parent.iterdir() if item.is_file()),
            default=name_file.stat().st_mtime,
        )
    except OSError:
        return _utc_now()
    return datetime.fromtimestamp(newest, timezone.utc).isoformat(timespec="seconds")


def _relative_key(workspace: Path, base: Path) -> str:
    """Key a run by its path under ``base``, falling back to its directory name."""

    try:
        relative = workspace.relative_to(base)
    except ValueError:
        return workspace.name
    return relative.as_posix() or workspace.name


@dataclass(frozen=True, slots=True)
class SaveReport:
    """What :meth:`Project.save` actually persisted, and what it could not.

    ``save()`` used to refuse the WHOLE project when any one simulation was
    un-persistable -- so a project with 14 saveable packages and one simulation
    carrying an inline array wrote nothing at all, and its library did not
    survive the session. The default is now to persist what it can and say what
    it skipped; ``strict=True`` restores the refuse-everything behaviour.

    ``simulations``/``packages``/``grids`` count what was persisted;
    ``pickled`` is how many of those packages and grids needed the pickle
    fallback rather than readable JSON, which is the number worth watching --
    a pickled entry is version-fragile. ``skipped`` pairs each unpersisted
    simulation with why.

    Truthy when anything was written, and usable anywhere a path is (``os.fspath``,
    ``open()``, ``Path(...)``) -- though it is NOT equal to a ``Path``.
    """

    project_spec: Path
    simulations: int = 0
    packages: int = 0
    grids: int = 0
    pickled: int = 0
    skipped: tuple[tuple[str, str], ...] = ()
    #: Files copied into ``inputs/`` by ``save(bundle_sources=True)``, as
    #: ``(from, to)`` pairs. Empty unless bundling was asked for.
    bundled: tuple[tuple[str, str], ...] = ()
    #: Context fields that could not be written, as ``(model, field)`` pairs. A grid,
    #: a domain and a surfaces table are not serializable (ledger 193), and were
    #: dropped in silence until ledger 257 -- so ``skipped=0`` read as "nothing was
    #: lost" while a model's build-time grid guard quietly stopped applying.
    context_dropped: tuple[tuple[str, str], ...] = ()

    def __fspath__(self) -> str:
        """The project spec path, so a report can stand in for one in ``open``/``Path``."""

        return str(self.project_spec)

    def __bool__(self) -> bool:
        """Whether anything was persisted."""

        return bool(self.simulations or self.packages or self.grids)

    def __repr__(self) -> str:
        """Compact summary naming the skip count, never the full skip list."""

        pieces = [
            f"project_spec={self.project_spec.name!r}",
            f"simulations={self.simulations}",
            f"packages={self.packages}",
            f"grids={self.grids}",
            f"pickled={self.pickled}",
            f"skipped={len(self.skipped)}",
        ]
        return f"SaveReport({', '.join(pieces)})"


@dataclass(slots=True)
class Run:
    """One materialized execution of a :class:`SimulationSpec`.

    A run owns its workspace and lifecycle. It can be created from a spec or
    reopened later from the JSON manifest and native MF6 input files.
    """

    name: str
    workspace: Path
    spec: SimulationSpec | None = None
    executable: str = "mf6"
    status: str = "created"
    created_at: str = field(default_factory=_utc_now)
    updated_at: str = field(default_factory=_utc_now)
    completed_at: str | None = None
    success: bool | None = None
    report: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    # The manifest's own ``spec`` block, kept verbatim when this run was reopened
    # rather than built. ``Run.reopen`` cannot rebuild a ``SimulationSpec`` (the
    # manifest stores a readable summary, not the recipe), so without this the
    # next ``save_manifest`` -- which ``_record_status`` fires on every state
    # change -- would overwrite the record of what the run contained with null.
    spec_block: dict[str, Any] | None = field(default=None, repr=False)
    # The full ``SimulationSpec.to_dict()`` recipe, when the spec could be
    # serialized. ``spec_block`` above is a READABLE SUMMARY -- names only -- so
    # it can say what a run contained but never rebuild it. This is what lets a
    # run reopened in a later session be rebuilt, re-run, or forked.
    spec_recipe: dict[str, Any] | None = field(default=None, repr=False)
    #: The project root a relative source path anchors to, when the run is reopened
    #: without its ``Project``.
    project_root: Path | None = field(default=None, repr=False)
    # What each library package MEANT when this run was prepared. The manifest
    # records package KEYS, so without this, editing 'npf/base' afterwards leaves
    # every run claiming to have used it -- true of the key, false of the model.
    package_digests: dict[str, str] | None = field(default=None, repr=False)
    # Where this run's OWN copy of each of those packages lives, relative to the
    # run folder (``_library/packages/...``). A rebuild resolves ``mf.ref`` from
    # these first, so editing the library later does not change what the run
    # rebuilds. None for a run prepared before the copies existed.
    package_snapshot: dict[str, str] | None = field(default=None, repr=False)
    # Each model's grid CRS as built: what a reopened run is drawn in. Unrecorded,
    # it fell back to EPSG:2927 whatever CRS the run had.
    model_crs: dict[str, str] | None = field(default=None, repr=False)
    built: BuiltSimulation | None = field(default=None, repr=False)
    simulation: Any | None = field(default=None, repr=False)
    build_context: SpecBuildContext | None = field(default=None, repr=False)
    _model_views: dict[str, SimulationBase] = field(
        default_factory=dict, init=False, repr=False
    )
    # `_record_status` re-summarizes the spec on every status change, so the
    # not-re-runnable warning would fire three or four times for one run.
    _warned_unserializable: bool = field(default=False, init=False, repr=False)
    # The copies once read -- or, for a run prepared in this session, the very
    # objects that were copied, so a build reads nothing back.
    _snapshot_specs: dict[str, PackageSpec] | None = field(
        default=None, init=False, repr=False
    )
    # Keys ``Project.reopen_run`` found no copy of: they resolve from the library
    # as it is NOW, and ``build`` checks them against the recorded fingerprints.
    _library_keys: tuple[str, ...] = field(default=(), init=False, repr=False)
    # What those keys resolved to at the last build, with fingerprints; ``write``
    # records them once the inputs on disk are rebuilt from them.
    _from_library: dict[str, tuple[PackageSpec, str | None]] = field(
        default_factory=dict, init=False, repr=False
    )

    def __post_init__(self) -> None:
        """Coerce ``workspace`` to a :class:`~pathlib.Path`."""

        self.workspace = Path(self.workspace)

    @property
    def manifest_path(self) -> Path:
        """Path to this run's lifecycle manifest."""

        return self.workspace / RUN_MANIFEST_NAME

    @property
    def pest_runs(self) -> list:
        """The PEST runs done on this run's models (anywhere under the run).

        Returns a list of :class:`~myflopy.modflow.mf6.pest.runs.PestRunHandle`;
        call ``.review()`` to open one as :class:`IesResults`. Empty until a
        calibration has been built/run with its default workspace under this run.

        The run's model is attached **lazily**, so ``review()`` can draw spatial
        maps (which need the grid) while merely listing resolves nothing: a
        listing property should not build a model view, cache it on the run, or
        raise on a run whose model is ambiguous. A run holding several models
        cannot pick one; there ``review(model=...)`` is the caller's job.

        Scans the whole run directory rather than one ``pest/`` subdirectory: a
        run gives each model its own folder, and a calibration now defaults to
        ``<model folder>.pest/<name>`` (ledger 107), so there is no single root
        that would hold every model's runs.
        """

        from myflopy.modflow.mf6.pest.runs import find_pest_runs

        return find_pest_runs(self.workspace,
                              model_factory=self._pest_review_model)

    def _pest_review_model(self):
        """This run's model for ``PestRunHandle.review()``, or None if ambiguous."""

        try:
            return self.model()
        except (KeyError, ValueError):
            # No discoverable model, or several and none nominated -- a
            # model-less review still opens, it just refuses spatial maps
            # with a message naming the fix.
            return None

    def build(self) -> Run:
        """Build this run's FloPy simulation from its specification.

        Library references resolve to the run's OWN copies first (written into
        ``<run>/_library/`` by ``prepare_run``), so a run reopened in a later
        session builds what it was prepared with even after the library was
        edited. A key it holds no copy of resolves from the build context's
        library; on a run reopened through ``Project.reopen_run`` that WARNS
        when the library's version is not the one the run recorded.
        """

        if self.spec is None:
            # A run whose packages carried inline arrays has no recipe -- `to_dict`
            # refuses an ndarray, so `spec_recipe` was written as None and this used
            # to be the end of the road. Its frozen snapshot holds the resolved
            # simulation, arrays and all, and reads no source file.
            thawed = thaw_simulation(self.workspace, run_name=self.name)
            if thawed is not None:
                logger.info("run %r rebuilt from its frozen snapshot", self.name)
            if thawed is None:
                raise ValueError(
                    f"Run {self.name!r} has neither a rebuildable recipe nor a frozen "
                    "snapshot, so it cannot be rebuilt -- it is openable and "
                    "inspectable, but summary-only. A run is in this state when it was "
                    "prepared before snapshots existed, when it was staged with "
                    "build=False and never built, or when freezing it failed. Prepare "
                    "the run again to make it rebuildable."
                )
            self.spec = thawed
        self.workspace.mkdir(parents=True, exist_ok=True)
        build_context = self._context_for_build()
        self.built = self.spec.build_flopy(
            self.workspace,
            build_context=build_context,
        )
        self.simulation = self.built.simulation
        patch_simulation_plot(self.simulation)
        self.simulation.exe_name = self.executable
        self._model_views.clear()
        self._from_library = self._check_library_fallback(build_context)
        self._freeze()
        self._record_status("built")
        return self

    def _freeze(self) -> None:
        """Snapshot the resolved simulation beside the run, best effort.

        At BUILD rather than at ``prepare_run``, which saves its manifest before
        anything is resolved. A failure here must never fail a build that worked:
        the snapshot is what makes the run rebuildable LATER, and losing it costs
        that and nothing else.
        """

        if self.spec is None:
            return
        try:
            freeze_simulation(
                self.spec,
                self.workspace,
                run_name=self.name,
                # The contexts the build RESOLVED, not the ones the spec declared: a
                # deferred grid only exists on the former.
                contexts={
                    name: built.context
                    for name, built in (getattr(self.built, "models", {}) or {}).items()
                },
                # What a lazy package actually built, so the snapshot holds arrays
                # rather than a recipe pointing back at source files (ledger 270).
                resolved={
                    name: dict(getattr(built, "resolved", {}) or {})
                    for name, built in (getattr(self.built, "models", {}) or {}).items()
                },
            )
        except Exception as error:  # noqa: BLE001 - see below; the tuple cannot be closed
            # Deliberately broad. This walks a spec built from arbitrary user
            # objects and writes to a filesystem, so it reaches numpy, json, and
            # every OSError a full disk or a read-only directory produces. The
            # model is already built and about to run; refusing that because a
            # convenience snapshot failed would be the worse outcome by far.
            # `warnings.warn`, not `logger.warning`: `_logging` attaches a
            # NullHandler to the myflopy root, which suppresses Python's fallback
            # stderr handler -- so a library log line prints NOTHING unless the
            # application configured logging. This failure is the moment a run stops
            # being rebuildable, and it hid a real defect through five commits while
            # `prepare_run` had already announced the opposite. It uses the same
            # channel `_warn_artifact_versions` uses for the same reason.
            logger.warning(
                "run %r built but could not be frozen, so it will not be rebuildable "
                "from its own snapshot", self.name, exc_info=True,
            )
            warnings.warn(
                f"Run {self.name!r} built but its frozen snapshot could not be "
                f"written ({error}), so it is not rebuildable from its own folder. "
                "The run itself is unaffected and will execute normally; "
                f"reopen_run({self.name!r}) will need the manifest's recipe instead, "
                "and has none if this model carries inline arrays.",
                stacklevel=3,
            )

    @property
    def frozen(self) -> dict[str, Any] | None:
        """What this run's frozen snapshot holds, or ``None`` when it has none."""

        return frozen_summary(self.workspace)

    @property
    def resolution(self) -> list[dict[str, Any]]:
        """One row per package: how it reached this run's mesh, and how big it is.

        The spec-side twin is ``SimulationSpec.move_report()``, which is computable
        BEFORE a build and is what you read to find out what the build will refuse.
        This one is the record: it is recomputed after the build, so a ``rebuilt``
        row can say what the package actually resolved to -- the lazy DISV's real
        ``ncpl`` exists nowhere else (ledger 271).

        `prepare_run` builds in memory, so the rows are there before `execute()`; a run
        reopened but not yet rebuilt reports the spec-side verdicts with no
        `resolved_ncpl`.
        """

        rows: list[dict[str, Any]] = []
        built_models = getattr(self.built, "models", {}) or {}
        for model in getattr(self.spec, "models", ()) or ():
            built = built_models.get(model.name)
            resolved = dict(getattr(built, "resolved", {}) or {})
            for package in model.packages:
                how, why = package_move(package)
                row: dict[str, Any] = {
                    "model": model.name,
                    "package": getattr(package, "name", "?"),
                    "how": how,
                    "why": why,
                }
                concrete = resolved.get(row["package"])
                if concrete is not None:
                    # `ncpl` is what a lazy DISV resolves TO; a lazy boundary condition
                    # resolves to RECORDS, and reading `options["ncpl"]` off one gave
                    # `None` -- a column that said nothing for every package except the
                    # discretization. Report whichever the package actually produced
                    # (ledger 280).
                    options = concrete.options or {}
                    if "ncpl" in options:
                        row["resolved_ncpl"] = options["ncpl"]
                    records = options.get("stress_period_data")
                    if isinstance(records, dict) and records:
                        first = records[min(records)]
                        row["resolved_records"] = len(first)
                rows.append(row)
        return rows

    def _context_for_build(self) -> SpecBuildContext | None:
        """The build context, with this run's own package copies taking precedence."""

        if not self.package_snapshot:
            if self.build_context is not None:
                return self.build_context
            # No Project in scope -- `Run.reopen` or `load_run`. The run still knows
            # which project it belongs to, and a RELATIVE grid source has to anchor
            # to that: without it the resolver falls back to the current working
            # directory, so the same run rebuilt from a different directory looked
            # for `inputs/domain.gpkg` beside wherever the interpreter happened to
            # be. Bundling makes relative paths the normal case (ledger 256).
            return (
                SpecBuildContext(project_root=self.project_root)
                if self.project_root is not None
                else None
            )
        if self._snapshot_specs is None:
            self._snapshot_specs = _read_package_snapshot(
                self.workspace, self.package_snapshot, run_name=self.name
            )
        base = self.build_context or SpecBuildContext(
            project_root=self.project_root,
            simulation_workspace=self.workspace,
            grid_workspace=self.workspace / "_grid",
        )
        return replace(
            base, package_specs={**base.package_specs, **self._snapshot_specs}
        )

    def _check_library_fallback(
        self, build_context: SpecBuildContext | None
    ) -> dict[str, tuple[PackageSpec, str | None]]:
        """Fingerprint the keys that resolved from the CURRENT library; warn if moved.

        Only the keys ``Project.reopen_run`` flagged -- those this run holds no
        copy of. A key whose fingerprint differs from the recorded one, or that
        was never fingerprinted, is not known to be what the run ran with.
        """

        if not self._library_keys or build_context is None:
            return {}
        recorded = self.package_digests or {}
        resolved = {
            key: (spec, package_digest(spec))
            for key, spec in build_context.package_specs.items()
            if key in self._library_keys
        }
        moved = [
            key
            for key, (_, digest) in resolved.items()
            if digest is None or recorded.get(key) != digest
        ]
        if moved:
            logger.warning(
                "Run '%s' was rebuilt with library package(s) %s as they are NOW. "
                "The run holds no copy of them -- it was prepared before runs "
                "kept one, or they could not be copied -- and they differ from "
                "what it was prepared with, or were never fingerprinted. Its "
                "results on disk came from the earlier versions; writing or "
                "executing it replaces them and records these. "
                "project.runs.drift() compares the two.",
                self.name,
                ", ".join(repr(key) for key in moved),
            )
        return resolved

    def _adopt_library_versions(self) -> None:
        """Record what the rewritten inputs came from, and keep a copy from now on."""

        packages = {key: spec for key, (spec, _) in self._from_library.items()}
        copied = _write_package_snapshot(
            self.workspace,
            packages,
            run_name=self.name,
            existing=self.package_snapshot,
        )
        self.package_snapshot = {**(self.package_snapshot or {}), **copied}
        self._snapshot_specs = {
            **(self._snapshot_specs or {}),
            **{key: packages[key] for key in copied},
        }
        self.package_digests = {
            **(self.package_digests or {}),
            **{key: digest for key, (_, digest) in self._from_library.items()},
        }
        self._library_keys = tuple(
            key for key in self._library_keys if key not in copied
        )
        self._from_library = {}

    def write(self, silent=False) -> Path:
        """Write MF6 input files, building first when necessary.
        silent defaults to False, so info will print to console"""

        if self.simulation is None:
            self.build()
        assert self.simulation is not None  # build() populates it
        self.simulation.write_simulation(silent=silent)
        if self._from_library:
            self._adopt_library_versions()
        self._record_status("written")
        return self.workspace

    def execute(self, *, silent: bool = True) -> tuple[bool, list[str]]:
        """Write the MF6 input files (if needed) and run the simulation.

        Records the run's success and MF6 report in the manifest.

        Parameters
        ----------
        silent : bool, default True
            Suppress MF6's console output during the run.

        Returns
        -------
        tuple[bool, list[str]]
            ``(success, report)`` -- whether MF6 converged, and its captured
            output lines.

        Examples
        --------
        >>> ok, report = run.execute()
        >>> if not ok: print("\\n".join(report[-20:]))
        """

        if (
            self.status not in {"written", "completed", "failed"}
            or self.simulation is None
        ):
            self.write()
        assert self.simulation is not None  # write()/build() populate it
        self.simulation.exe_name = self.executable
        success, report = self.simulation.run_simulation(silent=silent, report=True)
        self.success = bool(success)
        self.report = list(report)
        self.completed_at = _utc_now()
        self._record_status("completed" if success else "failed")
        return self.success, self.report

    def open(self, *, verbosity_level: int = 0):
        """Load this run's native MF6 simulation from disk."""

        if not (self.workspace / "mfsim.nam").exists():
            raise FileNotFoundError(
                f"No MF6 simulation found in run workspace: {self.workspace}"
            )
        self.simulation = flopy.mf6.MFSimulation.load(
            sim_ws=str(self.workspace),
            exe_name=self.executable,
            verbosity_level=verbosity_level,
        )
        # Writing this simulation writes what was on disk, not a library rebuild.
        self._from_library = {}
        patch_simulation_plot(self.simulation)
        self._model_views.clear()
        return self.simulation

    @property
    def model_names(self) -> tuple[str, ...]:
        """Return known model names for this run."""

        if self.built is not None:
            return tuple(self.built.models)
        if self.spec is not None:
            return self.spec.model_names
        if self.simulation is not None:
            try:
                return tuple(self.simulation.model_names)
            except (AttributeError, TypeError):
                # A duck-typed simulation with no `model_names`, or one whose
                # `model_names` is None. Fall through to reading the names off
                # the .nam files on disk.
                logger.debug(
                    "simulation exposes no model_names; reading them from the "
                    "workspace instead",
                )
        return tuple(
            path.stem
            for path in sorted(self.workspace.glob("*.nam"))
            if path.name.lower() != "mfsim.nam"
        )

    def _default_model_name(self) -> str:
        """Return the only known model name or ask the caller to choose."""

        names = self.model_names
        if len(names) == 1:
            return names[0]
        if not names:
            raise KeyError(f"Run '{self.name}' has no discoverable models.")
        raise ValueError(
            "Model name is required because this run has multiple models: "
            f"{', '.join(names)}"
        )

    def model(self, name: str | None = None) -> SimulationBase:
        """Return the preferred myflopy model view for this run.

        Built models return a live :class:`ModelView` with myflopy helper methods
        (``.packages``, ``.budget``, ``.plot.map(...)``, ``.diff(...)``, plus the
        kind-appropriate dependent variable: ``.hds`` on GWF, ``.conc`` on GWT,
        ``.temp`` on GWE). Reopened/file-backed models return a ``LoadedMf6Run``,
        which answers the same grammar over the written files.

        Not GWF-only: §6.0 made the view kind-aware, so this returns a usable view
        for GWT/GWE/PRT models too. Use :meth:`flopy_model` when you need the raw
        FloPy model object.

        Parameters
        ----------
        name : str, optional
            Model name; may be omitted when the run has exactly one model.

        Returns
        -------
        SimulationBase
            A :class:`ModelView` (built) or ``LoadedMf6Run`` (reopened).

        Raises
        ------
        ValueError
            If ``name`` is omitted but the run has multiple models.

        Examples
        --------
        >>> model = run.model("valley")
        >>> heads = model.hds.array()
        """

        name = self._default_model_name() if name is None else name
        if name in self._model_views:
            return self._model_views[name]

        if self.built is not None and name in self.built.models:
            view = ModelView(run=self, model_name=name)
            self._model_views[name] = view
            return view

        if (self.workspace / "mfsim.nam").exists():
            view = load_mf6_run(
                self.workspace,
                model_name=name,
                verbosity_level=0,
            )
            self._model_views[name] = view
            return view

        raise KeyError(f"Run '{self.name}' has no model named '{name}'.")

    def flopy_model(self, name: str | None = None) -> Any:
        """Return a raw FloPy model object by name."""

        name = self._default_model_name() if name is None else name
        if self.built is not None and name in self.built.models:
            return self.built.built_model(name).model
        if self.simulation is not None:
            model = self.simulation.get_model(name)
            if model is not None:
                return model
        if (self.workspace / "mfsim.nam").exists():
            self.open()
            assert self.simulation is not None  # open() populates it
            model = self.simulation.get_model(name)
            if model is not None:
                return model
        raise KeyError(f"Run '{self.name}' has no raw FloPy model named '{name}'.")

    @property
    def spec_name(self) -> str | None:
        """The simulation this run was prepared from, per its manifest."""

        return (self.spec_block or {}).get("name")

    @property
    def package_keys(self) -> dict[str, tuple[str, ...]]:
        """What each model carried, as the LIBRARY KEYS recorded at run time.

        Read from the manifest's own summary block rather than from
        :attr:`spec`, so it answers for a run reopened in a later session whose
        recipe could not be revived -- that block is the authority on what a run
        WAS. Returns ``{}`` for a run with no manifest.

        These are labels, not contents: a key says which entry was used, and
        ``project.runs.drift()`` says whether that entry still means the same
        thing.
        """

        return {
            model.get("name", "?"): tuple(model.get("packages", ()))
            for model in (self.spec_block or {}).get("models", ())
        }

    @property
    def outcome(self) -> ListingOutcome:
        """What MF6 reported in ``mfsim.lst``, read fresh each time.

        ``status`` is myflopy's verb -- what the library did to this run. This is
        MF6's verdict, and the two are independent: a run written by one session
        and executed by another sits at ``status="written"`` forever while its
        listing says ``Normal termination``.

        Deliberately NOT stored in the manifest. The listing is the authority and
        is right there; caching a verdict beside it would create a second source
        of truth that can disagree with the file it came from.
        """

        return read_listing_outcome(self.workspace, since=self._prepared_at())

    def _prepared_at(self) -> float | None:
        """When this run was prepared, as a POSIX timestamp, if it is knowable.

        ``None`` for a run discovered on disk without a manifest: its
        ``created_at`` is whenever the scan happened to construct it, so using
        that would call every such run's listing superseded.
        """

        if self.spec_block is None and self.status == "discovered":
            return None
        try:
            return datetime.fromisoformat(self.created_at).timestamp()
        except ValueError:
            return None

    def reconcile(self, *, write: bool = True) -> ListingOutcome:
        """Adopt MF6's verdict into ``status``/``success``, persisting by default.

        A run MF6 never finished -- killed, crashed, still going, never started --
        has no verdict, so nothing is changed: reporting that as failure would be
        the same class of lie as reporting ``status="written"`` as success.
        """

        outcome = self.outcome
        if outcome.succeeded is None:
            return outcome
        self.success = outcome.succeeded
        self.status = "completed" if outcome.succeeded else "failed"
        self.completed_at = self.completed_at or _utc_now()
        if write:
            self._record_status(self.status)
        return outcome

    def to_dict(self) -> dict[str, Any]:
        """Return the readable manifest representation of this run."""

        return {
            "run": {
                "name": self.name,
                "status": self.status,
                "executable": self.executable,
                "created_at": self.created_at,
                "updated_at": self.updated_at,
                "completed_at": self.completed_at,
                "success": self.success,
            },
            # A live spec re-summarizes; otherwise carry the block forward
            # unchanged. Never write null over provenance this run did not author.
            "spec": _spec_summary(self.spec)
            if self.spec is not None
            else self.spec_block,
            # The project this run belongs to, so a relative grid source anchors to
            # it when the run is reopened WITHOUT its Project. Without this,
            # `Run.reopen`/`load_run` resolved a relative path against the current
            # working directory -- which bundling makes load-bearing, since a bundled
            # recipe stores `inputs/<file>` (ledger 256).
            "project_root": str(self.project_root) if self.project_root else None,
            "spec_recipe": self._spec_recipe_block(),
            "package_digests": self.package_digests,
            "package_snapshot": self.package_snapshot,
            "crs": self._model_crs_block(),
            "metadata": self.metadata,
        }

    def _model_crs_block(self) -> dict[str, str] | None:
        """Each built model's grid CRS, or the block already on disk."""

        if self.built is None:
            return self.model_crs
        recorded = {
            name: str(built.context.grid.crs)
            for name, built in self.built.models.items()
            if getattr(built.context.grid, "crs", None) is not None
        }
        return recorded or self.model_crs

    def _spec_recipe_block(self) -> dict[str, Any] | None:
        """The rebuildable recipe for this run, or the one already on disk.

        A spec whose options embed raw arrays cannot be JSON-encoded; that is not
        an error here, it just means this run is summary-only and stays as
        rebuildable as it was before. Everything expressed as ``mf.ref`` into the
        package library serializes, which is the shape that matters -- the arrays
        live in the library and are persisted there.
        """

        if self.spec is None:
            return self.spec_recipe
        try:
            return self.spec.to_dict()
        except ValueError as error:  # not JSON-encodable -- an embedded array
            # WARNING, not DEBUG. This is the moment the run stops being
            # re-runnable, and it was announced at a level nobody sees: the
            # manifest keeps a readable summary, `reopen_run(...).spec` comes
            # back None, and nothing said so. De-duplicated because
            # `_record_status` re-summarizes on every status change, so a single
            # run would otherwise repeat it three or four times.
            if not self._warned_unserializable:
                self._warned_unserializable = True
                logger.warning(
                    "Run '%s' has no rebuildable recipe: its spec could not be "
                    "stored as one (%s), so the manifest keeps a readable summary "
                    "only and reopen_run('%s').spec will be None. It is still "
                    "re-runnable from the frozen snapshot written at build time. "
                    "Register array-bearing packages with "
                    "project.add_package(key, spec) and reference them with "
                    "mf.ref(key) if you want the recipe itself to serialize.",
                    self.name,
                    error,
                    self.name,
                )
            return self.spec_recipe

    def save_manifest(self) -> Path:
        """Persist the current lifecycle state."""

        _write_json(self.manifest_path, self.to_dict())
        return self.manifest_path

    def _record_status(self, status: str) -> None:
        """Update and persist the run status."""

        self.status = status
        self.updated_at = _utc_now()
        self.save_manifest()

    @classmethod
    def reopen(cls, workspace: Path | str, *, load: bool = True) -> Run:
        """Reopen a run from its manifest and optionally load the MF6 simulation."""

        workspace = Path(workspace)
        manifest_path = workspace / RUN_MANIFEST_NAME
        if not manifest_path.exists():
            raise FileNotFoundError(f"Run manifest not found: {manifest_path}")
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        run_data = data["run"]
        run = cls(
            name=run_data["name"],
            workspace=workspace,
            executable=run_data.get("executable", "mf6"),
            status=run_data.get("status", "created"),
            created_at=run_data.get("created_at") or _utc_now(),
            updated_at=run_data.get("updated_at") or _utc_now(),
            completed_at=run_data.get("completed_at"),
            success=run_data.get("success"),
            metadata=dict(data.get("metadata", {})),
            spec_block=data.get("spec"),
            spec_recipe=data.get("spec_recipe"),
            package_digests=data.get("package_digests"),
            package_snapshot=data.get("package_snapshot"),
            model_crs=data.get("crs"),
            project_root=(
                Path(data["project_root"]) if data.get("project_root") else None
            ),
        )
        run.spec = _spec_from_recipe(run.spec_recipe, name=run.name)
        # The frozen snapshot WINS over the manifest's recipe. A recipe is a set of
        # instructions and re-running them re-resolves: a deferred grid goes back to
        # its boundary and refinement files, and a lazy package would go back to its
        # sources. The snapshot is what this run actually built, so a reopened run
        # rebuilds what it RAN. To re-resolve, prepare the run again (ledger 251).
        thawed = thaw_simulation(workspace, run_name=run.name)
        if thawed is not None:
            run.spec = thawed
        if load:
            run.open()
        return run

    @classmethod
    def load(
        cls, workspace: Path | str, *, executable: str = "mf6", load: bool = True
    ) -> Run:
        """Open a native MF6 workspace, with or without a myflopy manifest."""

        workspace = Path(workspace)
        manifest_path = workspace / RUN_MANIFEST_NAME
        if manifest_path.exists():
            return cls.reopen(workspace, load=load)
        run = cls(
            name=workspace.name,
            workspace=workspace,
            executable=executable,
            status="loaded",
        )
        if load:
            run.open()
        return run


class Project:
    """Durable workspace + in-memory spec registry for a family of related runs.

    The lifecycle wrapper of the package-first API. A project holds reusable
    ``grid``/``package``/``simulation`` libraries (via ``add_grid`` / ``add_package``
    / ``add_simulation``) and turns a simulation spec into a :class:`Run` with
    ``prepare_run(name, sim)`` -- built in memory so you can inspect
    ``run.model(...)`` before ``run.execute()`` writes and runs MF6. Run manifests
    and MF6 files are durable on disk; the Python spec registry is in-memory
    because its builders may be arbitrary callables. Geometry lives on the model
    (``ModelContext``), not the project.

    Parameters
    ----------
    root : Path or str
        Project root directory (created lazily; runs land under ``<root>/runs``).
    name : str, optional
        Project name (defaults to the root directory name).
    layout : ProjectLayout, optional
        Custom on-disk directory layout (defaults to one rooted at ``root``).

    Examples
    --------
    >>> project = mf.Project("runs/valley", name="valley")
    >>> sim = mf.simulation(flow)                 # a SimulationSpec
    >>> run = project.prepare_run("baseline", sim)
    >>> run.model("valley").hds                   # inspect before running
    >>> ok, report = run.execute()                # writes + runs MF6
    """

    def __init__(
        self,
        root: Path | str,
        *,
        name: str | None = None,
        layout: ProjectLayout | None = None,
    ):
        """Open (or define) a project at ``root`` with empty spec registries.

        Uses ``layout`` if given (else a default rooted at ``root``); ``name``
        defaults to the root directory name. Packages/grids/simulations start
        empty and are populated via ``add_*``.
        """

        self.layout = layout or ProjectLayout(Path(root))
        self.root = self.layout.root
        self.name = name or self.root.name
        # dict SUBCLASSES: every dict operation still works and the entry type is
        # still inferred, but each also answers names()/summary() and gives a
        # did-you-mean on a missing key. See myflopy.library.
        self.packages: PackageLibrary = PackageLibrary()
        self.grids: GridLibrary = GridLibrary()
        self.simulations: SimulationLibrary = SimulationLibrary()
        for library in (self.packages, self.grids, self.simulations):
            library._project_name = self.name
        self._runs: RunLibrary | None = None

    @property
    def runs_dir(self) -> Path:
        """Directory containing project-owned run workspaces."""

        return self.root / "runs"

    @property
    def manifest_path(self) -> Path:
        """Path to the project workspace manifest written by :meth:`save`."""

        return self.layout.manifest_path

    def add_package(
        self, key: str, package: PackageSpec, *, overwrite: bool = False
    ) -> PackageSpec:
        """Add or replace a reusable package spec under a project ``key``.

        The key is a free-form library address that models reference later, for
        example ``mf.ref("npf/base")`` or a semantic key like
        ``mf.ref("k/calibrated")``. The package keeps its own ``PackageSpec.name``
        (the MF6 package it builds); the key does not have to match it. Resolved
        packages are checked for name collisions when a run is built.
        """

        _refuse_wrong_library(package, PackageSpec, "add_package")
        _refuse_overwrite(self.packages, key, "add_package", overwrite, self.name)
        self.packages[key] = package
        return package

    def add_grid(
        self, key: str, grid: Any, *, overwrite: bool = False
    ) -> GridSpec | GridRef:
        """Add or replace a reusable grid under a project ``key``.

        The grid may be a ``GridSpec`` recipe or an already-built grid object
        (such as a ``VoronoiGridPlus``). Models reference it with
        ``mf.grid_ref(key)``. The grid is normalized to a ``GridSpec`` and
        resolved into a model only when a run is built.

        The return -- and ``project.grids[key]`` -- is annotated
        ``GridSpec | GridRef`` rather than ``Any`` so an editor can actually
        resolve ``project.grids[key].resolve(...)``; with ``Any`` it falls back
        to guessing among every ``resolve`` it has indexed. In practice a
        library entry is a ``GridSpec``: a ``GridRef`` points INTO this library,
        so storing one here is accepted by ``add_grid`` and rejected later by
        the build (:func:`_library_grid_spec`).
        """

        _refuse_wrong_library(grid, GridSpec, "add_grid")
        _refuse_overwrite(self.grids, key, "add_grid", overwrite, self.name)
        grid = _grid_entry(grid)
        self.grids[key] = grid
        return grid

    def add_simulation(
        self, simulation: SimulationSpec, *, overwrite: bool = False
    ) -> SimulationSpec:
        """Add or replace a reusable simulation specification."""

        _refuse_wrong_library(simulation, SimulationSpec, "add_simulation")
        _refuse_overwrite(
            self.simulations, simulation.name, "add_simulation", overwrite, self.name
        )
        reason = self._why_unserializable(simulation)
        if reason is not None:
            # Said HERE, at the line that registered it, rather than discovered
            # as `SaveReport(simulations=0, skipped=1)` -- which is truthy, raises
            # nothing, and is what left a real project's `specs/` directory empty.
            logger.warning(
                "Simulation '%s' will NOT persist: %s. project.save() will skip "
                "it and Project.load() will not see it. Register array-bearing "
                "packages with project.add_package(key, spec) and reference them "
                "with mf.ref(key).",
                simulation.name,
                reason,
            )
        self.simulations[simulation.name] = simulation
        return simulation

    def add_simulation_from_yaml(self, source: str | Path) -> SimulationSpec:
        """Load a simulation from a YAML file (or text) and register it (plan §5.6).

        Thin wrapper over :meth:`SimulationSpec.from_yaml` + :meth:`add_simulation`;
        ``source`` is a ``.yaml`` path (``str``/``Path``) or raw YAML text.
        """

        return self.add_simulation(SimulationSpec.from_yaml(source))

    @property
    def runs(self) -> RunLibrary:
        """Every run under this project, scanned once and cached.

        Runs live on disk rather than in the registries, so this is a scan --
        cheap (measured 0 ms on a 4-run project, 2 ms over 39 workspaces), but a
        scan. It is cached because ``project.runs["baseline"]`` must return the
        same object twice; call :meth:`Project.refresh_runs` after executing a run
        from somewhere else.

        Reopening a single run by name stays :meth:`reopen_run`, which also
        attaches the build context so a revived spec can resolve its library
        references. This is the browse surface, not a second way to open one.
        """

        if self._runs is None:
            self._runs = RunLibrary.discover(self)
        return self._runs

    def refresh_runs(self) -> RunLibrary:
        """Re-scan the project for runs and return the refreshed library."""

        self._runs = RunLibrary.discover(self)
        return self._runs

    def package(self, key: str) -> PackageSpec:
        """Return one registered package specification by library key.

        The symmetric partner of :meth:`simulation` and :meth:`grid`.
        Equivalent to ``project.packages[key]``, which raises the same
        did-you-mean ``KeyError``; this spelling exists because the retired
        ``ProjectSpec`` had ``package_spec(key)`` and the 2026-06-17 fold into
        ``Project`` kept ``simulation(name)`` while dropping it.
        """

        return self.packages[key]

    def grid(self, key: str) -> GridSpec | GridRef:
        """Return one registered grid by library key.

        Equivalent to ``project.grids[key]``. The union return is deliberate: a
        ``GridRef`` points INTO this library and is rejected at build time, so
        narrowing it is worth something (ledger 171).
        """

        return self.grids[key]

    def summary(self) -> pd.DataFrame:
        """Return one table covering all three spec libraries.

        A ``library`` column distinguishes them; per-library columns are collapsed
        into ``detail`` so the three shapes can share one frame. Use
        ``project.packages.summary()`` for the full per-tier columns.
        """

        import pandas as pd

        rows: list[dict[str, Any]] = []
        for library in (self.packages, self.grids, self.simulations):
            for key in library.names():
                row = library._row(key, library[key])
                identity = row.get("key") or row.get("name")
                detail = ", ".join(
                    f"{name}={value}"
                    for name, value in row.items()
                    if name not in {"key", "name"}
                )
                rows.append(
                    {"library": library._kind, "key": identity, "detail": detail}
                )
        return pd.DataFrame(rows, columns=["library", "key", "detail"])

    def __repr__(self) -> str:
        """Counts only -- no disk, no validate(), no build.

        ``len()`` on three in-memory dicts. Runs are deliberately NOT counted:
        they live on disk, and PyCharm calls ``__repr__`` on every debugger step.
        """

        return (
            f"Project({self.name!r}, root={str(self.root)!r}, "
            f"packages={len(self.packages)}, grids={len(self.grids)}, "
            f"simulations={len(self.simulations)})"
        )

    def _repr_html_(self) -> str:
        """Return the three library tables. Still no disk, still no validate()."""

        blocks = [
            f"<h4>Project: <code>{escape(self.name)}</code></h4>",
            f"<p><code>{escape(str(self.root))}</code></p>",
        ]
        for library in (self.packages, self.grids, self.simulations):
            blocks.append(library._repr_html_())
        return "<div>" + "".join(blocks) + "</div>"

    def simulation(self, name: str) -> SimulationSpec:
        """Return one registered simulation specification by name."""

        try:
            return self.simulations[name]
        except KeyError as error:
            raise KeyError(
                f"Project '{self.name}' has no simulation named '{name}'."
            ) from error

    def _build_context(self, workspace: Path) -> SpecBuildContext:
        """Return the build context that resolves project package references.

        The project's package library is passed through so that any
        ``PackageRef`` declared on a model resolves to its concrete
        ``PackageSpec`` at build time.
        """

        return SpecBuildContext(
            project_root=self.root,
            simulation_workspace=workspace,
            grid_workspace=workspace / "_grid",
            package_specs=dict(self.packages),
            grid_specs={key: _library_grid_spec(key, value) for key, value in self.grids.items()},
        )

    def prepare_run(
        self,
        name: str,
        simulation: SimulationSpec | str,
        *,
        executable: str = "mf6",
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
        build: bool = True,
    ) -> Run:
        """Create a run and build its FloPy simulation in memory, ready to inspect.

        The returned run is built but **not** written or executed -- no MF6 input
        files exist yet -- so you can inspect it with ``run.model(...)`` (and your
        plotting/explorer tools) before committing to ``run.execute()``.

        Each library package the simulation references by ``mf.ref`` is copied
        into ``<run>/_library/``, and the run rebuilds from those copies: editing
        the library afterwards does not change what reopening this run in a
        later session builds. ``project.runs.drift()`` still says when the
        library has moved on from it.

        Parameters
        ----------
        name : str
            Run name; its workspace is ``<project>/runs/<name>``.
        simulation : SimulationSpec or str
            The simulation spec to build, or the name of one registered via
            :meth:`add_simulation`.
        executable : str, default "mf6"
            MF6 executable name/path used when the run is executed.
        metadata : dict, optional
            Arbitrary metadata recorded in the run manifest.
        overwrite : bool, default False
            Allow reusing/overwriting an existing run workspace of the same name.
        build : bool, default True
            Build the FloPy simulation now; ``False`` defers building (e.g. when
            staging several runs to build/execute selectively later).

        Returns
        -------
        Run
            A built (unless ``build=False``), unexecuted run.

        Examples
        --------
        >>> run = project.prepare_run("baseline", sim)
        >>> run.model("valley")        # inspect
        >>> run.execute()              # write + run MF6
        """

        spec = (
            self.simulations[simulation] if isinstance(simulation, str) else simulation
        )
        spec = _materialize_models(spec)
        workspace = self.runs_dir / name
        if workspace.exists() and any(workspace.iterdir()) and not overwrite:
            raise FileExistsError(f"Run workspace already exists: {workspace}")
        workspace.mkdir(parents=True, exist_ok=True)
        # A reused directory keeps whatever the PREVIOUS simulation froze there, and a
        # snapshot that outlives its simulation rebuilds the wrong model in silence.
        # Dropping it here also leaves a run staged with build=False holding no answer
        # it was never asked for (ledger 251).
        if discard_snapshot(workspace):
            logger.debug("discarded run %r's previous frozen snapshot", name)
        references = self._library_references(spec)
        snapshot = _write_package_snapshot(workspace, references, run_name=name)
        run = Run(
            name=name,
            workspace=workspace,
            spec=spec,
            executable=executable,
            metadata={} if metadata is None else dict(metadata),
            build_context=self._build_context(workspace),
            package_digests=self._package_digests(spec),
            package_snapshot=snapshot,
            project_root=self.root,
        )
        # Build from the very objects just copied, rather than reading them back.
        run._snapshot_specs = {key: references[key] for key in snapshot}
        run.save_manifest()
        if build:
            run.build()
        return run

    def run(
        self,
        name: str,
        simulation: SimulationSpec | str,
        *,
        executable: str = "mf6",
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
        silent: bool = True,
    ) -> Run:
        """Create, build, write, and execute one named simulation run."""

        run = self.prepare_run(
            name,
            simulation,
            executable=executable,
            metadata=metadata,
            overwrite=overwrite,
        )
        run.execute(silent=silent)
        return run

    def reopen_run(self, name: str, *, load: bool = True) -> Run:
        """Reopen one project-owned run by name.

        The run is given this project's build context, so a spec revived from the
        manifest can resolve its ``mf.ref`` references and be rebuilt or re-run.
        Package references resolve to the run's OWN copies, written by
        :meth:`prepare_run`, so what it rebuilds does not follow later edits to
        the library. A run prepared before those copies existed resolves against
        the library as it is now, and warns on rebuild when that is not what it
        recorded. Reopening through :func:`load_run` or :meth:`Run.reopen`
        resolves package references from the run's copies and nothing else --
        there is no library in scope.
        """

        run = Run.reopen(self.runs_dir / name, load=load)
        # Unconditionally: a run whose manifest holds no recipe rebuilds from its
        # frozen snapshot instead, and that spec needs the same context to resolve
        # `mf.grid_ref` and any reference it holds no copy of. Guarding this on
        # `spec is not None` left exactly the runs the snapshot rescues without one.
        run.build_context = self._build_context(run.workspace)
        if run.spec is not None:
            copied = run.package_snapshot or {}
            run._library_keys = tuple(
                key for key in _package_ref_keys(run.spec) if key not in copied
            )
        return run

    def discover_runs(self, *, load: bool = False) -> list[Run]:
        """Discover project-owned runs that have lifecycle manifests."""

        return [
            Run.reopen(path.parent, load=load)
            for path in sorted(self.runs_dir.glob(f"*/{RUN_MANIFEST_NAME}"))
        ]

    def discover_native_runs(
        self,
        search_root: Path | str | None = None,
        *,
        load: bool = False,
    ) -> list[Run]:
        """Discover native MF6 workspaces, including runs without manifests."""

        root = self.root if search_root is None else Path(search_root)
        runs = []
        for name_file in sorted(root.rglob("mfsim.nam")):
            workspace = name_file.parent
            manifest_path = workspace / RUN_MANIFEST_NAME
            if manifest_path.exists():
                run = Run.reopen(workspace, load=load)
            else:
                run = Run(
                    name=workspace.name,
                    workspace=workspace,
                    status="discovered",
                )
                if load:
                    run.open()
            runs.append(run)
        return runs

    def _package_digests(self, spec: SimulationSpec) -> dict[str, str]:
        """Fingerprint what each LIBRARY package meant for this simulation.

        Only ``mf.ref`` entries are recorded. A package written inline is baked
        into the run's own recipe and cannot drift; a reference names a library
        entry that can be edited afterwards, which is the whole point of
        recording it. A ref whose key is not defined is skipped -- ``validate``
        already reports that, and a digest is not the place to raise it.

        A package that cannot be fingerprinted (a closure builder, an option that
        does not pickle) records ``None`` rather than being dropped, so drift
        says "unfingerprinted" instead of silently under-reporting. It must never
        make ``prepare_run`` fail: the package was buildable before drift existed
        and still is.
        """

        library = self.packages.digests()  # cached until the library changes
        return {key: library[key] for key in self._library_references(spec)}

    def _library_references(self, spec: SimulationSpec) -> dict[str, PackageSpec]:
        """The library packages ``spec`` references by ``mf.ref``, in order.

        A ref whose key is not defined is left out, as :meth:`_package_digests`
        does -- ``validate`` reports it.
        """

        return {
            key: self.packages[key]
            for key in _package_ref_keys(spec)
            if key in self.packages
        }

    def _library_digests(self) -> dict[str, str | None]:
        """Fingerprint the whole package library once, for a batch comparison."""

        return self.packages.digests()

    def _why_unserializable(self, simulation: SimulationSpec) -> str | None:
        """Return a readable reason a simulation cannot be persisted, or None."""

        try:
            simulation.to_dict()
        except ValueError as error:
            return str(error)
        return None

    def _unresolved_refs(self, simulation: SimulationSpec) -> list[str]:
        """Return issues for package/grid references the library cannot resolve.

        A saved project must be rebuildable, so every ``mf.ref`` / ``mf.grid_ref``
        in a simulation must point at a defined library entry.
        """

        issues: list[str] = []
        package_keys = set(self.packages)
        grid_keys = set(self.grids)

        def check_packages(packages, where: str) -> None:
            """Record an issue for any package ref whose key is not defined in the project."""

            for package in packages:
                if isinstance(package, PackageRef) and package.key not in package_keys:
                    issues.append(
                        f"{where} references undefined package '{package.key}'."
                    )

        check_packages(simulation.packages, f"Simulation '{simulation.name}'")
        for model in simulation.models:
            where = f"Model '{model.name}' in simulation '{simulation.name}'"
            check_packages(model.packages, where)
            if isinstance(model.grid, GridRef) and model.grid.key not in grid_keys:
                issues.append(f"{where} references undefined grid '{model.grid.key}'.")
        return issues

    def validate(self) -> list[str]:
        """Return human-readable persistence issues without touching disk.

        Persistence requires serializable, rebuildable specs. Simulations that
        embed raw arrays or exchanges, derive from unknown simulations, or
        reference undefined library packages/grids are reported here rather than
        failing silently at :meth:`save` (or later at build).
        """

        issues: list[str] = []
        names = set(self.simulations)
        for simulation in self.simulations.values():
            if (
                simulation.derived_from is not None
                and simulation.derived_from not in names
            ):
                issues.append(
                    f"Simulation '{simulation.name}' derives from unknown "
                    f"simulation '{simulation.derived_from}'."
                )
            reason = self._why_unserializable(simulation)
            if reason is not None:
                issues.append(
                    f"Simulation '{simulation.name}' is not serializable: {reason}"
                )
            issues.extend(self._unresolved_refs(simulation))
        return issues

    def _project_dict(self, *, simulations: Iterable[str] | None = None) -> dict[str, Any]:
        """Return the durable project document.

        ``simulations`` lists the names actually written. It must exclude anything
        :meth:`save` skipped: :meth:`load` raises ``FileNotFoundError`` for a
        listed entry whose file is missing, which is a contract worth keeping --
        so a skipped simulation is left OUT of the listing rather than listed and
        absent.
        """

        return {
            "kind": "Project",
            "name": self.name,
            "layout": self.layout.to_dict(),
            "simulations": sorted(
                self.simulations if simulations is None else simulations
            ),
            "packages": sorted(self.packages),
            "grids": sorted(self.grids),
        }

    def _save_grid(self, key: str, grid: Any) -> bool:
        """Persist one project-library grid, pickling built grid objects.

        Returns whether it fell back to pickling, so ``SaveReport.pickled``
        counts grids as its docstring has always claimed.
        """

        if isinstance(grid, GridSpec) and grid.method == "object":
            pickle_path = self.layout.grid_pickle_path(key)
            pickle_path.parent.mkdir(parents=True, exist_ok=True)
            with pickle_path.open("wb") as handle:
                pickle.dump(grid.obj, handle)
            _write_json(self.layout.grid_sidecar_path(key), _artifact_versions())
            doc = {
                "kind": "GridSpec",
                "name": grid.name,
                "grid_type": grid.grid_type,
                "method": "pickle",
                "pickle_path": pickle_path.relative_to(self.root).as_posix(),
            }
            _write_json(self.layout.grid_spec_path(key), doc)
            return True
        # Recipe GridSpec or GridRef: already serializable.
        _write_json(self.layout.grid_spec_path(key), grid.to_dict())
        return False

    def _save_package(self, key: str, package: PackageSpec) -> bool:
        """Persist one project package, pickling it when it carries arrays.

        Source-driven / scalar packages serialize to readable JSON. A package
        whose options embed arrays (a computed ``k`` field, explicit
        stress-period data, etc.) is pickled as an artifact instead, mirroring
        how built grids are persisted, so it can be reused across variants
        without recomputing. Returns whether it fell back to pickling.
        """

        return _persist_package(
            package,
            spec_path=self.layout.package_spec_path(key),
            pickle_path=self.layout.package_pickle_path(key),
            sidecar_path=self.layout.package_sidecar_path(key),
            root=self.root,
        )

    def save(self, *, strict: bool = False, bundle_sources: bool = False) -> SaveReport:
        """Persist the project document, its simulations, and its package library.

        Persistence is opt-in. Returns a :class:`SaveReport` naming what was
        written and what could not be.

        A simulation that embeds a raw array, derives from an unregistered
        simulation, or references an undefined library entry cannot be
        persisted. By default such a simulation is SKIPPED and everything else is
        written -- because refusing the whole project loses a library whose every
        package would have saved.

        A simulation embedding a raw array is skipped with a logged warning: there
        is nothing to fix at save time. A simulation with a DANGLING reference is
        skipped and then RAISES, because that is a mistake with a fix -- but the
        raise comes after the write, so the library persists either way.
        ``strict=True`` restores the previous behaviour and raises before writing
        anything.

        Packages are unaffected by a skipped simulation: :meth:`_save_package`
        already pickles a package whose ``to_dict`` raises, so an array-bearing
        package persists either way.

        ``bundle_sources=True`` copies every file the project's GRIDS are built from
        into ``inputs/`` and repoints the saved recipes at those copies, so the
        project is an archive rather than a set of references into the directory it
        was authored in. It repoints the LIVE registry too: otherwise a later plain
        ``save()`` would write the original paths back and quietly un-bundle the
        project. Package sources are deliberately not bundled -- they are resolved
        eagerly, so a saved package holds records rather than a path, and copying
        them would archive files nothing opens.
        """

        bundled: list[tuple[str, str]] = []
        if bundle_sources:
            bundled = self._bundle_grid_sources()

        dropped = self._report_dropped_contexts()

        issues = self.validate()
        if issues and strict:
            raise ValueError(
                "Project cannot be saved:\n"
                + "\n".join(f"- {issue}" for issue in issues)
            )

        # Two kinds of un-persistable simulation, handled differently on purpose.
        # An embedded array is nothing the caller can fix at save time, so it is
        # skipped and logged. A dangling reference is a MISTAKE -- the fix is to
        # define the entry -- so it still raises. Either way the write happens
        # first, so one bad simulation no longer costs the whole library.
        skipped: list[tuple[str, str]] = []
        integrity: list[str] = []
        for name, simulation in self.simulations.items():
            reason = self._why_unserializable(simulation)
            broken = self._unresolved_refs(simulation)
            if simulation.derived_from is not None and (
                simulation.derived_from not in self.simulations
            ):
                broken.append(
                    f"Simulation '{name}' derives from unknown simulation "
                    f"'{simulation.derived_from}'."
                )
            if broken:
                integrity.extend(broken)
                skipped.append((name, "; ".join(broken)))
            elif reason is not None:
                skipped.append((name, reason))

        skipped_names = {name for name, _ in skipped}

        # Checked BEFORE the write, unlike the dangling-ref case below: two keys
        # that name one file cannot be written harmlessly -- the second silently
        # overwrites the first, and the reload hands the survivor back under both
        # keys. Refusing costs a rename; writing costs a spec.
        clashes = [
            *_colliding_keys(self.packages, self.layout.package_spec_path),
            *_colliding_keys(self.grids, self.layout.grid_spec_path),
        ]
        if clashes:
            raise ValueError(
                "Project cannot be saved -- these library keys map to the same "
                "file, so one would silently overwrite the other. Rename one of "
                "each pair:\n"
                + "\n".join(
                    f"- '{first}' and '{second}' both -> "
                    f"{path.relative_to(self.root).as_posix()}"
                    for first, second, path in clashes
                )
            )

        self.layout.ensure()
        _write_json(
            self.manifest_path,
            {
                "project": {
                    "name": self.name,
                    "spec": f"{self.layout.specs_dir_name}/{PROJECT_SPEC_NAME}",
                }
            },
        )
        written_names = [
            name for name in self.simulations if name not in skipped_names
        ]
        _write_json(
            self.layout.project_spec_path,
            self._project_dict(simulations=written_names),
        )
        for name in written_names:
            _write_json(
                self.layout.simulation_spec_path(name),
                self.simulations[name].to_dict(),
            )
        pickled = 0
        for key, package in self.packages.items():
            if self._save_package(key, package):
                pickled += 1
        for key, grid in self.grids.items():
            if self._save_grid(key, grid):
                pickled += 1
        for name, reason in skipped:
            logger.warning(
                "Project '%s': simulation '%s' not persisted -- %s",
                self.name,
                name,
                reason,
            )
        if integrity:
            raise ValueError(
                "Project saved, but these simulations were left out and must be "
                "fixed:\n" + "\n".join(f"- {issue}" for issue in integrity)
            )
        return SaveReport(
            project_spec=self.layout.project_spec_path,
            simulations=len(written_names),
            packages=len(self.packages),
            grids=len(self.grids),
            pickled=pickled,
            skipped=tuple(skipped),
            bundled=tuple(bundled),
            context_dropped=tuple(dropped),
        )

    def _report_dropped_contexts(self) -> list[tuple[str, str]]:
        """Name the context fields this save cannot write, once per save.

        Warned HERE rather than in ``ModelSpec.to_dict``: ``to_dict`` is called by
        ``validate``, by ``_why_unserializable``, by ``to_yaml`` and by the save
        itself, so a warning there fires three or four times for one save -- the same
        de-duplication problem ledger 227 hit with its own warning.
        """

        dropped = [
            (f"{name}.{model.name}", field)
            for name, simulation in self.simulations.items()
            for model in simulation.models
            for field in dropped_context_fields(model.context)
        ]
        if dropped:
            fields = sorted({field for _model, field in dropped})
            models = sorted({model for model, _field in dropped})
            logger.warning(
                "Saving %r: %d model context(s) carry %s, which a project file cannot "
                "hold -- they reload as None (%s). Register the grid with "
                "project.add_grid(key, grid) and declare the model as "
                "mf.gwf(..., grid=mf.grid_ref(key)), which persists and resolves at "
                "build time, so the full build-time geometry check and the per-package "
                "grid pins still apply. Without that, a fingerprint of the cell "
                "NUMBERING is recorded and compared instead: a different mesh is still "
                "caught, but one with the same connectivity and different coordinates "
                "-- a rescale, a mirror, a unit conversion -- is not.",
                self.name, len(models), " and ".join(fields), ", ".join(models[:5]),
            )
        return dropped

    def _bundle_grid_sources(self) -> list[tuple[str, str]]:
        """Copy every grid source into ``inputs/`` and repoint the project at it.

        Both tiers: the grid library, and a ``GridSpec`` written inline on a model.
        A ``GridRef`` is followed to its library entry rather than copied twice, and
        a ``method='object'`` grid has no sources -- it is already pickled whole.
        """

        inputs = self.layout.inputs_dir
        inputs.mkdir(parents=True, exist_ok=True)
        copied: list[tuple[str, str]] = []

        def _bundle(spec: GridSpec) -> GridSpec:
            if spec.method == "object":
                return spec
            result, report = bundle_grid_sources(
                spec, root=self.root, inputs_dir=inputs, project_root=self.root
            )
            copied.extend((item["from"], item["to"]) for item in report["copied"])
            return result

        for key, grid in list(self.grids.items()):
            if isinstance(grid, GridSpec):
                self.grids[key] = _bundle(grid)

        for name, simulation in list(self.simulations.items()):
            models = tuple(
                replace(model, grid=_bundle(model.grid))
                if isinstance(model.grid, GridSpec)
                else model
                for model in simulation.models
            )
            if models != simulation.models:
                self.simulations[name] = replace(simulation, models=models)

        if copied:
            logger.info("bundled %d grid source file(s) into %s", len(copied), inputs)
        return copied

    @staticmethod
    def _restore_layout(root: Path) -> ProjectLayout | None:
        """Rebuild a saved custom layout via the fixed project manifest, if any."""

        manifest_path = root / PROJECT_MANIFEST_NAME
        if not manifest_path.exists():
            return None
        spec_rel = _read_json(manifest_path).get("project", {}).get("spec")
        spec_path = root / spec_rel if spec_rel else None
        if spec_path is None or not spec_path.exists():
            return None
        layout_block = _read_json(spec_path).get("layout")
        if not layout_block:
            return None
        return ProjectLayout.from_dict(root, layout_block)

    @classmethod
    def load(cls, root: Path | str, *, layout: ProjectLayout | None = None) -> Project:
        """Load a project previously written with :meth:`save`.

        A saved custom layout round-trips automatically. A listed package, grid,
        or simulation whose file is missing raises rather than loading a broken
        project silently.
        """

        root = Path(root)
        if layout is None:
            layout = cls._restore_layout(root)
        project = cls(root, layout=layout)
        spec_path = project.layout.project_spec_path
        if not spec_path.exists():
            raise FileNotFoundError(f"Project spec not found: {spec_path}")
        data = _read_json(spec_path)
        project.name = data.get("name", project.name)
        for key in data.get("packages", ()):
            package_path = project.layout.package_spec_path(key)
            if not package_path.exists():
                raise FileNotFoundError(
                    f"Project lists package '{key}' but its spec is missing: {package_path}"
                )
            project.packages[key] = project._load_package(key, _read_json(package_path))
        for name in data.get("simulations", ()):
            sim_path = project.layout.simulation_spec_path(name)
            if not sim_path.exists():
                raise FileNotFoundError(
                    f"Project lists simulation '{name}' but its spec is missing: {sim_path}"
                )
            project.add_simulation(SimulationSpec.from_dict(_read_json(sim_path)))
        for key in data.get("grids", ()):
            grid_path = project.layout.grid_spec_path(key)
            if not grid_path.exists():
                raise FileNotFoundError(
                    f"Project lists grid '{key}' but its spec is missing: {grid_path}"
                )
            project.grids[key] = project._load_grid(key, _read_json(grid_path))
        return project

    def _load_grid(self, key: str, doc: dict[str, Any]) -> GridSpec | GridRef:
        """Reconstruct one project-library grid from its document."""

        if doc.get("kind") == "GridRef":
            return GridRef.from_dict(doc)
        if doc.get("method") == "pickle":
            self._warn_version_mismatch("grid", key, self.layout.grid_sidecar_path(key))
        return GridSpec.from_dict(doc)

    def _load_package(self, key: str, doc: dict[str, Any]) -> PackageSpec:
        """Reconstruct one project package from its document."""

        return _revive_package(
            key,
            doc,
            root=self.root,
            sidecar_path=self.layout.package_sidecar_path(key),
            stacklevel=4,
        )

    def _warn_version_mismatch(self, kind: str, key: str, sidecar_path: Path) -> None:
        """Warn if a pickled artifact was written under a different library stack."""

        _warn_artifact_versions(kind, key, sidecar_path, stacklevel=4)
