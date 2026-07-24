"""Project and run lifecycle objects for spec-driven MODFLOW workflows."""

from __future__ import annotations

import json
import pickle
import re
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib import import_module
from pathlib import Path
from typing import Any

import flopy

from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.project.run_model import load_mf6_run, patch_simulation_plot
from myflopy.specs import (
    BuiltSimulation,
    GridRef,
    GridSpec,
    PackageRef,
    PackageSpec,
    SimulationSpec,
    SpecBuildContext,
    _grid_entry,
    _package_entry_label,
)

RUN_MANIFEST_NAME = "run.json"
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
        except Exception:
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

        return package.metadata if isinstance(package, PackageSpec) else {}

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
                "context_metadata": model.context.metadata,
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
    built: BuiltSimulation | None = field(default=None, repr=False)
    simulation: Any | None = field(default=None, repr=False)
    build_context: SpecBuildContext | None = field(default=None, repr=False)
    _model_views: dict[str, SimulationBase] = field(
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
        """The PEST runs done on this run's model (``<workspace>/pest``).

        Returns a list of :class:`~myflopy.modflow.mf6.pest.runs.PestRunHandle`;
        call ``.review()`` to open one as :class:`IesResults`. Empty until a
        calibration has been built/run with its default workspace under this run.
        """

        from myflopy.modflow.mf6.pest.runs import find_pest_runs

        return find_pest_runs(self.workspace / "pest")

    def build(self) -> Run:
        """Build this run's FloPy simulation from its specification."""

        if self.spec is None:
            raise ValueError("A reopened run without a spec cannot be rebuilt.")
        self.workspace.mkdir(parents=True, exist_ok=True)
        self.built = self.spec.build_flopy(
            self.workspace,
            build_context=self.build_context,
        )
        self.simulation = self.built.simulation
        patch_simulation_plot(self.simulation)
        self.simulation.exe_name = self.executable
        self._model_views.clear()
        self._record_status("built")
        return self

    def write(self) -> Path:
        """Write MF6 input files, building first when necessary."""

        if self.simulation is None:
            self.build()
        assert self.simulation is not None  # build() populates it
        self.simulation.write_simulation(silent=True)
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
            except Exception:
                pass
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
        """Return the preferred myflopy GWF model view for this run.

        Built GWF models return a live :class:`ModelView` with myflopy helper
        methods (``.hds``, ``.packages``, ``.cor(...)``, ``.diff(...)``, ...).
        Reopened/file-backed GWF models return a ``LoadedMf6Run``.

        Use :meth:`flopy_model` when you need a raw FloPy model object or are
        working with a non-GWF model type.

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
            "spec": None if self.spec is None else _spec_summary(self.spec),
            "metadata": self.metadata,
        }

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
        )
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
        self.packages: dict[str, PackageSpec] = {}
        self.grids: dict[str, Any] = {}
        self.simulations: dict[str, SimulationSpec] = {}

    @property
    def runs_dir(self) -> Path:
        """Directory containing project-owned run workspaces."""

        return self.root / "runs"

    @property
    def manifest_path(self) -> Path:
        """Path to the project workspace manifest written by :meth:`save`."""

        return self.layout.manifest_path

    def add_package(self, key: str, package: PackageSpec) -> PackageSpec:
        """Add or replace a reusable package spec under a project ``key``.

        The key is a free-form library address that models reference later, for
        example ``mf.ref("npf/base")`` or a semantic key like
        ``mf.ref("k/calibrated")``. The package keeps its own ``PackageSpec.name``
        (the MF6 package it builds); the key does not have to match it. Resolved
        packages are checked for name collisions when a run is built.
        """

        self.packages[key] = package
        return package

    def add_grid(self, key: str, grid: Any) -> Any:
        """Add or replace a reusable grid under a project ``key``.

        The grid may be a ``GridSpec`` recipe or an already-built grid object
        (such as a ``VoronoiGridPlus``). Models reference it with
        ``mf.grid_ref(key)``. The grid is normalized to a ``GridSpec`` and
        resolved into a model only when a run is built.
        """

        grid = _grid_entry(grid)
        self.grids[key] = grid
        return grid

    def add_simulation(self, simulation: SimulationSpec) -> SimulationSpec:
        """Add or replace a reusable simulation specification."""

        self.simulations[simulation.name] = simulation
        return simulation

    def add_simulation_from_yaml(self, source: str | Path) -> SimulationSpec:
        """Load a simulation from a YAML file (or text) and register it (plan §5.6).

        Thin wrapper over :meth:`SimulationSpec.from_yaml` + :meth:`add_simulation`;
        ``source`` is a ``.yaml`` path (``str``/``Path``) or raw YAML text.
        """

        return self.add_simulation(SimulationSpec.from_yaml(source))

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
        run = Run(
            name=name,
            workspace=workspace,
            spec=spec,
            executable=executable,
            metadata={} if metadata is None else dict(metadata),
            build_context=self._build_context(workspace),
        )
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
        """Reopen one project-owned run by name."""

        return Run.reopen(self.runs_dir / name, load=load)

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

    def _project_dict(self) -> dict[str, Any]:
        """Return the durable project document."""

        return {
            "kind": "Project",
            "name": self.name,
            "layout": self.layout.to_dict(),
            "simulations": sorted(self.simulations),
            "packages": sorted(self.packages),
            "grids": sorted(self.grids),
        }

    def _save_grid(self, key: str, grid: Any) -> None:
        """Persist one project-library grid, pickling built grid objects."""

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
        else:
            # Recipe GridSpec or GridRef: already serializable.
            _write_json(self.layout.grid_spec_path(key), grid.to_dict())

    def _save_package(self, key: str, package: PackageSpec) -> None:
        """Persist one project package, pickling it when it carries arrays.

        Source-driven / scalar packages serialize to readable JSON. A package
        whose options embed arrays (a computed ``k`` field, explicit
        stress-period data, etc.) is pickled as an artifact instead, mirroring
        how built grids are persisted, so it can be reused across variants
        without recomputing.
        """

        try:
            doc = package.to_dict()
        except ValueError:
            doc = None
        if doc is not None:
            _write_json(self.layout.package_spec_path(key), doc)
            return
        pickle_path = self.layout.package_pickle_path(key)
        pickle_path.parent.mkdir(parents=True, exist_ok=True)
        with pickle_path.open("wb") as handle:
            pickle.dump(package, handle)
        _write_json(self.layout.package_sidecar_path(key), _artifact_versions())
        _write_json(
            self.layout.package_spec_path(key),
            {
                "kind": "PackageSpec",
                "name": package.name,
                "method": "pickle",
                "pickle_path": pickle_path.relative_to(self.root).as_posix(),
            },
        )

    def save(self) -> Path:
        """Persist the project document, its simulations, and its package library.

        Persistence is opt-in. Returns the path to the written project spec.
        """

        issues = self.validate()
        if issues:
            raise ValueError(
                "Project cannot be saved:\n"
                + "\n".join(f"- {issue}" for issue in issues)
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
        _write_json(self.layout.project_spec_path, self._project_dict())
        for name, simulation in self.simulations.items():
            _write_json(self.layout.simulation_spec_path(name), simulation.to_dict())
        for key, package in self.packages.items():
            self._save_package(key, package)
        for key, grid in self.grids.items():
            self._save_grid(key, grid)
        return self.layout.project_spec_path

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

    def _load_grid(self, key: str, doc: dict[str, Any]) -> Any:
        """Reconstruct one project-library grid from its document."""

        if doc.get("kind") == "GridRef":
            return GridRef.from_dict(doc)
        if doc.get("method") == "pickle":
            self._warn_version_mismatch("grid", key, self.layout.grid_sidecar_path(key))
        return GridSpec.from_dict(doc)

    def _load_package(self, key: str, doc: dict[str, Any]) -> PackageSpec:
        """Reconstruct one project package from its document."""

        if doc.get("method") == "pickle":
            self._warn_version_mismatch(
                "package", key, self.layout.package_sidecar_path(key)
            )
            path = self.root / doc["pickle_path"]
            with path.open("rb") as handle:
                return pickle.load(handle)
        return PackageSpec.from_dict(doc)

    def _warn_version_mismatch(self, kind: str, key: str, sidecar_path: Path) -> None:
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
                stacklevel=3,
            )
