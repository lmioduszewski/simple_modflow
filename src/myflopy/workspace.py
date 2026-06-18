"""Project and run lifecycle objects for spec-driven MODFLOW workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any

import flopy

from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.project.run_model import load_mf6_run, patch_simulation_plot
from myflopy.specs import (
    BuiltSimulation,
    PackageRef,
    PackageSpec,
    SimulationSpec,
    SpecBuildContext,
    _grid_entry,
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

    def __init__(self, *, run: "Run", model_name: str):
        self._initialize_from_built_run(run, model_name)


def load_run(
    workspace: Path | str,
    *,
    executable: str = "mf6",
    load: bool = True,
) -> Run:
    """Load a run workspace and return the standard run object."""

    return Run.load(workspace, executable=executable, load=load)


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

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root))

    @classmethod
    def for_project(cls, name: str, root: Path | str | None = None) -> ProjectLayout:
        """Return the default layout for a project name (defaults under ~/mf6)."""

        project_root = Path.home() / "mf6" / _slug(name) if root is None else Path(root)
        return cls(project_root)

    @property
    def manifest_path(self) -> Path:
        return self.root / PROJECT_MANIFEST_NAME

    @property
    def specs_dir(self) -> Path:
        return self.root / self.specs_dir_name

    @property
    def project_spec_path(self) -> Path:
        return self.specs_dir / PROJECT_SPEC_NAME

    @property
    def simulation_specs_dir(self) -> Path:
        return self.specs_dir / self.simulations_dir_name

    @property
    def package_specs_dir(self) -> Path:
        return self.specs_dir / self.packages_dir_name

    @property
    def inputs_dir(self) -> Path:
        return self.root / self.inputs_dir_name

    def simulation_spec_path(self, name: str) -> Path:
        """Return the durable spec path for one simulation name."""

        return self.simulation_specs_dir / f"{_slug(name)}.json"

    def package_spec_path(self, key: str) -> Path:
        """Return the durable spec path for one project package key."""

        *parents, filename = _package_key_parts(key)
        return self.package_specs_dir.joinpath(*parents, f"{filename}.json")

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
        return {
            "root": str(self.root.as_posix()),
            "specs_dir": self.specs_dir_name,
            "simulations_dir": self.simulations_dir_name,
            "inputs_dir": self.inputs_dir_name,
            "packages_dir": self.packages_dir_name,
        }


def _spec_summary(spec: SimulationSpec) -> dict[str, Any]:
    """Return stable, readable provenance for a simulation specification."""

    def package_name(package: PackageSpec | PackageRef) -> str:
        return package.name

    def package_metadata(package: PackageSpec | PackageRef) -> dict[str, Any]:
        return package.metadata if isinstance(package, PackageSpec) else {}

    return {
        "name": spec.name,
        "models": [
            {
                "name": model.name,
                "type": model.model_type.value,
                "packages": [
                    package_name(package)
                    for package in model.packages
                    if package.enabled
                ],
                "package_metadata": {
                    package_name(package): metadata
                    for package in model.packages
                    if package.enabled and (metadata := package_metadata(package))
                },
                "context_metadata": model.context.metadata,
                "post_build_hooks": [hook.name for hook in model.hooks],
            }
            for model in spec.models
        ],
        "simulation_packages": [
            package_name(package) for package in spec.packages if package.enabled
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
        self.workspace = Path(self.workspace)

    @property
    def manifest_path(self) -> Path:
        """Path to this run's lifecycle manifest."""

        return self.workspace / RUN_MANIFEST_NAME

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
        self.simulation.write_simulation(silent=True)
        self._record_status("written")
        return self.workspace

    def execute(self, *, silent: bool = True) -> tuple[bool, list[str]]:
        """Write and execute the simulation with FloPy 3.10."""

        if (
            self.status not in {"written", "completed", "failed"}
            or self.simulation is None
        ):
            self.write()
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
        methods. Reopened/file-backed GWF models return ``LoadedMf6Run``.

        Use :meth:`flopy_model` when you need a raw FloPy model object or are
        working with a non-GWF model type.
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
    """Workspace and in-memory specification registry for related runs.

    Run manifests and MF6 files are durable. The Python specification registry
    is intentionally in-memory because its builders may be arbitrary callables.
    """

    def __init__(
        self,
        root: Path | str,
        *,
        name: str | None = None,
        layout: ProjectLayout | None = None,
    ):
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

        The key is how models reference the package later, for example
        ``mf.ref("npf/base")``. Packages live independently of any model and
        are resolved into a simulation only when a run is built.
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
            grid_specs={key: _grid_entry(value) for key, value in self.grids.items()},
        )

    def prepare_run(
        self,
        name: str,
        simulation: SimulationSpec | str,
        *,
        executable: str = "mf6",
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
    ) -> Run:
        """Create an unbuilt run from a simulation spec or registered spec name."""

        spec = (
            self.simulations[simulation] if isinstance(simulation, str) else simulation
        )
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

    def validate(self) -> list[str]:
        """Return human-readable persistence issues without touching disk.

        Persistence requires serializable specs. Simulations that embed raw
        arrays, or that contain exchanges, are reported here rather than failing
        silently at :meth:`save`.
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
        return issues

    def _project_dict(self) -> dict[str, Any]:
        """Return the durable project document."""

        return {
            "kind": "Project",
            "name": self.name,
            "layout": self.layout.to_dict(),
            "simulations": sorted(self.simulations),
            "packages": sorted(self.packages),
        }

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
            _write_json(self.layout.package_spec_path(key), package.to_dict())
        return self.layout.project_spec_path

    @classmethod
    def load(cls, root: Path | str, *, layout: ProjectLayout | None = None) -> Project:
        """Load a project previously written with :meth:`save`."""

        project = cls(root, layout=layout)
        spec_path = project.layout.project_spec_path
        if not spec_path.exists():
            raise FileNotFoundError(f"Project spec not found: {spec_path}")
        data = _read_json(spec_path)
        project.name = data.get("name", project.name)
        for key in data.get("packages", ()):
            package_path = project.layout.package_spec_path(key)
            if package_path.exists():
                project.packages[key] = PackageSpec.from_dict(_read_json(package_path))
        for name in data.get("simulations", ()):
            sim_path = project.layout.simulation_spec_path(name)
            if sim_path.exists():
                project.add_simulation(SimulationSpec.from_dict(_read_json(sim_path)))
        return project
