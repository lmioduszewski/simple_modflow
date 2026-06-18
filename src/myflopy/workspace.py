"""Project and run lifecycle objects for spec-driven MODFLOW workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import flopy

from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.project.run_model import load_mf6_run, patch_simulation_plot
from myflopy.specs import (
    BuiltSimulation,
    ModelSpec,
    PackageRef,
    PackageSpec,
    SimulationSpec,
    SpecBuildContext,
)


RUN_MANIFEST_NAME = "run.json"
PROJECT_MANIFEST_NAME = "project.json"


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
        self, root: Path | str, *, name: str | None = None, create: bool = True
    ):
        self.root = Path(root)
        self.name = name or self.root.name
        self.packages: dict[str, PackageSpec] = {}
        self.models: dict[str, ModelSpec] = {}
        self.simulations: dict[str, SimulationSpec] = {}

        if create:
            self.runs_dir.mkdir(parents=True, exist_ok=True)
            self.save_manifest()
        elif self.manifest_path.exists():
            data = json.loads(self.manifest_path.read_text(encoding="utf-8"))
            self.name = data.get("project", {}).get("name", self.name)
        else:
            raise FileNotFoundError(f"Project manifest not found: {self.manifest_path}")

    @property
    def runs_dir(self) -> Path:
        """Directory containing project-owned run workspaces."""

        return self.root / "runs"

    @property
    def manifest_path(self) -> Path:
        """Path to the project manifest."""

        return self.root / PROJECT_MANIFEST_NAME

    def save_manifest(self) -> Path:
        """Persist the small project workspace manifest."""

        _write_json(
            self.manifest_path, {"project": {"name": self.name, "runs_dir": "runs"}}
        )
        return self.manifest_path

    def add_package(self, package: PackageSpec) -> PackageSpec:
        """Add or replace a reusable package specification."""

        self.packages[package.name] = package
        return package

    def add_model(self, model: ModelSpec) -> ModelSpec:
        """Add or replace a reusable model specification."""

        self.models[model.name] = model
        return model

    def add_simulation(self, simulation: SimulationSpec) -> SimulationSpec:
        """Add or replace a reusable simulation specification."""

        self.simulations[simulation.name] = simulation
        return simulation

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

    @classmethod
    def reopen(cls, root: Path | str) -> Project:
        """Reopen an existing project workspace."""

        return cls(root, create=False)
