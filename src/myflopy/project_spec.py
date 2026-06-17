"""Durable project specifications for related MODFLOW 6 simulations."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import json
from pathlib import Path
import re
from typing import Any

from myflopy.specs import SimulationSpec
from myflopy.workspace import Run


PROJECT_MANIFEST_NAME = "project.json"
PROJECT_SPEC_NAME = "project_spec.json"


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._-")
    if not slug:
        raise ValueError(
            "Project and simulation names must contain at least one safe character."
        )
    return slug


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


@dataclass(frozen=True, slots=True)
class ProjectLayout:
    """Directory layout for a durable project specification."""

    root: Path
    specs_dir_name: str = "specs"
    simulations_dir_name: str = "simulations"
    inputs_dir_name: str = "inputs"
    packages_dir_name: str = "packages"

    def __post_init__(self) -> None:
        object.__setattr__(self, "root", Path(self.root))

    @classmethod
    def for_project(cls, name: str, root: Path | str | None = None) -> ProjectLayout:
        """Return the default layout for a project name."""

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

    @property
    def simulations_dir(self) -> Path:
        return self.root / self.simulations_dir_name

    def simulation_workspace(self, name: str) -> Path:
        """Return the materialized MF6 workspace for one simulation name."""

        return self.simulations_dir / _slug(name)

    def simulation_spec_path(self, name: str) -> Path:
        """Return the durable spec path for one simulation name."""

        return self.simulation_specs_dir / f"{_slug(name)}.json"

    def ensure(self) -> None:
        """Create the default project directories, excluding optional snapshots."""

        for path in (
            self.root,
            self.specs_dir,
            self.simulation_specs_dir,
            self.package_specs_dir,
            self.inputs_dir,
            self.simulations_dir,
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


@dataclass(frozen=True, slots=True)
class ProjectSpec:
    """Durable specification for a named collection of MF6 simulations."""

    name: str
    root: Path | str | None = None
    simulations: tuple[SimulationSpec, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "simulations", tuple(self.simulations))
        names = [simulation.name for simulation in self.simulations]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError(f"Duplicate simulation names: {', '.join(duplicates)}")

    @property
    def layout(self) -> ProjectLayout:
        """Return the resolved project directory layout."""

        return ProjectLayout.for_project(self.name, self.root)

    def with_simulation(self, simulation: SimulationSpec) -> ProjectSpec:
        """Return a copy with ``simulation`` added or replaced by name."""

        simulations = list(self.simulations)
        for index, existing in enumerate(simulations):
            if existing.name == simulation.name:
                simulations[index] = simulation
                break
        else:
            simulations.append(simulation)
        return replace(self, simulations=tuple(simulations))

    def with_simulations(self, *simulations: SimulationSpec) -> ProjectSpec:
        """Return a copy with each supplied simulation added or replaced."""

        spec = self
        for simulation in simulations:
            spec = spec.with_simulation(simulation)
        return spec

    def simulation(self, name: str) -> SimulationSpec:
        """Return one named simulation specification."""

        for simulation in self.simulations:
            if simulation.name == name:
                return simulation
        raise KeyError(f"Project '{self.name}' has no simulation named '{name}'.")

    def validate(self) -> list[str]:
        """Return human-readable validation issues without touching disk."""

        issues: list[str] = []
        names = {simulation.name for simulation in self.simulations}
        for simulation in self.simulations:
            if (
                simulation.derived_from is not None
                and simulation.derived_from not in names
            ):
                issues.append(
                    f"Simulation '{simulation.name}' derives from unknown "
                    f"simulation '{simulation.derived_from}'."
                )
            try:
                simulation.to_dict()
            except ValueError as error:
                issues.append(
                    f"Simulation '{simulation.name}' is not serializable: {error}"
                )
        return issues

    def to_dict(self) -> dict[str, Any]:
        """Return the durable project spec manifest."""

        return {
            "kind": "ProjectSpec",
            "name": self.name,
            "layout": self.layout.to_dict(),
            "simulations": [simulation.name for simulation in self.simulations],
            "metadata": self.metadata,
        }

    def save(self) -> Path:
        """Persist this project spec and each named simulation spec."""

        issues = self.validate()
        if issues:
            raise ValueError(
                "ProjectSpec is not valid:\n"
                + "\n".join(f"- {issue}" for issue in issues)
            )

        layout = self.layout
        layout.ensure()
        _write_json(
            layout.manifest_path,
            {
                "project": {
                    "name": self.name,
                    "spec": f"{layout.specs_dir_name}/{PROJECT_SPEC_NAME}",
                    "simulations_dir": layout.simulations_dir_name,
                }
            },
        )
        _write_json(layout.project_spec_path, self.to_dict())
        for simulation in self.simulations:
            _write_json(
                layout.simulation_spec_path(simulation.name), simulation.to_dict()
            )
        return layout.project_spec_path

    @classmethod
    def load(cls, root: Path | str) -> ProjectSpec:
        """Load a project spec from a project root directory."""

        root = Path(root)
        project_spec_path = root / "specs" / PROJECT_SPEC_NAME
        if not project_spec_path.exists():
            manifest_path = root / PROJECT_MANIFEST_NAME
            if not manifest_path.exists():
                raise FileNotFoundError(f"Project spec not found under: {root}")
            manifest = _read_json(manifest_path)
            project_spec_path = root / manifest["project"].get(
                "spec",
                f"specs/{PROJECT_SPEC_NAME}",
            )
        data = _read_json(project_spec_path)
        simulations = tuple(
            SimulationSpec.from_dict(
                _read_json(root / "specs" / "simulations" / f"{_slug(name)}.json")
            )
            for name in data.get("simulations", ())
        )
        return cls(
            name=data["name"],
            root=root,
            simulations=simulations,
            metadata=dict(data.get("metadata", {})),
        )

    def _simulation_for_materialization(self, name: str) -> SimulationSpec:
        simulation = self.simulation(name)
        models = []
        for model in simulation.models:
            if model.options.get("model_rel_path") in {None, "."}:
                models.append(model.with_options(model_rel_path=model.name))
            else:
                models.append(model)
        return simulation.with_models(*models).with_workspace(
            self.layout.simulation_workspace(name)
        )

    def build(self, name: str, *, overwrite: bool = False) -> Run:
        """Build a named simulation into ``simulations/<name>/``."""

        workspace = self.layout.simulation_workspace(name)
        if workspace.exists() and any(workspace.iterdir()) and not overwrite:
            raise FileExistsError(f"Simulation workspace already exists: {workspace}")
        run = Run(
            name=name,
            workspace=workspace,
            spec=self._simulation_for_materialization(name),
            executable=self.simulation(name).executable,
        )
        run.build()
        return run

    def run(
        self,
        name: str,
        *,
        overwrite: bool = False,
        silent: bool = True,
    ) -> Run:
        """Build, write, and execute a named simulation."""

        run = self.build(name, overwrite=overwrite)
        run.execute(silent=silent)
        return run


__all__ = [
    "ProjectLayout",
    "ProjectSpec",
]
