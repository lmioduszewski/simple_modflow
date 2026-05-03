"""Dataclasses for project-level model specs, run definitions, and run records.

These objects are intentionally lightweight. They hold stable identifiers,
workspace paths, reusable package references, and small bits of metadata so the
project layer can organize MF6 runs without depending on large pickled model
objects.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def utc_now_iso() -> str:
    """Return the current UTC timestamp in ISO-8601 format."""

    return datetime.now(timezone.utc).isoformat()


def default_run_artifact_paths(run_id: str) -> dict[str, str]:
    """Return the default relative file layout for a run workspace.

    Parameters
    ----------
    run_id
        Run identifier that will also be used as the default MF6 model name.
    """

    return {
        "workspace": ".",
        "run_manifest": "run.toml",
        "simulation_name_file": "mfsim.nam",
        "listing_file": "mfsim.lst",
        "model_name_file": f"{run_id}.nam",
        "model_object_file": f"{run_id}.model",
        "heads_file": f"{run_id}.hds",
        "budget_file": f"{run_id}.cbc",
    }


@dataclass(slots=True)
class ModelSpec:
    """Conceptual base-model definition shared by many related runs.

    Parameters
    ----------
    name
        Stable model-family identifier.
    description
        Optional human-readable description of the model family.
    grid_ref
        Optional identifier for the grid family or source grid asset.
    tags
        Optional grouping labels.
    default_packages
        Optional package/artifact ids that are typical defaults for this model
        family.
    default_regions
        Optional region/group names that conceptually belong to this model family.
    metadata
        Free-form extra metadata.
    """

    name: str
    description: str | None = None
    grid_ref: str | None = None
    tags: list[str] = field(default_factory=list)
    default_packages: dict[str, str] = field(default_factory=dict)
    default_regions: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class RunSpec:
    """Planned run definition used before or during model construction.

    Parameters
    ----------
    run_id
        Unique run identifier within a project catalog.
    model_spec
        Name of the parent :class:`ModelSpec`.
    scenario
        Optional scenario label such as ``baseline`` or ``high_recharge``.
    grid_ref
        Optional grid identifier when it differs from the parent model spec.
    package_versions
        Mapping from package keys like ``"chd"`` or ``"sfr"`` to registered
        reusable artifact ids.
    workspace
        Optional explicit workspace path. Relative paths are resolved against the
        project root.
    notes
        Optional notes for the run manifest.
    tags
        Optional grouping labels.
    parameter_overrides
        Optional run-specific numeric/logical overrides recorded for provenance.
    metadata
        Additional free-form metadata.
    """

    run_id: str
    model_spec: str
    scenario: str | None = None
    grid_ref: str | None = None
    package_versions: dict[str, str] = field(default_factory=dict)
    workspace: Path | None = None
    notes: str | None = None
    tags: list[str] = field(default_factory=list)
    parameter_overrides: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def resolve_workspace(self, project_root: Path, runs_dir_name: str = "runs") -> Path:
        """Resolve the run workspace from the project root and optional override."""

        if self.workspace is None:
            return project_root / runs_dir_name / self.run_id
        if self.workspace.is_absolute():
            return self.workspace
        return project_root / self.workspace

    def to_run_record(
        self,
        *,
        project_root: Path,
        runs_dir_name: str = "runs",
        status: str = "registered",
    ) -> RunRecord:
        """Materialize this planned run definition into a concrete run record.

        Parameters
        ----------
        project_root
            Root folder of the containing :class:`~simple_modflow.project.catalog.ProjectCatalog`.
        runs_dir_name
            Default subdirectory used when ``workspace`` is not set explicitly.
        status
            Initial run status to record in the manifest.
        """

        workspace = self.resolve_workspace(project_root, runs_dir_name=runs_dir_name)
        return RunRecord(
            run_id=self.run_id,
            model_spec=self.model_spec,
            workspace=workspace,
            scenario=self.scenario,
            grid_ref=self.grid_ref,
            package_versions=dict(self.package_versions),
            status=status,
            notes=self.notes,
            tags=list(self.tags),
            parameter_overrides=dict(self.parameter_overrides),
            metadata=dict(self.metadata),
            paths=default_run_artifact_paths(self.run_id),
        )


@dataclass(slots=True)
class RunRecord:
    """Concrete run manifest record stored alongside or within a run workspace.

    This object is what ties together:
    - the run id and model spec
    - the on-disk workspace
    - reusable package references
    - lightweight provenance metadata
    """

    run_id: str
    model_spec: str
    workspace: Path
    scenario: str | None = None
    grid_ref: str | None = None
    package_versions: dict[str, str] = field(default_factory=dict)
    status: str = "registered"
    created_at: str = field(default_factory=utc_now_iso)
    completed_at: str | None = None
    notes: str | None = None
    tags: list[str] = field(default_factory=list)
    parameter_overrides: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    paths: dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        self.workspace = Path(self.workspace)
        if not self.paths:
            self.paths = default_run_artifact_paths(self.run_id)

    @property
    def manifest_path(self) -> Path:
        """Path to the TOML manifest that represents this run record."""

        return self.workspace / self.paths.get("run_manifest", "run.toml")

    @property
    def workspace_parent(self) -> Path:
        """Parent folder passed to ``SimulationBase`` as ``mf_folder_path``."""

        return self.workspace.parent

    def simulation_kwargs(self) -> dict[str, Path | str]:
        """Return constructor kwargs for a catalog-managed ``SimulationBase``.

        The returned mapping aligns the MF6 model name and workspace layout with
        this run record.
        """

        return {
            "name": self.run_id,
            "mf_folder_path": self.workspace_parent,
        }

    def get_path(self, key: str) -> Path:
        """Resolve a named path entry from ``self.paths`` to an absolute path."""

        if key == "workspace":
            return self.workspace
        relative = self.paths[key]
        return self.workspace / relative
