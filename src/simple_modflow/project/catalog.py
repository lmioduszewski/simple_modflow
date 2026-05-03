"""Project-level registry for model specs, runs, and reusable package artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from simple_modflow.project.components import (
    PackageArtifact,
    PackageCompatibilityError,
    apply_package_artifact,
    build_package_artifact,
    derive_package_artifact,
    format_validation_report,
    package_artifact_apply_order,
    validate_package_artifact_reference,
)
from simple_modflow.project.compare import RunComparison
from simple_modflow.project.discovery import DiscoveredRun, discover_existing_runs
from simple_modflow.project.loaders import RunLoader
from simple_modflow.project.manifest_io import (
    PROJECT_MANIFEST_NAME,
    RUN_MANIFEST_NAME,
    load_package_artifact,
    load_model_spec,
    load_run_record,
    read_toml,
    save_package_artifact,
    save_model_spec,
    save_run_record,
    write_toml,
)
from simple_modflow.project.specs import ModelSpec, RunRecord, RunSpec


class ProjectCatalog:
    """Manage the on-disk structure and metadata for a modeling project."""

    def __init__(
        self,
        root: Path,
        *,
        name: str | None = None,
        runs_dir_name: str = "runs",
        model_specs_dir_name: str = "model_specs",
        run_records_dir_name: str = "run_records",
        package_artifacts_dir_name: str = "package_artifacts",
        create: bool = True,
    ):
        """Create or reopen a project catalog rooted at ``root``."""

        self.root = Path(root)
        self.name = name or self.root.name
        self.runs_dir_name = runs_dir_name
        self.model_specs_dir_name = model_specs_dir_name
        self.run_records_dir_name = run_records_dir_name
        self.package_artifacts_dir_name = package_artifacts_dir_name

        if create:
            self.ensure_structure()
            self.write_manifest()
        elif self.project_manifest_path.exists():
            self._load_manifest()
        self.loader = RunLoader(self)
        self.compare = RunComparison(self)

    @property
    def project_manifest_path(self) -> Path:
        """Path to the catalog's top-level project manifest."""

        return self.root / PROJECT_MANIFEST_NAME

    @property
    def runs_dir(self) -> Path:
        """Directory where catalog-managed run workspaces live by default."""

        return self.root / self.runs_dir_name

    @property
    def model_specs_dir(self) -> Path:
        """Directory holding model-spec manifest files."""

        return self.root / self.model_specs_dir_name

    @property
    def run_records_dir(self) -> Path:
        """Directory holding external/imported run manifests."""

        return self.root / self.run_records_dir_name

    @property
    def package_artifacts_dir(self) -> Path:
        """Directory holding reusable package artifact manifests."""

        return self.root / self.package_artifacts_dir_name

    def ensure_structure(self):
        """Create the expected folder structure for the catalog."""

        self.root.mkdir(parents=True, exist_ok=True)
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.model_specs_dir.mkdir(parents=True, exist_ok=True)
        self.run_records_dir.mkdir(parents=True, exist_ok=True)
        self.package_artifacts_dir.mkdir(parents=True, exist_ok=True)

    def _load_manifest(self):
        """Load catalog settings back from the top-level project manifest."""

        data = read_toml(self.project_manifest_path).get("project", {})
        self.name = data.get("name", self.name)
        self.runs_dir_name = data.get("runs_dir", self.runs_dir_name)
        self.model_specs_dir_name = data.get("model_specs_dir", self.model_specs_dir_name)
        self.run_records_dir_name = data.get("run_records_dir", self.run_records_dir_name)
        self.package_artifacts_dir_name = data.get("package_artifacts_dir", self.package_artifacts_dir_name)

    def write_manifest(self):
        """Write the top-level project manifest."""

        write_toml(
            {
                "project": {
                    "name": self.name,
                    "runs_dir": self.runs_dir_name,
                    "model_specs_dir": self.model_specs_dir_name,
                    "run_records_dir": self.run_records_dir_name,
                    "package_artifacts_dir": self.package_artifacts_dir_name,
                }
            },
            self.project_manifest_path,
        )

    def register_model_spec(self, spec: ModelSpec, *, overwrite: bool = False) -> Path:
        """Register a conceptual model spec in the catalog."""

        path = self.model_specs_dir / f"{spec.name}.toml"
        if path.exists() and not overwrite:
            raise FileExistsError(f"Model spec already exists: {spec.name}")
        save_model_spec(spec, path)
        return path

    def load_model_spec(self, name: str) -> ModelSpec:
        """Load one model spec by name."""

        path = self.model_specs_dir / f"{name}.toml"
        if not path.exists():
            raise FileNotFoundError(f"Model spec not found: {name}")
        return load_model_spec(path)

    def list_model_specs(self) -> pd.DataFrame:
        """Return a tabular summary of registered model specs."""

        rows = []
        for path in sorted(self.model_specs_dir.glob("*.toml")):
            spec = load_model_spec(path)
            rows.append(
                {
                    "name": spec.name,
                    "description": spec.description,
                    "grid_ref": spec.grid_ref,
                    "tags": spec.tags,
                }
            )
        return pd.DataFrame(rows)

    def create_run(self, run_spec: RunSpec, *, overwrite: bool = False) -> RunRecord:
        """Pre-register a run definition and create its workspace on disk.

        This is mainly for advanced workflows:
        - planning runs before building a model
        - creating manifests for runs that will be built later
        - tests, imports, or scripted batch setup

        For the normal interactive build flow, you can usually skip this and create the
        run automatically via ``SimulationBase(..., project_catalog=..., run_id=...,
        model_spec=..., package_versions=...)``.
        """
        record = run_spec.to_run_record(
            project_root=self.root,
            runs_dir_name=self.runs_dir_name,
            status="registered",
        )
        if record.workspace.exists() and any(record.workspace.iterdir()) and not overwrite:
            raise FileExistsError(f"Run workspace already exists: {record.workspace}")
        record.workspace.mkdir(parents=True, exist_ok=True)
        save_run_record(record)
        return record

    def register_run_record(
        self,
        record: RunRecord,
        *,
        overwrite: bool = False,
        manifest_path: Path | None = None,
    ) -> Path:
        """Persist an existing run record to a manifest path."""

        manifest_path = record.manifest_path if manifest_path is None else manifest_path
        if manifest_path.exists() and not overwrite:
            raise FileExistsError(f"Run manifest already exists: {manifest_path}")
        if manifest_path == record.manifest_path:
            record.workspace.mkdir(parents=True, exist_ok=True)
        else:
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
        save_run_record(record, path=manifest_path)
        return manifest_path

    def register_external_run(self, record: RunRecord, *, overwrite: bool = False) -> Path:
        """Register a run manifest for a workspace that lives outside ``runs/``."""

        manifest_path = self.run_records_dir / f"{record.run_id}.toml"
        if manifest_path.exists() and not overwrite:
            raise FileExistsError(f"Run manifest already exists: {manifest_path}")
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        save_run_record(record, path=manifest_path)
        return manifest_path

    def discover_existing_runs(self, search_root: Path) -> list[DiscoveredRun]:
        """Search ``search_root`` for MF6 runs that can be imported or explored."""

        return discover_existing_runs(search_root)

    def import_existing_runs(
        self,
        search_root: Path,
        *,
        model_spec: str = "imported_model",
        overwrite: bool = False,
        skip_existing: bool = True,
    ) -> list[RunRecord]:
        """Import an existing directory tree of MF6 workspaces into this catalog."""

        imported = []
        for discovered in self.discover_existing_runs(search_root):
            record = discovered.to_run_record(model_spec=model_spec)
            manifest_path = self.run_records_dir / f"{record.run_id}.toml"
            if manifest_path.exists():
                if overwrite:
                    self.register_external_run(record, overwrite=True)
                    imported.append(record)
                    continue
                if skip_existing:
                    imported.append(load_run_record(manifest_path))
                    continue
                raise FileExistsError(f"Run manifest already exists: {manifest_path}")
            self.register_external_run(record, overwrite=False)
            imported.append(record)
        return imported

    def discover_run_manifests(self) -> list[Path]:
        """Return all run-manifest paths currently discoverable in this catalog."""

        manifests = []
        if self.runs_dir.exists():
            manifests.extend(self.runs_dir.glob(f"*/{RUN_MANIFEST_NAME}"))
        if self.run_records_dir.exists():
            manifests.extend(self.run_records_dir.glob("*.toml"))
        return sorted(Path(path) for path in manifests)

    def discover_runs(self) -> list[RunRecord]:
        """Load every currently discoverable run record."""

        return [load_run_record(path) for path in self.discover_run_manifests()]

    def load_run(self, run_id: str) -> RunRecord:
        """Load one run record by id."""

        manifest_paths = [
            self.runs_dir / run_id / RUN_MANIFEST_NAME,
            self.run_records_dir / f"{run_id}.toml",
        ]
        for manifest_path in manifest_paths:
            if manifest_path.exists():
                return load_run_record(manifest_path)
        raise FileNotFoundError(f"Run not found: {run_id}")

    def load_run_model(
        self,
        run: str | RunRecord,
        *,
        crs: str = "EPSG:2927",
        verbosity_level: int = 0,
    ):
        """Open a catalog-managed run as a file-backed ``LoadedMf6Run``."""

        return self.loader.load_run_model(
            run,
            crs=crs,
            verbosity_level=verbosity_level,
        )

    def attach_run(self, model, run: str | RunRecord):
        """Associate an in-memory model with a catalog-managed run definition."""

        record = run if isinstance(run, RunRecord) else self.load_run(run)
        model.project_catalog = self
        model.run_record = record
        model.run_id = record.run_id
        return record

    def create_package_artifact(
        self,
        artifact_id: str,
        *,
        model,
        package_name: str,
        description: str | None = None,
        tags: list[str] | None = None,
        metadata: dict | None = None,
        overwrite: bool = False,
    ) -> PackageArtifact:
        """Capture a reusable package artifact from an existing model package."""

        artifact = build_package_artifact(
            model,
            artifact_id=artifact_id,
            package_name=package_name,
            description=description,
            tags=tags,
            metadata=metadata,
        )
        self.register_package_artifact(artifact, overwrite=overwrite)
        return artifact

    def register_package_artifact(self, artifact: PackageArtifact, *, overwrite: bool = False) -> Path:
        """Persist a reusable package artifact manifest in this catalog."""

        path = self.package_artifacts_dir / f"{artifact.artifact_id}.toml"
        if path.exists() and not overwrite:
            raise FileExistsError(f"Package artifact already exists: {artifact.artifact_id}")
        save_package_artifact(artifact, path)
        return path

    def load_package_artifact(self, artifact_id: str) -> PackageArtifact:
        """Load one reusable package artifact by id."""

        path = self.package_artifacts_dir / f"{artifact_id}.toml"
        if not path.exists():
            raise FileNotFoundError(f"Package artifact not found: {artifact_id}")
        return load_package_artifact(path)

    def list_package_artifacts(self) -> pd.DataFrame:
        """Return a tabular summary of registered reusable package artifacts."""

        rows = []
        for path in sorted(self.package_artifacts_dir.glob("*.toml")):
            artifact = load_package_artifact(path)
            rows.append(
                {
                    "artifact_id": artifact.artifact_id,
                    "package_type": artifact.package_type,
                    "source_run_id": artifact.source_run_id,
                    "source_workspace": str(artifact.source_workspace),
                    "source_package_name": artifact.source_package_name,
                    "storage_mode": artifact.storage_mode,
                    "derived_from_artifact_id": artifact.derived_from_artifact_id,
                    "derived_from_run_id": artifact.derived_from_run_id,
                    "tags": artifact.tags,
                }
            )
        return pd.DataFrame(rows)

    def derive_package_artifact(
        self,
        source_artifact_id: str,
        *,
        artifact_id: str,
        description: str | None = None,
        tags: list[str] | None = None,
        metadata_updates: dict | None = None,
        package_data_updates: dict | None = None,
        overwrite: bool = False,
    ) -> PackageArtifact:
        """Create a lightly modified snapshot artifact with recorded lineage.

        This is the current safe way to expand reusable artifacts: derive a new
        same-grid snapshot from an earlier artifact while recording lineage and
        any targeted data/metadata edits.
        """

        source_artifact = self.load_package_artifact(source_artifact_id)
        artifact = derive_package_artifact(
            source_artifact,
            artifact_id=artifact_id,
            description=description,
            tags=tags,
            metadata_updates=metadata_updates,
            package_data_updates=package_data_updates,
        )
        self.register_package_artifact(artifact, overwrite=overwrite)
        return artifact

    def _resolve_run_artifacts(
        self,
        run: str | RunRecord,
        *,
        package_keys: list[str] | None = None,
    ) -> list[tuple[str, str, PackageArtifact]]:
        """Resolve a run's package references into loaded artifact objects."""

        record = run if isinstance(run, RunRecord) else self.load_run(run)
        items = list(record.package_versions.items())
        if package_keys is not None:
            keys = set(package_keys)
            items = [(key, value) for key, value in items if key in keys]
        resolved = []
        for key, artifact_id in items:
            resolved.append((key, artifact_id, self.load_package_artifact(artifact_id)))
        return resolved

    def validate_run_package_artifacts(
        self,
        model,
        run: str | RunRecord,
        *,
        package_keys: list[str] | None = None,
        raise_on_error: bool = False,
    ) -> pd.DataFrame:
        """Validate a run's referenced package artifacts before attaching them.

        Parameters
        ----------
        model
            Target model that will receive the artifacts.
        run
            Run id or run record whose ``package_versions`` should be checked.
        package_keys
            Optional subset of package keys to validate.
        raise_on_error
            When ``True``, raise :class:`PackageCompatibilityError` with a
            readable multi-line report if any artifact is not ready.
        """

        resolved = self._resolve_run_artifacts(run, package_keys=package_keys)
        referenced_package_types = {artifact.package_type for _, _, artifact in resolved}
        available_package_types = {name.lower() for name in model.package_names}
        rows = []
        for package_key, artifact_id, artifact in resolved:
            row = validate_package_artifact_reference(
                artifact,
                model,
                available_package_types=available_package_types,
                referenced_package_types=referenced_package_types - {artifact.package_type},
            )
            row["package_key"] = package_key
            row["artifact_id"] = artifact_id
            rows.append(row)

        frame = pd.DataFrame(rows)
        if raise_on_error and not frame.empty and not bool(frame["valid"].all()):
            raise PackageCompatibilityError(
                "Referenced package artifacts are not ready to apply:\n"
                + format_validation_report(rows)
            )
        return frame

    def apply_package_artifact(
        self,
        model,
        artifact_id: str,
        *,
        validate: bool = True,
    ):
        """Attach a stored package artifact to a compatible target model."""

        artifact = self.load_package_artifact(artifact_id)
        return apply_package_artifact(model, artifact, validate=validate)

    def apply_run_package_artifacts(
        self,
        model,
        run: str | RunRecord,
        *,
        package_keys: list[str] | None = None,
        validate: bool = True,
    ) -> dict[str, object]:
        """Resolve a run's package artifact ids and attach them to a model.

        This turns `run.package_versions` into an operational input: each value
        is treated as a registered package artifact id and attached if the
        target model is compatible. Artifacts are applied in a dependency-aware
        order so packages like ``mvr`` are attached after the packages they
        depend on.
        """

        record = run if isinstance(run, RunRecord) else self.load_run(run)
        resolved = self._resolve_run_artifacts(record, package_keys=package_keys)
        if validate:
            self.validate_run_package_artifacts(
                model,
                record,
                package_keys=package_keys,
                raise_on_error=True,
            )
        attached = {}
        resolved.sort(key=lambda item: (package_artifact_apply_order(item[2].package_type), item[0]))
        for key, artifact_id, artifact in resolved:
            attached[key] = apply_package_artifact(
                model,
                artifact,
                validate=validate,
            )
        return attached

    def list_runs(self) -> pd.DataFrame:
        """Return a tabular summary of registered/discovered runs."""

        rows = []
        for record in self.discover_runs():
            rows.append(
                {
                    "run_id": record.run_id,
                    "model_spec": record.model_spec,
                    "scenario": record.scenario,
                    "status": record.status,
                    "workspace": str(record.workspace),
                    "tags": record.tags,
                }
            )
        return pd.DataFrame(rows)

    def run_summary(self) -> pd.DataFrame:
        """Return a notebook-friendly run summary with archive-style columns."""

        rows = []
        for record in self.discover_runs():
            rows.append(
                {
                    "run_id": record.run_id,
                    "model_spec": record.model_spec,
                    "status": record.status,
                    "workspace": str(record.workspace),
                    "family": record.metadata.get("family"),
                    "grid_type": record.metadata.get("grid_type"),
                    "has_heads": record.metadata.get("has_heads"),
                    "has_budget": record.metadata.get("has_budget"),
                    "tags": record.tags,
                }
            )
        return pd.DataFrame(rows)
