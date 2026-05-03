"""TOML serialization helpers for project catalogs, run records, and artifacts."""

from __future__ import annotations

from pathlib import Path
import tomllib

from simple_modflow.project.components import PackageArtifact
from simple_modflow.project.specs import ModelSpec, RunRecord


RUN_MANIFEST_NAME = "run.toml"
PROJECT_MANIFEST_NAME = "project.toml"


def _toml_quote(value: str) -> str:
    """Return a TOML-safe quoted string."""

    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _clean_dict(data: dict) -> dict:
    """Drop ``None`` values and empty nested dictionaries from a structure."""

    cleaned = {}
    for key, value in data.items():
        if value is None:
            continue
        if isinstance(value, dict):
            nested = _clean_dict(value)
            if nested:
                cleaned[key] = nested
            continue
        if isinstance(value, list):
            cleaned[key] = list(value)
            continue
        cleaned[key] = value
    return cleaned


def _format_value(value):
    """Format one scalar/list value for the limited TOML writer used here."""

    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, Path):
        return _toml_quote(value.as_posix())
    if isinstance(value, str):
        return _toml_quote(value)
    if isinstance(value, list):
        return "[" + ", ".join(_format_value(item) for item in value) + "]"
    raise TypeError(f"Unsupported TOML value type: {type(value)!r}")


def _write_table(lines: list[str], data: dict, prefix: str = ""):
    """Append one TOML table and any nested child tables to ``lines``."""

    scalars = []
    nested = []
    for key, value in data.items():
        if isinstance(value, dict):
            nested.append((key, value))
        else:
            scalars.append((key, value))

    if prefix:
        lines.append(f"[{prefix}]")
    for key, value in scalars:
        lines.append(f"{key} = {_format_value(value)}")
    if prefix and nested:
        lines.append("")

    for index, (key, value) in enumerate(nested):
        child_prefix = f"{prefix}.{key}" if prefix else key
        _write_table(lines, value, child_prefix)
        if index != len(nested) - 1:
            lines.append("")


def dump_toml(data: dict) -> str:
    """Serialize a nested dictionary into the limited TOML subset used here."""

    cleaned = _clean_dict(data)
    lines: list[str] = []
    _write_table(lines, cleaned)
    return "\n".join(lines).strip() + "\n"


def write_toml(data: dict, path: Path):
    """Write a TOML payload to disk, creating parent folders as needed."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dump_toml(data), encoding="utf-8")


def read_toml(path: Path) -> dict:
    """Read a TOML file into a Python dictionary."""

    with path.open("rb") as file:
        return tomllib.load(file)


def _serialize_path(path: Path, *, base_dir: Path) -> str:
    """Serialize a path relative to ``base_dir`` where possible."""

    if path == base_dir:
        return "."
    try:
        return path.relative_to(base_dir).as_posix()
    except ValueError:
        return path.as_posix()


def _deserialize_path(raw: str, *, base_dir: Path) -> Path:
    """Restore a serialized relative/absolute path back to a ``Path`` object."""

    path = Path(raw)
    if path.is_absolute():
        return path
    if raw == ".":
        return base_dir
    return base_dir / path


def model_spec_to_dict(spec: ModelSpec) -> dict:
    """Convert a :class:`ModelSpec` into a TOML-ready dictionary."""

    return {
        "model_spec": {
            "name": spec.name,
            "description": spec.description,
            "grid_ref": spec.grid_ref,
            "tags": spec.tags,
        },
        "packages": dict(spec.default_packages),
        "regions": {
            "default_regions": list(spec.default_regions),
        },
        "metadata": dict(spec.metadata),
    }


def save_model_spec(spec: ModelSpec, path: Path):
    """Write a :class:`ModelSpec` manifest to disk."""

    write_toml(model_spec_to_dict(spec), path)


def load_model_spec(path: Path) -> ModelSpec:
    """Load a :class:`ModelSpec` from a TOML manifest path."""

    data = read_toml(path)
    spec_data = data.get("model_spec", {})
    package_data = data.get("packages", {})
    region_data = data.get("regions", {})
    metadata = data.get("metadata", {})
    return ModelSpec(
        name=spec_data["name"],
        description=spec_data.get("description"),
        grid_ref=spec_data.get("grid_ref"),
        tags=list(spec_data.get("tags", [])),
        default_packages=dict(package_data),
        default_regions=list(region_data.get("default_regions", [])),
        metadata=dict(metadata),
    )


def run_record_to_dict(record: RunRecord, *, manifest_dir: Path | None = None) -> dict:
    """Convert a :class:`RunRecord` into a TOML-ready dictionary."""

    manifest_dir = record.workspace if manifest_dir is None else manifest_dir
    return {
        "run": {
            "id": record.run_id,
            "model_spec": record.model_spec,
            "scenario": record.scenario,
            "status": record.status,
            "created_at": record.created_at,
            "completed_at": record.completed_at,
            "workspace": _serialize_path(record.workspace, base_dir=manifest_dir),
            "notes": record.notes,
            "tags": list(record.tags),
            "grid_ref": record.grid_ref,
        },
        "packages": dict(record.package_versions),
        "parameter_overrides": dict(record.parameter_overrides),
        "paths": dict(record.paths),
        "metadata": dict(record.metadata),
    }


def save_run_record(record: RunRecord, path: Path | None = None):
    """Write a run manifest for a :class:`RunRecord`."""

    path = record.manifest_path if path is None else path
    write_toml(run_record_to_dict(record, manifest_dir=path.parent), path)


def load_run_record(path: Path) -> RunRecord:
    """Load a :class:`RunRecord` from a TOML manifest path."""

    data = read_toml(path)
    run_data = data.get("run", {})
    workspace = _deserialize_path(run_data.get("workspace", "."), base_dir=path.parent)
    return RunRecord(
        run_id=run_data["id"],
        model_spec=run_data["model_spec"],
        workspace=workspace,
        scenario=run_data.get("scenario"),
        grid_ref=run_data.get("grid_ref"),
        package_versions=dict(data.get("packages", {})),
        status=run_data.get("status", "registered"),
        created_at=run_data.get("created_at") or "",
        completed_at=run_data.get("completed_at"),
        notes=run_data.get("notes"),
        tags=list(run_data.get("tags", [])),
        parameter_overrides=dict(data.get("parameter_overrides", {})),
        metadata=dict(data.get("metadata", {})),
        paths=dict(data.get("paths", {})),
    )


def package_artifact_to_dict(artifact: PackageArtifact, *, manifest_dir: Path | None = None) -> dict:
    """Convert a :class:`PackageArtifact` into a TOML-ready dictionary."""

    manifest_dir = artifact.source_workspace if manifest_dir is None else manifest_dir
    return {
        "package_artifact": {
            "id": artifact.artifact_id,
            "package_type": artifact.package_type,
            "source_workspace": _serialize_path(artifact.source_workspace, base_dir=manifest_dir),
            "source_model_name": artifact.source_model_name,
            "source_package_name": artifact.source_package_name,
            "source_run_id": artifact.source_run_id,
            "source_file": artifact.source_file,
            "description": artifact.description,
            "tags": list(artifact.tags),
            "storage_mode": artifact.storage_mode,
            "derived_from_artifact_id": artifact.derived_from_artifact_id,
            "derived_from_run_id": artifact.derived_from_run_id,
        },
        "compatibility": dict(artifact.compatibility),
        "package_data": dict(artifact.package_data),
        "metadata": dict(artifact.metadata),
    }


def save_package_artifact(artifact: PackageArtifact, path: Path):
    """Write a reusable package artifact manifest to disk."""

    write_toml(package_artifact_to_dict(artifact, manifest_dir=path.parent), path)


def load_package_artifact(path: Path) -> PackageArtifact:
    """Load a :class:`PackageArtifact` from a TOML manifest path."""

    data = read_toml(path)
    artifact_data = data.get("package_artifact", {})
    source_workspace = _deserialize_path(artifact_data["source_workspace"], base_dir=path.parent)
    return PackageArtifact(
        artifact_id=artifact_data["id"],
        package_type=artifact_data["package_type"],
        source_workspace=source_workspace,
        source_model_name=artifact_data["source_model_name"],
        source_package_name=artifact_data["source_package_name"],
        source_run_id=artifact_data.get("source_run_id"),
        source_file=artifact_data.get("source_file"),
        description=artifact_data.get("description"),
        tags=list(artifact_data.get("tags", [])),
        metadata=dict(data.get("metadata", {})),
        storage_mode=artifact_data.get("storage_mode", "snapshot"),
        derived_from_artifact_id=artifact_data.get("derived_from_artifact_id"),
        derived_from_run_id=artifact_data.get("derived_from_run_id"),
        compatibility=dict(data.get("compatibility", {})),
        package_data=dict(data.get("package_data", {})),
    )
