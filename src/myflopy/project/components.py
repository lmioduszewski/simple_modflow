"""Reusable package-artifact capture, validation, derivation, and application."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from hashlib import sha1
import json
from pathlib import Path
from typing import Any

import flopy
import numpy as np


SUPPORTED_PACKAGE_ARTIFACT_TYPES = {"chd", "drn", "ghb", "ic", "npf", "rch", "uzf", "lak", "sfr", "mvr"}
PACKAGE_ARTIFACT_APPLY_ORDER = {
    "ic": 10,
    "npf": 20,
    "rch": 30,
    "chd": 40,
    "drn": 50,
    "ghb": 60,
    "uzf": 70,
    "lak": 80,
    "sfr": 90,
    "mvr": 100,
}


class PackageCompatibilityError(ValueError):
    """Raised when a reused :class:`PackageArtifact` does not fit the target model.

    A captured package definition carries assumptions about the model it came from
    (grid type, cell count, layer count, time discretization). This error is raised
    when one of those assumptions is violated as the artifact is applied to a
    different model -- for example reusing a DISV package built for one Voronoi grid
    on a model with a different ``ncpl``. A :class:`ValueError`, so existing
    value-error handling still catches it.
    """


@dataclass(slots=True)
class PackageArtifact:
    """Reusable package definition captured from a source run or in-memory model.

    The artifact stores enough information to:
    - trace where the package originally came from
    - validate whether it can be reused on a new model
    - recreate supported package types on a compatible target model
    """

    artifact_id: str
    package_type: str
    source_workspace: Path
    source_model_name: str
    source_package_name: str
    compatibility: dict[str, Any]
    package_data: dict[str, Any]
    source_run_id: str | None = None
    source_file: str | None = None
    description: str | None = None
    tags: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    storage_mode: str = "snapshot"
    derived_from_artifact_id: str | None = None
    derived_from_run_id: str | None = None

    def __post_init__(self):
        """Coerce ``source_workspace`` to a ``Path`` and lowercase the package type."""

        self.source_workspace = Path(self.source_workspace)
        self.package_type = self.package_type.lower()


def _grid_ncpl_array(modelgrid) -> list[int]:
    """Return per-layer cell counts from a model grid as a plain Python list."""

    ncpl = getattr(modelgrid, "ncpl", None)
    if isinstance(ncpl, (int, np.integer)):
        return [int(ncpl)] * int(modelgrid.nlay)
    return np.asarray(ncpl, dtype=int).reshape(-1).astype(int).tolist()


def build_grid_compatibility_snapshot(model) -> dict[str, Any]:
    """Build a compact exact-grid compatibility fingerprint for a model."""

    modelgrid = model.modelgrid
    x = np.asarray(modelgrid.xcellcenters, dtype=float).reshape(-1)
    y = np.asarray(modelgrid.ycellcenters, dtype=float).reshape(-1)
    payload: dict[str, Any] = {
        "grid_type": model.grid_type,
        "nlay": int(modelgrid.nlay),
        "ncpl_array": _grid_ncpl_array(modelgrid),
        "node_count": int(sum(_grid_ncpl_array(modelgrid))),
        "xcellcenters": np.round(x, 6).tolist(),
        "ycellcenters": np.round(y, 6).tolist(),
    }

    iverts = getattr(modelgrid, "iverts", None)
    if iverts is not None:
        payload["iverts"] = [[int(vertex) for vertex in cell] for cell in iverts]

    payload["grid_hash"] = sha1(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return payload


def package_artifact_apply_order(package_type: str) -> int:
    """Return the preferred application order for a reusable package type."""

    return PACKAGE_ARTIFACT_APPLY_ORDER.get(package_type.lower(), 999)


def _build_grid_compatibility_issues(artifact: PackageArtifact, model) -> list[str]:
    """Return human-readable grid compatibility issues instead of raising."""

    target = build_grid_compatibility_snapshot(model)
    source = artifact.compatibility
    issues: list[str] = []

    if source.get("grid_type") != target.get("grid_type"):
        issues.append(
            f"expects grid_type={source.get('grid_type')}, got {target.get('grid_type')}"
        )
    if source.get("nlay") != target.get("nlay"):
        issues.append(
            f"expects nlay={source.get('nlay')}, got {target.get('nlay')}"
        )
    if source.get("node_count") != target.get("node_count"):
        issues.append(
            f"expects node_count={source.get('node_count')}, got {target.get('node_count')}"
        )
    if source.get("grid_hash") != target.get("grid_hash"):
        issues.append("was captured from a different grid layout")
    return issues


def validate_package_artifact_compatibility(artifact: PackageArtifact, model):
    """Raise if an artifact cannot be safely reused on the target model.

    This is the strict exact-grid compatibility gate used for current snapshot
    artifacts.
    """

    issues = _build_grid_compatibility_issues(artifact, model)
    if issues:
        raise PackageCompatibilityError(
            f"Package artifact '{artifact.artifact_id}' is incompatible with model '{model.name}': "
            + "; ".join(issues)
            + "."
        )


def _normalize_record_item(value):
    """Normalize NumPy/FloPy values into TOML/JSON-friendly Python objects."""

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_normalize_record_item(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): _normalize_record_item(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_normalize_record_item(item) for item in value]
    if isinstance(value, list):
        return [_normalize_record_item(item) for item in value]
    return value


def _extract_mf6_value(package, attribute: str, default=None):
    """Read a package attribute, unwrapping FloPy data containers where needed."""

    value = getattr(package, attribute, None)
    if value is None:
        return default

    getter = getattr(value, "get_data", None)
    if callable(getter):
        try:
            value = getter()
        except TypeError:
            pass
    elif hasattr(value, "array"):
        value = value.array

    return _normalize_record_item(value)


def _extract_stress_period_data(package) -> dict[str, list[list[Any]]]:
    """Serialize FloPy stress-period data into a plain nested Python structure."""

    raw = package.stress_period_data.get_data()
    stress_period_data: dict[str, list[list[Any]]] = {}
    for period, rows in raw.items():
        normalized_rows = []
        for row in rows.tolist():
            if not isinstance(row, (list, tuple)):
                row = [row]
            normalized_rows.append([_normalize_record_item(item) for item in row])
        stress_period_data[str(int(period))] = normalized_rows
    return stress_period_data


def _restore_rows(rows: list[list[Any]], *, tuple_indices: tuple[int, ...] = ()) -> list[list[Any]]:
    """Restore tuple-valued fields when rebuilding package rows from manifests."""

    restored = []
    for row in rows:
        normalized = list(row)
        for index in tuple_indices:
            if index < len(normalized) and isinstance(normalized[index], list):
                normalized[index] = tuple(normalized[index])
        restored.append(normalized)
    return restored


def _resolve_package(model, package_name: str):
    """Resolve a package from a model or simulation by common name variants."""

    package = model.gwf.get_package(package_name)
    if package is None:
        package = getattr(model.gwf, package_name.lower(), None)
    if package is None and getattr(model, "sim", None) is not None:
        package = model.sim.get_package(package_name)
    if package is None and getattr(model, "sim", None) is not None:
        package = getattr(model.sim, package_name.lower(), None)
    if package is None:
        raise FileNotFoundError(f"Package '{package_name}' not found on model '{model.name}'.")
    return package


def _deep_merge_dicts(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge nested dictionaries without mutating the inputs."""

    merged = deepcopy(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge_dicts(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def build_package_artifact(
    model,
    *,
    artifact_id: str,
    package_name: str,
    description: str | None = None,
    tags: list[str] | None = None,
    metadata: dict[str, Any] | None = None,
) -> PackageArtifact:
    """Create a reusable package artifact from an existing model package.

    Parameters
    ----------
    model
        Source model containing the package to capture.
    artifact_id
        Stable id used to store and reference the captured artifact later.
    package_name
        Model package key/name to capture, such as ``"chd"`` or ``"sfr"``.
    description, tags, metadata
        Optional human-facing provenance fields stored with the artifact.
    """

    package = _resolve_package(model, package_name)
    package_type = getattr(package, "package_type", "").lower()
    if package_type not in SUPPORTED_PACKAGE_ARTIFACT_TYPES:
        raise NotImplementedError(
            f"Package artifact support is not implemented for package type '{package_type}'."
        )

    if package_type in {"chd", "drn", "ghb", "rch"}:
        package_data = {
            "stress_period_data": _extract_stress_period_data(package),
        }
        auxiliary = _extract_mf6_value(package, "auxiliary")
        if auxiliary is not None:
            package_data["auxiliary"] = list(auxiliary)
    elif package_type == "ic":
        package_data = {
            "strt": _extract_mf6_value(package, "strt"),
        }
    elif package_type == "npf":
        package_data = {
            "icelltype": _extract_mf6_value(package, "icelltype"),
            "k": _extract_mf6_value(package, "k"),
            "k33": _extract_mf6_value(package, "k33"),
            "perched": _extract_mf6_value(package, "perched", False),
            "save_specific_discharge": _extract_mf6_value(package, "save_specific_discharge", True),
        }
    elif package_type == "uzf":
        package_data = {
            "mover": _extract_mf6_value(package, "mover", False),
            "simulate_et": _extract_mf6_value(package, "simulate_et", False),
            "linear_gwet": _extract_mf6_value(package, "linear_gwet", False),
            "square_gwet": _extract_mf6_value(package, "square_gwet", False),
            "simulate_gwseep": _extract_mf6_value(package, "simulate_gwseep", False),
            "unsat_etwc": _extract_mf6_value(package, "unsat_etwc", False),
            "unsat_etae": _extract_mf6_value(package, "unsat_etae", False),
            "nuzfcells": _extract_mf6_value(package, "nuzfcells"),
            "ntrailwaves": _extract_mf6_value(package, "ntrailwaves", 7),
            "nwavesets": _extract_mf6_value(package, "nwavesets", 40),
            "packagedata": _extract_mf6_value(package, "packagedata"),
            "perioddata": _normalize_record_item(package.perioddata.get_data()),
        }
    elif package_type == "lak":
        package_data = {
            "mover": _extract_mf6_value(package, "mover", False),
            "surfdep": _extract_mf6_value(package, "surfdep"),
            "maximum_iterations": _extract_mf6_value(package, "maximum_iterations"),
            "maximum_stage_change": _extract_mf6_value(package, "maximum_stage_change"),
            "time_conversion": _extract_mf6_value(package, "time_conversion"),
            "length_conversion": _extract_mf6_value(package, "length_conversion"),
            "nlakes": _extract_mf6_value(package, "nlakes"),
            "noutlets": _extract_mf6_value(package, "noutlets"),
            "ntables": _extract_mf6_value(package, "ntables"),
            "packagedata": _extract_mf6_value(package, "packagedata"),
            "connectiondata": _extract_mf6_value(package, "connectiondata"),
            "tables": _extract_mf6_value(package, "tables"),
            "outlets": _extract_mf6_value(package, "outlets"),
            "perioddata": _normalize_record_item(package.perioddata.get_data()),
        }
    elif package_type == "sfr":
        package_data = {
            "mover": _extract_mf6_value(package, "mover", False),
            "nreaches": _extract_mf6_value(package, "nreaches"),
            "maximum_picard_iterations": _extract_mf6_value(package, "maximum_picard_iterations"),
            "maximum_iterations": _extract_mf6_value(package, "maximum_iterations"),
            "maximum_depth_change": _extract_mf6_value(package, "maximum_depth_change"),
            "length_conversion": _extract_mf6_value(package, "length_conversion"),
            "time_conversion": _extract_mf6_value(package, "time_conversion"),
            "packagedata": _extract_mf6_value(package, "packagedata"),
            "connectiondata": _extract_mf6_value(package, "connectiondata"),
            "diversions": _extract_mf6_value(package, "diversions"),
            "perioddata": _normalize_record_item(package.perioddata.get_data()),
        }
    elif package_type == "mvr":
        package_data = {
            "modelnames": _extract_mf6_value(package, "modelnames", False),
            "maxmvr": _extract_mf6_value(package, "maxmvr"),
            "maxpackages": _extract_mf6_value(package, "maxpackages"),
            "packages": _extract_mf6_value(package, "packages"),
            "perioddata": _normalize_record_item(package.perioddata.get_data()),
        }
    else:
        raise NotImplementedError(
            f"Package artifact support is not implemented for package type '{package_type}'."
        )

    compatibility = build_grid_compatibility_snapshot(model)
    source_file = getattr(package, "filename", None)
    return PackageArtifact(
        artifact_id=artifact_id,
        package_type=package_type,
        source_workspace=model.workspace,
        source_model_name=model.name,
        source_package_name=package_name.upper(),
        source_run_id=None,
        source_file=source_file,
        description=description,
        tags=[] if tags is None else list(tags),
        metadata={} if metadata is None else dict(metadata),
        compatibility=compatibility,
        package_data=package_data,
    )


def derive_package_artifact(
    artifact: PackageArtifact,
    *,
    artifact_id: str,
    description: str | None = None,
    tags: list[str] | None = None,
    metadata_updates: dict[str, Any] | None = None,
    package_data_updates: dict[str, Any] | None = None,
) -> PackageArtifact:
    """Create a new artifact derived from an existing snapshot artifact.

    This is a cautious first step toward more flexible reuse: the derived artifact keeps
    the same exact-grid compatibility fingerprint, records its lineage, and allows small
    data/metadata edits without pretending to be a cross-grid recipe.
    """

    metadata = dict(artifact.metadata)
    if metadata_updates:
        metadata = _deep_merge_dicts(metadata, metadata_updates)

    package_data = dict(artifact.package_data)
    if package_data_updates:
        package_data = _deep_merge_dicts(package_data, package_data_updates)

    return PackageArtifact(
        artifact_id=artifact_id,
        package_type=artifact.package_type,
        source_workspace=artifact.source_workspace,
        source_model_name=artifact.source_model_name,
        source_package_name=artifact.source_package_name,
        compatibility=deepcopy(artifact.compatibility),
        package_data=package_data,
        source_run_id=artifact.source_run_id,
        source_file=artifact.source_file,
        description=artifact.description if description is None else description,
        tags=list(artifact.tags) if tags is None else list(tags),
        metadata=metadata,
        storage_mode=artifact.storage_mode,
        derived_from_artifact_id=artifact.artifact_id,
        derived_from_run_id=artifact.source_run_id or artifact.derived_from_run_id,
    )


def _artifact_dependencies(artifact: PackageArtifact) -> set[str]:
    """Infer package-type dependencies implied by a stored artifact."""

    dependencies: set[str] = set()
    if artifact.package_type == "mvr":
        for row in artifact.package_data.get("packages", []) or []:
            for item in row if isinstance(row, list) else [row]:
                if isinstance(item, str) and item.lower() in SUPPORTED_PACKAGE_ARTIFACT_TYPES:
                    dependencies.add(item.lower())
    if artifact.package_type in {"uzf", "lak", "sfr"} and artifact.package_data.get("mover", False):
        dependencies.add("mvr")
    dependencies.discard(artifact.package_type)
    return dependencies


def _artifact_structural_issues(artifact: PackageArtifact) -> list[str]:
    """Return structural data issues found within one stored artifact snapshot."""

    data = artifact.package_data
    issues: list[str] = []

    if artifact.package_type == "uzf":
        packagedata = data.get("packagedata") or []
        perioddata = data.get("perioddata") or {}
        if int(data.get("nuzfcells") or 0) != len(packagedata):
            issues.append(
                f"stores nuzfcells={data.get('nuzfcells')} but packagedata has {len(packagedata)} records"
            )
        if not perioddata:
            issues.append("has no perioddata records")
    elif artifact.package_type == "lak":
        packagedata = data.get("packagedata") or []
        connectiondata = data.get("connectiondata") or []
        if int(data.get("nlakes") or 0) != len(packagedata):
            issues.append(
                f"stores nlakes={data.get('nlakes')} but packagedata has {len(packagedata)} records"
            )
        if not connectiondata:
            issues.append("has no connectiondata records")
    elif artifact.package_type == "sfr":
        packagedata = data.get("packagedata") or []
        connectiondata = data.get("connectiondata") or []
        if int(data.get("nreaches") or 0) != len(packagedata):
            issues.append(
                f"stores nreaches={data.get('nreaches')} but packagedata has {len(packagedata)} records"
            )
        if not connectiondata and len(packagedata) > 1:
            issues.append("has no connectiondata records for a multi-reach network")
    elif artifact.package_type == "mvr":
        if not (data.get("packages") or []):
            issues.append("defines no source/target packages")
        if not (data.get("perioddata") or {}):
            issues.append("has no perioddata records")

    return issues


def validate_package_artifact_reference(
    artifact: PackageArtifact,
    model,
    *,
    available_package_types: set[str] | None = None,
    referenced_package_types: set[str] | None = None,
) -> dict[str, Any]:
    """Validate one artifact against a target model and package context.

    Unlike the strict compatibility check, this helper also reports dependency
    issues such as an ``mvr`` artifact being referenced without ``lak`` or
    ``uzf`` packages available.
    """

    available_package_types = set() if available_package_types is None else set(available_package_types)
    referenced_package_types = set() if referenced_package_types is None else set(referenced_package_types)

    issues: list[str] = []
    issues.extend(_build_grid_compatibility_issues(artifact, model))
    issues.extend(_artifact_structural_issues(artifact))

    dependencies = sorted(_artifact_dependencies(artifact))
    missing_dependencies = sorted(
        dependency
        for dependency in dependencies
        if dependency not in available_package_types and dependency not in referenced_package_types
    )
    if missing_dependencies:
        issues.append(
            "missing dependent packages: "
            + ", ".join(missing_dependencies)
            + f" (required by artifact type '{artifact.package_type}')"
        )

    return {
        "artifact_id": artifact.artifact_id,
        "package_type": artifact.package_type,
        "valid": not issues,
        "dependencies": dependencies,
        "missing_dependencies": missing_dependencies,
        "issue_count": len(issues),
        "message": "ok" if not issues else "; ".join(issues),
    }


def format_validation_report(report_rows: list[dict[str, Any]]) -> str:
    """Format a validation report into a readable multi-line error message."""

    lines: list[str] = []
    for row in report_rows:
        if row.get("valid", False):
            continue
        lines.append(
            f"- {row['package_key']} -> {row['artifact_id']} ({row['package_type']}): {row['message']}"
        )
    return "\n".join(lines)


def _restore_stress_period_data(package_data: dict[str, Any]) -> dict[int, list[list[Any]]]:
    """Restore serialized stress-period data to the shape FloPy expects."""

    restored: dict[int, list[list[Any]]] = {}
    for period, rows in package_data["stress_period_data"].items():
        restored_rows = []
        for row in rows:
            row = list(row)
            if row and isinstance(row[0], list):
                row[0] = tuple(row[0])
            restored_rows.append(row)
        restored[int(period)] = restored_rows
    return restored


def apply_package_artifact(model, artifact: PackageArtifact, *, validate: bool = True):
    """Attach a supported stored package artifact to a target model.

    Parameters
    ----------
    model
        Target model that will receive the recreated package.
    artifact
        Stored package artifact to recreate.
    validate
        When ``True``, enforce exact-grid compatibility before attaching.
    """

    if validate:
        validate_package_artifact_compatibility(artifact, model)

    if artifact.package_type in {"chd", "drn", "ghb", "rch"}:
        stress_period_data = _restore_stress_period_data(artifact.package_data)
        auxiliary = artifact.package_data.get("auxiliary")

    if artifact.package_type == "chd":
        package = flopy.mf6.ModflowGwfchd(
            model.gwf,
            print_input=False,
            print_flows=False,
            save_flows=True,
            filename=f"{model.name}.chd",
            pname="chd",
            stress_period_data=stress_period_data,
        )
    elif artifact.package_type == "drn":
        package = flopy.mf6.ModflowGwfdrn(
            model.gwf,
            print_input=False,
            print_flows=False,
            save_flows=True,
            filename=f"{model.name}.drn",
            pname="drn",
            stress_period_data=stress_period_data,
        )
    elif artifact.package_type == "ghb":
        package = flopy.mf6.ModflowGwfghb(
            model.gwf,
            print_input=False,
            print_flows=False,
            save_flows=True,
            filename=f"{model.name}.ghb",
            pname="ghb",
            stress_period_data=stress_period_data,
            auxiliary=auxiliary,
        )
    elif artifact.package_type == "ic":
        package = flopy.mf6.modflow.mfgwfic.ModflowGwfic(
            model.gwf,
            pname="ic",
            strt=artifact.package_data["strt"],
            filename=f"{model.name}.ic",
        )
    elif artifact.package_type == "npf":
        package = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            model.gwf,
            pname="npf",
            icelltype=artifact.package_data["icelltype"],
            k=artifact.package_data["k"],
            perched=artifact.package_data.get("perched", False),
            k33=artifact.package_data.get("k33"),
            save_flows=True,
            save_saturation=True,
            save_specific_discharge=artifact.package_data.get("save_specific_discharge", True),
            filename=f"{model.name}.npf",
        )
    elif artifact.package_type == "rch":
        maxbound = max((len(rows) for rows in stress_period_data.values()), default=0)
        package = flopy.mf6.ModflowGwfrch(
            model.gwf,
            pname="rch",
            print_input=False,
            print_flows=False,
            save_flows=True,
            maxbound=maxbound,
            stress_period_data=stress_period_data,
            filename=f"{model.name}.rch",
            auxiliary=auxiliary,
        )
    elif artifact.package_type == "uzf":
        package = flopy.mf6.ModflowGwfuzf(
            model.gwf,
            save_flows=True,
            print_flows=False,
            pname="uzf",
            nuzfcells=artifact.package_data["nuzfcells"],
            packagedata=_restore_rows(artifact.package_data["packagedata"], tuple_indices=(1,)),
            perioddata={int(per): _restore_rows(rows) for per, rows in artifact.package_data["perioddata"].items()},
            mover=artifact.package_data.get("mover", False),
            simulate_et=artifact.package_data.get("simulate_et", False),
            linear_gwet=artifact.package_data.get("linear_gwet", False),
            square_gwet=artifact.package_data.get("square_gwet", False),
            simulate_gwseep=artifact.package_data.get("simulate_gwseep", False),
            unsat_etwc=artifact.package_data.get("unsat_etwc", False),
            unsat_etae=artifact.package_data.get("unsat_etae", False),
            budget_filerecord=f"{model.name}_budget.uzf",
            budgetcsv_filerecord=f"{model.name}_uzf_budget.csv",
            package_convergence_filerecord=f"{model.name}_uzf_package_convergence.csv",
            ntrailwaves=artifact.package_data.get("ntrailwaves", 7),
            nwavesets=artifact.package_data.get("nwavesets", 40),
            filename=f"{model.name}.uzf",
        )
    elif artifact.package_type == "lak":
        package = flopy.mf6.ModflowGwflak(
            model.gwf,
            save_flows=True,
            print_input=False,
            print_flows=False,
            print_stage=True,
            stage_filerecord=f"{model.name}_stage.lak",
            budget_filerecord=f"{model.name}_budget.lak",
            budgetcsv_filerecord=f"{model.name}_lake_budget.csv",
            package_convergence_filerecord=f"{model.name}_lake_convergence.csv",
            mover=artifact.package_data.get("mover", False),
            surfdep=artifact.package_data.get("surfdep"),
            time_conversion=artifact.package_data.get("time_conversion"),
            length_conversion=artifact.package_data.get("length_conversion"),
            nlakes=artifact.package_data["nlakes"],
            noutlets=artifact.package_data["noutlets"],
            ntables=artifact.package_data["ntables"],
            packagedata=_restore_rows(artifact.package_data["packagedata"]),
            connectiondata=_restore_rows(artifact.package_data["connectiondata"], tuple_indices=(2,)),
            tables=artifact.package_data.get("tables"),
            outlets=artifact.package_data.get("outlets"),
            perioddata={int(per): _restore_rows(rows) for per, rows in artifact.package_data["perioddata"].items()},
            filename=f"{model.name}.lak",
            pname="lak",
            maximum_iterations=artifact.package_data.get("maximum_iterations"),
            maximum_stage_change=artifact.package_data.get("maximum_stage_change"),
        )
    elif artifact.package_type == "sfr":
        package = flopy.mf6.ModflowGwfsfr(
            model.gwf,
            save_flows=True,
            print_input=True,
            print_flows=True,
            pname="sfr",
            nreaches=artifact.package_data["nreaches"],
            packagedata=_restore_rows(artifact.package_data["packagedata"], tuple_indices=(1,)),
            connectiondata=_restore_rows(artifact.package_data["connectiondata"]),
            perioddata={int(per): _restore_rows(rows) for per, rows in artifact.package_data["perioddata"].items()},
            maximum_picard_iterations=artifact.package_data.get("maximum_picard_iterations"),
            maximum_iterations=artifact.package_data.get("maximum_iterations"),
            maximum_depth_change=artifact.package_data.get("maximum_depth_change"),
            budget_filerecord="sfr_budget.sfr",
            stage_filerecord="sfr_stage.sfr",
            length_conversion=artifact.package_data.get("length_conversion"),
            time_conversion=artifact.package_data.get("time_conversion"),
            mover=artifact.package_data.get("mover", False),
            diversions=None if artifact.package_data.get("diversions") is None else _restore_rows(artifact.package_data["diversions"]),
        )
    elif artifact.package_type == "mvr":
        package = flopy.mf6.ModflowGwfmvr(
            model.gwf,
            modelnames=artifact.package_data.get("modelnames", False),
            budget_filerecord=f"{model.name}.mvr.bud",
            budgetcsv_filerecord=f"{model.name}.mvr.csv",
            maxmvr=artifact.package_data["maxmvr"],
            maxpackages=artifact.package_data["maxpackages"],
            packages=_restore_rows(artifact.package_data["packages"]),
            perioddata={int(per): _restore_rows(rows) for per, rows in artifact.package_data["perioddata"].items()},
            filename=f"{model.name}.mvr",
            pname="mvr",
        )
    else:
        raise NotImplementedError(
            f"Package artifact support is not implemented for package type '{artifact.package_type}'."
        )

    attached = getattr(model, "attached_package_artifacts", {})
    attached[artifact.package_type] = artifact.artifact_id
    model.attached_package_artifacts = attached
    return package
