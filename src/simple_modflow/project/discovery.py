from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

from simple_modflow.project.specs import RunRecord, default_run_artifact_paths


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def _infer_family(run_name: str) -> str | None:
    match = re.match(r"^[A-Za-z]+", run_name)
    if match is None:
        return None
    return match.group(0).lower()


@dataclass(slots=True)
class DiscoveredRun:
    run_id: str
    workspace: Path
    model_name: str
    model_file: Path | None = None
    tags: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    paths: dict[str, str] = field(default_factory=dict)

    def to_run_record(self, *, model_spec: str = "imported_model", status: str = "imported") -> RunRecord:
        return RunRecord(
            run_id=self.run_id,
            model_spec=model_spec,
            workspace=self.workspace,
            status=status,
            tags=list(self.tags),
            metadata=dict(self.metadata),
            paths=dict(self.paths),
        )


def _infer_model_name_from_workspace(workspace: Path) -> str | None:
    model_files = sorted(workspace.glob("*.model"))
    if model_files:
        return model_files[0].stem

    name_files = sorted(path for path in workspace.glob("*.nam") if path.name.lower() != "mfsim.nam")
    if name_files:
        return name_files[0].stem

    for suffix in (".hds", ".cbc", ".disv", ".disu"):
        files = sorted(workspace.glob(f"*{suffix}"))
        if files:
            return files[0].stem
    return None


def discover_existing_runs(search_root: Path) -> list[DiscoveredRun]:
    search_root = Path(search_root)
    workspaces = sorted({path.parent for path in search_root.rglob("mfsim.nam")})
    model_file_workspaces = sorted(path.parent for path in search_root.rglob("*.model"))
    for workspace in model_file_workspaces:
        if workspace not in workspaces:
            workspaces.append(workspace)
    workspaces = sorted(set(workspaces))

    counts: dict[str, int] = {}
    discovered: list[DiscoveredRun] = []

    for workspace in workspaces:
        model_name = _infer_model_name_from_workspace(workspace)
        if model_name is None:
            continue
        model_file = next(iter(sorted(workspace.glob(f"{model_name}.model"))), None)
        relative_workspace = workspace.relative_to(search_root)

        run_id = model_name
        counts.setdefault(run_id, 0)
        counts[run_id] += 1
        if counts[run_id] > 1:
            suffix = _slug(relative_workspace.as_posix())
            run_id = f"{model_name}__{suffix}"

        paths = default_run_artifact_paths(model_name)
        paths["workspace"] = workspace.as_posix()
        if model_file is not None:
            paths["model_object_file"] = model_file.name

        for key, pattern in [
            ("simulation_name_file", "mfsim.nam"),
            ("listing_file", "mfsim.lst"),
            ("heads_file", f"{model_name}.hds"),
            ("budget_file", f"{model_name}.cbc"),
            ("model_name_file", f"{model_name}.nam"),
            ("disv_file", f"{model_name}.disv"),
            ("disu_file", f"{model_name}.disu"),
        ]:
            path = workspace / pattern
            if path.exists():
                paths[key] = path.name

        family = _infer_family(model_name)
        tags = []
        if family is not None:
            tags.append(family)
        if len(relative_workspace.parts) > 1:
            tags.extend(_slug(part).lower() for part in relative_workspace.parts[:-1])

        metadata = {
            "source_root": str(search_root),
            "relative_workspace": relative_workspace.as_posix(),
            "family": family,
            "has_mfsim": (workspace / "mfsim.nam").exists(),
            "has_heads": "heads_file" in paths,
            "has_budget": "budget_file" in paths,
            "grid_type": "disu" if "disu_file" in paths else ("disv" if "disv_file" in paths else None),
        }

        discovered.append(
            DiscoveredRun(
                run_id=run_id,
                workspace=workspace,
                model_name=model_name,
                model_file=model_file,
                tags=sorted(set(tags)),
                metadata=metadata,
                paths=paths,
            )
        )

    return discovered
