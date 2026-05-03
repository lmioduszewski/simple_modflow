"""Small convenience helpers for discovery/import workflows."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from simple_modflow.project.catalog import ProjectCatalog
from simple_modflow.project.discovery import discover_existing_runs
from simple_modflow.project.specs import ModelSpec


def summarize_discovered_runs(search_root: Path) -> pd.DataFrame:
    """Discover legacy/native MF6 runs below ``search_root`` and summarize them."""

    rows = []
    for run in discover_existing_runs(search_root):
        rows.append(
            {
                "run_id": run.run_id,
                "model_name": run.model_name,
                "workspace": str(run.workspace),
                "family": run.metadata.get("family"),
                "grid_type": run.metadata.get("grid_type"),
                "has_mfsim": run.metadata.get("has_mfsim", False),
                "has_heads": run.metadata.get("has_heads", False),
                "has_budget": run.metadata.get("has_budget", False),
                "tags": list(run.tags),
            }
        )
    return pd.DataFrame(rows)


def import_run_archive(
    search_root: Path,
    catalog_root: Path,
    *,
    project_name: str | None = None,
    model_spec: str = "imported_model",
    register_model_spec: bool = True,
    overwrite: bool = False,
    skip_existing: bool = True,
) -> ProjectCatalog:
    """Import a directory of existing MF6 runs into a lightweight catalog.

    Parameters
    ----------
    search_root
        Directory tree containing existing MF6 workspaces.
    catalog_root
        Destination root for the lightweight project catalog.
    project_name
        Optional display name for the new catalog.
    model_spec
        Model spec id assigned to imported runs.
    register_model_spec
        Whether to create/update the corresponding model spec automatically.
    overwrite
        Whether existing imported manifests may be replaced.
    skip_existing
        When ``True``, existing manifests are reused rather than raising.
    """

    catalog = ProjectCatalog(catalog_root, name=project_name)
    if register_model_spec:
        try:
            catalog.register_model_spec(
                ModelSpec(
                    name=model_spec,
                    description=f"Imported legacy runs discovered under {Path(search_root)}",
                    grid_ref="imported",
                    tags=["imported", "legacy"],
                ),
                overwrite=overwrite,
            )
        except FileExistsError:
            if overwrite:
                raise
    catalog.import_existing_runs(
        search_root,
        model_spec=model_spec,
        overwrite=overwrite,
        skip_existing=skip_existing,
    )
    return catalog
