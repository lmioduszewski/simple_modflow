# myflopy

## Developer Guide

- Codebase structure map: [docs/codebase_structure.md](docs/codebase_structure.md)
- Core specification API: [docs/core_specs.md](docs/core_specs.md)
- Preferred API guide: [docs/preferred_api.md](docs/preferred_api.md)
- Projects, reusable packages, and swappable grids: [docs/project_workflow.md](docs/project_workflow.md)
- Refactor review and consolidation strategy: [docs/refactor_review_report.md](docs/refactor_review_report.md)
- Local run clutter guide: [docs/local_run_clutter.md](docs/local_run_clutter.md)

Utilities for building, running, and post-processing MODFLOW 6 models with FloPy.

Use `import myflopy`.

Preferred vector-builder names for GIS-driven package/material inputs are:
`DRNFromVector`, `GHBFromVector`, `CHDFromVector`, `RCHFromVector`, and
`KFromVector`.

## Status

This project is undergoing a clean architecture refactor. Existing APIs,
aliases, and internal layouts may be removed when they do not fit the new
design; backward compatibility with `simple_modflow` is not a goal.

## Installation

Core install:

```bash
pip install .
```

Optional visualization extras:

```bash
pip install .[viz]
```

Developer tooling:

```bash
pip install .[dev]
```

## Notes

- Example notebooks, sample outputs, and ad hoc artifacts have been moved under `examples/` so they no longer live
  inside the importable package tree.
- Some integrations in this repo are environment-specific and are not declared as install dependencies here,
  especially custom `figs` wrappers plus GRASS/QGIS-related tooling.

## Test Helper

If Windows keeps locking old pytest temp folders, use the repo helper below. It keeps only one temp test
workspace at a time under `.pytest-work` and recreates it fresh for each run.

```powershell
.\scripts\pytest_local.cmd tests\test_mp3du_particles.py
```
