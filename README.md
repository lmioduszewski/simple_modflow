# simple_modflow

## Developer Guide

- Codebase structure map: [docs/codebase_structure.md](docs/codebase_structure.md)
- Preferred API guide: [docs/preferred_api.md](docs/preferred_api.md)

Utilities for building, running, and post-processing MODFLOW 6 models with FloPy.

Preferred vector-builder names for GIS-driven package/material inputs are:
`DRNFromVector`, `GHBFromVector`, `CHDFromVector`, `RCHFromVector`, and
`KFromVector`. Older names remain supported for compatibility.

## Status

This project is in the middle of a packaging and architecture cleanup. The codebase now has clearer internal
subpackages for MF6 grid helpers and simulation helpers, but some legacy façade modules are still kept for
compatibility while the refactor continues.

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
