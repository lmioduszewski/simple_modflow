# `simple_modflow` Codebase Structure

This document is the maintained map of the project layout.

Use it when you want to answer questions like:
- "Where does run loading happen?"
- "Where is the high-level model API?"
- "Where do grouped comparisons live?"
- "Which file should I edit for UZF artifacts or lazy budget loading?"

The guide focuses on the active, maintained code paths. Legacy or archived code
is called out explicitly so it is easy to avoid unless you are intentionally
working on it.

## Top Level

- [README.md](C:/Users/lukem/Python/Projects/simple_modflow/README.md)
  Main project readme.
- [pyproject.toml](C:/Users/lukem/Python/Projects/simple_modflow/pyproject.toml)
  Packaging and build metadata.
- [docs/codebase_structure.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/codebase_structure.md)
  This structure guide.
- [docs/preferred_api.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/preferred_api.md)
  Preferred usage guide for the modern public API.
- [docs/mf6io_reference.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/mf6io_reference.md)
  Internal summary of the external `mf6io.pdf` reference, especially the MF6
  package-output budget semantics that `simple_modflow` relies on.
- [src](C:/Users/lukem/Python/Projects/simple_modflow/src)
  All importable library code.
- [tests](C:/Users/lukem/Python/Projects/simple_modflow/tests)
  Automated regression and workflow tests.
- [examples](C:/Users/lukem/Python/Projects/simple_modflow/examples)
  Example notebooks and walkthrough material.

## `src/simple_modflow`

- [src/simple_modflow/__init__.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/__init__.py)
  Top-level lazy public API. This is where names like `SimulationBase`,
  `ProjectCatalog`, `ModelGroup`, `TriangleGrid`, and `explore_runs` are
  exported.

### `modflow`

Core MODFLOW-oriented functionality. Most day-to-day work is in the `mf6`
subdirectory.

#### `modflow/calcs`

- [src/simple_modflow/modflow/calcs/calibration.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/calcs/calibration.py)
  Calibration and comparison helpers built on top of model outputs.

#### `modflow/mf6`

This is the main MF6 implementation area.

- [src/simple_modflow/modflow/mf6/simplemodel.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simplemodel.py)
  Higher-level model assembly helpers from the older workflow surface.
- [src/simple_modflow/modflow/mf6/boundaries.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/boundaries.py)
  Boundary geometry processing and cell-intersection helpers.
- [src/simple_modflow/modflow/mf6/boundary_support.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/boundary_support.py)
  Shared support functions for boundary packages.
- [src/simple_modflow/modflow/mf6/drn.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/drn.py)
  Drain-specific workflow helpers. Preferred vector-builder name:
  `DRNFromVector` (compatibility name: `DRN`).
- [src/simple_modflow/modflow/mf6/ghb.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/ghb.py)
  General-head-boundary workflow helpers. Preferred vector-builder name:
  `GHBFromVector` (compatibility name: `GHB`).
- [src/simple_modflow/modflow/mf6/chd.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/chd.py)
  Vector-driven constant-head builders. Preferred builder name:
  `CHDFromVector`.
- [src/simple_modflow/modflow/mf6/recharge.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/recharge.py)
  Recharge-specific workflow helpers. Preferred vector-builder name:
  `RCHFromVector` (compatibility name: `RechargeFromShp`).
- [src/simple_modflow/modflow/mf6/kflow.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/kflow.py)
  Vector-driven hydraulic-conductivity/material-property builders. Preferred
  builder name: `KFromVector`.
- [src/simple_modflow/modflow/mf6/lakes.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/lakes.py)
  Lake geometry, connection, and package-data helpers.
- [src/simple_modflow/modflow/mf6/sfr.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/sfr.py)
  SFR workflow helpers and input builders.
- [src/simple_modflow/modflow/mf6/surface_water_validation.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/surface_water_validation.py)
  Shared validation/reporting layer for coupled `SFR`, `LAK`, and `MVR`
  workflows.
- [src/simple_modflow/modflow/mf6/uzf.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/uzf.py)
  UZF workflow helpers and package-data preparation.
- [src/simple_modflow/modflow/mf6/budget.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/budget.py)
  Budget-access layer. This is the place to look for single-model budget
  dataframe logic and package budget normalization.
- [src/simple_modflow/modflow/mf6/budget_tables.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/budget_tables.py)
  Budget dataframe shaping, package-output reshaping, and observation-area
  budget table helpers.
- [src/simple_modflow/modflow/mf6/budget_plotting.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/budget_plotting.py)
  Budget plotting and observation-area visualization helpers.
- [src/simple_modflow/modflow/mf6/headsplus.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/headsplus.py)
  Core heads reader/accessor. This is now the reader-focused shell that
  delegates observation and plotting work to the helper modules below.
- [src/simple_modflow/modflow/mf6/heads_observations.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/heads_observations.py)
  Observation lookup and observation-head table shaping for `HeadsPlus`.
- [src/simple_modflow/modflow/mf6/heads_plotting.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/heads_plotting.py)
  Heads plotting and choropleth presentation helpers for `HeadsPlus`.
- [src/simple_modflow/modflow/mf6/observations.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/observations.py)
  Reusable observation-target layer, currently centered on `HeadTargets` so the
  same target dataset can drive plotting, residual statistics, and PEST setup.
- [src/simple_modflow/modflow/mf6/headsplus.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/headsplus.py)
  Backward-compatible entry point for the heads API.
- [src/simple_modflow/modflow/mf6/archive](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/archive)
  Legacy or archived MF6 code. Avoid unless you intentionally need older logic.
- [src/simple_modflow/modflow/mf6/maybe junk](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/maybe%20junk)
  Explicitly non-primary code. Treat as experimental or deprecated.

#### `modflow/mf6/grid`

Grid construction, Voronoi handling, plotting, and optimization.

- [src/simple_modflow/modflow/mf6/grid/triangle.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/triangle.py)
  Main `TriangleGrid` implementation. Domain setup, feature refinement, mesh
  cleanup, optimization, and high-level `build_mesh(...)` profiles live here.
- [src/simple_modflow/modflow/mf6/grid/voronoi.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/voronoi.py)
  `VoronoiGridPlus` implementation and grid utilities.
- [src/simple_modflow/modflow/mf6/grid/plotting.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/plotting.py)
  Grid-specific plotting helpers.
- [src/simple_modflow/modflow/mf6/grid/mesh_quality.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/mesh_quality.py)
  Mesh quality metrics and reporting.
- [src/simple_modflow/modflow/mf6/grid/geometry_cleanup.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/geometry_cleanup.py)
  Geometry cleanup and selective resampling helpers.
- [src/simple_modflow/modflow/mf6/grid/seed_optimization.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/seed_optimization.py)
  Seed-relaxation and CVT-style optimization helpers used by `TriangleGrid`.

#### `modflow/mf6/simulation`

Shared model-building and model-access infrastructure.

- [src/simple_modflow/modflow/mf6/simulation/base.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/base.py)
  `SimulationBase`, the main model object. If you are looking for the central
  model API, start here.
- [src/simple_modflow/modflow/mf6/simulation/packages.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/packages.py)
  High-level package wrappers that apply `simple_modflow` defaults and
  optionally capture reusable package artifacts.
- [src/simple_modflow/modflow/mf6/simulation/discretization.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/discretization.py)
  Discretization helpers like `DisvGrid`, `DisuGrid`, and temporal
  discretization wrappers.
- [src/simple_modflow/modflow/mf6/simulation/accessors.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/accessors.py)
  Shared accessor constructors used by `SimulationBase`, including heads,
  budgets, inputs, outputs, and plotting access points.
- [src/simple_modflow/modflow/mf6/package_explorer.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/package_explorer.py)
  Preferred normalized package-input exploration layer behind
  `model.packages...`, including normalized input tables, budget-backed result
  tables, and stage-result explorers.
- [src/simple_modflow/modflow/mf6/simulation/indexing.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/indexing.py)
  Cell, node, layer, and time indexing helpers.
- [src/simple_modflow/modflow/mf6/simulation/regions.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/regions.py)
  Region and region-group registry support.
- [src/simple_modflow/modflow/mf6/simulation/runtime.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/runtime.py)
  Shared model write/run execution helper.
- [src/simple_modflow/modflow/mf6/pest](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/pest)
  The first reusable pyEMU/PEST layer. This is where `PestProject`,
  calibration parameter specs, support-file builders, and forward-run helpers
  now live.

#### `modflow/utils`

Low-level support utilities used across the package.

- [src/simple_modflow/modflow/utils/outputs.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/utils/outputs.py)
  Output-specific helper objects such as UZF, LAK, and SFR output accessors.
- [src/simple_modflow/modflow/utils/datatypes](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/utils/datatypes)
  Reusable datatype and reader helpers.
- [src/simple_modflow/modflow/utils/validators.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/utils/validators.py)
  Common validation helpers.

### `project`

Project organization, reusable package artifacts, run loading, run exploration,
and grouped-model comparison live here.

- [src/simple_modflow/project/catalog.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/catalog.py)
  `ProjectCatalog`: project registry for model specs, runs, and package
  artifacts.
- [src/simple_modflow/project/specs.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/specs.py)
  Data models for `ModelSpec`, `RunSpec`, and `RunRecord`.
- [src/simple_modflow/project/components.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/components.py)
  Reusable package artifact capture, validation, derivation, and application.
- [src/simple_modflow/project/run_model.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/run_model.py)
  `LoadedMf6Run`: lazy file-backed model loader for existing MF6 directories.
- [src/simple_modflow/project/explorer.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/explorer.py)
  `RunExplorer`: notebook-friendly browsing over a directory of MF6 runs.
- [src/simple_modflow/project/model_group.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/model_group.py)
  `ModelGroup`: grouped multi-model API for comparing heads, budgets, inputs,
  package tables, and selected outputs.
- [src/simple_modflow/project/loaders.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/loaders.py)
  Load helpers used by catalogs and explorers.
- [src/simple_modflow/project/discovery.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/discovery.py)
  Run-discovery logic for plain MF6 directory trees.
- [src/simple_modflow/project/compare.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/compare.py)
  Catalog-level run comparison helpers.
- [src/simple_modflow/project/manifest_io.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/manifest_io.py)
  TOML manifest read/write helpers for catalogs, specs, runs, and artifacts.
- [src/simple_modflow/project/helpers.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/helpers.py)
  Miscellaneous project-layer support functions.

## `examples`

- [examples/mf6](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6)
  Current MF6 examples and notebooks.
- [examples/mf6/notebooks](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks)
  Practical workflow notebooks, including:
  - Cumberland grid workflows
  - first-slice PEST/pyEMU calibration workflows
  - run explorer workflows
  - package artifact workflows
  - optimized related-run workflows

## `tests`

- [tests/test_project_catalog.py](C:/Users/lukem/Python/Projects/simple_modflow/tests/test_project_catalog.py)
  Main workflow regression suite for catalogs, loaded runs, lazy loading,
  artifacts, model groups, and related-run workflows.
- [tests/test_mf6_refactor_smoke.py](C:/Users/lukem/Python/Projects/simple_modflow/tests/test_mf6_refactor_smoke.py)
  Broad MF6 smoke/integration coverage, including notebook presence and larger
  workflow tests.
- [tests/test_mesh_optimization.py](C:/Users/lukem/Python/Projects/simple_modflow/tests/test_mesh_optimization.py)
  Mesh cleanup, quality reporting, and optimization regression coverage.
- [tests/test_mf6_pest.py](C:/Users/lukem/Python/Projects/simple_modflow/tests/test_mf6_pest.py)
  First-slice calibration coverage for `HeadTargets`, bounds logic,
  `PestProject`, and forward-run reapplication of `K` and `DRN` parameters.

## How To Find Code Quickly

- If you are building a new model: start in [simulation/base.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/base.py:1) and [simulation/packages.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/packages.py:1).
- If you are building package data from shapefiles or geopackages: prefer the
  vector-builder classes `DRNFromVector`, `GHBFromVector`, `CHDFromVector`,
  `RCHFromVector`, and `KFromVector`.
- If you are working on Voronoi or Triangle grids: start in [grid/triangle.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/triangle.py:1) and [grid/voronoi.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/grid/voronoi.py:1).
- If you are diagnosing stream/lake/mover coupling problems: start in [surface_water_validation.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/surface_water_validation.py:1), then check [sfr.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/sfr.py:1) and [lakes.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/lakes.py:1).
- If you are reopening existing runs: start in [project/run_model.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/run_model.py:1).
- If you are comparing runs: start in [project/model_group.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/model_group.py:1).
- If you are working on the preferred package input tables and choropleths:
  start in [package_explorer.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/package_explorer.py:1),
  then check [simulation/accessors.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/simulation/accessors.py:1)
  and [project/model_group.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/model_group.py:1).
- If you are working on reusable package snapshots: start in [project/components.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/project/components.py:1).
- If you are working on heads logic: start in [headsplus.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/headsplus.py:1), then look at [heads_observations.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/heads_observations.py:1) and [heads_plotting.py](C:/Users/lukem/Python/Projects/simple_modflow/src/simple_modflow/modflow/mf6/heads_plotting.py:1).

## Notes On Legacy Areas

- `archive/` and `maybe junk/` under `mf6` are not part of the preferred
  workflow path.
- Older workflow helpers still exist in some modules because backward
  compatibility matters, but new development should prefer the catalog/run/group
  architecture plus `SimulationBase`.
