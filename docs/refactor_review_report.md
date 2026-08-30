# Refactor Review Report

Date: 2026-06-12

## Executive Verdict

The refactor has **not gone off the rails**. The new architecture has a coherent
center:

- `SimulationBase` is the live-model facade.
- `grid`, `simulation`, `pest`, and `project` are meaningful ownership
  boundaries.
- `LoadedMf6Run`, `RunExplorer`, and `ModelGroup` provide a consistent path from
  one live model to reopened runs and multi-run comparison.
- Model-bound namespaces such as `model.packages`, `model.targets`,
  `model.visualize`, `model.particle_tracking`, and `model.parallel` give the
  API an understandable overall shape.
- The test suite provides broad coverage of the new behavior.

The main risk is no longer the overall direction. The main risk is that several
successful feature areas have accumulated too much code behind a small number
of facades. That makes ownership harder to see, increases the cost of review,
and encourages further additions to already oversized files.

The recommended next phase is therefore **consolidation, not another redesign**:

1. Freeze the preferred public API long enough to document and test it.
2. Split the four largest implementation modules by responsibility without
   changing their user-facing namespaces.
3. make compatibility names visibly secondary.
4. Remove generated outputs and executed notebook payloads from version control.
5. Split the largest test files along the same ownership boundaries as the
   production code.

## Scope Of This Review

This review covers:

- The current `codex-refactor-mf6-core` branch relative to `main`.
- The current working tree as of 2026-06-12.
- Importable code under `src/myflopy`.
- The new tests, docs, examples, and notebook organization.
- The actual lazy exports in `myflopy.__init__` and
  `myflopy.modflow.mf6.__init__`.

The branch contains four commits after `main`:

- `01fa7cb` - clean commit without large files
- `784bd1e` - PEST demo working
- `6250b5f` - repository cleanup and MP3DU fixes
- `fb24b54` - example cleanup

There are also current uncommitted changes in canonical examples, observation
comparison normalization, PRT stopping controls, tests, and notebooks.

## Scale

Approximate current scale:

| Area | Size |
|---|---:|
| `src/myflopy` | about 45,300 Python lines |
| `src/myflopy/modflow` | about 39,600 Python lines |
| `src/myflopy/project` | about 5,400 Python lines |
| `tests` | about 12,000 Python lines |
| top-level `myflopy` lazy exports | 90 names |
| `myflopy.modflow.mf6` exports | 94 names |
| collected tests | 235 |

The branch adds a large amount of useful implementation and test code, but the
largest raw change is examples and notebook output. Relative to `main`,
`examples` accounts for more than 700,000 added text lines, mostly because
notebook outputs are stored inline.

## What Is Working Well

### 1. The architecture now has a recognizable center

The preferred path is understandable:

```text
TriangleGrid / VoronoiGridPlus
              |
              v
        SimulationBase
     /       |        \
packages   targets   outputs/visualization
     \       |        /
        ProjectCatalog
              |
      LoadedMf6Run / RunExplorer
              |
          ModelGroup
```

This is substantially clearer than a flat collection of MF6 helper modules.

### 2. Live models and reopened models share useful concepts

`SimulationBase`, `LoadedMf6Run`, and `ModelGroup` converge on heads, budgets,
package exploration, outputs, and summaries. That is a strong design choice
because notebook workflows do not need a completely different mental model
after a run has been written to disk.

### 3. The model-bound namespaces are a good API direction

These access points are memorable and leave room for internal reorganization:

- `model.packages`
- `model.targets`
- `model.visualize`
- `model.particle_tracking`
- `model.parallel`
- `model.outputs`

The same is true for grouped access such as `group.packages` and `group.hds`.

### 4. The refactor created real regression coverage

The suite exercises:

- public and compatibility imports
- grid construction and optimization
- live models, loaded runs, catalogs, artifacts, and model groups
- normalized package input and result exploration
- targets and PEST integration
- PRT and MP3DU particle tracking
- interactive plotting
- canonical-model and parallel workflows

On 2026-06-12, the full suite in `mf-env` produced:

```text
234 passed, 1 failed, 125 warnings
```

The one failure was the canonical MPI test. Intel MPI reported that it could
not create its bootstrap stdout pipe. The serial canonical model run succeeded,
so this is an environment/launcher failure rather than evidence of a model
calculation regression.

## Main Findings

### Priority 1: Large facade modules now hide too many responsibilities

The clearest readability risk is concentrated in four files:

| File | Approximate lines | Concern |
|---|---:|---|
| `modflow/mf6/package_explorer.py` | 5,084 | Tables, normalization, maps, plots, namespaces, registry, and package-specific behavior in one module |
| `modflow/mf6/observations.py` | 2,608 | Five target types, model-bound wrappers, plotting facades, registry, and shared normalization |
| `project/model_group.py` | 2,416 | Group alignment, comparisons, maps, plots, and every grouped namespace |
| `modflow/mf6/grid/triangle.py` | 1,630 | Domain definitions, feature registration, cleanup, build, quality, and optimization |

These files are not unreadable because the code is inherently bad. They are
hard to review because they each contain several distinct ownership areas.

Recommended target boundaries:

- Keep `ModelPackages` in `package_explorer.py`, but move table builders,
  package-specific explorers, and visualization helpers into focused modules.
- Keep target registry and public target classes in `observations.py`, but move
  shared series logic, plotting wrappers, and each target family into focused
  modules.
- Keep `ModelGroup` and its public namespace graph in `model_group.py`, but move
  group alignment, map rendering, and package-specific grouped accessors out.
- Keep `TriangleGrid` as the public facade, but move domain/feature specs,
  cleanup orchestration, and optimization orchestration into collaborators.

Do not change the user-facing namespace while making these splits.

### Priority 1: The public API is broad and has one important naming collision

The top-level package exports 90 names and the MF6 package exports 94. All
currently resolve successfully in `mf-env`, but the surface is large enough
that it is difficult to tell which names are preferred, advanced, or retained
only for compatibility.

The clearest ambiguity is `GHB`:

```python
from myflopy.modflow.mf6 import GHB
# legacy/vector builder from myflopy.modflow.mf6.ghb

from myflopy.modflow.mf6.simulation import GHB
# package wrapper from myflopy.modflow.mf6.simulation.packages
```

The preferred vector builder is `GHBFromVector`, but the short name `GHB`
means different things at different import levels. `CHD` and `Drains` avoid
this exact collision because their preferred vector builders have distinct
names.

Recommendation:

- Reserve short package names such as `GHB`, `CHD`, `LAK`, and `UZF` for
  package attachment wrappers.
- Keep vector builders named `GHBFromVector`, `CHDFromVector`, and so on.
- Treat `myflopy.modflow.mf6.ghb.GHB` as a compatibility import and
  remove it from the preferred `mf6` export surface in a future deprecation
  cycle.

> **RESOLVED (2026-07-16, implementation plan 3.2):** `GHBFromVector` /
> `DRNFromVector` are now the real class names; the bare `GHB`/`DRN`
> spellings warn and are hidden from completion (`docs/deprecation_policy.md`).
> The simulation-layer `GHB` package wrapper keeps its name.

### Priority 1: Generated examples and outputs dominate the repository

The branch tracks:

- 70 files under `examples/mf6/notebooks`
- 40 files under `examples/mf6/notebooks/split_profile`
- executed notebooks
- a built wheel under `.build-smoke`
- generated `src/myflopy.egg-info`
- IDE metadata under `.idea`

This makes review diffs noisy and gives examples more visual weight than the
library itself. It also makes it harder to tell which notebook is authoritative.

Recommendation:

- Keep source notebooks, small audit JSON, and intentionally curated fixtures.
- Do not track executed notebooks unless a specific rendered output is a
  required artifact.
- Do not track split model workspaces or profiling outputs.
- Do not track wheels, `egg-info`, logs, or IDE metadata.
- Store large demonstration outputs outside the repository or publish them as
  release artifacts.

### Priority 2: Documentation has drifted as features were appended

`docs/codebase_structure.md` is a useful ownership map and
`docs/preferred_api.md` contains substantial workflow guidance. However,
`preferred_api.md` has become a long append-only document.

Concrete issues:

- The interactive visualization code fence opened around line 940 is not closed
  before the parallel-model heading, so the parallel section is structurally
  inside the code example.
- Additional visualization calls appear after the parallel section instead of
  beside the earlier visualization calls.
- Surface-water mapping notes appear after `Where To Look Next`, which should be
  the end of the document.
- The distinction between top-level preferred imports, MF6 imports, and direct
  implementation-module imports is not explicit enough.

Recommendation:

- Make one short API reference the source of truth.
- Keep `preferred_api.md` workflow-oriented and move detailed method inventories
  to reference pages.
- Add a docs check for unbalanced fenced code blocks and invalid local links.

### Priority 2: Tests are broad but repeat the production-code concentration

The largest test files are:

| File | Approximate lines |
|---|---:|
| `tests/test_project_catalog.py` | 4,656 |
| `tests/test_mf6_refactor_smoke.py` | 2,801 |
| `tests/test_mf6_pest.py` | 1,383 |
| `tests/test_interactive_plotting.py` | 818 |

This gives strong regression protection, but it makes ownership and failure
triage harder. `test_mf6_refactor_smoke.py` in particular now tests many
unrelated architectural areas.

Recommendation:

- Split tests by the same boundaries proposed for production modules.
- Keep one small public-import smoke file.
- Mark external executable, canonical, MPI, and slow integration tests
  separately.
- Make the default local suite skip MPI unless the launcher passes a small
  executable health check, not merely an executable-presence check.

### Priority 2: Compatibility code is still mixed into active code paths

The package correctly retains compatibility surfaces, but the boundary is not
always visually obvious:

- `mf6/archive` is clearly marked.
- `mf6/maybe junk` is clearly non-primary but should not remain in an importable
  package long term.
- `mfsimbase.py` and `voronoiplus.py` are compatibility facades.
- Old builder names coexist with preferred `*FromVector` names.
- `budget.py` contains compatibility-oriented wrappers beside newer package
  exploration behavior.

Recommendation:

- Introduce a documented compatibility policy.
- Mark compatibility exports in docstrings and API docs.
- Add deprecation warnings only when there is a clear replacement and migration
  path.
- Move `maybe junk` out of `src` immediately; archive or delete it after review.

### Priority 2: Warning output identifies a few maintenance tasks

The full suite emitted 125 warnings:

- FloPy internal/deprecation warnings.
- `trame_vtk` deprecation warnings.
- 50 pandas fragmentation warnings from observation-frame construction.

The observation warning is locally actionable. Repeated column insertion should
be replaced with one construction or concatenation step. The dependency
warnings should be tracked, but they do not currently indicate failing behavior.

## Recommended Review Strategy

### Phase 0: Establish the review baseline

Before more feature work:

- Commit or intentionally set aside the current working-tree changes.
- Define the exact branch comparison point.
- Remove generated repository artifacts so future diffs show code, not output.
- Record the current full-suite result and the known MPI launcher failure.

Exit condition: reviewers can inspect a clean diff containing only intentional
source, tests, docs, and curated examples.

### Phase 1: Freeze and classify the API

Classify every exported name as one of:

- **Preferred**: documented for ordinary users.
- **Advanced**: supported but normally imported from a subpackage.
- **Compatibility**: retained for older code.
- **Internal**: not exported and free to change.

Recommended import policy:

```python
import myflopy as mf
```

Use top-level `mf` for common workflow objects. Use focused subpackages for
package wrappers and advanced utilities:

```python
from myflopy.modflow.mf6.simulation import (
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.packages import (
    InitialConditions,
    KFlow,
    Storage,
    OutputControl,
)
```

Exit condition: every public export has a category, a one-line description, and
at least one import test.

### Phase 2: Split large modules behind stable facades

Suggested sequence:

1. `package_explorer.py`
2. `observations.py`
3. `model_group.py`
4. `grid/triangle.py`

For each split:

- Preserve public class names and user-facing property paths.
- Move behavior in small mechanical steps.
- Add focused tests before moving behavior.
- Avoid simultaneous API renaming.

Exit condition: no active module exceeds roughly 1,500 lines without a written
reason, and each module has one clear responsibility.

### Phase 3: Align tests and docs with ownership

- Split monolithic test files.
- Add API-surface tests for preferred exports.
- Add documentation structure checks.
- Ensure every preferred workflow has one canonical, unexecuted example.
- Keep environment-dependent integration tests explicitly marked.

Exit condition: a developer can identify the correct source file, test file, and
example from the same subsystem name.

### Phase 4: Retire or isolate legacy code

- Move `maybe junk` out of the package.
- Decide whether `archive` belongs in source distribution.
- Document or deprecate old builder names.
- Remove compatibility paths only after measuring internal/example usage.

Exit condition: ordinary source browsing does not mix preferred and abandoned
implementation paths.

## Folder And File Map

### Repository root

| Path | Responsibility |
|---|---|
| `README.md` | Short project introduction and links to developer docs |
| `pyproject.toml` | Packaging, dependencies, optional extras, and pytest configuration |
| `docs/` | Architecture, API, workflow, and integration documentation |
| `examples/mf6/` | Scripts, source notebooks, and currently some generated model artifacts |
| `scripts/` | Local test helpers |
| `src/myflopy/` | Importable package |
| `tests/` | Unit, smoke, workflow, canonical, and integration tests |

### `src/myflopy`

| File/folder | Responsibility |
|---|---|
| `__init__.py` | Top-level lazy public API |
| `py.typed` | Marks the package as typed |
| `project/` | Catalogs, manifests, run loading, artifacts, and grouped comparison |
| `modflow/` | MODFLOW-oriented implementation |

### `src/myflopy/project`

| File | Responsibility |
|---|---|
| `__init__.py` | Project-layer public exports |
| `specs.py` | Internal catalog metadata plus `RunSpec` and `RunRecord` |
| `catalog.py` | `ProjectCatalog`; registers models, runs, and package artifacts |
| `components.py` | Package artifact capture, compatibility checks, derivation, and application |
| `manifest_io.py` | TOML serialization and deserialization |
| `discovery.py` | Finds existing MF6 run directories |
| `loaders.py` | `RunLoader` construction and loading |
| `run_model.py` | `LoadedMf6Run`, the lazy file-backed model facade |
| `explorer.py` | `RunExplorer` and `explore_runs` |
| `model_group.py` | `ModelGroup` and all grouped result/input namespace facades |
| `compare.py` | Catalog-level `RunComparison` |
| `helpers.py` | Archive import and discovery-summary helpers |

### `src/myflopy/modflow/mf6/simulation`

| File | Responsibility |
|---|---|
| `__init__.py` | Simulation-layer exports |
| `base.py` | `SimulationBase`, the central live-model facade |
| `discretization.py` | `DisvGrid`, `DisuGrid`, and `TemporalDiscretization` |
| `packages.py` | Thin package-attachment wrappers and artifact hooks |
| `accessors.py` | Heads, budgets, inputs, outputs, choropleth, and cross-section constructors |
| `indexing.py` | Layer/cell/node/time indexing helpers |
| `regions.py` | Model regions, region groups, and registry |
| `runtime.py` | Shared simulation write/run helper |

### `src/myflopy/modflow/mf6/grid`

| File | Responsibility |
|---|---|
| `__init__.py` | Grid-layer exports |
| `triangle.py` | `TriangleGrid`, mesh profiles, build orchestration, cleanup, and optimization facade |
| `voronoi.py` | `VoronoiGridPlus` and Voronoi-grid facade |
| `connectivity.py` | DISU connectivity construction |
| `geometry.py` | Geometry, adjacency, centroid, overlap, and refinement helpers |
| `geometry_cleanup.py` | Geometry cleanup and resampling |
| `selection.py` | Grid cell selection and boundary-cell helpers |
| `surfaces.py` | Grid surface and raster helpers |
| `plotting.py` | Grid plots, maps, and `GridSection` |
| `mesh_quality.py` | Mesh quality metrics and report |
| `seed_optimization.py` | CVT-style seed relaxation |
| `helpers.py` | Small shared grid helpers |

### `src/myflopy/modflow/mf6/pest`

| File | Responsibility |
|---|---|
| `__init__.py` | Supported PEST public exports |
| `specs.py` | Parameter and observation specification data classes |
| `project.py` | `PestProject` orchestration |
| `parameters.py` | Parameter support files and registration |
| `observations.py` | PEST observation-file preparation |
| `forward_run.py` | Forward-run parameter reapplication and execution |
| `results.py` | `PestRunResults`, `PestRunReview`, and reopening results |
| `gis.py` | GIS parameter-source loading and bounds |
| `geostats.py` | pyEMU geostatistical structure helper |
| `synthetic_demo.py` | Compact synthetic end-to-end demonstration |
| `gold_standard_demo.py` | Full canonical PEST demonstration |

### Active MF6 root modules: model building and inputs

| File | Responsibility |
|---|---|
| `simplemodel.py` | Compact higher-level model builder |
| `boundaries.py` | Shared vector boundary intersection behavior |
| `boundary_support.py` | Normalization and stress-period-data helpers |
| `drn.py` | DRN vector builder; preferred name `DRNFromVector` |
| `ghb.py` | GHB vector builder; preferred name `GHBFromVector` |
| `chd.py` | CHD vector builder |
| `recharge.py` | Recharge vector builder and compatibility names |
| `kflow.py` | Vector-driven hydraulic-conductivity builder |
| `lakes.py` | Lake geometry, connection, package, and period data |
| `sfr.py` | SFR package-data construction and attachment |
| `uzf.py` | UZF package-data preparation |
| `surface_water_validation.py` | LAK/SFR/MVR validation and reports |

### Active MF6 root modules: results and exploration

| File | Responsibility |
|---|---|
| `package_explorer.py` | Normalized package inputs/results and model-bound package namespace |
| `observations.py` | Reusable targets, comparisons, model-bound targets, and target registry |
| `headsplus.py` | Heads reader facade |
| `heads_observations.py` | Heads observation lookup and shaping |
| `heads_plotting.py` | Heads plots and choropleths |
| `budget.py` | Budget reader and compatibility wrappers |
| `budget_tables.py` | Budget table normalization |
| `budget_plotting.py` | Budget plotting |
| `interactive_plotting.py` | Interactive exports and `ModelVisualization` |
| `contour_plotting.py` | Contour generation and plotting |
| `cross_section_plotting.py` | Model cross-section plotting |

### Active MF6 root modules: workflows

| File | Responsibility |
|---|---|
| `parallel.py` | Model splitting, serial validation, MPI execution, and result reconstruction |
| `prt.py` | Native MF6 PRT setup, execution, reopening, and result presentation |
| `canonical.py` | Canonical model contract and diagnostic signals |
| `canonical_example.py` | Canonical model construction |
| `package_explorer.py` | Package-aware input/result workflow facade |

### Other `modflow` areas

| Path | Responsibility |
|---|---|
| `modflow/mp3du/particles.py` | Supported MP3DU input, execution, and results |
| `modflow/mp3du/legacy_prt.py` | Legacy experimental PRT wrappers |
| `modflow/gwt/gwt.py` | Early GWT wrapper classes |
| `modflow/calcs/calibration.py` | Calibration statistics and plotting |
| `modflow/utils/` | Raster, surfaces, inputs, outputs, validators, readers, choropleths, and cross sections |
| `modflow/mf6/archive/` | Archived MF6 implementations |
| `modflow/mf6/maybe junk/` | Non-primary experimental files; should leave `src` |
| `modflow/mf6/mfsimbase.py` | Compatibility import for `SimulationBase` |
| `modflow/mf6/voronoiplus.py` | Compatibility imports for grid classes |

### Tests

| File | Primary coverage |
|---|---|
| `test_mf6_refactor_smoke.py` | Broad refactor imports, compatibility, model and grid smoke behavior |
| `test_project_catalog.py` | Catalogs, runs, artifacts, loaded runs, groups, and package exploration |
| `test_mf6_pest.py` | Targets, PEST specs, support files, results, and forward runs |
| `test_mf6_prt.py` | MF6 PRT setup and results |
| `test_parallel_model.py` | Split workflow, serial validation, MPI, and reconstruction |
| `test_interactive_plotting.py` | Interactive visualization and HTML exports |
| `test_mesh_optimization.py` | Mesh cleanup, quality, and optimization |
| `test_mp3du_particles.py` | Supported MP3DU API and compatibility |
| `test_master_visualization_prt_example.py` | Canonical example and notebook contract |
| `test_package_explorer_registry.py` | Package-explorer registration and namespace consistency |
| `test_grid_array_shapes.py` | Array-shape normalization |
| `test_choro_period_selection.py` | Choropleth period selection |
| `test_manifest_io.py` | TOML manifest behavior |
| `test_gold_standard_pest_demo.py` | Full PEST demo |
| `test_synthetic_pest_demo.py` | Synthetic PEST demo |
| `test_typing_surface.py` | Typing marker and public typing surface |
| `conftest.py` | Shared fixtures and canonical-run setup |

## New API: Brief But Complete Inventory

### Preferred top-level import

```python
import myflopy as mf
```

The top-level public names are grouped below. This is the complete current
top-level surface, excluding `__version__`.

### Live model, compact builders, and grids

- `SimulationBase`
- `SimpleModel`
- `SimpleModelConfig`
- `build_simple_model`
- `TriangleGrid`
- `MeshBuildProfile`
- `VoronoiGridPlus`

### Vector and material builders

- `CHDFromVector`
- `DRNFromVector`
- `GHBFromVector`
- `RCHFromVector`
- `KFromVector`

### Core specification API

- `PackageSpec`
- `ModelSpec`
- `ExchangeSpec`
- `SimulationSpec`

### Project, run, artifact, and group API

- `ProjectCatalog`
- `RunSpec`
- `RunRecord`
- `PackageArtifact`
- `PackageCompatibilityError`
- `LoadedMf6Run`
- `RunLoader`
- `RunExplorer`
- `RunComparison`
- `ModelGroup`
- `DiscoveredRun`
- `discover_existing_runs`
- `summarize_discovered_runs`
- `explore_runs`
- `load_mf6_run`
- `import_run_archive`
- `patch_simulation_plot`

### Observation target API

- `HeadTargets`
- `LakeStageTargets`
- `SfrStageTargets`
- `SfrFlowTargets`
- `DrnFlowTargets`

### PEST API

- `PestProject`
- `PestRunResults`
- `PestRunReview`
- `open_pest_run`
- `VectorParameterSource`
- `ExpGeoStruct`
- `KPilotPointParameter`
- `DrainElevationParameter`
- `DrainConductanceParameter`
- `HeadTargetObservationSpec`
- `LakeStageObservationSpec`
- `SfrStageObservationSpec`
- `SfrFlowObservationSpec`
- `DrnFlowObservationSpec`

### Native PRT and MP3DU API

- `ParticleTracking`
- `PRTProject`
- `PRTReleasePoints`
- `PRTRunResults`
- `open_prt_run`
- `ParticleTrackingInput`
- `prepare_particle_tracking`
- `run_particle_tracking`

### Parallel API

- `ParallelModelWorkflow`
- `ParallelSplitRun`
- `ParallelSplitResults`
- `ParallelEnvironment`
- `ParallelCompatibilityError`

### Visualization API

- `ModelVisualization`
- `ModelMapStyle`
- `FrameExportProgress`
- `StandaloneHtmlSlider`
- `ParticleTrackingScene`
- `plot_model_head_map`
- `plot_particle_pathlines`
- `build_particle_tracking_scene`
- `export_matplotlib_slider_html`
- `export_cross_section_slider_html`
- `export_head_map_slider_html`
- `export_head_layer_mosaic_slider_html`
- `export_particle_tracking_html`

### Surface-water validation API

- `SurfaceWaterValidationIssue`
- `SurfaceWaterValidationReport`
- `validate_surface_water_configuration`

### Canonical example and diagnostics API

- `CANONICAL_MODEL_CONTRACT`
- `CanonicalModelContract`
- `CanonicalModelConfig`
- `build_canonical_model`
- `canonical_head_signals`
- `canonical_feature_signals`
- `canonical_sfr_signals`
- `canonical_partition_mask`

### General readers and package namespaces

- `read_gpkg`
- `read_shp_gpkg`
- `modflow`
- `project`

### Focused simulation API

Import from `myflopy.modflow.mf6.simulation` or
`myflopy.modflow.mf6.simulation.packages`:

- Discretization: `DisvGrid`, `DisuGrid`, `TemporalDiscretization`
- Common package wrappers: `InitialConditions`, `KFlow`, `Storage`,
  `OutputControl`, `Recharge`, `Wells`, `Drains`, `GHB`, `CHD`, `LAK`, `UZF`
- Semantic advanced builders: `LAKBuilder`, `SFRBuilder`, `UZFBuilder`,
  `MVRBuilder`
- Regions: `ModelRegion`, `RegionGroup`, `RegionRegistry`

The old `MVR` wrapper has been replaced by `MVRBuilder`, `Move`, and
`MoverConnection` so examples can reference stable stream/lake endpoints instead
of hand-authored package rows.

### Important compatibility/advanced MF6 names

The `myflopy.modflow.mf6` namespace additionally exposes package
wrappers, older builders, lower-level grid/discretization classes, region
classes, lake/SFR/UZF helpers, and plotting functions.

Notable compatibility or advanced names include:

- `Boundaries`
- `DRN`
- `GHB` from the vector-builder module
- `RechargeFromShp`
- `SFRBuilder`
- `UZFBuilder`
- `LakeTableBuilder`
- `ModelCrossSectionStyle`
- `plot_model_cross_section`
- `irregular_voronoi_grid`

These are supported today, but they should not all carry equal prominence in
the preferred API documentation.

## Model-Bound API

The most important new API is not the flat export list. It is the namespace
graph attached to a model.

### `SimulationBase`

Core lifecycle and summary:

- `model.workspace`
- `model.summary()`
- `model.file_summary()`
- `model.packages.summary()` *(was `model.package_summary()`, deleted 2026-08-30)*
- ~~`model.output_summary()`~~ *(deleted 2026-08-30, no replacement needed)*
- `model.grid_summary()`
- `model.hds.summary()` *(was `model.result_summary()`, deleted 2026-08-30)*
- `model.list_input_files()`
- `model.list_output_files()`
- `model.run_simulation()`
- `model.load_all()`

Grid, time, heads, and budgets:

- `model.modelgrid`
- `model.vor`
- `model.grid_type`
- `model.idomain`
- `model.times`
- `model.per_dates`
- `model.hds`
- `model.all_heads`
- `model.bud(package)`
- `model.budget_cumulative`
- `model.budget_incremental`

Namespaces:

- `model.packages`
- `model.targets`
- `model.outputs`
- `model.visualize`
- `model.particle_tracking`
- `model.parallel`

Regions and groups:

- `model.add_region(...)`
- `model.add_region_from_cells(...)`
- `model.add_region_from_geometry(...)`
- `model.get_region(...)`
- `model.list_regions()`
- `model.add_group(...)`
- `model.list_groups()`
- `model.region_heads(...)`

Artifacts and validation:

- `model.register_package_reference(...)`
- `model.apply_registered_package_artifacts()`
- `model.validate_referenced_package_artifacts()`
- `model.validate_surface_water(...)`

### `model.packages`

Available package namespaces:

- `rch`
- `chd`
- `drn`
- `ghb`
- `wel`
- `uzf`
- `ic`
- `npf`
- `sto`
- `lak`
- `sfr`
- `surface_water`

Common input operations:

- `.inputs.get()`
- `.inputs.summary()`
- `.inputs.map(...)`
- field access such as `.inputs.finf`, `.inputs.pet`, or `.inputs.q`

Common result operations:

- `.results.get()`
- `.results.summary()`
- `.results.wide()`
- `.results.long()`
- `.results.stack()`
- `.results.map(...)`
- `.results.plot_timeseries(...)`

Package-specific additions:

- `model.packages.lak.connections`
- `model.packages.lak.budget`
- `model.packages.lak.results.stage`
- `model.packages.lak.results.stage_change`
- `model.packages.lak.results.q`
- `model.packages.sfr.budget`
- `model.packages.sfr.results.stage`
- `model.packages.sfr.results.q`
- `model.packages.sfr.results.long_profile(...)`
- `model.packages.surface_water.results.q`

### `model.targets`

Registry access:

- `model.targets.keys()`
- `model.targets.summary()`
- `model.targets.heads`
- `model.targets.lake_stage`
- `model.targets.sfr_stage`
- `model.targets.sfr_flow`
- `model.targets.drn_flow`

Common target operations:

- `.get()`
- `.to_long()`
- `.to_wide()`
- `.summary()`
- `.compare()`
- `.residuals()`
- `.stats()`
- `.to_flopy_obs()`
- `.attach_flopy_obs()`

Bound targets add model-aware simulated-series and plotting behavior.

### `model.visualize`

The visualization facade owns:

- head maps
- cross sections
- standalone HTML sliders
- Plotly animation exports
- particle pathline presentation

### `model.particle_tracking`

- `model.particle_tracking.prt(...)`
- `model.particle_tracking.open_prt(...)`
- `model.particle_tracking.mp3du(...)`

### `model.parallel`

- `model.parallel.environment`
- `model.parallel.topology()`
- `model.parallel.split_model(...)`
- `model.parallel.prepare(...)`
- `model.parallel.run(...)`

The returned split run supports validation, serial execution, MPI execution,
head comparison, partition plotting, and reconstructed results.

## Grouped API

`ModelGroup` mirrors the single-model exploration shape where practical:

- `group.summary()`
- `group.hds.get()`
- `group.hds.compare()`
- `group.bud(package).get()`
- `group.bud(package).compare()`
- `group.packages`
- `group.outputs`

Grouped package operations generally provide:

- `get()`
- `compare()`
- `map(...)`
- `compare_map(...)`
- `subplot_map(...)` where appropriate

Available grouped package namespaces include:

- `rch`
- `chd`
- `drn`
- `ghb`
- `wel`
- `uzf`
- `sfr`
- `lak`
- `surface_water`

## API Decisions Still Needed

The following should be decided before calling the refactor complete:

1. Is the top-level package intended to expose package wrappers, or should those
   always come from `mf6.simulation.packages`?
2. Should `myflopy.modflow.mf6.GHB` mean the package wrapper or remain the
   old builder?
3. Which of the 90 top-level names are truly preferred?
4. Are canonical example helpers part of the user API or developer/test API?
5. Are plotting export functions preferred public functions, or implementation
   details behind `model.visualize`?
6. What is the deprecation policy for `DRN`, `GHB`, `RechargeFromShp`,
   `mfsimbase.py`, and `voronoiplus.py`?

## Proposed Definition Of Done

The refactor is ready to leave the cleanup phase when:

- The preferred API fits on one short reference page.
- Every preferred export has a stable import path and import test.
- Compatibility names are labeled and do not collide with preferred meanings.
- `package_explorer.py`, `observations.py`, and `model_group.py` have been split
  behind unchanged public facades.
- No generated wheels, `egg-info`, IDE metadata, executed notebooks, or model
  run workspaces are tracked.
- The default suite is fast enough for ordinary development and excludes
  environment-dependent MPI execution by default.
- Slow, canonical, executable-dependent, and MPI tests are separately
  selectable.
- Documentation has no malformed code fences and points to one canonical
  example for each preferred workflow.
- `archive` and `maybe junk` have explicit disposition.

## Immediate Next Actions

Recommended order:

1. Clean tracked generated artifacts and notebook outputs.
2. Repair and shorten `docs/preferred_api.md`.
3. Create an explicit preferred/advanced/compatibility export manifest.
4. Resolve or document the `GHB` name collision.
5. Split `package_explorer.py` while preserving `model.packages`.
6. Split `observations.py` while preserving `model.targets`.
7. Split `model_group.py` while preserving `group.packages`.
8. Split the large test files to match those ownership boundaries.
9. Add an MPI launcher health check or change MPI execution to an opt-in marker.
10. Address the pandas fragmentation warning in observation-frame construction.

The architectural destination should remain the same. The next work should make
that destination easier to see in the source tree.
> Historical note: this report describes the pre-`Project`/`Run` catalog
> architecture. `ProjectCatalog`, `RunSpec`, `RunRecord`, `RunExplorer`, and
> their TOML manifest layer were removed during the clean `myflopy` refactor.
