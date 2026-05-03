# MF6 Examples

This folder holds notebooks, sample model outputs, and one-off artifacts that were previously stored under
`src/simple_modflow/modflow/mf6/`.

Contents:

- `sample_model_output/`: sample rasters and triangle outputs used by the `simplemodel.py` demo block
- `artifacts/`: ad hoc outputs, temporary plots, logs, and other non-library files
- `notebooks/`: exploratory notebooks moved out of the importable package tree

Suggested starting points:

Preferred GIS/vector builder names in the current API:
- `DRNFromVector`
- `GHBFromVector`
- `CHDFromVector`
- `RCHFromVector`
- `KFromVector`

Older builder names like `DRN`, `GHB`, and `RechargeFromShp` still work for
compatibility, but the names above are the intended long-term surface.

- `notebooks/simple_model_workflow.ipynb`: minimal config-first `SimpleModel` example
- `notebooks/triangle_voronoi_simplemodel_workflow.ipynb`: build a mesh with `TriangleGrid`, convert to `VoronoiGridPlus`, then run a simple model
- `notebooks/feature_rich_model_workflow.ipynb`: build and run a small but more realistic MF6 model with drains, GHB, recharge/UZF, lake, stream, and region groups
- `notebooks/refined_feature_rich_model_workflow.ipynb`: build and run a denser `TriangleGrid -> VoronoiGridPlus` model with refinement zones, then add the same boundary and package workflow on top of a more realistic mesh
- `notebooks/code_geometry_refined_model_workflow.ipynb`: define the domain, lake, stream, drain, GHB, and recharge zones directly in code with `shapely`, use those same geometries for grid refinement and model regions, and write temporary GIS layers only where today’s package-builder APIs still expect files
- `notebooks/cumberland_mesh_optimization_workflow.ipynb`: use the actual Cumberland mesh-input pattern from the pre-AlgoMesh workflow, compare a baseline Triangle/Voronoi mesh with a conservative optimized build, and optionally compare both against the imported AlgoMesh `.disu` grid
- `notebooks/cumberland_cvt_diagnostics_workflow.ipynb`: focused diagnostics notebook for the new constrained full-seed CVT/Lloyd implementation, comparing `baseline`, `cleanup_only`, `experimental_cvt`, and optional `AlgoMesh DISU`
- `notebooks/import_existing_runs_workflow.ipynb`: scan an archive of existing MF6 runs, browse it naturally with `RunExplorer`, inspect lightweight run summaries, optionally import into a `ProjectCatalog`, and compare runs lazily
- `notebooks/project_artifact_workflow.ipynb`: the modern project-managed workflow for related runs, showing inline artifact creation, preflight validation, auto-apply on `run_simulation()`, and derived artifacts with lineage
- `notebooks/optimized_related_run_artifact_workflow.ipynb`: a fuller related-run example on a modest optimized Voronoi model, showing source-run artifact capture, derived artifacts for scenario changes, validation, auto-apply, and result comparison on the same mesh
- `notebooks/refined_end_to_end_preferred_api_workflow.ipynb`: a robust end-to-end example that writes real geopackage inputs on the fly, builds an optimized refined Voronoi grid, validates the coupled surface-water inputs, and uses the preferred `*FromVector` builders plus UZF/LAK/SFR with MVR routing on a stable two-period model
- `import_existing_runs_workflow.py`: command-line version of the archive import and comparison workflow
