# MF6 Examples

This folder holds notebooks, sample model outputs, and one-off artifacts that were previously stored under
`src/myflopy/modflow/mf6/`.

Contents:

- `sample_model_output/`: sample rasters and triangle outputs used by the `simplemodel.py` demo block
- `artifacts/`: ad hoc outputs, temporary plots, logs, and other non-library files
- `notebooks/`: exploratory notebooks moved out of the importable package tree

Suggested starting points:

The preferred GeoPackage-first API is `GeoPackageSource`. It maps one
GeoPackage layer onto a `ModelContext` grid and returns reusable package specs
through `.chd()`, `.ghb()`, `.drn()`, `.wel()`, and `.rch()`, or conductivity
data through `.k_array()`.

- `notebooks/simple_model_workflow.ipynb`: minimal `SimpleModelConfig` to `SimulationSpec` and `Project.run()` example
- `notebooks/triangle_voronoi_simplemodel_workflow.ipynb`: build a mesh with `TriangleGrid`, convert to `VoronoiGridPlus`, then run a simple model
- `notebooks/feature_rich_model_workflow.ipynb`: build and run a small but more realistic MF6 model with drains, GHB, recharge/UZF, lake, stream, and region groups
- `notebooks/refined_feature_rich_model_workflow.ipynb`: build and run a denser `TriangleGrid -> VoronoiGridPlus` model with refinement zones, then add the same boundary and package workflow on top of a more realistic mesh
- `notebooks/code_geometry_refined_model_workflow.ipynb`: define the domain, lake, stream, drain, GHB, and recharge zones directly in code with `shapely`, use those same geometries for grid refinement and model regions, and write temporary GIS layers only where today’s package-builder APIs still expect files
- `notebooks/cumberland_mesh_optimization_workflow.ipynb`: use the actual Cumberland mesh-input pattern from the pre-AlgoMesh workflow, compare a baseline Triangle/Voronoi mesh with a conservative optimized build, and optionally compare both against the imported AlgoMesh `.disu` grid
- `notebooks/cumberland_cvt_diagnostics_workflow.ipynb`: focused diagnostics notebook for the new constrained full-seed CVT/Lloyd implementation, comparing `baseline`, `cleanup_only`, `experimental_cvt`, and optional `AlgoMesh DISU`
- `import_existing_runs_workflow.py`: discover native MF6 workspaces through `Project.discover_native_runs()`
- `notebooks/refined_end_to_end_preferred_api_workflow.ipynb`: historical end-to-end workflow that will be migrated to `GeoPackageSource`
- `notebooks/pest_first_slice_workflow.ipynb`: the first reusable PEST/pyEMU workflow example, showing GIS-defined `K` and `DRN` inputs, `HeadTargets`, `PestProject` template generation, and a forward run that reapplies pilot-point `K` multipliers plus drain elevation/conductance changes
- `notebooks/targets_api_quickstart.ipynb`: short focused notebook showing the flexible target API, including `HeadTargets`, `LakeStageTargets`, `model.targets...`, plotting helpers, FloPy observation object generation, and the PEST/reopen touchpoints
- `notebooks/master_large_model_visualization_prt_workflow.ipynb`: master integration notebook that defaults to a 10,000-cell, four-layer DISV/Voronoi model with CHD, GHB, DRN, UZF, SFR, two lakes, MVR, canonical observations, Matplotlib/Plotly standalone result animations, and MF6 PRT/PyVista review
- `master_large_model_review.py`: bounded command-line review harness for the same 10,000-cell model, with timed model execution, resumable selected-frame exports, optional bounded Plotly and PRT/PyVista review, and a JSON performance/result summary
- `../../docs/interactive_visualization_and_prt.md`: canonical standalone Matplotlib slider, Plotly animation, MF6 PRT, and PyVista 3D particle-viewing workflow
- `../../docs/princeton_2026_visualization_review.md`: detailed mapping from the Princeton 2026 FloPy training notebooks into the preferred `myflopy` API
- `notebooks/cumberland_predev_pest_first_slice_workflow.ipynb`: Cumberland-specific version of the first-slice calibration workflow, following the pre-development geometry adjustments from `model_pit_pre_excavation.py` but using one steady-state period and the current `ks_v5.gpkg` / `drn.gpkg` inputs
- `notebooks/cumberland_observed_snapshot_pest_first_slice_workflow.ipynb`: Cumberland-specific first-slice workflow that swaps the pseudo-targets for a real observed-head snapshot from `calib_observations.xlsx`, while keeping the model side to a one-period steady-state calibration example
- `cumberland_forward_mp3du.py`: command-line Cumberland MP3DU forward-tracking harness that exercises the stable `ParticleTrackingInput` / `prepare_particle_tracking` / `run_particle_tracking` API, writes a timestamped MP3DU workspace under `artifacts`, and captures start-cell plus endpoint diagnostics
- `cumberland_steady_snapshot_pest_first_slice.py`: command-line Cumberland steady-state snapshot calibration example using the current GIS inputs, one representative observed-head snapshot plus low-weight supplemental wells, timestamped artifact workspaces, and the current first-slice `PestProject` API, including `--fast-mode`, `--k-only`, optional local parallel workers, and any desired `--noptmax`
- `cumberland_transient_observed_pest_first_slice.py`: command-line Cumberland transient calibration example with one initial steady-state period plus 12 monthly transient periods, real observed heads, low-weight one-time supplemental wells, timestamped workspaces, and the current first-slice `PestProject` API, including `--fast-mode`, `--k-only`, and local parallel worker options
- `import_existing_runs_workflow.py`: command-line version of the archive import and comparison workflow

Reference handout:

- `../docs/myflopy_api_pamphlet.pdf`: a short visually organized API pamphlet covering the main workflows by topic/page
- `../docs/mp3du_quickstart.md`: focused quickstart for the supported MP3DU particle-tracking API
- `../../docs/local_run_clutter.md`: where example runs, pytest temp workspaces, and scratch MP3DU folders are stored locally
# Canonical master notebooks

The preferred end-to-end learning path starts at
`../../docs/CANONICAL_MASTER_NOTEBOOKS.md`. Its five notebooks all use
`myflopy.build_canonical_model()` so model construction, observations,
visualization, PRT, parallel splitting, PEST, and completed-run review remain
aligned.
