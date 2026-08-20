# The myflopy Manual

A Python-first toolkit for building, running, visualizing, and calibrating
MODFLOW 6 models on unstructured (Voronoi) grids, built on top of FloPy.

> **Applies to:** the `myflopy` package on the `myflopy` branch of
> `lmioduszewski/simple_modflow` (package path `src/myflopy/`).
> **Last revised:** 2026-06-24.

---

## About this manual

This is the reference manual for **myflopy**. It is written for hydrogeologists
and modelers who already know MODFLOW 6 and groundwater-flow concepts, but who do
**not** need prior FloPy internals — the manual explains where myflopy sits on top
of FloPy and what it adds.

The manual is organized as a progression: **Foundations → Core Concepts →
Building a Model → Running & Exploring → Calibration → Reference.** If you are new,
read [Chapter 3, *Getting Started*](03_getting_started.md) first — it is a complete
end-to-end overview you can work through in about fifteen minutes. Then dip into
the chapter you need.

Throughout, two example models recur so the code stays concrete:

- a **"from-scratch" valley** — a few-layer model with SFR, LAK, UZF, DRN, GHB,
  and RCH, assembled with the package-first API; and
- the **canonical valley model** — a shared, tested reference model used in the
  results, visualization, and calibration chapters.

### Conventions

| Convention | Meaning |
|---|---|
| `import myflopy as mf` | Assumed at the top of every example; the package-first surface is `mf.*`. |
| `python` code blocks | Runnable as shown unless a comment marks an ellipsis or a placeholder path. |
| **Note** | Background or clarification. |
| **Tip** | A recommended practice or shortcut. |
| **Warning** | A sharp edge — something that will bite if ignored. |
| `file_path:line` | A pointer into the source tree. |

> **Note — the preferred API.** myflopy has a layered design. The **package-first
> API** (`mf.gwf`, `mf.disv`, `mf.ghb`, `mf.sfr`, …) is the canonical, preferred
> way to build models and is what this manual teaches. The older object-oriented
> builder classes (`SFRBuilder`, `LAKBuilder`, …) are the **engine underneath**
> that facade — documented in [Chapter 20](20_api_reference.md) for completeness,
> not as a competing path. See [Chapter 4](04_architecture.md).

---

## Master table of contents

### Part I — Foundations

**1. Introduction** — *[01_introduction.md](01_introduction.md)*
- 1.1 What myflopy is, in one paragraph
- 1.2 Design philosophy: Python-first, declarative, GIS-aware, unstructured-first
- 1.3 Headline capabilities (Voronoi/DISV grids, specs, GIS-driven BCs, multi-physics, visualization, PEST-on-Voronoi, parallel runs)
- 1.4 When to use myflopy vs. raw FloPy vs. modflow-setup
- 1.5 A taste: the smallest complete model (forward reference to Ch. 3)
- 1.6 How the library is organized (the `mf.*` surface; the facade/engine layering)

**2. Installation & Environment** — *[02_installation.md](02_installation.md)*
- 2.1 Python version and core dependencies (FloPy, numpy, pandas, geopandas, shapely, rasterio)
- 2.2 The MODFLOW 6 executable (`mf6`) — install, locate, put on `PATH`
- 2.3 Optional extras: PEST++ (`pestpp-ies`, `pestpp-glm`) + pyEMU for calibration; PyVista + trame for 3-D scenes; Plotly for interactive maps
- 2.4 Installing myflopy (editable from source; `PYTHONPATH=src`)
- 2.5 Verifying the install (`import myflopy as mf`; build-and-run smoke test)
- 2.6 Working directories, units, and the on-disk run layout at a glance

**3. Getting Started (the cheat sheet)** — *[03_getting_started.md](03_getting_started.md)*
- 3.1 The mental model in one diagram
- 3.2 Hello, model — the smallest thing that runs
- 3.3 A real model end-to-end: grid → layers → context → packages → simulation → run
- 3.4 Reading results, plotting, and a first calibration — at a glance
- 3.5 The one-screen API map and "where to go next"

### Part II — Core Concepts

**4. Architecture & the Mental Model** — *[04_architecture.md](04_architecture.md)*
- 4.1 The package-first API as the canonical surface
- 4.2 The facade → engine → atom pattern (`mf.uzf` → `UZFBuilder`; `LayerStack` → `LayerSurfaces` → `Surface`)
- 4.3 The object hierarchy: `Project → SimulationSpec → ModelSpec → ModelContext → packages`
- 4.4 Geometry rides on the *model*, not the project (`model.myflopy_context`)
- 4.5 Immutability and composition (frozen specs; `with_*` copy-on-write)
- 4.6 Model variants & reusable, swappable packages — package/grid libraries, `mf.ref`/`mf.grid_ref`, `derive` + `replace_package`/`replace_grid`, lineage, and what re-maps on a grid swap
- 4.7 The build lifecycle: spec → `build(context)` → FloPy objects → write → run
- 4.8 The `mf.*` import surface (preferred vs. compatibility exports)

**5. The Spec System** — *[05_specs.md](05_specs.md)*
- 5.1 `PackageSpec` — a serializable description of one FloPy package; the three forms (`()`, `.gpkg()`, `.flopy()`)
- 5.2 `ModelSpec` — one GWF/GWT/GWE/PRT model; `mf.gwf/gwt/gwe/prt`
- 5.3 `SimulationSpec` — tdis + ims + models + exchanges; `.model()`, `.with_model()`, `.with_package()`
- 5.4 `ModelContext` — grid, surfaces, dates, domain, metadata; `with_metadata`
- 5.5 `PackageRef` / `GridRef` and the `ref` / `grid_ref` helpers (library reuse)
- 5.6 `PostBuildHook` — attach observations, validation, or extra FloPy objects after the build
- 5.7 `BuiltModel` / `BuiltSimulation` — the products of a build
- 5.8 `SpecBuildContext` — what the build threads through (workspaces, libraries)
- 5.9 Serialization round-trip: `to_dict` / `from_dict`

### Part III — Building a Model

**6. Grids & Discretization** — *[06_grids.md](06_grids.md)*
- 6.1 Why unstructured: Voronoi/DISV for myflopy
- 6.2 `VoronoiGridPlus` — the grid object (`get_disv_gridprops`, `gdf_vorPolys`, `gdf_topbtm`, `ncpl`)
- 6.3 `TriangleGrid` and `MeshBuildProfile` — the triangulation underneath (presets, refinement)
- 6.4 `GridSpec` — declarative grid recipes from GIS (`GridSpec.voronoi(boundary=, refinement=, breaklines=, points=)`)
- 6.5 **Eager vs. deferred** grids — the key seam (`.resolve(workspace)` vs. `mf.gwf().with_grid(...)`), and why the GIS package helpers require eager
- 6.6 CRS handling and reprojection
- 6.7 Inspecting and plotting a grid before you build on it

**7. Layers & Surfaces** — *[07_layers.md](07_layers.md)*
- 7.1 `LayerStack` — the facade: declare a top, then `.add()` named layers by thickness or bottom
- 7.2 `Surface` atoms — `Raster` / `Contours` / `Points` / `Array` / `Flat` (constant) and surface algebra
- 7.3 Sampling onto the grid (area-weighted vs. point-at-centroid)
- 7.4 Top-down reconciliation and pinch-out → `idomain` (active / pass-through / inactive)
- 7.5 `build()` → `LayerBuildResult` (top/botm/idomain/thickness); `build(attach=True)` and `vor.gdf_topbtm`
- 7.6 Quality control: `.qc()` → `LayerQCReport`, `.report()`
- 7.7 Visual checks: cross-sections, 3-D surfaces, VTK export
- 7.8 The engine beneath: `LayerSurfaces` (when and why to drop down)

**8. Flow Packages** — *[08_flow_packages.md](08_flow_packages.md)*
- 8.1 `mf.disv` — discretization from `get_disv_gridprops()` + layer arrays
- 8.2 `mf.ic` — initial heads
- 8.3 `mf.npf` — `k`, `k33`, icelltype, Newton options
- 8.4 `mf.sto` — steady-state vs. transient flags, storage
- 8.5 `mf.oc` — head/budget filerecords and save records
- 8.6 `mf.tdis` — stress periods and time stepping
- 8.7 `mf.ims` — solver complexity, inner/outer iterations, linear acceleration; one solver per model
- 8.8 Common patterns: arrays vs. constants; per-layer vs. per-cell

**9. Boundary Conditions** — *[09_boundary_conditions.md](09_boundary_conditions.md)*
- 9.1 The three forms of every list BC: direct `()`, GIS-driven `.gpkg()`, raw `.flopy()`
- 9.2 `mf.chd`, `mf.ghb`, `mf.drn`, `mf.riv`, `mf.wel`, `mf.rch`, `mf.evt`
- 9.3 GIS-driven BCs in depth: `GeoPackageSource` and `mf.X.gpkg(path, layer=, context=, nper=)`
- 9.4 Field mapping, layer/period base conventions, single-field vs. per-period sequences
- 9.5 Recharge as array vs. list; from rasters/PRISM via `RCHBuilder`
- 9.6 The raw escape hatch: `chd_spec`/`drn_spec`/… and when to use it

**10. Advanced Packages (the surface-water network)** — *[10_advanced_packages.md](10_advanced_packages.md)*
- 10.1 `mf.uzf` — unsaturated zone + ET; cells, hydraulic properties, infiltration/PET
- 10.2 `mf.sfr` — streams from a centerline: reaches, automatic connectivity, `StreamConnection`/`StreamDiversion`, inflow, unit conversions
- 10.3 `mf.lak` — lakes from polygons: stage/bottom/leakance, `LakeOutlet`, `LakeTable`/`LakeTableBuilder` bathymetry, connection modes
- 10.4 `mf.mvr` — the water mover: `mf.Move` + `MoverConnection`
- 10.5 **Semantic movers** — `mf.sfr_connection(sfr, "main_stem")` / `mf.lak_connection(lak, "valley_lake")` instead of raw indices; coordinate-based endpoints
- 10.6 MVR rules: moved packages must be declared and ordered before the mover (validation)
- 10.7 Surface-water validation: `validate_surface_water_configuration`

**11. Multi-Physics: Transport, Energy, Particles** — *[11_multiphysics.md](11_multiphysics.md)*
- 11.1 `mf.gwt` — solute transport; the transport package set
- 11.2 `mf.gwe` — heat/energy transport
- 11.3 `mf.prt` — particle tracking as a model (see also Ch. 15)
- 11.4 Exchanges: `build_gwf_gwt_exchange`, `build_gwf_gwe_exchange`, `build_gwf_prt_exchange`, `build_gwf_gwf_exchange`, and `ExchangeSpec`
- 11.5 Putting coupled models in one `SimulationSpec`; shared `ModelContext`
- 11.6 Solver assignment across coupled models

### Part IV — Running & Exploring

**12. Projects, Runs & Scenarios** — *[12_projects_runs.md](12_projects_runs.md)*
- 12.1 `Project(root, name=)` — the durable workspace and run lifecycle
- 12.2 Libraries: `add_package` / `add_grid` / `add_simulation` for reuse
- 12.3 `prepare_run(name, sim)` → `Run`; building in memory and inspecting `run.model(...)`
- 12.4 `run.execute()` — write + run MF6; `success`, `report`, `workspace`
- 12.5 Scenarios: copy-on-write spec edits, registering and running variants
- 12.6 The on-disk layout (`runs/<name>/`, `_grid/`, manifests) and customizing it
- 12.7 Persistence: `Project.save` / `load`, `Run.reopen`
- 12.8 Importing existing MF6 runs: `load_run` / `Run.load`, `LoadedMf6Run`, `discover_native_runs`

**13. Reading Results** — *[13_results.md](13_results.md)*
- 13.1 The model view returned by `run.model(name)` (live vs. file-backed)
- 13.2 Heads: `model.hds.array(layer=…)`, `model.hds.kstpkper`, `model.all_heads`
- 13.3 Budgets and cell-by-cell flows
- 13.4 Per-package results: SFR/LAK/UZF/DRN output series
- 13.5 Observation **targets** for review: `HeadTargets`, `LakeStageTargets`, `Sfr*Targets`, `DrnFlowTargets` — `compare`, `residuals`, `stats`
- 13.6 The raw FloPy objects when you need them: `run.flopy_model(name)`

**14. Visualization** — *[14_visualization.md](14_visualization.md)*
- 14.1 The `myflopy.viz` front door (`Fig`, `subplots`, `mpl_axes`, `PALETTE`, themes)
- 14.2 Choropleth maps over Voronoi cells (matplotlib `.plot_mpl()` and Plotly)
- 14.3 Cross-sections through the grid
- 14.4 Standalone interactive HTML sliders: head maps, layer mosaics, cross-sections through time (`export_*_slider_html`, `StandaloneHtmlSlider`, `ModelMapStyle`)
- 14.5 The model-bound entry point `ModelVisualization`
- 14.6 Progress callbacks and frame caching for large exports

**15. Particle Tracking** — *[15_particle_tracking.md](15_particle_tracking.md)*
- 15.1 Native MF6 PRT via `PRTProject` (`model.particles.prt(...)`)
- 15.2 Release points: `PRTReleasePoints.from_cells(...)`
- 15.3 Tracking controls: porosity, stop times, weak sinks, tracking times
- 15.4 Results: `PRTRunResults`, pathlines; reopening with `open_prt_run`
- 15.5 Plotting: `model.plot.map(pathlines=...)` (2-D) and `model.plot.grid(pathlines=..., backend="vtk")` (3-D)
- 15.6 The legacy mod-PATH3DU path (`ParticleTrackingInput`) — when and why

**16. Parallel Workflows** — *[16_parallel.md](16_parallel.md)*
- 16.1 `ParallelModelWorkflow` — the entry point; supported topologies
- 16.2 Inspecting topology and the run `environment` (mf6 / mpiexec)
- 16.3 Splitting: `split_model(workspace=, nparts=, mask=)` → `ParallelSplitRun`
- 16.4 Writing, executing serially or under MPI
- 16.5 Reassembling outputs on the original grid (`ParallelSplitResults`)
- 16.6 `ParallelCompatibilityError` and the single-GWF rule

### Part V — Calibration & Uncertainty

**17. PEST / pyEMU Integration** — *[17_pest.md](17_pest.md)*
- 17.1 The big picture: a declarative facade that compiles to native pyEMU `PstFrom`
- 17.2 The front door: `cal = model.pest(name, start_datetime=…)`
- 17.3 Parameterization: `cal.parameterize(target, style=…, bounds=…)` — targets (`k`, `k33`, `recharge`, `chd`, `ghb.*`, `drn.*`, `wel`, `uzf.vks`, `porosity`) and styles (`constant` / `zone` / `grid` / `pilotpoints`)
- 17.4 Geostatistics: `ExpGeoStruct`, anisotropy/bearing/nugget, the grid and pilot-point styles on Voronoi
- 17.5 Observations & forecasts: `cal.observe(...)` / `cal.forecast(...)` with the target classes
- 17.6 Building: `cal.build(pst_name, noptmax=…)`; the injected forward run
- 17.7 History matching with ensembles: `cal.run_ies(reals=, iterations=, workers=)`; `cal.prior` / prior Monte Carlo
- 17.8 Reviewing runs: discovery (`model.pest_runs`), `PestRunHandle.review()` → `IesResults` (phi, ensembles, forecasts, parameter-field maps)
- 17.9 The DBTL loop and "a good fit ≠ a good forecast"
- 17.10 Current parameter/observation coverage and known gaps

### Part VI — Reference & Appendices

**18. Serialization & File Formats** — *[18_serialization.md](18_serialization.md)*
- 18.1 The serialization boundary: what round-trips and what does not
- 18.2 `to_dict` / `from_dict` on specs and sources
- 18.3 Durable source references (`DataSourceSpec` and subclasses; `LiteralSource`)
- 18.4 Recipes on disk and the path to YAML/TOML

**19. Troubleshooting & FAQ** — *[19_troubleshooting.md](19_troubleshooting.md)*
- 19.1 Convergence: Newton, IMS complexity, stiff movers
- 19.2 The eager-vs-deferred grid gotcha (GIS helpers need a concrete grid)
- 19.3 "reach_top is required…" and other surface-elevation errors (`build(attach=True)`)
- 19.4 MVR ordering and undeclared-package errors
- 19.5 Units and unit conversions (length/time conversions on SFR/LAK)
- 19.6 OC filerecords ("HEAD SAVE FILE NOT SPECIFIED")
- 19.7 Missing executables (`mf6`, `pestpp-ies`, `mpiexec`)

**20. API Reference** — *[20_api_reference.md](20_api_reference.md)*
- 20.1 How to read the reference (every public export carries a full docstring)
- 20.2 Package-first surface (`mf.gwf` … `mf.mvr`, `sfr_connection`, `lak_connection`)
- 20.3 Specs (`specs.py`)
- 20.4 Grids & layers (`voronoi.py`, `triangle.py`, `layers.py`, `surfaces.py`)
- 20.5 Sources & GIS (`sources.py`, `geopackage.py`)
- 20.6 Observations & PEST (`observations/`, `pest/`)
- 20.7 Visualization, PRT, parallel
- 20.8 The legacy OO engine (`SFRBuilder`, `LAKBuilder`, `UZFBuilder`, `RCHBuilder`, `MVRBuilder`, `SimulationBase`)

**21. Comparison: myflopy vs. FloPy vs. modflow-setup** — *[21_comparison.md](21_comparison.md)*
- 21.1 What each tool is for
- 21.2 Where myflopy wins (unstructured, transport/energy/PRT, viz, parallel, GIS-driven BCs, SFR-from-centerline, PEST on Voronoi)
- 21.3 Where modflow-setup wins (single-file YAML, NHDPlus SFR, LGR)
- 21.4 Interop: myflopy is FloPy all the way down

**Appendix A — The Canonical Valley Model** — *[A_canonical_model.md](A_canonical_model.md)*
- A.1 What it is and why it is the shared fixture
- A.2 `CanonicalModelConfig` and `build_canonical_model`
- A.3 The model contract (`CanonicalModelContract`) and the diagnostic signals

**Appendix B — Glossary** — *[B_glossary.md](B_glossary.md)*
- DISV, DISU, Voronoi, idomain, pinch-out, ModelContext, spec, facade/engine, IES, geostruct, …

---

## Status

This manual is being authored chapter by chapter. Chapters not yet written are
listed above for structure; the **Getting Started** chapter is complete and is the
best entry point today.

| Chapter | Status |
|---|---|
| 3. Getting Started | ✅ Written |
| 4. Architecture & the Mental Model | ✅ Written |
| All others | ⬜ Outlined |
