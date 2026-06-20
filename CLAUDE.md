# myflopy / simple_modflow — Project Context for Claude

## Repository
- GitHub: `lmioduszewski/simple_modflow`
- Active development branch: `myflopy`
- Package lives at `src/myflopy/` on that branch

## What this project is
A Python-first MODFLOW 6 toolkit built on top of FloPy. Key strengths:
- **Voronoi/unstructured grids** (DISV/DISU) via `VoronoiGridPlus` + `TriangleGrid`
- **Declarative spec API** — `SimulationSpec`, `ModelSpec`, `PackageSpec`, `GridSpec` dataclasses in `specs.py`
- **Package-first API** — `gwf()`, `gwt()`, `gwe()`, `prt()`, `lak()`, `sfr()`, etc. in `package_api.py`
- **Multi-physics** — GWF, GWT (transport), GWE (energy), PRT (particle tracking) with exchange packages
- **Rich interactive visualization** — Plotly choropleth maps, HTML sliders, cross-sections, animations
- **Workspace/project management** — `Project`, `Run`, `load_run` in `workspace.py`
- **Parallel model workflows** — `ParallelModelWorkflow` in `parallel.py`

## Planned work: PEST / pyemu integration

### What's already built (`src/myflopy/modflow/mf6/pest/`)
- `PestProject` (`project.py`) — orchestrates pyEMU/PstFrom; builds `.pst`, manages forward run
- `KPilotPointParameter` — pilot point K on Voronoi cells, IDW interpolation, zone-aware (polygon source)
- `DrainElevationParameter` — adjustable drain elevations (absolute/additive, by feature or group)
- `DrainConductanceParameter` — adjustable drain conductance (multiplier/absolute, by feature or group)
- `HeadTargetObservationSpec` — water level observations at monitoring wells
- `LakeStageObservationSpec` — LAK package stage observations
- `SfrStageObservationSpec` / `SfrFlowObservationSpec` — SFR stage and flow observations
- `DrnFlowObservationSpec` — drain zone seepage observations
- `forward_run.py` — injected forward run: applies K/drain params, regenerates output CSVs
- `results.py` — `PestRunResults`, `PestRunReview`, `open_pest_run` for post-run analysis
- `geostats.py` — `ExpGeoStruct` / `build_geostruct` for pilot point kriging
- Demo files: `gold_standard_demo.py`, `synthetic_demo.py`

### What's missing / next steps for PEST
1. **`RechargeMultiplierParameter`** — recharge is a primary calibration target, not yet parameterized
2. **`KPilotPointParameter` from raster** — currently polygon-zone only; need raster zone support
3. **Well package parameters** (`WelSpec`) — no adjustable pumping rates
4. **UZF / GHB / CHD parameters** — not yet covered
5. **Regularization helpers** — no Tikhonov or preferred-value regularization setup
6. **Ensemble methods** — `draw_prior` exists but no IES/GLM ensemble smoother wrappers
7. **Sensitivity / identifiability analysis** helpers
8. **Parallel PEST++ workers** — not wired up

### Broader things myflopy could learn from modflow-setup (DOI-USGS)
- Automated raster/shapefile ingestion with CRS reprojection + resampling strategies
- Systematic layer data setup (tops/bottoms/K from rasters with pinch-out reconciliation)
- Boundary conditions (GHB, RIV, CHD) auto-built from GIS boundary features
- SFR network from stream centerline shapefile (snap reaches, compute lengths/slopes)
- Spec YAML/TOML serialization for reproducibility (`SimulationSpec.to_yaml()`)
- `idomain` auto-derived from model boundary polygon
- LGR parent-child model pairs with GWF-GWF exchange

## Comparison: myflopy vs modflow-setup
- myflopy wins on: unstructured grids, transport/energy/PRT models, visualization, parallel runs, Python-first API
- modflow-setup wins on: YAML no-code setup, NHDPlus SFR, LGR, automatic GIS data ingestion
