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

## Two API layers + capability map
This codebase has **two API generations** — check both before adding code, or you risk
duplicating an existing capability:
- **Modern declarative spec API** (preferred): `src/myflopy/*.py` — `package_api`, `specs`,
  `sources`, `geopackage`, `surfaces`, `advanced`, `builders`, `workspace`. Authoritative
  export list: `src/myflopy/__init__.py`.
- **Legacy OO API**: `src/myflopy/modflow/mf6/*.py` — `simplemodel`, `boundaries`, `sfr`,
  `lakes`, `recharge`, and builder classes (`SFRBuilder`, `LAKBuilder`, `RCHBuilder`, …).

**`myflopy_context.md` is the accurate, code-derived capability map** (rebuilt 2026-06-20).
Treat any "gap" as a hypothesis to re-verify against the code before building.

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
> Re-verified against the code on 2026-06-20. **Most items previously listed here are
> already built** (see `myflopy_context.md`). Genuine remaining gaps only:
- **YAML/TOML spec serialization** — thin wrapper over the existing
  `SimulationSpec.to_dict()` / `from_dict()` (round-trip already implemented in `specs.py`)
- **NHDPlus direct SFR reader** — `SFRBuilder` already builds reaches from any stream
  centerline LineString table; only national NHDPlus ingestion is missing
- **Area-weighted raster resampling** — sampling is point-at-centroid today
- **Reading existing MODFLOW array files** as source data
- **LGR parent-child model pairs** — absent (niche for a Voronoi-first toolkit)

Already built — do NOT rebuild: GIS-driven BCs (`GeoPackageSource.chd/ghb/drn/wel/rch`,
`mf.ghb.gpkg`, legacy `Boundaries`), CRS reprojection (vector + raster), layer surfaces +
reconcile + **pinch-out/idomain** (`surfaces.py`), SFR from centerline (`SFRBuilder`),
recharge from GIS/PRISM (`RCHBuilder`).

## Comparison: myflopy vs modflow-setup
- myflopy wins on: unstructured grids, transport/energy/PRT models, visualization, parallel
  runs, Python-first API, GIS-driven BCs, SFR-from-centerline, PEST on Voronoi grids
- modflow-setup wins on: single-file YAML no-code setup, NHDPlus SFR, LGR
