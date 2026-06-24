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

**`docs/myflopy_context.md` is the accurate, code-derived capability map** (rebuilt 2026-06-20).
Treat any "gap" as a hypothesis to re-verify against the code before building.

## Planned work: PEST / pyemu integration

### What's already built (`src/myflopy/modflow/mf6/pest/`)
- **One unified parameterization API**: `cal.parameterize(target, style=...)` compiling
  to native `pyemu.utils.PstFrom`. Targets: `k`, `k33`, `recharge`, `chd`, `ghb.cond`/
  `ghb.bhead`, `drn.cond`/`drn.elev`, `wel`. Styles: `constant`, `zone`, `grid` (one
  geostat-correlated multiplier per Voronoi cell + correlated prior), and `pilotpoints`
  (IDW from a `pp_space` net or explicit `pp_points` — `pilot_points.py`; pyEMU's own
  pilot points are unusable on unstructured grids). The legacy `add_parameter` +
  `build_pst` + `KPilotPointParameter`/`DrainElevation`/`DrainConductance` specs were
  RETIRED (deleted 2026-06-22) — do not reintroduce them.
- `PestProject` (`project.py`) — orchestrates pyEMU/PstFrom; native `build()` + `run_ies`/
  `prior` (`workers=` runs parallel PESTPP-IES agents). **Construct it with
  `model.pest(name, start_datetime=...)`** (the front door — defaults the workspace to
  `<model workspace>/pest/<name>`), not by importing `PestProject` directly.
- **Observations** accepted directly by `cal.observe(...)`/`cal.forecast(...)`: `HeadTargets`,
  `LakeStageTargets`, `SfrStageTargets`, `SfrFlowTargets`, `DrnFlowTargets` (or their
  pre-built `*ObservationSpec`). All are wired into the native build + forward run.
- `forward_run.py` — injected forward run: pyEMU's `apply_list_and_array_pars` for params,
  plus post-processors that regenerate head + named-series (lake/SFR/DRN) simulated CSVs.
- **Run discovery + review**: `model.pest_runs` / `run.pest_runs` (`runs.py`,
  `find_pest_runs`/`PestRunHandle`) list calibrations done on a model; `.review()` reopens
  one as `IesResults` (`ies.py`, `open_ies_run`) — phi, ensemble-vs-obs, forecasts,
  parameter-field maps. (The old `results.py` deterministic-review layer was deleted.)
- `geostats.py` — `ExpGeoStruct` / `build_geostruct`: the **single** geostruct builder
  (`project._geostruct_for` delegates to it); `grid` style takes `anisotropy`/`bearing`/
  `nugget` via `parameterize`.
- The PEST notebooks (`canonical_04/05/06`) calibrate the **canonical valley
  model itself** (perturb its K/recharge as "truth", then calibrate back); 06 uses
  pilot-point K. The old standalone demo models (`gold_standard_demo`, `synthetic_demo`,
  `modern_pest_demo.build_calibration_demo`) were retired in favor of one model everywhere.

### What's missing / next steps for PEST
> Verify against the code before building — most of the old list is now done
> (recharge/wel/ghb/chd are `parameterize` targets; IES + parallel workers ship via
> `run_ies(workers=)`; named-series obs are wired).
1. **Pilot points / zones from raster** — pilot nets + polygon zones exist; raster-driven
   zone arrays do not.
2. **UZF parameters** — not yet a `parameterize` target (GHB/CHD/WEL/DRN already are).
3. **Regularization helpers** — no Tikhonov or preferred-value regularization setup.
4. **Sensitivity / identifiability analysis** helpers.

### Broader things myflopy could learn from modflow-setup (DOI-USGS)
> Re-verified against the code on 2026-06-20. **Most items previously listed here are
> already built** (see `docs/myflopy_context.md`). Genuine remaining gaps only:
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
