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

## Canonical model-building API: package-first (use this for new work)
**The package-first API in `package_api.py` is THE preferred, canonical way to build models.**
Prefer it for all new model-building work; do not reach for the legacy OO path below unless
you have a specific reason. The OO builder classes (`SFRBuilder`, `LAKBuilder`, `UZFBuilder`,
`RCHBuilder`, `MVRBuilder`, `GHB`, `Recharge`, …) are the **engine underneath** the
package-first facade, not a competing API — e.g. `mf.uzf(...)` literally calls
`UZFBuilder(...).build()`. (The canonical model in `canonical_example.py` is written in the
older imperative builder style for historical/computed-cell reasons; that does NOT make the
imperative path "more proven" — package-first is the same engine and is the one to teach.)

### How it fits together (the mental model)
```
Project            ← durable workspace + run/scenario lifecycle (workspace.py); holds NO geometry
  └─ SimulationSpec ← one MF6 simulation: tdis + ims solver + the model(s)
       └─ ModelSpec = mf.gwf(name, context=…, packages=[…])   ← one model
            ├─ ModelContext(grid=, domain=, surfaces=, dates=)  ← geometry; rides on the MODEL, not the project
            └─ packages: mf.disv, mf.npf/ic/sto/oc,
                         mf.chd/ghb/drn/wel/rch (+ .gpkg / .flopy),
                         mf.uzf/sfr/lak (+ .flopy), mf.mvr
```
- **`Project(root, name=)`** (`workspace.py`): `add_grid/add_package/add_simulation` (reusable
  libraries), `prepare_run(name, sim)` → `Run` (built in memory; inspect `run.model(...)`),
  `run.execute()` writes + runs MF6. The Project is the lifecycle wrapper; geometry lives on the model.
- **`ModelContext`** (`specs.py`) attaches to the **model** via `mf.gwf(..., context=ctx)`, NOT to
  the project. It carries `grid`, `domain` (idomain), `surfaces`, `dates`. The built model exposes
  it as `model.myflopy_context`. The GIS-aware package helpers (`mf.uzf`, `mf.sfr`, `mf.X.gpkg`)
  take `context=` so they can map features/cells onto the grid.
- **Package-first surface**: `mf.gwf/gwt/gwe/prt`; `mf.disv`; `mf.ic/npf/sto/oc/tdis/ims`;
  list BCs `mf.chd/ghb/drn/wel/rch` each with `()` (direct data), `.gpkg(path, context=, nper=)`
  (from GeoPackage), `.flopy(...)` (raw FloPy escape hatch); advanced `mf.uzf/sfr/lak` (`()` =
  high-level builder, `.flopy(...)` = raw); `mf.mvr` with `mf.Move(mf.MoverConnection("sfr",0),
  mf.MoverConnection("lak",0))`. **MVR is validated**: moved packages must be declared in the
  model AND ordered before the mover (`test_advanced_specs`).
- **Layers** (same facade/engine pattern): use the facade **`mf.LayerStack`** (`layers.py`) —
  `LayerStack(vor, top=Raster("ground.tif")).add("sand", thickness=20, pinch="inactive")
  .add("clay", bottom=Contours(...)).build()` → disv-ready top/botm/idomain (plus `.qc()`,
  `.cross_section()/.surface_3d()/.vtk_3d()`, `from_modflow`). It **compiles to** the
  `LayerSurfaces` engine (`surfaces.py`: area-weighted sampling, top-down reconcile, pinch-out
  → idomain), which uses atomic `Surface` objects (`raster`/`from_contours`/`from_points`/
  `from_array`/constant + algebra). Use `LayerStack`; `LayerSurfaces`/`Surface` are the engine/atoms.

### Grid: eager (works now) vs deferred (a known seam)
- **Eager** (use for the full GIS stack): build the grid object first (`VoronoiGridPlus`, or
  `mf.GridSpec.voronoi(...).resolve(workspace)`), put it in `ModelContext(grid=vor, domain=idomain)`,
  then declare packages. **Required today** for `mf.uzf/sfr/lak/.gpkg` because they resolve cells
  eagerly at declaration (e.g. `UZFBuilder.build()` bakes resolved data into the `PackageSpec`).
- **Deferred** (`mf.gwf().with_grid(mf.GridSpec.voronoi(boundary=<gpkg>, refinement=<gpkg>))`):
  the project builds the grid at run time into `run.workspace/_grid/<model>` and populates
  `context.grid`. In `ModelSpec.build` the grid resolves BEFORE packages build, and the model
  carries `model.myflopy_context` — so deferred GIS packages are *feasible*, just not wired: the
  helpers would need a grid-lazy mode emitting a `PackageSpec` whose `build(model)` resolves against
  `model.myflopy_context` instead of resolving eagerly. Until then, deferred GridSpec composes only
  with disv + simple BCs, not with `mf.uzf/sfr/lak`.

Reference: `examples/mf6/package_first_full_stack.py` (full package-first stack on a Project).

## Legacy OO API (the engine; avoid for new model assembly)
`src/myflopy/modflow/mf6/*.py` — `simplemodel`, `boundaries`, `sfr`, `lakes`, `recharge`, and the
builder classes. Still the engine under the facade and still used by `canonical_example.py`. When
adding a capability, **grep both layers first** (package-first + these) to avoid duplicating one.
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
