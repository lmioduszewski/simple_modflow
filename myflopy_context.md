# myflopy — Project Context & Roadmap

Generated from a Claude Code web session (June 2026). Open this file in any session to restore context quickly.

---

## What myflopy is

A Python-first MODFLOW 6 toolkit built on top of FloPy, living on the `myflopy` branch of `lmioduszewski/simple_modflow`. It is **not** a thin wrapper — it is a full model-authoring, visualization, and calibration toolkit aimed at unstructured-grid groundwater models.

### Core architecture (src/myflopy/)

| File / Dir | What it does |
|---|---|
| `specs.py` (61 KB) | Composable `SimulationSpec`, `ModelSpec`, `PackageSpec`, `GridSpec` dataclasses — build a model by assembling specs |
| `package_api.py` (35 KB) | Top-level factory functions: `gwf()`, `gwt()`, `gwe()`, `prt()`, `lak()`, `sfr()`, `rch()`, `uzf()`, `mvr()`, `wel()`, etc. |
| `workspace.py` | `Project`, `ProjectLayout`, `ModelView`, `Run`, `load_run` — project/run management |
| `builders.py` | `build_gwf_gwt_exchange`, `build_gwf_gwe_exchange`, `build_gwf_prt_exchange`, `build_gwf_gwf_exchange`, `build_ims` |
| `advanced.py` | Higher-level spec factories: `lak_spec`, `sfr_spec`, `rch_spec`, etc. |
| `sources.py` | `RasterSource`, `ShapeSource`, `TableSource`, `GeoPackageSourceSpec`, `LiteralSource` |
| `geopackage.py` | `GeoPackageSource`, `CellSurfaceOffset` |
| `grid_spec_resolver.py` | Resolves `GridSpec` into FloPy DISV grid props |
| `modflow/mf6/grid/` | `VoronoiGridPlus`, `TriangleGrid`, `MeshBuildProfile` |
| `modflow/mf6/pest/` | Full PEST/pyemu calibration layer (see below) |
| `modflow/mf6/interactive_plotting.py` | HTML slider exports, `ModelVisualization`, `ParticleTrackingScene` |
| `modflow/mf6/parallel.py` | `ParallelModelWorkflow`, `ParallelEnvironment`, `ParallelSplitRun` |
| `modflow/mf6/canonical.py` | `CanonicalModelContract`, standardized output signals |
| `modflow/mf6/observations.py` | `HeadTargets`, `DrnFlowTargets`, `LakeStageTargets`, `SfrFlowTargets`, `SfrStageTargets` |

### Model types supported
- **GWF** — groundwater flow
- **GWT** — solute transport + GWF-GWT exchange
- **GWE** — energy transport + GWF-GWE exchange
- **PRT** — native MODFLOW 6 particle tracking + GWF-PRT exchange

---

## PEST / pyemu module — what's built

All code lives in `src/myflopy/modflow/mf6/pest/`.

### `PestProject` (project.py)
The main orchestrator. Wraps `pyemu.utils.PstFrom`.

```python
pest = PestProject(model, name="cal", workspace=Path("pest_ws"), start_datetime="2020-01-01")
pest.add_parameter(KPilotPointParameter(...))
pest.add_observation(HeadTargetObservationSpec(...))
pst = pest.build_pst("cal.pst")
```

Full pipeline: writes MF6 input → builds PstFrom workspace → prepares observation CSVs → prepares parameter template files → writes `pest_forward_config.json` → builds `.pst` → registers templates → finalizes observation values/weights.

Also: `pest.draw_prior(num_reals=1000)` for prior ensemble sampling.

### Parameters implemented

| Class | What it adjusts | Style |
|---|---|---|
| `KPilotPointParameter` | Hydraulic conductivity on Voronoi cells | Pilot points, IDW interpolation, zone-aware via polygon GIS source |
| `DrainElevationParameter` | DRN package cell elevations | By feature or grouped, absolute or additive |
| `DrainConductanceParameter` | DRN package conductance | By feature or grouped, multiplier or absolute |

Supports:
- Bounds modes: `absolute`, `multiplier`, `from_columns`, `multiplier_from_columns`
- Parameter space: `absolute` or `multiplier`
- PEST template file (`.tpl`) + CSV generation
- `ExpGeoStruct` / `build_geostruct` for pilot point kriging

### Observations implemented

| Class | Targets |
|---|---|
| `HeadTargetObservationSpec` | Water level at monitoring wells |
| `LakeStageObservationSpec` | LAK package stage |
| `SfrStageObservationSpec` | SFR package stream stage |
| `SfrFlowObservationSpec` | SFR package streamflow |
| `DrnFlowObservationSpec` | DRN zone seepage totals |

Each writes: simulated output CSV, `.ins` instruction file, target values CSV, location GeoPackage.

### Forward run (forward_run.py)
Injected script that PEST calls each iteration:
- Reads `pest_forward_config.json`
- IDW-interpolates pilot point K values to Voronoi cell centroids
- Applies drain elevation/conductance updates from parameter CSVs
- Reruns MF6
- Regenerates head and named-series (lake, SFR, DRN) output CSVs

### Results (results.py)
`PestRunResults`, `PestRunReview`, `open_pest_run` — post-run analysis and review.

### Demo files
- `gold_standard_demo.py` — full workflow with a real model
- `synthetic_demo.py` — synthetic test case

---

## PEST — what's missing / roadmap

| Item | Notes |
|---|---|
| `RechargeMultiplierParameter` | Recharge is a primary calibration target — not yet parameterized. Highest priority. |
| `KPilotPointParameter` from raster | Currently polygon-zone only. Need raster zone support (intersect cell centroids with a lithology raster). |
| `WelSpec` / well rate parameters | No adjustable pumping rates yet |
| UZF / GHB / CHD parameters | Not covered |
| Regularization helpers | No Tikhonov or preferred-value regularization setup |
| Ensemble methods | `draw_prior` exists; no IES/GLM ensemble smoother wrappers |
| Sensitivity / identifiability analysis | Not started |
| Parallel PEST++ workers | Not wired up |

---

## myflopy vs modflow-setup (DOI-USGS)

### What myflopy does that modflow-setup cannot

- Voronoi/unstructured DISV grids with programmatic refinement zones
- Solute transport (GWT), energy (GWE), particle tracking (PRT) models
- Native MODFLOW 6 PRT with interactive 3D visualization
- Parallel model workflows (`ParallelModelWorkflow`)
- Interactive HTML outputs — head map sliders, cross-sections, particle tracking
- PEST calibration on Voronoi grids (unique — modflow-setup PEST is structured-grid only)
- Declarative Python spec system
- Workspace/project/run management
- Canonical model contract / standardized output signals
- Surface water validation (`validate_surface_water_configuration`)

### What modflow-setup does that myflopy cannot (yet)

- Single YAML config file drives the entire model — no Python required
- Automatic CRS reprojection + resampling of all source data to the model grid
- Automated SFR from NHDPlus stream network data
- Local Grid Refinement (LGR) — parent + child structured model pairs
- GHB / RIV / CHD auto-built from GIS boundary shapefiles
- Reads existing MODFLOW array files as source data

### Priority borrowings from modflow-setup

| Feature | Why | Effort |
|---|---|---|
| `RechargeMultiplierParameter` | Completes PEST parameter coverage | Low |
| `KPilotPointParameter` from raster | Common lithology-zone workflow | Low–Medium |
| Idomain from boundary polygon | Removes boilerplate every model | Low |
| Boundary BCs from GIS features | Eliminates manual SPD construction | Medium |
| Spec YAML serialization (`SimulationSpec.to_yaml()`) | Reproducibility, shareability | Medium |
| Raster data ingestion pipeline | Core usability gap for layer data | Medium |
| SFR from stream centerline shapefile | Currently fully manual | High |
| LGR parent-child models | Niche but powerful | High |

---

## Why myflopy is worth using

**Real advantages over bare FloPy + nothing:**
1. Voronoi grid construction (`VoronoiGridPlus` + `TriangleGrid`) is significantly easier than raw FloPy triangle utilities
2. PEST calibration on Voronoi grids — no other public tool does this
3. Interactive HTML deliverables for clients (sliders, particle tracking, cross-sections)
4. Multi-physics (GWF+GWT+GWE+PRT) in one coherent spec system
5. Run/project management keeps multi-run workflows organized

**Honest caveats:**
- For structured-grid GWF models with lots of GIS data, modflow-setup is faster to set up
- No public docs or community — onboarding collaborators takes effort
- YAML serialization not yet implemented, so models aren't easily shareable as config files

---

## ACTIVE WORK: Layer / raster pipeline overhaul

> This is the next thing to implement. Decided in the June 2026 session.
> User wants this done **locally**, not pushed to the remote branch directly.

### Where the current pipeline lives

| File | Relevant code |
|---|---|
| `src/myflopy/modflow/mf6/grid/surfaces.py` | `get_raster_vals_at_centroids(vor, raster_files, labels)`, `get_gdf_topbtm_multilyr`, `get_raster_from_strike_dip` |
| `src/myflopy/modflow/mf6/grid/geometry.py` (≈312–378) | `reconcile_surfaces(vor, df, min_sep, trigger_sep, which)`, `adjust_cells_by_id`, `adjust_top_btm_overlaps` |
| `src/myflopy/modflow/mf6/grid/voronoi.py` | `gdf_topbtm` property, `nlay` property, `reconcile_surfaces` method, `idomain` setter |
| `src/myflopy/package_api.py` (≈437–465) | `disv(*, nlay, ncpl, ..., top, botm, idomain=None, ...)` |

### Gaps vs. modflow-setup's `discretization.py`

modflow-setup has: `make_idomain`, `create_vertical_pass_through_cells`,
`fix_model_layer_conflicts`, `fill_cells_vertically`, `fill_empty_layers`,
`verify_minimum_layer_thickness`. myflopy's gaps:

1. **No `IDOMAIN = -1` pass-through.** `reconcile_surfaces` only forces `min_sep`
   spacing — it never produces true pinch-outs. A thin/zero-thickness middle layer
   stays in the solution as a paper-thin cell instead of being removed with vertical
   flow passing through.
2. **No automatic idomain** from thickness or NaN coverage.
3. **Raster sampler does not reproject** — centroids are sampled in the grid CRS with
   no transform to the raster CRS. Silent garbage if CRS differ.
4. **nodata not converted to NaN.** Bug in `get_raster_vals_at_centroids`:
   ```python
   sampled = np.array(
       [val[0] if val is not None else np.nan for val in src.sample(zip(xs, ys))]
   )
   ```
   `rasterio.sample()` returns the nodata value (e.g. `-9999`), never `None`, so the
   guard never fires and nodata leaks in as a real elevation.
5. **No convenience layer constructors** — flat layers, constant-thickness, offset-below.
6. **No minimum-thickness validation / reporting.**

### Key design fact (confirmed)

For Voronoi/DISV grids the **grid IS the domain**, so "idomain from boundary polygon"
is irrelevant. Only two things matter: interior holes (carved by polygons) and
**layer pinch-outs**. MODFLOW 6 `IDOMAIN = -1` = vertical pass-through: the cell is
removed from the solution but vertical flow passes through, connecting the layer above
to the layer below. That is exactly the mechanism for "make a middle layer inactive in
an area and have the other layers act like it's not there."

### THE FIRST SLICE (what to implement)

New file: `src/myflopy/modflow/mf6/grid/layer_stack.py`

Source classes (each resolves to a per-cell elevation array):
- `Raster(path, fill=None)` — sample a GeoTIFF at centroids. Must reproject centroids
  to the raster CRS and convert nodata → NaN. `fill="propagate"` carries the previous
  surface down where NaN.
- `Flat(elevation)` — constant elevation everywhere.
- `OffsetBelow(distance)` — previous surface minus a fixed distance.
- `ConstantThickness(thickness)` — previous surface minus thickness (alias-ish of
  OffsetBelow but named for layer intent).
- (stretch) `StrikeDip(...)` — wrap existing `get_raster_from_strike_dip`.

`LayerStack` class:
```python
from myflopy.layers import LayerStack, Raster, Flat, OffsetBelow, ConstantThickness

stack = LayerStack(vor=vor)
stack.top(Raster("ground.tif"))
stack.add(Raster("l1_bot.tif"), name="alluvium")
stack.add(Flat(550.0),          name="clay")
stack.add(ConstantThickness(30),name="sand")
stack.add(OffsetBelow(20),      name="weathered")
stack.add(Raster("bedrock.tif"),name="bedrock", fill="propagate")

result = stack.build(
    minimum_thickness=1.0,
    pinch_out=True,      # thin cells -> idomain = -1 (pass-through). OPT-IN.
    reconcile="bottom",  # reuse existing reconcile_surfaces
)
# result.top, result.botm, result.idomain  -> feed mf.disv(...)
```

`build()` order of operations:
1. Resolve each source to a per-cell elevation array (top + each bottom).
2. `reconcile` (reuse `geometry.reconcile_surfaces`) to enforce ordering.
3. Thickness check against `minimum_thickness`.
4. **Pinch-out (opt-in):** where thickness < `minimum_thickness`, set
   `idomain = -1` (NOT `0` — `0` would block vertical flow). Default `pinch_out=False`
   keeps existing `min_sep` behavior so nothing breaks.
5. Return a small result object exposing `.top`, `.botm`, `.idomain`.

Also fix `get_raster_vals_at_centroids` in `surfaces.py`:
- Reproject centroid coords from `vor.crs` to the raster CRS before sampling.
- Read `src.nodata` and convert matches → NaN (don't rely on `is not None`).

Export from `src/myflopy/__init__.py` (or `grid/__init__.py`) so
`from myflopy.layers import ...` resolves.

Pinch-out is **opt-in**; the default path must reproduce today's behavior.
