# myflopy — Capability Map (code-derived)

> **Read this first.** This file was rebuilt on **2026-06-20** by scanning the actual
> `src/myflopy/` tree, not from memory or a prior plan. The previous version of this
> file listed features as "missing" that were already implemented, which caused real
> duplicated work (a standalone `LayerStack` that duplicated `LayerSurfaces`).
> The **PEST**, **serialization**, and **area-weighted sampling** rows were re-verified
> against the code and refreshed on **2026-07-03** (post the 2026-06-22 PEST unification):
> the legacy `KPilotPointParameter`/`build_pst` API is gone, IES + workers ship, and
> exchanges + hooks now serialize.
>
> **Discipline to avoid that again:** before building anything, `grep` the codebase for
> the capability **in both API layers** (see below). Treat every "gap" in this file as a
> hypothesis to re-verify against the code, not a fact. When in doubt, read the module.

---

## Two API generations (the core duplication trap)

A capability often exists in one or both of these. Always check both before adding code.

| Layer | Where | Style |
|---|---|---|
| **Modern declarative spec API** (preferred) | `src/myflopy/*.py` — `package_api.py`, `specs.py`, `sources.py`, `geopackage.py`, `surfaces.py`, `advanced.py`, `builders.py`, `workspace.py`, `grid_spec_resolver.py` | `mf.gwf(...)`, `mf.disv(...)`, `mf.ghb.gpkg(...)`, dataclass specs, `to_dict`/`from_dict` |
| **Legacy OO API** | `src/myflopy/modflow/mf6/*.py` — `simplemodel.py`, `simulation/`, `boundaries.py`, `drn.py`/`ghb.py`/`chd.py`, `recharge.py`, `lakes.py`, `sfr.py`, `mvr.py`, `uzf.py`, `ModelSurface` | `SimulationBase`, `Boundaries`, builder classes (`SFRBuilder`, `LAKBuilder`, `RCHBuilder`, …) |

The authoritative public surface is `src/myflopy/__init__.py` (`_EXPORTS`).

---

## Capability map

Legend: ✅ built · 🟡 partial / has primitives · ❌ missing

### Grids & discretization
| Capability | Status | Where |
|---|---|---|
| Voronoi / unstructured grids (DISV, DISU) | ✅ | `grid/voronoi.py` (`VoronoiGridPlus`), `grid/triangle.py` (`TriangleGrid`), `simulation/discretization.py` |
| Programmatic refinement regions | ✅ | `TriangleGrid`, `MeshBuildProfile`, `grid/seed_optimization.py`, `grid/mesh_quality.py` |
| Layer surfaces from raster / flat / points / contours(GRASS) | ✅ | `surfaces.py` (`Surface`, `LayerSurfaces`) |
| Relative surfaces (`offset_below`, `constant_thickness`) | ✅ | `surfaces.py` (added 2026-06) |
| Surface reconcile (enforce top-down ordering) | ✅ | `grid/geometry.py` `reconcile_surfaces` |
| Pinch-out → `idomain = -1` (vertical pass-through) | ✅ | `surfaces.py` `to_disv(pinch_out=...)`, `top_botm_idomain`, `thickness_report` (added 2026-06) |
| Raster sampling at centroids (reproject + nodata→NaN) | ✅ | `grid/surfaces.py` `get_raster_vals_at_centroids` (fixed 2026-06) |
| idomain from boundary polygon | N/A | For Voronoi the grid *is* the domain; only interior holes + layer pinch-outs matter |
| Area-weighted raster **resampling** (vs point-at-centroid) | ✅ | `grid/surfaces.py` `_area_weighted_sample` — per-cell mean of pixels within each Voronoi polygon; **default** (`method="area"`), centroid fallback for uncovered cells |
| `GridSpec.resolve` wiring | 🟡 | `grid_spec_resolver.py` resolves `voronoi` + `python`; `GridSpec.structured`/`from_geopackage` **fail fast** (unwired) — use `voronoi`/`python`/`from_object` |
| LGR parent/child structured grids | ❌ | no `Lgr`/`ModflowLgr` anywhere (niche for an unstructured-first tool) |

### Sources & GIS ingestion
| Capability | Status | Where |
|---|---|---|
| Declarative sources (raster, shape, table, geopackage, literal) | ✅ | `sources.py` (`RasterSource`, `ShapeSource`, `TableSource`, `GeoPackageSourceSpec`, `LiteralSource`) |
| GeoPackage feature → cell pipeline | ✅ | `geopackage.py` `GeoPackageSource` (`.chd/.ghb/.drn/.wel/.rch`) |
| CRS reprojection of vector sources | ✅ | `boundaries.py` `Boundaries.gdf` (`to_crs`) |
| CRS reprojection of rasters | ✅ | `grid/surfaces.py` (fixed 2026-06) |
| Import existing MODFLOW layer surfaces | 🟡 | `LayerStack.from_modflow` / `LayerSurfaces.from_modflow` read top/botm/idomain from a built MF6 model; a raw array-file *source* reader is still missing |

### Boundary conditions
| Capability | Status | Where |
|---|---|---|
| CHD/GHB/DRN/WEL/RCH/UZF package builders | ✅ | `package_api.py` (`_CHDPackage`…`_MVRPackage`), `advanced.py` `*_spec` |
| **BCs auto-built from GIS** (polygon/line) | ✅ | `mf.ghb.gpkg(...)`, `mf.drn.gpkg(...)`; `GeoPackageSource`; legacy `Boundaries` |
| Perimeter / edge-cell BCs (CHD/GHB on grid edge) | ✅ | `Boundaries.edge_intersections`, `edges_only=True` on `.gpkg()` builders |
| Cells ordered along a line (for line BCs / SFR) | ✅ | `Boundaries.sorted_cells_along_line` |
| Recharge from GIS + PRISM precip scaling + area scaling | ✅ | `recharge.py` `RCHBuilder`, `PrismPrecipScaling`, `Boundaries.shp_to_vor_poly_scale` |
| Dedicated RIV builder | ❌ | no `_RIV`/`riv`; represent via GHB/DRN |

### Advanced packages
| Capability | Status | Where |
|---|---|---|
| LAK (lakes, connections, outlets, lake tables) | ✅ | `lakes.py` (`LAKBuilder`, `LakeTableBuilder`), `advanced.py` `lak_spec` |
| **SFR from stream centerline** (LineString → reaches by cell intersection) | ✅ | `sfr.py` `SFRBuilder`/`StreamNetwork`; `rlen` from geometry; reach-top from grid surfaces |
| SFR connectivity inference | ✅ | `connection_mode` = `automatic` (geometric) / `nodes` / `explicit` |
| SFR direct from **NHDPlus** national dataset (sfrmaker-style) | ❌ | only custom centerline tables; no NHDPlus reader |
| MVR (mover), UZF | ✅ | `mvr.py` (`MVRBuilder`), `uzf.py` (`UZFBuilder`) |
| Surface-water configuration validation | ✅ | `surface_water_validation.py` `validate_surface_water_configuration` |

### Multi-physics, runs, viz
| Capability | Status | Where |
|---|---|---|
| GWF / GWT / GWE / PRT models + exchanges | ✅ | `package_api.py` `gwf/gwt/gwe/prt`; `builders.py` `build_gwf_{gwt,gwe,prt,gwf}_exchange` |
| MODPATH-style particle tracking (mp3du) | ✅ | `modflow/mp3du/` (`ParticleTrackingInput`, `run_particle_tracking`) |
| Workspace / project / run management | ✅ | `workspace.py` (`Project`, `Run`, `load_run`), `project/` |
| Parallel model split (partition, MPI) | ✅ | `parallel.py` (`ParallelModelWorkflow`, `ParallelSplitRun`) |
| Interactive HTML viz (sliders, cross-sections, particle scenes) | ✅ | `interactive_plotting.py` |
| Canonical model contract (standardized signals) | ✅ | `canonical.py`, `canonical_example.py` |

### Serialization
| Capability | Status | Where |
|---|---|---|
| `to_dict`/`from_dict` round-trip for all specs | ✅ | `specs.py` — Package/Grid/Model/SimulationSpec **incl. inter-model exchanges + post-build hooks** (importable builders; lambdas/closures fail loud) |
| Workspace pickling / native MF6 reload | ✅ | `workspace.py`, `project/` (`load_mf6_run`) |
| **YAML/TOML** config file wrapper | ❌ | no `to_yaml`; would be a thin wrapper over existing `to_dict`/`from_dict` |

### PEST / pyemu (`modflow/mf6/pest/`)
| Capability | Status | Where |
|---|---|---|
| `PestProject` orchestration (native PstFrom build + forward run) | ✅ | `pest/project.py`, `pest/forward_run.py` — construct via `model.pest(name, ...)` |
| **Unified** `cal.parameterize(target, style=...)` | ✅ | `pest/project.py` + `pest/native_parameters.py`; targets `k/k33/recharge/chd/ghb/drn/wel`, styles `constant/zone/grid/pilotpoints` |
| Pilot points on Voronoi (IDW; pyEMU's own are unusable on DISU) | ✅ | `pest/pilot_points.py`; `pp_space` net or explicit `pp_points` |
| Observations: head, lake stage, SFR stage/flow, DRN flow | ✅ | `pest/observations.py` (`cal.observe` / `cal.forecast`) |
| PESTPP-IES + prior Monte Carlo + parallel workers | ✅ | `pest/project.py` `run_ies(workers=)`, `prior`, `draw_prior` |
| Results / review / reopen a run | ✅ | `pest/ies.py` (`open_ies_run`, `IesResults`), `pest/runs.py` (`find_pest_runs`, `model.pest_runs`) |
| Geostats for grid / pilot-point priors | ✅ | `pest/geostats.py` (`ExpGeoStruct`, `build_geostruct`) |
| Pilot points / zones from **raster** (vs pp-net / polygon) | 🟡 | pp-net + polygon zones exist; raster-driven zone arrays do not |
| UZF parameters; Tikhonov/preferred-value **regularization**; identifiability | ❌ | not yet a `parameterize` target / not built |

---

## The genuinely short gap list (verify each before building)

1. **YAML/TOML serialization** — thin `yaml.dump(spec.to_dict())` / `from_dict(yaml.load(...))` wrapper. `to_dict`/`from_dict` are now complete (incl. exchanges + hooks); only the file-format wrapper is missing. *(Low effort.)*
2. **PEST raster-driven zones/pilot points**, Tikhonov/preferred-value **regularization**, **identifiability**, **UZF parameters** — the real PEST roadmap. (IES + prior MC + PEST++ workers already ship via `run_ies(workers=)`.)
3. **NHDPlus SFR reader** — only piece of "SFR from streams" not already covered by `SFRBuilder`.
4. **`GridSpec.structured`/`from_geopackage` resolution** and a **raw MF6 array-file source reader** — both currently fail-fast / missing; niche for a Voronoi-first tool.
5. **LGR** — absent; likely not worth it for a Voronoi-first toolkit.

## Already done — do NOT rebuild
GIS-driven BCs (CHD/GHB/DRN/WEL/RCH via `GeoPackageSource` + `mf.*.gpkg`), edge/perimeter
BCs, recharge from GIS/PRISM, CRS reprojection (vector + raster), **area-weighted raster
sampling** (default), layer surfaces + reconcile + **pinch-out/idomain**, **SFR from
centerline** (`SFRBuilder`), LAK/MVR/UZF builders, full spec `to_dict`/`from_dict`
(**incl. exchanges + hooks**), multi-physics + exchanges, parallel split, viz, and the
**unified PEST** stack — `cal.parameterize` (constant/zone/grid/pilotpoints), pilot points
on Voronoi, PESTPP-IES + prior MC + parallel workers, and run review (`open_ies_run`).
(The old version of this file wrongly listed several of these as missing.)
