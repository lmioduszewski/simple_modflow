# myflopy — Capability Map (code-derived)

> **Read this first.** This file was rebuilt on **2026-06-20** by scanning the actual
> `src/myflopy/` tree, not from memory or a prior plan. The previous version of this
> file listed features as "missing" that were already implemented, which caused real
> duplicated work (a standalone `LayerStack` that duplicated `LayerSurfaces`).
> The **PEST**, **serialization**, and **area-weighted sampling** rows were re-verified
> against the code and refreshed on **2026-07-03** (post the 2026-06-22 PEST unification):
> the legacy `KPilotPointParameter`/`build_pst` API is gone, IES + workers ship, and
> exchanges + hooks now serialize.
> **Refreshed 2026-07-07:** added the inspection/diff stack (`diff()`, `model.config`,
> ModelDiff tiers), the viz front door + unified view grammar, the sectioned hover
> system, the colorscale policy, map-mosaic view sync, ATS, and `to_xugrid`. The
> forward roadmap now lives in `docs/implementation_plan_2026-07.md` (rev. 3).
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
| GeoPackage feature → cell pipeline | ✅ | `geopackage.py` `GeoPackageSource` (`.chd/.ghb/.drn/.riv/.wel/.rch`) |
| CRS reprojection of vector sources | ✅ | `boundaries.py` `Boundaries.gdf` (`to_crs`) |
| CRS reprojection of rasters | ✅ | `grid/surfaces.py` (fixed 2026-06) |
| Import existing MODFLOW layer surfaces | 🟡 | `LayerStack.from_modflow` / `LayerSurfaces.from_modflow` read top/botm/idomain from a built MF6 model; a raw array-file *source* reader is still missing |

### Boundary conditions
| Capability | Status | Where |
|---|---|---|
| CHD/GHB/DRN/RIV/WEL/RCH/UZF package builders | ✅ | `package_api.py` (`_CHDPackage`…`_MVRPackage`), `advanced.py` `*_spec` |
| **BCs auto-built from GIS** (polygon/line) | ✅ | `mf.ghb.gpkg(...)`, `mf.drn.gpkg(...)`; `GeoPackageSource`; legacy `Boundaries` |
| Perimeter / edge-cell BCs (CHD/GHB on grid edge) | ✅ | `Boundaries.edge_intersections`, `edges_only=True` on `.gpkg()` builders |
| Cells ordered along a line (for line BCs / SFR) | ✅ | `Boundaries.sorted_cells_along_line` |
| Recharge from GIS + PRISM precip scaling + area scaling | ✅ | `recharge.py` `RCHBuilder`, `PrismPrecipScaling`, `Boundaries.shp_to_vor_poly_scale` |
| Dedicated RIV builder | ✅ | `mf.riv` (`_RIVPackage`: `()`/`.gpkg`/`.flopy`), `riv_spec`, `GeoPackageSource.riv`, registry entry (`stage/cond/rbot` earth, q RdBu) — added 2026-07-17 |

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
| GWT/GWE **package helpers** (adv/dsp/mst/ssm/cnc…; est/cnd/ctp…) | ❌ | raw `PackageSpec(flopy.mf6.ModflowGwt*, ...)` only; `mf.ic/oc/disv` hardcode GWF classes — plan §5.3 |
| GWT/GWE **results tier** (conc/temp reader, maps, group/diff, obs) | ❌ | no concentration/temperature reader exists anywhere — plan §6.0–6.2 |
| MF6 PRT (release points, run, pathlines, 3-D scenes) | ✅ | `prt.py` (`PRTProject`, `PRTRunResults`, `open_prt_run`), `model.particle_tracking` |
| PRT derived cell maps (travel-time / endpoints choropleths) | ❌ | pathlines/scene/map exist; grammar-integrated cell maps are plan §6.3 |
| MODPATH-style particle tracking (mp3du) | ✅ | `modflow/mp3du/` (`ParticleTrackingInput`, `run_particle_tracking`) |
| Workspace / project / run management | ✅ | `workspace.py` (`Project`, `Run`, `load_run`), `project/` |
| Parallel model split (partition, MPI) | ✅ | `parallel.py` (`ParallelModelWorkflow`, `ParallelSplitRun`) |
| Interactive HTML viz (sliders, cross-sections, particle scenes) | ✅ | `interactive_plotting.py` |
| Adaptive time stepping (ATS; zero-based `ats=` periods) | ✅ | `simulation/discretization.py` (`TemporalDiscretization(ats=...)`) |
| xarray/ugrid interchange (`to_xugrid`, NetCDF-ready) | ✅ | `grid/voronoi.py` `VoronoiGridPlus.to_xugrid`, `headsplus.py` `HeadsPlus.to_xugrid` (extra: `xugrid`) |
| Canonical model contract (standardized signals) | ✅ | `canonical.py`, `canonical_example.py`; profiles: full 100×100, `validation()` 50×50, `testing()` 21×21 (smallest contract-complete; the test suite's shared session fixture + `canonical_fast_tour.ipynb`) |

### Visualization front door, grammar, hover, colors
| Capability | Status | Where |
|---|---|---|
| Plotting front door (`viz.Fig/subplots/mpl_axes/PALETTE`, figs backend) | ✅ | `viz.py` — import every figure from here, not raw plotly/matplotlib. figs is vendored (`myflopy/_vendor/figs`, synced via `scripts/sync_vendored_figs.py`): viz.py imports the live figs first and falls back to the snapshot, so pip installs work without the local figs project |
| Unified view grammar: `map/plot/xs` + `mosaic/animate` on every leaf (model / group / diff), `backend="plotly"\|"mpl"` | ✅ | `package_plotting.py` (`SpatialView`, composers), `viz.mosaic` |
| Map mosaics framed to data + **synced pan/zoom** (`sync_views=`) | ✅ | `viz.py` (`shared_map_view`, `_map_sync_post_script`), `Choro.map_view` |
| **Sectioned hover** (`HoverSpec`/`HoverStyle`; layer/surface tables, dry marking; `hover_*` sugar on every map verb) | ✅ | `utils/datatypes/hover.py`; defaults wired on all choropleth paths; LAK/SFR maps join feature stage |
| **Colorscale policy** (diverging only for signed q-like + diff maps; `'earth'` otherwise) | ✅ | `package_registry.py` + call sites; pinned by `tests/test_colorscale_policy.py` |
| Choropleth engine (plotly + mpl backends, contours, hillshade) | ✅ | `utils/datatypes/choros.py` (`Choro`), `simulation/accessors.py` `build_choro` |

### Inspection & comparison
| Capability | Status | Where |
|---|---|---|
| Package explorers (`model.packages.<pkg>.inputs/.results`, normalized tables + maps) | ✅ | `package_explorer.py` facade over the `package_*` family; registry in `package_registry.py` |
| Single-model config inspector (`model.config.settings/ims/tdis/section`) | ✅ | `project/model_config.py`, `simulation/base.py` |
| **ONE `diff()` verb** — `model.diff(other)` / `group.diff()`, reference-star N-way | ✅ | `project/model_diff.py` (setup: packages + config + LAK/SFR connections), `project/model_results_diff.py` (heads, budget, per-package q, UZF, stage, MVR) |
| Grouped multi-model access (`group.hds/bud/packages`, member + Δ maps) | ✅ | `project/group/` (`ModelGroup`; facade `model_group.py`) — full single↔group symmetry |
| Diff usage guide | ✅ | `docs/model_diff_cheatsheet.md` |

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

> The sequenced roadmap for all of these is `docs/implementation_plan_2026-07.md` (rev. 4,
> 2026-07-14 — every claim re-verified against the code; read its "Resolved decisions"
> D8–D11 before starting).

1. **EVT / MAW / HFB packages** — no support yet (plan §5.2–5.4). DONE 2026-07-17:
   D8 (`mf.chd/ghb/drn/wel` each carry a real `.flopy(...)` escape hatch,
   `test_simple_list_bcs_have_a_real_flopy_escape_hatch`) and §5.1 `mf.riv` (all four
   pieces: helper + `GeoPackageSource.riv` + `riv_spec` + registry/explorers).
   CSUB is an explicit non-goal (plan §5.9). Also `mf.dis`/`mf.disu` passthroughs (§5.5).
2. **GWT/GWE integration** — package helpers (§5.3: `mf.adv/dsp/mst/ssm/cnc/...`,
   `mf.est/cnd/ctp/...`, model-type dispatch for `mf.ic/oc/disv`) AND the results tier
   (§6.0–6.2: `model.conc`/`model.temp`, grammar + hover + colors, budget terms,
   group/diff, `ConcTargets` → PEST).
3. **PRT derived cell maps + hover** (§6.3) and **PEST-IES viz upgrades**
   (field-map hover/colors, uncertainty maps, prior/posterior mosaics — §6.4).
4. **YAML/TOML serialization** — thin wrapper over the complete `to_dict`/`from_dict`
   round-trip (§5.6). *(Low effort.)*
5. **PEST raster-driven zones/pilot points**, Tikhonov/preferred-value **regularization**,
   **identifiability**, **UZF parameters** (§5.8). (IES + prior MC + PEST++ workers
   already ship via `run_ies(workers=)`.)
6. **NHDPlus SFR reader** — only piece of "SFR from streams" not already covered by `SFRBuilder`.
7. **`GridSpec.structured`/`from_geopackage` resolution** and a **raw MF6 array-file source
   reader** — both currently fail-fast / missing; niche for a Voronoi-first tool.
8. **LGR** — absent; likely not worth it for a Voronoi-first toolkit.

## Already done — do NOT rebuild
GIS-driven BCs (CHD/GHB/DRN/WEL/RCH via `GeoPackageSource` + `mf.*.gpkg`), edge/perimeter
BCs, recharge from GIS/PRISM, CRS reprojection (vector + raster), **area-weighted raster
sampling** (default), layer surfaces + reconcile + **pinch-out/idomain**, **SFR from
centerline** (`SFRBuilder`), LAK/MVR/UZF builders, full spec `to_dict`/`from_dict`
(**incl. exchanges + hooks**), multi-physics model shells + exchanges, parallel split,
ATS, `to_xugrid`, and the **unified PEST** stack — `cal.parameterize`
(constant/zone/grid/pilotpoints), pilot points on Voronoi, PESTPP-IES + prior MC +
parallel workers, and run review (`open_ies_run`).
Also done (2026-06/07): the **viz front door + unified view grammar** (`map/plot/xs` +
`mosaic/animate` everywhere, synced map mosaics), the **sectioned hover system**
(`HoverSpec`, defaults on every choropleth, LAK/SFR stage joins), the **colorscale
policy**, the **inspection/diff stack** (ONE `diff()` verb, `model.config`, ModelDiff
setup + results tiers), and full single↔group explorer symmetry.
(The old version of this file wrongly listed several of these as missing.)
