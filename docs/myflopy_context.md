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
| **Legacy OO API** | `src/myflopy/modflow/mf6/*.py` — `simplemodel.py`, `simulation/`, `boundaries.py`, `drn.py`/`ghb.py`/`chd.py`, `recharge.py`, `lakes.py`, `sfr.py`, `mvr.py`, `uzf.py` | `SimulationBase`, `Boundaries`, builder classes (`SFRBuilder`, `LAKBuilder`, `RCHBuilder`, …) |

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
| GeoPackage feature → cell pipeline | ✅ | `geopackage.py` `GeoPackageSource` (`.chd/.ghb/.drn/.riv/.wel/.rch/.evt`) |
| CRS reprojection of vector sources | ✅ | `boundaries.py` `Boundaries.gdf` (`to_crs`) |
| CRS reprojection of rasters | ✅ | `grid/surfaces.py` (fixed 2026-06) |
| Import existing MODFLOW layer surfaces | 🟡 | `LayerStack.from_modflow` / `LayerSurfaces.from_modflow` read top/botm/idomain from a built MF6 model; a raw array-file *source* reader is still missing |

### Boundary conditions
| Capability | Status | Where |
|---|---|---|
| CHD/GHB/DRN/RIV/WEL/RCH/EVT/UZF package builders | ✅ | `package_api.py` (`_CHDPackage`…`_MVRPackage`), `advanced.py` `*_spec` |
| **BCs auto-built from GIS** (polygon/line) | ✅ | `mf.ghb.gpkg(...)`, `mf.drn.gpkg(...)`; `GeoPackageSource`; legacy `Boundaries` |
| Perimeter / edge-cell BCs (CHD/GHB on grid edge) | ✅ | `Boundaries.edge_intersections`, `edges_only=True` on `.gpkg()` builders |
| Cells ordered along a line (for line BCs / SFR) | ✅ | `Boundaries.sorted_cells_along_line` |
| Recharge from GIS + PRISM precip scaling + area scaling | ✅ | `recharge.py` `RCHBuilder`, `PrismPrecipScaling`, `Boundaries.shp_to_vor_poly_scale` |
| Dedicated RIV builder | ✅ | `mf.riv` (`_RIVPackage`: `()`/`.gpkg`/`.flopy`), `riv_spec`, `GeoPackageSource.riv`, registry entry (`stage/cond/rbot` earth, q RdBu) — added 2026-07-17 — in the canonical model as the un-routed lake-outlet river (2026-07-18) |
| Dedicated EVT builder | ✅ | `mf.evt` (`_EVTPackage`: `()`/`.gpkg`/`.flopy`, `nseg=1` default), `evt_spec`, `GeoPackageSource.evt`, registry entry (`surface/rate/depth` earth, q RdBu) — added 2026-07-17 — in the canonical model on the valley walls, disjoint from UZF (2026-07-18). File-less whole-domain builder `mf.evt(context=, nper=, rate=, depth=)` via `EVTBuilder` (shared `_ArealBuilder` base with `RCHBuilder`), surface default `model_top` — added 2026-07-21 (§4.7.5) |

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
| GWT/GWE **package helpers** (adv/dsp/mst/ist/ssm/cnc/src; est/cnd/ctp/esl) | ✅ | `package_api.py` factories build `ModflowGwt*`/`ModflowGwe*` (plan §5.3B); `mf.ic/oc/disv/dis/disu` dispatch the FloPy class on the model kind (§5.3A/§5.5) |
| GWT/GWE **results tier** — reader + maps | ✅ | `model.conc`/`model.temp` (`headsplus.py` `DependentVariableFile` base + `ConcResults`/`TempResults`); full grammar (`get/summary/array/map/section/mosaic/animate`), field hover (`conc_hover`/`temp_hover`), `'earth'` colorscale; kind-gated (§6.0/6.1/6.2, 2026-07-24) |
| GWT/GWE results — budget **tables** | ✅ | Transport budget terms read correctly through `build_budget_result_table` and `model.bud(...)` as of 2026-07-27: transport models save flows by default, imeth=1 full-array terms (`STORAGE-AQUEOUS`/`STORAGE-CELLBLK`) build a real per-cell table, both paths agree on zero-based nodes, hover units follow the model kind (`M/T`/`E/T`), and `model.bud("ssm")` resolves (ledger 90–92, 95) |
| GWT/GWE results — budget **noun** | ✅ | `model.budget.<term>` (`package_results.py` `ModelBudgetNamespace` → `CellBudgetResultsExplorer`): every term in the model's own budget file as a spatial noun with the full verb set. Terms **discovered** from the file, so they follow kind + packages (`storage_aqueous` on GWT, `storage_cellblk` on GWE). Available on all kinds, incl. GWF terms with no package accessor (`sto_ss`, `data_spdis`) — §6.1/6.2 item 3, 2026-07-27 |
| GWT/GWE results — **group + diff** | ✅ | `group.conc`/`group.temp` + `group.diff().conc`/`.temp` (`project/group/conc.py`, `temp.py`). NOT clones: `GroupHeads` was refactored onto a shared `_GroupFieldView`, so all three kinds are one implementation with five class attributes each (§6.1/6.2 item 4, 2026-07-27). Needed `SimulationBase.all_conc`/`.all_temp` (kind-gated, unlike `all_heads`) |
| GWT/GWE **PEST calibration** | ✅ | `ConcTargets` + `ConcObservationSpec` + `prepare_conc_observations` + a `.ucn` forward-run post-processor: a control file history-matched against CONCENTRATIONS builds and runs end to end, recovering K (§6.1/6.2 item 5, 2026-07-28). `TempTargets` (GWE) still deferred |
| Canonical **transport fixture** | ✅ | `canonical_transport.py` — a GWT sibling on the canonical simulation (flat, so PEST works), `build_canonical_transport_calibration_demo`, and `examples/mf6/notebooks/canonical_07_transport_and_prt.ipynb` (§6.1/6.2 item 6, ledger 101/103) |
| MF6 PRT (release points, run, pathlines, 3-D scenes) | ✅ | `prt.py` (`PRTProject`, `PRTRunResults`, `open_prt_run`), `model.particle_tracking`; raw CSV on `results.track_records` |
| PRT **spec-first declarability** (mip/prp/ems + exchange, no FMI) | ✅ | `mf.mip`/`mf.prp`/`mf.ems` factories + prt6 dis/disv/oc dispatch; `mf.simulation` solver default is kind-aware (IMS vs EMS) — §6.3A, 2026-07-24 |
| PRT derived cell maps (travel-time / endpoints / capture choropleths) | ✅ | `prt_maps.py` — `results.travel_time` / `.endpoints` / `.capture` nouns (get/summary/plot/map/mosaic), release groups via PRP boundnames — §6.3B, 2026-07-25 |
| PRT plotly pathline map + `pathline_hover` | ✅ | `results.pathlines` view (`prt_maps.PRTPathlineView`): one `Scattermap` polyline per particle over `base="heads"`/`None`/a `Choro`, per-vertex `pathline_hover`, `plot()` = elevation vs travel time, `mosaic()` per release group, `backend="mpl"` = the FloPy plan view — §6.3C, 2026-07-25 |
| Map overlays composed into mosaics | ✅ | `Choro.add_overlay`/`overlay_traces` + `viz.mosaic` copying them (was silently dropping contours/locs) — §6.3C |
| Category colors (release groups, zones) | ✅ | `viz.PALETTE.categorical` + memoized `viz.category_colors` — one color per name across every figure |
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
| Plotting front door (`viz.Fig/subplots/mosaic/mpl_axes/PALETTE/category_colors`, figs backend) | ✅ | `viz.py` — import every figure from here, not raw plotly/matplotlib. figs is vendored (`myflopy/_vendor/figs`, synced via `scripts/sync_vendored_figs.py`): viz.py imports the live figs first and falls back to the snapshot, so pip installs work without the local figs project |
| Unified view grammar: `map/plot/section` + `mosaic/animate` on every leaf (model / group / diff), `backend="plotly"\|"mpl"` | ✅ | `package_plotting.py` (`SpatialView`, composers), `viz.mosaic` |
| **Derived-table views** (`<pkg>.<inputs\|results>.<noun>.<verb>` for merged tables, not just mappable fields) | ✅ | **`docs/view_layer_conventions.md` is the normative rule.** Reference impl: `SfrProfileView` (`sfr.results.profile`). Before 2026-07-18 derived tables sat OUTSIDE the grammar as `foo()`/`plot_foo()` pairs — that gap is what let two notebooks hand-roll matplotlib over a built-in |
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
| Grouped multi-model access (`group.hds/conc/temp/bud/packages`, member + Δ maps) | ✅ | `project/group/` (`ModelGroup`; facade `model_group.py`) — full single↔group symmetry, transport fields included as of 2026-07-27 |
| Diff usage guide | ✅ | `docs/model_diff_cheatsheet.md` |

### Serialization
| Capability | Status | Where |
|---|---|---|
| `to_dict`/`from_dict` round-trip for all specs | ✅ | `specs.py` — Package/Grid/Model/SimulationSpec **incl. inter-model exchanges, post-build hooks, and list-BC `functools.partial` builders** (§5.6A; importable builders; lambdas/closures fail loud) |
| Workspace pickling / native MF6 reload | ✅ | `workspace.py`, `project/` (`load_mf6_run`) |
| **YAML** config file wrapper | ✅ | `specs_io.py` — `SimulationSpec.to_yaml`/`from_yaml` + `Project.add_simulation_from_yaml` (§5.6, PyYAML safe mode); example `examples/mf6/yaml_spec/`. **TOML** deferred (no null type; `tomllib` is 3.11+) |

### PEST / pyemu (`modflow/mf6/pest/`)
| Capability | Status | Where |
|---|---|---|
| `PestProject` orchestration (native PstFrom build + forward run) | ✅ | `pest/project.py`, `pest/forward_run.py` — **front door is `model.pest(name, start_datetime=...)`**, which fills `start_datetime` from TDIS and defaults the workspace to `<model ws>.pest/<name>` — a SIBLING of the model directory, outside the tree `PstFrom` copies (ledger 107); the legacy `<model ws>/pest` root is still discovered. Direct `PestProject(...)` construction is the **advanced** path: `start_datetime` is mandatory and you own where it lands |
| **Unified** `cal.parameterize(target, style=...)` | ✅ | `pest/project.py` + `pest/native_parameters.py`; targets `k/k33/recharge/chd/ghb/drn/wel/uzf.vks` plus transport `porosity` (resolves the GWT sibling via `pest/model_lookup.py`), styles `constant/zone/grid/pilotpoints` (`pilotpoints` is flow-only — it interpolates against NPF K) |
| Pilot points on Voronoi (IDW; pyEMU's own are unusable on DISU) | ✅ | `pest/pilot_points.py`; `pp_space` net or explicit `pp_points` |
| Observations: head, lake stage, SFR stage/flow, DRN flow | ✅ | `pest/observations.py` (`cal.observe` / `cal.forecast`) |
| PESTPP-IES + prior Monte Carlo + parallel workers | ✅ | `pest/project.py` `run_ies(workers=)`, `prior`, `draw_prior` |
| Results / review / reopen a run | ✅ | `pest/ies.py` (`open_ies_run`, `IesResults`), `pest/runs.py` (`find_pest_runs`) — `model.pest_runs` / `run.pest_runs` → `PestRunHandle.review()`. Both attach the model so spatial maps work (`run.pest_runs` lazily, so listing does not load the simulation); a run with several models cannot pick one — pass `review(model=...)` |
| Parameter-field maps with hover + per-stat policy colors | ✅ | `pest/ies.py` `plot_field` (`_field_map_policy`: `change` → log-centered diverging, red = reduced; `mean`/`base` → earth+log; `std` → earth linear) + `parameter_field_hover` (§6.4A) |
| Uncertainty-reduction map (`did the data inform this region`) | ✅ | `pest/ies.py` `plot_field(stat="reduction")` — `1 - post_sd/prior_sd`, anchored 0–1, falls back to the data range when a posterior spread grew (§6.4B) |
| Prior-vs-posterior field mosaic (one shared scale, synced views) | ✅ | `pest/ies.py` `plot_field_mosaic` over `viz.mosaic(colorbar=)`; `mean`/`std` only — the only stats with a separate prior and posterior form (§6.4B). Plotly only |
| Observation residual map (heads + conc as points, DRN zones as cells) | ✅ | `pest/ies.py` `obs_residuals` / `plot_obs_residuals`, joined to the saved `*_target_locations.gpkg`/`.csv` + the set's own `mapping_file`; dispatch is on the set's declared `geometry` (`points`/`zones`), not its kind; `prefix=` selects one family so mixed units do not share a color scale (warns otherwise); excludes forecasts and unmeasured times; both backends (§6.4B). Lake/SFR record only a lake/reach number — ledger 76 |
| Geostats for grid / pilot-point priors | ✅ | `pest/geostats.py` (`ExpGeoStruct`, `build_geostruct`) |
| Zones for `style="zone"` from raster / polygon / array | ✅ | `pest/zones.py` `ZoneSpec` (`from_raster` majority-vote, `from_polygons`, `from_array`) + `resolve_zone_array`, which normalizes per family — pyEMU wants per-cell for array targets and `(nlay, ncpl)` for list targets and rejects the other cryptically (ledger 118). The earlier "polygon zones exist" claim here was FALSE: `zone_array` had one occurrence in all of `src/`, an unvalidated pass-through |
| Ensemble sensitivity / data worth | ✅ | `pest/ies.py` `sensitivity()` / `plot_sensitivity()` — `learned` (1 − post sd/prior sd, per group) and, with `forecast=`, the ensemble correlation that drives it. No jacobian, so no new run mode; NOT composite scaled sensitivity, and noisy below ~1/√n_reals (ledger 120) |
| Jacobian-based **identifiability**/CSS/FOSM (`Schur`/`ErrVar`); Tikhonov **regularization** | ❌ | regularization is a PESTPP-GLM concept and myflopy ships no GLM runner — PESTPP-IES ignores prior-information equations and rejects a version=2 regularized pst outright (ledger 117); identifiability needs a jacobian IES does not produce |

---

## The genuinely short gap list (verify each before building)

> The sequenced roadmap for all of these is `docs/implementation_plan_2026-07.md` (rev. 4,
> 2026-07-14 — every claim re-verified against the code; read its "Resolved decisions"
> D8–D11 before starting).

1. **MAW / HFB packages** — no support yet (plan §5.4). DONE 2026-07-17: D8
   (`mf.chd/ghb/drn/wel` each carry a real `.flopy(...)` escape hatch,
   `test_simple_list_bcs_have_a_real_flopy_escape_hatch`), §5.1 `mf.riv`, and §5.2
   `mf.evt` (each with all four pieces: helper + `GeoPackageSource` resolver +
   `*_spec` + registry/explorers). CSUB is an explicit non-goal (plan §5.9).
   `mf.dis`/`mf.disu` passthroughs **DONE 2026-07-23** (§5.5: model-kind dispatch,
   PRT excluded, viz stays DISV-only).
2. **GWT/GWE integration** — package helpers (§5.3 **DONE 2026-07-22**: `mf.adv/dsp/
   mst/ist/ssm/cnc/src` (GWT) + `mf.est/cnd/ctp/esl` (GWE); model-type dispatch for
   `mf.ic/oc/disv/dis/disu` builds the gwf/gwt/gwe FloPy class off the model kind).
   Results tier: reader + maps **DONE 2026-07-24** (§6.0/6.1/6.2: `model.conc`/
   `model.temp` grammar + hover + `'earth'` colorscale, kind-gated). Transport
   budget **tables** read correctly as of 2026-07-27 (ledger 90–92, 95 — five
   defects, from a zero-byte `.cbc` to silently empty tables), and the
   **`model.budget.<term>` noun** shipped the same day (§6.1/6.2 item 3:
   runtime-discovered terms, full spatial verb set, every model kind).
   `GroupConc`/`GroupTemp` + `diff().conc`/`.temp` shipped 2026-07-27 (§6.1/6.2
   item 4), on a shared `_GroupFieldView` rather than as clones — which also
   fixed `section` on the transport readers (ledger 99). Still deferred:
   `ConcTargets`/`TempTargets` → PEST.
3. **PRT**: **DONE 2026-07-25** — derived cell maps (§6.3B: `results.travel_time`/
   `.endpoints`/`.capture` nouns + release groups) and the plotly pathline map +
   `pathline_hover` (§6.3C: `results.pathlines`). **PEST-IES viz upgrades**: field-map
   hover + per-stat policy colors **DONE 2026-07-25 (§6.4A)**; uncertainty-reduction
   map, prior/posterior mosaic and the residual map **DONE 2026-07-26 (§6.4B)**;
   forecast-less `settings`/`report()` fix (ledger 73), the last `ies.py` color
   literals, the field-fixture buildout and the preferred-API docs
   **DONE 2026-07-26 (§6.4C)** — **§6.4 is complete**.
4. **YAML serialization DONE 2026-07-23 (§5.6)** — `SimulationSpec.to_yaml`/`from_yaml`
   + `Project.add_simulation_from_yaml` in `specs_io.py`; §5.6A first made the list-BC
   `functools.partial` builders round-trip. **TOML still deferred** (no null type;
   `tomllib` is 3.11+ vs the `>=3.10` floor — ledger 54).
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
Also done (2026-06/07): the **viz front door + unified view grammar** (`map/plot/section` +
`mosaic/animate` everywhere, synced map mosaics), the **sectioned hover system**
(`HoverSpec`, defaults on every choropleth, LAK/SFR stage joins), the **colorscale
policy**, the **inspection/diff stack** (ONE `diff()` verb, `model.config`, ModelDiff
setup + results tiers), and full single↔group explorer symmetry.
(The old version of this file wrongly listed several of these as missing.)
