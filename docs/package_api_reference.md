# Package API reference

The vocabulary of the **package-first API** — the canonical way to build and read
MODFLOW 6 models in myflopy. Two halves:

- **Build side** — `mf.<pkg>(...)` helpers that assemble a model (`package_api.py`).
- **Read side** — the `model.packages.<pkg>.<inputs|results>.<noun>.<verb>()`
  view grammar for tables, plots, and maps.
- **Drawing** — the six plotting verbs, free in `myflopy.plot` and bound as
  `model.plot` / `vor.plot` / `stack.plot` (section C).

> **Authoritative sources (this doc is the human-readable overview of them):**
> - Build-side signatures are pinned in `tests/api_snapshot.json` (regenerate with
>   `python scripts/derive_api_snapshot.py`; the diff is the public-API change log).
> - Per-package fields/results/capabilities live in
>   `src/myflopy/modflow/mf6/package_registry.py` (the single per-package descriptor).
> - The read-side grammar is normative in `docs/view_layer_conventions.md`.
> - Deliberate gaps/deviations are in `docs/compromises_and_deferrals.md`.
>
> When these disagree with this file, they win — update this file to match.

> **Starting a new model?** `docs/model_building_cheatsheet.md` is the ordered
> walkthrough — contours → surfaces → layer stack → packages → run → results —
> with every call verified end to end. Runnable:
> `examples/mf6/contours_to_model.py`. This file is the per-symbol reference; the
> cheat sheet is the path through it.

---

## A · Build side — `mf.*`

Everything below is an attribute of the top-level `myflopy` module (imported as
`mf`). The full public surface is `mf.__all__`; lower-level builders/factories are
`mf.__engine__` (imported but omitted from `dir(mf)` so discovery points here).

### Models & simulation

| helper | builds |
|---|---|
| `mf.gwf(name, context=, packages=[...])` | a groundwater-flow model |
| `mf.gwt(...)` / `mf.gwe(...)` / `mf.prt(...)` | transport / energy / particle-tracking models |
| `mf.simulation(...)` / `mf.SimulationSpec(...)` | one MF6 simulation (tdis + solver + models) |

### Discretization, core, solver

`mf.disv`, `mf.dis`, `mf.disu`, `mf.ic`, `mf.npf`, `mf.sto`, `mf.oc`, `mf.tdis`,
`mf.ims`, `mf.ems` — each `mf.<name>(...)` returns a `PackageSpec`.

`mf.ic` / `mf.oc` / `mf.disv` / `mf.dis` / `mf.disu` are **model-kind-aware**: they
build the GWF/GWT/GWE (and, for `oc`/`disv`/`dis`, PRT) FloPy class off the model's
type, so the same helper serves a transport, energy, or particle-tracking model
(for GWT/GWE `mf.oc`, pass `concentration_filerecord=`/`temperature_filerecord=`
via `**options`; for PRT, `trackcsv_filerecord=`). `mf.npf` and `mf.sto` are
GWF-only (MF6 has no GWT/GWE variant); `mf.ic` and `mf.disu` reject PRT (MF6 has
no `ModflowPrtic`/`ModflowPrtdisu`).

**Solvers**: `mf.ims(models=…)` for GWF/GWT/GWE; **`mf.ems(models=…)`** for PRT —
PRT models are *explicit* in MF6 and are rejected under IMS6 ("Explicit models
require EMS6"). `mf.simulation(...)`'s solver default is kind-aware: one IMS per
GWF/GWT/GWE model, one EMS per PRT model.

**Grid helpers** (`mf.disv` is the default — myflopy is Voronoi/DISV-first):

- `mf.disv(nlay=, ncpl=, nvert=, vertices=, cell2d=, top=, botm=, idomain=)` —
  vertex/unstructured grid; pull the mesh from a `VoronoiGridPlus` /
  `vor.get_disv_gridprops()` and the layer arrays from a `LayerStack`.
- `mf.dis(nlay=, nrow=, ncol=, delr=, delc=, top=, botm=, idomain=)` — structured
  (row/column) grid; a passthrough for rectilinear or externally supplied models.
- `mf.disu(nodes=, nja=, top=, bot=, area=, iac=, ja=, idomain=)` — fully
  unstructured grid with connectivity you supply yourself (`ihc`/`cl12`/`hwva`,
  `nvert`/`vertices`/`cell2d` via `**options`); prefer `mf.disv` unless ingesting an
  existing DISU mesh.

`mf.dis`/`mf.disu` build and run, but the choropleth-map / cross-section / animation
viz targets **DISV/Voronoi** meshes and is not wired for structured or raw-DISU
grids (use FloPy's own plotting, or DISV).

### List boundary conditions

Each list BC offers up to three forms; the record fields are its input nouns:

| helper | record fields | `()` direct | `.gpkg(path, context=, nper=)` | `.flopy(...)` |
|---|---|:---:|:---:|:---:|
| `mf.chd` | `head` | ✓ | ✓ | ✓ |
| `mf.ghb` | `bhead`, `cond` | ✓ | ✓ | ✓ |
| `mf.drn` | `elev`, `cond` | ✓ | ✓ | ✓ |
| `mf.riv` | `stage`, `cond`, `rbot` | ✓ | ✓ | ✓ |
| `mf.wel` | `q` | ✓ | ✓ | ✓ |
| `mf.rch` | `recharge` | ✓ | ✓ | ✓ |
| `mf.evt` | `surface`, `rate`, `depth` | ✓ | ✓ | ✓ |

- `()` — direct records: `mf.riv(stress_period_data={0: [[(0, 3), 100.0, 5.0, 98.0]]})`.
- `.gpkg` — map GeoPackage features onto cells; field args accept a column name or a
  `CellSurfaceOffset` (surface-relative), e.g. `mf.drn.gpkg("seeps.gpkg", context=ctx,
  nper=12, elevation=mf.CellSurfaceOffset("cell_top", offset=-2.0))`.
- `.flopy` — raw FloPy escape hatch (the same records, no myflopy processing).

**Domain-aware builder forms** (top-active cells picked from the model domain):

- `mf.rch(context=, nper=, recharge=)` — one value per selected cell/period.
- `mf.evt(context=, nper=, rate=, depth=, surface=)` — `surface` defaults to
  `CellSurfaceOffset("model_top")` (land surface) resolved through the surface engine.

`recharge`/`rate`/`depth` accept a scalar, a per-cell sequence, a `{cellid: value}`
mapping, or a `{period: ...}` mapping.

### Advanced packages

| helper | `()` (high-level builder) | `.flopy(...)` (raw) |
|---|:---:|:---:|
| `mf.uzf` | ✓ | ✓ |
| `mf.sfr` | ✓ | ✓ |
| `mf.lak` | ✓ | ✓ |
| `mf.mvr` | ✓ | ✓ |

Movers: `mf.mvr(moves=[mf.Move(mf.MoverConnection("sfr", 0), mf.MoverConnection("lak", 0))])`.
Moved packages must be declared in the model and ordered before the mover.

### GWT / GWE transport packages

Thin factories (same shape as `mf.ic`), for solute-transport (GWT) and
energy-transport (GWE) models. Build them on a `mf.gwt(...)` / `mf.gwe(...)` model
alongside the kind-aware `mf.ic`/`mf.oc`/`mf.disv`.

| factory | package | required args |
|---|---|---|
| `mf.adv` | advection (GWT) | `scheme=` (central/upstream/tvd/utvd) |
| `mf.dsp` | dispersion (GWT) | `alh=`, `ath1=` |
| `mf.mst` | mobile storage/transfer (GWT) | `porosity=` |
| `mf.ist` | immobile storage/transfer (GWT) | `porosity=`, `volfrac=`, `zetaim=` |
| `mf.ssm` | source-sink mixing (GWT) | `sources=` (`[[pname, srctype, auxname], …]`) |
| `mf.cnc` | constant concentration (GWT, list BC) | `stress_period_data=` → `(cellid, conc)` |
| `mf.src` | mass-source loading (GWT, list BC) | `stress_period_data=` → `(cellid, smassrate)` |
| `mf.est` | energy storage/transfer (GWE) | `porosity=`, `heat_capacity_solid=`, `density_solid=` |
| `mf.cnd` | conduction/dispersion (GWE) | `ktw=`, `kts=` |
| `mf.ctp` | constant temperature (GWE, list BC) | `stress_period_data=` → `(cellid, temp)` |
| `mf.esl` | energy source loading (GWE, list BC) | `stress_period_data=` → `(cellid, senerrate)` |

The transport list BCs (`cnc`/`src`/`ctp`/`esl`) are `package_api` factories only —
no `.gpkg` and no registry entry yet. Concentration/temperature **input maps +
hover** are still Phase 6. The **results tier** (`model.conc`/`model.temp`) is
**built** as of §6.0/6.1/6.2 (2026-07-24) — see the read-side model-level reads
below. Transport budget **tables** read correctly as of 2026-07-27 — `model.bud()`
and the package result tables handle both MF6 record shapes, agree on zero-based
node ids, and label hovers `M/T` (GWT) / `E/T` (GWE) rather than `ft³/d`
(ledger 90–92, 95). The `model.budget.<term>` noun shipped 2026-07-27 (§6.1/6.2
item 3) — see "Model-level reads" below. The grouped transport fields
(`group.conc`/`group.temp`, `group.diff().conc`/`.temp`) shipped the same day
(§6.1/6.2 item 4), which also fixed `section` on `model.conc`/`model.temp` — it had
raised since those readers shipped (ledger 99). **Transport calibration** ships as of
2026-07-28: `ConcTargets` observations go straight into `cal.observe(...)`, and the
canonical model has a GWT sibling (`build_canonical_transport_calibration_demo`) whose
worked example is `examples/mf6/notebooks/canonical_07_transport_and_prt.ipynb`.
Zones for `style="zone"` come from `mf.ZoneSpec` (`from_raster` — majority vote, not mean — `from_polygons`, `from_array`); the shape each target family needs is resolved for you.

Transport properties are `parameterize` targets too: `cal.parameterize("porosity")`
(aliases `mst.porosity`, `n`) resolves the GWT sibling automatically. Calibrate it
alongside heads — porosity is absent from the flow equation, and `v = Ki/n` makes it
near-collinear with K from concentration alone. Still deferred: `TempTargets` (GWE
observations) and `dsp.alh` (ledger 103).

### PRT particle-tracking packages

A PRT model is fully declarable spec-first (§6.3A, 2026-07-24) — a coupled
GWF+PRT simulation needs **no FMI package** (flows pass through the GWF-PRT
exchange) and no raw PackageSpecs:

| factory | package | notes |
|---|---|---|
| `mf.mip` | model input (PRT) | `porosity=` required; optional `retfactor=`/`izone=` |
| `mf.prp` | particle release points (PRT) | `packagedata=` rows `(irpt, (layer, cell), x, y, z[, boundname])`; `nreleasepts`/`perioddata`/`boundnames` defaulted |
| `mf.ems` | explicit solver | `models=` — required for PRT (IMS6 rejects explicit models) |

```python
particles = mf.prt("particles", packages=[
    mf.disv(**grid), mf.mip(porosity=0.25),
    mf.prp(packagedata=[(0, (0, 0), x, y, z, "west_wells")]),   # boundname → release group
    mf.oc(trackcsv_filerecord="particles.trk.csv",
          budget_filerecord="particles.bud", saverecord=[("BUDGET", "ALL")]),
])
sim = mf.simulation(flow, particles,        # solver default: IMS for flow, EMS for particles
    exchanges=[mf.ExchangeSpec("gwfprt", mf.build_gwf_prt_exchange,
                               models=("flow", "particles"))])
```

Release-point `boundname`s are echoed (uppercased) into the track CSV's `name`
column — the release-group key the capture map reads. The high-level
`model.particle_tracking` / `PRTProject` runtime remains the sanctioned path for
*post-hoc* tracking over an already-run flow model (a separate simulation, where
FMI + grid copying ARE needed); on that path the groups come from
`PRTReleasePoints.from_cells(..., group="west_wells")` (or `from_points`), and
`PRTReleasePoints.merge(west, east)` combines grouped sets into the one PRP MF6
wants, renumbering `irpt`. Merging a **named** set with an **un-named** one is
allowed and turns BOUNDNAMES on for the whole package, so MF6 synthesizes a label
for the points that had none (`PRP000000002`) — it appears as an ordinary group in
`capture` and on the pathline legend. Name every set if you do not want that.
A release with no groups at all leaves `name` blank throughout, and the views fall
back to per-particle coloring.

**Reading a finished PRT run (§6.3B/§6.3C, 2026-07-25).** Trajectories are not a
per-cell field, so `PRTRunResults` exposes four view nouns that answer the normal
verbs (`get`/`summary`/`plot`/`map`/`mosaic`, `prt_maps.py`) — one keeping the
paths, three rolling them onto cells:

| noun | `get()` rows | `map()` draws |
|---|---|---|
| `results.pathlines` | the **normalized track records**: every raw track-CSV column plus `cell, layer, travel_time, release_group, particle` | one map polyline per particle over a base map (`base="heads"`/`None`/a `Choro`), colored per release group |
| `results.travel_time` | one per `(layer, cell)`: `layer, cell, travel_time, particle_count, min_time, max_time, release_groups` | time-of-travel choropleth (`stat="median"`, `'earth'`, optional `logscale`) |
| `results.endpoints` | the same columns (the count is the mapped one) | particle-termination counts per cell |
| `results.capture` | those columns prefixed by `group` (one per `(group, layer, cell)`) | **mosaic**, one panel per release group on a shared scale; `group=` gives one `Choro` |

```python
results = model.particle_tracking.prt(workspace=ws, release_points=rp).run()
results.pathlines.map()                                  # tracks over the water table
results.travel_time.map(stat="median", logscale=True)   # time-of-travel figure
results.capture.map()                                    # capture zones, side by side
results.travel_time.plot()                               # cumulative arrival curve
```

`pathlines` is a **view, not a frame**: `results.pathlines.get()` is the record
table and `results.track_records` the untouched CSV (what FloPy/PyVista consume).
Its `map()` returns the `Choro` carrying one `Scattermap` line per particle —
hovering a vertex reports that particle's cell, layer, elevation, release
group, and elapsed time — `plot()` draws elevation against travel time, `mosaic()` gives one panel
per release group, and `backend="mpl"` returns the existing FloPy plan view.
`max_particles=` (default 250) caps a large run by sampling **stratified across
release groups** (so a cap at or above the group count cannot drop a whole capture
zone), and says so in the title and a warning. Selectors the call cannot honor
raise rather than being dropped: base-map options on `backend="mpl"`, `per`/`layer`
on an already-built `base`, and a `base` shared across `mosaic()` panels. Group colors come from
`viz.category_colors` / `PALETTE.categorical`, so a group reads the same on the
pathline map, the arrival curve, and the capture bars.

The three per-cell maps are **time-integrated** over the whole run: they reject `per=` (and
`animate()`) rather than accept an axis they would ignore, and carry no period
footer (`result_hover(..., footer=())`). `layer=None` pools every layer into one
plan-view panel, recomputing the statistic over the pooled particles rather than
averaging per-layer values. Cells no particle reached stay blank rather than
reading as a zero travel time. `plot()` is a distribution across particles
(cumulative arrivals / count bars), not a stress-period series. `capture` needs
release-group boundnames; without them it names the fix and offers
`by="release_point"`.

### Geometry, layers, context

- `mf.ModelContext(grid=, domain=, surfaces=, dates=)` — attaches to the **model**
  (`mf.gwf(..., context=ctx)`); carries the grid/idomain/surfaces the GIS-aware helpers use.
- `mf.LayerStack(vor, top=...).add(...).build()` — declarative layer top/botm/idomain
  (compiles to `mf.LayerSurfaces` / `mf.Surface`). `vor` is optional: the layering
  can be declared before any grid exists and bound at `build(vor)`.
  `.add(name, ..., split=3)` cuts one geologic **unit** into N model layers
  (`name_1..name_N`) without moving its base; the result carries
  `units` (unit → layer indices) and `per_layer({unit: value})` for the positional
  per-layer lists `mf.npf`/`sto`/`ic` take.
- `mf.Toward(target, fraction)` — the contact a share of the way from the surface
  above down to `target`; the primitive `split=` is built on.
- `mf.CellSurfaceOffset("model_top"|"cell_top"|"cell_bottom", offset=, minimum=)`.
- Grids: `mf.VoronoiGridPlus`, `mf.TriangleGrid`, `mf.GridSpec`; sources:
  `mf.GeoPackageSource`, `mf.read_gpkg`, `mf.read_shp_gpkg`.

### Specs & lifecycle

`mf.SimulationSpec`, `mf.ModelSpec`, `mf.PackageSpec`, `mf.ExchangeSpec`;
`mf.Project(root, name=)` → `prepare_run`/`run` → `mf.Run`; `mf.load_run` /
`mf.load_mf6_run`; `mf.ModelGroup`; `mf.ParallelModelWorkflow`.

**Serialization.** Any spec round-trips to a JSON-safe dict (`to_dict()`/
`from_dict()`) — builders as importable references (including the list-BC
`functools.partial` builders), sources tagged, paths POSIX-encoded. A whole
simulation also reads/writes **YAML**:

- `sim.to_yaml()` → YAML string; `sim.to_yaml(path)` also writes the file.
- `mf.SimulationSpec.from_yaml(source)` — `source` is YAML text, a `Path`, or a
  filename `str`.
- `project.add_simulation_from_yaml(path)` — load + register in one call.

A spec whose values aren't JSON-representable (e.g. a computed numpy array baked
into a package's options) can't serialize — that raises in `to_dict` before YAML
is involved. TOML is not supported (no null type; `tomllib` is 3.11+). Example:
`examples/mf6/yaml_spec/`.

### Observations & PEST

Targets: `mf.HeadTargets`, `mf.SfrStageTargets`, `mf.SfrFlowTargets`,
`mf.LakeStageTargets`, `mf.DrnFlowTargets`. Calibration: `model.pest(name,
start_datetime=)` → `PestProject` → `cal.parameterize(...)` / `cal.observe(...)` /
`cal.forecast(...)` / `cal.build(...)` / `cal.run_ies(...)`; discovery via
`model.pest_runs` or `run.pest_runs`. **`model.pest(...)` is the front door** — it
fills `start_datetime` from TDIS and puts the run where `pest_runs` will find it.
`mf.PestProject(...)` can be constructed directly, but that is the advanced path.

**Reviewing a run** (`IesResults`, from `cal.run_ies(...)` or
`model.pest_runs[i].review()`) — `IesResults` is a flat class, not the
`noun.verb` grammar, so its figures are `plot_*` and its tables are bare nouns:

| data | figure | what it answers |
| --- | --- | --- |
| `phi`, `phi_contributions()` | `plot_phi()`, `plot_phi_distribution()`, `plot_phi_contributions()` | did misfit drop, and which group owns it |
| `conflict()` | `plot_conflict()`, `plot_prior_vs_obs()`, `plot_vs_obs()` | is the prior wide enough; does the posterior bracket the data |
| `parameters_at_bounds()` | `plot_parameters_at_bounds()` | is the prior too tight |
| `sensitivity(forecast=…)` | `plot_sensitivity(forecast=…)` | what the data informed, and what drives a forecast — ensemble-based, **not** CSS |
| `forecasts()`, `forecast(name)` | `forecast(name).plot()` | the payoff: posterior prediction spread |
| `field(target)` | `plot_field(target, stat=…)` | property patterns — `mean`/`std`/`base`/`change`/`reduction` |
| — | `plot_field_mosaic(target, stat=…)` | prior vs posterior on one shared scale (`mean`/`std` only) |
| `obs_residuals(prefix=…)` | `plot_obs_residuals(prefix=…)` | *where* the model is biased (heads/conc as points, DRN zones as cells; `prefix=` picks one family so mixed units do not share a color scale) |
| `capture_fields`, `observation_sets` | `report(path)` | what the run recorded; the headline plots bundled to HTML |

Every figure takes `backend="matplotlib"` except `plot_field_mosaic`, which composes
Plotly subplots and says so. Spatial methods need the model grid, which `run_ies()` and
`review()` attach automatically — including from `run.pest_runs`, which resolves the
model on first `review()` so that listing runs stays cheap. (A run holding several
models cannot pick one; pass `review(model=...)` there.) `report(path)` covers phi, obs,
bounds, forecasts and `plot_field` mean/change — not the reduction, mosaic or residual
maps, which are called directly.

### Engine layer (`mf.__engine__`)

The builders/factories under the facade — reach for them only with a specific reason:
`RCHBuilder`, `EVTBuilder`, `SFRBuilder`, `LAKBuilder`, `UZFBuilder`, `MVRBuilder`,
and the `<pkg>_spec` factories (`chd_spec`, `rch_spec`, `evt_spec`, …). E.g.
`mf.rch(context=...)` is `RCHBuilder(...).build()`, and `mf.evt(context=...)` is
`EVTBuilder(...).build()` (both subclasses of the shared `_ArealBuilder`).

---

### Recharge and ET, array form — `mf.rch.array` / `mf.evt.array`

| call | writes |
|---|---|
| `mf.rch(...)` / `.flopy(...)` / `.gpkg(...)` | list form, one record per cell |
| `mf.rch.array(recharge=…, context=…)` | `READASARRAYS`, one array per period |
| `mf.evt.array(rate=…, depth=…, surface=…, context=…)` | the same for ET |

Array when the boundary covers the model; list when it is sparse or needs
boundnames. `irch`/`ievt` default to `"top_active"`, derived from `context.domain`
— without one MODFLOW 6 applies the flux to layer 1 and silently skips any column
whose layer 1 is inactive. Segmented ET has no array form and `nseg` raises.

### Horizontal flow barriers — `mf.hfb`

A barrier sits on the **face between two cells**, so MODFLOW 6 addresses it as a
cell *pair*. FloPy takes any pair at all and never checks it is a real connection;
the run aborts instead. These resolve geometry and validate first.

| call | from |
|---|---|
| `mf.hfb(pairs=…, hydchr=…, layers=…, context=…)` | explicit cell pairs |
| `mf.hfb.line(trace, context=…, hydchr=…)` | a fault trace → the faces it crosses |
| `mf.hfb.gpkg(path, context=…, hydchr="hydchr")` | a GeoPackage line layer |
| `mf.hfb.enclose(polygon, context=…, hydchr=…)` | a watertight cutoff wall around a region |
| `mf.hfb.flopy(stress_period_data=…)` | the FloPy-native form, unvalidated |

`hydchr` is a hydraulic **characteristic**, 1/T — barrier K divided by barrier
thickness, not K. `0` is impermeable; a negative value is a conductance multiplier.

`.enclose` is separate from `.line` on purpose: a ring passes *through* cells, so
the faces it crosses are **not** watertight. The cut between the enclosed cells and
their neighbours is.

Grid verbs underneath: `vor.barrier_faces(geom)`, `vor.enclosed_faces(polygon)`,
`vor.shared_face(a, b)`.

### Importing an existing MODFLOW-USG model

| call | returns | notes |
|---|---|---|
| `mf.read_usg(nam, gsf=, crs=, nper=, require_layered=, read_boundaries=)` | `UsgModel` | A USG `DISU` has no coordinates, so the `.gsf` is required for anything but inspection. |
| `usg.report()` | `str` | What converted, what was approximated, what was omitted — with counts. |
| `usg.validate()` | `list[Finding]` | What MODFLOW 6 will **reject**, before running. `finding.blocks_run`. |
| `usg.to_mf6(name, crs=, newton=, start_date_time=, complexity=, include=, fix_for_mf6=, local_origin=)` | `SimulationSpec` | DISV. Keep `local_origin=True`. |
| `usg.surfaces()` / `usg.boundary_frame(ftype)` / `usg.cln_polygons()` | `Surface`s / GeoDataFrame / GeoDataFrame | Grid-independent — these survive a change of grid, node numbers do not. |

`UsgModel` also exposes the model reshaped to `(nlay, ncpl)`: `top`, `botm`, `idomain`,
`strt`, `k`, `k33`, `ss`, `sy`, `icelltype`, `thickness`, `uppermost_active`.

## B · Read side — the `noun.verb` grammar

```
model.packages.<pkg>.<inputs|results>.<noun>.<verb>()
```

Every **noun** is an object; every object answers the same **verbs**. Different plot
*types* are different verbs — `plot()` is not a dispatcher.

### The verbs

| verb | available on | returns |
|---|---|---|
| `get()` | every noun | the DataFrame |
| `summary()` | every noun | a compact digest frame |
| `plot()` | every noun | the noun's **one** canonical chart (`viz.Fig`) † |
| `map()` | **spatial** nouns | choropleth on the grid |
| `section()` | spatial nouns | cross-section |
| `mosaic()` | spatial nouns | small-multiples panel |
| `animate()` | spatial nouns | animation over stress periods |
| `long()` / `wide()` / `stack()` | cell-budget fields | frame reshapes |
| `__call__(per=)` | period-bound nouns | a fresh view rebound to `per` |

† Exception: `lak.results.q.budget.plot()` currently returns a **matplotlib** figure
(the `viz.Fig` conversion is deferred — ledger 52). Derived-table nouns (`profile`,
`budget`) are not spatial, so they have no `map`/`section` (ledger 17); use the field
noun's own `map()` (e.g. `sfr.results.q.map()`).

### Nouns per package

**List BCs** (`chd/ghb/drn/riv/wel/rch/evt`):
- `inputs.<field>` — the record fields above (e.g. `riv.inputs.stage`,
  `evt.inputs.rate`, `rch.inputs.recharge`) — each a spatial field noun.
- `results.q` — the groundwater exchange (spatial). The `.get()` column names its
  reference frame: **`q_gwf`** (aquifer `.cbc` record; negative = the feature gains).

**SFR** (`model.packages.sfr`):
- `results.q` — stream↔aquifer exchange (spatial; `q_gwf`).
- `results.stage` — reach stage (spatial).
- `results.profile` — the merged multi-field reach profile (`get`/`summary`/`plot`).
- `results.q.profile` / `results.stage.profile` — single-field reach profiles
  (`get`/`summary`/`plot`). *The frame is `…profile.get()`, not `…profile()`.*
- `budget.<term>` — budget terms: `flow_ja_face`, `gwf`, `ext_inflow`, `ext_outflow`,
  `rain`, `runoff`, `evaporation`, `storage`, `from_mvr`, `to_mvr`, `mvr`,
  `stream_fluxes`, `auxiliary`, `types`.

**LAK** (`model.packages.lak`):
- `results.q` — lake↔aquifer exchange (spatial; column **`q_lake`**, feature-referenced,
  negative = the lake loses).
- `results.q.budget` — exchange summarized by connection type (`get`/`summary`/`plot`).
- `results.stage`, `results.stage_change`, `results.connections`.
- `budget.<term>` — `gwf`, `flow_ja_face`, `ext_inflow`, `ext_outflow`, `rainfall`,
  `runoff`, `evaporation`, `storage`, `from_mvr`, `to_mvr`, `mvr`, `withdrawal`,
  `constant`, `lake_fluxes`, `auxiliary`, `types`.

**UZF** (`model.packages.uzf`): `results.gwrch`, `results.sat`, `results.fields`.

**PRT** (a finished run's `PRTRunResults`, not `model.packages`): `results.pathlines`
(trajectories — `map` draws polylines, not cells), `results.travel_time`,
`results.endpoints`, `results.capture`. These are derived, time-integrated views, so
they reject `per=`; `pathlines` also refuses `section`/`animate` (a trajectory is not a
per-cell field to slice) and `travel_time`/`endpoints`/`capture` refuse `animate`.
The raw MF6 track table is `results.track_records`. Full detail in §A above.

**Model-level reads:** `model.hds` (GWF heads explorer — `get`/`summary`/`array`/
`map`/`section`/`mosaic`/`animate`), and its transport twins **`model.conc`** (GWT
concentration) and **`model.temp`** (GWE temperature) — the same dependent-variable
grammar, reading the `.ucn` binary via the shared `DependentVariableFile` base. The
readers are **kind-gated**: `.hds` on a transport model (or `.conc` on a flow model)
raises a clear error. `model.outputs.<pkg>.bud` (raw budget accessor),
`model.targets.<family>` (`compare`/`stats`/`calibration_plot`), `model.pest_runs`.

**`model.budget.<term>`** — every term in the model's *own* budget file, as a spatial
noun (`get`/`summary`/`plot`/`map`/`section`/`mosaic`/`animate`). Terms are **discovered
from the file**, not declared, so they follow the model's kind and packages:

```python
model.budget.types                     # the MF6 record names actually present
model.budget.source_sink_mix.get()     # GWT/GWE: the SSM term, per cell
model.budget.storage_aqueous.map()     # GWT storage (GWE spells it storage_cellblk)
model.budget.drn.summary()             # GWF: any boundary package
model.budget.sto_ss.get(per=0)         # ...and terms with no package accessor
model.budget["SOURCE-SINK MIX"]        # indexing takes either spelling
```

Available on **every** model kind, not just transport: the plumbing is kind-neutral,
and `STO-SS`/`STO-SY`/`DATA-SPDIS`/`DATA-SAT` have no package accessor at all.
`FLOW-JA-FACE` appears in the namespace but *refuses* `get()` — it is indexed by cell
connection, not by cell (ledger 91/96).

Distinguish it from its three neighbours, which are genuinely different things:

| spelling | what it reads |
|---|---|
| `model.budget.<term>` | the model budget file, per cell, one noun per term |
| `model.bud(pkg)` | the same file, legacy compatibility wrapper, raw frames |
| `model.budget_cumulative` / `_incremental` | the **listing** file: whole-model totals |
| `model.packages.<pkg>.budget.<term>` | the **package-output** file — a different file, feature-first node layout |

The column is MF6's raw `q`, unlike `model.packages.<pkg>.results.q`, which renames it
`q_gwf`/`q_lake` to carry the reference frame. The values are identical; only the name
differs (ledger 97).

**Grouped reads:** `group.hds` / `group.conc` / `group.temp` are the multi-model twins
of the model-level readers, with the same verbs plus `compare()`:

```python
group.conc.get()                     # aligned concentration, one block per member
group.conc.compare()                 # vs the reference: conc / reference_conc / diff
group.conc.map(); group.conc.plot(); group.conc.section(line=line)
group.diff().conc.map("variant")     # Δconc choropleth -- the ONE public diff route
group.diff().temp.summary()          # max/mean |Δtemp|, RMSE, argmax cell/layer
```

All three are the same class configured by five attributes (`_GroupFieldView`), so a
verb added to one is available to all. `diff` is always **model − reference**. Every
member of a `group.conc` must be a GWT model; the kind-gated reader supplies the error
otherwise (ledger 100).

### Deprecated spellings (D12)

The following resolve through `__getattr__` with a `DeprecationWarning` and stay out
of `dir()`/completion; they return exactly what they used to. Prefer the noun form:

| deprecated | use instead |
|---|---|
| `sfr.results.long_profile()` / `plot_long_profile()` | `sfr.results.profile.get()` / `.plot()` |
| `sfr.results.q.plot_profile()` / `sfr.results.stage.plot_profile()` | `…profile.plot()` |
| `lak.results.q.budget_summary()` / `plot_budget()` | `lak.results.q.budget.get()` / `.plot()` |

---

---

## C · Drawing — the plotting verbs

**Highlighting a subset of cells** — `select=` is an option on `map` at every
scope (free `mf.plot.map`, `model.plot.map`, `vor.plot.map`, and every
`…inputs/results.map` leaf):

```python
model.plot.map(layer=0, select=cells)                    # outline (default)
vor.plot.map(select="all_streams")                       # a registered region
model.plot.map(select="wetland.gpkg")                    # a vector file
vor.plot.map(values, select=cells, select_style="dim")   # fade the rest instead
vor.plot.map(select=cells, select_color="#0072B2")       # avoid a red colorscale
```


Six verbs cover every picture. They exist as free functions and bound to the
objects, and the bound form calls the free one, so the two cannot diverge.

```python
from myflopy import plot

plot.map(model, layer=0)      # or  model.plot.map(layer=0)
plot.grid(vor)                # or  vor.plot.grid()   -- or just vor.plot()
```

### The verbs

| verb | draws | notes |
|---|---|---|
| `map(source, values=, ...)` | plan view | contours, locations, hillshade, pathlines are **options**, not verbs |
| `section(source, ...)` | vertical slice | results through a model, geometry through a grid |
| `surface(source, ...)` | 3-D height field `z(x, y)` | Plotly |
| `grid(source, backend=)` | the mesh itself | `"plotly"` flat 2-D (no CRS needed), `"vtk"` the 3-D volume |
| `mosaic(panels, ncols=)` | many pictures, one figure | **combinator** — takes pictures, not a subject |
| `animate(frames, backend=)` | frames in sequence | **combinator**; `"plotly"` live, `"png"` rasterized |

### Which scope answers which

A scope answers the verbs it can *know*: a bare grid has no results, a layer
stack has no time.

| scope | verbs |
|---|---|
| `myflopy.plot` | all six |
| `model.plot` | all six |
| `vor.plot` | `map` `section` `grid` |
| `stack.plot` | `map` `section` `surface` `grid` |

Pinned by `tests/test_plot_vocabulary.py`, which fails naming both the scope and
the stray verb if one drifts.

### The parameters are in the signature, not just the docstring

Every verb spells its parameters out in the `def`. `plot.map` names 35 of them;
`model.plot.map` names the same 35. That duplication is deliberate and is the
only thing that works: **PyCharm and Pylance are static analyzers** — they read
the `def` line and never execute the module, so `__signature__`, `functools.wraps`
and a beautifully written docstring reach `help()` and reach no editor at all.
A `**kwargs` forwarder shows the caller nothing.

Each verb still keeps a `**kwargs` tail, but only where the tail is *genuinely*
open — `map`'s rides through to the `go.Choroplethmap` trace, whose parameter
names belong to Plotly and are validated at render time.

Three tests keep the copies honest, and all three fail by name:

| test | catches |
|---|---|
| `test_a_named_default_matches_the_link_that_owns_it` | a restated default that drifted from the function that owns it |
| `test_the_bound_model_verb_mirrors_the_free_one` | a parameter added to `plot.map` and forgotten on `model.plot.map` |
| `test_every_signature_parameter_is_documented` | accepted-but-undocumented **and** documented-but-not-accepted |

Narrower scopes may offer *fewer* parameters, never other ones: `vor.plot.map`
drops `per`/`layer`/`type` because a bare grid has no results to read. Anything
dropped still rides the `**kwargs` tail, so narrowing a signature never narrowed
what already worked.

**One deliberate exception.** The namespace-level `field=` sugar —
`packages.lak.results.map(field="stage")` — stays `(*args, field=None, **kwargs)`.
It dispatches to accessors with genuinely incompatible signatures (a lake's
`connections` map has no `per=`; `DrnInput.map` takes `per` positionally), so a
single merged signature could only be produced by lying about one of them. Call
the leaf — `packages.lak.results.stage.map(` — when you want completion; every
leaf verb is explicit.

### Everything returned is a Picture

```python
picture = model.plot.map(layer=0)
picture                       # renders itself in Jupyter
picture.fig                   # the figure, to adjust before display
picture.show()
picture.save("heads.png")     # suffix picks the format
picture.html("heads.html")    # standalone page
```

There is never a trailing `.plot()`. Three renderers back this: most pictures are
Plotly; `MplPicture` (Matplotlib) offers `.axes` and `VtkScene` (PyVista) offers
`.scene` instead of `.fig`, and their `.fig` raises saying so.

### Standalone HTML

Two artifacts for two jobs. `plot.animate(frames).html(path)` writes an
interactive page — for a choropleth it ships the cell geometry once and restyles
per frame rather than re-embedding the grid. The `mf.export_*_slider_html`
functions instead pre-render frames to PNG and page through them: their size is
independent of cell count, and they carry the machinery a long export needs
(`resume=`, `progress=`, external frame directories, a concurrent-writer lock).
They render through FloPy's `PlotMapView`, so they are a *different picture*,
not a second spelling.

> Runnable tour: `examples/mf6/notebooks/plotting_vocabulary_tour.ipynb`.
> Normative rules: `docs/view_layer_conventions.md`.

---

## Worked example

```python
import myflopy as mf

# build
ctx = mf.ModelContext(grid=vor, domain=idomain, surfaces=vor.gdf_topbtm)
gwf = mf.gwf("valley", context=ctx, packages=[
    mf.disv(...), mf.ic(strt=...), mf.npf(k=...), mf.sto(...),
    mf.rch(context=ctx, nper=6, recharge=6.0e-4),         # domain-aware builder
    mf.evt(context=ctx, nper=6, rate=2.0e-3, depth=2.5),  # surface = model top
    mf.riv.gpkg("river.gpkg", context=ctx, nper=6),       # from GeoPackage
    mf.lak(context=ctx, ...), mf.sfr(context=ctx, ...), mf.mvr(...),
    mf.oc(...),
])
sim = mf.SimulationSpec("base", models=(gwf,), packages=(mf.tdis(...), mf.ims(...)))
run = mf.Project(root).run("base", sim)               # build + run MF6

# read
model = run.model()
model.packages.sfr.results.profile.plot(per=5)      # merged reach profile
model.packages.sfr.results.q.profile.get()          # q_gwf along the stream
model.packages.lak.results.q.budget.get(per=5)      # exchange by connection type
model.packages.riv.results.q.map(per=5)             # exchange choropleth
model.packages.sfr.budget.flow_ja_face.get()        # in-channel routing flow
```
