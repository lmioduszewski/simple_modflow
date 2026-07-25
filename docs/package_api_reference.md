# Package API reference

The vocabulary of the **package-first API** — the canonical way to build and read
MODFLOW 6 models in myflopy. Two halves:

- **Build side** — `mf.<pkg>(...)` helpers that assemble a model (`package_api.py`).
- **Read side** — the `model.packages.<pkg>.<inputs|results>.<noun>.<verb>()`
  view grammar for tables, plots, and maps.

> **Authoritative sources (this doc is the human-readable overview of them):**
> - Build-side signatures are pinned in `tests/api_snapshot.json` (regenerate with
>   `python scripts/derive_api_snapshot.py`; the diff is the public-API change log).
> - Per-package fields/results/capabilities live in
>   `src/myflopy/modflow/mf6/package_registry.py` (the single per-package descriptor).
> - The read-side grammar is normative in `docs/view_layer_conventions.md`.
> - Deliberate gaps/deviations are in `docs/compromises_and_deferrals.md`.
>
> When these disagree with this file, they win — update this file to match.

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
below; still deferred: GWT/GWE budget views, `GroupConc`/`GroupTemp` group/diff,
and `ConcTargets`/`TempTargets` PEST observations (ledger 56).

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
  (compiles to `mf.LayerSurfaces` / `mf.Surface`).
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
`model.pest_runs`.

### Engine layer (`mf.__engine__`)

The builders/factories under the facade — reach for them only with a specific reason:
`RCHBuilder`, `EVTBuilder`, `SFRBuilder`, `LAKBuilder`, `UZFBuilder`, `MVRBuilder`,
and the `<pkg>_spec` factories (`chd_spec`, `rch_spec`, `evt_spec`, …). E.g.
`mf.rch(context=...)` is `RCHBuilder(...).build()`, and `mf.evt(context=...)` is
`EVTBuilder(...).build()` (both subclasses of the shared `_ArealBuilder`).

---

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
| `xs()` | spatial nouns | cross-section |
| `mosaic()` | spatial nouns | small-multiples panel |
| `animate()` | spatial nouns | animation over stress periods |
| `long()` / `wide()` / `stack()` | cell-budget fields | frame reshapes |
| `__call__(per=)` | period-bound nouns | a fresh view rebound to `per` |

† Exception: `lak.results.q.budget.plot()` currently returns a **matplotlib** figure
(the `viz.Fig` conversion is deferred — ledger 52). Derived-table nouns (`profile`,
`budget`) are not spatial, so they have no `map`/`xs` (ledger 17); use the field
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
they reject `per=`; `pathlines` also refuses `xs`/`animate` (a trajectory is not a
per-cell field to slice) and `travel_time`/`endpoints`/`capture` refuse `animate`.
The raw MF6 track table is `results.track_records`. Full detail in §A above.

**Model-level reads:** `model.hds` (GWF heads explorer — `get`/`summary`/`array`/
`map`/`xs`/`mosaic`/`animate`), and its transport twins **`model.conc`** (GWT
concentration) and **`model.temp`** (GWE temperature) — the same dependent-variable
grammar, reading the `.ucn` binary via the shared `DependentVariableFile` base. The
readers are **kind-gated**: `.hds` on a transport model (or `.conc` on a flow model)
raises a clear error. `model.outputs.<pkg>.bud` (raw budget accessor),
`model.targets.<family>` (`compare`/`stats`/`calibration_plot`), `model.pest_runs`.

### Deprecated spellings (D12)

The following resolve through `__getattr__` with a `DeprecationWarning` and stay out
of `dir()`/completion; they return exactly what they used to. Prefer the noun form:

| deprecated | use instead |
|---|---|
| `sfr.results.long_profile()` / `plot_long_profile()` | `sfr.results.profile.get()` / `.plot()` |
| `sfr.results.q.plot_profile()` / `sfr.results.stage.plot_profile()` | `…profile.plot()` |
| `lak.results.q.budget_summary()` / `plot_budget()` | `lak.results.q.budget.get()` / `.plot()` |

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
