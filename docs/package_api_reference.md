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

`mf.disv`, `mf.dis`, `mf.disu`, `mf.ic`, `mf.npf`, `mf.sto`, `mf.oc`, `mf.tdis`, `mf.ims`
— each `mf.<name>(...)` returns a `PackageSpec`.

`mf.ic` / `mf.oc` / `mf.disv` / `mf.dis` / `mf.disu` are **model-kind-aware**: they
build the GWF/GWT/GWE FloPy class off the model's type, so the same helper serves a
transport or energy model (for GWT/GWE `mf.oc`, pass
`concentration_filerecord=`/`temperature_filerecord=` via `**options`). `mf.npf`
and `mf.sto` are GWF-only (MF6 has no GWT/GWE variant).

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
grids (use FloPy's own plotting, or DISV). PRT is not dispatched by any grid helper —
`mf.prt` builds its own dis/disv internally (and MF6 has no `ModflowPrtdisu`).

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
hover** and the **results tier** (`model.conc`/`model.temp`) are Phase 6.

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

**Model-level reads:** `model.hds` (heads explorer — `get`/`array`/`map`/`xs`/
`mosaic`/`animate`), `model.outputs.<pkg>.bud` (raw budget accessor),
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
