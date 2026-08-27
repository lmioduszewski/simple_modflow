# Building a model, package-first — cheat sheet

Every call below was executed end to end on 2026-08-26, in a shell with **no
GRASS environment variables set**: DEM raster + digitized contours + a second
raster → 3-layer stack → DISV → MF6 converged → heads read → maps drawn. The
runnable version is **`examples/mf6/contours_to_model.py`**; it writes its own
raster and contour fixtures, so it runs anywhere.

Source kinds mix freely — §3 is the case where some contacts are rasters and some
are contours. §6 covers project workspace layout, reusable libraries, scenarios,
and what persistence will and won't accept.

> Reference for the read-side grammar: `docs/package_api_reference.md`.
> Normative view-layer rules: `docs/view_layer_conventions.md`.

---

## The mental model

```
Project                 durable workspace + run lifecycle  — holds NO geometry
  └─ SimulationSpec     one MF6 simulation: tdis + solver + model(s)
       └─ ModelSpec     = mf.gwf(name, context=ctx, packages=[...])
            ├─ ModelContext(grid=, domain=, surfaces=)   ← geometry, on the MODEL
            └─ packages  disv · npf/ic/sto/oc · chd/ghb/drn/riv/wel/rch/evt
                         · uzf/sfr/lak · mvr
```

Geometry rides on the **model**, not the project. The project is the lifecycle
wrapper — runs, scenarios, reuse.

---

## 0 · Setup

```python
import numpy as np
import myflopy as mf
```

**Using `mf.Contours`?** It shells out to GRASS, which is a **system**
dependency — install it with your package manager (Linux/macOS) or via OSGeo4W
(Windows). Nothing else is needed: the launcher is found on `PATH` (or by
globbing the OSGeo4W/QGIS bundles on Windows), and the GRASS Python bindings —
which live inside the install at `<prefix>/etc/python` and are *not* on
`sys.path` — are located by asking that launcher, so no `PYTHONPATH` either.

**It is pure Python from your side.** `mf.Contours` starts and tears down its own
GRASS session per interpolation — you never open a GRASS shell, never set
`GRASS_BIN`, never set `PYTHONPATH`, and nothing is left running afterwards.
Re-verified 2026-08-26 with both variables explicitly unset:

```
launcher : /usr/bin/grass                 (shutil.which)
bindings : /usr/lib/grass84/etc/python    (grass --config python_path)
```

Only if GRASS is installed somewhere undiscoverable, point at the launcher:

```bash
export GRASS_BIN=/usr/bin/grass     # or .../grass84.bat on Windows
```

Skip all of this if you're starting from rasters.

---

## 1 · Grid

Build it **eagerly**. The GIS-aware helpers (`mf.uzf/sfr/lak`, `mf.X.gpkg`)
resolve cells the moment you call them, so they need a real grid, not a recipe.

```python
grid = mf.GridSpec.voronoi(
    name="valley",
    boundary=mf.ShapeSource("boundary.gpkg", crs="EPSG:2927"),
    refinement=mf.ShapeSource("refine_area.gpkg", crs="EPSG:2927"),  # optional
    crs="EPSG:2927",
    boundary_max_area=40_000.0,      # target cell area in CRS units²
)
vor = grid.resolve(workspace="_grid")     # -> VoronoiGridPlus
```

`boundary=` takes a **`DataSourceSpec`/`ShapeSource`, not a bare path string.**

<details><summary>Direct construction, if you'd rather drive Triangle yourself</summary>

```python
tri = mf.TriangleGrid(model_ws="_grid", angle=30)
tri.set_domain_rectangle(x_dist=3000, y_dist=2000, origin=(0, 0), max_area=40_000)
tri.build()
vor = mf.VoronoiGridPlus(tri, crs="EPSG:2927")
```
</details>

---

## 2 · Contours → surfaces

This is the atom layer. A `Surface` is lazy — it interpolates **once**, caches the
raster next to the source, then samples many times.

```python
ground   = mf.Contours("ground_contours.gpkg", z="Elev", epsg="2927",
                       resolution=25, region_vector="boundary.gpkg")
clay_top = mf.Contours("clay_top.gpkg", z="Elev", epsg="2927",
                       resolution=25, region_vector="boundary.gpkg")
```

**A region is required** — pass `region_vector=` (any polygon; your boundary is
the natural choice) or `region_raster=`. Without one GRASS has no extent and
raises. `z=` names the elevation attribute on the contour lines; it defaults to
`"Elev"`.

The interpolated GeoTIFF lands beside the contours as `<name>.interp.tif` and is
reused on later runs. Pass `out=` to place it elsewhere.

### Every surface source

| source | use |
|---|---|
| `mf.Contours(path, z="Elev", region_vector=...)` | digitized elevation contours (GRASS) |
| `mf.Raster(path, fill=...)` | a DEM or any GeoTIFF |
| `mf.Points(xs, ys, zs, method="linear")` | scattered control points — **the GRASS-free fallback** |
| `Array(values)` | one value per cell, already computed (`from myflopy.layers import Array`) |
| `mf.Surface.from_array(values)` | same, via the atom class |

### Surface algebra

Surfaces compose, so a contact you have no contours for can be *derived*:

```python
clay_base = clay_top.below(25)                       # 25 ft below, parallel
bedrock   = ground.below(120)
capped    = clay_top.capped_at(ground - 2)           # never within 2 ft of ground
floored   = clay_top.floored_at(bedrock + 5)
blended   = mf.Surface.maximum(clay_top, bedrock)    # elementwise, BOTH surfaces
```

**Watch the binding.** `maximum` / `minimum` / `shift` / `clamp` / `isopach` /
`where` are **classmethods**, so calling one on an instance binds the instance to
`cls` and drops it. Every such spelling now raises `TypeError` naming the fix,
so this is a stop, not a silent wrong answer:

```python
mf.Surface.maximum(a, b)     # correct — max of a and b
a.floored_at(b)              # correct — the same thing, instance-bound
a.maximum(b)                 # TypeError: takes two or more surfaces, got 1
```

The instance-bound pair `capped_at` / `floored_at` do the two-surface job safely
(`a.capped_at(b)` is `min(a, b)`), so prefer them. Reach for `mf.Surface.minimum`
/ `maximum` when you have three or more surfaces to blend.

Guards, and why each one is where it is: `maximum`/`minimum` take `*surfaces`, so
a misbound call arrives as one operand rather than an arity error — they reject
fewer than two. `clamp`'s bounds are both optional, so a misbound `a.clamp(b)`
arrives with neither — it requires `lower=` or `upper=`. `shift`, `isopach` and
`where` already raise for the now-missing positional argument, so they need no
guard of their own.

| instance methods (safe to chain) | classmethods (call on `mf.Surface`) |
|---|---|
| `.above(d)` `.below(d)` `.between(lower=, upper=)` | `mf.Surface.shift(s, d)` |
| `.capped_at(o)` `.floored_at(o)` | `mf.Surface.maximum(a, b, ...)` `mf.Surface.minimum(a, b, ...)` — two or more |
| `.thickness()` `.within(zone, outside=)` | `mf.Surface.clamp(s, lower=, upper=)` |
| `+` `-` | `mf.Surface.isopach(src)` `mf.Surface.where(zone, inside, outside)` |

---

## 3 · Layer stack

Declare **top once**, then each layer by its *bottom* or its *thickness*.

### Mixing rasters and contours

This is the normal case, and it needs no special handling: **`Surface` is one
type whatever produced it**, so a raster contact, a contoured contact and a
derived contact are interchangeable arguments to `top=` and `bottom=`. Nothing
below cares where a surface came from.

```python
ground   = mf.Raster("ground_dem.tif")                       # LiDAR DEM
clay_top = mf.Contours("clay_top.gpkg", z="Elev",            # digitized contours
                       region_vector="boundary.gpkg")
bedrock  = mf.Raster("bedrock.tif")                          # another raster

stack = (
    mf.LayerStack(vor, top=ground, length_units="feet")
      .add("upper_sand",    bottom=clay_top,           min_thickness=2.0, pinch="inactive")
      .add("clay",          bottom=clay_top.below(20), pinch="passthrough")
      .add("lower_aquifer", bottom=bedrock,            min_thickness=5.0, pinch="inactive")
)

print(stack.qc())          # read this BEFORE build — it is a report, not a picture
layers = stack.build(attach=True)
```

Verified 2026-08-26 on a 293-cell grid, three layers from three different source
kinds:

```
[0] upper_sand      thickness   28.89 ..   61.83   active 293/293
[1] clay            thickness   20.00 ..   20.00   active 293/293
[2] lower_aquifer   thickness   51.31 ..   71.69   active 293/293
```

Three things worth knowing when the sources are mixed:

* **Resolution differences don't matter.** Every surface is sampled onto the same
  Voronoi cells by area-weighted averaging (`grid/surfaces.py`), so a 25-ft DEM
  and a coarsely-contoured contact land on the same footing. You do **not** need
  to pre-align, resample, or clip them to each other.
* **Only `mf.Contours` costs anything.** Its GRASS interpolation runs once, caches
  to `<name>.interp.tif`, and is reused on every later build. Rasters are read
  directly. So a stack with one contoured contact is slow on the *first* build
  and fast thereafter.
* **Mismatched CRS is the one thing to check yourself.** Sources are sampled in
  the grid's CRS; a contact in a different projection samples in the wrong place
  and shows up as a nonsense `qc()` thickness rather than an error.

A contact you have *neither* a raster nor contours for is derived (`clay_top
.below(20)` above) — see the algebra table in §2.

### Where a layer only exists in part of the domain

Give it a source that goes above the layer above it, and let `pinch="inactive"`
cut it out — that is the whole mechanism. `qc()` reports the count, and
`layers.idomain[k]` is what MODFLOW sees.

`qc()` reports per layer: NaN bottoms, cells below `min_thickness`, pinch-outs,
non-positive thickness in active cells, and how far the top-down reconcile had to
move anything. A clean stack looks like:

```
LayerStack QC [OK]: 3 layers, 293 cells, 1 active component(s)
  [0] upper_sand     nan_bottom=0 thin=0 pinched=0 thickness<=0(active)=0  reconcile_moved=0
```

**`attach=True` also publishes elevations onto `vor.gdf_topbtm`**, which the
grid-aware builders need — SFR reach tops and LAK lake-cell layering read it. Set
it whenever you use `mf.sfr`/`mf.lak`. It is also what gives maps their
layer-elevation hover rows.

### What `build()` returns

`LayerBuildResult` — fields `top`, `botm`, `idomain`, `thickness`, `names`,
`nlay`, `vor`; methods `qc()`, `report()`, `validate()`, `prune_isolated()`,
`attach_to_grid()`, and the `plot` namespace.

### `pinch=` — what happens where a layer goes thin

| value | effect |
|---|---|
| `"passthrough"` *(default)* | keep the cell active; layer stays vanishingly thin |
| `"inactive"` | `idomain = 0` — a true pinch-out |
| `"floor"` | clamp thickness to `min_thickness` |

### Look at it before you trust it

```python
layers.plot.section(y=1000)          # filled, layer-coloured cross-section
layers.plot.map("clay")              # per-layer thickness
layers.plot.surface("all")           # every contact as 3-D height fields
layers.plot.grid(scale=12)           # the layered cell VOLUME (needs `viz3d` extra)
```

Starting from an existing model instead? `mf.LayerStack.from_modflow(vor, source)`.

---

## 4 · Context

```python
ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)
```

Attaches to the **model** via `mf.gwf(..., context=ctx)`, and the built model
exposes it as `model.myflopy_context`. Every GIS-aware helper takes `context=`.

---

## 5 · Packages

One flat, declarative list — the readable payoff.

```python
gp = vor.get_disv_gridprops()

flow = mf.gwf(
    "valley",
    context=ctx,
    save_flows=True,
    newtonoptions="UNDER_RELAXATION",       # for drying/rewetting
    packages=[
        mf.disv(nlay=layers.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                vertices=gp["vertices"], cell2d=gp["cell2d"],
                top=layers.top, botm=layers.botm, idomain=layers.idomain,
                length_units="FEET"),
        mf.ic(strt=150.0),
        mf.npf(k=[25.0, 0.05, 40.0], k33=..., icelltype=1, save_flows=True),
        mf.sto(steady_state={0: True}),

        mf.ghb.gpkg("bcs.gpkg", layer="underflow", context=ctx, nper=nper),
        mf.drn(stress_period_data=drn_data),
        mf.rch(context=ctx, nper=nper, recharge=rch_array),

        mf.uzf(context=ctx, nper=nper, cells=uzf_cells, vks=0.25,
               thtr=0.08, thts=0.34, thti=0.17, eps=4.0, finf=finf, pet=pet),
        sfr, lak,                            # see below
        mf.oc(head_filerecord="valley.hds", budget_filerecord="valley.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ],
)
```

`k` accepts a scalar, one value per layer, or a full `(nlay, ncpl)` array.

### Three spellings for every list BC

`chd` `ghb` `drn` `riv` `wel` `rch` `evt` each answer all three:

| spelling | when |
|---|---|
| `mf.ghb(stress_period_data=...)` | you already have the records |
| `mf.ghb.gpkg(path, context=ctx, nper=n, head="head", conductance="conductance")` | driven from GIS |
| `mf.ghb.flopy(...)` | raw FloPy escape hatch |

### Surface water

```python
sfr = mf.sfr(context=ctx, nper=nper, streams="streams.gpkg",
             connection_mode="automatic", width=18.0, gradient=0.0012,
             roughness=0.030, streambed_k=0.05, streambed_thickness=1.5,
             inflow={0: {"north_trib": 15000.0}}, mover=True)

lak = mf.lak(context=ctx, nper=nper, lakes="lakes.gpkg", lake_id_field="name",
             starting_stage={"valley_lake": 101.0}, bed_leakance=0.11,
             connection_modes="automatic", mover=True)

mvr = mf.mvr(nper=nper, moves=(
    mf.Move(mf.sfr_connection(sfr, "main_stem"),
            mf.lak_connection(lak, "valley_lake"), value=0.5),
))
```

Pull `sfr`/`lak` out as handles so the mover references them **semantically**
rather than by raw reach index. MVR is validated: moved packages must be declared
in the model *and* ordered before the mover.

---

## 6 · Simulation, project, run

```python
sim = mf.SimulationSpec(
    "baseline",
    models=(flow,),
    packages=(
        mf.tdis(nper=nper, perioddata=[(1.0, 1, 1.0)] * nper),
        mf.ims(models=("valley",), complexity="COMPLEX",
               outer_maximum=250, inner_maximum=250,
               linear_acceleration="BICGSTAB"),
    ),
)

project = mf.Project("./project", name="valley")
project.add_simulation(sim)

run = project.prepare_run("baseline", "baseline")   # builds in memory
success, report = run.execute()                     # writes + runs MF6
if not success:
    print("\n".join(report[-25:]))
```

### Setting up the workspace

`mf.Project(root, name=)` creates nothing on disk until you ask it to. Runs land
under `<root>/runs/<name>`; persistence (opt-in, see below) writes under
`<root>/specs/`:

```
project/
├─ project.json                     # the manifest        (project.manifest_path)
├─ specs/
│   ├─ project_spec.json            # written by .save()
│   ├─ simulations/<name>.json
│   ├─ packages/<key>.json          # the reusable library
│   └─ grids/<key>.…
└─ runs/
    ├─ baseline/                    # MF6 input + output  (run.workspace)
    └─ dry/
```

Rename any of those directories with `ProjectLayout(root, specs_dir_name=...,
simulations_dir_name=..., inputs_dir_name=..., packages_dir_name=...,
grids_dir_name=...)` passed as `layout=`. A custom layout round-trips through
`save()`/`load()`.

### Reusable libraries

Register a grid or a package once and reference it from several simulations, so a
change lands in one place:

```python
project.add_grid("valley_voronoi", vor)
project.add_package("npf/base", mf.npf(k=[25.0, 0.05, 40.0], icelltype=1))

flow = mf.gwf("valley", context=ctx, packages=[
    mf.disv(...),
    mf.ref("npf/base"),          # ← resolved from the library when the run builds
    ...
])
```

### Scenarios

A scenario is just another simulation on the same project. The runs are
independent workspaces, so you can compare them directly:

```python
project.add_simulation(make_sim("baseline", rch=4.0e-4))
project.add_simulation(make_sim("dry",      rch=1.0e-4))

base = project.prepare_run("baseline", "baseline"); base.execute()
dry  = project.prepare_run("dry",      "dry");      dry.execute()

drawdown = base.model("valley").hds.array(layer=0) - dry.model("valley").hds.array(layer=0)
```

| call | does |
|---|---|
| `project.prepare_run(name, sim)` | build **in memory** — nothing written yet |
| `run.execute()` | write MF6 input, run it, return `(success, report)` |
| `project.run(name, sim)` | both at once |
| `project.discover_runs()` | every run with a lifecycle manifest |
| `project.reopen_run(name)` | reopen one by name |
| `mf.load_run(path)` | reopen a run from a bare path, no project needed |

**`prepare_run` before `execute` is the useful habit** — the model is fully built
and inspectable (`run.model("valley")`, and every plotting verb) while no input
files exist yet, so a bad layer stack or a mis-mapped BC is visible before MF6
ever sees it.

### Persistence is opt-in, and refuses numpy

`project.save()` writes the project + its simulations + its package library;
`mf.Project.load(root)` reads them back. But a spec can only be persisted if it is
**serializable**, and *numpy is not*:

```python
mf.npf(k=[1.0] * ncpl)          # list      -> serializable
mf.npf(k=np.ones(ncpl))         # ndarray   -> NOT serializable
mf.oc(saverecord=[("HEAD", np.int64(1))])   # even a numpy SCALAR blocks it
```

Which means an eagerly-built stack — whose `layers.top/botm/idomain` are arrays,
and whose `vor.get_disv_gridprops()["ncpl"]` is a `numpy.int64` — cannot be saved
as-is. Two honest options:

* **Don't save.** Runs, scenarios, discovery and reopening all work without it.
  This is the normal choice for an array-heavy model.
* **Convert at the boundary** — `top=layers.top.tolist()`, `ncpl=int(gp["ncpl"])`,
  and so on. Verified to round-trip through `save()` → `load()`.

Call **`project.validate()`** to find out which it is; it lists the offending
simulations without touching disk, and `save()` raises with the same message
rather than writing a project that cannot be reloaded.

For a single simulation the lighter option is `sim.to_yaml(path)` /
`mf.SimulationSpec.from_yaml(path)` (`Project.add_simulation_from_yaml` too),
which carries the same serializability rule.

---

## 7 · Read and draw

```python
model = run.model("valley")

model.hds.array(layer=0)                 # numpy
model.packages.ghb.results.q.get()       # a DataFrame

model.plot.map(layer=0, contours=True)   # six verbs, four scopes
model.plot.section(cells=[12, 40, 88])
model.plot.grid()
model.plot.animate([model.plot.map(layer=0, per=p) for p in range(model.nper)])
```

Read grammar is always `model.packages.<pkg>.<inputs|results>.<noun>.<verb>()`.
Everything drawn returns a `Picture`: `.fig` / `.show()` / `.save(p)` / `.html(p)`,
never a trailing `.plot()`.

---

## Gotchas that cost real time

1. **`mf.Contours` needs GRASS *installed*** — it is a system dependency, not a
   pip one, so `pip install` will never supply it. Discovery itself is automatic
   on every platform (fixed 2026-08-26; it used to glob `grass*.bat` only and so
   never found `/usr/bin/grass`). `GRASS_BIN` remains the override for an install
   that isn't on `PATH`.
2. **`mf.Contours` needs a region** — `region_vector=` or `region_raster=`.
3. **`boundary=` wants a `ShapeSource`,** not a path string.
4. **`GridSpec.resolve` is keyword-only**: `resolve(workspace=...)`, not
   `resolve(path)`.
5. **Build the grid eagerly.** Deferred `GridSpec` composes with disv + simple BCs
   only — `mf.uzf/sfr/lak` and `.gpkg` helpers resolve cells at declaration.
6. **`attach=True`** on `stack.build()` whenever SFR or LAK is in the model.
7. **Read `stack.qc()` before `build()`.** A silently pinched-out layer becomes an
   inactive region you'll spend an afternoon chasing in the head field.
8. `GridSpec.structured` / `from_geopackage` exist but **fail fast** — only
   `voronoi` and `python` resolve. Use `GridSpec.from_object` for a built grid.
9. **`project.save()` refuses numpy** — arrays *and* scalars. An eagerly-built
   stack cannot be persisted without `.tolist()` at the boundary. Call
   `project.validate()` to see this without touching disk; running and reopening
   work fine either way.
10. **Mixed source resolutions need no pre-alignment.** Area-weighted sampling
    puts a 25-ft DEM and a coarse contoured contact on the same footing. A
    mismatched **CRS**, though, samples in the wrong place and surfaces as a
    nonsense `qc()` thickness rather than an error — check that yourself.
