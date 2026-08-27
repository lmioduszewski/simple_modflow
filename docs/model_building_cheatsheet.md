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
Project                 the workspace everything is declared against
  ├─ inputs/            your reference data (relative paths resolve HERE)
  ├─ grids     ─────────  add_grid("base", GridSpec)     → mf.grid_ref("base")
  ├─ packages  ─────────  add_package("npf/base", spec)  → mf.ref("npf/base")
  ├─ simulations ───────  add_simulation(sim)
  │    └─ SimulationSpec   tdis + solver + model(s)
  │         └─ ModelSpec   = mf.gwf(name, context=ctx, packages=[...])
  │              ├─ ModelContext(grid=, domain=, surfaces=)   ← geometry
  │              └─ packages  disv · npf/ic/sto/oc · chd/ghb/… · uzf/sfr/lak
  └─ runs/<name>/       one workspace per run or scenario
```

**Make the Project first.** It is not just a run wrapper: it is the root that
reference-data paths resolve against, and the registry that `mf.grid_ref` /
`mf.ref` resolve through. Declare against it from the start and the whole model
is a portable recipe — `save()` it, `load()` it elsewhere, rebuild. Build
geometry first and hand the project finished arrays, and you get a project that
runs but cannot be persisted (§7).

Two things that follow, and often trip people:

* **Geometry attaches to the MODEL, not the project.** `ModelContext` rides on
  `mf.gwf(..., context=ctx)`. The project holds the *declarations*; the model
  holds the resolved geometry.
* **The layering does not need the grid to exist.** Surfaces are lazy, and
  `LayerSurfaces` is just an ordered list of them, so §3 comes before §4 on
  purpose. Binding it to a grid (§5) produces a `PackageSpec`, which registers on
  the project like anything else and is referenced with `mf.ref`.

There is no `add_surface`/`add_layers` registry — the layering enters the project
as its **bound DISV package**, which is the same door every other package uses.

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

## 1 · Project — make this first

```python
project = mf.Project("./valley", name="valley")
project.layout.ensure()          # creates specs/, inputs/, runs/ …
```

`mf.Project(root, name=)` touches no disk until asked. The layout:

```
valley/
├─ project.json                  # manifest              project.manifest_path
├─ inputs/                       # YOUR reference data   project.layout.inputs_dir
│   ├─ boundary.gpkg
│   ├─ ground_dem.tif
│   └─ clay_top.gpkg
├─ specs/                        # written by .save()
│   ├─ project_spec.json
│   ├─ simulations/<name>.json
│   ├─ packages/<key>.json
│   └─ grids/<key>.json
├─ grids/                        # resolved grid artifacts (your choice of name)
└─ runs/<name>/                  # MF6 input + output    run.workspace
```

Rename any of those with `ProjectLayout(root, specs_dir_name=…,
simulations_dir_name=…, inputs_dir_name=…, packages_dir_name=…,
grids_dir_name=…)` passed as `layout=`; a custom layout round-trips through
`save()`/`load()`.

### Put reference data under the project and refer to it *relatively*

This is the payoff, and it is worth being deliberate about:

```python
mf.ShapeSource("inputs/boundary.gpkg", crs=CRS)      # ✅ relative to project root
mf.ShapeSource("/home/me/gis/boundary.gpkg")         # ⚠️ absolute — pins the machine
mf.ShapeSource("/mnt/shared/regional.gpkg", external=True)   # deliberate, marked
```

A relative `DataSourceSpec` resolves against `project.root` — **not** the working
directory. Verified 2026-08-26 by building a run after `os.chdir("/")`: the grid
still resolved, `ncpl=293`. `external=True` opts a path out of that, for genuinely
shared data that should not be copied into the project.

### The libraries

Register once, reference from many simulations, change in one place:

```python
project.add_grid("base", grid_spec)                       # → mf.grid_ref("base")
project.add_package("npf/base", mf.npf(k=10.0, icelltype=1))   # → mf.ref("npf/base")
project.add_simulation(sim)
project.add_simulation_from_yaml("specs/simulations/alt.yaml")
```

`project.grids`, `project.packages` and `project.simulation(name)` read them back.

---

## 2 · Contours and rasters → surfaces

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

## 3 · Layering — declare it before the grid

A `Surface` is lazy and knows nothing about any grid, and **`LayerSurfaces` is
just an ordered list of them**. So the layering — which contacts, in what order,
under what names — is declarable with no grid in existence:

```python
layering = mf.LayerSurfaces(
    [ground,                 # [0] is the model TOP
     clay_top,               # each subsequent surface is a layer BOTTOM
     clay_top.below(20),
     bedrock],
    labels=["ground", "upper_sand", "clay", "lower_aquifer"],
)
```

Nothing is read, interpolated or sampled here. `mf.Contours` has not touched
GRASS yet; `mf.Raster` has not opened a file. The grid arrives later, in §5.

`labels` name the columns that land on `vor.gdf_topbtm`, so pick the names you
want to see in section plots and hover rows. There is one more label than there
are layers — the first is the top.

### `LayerStack` defers its grid too — and keeps per-layer control

`LayerSurfaces` is the engine; **`LayerStack` is the facade, and as of 2026-08-27
its `vor` is optional.** So you no longer choose between declaring early and
having per-layer rules — take both:

```python
stack = (
    mf.LayerStack(top=ground, length_units="feet")     # <- no grid
      .add("upper_sand",    bottom=clay_top,           min_thickness=2.0, pinch="inactive")
      .add("clay",          thickness=20.0,            pinch="passthrough")
      .add("lower_aquifer", bottom=bedrock,            min_thickness=5.0, pinch="inactive")
)

layers = stack.build(vor)        # the grid arrives here
print(stack.qc(vor))
disv = stack.to_disv(vor)        # ...or straight to a registerable PackageSpec
```

Bind it once instead, when you want to keep it — or reuse one declaration across
grids, which is the point of returning a copy:

```python
bound = stack.for_grid(vor)
bound.build(); bound.plot.section(y=1000)

coarse = stack.for_grid(vor_coarse).build()    # same layering,
fine   = stack.for_grid(vor_fine).build()      # two grids
```

`LayerStack(vor, top)` is unchanged and still correct whenever the grid already
exists. Using a deferred stack without a grid raises and names the fix, rather
than failing as `NoneType has no attribute ncpl` three frames down; so does
`LayerStack(ground)`, which reads as the deferred form but binds the surface to
`vor`.

| | `LayerSurfaces` | `LayerStack` |
|---|---|---|
| declare before the grid | ✅ | ✅ (since 2026-08-27) |
| `thickness=20.0` sugar | ✗ — derive it (`s.below(20)`) | ✅ |
| `min_thickness` / `pinch` | one **global** rule | **per layer** |
| QC report | `.thickness_report(vor)` | `.qc(vor)` |
| `.plot` namespace | ✗ | ✅ (bind first, or `build(vor).plot`) |
| one call to a DISV spec | `.to_disv(vor)` | `.to_disv(vor)` |

**Prefer `LayerStack`** — it is the documented facade and the per-layer rules are
usually what a real stack needs. Reach for `LayerSurfaces` when you already hold
an explicit surface for every contact and one global pinch rule is right, which
makes the list form the shorter spelling.

---

## 4 · Grid — declare it in the project, then resolve it

Declare the **recipe** on the project, so the grid is part of the saved document
rather than something you rebuilt by hand:

```python
grid_spec = mf.GridSpec.voronoi(
    name="valley",
    boundary=mf.ShapeSource("inputs/boundary.gpkg", crs=CRS),        # relative
    refinement=mf.ShapeSource("inputs/refine_area.gpkg", crs=CRS),   # optional
    crs=CRS,
    boundary_max_area=40_000.0,          # target cell area in CRS units²
)
project.add_grid("base", grid_spec)
```

`boundary=` takes a **`DataSourceSpec`/`ShapeSource`, not a bare path string.**

Then get a real grid object out of it — **still declared in the project**:

```python
vor = project.grids["base"].resolve(
    project_root=project.root,                     # relative paths anchor here
    workspace=project.root / "grids" / "base",     # where Triangle's files land
)
```

**Why resolve now rather than defer?** Because a `LayerStack` and the GIS package
helpers (`mf.uzf/sfr/lak`, `mf.X.gpkg`) resolve cells the moment you call them —
they need a grid object, not a recipe. Resolving from the library gets you both:
the declaration lives in the project and survives `save()`, and you hold a real
grid for the next four sections.

<details><summary>Fully deferred, when the model needs no layer stack or GIS packages</summary>

Hand the model a reference and let the run build the grid:

```python
flow = mf.gwf("valley", grid=mf.grid_ref("base"), packages=[...])
```

The grid is built into `run.workspace/_grid/<model>` at `prepare_run` time.
Verified to round-trip: `save()` → `load()` → `prepare_run` rebuilt the same
293-cell grid. This composes with `disv` + simple BCs only — see gotcha 5.
</details>

<details><summary>Direct construction, if you'd rather drive Triangle yourself</summary>

```python
tri = mf.TriangleGrid(model_ws="grids/base", angle=30)
tri.set_domain_rectangle(x_dist=3000, y_dist=2000, origin=(0, 0), max_area=40_000)
tri.build()
vor = mf.VoronoiGridPlus(tri, crs=CRS)
project.add_grid("base", vor)        # a built grid can go in the library too
```
</details>

---

## 5 · Bind the layering to the grid, and register it

`to_disv` samples the surfaces onto the cells and hands back a **`PackageSpec`** —
so the layer stack becomes an ordinary project package:

```python
disv = layering.to_disv(
    vor,
    pinch_out=True, minimum_thickness=2.0,     # thin cells -> idomain
    length_units="FEET",
    attach=True,                               # publish onto vor.gdf_topbtm
)
project.add_package("disv/base", disv)         # <- registered on the project
```

and the model refers to it by name, alongside every other library entry:

```python
flow = mf.gwf("valley", context=ctx, packages=[
    mf.ref("disv/base"),          # the layer stack
    mf.ref("npf/base"),
    mf.ic(strt=150.0),
    ...
])
```

Verified end to end 2026-08-26: layering declared before the grid, grid resolved
from the project library, `to_disv` → `add_package` → `mf.ref`, run converged at
`nlay=3`, and the whole project `save()`/`load()`-ed and rebuilt.

**Check it before you trust it.** `thickness_report` is the `LayerSurfaces`
equivalent of `qc()` and takes the grid, since that is what thickness needs:

```python
print(layering.thickness_report(vor, minimum_thickness=2.0))
```

### Or bind through the `LayerStack` facade

The facade gives `.add(thickness=…)`, per-layer `min_thickness`/`pinch`, `.qc()`
and `.plot`. It can be declared back in §3 with no grid and bound here, or
constructed with the grid directly when the grid already exists:

Sources mix here exactly as they do in `LayerSurfaces` (§2) — `top=` and
`bottom=` take a raster, a contoured contact or a derived one interchangeably:

```python
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

`stack.build()` returns arrays; feed them to `mf.disv(top=layers.top, …)` as in §7
— or use `LayerSurfaces.to_disv` above, which does it in one call and is the
registerable form.

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

#### Where a layer only exists in part of the domain

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

#### What `build()` returns

`LayerBuildResult` — fields `top`, `botm`, `idomain`, `thickness`, `names`,
`nlay`, `vor`; methods `qc()`, `report()`, `validate()`, `prune_isolated()`,
`attach_to_grid()`, and the `plot` namespace.

#### `pinch=` — what happens where a layer goes thin

| value | effect |
|---|---|
| `"passthrough"` *(default)* | `idomain = -1` — cell inactive but vertically transmissive |
| `"inactive"` | `idomain = 0` — a true pinch-out, no flow through |
| `"floor"` | clamp thickness to `min_thickness`, cell stays active |

#### Erosion — a surface that cuts down through the layers below

An incised channel, a buried valley, a scoured contact: the overlying surface
drops *below* one or more contacts underneath it. **There is no separate "cut"
operation** — `reconcile="bottom"` already does it. It pushes every contact the
cut passes through down to `cut - min_sep`, cascading, and the layers that end up
thinner than `min_thickness` are handled by their `pinch`.

Measured on a 441-cell grid, ground at 200 with a channel incised to 115, through
`sand` (base 150) and `clay` (base 120) into `till`:

```python
stack = (mf.LayerStack(top=ground, length_units="feet")
         .add("sand", bottom=Flat(150), min_thickness=2.0, pinch="inactive")
         .add("clay", bottom=Flat(120), min_thickness=2.0, pinch="inactive")
         .add("till", bottom=Flat(60),  min_thickness=2.0, pinch="inactive"))
```

```
              thickness in channel      idomain in channel
  sand              0.10                  0    ← cut out
  clay              0.10                  0    ← cut out
  till             54.80                  1    ← now directly under the channel floor
```

with `pinch="passthrough"` the cut-out layers read `idomain = -1` instead — same
geometry, but flow still passes vertically through the gap. **That is usually the
one you want for an erosional cut**, since a column of `0`s can disconnect the
till from the channel above it.

`qc()` names the cut rather than hiding it:

```
[0] sand   nan_bottom=0 thin=139 pinched=139 ... reconcile_moved=139 (max 10.10)
```

##### Or declare the cut explicitly

Relying on reconcile means the geometry is a side effect of a repair pass. To
make it intentional, cap the contact against the cutting surface:

```python
.add("sand", bottom=Flat(150).capped_at(ground - 2.0), min_thickness=2.0)
```

Same channel, different answer — and the difference is hydrogeology, not style:

| | in the channel | `qc()` |
|---|---|---|
| reconcile does it | sand `0.10` thick, **idomain 0** — unit absent | `reconcile_moved=139` |
| `capped_at(ground - 2)` | sand `2.00` thick, **idomain 1** — thin veneer kept | `reconcile_moved=0` |

Use reconcile when the unit is genuinely truncated and gone; use `capped_at` when
a remnant survives, or when you want the cell to stay active. The `reconcile_moved=0`
in the second row is the tell that you described the geometry rather than
discovering it.

#### Look at it before you trust it

```python
layers.plot.section(y=1000)          # filled, layer-coloured cross-section
layers.plot.map("clay")              # per-layer thickness
layers.plot.surface("all")           # every contact as 3-D height fields
layers.plot.grid(scale=12)           # the layered cell VOLUME (needs `viz3d` extra)
```

Starting from an existing model instead? `mf.LayerStack.from_modflow(vor, source)`.

---


### Why registering matters for persistence

Putting the DISV in the **package library** rather than inline in the model is
what lets an array-heavy model be saved at all. A simulation must be JSON-clean,
and `to_disv` emits numpy arrays for `top`/`botm` — but a *library* entry falls
back to a pickle sidecar beside its JSON:

```
specs/packages/disv/base.json     +  base.pkl  +  base.versions.json
```

So the simulation holds only `mf.ref("disv/base")`, stays serializable, and the
arrays ride along in the sidecar. Verified: `save()` → `load()` → `prepare_run`
rebuilt the same 3-layer model. §8 has the full persistence rules.

---

## 6 · Context

```python
ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)
```

Attaches to the **model** via `mf.gwf(..., context=ctx)`, and the built model
exposes it as `model.myflopy_context`. Every GIS-aware helper takes `context=`.

---

## 7 · Packages

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

## 8 · Simulation, runs, scenarios, persistence

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

### Persistence is opt-in, and the library is the escape hatch

`project.save()` writes the project + its simulations + its package library;
`mf.Project.load(root)` reads them back.

The rule that matters: **a SIMULATION must be JSON-clean, and numpy is not.**

```python
mf.npf(k=[1.0] * ncpl)                      # list     -> serializable
mf.npf(k=np.ones(ncpl))                     # ndarray  -> blocks save()
mf.oc(saverecord=[("HEAD", np.int64(1))])   # even a numpy SCALAR blocks it
```

So an array-laden package written **inline** in the model — `mf.disv(top=layers
.top, ncpl=gp["ncpl"], …)` — cannot be saved. Three ways out, best first:

1. **Register it and reference it (§5).** `project.add_package("disv/base", disv)`
   + `mf.ref("disv/base")`. The simulation then holds only a reference and stays
   JSON-clean, while the library entry falls back to a **pickle sidecar**:
   `specs/packages/disv/base.json` + `base.pkl` + `base.versions.json`. This is
   the intended shape, and it is why declaring against the project pays off.
   Verified to round-trip and rebuild.
2. **Convert at the boundary** — `top=layers.top.tolist()`, `ncpl=int(gp["ncpl"])`.
   Also verified, but you must remember it at every call site.
3. **Don't save.** Runs, scenarios, discovery and reopening never needed it.

Call **`project.validate()`** to see which case you are in; it lists the offending
simulations without touching disk, and `save()` raises the same message rather
than writing a project that cannot be reloaded. Note it inspects *simulations*, so
a model built entirely from `mf.ref`/`mf.grid_ref` reports clean — correctly, since
the arrays are the library's problem and the library can hold them.

For a single simulation the lighter option is `sim.to_yaml(path)` /
`mf.SimulationSpec.from_yaml(path)` (`Project.add_simulation_from_yaml` too),
which carries the same serializability rule.

---

## 9 · Read and draw

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
