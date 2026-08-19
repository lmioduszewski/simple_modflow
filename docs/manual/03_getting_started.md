# 3. Getting Started

This chapter is a **complete, end-to-end overview** of building, running,
visualizing, and calibrating a model with myflopy. It is deliberately dense — a
cheat sheet you can return to — but every concept here is expanded in a later
chapter, and cross-references point the way.

If you read nothing else, read this. By the end you will have:

- the **mental model** of how the pieces fit together;
- a **minimal model** that builds and runs (verified — it really runs);
- a **realistic, multi-package model** assembled the preferred way;
- and a tour of **results, visualization, and calibration**.

> **Prerequisites.** myflopy importable (`import myflopy as mf`), the `mf6`
> executable on your `PATH`, and — for the calibration section only — `pestpp-ies`
> and `pyemu`. See [Chapter 2, *Installation*](02_installation.md).

---

## 3.1 The mental model in one diagram

myflopy describes a model as a set of **immutable specifications**, then a
**`Project`** materializes them into FloPy objects, writes MODFLOW 6 input, and
runs it.

```
Project                  durable workspace + run/scenario lifecycle (holds NO geometry)
  └─ SimulationSpec      one MF6 simulation: tdis + ims solver + the model(s)
       └─ ModelSpec      one model = mf.gwf(name, context=ctx, packages=[...])
            ├─ ModelContext   grid + domain(idomain) + surfaces + dates  ← rides on the MODEL
            └─ packages       mf.disv, mf.ic/npf/sto/oc,
                              mf.chd/ghb/drn/riv/wel/rch/evt,
                              mf.uzf/sfr/lak, mf.mvr
```

Five terms unlock everything else:

| Term | What it is |
|---|---|
| **`Project`** | A durable workspace that owns runs and scenarios. It holds **no geometry**. |
| **`SimulationSpec`** | One MF6 simulation: time discretization (`tdis`), a solver (`ims`), and one or more models. |
| **`ModelSpec`** | One GWF / GWT / GWE / PRT model, built with `mf.gwf(...)` etc. and a flat list of package specs. |
| **`ModelContext`** | The **geometry** — grid, active domain (`idomain`), surfaces, dates. It attaches to the *model*, and the built model exposes it as `model.myflopy_context`. |
| **package spec** | A serializable description of one MODFLOW package, produced by `mf.disv(...)`, `mf.ghb(...)`, `mf.sfr(...)`, … |

> **Note — the package-first API is the canonical way to build models.**
> `mf.gwf`, `mf.disv`, `mf.ghb`, `mf.sfr`, … are the preferred surface. Behind each
> high-level helper sits an *engine* (`mf.uzf` runs `UZFBuilder`, `LayerStack`
> compiles to `LayerSurfaces`), but you rarely touch the engine directly. See
> [Chapter 4](04_architecture.md).

---

## 3.2 Hello, model — the smallest thing that runs

A one-layer, steady-state model on a small Voronoi grid, with a constant-head
gradient across it. This is complete and **verified to build, run, and converge**.

```python
import numpy as np
import myflopy as mf
from myflopy.layers import Array
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi

# 1) GRID — a small Voronoi grid. (Real models build this from a boundary
#    polygon; see §3.3 and Chapter 6. Here we use a demo grid generator.)
vor = rectangular_voronoi(CanonicalModelConfig(nrow=12, ncol=12, nlay=1, nper=1))
gp = vor.get_disv_gridprops()
ncpl = int(vor.ncpl)

# 2) LAYERS — one 20-ft-thick aquifer under a flat ground surface at 100 ft.
layers = (
    mf.LayerStack(vor, top=Array(np.full(ncpl, 100.0)), length_units="feet")
    .add("aquifer", thickness=20.0)
    .build(attach=True)                 # -> disv-ready top / botm / idomain
)

# 3) CONTEXT — the geometry the package builders read; rides on the model.
ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)

# 4) A constant-head gradient: 95 ft on the first cell, 90 ft on the last.
chd = {0: [[(0, 0), 95.0], [(0, ncpl - 1), 90.0]]}

# 5) MODEL — one flat, declarative package list.
flow = mf.gwf("hello", context=ctx, packages=[
    mf.disv(nlay=1, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
            vertices=gp["vertices"], cell2d=gp["cell2d"],
            top=layers.top, botm=layers.botm, idomain=layers.idomain,
            length_units="FEET"),
    mf.ic(strt=92.0),
    mf.npf(k=10.0),
    mf.sto(steady_state={0: True}),
    mf.chd(stress_period_data=chd),
    mf.oc(head_filerecord="hello.hds", budget_filerecord="hello.cbc",
          saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
])

# 6) SIMULATION — time discretization + solver + the model.
sim = mf.SimulationSpec("baseline", models=(flow,), packages=(
    mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
    mf.ims(models=("hello",)),
))

# 7) RUN — a Project materializes, writes, and runs it.
project = mf.Project("runs/hello")
project.add_simulation(sim)
run = project.run("baseline", "baseline")          # build + write + run MF6

print("converged:", run.success)                    # -> True
heads = run.model("hello").hds.array(layer=0)
print("head range: %.2f -> %.2f" % (heads.min(), heads.max()))   # -> 90.00 .. 95.00
```

That is the whole shape of a myflopy program: **grid → layers → context →
packages → simulation → run → results.** Everything below adds detail to those
seven steps.

> **Tip.** `project.run(name, sim)` is the one-shot "build and execute" call. When
> you want to inspect the in-memory model *before* running, use
> `run = project.prepare_run(name, sim)` then `run.execute()` — see §3.5.

---

## 3.3 A real model, end-to-end

Real models have **layers from rasters/contours**, **boundary conditions from
GeoPackages**, and a **surface-water network** (SFR + LAK + UZF) tied together
with a **mover**. The structure is identical to §3.2 — just richer. The example
below is condensed from the tested reference template at
[`examples/mf6/package_first_full_stack.py`](../../examples/mf6/package_first_full_stack.py),
which builds a two-layer valley with **GHB, RCH, DRN, UZF, SFR, LAK, and MVR** and
converges.

### Step 1 — a Project

```python
import myflopy as mf

project = mf.Project("modeling/valley", name="valley_demo")
```

The `Project` is the durable home for runs and scenarios. It holds no geometry.

### Step 2 — a grid (eager, from GIS)

The GIS-aware package helpers (`mf.uzf`, `mf.sfr`, `mf.lak`, and the `.gpkg`
forms) resolve grid cells the *moment you call them*, so the grid must exist
first. Build it **eagerly** from a boundary polygon (and optional refinement):

```python
grid_spec = mf.GridSpec.voronoi(
    boundary=mf.ShapeSource("gis/domain.shp"),
    refinement=mf.ShapeSource("gis/streams.shp"),   # finer cells near streams
)
vor = grid_spec.resolve(workspace="modeling/valley/_grid")   # -> VoronoiGridPlus
gp = vor.get_disv_gridprops()
```

> **Warning — eager vs. deferred.** You *can* defer grid construction to run time
> with `mf.gwf(...).with_grid(grid_spec)`, but deferred grids do **not** yet
> compose with the GIS package helpers (`mf.uzf/sfr/lak/.gpkg`). For the full
> surface-water stack, build the grid eagerly as above. See
> [Chapter 6, §6.5](06_grids.md).

### Step 3 — layers from surfaces

Author layers with the **`LayerStack`** facade: a top surface, then named layers
added by thickness or bottom surface. Surfaces come from rasters, contours,
points, or arrays.

```python
from myflopy.layers import Raster, Contours

stack = (
    mf.LayerStack(vor, top=Raster("gis/ground.tif"), length_units="feet")
    .add("alluvium", thickness=40.0)                       # upper aquifer
    .add("basin_fill", bottom=Contours("gis/bedrock.shp")) # to a bedrock surface
)
layers = stack.build(attach=True)     # top/botm/idomain; publishes vor.gdf_topbtm
print(layers.report())                # per-layer thickness + pinch summary
print(layers.qc())                    # problems to fix before MF6 runs
```

`build(attach=True)` also writes the surface elevations onto `vor.gdf_topbtm`, so
the grid-aware builders (SFR reach tops, LAK lake-cell layering) can read them.
See [Chapter 7, *Layers & Surfaces*](07_layers.md).

### Step 4 — the model context

```python
ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)
```

This is the geometry every builder reads. It rides on the model.

### Step 5 — boundary conditions (three forms)

Every list boundary condition (`chd`, `ghb`, `drn`, `riv`, `wel`, `rch`, `evt`) comes in three
forms:

```python
# (a) direct — you supply FloPy stress-period data:
ghb = mf.ghb(stress_period_data={0: [[(1, c), 86.0, 50.0] for c in mouth_cells]})

# (b) GIS-driven — read features from a GeoPackage layer:
drn = mf.drn.gpkg("gis/bcs.gpkg", layer="springs", context=ctx, nper=1)

# (c) raw FloPy escape hatch — when you already have prepared data:
rch = mf.rch.flopy(stress_period_data=rch_data)
```

See [Chapter 9, *Boundary Conditions*](09_boundary_conditions.md).

### Step 6 — the surface-water network

Streams from a centerline, lakes from polygons, an unsaturated zone — each a
single call against the context. Pull the SFR and LAK specs out as handles so the
mover can reference them **semantically**:

```python
sfr = mf.sfr(
    context=ctx, nper=1, streams="gis/streams.gpkg",
    connection_mode="automatic",
    connections=(mf.StreamConnection("north_trib", "main_stem"),
                 mf.StreamConnection("south_trib", "main_stem")),
    width=18.0, gradient=0.0012, roughness=0.030,
    streambed_k=0.05, streambed_thickness=1.5,
    inflow={0: {"north_trib": 15000.0, "south_trib": 10000.0}},
    length_conversion=3.28081, time_conversion=86_400.0, mover=True,
)

lak = mf.lak(
    context=ctx, nper=1, lakes="gis/lakes.gpkg", lake_id_field="name",
    starting_stage={"valley_lake": 101.0}, lake_bottom={"valley_lake": 96.0},
    bed_leakance=0.11, connection_modes="automatic",
    status={"valley_lake": ["ACTIVE"]}, mover=True,
    length_conversion=3.28081, time_conversion=86_400.0,
)

uzf = mf.uzf(context=ctx, nper=1, cells=floor_cells,
             vks=0.25, thtr=0.08, thts=0.34, thti=0.17, eps=4.0,
             finf={0: [3.0e-5] * len(floor_cells)},
             pet={0: [1.0e-4] * len(floor_cells)}, extdp=7.0)
```

See [Chapter 10, *Advanced Packages*](10_advanced_packages.md).

### Step 7 — assemble the model

One flat, declarative package list — the readable payoff:

```python
flow = mf.gwf(
    "valley",
    context=ctx,
    newtonoptions="UNDER_RELAXATION",
    save_flows=True,
    packages=[
        mf.disv(nlay=2, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                vertices=gp["vertices"], cell2d=gp["cell2d"],
                top=layers.top, botm=layers.botm, idomain=layers.idomain,
                length_units="FEET"),
        mf.ic(strt=strt),
        mf.npf(k=k, k33=k * 0.1, save_flows=True),
        mf.sto(steady_state={0: True}),
        ghb, drn, rch, uzf, sfr, lak,
        # MVR: move half the main stem's OUTFLOW (its final reach, resolved from
        # geometry — no hard-coded reach numbers) into the lake.
        mf.mvr(nper=1, moves=(
            mf.Move(mf.sfr_connection(sfr, "main_stem"),
                    mf.lak_connection(lak, "valley_lake"),
                    value=0.5),
        )),
        mf.oc(head_filerecord="valley.hds", budget_filerecord="valley.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ],
)
```

> **Warning — MVR ordering.** Packages moved by MVR must be **declared in the
> model and ordered before the mover**. myflopy validates this. The semantic
> endpoints — `mf.sfr_connection(sfr, "main_stem")` and
> `mf.lak_connection(lak, "valley_lake")` — resolve a stream's outlet reach and a
> lake's number from names/geometry, so you never hand-code indices. Pass
> `at=(x, y)` to target the reach nearest a coordinate.

### Step 8 — simulation and run

```python
sim = mf.SimulationSpec(
    "baseline",
    models=(flow,),
    packages=(
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=("valley",), complexity="COMPLEX",
               outer_maximum=250, inner_maximum=250,
               linear_acceleration="BICGSTAB"),
    ),
)
project.add_simulation(sim)

run = project.prepare_run("baseline", "baseline")   # build in memory
success, report = run.execute()                     # write + run MF6
print("converged:", success, "->", run.workspace)
```

That is a complete, GIS-driven, surface-water-coupled model in one screen of
declarative code.

---

## 3.4 Reading results, plotting, and calibrating — at a glance

### Results

`run.model(name)` returns a model **view** with myflopy helpers on top of the
FloPy output:

```python
model = run.model("valley")

heads = model.hds.array(layer=0)        # final-step heads for layer 0 (ncpl,)
steps = model.hds.kstpkper             # available (kstp, kper) output times
```

To compare against field measurements, wrap them in an observation **target** and
let it align with the simulation:

```python
from myflopy.modflow.mf6.observations import HeadTargets

obs = HeadTargets(locations="gis/wells.gpkg", values="data/observed_heads.csv")
obs.compare(model)     # per-target observed vs. simulated + residuals
obs.stats(model)       # n, mean error, MAE, RMSE
```

See [Chapter 13, *Reading Results*](13_results.md).

### Visualization

Every picture comes from one of a handful of verbs, on the model itself or on
`myflopy.plot`; figure/axis plumbing goes through the `myflopy.viz` front door
(figures, subplots, a shared palette). A layer head map is one call:

```python
model.plot.map(layer=0)                     # interactive Plotly map
model.plot.map(layer=0).plot_mpl(title="Layer 1 head (ft)")   # static matplotlib
```

The verbs are `map` (plan view), `section` (vertical slice), `surface` (3-D),
`grid` (the bare mesh), plus `mosaic` and `animate`. They exist both bound —
`model.plot.map(...)`, `vor.plot.grid()` — and free, as
`plot.map(model, ...)`; they are the same functions either way. Everything they
return is a *Picture*: it renders itself in Jupyter and answers `.fig`,
`.show()`, `.save(path)` and `.html(path)`.

For something you can share without a kernel, export a **standalone HTML slider**
that pages through time:

```python
from myflopy.modflow.mf6.interactive_plotting import export_head_map_slider_html

export_head_map_slider_html(model, "valley_heads.html", layer=0)
```

See [Chapter 14, *Visualization*](14_visualization.md).

### Calibration (a peek)

myflopy's PEST layer is a declarative facade that **compiles to native pyEMU
`PstFrom`**. The whole loop — parameterize, observe, forecast, build, run — is a
handful of calls off `model.pest(...)`:

```python
cal = model.pest("valley_pest", start_datetime="2024-01-01")

cal.parameterize("k",        style="pilotpoints", bounds=(0.1, 10.0))
cal.parameterize("recharge", style="constant",    bounds=(0.5, 2.0))

cal.observe(obs)               # the HeadTargets from above
cal.forecast(forecast_targets) # the predictions you actually care about

cal.build("valley.pst", noptmax=0)            # compile to PstFrom + forward run
cal.run_ies(reals=50, iterations=3, workers=12)   # PESTPP-IES, 12 parallel agents
```

Then review the ensemble — phi, observed-vs-simulated, forecasts, and parameter
field maps:

```python
handle = model.pest_runs[-1]   # discover calibrations done on this model
results = handle.review()      # -> IesResults
results.plot_phi()
```

Parameter **styles** are `constant`, `zone`, `grid` (one geostatistically
correlated multiplier per Voronoi cell), and `pilotpoints` (IDW from a pilot net).
See [Chapter 17, *PEST / pyEMU*](17_pest.md).

---

## 3.5 The one-screen API map

Everything you reach for, by job. All are attributes of `import myflopy as mf`
unless noted.

**Project & run lifecycle**

| Call | Purpose |
|---|---|
| `mf.Project(root, name=)` | Durable workspace; owns runs/scenarios. |
| `project.add_simulation/add_package/add_grid` | Register reusable specs/libraries. |
| `project.prepare_run(name, sim)` → `Run` | Build in memory (inspect before running). |
| `run.execute()` → `(success, report)` | Write + run MF6. |
| `project.run(name, sim)` → `Run` | One-shot build + execute. |
| `run.model(name)`, `run.flopy_model(name)` | Model view (myflopy helpers) / raw FloPy. |
| `mf.load_run(workspace)` | Reopen an existing run from disk. |

**Geometry**

| Call | Purpose |
|---|---|
| `mf.GridSpec.voronoi(boundary=, refinement=)` | Declarative grid recipe from GIS. |
| `grid_spec.resolve(workspace=)` → `VoronoiGridPlus` | Build the grid eagerly. |
| `mf.LayerStack(vor, top=).add(...).build(attach=True)` | Layered surfaces → top/botm/idomain. |
| `mf.ModelContext(grid=, domain=, surfaces=, dates=)` | The geometry bundle (rides on the model). |
| `mf.Raster / mf.Contours / mf.Array` (from `myflopy.layers`) | Surface atoms for `LayerStack`. |

**Models & simulation**

| Call | Purpose |
|---|---|
| `mf.gwf / mf.gwt / mf.gwe / mf.prt` | One flow / transport / energy / particle model. |
| `mf.SimulationSpec(name, models=, packages=)` | tdis + ims + models (+ exchanges). |
| `mf.tdis / mf.ims` | Time discretization / solver. |

**Packages**

| Call | Purpose |
|---|---|
| `mf.disv / mf.ic / mf.npf / mf.sto / mf.oc` | Core flow packages. |
| `mf.chd / mf.ghb / mf.drn / mf.riv / mf.wel / mf.rch / mf.evt` | List BCs — each `()`, `.gpkg()`, or `.flopy()`. |
| `mf.uzf / mf.sfr / mf.lak` | Unsaturated zone / streams / lakes. |
| `mf.mvr` + `mf.Move` + `mf.sfr_connection / mf.lak_connection` | Water mover with semantic endpoints. |

**Results, viz, calibration**

| Call | Purpose |
|---|---|
| `model.hds.array(layer=)`, `model.hds.kstpkper` | Heads. |
| `HeadTargets / LakeStageTargets / Sfr*Targets / DrnFlowTargets` | Observation targets (`compare/stats`). |
| `build_choropleth(model.vor, custom_zs=).plot_mpl()`, `export_*_slider_html(...)` | Maps and interactive sliders. |
| `model.pest(name).parameterize / observe / forecast / build / run_ies` | Calibration. |
| `model.particle_tracking.prt(workspace=, release_points=)` | Particle tracking (MF6 PRT). |
| `mf.ParallelModelWorkflow(model).split_model(...)` | Domain-decomposed parallel run. |

---

## 3.6 Where to go next

- **Understand the design** → [Chapter 4, *Architecture*](04_architecture.md) and
  [Chapter 5, *The Spec System*](05_specs.md).
- **Build from your own GIS** → [Chapter 6, *Grids*](06_grids.md),
  [Chapter 7, *Layers*](07_layers.md), and
  [Chapter 9, *Boundary Conditions*](09_boundary_conditions.md).
- **Add streams, lakes, and movers** →
  [Chapter 10, *Advanced Packages*](10_advanced_packages.md).
- **Organize work and scenarios** →
  [Chapter 12, *Projects, Runs & Scenarios*](12_projects_runs.md).
- **Visualize and calibrate** → [Chapter 14, *Visualization*](14_visualization.md)
  and [Chapter 17, *PEST / pyEMU*](17_pest.md).
- **When something breaks** →
  [Chapter 19, *Troubleshooting & FAQ*](19_troubleshooting.md).

> **The copy-me template.** The full, tested version of the §3.3 model lives at
> [`examples/mf6/package_first_full_stack.py`](../../examples/mf6/package_first_full_stack.py).
> Start there for a real build.
