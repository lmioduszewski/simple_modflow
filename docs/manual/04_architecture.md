# 4. Architecture & the Mental Model

[Chapter 3](03_getting_started.md) showed you *what* a myflopy program looks
like. This chapter explains *why* it is shaped that way. Understanding these five
ideas will make the rest of the manual feel obvious:

1. The **package-first API** is the one canonical way to build models.
2. Every high-level helper is a thin **facade over an engine over atoms**.
3. A model is a tree of **immutable specifications** that a `Project` materializes.
4. **Geometry rides on the model**, not on the project.
5. Because specs are immutable and composable, **variants and reusable packages**
   are nearly free.

---

## 4.1 The package-first API is the canonical surface

myflopy grew two generations of model-building code. The one you should use — and
the only one this manual teaches for assembling models — is the **package-first
API**: the flat namespace of helpers on `import myflopy as mf`.

```python
import myflopy as mf

flow = mf.gwf("valley", context=ctx, packages=[
    mf.disv(...), mf.ic(...), mf.npf(...), mf.sto(...),
    mf.ghb(...), mf.drn(...), mf.rch(...),
    mf.uzf(...), mf.sfr(...), mf.lak(...), mf.mvr(...),
    mf.oc(...),
])
```

Each `mf.<package>(...)` call returns a **package specification** — a small,
immutable, serializable description of one MODFLOW 6 package. You assemble a model
by handing `mf.gwf` a flat list of them. That is the whole idea: *declare what the
model is made of; let the library build it.*

> **Note — what about the builder classes?** You will see `SFRBuilder`,
> `LAKBuilder`, `UZFBuilder`, `RCHBuilder`, `MVRBuilder`, and friends in the
> source and in [Chapter 20](20_api_reference.md). These are **not** a competing
> API — they are the *engine underneath* the package-first facade (see §4.2). New
> work uses `mf.sfr(...)`, not `SFRBuilder(...)`.

---

## 4.2 The facade → engine → atom pattern

The package-first helpers are deliberately thin. Underneath almost every one is a
richer object that does the real work, and sometimes a third layer of primitives.
The same three-layer pattern repeats across the library:

| Facade (you use this) | Engine (does the work) | Atoms (primitives) |
|---|---|---|
| `mf.uzf(...)` | `UZFBuilder` | — |
| `mf.sfr(...)` | `SFRBuilder` | `StreamConnection`, `StreamDiversion` |
| `mf.lak(...)` | `LAKBuilder` | `LakeOutlet`, `LakeTable` |
| `mf.mvr(...)` | `MVRBuilder` | `Move`, `MoverConnection` |
| `mf.LayerStack(...)` | `LayerSurfaces` | `Surface` (`Raster`/`Contours`/`Array`/…) |

Concretely, `mf.uzf(context=ctx, ...)` literally constructs a `UZFBuilder`,
resolves the unsaturated cells against the grid in the context, and returns the
`PackageSpec` that `UZFBuilder.build()` produces. The facade exists to give you a
**readable, declarative call**; the engine exists to hold the **resolution logic**;
the atoms exist for the cases where you need to describe one stream connection or
one surface precisely.

You normally live entirely at the facade layer. You drop to the engine only when
you need something the facade doesn't surface — which is rare, and called out
explicitly where it applies.

```
mf.sfr(context=ctx, streams="streams.gpkg", ...)   ← FACADE: one readable call
        │
        ▼
   SFRBuilder(...).build()                          ← ENGINE: resolves reaches → cells
        │            ▲
        ▼            │
   PackageSpec   StreamConnection(...)              ← ATOMS: precise pieces
```

This is why [Chapter 7](07_layers.md) tells you to use `LayerStack` and treats
`LayerSurfaces`/`Surface` as the engine/atoms: same pattern, different package.

---

## 4.3 The object hierarchy

A myflopy model is a nested tree of specifications. From the outside in:

```
Project                  durable workspace + run/scenario lifecycle (holds NO geometry)
  └─ SimulationSpec      one MF6 simulation: tdis + ims solver + the model(s) + exchanges
       └─ ModelSpec      one GWF / GWT / GWE / PRT model
            ├─ ModelContext   grid + domain(idomain) + surfaces + dates   ← geometry
            └─ packages       a flat tuple of PackageSpec
```

Each layer has a single job:

- **`Project`** ([Chapter 12](12_projects_runs.md)) is the durable workspace. It
  owns runs and the reusable package/grid libraries, and it manages the lifecycle
  (`prepare_run` → `Run.execute`). It holds **no geometry and no model data** — it
  is the filing cabinet, not the model.
- **`SimulationSpec`** ([§5.3](05_specs.md)) is one MODFLOW 6 simulation: time
  discretization (`tdis`), a solver (`ims`), the model(s), and any inter-model
  exchanges. It is what you register and run.
- **`ModelSpec`** ([§5.2](05_specs.md)) is one model. You create it with
  `mf.gwf/gwt/gwe/prt(...)`; it carries the package list, the context, and
  post-build hooks.
- **`ModelContext`** ([§5.4](05_specs.md)) is the geometry bundle — grid, active
  domain, surfaces, dates — that the package builders read.
- **`PackageSpec`** ([§5.1](05_specs.md)) is one package.

Everything in this tree except the built FloPy objects is a frozen dataclass, which
is what makes §4.5 and §4.6 possible.

---

## 4.4 Geometry rides on the model, not the project

This is the single most important structural decision in myflopy, and the one most
likely to trip up someone coming from a project-centric tool.

**The `Project` holds no geometry.** Grids, layers, the active domain, and surfaces
all live in a **`ModelContext`** that you attach to the *model*:

```python
ctx = mf.ModelContext(grid=vor, domain=idomain, surfaces=vor.gdf_topbtm)
flow = mf.gwf("valley", context=ctx, packages=[...])
```

After a build, the model exposes its context as `model.myflopy_context`, and the
GIS-aware package helpers (`mf.uzf`, `mf.sfr`, `mf.lak`, and the `.gpkg` boundary
forms) read the grid from `context=` to map features onto cells.

Why put geometry on the model rather than the project?

- A project can hold **several models on different grids** (a coarse regional model
  and a refined local one) without a single "project grid" forcing them together.
- A coupled simulation's GWF, GWT, and PRT models can **share one context** by
  passing the same `ctx` to each — the geometry is stated once and reused, but it
  still belongs to the models, not to an ambient project state.
- It keeps the `Project` a pure lifecycle object, which is what makes runs and
  scenarios cleanly reproducible.

> **Tip.** When you build a multi-physics simulation ([Chapter 11](11_multiphysics.md)),
> construct one `ModelContext` and pass it to every `mf.gwf/gwt/...` call. They all
> resolve against the same grid and domain.

---

## 4.5 Immutability and composition

Every spec is a **frozen dataclass**. You never mutate one; you derive a new one
with a `with_*` method. This *copy-on-write* style is the foundation of variants
(§4.6) and safe reuse.

`PackageSpec` carries four such methods:

```python
npf = mf.npf(k=10.0)

npf_high = npf.with_options(k=100.0)          # change a FloPy option
npf_src  = npf.with_inputs(k=mf.RasterSource("k.tif"))   # swap an input source
npf_tag  = npf.with_metadata(note="calibrated")          # annotate (non-FloPy)
# npf itself is unchanged; each call returns a NEW PackageSpec.
```

The model and simulation specs compose the same way:

```python
# ModelSpec: add or replace a package in its slot, swap the grid or context.
flow2 = flow.with_package(npf_high)            # replace npf; everything else reused
flow3 = flow.with_grid(other_grid)
flow4 = flow.with_context(other_ctx)

# SimulationSpec: swap a whole model, or read one out and put a new one back.
sim2 = sim.with_model(flow2)
npf  = sim.model("valley").package("npf")      # read a package out of the tree
```

Because each operation returns a fresh, independent spec and shares the unchanged
parts, building a family of related models is cheap and side-effect-free — the
baseline is never disturbed by a variant.

> **Note — when immutability has a cost.** Specs that embed **raw NumPy arrays**
> (an `npf` with a `k` array, a `wel` with explicit per-cell data) work perfectly
> in memory but cannot be serialized to disk. For durable, swappable inputs,
> express data as **source specs** (`mf.RasterSource`, `mf.GeoPackageSourceSpec`,
> …) instead of materialized arrays. See [Chapter 18](18_serialization.md).

---

## 4.6 Model variants & reusable, swappable packages

Real projects are not one model — they are a **baseline plus a family of
variants**: high-K vs. low-K, a refined grid, a drought recharge scenario, a
calibrated parameter field. myflopy makes this a first-class workflow built on two
ideas: **libraries** (define a thing once) and **copy-on-write derivation** (swap
one thing, reuse the rest).

### Libraries: define packages and grids once

A `Project` owns two reusable libraries. You register entries under string keys,
independent of any model:

```python
project = mf.Project("modeling/valley", name="valley")

# Package library — each key maps to a PackageSpec.
project.add_package("npf/base",   mf.npf(k=10.0))
project.add_package("npf/high_k", mf.npf(k=100.0))
project.add_package("rch/2020",   mf.rch.gpkg("gis/rch.gpkg", context=ctx, nper=12))

# Grid library — register a grid you built and liked (or a recipe).
project.add_grid("base",    vor)
project.add_grid("refined", vor_refined)
```

### Reference library entries by key

A model refers to library entries with `mf.ref(key)` (a `PackageRef`) and
`mf.grid_ref(key)` (a `GridRef`). **Nothing is materialized at declaration** — a
reference is resolved against the project libraries when the run is built:

```python
flow = mf.gwf(
    "flow",
    grid=mf.grid_ref("base"),
    packages=[mf.ref("npf/base"), mf.ref("rch/2020"),
              mf.oc(saverecord=[("HEAD", "ALL")])],
)
baseline = mf.SimulationSpec(
    "baseline", models=[flow],
    packages=[mf.tdis(nper=12, perioddata=perioddata), mf.ims(models=["flow"])],
)
project.add_simulation(baseline)
```

This keeps a single source of truth: change `npf/base` once and every simulation
that references it picks up the change at the next build.

### Derive a variant and swap one piece

To make a scenario, **derive** a named copy of a simulation and `replace_package`
and/or `replace_grid`. Everything you don't touch is reused, and the edit is
recorded in `simulation.lineage`:

```python
refined_hk = (
    project.simulation("baseline")
    .derive("refined_hk")
    .replace_package("flow", mf.ref("npf/high_k"))   # swap a package by key
    .replace_grid("flow",   mf.grid_ref("refined"))  # swap the grid
)
project.add_simulation(refined_hk)
project.run("refined_hk", "refined_hk")
```

`derive`/`replace_*` are the high-level, library-aware idioms. They are sugar over
the copy-on-write methods from §4.5 — under the hood `replace_package` does
`sim.with_model(model.with_package(...))` and appends a lineage entry. When you are
working with concrete specs rather than library keys, the low-level form is just as
valid:

```python
high_k_flow = baseline.model("flow").with_package(
    baseline.model("flow").package("npf").with_options(k=100.0)
)
high_k = baseline.with_model(high_k_flow)
```

### What re-maps on a grid swap — and what doesn't

Swapping the grid in a variant is powerful, but only **source-driven** inputs
follow the new grid automatically:

- ✅ **Re-maps.** GIS-driven packages (`mf.rch.gpkg`, `mf.drn.gpkg`,
  K-from-vector, …) are evaluated against the model grid *at build time*, so they
  re-resolve onto whatever grid the variant uses.
- ❌ **Does not re-map.** Inputs tied to concrete cell IDs or raw per-cell arrays
  (a `chd` with explicit `(layer, cell)` records, an `npf` carrying a NumPy `k`
  array sized to the old grid) are bound to the grid they were written for.

> **Tip.** If you expect to swap grids, keep swappable inputs **source-driven**
> (GeoPackage/raster sources) rather than materialized arrays. That single habit is
> what makes grid-refinement scenarios painless. See
> [Chapter 12, §12.5](12_projects_runs.md).

---

## 4.7 The build lifecycle

Nothing touches the disk or FloPy until you build. The lifecycle is:

```
spec tree  ──build(context)──▶  FloPy objects  ──write──▶  MF6 input files  ──run──▶  outputs
(immutable)                     (a Simulation)             (on disk)                  (.hds/.cbc)
```

Step by step:

1. **Declare.** You compose `PackageSpec` / `ModelSpec` / `SimulationSpec`
   in memory. References (`mf.ref`, `mf.grid_ref`) and deferred grids are still
   unresolved.
2. **Prepare.** `run = project.prepare_run(name, sim)` resolves library
   references, resolves the grid into the context, and calls each spec's
   `build(...)` to produce the FloPy objects — all **in memory**. You can inspect
   `run.model(name)` before anything is written.
3. **Execute.** `run.execute()` (or the one-shot `project.run(name, sim)`) writes
   the MF6 input files into `runs/<name>/` and runs the `mf6` executable, returning
   `(success, report)`.
4. **Read.** `run.model(name)` now exposes results (heads, budgets, package output)
   — see [Chapter 13](13_results.md).

> **Note — the eager/deferred seam.** Most grids resolve at step 2. But the
> GIS-aware package helpers resolve their cells **eagerly, at declaration** (step
> 1), because `UZFBuilder.build()` bakes resolved cell data into the `PackageSpec`.
> That is why the full surface-water stack needs a concrete grid in the context
> *before* you call `mf.uzf/sfr/lak`. The mechanics and the workaround are in
> [Chapter 6, §6.5](06_grids.md).

---

## 4.8 The `mf.*` import surface

`import myflopy as mf` exposes a curated, lazily-imported namespace. It is
intentionally split:

- **Preferred surface** (`myflopy.__preferred__`) — the package-first helpers and
  the core spec/grid/context types. This is what `dir(myflopy)` and tab-completion
  show, and what this manual teaches: `mf.gwf`, `mf.disv`, `mf.ghb`, `mf.sfr`,
  `mf.Project`, `mf.SimulationSpec`, `mf.ModelContext`, `mf.LayerStack`, …
- **Compatibility surface** (`myflopy.__compatibility__`) — the engine classes and
  raw `*_spec` factories (`SFRBuilder`, `LAKBuilder`, `ghb_spec`, …). These remain
  importable for power users and internal use but are kept out of the discovery
  surface so newcomers are steered to the facade.

Imports are lazy: `myflopy/__init__.py` maps each public name to its module and
resolves it on first access through `__getattr__`. The practical upshot — you pay
the import cost of (say) the visualization stack only if you actually touch
`mf.plot_model_head_map`, which keeps `import myflopy as mf` fast.

```python
import myflopy as mf

mf.gwf            # preferred — resolves myflopy.package_api.gwf
mf.SFRBuilder     # compatibility — the engine under mf.sfr, still importable
```

---

## Where this leaves you

You now have the conceptual scaffolding for everything that follows:

- the **specs** in detail → [Chapter 5](05_specs.md);
- building each part of a model → Chapters [6](06_grids.md)–[11](11_multiphysics.md);
- the project, run, and **scenario** lifecycle in practice →
  [Chapter 12](12_projects_runs.md).
