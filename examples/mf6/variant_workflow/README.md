# Variant Workflow Examples

Four runnable examples showing how to keep model variants reusable, trackable,
and tidy with `myflopy`. Examples 01–03 use a small structured grid as a stand-in
mesh; example 04 builds a real Voronoi/DISV grid from layer surfaces. The
patterns are identical for unstructured models.

```
python 01_package_library_and_variants.py
python 02_grid_aware_packages.py
python 03_coupled_gwf_gwt.py
python 04_layer_surfaces_to_disv.py
```

- **`concerns.py`** — the reusable building-block functions. The "messy" array
  logic (your `get_ks(vor)`, raster sampling, zonal masks) lives here, so the
  model/variant files stay short.

## The one rule

> The reuse unit follows what is fixed. Anything **computed on a grid** (a `k`
> array, cell-mapped boundaries, layer top/botm) is **bound to that grid**.

| You vary… | Reuse unit | Pattern |
|---|---|---|
| Stresses, properties | a **package** | pickle once, `mf.ref` (example 01) |
| Layer top/bottom (same mesh) | the **DIS package** | a package variant, `mf.ref` |
| The **mesh** (refine) | grid + recomputed packages | grid-aware builder (example 02) |

## 01 — package library + variants (fixed mesh)

Compute the expensive packages once, `project.save()` (array packages are
pickled to `specs/packages/*.pkl`, scalar/source ones to JSON), then define
variants that reuse most packages and swap one `mf.ref(...)`. Reload in a fresh
session with `mf.Project.load(...)` and reuse without recomputing.

A variant is one readable line: which package it references.

## 02 — grid-aware packages (varying mesh)

When the mesh varies you cannot reuse a pickled array (wrong size). Instead the
package is a **builder function** that computes from the grid at build time, so
it recomputes for whatever mesh the model is built on. This is the right tool
for mesh-refinement variants and for source-driven inputs.

## 03 — coupled GWF–GWT

The exchange is cheap assembly that lives in the variant function; the reusable
**flow** packages come from the library. Each model is written into its own
subdirectory automatically (`runs/<name>/gwf/`, `runs/<name>/gwt/`).

Coupling notes baked into the example: each model needs its **own IMS**
(transport uses `BICGSTAB`); a flow boundary requires an **SSM** package on the
transport model with an auxiliary concentration on the boundary.

## 04 — layer surfaces → a real DISV grid

Build a Voronoi mesh once, describe the layering as a stack of `Surface`s
(`flat` / `raster` / `from_contours` / `from_points`), and turn it into
discretization with one call:

```python
layers = mf.LayerSurfaces([top, layer1_botm, layer2_botm])
disv   = layers.to_disv(vor)   # samples + reconciles overlaps -> mf.disv spec
```

Each surface is an independent, swappable ingredient: the example builds two
models on the **same mesh** that differ only in one layer surface (flat vs.
point-interpolated dip), and both solve. Swap in `Surface.raster("ground.tif")`
or `Surface.from_contours("aq_botm.gpkg", z="elev", region_raster="dom.tif")`
unchanged. `to_disv` also writes `vor.gdf_topbtm`, so choropleth/mapping code
keeps working.

## What persists vs. what stays in code

- **Persisted (reusable artifacts):** computed packages and grids — pickled by
  `save()`, reused across variants and sessions via `mf.ref` / `mf.grid_ref`.
- **In code (cheap assembly):** which models, which references, which exchange —
  your small `variant()` / `build()` function, version-controlled.

You never serialize the arrays themselves to JSON; you pickle them as artifacts
and keep the assembly in readable Python.
