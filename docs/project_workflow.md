# Projects, Reusable Packages, and Swappable Grids

A `Project` is the single workspace for related MODFLOW 6 simulations. It owns
two reusable libraries — **packages** and **grids** — that are defined once,
independently of any model, and referenced (or swapped) across simulations.
Persistence is opt-in via `save()`.

## Quickstart

```python
import myflopy as mf

project = mf.Project("modeling/cumberland", name="cumberland")

# 1. Package library: define packages once, independent of any model.
project.add_package("npf/base",   mf.npf(k=10.0))
project.add_package("npf/high_k", mf.npf(k=100.0))
project.add_package("rch/2020",   mf.rch.gpkg("inputs/rch.gpkg", context=ctx, nper=12))

# 2. Grid library: register a grid recipe or an already-built grid object.
project.add_grid("base", vor)              # a VoronoiGridPlus you built and liked

# 3. A model references packages and a grid by key.
flow = mf.gwf("flow",
    grid=mf.grid_ref("base"),
    packages=[mf.ref("npf/base"), mf.ref("rch/2020"), mf.oc(saverecord=[("HEAD", "ALL")])],
)
baseline = mf.SimulationSpec("baseline",
    models=[flow],
    packages=[mf.tdis(nper=12, perioddata=perioddata), mf.ims(models=["flow"])],
)
project.add_simulation(baseline)

# 4. Build and run. `run(run_name, simulation)` writes + executes into runs/<run_name>/.
#    `simulation` may be a SimulationSpec or a registered simulation name.
run = project.run("baseline", "baseline")
m = run.model("flow")                      # a ModelView: m.hds, m.packages, m.visualize ...
```

`mf.ref("npf/base")` and `mf.grid_ref("base")` are the reuse idioms. A
`PackageRef` / `GridRef` is resolved against the project libraries when the run
is built, so nothing is materialized until you build.

## Scenarios: reuse and swap

Derive a named variant and swap a package and/or a grid; everything else is
reused. The change is recorded in `simulation.lineage`.

```python
refined_hk = (
    project.simulation("baseline")
    .derive("refined_hk")
    .replace_package("flow", mf.ref("npf/high_k"))   # swap a package
    .replace_grid("flow", mf.grid_ref("refined"))    # swap the grid
)
project.add_simulation(refined_hk)
project.run("refined_hk", "refined_hk")
```

GIS-driven packages (`mf.rch.gpkg(...)`, `mf.drn.gpkg(...)`, `K`-from-vector,
etc.) re-map onto whatever grid resolves, because they are evaluated against
the model grid at build time. Inputs tied to concrete cell IDs or raw per-cell
arrays do **not** re-map — keep swappable inputs source-driven.

## Grids: build, look, then register

Grid generation with `TriangleGrid` / `VoronoiGridPlus` is iterative and
visual, so the primary path is to build the grid yourself and register the
finished object.

```python
tg = mf.TriangleGrid(...)
tg.set_domain_file("inputs/domain.gpkg")
tg.add_region_file("inputs/wellfield.gpkg", maximum_area=400)   # refine an area
vor = tg.build_mesh(profile="balanced")
vor.show()                                                      # look at it
project.add_grid("refined", vor)                                # register when happy
```

To **refine an area**, build a second grid with the refinement, look at it, and
register it under a new key — then swap to it with `replace_grid`.

### Grid definition tiers

| How you define the grid | When to use it | Persisted by `save()` as |
|---|---|---|
| Built object — `project.add_grid(key, vor)` / `mf.GridSpec.from_object(vor)` | Your default: build, look, use | a **pickle** (a regenerable cache) plus a version sidecar |
| `mf.GridSpec.python("grids.py", function="build_grid")` | A settled, reproducible build recipe | JSON pointing at your builder |
| `mf.GridSpec.voronoi(...)` / `.structured(...)` / `.from_geopackage(...)` | Simple or already-existing grids | JSON recipe |

A built-object grid is held in memory and is **not** serializable on its own
(`GridSpec.from_object(...).to_dict()` raises). It is pickled only when you call
`project.save()`. Pickle is a cache, not an archive: `GridSpec.python` is the
reproducible, version-robust path. `Project.load()` warns if a pickle's recorded
`flopy` / `geopandas` / `shapely` / `numpy` versions differ from the current
stack.

## Persistence is opt-in

In memory, everything works immediately — including raw arrays and coupled
models. Nothing is written until you call `save()`.

```python
project.save()                       # writes specs/, packages/, grids/ under the root
later = mf.Project.load("modeling/cumberland")
later.run("refined_hk", "refined_hk", overwrite=True)
```

### Serialization boundary

`save()` requires serializable specs. It can persist simulations whose package
data is JSON-friendly or expressed as data sources, the package library, and
the grid library. It **cannot** yet persist:

- **Simulations with exchanges** (coupled GWF↔GWT, GWF↔PRT, …).
- **Specs that embed raw arrays** — e.g. an `npf` carrying a NumPy `k` array, or
  a `wel` with explicit per-cell stress data.

`project.validate()` reports exactly which simulation fails and why, instead of
failing silently. For durable projects, express data as source specs
(`mf.RasterSource`, `mf.GeoPackageSourceSpec`, …) rather than materialized
arrays; otherwise keep the project in memory.

## Customizing the on-disk layout

`Project` writes under `root` using an internal `ProjectLayout`. To customize
directory names, pass one:

```python
from myflopy.workspace import ProjectLayout

layout = ProjectLayout(root, specs_dir_name="recipes", packages_dir_name="pkgs")
project = mf.Project(root, name="cumberland", layout=layout)
```
