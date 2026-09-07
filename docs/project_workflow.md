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
m = run.model("flow")                      # a ModelView: m.hds, m.packages, m.plot ...
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

## Coupled (multi-model) runs

When a simulation has more than one model — for example GWF coupled to GWT,
GWE, or PRT — each model's input files are written into its own subdirectory
named after the model, so it is clear which files belong to which model:

```
project/runs/high_k/
├── mfsim.nam
├── gwf/            ← flow model files
└── gwt/            ← transport model files
```

This is automatic for multi-model simulations. Single-model runs stay flat
(files at the run root). A model that sets its own `model_rel_path` is left as
declared. Building and running coupled models works today; note that a
simulation containing an exchange cannot yet be persisted with `save()` (see
the serialization boundary below) — though its built MF6 files remain
reopenable with `mf.load_run(...)`.

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
| `mf.GridSpec.voronoi(...)` | A Voronoi mesh from GIS boundary/refinement | JSON recipe |
| `.structured(...)` / `.from_geopackage(...)` | **Constructors exist but FAIL FAST** — only `voronoi` and `python` resolve. Use `GridSpec.from_object` for a built or structured grid | — |

A built-object grid is held in memory and is **not** serializable on its own
(`GridSpec.from_object(...).to_dict()` raises). It is pickled only when you call
`project.save()`. Pickle is a cache, not an archive: `GridSpec.python` is the
reproducible, version-robust path. `Project.load()` warns if a pickle's recorded
`flopy` / `geopandas` / `shapely` / `numpy` versions differ from the current
stack.

## Did it actually work?

A run manifest's `status` is **myflopy's verb** — what the library did to the run
(`created` → `built` → `written` → `completed`). It says nothing about what MF6
did afterwards, so a run written by one session and executed by another sits at
`status="written", success=None` forever.

`run.outcome` is **MF6's verdict**, read fresh from `mfsim.lst`:

```python
run.outcome
# ListingOutcome(outcome='normal', elapsed='3 Hours, 53 Minutes,  1 Seconds', stale=True)

run.outcome.succeeded        # True / False / None
run.reconcile()              # adopt the verdict into status/success and persist
```

Four outcomes: `normal`, `failed` (MF6 reported premature termination),
`incomplete` (a listing exists but reports neither — killed, crashed, or still
running), and `absent` (no listing; MF6 never started here).

**`succeeded` is `None`, not `False`, for an unfinished run.** Calling an
interrupted run a failure would be the same class of lie as calling
`status="written"` a success — so `reconcile()` leaves such a run alone rather
than recording a verdict MF6 never gave.

**`stale=True` means the listing predates its own inputs** — the verdict is real,
but describes a superseded run. Only MF6 *inputs* count; results are legitimately
newer than the listing that announced them.

The outcome is deliberately **never cached into the manifest**. The listing is
the authority and it is right there; storing a copy beside it would create a
second source of truth that can disagree with the file it came from.

## Seeing what a project holds

`project.packages`, `.grids` and `.simulations` are **`dict` subclasses** — every
dict operation works exactly as before, and each also answers the same three
questions:

```python
project                                  # Project('tentrails', root='…',
                                         #         packages=14, grids=1, simulations=3)
project.packages.names()                 # ['chd/base', 'crisp/base', 'dis/base', …]
project.packages.summary()               # key | package | builder | enabled | options | …
project.package("npf/base")              # == project.packages["npf/base"]
project.summary()                        # all three libraries, one frame
```

A missing key names the near miss and everything defined, on the spelling you
already type:

```python
project.packages["npf/bass"]
# KeyError: Project 'tentrails' has no package 'npf/bass'. Did you mean: 'npf/base'?
#           14 packages: chd/base, crisp/base, dis/base, …
```

`project.packages.persistence()` reports how each package *would* be persisted
(`json` or `pickle`, with the reason). It is a method rather than a `summary()`
column because deciding it means attempting `to_dict()` — 119 ms on a
100 000-record list BC — and Jupyter calls `_repr_html_` on every bare `project`.

**Library keys are not attribute-completable and never will be**: they contain
`/`, so they are not Python identifiers. They are answered at runtime instead —
by `names()`, `summary()`, the repr, and that `KeyError`.

`repr(project)` touches no disk and counts no runs: runs live on the filesystem,
and PyCharm calls `__repr__` on every debugger step. Use `project.discover_runs()`
for those.

## Forking and editing packages

Two verbs, spelled apart because they differ in one important way: `fork` leaves
the original alone, `edit` changes what an existing key **means**.

```python
project.packages.fork("rch_pond", "rch_pond_vol",
                      name="rch_pond_vol", recharge=pond_flow_vol,
                      note="max-volume pond")

project.packages.edit("npf/base", k=res.per_layer(K_revised),
                      note="post-calibration K")
```

Both are **built, not patched** — the new spec goes through the same constructor
`mf.npf(...)` does, so a typo raises here rather than at MF6 build time:

```python
project.packages.fork("npf/base", "npf/high", kk_typo=1.0)
# TypeError: npf: ModflowGwfnpf does not accept 'kk_typo'. Did you mean 'k'?
```

Both **deep-copy by default**, including options the fork never touched. That is
not caution for its own sake: `with_options` rebuilds the options dict but not its
values, so without the copy, editing a "copy" writes through into the original's
arrays.

Cost scales with the number of **Python objects** copied, not with cell count.
Arrays copy at memcpy speed — a 9,405-cell `npf` forks in 0.1 ms, a 72×9,405
recharge array in 1.7 ms — but a record list does not: the USG-derived
72-period, 9,090-record EVT takes **808 ms** for its 654,480 records.
`deep=False` is 0.1 ms on that same package and opts back into sharing.

`fork` refuses an existing key unless `overwrite=True`; `edit` is never refused.

Simulations fork the same way, swapping packages by library key:

```python
project.simulations.fork("baseline", "mounding",
                         swap={"rch/series": "rch_pond"})
```

The flat form applies to **every model carrying that label** — right for a
single-model simulation and for genuinely shared packages. A coupled GWF/GWT
simulation has `ic`, `oc` and `dis` on both models, so name the model when you
mean one:

```python
project.simulations.fork("coupled", "hot_start",
                         swap={"tport": {"ic": "ic/warm"}})
```

A label matching nothing raises, naming what the source actually carries — a typo
that silently swapped nothing would hand back a "variant" identical to its parent.
The fork records `derived_from` and appends to `lineage`, one entry per swap.

## The run ledger

`project.runs` is every run under the project — with a manifest, without one, or
broken:

```python
project.runs.summary()
#  key                  status      outcome     stale  elapsed                          simulation
#  baseline             completed   normal      True   3 Hours, 53 Minutes,  1 Seconds  baseline
#  from_usg             written     incomplete  False  —                                tentrails
#  mounding             completed   normal      False  15 Minutes,  4.265 Seconds       mounding
#  mounding_lower_L3_K  discovered  normal      False  14 Minutes, 32.796 Seconds       —
#  mounding_vol         completed   normal      False  16 Minutes, 54.890 Seconds       mounding_vol
```

Read `status` and `outcome` as two axes: `status` is what **myflopy** did to the
run, `outcome` is what **MF6** did afterwards. A run written by one session and
executed by another sits at `written` forever while its listing says
`Normal termination`; `run.reconcile()` adopts the verdict.

`stale=True` means an MF6 input file is newer than the listing — the verdict is
real but describes superseded inputs.

Runs with no `run.json` are listed as `discovered` rather than hidden: there are
far more bare MF6 workspaces in practice than manifests, and a ledger that
omitted them would be worse than none. A workspace under `runs/` is keyed by its
directory name either way; one elsewhere keeps its path, so `runs/from_usg` and
`scratch/from_usg` stay distinct. A truncated manifest is a row reading
`unreadable`, never an exception — one bad file must not hide every other run.

### What an edit invalidated

The manifest records package **keys**, not their contents, so
`project.packages.edit("npf/base", ...)` leaves every earlier run still naming
`npf/base` — true of the key, false of the model it actually used. Each run also
records a fingerprint of what those keys meant at the time:

```python
project.runs.drift()
#  run       key         ran_with          library_now       state
#  baseline  npf/base    8b6bcfeb950cb270  0e83cf9f816b6102  edited
#  baseline  rch/series  333df8c89afe0a37  333df8c89afe0a37  unchanged
#  mounding  rch_pond    8a2aacef0f6bca5d  —                 removed
```

`summary()` carries the one-cell version as a `drift` column naming the moved
keys. A run reading `unchanged` throughout is reproducible from the library as it
stands; anything else names exactly what to look at.

Only `mf.ref` entries are fingerprinted — a package written inline is baked into
the run's own recipe and cannot drift. Drift only reports; it never rewrites a
run's stored history.

Five states, and the four that are not `unchanged` all mean "do not assume":

| state | meaning |
|---|---|
| `unchanged` | reproducible from the library as it stands |
| `edited` | the key still exists and means something else |
| `removed` | the key is gone from the library |
| `unfingerprinted` | the package cannot be fingerprinted (a closure builder, an option that does not pickle) — it is still perfectly buildable |
| `unknown` | the run predates digests; a missing fingerprint is not evidence of sameness |

An `outcome` of `superseded` means the listing in that directory is OLDER than
the run's own preparation, so its verdict belongs to a previous run.
`prepare_run(overwrite=True)` reuses the directory and leaves the old listing and
results behind, so without this a re-prepared, never-executed run reported the
previous execution's `normal` and its elapsed time. `reconcile()` refuses to
adopt a superseded verdict.

The fingerprint covers the option and input VALUES a spec holds, its name, its
builder REFERENCE, `enabled` and `requires`. It does not see the body of an
importable builder, or the contents of a file a source points at — it says the
recipe is unchanged, not that the world is. Dict key order, `requires` order and
list-vs-tuple are normalised away, so a `save()`/`load()` round trip does not
read as an edit.

The library is **cached**, so `project.runs["baseline"]` returns the same object
twice; call `project.refresh_runs()` after executing a run from elsewhere.
Opening one run by name stays `project.reopen_run(name)`, which also attaches the
build context so a revived spec resolves its library references.

## Persistence is opt-in

In memory, everything works immediately — including raw arrays and coupled
models. Nothing is written until you call `save()`.

```python
report = project.save()              # writes specs/, packages/, grids/ under the root
report                               # SaveReport(project_spec='project_spec.json',
                                     #            simulations=2, packages=14, grids=1,
                                     #            pickled=9, skipped=0)
later = mf.Project.load("modeling/cumberland")
later.run("refined_hk", "refined_hk", overwrite=True)
```

`save()` returns a **`SaveReport`**, not a path. It is `os.fspath`-able — so
`Path(report)` and `open(report)` work — but it is deliberately not *equal* to a
`Path`. `pickled` is the number worth watching: a pickled entry is a cache, not
an archive. `skipped` pairs each unpersisted simulation with the reason.

### Serialization boundary

`save()` persists what it can and reports the rest. A simulation it cannot
serialize is **skipped** — written into neither the project document nor its own
file — and named in `report.skipped`; the package and grid libraries are
persisted regardless. Before 2026-09-04 a single such simulation refused the
whole project, so a library of otherwise-saveable packages wrote nothing at all.
Pass `strict=True` for that older behaviour (raise before writing anything).

It can persist simulations whose package data is JSON-friendly or expressed as
data sources, the package library, and the grid library. It **cannot** yet
persist:

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
