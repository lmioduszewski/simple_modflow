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
    .with_packages(mf.ref("npf/high_k"), slot="npf", on="flow")  # swap a package
    .replace_grid("flow", mf.grid_ref("refined"))                # swap the grid
)
project.add_simulation(refined_hk)
project.run("refined_hk", "refined_hk")
```

**Swapping the grid re-maps only what was declared LAZILY.** A GIS-driven helper given
a `context=` (`mf.drn.gpkg(path, context=ctx, …)`, `mf.rch(context=ctx, …)`,
`K`-from-vector, …) resolves its cells when you **declare** it, not when the model
builds, so its records carry the cell numbers of the grid that was on the context at
declaration. Omit `context=` and the same declaration becomes a RECIPE that resolves
against whatever grid the build supplies: the seven `.gpkg` forms (ledger 279) and
`mf.rch(...)` / `mf.evt(...)` (ledger 284). `mf.uzf/sfr/lak` are still eager-only.

Since 2026-09-13 that is refused rather than run: `ModelSpec.build` compares the
built DISV, the context domain and every package's cellids against the grid, and
raises naming what is stale. It used to build and run — MF6 applied the old cell
numbers to the new geometry and reported Normal termination. **Read
`sim.move_report()` before `prepare_run`**: a package declared lazily reads `rebuilt`
and moves itself, while one declared with `context=` reads `stale` and must be
re-declared against the new grid.

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

> **Runnable tour:** `examples/mf6/notebooks/project_libraries_tour.ipynb` walks
> this whole surface on a two-cell model in a temp directory — define, browse,
> fork, edit, run, reconcile, drift, save, reload. It executes in seconds and a
> test runs it, so it cannot rot silently.

## Seeing what a project holds

`project.packages`, `.grids` and `.simulations` are **`dict` subclasses** — every
dict operation works exactly as before, and each also answers the same three
questions:

```python
project                                  # Project('tentrails', root='…',
                                         #         packages=14, grids=1, simulations=3, aoi=2)
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

## Areas of interest: places to look at

`project.aoi` -- and `model.aoi` on any model of the project -- is a fourth library:
named **areas of interest** (ledger 310). Not model geometry; places you keep going
back to. A name then works wherever a place does:

```python
project.add_aoi("site", (500, 600, 1000, 1100))            # (xmin, ymin, xmax, ymax), grid CRS
project.add_aoi("wetland", "gis/wetland.gpkg", layer="aoi")  # any geometry, a file, a GeoDataFrame
model = run.model("flow")
model.add_aoi("north", north_polygon)                      # defined on the model...
sorted(project.aoi)                                        # ...kept by the project: ['north', 'site', 'wetland']
model.plot.map(layer=0, extent="site", backend="mpl")      # the window a map opens on
model.packages.npf.transmissivity.mosaic(per=3, extent="site")
model.plot.map(layer=0, select="wetland")                  # the cells under it
model.plot.map(layer=0, features=model.aoi)                # every area, one layer each
```

**Defined on the model, kept by the project.** A model that came from this project --
`prepare_run`, `reopen_run`, `project.runs[...]`, or `mf.load_run` of one of its runs
while the `Project` is alive -- answers `model.aoi` with the project's own library, so
`model.add_aoi(...)` IS `project.add_aoi(...)` and `project.save()` keeps it. With no
live `Project` object, a model reads the project's SAVED areas, and anything it adds
lasts for the session only: only `project.save()` writes a project's files. A model in
no project keeps its areas in memory, like `model.regions`.

Saved INLINE in `project_spec.json` (WKT + CRS), so an archive needs no file to restore
them. `crs=None` means the grid's frame. A name that is both an area and a model
region is refused by `select=` rather than guessed. `project.aoi["x"] = geom` replaces
as a dict does; `add_aoi` refuses an existing name unless `overwrite=True`.

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

**An edit VOIDS the grid pin it speaks for** (ledger 275). A resolved `.gpkg`
package carries a certificate saying the library resolved its records on a named mesh;
an override replaces that option wholesale, so after an edit the records are yours and
the certificate is a lie. The pin is voided and the package reads `unverifiable` rather
than `stale` — and if you move the model, the refusal names the edit instead of telling
you to re-declare records you just typed yourself. Only options the pin *speaks for*
void it: editing `stress_period_data` does, renaming the file or toggling `save_flows`
does not. There is no way back to a pin through the edit route — you cannot pin cell
numbers the caller typed (ledger 250) — so declare the package fresh against the new
grid when you want one.

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
#  run       key         ran_with          library_now       state      source
#  baseline  npf/base    8b6bcfeb950cb270  0e83cf9f816b6102  edited     —
#  baseline  rch/series  333df8c89afe0a37  333df8c89afe0a37  unchanged  ok
#  baseline  chd         —                 —                 —          edited
#  mounding  rch_pond    8a2aacef0f6bca5d  —                 removed    —
```

The third row is a package declared **inline** rather than through `mf.ref`. It has
no library key, so `state` — *does this key still mean what the run used?* — has
nothing to answer about and says so; `source` is the column that answers for it.

`summary()` carries the one-cell version as a `drift` column naming the moved
keys. A run reading `unchanged` throughout is reproducible from the library as it
stands; anything else names exactly what to look at.

Its `grid` column says what each model's grid was. A `GridSpec` is a recipe, and the
grid it resolves to changed on 2026-09-25 without the recipe changing: myflopy now
repairs the cells FloPy draws wrong at reflex corners and on the boundary (ledger 316).
So a run records each model's grid pin and how its cells were made (`run.grids`), and
the column reads `unrecorded` for a run from before that (it may sit on FloPy's cells,
holes and all), `m: FloPy's cells` for one built with `tessellation="flopy"`, and
`m: 4 cells repaired` where the repair changed something. `—` means there is nothing
to say.

**A run rebuilds from its own copy, not from the library** (ledger 230).
`prepare_run` copies each library package the simulation references into
`<run>/_library/`, and a rebuild resolves `mf.ref` from those copies first. So
`edited` means the library has moved on -- preparing the run again from it would
give a different model -- not that this run can no longer be rebuilt as it ran.
A run prepared before the copies existed has nothing else to build from: it
resolves against the library as it is now, warns naming each key that differs
from what it recorded, and records the versions it used once its inputs are
rewritten.

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

**`drift()` also reports a `source` column** (ledger 276) — whether the FILE a
package was resolved from still says what it said. The two answer different
questions and a package is routinely both at once: `state` is *does this key still
mean what the run used?*, `source` is *is the file it was built from still the one
it was built from?* Values are `ok`, `edited`, `missing`, `unfingerprinted` (names a
source but records no hash — anything declared before this landed) and `—` for the
majority that name no file. It hashes BYTES, so rewriting a GeoPackage without
changing a feature reads as `edited`: a prompt to look, which is the safe direction.

**It answers about what the RUN read, not about the library's current declaration**
(ledger 294). Each run records `metadata['built_from']` — the content hash of every
source it actually opened — and `source` compares against that. Without it, two runs
of one simulation prepared either side of an edit were indistinguishable: they share
one library declaration, so they recorded the identical fingerprint while building
different models, and both got the same verdict, including the run that AGREES with
the file. It also covers packages declared **inline**, which had no row at all,
because `package_digests` records only `mf.ref` entries.

For an **eager** package `source` is PROVENANCE, not buildability: it resolved at
declaration and holds records, so it builds fine with its GeoPackage deleted —
`missing` there means the file you built from has moved, not that anything is
broken. Only a lazy package is re-hashed at build, because only a lazy package reads
a file there. A package naming several files (a layer stack names one per surface
plus the top, and an excavation its footprint file and bed) gets one cell, reporting
the WORST with a count: `edited (1 of 3)`.

The fingerprint covers the option and input VALUES a spec holds, its name, its
builder REFERENCE, `enabled` and `requires`. It does not see the body of an
importable builder, or the contents of a file a source points at — it says the
recipe is unchanged, not that the world is. Dict key order, `requires` order,
list-vs-tuple and a `Path`-vs-string spelling are normalised away, so a
`save()`/`load()` round trip does not read as an edit.

A `Surface`, a `LayerStack` or a source is fingerprinted by its **declaration** —
what its `to_dict()` writes — because that is the only definition under which the
live spec and the reloaded one can agree (ledger 273). Two consequences worth
knowing: what `to_dict` omits is invisible, so previewing a stack with `.plot` no
longer marks its package edited, and two stacks bound to different grids
fingerprint alike; and a stack whose top is a per-cell `Surface.from_array` cannot
be written down at all, so it reads `unfingerprinted` rather than carrying a
digest that would drift on its own. An excavation (`stack.excavate`) is part of that
declaration (ledger 323): cutting a pond into a stack moves the digest of a lazy
`mf.disv(stack)`, and the footprint file is hashed by content like the stack's other
sources, so a redrawn pond reads `edited` under `source`.

**PEST calibrations are reported, not listed.** A calibration builds into
`<model workspace>.pest/` — inside the project — holding the PstFrom template,
the IES and prior masters, and every pyEMU worker. Each is a wholesale copy of
the model workspace, so each carries an `mfsim.nam`; listing them would add a
dozen rows per calibration to a ledger whose job is telling you what you ran.
They are excluded and counted instead:

```python
project.runs
# RunLibrary(1 runs: 'baseline'; 1 calibration excluded -- see .calibrations)

project.runs.calibrations
#  calibration          workspaces  open_with
#  runs/baseline.pest            6  model.pest_runs
```

Nothing is hidden — the repr says how many, and `model.pest_runs` /
`mf.find_pest_runs` open them properly, returning handles that understand phi,
ensembles and forecasts. `project.discover_native_runs()` is unchanged and still
returns every native MF6 workspace, calibration copies included; it is the raw
sweep, where `project.runs` is the curated ledger.

The library is **cached**, so `project.runs["baseline"]` returns the same object
twice; call `project.refresh_runs()` after executing a run from elsewhere.
Opening one run by name stays `project.reopen_run(name)`, which also attaches the
build context so a revived spec resolves its library references -- from the
run's own copies first. `mf.load_run` has no library in scope, so there package
references resolve from those copies and nothing else.

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

`context_dropped` pairs each model with the context fields a project file cannot
hold -- `grid`, `domain`, `surfaces` -- which reload as `None`. That is not new
behaviour, but it used to be silent, and it matters: the build-time grid check and
the per-package grid pins have nothing to compare against once the grid is gone. A
16-character fingerprint of the cell NUMBERING is recorded and compared instead, so
a different mesh is still caught; one with the same connectivity and different
coordinates (a rescale, a mirror, a unit conversion) is not. To keep the full check,
register the grid and reference it:

```python
project.add_grid("coarse", vor)
flow = mf.gwf("flow", grid=mf.grid_ref("coarse"), packages=[...])
```

A `GridRef` is serialized on the model itself and resolved at build time, so it
survives the round trip where a grid object on the context does not.

`bundled` lists the files `save(bundle_sources=True)` copied into `inputs/`.

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
