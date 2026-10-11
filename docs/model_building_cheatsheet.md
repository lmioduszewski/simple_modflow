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

**Exactly one region is required** — `region_vector=` (a polygon, normally your
model domain) or `region_raster=` (a raster whose extent to use). Giving neither
or both raises at construction. `z=` names the elevation attribute on the contour
lines; it defaults to `"Elev"`.

**Prefer `region_vector`.** `region_raster` uses the WHOLE raster, and a DEM
usually covers far more ground than your model. Measured on a real project: a
47,105 × 46,525 ft LiDAR as the region gave 3.5 M cells at `resolution=25` — 7.3×
the area of the domain, and 6.4× the area the contours covered. The extra area is
the worst case for `r.surf.contour`, which flood-fills outward from contour
lines and has nothing to work from out there. Switching to the domain polygon
took each layer to ~30 s.

**The output is clipped to the region already.** `clip=True` is the default, and
with `region_vector` it nulls everything outside the POLYGON, not merely its
bounding box — so a domain-shaped raster with nodata in the bbox corners. Verified
on a real domain: raster bounds equal to the domain's, 64% of pixels non-null.
Pass `clip=False` only when a neighbouring model needs the overlap.

**myflopy does not reproject for you.** GRASS is given `-o` (override projection
check), so contours and region are placed at their raw coordinates and a CRS
mismatch is silently accepted, not corrected. Contours in EPSG:2926 with a
region raster in EPSG:2927 landed 597,000 ft apart — `v.to.rast` then rasterized
nothing inside the region, `r.surf.contour` had no seed cells, and it scanned
forever at 0%. Reproject first:

```python
gpd.read_file(src).to_crs("EPSG:2927").to_file(dst, driver="GPKG")
```

The interpolated GeoTIFF lands beside the contours as `<name>.interp.tif` and is
reused on later runs. Pass `out=` to place it elsewhere.

#### When does the interpolation actually run?

Not when you declare the surface. Measured by counting calls to the interpolator:

| | interpolations |
|---|---|
| `mf.Contours(...)` + `.add(...)` the stack | **0** — fully lazy |
| `stack.draft_grid()` / a gridless `stack.plot` | **1** — it needs the raster's *bounds* |
| `build()` · `qc()` · `to_disv()` · later plots | still 1 — cached |

So a gridless preview is **not** cheaper than a build: the only way to learn where
a contour surface is, is to interpolate it. After that the GeoTIFF is reused
across builds, sessions and processes. `surface.cache_status()` reports
`fresh`/stale/absent; `refresh=True` on `build`/`qc`/`to_disv` forces a rebuild.

#### Watching it work

```python
ground = mf.Contours("ground.gpkg", z="Elev", region_vector="domain.gpkg",
                     progress=True)
```

```
Reading features...
  11%  22%  33%  44%  55%  66%  77%  88% 100%
Writing raster map...
   0%   4%   8%  12%  16% ...
```

**In a terminal you already see this** without the flag — GRASS modules are
subprocesses writing to the inherited file descriptors. **In Jupyter you do not**:
those descriptors belong to the kernel's console, not the cell, and
`contextlib.redirect_stdout` cannot catch them because the writing happens below
Python entirely. `progress=True` swaps the descriptors for a pipe and re-emits
through `sys.stdout`, which is what the cell displays.

It is a **named parameter, deliberately not part of `**grass_kwargs`** — those
feed the derived-raster cache signature, so a display flag riding through them
would silently re-interpolate the first time you asked to watch. Verified:
toggling `progress` leaves `cache_status()` at `fresh`.

### Every surface source

| source | use |
|---|---|
| `mf.Contours(path, z="Elev", region_vector=...)` | digitized elevation contours (GRASS) |
| `mf.Raster(path, fill=…, nodata=…)` | a DEM or any GeoTIFF |
| `mf.Points(xs, ys, zs, method="linear")` | scattered control points — **the GRASS-free fallback** |
| `Array(values)` | one value per cell, already computed (`from myflopy.layers import Array`) |
| `mf.Surface.from_array(values)` | same, via the atom class |

#### If your DEM was clipped but not labelled

A raster clipped to a boundary and exported **without a no-data value in its
header** is a very common GIS product — the fill (usually `0`) is then read as a
real elevation. It shows up as a rectangular surface instead of a domain-shaped
one, ground at elevation zero outside the boundary, and every layer beneath
collapsing to `min_sep` once its contact gets capped to that zero.

```python
mf.Raster("clipped.tif", nodata=0)      # declare it at read time
```

Measured on a real project: 4.3 M pixels of exact `0.0` outside the domain, and
median layer thicknesses of `[0.10, 0.10, 0.10]`. With `nodata=0` they became
`[6.62, 21.33, 58.60]` — identical to rewriting the file properly.

Masking happens *before* the per-cell average, so a cell straddling the boundary
is the mean of its real pixels rather than a blend with the sentinel.

The durable fix is still to label the raster, which helps every other tool too —
in QGIS, **Clip Raster by Mask Layer** with *"Assign a specified nodata value"*
set to `-9999` (never `0`; it is a plausible elevation). Check it took with
`gdalinfo out.tif | grep -i nodata`.

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

#### Looking at it before any grid exists

`.plot` needs no grid at all. On a gridless stack it draws on a coarse **draft
grid** covering the surfaces' own extent, so you can check the layering before
committing to a mesh:

```python
stack.plot.section(y=1000)     # works with vor is None
stack.plot.map()
stack.plot.surface("all")
```

```
INFO  LayerStack.plot: no grid on this stack, drawing on a draft grid of 574
      cells over (0, 0)-(3000, 2000). Pictures only -- build()/qc()/to_disv()
      still need a real grid.
```

The extent comes from the first `Raster` or `Contours` surface in the stack —
those are the only kinds that know where they are, and the search descends
through algebra, so a raster buried in `Flat(150).capped_at(ground - 2)` is still
found. The draft is built once (~0.05 s for a few hundred cells) and cached.

**Pictures only, deliberately.** `build()`, `qc()` and `to_disv()` still demand a
real grid: an approximate picture is useful, invented model geometry is not. A
draft has no boundary polygon, no refinement, and whatever extent your DEM
happens to cover rather than your model domain — which is exactly why it must
never become a DISV.

Control it, or supply an extent when nothing carries one:

```python
stack.for_grid(stack.draft_grid(cells=2000)).plot.section(y=1000)   # finer
stack.draft_grid(extent=(0, 0, 500, 400))                            # all-Flat stack
```

With no `crs=` the draft takes its raster's CRS, else the package default
`EPSG:2927` (ledger 391).

A stack of only `Flat`/`Array` surfaces describes thicknesses with no notion of
where they are, so `.plot` raises pointing at `extent=` rather than inventing
one.

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

#### Refinement: three channels, split by geometry type

Refinement is declared on the **spec**, not at resolve time, and there are three
separate inputs — they are *not* symmetric, which is the part that trips people up:

| | how many sources | geometry | default buffer | default priority | sizing option |
|---|---|---|---|---|---|
| `refinement=` | **one** | polygons | `0` | `0` | `refinement_max_area` |
| `breaklines=` | **a list** | lines → buffered into corridors | `10` | `1` | `breakline_max_area` |
| `points=` | **a list** | points | — | — | none — locations only |

Within any source, **every feature becomes its own region.** So "several refinement
areas" means several polygons in one layer, not several files. Give each its own
cell size with a column:

```python
refinement=mf.GeoPackageSourceSpec(
    refine, layer="areas",
    fields={"area": "max_area", "label": "name", "priority": "rank"},
),
refinement_max_area=5_000.0,     # fallback for rows with no max_area
```

`fields` reads **`{logical_key: your_column}`** — the key is what myflopy looks up,
the value is the column in your layer. Only `area`, `label` and `priority` are
understood, and an unrecognised key is ignored silently, so check the spelling
first if a per-feature size seems not to apply. Only `GeoPackageSourceSpec` carries
`fields`; a `ShapeSource` has no field map, so a shapefile can only take one size
for all its polygons via `refinement_max_area`.

**Polygons and lines in the same layer will not work.** A line reaching
`refinement=` raises `A LineString region has no area; pass buffer>0 to give it
width`. Setting `refinement_buffer` fixes it but buffers *every* feature in the
source — measured, a 400×400 ft polygon in a buffered layer came out 67% larger.
Split by geometry type instead, which can still be one file:

```python
mf.GridSpec.voronoi(
    boundary=mf.ShapeSource(domain, crs=CRS),
    refinement=mf.GeoPackageSourceSpec(refine, layer="areas",
                                       fields={"area": "max_area", "label": "name"}),
    breaklines=[mf.GeoPackageSourceSpec(refine, layer="creeks")],
    breakline_max_area=800.0,
    breakline_buffer=40,          # half-width of the corridor, CRS units
    crs=CRS, boundary_max_area=20_000.0,
)
```

**Several polygon geopackages.** `refinement=` is singular — a list fails with
`AttributeError: 'list' object has no attribute 'path'`. Either merge them into one
layer with an `area` column, or use `breaklines=` as the multi-source channel with
the buffer off, since `buffer(0)` leaves a polygon unchanged:

```python
refinement=mf.GeoPackageSourceSpec(polys_a, fields={"area": "max_area"}),
breaklines=[mf.GeoPackageSourceSpec(polys_b, fields={"area": "max_area"}),
            mf.GeoPackageSourceSpec(polys_c, fields={"area": "max_area"})],
breakline_buffer=0,
```

Buffer and the priority *default* are global per channel, so all breakline sources
share them; per-feature `area`/`priority` columns still override.

**All of these are named parameters**, so your editor completes them — `boundary_*`,
`refinement_*`, `breakline_*` and `region_point_tolerance`. An unrecognised keyword
now raises with the nearest match:

```
GridSpec.voronoi() got an unknown option 'breakline_bufer'.
Did you mean 'breakline_buffer'?
```

That used to be accepted and never read, so the typo silently meshed at the default
— measured, `breakline_bufer=40` gave corridors a quarter of the intended width
under a clean run. The wrong *shape* of source is refused the same way: a list in
`refinement=`, a single source in `breaklines=`, or a live GeoDataFrame anywhere all
name the fix instead of failing later as `'X' object has no attribute 'path'`.

**The cells tile the boundary.** FloPy's `tri2vor` left a hole at a reflex (inside)
corner of the boundary -- every 270° corner measured -- sometimes dropped a vertex lying
exactly on a straight edge, and kept cells in only the first part of a boundary drawn as
two features; a missing cell at a corner was the visible sign, and nudging the boundary
until it went away was the usual dodge. `VoronoiGridPlus` now checks every cell against
Triangle's own mesh and repairs the wrong ones (ledger 316), and a grid with nothing wrong
is unchanged byte for byte. `vor.tessellation_report.summary()` says what was repaired;
`GridSpec.voronoi(..., tessellation="flopy")` reproduces a grid made before the repair.
A boundary with a hole in it (an island, a lake cut out) builds now (ledger 317); one
that touches itself -- a hole meeting the edge, two parts meeting at a corner -- is
refused by name, because the cell there would be two pieces.

**Keep refinement polygons clear of the boundary.** A refinement polygon that crosses the
boundary -- or ends exactly on it, which is what clipping it to the boundary yourself
produces -- is clipped to a slightly inset domain, and Triangle meshes the thin sliver
between: 6,903 nodes against 493 for the same box kept 50 ft clear (the L, `max_area`
5000, refinement 1000), and 44,324 around a hole. Stop refinement polygons about a cell
width short of the edge.

**A point on the boundary becomes a boundary vertex.** A `points=` point strictly inside
and nearer an edge than twice what the file writes (0.02 ft at state plane) -- a well
snapped onto the boundary in GIS, or reprojected through latitude and longitude -- is
placed ON that edge; one that close to a corner or another point is merged into it, and
a warning says how many moved and how far (ledger 316). Left alone, Triangle bisected the
edge toward it: it crashed, hung or ran out of memory, or built a grid that could not be
written. A point OUTSIDE the boundary still gets no cell -- it is warned about now -- and
a point on the drawn boundary can round outside it, or sit outside the 2 ft simplified
one: move a well you need inside by more than that. Farther from an edge a point still
grades small cells toward it -- a 1000 ft square at state plane, max_area 2000, is 684
cells alone and 853 with a well 0.05 ft inside, 788 at 0.1 ft, 1,005 at 0.5 ft
(`docs/reproducers/gridspec_build_failures.py`); that is legal, and kept.

**Sources must be file-backed.** A `GridSpec` is a serializable recipe — it has to
survive `to_dict()` / `from_yaml()` — so an in-memory GeoDataFrame cannot be one.
If you built your refinement polygons programmatically, write them out first:

```python
gdf.to_file(proj.root / "inputs" / "refine.gpkg", layer="areas", driver="GPKG")
```

Writing under the project root also means the relative path anchors via
`project_root`. If you truly want no file, drop to `TriangleGrid.add_region_polygon`
in a loop over the frame and give up the durable spec.

**Trap: overlapping regions resolve by priority, and the defaults differ.** Regions
sort by `(-priority, area)` and each claims only what higher-priority regions left.
At equal priority the smaller claims first, so nesting works — inner polygon takes
its area, outer gets the donut. But breaklines default to priority `1` and
refinement to `0`, so a corridor takes the overlap from a refinement polygon:

```
('wellfield', 150.0, 'line',       geom 160000, claim 160000)   # took the whole overlap
('a',         400.0, 'refinement', geom 250000, claim 209995)   # ceded 40005 ft²
```

That is usually right for a stream and wrong when you used `breaklines=` only to
pass a second polygon file — set `priority` explicitly there. A region left with no
unique area raises `does not have any unique interior area left after resolving
overlaps`, which is what a duplicated polygon gets you -- except an EXACT duplicate
breakline row (the same line, `max_area` and priority), which is dropped with a warning.

**Check it before meshing.** `build=False` reads the sources and resolves regions
without running Triangle:

```python
tri = grid_spec.resolve(project_root=project.root, build=False)
[(r.label, r.max_area, r.source, round(r.claim_geometry.area))
 for r in tri._prepared_regions]
```

`source` says which channel each region came from (`refinement` / `line` /
`domain`), and `geometry` vs `claim_geometry` is exactly what overlap resolution
removed. After a real build, `_triangle.0.poly` in the workspace is the file
Triangle was actually given — its region section lists each point and its max area.

Then get a real grid object out of it — **still declared in the project**:

```python
vor = project.grids["base"].resolve(
    project_root=project.root,                     # relative paths anchor here
    workspace=project.root / "grids" / "base",     # where Triangle's files land
)
```

`project_root` anchors **every relative path in the recipe** — source paths,
`idomain_path`, `rasters` — so a spec written with project-relative paths stays
portable. Omit it and they resolve against the current working directory instead.
`workspace` is where Triangle's `_triangle.*` files are written; give each grid
its own directory, because the file names are fixed and two grids sharing one
overwrite each other's meshes. Reading `_triangle.0.poly` there is the honest way
to see what Triangle was actually asked for — its region section lists every
refinement point and the maximum area attached to it.

Two more arguments, both for inspection:

```python
tri = grid_spec.resolve(project_root=project.root, build=False)   # no meshing
[(r.label, r.max_area) for r in tri._prepared_regions]            # what Triangle will get
vor, tri = grid_spec.resolve(..., return_triangle=True)           # keep both
```

`build=False` reads the sources, resolves regions against one another and picks
their interior points, then stops — the cheap way to check a recipe before paying
for a mesh. Full argument-by-argument reference is in `GridSpec.resolve.__doc__`.

**Trap: the CRS has a hardcoded fallback.** If neither the spec nor its options
name one, the grid is built as `"EPSG:2927"` (Washington State Plane South, feet).
Always pass `crs=` to `GridSpec.voronoi`.

**Why resolve now rather than defer?** Because a `LayerStack` and the GIS package
helpers (`mf.sfr`/`mf.lak`, and `mf.X.gpkg`/`mf.uzf` when given `context=`) resolve cells
the moment you call them —
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
`units`, `nlay`, `vor`; methods `per_layer()`, `qc()`, `report()`, `validate()`,
`prune_isolated()`, `attach_to_grid()`, and the `plot` namespace. A stack with an
excavation (§5, "Excavating a pond") also fills `excavations` (name →
`mf.ExcavationRecord`: its cells and sampled bed), `ground` (the top before the cut),
`excavated` and `absorbed`, read by `excavation_report()`; `excavate()` cuts a result
that has no stack behind it.

#### `split=` — one geologic unit, several model layers

Geology and discretization are different decisions. `.add()` states the geology;
`split=` says how finely MODFLOW should slice it. The unit's **geometry does not
change** — the last cut lands exactly on the declared bottom.

```python
stack.add("sand", bottom=Raster("sand_base.tif"), split=3)      # three equal slices
stack.add("till", thickness=60, split=[0.25, 0.75])             # 15 ft, then 45 ft
stack.replace("sand", split=4)                                  # or decide later
```

Shares must sum to 1 (`[30, 70]` is rejected, naming the fix). The layers are
named `sand_1 … sand_N`; a unit with no split — or `split=1` — keeps its bare
name, so adding a split later never renames anything else.

**`names=` if the generated ones don't read well.** One name per split layer,
top to bottom:

```python
stack.add("sand", bottom=Raster("sand_base.tif"), split=3,
          names=["upper sand", "mid sand", "lower sand"])
```

It is a labelling choice and nothing else: same contacts, same thicknesses, same
`units` map — the UNIT is still `"sand"`, so `per_layer({"sand": …})` and
`plot.grid(layers="sand")` are unaffected by what you called the slices. The list
must cover **every** layer, because naming only some of them would leave which
layer a name refers to depending on where you stopped counting. `names=` without
a `split=` of 2 or more is refused (an unsplit unit already has a name — rename
the unit). `split` and `names` are validated as a pair, including inside
`replace`, so re-splitting a named unit means restating the names and dropping
the split reads `replace("sand", split=None, names=None)`.

**The unit stays addressable.** This is the part that matters, because splitting
changes `nlay` and every per-layer argument downstream is a *positional list*:

```python
layers = stack.build(vor)
layers.units                       # {'fill': [0], 'sand': [1,2,3], 'clay': [4]}
mf.npf(k=layers.per_layer({"fill": 30.0, "sand": 25.0, "clay": 0.05}))
layers.plot.grid(layers="sand")    # a unit name still selects its layers
```

Use `per_layer()` rather than writing `k=[30, 25, 25, 25, 0.05]` by hand. A list
of the *wrong* length is at least rejected by FloPy; a stale list of the **right**
length after a split is accepted silently with new meaning. A layer name as a key
overrides its unit, so one slice can differ from its siblings.

**`min_thickness` and `pinch` stay unit-scoped.** A unit pinches out whole or not
at all. Judged per slice, a 2.5 ft unit split three ways came back
`idomain [0, 1, 0]` — inactive cells inside a unit that is fully present, blocking
vertical flow. Splitting is a numerical choice and must not change that verdict.

Two limits worth knowing:

- A unit whose bottom is an `Isopach` **cannot be split** — it is measured from
  the layer above, so each cut would measure from the previous *cut* and the
  contacts would drift downward instead of dividing the unit. Give it an absolute
  bottom, or declare it with `thickness=` and split that.
- Fixed lifts (*"20 ft layers"*) are not a `split=` mode, because a unit of
  varying thickness would need a varying number of layers and `nlay` is global.
  For a `thickness=`-declared unit just divide: a 60 ft unit in 20 ft lifts is
  `split=3`. For a surface-bounded unit, repeat a constant cut instead:
  ```python
  lift = Surface.constant_thickness(20).floored_at(sand_base)
  stack.add("sand_1", bottom=lift).add("sand_2", bottom=lift).add("sand_3", bottom=sand_base)
  ```

Underneath it is one new surface operation, `mf.Toward(target, fraction)` —
the contact a share of the way down from the surface above to `target`. You rarely
need it directly; `split=` also handles the renormalization that makes equal
shares actually equal.

#### `pinch=` — what happens where a layer goes thin

| value | effect |
|---|---|
| `"floor"` *(default)* | cell stays active, held open at its `min_thickness` |
| `"passthrough"` | `idomain = -1` — cell inactive but vertically transmissive, and it can carry **no boundary condition** |
| `"inactive"` | `idomain = 0` — a true pinch-out, no flow through |
| `"floor"` | clamp thickness to `min_thickness`, cell stays active |

#### Erosion — a surface that cuts down through the layers below

An incised channel, a buried valley, a scoured contact: the overlying surface
drops *below* one or more contacts underneath it. **Erosion needs no separate "cut"
operation** — `reconcile="bottom"` already does it. It pushes every contact the
cut passes through down to `cut - min_sep`, cascading, and the layers that end up
thinner than `min_thickness` are handled by their `pinch`. (A VOID — a pond, pit or
basin that will hold a lake — is different, and is `stack.excavate(bottom, within)`:
see "Excavating a pond" below.)

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

`min_thickness` is a real minimum under the default `pinch="floor"`: reconcile holds
the layer open at it, so `.add("sand", thickness=0.4, min_thickness=5)` builds 5 ft,
not 0.4. `min_sep` is only the fallback for layers that pinch. A split unit divides
its minimum among its slices by their shares, so the unit keeps the thickness you
declared. Reconcile never reaches a PINCHING unit at or above its own `min_thickness`
-- a 0.8 ft clay with `min_thickness=0.5, pinch="inactive"` builds 0.8 ft and active --
and a unit below it is pinched, crushed to `min_sep` only where it is also within
`trigger_sep` (1.0 by default): a 1.5 ft clay with `min_thickness=2` stays 1.5 ft and
reads idomain 0 (ledger 377).

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

#### Excavating a pond — `stack.excavate(bottom, within)`

A lake's cells are WATER, not aquifer: without the cut, a `bathy` lake on a stack whose
top is the ground connects its floor into the layer CONTAINING the bed, and MF6 moves
that vertical connection up to the cell's top — on a real pond 411 of 612 vertical
connections exchanged nothing all run (ledger 327). The cut makes the first active cell under the
bed start exactly at the bed:

```python
bed = mf.Contours("pond_contours.gpkg", z="elevation", region_vector="pond.gpkg")
stack = (mf.LayerStack(top=ground)
         .add("sand", bottom=lyr1, min_thickness=3)
         .add("clay", bottom=lyr2, split=[0.5, 0.5]))
baseline = stack.build(vor)                                   # before: it mutates the stack
pond = stack.excavate(bottom=bed, within="pond.gpkg").build(vor, attach=True)
ctx = mf.ModelContext(grid=vor, domain=pond.idomain)
print(pond.excavation_report())      # cut / rim / absorbed / connected layer and thickness
```

Pass `within` the same polygons you give `mf.lak(lakes=)`. The cut runs after
reconcile and the pinch policy; a cell left under the bed thinner than `min_remnant`
(0.1 ft) is absorbed into the one below rather than moving the bed. With a transparent
lakebed the connected cell's thickness sets the infiltration rate, so read the report's
`connected_min`/`connected_max` — a larger `min_remnant` changes the geology under the
pond, not just the solver. After the cut, "model top" inside the pond is the bed:
`irch="top_active"` recharge lands under it (MF6 zeroes it while a lake is wet there),
a `.gpkg` boundary on a dug cell is dropped with a warning, and a package declared
against a context built BEFORE the cut -- which still names the removed cells -- is
refused at build, naming the excavation: a boundary or observation on a dug cell, or
an `irch`/`ievt` array pointing at one (ledger 321; `irch=None` warns instead, with the
rate it loses, ledger 322). The excavated stack saves like any other (ledger 323):
`mf.disv(stack)` on the gridless stack is lazy and re-covers the pond's cells on
whatever mesh the model is built on, the footprint file and bed are carried by
`project.save(bundle_sources=True)` and hashed, so `runs.drift()` reports an edited
footprint, and a reopened run's frozen snapshot keeps the records.

`mf.lak` then reads the excavation's RECORD (ledger 325) -- `attach=True` puts it on the
grid, or pass `surfaces=pond` -- and a lake whose cells ARE the excavation's is written
by the excavated rule: its floor connects to the cell under the bed (MF6 counts from
that cell's top, which is now the bed) and its walls go into the NEIGHBOURS' cells, with
MF6's own face geometry. Leave `lake_bottom` out -- the lake takes the excavation's bed
-- and pass `unexcavated="raise"` so a record that never arrived is an error, not a lake
quietly written by the old rule:

```python
lak = mf.lak(context=ctx, nper=nper, lakes="pond.gpkg", lake_id_field="name",
             starting_stage={"pond": 526.0}, bed_leakance=3000.0, unexcavated="raise")
lak.metadata["lake_rules"]           # {'pond': 'excavated'}
```

Declare the model with NEWTON (`mf.gwf(..., newtonoptions="NEWTON")`): without it MF6
carries nothing through a wall whose cell is dry, and the build warns. It also warns
about walls into `icelltype` 0 cells, recharge on the pond's floor with no `rainfall=`
(MF6 switches it off while the pond is wet), and a bed so transparent that the mesh
sets the infiltration rate. A lake NOT excavated keeps the old rule and is warned about
by default: MF6 counts its floor from each cell's top, not from the bed.

An excavated pond couples tightly to the aquifer, and MF6's default lake solver pays for
it in outer iterations (+81% and 16 failed time steps on a real pond with its contour
table, +80% and a 2.16% worst lake-budget step without one). If your MF6 is a develop /
6.8 build, `mf.lak(..., implicit=True)` solves the stage in the flow matrix instead
(+11% with the table, +14% without, none failed, the lake budget closed; ledger 327
measured it end to end through this API) -- give the IMS `linear_acceleration="BICGSTAB"`
(CG is refused at build) and keep NEWTON on. An MF6 without the option is warned about
when the run starts (ledger 326).

**A stage-volume table is optional on an excavated lake.** Without one MF6 builds the
lake's volume from its cells' tops, and after the cut those tops ARE the bed: on the real
pond MF6's own curve came within 1% of the contour table from 530 to 533 ft (-3.1% at
527, +1.3% at 536), where the un-excavated geometry held half the storage (ledger 327).
Keep `LakeTableBuilder` for sub-cell bathymetry or a rim above the stages you reach.

**No stack?** A DISV built from a loaded run's arrays has no `LayerStack` to call
`.excavate` on. Put the arrays in a `LayerBuildResult` and cut that --
`dug = result.excavate(bottom=bed, within="pond.gpkg")` returns a NEW result -- then
build the DISV from `dug` and pass `ModelContext(grid=vor, domain=dug.idomain,
surfaces=dug)`, which is what carries the record to `mf.lak`.

`docs/reproducers/lak_vertical_floor.py` runs MF6 on the premise -- a VERTICAL connection
counts from its cell's TOP whatever bed is written for it -- and then the same nine-cell
bench through `stack.excavate` + `mf.lak`, where the floor carries MF6's formula from the
bed exactly.

#### Look at it before you trust it

```python
layers.plot.section(y=1000)          # filled, layer-coloured cross-section
layers.plot.map("clay")              # per-layer thickness
layers.plot.surface("all")           # every contact as 3-D height fields (plotly)
layers.plot.grid(scale=12)           # the layered cell VOLUME (needs `viz3d` extra)
```

`section` takes the line however you have it — a shapely geometry, a vector file
to read it from, the points themselves, or `x=`/`y=` for an axis-aligned slice.
With none of them you get a West–East line through the grid centre.

```python
layers.plot.section(line=LineString([(x0, y0), (x1, y1)]))
layers.plot.section(line="sections/a_a_prime.gpkg")
layers.plot.section(color_by="thickness")     # 'layer' (default) or 'thickness'
```

**To look at contacts individually, render them as VTK sheets:**

```python
layers.plot.surface("all", backend="vtk").show()             # one actor per contact
layers.plot.surface("clay", backend="vtk", scale=12).show()  # just this one
```

**In Jupyter, call `.show()` on a VTK scene** — unlike every other Picture, where
the bare expression is the idiomatic form. `.show()` hands you PyVista's trame
widget: interactive, server-backed, and nothing stored in the notebook file. The
bare expression instead emits a **self-contained 1–1.5 MB HTML blob**, which
JupyterLab commonly drops on the floor (the iopub data-rate limit is 1 MB/s by
default, and embedded vtk.js scripts can be sanitized) — and on the occasions it
*does* render, it writes those megabytes into your `.ipynb` for every scene.

For a file to keep or send, ask for one explicitly:

```python
layers.plot.surface("all", backend="vtk").html("surfaces.html")
```

That page opens by double-clicking it, not only over a web server. PyVista emits
its vtk.js bundle as `<script type="module">`, which a browser will not execute
from `file://` — the viewer shell loads, its loader never runs, and you get
vtk.js's *"Drop File / Explore Scene"* placeholder with
`ReferenceError: OfflineLocalView is not defined` in the console. `.html()`
rewrites that one marker on the way out; the bundle uses no module syntax, and a
test fails if a future PyVista ships one that does.

`grid(backend="vtk")` draws the cell VOLUME and `surface(backend="vtk")` draws
the CONTACTS as height fields — the same subject as their Plotly forms, drawn by
a different renderer, per the `backend=` rule. Cells with no data are dropped
rather than drawn at z=0, where they would read as a real contact at sea level.

##### Moving around in a 3-D scene

VTK and vtk.js do not share an interactor style, so the widget and the written
page answer to different gestures. **Shift + left-drag pans in both**; the
middle-drag habit works only in the widget:

| gesture | `.show()` widget | exported `.html` |
|---|---|---|
| left-drag | rotate | rotate |
| shift + left-drag | pan | pan |
| middle-drag | pan | *nothing* |
| right-drag | dolly | *nothing* |
| ctrl + left-drag | roll | roll |
| wheel | zoom | zoom |
| keys | *none* | `r` reset · `w` wireframe · `s` surface · `v` points |

The written page has no UI at all; click the canvas once to give it focus before
the keys respond. The widget has a toolbar at top-left (the dots button collapses
it) with reset, isometric, **+X/+Y/+Z views** and **parallel projection** — that
pair is how you get an orthographic slab view of the volume without writing code.
Both of those buttons disappear if you flip the toolbar to local rendering, and
parallel projection is *not* saved into an exported page: trame serializes camera
position, focal point, view-up and clipping range, and nothing else.

##### Showing, hiding, and reading the cells

Every layer is its own named actor, so a scene can be taken apart after it is
built:

```python
scene = layers.plot.grid()
scene.scene.actors["clay"].visibility = False   # hide one layer
scene.show()
```

Each cell carries `layer`, `thickness`, `top`, `botm` and `cellid` whichever one
is drawn, so `color_by` picks the picture without rebuilding the data:

```python
layers.plot.grid(color_by="thickness")   # viridis, per-cell thickness
layers.plot.grid(color_by="top")         # the layer's own top contact
layers.plot.grid(color_by="cellid")      # which cell is which
```

`top`/`botm` are the flat per-cell contacts; `elevation` is vertex z, so it ramps
*within* a cell. `cmap` defaults to whatever suits the chosen scalar.

**A written page keeps only the scalar you drew** — vtk.js serializes the active
array and drops the rest — and it has no picker, no keyboard and no widgets. For
a scene you want to interrogate rather than look at, write the mesh instead and
open it in ParaView, which gives per-block visibility, a cell inspector and
interactive slicing:

```python
layers.plot.grid().meshes[0].save("stack.vtu")
```

For per-cell numbers on a map instead of in 3-D, the choropleth hover already
reports **cell number, area, x, y, model top and every layer bottom** — but only
if the result was attached to the grid. Without that the elevation rows silently
vanish and nothing warns you:

```python
layers = stack.build(attach=True)     # or layers.attach_to_grid(vor) after the fact
layers.plot.map(basemap=True)
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

lake_bed = mf.Contours("lake_contours.gpkg", z="Elev", epsg="2927",
                       region_vector="lakes.gpkg")   # bathymetry, clipped to the lake
lak = mf.lak(context=ctx, nper=nper, lakes="lakes.gpkg", lake_id_field="name",
             starting_stage={"valley_lake": 101.0}, lake_bottom=lake_bed,
             bed_leakance=0.11, connection_modes="bathy", mover=True)

mvr = mf.mvr(nper=nper, moves=(
    mf.Move(mf.sfr_connection(sfr, "main_stem"),
            mf.lak_connection(lak, "valley_lake"), value=0.5),
))
```

Pull `sfr`/`lak` out as handles so the mover references them **semantically**
rather than by raw reach index. MVR is validated: moved packages must be declared
in the model *and* ordered before the mover.

This lake is declared on a stack whose top is the ground, so `mf.lak` warns that MF6
will count its floor from each cell's top. For a pond, cut it into the stack first
(§5, "Excavating a pond") and leave `lake_bottom` out.

#### A lakebed from bathymetry

`connection_modes="bathy"` (the default; `"automatic"` is its old name) connects
each lake cell through its bed *and* through every step in it, so `lake_bottom`
has to vary **per cell**. It is required — leave it out and `mf.lak` raises
`lake_bottom is required` — **unless the layer stack was excavated for the lake**
(§5, "Excavating a pond"): then the bed is the excavation's, declared once, and a
`lake_bottom` you do give must equal it. On a stack whose top is the ground, cut the
pond in first: MF6 counts each VERTICAL connection from its cell's TOP, so a bed
inside a cell does not lower it, and `mf.lak` warns about such a lake
(`unexcavated="warn"`; `"raise"` refuses it). Any of these works:

```python
lake_bottom=lake_bed                                # any Surface: mf.Contours, mf.Raster, algebra
lake_bottom="bathy.tif"                             # a raster path, Path or str
lake_bottom=pd.Series({812: 94.5, 813: 92.0})       # per cell, already computed
lake_bottom={"valley_lake": lake_bed, "pond": "pond.tif"}   # one per lake
lake_bottom=mf.Surface.offset_below(6)              # 6 ft below the model top
```

A `str` that names a column of `lakes` is still read as the column. A Surface or
raster is sampled onto the grid **once** — it used to be once per cell read, the
whole grid each time: 5 samples for a 2-cell pond. `units="meters"` on a Surface
converts to the model's units, and a relative Surface measures from the model top.

**A lake covers its cells by AREA** (ledger 241) -- a polygon drawn along cell
faces reaches only the cells inside it, and a sliver of a cell, which a CRS round
trip leaves, is dropped unless the whole lake is that thin (a detached part of a
multipart lake is judged on its own, so one narrower than 0.01 ft goes too). Every
lake cell gets its own vertical connection, so a spurious one moves the budget. A
lake covering no cell raises, with both extents and both CRSs.

**Only the lake's own cells need a value.** The neighbours each step is measured
against are land: where the source has nothing for one — a raster clipped to the
lake, a Series of lake cells — it sits at the model top. Two adjustments WARN: a
bed above its cell's model top is capped at the top, and a lake cell the raster
misses entirely (a sliver the polygon clips at its edge) gets the top too, zero
depth. A raster that misses the whole lake raises, and so does a Series missing a
lake cell.

A **flat**-bottomed basin is the other mode: `connection_modes="rectangular"`
with one `lake_bottom` and a `lake_top` per lake. A flat bed under `"bathy"` has
no steps, so on a lake that was NOT excavated it builds no sidewall connections and
warns; an excavated lake's walls go into its neighbours' cells whatever its bed.

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

**And for an array that comes from a POLYGON MAP, there is a fourth way that keeps
the recipe instead: `mf.Polygons`** (ledger 299).

```python
mf.npf(k=mf.Polygons("geology.gpkg", column="k", layer="units", default=1.0))
```

It records WHERE the values come from rather than what they resolved to, so the
declaration is JSON-clean, survives `save()`/`load()`, and resolves against whatever
grid the build supplies — one declaration builds at ncpl 36 and ncpl 100. `move_report`
reads it `rebuilt` where a plain array reads `unverifiable`, because an array names its
cells by POSITION and carries nothing to check against a new mesh.

Available on `mf.npf(k=, k22=, k33=)`, `mf.ic(strt=)` and `mf.sto(ss=, sy=)`.
`icelltype` / `iconvert` are REFUSED: they are integer flags, and every polygon rule
returns a float — there is no cell type 0.4.

**Name the rule if the default is not what you mean.** `rule="majority"` (default)
gives a cell the value of whichever polygon covers most of it and so invents nothing;
`rule="mean"` is the area-weighted mean, right for a continuous field and wrong for a
geology map — measured on a 64-cell grid with one contact it wrote 12 distinct
conductivities where the map holds 2. `rule="centroid"` matches `ZoneSpec.from_polygons`.

**A cell no polygon reaches is REFUSED unless you pass `default=`.** An unfilled cell
is a NaN, and MF6 answers a NaN budget while reporting normal termination.

For a BOUNDARY condition — a conductance, a stage, an invert — this is the wrong tool:
`mf.drn.gpkg(...)` and its six siblings already read a polygon's attribute columns and
have been lazy since ledger 279. `GeoPackageSource.k_array` still works and is not
deprecated; it is the EAGER form, and it is what cannot be saved.

Call **`project.validate()`** to see which case you are in; it lists the offending
simulations without touching disk. Since 2026-09-04 `save()` no longer raises on
them — it writes the project document, every package and every grid, SKIPS the
simulation it cannot serialize, and names it in the returned `SaveReport.skipped`
(a dangling `mf.ref` still raises, but only after the write). Pass `strict=True`
for the old all-or-nothing behaviour. Note `validate()` inspects *simulations*, so
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
model.ghb.results.q.get()       # a DataFrame

model.plot.map(layer=0, contours=True)   # six verbs, four scopes
model.plot.section(cells=[12, 40, 88])   # cross-section: cells by layer, water level on top
model.plot.grid()
model.plot.animate([model.plot.map(layer=0, per=p) for p in range(model.nper)])
```

Read grammar is always `model.<pkg>.<inputs|results>.<noun>.<verb>()`.
Everything drawn returns a `Picture`: `.fig` / `.show()` / `.save(p)` / `.html(p)`,
never a trailing `.plot()` -- except a `backend="mpl"` figure, which is `section`'s
default (`.savefig(p)`; `backend="plotly"` for the Picture, `fill=None` for the head
line profile).

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
   only — `mf.sfr`/`mf.lak` resolve cells at declaration; `mf.uzf` and the `.gpkg`
   helpers do too WITH `context=`, and resolve at build without it.
6. **`attach=True`** on `stack.build()` whenever SFR or LAK is in the model.
7. **Read `stack.qc()` before `build()`.** A silently pinched-out layer becomes an
   inactive region you'll spend an afternoon chasing in the head field.
8. `GridSpec.structured` / `from_geopackage` exist but **fail fast** — only
   `voronoi` and `python` resolve. Use `GridSpec.from_object` for a built grid.
9. **`project.save()` skips a simulation that embeds numpy** — arrays *and*
   scalars — but persists everything else and names what it left out
   (`report.skipped`). It used to refuse the whole project, so one such
   simulation cost you the entire package library. A package carrying arrays is
   fine either way: it is pickled rather than JSON-encoded. Call
   `project.validate()` to see the issues without touching disk. A **dangling
   `mf.ref`** still raises — that is a mistake with a fix — but only after the
   write, so the library persists regardless.
10. **Mixed source resolutions need no pre-alignment.** Area-weighted sampling
    puts a 25-ft DEM and a coarse contoured contact on the same footing. A
    mismatched **CRS**, though, samples in the wrong place and surfaces as a
    nonsense `qc()` thickness rather than an error — check that yourself.
11. **A lake on a stack whose top is the ground exchanges nothing through most of
    its floor.** MF6 counts a VERTICAL lake connection from its cell's TOP, not
    from the bed written for it, so a bed inside a cell leaks nothing until the
    stage passes that top — 411 of 612 on a real pond carried exactly 0 all run
    (ledger 327).
    Cut the pond in with `stack.excavate(bottom=, within=)` before `mf.lak`, leave
    `lake_bottom` out, and pass `unexcavated="raise"` so a record that never
    reached the lake is an error.
