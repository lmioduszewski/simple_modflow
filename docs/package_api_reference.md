# Package API reference

The vocabulary of the **package-first API** — the canonical way to build and read
MODFLOW 6 models in myflopy. Two halves:

- **Build side** — `mf.<pkg>(...)` helpers that assemble a model (`package_api.py`).
- **Read side** — the `model.<pkg>.<inputs|results>.<noun>.<verb>()`
  view grammar for tables, plots, and maps.
- **Drawing** — the six plotting verbs, free in `myflopy.plot` and bound as
  `model.plot` / `vor.plot` / `stack.plot` (section C).

> **Authoritative sources (this doc is the human-readable overview of them):**
> - Build-side signatures are pinned in `tests/api_snapshot.json` (regenerate with
>   `python scripts/derive_api_snapshot.py`; the diff is the public-API change log).
> - Per-package fields/results/capabilities live in
>   `src/myflopy/modflow/mf6/package_registry.py` (the single per-package descriptor).
> - The read-side grammar is normative in `docs/view_layer_conventions.md`.
> - Deliberate gaps/deviations are in `docs/compromises_and_deferrals.md`.
>
> When these disagree with this file, they win — update this file to match.

> **Starting a new model?** `docs/model_building_cheatsheet.md` is the ordered
> walkthrough — contours → surfaces → layer stack → packages → run → results —
> with every call verified end to end. Runnable:
> `examples/mf6/contours_to_model.py`. This file is the per-symbol reference; the
> cheat sheet is the path through it.

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
| `mf.simulation(*models, packages=[tdis, ims])` / `mf.SimulationSpec(...)` | one MF6 simulation; `packages=` carries timing and solvers, defaults fill the rest |

### Discretization, core, solver

`mf.disv`, `mf.dis`, `mf.disu`, `mf.ic`, `mf.npf`, `mf.sto`, `mf.oc`, `mf.tdis`,
`mf.ims`, `mf.ems` — each `mf.<name>(...)` returns a `PackageSpec`.

`mf.ic` / `mf.oc` / `mf.disv` / `mf.dis` / `mf.disu` are **model-kind-aware**: they
build the GWF/GWT/GWE (and, for `oc`/`disv`/`dis`, PRT) FloPy class off the model's
type, so the same helper serves a transport, energy, or particle-tracking model
(for GWT/GWE `mf.oc`, pass `concentration_filerecord=`/`temperature_filerecord=`
via `**options`; for PRT, `trackcsv_filerecord=`). `mf.npf` and `mf.sto` are
GWF-only (MF6 has no GWT/GWE variant); `mf.ic` and `mf.disu` reject PRT (MF6 has
no `ModflowPrtic`/`ModflowPrtdisu`).

**Solvers**: `mf.ims(...)` for GWF/GWT/GWE; **`mf.ems(...)`** for PRT — PRT models
are *explicit* in MF6 and are rejected under IMS6 ("Explicit models require
EMS6"). Everything simulation-level goes through ONE list,
`mf.simulation(*models, packages=[tdis, ims, ims2])`, the way model packages go
through `mf.gwf(packages=[...])`; myflopy fills in whatever the list omits — a
steady `mf.tdis()` and a solver for every model no supplied solver claims, one
IMS per GWF/GWT/GWE and one EMS per PRT. **`models=` is optional on a solver**:
one that declares none claims every model the others left, so a single-model
simulation writes `packages=[mf.ims(complexity="COMPLEX")]` and names nothing.
Solvers are ordered by MODEL DECLARATION order however you list them — MF6
requires a GWF solution in `mfsim.nam` before the solution of any model it feeds.
A **model** package handed to `packages=` raises (ledger 200).

**Exchanges**: `mf.exchange(flow, transport)` reads the two models' kinds and
picks the MF6 class; `mf.exchange.gwfgwt/.gwfgwe/.gwfprt/.gwfgwf(a, b)` say it
outright and accept plain names. **Nothing wires one for you** -- a coupled
simulation without an exchange builds, runs, and moves nothing, under normal
termination. Order does not matter; the GWF model is put first because that is
what MF6 wants. `gwfgwf` is the one carrying real content and refuses a call
without `nexg`/`exchangedata`, which MF6 rejects as *Required block "DIMENSIONS"
not found*. The four `build_*_exchange` functions underneath stay engine-tier.

**Grid helpers** (`mf.disv` is the default — myflopy is Voronoi/DISV-first):

- `mf.disv(nlay=, ncpl=, nvert=, vertices=, cell2d=, top=, botm=, idomain=)` —
  vertex/unstructured grid; pull the mesh from a `VoronoiGridPlus` /
  `vor.get_disv_gridprops()` and the layer arrays from a `LayerStack`.
- `mf.dis(nlay=, nrow=, ncol=, delr=, delc=, top=, botm=, idomain=)` — structured
  (row/column) grid; a passthrough for rectilinear or externally supplied models.
- `mf.disu(nodes=, nja=, top=, bot=, area=, iac=, ja=, idomain=)` — fully
  unstructured grid with connectivity you supply yourself (`ihc`/`cl12`/`hwva`,
  `nvert`/`vertices`/`cell2d` via `**options`); prefer `mf.disv` unless ingesting an
  existing DISU mesh. A `VoronoiGridPlus`'s `iac`/`ja`/`cl12`/`hwva` are MF6's own DISV
  geometry -- `cl12` PER SIDE, each cell's DISV centre to the face -- so a DISU written
  from them reproduces the DISV (ledger 328).

`mf.dis`/`mf.disu` build and run, but the choropleth-map / cross-section / animation
viz targets **DISV/Voronoi** meshes and is not wired for structured or raw-DISU
grids (use FloPy's own plotting, or DISV).

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

> **`gaps=` and `snap_within=` on every `.gpkg` helper** (ledger 288). `gaps="warn"` (default) names a feature that reached no cell and says what was not written; `"nearest"` moves a POINT that only just missed, within `snap_within` (default: one cell width, measured on the cell it would snap to); `"raise"` refuses; `"ignore"` is silent. A source that lands NOWHERE is refused in every mode. An areal feature is clipped to the part inside the grid, never moved, and the clipped-away fraction is reported.


- `()` — direct records: `mf.riv(stress_period_data={0: [[(0, 3), 100.0, 5.0, 98.0]]})`.
- `.gpkg` — map GeoPackage features onto cells; field args accept a column name or a
  `CellSurfaceOffset` (surface-relative), e.g. `mf.drn.gpkg("seeps.gpkg", context=ctx,
  nper=12, elevation=mf.CellSurfaceOffset("cell_top", offset=-2.0))`. A source whose
  features land on no active cell RAISES, naming both extents and both CRSs (the file
  may have none); a feature that misses while others land is WARNED by name (ledger 211).
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

**A lakebed from bathymetry.** `mf.lak(connection_modes="bathy")` (the default)
needs `lake_bottom` per cell: a `Surface` (`mf.Contours(...)`, `mf.Raster(...)`,
or algebra on one), a raster path (`Path` or `str`), or a cell-indexed
`pd.Series` — bare for every lake, or `{lake_id: ...}` for each. A Surface or
raster is sampled onto the grid once. Only the lake's cells need a value: the
neighbours its sidewall steps are measured against sit at the model top where the
source has none, and a bed above its cell's top is capped there, with a warning.
**A pond cut into the layer stack takes its bed from the cut** — `lake_bottom` is
then optional and, given, must equal it (below); on a stack whose top is the ground,
excavate first.

```python
bed = mf.Contours("pond_contours.gpkg", z="Elev", epsg="2927", region_vector="pond.gpkg")
lak = mf.lak(context=ctx, nper=nper, lakes="pond.gpkg", lake_id_field="name",
             starting_stage={"pond": 101.0}, lake_bottom=bed, bed_leakance=0.1)
```

A flat-bottomed basin is `connection_modes="rectangular"`, with one `lake_bottom`
and a `lake_top` per lake.

**Which layer of a GeoPackage.** `mf.lak(lakes="site.gpkg", layer="ponds")` reads that
layer, by NAME; omitted, the first, which was all `mf.lak` could read before. A name
the file does not hold raises, listing the layers present, and an integer is refused
-- beside `only_layer` it reads as a model layer (`geopandas.list_layers(path)` gives
the names). `stack.excavate(within=, layer=)` takes the same name. The
lakes are read by `grid.geometry.read_footprints`, their cells chosen by
`footprint_cells` and a Surface/raster bed sampled by `surfaces.sample_bed` -- the
shared reader, cell rule and bed sampler an excavation of the layer stack
(`stack.excavate`) is built on.

**What the bed does in MF6.** MF6 ignores the elevations written on a VERTICAL
lake connection and uses the connected cell's TOP instead, so a bed lying inside a
cell, below its top, does not lower that connection: it picks the layer each
VERTICAL connection lands in. It is also the bottom of each HORIZONTAL face
`"bathy"` writes on the lake cell itself, and MF6 does read those, exchanging water
through them from the bed up. Measured: a VERTICAL connection under a cell topped
at 8.0 with the bed written as 3.0, 6.0 or 0.0 held the same 400.0 ft3 at stage
12.0; held at 6.0, below the top, it carried nothing and the lake held 0.0 ft3,
while a HORIZONTAL face on the same cell from the bed (3.0) up carried 3.46 ft3/d
out of the lake in a steady period. In a transient one MF6 lets no more out of a
lake than it holds plus what arrives, so the same face carried 0.0 (2.0 with an
inflow of 2.0), while water from the aquifer still came IN through it.
**An excavated lake is the way out** (ledger 325). When the lake's cells ARE an
excavation's (`stack.excavate(within=<the lakes= polygons>)`, its record reaching
`mf.lak` through `attach=True` or `surfaces=<the result>`), `mf.lak` writes it by the
EXCAVATED rule in either mode: one VERTICAL connection per lake cell to the top active
cell, whose top IS the bed, and HORIZONTAL walls INTO each neighbour whose top is higher,
`CONNLEN` from the neighbour's DISV centre to the shared face and `CONNWIDTH` the face --
MF6's own geometry. `lake_bottom` defaults from the record; a partial match raises. A
pit lake in uniform flow carries 1.005 of the through-flow of the same pit as a
K = 1e4 inclusion, against 0.52 for `"bathy"` on the un-excavated stack (0.39 while its
`CONNLEN` was the centre-to-centre spacing, before ledger 328). Every other lake
keeps its rule and, by default, is warned about (`unexcavated="warn"`/`"raise"`/
`"ignore"`); `metadata["lake_rules"]` says which rule each lake got. A stage-volume
table is optional on an excavated lake: MF6's own curve, from cell tops that now ARE the
bed, came within 1% of a contour table from 530 to 533 ft on a real pond (ledger 327).
`docs/reproducers/lak_vertical_floor.py` runs MF6 on the premise (a VERTICAL connection
counts from its cell's top, whatever bed is written) and on the same bench excavated.

**Solving the stage in the matrix: `mf.lak(implicit=True)`** (ledger 326). MF6's LAK
`IMPLICIT` makes each lake stage one more unknown of the flow matrix, solved WITH the
heads instead of by turns -- far fewer outer iterations for a strongly coupled lake, the
same converged answer (a transient pit matched the default formulation's stages to
1e-8 ft). It needs a BICGSTAB solution (a CG one is refused at build; the IMS
`mf.simulation` makes picks BICGSTAB -- for a library lake reached through `mf.ref` too,
once the build resolves it), every lake connected, and ideally NEWTON (warned).
It exists in MF6's develop branch / 6.8 only and FloPy 3.10 cannot write it, so
myflopy writes it into the LAK file itself (`lak_implicit.ModflowGwflakImplicit`); the
executable is probed once at run time and one without the option is warned about. A
run reopened with `Run.open` or `mf.load_mf6_run` keeps it, and so do
`ParallelModelWorkflow.split_model` and a package artifact; a bare
`flopy.mf6.MFSimulation.load` drops it silently. `model.config` shows it as the lake's
`implicit` setting, so `diff().config` tells two runs apart by it.
The `mf.lak` docstring (`help(mf.lak)`) spells out what each mode writes, and every
parameter; a `LakeOutlet`'s `rate` is negative to remove water and positive to add
it, and per-period values are zero-based. What it used to accept without checking
(ledger 318) is refused at the call now: a period past the last, a `status` MF6
does not know, a negative forcing, `nlakes`/`pname` in `**options`, a rectangular
bed no layer holds; an array or Series is read by period.

The same bathymetry can also give the lake a stage/volume/area table:
`tables={"pond": mf.LakeTableBuilder(dem=bed.resolve_source(), footprint=Path("pond.gpkg"), stages=...)}`.
`dem=` takes the raster PATH that `resolve_source()` returns, not the Surface.
`footprint=` is any form `mf.lak(lakes=)` reads -- a path (`footprint_layer=` picks
the layer), a frame, an `mf.Area`, a geometry -- reprojected to the raster's CRS,
and refused when it reaches no pixel with a bed. Start `stages` at or below the
lowest bed: extra zero-volume rows are reduced to one, the same curve to MF6. The
class docstring lists every field, its mode and its errors; the prism
(`area`/`top`/`bottom`) and explicit-`rows` forms (a DataFrame too) are there.

**Which cells a lake covers.** By AREA, through the shared rule (ledger 241 -- the
same `cell_shares` behind `.gpkg` boundaries and head targets): a lake polygon
drawn along cell faces reaches only the cells inside it, and an overlap that is a
sliver of a cell -- what a CRS round trip leaves -- is noise and is dropped,
unless every piece of the lake is that thin. A MULTIPART lake is measured per
part, so a detached part narrower than the floor (0.01 ft) goes with the slivers
while a 1 ft pond stays. Each lake cell becomes a full LAK cell with its own
vertical connection, so a spurious one changes the lake's sidewall topology and
its budget. A lake covering no cell RAISES, naming why and giving both extents and
both CRSs (ledger 211) -- from `validate()`, so reading `LAKBuilder.lake_cells` to
see WHICH lake missed still answers `{"pond": []}`. (A lake polygon with no area
at all -- a collapsed outline -- is refused earlier, on the first read of
`lake_table`.)

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
hover** are still Phase 6. The **results tier** (`model.conc`/`model.temp`) is
**built** as of §6.0/6.1/6.2 (2026-07-24) — see the read-side model-level reads
below. Transport budget **tables** read correctly as of 2026-07-27 — `model.bud()`
and the package result tables handle both MF6 record shapes, agree on zero-based
node ids, and label hovers `M/T` (GWT) / `E/T` (GWE) rather than `ft³/d`
(ledger 90–92, 95). The `model.budget.<term>` noun shipped 2026-07-27 (§6.1/6.2
item 3) — see "Model-level reads" below. The grouped transport fields
(`group.conc`/`group.temp`, `group.diff().conc`/`.temp`) shipped the same day
(§6.1/6.2 item 4), which also fixed `section` on `model.conc`/`model.temp` — it had
raised since those readers shipped (ledger 99). **Transport calibration** ships as of
2026-07-28: `ConcTargets` observations go straight into `cal.observe(...)`, and the
canonical model has a GWT sibling (`build_canonical_transport_calibration_demo`) whose
worked example is `examples/mf6/notebooks/canonical_07_transport_and_prt.ipynb`.
Zones for `style="zone"` come from `mf.ZoneSpec` (`from_raster` — majority vote, not mean — `from_polygons` — by cell **centroid**, and a centroid on a zone's edge belongs to the first zone covering it — `from_array`); the shape each target family needs is resolved for you.

Transport properties are `parameterize` targets too: `cal.parameterize("porosity")`
(aliases `mst.porosity`, `n`) resolves the GWT sibling automatically. Calibrate it
alongside heads — porosity is absent from the flow equation, and `v = Ki/n` makes it
near-collinear with K from concentration alone. Still deferred: `TempTargets` (GWE
observations) and `dsp.alh` (ledger 103).

### PRT particle-tracking packages

A PRT model is fully declarable spec-first (§6.3A, 2026-07-24) — a coupled
GWF+PRT simulation needs **no FMI package** (flows pass through the GWF-PRT
exchange) and no raw PackageSpecs:

| factory | package | notes |
|---|---|---|
| `mf.mip` | model input (PRT) | `porosity=` required; optional `retfactor=`/`izone=` |
| `mf.prp` | particle release points (PRT) | `packagedata=` rows `(irpt, (layer, cell), x, y, z[, boundname])`; `nreleasepts`/`perioddata`/`boundnames` defaulted |
| `mf.ems` | explicit solver | `models=` — required for PRT (IMS6 rejects explicit models) |

```python
particles = mf.prt("particles", packages=[
    mf.disv(**grid), mf.mip(porosity=0.25),
    mf.prp(packagedata=[(0, (0, 0), x, y, z, "west_wells")]),   # boundname → release group
    mf.oc(trackcsv_filerecord="particles.trk.csv",
          budget_filerecord="particles.bud", saverecord=[("BUDGET", "ALL")]),
])
sim = mf.simulation(flow, particles,        # solver default: IMS for flow, EMS for particles
    exchanges=[mf.exchange(flow, particles)])   # -> ModflowGwfprt, chosen from the kinds
```

Release-point `boundname`s are echoed (uppercased) into the track CSV's `name`
column — the release-group key the capture map reads. The high-level
`model.particle_tracking` / `PRTProject` runtime remains the sanctioned path for
*post-hoc* tracking over an already-run flow model (a separate simulation, where
FMI + grid copying ARE needed); on that path the groups come from
`PRTReleasePoints.from_cells(..., group="west_wells")` (or `from_points`), and
`PRTReleasePoints.merge(west, east)` combines grouped sets into the one PRP MF6
wants, renumbering `irpt`. Merging a **named** set with an **un-named** one is
allowed and turns BOUNDNAMES on for the whole package, so MF6 synthesizes a label
for the points that had none (`PRP000000002`) — it appears as an ordinary group in
`capture` and on the pathline legend. Name every set if you do not want that.
A release with no groups at all leaves `name` blank throughout, and the views fall
back to per-particle coloring.

**`from_points` puts each point in the cell ACTIVE in its release layer** (ledger
242, 2026-09-12). It located points with FloPy's `modelgrid.intersect`, which
answers with the lowest CELL2D and reads no idomain: with cells 0 and 4 inactive
in layer 0, a point on the 4|5 face resolved to cell 4 and one on the 0/1/4/5
vertex to cell 0 — and FloPy's PRP then refused the whole package, `Cellid (0, 4)
is outside of the active model grid`, naming an id the caller never wrote (the
active neighbour the point also touches releases and tracks). It now uses the
shared coverage rule — the same one a head target and a `.gpkg` feature resolve
by — so a point on a face or vertex takes the lowest-numbered cell active in
`layer=`, a frame in another CRS is reprojected to the grid's first, a point
landing on no active cell is left out and named in a warning (`irpt` renumbered
over the survivors, as `merge` does), and none of them landing raises with both
extents and both CRSs (ledger 211) where FloPy named only the point. It is also
much faster, since one spatial-index query replaces FloPy's per-point scan of
every cell: 200 points on a 10,000-cell grid, 7.5 s → 0.016 s. The coordinates
WRITTEN are still the model's own frame, which on an offset grid is not the
world's (ledger 218). `from_cells` trusts the cell ids it is given — it neither
range-checks them nor reads idomain. A **structured (DIS) grid keeps FloPy's
`intersect`, and with it FloPy's own refusal**: the shared rule is indexed by cell
id and a DIS cell is `(row, col)`, so neither half of this — not the activity, not
the ledger-211 miss message — reaches a structured model, where one off-grid point
still raises `point given is outside of the model area` and ends the call.

**Reading a finished PRT run (§6.3B/§6.3C, 2026-07-25).** Trajectories are not a
per-cell field, so `PRTRunResults` exposes four view nouns that answer the normal
verbs (`get`/`summary`/`plot`/`map`/`mosaic`, `prt_maps.py`) — one keeping the
paths, three rolling them onto cells:

| noun | `get()` rows | `map()` draws |
|---|---|---|
| `results.pathlines` | the **normalized track records**: every raw track-CSV column plus `cell, layer, travel_time, release_group, particle` | one map polyline per particle over a base map (`base="heads"`/`None`/a `Choro`), colored per release group |
| `results.travel_time` | one per `(layer, cell)`: `layer, cell, travel_time, particle_count, min_time, max_time, release_groups` | time-of-travel choropleth (`stat="median"`, `'earth'`, optional `logscale`) |
| `results.endpoints` | the same columns (the count is the mapped one) | particle-termination counts per cell |
| `results.capture` | those columns prefixed by `group` (one per `(group, layer, cell)`) | **mosaic**, one panel per release group on a shared scale; `group=` gives one `Choro` |

```python
results = model.particle_tracking.prt(workspace=ws, release_points=rp).run()
results.pathlines.map()                                  # tracks over the water table
results.travel_time.map(stat="median", logscale=True)   # time-of-travel figure
results.capture.map()                                    # capture zones, side by side
results.travel_time.plot()                               # cumulative arrival curve
```

`pathlines` is a **view, not a frame**: `results.pathlines.get()` is the record
table and `results.track_records` the untouched CSV (what FloPy/PyVista consume).
Its `map()` returns the `Choro` carrying one `Scattermap` line per particle —
hovering a vertex reports that particle's cell, layer, elevation, release
group, and elapsed time — `plot()` draws elevation against travel time, `mosaic()` gives one panel
per release group, and `backend="mpl"` returns the existing FloPy plan view.
`max_particles=` (default 250) caps a large run by sampling **stratified across
release groups** (so a cap at or above the group count cannot drop a whole capture
zone), and says so in the title and a warning. Selectors the call cannot honor
raise rather than being dropped: base-map options on `backend="mpl"`, `per`/`layer`
on an already-built `base`, and a `base` shared across `mosaic()` panels. Group colors come from
`viz.category_colors` / `PALETTE.categorical`, so a group reads the same on the
pathline map, the arrival curve, and the capture bars.

The three per-cell maps are **time-integrated** over the whole run: they reject `per=` (and
`animate()`) rather than accept an axis they would ignore, and carry no period
footer (`result_hover(..., footer=())`). `layer=None` pools every layer into one
plan-view panel, recomputing the statistic over the pooled particles rather than
averaging per-layer values. Cells no particle reached stay blank rather than
reading as a zero travel time. `plot()` is a distribution across particles
(cumulative arrivals / count bars), not a stress-period series. `capture` needs
release-group boundnames; without them it names the fix and offers
`by="release_point"`.

### Geometry, layers, context

- `model.packages.summary()` — one row per package: its `package` name and `type`
  (two packages of one type are two rows), `kind`, `mappable`, the input
  `fields` and `results` you can draw. Cheap by default (registry + package list,
  no package data read); `summary(detail="data")` adds `records`, `periods` and
  `layers` — the last of which tells you a package occupies no layer 0 *before*
  you draw an empty map of it.
- `model.packages.mosaic()` — every mappable package drawn as one grid;
  `packages=[...]` to restrict it.

- `mf.ModelContext(grid=, domain=, surfaces=, dates=)` — attaches to the **model**
  (`mf.gwf(..., context=ctx)`); carries the grid/idomain/surfaces the GIS-aware helpers use.
- `mf.LayerStack(vor, top=...).add(...).build()` — declarative layer top/botm/idomain
  (compiles to `mf.LayerSurfaces` / `mf.Surface`). `vor` is optional: the layering
  can be declared before any grid exists and bound at `build(vor)`.
  `.add(name, ..., split=3)` cuts one geologic **unit** into N model layers
  (`name_1..name_N`, or `names=[...]` for your own, one per layer) without moving
  its base; the result carries
  `units` (unit → layer indices) and `per_layer({unit: value})` for the positional
  per-layer lists `mf.npf`/`sto`/`ic` take.
- `stack.excavate(bottom, within, *, name=, layer=, min_remnant=0.1)` — a void cut
  into the stack: a pond, basin, pit or trench (ledger 320). Inside `within` — the
  SAME polygons you pass to `mf.lak(lakes=)`, read by the same reader and covering
  the same cells — every surface above `bottom` is lowered to it, so the first active
  cell under the bed has its top exactly at the bed and the cells above go idomain 0.
  It runs in `build()` AFTER reconcile and the pinch policy, and absorbs rather than
  deepens: a cell under the bed thinner than `min_remnant` (or inactive) is folded
  into the cell below, and the bed never moves. `bottom` is an absolute elevation (a
  Surface, a number, a raster path, or a column of the footprint file); a surface
  measured from the one above is refused (write `ground - 9`), and so is algebra over
  a `fill="propagate"` source (pass it bare: it takes the ground where it has no data,
  and those cells are not cut). An excavation's name may not be a unit's, nor a layer
  a unit builds later (`<unit>_i` of a `split`). Footprints of two
  excavations may not share a cell. Mutates and returns the stack, like `add` — build
  the un-excavated baseline first. The result gains `excavations` (name →
  `mf.ExcavationRecord`: cells + sampled bed), `ground` (the pre-cut top), `excavated`
  and `absorbed`; `result.excavation_report()` says what was removed, absorbed and
  connected; `build(attach=True)` publishes the records as `vor.excavations` beside
  `gdf_topbtm`. `result.excavate(...)` cuts a BUILT result (a model with no stack).
  `stack.to_disv(vor)` is now `stack.build(vor).to_disv()`. A package declared
  against a context from BEFORE the cut is refused at build, naming the excavation:
  a list boundary or an observation on a dug cell, an array recharge/ET index
  pointing at one (ledger 321). An excavated stack SAVES (ledger 323): its excavations
  ride in the `$layer_stack` payload — the footprint file by path (under `"zone"`, so
  `save(bundle_sources=True)` carries it, a lazy `mf.disv` digests its content and
  `runs.drift()` reports it edited), an `mf.Area` or geometry inline as WKT + CRS —
  so `mf.disv(<gridless excavated stack>)` is lazy and re-covers the cells on whatever
  mesh it is built on, `with_grid` included. A payload is refused at load wherever
  `excavate()` would refuse the call (a repeated or empty name, a relative bed, `layer=`
  beside an inline area); only the footprint file is left unread until build. A run's
  frozen snapshot keeps the records, so a reopened run still names its excavation.
  `mf.lak` reads the records (ledger 325): a lake whose cells ARE an excavation's is
  written by the excavated rule -- floor at the bed, walls into the neighbours' cells --
  and `lake_bottom` defaults from the record.
- `mf.Toward(target, fraction)` — the contact a share of the way from the surface
  above down to `target`; the primitive `split=` is built on.
- `mf.CellSurfaceOffset("model_top"|"cell_top"|"cell_bottom", offset=, minimum=)`.
- `GridSpec.voronoi(boundary=, refinement=, breaklines=, points=, ...)` — refinement
  arrives on **three asymmetric channels**: `refinement=` takes ONE polygon source
  (every feature = one region, sized by `fields={"area": "<column>"}` or
  `refinement_max_area`), `breaklines=` takes a LIST of line sources buffered into
  corridors (`breakline_buffer`, default 10; priority **1**, vs 0 for refinement),
  `points=` takes a LIST of point sources pinned as vertices. Mixing polygons and
  lines in one `refinement=` layer raises. Every sizing option is a **named
  parameter** (`boundary_*`/`refinement_*`/`breakline_*`), and an unknown keyword
  raises with the nearest match rather than being silently ignored. Sources must be
  file-backed specs — a GeoDataFrame, a bare path, or a list is refused by name.
  The cells are repaired to tile the boundary (ledger 316); `tessellation="flopy"`
  keeps FloPy's `tri2vor` cells as they were.
  Full table + traps in `docs/model_building_cheatsheet.md`.
- `GridSpec.resolve(project_root=, workspace=, build=, return_triangle=)` — turns a
  grid *recipe* into a grid *object*. `project_root` anchors relative source paths;
  `workspace` is where Triangle's `_triangle.*` files land (one directory per grid);
  `build=False` returns the prepared-but-unmeshed `TriangleGrid` for inspection;
  `return_triangle=True` returns `(vor, tri)`. Full reference in the docstring.
  Note the CRS falls back to `"EPSG:2927"` when the spec names none.
- Grids: `mf.VoronoiGridPlus`, `mf.TriangleGrid`, `mf.GridSpec`; sources:
  `mf.GeoPackageSource`, `mf.read_gpkg`, `mf.read_shp_gpkg`.

### Specs & lifecycle

`mf.SimulationSpec`, `mf.ModelSpec`, `mf.PackageSpec`, `mf.ExchangeSpec`;
`mf.Project(root, name=)` → `prepare_run`/`run` → `mf.Run`; `mf.load_run` /
`mf.load_mf6_run`; `mf.ModelGroup`; `mf.ParallelModelWorkflow`.

**The project's libraries.** `project.packages` / `.grids` / `.simulations` are
`dict` **subclasses** — every dict operation works unchanged, and each also
answers the same questions. Full walkthrough in
[project_workflow.md](project_workflow.md).

| | |
|---|---|
| what is defined | `project.packages.names()`, `.summary()`, `repr(project)`, `project.summary()` |
| give me that one | `project.package(key)` / `.grid(key)` / `.simulation(name)`, or subscript (a miss names the near miss) |
| make a variant | `project.packages.fork(source, key, **overrides)`, `project.simulations.fork(source, name, swap={...})` |
| change one in place | `project.packages.edit(key, **overrides)` |
| how it would persist | `project.packages.persistence()` |

Forks and edits are **built, not patched** — they go through the same
constructor `mf.npf(...)` does, so a typo raises there rather than at MF6 build
time, and they deep-copy by default so editing a fork cannot write through into
the original's arrays.

**Runs.** `project.runs` is the ledger: one row per run under the project,
including bare MF6 workspaces with no manifest. `run.model_options` returns the MODEL-level options (`NEWTON`, `print_flows`, ...) as `mf.gwf(...)` keyword arguments, so rebuilding specs from a workspace does not silently drop them (ledger 286). `run.outcome` is what **MF6**
reported in `mfsim.lst` (`normal`/`failed`/`incomplete`/`absent`/`superseded`),
which is a different axis from `run.status`, what **myflopy** did;
`run.reconcile()` adopts the verdict. `project.runs.drift()` reports which
library packages a run used have changed since it ran. `project.save()` returns a
`SaveReport`, not a path.

**Serialization.** Any spec round-trips to a JSON-safe dict (`to_dict()`/
`from_dict()`) — builders as importable references (including the list-BC
`functools.partial` builders), sources tagged, paths POSIX-encoded, **tuples
tagged** (`{"$tuple": [...]}`, so a cellid stays a cellid — `observations=`
writes its records with `repr` and emitted `[0, 0, 0]` without it), and
**integer stress-period keys restored on read** (`to_dict()` must stringify them
to reach JSON; `from_dict()` turns `'0'` back into `0` for package
options/inputs, so a saved list BC keeps its rows — ledger 179). A whole
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

A target given as a POINT resolves to the one cell holding it; on a shared face
or vertex, the lowest-numbered cell active in that target's own layer. Targets
are reprojected to the grid's CRS first, a target on no active cell is left out
and named in a warning -- a Python `UserWarning` as well as the log, since a log line
alone reaches nobody by default (ledger 342, ledger 288's rule) -- and none of them
landing raises. A target given as a
cell id is checked instead: an id the grid does not have, or one inactive in its
layer, raises (ledger 233). A name used in TWO layers is two observations,
`<name>_lay<layer>`: `targets.observation_labels(model)` is the one definition,
and the MF6 records, the simulated CSV, the forward-run mapping and the pst all
carry it (ledger 234).

A DRN seepage ZONE reaches the cells it COVERS -- a polygon by area, a trace by
the cells it runs through, a point the one cell holding it -- active in the
zone's own layer, and comes back in the grid's CRS whatever the zone file's was.
It used to take every cell it touched, so a zone drawn on one cell claimed nine
and two zones abutting on a column of faces both claimed that column, summing
every drain in it into two observations; a zone off the grid came back empty and
its flow as zero, and now raises (ledger 238). A zone given as several rows is
one zone: its cells are the union -- on the observation route as on the budget
file, which until ledger 336 summed only the LAST row's cells there -- and it is
reported missing only when every part of it misses. Zone membership is geometry
only: a covered cell that carries no drain is still in the zone and adds zero.

**A zone's cell is a `(layer, cell)` pair, on every route** (ledger 342). MF6's
cell-budget file names a drain by its model-wide USER node -- `layer * ncpl + cell` on
a DISV, whatever the idomain removes -- and the budget route matched that against a
zone's in-layer cell numbers, so a zone in layer 1 read whichever LAYER-0 drain shared
a number (-6.77 for a zone whose drains held -5975.73); it now matches the pair through
`model.node_to_lni`, and the zone snapshot a PEST build writes carries `layer`, so the
forward run reads the same drains. A zone's rows may sit in different layers; it is one
zone of their union. A row gives explicit `cells` or is resolved from its geometry, and
one table may mix the two (a missing `cells` entry means "resolve me"); a dict or a list
of record dicts with a `geometry` key is a geometry zone in the grid's coordinates.
Explicit cells are taken as given, except that reading the budget file refuses a
`(layer, cell)` the model's grid does not have (ledger 343): it matched no row, so a zone
of cells `[20, 99]` on 36-cell layers read cell 20 alone, where MF6 refuses the same
cellid on the observation route.
`to_flopy_obs` writes ONE record per zone cell -- `<name>_c<cell>`, or
`<name>_lay<layer>_c<cell>` for a zone in several layers -- where it wrote one per
locations ROW (every cell of a zone drawn as two rows twice, under one name, which MF6
accepts), and refuses two zones whose labels coincide once written (`"a b"`, `"a_b"`).
`compare()` gives one set of zone columns per zone on both routes -- its `cells` the
`(layer, cell)` pairs -- and a set given no values compares with `time` equal to `per`,
so `residuals()` of a capture template works (it raised `KeyError: 'per'`).

**A target keyed by stress period reads the value at the period's END** -- head,
concentration, reach stage, reach flow and zone flow alike (ledger 336; the user's
rule, "a period's flow should be end of period value"). The record read is the one
that ENDS the period on the clock of the file the value comes from -- the head and
stage files follow `SAVE HEAD`, the budget files `SAVE BUDGET` (ledger 301) -- found
by time, so an ATS period that takes more steps than it declared is read at its real
end. One rule, `myflopy.modflow.mf6.period_end`, serves every family, PEST's
forward run and `model.hds`/`conc`/`temp.get()`. The end is where a record's time
into the period meets the period's length, to 1e-8 of it: FloPy writes that length
at 8 decimals, so an hour held as `1/24` ends MF6's period at `0.04166667`.
**A period whose end the file does not hold is REFUSED**, naming the
period, the steps output control saved and the fix (`SAVE BUDGET LAST` or `ALL`):
under `FIRST`, `STEPS` short of the end or `FREQUENCY` the value does not exist, and
the step before it is a different number. A reach's flow is what it sends downstream
-- to the next reach, a mover, or out of the network (`EXT-OUTFLOW`, so an outlet
reach reads its outflow) -- the quantity MF6's `DOWNSTREAM-FLOW` observation reports;
a zone's flow is its cells ADDED within that one step. Both flows used to add every
step a period saved as well: 3.0x at `nstp=3` under `SAVE BUDGET ALL`.

A lake, reach or zone target reads one of two routes: MF6's own observation, when
EVERY target of the set is DECLARED in the owning package's OBS, else the output file.
**Declared means by identity** (ledger 344): a target's label -- its name as
`to_flopy_obs` writes it -- is declared in the LAK/SFR/DRN package's own OBS as the
family's quantity (`STAGE`, `DOWNSTREAM-FLOW`, `DRN`) and for the target's own lake,
reach or cell (a number from one, a boundname resolved through the package's
`packagedata`, or the cellid). Only the files that OBS names are read, each record
from its own column -- a CSV, or MF6's binary form for a block declared `BINARY` (ledger
345). It used to take the first `*.csv` anywhere in the workspace, in
sorted order, whose header held every NAME: a lake target read another lake's stage, a
GWF head (`heads.csv` sorts first) or a LAK `VOLUME`; an SFR stage target the
`DOWNSTREAM-FLOW` of a same-named flow target; and a stale copy -- a backup, or the CSV
a previous run left in a directory `prepare_run(overwrite=True)` reused -- froze every
PEST realization at its values. A label declared as another quantity only, for another
lake, reach or cell, or twice, is REFUSED naming the target and the record; so is a
set declared in PART (it silently moved the whole set to the file route, which reads
other numbers -- a dry lake's stage, an SFR flow's sign, the LAK stage file's record
clock, and until ledger 346 another clock altogether); a set with nothing declared reads the file, as before. A model carrying several
LAK, SFR or DRN packages needs `package=<name>` on the targets (their lake and reach
numbers each start at 0), carried into PEST's forward run -- the packages are found by
TYPE, so one NAMED `lak` does not hide a `lak_2` beside it (ledger 345); with none it is refused
wherever the targets are declared -- and for LAK and SFR always, since their file
belongs to one package -- while an undeclared DRN zone still sums every package's
drains. `to_flopy_obs` refuses a `kind=` other than the family's quantity and two
names MF6 would upper-case to one label (`"Pond A"`, `"pond-a"`). An SFR stage's
no-data value (3e30) reads `NaN` there, as a lake's does. On the file it reads ONE value per lake per saved
record and ONE per reach or zone per stress period, at the period's end (ledgers 331,
336; it had added a lake's copies, one per connected cell, and a period's steps).
**The observation route is one row per stress period too, at its end, keyed `0 ..
nper - 1`** (ledger 346, the user's rule). MF6 writes a continuous observation at EVERY
time step, labelled by simulated time and nothing else, and the route handed every
step on: a target keyed by period matched nothing (`compare()` gave `NaN`, the PEST
build refused it -- and before ledger 340 dropped it, leaving the STARTING model's value
at weight 1), and under ATS the number of rows moved with the parameters (12, 13, 14 at
K = 10, 1, 2) while PEST reads every realization by the BASE run's line positions -- a
K = 2 realization read t = 23 as the end of the last period, in silence. Which row ends
each period is `period_end.period_end_rows`: each period's time is measured from MF6's
OWN previous period end, because MF6 accumulates the 8-decimal lengths FloPy wrote (an
hourly live model's ends run 3.3e-9 further ahead of the lengths the model holds every
period, past the 1e-8 tolerance from the fourth), and each end is checked by the rule
every route reads. A file short of a period's end -- a run that stopped, a truncated
copy -- is REFUSED naming the period. Every route now labels `time` by stress period,
bar the LAK stage FILE, keyed by saved record (below). **A target keyed by MF6's
simulated time is refused there**, in `compare()` and the build alike, naming the period
each key ends (`10 is the end of period 0 in simulated time -- key it 0`) -- when a key
names no period. A time key BELOW `nper` is also a period number and reads that period
in silence: `1` on daily periods, and on twelve 2-day periods the old keys 2, 4, 6 (the
ends of periods 0, 1, 2) read periods 2, 4, 6 (ledger 347). A set written for the old
spelling is refused only when a key reaches past the last period -- re-key it. PEST templates built on this route at more than one step per period expect a
line per step and must be rebuilt. The
observation route's SFR flow also keeps MF6's sign (`DOWNSTREAM-FLOW` is negative)
where the file route reports the magnitude. On the LAK stage-file route the record
index IS the period only when output control saves one step per period, the step that
ends it. Saving several (`SAVE HEAD ALL` at `nstp=3`), a lake target keyed by period
silently reads a step inside an EARLIER period, while the SFR and head targets beside
it read period end -- key it by saved-record index there, until
`build_lak_stage_result_table` carries the stage file's own `kstpkper` (ledger 301).

**A target's `time` is a zero-based stress period, read by NUMBER, one rule for
`compare()` and the PEST build** (`_shared._numeric_periods`, ledger 340): `1.0` is
period 1, and a key with a fraction (`1.5`, which was truncated to 1) is refused. Two
values for one target in one period are refused in both, naming them -- the build kept
the last. A blank value is unmeasured; a value that is not a finite number (`"12,5"`,
`inf`) is refused when the targets are built, naming the target and its column.

**What PEST is handed is exactly the targets** (ledger 340). `cal.build()` matches each
target to the simulated row whose index EQUALS its period and to its column by name,
and takes the observation pyEMU's instruction file reads there -- never a name rebuilt
from the key's spelling, which missed `per:1` for a key of `1.0`. A target nothing
holds is refused -- a head or concentration target whose location lands on no active
cell included (ledger 341: `compare()` only warns about one, since one target file may
serve several grids). Every simulated row no target measures stays at weight 0 with the base
run's value (a forecast set's other periods too), so nothing unmeasured is
history-matched. A set given no values registers: observed, it only captures; forecast,
it forecasts every row. The head and concentration tables are written by ONE writer for
the base run and every forward run, and two sets under one prefix are refused. A PEST
template built before this bakes the old forward run and must be rebuilt.

Targets: `mf.HeadTargets`, `mf.SfrStageTargets`, `mf.SfrFlowTargets`,
`mf.LakeStageTargets`, `mf.DrnFlowTargets`. Calibration: `model.pest(name,
start_datetime=)` → `PestProject` → `cal.parameterize(...)` / `cal.observe(...)` /
`cal.forecast(...)` / `cal.build(...)` / `cal.run_ies(...)`; discovery via
`model.pest_runs` or `run.pest_runs`. **`model.pest(...)` is the front door** — it
fills `start_datetime` from TDIS and puts the run where `pest_runs` will find it.
`mf.PestProject(...)` can be constructed directly, but that is the advanced path.

**Reviewing a run** (`IesResults`, from `cal.run_ies(...)` or
`model.pest_runs[i].review()`) — `IesResults` is a flat class, not the
`noun.verb` grammar, so its figures are `plot_*` and its tables are bare nouns:

| data | figure | what it answers |
| --- | --- | --- |
| `phi`, `phi_contributions()` | `plot_phi()`, `plot_phi_distribution()`, `plot_phi_contributions()` | did misfit drop, and which group owns it |
| `conflict()` | `plot_conflict()`, `plot_prior_vs_obs()`, `plot_vs_obs()` | is the prior wide enough; does the posterior bracket the data |
| `parameters_at_bounds()` | `plot_parameters_at_bounds()` | is the prior too tight |
| `sensitivity(forecast=…)` | `plot_sensitivity(forecast=…)` | what the data informed, and what drives a forecast — ensemble-based, **not** CSS |
| `forecasts()`, `forecast(name)` | `forecast(name).plot()` | the payoff: posterior prediction spread |
| `field(target)` | `plot_field(target, stat=…)` | property patterns — `mean`/`std`/`base`/`change`/`reduction` |
| — | `plot_field_mosaic(target, stat=…)` | prior vs posterior on one shared scale (`mean`/`std` only) |
| `obs_residuals(prefix=…)` | `plot_obs_residuals(prefix=…)` | *where* the model is biased (heads/conc as points, DRN zones as cells; `prefix=` picks one family so mixed units do not share a color scale) |
| `capture_fields`, `observation_sets` | `report(path)` | what the run recorded; the headline plots bundled to HTML |

Every figure takes `backend="matplotlib"` except `plot_field_mosaic`, which composes
Plotly subplots and says so. Spatial methods need the model grid, which `run_ies()` and
`review()` attach automatically — including from `run.pest_runs`, which resolves the
model on first `review()` so that listing runs stays cheap. (A run holding several
models cannot pick one; pass `review(model=...)` there.) `report(path)` covers phi, obs,
bounds, forecasts and `plot_field` mean/change — not the reduction, mosaic or residual
maps, which are called directly.

### Engine layer (`mf.__engine__`)

The builders/factories under the facade — reach for them only with a specific reason:
`RCHBuilder`, `EVTBuilder`, `SFRBuilder`, `LAKBuilder`, `UZFBuilder`, `MVRBuilder`,
and the `<pkg>_spec` factories (`chd_spec`, `rch_spec`, `evt_spec`, …). E.g.
`mf.rch(context=...)` is `RCHBuilder(...).build()`, and `mf.evt(context=...)` is
`EVTBuilder(...).build()` (both subclasses of the shared `_ArealBuilder`).

---

### UZF without a grid — `mf.uzf` with no `context=`

```python
mf.uzf(nper=12, vks=0.25, thtr=0.08, thts=0.34, thti=0.17, finf=3.0e-5, pet=1.0e-4)
```

Omit `context=` and the package resolves at BUILD time from the model's DOMAIN, so one
declaration builds on any mesh and `with_grid` reads it `rebuilt` (ledger 304). Measured:
75 packagedata rows at ncpl 25 and 243 at ncpl 81, each byte-identical to a fresh eager
declaration on that mesh.

It needs no grid at all — the cell set, the `iuzno` numbering and `ivertcon` come from the
idomain — so it also accepts a bare `ModelContext(domain=...)` or a model with no context,
where the mesh arrives through a hand-written discretization. And unlike every other lazy
form it names no file, so there is nothing to bundle or fingerprint; a recipe of scalars
also means a model carrying a lazy UZF **saves**, where the eager form (records of cellid
tuples) persists nothing.

Refused at the declaration, each still available with `context=`:

| refused lazily | why |
| --- | --- |
| `cells=[(0, 1), …]` | you cannot pin cell numbers the caller typed |
| a per-cell `vks`/`finf`/… | its length is the cell count of the mesh it was written for |
| `mover=True` | a UZF mover endpoint is a hand-typed `iuzno`, an enumeration over ACTIVE CELLS — in range on every mesh and a different column on each |

`mf.sfr` and `mf.lak` still require `context=`.

### Moving a calibrated array to a new mesh — `with_grid(remap=True)`

```python
moved = sim.with_grid(finer_grid, remap=True)          # rule="geometric", uncovered="nearest"
moved = sim.with_grid(finer_grid, remap=True, rule="harmonic", uncovered="raise")
```

Carries per-cell arrays across by area overlay, so `move_report()` reads **`remapped`**
instead of `unverifiable` and the build accepts them. Serves `npf` `k`/`k22`/`k33` and
`ic.strt`/`sto.ss`/`sto.sy` (ledger 300).

| argument | meaning |
| --- | --- |
| `remap` | opt-in. The default does not move — a remap is the inferior route wherever a source exists |
| `rule` | the CONDUCTIVITY rule: `"geometric"` (default), `"mean"`, `"harmonic"`, `"majority"` |
| `uncovered` | `"nearest"` (default, fills and warns in detail) or `"raise"` |

**Use it for an array with no recipe** — a PEST-calibrated K field, an array read out of an
imported model. Where the package has a source, re-declare it instead: that measured 7.1%
worst-case boundary-flux error against this engine's 12.5% at best.

**`rule=` exists because K has no safe default.** The correct average depends on which way
the water flows, which the K field does not know: measured in MF6, the arithmetic mean gave
264× the right flux with the contrast across the flow and was exact with it along.
`geometric` is never best and never catastrophic, and the engine warns naming the cells
where arithmetic and harmonic diverge by more than 2×.

Integer flags (`icelltype`, `iconvert`, `irch`, `ievt`) and everything extensive are
refused by name. Two meshes in different coordinate systems are refused before anything is
carried; a unit rescale under the same declared CRS is the known blind spot.

### A per-cell array from a polygon map — `mf.Polygons`

| call | fills |
| --- | --- |
| `mf.npf(k=mf.Polygons(p, column="k"))` | horizontal K, one value per mapped unit |
| `mf.npf(k22=…, k33=…)` | the other two conductivities |
| `mf.ic(strt=mf.Polygons(p, column="head"))` | starting heads |
| `mf.sto(ss=…, sy=…)` | storage |

```python
mf.npf(k=mf.Polygons("geology.gpkg", column="k", layer="units", default=1.0))
```

A serializable RECIPE, resolved at BUILD against whatever grid the model supplies —
so one declaration builds on two meshes, `move_report` reads it `rebuilt`, and
`Project.save()` keeps it. An `ndarray` in the same slot can do none of those
(ledger 299).

| argument | meaning |
| --- | --- |
| `column` | attribute holding the value. A NULL is REFUSED naming the feature |
| `rule` | `"majority"` (default), `"mean"`, `"centroid"` — see below |
| `default` | value for a cell no polygon reaches. Without it, REFUSED |
| `layer` | layer/table name inside a multi-layer file |
| `layer_field` | column naming the model layer. Omitted, the map is areal and applies to every layer |

**The rule is named, not inherited, because the three disagree.** `majority` gives a
cell whichever polygon covers most of it and so writes only values the map contains;
`mean` is the area-weighted mean and blends across a contact — on a 64-cell grid with
one contact it wrote 12 distinct conductivities where the map holds 2, and its
area-integral ran 0–60% from `majority`'s across 19 contact positions. `centroid`
matches `ZoneSpec.from_polygons`; it differs from `majority` on at most one cell in 64.

`icelltype` and `iconvert` are **refused** — integer flags, and every rule returns a
float. So is a `Surface` or `LayerStack` in any array slot: those describe an elevation
field and used to be accepted, then die inside FloPy naming neither.

For a boundary condition's conductance or elevation this is the wrong tool — the seven
`mf.<pkg>.gpkg` forms already do it. `GeoPackageSource.k_array` is the EAGER sibling:
still supported, but it needs the grid at declaration and its result cannot be saved.

---

### Recharge and ET, array form — `mf.rch.array` / `mf.evt.array`

| call | writes |
|---|---|
| `mf.rch(...)` / `.flopy(...)` / `.gpkg(...)` | list form, one record per cell |
| `mf.rch.array(recharge=…, context=…)` | `READASARRAYS`, one array per period |
| `mf.evt.array(rate=…, depth=…, surface=…, context=…)` | the same for ET |

Array when the boundary covers the model; list when it is sparse or needs
boundnames. `irch`/`ievt` default to `"top_active"`, derived from `context.domain`
— without one MODFLOW 6 applies the flux to layer 1 and silently skips any column
whose layer 1 is inactive. "Top active" is the highest cell with `idomain > 0`, so
a vertical pass-through cell (`-1`) is skipped; it read `!= 0` and lost 42% of a
pass-through stack's recharge until ledger 319. Every "top active" in myflopy is
`grid.geometry.top_active_layer`. The index is zero-based, and a typed one is written
as given in either mode (an `ievt` beside an omitted ET `surface` was dropped until
ledger 324). Segmented ET has no
array form and `nseg` raises. The build checks the index against the DISV it is
built beside: a per-column index naming a REMOVED cell (`idomain <= 0`) in a column
that has an active one is refused (ledger 321 -- most often a `context=` from before
an excavation, or, for a package declared without one, the MODEL's context `domain=`),
and `irch=None` / one layer everywhere, the documented opt-out, is WARNED with the rate
it loses (ledger 322) -- exact for recharge, "up to" at the maximum rate for ET.

### Horizontal flow barriers — `mf.hfb`

A barrier sits on the **face between two cells**, so MODFLOW 6 addresses it as a
cell *pair*. FloPy takes any pair at all and never checks it is a real connection;
the run aborts instead. These resolve geometry and validate first.

| call | from |
|---|---|
| `mf.hfb(pairs=…, hydchr=…, layers=…, context=…)` | explicit cell pairs |
| `mf.hfb.line(trace, context=…, hydchr=…)` | a fault trace → the connections it cuts |
| `mf.hfb.gpkg(path, context=…, hydchr="hydchr")` | a GeoPackage line layer |
| `mf.hfb.enclose(polygon, context=…, hydchr=…)` | a watertight cutoff wall around a region |
| `mf.hfb.flopy(stress_period_data=…)` | the FloPy-native form, unvalidated |

`hydchr` is a hydraulic **characteristic**, 1/T — barrier K divided by barrier
thickness, not K. `0` is impermeable; a negative value is a conductance multiplier.

A trace cuts the **connections** whose centre-to-centre segment crosses it -- not the
faces it passes through, which run between cells ALONG the fault and block nothing
across it (68-100% of the flow still crossed until ledger 216). A cell centre lying
on the trace belongs to its left, looking along the trace as drawn. `.enclose` seals
a REGION: cells inside by centroid, and it refuses a cut that leaks.

Grid verbs underneath: `vor.barrier_faces(geom)`, `vor.enclosed_faces(polygon)`,
`vor.shared_face(a, b)`.

### Importing an existing MODFLOW-USG model

| call | returns | notes |
|---|---|---|
| `mf.read_usg(nam, gsf=, crs=, nper=, require_layered=, read_boundaries=)` | `UsgModel` | A USG `DISU` has no coordinates, so the `.gsf` is required for anything but inspection. |
| `usg.report()` | `str` | What converted, what was approximated, what was omitted — with counts. |
| `usg.validate()` | `list[Finding]` | What MODFLOW 6 will object to, before running. `severity` is `reject` (will not start) or `allow` (starts, and the boundary behaves differently) — which one depends on the layer's `icelltype` (ledger 291). |
| `usg.to_mf6(name, crs=, newton=, start_date_time=, complexity=, include=, fix_for_mf6=, local_origin=)` | `SimulationSpec` | DISV. Keep `local_origin=True`. `fix_for_mf6=True` lifts a GHB head / DRN invert / RIV **bed** back inside its cell, a river as a pair so its depth survives. |
| `usg.surfaces()` / `usg.boundary_frame(ftype)` / `usg.cln_polygons()` | `Surface`s / GeoDataFrame / GeoDataFrame | Grid-independent — these survive a change of grid, node numbers do not. |
| `usg.export_gis(dir, crs=, stream_lines=, clip_to=, resolution=, bed_thickness=, overwrite=)` | `ExportManifest` | Writes the whole model as GIS files that outlive **any** grid. See below. |

`UsgModel` also exposes the model reshaped to `(nlay, ncpl)`: `top`, `botm`, `idomain`,
`strt`, `k`, `k33`, `ss`, `sy`, `icelltype`, `thickness`, `uppermost_active`.

#### `export_gis` — the form that survives re-gridding

`to_mf6()` gives a simulation on the grid it was converted on. `export_gis()` gives the
model's *content*, so you can rebuild it on whatever mesh you end up with — which is
usually why a USG model is imported at all.

| written | what it is | feeds |
|---|---|---|
| `boundaries.gpkg:<pkg>_segments` | each list BC as **one line per source cell** | `mf.drn/ghb/chd/riv` |
| `boundaries.gpkg:source_cells`, `:domain` | the old mesh and its outline, so values stay auditable | — |
| `streams.gpkg:streams` | one `LineString Z` per CLN stream, Z = bed elevation | `mf.sfr` |
| `streams.gpkg:stream_nodes` | a point per old cell, with `station` for interpolation | `mf.sfr` |
| `lakes.gpkg:lakes` / `:lake_nodes` | a polygon per CLN waterbody, and its per-node bed | `mf.lak` |
| `arrays/lakebed_bottom.tif` | per-cell lakebed — pass this as `lake_bottom=`, not the column | `mf.lak` |
| `lake_forcing.csv` | the CLN P−ET wells as a rate over each lake | `mf.lak(rainfall=, evaporation=)` |
| `arrays/*.tif` | every per-cell array, **one file per distinct array** | `mf.Raster` |
| `cell_values.gpkg` | the same arrays at the cell centres — exact, unlike the rasters | `Surface.from_points` |
| `periods.csv` | timing, and which array each package uses in each period | — |
| `manifest.json`, `README.md` | every column, its units, and every assumption made | — |

**A model that fails to converge mid-run? `mf.tdis(ats=True)`.** Adaptive time
stepping lets MODFLOW 6 choose its own step within a period -- and, crucially,
**retry a failed step at a smaller one instead of ending the run**
(`ats_dtfailadj`, default 5.0; `0` restores the stop-on-failure behaviour). Periods
are selected zero-based: `ats=True`, `ats=[7, 9]`, or `ats={7: {"dtmin": 1e-4}}`.

**`mf.Spread` and `mf.CellSurfaceOffset` both SERIALIZE** (ledger 274) — a package
holding one saves as reviewable JSON under a `$spread` / `$cell_surface_offset` token
instead of falling through to the pickle fallback, and a hand-edited payload is refused
where you can act on it. They are **refused on the ARRAY forms** (`mf.rch.array`,
`mf.evt.array`): both describe a per-FEATURE resolution and the array form has no
features, so neither could ever act — use the `.gpkg` list form, where each feature
resolves its own value.

**Building from a line BC? Wrap the conductance in `mf.Spread`.** A feature's value is
otherwise written to every cell it covers, which is right for an elevation and
multiplies a conductance by the cell count (measured: 2.02x, on an unchanged grid).
`conductance=mf.Spread("conductance")` gives each cell its share — by length for a
line, area for a polygon, and the whole value for a point, which has nothing to divide.
A polygon covers a cell by AREA: one drawn along cell faces reaches only the cells
inside it, not the neighbours it shares an edge with (ledger 210). A line lying ON a
face is split evenly between its two cells, and a point on a face or vertex goes whole
to the lowest-numbered active cell it touches (ledger 212). A line PIECE that is noise
goes the way an areal one does (ledger 237): noise is under the 0.01 tolerance, in the
grid's length unit, AND under 1% of the cell's own mean width — so a sub-unit lab-scale
grid keeps its coverage, and a feature shorter than the tolerance in TOTAL is what was
drawn and keeps every piece of itself. That is what stops a line ending ON a face and
overshooting it by 1e-6 ft — what a CRS round trip does — from writing a second record
for the cell beyond, at the same elevation and conductance; `mf.sfr` has dropped such a
reach since ledger 214. The same rule serves the cell selector, so `map(select=<line>)`,
a region and a line head target move with it, not only a `.gpkg` source.

Runnable walkthrough: **`examples/mf6/notebooks/usg_export_to_new_model.ipynb`** —
each written file, and the call that consumes it. It builds no model; it lays out the
workflow, and runs end to end against a demo grid so the shapes are visible first.

Three things worth knowing, each measured rather than assumed:

- **Column names are the builders' parameter names.** `streams` carries `rwid`/`rhk`/
  `rgrd`/`rbth`/`man`; `lakes` carries `strt`/`lake_bottom`/`bedleak`. So
  `mf.sfr(streams=..., width="rwid", streambed_k="rhk", ...)` works with no renaming. The one
exception is `lake_bottom`: pass `arrays/lakebed_bottom.tif`, since the per-lake column is a
single number for a bed that varies and `mf.lak` then rejects the cells it misses.
- **Every conductance has a `*_per_ft` twin.** A list BC here is a *line* (measured: mean
  in-set neighbour counts of 1.87–2.52, against 4–6 for anything areal), and splitting a
  cell in two leaves the elevation alone while doubling the conductance. Rebuild it as
  `cond_per_ft × length in the new cell` and the feature total is preserved. Do **not**
  treat it as a leakance — `corr(cond, cell area) ≈ 0` over 205 DRN records. The chain
  those lengths come from is joined by a shared **face** — `grid.shared_face`, the rule the
  HFB views resolve a barrier's own face against, and pair-for-pair the older
  `adjacent_cells_idx` on every mesh measured (ledger 244). A shared *vertex* is
  not a connection, so on a structured or quadtree mesh a diagonal staircase of six 25 ft
  cells measures the 125.0 ft it occupies rather than 145.71, and a 2x2 block is a ring
  (75.0 ft, `conductance_per_ft` 4.00) rather than a clique whose walk crossed the corner
  (85.36 ft, 3.31); a Voronoi mesh is unaffected, its vertices being 3-valent. A cell left
  alone on its chain — nothing shares a face with it, or it is a branch arm the open walk
  could not reach — has `length_ft` 0 and a NaN twin, counted in the manifest as
  `cells_alone_on_a_chain` with a note and a logged warning; take those cells from the
  point layer. Reading the layer with `conductance=mf.Spread("conductance")` recovers the
  feature total exactly either way — it is the per-foot hand recipe that goes NaN.
- **A repeated period is written once.** Ten Trails is 72 periods and a strict 12-month
  cycle, so recharge is 11 distinct arrays, not 72; `periods.csv` is the index.

What CLN cannot supply is written as a documented default rather than omitted, so a
missing column never becomes a silent zero: SFR's `rbth` and `man` (CLN routes with a pipe
conductivity, which has no MF6 counterpart), and `bedleak`, which is `FSKIN / bed_thickness`
— CLN has a *skin*, not a bed, so the thickness is a choice. `fskin` is written unmodified
beside it, and `mf.lak` takes one `bed_leakance` per lake while `lake_nodes` keeps the
per-node variation.

## B · Read side — the `noun.verb` grammar

```
model.<pkg>.<inputs|results>.<noun>.<verb>()
```

`model.<pkg>` is myflopy's view of one package; **FloPy's own package object is
`model.fl.<pkg>`** (`model.fl.drn.stress_period_data`; any other package by name --
`model.fl.npf`, `model.fl.disv`). `model.packages` holds what spans packages --
`summary()` and `mosaic()`. Group and diff scopes still read
`group.packages.<pkg>` / `diff.packages.<pkg>` (ledger 339).

Every **noun** is an object; every object answers the same **verbs**. Different plot
*types* are different verbs — `plot()` is not a dispatcher.

### The verbs

| verb | available on | returns |
|---|---|---|
| `get()` | every noun | the DataFrame |
| `summary()` | every noun | a compact digest frame |
| `plot()` | every noun | the noun's **one** canonical chart (`viz.Fig`) † |
| `map()` | **spatial** nouns | choropleth on the grid |
| `section()` | spatial nouns | cross-section |
| `mosaic()` | spatial nouns | small-multiples panel |
| `animate()` | spatial nouns | animation over stress periods |
| `long()` / `wide()` / `stack()` | cell-budget fields | frame reshapes |
| `__call__(per=)` | period-bound nouns | a fresh view rebound to `per` |

† Exception: `lak.results.q.budget.plot()` currently returns a **matplotlib** figure
(the `viz.Fig` conversion is deferred — ledger 52). Derived-table nouns (`profile`,
`budget`) are not spatial, so they have no `map`/`section` (ledger 17); use the field
noun's own `map()` (e.g. `sfr.results.q.map()`).

### Nouns per package

**List BCs** (`chd/ghb/drn/riv/wel/rch/evt`):
- `inputs.<field>` — the record fields above (e.g. `riv.inputs.stage`,
  `evt.inputs.rate`, `rch.inputs.recharge`) — each a spatial field noun.
- `results.q` — the groundwater exchange (spatial). The `.get()` column names its
  reference frame: **`q_gwf`** (aquifer `.cbc` record; negative = the feature gains).

**SFR** (`model.sfr`):
- `results.q` — stream↔aquifer exchange (spatial; `q_gwf`).
- `results.stage` — reach stage (spatial).
- `results.profile` — the merged multi-field reach profile (`get`/`summary`/`plot`).
- `results.q.profile` / `results.stage.profile` — single-field reach profiles
  (`get`/`summary`/`plot`). *The frame is `…profile.get()`, not `…profile()`.*
- `budget.<term>` — budget terms: `flow_ja_face`, `gwf`, `ext_inflow`, `ext_outflow`,
  `rain`, `runoff`, `evaporation`, `storage`, `from_mvr`, `to_mvr`, `mvr`,
  `stream_fluxes`, `auxiliary`, `types`.

**LAK** (`model.lak`):
- `results.q` — lake↔aquifer exchange (spatial; column **`q_lake`**, feature-referenced,
  negative = the lake loses).
- `results.q.budget` — exchange summarized by connection type (`get`/`summary`/`plot`).
- `results.stage`, `results.stage_change`, `results.connections`. A lake with no
  stage -- DRY (MF6 writes -1e30) or INACTIVE (+1e30) -- is `NaN` in `results.stage`,
  a gap in its series, blank on its maps and an undefined `stage_change`, never a
  number 1e30 ft away; `model.outputs.lak.stage.get(masked=False)` is the raw file
  and the one place dry and inactive differ (ledger 329).
- `budget.<term>` — `gwf`, `flow_ja_face`, `ext_inflow`, `ext_outflow`, `rainfall`,
  `runoff`, `evaporation`, `storage`, `from_mvr`, `to_mvr`, `mvr`, `withdrawal`,
  `constant`, `lake_fluxes`, `auxiliary`, `types`.

**UZF** (`model.uzf`): `results.gwrch`, `results.sat`, `results.fields`.

**Static arrays** (`model.npf` / `.ic` / `.sto`): one noun per STORED array --
`npf.k`/`k22`/`k33`, `ic.strt`, `sto.ss`/`sy`. Plus two DERIVED nouns on `npf` (ledger
308): **`npf.transmissivity`** -- K times the SATURATED thickness at a stress period
(`per=`, its last saved time step), by MF6's own ICELLTYPE rule: `min(head, top) -
bottom` in a convertible layer, the full thickness in a confined one, matching MF6's
`DATA-SAT` -- and **`npf.transmissivity_max`**, K times the full thickness, which needs no
run. Both answer `get` (with `k` and `thickness` columns; `per` for the saturated one) /
`summary` / `map` / `mosaic(by="layer")` / `animate()`, and draw THROUGH
`model.plot.map(type="trans" | "trans_max")`, so the noun and the free map are one
picture. They are not in `npf.fields`, which lists what the package file stores.
Transmissivity is WRITTEN as a whole number with thousands separators -- `112,000 ft²/d`,
never `1.12e+05` -- in the hover and on every colorbar, both backends (ledger 313); a
value below 1 keeps its significant digits.

**PRT** (a finished run's `PRTRunResults`, not `model.packages`): `results.pathlines`
(trajectories — `map` draws polylines, not cells), `results.travel_time`,
`results.endpoints`, `results.capture`. These are derived, time-integrated views, so
they reject `per=`; `pathlines` also refuses `section`/`animate` (a trajectory is not a
per-cell field to slice) and `travel_time`/`endpoints`/`capture` refuse `animate`.
The raw MF6 track table is `results.track_records`. Full detail in §A above.

**Model-level reads:** `model.hds` (GWF heads explorer — `get`/`summary`/`array`/
`map`/`section`/`mosaic`/`animate`), and its transport twins **`model.conc`** (GWT
concentration) and **`model.temp`** (GWE temperature) — the same dependent-variable
grammar, reading the `.ucn` binary via the shared `DependentVariableFile` base. The
readers are **kind-gated**: `.hds` on a transport model (or `.conc` on a flow model)
raises a clear error. `model.outputs.<pkg>.bud` (raw budget accessor),
`model.targets.<family>` (`compare`/`stats`/`calibration_plot`), `model.pest_runs`.

**`model.budget.<term>`** — every term in the model's *own* budget file, as a spatial
noun (`get`/`summary`/`plot`/`map`/`section`/`mosaic`/`animate`). Terms are **discovered
from the file**, not declared, so they follow the model's kind and packages:

```python
model.budget.types                     # the MF6 record names actually present
model.budget.source_sink_mix.get()     # GWT/GWE: the SSM term, per cell
model.budget.storage_aqueous.map()     # GWT storage (GWE spells it storage_cellblk)
model.budget.drn.summary()             # GWF: any boundary package
model.budget.sto_ss.get(per=0)         # ...and terms with no package accessor
model.budget["SOURCE-SINK MIX"]        # indexing takes either spelling
```

Available on **every** model kind, not just transport: the plumbing is kind-neutral,
and `STO-SS`/`STO-SY`/`DATA-SPDIS`/`DATA-SAT` have no package accessor at all.
`FLOW-JA-FACE` appears in the namespace but *refuses* `get()` — it is indexed by cell
connection, not by cell (ledger 91/96).

Distinguish it from its three neighbours, which are genuinely different things:

| spelling | what it reads |
|---|---|
| `model.budget.<term>` | the model budget file, per cell, one noun per term |
| `model.bud(pkg)` | the same file, legacy compatibility wrapper, raw frames |
| `model.budget_cumulative` / `_incremental` | the **listing** file: whole-model totals |
| `model.<pkg>.budget.<term>` | the **package-output** file — a different file, feature-first node layout |

The column is MF6's raw `q`, unlike `model.<pkg>.results.q`, which renames it
`q_gwf`/`q_lake` to carry the reference frame. The values are identical; only the name
differs (ledger 97).

**Each of those files also has its own SAVE CLOCK, and `model.kstpkper` is the head
file's** (ledger 301). One MF6 `OC` `saverecord` sets two independent flags — `SAVE HEAD`
drives `idvsave`, `SAVE BUDGET` drives `icbcfl` — and LAK/SFR gate their *stage* output on
the first while gating their *budget* output on the second. So `SAVE HEAD LAST` +
`SAVE BUDGET ALL` at `nstp=5` writes 2 head records against 10 LAK budget records, and
`STEPS 1` vs `STEPS 3` writes two of each naming different time steps. Every table built
here labels its rows from the file it actually read; **never pair records from one file
with `model.kstpkper`**. `modflow/utils/outputs.py`'s `SFRStageOutput.get` is the pattern
to copy. The budget filename comes from the OC record too — it is not always
`<name>.cbc`.

**Grouped reads:** `group.hds` / `group.conc` / `group.temp` are the multi-model twins
of the model-level readers, with the same verbs plus `compare()`:

```python
group.conc.get()                     # aligned concentration, one block per member
group.conc.compare()                 # vs the reference: conc / reference_conc / diff
group.conc.map(); group.conc.plot(); group.conc.section(line=line)
group.diff().conc.map("variant")     # Δconc choropleth -- the ONE public diff route
group.diff().temp.summary()          # max/mean |Δtemp|, RMSE, argmax cell/layer
```

All three are the same class configured by five attributes (`_GroupFieldView`), so a
verb added to one is available to all. `diff` is always **model − reference**. Every
member of a `group.conc` must be a GWT model; the kind-gated reader supplies the error
otherwise (ledger 100).

### Deprecated spellings (D12)

The following resolve through `__getattr__` with a `DeprecationWarning` and stay out
of `dir()`/completion; they return exactly what they used to. Prefer the noun form:

| deprecated | use instead |
|---|---|
| `sfr.results.long_profile()` / `plot_long_profile()` | `sfr.results.profile.get()` / `.plot()` |
| `sfr.results.q.plot_profile()` / `sfr.results.stage.plot_profile()` | `…profile.plot()` |
| `lak.results.q.budget_summary()` / `plot_budget()` | `lak.results.q.budget.get()` / `.plot()` |

---

---

## C · Drawing — the plotting verbs

**Highlighting a subset of cells** — `select=` is an option on `map` at every
scope (free `mf.plot.map`, `model.plot.map`, `vor.plot.map`, and every
`…inputs/results.map` leaf):

```python
model.plot.map(layer=0, select=cells)                    # outline (default)
vor.plot.map(select="all_streams")                       # a registered region
model.plot.map(select="wetland.gpkg")                    # a vector file
vor.plot.map(values, select=cells, select_style="dim")   # fade the rest instead
vor.plot.map(select=cells, select_color="#0072B2")       # avoid a red colorscale
```

**The window a map opens on** — `extent=` is an option on `map` at every scope and on
both backends (ledger 309), and the one way to frame a static map, which cannot be
zoomed. It sets the VIEW only: no cell is cropped and the colour scale still spans
the whole field (pass `zmin`/`zmax` for a close-up). Four numbers are
`(xmin, ymin, xmax, ymax)` in the grid's CRS -- shapely's order, not FloPy's -- and a
box off the grid raises printing both boxes and both CRSs. In a mosaic, one extent
frames every panel and a list frames each on its own (a list is never synced).

```python
model.plot.map(layer=0, extent=(500, 600, 1000, 1100), backend="mpl")   # a static close-up
model.plot.map(layer=0, extent=[(500, 600), (1000, 1100)])              # points -> their box
model.plot.map(layer=0, extent=transect_line)                           # any geometry's bounds
model.npf.transmissivity.mosaic(per=3, extent=site_gdf)        # every panel
model.hds.mosaic(by="layer", extent=[site, None, north, None])          # one per panel
model.add_aoi("site", (500, 600, 1000, 1100))                           # a named area...
model.plot.map(layer=0, extent="site", backend="mpl")                   # ...used by name
```

Named areas live in `model.aoi` / `project.aoi` (ledger 310): defined on a model,
kept by its project, and usable as `extent=`, `select=` and `features=`. See
`docs/project_workflow.md`, "Areas of interest".

**Custom geometry over a map** — `features=` is an option on `map` at the same
scopes (ledger 254): a vector file path, a `GeoDataFrame`, a shapely geometry,
a `mf.plot.Feature(...)` carrying its styling, or a list/dict of those. One
legend entry and an attribute hover per layer, drawn by both backends and kept
by `mosaic`/`animate`. `locs=` is deprecated in its favour.

```python
model.plot.map(layer=0, features="wells.gpkg")                       # first layer of the file
model.plot.map(layer=0, features=mf.plot.Feature("site.gpkg", layer="ponds", fill=True))
model.npf.k.map(logscale=True, features={"Wells": wells_gdf, "Boundary": boundary})
```

**Log colour scales** — `logscale=True` on any map colours `log10(value)` but
labels the colorbar in real units at whole decades on both backends; the hover
keeps the real value and `zmin`/`zmax`/`zmid` are given in DATA units
(`model.npf.k.map(logscale=True, zmin=1e-4, zmax=10)`). `colorscale=`
takes any Plotly name (`_r` reverses it) or explicit stops (ledger 252/253).


Six verbs cover every picture. They exist as free functions and bound to the
objects, and the bound form calls the free one, so the two cannot diverge.

```python
from myflopy import plot

plot.map(model, layer=0)      # or  model.plot.map(layer=0)
plot.grid(vor)                # or  vor.plot.grid()   -- or just vor.plot()
```

### The verbs

| verb | draws | notes |
|---|---|---|
| `map(source, values=, backend=)` | plan view | contours, locations, hillshade, pathlines are **options**, not verbs; `"plotly"` or `"mpl"` |
| `section(source, backend=, fill=)` | vertical slice | line profile by default; `fill="layer"`/`"results"`/array draws the CELLS (grid + layers + a field), `backend="mpl"` only |
| `surface(source, ...)` | 3-D height field `z(x, y)` | Plotly only — there is no matplotlib rendering of a height field here |
| `grid(source, backend=)` | the mesh itself | `"plotly"` flat 2-D (no CRS needed), `"mpl"` FloPy's patch renderer, `"vtk"` the 3-D volume |
| `mosaic(panels, ncols=)` | many pictures, one figure | **combinator** — takes pictures, not a subject |
| `animate(frames, backend=)` | frames in sequence | **combinator**; `"plotly"` live, `"png"` rasterized |

### Which scope answers which

A scope answers the verbs it can *know*: a bare grid has no results, a layer
stack has no time.

| scope | verbs |
|---|---|
| `myflopy.plot` | all six |
| `model.plot` | all six |
| `vor.plot` | `map` `section` `grid` |
| `stack.plot` | `map` `section` `surface` `grid` |

Pinned by `tests/test_plot_vocabulary.py`, which fails naming both the scope and
the stray verb if one drifts.

### The parameters are in the signature, not just the docstring

Every verb spells its parameters out in the `def`. `plot.map` names 35 of them;
`model.plot.map` names the same 35. That duplication is deliberate and is the
only thing that works: **PyCharm and Pylance are static analyzers** — they read
the `def` line and never execute the module, so `__signature__`, `functools.wraps`
and a beautifully written docstring reach `help()` and reach no editor at all.
A `**kwargs` forwarder shows the caller nothing.

Each verb still keeps a `**kwargs` tail, but only where the tail is *genuinely*
open — `map`'s rides through to the `go.Choroplethmap` trace, whose parameter
names belong to Plotly and are validated at render time.

Three tests keep the copies honest, and all three fail by name:

| test | catches |
|---|---|
| `test_a_named_default_matches_the_link_that_owns_it` | a restated default that drifted from the function that owns it |
| `test_the_bound_model_verb_mirrors_the_free_one` | a parameter added to `plot.map` and forgotten on `model.plot.map` |
| `test_every_signature_parameter_is_documented` | accepted-but-undocumented **and** documented-but-not-accepted |

Narrower scopes may offer *fewer* parameters, never other ones: `vor.plot.map`
drops `per`/`layer`/`type` because a bare grid has no results to read. Anything
dropped still rides the `**kwargs` tail, so narrowing a signature never narrowed
what already worked.

**One deliberate exception.** The namespace-level `field=` sugar —
`packages.lak.results.map(field="stage")` — stays `(*args, field=None, **kwargs)`.
It dispatches to accessors with genuinely incompatible signatures (a lake's
`connections` map has no `per=`; `DrnInput.map` takes `per` positionally), so a
single merged signature could only be produced by lying about one of them. Call
the leaf — `packages.lak.results.stage.map(` — when you want completion; every
leaf verb is explicit.

### Everything returned is a Picture

```python
picture = model.plot.map(layer=0)
picture                       # renders itself in Jupyter
picture.fig                   # the figure, to adjust before display
picture.show()
picture.save("heads.png")     # suffix picks the format
picture.html("heads.html")    # standalone page
```

There is never a trailing `.plot()`. Three renderers back this: most pictures are
Plotly; `MplPicture` (Matplotlib) offers `.axes` and `VtkScene` (PyVista) offers
`.scene` instead of `.fig`, and their `.fig` raises saying so. A `backend="mpl"`
picture is the exception -- a bare Matplotlib `Figure` -- and that is
`model.plot.section()`'s default (the cross-section by layer; `fill=None` is the
head line profile, `backend="plotly"` the interactive picture).

### Standalone HTML

Two artifacts for two jobs. `plot.animate(frames).html(path)` writes an
interactive page — for a choropleth it ships the cell geometry once and restyles
per frame rather than re-embedding the grid. The `mf.export_*_slider_html`
functions instead pre-render frames to PNG and page through them: their size is
independent of cell count, and they carry the machinery a long export needs
(`resume=`, `progress=`, external frame directories, a concurrent-writer lock).
They render through FloPy's `PlotMapView`, so they are a *different picture*,
not a second spelling.

> Runnable tour: `examples/mf6/notebooks/plotting_vocabulary_tour.ipynb`.
> Normative rules: `docs/view_layer_conventions.md`.

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
model.sfr.results.profile.plot(per=5)      # merged reach profile
model.sfr.results.q.profile.get()          # q_gwf along the stream
model.lak.results.q.budget.get(per=5)      # exchange by connection type
model.riv.results.q.map(per=5)             # exchange choropleth
model.sfr.budget.flow_ja_face.get()        # in-channel routing flow
```
