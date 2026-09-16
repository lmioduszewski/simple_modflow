# myflopy / simple_modflow — Project Context for Claude

## Repository
- GitHub: `lmioduszewski/simple_modflow`
- Active development branch: `myflopy`
- Package lives at `src/myflopy/` on that branch

## Adding a package? Read plan §4.7 FIRST (2026-07-18)
A 6-angle agent sweep measured **92 sites in `src/` that encode per-package knowledge**,
carrying **33 unintended gaps** today. `mf.riv`/`mf.evt` followed the documented
four-piece checklist, under two adversarial reviews, and still left gaps; `wel` is
missing from the artifact subsystem, and `model.chd/drn/ghb/riv/wel/evt` accessors
didn't exist while `model.rch/uzf/sfr/lak` did — rot that predates riv/evt.
(Both are now FIXED: artifacts in 4.7.3, the six accessors in 4.7.7.) Plan §4.7 consolidates
these onto the existing `package_registry.py` descriptor and lists the sites. Note 6 of
the gaps are INTENTIONAL (deprecated/frozen legacy tiers) and must not be "fixed".

**What adding a package actually costs (measured 2026-07-18, plan 4.7.6):** one
descriptor entry reaches **9 surfaces automatically** (diff tiers, mover list,
budget basing, suffix map, all three artifact sets, budget term). **7 still need
hand-writing**: the `<pkg>_spec` factory, the `package_api` helper class and the
`GeoPackageSource` resolver — all deliberate, they carry prose and explicit
signatures — plus four namespace properties (`ModelPackages`, `SimulationBase`,
`GroupPackages`, `_PackageDiffNamespace`). Those four stay hand-written **on
purpose** (4.7.7): a `.pyi` stub replaces a whole module for type checkers, so
generating them would mean hand-maintaining stubs for every other public name in
four large modules. Instead `test_package_descriptor.py` asserts all four
namespaces cover the whole registry, so a half-wired package fails immediately. `tests/test_package_descriptor_payoff.py` is the
machine-checked version of this list and fails in BOTH directions: adding a
hand-written site, or closing one without recording it.

**As of 4.7.2 (2026-07-18) `package_registry.py` is the single source of per-package
truth** — FloPy class, record fields, `.gpkg` defaults, capabilities, tier flags,
file suffix, budget `node2` basing, prose blurb, for all 10 packages. Adding a package
means adding ONE descriptor entry there. The ~92 hardcoded sites still exist and are
still authoritative until 4.7.3 deletes them one at a time; until then
`tests/test_package_descriptor.py` proves the descriptor equals them, asserting
against the live FloPy signature/dfn and the real lists rather than repeating
literals. Its `RETIRE WITH 4.7.3` tests must be deleted as each list goes, or they
become circular and stop testing anything.

**Changing the package API?** `tests/api_snapshot.json` pins every helper signature,
`*_spec` output, and per-surface package list. If a change is intentional, rerun
`python scripts/derive_api_snapshot.py` and **review the diff** — it is the public-API
change log. An un-regenerated snapshot fails CI.

## Adding a table, plot, or map to a package? Read `docs/view_layer_conventions.md`
The rule is `model.packages.<pkg>.<inputs|results>.<noun>.<verb>()` — **every noun
is an object, every object answers the same verbs** (`get`/`summary`/`plot`, plus
`map`/`section`/`mosaic`/`animate` for spatial nouns). Never add a `foo()` + `plot_foo()`
method pair; add a view class exposed as a noun property.
Figures are **always `viz.Fig`**, never a bare `go.Figure` — a bare figure silently
drops `scrollZoom`, `dragmode='pan'`, and the house template. Colors come from the
policy helper, not hex literals at the call site.
**Drawing anything at all? Six verbs, and only six** (Phase 8, 2026-08-20):
`map` `section` `surface` `grid` (pictures) + `mosaic` `animate` (combinators over
pictures). They are free functions in `myflopy.plot` AND bound as `model.plot` /
`vor.plot` / `stack.plot`; the bound form calls the free one. **Geometry chooses
the verb** — contours, locations, hillshade and pathlines are OPTIONS on `map`,
never verbs; **`backend=` switches the renderer, never the subject** (which is why
the 3-D volume is `grid(backend="vtk")`, not `surface(...)` — `surface` means a
height field; `surface(backend="vtk")` exists and draws those height fields as
separate PyVista sheets, which is the rule holding rather than bending). Everything returned is a `viz.Picture`: `.fig` / `.show()` /
`.save(path)` / `.html(path)`, and never a trailing `.plot()`.
**A cross-section of the MODEL (grid + layers + results) is
`section(fill=..., backend="mpl")`** (ledger 177, 2026-09-02) -- `fill="layer"`
colours cells by layer and draws the head as a water surface; `fill="results"`
colours them by head/conc/temp with a colorbar; a `(nlay, ncpl)` array paints
anything else. Plain `section()` is a LINE PROFILE (head vs distance) and draws
no cells -- that difference is what sent this looking twice. `fill=` is mpl-only
and raises rather than degrading: on the plotly backend, on any profile-only
argument, on `"results"` for a bare grid, and on a flat `(ncpl,)` array.
Pick what it contains with `layers=` (which layers' cells; masks and CROPS the
view), `head_layers=` (whose water levels ride on top -- default `0`, a list
draws one labelled surface each, `None` none) and `layer_labels=` (unit names for
the legend, defaulting to `model.myflopy_context.surfaces.names` when the model
was declared with a `LayerBuildResult` as its `surfaces`). **`layer=` RAISES
here** -- it means "overlay these layers' head profiles" on the profile branch.
RETIRED as the spelling to teach: `mf.plot_model_cross_section` (still there,
still the engine, but on no scope and not exported top-level).
**Want it static? `backend="mpl"`, on `map` / `section` / `grid`, at EVERY scope**
(ledger 176, 2026-09-02) — free verb, bound namespace, and noun. The NOUN-tier
`mosaic`/`animate` take it too, which is where it earns its keep because there is
no intermediate picture to call `.plot_mpl()` on; the free `mosaic` has no
`backend` at all and the free `animate`'s is `{"plotly", "png"}`, a rasteriser
over frames rather than this switch. It returns a bare `matplotlib.figure.Figure`
everywhere (`normalize_backend` resolves the aliases, `as_mpl_figure` the three
renderers that each returned something different). **`.plot_mpl()` still works
and is NOT the spelling to teach** — it is what `backend=` calls underneath.
`surface` deliberately has none: no matplotlib rendering of a 3-D height field
exists here, and a parameter that cannot act is the same defect as a hidden one.
Two traps: `Choro.plot_mpl` **ignores overlays**, so a picture whose subject IS
an overlay (both HFB verbs) must draw it onto the axes itself or return a valid
figure with its subject missing; and `PALETTE.highlight` is Plotly's
`rgb(...)`, which matplotlib rejects — use `PALETTE.mpl_highlight`.
**3-D scenes: build-time options on the verb, post-render handles on the Picture**
(ledger 157, 2026-08-28). Both layer scenes draw ONE NAMED ACTOR PER LAYER —
`scene.scene.actors["clay"].visibility = False` — and every cell carries
`layer`/`thickness`/`top`/`botm`/`cellid` whichever `color_by` draws, so
`scene.meshes[0].save("stack.vtu")` hands ParaView the lot. Two traps worth
carrying: flopy's `Vtk.add_array` **NaN-masks FLOAT arrays at `idomain == 0`**
(so attach floats after `to_pyvista()`, never through it), and its own `top`
array is NaN for every layer below 0. What a WRITTEN page cannot do, measured:
no widgets, no picker, no keyboard, and only the ACTIVE scalar is serialized.
`tests/test_plot_vocabulary.py` pins the verb set at every scope and fails naming
the stray verb, so adding one in one place and forgetting another is caught.
**Highlighting cells is `map(select=...)`, on EVERY scope** (ledger 160,
2026-08-29) -- cell indices, a boolean mask, a registered region name, a vector
file, or a geometry. It draws a **dissolved-boundary outline** by default, NOT
Plotly's `selectedpoints`: selection styling on a choropleth exposes only
`opacity`, so it can only dim the rest, which measured a **6.9x loss of readable
contrast** across the unselected cells on a real head field. `select_style="dim"`
restores the old picture, `"both"` draws each. Applied on `Choro` (not in the
verb) so `mosaic`/`animate` keep it, and `mode="lines"` so a user's lasso cannot
clobber it.
**Log colour scales and custom geometry are OPTIONS on `map`, at EVERY scope**
(ledger 252–254, 2026-09-13). `logscale=True` colours log10 but LABELS the
colorbar in real units at whole decades on both backends (`viz.log_decade_colorbar`;
a `LogNorm` on mpl); `zmin`/`zmax`/`zmid` are then DATA units (`zmin=1e-4,
zmax=10`), a non-positive one raises, contours use the real field, and `mosaic`
relabels a grid of log panels itself while refusing to pool log with linear.
`colorscale=` takes ANY Plotly name (`_r` reverses) or stops, and the static
backend builds its colormap from Plotly's own stops, so a name renders the same
on both (the mirrored-`RdBu` hand table is gone). `features=` draws
wells/ponds/boundaries over the cells: a path, GeoDataFrame, geometry,
`mf.plot.Feature(...)` (`layer`/`color`/`width`/`size`/`fill`/`label`/`hover`/`crs`),
or a list/dict of those — one legend entry and an attribute hover per layer,
drawn by both backends, kept by `mosaic`/`animate`, resolved eagerly so a bad
path raises at the call. Engine: `utils/datatypes/features.py`. **`locs=` is
DEPRECATED** (warns, drawn as a layer named "locations") and is named on no
signature — do not reintroduce it.
**The COMPOSERS name the map options too** (ledger 258, 2026-09-14).
`mosaic`/`animate` are defined ONCE on `SpatialView`, and they now name all 23
per-panel drawing options plus `backend`, with descriptions spliced from
`plot.map` — so `packages.npf.k.mosaic(logscale=True, features=..., show_layer_elevs=True)`
completes in an editor. Only options the caller actually SET are forwarded
(`_explicit_options`), and they are REFUSED on `kind="plot"`/`"section"`. They
rotted because `test_noun_signatures.py` discovers subjects by a body mentioning
`plot.map(`/`_choropleth_factory`, which a composer's never does; it sweeps them
explicitly now.
**`show_layer_elevs`/`hover_surfaces` are Tier 1** (same ledger): the layer frame
is the GRID's, so they act on a K or conductance map, not just on heads. They
were inert everywhere but heads because `Choro._build_hover_context` fetched the
elevations inside its heads branch and every noun passes `hover_heads=False`;
and `HoverSpec._effective_blocks` answered a surfaces-only request with a
`LayerTable(fields=("__surface__",))` whose key nothing populated. Both fixed.
The elevations are still fetched ONLY when a map asks — `gdf_topbtm` is lazy.
**Never give a picture verb a bare `**kwargs`** (8.8, 2026-08-24): PyCharm and
Pylance are STATIC — they read the `def` line and never run the module, so
`__doc__`, `__signature__` and `functools.wraps` reach `help()` and reach no
editor. Name the parameters, mirroring the default of whichever link in the
forwarding chain owns each one, and keep a `**kwargs` tail only where it is
genuinely open (Plotly trace names). The bound namespace method repeats the free
verb's list; three tests fail by name if a default drifts, a bound copy falls
behind, or a parameter is documented-but-unaccepted (or the reverse). A narrower
scope may offer FEWER parameters, never other ones. The one exemption is the
namespace-level `field=` sugar, which dispatches to accessors with incompatible
signatures — it stays `*args, **kwargs` and points at the leaf.
**This reaches the NOUN tier too, and as of 2026-09-02 all 15 are converted**
(`tests/test_noun_signatures.py`, which DISCOVERS its subjects by AST sweep
rather than listing them; its `UNCONVERTED` set is now empty and asserted empty).
A noun names Tier 1 (`NOUN_MAP_PARAMS`, 22 including `backend`) plus its own
selectors, REFUSES `values`/`type`/the legacy hover trio, and omits the Tier 2
per-layer names that were measured byte-identical on a record noun. **Docstrings are STATIC, and generated** (ledger 178/180, 2026-09-05): the free
verb is the single author, and `scripts/derive_docstrings.py` writes the spliced
result back into each source file, because PyCharm and Pylance never run the
module and were showing ONE line of `model.plot.section` against 169. Change a
reference on `plot.map`, rerun the script, review the diff.
**Text between `.. generated from myflopy.plot.<verb> ...` and `.. end
generated` is NOT yours** -- it is replaced wholesale on the next run. Everything
outside a fence is the method's own: edit it in place, no regeneration needed.
The fences are what let the splice recover its own input (`strip_generated`), so
idempotency and the ratchet are structural rather than guarded -- the first
version cached the pre-splice text instead, and once its output was written back
the cache held that, every inherited entry was shadowed, and **all 29 methods
silently stopped propagating while `--check` stayed green**. Editing a reference
without regenerating fails; so does hand-editing inside a fence. A name a method
documents that its reference ALSO documents is a *shadow* -- permanent, because
the method's entry wins -- and every one is pinned in
`tests/docstring_shadows.json`; read that diff, it is the list of places your
edit will not reach.
**Converting one? The hazard is a `**kwargs` tail something READ.** PRT's `map`
consulted `base_kwargs` for three separate refusals; naming the parameters
emptied it and would have disabled all three silently (`_explicit_options`
reconstructs it from the signature's own defaults). And `GridPlots.map` /
`StackPlots.section` are NOT nouns — they are verb scopes, where `values=` and
`custom_hover=` are legitimate; `VERB_SCOPES` splits them out.
Runnable tour: `examples/mf6/notebooks/plotting_vocabulary_tour.ipynb`.
RETIRED, do not reintroduce: `model.cor`, `model.srf`, `model.visualize`, the
eleven `vor.*` plotting aliases, `stack.preview/views/vtk_3d/surface_3d`.

Only the *spatial* half of this was documented before 2026-07-18; derived tables
had no rule, and the gap produced `sfr.results.long_profile`/`plot_long_profile`
plus ~20 lines of hand-rolled matplotlib in two notebooks that redrew — off-color —
a figure the library already built.

## Browsing a Project: libraries, runs, forks (2026-09-04/06, ledger 179-191)
**`docs/project_workflow.md` is the walkthrough**; `docs/project_libraries_design.md`
is the reasoning record. `project.packages` / `.grids` / `.simulations` are `dict`
**subclasses** — every dict operation works unchanged, plus `names()`, `summary()`,
a did-you-mean `KeyError`, and `_repr_html_`. `project.package(key)` / `.grid(key)`
mirror the long-standing `simulation(name)`; `package_spec(key)` had existed on the
retired `ProjectSpec` and was dropped in the 2026-06-17 fold, unnoticed because
**`api_snapshot.json` covered no part of the `Project` surface** — it does now
(`Project.members`, `Run.members`, one row per library), so a verb vanishing fails CI.
**Variants: `packages.fork(source, key, **overrides)` / `packages.edit(key, ...)` /
`simulations.fork(source, name, swap={...})`.** Forks are **built, not patched** —
`copy.deepcopy` reconstructs a frozen+slots dataclass through `__reduce_ex__`, so
`__post_init__` never runs and an in-place option write accepts a typo silently;
they go through `dataclasses.replace` instead, and deep-copy by default so editing a
fork cannot write through into the parent's arrays. A rename carries into
`pname`/`filename` and REFUSES a `*_filerecord` still naming the old package. A flat
`swap=` reaches every model carrying the label; `{"model": {...}}` names one.
**Runs are a ledger, not a list: `project.runs`.** It includes bare MF6 workspaces
(`discovered`) and truncated manifests (`unreadable`) — measured 39 workspaces
against 14 manifests under `examples/`. **`status` is myflopy's verb, `outcome` is
MF6's verdict** read from `mfsim.lst`: `normal`/`failed`/`incomplete`/`absent`/
`superseded`, the last meaning the listing predates this run's own preparation
(`prepare_run(overwrite=True)` reuses the directory — design D5, still open, ledger
191). `run.reconcile()` adopts a verdict; it will not adopt a superseded one.
**`project.runs.drift()`** says which library packages a run used have changed since,
via a content fingerprint stored per run. A reopened run BUILDS from its own copy
of those packages (`<run>/_library/`, ledger 230), so drift says the library moved
on, not that the run can no longer be rebuilt as it ran. `project.save()` returns a **`SaveReport`**,
not a path, and persists what it can rather than refusing the whole project.
Runnable tour: `examples/mf6/notebooks/project_libraries_tour.ipynb` (executed by
`tests/test_project_libraries_tour.py`, not merely source-checked).
Two costs that surprise: fork and digest scale with the number of PYTHON OBJECTS,
not cells — a 72×9,405 array package digests in ~5 ms, the 654,480-record EVT in
**1.3 s**, so `PackageLibrary.digests()` caches until the library is mutated.

## Compromise ledger (standing user rule, 2026-07-17)
**`docs/compromises_and_deferrals.md`** records every deliberate scope cut,
test-fidelity trade, judgment call, or considered-but-omitted optional inside
delivered work. Any change that makes such a call MUST add/update its ledger
entry **in the same pass** (it is part of the docs-always-in-sync rule).
Entries are removed only when the compromise is actually undone.
**An entry records the DECISION, not the debugging** -- what was deliberately not done,
the measurement that settled a judgment call, the trap a later reader would re-discover.
Bug narrative belongs in the commit message; entries doubled in length (22.1 lines each
for 1-132, 50.9 for 133-264) and that growth is what they filled with. The labelled
`What:`/`Why:` format the preamble used to prescribe is carried by 42 of entries 1-132
and **none** of 133-264.
**A number is an address**: `ledger N` is cited **689 times over tracked files** --
`docs/` 317, `tests/` 173, `src/` 138, this file 59, and none at all in the notebooks --
naming 159 distinct entries, so entries are never renumbered and a number is never reused
(a resolved entry is restated below, struck through, keeping its number). Adding an entry
means rerunning `python scripts/derive_ledger_index.py`;
`tests/test_ledger_format.py` fails if the index is stale, if the numbering breaks, if a
title never closes its bold, or if any citation names an entry that does not exist. Its
sweep reads **git-tracked files only** -- an `rglob` walk also reads the stale repo copies
under `.claude/worktrees/`, which were 62% of the occurrences and made the count differ
per machine.

## Printing? Only if a human asked for it (plan 7.2, 2026-08-01)
Library code logs; it does not print. Use `myflopy._logging.get_logger(__name__)` —
`logger.info` for progress on slow work, `logger.warning` for a degraded or ignored
input, `logger.debug` for per-row detail. **The one exemption: output whose job is to
report something a human explicitly asked for.** In practice that means output behind a
flag the caller passed (`verbose=`, `progress=`, `verbosity_level`) — the flag IS the
human asking — plus `runtime.run_simulation`'s "Success is:" line.
`tests/test_no_library_prints.py` enforces it: gated prints are exempt by AST
inspection without being listed, everything else must be in an allowlist exact in both
directions, and `import myflopy` must be silent.

**Also don't call flopy APIs it has deprecated.** `tests/test_no_deprecated_flopy_calls.py`
names them (`gwf.package_names`, `gwf.package_name_dict` — use `get_package_list()`, or
myflopy's own `model.package_names` -- the real names, one per package;
`model.package_types` lists each type once, ledger 228). Note `getattr(gwf, "package_name_dict", {})` still
warns: reading the attribute is what warns, so the default only covers absence.

## Catching exceptions? Name them (plan 7.3, 2026-07-31)
**No bare `except:` anywhere** — it swallows `KeyboardInterrupt`. **`except Exception`
needs a `# noqa: BLE001` and a comment saying why the set cannot be closed**; eighteen
survive on that basis (flopy's `utils/voronoi.py` has a literal `raise Exception(...)`;
`pickle.dump` walks an unbounded graph; `__dir__`/`__repr__` must never raise).
`tests/test_exception_narrowness.py` pins the allowlist EXACTLY — adding a broad catch
fails, and so does closing one without recording it.

Every narrowed handler logs the swallow at DEBUG via `myflopy._logging.get_logger`; the
message convention (name the DEGRADATION, not just the error) is in that module's
docstring. Two lessons worth carrying: these `try` blocks usually wrap a **pipeline**,
so the honest tuple is wider than it first looks (a lazy `MFSimulation.load` reaches
flopy's `MFDataException`/`FlopyException`, which subclass `Exception` directly); and a
broad handler around a block that **validates its own arguments** will eat the
validation — that bug was found twice, in `read_gpkg` and `contour_line_segments`.

## Test suite (fast by design)
- Full suite (**3472 passed / 4 skipped**, 2026-09-15): `pytest -n 10` ≈ **85 s** (worksteal dist is in
  addopts); serial (`-n0`) ≈ **6m28s**; inner loop `pytest -m "not slow"` ≈ 32 s.
  **A parallel run competing with another MF6 workload fails intermittently** —
  a different test each time, all passing in isolation, wall time 85 s → 245 s.
  That is contention, not a regression; re-run alone before believing a failure.
  Serial is ~3.7x the parallel run — that gap is the price of the `-n0`
  pre-commit check below, not a regression.
  **Verify sign/column changes with `-n0`** — a session-fixture/xdist interaction
  can report green while serial catches real failures (see ledger 48).
- Tests share ONE session-scoped canonical model on
  `CanonicalModelConfig.testing()` (21×21, smallest contract-complete profile —
  all 15 packages + 5 obs families). Weekly CI re-runs everything on the 50×50
  profile via `SIMPLE_MODFLOW_CANONICAL_PROFILE=validation`.
- Tour/verify interactively: `examples/mf6/notebooks/canonical_fast_tour.ipynb`.
- **A test grid is an ASSUMPTION, and the one difference that matters is the vertex
  DTYPE** (ledger 263, corrected in 265). `rectangular_voronoi` is a real SciPy Voronoi
  with ragged, unequal cells -- what it cannot show is `numpy.int64`: its cell2d vertex
  indices are all python `int` (82 of 82 at 4x4), where a `GridSpec.voronoi` grid mixes
  them (114 `int64` + 36 `int`). `_json_value` passed `np.float64` (a `float` subclass)
  while refusing `np.int64` (not an `int` subclass), so real models could not be
  serialized, frozen or rebuilt while 3,315 tests passed. Use the session-scoped
  **`triangulated_grid`** fixture in `tests/conftest.py` for that class of test -- and
  note that **raggedness is NOT the discriminator**, so asserting it proves nothing.
  Stronger still, and fixture-independent: push the numpy scalar types through the
  boundary directly (`tests/test_real_grid_round_trip.py::test_a_numpy_scalar_serializes`).
- **A numpy scalar at a serialization or hashing boundary is a recurring bug class**
  (ledger 265). `tests/test_numpy_at_the_boundary.py` is the guard; most of it needs no
  grid, which is why it finds what a grid fixture cannot.
  The trap is numpy's inconsistent subclassing: `np.float64` IS a `float` and `np.str_`
  IS a `str`, so both slip through an `isinstance` passthrough unconverted, while
  `np.int64` and `np.bool_` subclass nothing and raise. Converting only the half that
  raises leaves the other half in the output -- `json.dumps` accepts it and **PyYAML,
  which dispatches on the exact type, does not**, so `to_yaml` broke for every model
  carrying a DISV. **Test with `type(x) is`, never `isinstance`**, or the check passes on
  the unconverted value. Same class, three more: `_canonical` must NORMALIZE then fall
  through (returning early skipped the string interning, so a `np.str_` boundname read as
  EDITED through `runs.drift()`); it and `_json_value` must agree about a 0-d `ndarray`;
  and `.item()` on `np.longdouble` returns its own type, so recursing never ends.

## What this project is
A Python-first MODFLOW 6 toolkit built on top of FloPy. Key strengths:
- **Voronoi/unstructured grids** (DISV/DISU) via `VoronoiGridPlus` + `TriangleGrid`
- **Declarative spec API** — `SimulationSpec`, `ModelSpec`, `PackageSpec`, `GridSpec` dataclasses in `specs.py`
- **Package-first API** — `gwf()`, `gwt()`, `gwe()`, `prt()`, `lak()`, `sfr()`, etc. in `package_api.py`
- **Multi-physics** — GWF, GWT (transport), GWE (energy), PRT (particle tracking) with exchange packages
- **Rich interactive visualization** — Plotly choropleth maps, HTML sliders, cross-sections, animations
- **Workspace/project management** — `Project`, `Run`, `load_run` in `workspace.py`
- **Parallel model workflows** — `ParallelModelWorkflow` in `parallel.py`

## Canonical model-building API: package-first (use this for new work)
**Starting a new model? `docs/model_building_cheatsheet.md`** is the ordered path —
contours → `Surface` → `LayerStack` → `ModelContext` → packages → `Project`/run →
results — every call verified end to end 2026-08-26, runnable as
`examples/mf6/contours_to_model.py`. The two library gaps it used to document
rather than fix were both **closed 2026-08-26**: GRASS launcher discovery globbed
`grass*.bat` only, so `mf.Contours` needed `GRASS_BIN` plus a hand-set
`PYTHONPATH` on Linux — `contour_interp.py` now falls back to `shutil.which` on
non-Windows and asks the resolved launcher for its bindings path (ledger 146);
and `Surface.maximum/minimum/clamp` are classmethods, so a misbound `a.maximum(b)`
dropped `a` silently — they now raise `TypeError` naming `Surface.maximum(a, b)`
and the fluent `a.floored_at(b)` (ledger 147).
**`docs/package_api_reference.md` is the human-readable map of the whole API** — the
`mf.*` build helpers and the `model.packages.<pkg>.<inputs|results>.<noun>.<verb>()`
read grammar (pinned by `tests/api_snapshot.json` + `package_registry.py`).
**The package-first API in `package_api.py` is THE preferred, canonical way to build models.**
Prefer it for all new model-building work; do not reach for the legacy OO path below unless
you have a specific reason. The OO builder classes (`SFRBuilder`, `LAKBuilder`, `UZFBuilder`,
`RCHBuilder`, `EVTBuilder`, `MVRBuilder`, `GHBFromVector`, `Recharge`, …) are the **engine underneath** the
package-first facade, not a competing API — e.g. `mf.uzf(...)` literally calls
`UZFBuilder(...).build()`. (The canonical model in `canonical_example.py` is written in the
older imperative builder style for historical/computed-cell reasons; that does NOT make the
imperative path "more proven" — package-first is the same engine and is the one to teach.)

### How it fits together (the mental model)
```
Project            ← durable workspace + run/scenario lifecycle (workspace.py); holds NO geometry
  └─ SimulationSpec ← one MF6 simulation: tdis + ims solver + the model(s)
       └─ ModelSpec = mf.gwf(name, context=…, packages=[…])   ← one model
            ├─ ModelContext(grid=, domain=, surfaces=, dates=)  ← geometry; rides on the MODEL, not the project
            └─ packages: mf.disv, mf.npf/ic/sto/oc,
                         mf.chd/ghb/drn/riv/wel/rch/evt (+ .gpkg / .flopy),
                         mf.uzf/sfr/lak (+ .flopy), mf.mvr
```
- **`Project(root, name=)`** (`workspace.py`): `add_grid/add_package/add_simulation` (reusable
  libraries), `prepare_run(name, sim)` → `Run` (built in memory; inspect `run.model(...)`),
  `run.execute()` writes + runs MF6. The Project is the lifecycle wrapper; geometry lives on the model.
- **`ModelContext`** (`specs.py`) attaches to the **model** via `mf.gwf(..., context=ctx)`, NOT to
  the project. It carries `grid`, `domain` (idomain), `surfaces`, `dates`. The built model exposes
  it as `model.myflopy_context`. The GIS-aware package helpers (`mf.uzf`, `mf.sfr`, `mf.X.gpkg`)
  take `context=` so they can map features/cells onto the grid.
- **Moving a model is `sim.with_grid(grid, on=, name=, grid_free=())`** (ledger 271,
  2026-09-15). It drops everything derived from the old mesh and CLASSIFIES every
  package -- `rebuilt` (carries a recipe), `stale` (pinned elsewhere, or a
  hand-written DISV, which IS the old mesh), `unverifiable` (names cells with no pin),
  `unchanged`. **It reports and stamps; the BUILD refuses**, so read
  `sim.move_report()` before `prepare_run`, and `run.resolution` after. `on=None`
  moves EVERY model with a grid, because coupled models must share a mesh (ledger
  268); `on=` is the spelling that can strand a sibling. `grid_free=` is the caller
  asserting mesh-independence. The stamp is the only thing behind a hand-typed cellid
  across a move -- every number stays in range on a larger mesh. **A DEFERRED-grid
  model has no old pin to stamp**, so its hand-written cellids move unrefused; they
  were never verifiable anyway. `replace_grid` is a warned `__getattr__`-only alias.
- **A grid swap re-maps NOTHING outside the discretization** (ledger 249,
  2026-09-13). Every `.gpkg` helper resolves its cells at DECLARATION, so
  moving the grid leaves every package's cellids,
  connection data and dimension counts on the old mesh. Nothing looked, and MF6 ran the
  old cell numbers against the new geometry under **Normal termination** -- a lazy CHD
  raised `IndexError` on a finer grid and resolved to ZERO cells on a coarser one.
  `ModelSpec.build` now compares the built DISV's `ncpl`/`nvert`/cell centres, the
  context `domain` (cells AND layers), and every `stress_period_data`
  cellid against the grid (NOT `surfaces` -- that check was removed during review,
  because `surfaces` has no reliable "rows are cells" invariant), and raises listing all of it at once. **The cell-centre test
  is index-aligned after removing each cloud's own centroid**, with a tolerance of a
  thousandth of a cell taken from the data. Three traps, each of which it got wrong
  first: a Voronoi cell's DISV centre is **not** its polygon centroid (38.9 length units
  apart on one 4x4 grid -- compare against the grid's own `get_disv_gridprops()`); a
  converted model's frame is deliberately not the world's (1.34e+06 apart on the USG
  fixture, the State Plane offset), so the test must be translation-invariant; and
  **sorting to gain rotation-invariance threw away cell ORDER**, which is the only thing
  a cellid means -- a renumbered mesh moved the mean cell 509 length units and built
  clean. An absolute tolerance also refused the same mesh written to a hundredth of a
  foot, and float32 at State Plane. **It is a range check, not a proof**: a stale cellid
  that still falls inside a LARGER new grid passes, which the grid pin closes; advanced
  packages are not cellid-checked because the cellid's position differs per package
  (SFR's `packagedata[1]` is a cellid, LAK's is a stage). To move a model to a new grid
  today, re-declare its GIS packages; lazy resolution is planned work.
  **A resolved package also carries a PIN** (ledger 250, 2026-09-13) -- `grid_pin(grid)`,
  a hash of the mesh's CONNECTIVITY (the cell count and each cell's canonical vertex
  ring; NOT the vertex count, which FloPy computes at write time) and nothing continuous. That closes the range check's gap: a
  package resolved on a small grid and rebuilt on a larger one has every cellid in range
  and every one naming a different cell. Being connectivity-only makes it exact where
  geometry needs a tolerance -- the same mesh at fewer decimals, in float32, or shifted
  to a local origin pins the SAME, a renumbering does not, and a mesh stretched without
  renumbering pins the same (the geometry test catches that; the two are complementary).
  It rides in `metadata`, so it is serialized and forked but **excluded from
  `package_digest`** -- the records already encode the grid, and digesting the pin would
  report every stored run as drifted on the day it landed. Stamped once in
  `GeoPackageSource._metadata`, which all seven list BCs reach. **A hand-written package
  stays unpinned**: you cannot pin cell numbers the caller typed.
  **Every grid-bound helper the library resolves itself now pins** (ledger 267):
  the seven `.gpkg` forms, `rch`/`evt` in BOTH list and array shapes, `uzf`, `sfr`,
  `lak` and `hfb`. `lak` and `hfb` are the two the range check is blind to in the
  direction that matters -- LAK's cellids live in `connectiondata`, which `_cellids_in`
  deliberately will not walk, and a barrier's pairs are all still in range on a LARGER
  mesh -- and both BUILT CLEAN across 28 -> 92 before the pin. **The array form pins
  CONDITIONALLY, on whether the resolved payload names cells**: `mf.rch.array(irch=0)`
  means the same layer on every mesh and must stay movable, while the default
  `irch="top_active"` resolves to an `(ncpl,)` array and must not.
  `tests/test_grid_pin.py` BUILDS each helper and reads its metadata, in both
  directions -- reading the source is what missed `rch.array` twice.
- **A run keeps a FROZEN SNAPSHOT, so it rebuilds without its sources** (ledger 251,
  2026-09-13). `<run>/.myflopy/frozen.npz`, written at the end of `Run.build`: the
  resolved spec as JSON, every array beside it, plus the grid's vertices, cells, CRS
  and layer elevations. It exists because a recipe cannot hold an array -- `to_dict`
  refuses `ndarray`, so a model with an inline `k` wrote `spec_recipe=None` and could
  not be rebuilt AT ALL. Measured on a real Ten Trails run: **2.78 MB against a 1,023
  MB run, 0.27%** -- binary beats MF6's ~15 characters per number. Unconditional.
  **The snapshot WINS over the manifest's recipe on reopen**, and it holds the grid the
  build RESOLVED: a declared `ModelSpec` carries the grid RECIPE, so freezing the spec
  captured no grid for a DEFERRED one and a reopened run went back to its boundary
  GeoPackage and died there. Prepare the run again to re-resolve.
  **ONE FILE, not a directory**, because `os.replace` of a file is atomic and a
  directory needed two steps and a shared staging name. **Stamped with its run name
  and discarded by `prepare_run` when a directory is reused** -- without that, a
  re-prepared run thawed and rebuilt the PREVIOUS simulation in silence.
  **An object-dtype array is REFUSED at freeze time**, naming the package and the
  offending fields: `savez` stores one by pickling and `load` refuses a pickle, and
  FloPy's own list-BC recarray is that shape (a cellid is a tuple, a boundname a
  string), so writing it produced a snapshot nothing could ever read. Not carried:
  `ModelContext.dates`.
- **`project.save(bundle_sources=True)` makes a project an ARCHIVE** (ledger 256,
  2026-09-14) -- it copies every file the project's GRIDS are built from into
  `inputs/`, repoints the saved recipes AND the live registry (a later plain `save()`
  would otherwise write the original paths back), and reports the pairs in
  `SaveReport.bundled`. **Grids only, measured not assumed**: package sources resolve
  EAGERLY, so a saved package holds records rather than a path -- zero file reads
  inside `PackageSpec.build` across 306 tests, a project rebuilds with its GeoPackage
  deleted, and `metadata['source_path']` is read by nothing. A `GridSpec` is the
  opposite: `ModelSpec.build` re-resolves it on EVERY build. `rasters=` and
  `idomain_path=` are bundled too even though they are read later still (the first
  `gdf_topbtm`, the first edge query). An `external=True` source IS bundled and the
  flag cleared -- it meant "never anchor", which stops being true once a copy is
  inside. **A run now records its project root**, because bundling stores RELATIVE
  paths and `Run.reopen`/`load_run` had no `project_root`, so a relative source
  resolved against the current working directory.
- **A saved project keeps a 16-character GRID FINGERPRINT, because it cannot keep the
  grid** (ledger 257, 2026-09-14). `ModelSpec.to_dict` cannot carry `context.grid`,
  `domain` or `surfaces` -- inlining them was refused in ledger 193 (floods
  `Project.validate`, breaks `to_yaml`, ~3.3 MB per model at 9,405 cells) -- and it
  dropped them SILENTLY, which undid both new guards: a model that raised live BUILT
  CLEAN after `save()` + `load()`, and a package pin survived the save and was never
  compared. `to_dict` records `grid_pin` into the context metadata, and
  `disv_pin(disv)` takes the same fingerprint off a BUILT package, so the guard still
  compares a reloaded model against its mesh. **What the fingerprint cannot see:**
  connectivity-only means a rescale, a mirror or a unit conversion passes -- the save
  warning and the refusal both say so. `SaveReport.context_dropped` names what was
  lost. **To keep the FULL check across a save**, register the grid and use
  `mf.gwf(..., grid=mf.grid_ref(key))`, which is serialized on the model and resolved
  at build time.
- **Registering an existing key RAISES** (ledger 208, 2026-09-10):
  `add_package`/`add_grid`/`add_simulation` take `overwrite=False`. They used to
  replace silently while `packages.fork` refused the same collision — the guarded
  verb was the deliberate one and the unguarded verb was the accident (a second
  script run, a cell executed twice). **`packages.edit(key, ...)` stays unguarded**:
  changing what a key MEANS is supported, and `runs.drift()` shows the cost.
- **`mf.disv(layers)`** (ledger 206, 2026-09-10) — the discretization was the ONE
  package you could not get from a `mf.<name>(...)` helper; it came only from
  `result.to_disv()`. It delegates to that same method rather than re-deriving,
  because a second resolution of the same stack measured 8.5 ft of drift. Mixing
  it with hand-written `nlay=`/`vertices=` raises. The hand-written form still works.
- **One verb per operation, both tiers** (ledger 204, 2026-09-10):
  `with_packages(*pkgs, slot=, on=, name=)` / `without_packages(*names, on=)`. One or
  many, varargs or a list. GONE: `ModelSpec.with_package`/`without_package`,
  `SimulationSpec.with_package`/`add_package`/`replace_package` — four spellings for
  add-or-replace, of which two differed only in the absent-slot case and one
  recorded its own lineage as `'replace_package'`. **`slot=` carries the must-exist
  assertion** `replace_package` held, and that is load-bearing: without it a missed
  swap ADDS, and MF6 sums two RCH packages under normal termination (0.8000 against
  an intended 0.4000). `on=<model>` picks the tier; `name=` derives a named variant
  in one step. **`without_packages` at the simulation tier is NEW** — there was no
  way to drop a tdis or solver short of rebuilding the tuple.
  Lineage now reads `{"operation": "with_packages", "model":, "packages": [...]}`.
- **Coupled models must share a mesh, and that is now REFUSED at build** (ledger 268,
  2026-09-15). MF6 requires it for `gwfgwt`/`gwfgwe`/`gwfprt` — `exg-gwfgwt.f90:224`
  terminates on a node-count mismatch and `:239-255` on an idomain one — but it compares
  **counts and idomain only, never geometry**, so a mesh of the SAME size passes its
  checks and it maps node-for-node onto the wrong one. Measured: GWF on a 3x2 Voronoi,
  GWT on a 2x3, both `ncpl=6`, built clean with zero warnings. **`gwfgwf` is excluded** —
  a coarse/fine pair is what that exchange is for. DISV only: a structured pair is not
  compared, and MF6's node count is all that stands behind it.
- **Coupling models is `mf.exchange(a, b)`** (ledger 203, 2026-09-10) — it reads the
  two kinds and picks the MF6 class; `.gwfgwt/.gwfgwe/.gwfprt/.gwfgwf` are explicit
  and take plain names. The four `build_*_exchange` functions resolved on `mf` but
  were in NEITHER `__all__` NOR `dir(mf)`, so autocomplete never offered the one
  name a coupled model needs — and the spelling four docstrings taught raised
  ValueError. **Nothing wires an exchange for you**: a coupled sim without one
  builds, runs and moves nothing under normal termination.
- **Simulation assembly is ONE list** (ledger 202, 2026-09-10):
  `mf.simulation(*models, name=, packages=[tdis, ims, ims2], exchanges=)`. `tdis=`,
  `solver=` and `complexity=` are GONE — `complexity=` was ignored the moment
  `solver=` was passed. myflopy fills what the list omits; a solver with no
  `models=` claims the rest, so one model needs `packages=[mf.ims(complexity=...)]`
  and no names. Solvers are re-ordered by model declaration (MF6 demands GWF's
  solution first). A model package in `packages=` RAISES.
- **A `.gpkg` source that lands nowhere RAISES** (ledger 211, 2026-09-11). Every
  feature off the grid, on inactive cells, or only bordering cells used to build an
  EMPTY package -- MF6 rejects it naming the file, or under `maxbound=` runs without
  it. The error prints both extents and both CRSs, because a location mismatch is two
  boxes; a partial miss WARNS by name. A polygon covers a cell by AREA (ledger 210):
  drawn along faces it reaches only the cells inside it, never its neighbours. A LINE
  piece under the same 0.01 tolerance is noise too (ledger 237) -- and this is the
  SHARED rule, so it moves `map(select=<line>)`, a region and a line head target as
  well as a `.gpkg` source: a line ending on a face and overshooting by 1e-6 ft -- a CRS
  round trip moves an end by 3e-7 -- wrote a full second record, where `mf.sfr` had
  dropped the reach since ledger 214.
- **Sampling a raster onto cells: `gaps=` says what a MISSING cell means** (ledger
  248, 2026-09-13). `get_raster_vals_at_centroids` handed back NaN for a cell off
  the raster's data with nothing said, and a NaN in `strt` or `gdf_topbtm` makes MF6
  report success and return a NaN budget. Default `"warn"` names the cells;
  `"nearest"` takes each from its nearest valued cell; `"raise"` refuses; `"ignore"`
  is for a caller whose next stage fills them, which is why `Surface.sample` passes
  it -- a partial raster is NORMAL for a layer surface (`fill="propagate"`,
  pinch-out), and warning there would have drowned the default.
  **The motivating trap: `InterpolatedSurface` writes a hull-limited raster.**
  `griddata` fills only the CONVEX HULL of its input points, and those points are the
  source grid's CELL CENTROIDS -- so a raster interpolated from one model is blank in
  a rim half a cell wide around the whole footprint, and a different grid's boundary
  cells land in it. For remapping between grids the real answer is an area overlay on
  the source cells' POLYGONS (Phase 2), not a raster round trip. A raster covering NO
  cell raises with both extents and both CRSs -- but only when they are actually
  disjoint; all-no-data where they overlap is a blank raster, not a misplaced one.
  `vor.get_raster_vals_at_centroids` forwards `method=`/`nodata=`/`gaps=` now; it
  forwarded none of them.
- **A `.gpkg` NULL raises too, naming the row** (ledger 247, 2026-09-13). A null
  attribute reached a written package as a NON-FINITE number by two separate
  routes: `k_array` did `float(row[value])`, which makes NaN, and a boundary record
  kept the `None`, which numpy casts to NaN when FloPy builds the recarray. MF6 then
  answers a NaN budget while reporting normal termination. Every resolved value is
  checked now -- one funnel for chd/ghb/drn/riv/wel/rch/evt and `k_array` -- and the
  error names the feature, the field and what it found (a NaN offset column on a
  `CellSurfaceOffset` names the COLUMN, not the surface). The legacy `KFromVector`
  filled a null with the layer default; substituting a number the user never gave is
  its own silent guess, so this refuses instead. A feature that ALSO lands nowhere
  reports the miss (ledger 211), which explains more.
- **`edges_only=` is the edge of the ACTIVE domain, per layer** (ledger 246,
  2026-09-13). `_edge_set` called `get_grid_edge()` with no argument, so only a grid
  carrying its own `idomain_path` was ever honoured -- a `ModelContext(domain=...)`
  array, the documented route and what `LayerStack` produces, never reached the edge
  test, though the same class already honours it per cell. On a 5x5 grid with one
  inactive cell the package wrote 16 cells where the active edge is 20, dropping the
  hole's four face-neighbours; with the outer ring inactive it refused the source
  outright. The context's domain WINS over the grid's file (ledger 209's precedence);
  with none given the file still decides.
- **A lake COVERS its cells, by area** (ledger 241, 2026-09-12) -- `mf.lak` kept
  every cell its polygon grazed and had its own activity test. A one-cell lake
  round-tripped 2927 -> 4326 -> 2927 became 4 cells (extras overlapping by 7.41e-6,
  9.07e-14 and 7.41e-6 ft2) and one grown by 1e-6 ft became 9, each a full LAK cell
  with its own vertical connection. `LAKBuilder.lake_cells` is `cell_shares` now and
  `_active` is `domain_active` -- a 1-D idomain masks the lake per cell where it was
  ignored, and a layer past the domain names itself. A lake on no cell follows 211.
- **A PRT release point goes in the cell ACTIVE in its layer** (ledger 242,
  2026-09-12). `PRTReleasePoints.from_points` located points with flopy's
  `modelgrid.intersect`, which answers with the lowest CELL2D and reads no idomain:
  with cells 0 and 4 inactive in layer 0, a point on the 4|5 face resolved to cell 4
  and one on the 0/1/4/5 vertex to cell 0. Nothing silent came of it -- flopy's PRP
  checks every cellid, so `PRTProject` died with `Cellid (0, 4) is outside of the
  active model grid`, an id the caller never wrote, where the active neighbour the
  point also touches releases and tracks. It is `cell_shares` + `domain_active` now
  -- the same rule as a head target -- so a tie goes to the lowest cell ACTIVE in
  `layer=`; a point off the grid raises with both extents and both CRSs while a
  partial miss warns by name (ledger 211, where flopy's `intersect` named only the
  point and ended the whole call), and a frame in another CRS is reprojected first.
  What is WRITTEN is still the model's own frame (ledger 218). **Both halves are
  DISV-only**: a structured (DIS) grid keeps flopy's `intersect` AND flopy's
  all-or-nothing refusal, because the shared rule is indexed by cell id and a DIS
  cell is `(row, col)`. `from_cells` trusts the ids it is given. One sindex query
  per point also replaces flopy's per-point scan of every cell -- 200 points on a
  10,000-cell grid, 7.5 s -> 0.016 s.
- **A lakebed from bathymetry is `mf.lak(lake_bottom=<Surface>)`** (ledger 231,
  2026-09-11) -- `mf.Contours`/`mf.Raster`/algebra, a raster `Path` or `str`, or a
  cell-indexed Series. Sampled ONCE (it was one full-grid sample per cell read).
  `bathy` reads the bed at the lake's NEIGHBOURS too; where the source has none they
  are the model top, and a per-cell bed above its model top is capped with a warning.
  A scalar bed is never capped -- `rectangular` facilities keep their flat bottom.
- **Package-first surface**: `mf.gwf/gwt/gwe/prt`; `mf.disv` (+ `mf.dis`/`mf.disu`
  structured/unstructured passthroughs, viz stays DISV-only); `mf.ic/npf/sto/oc/tdis/ims`;
  list BCs `mf.chd/ghb/drn/riv/wel/rch/evt` each with `()` (direct data), `.gpkg(path, context=, nper=)`
  (from GeoPackage), `.flopy(...)` (raw FloPy escape hatch); advanced `mf.uzf/sfr/lak` (`()` =
  high-level builder, `.flopy(...)` = raw); `mf.mvr` with `mf.Move(mf.MoverConnection("sfr",0),
  mf.MoverConnection("lak",0))`. **MVR is validated**: moved packages must be declared in the
  model AND ordered before the mover (`test_advanced_specs`). PRT is spec-first
  declarable too (§6.3A): `mf.mip/prp` + `mf.ems` (PRT is explicit — EMS, never IMS;
  `mf.simulation` defaults kind-aware) + prt6-dispatched `mf.disv/dis/oc`; coupled
  GWF+PRT needs NO FMI (the exchange passes flows) — FMI/grid copying is only for
  `PRTProject`'s post-hoc separate-simulation path.
- **Layers** (same facade/engine pattern): use the facade **`mf.LayerStack`** (`layers.py`) —
  `LayerStack(vor, top=Raster("ground.tif")).add("sand", thickness=20, pinch="inactive")
  .add("clay", bottom=Contours(...)).build()` → disv-ready top/botm/idomain (plus `.qc()`,
  **`vor` is OPTIONAL since 2026-08-27** — `LayerStack(top=...)` declares the layering
  before any grid exists, and `build(vor)`/`qc(vor)`/`to_disv(vor)` take it later, or
  `for_grid(vor)` returns a bound copy (ledger 148). That is what lets the project-first
  order (declare layering → build grid) keep the per-layer `thickness=`/`min_thickness`/
  `pinch` control that only the facade has; `LayerSurfaces` has one global rule.
  **A UNIT is not a LAYER as of 2026-08-28 (ledger 158).** `.add(name, ..., split=3)`
  or `split=[0.3, 0.7]` cuts one geologic unit into N model layers `name_1..name_N`
  WITHOUT moving its base (`names=[...]` names them yourself -- one per layer, needs
  a split of 2+, and `replace` validates `split`/`names` as a PAIR rather than
  dropping names silently; the UNIT key is unaffected either way); `res.units` maps unit → layer indices and
  `res.per_layer({unit: value})` builds the positional lists `mf.npf`/`sto`/`ic`
  want — never hand-write those, a stale list of the RIGHT length after a split is
  accepted silently with new meaning. `min_thickness`/`pinch` stay UNIT-scoped (per
  slice, a present unit comes back with holes in it). Engine side this is ONE new
  primitive, `Surface.toward(target, f)` — general `Surface - Surface` stays
  refused. There is deliberately no `.split()` verb: `replace(name, split=N)`.
  Not splittable: an `Isopach` bottom (measured from above, so cuts would drift).
  **A stack SERIALIZES as of ledger 269** — `Surface`/`LayerStack` carry tagged
  `to_dict`/`from_dict` under `$surface`/`$layer_stack` tokens (NOT a `kind` tag:
  `Surface` has a `kind` field of its own and the source dispatcher keys on `kind`,
  so one tagged as the other comes back a plain dict, silently). The payload is the
  DECLARATION — top, units, `min_thickness`/`pinch`/`split`/`names` — never `vor` and
  never the built arrays. **`Surface.from_array`/`from_points` are REFUSED by name**:
  their content is per-cell data, which ledger 193 refused to inline. And a stack does
  NOT determine its own arrays — `reconcile`/`min_sep`/`trigger_sep`/`method` live on
  the `build()`/`to_disv()` call and move the answer, so a round-tripped stack is a
  recipe that still needs them.
  Also fixed there: `thickness=<Surface>` was silently an ELEVATION, and
  `build`/`to_disv` disagreed on a split stack's idomain.
  `.plot.map()/.plot.section()/.plot.surface()` and `.vtk_3d()`, `from_modflow`).
  **Pictures come from `stack.plot`** as of 8.5a — `preview`/`thickness_map`/
  `cross_section`/`surface_3d`/`views` are gone. `qc()` stays a method: it is a
  report you read, not a picture. It **compiles to** the
  `LayerSurfaces` engine (`surfaces.py`: area-weighted sampling, top-down reconcile, pinch-out
  → idomain), which uses atomic `Surface` objects (`raster`/`from_contours`/`from_points`/
  `from_array`/constant + algebra). Use `LayerStack`; `LayerSurfaces`/`Surface` are the engine/atoms.
  **`vor.gdf_topbtm` is GEOMETRY-FIRST, and that is a contract** (ledger 196, 2026-09-10):
  several consumers decode it POSITIONALLY (`iloc[:, 1]` = model top, `iloc[:, 2:]` =
  bottoms), so where the geometry column sits is load-bearing. `reconcile=True` — the
  DEFAULT — used to append geometry at the end, which made the hover read the first
  bottom as the top (5.0 against a truth of 10.0) and hand back a Polygon as the last
  bottom, silently. Normalization lives on the **getter** (`grid.geometry_first`), not
  the setter, because `__setstate__` and in-place edits both bypass the setter — and the
  getter STORES what it normalized and returns that object, because returning a copy
  loses `vor.gdf_topbtm[2] = ...` silently. Geometry is found by **dtype**, not the name.
  Column *names* (int vs unit strings) were the wrong diagnosis: those raise loudly.
  `ModelContext.surfaces` is a captured reference the getter cannot reach, so SFR, LAK
  and `CellSurfaceOffset` all resolve through ONE function, `specs._layer_elevations`
  (ledger 209): the context wins when given, and a disagreement with the grid's live
  table WARNS. SFR used to read the grid while LAK read the context -- two packages in
  one model on different elevations after a notebook cell re-run (70.0 vs 60.0).

### The discretization can be LAZY (ledger 270, 2026-09-15)
`mf.disv(stack)` on a **gridless** `LayerStack` emits a lazy spec that resolves the
stack against whatever grid the build supplies, publishes the `idomain` and sampled
elevations onto the model context, and records what it built in `BuiltModel.resolved`
so the run freezes **arrays, not the recipe** (a recipe sends a reopened run back to
`gis/ground.tif` — ledger 251's failure). One declaration builds at ncpl 28 and 92.
**This is what makes a LAYERED MODEL SAVEABLE AT ALL**: `mf.disv(result)` holds
resolved arrays, `to_dict` refuses an ndarray (ledger 193), so `Project.save()` used
to persist **zero** simulations for every layered model and `load()` returned an empty
project. `mf.disv(result)` stays EAGER on purpose (24 option-read assertions).
**`mf.disv` forwards all eight of `to_disv`'s build kwargs** — `attach`,
`default_min_thickness`, `default_pinch`, `reconcile`, `min_sep`, `trigger_sep`,
`method`, `refresh` — because **a stack does not determine its own arrays**:
`default_pinch="inactive"` takes idomain layer 0 from 28 to 0, `reconcile=False` moves
`botm.sum` by 26.6. On an already-built `LayerBuildResult` none of them can act, so a
non-default one **raises** rather than being ignored. A declared `domain` that
disagrees with what the stack resolves is refused, naming both shapes.

### Grid: eager (works now) vs deferred (a known seam)
- **Eager** (use for the full GIS stack): build the grid object first (`VoronoiGridPlus`, or
  `mf.GridSpec.voronoi(...).resolve(workspace)`), put it in `ModelContext(grid=vor, domain=idomain)`,
  then declare packages. **Required today** for `mf.uzf/sfr/lak/.gpkg` because they resolve cells
  eagerly at declaration (e.g. `UZFBuilder.build()` bakes resolved data into the `PackageSpec`).
- **Deferred** (`mf.gwf().with_grid(mf.GridSpec.voronoi(boundary=<gpkg>, refinement=<gpkg>))`):
  the project builds the grid at run time into `run.workspace/_grid/<model>` and populates
  `context.grid`. In `ModelSpec.build` the grid resolves BEFORE packages build, and the model
  carries `model.myflopy_context` — so deferred GIS packages are *feasible*, just not wired: the
  helpers would need a grid-lazy mode emitting a `PackageSpec` whose `build(model)` resolves against
  `model.myflopy_context` instead of resolving eagerly. Until then, deferred GridSpec composes only
  with disv + simple BCs, not with `mf.uzf/sfr/lak`.

Reference: `examples/mf6/package_first_full_stack.py` (full package-first stack on a Project).

## Recharge/ET in ARRAY form? `mf.rch.array` / `mf.evt.array` (2026-08-30, ledger 163)
Same physics as `mf.rch(...)`, one array per period instead of one record per cell
— **6.1x faster end to end** on a 72-period 9405-cell model. A fourth method on the
existing helper, NOT `mf.rcha`: MF6 has one recharge package and `READASARRAYS` is
an option inside its file; FloPy's two classes are a FloPy artifact.
**`irch="top_active"` is the default and must stay one.** MF6 with no IRCH applies
recharge to layer 1 unconditionally and SILENTLY skips any column whose layer 1 is
inactive — measured 40% loss under "Normal termination". `irch` is **zero-based**
(FloPy adds one on write; passing 1-based gives `Invalid layer number`).
`nseg > 1` and `boundnames` RAISE rather than being silently dropped — segmented ET
has no array form, because NSEG lives in DIMENSIONS which MF6 never reads under
READASARRAYS. Arrays are opt-in: for a sparse BC they are a pessimization.

## Barriers? `mf.hfb`, and it is NOT in the package registry (2026-08-29, ledger 162)
`mf.hfb()` / `.line` / `.gpkg` / `.enclose` / `.flopy`. A barrier sits on the FACE
between two cells, so MF6 wants a cell PAIR — and FloPy validates each cellid
against idomain but **never checks the two are connected**, so a bad pair kills the
run mid-way. `.line`/`.gpkg` resolve a fault trace to the CONNECTIONS it cuts --
centre-to-centre segments crossing it (`vor.barrier_faces`, ledger 216) -- and every
route validates first. They took the faces the trace crosses until 2026-09-11: those
run between cells ALONG the fault, and 68-100% of the cross-fault flow still crossed
while `results.q.summary()` read as sealed (it sums |q| on the barrier's own faces).
**Deliberately not a `package_registry` entry**: HFB is face-indexed (no cellid, so
`record_fields` is undefinable) and MF6 writes **no HFB budget record** at all. A
registry entry would have to lie in the one file whose premise is that it cannot. It
rides `run_model._NON_REGISTRY_SUFFIXES` instead, one line, like `mvr`.
**But the FLOW is still queryable** — `model.packages.hfb.results.q.get()`. A barrier
sits ON a connection and `FLOW-JA-FACE` carries every connection's flow, looked up via
the model's own `IA`/`JA` from the `.grb` (NOT reconstructed — idomain makes MF6
renumber). Inputs side: `get`/`summary`/`segments`/`map`, all hand-written, because
the registry view is cell-keyed and a barrier is an EDGE.
Three measured traps: a **duplicated face** is applied TWICE by MF6 and left
permanently wrong (`condsat_reset` restores the already-modified value); a closed
wall through `.line` is watertight now too (82 cuts on the canonical ring against
`.enclose`'s 83), but only `.enclose` refuses a cut that leaks; and `hydchr` is
K/thickness, 1/T, not K.

## Importing a MODFLOW-USG model? `mf.read_usg` (2026-08-28, ledger 159)
`mf.read_usg(nam, gsf=)` -> `UsgModel` -> `.to_mf6()` -> `SimulationSpec` on **DISV**.
Module `src/myflopy/modflow/usg/`; tests `tests/test_usg_import.py` build their own
synthetic 2-layer USG model, so CI never needs a real one.
**Read the report, not just the spec:** `usg.report()` states every approximation and
omission with counts, and `usg.validate()` names what MODFLOW 6 will REJECT before you
run (USG accepts a head boundary below its cell bottom; MF6 refuses to start).
Three traps that cost a day and will cost it again:
- **Coordinates must be LOCAL.** MF6 builds DISV conductances from raw vertex
  coordinates; on State Plane (~1.34e6 ft) it loses the precision and returns a **NaN
  budget while printing "Normal termination"**. Measured: 0/9405 cells finite as-is,
  9405/9405 shifted. `to_mf6(local_origin=True)` is the default -- do not turn it off.
  So a converted model's own frame is NOT the world's: PRT release points and pathline
  maps convert between the two (ledgers 218, 219), and the converted NPF saves the
  specific discharge and saturation PRT reads (217).
- **`complexity="MODERATE"` kills MF6 with SIGFPE** on a converted USG model; the
  default is `COMPLEX`. And SMS's delta-bar-delta/backtracking tuning IS carried over
  (it maps field-for-field) -- it is what made the original converge.
- **An empty period dict is not an absent one.** `steady_state={}` makes FloPy write
  empty period blocks, dropping every TRANSIENT flag; MF6 then solves steady and
  returns NaN. Pass `or None`.
**CLN has no MF6 counterpart and is NOT converted** -- read, segmented into waterbodies
vs streams, and reported (`model.cln_polygons()` gives them as polygons for a later
LAK/SFR rebuild). On the Ten Trails model that costs no pumping (the whole WEL package
is CLN-local P-ET) but removes the lake/stream stage feedback.
`ETS` becomes a **list-based** EVT because MF6 cannot combine segments with
READASARRAYS -- 72 x 9,090 records, 61 MB. That is the physics, not a format choice.

## Re-gridding an imported USG model? `usg.export_gis` (2026-08-30, ledger 169)
`to_mf6()` binds the model to the grid it converted on. `usg.export_gis(dir)` writes
its **content** -- vectors, rasters, tables -- so it can be rebuilt on any later mesh,
which is usually the reason a USG model was imported. Column names ARE the builders'
parameter names (`streams` carries `rwid`/`rhk`/`rgrd`/`rbth`/`man`, `lakes` carries
`strt`/`lake_bottom`/`bedleak`), so `mf.sfr(streams=..., width="rwid", ...)` needs no
renaming. Full file-by-file map in `docs/package_api_reference.md`; runnable
workflow (export -> new model, one package at a time) in
`examples/mf6/notebooks/usg_export_to_new_model.ipynb`.
Four measured facts drive the design, and each is a trap if forgotten:
- **A list BC is a LINE, not a patch** -- mean in-set neighbour counts 1.87 (DRN) /
  1.94 (GHB) / 2.52 (CHD), against 4-6 for anything areal. So they are written as one
  segment per source cell, and every conductance gets a **`*_per_ft` twin**: rebuild as
  `cond_per_ft x new length`. It is NOT a leakance -- `corr(cond, cell area) ~ 0` over
  205 DRN records, and GHB is one constant per family across a 249,000x area spread.
  **A chain is joined by a shared FACE** (ledger 244, 2026-09-12) -- `grid.shared_face`,
  the rule the HFB views resolve a barrier's own face against (and pair-for-pair the older
  `adjacent_cells_idx` on every mesh measured), not `touches`, which a single shared
  VERTEX satisfies. On a structured or quadtree mesh that measured a diagonal
  staircase of six 25 ft cells at 145.71 ft where it occupies 125.0 -- three of its
  segments 30.18 ft long inside 25 ft cells -- and made a 2x2 block a clique whose walk
  stepped across the corner (85.36 ft, `conductance_per_ft` 3.31 against 4.00). It changes
  nothing on a Voronoi mesh, whose vertices are 3-valent, which is why ledger 169 did not
  see it. A cell left ALONE ON ITS CHAIN -- nothing shares a face with it, or it is a
  branch arm the open walk could not reach -- keeps the degenerate segment a lone cell
  always had (`length_ft` 0, and a NaN `*_per_ft` where the package HAS one -- CHD carries
  no extensive column), so the export counts them
  (`cells_alone_on_a_chain`), warns, and points at the POINT layer. `mf.Spread` is
  unaffected -- the feature total comes back exact either way, on the same mesh and on a
  4x finer one; it is the README's per-foot hand recipe that goes NaN.
- **A repeated period is written once** -- 72 periods, 12-month cycle, 11 distinct
  recharge arrays. `periods.csv` is the index.
- **A CLN stream bed IS the old model top** (max |diff| 0.00003 ft on all three creeks),
  and that top was CONDITIONED monotonic along the channels: 0 rises in 97 CLN nodes
  where the raw new DEM has 17% at 300 ft sampling and 38% at 50 ft. Clamping a new
  model's streams to its top reproduces the old model, but needs the same conditioning
  -- and finer grids make it worse, not better.
- **A centerline drawn independently runs mouth-first.** Both Ten Trails creeks do, so
  projecting nodes onto them raw gives a station that DECREASES downstream and a profile
  that reads as climbing. `_oriented_station` flips on the sign of `corr(order, station)`.
- **The CLN carries NO bathymetry.** Every feature's bed is the model top, lakes
  included (Keevie 313/313 nodes, Marjorie 104/104, to 0.00003 ft). So a LAK built
  from the exported bottoms has zero depth, and 312 of Keevie's round just ABOVE
  their cell top. Through the raster, `mf.lak` now caps those at the top and WARNS
  (ledger 231); it used to reject them as "bottom does not intersect active cell",
  so the warning is the signal now, not an error. Supply your own depth. Also: 3 of RockCreek's 39 cells are inactive in
  every layer (USG allowed the connection, MF6's SFR will not) -- see `cell_active`.
**`free_row` (`usg/_io.py`), not `free_floats`, for any row with a trailing label** --
`Wlnd217` contributes a phantom `217.0` to a regex number scan. Same class of bug as the
`IPRN`-as-multiplier one that negated an ET surface.
CLN features are now grouped by **name**, not by connected component: components merge a
tributary into its trunk (6 features found where the file names 7). `FSKIN`, the conduit
table and the node labels are kept -- `FSKIN` was read and discarded before.

## Legacy OO API (the engine; avoid for new model assembly)
`src/myflopy/modflow/mf6/*.py` — `simplemodel`, `boundaries`, `sfr`, `lakes`, `recharge`, and the
builder classes. Still the engine under the facade and still used by `canonical_example.py`. When
adding a capability, **grep both layers first** (package-first + these) to avoid duplicating one.
**`docs/myflopy_context.md` is the accurate, code-derived capability map** (rebuilt 2026-06-20).
Treat any "gap" as a hypothesis to re-verify against the code before building.

## Planned work: PEST / pyemu integration

### What's already built (`src/myflopy/modflow/mf6/pest/`)
- **One unified parameterization API**: `cal.parameterize(target, style=...)` compiling
  to native `pyemu.utils.PstFrom`. Targets: `k`, `k33`, `recharge`, `chd`, `ghb.cond`/
  `ghb.bhead`, `drn.cond`/`drn.elev`, `wel`, `uzf.vks` (`vks`), and transport `porosity` (`mst.porosity`/
  `n`), whose GWT sibling is resolved via `pest/model_lookup.py` — calibrate it with
  heads in the mix, since `v = Ki/n` makes K and porosity near-collinear from
  concentration alone. Styles: `constant`, `zone`, `grid` (one
  geostat-correlated multiplier per Voronoi cell + correlated prior), and `pilotpoints`
  (IDW from a `pp_space` net or explicit `pp_points` — `pilot_points.py`; pyEMU's own
  pilot points are unusable on unstructured grids). The legacy `add_parameter` +
  `build_pst` + `KPilotPointParameter`/`DrainElevation`/`DrainConductance` specs were
  RETIRED (deleted 2026-06-22) — do not reintroduce them.
- `PestProject` (`project.py`) — orchestrates pyEMU/PstFrom; native `build()` + `run_ies`/
  `prior` (`workers=` runs parallel PESTPP-IES agents). **Construct it with
  `model.pest(name, start_datetime=...)`** (the front door — defaults the workspace to
  `<model workspace>.pest/<name>` — a SIBLING of the model dir, ledger 107), not by
  importing `PestProject` directly.
- **Observations** accepted directly by `cal.observe(...)`/`cal.forecast(...)`: `HeadTargets`,
  `LakeStageTargets`, `SfrStageTargets`, `SfrFlowTargets`, `DrnFlowTargets` (or their
  pre-built `*ObservationSpec`). All are wired into the native build + forward run.
- `forward_run.py` — injected forward run: pyEMU's `apply_list_and_array_pars` for params,
  plus post-processors that regenerate head + named-series (lake/SFR/DRN) simulated CSVs.
- **Run discovery + review**: `model.pest_runs` / `run.pest_runs` (`runs.py`,
  `find_pest_runs`/`PestRunHandle`) list calibrations done on a model; `.review()` reopens
  one as `IesResults` (`ies.py`, `open_ies_run`) — phi, ensemble-vs-obs, forecasts,
  parameter-field maps. (The old `results.py` deterministic-review layer was deleted.)
- `geostats.py` — `ExpGeoStruct` / `build_geostruct`: the **single** geostruct builder
  (`project._geostruct_for` delegates to it); `grid` style takes `anisotropy`/`bearing`/
  `nugget` via `parameterize`.
- The PEST notebooks (`canonical_04/05/06`) calibrate the **canonical valley
  model itself** (perturb its K/recharge as "truth", then calibrate back); 06 uses
  pilot-point K. The old standalone demo models (`gold_standard_demo`, `synthetic_demo`,
  `modern_pest_demo.build_calibration_demo`) were retired in favor of one model everywhere.

### What's missing / next steps for PEST
> Verify against the code before building — most of the old list is now done
> (recharge/wel/ghb/chd are `parameterize` targets; IES + parallel workers ship via
> `run_ies(workers=)`; named-series obs are wired).
1. ~~Zones from raster~~ **DONE 2026-07-30** — `mf.ZoneSpec` (raster majority-vote /
   polygon / array), normalized per target family (`pest/zones.py`, ledger 118).
2. ~~UZF parameters~~ **DONE 2026-07-30** — `uzf.vks` (ledger 117).
3. ~~Regularization helpers~~ **RETIRED, do not build** — PESTPP-IES ignores prior
   information equations and a version-2 regularized pst fails to PARSE, so it would
   break `run_ies` (measured). Replaced by the correlated **pilot-point prior**
   (ledger 119); for IES the knob is `run_ies(ies_reg_factor=…)`.
4. **Sensitivity / identifiability analysis** — the one item still open. Splits by
   jacobian: ensemble-based needs no new run; `Schur`/`ErrVar` need a PESTPP-GLM run
   mode (IES produces no jacobian).

### Broader things myflopy could learn from modflow-setup (DOI-USGS)
> Re-verified against the code on 2026-06-20; PEST/serialization/sampling rows refreshed
> 2026-07-03. **Most items previously listed here are already built** (see
> `docs/myflopy_context.md`). Genuine remaining gaps only:
- **YAML spec serialization — DONE 2026-07-23 (§5.6).** `SimulationSpec.to_yaml()`
  /`from_yaml()` (+ `Project.add_simulation_from_yaml`) wrap `to_dict`/`from_dict`
  in `specs_io.py` (PyYAML core dep, safe mode). §5.6A first taught `_callable_ref`
  to serialize the list-BC `functools.partial(_build_named, cls)` builders, so
  chd/ghb/drn/riv/wel/rch/evt now round-trip too. TOML is deferred (no null type;
  `tomllib` is 3.11+ vs the `>=3.10` floor) — ledger 54. Example:
  `examples/mf6/yaml_spec/`.
- **NHDPlus direct SFR reader** — `SFRBuilder` already builds reaches from any stream
  centerline LineString table; only national NHDPlus ingestion is missing
- **Reading existing MODFLOW array files** as source data — a raw array-file *source* is
  missing (existing-model *surfaces* import via `LayerStack.from_modflow`). Note
  **area-weighted raster sampling is done and the default** (`grid/surfaces.py`
  `_area_weighted_sample`) — not a gap.
- **`GridSpec.structured` / `from_geopackage` resolution** — the constructors exist but
  **fail fast** (unwired); myflopy resolves `voronoi` + `python` only. Use
  `GridSpec.from_object` for a built structured/existing grid.
- **LGR parent-child model pairs** — absent (niche for a Voronoi-first toolkit)

Already built — do NOT rebuild: GIS-driven BCs (`GeoPackageSource.chd/ghb/drn/riv/wel/rch/evt`,
`mf.ghb.gpkg`, legacy `Boundaries`), CRS reprojection (vector + raster), layer surfaces +
reconcile + **pinch-out/idomain** (`surfaces.py`), SFR from centerline (`SFRBuilder`),
recharge from GIS/PRISM (`RCHBuilder`).

## Comparison: myflopy vs modflow-setup
- myflopy wins on: unstructured grids, transport/energy/PRT models, visualization, parallel
  runs, Python-first API, GIS-driven BCs, SFR-from-centerline, PEST on Voronoi grids
- modflow-setup wins on: NHDPlus SFR, LGR (myflopy now also loads a full simulation
  from a single YAML file via `SimulationSpec.from_yaml`, §5.6)
