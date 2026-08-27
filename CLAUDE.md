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
height field). Everything returned is a `viz.Picture`: `.fig` / `.show()` /
`.save(path)` / `.html(path)`, and never a trailing `.plot()`.
`tests/test_plot_vocabulary.py` pins the verb set at every scope and fails naming
the stray verb, so adding one in one place and forgetting another is caught.
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
Runnable tour: `examples/mf6/notebooks/plotting_vocabulary_tour.ipynb`.
RETIRED, do not reintroduce: `model.cor`, `model.srf`, `model.visualize`, the
eleven `vor.*` plotting aliases, `stack.preview/views/vtk_3d/surface_3d`.

Only the *spatial* half of this was documented before 2026-07-18; derived tables
had no rule, and the gap produced `sfr.results.long_profile`/`plot_long_profile`
plus ~20 lines of hand-rolled matplotlib in two notebooks that redrew — off-color —
a figure the library already built.

## Compromise ledger (standing user rule, 2026-07-17)
**`docs/compromises_and_deferrals.md`** records every deliberate scope cut,
test-fidelity trade, judgment call, or considered-but-omitted optional inside
delivered work. Any change that makes such a call MUST add/update its ledger
entry **in the same pass** (it is part of the docs-always-in-sync rule).
Entries are removed only when the compromise is actually undone.

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
myflopy's own `model.package_names`). Note `getattr(gwf, "package_name_dict", {})` still
warns: reading the attribute is what warns, so the default only covers absence.

## Catching exceptions? Name them (plan 7.3, 2026-07-31)
**No bare `except:` anywhere** — it swallows `KeyboardInterrupt`. **`except Exception`
needs a `# noqa: BLE001` and a comment saying why the set cannot be closed**; ten
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
- Full suite (**1609 passed / 1 skipped**, 2026-08-24): `pytest -n 10` ≈ **75 s** (worksteal dist is in
  addopts); serial (`-n0`) ≈ **4m30s**; inner loop `pytest -m "not slow"` ≈ 32 s.
  Serial is ~3.7x the parallel run — that gap is the price of the `-n0`
  pre-commit check below, not a regression.
  **Verify sign/column changes with `-n0`** — a session-fixture/xdist interaction
  can report green while serial catches real failures (see ledger 48).
- Tests share ONE session-scoped canonical model on
  `CanonicalModelConfig.testing()` (21×21, smallest contract-complete profile —
  all 15 packages + 5 obs families). Weekly CI re-runs everything on the 50×50
  profile via `SIMPLE_MODFLOW_CANONICAL_PROFILE=validation`.
- Tour/verify interactively: `examples/mf6/notebooks/canonical_fast_tour.ipynb`.

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
  `.plot.map()/.plot.section()/.plot.surface()` and `.vtk_3d()`, `from_modflow`).
  **Pictures come from `stack.plot`** as of 8.5a — `preview`/`thickness_map`/
  `cross_section`/`surface_3d`/`views` are gone. `qc()` stays a method: it is a
  report you read, not a picture. It **compiles to** the
  `LayerSurfaces` engine (`surfaces.py`: area-weighted sampling, top-down reconcile, pinch-out
  → idomain), which uses atomic `Surface` objects (`raster`/`from_contours`/`from_points`/
  `from_array`/constant + algebra). Use `LayerStack`; `LayerSurfaces`/`Surface` are the engine/atoms.

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
