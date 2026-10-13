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
budget basing, suffix map, all three artifact sets, budget term). **8 still need
hand-writing**: the `<pkg>_spec` factory, the `package_api` helper class and the
`GeoPackageSource` resolver — all deliberate, they carry prose and explicit
signatures — plus five namespace properties (`ModelPackages`, `SimulationBase`
for `model.<pkg>`, `FlopyPackages` for `model.fl.<pkg>`, `GroupPackages`,
`_PackageDiffNamespace`). Those five stay hand-written **on
purpose** (4.7.7): a `.pyi` stub replaces a whole module for type checkers, so
generating them would mean hand-maintaining stubs for every other public name in
four large modules. Instead `test_package_descriptor.py` asserts all five
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
`*_spec` output, and per-surface package list -- and the parameters (names, kinds,
defaults; NOT annotations) of EVERY public callable (ledger 443): the set
`scripts/public_surface.py` discovers, the same one `test_public_parameters_documented`
measures, so a default changing on `Run.reopen` or `with_grid` shows in the diff and
a type-only change does not. If a change is intentional, rerun
`python scripts/derive_api_snapshot.py` and **review the diff** — it is the public-API
change log. An un-regenerated snapshot fails CI.

## Adding a table, plot, or map to a package? Read `docs/view_layer_conventions.md`
**A binary output file's rows are labelled from ITS OWN clock** (ledger 301,
2026-09-23) -- never from `model.kstpkper`, which is the HEAD file's. One MF6
`OC` `saverecord` sets two INDEPENDENT flags, `idvsave` from `SAVE HEAD` and
`icbcfl` from `SAVE BUDGET`, and LAK/SFR gate their STAGE file on the first
while gating their BUDGET file on the second -- so within one package the two
outputs follow different schedules. `SAVE HEAD LAST` + `SAVE BUDGET ALL` at
`nstp=5` writes **2 head records against 10 LAK budget records**;
`STEPS 1` vs `STEPS 3` writes **two of each naming different time steps**, which
no length check can see. Both were live: the first crashed
`packages.lak.budget.get()` and the second stamped every flow with the wrong
step, silently. `SFRStageOutput.get` is the reference pattern. **The budget
FILENAME comes from the OC record too** -- `<name>.cbc` is a default, not a
rule, and assuming it left every cell-budget surface raising `FileNotFoundError`
on a real model whose OC said `BUDGET FILEOUT tentrails.cbb`. **A reader whose
records do not line up must RAISE** -- the removed fallback warned once (Python
dedupes by line, and IPython never resets `__warningregistry__`) and returned a
raw list from a method documenting a DataFrame, and the `isinstance` guard that
caught it downstream `continue`d, returning **0 rows against a file holding 38**.
**A model with TWO packages of one type** (two `RIV6` entries, a main stem plus
a named tributary) writes two budget records per saved time, which overran the
model-budget clock the same way -- `model.bud("riv").df` raised `IndexError`.
**The STAGE clock too** (ledger 348): `build_lak_stage_result_table` labels each
record from the stage file's own headers (`model.outputs.lak.stage.kstpkper()`) and
`per` is the stress period -- it was the saved record's INDEX, so at `nstp=3` a lake
target keyed by period read steps 1-3 of period 0 and PEST got nine rows for three
periods -- and a record count that disagrees with the headers RAISES. `group/sfr.py`
decoded a reach-by-time frame as time-by-reach (`per` held the reach); both group
stage tables are rebuilt on the single-model tables. STILL OPEN: `_kstpkper` is never
invalidated.
The rule is `model.<pkg>.<inputs|results>.<noun>.<verb>()` (ledger 339: `model.drn` is
myflopy's view, **FloPy's package object is `model.fl.drn`**, and `model.packages` keeps
only the all-packages `summary()`/`mosaic()`; group and diff still read
`group.packages.<pkg>` / `diff.packages.<pkg>`) — **every noun
is an object, every object answers the same verbs** (`get`/`summary`/`plot`, plus
`map`/`mosaic`/`animate` for spatial nouns, and `section` -- a line profile -- only on
the dependent-variable nouns, `LineProfileView`: a verb that can never act is not
offered, ledger 462). Never add a `foo()` + `plot_foo()`
method pair; add a view class exposed as a noun property.
**Several packages of one type are found by TYPE and never picked in silence** (ledger
352; Ten Trails carries RIV `riv` + `crisp`). FloPy's `get_package("lak")` answers a
package NAMED `lak` first, so `model.fl.lak`/`model.package("lak")` returned `lak` alone
beside `lak_2` and every LAK view read its lake as the model's. ONE rule now:
`package_registry.model_packages_of_type` (RCHA/EVTA count as RCH/EVT) -- `model.fl.<pkg>`
is a LIST of several, `model.package(s)` a type then an exact name (no filename or
prefix fallback). Views whose rows are CELLS read every package, each row's `pname`
naming it -- the budget-file reader kept each time's FIRST record, so `crisp` (53% of the
river exchange) was absent from `model.riv.results.q` -- and match a term EXACTLY by
`paknam2` (the `paknam=` it tried never matched: MF6 writes the model's name there).
Views of one package's FEATURES (LAK lakes, SFR reaches, UZF `ifno`) go through
`one_package`: several raise `mf.SeveralPackagesError` naming them and the ONE spelling
`model.lak(package="lak_2")` (`model.sfr`/`model.uzf`/`model.outputs.<pkg>` alike), and
each explorer passes `package` to every builder below it. It is a `LookupError` ON
PURPOSE: the group/diff tiers catch `ValueError` as "package missing". Trap: a property
raising `AttributeError` makes Python retry `__getattr__` -- `FlopyPackages` re-raises
the type's refusal there, or `model.fl.lak` answered a GHB NAMED `lak`.
**Several records on ONE cell combine by what the number IS** (ledger 297,
2026-09-22): `FieldSpec.kind` -- `extensive` (a flux or conductance; MF6 superposes
them, so they ADD), `elevation`, `intensive` -- and `agg` is DERIVED from it, never
declared beside it. Every entry read `"sum"` before, elevations included, so two CHD
records of 100.0 and 102.0 rendered **202.0**. Two namespaces exist and three of seven
packages differ between them: `_ELEVATION_FIELDS`/`_EXTENSIVE_FIELDS` are keyed by the
HELPER's argument names (`drn.elevation`, `ghb.head`, `wel.rate`), NOT the record names
(`drn.elev`, `ghb.bhead`, `wel.q`) -- a set rebuilt against `FieldSpec` names matches
NOTHING for those three and raises nothing. `_EXTENSIVE_FIELDS` is deliberately a
STRICT SUBSET of `kind == "extensive"`: it asks whether the declared column holds a
TOTAL, which an unwrapped conductance does not.
**The rule is ONE function, `package_registry.agg_for_kind`, and the observation targets
ask it too** (ledger 331): `_aggregate_named_series` takes a REQUIRED `kind=`. A stage
table repeats a lake's one stage on every connected cell, and the targets ADDED the
copies -- the user's pond, 763 cells, read **401,339.97 ft for 526.00**, in the CSVs PEST
read. Copies that disagree are REFUSED, not averaged. Every INPUT view -- map, series,
group, diff -- defaults through `package_registry.input_field_agg` (ledger 472; the
series and group tier had summed elevations), and records pair across models one to one.
**A period's value is its END-of-period value, for a flow too** (ledger 336, the user's
rule): every family -- heads, conc, SFR stage, SFR flow, DRN flow, PEST's forward run,
`model.hds.get()` -- picks the record through ONE function, `period_end.period_end_steps`,
on the file's OWN clock and by `pertim == perlen` (ATS breaks `kstp == nstp`) to 1e-8,
because FloPy WRITES `perlen` at 8 decimals (at 1e-9 a live model's hourly periods were
all refused), and a period saved without its end is REFUSED naming `SAVE ... LAST`. The
SFR/DRN flow targets had ADDED a period's steps: the canonical outlet at `nstp=3` read
**113,453.29 for 37,817.76**. Do not add another "last saved step".
**On a picture or a group table, `per=N` is the LAST SAVED step of period N** (ledger
337, the user's words), through ONE function, `period_end.last_saved_step(s)`: the map,
both sections, `surface`, the section composers, `HeadsPlus.array`, the group's
`_coerce_kstpkper`/`_reduce_to_period_end` and `canonical_sfr_signals` call it, and
`period_end_steps` is it plus the end check. Deliberately NOT the end-checked rule, so a
`SAVE HEAD FIRST` run still draws. The profile, the noun, `surface` and a group read `per`
as a saved-RECORD index until then: `per=1` at `nstp=3` drew (1, 0) beside the map's
(2, 1). **A view that reduces a results table to one value per period reads that step
through `period_end.at_last_saved_step(s)`** (ledger 348) -- before ANY rows combine:
the results maps chose the step in the payload builders, every series in
`SpatialView._series_selection`, and `stage_change`, the SFR long profile, the LAK
`q.budget`, the surface-water table, the group stage tables and diff maps and
`region_heads(per=)` likewise -- and so do the TABLES keyed by period (ledger 349):
`q.wide()`, `model.lak.budget`/`model.sfr.budget` `summary()`/`wide()` (grouping by
`kstpkper` keeps every step) and the SFR field profiles. They had ADDED a period's
saved steps (a LAK exchange map at `nstp=3` read -2793 ft3/d for -587; once the map
moved, a CHD `wide(per=1)` read -2152.678 beside its map's -452.712) or painted its
FIRST (`agg="first"`). `get(per=N)` stays a FILTER: every step of period N; `long()`
keeps every step.
**The OBSERVATION route is one row per period at its end too** (ledger 346): an MF6 OBS
CSV is every time step labelled by `totim` alone, so `period_end.period_end_rows` finds
each period's END row -- measuring each period from MF6's OWN previous end, since MF6
adds FloPy's 8-decimal lengths and a live hourly model outgrows 1e-8 of the cumulative
`perlen` from the fourth period -- through `period_end_steps`' check; a file short of an
end is refused, a target keyed by MF6's time (`10, 20, 30`) is refused naming the period
WHEN A KEY NAMES NO PERIOD, and ATS no longer moves the row count PEST's positional `.ins`
reads (12/13/14 rows at K = 10/1/2 read t = 23 as the last period's end). **A time key
below `nper` IS a period and reads that period in silence** (ledger 347) -- not a daily
corner: on twelve 2-day periods the old keys 2, 4, 6 read periods 2, 4, 6; only a set
reaching past the last period is refused, so audit a set written before 346. A LAK
STAGE-FILE target is keyed by period too since ledger 348 (it was keyed
by SAVED RECORD: at `nstp=3` periods 0-2 read steps 1-3 of period 0).
**A flow target is MF6's own sign on BOTH routes** (ledger 350, the user's call): an SFR
reach reads MF6's `DOWNSTREAM-FLOW` and a DRN zone its `DRN` flow, zero or NEGATIVE
(water leaving), so a gauged `Q` is entered `-Q`, and a target above zero is REFUSED at
construction naming `-Q` (`_shared._refuse_flows_leaving_the_frame`). The SFR file route
took the MAGNITUDE, so one reach read with opposite signs under two names and a positive
`Q` matched exactly scored `-2Q` (phi 327,106 on three exact observations), a misfit only
a drier reach lowers. Do not normalize a flow's sign anywhere on this path.
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
`.save(path)` / `.html(path)`, and never a trailing `.plot()` -- bar a `backend="mpl"`
figure, which is `section`'s DEFAULT (next paragraph).
**`section()` IS the model's cross-section, drawn by Matplotlib** (ledger 333, the
user's call, 2026-09-26; the picture since ledger 177): the defaults are
`fill="layer"` -- cells coloured by layer, the head as a water surface on top -- and
`backend="mpl"`, the ONE verb whose default is static, so a bare call returns a bare
`Figure` (`.savefig`), not a Picture. `fill="results"` colours the cells by
head/conc/temp with a colorbar; a `(nlay, ncpl)` array paints anything else.
**`fill=None` is the LINE PROFILE** (head vs distance, no cells -- the old default)
and the ONE spelling of it (both pictures open on the FIRST saved time, ledger 337;
its interpolation-only arguments -- `num_points`, `use_rbf`, `clip`, `surf_type`, the
surface tail -- RAISE without interpolation, naming `interpolate=True`; `spacing` and
`extrapolate_beyond_section_ends` never acted and are GONE; it reads its model top from
the model it draws, where a FloPy reload printed and took the FIRST flow model's): a profile argument (`layer=`, `interpolate=`, ...) given
without it RAISES naming `fill=None`, so the picture never depends on which other
arguments were passed. Both pictures draw on BOTH backends, the same cells (ledger
332): one `FilledSection` reads one FloPy `PlotCrossSection`, so Plotly is not a
second derivation of the geometry. `.ani` stays the Plotly profile's
(`fill=None, backend="plotly"`). It raises rather than degrading: on the other
picture's arguments, on `"results"` or a results argument for a bare grid, and on
a flat `(ncpl,)` array. **The NOUN `model.hds.section()` stays the head LINE
PROFILE** (ledger 334, the user's call): it names what measured to act
(`NOUN_SECTION_PARAMS`), and `fill=` and the cross-section's other arguments RAISE
there, naming `model.plot.section(fill=...)`; its `**refused` tail forwards nothing.
`XSection` checks its own tail when built, so a typo raises on every section.
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
**A `surface` is MASKED to `clip=` or the grid's footprint** (ledger 314). The
interpolator is not a mask: `griddata` is NaN outside the convex hull of the cell
centres, which looks like a clip and spans every concavity of a non-convex domain --
`.show()` drew exactly that while only `plot_heatmap` masked. Draw
`InterpolatedSurface.clipped_surface`, never `.surface`. One mask rule (pixel centre
inside, `geometry_mask` on `transform`), and `transform` is registered to
`xy_meshgrid`'s pixel CENTRES -- it sat half a pixel off, so a saved raster and
`XSection(interpolate=True)` moved when that was fixed. The rim past the outermost
centres stays empty (ledger 248).
**A no-value cell is `|value| >= 1e29`, both signs: `headsplus._is_sentinel`**
(ledger 315). MF6 writes +1e30 inactive and -1e30 DRY. The head surface cut at
`> 10_000`, so a dry cell cratered it to -1e30 and real heads above 10,000 ft were
erased. The map and the section tested their own sign or none, until every read in
`choros.py` went through `Choro._layer_field`. Before that, a dry cell set every
colour range to -1e30 and made `contours=True` draw NOTHING. **The user's calls:**
on a mounding map a DRY cell reads 0 and an inactive one is blank; on a
depth-below-ground (`bgs`) map both are blank.
**A lake with no stage is the same two numbers** (ledger 329): LAK writes -1e30 DRY and
+1e30 INACTIVE, and every stage view drew them. Masked ONCE where the file is read,
`model.outputs.lak.stage.get()` (`masked=False` for the raw file), so the table, maps,
series, exchange hover, group and diff inherit it. A table keeps the ROW with a `NaN`
(dropping on it deleted the transitions INTO a dry period, and `LakeStageTargets`' rows);
nothing may SUM an all-NaN group, which pandas makes 0.0 -- `_aggregate_leaving_gaps`
serves the series and both map payloads (`agg="sum"` painted a dry lake 0 ft); a diff
COUNTS a one-sided dry lake (`n_one_sided`). **A lake TARGET is not a view** (ledger 351,
the user's rule): a DRY lake is calibrated against its BOTTOM on both routes -- MF6's own
STAGE observation of it, since the solver clamps the stage there -- so the stage-file
route reads the file RAW (`build_lak_stage_result_table(masked=False)`) and puts
`package_tables.lake_bottoms` in place of -1e30, MF6's `lakebot` rule (VERTICAL: the cell
TOP, never the BELEV written; HORIZONTAL: BELEV, or the cell BOTTOM when BELEV == TELEV;
the lowest; EMBEDDED refused -- its bottom is its table's). 1e30 there stalled PESTPP-IES
(lambda stuck near 1e57, phi 1.11 against 0.085; `docs/reproducers/dry_lake_bottom.py`).
Exact on a model LOADED from its files; a model still in memory is off by FloPy's
8-decimal write (<= 5e-9). What PEST reads for an INACTIVE lake is 1e30, every row kept:
PEST++ fails a run on an empty field, on `nan` and on a missing row, and reads an extra
column as the wrong one (`docs/reproducers/pest_no_value_fields.py`). **A reach with no
stage is a gap too** (ledger 353, the user's call): `model.outputs.sfr.stage.get()` masks
SFR's -1e30 (dry) and +1e30 (inactive) the same way, so every SFR stage VIEW -- maps,
series, both profiles, the q hover, `canonical_sfr_signals`, group, diff -- leaves it
blank. **SFR TARGETS were left as they were**: `SfrStageTargets`' file route reads the file
RAW (`build_sfr_stage_result_table(masked=False)`), so a dry reach still reaches PEST as
-1e30, while MF6's STAGE observation reports the streambed top there (dry AND, measured,
inactive) -- the two routes differ at such a reach. **A target's PICTURES are stage views**
(ledger 355): `CalibrationPlot`'s `timeseries`/`obs_vs_sim`/`residuals_by_period` blank a
simulated `|value| >= 1e29` (`calibration._blank_no_value`) while `compare()` keeps it --
its 1:1 line ran from -1e30 to +1e30.
**Highlighting cells is `map(select=...)`, on EVERY scope** (ledger 160,
2026-08-29) -- cell indices, a boolean mask, a registered region name, a vector
file, or a geometry. It draws a **dissolved-boundary outline** by default, NOT
Plotly's `selectedpoints`: selection styling on a choropleth exposes only
`opacity`, so it can only dim the rest, which measured a **6.9x loss of readable
contrast** across the unselected cells on a real head field. `select_style="dim"`
restores the old picture, `"both"` draws each. Applied on `Choro` (not in the
verb) so `mosaic`/`animate` keep it, and `mode="lines"` so a user's lasso cannot
clobber it.
**Named areas of interest: `model.aoi` / `project.aoi`** (ledger 310) -- `mf.Area`s by
name, usable as `extent="site"`, `select="site"`, `features=model.aoi`. DEFINED ON THE
MODEL, KEPT BY THE PROJECT: a model finds its live project through
`myflopy.areas`' root-keyed WEAK registry (a live `Project` registers there, every live
one per root kept, newest wins -- ledger 312), asked for the NEAREST ancestor of the
model's workspace FIRST (a copied project's runs record the ORIGINAL root) and then the
recorded `view._project_root` -- NOT through the five Run construction paths. No live
project -> the SAVED areas, additions session-only (the user's call). Saved INLINE in
`project_spec.json`. Trap that cost a round: a pytest fixture holding the `Project`
keeps it in the weak registry, so a "no live project" test passes for the wrong
reason -- build and drop it inside the test and assert `live_areas(root) is None`.
**`extent=` is the window a map opens on, at EVERY scope and backend** (ledger 309):
VIEW ONLY (the user's rule -- no crop, no colour re-range; `zmin`/`zmax` tighten a
close-up). One resolver, `utils/datatypes/extent.resolve_extent`, run at the call;
plotly frames it through `latlon_bounds`/`map_view`, Matplotlib sets padded limits
LAST in `plot_mpl` and only when given. Four numbers are SHAPELY order
`(xmin, ymin, xmax, ymax)` -- FloPy's `(xmin, xmax, ymin, ymax)` is undetectable when
it happens to be a valid box. A mosaic takes one extent for every panel or a LIST,
one per panel, told apart by the ITEMS (`viz.is_extent_list`); different windows are
never synced (`sync_views=None`; `True` raises). NOT on `grid`/`section`/`surface`.
**Transmissivity: `map(type="trans")` is SATURATED, `"trans_max"` is K x full
thickness** (ledger 308, the user's call), and `packages.npf.transmissivity`/`_max` are
the nouns -- they draw THROUGH the free map (`StaticArrayFieldExplorer._map_type`), so
one picture. Saturated follows MF6's ICELLTYPE rule (confined = full thickness;
THICKSTRT honoured) and matches `DATA-SAT` to 5e-7. One computation:
`modflow/mf6/hydraulics.model_transmissivity` (import layer 1). Total-over-layers is
DEFERRED (needs a non-int `layer=` at every scope). T is WRITTEN `112,000 ft²/d`, not
`1.12e+05` (user, ledger 313): `HoverSpec.formats={"T": "thousands"}` for the hover,
d3 `","` / comma decade labels on plotly bars (composers read `Choro.number_style`),
`hover.thousands_formatter()` on static bars -- ONE rule, `format_number(style=)`.
**A plotly map is fitted for the canvas plotly will draw** (ledger 307, amended by 311):
`viz.map_canvas` (an export size, else the figure's, else plotly's 700x450, less the
margins plotly DRAWS -- it grows them for a colorbar and animation controls, so pass
`colorbar=`/`controls=` -- times the subplot domain) into `Choro.map_view(size=)`.
Measure framing in a BROWSER (`docs/reproducers/map_canvas_margins.py`); a test that
takes `map_canvas` as the truth certified the colorbar bug. The old fit assumed 1000x700 and cropped single maps, mosaic
panels (~10% of the grid) and PNGs; a hillshade reset the view after it. Pass a SIZE
whenever you frame a map, and never set `center`/`zoom` after the fit.
**A map of an UNRUN model reads its inputs, and `type`/`layer` are checked at the call**
(ledger 305, 2026-09-23). `_grid_of` binds a model only when it has RESULTS, so an
unrun model was drawn as a bare grid: `type='ks'`, a typo'd `type`, and the documented
`'rch'` all drew NODE NUMBERS 0..ncpl-1, and `select="<region>"` was opened as a file.
`Choro.model` still means "whose results"; `Choro._inputs` means "whose inputs and
names" and falls back to the input-only model the free verb passes. Do not fake the
results binding for an unrun model -- the heads hover, timing reader and animation
clock all assume results. A results `type` on an unrun model RAISES; a new map type
goes in `choros.MAP_TYPES` with whether it needs results. `'input_rch'` is retired.
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
`mosaic`/`animate` are defined ONCE on `SpatialView`, and they now name all 24
per-panel drawing options plus `backend`, with descriptions spliced from
`plot.map` — so `packages.npf.k.mosaic(logscale=True, features=..., show_layer_elevs=True)`
completes in an editor. Only options the caller actually SET are forwarded
(`_explicit_options`), and they are REFUSED on `kind="plot"`/`"section"`. They
rotted because `test_noun_signatures.py` discovers subjects by a body mentioning
`plot.map(`/`_choropleth_factory`, which a composer's never does; it sweeps them
explicitly now.
**...and HONOUR them** (ledger 306): a panel's `zmin`/`zmax` beat the pooled range
(widest across panels) on both backends -- they were forwarded and pooled over. The
static composers draw each panel through `Choro.plot_mpl(ax=, norm=)`, never by
re-plotting cells, which is what dropped `features=`/`select=` from every static
mosaic and animation.
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
The sweep keys on a body mentioning `plot.map(`/`plot.section(`, so the noun's
`section` (it renders through `render_xsections`) sat outside it with a `**kwargs`
tail that swallowed `fill=` until ledger 334; it is swept explicitly now, and every
`Class.section` with a tail is pinned.
A noun names Tier 1 (`NOUN_MAP_PARAMS`, 25 including `backend`) plus its own
selectors, REFUSES `values`/`type`/the legacy hover trio, and omits the Tier 2
per-layer names that were measured byte-identical on a record noun -- and REFUSES
them when they ride its tail (`refuse_noun_parameters(honours=)`, ledger 465;
`hover_layers` passes where the picture keeps a per-layer hover). **Docstrings are STATIC, and generated** (ledger 178/180, 2026-09-05): the free
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
is the reasoning record. `project.aoi` (areas of interest, ledger 310) is a fourth,
coercing library. `project.packages` / `.grids` / `.simulations` are `dict`
**subclasses** — every dict operation works unchanged, plus `names()`, `summary()`,
a did-you-mean `KeyError`, and `_repr_html_`. `project.package(key)` / `.grid(key)`
mirror the long-standing `simulation(name)`; `package_spec(key)` had existed on the
retired `ProjectSpec` and was dropped in the 2026-06-17 fold, unnoticed because
**`api_snapshot.json` covered no part of the `Project` surface** — it does now
(`Project.members`, `Run.members`, one row per library), so a verb vanishing fails CI,
and every method's parameters and defaults (ledger 443), so a keyword changing does too.
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
on, not that the run can no longer be rebuilt as it ran.
**A second column, `source`, says whether the FILE under a package still says what
it said** (ledger 276) — `ok`/`edited`/`missing`/`unfingerprinted`/`—`. Separate from
`state` on purpose: a package is routinely `unchanged` (the key means the same) AND
`edited` (its GeoPackage moved on), and one column would hide one of them. Nothing
recorded a source's CONTENT before — only its path — so an edited input read
`unchanged` forever. It hashes bytes, so a rewritten `.gpkg` reads `edited` even with
no feature changed: the safe direction.
**It answers about what the RUN READ, not about the library's current declaration**
(ledger 294, 2026-09-21). A run records `metadata['built_from']` — the content hash of
every source it actually OPENED — and `source` compares against that. Without it two
runs prepared either side of an edit were INDISTINGUISHABLE: they share one library
declaration, so they recorded the identical digest while building CHD heads of 9.0 and
42.0, and `drift()` gave both the same verdict — including the run that AGREES with the
file. It also covers packages declared **INLINE**, which had no row at all:
`_package_digests` records only `mf.ref` entries, on the premise that an inline package
"cannot drift" — which lazy packages falsified, since one holds a PATH and opens it at
build. Such a row reads `state` `—`: there is no key, so that question has nothing to
answer about. **Only a LAZY package is re-hashed at build** (only it reads a file
there); an eager one carries its declaration hash, so `source` for an eager package is
PROVENANCE, not buildability. That split is also what keeps the cost honest — hashing is
**0.15%** of the work it describes for an ordinary GeoPackage, but a full read of a
3.7 GB DEM is **~1.9 s**, and the build only samples such a raster. A package naming
several files reports the WORST with a count (`edited (1 of 3)`).
**NOT built, and why** (same ledger): `unfrozen` is a per-RUN fact and `drift()` is per
(run, package key), so it belongs in `runs.summary()`. `legacy-digest` has no era marker to
read, and 0 of 84 real stored digests move under HEAD. **`runs.summary()` has a `grid`
column** (ledger 330): a run records each model's grid pin and how its cells were made
(`run.grids`), and the column reads `unrecorded` for a run from before that, `FloPy's
cells`, or `N cells repaired`. It does NOT re-resolve the recipe to ask whether the grid
changed -- that is a triangulation per grid per call.
`project.save()` returns a **`SaveReport`**,
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

## Writing a public function or class? Document EVERY parameter (ledger 338, 2026-09-27)
The user's rule, asked for repeatedly: one numpydoc entry per signature parameter, in
order, head `name : type` (what the code ACCEPTS -- run a check when the code is
ambiguous, never copy the annotation), a body giving the default, units, which
mode/branch it belongs to and how it interacts with the others, plus a `**options` entry
saying where a tail goes. `tests/test_public_parameters_documented.py` discovers
`__all__`, `__engine__`, the helper methods and every public method on an exported class
(715 callables), plus -- through `scripts/public_surface.py`'s `reached_classes`, a static
walk with no list -- the public methods of every class a user REACHES from those by
attribute, call or returned value: the read-grammar views, `model.plot`, `model.targets`
(400 more, ledger 436); and fails naming the callable -- no allowlist. `python
tests/test_public_parameters_documented.py <file>` lists what a file still lacks. It
checks that a type is WRITTEN, not that it is TRUE: the 2,004 entries written for it were
read against the code by a second pass, which corrected 175. A dataclass users RECEIVE may
use `Attributes`. A gap inside a generated fence is fixed in the free verb, then regenerated.

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

**Importing a bare MF6 workspace? Carry the MODEL's own options** (2026-09-18,
ledger 286). `mf.load_mf6_run(ws)` + a package-by-package rebuild drops them, because
**`NEWTON` is in the model's NAME FILE, not in any package** — measured on a real
18,107-cell model: 90,408 cells converted to DRY at ITER=0, period-1 recharge fell
2,854,736 → 1,645 ft3/d, every head-dependent boundary went to zero, 200% budget
discrepancy, and the solver diverged at 143,000 ft head changes. Spread
`**run.model_options` into `mf.gwf(...)`. `model.config` was blind in the same place —
it walked `gwf.packagelist`, so `ref.diff(new)` compared two runs differing ONLY in
`NEWTON` as EQUAL; `config.model` is the new section. Traps: FloPy returns
`newtonoptions` in three nested shapes, **`newtonoptions=""`/`[]` write NOTHING**
(so normalizing the bare case to empty drops the option while every round-trip test
passes), and the reader KEEPS `*filerecord` entries that the diff extractor skips.
**The SOLVER is still not carried** — the same model drifted `outer_maximum` 2000 → 300;
`config.ims` can read it today, it wants kwargs shaping (plan D2).

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
- Full suite (**8061 passed / 4 skipped**, 2026-10-12, the second and third rounds of the Ten
  Trails PEST preparation -- ledgers 543-563: date-keyed hydrographs, array recharge/storage/
  monthly multipliers, zones that read every layer, the by-value polygon vote, `mf.Basin`, screened
  head targets on one period-end reader, lake/stream feature parameters, the calibrated write-back):
  serial (`-n0`) **20m05s** and `-n 10` **4m11s**, NOT quiet clocks (load averages ~14.7 and ~7.2 at
  their starts; the ratio 4.8x is high for that reason). Before it, 7777 passed / 4 skipped
  (2026-10-12, the first round -- ledgers 499-542: target hygiene, the IES review and run policy,
  parameters by package TYPE with additive elevations, pilot points that stack and interpolate in
  log space, `mf.AverageYear`, RIV flow targets, `mf.RasterValues`): serial (`-n0`) 18m18s and
  `-n 10` 4m22s, NOT quiet clocks -- other lanes' tests and PESTPP-IES runs shared the machine (load
  averages ~14 and ~4.3 at their starts, ~20 during the parallel run). The first `-n 10` run
  failed one test, `test_sfr_stream_inputs.py`'s self-diversion case (MF6 6.8 segfaults on it
  intermittently; it retries 5 times and passed 5/5 alone three times, in the serial run and on
  the rerun). Before it, 7435 passed / 4 skipped (2026-10-12, the backlog's last batch -- ledgers
  489-498): serial (`-n0`) 11m49s and `-n 10` 2m31s, both in a checkout outside `.claude` with a
  RELATIVE `PYTHONPATH=src`, at load averages of ~4.7 and ~3.2.
  Before that, 7374 passed / 4 skipped (2026-10-11, the third round --
  ledgers 447-488: the read views' defects the docstring pass found -- a period read no
  longer mutates a package, `per=N` reads the block MF6 has in force, SFR q per reach,
  diffs by package TYPE and the mover records -- plus PEST on DISU and one unit table per
  dimension): `-n 10` **2m26s** at a load average of ~4.5 (2m35s again at ~2.9) and
  serial (`-n0`) **16m25s**, NOT a quiet clock -- it started at ~13.4 as the parallel run
  ended, with a game client running on the machine; the slowest tests are still the
  PEST++ end-to-end runs (25-26 s at `-n 10`), and no new single test dominates. Before
  it, 6696 passed / 4 skipped (ledgers 412-446, the second round, the docstring guard over
  the read-grammar views and the executed examples of `tests/test_examples_run.py`):
  serial 10m27s and `-n 10` 2m16s, at ~13.5 and ~1.0 (one earlier `-n 10` of that tree
  took 3m57s, unexplained). Before it, 5972 passed /
  4 skipped (ledgers 362-411, the six PEST leftovers and the docstring pass's 68 code
  defects): serial 9m09s and `-n 10` 2m11s, at ~11.9 and ~1.7. Before that, 5583 passed / 4 skipped
  (2026-10-02, the PEST-build series, ledgers 356-361 -- the build reads the run only,
  observation groups, a build deletes its workspace only when a build wrote it, and a
  parallel run its master and workers only when they are copies of its template): serial
  (`-n0`) 7m41s and `-n 10` 1m34s, the parallel one starting at a load average of ~0.4
  and the serial one at ~10 as the parallel run ended.
  Before it, 5569 passed / 4
  skipped (ledger 360): serial 7m41s and `-n 10` 1m40s, at ~8.9 and ~2.5. Before that,
  5516 passed / 4 skipped (ledgers 353-355 -- SFR stage views
  masked, the suite's `PYTHONPATH` made absolute, and their review): serial (`-n0`)
  **7m47s** and `-n 10` **1m36s**, both in a worktree outside `.claude` with a RELATIVE
  `PYTHONPATH=src` (ledger 354), the parallel one starting at a load average of ~1.1. Before
  it, 5514 passed / 4 skipped (ledgers 353-354): `-n 10` 1m42s at ~3.4, serial not
  re-measured. Before it, 5489 passed / 4 skipped (2026-09-30, ledger 352 --
  several packages of one type): serial (`-n0`) **7m37s** and `-n 10` **1m33s**, both in a
  detached checkout outside `.claude` with an absolute `PYTHONPATH`, the parallel one starting
  at a load average of ~3.4 and the serial one at ~7 as the parallel run ended. Before it,
  5474 passed / 4 skipped
  (2026-09-29, the PEST-targets series, ledgers 340-351 on top of 339): serial (`-n0`) 7m41s
  and `-n 10` 1m38s, the serial one starting at a load average of ~8. Before it, 5460 passed / 4 skipped (ledger 350 -- flow targets in
  MF6's sign): `-n 10` **1m38s** in a worktree outside `.claude`, at a load average of
  ~0.1 at its start (1m43s at ~2.7 after the last docs edits), absolute `PYTHONPATH` (a relative `PYTHONPATH=src` failed the two
  tests that run a PEST `forward_run.py` in its template then: the subprocess imported another
  checkout's myflopy -- fixed by ledger 354); serial not re-measured. Before it, 5454 passed / 4 skipped
  (ledger 349 -- the review of 348): `-n 10` 1m38s at ~0.5. Before it, 5451 passed / 4 skipped
  (ledger 348 -- the LAK stage clock): `-n 10` 1m34s at ~0.5. Before it, 5437 passed / 4 skipped
  (ledger 347 -- the review of 346): `-n 10` 1m38s at ~0.2. Before it, 5434 passed / 4 skipped
  (ledger 346 -- the observation route keyed by period): `-n 10` 1m40s at ~1.5. Before it,
  5424 passed / 4 skipped (ledger 345, the review of 344):
  `-n 10` 1m39s in a worktree outside `.claude`, at a load average of ~0.3 at its
  start, with an ABSOLUTE `PYTHONPATH` (the bullet below); serial not re-measured. Before it,
  5421 passed / 4 skipped (ledger 344 on the PEST-target series 340-343): `-n 10` 1m42s
  at ~0.8. Before it, 5371 passed / 4 skipped (2026-09-28, ledger 316's five build-failure
  fixes and their review, rebased onto ledger 339): `-n 10` **1m34s** in a detached checkout
  outside `.claude`, at a load average of ~0.6 at its start; serial not re-measured. The
  same fixes before the rebase counted the same 5371 in 1m30s at ~6.1, and ledger 339 alone
  5325 in 1m29s. Before them, 5333 passed / 4
  skipped (ledger 316's written-resolution amendment on top of 338): `-n 10` **1m32s** at
  ~3.6. Before that, 5325 passed / 4
  skipped (ledger 338 -- the docstring guard adds 717 cases): serial (`-n0`) **7m08s** and `-n 10` **1m38s**, both in a detached checkout
  outside `.claude`, the serial one starting at a load average of ~11 as the parallel run
  ended. Before it, 4608 passed / 4 skipped (ledger 337): serial (`-n0`) **6m57s**
  in a detached checkout outside `.claude` at a load average of ~0.9. Before it, 4575 passed /
  4 skipped (2026-09-26, ledger 316's revision on top of 336):
  `-n 10` **2m08s** in a detached checkout outside `.claude` at a load average of ~10.7;
  serial not re-measured. Before it, 4566 passed / 4 skipped (ledger 336 on 332-335):
  serial (`-n0`) **7m32s** and `-n 10` **2m46s**, both in a detached checkout outside
  `.claude` (load average ~2.5 at the end). Before it, 4542 passed / 4 skipped (ledgers 332-335): serial (`-n0`)
  **7m00s** in a detached checkout outside `.claude` (so `test_docs_structure` runs its
  cases), load average ~1.3 at its start -- the quietest serial clock since 4037. Ledger
  335's tree, the same count, took **14m49s** at a load of ~23.5 from other sessions'
  suites: the load, not its temp guard, which lists one near-empty directory twice per
  test (5.7 us, 26 ms a run, measured). The ledger-334 tree's `-n 6` inside the worktree
  took 5m18s at ~21-30. Before that, 4380 (ledger 331 on 330):
  serial (`-n0`) **9m22s** and `-n 10` **3m38s**, both under other sessions' load. Before
  that merge, ledger 331's branch alone (4351) took serial (`-n0`) **7m33s** (load average ~3.9 at its start, the 5-minute figure
  ~11 as another workload wound down, and two scripts reading a real run shared its
  first two minutes) and `pytest -n 10` **2m13s** beside another session's serial suite
  (~4 at start), a ratio of 3.4x -- 2m55s when that session's own `-n 10` overlapped
  the last 90 s. The same 4351 took 2m12s at ~14-15 the day before, 4350 1m37s at ~6.4.
  4336 (ledger 318's fix and its
  review): serial (`-n0`) **7m03s** and `pytest -n 10` **1m35s**, the branch's final
  check with no other MF6 workload on the machine (load average ~3.5 at the serial
  run's start), a ratio of 4.4x; an earlier `-n 10` of the same 4336 took 1m32s at ~7
  (4332 took 1m38s at ~1.6, 4297 with ledger 329 1m54s at ~3, 4286 2m50s at ~4, and
  4271 with ledger 328 2m21s). 4266 on
  the same tree took serial (`-n0`) **8m49s** and `-n 10` **1m37s**, both with other work
  on the machine (load average 1.5-4.8 through the serial run, ~10 at the parallel one's
  start), so neither is a quiet clock -- nor were 4261's 7m58s / 1m43s the same day (ledgers 318-327, MF6
  runs of a 20,323-cell model sharing the machine). The last quiet serial
  measurement is **4037 in ≈ 6m29s** (2026-09-24; 3853/6m24s on 2026-09-23). The
  inner loop `pytest -m "not slow"` ≈ 32 s is from 2026-09-19 and was NOT re-measured
  (serially, under the same load, it took 2m02s for 3773 tests).
  **A parallel run competing with another MF6 workload fails intermittently** —
  a different test each time, all passing in isolation, wall time 85 s → 245 s.
  That is contention, not a regression; re-run alone before believing a failure.
  **Read the two numbers as a RATIO, not a clock.** Both grew ~1.5x between
  2026-09-14 and 2026-09-19 on the same tree under load, while serial stayed
  ~3.7x the parallel run — so a wall time well off these is machine contention,
  and it is the ratio moving that would mean something. That 3.7x gap is the
  price of the `-n0` pre-commit check below, not a regression.
  **Verify sign/column changes with `-n0`** — a session-fixture/xdist interaction
  can report green while serial catches real failures (see ledger 48).
- **A relative `PYTHONPATH` is fine: the suite makes it absolute** (ledger 354, which
  retires 344's "give it as an ABSOLUTE path"). Python resolves a relative entry against
  the directory a process STARTS in, and the PEST end-to-end tests run `python
  forward_run.py` with the TEMPLATE as working directory, so `PYTHONPATH=src` resolved
  there, found nothing, and the forward run imported the venv's EDITABLE install -- the
  main checkout -- instead of the code under test (measured in 344: `TypeError` on the
  main checkout's signature). `tests/conftest.py`'s `pytest_configure` now rewrites every
  relative entry against the directory pytest was started in, before any xdist worker or
  test starts a process, and `tests/test_subprocess_imports_the_code_under_test.py` fails
  naming both files if a child imports another myflopy. A reproducer that runs PEST++
  makes it absolute itself, through the same function (`docs/reproducers/
  _absolute_pythonpath.py`; ledger 425) -- a new one calls
  `make_pythonpath_absolute()` at module level or `tests/test_reproducers_make_pythonpath_absolute.py`
  fails naming it. **The rule lives in the library, `src/myflopy/_pythonpath.py`**
  (ledger 454), dependency-free, loaded BY PATH by conftest and that file, and
  `PestProject.run_ies`/`prior` hand PEST++ the absolute value too -- refusing a
  relative entry the process has `chdir`ed away from. Any other script you run BY HAND,
  or a notebook kernel started elsewhere, still resolves its OWN relative `PYTHONPATH`
  against its own start directory.
- **A session writes its temp files in its OWN root, and a leak fails by name** (ledger
  335): `tests/conftest.py` points TMPDIR at `/tmp/myflopy-pytest-*` (one subdirectory
  per xdist worker) and removes it at the end, even after a failure. A test that leaves
  anything in `tempfile.gettempdir()` ERRORS at its teardown naming the entry; what a
  fixture or a process leaves fails the run, with a red section and "temp leaks" in the
  count. Write under `tmp_path`/`tmp_path_factory`, or clean up what you make. /tmp is
  a 16 GB tmpfs with a per-user quota: when it is full EVERY write the user makes there
  fails, Claude Code's own tool output included (commands return "Exit code 1" with
  nothing printed).
- **Nor in the WORKING directory** (ledger 356): the controller lists the directory
  pytest started in before any worker starts, and an entry the run adds there (bar
  `.pytest-work`, `.pytest_cache`, `__pycache__`) fails the run under a red "left in the
  working directory" section. Session-level, not per test: xdist workers share that
  directory. It caught pyEMU's `PstFrom.log` and an empty `s2/` from a
  `SimulationSpec.build()` called with no workspace (which makes `./<name>`, as
  documented) -- pass one under `tmp_path`.
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
`mf.*` build helpers and the `model.<pkg>.<inputs|results>.<noun>.<verb>()`
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
Project            ← durable workspace + run/scenario lifecycle (workspace.py); holds no MODEL geometry
                     (its `aoi` areas of interest are places to LOOK at, ledger 310)
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
  hand-written DISV, which IS the old mesh -- but **an eagerly declared `.gpkg` is
  now `rebuilt`, not `stale`**: it records the recipe it was built from, so the move
  RE-DECLARES it against the new mesh exactly, no overlay involved. Recorded rather
  than re-derived because `metadata['fields']` is three fields short -- `nper`,
  `edges_only`, `boundnames` -- and `edges_only` decides WHICH CELLS are kept, so a
  guess writes a different package in silence. **`mf.uzf` now has a recipe too** --
  declare it without `context=` and it resolves at BUILD from the DOMAIN, reading
  `rebuilt` (ledger 304); `sfr`/`lak` still have none and stay `stale`, ledger 296), `unverifiable` (names cells with no pin --
  including the ARRAY forms, which name cells by POSITION and so carry no cellid to
  find, ledger 282; **and any package holding a per-cell ARRAY** -- `npf`/`ic`/`sto`
  take no `context=`, so nothing stamped them and they read `unchanged` while a lazy
  disv rebuilt correctly around them and the build carried K onto the WRONG MESH,
  ledger 295 --
  or an EDIT voided the pin -- ledger 275: `packages.edit`/`fork`/`with_options`
  replacing a certified option makes the library's certificate a lie, so it is voided
  and the refusal names the edit instead of telling you to re-declare what you just
  typed),
  `unchanged`. **It reports and stamps; the BUILD refuses**, so read
  `sim.move_report()` before `prepare_run`, and `run.resolution` after.
  **What the new mesh would have NO DATA for is `model.coverage_on(grid)` /
  `sim.coverage_on(grid, on=)`** (ledger 298) -- pure geometry, so the same number
  whatever moves. Reported as **AREA IN MODEL UNITS**, not a coverage fraction: polygon
  intersection is exact to ~1e-15, so on a mesh refined over the SAME footprint
  `coverage < 1` flags 765 of 2,352 cells (32.5%) short by 4e-15 of themselves, where
  the same residual in area is ~1e-11 ft2 and reads as zero. A footprint grown 1,200 ft
  east reports 3,600,000 ft2, which is checkable by hand. **A DISJOINT pair is a
  separate verdict, not "100% uncovered"** -- `grid_pin` is connectivity-only, so the
  same mesh in another frame (what a USG local origin produces) pins IDENTICALLY while
  overlapping nothing, and calling that a gap is what would invite filling every cell
  from a source 1.5e6 ft away. A model with no grid REFUSES rather than answering. `on=None`
  moves EVERY model with a grid, because coupled models must share a mesh (ledger
  268); `on=` is the spelling that can strand a sibling. `grid_free=` is the caller
  asserting mesh-independence. The stamp is the only thing behind a hand-typed cellid
  across a move -- every number stays in range on a larger mesh. **A DEFERRED-grid
  model has no old pin to stamp**, so its hand-written cellids move unrefused; they
  were never verifiable anyway -- and a per-cell array cannot be measured there
  either. The array detection is ncpl-EXACT (a per-LAYER list `k=[10.0, 1.0]` means the
  same on any mesh, so guessing from "is it a sequence" would refuse a correct move),
  and `with_grid` supplies the count from the OLD grid it STILL HOLDS at the call site
  -- it does not mutate, so the mesh is gone from the RESULT and present at the CALL.
  The shared rule is `grid_identity.column_count` / `cellwise_options`, in import layer
  1 because three callers need it; it tells a scalar and a per-layer list from a flat
  `ncpl` list, an `(nlay, ncpl)` array and a NESTED per-layer list of per-cell lists.
  It gave the arrays a VERDICT and not a route; **`remap=True` is now the route**
  (ledger 300, below). `replace_grid` is a warned `__getattr__`-only alias.
  **A MOVER is classified from its SIBLINGS, not from itself** (ledger 303): `mf.mvr`'s
  records name a provider by package name plus an INTEGER -- a reach or lake number -- and
  none of the mover's own options change with the mesh, so it read `unchanged` beside a
  `stale` SFR while holding an index into the old reach network (measured: outlet at reach
  9 of 10 on one mesh, 20 of 21 on another). `package_moves(packages, ...)` is the
  two-pass classifier and **all three** callers use it -- `move_report`, `with_grid`'s
  stamping loop and `Run.resolution` -- because updating two leaves the run table and the
  pre-build report contradicting each other. A provider that is anything but `unchanged`
  invalidates the index, **`rebuilt` included**: a rebuilt SFR renumbers every reach, so a
  CORRECT move of the provider is what makes a typed index wrong. The marker is
  `INDEXED_PACKAGES`, stamped in `mvr_spec` (the funnel both entry points reach), and
  deliberately NOT `requires=` -- that means "build this first", so reading it as an index
  reference would refuse a package ordered after a discretization whenever the disv
  rebuilt (`mf.hfb` carried `("disv",)` until ledger 378). It is conservative about WHICH provider counts (a lake number is
  feature-derived and usually survives) because a false refusal costs a re-declaration
  where a false pass routes water to the wrong reach.
  **A gridless `mf.uzf` is marked, not pinned** (ledger 302): it resolves from the IDOMAIN
  and needs no grid, so there was nothing to pin and it read `unchanged` -- a 25 -> 81 move
  BUILT CLEAN with 75 UZF objects on columns 0..24 of 81. It carries `CELLWISE_OPTIONS`
  now, the same two-branch shape `_array_bc_spec` uses. The marker is read by
  `package_move`, so it closes the `with_grid` path and NOT a hand-assembled model on the
  wrong grid -- pinned by a test so the claim stays honest.
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
  (SFR's `packagedata[1]` is a cellid, LAK's is a stage).
- **`mf.rch(...)` / `mf.evt(...)` declared WITHOUT `context=` resolve at BUILD time too**
  (ledger 284, 2026-09-18) -- measured 28 records at ncpl 28 and 92 at ncpl 92,
  byte-identical to the eager form on each. The recipe is FLAT builder kwargs in
  `inputs`, recording only what the caller SET. **Its precondition is WIDER than the
  `.gpkg` one**: the areal builders touch only `context.domain` and `grid.ncpl`, both of
  which the built DISV carries, so `_lazy_resolution_context(needs_grid=False)` keeps the
  lazy form from being NARROWER than the eager one -- it accepts a bare
  `ModelContext(domain=...)` and a model with no context at all, and synthesizes an
  all-active domain when the DISV declares no `idomain=`. **Caller-typed cell identity is
  REFUSED in a lazy recipe** -- an explicit `cells=` list, `name_by_cell=`, a
  `{cellid: value}` mapping, a per-cell sequence -- all still available with `context=`.
  Two traps it introduced and the mutation check caught: binding the partial to
  `RCHBuilder` narrowed the option check so the lazy branch refused `print_input` (a
  `lazy != eager` divergence a typo test cannot see -- assert BOTH directions), and a
  `package=` KEYWORD on the partial moved the digest, because `package_digest` pickles
  `name` and the builder reference into one memo table and an interned string collides
  (168 bytes vs 172; `docs/reproducers/digest_pickle_memo.py`).
- **A `.gpkg` list BC declared WITHOUT `context=` resolves at BUILD time** (ledger 279,
  2026-09-17) -- one declaration builds correctly on two meshes (measured: chd 11
  records at ncpl 36, 13 at ncpl 100; drn 10 and 7 carrying a `Spread` and a
  `CellSurfaceOffset`), and `move_report` reads `rebuilt` rather than `stale`. **With
  `context=` nothing changes** -- the eager branch is untouched and still resolves at
  declaration, so the "re-declare its GIS packages to move a model" advice still holds
  for anything declared eagerly. The recipe is a `BoundarySource` in
  `PackageSpec.inputs` under `stress_period_data` (ledger 277); the builder is bound to
  its FloPy class with `functools.partial` so the typo check stays byte-identical to the
  eager form. **The discretization must be declared BEFORE it** (refused by name
  otherwise, found by `package_type` -- `model.get_package("dis")` prefix-matches a
  DISV), and the domain and layer elevations fall back to the BUILT discretization when
  the context carries none, which is the normal case after `with_grid` since that drops
  both. Reading the idomain off `model.modelgrid` instead is wrong: flopy synthesizes an
  all-ones float array where the package's own field says None.
  **The run freezes what it BUILT** (ledger 280): the snapshot holds records and the
  new mesh's pin, so a lazy run rebuilds with its GeoPackage deleted -- inline AND
  through `mf.ref`, which failed where inline succeeded because `freeze._for` keys on
  the DECLARED entry's name and a `PackageRef`'s name is its library key. The declared
  spec keeps its recipe and stays `rebuilt`; the snapshot is pinned and reads `stale`,
  which is correct for something already resolved. `Run.resolution` reports
  `resolved_ncpl` for a discretization and `resolved_records` for a boundary condition.
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
  `SaveReport.bundled`. **Grids AND lazy packages** (ledger 293, 2026-09-19): a
  package declared without `context=` holds a RECIPE, and a lazy model has **TWO
  carriers, of which the easily-forgotten one fails FIRST** -- the `.gpkg` recipes are
  `DataSourceSpec`s an existing walk would find, while the gridless `LayerStack`
  behind `mf.disv` is not a dataclass and names no source spec, and its raster is
  opened inside the DISCRETIZATION, before any BC is built. So a BC-only bundler
  reports a complete archive and ships one that cannot open. An EAGER package still
  carries nothing, and that is measured, not assumed: it resolves at declaration and
  holds records -- zero file reads inside `PackageSpec.build` across 306 tests, and a
  project rebuilds with its GeoPackage deleted. A `GridSpec` is the
  opposite: `ModelSpec.build` re-resolves it on EVERY build. `rasters=` and
  `idomain_path=` are bundled too even though they are read later still (the first
  `gdf_topbtm`, the first edge query). An `external=True` source IS bundled and the
  flag cleared -- it meant "never anchor", which stops being true once a copy is
  inside. **A run now records its project root**, because bundling stores RELATIVE
  paths and `Run.reopen`/`load_run` had no `project_root`, so a relative source
  resolved against the current working directory. The PACKAGE tier needed the same
  anchor, in `ModelSpec.build` and **after `resolved_packages`** -- before it, a
  `PackageRef` carries no `inputs` and every `mf.ref` package is silently skipped
  (ledger 280's miss, one function over). `Surface.out` IS carried: it is the GRASS
  cache, so an archive without it needs GRASS to rebuild, and with it the cache reads
  `stale` and is reused rather than regenerated. `SaveReport.not_bundled` names what
  was left behind (a stack holding per-cell arrays cannot be serialized at all).
- **A lazy package's digest reads its file's CONTENT, not its path** (ledger 293).
  Path-keying was wrong in BOTH directions at once: bundling moved every lazy digest
  on a byte-identical archive, so `runs.drift()` called every stored run edited the
  day it was archived; and a DIFFERENT file at the same path digested IDENTICALLY.
  `BoundarySource.source_sha` rides in `inputs`, so it is digest-VISIBLE, and it
  REPLACES the path rather than joining it -- keeping both reintroduces the churn
  while looking fixed. It falls back to the path when the file could not be read.
  **Both carriers do this** -- `BoundarySource.source_sha` and `Surface.source_shas`,
  the latter stamped at LAZY-PACKAGE declaration (`mf.disv`/`mf.evt`) so a bare
  `Surface.raster(p)` still opens nothing. `Surface.out` is dropped from the digest
  entirely: a WRITE target whose content is decided by what IS digested. No grid
  source moved -- `digest_payload()` defaults to `to_dict()`.
  **A lazy `.gpkg` RESOLVES its path at declaration**, because the declaration checked
  it against the cwd while the build anchored it to the project root -- two bases for
  one string, and with a same-named file under the root the model BUILT CLEAN from the
  wrong one (999.0 against a declared 9.0, no warning). The anchor now only ever sees
  bundler-written project-relative paths. This also makes lazy AGREE with eager, which stores records
  and is therefore location-blind and content-sensitive already -- the two had been
  drifting in OPPOSITE directions from the same declaration. `drift()`'s `source`
  column was inverted for the same reason (only the eager form stamped
  `source_path`), so the package that cannot build without its file read `—`.
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
  boxes; a partial miss WARNS by name -- through **`warnings.warn` as well as the
  logger** since ledger 288, because `import myflopy` installs a `NullHandler` and root
  has none, so `logging.lastResort` never fires and the log line reached NOBODY (zero
  bytes on stderr, measured in a bare session and in a real ipykernel). **`gaps=` says
  what a miss means** -- `"warn"` (default) / `"nearest"` / `"raise"` / `"ignore"`, the
  same four words as the raster sampler, validated through one `check_gaps`. The
  all-missed refusal is ORTHOGONAL to the mode, `"ignore"` included: the raster side
  lets `"ignore"` through because `Surface.sample` has a next stage, and a boundary
  package has none. **`gaps="nearest"` moves POINTS only** -- translation recovers a
  point exactly, recovers a LINE only when it happens to lie parallel to the boundary it
  missed, and leaves a polygon touching at 0.0001 of itself; orientation luck is not a
  rule. An areal feature is CLIPPED to the part inside the grid (which `cell_shares`
  already did exactly) and the clipped-away fraction is now REPORTED.
  **`snap_within=` is `None` in the recipe and resolved at BUILD** against the candidate
  cell's own `sqrt(area)` -- a number resolved at declaration carries the old mesh's cell
  width across `with_grid`, and cell areas span 3.9-18.5x within one mesh. Compared
  against POLYGON distance, not centroid (3.07x apart on one real miss; the raster
  sibling uses centroid). A total of what was NOT written is reported only where the
  library knows the field is extensive -- a well's rate, or a `Spread(mode="clip")` field
  -- because an unwrapped conductance is copied whole to every covered cell and
  declared-minus-written would call a clean build a -300% loss.
  **The blind spot to state out loud:** a unit/scale CRS error leaves near features ON
  the grid in the WRONG cells. Measured, EPSG:32149 against a 600 ft model at the origin:
  2 of 5 wells landed 161 ft and 484 ft from where they were drawn, 60% of the pumping
  gone. No tolerance can see that. A polygon covers a cell by AREA (ledger 210):
  drawn along faces it reaches only the cells inside it, never its neighbours. A LINE
  piece under the same 0.01 tolerance is noise too (ledger 237) -- and this is the
  SHARED rule, so it moves `map(select=<line>)`, a region and a line head target as
  well as a `.gpkg` source: a line ending on a face and overshooting by 1e-6 ft -- a CRS
  round trip moves an end by 3e-7 -- wrote a full second record, where `mf.sfr` had
  dropped the reach since ledger 214.
  **A multi-part POINT divides itself between its parts** (ledger 287): every part used
  to get a share of 1.0, so a MultiPoint of 3 wells declared at `-1000` wrote **-3000**
  under `mf.Spread` -- which exists to conserve a feature's total -- and a part off the
  grid vanished into a sum indistinguishable from a complete smaller feature, never
  entering `missed`. Lines and polygons were always right (they divide by their own
  extent), which is what isolates it to the zero-extent branch. It CHANGES WRITTEN
  VALUES for a multi-part point under `Spread`; an unwrapped value is untouched.
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
- **An elevation outside its own cell: `min_above_bottom=` on chd/ghb/drn/riv** (ledger
  290, 2026-09-18). **MF6's behaviour splits on `icelltype`, and BOTH branches are live
  here.** At `icelltype = 0` — what `mf.npf(k=10.0)` writes, naming no icelltype at all —
  all four ACCEPT a below-bottom value and over-drain under Normal termination (2.17x on a
  bench, 1.09x on a real 64-cell model). At `icelltype != 0` — `CanonicalModelConfig` and
  17 taught sites — MF6 REFUSES by name, NEWTON does not change it, and the refusal is a
  PER-PERIOD read, so an offender in period 3 dies after periods 1-2 have solved.
  **The DETECTION is default-on and the correction is opt-in**, because what is constant
  across models is the silence, not the magnitude: nothing reported this before, and the
  count/worst depth/lift now ride beside `coverage` in metadata.
  **Half the arithmetic already shipped** and the entry says so:
  `CellSurfaceOffset("cell_bottom", offset=F, minimum=<declared>)` IS `min_above_bottom=F`
  for a constant or single-column elevation. What it cannot floor is a per-period list or
  a surface-relative value (`minimum=` is per FEATURE, and `minimum=<CellSurfaceOffset>`
  raises) — and surface-relative is the idiom the docs teach (63 `model_top` + ~50
  `cell_top` against 5 uses of `minimum=` in the whole tree).
  **`min_above_bottom=0.0` cannot mean "at the bottom"**: FloPy writes `botm` fixed-point
  and a BC value in E format, so exactly-equal can round BELOW and MF6's strict compare
  trips — 13 of 24 real runs FAILED, round numbers being the immune cases — and a CHD
  exactly there aborts the solver. Nudged by `max(|bottom|*1e-6, 1e-6)`; a plain absolute
  1e-6 is not enough.
  **RIV is a PAIR**: MF6 checks `rbot` and never names `stage`, and flooring `rbot` alone
  trips `stage >= rbot`, turning a running model into a hard failure. Floor rbot, raise
  stage by the same delta (measured: the three readings are 3.5x apart in flux; the naive
  ones collapse a declared 1.0 ft depth to 0.000). A declared `stage < rbot` is refused.
  **The floor runs AFTER the `gaps="nearest"` snap** — a snapped point is in a different
  cell with a different bottom.
- **A `.gpkg` NULL raises too, naming the row** (ledger 247, 2026-09-13, amended by 289).
  **Two members of that family survived until 2026-09-18** (ledger 289): a NULL in a
  `CellSurfaceOffset`'s `minimum` column read back as `nan`, `max(value, nan)` returns
  `value`, and the result is FINITE -- so the guard passed and the FLOOR SILENTLY DID NOT
  EXIST (measured: output byte-identical to no floor, 11 of 11 records below their own cell
  bottom). And **`Spread` on an ELEVATION** was accepted and divided it by the cell's area
  share (a declared 75.0 written as 2.06-11.87) -- ledger 274 refused that on the ARRAY
  form and never reached the list form. `_ELEVATION_FIELDS` names the six fields, and the
  refusal is called from BOTH the eager and the lazy path because a lazy declaration never
  reaches `_bc_from_features`. A null
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
  **A pond is `stack.excavate(bottom, within)`** (ledger 320) -- a VOID, cut in
  `build()` AFTER reconcile and pinch (cutting before them was measured wrong under
  every policy). The connected cell's top is EXACTLY the bed; cells above go idomain
  0, never -1; a sliver under the bed is ABSORBED into the cell below (`min_remnant`,
  0.1), never deepened. `within` is read, covered and sampled through the SAME helpers
  as `mf.lak` (`geometry.read_footprints`/`footprint_cells`, `surfaces.sample_bed`), so
  record cells == lake cells; `LakeTableBuilder(footprint=)` reads through
  `read_footprints` too, reprojected to the RASTER's CRS (ledger 318). **A lake with NO
  bathymetry takes `bottom=mf.Basin(lakes, depth, water_surface=)`** (ledger
  554): 0 at the shore, `depth` at the pole of inaccessibility, `"bowl"` by
  default (1.5x a `"cone"`'s volume -- the user's call), each cell the AREA MEAN of the
  profile (a centroid read missed the volume by 10%), NaN off the footprint -- so as a
  stack top or unit bottom it RAISES (MF6 ran that NaN to normal termination, every head
  NaN; `mf.Min(<surface>, bed)` composes; ledger 555); the
  excavation footprint and `Surface.basin` share ONE reader (`_coerce_footprint`). The result carries `excavations` (`mf.ExcavationRecord`),
  `ground`, `excavated`, `absorbed`; `attach=True` publishes `grid.excavations`, which
  the `gdf_topbtm` SETTER clears. `stack.to_disv` is now `build().to_disv()` (0 digests
  moved, 12 stacks). **An excavated stack SAVES** (ledger 323): excavations ride in the
  `$layer_stack` payload only when there are any (16 of 16 un-excavated lazy-DISV
  digests unmoved), the footprint FILE under the `"zone"` key `Surface.where` uses -- so
  the bundler, the build anchor, `built_from` and the hash stamping find it with no
  wiring -- and an `mf.Area`/geometry inline as WKT. `from_dict` builds entries
  DIRECTLY (reading the file there checks a bundled relative path against the cwd) and
  REFUSES unknown keys -- but runs every OTHER refusal `excavate()` makes, through the
  same helpers (unchecked, a repeated name dropped one pond's RECORD while its cut
  stayed). The lazy `mf.disv` publishes a table built from the RESULT, not
  `grid.gdf_topbtm`, with the records in `attrs["excavations"]`, marked as the DISV's
  own so `_layer_elevations` does not blame the context when `attach=False` leaves the
  grid's table behind. The frozen snapshot carries the records too, onto the thawed
  table and grid. **`mf.lak` reads the record** (ledger 325) through the SAME resolver
  as the elevations (`specs._layer_elevation_source`: a `LayerBuildResult` as `surfaces`,
  else the grid's `excavations`), never from geometry -- a geometric test was measured
  wrong three ways. A lake whose cells ARE a record's gets the EXCAVATED rule: its
  vertical goes to the top cell with idomain > 0 AND thickness (activity alone picks a
  zero-thickness layer when the context has no `domain`), its walls go INTO the
  neighbours' columns with MF6's own DISV `CONNLEN`/`CONNWIDTH`
  (`grid/connectivity.DisvFaces` -- the centre-to-centre spacing measured 0.93 of a pit's
  through-flow, this 1.005, the old rule 0.39, 0.52 since ledger 328), and `lake_bottom`
  defaults from the record. A partial match RAISES whatever `unexcavated=` says. Every
  other lake keeps the old rule (`tests/data/lak_legacy_golden.json`, written by the old
  code) except faces under 1e-6 and CONNLEN (below), and is WARNED about
  (`unexcavated="warn"`; the canonical fixtures filter it, F8).
  **`grid.cl12` is MF6's PER-SIDE length** (ledger 328, the user's call reversing 325's
  byte-identical legacy): each cell's DISV centre to the LINE of the shared face, from
  the ONE rule `DisvFaces` (`DisvGeom.f90`'s `shared_edge`/`distance_normal`, the pairs
  and the edge from the LOWER cell's ring, as MF6 builds them). From 2025-07-02 it held
  the full centre spacing: every DISU myflopy wrote carried ~half the DISV's horizontal
  conductance (a DISU now matches its DISV to 2e-10 ft, `docs/reproducers/
  disu_cl12_matches_disv.py`), and a legacy lake wall's CONNLEN was ~2x MF6's (1.74-2.35
  on real ponds). DISV models were never affected -- MF6 computes its own. Built from the
  vertex rings, not polygons: 0.08 s against ~8 s at 18,107 cells, same pairs on every
  mesh measured. **A Voronoi mesh cannot test WHICH side a length is on**: its face is the
  generators' perpendicular bisector, so cl1 == cl2 to rounding (the canonical grid to
  5.1e-7 ft) and a table with the sides SWAPPED passes -- test per-side lengths on quads
  whose DISV centres sit off their centroids, as `tests/test_disu_matches_disv.py` does.
  `calculate_distance` is NOT this number (polygon centroid to the shared boundary,
  0.48-15.9x `cl12`). `grid_identity.lak_geometry_complaints` checks the
  connections against the BUILT DISV from the LAK builder AND
  `_discretization_complaints`; `lak_build_warnings` says no-NEWTON, `icelltype` 0
  walls, RCH/EVT under a wet lake, and a transparent bed -- for excavated lakes only.
  **`mf.lak(implicit=True)` writes MF6's LAK `IMPLICIT`, which FloPy 3.10 cannot** (ledger
  326): `lak_implicit.ModflowGwflakImplicit`'s own `write` inserts it after FloPy writes,
  so every write of THAT OBJECT gets it once. **The flag is the CLASS, so anything that
  REBUILDS a LAK package loses it** unless it carries the class: FloPy's splitter
  (`package_factory`), an artifact restore, FloPy's loader, and `model.config` (no FloPy
  dataset to read) each did, silently -- `split_model`, `apply_package_artifact`,
  `keep_implicit` and the config's `implicit` row carry it now; a bare FloPy load still
  cannot. Two traps: FloPy's `MFPackage.__init_subclass__`
  registers ANY subclass as its LOADER's class for that type -- defining one made every
  LAK loaded afterwards the subclass, IMPLICIT or not; the registry is restored right
  after the class. A CG solution is refused at build (MF6 stops; `complexity="SIMPLE"`
  or none defaults to CG, measured) -- except the IMS `mf.simulation` generated, which the
  build turns BICGSTAB once an `mf.ref` lake resolves (`GENERATED_SOLVER`) -- and the
  executable is probed once per file at run time.
  `docs/reproducers/lak_implicit_requirements.py` runs MF6 on each claim.
  **Why a CUT and not a better connection record** (ledger 327): MF6 overwrites a
  VERTICAL connection's BELEV/TELEV with its cell's TOP, so a bed inside a cell exchanges
  nothing through the floor until the stage passes that top -- with a table, NEWTON or
  `IMPLICIT` alike (`docs/reproducers/lak_vertical_floor.py`, 26 cases against MF6's own
  formula, then the same bench excavated). On the real pond 411 of 612 carried exactly 0
  all run; after the cut a stage-volume table is OPTIONAL (MF6's own curve within 1% of
  the contour table, 530-533 ft). **End to end through the API on that pond** (ledger 327)
  the default formulation reproduced the prototype's run to every digit the ledger shows
  at +80% outer iterations and a 2.16% worst LAK step; `implicit=True` took +14% with the
  lake budget closed, which is why the teaching passes it. MF6 writes the pond's stage as
  -1e30 on every dry recovery day (D10); its stage view leaves a gap there (ledger 329). **A model with NO stack** cuts a `LayerBuildResult`
  built from its loaded arrays (`result.excavate` returns a NEW one --
  `mounding_L3K_to_lak.ipynb`), and its loaded IRCH must become `"top_active"` or the
  build refuses the dug cells it names (ledger 321).
  **A package declared against a context from BEFORE the cut is refused at build**
  (ledger 321): `_discretization_complaints` checks every list-BC record, OBS6 cellid
  and per-column `irch`/`ievt` against the BUILT DISV's idomain (`<= 0` = removed) and
  names the excavation. FloPy does NOT catch a list BC on a removed cell (measured: a
  CHD at a dug cell built clean); HFB is exempt (MF6 warns and drops it). A spec is
  paired with the FloPy package it BUILT, not by name (a bare `PackageSpec` is `chd_0`
  to FloPy), and a LAZY package's refusal names the MODEL's stale `domain=`. `irch=None`
  / one layer everywhere stays the opt-out and WARNS with the rate lost, which matches
  MF6's own budget to 1e-6 for recharge and is an UPPER bound for ET (ledger 322).
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
**The lazy spec holds its OWN copy of the stack** (ledger 270, amended 2026-09-25) --
a SNAPSHOT: `stack.add(...)` after `mf.disv(stack)` does not reach the declared
package, a library copy or its digest; declare again to take an edit. It held the
caller's object whenever nothing was hash-stamped (an all-inline stack). The copy is
`from_dict(to_dict())` -- what a reload builds -- so it digests the same (124 of 124
across the suite); `stamp_source_hashes` always rebuilds now, `$tuple` tag removed first.
**So does the EAGER `mf.disv(result)`**: it held the result's own arrays (`np.asarray`
shares memory), so `result.idomain[...] = 0` after the call edited the package; it
copies now, and a result given a `grid=` of another cell count is refused at the call.

### Grid: eager (works now) vs deferred (a known seam)
- **Eager** (use for the full GIS stack): build the grid object first (`VoronoiGridPlus`, or
  `mf.GridSpec.voronoi(...).resolve(workspace)`), put it in `ModelContext(grid=vor, domain=idomain)`,
  then declare packages. **Required today** for `mf.sfr`/`mf.lak` because they resolve cells
  eagerly at declaration. **`mf.uzf` no longer needs it** (ledger 304): omit `context=` and
  it resolves at BUILD from the domain -- it never used the grid for anything but its own
  pin (measured, one attribute: `get_disv_gridprops`), so it needs strictly LESS than
  `mf.rch`. Refused lazily, each still fine with `context=`: an explicit `cells=` list, a
  per-cell value, and `mover=True` -- a UZF mover endpoint is a hand-typed `iuzno`, an
  enumeration over ACTIVE CELLS that is in range on every mesh and means a different column
  on each (measured 18-20 against 33-35), and no `uzf_connection` exists to make it
  symbolic.
  **The seven `.gpkg` forms no longer need it** (ledger 279), and neither does
  `mf.rch(...)`/`mf.evt(...)` (ledger 284) — omit `context=` and they resolve at BUILD
  time against whatever grid the model supplies. Passing `context=` still resolves and
  pins at declaration, unchanged.
- **Deferred** (`mf.gwf().with_grid(mf.GridSpec.voronoi(boundary=<gpkg>, refinement=<gpkg>))`):
  the project builds the grid at run time into `run.workspace/_grid/<model>` and populates
  `context.grid`. In `ModelSpec.build` the grid resolves BEFORE packages build, and the model
  carries `model.myflopy_context` — so deferred GIS packages are *feasible*, just not wired: the
  helpers would need a grid-lazy mode emitting a `PackageSpec` whose `build(model)` resolves against
  `model.myflopy_context` instead of resolving eagerly. **That mode now exists** for the seven
  `.gpkg` forms and for `mf.rch`/`mf.evt` — omit `context=`. Deferred GridSpec therefore composes
  with disv, the list BCs and areal recharge/ET; still not with `mf.uzf/sfr/lak`.
- **A grid built from a triangulation is REPAIRED to tile its boundary** (ledger 316).
  `VoronoiGridPlus(tri)` runs FloPy's `tri2vor`, checks every cell against the dual of
  Triangle's own mesh (`grid/tessellation.py`) and rewrites only the wrong ones, on
  FloPy's vertex ids: a grid with nothing wrong is byte-identical and keeps its pin; one
  that had a broken cell (every 270 deg reflex corner did) re-pins with no cell
  renumbered. **Every returned ring IS the dual's ring under one relabelling** -- each
  dual vertex gets FloPy's id at its location (to 1000 ulp) or an appended one, one to
  one -- so no face can be lost; the first version kept FloPy's ring unless a detector
  objected, and a fuzz pass found faces it missed (a Qhull-merged vertex accepted at 5% of
  the circumradius, FloPy's 180-degree filter re-applied). Circumcentres within 1000 ulp
  of the coordinates are ONE vertex (an edge of two ulps turned over under reprojection
  and GEOS refused the grid), so two NODES that close, or a node whose circumcentres all
  merge, have no cell left and are REFUSED by name rather than dropped -- a point
  constraint on a boundary edge, or within rounding of one, did that, and is now PLACED on
  the edge first (`TriangleGrid._place_point_constraints`, 0.02 ft at state plane; ledger
  316's five). The certificate
  stays at 1e-9 of the edge: widened with the merge it passed a boundary circumcentre
  outside the domain and returned a cell crossing itself. **A corner is never split
  closer than the FILE records** (2026-09-28): FloPy writes a coordinate above 1e5 to
  nine significant digits -- 0.01 ft at x = 1.26e6, `grid_identity.written_resolution` --
  and four nodes almost on one circle give two circumcentres a hair apart. The repair
  kept both, 7.4e-5 ft apart, where FloPy drew one corner; the face between them was
  written with ZERO length and a real 20,323-cell model ran to normal termination with
  every head NaN. Dual vertices joined by a face shorter than one written unit on both
  axes are now one corner -- FloPy's where it drew one -- counted in `unsplit`; a face
  FloPy drew ITSELF is joined only when it is a FACE two cells share and both its ends
  PRINT as one point (a plain four-corner boundary at state plane on the default profile
  had one: NaN heads on FloPy's own cells), so no grid that runs today moves. And a
  model whose DISV would still write a face MF6 cannot compute -- zero length, or both
  cell centres written on its line (CONNLEN 0, a division by zero) -- is REFUSED at build
  naming the cells (`unwritable_faces`, comparing the numbers MF6 reads back, not the
  text; a snapshot rebuild is the grid that ran). Measure what the
  FILE holds, not what float64 holds: every precision check in this module was right in
  memory. `tessellation_report` says what
  changed and carries `legacy_pin` (FloPy's own pin), and a package pinned to that is
  refused by a message naming the repair; one resolved under the FIRST version of the
  repair (9bd3651) has no pin on record, and its refusal says so; one resolved on a grid
  the repair split from bb56668 matches `split_pin`, and its refusal says the cells kept
  their numbers. Where FloPy's cell SET
  is wrong -- a boundary in two parts (ledger 317), a point outside the domain -- the
  dual's cells are taken whole, renumbered, logged at WARNING. `tessellation="flopy"`
  (on `GridSpec.voronoi` or `VoronoiGridPlus`) rebuilds FloPy's cells exactly, and a grid
  rebuilt from verts/iverts (snapshot, pickle, reloaded run) is never repaired. A
  boundary with an interior ring builds (ledger 317): `TriangleGrid` gives Triangle each
  hole as ring plus interior point, never a region's ring as a hole, and refuses a
  boundary that touches itself (a hole meeting its edge, parts meeting at a corner).

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
**"Top active" is `idomain > 0`, one rule: `grid/geometry.top_active_layer`** (ledger
319) -- array RCH/EVT, the areal builders, UZF, SFR's `reach_layer` and the USG
importer all call it. The array form read `!= 0`, so a pass-through `-1` counted as
active and MF6 dropped that column's flux: 42.2% of a `pinch="passthrough"` stack's
recharge, Normal termination. An EAGER array package on such a domain moves its digest
(0 of 25 stored, 0 of 47 in the suite); a lazy one does not. Do not add a sixth copy.
`nseg > 1` and `boundnames` RAISE rather than being silently dropped — segmented ET
has no array form, because NSEG lives in DIMENSIONS which MF6 never reads under
READASARRAYS. **`mf.Spread`/`mf.CellSurfaceOffset` raise here too** (ledger 274):
both describe a per-FEATURE resolution and the array form has no features, so
neither could act — they were accepted, pickled and digested before. That refusal
saw only the VALUE fields until ledger 283, so `irch=mf.Spread(...)` was accepted
and digested; it reads every keyword now.
**Both array helpers go LAZY when the part that needs a grid is left to them**
(ledger 285) — `irch`/`ievt` at `"top_active"`, or an omitted EVT `surface`, with no
`context=`. That is what finally makes an array-form model SAVEABLE: the default
resolves to an `(ncpl,)` array, `to_dict` refuses an ndarray, and `Project.save()`
persisted **zero simulations** (ledger 270's defect one package over). **The trigger is
narrower than the list form's on purpose**: `irch=0 | None | <array>` with no context is
legal and digestible today, so keying on `context is None` would move those digests and
make `runs.drift()` call every reloaded package edited. Measured 0 of 10 moved.
**`mf.evt.array` with `surface` omitted resolves the MODEL TOP**, and refuses when a
`context=` cannot answer (ledger 283). A typed `ievt` beside it was DROPPED until ledger
324 (`ievt=1` wrote no IEVT; MF6 used layer 1) -- the lazy branch builds its options from
the value fields and keywords, so anything else must be carried by hand. It used to leave the argument out entirely, so FloPy's own
`surface=0.0` default wrote `CONSTANT 0.00000000` — the datum — and MF6 then draws the
full rate everywhere, because every head is above the ET surface, making `depth=` inert:
measured EVT total **-151.13 against a correct 0.0**, Normal termination, no warning.
Same pass: **`mf.evt(context=...)` silently dropped `nseg` and `auxiliary`** — both are
NAMED parameters, so neither ever reached `**options`, and `mf.rch` kept `auxiliary` only
because it has no parameter of that name. `auxiliary` is forwarded now, `nseg > 1` raises. Arrays are opt-in: for a sparse BC they are a pessimization.

## Moving a CALIBRATED array to a new mesh? `with_grid(remap=True)` (ledger 300)
The overlay remap engine, B3's last commit. It intersects the OLD grid's own cell polygons
with the new grid's and carries per-cell arrays across on the resulting areas, so
`unverifiable` becomes a fifth verdict, **`remapped`**. Serves `npf` `k`/`k22`/`k33` and
`ic.strt`/`sto.ss`/`sto.sy`.
**OPT-IN, and the default does not move, for a measured reason: a remap is the WORSE route
wherever a source exists.** Re-declaring a polygon map on the new mesh (ledger 299) held to
7.1% worst-case boundary-flux error where the best reducer managed 12.5% and the worst four
orders of magnitude. This is for an array with NO recipe -- a PEST-calibrated K field, an
array read out of an imported model -- and the warning says so.
**K HAS NO SAFE DEFAULT, and no per-field declaration could ever give it one** (so
`FieldSpec.remap`, reserved in ledger 299, is the wrong shape). The right average depends
on FLOW DIRECTION, a property of the SOLUTION: measured in MF6 on one 2:1 coarsening, with
the contrast ACROSS the flow the arithmetic mean gave **264x** the fine-mesh flux and the
harmonic was nearly exact; rotated 90 degrees the arithmetic was EXACT and the harmonic
wrong ~250x. **The ordering reverses.** `rule="geometric"` is the default (never best,
never catastrophic; the classical answer for a heterogeneous 2-D field) and the engine
computes arithmetic and harmonic ALONGSIDE it, WARNING which cells diverge by more than 2x.
**TWO DENOMINATORS; only the intensive one is built.** Intensive divides by the new cell's
COVERED area -- which is what makes a partly covered cell right rather than diluted
(constant 7.0 field: 7.000000 against 0.0066-6.9878 for the cell's own area). Extensive
divides by the OLD cell's area and is DEFERRED, because `kind == "extensive"` holds a RATE
*and* a TOTAL needing opposite treatment (measured: "split and renormalize" halves RCH on a
2:1 refinement; the area-mean DOUBLES a WEL total). **`npf`/`ic`/`sto` are not in
`_PACKAGE_EXPLORER_SPECS` at all**, so no `FieldSpec` exists for `k` or `strt`;
`grid_identity.cellwise_options` is what enumerates them.
**Integer flags are refused by name** (`icelltype`, `iconvert`, `irch`, `ievt`): the mean
gave 5 of 25 cells values like 0.3851, and `irch`/`ievt` are zero-based LAYER INDICES.
**THE FRAME GUARD HAD TO BE BUILT -- `CoverageReport.disjoint` IS NOT IT.** Measured: two
identical footprints declared EPSG:2927 and EPSG:2926 give `disjoint=False,
uncovered_area=0.0` and `summary()` says "every cell of the new mesh has data under it".
`disjoint` catches only the zero-pieces (translated) case, so
`refuse_a_different_frame` compares CRS FIRST, then `disjoint`. A grid is REFUSED, never
reprojected -- unlike a `.gpkg` source -- because both meshes are the modeller's own.
**Blind spot, stated out loud:** a unit rescale keeping the same declared CRS passes both
tests and is indistinguishable from a legitimate SUBDOMAIN refinement. A test pins the
limit. Same class as ledger 211's.
**`CoverageReport.uncovered` IS NOT THE FILL LIST** -- it names cells that are SHORT, one
measured at 98.9% covered, so filling from it threw away 15 correct values of 51 named. The
predicate is `~reached`. And **mask the WEIGHTS per LAYER, never pass a NaN value through**:
ten reached cells came back NaN, three ACTIVE on the new mesh.
**Write the arrays with `replace`, NOT `with_options`** -- that voids the pin and
`package_move` then says "its records were edited", a lie about a correct remap. Stamp
`REMAPPED_FROM` INSTEAD of `MOVED_FROM` (`moved_package_complaints` fires on that alone, so
leaving it refuses the fix) and do NOT clear `PIN_VOIDED`; a voided package is not remapped
at all.
**Saving a remapped model takes three steps, and that is the whole cost** -- an earlier
draft of ledger 300 said it could not be saved at all, which was wrong. The arrays are
ndarrays (`to_dict` refuses them, ledger 193) and **the model was already unsaveable before
the move** -- which is why it needed remapping rather than re-declaring. The route, pinned
by a test: `with_grid(<built grid>, remap=True)`, then `with_grid(mf.grid_ref(key))` to
re-point (the `remapped` verdict and the arrays survive, because `REMAPPED_FROM` is tested
before the cellwise branch), then `add_package` the array and `mf.ref` it. A `grid_ref`
target is refused for the remap itself: it needs real polygons.
**NOT built:** collapsing that into one call, `runs.drift()` telling a remap from an edit
(the values move the digest, the stamp does not), hand-typed list BCs (no per-cell values
to weight), and `uzf`/`sfr`/`lak`/`hfb.line` plus the eager array-form `rch`/`evt`
(measured `stale`, not `unverifiable` -- they carry a pin, so re-declaration is their route).

## A per-cell array from a polygon map? `mf.Polygons` (2026-09-23, ledger 299)
`mf.npf(k=mf.Polygons("geology.gpkg", column="k", layer="units", default=1.0))`, and the
same in `k22`/`k33`, `mf.ic(strt=)`, `mf.sto(ss=, sy=)`. A serializable RECIPE resolved
at BUILD against whatever grid the model supplies -- one declaration builds at ncpl 36
and ncpl 100, and `move_report` reads it **`rebuilt`** where a plain ndarray reads
`unverifiable` (ledger 295: an array names cells by POSITION, so there is nothing to
check).
**`GeoPackageSource.k_array` already painted K from polygons and is NOT deprecated** --
do not "fix" it, and do not cite novelty as the reason this exists. What it cannot be is
LAZY or SAVED, and both are measured: it takes a live `ModelContext` so it resolves at
DECLARATION, and it returns an ndarray, which `to_dict` refuses (ledger 193) -- so
`Project.save()` reported **`simulations=0, skipped=1`** and `load()` returned an empty
project, silently, because `save()` reports rather than raises. The recipe form reports
`simulations=1`. For a BC's conductance or elevation this is the wrong tool: the seven
`mf.<pkg>.gpkg` forms already read attribute columns and are lazy since ledger 279.
**THE RULE IS NAMED, NOT INHERITED, AND THE PAIRING IS EASY TO GET WRONG.**
`rule="majority"` (default) / `"mean"` / `"centroid"`. The fork is **`mean` against the
other two**: on a 64-cell grid with one contact it wrote 12 distinct conductivities where
the map holds 2, and its area-integral ran 0.0%-60.5% from majority's across 19 contact
positions -- the arithmetic-vs-harmonic problem in miniature, which is why the default is
the rule that invents nothing. `majority` and `centroid` differ on **at most ONE cell of
64** over that sweep, so a single fixture distinguishes nothing and the test sweeps all
19. **`majority` votes by VALUE, not by feature** (ledger 553): the features carrying one
value are UNITED per cell before the values compete, so a unit drawn as two touching
features (30% + 30%) beats another unit's 40% with or without `default=`, as its
dissolved self does; a tie between values goes to the value of the lowest row among
their features in the cell. A map whose units are each ONE feature paints bit for bit
what the vote by feature painted (the sweep, asserted against that vote written out). An earlier "8 of 64 (12.5%)" figure was centroid vs the area-weighted MEAN, not vs
majority -- that mis-mapping reached two docstrings before being caught. **It is NOT in
ledger 169, and this file claimed it was**: entry 169 contains no such number and the
only occurrence in the tree was that claim. Ledger 299's own sweep is the citable
measurement; a number with no ledger behind it is not an address.
**`cell_shares` is the WRONG helper here.** It normalizes by the whole FEATURE (what
`mf.Spread` needs to divide a total) and drops slivers; filling an array needs share of
the CELL and must keep them. Use `gpd.overlay`, as `cell_coverage` does (ledger 298).
**That reason is right for a polygon map and only HALF right for a grid-to-grid remap**:
when the "feature" IS an old cell polygon, `cell_shares`' per-feature denominator is
exactly the EXTENSIVE remap weight (share of the OLD cell). The reasons to avoid it there
are the sliver rule and a measured 2x cost, not the denominator.
**Refusals are the feature.** A NULL names the feature (ledger 247); a cell no polygon
reaches is REFUSED unless `default=` is given (an unfilled cell is a NaN and MF6 answers
a NaN budget under normal termination); a map landing nowhere prints both extents and
both CRSs (ledger 211). `icelltype`/`iconvert` are refused -- integer flags, and every
rule returns a float, so there is no cell type 0.4; `majority` would be meaningful but
that needs a cast this commit did not measure. The parameter list is an ALLOWLIST, not
"anything that is a `Polygons`", because a helper's `**options` tail also carries
`steady_state`/`transient`/filerecords.
**Two things the array slots accepted and could not use now refuse at DECLARATION**: a
`Surface` or `LayerStack` (accepted, then died inside FloPy with an `MFDataException`
naming neither), and a recipe sitting beside an option of the same name (reached FloPy as
`got multiple values for keyword argument 'k'`). The second was INVERTED:
`_check_input_collisions` EXEMPTS the source shape on the recorded premise that "a future
resolver will consume one into the other", so the shape it called supported got the worse
error. The resolver REPLACES rather than merges, so `PackageSpec.build` refuses it now,
for every package. The `Surface` refusal is deliberately NARROW -- it names the two types
measured to fail rather than allowlisting valid ones, because these slots really do take
a scalar, a sequence, an ndarray, an `OPEN/CLOSE` path string (measured: it builds and
writes `OPEN/CLOSE ... FACTOR 1.0`, at an ABSOLUTE path) and FloPy's own layered dicts.
**NOT built, deliberately: layer ELEVATIONS from a polygon map.** Scoped out by the user
to a follow-up. Nothing does it today -- `Surface.where(zone, inside, outside)` is a
two-way inside/outside choice on centroids, not a per-unit lookup -- but it is a
DIFFERENT route (through `Surface`/`LayerStack`), not another keyword.
**Placement is load-bearing:** the class is in `sources.py` (import layer 0), so the
import ratchet does not move -- `specs.py` already imports that module. The resolver is
`geopackage.array_from_polygons` (layer 3), the rule is
`grid/geometry.polygon_cell_values` (layer 1), and `POLYGON_RULES` is stated in layer 0
with the resolver re-stating the three names (it cannot import upward); a test holds the
two together. **`package_api.py` already had `_lazy_array_factory`** for `mf.rch.array`
at ~line 5890 -- a second definition at line 1000 was silently shadowed by it. Grep
before adding a private helper to that file.
**A RASTER of values is `mf.RasterValues`, and `*` multiplies recipes** (ledger
541). Same slots, plus a per-layer LIST (`k=[mf.RasterValues("k_1.tif"),
5.0, ...]`, one entry per layer; one recipe alone fills every layer) and, AREAL,
`mf.rch.array(recharge=)`/`mf.evt.array(rate=, depth=, surface=)` bare or
`{period: recipe}` (each distinct recipe sampled once per build; `context=` refused
beside one). `recipe * factor` is an `ArrayProduct` of `RasterValues`, `Polygons` and
numbers; `Polygons(value=0.1, default=1.0)` is a constant multiplier (exactly one of
`column`/`value`). It samples through the ONE sampler, `get_raster_vals_at_centroids`
(`band=`/`mean=` added there), and `mean="arithmetic"` reproduces
`Surface.raster(p).values(vor)` BIT FOR BIT -- do not add a second sampler.
`gaps="raise"` by default (`"nearest"` or `default=` fill; `"warn"`/`"ignore"` refused);
a raster on NO cell raises even with `default=`. An `ndarray` factor is refused (it is
what makes a package unsaveable; `with_grid(remap=True)` is that route). **`mf.Raster`
(a `Surface`) stays refused in an array slot** -- it is an elevation with length units.
**`default=` COMPETES for the part of each cell no polygon covers** (ledger
551, the user's call): under `majority` the default is one more VALUE, its vote that
remainder PLUS any polygon carrying the default's own value (ledger 553), and a value
takes a cell only when its cover (united over its features, so a pit split into
touching features takes what the dissolved pit takes, ledger 552) exceeds that vote --
a tie, to 1e-9 of the cell, leaves the default,
so `Polygons(value=0.1, default=1.0)` multiplies the cells it covers MORE than half of --
6 on a 6x6 mesh where it touches 11, which the old rule all multiplied; under `mean` the
remainder is averaged in at the default. Without `default=` nothing changes. A map that
TILES the grid, none of its values equal to the default, paints bit for bit what it
painted before (the 19-contact sweep, asserted: a remainder under 1e-9 of the cell is
rounding and counts as none); the remainder is
measured against the UNION of a cell's pieces, so overlapping polygons are not counted
twice. The recipe did not change, so a stored package's DIGEST did not move while its
resolved array did -- `runs.drift()` cannot see it; a frozen run keeps its old arrays.
`ArrayProduct` is NOT a `_SPEC_TOKEN` object: `specs.py` handles it by type so its
factors are walked as SOURCES (bundler, anchor, `source_paths`, digest) with no wiring.
A recipe lives in `inputs` ONLY: `PackageSpec` refuses one in `options`, so
`packages.fork`/`edit`/`with_options` cannot carry one -- declare the variant through
the helper and `add_package(..., overwrite=True)` (ledger 542).

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
**But the FLOW is still queryable** — `model.hfb.results.q.get()`. A barrier
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
run. **Whether it rejects splits on `icelltype`, which a converted model inherits from
the LPF `LAYTYP`** (ledger 291): a convertible layer is refused by name, a confined one
runs and moves more water than the boundary drawn (measured 2.0x for a drain), with
nothing in MF6's listing file saying so -- a finding with `blocks_run` False is not a
harmless one, and `to_mf6` warns about those at conversion time as well. MF6 stops at
the FIRST package it objects to, which is the argument for reading the list over
running and seeing.
**`fix_for_mf6=True` lifts a GHB head, a DRN invert and a RIV BED back inside their
cells.** A river moves as a PAIR -- `rbot` to the floor, and `stage` up by the SAME
delta **only where the stage was below the bottom too**, so the declared depth
survives while a legal water surface is left exactly where it was drawn. Lifting it
unconditionally (ledger 290's rule, generalised from a case where both were below)
measured a stage 3 ft above land surface and a **reversed sign**, -272.73 -> +1254.55;
`river_pair_lift`/`bottom_clearance` in `modflow/mf6/boundary_support.py` are the one
shared rule and `mf.riv.gpkg` was corrected with it (ledger 291). It does NOT touch a
stage below its own bed: MF6 refuses that at any icelltype and a converter cannot know
which number is wrong. Until 2026-09-18 the clamp lifted RIV's `stage` -- the first
value in the record, so the one a positional read lands on -- and left `rbot`, the
field MF6 checks, exactly where it was.
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
  to native `pyemu.utils.PstFrom`. Targets: `k`, `k33`, `sto.ss`/`sto.sy` (refused where
  they cannot act: an all-steady run, SY on confined cells only), `recharge` (a list RCH's
  records, or an ARRAY RCH's per-period `recharge` files, whose periods come from FloPy's
  records -- files are numbered by block ORDINAL -- and never `irch`), `chd`, `ghb.cond`/
  `ghb.bhead`, `drn.cond`/`drn.elev`, `riv.cond`/`riv.stage`, `wel`, `uzf.vks` (`vks`),
  `lak.bedleak` (`bedleak`), `sfr.rhk` (`rhk`), and transport `porosity` (`mst.porosity`/
  `n`), whose GWT sibling is resolved via `pest/model_lookup.py` — calibrate it with
  heads in the mix, since `v = Ki/n` makes K and porosity near-collinear from
  concentration alone. **A target's files are found by package TYPE** from FloPy's own
  records of the staged model (`native_parameters.target_files`, ledger
  517) -- never a file-name glob, which scaled one of two GHBs in silence --
  and several packages of the type are REFUSED without `package=` (a name, or a
  sequence sharing ONE set; ledger 352's rule -- refused at build when the packages'
  records differ in width, BOUNDNAMES or AUXILIARY on one only, since pyEMU reads one
  set's files as one table; ledger 559). **An ELEVATION (`chd`, `ghb.bhead`,
  `drn.elev`, `riv.stage`, from `FieldSpec.kind`) is an ADDITIVE offset in model units**,
  default bounds (-2, 2), a multiplier on it refused; bounds that exclude the starting
  value (0 for an offset, 1 for a multiplier and every pilot point) are refused
  (ledger 518); `riv.stage` moves RBOT with STAGE (ledger 519).
  An argument that cannot act is refused (`layers=` on a list target -- naming
  `package=`, since `zones=` leaves no list record out (ledger 559) -- or recharge,
  `correlation=` on constant/zone/feature, `pp_spacing`/`pp_points`/`pp_nearest` on any
  style but array pilot points, a style other than the five or one the target does not
  take, an `**extra` key pyEMU's `add_parameters` does not take; ledger 520), and so is a
  `name=` another call holds (523). **`periods=` groups a target's per-period files**
  (ledger 549; it replaced `temporal=`, which nothing read): `'month'` is one
  parameter set per calendar month of the TRANSIENT periods, `m01`..`m12` by each
  period's MIDPOINT on the run's calendar (`period_calendar`, ledger 529), shared by
  every cycle; a STEADY period is in no seasonal group (the user's rule, Q3 -- stack a
  plain call for it); `{label: [periods]}` gives explicit groups. A period repeating
  another group's block (MF6's rule for a period with no block) is REFUSED, naming the
  fix: a block per period (`mf.AverageYear(...).by_month(...)`). **`lak.bedleak` and
  `sfr.rhk` are FEATURE targets** (ledger 558): a row is keyed by its lake
  and connection / reach NUMBER, never its cell (a lake meets one cell twice, two reaches
  share one), with `use_rows` computed from the
  staged package and the template checked to scale exactly those rows (pyEMU falls back to
  EVERY row on a miss). `style='feature'` is one multiplier per lake (bedleak's default)
  or per named stream (a reach BOUNDNAME); `features=` picks lakes/reaches/streams by
  ZERO-based number or boundname, slugged by the observation groups' ONE rule
  (`pest/naming.py`'s `pest_group_name`); the boundname match is
  `package_registry.feature_boundnames`, shared with the observation targets. Refused: a
  `NONE` bedleak on a lake it scales, ANY reach with cellid `NONE` (one column short --
  pyEMU multiplied its Manning's n), and `correlation=` on an rhk `grid` whose reaches share
  a column (pyEMU then decorrelates every reach). Styles: `constant`, `zone`, `grid` (one
  geostat-correlated multiplier per Voronoi cell + correlated prior), and `pilotpoints`
  (a `pp_spacing=` net in MODEL LENGTH UNITS or explicit `pp_points` — `pilot_points.py`;
  pyEMU's own pilot points are unusable on unstructured grids). **Pilot points MULTIPLY
  onto what pyEMU wrote** (ledger 525): ONE forward-run command per array
  file -- per group of files covered by the same sets, so every period file of an
  array RCH shares one interpolation (ledger 547) -- after `apply_list_and_array_pars`, so they stack with constant/zone/grid and
  with each other -- they used to rebuild the layer from the build-time base and erase
  every other multiplier on it. Each cell reads the `pp_nearest` (8) nearest points of
  its SET -- one call, one layer, one `zones=` label -- by IDW in log10 space, through
  weights stored at build (526); `zones=` is honoured (527), and
  the template's `pilot_point_sets.json` lists every set. The legacy `add_parameter` +
  `build_pst` + `KPilotPointParameter`/`DrainElevation`/`DrainConductance` specs were
  RETIRED (deleted 2026-06-22) — do not reintroduce them.
- `PestProject` (`project.py`) — orchestrates pyEMU/PstFrom; native `build()` + `run_ies`/
  `prior` (`workers=` runs parallel PESTPP-IES agents). **Construct it with
  `model.pest(name, start_datetime=...)`** (the front door — defaults the workspace to
  `<model workspace>.pest/<name>` — a SIBLING of the model dir, ledger 107), not by
  importing `PestProject` directly.
- **A build only READS the run and never touches the model object** (ledger 356, the
  user's call). `build()` copies the run's INPUT files to `<template>.staging` -- what
  `mfsim.nam` leads MF6 to read (`project._simulation_input_files`: name-file entries,
  `FILEIN` scanned in turn, `OPEN/CLOSE`), never its outputs, run metadata or an old
  PEST folder, which every template and IES worker carried (92 of the canonical full
  template's 107 MB; ledger 383) -- loads the
  simulation from THAT (`MFSimulation.load` + `keep_implicit` -- FloPy's loader drops LAK
  `IMPLICIT`, ledger 326), relayers (ledger 110) and externalizes it there, hands it to
  `PstFrom`, and removes it, failure included. It used to do all of that to the CALLER's
  model and write it into the run directory (every input rewritten, 22 files added on
  the test toy; 18 storage entries of the object changed). Everything the build reads off
  a FloPy model goes through ONE helper, `native_parameters.staged_model` -- the recipe
  arrays, the pilot points' cell centres (their base is the template's own array file,
  ledger 525), the default spatial reference;
  do not read `project.model.gwf` there. **The run's FILES win over the object**: the base
  observations come from its OUTPUT files, so an in-memory edit after the run is not in
  the template, and a model never written is REFUSED (no `mfsim.nam`). Two old-build
  crashes went with it: a lake with empty period blocks (FloPy died externalizing the
  built object) and an array target on a run opened with `mf.load_mf6_run` whose arrays
  were inline (loaded lazily: NPF never loaded, so nothing was externalized and the build
  found no file for `k`; with pilot points, no base array). pyEMU opens `PstFrom.log` in the CWD; it is built with
  the staging dir as cwd and the log moved into the template. **Resolve every path BEFORE
  that `chdir`** (ledger 357): the template's was resolved inside it, so a RELATIVE model
  workspace -- every calibrating notebook's `Path('../artifacts/...')` -- built its
  template inside the staging copy and died; tests under `tmp_path` (absolute) never saw
  it. A rebuild closes the previous build's log before deleting the template. Cost: canonical full
  (100x100x4) 31 s -> 9 s -- externalizing a LOADED simulation is far cheaper than a built
  one (23.9 s vs 2.0 s on the validation profile); the copy itself 0.03-0.09 s on tmpfs.
- **A build deletes its workspace ONLY when a build wrote it** (ledger 360).
  `model.pest(name, workspace=W)` makes `W` the template, and `build()` used to `rmtree` it
  WHATEVER it held, and `PstFrom(remove_existing=True)` again: `W` = a Project root
  deleted the model's own run. The template and `<W>.staging` are now deleted only when
  absent, empty, or carrying the `myflopy_pest_build.json` the build writes into each,
  naming the directory's RESOLVED path -- so a COPY (pyEMU's IES master and workers, a
  backup) is not taken for one; a template built before 360 counts only when its metadata
  names this calibration and this run and it IS where `model.pest(name)` puts it (ledger
  361: the NAME alone was compared, so a `cp -r` backup was deleted). **The marker is
  VISIBLE on purpose** (361): as a dot-file it survived `rm W/*`, and a folder emptied that
  way and refilled was deleted by the next build. Anything else is REFUSED
  (`FileExistsError`), checked before anything is deleted, and `PstFrom` gets
  `remove_existing=False`. A `W` that is, contains or sits inside the run directory is
  refused (`ValueError`) at `model.pest(...)` and at `build()` whatever markers it holds --
  ledger 107's allowance of a template inside the run is gone. **The parallel RUN is
  guarded the same way** (361, amended by ledger 532): a launch deletes `master_dir`
  (default `<template parent>/<name>_ies_master`) and the workers it uses, so
  `run_ies`/`prior(workers=N)` first requires each to be absent, empty, or a COPY of THIS
  template -- the marker, copied, names the template. **Each master has its own workers,
  `<master>.workers/worker_<n>`** (they shared `<model>.pest/worker_<n>` until then, so a
  second calibration's run deleted the first's), **a master holding results is refused
  unless `overwrite=True`**, and **a directory a run is still going in is refused whatever
  `overwrite` says**: every launch writes `myflopy_pest_run.lock` (host, pid, process start)
  into the template (serial) or master and worker directory (parallel); a live pid, another
  host, or a PEST++ process working there (`/proc` cwd) is a live run, a dead pid a STALE
  lock (warned, ignored). The launch copies the master itself (without the template's
  serial outputs) and hands pyEMU `reuse_master=True, cleanup=False`; workers go after a
  CLEAN finish unless `keep_workers=True`. **`cal.resume_ies(workers=N)` continues an
  interrupted history match** (ledger 533): PESTPP-IES refuses `/r`, so it restarts
  from the run's own ensembles and the outputs are spliced back into ONE run of
  iterations `0..N` (`pest/run_state.py`).
- **Observations** accepted directly by `cal.observe(...)`/`cal.forecast(...)`: `HeadTargets`,
  `LakeStageTargets`, `SfrStageTargets`, `SfrFlowTargets`, `DrnFlowTargets`, `RivFlowTargets`
  (or their pre-built `*ObservationSpec`). All are wired into the native build + forward run.
- **A screened well is ONE head observation** (ledger 556): `screen_top`/
  `screen_bottom` on a `HeadTargets` location (on when the columns are present) reads the
  ACTIVE layers the screen opens into, weighted by NPF `K` x saturated open length --
  `hydraulics.screened_head`, sharing `saturated_thickness` with `transmissivity()` -- and
  the deepest open layer's head when all are dry (under NEWTON a finite head; Ten Trails
  holds ~17% of active cells below their bottom). The forward run weights with the
  REALIZATION's K files (`k_sources` in the head config, the per-layer `npf_k_layer*.txt`
  the build re-stores and pyEMU/pilot points rewrite). A set with only `layer` columns keeps
  ledger 234 (two layers = two full-weight observations). `compare()`/`simulated_heads()` and
  the forward-run writers read through ONE period-end reader of the OC-declared file,
  `target_values.read_target_values` (ledger 557) -- do not read a target's
  heads through `HeadsPlus`/`all_heads` again (7 s and 1.9 GB per `compare()` on Ten Trails).
- **`RivFlowTargets` (creek baseflow) and `DrnFlowTargets` are ONE zone implementation**
  (ledger 534): `observations/_zones.py`, each family a thin subclass naming its
  MF6 quantity and sign rule -- add a zone family there, never copy `drn.py`. RIV is MF6's
  sign with NO sign refused (baseflow `-Q`, a losing reach positive); DRN keeps its `-Q`
  refusal. The budget route reads by `paknam2` through `build_budget_result_table`:
  `package=` reads one package, `None` sums every package of the type in the zone's cells,
  and a package in scope that wrote no record (`save_flows` off) is REFUSED (ledger
  536). A period with no boundary in a zone reads 0.0 (537); a zone with none in ANY
  saved period is REFUSED naming its layer(s) and `package=` (537 as amended: wrong
  layer or package, not zero flow). A test
  double reaching that route needs a budget READER, not `model.bud`:
  `tests.test_drn_zones_on_faces.attach_budget`.
- **Dated readings become period keys in ONE place** (ledger 543):
  `HeadTargets.from_hydrographs` (and `ConcTargets`) / `LakeStageTargets.from_hydrographs`
  bin a WIDE or LONG date-keyed table against `model.calendar` through
  `observations/hydrographs.py` -- history: the reading NEAREST each period's END (ledger
  336; on the Ten Trails TDIS, start 2017-09-30, a period ends 00:00 on its month's last
  day, so containment would move a last-day reading a month on); average year: per
  calendar month the median over years of each year's monthly median, on the LAST
  `mf.AverageYear` cycle (a model that is not one is refused). Binning is EAGER and the set
  carries the model's calendar as a CLOCK STAMP: `compare()` and the PEST build refuse it
  on a model with other dates (another cycle count). Do not bin dates anywhere else. The
  build records each set's `period_end_dates`, and the IES review draws dates (ledger
  544). **How two calendars differ is ONE rule, `PeriodCalendar.clock_difference`**, and
  a run's cycle ONE reading, `period_calendar.monthly_cycle` (ledger 546): the stamp check
  and `AverageYear.equilibrium` both ask the first and word it themselves; never compare
  two calendars by hand.
- **What PEST is handed is exactly the targets** (ledger 340, 2026-09-27). The build
  matches each target to the row pyEMU's INSTRUCTION FILE reads -- by POSITION, the key
  read as a period by NUMBER through `_shared._numeric_periods` (`1.0` is period 1,
  `1.5` refused, in `compare()` too) -- and never rebuilds pyEMU's observation NAME from
  a key's spelling: that is what silently dropped `per:1.0` against pyEMU's `per:1`, and
  `time:10` against `time:10.0`. A target nothing holds (a head/conc target on no
  active cell too, ledger 341 -- `compare()` only warns), two values for one (target,
  period), a value that is not a finite number and two sets under one prefix are
  REFUSED by name; a blank is unmeasured. The forward-run writers and the IES review read target
  names as TEXT, through `forward_run.read_target_table`: pandas read a well named `007`
  as `7` and `NA` as missing (ledger 341). Every row no target measures carries weight 0 (a forecast set's
  other periods too); pyEMU's default is weight 1 and the BASE value, which anchored the
  calibration to the starting model. Heads/conc base and forward CSVs come from ONE
  writer (the pivot SORTED its columns: PEST read W1's head as W10's in every
  realization); the named-series forward run reloads the calibrated model BY NAME.
  **A template built before this bakes the old forward run and must be rebuilt.**
- **A target's `group` label is its PEST observation group** (ledger 358, the user's
  call). `assign_observation_groups` (in `pest/observations.py`, run by
  `finalize_observations`) puts a MEASURED observation in its label made a PEST name
  (`pest_group_name`: lower-cased, every run of characters but `[a-z0-9_.:-]` one `_`),
  an unlabelled one in its set's PREFIX, every row no target measures in
  `<prefix>_unmeasured` and a forecast in its prefix -- so a user's group holds exactly
  what was measured under it. pyEMU named every column a group of its own, so a group was
  one well. REFUSED by name: a label that is not text (a number -- deliberately not "no
  label"; the advice is `astype("string")`, which keeps a blank blank), the TEXT
  `nan`/`None`/`<NA>` (ledger 359: `astype(str)`, the first advice, wrote a blank so and
  made a group `nan` of every unlabelled well), no letter or digit, a set's group name, two labels becoming one name
  (`North`/`north`), and a label OR PREFIX starting `l_`/`less`/`g_`/`greater`, which
  PESTPP-IES reads as an INEQUALITY on the bare start (`lessard`, `greaterville`;
  measured: phi 0 where the same data scored 424.864; PESTPP-GLM at NOPTMAX=0 does not).
  That refusal is of a USER's label: **a bound is `inequality="at_most"|"at_least"` on
  the target SET** (ledger 564, every family): the build alone writes
  `l_`/`g_` before each measured observation's ordinary group, the user-label refusal
  above stays and names it, a flow bound is in MF6's sign (`-2.0` + `"at_most"` = at
  least 2 leaving), `compare()` gives PEST++'s residual plus `raw_residual`, and the IES
  review scores one-sided through ONE helper on pyEMU's tag rule (`phi_contributions`
  reads PESTPP-IES's `phi.group.csv`). **Difference observations are an OBJECT,
  `targets.differences(kind="mean"|"previous")`** (ledger 565): each
  label minus its mean over ITS measured periods (source weight) or its previous
  measured period (`w/sqrt(2)`), built in myflopy's forward run from the ONE period-end
  reader (`target_values.target_differences` is the one rule), ONE group per set
  (`hdiff`/`cdiff`), undefined rows 0.0 at weight 0 -- never pyEMU's
  `setup_temporal_diff_obs`, which differences by group.
  No length limit (PEST++ 5.2.16 read 5000 characters); a comma or `"` breaks the
  version-2 file. Lake/SFR targets have no label column -- split them by prefix. **No
  default prefix is inside another's groups** (a phi-factor tag matches anywhere in a
  group name): the lake default is `lak_stage`, was `stage` -- inside `sfr_stage` --
  until ledger 372; `test_no_default_prefix_is_inside_another_sets_group`. **The
  review draws one time series per target LOCATION** (`IesResults._location_panels`, a
  column found by pyEMU's `<base>_<index>:<value>` name), never per group, and marks only
  MEASURED observations -- per group it drew a group of wells as one zig-zag line, and
  the old per-column review drew every zero-weight row as a red "measured" triangle at
  the base run's value. A run built before keeps one group per column and still reviews.
  `phi_contributions`/`plot_conflict` read the user's groups unchanged.
- **A DRN zone's cell is a `(layer, cell)` pair, on every route** (ledger 342). A
  cell-budget record names a drain by MF6's model-wide USER node (`layer * ncpl + cell`
  on a DISV, past an inactive cell too -- pinned by `tests/test_drn_zone_layers.py`);
  the budget route matched it against a zone's IN-LAYER numbers, so a zone in layer 1
  read the layer-0 drain sharing a number (-6.77 for drains holding -5975.73), and the
  PEST zone snapshot dropped `layer`. It matches the pair through `model.node_to_lni`
  now -- a test double reading the DRN budget must carry that map -- and the snapshot
  carries `layer`. A zone is ONE zone of its rows' pairs, compared once and written to
  MF6 once per cell (`<name>_lay<layer>_c<cell>` when it spans layers); a row with no
  `cells` entry is resolved from its geometry, so one table may mix the two, and an
  explicit pair the grid does not have is refused on the budget route (ledger 343: it
  matched nothing, so `[20, 99]` read cell 20 alone). A row with no `layer` reads EVERY
  layer of its columns, those a GEOMETRY covers (ledger 535, 482's area rule) and those a
  `cells` list names (ledger 550, the user's call) alike -- both read
  layer 0 alone, 1.3% of the Ten Trails `ss` south drains; `to_flopy_obs(None)` refuses
  such a zone, having no layers to read. A partial
  miss of a DRN zone or a head/conc target is a Python warning as well as a log line
  (ledger 288's rule).
- **The observation route reads what the package DECLARES, by identity** (ledger 344).
  Lake/SFR/DRN targets took the first `*.csv` in the workspace, in sorted order, whose
  header held their NAMES -- a GWF head, a LAK `VOLUME`, another lake's stage, a
  same-named SFR flow, or a STALE copy (a backup; the CSV a reused run directory kept),
  which froze every PEST realization. `_shared._declared_route` asks FloPy for the
  owning package's OBS records (`package.obs[i].continuous`): a label counts only when
  declared as the family's quantity for the target's OWN lake/reach/cell (a number from
  one, or a boundname via `packagedata`), and `_read_declared_observations` reads only
  those files, each record BY POSITION, as text or MF6's `BINARY` form (ledger 345; MF6
  cuts names to 40 characters; the header is checked, not matched). Refused by name:
  another quantity only, another id, a label declared twice, a set declared in PART
  (the route is chosen for the whole set, since the two routes did not read the same
  numbers -- a dry lake until ledger 351, an SFR flow's sign until ledger 350, the LAK
  stage file's record clock until ledger 348), a
  missing declared file, and a model with several LAK/SFR/DRN packages unless
  the targets name one (`package=`, carried into the forward run; packages found by
  TYPE -- FloPy's `get_package("lak")` answers a package NAMED `lak` first; the stage
  VIEWS read that one alone in silence until ledger 352). `to_flopy_obs` refuses
  a `kind=` other than the quantity and two names one label once upper-cased. A test
  double reaching the observation route must answer `model.package(type)` with a
  package whose `obs[i].continuous.get_data()` holds the records
  (`tests/test_mf6_pest.py::_declare_observations`); one without `package` declares
  nothing and reads the file route. A template built before this reads by
  identity without a rebuild -- its baked forward run calls the library's
  `simulated_series` -- but cannot pass `package=`, so a model of several packages of
  the type must be rebuilt.
- `forward_run.py` — injected forward run: pyEMU's `apply_list_and_array_pars` for params,
  plus post-processors that regenerate head + named-series (lake/SFR/DRN/RIV) simulated CSVs.
  Every named-series set goes through ONE `write_named_series_targets(sets=[...])` call that
  opens the run once (ledger 540: one full reload PER SET cost ~7 s each beside
  a 17-20 s Ten Trails `ss` run); do not add a per-set call back.
- **Run discovery + review**: `model.pest_runs` / `run.pest_runs` (`runs.py`,
  `find_pest_runs`/`PestRunHandle`) list calibrations done on a model; `.review()` reopens
  one as `IesResults` (`ies.py`, `open_ies_run`) — phi, ensemble-vs-obs, forecasts,
  parameter-field maps. (The old `results.py` deterministic-review layer was deleted.)
- **Write-back: `ies.calibrated(realization="base").to_library(project)`** (ledger
  560). It REPLAYS the realization's own forward run -- templates filled at
  full precision (pyEMU's filler keeps 4 significant digits), `forward_run.py` cut by AST
  before `pyemu.os_utils.run`, a subprocess with the absolute `PYTHONPATH` -- and reads
  back each dataset the build metadata's `parameters` records (version 2: `model`,
  `variable`, `fields`). Never re-apply a multiplier in Python. `base` reads
  `<case>.<it>.base.par` (16 digits), any other label the par ensemble (6). The base of
  the substitution is the RUN's frozen snapshot and library copies, never the library now;
  values-only, so a list package keeps its pin after every cell is checked and drops its
  `.gpkg` move recipe (a boundname is compared without regard to case: FloPy LOADS it
  lower-cased); forks to `f"{label}/{pest}"` keys plus a forked simulation, never
  edits -- the fork is refused when the library's other packages moved since the run
  (ledger 562). A fork/edit/`with_options` override over a recipe-filled option REPLACES the
  recipe (ledger 561). LAK `bedleak` and SFR `rhk` come back too, and FloPy reads their
  columns as TEXT (`bedleak` may say `NONE`, `man` may name a time series): a spelled
  number is a number, never identity, and a keyword the forward run left stays (ledger 563).
  **A review reads ONE run** (ledgers 507, 508): a serial launch deletes
  an earlier run's `<name>.*` outputs from the template (`runs.pestpp_ies_outputs`, the
  one list) and says so (`UserWarning`, ledger 532), a parallel master is copied
  from the template without them, and `IesResults.iterations` counts only what the run's phi file records. A handle's
  `executions` are the masters named EXACTLY `<name>_ies_master`/`<name>_prior_master`
  plus a serial run in the template (kind from its `<name>.pst`'s NOPTMAX), its case is
  `<name>` (what PESTPP-IES writes, whatever `build(filename)` wrote), and
  `review(kind)` REFUSES a kind it lacks -- it fell back silently to the prior.
  `field()` takes a log-transformed field's statistics in log10 (geometric mean), as
  `sensitivity()` does (ledger 514). `run_ies`/`prior` default
  `overdue_giveup_fac=5.0` (ledger 510) and build with `chunk_len=10**6`
  (ledger 511).
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
