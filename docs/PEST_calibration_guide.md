# A practical guide to PEST++ calibration & uncertainty with myflopy

A plain-language summary of the methods taught in the GWM² 2026 PEST workshop
(Jeremy White, Mike Fienen, Katie Markovich), mapped to the myflopy calibration
API. The goal is to give you the *why* behind each step and the one or two lines
of myflopy that do it.

The whole thing is one idea from Bayes' rule:

> **posterior  ∝  prior  ×  likelihood**
>
> What we believe *after* looking at the data = what we believed *before* (the
> **prior**), reshaped by how well the model can match the data (the
> **likelihood**). "Learning" is just reducing uncertainty about predictions we
> care about.

Everything below is a piece of that sentence.

---

## 0. Purpose first (and keep it simple)

Before any calibration: know your purpose. Are you *explaining* a system or
*predicting* an unseen outcome? Applied work is almost always prediction, which
means the deliverable is a **forecast with honest uncertainty**, not a single
"calibrated" model. Build the *minimum feasible complexity* — a model fast
enough that you can run it hundreds of times, because you will.

A forecast you can defend > a pretty model you ran once.

**Time setup for an average-year calibration.** Declare the steady period and the
monthly cycles with `mf.AverageYear` (ledger 530) -- its `tdis()`, STO flags
and `by_month()` forcing -- and before building PEST, run it once and read
`year.equilibrium(run.model(name), tolerance=...)`: every cycle is twelve more
periods in EVERY realization's forward run, so take the fewest cycles whose last two
agree to within your target noise. Targets are compared on the last cycle
(`year.months_of(-1)` maps each calendar month to its period). A run recovers its
cycle from its own files with `mf.AverageYear.from_model(model)`, and
`model.calendar` dates every period; `model.pest(name)` takes its `start_datetime`
from the same TDIS `START_DATE_TIME`.

---

## 1. The prior — what you're uncertain about (parameterization)

The prior is your expert knowledge expressed as adjustable parameters and their
plausible ranges. Use **lots** of parameters (heterogeneity is real), but
express their spatial/temporal correlation so they stay physically sensible.

Key parameters: hydraulic conductivity (K), recharge, storage, riverbed/GHB/drain
conductances, pumping. Two shapes:
- **array-like** (K, recharge, storage) — one value per cell, usually controlled
  by pilot points or multipliers at several scales;
- **list-like** (pumping, boundary heads, conductances) — one value per feature.

In myflopy you *declare* what to adjust; it compiles to native pyEMU:

```python
cal = model.pest("calib", start_datetime="2020-01-01")   # the front door; workspace
                                                          # defaults to <model ws>.pest/calib

# bounds = how far the multiplier can move; physical = hard limits on the final value
cal.parameterize("k",        style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
cal.parameterize("recharge", style="grid",     bounds=(0.5, 1.5), physical=(0, 1e-2))
cal.parameterize("ghb.cond", bounds=(0.1, 10))

# transport properties live on the GWT sibling; the target finds it for you
cal.parameterize("porosity", style="constant", bounds=(0.5, 2), physical=(0.05, 0.5))
cal.parameterize("uzf.vks", style="grid", correlation=400, bounds=(0.2, 5))

# zones: one multiplier per label, from a raster, a polygon layer, or an array
# (a second "k" parameter needs its own name: a repeated one is refused, ledger 523)
cal.parameterize("k", style="zone", zones=mf.ZoneSpec.from_raster("hsu.tif", model),
                 name="kz")

# elevations are OFFSETS in model units; bounds default to (-2, 2)
cal.parameterize("drn.elev")
# two RIV packages (riv + crisp): name the one, or share one set
cal.parameterize("riv.cond", package="crisp", bounds=(0.1, 10))
cal.parameterize("riv.stage", package=("riv", "crisp"), bounds=(-1, 1))
```

**A target scales the packages of its TYPE.** `ghb.cond` scales every GHB record
whatever the package is named, through the files FloPy's own records of the run's
packages name -- it used to guess the file name from the type, so a GHB named `edge`
was not found and a second DRN (`seeps` beside `drn`) was left out with no word (ledger
517). A model with SEVERAL packages of the type refuses a call with no
`package=`, naming each and the spellings: `package="crisp"` picks one, and
`package=("riv", "crisp")` shares one parameter set (its default name gains
`_<package>`, `rivcond_crisp`). Two creeks are usually two parameters.

**An elevation is an offset, never a multiplier.** `chd`, `ghb.bhead`, `drn.elev` and
`riv.stage` are additive in model length units -- `bounds=(-2, 2)` by default --
because a multiplier on a head scales with the datum (x1.1 moves a 400 ft head 40 ft and
a 40 ft one 4). `additive=False` there is refused, and so are bounds that EXCLUDE the
value a parameter starts at: 0 for an offset, 1.0 for a multiplier and every pilot point
(`bounds=(2, 5)` would start each one outside its own range, ledger 518).
To raise K by 2-5x, change the model and calibrate around it. **`riv.stage` moves each
record's river bottom by the same offset**, so its depth -- what sets the river's flux
-- survives every realization (ledgers 290/291); a record already below its bottom is
refused at build.

**An argument that cannot act is refused, not ignored** (ledger 520):
`layers=` on a list target (use `zones=`), `correlation=` on `constant`/`zone`, a
variogram shape (`anisotropy`/`bearing`/`nugget`) without `correlation=`,
`temporal=`, and any style but `constant`/`zone`/`grid`/`pilotpoints` -- pyEMU's own
`"pp"` names `"pilotpoints"` -- and `"pilotpoints"` on a list target (pilot points
interpolate onto an array; use `"grid"` with `correlation=`). A parameter name used twice is refused too: give each
stacked call its own `name=`.

**Zones** are per-cell integer labels; `style="zone"` gives one adjustable multiplier
per distinct label. `ZoneSpec.from_raster` samples a categorical raster by **majority
vote**, not by the area-weighted mean this library uses for surfaces — averaging zone
ids is meaningless, and a cell straddling zones 1 and 3 would come out as zone 2, which
it does not touch and may not exist. A cell smaller than a pixel takes the label under
its centroid; only a cell off the raster, or on its no-data, gets `nodata_zone`
(ledger 213). `from_polygons(gdf, column="hsu")` labels cells by
the polygon they fall in — by the cell's **centroid**, and a zone is closed, so a
centroid exactly on a zone's edge is in it and a boundary drawn along a row of
centroids does not leave that row unparameterized (ledger 243); a cell covered by
two zones takes the first row, a polygon whose label is NULL is not a zone at all,
and a layer that reaches no cell says so. A plain array still works.

You never have to think about SHAPE. pyEMU wants a per-cell array for array targets
(`k`, `k33`, `porosity`) and a `(layer, cell)` array for list targets (`recharge`,
`chd`, `ghb.*`, `drn.*`, `wel`, `uzf.vks`), and rejects the wrong one with an error
naming neither the target nor the shape. One `ZoneSpec` serves both.

**`uzf.vks`** is UZF's saturated vertical K, and the only UZF quantity worth
calibrating on its own: measured on the canonical model, `vks x3` moves heads by
1.026 ft while `finf x2` moves them 0.011 ft and `surfdep x5` moves them 0.0001 ft.
`thts` and `eps` are left out deliberately — they are 0.995-collinear with each other,
so a set containing both is ill-posed in the same way K and porosity are.

UZF is also the one list package whose external file does not open with the cell
identity (it numbers its own cells: `iuzno layer icell2d ...`), which the recipe
declares via `index_cols`. You never pass that, but it is why a UZF parameter lands on
its own cell rather than all 189 of them sharing one set of coordinates.

**Porosity needs concentration data — and heads.** It does not appear in the flow
equation at all, so heads alone cannot constrain it (measured on the canonical model:
a 3x porosity change moves heads by *exactly* 0.000000 ft). But transport velocity is
`v = Ki/n`, so raising K and lowering porosity move a plume almost identically —
the concentration responses to `K x3` and `porosity /3` have cosine similarity **0.98**.
Estimating both from concentration alone is therefore ill-posed. Observe **heads and
concentration together**: heads pin K, and concentration then pins porosity.
`build_canonical_transport_calibration_demo` ships both target families for this reason.

`layers=` is rejected on a target MODFLOW writes as one whole-grid array, rather than
being quietly ignored the way it used to be. A multi-layer `style="grid"` parameter is
one group per layer, named by the ZERO-BASED layer (`kl0`, `kl1`, as pilot points are);
a correlated one is refused above 5,000 parameters per group, because its prior is a
dense covariance of every pair (about 1 GB at 5,000, 12 GB for one 20,323-cell layer;
ledger 522) -- pilot points carry a correlated field for a fraction of that.

**`style="pilotpoints"` works on porosity too**, as of 2026-07-30. It was briefly
documented here as flow-only, because the pilot-point interpolation read a hardcoded
`gwf.npf.k` as its base array — which was not a transport limitation at all but a bug
that also made `parameterize("k33", style="pilotpoints")` overwrite K33 with horizontal
K, 30x wrong and silent. Pilot points now multiply the array their own target names, so
every array target (`k`, `k33`, `porosity`) supports them. See ledger 115.

### Pilot points (array targets)

```python
# a net 1,500 ft apart on layers 0 and 1 (pp_spacing is in MODEL LENGTH UNITS)
cal.parameterize("k", style="pilotpoints", pp_spacing=1500.0, layers=[0, 1],
                 bounds=(0.1, 10), physical=(1e-3, 300), correlation=4000.0)
# a broad multiplier under them: the two MULTIPLY
cal.parameterize("k", style="constant", layers=[0, 1], name="kc", bounds=(0.2, 5))
print(cal.settings())    # after build(): the point count of every set
```

- **Placement.** `pp_spacing=` lays a regular net -- squares that wide from the lowest
  cell centre, the cell centre nearest each occupied square's middle becoming a point
  -- or `pp_points=` gives the locations (a GeoDataFrame of Points or `(x, y)` pairs in
  the grid's frame). Give exactly one. `pp_spacing` replaced `pp_space`, which counted
  in median cell widths and laid 4,501 points per layer on Ten Trails' refined mesh
  (ledger 527); 1,500 ft lays about 100 there. `build()` records each set's
  point count on the spec (`spec.pilot_point_sets`) and in `cal.settings()`.
- **Interpolation.** Each cell takes the `pp_nearest` (default 8) nearest points of its
  set, weighted by inverse distance squared, in **log10 space** for a log-transformed
  parameter -- a cell between a point at 10 and one at 0.1 reads 1, not 5.05. The
  weights are computed once at build and stored in the template
  (`<set>_pp.factors.npz`), so a forward run does no distance arithmetic. The old rule
  (every point of the layer, linear) let a block of points set to x5 reach only 2.3-4.1
  between them and biased a lognormal field 13-60% high
  (`docs/reproducers/pilot_point_interpolation.py`, ledger 526).
- **Stacking.** Every pilot-point set on an array file is applied by one forward-run
  command after pyEMU's multipliers, multiplying the file as pyEMU wrote it: a constant,
  zone or grid multiplier and any number of pilot-point sets on one layer all reach
  MF6 as a product. Until ledger 525 a pilot-point set rebuilt the layer from
  the build-time base and erased them all (`docs/reproducers/pilot_points_stack_and_zone.py`).
  A covered cell is clamped to the tightest `physical=` of the sets covering it and of
  the pyEMU multipliers on its file, so a constant's `physical=` still holds under pilot
  points (ledger 528).
- **Zones.** `zones=` gives each label `> 0` its own set on each layer: its own net (or
  the explicit points inside it), interpolated onto its own cells only, with its own
  prior covariance; cells labelled `<= 0` keep their value. Different spacing, bounds
  or variogram per zone are separate calls, each with a masked `zones=` and a `name=`:

  ```python
  cut = mf.ZoneSpec.from_polygons("cut.gpkg", model)          # 1 inside, 0 outside
  outside = mf.ZoneSpec.from_array(1 - cut.values)
  cal.parameterize("k", style="pilotpoints", layers=[2], zones=outside, name="k3",
                   pp_spacing=1500.0, bounds=(0.1, 10), correlation=4000.0)
  cal.parameterize("k", style="pilotpoints", layers=[2], zones=cut, name="k3cut",
                   pp_spacing=500.0, bounds=(0.01, 100), physical=(1e-4, 10))
  ```

  No point outside the polygon moves a cell inside it, so the cut's contrast in the
  base field survives at multiplier 1.
- **Names.** A set is `<name>l<layer>` -- `<name>l<layer>z<zone>` with zones -- and its
  points `...pp0000`, lower-cased as PEST reads them. Every call needs its own name --
  `parameterize()` refuses a repeated one (any case, ledger 523), on one layer or two --
  so pilot points declared one layer per call take a name each. The build metadata's
  `parameters` record names the same sets as `pilot_point_sets.json`.
- **Refused at build**: an explicit point in no cell, in a cell labelled `<= 0`, or given
  twice; a zone `> 0` no explicit point lies in; zones labelling nothing; `physical=`
  limits on one cell -- of its pilot-point sets, or of them and a multiplier on the same
  file -- that share no value.
- **What a write-back reads.** `pilot_point_sets.json` in the template lists every set
  (layer, zone, its value, points and factor files, interpolation, limits) and, per
  array file, whether the forward run starts from pyEMU's output or from the captured
  `<file>.pp_base`.

**An array target needs its array in the run's files.** `parameterize("k33")` on a
model whose NPF was written with no `k33=` has nothing to multiply -- MF6 then takes
K33 from K -- and `build()` refuses it naming the target and the package. It used to
die inside the build with `cannot reshape array of size 1 into shape (2,newaxis)`.
Declare `k33=` on the NPF and run the model again, or calibrate `k` alone.

`physical=` is the safety rail: no matter what calibration does, the *final* K
stays inside those bounds. Always set it for multipliers.

**`correlation=` on pilot points is read as of 2026-07-30.** It was documented from
the start and had no effect: pilot VALUES were drawn independently from their bounds,
and the reasoning in the code was that IDW interpolation supplies the spatial
structure. It supplies smoothness, but not the correlation LENGTH you asked for — the
field's structure came from pilot spacing instead, so a 200 m and a 2000 m variogram
gave the same prior. Measured after the fix on the canonical model (nearest pilot pair
~372 m): a 150 m range leaves neighbours independent (+0.01), a 3000 m range correlates
them (+0.85) and decays with distance. It matters most where data are thin — IES then
moves neighbouring points together instead of letting each wander alone.

**The prior `run_ies()`/`prior()` hand PESTPP-IES holds every parameter** (ledger
385): a draw of the whole control file -- each correlated `grid` parameter through
its variogram, every other one on its own from its bounds -- with each correlated
pilot-point set's values replaced by its own correlated draw. Before, a correlated
`grid` parameter beside correlated pilot points, or correlated pilot points beside any
other parameter, stopped both calls before a model ran. `cal.draw_prior()` draws the
same way by default (`use_specsim=False`, ledger 384); pyEMU's spectral simulation
(`use_specsim=True`) is possible only on a regular structured grid and is refused on a
DISV/Voronoi one. PESTPP-IES holds the initial realizations to each parameter's bounds,
so a drawn value beyond a bound runs at the bound. **Every prior draw spreads as far as
the bounds say** -- six standard deviations between them (`PRIOR_SIGMA_RANGE`), pyEMU's
own rule -- pilot points included since ledger 424: their correlated draw took the
variogram's sill as a log10 variance of 1 whatever the bounds, so with bounds `(0.5, 2)`
77% of it lay outside them and PESTPP-IES ran those at a bound
(`docs/reproducers/pilot_point_prior_bounds.py`). Widen the bounds to widen the prior.

**On Tikhonov regularization: myflopy does not offer it, deliberately.** It is a
PESTPP-GLM concept, and the runner here is PESTPP-IES, which prints
`prior information equations not supported in ensemble methods, ignoring`. Worse, a
version-2 control file in `regularization` mode does not parse at all — `pestpp-ies`
exits with a control-file error before any model run — so `cal.run_ies()` refuses such
a control file up front rather than letting it fail obscurely. For IES, use
`cal.run_ies(ies_reg_factor=0.1)`, and put spatial structure in the prior via
`correlation=`.

> **Tip from the workshop:** build and inspect the prior *before* history
> matching ("prior Monte Carlo, early and often") — see §5.

---

## 2. The likelihood — how you measure fit (observations & weights)

The **objective function** (phi, Φ) is the weighted sum of squared residuals —
the model's misfit to the data. Minimizing it is "calibration." Each
observation's **weight** scales its influence.

The deceptively hard part is the weights:
- A good starting rule is **weight = 1 / (measurement noise std)** — better data
  counts more.
- But raw weights let high-magnitude or densely-sampled observations dominate.
  The fix is **visibility weighting**: group observations by type and set weights
  so each *group* contributes roughly equally (or by your chosen split). This
  focuses the fit on what matters for your purpose.

Observations are first-class myflopy objects; you hand them to the project:

```python
cal.observe(head_targets)          # history-matching targets (measured values + weights)
```

**Where a weight comes from** (ledger 504): a reading's OWN weight -- a `weight`
column in a LONG values table (`time`, `name`, value, `weight`) -- beats its location's
(the locations' `weight` column, for head, concentration and DRN/RIV-zone targets), which
beats the set's (`LakeStageObservationSpec(weight=)`, `SfrStageObservationSpec(weight=)`,
`SfrFlowObservationSpec(weight=)`), which beats 1.0. A blank weight means "none of its
own". A WIDE values table cannot carry one and is refused naming the LONG form.

**What the control file holds is exactly your targets** (ledger 340). Each target is
matched to the simulated row whose index equals its `time` read as a zero-based stress
period — `1.0` is period 1 — and `build()` REFUSES, naming the target, a key with a
fraction (`1.5`), a negative key, a DATE (ledger 499: pandas reads a datetime
as nanoseconds since 1970, so a table read with `parse_dates` was keyed by period
1706659200000000000 -- bin dated readings to stress periods first), two readings of one
target in one period, a value that is not a number, and a target at a period (or a
name) the simulated table does not hold. `compare()` and `stats()` refuse the same keys
on every route, the file routes included (ledger 500: there a period past the
run's last read `NaN`, and `stats()` left it out in silence). A blank value is
unmeasured. The base run's heads and concentrations are read from the file the model's OC
declares (`HEAD FILEOUT`/`CONCENTRATION FILEOUT`, ledger 505), as every forward
run's are; a declared file missing from the run directory is refused rather than another
read in its place. A base-run table with a cell holding no number is refused too,
whatever that observation's weight (ledger 506): PEST++ fails a run on reading an
empty field at weight 0 as at weight 1, and every realization would hand it the same one. Every simulated row no target measures stays in the control file
at weight 0, holding the base run's own value: reported in every ensemble, never
history-matched. (Until then each such row carried weight 1, pulling the calibration
back toward the uncalibrated model at every period nobody measured.)

**A flow target is in MF6's own sign** (ledger 350): an SFR reach's streamflow
(`SfrFlowTargets`) and a DRN zone's discharge (`DrnFlowTargets`) are water LEAVING the
reach or the aquifer, which MF6 writes negative, on both routes and in the tables PEST
reads — so a gauged `Q` is entered as `-Q`. A value above zero can never be matched and
is refused when the targets are built, naming `-Q`; zero (a dry reach or spring) is
accepted.

**Creek baseflow: `RivFlowTargets`** (ledger 534). A river boundary's MF6
`RIV` flow is the exchange between aquifer and river: NEGATIVE where the aquifer feeds
the river (a gaining reach — baseflow), POSITIVE where the river leaks in. So a gauged
baseflow gain `Q` is entered `-Q`, as for `SfrFlowTargets`, and a net-losing zone is a
positive number; no sign is refused. A zone is a set of cells, drawn or listed exactly
as for `DrnFlowTargets` — the two share one implementation:

```python
baseflow = mf.RivFlowTargets(
    locations=gauged_reaches,                      # one polygon per gauged reach, or {"name": cells}
    values={"time": [0], "rock_creek": [-208_000.0], "crisp_creek": [-240_000.0]},
)
cal.observe(baseflow)                              # group "riv_flow" unless the zones carry a label
```

| zone flow (`DrnFlowTargets`, `RivFlowTargets`) | what it reads |
|---|---|
| budget route (no OBS declared — the default) | the model's cell-budget file, by the package NAME MF6 writes on each record: `package="crisp"` reads that package alone; `package=None` sums EVERY package of the type in the zone's cells (`riv` + `crisp`). Every package read must be declared `save_flows=True`, or reading raises naming it (ledger 536). |
| observation route | the package's declared `DRN`/`RIV` observations (`to_flopy_obs` / `attach_flopy_obs`), one per zone cell; with several packages of the type, `package=` is required (ledger 352). |
| a zone drawn with no `layer` | every layer of the columns it covers (ledger 535; it read layer 0 alone). A `layer` picks one; an explicit `cells` list with no layer is layer 0. |
| a period with no boundary in the zone's cells | 0.0, on both routes (ledger 537). |

The budget route is the cheap one for a big model (Ten Trails `ss`: 0.02 s a set) and
needs no OBS declaration; declare RIV observations only to cross-check. A zone by RIV
BOUNDNAME is not offered (ledger 539): draw the creek or list its cells.
Every named-series set (lake, SFR, DRN, RIV) is regenerated by ONE forward-run call that
opens the run once (ledger 540: it reloaded the whole simulation once per
set, about 7 s of a 17-20 s Ten Trails realization each).

**A dry lake is calibrated against its bottom** (ledger 351). When a lake empties,
MF6 holds its stage at the lake's bottom -- the lowest cell top under a vertical
connection, or a horizontal connection's bottom elevation -- and `LakeStageTargets` hands
PEST that number, on both of its routes, so the simulated stage falls continuously to
the bottom and stays there while the lake is dry. (It was 1e30 on the route reading the
stage file, which PESTPP-IES could not recover from: its lambda stayed near 1e57 and it
stopped calibrating.) Only a lake set `STATUS INACTIVE` has no stage, and is written
1e30. So a lake you OBSERVED dry is a target AT its bottom: `lake_bottoms(model)` in
`myflopy.modflow.mf6.package_tables` gives the number MF6 uses, and a simulated stage can
never go below it. **A dry stream reach is NOT treated that way** (ledger 353): with
`SfrStageTargets` reading the stage FILE, a reach that goes dry reads MF6's -1e30 in the
tables PEST reads -- the magnitude that stalled PESTPP-IES for a lake (ledger 351) --
while on its observation route MF6's STAGE observation gives the streambed top. The SFR
stage PICTURES leave a dry reach blank -- the targets' own `.plot` pictures too (ledger
355) -- while the targets' numbers were left as they were.

> Weighting matters *a lot* and there is no one-size-fits-all. Pre-process data
> for information content, use more groups rather than fewer, and aim weights at
> the model's purpose.

### Observation groups: your `group` labels (ledger 358)

**A target's group is its PEST observation group.** `HeadTargets`, `ConcTargets`,
`DrnFlowTargets` and `RivFlowTargets` read a `group` column from their locations (`group_column=`; `groups=`
on `HeadTargets.from_cells`), and every MEASURED observation of a target is put in the
group its label names:

| what the observation is | its PEST group (`obgnme`) |
|---|---|
| measured, its target labelled `"North Wells"` | `north_wells` — lower-cased, every run of characters other than letters, digits, `_`, `.`, `:`, `-` one `_` |
| measured, no label (missing or blank) | the set's prefix: `hds`, `conc`, `lak_stage`, `sfr_stage`, `sfr_flow`, `drn_flow`, `riv_flow` (or the `prefix=` you passed) |
| a period no target measures, a blank target | `<prefix>_unmeasured` (weight 0) |
| a forecast | the forecast set's prefix (`fore1`, ...); the set's other periods `fore1_unmeasured` |

So PEST++ reports phi by YOUR groups — and, with no labels at all, by data type, since
each set is its own group. Zero-weight rows and forecasts never join a group of yours,
so a group's membership is exactly what was measured under its label. Lake and SFR
targets carry no label: their measured observations are in their set's prefix, and two
sets registered under two prefixes (`cal.observe(north_lakes, prefix="north_lakes")`)
are two groups. (pyEMU used to name every simulated column a group of its own,
`oname:hds_otype:lst_usecol:w1`, so a group was one well.)

`build()` refuses a label it cannot use, naming it: one that is **not text** (a group
column of numbers — use `frame["group"].astype("string")`, or
`.astype("Int64").astype("string")` for whole-number ids, or names like `"zone1"`); the
text **`nan`, `None` or `<NA>`**, which is what `astype(str)` writes for a MISSING label
and would make one group of every well you left unlabelled (ledger 359); one with no
letter or digit, one taking a group the build gives a set (`hds`,
`lak_stage_unmeasured`, ...), and **two different labels that become one name**
(`"North"` and `"north"`, `"north wells"` and `"north_wells"` — one is very likely a
typo). Surrounding whitespace is not part of a label. A DRN or RIV zone given as several
rows must carry one label. **A group whose name starts with `l_`, `less`, `g_` or
`greater` is refused too**, and so is a set prefix that does: PESTPP-IES reads such a
group as an INEQUALITY (`less...`: simulated at most the target; `greater...`: at
least), matching on the bare start, so `"Greater Seattle"` or `"lessard"` would score 0
wherever the model is on the allowed side of the target — measured,
`docs/reproducers/pest_observation_groups.py`. (PESTPP-GLM at NOPTMAX=0 reads them as
ordinary groups; PESTPP-IES does not.)

**Using the groups with PEST++.** Phi by group is in `pestpp-glm`'s `.rec`
("Contribution to phi from observation group ...") and in `pestpp-ies`'s
`<case>.phi.group.csv`, one column per group, and `ies.phi_contributions()` /
`plot_phi_contributions()` / `plot_conflict()` read the same groups. To balance data
types, give PESTPP-IES a **phi-factor file** — one `tag,fraction` line per group — and
it rescales the weights so each tag's share of phi is its fraction:

```python
wells = gpd.read_file("wells.gpkg")          # its `group` column: "shallow" / "deep"
cal.observe(mf.HeadTargets(locations=wells, values=heads))
cal.observe(mf.DrnFlowTargets(locations=springs, values=seepage))   # springs' group: "springs"
cal.build()
(cal.template_workspace / "phi_factors.csv").write_text("shallow,0.4\ndeep,0.4\nsprings,0.2\n")
ies = cal.run_ies(reals=50, iterations=3, ies_phi_factor_file="phi_factors.csv")
```

Measured on PEST++ 5.2.16: a tag matches ANYWHERE in a group's name (`orth` took
`north`), every NONZERO-weight group must be matched by exactly one tag — an untagged one
stops IES ("group 'SOUTH' not identified with any tags"), and so does a group two tags
match (`north` and `nor`: "mapped to multiple tags") — while the zero-weight
`*_unmeasured` groups need none. So name groups so that none contains another's name
(`north` and `north_deep` cannot be given separate fractions). **No default prefix
contains another**, so `hds`, `conc`, `lak_stage`, `sfr_stage`, `sfr_flow`, `drn_flow`
and `riv_flow` each take a tag of their own. A lake set's default was `stage` until ledger 372,
which is inside an SFR set's `sfr_stage`: every tag matching the lake's group matched the
stream's too and PESTPP-IES stopped ("group 'SFR_STAGE' mapped to multiple tags",
ledger 359). A template built before then keeps `stage` in its observation names and
groups -- and, like every template built before ledger 340, must be rebuilt anyway.

### Noise and structural error

Two error sources live in the residuals: (1) **measurement noise** (instruments,
human error) and (2) **structural error** — the model is a simplification, so
even "perfect" parameters can't match reality exactly. Pretending structural
error is zero causes **overfitting**: the model contorts parameters to chase
noise, which *looks* great in the calibration period and predicts badly.

PESTPP-IES carries noise into the answer by fitting each realization to a
slightly different, noise-perturbed copy of the data — so the uncertainty in the
data flows through to the uncertainty in the forecast.

---

## 3. Forecasts — the things you actually care about

A forecast is a model output you want to *predict* but have no data for (a future
flux, an unseen head, a travel time). Declare them so the uncertainty tools track
them. Two workshop habits:
- **Always add a forecast period** to the end of the history model and monitor
  "future" outputs.
- **Frame the question to play to your strengths**: "*Will* the spring flow drop
  below X?" (a probability) is far more answerable than "*What* will the flow be?"

```python
cal.forecast(prediction_targets)   # same target objects as observe(), but zero-weight
```

myflopy registers these as zero-weight PEST++ forecasts, so every uncertainty
method reports their prior→posterior distribution for free. A forecast is the rows its
targets name (a blank value is fine; a number is written as a placeholder); the set's
other periods are zero-weighted and are not forecasts. A set given no values at all
forecasts every row: its locations at every period.

---

## 4. Build and sanity-check

```python
print(cal.settings())              # review the resolved config — never a black box
pst = cal.build("calib.pst", noptmax=0)   # noptmax=0 = "run once, compute residuals"
```

`noptmax=0` is the universal first check: run the model once through PEST and make
sure the plumbing works before launching anything long.

**The build only reads your run** (ledger 356). It copies the simulation's INPUT
files -- what `mfsim.nam` leads MODFLOW 6 to read, and nothing else of the run
directory: no outputs, no run metadata, no old PEST folder (ledger 383) -- to a
staging directory beside the template, loads the simulation from that copy, stores
every array and list there in external files (the form pyEMU's multipliers act on)
and hands that to PEST; the run's files and the model object you called
`model.pest(...)` on are left exactly as they were, and the staging copy is removed
when the build ends. Two consequences: the model must have been **written and run**
(its directory holds `mfsim.nam`, and its output files are where the base
observations come from) -- a model never written is refused -- and the template is
made from the run's FILES, so an edit made to the model object in memory after the
run does not reach PEST: write and re-run first. pyEMU's `PstFrom.log` lands in the
template, not in your working directory.

**What a build deletes** (ledger 360). Every build deletes and rewrites the template
directory -- `<model workspace>.pest/<name>` by default, or the `workspace=` you pass --
and the staging directory beside it (`<template>.staging`). It deletes each **only when
it is absent, empty, or one a previous build wrote there**: a build writes a
`myflopy_pest_build.json` into it, naming the directory -- a visible file, so a folder
you empty (`rm <template>/*`) loses it with the rest. So re-running `cal.build()`, or
a fresh `model.pest(name)` on the same workspace, rebuilds the template -- the output of a
serial `run_ies`/`prior` in it included -- while any other folder (a Project root, a data
folder, a copy of a template such as an IES master directory) is refused with
`FileExistsError`, naming it and what it holds, before anything is deleted. A
`workspace=` that is, contains or sits inside the model's run directory is refused with
`ValueError` whatever it holds, at `model.pest(...)`. Leave `workspace=` out, or pass a
new or empty folder. A template built before this rule has no marker: it is rebuilt when
it is this calibration's, at the place `model.pest(name)` puts it (ledger 361: a copy of
it anywhere else, a backup under the same name included, is not); delete any other
yourself.

---

## 5. History matching with ensembles (PESTPP-IES)

Older calibration (PESTPP-GLM) finds *one* best parameter set using gradients
("one model run per parameter"). **Ensemble methods** (PESTPP-IES) instead evolve
a *collection* of randomly-sampled parameter sets (an **ensemble**), approximating
the gradients from the ensemble itself — so the cost is independent of the number
of parameters, and you get uncertainty for free.

### Prior Monte Carlo — "early and often"

Before history matching, run the *prior* ensemble once and look at it: does the
spread of simulated values **bracket** the observed data? If the observations
fall outside the prior, you have **prior-data conflict** — the model as built
can't reproduce the data with plausible parameters. That's a modeling problem to
fix (structure, bounds, weights), not something to calibrate away.

```python
prior = cal.prior(reals=50)      # run the prior ensemble once (NOPTMAX=-1), no iterations
prior.plot_prior_vs_obs()        # grey prior spaghetti vs measured (conflicts highlighted orange)
prior.conflict()                 # table: which observations the prior can't reach
prior.plot_conflict()            # % of observations in prior-data conflict, by group
```

If the grey envelope brackets the red markers, you have a green light to history
match. If not, fix the conflict first — calibrating into prior-data conflict
injects bias.

### Run it

```python
ies = cal.run_ies(reals=50, iterations=3, workers=10)   # one line; sensible defaults
print(ies.settings)
```

Defaults make `cal.run_ies()` a valid first call. Production knobs (localization,
covariance re-inflation, bad-phi filtering) pass straight through as keyword
arguments when you need them.

**Defaults set for long runs on a shared machine.**

- `overdue_giveup_fac=5.0` (ledger 510). In a parallel run PESTPP-IES kills a
  model run that takes this many times the average run, and drops its realization.
  Left unset, PESTPP-IES uses 2.0, and on a loaded machine that culls a NEWTON model's
  slow, hard-to-converge realizations selectively. Pass a larger value to be more
  tolerant, or `None` for PESTPP-IES's own 2.0.
- Every multiplier file is applied in the forward run's own process (`PstFrom`'s
  `chunk_len=10**6`, ledger 511). pyEMU's default of 50 starts a
  `cpu_count`-process pool in every realization once 50 or more files are
  parameterized. Monthly recharge arrays and per-period river files reach that
  quickly. Measured on a 61-file toy: 21 child processes and 2.9 s per apply,
  against none and 1.7 s.

**Reusing an ensemble** (ledger 509). An `ies_par_en=` you pass is used as
given. When a `grid` or array `pilotpoints` parameter carries `correlation=`, a call
that passes no `ies_par_en` draws a fresh correlated prior into it. Pass
`ies_restart_obs_en=` (a prior Monte Carlo's simulated values, which PESTPP-IES then
does not rerun) together with the `ies_par_en` that produced them, for example
`cal.run_ies(ies_par_en=prior_dir / "cal.0.par.csv", ies_restart_obs_en=prior_dir / "cal.0.obs.csv")`.
A restart ensemble set without one, while a prior would be drawn, is refused: it
would pair every realization's observations with other parameters. The option
persists on `cal.pst`, so `del cal.pst.pestpp_options["ies_restart_obs_en"]` drops it.

**One run's outputs per directory** (ledger 507). A serial run (`workers`
unset) writes `<name>.*` into the template. Before it starts, it deletes the outputs an
earlier run of the calibration left there: each iteration's ensembles and
`.base.par`, the phi files and the noise ensemble, except a file a PEST++ option
names. A rerun with fewer iterations, or a `prior()` after a run, is therefore not
reviewed with the earlier run's last iteration as its posterior. A serial run's
results last until the next serial run, and the deletion is said with a `UserWarning`
(`overwrite=True` says it is meant and silences it; ledger 532). A parallel
run's master is copied from the template without those outputs. The review also
counts only the iterations the run's phi file records, so a directory filled
before this fix also reads correctly.

**Where a parallel run writes, and what it deletes** (ledger 361, ledger 532). With
`workers=N` the results are collected in the master directory -- `<name>_ies_master` (or
`<name>_prior_master`) beside the template, or the `master_dir=` you pass -- and the
workers run in `<master>.workers/worker_0` .. `worker_<N-1>`. Each master has its own
worker directory, so two calibrations of one model, or a `prior()` in one kernel and a
`run_ies()` in another, never share one. A launch deletes the master and the workers it
uses, then copies the template into them, so myflopy checks them first:

- each must be absent, empty, or a COPY of this calibration's template; a folder of
  yours, or another calibration's copy, is refused with `FileExistsError`, and one that
  is, contains or sits inside the run directory or the template with `ValueError`;
- **a master holding an earlier run's results is refused unless `overwrite=True`**,
  naming the folder, whether its run finished, and the spellings -- `overwrite=True` to
  replace it, another `master_dir=` to keep it, `cal.resume_ies(workers=N)` to continue
  an interrupted history match;
- **a directory a run is still going in is refused whatever `overwrite` says**
  (`RuntimeError`). Every launch writes `myflopy_pest_run.lock` into the directories it
  runs in, naming the host, the process and when it started, and removes it when the run
  returns; a lock naming a process still running, or written on another host, or a PEST++
  process working in the folder (read from `/proc` on Linux), means a run is going. A
  lock whose process is gone -- the kernel was killed -- is STALE: it is said so with a
  `UserWarning` and ignored. `build()` refuses a template a serial run is going in.

When the run finishes cleanly its worker directories are deleted -- each is a full copy
of the template, disposable, and a resume makes them again; `keep_workers=True` keeps
them. A run that fails or is interrupted keeps them for diagnosis. The master and its
results are always kept.

**Resuming an interrupted history match** (ledger 533). A multi-day run stopped
part-way (a killed kernel, a reboot) keeps every iteration it completed. In a fresh
kernel, declare and `build()` the calibration again as before -- the template is rebuilt,
the master is not -- and continue it:

```python
ies = cal.resume_ies(workers=10)              # parallel: the master <name>_ies_master
ies = cal.resume_ies(workers=10, iterations=5) # or go further than the run asked
```

PESTPP-IES refuses the PEST `/r` restart switch, so this restarts it from the run's own
ensembles: the prior for regularization, the last completed iteration's parameters and
simulated values (not rerun), and the noise ensemble -- and from the Marquardt lambda the
run's `.rec` shows in force, without which the resumed run takes another path (measured on
the canonical testing model: the resumed iteration 2 read phi 1.715 against the
uninterrupted run's 1.714, and 1.651 without the lambda). The run's earlier outputs are moved
to `<name>.resumed_from_<n>/` while it goes and joined back after, the resumed
iterations renumbered after them, and the run's control file is put back as it was
launched (`NOPTMAX` the run's total), so the directory holds one run of iterations `0..N`
and reviews like any other -- interrupted again, it resumes again to the same total. The
archive keeps the interrupted run's `.rec` and `.log`; a later serial run in the template
deletes it with the other earlier outputs.
The resumed run's workers are made from the template, so a template rebuilt with other
parameters or observations is refused. A SERIAL run (`resume_ies()` with no `workers`)
resumes only while its template has not been rebuilt: `build()` deletes the template and
with it a serial run's results (ledger 360) -- run long calibrations in parallel. A prior
Monte Carlo cannot be resumed: its simulated values are written only at the end of the
sweep.

### The DBTL loop

Real calibration is iterative: **Design → Build → Test → Learn**, repeated. Each
cycle you revise parameters, prior, weights, or the model itself and re-run. The
recommended cadence:

| Stage | Settings |
|---|---|
| Sanity | `noptmax=0`, plot & check |
| First real run | 30–50 reals, `iterations=3`, `bad_phi_sigma=1.5` |
| Struggling to fit | 100–200 reals, more iterations, add multimodal/re-inflation knobs |

```python
ies = cal.run_ies(reals=50, iterations=3, bad_phi_sigma=1.5, workers=10)
```

---

## 6. Assess a run — the per-cycle checklist

Each cycle, look at these. myflopy gives one method per item:

```python
ies.plot_phi()                  # did misfit drop? quickly early but not TOO quickly (overfit)?
ies.plot_phi_distribution()     # spread of phi: prior vs posterior — collapse to ~0 = overfitting
ies.plot_vs_obs()               # prior(grey)/posterior(blue) vs measured(red), one panel per location
ies.plot_vs_obs(groups=["deep"])  # ...only the locations measured in group "deep" (the name the pst holds: "north_wells" for "North Wells")
ies.plot_phi_contributions()    # which observation groups make up the misfit (visibility weighting)
ies.plot_parameters_at_bounds() # % of parameters pinned at bounds — prior too tight / compensating?
ies.forecasts()                 # prior→posterior uncertainty table for every forecast
ies.forecast("spring").plot()   # the payoff: posterior forecast distribution
ies.plot_field("k", stat="mean")     # property patterns — plausible or laughable?
ies.plot_field("k", stat="std")      # where is K still uncertain?
ies.plot_field("k", stat="change")   # where did calibration move K? (posterior/prior)
ies.plot_field("k", stat="reduction")  # did the data inform this region at all?
ies.plot_field_mosaic("k", stat="mean")  # prior vs posterior, side by side, one scale
ies.plot_obs_residuals()        # WHERE is the model biased, and by how much?
ies.plot_sensitivity()          # which parameter groups did the data inform?
ies.plot_sensitivity(forecast="pred_1")   # ...and which ones drive this prediction?
ies.best()                      # the single realization to carry forward (the "base" / min-error-variance one)
ies.report("review.html")       # the headline plots bundled into one HTML
```

Every plot takes `backend="matplotlib"` for static matplotlib/seaborn output instead
of interactive Plotly — except `plot_field_mosaic`, which composes Plotly subplots and
raises rather than pretending. `report(...)` bundles the headline set — phi, obs,
bounds, forecast and `plot_field` mean/change — not every figure on the class: the
reduction map, the mosaic, the residual map, `plot_prior_vs_obs` and `plot_conflict`
are called directly.

The scoring calls (`phi_contributions`, `obs_residuals`, `plot_obs_residuals`) read
ONE posterior realization, `realization="base"` by default, and a label the ensemble
lacks RAISES naming some it has -- they scored the ensemble mean in silence until
ledger 480. A run with no `base` realization (`ies_include_base` off) needs a
label, or `best(criterion="min_phi")`; `best()` itself refuses `"base"` there. The
backend switch is the house one: `"plotly"`/`"interactive"` or
`"matplotlib"`/`"mpl"`/`"static"` (`"go"` and `"sns"` are no longer spellings).

**Reading the field maps.** They are colored by policy, and the policy depends on
the statistic, not the parameter:

- **Statistics in the parameters' space** (ledger 514). When every parameter
  of the field is log-transformed (the default) and every value is positive,
  `field()` takes its statistics in log10, the space PESTPP-IES estimates in. The
  mean is GEOMETRIC, the std is in decades, `change` is the ratio of geometric means,
  and `reduction` agrees with `sensitivity()`'s `learned`. Arithmetic statistics of a
  log-normal K read an unchanged median as roughly halved wherever the data
  narrowed the spread (measured: change 0.54 against a geometric 1.00).
  `field(...).attrs["space"]` names the space, and the map title says
  "geometric mean"/"std of log10".
- `stat="mean"` / `"base"` — absolute conductivity, which spans decades, so it is
  drawn on a **log** scale with the colorbar relabeled back into real units
  (`0.001`, `0.01`, …) on *both* backends. A cell whose K is exactly 0 has no
  logarithm and draws as a gap. (Since ledger 252 every log map reads this way,
  hand-built or by policy; `plot_field`'s only extra is the `x` suffix on the
  `change` map's ratios.)
- `stat="std"` — a spread (in decades for a log field, model units otherwise), drawn
  **linear**: a posterior standard deviation is legitimately zero wherever the ensemble collapsed, and those are
  precisely the cells this map exists to show.
- `stat="change"` — `posterior_mean / prior_mean` (geometric means for a log field), a **ratio** whose neutral value
  is 1, so it is mapped as `log10(ratio)` on a diverging red-white-blue scale with
  symmetric limits: white is "calibration left this alone", **red is reduced**,
  blue is increased, and a halving sits exactly as far from white as a doubling.
  The hover shows the raw ratio (`0.5x`), not the logarithm.

- `stat="reduction"` — `1 - posterior_sd / prior_sd`, the share of the prior spread
  the data removed. **This is the "did the data inform this region" map.** Its scale
  is anchored to 0–1 (0 = the data said nothing here, 1 = the ensemble collapsed)
  rather than autoscaled, so two layers are comparable and a field that only ever
  reduces 0.95–1.0 does not stretch into a dramatic-looking map of a trivial range.
  If any cell went *negative* — the posterior spread grew — the scale falls back to
  the data range, so those cells stay visible instead of clipping to the bottom
  color where they would read as "no reduction".

Hovering a cell gives the plotted statistic plus prior mean, posterior mean, and
posterior standard deviation; the title states which iterations and how many
realizations the figure is summarizing. A cell whose prior mean is zero has an
undefined ratio and says so, rather than silently reading as missing data.

**What did the data actually inform, and what drives the forecast?**
`sensitivity()` answers both from the ensembles the run already wrote — no jacobian,
so no extra model runs. `learned` is `1 - posterior_sd / prior_sd` per parameter group
(measured in the parameter's transform space): a group near zero was not constrained by
the observations you have. Passing `forecast=` adds the ensemble correlation between
each group and that prediction.

The pairing is what to read. A group with **high forecast correlation and low
learning** is where your forecast uncertainty comes from: it matters, and the data did
not pin it down. That is the argument for collecting a different kind of measurement.

Two honest caveats, both worth taking seriously:

- **This is not CSS.** Composite scaled sensitivity is a local derivative at one
  parameter set; this is a global measure conditioned on your prior. A parameter can be
  locally sensitive yet score near zero here because the prior never moved it far. The
  jacobian-based answers (`Schur`/`ErrVar` identifiability, FOSM) need a PESTPP-GLM run,
  which myflopy does not launch — an IES run produces no jacobian.
- **Small values are noise.** With `n` realizations, a correlation around `1/sqrt(n)` is
  indistinguishable from zero — at `reals=50` that is ~0.14. `plot_sensitivity` draws
  that floor on the chart, and `n_reals` is in the returned frame; do not rank bars
  below it against each other.

**Prior vs posterior, side by side.** `plot_field_mosaic("k", stat="mean")` draws one
panel per ensemble on **one shared color scale** — which is the whole point, since two
autoscaled panels would make any difference look identical. The panels also pan and
zoom together — though only in a saved/shown HTML page, not in inline notebook output,
which gets the shared *starting* view but no live linking. Only `stat="mean"` and
`stat="std"` are accepted: they are the only stats with a separate prior and posterior
form. `change` and `reduction` are already prior→posterior comparisons, and `base` is a
single posterior realization with no prior counterpart, so composing any of them "prior
vs posterior" would draw one map twice. Plotly only — for static output, call
`plot_field(..., backend="matplotlib")` once per panel.

**Where is the model actually wrong?** `plot_obs_residuals()` puts the misfit on the
grid rather than in a histogram: point targets — head *and concentration* — draw as
markers at their coordinates, DRN
zones color the cells they cover, and **both read on one diverging scale centered on
zero** — so a point and the cell under it mean the same thing at the same color. **Red
is under-simulated**, blue is over-simulated, white is on the money. The residual is
`simulated − measured` (the convention `phi_contributions` uses; note PEST's own `.res`
file reports the opposite sign, and every label on this figure says which one it is
in). A DRN zone's seepage is in MF6's own sign, negative -- water leaving the aquifer,
entered as `-Q` (ledger 350) -- so on a zone "over-simulated" (blue) means it drains
LESS than was measured. Each location is summarized over time, so a target measured in six stress periods
is one point showing its mean bias, with the count in the hover. `obs_residuals()`
returns the same numbers as a DataFrame.

Lake and SFR targets are *not* on this map: they record only a lake or reach number, with
no geometry, so they cannot be placed from the run directory alone. They are absent from
`obs_residuals()` too, so the figure never silently implies it covered them. What decides
this is the observation set's declared **geometry** (`"points"` / `"zones"`), not its
kind — a family that records no geometry is skipped.

**One scale, so mind the units.** A run history-matched against more than one *kind* of
measurement — heads in length, concentration in mass/volume, DRN seepage in
volume/time — has residuals in different units, and one shared color limit hands the
scale to whichever family has the bigger numbers. The others then render uniformly
white, which reads as a perfect fit rather than as an unreadable figure. Pass
`prefix=` to draw one family at a time (`plot_obs_residuals(prefix="conc")`); drawing
them together still works and warns. `obs_residuals()` takes the same `prefix=`, and
its `value_kind` column names which measurement each row is.

Two other things are deliberately excluded, both of which would otherwise be drawn as
misfit the model was never fitted to. **Forecasts**, because `cal.forecast(...)` registers
an ordinary observation set at weight 0 and nothing in the saved metadata distinguishes it
— a forecast's apparent "residual" is often the largest number in the run and would set
the whole map's color scale. And **unmeasured times**: pyEMU creates one observation per
row of the simulated output, so a target measured in one stress period out of six leaves
five rows carrying the model's own output -- at weight 0 in a run built since ledger 340,
at weight 1.0 in one built before it. Averaging those in drags every residual toward zero
and can flip its sign, so they are excluded by joining back to the
`*_target_values.csv` snapshot of what was actually measured (which, since ledger 340,
names each target's observation).

What you're checking for:
- **Phi**: should drop fast at first. If it collapses to near-zero, you're
  overfitting (see §7).
- **Obs vs sim**: the posterior spread should *cover* the measured values without
  being implausibly narrow.
- **Forecasts**: value the *spread*, not a tiny phi. Use `capture=True` on a
  parameter (e.g. `cal.parameterize("k", ..., capture=True)`) to enable the field
  maps. **Capture costs one zero-weight observation per cell of every layer captured.**
  On Ten Trails, `k33` over 5 x 20,323 cells added 101,615 observations and grew the
  template from 18 to 33 MB. Every observation ensemble then carries those columns:
  one 5-layer field at 200 realizations is a 352 MB CSV per iteration, ~2 s per read.
  A capture covers the layers its parameter covers (`layers=`). Consider
  `run_ies(ies_save_binary=True)`; the review reads `.jcb` ensembles (ledger
  512) and reads each file once (ledger 513).
- **Property maps**: if the calibrated K pattern looks physically absurd, adjust
  the geostatistics/prior — don't accept it just because phi dropped.
- **% of parameters at their bounds**: many parameters pinned at bounds means the
  prior is too tight or compensating for structural error.

---

## 7. The most important lesson: a good fit ≠ a good forecast

The workshop's "pepsi challenge": two runs, one with an *okay* fit and one with a
*crazy-good* fit. The crazy-good fit has a tiny, confident posterior — that
*misses the truth* for some forecasts. The okay fit keeps a wider, more honest
posterior that brackets the truth.

With an imperfect model (every model is imperfect), chasing the data too hard
injects **bias** into predictions that depend on parameter combinations the data
can't constrain. The cure is to **deliberately under-fit**: heavier noise/weights,
fewer iterations, filter signal you can't (or shouldn't) match. If prior
uncertainty is already good enough for the decision, you may not need to history
match at all.

Bias–variance in one line: fitting harder lowers variance but can raise bias, and
in groundwater the bias is usually *unknowable* — so err toward honest uncertainty.

---

## 8. Beyond history matching (where this goes next)

- **Linear uncertainty / data worth (FOSM)** — cheap prior→posterior uncertainty
  and "which observations are worth collecting" without an ensemble run.
- **Management optimization** — PESTPP-OPT (linear, chance-constrained) and
  PESTPP-MOU (multi-objective, Pareto trade-offs under uncertainty): "how much can
  we pump while keeping the spring above X at 95% confidence?"
- **Emulation** — fast surrogate models (DSI, GPR) that stand in for a slow model
  during uncertainty/optimization.

In myflopy today, the run-time knobs for the advanced IES features are available
through `run_ies(**pestpp_options)`; FOSM, optimization, and emulation facades are
on the roadmap.

---

## The whole workflow, start to finish

```python
# model.pest(...) is the front door (returns a PestProject); the workspace
# defaults to <model workspace>.pest/<name>
cal = model.pest("calib", start_datetime="2020-01-01")

# 1-2-3: prior, likelihood, forecasts
cal.parameterize("k",        style="constant", bounds=(0.2, 5), physical=(1e-3, 100), capture=True)
cal.parameterize("recharge", style="grid",     bounds=(0.5, 1.5), physical=(0, 1e-2))
cal.observe(head_targets)
cal.forecast(prediction_targets)

# 4: build & sanity-check
print(cal.settings())
cal.build("calib.pst", noptmax=0)

# 5: history match
ies = cal.run_ies(reals=50, iterations=3, bad_phi_sigma=1.5, workers=10)

# 6: assess (the DBTL checklist)
ies.plot_phi(); ies.plot_vs_obs()
ies.forecasts(); ies.forecast("spring").plot()
ies.plot_field("k", stat="change"); ies.plot_field("k", stat="reduction")
ies.plot_field_mosaic("k", stat="mean"); ies.plot_obs_residuals()
ies.report("review.html")

# 7: learn, revise, repeat
```

Every plot above takes `backend="matplotlib"` if you prefer static
matplotlib/seaborn output over interactive Plotly — except
`plot_field_mosaic`, which composes Plotly subplots and refuses.
