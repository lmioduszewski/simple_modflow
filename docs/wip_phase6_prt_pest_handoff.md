# WIP handoff — Phase 6 PRT (§6.3) + PEST-IES (§6.4) integrations

> **TEMPORARY session-handoff document** (2026-07-24). Delete when §6.3+§6.4 are
> complete and their plan banners are written. This captures verified scoping
> results, locked user decisions, and the exact next steps so a fresh session can
> continue mid-stream without re-scoping.

## Where things stand

- **Committed:** `fb79d99` **6.3A DONE** (PRT build-side declarability:
  `mf.mip`/`mf.prp`/`mf.ems`, prt6 dis/disv/oc dispatch, kind-aware
  `mf.simulation` solver default, docstring fixes, e2e MF6 test, docs synced —
  see the §6.3.3 banner in the implementation plan). Pushed.
- **6.3B DONE 2026-07-25** — release groups (`PRTReleasePoints.from_cells/
  from_points(group=)`, `merge`, `PRTProject` boundnames) + `prt_maps.py` with the
  `results.travel_time` / `.endpoints` / `.capture` nouns. NAMING: they are nouns,
  not the `capture_map()`/`travel_time_map()` spellings drafted below — the plan's
  original names were `verb_noun` pairs the view-layer conventions forbid. An
  adversarial review of the change also found and fixed a real pre-existing bug:
  `logscale` was silently dropped on every `custom_zs` choropleth (ledger 60).
- **6.3C DONE 2026-07-25** — `results.pathlines` is now a **view** (the raw CSV
  moved to `results.track_records`, clean break, ledger 61) whose `map()` draws one
  `Scattermap` polyline per particle over `base="heads"`/`None`/an existing `Choro`,
  with per-vertex `pathline_hover`, `plot()` = elevation vs travel time, `mosaic()`
  per release group, `backend="mpl"` = the old FloPy plan view. Two shared-surface
  changes came with it: `Choro.add_overlay`/`overlay_traces` + `viz.mosaic` copying
  overlays (fixes a pre-existing silent drop of contours/locs, ledger 62) and
  `PALETTE.categorical` + memoized `viz.category_colors` (ledger: policy colors for
  named categories), retrofitted into `travel_time.plot()`/`capture.plot()`. New
  grid helper `vor.points_to_latlon(x, y)` — **6.4B's residual points want it too**.
  Canonical notebook 03 cell 5's hand-rolled matplotlib is now the library figure
  (ledger 64). **§6.3 is complete**; the plan's acceptance checklist item is ticked.
- **6.4A DONE 2026-07-25** — `build_choropleth` + `vor.choropleth` widened,
  `_field_map_policy` per-stat colors, `parameter_field_hover`, run context on
  the title, ies.py color literals retired. See the 6.4A block below for the
  four code-established corrections; ledger 67–73.
- **Next:** 6.4B→C.
- Full serial suite after 6.4B: **940 passed / 1 skipped** (`pytest -n0`, ~162 s).
- Task list: #54 6.4A, #55 6.4B, #56 6.4C.

## Locked user decisions (AskUserQuestion, 2026-07-24)

1. **Full §6.3** — all three pieces: derived per-cell maps + NEW plotly pathline
   map + build-side declarability (6.3A done).
2. **Real release groups via boundnames** — `PRTReleasePoints`/`PRTProject` (and
   `mf.prp`, already done) gain group labels; `results.capture.map(by="release_group")`
   reads the track CSV `name` column. NOT the irpt fallback. (DONE, 6.3B.)
3. **Residual map scope: heads + DRN zones now**; lake/SFR residual mapping
   deferred (no coordinates without live-model resolution) — ledger it.

## Verified scoping facts (from a 5-agent sweep + live probes; trust these)

### PRT data layer
- Track CSV schema (`PRTRunResults.track_records` = raw `pd.read_csv`):
  `kper, kstp, imdl, iprp, irpt, ilay, icell, izone, istatus, ireason, trelease,
  t, x, y, z, name`. Particle key = `(imdl, iprp, irpt, trelease)`.
- **`icell` is a ONE-based whole-grid node number** (verified by live 2-layer
  DISV probe): zero-based `cell = (icell-1) % ncpl`, `layer = (icell-1) // ncpl`
  (`== ilay-1`). Nothing in the repo does this conversion yet. The synthetic
  fixture in tests/test_mf6_prt.py:180-199 has non-realistic basing — do not
  trust it for icell semantics.
- `terminal_points` = rows with `ireason == 3` (prt.py:284-291). ireason codes:
  0=release record, 1=cell-to-cell transition, 3=termination. Travel time per
  particle = `t - trelease` on the terminal row.
- Boundnames are echoed **uppercased** into the `name` column (verified:
  `WEST_WELLS`). `mf.prp` (done, package_api.py) auto-enables boundnames from
  6-tuple rows. `PRTProject`/`PRTReleasePoints` do NOT support groups yet —
  that is 6.3B work.
- `PRTRunResults.flow_model` is a `SimulationBase` with `.vor` (ncpl,
  gdf_vorPolys) and `.hds` — everything a Choro needs (prt.py:256-262).

### Derived-map recipe (the fixed house pattern — copy it)
Precedent: `SurfaceWaterExchangeResultsExplorer`
(package_surface_water.py:1832-1950) — a NON-registry SpatialView over a derived
table. Recipe:
1. Normalized per-cell DataFrame with zero-based `cell` column + value column +
   extra hover columns.
2. `build_cell_input_map_payload(frame, ncpl=, value_column=, per=, layer=,
   agg=)` (package_plotting.py:381-470) → `(values list[ncpl], hover dict)`;
   **passes ALL extra frame columns into the hover dict** (numeric agg="first").
3. `flow_model.cor(per=, layer=, type='custom', custom_zs=values,
   custom_hover=hover, hover_heads=False, hover_ks=False,
   hover_spec=result_hover(...), colorscale='earth', logscale=...)` then
   `_apply_backend(choro, backend)`.
- `type='custom'` → `Choro.zs` returns custom_zs verbatim (log10 when
  logscale=True). custom_zs must be a **plain Python list exactly ncpl long**.
- What custom maps lose: per-layer hover table + top/botm surfaces (fine);
  contours still work. Choro requires `model.hds` to exist (kstpkper default,
  choros.py:190) — PRT maps build against `flow_model`, which has it.
- Choro colorscale setter whitelists 17 names; `'earth_r'` etc. silently fall
  back — reversal needs explicit `[[pos,color],...]` stops or
  `reversescale=True` through kwargs (choros.py:693-716, 928-940).
- `result_hover` (hover.py:621-638) has **NO footer= param** (hardwired
  `('period','date')`) — ADD footer= (non-breaking, default unchanged); PRT
  time-integrated maps pass `footer=()`.
- `SpatialView` (package_plotting.py:788+): host provides tabular `get()`
  (long frame w/ per/layer/cell + value col) + atomic `map()`; mosaic facets
  ONLY by 'layer'/'model' → **capture map must hand-build (label, Choro) panels
  and call `viz.mosaic` directly** (viz.py:173-304; shared coloraxis comes from
  the FIRST panel; `sync_views=True` syncs pan/zoom).
- Import layering: prt.py is layer 0 with 4 deferred imports (allowlist);
  package_plotting is layer 1 → **the view classes go in a NEW module**
  (`src/myflopy/modflow/mf6/prt_maps.py`), and prt.py lazily imports it inside
  the noun properties via ONE shared helper (+1 to prt.py's deferred count —
  regen `scripts/derive_import_layers.py`).

### PEST-IES layer (for 6.4)
- ies.py viz audit **already passes**: every plot method + IesForecast.plot use
  viz.Fig/viz.subplots/viz.mpl_axes. **STALE as of 6.4A** — the `'darkorange'`,
  `'rgba(80,80,80,0.35)'` and `_POST_COLOR.replace('0.55','1.0')` nits listed here
  were all retired to `viz.PALETTE` in 6.4A. Only two bare `color="0.6"` greys
  remain (ies.py:689, :1139).
- `field(target, layer=0)` (ies.py:1062-1097) already returns per-cell
  `cell, prior_mean, prior_std, posterior_mean, posterior_std, base?(when a
  'base' realization exists), change(=posterior_mean/prior_mean)`.
- **The captured field is ALWAYS absolute resolved K** (base × multipliers,
  clamped, re-read from the rewritten array file) — never a multiplier. So the
  colorscale policy is **per-STAT**: `change` → RdBu log-centered at 1;
  `mean/std/base` → 'earth' (+logscale for K). Today: plotly branch always
  'earth' (unoverridable), mpl branch `RdBu_r` for change else `'viridis'`
  (ies.py:1144-1152) — reconcile both backends.
- **Load-bearing seam:** `build_choropleth` (grid/plotting.py:34-72) has a FIXED
  no-kwargs signature — `plot_field`'s `**choropleth_kwargs` raise TypeError for
  colorscale/logscale/hover_spec/zmid today. DECIDED FIX: extend
  build_choropleth with a `**choro_kwargs` passthrough to Choro. (Choro itself
  supports logscale/named-RdBu/explicit stops/zmid-through-kwargs already:
  choros.py:133,217,642-643,693-716,928-940.)
- Residual map (heads + DRN): head-target locations ARE on disk in the results
  workspace — `<prefix>_target_locations.gpkg` (Point geometry, name/layer/
  group/weight) + `<prefix>_head_target_map.csv` (name,layer,cell)
  (observations.py:185-208,222-235); `myflopy_pest_metadata.json`
  `observation_sets` records each set's kind + locations_file, and
  `IesResults._metadata` already parses that JSON but only exposes
  capture_fields (ies.py:243-259) — expose observation_sets. DRN zones snapshot
  explicit `cells` lists (observations.py:299-318). Lake/SFR have NO coords →
  deferred (ledger).
- Uncertainty: `plot_field(stat='std')` exists; variance reduction
  `1 - post_sd/prior_sd` is one derived column from field().
- `field_mosaic`: `viz.mosaic` composes Choro panels with a shared coloraxis
  (diff=True → zero-centered) + synced views — nearly free.
- Ensembles: `obs_ensemble(i)`/`par_ensemble(i)` read `<case>.<N>.obs.csv`/
  `.par.csv`; prior/posterior = lowest/highest iteration on disk; iterations
  globbed (ies.py:275-342). `.iterations`/`.settings` provide hover-footer info.
- `IesResults.model` is auto-attached on all three entry paths (run_ies/prior/
  pest_runs.review) — the choropleth path is safe.
- Tests: **no synthetic ensemble fixture exists anywhere.** Build a canned
  run-directory fixture: minimal pyemu-loadable `.pst` + `case.0.obs.csv`/
  `case.1.obs.csv` + `.par.csv` files + `phi.actual.csv` +
  `myflopy_pest_metadata.json`; then `open_ies_run(dir, model=fake)` with the
  existing `_two_cell_vor_clockwise` grid (test_mf6_pest.py:104-118) exercises
  field/plot_field/hover with zero pestpp runs. Slow e2e extension point:
  `test_ies_capture_field_and_spatial_maps_end_to_end`
  (test_mf6_pest.py:1397-1455) — already slow-marked + xdist-prioritized (slow
  marking is by NAME in `_SLOW_TESTS`, conftest.py:139-176, NOT decorators).
  pestpp-ies is on PATH (no skip guard; missing binary FAILS not skips).
  There is NO canonical IES e2e test (plan text was wrong).

## Next steps, in order

### 6.3B — release groups + derived cell maps — **DONE 2026-07-25** (kept for the design rationale; the `*_map()` names below became nouns)
1. `result_hover` (hover.py:621): add `footer: tuple = ("period", "date")`
   param, pass through to HoverSpec. Non-breaking.
2. `PRTReleasePoints` (prt.py:129-231): add `group: str | Sequence[str] | None =
   None` to `from_cells` and `from_points` → rows become 6-tuples
   `(irpt, cellid, x, y, z, group)` when given (single str applies to all
   points; sequence is per-cell). Consider a small `merge(*sets)` classmethod
   that renumbers irpt so multiple grouped sets combine into one PRP.
3. `PRTProject.__init__` (prt.py:406-419): pass `boundnames=True` to
   `ModflowPrtprp` when any packagedata row has ≥6 fields.
4. NEW `src/myflopy/modflow/mf6/prt_maps.py`:
   - `pathline_cell_table(pathlines, ncpl)` — vectorized icell→(cell, layer)
     normalization + `travel_time = t - trelease` on terminal rows; release
     group from `name` (may be empty).
   - `PRTTravelTimeView(SpatialView)`: `__init__(results)`; `self.model =
     results.flow_model`; `value_name = "travel_time"`.
     `get(stat="median", layer=None)` → columns `cell, layer, travel_time,
     particle_count, min_time, max_time` (stat picks the aggregation for the
     travel_time column; count/min/max always included → they feed hover
     extra_fields). `summary()`. `map(stat=, layer=None, logscale=False,
     backend=, hover=None, **kwargs)` → recipe above with
     `result_hover("travel_time", title="Travel time",
     units={"travel_time": "d"}, extra_fields=("particle_count", "min_time",
     "max_time"), footer=())`, colorscale 'earth'.
   - `PRTEndpointsView(SpatialView)`: `get()` → `cell, layer, particle_count`;
     `map()` → counts choropleth, 'earth', result_hover("particle_count",
     title="Endpoints", footer=()).
   - `PRTCaptureView`: `get(by="release_group")` → `group, cell, layer,
     particle_count` (group = name column; when all-empty raise a clear error
     naming `group=` on the release points and `by="release_point"` as the
     irpt fallback); `map(by=, ncols=2, sync_views=True)` → per-group endpoint
     Choros hand-composed via `viz.mosaic` (shared scale).
   - layer=None on get/map means aggregate across layers (a capture-zone map is
     layer-agnostic by default); explicit layer= filters.
5. `PRTRunResults` (prt.py): add `travel_time`/`endpoints`/`capture` properties
   via one deferred import of prt_maps (keep the ratchet delta to +1).
6. Tests: NEW `tests/test_prt_maps.py` — module-scoped tiny coupled GWF+PRT run
   (copy the builder from
   `test_package_api.py::test_prt_model_fully_declarable_and_runs`, clockwise
   `_disv_values` winding, 2 release points with groups). Assert: normalization
   math (incl. a synthetic multi-layer frame: icell=3 on ncpl=2 → layer 1 cell
   0), travel-time values sane, endpoint counts, capture groups
   == {WEST_WELLS, EAST_WELLS}, maps render with colorscale 'earth' +
   hover template contains the title, footer has no "Period". Cross-ref comment
   in test_colorscale_policy.py pointing at the PRT pins.
7. Regen import layers + snapshot (prp/PRTReleasePoints signature changes show
   in `mf.__engine__` only if exported — check `derive_api_snapshot.py`
   surfaces); docs sync (package_api_reference read-side, myflopy_context PRT
   derived-maps row ❌→✅, plan 6.3.1 banner, view_layer_conventions if the
   noun pattern needs a line); ledger: result_hover footer default stays
   ('period','date') for back-compat; capture-map error-not-fallback choice.
8. `pytest -n0` full; commit "6.3B: PRT release groups + derived cell maps".

### 6.3C — plotly pathline map (task #53) — **DONE 2026-07-25**
Built as a fourth **noun** (`results.pathlines`), not the `pathline_map()` method
drafted here — same naming rule that renamed 6.3B's maps. It returns the `Choro`
(not a bare `viz.Fig`) so `_apply_backend`, view syncing, and the mosaic composer
all keep working. Corrections to what was drafted, worth carrying into 6.4:
- Maps are `go.Choroplethmap` in **WGS84**, so any overlay of model x/y needs
  reprojection: `vor.points_to_latlon(x, y)` (new, `grid/geometry.py`). 6.4B's
  head-residual points and DRN zones need exactly this.
- `viz.mosaic` copied only `panel.get_choropleth()` — overlays were silently
  dropped. Now `Choro.add_overlay()`/`overlay_traces()` and mosaic copies both
  (ledger 62). Any 6.4 figure that draws over a map should use `add_overlay`.
- `PALETTE` had no qualitative sequence; `viz.category_colors` +
  `PALETTE.categorical` now memoize name → color across figures. Use it for IES
  realization groups/zones rather than new hex literals.
- `HoverSpec` renders fine on a non-choropleth trace: build a
  `HoverContext(ncpl=<points in the trace>)` — the assembler's `ncpl` is a row
  count, not a grid property. `pathline_hover()` is the worked example.

### 6.4A — IES plumbing/policy/hover (task #54) — **DONE 2026-07-25**
Delivered as scoped, with four corrections that a 12-agent scouting/design sweep
established against the code. Carry these into 6.4B/C:
- **Diverging scales ship as STOPS, never the name `'RdBu'`.** `_PLOTLY_TO_MPL_CMAP`
  maps `'rdbu'` → `'RdBu_r'`, the reversed colormap, so a named diverging scale
  renders **mirrored** between `plot()` and `plot_mpl()`. New helper
  `_red_white_blue_diverging_colorscale()` (package_plotting.py), derived from its
  sibling. NOTE the same trap still bites every signed-`q` map in the registry —
  ledger 70, deliberately not fixed here.
- **`zmid` does not center a matplotlib map.** `plot_mpl` reads `self._zmin`/`_zmax`
  and never consults the trace kwargs. Symmetric `zmin`/`zmax` **in log space** is
  what makes both backends agree. No `TwoSlopeNorm` is involved: deleting the
  hardcoded `cmap=` was the whole mpl fix, since `plot_mpl` derives its colormap
  from `self.colorscale`.
- **The hover footer cannot carry run context.** `_render_footer` recognizes only
  period/step/date/area/model with no `else`, so `footer=("iteration",)` renders
  nothing silently. Iteration + realization count ride the figure TITLE instead.
- **Never touch `IesResults.settings` for the realization count** — it does
  `len(self.forecast_names)` → `list(None)` and raises `TypeError` on any
  forecast-less run (ledger 73, still unfixed; fold into 6.4C). Use
  `self.prior._df.shape[0]` / `self.posterior._df.shape[0]`.
- Also widened `VoronoiGridPlus.choropleth` (voronoi.py), a hand-copied twin of
  `build_choropleth`'s signature — widening only the function would have left
  `vor.choropleth(colorscale=...)` raising `TypeError`.
- **A log map's colorbar is labeled in log10 units unless you relabel it** — nothing
  in `choros.py` sets `tickvals`/`ticktext`. 6.4A relabels inside `plot_field` for
  BOTH backends (`_log_decade_colorbar` + `_relabel_log_colorbar`, since `plot_mpl`
  has no tick hook). Ledger 74; a hand-built log `Choro` still reads in log10.
- **`viz.mosaic` discards panel `zmin`/`zmax` AND the colorbar** — it recomputes the
  shared range from the pooled data and only centers under `diff=True`. So 6.4B's
  `field_mosaic` MUST pass `diff=True` for a `change` mosaic or the white "no change"
  point moves off 1x (measured: fraction 0.231, not 0.5). Ledger 71.
- `api_snapshot.json` unchanged (`build_choropleth` was never pinned); import
  layers regenerated, `pest.ies` 1→2 (deferred count unchanged).
### 6.4B — new figures (task #55) — **DONE 2026-07-26**
Three of the four planned names changed, each for a reason established against the
code. Carry these into 6.4C:
- **The uncertainty map is `plot_field(stat="reduction")`, not a new method.**
  `plot_field(stat="std", which=...)` already WAS the posterior-sd map (6.4A), so
  `field_uncertainty_map(which="std"|"reduction")` would have been a second spelling
  of a shipped figure plus a second meaning for `which`. Only `1 - post_sd/prior_sd`
  was new, and `field()` already had `prior_std`. Ledger 75.
- **`plot_field_mosaic`, not `field_mosaic`** — every figure method on `IesResults` is
  `plot_*`. It only accepts `mean`/`std`, which makes ledger 71's `diff=True` trap
  **unreachable** instead of documented: `change` has no prior/posterior pair to
  compose. Panels come from a new `_field_choro` seam (`viz.mosaic` needs the `Choro`,
  `plot_field` returns a rendered figure).
- **`viz.mosaic` gained `colorbar=`** (dict, or a `(cmin, cmax) -> dict` callable —
  the pooled limits exist only inside `mosaic`). Additive: omit it and every existing
  mosaic is byte-identical. Ledger 71 half-retired.
- **`plot_obs_residuals()` has no `map=` flag** — nothing was specified for
  `map=False`, and a boolean with one honest value is worse than a named method.
  Ledger 77.
- **`usecol` from `pst.try_parse_name_metadata()` is unusable as a join key** — it
  truncates at the first underscore (`obs_00` → `obs`). Prefix and location are parsed
  out of the whole obs name by `_OBS_NAME_RE` instead. This is the single fact the
  residual map depends on.
- **`plot_mpl` draws in MODEL coordinates**, so the static backend scatters raw x/y and
  only the Plotly path needs `points_to_latlon`. That is why residuals are not
  plotly-only. `plot_mpl` ignores `Choro` overlays entirely.
- `mpl_colormap_for` extracted from `plot_mpl` into `choros.py` so overlaid points and
  the cells beneath them share one colormap. Ledger 79.
- Lake/SFR residuals deferred — they persist only a lake/reach number. Ledger 76.
- Import layers regenerated: `pest.ies` deferred 1→**0** (hoisted, per the ratchet's
  own advice) and layer 2→4. `api_snapshot.json` unchanged again.
### 6.4C — fixture/tests/docs (task #56)
- **Most of the planned test work landed in 6.4B** — `_stub_ies_results` in
  `test_mf6_pest.py` is a hand-built run directory (location snapshots + a fake
  `pst`/posterior) that drives the whole residual path fast, and the e2e was extended
  in place. What 6.4C still owes: a stub that also carries **capture fields**, so
  `field`/`plot_field`/`plot_field_mosaic` get fast coverage too (today only the slow
  e2e proves them); the last two `ies.py` color literals — bare `color="0.6"` greys at
  **ies.py:689 and :1139** (use `viz.PALETTE.mpl_ensemble`). NOTE: the `'darkorange'`
  literals and the `_POST_COLOR.replace('0.55','1.0')` alpha hack this file used to list
  here were already retired in 6.4A — verified gone, do not go looking for them.
  Also: ledger 73's `settings`/`report()`
  `TypeError` on a forecast-less run; docs (plan §6.4 banner — it still has none;
  snapshot regen if any exported signature changed; full `-n0` suite; commit;
  **push everything** (origin is currently at 7322af9); delete this file in the
  final 6.4 commit.

## Standing rules (apply to every step)
- Python: `/home/lukem/python/envs/gw/bin/python` (bare python/python3 won't work).
- Verify serially `pytest -n0 -q` before every commit (xdist can mask failures).
- Regen `scripts/derive_api_snapshot.py` on any public-signature change and
  `scripts/derive_import_layers.py` on any import-structure change; review diffs.
- Docs in sync IN THE SAME COMMIT: package_api_reference.md, myflopy_context.md,
  implementation-plan banners, CLAUDE.md, view_layer_conventions.md as touched.
- Every deliberate scope cut → docs/compromises_and_deferrals.md (same pass);
  next ledger number is 80.
- Figures always viz.Fig; colors from policy helpers, never hex at call sites.
- Commit messages end with: Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
