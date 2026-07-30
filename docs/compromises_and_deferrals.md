# Compromise & deferral ledger

> **Standing rule (user instruction, 2026-07-17):** every deliberate scope cut,
> test-fidelity trade, judgment call on ambiguous instructions, or optional
> capability that was considered and NOT built must be recorded here **in the
> same pass** as the work that made the call. This ledger is the honest delta
> between "what the plan/sketch implied" and "what shipped". It is NOT a
> roadmap — future work that was never in scope lives in
> `docs/implementation_plan_2026-07.md`; this file only tracks compromises
> *inside delivered work*.
>
> Format per entry: **What** was cut/traded · **Why** · **Impact** · **Revisit
> trigger**. Remove an entry only when the compromise is actually undone (note
> the commit), never because it got stale.

## Fast test suite (2026-07-17, merged `c78c8b9`)

1. **Default suite runs the smallest canonical profile.**
   - What: `pytest` uses `CanonicalModelConfig.testing()` (21×21×4, nper=6);
     the 50×50 `validation` profile runs only in the weekly CI lane
     (`SIMPLE_MODFLOW_CANONICAL_PROFILE=validation`), the 100×100 full profile
     only manually.
   - Why: full suite 11m44s → 47.5s; every profile is contract-complete (all
     15 packages + 5 obs families), so *feature* coverage is unchanged.
   - Impact: day-to-day tests exercise smaller matrices; resolution-dependent
     numerical regressions would surface weekly, not per-push.
   - Revisit: if a bug ever slips through that validation() would have caught,
     promote validation to per-push CI.

2. **IES e2e tests trim the ensemble exploration.**
   - What: `run_ies(reals=6, iterations=1, lambda_scale_fac=1.0, workers=4)`
     in `tests/test_mf6_pest.py` (was the pyEMU defaults; ~82 forward runs →
     ~20).
   - Why: the two IES e2e tests dominated suite wall time (~270s combined).
   - Impact: tests still cover the full master/agent path, phi/ensemble
     mechanics, and review layer, but with one lambda scale factor and a
     small ensemble. The shipped `run_ies` API is untouched.
   - Revisit: weekly validation lane, or if an IES option-handling bug slips
     through.

3. **Pilot-point test density sized to the small domain.**
   - What: `pp_space=4` (was 8) in the pilot-point tests, keeping
     `n_pp >= 10` on the 2,100 m testing() domain.
   - Why/Impact: geometry-appropriate; sparser nets than production use.
   - Revisit: with any pilot-point algorithm change, re-run on validation().

4. **Large-multilayer render test shrunk.**
   - What: `test_interactive_plotting` large-multilayer cases now use
     nlay=4 / nper=6, mosaic `ncols=2`.
   - Why: rendering cost; the assertions still cover external frames, the
     1e30 dry-value mask, and NaN layers.
   - Impact: smaller frame payloads than the original stress shape.
   - Revisit: if an external-frames size threshold bug appears.

5. **Phase 4.5 second-tier splits deferred.**
   - What: splitting `grid/triangle.py` and `package_plotting.py` (the two
     remaining god-modules) postponed to Phase 8 prep.
   - Why: Phase 4's goal (import hygiene + the three worst modules) was met;
     recorded in the plan banner.
   - Revisit: Phase 8.

## Phase 5 — D8 / mf.riv / mf.evt (2026-07-17, branch `phase-5-package-api`)

> Entries 6, 13 and 14 below were **scheduled** rather than open-ended: plan
> §4.7 (package declaration consolidation) absorbed them — 4.7.5 built the
> `_ArealBuilder`/`EVTBuilder` (entry 6, **RESOLVED 2026-07-21**), 4.7.1 fixed
> `rch_spec` (13) and the budget lists, and §4.7.2 carries the mover list (14).
>
> **Correction (2026-07-18):** this banner previously said 4.7.1 would fix
> entry 14 (`rch` in `_MOVER_PACKAGES`). 4.7.1 shipped and did NOT — entry 14
> itself always recorded that as a deliberate out-of-scope deferral, so the
> banner overstated the promise, not the entry. Moved to 4.7.2, where `mover`
> becomes a declared descriptor capability and the fix lands with its permanent
> home. Flagged by the 2026-07-18 next-step survey.

6. **No `EVTBuilder` (no file-less whole-domain EVT form). — RESOLVED 2026-07-21 (plan §4.7.5).**
   - What: `mf.evt` had `()` / `.gpkg` / `.flopy` but no builder form like
     `mf.rch(context=, nper=, recharge=)` that self-selects top-active cells.
   - Resolution: `mf.evt(context=, nper=, rate=, depth=)` now builds a file-less
     whole-domain EVT package via `EVTBuilder`, a thin subclass of the new shared
     `_ArealBuilder` base (`modflow/mf6/areal.py`, extracted from `RCHBuilder`).
     The per-cell ET `surface` -- the one piece with no RCH analogue -- resolves
     through the shared `SurfaceResolver`/`CellSurfaceOffset` engine (default
     `model_top`/land surface), never reading `gdf_topbtm` columns directly.
     `rate`/`depth` reuse the RCH value-broadcast machinery verbatim. Deferred
     scope for the builder is recorded in entry 51.

7. **EVT `.gpkg` emits only the `nseg=1` record shape.**
   - What: `GeoPackageSource.evt` maps `surface/rate/depth` (4-field records).
     Segmented ET (`nseg>1` with `pxdp`/`petm`/`petm0` values) works only via
     `()`/`.flopy` with explicitly assembled records.
   - Why: per-segment GeoPackage column conventions would be speculative API.
   - Impact: GIS-driven segmented ET is not one call.
   - UPDATE 2026-07-18 (post-review): the limit is now enforced, not implied —
     `GeoPackageSource.evt` raises a clear `ValueError` on `nseg != 1` instead
     of letting FloPy fail later with an opaque record-shape error, and the
     `.gpkg` docstring no longer advertises `nseg` as a usable option.
   - Revisit: if a real model needs segmented ET from GIS.

8. **RIV input hover is the generic per-field hover, not the sketched
   multi-field one.**
   - What: plan §5.1 sketched `cell_input_hover("stage",
     extra_fields=("cond", "rbot"))`; shipped the same generic per-field
     hover every other list BC uses (`inputs.stage.map()` hovers stage, etc.).
   - Why: no existing package passes `extra_fields`; a bespoke riv hover would
     have been the first, diverging from siblings for marginal value.
   - Impact: one field per input-map hover (consistent with chd/ghb/drn/wel).
   - Revisit: if multi-field input hover is wanted, add it for ALL list BCs
     at once, not riv alone.

9. ~~**RIV/EVT are not in the canonical model.**~~
   **RESOLVED 2026-07-18.** Both are now contract packages (13 → 15), so they
   run inside the full integration matrix (groups, diffs, PEST, parallel
   splits). Placement was chosen to avoid double-counting ET: **EVT sits on the
   valley walls alongside RCH** (the standard recharge/ET pairing) and is
   **disjoint from UZF**, which covers the floor and — verified — runs with
   `simulate_et` auto-enabled but no `linear_gwet`/`square_gwet`, i.e. vadose-zone
   ET only. **RIV is the un-routed outlet river below the lake** (xn 0.92–0.99),
   carved out of the UZF footprint exactly as the lake/stream cells are.
   Pinned by `test_canonical_evt_and_uzf_footprints_stay_disjoint` (mutation-tested:
   growing UZF onto the walls fails it) and
   `test_canonical_riv_and_evt_carry_real_flux`.

10. **RIV/EVT are not PEST `parameterize` targets.**
    - What: `cal.parameterize` supports `chd`/`ghb.cond`/`drn.cond`/`wel`
      etc., but not `riv.cond`/`riv.stage`/`evt.rate`.
    - Why: the parameterization backlog is plan §5.8; not silently in 5.1/5.2
      scope — listed here so the asymmetry with sibling BCs is visible.
    - Revisit: fold into §5.8 when it runs.

## Phase 5 adversarial review follow-up (2026-07-18)

An ultracode review of the D8/riv/evt diff confirmed 9 findings; all were fixed
(see the review-fixes commit). The residue below is what was deliberately NOT
done, plus one assurance gap in the review itself.

12. ~~**The review that vetted this work was only ~50% complete.**~~
    **RESOLVED 2026-07-18** — review re-run to completion (18/18 agents,
    worktree-isolated; 4 confirmed of 14 claims). The **budget off-by-one is now
    FIXED** (commit below): it was worse than reported — `model.bud()` returned
    1-based nodes for **chd, riv, wel AND evt** (chd/wel predate riv/evt by
    years), while `group.bud('drn'|'ghb')` came back one cell low. Both halves
    were fixed together as required. **All four confirmed findings are now
    closed** (2026-07-18): the budget pair, plus `default_input` pinning (now
    table-driven over the registry's `cell_stress` set) and the
    `_PackageDiffNamespace` riv/evt properties. Nothing outstanding from this
    review; entry kept as the record of what it found. Delete at the next
    ledger tidy.

13. ~~**`rch_spec` keeps the fragile `maxbound` inference.**~~
    **RESOLVED 2026-07-18** — removed, matching `evt_spec` and every other list
    BC. Verified behavior-preserving: FloPy computes `MAXBOUND` at write time and
    the written `.rch` file is byte-identical (`MAXBOUND 3` for a 3-record
    period). `test_list_bc_specs_accept_every_native_flopy_input_shape` now
    covers all 7 list BCs against a bare list, a `None` period, and an empty
    dict, and asserts none of them pre-compute `maxbound`.

14. ~~**`_MOVER_PACKAGES` still lists `rch`, which has no mover support.**~~
    **RESOLVED 2026-07-18 in 4.7.2.** `rch` removed; the tuple is now exactly
    the seven packages whose FloPy constructor exposes `mover`
    (`drn ghb riv wel uzf sfr lak`). Unobservable — MF6 has no RCH mover, so the
    try-wrapped lookup always found nothing — and it upgrades a documented
    discrepancy into an equality invariant: `test_the_mover_list_equals_the_mover_capability`
    asserts the tuple equals `PackageCapabilities.mover` derived from live FloPy.
    The API snapshot caught the change and its diff is that one line.
    (A survey agent claimed `ghb` was wrong here too; checked against live
    FloPy — GHB **is** a valid MF6 mover provider, so only `rch` was wrong.)

15. **RIV/EVT input maps hover one field at a time (unchanged from entry 8).**
    - Confirmed still true post-review; the generic per-field hover is what all
      list BCs use. Not a regression, just not the plan's sketch.

16. **Mutation-testing residue is a real process hazard.**
    - What: the review's test-quality agents mutated `package_api.py` and
      `advanced.py` in the shared working tree and did not revert. One verifier
      then read the mutated tree and reported a false defect ("nseg forward
      deleted") as fact.
    - Why: workflow agents ran without `isolation: "worktree"`.
    - Impact: none persisted (mutations were uncommitted and reverted), but a
      review verdict was contaminated, and a careless `git restore` during
      cleanup briefly wiped uncommitted fixes.
    - Revisit: run review workflows that may mutate code with
      `isolation: "worktree"`, and always `git status` before trusting a
      subagent's file-based claim.

## Notebook session-state call (2026-07-17, `c578b53`)

11. **`RUN_SUITE` restored to `False` in `canonical_fast_tour.ipynb`.**
    - What: the user's session had flipped it to `True`; the commit keeps
      their three `.show()` edits but restores the documented `False` gate
      default and strips executed outputs (D2).
    - Why: `True` as the checked-in default would make every fresh notebook
      run launch the full pytest suite; judged to be leftover session state.
    - Impact: none if the judgment was right; user was flagged and can flip
      it back if `True` was intended.
    - Revisit: user says so.

## SFR profile view + view-layer convention (2026-07-18)

17. **The profile view has no `.map()`.**
    - What: `SfrProfileView` answers `get`/`summary`/`plot` but not `map`/`xs`/
      `mosaic`/`animate`, so it is not a full citizen of the spatial grammar.
    - Why: the profile is a merged multi-field table; `map()` would have to ask
      "which field?", and `sfr.results.q.map()` / `sfr.results.stage.map()`
      already answer that per field. Adding an ambiguous verb to satisfy
      symmetry would make the API worse, not better.
    - Impact: `docs/view_layer_conventions.md` says spatial verbs apply to
      "spatial nouns", which the profile is not. The asymmetry is documented
      rather than hidden.
    - Revisit: if a `profile.map(field=...)` need appears in practice.

18. **The loose field-level plot verbs — RESOLVED 2026-07-22 (plan §4.8).**
    - What was deferred: `sfr.results.q.plot_profile`, `sfr.results.stage.plot_profile`,
      and `lak.results.q.plot_budget` / `budget_summary` kept their loose spellings
      when only `long_profile`/`plot_long_profile` moved onto the view shape.
    - Resolution: each became a view-class noun answering `get`/`summary`/`plot` —
      `sfr.results.q.profile` and `sfr.results.stage.profile` (both
      `SfrReachProfileView`), and `lak.results.q.budget` (`LakBudgetView`). The old
      spellings are D12 warned aliases resolved via `__getattr__` (hidden from
      completion), preserving their exact old returns. Callers migrated (2 canonical
      notebooks, the API pamphlet, `test_colorscale_policy`); tests added to
      `test_sfr_profile_view.py` and `test_colorscale_policy.py`.
    - **Deviation from the plan's "keep old spellings as warned aliases":** the SFR
      field explorers' old *data* method was itself named `profile()`, which collides
      with the new `profile` noun property (a name cannot be both). So `profile`
      becomes the noun and the frame comes from `q.profile.get()`; there is no warned
      `profile()`-returns-a-frame alias. Safe because the old `profile()` df method had
      zero callers. `plot_profile`/`plot_budget`/`budget_summary` (no collision) are
      proper warned aliases.
    - **LAK budget figure stays matplotlib** (entry 52) — the `viz.Fig` conversion is
      a separate view-layer concern, deferred to keep this a pure verb migration with
      exact-return aliases.
    - Entry 17 (`profile` has no `.map()`) stands: spatial verbs live on the field
      explorers (`q.map()`), not on the derived-table nouns.

19. **`docs/refactor_review_report.md` still lists `long_profile(...)`.**
    - What: its API inventory (line ~851) names the now-deprecated spelling.
    - Why: that file is a **dated snapshot** ("Date: 2026-06-12") of a review,
      not a living reference; editing it would falsify the record of what the
      API looked like then.
    - Impact: a reader of that file could copy a deprecated name — but it warns
      on use and names its replacement.
    - Revisit: if the report is ever converted into a living document.

20. **The user's executed `canonical_fast_tour.ipynb` outputs were discarded.**
    - What: the working tree carried a 103,008-line output diff from a real
      run (13 of 22 cells). Editing cell 11 required touching the file, so the
      outputs were stripped per D2 and the notebook reset to HEAD before the
      source edit was reapplied.
    - Why: committing a source fix under a 103k-line output diff makes the
      change unreviewable, and D2 forbids tracked outputs regardless.
    - Impact: none that is not regenerable — the fast tour re-runs in ~1 min.
      The pre-strip file was backed up to the session scratchpad
      (`canonical_fast_tour.with_outputs.ipynb`) in case it is wanted.
    - Revisit: n/a.

## Calibration plot + canonical head observations (2026-07-18)

21. **The canonical model's "measured" heads are synthetic, not data.**
    - What: head targets read off a fitted analytic "observed water table"
      (`_OBSERVED_WATER_TABLE_COEF`), sampled at four named features plus a
      `monitor_1`…`monitor_8` valley-floor network.
    - Why: the model is synthetic — there are no measurements. Before this,
      every head target was `NaN`, which made `compare()`, `stats()`,
      `residuals()` and `calibration_plot()` all vacuous, and produced the
      reported blank calibration plot.
    - Impact: the surface is a seven-coefficient **trend**, not a copy of the
      simulated heads, so residuals stay real (a copy would score a perfect fit
      while teaching nothing). It reads as well calibrated — ME ≈ ±0.05 ft,
      RMSE 2.07 ft = 3.9% of range (testing) and 2.32 ft = 4.3% (validation).
    - **UPDATE 2026-07-18:** the earlier ~20%-of-range bias recorded here is
      GONE; entry 25 is resolved. What remains is only that the observations
      are synthetic at all, which is inherent to a synthetic model.
    - Revisit: if the model's physics changes, the coefficients must be refitted
      — `tests/test_canonical_head_observations.py` fails on both profiles if
      they drift, which is the intended tripwire.

22. **Production wells carry no measured value.**
    - What: `shallow_pumping` / `deep_pumping` remain observation *locations*
      (their simulated series and head-change signals still work) but have no
      "measured" head in any period.
    - Why: a regional water-level survey does not report a head measured inside
      a pumping well — it reads drawdown, and the deep well is screened in a
      different unit from the water table the survey describes. Including them
      planted +11 ft (testing) and +23 ft (validation) outliers that dominated
      every residual statistic.
    - Impact: 60 paired points from 10 measured wells; the two production wells
      contribute simulated series only. This is standard practice for a real
      calibration dataset, not a shortcut.
    - Revisit: if drawdown-aware synthetic observations are ever generated from
      a truth run (`canonical_calibration.py` already does that for PEST).

23. **The empty-cross-plot diagnosis warns rather than raising.**
    - What: `CalibrationPlot.from_obs_vs_sim` emits a `UserWarning` and draws an
      on-figure box when nothing pairs, instead of raising.
    - Why: an empty overlap is a legitimate intermediate state (targets declared
      before a run, a join still being built). Raising would break notebooks
      mid-tour; silence was the original defect.
    - Impact: a script that ignores warnings still gets a figure — but it now
      carries a visible annotation, so the failure cannot pass unnoticed in a
      notebook.
    - Revisit: if a strict mode is ever wanted, add `on_empty="raise"`.

24. **Only the `obs_vs_sim` path got the diagnosis.**
    - What: `from_timeseries`, `from_residuals_by_period`, and the `"heads"`
      MultiIndex path can still render empty without explanation.
    - Why: `obs_vs_sim` is the default path behind every
      `targets.<family>.calibration_plot()`, so it covers the reported defect
      and all five target families. The others take different frame shapes and
      would each need their own diagnosis.
    - Impact: an empty time-series plot is still silent.
    - Revisit: when the view-layer sweep in plan §4.8 lands.

25. ~~**OUTSTANDING: the fast tour must show a WELL-CALIBRATED model.**~~
    **RESOLVED 2026-07-18.** The scatter went from RMSE 6.6 ft (~20% of head
    range, every well biased the same direction) to **RMSE 2.07 ft = 3.9% of
    range with ME +0.05 ft** on the testing profile the fast tour uses, and
    2.32 ft / 4.3% / −0.03 ft on validation. Both approaches the user approved
    were used together:

    1. **A calibrated observation surface.** Refitting `regional` itself would
       have been circular — it drives CHD heads, GHB heads, initial conditions,
       SFR streambed tops and RIV stage, so changing it moves the very solution
       being fitted. Instead a *separate* trend surface was fitted offline by
       least squares against the simulated heads and frozen as literals
       (`_OBSERVED_WATER_TABLE_COEF`), so build time stays run-free. Its `tanh`
       term represents the mid-valley bedrock constriction, which a smooth
       polynomial cannot follow.
    2. **A wider network.** Eight `monitor_*` wells on the valley floor, sited
       by an a-priori rule fixed BEFORE any residual was inspected: spread along
       the valley, never inside the constriction band, and at least 1.2 cell
       widths clear of every lake/stream/river/drain/pond/pumping cell. Sited on
       normalized position, so the same physical locations are sampled at any
       resolution.

    Two findings made the difference and are worth remembering: the head field
    is nearly resolution-independent (1.6 ft RMS between profiles at matched
    normalized positions), so fixed fractional well positions transfer; and the
    residual statistics were dominated by the two production wells, which should
    never have carried survey values at all (entry 22).

    Guarded by `test_the_model_reads_as_well_calibrated` — RMSE < 5% of range
    AND |ME| < 0.75 ft, the second bound being what stops a surface that is
    precise but systematically offset, which is exactly how the first attempt
    failed. See also entries 21 and 22.

26. **The observation-surface coefficients have no regeneration script.**
    - What: `_OBSERVED_WATER_TABLE_COEF` in `canonical_example.py` is seven
      frozen float literals. They were produced by a least-squares fit against
      the canonical model's simulated heads in a **throwaway session**; that fit
      was never committed, so the numbers cannot be reproduced or updated from
      anything in the repo.
    - Why not done in the same pass: the fit needs a built *and run* canonical
      model on both profiles (~2 minutes), which makes it a `scripts/` tool
      rather than part of the build or the test suite. Landing the calibration
      result was the user's ask; the tooling around it was not.
    - Impact: **this is the fragile part of the well-calibrated demo.** If the
      model's physics changes — K field, boundaries, layer geometry, the
      regional gradient — the coefficients silently stop matching.
      `test_the_model_reads_as_well_calibrated` will *catch* the drift (RMSE
      < 5% of range, |ME| < 0.75 ft), so it fails loudly rather than quietly,
      but whoever hits that failure has to reconstruct the fitting procedure
      from the docstring and this entry instead of rerunning a script.
    - Shape of the fix: `scripts/fit_observation_surface.py`, the honest
      counterpart to `scripts/derive_api_snapshot.py` — build + run both
      profiles, fit the documented basis
      `[1, xn, xn^2, xn^3, cross_relief, tanh((xn - 0.45) / 0.08), layer]` over
      the valley-floor cells, print the coefficients for pasting, and report
      RMSE/ME per profile against the test's thresholds. A `--check` mode could
      verify the committed literals still fit, though that would cost two model
      runs and so belongs in the weekly lane, not per-push CI.
    - Revisit: next time the canonical model's physics is touched, or sooner if
      the calibration test ever fails.

## Package descriptor — plan 4.7.2 (2026-07-18)

27. ~~**The descriptor is written but not yet read.**~~
    **RESOLVED 2026-07-18 in 4.7.3** (5 commits). Six lists across five
    subsystems now derive from the descriptor: the results-diff pair, the
    model-diff pair, budget node basing, the suffix map, and the artifact sets.
    Each was verified identical to the literal it replaced before the literal
    was deleted.
    **Still hand-written by choice:** the `ModelGroup` `GroupPackageInputs`
    accessors and the `SimulationBase` package properties — code objects, not
    data, so they belong with 4.7.4.

28. **The `RETIRE WITH 4.7.3` tests become circular the moment they are wrong.**
    - What: the tier, suffix and node-basing assertions compare the descriptor
      against the hardcoded list. Once 4.7.3 deletes a list and derives it from
      the descriptor, that test compares the registry with itself and passes
      unconditionally.
    - Why: unavoidable for a "prove equality before replacing" step — the test
      is a migration scaffold, not a permanent invariant.
    - Impact: a stale scaffold would give false confidence, which is worse than
      no test. Each is marked `RETIRE WITH 4.7.3` in the file, the obligation is
      stated in the module docstring, in CLAUDE.md, and in the plan.
    - Revisit: delete each assertion in the same commit that deletes its list.

29. **`evt.record_fields` records the `nseg=1` prefix, not the full record.**
    - What: FloPy's dfn lists `surface rate depth pxdp petm petm0`; the
      descriptor carries `("surface", "rate", "depth")`.
    - Why: myflopy's GeoPackage path supports `nseg=1` only and raises on
      anything else (entry 7), so the segmented tail has no generator to feed.
    - Impact: 4.7.4 must not assume `record_fields` is the complete MF6 record
      for every package. The test asserts it IS a prefix of the dfn rather than
      hardcoding three names, so a reordering upstream still fails.
    - Revisit: with entry 7, if segmented ET from GIS is ever built.

30. **The registry still covers GWF packages only.**
    - What: 10 GWF packages are described; the GWT/GWE/PRT analogues are absent.
    - Why: plan 5.3 adds those packages and is explicitly gated behind 4.7.
      Describing packages that do not yet have a package-first surface would be
      speculative.
    - Impact: 5.3 must extend the descriptor rather than assume it is complete —
      which is precisely the discipline 4.7 exists to establish.
    - Revisit: 5.3.

31. **`gpkg_defaults` records parameter names, not MF6 record names.**
    - What: for DRN it is `{"elevation": "elevation", "conductance": "conductance"}`
      while `record_fields` is `("elev", "cond")`. The two vocabularies differ
      on purpose — GeoPackage columns are user-facing and verbose, MF6 record
      fields are terse.
    - Why: the resolver signature is the contract for `.gpkg` callers; renaming
      either side to match the other would be a breaking API change for no gain.
    - Impact: 4.7.4 needs an explicit mapping between the two when it collapses
      the resolvers, not a naive zip. Both are pinned against live sources, so
      the mapping can be derived rather than guessed.
    - Revisit: 4.7.4.

## 4.7.3 list deletion (2026-07-18)

32. **Three derived lists are unions, not pure derivations.**
    - What: `_PACKAGE_SUFFIX_TO_TYPE` keeps an explicit `_NON_REGISTRY_SUFFIXES`
      (`dis disu disv ic mvr npf oc sto`), and the artifact sets keep
      `_NON_REGISTRY_ARTIFACT_TYPES/_ORDER` (`ic npf mvr`).
    - Why: those are structural packages — discretization, initial conditions,
      node properties, output control, the mover — with no per-package
      behaviour (fields, results, capabilities) for the registry to describe.
      Forcing them in would make the registry a list of "packages that exist"
      rather than a descriptor of package behaviour.
    - Impact: two sources of truth remain for those lists, and the merge order
      matters — the non-registry half merges LAST, so a suffix or order claimed
      by both would silently shadow the registry. Both are now tested
      (collision, duplicate order, bracketing).
    - Revisit: if the registry ever grows a "structural" kind; not planned.

33. **Iteration order of the derived lists changed.**
    - What: `_CELL_BUDGET_PACKAGES`, `_MOVER_PACKAGES` and `_CONNECTION_PACKAGES`
      now follow the registry's declaration order rather than their previous
      hand-written order. Sets are identical.
    - Why: preserving the old order would have meant re-encoding it somewhere,
      which is the duplication being removed.
    - Impact: `MvrResultDiff.summary()` and the cell-budget diff build their
      rows in a different order. No test pins row order and the API snapshot
      stores these sorted, so nothing observable changed — but a consumer that
      relied on positional order rather than the `package` column would break.
    - Revisit: if row order ever becomes contractual, sort explicitly.

34. **Scaffold retirement is a manual obligation, not an enforced one.**
    - What: five `RETIRE WITH 4.7.3` tests were deleted as their lists went, and
      each retirement is documented in place. Nothing mechanically prevents a
      future refactor from leaving a scaffold behind once it goes circular.
    - Why: detecting "this assertion has become tautological" is not something
      pytest can do; it needs a human noticing that both sides now come from
      one source.
    - Impact: a stale scaffold passes unconditionally and reads as coverage.
      Mitigated by the in-place retirement notes explaining where the coverage
      moved, so the next reader can tell deletion from omission.
    - Revisit: 4.7.4, which will make several more assertions circular.

## 4.7.3 adversarial review follow-up (2026-07-18)

A review of the 4.7.2+4.7.3 diff raised 27 candidates; 6 survived independent
refutation, of which 2 were verified-clean negatives and 1 was cosmetic. The
three actionable ones are below — two were defects in the tests written the
same day, which is the useful part of the result.

35. **`_CELL_BUDGET_PACKAGES` has no consumer in `src/`.**
    - What: the constant is read only by `scripts/derive_api_snapshot.py` and
      two test modules. No production code uses it, so `tiers.results_diffable`
      currently drives nothing.
    - Why not removed: it is part of the pinned API surface, and
      `tests/test_model_group_symmetry.py` uses it as a genuine cross-check
      (hardcoded group accessors vs the registry). Deleting a snapshot-pinned
      constant is a bigger change than a review-fix pass should carry.
    - Impact: a reader may assume the flag gates per-cell results diffing when
      it gates nothing; the real diff path reaches packages another way.
    - Revisit: 4.7.4 — either wire it into the results-diff namespace or delete
      it and drop the tier flag.

36. **The review found a circular test in the commit that documented circular
    tests as the hazard.**
    - What: `test_the_results_diff_namespace_still_covers_every_flagged_package`
      asserted `set(_CELL_BUDGET_PACKAGES) == flagged` while
      `_CELL_BUDGET_PACKAGES` is derived from that same flag. Never able to
      fail. Written in the same pass as ledger 28, which warns about exactly
      this.
    - Why it happened: the retirement notes correctly identified WHICH tests
      went circular, then the replacement re-introduced the property by reading
      the derived constant instead of a consumer.
    - Impact: deleted, not rewritten — no honest consumer-side check exists
      while entry 35 stands. Also removed a tautological
      `assert orders == sorted(orders)` whose subject is built with `sorted()`.
    - Revisit: treat "does this assertion's other side trace back to the
      registry?" as a checklist item for every 4.7.4 replacement test.

37. ~~**Artifact capture/restore dropped MOVER from every list BC.**~~
    **FIXED 2026-07-18 in the same pass.** `components.py` hardcoded the mover
    set to `{uzf, lak, sfr}` in three places (capture, restore, dependency
    inference), so `mf.drn(spd, mover=True)` captured and restored with MOVER
    off — silent, and `_artifact_dependencies` reported no `mvr` dependency, so
    ordering validation could not catch it either. All three now read
    `capabilities.mover`. Guarded by three tests in
    `tests/test_package_artifact_round_trip.py`, verified to fail against each
    half of the old code independently. Note this survived 4.7.3(5/5), which
    edited the same file — deriving one list in a module does not find the
    others.

## 4.7.4 collapse (2026-07-18)

38. **The helper classes were NOT collapsed — the plan's premise was wrong.**
    - What: plan §4.7.4 called for `_ListBCPackage(descriptor)` replacing the 7
      `_CHDPackage`-style classes in `package_api.py`, on the stated grounds
      that the file is "3,038 lines, mostly repeated structure". Measured, the
      seven list-BC classes are 61% hand-written docstrings, 21% signatures and
      17% executable body.
    - Why declined (user decision, 2026-07-18): the collapse would delete ~747
      lines of package-specific prose in favour of generated text, and replace
      257 lines of explicit signatures with `**kwargs` — so
      `mf.riv.gpkg(stage=..., rbot=...)` would lose its named parameters,
      defaults and IDE completion. That trades the public surface's
      discoverability, which the user has twice said they care about, for ~150
      lines.
    - Impact: the seven classes remain, each ~175 lines of mostly documentation.
      Anyone adding a list BC still writes a helper class by hand — but it is a
      docstring-writing exercise, not logic duplication, since the bodies are
      now one-line delegations to the shared `*_spec`.
    - Revisit: only if the *bodies* diverge again. Prose duplication is not the
      thing 4.7 exists to remove.

39. **4.7.4 added lines rather than removing them.**
    - What: net +40 lines across `advanced.py` and `geopackage.py`, against a
      plan that billed this as the "biggest LOC win".
    - Why: the two shared bodies carry docstrings explaining the invariants they
      now own (why `maxbound` is never set; why record order is asserted). The
      deleted code was ~10-line dict literals with no explanation.
    - Impact: the win is real but categorical rather than numeric — a divergence
      class is gone. Do not use LOC as the success measure for 4.7.5/4.7.6.
    - Revisit: n/a.

40. **The record-order assertion is an `assert`, not a raise.**
    - What: `_bc_from_features` uses a bare `assert` to check field order
      against the descriptor, so it vanishes under `python -O`.
    - Why: it guards an internal call-site contract (the seven wrappers), not
      user input — the public signatures make the order unreachable from
      outside. A `ValueError` would imply callers can trigger it.
    - Impact: under `-O` a future wrapper written with swapped fields would
      build a wrong model silently. Mitigated by
      `test_record_field_order_is_checked_against_the_descriptor`, which also
      checks each resolver's declared signature order statically.
    - Revisit: if myflopy is ever run under `-O` in anger.

## 4.7.6 payoff test (2026-07-18)

41. **The payoff test asserts 9-of-16 automatic, not "every surface".**
    - What: plan §4.7.6 asked for a test that a synthetic descriptor "appears
      automatically in every surface". Measured, 9 surfaces are automatic and 7
      are not, so that assertion would have failed the day it was written.
    - Why: 3 of the 7 are deliberate (ledger 38 — prose and signatures), and 4
      are namespace properties deferred pending a registry-generated `.pyi`.
      Neither group is an oversight; both are recorded decisions.
    - Impact: the test asserts the true split and ratchets both ways, which is
      strictly more useful than a green light that would have required lying
      about scope. The plan text is corrected rather than the test weakened.
    - Revisit: each time a manual surface is closed.

42. **Four namespace properties remain hand-written, blocked on typing.**
    - What: `ModelPackages.<pkg>`, `SimulationBase.<pkg>`, `GroupPackages.<pkg>`
      and `_PackageDiffNamespace.<pkg>` must still be written per package.
    - Why: generating them at runtime degrades `diff.packages.riv` to `Any` for
      downstream users, because myflopy ships `py.typed` and static analysers
      cannot see runtime-generated properties. The plan's own "stays
      duplicated" list calls for emitting a checked-in `.pyi` from the registry
      with a test that the stub matches — that is the real unit of work, and it
      was never in 4.7.3's or 4.7.4's scope.
    - Impact: this is the largest remaining per-package cost, and it is what
      keeps 5.3 from being fully mechanical. Note `SimulationBase` accessors are
      also *incomplete* today (rch/uzf/sfr/lak exist; chd/drn/ghb/wel/riv/evt do
      not), so generating them would close a real gap as well as a DRY one.
    - Revisit: as its own step, before or alongside 5.3.

43. **The payoff probe runs in a subprocess.**
    - What: `tests/test_package_descriptor_payoff.py` shells out to a fresh
      interpreter rather than monkeypatching in-process.
    - Why: the derived lists are module-level constants evaluated at import, so
      an in-process injection would not propagate, and `importlib.reload` would
      leave other modules holding stale references — actively dangerous under
      `-n 10` with a session-scoped canonical model.
    - Impact: ~1.6 s per run and the probe body lives as a string, so it is not
      linted or type-checked. Kept short and assertion-free for that reason;
      all judgment lives in the test module.
    - Revisit: if the probe grows past ~60 lines, move it to a real file under
      `scripts/` so it can be linted.

## 4.7.7 namespace accessors (2026-07-18)

42. ~~**Four namespace properties remain hand-written, blocked on typing.**~~
    **RESOLVED 2026-07-18, by rejecting the premise.** The entry assumed the
    fix was "generate them + emit a registry-derived `.pyi`". Measuring first
    showed three of the four classes were already complete — only
    `SimulationBase` was missing anything, six accessors. And the `.pyi` half
    was worse than described: a stub file **replaces the whole module** for type
    checkers, so it would have meant hand-maintaining stubs for every other
    public name in four large modules.
    Shipped the six accessors plus
    `test_every_namespace_exposes_every_registry_package`, which asserts all
    four namespaces cover the registry. See entry 44 for what that leaves.

44. **The four namespace classes still write one property per package.**
    - What: adding a package still means four hand-written properties.
    - Why: generating them costs static typing (see 42), and the classes are not
      uniform anyway — `ModelPackages.uzf` returns `UzfPackageExplorer`,
      `GroupPackages.sfr` returns a parameterized
      `GroupResultsOnlyPackageAccessor[...]`, `_PackageDiffNamespace.lak`
      returns `_LakDiffNode`. Only the seven cell-stress packages share a shape;
      generating just those would leave the advanced ones hand-written anyway,
      for a third of the benefit and all of the typing cost.
    - Impact: four edits per new package — but a forgotten one now fails a test
      immediately rather than surfacing as an `AttributeError` months later,
      which is the failure that actually happened.
    - Revisit: if the advanced explorer types are ever unified, generation
      becomes worth reconsidering.

## Q sign convention — decided 2026-07-19, IMPLEMENTED 2026-07-20

45. **DONE: keep MF6's raw `q` signs; name the reference frame instead.**
    - **User decision, verbatim intent:** all `q` values must match MF6's own
      sign. Diverging from MF6 "is just asking for trouble". The ambiguity is
      resolved by naming the reference frame in the **column**, and by saying it
      in the docstrings — never by negating data.

      | column | reference | negative means |
      |---|---|---|
      | `q_gwf` | the GWF cell | discharge **out of the aquifer** (the feature *gains*) |
      | `q_lake` | the lake | the lake **loses** |

      Which frame applies follows from the FILE the record came from:
      `.cbc` cell records (SFR, DRN, RIV, CHD, WEL, GHB, RCH, EVT) are
      aquifer-referenced (`q_gwf`) and mutually consistent; package budget files
      (LAK's `GWF` record) are feature-referenced (`q_lake`). The trap is that
      LAK's record is *named* `GWF` but written from the lake's point of view —
      verified, the canonical perched lake is losing and reports ≈ −5888.
    - **SFR is `q_gwf`, not `q_sfr`.** The decision doc floated `q_sfr` with
      "negative = losing", but SFR's exchange is the aquifer's `.cbc` `SFR`
      record — *aquifer*-referenced, so negative = gaining, identical in frame to
      every list BC. Confirmed on the canonical model (41 gaining reaches, all
      `q < 0`). The user chose `q_gwf` (honest to the data) over `q_sfr`. Only
      LAK, which really is feature-referenced, gets a per-feature name (`q_lake`).
    - **Accessor stays `q`** (`model.packages.lak.results.q`) so the
      `results.<noun>.<verb>` grammar stays uniform (`docs/view_layer_conventions.md`).
      The explorer now carries a `result_name` (`"q"`, the accessor) distinct
      from `value_name` (`"q_gwf"`/`"q_lake"`, the emitted column); only the
      DataFrame column and its docstrings became frame-explicit.
    - **Scope: all user-facing q tables** (user's choice) — the result
      explorers *and* the group-compare + model-diff tables. Diff columns follow
      the raw name: `q_gwf` → `reference_q_gwf`, `q_gwf_diff`. Generic `.budget`
      term tables and the MVR mover-flow diff stay `q` (multi-term budget values
      with no gaining/losing ambiguity). The combined `surface_water` map draws
      the derived `exchange_intensity` field (normalized, positive = gaining).
    - **Plots/maps** keep gaining = blue, losing = red — reached by
      `_exchange_colorscale(frame)`, oriented from `ResultSpec.reference_frame`,
      not from pre-flipped data.
    - **What the implementation did (all landed):**
      1. Reverted the normalization in `build_sfr_budget_result_table` (commit
         `382530e`) and every downstream flip — commit 1 of 2.
      2. Replaced `ResultSpec.raw_gaining_sign` with `reference_frame`
         (`"gwf"`/`"feature"`); it describes the source file, so it can't go
         stale the way the old sign did.
      3. Renamed the columns via `value_name`, a TRUE rename at the read boundary
         so any missed consumer fails loudly; `result_name` keeps the accessor.
         `test_package_descriptor.py` now asserts value_name matches the frame.
      4. Re-pointed the three single-model maps and the two group maps.
      5. PEST is untouched: SFR flow obs take `.abs()`; DRN reads the separate
         legacy budget path. DRN resolves by naming (entry 46) — no run changes.
    - **Two bugs surfaced and fixed while implementing** (see entry 48).
    - Done: 2026-07-20. Removed from the outstanding queue.

46. **DRN's PEST series is sign-inconsistent with SFR/LAK (audit finding).**
    - What: the canonical springs discharge — the feature GAINS — and
      `targets.drn_flow.simulated_series()` reports −1274.9, while a gaining SFR
      reach reported +864 after the (now-to-be-reverted) normalization.
      `pest/forward_run.py` writes that series straight into the PEST simulated
      CSV, so observed discharge entered as a positive number yields a residual
      of ~−2550 with no warning.
    - Why not fixed yet: DRN was never in the surface-water convention, and its
      display layer and targets are currently consistent *with each other*.
      Normalizing only the targets would CREATE a plot-vs-PEST disagreement.
    - Resolution under entry 45: DRN is `.cbc`-sourced, so it is `q_gwf` —
      aquifer-referenced, negative = out of the aquifer. Naming makes it correct
      and consistent without flipping anything or breaking existing PEST runs.
      **Resolved 2026-07-20** with entry 45: DRN's exchange column is now `q_gwf`
      and its docstring states the frame; the PEST target series is unchanged
      (still MF6's raw sign), so no calibration run shifts.

47. **The 2026-07-19 sign audit (64 read paths, 63 locations) — TRIAGED 2026-07-21.**
    - The audit traced every SFR/LAK exchange read. Recovered its findings from
      the workflow journal and triaged all against post-rename code:
    - **Resolved by the revert (entry 45).** Every "raw-vs-normalized disagree"
      finding — `group.bud()`/`model.bud("sfr"|"lak").df`, the raw `budget_tables`
      reads, the `GroupBudget` split — is moot: nothing is normalized any more, so
      the raw escape-hatch reads and the `results.q` reads now carry the *same*
      MF6 sign (they differ only in column name: `q` vs `q_gwf`/`q_lake`).
    - **Resolved by the rename + its fixes.** The `gaining_sign` dual-meaning
      (→ `reference_frame`); all six map inversions (single + group × SFR/LAK/SW);
      the `mf6io_reference.md` LAK/SFR contradiction; the `_signed_exchange_colors`
      docstring; the circular colorscale test (rebuilt to assert on real maps).
    - **Four sites the audit flagged that the main rename pass MISSED, now fixed
      2026-07-21:** (a) `group/surface_water.py` combined map was still inverted
      (drew gaining red) — the only map I hadn't re-pointed; (b) the stale policy
      bullet in `view_layer_conventions.md` and (c) in `implementation_plan §2`
      that instructed the old blue-at-negative rule; (d) `canonical_02` cell 7 read
      `profile['q']` (now `q_gwf`). A new `test_every_group_exchange_map_puts_blue_
      on_gaining` covers the group tier that had no test.
    - **Surviving items — separable, pre-existing, non-sign (see entry 49).**
    - Done: 2026-07-21.

49. **Two pre-existing, non-sign SFR observation/plot issues (audit tail).**
    A 6-agent adversarial review (2026-07-21) verified what each really needed:
    49a warranted a small fix; 49b stays deferred.
    - **49a — FIXED (`observations/sfr.py` `SfrFlowTargets._package_flow_table`).**
      When in-channel FLOW-JA-FACE *routing* flow was unavailable, the fallback
      relabeled the stream-aquifer *exchange* (`sfr.results.q`, ~1e2 leakage) as
      `sim_flow` and fed it to PEST — two different physical quantities under one
      name, silently. The original ledger rationale ("only fires when FLOW-JA-FACE
      is unavailable, not a normal run") was **wrong in both directions**: a
      genuinely-absent record *crashes* on the unwrapped `.get()` at
      `sfr.py:329`, and the fallback actually fires on a **real** model when
      FLOW-JA-FACE is present but the SFR budget was saved on a different cadence
      than heads — `SFRBudget.get` swallows the length-mismatch `ValueError`
      (`budget.py:343-351`) and returns the raw record list (a non-DataFrame), so
      `not isinstance(flow, pd.DataFrame)` is True. Fix: the routing-missing branch
      now **warns and returns an empty frame** instead of substituting the
      exchange, so PEST/`compare` surface a *missing* value rather than fitting the
      wrong quantity. Two tests added (`test_mf6_pest.py`): one covers the primary
      FLOW-JA-FACE aggregation path (outbound-sum + TO-MVR), which had **no
      coverage** before; one asserts the refusal-to-substitute + warning. The
      mock-only branch (`not hasattr(model, "outputs")`) is left as-is on purpose —
      a real model always exposes `.outputs`, so it never runs in production; only
      test doubles reach it, and one existing test depends on it. Done: 2026-07-21.
    - **49b — DEFERRED (`budget.py` `SFRBudget.plot_flows`).** Names an SFR
      reach-to-reach series `riv_flows` (a misnomer — no RIV package is read); the
      unstated "method-shape concern" is that it is the `plot_<noun>()` anti-pattern
      the view-layer conventions forbid — it `fig.show()`s and returns `None`
      instead of a `viz.Fig`. Verified **zero callers** (not in any notebook, doc,
      `__all__`, or `api_snapshot.json`; reachable only via the legacy
      `model.outputs.sfr.bud` accessor). The misnomer is a *local variable* — the
      plotted trace is `stream {i}`, so nothing mislabeled reaches output; it is
      cosmetic. Its real correctness gaps (fragile `get('flow')` substring term
      selection; no per-reach groupby / TO-MVR aggregation; unnormalized single
      `(kstp,kper)` selector) only bite *if the dead method is called*, which
      nothing does. A bare delete was rejected: the modern routing noun
      `model.packages.sfr.budget.flow_ja_face` has **no `.plot()`** equivalent, so
      deletion is a small capability regression that belongs with adding the
      replacement verb — scope, not triage.
    - Revisit 49b: when SFR observations / the SFR budget view layer are next
      revised, retire `plot_flows` into `flow_ja_face.plot()` (returning `viz.Fig`)
      and carry over its three correctness gaps, rather than patching a dead method.

48. **Two bugs surfaced while implementing the reference-frame rename (entry 45).**
    - **Group LAK map was inverted** (pre-existing). The single-model LAK map was
      fixed in `c55ea23` to put gaining on blue, but that fix never reached the
      group layer: `project/group/lak.py` kept the raw blue-at-negative scale, so
      a gaining lake drew red on `group.packages.lak.results.q.map(...)` while the
      single-model map drew it blue. Fixed by orienting both group maps from the
      declared frame (`_exchange_colorscale("feature")` for LAK,
      `"gwf"` for SFR); verified the two maps now agree.
    - **xdist masked a silent-skip.** After the rename, the SFR profile plot's
      guard `if "q" in frame.columns` silently dropped the exchange trace (the
      column is now `q_gwf`), and `build_sfr_q_map_payload` read a missing `"q"`.
      The full suite under `pytest -n 10` reported **green** while the same tests
      **failed serially** (`-n0`) — a session-fixture/xdist interaction hid real
      failures behind a passing count. Both were caught only by running serially
      and by a direct scripted check of `results.q.map()`.
    - **Standing lesson (added to the q-sign memory):** a green xdist count is not
      proof. Verify a sign/column change **serially** and by observing the real
      call site (`.map()`, `.plot()`), not just via a hand-built frame — the same
      failure mode that let the 2026-07-19 map inversion ship through a green test.

50. **Flat-bottom `bathy` lake: WARN, not raise (judgment call).**
    - `LAKBuilder`'s `bathy` mode connects a cell to a neighbor only where its
      lake bottom is BELOW the neighbor's (an exposed step). A flat `lake_bottom`
      has no steps, so a multi-cell bathy lake produces ZERO horizontal
      connections — a lake with no lateral aquifer exchange (an earlier config
      dropped from 39 horizontal connections to 0 this way).
    - **Chose a warning over a hard error.** The geometry is internally
      consistent, not malformed, and the codebase already had an intentional test
      (`test_bathy_scalar_bottom_makes_no_exposed_steps`) plus fixtures relying on
      it building. Raising would break internally-consistent configs; a
      `UserWarning` surfaces the likely mistake and names the three real choices
      (real bathymetry / `rectangular` / `only_vertical`) while letting it build.
    - **Trade recorded:** a warning can be missed (a hard error cannot), so a
      genuinely-degenerate lake can still reach a run. Accepted because the config
      is valid MF6 and the user's guidance was "internally consistent — needs a
      guard, not a fix." Revisit only if silent degenerate lakes recur in
      practice.
    - Two incidental flat-`bathy` test fixtures were cleaned up rather than left
      warning: `_wide_lake_builder` became `rectangular` (its docstring promised
      per-cell sidewall building, which flat bathy never did), and
      `test_replace_lak_on_loaded_run` (a replace-mechanism test, not a geometry
      one) filters the expected warning.

52. **LAK budget bar (`lak.results.q.budget.plot`) is still matplotlib (plan §4.8, 2026-07-22).**
    - What: the migrated `budget.plot()` renders the connection-type bar chart with
      matplotlib (`mpl_axes`, `ax.bar`, hex colors), carried over verbatim from the
      old `plot_budget`. The house rule is "figures are always `viz.Fig`" with colors
      from the policy helper, not hex literals.
    - Why deferred: §4.8 was scoped as a pure *verb* migration (loose `plot_foo` →
      `noun.plot`) with D12 aliases preserving *exact* returns. Converting the backend
      mpl→plotly changes the return type and pulls in the discrete signed-color policy
      plumbing (blue = the lake gains / inflow, red = loses / outflow) — a separable
      concern that would have bloated the verb sweep and broken the "exact return"
      alias contract.
    - Impact: `budget.plot()` and its `plot_budget` alias both return a matplotlib
      figure; the SFR profile nouns already return `viz.Fig`. One node out of uniform.
    - Revisit: when the discrete-signed-bar color policy is factored out (the same
      helper the SFR signed-exchange bars use), convert the LAK budget bar to `viz.Fig`.

51. **`_ArealBuilder`/`EVTBuilder` deferred scope (plan §4.7.5, 2026-07-21).**
    The 4.7.5 extraction was kept minimal by explicit user decision ("do the
    minimal but record the quirks as future things on the ledger"). Deferred:
    - **Segmented ET (`nseg>1`)** — the builder emits only the single-segment
      `(cellid, surface, rate, depth)` record; `pxdp`/`petm`/`petm0`/
      `surf_rate_specified` and `auxiliary` are direct/`.flopy`-form only (mirrors
      entry 7's `.gpkg` `nseg=1` limit). RCH has no segmentation analogue.
    - **`cell_top` surface for a cell top-active in a layer below 0** — needs
      per-layer `layer_N_top` surface columns (the same requirement the `.gpkg`
      path already has). The builder default is `model_top` (land surface), which
      resolves for every column regardless of top-active layer, so this only bites
      a caller who explicitly passes `CellSurfaceOffset("cell_top")` on a deep
      top-active cell.
    - **`surface_only` cells keyword is a `top_active` alias** — accepted by
      `__post_init__` but `_ArealBuilder.selected_cells` takes the first active
      layer for both (inherited verbatim from `RCHBuilder`; `UZFBuilder`
      distinguishes `surface_only`). Not "fixed" to avoid changing shipped RCH.
    - **No-domain `top_active` == `all_active`** — with `context.domain is None`
      both selectors return every cell on `layer` (inherited RCH behavior).
    - **`_ArealBuilder` is not yet consumed by `UZFBuilder`**, a third
      independent reimplementation of the same idomain selection; it diverges on
      `surface_only` semantics, so folding it in is its own follow-up.
    Revisit any of these on real user demand.

53. **`mf.dis`/`mf.disu` build+run but get no myflopy viz (plan §5.5, 2026-07-23).**
    - What: the new `mf.dis` (structured) and `mf.disu` (fully unstructured)
      passthroughs produce valid, runnable MF6 models, but the choropleth-map /
      cross-section / animation view layer targets **DISV/Voronoi** meshes only. A
      DIS or raw-DISU model has no `.map()`/`.xs()`/`.animate()` rendering from
      myflopy — a caller falls back to FloPy's own plotting, or uses DISV.
    - Why: myflopy is Voronoi/DISV-first by design; the whole spatial viz stack
      (`VoronoiGridPlus`, cell-polygon choropleths) is built on the vertex mesh.
      Wiring structured/DISU rendering is a large, separate effort with little
      demand for the Voronoi-first use cases §5.5 targets. The passthroughs exist so
      structured/externally-supplied grids can be *built* with the same package-first
      grammar, not to make them first-class in the viz layer.
    - Also deferred (mirrors `mf.disv`/`ic`/`oc`): **no `package_registry.py` entry
      and no results tier** for `dis`/`disu` — they are grid packages, not BCs. And
      **PRT is not dispatched** by any grid helper (`mf.prt` builds its own dis/disv
      internally; MF6 has no `ModflowPrtdisu`) — the shared dispatch raises a clear
      `ValueError` for `prt6` rather than growing a fourth column.
    - Revisit: if structured-grid choropleth/xs rendering is ever requested, it
      belongs with the Phase-6 generic dependent-variable surface work, not here.

54. **YAML serializer ships without TOML (plan §5.6, 2026-07-23).**
    - What: `SimulationSpec.to_yaml`/`from_yaml` + `Project.add_simulation_from_yaml`
      (in `specs_io.py`) support **YAML only**. The plan listed "optional TOML via
      stdlib `tomllib` + `tomli-w` extra"; that was deferred (user-approved).
    - Why deferred: two concrete frictions for a rarely-used format. (a) **TOML has
      no null type** — `to_dict()` legitimately emits `None` (e.g. optional builder/
      workspace fields), which `tomli-w` cannot write; supporting TOML means stripping
      `None` before dump and tolerating the asymmetry on read. (b) **`tomllib` is
      stdlib only on 3.11+** while the package targets `requires-python >=3.10`, so
      the read path would need a `tomli` backport dependency and a version gate. YAML
      has neither problem (PyYAML is one small pure-Python core dep, safe mode, and
      represents null natively).
    - Impact: no `to_toml`/`from_toml`; YAML is the one file format. The dict form
      (`to_dict`/`from_dict`) is unaffected and remains the substrate any future TOML
      wrapper would sit on.
    - Revisit: add TOML if a user wants it AND the floor moves to 3.11 (or a backport
      dep is acceptable); the null-stripping is the only real design work left.

55. **List-BC builders serialize via a bespoke partial encoding (plan §5.6A, 2026-07-23).**
    - What: the 7 list BCs (chd/ghb/drn/riv/wel/rch/evt) build through
      `functools.partial(_build_named, <FloPy class>)` (`advanced.py:_factory`), which
      is picklable but not importable, so `_callable_ref` used to reject them and they
      could not `to_dict()`. `_callable_ref`/`_resolve_callable` now special-case
      `functools.partial`, emitting `{"partial": <func ref>, "args": [{"$callable":
      <class ref>}], "keywords": {...}}`.
    - Trade: this couples the serializer to the *shape* of a partial-based builder. It
      is fully general for any `partial` over importable pieces, but the only builder
      that uses it today is `_build_named`. An alternative — replacing the partial with
      a single module-level `build_list_bc` function carrying the class ref in options —
      was **not** taken: it would churn the artifact/build path and the options shape
      (and the seven `*_spec` factories' public signatures) for no user-visible gain.
    - Revisit: if a second, differently-shaped partial builder ever appears, reconsider
      whether the general partial encoding or a module-level-function refactor is cleaner.

56. **GWT/GWE results tier shipped reader+maps only (plan §6.0/6.1/6.2, 2026-07-24).**
    - What: the user scoped this pass to the **keystone + readable `model.conc`/
      `model.temp`** (grammar `get/summary/array/map/xs/mosaic/animate` + field hover +
      `'earth'` colorscale, via the full value-kind `Choro` hook). Four 6.1/6.2
      sub-items were **deliberately deferred**: (3) GWT/GWE **budget views** (`_get_budget_reader`
      is generic but `ModelView` budget plumbing is GWF-shaped)
      — **item 3 CLOSED 2026-07-27**: `model.budget.<term>` ships (ledger 97), and the
      stated blocker turned out to be wrong twice over. `_get_budget_reader` was generic
      *and* the `ModelView` plumbing was not the obstacle; the real obstacle was six
      correctness bugs in the shared budget-reading path (ledger 91–96), which had to be
      fixed before any view could be built on it; (4) **`GroupConc`/
      `GroupTemp`** member + Δ diff maps (near-mechanical `GroupHeads` clone with
      `elev`→`conc`, blocked on an `all_conc`/`all_temp` group table)
      — **item 4 CLOSED 2026-07-27** (ledger 100). The blocker was real and was
      built. "Near-mechanical clone" was NOT taken: §6.2 says to fix the
      abstraction rather than copy-paste, so `GroupHeads` moved onto a shared
      `_GroupFieldView` and all three kinds are now one implementation. The
      estimate also missed that `xs` was broken on the transport readers
      entirely (ledger 99); (5) **`ConcTargets`/
      `TempTargets`** + end-to-end transport **calibration** (~8 PEST wiring points + a
      new concentration forward-run post-processor + transport `parameterize` targets);
      (6) the **canonical `transport=True` fixture** (large — the canonical runs on the
      single-model `SimulationBase`, so a coupled GWT needs the multi-model spec path).
    - Why: each is additive on top of a readable reader, and the canonical fixture is a
      sub-project of its own; landing the reader first makes transport models *readable*
      (the stated goal) in a bounded change and lets the rest follow incrementally.
    - Two smaller judgment calls inside the shipped work:
      - **`self.gwf` is now the kind-neutral flopy handle** on a transport `ModelView`
        (a slight naming lie), NOT renamed to `flopy_model` everywhere — `model.gwf` is
        read in ~20 pervasive sites (components, xsections, inputs, choros, ...) whose
        calls (`modelgrid`/`output`/`get_package`) are kind-agnostic; a full rename was
        judged higher-risk than the naming imprecision. `.model_type` is the honest kind.
      - **`ConcResults`/`TempResults` keep the internal store column and `.ucn` default
        suffix**; the reader reads `<model>.ucn` (or an explicit path). A GWE model that
        names its temperature output something else must pass the path. Documented on the
        reader; revisit if a different convention shows up.
    - Revisit: build the four deferred sub-items when transport calibration or grouped
      transport comparison is actually needed; the reader/hook/colorscale substrate is
      in place for all of them.

57. **PRT derived maps leave unreached cells blank, not zero (plan §6.3B, 2026-07-25).**
    - What: `travel_time.map()` / `endpoints.map()` / `capture.map()` fill cells no
      particle reached with **NaN** (`fill_value=float("nan")`), so they are drawn as
      gaps. The alternative — `fill_value=0` — is available through `map(fill_value=0)`.
    - Why: on an `'earth'` scale a zero-filled grid renders every untouched cell as a
      solid low-end color, which reads as "0 days travel time" / "0 particles measured
      here" rather than "no particle went here". A capture-zone figure is about *where*
      particles are, so absence should be absence. For counts the zero reading is at
      least arguable, but one rule across the three maps beats a per-noun exception.
    - Trade: a NaN cell has no hover, so there is no way to hover-confirm "nothing here";
      the choropleth simply shows the grid outline. Accepted.
    - Revisit: if a user wants an explicit zero-count background, `fill_value=0` already
      does it — promote it to the default only if that turns out to be the common ask.

58. **`capture` raises on an ungrouped run instead of falling back (plan §6.3B, 2026-07-25).**
    - What: `results.capture.get()/map()` needs PRP boundnames. When the run has none,
      it raises a `ValueError` naming both fixes (`PRTReleasePoints.from_cells(...,
      group='west_wells')` and `by='release_point'`) rather than silently grouping by
      `irpt`.
    - Why: the silent fallback produces a *plausible* figure — one panel per particle —
      that answers a different question than the one asked, and on a run with hundreds
      of release points it would produce hundreds of mosaic panels. The user's locked
      decision for §6.3 was "real release groups", not the irpt fallback; keeping the
      fallback reachable but explicit honors both.
    - Trade: an extra keyword for anyone who genuinely wants per-release-point capture.
    - Revisit: not expected to change.

59. **PRT views replace the series verb rather than implementing it (plan §6.3B, 2026-07-25).**
    - What: `SpatialView.plot()` draws a value-by-stress-period series. PRT results have
      no period axis, so the three PRT nouns override `plot()` with distribution figures
      (cumulative arrival curve; particles-per-cell and particles-per-group bars) and
      make `_series_table()` raise a message that says so. `mosaic(kind="plot")`
      therefore raises too. For the same reason `animate()` raises rather than emitting
      a one-frame animation labelled "Period 0", and `map()` rejects `per=` instead of
      swallowing it into `cor()`'s `**kwargs` — the two ways an ignored period selector
      could have quietly claimed a period the values do not belong to.
    - Why: the grammar's promise is "every noun answers the same verbs", and `plot()`
      returning a `viz.Fig` is kept — what changes is what the x-axis *means*, because a
      stress-period axis does not exist for a time-integrated tracking run. A raise with
      the reason beats a `KeyError: 'per'` from inside a groupby.
    - Trade: `plot()` is not uniform in meaning across nouns (period series elsewhere,
      distribution here). Documented on each method.
    - Revisit: if transient PRT tracking (`track_times`) grows a per-period story, a real
      series verb could be added alongside.

60. **`logscale` now applies to custom-`zs` choropleths, blanking non-positive values
    (plan §6.3B, 2026-07-25).**
    - What: `Choro.zs` returned `custom_zs` verbatim *before* reaching its
      `if self.logscale: zs = np.log10(zs)` branch, so `logscale=True` was silently a
      no-op on every `type="custom"` map — including the PRT travel-time map, whose
      docstring advertises it. Found by an adversarial review of this change (it was the
      one finding of 32 that survived verification). Fixed at the root, in `Choro.zs`.
    - Judgment call inside the fix: non-positive values become **NaN**, not `-inf`.
      A travel time of 0 is ordinary (a particle terminating at release), and `-inf`
      would drag the shared color range to negative infinity and blank the entire map.
      NaN draws as a gap, matching ledger 57's "absence is absence".
    - Trade: the pre-existing non-custom branch (heads/K/rch) still lets `log10(0)` through
      as `-inf`. It was left alone deliberately — changing how existing head and K maps
      render is not this change's business, and no in-repo caller pairs `custom_zs` with
      `logscale`, so the fix has no blast radius today.
    - Revisit: fold the non-custom branch onto the same `_logscaled` helper when a map
      there is actually seen misbehaving on zeros.

61. **`PRTRunResults.pathlines` became a view; the raw frame moved to `track_records`
    (plan §6.3C, 2026-07-25).**
    - What: adding the pathline map meant `pathlines` had to be a noun answering
      `get/summary/map/plot/mosaic`, but the name was already a `cached_property`
      returning the raw track CSV (~20 call sites across 6 files, incl. canonical notebook 03's
      `prt_results.pathlines.empty`). The raw frame is now `results.track_records`;
      `pathlines.get()` returns the normalized **superset** (every raw column plus
      `cell`/`layer`/`travel_time`/`release_group`/`particle`).
    - Judgment call: a **clean break**, chosen over a DataFrame-proxy shim that would
      have forwarded `.empty`/`len()`/`[...]` to `.get()`. The proxy would have kept
      unmigrated notebooks running, but it blurs "a noun is an object" and would silently
      absorb typos as frame lookups. `cached_property` also cannot carry a deprecation
      shim the way a renamed method can (D12 is about `__getattr__` aliases).
    - Trade: user notebooks outside the repo that treat `results.pathlines` as a frame
      raise `AttributeError` until they add `.get()`. Every in-repo site (src, tests,
      example script, notebook 03 cell 3) was migrated in the same commit.
    - Revisit: not planned. If the break proves noisy in practice, a `__getattr__` on the
      view that raises a *pointed* message (rather than forwarding) is the middle ground.

62. **`viz.mosaic` now copies map overlays — a pre-existing silent drop (plan §6.3C,
    2026-07-25).**
    - What: the composer copied exactly one trace per map panel
      (`cell_traces.append([panel.get_choropleth()])`), so contour lines and location
      markers have been missing from **every** mosaic since they were added — the
      figure looked complete, just without them. Rather than special-casing pathlines,
      `Choro` gained `add_overlay()`/`overlay_traces()` and `add_contours`/`add_locs`
      were refactored into non-mutating `_contour_traces()`/`_locs_traces()` builders,
      so one accessor answers for all three sources (contours, locs, registered
      overlays); mosaic copies `[get_choropleth(), *overlay_traces()]`.
    - Trade: an existing mosaic of contoured maps now *renders its contours*, which
      changes figures a user may already have. Taken deliberately (locked decision): it
      is the same family of defect as ledger 60 — a documented feature silently doing
      nothing — and the shared-color-scale logic is now guarded to the cell trace
      (`z is not None`), since `Scattermap` has neither `z` nor a top-level `coloraxis`.
    - Small fix folded in: `_locs_traces` skips geometry that is neither Polygon nor
      Point. The old loop left `coords`/`mode` loop-local, so such a feature raised
      `NameError` on the first row and — worse — on any later row silently redrew the
      *previous* feature's geometry. Mosaic now asks every panel for its overlays, so
      that path is reachable from more places.
    - Revisit: three sibling paths still drop overlays and are now asymmetric with
      mosaic — `SpatialView`'s **animation** and `Choro.ani` rebuild frames from the
      cell trace alone, and the **matplotlib** mosaic reads only `choro.zs`.
      `add_hillshade` is a layout image, not a trace, so it is not carried either;
      left alone — layout images are per-subplot-anchored and need their own
      coordinate handling. Also pre-existing and now inherited by overlays: calling
      `.choropleth`/`.plot()` twice appends every trace again (it never clears
      `self.fig.data`), so a re-plotted map doubles its traces.

63. **The pathline map caps particles at 250 rather than drawing every trace (plan
    §6.3C, 2026-07-25).**
    - What: one `Scattermap` per particle is what makes per-particle hover and legend
      toggling work, but a release of thousands would produce a figure that is unusable
      to render. `map(max_particles=250)` draws a sample; `None` draws all.
    - Judgment call: the sample is **stratified by release group** (round-robin quota,
      then evenly spaced within each group), not "the first N" or a flat even sample
      across sorted particles — either of those can drop an entire capture zone and
      change what the figure appears to say. The cap is announced twice: in the figure
      title (`"(2 of 4 particles)"`) and a `UserWarning`.
    - Trade: a capped figure is not the whole run, and per-particle hover on a sampled
      figure can mislead someone who does not read the title. The no-group-dropped
      promise holds only while `max_particles` is at least the number of groups; below
      that the round-robin keeps the first N groups in sorted order. The alternative — one
      NaN-separated trace per group — renders far more particles but loses per-particle
      legend toggling and complicates the hover, so it was not taken.
    - Revisit: if real users routinely release thousands, add `render="grouped"` for the
      NaN-separated fast path alongside the per-particle default.

64. **Canonical notebook 03 lost its stream/lake outlines when it moved onto the library
    figure (plan §6.3C, 2026-07-25).**
    - What: cell 5 hand-rolled ~20 lines of matplotlib — a head choropleth, stream and
      lake cell outlines, and pathlines — exactly the drift `docs/view_layer_conventions.md`
      exists to prevent. It is now `prt_results.pathlines.map(...).plot()`. The head field
      and the paths are the same or better (interactive, hover, house styling); the
      colored region outlines are gone, because `outline_regions=` exists only on the
      matplotlib branch (`plot_mpl`).
    - Trade: the notebook's "what to look for" prose points at the stream corridor and
      terminal lake, which are now read from the head field rather than outlined.
    - Revisit: region outlines are a natural first user of `Choro.add_overlay` (62) —
      a `regions=` argument on the plotly map would restore them everywhere, not just
      in this notebook.

65. **`viz.category_colors` remembers name → color in a process-global memo, and a
    7-color palette can still collide (plan §6.3C, 2026-07-25).**
    - What: the helper exists so a release group keeps one color across its map, its
      arrival curve, and its capture bars. That requires memory outside any one figure,
      so the mapping lives in a module-level dict for the life of the process.
    - Judgment call: a name new to a call takes the least-used color that **none of the
      other names in the same call** already hold. The first implementation indexed
      `PALETTE.categorical` by the global counter, which an adversarial review showed
      hands two groups the *same* color once the palette wraps (register a group, then
      seven others, then its sibling). Per-call collision avoidance fixes the case that
      matters — the categories of one figure.
    - Residual, accepted: two names first seen in **separate** calls can still collide
      once more than `len(PALETTE.categorical)` names exist. Seven colors cannot promise
      otherwise, and re-assigning an old name to fix it would break the stability the
      helper exists for. Colors are therefore also not reproducible across processes:
      they depend on what was registered first.
    - Mitigation: `memoize=False` (used by `pathlines.map(color="particle")`) colors one
      figure without touching the memo — otherwise a single 250-particle map would
      register 250 names and shift every later figure's group colors.
    - Testing: the memo is global, so any test asserting a specific hex is order- and
      worker-dependent. `tests/conftest.py::isolated_category_colors` saves, clears, and
      restores it; use it rather than asserting colors against whatever ran first.
    - Revisit: if collisions show up in practice, the options are a longer palette or
      scoping the memo to a figure/session object rather than the process.

66. **PRT view tests split between one real MF6 run and hand-written track CSVs (plan
    §6.3C, 2026-07-25).**
    - What: `tests/test_prt_maps.py` runs real two-cell GWF+PRT simulations for the
      fixtures whose *fidelity to MF6* is the point (grouped, un-grouped, and mixed
      named/un-named releases), and feeds hand-written track CSVs through
      `open_prt_run` for the cases MF6 cannot easily be made to produce on a two-cell
      grid: multi-layer endpoints, a pooled-vs-per-layer statistic that differs, seven
      particles arranged so a cap must stratify, rows out of time order, and an empty
      run. Roughly a third of the file is synthetic.
    - Why: a two-cell, one-layer model cannot exhibit multi-layer or many-particle
      behavior, and forcing it would need a fixture heavy enough to slow the suite.
    - Risk, and it has already bitten: a synthetic CSV encodes an *assumption* about
      what MF6 writes. The first version of the mixed named/un-named test asserted that
      MF6 leaves `name` blank for un-named points. Checking against a real run showed
      it **synthesizes** `PRP000000002` instead; the test now runs MF6 and the blank
      case is kept only as an explicitly-labelled external-CSV defense. `icell` basing
      is pinned the same way (real run + a unit).
    - Rule of thumb going forward: anything asserting *what MF6 produces* must run MF6;
      synthetic CSVs are for exercising *our* transformations on shapes the tiny model
      cannot reach, and their docstrings must not claim MF6 provenance.

## 6.4A — IES field-map policy, hover, and plumbing (2026-07-25)

67. **The per-stat field-map color policy is module-private to `ies.py`, not a shared
    helper (plan §6.4 item 1).**
    - What: `_field_map_policy(stat, values)` lives in `src/myflopy/modflow/mf6/pest/ies.py`
      and returns `Choro` kwargs, rather than becoming a public `viz` or
      `package_registry` surface.
    - Why: `package_registry` is keyed by package/field **name**, and a *statistic* of a
      captured array is not a package field — the same captured K array is a magnitude
      as a `mean` and a ratio as a `change`, so a stat axis would be a second keyspace
      for one caller. `viz.py` holds no continuous-colorscale policy today. The house
      precedent for a registry-free derived-map scale is a constant in the consumer
      module (`prt_maps.PRT_COLORSCALE`).
    - Cost, accepted: `pest.ies` moves import layer 1 → 2, because the policy
      module-level-imports `package_plotting` for the two color helpers. Verified
      acyclic (nothing imports `pest.ies`); `prt_maps` already sits at layer 2 for the
      same reason.
    - Revisit: when 6.4B/6.4C give it a caller outside `IesResults`, promote it to
      `package_plotting` beside the color helpers. Promotion is mechanical; choosing a
      public signature before a second caller exists is not.

68. **`zmid` is emitted but is not what centers the `change` map.**
    - What: the policy returns `zmid=0.0` *and* symmetric `zmin`/`zmax` in log space.
      The symmetric limits are what actually center it.
    - Why: `zmid` reaches only the Plotly trace; `Choro.plot_mpl` reads `self._zmin`/
      `_zmax` and never consults the trace kwargs, so `zmid` alone would center the
      interactive map and leave the static one autoscaled — the two backends
      disagreeing about where "no change" sits.
    - Accepted: a redundant key on the plotly side, in exchange for one code path that
      is correct on both.

69. **Diverging parameter-field scales ship as explicit STOPS, never as the name
    `'RdBu'`.**
    - What: `_red_white_blue_diverging_colorscale()` (new, `package_plotting.py`) is
      derived by reversing its blue-white-red sibling.
    - Why: `Choro`'s plotly→matplotlib table maps `'rdbu'` to `'RdBu_r'`, the REVERSED
      colormap, so a diverging scale passed by name renders **mirrored** between
      `plot()` and `plot_mpl()` — red meaning "reduced" in one and "increased" in the
      other. Stops survive both backends. Pinned by
      `test_field_map_policy_survives_the_real_choropleth_front_door`.
    - Not fixed here: the `_PLOTLY_TO_MPL_CMAP['rdbu']` entry itself is left alone.
      `'RdBu'` is the registry's declared scale for signed `q` on ten packages, so
      flipping the table would silently re-color every one of those static maps. That
      is its own change with its own review — see 70.

70. **The signed-`q` maps still render mirrored between backends (pre-existing, NOT
    fixed in 6.4A).**
    - What: `package_registry` declares `colorscale="RdBu"` for `q` on rch/chd/drn/ghb/
      riv/wel/evt/sfr/lak. Because of the table entry in 69, `results.q.map(backend='mpl')`
      puts red where `plot()` puts blue.
    - Why not fixed: out of 6.4A's scope, and the fix is a judgment call about which
      orientation is *correct* for signed q (which is per-package — see the q-sign
      convention section) rather than a mechanical flip. 6.4A avoids the trap for its
      own maps by shipping stops.
    - Cost: anyone reading a static signed-q map today reads it backwards relative to
      the interactive one.
    - **Correction (2026-07-27):** this entry understated the blast radius. The
      mirroring was NOT confined to the registry's per-package `q` scales — the
      **grouped difference maps** were affected too, via
      `get_default_group_compare_colorscale()` returning the same bare name to six
      call sites (`group/spatial.py`, `inputs.py`, `results.py`, `sfr.py`, `uzf.py`,
      `lak.py`). That half is **fixed** (ledger 87); the per-package `q` scales listed
      above are still open.

71. **Neither the decade colorbar NOR the symmetric centering survives `viz.mosaic`.**
    - What: `_log_decade_colorbar` puts real-unit ticks (`0.1x`, `1x`, `10x`) on a log
      map, and the `change` policy centers it with symmetric `zmin`/`zmax`.
      `viz.mosaic` moves every panel onto a shared `coloraxis` and **discards both**:
      it copies only `trace.colorscale` and then recomputes `cmin`/`cmax` from the
      pooled `z` data, adding `cmid=0.0` only when the caller passes `diff=True`
      (viz.py:301-313). The per-trace `zmin`/`zmax` are ignored once
      `trace.coloraxis` is set, and the coloraxis has no `colorbar` key at all.
    - Measured, two `change` panels spanning 0.5x–10x: `mosaic(diff=False)` gives
      `cmin=-0.301, cmax=1.0, cmid=None`, putting the white "no change" point at
      fraction **0.231** of the diverging scale instead of 0.5. Since stop 0 is red,
      cells that *increased* by up to ~2.2x then render on the red half and read as
      decreases. `mosaic(diff=True)` gives `cmin=-1.0, cmax=1.0, cmid=0.0` — centered.
      Tick labels are lost in both modes.
    - Why deferred: teaching `mosaic` to carry a colorbar, or to honor an
      already-symmetric panel range, is a shared-surface change affecting every mosaic
      in the project; 6.4B builds `field_mosaic` and is the right place to decide it.
    - **Rule for 6.4B:** a `change` mosaic MUST pass `diff=True`, or its white point
      silently moves off 1x. A `mean`/`std` mosaic is sequential, so centering is moot
      there and only the tick labels are at stake.
    - No runtime impact in 6.4A: `plot_field` returns a `viz.Fig`, and mosaicking those
      raises (`choroplethmap` is not compatible with an `xy` subplot), so only a future
      helper passing `Choro` panels can reach this.
    - **HALF-RETIRED 2026-07-26 (6.4B).** The colorbar half is fixed at the source:
      `viz.mosaic` now takes `colorbar=`, a dict **or** a `(cmin, cmax) -> dict`
      callable. The callable form is what the problem actually needs — the useful ticks
      depend on the pooled limits, and those exist only inside `mosaic`. Purely
      additive: the whole block is guarded on `colorbar is not None`, so omitting it
      leaves every pre-existing mosaic untouched. (The pin,
      `test_a_mosaic_without_a_colorbar_argument_is_untouched`, asserts the two
      colorbar properties are unset — the byte-identity is by inspection of the guard,
      not by the test.)
      The centering half is **unreachable rather than fixed**: `plot_field_mosaic`
      accepts only `mean`/`std`, since `change`/`reduction`/`base` have no separate
      prior and posterior form to compose, so no `change` mosaic can be built through
      the front door and the `diff=True` rule above never has to be remembered. A
      hand-rolled `viz.mosaic` of `change` panels still needs it.

72. **Run context (iterations, realization count) rides the figure title, not the hover
    footer — the plan asked for the footer.**
    - What: plan §6.4 item 1 and the handoff both specify a hover footer carrying
      iteration and realization count. `plot_field` puts them in the title instead.
    - Why: `HoverSpec._render_footer` has `if/elif` branches for exactly
      `period`/`step`/`date`/`area`/`model` and no `else`, so `footer=("iteration",)`
      renders **nothing** and raises nothing. Extending that vocabulary for one caller
      would add two customdata columns per cell for a value identical in every cell,
      and would still be invisible on the matplotlib backend, which has no hover at all.
      The title is free, visible on both backends, and survives `report()`'s HTML.
    - Also: the counts come from `self.prior._df.shape[0]` / `posterior`, NOT from
      `settings.num_reals` — that is the *requested* count and is `None` for a run
      myflopy did not launch. (This entry also used to cite the forecast-less
      `TypeError` as a second reason; that was entry 73, **fixed in 6.4C**.)

73. **`IesResults.settings` and `report()` raised on a forecast-less run —
    RESOLVED in 6.4C (2026-07-26).**
    - What was wrong: pyEMU's `Pst.forecast_names` returns `None` when the `forecasts`
      option was never written, so `list(None)` raised `TypeError` in `settings`,
      `report()`, `forecasts()` and `forecast()`. A run with no `cal.forecast(...)` is
      the ordinary case, not an exotic one — `PestProject._apply_forecasts` early-returns
      without writing the option — so the whole review surface died on it.
    - Fixed at the single source (`IesResults.forecast_names`) rather than the four
      call sites: the declared return type was already `list[str]`, nothing anywhere
      branched on the `None`, and `forecasts()` already carried an `if not rows` guard
      that was simply unreachable. The build side already did the same thing
      (`PestProject.settings` uses `.get("forecasts", "")` → 0).
    - The fix filters falsy names rather than only guarding `None`: pyEMU splits an
      empty option string into `[""]`, which would report one forecast and send
      `report()` looking up an ensemble column named `""`.
    - Pinned by six tests over a real `pyemu.Pst` (the existing `_stub_ies_results` has
      no `pestpp_options`, so it raises `AttributeError` before reaching the defect).

74. **`logscale` maps read in log10 units unless the caller relabels the colorbar
    (pre-existing, fixed only for `plot_field`).**
    - What: nothing in `choros.py` sets `tickvals`/`ticktext`, so **every** `Choro`
      built with `logscale=True` labels its colorbar `-3 … 2` for a field that runs
      0.001 … 100 — on both backends. 6.4A adds `_log_decade_colorbar` (Plotly) plus
      `_relabel_log_colorbar` (matplotlib, since `plot_mpl` has no tick hook) and wires
      both into `plot_field` only.
    - Why not fixed generally: the right home is `Choro` itself — it knows it is
      log-scaled — but that changes every existing log map in the project, including
      stored notebook output, and belongs with the `viz.mosaic` colorbar question
      (ledger 71) rather than inside a PEST change.
    - Caught by review, not by design: the first cut of 6.4A put the decades on the
      Plotly branch only, which *introduced* a regression — the matplotlib K map had
      been linear before, so its colorbar used to read real K. Pinned now by
      `test_a_log_map_reads_in_real_units_on_both_backends`.
    - Consequence today: a hand-built log `Choro` (e.g. the truth-K map in
      `canonical_06` cell 18) still reads in log10. That notebook says so in prose.

## 6.4B — IES uncertainty reduction, mosaics, residual maps (2026-07-26)

75. **The uncertainty map became a `stat`, not the planned `field_uncertainty_map`
    method.**
    - What: the plan asked for `field_uncertainty_map(target, which="std"|"reduction")`.
      Shipped instead as `plot_field(target, stat="reduction")`, plus a `reduction`
      column on `field()`.
    - Why: `plot_field(stat="std", which="prior"|"posterior")` **already shipped in
      6.4A** and is exactly the posterior-sd map, so the planned method would have been
      a second spelling of an existing figure — the `foo()`/`plot_foo()` duplication
      `docs/view_layer_conventions.md` forbids. Worse, `which="std"` would have given
      `which` a second meaning inside a class where it means prior-vs-posterior in
      `plot_field`, `parameters_at_bounds` and `plot_parameters_at_bounds`. The only
      genuinely new content was `1 - posterior_std / prior_std`, one derived column
      (`field()` already returned `prior_std`).
    - Cost: the plan's literal name does not exist. Anyone following the old plan text
      finds `stat="reduction"` named in the same item and in the error message that
      lists every valid stat.
    - Judgment call inside it: the reduction scale is anchored to `[0, 1]` — an
      absolute frame — rather than autoscaled like `std`, so two layers are comparable
      and a field that only reduces 0.95–1.0 does not stretch into a dramatic-looking
      map of a trivial range. It **falls back to the data range** when any cell is
      negative, because a posterior spread that *grew* must not clip to the bottom
      color, where it would read as "the data said nothing" instead of "worse than the
      prior".

76. **Lake and SFR observations cannot be placed on the residual map (deferred).**
    - What: `plot_obs_residuals` covers head targets (points) and DRN zones (cells).
      `LakeStageTargets`, `SfrStageTargets` and `SfrFlowTargets` are silently absent
      from the map — but `obs_residuals()` is likewise silent, so the figure never
      claims to be complete, and a run with *only* those targets raises a `ValueError`
      naming the reason rather than drawing an empty grid that would read as zero
      residuals everywhere.
    - Why: they persist only a lake or reach **number** (`observations.py` snapshots
      `["name", "lake"]` / `["name", "reach"]`) — no geometry, no cells. Heads and DRN
      zones resolve from the snapshot files **alone**; lakes and reaches would have to
      resolve id → cellid against the live model at review time. (The *map* already
      needs `IesResults.model` for the grid, so that part is not new — what is new is
      depending on the model's LAK/SFR package contents, and on lake/reach numbering
      being unchanged since the build. The `obs_residuals()` table needs no model
      today, and would start to.)
    - Undeferring: give the lake/SFR `prepare_*` functions a cell snapshot at build
      time, the way DRN zones already have one. That fixes it at the source and needs
      no model at review time — but it only helps runs built *after* the change.

77. **`plot_obs_residuals` has no `map=` flag, unlike the plan's `map=True`.**
    - What: the plan named the method `plot_obs_residuals(map=True)`. Shipped as
      `plot_obs_residuals()` — it *is* the map.
    - Why: nothing was specified for `map=False`, and a boolean with one honest value
      is a worse API than a named method. The natural `map=False` figure is a
      simulated-vs-measured 1:1 scatter, which is a different figure with a different
      argument list, and `plot_vs_obs` already covers most of that ground.
    - If a 1:1 scatter is wanted later it should be its own `plot_*` method, not a flag
      on this one.

78. **The residual sign convention differs from PEST's own `.res` file.**
    - What: `residual = simulated - measured`. PEST's `.res` reports
      `measured - modelled`, i.e. the opposite sign.
    - Why not matched: `phi_contributions` (shipped, ies.py) already computes
      `weight * (simulated - obsval)`, so matching PEST would have made the two
      surfaces of one class disagree. Per the standing q-sign rule, the frame is
      **named** rather than signed into agreement: the column is documented, the hover
      label reads `residual (sim − meas)`, and the figure title says
      `residuals (simulated − measured)`.
    - Consequence: a reader coming from a PEST `.res` file sees flipped colors. Every
      label on the figure says which convention it is in.

79. **`mpl_colormap_for` was extracted from `Choro.plot_mpl` rather than duplicated.**
    - What: the static residual map scatters points over the cells and must use the
      *same* colormap the cells were drawn with, or a point and the cell beneath it
      render different colors for the same value. The colorscale→colormap conversion
      moved out of `plot_mpl` into a module-level `mpl_colormap_for` in `choros.py`.
    - Why recorded: this is a shared-surface change made for one caller. It is a pure
      extraction (`plot_mpl`'s behavior is unchanged and its tests still pass), but it
      widens `choros.py`'s public surface by one name.
    - Note it inherits the `'rdbu'` → `'RdBu_r'` trap (ledger 69/70): passing a
      diverging scale **by name** still mirrors between backends. Its docstring says so.

80. **`obs_residuals` uses `weight > 0` as its "is this a forecast?" proxy.**
    - What: `cal.forecast(...)` sets are persisted as ordinary observation sets — the
      `is_forecast` flag lives on the prepared item, not on the metadata dict that gets
      written (`project.py`) — so a completed run cannot tell a forecast target from a
      calibration target. Zero weight is used instead.
    - Why acceptable: a residual map answers "where was the model fitted badly", and a
      zero-weight observation was not fitted at all. Excluding *deliberately silenced*
      targets along with forecasts is the same answer to the same question, and it
      matches `phi_contributions`, which already filters `weight > 0`.
    - Cost: someone who zero-weights an observation but still wants to see its misfit
      spatially cannot. `obs_residuals()` is the escape hatch only in the sense that it
      applies the same filter — there is no `include_zero_weight=` today.
    - Found by review, not by design: the first cut drew the canonical model's own
      forecast point as a calibration residual, where its 0.799 was the second-largest
      value on the map and could have set the color scale for everything else.

81. **Unmeasured-time rows are excluded by joining to `*_target_values.csv`, and the
    filter silently does not apply when that snapshot is missing.**
    - What: pyEMU creates one observation per row of the simulated output. Rows with no
      target keep pyEMU's defaults — weight 1.0 and an `obsval` equal to the base
      model's own simulated value — so they read as perfect fits and drag a location's
      mean residual toward zero (measured: with 4 phantom rows against 2 real ones, a
      +0.40 residual became −0.10, i.e. the marker changes color).
    - The fix joins back to the measured `(prefix, location, time)` triples in the
      values snapshot. `_measured_observation_keys()` returns `None` when **no**
      observation set recorded a `values_file`, and the caller then does not filter.
    - Why that degradation: dropping every observation would be worse than the bias,
      and a run old enough to lack the snapshot is exactly the run whose metadata we
      cannot reason about. A run built by current myflopy always writes it.
    - Note the underlying oddity is upstream and NOT fixed here: those phantom rows
      carry weight 1.0, so PESTPP is history-matching the base model's own output at
      times nobody measured. That is a build-side question (`observations.py`), with
      its own tests, and is out of scope for a review-layer change.

82. **`run.pest_runs` attaches its model lazily, and gives up silently on a
    multi-model run.**
    - What was wrong (found scoping 6.4C): `Run.pest_runs` called
      `find_pest_runs(root)` with no `model=`, so every handle carried `_model=None`
      and `run.pest_runs[i].review()` returned an `IesResults` that refused **every**
      spatial map — while `canonical_04` and `canonical_06` both told the reader it
      was interchangeable with `model.pest_runs`. (`package_api_reference.md` did not:
      it never mentioned `run.pest_runs` at all until this pass added it.)
    - Fixed by attaching a **`model_factory`** (a zero-argument callable) rather than a
      model, so discovery stays side-effect free: resolving eagerly would build a model
      view, cache it on the run, and raise on a run holding several models — all while
      the caller may only want to list what was calibrated. Note the cost argument is
      about side effects, **not** load time: `Run.model()` returns a lazy
      `LoadedMf6Run` for a reopened run and does not read the MF6 files itself.
    - The compromise: `Run._pest_review_model` swallows `KeyError`/`ValueError` from
      `Run.model()` — a run with no discoverable model, or several and none nominated,
      yields `model=None` and the "Spatial maps need the model grid" message (widened
      in this pass to name `review(model=...)`, which it did not before). Guessing one
      of several models would be worse than declining.
    - Not done, and a real limitation: the factory takes no arguments, so it cannot use
      the `model_name` the handle already carries from the build metadata. A
      multi-model run therefore declines even when the run's own metadata names which
      model was calibrated. Threading that through would remove the compromise above.
    - Fixing the code rather than the docs was chosen deliberately: it makes the two
      notebooks' existing claim TRUE without editing a canonical notebook.
    - Also inherited, not introduced: a reopened run's grid gets `load_mf6_run`'s
      default `crs="EPSG:2927"`, since `Run.model()` passes no CRS. `review()` on a
      reopened run of a model in another CRS will therefore place its spatial maps
      wrongly on the basemap. Pre-existing in `Run.model()`; this change is what makes
      it reachable from the review layer.

83. **`report()` still dies on a workspace with no phi output (deferred).**
    - What: `report()` calls `plot_phi` and `plot_phi_distribution` unguarded, while
      the diagnostics below them sit in `try/except`. `plot_phi` raises
      `FileNotFoundError` when `<case>.phi.actual.csv` is absent; `plot_phi_distribution`
      raises `ValueError` when every phi is non-finite (including all-zero phi, since
      `log10(0)` is filtered out).
    - Why not fixed in 6.4C: 6.4C's remit was the **forecast-less** `TypeError`, which
      is a contract bug. Making `report()` skip its two headline figures is a behavior
      change — a report that silently omits phi convergence looks like a report, and a
      missing phi file usually means the run did not finish, which the caller should
      hear about.
    - Explicitly do **not** widen the `except ValueError` guarding `plot_vs_obs`: that
      `ValueError` is the deliberate "no nonzero-weight groups" signal, and widening it
      would also swallow the `FileNotFoundError` from a missing ensemble.

84. **The phi-contribution pie charts still take the backend's default colors
    (deferred).**
    - What: `plot_phi_contributions(kind="pie")` passes no colors on either backend, so
      the same observation group is one color in Plotly and another in matplotlib.
      Observation groups are named categories and `viz.category_colors` /
      `PALETTE.categorical` exists for exactly this.
    - Why not fixed in 6.4C: it changes rendered output on both backends (not a
      refactor), and `category_colors` memoizes for the **process lifetime**, so the
      assignment depends on which figure a session drew first — acceptable for a
      standalone map, but it needs thinking about inside a `report()` bundle.
    - This is the last place in `ies.py` where a semantic color is left to a default;
      the audit in 6.4C found every other figure on the viz front door.

85. **`edgecolor="black"` in `plot_obs_residuals` stays a literal: PALETTE covers
    semantic colors, not structural ones.**
    - What: the residual scatter outlines its markers in black. The 6.4C color audit
      asked whether that belongs in `viz.PALETTE`.
    - Decision: no. The marker's *meaning* is carried by its fill (`c=residual` on the
      shared diverging cmap); the edge only separates a marker from the cell beneath
      it. `PALETTE` has no edge/outline member, and adding one for this call site would
      mean reconciling the other structural literals scattered across three modules in
      the same pass (`contour_plotting.py` `#222222`, `package_plotting.py` `#666666`
      ×2, and in `choros.py` both the `edgecolor` parameter default and the `"black"`
      region boundary).
    - Known asymmetry, accepted: the Plotly twin draws no marker outline at all.

86. **`mf.PestProject` is documented as the advanced path but still exported as
    preferred (deferred).**
    - What: 6.4C documents `model.pest(name, ...)` as the front door and direct
      `PestProject(...)` construction as advanced, in the capability map, the package
      API reference and the class docstring. But `PestProject` is still a first-tier
      name in `mf.__all__`, not in `mf.__engine__` alongside the other builders.
    - Why not moved: demoting it removes it from `__all__`/autocomplete and changes
      `tests/api_snapshot.json` — a public-API change, which plan §6.4 item 6 does not
      ask for and which should be an explicit decision rather than a side effect of a
      documentation pass.
    - Related asymmetry, also untouched: `find_pest_runs`/`PestRunHandle` are exported
      at top level while `open_ies_run`/`IesResults` are not.

87. **Every grouped difference map rendered mirrored on the static backend
    (pre-existing, FIXED 2026-07-27).**
    - What: `get_default_group_compare_colorscale()` returned the NAME `"RdBu"`, and
      `choros.py`'s `_PLOTLY_TO_MPL_CMAP` maps `"rdbu"` to matplotlib's **reversed**
      `"RdBu_r"`. So all six diff-map families — `group/spatial.py`, `inputs.py`,
      `results.py`, `sfr.py`, `uzf.py`, `lak.py` — put red where `plot()` put blue.
      One figure said "lower here" interactively and "higher here" statically.
    - Found while scoping §6.1 items 3–4 (Δconc/Δtemp need the same scale); the user
      chose to fix all six at once rather than ship stops only for the new maps and
      leave the transport diffs disagreeing with `diff().hds`.
    - **This is a visible change** to the static rendering of existing diff figures.
      It is a correction, not a restyle — but a report generated before today and
      regenerated after will not match.
    - The old test asserted `== "RdBu"`, i.e. it pinned the *name* and so passed the
      entire time the maps were wrong. Replaced with an assertion on the rendered ends
      of BOTH backends, per the front-door pattern established in 6.4A.
    - Judgment call inside the fix: the stops are **stated twice** — in
      `package_registry._GROUP_COMPARE_DEFAULT_COLORSCALE` and in
      `package_plotting._red_white_blue_diverging_colorscale()`. Importing the latter
      into the registry was tried and reverted: the registry is the low-level source of
      per-package truth, and depending on a drawing module inverted the import layering
      (measured: it pushed ~18 modules up a layer). `test_the_two_diverging_orientations
      _do_not_drift` asserts the two are equal, so the duplication cannot rot.

88. **Non-field maps borrowed `model.hds` for their time axis, which broke every
    package map on transport models (pre-existing, FIXED 2026-07-27).**
    - What: `Choro` resolves a dependent-variable reader from its `type`, and
      `_DEPVAR_READER_ATTR` knows only `hds`/`conc`/`temp`. Every other type — including
      `type='custom'`, which is what the whole package-explorer grammar uses — fell
      back to `self.model.hds` for `kstpkper` and the layer table. On a GWT/GWE model
      §6.0's kind guard turns that into `AttributeError: model 'x' is a GWT model;
      '.hds' (heads) is only available on GWF models`.
    - Effect: no package input/result map, no group diff map, and no budget map worked
      on a transport model at all — the failure was in the *timing* lookup, before
      anything was drawn, which is why it looked unrelated to color or data.
    - Fixed by borrowing the model's OWN dependent variable (`Choro._timing_reader` →
      `SimulationBase._field_reader` → the existing `accessors.field_reader`, which was
      written for exactly this and had one caller). `_field_reader` is deliberately
      private and **ungated**, unlike `.hds`/`.conc`/`.temp`: kind-agnostic machinery
      needs *a* field without asserting which physics it is.
    - Kept the `getattr(self.model, "_field_reader", None) or self.model.hds` form so
      duck-typed stand-ins that are not a `SimulationBase` keep working.

89. **Reopening any multi-model run raised `RecursionError` (pre-existing, FIXED
    2026-07-27).**
    - What: a simulation with several models gives each one its own subdirectory
      (`workspace._materialize_models`), so a coupled GWF+GWT/GWE/PRT run keeps nothing
      but `mfsim.nam` at the top. `LoadedMf6Run` discovered packages by globbing the
      workspace ROOT, found none, left the grid type `"unknown"`, and the `grid_type`
      property then fell back to inspecting `self.gwf` — whose loader consults
      `grid_type`. `load_run(ws).model(name)` recursed to death for **every** model of
      a coupled run, flow or transport.
    - Four distinct defects, all fixed: (a) package/grid discovery now looks in
      `<workspace>/<model>/` via `_model_dir`; (b) `_infer_model_name` falls back to
      one level of subdirectories; (c) `_ensure_core_loaded` uses the DISCOVERED
      `_grid_type_override` instead of the `grid_type` property, which makes the cycle
      structurally impossible rather than merely unreached; (d) the constructor's
      progress message interpolated `self.grid_type` — so a "lazy, file-backed" loader
      loaded the entire FloPy simulation just to build a string, **even at verbosity 0
      where the string is discarded**.
    - (c) needs its own test: once (a) works, `grid_type` answers from the override and
      never reaches `self.gwf`, so the guard is unexercised by the multi-model case.
      Pinned separately against a workspace with no recognizable grid package.
    - Also fixed alongside: `LoadedMf6Run` never set `model_type`, so a reopened GWT/GWE
      model reported the class default `'gwf6'` and offered `.hds` — silently the wrong
      physics under a familiar name, with the kind guard unable to help because it
      believed the lie. Now read from `mfsim.nam`'s `models` block.
    - Scope note: found while scoping §6.1 items 3–4, and fixed FIRST at the user's
      direction, because it meant those features would work only on a live built run
      and appear broken to anyone returning to a run in a later session — the
      documented purpose of `load_run`.

90. **`mf.gwt`/`mf.gwe` now default `save_flows=True` — a build-side behavior change
    (2026-07-27).**
    - What was wrong: MF6 writes a **zero-byte** `.cbc` for a transport model unless
      `SAVE_FLOWS` is set on the model itself, even when OC asks for a budget. FloPy then
      raises `ValueError: datafile error: file is empty`, which names neither the model
      nor the missing option. Every myflopy transport model was in this state, because
      `mf.gwt`/`mf.gwe` defaulted `save_flows=None` while `mf.npf`, `mf.sto` and every
      list BC already default it `True`.
    - The change: transport models now save flows by default, so asking OC for a budget
      yields one. **Existing transport models start writing a `.cbc` they did not
      before** — more disk, and a file appears where none used to.
    - Scoped to `mf.gwt`/`mf.gwe` only. `mf.gwf` keeps `save_flows=None`: a flow model
      already gets its budget from NPF and the list BCs, so flipping it there would be a
      behavior change to every existing model for no gain.
    - Measured, not assumed: model-level `SAVE_FLOWS` alone is sufficient — stripping it
      from MST/SSM and re-running still wrote all three terms. So the transport package
      helpers (`mf.mst`/`mf.ssm`/`mf.est`/...) were left alone.
    - The terms are not what plan §6.1 item 3 guessed. Measured on real runs:
      GWT = `STORAGE-AQUEOUS` / `FLOW-JA-FACE` / `SOURCE-SINK MIX`;
      GWE = `STORAGE-CELLBLK` / `FLOW-JA-FACE` / `SOURCE-SINK MIX`. There is **no term
      named `SSM`** (the SSM package's record is `SOURCE-SINK MIX`), and `DECAY` appears
      only when MST declares decay or sorption.

91. **Transport budgets were unreadable four different ways — all FIXED, with two
    deliberate limits (2026-07-27).**
    Found while scoping §6.1 item 3 (`model.budget.<term>`), reproduced against real
    coupled GWF+GWT and GWF+GWE runs, and fixed under the user's "fix everything
    found" decision. What was wrong:
    - **(a) Two of a transport model's three terms returned an EMPTY table, silently.**
      MF6 writes storage and `FLOW-JA-FACE` as imeth=1 *full arrays* with no `node`
      column. `pd.DataFrame.from_records` turns a `(nlay, 1, ncpl)` float array into a
      nonsense `(1, 1)` frame with an integer column name, which the builder's
      `"node" not in frame.columns` guard then skipped — `table=(0, 7)`, no error. A
      caller could not tell "no flow" from "not read".
    - **(b) The two budget paths disagreed about node basing.** The modern explorer
      zero-based with its own rule; `model.bud()` looked the package up in the
      registry, which no transport record matches — so the legacy path handed back
      MF6's raw 1-based ids. Measured against the CHD source cells: `0,6,7,8,14…`
      (right) vs `1,7,8,9,15…` (wrong).
    - **(c) Units were hardcoded `ft³/d`** at both hover sites — the wrong *dimension*
      for GWT (mass/time) and GWE (energy/time), and the wrong *units* for any GWF
      model not declared in feet and days.
    - **(d) `model.bud("ssm")` was unreachable on every transport model.** MF6 names
      the record for the process (`SOURCE-SINK MIX`), not the package, and the lookup
      substring-matched package names — so the obvious call raised while the
      non-obvious literal worked.
    - **The fix is one shared converter.** `_model_budget_record_frame` handles both
      record shapes and zero-bases `node` exactly once, and BOTH paths now call it, so
      they cannot drift apart again. Cross-checked physically, not just structurally:
      the storage and source-sink terms cancel to 1e-6 relative on both kinds
      (GWT ∓27.26, GWE ∓1.141e8), which pins the table's *values*, not only its shape.
    - **Deliberate limit 1 — `FLOW-JA-FACE` raises instead of returning a table.** It
      is indexed by cell *connection* (388 values on a 64-cell grid), so there is no
      honest mapping onto cells. Exposing it would need a connection-indexed noun,
      which is a different feature; the bug being fixed here was *silence*, so refusing
      loudly is the correct end state for this pass. **How that refusal is decided is
      itself load-bearing — see 96.**
    - **Deliberate limit 2 — the helper lives in `package_explorer_utils` (layer 0),
      not in `package_budget` (layer 2).** Importing `package_budget` from
      `budget_tables` worked but pushed `budget_tables`/`budget_plotting`/`budget`/
      `outputs` each up one layer. Moving the two shared functions down instead keeps
      the code shared with zero layering movement.

92. **A fifth off-by-one, found while fixing 91 and fixed with it: `model.bud("sfr"
    |"lak"|"uzf").df` returned node ids one cell high (2026-07-27).**
    - The registry marked these three `zero_base_budget_nodes=False`, and both it and
      `_cell_based_budget_packages`' docstring justified that by saying their records
      "carry feature ids, not model cells". **Measured false.** In the *model* budget
      file every record's `node` is a 1-based model cell, `SFR`/`LAK`/`UZF-GWRCH`
      included; it is `node2` that holds the feature id there. The feature-first layout
      the belief described is real but belongs to the separate *package-output* budget
      file that `model.outputs.<pkg>.bud` reads — the two were conflated.
    - This is exactly the off-by-one class `tests/test_budget_node_basing.py` exists to
      prevent; it escaped only because that module parametrizes `cell_stress` packages
      alone. It now covers the advanced three as well, against ground truth taken from
      `packagedata`/`connectiondata`.
    - Nothing in `src/`, `tests/`, or `examples/` consumed those ids, so no caller
      shifts under the fix.

93. **`zero_base_budget_nodes` is now misnamed — rename DEFERRED (2026-07-27).**
    - After 91/92, `node` is zero-based once and universally, in the shared record
      converter. The descriptor field and `_cell_based_budget_packages()` survive but
      now govern only `node2`, whose meaning genuinely is per-package (a boundary index
      for the list BCs, a feature id for SFR/LAK/UZF).
    - Why not renamed now: the field is pinned by `tests/api_snapshot.json` and counted
      by `tests/test_package_descriptor_payoff.py`, so a rename is a separate,
      mechanical pass that would otherwise be buried inside a correctness fix. Both
      docstrings state the narrowed scope in the meantime.
    - The registry keeps its 9 automatic surfaces: this narrows what one of them
      *means*, it does not remove it.

94. **`node2` basing is left asymmetric — across packages AND across the two paths —
    deliberately (2026-07-27).**
    - The asymmetry is real and larger than "one rule with exceptions". Measured on the
      canonical model, `model.packages.<pkg>.results.q.get()` vs `model.bud(<pkg>).df`:
      - `drn` (and the other cell BCs): **0..11 both paths — they agree.**
      - `sfr`: **modern 0..54, legacy 1..55 — the two paths differ by one.** The modern
        builder zero-bases `node2` with a min>=1 heuristic for every package; the legacy
        one consults the registry, which excludes SFR.
      - `lak`: **modern 162..292, legacy 1..1 — not even the same quantity.**
        `lak.results.q` reads the separate *package-output* budget file, where `node2`
        is a GWF cell, while `model.bud("lak")` reads the model budget, where it is the
        lake id.
    - All of it predates this pass and is preserved bit-for-bit rather than unified;
      this change touched only `node`.
    - Why: nothing reads `node2` as a meaningful id on the model-budget path — the only
      consumer is `group/_shared.py`, which uses it as an opaque diff *alignment key*,
      where both sides shift together and absolute basing is irrelevant. Zero-basing it
      for LAK would silently redefine `model.bud("lak").df.node2` from a 1-based lake id
      to a 0-based one for no benefit, and the honest fix is the rename in 93, not a
      second basing change smuggled in here.

95. **Budget hover units fall back to dimensional placeholders, and two package-output
    hovers still hardcode feet (2026-07-27).**
    - `budget_value_units` reads `dis.length_units`/`tdis.time_units` and returns
      `ft³/d` for the canonical model, which declares feet and days. A model that
      declares neither now reads `L³/T` (GWF), `M/T` (GWT), `E/T` (GWE) rather than
      guessing feet and days as before. **This changes the hover on any undeclared GWF
      model** from a confident wrong label to an honest dimensional one.
    - Not converted in this pass: the SFR and LAK *package-output* hovers
      (`hover.py:675,688`), which hardcode `ft`/`ft²`/`ft³/d` across several fields.
      Most of those are geometry (`rlen`, `stage`, `flow_area`), not budget values, so
      they need a length-unit helper rather than this one, and both packages are
      GWF-only — no transport model reaches them. Separable; left for the pass that
      needs it.

96. **The connection-indexed guard keys on a hardcoded record NAME, because size
    cannot decide it (2026-07-27).**
    - Caught by adversarial review of the 91 fix, before commit. The first version of
      the guard discriminated a connection-indexed full array from a cell-indexed one
      by **length alone** (`values.size != len(model.node_to_lni)`). That is unsound:
      under an idomain reduction MF6 expands cell arrays back to `nodesuser` while
      `FLOW-JA-FACE` stays at the reduced `nja`, so the two counts are independent.
    - **Reproduced, not theorised.** A 1-layer DISV with `ncpl=4` and two active
      adjacent cells gives `nja == nodesuser == 4`, and *both* records arrive with the
      identical shape `(1, 1, 4)` — so neither size nor shape can separate them.
      `model.bud("FLOW-JA-FACE").df` returned a four-row "cell" table reporting
      `q=+16.85` through a cell whose `idomain` is 0. Pre-fix, that same call raised
      loudly, so the size-only guard was **strictly worse than what it replaced**, on
      the exact failure mode the fix existed to remove.
    - Why a literal set and not a registry field: this is an MF6 **file-format** fact
      (which records are written imeth=1 over connections), not a per-package one, so
      it does not belong on a package descriptor. If MF6 adds a second such record it
      must be added to `_CONNECTION_INDEXED_BUDGET_RECORDS` — that is the known cost of
      this approach, accepted because the alternatives are worse: the header's shape is
      ambiguous on 1-layer models, and length is provably wrong.
    - The guard resolves the caller's string to the real record name first, so the
      substring route (`model.bud("flow")` → `FLOW-JA-FACE`) is caught too; a guard
      that tested the caller's own string would have missed exactly that alias.
    - Test-fidelity note: every other budget test uses an all-ones idomain, where
      `nja > nodesuser` is guaranteed, so none of them could see this. The regression
      test builds and runs the degenerate 4-cell model, and **asserts the collision
      still holds** before asserting the refusal — otherwise a future grid change would
      silently turn it into a test of nothing.

97. **`model.budget.<term>` — four judgment calls behind the new noun (2026-07-27).**
    Closes plan §6.1/6.2 item 3 (ledger 56 sub-item 3).
    - **Terms are DISCOVERED, not declared — diverging from the package namespaces.**
      `lak.budget.<term>` / `sfr.budget.<term>` hand-write one `@property` per term with
      the MF6 record name as a literal. That is untenable here: the model budget's term
      set varies with kind and configuration (GWT storage is `STORAGE-AQUEOUS`, GWE's is
      `STORAGE-CELLBLK`, `DECAY` appears only when MST declares it, and a flow model's
      terms are whichever packages it carries). So this namespace resolves through
      `__getattr__`/`__dir__` against the live file. The inconsistency with the package
      namespaces is accepted deliberately; the alternative was a declared list that is
      wrong for most models.
    - **Available on every model kind, not gated to transport** (user decision). The
      plumbing underneath is kind-neutral, so gating would have been extra work to
      *remove* capability — and `STO-SS`/`STO-SY`/`DATA-SPDIS`/`DATA-SAT` have no package
      accessor at all, so this is their only route that is not the legacy `model.bud()`
      wrapper. Accepted cost: for the boundary packages it is a third spelling alongside
      `model.packages.<pkg>.results.q` and `model.bud(pkg)`. They are pinned to agree
      numerically by `test_a_term_agrees_with_the_package_results_path`.
    - **The column stays MF6's raw `q`, not `q_gwf`.** `model.packages.drn.results.q`
      renames it `q_gwf` so the column carries its reference frame. This namespace is
      explicitly the raw model-budget view, and the frame-naming convention does not
      generalize: `q_gwf` is meaningless for `STORAGE-AQUEOUS`, which is a mass rate, not
      a gwf-referenced exchange. Values are identical either way — verified — and the
      hover unit label already states the dimension (`M/T`, `E/T`, `L³/T`).
    - **`FLOW-JA-FACE` is listed but refuses `get()`** (user decision). Omitting it would
      put `dir()` and `types` back in disagreement, which is exactly what makes the
      package namespaces confusing — a term present in the file but not declared there is
      silently unreachable by attribute. The error explains the reason (ledger 96).

98. **Two pre-existing defects in the PACKAGE-level `budget.<term>` namespaces, found
    while scoping 97 and NOT fixed (2026-07-27).**
    - **`budget.get(term="ext_inflow")` silently returns an empty frame.** The only term
      normalizer on that path (`_normalize_term_filter`) upcases but does not convert `_`
      back to `-`, so `"ext_inflow"` becomes `"EXT_INFLOW"` and never matches
      `"EXT-INFLOW"`. Only the attribute route (`budget.ext_inflow`) is safe. The new
      `budget_term_attribute` normalizer added in 97 is the missing piece, but wiring it
      into the package path changes behavior for `lak`/`sfr` term filters, which is a
      separate change with its own blast radius.
    - **`PackageBudgetTermExplorer` has no `plot()`.** It answers `types`/`get`/`summary`/
      `wide` only, so it violates `view_layer_conventions.md`'s "every noun answers the
      same verbs" rule (already noted at ledger entry 62). The new model-level noun does
      NOT inherit this — it returns `CellBudgetResultsExplorer`, a `SpatialView` with the
      full verb set — so the two tiers now differ in what verbs a budget term answers.
    - Why deferred: both are in the package-output subsystem, which reads a different
      file with different node semantics, and #57 was scoped to the model-level noun.
      Neither is a wrong-data bug; the first is a silent empty result and the second a
      missing verb. Also relevant: **no test anywhere covers `<pkg>.budget.<term>`** — its
      only consumer is `canonical.py:255` — so that tier would need test coverage built
      before it is safe to change.

99. **`xs` was broken on the transport readers from the day they shipped — FIXED
    2026-07-27.**
    - `model.conc.xs()` and `model.temp.xs()` raised `AttributeError: model 'trans' is
      a GWT model; '.hds' (heads) is only available on GWF models`, because
      `XSection.all_heads` cached its value table from `self.model.hds` — the very
      accessor the §6.0 kind guard refuses on a transport model. Meanwhile
      `myflopy_context.md` advertised the transport readers as shipping the "full
      grammar (`get/summary/array/map/xs/mosaic/animate`)". The docs were wrong.
    - Found while scoping §6.1/6.2 item 4, by testing the claim rather than reading it;
      the scoping agent had reported `XSection` as "heads-hardwired" and therefore a
      deep blocker for `GroupConc.xs`.
    - It was not deep. `XSection` reads its values POSITIONALLY (`y[0] for y in
      y_lyr.values`, xsections.py), so the per-kind column name (`elev`/`conc`/`temp`)
      never surfaces. The fix points the cache at `model._field_reader` — the ungated
      kind-neutral reader added 2026-07-26 (ledger 88) — with a `.hds` fallback for
      duck-typed stand-ins. One line of behavior change.
    - **The property is still named `all_heads`.** It now returns whatever field the
      model has. Renaming it would touch every `XSection` consumer for a cosmetic gain;
      the docstring states the mismatch instead. Deferred.
    - Fixing it here rather than separately was a user decision: it unblocked
      `group.conc.xs()` in the same pass, so the group tier did not ship with a hole.

100. **`GroupConc`/`GroupTemp` — the "near-mechanical clone" was refused (2026-07-27).**
    Closes plan §6.1/6.2 item 4 (ledger 56 sub-item 4).
    - Ledger 56 predicted a "near-mechanical `GroupHeads` clone with `elev`→`conc`".
      §6.2 says the opposite: *"if it is not mostly mechanical, stop and fix the 6.0
      abstraction instead of copy-pasting."* Cloning would have made three copies of
      `get`/`compare`/`compare_map` (~120 lines each) differing only in literals.
    - What shipped instead: `GroupHeads` moved onto a new `_GroupFieldView` base — the
      group-side counterpart of `DependentVariableFile` — and all three kinds are now
      **one implementation configured by five class attributes**
      (`reader_attribute`/`table_attribute`/`value_column`/`value_label`/`value_unit`).
      A test asserts none of the three re-implements the shared methods, so a future
      copy-paste fails rather than quietly accumulating.
    - **Cost of that choice, stated plainly:** it rewrites working GWF code that the
      group tier depends on. Mitigated by mutation testing (reversing the diff
      subtraction, repointing the reader/table attributes, and repointing the diff leaf
      are all caught), and the existing heads tests were left untouched as the
      regression net.
    - **`all_conc`/`all_temp` are kind-GATED; `all_heads` is not.** The asymmetry is
      deliberate: `all_heads` predates the §6.0 kind guard and un-gating is not worth a
      behavior change to every existing caller, but a new transport table should not
      inherit the looseness — an ungated `all_conc` on a flow model reads a missing
      `.ucn` and fails far from the cause. Recorded rather than unified.
    - **`GroupConc`/`GroupTemp` live in their own modules, not in `spatial.py`.** Beyond
      the plan's own recommendation, `spatial.py` sits at exactly 1 on the deferred-import
      ratchet (its lazy `XSection` import); adding a second field class there would have
      taken it to 2 and failed `test_deferred_import_ratchet`. Separate modules take a
      fresh allowlist entry instead.
    - **Diff orientation now has a test, for all three kinds.** Nothing pinned that
      `diff` is `model - reference`: the shared colorscale helper is tested in isolation,
      so reversing the subtraction would invert every difference map while
      `test_diff_maps_use_rdbu_negative_red_positive_blue` still passed. The new test
      recomputes the expected difference from `get()` and compares element-wise.
    - **Test fixtures build two models that DIFFER** (porosity 0.2 vs 0.05). A grouped
      comparison over identical models proves the plumbing runs and nothing else; the
      earlier scoping assumed no transport fixture was possible, but
      `irregular_voronoi_grid` is deterministic, so two coupled runs share a grid and
      group cleanly. `_coupled_run` gained a `porosity=` keyword for this.

101. **The canonical transport fixture: the recorded blocker was wrong, and the real
     obstacles were different (2026-07-28).** Closes plan §6.1/6.2 item 6 (ledger 56
     sub-item 6) for the MODEL; the calibration built on it is item 5, closed the same
     day (ledger 103, and the PEST wiring in ledger 104/105).
     - Ledger 56 deferred this as "large — the canonical runs on the single-model
       `SimulationBase`, so a coupled GWT needs the multi-model spec path." Measured:
       **`SimulationBase` already owns a real `MFSimulation`**, and a GWT model, a
       second IMS and a `GWF6-GWT6` exchange attach to it with no library change.
     - Worse, the spec path would have **broken** what it was supposed to enable.
       `Project.prepare_run` gives each model its own subdirectory, so external arrays
       land at `<ws>/flow/flow.npf_k.txt` — and PEST's parameter-file resolution globs
       the workspace ROOT. The flat sibling keeps `viz_prt_master.npf_k_layer1.txt`
       *and* `viz_prt_master_t.mst_porosity.txt` at the root, so calibration works and
       a future transport `parameterize` target is a recipe entry rather than a rework.
       (Binary outputs stay at the root either way — observations alone would have
       survived the spec path, which is what made the wrong choice look viable.)
     - **What actually blocked it, none of it predicted.** Three MF6 refusals, each
       found only by running: (a) `MODELNAME` is capped at 16 characters and the
       canonical name is already 14, so `viz_prt_master_trans` was rejected before
       anything solved; (b) a flow model with boundary packages *requires* an SSM
       package on the transport model; (c) every advanced flow package needs its
       transport counterpart — the canonical's LAK/SFR/UZF/MVR mean LKT/SFT/UZT/MVT or
       MF6 will not run.
     - **The subtlest one: the second solver needs an EXPLICIT `filename`.** Without it
       FloPy auto-names it `<sim>_0.ims`, which lands it BEFORE the flow solver in
       `mfsim.nam` — and MF6 rejects exactly that ("the IMS specified for GWF must be
       listed in mfsim.nam before the IMS for GWT"). Registration ORDER does not fix
       it, in either direction; only the filename does. Worth recording because the
       error names an ordering problem whose cause is a naming one.
     - **The source is a CNC package on the transport model**, not auxiliary
       concentrations threaded through the flow model's boundaries. The flow model is
       left untouched, so `CANONICAL_MODEL_CONTRACT.validate()` — which is entirely
       `model.gwf`-scoped — keeps passing verbatim, and every existing canonical test
       sees the model it always saw. A test asserts that rather than assuming it.

102. **`ConcTargets` keeps `head_target`/`sim_head` as its internal column names
     (2026-07-28).**
     - The normalized target frames name their value columns `head_target` and
       `sim_head` for EVERY field kind, concentration included. Read literally that
       is wrong for a conc target, and a caller doing `targets.to_long()` sees a
       `head_target` column holding concentrations.
     - Why it is kept: those names are woven through `calcs/calibration.py` (the
       `CalibrationPlot` machinery, ~8 sites) and through the PEST observation
       builder, whose value assignment reads
       `getattr(row, "value", getattr(row, "head_target", np.nan))`. Renaming
       per-kind would either fork that machinery or -- if the rename were done
       without also renaming in the PEST builder -- assign **`obsval = NaN` to every
       observation with no error at all**, producing a control file that runs and
       calibrates to nothing.
     - Sharing them buys something real: concentration targets get the existing
       calibration plots, residual statistics and PEST wiring unchanged.
     - Proper fix (deferred): rename to a neutral `target_value`/`simulated_value`
       across `observations/` and `calcs/calibration.py` in one pass, with the PEST
       builder updated in the same commit. Separable, and not safe to do halfway.

103. **The transport calibration demo estimates FLOW parameters from concentration
     data (2026-07-28).** Plan §6.1/6.2 item 5, observation half.
     - `build_canonical_transport_calibration_demo` perturbs K and asks the
       calibration to recover it from *concentrations*, not from heads. That needs no
       new `parameterize` target: `k`/`recharge` already exist, and the flat
       simulation layout (ledger 101) keeps their external arrays where PEST globs.
     - Measured first, not assumed: on the canonical testing profile a 3x K change
       moves concentration at 9 of 12 monitoring wells by more than 1% (max 0.175
       against well values of 0.02-0.25), and the spoiled model misfits all 66
       truth-derived targets (max 0.175, RMSE 0.057).
     - **Closed 2026-07-28: PESTPP-IES was run on this demo and converged** (reported
       by the user from a live notebook run; the build, the misfit and the forward run
       were measured here, the convergence was not). Until then the demo was only known
       to BUILD and misfit -- worth distinguishing, because a well-posed-looking
       calibration that cannot actually close its residuals is a common way for a
       synthetic demo to be quietly useless.
     - **Well placement is load-bearing.** CNC pins the source cells at the source
       concentration for the whole run, so a well there reads 1.0 regardless of K and
       constrains nothing. `monitoring_well_cells` keeps intermediate-concentration
       cells downgradient of the source; a test asserts no well sits on a source cell
       and that no target sits at the source value.
     - ~~**Transport parameters (`mst.porosity`, `dsp.alh`) are still not
       `parameterize` targets.**~~ **`porosity` CLOSED 2026-07-29** (see 110);
       `dsp.alh` deliberately NOT added: `ath1` is derived from `alh` at build time
       only, so a multiplier on `alh` alone silently breaks the 10:1 transverse ratio
       the model was built with. That needs a linked-parameter concept, not a recipe
       entry.

104. **`prepare_conc_observations` keeps a DEFERRED import of `load_mf6_run`, against
     the ratchet's default advice (2026-07-28).**
     - The deferred-import ratchet's failure message says to hoist an import to module
       level rather than raise the allowlist. Here hoisting is not available:
       `pest/observations.py` sits at import layer 6 and `project/run_model.py` at 9, so
       a module-level import would push `pest.observations` to 10 -- ABOVE
       `pest/project.py` (7), which imports it. That inverts the graph and fails
       `test_runtime_imports_never_point_upward`.
     - The import is needed because concentrations live on the GWT sibling while the
       calibration hangs off the FLOW model: `project.model.all_conc` is kind-gated and
       would refuse, so the transport view has to be reopened by name.
     - Same shape as the `package_registry` hoist reverted earlier this session
       (ledger 91): the ratchet's advice is right by default and wrong when the target
       module is higher in the graph. Recorded so the next reader does not "fix" it.

105. ~~**Concentration observations are not yet drawn on IES residual maps
     (2026-07-28).**~~ **CLOSED 2026-07-29.** `prepare_conc_observations` now
     snapshots the point locations `match_to_model` was already computing (and
     discarding), and `_observation_locations` dispatches on a declared
     `geometry` (`"points"`/`"zones"`) instead of a hardcoded kind tuple, reading
     each family's own `mapping_file` rather than assuming the head suffix. Runs
     already on disk keep reviewing via `_LEGACY_GEOMETRY_BY_KIND`.

106. **`canonical_07` is verified by hand, because nothing in this repo executes a
     notebook (2026-07-28).**
     - There is no nbval/nbmake/papermill/nbclient anywhere, and CI runs only
       `pytest -m "not slow"`. The canonical notebooks are checked by *grepping their
       JSON source* for required substrings. That catches a missing API name and
       nothing else — not a typo'd kwarg, not a stale signature, not a cell that
       raises.
     - So this notebook was executed end to end by hand on the full `validation()`
       (50x50) profile before commit: all 11 code cells ran clean. `nbconvert` is not
       installed either, so the run extracted the code cells and executed them in
       order in one namespace, which validates the code but not the notebook
       machinery around it.
     - `canonical_07` also falls outside BOTH existing notebook globs
       (`canonical_0[0-3]` asserts exactly 4, `canonical_0[4-6]` exactly 3), so it
       would have shipped with zero static checking. It has its own test now, which
       pins the ten API names the notebook exists to teach plus the no-committed-
       outputs rule.
     - Deferred: adding notebook execution to CI. It would need a jupyter dependency
       and would move the slow lane's runtime substantially (this notebook alone runs
       two coupled models and a PRT simulation), so it is a CI-policy change rather
       than part of this feature.

107. **PstFrom copied the PEST template into itself; guarded, with a real capability
     limit (2026-07-28).**
     - `PestProject` defaults its template to `<model workspace>/pest/<name>` — INSIDE
       the workspace `PstFrom` copies wholesale. Whenever a `pest/` directory exists
       when the copy starts, `shutil.copytree` descends into the destination it is
       currently filling: `pest/<name>/pest/<name>/pest/...`. A user re-running a
       calibration cell produced a **17 GB tree nested 40 levels deep** before it died
       of recursion.
     - Only a notebook comment documented this (`canonical_04`: "a stale pest/ template
       inside the model dir makes PstFrom recurse"), and the workaround was to
       `rmtree` the whole artifact root — which only works if the cleanup and the build
       are in the same cell. `canonical_07` split them across sections and reproduced
       the bug immediately.
     - **My first guard was wrong and testing caught it.** I assumed a *stale template*
       was required and deleted just that; it still recursed. The mechanism is that the
       destination is inside the source, so an EMPTY `pest/` is enough.
       `PstFrom`'s own `remove_existing=True` does not help — it clears the
       destination, which is not what recurses.
     - What ships: `_clear_stale_template` removes this run's template and, if that
       leaves `pest/` empty, removes `pest/` too — restoring exactly the state a
       first-ever build sees. Measured: three consecutive rebuilds stay at nesting 1
       and 15 MB.
     - ~~**The limit:** a second, differently-named calibration at the DEFAULT
       location could not be made safe this way, so `model.pest_runs` with several
       runs required explicit workspaces.~~ **CLOSED 2026-07-29** by the proper fix
       below; several named calibrations coexist again.
     - **The proper fix, 2026-07-29: the default moved OUT of the copied tree** to
       `<model workspace>.pest/<name>`, a sibling of the model directory. Excluding
       the subtree instead is impossible through supported API: pyEMU 1.4.0 copies
       with a bare `shutil.copytree(o_d, n_d, symlinks=True)` in the private
       `_try_copy_dir` (`pyemu/utils/os_utils.py`), called from `PstFrom.__init__`
       itself — no `ignore=`, no kwarg, no hook. Only monkeypatching a private
       function would work, which would break silently on any pyEMU upgrade and bring
       the recursion back.
     - Moving cost less than this entry assumed. `find_pest_runs` is
       location-agnostic (it `rglob`s whatever root it is handed), the forward run
       uses only basenames, and `original_workspace` in the metadata has no reader.
       Only two call sites hardcoded a root. The one hard constraint is that IES
       masters must stay SIBLINGS of the template (`runs.py` globs
       `template_dir.parent`); break that and a finished run silently lists as
       "built (not run)".
     - `_clear_stale_template` is kept for an explicit `workspace=` inside the model
       directory, and finally has a test — it had none.

108. **Two of my own `canonical_07` test assertions fought the way notebooks are
     actually used (2026-07-28).**
     - `assert "RUN_IES = False" in source` fails the moment someone flips the gate to
       True to run the calibration — which is the entire point of the gate. No other
       test in the repo asserts a gate VALUE; I invented it. Now asserts the flag
       exists.
     - `assert not any(cell.get("outputs"))` fails the moment someone runs the notebook
       locally, since outputs land in the working tree. The existing convention asserts
       no ERROR outputs, which tolerates a local run and still catches a committed
       failure. Now matches.
     - Recorded because both were written from "what should the committed file look
       like" without asking "what does this do to someone working in the file".

109. **`plot_obs_residuals` warns about mixed units rather than refusing them
     (2026-07-29).**
     - One diverging color scale spans every family on the figure, so a run
       history-matched against heads (length) *and* concentration (mass/volume)
       hands the limit to whichever has the bigger numbers. The other family then
       renders uniformly white — which reads as **a perfect fit**, not as a broken
       figure. That is the dangerous failure mode, and it is what closing ledger 105
       makes reachable in practice, since the transport demo now carries both.
     - `prefix=` (on `obs_residuals` and `plot_obs_residuals`) is the fix; the
       default still draws everything and emits a `UserWarning` naming the kinds and
       the way out.
     - **Why warn, not refuse.** The mixing predates this change: head targets (ft)
       and DRN seepage (ft³/d) have always shared the scale, so refusing would break
       reviewing runs already on disk. A warning is visible, keeps old runs working,
       and points at the parameter.
     - Not fixed: a genuinely correct multi-unit figure needs either one colorbar per
       family or normalized residuals, and normalizing changes what the number means.
       Deferred as a real design question, not an oversight.

110. **`parameterize("porosity")` re-stores the model's griddata one array per
     layer, mutating the FloPy model (2026-07-29).** Closes the porosity half of
     ledger 103.
     - MST porosity is supplied as a scalar (`ModflowGwtmst(gwt, porosity=0.25)`),
       so `set_all_data_external()` writes ONE file of `nlay * ncpl` values with no
       LAYERED keyword. pyEMU cannot parameterize that on a Voronoi grid:
       `write_array_tpl` calls `get_xy([i, j])` for **every** array row against a
       spatial reference holding `ncpl` entries, so it raises
       `IndexError: index 441 is out of bounds` from inside `add_parameters` --
       naming neither the target nor the cause. Measured, not assumed.
     - What ships: `relayer_array_target` calls `make_layered()` + a per-layer
       `set_data` on every declared ARRAY target before externalizing, giving
       porosity the shape `npf_k` already had (`..._layer1.txt` ...). Same numbers,
       and MF6 reads it identically; a test asserts the values survive.
     - **The compromise: this mutates the user's model object.** `_ensure_external_model`
       already calls `set_all_data_external()` on it, so the model was never left
       untouched by `build()` -- but this changes the STORAGE SHAPE of a package the
       user configured, which the previous call did not. Applied unconditionally to
       declared array targets rather than guessed at, because `make_layered()` is
       idempotent (verified on the already-layered `npf_k`) and detecting
       "needs relayering" from FloPy's API is not reliable.
     - Latent bug fixed in passing: a flow model built with a scalar `k` hit the same
       pyEMU `IndexError`. Not just a transport concern.
     - ~~Not done: `style="pilotpoints"` on transport (refused).~~ **That refusal was
       retired 2026-07-30 (ledger 115):** the hardwired `gwf.npf.k` base array was a
       BUG, not a transport limitation -- it also overwrote K33 with horizontal K on a
       flow target. With the base array resolved from the recipe, porosity pilot points
       are correct and supported. Still not done: `dsp.alh` (see 103).

111. **`layers=` was silently ignored, and is now an error (2026-07-29).**
     - `add_native_parameter` filtered per-layer files and then fell through with
       `if selected:` (and `... or files` in `pilot_points.py`) -- so `layers=[0]` on a
       target with no per-layer files, or `layers=[9]` on a 4-layer model, parameterized
       EVERY layer while the control file looked exactly as the caller intended.
     - Now `_select_layer_files` refuses, naming which case it is: a whole-grid array has
       no per-layer file to pick, and a bad layer index lists the ones that exist.
     - **This is a behavior change**, not purely additive: a call that used to "work"
       now raises. Every existing use in the repo passes valid layers, so nothing
       broke -- but a user script relying on the fallthrough would have been relying
       on a bug.

112. **The transport demo now perturbs porosity too, and ships head targets
     (2026-07-29).**
     - `build_canonical_transport_calibration_demo` gained `head_targets` and
       `start_porosity_factor`. Not decoration: measured on the canonical testing
       profile, the concentration responses to `K x3` and `porosity /3` have cosine
       similarity **0.98** (transport velocity is `v = Ki/n`), while heads respond to
       `K x3` by 1.098 ft and to porosity by exactly 0.000000 ft. Estimating both from
       concentration alone is ill-posed; heads break the tie.
     - `start_porosity_factor` is **1.6, not 2.0**: 0.25 -> 0.40 is a substantial error
       that is still a plausible porosity, whereas 0.50 sits at the practical maximum
       for unconsolidated sediment -- a starting value no `physical=` upper bound could
       sit above, so every multiplier above 1.0 would clamp. Found by a test that failed
       for exactly that reason.
     - This CHANGES the demo other tests and notebooks build. The concentration targets
       are unchanged (sampled from the truth run before either perturbation), so
       existing conc-only calibrations still behave as before, but the model handed to
       PEST now starts with the wrong porosity as well as the wrong K.

113. **The relayering of ledger 110 broke the second build on one model, and the
     fix reads a private FloPy attribute (2026-07-29).**
     - Found by the ledger-107 coexistence test, not by inspection.
       `_ensure_external_model` calls `set_all_data_external()` on every build,
       so on build #2 the arrays are EXTERNAL and FloPy refuses `make_layered()`
       (`Converting external file data into layered data currently not
       support`). Relayering unconditionally therefore broke exactly the
       capability ledger 107 was restoring. My "idempotent" claim in 110 was
       tested only on a fresh, never-externalized model.
     - Now skipped when the array is already stored per layer, and when the model
       has a single layer (a whole-grid file there already has `ncpl` values, so
       pyEMU indexes it fine). **FloPy exposes no public flag for "is this stored
       per layer"**, so this reads `array._get_storage_obj().layered` through a
       `getattr` fallback: if that private accessor moves, the fallback attempts
       the relayer, which is the safe direction.
     - `store_internal()` runs first so a model LOADED from a previous run — its
       griddata already external and unlayered — can still be calibrated. That
       leaves the old whole-grid file on disk, unreferenced by the package. So
       `_resolve_files` now matches per-layer files BEFORE the exact name;
       resolving the stale file instead would point every parameter at values
       MODFLOW no longer reads, and the calibration would run happily against
       them. The stale file is left in place rather than deleted — this code
       should not remove files from a user's model directory.

114. **`Run.pest_runs` now scans the whole run directory (2026-07-29).**
     - It scanned `<run workspace>/pest`. With the default at
       `<model workspace>.pest/<name>` and `prepare_run` giving each model its own
       folder, there is no single subdirectory that holds every model's runs.
     - Scanning the run root is a strict superset of what it found before, and
       `find_pest_runs` already skips `*_master` copies, so a run is still listed
       once. The cost is walking a larger tree for one filename.
     - A third copy problem, found while verifying 114: a template is a COPY of
       the model workspace, so a calibration that lived INSIDE it (the old
       default) is copied into every template built afterwards, and discovery
       listed it twice — `review()` on the copy would open a directory nothing
       ever ran in. `find_pest_runs` now skips any build whose metadata is nested
       inside another build's directory, the same reasoning as the existing
       `*_master` skip.

115. **Pilot points scaled the K field whatever target they were pointed at
     (2026-07-30).** Found while scoping §5.8; verified independently before acting.
     - `add_pilot_point_parameter` read `project.model.gwf.npf.k.get_data()` as the
       base array for EVERY target and wrote the interpolated result to the target's
       own file. Measured on the canonical model: `parameterize("k33",
       style="pilotpoints")` moved K33 from 2.863637 to 85.909123 — **30x wrong**,
       vertical anisotropy destroyed — with forward run exit 0 and MODFLOW reporting
       `Normal termination`. Nothing warned.
     - **The ledger-110 guard was too narrow and this is the cost.** It refused
       `recipe.model != "flow"`, which describes the symptom I happened to notice
       (transport) rather than the precondition (does the recipe name its base
       array?). `k33` is a flow target, so it sailed through, and every future flow
       array target would have inherited the bug. The replacement guard states the
       real precondition.
     - Fixing it made the old refusal's stated reason false, so transport pilot
       points are now SUPPORTED, not refused: verified that porosity pilot points
       interpolate against the MST array (0.40, the spoiled truth) and not NPF K.
       A capability gain arriving inside a bug fix, recorded here because the guide
       said the opposite for one day.
     - The forward-run parameters were renamed `k_file`/`base_k` -> `array_file`/
       `base_value`. The names were how the K assumption survived review; they are
       generated into each template at build time, so old templates keep their own
       frozen copy and nothing on disk breaks.

116. **A pilot-points-only calibration could not run (2026-07-30).**
     - `PstFrom.__init__` inserts `apply_list_and_array_pars` into `pre_py_cmds`
       unconditionally (pst_from.py:296-302) but only writes the
       `mult2model_info.csv` manifest it reads when `par_dfs` is non-empty
       (:757-762). Pilot points register through template files after `build_pst`
       and never call `pf.add_parameters`, so every forward run died with
       `FileNotFoundError: mult2model_info.csv` raised from inside pyEMU, naming
       nothing the user had written.
     - Hidden because every example and test pairs pilot points with a second
       parameter. "Calibrate K with pilot points" alone is an ordinary request.
     - An empty manifest is not a fix: `apply_list_and_array_pars` asserts
       `ddf.shape[0] > 0` and has no zero-parameter path. So the command is dropped
       instead, keyed on the same `par_dfs` signal pyEMU itself branches on.
     - **Compromise:** this reaches into `pf.pre_py_cmds`, a public list but one
       whose CONTENTS are pyEMU's. The match is on the substring
       `apply_list_and_array_pars`; if pyEMU renames that helper the filter silently
       stops matching and the bug returns. A test covers the behaviour, not the
       mechanism, so it would catch the regression.

117. **`uzf.vks` ships alone, and UZF needed a new recipe field to be placeable
     (2026-07-30).** Closes §5.8 item 1.
     - UZF is the first list package whose external file does not open with the cell
       identity: it numbers its own cells, so a packagedata row is
       `iuzno layer icell2d landflag ivertcon surfdep vks ...`. Every recipe until now
       assumed `(layer, cell)` at columns (0, 1), hardcoded.
     - **Measured consequence of getting it wrong** (run both ways): with the default,
       pyEMU reads the LAYER as the cell number and all 189 UZF parameters receive cell
       0's coordinates. The `.pst` builds, the multipliers reach the correct rows, and
       MODFLOW terminates normally — only the geostatistics are meaningless, surfacing
       later as `Exception: error inverting cov` from the prior draw, an error naming
       parameters but not the cause. With `index_cols=(1, 2)`: 189 distinct coordinates
       matching the true UZF cell centroids, and a clean prior.
     - **Only vks.** `vks x3` moves heads 1.026 ft on the canonical testing profile;
       `finf x2` moves them 0.011 ft and `surfdep x5` 0.00012 ft. `thts`/`eps` are
       0.995-collinear with each other and ~0.78-0.80 with vks from heads alone, so
       shipping the set would hand users the K/porosity trap of ledger 112 with no
       second data type available to break it.
     - **Not done: UZF PERIODDATA** (`finf`, `pet`, `extdp`). Those rows carry only
       `ifno` — no cellid at all — so pyEMU returns `x`/`y` as `None`, `correlation=`
       dies with a `TypeError` deep in a geostats helper, and `zones=` can never work.
       Adding them needs an explicit refusal for the spatial styles, which is worth
       doing when someone actually wants a recharge-like UZF flux parameter.
     - `use_col=6` is a positional promise about a file MODFLOW writes. Verified that
       enabling `boundnames` appends column 11 and leaves vks at 6 — the caveat an
       adversarial reviewer raised, closed by measurement rather than by reading the
       dfn ordering. A test pins it.
