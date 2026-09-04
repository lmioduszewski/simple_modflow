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
      DIS or raw-DISU model has no `.map()`/`.section()`/`.animate()` rendering from
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
    - `model.conc.section()` and `model.temp.section()` raised `AttributeError: model 'trans' is
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
      `group.conc.section()` in the same pass, so the group tier did not ship with a hole.

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

118. **`style="zone"` was never a working feature, so §5.8 item 2 built the whole
     front door, not just its raster half (2026-07-30).**
     - `zones=` was stored verbatim on the spec and handed to pyEMU as
       `zone_array` unchanged. `zone_array` had exactly ONE occurrence in all of
       `src/`. There was no resolver, no validation, no test and no doc.
     - **The contract is asymmetric between families**, measured on this stack:
       ARRAY targets accept `(ncpl,)`/`(ncpl,1)`/`(1,ncpl)` and reject
       `(nlay,ncpl)`; LIST targets accept only `(nlay,ncpl)` and reject the
       per-cell forms with `IndexError: index 1 is out of bounds for axis 1`,
       raised from inside pyEMU naming neither target nor shape. One zonation
       therefore could not serve both packages — but a zonation is a property of
       the MAP, not of the package, so `resolve_zone_array` normalizes per family
       and the caller never sees the difference.
     - **`docs/myflopy_context.md:150` claimed polygon zones already existed.**
       They did not, for `parameterize`. Corrected in the same pass.
     - **A single-layer list target could not be zoned by pyEMU at all.** The
       correct `(1, ncpl)` shape hits `checker2` (`pst_from.py:2121-2134`), which
       assumes a `(1, n)` array on a vertex grid is an idomain array and reshapes
       it to `(n, 1)`; every lookup then fails exactly like a wrong shape.
       **Compromise:** the resolver pads to `(2, ncpl)`. The padding row addresses
       no real cell so nothing indexes it, but it is a workaround for a pyEMU
       heuristic, not a fix — if pyEMU stops reshaping, the pad becomes harmless
       dead weight rather than breaking.
     - Raster sampling is by **majority vote**, deliberately not the
       area-weighted mean `grid/surfaces.py` uses for surfaces: zone ids are
       labels, and a cell straddling zones 1 and 3 would average to zone 2 — one
       it does not touch and which may not exist. A test asserts the id set is
       `{1, 3}`.
     - **Zone id 0 is left inconsistent, and warns.** pyEMU skips ids below 1 for
       array targets but makes a real adjustable parameter for zone 0 on list
       targets. Normalizing shape cannot normalize that, and silently remapping
       the caller's ids would be worse, so it is called out at resolve time.
     - Not done: overlapping zone polygons take the first match rather than being
       refused — that is a modelling mistake to catch upstream, and guessing an
       ordering here would hide it.

119. **The regularization item was retired and replaced by the gap it was hiding
     (2026-07-30).** Closes §5.8 item 3, not by building what it asked for.
     - **What was asked for would have broken the shipped runner.** Measured on
       PEST++ 5.2.16: `pyemu.helpers.zero_order_tikhonov` flips `pestmode` to
       `regularization` and adds prior-information equations; PESTPP-IES prints
       `prior information equations not supported in ensemble methods, ignoring`;
       and because myflopy writes version-2 control files, the `* regularization`
       keywords land in the control-data keyword block and `pestpp-ies` exits 1
       with `control file parsing error` before a single model run. A
       `cal.regularize()` would have been a loaded gun.
     - `_refuse_regularization_mode` therefore refuses such a control file at
       launch, naming `ies_reg_factor` — IES's own equivalent knob, which already
       works through `run_ies(**pestpp_options)` with no new code.
     - **The real gap: `correlation=` did nothing for pilot points.** It was
       documented at `pilot_points.py:14` from the start and read nowhere in that
       module. `_inject_geostatistical_prior` filtered on `style == "grid"`, with
       the stated reasoning that IDW interpolation supplies spatial smoothness.
       That reasoning is half right — IDW output IS smooth — but the pilot VALUES
       were independent draws, so the field's correlation length came from point
       spacing rather than from the variogram. A 200 m and a 2000 m range gave the
       same prior.
     - Measured after the fix (nearest pilot pair ~372 m on the canonical testing
       profile): range 150 m -> mean prior correlation +0.010 near / +0.006 far;
       range 3000 m -> +0.845 near / +0.495 far. A test pins both directions.
     - **Judgment call: one covariance per LAYER**, not one across all layers.
       Each layer's points are a separate parameter group, and correlating across
       layers would assert a vertical structure the variogram never described --
       `correlation=` is a map-plane range.
     - Not done: a PESTPP-GLM runner. It is the only way regularization proper
       becomes reachable, and it would need its own results-review layer
       (`IesResults` reads IES-specific files). Recorded as Tier B in §5.8's
       remaining item.

120. **Sensitivity ships ensemble-based only; the jacobian half is deferred with
     its cost measured (2026-07-30).** Closes §5.8 item 4 Tier A.
     - pyEMU's `Schur`/`ErrVar` take `jco` as a REQUIRED positional and give the
       textbook answers (CSS, identifiability, parameter contribution, FOSM). They
       are thin to wrap — but a completed PESTPP-IES run produces no jacobian, so
       none of it is reachable from `run_ies`. `pestpp-glm` with `NOPTMAX=-2`
       produces one in `npar_adj + 1` model runs. Tier B is therefore a NEW RUN
       MODE plus a second results-review layer (`IesResults` reads IES-specific
       files), not a wrapper. Deferred, with that cost recorded rather than
       rediscovered.
     - What ships instead reads the ensembles the run already wrote:
       `sensitivity()` / `plot_sensitivity()`, giving `learned`
       (`1 - post_sd/prior_sd` per group, in the parameter's transform space) and,
       with `forecast=`, the ensemble correlation driving that prediction.
     - **The main way this could ship silently wrong is by being read as CSS.**
       It is a global, prior-conditioned measure, not a local derivative; the two
       genuinely disagree when the prior never moved a locally-sensitive
       parameter. Said plainly in the method docstring, the guide, and the
       capability map rather than left for the reader to infer.
     - **Sampling noise is drawn, not just documented.** At `n` realizations a
       correlation near `1/sqrt(n)` is indistinguishable from zero (~0.14 at the
       common `reals=50`). `plot_sensitivity` draws that floor as a line on the
       chart, because the reader is looking at bars, and `n_reals` is returned in
       the frame.
     - **Judgment call: group correlations are means of ABSOLUTE values.** A group
       whose members push a forecast in opposite directions is still influential,
       and a signed mean cancels it to zero -- hiding exactly the group that
       matters. Caught by a mutation test that the first version of the test could
       not see, because that fixture had one parameter per group.
     - Not done: `pyemu.EnDS` data worth. It is the natural next ensemble-based
       addition, but at typical IES ensemble sizes the cross-covariance is
       rank-deficient, and shipping a number that looks like data worth without
       being able to say when it is trustworthy is worse than not shipping it.

121. **The parameterize surface had no completeness test; it does now, and it
     caught a third stale doc while being written (2026-07-30).**
     - Adding a `parameterize` target reaches ~12 code behaviours from one
       `_Recipe` entry but needs prose in four places, and NOTHING failed if you
       missed one. Measured cost in a single day: the `porosity` commit left
       `docs/manual/README.md` on the pre-porosity target set, the `uzf.vks`
       commit left the target out of `parameterize`'s own docstring, and both
       were found by hand rather than by CI.
     - `tests/test_parameterize_payoff.py` is the parameterize twin of
       `test_package_descriptor_payoff.py` and keeps its contract: AUTOMATIC
       sites may grow but never shrink; MANUAL sites must be EXACTLY as
       recorded, so closing one also fails and the win gets written down.
     - **A doc test that could not fail.** The first version searched whole
       files for each target name; deleting `porosity` from the manual's target
       list still passed, because the word appears elsewhere in that file. Found
       by mutating the test rather than by review. It now anchors on the
       enumeration itself (a bullet or a table row) and searches only that.
       Recorded because "the test is green" was, for one draft, worse than no
       test at all.
     - **Compromise: the doc check is per-target TOKENS, not exact names.**
       Short capability tables abbreviate `ghb.cond`/`ghb.bhead` to `ghb`, and
       forcing the long spellings everywhere would make the tables unreadable.
       So `DOC_TOKENS` maps each target to how it is spelled in prose, and is
       asserted to cover the registry exactly — adding a target without a token
       fails, which is the moment the author is meant to go and write the prose.
     - Not covered: that the enumeration is CORRECT, only that every target is
       mentioned in it. A doc could still describe a target wrongly. Catching
       that needs generated prose, which costs the per-target explanation these
       lists exist to carry.

122. **7.3: the 60 exception swallows, what narrowing them cost, and the five
     defects they were hiding (2026-07-31).**
     - **The plan's site list was stale in both directions.** It named six
       files; the real count was **60 handlers across 27 files** -- 49
       `except Exception` and **11 bare `except:`** the plan never mentioned.
       `utils/gdal.py` no longer exists (atticked in Phase 1, as the plan
       predicted it would be).
     - **`_logging.py` was pulled forward from 7.2.** 7.3 asks for a
       `logger.debug` on every narrowed handler, and 7.2 owns the logger. Doing
       7.3 first without it would have meant a second pass over all 60 sites.
       7.2 shrinks to the `print()` replacement, which is what it was really
       about. `_logging` is graph-EXTERNAL in the import-layer derivation, like
       `_vendor`: it is a stdlib-only leaf that modules at every depth import,
       so counting it would push every current L0 leaf to L1 and say nothing
       about the architecture.
     - **Ten handlers stay `except Exception` ON PURPOSE**, each with the
       measurement that justifies it, and `tests/test_exception_narrowness.py`
       pins the list exactly in both directions. The recurring reason is that
       the raisable set is genuinely open: flopy's `utils/voronoi.py` contains a
       literal `raise Exception(...)`; `pickle.dump(model)` walks an unbounded
       third-party graph; `grid_spec_resolver` does a `setattr` on an object a
       USER'S builder script returned (pydantic raises a plain `ValueError`
       there); `__dir__`/`__repr__` must never raise. This is the honest number:
       "narrow all 60" was not achievable, and pretending otherwise would have
       traded working fallbacks for crashes.
     - **Several tuples came out WIDER than a first pass would write**, and the
       reason generalizes: most of these `try` blocks wrap a PIPELINE, not a
       call. On a file-backed model `model.gwf`/`model.sim`/`model.vor` are LAZY
       properties, so the first access performs a whole `MFSimulation.load` --
       dragging in flopy's `MFDataException`/`FlopyException` (which subclass
       `Exception` directly, so no builtin tuple reaches them), a
       `UnicodeDecodeError` on a binary `mfsim.nam`, and myflopy's own
       post-load arithmetic. Three MRO facts were measured rather than assumed:
       `rasterio.errors.RasterioError` is NOT an `OSError` (only
       `RasterioIOError` is, by multiple inheritance), pyproj's `CRSError` is a
       `RuntimeError`, and `numpy.ma.MaskError` subclasses `Exception` directly.
     - **Five defects the broad handlers were hiding**, each now pinned by a
       test in `tests/test_narrowed_swallows.py` or `tests/test_model_diff.py`:
       1. **A group member without the package hid every other member's
          differences.** `GroupPackageInputs.get()` builds the WHOLE group
          before `compare()` filters to one model, unguarded -- so in a
          three-model group where one member has no GHB, the builder raised,
          `_value_cells_changed` caught it and returned 0, and `report()`
          printed "identical to reference" for a model whose conductances
          really did differ. `model.diff(a, b)` is a documented door, so a 3+
          model group is ordinary usage. HIGHEST severity of the five.
       2. **`_model_nlay` reported a cell count as a layer count** --
          `np.asarray(botm).reshape(-1).size` is the total number of botm
          ENTRIES. On a 3-layer, 4-cell grid it answered 12.
       3. **`read_gpkg` caught its own raise.** The bare `except:` swallowed the
          `TypeError: Unexpected geometry type` raised four lines above it,
          silently returning a PARTIAL feature set. Its comment described
          terminating a `while` loop that no longer existed.
       4. **`contour_line_segments` caught its own raise too** -- the same shape,
          found independently: a misspelled `method=` drew no contours and said
          nothing. Validation now runs before the `return []` early exits, so
          the error does not depend on whether the data would have contoured.
       5. **`_resample_timeseries_df` did `return print(...)`** -- printed a
          message and handed the caller None, so the real symptom arrived later
          as an `AttributeError` on None, far from the bad argument.
     - **`drn.py` now RAISES where it used to guess.** A failed `gdf_topbtm`
       lookup fell back to `drain_elevation = bottom_addition` -- a drain at the
       bare offset (often 0) instead of at the layer bottom, which drains the
       aquifer. There is no defensible default: the cell not being in the
       layer-surface frame means the grid and the drain footprint disagree, and
       MODFLOW would have run happily on the wrong answer. This is a **behaviour
       change**: a build that silently produced bad drains now fails with a
       message naming the cell and pointing at `bottoms=`.
     - **Compromise: one diff asymmetry is documented, not fixed.**
       `_package_present` matches by PREFIX, so an array-form RCHA/EVTA on an
       externally loaded model reports as present while the cell tier -- which
       cannot read array-form packages at all -- yields no keys; two such models
       compare as identical. Telling that apart needs a third state ("present,
       unreadable by this tier") threaded through the summary, not a wider
       `except`. Recorded in a comment at the site and deferred.
     - **Compromise: the mover-flow probe gets no log line.** It tries 7
       packages x 2 directions and most are EXPECTED to fail, so one debug line
       per absent term would be 14 lines of noise per call. Skipping is the
       documented behaviour there, not a degradation.
     - **A debug line introduced its own bug**, worth recording because it is
       the standing risk of this whole exercise: the new `logger.debug` in
       `_value_cells_changed` referenced `self._package_name`, which that class
       does not have. It only executes when the handler fires, so the whole
       suite stayed green -- it surfaced only when the fix was MUTATED to
       confirm the new test could fail. Log lines added to rarely-taken branches
       are untested code by construction.
     - **The ratchet keys on `path:function`, not line number.** Two broad
       handlers in one function collapse to one key, so a second one there would
       not fail. Deliberate: keying by line would turn every unrelated edit into
       a test failure, which is how ratchets end up deleted.
     - **`docs/myflopy_context.md` deliberately NOT updated.** It is the
       code-derived map of the model-building *capability* surface;
       `myflopy._logging` is private and narrowing an exception handler adds no
       capability. The convention lives in `CLAUDE.md` and the module docstring
       instead. Recorded so the omission reads as a decision rather than a miss.

123. **7.2 (prints → logging), 7.4 (docs debt), and the flopy deprecations —
     what was exempted and why (2026-08-01).**
     - **The count was 60, not 61.** A line-grep for `print(` over-counts:
       several hits are docstring `>>> print(...)` examples, and one was the
       word "foot**print(s)**". The AST count is the real one, and the same
       AST walk is what `tests/test_no_library_prints.py` uses.
     - **The exemption rule is the user's, stated 2026-08-01:** *anything whose
       job is to print a report a human asked for* keeps printing. Operationally
       that resolves to **output behind a flag the caller passed** —
       `verbose=True`, `progress=…`, `verbosity_level>0` — because the flag IS
       the human asking. 13 prints qualify (7 in `triangle.py`'s CVT optimizer,
       plus `ghb`, `recharge`, `calibration` x3, `run_model._progress`,
       `interactive_plotting._emit_progress`). They are exempt by AST inspection
       rather than by being listed, so adding a `verbose=`-gated print later
       needs no test edit.
     - **Exactly one UNCONDITIONAL print survives**: `runtime.run_simulation`'s
       `"Success is: …"`. It is the interactive "did my model run" answer,
       printed directly beneath flopy's own `silent=False` MF6 output — a
       success line that appeared only if you had configured logging first would
       be a worse answer to the question the caller just asked. The chatter
       AROUND it (saving the `.model` snapshot, and failing to) is a side effect
       and now logs.
     - **Level policy**: `info` for progress on slow work (sampling a raster,
       building the grid, reading a GeoPackage), `warning` for a degraded or
       ignored input (unrecognized colorscale, a location intersecting no cell,
       reprojecting mismatched CRS), `debug` for per-row or per-frame detail.
     - **Two `return print(...)` calls became raises** (`mf2Dplots`): they
       printed a message and handed the caller None, so the real symptom arrived
       later as an `AttributeError` on None. Same defect class as 7.3's
       `_resample_timeseries_df`. Two more prints that were immediately followed
       by a bare `raise ValueError` had their message moved INTO the exception.
     - **One leftover debugging print was deleted**, not converted:
       `calibration.py`'s `print(len(lak_obs['Deep Lake']))`. It named a specific
       lake from someone's model.
     - **Compromise: `InterpolatedSurface.interpolator` still silently ignores an
       unsupported value.** Its docstring documents "ignored with a warning", so
       7.2 made the warning real (a `logger.warning`) rather than changing the
       contract. The trap behind it is separate and NOT fixed: after the warning
       the code falls through to the `use_rbf` branch, so a typo'd interpolator
       returns an rbf surface rather than the requested one. The same is true of
       the `kstpkper` setter, which keeps the previous value. Both are behaviour
       changes beyond 7.2's scope; recorded here so the next reader does not
       assume the warning is the whole story.
     - **7.4 found 120 dead links.** `docs/codebase_structure.md` wrote every
       link as an absolute Windows path (`C:/Users/lukem/Python/Projects/...`),
       so the repo's own structure guide was fully broken for every reader on
       every platform. Rewritten as repo-relative. The tracked
       `myflopy_api_pamphlet.pdf` is now untracked and gitignored — it is fully
       regenerated by its script, so tracking it made every rebuild a binary
       diff and had made the generator effectively unrunnable all week.
     - **Compromise: the manual's 21 unwritten chapters are ALLOWLISTED, not
       fixed.** `docs/manual/README.md` is a complete table of contents for a
       manual whose text is chapters 3 and 4. Failing on those links would just
       get the test deleted; listing them makes the promise auditable, and a
       second test fails if a listed chapter is later written and left in the
       list. A typo'd chapter name still fails, because it will not be listed.
     - **Compromise: the link check ignores anchors and external URLs.** Anchors
       (`#section`) would mean slugifying every heading the way GitHub does;
       external URLs need the network and rot for reasons outside this repo.
     - **The flopy deprecation fix needed no `_flopy_compat.py`.** 7.1 proposed
       centralizing every flopy touchpoint; for `package_names`/`package_name_dict`
       the narrow version was free, because `SimulationBase.package_names`
       already read the supported `get_package_list()`. Three call sites simply
       stopped reaching past it. Worth noting before 7.1 is revisited: the
       payoff there may be smaller than the plan assumes.
     - **A "defensive" getattr still fires a DeprecationWarning.**
       `getattr(gwf, "package_name_dict", {})` looks safe but warns anyway —
       READING the attribute is what warns, and the default only covers absence.
       `tests/test_no_deprecated_flopy_calls.py` matches the AST attribute
       access, so it catches the getattr spelling too.
     - **That guard is deliberately a NAMED LIST, not a blanket rule.** A general
       "no deprecated calls" is unenforceable against dependencies we do not
       control, and `warnings.simplefilter("error")` would fail on
       pandas/geopandas/pyemu deprecations that are not ours to fix.

124. **Phase 8 was retracted and rewritten: the goal was a vocabulary, not a
     directory (2026-08-01).**
     - The original §8 was "move nine standalone plotting modules into
       `myflopy/plot/`". A 12-agent re-verification against current code found
       **six of its nine import-graph rows stale** (Phases 4-7 had moved the
       graph under it), and three structural problems the plan could not survive:
       1. **The move order was not leaf-first.** Measured by import closure:
          `contour_plotting`, `cross_section_plotting`, `budget_plotting` and
          `mf2Dplots` pull nothing; `xsections` pulls FIVE of the eight other
          movers — and it was scheduled third, ahead of `grid`, `budget` and
          `heads`, which it depends on.
       2. **Its acceptance test was impossible.** "Existing plotting tests must
          pass unedited" cannot hold, because `pyproject.toml:158-164` escalates
          `myflopy` DeprecationWarnings to errors, so the D12 facades the plan
          prescribes break every test importing an old path — 11 files, not the
          five it named (two of which import no moving module at all).
       3. **Its scope boundary was already violated.** §8 says the `package_*`
          explorer family stays put, but `headsplus.py:50` imports
          `package_plotting` at module level; importing `xsections` loads 37
          myflopy modules including `package_plotting` and `hover`.
     - **A third hazard nothing had reported:** `tests/test_interactive_plotting.py`
       monkeypatches by STRING (`monkeypatch.setattr("myflopy.modflow.mf6.
       interactive_plotting.plot_model_head_map", …)`, three sites). After a move
       the patch lands in the stub's namespace while the function under test
       resolves the name from its own globals — the patch silently becomes a
       no-op and the test passes while testing nothing.
     - **The cost/benefit was upside-down.** 5,518 lines would move; 7,367 lines
       of plotting would stay (`package_surface_water` 2256, `package_plotting`
       1990, `prt_maps` 1237, `hover` 824, `interpolated_surface` 525, `viz` 451)
       plus `pest/ies.py`'s twelve plot methods. `__compatibility__` would go
       23 → ~49 entries, frozen for two tagged releases when the repo has ONE tag.
       And the user-visible surface would not change at all.
     - **The target layout was also incomplete**: three plotting modules created
       after the plan's 2026-07-14 verification (`grid/interpolated_surface.py`,
       `prt_maps.py`, the plotting half of `pest/ies.py`) appear nowhere in it.
     - **What replaced it.** The user's stated goal was discoverability, so §8 is
       now a VOCABULARY phase: three layers (pictures / composition / output),
       four picture verbs chosen by GEOMETRY (`map`, `section`, `surface`,
       `timeseries`), with content and renderer as options. Files do not move.
     - **Two collapses the user identified that the code already supported.**
       `contours` and `pathlines` are not picture kinds — `plot_cell_contours`
       and `plot_particle_pathlines` both take `ax=`, i.e. they draw onto an
       existing map, and `Choro` already has eight content parameters and five
       `add_*` methods. And `mosaic`/`animate` are combinators over arbitrary
       pictures, not map verbs — `viz.mosaic`'s docstring already says "Compose
       arbitrary panel objects", it was just shadowed by `<node>.mosaic()` sugar.
     - **VTK is a BACKEND of `surface`, not a picture kind.** The 3-D code splits
       along renderer lines, not content: plotly (`InterpolatedSurface.
       surface_trace`, `surface_3d`, `plot3d`) vs VTK/PyVista (`vtk_3d`,
       `ParticleTrackingScene`, `export_particle_tracking_html`). The precedent
       already exists as `map(backend="plotly"|"mpl")`.
     - **Compromise: no deprecations at all.** The user confirmed there are no
       other consumers and that old NOTEBOOKS need not keep working, only old
       MODELS need to keep opening. Verified they do: no pickled object caches a
       plotting class (`workspace.py` pickles grids and array-bearing packages;
       no `Choro`/`GridSection`/`XSection` is cached on either). So each stage
       deletes what it replaces. This is a deliberate departure from the
       deprecation policy, justified by a single-user codebase at one tag.
     - **Compromise: the ~35 internal `build_*_table`/`build_*_payload` helpers
       keep their prefix.** Only three `build_*` names were ever public plotting
       API. Renaming the rest is churn with no user-visible gain.
     - **Known long pole: `animate`.** Every other stage renames or collapses
       existing behaviour; a frame-accepting `animate(frames)` is new code,
       because today's `Animation(model, periods)` is model-bound and redraws
       from the model rather than composing pictures.
     - **Notebooks are edited IN the stage that breaks them, not batched at the
       end (user's call, 2026-08-01).** 11 notebooks call the changing API. The
       first draft deferred all of them to 8.7; the user runs them live, so that
       would have left a broken notebook standing for days. Each stage now fixes
       what it breaks in the same commit, and re-executes the canonical set via
       `scripts/render_notebooks.py` -- a source edit that leaves a notebook
       unrunnable is exactly the failure a grep cannot see.
     - **Compromise: two notebooks are flagged, not edited.**
       `bearcreek_uncertainty.ipynb` and `demo_ies_uncertainty.ipynb` are
       UNTRACKED working files -- the user's uncommitted work. Between them they
       hold 2 `build_choropleth` calls and 8 `.plot()` calls that this phase
       breaks. They get a written list of what to change and an explicit ask;
       rewriting someone's in-progress analysis is not ours to do. Tracked
       canonical notebooks are edited surgically (cell-level, outputs intact) --
       never checked out or regenerated wholesale.

125. **8.1: one contract for every picture, and two things the survey got wrong
     (2026-08-18).**
     - `viz.Picture` gives `Choro`, `XSection`, `GridSection` and
       `InterpolatedSurface` one shape: `.fig`, `.show()`, `.save(path)`,
       `.html(path)`, and `_repr_mimebundle_` so a picture renders itself. Before
       it there were four names for the figure (`.choropleth`, `.fig`,
       `.figure`, and none) and THREE incompatible meanings for `.plot()`:
       return the figure (`Choro`), show it and return None (`GridSection`), or
       open a browser window and return None (`InterpolatedSurface`).
     - **Correction to the scoping discussion.** I told the user `Choro` had
       three spellings of "give me the figure". It had two (`.plot()` and
       `.choropleth`). `get_choropleth()` returns the choropleth **trace**, which
       `viz.mosaic` composes panels out of -- deleting it as a duplicate would
       have broken every mosaic. It survives, with a test saying why.
     - **A live bug the contract had to fix to exist.** `add_choropleth()` calls
       `fig.add_trace(...)` unconditionally, so the old `.choropleth` property
       was NOT idempotent: touching the figure twice drew every trace twice, and
       `.plot()` just returned it. `.fig` now assembles once and caches, which is
       what makes `picture.fig.update_layout(...)` then `picture.show()` safe.
     - **`InterpolatedSurface` was violating the house-figure rule.** Its
       `plot()` built a BARE `go.Figure` and called `.show(renderer="browser")`,
       so it returned None, could not be embedded, forced a browser window, and
       dropped the template, `scrollZoom` and `dragmode="pan"` that
       CLAUDE.md requires. Replaced by a `.fig` built on `viz.Fig`.
     - **Compromise: `save()` raises its kaleido ImportError inline instead of
       through `myflopy._optional.require`.** `viz` is deliberately
       externals-only (Layer 0), and importing any myflopy module there would
       either push every current L0 leaf up a layer (module-level) or trip the
       deferred-import ratchet, which only ever moves DOWN (function-level). The
       duplicated message is the cheaper of the two.
     - **`_assembled` is a CLASS attribute, not an `__init__` assignment.**
       `Choro` is constructed via `object.__new__` by several test doubles; a
       `.fig` that raises AttributeError on those is a contract that only
       half-holds. Found by the suite, not by review.
     - **A grep-based sweep missed a production call site**, exactly as the plan
       warned. `pest/ies.py:1686` and `:2234` call `.plot()` on a variable named
       `choro`, but my first regex keyed on names containing "choro"/"map" and
       still missed them because of the surrounding expression. Receivers were
       then resolved by reading every `.plot()` in `src/`, which is what the plan
       says to do and what I should have done first.
     - **Notebooks: 6 tracked files edited surgically** (12 edits, text-level, so
       nothing but the changed strings moved and all outputs stayed intact).
       `bearcreek_uncertainty.ipynb` and `demo_ies_uncertainty.ipynb` are
       UNTRACKED working files and were deliberately NOT touched -- they hold 6
       picture `.plot()` calls between them, listed for the user to apply.
     - **Left for 8.2 on purpose:** `ies.forecast(name).plot()` in two notebooks
       and `self.forecast(name).plot()` in `ies.py:2318` are the grammar's
       TIMESERIES verb, not a picture accessor. `model.gwf.plot()` and
       `model.gwf.chd.plot()` are flopy's own and are never ours.
     - **Collision found for 8.4:** `vor.plot` already exists and is **flopy's**
       `VoronoiGrid.plot` (verified via `__qualname__`), called as
       `model.vor.plot()` in two notebooks. Making `vor.plot` a namespace would
       shadow it. Recorded in the plan at 8.4 rather than decided here.

126. **8.2: `xs` became `section`; `plot` stayed, because the rename it was
     given rested on a false premise (2026-08-18).**
     - `xs` → `section` across 25 call sites and 4 definitions, plus the `kind=`
       string values. `"section"` was ALREADY a half-accepted alias
       (`package_plotting.py` matched `("xs", "section")`), so half the plumbing
       was in place; it is now the only spelling.
     - **`xs` means two different things and only one is the verb.**
       `XSection.xs` and `InterpolatedSurface.xs` are x-COORDINATE properties --
       a blind rename would have corrupted the geometry code. The separating
       rule, verified against every occurrence in the tree: `.xs(` with a paren
       is always the verb; `.xs` without a paren is always coordinates.
     - **`plot` → `timeseries` was planned, scoped, and DROPPED on evidence.**
       An inventory of every `plot()` on a grammar node found that it does not
       draw time series: distance profiles (`SfrProfileView`,
       `SfrReachProfileView`), bar charts (`LakBudgetView`, `PRTEndpointsView`,
       `PRTCaptureView`), a cumulative arrival curve (`PRTTravelTimeView`), a
       histogram (`IesForecast`). Renaming would have produced
       `endpoints.timeseries()` returning a bar chart of cells.
     - **The defect was the DOCUMENTATION, not the code.**
       `docs/view_layer_conventions.md` and `myflopy_context.md` both described
       `plot` as the series verb, and the code never honoured it. That stale
       sentence is what led the user to ask for the rename in the first place --
       a doc error that nearly became a 91-site API change. The doc now defines
       `plot()` as *the node's non-spatial chart*, lists what shape each node
       actually returns, and states the corollary: a `plot()` that returns a MAP
       is misnamed.
     - **Compromise: one verb whose name does not change with the chart type.**
       Splitting into `timeseries`/`profile`/`bars`/`histogram` would be more
       accurate per node and was rejected: it breaks "every noun answers the same
       verbs", which is the property that makes the grammar guessable without
       reading source. `chart` was considered as a more literal name and rejected
       as ~91 call sites for a marginal gain, now that 8.1 has removed the
       `Choro.plot()` collision that made `plot` ambiguous.
     - `RchInput.plot`/`UzfInput.plot`/`DrnInput.plot` → `map()`: they return
       CHOROPLETHS, so they were misnamed under every option. No callers existed.
     - The user's two UNTRACKED notebooks (`bearcreek_uncertainty`,
       `demo_ies_uncertainty`) were fixed on request in this pass -- 4 edits,
       picture receivers only; `model.vor.plot()` (flopy's own) and
       `ies.forecast(...).plot()` (the grammar verb, which survives) left alone.

127. **8.3: `myflopy.plot`, and the legacy choropleth layer removed (2026-08-18).**
     - **8.3a** added `src/myflopy/plot/` with four verbs chosen by GEOMETRY --
       `map` (plan view), `section` (vertical slice), `surface` (3-D) -- plus
       `mosaic` (re-exported from `viz`, not reimplemented) and `animate`. That
       one rule retired three names on its own: `plot3d` is a `surface`,
       `map_nodes`/`plot2d` are `map(values=...)`, and contours are an option on
       `map`, not a verb. A test pins each retirement.
     - **`mf.plot` needed a line in the root `__getattr__`**, which special-cases
       subpackages one at a time. Without it `import myflopy.plot` works while
       `mf.plot` raises -- a half-wired export that no import test would catch.
     - **No `timeseries` at module scope, deliberately.** A chart belongs to a
       node that knows the model's periods; there is no useful "chart these bare
       arrays" that plotly does not already do better.
     - **8.3b deleted the legacy choropleth layer**, which ran deeper than the
       plan's "five entry points": the whole `heads_plotting.py` module (its four
       public names were all compatibility wrappers), the three `HeadsPlus`
       methods that fronted them (`choropleth`, `plot_choropleth`, `plot_heads`),
       `mf2Dplots.ChoroplethPlot` (orphaned -- its only caller was
       `heads_plotting.choropleth`), and `plot_drn_choropleth` with its
       `budget.plot_choro` wrapper. Verified beforehand: NO external callers.
     - **The user drew the line explicitly (2026-08-18):** remove the legacy trio
       from `model.hds`, keep the grammar verbs (`map`/`section`/`plot`/`mosaic`/
       `animate`). `hds` is for reading and querying head data; the pictures come
       from the grammar or from `myflopy.plot`.
     - **Compromise: a small capability delta on `plot_heads`.** It accepted a
       PATH to a locations file plus a `crs=`, resolving features to cells
       itself. Its replacement, `model.hds.plot(cells=...)`, takes cell ids. The
       feature-to-cell resolution still exists (`vor.get_vor_cells_as_series`,
       `datatypes/locs.py`) but the caller now does that step. Recorded rather
       than hidden, because "deleted, replaced by" is not quite true here.
     - **`build_choropleth`/`build_grid_section` became PRIVATE factories**
       (`_choropleth_factory`/`_grid_section_factory`) rather than being deleted:
       they are the implementation behind `plot.map`/`plot.section` and the
       `vor.choropleth`/`vor.cross_section` aliases that 8.4 will fold into
       `vor.plot`. Dropped from `grid/__init__.__all__`.
     - **A blanket rename hit a test's EXPECTED STRING.** The `build_choropleth`
       -> `_choropleth_factory` sweep rewrote an assertion in
       `test_master_visualization_prt_example.py` that describes what a NOTEBOOK
       must contain -- so the test started demanding a private name inside a
       notebook. Caught by the suite. The lesson is the same one 8.1 recorded:
       a rename over `tests/` can change what a test is asserting ABOUT, not just
       how it spells it.
     - **All four notebooks migrated off the deep import.**
       `from myflopy.modflow.mf6.grid.plotting import build_choropleth` ->
       `from myflopy import plot`, and the calls to `plot.map(...)`. This was the
       stated goal of the whole phase: nothing outside the package names an
       internal module path any more.
     - **Notebooks are not executed by the suite**, so a migration that reads
       fine and dies on first run is the live failure mode. The exact rewritten
       call shapes are now pinned by a canonical-model test. That test caught one
       real error: a comment claiming `plot.map(..., backend="mpl")`, which the
       signature does not accept -- the mpl backend is `.plot_mpl()`.

128. **8.4a: `model.plot`, and the 25-site workaround it let us delete (2026-08-18).**
     The stage as planned was "bind the verbs, delete `model.cor`/`section`/`srf`".
     Measuring first turned up more than a rename, and three things the plan had
     wrong.
     - **`SimulationBase.plot` ALREADY EXISTED** (`base.py:1419`), a one-line
       delegate to FloPy's `MFSimulation.plot`. The plan flagged only the
       `vor.plot` collision, and the user's shadowing approval was given for that
       one. Shadowing it costs nothing measurable: **zero callers** anywhere in
       `src/`, `tests/`, `docs/` or `examples/`, and the renderer is still
       reachable as `model.sim.plot(...)`, which is where it actually lives.
       `patch_simulation_plot` (`run_model.py`) patches `sim.plot`, not
       `model.plot`, so it is unaffected.
     - **The real find: `show_layer_elevs`.** `_choropleth_factory` hardcoded
       `False` while `cor` defaulted `True`, so from 8.3 until this commit
       `plot.map(model)` silently dropped five hover rows ('Top of Model' plus
       each 'Layer N Bottom') that `model.cor()` showed. The codebase was already
       working around it: **25 call sites** in 11 files repeated
       `kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(model))`,
       and that helper was **duplicated verbatim** in `project/group/_shared.py`
       and `package_explorer_utils.py`. The rule now lives once, in the factory,
       resolved from `vor.gdf_topbtm` -- the exact condition both copies tested.
       All 25 rituals, both helper copies, and their re-exports are gone. Pinned
       in both directions (model keeps the hover; a grid without `gdf_topbtm`
       still opts out) and mutation-tested.
     - **~90 lines of signature restatement deleted.** `SimulationBase.cor` (31
       params) -> `accessors.build_choro` (31 params) -> `Choro.__init__` was
       three hand-maintained copies of one signature; `section` -> `build_xsection`
       -> `XSection` was the same shape in 15. Both middle layers had exactly one
       caller. **They had already drifted**, which is the argument against them:
       `build_xsection` disagreed with `XSection.__init__` on `use_rbf` (False vs
       True) and `x_or_y` (None vs `'x'`), and four `XSection` parameters
       (`surf_type`, `interpolator`, `section_name`, `clip`) were simply
       unreachable through `model.section()`.
     - **COMPROMISE — the drifted defaults resolve toward the class.**
       `model.plot.section()` uses `XSection.__init__`'s own defaults, so
       `use_rbf` is now True and `x_or_y` None where `model.section()` gave False
       and `'x'`. Not preserved, because the class is the source of truth and the
       restatement was the bug; the blast radius is one internal caller
       (`interactive_plotting.py:1317`) and zero notebooks. Same call for
       `model.plot.surface()`, which uses `InterpolatedSurface`'s
       `kstpkper[0]` rather than `ModelSurface.hds()`'s `(nstp[0]-1, 0)` --
       identical on a steady-state model, different on a transient one.
     - **COMPROMISE — discoverability, not behaviour, is what the collapse
       costs.** `plot.map`'s signature is `(source, /, values=None, **kwargs)`, so
       `help()` and IDE completion no longer advertise the 15 contour/hillshade/
       bounds parameters `model.cor` spelled out. Verified they all still ARRIVE:
       each is named in `Choro.__init__` with an identical default and rides
       `**choro_kwargs` untouched. Deferred rather than fixed -- naming them on
       `_choropleth_factory` would rebuild the restatement this entry deletes.
     - **`ModelSurface` was already dead.** Its only apparent caller
       (`xsections.py:273`) was inside a string used as a commented-out block, so
       `surface_data.py` went with `model.srf` rather than surviving as an
       internal helper. That also removes its `plot: bool = False` flag, a
       Layer-3 concern 8.1 had left embedded in a Layer-1 constructor.
     - **Two monkeypatches had to move, not just be renamed.** Tests that patched
       `model.cor`/`model.section` now patch `ModelPlots.map`/`.section` on the
       CLASS: `model.plot` returns a fresh namespace per access (deliberately --
       `from_built_run` rehydrates via `cls.__new__`, so a memoized namespace
       could outlive its state), and patching an attribute of a throwaway object
       does nothing.
     - **Import layering cost: zero.** `simulation.base` is L8 and `myflopy.plot`
       is L5, so the namespace imports downward at module level; the derived layer
       map is unchanged and `deferred_total` stays at 60. This is why the GRID
       namespace (8.4b) must be built in `grid/plotting.py` (L3) instead --
       `voronoi.py` is L4, and reaching for `myflopy.plot` there would point
       upward and force a deferred import, which the exact-match ratchet forbids.

129. **8.4b: the eleven grid aliases, and what only five of them were (2026-08-18).**
     `vor.plot` is now a namespace with three verbs -- `map`, `section`, `grid`
     -- **deliberately shadowing FloPy's inherited `VoronoiGrid.plot`**, which
     the user approved. Reading the eleven aliases before deleting them changed
     what the stage was:
     - **Only five were distinct pictures.** `choropleth`/`show` were the same
       Choro (`show` returned it rather than showing anything); `cross_section`
       is now `plot.section`; `plot2d` is now `plot.grid`;
       `show_selected_cells` became the `select=` option on `map`.
     - **`show_overlapping_geometry` folded into that same option.** It was a
       two-liner over `get_vor_cells_as_series` + `show_selected_cells`, so
       `map(select=<gpkg or geometry>)` now accepts either cell indices or
       something to intersect. One option, both spellings.
     - **`mapit` was dead code** -- it calls GeoPandas `.explore()`, which needs
       `folium`/`mapclassify`; neither is installed nor declared in
       `pyproject.toml`, so every call raised. **`plottri` is not a picture of
       the Voronoi grid at all** -- it draws the *triangulation* the grid was
       derived from, returns None, opens a browser, and raises on any grid not
       built from a Triangle mesh. Both deleted outright.
     - **`dash_selector` is a LIVE user workflow, and only the alias went.**
       (User, after the fact: it is how they box/lasso-select Voronoi cells and
       get the cell-id list back, for querying a batch of cells quickly.) The
       capability is untouched -- `Choro.dash_selector()` is where it always
       lived; `vor.dash_selector` was a property wrapping a default Choro. The
       spelling is now `vor.plot.map().dash_selector()`, which also accepts the
       map's own options (`values=`, `zmin=`, ...) instead of only the default.
       Being a property was in fact a bug for this use: merely tab-completing or
       `repr`-ing `vor.dash_selector` launched a Dash server on port 8050. It
       only starts when called now.
     - **`plot3d` is retired with no replacement in this stage.** It drew 3-D
       Voronoi edges as a bare `go.Figure`. `plot.surface` is the 3-D verb, but
       it interpolates a field rather than drawing mesh geometry, so this is a
       genuine (small) capability gap until 8.5 does the 3-D work. Zero callers
       anywhere, which is why it is recorded rather than rebuilt.
     - **COMPROMISE — `grid` is a fifth verb, not a `map` option.** It breaks the
       "geometry chooses the verb" rule on its face: a mesh view IS a plan view.
       It earns the exception on a measured constraint -- **`Choro` hard-requires
       a CRS** (a naive-geometry grid raises "Cannot transform naive geometries"),
       because its basemap does, while `plot2d` and FloPy's `plot` are both
       CRS-free. Deleting both without a CRS-free verb would remove the only view
       available for a grid still being refined, before a projection exists.
     - **The namespace defines `__call__`.** `vor.plot()` -> `vor.plot.grid()`,
       so the two UNTRACKED notebooks that call `model.vor.plot()` keep working
       untouched -- they use it purely as "show me the grid", which is exactly
       what `grid()` draws. That meant this stage needed **no edits to the user's
       work-in-progress files at all**. The return type does change (a plotly
       Picture, not an mpl Axes); `plot.grid().plot_mpl()` is FloPy's renderer,
       called unbound, and is pinned by a test that reads its source.
     - **`vor.choropleth` was a FOURTH copy of the choropleth signature** (after
       `cor`, `build_choro`, and the factory), and its own docstring admitted it:
       "A thin restatement ... keep the two in step". Gone with the rest.
     - **A migration bug fixed in passing.** `inputs.py:156` called
       `vor.show_selected_cells(drn_cells)` and **discarded the result** -- a
       method named `map()` that returned None. It now returns the picture.
     - **`grid/__init__.py` exported nine of these names with zero consumers**
       repo-wide. Dropped.
     - **Layer cost: still zero.** `GridPlots` lives in `grid/plotting.py` (L3),
       which `voronoi.py` (L4) already imported. Building it in `myflopy.plot`
       (L5) instead would have pointed upward and forced a deferred import; the
       derived layer map is byte-identical and `deferred_total` stays at 60.
       A test now pins that ordering so the next person does not "tidy" it.
     - **`import myflopy` is not a syntax check.** Removing the last name from
       `voronoi.py`'s `if TYPE_CHECKING:` block left it empty -- a hard
       `IndentationError` -- and `python -c "import myflopy"` still printed "ok",
       because the top-level package resolves subpackages lazily via `__getattr__`
       and never touched the file. `ruff` caught it as `invalid-syntax`. The habit
       worth keeping: after deleting imports, lint the FILES rather than importing
       the package, and diff the lint output against the pre-change baseline so
       one new error is not lost among the pre-existing ones.

130. **`Choro.ani` broke the idempotence the Picture contract requires (2026-08-20).**
     Found while scoping 8.5, fixed on the spot because it is one line and live.
     `.ani` builds the frames figure and assigns it to `_fig`, but left
     `_assembled` False -- so the next `.fig` access re-ran the assembly and
     appended the static choropleth, contours, locs and overlays **on top of the
     animation**. Since `.show()`, `.html()` and `.save()` all route through
     `.fig`, `choro.ani` followed by any of them produced a doubled figure.
     - **Why nothing caught it.** The shipping exporters
       (`ModelVisualization.plotly_*_animation`) read `.ani` and never touch
       `.fig`, so the whole tested path missed it. It is the *natural* new
       spelling -- `model.plot.map(...).ani` then `.show()` -- that was broken,
       which is exactly what 8.6 will make canonical.
     - This is the same defect class 8.1 fixed for the static figure (ledger
       125): assembly that mutates rather than rebuilds. The `Picture` docstring
       already required `fig` to be idempotent; `.ani` was simply not honoring it.
     - Pinned by `test_an_animation_is_the_choros_figure_from_then_on` and
       mutation-tested. **8.6 still has to decide what `.ani` IS** -- it returns a
       bare `viz.Fig` rather than a Picture, so an animation cannot answer
       `.html(path)` today. That is the redesign, not this fix.

131. **`plot._grid_of` read `.vor` as "this is a model" (2026-08-20).**
     Found while scoping 8.5a; fixed before `stack.plot` could depend on it.
     A `LayerStack` and a `LayerBuildResult` both carry `.vor`, so the duck-typed
     dispatcher classified them as MODELS and sent them down the results path,
     where they died on `.hds` (verified: `AttributeError: ... no attribute
     'hds'`). The discriminator is now RESULTS -- `hds`/`conc`/`temp` -- because
     that is what actually distinguishes a model from geometry that happens to
     know its grid.
     - The 8.4 sweep flagged this as a hazard for the model namespace and
       `ModelPlots` sidesteps it by binding its subject explicitly (ledger 128).
       That did not help the FREE functions, which is where a layer stack enters.
     - **NOTE for 8.5a's `stack.plot.map()`:** delegating to `plot.map(vor,
       values=thickness)` returns a `Choro`, which draws on a web basemap.
       `LayerBuildResult.thickness_map` uses plain geopandas plotting with no
       basemap. For a real georeferenced stack the Choro is better; for
       `layer_management_workflow.ipynb`, whose synthetic (0,0)-(1000,600)
       EPSG:2927 domain reprojects to open ocean off Oregon (measured:
       -126.86, 45.15), it is worse. `.plot_mpl()` is basemap-free, so
       `stack.plot.map().plot_mpl()` reproduces today's picture exactly while
       `stack.plot.map()` gives the georeferenced one. Both spellings kept.

132. **8.5a: `stack.plot`, and a Picture that is not Plotly (2026-08-20).**
     Deferred out of 8.4 because binding these verbs is a RETURN-TYPE change, not
     a rename. Three views, three different violations: `thickness_map`/`preview`
     returned a matplotlib Axes, `cross_section` returned an Axes, and
     `surface_3d` returned a **bare `go.Figure`** (not `viz.Fig`, so no house
     template) with `html_path=`/`browser=` baked in.
     - **`viz.MplPicture` is the new idea.** Some drawings genuinely ARE
       matplotlib: the filled, layer-coloured section comes from FloPy's
       `PlotCrossSection` and has no Plotly equivalent to defer to. Rather than
       exempt it from the grammar or fake a `fig`, `MplPicture` answers the same
       four verbs over an Axes. `.fig` deliberately RAISES with a message naming
       `.axes`, because `fig` is documented package-wide as the Plotly figure and
       returning an `mpl.Figure` would break every caller expecting `.add_trace`.
       `.html(path)` writes a self-contained page with the PNG inlined — no
       network, no plotly.js. **No change to `Picture` was needed**: overriding
       all four methods was always legal; only its docstring claimed otherwise.
     - **COMPROMISE — `plot.map()` is basemap-FREE by default**, unlike every
       other `map` in the vocabulary. A `Choro` draws on a web basemap, and a
       stack under construction is usually on synthetic or local coordinates:
       `layer_management_workflow.ipynb`'s (0,0)-(1000,600) EPSG:2927 domain
       reprojects to open ocean off Oregon (measured -126.86, 45.15), so the
       shared choropleth would be strictly worse there. `basemap=True` opts into
       it for a genuinely georeferenced stack. Both spellings tested.
     - **`views()` deleted, not replaced.** It rendered four fixed pictures and
       returned None — Layer-3 display fused into Layer 1. The notebook cell now
       calls the four verbs, which is longer and honest; a general composer
       (`plot.mosaic`) takes plotly panels and cannot hold the matplotlib two.
       Recorded as a real, if small, ergonomic loss.
     - **Laziness moved an error.** `plot.surface("nope")` no longer raises at
       call time; the `KeyError` now surfaces when `.fig` assembles. That is the
       price of a Picture that builds on demand, and the test was updated to
       assert it at the new point rather than papered over.
     - **`layers.py` moved L3 -> L4** by importing `grid/plotting.py` for the
       `basemap=True` path. The first attempt reached `myflopy.plot` (L5) with a
       function-local import and tripped the ratchet 2 -> 3; the fix is the same
       one 8.4b used for `GridPlots` — call the L3 factory, not the L5 front door.
       `deferred_total` stays 60.
     - **The notebook was executed, not just edited.** `scripts/render_notebooks.py`
       needs nbconvert, which is not installed, so the code cells were run
       directly: all clean except a pre-existing cell hardcoding a **Windows**
       GRASS path (`AppData/Local/.../grass84.bat`), which cannot run on Linux
       and is unrelated to this stage.
     - Still on 8.5b: `vtk_3d` keeps its `IFrame` return and its `Path.cwd()`
       write, and becomes `plot.grid(backend="vtk")`.

133. **8.5b: the 3-D scenes get a contract, and `backend=` earns its name (2026-08-20).**
     `plot.grid(backend="vtk")` is the 3-D layered mesh; `pathlines=` draws
     particle tracks as tubes over it. `LayerBuildResult.vtk_3d`,
     `ParticleTrackingScene`, `export_particle_tracking_html` and
     `ModelVisualization`'s two particle wrappers are gone.
     - **The plan said `surface(backend="vtk")`; that would have been dishonest.**
       Measured, the three things draw three different SHAPES:
       `InterpolatedSurface` emits one `go.Surface(x, y, z)` -- a height field;
       `vtk_3d` renders a DISV cell VOLUME coloured by layer; the particle scene
       renders POLYLINE tubes. A `backend=` switch that turns a height field into
       a volume changes the SUBJECT, which `myflopy/plot/__init__.py` explicitly
       forbids ("Geometry chooses the verb. Not content, and not renderer.").
       `grid` already means "the mesh itself", so a 3-D layered mesh is that same
       subject redrawn -- and `backend=` is then a true renderer switch, matching
       the `map(backend="mpl")` precedent. Pathlines stay an OPTION, exactly as
       they already are on the 2-D map.
       **My own first framing of this was wrong** ("field vs geometry"):
       `surface_3d` already draws layer GEOMETRY via
       `InterpolatedSurface(surf_type="lyr")`. The real axis is height-field vs
       volume vs polyline.
     - **`viz.VtkScene` is the third renderer.** `Picture` is Plotly-shaped
       (`html` -> `fig.write_html`, `save` -> `fig.write_image`), and a
       `pv.Plotter` answers none of it -- measured: no `write_html`, no
       `write_image`, no `_repr_html_`. It has `export_html` and `screenshot`,
       so `VtkScene` overrides all four verbs, exactly as `MplPicture` (8.5a)
       does for Matplotlib. `.fig` raises and names `.scene`.
     - **Nothing is written to disk just by building a picture any more.**
       `vtk_3d` exported an HTML file on EVERY call -- into `Path.cwd()` when
       given no path -- and returned an `IFrame` pointing at it. That litter was
       real enough that **`.gitignore` carried a line naming the file**
       (`layer_vtk_3d_all.html`); the line is deleted with the behaviour, and a
       test asserts the working directory stays clean. Jupyter display now
       exports to a string, not a file.
     - **The optional dependencies were never actually optional.**
       `_optional._EXTRA_FOR_MODULE` had never heard of `pyvista` or `trame`
       despite the `viz3d` extra existing, and `layers.py` imported pyvista with
       NO guard at all. Both now route through `require(...)`, so a missing
       dependency names the extra. `dash` was added at the same time (the
       `Choro.dash_selector` path had the same gap).
     - **COMPROMISE — `VtkScene.html` inlines ~1 MB.** Measured 1,080,568 bytes
       for a trivial scene, self-contained with no remote `<script src>`. That is
       inherent to vtk.js offline export, not something this stage introduced;
       recorded because `_repr_mimebundle_` now embeds that per inline display.
     - **`PRTRunResults.export_3d_html` KEPT** rather than folded away. It is
       `self.scene(**kwargs).html(path)` with the plotter closed afterwards --
       "give me a file I can send" is a distinct intent from "show me", and the
       close matters for a long notebook session.
     - Ratchet moved DOWN, 60 -> 59: `prt.py` lost a deferred import when
       `export_3d_html` stopped reaching for the module-level exporter.

134. **8.2's `xs` -> `section` rename left a canonical notebook broken (2026-08-20).**
     Found while scoping 8.6. `canonical_fast_tour.ipynb` cell[9] called
     `model.hds.animate(kind='xs', line=line)`, and `kind='xs'` has raised
     `ValueError: kind must be 'map', 'plot', or 'section'` since 8.2 --
     `package_plotting.py:1614`. `docs/model_diff_cheatsheet.md:139` carried the
     same stale spelling. Both fixed.
     - **The tests were NOT the gap.** `test_view_grammar_composers.py:354`
       exercises `animate(kind="section")` directly and has passed throughout;
       the code was covered. What was missed is that **the suite does not execute
       notebooks**, so a notebook can name an API that no longer exists and stay
       green forever. `git log` confirms 8.2 (`35ddc73`) never touched
       `canonical_fast_tour.ipynb` despite renaming the verb it calls.
     - The Phase 8 notebook policy ("every stage updates the notebooks it
       breaks, in the same commit") was written precisely for this and was not
       followed in 8.2 -- the stage's own table listed `canonical_fast_tour`
       under 8.6 for the `animate` SIGNATURE change, which is a different break,
       and the `kind=` argument slipped between the two rows.
     - Corrective habit, applied from 8.5a onward and worth keeping: **execute
       the notebooks a stage touches**, do not just edit them. That is what
       caught this class of error in 8.5a and 8.5b.

135. **The Picture idempotence contract was proven against a STUB, and four real
     classes had drifted off it (2026-08-20).**
     8.1 wrote `test_assembling_twice_does_not_draw_twice` against
     `_TwoTraceStub` -- a class written in the test file to satisfy the contract.
     It proved the rule about itself and checked no real subclass. Applying the
     same assertions to the live classes found **four violators immediately**:
     - `XSection.fig` rebuilt from scratch on every access, so
       `xs.fig.update_layout(...)` then `xs.show()` silently dropped the edit.
     - `GridSection.fig` likewise.
     - `InterpolatedSurface.fig` delegated to `clipped_fig(clip=False)`, a fresh
       build each time.
     - **`GridMesh.fig` -- written by me in 8.4b**, with the same defect, three
       stages after the contract was introduced.
     Only `Choro` (fixed in 8.1) and `LayerSurface` (8.5a) were correct.
     - **`XSection.ani` was the `Choro.ani` bug again** (ledger 130): it built the
       animation figure, returned it, and never assigned `self._fig`, so
       `plot.section(...).ani` then `.show()`/`.html()` wrote the STATIC section.
       I fixed one sibling in `98814ec` and did not check the other.
     - The lesson is about the SHAPE of the test, not the bug: a contract test
       that constructs its own conformant subject is a tautology. The test now
       enumerates `Picture.__subclasses__()` and asserts every myflopy plotly
       picture is exercised, so a new one cannot silently opt out. Both guards
       are mutation-tested.
     - `MplPicture`/`VtkScene` are excluded by name -- their `.fig` raises on
       purpose -- and `LayerSurface` needs a built stack to instantiate.

136. **8.6a: `animate` becomes a combinator, and the plan's "only new code" stage
     turned out to be mostly promotion (2026-08-20).**
     `plot.animate(model, periods=)` returned `Animation` -- a config bag with
     `.sliders`/`.updatemenus` and nothing else. It was in `__all__`, documented
     as a verb, and could not `.show()`, `.save()` or `.html()`. It had **zero
     callers** outside its own docstring, so it is deleted rather than adapted.
     - **The Layer-2 combinator already existed, privately.**
       `SpatialView._plotly_animation(frames, title=)` took exactly
       `[(label, Picture), ...]`, and `_frame_panels` already produced that
       shape. The node grammar also already had a `backend="plotly"|"mpl"`
       switch. So the stage is: promote the private builder to
       `viz._build_frame_figure`, wrap it in a Picture, and add the one genuinely
       missing piece -- a raster backend.
     - **Delegating fixed a real bug.** `_plotly_animation` passed only
       `choro.get_choropleth()`, so contours, location markers and pathlines
       drawn over a map silently vanished between the static picture and its
       animation. The shared builder includes `overlay_traces()`.
     - **Two Picture classes, matching the renderer split viz.py already has.**
       `FrameAnimation` (plotly, `.fig` is the frames figure) and
       `SliderAnimation` (rasterized frames, `.fig` RAISES and names `.frames`,
       exactly as `MplPicture`/`VtkScene` do). Heterogeneous frames are a
       `SliderAnimation`, so they need no third concept.
     - **Mixed frames + `backend="plotly"` RAISE rather than falling back.**
       Silently returning a raster page instead of an interactive figure changes
       WHAT you get, which is the same class of lie that made 8.5 rename
       `surface(backend="vtk")` to `grid(backend="vtk")` (ledger 133). Plotly
       itself does not validate this -- it accepts mismatched frame traces and
       fails in the browser -- so the builder checks trace counts explicitly.
     - **The PNG bridge does NOT go through Plotly's static export.** Measured:
       `kaleido` v1 imports fine (so viz.py's friendly ImportError never fires)
       and then raises `RuntimeError: Kaleido requires Google Chrome`; kaleido is
       not in `pyproject.toml` at all. Every plotly picture here already carries
       a `plot_mpl()`, so `_mpl_figure_of` routes through those -- normalizing
       the three return shapes they use (Figure, Axes, `(fig, ax)`).
     - **CORRECTION to my own earlier claim.** I told the user the raster slider
       was ~0.7 MB against 17.7 MB for plotly -- a 25x win. At shipped defaults
       (dpi=140) it is 8.65 MB for 20 frames. Measured on the canonical model
       (441 cells x 6 frames): plotly 2.05 MB, png 1.10 MB. The honest statement
       is that **png size is independent of cell count while plotly grows with
       cells x frames**, so png wins on big grids, not universally.
     - **`export_matplotlib_slider_html` KEPT as the engine.** The user pushed
       back on deleting it, correctly: it takes a `render_frame` callable plus
       values and knows nothing about MODFLOW, so it is a general primitive, not
       a redundant wrapper. `SliderAnimation.export()` calls it and returns the
       full `StandaloneHtmlSlider` handle; `.html(path)` returns a `Path` like
       every other picture. Only the three MODEL-SPECIFIC wrappers are
       duplicative, and those go in 8.6b with `model.visualize`.
     - **COMPROMISE — `_build_frame_figure` lives in `viz.py`, not next to the
       slider machinery.** `package_plotting` is L1 and `interactive_plotting` is
       L2, so the node grammar could not import the builder from there without a
       deferred import, which the exact-match ratchet forbids. It is pure figure
       assembly and needs only what viz already has, so viz (L0) is its correct
       home anyway; `SliderAnimation` stays at L2 because it needs the exporter.
       `deferred_total` unchanged at 59.

137. **8.6b: `model.visualize` deleted -- but not for the reason the plan gave
     (2026-08-20).**
     The plan said to delete it "once its four slider exports are reachable as
     `<picture>.html(path)`". **They are not, and I verified it before acting.**
     `plot_model_head_map` renders through FloPy's `PlotMapView` --
     `plot_array` + `plot_grid` + `contour_array`, styled by a `ModelMapStyle`
     -- while `animate(backend="png")` rasterizes via `Choro.plot_mpl`, a
     geopandas cell fill. Different drawings, not two spellings of one.
     - **I told the user the three exporters were "redundant wrappers" and that
       was wrong.** They pushed back ("seems like the export_* functions are
       useful right?"), which is what prompted the check. Deleting them would
       have removed a picture. All four `export_*_slider_html` functions stay.
     - What IS redundant is the CLASS: five methods, each a one-line delegation
       to a module-level function already exported in `mf.__all__`. So the
       namespace goes and the functions stay --
       `model.visualize.head_map_slider_html(path)` becomes
       `mf.export_head_map_slider_html(model, path)`.
     - **COMPROMISE — ergonomics for consistency.** A method on the model reads
       better than a function you pass the model to. The trade is that Phase 8's
       rule is "pictures come from verbs on objects", and these are neither
       pictures nor verbs: they are exporters that write a file and return a
       handle. Keeping the namespace would have left a second way to make an
       animation sitting outside the grammar.
     - **The two plotly-animation methods folded into the verb.** They were
       `model.plot.map(...).ani` plus frame subselection plus an HTML write, all
       of which `plot.animate(frames).html(path)` now does -- with the frame
       selection explicit at the call site, which is honest about the cost (one
       rendered map per frame).
     - **`_write_plotly_choropleth_restyle_html` was re-homed, not deleted.** It
       had exactly one caller, `plotly_head_map_animation`, so deleting the
       method would have orphaned the real size fix: it ships the cell geometry
       ONCE and restyles values per frame. It is now what `FrameAnimation.html()`
       writes for a choropleth, with plotly's own writer as the fallback for
       anything else. Moved to `viz.py` (L0) with `_build_frame_figure`, since
       `package_plotting` (L1) cannot reach `interactive_plotting` (L2).
     - **A guard was moved DOWN rather than lost.** `zmin >= zmax` was validated
       inside the deleted method; it now lives in `_choropleth_factory`, so it
       covers every map instead of only the animated head map. An inverted range
       renders an all-one-colour picture with no error, which reads as a broken
       model.
     - Ratchet 59 -> 58 (`base.py` lost the deferred import behind the property);
       api_snapshot dropped `visualize` and `ModelVisualization`.

138. **8.7: the vocabulary pinned, and an allowlist that had gone stale in
     silence (2026-08-20).**
     `tests/test_plot_vocabulary.py` pins the verb set at all four scopes, that
     no scope invents a spelling outside it, and that every retired name stays
     gone. The conventions doc gained the three-layer frame.
     - **The plan's acceptance criterion was wrong twice over:** "all three
       scopes expose the same verb set". There are FOUR scopes, and they do not
       and should not match -- a bare grid has no results (no `surface`, no
       `animate`), a layer stack has no time (no `animate`). The test declares
       each subset and asserts `scope_verbs <= VOCABULARY`, which is the property
       worth having: a scope may answer fewer verbs, never different ones.
     - **`conftest._SLOW_TESTS` named two tests that 8.5b had renamed.** It marks
       heavy tests by STRING, so a rename silently un-marks one. Nothing failed
       -- the tests kept passing, just in the fast lane -- so two pyvista VTK
       exports had been running in the inner loop since 8.5b -- measured 3.1 s
       of pyvista work against a ~32 s fast loop, so roughly a tenth of it. A
       guard now asserts
       every name in the set resolves to a real test, plus a `_RETIRED_SLOW_TESTS`
       register asserted disjoint from the collected names. Mutation-tested.
     - **COMPROMISE — the conventions doc was EXTENDED, not rewritten.** The plan
       said "rewritten around the three layers". Its existing grammar,
       colorscale-policy, signed-exchange and naming sections were accurate and
       hard-won; a rewrite would have destroyed the specifics the file exists to
       record. The three-layer model is now the frame at the top and the grammar
       reads as one scope within it.
     - **COMPROMISE — `mosaic`/`animate` sit on `model.plot` but not on
       `vor.plot`/`stack.plot`.** They are subject-free combinators, so strictly
       they belong at module level only. They are on the model because that is
       the common entry point and discoverability was the phase's whole goal.
       Recorded as a deliberate asymmetry rather than left to look accidental,
       and pinned by the scope table so it cannot drift further.

139. **The master notebook's missing output file was a typo, not a regression
     (2026-08-20).**
     Cell 10 asserted ten MF6 output files exist and one did not:
     `{model.name}_lake_budget.csv`. The file MF6 writes is
     `{model.name}_lak_budget.csv` -- the package abbreviation, matching
     `_uzf_budget.csv` and `_sfr_budget.sfr` beside it. `git log -S` dates the
     typo to `fb24b54`, long before Phase 8; nothing else in the tree shares it.
     - Worth recording because I hunted it as a Phase 8 regression first. The
       notebook had presumably never been run end-to-end since that commit --
     which is the same gap ledger 134 records: the suite does not execute
       notebooks, so a broken cell stays green indefinitely.

140. **`mf.LayerStack` was exported without any of its arguments (2026-08-20).**
     Found writing the plotting tour: `mf.LayerStack` is on the facade but
     `Raster`, `Flat`, `Contours`, `Points`, `Array`, `Isopach`, `Min`, `Max`,
     `Clamp` and `Where` were not -- all ten are in `myflopy.layers.__all__` and
     none reached `mf`. So CLAUDE.md's own documented one-liner,
     `LayerStack(vor, top=Raster("ground.tif"))`, did not work from
     `import myflopy as mf`; it needed a second, deeper import for the argument
     types. All ten added to the lazy export map; api_snapshot diff is exactly
     those names.
     - The lesson is about what a facade owes: exporting a class without the
       types its constructor takes is only half an export. Writing a notebook
       that used the documented spelling is what surfaced it -- no test did,
       because the tests import from `myflopy.layers` directly.

141. **My notebook runner is not a Jupyter kernel (2026-08-20).**
     The master notebook "failed" at cell 25 with `NameError: display`. Jupyter
     injects `display` as a builtin; seven notebooks in this repo rely on that,
     which is normal practice. The failure was in my `exec`-based runner, not the
     notebook. Recorded so the next person checking a notebook this way injects
     `display` (and remembers that a bare `exec` loop differs from a kernel in
     other ways too -- no `_`/`__` history, no rich reprs, no cell ordering
     guarantees beyond what the loop imposes).

142. **The sectioned hover reached the grammar but not the map verb (2026-08-20).**
     User-reported: `model.plot.map()` showed the flat `Cell No. / Area / x / y`
     dump while `model.hds.map()` showed the sectioned hover the 6.x work built.
     Measured side by side -- the first resolved `hover_spec` to `None`, the
     second to a `HoverSpec`.
     - **NOT a Phase 8 regression.** `model.cor()` never passed a `hover_spec`
       either, so the flat hover is what the direct map verb has always given.
       What changed is that `model.plot.map` is now THE map verb, so "the
       documented way to draw a map" and "the way that gets the good hover"
       stopped being different things by accident and started being different
       things visibly.
     - Fixed where `show_layer_elevs` was fixed (ledger 128): the factory
       resolves a default from the map's `type` -- `hds`/`conc`/`temp` to
       `head_hover`/`conc_hover`/`temp_hover`, the same builders the grammar
       uses. `type="custom"` (the package/group path, which supplies
       `custom_hover`) and `rch`/`ks` are left alone, so their behaviour is
       unchanged.
     - **A bare grid still gets the flat hover, deliberately.** The sectioned
       one needs model context -- layers, periods, dates -- so the default is
       gated on `model is not None`. Both halves are pinned.
     - Lesson worth keeping: the two entry points to one picture had drifted
       apart with nothing asserting they agree. The test now compares the
       rendered `hovertemplate` of `model.plot.map` against `model.hds.map`
       directly, which is the only form that could have caught this.

143. **Docstrings are the API surface for a `**kwargs` facade (2026-08-20).**
     The plotting verbs forward to picture classes with 15-39 constructor
     parameters, so there is no signature to read -- the docstring is the only
     place a caller can learn what is accepted. They were prose-only: no
     `Parameters`, no `Returns`, no examples, which is what an editor renders as
     a structured tooltip.
     - All six verbs rewritten in NumPy style with every forwarded parameter
       documented, including the traps: pass a diverging colorscale as STOPS not
       a name (the plotly-to-matplotlib table reverses `rdbu`, ledger 69/70), and
       `zmin >= zmax` raises rather than rendering one flat colour.
     - **The BOUND methods are what an editor shows** on `model.plot.map(`, and
       they were 57 characters. Duplicating thirteen docstrings would drift
       within a release, so `_inherit_verb_docs` appends the free function's full
       reference to each at import, from the single source. Runtime
       introspection (`help()`, Jupyter `?`, most editor hovers) reads `__doc__`
       and gets it; a purely static reader still sees the short source
       docstring, which is written to stand alone.
     - **COMPROMISE — `StackPlots` is bound from `myflopy.plot`, not `layers`.**
       `layers` is L4 and `plot` is L5, so the binding cannot live next to the
       class; it happens where both are in scope. Recorded because it looks
       misplaced until you check the layering.
     - Writing the guard immediately caught `viz.mosaic` with no `Returns`
       section. Fixed rather than excluded from the check.

144. **Hover audit: every model-backed map is sectioned (2026-08-20).**
     Follow-up to ledger 142, asked directly: is it ALL maps now, or is more
     work needed? Built one map of each kind on the canonical model and read the
     resolved spec off each -- the verb, the dependent-variable readers, a
     static-array package (`npf.k`, `sto.ss`), a per-period array package
     (`rch.inputs`) and four cell-stress packages. All sectioned.
     - **`vor.plot.map()` is the one deliberate exception** and is not a gap: a
       bare grid has no layers, periods or dates, so there is nothing to section
       BY. It keeps the flat `Cell No. / Area / x / y` hover, which is the right
       answer for geometry.
     - Pinned by `test_every_model_backed_map_gets_the_sectioned_hover`, which
       builds ten real maps rather than checking one and generalizing. The
       paths differ -- some pass `hover_spec`, some `custom_hover`, some
       neither -- so breadth is the only honest check.
     - Noted in passing, NOT fixed: `model.packages.npf` exposes its nouns
       directly (`packages.npf.k`) with no `.inputs` tier, unlike the
       cell-stress packages (`packages.ghb.inputs.head`). The documented grammar
       in CLAUDE.md is `<pkg>.<inputs|results>.<noun>.<verb>()`, so the
       static-array packages deviate from it. Pre-existing and out of scope
       here; recorded so it is not re-discovered as new.

145. **Explicit signatures on the plotting verbs; a `**kwargs` tail kept (2026-08-24).**
     The 8.7 docstring pass documented ~30 parameters for `plot.map` and the
     user reported the IDE still showed `**kwargs`. Correct: PyCharm and Pylance
     are static analyzers, so `__doc__`, `__signature__` and `functools.wraps`
     all miss them. Only the `def` line reaches an editor. The verbs now name
     their parameters; the bound namespace methods repeat them.
     - **The duplication is the compromise.** `plot.map` and `ModelPlots.map`
       carry the same 35 names. Generating one from the other would put them
       back out of static reach, and a `.pyi` stub is ruled out by CLAUDE.md
       §4.7.7. `test_the_bound_model_verb_mirrors_the_free_one` fails naming the
       parameter that drifted, which is what makes the copies safe rather than
       merely duplicated.
     - **A `**kwargs` tail survives on every verb, deliberately.** Explicit
       parameters are ADDITIVE: naming the common ones without keeping the tail
       would break any call passing something unlisted. `map`'s tail is also
       genuinely open — it reaches the `go.Choroplethmap` trace, whose names
       Plotly owns and validates at render time.
     - **Restated defaults are a real hazard**, so they are machine-checked.
       `test_a_named_default_matches_the_link_that_owns_it` reads the expected
       value off the forwarding target at runtime; nothing is hand-written. It
       caught `show_layer_elevs` immediately: `_choropleth_factory` owns it with
       `None` (resolve from the grid) while `Choro` declares `True`, and the
       nearest link is the one a caller sees.
     - **`FieldMappable`'s `field=` sugar stays `(*args, **kwargs)` on purpose.**
       Measured: it dispatches to accessors with incompatible signatures —
       `LakConnectionsExplorer.map` has no `per`, `DrnInput.map` takes `per`
       positionally, `PRTCaptureView.map` has neither `per` nor `model`. One
       merged signature could only be written by lying about some of them. Its
       docstring now says so and points at the leaf. NOT a deferral: unifying
       those accessors would flatten distinctions that carry meaning.
     - Scoping correction worth keeping: the leaf verbs were assumed bare and
       measured otherwise. `packages.npf.k.map`, `ghb.results.q.map`,
       `model.hds.map` and the rest already resolve to concrete Explorer classes
       with explicit typed parameters. Only the four namespace-level dispatchers
       were bare.

     Defects the work surfaced, all fixed in the same pass:
     - `plot.grid` documented `layers`/`scale`/`color_by`/`cmap` under "Other
       Parameters". Those live on `LayerBuildResult._vtk_plotter`, reachable
       only via `stack.plot.grid()`; through `plot.grid` they raised. Meanwhile
       the seven parameters it *does* accept (`vertical_exaggeration`,
       `model_style`, ...) were undocumented, and `_vtk_plotter`'s
       `width`/`height` were undocumented anywhere.
     - `SimulationBase.plot` and `LayerStack.plot` property docstrings both
       omitted `grid` from the verb list — added in 8.5b, never documented.
     - `stack.plot.map(**kwargs)` silently DISCARDED its kwargs unless
       `basemap=True`: the Matplotlib branch ignored them. Now raises naming them.
     - `GridPlots.grid(**kwargs)` was a dead passthrough — `GridMesh` takes only
       the grid, so every kwarg raised. Signature is now `grid(self)`.
     - `_inherit_verb_docs` appended the free docstring behind a 70-dash rule,
       which is malformed NumPy (a section underline with no title above it) and
       made PyCharm drop structured rendering for the whole docstring. It now
       splices sections into one well-formed docstring and drops the `source`
       entry, which a bound verb does not take.

146. **GRASS auto-discovery reaches POSIX, and finds its own bindings (2026-08-26).**
     `_find_grass_launcher` globbed `grass*.bat` only — an OSGeo4W shape — so
     `/usr/bin/grass` was invisible and `mf.Contours` raised "Could not find a
     GRASS launcher … install GRASS via OSGeo4W" on a machine where GRASS worked
     fine. `shutil.which` over `_GRASS_EXECUTABLES` now runs after the `.bat`
     glob, and `_grass_modules` asks the resolved launcher for
     `--config python_path` and appends it to `sys.path`. Measured on this
     machine: launcher `/usr/bin/grass`, bindings `/usr/lib/grass84/etc/python`,
     both found with `GRASS_BIN` **and** `PYTHONPATH` unset — the two env vars
     the cheat sheet used to require are now genuinely optional.
     - **The `which` fallback is gated to non-Windows, deliberately.** A GRASS on
       Windows `PATH` but outside the OSGeo4W/QGIS bundle directories is still
       not auto-discovered. The gate keeps Windows resolution byte-identical
       (the platform the `.bat` search was written for and the one that can't be
       re-measured here), and `GRASS_BIN` covers the gap. Revisit if a Windows
       user reports a PATH-only install.
     - **`_GRASS_EXECUTABLES` is a fixed tuple**, so a future `grass9` binary
       installed *without* the unversioned `grass` symlink would need `GRASS_BIN`.
       Nearly every distribution ships the plain name, so the list is a fallback
       for side-by-side installs rather than the primary path. Revisit: add names
       when GRASS 9 ships.
     - **The launcher is re-queried per `run()`**, not cached. It is one
       subprocess against an interpolation that takes seconds to minutes, and a
       module-level cache would go stale across a `GRASS_BIN` change inside one
       session. Not worth the state.
     - **The retry path cannot be unit-tested without GRASS installed**, so it
       carries `# pragma: no cover`. What *is* tested is everything around it:
       PATH discovery, env-var precedence, both `--config` shapes, the
       launcher-won't-run and no-launcher-at-all branches, the `sys.path`
       de-duplication, and — in both directions — that the error message says
       OSGeo4W on Windows and does not on POSIX. The end-to-end import was
       verified by hand instead, recorded above.
     - **The notebook preflight reaches into private functions.**
       `layer_management_workflow.ipynb` cell 25 hardcoded an
       `AppData/…/grass84.bat` path (flagged in entry 132 as un-runnable on
       Linux); it now calls `_default_grass_bin` / `_grass_modules` to report
       what discovery found. There is no public API for "is GRASS available",
       and inventing one for a preflight cell was out of scope — the cell says
       so in a comment.
147. **The classmethod-binding tripwire is arity-based, so a 2+-argument
     misbinding still slips through (2026-08-26).**
     `Surface.maximum/minimum/clamp/where/isopach/shift` are classmethods, so an
     instance call binds the instance to `cls` and drops it. Measured before the
     fix: `a.maximum(b).operands` was `(b,)` — no error, no warning, a
     plausible-looking Surface that ignores one input, and downstream a wrong
     contact elevation that surfaces only as odd heads. `maximum`/`minimum` now
     reject fewer than two surfaces and `clamp` requires `lower=` or `upper=`,
     each naming both fixes (`Surface.maximum(a, b)` and the fluent
     `a.floored_at(b)`). `shift`/`isopach`/`where` already raised for the
     now-missing positional argument and are pinned by test rather than guarded.
     - **The compromise: `a.maximum(b, c)` still drops `a` silently.** Two
       operands arrive, the arity guard is satisfied, and nothing distinguishes
       it from a legitimate `Surface.maximum(b, c)`. Judged acceptable because
       the fluent mental model that produces the bug is binary — someone who
       writes three surfaces is already thinking variadically and spells it on
       the class. Same residue on `a.clamp(b, upper=c)`. Revisit if either
       spelling turns up in a real notebook.
     - **A class-only descriptor was considered and rejected.** Replacing
       `@classmethod` with a descriptor that refuses instance binding would close
       the hole exactly, at every arity. It was not taken because it puts the
       *correct* spelling out of static reach: PyCharm and Pylance read the `def`
       line and the decorator, so a custom descriptor loses the "this is a
       classmethod" signal and `Surface.maximum(a, b)` starts drawing a
       wrong-argument squiggle. Making the right call look broken in the editor
       is a worse trade than the residual hole (see ledger 145 — only the `def`
       line reaches an IDE). Raising on attribute access would also break
       `hasattr`/`inspect.getmembers` over an instance.
     - **`where` gets no guard even in principle.** `Surface.where(zone, inside,
       outside)` takes three required positionals: `a.where(zone, b)` raises for
       the missing one, and `a.where(zone, b, c)` is semantically the call the
       user wanted — `a` plays no role in the correct spelling, so nothing is
       lost. The test pins the loud failure so a later default value cannot
       quietly reopen it.
     - **One-surface envelopes are now an error even when deliberate.**
       `Surface.minimum(x)` returned `x`, a no-op; making the arity the tripwire
       costs that. No call site in the repo splats a variable-length list into
       either constructor (checked across `src`, `tests`, `docs`, `examples`), and
       a splat that collapses to one element is itself worth hearing about.

148. **`LayerStack` accepts its grid late; `vor` stays first-positional (2026-08-27).**
     `LayerSurfaces` could always be declared before the grid; `LayerStack` could
     not, because `vor` was a required constructor argument. That forced a choice
     between the project-first ORDERING (declare the layering, then build the
     grid) and the per-layer CONTROL only the facade offers (`thickness=`,
     per-layer `min_thickness`/`pinch`). `vor` is now optional and
     `build`/`qc`/`to_disv` accept one, mirroring the override `to_disv` already
     had. `for_grid(vor)` returns a BOUND COPY.
     - **`vor` stays the first positional parameter** rather than moving to
       keyword-only, so `LayerStack(vor, top)` — which is in every notebook, doc
       and example — keeps working untouched. The cost is that `LayerStack(ground)`
       reads like the deferred form and silently binds a Surface to `vor`; that is
       guarded with a TypeError naming both spellings, so it is a stop rather than
       a failure three frames down in `sample()`.
     - **`for_grid` copies rather than mutates.** The layering is the expensive
       thing to write and the grid is the thing you vary, so one declaration can
       serve a coarse test grid and a fine production grid. Mutating would make
       the second binding silently reuse the first.
     - **`.plot` cannot take a grid** — it is a property — so on a deferred stack
       it raises naming `build(vor).plot` and `for_grid(vor).plot`. Considered
       making it a method for symmetry; rejected, because `stack.plot.map()` is
       documented in CLAUDE.md and four notebooks, and the whole plotting
       vocabulary treats `.plot` as a namespace rather than a call.
     - NOT done: `LayerSurfaces` still takes only a global `minimum_thickness` /
       `pinch` in `to_disv`. Per-layer rules there would duplicate what the facade
       exists to provide. The cheat-sheet table records this as the real remaining
       difference between the two.

149. **`region_vector` promoted from `**grass_kwargs` to a real parameter (2026-08-27).**
     Reported as "region_vector doesn't appear to be an arg of mf.Contours". It
     worked, but only as a `**grass_kwargs` passenger, so it never appeared in the
     signature -- the same static-visibility class as ledger 145 and the lazy
     exports, on a parameter that is REQUIRED half the time.
     - **The real defect was the cache, not the naming.** `sources` is
       content-tracked and `params` holds values only. `region_raster` was a field
       and went into `sources`; `region_vector` rode kwargs, so only its PATH
       STRING was recorded. Editing the domain polygon left `cache_status()` at
       `fresh` and returned a raster interpolated over the old extent. Both region
       kinds are now in `sources`.
     - **"Exactly one region" is now checked at construction**, not in
       `_set_region`, which runs only after a GRASS session has started -- late,
       and buried in GRASS's own output.
     - **Behaviour change, accepted:** constructing a region-less `Contours` used
       to succeed and fail later at interpolation; it now raises immediately. Two
       test fixtures relied on the old laxity (they check `cache_status` without
       ever interpolating) and were updated. Existing caches invalidate once,
       because `params` changes shape when `region_vector` leaves `grass_kwargs`.
     - NOT done: myflopy still does not reproject contours to the region's CRS.
       GRASS is invoked with `-o`, so a mismatch is accepted silently rather than
       corrected -- measured cost, a 597,000 ft offset that hangs `r.surf.contour`
       at 0% forever. A CRS check and a no-overlap check are filed as a task; the
       cheat sheet documents the manual `to_crs` in the meantime.

150. **`clip` promoted alongside `region_vector`; it stays in the cache signature (2026-08-27).**
     Same invisibility as ledger 149 -- reachable only through `**grass_kwargs` --
     on a flag that defaults to ON and shapes the output (with a vector region it
     clips to the POLYGON, not its bounding box).
     - **The care needed was the opposite of `progress`'s.** `progress` is
       deliberately EXCLUDED from the derived-raster signature: it changes what
       you see, not what is written. `clip` changes what is written, so moving it
       out of `grass_kwargs` would have silently dropped it from the signature
       and let a clipped and an unclipped surface share one cached file. It is
       added to `params` explicitly, with a test that fails if it leaves.

151. **`mf.Raster(..., nodata=)` overrides the header (2026-08-27).**
     A raster clipped to a boundary but exported with no no-data value in its
     header is a common GIS product; the fill (usually 0) is then read as real
     ground. Measured: 4.3M pixels of exact 0.0 outside a domain, collapsing all
     three layers to `min_sep` once contacts were capped to it.
     - **Masking is threaded into both samplers rather than applied after.**
       `method="area"` averages each cell's pixels, so masking afterwards leaves
       a boundary-straddling cell already blended -- half 0, half real ground --
       reading as a plausible elevation with nothing left to detect. Pinned by a
       test asserting a straddling cell returns 500, not 250.
     - Writing it produced a bug now covered by its own test: a cell whose every
       pixel is masked counts as UNCOVERED and falls to the centroid fallback,
       which did not forward the override and handed the sentinel back.
     - NOT done: no heuristic warning for "this raster declares no nodata and is
       suspiciously full of zeros". Too clever -- 0 is a legitimate elevation.
       The docs name the durable fix (`gdalwarp -dstnodata`) instead.

152. **Interpolated surfaces drop no-value samples (2026-08-27).**
     `griddata` triangulates every point it is given, so ONE NaN vertex makes
     every output cell of every triangle touching it NaN; scattered NaNs erase
     the whole surface. Measured: 300 valid cells of 582 gave 0 finite pixels of
     14,400, which renders as nothing and reads as "plotting is broken".
     - Pre-existing, not introduced by 151 -- any properly declared no-data
       raster reached it. `nodata=` merely made it easy to reach.
     - Fixed at `_finite_samples`, which griddata AND rbf both draw from, rather
       than at either call site.
     - An all-NaN surface now RAISES naming the likely causes. Returning an empty
       picture is what made the original failure so hard to place.

153. **`surface(backend="vtk")` renders contacts as separate sheets (2026-08-27).**
     `grid(backend="vtk")` fuses layers into one volume, so no contact can be
     isolated in it -- which is the actual question when someone wants to see
     surfaces individually. Separate actors let a viewer hide the sheets above.
     - **Does not weaken the `backend=` rule**: `surface` still means a height
       field and only the renderer changes. Both paths interpolate through the
       shared `InterpolatedSurface`, so a VTK sheet and its Plotly counterpart
       are the same numbers.
     - NaN cells are thresholded out. Drawn, they would be a sheet at z=0 --
       a contact at sea level, exactly the failure 151 exists to prevent.
     - `show_edges` is recorded in `test_plot_vocabulary`'s stack-local allowlist
       rather than added to the free verb: it is a vtk sheet detail, like
       `scale`/`cmap`/`width` beside it. The guard rejected it first; the
       allowlist entry is the reviewed answer, not a workaround.
     - The import ratchet rejected the new deferred import; hoisting it removed
       the duplicate in `_surface_fig` too, so the deferred total went 58 -> 57.
     - **DEBT REPAID 2026-08-28 (see 157).** "Separate actors let a viewer hide
       the sheets above" justified this entry and was not actually reachable:
       `add_mesh` was called with `label=` but not `name=`, so `plotter.actors`
       was keyed by address strings like `UnstructuredGrid(Addr=0x1f3bd750)` --
       and on the volume path that dict also holds the scalar bar, so indexing
       by position hit it. The sheets are now named for their contacts.

154. **A bare EPSG number is normalized to `EPSG:<n>` (2026-08-27).**
     `crs="2927"` failed with `CRSError: The WKT could not be parsed`, nine frames
     deep in raster sampling and several steps from the line that was wrong.
     - **Why it was so far from its cause:** pyproj accepts a bare `"2927"`, so
       geopandas reprojection, the boundary clip and the whole grid build
       succeed. Rasterio's `warp_transform` is the only thing that rejects it, so
       the first symptom is the first RASTER sample.
     - Normalized at `VoronoiGridPlus.__init__`. A bare integer string is never
       valid WKT or PROJ, so reading it as EPSG cannot be wrong. It also lines up
       with `mf.Contours(epsg="2927")`, which has always taken the bare form --
       the inconsistency that made `"2927"` natural to write.
     - NOT done: no general CRS validation at construction. Only the unambiguous
       case is handled; anything else still passes through to pyproj/rasterio to
       interpret and to complain about in their own words.

155. **`VtkScene.html()` rewrites PyVista's module script so the page opens from disk (2026-08-27).**
     Reported as "grid.html doesn't display correctly in Firefox" -- it showed
     vtk.js's "Drop File / Explore Scene" placeholder.
     - **Not a Firefox bug and not a rendering bug.** PyVista emits its bundle as
       `<script type="module">`; a module script is not executed from `file://`.
       The bundle's last act is `window.OfflineLocalView = {...}` and a following
       CLASSIC script calls `OfflineLocalView.load(...)`, so the shell runs, the
       loader does not, and the console says
       `ReferenceError: OfflineLocalView is not defined`. Served over HTTP the
       identical file is fine -- which is exactly why it looked browser-specific,
       and why my first investigation (served over localhost) found nothing.
     - **Post-processing someone else's bundle is the compromise.** It is one
       string replacement of a marker myflopy does not own, so
       `test_vtk_html_export.py` asserts the assumption it rests on: no
       `import`/`export`/`import.meta`/dynamic import and no top-level `await`
       anywhere in the emitted scripts. If a future PyVista ships module syntax,
       that test fails rather than writing a page that silently throws.
     - Ordering improves rather than degrades: a classic script runs at parse
       time, BEFORE the consumer below it, where the module was merely deferred
       and happened to win a `setTimeout(..., 0)` race.
     - Applied to `_repr_mimebundle_` too, since notebook output has the same
       constraint, and pinned as still self-contained (no network fetches).

156. **A section line may be what its docstring said, and `color_by` is checked (2026-08-28).**
     Both `StackPlots.section` and the free `plot.section` documented a
     `LineString` or a `Path`, and `LayerBuildResult._resolve_line` did a bare
     `[tuple(pt) for pt in line]` -- so the only form that worked was the one
     form neither docstring mentioned.
     - **It failed LAZILY**, which was worse than failing. A `LayerSection` is a
       Picture, so `section(line=LineString(...))` returned fine and blew up
       later on `.axes` with `TypeError: 'LineString' object is not iterable`,
       pointing into FloPy rather than at the argument. The line is now coerced
       in `LayerSection.__init__` -- at the call -- not at draw time.
     - Reuses `_as_linestring` and `read_shp_gpkg`, the same pair
       `GridPlots.section` already reads the same inputs with, rather than a
       second reader that could drift from it.
     - **No broad `except` around the reader**, per the standing rule: a broad
       handler there would eat the point/length validation below it, which is
       the exact bug already found twice in `read_gpkg` and
       `contour_line_segments`.
     - `color_by` is validated on `section` (`layer`/`thickness`) and on
       `grid`. Both silently drew the wrong picture for an unknown value --
       `section(color_by='banana')` fell through to the layer-coloured branch,
       and on `grid` anything that was not `'layer'` MEANT `'elevation'`.

157. **The 3-D volume carries its numbers, and comes apart by layer (2026-08-28).**
     Asked for layer toggles, cell readouts and quick sections in the VTK scene.
     Three of those four wants were already answered elsewhere (see the notes at
     the end); what was genuinely missing was that the mesh carried no data and
     the scene had no handles.
     - **Float arrays must NOT go through flopy's `Vtk.add_array`.** It masks
       float arrays to NaN wherever `idomain == 0` -- i.e. on exactly the cells
       a `pinch="inactive"` layer creates, the ones worth inspecting. Measured:
       138 of 843 cells on a stack with one pinching layer; integer arrays come
       through intact. `thickness`/`top`/`botm`/`cellid` are attached to the
       PyVista mesh after `to_pyvista()` instead.
     - **`top` and `botm` are overwritten, not added.** FloPy's own `top` cell
       array is NaN for EVERY layer below 0, so the per-layer top is rebuilt as
       `vstack([top, botm[:-1]])`.
     - `color_by` widens to `layer`/`thickness`/`top`/`botm`/`cellid`/
       `elevation`. **`elevation` is kept although it overlaps `top`**: it is a
       POINT array of vertex z, so it ramps within a cell where `top` is a flat
       per-cell contact. Two names, genuinely two pictures -- but close enough
       that the docstring says which is which.
     - `cmap` defaults to `None` = "whatever suits this scalar" rather than
       `tab10`, which is a qualitative map and wrong for a continuous field.
       This makes `cmap` honest for every scalar for the first time; previously
       the elevation branch hardcoded `terrain` and ignored the argument.
       `_VTK_GRID_DEFAULTS` is mirrored, so it moves with the signature.
     - **One actor per layer instead of one fused mesh.** Costs `nlay` draw
       calls; buys the only thing that makes a toggle expressible. Actors on
       both paths are now `name=`d, so `scene.scene.actors["clay"].visibility
       = False` works.
     - **The scalar bar is not deduplicated by hand.** Every actor asks for one
       and PyVista keys bars by TITLE, so the scene gets exactly one. An
       explicit `clim` is passed anyway: PyVista *also* syncs each mapper to the
       shared bar's LUT, which would paper over a divergence, and relying on
       that is relying on undocumented behaviour. **The consequence is that the
       clim cannot be pinned by a test** -- `mapper.scalar_range` returns the
       same value with or without it. The test says so rather than asserting
       something that cannot fail.
     - `scene.meshes` is populated at last (it was documented on `VtkScene` and
       left `()` by both layer scenes). For the volume it is the ASSEMBLED mesh,
       not the per-layer pieces the actors draw -- a deliberate looseness,
       because `scene.meshes[0].save("stack.vtu")` for ParaView wants the whole
       selection with every array on it.
     - `_vtk_surface_plotter` omitted `off_screen`, so `surface(...).save(...)`
       raised "Nothing to screenshot" in a plain script while `grid(...)` wrote
       a PNG. **`tests/conftest.py` sets `PYVISTA_OFF_SCREEN`, so the suite
       structurally could not reproduce it** -- which is how it survived. The
       test neutralises the global and pins the constructor argument.
     - A layer may no longer be named `top` (it collided with the model-top
       contact and with the actor name, and used to fail much later inside
       pandas as "cannot reindex on an axis with duplicate labels"), and asking
       for the same contact twice draws it once.
     - **NOT done, and deliberately.** Interactive toggling/picking/slicing on
       the live scene, and visibility checkboxes injected into the exported
       page, were both scoped and deferred pending which viewing path matters.
       Measured constraints if they are ever picked up: the exported vtk.js page
       has NO widget manager, NO picker and NO keyboard, and serializes ONLY the
       active scalar (7 of 8 attached arrays appear zero times in its
       `index.json`); trame has no keyboard channel at all; `vtkButtonWidget`
       needs a hover event trame's client never sends; and `add_mesh_slice`
       bakes eight dead actors into any export. Cross-sections stay
       `stack.plot.section()`, which already answers the question.

158. **A unit splits into model layers; geology and discretization separate (2026-08-28).**
     Asked for "additional surface algebra" to split declared layers by amount or
     percentage. Half the request already worked and half was inexpressible.
     - **Fixed thickness always worked**; fractional splits did not, and could
       not. `B - A`, `A * (1/3)`, `A / 3` and `-A` all raise, deliberately
       (`test_subtracting_a_surface_from_a_surface_is_unsupported`), and the only
       subtraction in the algebra is `isopach` with one operand hardwired to
       `previous`. So the per-cell interval `top - bottom` had no spelling.
     - **ONE new primitive, not general arithmetic.** `Surface.toward(target, f)`
       = `previous + f*(target - previous)`. The recurrence form is what keeps it
       lazy: `values` carries a single `previous` slot, so a cut measured from
       the contact above needs no anchor plumbing and no grid. General
       `Surface - Surface` was reconsidered and still refused -- it would widen
       the elevation-vs-thickness ambiguity that this same pass had to fix
       elsewhere (see below).
     - **`split=` lives on `.add()`, and there is no `.split()` verb.** Not
       ergonomics: `split` and `pinch`/`min_thickness` are ONE decision. Judged
       per slice, a 2.5 ft unit split three ways returns `idomain [0, 1, 0]` --
       inactive cells inside a unit that is fully present, and 0 blocks vertical
       flow. Whoever expands the split must hold the unit's pinch policy, which
       is `.add()`. `replace(name, split=N)` covers editing afterwards; it
       already mutates by name, returns self, and keeps unmentioned fields.
     - **The facade renormalises, because the naive form is silently wrong.**
       Cumulative fractions 1/3, 2/3, 1 fed to `toward` give thicknesses
       [20, 26.67, 13.33], not thirds -- each cut resolves against the previous
       CUT. `g_i = share_i / (1 - sum(shares before i))`, with the last step
       forced to exactly 1.0 so the base lands ON the declared bottom rather than
       a float's width above it. That last line is the difference between
       conservative and nearly-conservative, and the test asserts it with
       `array_equal`, not `allclose`, or it would pass either way.
     - **`trigger_sep` is now exposed and scales as `1.0 / max split`.** It was
       accepted by `LayerSurfaces.sample` and never passed by `build`, so it was
       always 1.0 -- sized for units and destructive to their slices: a 2.4 ft
       unit split three ways came back `[0.1, 1.5, 0.1]` = 1.7 ft with the base
       moved 0.7 ft and the layer below silently absorbing it. Scaled, it is
       exactly `[0.8, 0.8, 0.8]`. An unsplit stack gets 1.0, so nothing already
       built changes shape.
     - **COMPROMISE, measured and not fixed:** a unit that pinches out ENTIRELY
       grows by `(N-1) * min_sep` when split -- 0.1 ft at N=1 up to 0.5 ft at
       N=5 -- because reconcile floors each sub-contact independently, and
       everything below shifts down with it. Fixing it needs either a per-layer
       `min_sep` (a signature change to `reconcile_surfaces`, which
       `vor.reconcile_surfaces` also exposes) or collapsing a pinched unit's
       sub-contacts before reconcile rather than after. Both are wider than this
       change; the cells involved are `idomain = 0` by construction, and a
       smaller `min_sep` shrinks the drift proportionally. **Promoted to its own
       deferral as 165**, with the measurements and both candidate fixes.
     - **FOLLOW-UP 2026-08-31: `names=` added, generated names kept as the
       default.** Asked for on the grounds that `bottom layer 1_1` reads badly.
       Three judgment calls, all in the direction of refusing rather than
       guessing:
       - **All or nothing.** A partial list (name slices 1 and 3, generate 2)
         was considered and refused: which model layer a given name refers to
         would then depend on where the caller stopped counting, and a layer
         name is a DataFrame column and a 3-D actor, not a comment.
       - **`names=` requires a split of 2+.** `names=["x"]` on an unsplit unit
         is a rename of the unit, and the unit already has a name; accepting it
         would give two spellings for one thing and quietly break the `units`
         key. `split=1` is a documented no-op, so it is refused there too.
       - **`replace` validates `split` and `names` as a pair.** Changing one
         without the other raises rather than dropping or recycling names the
         caller wrote — dropping a split from a named unit is
         `replace(name, split=None, names=None)`. Slightly more typing, and the
         alternative is exactly the silent renaming the rest of this entry
         exists to prevent.
       The UNIT key is untouched by naming, so `units` / `per_layer` /
       `plot.grid(layers=<unit>)` behave identically either way — asserted by
       building the same stack twice and comparing `botm`/`thickness` with
       `array_equal`. Not pinned in `api_snapshot.json`: it covers the `mf.*`
       helpers and export lists, and `LayerStack` methods were never in it.
     - **NOT a `split=` mode: fixed lifts.** "20 ft layers" on a unit of varying
       thickness needs a varying number of layers, and `nlay` is global. A
       `thickness=`-declared unit just divides (60 ft in 20 ft lifts is
       `split=3`); a surface-bounded one repeats a constant cut, which is exact
       today and is documented in the cheat sheet instead.
     - **Not splittable: an `Isopach` bottom.** It is measured from the layer
       above, so each cut would measure from the previous cut. Refused at build
       with a message naming both fixes. The guard recurses into operands,
       because a composite (`Flat(150).capped_at(Isopach(map))`) is a `min` node
       that looks absolute over an operand that is not.
     - **`LayerBuildResult.units` + `per_layer()` instead of full unit-first
       plumbing.** Nothing downstream of `layers.py` reads layer NAMES -- npf/
       sto/ic take positional lists, PEST takes `tuple[int,...]`, `.gpkg` reads a
       1-based integer column -- so a split cannot break a consumer, but it does
       silently RETARGET them. `per_layer()` closes that for the code-side cases
       at the cost of one field and one method. **Deferred:** unit names accepted
       directly by `mf.npf`/PEST `layers=`/`.gpkg` `layer_field`. The `.gpkg`
       integer attribute is the one case the library cannot protect, because the
       number lives in a GIS file outside the code.
     - **Fixed in the same pass: `thickness=<Surface>` was silently an
       ELEVATION.** `thickness=Flat(20)` under a top of 100 set the bottom to 20
       -- an 80-thick layer -- while the docstring promised a thickness; a raw
       ndarray at least failed loudly. Now refused, EXCEPT for the three kinds
       that really do measure from above (`isopach`, `constant_thickness`,
       `offset_below`), which mean there exactly what they say. The first attempt
       rejected all Surfaces and broke `test_facade_isopach_layer`, which was
       testing correct behaviour -- the narrower rule is the right one.
     - `to_disv` now takes its idomain from `build()` rather than letting the
       engine judge each layer alone; otherwise the two disagreed on any split
       stack, which is a worse bug than either answer.

159. **A MODFLOW-USG model imports, and says what it could not bring (2026-08-28).**
     `mf.read_usg(nam, gsf=)` reads a MODFLOW-USG model into a `UsgModel` and
     `.to_mf6()` converts it to a `SimulationSpec` on a DISV grid. The reference
     model is the Ten Trails valley model (47,025 nodes = 5 x 9,405, 72 stress
     periods). Everything below is deliberate; `report()` and `validate()` state
     each one at runtime so no omission is silent.
     - **CLN is NOT converted.** MODFLOW 6 has no Connected Linear Network. The
       804 CLN nodes are read, segmented by graph shape into 4 waterbodies (707
       nodes, 2-D meshes) and 2 streams (97 nodes, chains), and reported with the
       layer-1 cells each touches -- but nothing is written. Dropping them costs
       no pumping at all (the entire WEL package is CLN-local P-ET, `ITMP = 0`
       GWF wells in every period) but does remove the lake/stream stage feedback
       and ~15,000 ft3/d of net surface-water flux. `model.cln_polygons()` returns
       the features as polygons, which is what a later LAK/SFR rebuild starts from.
     - **ETS becomes a LIST-based EVT, and it is large.** `NETSEG = 2` and
       MODFLOW 6 cannot combine segments with `READASARRAYS`, so the array
       package becomes one record per active column per period: 72 x 9,090 =
       654,480 records, a 61 MB EVT file. Using `EVTA` instead would silently
       drop the segment shape, which is a change to the physics, not the format.
     - **`NETSOP = 3` is resolved once, and that is exact.** USG applies ET to
       the highest ACTIVE cell; MODFLOW 6 needs an explicit cell. IBOUND is
       static for the whole run, so resolving it at conversion is exact rather
       than an approximation.
     - **CHD loses its within-period ramp.** USG interpolates `shead` -> `ehead`
       across a stress period; MODFLOW 6 holds one value. The end-of-period head
       is taken (max |ehead - shead| here is 1.40 ft).
     - **SMS -> IMS is partial, by necessity, but keeps the tuning.** The
       delta-bar-delta under-relaxation and backtracking controls map field for
       field and ARE carried -- they are how the original model was made to
       converge, and discarding them leaves MODFLOW 6 taking Newton steps of
       100,000+ ft on a model whose heads span 400. What does NOT carry: SMS's
       `HICLOSE` is not MODFLOW 6's `INNER_DVCLOSE` (different inner solvers), and
       `IACL`/`NORDER`/`LEVEL`/`NORTH`/`RCLOSEPCGU` have no counterpart, so the
       complexity preset owns the inner solve outright. Mixing half of SMS's
       numbers into half a preset produces a solver that is neither.
     - **`COMPLEX` is the default preset, not `MODERATE`.** On this model
       `MODERATE` does not converge slowly -- MODFLOW 6 dies with SIGFPE during
       the first solve. Measured; `COMPLEX` runs the same model.
     - **The mesh is written to a LOCAL origin.** MODFLOW 6 builds DISV
       conductances from raw vertex coordinates, and on a State Plane grid
       (~1.34 million ft here) that arithmetic loses enough precision to return a
       NaN budget while still reporting "Normal termination". Measured on identical
       input differing only in the shift: **0 of 9,405 cells finite as-is,
       9,405 of 9,405 shifted**. `to_mf6(local_origin=True)` (the default) writes
       the mesh relative to its own corner and declares that corner as
       `xorigin`/`yorigin`, so the model stays georeferenced.
     - **`fix_for_mf6` is OFF by default.** MODFLOW-USG accepts a head boundary
       below its cell bottom and MODFLOW 6 refuses to run. The opt-in raises GHB
       heads to just above the bottom (a relative nudge -- exact equality still
       trips MODFLOW 6 once the written text rounds) and omits CHD records in the
       periods where they sit below, which is inert in USG anyway. Off by default
       because a boundary head is not something to change unasked.
     - **No DRAWDOWN output.** The USG OC asks for it; MODFLOW 6 does not produce
       drawdown.
     - **Boundary files are truncated to `NPER`.** They routinely carry more
       period blocks than the model runs (this WEL holds 612 and its CHD 792 for a
       72-period run) and MODFLOW reads the first `NPER`, so the reader does too.
     - **A start date is refused rather than invented.** This DISU parks a date on
       each stress-period line whose year is the literal constant 2023 on all 612
       lines, so the sequence jumps backwards every January. Dates are used only
       when strictly increasing; otherwise `to_mf6(start_date_time=)` must supply one.
     - **DISU is read through FloPy; BAS6 and LPF are not.** FloPy 3.10 cannot load
       a BAS whose `STRT` comes from a binary unit (`EXTERNAL -61`): it routes the
       record to `Util2d.load_txt`, which tests `"," in line` against `bytes`. LPF
       then fails too because it reaches for the BAS that never loaded. Both are
       read here instead; DISU, which FloPy handles correctly, is delegated.
     - **Non-layered USG grids are refused, not approximated.** A nested or
       ghost-node-refined DISU has no DISV equivalent, so `require_layered=True`
       raises rather than misplacing cells. The test is on the connectivity, not
       on `IVSD`, because `IVSD` states intent while the connections state fact.

160. **A highlight is a drawing, not a dimming; and it works on every map (2026-08-29).**
     `map(select=...)` moved from grid-only to every scope, and changed mechanism.
     Asked for: "add selections to any map... is there a better way to highlight
     certain cells?" Both halves were answered by measurement.
     - **`selectedpoints` is the wrong mechanism for a field.** Plotly's selection
       styling on a choroplethmap exposes exactly one property --
       `go.choroplethmap.selected.Marker()` is `['opacity']` -- so "highlight"
       there can only mean *erase everything else*; unselected cells render at
       `0.2 x opacity` (plotly's `DESELECTDIM`). Measured on a real head field,
       200,000 random pairs of UNSELECTED cells in CIE Lab: perceptually
       distinguishable pairs fall from **67.2% to 9.7%**, a 6.9x loss of contrast
       across the cells you did not select, while the colorbar still advertises
       the full range. Six of ten adjacent `earth` deciles become
       indistinguishable.
     - **The default is now a dissolved-boundary outline** drawn as one
       `go.Scattermap`. The field keeps full opacity, and 804 contiguous cells
       dissolve in 17 ms for +0.034 MB on a 5.3 MB page. A scattered 800-cell
       selection costs 30 ms and +0.197 MB, which is why the surveys' proposed
       per-cell-`marker.line` fallback was NOT built: one geometry route, no
       heuristic threshold.
     - **`select_style="dim"` reproduces the old picture verbatim** (the 0.2 is
       now the named constant `viz.HIGHLIGHT_DIM_OPACITY` rather than plotly's
       implicit default). `"both"` draws each. **The grid-scope default DID
       change** -- a deliberate call, approved: three call sites in the repo
       (`utils/inputs.py:156`, one test, the USG workbook) and all three read
       better as outlines.
     - **`mode="lines"` is load-bearing.** Plotly treats a scatter-like trace as
       selectable only when it has markers or text, so the outline is immune to
       the user's own box/lasso -- which writes the very `selectedpoints` the dim
       uses, and which `Choro.dash_selector()` reads back. A programmatic dim and
       a genuine interactive selection cannot coexist on one trace; an outline
       and a selection can.
     - **Applied on `Choro`, not in the verb.** `mosaic` and `plot.animate`
       rebuild from `overlay_traces()` and `get_choropleth()`, never from `.fig`,
       so the old post-render patch was silently dropped by both -- a user got
       four unhighlighted maps and no warning. Both now carry it.
     - **Three inherited bugs fixed rather than preserved.** `select=[]` dimmed
       the whole map (an empty tuple serialises to `[]`, truthy in JS, so the
       branch fired with nothing selected); `select=<geometry>` highlighted
       nothing while dimming everything (`get_vor_cells_as_series(...).to_list()`
       is a list *of lists*, one per feature); out-of-range and negative indices
       were accepted silently. One `_resolve_cells` now serves the plotly and
       matplotlib backends and `plot_mpl(outline_regions=)`.
     - **Resolved eagerly in `__init__`**, so a bad `select=` raises at the call
       rather than from a notebook's display hook cells later.
     - **`select_width` is a module constant, not a 41st parameter.** A contour
       SET has variable density and earns `contour_width`; a highlight boundary
       has one job. `select_name` is derived (region name, file stem, else
       "selection"), which labels the legend correctly for free.
     - **Deferred, and purely additive later:** `select={"streams": [...]}` for
       named groups with per-group colour. `get_vor_cells_as_series` already
       returns that shape, but per-group legend/colour multiplies the surface and
       cannot collide with any v1 spelling.
     - **NOT given to the other picture verbs.** `surface` draws an interpolated
       height field with no cells, `grid` is one trace per cell, and `section`
       already has a `cells=` parameter meaning something else. Recorded because
       the outline mechanism *would* generalise where `selectedpoints` never could.
     - **Also fixed here:** `plot._grid_of` returned a `UsgModel` (`.grid`) and a
       `BuiltModel` (`.context.grid`) *as if they were grids*, so `plot.map(usg)`
       died with `AttributeError: no attribute 'gdf_vorPolys'` several frames from
       the call -- reachable the moment `mf.read_usg` shipped. Non-drawable
       sources now raise `TypeError` naming what was expected.
     - **Left alone, deliberately:** `grid/selection.py:44` *returns* a
       `ValueError` instead of raising it, shared by 20+ call sites across
       readers/geometry/boundaries/particles. `_resolve_cells` sidesteps it by
       coercing `str` to `Path` first. Its blast radius is far wider than this
       change and it is not made worse by shipping the highlight.

161. **Mounding reaches the hover, and says what it is measured from (2026-08-29).**
     `show_mounding=True` coloured the map by mounding but showed no mounding
     number in the hover, and where a number did appear its label did not say
     above WHAT.
     - **Two hover renderers, one of them unfed.** `Choro.hover_dict` carried the
       mounding, but a model map renders through `hover_spec`, which builds from
       the context PAYLOAD and never reads `hover_dict`. And a spec only renders
       fields it NAMES, so feeding the payload was necessary and not sufficient --
       `_resolved_hover_spec` now appends the field with `with_fields`, the same
       mechanism `hover_fields=` uses.
     - **The label now names the datum.** `"Layer 1 Mounding"` became
       `"Mounding above layer 1 bottom"`, and above ground
       `"Mounding above ground surface"`. The old label was not merely vague: at
       `layer=-1` the code measures from the model top but sets `self.layer = 0`,
       so it claimed layer 1's bottom as the datum when the top was used.
     - **`layer=-1` is resolved at construction**, not as a side effect of
       computing the series. Reading the label before the numbers used to report
       the wrong datum, because the above-ground flag was set inside
       `_mounding_series`.
     - **One `_mounding_series` serves both renderers.** They are two renderings
       of one number; computing it in each is how they drift.
     - Adjacent dead code removed while here: `hover_dict`'s layer-elevation
       branch opened by reading `self.vor.gdf_topbtm.columns` into `layer_nums`,
       which was immediately overwritten and never read -- so
       `show_layer_elevs=True` crashed with `'NoneType' has no attribute
       'columns'` on any grid without a layer frame, for a value that was thrown
       away. The model branch takes its elevations from the model, which
       necessarily has them or it could not have run.
     - Related: `read_usg` now publishes the imported layer elevations onto
       `grid.gdf_topbtm` (`UsgModel.attach_layers_to_grid`), so an imported grid
       also feeds the layer hover, the mounding colorscale and the surface-aware
       SFR/LAK builders -- a `.gsf` carries geometry and nothing else.

162. **Horizontal flow barriers, resolved from geometry and validated first (2026-08-29).**
     `mf.hfb` with five entry points. A barrier sits on the FACE between two cells,
     so MODFLOW 6 addresses it as a cell pair — and that is the whole difficulty.
     - **NOT a `package_registry` entry, deliberately.** The registry describes
       cell-indexed boundary conditions; HFB is face-indexed. `record_fields` is
       documented as "the fields AFTER the cellid" and HFB has no cellid — measured,
       `_record_fields_from_dfn("hfb")` raises `ValueError`. `results` would be
       empty because **MODFLOW 6 writes no HFB budget record**: a run with HFB
       produced cbc texts `['CHD','EVTA','FLOW-JA-FACE','RCHA']`, and HFB's only
       trace is a correction to `FLOW-JA-FACE`. A registry entry would have to
       state falsehoods in the one file whose entire premise is that it cannot
       drift. Measured cost of registering it anyway: **12 test failures across 7
       files**; and a *consistent lie* passes the descriptor suite and then dies
       three layers down (`KeyError: 'cell'`, a budget-text `ValueError`, and
       silent `cellid2` corruption in the artifact tier).
     - **The one thing registry membership would have bought** — package discovery
       on reopen — costs ONE line in `run_model._NON_REGISTRY_SUFFIXES`, the same
       route `mvr` takes for the same reason (it describes a relationship rather
       than a cell).
     - **A results tier DOES exist, and this entry's first draft was wrong to say
       otherwise.** MODFLOW 6 writes no HFB *budget record* -- that part holds --
       but the flow across a barrier is not lost: a barrier sits ON a
       cell-to-cell connection, and `FLOW-JA-FACE` carries the flow on every
       connection, so a barrier's flow is simply that entry.
       `packages.hfb.results.q.{get,summary,map}` looks each pair up via the
       model's own `IA`/`JA`, read from the binary grid file MODFLOW 6 writes
       beside its output -- NOT reconstructed from the grid, because `idomain`
       removal makes MODFLOW 6 renumber and a reconstruction would silently index
       the wrong connection. The column is `q_cell1_to_cell2`: MODFLOW 6's own
       sign, with the frame named rather than negated.
     - **The view layer is bespoke, not registry-backed.** The registry view is
       cell-keyed (`split_cellid_columns` knows exactly one cellid) and a
       barrier's geometry is a shared EDGE, so `HfbPackageExplorer` and
       `HfbResultsExplorer` are hand-written, in the same spirit as the
       `StaticArrayPackageExplorer` that serves `ic`/`npf`/`sto`. `map` draws
       barriers as LINES over the cell field rather than filling cells, riding
       the same `add_overlay` path as the `select=` outline of ledger 160, and
       `segments()` returns the faces as a GeoDataFrame -- the shape that
       survives a change of grid, where a cell-pair index does not.
     - **No artifact tier.** The `LIST_BC` branch is measured fatal on HFB twice
       over. Deferred with this entry.
     - **`.enclose` is a separate entry point, not `closed=True` on `.line`.** A
       ring passes THROUGH cells, so the crossed-face set has a gap wherever it
       enters and leaves one: measured on the canonical grid, a square ring crossed
       53 faces and left all 441 cells hydraulically connected, where the cut
       between the enclosed set and its neighbours is 83 faces and seals 105.
       Buffering does not fix it; it is structural, not a tolerance.
     - **Duplicate faces are collapsed, with a warning.** The most dangerous input
       found: MODFLOW 6 accepts a repeated face silently, applies its series
       formula twice (conductance 14.4727 -> 0.8418 once -> 0.4335 twice), and
       `condsat_reset` then restores the ALREADY-MODIFIED value — so the face stays
       wrong for the rest of the run even after an empty period removes the barrier.
     - **`maxhfb` is not exposed.** FloPy computes it from the data at write time;
       a hand-set value that disagrees is the `maxbound` defect class again.
     - **`hydchr` accepted raw only.** A `k=` + `thickness=` convenience pair is a
       natural addition and is deliberately not in v1 — one number, one meaning.
     - **Vertical barriers are accepted, not generated.** Legal on DISV since
       MODFLOW 6 6.7.0; FloPy 3.10's embedded dfn still asserts the old same-layer
       rule and is out of date. The validator accepts an adjacent-layer pair and
       rejects one that skips a layer; the geometry helpers only ever produce
       lateral pairs, because a line is lateral.
     - **Amends ledger 159**: the USG importer's `_hfb_spec` no longer hand-rolls a
       bare `PackageSpec` and no longer drops cross-layer pairs on sight — it calls
       `mf.hfb`, so the Ten Trails model's 28 barriers all convert and are validated.
     - **Left alone**: `grid/selection.py`'s `MultiLineString` gap and its
       `return ValueError`. `barrier_faces` reaches the grid through `sindex.query`
       rather than `get_vor_cells_as_series`, so it never touches them; they remain
       the separate task already filed.

163. **Array-form recharge and ET: `mf.rch.array` / `mf.evt.array` (2026-08-30).**
     MODFLOW 6 offers RCH and EVT in two input shapes and myflopy exposed only the
     list one, so a whole-grid recharge field was written one record per cell.
     Measured on the Ten Trails model (72 periods, 9405 columns): the array form
     is **6.1x faster end to end** (0.88 s vs 5.35 s) and the FloPy constructor
     alone **413x**.
     - **A fourth method on the existing helper, not `mf.rcha`.** The registry's
       one-FloPy-class-per-package assumption does not bind: `flopy_class` has a
       single production consumer, inside `_list_bc_spec`, which `.array()` does
       not route through. An RCHA written as `pname="rch"` runs, reloads and is
       discovered identically -- MODFLOW 6 has ONE recharge package and
       `READASARRAYS` is an option inside its file; FloPy's two classes are a
       FloPy artifact. Siblings would have forced four hand-written namespace
       properties each plus an incoherent `GeoPackageSource.rcha`.
     - **`irch="top_active"` is the default, and that is the whole safety story.**
       MODFLOW 6 with no IRCH array applies recharge to layer 1 unconditionally,
       and a column whose layer 1 is inactive is outside the reduced node
       numbering and skipped SILENTLY. Measured on a 10-column model with 4 such
       holes: list 1000.0, array without IRCH **600.0**, array with IRCH 1000.0 --
       a 40% loss under "Normal termination". `irch=None` is the explicit opt-out.
     - **`irch` is ZERO-based.** FloPy adds one on write. Found the hard way:
       passing a 1-based layer produced `2` in the file and MODFLOW 6 died with
       `Invalid layer number: 3`. The design brief had this backwards.
     - **`nseg > 1` raises rather than silently flattening.** NSEG lives in the
       DIMENSIONS block, which MODFLOW 6 never reads under `READASARRAYS`, so
       segmented ET has no array form at all. Dropping the segments is a change to
       the physics, not the format -- ledger 159 already says so about the Ten
       Trails ETS, which is why it remains a 61 MB list EVT.
     - **`boundnames` raises**: MODFLOW 6 refuses them under `READASARRAYS`.
     - **The inputs tier now reads array packages.** It raised `AttributeError`
       before, because an array package has no `stress_period_data` at all. One
       branch in `build_cell_package_input_table` expands `{per: ncpl-array}` into
       the same tidy frame, which closes the inputs, diff and group tiers together
       -- so **the existing entry at `:2397`** (the array diff asymmetry, accepted
       when RCHA could only arrive from an externally loaded model) is now closed.
     - **The budget record is matched by package NAME, not by substring.** FloPy
       resolves `text=` by first-hit substring and `"RCH"` is inside both `"RCHA"`
       and `"UZF-GWRCH"`, so a model declaring UZF before RCH would have had
       `packages.rch.results.q` quietly return UZF's recharge term. `paknam=` now
       leads, with the text-only match kept as a fallback because FloPy RAISES
       rather than returning empty when a name is absent. Pre-existing hazard;
       arrays made the namespace dense enough to be worth closing.
     - **Deliberately absent**: no `.gpkg` sibling (feature-to-cell mapping is
       inherently list-shaped), no artifact tier (`components.py` raises
       `NotImplementedError` on `package_type == "rcha"` -- deferred), and DISU is
       excluded by MODFLOW 6 itself.
     - **Not a default, and never should be.** The array form needs a value for
       every column including inactive ones; for a sparse boundary (40 cells of
       9405) it is a pessimization -- ~1 kB as a list against ~160 kB as an array.
     - **Amends ledger 159**: the USG importer's `_recharge_spec` calls
       `mf.rch.array` instead of hand-rolling a bare `PackageSpec`, and resolves
       `irch` from the model's own idomain, so `NRCHOP = 1` and `NRCHOP = 3` are
       both honoured rather than approximated by MODFLOW 6's layer-1 default.
     - **Still open**: the PEST `_ALIASES["rcha"] -> "recharge"` entry is stale and
       an array recipe for `parameterize("recharge")` is not written.

164. **`VoronoiGridPlus.from_gsf` reads the `.gsf` itself rather than reusing
     flopy's `UnstructuredGrid.from_gridspec` (2026-08-28).**
     A MODFLOW-USG `DISU` stores connectivity and geometric measures but no
     coordinates, so the `.gsf` is the only file that knows where the cells are.
     - **Why not delegate.** flopy's reader rejects the standard
       `UNSTRUCTURED GWF` header outright -- `not (A) or (B)` where `not (A or B)`
       was meant, so the two-word form always raises "Invalid GSF file, no
       header" -- and it returns all layers stacked with the raw 3-D vertex
       lists, which is precisely what must not reach a plan grid. Wrapping it
       would mean patching a third-party parser and then undoing its output.
       ~90 lines of our own parser buys the header tolerance, the layer
       selection, and clear errors that name the offending file and node.
     - **Why the collapse matters.** A `.gsf` cell is a hexahedron: 8 vertices,
       the bottom four under the top four. Passed through as-is, every polygon
       traces its outline twice -- `is_valid` False, DOUBLE the true area -- and
       the map still *draws correctly*, so the only visible symptom is that
       hover stops working, because maplibre abandons hit-testing on a
       self-overlapping ring. Measured in headless Chromium: 0/5 hover probes on
       the degenerate grid, 5/5 on the collapsed one, with the two renderings
       within 0.4% of each other on colored pixel count.
     - **`from_gsf`, not `vor_from_gsf`.** Its neighbour is `vor_from_disu`, so
       the two constructors on this class now disagree about the `vor_` prefix.
       `from_gsf` matches the `from_*` convention used everywhere else in the
       codebase (`Surface.from_contours`, `LayerStack.from_modflow`,
       `GridSpec.from_object`); renaming `vor_from_disu` to match would break the
       two Cumberland notebooks that call it, so the inconsistency stays until
       something else touches that method.
     - **NOT done: layer elevations.** A `.gsf` carries z per vertex, which is
       enough to populate `gdf_topbtm` and light up the layer-elevation rows in
       the map hover. `from_gsf` drops z entirely and returns a plan view only.
       A USG model's tops and bottoms are also in its DISU, so the better home
       for that is a reader that takes both files, not this one.
     - **NOT done: idomain from the USG BAS.** `from_gsf` takes `idomain` by
       hand like every other constructor here. Note flopy's `MfUsgBas.load` sizes
       IBOUND from DISU alone and has no CLN awareness, so it under-reads the
       array on any model with CLN nodes -- reading idomain automatically would
       mean parsing the BAS ourselves.

165. **Splitting a pinched-out unit inflates it by `(N-1) * min_sep` — fix DEFERRED
     (2026-08-30, from 158).**
     Splitting is supposed to change discretization, not geometry, and for a unit
     with real thickness it does (158 scaled `trigger_sep` to make it exact).
     A unit that pinches out **entirely** is the one case left: reconcile floors
     each sub-contact at `min_sep` independently, so the unit's total grows
     linearly with the number of slices and everything below it shifts down.
     - Re-measured 2026-08-30, still live, on a zero-thickness unit with
       `pinch="inactive"` and the default `min_sep=0.1`:

       | split | unit total | base |
       |---|---|---|
       | none | 0.100 | 99.900 |
       | 2 | 0.200 | 99.800 |
       | 3 | 0.300 | 99.700 |
       | 5 | 0.500 | 99.500 |

     - **Why it is tolerable:** the cells are `idomain = 0` by construction, so
       nothing is solved in them; the drift is bounded and known; and it scales
       with `min_sep`, so a caller who cares can shrink it. It cannot silently
       grow — the unit either pinches whole or not at all (158).
     - **Why it is not fixed:** both routes are wider than the change that
       surfaced it.
       1. **Per-layer `min_sep`.** `reconcile_surfaces` (`grid/geometry.py:332`)
          takes a scalar and applies it column by column; it is also exposed as
          `vor.reconcile_surfaces`, so widening the parameter is a public
          signature change with its own snapshot review.
       2. **Collapse a pinched unit's sub-contacts before reconcile.** Cleaner in
          principle — the split provably could not move geometry — but reconcile
          is single-pass by design, and the pinch verdict is currently computed
          *after* it, so this reorders the engine and re-opens
          `_validate_pinch_invariant`.
     - **If picked up:** route 2 is the better answer and subsumes route 1 for
       this case. Do it with the two-pass shape the 158 design discussion
       rejected on scope grounds (reconcile the units as declared, then cut each
       unit's interval), not by special-casing zero thickness. Test it against
       the table above: every row should read `0.100 / 99.900`.

166. **`model.packages.summary()` / `.mosaic()`, and three bugs found by asking
     for them (2026-08-30).**
     Asked for a quick way to see every mappable package and its key facts.
     There was none: `model.package_names` gives names, `model.summary()` is one
     model-level row, and `ModelPackages` had no `__iter__`, `__dir__` or
     `__repr__`. The registry already knew every package's mappable fields;
     nothing exposed them together.
     - **`summary()` is CHEAP by default** -- registry plus the package list, no
       package data read -- because `LoadedMf6Run` overrides package discovery
       precisely to avoid `_ensure_core_loaded()`, and a rich default would have
       thrown that away. `detail="data"` opts into `records`/`periods`/`layers`.
       The `layers` column is the one that pays for itself: it reports `WEL` as
       `1, 3`, which is exactly why `wel.inputs.map()` drew nothing at the
       default `layer=0`.
     - **`mosaic()` is a verb, not a flag on `summary()`.** A `summary(mosaic=True)`
       would be the table-and-picture conflation `view_layer_conventions` exists
       to prevent. `ModelPackages` is not a policed scope in
       `test_plot_vocabulary` (only module/model/grid/stack are), so neither
       method has a vocabulary cost.
     - **`mappable` is PROBED, not inferred from the registry.** Mutation testing
       caught this: HFB is deliberately absent from `package_registry` (ledger
       162) yet `packages.hfb.inputs.map()` works, so a registry-derived
       predicate reported the one hand-written package as undrawable -- exactly
       backwards. The column now looks for a callable `map()`.
     - **Bug: every map died on an unrun model.** `_grid_of` probes
       `hasattr(source, "hds")`; `hasattr` swallows only `AttributeError`, so
       FloPy's `FileNotFoundError` escaped a CAPABILITY PROBE and killed even
       INPUT maps, which need no results at all. Now `_has_results` catches the
       closed set (`OSError`/`ValueError`/`KeyError`) and logs the degradation.
     - **Bug: `Choro.per` dereferenced a model it does not always have.** With
       the probe fixed, input maps arrived with `model=None` and died in the
       `per` setter. For a model-less choropleth `per` is a LABEL -- which
       period's records were selected is already baked into the values -- so it
       is stored rather than resolved to a `kstpkper`.
     - **Bug: an empty selection reported the wrong thing.** A selection matching
       no records comes back WITHOUT its value column, so guarding the column
       before the emptiness turned "no records for this period/layer" into
       `KeyError: Value column 'q' was not found`, naming a column the package
       certainly has -- and made the all-fill branch below it unreachable. Found
       in BOTH `build_cell_input_map_payload` and
       `build_group_input_compare_map_payload`; fixed in both.
     - **`package_summary()` deleted outright, not deprecated.** Both bodies
       (`SimulationBase` and the `LoadedMf6Run` override). It is absent from
       `tests/api_snapshot.json` -- which DOES cover this class, so the absence
       is meaningful -- absent from `__all__` and `__compatibility__`, and had
       zero callers anywhere. `docs/deprecation_policy.md` exists to retire a
       name "without breaking existing scripts"; there were none, so the
       machinery would have cost an alias, a compatibility entry, a policy row
       and a guard test to protect nobody. Recorded as a deliberate departure
       from the two-release rule, and as the precedent 167 relies on.

167. **The rest of the `*_summary()` family, and `model.outputs` -- PART DONE,
     rest DEFERRED (2026-08-30).**
     **Update, same day:** `output_summary` and `result_summary` are now DELETED
     too -- the zero-caller half of the family, taken because it cost nothing and
     removed a GWF-only trap (both hardcoded `f"{self.name}.hds"`, so on a
     GWT/GWE model they reported `has_heads=False` and null statistics in
     silence). `result_summary` needed no replacement built:
     `model.hds.summary()` already existed, is already tested, is kind-dispatched
     through `model._field_reader`, and reports MORE (`records`, `periods`,
     `layers` on top of the same statistics). `output_summary` needed none at
     all -- it welded "which output files exist" to "how far did the run get",
     which is why it had no single natural home.
     **Still open: `file_summary` and `grid_summary`**, because those two are the
     ones whose destinations do not exist. Everything below still stands for
     them.

     `package_summary` (166) was one of six. The other five are the last
     surviving instance of the prefix shape `view_layer_conventions.md:272`
     forbids, and the intended end state is nouns:
     `model.results.summary()`, `model.files.summary()`, `model.vor.summary()`.
     - **Measured cost of finishing it: three lines in two notebooks.**
       `output_summary` and `result_summary` had ZERO callers anywhere (DONE);
       `grid_summary` has one (`usg_import_workbook.ipynb:620`); `file_summary`
       has two (`usg_import_workbook.ipynb:610`, `usg_mf6_model.ipynb:443`).
       None is in `api_snapshot.json`; none has a test. Both notebooks are the
       user's live USG workbooks, so they are edited surgically, never
       regenerated.
     - **Why not now: the destinations do not exist.** `model.results` and
       `model.files` are both absent (`hasattr` False; the names are free).
       `model.vor` exists as a plot scope but has no `summary()`. And a
       recursive `model.results.summary()` has a hole today --
       `lak.results.summary()`, `sfr.results.summary()` and
       `surface_water.results.summary()` are all `AttributeError`, where the
       eight registry-backed namespaces have it.
     - **Three defects the replacements must not inherit.** `_OUTPUT_SUFFIXES`
       lists `".obs.csv"` but the test is `path.suffix`, and
       `Path("a.obs.csv").suffix == ".csv"` -- so that member has never matched,
       and all 12 output CSVs on the canonical model are labelled
       `category="input"` (the honest split is 35/12, not 41/6);
       `list_input_files()`/`list_output_files()` share it. `output_summary` and
       `result_summary` both hardcode `f"{self.name}.hds"`, so on a GWT/GWE model
       they silently return `has_heads=False`, where `model._field_reader`
       already dispatches correctly. And `result_summary` reads the LAST time
       step only, so its head range disagrees with `model.hds.summary()` over all
       periods for reasons no column explains.
     - **`model.outputs` is NOT what it was assumed to be.** The belief that it
       was being replicated under `model.packages.<pkg>.results` is half right in
       a misleading way: the `.stage` half was superseded by `.results` (and is
       richer there -- tidy frames rather than raw ndarrays), but the `.bud` half
       went to the SIBLING namespace `packages.<pkg>.budget.<term>`.
       `packages.lak.results.bud` is an `AttributeError`. A cleanup driven by the
       `.results` framing would conclude the raw budget accessor has no
       replacement and either keep `outputs` or duplicate `.budget` under
       `.results`.
     - **Nothing in the repo records that intent.** Searched all of `docs/`,
       CLAUDE.md, the ledger, the plan, every docstring and the full git log:
       no entry, commit or docstring says `model.outputs` was to be superseded.
       What IS written points the other way -- `budget.py:312` and `:426` both
       call `model.outputs.lak.stage` "the canonical home for this behavior",
       and those two shims (`LakStage`, `SFRStage`) are themselves dead.
     - **It cannot be deleted, only demoted.** 18 call sites inside `src/`, and
       the load-bearing ones BUILD the tier that replaced it --
       `package_budget.py` constructs `packages.<lak|sfr>.results.stage/q`
       through `outputs.<pkg>.stage`, and `package_surface_water.py` calls
       `model.outputs.lak.bud.types` from inside `LakResultsNamespace` itself.
       It is also pinned live in `api_snapshot.json` and `test_typing_surface.py`
       and taught without a deprecation marker in
       `docs/package_api_reference.md`.
     - **One capability has no replacement at all:**
       `outputs.uzf.ifno_to_cellid`. It is the only member that works on an
       UNRUN model -- it reads `uzf.packagedata` plus the modelgrid -- so it is
       an INPUT index map misfiled under "outputs", and its home is
       `packages.uzf.inputs`. It is also the only part of the namespace any test
       asserts on.
     - **Free wins found alongside, not taken here:** `model.lak_output`,
       `model.sfr_output`, `model.uzf_output` and `group.outputs` have zero
       readers and exist in the snapshot only because the derive script sweeps
       every public property. And `calibration.py:740` calls
       `model.outputs.lak.stage.nlakes`, which already raises -- `nlakes` lives
       only on the dead `LakStage` shim.
     - **When picked up:** the two zero-caller summaries first (free, and it
       removes the GWF-only trap), then `vor.summary()` -- deciding on the way
       whether `nlay`/`node_count` belong on it at all, since a 2-D grid cannot
       know them and the narrower-scope rule says fewer members, never other
       ones. Then `model.files`, name-based rather than suffix-based. For
       `outputs`: delete the four zero-caller names, rehome `ifno_to_cellid`,
       correct the live doc, and only then demote the property behind
       `__compatibility__` -- do NOT reuse the word `outputs` for the file
       namespace.

168. **An input map's hover shows the whole record, not just the coloured field
     (2026-08-30).**
     Reported against `ghb.inputs.map()` and `drn.inputs.map()`: the hover gave
     the cell number, the one coloured field, and the period. Reading a GHB map
     means asking "what head, against what conductance", and that answered half
     of it.
     - **The data was already there.** `build_cell_input_map_payload` aggregates
       every field in the selection -- the payload for a GHB map already carried
       `bhead`, `cond`, `Layer`, `Record Count` and `Package`. `cell_input_hover`
       already took an `extra_fields` argument and built a `Fields` block from
       it. The three call sites in `package_inputs.py` simply passed
       `cell_input_hover(self.field_name)` and nothing else, so the rest was
       computed and dropped.
     - The sibling fields render as an untitled block, matching `lak_hover` /
       `sfr_hover`, which list a feature's fields the same way.
     - **`Layer` is always shown; `Record Count` only when records merged.** The
       layer is worth saying on every map -- the commonest reason a boundary map
       comes back empty is that the package has no records in the default layer 0
       (see 166). The record count is only informative when some cell aggregated
       more than one record, which is the one thing that explains a summed value;
       measured on the canonical model every package is 1 per cell, so showing it
       unconditionally would read "records 1" forever. The call site inspects the
       payload and asks for it only when `max > 1` -- data-dependent composition
       at the CALL site, not in the spec.
     - **No change to what drives the colour**, which was already right:
       `field=` selects it and the default comes from the registry's
       `default_input` (`bhead` for GHB, `elev` for DRN, `stage` for RIV). Only
       the hover was wrong.
     - Wired at all three input explorers -- cell-stress, UZF and static-array.
       `UzfFieldInputsExplorer` had no `package_name`, so it gained one as a class
       attribute. UZF and NPF payloads carry a single array, and `Fields.build`
       skips names the payload lacks, so no empty block renders there.

169. **`usg.export_gis` — the USG model as grid-independent GIS (2026-08-30).**
     `to_mf6()` binds an imported model to the grid it converted on, which is the
     wrong unit of reuse: the reason to import a USG model is usually to rebuild it
     on a better mesh, and the mesh changes again. `export_gis` writes the content
     instead. The scope cuts, and the judgment calls, in order:
     - **It builds no packages.** It writes files whose columns are the builders'
       parameter names and stops. An earlier scope (a `CellOverlay` primitive plus
       per-field reducers, ~7 days) was designed and then dropped in favour of this
       on the user's direction, and the direction was right: a transfer has to be
       re-run and re-argued on every re-grid, an extract does not. The overlay work
       is NOT deferred-with-intent — it is unnecessary unless someone wants an
       automatic grid-to-grid path later.
     - **Conductance conserves the feature total, not a leakance.** Measured over
       205 DRN records `corr(cond, cell area) ~ 0`, and GHB is a single constant per
       family (172,357.9 in layer 3, 0.15 in layer 5) across cells spanning 2.8 to
       697,010 ft2 — a 249,000x range. These are calibrated numbers with no geometry
       in them, so `sum C` is the only thing worth preserving. The `*_per_ft` twin
       does that and nothing more; it is deliberately NOT presented as a leakance.
     - **Rasters lose the smallest cells, and the loss is measured rather than
       warned about.** Default resolution is a quarter of the median cell width
       (17.2 ft on Ten Trails); cells below ~295 ft2 are absorbed. Sampling `top.tif`
       back at the cell centres — the same operation a consumer performs — costs at
       most 5.04 ft and typically 0.240 ft. `cell_values.gpkg` carries every array
       exactly, and the manifest says so with the numbers.
     - **`bedleak = FSKIN / 1.0 ft` is a choice, not a conversion.** MF6's `bedleak`
       is a leakance (1/T); `FSKIN` is a conductivity (L/T). CLN has a *skin*, not a
       bed, so no thickness exists to divide by. `fskin` is written unmodified beside
       the derived column and `bed_thickness=` is a parameter.
     - **`mf.lak` takes one `bed_leakance` per lake**, and `FSKIN` varies 9x within
       Horseshoe Lake (0.054–0.498). The median is written to `lakes`, the variation
       to `lake_nodes`, and the manifest names every lake whose spread exceeds 2x.
       Teaching `LAKBuilder._leakance` to accept a cell-indexed Series — which
       `_bottom` right beside it already does — is a real asymmetry, NOT fixed here.
     - **SFR `rbth` and `man` are invented** (1.0 ft, 0.03). CLN routes with a pipe
       conductivity, which has no MF6 counterpart at all. Written as documented
       defaults rather than omitted, so a missing column never becomes a silent zero.
     - **`mf.sfr`'s `reach_top` still cannot take a raster or a line's own Z**, so a
       streambed profile cannot cross a re-grid through the builder. `_reach_values`
       accepts a scalar, a column, or a mapping keyed by reach number — and reach
       numbers change with the grid. `mf.lak`'s `lake_bottom` already samples a
       raster. The export routes around it (Z on the line, plus `station` on the
       nodes for `np.interp`) but the asymmetry is real and NOT closed. Note the
       fallback is actively wrong here: `reach_top=None` samples the *model top*,
       which on Ten Trails' CrispCreek is up to 160 ft above the bed.
     - **A single-cell stream gets no line.** Its node still reaches `stream_nodes`
       with its real bed elevation and `FSKIN`; inventing a centerline through one
       cell would put a reach somewhere nobody chose. Named in the manifest notes.
     - **The stream profile is nearly redundant, and is written anyway.** Measured,
       the CLN bed IS the old model top to 0.00003 ft on all three creeks — so it
       carries no independent elevation data. But that top was *conditioned* monotonic
       along the channels (0 rises in 97 nodes, where the raw new DEM has 17% at
       300 ft sampling and 38% at 50 ft; p ~ 0.001 that 0/38 is chance). So the
       profile's value is as the conditioned shape to aim at, not as bed elevations.
       No conditioning helper ships — that is the user's workflow, and a swath search
       confirmed there is nothing to snap to (a 200 ft-wide minimum finds ground only
       0.24–0.45 ft lower, so the rises are DEM roughness, not misplacement).
     - **Features are grouped by name, breaking the old behaviour on purpose.**
       `cln_polygons()` previously grouped by connected component, which merges a
       tributary into its trunk — 6 features where the file names 7. Names win when
       present; shape is the fallback. `ClnFeature.kind` still comes from mean degree,
       so a 2-node "lake" is correctly classified as a chain.
     - **`free_row` replaces positional number-scanning for labelled rows.**
       `Wlnd217` contributes a phantom `217.0` to `free_floats`, giving that row one
       more number than its neighbours. Harmless while only indices 1/3/4 were read;
       fatal the moment anything indexes from the end. Same class as the
       `IPRN`-as-multiplier bug that once negated an entire ET surface.
     - **The CLN carries no bathymetry, and saying so is the deliverable.** Measured,
       EVERY feature's bed is the model top -- lakes as well as streams, to 0.00003 ft
       (Keevie 313/313 nodes at top, Black Diamond 150/150, Horseshoe 140/140,
       Marjorie 104/104). So a LAK built from these bottoms has zero depth, and worse,
       312 of Keevie's nodes round to just ABOVE their cell top, which `mf.lak` rejects
       as "bottom does not intersect active cell". This is reported in the manifest
       rather than papered over: inventing a depth is the user's call, not the
       exporter's. Found by round-tripping the export back through `mf.lak`, which is
       why that round trip is worth doing on any future exporter.
     - **A dissolved footprint touches cells the CLN never had.** `LAKBuilder` resolves
       lake cells from the polygon, so a lakebed raster covering only the CLN's own
       cells leaves it asking for a bottom that is not there (measured: 249 fringe
       cells over four lakes). They are filled from the nearest node. The alternative
       -- shrinking the polygon -- would drop real lake area.
     - **Three of RockCreek's 39 cells are inactive in every layer.** MODFLOW-USG
       accepted the CLN-GWF connection; MF6's SFR refuses to place a reach there. The
       `cell_active` column on `stream_nodes`/`lake_nodes` and a manifest note carry
       it. NOT filtered automatically -- which reaches to drop is a modelling choice.
     - **The workflow notebook builds no model, on instruction.**
       `examples/mf6/notebooks/usg_export_to_new_model.ipynb` shows the call that
       consumes each written file and stops there. It defaults to `USE_DEMO_GRID =
       True`, loading the OLD USG mesh, so it runs end to end out of the box and the
       shapes are visible before a real grid is committed; the `mf.sfr`/`mf.lak`
       calls are left commented because they need an `nper` and a layer decision the
       notebook must not make. Its `line_bc_records` helper PRINTS what the layer map
       drops -- with old layer 5 out, the whole Qpon GHB family (67 records, the
       model's only deep boundary) goes with it, and a silent drop there is exactly
       the failure this export exists to prevent.
     - **`mf.Spread` (2026-08-31), because the canonical path was wrong for a
       conductance.** `GeoPackageSource._boundary_data` writes a feature's value to
       EVERY cell it intersects. Right for an elevation or a head; for an extensive
       field it multiplies by the cell count. Measured through `mf.drn.gpkg` on the
       grid the values came from -- not a regrid -- 205 line features became 424
       records and 66,007.43 ft2/d became 133,631.25, **a factor of 2.02**. Wrapping
       the field (`conductance=mf.Spread("conductance")`) hands each cell its share of
       the feature, by length for a line and area for a polygon. On the user's real
       18,107-cell grid: 864 records over 209 cells, 65,997.5 ft2/d -- 99.98% of
       source, continuous coverage where the original had 68 cells.
       **Opt-in, deliberately**: the default is unchanged, so no existing `.gpkg`
       caller's results move. `mode="clip"` (default) lets the part outside the grid
       go, `"retained"` renormalizes; `min_share=0.01` drops corner clips (measured
       11 of 79 cells holding 0.90% between them). A point keeps its whole value --
       nothing to divide -- so the point layers stay exact without it.
     - **Boundaries are written twice, as points AND lines.** The lines are the build
       path (continuous cells, needs `Spread`); the points resolve to one cell each so
       an unwrapped conductance is already exact there. Keeping both is ~40 KB and
       removes a footgun; naming one "BUILD FROM THIS" in the manifest is what makes
       the choice visible rather than a coin flip.
     - **Column names are the package registry's, not USG's** -- `elevation`,
       `conductance`, `head`, so `mf.drn.gpkg(path, layer="drn_lines", ...)` needs no
       field arguments. CHD's `ehead` has no MF6 counterpart and is written as
       `ehead_usg_only` rather than dropped.
     - **Three defects found by actually writing and running the rebuilt model
       (2026-08-31), each of which produced input MF6 rejects or mis-solves:**
       (a) `Spread` had no `_metadata_value` branch, so `prepare_run` died with
       `TypeError: Object of type Spread is not JSON serializable` -- long after the
       package built cleanly. Every dataclass in the `RowValue` union needs one, and
       a test now asserts the whole union round-trips.
       (b) `GeoPackageSource._active` used `bool(idomain[layer, cell])`, and
       **`bool(-1)` is True**. MF6 idomain is three-valued -- `>0` active, `0`
       inactive, `<0` vertical passthrough -- and a passthrough cell holds no
       boundary. Measured: 347 of 864 DRN records landed in passthrough cells on a
       stack using `pinch="passthrough"`, failing the run at read time. Now `> 0`.
       (c) `SFRBuilder` kept any reach with `length > 0`, so a stream clipping a cell
       corner produced a 0.005 ft reach; MF6 divides by reach length. `mf.sfr` gains
       `min_reach_length` (default 0.0, no behaviour change). NOTE: this did NOT fix
       the SIGFPE it was suspected of -- that remains open, see below.
     - **`mf.tdis(nper=)` does not infer `nper` from `perioddata`** and defaults to 1.
       A 72-row perioddata with the default wrote `NPER 1`; MF6 would have solved one
       31-day period and terminated normally. Not changed -- a mismatch is arguably a
       caller error -- but it is the single easiest way to silently get a steady
       answer from a transient model, and it deserves a guard.
     - **A raster export does not cover a larger new domain.** 2,287 of 18,107 cells
       fell outside the old model, and `Surface.raster(...).values(vor)` returns NaN
       there; one NaN makes every budget term NaN with no warning. The consuming
       notebook fills gaps explicitly rather than the export inventing values.
     - **OPEN: the rebuilt model SIGFPEs at the first timestep.** All 12 packages
       read; the crash is arithmetic during the solve, after LAK setup. Ruled out by
       measurement: coordinate magnitude (shifting DISV to a local origin changed
       nothing), NaN in K/K33/Ss/Sy, zero or negative layer thickness (min 0.100 ft),
       and sliver SFR reaches. Not diagnosed further -- it is model physics, not the
       transfer.
     - **`min_thickness` now sets geometry, and `pinch` defaults to `"floor"`
       (2026-09-01, approved by the user, who confirmed no existing models are
       affected).** Before this, `min_thickness` only ever decided an idomain value;
       the built thickness came from the global `min_sep`. So
       `.add(name, min_thickness=5)` produced a 0.1 ft layer, and `LayerStack.add`'s
       own docstring described a clamp the code never performed (`layers.py:2000`
       and `:2268` said `"floor"` clamps; `_idomain_from_thickness` left the
       geometry untouched and only skipped the idomain change).
       Now reconcile takes a **per-layer** separation: a `"floor"` layer is spaced at
       its own `min_thickness`, and `min_sep` is the fallback for layers that pinch,
       which must still be allowed to come out thin -- being thin is the signal that
       they pinch out. `reconcile_surfaces` accepts a scalar or one value per bottom
       contact.
       Three consequences worth carrying: the split-aware division (a unit's minimum
       is divided among its slices, or a 3 ft unit split three ways builds 9 ft --
       the geometry twin of the ledger-158 pinch bug); the guard on
       `min_sep < min_thickness` still stands for pinching layers and now names all
       three ways out; and the default flip means thin cells stay ACTIVE, which on
       the Ten Trails rebuild took 76,738 active / 13,797 passthrough to 90,535 / 0.
       That matters beyond tidiness -- a passthrough cell carries no boundary, and
       it had silently swallowed 347 DRN records.
       Five tests were updated to the new default and five added for the new
       behaviour; `test_per_layer_min_thickness_override` now passes
       `pinch="passthrough"` explicitly, because a threshold only reaches idomain
       under a pinching policy.
     - **`to_disv` resolved the geometry a SECOND time and disagreed with `build`
       by 8.5 ft** -- found by the user, whose drains landed below their cell
       bottoms. `to_disv` already called `build()` for the idomain (the ledger-158
       fix) but then passed the SCALAR `min_sep` to `ls.to_disv`, so the per-layer
       separations reached one resolution and not the other. The failure is
       maximally quiet: `attach_to_grid()` publishes `build`'s surfaces, so
       `mf.CellSurfaceOffset("cell_bottom", offset=2)` placed boundaries two feet
       above bottoms **MODFLOW never saw**, and MF6 then rejected them against the
       bottoms it did. Same class as ledger 158, one field over -- the lesson is
       that ANY resolution argument `to_disv` does not share with `build` writes a
       different model than the one you inspected.
       `test_to_disv_and_build_resolve_identical_geometry` pins top, botm AND
       idomain; `test_cell_surface_offset_sees_the_reconciled_bottoms` pins the
       consequence end to end.
     - **`mf.tdis(ats=...)` (2026-09-01)** -- adaptive time stepping reached the
       package-first path. The engine (`_resolve_ats_periods`, `_build_ats_records`
       in `mf6/simulation/discretization.py`) and its tests already existed; only
       the legacy `TemporalDiscretization` could reach them, so the canonical
       `mf.tdis` accepted raw FloPy `ats_perioddata` and nothing else. Same shape of
       gap CLAUDE.md records for other packages: engine complete, facade missing.
       `ats=` takes `True`, an iterable of **zero-based** period indices, or a
       mapping of per-period overrides; `ats_perioddata=` stays as the escape hatch
       and passing both raises rather than silently choosing.
       **It also fixes a live defect for existing `ats_perioddata` users:** FloPy
       does not size `MAXATS` from the record list, so the written file declared
       `MAXATS 1` however many records were supplied and MF6 read only the first.
       Measured: 4 records in, `MAXATS 1` out. The legacy path corrected this via
       `_set_maxats`; the spec path did not. `SimulationSpec.build_flopy` now sizes
       it for every tdis carrying records, whichever way they were built.
       **The two ATS helpers stay function-level imports.** Hoisting them to module
       level in `package_api`/`specs` imports cleanly but leaves a partially
       initialised graph -- 10 unrelated tests fail, including the noun-signature
       suite. The deferred-import ratchet was regenerated instead (`package_api`
       2 -> 3, `specs` 2 -> 3), which is the correct call when hoisting is what
       breaks rather than what fixes.
     - **Not exported:** HFB (face-indexed, no line yet), the `pxdp`/`petm` ET segment
       arrays (uniform scalars on this model — 0.3 and 1.0 — so they belong in the
       call, not a file), and observation/PEST scaffolding. `icelltype` is constant so
       it is not written as a raster.

170. **`GridSpec.resolve` documented rather than repaired; the CRS fallback left in
     place (2026-08-31).**
     Asked for thorough docs on `resolve()` -- "it's not clear what the args are".
     The four arguments are now a full numpydoc block on the method, plus a short
     section in `model_building_cheatsheet.md` and a line in
     `package_api_reference.md`. Two findings surfaced while tracing them, and
     both were documented rather than changed:
     - **`_voronoi_options` falls back to `crs="EPSG:2927"`** when neither the
       spec nor its options name one -- Washington State Plane South, feet, a
       silent regional default for every user everywhere else. Changing it now
       (to `None`, or to raising) would break any existing spec relying on it,
       and this pass was scoped to documentation. It is called out as a trap in
       all three places instead. **A real fix is a candidate deferral**: require
       `crs=` on `GridSpec.voronoi`, with a deprecation cycle for the default.
     - **`resolve()` is not cached** and re-runs Triangle on every call,
       overwriting `workspace`. Fixed file names (`_triangle.0.poly` and
       friends) mean two grids sharing a workspace clobber each other's meshes.
       Documented as "one directory per grid" rather than fixed, because
       namespacing the files is a `TriangleGrid` change with its own blast
       radius and the failure is loud, not silent.
     - **`return_triangle` is silently ignored** by the `object` and `pickle`
       methods (they return before the check that raises for `python` specs).
       Left as-is: the argument is meaningless there and raising would make
       `from_object` specs harder to swap in as a drop-in replacement, which is
       their whole point.

171. **`Project.grids` typed so an editor can find `resolve`; `resolve() -> Any` left
     alone (2026-08-31).**
     Reported as "PyCharm shows me `async def _resolve_with_query(...)` when I hover
     over `resolve`" -- aiohttp's DNS resolver. Not a myflopy method at all: with
     `self.grids: dict[str, Any]`, the receiver has no type, so an editor falls back
     to guessing among every `resolve` in its index. Now
     `dict[str, GridSpec | GridRef]`; `add_grid` and `_load_grid` were `-> Any` too.
     Verified with `mypy`, which now reveals `GridSpec | GridRef` where it revealed
     `Any`.
     - **The union forced a decision about `GridRef`.** `GridRef` had no `resolve`,
       so the honest union produced a *new* warning on the user's working line
       (`Item "GridRef" ... has no attribute "resolve"`). Three options: annotate
       `dict[str, GridSpec]` (a lie -- `_save_grid`/`_load_grid` both handle a
       `GridRef` entry and it round-trips); leave the warning; or give `GridRef` a
       `resolve` that raises. Took the third. A method that always raises is
       unusual, but a reference genuinely cannot build anything, and the message
       now names the entry to look up instead of arriving as a bare
       `AttributeError`. Its signature is asserted equal to `GridSpec.resolve`'s,
       or narrowing the union would buy nothing.
     - **`crs: str | None` widened to `str | int | None`** (5 sites in `specs.py`,
       3 in `sources.py`). `crs=2927` is the natural spelling for an EPSG code and
       works at runtime everywhere it is used -- it was simply a type error nobody
       saw, because the `Any` above meant no checker ever reached the call.
     - **REFUSED: `@overload`s to narrow `resolve() -> Any`.** The return depends on
       the spec's *method* -- `object` hands back its grid, `pickle` unpickles one,
       `python` returns whatever a user's builder made, `voronoi` returns a
       `VoronoiGridPlus` (or a `TriangleGrid` under `build=False`, or a tuple under
       `return_triangle=True`). Method is runtime state, not in the type, so
       overloads on `build`/`return_triangle` would be accurate for `voronoi` specs
       and would LIE for the other three. `Any` is the honest annotation; the
       docstring's Returns section carries all five shapes instead. A real fix means
       splitting `GridSpec` by method, or a separate `resolve_voronoi()` narrow
       enough to promise a type -- both wider than this change.

172. **Grid refinement documented across three channels; `GeoPackageSourceSpec.fields`
     was documented backwards (2026-08-31).**
     Asked twice in one session how refinement works with several geopackages and
     with mixed polygon/line layers -- the behaviour was correct and entirely
     undocumented. Now a table + section in `model_building_cheatsheet.md` and a
     rewritten `GridSpec.voronoi` docstring. Every claim was measured on synthetic
     fixtures, not read off the source.
     - **The `fields` docstring described a class that does not exist.** It claimed
       a "source column name -> target field name" mapping recorded by
       `mf.ghb.gpkg(...)`, with the example `{"head": "bhead", "k": "cond"}`.
       Grepped: `DataSourceSpec.fields` is consumed ONLY by
       `grid_spec_resolver._source_value`, and only as `{logical_key: column}` --
       the opposite direction. The runtime `GeoPackageSource` behind `mf.ghb.gpkg`
       has no `fields` at all; it uses `name_field`/`layer_field`/`period_field`.
       Corrected, with the direction called out explicitly, because the refinement
       docs tell users to write `fields={"area": "max_area"}` and the class doc
       said to write it the other way round.
     - **The three channels are asymmetric and stay that way.** `refinement=` takes
       ONE source, `breaklines=` a LIST, `points=` a LIST; buffers default 0 / 10 /
       n-a and priorities 0 / 1 / n-a. Documented rather than harmonized: the
       asymmetry is load-bearing (a line has no area, so it MUST be buffered; a
       polygon must not be), and `breaklines` with `breakline_buffer=0` is already
       a working multi-source polygon channel. Measured: a polygon routed that way
       arrives at its exact area with its own `max_area`.
     - **DOCUMENTED, NOT FIXED: `refinement=[a, b]` fails as
       `AttributeError: 'list' object has no attribute 'path'`.** A type check
       naming the two real options (merge the layers, or use `breaklines` with
       buffer 0) would be three lines and strictly better. Left out because this
       pass was scoped to documentation and the user asked for docs; the cheatsheet
       quotes the raw error so it is at least searchable. **Worth doing.**
     - **Two silent behaviours worth knowing, both now documented:** an unrecognised
       `fields` key is ignored (a typo'd `"aera"` silently falls back to the global
       option -- verified), and `refinement_buffer` applies to every feature in the
       source, not only the lines that need it (a 400x400 polygon came out 67%
       larger). Neither is a bug; both are invisible without being told.

173. **`GeoPackageSourceSpec` given a full docstring; two more descriptions were wrong
     (2026-08-31).**
     Follow-on from 172 -- the class the refinement docs point at listed its
     arguments thinly and defined none of them. Now a complete numpydoc block
     (every field including the inherited `external`/`metadata`, a keys table for
     `fields`, Notes, See Also, four runnable Examples). Verified by constructing
     and reading every Example and round-tripping each through `to_dict`.
     - **`query` is `pandas.DataFrame.query`, not SQL.** Documented as "an OGR/SQL
       `WHERE` expression"; `_read_source` reads the whole layer then calls
       `gdf.query(...)`. Measured: `"active = 1"` raises `ValueError: cannot assign
       without a target object`, `"active == 1"` returns 1 of 2 features. The old
       docstring's own example was the broken spelling -- so was the replacement
       example added earlier in this same session, which is a good argument for
       running docstring examples rather than eyeballing them. Also noted that it
       filters after reading, so it does not reduce IO.
     - **`crs` is a FALLBACK, not an override.** Both this class and `ShapeSource`
       said "override"; `_read_source` applies it only `if gdf.crs is None`. A file
       carrying its own CRS keeps it and the argument silently does nothing -- which
       is correct behaviour (it cannot reinterpret coordinates by accident) and the
       opposite of what the word "override" promises. This bit the user earlier in
       the session on a 4326 refinement layer passed `crs=2927`. The reprojection
       target is the CRS on the `GridSpec`, which is now said explicitly in both.
     - **Positional order is a trap, documented rather than changed.** These are
       dataclasses extending `DataSourceSpec`, so the real order is
       `(path, external, metadata, layer, query, crs, fields)`:
       `GeoPackageSourceSpec("x.gpkg", "areas")` sets `external="areas"`, truthy,
       and names no layer. Making the subclass fields keyword-only (`kw_only`) would
       fix it properly but changes the signature of every source class and their
       `_from_dict` callers; a `.. warning::` in each docstring was the proportionate
       call for a documentation pass. **Candidate deferral.**

174. **`GridSpec.voronoi` stops swallowing unknown keywords, and names its sizing
     options (2026-08-31).**
     Noticed that `breakline_buffer` worked but was not in the signature. It was
     reaching `**engine_options` and landing in a bucket nothing validates.
     Measured before changing anything: `breakline_buffer=40` gives corridors of
     239,869 ft2, the typo `breakline_bufer=40` gives 59,967 ft2 -- the default 10,
     silently, under a clean run and a "Normal termination". `compleletly_made_up=123`
     was accepted, stored in `options`, and never read. This is CLAUDE.md 8.8's
     bare-`**kwargs` rule (written for picture verbs) biting a module the rule had
     never been applied to.
     - **17 options are now named parameters**, grouped `boundary_*` / `refinement_*`
       / `breakline_*` plus `region_point_tolerance`. Static analysers read the `def`
       line and nothing else, so this is the only thing that makes them complete.
     - **DEVIATION FROM 8.8, deliberate: the named parameters default to `None`, not
       to the resolver's real default.** 8.8 says mirror the owning link's default.
       Mirroring here would mean always inserting a value, and a present-but-`None`
       `boundary_max_area` SHADOWS its own `default_cell_area`/`max_area` aliases in
       `_option`, which returns the first key PRESENT rather than the first non-None.
       So `None` means "not given" and the effective defaults live in the docstring's
       `Other Parameters` block instead. This also removes the default-drift risk
       that 8.8's mirroring requires a test to police.
     - **The allowlist is derived, and ratchets both ways.** `_VORONOI_OPTION_KEYS`
       admits the 5 legacy aliases and the two dict-read keys as well as the 17;
       `test_the_option_allowlist_matches_what_the_resolver_reads` re-greps
       `grid_spec_resolver.py` for every `_option(options, ...)` and `options.get(...)`
       and fails if the sets differ in EITHER direction -- a new resolver read the
       allowlist would reject, or an allowlist entry nothing reads (which would
       silently accept a typo again).
     - **Source shape is validated at the call too.** `refinement=[a, b]`, a single
       source in `breaklines=`, and a live `GeoDataFrame` all used to die inside the
       resolver as `AttributeError: 'X' object has no attribute 'path'`, naming
       neither the argument nor the fix. Each now raises naming both; the
       GeoDataFrame message explains that a `GridSpec` is a serializable recipe and
       gives the one-line `to_file` fix. Closes the deferral recorded in 172.
     - **Technically breaking, knowingly.** Any caller passing a keyword that was
       silently ignored now gets a `TypeError` at construction. That is a bug being
       surfaced rather than behaviour removed -- the option never did anything -- and
       the full suite (1968 tests) passes unchanged, so nothing in-repo relied on the
       silence.

175. **A NaN stack refuses to draw, instead of failing inside Matplotlib (2026-08-31).**
     `stack.plot.section(...)` on a stack with a NaN top or bottom died eight frames
     deep as `ValueError: Axis limits cannot be NaN or Inf`, from flopy's
     `PlotCrossSection._set_axes_limits` -- naming neither the layer at fault nor the
     cause. `LayerBuildResult.qc()` had detected the condition since it was written
     (`nan_top`, `nan_botm`, `nan_active_cells`); nothing on the picture path
     consulted it.
     - `_require_finite_geometry(verb)` now runs before the renderer, names each
       layer with its NaN count, gives the three real causes (a raster not covering
       the grid, contours interpolated inside a smaller hull, a nodata value read as
       elevation) and the three fixes (extend the source, `fill='propagate'`,
       idomain).
     - **Scoped to the section verb**, the one that failed. `vertex_grid()` is public
       and returns a flopy grid; raising there would break inspecting a
       partially-built stack, which is a legitimate thing to do while debugging
       exactly this. The helper is on `LayerBuildResult` so `map`/`surface` can adopt
       it if they turn out to fail the same way -- **not yet verified that they do**,
       and it was not worth building a NaN fixture for each on a report of one.
     - Worth carrying: the NaN **propagates downward**. A hole in one unit's bottom
       makes every `thickness=`-declared unit beneath it NaN too, so the report names
       several layers when one source is at fault. The measured fixture shows
       `sand` 3 cells and `clay` 3 cells from a single 3-cell hole.


169b. **The 8.8 no-bare-kwargs rule reaches the NOUN tier -- one method converted,
      17 recorded as debt (2026-08-30).**
      Reported as "`model.hds.map()` docstring doesn't show all arguments...
      `show_layer_elevs`, `show_mounding` aren't even mentioned". Measured: 29 of
      its 39 reachable parameters were invisible, and `show_mounding` does not
      decorate the map -- it REWRITES it, z going
      `153.76, 147.67, 144.27` to `22.00, 16.28, 15.79`.
      - **It is 18 methods, not one.** An AST sweep finds 18 classes defining a
        picture verb that takes a bare `**kwargs` and forwards to a plot function;
        13 have no `Parameters` section at all and 8 have a one-line docstring.
        Across the family that is 516 invisible parameter-slots and 110 named
        parameters an editor autocompletes and nothing explains.
      - **Why no test saw it.** Every enumeration in `test_plot_vocabulary.py` is a
        hand-written four-entry whitelist -- `SCOPES` is
        `{module, model, grid, stack}`, `_namespace` a four-entry dict. There is no
        discovery step, so the noun tier was unreachable by construction, in both
        directions. All 64 of its tests pass over this defect.
      - **It had been excused once already.** Ledger 145 recorded the leaf verbs as
        "assumed bare and measured otherwise". That check asked ONE question -- is
        the signature `(*args, **kwargs)`? -- where 8.8 states four requirements.
        Pointing the six tests at `DependentVariableFile.map` today, four fail.
        That is why the new test DISCOVERS its subjects instead of listing them.
      - ~~**The contract is signature-static, docstring-runtime**~~ -- **REVERSED
        2026-09-02, see 178.** This was written up as the contract; it was never
        a decision anyone made, only a description of what the code happened to
        do. The claim that "what a purely static reader gets is the SIGNATURE,
        which is the load-bearing assertion" is wrong: an editor is where these
        docstrings are actually read, and PyCharm showed ONE line of
        `model.plot.section` against 169 assembled at import. The splice is still
        the single author; its output is now written into the source.
      - **The tiers are MEASURED, not stylistic.** Every Tier 1 name changes the
        figure on both a results noun and a record noun. Tier 2
        (`kstpkper`/`per_timestep`/`bgs`/`hover_layers`/`hover_surfaces`/
        `show_layer_elevs`/`show_mounding`) changes it only on a per-layer field;
        on a record noun those produce a byte-identical figure, and `show_mounding`
        is worse than inert -- it injects a head-derived row into an elevation
        tooltip on every cell. Naming a parameter that provably cannot act is the
        same defect as hiding one that can. Four more (`hover_fields`, `zoom`,
        `rch_scale`, `animation_kstpkpers`) were measured inert on nouns and are
        recorded in `NOUN_INERT_PARAMS` rather than silently omitted.
      - **Five parameters now RAISE from the noun** rather than being accepted.
        `values=` was the reason: measured on three classes, it repainted every
        cell from the supplied array while the hover, title and colorscale went on
        reporting the noun's real field -- `hds.map(values=[1.0]*ncpl)` drew
        all-1.0 cells whose tooltips read 131.5 ft. That is a wrong picture, not a
        missing docstring, and the report would never have surfaced it. **This is a
        breaking change** for any notebook passing `values=`/`type=`/the legacy
        hover trio to a noun -- every such call was producing a mislabelled figure.
      - **Prerequisite fixes to the splice machinery**, all measured first:
        `_drop_parameter` walked the whole sectioned remainder, so dropping a
        parameter named `grid` or `section` also deleted the See Also entry sharing
        its name; it could not split a grouped head (`zmin, zmax :` was a no-op for
        both names); and `_inherit_verb_docs` was not idempotent -- a second call
        took `ModelPlots.map` from 149 lines to 297. `GridPlots.map` carried two
        `Parameters` headings, which is malformed NumPy: a reader stops at the
        first block.
      - **DEFERRED: the other 17**, listed exactly in
        `tests/test_noun_signatures.py::UNCONVERTED` and xfailed there. The list is
        exact in both directions, so a class fixed without being removed fails as
        loudly as a new one added -- it can only shrink. That shape is deliberate:
        ledger 145 is the precedent for what an unchecked "measured clean" note
        becomes. The remaining work is 17 signature rewrites (~30 named parameters
        each, from ~7) plus their `refuse_noun_parameters` call; the docstrings
        splice from `plot.map` and need no authorship.
      - Also deferred, found while measuring and NOT fixed here: `logscale=True`
        raises on heads (an object-dtype `zs` path, a numeric bug rather than a
        signature one); `contour_method='nearest'` is documented on `plot.map` and
        rejected by it; `contour_resolution` is inert under the default linear
        method. The last two are "a parameter you document must be one you accept"
        failing one level down, at the value set.


176. **`backend="mpl"` reaches the verb tier, and the noun tier is finished
      (2026-09-02).**
      Reported as "I'm having trouble understanding how I plot a cross section
      with the mpl backend -- `model.plot.section(backend='mpl')` rendered a
      plotly plot. I know I've created mpl cross sections before." Both halves
      were right, and the diagnosis was the interesting part.
      - **It was never a regression, and the check mattered.** `plot.section` was
        `def section(source, /, **kwargs)` before the 8.8 conversion and
        `XSection.__init__` parks unknown kwargs in `self._kwargs`, so
        `backend="mpl"` was swallowed there too -- identical behaviour before and
        after. Verified over every commit that has touched
        `src/myflopy/plot/__init__.py` since it was created: **0 occurrences of
        an mpl backend value at all nine**. The verb tier was born without the
        switch at 8.3a (2026-08-18), on the same day 8.2 gave the noun tier one.
      - **Three different answers to one typo**, measured before anything was
        changed: `model.plot.section(..., backend="mpl")` returned an `XSection`
        with `{'backend': 'mpl'}` stashed on it; `plot.map(..., backend="mpl")`
        returned a `Choro`, swallowed by `**trace_kwargs`;
        `vor.plot.section(line, backend="mpl")` raised `TypeError`. Meanwhile
        `model.hds.section(line=..., backend="mpl")` had worked all along.
      - **Why no test saw it.** `test_map_select.py` carries a test *titled*
        "Every grammar leaf offers `backend="mpl"`" whose body calls
        `picture.plot_mpl()` and never passes a backend -- the docstring asserts
        more than the test does. And `test_plot_vocabulary._chains()` anchors
        `"section": [XSection.__init__]`, the Plotly constructor;
        `render_xsections`, which is where `backend=` actually lives, is not in
        the chain, so a verb missing it was unmeasurable by construction.
      - **`backend` is now in `NOUN_MAP_PARAMS`, and that is the load-bearing
        decision.** It is the one name in that tier that is not a drawing option.
        It had always worked on the unconverted nouns *through their `**kwargs`
        tail*, so converting one to a named signature without it would have
        DELETED the matplotlib backend from that noun -- silently, with no test
        objecting. Listing it makes `test_noun_signatures` refuse the conversion
        instead. This is why the two jobs were done in one pass rather than in
        sequence: fifteen chances to quietly undo the thing being built.
      - **COMPROMISE -- `surface` gets no backend.** `InterpolatedSurface` has no
        `plot_mpl` and there is no matplotlib rendering of a 3-D height field
        here. Naming a parameter that cannot act is the same defect as hiding one
        that can (the rule 169b established), so it is documented as absent
        rather than added and made to raise. `animate`'s
        `backend={"plotly", "png"}` is left alone for a different reason: a
        rasteriser over frames is not a second renderer of one subject.
      - **`Choro.plot_mpl` ignores overlays, which makes HFB a special case.**
        A barrier is an edge drawn OVER the cells, so the ordinary
        `_apply_backend` route -- correct for every other noun -- returns a
        perfectly valid figure of the grid with every barrier missing: the one
        thing the picture is of. Both HFB verbs now draw the segments onto the
        axes themselves. Measured on the canonical wall: **32 barriers, 32
        `Line2D` artists**, on the inputs and results side both. That count is
        the assertion, because a type check passes straight over the defect.
        Two details worth carrying: `plot_mpl` renders in MODEL coordinates
        (measured, x limits -105..2205 against grid bounds 0..2100 on EPSG:2927)
        so the segments go on unprojected, where the Plotly overlay needs
        EPSG:4326; and `PALETTE.highlight` is Plotly's `rgb(214,39,40)`, which
        matplotlib REJECTS -- `PALETTE.mpl_highlight` is the twin the palette
        already carried for exactly this.
      - **BREAKING (small): `PRTPathlineView.map(backend="mpl")` returns a bare
        `Figure`, not FloPy's `(fig, ax)`.** Every other `backend="mpl"` in the
        grammar returns one figure, and a lone verb handing back a tuple means a
        caller looping over several pictures gets a `TypeError` from the odd one
        out. The axes are still `figure.axes[0]`. One test updated; no other
        caller anywhere in the repo.
      - **Naming PRT's parameters emptied a tail three checks were reading.**
        `map` consulted `**base_kwargs` to refuse base-map options on the
        Matplotlib branch, to refuse them on an already-built base, and to
        forward them. Naming them (8.8) would have silently disabled all three.
        `_explicit_options` reconstructs the same information from the named
        parameters by comparing each against the default in the method's OWN
        signature, so the two cannot drift. This is the general hazard of the
        8.8 conversion and the one to look for next time, so it now has a test
        rather than a comment. Note what that test canNOT probe with: the PRT
        fixture's grid is TWO cells, and two points produce no contour traces at
        all, so `contours=True` leaves the trace count unmoved and reads as
        "never forwarded" when the option arrived correctly -- instrumenting
        `_explicit_options` shows `{'zmin', 'contours', 'contour_levels'}`
        reaching the base map, and the same call on the canonical grid goes 1
        trace to 15. `zmin` is the probe instead.
      - **`values=` survives on `HfbPackageExplorer.map` where every other noun
        refuses it**, and the difference is real rather than an oversight. Those
        nouns fix a field, so an override repaints the cells while the hover and
        colorbar go on describing the real one. This noun's subject is the
        BARRIERS; the cells are a backdrop, and choosing what the backdrop shows
        contradicts nothing -- exactly as on `vor.plot.map`, which it delegates
        to. The rule as written already permits this (it tests keyword-only
        names, and `values` is positional here), which is fortunate rather than
        designed; it is recorded so the next reader does not "fix" it.
      - **`GridPlots.map` and `StackPlots.section` were NOT nouns**, and treating
        them as such was a category error inherited from 169b. They are the verb
        tier wearing a noun's shape. A noun fixes what is drawn, so
        `values`/`type`/`custom_hover` contradict it; a scope fixes only the
        SUBJECT, and `vor.plot.map(values=node_ids)` is the whole point of a grid
        map while `custom_hover` is the only hover a bare grid HAS, the sectioned
        one needing a model. Applying the noun rules would have deleted three
        working parameters to satisfy a rule written about something else. They
        are now split out as `VERB_SCOPES`, still covered by the docstring test
        and by `test_plot_vocabulary`'s mirror/subset/defaults rules.
      - **`GroupLakConnections` splices its own docstring**, unlike the other
        fourteen. It is layer 9 and `myflopy.plot` is layer 7, so importing it
        from the splicer pulls in `myflopy.project` (layer 13) while
        `simulation.base` is still half-built and the import fails outright --
        found by hitting it, not by reading the layer map. It reaches down for
        the reference from the bottom of its own module instead, which is the
        mirror image of what `StackPlots` needs.
      - **Three normalizations became one, each after a second call site wanted
        it.** `normalize_backend` (the alias set was written out three times) and
        `as_mpl_figure` -- because the renderers disagree about what they return:
        `Choro.plot_mpl` a `Figure`, the `figs` cross-section helper a
        `(fig, ax)` tuple, FloPy's patch renderer an `Axes`. Measured, not
        assumed: `vor.plot.grid(backend="mpl")` returned an `Axes` on the first
        pass and the test caught it.
      - **The docstring contract is unchanged and deliberate: signature-static,
        docstring-runtime.** Each noun's SOURCE carries its own local parameters,
        its Returns/See Also, and 5-6 Examples written for that noun; the ~22
        shared entries are spliced from `plot.map` at import. Every parameter
        therefore appears twice at runtime -- in the signature an editor reads,
        and as a full entry in `help()` -- without twenty-two entries being
        hand-copied into fifteen files, which is the drift the splice exists to
        prevent. `GroupLakConnections.map` comes out at 210 lines and 28
        documented entries (38 before the review below removed a spliced entry
        for a parameter it does not have, and a duplicate one).
      - **The three items 169b deferred are now CLOSED, and two of them this
        pass had made worse.** Splicing `plot.map`'s reference onto eighteen
        nouns propagates its errors as well as its text, so a single wrong entry
        became eighteen wrong entries -- which is the cost side of the splice and
        is worth stating plainly.
        - `logscale=True` raised on heads:
          `TypeError: loop of ufunc does not support argument 0 of type float
          which has no callable log10 method`. Two branches of `Choro.zs` offer
          it; the `custom_zs` one routed through `_logscaled` (which coerces to
          float and blanks non-positive values to NaN) while the head branch
          called `np.log10` raw. So every RECORD noun honoured `logscale=` and
          the one field people actually log-scale did not. Fixed by using the
          helper that was written for this and already sat beside it. Guarded by
          comparing the log render against the map's OWN linear render -- not
          against `hds.array(...)`, which resolves its default output time
          separately and differs in one cell.
        - `contour_method='nearest'` was documented and has never been
          implemented; the accepted set is `linear`/`tri`/`tricontour` and
          `cubic`/`clough`/`clough_tocher`/`cloughtocher`. Documentation
          corrected to what the code accepts, aliases named.
        - `contour_resolution` is not inert in general -- it is inert under the
          DEFAULT method. `_linear_contour_segments` triangulates the cell
          centres and does not take a resolution at all; only
          `_cubic_contour_segments` interpolates onto a resolution-square grid.
          Measured 338 contour points at 40, 150 and 400 alike under `linear`.
          The entry now says which method it applies to.
        - **A new guard for the whole class**: `test_plot_vocabulary` already
          asserted that every documented NAME is accepted. It now asserts the
          same of the VALUES inside a documented `{...}` set, which is where this
          defect hid one level down. Mutation-checked -- putting `'nearest'` back
          in the docstring fails naming the value and the error it raises.
      - **An Examples block is prose to every other check, and two of them
        lied.** `LakConnectionsExplorer` was documented at
        `model.packages.lak.inputs.connections` -- the explorer hangs off the
        PACKAGE, not off `.inputs` -- and `SurfaceWaterExchangeResultsExplorer`
        at `model.surface_water.results.exchange`, where it is
        `model.packages.surface_water.results.q`. Both plausible, both wrong,
        and invisible to the signature tests and to a
        documents-what-it-names check alike. A new parametrized test walks the
        ATTRIBUTE CHAIN of every `>>>` line against the canonical model; it fails
        naming the noun and the path. Not the call, only the chain -- drawing
        eighteen maps in every documented combination would make this the slowest
        file in the suite for much less. Mutation-checked: retyping one example
        as `model.packages.lake...` fails by name.
      - **A noun's own hover spec is the BASE, not an override -- caught by a
        test that had been adapted rather than fixed.** The conversion first
        passed each noun's default through `hover=`, the call-site override slot,
        because that is the name the verb exposes. It renders identically, so
        the only symptom was `choro.hover_spec` coming back `None` where
        `model.hds.map().hover_spec` carries a spec -- and
        `test_model_budget_namespace` was edited to read
        `_resolved_hover_spec()` instead, which passes over the change rather
        than questioning it. Reverted at all 11 sites: the default now goes to
        `hover_spec=` (the base, which is where it went before 8.8) and `hover=`
        stays the caller's. The adapted test is back to its original assertion.
        The lesson is the one worth carrying: a green suite after a test edit is
        not evidence, and "the assertion moved" is the signal to look at.
      - **The debt list is now EMPTY, and that is the point of the pass.**
        `UNCONVERTED` is kept rather than deleted, and an assertion says so, so
        that the three parametrized tests go from "xfail on the debt list" to
        "assert on every noun in `src/`". Anything added back is debt with a
        ledger entry, not an exemption. 18 discovered methods: 16 nouns and the
        2 verb scopes, 51 xfails to 0.
      - **Also fixed in passing**: the ledger had TWO entries numbered 169. The
        noun-tier one is now 169b, on the user's instruction; `usg.export_gis`
        keeps 169, which is what CLAUDE.md references.
      - **`hover_spec=` was the pre-8.8 spelling of `hover=`, and naming `hover`
        silently took it away.** Every noun used to end in
        `kwargs.setdefault("hover_spec", <default>)`, so a caller passing
        `hover_spec=` kept it. Once `hover=` is a named parameter the noun always
        supplies one, and `Choro` stores it as `_hover_override`, which
        `_resolved_hover_spec` prefers -- so the caller's `hover_spec=` still
        reached the constructor and lost, with no error. Measured on the
        canonical model, HEAD vs converted: `drn.results.q`, `budget.sto_ss`,
        `lak.results.stage`, `sfr.results.stage` all went honoured -> dropped.
        `refuse_noun_parameters` could not cover it -- `hover_spec` is not a
        `plot.map` parameter, and `test_the_tiers_partition_the_free_verb`
        forbids naming one that is not -- so the resolution is a sibling helper,
        `resolve_noun_hover(hover, trace_kwargs)`: it pops the legacy key off the
        tail and returns it when `hover=` is None, `hover=` winning when both are
        given. **Judgment call: `hover_spec` stays UNNAMED.** Naming it would put
        two spellings of one argument in an editor's autocomplete, and the rule
        is that a narrower scope offers fewer parameters, never other ones.
        Applied in `package_results.py` (both methods) and, in the adversarial
        review below, `project/group/lak.py`. **The same defect is still live in
        `package_surface_water.py` (5 sites) and `package_inputs.py` (3 sites)**
        -- each needs the one-line `hover = resolve_noun_hover(hover,
        trace_kwargs)` after its `refuse_noun_parameters` call;
        `test_package_map_hover_spec_is_overridable` fails until they get it
        (it routes through `package_surface_water.py`).
      - **A test that read `Choro.hover_spec` had to move to
        `_resolved_hover_spec()`.** `test_terms_with_no_package_accessor_are_reachable`
        asserted `budget.sto_ss.map().hover_spec.title == "STO-SS q"`; the
        attribute is now `None` for every noun, because the default arrives as
        `hover=`. The attribute is one of two inputs, the method is the answer,
        and the question the line asks is "which spec will this map draw with".

  - **Adversarial review of `GroupLakConnections.map` (2026-09-02).** Three
    defects, all introduced by the conversion, all found by running code rather
    than reading it. The mechanical part of the conversion was clean: all 22
    `NOUN_MAP_PARAMS` keyword-only, none of `NOUN_REFUSED_PARAMS`/
    `NOUN_INERT_PARAMS` offered, every one of the 21 forwardable names actually
    forwarded, every default equal to `ModelPlots.map`'s, nothing dropped or
    renamed against HEAD (only `colorscale`'s annotation widened). The defects
    were all in the parts a signature diff does not cover:
      - **The docstring splice documented a parameter the method does not
        have.** `inherit_map_docs(..., keep=set(NOUN_MAP_PARAMS) | {"per",
        "layer"})` inherited `plot.map`'s `per` entry, but this noun has no
        `per` -- connection geometry is static, so the forward hardcodes
        `per=0`. `map(per=3)` therefore reached the tail and died as
        `ModelPlots.map() got multiple values for keyword argument 'per'` --
        naming a class the caller never typed, which is verbatim the defect
        `refuse_noun_parameters` was written against, reintroduced through the
        docs. `keep` is now exactly `NOUN_MAP_PARAMS`; `layer` came out with it
        because the method documents `layer` locally and inheriting it too put
        **two `layer` entries in one Parameters block**. 38 entries -> 28, no
        duplicates, no lies. **The three signature tests cannot catch either
        one**: they check named -> documented, never documented -> named, and
        the splice runs at import so the source they parse looks fine.
      - **`per=` is now refused explicitly**, since every sibling noun takes one
        and reaching for it here is the natural mistake. A local `TypeError`
        naming the noun, not a collision inside a class the caller never named.
      - **The `hover_spec=` regression, found independently and then handed to
        the shared fix.** Measured on the canonical model:
        `map(hover_spec=cell_input_hover("belev"))` set `Choro.hover_spec` and
        left the figure **byte-identical** to `map()`. First fixed here by
        raising; that was **wrong and was reverted** -- `test_hover_integration.py::
        test_package_map_hover_spec_is_overridable` is a HEAD test pinning
        `hover_spec=` as a *working override*, so refusing it breaks a contract
        rather than restoring one. Now `resolve_noun_hover(hover, trace_kwargs)`,
        the shared helper, with `hover=` winning when both are given.
    Verified by running: `backend="mpl"`/`"matplotlib"` return a
    `matplotlib.figure.Figure` and `backend="nope"` raises `ValueError`;
    all 6 docstring Examples execute; each of the 23 named selectors and drawing
    parameters changes `fig.to_dict()`. Two apparent inertias were **the test
    model, not the code**: the canonical model has ONE lake, so `lake=0` is the
    default (`lake=1` selects 0 of 72 rows and does differ), and the group has
    one member, so `model_name="base"` is the default. `contour_method` acts
    (`"linear"` vs `"cubic"`); `"nearest"` is documented on `plot.map` and
    rejected by it, which is the deferral already recorded above.
    **Caveat on the parallel sweep**: the shared tree was being edited by other
    agents throughout, so suite results moved between runs -- a transient
    self-import in `prt_maps.py` broke `import myflopy` outright at one point,
    and `test_deferred_import_ratchet` and
    `test_terms_with_no_package_accessor_are_reachable` each failed and then
    passed with no change to this file. Attribution here was done by
    `git stash push` on `project/group/lak.py` alone and re-running.

177. **`layer=` on the grid's feature-to-cell selectors (2026-09-03).**
    `vor.get_vor_cells_as_dict` and `vor.get_vor_cells_as_series` both reached a
    vector file through a bare `gpd.read_file(path)`, which reads the **first**
    layer and, for a multi-layer GeoPackage, emits only a `UserWarning`. So there
    was no way to ask for any other layer, and reading the wrong table was silent
    in any caller that does not surface warnings. Both now take
    `layer: str | int | None = None`, threaded through a shared `_read_locs`
    helper in `grid/selection.py`; the wrapper methods on `VoronoiGridPlus` mirror
    it. Reported from a real pond footprint stored beside other layers in one
    `.gpkg`.
    - **`layer=None` is GDAL's own default**, so every call written before the
      parameter is byte-identical -- pinned by
      `test_vector_layer_selection.py::test_default_is_unchanged`, which asserts
      the no-argument call equals `layer=<first layer name>` rather than
      restating a literal.
    - **Validated up front rather than by catching.** GDAL's message is
      `Layer 'x' could not be opened`, which does not say what IS in the file --
      exactly the question the caller has. `_layer_names` uses
      `geopandas.list_layers` to raise a `ValueError` naming the available
      layers. That function is geopandas **1.0+** while `pyproject.toml` declares
      `geopandas>=0.14`, so it is fetched with `getattr` and the check is
      **skipped** on an older geopandas, which then gets GDAL's poorer message.
      Deliberate: raising the floor for one error string is not worth it, and
      validating up front avoids a broad `except` that the narrowness rule would
      otherwise force onto the allowlist (pyogrio's `DataLayerError` subclasses
      `RuntimeError`, fiona's `DriverError` does not, so no narrow tuple covers
      both engines).
    - **`get_vor_cells_as_series` REFUSES `layer=` for an in-memory geometry**,
      where there is no file to pick a layer from -- a parameter that cannot act
      is the same defect as a missing one (the rule `surface` follows for
      `backend=`).
    - **Not extended to `read_gpkg` / `read_shp_gpkg`** in
      `modflow/utils/datatypes/readers.py`, which have the same first-layer-only
      behaviour. `read_gpkg`'s docstring claims it is "useful when you have
      multiple layers", which its `gpd.read_file(gpkg_path)` does not deliver;
      that is a separate defect, on a reader with more callers, and was out of
      scope here.
    - **The nested return shape was left alone.** `get_vor_cells_as_dict` returns
      `{name: [[cell, ...]]}` -- one inner list per matched feature, not a flat
      cell list. `utils/datatypes/locs.py` and `mf6/heads_observations.py` both
      depend on it (the latter's `assert len(cell_ids) == 1`, whose message says
      "more than one cell found", is in fact counting groups). Changing it is a
      breaking change worth doing separately; the new test flattens twice and
      says why.


177. **`section(fill=...)`: the cross-section of the MODEL, not of a line
      (2026-09-02).**
      Reported as "so we really have no way to easily draw an mpl cross section
      showing the grid, layers, and results? Please double check. It's not on
      `model.plot.grid`?" -- after the same user had gone looking for it once
      already the same day. Audited every route before building anything:
      | route | drew |
      |---|---|
      | `model.plot.grid(backend="mpl")` | a PLAN view. Not a section at all |
      | `model.plot.section(line, backend="mpl")` | 2 lines, 0 filled collections |
      | `model.hds.section(line, backend="mpl")` | the same 2 lines |
      | `vor.plot.section(line, backend="mpl")` | 124 lines -- cell OUTLINES, no results |
      | `plot_model_cross_section(m, line, kstpkper=)` | 2 collections + 31 lines |
      So the picture existed, in exactly one place, reachable only by deep import
      from `myflopy.modflow.mf6` and bound to no scope. That is why it was lost
      twice.
      - **One verb, not a sixth.** Geometry chooses the verb and both pictures
        are vertical slices, so `fill=` selects the subject within `section`
        rather than earning a `filled_section`. The two answer different
        questions -- *what is the head along this line* against *what does the
        model look like through here* -- but they answer them about the same
        geometry.
      - **`fill=` is Matplotlib ONLY, and raises rather than degrading.** FloPy's
        `PlotCrossSection` is the renderer and there is no Plotly counterpart to
        defer to. `fill=` with the default backend raises naming `backend="mpl"`,
        which is the same shape `plot.grid` already uses for `pathlines=` with
        the plotly backend. Silently returning the line profile instead would
        have been this session's own defect, one verb along.
      - **The profile-only arguments RAISE too.** `interpolate`, `spacing`,
        `num_points`, `use_rbf`, `interpolator`, `x_or_y`, `show_model_top`,
        `show_model_btm`, `surf_type`, `clip`,
        `extrapolate_beyond_section_ends`, `animation_kstpkpers` shape a sampled
        line; the filled branch walks the grid cell by cell and has nothing to
        sample. `_PROFILE_ONLY_SECTION_ARGS` names them so the message can list
        what was given, instead of turning a knob that does nothing.
      - **A COMPROMISE caught by its own test: the grid branch silently
        downgraded `fill="results"` to `"layer"`.** First draft coerced the
        unanswerable value because a bare grid has no results -- handing back a
        picture that looks right and is not the one asked for. It raises now,
        naming what a grid does not have. Written into the test before the fix,
        which is the only reason it did not ship.
      - **A flat `(ncpl,)` array is refused.** It would colour every layer the
        same and read as a real answer; the message says to stack the layers or
        use `fill="results"`. `model.hds.array()` returns exactly that flat
        shape, so this is the mistake the API invites.
      - **`fill="results"` is kind-neutral**, resolved through the model's own
        `_field_reader` rather than `.hds` -- heads on GWF, concentration on GWT,
        temperature on GWE, the same rule ledger 99 established for `XSection`.
      - **`cells=` keeps working**, converted to the polyline through those
        cells' centroids, rather than being refused for a reason a caller would
        find arbitrary.
      - **The layer fill still draws the RESULTS.** `fill="layer"` is not the
        geology alone: the simulated head goes on as a water surface over it,
        because that is the picture people mean. Measured: 2 filled collections
        + 31 head-surface artists on the canonical model. `fill="results"`
        replaces the layer colouring instead of adding to it -- a cell has one
        fill, and drawing both would leave the legend and the colorbar
        describing the same patch differently -- so its layer legend is
        suppressed and it earns a colorbar (the observable difference the test
        asserts on).
      - **Fixed in passing:** `docs/build_myflopy_api_pamphlet.py` documented
        `mf.plot_model_cross_section(model, line=section.line)`, which raises
        `AttributeError` under the documented `import myflopy as mf` -- it is
        exported from `myflopy.modflow.mf6`, never top-level. It now shows both
        section spellings on the grammar instead. `model.xsect(...)` in the same
        block was also long retired.
      - **Three follow-ups from the first look at it, same day.** "The labelling
        of the units didn't show up in the legend. And the `layer` arg doesn't
        seem to do anything for the mpl pic. I'd like to be able to define what
        layers to draw. And really I should be able to decide what water levels
        to draw."
        - **`layer=` was accepted and IGNORED** -- the defect this family of work
          exists to end, reintroduced by the very change that closed it
          elsewhere. `_PROFILE_ONLY_SECTION_ARGS` listed twelve names and not
          this one, because `layer` is real on the profile branch. The filled
          branch takes `layers=` (plural) and `layer=` now raises naming it: two
          spellings for two meanings beats one that silently means neither.
        - **`layers=` also had to stop the OUTLINES.** Reported after the first
          cut: "only the layers I list are filled, but cells for all layers
          still draw". FloPy's `plot_grid` outlines every cell and takes no
          layer filter, so masking the fill left the excluded layer drawn.
          Cropping cannot cover for it -- a layer's elevation range OVERLAPS its
          neighbours', and the canonical bottom layer spans -1.3..59.8 straight
          through a retained band of 28.5..166.0. With a subset the edges now
          come from the filled collection itself, which is already masked; with
          every layer the two-collection picture is untouched.
        - **`legend=` places it, with words that exist.** `"bottom"` is the
          obvious one and is not a matplotlib `loc` value, so a friendly set
          (`top`/`bottom`/`left`/`right`, the four corners, `auto`) maps onto
          what `loc` accepts and matplotlib's own strings pass through
          unchanged. `"outside <side>"` is not a matplotlib placement at all and
          is the one that matters here: a layer legend anywhere INSIDE the axes
          sits on top of the geology it describes. An unknown word RAISES and
          lists the accepted ones -- falling back to `"best"` would be a legend
          that ignored where it was told to go.
        - **`show_grid=False` was silently dropped**, found by the probe written
          to check the above. The filled branch builds no class that takes an
          open tail -- that is the profile branch's `XSection` -- so anything
          left in `**kwargs` was accepted and discarded. It is a named parameter
          now, and leftovers RAISE.
        - **`layers=` masks AND crops.** FloPy's `plot_array` skips NaN cells, so
          masking hides them without moving any other cell's geometry -- but the
          axis stays full height, which puts the two layers you asked for in a
          thin band. Measured: `layers=[0, 1]` goes 124 drawn cells to 62. The
          legend narrows to match.
        - **The crop must come from the DRAWING, not from the grid.** Reported as
          "the axes need to adjust to the new extent; the y-axis scale still
          assumes all layers are drawn". The first cut used the modelgrid's own
          `top`/`botm`, which are whole-GRID statistics -- so the limits came
          from cells the section line never crosses. Measured, `layers=[0]` gave
          an axis of 66.3..164.2 for content spanning 72.0..146.2: a quarter of
          the height empty. Now measured off the drawn polygons, which is exact
          and needs no grid introspection. Water surfaces are `Line2D` rather
          than collections and had to be added explicitly -- without them
          `layers=[3]` with the default `head_layers=0` cropped out the very line
          its own legend advertised. Full-layer sections keep FloPy's framing:
          it is what `plot_model_cross_section` and everything built on it
          already use, and tightening it unasked would move every one.
        - **`head_layers=` chooses whose water levels are drawn**, defaulting to
          `0` -- the single water table every earlier figure got, so nothing
          moves. A list draws one surface per layer, each LABELLED, because
          several unlabelled lines on one section cannot be told apart; colours
          come from `category_colors`, keyed by label, so a layer's water level
          is the same colour on every section it appears in. `None` draws the
          geology alone. Ignored under `fill="results"`, which already paints
          that quantity onto the cells.
        - **Layer names come from the SPECS, which was the user's own
          suggestion.** `ModelSpec.build` stores the `ModelContext` on the model
          (`specs.py:2141`), so a model declared with
          `ModelContext(surfaces=stack.build(vor))` carries `.names` -- "sand",
          "clay" -- and the legend uses them with no argument at all.
          `layer_labels=` overrides. The fallback matters and is tested: the
          canonical model is built imperatively, has `myflopy_context is None`,
          and correctly gets `Layer N`. A `surfaces` that is a plain frame (the
          shape `ModelContext`'s own docstring shows) also has no names.
      - **NOT done, and deliberately:** `plot_model_cross_section` and
        `plot_layered_cross_section` stay where they are and keep their callers.
        They are the ENGINE under `fill=` (`LayerStack.plot.section` uses the
        same renderer), so deleting them would mean inlining a renderer into a
        verb. What changed is which spelling the docs teach.


178. **Spliced docstrings written into the source, because an editor cannot run
      your module (2026-09-02).**
      Reported as "`model.plot.section` docstring doesn't show in PyCharm. It
      shows the args but the docstring itself is only one line. Needs complete
      docstring. This is a reoccurring issue" -- and then, decisively, "I don't
      remember deciding the docstrings were runtime only."
      - **The user was right that they never decided it.** Ledger 169b recorded
        "signature-static, docstring-runtime" as the contract earlier the same
        day, reasoning that the SIGNATURE is what a static reader needs. That is
        the same argument plan 8.8 demolished for `**kwargs`, one field along:
        `__doc__` assigned at import reaches `help()` and reaches no editor.
        Measured across the bound tier: `ModelPlots.section` 1 source line
        against 169 runtime, `.map` 1 against 188, `.surface` 1 against 51,
        `.grid` 3 against 82 -- and every noun the same shape.
      - **Generate, then pin** -- the contract this repo already runs on for
        `api_snapshot.json` and `import_layers.json`.
        `scripts/derive_docstrings.py` writes each assembled docstring back into
        its own file; `tests/test_docstrings_are_static.py` fails in BOTH
        directions, so editing `plot.map`'s reference without regenerating fails,
        and hand-editing a generated docstring fails too. 29 methods, discovered
        by the `_myflopy_own_doc` marker the splicers already stamp rather than
        listed, so a noun added to the splice cannot be left behind.
      - **The prerequisite was idempotency, and the splice did not have it.**
        Writing the text into the source means the splicers re-run over their own
        output at every import; a non-fixed-point would grow the docstring each
        time. Measured: a second pass took `ModelPlots.section` from 169 lines to
        179. Three duplications, each its own small bug:
        - the free verb's extended PROSE was appended unconditionally;
        - the `Bound form of ...` footer likewise;
        - and `**kwargs` was duplicated because `_parameter_name` requires a
          `name : type` head, while a VAR-ARG entry legally has no type -- so the
          name filter could not see it to drop the inherited copy.
      - **A fourth needed care: an entry head with NO type at all.**
        `viz.mosaic` documents every parameter as a bare `panels` with the
        description indented under it -- legal NumPy, and equally invisible to
        the filter. The obvious fix (treat a bare identifier at column 0 as a
        head) is WRONG and was caught immediately: `Returns` is a valid
        identifier, so section headings started parsing as parameters and all
        thirteen methods drifted. It is only unambiguous INSIDE a `Parameters`
        block, so `_bare_parameter_name` is scoped there and `_documented_names`
        walks the blocks rather than the whole remainder.
      - **`mosaic` needed two generation passes, and that is not drift.** Its
        first pass applies the canonical NumPy section ORDER (177 put `Returns`
        before `Examples`), which changes the text once; the second is a fixed
        point. Verified to three passes: 66, 68, 68, 68.
      - **Found in passing: `GridPlots.map` documented four of its own
        parameters with bare heads**, so once the filter could see them it
        dropped the typed inherited copies and left entries the noun docstring
        test rejects. They now carry real `name : type` heads -- which they
        should have had regardless.
      - **NOT changed: the splicers still run at import.** They remain the single
        author, so the reference is still written once on the free verb. What
        changed is that their output is committed rather than assembled only in
        memory.
