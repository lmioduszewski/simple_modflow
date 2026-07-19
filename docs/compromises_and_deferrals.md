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

> Entries 6, 13 and 14 below are now **scheduled** rather than open-ended: plan
> §4.7 (package declaration consolidation) absorbs them — 4.7.5 builds the
> `_ArealBuilder`/`EVTBuilder` (entry 6), 4.7.1 fixed `rch_spec` (13) and the
> budget lists, and §4.7.2 carries the mover list (14). They stay listed until
> the code actually lands.
>
> **Correction (2026-07-18):** this banner previously said 4.7.1 would fix
> entry 14 (`rch` in `_MOVER_PACKAGES`). 4.7.1 shipped and did NOT — entry 14
> itself always recorded that as a deliberate out-of-scope deferral, so the
> banner overstated the promise, not the entry. Moved to 4.7.2, where `mover`
> becomes a declared descriptor capability and the fix lands with its permanent
> home. Flagged by the 2026-07-18 next-step survey.

6. **No `EVTBuilder` (no file-less whole-domain EVT form).**
   - What: `mf.evt` has `()` / `.gpkg` / `.flopy` but no builder form like
     `mf.rch(context=, nper=, recharge=)` that self-selects top-active cells.
   - Why: `RCHBuilder` broadcasts one value per cell; EVT needs a per-cell
     *elevation* (`surface`), which is new design surface, not reuse.
     `mf.evt.gpkg(..., surface=CellSurfaceOffset("cell_top"))` covers the
     areal land-surface case through the proven mapping engine.
   - Impact: whole-domain ET requires a domain-covering polygon or direct
     `stress_period_data`.
   - Revisit: on real user demand (recorded in plan §5.2 banner).

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

18. **Only SFR's derived table was converted; the other loose plot verbs stay.**
    - What: `sfr.results.q.plot_profile`, `sfr.results.stage.plot_profile`,
      `lak.…plot_budget`, and the bare `plot` on the field explorers keep their
      current spellings. Only `long_profile`/`plot_long_profile` moved onto the
      view shape.
    - Why: those are *field*-level verbs already inside the documented spatial
      grammar, and they were not what the user reported. Converting them is a
      mechanical follow-on best done as one sweep with its own snapshot diff.
    - Impact: the convention doc is normative for NEW nouns; existing
      field-level `plot_*` verbs are not yet uniform with it.
    - Revisit: as a view-layer pass alongside plan §4.7's package consolidation.

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
