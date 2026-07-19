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
> `_ArealBuilder`/`EVTBuilder` (entry 6), 4.7.1 fixes `rch_spec` (13) and the
> budget/mover lists (14). They stay listed until the code actually lands.

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

14. **`_MOVER_PACKAGES` still lists `rch`, which has no mover support.**
    - What: added `riv` (FloPy confirms `mover=True`); left the pre-existing
      `rch` entry, which FloPy shows has no `mover` option.
    - Why: harmless — the budget-term lookup is `try`-wrapped and finds
      nothing — and removing it is an out-of-scope behavior change.
    - Revisit: whenever the mover diff is next touched.

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
    - What: head targets sample the regional water table (the surface that also
      seeds initial conditions) plus a small deterministic per-well offset,
      rather than carrying real measurements or `NaN`.
    - Why: the model is synthetic — there are no measurements. Before this,
      every head target was `NaN`, which made `compare()`, `stats()`,
      `residuals()` and `calibration_plot()` all vacuous, and produced the
      reported blank calibration plot. Sampling the *conceptual* surface rather
      than the simulated heads keeps residuals real; copying simulated heads
      would make the calibration plot a tautology.
    - Impact: the demo scatter shows a **systematic bias** (RMSE ≈ 6.6 ft on a
      ~33 ft head range, ~20%). That is honest — the regional surface is an
      approximation, and the model departs from it near pumping, the lake and
      the streams — but it is not a *well-calibrated* looking demo.
    - Revisit: **superseded by entry 25** — the user has asked for a
      well-calibrated demo, so this bias is a shortfall to close, not a
      compromise to keep.

22. **Pumping wells are observed only at period 0.**
    - What: `shallow_pumping` / `deep_pumping` carry a measured value in the
      first stress period and `NaN` afterwards.
    - Why: a regional-survey value is not a valid measurement inside a well that
      is actively drawing down (the deep well swings ~45 ft once pumping ramps).
      Carrying it forward would plant one meaningless outlier that dominates
      every residual statistic.
    - Impact: 14 paired points instead of 24; the pumping wells contribute one
      pre-development baseline each.
    - Revisit: if time-varying synthetic observations are ever generated from a
      truth run (that is what `canonical_calibration.py` already does for PEST).

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

25. **OUTSTANDING (user requirement, 2026-07-18): the fast tour must show a
    WELL-CALIBRATED model.**
    - What: `canonical_fast_tour.ipynb` currently renders a calibration scatter
      with a visible systematic bias — RMSE ≈ 6.6 ft on a ~33 ft head range
      (~20%), with `regional_center` biased ≈ −9 ft and `pond_mound` ≈ +3.7 ft
      in every period. The user's requirement is that the flagship notebook
      show a *well-calibrated* model, not merely a non-empty one.
    - Why it is not yet met: the "measured" heads are sampled from the regional
      water table, which is the **conceptual** surface used to seed initial
      conditions, not the model's converged solution. The two differ by several
      feet wherever the stresses bite. Fixing the plot (entry 23) and populating
      the targets (entry 21) made the calibration tier *work*; neither was aimed
      at making it look *good*.
    - Constraint that shapes the fix: the observations must not simply be copies
      of the simulated heads, or the scatter becomes a tautology that teaches
      nothing (see entry 21).
    - Candidate approaches, in rough order of preference:
      1. Sample observations from the regional surface but **calibrate the
         surface itself** — fit the `146.0 - 60.0*xn + 0.45*cross_relief`
         coefficients (and the 0.6 ft/layer decline) so the conceptual table
         tracks the converged solution. Keeps observations independent of the
         run while shrinking the bias; also improves initial conditions.
      2. **Widen the observation network.** Only 4 head targets exist and 2 are
         pumping wells — a thin, stress-dominated network by construction. Add
         valley-floor monitoring wells away from the pumping/lake/stream
         influence, where the regional surface is already a good approximation.
      3. Derive observations from a **truth run plus small noise**, the pattern
         `canonical_calibration.py::_truth_head_targets` already uses for the
         PEST notebooks. Gives an excellent-looking scatter, but makes the
         canonical targets depend on having run, and residuals become pure
         injected noise.
    - Revisit: next session. Entry 21 stays until this lands.
