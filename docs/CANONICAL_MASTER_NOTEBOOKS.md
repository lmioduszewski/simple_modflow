# Canonical Master Notebook Set

These notebooks are the authoritative end-to-end examples for `myflopy`, all
built on the same four-layer refined DISV/Voronoi canonical valley model.
Notebooks 00–03 build and explore it with `mf.build_canonical_model()`; the PEST
notebooks 04–06 calibrate that same model via `build_canonical_calibration_demo`
(it builds the canonical model as the synthetic truth, then perturbs K).

`canonical_model_template.ipynb` is the editable scaffold for creating or
reworking that canonical model from real GeoPackage inputs. It shows the
preferred builder API directly rather than calling `mf.build_canonical_model()`.

The set shares `canonical_notebook_style.py` for a consistent publication-style
header and visual language. Each notebook combines explanatory Markdown,
annotated executable code, interpretation guidance, and a clear handoff to the
next topic.

0. `canonical_fast_tour.ipynb`: the whole surface in about a minute — builds
   the smallest contract-complete profile (`CanonicalModelConfig.testing()`,
   the same one the test suite shares), optionally fires the entire myflopy
   test suite from a cell, then tours heads/maps/animations, package results,
   observation scatterplots, and a real 4-worker PESTPP-IES ensemble.
1. `canonical_00_model_and_build.ipynb`: conceptual model, mesh, build, and run.
2. `canonical_01_packages_and_observations.ipynb`: packages, named regions, and targets.
3. `canonical_02_visual_diagnostics.ipynb`: pronounced visual signals and standalone exports.
4. `canonical_03_prt_and_parallel.ipynb`: MF6 PRT and two-to-eight-way model splitting.
5. `canonical_04_pest_and_results.ipynb`: calibrate the canonical model with the
   declarative `PestProject` facade + PESTPP-IES, then review the fit.
6. `canonical_05_pest_calibration_setup.ipynb`: the modern PEST setup in detail —
   `parameterize` / `observe` / `forecast` / `build` over the canonical model.
7. `canonical_06_pest_ies_uncertainty.ipynb`: ensemble uncertainty + the per-cycle
   IES diagnostics (prior Monte Carlo, prior-data conflict, phi distribution &
   contributions, parameters at bounds, forecast uncertainty) and the spatial
   review layer — parameter-field maps, the uncertainty-reduction map, the
   prior-vs-posterior mosaic and the observation residual map -- then takes the
   calibrated base realization home (`ies.calibrated()`, ledger 560).
8. `canonical_07_transport_and_prt.ipynb`: solute transport on the canonical model: the
   GWT sibling, the concentration reader and transport budget, PRT particle tracking
   from the contaminant source, and a PEST calibration whose observations are
   **concentrations** (it recovers K from plume shape and timing, which heads alone
   cannot pin down).
9. `canonical_08_recalibration_rehearsal.ipynb`: the dress rehearsal of the Ten Trails
   recalibration and the acceptance test of the PEST preparation (ledger
   589). Not the four-layer canonical model: the canonical valley's 21x21
   grid carrying a model of Ten Trails' SHAPE built by
   `myflopy.modflow.mf6.canonical_rehearsal.build_rehearsal` -- five layers, NEWTON, an
   array RCH on an `mf.AverageYear` cycle, two RIV packages (`riv`, `crisp`), layered
   storage, two excavated `mf.Basin` lakes under one IMPLICIT LAK -- as a START and a
   TRUTH in one `mf.Project`, the truth DECLARED through the package-first API. The
   workflow, written once in the API to copy for Ten Trails: dated hydrographs binned in
   average-year mode, the equilibrium check, zoned pilot-point K + month-group recharge
   + Sy + crisp's conductance + per-lake bed leakance, heads + differences (each reading
   counted once: the levels carry a datum error, the differences the rest of the
   reading's weight) + lake stages
   + crisp's baseflow + a one-sided bound, **the truth round trip before any IES**
   (`run_template` at `truth_parameter_values` reproduces the noise-free data), PESTPP-IES,
   recovery against the truth, and the write-back (`ies.calibrated().to_library()`)
   rebuilt, run and checked against what IES ran. It ends in a `rehearsal_report`.
   **It is EXECUTED**, unlike 04-06: `tests/test_recalibration_rehearsal.py` runs it
   to its round trip and end to end and asserts the report's per-run gates (the knobs in
   its first cell set the run time; `docs/reproducers/rehearsal_recovery_budget.py` sweeps
   them). Whether the posterior covers the truth is judged POOLED over six converged runs
   (`canonical_rehearsal.pooled_recovery`, ledger 593), by the reproducer and
   an opt-in test (`SIMPLE_MODFLOW_REHEARSAL_COVERAGE=1`).

### Reference tours (not part of the numbered series)

- `plotting_vocabulary_tour.ipynb` — the whole drawing API in one place: the six
  verbs, the four scopes that carry them, the output contract every picture
  answers, and the two `backend=` switches. Runs on the smallest canonical
  profile, so it is fast. Read `docs/view_layer_conventions.md` alongside it —
  this notebook is the runnable version of that file.

Use `CanonicalModelConfig.validation()` for quick work and
`CanonicalModelConfig()` for the full model with at least 10,000 horizontal cells.
The PEST notebooks run PESTPP-IES, which fires hundreds of forward solves — use
parallel `workers=N` and a modest ensemble (they are "go get coffee" cells).
