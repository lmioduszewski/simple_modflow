# Archived PEST notebooks

These notebooks demonstrate the **first-generation, hand-rolled PEST slice**
(`KPilotPointParameter`, drain parameter specs, and a custom forward run). They
are kept as a record and **do not run**: that path (`add_parameter`/`build_pst`
and its parameter specs) was deleted on 2026-06-22, along with `mf.DRNFromVector`
and `mf.open_pest_run`, and the Cumberland ones read data that is not in the
repository (executed 2026-10-11, ledger 444). The workflow that replaced
them is the declarative facade:

- **`canonical_05_pest_calibration_setup.ipynb`** — `parameterize` / `observe` /
  `forecast` / `build`, compiling to native `pyemu.utils.PstFrom`.
- **`canonical_06_pest_ies_uncertainty.ipynb`** — `run_ies` + `IesResults`
  (phi convergence, ensembles vs observations, posterior forecast uncertainty).
- **`canonical_04_pest_and_results.ipynb`** — the same facade end to end on the
  canonical model: calibrate, run PESTPP-IES, review the ensemble. Lake, SFR and DRN
  targets (`LakeStageTargets`, `SfrStageTargets`, `SfrFlowTargets`,
  `DrnFlowTargets`) are accepted by `observe(...)` like `HeadTargets`; nothing in
  this folder covers an observation the facade lacks.

`targets_api_quickstart.ipynb` is not part of the slice and does not run either: it
reads a `model` and target files it never builds.
