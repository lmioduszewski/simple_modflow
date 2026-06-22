# Archived PEST notebooks

These notebooks demonstrate the **first-generation, hand-rolled PEST slice**
(`KPilotPointParameter`, drain parameter specs, and a custom forward run). They
still run against the legacy `PestProject.build_pst` path, but the recommended
workflow is now the declarative facade:

- **`canonical_05_pest_calibration_setup.ipynb`** — `parameterize` / `observe` /
  `forecast` / `build`, compiling to native `pyemu.utils.PstFrom`.
- **`canonical_06_pest_ies_uncertainty.ipynb`** — `run_ies` + `IesResults`
  (phi convergence, ensembles vs observations, posterior forecast uncertainty).

The comprehensive legacy multi-target demo (`canonical_04_pest_and_results.ipynb`)
remains in the parent folder because it still covers lake/SFR/DRN named-series
observations that the native facade has not yet absorbed.
