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
   contributions, parameters at bounds, forecast uncertainty).

Use `CanonicalModelConfig.validation()` for quick work and
`CanonicalModelConfig()` for the full model with at least 10,000 horizontal cells.
The PEST notebooks run PESTPP-IES, which fires hundreds of forward solves — use
parallel `workers=N` and a modest ensemble (they are "go get coffee" cells).
