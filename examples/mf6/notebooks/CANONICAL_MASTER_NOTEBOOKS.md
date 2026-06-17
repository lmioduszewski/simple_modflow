# Canonical Master Notebook Set

These notebooks are the authoritative end-to-end examples for `myflopy`.
Every notebook uses `mf.build_canonical_model()` and the same four-layer refined
DISV/Voronoi conceptual model.

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
5. `canonical_04_pest_and_results.ipynb`: full compact truth-to-PEST calibration
   analogue, forward-run validation, optional multi-worker solve, and completed-run review.

Use `CanonicalModelConfig.validation()` for quick work and
`CanonicalModelConfig()` for the full model with at least 10,000 horizontal cells.
