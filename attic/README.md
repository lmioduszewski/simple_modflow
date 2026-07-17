# Attic

Retired modules kept for reference (implementation plan Phase 2.1). Nothing
in `src/`, `tests/`, or `examples/` may import from here; ruff and packaging
exclude this directory. Each file notes why it was retired.

| File | Retired | Why |
|------|---------|-----|
| `grasspyV.py` | 2026-07-16 | hardcoded `grass84.bat`; live GRASS path is `modflow/utils/contour_interp.py` |
| `openET.py` | 2026-07-16 | hardcoded Cumberland paths; unimported |
| `mf3dplots.py` | 2026-07-16 | orphaned 3-D plot scratch; superseded by `interactive_plotting`/`layers` |
| `recovery_analysis.py` | 2026-07-16 | orphaned |
| `cj_approximation.py` | 2026-07-16 | orphaned Cooper-Jacob scratch |
| `gdal.py` | 2026-07-16 | zero importers; broad-except-laden GDAL helpers |
| `gis_functions.py` | 2026-07-16 | zero importers; top-level osgeo import |
| `gwt.py` | 2026-07-16 | transport scratch with personal paths; the real transport home is plan Phases 5.3 + 6.1 |
| `rasterize_points.py` | 2026-07-17 | extracted from `modflow/utils/surfaces.py` (plan 4.3); zero importers — `InterpolatedSurface` moved to `mf6/grid/interpolated_surface.py` |
