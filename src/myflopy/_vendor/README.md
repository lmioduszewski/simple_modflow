# Vendored dependencies

## figs

Vendored snapshot of the local `figs` project @ commit `bb1526b`, synced
2026-07-16. Do not edit by hand; rerun `python scripts/sync_vendored_figs.py`
(against the figs checkout) after figs changes and commit the result.

This is a TRIMMED closure: the full `mpl/` subpackage (myflopy uses
`plot_cross_section` at runtime) plus the modules providing `Fig`,
`Subplot`, `Template`, `create_hover`. The plotly-side aq-test and
scaling/export stacks are deliberately omitted — they would add
bokeh/cairosvg/reportlab/svglib/svgpathtools as runtime deps. The trimmed
closure requires only myflopy's existing core dependencies (plotly, pandas,
numpy, matplotlib, seaborn, geopandas, shapely). The AESI logo is embedded
as base64 in `figs/_logo_data.py` instead of shipping a PNG.

`myflopy.viz` imports the REAL figs first and falls back to this snapshot,
so the author's machine exercises live figs while installed environments
(and CI) exercise the vendored copy.
