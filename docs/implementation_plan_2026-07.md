# myflopy Consolidation & Completion Plan

**Date:** 2026-07-06 (rev. 3 — adds Phase 6: full GWT / GWE / PRT / PEST-IES integration
into the preferred API incl. visualization + hover; folds in work completed 2026-07-06)
**Source:** Full structural/code review of `src/myflopy` on branch `myflopy` (follow-up to
`docs/refactor_review_report.md`, 2026-06-12). Every finding below was verified against the
code; line numbers are anchors and WILL drift — always locate by symbol name, not line
number.
**Audience:** an implementing agent (e.g. Opus 4.8) working phase-by-phase. Read the whole
"Ground rules" section before touching anything.

---

## Completed since rev. 2 (2026-07-06, on branch `myflopy`, uncommitted) — do NOT redo

These subsystems now exist and are the *patterns to extend* in later phases:

1. **Sectioned hover system** — `src/myflopy/modflow/utils/datatypes/hover.py`:
   `HoverSpec` / `HoverStyle` / `Fields` / `LayerTable` / `HoverContext`, plus factories
   `head_hover`, `lak_hover`, `sfr_hover`, `cell_input_hover` (blue accent `#185FA5`),
   `result_hover` (teal `#0f6e56`), `compare_hover` (purple `#534AB7`),
   `surface_water_hover`. Per-cell strings are precomputed into `customdata`; the
   hovertemplate is a styled skeleton of `%{customdata[i]}` refs. Every choropleth path
   (single-model heads/inputs/results/stage, group member maps, group `diff()` maps)
   sets a default via `kwargs.setdefault("hover_spec", ...)`; call-site sugar
   (`hover=`, `hover_layers=`, `hover_surfaces=`, `hover_fields=`) resolves centrally in
   `Choro._resolved_hover_spec()`. LAK/SFR q-maps join per-feature **stage** via
   `join_lak_stage`/`join_sfr_stage` (package_surface_water.py). Top-level exports:
   `mf.HoverSpec`, `mf.HoverStyle`. Tests: `tests/test_hover_spec.py` (24),
   `tests/test_hover_integration.py` (11 slow, canonical).
2. **Colorscale policy** (user decision, pinned in `tests/test_colorscale_policy.py`):
   red/white/blue diverging ONLY for (a) signed gaining/losing "q"-like exchange fields —
   SFR/LAK/combined-SW use `_blue_white_red_diverging_colorscale()` (gaining/negative =
   BLUE, losing/positive = RED); per-package `q` results keep `"RdBu"` — and (b) ALL
   diff/compare maps (`"RdBu"` + `zmid=0` → negative RED, positive BLUE). **Everything
   else defaults to `'earth'`** (Plotly brown→cream→blue; the mounding-figure look).
   Centralized in `package_registry.py` specs + call-site fallbacks. `Choro.colorscale`
   setter now accepts explicit `[[pos, color], ...]` lists (was silently dropping them —
   SFR's documented scale had never rendered) and `plot_mpl` builds a
   `LinearSegmentedColormap` from list scales.
3. **Map view fitting + live pan/zoom sync** — `viz.mosaic` map subplots and
   `animate(kind="map")` are framed to the data (union extent) via `Choro.map_view()` /
   `Choro.latlon_bounds` / `viz.shared_map_view()`; `sync_views=True` (default) also
   injects a `plotly_relayout` JS handler so panels pan/zoom together on
   `show()`/`write_html()` (not notebook-inline). `sync_views=False` = same shared start
   view, independent interaction.
4. **figs backend changes (in the SEPARATE figs repo, `_fig.py`, UNCOMMITTED):**
   `Fig.add_post_script(js)` + threading into `show()`/`write_html()`; `write_html`
   signature fixed to plotly-compatible `write_html(self, file=None, config=None, ...)`
   (was mis-binding the path to `config`). **Phase 1.1's vendored snapshot must be taken
   AFTER these are committed in figs** — they are load-bearing for map sync.
5. The public diff surface is the ONE **`diff()` verb** (`group.diff().hds.map("b")`,
   `group.diff().packages.ghb.results.q.map("b")`); `compare_map` is internal plumbing —
   never document, teach, or test it directly.

---

## Resolved decisions (locked by the user)

Do not re-litigate these or present alternatives; implement as stated. If a *new*
ambiguity arises that this plan does not cover, stop and ask the user — do not guess.

| ID | Decision | Resolution |
|----|----------|------------|
| D1 | `figs` dependency | **Vendor a snapshot** into `src/myflopy/_vendor/figs/` with external-first import. Do NOT publish figs to PyPI; do NOT build a no-figs degraded mode. Details in Phase 1.1. |
| D2 | Executed notebooks | **Strip all outputs from every tracked notebook.** Rendered copies live outside git. Details in Phase 2.3. |
| D3 | Git-history rewrite | **Not part of this plan.** The agent must never rewrite history or force-push (Appendix D is user-only reference). |
| D4 | Structured-grid stance | **Add thin `mf.dis` / `mf.disu` passthroughs** (model-type-aware). `GridSpec.structured` stays fail-fast. Details in Phase 5.5. |
| D5 | Plotting consolidation | **Do it**, as its own late phase (Phase 8), after the god-module splits (Phase 4), API completion (Phases 5–6), and the `_flopy_compat` boundary (Phase 7.1). Details in Phase 8. |
| D6 | Hover defaults | DONE (see "Completed"). Heads default `layers="active+strip"`; sectioned/styled hover is the default everywhere; `custom_hover` stays as the raw escape hatch. |
| D7 | Colorscales | DONE (see "Completed"). Diverging only for signed q-like + diff maps; `'earth'` for everything else. New surfaces added by this plan MUST follow this policy. |

---

## Ground rules for the implementing agent

Non-negotiable project conventions. Violating them is a review failure even if tests pass.

1. **Environment.** mf-env venv: `C:\Users\lukem\Python\mf-env\.venv\Scripts\python.exe`
   with `PYTHONPATH=src`. Fast loop: `pytest -m "not slow"` (~60 s, 540 tests as of
   2026-07-06 evening). Full suite: ~14 min. Fast suite after every change-set; FULL
   suite at the end of each phase.
2. **Tests are mandatory.** Real, runnable pytest tests covering behavior (not smoke
   tests) for every capability/refactor; run pytest to prove they pass before declaring
   done. Refactors need existing tests green AND a new test pinning the seam.
3. **Never change public import paths.** Splits/moves happen behind facades (the proven
   pattern: `package_explorer.py` fronting the `package_*` family). The lazy-export
   tables in `myflopy/__init__.py` and `modflow/mf6/__init__.py` must keep resolving.
4. **Grep both API layers before building.** Package-first (`package_api.py`, `specs.py`,
   `advanced.py`, `builders.py`, `geopackage.py`) sits on the OO engine
   (`modflow/mf6/*.py` + `simulation/`). `docs/myflopy_context.md` is the capability map —
   update it (and `CLAUDE.md`) whenever a phase changes what exists.
5. **Zero-based stress periods** in every myflopy-facing API.
6. **ASCII-safe console output** (Windows console; see commit `e61b973`).
7. **Do not reintroduce retired APIs** (legacy PEST `add_parameter`/`build_pst` path,
   retired PEST demo models).
8. **One phase per branch/PR.** A phase is done when all its acceptance criteria pass
   plus a full-suite run.
9. **Check the working tree first.** As of rev. 3 the tree carries the completed-but-
   uncommitted hover/colorscale/map-sync work listed above PLUS the earlier in-flight
   diff/group work. Land it with the user before starting any phase.
10. **Verify before acting.** Re-verify every claim (orphan file, symbol, line) with grep
    before deleting or renaming — the codebase moves fast.
11. **No history rewrites, no force pushes, ever** (D3).
12. **New map surfaces follow D6/D7:** every new choropleth gets a `HoverSpec` default
    (via `kwargs.setdefault("hover_spec", ...)`) and a policy-compliant colorscale, with
    tests in the style of `test_hover_spec.py` / `test_colorscale_policy.py`.
13. **Plotting goes through `myflopy.viz`** (Fig/subplots/mpl_axes/PALETTE); the
    documented raw-plotly exceptions are 3-D scenes, mapbox debug maps, and animation
    re-wraps.

---

## Phase 0 — Baseline (prerequisite, ~0.5 day)

**Steps:**
1. Land or shelve ALL in-flight work with the user: the diff/group feature files
   (`model_diff.py`, `model_group.py`, `model_results_diff.py`, `run_model.py`,
   `package_*.py`, `simulation/*.py`, `headsplus.py`, `xsections.py`, their tests) AND
   the rev-3 completed work (hover.py + wiring, colorscale policy, viz.py map sync,
   `test_hover_*.py`, `test_colorscale_policy.py`, `test_view_grammar_composers.py`
   updates). Commit the **figs repo** changes (`_fig.py`: `add_post_script`,
   `write_html` fix) in figs itself.
2. Record the baseline fast/full suite results next to the PR.
3. Note the starting commit.

**Acceptance:** clean `git status` in BOTH repos; recorded baseline test results.

---

## Phase 1 — Distribution blockers (highest priority, ~1–2 days)

### 1.1 Vendor `figs` with external-first import (D1)

**Problem (verified):** `src/myflopy/viz.py` hard-imports
`from figs import Fig, Subplot, Template, create_hover` and
`from figs.mpl import REPORT, Theme, get_mplfig`. `figs` resolves to the user's local
project (`C:\Users\lukem\Python\Projects\figs`, package root `src/figs/`) and is not
declared in `pyproject.toml`, not on PyPI. 21+ modules import `myflopy.viz`, so a pip
install of myflopy is broken for anyone but the author.

**Implementation:**
1. **Snapshot the source** into `src/myflopy/_vendor/figs/` — ONLY the modules providing
   the symbols above plus their figs-internal deps, traced transitively; rewrite absolute
   `figs.` imports to relative so the subtree is self-contained. **Snapshot only after
   the figs repo has committed** `add_post_script`, the `write_html(file=None, ...)`
   signature fix, `_repr_mimebundle_` config merge, and `add_logo=False` default — all
   load-bearing for myflopy's map sync + inline scroll-zoom. Add
   `_vendor/__init__.py` and `_vendor/README.md` ("Vendored snapshot of figs, synced
   <date>. Do not edit by hand; run `scripts/sync_vendored_figs.py`.").
2. **External-first import in `viz.py`:**
   ```python
   try:
       from figs import Fig, Subplot, Template, create_hover
       from figs.mpl import REPORT, Theme, get_mplfig
   except ImportError:  # vendored fallback for installed environments
       from myflopy._vendor.figs import Fig, Subplot, Template, create_hover
       from myflopy._vendor.figs.mpl import REPORT, Theme, get_mplfig
   ```
   `viz.py` stays the ONLY runtime importer of figs (TYPE_CHECKING imports exempt);
   enforce with an AST/grep test.
3. **Sync script** `scripts/sync_vendored_figs.py`: copies the traced module set,
   rewrites imports, stamps the sync date. Idempotent.
4. **Packaging:** verify `myflopy._vendor.figs` lands in the wheel; add any figs data
   files to package-data if they exist.

**Tests:** `tests/test_vendored_figs.py` — vendored import + `Fig()` construction; a
subprocess with figs hidden proving `import myflopy.viz` works and map-sync post-scripts
still inject into `write_html` output; single-importer rule. CI (1.3) has no external
figs, so every CI run exercises the vendored path — the permanent drift alarm.

**Acceptance:** fresh venv, `pip install -e .`, NO figs on path → `import myflopy.viz`
works, `viz.Fig()` constructs, viz-touching fast tests pass.

### 1.2 Declare `matplotlib` (and audit the rest)

Add `"matplotlib>=3.8"` to `[project] dependencies` (imported top-level by 10+ live
modules, currently only transitive). Audit the rest against a fresh-venv install:
`seaborn` (used by `viz.mpl_axes`!), `openpyxl` (xlsx reads), and document that FloPy's
`Triangle` and the MODFLOW/PEST++ executables are runtime binaries, not pip deps.

**Acceptance:** fresh-venv install + import-surface smoke (1.3) green.

### 1.3 Add CI (GitHub Actions)

`.github/workflows/ci.yml`: matrix ubuntu-latest + windows-latest × Python 3.10/3.12;
`pip install -e .[dev]` → `pytest -m "not slow" -q` → (after 1.4) `ruff check src tests`.
Fast suite must not need MF6 binaries (that is the `slow` contract). First run will
surface optional-dep assumptions (pyvista/trame, pymetis/h5py, pyemu) — prefer
`pytest.importorskip` over fattening the CI install. Add the import-surface smoke step:
`python -c "import myflopy; [getattr(myflopy, n) for n in myflopy.__all__]"`.

**Acceptance:** green CI on both OSes on a no-op PR.

### 1.4 Lint + (scoped) type checking

`[tool.ruff]`: line-length 100, target py310, `extend-exclude = ["attic", "examples",
"scripts", "src/myflopy/_vendor"]`, `select = ["E","F","W","I","UP","B"]`. First pass:
safe autofixes only; triage the rest; never reformat in the same commit as logic changes.
Scoped mypy (or pyright basic) on the declarative core (`specs.py`, `package_api.py`,
`sources.py`, `workspace.py`, `advanced.py`, `builders.py`) — `py.typed` ships today with
nothing verifying it. Add tools to the `dev` extra.

**Acceptance:** `ruff check src tests` clean in CI; type check green on the scoped list.

---

## Phase 2 — Junk removal & repo hygiene (~1 day)

### 2.1 Delete/attic the orphaned modules

Verified orphans (re-verify each with grep before removing):

| File | Action | Notes |
|------|--------|-------|
| `modflow/utils/shiny_app.py` | delete | palmerpenguins Shiny demo |
| `modflow/mf6/pyqgis.py` | delete | mutates `sys.path` with hardcoded OSGeo4W path |
| `modflow/mf6/grasspyV.py` | attic | hardcoded `grass84.bat`; live GRASS path is `utils/contour_interp.py` — confirm |
| `modflow/utils/openET.py` | attic | hardcoded Cumberland paths |
| `modflow/utils/iheads.py` | delete | superseded by `HeadsPlus` |
| `modflow/mf6/mf3dplots.py` | attic | update the `viz.py` docstring that names it |
| `modflow/mf6/recovery_analysis.py` | attic | orphaned |
| `modflow/calcs/cj_approximation.py` | attic | orphaned |
| `modflow/gwt/` (whole dir) | attic `gwt.py`, delete dir | scratch with personal paths; the REAL transport home is Phases 5.3 + 6.1 — nothing may import `myflopy.modflow.gwt` |

`prism_ppt.py` (1 importer): cascade-check. `mf2Dplots.py` (1 importer) is deferred to
Phase 8's fold — do not touch here.

### 2.2 Strip `__main__` scratch blocks with personal paths

Sites: `ghb.py`, `utils/raster.py`, `utils/prism_ppt.py`, `utils/surfaces.py`. Remove the
`C:\Users\lukem\...` blocks; convert anything worth keeping into a test or `examples/`
script. **Acceptance:** `grep -rn "lukem" src/myflopy` empty.

### 2.3 Strip all notebook outputs from git (D2)

`nbstripout` in dev extra + `.pre-commit-config.yaml` (+ fallback
`scripts/strip_notebooks.py`); one-time strip-only commit over ALL tracked notebooks
(incl. `notebooks/archive/`); gitignored `docs/_rendered/` + `scripts/render_notebooks.py`
for browsable executed copies; evaluate `starting_heads.csv` (604 KB fixture?); gitignore
`scripts/scratch*`.

**Acceptance:** `nbstripout --dry-run` clean over tracked notebooks; none over ~100 KB.

### 2.4 Move MP3DU executables out of `src/`

`src/myflopy/modflow/mp3du/*.exe` (10.7 MB) → untracked `tools/mp3du/` (tracked README
tells the user what goes there). Resolution order in `particles.py`: explicit arg →
`MYFLOPY_MP3DU_DIR` env var → repo `tools/mp3du/` → deprecated `_MODULE_DIR` fallback
(warns via 3.1). gitignore `src/**/*.exe`. Unit-test the resolution order (monkeypatched
env, tmp dirs).

---

## Phase 3 — Deprecation mechanics & API-surface cleanup (~1 day)

### 3.1 Deprecation helper + policy

`src/myflopy/_deprecation.py` with `warn_deprecated(old, new, *, since)` and
`deprecated_module_getattr(mapping, module)`; `docs/deprecation_policy.md` (short):
every compatibility name warns with its replacement; names live ≥ 2 tagged releases;
`__compatibility__` in `myflopy/__init__.py` is authoritative. Tests: `pytest.warns`
per alias; fast suite passes with `-W error::DeprecationWarning` filtered to `myflopy.*`.

### 3.2 Resolve the GHB/DRN naming collision + FromVector consistency

Rename so `*FromVector` is the real class name in all four modules (`ghb.py` GHB→
GHBFromVector, `drn.py` DRN→DRNFromVector; chd/kflow already correct); keep warned
aliases via 3.1. Repoint `modflow/mf6/__init__.py`'s `"GHB"` export; update internal
import sites and `mfsimbase.py`; update docs.

---

## Phase 4 — Structural splits & import hygiene (~4–6 days)

**Mechanic:** create the new package, move code in small mechanical commits, keep the old
module as a re-export facade, change zero user-facing paths, run the fast suite after
each move. No renames, no "improvements" during a split.

### 4.1 Split `modflow/mf6/observations.py` (2,744 lines, 13 classes)

Target package `modflow/mf6/observations/`: `_shared.py` (normalization/series helpers),
`heads.py` (`HeadTargets`, `BoundHeadTargets`, `BoundHeadTargetPlots`), `lake.py`,
`sfr.py` (stage+flow, bound), `drn.py`, `plots.py` (`BoundNamedSeriesTargetPlots`),
`registry.py` (`TargetRegistry`), `__init__.py` re-exporting everything (the package
shadows the old module path automatically). `heads_observations.py` is NOT part of this
split. Consumers (top-level `_EXPORTS`, pest package, `model.targets`) keep working —
verify with grep. Watch pickling references. Tests: existing target tests unchanged +
`tests/test_observations_layout.py` (both import paths).

### 4.2 Split `project/model_group.py` (2,572+ lines, 28+ classes)

Target package `project/group/`: `spatial.py` (`_GroupSpatialView`, `GroupHeads`),
`budget.py`, `inputs.py`, `results.py` (`GroupCellPackageResults` + namespaces +
`GroupOutputs`), `sfr.py`, `lak.py` (incl. connections), `uzf.py`, `surface_water.py`,
`packages.py`, `core.py` (`ModelGroup`), `__init__.py`. Keep `model_group.py` as facade.
**Re-verify the class inventory after Phase 0** (in-flight work adds classes). Tests:
group/diff test files green unchanged + layout test.

### 4.3 Consolidate the three `surfaces` modules

`myflopy/surfaces.py` (engine — KEEP), `grid/surfaces.py` (sampler — KEEP),
`modflow/utils/surfaces.py` (legacy; `InterpolatedSurface` still imported by
`layers.py`, `surface_data.py`, `xsections.py`, `test_layers.py`). Extract
`InterpolatedSurface` to a proper home (`grid/interpolated_surface.py`, or into the
engine ONLY if semantics align — read first); repoint the four sites; leave a warned
facade; attic the rest symbol-by-symbol.

### 4.4 Import layering + deferred-import ratchet

~195 function-level `from myflopy...` imports (≈half of internal imports) are circular-
dependency scar tissue centered on `simulation/base.py`. Document the layer map in
`docs/architecture_layers.md`:
- **L0** `_optional`, `_deprecation`, `_vendor`, `viz`, `modflow/utils/*` (incl.
  `datatypes/hover.py` — it is a leaf)
- **L1** `modflow/mf6/grid/*`
- **L2** output readers: `headsplus`, `heads_observations`, `budget`, `budget_tables`
- **L3** explorers + observations + `*_plotting` (→ `myflopy/plot/` after Phase 8)
- **L4** engine facade: `simulation/*`, OO builders
- **L5** declarative: `specs`, `sources`, `geopackage`, `package_api`, `advanced`,
  `builders`, `layers`, `surfaces`, `grid_spec_resolver`
- **L6** lifecycle: `workspace`, `project/*`, `pest/*`, `parallel`, `prt`, `mp3du`
Runtime imports only point at same-or-lower layers; upward = TYPE_CHECKING or injection.
`tests/test_import_layering.py` (AST walk vs the map, `_vendor` excluded) + a deferred-
import **ratchet** (`tests/deferred_import_allowlist.json`, counts may only go DOWN).
Burn down ~30 easy conversions to prove the mechanism.

### 4.5 Second-tier splits (after 4.1/4.2)

`grid/triangle.py` (1,742) → move domain/feature specs + cleanup/optimization
orchestration into existing `grid/` collaborators behind the `TriangleGrid` facade.
`package_plotting.py` / `package_surface_water.py` stay in the explorer family unless
they keep growing. `specs.py` is cohesive — leave below ~2,500.

---

## Phase 5 — Package API completion (~5–8 days, independent of Phase 4)

### 5.0 The pattern to replicate (read first)

Every list-style boundary package = four pieces (DRN is the reference):
1. `package_api.py` helper class (`_DRNPackage`: `__call__` direct data / `.gpkg` GIS /
   `.flopy` escape hatch) + module-level singleton.
2. `geopackage.py` resolver method (`GeoPackageSource.drn`) — shared plumbing
   (`_boundary_data`, `_cells`, `_periods`, surface references, `CellSurfaceOffset`).
3. `advanced.py` spec factory (`drn_spec`).
4. Registry + exports + docs: `FieldSpec`/`ResultSpec` in `package_registry.py` (this is
   what feeds explorers, ModelDiff, **hover defaults, and the colorscale policy** — new
   entries MUST set `colorscale="earth"` for non-signed fields, `"RdBu"` for signed q),
   lazy-export entries, `docs/myflopy_context.md` + `CLAUDE.md` updates.

**Tests per package:** `test_package_api.py` / `test_geopackage_source.py` /
`test_advanced_specs.py` patterns + a slow e2e on a small grid where feasible + hover/
colorscale assertions in the style of `test_hover_integration.py`.

### 5.1 `mf.riv` — do first
RIV records `(cellid, stage, cond, rbot)`. All four pieces; `.gpkg` with surface
references for stage/rbot. Input hover: `cell_input_hover("stage", extra_fields=("cond",
"rbot"))`; colorscale `earth`; result q joins the RdBu q family.

### 5.2 `mf.evt`
EVT `(cellid, surface, rate, depth [, pxdp, petm, petm0])`; `nseg` via options. Reuse the
RCH areal-mapping prior art (`RCHBuilder`, `GeoPackageSource.rch`). Earth colorscale,
blue input hover.

### 5.3 GWT/GWE package-first helpers + model-type-aware core helpers

**Part A — model-type dispatch (prerequisite for 5.5 and 6.1/6.2):** `mf.ic`/`mf.oc`/
`mf.disv` hardcode `ModflowGwf*` classes (verified) so they cannot serve GWT/GWE models.
Add module-level dispatch builders (serialization-safe — `_callable_ref` resolves
module-level functions):
```python
_IC_CLASSES = {"gwf": flopy.mf6.ModflowGwfic, "gwt": flopy.mf6.ModflowGwtic,
               "gwe": flopy.mf6.ModflowGweic}
def _build_ic(model, **values): ...
```
Determine model kind from the flopy instance (check what flopy 3.10/3.11 exposes; prior
art: `specs.py` `ModelType`, `_require_model_type`). Round-trip serialization tests.

**Part B — thin package factories** in `package_api.py` (same shape as `mf.ic`):
- GWT: `mf.adv(scheme=)`, `mf.dsp(alh=, ath1=, ...)`, `mf.mst(porosity=, ...)`,
  `mf.ssm(sources=)`, `mf.cnc(...)` (list BC: `()`/`.flopy`; `.gpkg` later),
  `mf.src(...)`, `mf.ist(...)`.
- GWE: `mf.est(...)`, `mf.cnd(...)`, `mf.ctp(...)` (list BC), `mf.esl(...)`.
Update `examples/mf6/variant_workflow/03_coupled_gwf_gwt.py` + hand-rolled test specs to
use them (keep one raw-PackageSpec regression test). Registry `FieldSpec` entries for
`cnc`/`ctp` (blue input hover, earth). Do NOT resurrect `modflow/gwt/`.

### 5.4 `mf.maw` and `mf.hfb`
- **MAW:** high-level `mf.maw(wells=[...], context=)` computing connections from point +
  screen interval against layer elevations (new `MAWBuilder` in `modflow/mf6/maw.py`,
  UZFBuilder-shaped) + `.flopy` raw form. Slow e2e on the canonical grid.
- **HFB:** cell-*pair* records; on DISV pairs must share a face crossed by the barrier
  line. RESEARCH STEP first: read `grid/voronoi.py` + `grid/connectivity.py` for a
  shared-edge-with-geometry API; if absent, build the edge lookup there with unit tests,
  then `mf.hfb.gpkg(path, hydchr=...)`.

### 5.5 `mf.dis` / `mf.disu` passthroughs (D4)
Thin factories mirroring `mf.disv`, using 5.3A dispatch (Gwf/Gwt/Gwe/Prt dis classes;
note `prt.py` already builds `ModflowPrtdis(v)` internally — don't duplicate).
`GridSpec.structured` stays fail-fast (`test_gridspec_fail_fast.py` untouched).
Docstrings state the Voronoi-first stance. Explorer/registry wiring only if the
explorers can render structured grids — verify, else document DISV-only.

### 5.6 YAML/TOML spec serialization
`SimulationSpec.to_yaml/from_yaml` (+ optional TOML via stdlib `tomllib` + `tomli-w`
extra) in a small `specs_io.py`; PyYAML (`safe_load`/`safe_dump` only) as a core dep;
`Project.add_simulation_from_yaml`. Normalize through `_json_value` first. Round-trip
tests incl. exchanges/hooks/GridSpec/refs + a hand-written minimal YAML that builds +
an example file in `examples/`.

### 5.7 (Forward-looking, larger) grid-lazy GIS packages for deferred GridSpec
`.gpkg` helpers first: when the grid is deferred, return a `PackageSpec` whose
module-level builder resolves against `model.myflopy_context` at build time, with the
declaration payload (paths/params) in `PackageSpec.options` (serializable → composes
with 5.6: a full GIS model from one YAML). Then UZF (hardest — `UZFBuilder.build()`
bakes arrays), then SFR/LAK. Clear error when a lazy GIS package has no grid.

### 5.8 PEST parameterization backlog (build-side; the viz side is 6.4)
In priority order: UZF as a `parameterize` target; raster-driven zone arrays for
`style="zone"`; Tikhonov/preferred-value regularization helpers; sensitivity/
identifiability helpers. Facade COMPILES to native pyEMU PstFrom — never a parallel
engine.

---

## Phase 6 — Full model-type integration: GWT, GWE, PRT, PEST-IES (~6–10 days)

**Goal:** every model type myflopy claims (`mf.gwf/gwt/gwe/prt`) and the calibration
layer get the SAME first-class treatment heads already have: a results accessor in the
unified grammar (`map/plot/xs/mosaic/animate`), sectioned hover (D6), policy colorscales
(D7), group/diff support, and observation targets. **Verified starting state:** GWT/GWE
have NO results tier at all (no concentration/temperature reader exists anywhere in
`src/`); PRT has `PRTRunResults` (pathlines df, 3-D `scene`, `plot_map`,
`export_3d_html`) + `model.particle_tracking` but no grammar/hover/choropleth
integration; PEST-IES has a rich review layer (`IesResults`: `plot_phi`, `plot_vs_obs`,
`plot_phi_distribution`, `plot_parameters_at_bounds`, `plot_phi_contributions`,
`plot_prior_vs_obs`, `plot_conflict`, `field`/`plot_field`, `report`, backend switch)
that predates the hover/colorscale systems.

### 6.0 Design keystone: ONE generic dependent-variable surface

Heads, concentration, and temperature are the same shape: a binary output file read per
(kstpkper, layer, cell), explored through the same grammar. **Do not clone
`HeadsPlus` twice.** Instead:
1. Read `HeadsPlus` (`modflow/mf6/headsplus.py`) and factor the file-agnostic core into
   a parameterized base (or make `HeadsPlus` itself take `text`/`label`/`unit`
   parameters). FloPy's `HeadFile` reads concentration/temperature binaries via
   `text="concentration"` / `"temperature"` (verify against the installed flopy;
   concentration output is declared by GWT OC `concentration_filerecord`, `.ucn` by
   convention).
2. Instantiate: `model.hds` (existing, unchanged public behavior), `model.conc` (GWT),
   `model.temp` (GWE). Each carries: `.array()`, `.wide()/.long()`, the grammar verbs
   `map/xs/plot/mosaic/animate`, contours, and the hover sugar
   (`hover_layers`/`hover_surfaces` — `LayerTable` is already field-generic, so
   per-layer concentration tables + surfaces merge work for free).
3. `SimulationBase`/`ModelView`/`Run.model()` must recognize GWT/GWE models: today
   `Run.model()` documents "the preferred myflopy **GWF** model view" — extend model-kind
   detection so a transport/energy model view exposes `conc`/`temp` (and NOT `hds`), and
   `run.model("transport")` works. Read `workspace.py:Run.model` + `ModelView` first.

### 6.1 GWT results tier (concentration)

Building on 6.0:
1. **Hover:** `conc_hover(unit="mg/L")` factory — primary `conc`, title
   "Concentration", `layers="active+strip"` default, surfaces merge supported; unit
   configurable (model-dependent; do NOT hardcode — thread from a `units=` param with a
   documented default). Registry-style defaults on every conc map.
2. **Colorscale:** per D7, `'earth'` default for concentration; Δconc (diff maps) =
   `RdBu`, negative red / positive blue. Add to `test_colorscale_policy.py`.
3. **Budget:** GWT budget terms (STORAGE-AQUEOUS, SSM, DECAY, …) through the budget
   explorers — read `budget.py`/`package_budget.py`/`budget_tables.py` and extend the
   term discovery to transport CBC files; `model.bud` equivalent on the GWT view.
4. **Group/diff:** `GroupConc` mirroring `GroupHeads` (member maps + `diff()` Δconc maps
   via the existing compare payload machinery — it is value-column generic) +
   `compare_hover("conc", "diff", ...)`. Slot into Phase 4.2's `project/group/` layout
   if that split has landed (`group/conc.py`).
5. **Observations:** `ConcTargets` (measured vs simulated concentration at points)
   mirroring `HeadTargets` — reuses the 4.1 observations layout; wire into
   `cal.observe(...)`/`cal.forecast(...)` + `forward_run.py` post-processors so
   **transport calibration** works end-to-end (this is the PEST tie-in).
6. **Canonical transport fixture:** tests need a real GWT model. Add a small coupled
   GWF+GWT variant fixture — promote `examples/mf6/variant_workflow/03_coupled_gwf_gwt.py`
   into a conftest fixture (slow-marked) or add a `transport=True` option to a small
   canonical config. Every 6.1 feature gets a slow e2e against it (map renders, hover
   template asserts, budget table, group diff, obs round-trip).

### 6.2 GWE results tier (temperature)

A second instantiation of 6.0 — deliberately mechanical after 6.1:
`model.temp` (`text="temperature"`), `temp_hover(unit=...)` (thread model units;
document default), `'earth'` scale + `RdBu` diffs, GWE budget terms, `GroupTemp` +
diff, optional `TempTargets`, and a small GWF+GWE fixture (mirror the GWT one; the
exchange builder `build_gwf_gwe_exchange` already exists). If 6.1's factoring was done
right, this phase is mostly registrations + tests; if it is not mostly mechanical, stop
and fix the 6.0 abstraction instead of copy-pasting.

### 6.3 PRT integration

PRT results are trajectories, not per-cell fields — integrate where the grammar fits and
style the rest, rather than forcing everything into choropleths:
1. **Cell-based derived maps (the real win, grammar-compatible):** from
   `PRTRunResults.pathlines` / `terminal_points`, build per-cell summaries and render
   them as first-class `Choro` maps with hover + policy colorscales:
   - `results.travel_time_map(stat="median")` — travel time per terminating cell
     (classic capture-zone / time-of-travel figure); `'earth'` scale (or reversed —
     verify readability), `result_hover("travel_time", units={"travel_time": "d"})`
     with extra fields (particle count, min/max time).
   - `results.endpoints_map()` — particle-termination counts per cell.
   - `results.capture_map(by="release_group")` — which release group's particles end
     where (categorical → needs a small categorical-colorscale story; keep simple:
     one map per group via `viz.mosaic`, synced views).
   These are per-cell payloads → they inherit mosaic/animate/hover for free once built
   through `build_cell_input_map_payload`-style tables.
2. **Pathline map hover:** `plot_map`/pathline traces are scatter, not choropleth — the
   hover engine is trace-agnostic (customdata + template), so add a
   `pathline_hover()` spec rendering particle id, release point/group, current time,
   layer; apply `HoverStyle.to_hoverlabel()` to the scatter traces. Read
   `prt.py:PRTRunResults.plot_map` + `interactive_plotting.py` scene builders first.
3. **Build-side parity:** add `mf.mip(porosity=)` and `mf.prp(...)` thin factories
   (5.3A dispatch pattern) so a PRT model is fully declarable in a `SimulationSpec`
   without raw PackageSpecs; keep `model.particle_tracking.prt(...)`/`PRTProject` as the
   sanctioned high-level runtime path (it manages FMI/grid copying correctly — do not
   duplicate that logic in specs).
4. **3-D scene** stays raw plotly (documented viz.py exception). MP3DU untouched except
   Phase 2.4.
5. Group story: comparing PRT runs across models = diff of travel-time maps (mechanical
   via the compare payload once 6.3.1 exists) — mark forward-looking, implement only the
   single-model maps + hover now.

### 6.4 PEST-IES visualization & preferred-API integration

The review layer exists and is good; bring it up to the house standards and close the
workshop gaps (memory: spatial parameter-field maps were the biggest gap; DBTL framing +
input-capture-as-observations is the driving workflow):
1. **Parameter-field maps get hover + policy colors.** `IesResults.plot_field(target)`
   renders per-cell parameter stats — read its implementation (it routes through
   `build_choropleth` per the viz memory). Add `parameter_field_hover()`: primary = the
   plotted stat (e.g. posterior mean multiplier or K), block: prior mean, posterior
   mean, posterior sd (per cell, from the ensembles via `field()`), footer: iteration,
   realization count. Colorscale policy: **multiplier fields are relative → diverging
   `RdBu` centered at 1 (log-centered)** — negative/red = reduced, blue = increased —
   matching the user's diff-map rule; **absolute-K fields → `'earth'` + logscale**.
   Pin in `test_colorscale_policy.py`.
2. **Uncertainty maps:** `field_uncertainty_map(target)` — posterior sd or
   prior→posterior variance-reduction per cell (sequential `'earth'`), with hover.
   This is the "did the data inform this region" figure.
3. **Prior-vs-posterior mosaics:** `IesResults.field_mosaic(target, which=("prior",
   "posterior"), stat="mean")` — sugar over `viz.mosaic` (shared color scale, synced
   views come free now). One line of real code per panel; mostly tests + docstring.
4. **Residual map:** posterior mean residual per observation location on the grid
   (bubble scattermap or nearest-cell choropleth) — `plot_obs_residuals(map=True)`.
5. **Ensure every IES plot uses the viz front door** (`viz.Fig`/`mpl_axes`) — audit
   `ies.py`'s figure construction; it already has `backend=` switches, so this is
   verification + spot fixes, not a rewrite.
6. **Preferred-API surface check:** `model.pest(name, ...)` front door,
   `model.pest_runs` / `run.pest_runs` → `.review()` — already the intended path;
   document in the capability map that `PestProject` direct construction is advanced.
   The parameterization backlog itself is Phase 5.8.
7. **Tests:** fast tests with a synthetic ensemble fixture (small DataFrames standing in
   for pyemu ensembles — `ies.py` methods are DataFrame-driven; check what
   `test_mf6_pest.py` already fakes) + one slow end-to-end extension of the existing
   IES canonical test asserting the new maps render with hover + policy colors.

### 6.5 Capability-map + manual updates

After each of 6.1–6.4: update `docs/myflopy_context.md`, `CLAUDE.md` (the "multi-physics"
claim becomes real; remove stale gap lines), and add one canonical notebook or example
per surface (a transport-results walkthrough; a PRT travel-time map example; an IES
field-map review section in the existing PEST notebooks).

---

## Phase 7 — FloPy 4 readiness & robustness (~1–2 days)

### 7.1 `_flopy_compat.py` — centralize flopy-internal touchpoints

Verified touchpoints: `flopy.utils.binaryfile` (5 sites — grows with 6.0's
concentration/temperature readers: route them through the compat module from day one),
`flopy.plot(+.crosssection/.plotutil)`, `flopy.utils.triangle`, `flopy.export.vtk`,
`flopy.utils.voronoi`, `flopy.utils.geospatial_utils`, `flopy.discretization.vertexgrid`,
the `patch_simulation_plot` monkeypatch (`project/run_model.py`), and one
internal-deprecated `mfmodel` call (warns in test output — find and fix/isolate).
Create `src/myflopy/_flopy_compat.py` re-exporting every non-`flopy.mf6` symbol with a
one-line "used for" comment; mechanically repoint importers; AST test
(`tests/test_flopy_compat_boundary.py`): only `_flopy_compat.py` may import
`flopy.utils/plot/export/discretization`. Keep the `flopy>=3.10,<4` pin. **7.1 lands
BEFORE Phase 8** so plotting files move exactly once with clean imports.

### 7.2 Logging

`src/myflopy/_logging.py` (`get_logger`, NullHandler on the `myflopy` root — libraries
never configure handlers). Replace `print(...)` in live modules (grid `selection.py`,
`triangle.py`, `pest/project.py`, `simulation/runtime.py`, `choros.py` — incl. the
colorscale-fallback print, `readers.py`, `grid/surfaces.py`). Exemptions: deliberate
console **reports** (`ModelDiff.report()` and friends — keep printing, keep ASCII) and
`_vendor`. Test: `import myflopy` emits nothing to stdout.

### 7.3 Narrow the silent exception swallows

Sites: `grid_spec_resolver.py`, `workspace.py` (`model_names` fallback),
`utils/gdal.py`, `surface_water_validation.py`, `simulation/base.py`, plus the two added
knowingly in the hover work (`Choro._build_hover_context`'s per-dates and topbtm
guards — narrow to the concrete exceptions). Narrow each; `logger.debug` the swallow;
keep load-bearing fallbacks visible in debug logs.

### 7.4 Docs debt sweep

Verify/fix the June review's `preferred_api.md` structural issues (unclosed fence,
content after "Where To Look Next"); add `tests/test_docs_structure.py` (balanced
fences, resolving local links); gitignore generated PDFs (`myflopy_api_pamphlet.pdf`);
coordinate any `docs/` reorganization with the manual effort (don't do both).

---

## Phase 8 — Plotting consolidation into `myflopy/plot/` (D5; ~2–3 days)

**Position:** last structural phase — after 4.1 (observations split), Phases 5–6 (no new
API work targeting moving files), and 7.1 (imports already normalized). Only 7.4's docs
sweep and Phase 9 come after.

**Scope decision (explicit):** the `package_*` explorer family stays put
(`package_plotting.py`, `package_surface_water.py` are explorer components).
`datatypes/hover.py` also stays put — it is a spec engine (L0 leaf), not plotting.
This phase consolidates the *standalone* plotting modules only.

**Verified import graph (re-verify before moving — Phases 4–7 will have touched it):**

| Module (current) | Imported by |
|------------------|-------------|
| `utils/datatypes/choros.py` (`Choro`) | `grid/plotting.py`, `package_plotting.py`, `simulation/accessors.py`, `utils/inputs.py`, tests |
| `utils/datatypes/xsections.py` | `headsplus.py`, `package_plotting.py`, `simulation/accessors.py`, `project/model_group.py`, `model_results_diff.py`, tests |
| `mf6/contour_plotting.py` | `choros.py`, 1 test |
| `mf6/cross_section_plotting.py` | `layers.py`, `interactive_plotting.py`, `mf6/__init__.py` (exported), tests |
| `mf6/heads_plotting.py` | `headsplus.py` only |
| `mf6/mf2Dplots.py` | `heads_plotting.py` only |
| `mf6/budget_plotting.py` | `budget.py` only |
| `mf6/interactive_plotting.py` | `prt.py`, `simulation/base.py`, `mf6/__init__.py`, `viz.py`, top-level `__init__.py` (12 exported names), tests |
| `mf6/grid/plotting.py` | `grid/voronoi.py`, `pest/ies.py`, `viz.py`, 1 test |

**Target layout** — `src/myflopy/plot/`: `choropleth.py` (Choro), `xsections.py`,
`contours.py`, `cross_sections.py`, `heads.py` (+ fold the used parts of `mf2Dplots.py`,
attic the rest — completes 2.1's deferral), `budget.py` (reader `mf6/budget.py` stays),
`interactive.py`, `grid.py`, `__init__.py` re-exporting all public names. `viz.py`
remains the front door.

**Mechanic:** leaf-first move order (`contours` → `choropleth` → `xsections` → `grid` →
`budget` → `heads` → `cross_sections` → `interactive`); per move: `git mv`, fix internal
imports, old-path facade with 3.1's warning `__getattr__`, repoint internal importers
(so internal code never triggers its own deprecation warnings), fast suite, commit.
Update `_EXPORTS` targets (the 12 interactive names) and `mf6/__init__.py`; grep for
string references (pickles, docs, notebooks); update the layer map (plot/ = L3) and the
viz.py docstring. **Existing plotting tests must pass unedited** — if they break, the
move broke something. New `tests/test_plot_layout.py`: new path clean, old path
warns-and-works, top-level exports resolve.

---

## Phase 9 — Final wrap-up (~0.5 day)

1. Full-suite run recorded next to the Phase 0 baseline.
2. Final capability-doc alignment: `docs/myflopy_context.md`, `CLAUDE.md` (gap lists,
   module map, the now-real multi-physics claims), `docs/codebase_structure.md`.
3. Regenerate the API pamphlet if still the process (PDF gitignored per 7.4).
4. Walk Appendix C; link every box to its PR/commit.

---

## Sequencing, dependencies, effort

```
Phase 0 (baseline)      ── prerequisite for everything (BOTH repos clean)
Phase 1 (distribution)  ── independent; DO FIRST
Phase 2 (junk/hygiene)  ── independent; cheap; second
Phase 3 (deprecation)   ── before Phases 4 and 8
Phase 4 (splits/layering) ── 4.2 blocked on Phase 0
Phase 5 (package API)   ── independent of Phase 4;
                           order: 5.3A → 5.1 → 5.2 → 5.3B → 5.6 → 5.4 → 5.5 → 5.7
Phase 6 (model types)   ── 6.0/6.1 need 5.3 (buildable GWT models for fixtures);
                           6.2 after 6.1; 6.3 and 6.4 independent of 6.1/6.2;
                           6.1.5 (ConcTargets) composes with 4.1's observations split —
                           land whichever comes second on top of the first
Phase 7 (flopy4/robustness) ── after Phase 2; 7.1 MUST land before Phase 8;
                           6.0's new readers route through 7.1's compat module if 7.1
                           landed first (else repoint during 7.1)
Phase 8 (plotting)      ── after 4.1, Phase 5, Phase 6, and 7.1
Phase 9 (wrap-up)       ── last
```

| Item | Effort | Risk |
|------|--------|------|
| 1.1 figs vendoring | M | Low (CI exercises it) |
| 1.2 deps / 1.3 CI / 1.4 lint | S / M / M | Low–Med |
| 2.1–2.4 junk/hygiene | S–M | Low |
| 3.1–3.2 deprecation | S | Low |
| 4.1 observations split | M–L | Med |
| 4.2 model_group split | L | Med-high (in-flight work) |
| 4.3 surfaces / 4.4 layering | S–M / M | Low |
| 5.1 riv / 5.2 evt | M each | Low (proven pattern) |
| 5.3 gwt/gwe helpers + dispatch | M–L | Med (serialization round-trip) |
| 5.4 maw/hfb | L | Med-high (HFB geometry) |
| 5.5 dis/disu / 5.6 YAML | S–M / M | Low |
| 5.7 grid-lazy GIS | L–XL | High — forward-looking |
| 5.8 PEST build-side backlog | M–L | Med |
| 6.0 generic dependent-variable surface | M–L | Med (the keystone — get this right) |
| 6.1 GWT results tier | L | Med (fixture + budget terms) |
| 6.2 GWE results tier | M | Low if 6.0 held; STOP and fix 6.0 if not |
| 6.3 PRT maps + hover | M–L | Med (scatter-hover extension) |
| 6.4 PEST-IES viz/hover | M–L | Med |
| 7.1 flopy compat | M | Low |
| 7.2–7.4 logging/exceptions/docs | S–M | Low |
| 8 plotting consolidation | M–L | Med (9 modules, 12 exports move) |
| 9 wrap-up | S | Low |

S ≈ ≤half day, M ≈ 1–2 days, L ≈ 3–5 days, XL ≈ 1–2 weeks.

---

## Appendix A — Verification commands

```powershell
$env:PYTHONPATH = "src"
$py = "C:\Users\lukem\Python\mf-env\.venv\Scripts\python.exe"

& $py -m pytest -m "not slow" -q          # fast suite (after every change-set)
& $py -m pytest -q                        # full suite (end of each phase; ~14 min)

# import surface smoke
& $py -c "import myflopy; [getattr(myflopy, n) for n in myflopy.__all__]; print('ok')"

grep -rn "lukem" src/myflopy              # must be empty after Phase 2
grep -rn "X" src tests examples docs --include="*.py"   # orphan check before deleting X

# deferred-import ratchet baseline (~195 on 2026-07-06)
grep -rn "^\s\+from myflopy" src/myflopy --include="*.py" | wc -l

# hover/colorscale policy quick checks
& $py -m pytest tests/test_hover_spec.py tests/test_colorscale_policy.py -q
& $py -m nbstripout --dry-run (git ls-files "*.ipynb")
```

## Appendix B — Size snapshot (2026-07-06, for drift detection)

`src/myflopy` ~57.5k lines / ~90 files. Largest: `observations.py` 2,744 ·
`model_group.py` 2,572+ · `specs.py` 2,173 · `package_api.py` 1,823 ·
`package_plotting.py` 1,783+ · `grid/triangle.py` 1,742 · `package_surface_water.py`
1,740+ · `interactive_plotting.py` 1,344 · `simulation/base.py` 1,262 · `pest/ies.py`
1,235. Tests: 52 files, ~16k lines, 540 fast-passing (incl. `test_hover_spec.py`,
`test_hover_integration.py`, `test_colorscale_policy.py`). Tracked repo ~48 MB
(executed notebooks + mp3du exes — Phase 2 targets).

## Appendix C — Master acceptance checklist

- [ ] figs committed in its repo, then vendored; fresh venv `pip install -e .` →
      `import myflopy.viz` + map-sync post-scripts work without local figs (0, 1.1)
- [ ] `matplotlib` (+ `seaborn`, `pyyaml`) declared; import-surface smoke green (1.2, 5.6)
- [ ] CI green on ubuntu + windows (fast suite + ruff); vendored-figs path exercised (1.3, 1.4)
- [ ] `grep -rn "lukem" src/myflopy` empty; junk modules gone (2.1, 2.2)
- [ ] Notebooks stripped + pre-commit hook; no `.exe` under `src/`; `tools/mp3du/`
      resolution tested (2.3, 2.4)
- [ ] Deprecated aliases warn; GHB/DRN collision resolved; policy doc exists (3.x)
- [ ] `observations/` + `project/group/` packages; old paths work; no active module
      > ~1,800 lines without a written reason (4.1, 4.2)
- [ ] Layering test + deferred-import ratchet in CI (4.4)
- [ ] `mf.riv`, `mf.evt` with `()/.gpkg/.flopy` + registry + hover + earth/RdBu policy (5.1, 5.2)
- [ ] `mf.ic/oc/disv` dispatch on GWT/GWE; `mf.adv/dsp/mst/ssm/cnc/src/ist` +
      `mf.est/cnd/ctp/esl` exist; example 03 uses them (5.3)
- [ ] `mf.dis`/`mf.disu` passthroughs; fail-fast GridSpec tests untouched (5.5)
- [ ] `SimulationSpec.to_yaml/from_yaml` round-trips incl. exchanges/hooks (5.6)
- [ ] ONE generic dependent-variable surface; `model.conc` + `model.temp` with full
      grammar (map/xs/plot/mosaic/animate), `conc_hover`/`temp_hover`, earth + RdBu-diff
      colors, GWT/GWE budget terms, `GroupConc`/`GroupTemp` + `diff()` maps,
      `ConcTargets` wired into PEST, GWT + GWE test fixtures (6.0–6.2)
- [ ] PRT: `travel_time_map`/`endpoints_map` choropleths with hover + policy colors;
      pathline scatter hover styled; `mf.mip`/`mf.prp` factories (6.3)
- [ ] IES: `plot_field` hover + multiplier-diverging/absolute-earth colors;
      `field_uncertainty_map`; `field_mosaic` (synced prior/posterior); residual map;
      all IES figures on the viz front door (6.4)
- [ ] `_flopy_compat.py` is the only importer of flopy non-mf6 internals — including the
      6.0 readers (7.1)
- [ ] `myflopy` logger wired; no stray prints; exception swallows narrowed (7.2, 7.3)
- [ ] `myflopy/plot/` holds all standalone plotting; old paths warn-and-work;
      `_EXPORTS` repointed; `mf2Dplots` folded; existing plotting tests unedited (8)
- [ ] Final full suite recorded; capability docs match reality (9)

## Appendix D — (User-only) optional git-history purge (D3: out of scope for the agent)

For reference only — the agent must never do this. `git filter-repo` on `*.exe` + the
fat notebook paths, coordinate clones, force-push, re-clone. All content-level cleanup
in this plan works without it; this only shrinks `.git`.
