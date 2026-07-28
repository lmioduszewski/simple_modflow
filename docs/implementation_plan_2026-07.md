# myflopy Consolidation & Completion Plan

**Date:** 2026-07-14 (rev. 4 — full adversarial re-verification of every rev-3 claim
against the committed code. Reconciles the plan with the now-committed baseline (myflopy
`985eb9a`/`76c5f9f`/`962654d`, figs `bb1526b`); corrects claims that were wrong or had
drifted (two Phase 2 "orphan" rows, the `.flopy` pattern, `__compatibility__`, the
Phase 8 import graph, Appendix B numbers); reconciles the effort math and the 4.4 layer
map; adds the Linux environment, test-runtime budget, and release/tagging story; locks
new decisions D8–D12.)
**Prior:** rev. 3, 2026-07-06 (added Phase 6). **Source:** full structural/code review
of `src/myflopy` on branch `myflopy` (follow-up to `docs/refactor_review_report.md`,
2026-06-12), re-verified claim-by-claim on 2026-07-14. Line numbers are anchors and
WILL drift — always locate by symbol name, not line number.
**Audience:** an implementing agent working phase-by-phase. Read the whole
"Ground rules" section before touching anything.

---

## Completed since rev. 2 — COMMITTED as `985eb9a` + `76c5f9f` + `962654d` (2026-07-07) — do NOT redo

These subsystems exist, are committed, and are the *patterns to extend* in later phases.
Both repos are clean as of rev. 4 (2026-07-14); nothing below is in-flight anymore:

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
   SFR/LAK/combined-SW use `_exchange_colorscale(frame)`, which orients BLUE onto the
   gaining end per the declared reference frame (negative for `gwf` = SFR + list BCs,
   positive for `feature` = LAK / normalized `exchange_intensity`); per-package `q`
   results keep `"RdBu"` — and (b) ALL diff/compare maps (`"RdBu"` + `zmid=0` → negative
   RED, positive BLUE). **Everything else defaults to `'earth'`** (Plotly
   brown→cream→blue; the mounding-figure look). (The frame-naming replaced an earlier
   normalization; see `docs/compromises_and_deferrals.md` §45.)
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
4. **figs backend changes — COMMITTED in the separate figs repo as `bb1526b`
   (2026-07-07):** `Fig.add_post_script(js)` + threading into `show()`/`write_html()`;
   `write_html` signature fixed to plotly-compatible
   `write_html(self, file=None, config=None, ...)`; `_repr_mimebundle_` config merge;
   `add_logo=False` default. All four of Phase 1.1's snapshot preconditions are
   satisfied — **the vendoring gate is OPEN**. NOTE: the same commit restructured figs
   into `plotly/`/`mpl/`/`apps/` subpackages under `src/figs/`; Phase 1.1's module
   trace must be done against figs HEAD, not rev-3's description. Linux location of
   the figs checkout: `/home/lukem/python/Figures%20-%20Templates` (editable-installed
   in the `gw` env).
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
| D4 | Structured-grid stance | **DONE (5.5, 2026-07-23):** thin model-type-aware `mf.dis` / `mf.disu` passthroughs added (dispatch on model kind; PRT excluded; viz stays DISV-only). `GridSpec.structured` stays fail-fast. |
| D5 | Plotting consolidation | **Do it**, as its own late phase (Phase 8), after the god-module splits (Phase 4), API completion (Phases 5–6), and the `_flopy_compat` boundary (Phase 7.1). Details in Phase 8. |
| D6 | Hover defaults | DONE (see "Completed"). Heads default `layers="active+strip"`; sectioned/styled hover is the default everywhere; `custom_hover` stays as the raw escape hatch. |
| D7 | Colorscales | DONE (see "Completed"). Diverging only for signed q-like + diff maps; `'earth'` for everything else. New surfaces added by this plan MUST follow this policy. |
| D8 | `.flopy` on simple list BCs | First half DONE 2026-07-17 (branch `phase-5-package-api`): real `.flopy(...)` on `mf.chd/ghb/drn/wel` (thin passthrough shaped like `_RCHPackage.flopy`), GHB docstring aligned, pinned by `test_simple_list_bcs_have_a_real_flopy_escape_hatch`. Remaining: ship `mf.riv`/`mf.evt` with all three entry points (Phases 5.1–5.2). |
| D9 | GWT test fixture | **`transport=True` option on the small canonical config** (one-model-everywhere), NOT a promotion of example 03 into a second model family. GWE mirrors it in 6.2. Details in Phase 6.1.6. |
| D10 | `iheads.py` + `modflow/gwt/` removal | **Hard-delete now** (Phase 2.1) including their lazy-export map entries and smoke-test imports — they are public exports, not orphans. Justified pre-1.0 with zero tagged releases; note in the changelog. |
| D11 | Release/tagging | **Tag `v0.1.0` at the Phase 0 baseline commit**, then bump+tag at milestones (`v0.2.0` after Phase 1, `v0.3.0` after Phases 4–5). This starts 3.1's "≥ 2 tagged releases" deprecation clock and anchors the CI wheel build. |
| D12 | Alias visibility | Deprecated aliases/methods/classes are **hidden from autocompletion**: resolved only via `__getattr__` (module- or class-level), excluded from `__all__`, `dir()`, `TYPE_CHECKING` import blocks, and `.pyi` stubs. Runtime-compatible with a warning, invisible to IDEs and new users. Details in Phase 3.1. |

---

## Ground rules for the implementing agent

Non-negotiable project conventions. Violating them is a review failure even if tests pass.

1. **Environment (per-machine).** Windows: mf-env venv
   `C:\Users\lukem\Python\mf-env\.venv\Scripts\python.exe` with `PYTHONPATH=src`.
   Linux (the repo's current home, `/home/lukem/python/simple_modflow`): the
   provisioned env is `/home/lukem/python/envs/gw` (Python 3.14, flopy 3.10, figs
   editable). FULLY PROVISIONED as of 2026-07-16 (Phase 0/1): pytest/ruff/
   nbstripout/mypy/build, every optional extra (viz, viz3d, pest, parallel,
   xugrid), GDAL via system-bindings symlink, and **myflopy installed editable**
   — the old shadowing footgun (a stale non-editable copy in site-packages) is
   gone; `import myflopy` resolves to the repo with or without `PYTHONPATH=src`.
   Fast loop:
   `pytest -m "not slow"` (~60 s, ~540 tests — conftest auto-marks ~47 slow). Full
   suite: ~14 min. Fast suite after every change-set; FULL suite at the end of each
   phase.
2. **Tests are mandatory.** Real, runnable pytest tests covering behavior (not smoke
   tests) for every capability/refactor; run pytest to prove they pass before declaring
   done. Refactors need existing tests green AND a new test pinning the seam.
3. **Never change public import paths.** Splits/moves happen behind facades (the proven
   pattern: `package_explorer.py` fronting the `package_*` family). The lazy-export
   tables in `myflopy/__init__.py` and `modflow/mf6/__init__.py` must keep resolving.
4. **Grep both API layers before building.** Package-first (`package_api.py`, `specs.py`,
   `advanced.py`, `builders.py`, `geopackage.py`) sits on the OO engine
   (`modflow/mf6/*.py` + `simulation/`). **CANONICAL RULE (user, 2026-07-14): every
   change-set updates ALL documentation artifacts in the same pass** — every affected
   guide under `docs/` (including the capability map `docs/myflopy_context.md`),
   `CLAUDE.md`, tests, and any notebook/example demonstrating the touched surface.
   A missing doc/test/notebook update means the task is incomplete.
5. **Zero-based stress periods** in every myflopy-facing API.
6. **ASCII-safe console output** (Windows console; see commit `e61b973`).
7. **Do not reintroduce retired APIs** (legacy PEST `add_parameter`/`build_pst` path,
   retired PEST demo models).
8. **One phase per branch/PR.** A phase is done when all its acceptance criteria pass
   plus a full-suite run.
9. **Check the working tree first.** As of rev. 4 BOTH repos are clean: the rev-3
   in-flight work landed as myflopy `985eb9a`/`76c5f9f`/`962654d` and figs `bb1526b`
   (all 2026-07-07). Phase 0 reduces to environment provisioning, baseline recording,
   and tagging (D11). Still verify `git status` in both repos before starting any
   phase.
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

## Phase 0 — Baseline (prerequisite, ~0.5 day; step 1 ALREADY DONE)

**Steps:**
1. ~~Land or shelve ALL in-flight work~~ **DONE 2026-07-07**: the diff/group feature
   files and the rev-3 hover/colorscale/map-sync work landed as myflopy `985eb9a`
   (+ docs `76c5f9f`, docstrings `962654d`); the figs changes — including the
   `_repr_mimebundle_` merge and `add_logo=False` default that rev-3's Phase 0 text
   omitted from its figs change-set — landed as figs `bb1526b`. Verify both repos are
   still clean.
2. Provision the Linux environment (ground rule 1): install pytest, ruff, nbstripout
   into `/home/lukem/python/envs/gw`; confirm `PYTHONPATH=src` resolves the repo (not
   the shadowing site-packages copy).
3. Record the baseline fast/full suite results (and the full-suite wall time — the
   test-runtime budget in "Sequencing" tripwires at 2× this number) in
   `docs/phase_baselines.md`; later phases append their close-out runs there.
4. Note the starting commit and **tag it `v0.1.0`** (D11 — starts the deprecation
   clock; an annotated tag on the existing commit, no history rewrite).

**Acceptance:** clean `git status` in BOTH repos; a working local test loop;
recorded baseline test results + wall time; `v0.1.0` tag exists.

---

## Phase 1 — Distribution blockers (highest priority, ~3–6 days)

### 1.1 Vendor `figs` with external-first import (D1)

**Problem (re-verified 2026-07-14):** `src/myflopy/viz.py:39-40` hard-imports
`from figs import Fig, Subplot, Template, create_hover` and
`from figs.mpl import REPORT, Theme, get_mplfig`. `figs` resolves to the user's local
editable project (Linux: `/home/lukem/python/Figures%20-%20Templates`; Windows:
`C:\Users\lukem\Python\Projects\figs`; package root `src/figs/` on both) and is not
declared in `pyproject.toml`, not on PyPI. 21+ modules import `myflopy.viz`, so a pip
install of myflopy is broken for anyone but the author. Since rev. 3, figs `bb1526b`
restructured the package into `plotly/`/`mpl/`/`apps/` subpackages — trace against
figs HEAD.

**Implementation:**
1. **Snapshot the source** into `src/myflopy/_vendor/figs/` — ONLY the modules providing
   the symbols above plus their figs-internal deps, traced transitively **against figs
   HEAD** (the `bb1526b` restructure makes `figs.mpl` a subpackage sitting next to
   `plotly/`/`apps/` modules the trace must exclude); rewrite absolute `figs.` imports
   to relative so the subtree is self-contained. The snapshot gate is OPEN: all four
   preconditions (`add_post_script`, the `write_html(file=None, ...)` signature fix,
   `_repr_mimebundle_` config merge, `add_logo=False` default — all load-bearing for
   myflopy's map sync + inline scroll-zoom) are committed in figs `bb1526b`. Add
   `_vendor/__init__.py` and `_vendor/README.md` ("Vendored snapshot of figs @
   <figs-commit-hash>, synced <date>. Do not edit by hand; run
   `scripts/sync_vendored_figs.py`.") — **record the figs commit hash** so drift
   between the live and vendored copies is diagnosable.
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
   rewrites imports, stamps the sync date + figs commit hash. Idempotent. The figs
   checkout path is per-machine — take it from an env var or CLI arg, do NOT hardcode
   (Linux: `/home/lukem/python/Figures%20-%20Templates`, note the percent-encoded
   directory name).
4. **Packaging:** verify `myflopy._vendor.figs` lands in the wheel; add any figs data
   files to package-data if they exist.

**Tests:** `tests/test_vendored_figs.py` — vendored import + `Fig()` construction; a
subprocess with figs hidden proving `import myflopy.viz` works and map-sync post-scripts
still inject into `write_html` output; single-importer rule. CI (1.3) has no external
figs, so every CI run exercises the vendored path — the permanent drift alarm.

**Acceptance:** fresh venv, `pip install -e .`, NO figs on path → `import myflopy.viz`
works, `viz.Fig()` constructs, viz-touching fast tests pass.

### 1.2 Declare `matplotlib` (and audit the rest)

Add `"matplotlib>=3.8"` to `[project] dependencies` (COUNT CORRECTED rev. 4: 7 modules
import it top-level, 13 counting lazy function-level imports — rev-3's "10+ top-level"
overstated; the action is unchanged since it currently arrives only transitively via
figs). Audit the rest against a fresh-venv install: `seaborn` (lazy-imported by
`viz.mpl_axes`, `calibration.py`, and `pest/ies.py` — fails only when an mpl-backend
plot is requested, never at import time), `openpyxl` (no direct import anywhere — it is
pandas' implicit xlsx engine for the `pd.read_excel` calls in `recharge.py`/`ghb.py`/
`calibration.py`, so a fresh venv fails at runtime on xlsx reads; declare it or
document it as an extra), and document that FloPy's `Triangle` and the MODFLOW/PEST++
executables are runtime binaries, not pip deps.

**Acceptance:** fresh-venv install + import-surface smoke (1.3) green.

### 1.3 Add CI (GitHub Actions)

`.github/workflows/ci.yml`: matrix ubuntu-latest + windows-latest × Python 3.10/3.12
**+ 3.14** (the local dev interpreter is 3.14 — catch its breakage in CI, not locally);
`pip install -e .[dev]` → `pytest -m "not slow" -q` → (after 1.4) `ruff check src tests`.
CONTRACT CLARIFIED (2026-07-16, from the first real CI run): the fast suite never
RUNS models, but it does BUILD them — flopy validates exe paths at build time and
the Voronoi fixtures invoke `triangle` — so the fast job installs the `mf6` +
`triangle` binaries (`get-modflow --subset mf6,triangle`, seconds). First run will
surface optional-dep assumptions (pyvista/trame, pymetis/h5py, pyemu) — prefer
`pytest.importorskip` over fattening the CI install. Add the import-surface smoke step:
`python -c "import myflopy; [getattr(myflopy, n) for n in myflopy.__all__]"`.
Two more jobs (rev. 4): (a) a **wheel build** (`python -m build` + install the wheel in
a clean venv + import smoke) — this is also the automated check that `_vendor` lands in
the wheel (1.1 acceptance, previously verified by no automated step); (b) a
**scheduled slow-lane job** (nightly or weekly, ubuntu only) that installs MF6 binaries
(e.g. `modflow-devtools`' `get-modflow`) and runs `pytest -m slow` — without it the
slow e2e tests Phases 5–6 keep adding are never run automatically.

**Acceptance:** green CI on both OSes on a no-op PR; wheel job green; slow lane
scheduled (its first green run may trail the phase).

### 1.4 Lint + (scoped) type checking

`[tool.ruff]`: line-length 100, target py310, `extend-exclude = ["attic", "examples",
"scripts", "src/myflopy/_vendor"]`, `select = ["E","F","W","I","UP","B"]`. First pass:
safe autofixes only; triage the rest; never reformat in the same commit as logic changes.
Scoped mypy (or pyright basic) on the declarative core (`specs.py`, `package_api.py`,
`sources.py`, `workspace.py`, `advanced.py`, `builders.py`) — `py.typed` ships today with
nothing verifying it. Add tools to the `dev` extra.

**Acceptance:** `ruff check src tests` clean in CI; type check green on the scoped list.

---

## Phase 2 — Junk removal & repo hygiene (~2–3 days)

### 2.1 Delete/attic the orphaned modules

Re-verified 2026-07-14: 7 of the 9 rows are true orphans; **`iheads.py` and
`modflow/gwt/` are NOT** — both are public lazy exports with smoke-test importers
(D10 = hard-delete with the extra edit sites listed in their rows). Re-verify each
with grep before removing:

| File | Action | Notes |
|------|--------|-------|
| `modflow/utils/shiny_app.py` | delete | palmerpenguins Shiny demo |
| `modflow/mf6/pyqgis.py` | delete | mutates `sys.path` with hardcoded OSGeo4W path |
| `modflow/mf6/grasspyV.py` | attic | hardcoded `grass84.bat`; live GRASS path is `utils/contour_interp.py` (confirmed); ALSO update the markdown cell naming `grasspyV.SurfaceInterpFromShp` in `examples/mf6/notebooks/layer_management_workflow.ipynb` (docs-in-sync, ground rule 4) |
| `modflow/utils/openET.py` | attic | hardcoded Cumberland paths |
| `modflow/utils/iheads.py` | delete (D10 — NOT an orphan) | public lazy export: ALSO remove the `get_iheads` map entries in `modflow/__init__.py` and `modflow/utils/__init__.py`, and the import/assert in `tests/test_mf6_refactor_smoke.py`; superseded by `HeadsPlus` (it is a thin wrapper over it) |
| `modflow/mf6/mf3dplots.py` | attic | update the `viz.py` docstring that names it |
| `modflow/mf6/recovery_analysis.py` | attic | orphaned |
| `modflow/utils/gdal.py` | attic (found during Phase 1) | zero importers anywhere (re-verified 2026-07-16); broad excepts noted in 7.3 become moot |
| `modflow/mf6/gis_functions.py` | attic (found during Phase 1) | zero importers anywhere (re-verified 2026-07-16); top-level `osgeo` import |
| `modflow/calcs/cj_approximation.py` | attic | orphaned |
| `modflow/gwt/` (whole dir) | attic `gwt.py`; delete the dir — the only remainder is an empty `__init__.py` (D10 — NOT unimported) | scratch with personal paths; ALSO remove the `"gwt"` lazy submodule map entry in `modflow/__init__.py` and the `GWT` import in `tests/test_mf6_refactor_smoke.py`; the REAL transport home is Phases 5.3 + 6.1; afterwards add a permanent check to the layout test that nothing imports `myflopy.modflow.gwt` |

`prism_ppt.py` (1 importer): cascade-check. `mf2Dplots.py` (1 importer) is deferred to
Phase 8's fold — do not touch here.

### 2.2 Strip `__main__` scratch blocks with personal paths

Sites: `ghb.py`, `utils/raster.py`, `utils/prism_ppt.py`, `utils/surfaces.py`. Remove the
`C:\Users\lukem\...` blocks; convert anything worth keeping into a test or `examples/`
script. The other 7 `lukem` lines in `src/myflopy` live in files 2.1 already
deletes/attics (`openET.py`, `grasspyV.py`, `gwt/gwt.py`) — the acceptance grep passes
only if 2.1 lands too. Personal paths ALSO exist outside `src/` (rev. 4):
`examples/mf6/archive/cumberland_transient_observed_pest_first_slice.py`,
`examples/mf6/notebooks/canonical_model_template.ipynb` (loads `starting_heads.npy`
from a hardcoded `C:\Users\lukem` path), and two archive notebooks — fix or archive
these in the same pass. **Acceptance:** `grep -rn "lukem" src/myflopy` empty AND
`git grep -l "Users.lukem" -- examples` empty outside `examples/**/archive/`.

### 2.3 Strip all notebook outputs from git (D2)

`nbstripout` in dev extra + `.pre-commit-config.yaml` (+ fallback
`scripts/strip_notebooks.py`); one-time strip-only commit over ALL tracked notebooks
(incl. `notebooks/archive/`); gitignored `docs/_rendered/` + `scripts/render_notebooks.py`
for browsable executed copies; DELETE `examples/mf6/notebooks/starting_heads.csv`
(602 KB — re-verified 2026-07-14: nothing in the repo references it; the related
notebook loads `starting_heads.npy` from a personal path instead, see 2.2); gitignore
`scripts/scratch*`. Note: stripping shrinks the CHECKOUT (~31.5 MB of tracked
notebooks), not the ~170 MiB `.git` history (D3 keeps that out of scope).

**Acceptance:** `nbstripout --dry-run` clean over tracked notebooks; none over ~100 KB.

### 2.4 Move MP3DU executables out of `src/`

`src/myflopy/modflow/mp3du/*.exe` (10.7 MB) → untracked `tools/mp3du/` (tracked README
tells the user what goes there). Resolution order in `particles.py`: explicit arg →
`MYFLOPY_MP3DU_DIR` env var → repo `tools/mp3du/` → deprecated `_MODULE_DIR` fallback.
SEQUENCING FIX (rev. 4): the fallback's deprecation warning depends on 3.1, which lands
AFTER Phase 2 — use a plain `warnings.warn(..., DeprecationWarning)` here and convert
it to the 3.1 helper when 3.1 lands. gitignore `src/**/*.exe`. Unit-test the resolution
order (monkeypatched env, tmp dirs).

---

## Phase 3 — Deprecation mechanics & API-surface cleanup (~1 day)

> **DONE 2026-07-16** (branch `phase-3-deprecation`; see `docs/phase_baselines.md`
> and `docs/deprecation_policy.md`). Notes: the repurposed `__compatibility__`
> holds fully-qualified warned-alias names; the old second-tier meaning moved to
> `__engine__`; the mp3du legacy-PRT `__getattr__` and exe-path warning were
> absorbed into the helper alongside model_group's `_warn_deprecated`.

### 3.1 Deprecation helper + policy

`src/myflopy/_deprecation.py` with `warn_deprecated(old, new, *, since)` and
`deprecated_module_getattr(mapping, module)`; `docs/deprecation_policy.md` (short):
every compatibility name warns with its replacement; names live ≥ 2 tagged releases
(the clock is real now — D11 tags `v0.1.0` at Phase 0 and bumps at milestones);
`__compatibility__` in `myflopy/__init__.py` is authoritative. CORRECTIONS (rev. 4):
`__compatibility__` ALREADY EXISTS (since 2026-06-16) as an unwarned export-tier
marker over `_SECOND_TIER_EXPORTS` (21 names) — 3.1 REPURPOSES it as the warned-alias
registry rather than introducing it; and `project/model_group.py` already carries a
private `_warn_deprecated` (used 6× for `ModelGroup.rch/chd/drn/ghb/wel/uzf`) — absorb
it into the new module, do not leave two mechanisms.

**D12 — deprecated names are runtime-only, hidden from autocompletion.** Every warned
alias resolves ONLY through module-level `__getattr__` (PEP 562) — or class-level
`__getattr__` for deprecated methods/attributes (e.g. the `ModelGroup.rch/chd/...`
accessors) — never as a real module/class attribute. Deprecated names stay OUT of
`__all__`, OUT of the `TYPE_CHECKING` import blocks and any `.pyi` stubs (that is what
IDEs/static completers read), and module `__dir__()` must exclude them — implement
`__dir__` alongside `__getattr__` inside `deprecated_module_getattr` so hiding comes
free at every facade. Old code keeps working with a warning; new users never see the
old names suggested. Tests: `pytest.warns` per alias; per alias also assert
`name not in dir(module)` and `name not in __all__`; fast suite passes with
`-W error::DeprecationWarning` filtered to `myflopy.*`.

### 3.2 Resolve the GHB/DRN naming collision + FromVector consistency

Rename so `*FromVector` is the real class name in all four modules (`ghb.py` GHB→
GHBFromVector, `drn.py` DRN→DRNFromVector; chd/kflow already correct); keep warned
aliases via 3.1. SCOPE CORRECTION (rev. 4): smaller than rev. 3 implied — unwarned
aliases `GHBFromVector = GHB` / `DRNFromVector = DRN` already exist (since 2026-06-16)
and `mfsimbase.py` already imports the `*FromVector` names (its `__all__` lists both
spellings). Remaining work: flip which name is the real class, make the old
`GHB`/`DRN` names warn via 3.1 (hidden per D12 — the old names drop out of
`__all__`/`dir()`/TYPE_CHECKING; only `*FromVector` stays completion-visible),
repoint `modflow/mf6/__init__.py`'s `"GHB"`/`"DRN"` export entries, update docs.

---

## Phase 4 — Structural splits & import hygiene (~6–12 days)

> **4.1–4.4 DONE 2026-07-17** (branch `phase-4-splits`; see
> `docs/phase_baselines.md`). Deviations from the sketches below, decided
> against the real graph: `GroupOutputs` lives in `group/packages.py` (not
> `results.py`) to avoid a results↔lak cycle; the layering pin is the
> derived depth map (`docs/import_layering.md`, an output of
> `scripts/derive_import_layers.py`) rather than the L0–L6 sketch — the
> module graph was found ALREADY acyclic. §4.5 second-tier splits
> (triangle/package_plotting) deliberately deferred: not blocking, and
> `package_plotting` is better split with Phase 8 prep as §4.5 itself notes.

**Mechanic:** create the new package, move code in small mechanical commits, keep the old
module as a re-export facade, change zero user-facing paths, run the fast suite after
each move. No renames, no "improvements" during a split.

### 4.1 Split `modflow/mf6/observations.py` (3,014 lines, 13 classes)

Target package `modflow/mf6/observations/`: `_shared.py` (normalization/series helpers),
`heads.py` (`HeadTargets`, `BoundHeadTargets`, `BoundHeadTargetPlots`), `lake.py`,
`sfr.py` (stage+flow, bound), `drn.py`, `plots.py` (`BoundNamedSeriesTargetPlots`),
`registry.py` (`TargetRegistry`), `__init__.py` re-exporting everything (the package
shadows the old module path automatically). `heads_observations.py` is NOT part of this
split. Consumers (top-level `_EXPORTS`, pest package, `model.targets`) keep working —
verify with grep. Watch pickling references. NOTE (rev. 4): `pest/observations.py`
imports the private `_normalize_row_labels` from this module — `_shared.py` must keep
that cross-package import alive. Class inventory re-verified 2026-07-14: still exactly
13 (the +270 lines since rev. 3 are docstrings only). Tests: existing target tests
unchanged + `tests/test_observations_layout.py` (both import paths).

### 4.2 Split `project/model_group.py` (2,761 lines, 27 classes)

Target package `project/group/`: `spatial.py` (`_GroupSpatialView`, `GroupHeads`),
`budget.py`, `inputs.py` (+ `GroupPackageInputField`), `results.py`
(`GroupCellPackageResults` + namespaces + `GroupOutputs`), `sfr.py`, `lak.py` (incl.
connections), `uzf.py`, `surface_water.py`, `packages.py`, `core.py` (`ModelGroup`),
`__init__.py`. Keep `model_group.py` as facade. INVENTORY (rev. 4): 27 top-level
classes as of 2026-07-14 — `985eb9a` added `_GroupSpatialView` and
`GroupPackageInputField` and re-based views on `FieldMappable`/`LeafFieldSugar`
(rev-3's "28+" matched no committed tree). Re-verify again if any Phase 6 group work
(`GroupConc`/`GroupTemp`) lands first — the PREFERRED order is 4.2 before 6.1's group
work so the new classes get a proper home from day one. Tests: group/diff test files
green unchanged (incl. `test_model_group_diff_alignment.py`, new in `985eb9a`) +
layout test.

### 4.3 Consolidate the three `surfaces` modules

`myflopy/surfaces.py` (engine — KEEP), `grid/surfaces.py` (sampler — KEEP),
`modflow/utils/surfaces.py` (legacy; `InterpolatedSurface` still imported by
`layers.py`, `surface_data.py`, `xsections.py`, `test_layers.py`). Extract
`InterpolatedSurface` to a proper home (`grid/interpolated_surface.py`, or into the
engine ONLY if semantics align — read first); repoint the four sites; leave a warned
facade; attic the rest symbol-by-symbol.

### 4.4 Import layering + deferred-import ratchet

198 function-level `from myflopy...` imports (≈half of internal imports; 31 of them are
the intentional lazy-export machinery in `myflopy/__init__.py`) are circular-dependency
scar tissue centered on `simulation/base.py` (16 sites). LAYER-MAP CORRECTION (rev. 4):
the rev-3 map contradicted real edges — `headsplus` (L2) imports `heads_plotting`,
`budget.py` (L2) imports `budget_plotting`, `choros.py` imports `contour_plotting`, and
`utils/inputs.py` imports `choros` — so "L3 `*_plotting`" and a uniform-L0
`modflow/utils/*` are untenable; a literal AST test against that map fails immediately.
Fix: **derive the map FROM the actual AST graph first** (the layering doc is an output
of this step, not an input), then pin it. Corrected sketch to validate against the
graph before pinning:
- **L0** `_optional`, `_deprecation`, `_vendor`, `viz` (verified: imports externals
  only — rev-3's Phase 8 rows claiming viz imports plotting modules were phantom),
  `modflow/utils/*` leaves (`datatypes/hover.py` IS a leaf; `choros`/`xsections`/
  `inputs` are NOT — see L1)
- **L1** presentation leaves — the standalone plotting modules (today: `choros`,
  `xsections`, `contour_plotting`, `heads_plotting`+`mf2Dplots`, `budget_plotting`,
  `grid/plotting`; = `myflopy/plot/` after Phase 8) + `utils/inputs`; these import
  only L0 + externals
- **L2** `modflow/mf6/grid/*` + output readers: `headsplus`, `heads_observations`,
  `budget`, `budget_tables`
- **L3** explorers + observations + `interactive_plotting`
- **L4** engine facade: `simulation/*`, OO builders
- **L5** declarative: `specs`, `sources`, `geopackage`, `package_api`, `advanced`,
  `builders`, `layers`, `surfaces`, `grid_spec_resolver`
- **L6** lifecycle: `workspace`, `project/*`, `pest/*`, `parallel`, `prt`, `mp3du`
Runtime imports only point at same-or-lower layers; upward = TYPE_CHECKING, injection,
or an allowlisted deferred import — **deferred (function-level) imports are exempt from
the layering test** but counted by the ratchet. `tests/test_import_layering.py` (AST
walk vs the map, `_vendor` excluded) + a deferred-import **ratchet**
(`tests/deferred_import_allowlist.json`, baseline **198** as of 2026-07-14, counts may
only go DOWN; list the two lazy-export `__init__.py` files as intentional). Burn down
~30 easy conversions to prove the mechanism.

### 4.5 Second-tier splits (after 4.1/4.2)

`grid/triangle.py` (1,822 — +80 docstring lines since rev. 3) → move domain/feature
specs + cleanup/optimization orchestration into existing `grid/` collaborators behind
the `TriangleGrid` facade. GROWTH TRIGGER FIRED (rev. 4): `package_plotting.py` grew
771 → 1,914 lines (+148%, the `985eb9a` hover/colorscale/map-sync work) and is now the
largest file in the explorer family; `package_surface_water.py` is at 1,852. Rev-3's
"stay unless they keep growing" condition is met — plan a split of
`package_plotting.py` along its natural seams (payload builders / mosaic-animate /
colorscale helpers) behind the explorer facade, here or as Phase 8 prep. `specs.py` is
cohesive — leave below ~2,500 (currently 2,315).

### 4.6 Suite speedup — DONE 2026-07-17

Not originally a plan phase; recorded in `docs/phase_baselines.md` (branch
`phase-4.6-fast-suite`). Full suite 11m44s → 47.5s at `-n 10` via the
`CanonicalModelConfig.testing()` profile + session-scoped fixtures. Numbered here so
4.7 below is continuous.

### 4.7 Package declaration consolidation (the DRY keystone) — ~8–12 days

> **DO THIS BEFORE 5.3.** 5.3 adds ~11 GWT/GWE packages. Against **92 sites in
> `src/` that encode per-package knowledge**, that is a very large number of
> hand-edits at a **measured miss rate that already left 33 unintended gaps** in the
> current code. Consolidating first makes 5.3 dramatically cheaper and is the only
> way to stop paying this tax per package.

**The thesis: do not build a new abstraction — finish the one that exists.**
`modflow/mf6/package_registry.py` is already data-driven and it *works*: colorscales,
field names, budget terms and explorer wiring all flow from it, for every package,
with zero misses. The defect is that ~20 OTHER sites re-declare package knowledge
instead of reading from it. This phase extends the registry into a complete package
descriptor and deletes its competitors.

**Evidence this is not speculative (all from the 2026-07-17/18 riv+evt work):**
- `mf.riv`/`mf.evt` shipped with a documented four-piece checklist, the canonical
  docs rule, and TWO adversarial reviews — and **33 unintended gaps remain**
  (6–8 were caught and fixed post-review: `package_api.__all__`, `_DIFF_PACKAGES`,
  `_CELL_BUDGET_PACKAGES`, `_PACKAGE_SUFFIX_TO_TYPE`, plus the artifact set).
- **The rot is not new-package-specific.** `wel` is missing from the artifact
  subsystem entirely; `model.chd/drn/ghb/wel` convenience accessors don't exist
  while `model.rch/uzf/sfr/lak` do; `components.py:423` omits `drn/ghb/rch/wel`.
  These predate riv/evt by a long way and nobody noticed — the lists silently rot
  whether or not packages are being added.
- The **budget off-by-one** (riv/evt 1-based vs everything else 0-based) is a direct
  symptom: a hardcoded set in `budget_tables.py:63` omitted riv/evt, and
  `project/group/budget.py:33` "compensates" with an `if node.min() >= 1` heuristic
  that *guesses* what the set should have told it. Two bugs that cancel for riv/evt
  and leave `group.bud('drn'|'ghb')` **wrong today**.

**Measured site inventory (6-angle agent sweep, 2026-07-18, deduped by file:line).**
An earlier hand estimate of "~20" was wrong by 4.6×; the real figures are:

| Bucket | Count | Meaning |
|---|---|---|
| Sites in `src/` encoding per-package knowledge | **92** | the consolidation surface |
| — already complete (list every cell_stress pkg) | 34 | someone remembered; still boilerplate |
| — **true unintended gaps** | **33** | a package is silently absent — latent bugs |
| — feature gaps (PEST targets, obs families, canonical model) | 19 | real per-package work, NOT boilerplate |
| — intentional absences (deprecated / frozen legacy tiers) | 6 | **must NOT be extended** |
| Sites outside `src/` (tests, docs, examples) | 28 | also need updating per package |

Do not "fix" all 58 non-complete sites: 6 are intentional (deprecated-alias maps and
the frozen legacy OO tier correctly exclude riv/evt — D12), and 19 are genuine
feature decisions (e.g. `riv` as a PEST `parameterize` target is ledger entry 10, not
a missed edit). **The 33 are the bug surface.** Notable examples, all verified:
`budget_tables.py:63` (the live off-by-one), `components.py` ×6 (artifact capture /
restore — `wel` missing too), `simulation/base.py:736-756` (convenience accessors
exist for only 4 of 10 packages — `model.rch/uzf/sfr/lak` yes, `model.chd/drn/ghb/
wel/riv/evt` no), `model_diff.py:629-632` (static-typing properties),
`simulation/__init__.py` ×3 (export maps), `simplemodel.py` ×3 (per-package
dispatch), `choros.py:470`, `budget.py:196` + `budget_plotting.py:88` (hardcoded
`"drn"` defaults), `mfsimbase.py:49`, `boundaries.py:81`, `utils/inputs.py:19`,
`mp3du/particles.py` ×2, `parallel.py:229`.

Full machine-readable inventory: re-run the sweep workflow
(`docs/` has no copy — regenerate with the `package-declaration-inventory` workflow;
6 agents, ~13 min) or read `src/` for the symbols above.

#### 4.7.0 Make equivalence provable FIRST (do not skip)
> **DONE 2026-07-19.** `tests/api_snapshot.json` (1,162 lines) + regeneration script
> `scripts/derive_api_snapshot.py` + `tests/test_api_snapshot.py`. Covers all three
> dimensions: **23** helper signatures (every public method, parameter kind and
> default), **11** `*_spec` outputs from fixed literal inputs (name, builder identity,
> `requires`, and `option_summary` so data payloads are summarized not dumped), and
> **17** surfaces. Mutation-verified against the two failure modes 4.7.4 actually
> risks: silently changing a helper default, and silently dropping an option from a
> spec factory — both fail the snapshot. The script's `--check` mode runs inside the
> test, so a snapshot nobody regenerated fails CI. Intentional API changes: rerun the
> script and REVIEW THE DIFF, which is then the public-API change log for the PR.
> A `test_snapshot_covers_the_whole_public_package_surface` guard stops the net
> silently shrinking (a snapshot that stops covering something would still pass a
> plain equality check).

Golden-snapshot tests, landed before any refactor: (a) every `mf.<pkg>` public
signature (methods, params, defaults); (b) every `*_spec` output (`options` dict +
builder func/args) from fixed inputs; (c) the package set visible in each of the ~20
surfaces. Without this, "behavior-preserving" is a claim, not a fact. The snapshot
doubles as a table of every current inconsistency.

#### 4.7.1 Fix the live bugs first — do NOT refactor onto a broken baseline
> **Budget-basing pair DONE 2026-07-18.** Scope was larger than the review reported:
> `model.bud()` returned MF6's raw 1-based nodes for **chd, riv, wel and evt** (chd
> and wel long predate riv/evt), while `group.bud('drn'|'ghb')` came back one cell
> low from a double-subtraction. Fixed together: `_zero_base_budget_frame` now derives
> its package set from the registry's `cell_stress` entries instead of a hardcoded
> `{"drn","ghb","rch"}` (a preview of 4.7.3), and `GroupBudget._normalize_nodes` is
> deleted — `budget.df` is already normalized. New `tests/test_budget_node_basing.py`
> pins all three paths (legacy `bud()`, group, and the registry explorer) against
> ground truth taken from each package's own stress-period data; it failed 14 ways
> before the fix.

> **`default_input` pinning and `_PackageDiffNamespace` properties DONE 2026-07-18**
> — both mutation-verified (swapping riv `stage`->`cond` and dropping the riv diff
> property each now fail). The default_input test is table-driven over the registry's
> `cell_stress` set, so a new BC is covered automatically and an unpinned one fails
> loudly; the diff-property guard lives beside the `_DIFF_PACKAGES` invariant test
> because `py.typed` makes a `__getattr__`-only package type as `Any` downstream.

> **`rch_spec` maxbound DONE 2026-07-18** — removed (FloPy computes MAXBOUND at
> write time; written file verified byte-identical). All 7 list BCs are now covered
> against bare-list / `None`-period / empty-dict inputs.

> **Artifact subsystem DONE 2026-07-18 — 4.7.1 is now COMPLETE.** `wel`/`riv`/`evt`
> are supported; the type set and the capture/restore branches are registry-driven.
> Two latent bugs surfaced while fixing it: chd/drn CAPTURED `auxiliary` and dropped
> it on restore, and `boundnames` never round-tripped at all — it changes the record
> SHAPE (a trailing boundname field), so restoring without it made FloPy misparse
> every row and leave `NaN` cellids. That would have silently corrupted `wel`
> artifacts (canonical `wel` has `boundnames=True`) even if wel had simply been added
> to the old hardcoded set. `tests/test_package_artifact_round_trip.py` covers all 7
> list BCs end to end; mutation-verified.

#### 4.7.2 Grow the registry into a complete package descriptor — DONE 2026-07-18
Added to `PackageExplorerSpec`, for all 10 packages: `flopy_class` (a **string**,
so the registry keeps zero imports and the deferred-import ratchet is
undisturbed), `record_fields`, `gpkg_defaults`, `capabilities`
(`mover`/`auxiliary`/`boundnames`/`observations`/`edges_only`), `tiers`
(`diffable`/`connection_diffable`/`results_diffable`/`artifact_serializable`/
`artifact_apply_order`/`model_accessor`), `file_suffix`, `zero_base_budget_nodes`,
and a hand-written `blurb`. Budget text already existed via `ResultSpec`.

**Nothing consumes it yet** — that is 4.7.3. The value is entirely in
`tests/test_package_descriptor.py` (80 tests), where **every assertion reads the
ORIGINAL source rather than repeating the literal**: capabilities against
`inspect.signature` of the live FloPy class, `record_fields` against FloPy's own
**dfn**, `gpkg_defaults` against the live `GeoPackageSource` signatures, tiers
against the real tuples, node basing against `budget_tables`. A test that merely
restated the registry would pass while the descriptor was wrong. Mutation-verified
across five assertion families.

**One behaviour change, deliberate:** `rch` was removed from
`MvrResultDiff._MOVER_PACKAGES` (ledger entry 14, open since 5.1). MF6 has no RCH
mover, the lookup was try-wrapped and always found nothing, so this is
unobservable — and it turns a documented discrepancy into an equality invariant
against `capabilities.mover`. The API snapshot caught it, was regenerated, and the
diff is exactly that one line.

> **Note for 4.7.3:** the tier/suffix/node-basing tests are marked
> `RETIRE WITH 4.7.3` in the file. When a hardcoded list is deleted, its
> assertion becomes circular and must be deleted in the same commit — otherwise
> it silently stops testing anything.

#### 4.7.3 Delete the hardcoded lists — DONE 2026-07-18 (5 commits)
Five commits, lowest risk first, each verified identical to the literal it
replaced and each retiring its own now-circular scaffold test:

1. `_CELL_BUDGET_PACKAGES` + `_MOVER_PACKAGES` ← `tiers.results_diffable`,
   `capabilities.mover`
2. `_DIFF_PACKAGES` + `_CONNECTION_PACKAGES` ← `tiers.diffable`,
   `tiers.connection_diffable`
3. budget node basing ← `zero_base_budget_nodes` (was inferring from `kind`)
   — **narrowed 2026-07-27 to `node2` only.** The flag's `False` values for
   sfr/lak/uzf were measured wrong: in the *model* budget every record's `node`
   is a 1-based model cell, so `node` is now zero-based once, universally, in the
   shared record converter. Rename deferred (ledger 93).
4. `_PACKAGE_SUFFIX_TO_TYPE` ← `file_suffix` ∪ structural suffixes
5. artifact sets ← `tiers.artifact_serializable` / `artifact_apply_order` ∪
   `{ic, npf, mvr}`

**The union matters.** Three of these lists contain entries that are *not*
registry packages — `dis/disu/disv/ic/mvr/npf/oc/sto` in the suffix map,
`ic/npf/mvr` in the artifacts. They are structural packages with no per-package
behaviour to describe, so they stay explicit. Replacing any of those lists
wholesale would have silently dropped them and broken package discovery or
artifact restore.

**Each retirement got a replacement aimed at the NEW failure mode**, not a
weaker version of the old one — suffix collision between the two sources; a
serializable package with no apply order (falls back to 999, i.e. after `mvr`,
which references packages that must already exist); duplicate order numbers;
diff tiers overlapping. All mutation-verified.

Still hand-written, deliberately: the `ModelGroup` `GroupPackageInputs`
accessors and the `SimulationBase` package properties. Those are code objects
rather than data — they belong with 4.7.4's collapse. Their `model_accessor`
flag stays scaffold-tested against the real classes, which remains a genuine
cross-check rather than a circular one.

> The `test_diff_tier_covers_exactly_the_group_cell_bc_accessors` invariant
> inverted usefully: it used to compare two hand-typed lists, and now
> cross-checks the registry against the still-hardcoded group accessors. It is
> what catches the next package added to the descriptor but not the group.

#### 4.7.4 Collapse the per-package code — DONE 2026-07-18 (scoped down)
Shipped: one `_list_bc_spec(package, data, **opts)` behind the 7 named `*_spec`
wrappers, and one `GeoPackageSource._bc_from_features(package, **fields)` behind
the 7 resolvers. Every public name, signature and docstring is unchanged —
**`api_snapshot.json` came out byte-identical, which is the proof.**

**NOT shipped, and the plan text above was wrong to ask for it:**
`_ListBCPackage(descriptor)` replacing the 7 helper classes. The premise
("3,038 lines, mostly repeated structure") does not survive measurement. The
seven list-BC classes are 1,230 lines:

| | lines | share |
|---|---|---|
| docstrings | 747 | **61%** |
| signatures | 257 | 21% |
| executable bodies | 205 | 17% |

It is repeated *shape around hand-written prose*, not repeated logic. Each
docstring documents a genuinely different boundary (RIV's `rbot` cutover, DRN's
one-way flow, EVT's extinction depth) with worked record layouts. Collapsing
them would delete ~747 lines of prose in favour of generated text and replace
257 lines of explicit signatures with `**kwargs` — losing named parameters,
defaults and IDE completion on the public surface — to save ~150 lines of thin
delegation. Declined by the user 2026-07-18 (ledger 38).

**"Biggest LOC win" also did not hold.** The shipped collapse is +40 lines net:
the two shared bodies carry real docstrings explaining the invariants, and that
is worth more than the lines saved. The value delivered is the divergence class,
not the line count — "we never set `maxbound`" is now a property of one shared
implementation rather than something seven functions each have to remember,
which is exactly the bug that crashed `evt_spec`/`rch_spec`.

Two guards came free with the shared bodies, both mutation-verified: the record
field order is checked against `gpkg_defaults` (a silently reordered RIV record
builds a model that runs and is wrong), and `edges_only` is rejected for the
packages whose descriptor says they cannot support it.

#### 4.7.5 Extract `_ArealBuilder`; add `EVTBuilder` — DONE 2026-07-21
Lifted `RCHBuilder`'s cell selection + value broadcasting into a shared
`_ArealBuilder` base (`modflow/mf6/areal.py`, layer 2) with three subclass hooks
(`_row_for`, `_make_spec`, `_build_metadata`) plus `_prepared` for per-build
precompute; `RCHBuilder` is now a thin subclass with behavior unchanged. `EVTBuilder`
(`modflow/mf6/evapotranspiration.py`) adds `rate`/`depth` (reusing the base broadcast)
and the one genuinely new piece, `_resolve_surface`, which delegates to the shared
`SurfaceResolver`/`CellSurfaceOffset` engine (extracted from `GeoPackageSource` in
4.7.5a) — never reading `gdf_topbtm` directly. `mf.evt(context=, nper=, rate=, depth=)`
overloads the existing `__call__` exactly like `mf.rch`; `mf.EVTBuilder` is exported
in `__engine__`. Default surface is `model_top` (land surface), which resolves for
every column. Resolves ledger entry 6; deferred scope in ledger entry 51. Kept
minimal per user decision (three commits: 4.7.5a resolver, 4.7.5b base, 4.7.5c EVT).

#### 4.7.6 The payoff test — DONE 2026-07-18
`tests/test_package_descriptor_payoff.py` registers a synthetic descriptor in a
**fresh interpreter** (the derived lists are module-level constants, so it must
inject before any consumer imports) and reports where it lands.

**Measured answer: 9 surfaces automatic, 7 still hand-written.**

| automatic (no second edit) | still hand-written |
|---|---|
| `_DIFF_PACKAGES` | `advanced.<pkg>_spec` ¹ |
| `_CELL_BUDGET_PACKAGES` | `package_api.<pkg>` helper ¹ |
| `_MOVER_PACKAGES` | `GeoPackageSource.<pkg>` ¹ |
| `budget_tables.cell_based` | `ModelPackages.<pkg>` ² |
| `run_model.suffix_map` | `SimulationBase.<pkg>` ² |
| `components.LIST_BC` | `GroupPackages.<pkg>` ² |
| `components.SUPPORTED` | `_PackageDiffNamespace.<pkg>` ² |
| `components.APPLY_ORDER` | |
| `registry.budget_term` | |

¹ **deliberate** (ledger 38): named factories, helper classes and resolvers carry
per-package prose and explicit signatures. Generating them costs documentation
quality and IDE completion.
² **deferred**: namespace properties. The "stays duplicated" note below says
generating these degrades `diff.packages.riv` to `Any` downstream (the package
ships `py.typed`), so closing them needs a checked-in `.pyi` generated from the
registry — runtime generation alone is a regression.

So the plan's "assert it appears in EVERY surface" was never achievable as
written, and would have failed on day one. The test asserts the honest thing
instead, and is a **ratchet in both directions**: a new hand-written site fails
it, and so does closing one without moving the name out of `MANUAL_SURFACES` —
progress must be recorded, not absorbed. All three failure modes are
mutation-verified.

**Is 5.3 mechanical now?** Partly, and measurably so. The artifact subsystem —
where `wel` was missing for its entire life and `riv`/`evt` from birth — is now
automatic, as are the diff tiers, suffix map and budget basing. What remains per
new package is 3 documented declarations plus 4 namespace properties. The
category that produced the original six-site miss is closed; the category that
remains is visible, enumerated, and machine-checked.

#### 4.7.7 Complete the namespace accessors — DONE 2026-07-18
Scoped down from "generate the properties + emit a `.pyi`" after measuring.

**The measurement:** of the four namespace classes, `ModelPackages`,
`GroupPackages` and `_PackageDiffNamespace` were already **complete**. Only
`SimulationBase` was missing anything — six accessors, so `model.rch` worked
while `model.chd` raised `AttributeError`. Not a design decision; someone wrote
four and stopped.

**Why generation was rejected:** a `.pyi` stub **replaces the whole module** for
type checkers. Generating six properties would have meant hand-maintaining stubs
for every other public name in four large modules — a large, brittle artifact to
solve a one-class problem. Runtime generation without stubs would degrade
`diff.packages.riv` to `Any` for everyone downstream, since the package ships
`py.typed`.

**Shipped instead:** the six accessors written out (uniform
`self.package("<pkg>")`, same as the four that existed), all ten
`model_accessor` flags flipped to `True`, and the descriptor test upgraded from
"does the flag match reality" — which merely *documented* the gap — to
`test_every_namespace_exposes_every_registry_package`, which asserts all four
namespaces cover the whole registry. Mutation-verified: deleting one accessor
fails it. The API snapshot recorded exactly six additions.

So the per-package cost of these four surfaces is unchanged in lines but changed
in kind: forgetting one is now impossible to miss rather than invisible until a
user hits `AttributeError`.

#### What deliberately STAYS duplicated (the "good reason not to" list)
- **Static-typing surface.** Generated properties degrade `diff.packages.riv` to `Any`
  downstream (the package ships `py.typed`) — proven by review 2026-07-18. Generate at
  runtime BUT emit a **checked-in `.pyi` from the registry**, with a test that the stub
  matches. DRY must not cost type safety.
- **Per-package prose docstrings.** Template the *structure* from the descriptor; keep
  the prose as descriptor data. Do not trade doc quality for DRY.
- **Advanced packages (SFR/LAK/UZF/MVR).** Real data-model differences (reaches,
  connections, outlets). Share the spec wrapper + registry entry, NOT the record
  machinery.
- **Genuinely package-specific logic**, e.g. MVR's package-ordering validation.
- CHD having no conductance, WEL being signed, EVT being areal are **data** in the
  descriptor, not code branches.

### 4.8 View-layer conventions — PARTIALLY DONE 2026-07-18

Not originally a plan phase. Opened by a user report: a notebook cell hand-rolled
~20 lines of matplotlib to redraw `sfr.results.plot_long_profile`, losing
`scrollZoom`/pan/the house template and coloring gaining reaches **green** against
the documented blue (`docs/mf6io_reference.md`).

Root cause: the **spatial** half of the view grammar was documented
(`map/plot/xs` + `mosaic/animate` on every leaf) but **derived tables** — merged
frames that are not one mappable field — had no rule, so they accreted as
`foo()` + `plot_foo()` method pairs in five different verb spellings.

**Done:**
- `docs/view_layer_conventions.md` — the normative rule, both halves:
  `model.packages.<pkg>.<inputs|results>.<noun>.<verb>()`, every noun an object,
  every object answering the same verbs. Linked from `CLAUDE.md`.
- `SfrProfileView` (`sfr.results.profile`) as the reference implementation:
  `get`/`summary`/`plot`, `__call__(per=)` rebinding, and `signed_exchange=True`
  drawing per-reach bars colored off the **shared** diverging scale so discrete
  and continuous renderings cannot drift.
- `long_profile` / `plot_long_profile` → D12 warned aliases preserving their old
  return values exactly (the plot alias pins `signed_exchange=False`).
- Both notebooks collapsed to the one-liner; `canonical.py` migrated.
- `tests/test_sfr_profile_view.py` (15 tests; the sign-convention assertion is
  mutation-verified).

**Field-level verb sweep — DONE 2026-07-22 (ledger 18).** `sfr.results.q.profile`
and `sfr.results.stage.profile` (`SfrReachProfileView`) and `lak.results.q.budget`
(`LakBudgetView`) are now view nouns answering `get`/`summary`/`plot`; the loose
`plot_profile`/`plot_budget`/`budget_summary` spellings are D12 warned aliases
(hidden from completion, exact old returns). One deviation: the SFR field data
method was itself named `profile()`, so the noun replaces it and the frame comes
from `q.profile.get()` (no `profile()` df-alias; it had no callers). The LAK budget
bar stays matplotlib for now (ledger 52). Entry 17 (`profile` has no `map()`)
stands — spatial verbs live on the field explorers. These verbs are not in
`api_snapshot.json` (it tracks the `mf.<pkg>` facade), so the snapshot is unchanged;
callers migrated (2 canonical notebooks, API pamphlet, `test_colorscale_policy`),
tests added.

**Also landed 2026-07-18 (same user report, second defect).**
`targets.heads.calibration_plot()` drew fully formatted but empty axes. Two
causes, both fixed:

- **The plot could not say it was empty.** `CalibrationPlot.from_obs_vs_sim`
  dropped unpaired rows and rendered whatever survived — nothing. It now
  diagnoses *which side* is missing (all-NaN observed / all-NaN simulated /
  no overlapping rows / empty table), in the caller's own column names, as both
  a `UserWarning` and an on-figure annotation, on both backends. This is the
  default path behind **all five** target families, not just heads.
- **The canonical model had no measured heads.** Every head target carried
  `head: np.nan`, so `stats()` and `residuals()` were vacuous too. Targets now
  sample the regional water table (§`_synthetic_head_observations`), giving 14
  paired points and real statistics. Ledger entries 21–24 record the trade-offs
  (warn rather than raise; only the `obs_vs_sim` path diagnosed).
- **Follow-on, same day: the demo now reads as well calibrated.** The first fix
  left a ~20%-of-range systematic bias; the user's requirement is that the fast
  tour show a *well-calibrated* model. Closed by combining a fitted observation
  surface with a wider network — **RMSE 2.07 ft = 3.9% of head range, ME
  +0.05 ft** (testing, which the fast tour uses); 4.3% / −0.03 ft on validation.
  Note `regional` itself could NOT be refitted: it drives CHD/GHB/initial
  conditions/SFR/RIV, so fitting it to the solution is circular. Ledger entry
  25 records the method and the two findings that mattered.

Tests: `tests/test_calibration_plot_empty.py` (9),
`tests/test_canonical_head_observations.py` (11, including the
well-calibrated contract: RMSE < 5% of range and |ME| < 0.75 ft).

---

## Phase 5 — Package API completion (~8–15 days excluding 5.7, independent of Phase 4)

### 5.0 The pattern to replicate (read first)

> **SUPERSEDED IN SPIRIT BY 4.7 (2026-07-18).** The "four pieces" below is the
> minimum, not the whole job: riv/evt followed it exactly and still left 33
> unintended gaps, because **92 places** in `src/` encode package knowledge (see the
> measured inventory in 4.7). If 4.7 has landed,
> add the registry descriptor and let the surfaces derive — the list below is then
> historical. If 4.7 has NOT landed, follow it AND grep for every hardcoded package
> tuple/set first (4.7 lists the confirmed sites).

Every list-style boundary package = four pieces (DRN is the reference):
1. `package_api.py` helper class (`_DRNPackage`: `__call__` direct data / `.gpkg` GIS /
   `.flopy` escape hatch) + module-level singleton. CORRECTION + D8 (rev. 4): the
   simple list-BC helpers (`_CHDPackage`/`_GHBPackage`/`_DRNPackage`/`_WELPackage`)
   have NEVER had a `.flopy` method (only rch/uzf/sfr/lak/mvr do), even though the
   chd/drn/wel docstrings — and `advanced.py`'s — advertise one. **D8: as the first
   task of this phase, add a real `.flopy(...)` to all four existing helpers** (thin
   raw-FloPy passthrough, same shape as `_RCHPackage.flopy`), align the GHB docstring
   (currently "same as ()"), and ship `mf.riv`/`mf.evt` with all three entry points as
   Appendix C promises.
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
> **DONE 2026-07-17** (branch `phase-5-package-api`): all four pieces — `_RIVPackage`
> (`()`/`.gpkg` incl. `CellSurfaceOffset` stage/rbot/`.flopy`), `GeoPackageSource.riv`,
> `riv_spec`, registry entry (stage/cond/rbot `earth`, result q `RdBu`) + explorer/
> group wiring via the registry-generic paths + `mf.riv`/`riv_spec` lazy exports.
> CORRECTION (2026-07-18, post-review): the **diff** tier is NOT registry-generic —
> `_DIFF_PACKAGES` (`model_diff.py`) and `_CELL_BUDGET_PACKAGES`
> (`model_results_diff.py`) are hardcoded and were extended explicitly; the
> original DONE note overclaimed this. Any future list BC must extend both.
> Input hover rides the generic per-field `cell_input_hover` path like every other
> list BC (no bespoke extra-fields hover was needed). MF6 e2e:
> `test_geopackage_specs_run_through_project` runs RIV and asserts explorer
> table/colorscales on the completed run. **Canonical integration DONE 2026-07-18:**
> RIV is the un-routed outlet river below the lake (xn 0.92-0.99), carved out of
> the UZF footprint like the lake/stream cells; contract now 15 packages.

RIV records `(cellid, stage, cond, rbot)`. All four pieces; `.gpkg` with surface
references for stage/rbot. Input hover: `cell_input_hover("stage", extra_fields=("cond",
"rbot"))`; colorscale `earth`; result q joins the RdBu q family.

### 5.2 `mf.evt`
> **DONE 2026-07-17** (branch `phase-5-package-api`): `_EVTPackage` (`()`/`.gpkg`/
> `.flopy`, `nseg=1` default param, segments via `nseg` + record values),
> `GeoPackageSource.evt` (areal polygon→cell mapping on the shared `_boundary_data`
> plumbing — the same machinery `.rch` uses), `evt_spec` (maxbound inferred like
> `rch_spec`), registry entry (`surface/rate/depth` earth, q `RdBu`) + explorer/group
> wiring + lazy exports. MF6 e2e: `test_geopackage_specs_run_through_project`.
> **UPDATE 2026-07-21 (§4.7.5):** the earlier "no separate `EVTBuilder`" deviation
> is now reversed — `mf.evt(context=, nper=, rate=, depth=)` ships a file-less
> builder form via `EVTBuilder`, a thin subclass of the shared `_ArealBuilder`
> base extracted from `RCHBuilder`. Ledger entry 6 RESOLVED.
> **Canonical integration DONE 2026-07-18:** EVT sits on the valley walls beside
> RCH (standard recharge/ET pairing) and is DISJOINT from UZF, which does
> vadose-zone ET only (`simulate_et` auto-on via pet/extdp, no linear/square_gwet).
> Overlapping them would double-count one PET demand; the disjointness is pinned
> by `test_canonical_evt_and_uzf_footprints_stay_disjoint`.

EVT `(cellid, surface, rate, depth [, pxdp, petm, petm0])`; `nseg` via options. Reuse the
RCH areal-mapping prior art (`RCHBuilder`, `GeoPackageSource.rch`). Earth colorscale,
blue input hover.

### 5.3 GWT/GWE package-first helpers + model-type-aware core helpers

> **PREFER 4.7 FIRST.** This sub-phase adds ~11 packages against a measured **92
> per-package declaration sites** in `src/` (which today carry 33 unintended gaps).
> After 4.7 it is one descriptor per package plus genuinely new physics code.

**Part A — model-type dispatch — DONE 2026-07-22 (prerequisite for 5.5 and 6.1/6.2):**
`mf.ic`/`mf.oc`/`mf.disv` hardcoded `ModflowGwf*` classes so they could not serve
GWT/GWE models. Now they pass module-level dispatch builders (`build_ic`/`build_oc`/
`build_disv` in `builders.py`, beside `build_ims`) that resolve the FloPy class from the
built model's `model_type` (`"gwf6"`/`"gwt6"`/`"gwe6"`) via module-level class dicts —
serialization-safe (`_callable_ref` → `myflopy.builders:build_ic`; the class dicts are
never serialized). PRT raises a clear error (no `ModflowPrtic`); `npf`/`sto` stay GWF-only
(MF6 has no GWT/GWE variant). No registry or api_snapshot change (the `__call__`
signatures are unchanged). Tests in `test_package_api.py`. The original sketch:
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
use them (keep one raw-PackageSpec regression test; note the example exercises
adv/mst/ssm/ic/oc/dis but NOT dsp/cnc/src/ist — those factories need fresh test
fixtures rather than example retrofits). Registry `FieldSpec` entries for `cnc`/`ctp`
(blue input hover, earth). Do NOT resurrect `modflow/gwt/`.

### 5.4 `mf.maw` and `mf.hfb`
- **MAW:** high-level `mf.maw(wells=[...], context=)` computing connections from point +
  screen interval against layer elevations (new `MAWBuilder` in `modflow/mf6/maw.py`,
  UZFBuilder-shaped) + `.flopy` raw form. Slow e2e on the canonical grid.
- **HFB:** cell-*pair* records; on DISV pairs must share a face crossed by the barrier
  line. RESEARCH RESOLVED (rev. 4): no public shared-edge-with-geometry API exists —
  `shared_face_length` returns only a float and `find_adjacent_*` return ids — **but
  `build_disu_connectivity` (grid/connectivity.py) already computes the per-pair
  shared-face LineStrings internally and discards them.** Lift that segment logic into
  a public edge-lookup API (keyed by cell pair, intersectable with a barrier
  LineString; handle lines clipping cell corners or running coincident with faces),
  with unit tests, then `mf.hfb.gpkg(path, hydchr=...)`. Grid modules live at
  `modflow/mf6/grid/` (there is no top-level `myflopy/grid/`).

### 5.5 `mf.dis` / `mf.disu` passthroughs (D4) — DONE 2026-07-23
Thin factories mirroring `mf.disv`, using 5.3A dispatch (Gwf/Gwt/Gwe/Prt dis classes;
note `prt.py` already builds `ModflowPrtdis(v)` internally — don't duplicate).
`GridSpec.structured` stays fail-fast (`test_gridspec_fail_fast.py` untouched).
Docstrings state the Voronoi-first stance. Explorer/registry wiring only if the
explorers can render structured grids — verify, else document DISV-only.

**Delivered:** `_DIS_CLASSES`/`_DISU_CLASSES` + `build_dis`/`build_disu` in
`builders.py` (gwf6/gwt6/gwe6, keyed on `model_type`, serialization-safe); `mf.dis`
(nlay/nrow/ncol/delr/delc/top/botm/idomain) and `mf.disu` (nodes/nja/top/bot/area/
iac/ja/idomain + `**options`) factories in `package_api.py`, both mirroring `mf.disv`
and dispatching on the model kind. PRT is **excluded** exactly as `mf.disv` is —
`mf.prt` builds its own dis/disv internally, and MF6 has no `ModflowPrtdisu`; the
shared dispatch raises a clear `ValueError` for `prt6`. **No registry/explorer
wiring** — `dis`/`disu` follow `disv`/`ic`/`oc` (not in `package_registry.py`, no
results tier); the choropleth/xs/animate viz stays DISV/Voronoi-only, now stated in
the docstrings, `docs/package_api_reference.md`, and ledger. `test_gridspec_fail_fast.py`
untouched. Tests: `test_dis_dispatch_builds_structured_class_per_kind` + `dis`/`disu`
added to the round-trip/PRT-rejection test; snapshot regenerated (+2 signatures);
example 03 + `concerns.dis()` retrofitted to `mf.dis` (both run MF6 green).

### 5.6 YAML/TOML spec serialization — DONE 2026-07-23 (YAML; TOML deferred)
`SimulationSpec.to_yaml/from_yaml` (+ optional TOML via stdlib `tomllib` + `tomli-w`
extra) in a small `specs_io.py`; PyYAML (`safe_load`/`safe_dump` only) as a core dep;
`Project.add_simulation_from_yaml`. Normalize through `_json_value` first. Round-trip
tests incl. exchanges/hooks/GridSpec/refs + a hand-written minimal YAML that builds +
an example file in `examples/`.

**Delivered:** `specs_io.py` (`simulation_to_yaml`/`simulation_from_yaml`, PyYAML
safe mode, `yaml` imported lazily so `import myflopy` never pulls it in);
`SimulationSpec.to_yaml(path=None)` / `from_yaml(str|Path)` +
`Project.add_simulation_from_yaml`. PyYAML added as a core dep. Layer 0 (no
module-level myflopy deps) so `specs.py` imports it downward with **no** deferred-
import ratchet change. **Scoping found the round-trip was NOT actually complete for
list BCs** — chd/ghb/drn/riv/wel/rch/evt build via a non-importable
`functools.partial(_build_named, cls)`, so `to_dict` rejected them. **§5.6A** (done
first, user-approved) generalized `_callable_ref`/`_resolve_callable` to serialize
that partial (`{"partial": ref, "args": [{"$callable": class ref}], ...}`), so all 7
list BCs now round-trip and a real BC-carrying flow model survives YAML. **TOML
deferred** (no null type; `tomllib` is 3.11+ vs the `>=3.10` floor) — ledger 54; the
partial-encoding trade is ledger 55. Tests: `tests/test_specs_io.py` (string/file/
Project round-trip across BCs/refs/exchanges/hooks + a hand-written minimal YAML that
builds MF6) and `test_list_bc_partial_builders_round_trip` in `test_project_spec.py`.
Example: `examples/mf6/yaml_spec/` (`model.yaml` + `run.py`, runs MF6 green).

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
engine. (All four items re-verified still open 2026-07-14. Schedule: anytime after
5.3; its visualization counterpart is 6.4.)

### 5.9 Explicit non-goal: CSUB

`mf.csub` (compaction/subsidence) is deliberately OUT of this plan — rarely requested
for the Voronoi-first use cases and large (interbed data model). Recorded here because
`docs/myflopy_context.md` cites "plan §5.1–5.4" for the CSUB gap; this section makes
that pointer resolve. Revisit on user demand.

---

## Phase 6 — Full model-type integration: GWT, GWE, PRT, PEST-IES (~7–15 days)

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

### 6.0 Design keystone: ONE generic dependent-variable surface — DONE 2026-07-24

**Delivered (keystone + readable conc/temp; user-scoped increment).** The
`text=`/`.ucn` flopy claim verified true. Factored `HeadsPlus` into a
file-agnostic `DependentVariableFile(SpatialView, bf.HeadFile)` base
(`headsplus.py`) parameterized by `value_name`/`store_column`/`_output_suffix`/
`_binary_text`/`_choro_type`/default-hover; `HeadsPlus` keeps the head-only extras
(obs, mounding, legacy `choropleth`), and `ConcResults`/`TempResults` are ~15-line
subclasses. Model view is now **kind-aware**: `_initialize_from_built_run` stores a
kind-neutral flopy handle + `model_type`; `model.hds`/`.conc`/`.temp` are gated by
kind (a clear `AttributeError` on the wrong kind); `get_kstpkper`/`field_reader`
route the time axis to the right reader. The **full value-kind hook** (chosen over
the bypass) generalized `Choro`'s `type=='hds'`/`'elev'`/`model.hds` triple to a
`_depvar_reader`/`_value_column`/`_value_name` resolution (`_DEPVAR_READER_ATTR`),
so conc/temp maps get the same per-layer hover table as heads; `hover.py`
generalized the `LayerTable` header literal + `mark_dry` (head-only) and added
`conc_hover`/`temp_hover`. `'earth'` colorscale is already the `Choro` default.
Tests: `tests/test_gwt_gwe_results.py` (coupled GWF+GWT and GWF+GWE run MF6; read
`get/summary/array/map` + hover + colorscale + kind guard). **Deferred to
follow-ups** (per scope): GWT/GWE budget views (6.1.3), `GroupConc`/`GroupTemp`
group/diff (6.1.4), `ConcTargets`/`TempTargets` + transport calibration (6.1.5),
and the canonical `transport=True` fixture (6.1.6) — ledger 56.

#### Original design notes


Heads, concentration, and temperature are the same shape: a binary output file read per
(kstpkper, layer, cell), explored through the same grammar. **Do not clone
`HeadsPlus` twice.** Instead:
1. Read `HeadsPlus` (`modflow/mf6/headsplus.py`) and factor the file-agnostic core into
   a parameterized base (or make `HeadsPlus` itself take `text`/`label`/`unit`
   parameters). SCOPE NOTE (rev. 4): slightly bigger than text/label/unit — the single
   file-open is one `super().__init__` call, but `map()` couples to
   `model.cor(type='hds')` and the `head_hover` default, so the cor/SpatialView
   plumbing needs a value-kind hook too; other head-isms to parameterize:
   `value_name='head'`, the `.hds` default path, elev→head renames, the 1e29/1e30 dry
   sentinels. FloPy's `HeadFile` reads concentration/temperature binaries via
   `text="concentration"` / `"temperature"` (verify against the installed flopy — this
   sub-claim is still UNVERIFIED as of rev. 4, flopy was absent from the review
   machine; concentration output is declared by GWT OC `concentration_filerecord`,
   `.ucn` by convention).
2. Instantiate: `model.hds` (existing, unchanged public behavior), `model.conc` (GWT),
   `model.temp` (GWE). Each carries: `.array()`, `.wide()/.long()`, the grammar verbs
   `map/xs/plot/mosaic/animate`, contours, and the hover sugar
   (`hover_layers`/`hover_surfaces` — `LayerTable` is field-generic as claimed, with
   two head-isms to fix for conc/temp: the hardcoded `head`/`bot` column-header
   literal in the surfaces variant, and `_mark_dry` (head-below-bottom dagger), which
   must default off for non-head fields).
3. `SimulationBase`/`ModelView`/`Run.model()` must recognize GWT/GWE models: today
   `Run.model()` documents "the preferred myflopy **GWF** model view" — extend model-kind
   detection so a transport/energy model view exposes `conc`/`temp` (and NOT `hds`), and
   `run.model("transport")` works. Read `workspace.py:Run.model` + `ModelView` first.

### 6.1 GWT results tier (concentration)

> **Partially DONE 2026-07-24 (see §6.0 banner):** items 1 (hover) + 2 (colorscale)
> shipped, and `model.conc` is readable through the full grammar. **Still open:**
> 3 (budget view), 4 (`GroupConc`/diff), 5 (`ConcTargets`/PEST), 6 (canonical
> transport fixture) — ledger 56.
>
> **Item 3 prerequisites DONE 2026-07-27 — the budget path was BROKEN, not merely
> unplumbed.** Scoping item 3 against real coupled GWT and GWE runs found five
> defects, all now fixed (ledger 90–92, 95): transport models wrote a zero-byte
> `.cbc`; imeth=1 full-array terms returned an EMPTY table with no error; the two
> budget paths disagreed about node basing; units were hardcoded `ft³/d`; and
> `model.bud("ssm")` was unreachable. A fifth, `model.bud("sfr"|"lak"|"uzf")`
> returning node ids one cell high on GWF, was found and fixed with them.
> **Correction to item 3's premise below:** term discovery being "text-generic" was
> true but not sufficient — `build_budget_result_table` handled only ONE of MF6's
> two record shapes. Reader acquisition was also already fine (§6.0 made
> `_get_budget_reader` kind-neutral).
>
> **Item 3 DONE 2026-07-27 — `model.budget.<term>` ships.** A namespace on
> `SimulationBase` (so it lands on live *and* reopened views) whose terms are
> **discovered from the budget file**, each returning a `CellBudgetResultsExplorer`
> — the full spatial verb set, not the reduced one the package-level
> `<pkg>.budget.<term>` offers. Two scoping corrections worth keeping: the package
> namespaces hand-write every term as a literal (there is no normalizer to reuse),
> and ledger 52 is about `lak.results.q.budget`, a different object. Exposed on
> **every** model kind rather than gated to transport, since the plumbing is
> kind-neutral and `sto_ss`/`data_spdis` have no package accessor. See ledger 97
> for the four judgment calls and 98 for two package-tier defects left alone.

Building on 6.0:
1. **Hover:** `conc_hover(unit="mg/L")` factory — primary `conc`, title
   "Concentration", `layers="active+strip"` default, surfaces merge supported; unit
   configurable (model-dependent; do NOT hardcode — thread from a `units=` param with a
   documented default). Registry-style defaults on every conc map.
2. **Colorscale:** per D7, `'earth'` default for concentration; Δconc (diff maps) =
   `RdBu`, negative red / positive blue. Add to `test_colorscale_policy.py`.
3. **Budget:** GWT budget terms (STORAGE-AQUEOUS, SSM, DECAY, …) through the budget
   explorers. GOOD NEWS (rev. 4): term discovery is already text-generic
   (`get_unique_record_names`, `build_budget_result_table(budget_text=...)`); the real
   gap is reader acquisition — `_get_budget_reader` is literally
   `self.gwf.output.budget()` and `ModelView` is GWF-only (see 6.0.3) — so this is
   plumbing a transport budget reader/view, not new term tables; `model.bud`
   equivalent on the GWT view.
4. **Group/diff:** `GroupConc` mirroring `GroupHeads` (member maps + `diff()` Δconc maps
   via the existing compare payload machinery — it is value-column generic) +
   `compare_hover("conc", "diff", ...)`. Slot into Phase 4.2's `project/group/` layout
   if that split has landed (`group/conc.py`).
5. **Observations:** `ConcTargets` (measured vs simulated concentration at points)
   mirroring `HeadTargets` — home: `observations/conc.py` in the 4.1 layout (if 6.1
   lands before 4.1, add `conc.py` to 4.1's target file list so neither goes stale);
   wire into `cal.observe(...)`/`cal.forecast(...)` + `forward_run.py` post-processors
   so **transport calibration** works end-to-end (this is the PEST tie-in).
6. **Canonical transport fixture (D9, locked):** tests need a real GWT model. Add a
   `transport=True` option to a small canonical config (one-model-everywhere; do NOT
   promote example 03 into a second model family). Slow-marked, and **scoped session-
   or module-wide** so ONE MF6 run serves all 6.1 e2e assertions (see the test-runtime
   budget in "Sequencing"). Every 6.1 feature gets a slow e2e against it (map renders,
   hover template asserts, budget table, group diff, obs round-trip).

### 6.2 GWE results tier (temperature)

> **Partially DONE 2026-07-24 (see §6.0 banner):** `model.temp` reader + map +
> `temp_hover` + `'earth'` colorscale shipped alongside 6.1 (the keystone made it
> ~free). **Still open:** GWE budget view, `GroupTemp`/diff, `TempTargets` — ledger 56.
>
> **Budget prerequisites DONE 2026-07-27 alongside 6.1's** — every fix was verified
> against a real GWF+GWE run as well as GWF+GWT, because the two kinds do NOT share
> term names (GWE's storage term is `STORAGE-CELLBLK`, GWT's is `STORAGE-AQUEOUS`)
> and the hover unit differs (`E/T` vs `M/T`). See the §6.1 banner and ledger 90–92, 95.

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
1. **Cell-based derived maps — DONE 2026-07-25 (§6.3B):** `prt_maps.py` turns track
   records into three per-cell views exposed as **nouns** on `PRTRunResults` —
   `results.travel_time` / `.endpoints` / `.capture`, each answering
   `get/summary/plot/map/mosaic`. NAMING CORRECTION: the `*_map()` spellings drafted
   here are exactly the `verb_noun` pairs `docs/view_layer_conventions.md` forbids, so
   they became nouns (`travel_time.map(stat="median")`, `capture.map(by=…)`).
   Verified/decided while building:
   - **`icell` is a ONE-based whole-grid node number** — `cell=(icell-1)%ncpl`,
     `layer=(icell-1)//ncpl` (`== ilay-1`). Nothing in the repo did this conversion;
     `pathline_cell_table` is now its single home, pinned by a multi-layer unit.
   - **Real release groups** (locked user decision): `PRTReleasePoints.from_cells/
     from_points(group=…)` write PRP boundnames (+`merge()` renumbering `irpt`),
     `PRTProject` sets `boundnames` from the row width, and MF6 echoes the labels
     **uppercased** into the track CSV `name` column — the key `capture` groups by.
     Without groups `capture` raises naming `group=` and the `by="release_point"`
     fallback rather than silently degrading (ledger 58).
   - Time-integrated, so **no period axis**: `result_hover` gained `footer=` (default
     unchanged) and these maps pass `footer=()`; the series verb raises instead of
     failing inside a groupby on a missing `per` (`plot()` is a distribution figure —
     cumulative arrival curve / count bars). `layer=None` pools layers and
     **recomputes** the statistic over the pooled particles (never a stat of stats).
   - Unreached cells map to NaN, not 0 (ledger 57); `'earth'` for both counts and
     elapsed times (one-sided magnitudes, never signed) — `PRT_COLORSCALE`.
   - `capture.map()` returns the per-group **mosaic** (shared scale, synced views);
     `group=` returns one `Choro`, which is also what the facet composers use.
2. **Pathline map — DONE 2026-07-25 (§6.3C):** CORRECTION (rev. 4) was right that the
   existing pathline map (`plot_particle_pathlines`) is pure **matplotlib** (flopy
   PlotMapView), so this was a NEW plotly map, not a restyling. Built as a **fourth
   noun**, not a `pathline_map()` method (same naming rule as 6.3B):
   `results.pathlines` — a view whose `get()` is the normalized record table.
   Decided while building:
   - **`pathlines` was already taken** by the raw cached track frame (12 call sites
     incl. canonical notebook 03). Locked decision: **clean break** — the raw CSV moved
     to `results.track_records` (what FloPy/PyVista consume) and every site migrated;
     no proxy shim (ledger 61).
   - Maps are `go.Choroplethmap` over `vor.latlon`, i.e. **WGS84** — track `x`/`y` are
     model coordinates, so overlays need reprojection. New reusable grid helper
     `vor.points_to_latlon(x, y)` (`grid/geometry.py`, next to `get_gdf_latlon`);
     6.4B's residual points want the same one.
   - `viz.mosaic` copied **only** `panel.get_choropleth()`, silently dropping every
     overlay — contours and location markers have been missing from mosaics all along.
     Locked decision: fix it wholesale. `Choro` gained `add_overlay`/`overlay_traces`
     (one accessor now answers for contours, locs, and registered overlays alike) and
     mosaic copies them,
     guarding the coloraxis/`z` logic to the cell trace only (ledger 62).
   - `PALETTE` had **no qualitative sequence**; group colors were plotly's default
     colorway, so a group could change color between figures. Added
     `PALETTE.categorical` (Okabe-Ito, colorblind-safe) + memoized
     `viz.category_colors`, retrofitted into `travel_time.plot()`/`capture.plot()`.
   - One trace per particle (what makes per-particle hover work) is capped by
     `max_particles=250`, sampled **stratified by release group** so a cap cannot
     drop a whole capture zone, and announced in the title + a warning (ledger 63).
   - `pathline_hover()` renders **per trajectory vertex**, not per cell — the first
     spec to do so; `HoverContext(ncpl=<vertices>)` (the assembler's `ncpl` is really
     a row count). `xs()` raises (no per-cell field to slice); `plot()` is elevation
     vs travel time.
   - A 3-lens adversarial review (correctness / docs-grammar / mutation-tested tests)
     found and fixed, before commit: a **`KeyError` on any run mixing named and
     un-named release points** (MF6 writes a blank `name`; the label has to be
     resolved before the palette is keyed); `category_colors` handing two visible
     groups the **same color** once the 7-color palette wrapped (now per-call
     collision avoidance, ledger 65); a `color="particle"` map registering one memo
     entry per particle (`memoize=False`); `mosaic(base=<Choro>)` drawing every group
     onto **one shared map**; and four selectors silently dropped rather than raising
     (`**base_kwargs` on `backend="mpl"`, `per`/`layer` on a borrowed base, blank-base
     option clashes, late `color=` validation). Mutation testing also showed the first
     stratification and colour-memoization tests passed against broken code — both
     were rebuilt around fixtures where the mutation actually changes the result.
3. **Build-side parity — DONE 2026-07-24 (§6.3A):** `mf.mip`/`mf.prp` are plain
   5.3B-style factories, NOT 5.3A dispatch (correction: MIP/PRP exist only on PRT —
   no cross-kind ambiguity). Full declarability additionally required (scoping
   verified by a live spec-declared GWF+PRT MF6 run): **`mf.ems` + `build_ems`**
   (EMS must be *registered* via `register_solution_package`; a bare EMS PackageSpec
   is silently dropped and MF6 aborts "Explicit models require EMS6"), **prt6
   entries** in the dis/disv/oc dispatch tables (`ModflowPrtdis/Prtdisv/Prtoc`;
   ic/disu correctly stay rejected — no `ModflowPrtic`/`ModflowPrtdisu`), and a
   **kind-aware `mf.simulation` solver default** (EMS for prt models — the old
   one-IMS-per-model default was a runtime landmine). CORRECTION: the spec path
   needs **no FMI at all** — the GWF-PRT exchange passes flows directly; FMI/grid
   copying is exclusively `PRTProject`'s post-hoc separate-simulation concern, so
   there was nothing to duplicate. `mf.prp` derives `nreleasepts`, defaults
   `perioddata={0: ["FIRST"]}`, auto-enables `boundnames` from the row width, and
   pins `pname` (flopy otherwise numbers instances "prp_0"). Misleading docstrings
   fixed (mf.prt/mf.ims/mf.simulation + stale "PRT builds its own dis" comments).
   Pinned by `test_prt_model_fully_declarable_and_runs` (runs MF6; boundnames echo
   uppercased into the track CSV `name` column). `PRTProject` stays the sanctioned
   post-hoc runtime path.
4. **3-D scene** stays as-is — CORRECTION (rev. 4): it is **PyVista/trame** (not raw
   plotly); either way it remains a documented viz.py exception. MP3DU untouched
   except Phase 2.4.
5. Group story: comparing PRT runs across models = diff of travel-time maps (mechanical
   via the compare payload once 6.3.1 exists) — mark forward-looking, implement only the
   single-model maps + hover now.

### 6.4 PEST-IES visualization & preferred-API integration — DONE 2026-07-26

> **All seven items landed across 6.4A (2026-07-25), 6.4B and 6.4C (2026-07-26).**
> Three of the four planned figure names changed and one planned method was not built
> at all, each for a reason established against the code rather than assumed — see the
> per-item notes below and ledger 67–86. The recurring lesson: **the review layer was
> further along than the plan assumed**, so every item here started with "verify what
> already ships" and two of them shrank to nothing on contact. 6.4C's own scoping found
> the same thing a third time (the capture-field fixture it was told to build already
> existed), plus two wrong instructions in the handoff file it inherited.

The review layer exists and is good; bring it up to the house standards and close the
workshop gaps (memory: spatial parameter-field maps were the biggest gap; DBTL framing +
input-capture-as-observations is the driving workflow):
1. **Parameter-field maps get hover + policy colors. — DONE 2026-07-25 (6.4A).**
   `IesResults.plot_field(target)` now carries `parameter_field_hover()` (primary = the
   plotted stat, block = prior mean / posterior mean / posterior sd) and a per-**stat**
   colorscale policy `_field_map_policy`.
   > Two corrections to what this item originally said, both established against the
   > code rather than assumed. **(a)** The captured field is *always absolute resolved
   > K*, never a multiplier — pyEMU rewrites `org_array × multipliers` before each run
   > and we capture the rewritten array — so the policy keys on the STAT, not on
   > "multiplier vs absolute": `change` → diverging red-white-blue, log-centered at 1
   > with symmetric limits (red = reduced); `mean`/`base` → `'earth'` + logscale;
   > `std` → `'earth'`, linear (a spread is legitimately 0 where the ensemble
   > collapsed). **(b)** Iteration and realization count *cannot* go in the hover
   > footer — `HoverSpec._render_footer` recognizes only period/step/date/area/model
   > and silently drops anything else — so they ride the figure title, visible on the
   > matplotlib backend too. Ledger 67–73.
   Pinned in `test_colorscale_policy.py` (including a call-site pin proving both
   backends put red on the same end) and `test_hover_spec.py`.
2. **Uncertainty maps: `plot_field(target, stat="reduction")`. — DONE 2026-07-26 (6.4B).**
   > This item originally asked for a `field_uncertainty_map(target, which="std"|
   > "reduction")` method. Scoping found **half of it already shipped in 6.4A**:
   > `plot_field(stat="std", which="prior"|"posterior")` is exactly the posterior-sd
   > map, so a new method would have been a second spelling of an existing figure —
   > the `foo()`/`plot_foo()` duplication `view_layer_conventions.md` forbids — and
   > `which="std"` would have given `which` a second meaning in a class where it means
   > prior-vs-posterior everywhere else. Only the variance reduction was new, and it is
   > one derived column (`field()` already returned `prior_std`). So it landed as
   > `stat="reduction"`, not a method. Ledger 75.
   > The scale is **anchored to [0, 1]** — an absolute frame, unlike a spread in model
   > units — but falls back to the data range when any cell is negative, so a posterior
   > that *grew* is not clipped to the bottom color where it would read as "the data
   > said nothing".
3. **Prior-vs-posterior mosaics: `IesResults.plot_field_mosaic(...)`. — DONE 2026-07-26 (6.4B).**
   Panels come from a new `_field_choro` seam (`viz.mosaic` needs the `Choro`; `plot_field`
   returns a rendered figure). `stat` is restricted to `mean`/`std` — the only stats with
   a separate prior and posterior form — which is also what makes ledger 71's `diff=True`
   trap unreachable rather than merely documented. **`viz.mosaic` gained a `colorbar=`
   passthrough** accepting a dict or a `(cmin, cmax) -> dict` callable, because the pooled
   limits are known only inside `mosaic`; without it a log-scaled `mean` mosaic read
   `−3 … 2`. Ledger 71 half-retired.
4. **Residual map: `obs_residuals()` + `plot_obs_residuals()`. — DONE 2026-07-26 (6.4B).**
   Heads draw as points, DRN zones color their cells, **both on one symmetric scale** so a
   point and the cell under it mean the same thing at the same color. Shipped **without**
   the planned `map=True` flag: there is no `map=False` form, and a flag with one honest
   value is a worse API than a named method (ledger 77).
   > Two things scoping settled. **(a)** `pst.try_parse_name_metadata()`'s `usecol`
   > column truncates at the first underscore (`obs_00` → `obs`), so it is useless as a
   > join key for real location names — the prefix and location are parsed out of the
   > whole obs name instead. **(b)** `plot_mpl` draws in **model coordinates**, so the
   > static backend scatters raw x/y and only the Plotly path needs `points_to_latlon`.
   > That is why this is not plotly-only. Lake/SFR targets record only a lake or reach
   > number and are deferred — ledger 76.
5. **Ensure every IES plot uses the viz front door. — DONE 2026-07-26 (6.4C).**
   The audit found it **already clean**: every figure-producing entry point on
   `IesResults` (plus `IesForecast.plot` and the private `_field_choro` seam) routes
   through `viz.Fig`/`viz.subplots`/`viz.mpl_axes`/`viz.mosaic` — or a `Choro`, which is
   built on them — on both backends; `ies.py` constructs no raw `go.Figure`/
   `plt.subplots`; and the other eleven `pest/` modules build no figures at all. Two
   defects, both fixed:
   the last bare `color="0.6"` greys became `_MPL_PRIOR` — they draw the **prior**
   ensemble, whose Plotly twins already used `_PRIOR_COLOR`, so the handoff's instruction
   to use `mpl_ensemble` ("0.5") would have changed the render and inverted the semantics
   `viz.PALETTE` documents — and two `with sns.axes_style("whitegrid")` wrappers around
   `viz.mpl_axes`, which already does exactly that, were removed. `edgecolor="black"`
   stays a literal (ledger 85); the pie charts' default category colors are deferred
   (ledger 84).
6. **Preferred-API surface check. — DONE 2026-07-26 (6.4C).**
   `model.pest(name, ...)` is now documented as the front door in the capability map,
   the package API reference and the `PestProject` docstring, with direct construction
   called out as advanced (its docstring previously called *itself* "the single front
   door" and had an "or directly::" clause whose example was `model.pest(...)` again).
   > Scoping turned up a real defect behind the docs: **`run.pest_runs` attached no
   > model**, so `run.pest_runs[i].review()` refused every spatial map, while both PEST
   > notebooks told the reader it was interchangeable with
   > `model.pest_runs`. Fixed in the code rather than the docs — via a lazy
   > `model_factory`, since listing calibrations must not load a MF6 simulation — which
   > makes the notebooks' existing claim true without editing a canonical notebook.
   > Ledger 82. Moving `mf.PestProject` out of the preferred export tier is deferred
   > (ledger 86): that is a public-API change this item does not ask for.
   The parameterization backlog itself is Phase 5.8.
7. **Tests. — DONE 2026-07-26 (6.4B + 6.4C).**
   > The synthetic fixture this item asks for was **already built in 6.4B**
   > (`_stub_ies_with_capture`, plus `_stub_ies_results` for the residual path) — the
   > handoff's claim that "only the slow e2e proves the field maps" was stale on
   > arrival. The real deficit was that the stub was too *degenerate* to see the
   > behavior: prior and posterior shared one `means` tuple, so `change` was identically
   > 1.0, both mosaic panels were pixel-identical, and **a swapped `which=` would have
   > passed the whole suite**. 6.4C parameterized it (separate prior/posterior means,
   > a cell subset, the realization index, the capture family, the model) and added
   > ten fast field tests — including the NaN padding onto uncaptured cells, which ran
   > in no test at all because both the stub and the e2e captured every cell.
   Also added: six forecast-less tests over a real `pyemu.Pst` (ledger 73), two for
   `run.pest_runs` (ledger 82), and four backend assertions the audit found missing.

### 6.5 Capability-map + manual updates

After each of 6.1–6.4: update `docs/myflopy_context.md`, `CLAUDE.md` (the "multi-physics"
claim becomes real; remove stale gap lines), and add one canonical notebook or example
per surface (a transport-results walkthrough; a PRT travel-time map example; an IES
field-map review section in the existing PEST notebooks).

---

## Phase 7 — FloPy 4 readiness & robustness (~2–4 days)

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
BEFORE Phase 8** so plotting files move exactly once with clean imports. (Relative to
Phase 6: either order works — 6.0's readers route through the compat module if 7.1
landed first, else get repointed during 7.1; the sequencing block is authoritative.)

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

RESOLVED (rev. 4): the June review's `preferred_api.md` issues (unclosed fence,
content after "Where To Look Next") were already fixed by the 2026-06-16 rewrite —
nothing to do there. Remaining: add `tests/test_docs_structure.py` (balanced fences,
resolving local links); gitignore generated PDFs (`myflopy_api_pamphlet.pdf` is
tracked today); coordinate any `docs/` reorganization with the manual effort (don't do
both). 7.4 ships inside Phase 7's PR — rev-3's "7.4 comes after Phase 8" ordering is
RETRACTED (it split one phase across a PR boundary, violating ground rule 8, and
nothing in 7.4 depends on Phase 8).

---

## Phase 8 — Plotting consolidation into `myflopy/plot/` (D5; ~2–5 days)

**Position:** last structural phase — after 4.1 (observations split), 5.1–5.6 + the
*landed* parts of Phase 6 (no new API work targeting moving files), and 7.1 (imports
already normalized). CLARIFIED (rev. 4): **5.7 (grid-lazy GIS, forward-looking) does
NOT gate Phase 8** — it targets spec/geopackage code, not plotting files; and D5's
"after the god-module splits" means 4.1 specifically (4.2–4.5 don't touch the moving
files, though they land earlier anyway in the stated order). Only Phase 9 comes after.

**Scope decision (explicit):** the `package_*` explorer family stays put
(`package_plotting.py`, `package_surface_water.py` are explorer components).
`datatypes/hover.py` also stays put — it is a spec engine (L0 leaf), not plotting.
This phase consolidates the *standalone* plotting modules only.

**Verified import graph (CORRECTED 2026-07-14; re-verify before moving — Phases 4–7
will have touched it):**

| Module (current) | Imported by |
|------------------|-------------|
| `utils/datatypes/choros.py` (`Choro`) | `grid/plotting.py`, `simulation/accessors.py`, `utils/inputs.py`, 6 test files — NOT `package_plotting.py` (rev-3's edge was phantom; only a local variable name) |
| `utils/datatypes/xsections.py` | `headsplus.py`, `package_plotting.py`, `simulation/accessors.py`, `project/model_group.py`, `model_results_diff.py`, tests (+221 lines in `985eb9a`) |
| `mf6/contour_plotting.py` | `choros.py`, 1 test |
| `mf6/cross_section_plotting.py` | `layers.py`, `interactive_plotting.py`, `mf6/__init__.py` (exported), tests |
| `mf6/heads_plotting.py` | `headsplus.py` only |
| `mf6/mf2Dplots.py` | `heads_plotting.py` only |
| `mf6/budget_plotting.py` | `budget.py` only |
| `mf6/interactive_plotting.py` | `prt.py`, `simulation/base.py`, `mf6/__init__.py`, top-level `__init__.py` (**13** exported names — rev-3's "12" was never right), tests — NOT `viz.py` (docstring mention only) |
| `mf6/grid/plotting.py` | `grid/voronoi.py`, **`grid/__init__.py` (re-export — rev. 3 missed it)**, `pest/ies.py`, 1 test, + the `canonical_02`/`canonical_06` notebooks import `build_choropleth` — NOT `viz.py` (docstring mention only) |

**Target layout** — `src/myflopy/plot/`: `choropleth.py` (Choro), `xsections.py`,
`contours.py`, `cross_sections.py`, `heads.py` (+ fold the used parts of `mf2Dplots.py`,
attic the rest — completes 2.1's deferral), `budget.py` (reader `mf6/budget.py` stays),
`interactive.py`, `grid.py`, `__init__.py` re-exporting all public names. `viz.py`
remains the front door.

**Mechanic:** leaf-first move order (`contours` → `choropleth` → `xsections` → `grid` →
`budget` → `heads` → `cross_sections` → `interactive`); per move: `git mv`, fix internal
imports, old-path facade with 3.1's warning `__getattr__` (per D12: the old paths
warn-and-work at runtime but are absent from `__all__`/`dir()`/TYPE_CHECKING — only
the new `myflopy/plot/` names are completion-visible), repoint internal importers
(so internal code never triggers its own deprecation warnings), fast suite, commit.
Update `_EXPORTS` targets (the 13 interactive names) and `mf6/__init__.py` +
`grid/__init__.py`; grep for string references (pickles, docs, notebooks — incl. the
`build_choropleth` imports in the canonical notebooks); update the layer map (plot/ =
L1 presentation leaves per the corrected 4.4 map) and the viz.py docstring. **Existing
plotting tests must pass unedited** — if they break, the move broke something
(`985eb9a` added five more plotting-test importers this rule now covers:
`test_hover_spec`, `test_colorscale_policy`, `test_view_grammar_composers`,
`test_group_map_api`, `test_hover_integration`). Note `choros.py` now lazily imports
`datatypes/hover.py` inside `_build_hover_context` — hover stays put (L0 leaf), the
edge just needs to survive the move. New `tests/test_plot_layout.py`: new path clean,
old path warns-and-works and is hidden per D12 (not in `dir()`/`__all__`), top-level
exports resolve.

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
Phase 0 (baseline)      ── step 1 DONE (both repos committed + clean); remaining:
                           env provisioning, baseline record + wall time, v0.1.0 tag
Phase 1 (distribution)  ── independent; DO FIRST (1.1's figs gate is already open)
Phase 2 (junk/hygiene)  ── independent; cheap; second. 2.4's deprecation warning
                           depends on 3.1 — land a plain warnings.warn, convert in 3.1
Phase 3 (deprecation)   ── before Phases 4 and 8; tagging (D11) makes the policy real
Phase 4 (splits/layering) ── 4.2 PREFERRED before 6.1's group work (else re-derive
                           the class inventory); 4.5 after 4.1/4.2, no other blockers;
                           4.6 DONE; 4.7 (declaration consolidation) STRONGLY PREFERRED
                           before 5.3 — and its 4.7.1 bug-fix step is independent and
                           should land regardless (fixes a live off-by-one)
Phase 5 (package API)   ── independent of Phase 4 EXCEPT 5.3, which should follow 4.7;
                           order: D8 .flopy backfill → 5.1 → 5.2 (all DONE) → [4.7] →
                           5.3A → 5.3B → 5.6 → 5.4 → 5.5 → 5.8 (anytime after 5.3;
                           its viz side is 6.4)
                           → 5.7 (optional/forward-looking; gates NOTHING)
Phase 6 (model types)   ── 6.0/6.1 need 5.3 (buildable GWT models for fixtures);
                           6.2 after 6.1; 6.3 and 6.4 independent of 6.1/6.2;
                           6.1.5 (ConcTargets) composes with 4.1's observations split —
                           land whichever comes second on top of the first
Phase 7 (flopy4/robustness) ── after Phase 2; 7.1 MUST land before Phase 8; 7.1/6.0
                           in either order (see 7.1); 7.4 ships inside Phase 7's PR
Phase 8 (plotting)      ── after 4.1, 5.1–5.6, landed Phase 6 work, and 7.1
                           (5.7 does NOT gate it)
Phase 9 (wrap-up)       ── last
```

**Test-runtime budget (rev. 4):** the canonical fixtures are function-scoped today, so
every slow e2e re-runs MF6 from scratch; Phases 5–6 mandate many new slow e2e tests
plus new GWT/GWE/PRT fixtures, so the ~14-min full suite plausibly doubles or triples
unmanaged. Before 6.1: (a) make the canonical + transport run fixtures session- or
module-scoped (one MF6 run, many assertions); (b) record the full-suite wall time at
each phase end next to the Phase 0 baseline — tripwire at 2× the baseline; (c) the
slow tests run automatically in 1.3's scheduled slow lane.

| Item | Effort | Risk |
|------|--------|------|
| 1.1 figs vendoring | M | Med (figs restructured at `bb1526b` — re-trace at HEAD; CI exercises the snapshot) |
| 1.2 deps / 1.3 CI / 1.4 lint | S / M / M | Low–Med |
| 2.1–2.4 junk/hygiene | S–M | Low |
| 3.1–3.2 deprecation | S | Low |
| 4.1 observations split | M–L | Med |
| 4.2 model_group split | L | Med-high (Phase 6 adds Group* classes here — do 4.2 first) |
| 4.3 surfaces / 4.4 layering | S–M / M | Low |
| 4.5 second-tier splits (triangle, package_plotting) | M | Low–Med |
| 4.7.0 golden snapshots / 4.7.1 live-bug fixes | S / S–M | Low (pure safety net; 4.7.1 fixes wrong numbers) |
| 4.7.2 registry descriptor | M | Low (data only, zero behavior change) |
| 4.7.3 derive the ~20 lists | M–L | Med (one commit each, snapshot-guarded) |
| 4.7.4 collapse spec/resolver/helper triplication | L | Med–High (most public API touched — do last) |
| 4.7.5 `_ArealBuilder` + EVTBuilder — DONE 2026-07-21 | M | Med (refactors shipped RCH behavior; landed with RCH behavior unchanged) |
| D8 .flopy backfill (chd/ghb/drn/wel) | S | Low |
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

S ≈ ≤half day, M ≈ 1–2 days, L ≈ 3–5 days, XL ≈ 1–2 weeks. (Rev. 4: phase headers now
equal the min–max sum of their item rows — rev-3's headers understated Phases 1/4/5/6/7
by up to ~2×. Phase 5's header excludes 5.7, which is optional/forward-looking.)

---

## Appendix A — Verification commands

Linux (the repo's current home — primary):

```bash
py=/home/lukem/python/envs/gw/bin/python   # fully provisioned (Python 3.14)
# myflopy is installed EDITABLE in gw (2026-07-16): imports resolve to the repo
# with or without PYTHONPATH=src (exporting it stays harmless).
# Sanity check: $py -c "import myflopy; print(myflopy.__file__)" → repo path

$py -m pytest -m "not slow" -q            # fast suite (after every change-set)
$py -m pytest -q                          # full suite (end of each phase; ~14 min)

# import surface smoke
$py -c "import myflopy; [getattr(myflopy, n) for n in myflopy.__all__]; print('ok')"

grep -rn "lukem" src/myflopy              # must be empty after Phase 2
grep -rn "X" src tests examples docs --include="*.py"   # orphan check before deleting X

# deferred-import ratchet baseline (198 on 2026-07-14)
grep -rn "^\s\+from myflopy" src/myflopy --include="*.py" | wc -l

# hover/colorscale policy quick checks
$py -m pytest tests/test_hover_spec.py tests/test_colorscale_policy.py -q
git ls-files "*.ipynb" | xargs $py -m nbstripout --dry-run
```

Windows (mf-env box — same commands via PowerShell):

```powershell
$env:PYTHONPATH = "src"
$py = "C:\Users\lukem\Python\mf-env\.venv\Scripts\python.exe"
& $py -m pytest -m "not slow" -q
& $py -m pytest -q
& $py -m nbstripout --dry-run (git ls-files "*.ipynb")
```

## Appendix B — Size snapshot (2026-07-14, for drift detection)

`src/myflopy` 61,682 lines / 134 .py files (rev-3's "~90 files" was wrong even against
its own tree; the +4.2k lines since rev. 3 are dominated by the `962654d` docstring
commit). Largest: `observations.py` 3,014 · `project/model_group.py` 2,761 ·
`package_api.py` 2,500 (growth is docstrings, not logic) · `specs.py` 2,315 ·
`package_plotting.py` 1,914 · `package_surface_water.py` 1,852 · `grid/triangle.py`
1,822 · `interactive_plotting.py` 1,394 · `datatypes/choros.py` 1,322 ·
`simulation/base.py` 1,262 (`pest/ies.py` 1,255 is #11). Tests: 51 files, ~15.9k
lines, ~540 fast-passing (conftest auto-marks ~47 slow: 25 decorators + `_SLOW_TESTS`
+ `canonical_run` users — the raw `def test_` count is 587). Tracked repo 48.6 MB
(notebooks 31.5 MB + mp3du exes 11.2 MB — Phase 2 targets); `.git` history is
~170 MiB (out of scope, D3).

## Appendix C — Master acceptance checklist

- [x] figs committed in its repo (`bb1526b`, 2026-07-07) AND vendored (2026-07-16:
      `_vendor/figs` @ bb1526b via `scripts/sync_vendored_figs.py`; figs-hidden
      subprocess tests prove `import myflopy.viz` + post-script injection work
      without local figs; wheel verified to contain the snapshot) (0, 1.1)
- [x] `v0.1.0` tagged at baseline (2026-07-16); wheel-build CI job written and
      verified locally; scheduled slow lane exists (0, 1.3, D11)
- [x] `matplotlib` + `seaborn` + `openpyxl` declared; dash/osgeo lazified with
      subprocess isolation tests; import-surface smoke green (`pyyaml` lands with
      5.6) (1.2)
- [x] CI green on ubuntu + windows (fast suite + ruff + scoped mypy); vendored-figs
      path exercised — all 8 jobs green 2026-07-17 (`cc62159`) after 7 fix
      iterations; the final Windows blocker was backslash-vs-slash paths in the
      figs AST-guard test itself (1.3, 1.4)
- [x] `grep -rn "lukem" src/myflopy` empty (2026-07-16); junk modules gone incl. the
      `get_iheads`/`gwt` lazy-export + smoke-test removals (D10);
      `test_retired_modules_stay_retired` pins it; live notebooks/examples use
      home-relative paths (2.1, 2.2)
- [x] Notebooks stripped (26, −776k lines) + pre-commit hook + strip/render
      scripts; no `.exe` under `src/` (gitignored forever); `tools/mp3du/`
      resolution chain tested (6 tests) (2.3, 2.4)
- [x] Deprecated aliases warn AND are hidden from autocompletion per D12 (not in
      `__all__`/`dir()`/TYPE_CHECKING/stubs); GHB/DRN collision resolved; policy doc
      exists (3.x) — done 2026-07-16: `_deprecation.py` (module + instance
      `__getattr__` helpers, registry), `__compatibility__` repurposed (old
      second-tier meaning renamed `__engine__`), model_group/mp3du ad-hoc
      mechanisms absorbed, `error:myflopy:DeprecationWarning` pytest filter
- [x] `observations/` + `project/group/` packages; old paths work; no active module
      > ~1,800 lines without a written reason (4.1, 4.2) — done 2026-07-17; the
      four modules still over the line (`specs.py` 2,333, `package_plotting.py`
      1,906, `package_surface_water.py` 1,852, `grid/triangle.py` 1,822) carry
      their written reasons in §4.5 / §5 (cohesive or split planned there)
- [x] Layering test + deferred-import ratchet in CI (4.4) — done 2026-07-17:
      map derived from the real AST graph (`scripts/derive_import_layers.py`;
      143 modules, ALREADY acyclic, depths 0..13); honest ratchet baseline 75
      by AST count (the old 198 grep also matched TYPE_CHECKING imports),
      burned down to 55 (20 hoists)
- [ ] `.flopy` backfilled on chd/ghb/drn/wel with docstrings fixed (D8); `mf.riv`,
      `mf.evt` with `()/.gpkg/.flopy` + registry + hover + earth/RdBu policy (5.1, 5.2)
- [ ] `mf.ic/oc/disv` dispatch on GWT/GWE; `mf.adv/dsp/mst/ssm/cnc/src/ist` +
      `mf.est/cnd/ctp/esl` exist; example 03 uses them (5.3)
- [x] `mf.dis`/`mf.disu` passthroughs; fail-fast GridSpec tests untouched (5.5) — DONE 2026-07-23
- [x] `SimulationSpec.to_yaml/from_yaml` round-trips incl. exchanges/hooks/list BCs (5.6) — DONE 2026-07-23 (YAML; TOML deferred, ledger 54)
- [ ] ONE generic dependent-variable surface; `model.conc` + `model.temp` with full
      grammar (map/xs/plot/mosaic/animate), `conc_hover`/`temp_hover`, earth + RdBu-diff
      colors, GWT/GWE budget terms, `GroupConc`/`GroupTemp` + `diff()` maps,
      `ConcTargets` wired into PEST, GWT + GWE test fixtures (6.0–6.2)
- [x] PRT — DONE: `mf.mip`/`mf.prp`/`mf.ems` factories + prt6 dispatch (6.3A, 2026-07-24);
      `results.travel_time`/`.endpoints`/`.capture` choropleth nouns with hover + policy
      colors (6.3B, 2026-07-25); `results.pathlines` plotly map + `pathline_hover`,
      `Choro` overlays carried through `viz.mosaic`, `PALETTE.categorical` (6.3C, 2026-07-25)
- [x] IES **(§6.4 COMPLETE 2026-07-26)**: `plot_field` hover + per-stat diverging/earth
      colors (`parameter_field_hover`, `_field_map_policy`, `build_choropleth` widened —
      6.4A); the uncertainty map as `plot_field(stat="reduction")`, NOT the originally
      planned `field_uncertainty_map` (half of it already shipped in 6.4A — ledger 75);
      `plot_field_mosaic`, not `field_mosaic`, over a new `viz.mosaic(colorbar=)`;
      `obs_residuals`/`plot_obs_residuals` (6.4B); every IES figure verified on the viz
      front door, the last color literals retired, and the preferred-API docs (6.4C)
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
