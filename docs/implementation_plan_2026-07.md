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
| D4 | Structured-grid stance | **Add thin `mf.dis` / `mf.disu` passthroughs** (model-type-aware). `GridSpec.structured` stays fail-fast. Details in Phase 5.5. |
| D5 | Plotting consolidation | **Do it**, as its own late phase (Phase 8), after the god-module splits (Phase 4), API completion (Phases 5–6), and the `_flopy_compat` boundary (Phase 7.1). Details in Phase 8. |
| D6 | Hover defaults | DONE (see "Completed"). Heads default `layers="active+strip"`; sectioned/styled hover is the default everywhere; `custom_hover` stays as the raw escape hatch. |
| D7 | Colorscales | DONE (see "Completed"). Diverging only for signed q-like + diff maps; `'earth'` for everything else. New surfaces added by this plan MUST follow this policy. |
| D8 | `.flopy` on simple list BCs | **Add a real `.flopy(...)` escape hatch** to `mf.chd/ghb/drn/wel` (they never had one — only rch/uzf/sfr/lak/mvr do, though the chd/drn/wel docstrings advertise it) and ship `mf.riv`/`mf.evt` with all three entry points. Details in Phase 5.0. |
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

---

## Phase 5 — Package API completion (~8–15 days excluding 5.7, independent of Phase 4)

### 5.0 The pattern to replicate (read first)

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

### 6.0 Design keystone: ONE generic dependent-variable surface

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
2. **Pathline map hover:** CORRECTION (rev. 4): the existing pathline map
   (`plot_particle_pathlines`, interactive_plotting.py) is pure **matplotlib** (flopy
   PlotMapView), not plotly scatter — so this step means building a NEW plotly
   pathline map (scatter traces over the grid choropleth), not styling existing
   traces. The hover engine is trace-agnostic (customdata + template), so add a
   `pathline_hover()` spec rendering particle id, release point/group, current time,
   layer; apply `HoverStyle.to_hoverlabel()` to the new scatter traces. Read
   `prt.py:PRTRunResults.plot_map` + `interactive_plotting.py` scene builders first.
3. **Build-side parity:** add `mf.mip(porosity=)` and `mf.prp(...)` thin factories
   (5.3A dispatch pattern) so a PRT model is fully declarable in a `SimulationSpec`
   without raw PackageSpecs; keep `model.particle_tracking.prt(...)`/`PRTProject` as the
   sanctioned high-level runtime path (it manages FMI/grid copying correctly — do not
   duplicate that logic in specs).
4. **3-D scene** stays as-is — CORRECTION (rev. 4): it is **PyVista/trame** (not raw
   plotly); either way it remains a documented viz.py exception. MP3DU untouched
   except Phase 2.4.
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
                           the class inventory); 4.5 after 4.1/4.2, no other blockers
Phase 5 (package API)   ── independent of Phase 4;
                           order: D8 .flopy backfill → 5.3A → 5.1 → 5.2 → 5.3B → 5.6
                           → 5.4 → 5.5 → 5.8 (anytime after 5.3; its viz side is 6.4)
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
