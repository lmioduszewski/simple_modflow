# Phase baselines — test-suite results per implementation-plan phase

Companion to `docs/implementation_plan_2026-07.md` (rev. 4). One section per phase
close-out; the Phase 0 entry is the reference every later phase compares against
(test-runtime tripwire = 2× the Phase 0 full-suite wall time).

## Phase 0 — Baseline (2026-07-16)

- **Commit:** the commit tagged `v0.1.0` (D11) on branch `myflopy` — plan rev-4
  (`d7a39e1`) + this baseline record + the PEST++ POSIX command fix. figs repo at
  `bb1526b`, clean.
- **Machine/env:** Linux, `/home/lukem/python/envs/gw` (Python 3.14.4, flopy 3.10.0),
  `PYTHONPATH=src`. MODFLOW executables installed to `~/.local/bin` via
  `python -m flopy.utils.get_modflow` (mf6 6.7.0, triangle 1.6, + full bundle).
- **Fast suite** (`pytest -m "not slow" -q`): **532 passed, 1 skipped, 1 failed,
  ~10 s wall** (the Windows box's ~60 s figure does not apply here).
- **Full suite** (`pytest -q`): **580 passed, 1 skipped, 2 failed, 9m51s wall**
  (canonical run, fully provisioned, incl. the PEST++ command fix). Test-runtime
  tripwire for later phases: **20 min** (2× baseline).

### Known failures (environment/repo defects, NOT regressions)

1. `test_mf6_refactor_smoke.py::test_examples_layout_smoke` — **repo defect**: asserts
   the gitignored `examples/mf6/artifacts/` dir and the **untracked**
   `examples/mf6/sample_model_output/top_raster.tif` exist. Fails on any fresh clone
   and will fail in Phase 1.3 CI. Fix (with 1.3 or Phase 2): track or generate the
   sample data, mkdir-or-skip the artifacts assert.
2. `test_parallel_model.py::test_canonical_eight_part_split_runs_with_mpi` — the
   standard MODFLOW executables bundle ships a **serial** mf6 ("Can not run parallel
   mode with this executable: no MPI"). Needs an MPI-enabled mf6 build to pass here.
   Bonus finding: `parallel.environment.parallel_ready` is a **false positive** — it
   detects `mpiexec` but not whether mf6 is an MPI build, so the test's skip guard
   doesn't fire; the readiness check should probe the binary (e.g. parse
   `mf6 -v`/trial run) and the test would then skip cleanly.

### Environment provisioning performed (reproduce on a fresh Linux box)

```bash
gw=/home/lukem/python/envs/gw/bin/python
$gw -m pip install pytest nbstripout ruff mypy build   # plan ground rule 1 / 1.4
$gw -m pip install "dash>=2.15" "dash-bootstrap-components>=1.5"   # viz extra
$gw -m pip install pyvista trame trame-vtk trame-vuetify           # viz3d extra
$gw -m pip install pyemu                                            # pest extra
$gw -m pip install scikit-learn pymetis xugrid          # parallel + xugrid extras
# (h5py was already present)
# GDAL: no pip wheel builds here; system GDAL 3.12.2 symlinked into the env:
#   ln -s /usr/lib/python3/dist-packages/osgeo <gw site-packages>/osgeo
$gw -m pip install -e .    # editable — replaces the stale site-packages copy
                           # that used to shadow the repo without PYTHONPATH
$gw -m flopy.utils.get_modflow ~/.local/bin            # mf6, triangle, ...
# PEST++ 5.2.16 (pestpp-ies etc.) from github.com/usgs/pestpp releases
# (linux tar.gz → bin/pestpp-* → ~/.local/bin)
mkdir -p examples/mf6/artifacts                        # gitignored, asserted by a test
```

### Findings feeding Phase 1.2 (dependency audit) and test hygiene

- `dash` is a **hard top-level import of core code** (`utils/datatypes/choros.py`,
  reached via `grid/__init__.py` → `grid/plotting.py`) but is declared only in the
  `viz` extra — promote it, lazify the import, or guard it (1.2).
- `osgeo`/GDAL is a hard top-level import in `utils/raster.py` with no pip story —
  document as a system dependency (1.2).
- Three `test_mf6_pest.py` tests (`...builds_native_pst_with_multilayer_k`,
  `...grid_k_parameterization_on_voronoi_with_capture`,
  `...pilot_point_k_parameterization_on_voronoi`) **raise** on missing `pyemu`
  instead of `pytest.importorskip("pyemu")` like their siblings — add the guard.
- PEST++ 5.2.16 installed. **Linux bug found and fixed during baselining**:
  `pest/project.py` wrote the PST model command as `"{sys.executable}"
  forward_run.py` — POSIX PEST++ neither shell-parses quotes nor accepts absolute
  command paths (its run manager mangles the leading `/` and `execv` fails), so ALL
  IES realizations failed on Linux. Fixed via `build_forward_run_command()`: POSIX
  gets a `./run_forward.sh` wrapper written into the template workspace (travels
  with worker-dir copies); Windows keeps the quoted form. Pinned by
  `test_build_forward_run_command_is_pestpp_compatible` (fast).
- `tests/test_mf6_pest.py` (~line 22) hardcodes a `C:\Users\lukem\...` pyemu
  fallback path — in `tests/`, outside Phase 2.2's `src/myflopy`-scoped acceptance
  grep; remove it with Phase 2.2.
- Deferred-import ratchet baseline: **198** (`grep -rn "^\s\+from myflopy"
  src/myflopy --include="*.py" | wc -l`).

## Phase 1 — Distribution blockers (2026-07-16, branch `phase-1-distribution`)

- **Delivered:** figs vendored (`_vendor/figs` @ figs `bb1526b`, trimmed closure,
  sync script, external-first import in viz.py, single-importer AST rule);
  matplotlib/seaborn/openpyxl declared, dash/osgeo lazified (+ subprocess
  isolation tests); GitHub Actions (fast matrix ubuntu+windows ×
  py3.10/3.12/3.14, ruff+mypy lint job, wheel job, weekly slow lane);
  `ruff check src tests` green (577 safe autofixes + triage; documented
  burn-down ignore list); scoped mypy green on the declarative core (26 errors
  fixed, incl. one latent NameError bug in `budget.py` stream-flow plotting).
- **Fast suite:** 540 passed / 2 skipped / 0 failed, ~18 s.
- **Repo defects fixed:** layout smoke test no longer requires untracked local
  data (mkdir + skip); pyvista tests importorskip; pest tests' pyemu guard
  gap remains for 3 tests (they raise a clear ModuleNotFoundError — acceptable).
- **Still open for Phase 1 acceptance:** first push → remote CI green on both
  OSes (workflows are verified locally: wheel contains `_vendor`, py310
  grammar parse clean, figs-blocked smokes pass).

## Phase 2 — Junk removal & repo hygiene (2026-07-16, branch `phase-2-hygiene`)

- **Delivered:** 11 modules retired (8 to `attic/` with a README ledger, 3
  deleted incl. the D10 public-export removals for `get_iheads`/`modflow.gwt`);
  `test_retired_modules_stay_retired` guard; all `C:\Users\lukem` paths gone
  from `src/` and live examples (home-relative replacements resolve identically
  on the author's Windows boxes); all 26 notebooks output-stripped (−776k
  lines, every notebook now <100 KB) with a pre-commit nbstripout hook +
  `scripts/strip_notebooks.py` / `render_notebooks.py`; unreferenced
  `starting_heads.csv` deleted; MP3DU executables moved out of `src/` to the
  untracked `tools/mp3du/` with a tested 4-step resolution chain
  (`MYFLOPY_MP3DU_DIR` supported); `src/**/*.exe` gitignored.
- **Fast suite:** 555 passed / 1 skipped, ~20 s. Ruff + scoped mypy green.
- **Note:** the tracked checkout shed ~42 MB (31.5 MB notebook outputs +
  11.2 MB executables); `.git` history keeps the old blobs (D3 — no rewrite).

### Full-suite runs (append per phase)

| Date | Commit | Result | Wall time |
|------|--------|--------|-----------|
| 2026-07-16 | `d7a39e1` | 570 passed / 7 skipped / 5 failed (pre-pyemu, pre-PEST++) | 7m59s |
| 2026-07-16 | `d7a39e1` | 576 passed / 1 skipped / 5 failed (pyemu installed; 3 IES fails = PEST++ POSIX command bug) | 8m10s |
| 2026-07-16 | `d7a39e1`+fix | **BASELINE: 580 passed / 1 skipped / 2 failed (both known: layout smoke, MPI)** | **9m51s** |
| 2026-07-16 | Phase 1 end (`phase-1-distribution`) | 588 passed / 2 skipped / 1 failed (only the MPI environment limitation; layout smoke now skips; +12 new Phase-1 tests; parallel/xugrid extras now exercised) | 11m23s (within the 20 min tripwire) |
