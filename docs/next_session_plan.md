# Next session: close ledger 103, 105, 107

Scoped 2026-07-28. Working tree clean at `7b80775`. Nothing started — this is the
plan only. Three independent items; A and B can be done in either order, C is the
one with the widest blast radius.

All three forks below were **decided by the user**; they are not open questions.

---

## A · Concentration observations on IES residual maps (closes ledger 105)

Smallest item. Genuinely wiring — no new machinery.

**Key finding:** `_observation_locations` (`pest/ies.py:1588-1657`) is the ONLY
kind-specific step in the whole residual pipeline. Everything downstream
(`obs_residuals` `:1702`, `plot_obs_residuals` `:1781`) splits on **column
nullity** — rows with non-null `x` become markers, rows with non-null `cells`
become choropleth patches — and never branches on `kind` again.

**The crux:** `prepare_conc_observations` (`pest/observations.py:299`) calls
`match_to_model`, which already returns `x`/`y`, then slices only
`["name","layer","cell"]` and throws the coordinates away. We are discarding
exactly what the map needs.

Change surface:

| file:line | change |
|---|---|
| `observations.py:299` | keep `x`/`y` from `match_to_model` |
| `observations.py` (conc) | write `<prefix>_target_locations.gpkg` (reuse `_write_locations_snapshot`, as head does at `:203`) |
| `observations.py:322-335` | metadata gains `locations_file`, `geometry: "points"`, and the mapping filename |
| `ies.py:1604` | replace the `kind in ("head_targets","drn_flow")` tuple |
| `ies.py:1623` | replace `if kind == "head_targets"` with `if geometry == "points"` |
| `ies.py:1629` | un-hardcode the `_head_target_map.csv` suffix — read it from metadata |
| `ies.py:1820-1826` | error text claims only heads/DRN have geometry |

**Dispatch design:** add a `geometry` key to the metadata (`"points"` / `"zones"`),
with a **fallback keyed on `kind`** for run directories already written to disk
(`head_targets`→points, `drn_flow`→zones). Old runs must keep reviewing.

Tests: `_stub_ies_results` (`tests/test_mf6_pest.py:1665`) needs a conc branch —
today every entry it builds is head or drn. `tests/test_mf6_pest.py:1891` pins the
current refusal message and will need updating. Add a real end-to-end assertion on
the transport demo.

Docs: `PEST_calibration_guide.md:260-274`, `package_api_reference.md:284`,
`myflopy_context.md:148`, ledger 105 → closed.

**Constraint:** `prepare_conc_observations` must keep its deferred `load_mf6_run`
import (ledger 104) or `test_runtime_imports_never_point_upward` fails.

---

## B · Transport `parameterize` targets (closes ledger 103)

**Decided: `mst.porosity` only.** NOT `dsp.alh` — `ath1` is derived from `alh` at
build time only, so once externalized a multiplier on `alh` alone silently breaks
the 10:1 transverse ratio the model was built with.

**The easy half.** pyEMU's `apply_list_and_array_pars` is entirely filename-driven
(`pyemu/utils/helpers.py:1940`, `_process_array_file` `:2013`) — it never parses
model names or package types, so a GWT array works in the forward run with zero
changes. The model name enters at exactly one line:

- `native_parameters.py:206` — `stub = recipe.file_stub.format(model=model_name)`
- fed from `project.model.name` (the GWF model) at `native_parameters.py:239` and
  `pilot_points.py:96`

So: `_Recipe` gains a `model: "flow"|"transport"` field (default `"flow"`), and
`_resolve_files` resolves the transport name through the **existing**
`resolve_transport_model_name` (`pest/observations.py:240`) — the precedent the
observation side already set.

Verified file layout (flat coupled sim, `set_all_data_external()`):
```
viz_prt_master.npf_k_layer1..4.txt      <- flow, PER-LAYER
viz_prt_master_t.mst_porosity.txt       <- transport, ONE FILE (no LAYERED keyword)
```
`_resolve_files`' exact-name branch (`:210-213`) hits first, so single-file
porosity needs no change there.

**The hard half — three guards, all for silent-wrongness:**

1. **`style='pilotpoints'` on a transport target must be refused.**
   `pilot_points.py:109` hardcodes `project.model.gwf.npf.k` as the base array
   regardless of target, and calls it `base_k` through the forward run
   (`forward_run.py:40`). Porosity pilot points would interpolate against K.
2. **`layers=` on a single-file target must be refused.** `native_parameters.py:245`
   filters on `_layer(\d+)\.txt$`; with one file `selected` is empty and `:248`
   (`if selected:`) **silently keeps all files** — so `layers=[0]` appears to work
   and does nothing.
3. **`capture=True` on a non-layered file.** `_add_capture_field_observations`
   (`project.py:763-790`) leaves `layer_prefixes` empty, so the captured field has
   no layer identity for `plot_field`. Handle it or refuse it.

`parameterize` is NOT in `tests/api_snapshot.json`, so adding a target will not
move the snapshot.

Docs: `project.py:255-262` target list, `CLAUDE.md:161-165`, ledger 103 → closed.

### B2 · The demo must use heads AND concentration

**Decided.** Measured on the canonical transport model (testing profile):

| perturbation | head response | conc response |
|---|---|---|
| K x3 | **1.098 ft** | max dC 0.157 |
| porosity /3 | **exactly 0.000000 ft** | max dC 0.154 |

Cosine similarity of the two concentration responses at the monitoring wells:
**0.98** — near-collinear. Physically expected, since transport velocity is
`v = Ki/n`, so tripling K and dividing porosity by three do nearly the same thing.

So estimating K and porosity **jointly from concentration alone is ill-posed**.
Porosity is absent from the flow equation, so heads pin K and concentration then
pins porosity. `build_canonical_transport_calibration_demo` needs a `HeadTargets`
family added, and `canonical_07` a section explaining why two data types are
required. This is the scientifically correct story, not just a workaround.

---

## C · Move the PEST template out of the copied tree (closes ledger 107)

**Decided: `<model workspace>.pest/<name>`** — a sibling of the model directory
(model at `.../viz_prt_master`, template at `.../viz_prt_master.pest/<name>`).
Discovery searches the new root AND the legacy `<ws>/pest` so runs already on disk
stay visible.

**Why not exclude the subtree instead: it is impossible via supported API.**
pyEMU 1.4.0 copies with `shutil.copytree(o_d, n_d, symlinks=True)` in the private
`_try_copy_dir` (`pyemu/utils/os_utils.py:256-268`), called from
`PstFrom.__init__` itself (`pst_from.py:285` -> `_setup_dirs` `:956-980`). There is
no `ignore=`, no kwarg on `PstFrom.__init__` (`:213-227`), no hook. Only a
monkeypatch of a private function would work — rejected: a silent break on any
pyEMU upgrade brings the 17 GB recursion back.

**Moving is cheaper than ledger 107 assumed.** Three things are NOT coupled:

- `find_pest_runs` (`runs.py:107-175`) is **location-agnostic** — it `rglob`s
  whatever root it is handed. Only the two CALL SITES hardcode `<ws>/pest`:
  `workspace.py:420` and `simulation/base.py:284`.
- The forward run has **zero** dependence on template location — every path is a
  bare basename (`build_forward_run_command` `project.py:84-104`; `model_command`
  is `./run_forward.sh`, relative to the run dir which PEST++ sets to the template).
- `original_workspace` in the metadata (`project.py:622`) is **write-only** —
  grep confirms no reader.

**The one hard constraint:** masters must remain SIBLINGS of the template.
`runs.py:160` globs `template_dir.parent / f"{name}_*master*"`, pairing with
`project.py:1182` (`self.template_workspace.parent / f"{self.name}_{master_suffix}"`).
Break that and discovery silently reports a run as "built (not run)", and
`review()` falls back to the never-run template (`runs.py:79-80`).

After the move, `_clear_stale_template`'s refuse-branch becomes unreachable for the
default (the `inside` check at `project.py:485` early-returns) and remains only for
a user-supplied inside-the-model workspace. Multiple named runs work again.

Change surface: 2 discovery call sites; test assertions at
`tests/test_mf6_pest.py:932`, `:937`, `:1320` (+ discovery at `:1390`); prose in
`CLAUDE.md:172`, `PEST_calibration_guide.md:49,340`, `myflopy_context.md:139`;
and the `rmtree` workaround comments in 5 canonical notebooks — including
`canonical_07:262-270`, which exists only because of this bug and should go.

**Note:** there is currently NO test for `_clear_stale_template` at all
(grep for `clear_stale` / "Cannot build this calibration" in `tests/` finds
nothing). Add one while here.

---

## Verification checklist (applies to all three)

- `pytest -n0 -q` serially before each commit — xdist can mask failures (ledger 48)
- `ruff check src/ tests/`
- regenerate `scripts/derive_api_snapshot.py` and `scripts/derive_import_layers.py`,
  **review both diffs**
- ledger entry in the same pass for every judgment call
- mutation-test each fix: revert it, confirm a test fails
