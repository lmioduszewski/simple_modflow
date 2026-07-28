# ModelDiff cheat sheet

Compare two (or more) MODFLOW 6 models with one verb: `diff()`. It answers two
questions — **"are these the same model?"** (setup) and **"do they behave the
same?"** (results) — across every package, the config, connection geometry, and
the computed outputs.

- `diff()` compares each model against a **reference** (baseline). N‑way is
  reference‑star: reference vs each other model.
- Everything is **zero‑based** (period 0 = first stress period).
- Setup tiers work **pre‑run**; the results tier needs **completed runs**
  (`.hds` / `.lst` / `.cbc` on disk).

---

## 0. Load and build the diff

```python
import myflopy as mf

ref = mf.load_mf6_run(r"C:\Users\lukem\mf6\Lkpt_F9b")     # baseline / reference
new = mf.load_mf6_run(r"C:\Users\lukem\mf6\F9_maxflow")   # compared

# Two doors — both return the same ModelDiff:
diff = ref.diff(new)                                       # quick; ref is the reference
# diff = ref.diff([new, other])                            # N-way (reference-star)
# diff = ref.diff(r"C:\Users\lukem\mf6\F9_maxflow")        # accepts a path directly
diff = mf.ModelGroup(                                       # named; best for >2 models
    {"Lkpt_F9b":   r"C:\Users\lukem\mf6\Lkpt_F9b",
     "F9_maxflow": r"C:\Users\lukem\mf6\F9_maxflow"},
    reference="Lkpt_F9b",
).diff()
```

## 1. Headline reports

```python
print(diff.report())               # setup faithful-copy report (packages + config + connections)
print(diff.report(results=True))   # ^ plus heads + budget results (reads output files)
diff.summary()                     # package structural/value matrix
diff.model("F9_maxflow").report()  # focus one model
diff.reference, diff.model_names
```

---

## 2. Setup diff — package inputs (structural + value)

Every package is a node with the same shape as single models and groups:
`diff.packages.<pkg>.inputs` (declared data) and `.results` (computed outputs).
Cell-based BC packages: `ghb`, `chd`, `drn`, `wel`, `rch`.

```python
diff.packages.ghb.inputs.cells()      # STRUCTURAL: only_in_reference / only_in_model / shared
diff.packages.ghb.inputs.values()     # VALUE: bhead/cond + *_diff on shared cells
diff.packages.ghb.inputs.summary()    # per-package counts
diff.packages.ghb.inputs.map("F9b", value_column="cond")  # diverging Δ choropleth
# a choropleth is per-model: pass the model name (positional or model_name=) when the
# group has >1 non-reference model; it's inferred when there's exactly one.
```

## 3. Setup diff — configuration (TDIS / IMS / OC / package options)

```python
diff.config.settings()                 # every differing setting
diff.config.settings(section="tdis")   # timing / ATS / per-period nstp, tsmult
diff.config.settings(section="ims")    # solver block
diff.config.summary()                  # per-model count of differing settings

# inspect ONE model's config (no diff):
ref.config.settings()
ref.config.ims            # solver block as a dict
ref.config.tdis           # timing, incl. per-period rows
ref.config.section("lak") # any package's OPTIONS
ref.config.sections       # which sections exist
```

## 4. Setup diff — connection geometry (LAK / SFR)

For advanced packages the *inputs* tier is the connection/reach geometry.

```python
diff.packages.lak.inputs.connections()                    # belev/telev/connlen/connwidth diff
diff.packages.lak.inputs.connections().query("lake == 0") # confirm Lake 0 is identical
diff.packages.sfr.inputs.reaches()                        # reach geometry (reach/layer/cell/rlen)
diff.packages.lak.inputs.summary()                        # per-model connection counts
diff.packages.sfr.inputs.summary()
```

---

## 5. Results diff — behavior (needs completed runs)

Every `.summary()` carries a `within_tolerance` flag (`|Δ| ≤ atol + rtol·|ref|`,
defaults `1e-3`, overridable) and `argmax_*` columns for **where / when** the
biggest difference is.

```python
# heads -- a surface-root leaf, exactly like model.hds / group.hds
# (diff.conc / diff.temp are the transport twins, same verbs)
diff.hds.summary()                # max/mean|Δ|, RMSE, argmax cell/layer/kstpkper, within_tolerance
diff.hds.get(per=8, layer=1)      # aligned elev / reference_elev / diff per cell
diff.hds.map("F9b", per=8, layer=1)   # Δhead choropleth (model - reference)
diff.hds.plot(layer=0)            # mean Δhead by stress period, line per model
diff.hds.xs(line=line)            # reference vs model head profiles along a section line

# overall (volumetric) budget — term by term
diff.bud.summary()                # per term: totals, % change, within_tolerance
diff.bud.get()                    # per-timestep per-term diffs

# per-package cell budget (leakage) + UZF
diff.packages.ghb.results.summary()        # where did GHB leakage change (field q)
diff.packages.uzf.results.gwrch.summary()  # UZF recharge; also .sat

# lakes & streams — q + stage  ← the facility-stage question
diff.packages.lak.results.stage.summary()  # Δstage: worst lake + period, within_tolerance
diff.packages.lak.results.stage.timeseries()   # full per / lake / stage_diff table
diff.packages.lak.results.q.summary()      # Δ lake-groundwater exchange
diff.packages.sfr.results.stage.summary()  # Δreach stage
diff.packages.sfr.results.q.summary()      # Δ stream-groundwater exchange

# mover
diff.packages.mvr.results.summary()        # FROM-MVR/TO-MVR per package (empty if no mover)

# spatial Δ maps (per-cell flux difference choropleth, per model):
# all diff maps use the diverging red/white/blue scale (negative = red,
# positive = blue) and a purple-accented sectioned hover showing the model
# value, the reference value, and Δ (hover_* sugar works here too).
diff.packages.sfr.results.q.map("F9b")     # Δ stream leakage (interactive Choro)
diff.packages.lak.results.q.map("F9b")     # Δ lake exchange
diff.packages.ghb.results.map("F9b")       # Δ GHB leakage (also drn/chd/wel/rch)
diff.packages.uzf.results.gwrch.map("F9b") # Δ UZF recharge

# UNIFIED GRAMMAR -- every leaf (single-model / group / diff, inputs & results)
# has the same verbs: get/summary + map/plot/xs + mosaic/animate, all with backend=:
diff.hds.mosaic(by="model")                 # synced Δ mosaic, one panel per model
diff.hds.mosaic(by="model", backend="mpl")  # static (matplotlib) instead of Plotly
diff.hds.animate(over="period")             # Δ animated across stress periods (play/slider)
diff.hds.animate(kind="xs", line=line)      # animated cross section (ref vs models)
diff.packages.ghb.inputs.mosaic(by="model") # setup-value Δ mosaic
diff.packages.ghb.results.map(field="q")    # field= sugar == .results.q.map()

# NOTE: diff.packages.lak/.sfr INPUTS (connection geometry) have NO .map() --
# that tier is a set-diff (added/removed/changed), not a scalar field.
# The per-cell Δ map lives under .results.

# tighten / loosen the pass/fail on any results accessor:
diff.hds.summary(atol=0.01, rtol=0.0)
```

---

## 6. Reading the output

| Situation | Meaning |
|---|---|
| `within_tolerance = True` everywhere | identical within tolerance (regression pass) |
| `within_tolerance = False` | a real difference — check `argmax_*` for where/when |
| setup identical **+** results differ | same model, different answer (solver / numeric drift) |
| setup differs **+** results differ | your change moved the output (scenario) |

For the facility question: `diff.packages.lak.results.stage.summary()` → the
lake and period where stage rises most vs baseline.

## 7. Caveats

- **Results need runs.** Loaded completed runs are fine; unrun models raise.
- **Cell-level results** (heads, cell budget) assume the **same grid / cell
  numbering**. If two runs use different Voronoi meshes, lean on the
  grid-independent views: `lak/sfr .results.stage`, overall `bud`, `config`.
- One tree everywhere: `diff.packages.<pkg>.inputs` (`ghb/chd/drn/wel/rch`
  cells+values; `lak/sfr` connections) and `.results` (cell budgets, `uzf`
  gwrch/sat, `lak/sfr` q+stage, `mvr`); heads live at `diff.hds`, budget at
  `diff.bud` — mirroring `model.hds` / `group.hds`.
- `group.rch` (flat shortcut) still works but warns — use
  `group.packages.rch.inputs`, which mirrors `model.packages.rch.inputs`.
- Reference = baseline; everything zero-based.

---

## Quick reference

| Want | Call |
|---|---|
| Full setup report | `diff.report()` |
| Report incl. results | `diff.report(results=True)` |
| Which cells added/removed | `diff.packages.<pkg>.inputs.cells()` |
| Value differences per cell | `diff.packages.<pkg>.inputs.values()` |
| Config differences | `diff.config.settings()` |
| Lake/stream connection geometry | `diff.packages.lak.inputs.connections()` / `diff.packages.sfr.inputs.reaches()` |
| Head differences | `diff.hds.summary()` |
| Water-balance differences | `diff.bud.summary()` |
| Package leakage differences | `diff.packages.<pkg>.results.summary()` |
| Lake / reach stage differences | `diff.packages.lak.results.stage.summary()` / `diff.packages.sfr.results.stage.summary()` |
| Mover differences | `diff.packages.mvr.results.summary()` |
| One model's config | `model.config.settings()` |
