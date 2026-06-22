# A practical guide to PEST++ calibration & uncertainty with myflopy

A plain-language summary of the methods taught in the GWM² 2026 PEST workshop
(Jeremy White, Mike Fienen, Katie Markovich), mapped to the myflopy calibration
API. The goal is to give you the *why* behind each step and the one or two lines
of myflopy that do it.

The whole thing is one idea from Bayes' rule:

> **posterior  ∝  prior  ×  likelihood**
>
> What we believe *after* looking at the data = what we believed *before* (the
> **prior**), reshaped by how well the model can match the data (the
> **likelihood**). "Learning" is just reducing uncertainty about predictions we
> care about.

Everything below is a piece of that sentence.

---

## 0. Purpose first (and keep it simple)

Before any calibration: know your purpose. Are you *explaining* a system or
*predicting* an unseen outcome? Applied work is almost always prediction, which
means the deliverable is a **forecast with honest uncertainty**, not a single
"calibrated" model. Build the *minimum feasible complexity* — a model fast
enough that you can run it hundreds of times, because you will.

A forecast you can defend > a pretty model you ran once.

---

## 1. The prior — what you're uncertain about (parameterization)

The prior is your expert knowledge expressed as adjustable parameters and their
plausible ranges. Use **lots** of parameters (heterogeneity is real), but
express their spatial/temporal correlation so they stay physically sensible.

Key parameters: hydraulic conductivity (K), recharge, storage, riverbed/GHB/drain
conductances, pumping. Two shapes:
- **array-like** (K, recharge, storage) — one value per cell, usually controlled
  by pilot points or multipliers at several scales;
- **list-like** (pumping, boundary heads, conductances) — one value per feature.

In myflopy you *declare* what to adjust; it compiles to native pyEMU:

```python
cal = PestProject(model, "calib", workspace="calib_template", start_datetime="2020-01-01")

# bounds = how far the multiplier can move; physical = hard limits on the final value
cal.parameterize("k",        style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
cal.parameterize("recharge", style="grid",     bounds=(0.5, 1.5), physical=(0, 1e-2))
cal.parameterize("ghb.cond", bounds=(0.1, 10))
```

`physical=` is the safety rail: no matter what calibration does, the *final* K
stays inside those bounds. Always set it for multipliers.

> **Tip from the workshop:** build and inspect the prior *before* history
> matching ("prior Monte Carlo, early and often") — see §5.

---

## 2. The likelihood — how you measure fit (observations & weights)

The **objective function** (phi, Φ) is the weighted sum of squared residuals —
the model's misfit to the data. Minimizing it is "calibration." Each
observation's **weight** scales its influence.

The deceptively hard part is the weights:
- A good starting rule is **weight = 1 / (measurement noise std)** — better data
  counts more.
- But raw weights let high-magnitude or densely-sampled observations dominate.
  The fix is **visibility weighting**: group observations by type and set weights
  so each *group* contributes roughly equally (or by your chosen split). This
  focuses the fit on what matters for your purpose.

Observations are first-class myflopy objects; you hand them to the project:

```python
cal.observe(head_targets)          # history-matching targets (measured values + weights)
```

> Weighting matters *a lot* and there is no one-size-fits-all. Pre-process data
> for information content, use more groups rather than fewer, and aim weights at
> the model's purpose.

### Noise and structural error

Two error sources live in the residuals: (1) **measurement noise** (instruments,
human error) and (2) **structural error** — the model is a simplification, so
even "perfect" parameters can't match reality exactly. Pretending structural
error is zero causes **overfitting**: the model contorts parameters to chase
noise, which *looks* great in the calibration period and predicts badly.

PESTPP-IES carries noise into the answer by fitting each realization to a
slightly different, noise-perturbed copy of the data — so the uncertainty in the
data flows through to the uncertainty in the forecast.

---

## 3. Forecasts — the things you actually care about

A forecast is a model output you want to *predict* but have no data for (a future
flux, an unseen head, a travel time). Declare them so the uncertainty tools track
them. Two workshop habits:
- **Always add a forecast period** to the end of the history model and monitor
  "future" outputs.
- **Frame the question to play to your strengths**: "*Will* the spring flow drop
  below X?" (a probability) is far more answerable than "*What* will the flow be?"

```python
cal.forecast(prediction_targets)   # same target objects as observe(), but zero-weight
```

myflopy registers these as zero-weight PEST++ forecasts, so every uncertainty
method reports their prior→posterior distribution for free.

---

## 4. Build and sanity-check

```python
print(cal.settings())              # review the resolved config — never a black box
pst = cal.build("calib.pst", noptmax=0)   # noptmax=0 = "run once, compute residuals"
```

`noptmax=0` is the universal first check: run the model once through PEST and make
sure the plumbing works before launching anything long.

---

## 5. History matching with ensembles (PESTPP-IES)

Older calibration (PESTPP-GLM) finds *one* best parameter set using gradients
("one model run per parameter"). **Ensemble methods** (PESTPP-IES) instead evolve
a *collection* of randomly-sampled parameter sets (an **ensemble**), approximating
the gradients from the ensemble itself — so the cost is independent of the number
of parameters, and you get uncertainty for free.

### Prior Monte Carlo — "early and often"

Before history matching, run the *prior* ensemble once and look at it: does the
spread of simulated values **bracket** the observed data? If the observations
fall outside the prior, you have **prior-data conflict** — the model as built
can't reproduce the data with plausible parameters. That's a modeling problem to
fix (structure, bounds, weights), not something to calibrate away.

### Run it

```python
ies = cal.run_ies(reals=50, iterations=3, workers=10)   # one line; sensible defaults
print(ies.settings)
```

Defaults make `cal.run_ies()` a valid first call. Production knobs (localization,
covariance re-inflation, bad-phi filtering) pass straight through as keyword
arguments when you need them.

### The DBTL loop

Real calibration is iterative: **Design → Build → Test → Learn**, repeated. Each
cycle you revise parameters, prior, weights, or the model itself and re-run. The
recommended cadence:

| Stage | Settings |
|---|---|
| Sanity | `noptmax=0`, plot & check |
| First real run | 30–50 reals, `iterations=3`, `bad_phi_sigma=1.5` |
| Struggling to fit | 100–200 reals, more iterations, add multimodal/re-inflation knobs |

```python
ies = cal.run_ies(reals=50, iterations=3, bad_phi_sigma=1.5, workers=10)
```

---

## 6. Assess a run — the per-cycle checklist

Each cycle, look at these. myflopy gives one method per item:

```python
ies.plot_phi()                  # did misfit drop? quickly early but not TOO quickly (overfit)?
ies.plot_phi_distribution()     # spread of phi: prior vs posterior — collapse to ~0 = overfitting
ies.plot_vs_obs()               # prior(grey)/posterior(blue) vs measured(red) — does blue bracket red?
ies.plot_phi_contributions()    # which observation groups make up the misfit (visibility weighting)
ies.plot_parameters_at_bounds() # % of parameters pinned at bounds — prior too tight / compensating?
ies.forecasts()                 # prior→posterior uncertainty table for every forecast
ies.forecast("spring").plot()   # the payoff: posterior forecast distribution
ies.plot_field("k", stat="mean")     # property patterns — plausible or laughable?
ies.plot_field("k", stat="std")      # where is K still uncertain?
ies.plot_field("k", stat="change")   # where did calibration move K? (posterior/prior)
ies.best()                      # the single realization to carry forward (the "base" / min-error-variance one)
ies.report("review.html")       # all of the above bundled into one HTML
```

Every plot takes `backend="matplotlib"` for static matplotlib/seaborn output instead of interactive Plotly.

What you're checking for:
- **Phi**: should drop fast at first. If it collapses to near-zero, you're
  overfitting (see §7).
- **Obs vs sim**: the posterior spread should *cover* the measured values without
  being implausibly narrow.
- **Forecasts**: value the *spread*, not a tiny phi. Use `capture=True` on a
  parameter (e.g. `cal.parameterize("k", ..., capture=True)`) to enable the field
  maps.
- **Property maps**: if the calibrated K pattern looks physically absurd, adjust
  the geostatistics/prior — don't accept it just because phi dropped.
- **% of parameters at their bounds**: many parameters pinned at bounds means the
  prior is too tight or compensating for structural error.

---

## 7. The most important lesson: a good fit ≠ a good forecast

The workshop's "pepsi challenge": two runs, one with an *okay* fit and one with a
*crazy-good* fit. The crazy-good fit has a tiny, confident posterior — that
*misses the truth* for some forecasts. The okay fit keeps a wider, more honest
posterior that brackets the truth.

With an imperfect model (every model is imperfect), chasing the data too hard
injects **bias** into predictions that depend on parameter combinations the data
can't constrain. The cure is to **deliberately under-fit**: heavier noise/weights,
fewer iterations, filter signal you can't (or shouldn't) match. If prior
uncertainty is already good enough for the decision, you may not need to history
match at all.

Bias–variance in one line: fitting harder lowers variance but can raise bias, and
in groundwater the bias is usually *unknowable* — so err toward honest uncertainty.

---

## 8. Beyond history matching (where this goes next)

- **Linear uncertainty / data worth (FOSM)** — cheap prior→posterior uncertainty
  and "which observations are worth collecting" without an ensemble run.
- **Management optimization** — PESTPP-OPT (linear, chance-constrained) and
  PESTPP-MOU (multi-objective, Pareto trade-offs under uncertainty): "how much can
  we pump while keeping the spring above X at 95% confidence?"
- **Emulation** — fast surrogate models (DSI, GPR) that stand in for a slow model
  during uncertainty/optimization.

In myflopy today, the run-time knobs for the advanced IES features are available
through `run_ies(**pestpp_options)`; FOSM, optimization, and emulation facades are
on the roadmap.

---

## The whole workflow, start to finish

```python
from myflopy.modflow.mf6.pest import PestProject

cal = PestProject(model, "calib", workspace="calib_template", start_datetime="2020-01-01")

# 1-2-3: prior, likelihood, forecasts
cal.parameterize("k",        style="constant", bounds=(0.2, 5), physical=(1e-3, 100), capture=True)
cal.parameterize("recharge", style="grid",     bounds=(0.5, 1.5), physical=(0, 1e-2))
cal.observe(head_targets)
cal.forecast(prediction_targets)

# 4: build & sanity-check
print(cal.settings())
cal.build("calib.pst", noptmax=0)

# 5: history match
ies = cal.run_ies(reals=50, iterations=3, bad_phi_sigma=1.5, workers=10)

# 6: assess (the DBTL checklist)
ies.plot_phi(); ies.plot_vs_obs()
ies.forecasts(); ies.forecast("spring").plot()
ies.plot_field("k", stat="change")
ies.report("review.html")

# 7: learn, revise, repeat
```

Every plot above takes `backend="matplotlib"` if you prefer static
matplotlib/seaborn output over interactive Plotly.
