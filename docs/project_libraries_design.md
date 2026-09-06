# Project spec libraries — design

*One grammar over four libraries; run records that say what happened; forks that go
through the front door.*

**Scope: inputs and model structure only.** Loading a run, opening its model and
reaching results already works and is untouched by this.

Status: **BUILT** — all four steps landed 2026-09-04/05; see the per-step notes
and ledger 179-187. Kept as the reasoning record.

Originally: **design, not built.** Every defect cited below was reproduced against
`/home/lukem/python/simple_modflow` @ `myflopy` and the real project at
`/home/lukem/models/mf6/tentrails` on 2026-09-03.

---

## 1. Why this, and why in this order

The premise is already proven by the user's own project. Every package in every run
carries a library key:

```
baseline      dis/base npf/base sto/base ic/steady chd/base ghb/base drn/base
              riv/rock crisp/base oc/all rch/series
mounding      …                                     rch/background rch_pond
mounding_vol  …                                     rch/background rch_pond_vol
```

Three scenarios differing in exactly the recharge key. The discipline this design is
meant to support is one the user already practises by hand. What is missing is that
**none of it survives the session, and the run records are wrong.**

| # | Defect | Site | Measured |
|---|---|---|---|
| D1 | Reopening a run **erases** its provenance | `workspace.py:665`, `:634` | `Run.reopen` never sets `spec=`, so `to_dict()` writes `"spec": null` over the 11-package block. `_record_status` (`:657`) saves on every transition, so `reopen_run("baseline").execute()` destroys the record of what `baseline` was. Reproduced on a copy of the real manifest. |
| D2 | `status`/`success` cannot tell success from failure | `workspace.py:657` | All four runs: `status="written"`, `success=null`. Three reached `Normal termination`; `from_usg` did not. |
| D3 | The library does not persist | `workspace.py:1105` | `specs/packages/` is empty and `specs/project_spec.json` absent, so `Project.load` raises. The 11–12 keys per run are unresolvable in any later session. |
| D4 | `package_spec(key)` was dropped | git `10d7a87` | The retired `ProjectSpec` had it; the 2026-06-17 fold kept `simulation(name)` and dropped the package accessor. Restoring it undoes a regression. |
| D5 | `prepare_run(overwrite=True)` clobbers the record | `workspace.py:848` | Resets `status`/`success`/`metadata` while leaving the previous run's outputs in place — including orphan files only the old spec wrote. |

Build order follows severity, not API tidiness: **truth → persistence → discovery →
fork.** Fork is the feature that was asked for; it is last because a fork you cannot
reload and cannot trace is a worse version of writing `mf.rch(...)` twice.

---

## 2. The whole thing in use

```python
import myflopy as mf

# ── session 1 ────────────────────────────────────────  build the library once ──
proj = mf.Project("/home/lukem/models/mf6/tentrails", name="tentrails")

proj.add_grid("grid/base", vor)
proj.add_package("dis/base",   disv)
proj.add_package("npf/base",   mf.npf(k=res.per_layer(K), icelltype=1))
proj.add_package("sto/base",   sto);    proj.add_package("ic/steady", ic)
proj.add_package("chd/base",   chd);    proj.add_package("ghb/base",  ghb)
proj.add_package("drn/base",   drn);    proj.add_package("riv/rock",  riv)
proj.add_package("crisp/base", crisp);  proj.add_package("oc/all",    oc)
proj.add_package("rch/series",     mf.rch.array(context=ctx, recharge=series))
proj.add_package("rch/background", mf.rch.array(context=ctx, recharge={0: mean_rch}))
proj.add_package("rch_pond",       mf.rch.array(context=ctx_pond, name="rch_pond",
                                                recharge=pond_flow))

# fork the recharge package, not the whole model
proj.packages.fork("rch_pond", "rch_pond_vol",
                   name="rch_pond_vol", recharge=pond_flow_vol,
                   note="max-volume pond")

proj.save()
# SaveReport(project_spec='specs/project_spec.json', written=6, pickled=8, skipped=0)

# ── session 2, next morning ──────────────  nothing rebuilt, nothing remembered ──
proj = mf.Project.load("/home/lukem/models/mf6/tentrails")
proj
# Project('tentrails', root='…/tentrails', packages=14, grids=1, simulations=3)

proj.packages.names()
# ['chd/base','crisp/base','dis/base','drn/base','ghb/base','ic/steady','npf/base',
#  'oc/all','rch/background','rch/series','rch_pond','rch_pond_vol','riv/rock','sto/base']

proj.package("rch_pond_vol")
# PackageSpec(name='rch_pond_vol', builder='ModflowGwfrcha', status='enabled',
#             options={'recharge': '<ndarray (72, 9405)>', 'irch': 'top_active'},
#             metadata_keys=('forked_from',))

proj.runs.summary()
#  key           status     outcome     manifest  simulation    packages  updated
#  baseline      written    normal      behind    baseline      11        2026-09-03T…
#  mounding      written    normal      behind    mounding      12        2026-09-03T…
#  mounding_vol  written    normal      behind    mounding_vol  12        2026-09-03T…
#  from_usg      written    incomplete  behind    tentrails     11        2026-09-02T…

proj.runs["mounding_vol"].record.package_keys["tentrails"]
# ('dis/base','npf/base','sto/base','ic/steady','chd/base','ghb/base','drn/base',
#  'riv/rock','crisp/base','oc/all','rch/background','rch_pond_vol')
```

That `runs.summary()` is the whole argument for Step 1: `status` is what *myflopy*
did, `outcome` is what *MF6* did, and today only the first exists — which is why all
four of your runs look identical and none of them looks finished.

---

## 3. Step 1 — tracking truth

**BUILT 2026-09-04** — `265b897` (provenance), `872ff7d` (recipe), ledger 184 (outcome).

> **`status` is myflopy's verb. `outcome` is MF6's verdict. Two axes; never collapse
> them.** `status="written"` is not a lie; it reads as one because nothing reports
> what MF6 did afterwards.

New leaf module `src/myflopy/run_record.py` — **stdlib only, layer 0**. Not under
`myflopy/project/`: importing any `myflopy.project.*` submodule executes
`myflopy/project/__init__.py` (layer 13) at runtime. Zero new deferred imports, so
`tests/deferred_import_allowlist.json`'s `"myflopy.workspace": 1` is untouched.

```python
@dataclass(frozen=True, slots=True)
class RunRecord:
    """A run's manifest read back as a typed object — the authority on what a run
    WAS whenever its live SimulationSpec is gone. Nothing here is inferred."""

    name: str
    status: str
    executable: str
    created_at: str
    updated_at: str
    completed_at: str | None
    success: bool | None
    spec_block: dict[str, Any] | None          # the _spec_summary payload, verbatim
    metadata: dict[str, Any]

    @property
    def package_keys(self) -> dict[str, tuple[str, ...]]: ...   # model -> library keys
    @property
    def spec_name(self) -> str | None: ...


@dataclass(frozen=True, slots=True)
class ListingOutcome:
    """What mfsim.lst says MF6 did, and whether it still describes the inputs."""

    outcome: str                  # "normal" | "failed" | "incomplete" | "absent"
    finished_at: str | None
    elapsed: str | None
    manifest: str                 # "current" | "behind" | "rewritten"


def read_listing_outcome(workspace: Path) -> ListingOutcome: ...
```

**Three changes in `workspace.py`:**

1. **`Run.to_dict` must stop destroying blocks it did not author.** **BUILT — `265b897`.** Today it writes
   `"spec": None if self.spec is None`, which is how D1 erases provenance. It must
   carry forward *every* block it did not compute this call — `spec`, and any block a
   later step adds. This is the fix that must land first and alone; it is a live
   data-loss bug reachable with the code as it stands.

2. **The manifest stores the rebuildable recipe.** *(Revised 2026-09-04 — the
   original plan resolved the parent spec by name to avoid a format change; that
   constraint was lifted, and storing the recipe is strictly better.)* `run.json`
   gains `spec_recipe` = `SimulationSpec.to_dict()` beside the readable `spec`
   summary. **348 → 3,277 bytes** on a Ten Trails-shaped simulation, because the
   packages are `mf.ref(...)` and the arrays live in the library. `reopen_run`
   attaches the project's build context so those refs resolve; `Run.reopen` and
   `mf.load_run` do not, having no library in scope. A spec embedding raw arrays
   gets no recipe and stays summary-only. **BUILT — ledger 181.**

3. **`Run.outcome` reads `mfsim.lst`, view-only.** **BUILT — ledger 184.** It is
   never cached into the manifest at all: the listing is the authority and sits in
   the same directory, so a stored verdict would be a second source of truth that
   can disagree with its own source. `run.reconcile()` adopts it into
   `status`/`success` and persists those (which ARE manifest fields, and were
   wrong).

**Not fixed here:** `prepare_run(overwrite=True)` (D5). It needs an archive-or-refuse
policy and `examples/mf6/variant_workflow/01_package_library_and_variants.py` depends
on today's behaviour. Recorded in §7.

---

## 4. Step 2 — persistence

`_save_package` **already** pickles a package whose `to_dict()` raises
(`workspace.py:1073`), so the package library already survives arrays. The blocker is
narrower than it looks: `Project.save()` calls `validate()` and raises on *any* issue
(`workspace.py:1110`), and `validate()` reports an inline array on a `SimulationSpec`
as an issue. One un-persistable simulation therefore refuses the entire project —
which is why `specs/` is empty on a project whose 14 packages would each save fine.

```python
@dataclass(frozen=True, slots=True)
class SaveReport:
    """What save() actually persisted, and what it could not."""

    project_spec: Path
    simulations: int = 0
    packages: int = 0
    grids: int = 0
    pickled: int = 0                            # of packages+grids, how many took pickle
    skipped: tuple[tuple[str, str], ...] = ()   # (simulation name, reason)

    def __fspath__(self) -> str: ...
    def __bool__(self) -> bool: ...


class Project:
    def save(self, *, strict: bool = False) -> SaveReport: ...
```

**BUILT 2026-09-04** (ledger 180). One refinement found in implementation: the two
kinds of failure want different treatment. An embedded array is nothing the caller
can fix at save time — skip and log. A **dangling `mf.ref`** is a mistake with a
fix, so it still raises — but *after* the write, so the error costs a message
rather than the library. `test_project_save_rejects_unresolved_references` passes
unchanged, which is the evidence the contract was kept rather than weakened.

`strict=True` keeps today's refuse-everything behaviour. The default persists what it
can and names what it could not, because a project that saves 13 of 14 packages and
tells you which one it dropped is strictly more useful than one that saves nothing.

**A skipped simulation must also be removed from the key list `_project_dict` writes**
(`workspace.py:1040`), or `Project.load` raises `FileNotFoundError` on the next
session — that fail-loud contract is pinned by
`test_project_load_fails_loud_on_missing_file`.

**`save()` changing return type is a breaking change.** `tests/test_project_spec.py:468`
does `saved_path = project.save()` then compares it to a `Path`; `__fspath__` does not
participate in `==`, so that assertion must be rewritten. Recorded in §7.

---

## 5. Step 3 — discovery and inspection

**BUILT 2026-09-04** (ledger 182), except `RunLibrary` — see the note at the end
of this section.

One base, four subclasses, one override point. Four hand-written near-identical
classes is the cost CLAUDE.md §4.7.6 measured and rejected.

New module `src/myflopy/library.py` (imports `myflopy.specs` + stdlib + pandas → L2),
imported at module level from `workspace.py`. `RunLibrary` needs `Run` (L10) so it is
defined in `workspace.py` and imports the base downward — forced by layering, not
taste.

```python
_T = TypeVar("_T")

class SpecLibrary(dict[str, _T], Generic[_T]):
    """A named library of specs: a dict, plus the verbs every library answers.

    Subclasses override _row and declare _kind/_columns. Nothing else. All dict
    behaviour is inherited unchanged — in, len, iteration, .get, .items, .pop.
    """

    _kind: ClassVar[str] = "entry"
    _columns: ClassVar[tuple[str, ...]] = ("key",)
    _project_name: str | None = None      # class-level default; never a ctor argument

    def _row(self, key: str, entry: _T) -> dict[str, object]:
        """One summary() row. THE single override point."""

    def names(self) -> list[str]: ...
    def summary(self) -> pd.DataFrame: ...
    def copy(self) -> Self: ...            # dict.copy returns a plain dict — measured
    def __or__(self, other: Mapping[str, _T]) -> Self: ...
    def __missing__(self, key: str) -> NoReturn: ...
    def __repr__(self) -> str: ...
    def _repr_html_(self) -> str: ...


class PackageLibrary(SpecLibrary[PackageSpec]):
    _kind = "package"
    _columns = ("key", "package", "builder", "options", "inputs", "requires", "enabled")

    def fork(self, source: str, key: str, *, name: str | None = None,
             deep: bool = True, note: str | None = None, **overrides: Any) -> PackageSpec: ...

class GridLibrary(SpecLibrary["GridSpec | GridRef"]):
    _kind = "grid"
    _columns = ("key", "grid", "type", "method", "crs", "built", "cells")

class SimulationLibrary(SpecLibrary[SimulationSpec]):
    _kind = "simulation"
    _columns = ("name", "models", "model_packages", "sim_packages", "exchanges",
                "derived_from", "lineage")
```

**`_project_name` is a class-level default, not a constructor argument.** Measured:
`dict.copy()` and `dict.__or__` return a **plain `dict`**, silently dropping both the
subclass and any instance state. A required `__init__` argument would make the library
degrade to a bare dict on any of those. `copy`/`__or__` are overridden to preserve the
type, and the name is set by `Project` after construction — absent, the `KeyError`
message just omits the project name.

```python
proj.packages['npf/bass']
# KeyError: Project 'tentrails' has no package 'npf/bass'.
#           Did you mean: 'npf/base'?
#           14 packages: chd/base, crisp/base, dis/base, drn/base, ghb/base,
#           ic/steady, npf/base, oc/all, rch/background, rch/series, rch_pond,
#           rch_pond_vol, riv/rock, sto/base
```

`__missing__` is the hook rather than a new `require()` verb: it upgrades the message
on the spelling everyone already types, leaves `.get()` returning `None`
(`dict.get` never calls it), and adds no name to the grammar.

**There is no `persists` column.** It was in the first draft and measured **119 ms per
100 000-record list BC** — Jupyter calls `_repr_html_` on every bare `proj`, and the
user's USG-derived list EVT is ~650 000 records. It is an explicit call instead:

```python
    def persistence(self) -> pd.DataFrame: ...   # key, persists ("json"/"pickle"), reason
```

**Registration keeps its existing front doors** — `add_package` / `add_grid` /
`add_simulation`. No `.add()` on the libraries: a second spelling of an operation that
already has one, and the fourth library could not have it anyway (a run is created by
`prepare_run`, never added). `project.packages[key] = spec` keeps working (D1).

**Keyed accessors**, restoring what the fold dropped:

```python
    def package(self, key: str) -> PackageSpec: ...        # = self.packages[key]
    def grid(self, key: str) -> GridSpec | GridRef: ...    # = self.grids[key]
    def simulation(self, name: str) -> SimulationSpec: ... # exists already
```

There is **no `Project.run(name)` accessor** — `Project.run(name, simulation, …)`
already exists and *executes*. The run accessor is `project.runs[name]`, full stop.

**`__repr__` touches no disk.** `len()` on three in-memory dicts, and runs are not
counted. PyCharm calls `__repr__` on every debugger step, and `Run.reopen` does an
unguarded `json.loads(...)["run"]`, so one truncated manifest would otherwise make
`repr(project)` raise.

### The runs library

```python
class RunLibrary(SpecLibrary[Run]):
    """Every run under a project — with a manifest, without one, or broken."""

    _kind = "run"
    _columns = ("key", "status", "outcome", "manifest", "simulation",
                "models", "packages", "updated", "path")

    def refresh(self, *, load: bool = False) -> Self: ...
```

| condition | key | status |
|---|---|---|
| dir under `runs_dir` with `run.json` | dir name | the manifest's |
| dir under `root` with `mfsim.nam`, no `run.json` | path relative to `root` | `discovered` |
| `run.json` unparseable or missing `"run"` | dir name | `unreadable` |

Manifest-less runs are included deliberately: there are **39 MF6 workspaces against 14
`run.json` files** under `examples/`, and most notebooks call `sim.run_simulation`
directly. A ledger that silently omitted two thirds of real runs would be worse than
none. Path-keying them keeps `runs/from_usg` and `scratch/from_usg` distinct.
`unreadable` is a **row, never an exception** — `json.JSONDecodeError`, `KeyError`
and `OSError` caught by name, so the broad-catch allowlist does not grow.

Scan cost measured: `rglob("mfsim.nam")` is 0 ms on the real project (4 hits), 2 ms on
`examples/` (39 hits). Cached, with an explicit `refresh()`.

**`RunLibrary` — BUILT 2026-09-04 (ledger 185).** Original note kept below for
the reasoning that deferred it; the outcome reader it waited on landed in 184.

~~**`RunLibrary` is NOT built yet.**~~ Its two columns worth having — `outcome` and
`manifest` — both come from `read_listing_outcome`, the other half of Step 1 that
has not landed. A runs table whose only column is the `status` field that already
reads `written` for every one of your runs would add a surface without adding an
answer. It follows the outcome reader.

---

## 6. Step 4 — fork

**BUILT 2026-09-05** — `0dab9fd` (fork/edit), `043f55c` (rename carry), ledger 186.

> **A fork is built, never patched.** A forked spec leaves through the same
> constructor an ordinary `mf.rch(...)` leaves through, so every check that guards a
> hand-written package guards a forked one.

The obvious implementation is silently wrong. `PackageSpec` is
`@dataclass(frozen=True, slots=True)` and validates in `__post_init__` →
`_check_options`. `copy.deepcopy` reconstructs a slots dataclass through
`__reduce_ex__`, so **`__post_init__` never runs**. Measured:

```
deepcopy + in-place option write   -> typo 'kk_typo' ACCEPTED, surfaces at MF6 build
replace(deep-copied options)       -> TypeError: rch_pond: ModflowGwfrcha does not
                                      accept 'kk_typo'. Did you mean 'recharge'?
```

The mechanism is one line: deep-copy the option **values** into a fresh dict, then
hand it to `dataclasses.replace`. Never `object.__setattr__` on a copy.

```python
def fork_package(package: PackageSpec, /, *, name: str | None = None,
                 deep: bool = True, **overrides: Any) -> PackageSpec:
    options = copy.deepcopy(package.options) if deep else dict(package.options)
    return replace(package, options={**options, **overrides},
                   name=name or package.name)
```

Measured on this: the fork's arrays are **not** aliased to the parent's, including for
options the fork never touched — which is where the real hazard lives, since
`with_options` today rebuilds the dict but not its values and lets a "copy" write
through into the baseline's `k`.

`deep=True` costs **15.2 ms / 0.50 MB** on a Ten Trails-scale `SimulationSpec`
(9 405 cells, 72 periods, ~14 400 CHD records). `deep=False` is the opt-out for
someone who wants the aliasing on purpose.

**CORRECTED 2026-09-06 (ledger 188/189):** that figure is a small-model one and
misleads about the shape. Cost tracks the number of Python objects copied, not
cells — a 9 405-cell `npf` forks in 0.1 ms and a 72×9 405 array in 1.7 ms, while
the real 72-period × 9 090-record EVT takes **808 ms**.

Every fork stamps `metadata["forked_from"] = {"key", "at", "note"}`, riding the
existing `PackageSpec.metadata` slot that `_spec_summary` already copies into run
manifests. **No new dataclass field.**

**Simulation forking wraps `SimulationSpec.derive`**, which already exists with real
`derived_from`/`lineage`. `swap=` must be model-qualified — `{"model": {"old": "new"}}`
as well as the flat form — because a coupled GWF/GWT simulation has `ic`, `oc` and
`dis` on both models and a flat key matches two entries with no way to disambiguate.
And the fork must register its parent, or `Project.validate` refuses to save a project
whose `derived_from` names an unregistered simulation.

---

## 6a. Editing a package that a run already used

**BUILT 2026-09-05** — `3987866`, ledger 187. One correction: digests are
computed at `prepare_run`, not in `_spec_summary` as designed below — that
function sees `_package_entry_label` only, and a `PackageRef` is stored
unresolved, so digesting it needs the project library a module function cannot
reach.

Fork answers "I want a variant". It does not answer "`npf/base` was wrong, fix it" —
and that is the case with teeth, because **the run manifest records package *keys*,
not values.** Editing `npf/base` after `baseline` ran leaves the manifest saying
`baseline` used `npf/base`, while `npf/base` now means something else. Measured on the
current code: two models differing **50× in K** produce byte-identical manifests, and
`ModelGroup(...).diff().report()` calls them "identical to reference".

So editing is not the hard part. Knowing what an edit invalidated is.

### The edit verb

```python
class PackageLibrary(SpecLibrary[PackageSpec]):
    def edit(self, key: str, *, name: str | None = None, deep: bool = True,
             note: str | None = None, **overrides: Any) -> PackageSpec:
        """Replace the package at ``key`` with an edited copy of itself.

        Same key, same mechanism as :meth:`fork` — the edited spec is BUILT through
        ``dataclasses.replace``, so a typo raises here exactly as it would in
        ``mf.npf(...)``. Mutates the library (D5) and returns the new spec.
        """
```

`edit(key, …)` and `fork(source, key, …)` are the same operation with different
destinations, deliberately spelled apart: one keeps the key and changes what it means,
the other leaves the original alone. Plain `project.packages[key] = spec` also still
works — it is a dict — but skips the drift stamp below, so `edit` is the taught
spelling.

```python
proj.packages.edit("npf/base", k=res.per_layer(K_revised), note="post-calibration K")
# PackageSpec(name='npf', builder='ModflowGwfnpf', …)
```

### How it is tracked

One new field in `_spec_summary`, written per model alongside the keys it already
records:

```python
"packages":        ["dis/base", "npf/base", …],              # exists today
"package_digests": {"dis/base": "9f2c…", "npf/base": "24e5…"} # new
```

The digest is a content fingerprint of what the package *builds* — name, builder
reference, options, inputs, enabled, requires — **not** the pickle of the spec object.
That distinction is load-bearing: a `functools.partial` builder has no `__eq__`, so
hashing the object graph gives a different answer after a JSON round trip even when
every field compares equal. Measured on the field-based digest:

| property | result |
|---|---|
| stable across processes and `PYTHONHASHSEED` | yes |
| stable across `to_dict`/`from_dict` | yes |
| stable across `to_yaml`/`from_yaml` | yes |
| stable across `Project.save()`/`load()` | yes |
| distinguishes `k=1.0` from `k=50.0` | yes |
| distinguishes one cell changed by `1e-9` in a 9 405-cell array | yes |
| cost | **0.04 ms** scalar, ~5 ms for a 72×9 405 array, **1.3 s** for the 654 480-record EVT — hence the cache (ledger 189) |

That turns the question into a column:

```python
proj.runs.summary()
#  key           status     outcome  drift
#  baseline      completed  normal   npf/base
#  mounding      completed  normal   npf/base
#  mounding_vol  completed  normal   —

proj.runs["baseline"].drift()
#  key        ran_with  library_now  state
#  npf/base   24e5b755  8a1f0c93     edited
#  rch/pond   4bb359e0  —            removed
#  dis/base   9f2c1d40  9f2c1d40     unchanged
```

Three states, and the third matters as much as the first: `edited` (the key still
exists and means something else), `removed` (the key is gone from the library), and
`unchanged`. A run whose every package reads `unchanged` is reproducible from the
library as it stands; anything else tells you exactly which packages to look at.

**Editing is never refused.** It is your library, and iterating on `npf/base` is the
normal case — the point is that the consequence becomes visible rather than silent.
`edit` stamps `metadata["edited"] = {"at", "note"}` the way `fork` stamps
`forked_from`, so the spec itself also carries the trace.

**Legacy manifests have no digests**, including all four existing runs. Those rows
report `drift = unknown` rather than `—`; a missing fingerprint is not evidence of
sameness, and claiming otherwise would be the same class of lie as `success=None`
reading as "fine".

---

## 7. Explicitly not building

1. **`rename` beyond `pname`/`filename`.** A package name is embedded in **four**
   places, not three: `options["pname"]`, `options["filename"]`, the uzf/lak/sfr
   `*_filerecord` (`package_registry.py`), and `mvr_spec`'s `budget_filerecord` /
   `budgetcsv_filerecord` (`advanced.py:631`) — whose template
   (`{model_name}.{name}.bud`) is a different convention from the registry's
   (`{model_name}_{name}_budget.uzf`), and `lak`'s `stage_filerecord` is a third.
   `PackageSpec` carries **no package-type field**, so nothing can dispatch between
   them. Renaming also invalidates sibling `requires=` tuples, MVR
   `options["packages"]` and every `perioddata` source/destination string, and
   `MoverConnection.package`. `fork(name=…)` therefore rewrites `pname`/`filename`
   only and **raises** when the old name appears anywhere else.
2. **`GridLibrary.fork`.** `GridSpec.voronoi` routes 17 named parameters into
   `options` then dispatches `**engine_options` through an explicit key-set table into
   `triangle`/`mesh`/`voronoi_options`, with `min_angle` aliased to `angle`.
   Reproducing that is a second hand-written per-key site — the §4.7.6 hazard — and
   `GridSpec.obj` holds a *live mutable* built mesh besides. The supported moves stay
   `mf.grid_ref(key)` and `add_grid(new_key, …)`.
3. **`prepare_run(overwrite=True)` history loss (D5).** Needs an archive-or-refuse
   policy; the flagship variant example depends on today's behaviour.
4. **Attribute completion of library keys.** Keys contain `/` (`npf/base`), so they
   are not identifiers. No `.pyi`, no `__getattr__` tier. Keys are answered at runtime
   by `names()`, `summary()`, the repr, and the `__missing__` message.

**Breaking changes to record:** `save()` returns `SaveReport` rather than `Path`
(`tests/test_project_spec.py:468` must be rewritten); `project.packages` becoming a
read-only property makes `project.packages = {...}` raise (no call site does this —
grepped — but it is a change).

---

## 8. Implementation checklist

| Area | Work |
|---|---|
| `src/myflopy/run_record.py` | new, ~120 lines: `RunRecord`, `ListingOutcome`, `read_listing_outcome` |
| `src/myflopy/library.py` | new, ~200 lines: `SpecLibrary` + 3 subclasses |
| `src/myflopy/workspace.py` | `Run.to_dict` carry-forward; `Run.reopen(project=)`; `Run.outcome`/`record`/`reconcile`; `RunLibrary`; `Project.runs`/`package`/`grid`/`__repr__`/`_repr_html_`; `save() -> SaveReport` |
| `src/myflopy/specs.py` | `fork_package` helper; `package_digest`; `SimulationSpec` fork support |
| `_spec_summary` (`workspace.py:330`) | add `package_digests` per model — the one change that makes drift detectable |
| Tests | new `test_run_record.py`, `test_spec_libraries.py`; edit `test_workspace.py` (`save()` return type, reopen spec), `test_project_spec.py:468` |
| Pinned/derived | `scripts/derive_import_layers.py` + `docs/import_layering.md` (two new modules; the doc carries the ratchet count and must be regenerated with the JSON). **`tests/api_snapshot.json` does NOT cover `Project`/`Run`** — verified, it iterates `package_api.__all__` only — so either add a guard for this surface or record that it is unguarded |
| Docs (hard rule) | `docs/core_specs.md:436-453`, `docs/manual/03_getting_started.md:121,303,413`, `docs/manual/04_architecture.md`, `docs/manual/README.md`, `docs/myflopy_context.md:103,135`, `docs/codebase_structure.md:51-54`, `docs/package_api_reference.md`, `docs/model_building_cheatsheet.md` |
| Examples | `examples/mf6/package_first_full_stack.py`, `contours_to_model.py:303`, `demo_site/build_model.py`, `variant_workflow/01`, `variant_workflow/03_coupled_gwf_gwt.py:35` |
| Ledger | entries for each item in §7 plus the two breaking changes |

Steps 1 and 2 are independently valuable and land first; 3 and 4 are the API the
request was about, and are worth little without them.
