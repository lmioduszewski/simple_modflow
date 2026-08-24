# View-layer conventions (`model.packages.…`)

> **Status:** normative. Added 2026-07-18 after `sfr.results.long_profile` /
> `plot_long_profile` drifted outside the grammar and two notebooks
> hand-rolled ~20 lines of matplotlib to redraw a figure the library already
> built — losing `scrollZoom`, pan, the house template, *and* the documented
> gaining/losing colors in the process. Extended 2026-08-20 (plan 8.7) with the
> three-layer model Phase 8 built.

## The three layers

Everything that draws in myflopy sits in one of three layers. Knowing which
layer a thing belongs to answers most design questions about it.

| Layer | What it is | Verbs |
|---|---|---|
| **1 — Pictures** | one drawing of one subject | `map` `section` `surface` `grid` |
| **2 — Composition** | combines finished pictures | `mosaic` `animate` |
| **3 — Output** | what you do with any picture | `.fig` `.show()` `.save(p)` `.html(p)` |

Three rules follow from the table, and they are the whole design:

1. **Geometry chooses the Layer-1 verb — not content, and not renderer.** A map
   is a plan view whatever is drawn on it, so contours, well markers, a
   hillshade and particle pathlines are *options* on `map`, never verbs. This is
   what retired `plot3d` (a 3-D view is `surface`), `map_nodes` (a map whose
   values are node ids) and a top-level `contours`.
2. **Layer 2 takes a collection, Layer 1 takes a subject.** `map(model)` draws
   the model; `mosaic(panels)` and `animate(frames)` draw what you hand them.
   That is why the combinators do not fit the "first argument is the subject"
   rule the picture verbs follow — they have no subject.
3. **Layer 3 is uniform.** Every picture answers the same four things, so you
   never learn per-class output methods. There is no trailing `.plot()`.

### `backend=` switches the renderer, never the subject

`grid(backend="vtk")` draws the same mesh in 3-D. `animate(backend="png")`
rasterizes the same frames. Both are legal because the *subject* is unchanged.

The counter-example is instructive: plan 8.5 was written as
`surface(backend="vtk")`, and that would have been a lie — `surface` means a
height field `z(x, y)`, while the VTK scene draws a cell **volume** and the
particle scene draws **polyline tubes**. Three shapes, not one shape three ways.
Hence `grid(backend="vtk")`: `grid` already means "the mesh itself", so the
3-D layered mesh is that same subject, redrawn. See ledger 133.

### Three renderers, one contract

Most pictures are Plotly, and `Picture` is written in terms of `.fig`. That is a
*default*, not the contract. A picture whose native renderer is something else
answers the same four verbs by overriding them, and its `.fig` raises and names
what to use instead:

| Class | Native renderer | Instead of `.fig` |
|---|---|---|
| `MplPicture` | Matplotlib | `.axes` |
| `VtkScene` | PyVista | `.scene` |
| `SliderAnimation` | rasterized frames | `.frames` |

Returning an `mpl.Figure` from `.fig` would satisfy the letter and break every
caller reaching for `.add_trace`; raising is the honest answer.

### Two spellings, one implementation

```python
plot.map(model, layer=0)     # free function, from `myflopy.plot`
model.plot.map(layer=0)      # bound to the object
```

The bound form *calls* the free one, so they cannot diverge. Which verbs a scope
answers depends on what that scope can know: a bare grid has no results, so no
`surface` or `animate`; a layer stack has no time, so no `animate`.
`tests/test_plot_vocabulary.py` pins each scope's set, and fails naming both the
scope and the stray verb if one drifts.

### Write the parameters into the `def`, not only the docstring

**Never give a picture verb a bare `**kwargs`.** PyCharm and Pylance are static
— they read the `def` line and never execute the module — so a forwarder shows
the caller nothing no matter how thorough its docstring. `__signature__` and
`functools.wraps` reach `help()` and neither reaches an editor, and a `.pyi` stub
is ruled out here for the reason in CLAUDE.md §4.7.7 (a stub replaces the whole
module for type checkers).

So each verb names its parameters, with the exact default of whichever link in
its forwarding chain owns that argument. The bound form repeats them. Keep a
`**kwargs` tail **only where the tail is genuinely open** — `map`'s reaches the
`go.Choroplethmap` trace, whose names Plotly owns.

Two rules follow, both enforced:

* A narrower scope may offer **fewer** parameters, never other ones. `vor.plot.map`
  drops `per`/`layer`/`type` — a bare grid has no results — and what it drops
  still rides the tail, so narrowing a signature never narrows what worked.
* A parameter you document must be one you accept. `plot.grid` documented
  `layers`/`scale`/`color_by`/`cmap` for three releases; those live on the
  *layer-stack* builder and were never reachable through it.

The exception is a verb that is genuinely polymorphic. The namespace-level
`field=` sugar dispatches to accessors with incompatible signatures, so it keeps
`(*args, field=None, **kwargs)` and its docstring says to call the leaf instead.
That is a real property of the design, not a shortcut.

## Why this file exists

The **spatial** half of the grammar was already decided and documented
(`docs/myflopy_context.md`: *"`map/plot/section` + `mosaic/animate` on every leaf"*,
implemented by `SpatialView` / `FieldMappable` in `package_plotting.py`).

What was never written down is the other half: what to do with a **derived
table** — a purpose-built frame like the SFR long profile that merges several
fields and is not a single mappable field. Because no rule covered them, they
became loose method pairs on the namespace (`long_profile()` returning a frame,
`plot_long_profile()` returning a figure). Five different verb spellings
accumulated (`map`, `plot`, `plot_profile`, `plot_budget`, `plot_long_profile`),
none of them wrong, none of them predictable.

This file states both halves as one rule.

## The shape

```
model.packages.<pkg>.<inputs|results>.<noun>.<verb>()
```

- **`<pkg>`** — the package: `sfr`, `lak`, `chd`, `riv`, `evt`, …
- **`<inputs|results>`** — declared inputs, or simulated results.
- **`<noun>`** — *what you are looking at*. Either a **field** (`q`, `stage`,
  `head`) or a **derived table** (`profile`). Always a noun, never a verb, and
  never prefixed with the verb that consumes it.
- **`<verb>`** — *what you want out of it*, from the fixed set below.

Every noun is an object. Every object answers the same verbs. That is the whole
rule; it is what makes the API guessable without reading source.

## The verbs

| Verb | Returns | On |
|---|---|---|
| `get(...)` | `DataFrame` — the normalized rows | every noun |
| `summary()` | `DataFrame` — a compact digest | every noun |
| `plot(...)` | `viz.Fig` — the node's NON-SPATIAL chart | every noun |
| `map(...)` | choropleth | spatial nouns |
| `section(...)` | cross-section | spatial nouns |
| `mosaic(...)` / `animate(...)` | multi-panel / animated | spatial nouns |

These are the same Layer-1 and Layer-2 verbs as the table at the top of this
file, scoped to one noun. `surface` and `grid` are not in the grammar because a
package field has no 3-D contact and no mesh of its own — those belong to the
model, the grid, and the layer stack (`model.plot`, `vor.plot`, `stack.plot`).
Everything returned obeys the Layer-3 contract, so `…​.q.map(per=0).save("q.png")`
works without knowing which class produced it.

> **What `plot` means, precisely.** Not "time series" -- the codebase never
> honoured that, and saying so misled a reader as recently as 2026-08-18.
> `plot()` is *the chart of this node that is not a map and not a section*, and
> its shape follows the node: a series by stress period on a field
> (`SpatialView`, `FieldMappable`, `LakStageChangeExplorer`), a **distance
> profile** on the SFR profile views, a **bar chart** on `LakBudgetView` /
> `PRTEndpointsView` / `PRTCaptureView`, a **cumulative arrival curve** on
> `PRTTravelTimeView`, a **histogram** on `IesForecast`. One verb, one question
> -- "chart this node" -- and the node decides what chart answers it. That is
> what keeps the grammar guessable; a verb that changed name with the chart type
> (`timeseries`/`profile`/`bars`) would not.
>
> A corollary worth stating: if a `plot()` returns a MAP, it is misnamed. Three
> did -- `RchInput`/`UzfInput`/`DrnInput` returned choropleths -- and are now
> `map()` (plan 8.2).

The model-level **dependent-variable readers** follow the same spatial-noun verb
set: `model.hds` (GWF heads), `model.conc` (GWT concentration), `model.temp` (GWE
temperature) all share `get/summary/array/map/section/mosaic/animate`, built from one
`DependentVariableFile` base (`headsplus.py`) — add a new dependent variable by
subclassing it and setting `value_name`/`store_column`/`_binary_text`/`_choro_type`,
not by cloning the reader. The choropleth reads whichever field the map's `type`
selects (`_DEPVAR_READER_ATTR` in `choros.py`), so the per-layer hover table and
value labels stay field-generic.

Results that are not per-cell fields to begin with join the grammar the same way,
by *deriving* one: a finished PRT run exposes `results.travel_time`,
`results.endpoints`, and `results.capture` (`prt_maps.py`) rather than
`travel_time_map()` / `endpoints_map()` / `capture_map()`. Those maps are
**time-integrated** over the run, which is the one documented departure: they take
no `per=`, carry no period footer (`result_hover(..., footer=())`), and `plot()`
draws a distribution across particles instead of a series by stress period (see
the compromise ledger, entry 59).

A noun does **not** have to be a per-cell field at all: `results.pathlines` keeps
the trajectories, and its `map()` draws one polyline per particle over a base map
instead of coloring cells. It still answers the same verbs (`get` = the
normalized records, `plot` = elevation vs travel time, `mosaic` = one panel per
release group) and refuses the ones that make no sense (`section`, `animate`) with a
message naming the alternative. Overlay traces reach a mosaic through
`Choro.add_overlay` / `overlay_traces()`; a map panel that hand-adds traces to
`choro.fig` will have them dropped when composed.

**Categories are colored by policy too.** Named, unordered things — release
groups, zones, scenarios — take `viz.category_colors(names)`
(`PALETTE.categorical`, colorblind-safe), which *memoizes* name → color so a
group keeps one color across every figure it appears in. Never index the palette
at the call site, and never fall back to plotly's default colorway: that is what
let one release group be blue on its map and orange on its curve.

Rules that hold for all of them:

- **Figures are always `viz.Fig`**, never a bare `go.Figure`. `Fig` carries
  `dragmode='pan'`, `scrollZoom`, and the house template. Returning a bare
  figure silently drops all three — it is the defect this convention exists to
  prevent, and it has now happened twice (`hds.animate`, and these notebooks).
- **Colors follow the policy**, not the call site. Signed exchange fields use
  `_exchange_colorscale(frame)`, which puts BLUE on whichever end is gaining for
  the declared reference frame — the negative end for `gwf` (SFR + list BCs),
  the positive end for `feature` (LAK, and the normalized `exchange_intensity`).
  See "Signed exchange columns name their reference frame" below and
  `docs/mf6io_reference.md`. Read the endpoints off the helper rather than
  writing hex literals or hardcoding one orientation, so discrete and continuous
  renderings cannot drift apart — and so a package's colours cannot disagree
  with its declared sign.
- **A diverging scale ships as STOPS, never as a name.** `Choro`'s
  plotly→matplotlib table maps `'rdbu'` to `'RdBu_r'` — the *reversed* colormap —
  so a map built with `colorscale='RdBu'` renders mirrored between `plot()` and
  `plot_mpl()`: red means "less" in one and "more" in the other. Pass the stop
  list from a helper (`_blue_white_red_diverging_colorscale`,
  `_red_white_blue_diverging_colorscale`) and both backends agree. Also note that
  `zmid` reaches only the Plotly trace; what centers a diverging map on *both*
  backends is symmetric `zmin`/`zmax`.
- **Derived maps that no registry describes pin their scale at the source.**
  `package_registry` is keyed by package/field *name*, so a quantity it does not
  name — a PRT travel time, a *statistic* of a calibrated parameter field — keeps
  its policy as a constant or a small helper in the module that draws it
  (`prt_maps.PRT_COLORSCALE`, `ies._field_map_policy`), pinned in
  `tests/test_colorscale_policy.py`. Do not add a second keyspace to the registry
  for one consumer.
- **A mosaic pools its panels onto one color axis, which is the point — but it
  discards each panel's limits AND its colorbar.** That pooling is what makes
  small multiples comparable, so do not fight it per panel: pass
  `viz.mosaic(..., colorbar=...)`, either a dict or a `(cmin, cmax) -> dict`
  callable when the labels depend on the pooled range (a log mosaic's real-unit
  decades can only be computed once the shared limits are known — see
  `ies._log_decade_colorbar_for_mosaic`). For a *diverging* quantity also pass
  `diff=True`, or the neutral point lands wherever the pooled data happens to put
  it. Better still, refuse the combination that needs it: `plot_field_mosaic`
  accepts only stats that *have* a prior and a posterior form, which makes the
  trap unreachable rather than merely documented. Ledger 71.
- **Anything drawn OVER a static map takes its colormap from the same place the
  cells did** — `choros.mpl_colormap_for(colorscale)`. Building a second
  colormap at the call site is how an overlaid point and the cell beneath it end
  up different colors for the same value. Note the Plotly and matplotlib overlay
  paths are genuinely different: `plot_mpl` ignores `Choro` overlays entirely and
  draws in **model coordinates**, while the Plotly path needs
  `vor.points_to_latlon`. Supporting both is ~10 lines; dropping the static one
  leaves a figure that silently shows the base map and none of its data.
- **Period selection is `per=`** on the verb, and a noun may also be *called*
  to bind a period once: `results.profile(per=3).plot()` ==
  `results.profile.plot(per=3)`. Calling a noun returns a new bound view; it
  never mutates the original.

## What NOT to write

```python
results.long_profile(per=3)        # noun carrying no verb -> what does it return?
results.plot_long_profile(per=3)   # verb_noun -> a fourth spelling to memorize
results.get_profile_plot(per=3)    # ditto
```

```python
results.profile.get(per=3)         # the frame
results.profile.plot(per=3)        # the figure
```

If a new derived table needs a home, add a **view class** with `get`/`plot`
(plus `summary` where meaningful) and expose it as a noun property. Do not add
a `foo()` + `plot_foo()` pair.

## Worked example

`SfrProfileView` (`package_surface_water.py`) is the reference implementation:

```python
profile = model.packages.sfr.results.profile

profile.get()                       # merged reach table: geometry + stage + q
profile.summary()                   # compact digest
profile.plot()                      # Fig: streambed, stage, signed exchange bars
profile(per=11).plot()              # same figure, last stress period
profile.plot(signed_exchange=False) # the older unsigned single line
```

`signed_exchange=True` (the default) draws one bar per reach colored by sign,
blue where the reach gains and red where it loses, reading both colors off the
shared diverging scale.

## Naming

Prefer the **shortest noun that is unambiguous in its namespace**. Inside
`sfr.results` a profile can only be longitudinal, so the noun is `profile`, not
`long_profile` — and it matches `sfr.results.stage.profile` one level down,
where the same noun means the same thing scoped to one field:

- `sfr.results.profile` — every field along the stream
- `sfr.results.stage.profile` — one field along the stream
- `sfr.results.q.profile` — the exchange field along the stream
- `lak.results.q.budget` — the LAK exchange summarized by connection type

Namespace-level nouns merge fields; field-level nouns cover a single field.

### `budget` is spelled at two tiers, over two different files

- `model.budget.<term>` — the **model** budget file (`.cbc`), one noun per MF6 record,
  terms **discovered from the file** so they follow the model's kind and packages.
  Each term is a spatial noun with the full verb set.
- `model.packages.<pkg>.budget.<term>` — the **package-output** budget file, a
  genuinely different file whose node layout is feature-first. Its terms are
  hand-declared, and they answer a reduced verb set (`get`/`summary`/`wide`, no
  `plot`) — a known inconsistency, ledger 98.

Treating these two files as interchangeable is what produced the node off-by-one in
ledger 92. When adding anything here, name which file you are reading.

The loose field-level spellings `q.plot_profile` / `stage.plot_profile` /
`q.plot_budget` / `q.budget_summary` were retired onto these nouns in plan §4.8
(D12 warned aliases preserve their old returns exactly). One deliberate
exception: because the SFR field explorers' old data method was itself named
`profile()`, the noun `profile` replaces it directly — the frame now comes from
`q.profile.get()`, and there is no warned `profile()`-returns-a-frame alias (the
name is the noun). The `lak.results.q.budget` figure is still matplotlib pending
a separate `viz.Fig` conversion (recorded in the compromise ledger).

## Signed exchange columns name their reference frame

The **accessor** for a package's groundwater exchange is always `results.q`
(one noun, every package). The **DataFrame column** that `.get()` emits, however,
names the reference frame the sign is measured in — because myflopy keeps MF6's
raw sign and never negates data to force one convention:

| column | source record | reference frame | negative means |
|---|---|---|---|
| `q_gwf` | aquifer `.cbc` cell record (SFR + every list BC) | the GWF cell | discharge **out of the aquifer** (the feature *gains*) |
| `q_lake` | LAK package budget file | the lake | the lake **loses** |

So `model.packages.sfr.results.q.get()` has a `q_gwf` column and
`model.packages.lak.results.q.get()` has a `q_lake` column; both carry MF6's
exact sign. Which frame a package is in is declared once, on
`ResultSpec.reference_frame` (`package_registry.py`), and drives the column name,
the map's colour orientation, and the docstrings together — never a literal at a
call site. Plots always show gaining blue / losing red by consulting that frame,
so a package's colours and its declared sign cannot drift apart. Derived columns
(`q_per_length`, `q_per_area`, `q_diff`) keep their own names; the combined
`surface_water` map draws myflopy's own normalized `exchange_intensity` field
(positive = gaining).

Never resolve the SFR-vs-LAK sign difference by negating one package's data:
that diverges from every MF6 file and, once a downstream consumer flips it back,
silently inverts a map (it did, on 2026-07-19). Name the frame instead.

## Changing a noun's name

Renames go through the standard deprecation mechanism (`_deprecation.py`,
decision **D12**): the old spelling keeps working with a `DeprecationWarning`
naming its replacement, resolves **only** via `__getattr__` so it stays out of
`__all__` / `dir()` / completion, and is listed in `myflopy.__compatibility__`.

A deprecated alias must return **exactly what it returned before** — if the new
verb changed a default, pin the old default inside the alias (the retired plot
spelling does this for `signed_exchange=False`).

The retired spelling must not survive anywhere an IDE reads — not in `dir()`,
`__all__`, a docstring tooltip, or the **name of the private method backing it**
(`_legacy_profile_plot`, not `_plot_long_profile`). See
`docs/deprecation_policy.md` §2 for the full surface list and the guard test.

## Enforcement

- `tests/test_sfr_profile_view.py` pins the view shape, the signed-exchange
  colors against the shared scale, and the deprecated spellings.
- `tests/test_colorscale_policy.py` pins the colorscale policy.
- `tests/test_deprecation.py` pins D12 hiding and the compatibility registry.
- `tests/api_snapshot.json` records the public surface; regenerate with
  `scripts/derive_api_snapshot.py` and review the diff.
- `tests/test_plot_vocabulary.py` pins the verb set at every scope (module,
  model, grid, stack), that no scope invents a verb, and that every retired
  spelling stays gone. A verb added in one place and forgotten in another fails
  here, naming both.
- `tests/test_picture_contract.py` pins Layer 3 against the REAL picture
  classes, enumerated from `Picture.__subclasses__()`. It is written that way on
  purpose: the original version proved the contract about a stub defined in the
  test file, which let four classes drift off it for three stages (ledger 135).
- `tests/test_plot_front_door.py` pins `myflopy.plot` itself — the verbs, the
  dispatch, and what is deliberately absent.
