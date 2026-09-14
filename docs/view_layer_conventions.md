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

**`backend="mpl"` is the static renderer, at every scope that has one** — the
free verbs, the bound namespaces, and the nouns. `plot.map`, `plot.section` and
`plot.grid` take it; `model.hds.map`, `packages.<pkg>.inputs.<field>.map` and
the rest take it. The NOUN-tier `mosaic`/`animate` take it too, which is where it
earns its keep, because a composer has no intermediate picture to call
`.plot_mpl()` on. It returns a bare `matplotlib.figure.Figure` everywhere.

Measured, so the sentence above stays honest: the FREE `mosaic` has no `backend`
parameter at all, and the free `animate`'s is `{"plotly", "png"}` — see below. Aliases `"matplotlib"`/`"static"`
and `"interactive"` resolve through one function, `normalize_backend`.

`surface` deliberately does NOT take it: there is no matplotlib rendering of a
3-D height field here, and a parameter that cannot act is the same defect as one
that is hidden. `animate`'s `backend={"plotly", "png"}` is a different axis — a
rasteriser over frames, not a second renderer of one subject.

**A section has two subjects, and `fill=` picks which** (ledger 177, 2026-09-02).
Without it, `section` answers *what is the head along this line* -- a profile,
two traces, no geometry. With it, *what does the model look like through here* --
the cells, their layers, and a field painted on them. `fill="layer"` colours by
layer and draws the simulated head as a water surface over the geology;
`fill="results"` colours the cells by the model's own dependent variable (head /
concentration / temperature, resolved through `_field_reader`) with a colorbar;
any `(nlay, ncpl)` array does the same for a field of your own. It stays ONE
verb because both are vertical slices -- geometry chooses the verb, and this is
the same geometry drawn two ways.

**Choosing what goes in it: `layers=`, `head_layers=`, `layer_labels=`.**
`layers=` picks which layers' CELLS are drawn — masked, and the vertical extent
cropped to what remains, because a full-height axis puts two layers in a thin
band. `head_layers=` picks whose water levels go on top (default `0`, the single
water table; a list draws one labelled surface each through `category_colors`, so
a layer keeps its colour across figures; `None` draws none). `layer_labels=`
names the units in the legend, defaulting to the names the model's own build
context carries — `ModelContext(surfaces=stack.build(vor))` keeps them through
`ModelSpec.build`, so a spec-declared model already knows its layers are "sand"
and "clay". **`layer=` RAISES in the filled branch**: it means "overlay these
layers' head profiles" on the other branch, and it silently did nothing here on
the first cut.

`fill=` is **matplotlib only** (FloPy's `PlotCrossSection` is the renderer and
has no Plotly counterpart), and it RAISES rather than falling back: with the
default backend, with any line-profile argument that cannot act on it
(`interpolate`, `spacing`, `num_points`, ...), with `fill="results"` on a bare
grid, and on a flat `(ncpl,)` array, which would paint every layer the same and
look like an answer. Before this existed the picture was
`myflopy.modflow.mf6.plot_model_cross_section` -- a deep import, on no scope,
and lost twice by the person who asked for it.

**`.plot_mpl()` on the returned picture is the same renderer reached a second
way, and is not the spelling to teach.** It still exists, still works, and is
what `backend="mpl"` calls underneath. But it only exists once a picture has been
built, which is why it cannot serve `mosaic`/`animate`, and reaching for it was
how `canonical_02` ended up dropping out of the grammar to get cell outlines.
Prefer the parameter. Until 2026-09-02 the verb tier had no `backend` at all
while the noun tier had had one since 8.2, so `model.plot.section(backend="mpl")`
was accepted, stored on the section object, and ignored — a Plotly figure back,
no error (ledger 176).

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

### Colour scales and custom geometry are options on `map`

Rule 1 in practice (ledger 252–254, 2026-09-13). Three options every map
takes, at every scope, verb and noun alike:

- **`logscale=True`** colours `log10(value)` — Plotly has no log colour axis —
  but the picture READS in real units: the colorbar is labelled at whole
  decades (`0.001 … 100`) on both backends, the hover keeps the real value,
  contours are drawn on the real field, and `zmin`/`zmax`/`zmid` are given in
  DATA units (`zmin=1e-4, zmax=10`; a non-positive one raises). The static
  backend draws through a `LogNorm`, so its colorbar is a true log axis. The
  one helper behind the labels is `viz.log_decade_colorbar`; `Choro`,
  `mosaic` and the PEST field maps all call it. Caller `colorbar=` keys ride on
  top, which is how the PEST change map keeps its `0.1x / 1x / 10x`.
- **`colorscale=`** takes any Plotly name (`'earth'`, `'Viridis'`, `'Turbo'`,
  `'RdBu_r'`, …) or explicit stops; the matplotlib colormap is built from the
  same stops. An unknown name warns and draws `earth` — the deliberate inverse
  of `_field_map_policy`, which raises, because a typo in a notebook cell
  should cost a warning and not the figure.
- **`show_layer_elevs=` / `hover_surfaces=`** put the grid's layer frame in the
  hover — the model top and each layer's bottom at that cell. They are **Tier 1**
  as of ledger 258: the elevations come from the GRID, not from the mapped
  field, so they act on a K map, a conductance map and a head map alike. They
  used to be Tier 2 and inert everywhere but heads, because
  `Choro._build_hover_context` fetched the elevations inside its heads branch and
  every noun passes `hover_heads=False`. Still **fetched only when asked**:
  `gdf_topbtm` is lazy and samples rasters. Where the map DOES have a per-layer
  field the elevations merge into that table as a second column; where it does
  not, they render alone. `hover_layers` stayed Tier 2 — it needs a per-layer
  field, which a record noun does not carry.
- **`features=`** draws custom geometry over the cells — wells, pond
  footprints, a boundary. A vector file path, a `GeoDataFrame`/`GeoSeries`, a
  shapely geometry, a `myflopy.plot.Feature` carrying its styling
  (`name`/`layer`/`color`/`width`/`size`/`fill`/`opacity`/`label`/`hover`/`crs`),
  or a list or dict of those. Points are markers, lines polylines, polygons
  outlines or fills; each layer is ONE legend entry (one `legendgroup`, however
  many geometry kinds it holds) with a stable colour from `category_colors` and
  a hover listing its attribute columns. A source with no CRS is taken to be in
  the grid's. Resolved eagerly, like `select=`, so a missing file raises at the
  call; rides `overlay_traces()`, so `mosaic` and `animate` keep it; drawn by
  `plot_mpl` too. The engine is `utils/datatypes/features.py` — one resolver,
  two renderers, so the backends cannot disagree about what a layer is.
  `locs=` is DEPRECATED: it warns, is named on no signature, and is drawn as a
  layer called "locations" through the same engine.

### The composers name the map options too

`mosaic` and `animate` forward every per-panel drawing option to each panel, so
(8.8, ledger 258) they NAME them: `packages.npf.k.mosaic(logscale=True,
features=wells, show_layer_elevs=True)` completes in an editor and carries a
description for each. They are defined ONCE, on `SpatialView`, so every noun
inherits one signature and one spliced docstring.

Two rules keep that honest. Only the options a caller actually **set** are
forwarded (`_explicit_options` compares against the composer's own signature),
so a panel otherwise keeps its noun's colour policy. And they are **refused, not
ignored**, on `kind="plot"` / `kind="section"`, which draw no choropleth — the
`**kwargs` tail still rides to the section branch, which has drawing options of
its own.

Why this rotted: `tests/test_noun_signatures.py` discovers its subjects by
looking for a body that mentions `plot.map(` or `_choropleth_factory`. A
composer's body mentions neither — it goes through `_facet_panels` — so while
every noun's `map` was converted in 8.8, `mosaic` named eleven parameters, none
of them drawing options. The file now sweeps the composers explicitly.

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

### The docstring is generated, and the fence says which half

The same argument reaches the docstring: a `__doc__` assembled at import reaches
`help()` and reaches no editor. So the spliced text is written into the source by
`scripts/derive_docstrings.py`, and the boundary is marked:

```
Own summary line -- yours, edit it here.

.. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
Bound form of :func:`myflopy.plot.section`.
...
.. end generated
```

Inside a fence is the free verb's reference, replaced wholesale on the next run.
Outside is yours, edited in place with no tool run.

**Adding a noun?** Three steps, and the middle one is the one that gets missed:

1. Write its own summary, prose, and only the parameters it documents
   DIFFERENTLY from the reference.
2. **Register it** — add the method to `_inherit_noun_docs()` in
   `myflopy/plot/__init__.py` (or, for a class ABOVE that module in the import
   graph, call `inherit_map_docs` at the bottom of its own module, as
   `project/group/lak.py` does). Nothing discovers a noun that is not wired in.
3. Run `python scripts/derive_docstrings.py` — the shared entries arrive fenced.

Skip step 2 and the script reports success while the docstring stays as you
wrote it. It is not silent forever — `test_the_sweep_finds_every_spliced_method_and_no_others`
pins the count — but it is silent at the moment you would want to hear about it.

Two things worth knowing before you write one. Prose belongs **above**
`Parameters`: a paragraph at column 0 inside the list is read as an entry head by
every NumPy parser, and `_refuse_malformed_entries` now raises rather than
propagating it. And a parameter you document that the reference documents too is
a **shadow** — your entry wins, permanently, so reference edits will never reach
it. That is often right (`layer` is narrower on a record noun), but it is how
three `reverse` entries came to say the opposite of what the verb did. Every
shadow is pinned in `tests/docstring_shadows.json`; read that diff.

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
- **The exchange scales ship as STOPS from their helpers, because the helpers
  own the ORIENTATION** (`_blue_white_red_diverging_colorscale`,
  `_red_white_blue_diverging_colorscale`): which end is gaining is a policy
  decision, and a bare name cannot carry it. Names themselves are safe since
  2026-09-13 (ledger 253): `colorscale=` accepts every Plotly name (`_r`
  reverses one) and `plot_mpl` builds its colormap from Plotly's own stops for
  that name, so `'RdBu'` renders the same way on both backends. Before that a
  hand table mapped `'rdbu'` to matplotlib's *reversed* colormap and every
  other name to `gist_earth` (ledger 69/70). Also note that `zmid` reaches only
  the Plotly trace; what centers a diverging map on *both* backends is
  symmetric `zmin`/`zmax`.
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
  callable when the labels depend on the pooled range. A grid of `logscale`
  panels no longer needs it — `mosaic` rebuilds the real-unit decades from the
  pooled limits itself (`viz.log_decade_colorbar`), and REFUSES to pool a log
  panel with a linear one, since one axis cannot mean both (ledger 252). For a
  *diverging* quantity also pass
  `diff=True`, or the neutral point lands wherever the pooled data happens to put
  it. Better still, refuse the combination that needs it: `plot_field_mosaic`
  accepts only stats that *have* a prior and a posterior form, which makes the
  trap unreachable rather than merely documented. Ledger 71.
- **Anything drawn OVER a static map takes its colormap from the same place the
  cells did** — `viz.mpl_colormap_for(colorscale)` (re-exported by `choros`,
  built from Plotly's own stops for names and explicit stops alike). Building a second
  colormap at the call site is how an overlaid point and the cell beneath it end
  up different colors for the same value. Note the Plotly and matplotlib overlay
  paths are genuinely different: `plot_mpl` draws `features=` layers and the
  `select=` outline but ignores raw `add_overlay` traces (they are Plotly
  objects), and it draws in **model coordinates**, while the Plotly path needs
  `vor.points_to_latlon`. Supporting both is ~10 lines; dropping the static one
  leaves a figure that silently shows the base map and none of its data.
  **A picture whose subject IS an overlay must therefore draw it itself under
  `backend="mpl"`**, not defer to `_apply_backend`. The two HFB verbs are the
  worked example: a barrier is an edge drawn over the cells, so a `plot_mpl`
  passthrough returns a perfectly valid figure of the grid with every barrier
  missing. Measured on the canonical wall — 32 barriers, 32 `Line2D` artists —
  which is the assertion, because a type check passes over the defect.
  Use `PALETTE.mpl_highlight`, not `PALETTE.highlight`: Plotly's
  `rgb(214,39,40)` spelling raises in matplotlib.
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
