# View-layer conventions (`model.packages.…`)

> **Status:** normative. Added 2026-07-18 after `sfr.results.long_profile` /
> `plot_long_profile` drifted outside the grammar and two notebooks
> hand-rolled ~20 lines of matplotlib to redraw a figure the library already
> built — losing `scrollZoom`, pan, the house template, *and* the documented
> gaining/losing colors in the process.

## Why this file exists

The **spatial** half of the grammar was already decided and documented
(`docs/myflopy_context.md`: *"`map/plot/xs` + `mosaic/animate` on every leaf"*,
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
| `plot(...)` | `viz.Fig` | every noun |
| `map(...)` | choropleth | spatial nouns |
| `xs(...)` | cross-section | spatial nouns |
| `mosaic(...)` / `animate(...)` | multi-panel / animated | spatial nouns |

Rules that hold for all of them:

- **Figures are always `viz.Fig`**, never a bare `go.Figure`. `Fig` carries
  `dragmode='pan'`, `scrollZoom`, and the house template. Returning a bare
  figure silently drops all three — it is the defect this convention exists to
  prevent, and it has now happened twice (`hds.animate`, and these notebooks).
- **Colors follow the policy**, not the call site. Signed q-like fields use
  `_blue_white_red_diverging_colorscale()` — gaining/negative BLUE, losing/
  positive RED (`docs/mf6io_reference.md`). Read the endpoints off that helper
  rather than writing hex literals, so discrete and continuous renderings
  cannot drift apart.
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

Namespace-level nouns merge fields; field-level nouns cover a single field.

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
