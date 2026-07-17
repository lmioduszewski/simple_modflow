# Deprecation policy

How myflopy retires a public name without breaking existing scripts
(implementation plan 3.1; decisions D11/D12).

## The rules

1. **Every deprecated name keeps working** and emits a `DeprecationWarning`
   that names its replacement and the release that deprecated it (`since`).
2. **Deprecated names are runtime-only (D12).** They resolve *exclusively*
   through `__getattr__` — module-level (PEP 562) or class-level — and stay
   out of `__all__`, `dir()`, `TYPE_CHECKING` import blocks, and any stubs.
   Old code runs with a warning; IDEs and tab-completion never suggest the
   old names to anyone.
3. **A deprecated name lives at least two tagged releases** after its
   `since` version before it may be removed. The clock is real: `v0.1.0`
   was tagged at the Phase 0 baseline and releases bump at plan milestones
   (D11). Example: a name with `since = "0.2.0"` is removable in `0.4.0` at
   the earliest.
4. **`myflopy.__compatibility__` is authoritative.** It is the static tuple
   (in `src/myflopy/__init__.py`) of every warned alias, as fully-qualified
   old names. `tests/test_deprecation.py` cross-checks it against the
   registry that the helpers populate — adding a warned alias without
   listing it there (or vice versa) fails the suite.

## The mechanism — `myflopy/_deprecation.py`

One module provides the whole mechanism; do not hand-roll warnings or
ad-hoc `__getattr__` redirects elsewhere.

- `warn_deprecated(old, new, *, since)` — the standard warning. `old` must
  be fully qualified (`"myflopy...."`); messages start with that name (see
  "Enforcement" below). Also used directly for *behavior* deprecations
  (e.g. the MP3DU package-directory executable fallback) — those are not
  aliases and are not listed in `__compatibility__`.
- `deprecated_module_getattr(mapping, __name__, *, fallback_getattr=None,
  fallback_dir=None)` — returns a `(__getattr__, __dir__)` pair to assign
  at module bottom. `mapping` is `{old: (target, replacement, since)}` with
  `target` `"pkg.mod:attr"` or `"pkg.mod"`. Facades with existing
  lazy-export machinery chain it via the fallbacks; the returned `__dir__`
  excludes the deprecated names either way.
- `deprecated_instance_getattr(mapping, owner)` — returns a class-level
  `__getattr__` for deprecated methods/attributes; `mapping` is
  `{old: (instance_attr, replacement, since)}`. The old names must not also
  exist as real class attributes (that would defeat the hiding).

## Adding a deprecation (checklist)

1. Rename/move the real thing; wire the old name through one of the two
   `__getattr__` helpers (never a plain assignment or a bare
   `warnings.warn`).
2. Remove the old name from `__all__`, `TYPE_CHECKING` blocks, and stubs;
   add the new name there instead.
3. Add the fully-qualified old name to `myflopy.__compatibility__`.
4. Add tests: `pytest.warns(DeprecationWarning)` resolves the old name to
   the same object as the new one; `old not in dir(module)`;
   `old not in __all__`.
5. Update every doc/guide/notebook that teaches the old name (canonical
   rule: docs move in the same commit as the code).

## Enforcement in the test suite

`pyproject.toml` sets `filterwarnings = ["error:myflopy:DeprecationWarning"]`:
any `DeprecationWarning` whose message starts with `myflopy` (ours all do —
the message begins with the fully-qualified old name) becomes an **error**
unless a test captures it with `pytest.warns`. Internal code can therefore
never call a deprecated alias, and tests must exercise aliases explicitly.

## Current deprecations

| Old name | Use instead | Since | Removable in |
|---|---|---|---|
| `myflopy.modflow.mf6.GHB` (facade) and `myflopy.modflow.mf6.ghb.GHB` | `GHBFromVector` (same modules) | 0.2.0 | 0.4.0 |
| `myflopy.modflow.mf6.DRN` (facade) and `myflopy.modflow.mf6.drn.DRN` | `DRNFromVector` (same modules) | 0.2.0 | 0.4.0 |
| `ModelGroup.rch/.chd/.drn/.ghb/.wel/.uzf` flat shortcuts | `group.packages.<pkg>.inputs` | 0.1.0 | 0.3.0 |
| `myflopy.modflow.mp3du.particles.PRT/PrtMip/PrtOc/PrtPrp/PrtDisv/PrtFmi` | `myflopy.modflow.mp3du.legacy_prt.<name>` | 0.1.0 | 0.3.0 |

Behavior deprecations (warned, not aliases): loading MP3DU executables from
inside the package (`since` 0.1.0) — use `tools/mp3du/` or
`MYFLOPY_MP3DU_DIR`.

**Not deprecated, despite the similar name:** the `GHB` class in
`modflow/mf6/simulation/packages.py` (re-exported by
`myflopy.modflow.mf6.simulation` and `mfsimbase`) is the OO package
wrapper, a different class from the retired bare-`GHB` spelling of the
vector builder. It keeps its name.

## Related tiers (not deprecation)

`myflopy.__engine__` (formerly the meaning of `__compatibility__`, renamed
when 3.1 repurposed that dunder) lists the *unwarned* second-tier exports —
OO builders and `*_spec` factories that remain fully supported engine API
but are kept out of `__all__`/`dir()` so discovery points at the
package-first facade. See `docs/manual/04_architecture.md` §4.8.
