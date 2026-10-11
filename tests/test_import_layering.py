"""Import layering + deferred-import ratchet (implementation plan 4.4).

Three pinned invariants over ``src/myflopy`` (``_vendor`` excluded):

1. **Acyclic** — the module-level runtime import graph has no cycles.
2. **Layered** — every module-level runtime import points at a module whose
   layer (``tests/import_layers.json``) is the same or lower. Upward
   references belong in TYPE_CHECKING blocks or deferred imports.
3. **Current** — the three pinned files equal what the derivation writes now,
   byte for byte. The first two checks pass on a STALE map whenever a module
   moved to a higher layer without crossing a lower neighbour's -- which is how
   ``grid.tessellation`` sat at layer 0 while the graph put it at 1.
4. **Ratchet** — function-level ``myflopy`` imports per module match
   ``tests/deferred_import_allowlist.json`` EXACTLY. Hoisted one? Lower the
   allowlist. Added one? That is a regression — import at module level (the
   graph is a DAG; there is almost always a downward path) or justify a new
   allowlist entry in review.

After legitimate structural changes regenerate all pinned files with
``python scripts/derive_import_layers.py`` and review the diff.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from derive_import_layers import build_graph, sccs_and_depths, stale_files  # noqa: E402

LAYERS = json.loads((ROOT / "tests/import_layers.json").read_text())
ALLOWLIST = json.loads((ROOT / "tests/deferred_import_allowlist.json").read_text())

MODULES, EDGES, DEFERRED = build_graph()


def test_every_module_has_a_layer():
    missing = sorted(set(MODULES) - set(LAYERS))
    stale = sorted(set(LAYERS) - set(MODULES))
    assert not missing and not stale, (
        f"layer map out of date (missing={missing}, stale={stale}) — run "
        "`python scripts/derive_import_layers.py` and review the diff"
    )


def test_module_level_import_graph_is_acyclic():
    sccs, _, _ = sccs_and_depths(EDGES)
    cycles = [sorted(c) for c in sccs if len(c) > 1]
    assert not cycles, f"module-level runtime import cycles: {cycles}"


def test_runtime_imports_never_point_upward():
    violations = [
        f"{mod} (L{LAYERS[mod]}) -> {target} (L{LAYERS[target]})"
        for mod, targets in EDGES.items()
        for target in targets
        if LAYERS[target] > LAYERS[mod]
    ]
    assert not violations, (
        "module-level runtime imports may only point at the same or a lower "
        f"layer; move these to TYPE_CHECKING or restructure: {violations}"
    )


def test_deferred_import_ratchet():
    regressions = {
        m: (DEFERRED.get(m, 0), ALLOWLIST.get(m, 0))
        for m in set(DEFERRED) | set(ALLOWLIST)
        if DEFERRED.get(m, 0) != ALLOWLIST.get(m, 0)
    }
    assert not regressions, (
        "deferred (function-level) myflopy imports changed — module: (now, "
        f"allowed): {regressions}. If you hoisted imports, ratchet DOWN by "
        "regenerating with `python scripts/derive_import_layers.py`; if a "
        "count went UP, hoist the import to module level instead."
    )


def test_the_pinned_layer_files_are_what_the_derivation_writes():
    """A layer map that has fallen behind the code fails here, naming the file.

    ``test_every_module_has_a_layer`` compares the module SET and
    ``test_runtime_imports_never_point_upward`` only that no edge climbs, so a
    module whose layer ROSE (``grid.tessellation`` 0 -> 1, ``sfr`` 3 -> 4) left
    ``tests/import_layers.json`` and ``docs/import_layering.md`` wrong while all
    of them passed.
    """

    stale = stale_files()
    assert not stale, (
        f"{stale} differ from what `python scripts/derive_import_layers.py` writes; "
        "rerun it and review the diff"
    )
