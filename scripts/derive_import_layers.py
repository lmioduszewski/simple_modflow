"""Regenerate the import-layering pin from the real AST graph (plan 4.4).

Writes:
- ``tests/import_layers.json``            {module: layer} — layer = longest-path
  depth of the module's SCC in the module-level runtime import graph
- ``tests/deferred_import_allowlist.json`` {module: count} — function-level
  ``myflopy`` imports (the deferred-import ratchet; counts may only go DOWN)
- ``docs/import_layering.md``             the human-readable layer map (an
  OUTPUT of this derivation, per the plan — never hand-edit it)

Run it after adding a module or hoisting a deferred import, review the diff,
and commit the three files together. ``tests/test_import_layering.py`` pins
the results: no module-level runtime import may point at a higher layer, the
graph must stay acyclic, and deferred-import counts must match the allowlist
exactly (lower the allowlist when you hoist; never raise it).

Rules of the graph:
- module-level runtime imports only — ``TYPE_CHECKING`` blocks and
  function-level (deferred) imports are excluded
- ``myflopy._vendor`` is external (its parent-package side effects are
  import machinery, not architecture)
- edges to ancestor packages (their ``__init__`` runs first anyway) are
  import machinery, not architecture — the root lazy-export table stays L0
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
PKG = SRC / "myflopy"


def module_name(path: Path) -> str:
    rel = path.relative_to(SRC).with_suffix("")
    parts = list(rel.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def collect_modules() -> dict[str, Path]:
    return {
        module_name(p): p
        for p in sorted(PKG.rglob("*.py"))
        if "_vendor" not in p.parts
    }


def resolve_target(raw: str, module_names: set[str]) -> str | None:
    if raw.startswith("myflopy._vendor"):
        return None
    parts = raw.split(".")
    for i in range(len(parts), 0, -1):
        cand = ".".join(parts[:i])
        if cand in module_names:
            return cand
    return None


def top_level_runtime_imports(tree: ast.Module) -> list[str]:
    """myflopy import targets at module level, outside TYPE_CHECKING blocks."""

    out: list[str] = []

    def visit_body(body):
        for node in body:
            if isinstance(node, ast.If):
                test = node.test
                name = (
                    test.id if isinstance(test, ast.Name)
                    else test.attr if isinstance(test, ast.Attribute)
                    else None
                )
                if name == "TYPE_CHECKING":
                    continue
                visit_body(node.body)
                visit_body(node.orelse)
            elif isinstance(node, ast.Try):
                visit_body(node.body)
                for handler in node.handlers:
                    visit_body(handler.body)
                visit_body(node.orelse)
                visit_body(node.finalbody)
            elif isinstance(node, ast.Import):
                out.extend(a.name for a in node.names if a.name.startswith("myflopy"))
            elif isinstance(node, ast.ImportFrom):
                if node.level == 0 and node.module and node.module.startswith("myflopy"):
                    out.extend(f"{node.module}.{a.name}" for a in node.names)

    visit_body(tree.body)
    return out


def deferred_import_count(tree: ast.Module) -> int:
    """Function-level myflopy import statements (the ratchet metric)."""

    count = 0
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            for sub in ast.walk(node):
                if isinstance(sub, ast.ImportFrom):
                    if sub.level == 0 and sub.module and sub.module.startswith("myflopy"):
                        count += 1
                elif isinstance(sub, ast.Import):
                    count += sum(1 for a in sub.names if a.name.startswith("myflopy"))
    return count


def build_graph():
    modules = collect_modules()
    names = set(modules)
    edges: dict[str, set[str]] = {m: set() for m in modules}
    deferred: dict[str, int] = {}
    for mod, path in modules.items():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for raw in top_level_runtime_imports(tree):
            target = resolve_target(raw, names)
            if target is None or target == mod:
                continue
            if mod.startswith(target + "."):
                continue  # ancestor package: import machinery, not architecture
            edges[mod].add(target)
        n = deferred_import_count(tree)
        if n:
            deferred[mod] = n
    return modules, edges, deferred


def sccs_and_depths(edges: dict[str, set[str]]):
    sys.setrecursionlimit(100_000)
    counter = [0]
    index: dict[str, int] = {}
    lowlink: dict[str, int] = {}
    on_stack: dict[str, bool] = {}
    stack: list[str] = []
    sccs: list[list[str]] = []

    def strongconnect(v: str):
        index[v] = lowlink[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        on_stack[v] = True
        for w in edges[v]:
            if w not in index:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif on_stack.get(w):
                lowlink[v] = min(lowlink[v], index[w])
        if lowlink[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            sccs.append(comp)

    for v in edges:
        if v not in index:
            strongconnect(v)

    comp_of = {m: i for i, comp in enumerate(sccs) for m in comp}
    depth: dict[int, int] = {}
    for i, comp in enumerate(sccs):  # Tarjan emits reverse topological order
        d = 0
        for m in comp:
            for t in edges[m]:
                if comp_of[t] != i:
                    d = max(d, depth[comp_of[t]] + 1)
        depth[i] = d
    return sccs, comp_of, depth


def main() -> None:
    modules, edges, deferred = build_graph()
    sccs, comp_of, depth = sccs_and_depths(edges)
    cycles = [sorted(c) for c in sccs if len(c) > 1]
    if cycles:
        for c in cycles:
            print("CYCLE:", " <-> ".join(c))
        raise SystemExit("module-level import cycles found — fix before pinning")

    layers = {m: depth[comp_of[m]] for m in sorted(modules)}
    (ROOT / "tests/import_layers.json").write_text(
        json.dumps(layers, indent=1, sort_keys=True) + "\n"
    )
    (ROOT / "tests/deferred_import_allowlist.json").write_text(
        json.dumps(deferred, indent=1, sort_keys=True) + "\n"
    )

    by_layer: dict[int, list[str]] = {}
    for m, layer in layers.items():
        by_layer.setdefault(layer, []).append(m)
    total = sum(deferred.values())
    doc = [
        "# Import layering — derived map (do not hand-edit)",
        "",
        "Generated by `scripts/derive_import_layers.py` (implementation plan 4.4)",
        "and pinned by `tests/test_import_layering.py`. A module's layer is its",
        "longest-path depth in the module-level runtime import graph",
        "(TYPE_CHECKING and function-level imports excluded; `_vendor` external).",
        "",
        "**Rules** — module-level runtime imports may only point at the same or a",
        "lower layer, and the graph must stay acyclic. Upward references belong in",
        "`TYPE_CHECKING` blocks or (sparingly) function-level deferred imports —",
        "the deferred-import ratchet (`tests/deferred_import_allowlist.json`,",
        f"currently **{total}** function-level `myflopy` imports) counts those and",
        "only ever goes down. To update after a legitimate change: run the script,",
        "review the diff, commit the regenerated files together.",
        "",
    ]
    for layer in sorted(by_layer):
        doc.append(f"## Layer {layer}")
        doc.append("")
        for m in by_layer[layer]:
            suffix = f" *(deferred imports: {deferred[m]})*" if m in deferred else ""
            doc.append(f"- `{m}`{suffix}")
        doc.append("")
    (ROOT / "docs/import_layering.md").write_text("\n".join(doc))
    print(
        f"modules={len(modules)} layers=0..{max(by_layer)} "
        f"deferred_total={total} — wrote tests/import_layers.json, "
        "tests/deferred_import_allowlist.json, docs/import_layering.md"
    )


if __name__ == "__main__":
    main()
