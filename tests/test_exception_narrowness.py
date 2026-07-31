"""The 7.3 ratchet: no bare `except:`, and every broad catch is argued for.

Plan 7.3 narrowed ~60 handlers. Without a test, the next broad handler added in
a hurry costs nothing to write and is invisible in review -- which is how the
first sixty arrived. So this file makes the two rules mechanical:

1. **No bare `except:` anywhere.** It catches `KeyboardInterrupt` and
   `SystemExit`, so a bare handler inside a loop makes a long build
   un-interruptible. There were eleven; there are none.
2. **`except Exception` must carry `# noqa: BLE001`.** Not because the marker
   proves anything on its own, but because a reviewer typing it has to decide
   the catch is genuinely unenumerable and (by convention) write the sentence
   saying why. Ten survive that bar today, each for a stated reason: flopy's
   `voronoi` module contains a literal `raise Exception(...)`, `__dir__` and
   `__repr__` must never raise, `pickle.dump` walks an unbounded object graph,
   and a `setattr` on an object a user's own builder returned can raise
   anything at all.

The allowlist is EXACT in both directions, like `tests/api_snapshot.json`.
Adding a broad catch fails; removing one also fails, so the win gets recorded
rather than absorbed.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "myflopy"

#: Broad `except Exception` handlers that survive, as `module:function`.
#:
#: `_vendor` is excluded from the scan entirely -- it is third-party code we do
#: not edit.
ALLOWED_BROAD = {
    # flopy's `utils/voronoi.py` contains a literal `raise Exception(...)` for a
    # degenerate triangulation, so `Exception` is the only clause that catches
    # the case these guards exist for. None of the three swallows: each records
    # the failure in a report the caller reads.
    "modflow/mf6/grid/mesh_quality.py:triangle_quality_report",
    "modflow/mf6/grid/triangle.py:optimize_seeds",
    # An autocomplete that raises, and a repr that breaks the traceback you were
    # reading, are worse than a short answer. The budget-file chain behind them
    # was measured to raise NotImplementedError (from flopy's base `Grid.shape`)
    # and MFDataException, which no builtin tuple covers.
    "modflow/mf6/package_results.py:__dir__",
    "modflow/mf6/package_results.py:__repr__",
    # `pickle.dump(model)` serializes an unbounded third-party object graph, and
    # an escaped exception here would kill the MF6 run on the next line.
    "modflow/mf6/simulation/runtime.py:run_simulation",
    # Persistence of the splitter node mapping is best-effort; the run's results
    # are reconstructed from the in-memory splitter regardless.
    "modflow/mf6/parallel.py:write",
    # `grid` is whatever a USER'S builder script returned, so the types a
    # `setattr` on it can raise are genuinely open (pydantic raises a plain
    # ValueError). All that is lost is a provenance tag.
    "grid_spec_resolver.py:_resolve_python",
    # `model.vor` is a lazy property on file-backed models that loads flopy and
    # rebuilds the grid, reaching flopy, geopandas, pyproj, shapely and the
    # filesystem -- and its crs comes from user input.
    "modflow/mf6/grid/interpolated_surface.py:__init__",
    # Does not swallow: re-raises every failure to load the optional GRASS
    # SYSTEM dependency as one actionable ImportError, cause attached.
    "modflow/utils/contour_interp.py:_grass_modules",
}


def _python_files():
    return [p for p in sorted(SRC.rglob("*.py")) if "_vendor" not in p.parts]


def _enclosing_name(tree: ast.Module, node: ast.ExceptHandler) -> str:
    """The innermost enclosing def/class name, for a readable allowlist key.

    Keyed by NAME rather than line number so the allowlist survives ordinary
    edits above the handler. The trade: two broad handlers in one function
    collapse to a single key, so adding a second to an already-allowed function
    would not fail. That is deliberate -- keying by line would turn every
    unrelated edit into a test failure, which is how ratchets get deleted.
    """

    innermost, innermost_line = "", -1
    for parent in ast.walk(tree):
        if not isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        if parent.lineno <= node.lineno <= getattr(parent, "end_lineno", parent.lineno):
            if parent.lineno > innermost_line:
                innermost, innermost_line = parent.name, parent.lineno
    return innermost or "<module>"


def _handlers():
    """Yield ``(relative_path, lineno, kind, key, source_line)`` per handler."""

    for path in _python_files():
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        tree = ast.parse(text)
        for node in ast.walk(tree):
            if not isinstance(node, ast.ExceptHandler):
                continue
            relative = path.relative_to(SRC).as_posix()
            source = lines[node.lineno - 1]
            if node.type is None:
                kind = "bare"
            elif isinstance(node.type, ast.Name) and node.type.id == "Exception":
                kind = "broad"
            else:
                continue
            yield relative, node.lineno, kind, f"{relative}:{_enclosing_name(tree, node)}", source


def test_no_bare_except_anywhere():
    """`except:` also catches KeyboardInterrupt and SystemExit. One of the
    eleven this pass removed sat in the innermost cell loop of a GHB build, so
    Ctrl-C during a long boundary build did nothing at all."""

    bare = [f"{path}:{line}" for path, line, kind, _, _ in _handlers() if kind == "bare"]
    assert not bare, (
        f"bare `except:` at {bare}. Name the exceptions, or use "
        "`except Exception:` with a `# noqa: BLE001` and a reason if the set "
        "genuinely cannot be closed."
    )


def test_every_broad_except_is_marked():
    """The marker is the tripwire that forces the decision to be conscious."""

    unmarked = [
        f"{path}:{line}"
        for path, line, kind, _, source in _handlers()
        if kind == "broad" and "BLE001" not in source
    ]
    assert not unmarked, (
        f"`except Exception` without a `# noqa: BLE001` marker at {unmarked}. "
        "Narrow it to the exceptions the guarded call can actually raise -- or "
        "if the set is genuinely open, add the marker AND a comment saying why "
        "(see docs/compromises_and_deferrals.md, plan 7.3)."
    )


def test_the_broad_handlers_are_exactly_the_ones_argued_for():
    """Fails in BOTH directions. A new broad catch fails; closing one also
    fails, so the win is written down rather than silently absorbed."""

    found = {key for _, _, kind, key, _ in _handlers() if kind == "broad"}
    added = sorted(found - ALLOWED_BROAD)
    closed = sorted(ALLOWED_BROAD - found)
    assert not added, (
        f"new broad exception handlers: {added}. If the exception set really "
        "cannot be enumerated, add it here with the sentence explaining why."
    )
    assert not closed, (
        f"{closed} no longer catches broadly -- good. Delete it from "
        "ALLOWED_BROAD and record the narrowing in the compromise ledger."
    )


@pytest.mark.parametrize("path", [p.relative_to(SRC).as_posix() for p in _python_files()])
def test_every_file_parses(path):
    """A cheap syntax guard: this suite edited 27 files by script."""

    ast.parse((SRC / path).read_text(encoding="utf-8"))
