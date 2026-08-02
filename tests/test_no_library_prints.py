"""A library talks through logging, not stdout (plan 7.2).

`print()` in library code cannot be silenced, redirected, levelled, or turned
*up* when you actually want the detail. It is also invisible in review: nothing
fails, the line just appears in somebody's notebook forever. `read_gpkg`
announced "Imported N features" on every single call.

**The exemption, stated once.** Output whose job is to report something a human
explicitly asked for keeps printing. In practice that means output behind a
flag the caller passed -- `verbose=True`, `progress=...`, `verbosity_level>0` --
plus one unconditional case (see ALLOWED_UNGATED_PRINTS). Everything else goes
to `myflopy._logging`, where an application decides whether it is shown.

The gate check is AST-based: a `print` inside an `if verbose:` (or a
`progress`/`verbosity` test) is exempt without needing to be listed, because
the flag IS the human asking. Everything else must be named.
"""

from __future__ import annotations

import ast
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "myflopy"

#: Names in an `if` test that mean "the caller asked for this output".
REQUEST_FLAGS = ("verbose", "progress", "verbosity")

#: Unconditional prints that stay, as `module:function`.
#:
#: Exactly one. `run_simulation` IS the interactive "did my model run" moment --
#: flopy is already printing MF6's own output beneath it (`silent=False`), and a
#: success line that only appears if you configured logging first would be a
#: worse answer to the question the caller just asked. The chatter AROUND it
#: (saving the .model snapshot, and failing to) is a side effect, and logs.
ALLOWED_UNGATED_PRINTS = {
    "modflow/mf6/simulation/runtime.py:run_simulation",
}


def _python_files() -> list[Path]:
    return [p for p in sorted(SRC.rglob("*.py")) if "_vendor" not in p.parts]


def _enclosing_name(tree: ast.Module, lineno: int) -> str:
    innermost, innermost_line = "", -1
    for parent in ast.walk(tree):
        if not isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        if parent.lineno <= lineno <= getattr(parent, "end_lineno", parent.lineno):
            if parent.lineno > innermost_line:
                innermost, innermost_line = parent.name, parent.lineno
    return innermost or "<module>"


def _print_calls():
    """Yield `(relative_path, lineno, key, requested)` for every `print(...)`."""

    for path in _python_files():
        tree = ast.parse(path.read_text(encoding="utf-8"))
        requested: set[int] = set()
        for node in ast.walk(tree):
            if not isinstance(node, ast.If):
                continue
            if not any(flag in ast.unparse(node.test) for flag in REQUEST_FLAGS):
                continue
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name):
                    if sub.func.id == "print":
                        requested.add(sub.lineno)
        relative = path.relative_to(SRC).as_posix()
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id == "print":
                    yield (
                        relative,
                        node.lineno,
                        f"{relative}:{_enclosing_name(tree, node.lineno)}",
                        node.lineno in requested,
                    )


def test_unconditional_prints_are_exactly_the_ones_allowed():
    """Fails in BOTH directions, so removing the last one is recorded too."""

    found = {key for _, _, key, requested in _print_calls() if not requested}
    added = sorted(found - ALLOWED_UNGATED_PRINTS)
    removed = sorted(ALLOWED_UNGATED_PRINTS - found)
    assert not added, (
        f"unconditional print() in library code: {added}. Use "
        "`myflopy._logging.get_logger(__name__)` -- logger.info for progress, "
        "logger.warning for a degraded result, logger.debug for detail. Put it "
        "behind `verbose=`/`progress=` only if the caller asks for it by name."
    )
    assert not removed, (
        f"{removed} no longer prints unconditionally. Delete it from "
        "ALLOWED_UNGATED_PRINTS and say so in the compromise ledger."
    )


def test_importing_myflopy_is_silent():
    """The plan's own acceptance test. An import that prints makes every
    downstream script's output unparseable and every notebook noisier."""

    result = subprocess.run(
        [sys.executable, "-c", "import myflopy"],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr[-2000:]
    assert result.stdout == "", f"`import myflopy` wrote to stdout: {result.stdout!r}"


def test_reading_a_geopackage_is_silent_but_loggable():
    """The concrete regression: `read_gpkg` printed "Imported N features" on
    every call, into every notebook. It still says so -- at DEBUG/INFO, where
    the reader can ask for it."""

    probe = textwrap.dedent(
        """
        import io, logging, contextlib, tempfile, pathlib
        import geopandas as gpd, shapely as shp
        from myflopy.modflow.utils.datatypes.readers import read_gpkg

        path = pathlib.Path(tempfile.mkdtemp()) / "p.gpkg"
        gpd.GeoDataFrame({"n": ["a"]}, geometry=[shp.Point(0, 0)],
                         crs="EPSG:2927").to_file(path, driver="GPKG")

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            read_gpkg(path)
        print("SILENT" if out.getvalue() == "" else f"PRINTED {out.getvalue()!r}",
              file=__import__("sys").stderr)

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        log = logging.getLogger("myflopy")
        log.addHandler(handler); log.setLevel(logging.DEBUG)
        read_gpkg(path)
        print("LOGGED" if any("1 features" in r.getMessage() for r in records)
              else f"NOTLOGGED {[r.getMessage() for r in records]}",
              file=__import__("sys").stderr)
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", probe], cwd=ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr[-3000:]
    assert "SILENT" in result.stderr, result.stderr[-3000:]
    assert "LOGGED" in result.stderr, result.stderr[-3000:]


@pytest.mark.parametrize(
    "relative", [p.relative_to(SRC).as_posix() for p in _python_files()]
)
def test_every_module_that_logs_defines_its_logger(relative):
    """`logger.info(...)` with no `logger` in scope is a NameError that only
    fires on the branch that logs -- which is exactly the branch nobody runs."""

    tree = ast.parse((SRC / relative).read_text(encoding="utf-8"))
    logs = any(
        isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "logger"
        for node in ast.walk(tree)
    )
    if not logs:
        return
    defines = any(
        isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "logger" for t in node.targets)
        for node in ast.walk(tree)
    )
    assert defines, f"{relative} calls logger.* but never assigns `logger`"
