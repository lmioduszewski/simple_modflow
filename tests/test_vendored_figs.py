"""The vendored figs snapshot must fully substitute the external figs.

Three contracts (implementation plan §1.1):
1. The vendored subtree imports standalone and provides the full surface
   myflopy uses (Fig/Subplot/Template/create_hover + mpl REPORT/Theme/
   get_mplfig), including the map-sync hooks (add_post_script, plotly-
   compatible write_html).
2. With the external figs hidden, ``import myflopy.viz`` works, resolves to
   the vendored classes, and post-scripts still inject into write_html
   output (the map-sync mechanism).
3. Single-importer rule: outside ``_vendor``, only ``myflopy/viz.py`` may
   import figs at runtime (TYPE_CHECKING imports exempt).
"""

from __future__ import annotations

import ast
import inspect
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def test_vendored_figs_surface():
    from myflopy._vendor.figs import Fig, Subplot, Template, create_hover
    from myflopy._vendor.figs.mpl import REPORT, Theme, get_mplfig

    fig = Fig()
    assert hasattr(fig, "add_post_script")
    assert "file" in inspect.signature(fig.write_html).parameters
    assert callable(create_hover) and callable(get_mplfig)
    assert Subplot is not None and Template is not None
    assert isinstance(REPORT, Theme)


def test_vendored_figs_has_no_heavy_extra_deps():
    """The trimmed closure must not import figs' export/scaling stack.

    Runs in a subprocess (with the external figs blocked) because in-process
    the external figs — whose closure legitimately pulls these libraries —
    may already be in ``sys.modules`` from other tests.
    """

    code = """
import importlib.abc
import sys


class _BlockFigs(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] == "figs":
            raise ImportError("figs is blocked for this test")
        return None


sys.meta_path.insert(0, _BlockFigs())

import myflopy._vendor.figs
import myflopy._vendor.figs.mpl

for heavy in ("bokeh", "cairosvg", "reportlab", "svglib", "svgpathtools"):
    assert heavy not in sys.modules, f"vendored figs pulled in {heavy}"
print("ok")
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(SRC), "MPLBACKEND": "Agg"},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0 and "ok" in result.stdout, (
        result.stdout[-1000:] + "\n" + result.stderr[-2000:]
    )


def test_viz_works_without_external_figs(tmp_path):
    """Hide figs in a subprocess; myflopy.viz must fall back to the vendored copy."""

    out_html = tmp_path / "post_script_check.html"
    code = f"""
import importlib.abc
import sys


class _BlockFigs(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] == "figs":
            raise ImportError("figs is blocked for this test")
        return None


sys.meta_path.insert(0, _BlockFigs())

from myflopy import viz

assert viz.Fig.__module__.startswith("myflopy._vendor.figs"), viz.Fig.__module__

fig = viz.Fig()
fig.add_scattergl(x=[0, 1], y=[0, 1])
fig.add_post_script("/*MAPSYNC-VENDOR-TEST*/")
fig.write_html({str(out_html)!r})
html = open({str(out_html)!r}, encoding="utf-8").read()
assert "/*MAPSYNC-VENDOR-TEST*/" in html, "post-script did not inject"
print("ok")
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(SRC), "MPLBACKEND": "Agg"},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0 and "ok" in result.stdout, (
        result.stdout[-1000:] + "\n" + result.stderr[-2000:]
    )


def _scan_paths():
    # tests/ must not import figs either: external figs does not exist on
    # CI/installed machines — everything goes through myflopy.viz.
    yield from (SRC / "myflopy").rglob("*.py")
    yield from (ROOT / "tests").rglob("*.py")


def _runtime_figs_importers() -> list[str]:
    offenders = []
    for path in _scan_paths():
        if "_vendor" in path.parts:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))

        type_checking_spans: list[tuple[int, int]] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.If):
                test = node.test
                name = (
                    test.id
                    if isinstance(test, ast.Name)
                    else test.attr
                    if isinstance(test, ast.Attribute)
                    else None
                )
                if name == "TYPE_CHECKING":
                    type_checking_spans.append((node.lineno, node.end_lineno or node.lineno))

        def _is_type_checking(lineno: int) -> bool:
            return any(start <= lineno <= end for start, end in type_checking_spans)

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                roots = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                roots = [node.module.split(".")[0]]
            else:
                continue
            if "figs" in roots and not _is_type_checking(node.lineno):
                offenders.append(f"{path.relative_to(ROOT)}:{node.lineno}")
    return offenders


def test_viz_is_the_only_runtime_figs_importer():
    offenders = _runtime_figs_importers()
    outside_viz = [o for o in offenders if not o.startswith("src/myflopy/viz.py:")]
    assert not outside_viz, (
        "figs may only be imported by myflopy/viz.py at runtime "
        f"(TYPE_CHECKING exempt); offenders: {outside_viz}"
    )
    assert offenders, "viz.py should import figs (external-first) — AST scan is broken"
