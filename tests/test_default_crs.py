"""The default CRS is defined once (ledger 225).

EPSG:2927 was spelled three ways at 19 code sites, so changing or even finding the
default meant a grep, and a new site could pick any spelling. A raw 2927 anywhere
in ``src/myflopy`` outside ``_defaults.py`` fails here, naming the site.
"""

from __future__ import annotations

import ast
from pathlib import Path

import myflopy
from myflopy._defaults import DEFAULT_CRS, DEFAULT_EPSG

SRC = Path(myflopy.__file__).parent


def _literals(path: Path):
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
                docstrings.add(id(first.value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and id(node) not in docstrings:
            value = node.value
            if (value == 2927 and not isinstance(value, bool)) or (
                isinstance(value, str) and value.strip().upper() in {"2927", "EPSG:2927"}
            ):
                yield node.lineno


def test_no_module_spells_the_default_crs_itself():
    stray = [
        f"{path.relative_to(SRC)}:{line}"
        for path in sorted(SRC.rglob("*.py"))
        if path.name != "_defaults.py" and "_vendor" not in path.parts
        for line in _literals(path)
    ]
    assert not stray, f"use myflopy._defaults.DEFAULT_CRS / DEFAULT_EPSG at: {stray}"


def test_the_default_is_still_2927():
    """Location-agnostic models sit at its origin on purpose; this must not drift."""

    assert (DEFAULT_EPSG, DEFAULT_CRS) == (2927, "EPSG:2927")
