"""Every reproducer and example script removes the scratch directory it makes (ledger 460).

``/tmp`` is a 16 GB tmpfs with a per-user quota, and when it is full every write the
user makes there fails (ledger 335). Fifteen of the scripts under
``docs/reproducers/`` and ``examples/`` called ``tempfile.mkdtemp()`` and never
removed what it made -- 12 MB a run for ``flow_targets_period_end.py``.

The rule this checks, statically, because running every script takes ~40 minutes of
MF6 and PEST++ (ledger 444): a ``mkdtemp()`` result is bound to a name, and that
name is removed by ``shutil.rmtree`` inside a ``finally:`` block or registered with
``atexit.register(shutil.rmtree, <name>, ...)``; ``mkstemp()`` and
``NamedTemporaryFile(delete=False)`` are not used at all, and a
``TemporaryDirectory`` cleans itself. ``tests/test_examples_run.py`` runs the
scripts it executes with ``TMPDIR`` in a directory of its own and fails naming
whatever a run leaves there, which catches what the library leaves as well.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted(
    [*(ROOT / "docs" / "reproducers").glob("*.py"), *(ROOT / "examples").rglob("*.py")]
)


def _call_name(node: ast.AST) -> str | None:
    """``"tempfile.mkdtemp"`` for ``tempfile.mkdtemp(...)``, ``"mkdtemp"`` for a bare call."""

    if not isinstance(node, ast.Call):
        return None
    func = node.func
    if isinstance(func, ast.Attribute):
        base = func.value.id if isinstance(func.value, ast.Name) else "?"
        return f"{base}.{func.attr}"
    return func.id if isinstance(func, ast.Name) else None


def _is_mkdtemp(node: ast.AST) -> bool:
    return (_call_name(node) or "").split(".")[-1] == "mkdtemp"


def _first_name(call: ast.Call) -> str | None:
    """The name a cleanup call removes: its first positional argument, or the one after
    ``shutil.rmtree`` in ``atexit.register(shutil.rmtree, NAME, ...)``."""

    args = call.args
    if (_call_name(call) or "").endswith("register") and len(args) >= 2:
        args = args[1:]
    return args[0].id if args and isinstance(args[0], ast.Name) else None


def scratch_findings(source: str) -> list[str]:
    """What in ``source`` makes scratch it does not remove, one message per site."""

    tree = ast.parse(source)
    findings: list[str] = []
    removed: set[str] = set()
    for node in ast.walk(tree):
        # `finally:` bodies: shutil.rmtree(NAME, ...)
        if isinstance(node, ast.Try):
            for inner in (n for stmt in node.finalbody for n in ast.walk(stmt)):
                if (_call_name(inner) or "").endswith("rmtree"):
                    name = _first_name(inner)
                    if name:
                        removed.add(name)
        # atexit.register(shutil.rmtree, NAME, ...)
        if (_call_name(node) or "") == "atexit.register" and node.args and (
            isinstance(node.args[0], ast.Attribute) and node.args[0].attr == "rmtree"
        ):
            name = _first_name(node)
            if name:
                removed.add(name)
        name = (_call_name(node) or "").split(".")[-1]
        if name == "mkstemp":
            findings.append(f"line {node.lineno}: mkstemp() -- use a TemporaryDirectory")
        if name == "NamedTemporaryFile" and any(
            k.arg == "delete" and isinstance(k.value, ast.Constant) and k.value.value is False
            for k in node.keywords
        ):
            findings.append(f"line {node.lineno}: NamedTemporaryFile(delete=False)")
    bound: dict[int, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(
            node.targets[0], ast.Name
        ):
            for inner in ast.walk(node.value):
                if _is_mkdtemp(inner):
                    bound[id(inner)] = node.targets[0].id
    for node in ast.walk(tree):
        if _is_mkdtemp(node):
            target = bound.get(id(node))
            if target is None:
                findings.append(f"line {node.lineno}: mkdtemp() not bound to a name")
            elif target not in removed:
                findings.append(
                    f"line {node.lineno}: {target} = mkdtemp() is never removed "
                    "(shutil.rmtree in a finally:, or atexit.register(shutil.rmtree, ...))"
                )
    return findings


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: str(p.relative_to(ROOT)))
def test_the_script_removes_the_scratch_it_makes(script):
    assert scratch_findings(script.read_text(encoding="utf-8")) == []


@pytest.mark.parametrize(("source", "expected"), [
    ("import tempfile\nroot = tempfile.mkdtemp()\n", ["never removed"]),
    ("import tempfile\nf(tempfile.mkdtemp())\n", ["not bound to a name"]),
    ("import tempfile, shutil\nroot = tempfile.mkdtemp()\ntry:\n    pass\n"
     "finally:\n    shutil.rmtree(root)\n", []),
    ("import tempfile, shutil, atexit\nroot = Path(tempfile.mkdtemp())\n"
     "atexit.register(shutil.rmtree, root, ignore_errors=True)\n", []),
    ("import tempfile, shutil\nroot = tempfile.mkdtemp()\nshutil.rmtree(root)\n",
     ["never removed"]),
    ("import tempfile\nfd, p = tempfile.mkstemp()\n", ["mkstemp"]),
    ("import tempfile\nf = tempfile.NamedTemporaryFile(delete=False)\n",
     ["delete=False"]),
])
def test_the_rule_tells_a_cleaned_script_from_a_leaking_one(source, expected):
    """An rmtree OUTSIDE a finally is skipped by an exception, so it does not count."""

    found = scratch_findings(source)
    assert len(found) == len(expected)
    for message, fragment in zip(found, expected, strict=True):
        assert fragment in message
