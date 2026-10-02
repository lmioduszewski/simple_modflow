"""A subprocess a test starts in ANOTHER directory imports the myflopy under test (ledger 354).

Python resolves a relative ``PYTHONPATH`` entry against the directory a process starts
in. Run from a worktree with ``PYTHONPATH=src``, the session imports the worktree's
``src``, but a child started elsewhere -- PEST's ``forward_run.py`` runs with its
TEMPLATE as working directory -- resolved ``src`` there, found nothing, and imported
the venv's EDITABLE install, the main checkout: a change inside the library was
exercised, in every forward run, by code that did not contain it (ledger 344).
``tests/conftest.py`` now rewrites every relative entry to an absolute path when the
session starts, so the absolute ``PYTHONPATH`` CLAUDE.md asked for is no longer the
caller's job. These tests fail without that rewrite when the suite is run from a
worktree with a relative ``PYTHONPATH``; with an absolute one, or none, they pass
either way, which is the point: the session and its children agree.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import myflopy
from tests.conftest import _absolute_pythonpath


def test_a_subprocess_started_elsewhere_imports_the_myflopy_under_test(tmp_path):
    """A child started in an empty directory imports the session's own myflopy."""

    result = subprocess.run(
        [sys.executable, "-c", "import myflopy; print(myflopy.__file__)"],
        cwd=tmp_path, capture_output=True, text=True, timeout=300,
    )
    assert result.returncode == 0, result.stderr
    child = Path(result.stdout.strip().splitlines()[-1]).resolve()
    assert child == Path(myflopy.__file__).resolve(), (
        f"a subprocess started in {tmp_path} imported {child}, the session imports "
        f"{Path(myflopy.__file__).resolve()} -- PYTHONPATH is "
        f"{os.environ.get('PYTHONPATH')!r} (ledger 354)"
    )


def test_every_pythonpath_entry_the_session_hands_on_is_absolute():
    """What every child of the session inherits names no directory relative to its cwd."""

    value = os.environ.get("PYTHONPATH")
    entries = value.split(os.pathsep) if value else []
    assert [entry for entry in entries if not os.path.isabs(entry)] == []


def test_a_relative_entry_is_resolved_as_python_resolves_it(tmp_path):
    """Against the start directory, an empty entry meaning that directory, the order
    and every absolute entry kept -- what the interpreter itself puts on ``sys.path``."""

    start = tmp_path / "start"
    (start / "src").mkdir(parents=True)
    value = os.pathsep.join(["src", "/opt/abs", "", "../up"])
    assert _absolute_pythonpath(value, str(start)).split(os.pathsep) == [
        str(start / "src"), "/opt/abs", str(start), str(tmp_path / "up")]

    shown = subprocess.run(
        [sys.executable, "-c", "import sys; print('\\n'.join(sys.path))"],
        cwd=start, capture_output=True, text=True, timeout=300,
        env={**os.environ, "PYTHONPATH": value},
    )
    resolved = shown.stdout.splitlines()
    for entry in (str(start / "src"), "/opt/abs", str(tmp_path / "up")):
        assert entry in resolved, (entry, resolved[:6])
