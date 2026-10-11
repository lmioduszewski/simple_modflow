"""A reproducer that runs PEST++ makes ``PYTHONPATH`` absolute first (ledger 425).

PEST++ starts ``python forward_run.py`` in the template. Under a relative
``PYTHONPATH`` (``PYTHONPATH=src``, the way every reproducer documents its "before"
column) that child resolved ``src`` against the template, found nothing, and imported
the venv's editable install -- the main checkout -- so the reproducer measured the
wrong tree; a first run of ``lake_default_prefix.py`` exited 1 for that reason. Ledger
354 fixed it for pytest only. The one rewrite now lives in
``docs/reproducers/_absolute_pythonpath.py``, which ``tests/conftest.py`` loads too.
"""

from __future__ import annotations

import ast
import os
import re
from pathlib import Path

import pytest

from tests.conftest import _absolute_pythonpath

REPRODUCERS = Path(__file__).resolve().parents[1] / "docs" / "reproducers"

#: Starts PEST++ -- a ``pestpp-*`` executable, or the calibration's own runners.
_RUNS_PEST = re.compile(r"\"pestpp|pestpp-(ies|glm)\"|\.prior\(|\.run_ies\(")
_IMPORTS_MYFLOPY = re.compile(r"^(import myflopy|from myflopy)", re.MULTILINE)


def _runs_pest_with_myflopy(source: str) -> bool:
    return bool(_RUNS_PEST.search(source)) and bool(_IMPORTS_MYFLOPY.search(source))


def _calls_at_module_level(source: str, name: str) -> bool:
    return any(
        isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Name)
        and node.value.func.id == name
        for node in ast.parse(source).body
    )


def test_every_reproducer_running_pest_makes_pythonpath_absolute():
    scripts = {
        path.name: path.read_text(encoding="utf-8")
        for path in sorted(REPRODUCERS.glob("*.py"))
        if not path.name.startswith("_")
    }
    running = sorted(name for name, source in scripts.items() if _runs_pest_with_myflopy(source))
    assert {"lake_default_prefix.py", "pest_handed_targets.py", "pest_observation_groups.py",
            "obs_route_period_end.py"} <= set(running)
    missing = [name for name in running
               if not _calls_at_module_level(scripts[name], "make_pythonpath_absolute")]
    assert not missing, (
        f"These reproducers run PEST++ and never call make_pythonpath_absolute() at module "
        f"level, so a relative PYTHONPATH imports another myflopy in the forward run: {missing}"
    )


def test_the_rewrite_keeps_absolute_and_order_and_reads_empty_as_the_start(monkeypatch, tmp_path):
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("_absolute_pythonpath", REPRODUCERS / "_absolute_pythonpath.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    value = os.pathsep.join(["src", "/abs/lib", "", "../other"])
    monkeypatch.setenv("PYTHONPATH", value)
    rewritten = module.make_pythonpath_absolute(str(tmp_path))
    assert rewritten.split(os.pathsep) == [
        str(tmp_path / "src"), "/abs/lib", str(tmp_path), str(tmp_path.parent / "other")]
    assert os.environ["PYTHONPATH"] == rewritten
    monkeypatch.delenv("PYTHONPATH")
    assert module.make_pythonpath_absolute(str(tmp_path)) is None
    assert "PYTHONPATH" not in os.environ
    # tests/conftest.py rewrites with the same function.
    assert _absolute_pythonpath(value, str(tmp_path)) == rewritten


@pytest.mark.parametrize("name", ["lake_default_prefix.py", "pest_handed_targets.py"])
def test_the_guard_reads_a_reproducer_without_the_call_as_missing(name):
    source = (REPRODUCERS / name).read_text(encoding="utf-8")
    assert _calls_at_module_level(source, "make_pythonpath_absolute")
    stripped = source.replace("\nmake_pythonpath_absolute()\n", "\n")
    assert not _calls_at_module_level(stripped, "make_pythonpath_absolute")
