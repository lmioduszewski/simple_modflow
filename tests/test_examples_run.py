"""Examples, reproducers and notebooks that the library broke once, executed.

Each of these stopped running at some HEAD while every test passed, because the
suite read their SOURCE or not at all (ledger 444):

- ``examples/mf6/contours_to_model.py`` -- the runnable form of
  ``docs/model_building_cheatsheet.md`` -- read ``vor.gdf_topbtm["top"]``, while
  ``LayerStack``'s ``attach=True`` publishes INTEGER columns (0 = model top).
- ``docs/reproducers/top_active_passthrough.py`` declared the old ``irch`` through
  myflopy, which ledger 321's build guard now refuses; its "before" is measured on
  files FloPy writes instead.
- ``docs/reproducers/obs_route_declared_identity.py`` read ``model.outputs.lak`` on
  a model of two LAK packages, which ledger 352 refuses, naming
  ``model.outputs.lak(package=...)``.
- ``layer_management_workflow.ipynb`` built a ``Contours`` surface with no region,
  which ``Surface.contours`` has refused since the region became required.

Each runs as a reader would: the scripts in a subprocess, the notebook's code cells
in order in one namespace. The scripts write under ``tempfile.mkdtemp()``, so TMPDIR
points at the test's own directory and its teardown removes what they leave.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _run_script(script: Path, tmp_path: Path, *, timeout: float = 600.0) -> str:
    """Run ``script`` with this checkout's ``src`` first and its temp root in ``tmp_path``."""

    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    env["TMPDIR"] = str(tmp_path)
    env["MPLBACKEND"] = "Agg"
    done = subprocess.run(
        [sys.executable, str(script)], cwd=tmp_path, env=env, capture_output=True,
        text=True, timeout=timeout,
    )
    assert done.returncode == 0, (
        f"{script.relative_to(ROOT)} exited {done.returncode}\n"
        f"--- stdout (tail)\n{done.stdout[-3000:]}\n--- stderr (tail)\n{done.stderr[-3000:]}"
    )
    return done.stdout


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("grass") is None and not os.environ.get("GRASS_BIN"),
                    reason="mf.Contours needs GRASS GIS")
@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
def test_the_cheatsheet_example_runs_to_a_reloaded_model(tmp_path):
    out = _run_script(ROOT / "examples" / "mf6" / "contours_to_model.py", tmp_path)
    assert "converged: True" in out
    assert "reloaded from disk and rebuilt: nlay=3" in out


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
def test_the_top_active_reproducer_still_measures_the_old_rule(tmp_path):
    """The OLD irch is refused at build, and MF6 still drops 42.2% of it when written."""

    out = _run_script(ROOT / "docs" / "reproducers" / "top_active_passthrough.py", tmp_path)

    def line(start: str) -> str:
        return next(text for text in out.splitlines() if text.startswith(start))

    refused = line("OLD  ")
    assert "REFUSED at build" in refused
    assert "152.017 of 360 ft3/d (42.2%)" in refused
    assert "applied   207.98 ( 57.8%)" in line("OLD (FloPy, MF6)")
    assert "applied   360.00 (100.0%)" in line("LAZY ")


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
def test_the_declared_identity_reproducer_names_the_lak_package(tmp_path):
    out = _run_script(
        ROOT / "docs" / "reproducers" / "obs_route_declared_identity.py", tmp_path
    )
    assert "m.lak.stg (package 'lak'): [111.693, 110.705, 110.281]" in out


def test_the_layer_management_notebook_runs_top_to_bottom(tmp_path, monkeypatch):
    """Every code cell, in order, in one namespace; no GRASS (cell 18 only asks the cache)."""

    notebook = ROOT / "examples" / "mf6" / "notebooks" / "layer_management_workflow.ipynb"
    cells = [
        "".join(cell["source"])
        for cell in json.loads(notebook.read_text(encoding="utf-8"))["cells"]
        if cell["cell_type"] == "code"
    ]
    monkeypatch.chdir(tmp_path)                 # it writes under ./artifacts
    namespace: dict[str, object] = {"__name__": "__main__"}
    for index, source in enumerate(cells, start=1):
        try:
            exec(compile(source, f"<layer_management cell {index}>", "exec"), namespace)
        except Exception as error:  # noqa: BLE001 -- re-raised naming the cell and its text
            raise AssertionError(
                f"code cell {index} failed: {type(error).__name__}: {error}\n\n{source}"
            ) from error
    assert namespace["demo"].cache_status() == {"bedrock": "missing"}
