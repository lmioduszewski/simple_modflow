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
- ``docs/reproducers/pest_template_holds_inputs_only.py`` lifted a fixed list of
  names out of ``pest/project.py``, which stopped matching when a helper was split
  out of ``_simulation_input_files``.
- ``layer_management_workflow.ipynb`` built a ``Contours`` surface with no region,
  which ``Surface.contours`` has refused since the region became required.
- ``plotting_vocabulary_tour.ipynb`` carried a mangled line naming ``built`` before
  the cell defined it; ``refined_end_to_end_preferred_api_workflow.ipynb`` printed a
  ``connectiondata`` it no longer defines; ``simple_model_workflow.ipynb`` and
  ``triangle_voronoi_simplemodel_workflow.ipynb`` read a ``geometry`` column the head
  table no longer has, and the second built through the removed
  ``build_simple_model``/``mf_folder_path``; ``feature_rich_model_workflow.ipynb``
  called FloPy's ``model.gwf.plot()``, which FloPy 3.10 cannot draw for an SFR
  period with no settings.

Each runs as a reader would: the scripts in a subprocess, the notebook's code cells
in order in one namespace, with each cell's last expression and every ``display``
argument RENDERED as Jupyter would -- a picture assembles its figure lazily, so
executing the code alone draws nothing. A script's TMPDIR is a directory of its own,
and a run that leaves anything there FAILS naming it (ledger 460): /tmp is a
quota'd tmpfs, and fifteen scripts left their ``mkdtemp()`` behind (12 MB a run for
``flow_targets_period_end.py``). ``tests/test_scripts_clean_up_their_scratch.py`` holds every script
to the same rule statically.
"""

from __future__ import annotations

import ast
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

ROOT = Path(__file__).resolve().parents[1]


#: What a run may leave in its TMPDIR: pyEMU's ``PortManager`` makes ``port_locks``
#: once, under a fixed name (ledger 335's allowance).
_TMP_ALLOWED = {"port_locks"}


def _run_script(script: Path, tmp_path: Path, *, timeout: float = 600.0) -> str:
    """Run ``script`` with this checkout's ``src`` first, failing if it leaves temp files.

    The working directory is ``tmp_path`` and TMPDIR a directory of its own beside
    it, so anything the script -- or the library under it -- leaves in the temp root
    is a named failure, not a quiet use of the user's quota'd /tmp.
    """

    temp_root = tmp_path / "TMPDIR"
    temp_root.mkdir()
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    env["TMPDIR"] = str(temp_root)
    env["MPLBACKEND"] = "Agg"
    done = subprocess.run(
        [sys.executable, str(script)], cwd=tmp_path, env=env, capture_output=True,
        text=True, timeout=timeout,
    )
    assert done.returncode == 0, (
        f"{script.relative_to(ROOT)} exited {done.returncode}\n"
        f"--- stdout (tail)\n{done.stdout[-3000:]}\n--- stderr (tail)\n{done.stderr[-3000:]}"
    )
    left = sorted(entry.name for entry in temp_root.iterdir() if entry.name not in _TMP_ALLOWED)
    assert left == [], f"{script.relative_to(ROOT)} left {left} in its TMPDIR"
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


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
def test_the_template_reproducer_reads_this_checkouts_input_rule(tmp_path):
    """It lifts ``_simulation_input_files`` out of the source by what that function uses.

    A fixed list of names stopped it with a ``NameError`` once the walk was split out
    into ``_scan_simulation_inputs``.
    """

    pytest.importorskip("pyemu")
    out = _run_script(
        ROOT / "docs" / "reproducers" / "pest_template_holds_inputs_only.py", tmp_path
    )
    assert "  input          134 / 0.52 MB" in out
    assert "  run's other      0 / 0.00 MB" in out


NOTEBOOKS = ROOT / "examples" / "mf6" / "notebooks"

#: Notebook -> a check on its final namespace. Each is fast (3-15 s) and needs only
#: ``mf6``; the canonical and PEST notebooks are left out for their run time (ledger
#: 444), and the archive ones do not run by design.
EXECUTED_NOTEBOOKS = {
    "layer_management_workflow": lambda ns: ns["demo"].cache_status() == {"bedrock": "missing"},
    "plotting_vocabulary_tour": lambda ns: ns["built"].nlay == 2,
    "feature_rich_model_workflow": lambda ns: True,
    "refined_end_to_end_preferred_api_workflow": lambda ns: len(ns["lak"].connectiondata) > 0,
    "simple_model_workflow": lambda ns: list(ns["completed"].all_heads.columns) == ["elev"],
    "triangle_voronoi_simplemodel_workflow": lambda ns: ns["run"].success,
}


def _render(obj) -> None:
    """Draw ``obj`` the way Jupyter would display it, letting any failure raise.

    A myflopy picture assembles its figure LAZILY (``Picture.fig``), so a cell whose
    map is its last expression -- or is handed to ``display`` -- draws nothing until
    the front end asks for a representation. Executing the code alone left every such
    cell untested: with every map's assembly made to raise, five of the six notebooks
    still passed, ``feature_rich_model_workflow`` among them, whose last cell IS a map
    (review of ledger 444). IPython's own formatters print and swallow an
    exception in a ``_repr_*_`` method, so the methods are called directly here.
    """

    if obj is None:
        return
    for method in ("_repr_mimebundle_", "_repr_html_", "_repr_png_", "_repr_svg_"):
        if callable(getattr(obj, method, None)):
            getattr(obj, method)()
            return


def _run_cell(source: str, filename: str, namespace: dict) -> None:
    """Run one code cell, then render its last expression as Jupyter's display would."""

    tree = ast.parse(source, filename)
    last = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
    exec(compile(tree, filename, "exec"), namespace)
    if last is not None:
        _render(eval(compile(ast.Expression(last.value), filename, "eval"), namespace))


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
@pytest.mark.parametrize("name", sorted(EXECUTED_NOTEBOOKS))
def test_the_notebook_runs_top_to_bottom(name, tmp_path, monkeypatch):
    """Every code cell, in order, in one namespace, as a reader would run it.

    The working directory is a stand-in ``examples/mf6/notebooks`` under ``tmp_path``,
    so what a notebook writes beside itself (``../artifacts``, ``./artifacts``) lands
    there; the real notebooks directory is on ``sys.path`` for its local modules.
    """

    cells = [
        "".join(cell["source"])
        for cell in json.loads((NOTEBOOKS / f"{name}.ipynb").read_text(encoding="utf-8"))["cells"]
        if cell["cell_type"] == "code"
    ]
    here = tmp_path / "examples" / "mf6" / "notebooks"
    here.mkdir(parents=True)
    monkeypatch.chdir(here)
    monkeypatch.syspath_prepend(str(NOTEBOOKS))

    def display(*objs, **_options):
        for obj in objs:
            _render(obj)

    import IPython.display
    import plotly.io

    # `from IPython.display import display` inside a notebook must get the renderer
    # too (outside a kernel IPython's own just prints the repr); and `fig.show()`
    # serializes the figure instead of opening a browser tab -- plotly's default
    # renderer outside Jupyter is "browser".
    monkeypatch.setattr(IPython.display, "display", display)
    monkeypatch.setattr(plotly.io, "show", lambda fig, *a, **k: fig.to_plotly_json())
    namespace: dict[str, object] = {"__name__": "__main__", "display": display}
    for index, source in enumerate(cells, start=1):
        try:
            _run_cell(source, f"<{name} cell {index}>", namespace)
        except Exception as error:  # noqa: BLE001 -- re-raised naming the cell and its text
            raise AssertionError(
                f"{name}: code cell {index} failed: {type(error).__name__}: {error}\n\n{source}"
            ) from error
        finally:
            plt.close("all")
    assert EXECUTED_NOTEBOOKS[name](namespace)
