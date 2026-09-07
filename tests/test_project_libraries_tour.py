"""The runnable tour of the Project libraries surface.

The repo's other notebook tests inspect SOURCE. This one executes, because the
notebook's whole claim is that the surface works end to end -- a tour that does
not run is worse than no tour. It is small on purpose: a two-cell model, two MF6
runs, a few seconds.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "examples" / "mf6" / "notebooks" / "project_libraries_tour.ipynb"


def _code_cells() -> list[str]:
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return [
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    ]


def test_the_tour_runs_top_to_bottom(tmp_path, monkeypatch):
    """Execute every code cell in order, in one namespace, as a reader would."""

    monkeypatch.chdir(NOTEBOOK.parent)
    namespace: dict[str, object] = {"__name__": "__main__"}

    for index, source in enumerate(_code_cells(), start=1):
        try:
            exec(compile(source, f"<tour cell {index}>", "exec"), namespace)
        except Exception as error:  # noqa: BLE001 -- re-raised with the cell text
            raise AssertionError(
                f"tour cell {index} failed: {type(error).__name__}: {error}\n\n{source}"
            ) from error

    # The tour ends having exercised the whole grammar, not just imported it.
    project = namespace["project"]
    assert sorted(project.packages) == [
        "chd/base", "dis/base", "ic/steady", "npf/base", "npf/high_k", "oc/all",
    ]
    assert sorted(project.simulations) == ["baseline", "high_k"]
    assert sorted(project.refresh_runs()) == ["baseline", "high_k"]


def test_the_tour_teaches_the_distinctions_that_cost_the_most(tmp_path, monkeypatch):
    """The three that are easy to get wrong, and each is shown rather than told."""

    monkeypatch.chdir(NOTEBOOK.parent)
    namespace: dict[str, object] = {"__name__": "__main__"}
    for source in _code_cells():
        exec(compile(source, "<tour>", "exec"), namespace)

    project = namespace["project"]
    drift = project.refresh_runs().drift().set_index(["run", "key"])

    # 1. drift is per-key AND per-run: editing npf/base reaches baseline, and
    #    leaves high_k alone because it references npf/high_k.
    assert drift.loc[("baseline", "npf/base"), "state"] == "edited"
    assert drift.loc[("high_k", "npf/high_k"), "state"] == "unchanged"

    # 2. status is myflopy's verb, outcome is MF6's verdict -- two columns.
    summary = project.runs.summary()
    assert {"status", "outcome"} <= set(summary.columns)

    # 3. a calibration is reported, not listed as runs.
    assert len(project.runs.calibrations) == 1
    assert "calibration excluded" in repr(project.runs)


def test_the_tour_leaves_nothing_behind():
    """It builds in a temp directory, so running it does not litter the repo."""

    source = "\n".join(_code_cells())

    assert "tempfile.mkdtemp()" in source
    assert not (NOTEBOOK.parent / "tour").exists()


def test_the_notebook_ships_without_stored_output():
    """Outputs are stripped, like the other tours -- a diff of a notebook should
    be a diff of its prose and code."""

    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))

    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            assert cell.get("outputs") == []
            assert cell.get("execution_count") is None
