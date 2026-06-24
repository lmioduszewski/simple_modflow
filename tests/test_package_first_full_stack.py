"""Integration test for the package-first full-stack reference example.

Proves the package-first API assembles GHB + RCH + DRN + UZF + SFR + LAK + MVR
together on a multi-layer DISV grid inside a Project, writes every package, and
runs MF6 to convergence -- the end-to-end "does it actually run" check.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
EXAMPLES = ROOT / "examples" / "mf6"
if str(EXAMPLES) not in sys.path:
    sys.path.insert(0, str(EXAMPLES))


@pytest.mark.slow
def test_package_first_full_stack_builds_and_runs(tmp_path):
    pytest.importorskip("flopy")
    from package_first_full_stack import build_full_stack_project

    project, _ = build_full_stack_project(tmp_path / "proj")
    run = project.prepare_run("baseline", "baseline")

    # Every package assembled through the package-first API, in valid MVR order
    # (the mover is declared after the sfr/lak packages it moves between).
    packages = list(run.built.models["valley"].packages)
    for name in ("disv", "npf", "ic", "sto", "ghb", "drn", "rch", "uzf", "sfr", "lak", "mvr", "oc"):
        assert name in packages, f"missing package {name!r} (have {packages})"
    assert packages.index("mvr") > packages.index("sfr")
    assert packages.index("mvr") > packages.index("lak")

    # Writes + runs MF6 to convergence.
    success, report = run.execute()
    assert success, "MF6 did not converge:\n" + "\n".join(report[-25:])

    # Heads, budget, and the advanced-package input files are on disk.
    assert (run.workspace / "valley.hds").exists()
    assert (run.workspace / "valley.cbc").exists()
    for suffix in ("ghb", "drn", "uzf", "sfr", "lak", "mvr"):
        assert (run.workspace / f"valley.{suffix}").exists(), f"missing valley.{suffix}"
