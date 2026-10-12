"""PEST++ is handed an ABSOLUTE ``PYTHONPATH`` by ``run_ies``/``prior`` (ledger 454).

Ledger 354 made the test suite's ``PYTHONPATH`` absolute and ledger 425 the
reproducers', but ``PestProject.run_ies``/``prior`` handed PEST++ the environment as
given: a notebook or script started under ``PYTHONPATH=src`` had PEST++ start
``forward_run.py`` in the template, where ``src`` names nothing, and the forward run
imported the venv's editable install -- another checkout -- instead of the code the
calibration was built with. The rewrite is the one ``tests/conftest.py`` and the
reproducers use (``myflopy._pythonpath``); a relative entry the working directory has
moved away from since Python started is refused rather than guessed.
"""

from __future__ import annotations

import os
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

import myflopy as mf
from myflopy._pythonpath import pythonpath_for_children

pyemu = pytest.importorskip("pyemu")

#: The directory this session's ``src`` sits in, and ``src``'s name relative to it.
SRC = Path(mf.__file__).resolve().parents[1]


class _Launched(Exception):
    """Raised by the stand-in for PEST++ once it has looked at its environment."""


@pytest.fixture
def calibration(tmp_path):
    from tests.test_pest_structured_grid import _dis_run

    model = mf.load_mf6_run(_dis_run(tmp_path / "run"))
    cal = model.pest("pp", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("chd", style="constant", bounds=(0.9, 1.1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("pp.pst", noptmax=0)
    return cal


def _stand_in_for_pestpp(seen):
    """``pyemu.os_utils.run``'s stand-in: record the environment, and what a Python
    started in the template -- as PEST++ starts the forward run -- imports."""

    def run(command, cwd=".", **_kwargs):
        child = subprocess.run(
            [sys.executable, "-c", "import myflopy; print(myflopy.__file__)"],
            cwd=cwd, capture_output=True, text=True, timeout=300,
        )
        seen.update(pythonpath=os.environ.get("PYTHONPATH"), child=child)
        raise _Launched

    return run


def test_pestpp_is_handed_an_absolute_pythonpath_and_imports_this_myflopy(
    calibration, monkeypatch
):
    """Under ``PYTHONPATH=src`` from the directory holding ``src``, PEST++ gets the
    absolute path, a Python started in the template imports the session's own myflopy,
    and the caller's environment is restored afterwards."""

    monkeypatch.chdir(SRC.parent)
    monkeypatch.setenv("PYTHONPATH", SRC.name)
    seen = {}
    monkeypatch.setattr(pyemu.os_utils, "run", _stand_in_for_pestpp(seen))
    with pytest.raises(_Launched):
        calibration.prior(reals=2)

    assert seen["pythonpath"] == str(SRC)
    assert seen["child"].returncode == 0, seen["child"].stderr
    child = Path(seen["child"].stdout.strip().splitlines()[-1]).resolve()
    assert child == Path(mf.__file__).resolve()
    assert os.environ["PYTHONPATH"] == SRC.name


def test_a_relative_entry_the_working_directory_has_left_is_refused(
    calibration, monkeypatch, tmp_path
):
    """Read against a directory the process did not start in, ``src`` names a
    directory Python never put on ``sys.path``: refused, naming the entry, before
    PEST++ starts."""

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    monkeypatch.setenv("PYTHONPATH", "src")
    seen = {}
    monkeypatch.setattr(pyemu.os_utils, "run", _stand_in_for_pestpp(seen))
    with pytest.raises(RuntimeError, match=r"(?s)PYTHONPATH entry 'src' is relative.*"
                                           r"not on sys.path.*ABSOLUTE"):
        calibration.prior(reals=2)
    assert seen == {}


def test_an_absolute_or_absent_pythonpath_is_handed_on_as_it_is(monkeypatch, tmp_path):
    monkeypatch.setenv("PYTHONPATH", os.pathsep.join([str(SRC), "/opt/abs"]))
    assert pythonpath_for_children() == os.pathsep.join([str(SRC), "/opt/abs"])
    monkeypatch.delenv("PYTHONPATH")
    assert pythonpath_for_children() is None
