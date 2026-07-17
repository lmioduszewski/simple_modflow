"""Optional dependencies must not be required by the core import surface.

``dash``/``dash_bootstrap_components`` (the ``viz`` extra) and ``osgeo``/GDAL
(a system library) are import-time optional: ``import myflopy`` and resolving
every lazy export in ``myflopy.__all__`` must succeed with them absent. Each
test runs a subprocess whose import machinery blocks the module, simulating a
fresh install without the extra.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

_BLOCKER_TEMPLATE = """
import importlib.abc
import sys

BLOCKED = {blocked!r}


class _Blocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        root = fullname.split(".")[0]
        if root in BLOCKED:
            raise ImportError(f"{{fullname}} is blocked for this test")
        return None


sys.meta_path.insert(0, _Blocker())

import myflopy

for name in myflopy.__all__:
    getattr(myflopy, name)

# The lazy-export surface beyond __all__ (the unwarned engine tier) must
# also stay importable.
for name in getattr(myflopy, "__engine__", ()):
    getattr(myflopy, name)

print("ok")
"""


def _run_with_blocked_modules(*blocked: str) -> subprocess.CompletedProcess:
    code = _BLOCKER_TEMPLATE.format(blocked=tuple(blocked))
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env={**os.environ, "PYTHONPATH": str(SRC), "MPLBACKEND": "Agg"},
        capture_output=True,
        text=True,
    )


def _assert_ok(result: subprocess.CompletedProcess, blocked: str) -> None:
    assert result.returncode == 0 and "ok" in result.stdout, (
        f"import surface requires optional dependency {blocked!r}:\n"
        f"{result.stdout[-1000:]}\n{result.stderr[-2000:]}"
    )


def test_import_surface_works_without_dash():
    _assert_ok(_run_with_blocked_modules("dash", "dash_bootstrap_components"), "dash")


def test_import_surface_works_without_osgeo():
    _assert_ok(_run_with_blocked_modules("osgeo"), "osgeo")


def test_import_surface_works_without_viz3d():
    _assert_ok(
        _run_with_blocked_modules("pyvista", "trame", "trame_vtk", "trame_vuetify"),
        "pyvista/trame",
    )


def test_import_surface_works_without_pyemu():
    _assert_ok(_run_with_blocked_modules("pyemu"), "pyemu")
