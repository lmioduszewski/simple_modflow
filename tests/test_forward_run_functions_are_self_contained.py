"""Every function the PEST build hands pyEMU binds every module it uses INSIDE its own body.

pyEMU's ``PstFrom.add_py_function`` copies one function's SOURCE into the generated
``forward_run.py``, so a module-level name of ``pest/forward_run.py`` it reads is not
there (``forward_run.read_target_table``'s docstring states the rule). The R2-screens
batch's lint pass removed ``apply_pilotpoints_to_array``'s own ``import numpy as np``
in favour of the module's, which the copied function only found because pyEMU's
header happens to import numpy as ``np`` too.

The functions are DISCOVERED from the build's ``add_py_function`` calls in the source,
not listed, so a new forward-run function is covered the day it is wired.
"""

from __future__ import annotations

import builtins
import dis
import re
import types
from pathlib import Path

import pytest

from myflopy.modflow.mf6.pest import forward_run

SRC = Path(forward_run.__file__).resolve().parents[3]
CALL = re.compile(r"add_py_function\(\s*[^,]+,\s*f?[\"'](\w+)\(", re.S)


def _wired() -> list[str]:
    names = set()
    for path in SRC.rglob("*.py"):
        names.update(CALL.findall(path.read_text(encoding="utf-8")))
    return sorted(names)


def _globals_read(code: types.CodeType) -> set[str]:
    found = {ins.argval for ins in dis.get_instructions(code)
             if ins.opname in {"LOAD_GLOBAL", "LOAD_NAME"}}
    for const in code.co_consts:
        if isinstance(const, types.CodeType):
            found |= _globals_read(const)
    return found


def test_the_wired_functions_are_found():
    assert {"apply_pilotpoints_to_array", "_write_head_target_csv",
            "_write_conc_target_csv", "write_named_series_targets"} <= set(_wired())


@pytest.mark.parametrize("name", _wired())
def test_a_forward_run_function_reads_no_module_level_name(name):
    function = getattr(forward_run, name)
    stray = sorted(n for n in _globals_read(function.__code__) if not hasattr(builtins, n))
    assert not stray, (
        f"forward_run.{name} reads module-level {stray}; pyEMU copies only its source "
        "into the generated forward run, so import them inside its body."
    )
