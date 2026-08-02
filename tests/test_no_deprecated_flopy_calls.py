"""FloPy APIs we have been told are going away, and must stop calling.

FloPy marks a method deprecated long before it removes it, and the warning is
easy to live with -- which is exactly the problem. `model.gwf.package_names`
warned on **every canonical model build** for long enough that the line became
part of the expected test output, invisible.

The rule is narrow on purpose: this is not "no deprecated calls anywhere", which
would be unenforceable against a dependency we do not control. It is a list of
specific APIs, each with the supported replacement, added when we actually hit
one.

Two implementation notes, both learned by getting it wrong first:

* It reads the **AST**, not the text. A regex over source lines flagged the
  comment that explains the rule, in the very file the rule fixed.
* It scans for the access on a FLOPY object specifically (``<x>.gwf.…`` /
  ``<x>.sim.…``). myflopy's own ``SimulationBase.package_names`` is fine -- it
  *is* the replacement -- and the last test here checks that claim rather than
  trusting it.
"""

from __future__ import annotations

import ast
import inspect
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCANNED = ("src", "tests", "examples", "scripts")

#: ``(flopy handle, attribute) -> (what is wrong, what to use instead)``.
FORBIDDEN = {
    ("gwf", "package_names"): (
        "flopy deprecated MFModel.package_names in 3.9 and warns on every access",
        "use myflopy's own `model.package_names` (uppercased names), or "
        "`model.gwf.get_package_list()` directly",
    ),
    ("sim", "package_names"): (
        "the same deprecation, reached through the simulation handle",
        "use myflopy's own `model.package_names`",
    ),
    ("gwf", "package_name_dict"): (
        "flopy deprecated MFModel.package_name_dict in 3.9 -- and note that even "
        "`getattr(gwf, 'package_name_dict', {})` warns, because reading the "
        "attribute is what warns",
        "use `model.gwf.get_package_list()` (same names, uppercased)",
    ),
    ("gwf", "package_type_dict"): (
        "flopy deprecated MFModel.package_type_dict in 3.9",
        "use `model.gwf.get_package_list(ftype=...)`",
    ),
}


def _source_files():
    for folder in SCANNED:
        base = ROOT / folder
        if not base.exists():
            continue
        for path in sorted(base.rglob("*.py")):
            if "_vendor" not in path.parts:
                yield path


def _accesses(path: Path):
    """Yield ``(handle, attribute, lineno)`` for every ``<x>.<handle>.<attr>``."""

    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Attribute):
            continue
        if not isinstance(node.value, ast.Attribute):
            continue
        yield node.value.attr, node.attr, node.lineno


@pytest.mark.parametrize(
    "handle,attribute", sorted(FORBIDDEN), ids=lambda part: str(part)
)
def test_a_deprecated_flopy_api_is_not_called(handle, attribute):
    reason, replacement = FORBIDDEN[(handle, attribute)]

    hits = [
        f"{path.relative_to(ROOT)}:{line}"
        for path in _source_files()
        for found_handle, found_attribute, line in _accesses(path)
        if (found_handle, found_attribute) == (handle, attribute)
    ]
    assert not hits, f"{hits}: {reason}. Instead, {replacement}."


def test_our_own_package_names_uses_the_supported_call():
    """The replacement has to actually BE a replacement. If this property were
    itself implemented on the deprecated one, the scan above would be theatre --
    every call site would look clean while the warning still fired."""

    from myflopy.modflow.mf6.simulation.base import SimulationBase

    source = textwrap.dedent(inspect.getsource(SimulationBase.package_names.fget))
    body = ast.parse(source)
    attributes = {
        node.attr for node in ast.walk(body) if isinstance(node, ast.Attribute)
    }
    assert "get_package_list" in attributes
    assert "package_names" not in attributes
