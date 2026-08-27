"""Every lazily-exported name must also be statically visible.

`myflopy/__init__.py` resolves its public surface through a module-level
`__getattr__` over `_EXPORTS`, so `import myflopy as mf; mf.Raster(...)` works at
runtime without importing the world. But **PyCharm and Pylance never run
`__getattr__`** -- they are static, so a name that exists only in `_EXPORTS` is
invisible: no completion, no go-to-definition, no type checking. It looks to the
user exactly like the export is missing.

The `if TYPE_CHECKING:` block exists to close that gap, and it is hand-written,
so it drifts. It had drifted: `Raster`, `Contours`, `Points` and the rest of the
surface-constructor family were added to `_EXPORTS` on 2026-08-20 and never added
here, along with `ZoneSpec`, `ConcObservationSpec`, `lak_connection` and
`sfr_connection` -- 16 names in all, reported by a user who could not find
`mf.Raster` in PyCharm (2026-08-27).

These tests parse the source rather than importing it: `TYPE_CHECKING` is False at
runtime, so the block's names never exist to introspect.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

INIT = Path(__file__).resolve().parents[1] / "src" / "myflopy" / "__init__.py"

#: Handled by `__getattr__`'s explicit subpackage branch rather than `_EXPORTS`,
#: but still needing static declaration -- `from myflopy import plot` is the
#: documented plotting front door.
SUBPACKAGES = {"modflow", "project", "plot"}


@pytest.fixture(scope="module")
def tree() -> ast.Module:
    return ast.parse(INIT.read_text())


def _type_checking_names(tree: ast.Module) -> set[str]:
    """Names imported under `if TYPE_CHECKING:` at module level."""

    names: set[str] = set()
    for node in tree.body:
        if not (isinstance(node, ast.If) and ast.unparse(node.test) == "TYPE_CHECKING"):
            continue
        for sub in node.body:
            if isinstance(sub, (ast.Import, ast.ImportFrom)):
                names |= {a.asname or a.name for a in sub.names}
    return names


def _export_names(tree: ast.Module) -> set[str]:
    """Keys of the `_EXPORTS` mapping."""

    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "_EXPORTS" for t in node.targets
        ):
            return {
                key.value
                for key in node.value.keys
                if isinstance(key, ast.Constant) and isinstance(key.value, str)
            }
    raise AssertionError("no _EXPORTS mapping found in myflopy/__init__.py")


def test_every_lazy_export_is_statically_visible(tree):
    """The failure this catches is silent: the name works, the IDE just can't see it."""

    missing = sorted((_export_names(tree) | SUBPACKAGES) - _type_checking_names(tree))
    assert not missing, (
        "these names resolve at runtime but are INVISIBLE to PyCharm/Pylance -- "
        f"add them to the `if TYPE_CHECKING:` block in myflopy/__init__.py: {missing}"
    )


def test_the_type_checking_block_declares_nothing_extra(tree):
    """The reverse direction: a static import for a name that no longer exports.

    Harmless at runtime -- `TYPE_CHECKING` is False -- which is exactly why it
    would sit there indefinitely, promising an attribute that raises
    AttributeError when anyone follows the completion.
    """

    stale = sorted(_type_checking_names(tree) - _export_names(tree) - SUBPACKAGES)
    assert not stale, (
        "the TYPE_CHECKING block imports names that are not in `_EXPORTS`, so an "
        f"editor offers them and `mf.<name>` then raises AttributeError: {stale}"
    )


@pytest.mark.parametrize("name", sorted(SUBPACKAGES))
def test_the_declared_subpackages_actually_import(name):
    """`from myflopy import plot` must work, not merely type-check."""

    import myflopy

    assert getattr(myflopy, name) is not None


def test_the_statically_declared_names_resolve_at_runtime(tree):
    """A TYPE_CHECKING import pointing at the wrong module is invisible until
    someone actually reaches for the attribute. Resolve them all here instead."""

    import myflopy

    broken = []
    for name in sorted(_export_names(tree)):
        try:
            getattr(myflopy, name)
        except Exception as exc:  # noqa: BLE001 - reporting, not handling
            broken.append(f"{name}: {type(exc).__name__}: {exc}")
    assert not broken, "exported names that fail to resolve:\n  " + "\n  ".join(broken)
