"""Every ``_SPEC_TOKEN`` class must have a revive arm, or it round-trips to a dict.

``_json_value`` writes a token-bearing object as ``{token: payload}`` by duck-typing
on the class attribute, so declaring ``_SPEC_TOKEN`` is enough to get WRITTEN. Being
READ is a separate, hand-maintained if-chain in ``specs._spec_value``. Forget the arm
and the object comes back a **plain dict, with nothing raised** -- ledger 269 recorded
that trap for ``Surface``, and ledger 274 measured it again for ``Spread``:

    written: {"c": {"$spread": {"field": "cond"}}}
    read back: type=builtins.dict, equal to the original: False, raised: nothing

A half-registered package builds, saves and reloads under normal termination with
every such object in it turned inert. This module makes that impossible to ship.

It reads the tree by AST rather than importing it -- 0.4 s against 4.6 s for a
``pkgutil`` sweep, which is material against a 32 s inner loop. Two constraints follow
from that and are asserted rather than assumed: a ``_SPEC_TOKEN`` must be a string
LITERAL in a class body, and an arm must be spelled ``"$token" in value``. The
``$tuple`` form (``set(value) == {"$tuple"}``) is invisible to the sweep, which is why
it is exempt by name below rather than silently absent.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src" / "myflopy"

#: Tokens handled by `_spec_value` in a shape this sweep cannot see, or with no class.
_ARMS_WITHOUT_A_CLASS = {"$tuple"}


def _declared_tokens() -> dict[str, str]:
    """Every ``_SPEC_TOKEN = "..."`` in a class body under ``src/myflopy``."""

    found: dict[str, str] = {}
    for path in sorted(SRC.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for stmt in node.body:
                targets = (
                    stmt.targets
                    if isinstance(stmt, ast.Assign)
                    else [stmt.target]
                    if isinstance(stmt, ast.AnnAssign)
                    else []
                )
                for target in targets:
                    if not (isinstance(target, ast.Name) and target.id == "_SPEC_TOKEN"):
                        continue
                    assert isinstance(stmt.value, ast.Constant) and isinstance(
                        stmt.value.value, str
                    ), (
                        f"{node.name}._SPEC_TOKEN must be a string literal in the "
                        "class body so this test can read it without importing."
                    )
                    found[stmt.value.value] = f"{path.name}:{node.name}"
    return found


def _revive_arms() -> set[str]:
    """Every ``"$..." in value`` test inside ``myflopy.specs._spec_value``."""

    module = ast.parse((SRC / "specs.py").read_text())
    function = next(
        node
        for node in ast.walk(module)
        if isinstance(node, ast.FunctionDef) and node.name == "_spec_value"
    )
    return {
        node.left.value
        for node in ast.walk(function)
        if isinstance(node, ast.Compare)
        and isinstance(node.ops[0], ast.In)
        and isinstance(node.left, ast.Constant)
        and isinstance(node.left.value, str)
        and node.left.value.startswith("$")
    }


def test_every_spec_token_has_a_revive_arm():
    """Declaring a token gets you WRITTEN. Being read is a separate edit."""

    declared = _declared_tokens()
    missing = sorted(set(declared) - _revive_arms())

    assert not missing, (
        "declared with no revive arm: "
        + ", ".join(f"{token} ({declared[token]})" for token in missing)
        + " -- `_json_value` will WRITE it and `_spec_value` will hand it back as a "
        "PLAIN DICT with nothing raised. Add an arm in `specs._spec_value`."
    )


def test_every_revive_arm_has_a_class():
    """The other direction: a deleted class must not leave a dead arm behind."""

    stray = sorted(_revive_arms() - set(_declared_tokens()) - _ARMS_WITHOUT_A_CLASS)

    assert not stray, (
        f"`_spec_value` revives {stray}, which no class declares. Delete the arm, or "
        "add the token to _ARMS_WITHOUT_A_CLASS if it is handled without one."
    )


def test_the_sweep_can_actually_see_the_tokens():
    """A guard on the guard: an AST sweep that finds nothing passes both tests above."""

    declared = _declared_tokens()

    assert len(declared) >= 4, f"the sweep found only {declared} -- it has gone blind"
    assert {"$surface", "$layer_stack", "$spread", "$cell_surface_offset"} <= set(declared)


@pytest.mark.parametrize(
    "token, build",
    [
        ("$spread", lambda: __import__("myflopy").Spread(
            "cond", mode="retained", min_share=0.05, warn_below=0.5)),
        ("$cell_surface_offset", lambda: __import__("myflopy").CellSurfaceOffset(
            "cell_bottom", offset=2.0, minimum=1.0)),
        ("$surface", lambda: __import__("myflopy").Surface.flat(10.0)),
    ],
)
def test_a_token_object_survives_a_round_trip_as_itself(token, build):
    """Every field NON-DEFAULT, because `to_dict` omits defaults.

    A sample built from defaults round-trips through a `to_dict` that writes nothing
    at all, so it would certify a converter that drops every field.
    """

    from myflopy.specs import _json_value, _spec_value

    sample = build()
    written = _json_value(sample)
    assert set(written) == {token}, written

    revived = _spec_value(written)

    assert type(revived) is type(sample), f"{token} revived as {type(revived).__name__}"
    assert revived == sample
