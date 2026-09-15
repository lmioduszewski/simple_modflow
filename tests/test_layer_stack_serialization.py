"""A layer stack you can write down.

A `LayerStack` is a DECLARATION -- a top surface and a list of named units -- and
until now it could not be serialized at all, so a model carrying one could not ride
in a spec. The consequence was larger than it sounds: `mf.disv(result)` holds the
RESOLVED `top`/`botm`/`idomain` arrays, `to_dict` refuses an ndarray (ledger 193),
and so `Project.save()` persisted **zero** simulations for every layered model --
the canonical way to build one here -- naming the skip in `SaveReport.skipped` and
leaving `load()` to return an empty project.

This commit makes the declaration writable. It does NOT yet make `mf.disv` carry one
instead of arrays; that is the next commit, and it is what actually closes the save
(ledger 269).

**A stack does not determine its own arrays.** `reconcile`, `min_sep`, `trigger_sep`
and `method` live on the `build()`/`to_disv()` call, so a round-tripped stack is a
recipe that still needs them -- which is why they are forwarded in the next commit.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

import myflopy as mf
from myflopy.layers import LayerStack
from myflopy.specs import _json_value, _spec_value
from myflopy.surfaces import Surface


def _declared_stack():
    """A stack exercising every per-unit control the payload has to carry."""

    stack = mf.LayerStack(top=mf.Surface.flat(100.0), length_units="feet", time_units="days")
    stack.add("sand", thickness=20.0, pinch="inactive", min_thickness=2.0)
    stack.add("clay", thickness=10.0, split=3)
    stack.add("gravel", thickness=5.0, split=2, names=["gravel_upper", "gravel_lower"])
    return stack


def test_a_surface_round_trips_through_its_own_payload():
    for surface in (
        Surface.flat(10.0),
        Surface.raster("dem.tif"),
        Surface.flat(10.0) - 5.0,
    ):
        assert Surface.from_dict(surface.to_dict()) == surface


def test_an_algebraic_surface_keeps_its_operands():
    """`operands` is the one recursive field, so it is the one that can silently flatten."""

    derived = Surface.flat(100.0) - 12.5

    revived = Surface.from_dict(derived.to_dict())

    assert revived == derived
    assert len(revived.operands) == 1 and revived.operands[0].value == 100.0


@pytest.mark.parametrize("maker, kind", [
    (lambda: Surface.from_array(np.zeros(4)), "array"),
    (lambda: Surface.from_points(np.zeros(3), np.zeros(3), np.zeros(3)), "points"),
])
def test_a_surface_holding_arrays_is_refused_by_name(maker, kind):
    """Inlining per-cell arrays into a spec is what ledger 193 refused.

    Refusing by name beats degrading: a stack that quietly dropped its array surface
    would round-trip into a DIFFERENT model, which is the failure mode this whole
    phase exists to remove.
    """

    with pytest.raises(ValueError, match=rf"'{kind}' surface holds per-cell arrays"):
        maker().to_dict()


def test_a_declared_stack_round_trips_with_every_per_unit_control():
    stack = _declared_stack()

    revived = LayerStack.from_dict(stack.to_dict())

    assert [layer.name for layer in revived._layers] == ["sand", "clay", "gravel"]
    assert revived._layers[0].pinch == "inactive"
    assert revived._layers[0].min_thickness == 2.0
    assert revived._layers[1].split == 3
    assert revived._layers[2].names == ["gravel_upper", "gravel_lower"]
    assert revived.length_units == "feet" and revived.time_units == "days"
    assert revived._top == stack._top


def test_the_payload_carries_the_declaration_not_the_grid(triangulated_grid):
    """`vor` and the cached draft grid are deliberately absent.

    A grid cannot be inlined (ledger 193), and a stack is declarable without one --
    that is what lets the project-first order work. Tested on a GRID-BOUND stack,
    because a gridless one would pass whether or not `vor` is excluded.
    """

    bound = mf.LayerStack(triangulated_grid, top=mf.Surface.flat(10.0))
    bound.add("sand", thickness=4.0)
    assert bound.vor is not None, "a gridless stack proves nothing here"

    payload = bound.to_dict()

    assert set(payload) == {"top", "layers", "length_units", "time_units"}
    assert LayerStack.from_dict(payload).vor is None
    assert LayerStack.from_dict(payload, vor=triangulated_grid).vor is triangulated_grid


def test_a_stack_survives_the_spec_serializer_and_json():
    """The `$layer_stack` token, through the same funnel a package option takes."""

    stack = _declared_stack()

    written = _json_value(stack)
    revived = _spec_value(json.loads(json.dumps(written)))

    assert isinstance(revived, LayerStack)
    assert [layer.name for layer in revived._layers] == ["sand", "clay", "gravel"]


def test_the_token_does_not_collide_with_a_source_kind():
    """`Surface` HAS a `kind` field, and the source dispatcher keys on `kind`.

    Tagging a surface as `{"kind": "flat"}` made it come back a plain dict -- no
    error, just a different object. The `$`-token is what keeps the two apart.
    """

    written = _json_value(Surface.raster("dem.tif"))

    assert set(written) == {"$surface"}
    assert isinstance(_spec_value(written), Surface)


def test_a_round_tripped_stack_builds_the_same_arrays(triangulated_grid):
    """The point of a recipe: it has to reproduce what it described."""

    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)

    revived = LayerStack.from_dict(stack.to_dict())

    original = stack.build(triangulated_grid, attach=False)
    rebuilt = revived.build(triangulated_grid, attach=False)
    assert np.array_equal(original.top, rebuilt.top)
    assert np.array_equal(original.botm, rebuilt.botm)
    assert np.array_equal(original.idomain, rebuilt.idomain)
