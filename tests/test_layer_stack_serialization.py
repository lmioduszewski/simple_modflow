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


# -- a tuple INSIDE a token's payload (D1 of the excavate scoping) ---------------------
#
# `_json_value` tags every tuple it meets as `{"$tuple": [...]}` -- inside a token's
# payload too -- but `_spec_value` handed the payload to `from_dict` raw, so the tag
# was never taken off. `LayerStack.add` takes "a sequence of shares" for `split=` and
# one name per share for `names=`, and a tuple is a sequence -- but it came back as a
# DICT: the saved stack died building with `could not convert string to float:
# '$tuple'`, and its digest no longer matched the live declaration's, so
# `runs.drift()` would have called it edited.


def _tuple_split_stack():
    """A stack whose only tuple is the unequal split -- nothing else could be tagged."""

    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=30.0, split=(0.3, 0.7), names=("sand_top", "sand_base"))
    stack.add("clay", thickness=10.0)
    return stack


def test_a_tuple_split_comes_back_a_tuple(triangulated_grid):
    """Through the funnel a package input takes: `_json_value` -> JSON -> `_spec_value`."""

    stack = _tuple_split_stack()
    written = _json_value(stack)
    assert written["$layer_stack"]["layers"][0]["split"] == {"$tuple": [0.3, 0.7]}, (
        "the tag is what makes this case worth testing; if the writer stops "
        "tagging, this test proves nothing"
    )

    revived = _spec_value(json.loads(json.dumps(written)))

    split, names = revived._layers[0].split, revived._layers[0].names
    assert type(split) is tuple and split == (0.3, 0.7)
    assert type(names) is tuple and names == ("sand_top", "sand_base")
    original = stack.build(triangulated_grid, attach=False)
    rebuilt = revived.build(triangulated_grid, attach=False)
    assert rebuilt.nlay == original.nlay == 3
    assert rebuilt.names == original.names
    assert np.array_equal(original.botm, rebuilt.botm)
    assert np.array_equal(original.idomain, rebuilt.idomain)


def test_a_saved_tuple_split_rebuilds_and_digests_the_same(tmp_path):
    """Through `Project.save`/`load`, where a lazy `mf.disv` carries the stack.

    The digest is the half a build test cannot see: the reloaded package digested a
    DICT where the live one digested a tuple, so the same declaration read as edited.
    """

    from myflopy.specs import package_digest
    from myflopy.workspace import Project

    project = Project(tmp_path / "p", name="p")
    project.add_package("disv/split", mf.disv(_tuple_split_stack()))
    live = package_digest(project.package("disv/split"))
    project.save()

    reloaded = Project.load(tmp_path / "p").package("disv/split")

    split = reloaded.inputs["layer_stack"]._layers[0].split
    assert type(split) is tuple and split == (0.3, 0.7)
    assert package_digest(reloaded) == live


def test_a_tuple_split_survives_yaml():
    """YAML renders the tag in flow style, then reads it back through the same funnel."""

    simulation = mf.simulation(mf.gwf("m", packages=[mf.disv(_tuple_split_stack())]), name="s")

    revived = mf.SimulationSpec.from_yaml(simulation.to_yaml())

    disv = next(p for p in revived.models[0].packages if p.name == "disv")
    split = disv.inputs["layer_stack"]._layers[0].split
    assert type(split) is tuple and split == (0.3, 0.7)


def test_a_tuple_grass_kwarg_comes_back_a_tuple():
    """The `$surface` token had the same gap: `grass_kwargs` is free-form.

    Asserted on the GRASS cache signature, which is where it shows: it is built from
    `str(v)`, so a revived `{'$tuple': [1, 2]}` read `stale` against a cache written
    for `(1, 2)` and would have re-run the interpolation.
    """

    surface = Surface.from_contours(
        "contours.gpkg", region_vector="domain.gpkg", out="ground.tif", extra=(1, 2)
    )

    revived = _spec_value(json.loads(json.dumps(_json_value(surface))))

    assert type(revived.grass_kwargs["extra"]) is tuple
    assert revived.to_dict() == surface.to_dict()
    assert revived._derived_raster().params == surface._derived_raster().params


def test_the_untag_leaves_the_rest_of_a_payload_alone():
    """Only the tag comes off: a nested surface stays the plain dict `from_dict` wants.

    Running the whole of `_spec_value` on a payload would be wrong -- it revives any
    dict whose `kind` names a source and promotes digit-string keys to ints -- so the
    untag is deliberately narrower than that.
    """

    from myflopy.specs import _untag_tuples

    payload = {
        "top": {"kind": "flat", "value": 1.0},
        "keys": {"0": "a"},
        "nested": [{"$tuple": [1, {"$tuple": [2, 3]}]}],
    }

    untagged = _untag_tuples(payload)

    assert untagged["top"] == {"kind": "flat", "value": 1.0}
    assert untagged["keys"] == {"0": "a"}
    assert untagged["nested"] == [(1, (2, 3))]
    assert type(untagged["nested"][0][1]) is tuple
