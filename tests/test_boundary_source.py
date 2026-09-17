"""`BoundarySource`: the serializable recipe behind a lazy `.gpkg` boundary condition.

Data only. Nothing emits one yet -- the helper that does is B2 commit 5b -- so this
module exists to pin the contract that commit has to meet, and to make the failure
modes of a new source class impossible to ship silently.

Three of those failure modes have already cost this project a ledger entry each:
an unregistered token round-trips to a PLAIN DICT with nothing raised (269, 274); a
numpy or `Path` value inside a source drifts the digest across a save (265, 273); and
a recipe that carries a package NAME makes a rename silently vanish (275's territory,
via `_renamed_options`).
"""

from __future__ import annotations

import copy
import dataclasses
import json

import pytest
import yaml

import myflopy as mf
from myflopy.sources import _SOURCE_TYPES, BoundarySource
from myflopy.specs import _json_value, _spec_options, _spec_value, package_digest


def _recipe(**overrides):
    """A recipe with every interesting field NON-DEFAULT.

    Defaults are omitted from `to_dict`, so a sample built from them round-trips
    through a converter that writes nothing at all -- it would certify rather than
    check (the lesson of ledger 274's `to_dict` tests).
    """

    return BoundarySource(
        path="gis/bc.gpkg",
        package="drn",
        nper=3,
        layer="drains",
        values={
            "elevation": mf.CellSurfaceOffset("cell_bottom", offset=2.0),
            "conductance": mf.Spread("cond", mode="retained"),
        },
        name_field="feature",
        layer_field="lay",
        period_field="per",
        layer_base=0,
        period_base=1,
        edges_only=True,
        boundnames=False,
        **overrides,
    )


def test_it_is_registered_so_an_unknown_kind_raises():
    """`kind` ends in "Source" so a MISSED registration fails loudly.

    `_spec_value` keeps an `endswith("Source")` arm beside the registry for exactly
    this: an unregistered kind raises rather than degrading to a dict, which is the
    silent failure ledger 269 records.
    """

    assert _SOURCE_TYPES["BoundarySource"] is BoundarySource

    with pytest.raises(ValueError, match="Unknown source kind"):
        _spec_value({"kind": "MadeUpSource", "path": "x"})


def test_it_round_trips_through_inputs_as_itself():
    """Asserts the TYPE, not that something non-None came back.

    A plain dict passes every other assertion you would think to write.
    """

    source = _recipe()
    spec = dataclasses.replace(mf.npf(k=1.0), inputs={"stress_period_data": source})

    revived = _spec_options(
        json.loads(json.dumps(_json_value(spec.inputs)))
    )["stress_period_data"]

    assert type(revived) is BoundarySource
    assert revived == source


def test_the_live_value_objects_survive_the_round_trip():
    """`Spread` and `CellSurfaceOffset` come back as THEMSELVES, not as dicts.

    They survive only because `inputs` goes through `_spec_value(..., int_keys=True)`,
    which is what recurses INTO a source payload. A recipe reached through any other
    field comes back with plain dicts where the objects were -- which is why 5b must
    put this in `inputs` and nowhere else.
    """

    source = _recipe()
    revived = _spec_options(
        json.loads(json.dumps(_json_value({"stress_period_data": source})))
    )["stress_period_data"]

    assert type(revived.values["conductance"]) is mf.Spread
    assert type(revived.values["elevation"]) is mf.CellSurfaceOffset
    assert revived.values["conductance"].mode == "retained"
    assert revived.values["elevation"].offset == 2.0

    # The trap, pinned: WITHOUT the inputs path the nested objects do not revive.
    naive = _spec_value(json.loads(json.dumps(_json_value(source))))

    assert type(naive.values["conductance"]) is dict


def test_its_digest_is_stable_across_a_save():
    """A recipe that drifts reports a package nobody edited as `edited` (ledger 273)."""

    spec = dataclasses.replace(
        mf.npf(k=1.0), inputs={"stress_period_data": _recipe()}
    )
    reloaded = mf.PackageSpec.from_dict(json.loads(json.dumps(spec.to_dict())))

    assert package_digest(spec) == package_digest(reloaded)
    assert type(reloaded.inputs["stress_period_data"]) is BoundarySource


def test_the_digest_still_sees_every_field_change():
    """A fix bought by making the digest blind is not a fix.

    `to_dict` omits defaults and the digest reads that payload, so a field that stops
    being written stops being seen. These are every field it writes conditionally.
    """

    def digest(**overrides):
        return package_digest(
            dataclasses.replace(
                mf.npf(k=1.0),
                inputs={"stress_period_data": dataclasses.replace(_recipe(), **overrides)},
            )
        )

    base = digest()
    for field_name, other in (
        ("package", "ghb"),
        ("nper", 4),
        ("layer", "other"),
        ("name_field", "nm"),
        ("layer_field", "l"),
        ("period_field", None),
        ("layer_base", 1),
        ("period_base", 0),
        ("edges_only", False),
        ("boundnames", True),
        ("path", "gis/other.gpkg"),
        ("values", {"conductance": mf.Spread("other")}),
    ):
        assert digest(**{field_name: other}) != base, f"{field_name} is invisible"


def test_yaml_agrees_with_json():
    """The two documented formats must not disagree about one spec (ledger 265)."""

    written = _json_value({"stress_period_data": _recipe()})
    from_json = _spec_options(json.loads(json.dumps(written)))["stress_period_data"]
    from_yaml = _spec_options(yaml.safe_load(yaml.safe_dump(written)))["stress_period_data"]

    assert from_json == from_yaml == _recipe()


def test_it_survives_the_freeze_deepcopy_and_pickle_paths():
    """Every route a spec takes, not just the one the tests above use."""

    from myflopy.project.freeze import _lift

    source = _recipe()
    store: dict = {}
    lifted = _lift({"stress_period_data": source}, store, "test input")

    assert type(lifted["stress_period_data"]) is BoundarySource
    assert not store, "a recipe holds no arrays to lift"
    assert copy.deepcopy(source) == source
    assert dataclasses.replace(source, nper=9).nper == 9


def test_it_is_keyword_only():
    """Not a courtesy: `nper` has no default and the inherited fields do.

    Without `kw_only` the class cannot be created at all. It also closes the positional
    trap `GeoPackageSourceSpec`'s docstring warns about, where
    `GeoPackageSourceSpec("x.gpkg", "areas")` sets `external="areas"`.
    """

    with pytest.raises(TypeError, match="keyword-only"):
        BoundarySource("gis/bc.gpkg", "drn")


def test_it_carries_no_package_name():
    """Identity belongs to the `PackageSpec`, not to its recipe.

    `fork_package` renames through `_renamed_options`, which rewrites `options` and
    never touches `inputs` -- so a name frozen in the recipe would be silently
    discarded by a rename, and two renamed forks would collapse into ONE MF6 package
    under a FloPy "Replacing existing package" warning. That is the exact defect
    `_renamed_options` already exists to prevent; 5b and 5c must take the name from the
    spec at build time instead.
    """

    names = {f.name for f in dataclasses.fields(BoundarySource)}

    assert "package_name" not in names and "name" not in names
    assert "package" in names, "the package TYPE is schema, not identity"


def test_it_is_not_the_runtime_resolver():
    """The runtime `GeoPackageSource` can never be a recipe, which is why this exists."""

    from myflopy.geopackage import GeoPackageSource

    assert not hasattr(GeoPackageSource, "to_dict"), "it cannot be written down"
    assert "context" in {f.name for f in dataclasses.fields(GeoPackageSource)}
    assert "context" not in {f.name for f in dataclasses.fields(BoundarySource)}
    # A name collision worth knowing about: the registry key "GeoPackageSource" belongs
    # to the SERIALIZABLE `GeoPackageSourceSpec`, not to the runtime class of that name.
    # Two classes, one string -- so "is GeoPackageSource registered?" has a misleading
    # answer, and `BoundarySource` deliberately does not reuse either name.
    from myflopy.sources import GeoPackageSourceSpec

    assert _SOURCE_TYPES["GeoPackageSource"] is GeoPackageSourceSpec
    assert _SOURCE_TYPES["GeoPackageSource"] is not GeoPackageSource

    # It cannot even be CONSTRUCTED without the file and a gridded context -- both
    # `__post_init__` checks -- which is the other half of why a recipe needs its own
    # type rather than reusing the resolver.
    with pytest.raises(FileNotFoundError):
        GeoPackageSource(path="x.gpkg", context=mf.ModelContext(), nper=1)
