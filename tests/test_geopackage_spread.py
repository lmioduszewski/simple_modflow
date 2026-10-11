"""``mf.Spread``: splitting an extensive GeoPackage field across the cells it covers.

The behaviour being pinned is that a conductance belongs to the *feature*. Without
``Spread`` a feature's value is written to every cell it intersects, which is right
for an elevation and multiplies a conductance by the cell count -- on real data, by
2.02 on the very grid the values came from.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point, Polygon

import myflopy as mf
from myflopy.geopackage import CellSurfaceOffset, Spread, _extent


@pytest.fixture
def grid_and_context(tmp_path):
    """A 4x1 strip of unit square cells, and a context over it."""

    from myflopy.specs import ModelContext

    polys = gpd.GeoDataFrame(
        {"geometry": [Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)]) for i in range(4)]},
        crs="EPSG:2927",
    )

    class _Grid:
        gdf_vorPolys = polys
        ncpl = 4

        def get_grid_edge(self):
            return [0, 3]

    grid = _Grid()
    return grid, ModelContext(grid=grid, domain=np.ones((1, 4), dtype=int))


def _write(tmp_path, geometry, **columns):
    """Write a one-row GeoPackage and return its path."""

    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {**{k: [v] for k, v in columns.items()}, "geometry": [geometry]}, crs="EPSG:2927"
    ).to_file(path, layer="bc", driver="GPKG")
    return path


def _records(spec):
    """The period-0 records from a built PackageSpec."""

    return spec.options["stress_period_data"][0]


def test_extent_prefers_area_then_length():
    """A polygon measures by area, a line by length, a point by neither."""

    assert _extent(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])) == pytest.approx(4.0)
    assert _extent(LineString([(0, 0), (3, 0)])) == pytest.approx(3.0)
    assert _extent(Point(1, 1)) == 0.0


def test_without_spread_the_value_is_broadcast(tmp_path, grid_and_context):
    """The documented prior behaviour: every intersected cell gets the whole value."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1,
                       elevation="elevation", conductance="cond", boundnames=False)
    records = _records(spec)
    assert len(records) == 4
    assert sum(r[-1] for r in records) == pytest.approx(400.0)


def test_spread_conserves_the_total_across_cells(tmp_path, grid_and_context):
    """A line crossing four cells hands each its own share, summing to the original."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    records = _records(spec)
    assert sum(r[-1] for r in records) == pytest.approx(100.0)
    # 0.5 + 1 + 1 + 0.5 of a 3-unit line
    assert sorted(round(r[-1], 6) for r in records) == [
        pytest.approx(100 / 6), pytest.approx(100 / 6),
        pytest.approx(100 / 3), pytest.approx(100 / 3),
    ]


def test_spread_leaves_intensive_fields_alone(tmp_path, grid_and_context):
    """Only the wrapped field is split; an elevation stays whole on every cell."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    assert {r[1] for r in _records(spec)} == {10.0}


def test_a_point_keeps_its_whole_value(tmp_path, grid_and_context):
    """A point has no extent to divide, so Spread is a no-op rather than a zero."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, Point(1.5, 0.5), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    records = _records(spec)
    assert len(records) == 1
    assert records[0][-1] == pytest.approx(100.0)


def test_clip_loses_the_part_outside_the_grid(tmp_path, grid_and_context):
    """Half the line hangs past the mesh, so half the conductance is not applied."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(2.0, 0.5), (6.0, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond", mode="clip"), boundnames=False)
    assert sum(r[-1] for r in _records(spec)) == pytest.approx(50.0)


def test_retained_keeps_the_total_on_what_is_left(tmp_path, grid_and_context):
    """The same line under mode='retained' concentrates its whole total inside."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(2.0, 0.5), (6.0, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond", mode="retained"), boundnames=False)
    assert sum(r[-1] for r in _records(spec)) == pytest.approx(100.0)


def test_min_share_drops_a_corner_clip(tmp_path, grid_and_context):
    """A cell the feature barely grazes is dropped rather than given a sliver."""

    grid, ctx = grid_and_context
    # Two whole cells, then a 0.02 nick of cell 2: 0.99% of the feature, so
    # ``min_share`` drops it -- and 2% of the cell it nicks, twice the noise floor on
    # a line piece, so that does not. The old nick was 0.005, which is under the floor
    # and (measured: 0.004999999999999893 against a share guard of exactly 0.005, a
    # margin of 1e-16) never reached ``min_share`` at all.
    path = _write(tmp_path, LineString([(0.0, 0.5), (2.02, 0.5)]), layer=1,
                  elevation=10.0, cond=100.0)
    loose = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                        conductance=Spread("cond", min_share=0.0), boundnames=False)
    tight = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                        conductance=Spread("cond", min_share=0.01), boundnames=False)
    assert len(_records(loose)) == 3
    assert len(_records(tight)) == 2


def test_spread_rejects_an_unknown_mode():
    """A typo in mode is refused at construction, not silently treated as clip."""

    with pytest.raises(ValueError, match="clip.*retained"):
        Spread("cond", mode="area")


def test_spread_rejects_an_out_of_range_min_share():
    """min_share is a fraction; 1.0 would drop every cell."""

    with pytest.raises(ValueError, match="min_share"):
        Spread("cond", min_share=1.0)


def test_ghb_conductance_also_spreads(tmp_path, grid_and_context):
    """Spread is a value spec, not a DRN special case -- it rides any package's field."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, head=12.0, cond=90.0)
    spec = mf.ghb.gpkg(path, layer="bc", context=ctx, nper=1, head="head",
                       conductance=Spread("cond"), boundnames=False)
    records = _records(spec)
    assert sum(r[-1] for r in records) == pytest.approx(90.0)
    assert {r[1] for r in records} == {12.0}


def test_every_row_value_dataclass_survives_a_run_manifest():
    """Each dataclass in the RowValue union must convert to plain JSON.

    A run manifest is written with ``json.dumps`` and no ``default=``, so a spec
    that reaches it unconverted raises ``TypeError`` at ``prepare_run`` -- long
    after the package itself built cleanly. This fails the moment a new variant
    is added to the union without a ``_metadata_value`` branch.
    """

    import dataclasses
    import json
    import typing

    from myflopy.geopackage import CellSurfaceOffset, RowValue, _metadata_value

    dataclass_variants = [
        arg for arg in typing.get_args(RowValue)
        if isinstance(arg, type) and dataclasses.is_dataclass(arg)
    ]
    assert dataclass_variants, "RowValue should carry at least one dataclass variant"

    samples = {
        CellSurfaceOffset: CellSurfaceOffset("cell_top", offset=-1.0),
        Spread: Spread("conductance"),
    }
    missing = [cls.__name__ for cls in dataclass_variants if cls not in samples]
    assert not missing, f"add a sample for {missing} so this test still covers the union"

    for cls in dataclass_variants:
        converted = _metadata_value(samples[cls])
        assert isinstance(converted, dict), f"{cls.__name__} did not convert"
        assert converted["type"] == cls.__name__
        json.dumps(converted)  # the actual failure mode


def test_spread_metadata_records_its_settings():
    """The manifest keeps enough to reproduce the split, not just the field name.

    It did not, until ledger 274: `warn_below` was dropped, so whether a feature
    losing 40% of itself warned depended on a number no persisted artifact recorded,
    under a test whose own docstring claimed otherwise.
    """

    from myflopy.geopackage import _metadata_value

    assert _metadata_value(
        Spread("cond", mode="retained", min_share=0.05, warn_below=0.5)
    ) == {
        "type": "Spread",
        "field": "cond",
        "mode": "retained",
        "min_share": 0.05,
        "warn_below": 0.5,
    }


def test_the_two_converters_cover_the_same_fields():
    """`_metadata_value` and `to_dict` are separate on purpose, and must not diverge.

    They answer different questions -- provenance a human reads against a payload
    that round-trips -- so they are deliberately not one function (ledger 274). What
    ties them is this: neither may quietly stop recording a field. `to_dict` omits
    DEFAULTS, so it is compared against the dataclass's own field list.
    """

    import dataclasses

    from myflopy.geopackage import CellSurfaceOffset, _metadata_value

    for sample in (
        Spread("cond", mode="retained", min_share=0.05, warn_below=0.5),
        CellSurfaceOffset("cell_bottom", offset=2.0, minimum=1.0),
    ):
        declared = {f.name for f in dataclasses.fields(sample)}
        assert set(_metadata_value(sample)) - {"type"} == declared
        assert set(sample.to_dict()) == declared, "every field is non-default here"


def test_passthrough_cells_do_not_receive_boundaries(tmp_path):
    """idomain < 0 is vertical passthrough, not active -- MF6 refuses a boundary there.

    ``bool(-1)`` is ``True``, so a truthiness test lets these through and the run
    fails at read time with "Cell is outside active grid domain". A ``LayerStack``
    with ``pinch="passthrough"`` produces them in quantity.
    """

    from myflopy.specs import ModelContext

    polys = gpd.GeoDataFrame(
        {"geometry": [Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)]) for i in range(4)]},
        crs="EPSG:2927",
    )

    class _Grid:
        gdf_vorPolys = polys
        ncpl = 4

        def get_grid_edge(self):
            return [0, 3]

    # cell 1 active, cell 2 passthrough, cell 3 inactive, cell 4 active
    domain = np.array([[1, -1, 0, 1]])
    ctx = ModelContext(grid=_Grid(), domain=domain)

    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1,
                  elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    cells = sorted(r[0][1] for r in _records(spec))
    assert cells == [0, 3], "only the two truly active cells may carry a record"


def test_one_dimensional_domain_also_excludes_passthrough(tmp_path):
    """The single-layer branch of _active needs the same test, not truthiness."""

    from myflopy.specs import ModelContext

    polys = gpd.GeoDataFrame(
        {"geometry": [Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)]) for i in range(4)]},
        crs="EPSG:2927",
    )

    class _Grid:
        gdf_vorPolys = polys
        ncpl = 4

        def get_grid_edge(self):
            return [0, 3]

    ctx = ModelContext(grid=_Grid(), domain=np.array([1, -1, 0, 1]))
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1,
                  elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance="cond", boundnames=False)
    assert sorted(r[0][1] for r in _records(spec)) == [0, 3]


def test_a_spread_and_an_offset_round_trip_as_themselves():
    """The failure this guards is SILENT: without a revive arm they come back dicts.

    Asserts the TYPE, not that something non-None came back -- a plain dict passes
    every other assertion you would think to write (ledger 269, 274).
    """

    import json

    from myflopy.specs import _json_value, _spec_value

    for sample in (
        Spread("cond", mode="retained", min_share=0.05, warn_below=0.5),
        CellSurfaceOffset("cell_bottom", offset=2.0, minimum="floor_col"),
    ):
        revived = _spec_value(json.loads(json.dumps(_json_value(sample))))

        assert type(revived) is type(sample)
        assert revived == sample


@pytest.mark.parametrize(
    "offset", [-2.0, -2, "depth", "2.0", "1e5", "no", "null"]
)
def test_the_offset_union_keeps_its_side(offset):
    """`offset` is `float | str` -- a number, or the NAME of a GeoPackage column.

    A column literally called "2.0" is the case that would cross over and silently
    change every elevation it resolves. The YAML-hostile names ("no", "null") are
    here because PyYAML would otherwise read them as booleans and None.
    """

    import json

    import yaml

    sample = CellSurfaceOffset("model_top", offset=offset)
    for load, dump in ((json.loads, json.dumps), (yaml.safe_load, yaml.safe_dump)):
        revived = CellSurfaceOffset.from_dict(load(dump(sample.to_dict())))

        assert type(revived.offset) is type(offset)
        assert revived.offset == offset


def test_a_numpy_setting_does_not_look_like_an_edit():
    """`to_dict` normalizes through `_json_value` where `Surface.to_dict` does not.

    `type(np.float64(0.01)) is not float`, so an un-normalized numpy scalar fails the
    omit-defaults test, gets written, and returns a plain float that IS omitted --
    two payloads for one object, and a digest that moves across a save. The digest
    reads this payload since ledger 273, so the normalization is load-bearing.
    """

    import json

    from myflopy.specs import package_digest

    spec = dataclasses.replace(
        mf.npf(k=1.0), inputs={"c": Spread("cond", min_share=np.float64(0.01))}
    )
    reloaded = mf.PackageSpec.from_dict(json.loads(json.dumps(spec.to_dict())))

    assert package_digest(spec) == package_digest(reloaded)
    assert type(Spread("c", min_share=np.float64(0.05)).to_dict()["min_share"]) is float


def test_a_bare_payload_is_yaml_safe():
    """PyYAML dispatches on the EXACT type and refuses a numpy scalar json takes."""

    import yaml

    yaml.safe_dump(Spread("cond", min_share=np.float64(0.05)).to_dict())


def test_a_hand_edited_payload_is_refused_where_it_can_be_fixed():
    """Serializing these is what puts them in front of a text editor.

    Revival goes through `__init__`, so `__post_init__` re-runs and a bad payload
    raises here rather than dying inside a resolver several frames later.
    """

    with pytest.raises(ValueError, match="mode"):
        Spread.from_dict({"field": "cond", "mode": "nonsense"})
    with pytest.raises(ValueError, match="reference"):
        CellSurfaceOffset.from_dict({"reference": "sea_level"})
    with pytest.raises(ValueError, match="not usable"):
        Spread.from_dict({})
    with pytest.raises(ValueError, match="not usable"):
        Spread.from_dict({"field": "cond", "bogus": 1})
    with pytest.raises(ValueError, match="mapping"):
        Spread.from_dict("cond")


def test_a_package_holding_one_persists_as_json_not_a_pickle():
    """The user-visible win: a reviewable archive instead of a version-fragile cache.

    `_persist_package` falls back to pickling whatever `to_dict` refuses, and
    CLAUDE.md's rule is that a pickled entry is a cache, not an archive.
    """

    import tempfile

    from myflopy.workspace import Project

    with tempfile.TemporaryDirectory() as tmp:
        project = Project(Path(tmp) / "p", name="p")
        project.add_package(
            "npf/spread", dataclasses.replace(mf.npf(k=1.0), inputs={"c": Spread("cond")})
        )
        report = project.save()

        assert report.pickled == 0, "a spread no longer forces the pickle fallback"
        written = (Path(tmp) / "p" / "specs" / "packages" / "npf" / "spread.json").read_text()
        assert "$spread" in written
        assert type(Project.load(Path(tmp) / "p").package("npf/spread").inputs["c"]) is Spread


@pytest.mark.parametrize(
    "verb, kwargs",
    [
        ("mf.rch.array", {"recharge": Spread("c"), "irch": 0}),
        ("mf.evt.array", {"surface": 10.0, "rate": Spread("c"), "depth": 1.0, "ievt": 0}),
        (
            "mf.evt.array",
            {"surface": CellSurfaceOffset("model_top"), "rate": 1.0, "depth": 1.0, "ievt": 0},
        ),
    ],
)
def test_an_array_form_refuses_a_per_feature_value(verb, kwargs):
    """Neither can act without features, and a parameter that cannot act is a defect.

    It is also what makes ledger 274's "no stored digest moves" structural rather
    than merely unobserved: `prepare_run` writes the manifest BEFORE it builds
    (and `build=False` is supported), so a run could otherwise record a digest for
    a package holding one of these, which tokenizing would then move.
    """

    helper = mf.rch.array if verb == "mf.rch.array" else mf.evt.array

    with pytest.raises(TypeError, match="per-FEATURE resolution"):
        helper(**kwargs)


def test_the_array_form_still_takes_an_ordinary_value():
    """The control: the refusal must not be bought with a false positive."""

    assert mf.rch.array(recharge=0.001, irch=0) is not None


def _one_feature_file(tmp_path):
    """A one-row GeoPackage, so the refusal under test is reached rather than the file check."""

    import geopandas as gpd
    from shapely.geometry import box

    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["f"], "layer": [1], "x": [1.0]},
        geometry=[box(0.0, 0.0, 1.0, 1.0)],
        crs="EPSG:2927",
    ).to_file(path)
    return path


@pytest.mark.parametrize(
    ("helper_name", "field"),
    [("chd", "head"), ("ghb", "head"), ("drn", "elevation"),
     ("riv", "stage"), ("riv", "rbot"), ("evt", "surface")],
)
def test_spread_is_refused_on_an_elevation_field(tmp_path, helper_name, field):
    """`Spread` divides a feature's total across the cells it covers. An elevation has no
    total to divide, so each cell would get a FRACTION of the elevation declared.

    Measured before this: `elevation=mf.Spread("elev")` on a declared 75.0 wrote 2.06 to
    11.87 -- the elevation times each cell's area share -- with all 11 records below their
    own cell bottom, a clean build and no warning (ledger 289). Ledger 274 refused exactly
    this on the ARRAY form and never reached the list form.

    Asserted on the LAZY declaration (no `context=`), because that path never reaches
    `_bc_from_features` -- a single check site would have covered only half the API.
    """

    import myflopy as mf
    from myflopy.geopackage import Spread

    helper = getattr(mf, helper_name)
    with pytest.raises(TypeError, match=f"{field}=mf.Spread"):
        helper.gpkg(_one_feature_file(tmp_path), nper=1, **{field: Spread("x")})


def test_spread_is_still_allowed_where_there_IS_a_total_to_divide(tmp_path):
    """The refusing direction. A conductance, and a recharge column that really holds a
    volume, ARE spreadable -- that is what `Spread`'s own docstring teaches -- so the
    refusal must name elevations only."""

    import myflopy as mf
    from myflopy.geopackage import Spread

    path = _one_feature_file(tmp_path)
    for helper_name, field in (("drn", "conductance"), ("rch", "recharge"), ("wel", "rate")):
        spec = getattr(mf, helper_name).gpkg(path, nper=1, **{field: Spread("x")})
        assert spec.inputs["stress_period_data"].values[field].field == "x"


def test_a_misspelled_spread_column_is_named_before_any_record(tmp_path, grid_and_context):
    """A column named only inside `Spread(...)` is in the up-front column check.

    It was left out, so a typo there died mid-build as a bare `KeyError('typo')` naming
    neither the file nor the feature, where the same typo unwrapped got the named error.
    Both spellings now say the same thing.
    """

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    for conductance in ("typo", Spread("typo")):
        with pytest.raises(ValueError, match=r"GeoPackage 'bc.gpkg' is missing fields: typo"):
            mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                        conductance=conductance, boundnames=False)


def _evt_file(tmp_path):
    path = tmp_path / "evt.gpkg"
    gpd.GeoDataFrame(
        {"name": ["et"], "layer": [1], "surface": [10.0], "rate": [1.0], "depth": [2.0]},
        geometry=[LineString([(2.0, 0.5), (6.0, 0.5)])], crs="EPSG:2927",
    ).to_file(path, layer="bc", driver="GPKG")
    return path


@pytest.mark.parametrize("lazy", [False, True])
def test_two_spreads_with_different_modes_are_refused_naming_both(tmp_path, grid_and_context, lazy):
    """A feature has ONE set of cell shares, so one mode; the FIRST was used for both.

    Measured before the fix: `rate=Spread(mode="clip")` beside
    `depth=Spread(mode="retained")` on a line half off the grid wrote the depth CLIPPED
    -- half of what `retained` keeps -- with nothing said (ledger 396). Refused at
    the call, eager and lazy.
    """

    _, ctx = grid_and_context
    kwargs = {} if lazy else {"context": ctx}
    with pytest.raises(ValueError, match=r"rate=Spread\('rate', mode='clip'\), "
                                         r"depth=Spread\('depth', mode='retained'\)"):
        mf.evt.gpkg(_evt_file(tmp_path), layer="bc", nper=1, **kwargs,
                    rate=Spread("rate"), depth=Spread("depth", mode="retained"))


def test_two_spreads_with_one_mode_still_build(tmp_path, grid_and_context):
    """The control: agreeing modes are what one set of shares can honour."""

    _, ctx = grid_and_context
    spec = mf.evt.gpkg(_evt_file(tmp_path), layer="bc", context=ctx, nper=1,
                       rate=Spread("rate", mode="retained"), depth=Spread("depth", mode="retained"))
    records = _records(spec)
    assert sum(r[2] for r in records) == pytest.approx(1.0)
    assert sum(r[3] for r in records) == pytest.approx(2.0)


@pytest.mark.parametrize("route", ["source", "eager", "lazy"])
def test_a_spread_refusal_names_the_call_that_reached_it(tmp_path, grid_and_context, route):
    """Both declaration-time Spread refusals name the route actually taken.

    They said `mf.<pkg>.gpkg` whichever call reached them, so a caller of
    `GeoPackageSource(...).drn(...)` was told about a helper they never called.
    """

    from myflopy.geopackage import GeoPackageSource

    _, ctx = grid_and_context
    path = _evt_file(tmp_path)
    named = {"source": "GeoPackageSource", "eager": "mf", "lazy": "mf"}[route]
    calls = {
        "source": lambda pkg, **kw: getattr(
            GeoPackageSource(path, ctx, nper=1, layer="bc"), pkg)(**kw),
        "eager": lambda pkg, **kw: getattr(mf, pkg).gpkg(path, layer="bc", context=ctx,
                                                         nper=1, **kw),
        "lazy": lambda pkg, **kw: getattr(mf, pkg).gpkg(path, layer="bc", nper=1, **kw),
    }[route]
    spelled = "GeoPackageSource.{}" if named == "GeoPackageSource" else "mf.{}.gpkg"
    with pytest.raises(TypeError, match=rf"^{spelled.format('drn')}: elevation=mf\.Spread"):
        calls("drn", elevation=Spread("surface"), conductance="rate")
    with pytest.raises(ValueError, match=rf"^{spelled.format('evt')}: rate=Spread"):
        calls("evt", rate=Spread("rate"), depth=Spread("depth", mode="retained"))
