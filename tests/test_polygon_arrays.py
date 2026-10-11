"""A polygon map fills a per-cell array parameter, lazily and saveably.

``GeoPackageSource.k_array`` has painted K from polygons for a long time and still
does. What it cannot be is LAZY or SAVED: it needs the grid at DECLARATION, so its
model cannot be moved to another mesh, and it hands back an ``ndarray``, which
``to_dict`` refuses -- so ``Project.save()`` persisted ZERO simulations for such a
model. Both halves are asserted here against the recipe form (ledger 299).

The rules are the other half. Three of them exist in this repo and they DISAGREE; the
tests that matter are the ones proving they are three rules rather than one wearing
three names, and that the default invents no value the map does not contain.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.geopackage import array_from_polygons
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.modflow.mf6.grid.geometry import polygon_cell_values
from myflopy.sources import POLYGON_RULES, Polygons, source_from_dict
from myflopy.specs import ModelContext, package_digest


def _grid(n):
    return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1))


@pytest.fixture(scope="module")
def meshes():
    """Two DIFFERENT meshes -- the pair is what a recipe has to survive."""

    return _grid(6), _grid(10)


def _write(path, vor, *, values=((0.0, 0.5, 1.0), (0.5, 1.0, 100.0)), layer=None):
    """A geology map: one row per (x-fraction span, k), optionally layer-tagged."""

    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    rows, geoms = [], []
    for index, (start, stop, k) in enumerate(values):
        rows.append({"name": f"unit{index}", "k": float(k), "layer": 1})
        geoms.append(box(xmin + (xmax - xmin) * start, ymin, xmin + (xmax - xmin) * stop, ymax))
    frame = gpd.GeoDataFrame(rows, geometry=geoms, crs=vor.crs)
    if layer is not None:
        frame["layer"] = layer
    frame.to_file(path, layer="geology", driver="GPKG")
    return path


@pytest.fixture(scope="module")
def geology(tmp_path_factory, meshes):
    """Written against the LARGER mesh, so it covers both.

    A map covering only the first mesh is not a fixture bug to paper over with
    `default=`: a move onto a footprint the geology does not reach IS uncovered, and
    the refusal that catches it is tested on its own below.
    """

    _, large = meshes
    return _write(tmp_path_factory.mktemp("poly") / "geology.gpkg", large)


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    return stack


def _model(grid, package, tmp_path, name):
    """Build a one-array-package model on `grid` and return the built simulation."""

    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid),
               packages=[mf.disv(_stack()), package, mf.ic(strt=5.0)]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    return sim.build(tmp_path / name)


# --------------------------------------------------------------------------- the point


def test_one_declaration_fills_k_on_two_different_meshes(meshes, geology, tmp_path):
    """The reason this exists: the recipe resolves against whatever grid is supplied."""

    small, large = meshes
    recipe = mf.Polygons(path=geology, column="k", layer="geology")
    written = []
    for index, grid in enumerate((small, large)):
        built = _model(grid, mf.npf(k=recipe), tmp_path, f"m{index}")
        written.append(built.simulation.get_model("m").get_package("npf").k.get_data())

    assert written[0].shape == (1, small.ncpl)
    assert written[1].shape == (1, large.ncpl)
    for array in written:
        assert set(np.unique(array)) == {1.0, 100.0}


def test_a_polygon_backed_array_moves_rebuilt_where_a_plain_one_cannot(meshes, geology, tmp_path):
    """`unverifiable` -> `rebuilt`, which is the whole justification (ledgers 295, 299)."""

    small, large = meshes
    plain = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small),
               packages=[mf.disv(_stack()), mf.npf(k=np.full(small.ncpl, 10.0))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    assert _verdict(plain.with_grid(large), "npf") == "unverifiable"

    lazy = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small),
               packages=[mf.disv(_stack()),
                         mf.npf(k=mf.Polygons(path=geology, column="k", layer="geology"))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    moved = lazy.with_grid(large)
    assert _verdict(moved, "npf") == "rebuilt"
    built = moved.build(tmp_path / "moved")
    array = built.simulation.get_model("m").get_package("npf").k.get_data()
    assert array.shape == (1, large.ncpl)
    assert set(np.unique(array)) == {1.0, 100.0}


def _verdict(sim, package):
    return next(row["how"] for row in sim.move_report() if row["package"] == package)


def test_a_project_carrying_a_polygon_map_actually_saves(meshes, geology, tmp_path):
    """The `k_array` failure, in the one place it is visible: `simulations=0`.

    An ndarray cannot be serialized (ledger 193), so a model painted by ``k_array``
    persisted NOTHING and ``load()`` returned an empty project -- silently, since
    ``save()`` reports rather than raises.
    """

    small, _ = meshes
    resolver = mf.GeoPackageSource(geology, ModelContext(grid=small), 1, layer="geology")
    eager = mf.npf(k=resolver.k_array(value="k", nlay=1, defaults=5.0))
    lazy = mf.npf(k=mf.Polygons(path=geology, column="k", layer="geology"))

    for index, (label, package, expected) in enumerate(
        (("eager", eager, 0), ("lazy", lazy, 1))
    ):
        project = mf.Project(tmp_path / f"proj{index}", name=label)
        project.add_grid("mesh", small)
        project.add_simulation(
            mf.simulation(
                mf.gwf("m", grid=mf.grid_ref("mesh"),
                       packages=[mf.disv(_stack()), package, mf.ic(strt=5.0)]),
                name="s", packages=[mf.tdis(nper=1)],
            )
        )
        report = project.save()
        assert report.simulations == expected, f"{label}: {report}"
        assert len(mf.Project.load(tmp_path / f"proj{index}").simulations) == expected


# ------------------------------------------------------------------------------ rules


def test_the_rules_are_three_rules_and_not_one(meshes, tmp_path):
    """`mean` invents values the map does not contain; the other two never do.

    This is the reason the rule is NAMED on the recipe rather than inherited from
    whichever helper happened to resolve it first.
    """

    small, _ = meshes
    path = _write(tmp_path / "contact.gpkg", small, values=((0.0, 0.47, 1.0), (0.47, 1.0, 100.0)))
    zones = gpd.read_file(path, layer="geology")
    cells = small.gdf_vorPolys
    painted = {
        rule: polygon_cell_values(cells, zones, zones["k"].to_numpy(), rule=rule)[0]
        for rule in POLYGON_RULES
    }

    declared = {1.0, 100.0}
    assert set(np.unique(painted["majority"])) <= declared
    assert set(np.unique(painted["centroid"])) <= declared
    invented = set(np.unique(painted["mean"])) - declared
    assert invented, "an area-weighted mean across a contact must blend, or it is not a mean"


def test_majority_and_centroid_are_not_the_same_rule(meshes, tmp_path):
    """They agree almost everywhere -- which is exactly why this has to be asserted.

    Over 19 contact positions on a 64-cell grid they differed on at most ONE cell, so a
    single fixture proves nothing and a refactor collapsing one branch into the other
    would pass every other test here. The disagreeing cell is the one whose CENTRE and
    whose BULK fall on opposite sides of the contact.
    """

    small, _ = meshes
    cells = small.gdf_vorPolys
    for index, split in enumerate(np.linspace(0.05, 0.95, 19)):
        path = _write(
            tmp_path / f"c{index}.gpkg", small,
            values=((0.0, float(split), 1.0), (float(split), 1.0, 100.0)),
        )
        zones = gpd.read_file(path, layer="geology")
        values = zones["k"].to_numpy()
        majority, _ = polygon_cell_values(cells, zones, values, rule="majority")
        centroid, _ = polygon_cell_values(cells, zones, values, rule="centroid")
        if (majority != centroid).any():
            return
    pytest.fail(
        "majority and centroid agreed on every cell of 19 contact positions -- one "
        "branch is answering for both"
    )


def test_the_rule_names_agree_across_the_two_modules():
    """``POLYGON_RULES`` (layer 0) and the resolver's own dispatch (layer 1) must match.

    The resolver cannot import the tuple -- it sits a layer ABOVE ``sources`` -- so the
    three names are stated twice on purpose, and this is what keeps the restatement
    from drifting into a rule that validates and then raises.
    """

    cells = _grid(3).gdf_vorPolys
    empty = gpd.GeoDataFrame({"k": []}, geometry=[], crs=cells.crs)
    for rule in POLYGON_RULES:
        polygon_cell_values(cells, empty, np.array([]), rule=rule)
    with pytest.raises(ValueError, match="is not a rule"):
        polygon_cell_values(cells, empty, np.array([1.0]), rule="nearest")


def test_an_unknown_rule_is_refused_at_the_declaration():
    with pytest.raises(ValueError, match="is not a rule"):
        mf.Polygons(path="x.gpkg", column="k", rule="closest")


def test_a_tie_goes_to_the_lowest_row(meshes, tmp_path):
    """Two polygons covering a cell equally resolve to the FIRST, as `centroid_zones` does."""

    small, _ = meshes
    cells = small.gdf_vorPolys
    xmin, ymin, xmax, ymax = cells.total_bounds
    whole = box(xmin, ymin, xmax, ymax)
    zones = gpd.GeoDataFrame({"k": [1.0, 2.0]}, geometry=[whole, whole], crs=cells.crs)
    painted, reached = polygon_cell_values(cells, zones, zones["k"].to_numpy(), rule="majority")
    assert reached.all()
    assert set(np.unique(painted)) == {1.0}


# ------------------------------------------------------------------------- refusals


def test_a_null_is_refused_naming_the_feature(meshes, tmp_path):
    """Ledger 247's rule: never substitute a number the user did not give."""

    small, _ = meshes
    path = tmp_path / "null.gpkg"
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["west", "east"], "k": [1.0, None]},
        geometry=[box(xmin, ymin, (xmin + xmax) / 2, ymax),
                  box((xmin + xmax) / 2, ymin, xmax, ymax)],
        crs=small.crs,
    ).to_file(path, layer="geology", driver="GPKG")

    with pytest.raises(ValueError, match="'east'"):
        array_from_polygons(
            Polygons(path=path, column="k", layer="geology"), small, nlay=1, parameter="npf.k"
        )


def test_a_cell_no_polygon_reaches_is_refused_without_a_default(meshes, tmp_path):
    """An unfilled cell is a NaN, and MF6 answers a NaN budget reporting success."""

    small, _ = meshes
    path = _write(tmp_path / "half.gpkg", small, values=((0.0, 0.5, 5.0),))
    with pytest.raises(ValueError, match="reached by no polygon"):
        array_from_polygons(
            Polygons(path=path, column="k", layer="geology"), small, nlay=1, parameter="npf.k"
        )
    filled = array_from_polygons(
        Polygons(path=path, column="k", layer="geology", default=0.5),
        small, nlay=1, parameter="npf.k",
    )
    assert set(np.unique(filled)) == {0.5, 5.0}


def test_a_map_that_lands_nowhere_is_refused_with_both_extents(meshes, tmp_path):
    """Ledger 211's shape, so a CRS or location mismatch reads as one."""

    small, _ = meshes
    path = tmp_path / "far.gpkg"
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["far"], "k": [1.0]},
        geometry=[box(xmax + 9e5, ymax + 9e5, xmax + 1e6, ymax + 1e6)],
        crs=small.crs,
    ).to_file(path, layer="geology", driver="GPKG")

    with pytest.raises(ValueError, match="lands on the grid") as error:
        array_from_polygons(
            Polygons(path=path, column="k", layer="geology", default=1.0),
            small, nlay=1, parameter="npf.k",
        )
    assert "span" in str(error.value)


def test_a_partial_miss_warns_by_name(meshes, tmp_path):
    small, _ = meshes
    path = tmp_path / "partial.gpkg"
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["on", "off"], "k": [1.0, 2.0]},
        geometry=[box(xmin, ymin, xmax, ymax),
                  box(xmax + 9e5, ymax + 9e5, xmax + 1e6, ymax + 1e6)],
        crs=small.crs,
    ).to_file(path, layer="geology", driver="GPKG")

    with pytest.warns(UserWarning, match="'off'"):
        array_from_polygons(
            Polygons(path=path, column="k", layer="geology"), small, nlay=1, parameter="npf.k"
        )


def test_a_missing_column_says_which_columns_there_are(meshes, geology):
    small, _ = meshes
    with pytest.raises(ValueError, match="did you mean 'k'"):
        array_from_polygons(
            Polygons(path=geology, column="kk", layer="geology"),
            small, nlay=1, parameter="npf.k",
        )


@pytest.mark.parametrize(
    ("helper", "parameter"), [("npf", "icelltype"), ("sto", "iconvert")]
)
def test_an_integer_flag_refuses_a_polygon_map(helper, parameter, geology):
    """Every rule returns a float, and there is no cell type 0.4."""

    recipe = mf.Polygons(path=geology, column="k", layer="geology")
    kwargs = {parameter: recipe} | ({"k": 1.0} if helper == "npf" else {})
    with pytest.raises(TypeError, match="integer flag"):
        getattr(mf, helper)(**kwargs)


def test_a_parameter_that_is_not_an_array_refuses_a_polygon_map(geology):
    recipe = mf.Polygons(path=geology, column="k", layer="geology")
    with pytest.raises(TypeError, match="not a per-cell array"):
        mf.npf(k=1.0, alternative_cell_averaging=recipe)


@pytest.mark.parametrize("helper", ["npf", "ic", "sto"])
def test_a_surface_in_an_array_slot_is_refused_at_the_declaration(helper):
    """It was ACCEPTED and died at build inside FloPy naming neither the slot nor the type."""

    surface = mf.Surface.flat(10.0)
    kwargs = {"npf": {"k": surface}, "ic": {"strt": surface}, "sto": {"sy": surface}}[helper]
    with pytest.raises(TypeError, match="not supported"):
        getattr(mf, helper)(**kwargs)


def test_a_recipe_and_an_option_for_one_parameter_is_refused_clearly(geology):
    """It used to reach FloPy as `got multiple values for keyword argument 'k'`."""

    import flopy

    from myflopy.specs import PackageSpec

    spec = PackageSpec(
        "npf", flopy.mf6.ModflowGwfnpf, {"k": 10.0},
        inputs={"k": mf.Polygons(path=geology, column="k", layer="geology")},
    )
    with pytest.raises(TypeError, match="supplied twice"):
        spec.build(object())


def test_a_layer_beyond_the_model_is_refused(meshes, tmp_path):
    small, _ = meshes
    path = _write(tmp_path / "deep.gpkg", small, layer=4)
    with pytest.raises(ValueError, match="does not have"):
        array_from_polygons(
            Polygons(path=path, column="k", layer="geology", layer_field="layer"),
            small, nlay=2, parameter="npf.k",
        )


@pytest.mark.parametrize(
    ("layers", "why"),
    [
        ([0, 1], "0, a layer this model does not have -- with layer_base=1 they run 1 to 2"),
        ([1, 1.5], "1.5 is not a whole layer number"),
        ([1, None], "empty"),
    ],
)
def test_a_feature_naming_no_layer_is_refused_by_name(meshes, tmp_path, layers, why):
    """Below the base, a fraction or a NULL names no layer, and is refused naming the feature.

    Only a layer PAST the last was refused. One below `layer_base` matched no layer
    and painted nothing: measured, zones at `lay` 0 and 1 with `layer_base=1` and
    `default=1.0` gave the first zone's cells the default with no warning (ledger
    393). A fraction matched no layer the same way; a NULL died in numpy.
    """

    small, _ = meshes
    path = _write(tmp_path / "low.gpkg", small)
    frame = gpd.read_file(path, layer="geology")
    frame["layer"] = layers
    frame.to_file(path, layer="geology", driver="GPKG")
    position = next(i for i, one in enumerate(layers) if one != 1)
    with pytest.raises(ValueError, match="1 of 2 features name no model layer") as caught:
        array_from_polygons(
            Polygons(path=path, column="k", layer="geology", layer_field="layer", default=1.0),
            small, nlay=2, parameter="npf.k",
        )
    assert f"'unit{position}' (" in str(caught.value)
    assert why in str(caught.value)


def test_a_missing_file_is_refused_at_the_declaration(tmp_path):
    with pytest.raises(FileNotFoundError):
        mf.npf(k=mf.Polygons(path=tmp_path / "absent.gpkg", column="k"))


# ------------------------------------------------------------------- shape and layers


def test_an_areal_map_is_broadcast_to_every_layer(meshes, geology):
    """No `layer_field` means a map digitized in plan view, which is every layer."""

    small, _ = meshes
    array = array_from_polygons(
        Polygons(path=geology, column="k", layer="geology"), small, nlay=3, parameter="npf.k"
    )
    assert array.shape == (3, small.ncpl)
    assert np.array_equal(array[0], array[2])


def test_layer_field_paints_only_its_own_layer(meshes, tmp_path):
    small, _ = meshes
    path = _write(tmp_path / "bylayer.gpkg", small, values=((0.0, 1.0, 7.0),), layer=2)
    array = array_from_polygons(
        Polygons(path=path, column="k", layer="geology", layer_field="layer", default=0.25),
        small, nlay=3, parameter="npf.k",
    )
    assert set(np.unique(array[1])) == {7.0}
    assert set(np.unique(array[0])) == {0.25}
    assert set(np.unique(array[2])) == {0.25}


def test_ic_takes_a_polygon_map_too(meshes, geology, tmp_path):
    small, _ = meshes
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small),
               packages=[mf.disv(_stack()), mf.npf(k=1.0),
                         mf.ic(strt=mf.Polygons(path=geology, column="k", layer="geology"))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    strt = sim.build(tmp_path / "ic").simulation.get_model("m").get_package("ic").strt.get_data()
    assert set(np.unique(strt)) == {1.0, 100.0}


# ------------------------------------------------------------ serialization and digest


def test_the_recipe_round_trips_as_itself(geology):
    recipe = mf.Polygons(
        path=geology, column="k", layer="geology", rule="mean", default=2.0,
        layer_field="layer", layer_base=0,
    )
    revived = source_from_dict(recipe.to_dict())
    assert revived == recipe
    assert type(revived) is Polygons


def test_a_default_valued_field_is_omitted_from_the_payload(geology):
    """So adding a defaulted field later cannot move a stored digest (ledger 273)."""

    payload = mf.Polygons(path=geology, column="k").to_dict()
    assert set(payload) == {"kind", "path", "column"}


def test_the_digest_reads_content_not_location(meshes, tmp_path):
    """Ledger 293's rule: a copy digests the same, a different file does not."""

    small, _ = meshes
    first = _write(tmp_path / "a.gpkg", small)
    copied = tmp_path / "b.gpkg"
    copied.write_bytes(first.read_bytes())
    other = _write(tmp_path / "c.gpkg", small, values=((0.0, 1.0, 42.0),))

    digests = [
        package_digest(mf.npf(k=mf.Polygons(path=path, column="k", layer="geology")))
        for path in (first, copied, other)
    ]
    assert digests[0] == digests[1], "a byte-identical copy must not move the digest"
    assert digests[0] != digests[2], "a different file at a different path must move it"


def test_the_eager_declaration_did_not_move(meshes):
    """`mf.npf(k=10.0)` must be what it always was -- same builder, no inputs."""

    import flopy

    spec = mf.npf(k=10.0)
    assert spec.builder is flopy.mf6.ModflowGwfnpf
    assert spec.inputs == {}
    assert spec.options["k"] == 10.0


def test_every_source_class_is_registered(meshes):
    """The trap ``BoundarySource``'s docstring names, as a sweep rather than a spot check.

    ``_spec_value`` keeps an ``endswith("Source")`` arm beside the registry so an
    unregistered kind RAISES -- but nothing asserted the registry was COMPLETE, so a
    new source class could reach a released version registered nowhere and be found
    only by the first user who saved a project containing one.
    """

    from myflopy.sources import _SOURCE_TYPES, DataSourceSpec

    def descendants(cls):
        for child in cls.__subclasses__():
            yield child
            yield from descendants(child)

    missing = sorted(
        cls.__name__ for cls in descendants(DataSourceSpec) if cls.kind not in _SOURCE_TYPES
    )
    assert not missing, f"source classes missing from _SOURCE_TYPES: {missing}"


def test_a_multi_part_polygon_is_judged_by_its_TOTAL_area_in_a_cell():
    """Summing per (cell, unit) before picking the winner, not largest single PIECE.

    A unit digitized as two lobes reaches one cell as two overlay pieces. Picking the
    largest PIECE lets a one-part neighbour it actually outcovers win: here unit A holds
    60% of the cell in two pieces of 30%, unit B 40% in one, and the piece rule gives
    the cell to B. Nothing else in this file can see the difference, because every other
    fixture is single-part.
    """

    from shapely.geometry import MultiPolygon

    cells = gpd.GeoDataFrame({"_": [0]}, geometry=[box(0, 0, 10, 10)], crs="EPSG:2927")
    zones = gpd.GeoDataFrame(
        {"k": [1.0, 2.0]},
        geometry=[
            MultiPolygon([box(0, 0, 10, 3), box(0, 7, 10, 10)]),  # 60, in two pieces
            box(0, 3, 10, 7),                                     # 40, in one
        ],
        crs="EPSG:2927",
    )
    painted, reached = polygon_cell_values(cells, zones, zones["k"].to_numpy(), rule="majority")
    assert reached.all()
    assert painted[0] == 1.0


# ---------------------------------------------------------------- the frozen snapshot


def _built(grid, package, tmp_path, name):
    return mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid),
               packages=[mf.disv(_stack()), package, mf.ic(strt=5.0)]),
        name="s", packages=[mf.tdis(nper=1)],
    )


def test_the_build_records_the_RESOLVED_array_not_the_recipe(meshes, geology, tmp_path):
    """``BuiltModel.resolved`` is what the snapshot is written from (ledger 270).

    Recording this anywhere else looks equivalent and is dead: the freeze reads exactly
    this attribute, so a package recorded elsewhere is frozen as its RECIPE and a
    reopened run goes back to the polygon file. Found that way, in this commit.
    """

    small, _ = meshes
    sim = _built(small, mf.npf(k=mf.Polygons(path=geology, column="k", layer="geology")),
                 tmp_path, "resolved")
    resolved = sim.build(tmp_path / "resolved").built.models["m"].resolved

    assert "npf" in resolved, f"only {sorted(resolved)} were recorded"
    frozen = resolved["npf"]
    assert isinstance(frozen.options["k"], np.ndarray)
    assert frozen.options["k"].shape == (1, small.ncpl)
    assert frozen.inputs == {}, "the recipe must not survive into the resolved spec"


def test_a_frozen_run_thaws_with_its_polygon_file_deleted(meshes, geology, tmp_path):
    """The whole point of the snapshot: a run rebuilds without its sources (ledger 251)."""

    from myflopy.project.freeze import freeze_simulation, thaw_simulation

    small, _ = meshes
    copied = tmp_path / "throwaway.gpkg"
    copied.write_bytes(geology.read_bytes())
    sim = _built(small, mf.npf(k=mf.Polygons(path=copied, column="k", layer="geology")),
                 tmp_path, "freeze")
    built = sim.build(tmp_path / "freeze").built
    freeze_simulation(
        sim, tmp_path / "snap", run_name="s",
        contexts={name: model.context for name, model in built.models.items()},
        resolved={name: dict(model.resolved or {}) for name, model in built.models.items()},
    )

    copied.unlink()
    thawed = thaw_simulation(tmp_path / "snap")
    packages = {package.name: package for package in thawed.models[0].packages}
    assert set(np.unique(packages["npf"].options["k"])) == {1.0, 100.0}


# -------------------------------------------------- the subsystems it has to reach


def _project_with_a_polygon_map(root, grid, source, name="p"):
    project = mf.Project(root, name=name)
    project.add_grid("mesh", grid)
    sim = mf.simulation(
        mf.gwf("m", grid=mf.grid_ref("mesh"),
               packages=[mf.disv(_stack()),
                         mf.npf(k=mf.Polygons(path=source, column="k", layer="geology")),
                         mf.ic(strt=5.0)]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    project.add_simulation(sim)
    return project, sim


def test_bundling_archives_the_polygon_file(meshes, tmp_path):
    """An archive that cannot open is ledger 293's failure; this is the new carrier."""

    small, _ = meshes
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    source = _write(outside / "geology.gpkg", small)
    project, _ = _project_with_a_polygon_map(tmp_path / "proj", small, source)

    report = project.save(bundle_sources=True)
    assert report.not_bundled == ()
    assert any("geology.gpkg" in str(destination) for _, destination in report.bundled)
    assert (project.root / "inputs" / "geology.gpkg").exists()


def test_drift_reads_the_polygon_file_and_notices_an_edit(meshes, tmp_path):
    """Ledger 294's ``source`` column, on the new source type.

    It works because the recipe rides in ``inputs`` as a ``DataSourceSpec`` carrying a
    ``source_sha`` -- the existing walk finds it with no per-type wiring. Pinned because
    nothing here would fail if a refactor stopped it being found.
    """

    small, _ = meshes
    source = _write(tmp_path / "geology.gpkg", small)
    project, sim = _project_with_a_polygon_map(tmp_path / "proj", small, source)
    project.prepare_run("r1", sim).build()

    assert set(project.runs.drift()["source"]) == {"ok"}

    edited = gpd.read_file(source, layer="geology")
    edited.loc[0, "k"] = 999.0
    edited.to_file(source, layer="geology", driver="GPKG")

    assert set(project.runs.drift()["source"]) == {"edited"}
