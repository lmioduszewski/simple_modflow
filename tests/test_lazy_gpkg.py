"""A `.gpkg` boundary condition declared WITHOUT `context=` resolves at build time.

One declaration builds correctly on two different meshes -- the reason Phase 2 exists.
Until now every `.gpkg` helper resolved its cells at DECLARATION, so moving a model to a
new grid left every package's cellids on the old mesh (ledger 249), and `move_report`
called them `stale`. A lazy one reads `rebuilt` and resolves itself.

With `context=` nothing changes: the eager path is untouched, and every test here that
compares the two is asserting exactly that.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, box

import myflopy as mf
from myflopy.grid_identity import GRID_PIN, SOURCE_SHA, source_fingerprint
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext, package_move


def _grid(n):
    return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1))


@pytest.fixture(scope="module")
def meshes():
    """Two DIFFERENT meshes. The pair is the point -- see `test_the_two_meshes_differ`."""

    return _grid(6), _grid(10)


@pytest.fixture(scope="module")
def source(tmp_path_factory, meshes):
    """A CHD polygon and a DRN line, both spanning the western/middle of mesh A."""

    small, _ = meshes
    path = tmp_path_factory.mktemp("gpkg") / "bc.gpkg"
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [9.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)], crs=small.crs,
    ).to_file(path, layer="chd", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["mid"], "layer": [1], "cond": [10.0]},
        geometry=[LineString([(xmin, (ymin + ymax) / 2), (xmax, (ymin + ymax) / 2)])],
        crs=small.crs,
    ).to_file(path, layer="drn", driver="GPKG")
    return path


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    return stack


def _records(spec, grid, tmp_path, name):
    """Build a one-package model on `grid` and return its written records."""

    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid), packages=[mf.disv(_stack()), spec]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    built = sim.build(tmp_path / name)
    return built.simulation.get_model("m").get_package(spec.name).stress_period_data.get_data(0)


def test_a_lazy_gpkg_declares_a_recipe_and_reads_as_rebuilt(source):
    """`inputs` carries the recipe under the option it supplies (ledger 272 case b)."""

    from myflopy.sources import BoundarySource

    spec = mf.chd.gpkg(source, layer="chd", nper=1, head="head")

    assert type(spec.inputs["stress_period_data"]) is BoundarySource
    assert "stress_period_data" not in spec.options
    # NO PIN, which is what "nothing is resolved yet" actually means. This used to
    # assert the metadata was EMPTY, which is a wider claim than its own reason: the
    # spec carries provenance from ledger 293 -- where the file is and what it said --
    # and provenance is not a resolution. `metadata` is excluded from `package_digest`,
    # so neither kind of stamp moves anything.
    assert GRID_PIN not in spec.metadata, "nothing is resolved yet, so nothing is pinned"
    assert spec.metadata["source_path"] == str(source)
    assert spec.metadata[SOURCE_SHA] == source_fingerprint(source)
    assert package_move(spec)[0] == "rebuilt"


def test_the_two_meshes_differ(meshes):
    """The guard on the comparison below: 'equal on both' passes vacuously otherwise."""

    small, large = meshes

    assert small.ncpl != large.ncpl


@pytest.mark.parametrize("mesh", [0, 1])
def test_one_declaration_equals_the_eager_form_on_either_mesh(
    source, meshes, tmp_path, mesh
):
    """The contract: lazy == eager, on A and on B, for the SAME declaration.

    Compared as the WRITTEN records rather than as a digest -- a digest would differ for
    reasons that have nothing to do with the cells, and comparing counts alone would
    pass on a resolver that found the right number of wrong cells.
    """

    grid = meshes[mesh]
    lazy = mf.drn.gpkg(
        source, layer="drn", nper=1,
        elevation=mf.CellSurfaceOffset("cell_bottom", offset=2.0),
        conductance=mf.Spread("cond"),
    )
    written = _records(lazy, grid, tmp_path, f"lazy{mesh}")

    resolved = mf.disv(_stack().build(grid, attach=False))
    eager = mf.drn.gpkg(
        source, layer="drn", nper=1,
        context=ModelContext(
            grid=grid, domain=_stack().build(grid, attach=False).idomain
        ),
        elevation=mf.CellSurfaceOffset("cell_bottom", offset=2.0),
        conductance=mf.Spread("cond"),
    )
    assert resolved is not None
    expected = eager.options["stress_period_data"][0]

    assert len(written) == len(expected)
    assert [tuple(r)[0] for r in written] == [tuple(r)[0] for r in expected], "same cells"
    assert np.allclose(
        [tuple(r)[1] for r in written], [tuple(r)[1] for r in expected]
    ), "same elevations -- CellSurfaceOffset resolved off the built discretization"


def test_the_two_meshes_resolve_to_different_records(source, meshes, tmp_path):
    """Without this, 'lazy equals eager on both' would pass on a resolver doing nothing."""

    small, large = meshes
    spec = mf.chd.gpkg(source, layer="chd", nper=1, head="head")

    assert len(_records(spec, small, tmp_path, "a")) != len(
        _records(spec, large, tmp_path, "b")
    )


def test_the_domain_comes_from_the_built_discretization(source, meshes, tmp_path):
    """Correctness, not politeness.

    `with_grid` DROPS `context.domain`, so a lazy package's context routinely carries
    none. A resolver reading only the context writes boundaries into cells the
    discretization makes inactive -- measured here, 10 records against a correct 7 --
    and MF6 then refuses the run, except under `edges_only=True` where it runs normally
    on the wrong boundary.
    """

    small, _ = meshes
    result = _stack().build(small, attach=False)
    line = gpd.read_file(source, layer="drn").geometry.iloc[0]
    crossed = small.gdf_vorPolys[small.gdf_vorPolys.intersects(line)].index.tolist()

    idomain = result.idomain.copy()
    idomain[0, crossed[:3]] = 0
    disv = mf.disv(**{**result.to_disv(small).options, "idomain": idomain})

    lazy = mf.drn.gpkg(
        source, layer="drn", nper=1, elevation=10.0, conductance=mf.Spread("cond")
    )
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small), packages=[disv, lazy]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    built = sim.build(tmp_path / "dom")
    written = built.simulation.get_model("m").get_package("drn").stress_period_data.get_data(0)

    blind = mf.drn.gpkg(
        source, layer="drn", nper=1, elevation=10.0, conductance=mf.Spread("cond"),
        context=ModelContext(grid=small),
    )
    aware = mf.drn.gpkg(
        source, layer="drn", nper=1, elevation=10.0, conductance=mf.Spread("cond"),
        context=ModelContext(grid=small, domain=idomain),
    )

    assert len(written) == len(aware.options["stress_period_data"][0])
    assert len(written) < len(blind.options["stress_period_data"][0]), (
        "the inactive cells must actually have been dropped"
    )


@pytest.mark.parametrize(
    "package, kwargs",
    [
        ("chd", {"head": "head"}),
        ("ghb", {"head": "head", "conductance": 1.0}),
        ("drn", {"elevation": 10.0, "conductance": 1.0}),
        ("riv", {"stage": 12.0, "conductance": 1.0, "rbot": 8.0}),
        ("wel", {"rate": -1.0}),
    ],
)
def test_every_list_bc_has_a_lazy_form(source, meshes, tmp_path, package, kwargs):
    """All seven wrappers, not just the one the examples use."""

    layer = "chd" if package in ("chd", "ghb") else "drn"
    spec = getattr(mf, package).gpkg(source, layer=layer, nper=1, **kwargs)

    assert package_move(spec)[0] == "rebuilt"
    assert len(_records(spec, meshes[0], tmp_path, f"{package}")) > 0


def test_editing_a_lazy_package_is_refused_at_the_edit(source):
    """The source exemption had to learn a limit when a builder RESOLVES the recipe.

    `_check_input_collisions` exempts a source spec in `inputs` because shadowing its
    own option is what a source MEANS (ledger 272). But `_build_lazy_list_bc` takes
    `stress_period_data` as a required keyword-only parameter, so it is handed BOTH --
    and both edit routes constructed fine and died at BUILD with
    `functools.partial(...) got multiple values for keyword argument`, naming a partial
    object. That is the failure the check exists to prevent, let through by its own
    exemption.
    """

    from myflopy.specs import fork_package

    lazy = mf.chd.gpkg(source, layer="chd", nper=1, head="head")
    records = {0: [[(0, 0), 1.0]]}

    with pytest.raises(TypeError, match="part of this package's recipe"):
        lazy.with_options(stress_period_data=records)

    with pytest.raises(TypeError, match="part of this package's recipe"):
        fork_package(lazy, stress_period_data=records)

    # The exemption itself survives: a source shadowing an option whose builder does
    # NOT require it is still legal -- that is ledger 272's case (b), still unwired.
    import dataclasses

    from myflopy.sources import DataSourceSpec

    assert dataclasses.replace(
        mf.npf(k=1.0), inputs={"k": DataSourceSpec(path="x.tif")}
    ) is not None


def test_a_pass_through_edit_still_works_on_a_lazy_package(source):
    """The control: only the option the recipe SUPPLIES is refused."""

    from myflopy.specs import fork_package

    lazy = mf.chd.gpkg(source, layer="chd", nper=1, head="head")
    renamed = fork_package(lazy, name="chd_west", save_flows=True)

    assert renamed.name == "chd_west"
    assert renamed.options["save_flows"] is True
    assert type(renamed.inputs["stress_period_data"]).__name__ == "BoundarySource"


def test_a_rename_reaches_the_built_package(source, meshes, tmp_path):
    """Why the recipe carries no name (ledgers 277, 278).

    A name frozen in the recipe survives a rename, and two renamed forks then collapse
    into ONE MF6 package under a FloPy "Replacing existing package" warning while
    myflopy reports two.
    """

    from myflopy.specs import fork_package

    small, _ = meshes
    west = fork_package(
        mf.chd.gpkg(source, layer="chd", nper=1, head="head"), name="chd_west"
    )
    sim = mf.simulation(
        mf.gwf(
            "m", context=ModelContext(grid=small),
            packages=[mf.disv(_stack()), west],
        ),
        name="s", packages=[mf.tdis(nper=1)],
    )
    model = sim.build(tmp_path / "rename").simulation.get_model("m")

    assert "CHD_WEST" in model.get_package_list()


def _project_run(tmp_path, source, grid, *, use_ref):
    """Prepare and build a one-BC run, either inline or through `mf.ref`."""

    from myflopy.workspace import Project

    project = Project(tmp_path / "proj", name="p")
    lazy = mf.chd.gpkg(source, layer="chd", nper=1, head="head")
    if use_ref:
        project.add_package("chd/west", lazy)
        packages = [mf.disv(_stack()), mf.ref("chd/west")]
    else:
        packages = [mf.disv(_stack()), lazy]
    project.add_simulation(
        mf.simulation(
            mf.gwf("m", context=ModelContext(grid=grid), packages=packages),
            name="s", packages=[mf.tdis(nper=1)],
        )
    )
    run = project.prepare_run("r", "s")
    run.build()
    return run


@pytest.mark.parametrize("use_ref", [False, True], ids=["inline", "mf.ref"])
def test_a_lazy_run_rebuilds_with_its_geopackage_deleted(
    source, meshes, tmp_path, use_ref
):
    """The snapshot must hold RECORDS, not the recipe that points back at the file.

    `mf.ref` is parametrized because it is the DOCUMENTED library path and it failed
    where inline succeeded: `freeze._for` reads `name` off the DECLARED entry, and a
    `PackageRef`'s name is its library KEY ("chd/west"), while a builder only knows the
    resolved spec's own name ("chd"). The mismatch is completely silent -- the snapshot
    simply keeps the recipe, and nothing says so until the file moves.
    """

    import shutil

    from myflopy.project.freeze import thaw_simulation

    # A PRIVATE copy: the fixture is module-scoped, and deleting it broke every test
    # that ran after this one. A shared fixture is read-only by contract (ledger 48).
    private = tmp_path / "gis"
    private.mkdir()
    shutil.copy(source, private / source.name)

    small, _ = meshes
    run = _project_run(tmp_path, private / source.name, small, use_ref=use_ref)
    before = (
        run.built.built_model("m").model.get_package("chd")
        .stress_period_data.get_data(0)
    )

    shutil.rmtree(private)

    thawed = thaw_simulation(run.workspace)
    rebuilt = thawed.build_flopy(workspace=tmp_path / "again")
    after = rebuilt.built_model("m").model.get_package("chd").stress_period_data.get_data(0)

    assert len(after) == len(before) > 0
    assert [tuple(r)[0] for r in after] == [tuple(r)[0] for r in before]


def test_the_snapshot_holds_records_and_the_declared_spec_keeps_its_recipe(
    source, meshes, tmp_path
):
    """Two different things, and each has to stay itself.

    The snapshot is the resolved EAGER equivalent -- pinned, `stale`, no recipe -- so a
    reopened run rebuilds exactly what ran. The library's spec keeps the recipe and
    stays `rebuilt`, so it can still move to another mesh.
    """

    from myflopy.grid_identity import GRID_PIN, GRID_PIN_OPTIONS, grid_pin
    from myflopy.project.freeze import thaw_simulation

    small, _ = meshes
    declared = mf.chd.gpkg(source, layer="chd", nper=1, head="head")
    run = _project_run(tmp_path, source, small, use_ref=False)
    frozen = next(
        p for p in thaw_simulation(run.workspace).models[0].packages
        if getattr(p, "name", None) == "chd"
    )

    assert not frozen.inputs, "the recipe must not survive into the snapshot"
    assert len(frozen.options["stress_period_data"][0]) > 0
    assert frozen.metadata[GRID_PIN] == grid_pin(small)
    assert frozen.metadata[GRID_PIN_OPTIONS] == ["stress_period_data"]
    assert package_move(frozen)[0] == "stale", "a resolved snapshot IS pinned"

    assert declared.inputs, "the declared spec keeps its recipe"
    assert package_move(declared)[0] == "rebuilt"


def test_the_resolution_report_says_what_each_package_resolved_to(
    source, meshes, tmp_path
):
    """`ncpl` is what a lazy DISV resolves to; a lazy BC resolves to RECORDS.

    Reading `options["ncpl"]` off a boundary condition gave `None` -- a column that
    said nothing for every package except the discretization.
    """

    small, _ = meshes
    run = _project_run(tmp_path, source, small, use_ref=False)
    rows = {row["package"]: row for row in run.resolution}

    assert rows["disv"]["resolved_ncpl"] == small.ncpl
    assert rows["chd"]["resolved_records"] > 0
    assert "resolved_ncpl" not in rows["chd"], "a BC resolves to records, not a mesh"


def test_a_RELATIVE_declaration_names_the_file_the_caller_meant(tmp_path, monkeypatch):
    """Two bases for one string, and the build silently took the other one.

    `_lazy_bc_spec` validated the path against the CURRENT WORKING DIRECTORY while
    `_anchored_packages` re-based it against the PROJECT ROOT -- so a declaration that
    passed could fail at build, and where a same-named file happened to sit under the
    project root, the model BUILT CLEAN from the wrong one. Measured before the fix:
    heads of 999.0 against the declared 9.0, normal termination, no warning. Found by
    two independent adversarial reviews.

    The recipe resolves at declaration now, so the only relative paths reaching the
    anchor are the ones the bundler wrote, which are project-relative by construction.
    """

    import geopandas as gpd
    from shapely.geometry import box

    from myflopy.specs import SpecBuildContext, _anchored_packages

    here, root = tmp_path / "here", tmp_path / "proj"
    (here / "gis").mkdir(parents=True)
    (root / "gis").mkdir(parents=True)

    def write(path, head):
        gpd.GeoDataFrame(
            {"name": ["w"], "layer": [1], "head": [head]},
            geometry=[box(0, 0, 10, 10)], crs=2927,
        ).to_file(path, layer="chd", driver="GPKG")

    write(here / "gis" / "bc.gpkg", 9.0)        # what the caller is looking at
    write(root / "gis" / "bc.gpkg", 999.0)      # the decoy the anchor used to pick

    monkeypatch.chdir(here)
    spec = mf.chd.gpkg("gis/bc.gpkg", layer="chd", nper=1, head="head")

    assert spec.inputs["stress_period_data"].path == (here / "gis" / "bc.gpkg").resolve()

    anchored = _anchored_packages((spec,), SpecBuildContext(project_root=root))

    assert anchored[0].inputs["stress_period_data"].path == (
        here / "gis" / "bc.gpkg"
    ).resolve(), "the anchor must not repoint a path the declaration already resolved"

    # And the hash describes the file the build will actually open.
    assert spec.inputs["stress_period_data"].source_sha == source_fingerprint(
        here / "gis" / "bc.gpkg"
    )
