"""Move it, and read the report.

`sim.with_grid(new)` replaces `replace_grid`, which swapped the grid and re-mapped
NOTHING -- every package kept cellids resolved on the old mesh and nothing looked
(ledger 249). `with_grid` drops everything derived from the old grid, classifies each
package, and STAMPS the ones nothing can verify so the build refuses them.

**It reports and stamps; the BUILD refuses.** That is deliberate (user decision,
2026-09-15): the guard that refuses already exists and runs against what was actually
built, so putting the refusal there means one authority rather than two. `move_report`
is what lets you see the verdict without waiting for `prepare_run` to raise.

The four verdicts: `rebuilt` (a recipe resolves itself), `stale` (pinned to another
mesh, or a hand-written DISV that IS the old mesh), `unverifiable` (names cells with
no pin -- you cannot pin cell numbers the caller typed -- OR its pin was VOIDED
because an edit replaced the records the library resolved, ledger 275), `unchanged`
(names no cells, or asserted mesh-independent with `grid_free=`).
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.grid_identity import MOVED_FROM
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)
from myflopy.specs import ModelContext


@pytest.fixture
def boundary(tmp_path):
    path = tmp_path / "dom.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(path)
    return path


def _grid_spec(boundary, max_area):
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=max_area
    )


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack


def _lazy_simulation(boundary, max_area=40_000.0, name="s"):
    model = mf.gwf(
        "m",
        packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=90.0),
                  mf.oc(head_filerecord="m.hds")],
    ).with_grid(_grid_spec(boundary, max_area))
    return mf.simulation(model, name=name,
                         packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


# --- the move --------------------------------------------------------------- #


def test_a_lazy_model_moves_to_another_grid_and_builds(tmp_path, boundary):
    """The deliverable, stated as plainly as it can be."""

    spec = _lazy_simulation(boundary)
    before = spec.build_flopy(tmp_path / "a").built_model("m").model.get_package("disv")

    moved = spec.with_grid(_grid_spec(boundary, 12_000.0), name="s_fine")

    after = moved.build_flopy(tmp_path / "b").built_model("m").model.get_package("disv")
    assert int(before.ncpl.get_data()) == 28
    assert int(after.ncpl.get_data()) == 92
    assert moved.name == "s_fine"


def test_with_grid_drops_everything_derived_from_the_old_mesh(tmp_path, boundary):
    """Keeping them is how a swap used to build clean on numbers that meant something else."""

    grid = _grid_spec(boundary, 40_000.0).resolve(workspace=tmp_path / "g")
    context = ModelContext(grid=grid, domain=np.ones((2, grid.ncpl), dtype=int))
    spec = mf.simulation(
        mf.gwf("m", context=context,
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=90.0),
                         mf.oc(head_filerecord="m.hds")]),
        name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )

    moved = spec.with_grid(_grid_spec(boundary, 12_000.0))

    context = moved.models[0].context
    assert context.grid is None and context.domain is None and context.surfaces is None


def test_the_lineage_records_the_move(tmp_path, boundary):
    moved = _lazy_simulation(boundary).with_grid(_grid_spec(boundary, 12_000.0))

    assert moved.lineage[-1]["operation"] == "with_grid"
    assert moved.lineage[-1]["models"] == ["m"]


# --- the report ------------------------------------------------------------- #


def _mixed_simulation(tmp_path, boundary):
    """One package of each verdict, on an EAGER grid so there is a pin to move from."""

    grid = _grid_spec(boundary, 40_000.0).resolve(workspace=tmp_path / "g")
    source = tmp_path / "bc.gpkg"
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["w"], "layer": [1], "head": [9.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)], crs=grid.crs,
    ).to_file(source, layer="chd", driver="GPKG")
    context = ModelContext(grid=grid, domain=np.ones((2, grid.ncpl), dtype=int))
    model = mf.gwf(
        "m", context=context,
        packages=[
            mf.disv(_stack()),                                        # rebuilt
            mf.npf(k=1.0),                                            # unchanged
            mf.chd.gpkg(source, layer="chd", nper=1, head="head", context=context),
            mf.wel(stress_period_data={0: [[(0, 3), -50.0]]}),        # unverifiable
            mf.ic(strt=90.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    return mf.simulation(model, name="s",
                         packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def test_the_report_classifies_every_package(tmp_path, boundary):
    verdicts = {
        row["package"]: row["how"]
        for row in _mixed_simulation(tmp_path, boundary).move_report()
    }

    assert verdicts["disv"] == "rebuilt"
    assert verdicts["npf"] == "unchanged"
    assert verdicts["chd"] == "stale"
    assert verdicts["wel"] == "unverifiable"


def test_a_hand_written_disv_is_stale_not_unchanged(tmp_path, boundary):
    """It names no cellid and carries no pin, so both other tests miss it.

    It IS the old mesh -- it holds that grid's `cell2d` and `vertices` verbatim. The
    first version of the classifier called it "unchanged".
    """

    grid = _grid_spec(boundary, 40_000.0).resolve(workspace=tmp_path / "g")
    context = ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int))
    spec = mf.simulation(
        mf.gwf("m", context=context,
               packages=[mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
                         mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds")]),
        name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )

    verdicts = {row["package"]: row["how"] for row in spec.move_report()}

    assert verdicts["disv"] == "stale"


# --- the stamp, and what it refuses ----------------------------------------- #


def test_an_unverifiable_package_is_stamped_and_the_build_refuses_it(tmp_path, boundary):
    """The only thing standing behind a hand-typed cellid across a move.

    Every one of its numbers is still in range on the larger mesh and every one names
    a different cell, so neither the range check nor a pin can tell it apart from a
    correct one.
    """

    moved = _mixed_simulation(tmp_path, boundary).with_grid(_grid_spec(boundary, 12_000.0))

    wel = next(p for p in moved.models[0].packages if p.name == "wel")
    assert MOVED_FROM in (wel.metadata or {})
    with pytest.raises(ValueError, match=r"'wel' was moved to this grid by `with_grid`"):
        moved.build_flopy(tmp_path / "ws")


def test_grid_free_clears_the_stamp(tmp_path, boundary):
    """The caller's assertion that a package really is mesh-independent."""

    moved = _mixed_simulation(tmp_path, boundary).with_grid(
        _grid_spec(boundary, 12_000.0), grid_free=("wel",)
    )

    wel = next(p for p in moved.models[0].packages if p.name == "wel")
    assert MOVED_FROM not in (wel.metadata or {})
    assert moved.move_report(grid_free=("wel",))[3]["how"] == "unchanged"


# --- coupled models move together ------------------------------------------- #


def test_on_none_moves_every_model_that_has_a_grid(tmp_path, boundary):
    """Moving one and stranding its sibling is what ledger 268's check exists to catch."""

    flow = mf.gwf("flow", packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=90.0),
                                    mf.oc(head_filerecord="flow.hds")]
                  ).with_grid(_grid_spec(boundary, 40_000.0))
    transport = mf.gwt("tran", packages=[mf.disv(_stack()), mf.mst(porosity=0.3),
                                         mf.ic(strt=0.0), mf.adv(),
                                         mf.oc(concentration_filerecord="tran.ucn")]
                       ).with_grid(_grid_spec(boundary, 40_000.0))
    spec = mf.simulation(
        flow, transport, name="c",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
                  mf.ims(models=["flow"]), mf.ims(models=["tran"], name="ims_t")],
        exchanges=[mf.exchange(flow, transport)],
    )

    moved = spec.with_grid(_grid_spec(boundary, 12_000.0))

    assert moved.lineage[-1]["models"] == ["flow", "tran"]
    built = moved.build_flopy(tmp_path / "ws")
    for name in ("flow", "tran"):
        assert int(built.built_model(name).model.get_package("disv").ncpl.get_data()) == 92


def test_on_names_one_model_only(tmp_path, boundary):
    """`on=` is the escape hatch, and it is the one that can strand a sibling."""

    moved = _lazy_simulation(boundary).with_grid(_grid_spec(boundary, 12_000.0), on="m")

    assert moved.lineage[-1]["models"] == ["m"]


def test_an_edited_package_is_refused_with_a_sentence_you_can_act_on(tmp_path):
    """The old refusal was unanswerable, which is the point of ledger 275.

    An edited package kept the library's pin, so `with_grid` said "re-declare it
    against the new grid" -- to a user who may already have typed correct numbers FOR
    that grid. The message named records that no longer existed, and there was no
    route out: the edit route can never re-pin (ledger 250).

    It is ONE complaint, not two. A voided package classifies `unverifiable`, so it
    would also be stamped `MOVED_FROM` and listed a second time under a different
    explanation.
    """

    import geopandas as gpd
    from shapely.geometry import box

    from myflopy.grid_identity import GRID_PIN, PIN_VOIDED
    from myflopy.specs import fork_package

    def _grid(n):
        return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1))

    small, large = _grid(6), _grid(10)
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=20.0)
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [9.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)],
        crs=small.crs,
    ).to_file(path, layer="chd", driver="GPKG")

    resolved = mf.chd.gpkg(
        path, layer="chd", nper=1, head="head",
        context=ModelContext(grid=small, domain=np.ones((1, small.ncpl), dtype=int)),
    )
    edited = fork_package(resolved, stress_period_data={0: [[(0, 5), 9.0]]})

    assert GRID_PIN not in edited.metadata
    assert edited.metadata[PIN_VOIDED]["was"] == resolved.metadata[GRID_PIN]

    model = mf.gwf(
        "m",
        context=ModelContext(grid=large, domain=np.ones((1, large.ncpl), dtype=int)),
        packages=[mf.disv(stack), edited],
    )
    sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1)])
    moved = sim.with_grid(large, name="s2")

    rows = {row["package"]: row["how"] for row in moved.move_report()}
    assert rows["chd"] == "unverifiable", rows

    with pytest.raises(ValueError) as caught:
        moved.build()
    message = str(caught.value)

    assert "EDITED after the library resolved them" in message
    assert message.count("package 'chd'") == 1, f"one defect, one sentence:\n{message}"


def test_a_ref_declared_model_reports_what_its_packages_actually_are(tmp_path):
    """`Run.resolution` classified the `mf.ref`, not the package it names.

    A `PackageRef` carries no options, inputs or metadata, so `package_move` answered
    `unchanged` -- moves untouched -- for EVERY row of a model declared the documented
    library way. Same shape as ledger 280's freeze mis-key: the declared entry and the
    thing it resolves to are different objects, and only one of them can be classified.
    """

    from myflopy.workspace import Project

    grid = rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1))
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    result = stack.build(grid, attach=False)

    project = Project(tmp_path / "p", name="p")
    project.add_package("dis/hand", mf.disv(**result.to_disv(grid).options))
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_simulation(mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid),
               packages=[mf.ref("dis/hand"), mf.ref("npf/base")]),
        name="s", packages=[mf.tdis(nper=1)],
    ))

    run = project.prepare_run("r", "s")
    run.build()
    rows = {row["package"]: row["how"] for row in run.resolution}

    assert rows["dis/hand"] == "stale", "a hand-written DISV IS the old mesh"
    assert rows["npf/base"] == "unchanged", "a scalar npf genuinely names no cells"


# --- per-cell arrays, which nothing stamped (ledger 295) --------------------- #


def _array_simulation(grid, *, k):
    """A model carrying `k` on `grid`, declared the modern way (a LAZY disv)."""

    from myflopy.specs import ModelContext

    context = ModelContext(grid=grid, domain=np.ones((2, grid.ncpl), dtype=int))
    model = mf.gwf(
        "m",
        context=context,
        packages=[mf.disv(_stack()), mf.npf(k=k), mf.ic(strt=90.0),
                  mf.oc(head_filerecord="m.hds")],
    )
    return mf.simulation(
        model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )


def _two_meshes():
    """Two DIFFERENT meshes with the SAME cell count -- the dangerous shape.

    A range check cannot see this: every index is in range on both. It is the case
    ledger 249 was written about and the one a same-ncpl swap produces.
    """

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    b = rectangular_voronoi(CanonicalModelConfig(nrow=7, ncol=4, nlay=2, nper=1))
    assert a.ncpl == b.ncpl
    return a, b


def test_a_per_cell_array_does_not_move_silently(tmp_path):
    """The headline. `npf`/`ic`/`sto` take no `context=`, so NOTHING stamped them.

    Measured before the fix: the lazy disv rebuilt correctly on the new mesh, the
    report said `npf: unchanged -- names no cells`, the build COMPLETED, and the K
    array was byte-identical to the one written for the old mesh. Ledger 282 closed
    this for `mf.rch.array`/`mf.evt.array` only -- `CELLWISE_OPTIONS` is stamped in one
    place those two reach. B1/B2 made it MORE reachable, not less: the discretization
    used to be `stale` and refuse the build, and now it rebuilds and leaves nothing to
    object to.
    """

    a, b = _two_meshes()
    k = np.linspace(1.0, 100.0, a.ncpl * 2).reshape(2, a.ncpl)
    moved = _array_simulation(a, k=k).with_grid(mf.GridSpec.from_object(b))

    verdicts = {row["package"]: row["how"] for row in moved.move_report()}
    assert verdicts["npf"] == "unverifiable", verdicts
    assert verdicts["disv"] == "rebuilt", "the discretization still moves itself"

    with pytest.raises(ValueError) as caught:
        moved.build_flopy(tmp_path / "ws")

    assert "npf" in str(caught.value)


def test_a_grid_free_declaration_still_MOVES(tmp_path, boundary):
    """The other direction, and it is the one a careless fix breaks.

    `npf(k=[10.0, 1.0])` is a per-LAYER list and means the same on any mesh. Deciding
    from "is this a sequence?" would refuse a correct move, which is why the test is
    ncpl-EXACT and why this assertion exists beside the one above (ledger 284's
    lesson: assert BOTH directions).
    """

    a, b = _two_meshes()
    moved = _array_simulation(a, k=[10.0, 1.0]).with_grid(mf.GridSpec.from_object(b))

    verdicts = {row["package"]: row["how"] for row in moved.move_report()}
    assert verdicts["npf"] == "unchanged", verdicts

    built = moved.build_flopy(tmp_path / "ws")
    assert built.simulation.get_model("m").get_package("npf").k.array is not None


@pytest.mark.parametrize(
    "spelling",
    ["ndarray_2d", "flat_list", "nested_list", "ndarray_1d"],
    ids=["(nlay,ncpl) array", "flat ncpl list", "nested per-layer lists", "1-D ncpl array"],
)
def test_every_per_cell_SPELLING_is_seen(spelling):
    """FloPy accepts all four, so seeing only one is the same defect in miniature.

    The nested case is what a per-cell `k` on a multi-layer model looks like written
    by hand, and reading `.shape` alone missed it -- the same gap ledger 282 fixed for
    a flat list, one nesting level down.
    """

    from myflopy.specs import package_move

    a, _b = _two_meshes()
    n = a.ncpl
    payload = {
        "ndarray_2d": np.ones((2, n)),
        "flat_list": [1.0] * n,
        "nested_list": [[1.0] * n, [2.0] * n],
        "ndarray_1d": np.ones(n),
    }[spelling]

    assert package_move(mf.npf(k=payload), ncpl=n)[0] == "unverifiable"


def test_the_report_read_AFTER_a_move_agrees_with_the_build(tmp_path):
    """`with_grid` does not mutate, so it still holds the OLD mesh when it classifies.

    The spec it RETURNS does not -- its context grid is dropped. So a report read off
    the moved spec had no mesh to measure a payload against and said `unchanged` about
    the very package the build then refused. The stamp itself is the answer: once
    `with_grid` has decided, `package_move` reads its verdict rather than recomputing
    one it no longer has the inputs for.
    """

    a, b = _two_meshes()
    k = np.ones((2, a.ncpl))
    moved = _array_simulation(a, k=k).with_grid(mf.GridSpec.from_object(b))

    assert MOVED_FROM in (moved.models[0].packages[1].metadata or {})

    verdicts = {row["package"]: row["how"] for row in moved.move_report()}
    refused = [
        row["package"] for row in moved.move_report() if row["how"] == "unverifiable"
    ]

    with pytest.raises(ValueError) as caught:
        moved.build_flopy(tmp_path / "ws")

    for package in refused:
        assert package in str(caught.value), (
            f"{package} reads {verdicts[package]!r} but the build does not name it"
        )


# --- an EAGER .gpkg moves exactly, from its own recipe (ledger 296) ---------- #


def _eager_gpkg_simulation(tmp_path, grid, **overrides):
    """A model whose GHB was declared EAGERLY -- with a `context=`, as the docs teach."""

    from myflopy.specs import ModelContext

    context = ModelContext(grid=grid, domain=np.ones((2, grid.ncpl), dtype=int))
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [95.0], "cond": [2.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.3, ymax)], crs=grid.crs,
    ).to_file(path, layer="ghb", driver="GPKG")

    ghb = mf.ghb.gpkg(
        path, layer="ghb", nper=1, head="head",
        conductance=mf.Spread("cond", mode="retained"), context=context, **overrides
    )
    model = mf.gwf(
        "m", context=context,
        packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=95.0),
                  mf.oc(head_filerecord="m.hds"), ghb],
    )
    sim = mf.simulation(
        model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )
    return sim, ghb


def _ghb_records(built):
    return built.simulation.get_model("m").get_package("ghb").stress_period_data.get_data(0)


def test_an_EAGERLY_declared_gpkg_moves_EXACTLY(tmp_path):
    """The largest `stale` class there is, and it needs no overlay engine.

    Declaring with a `context=` is what every pre-B2 document teaches, and it resolves
    at declaration and holds finished records -- so a grid move used to classify it
    `stale` and the build refused it. Everything needed to resolve it AGAIN was present
    at the call, so it is recorded; the move re-declares it against the new mesh, which
    is EXACT where an overlay would only approximate.
    """

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    b = rectangular_voronoi(CanonicalModelConfig(nrow=9, ncol=9, nlay=2, nper=1))
    assert a.ncpl != b.ncpl

    sim, _ghb = _eager_gpkg_simulation(tmp_path, a)
    on_a = _ghb_records(sim.build_flopy(tmp_path / "a"))

    moved = sim.with_grid(mf.GridSpec.from_object(b))
    verdicts = {row["package"]: row["how"] for row in moved.move_report()}
    assert verdicts["ghb"] == "rebuilt", verdicts

    on_b = _ghb_records(moved.build_flopy(tmp_path / "b"))

    assert len(on_b) != len(on_a), "a different mesh resolves a different cell set"
    assert float(sum(r["cond"] for r in on_b)) == pytest.approx(
        float(sum(r["cond"] for r in on_a))
    ), "mf.Spread(mode='retained') conserves the feature total across meshes"


def test_the_recipe_records_what_METADATA_CANNOT(tmp_path):
    """`metadata['fields']` looks like a recipe and is three fields short.

    `nper`, `edges_only` and `boundnames` are absent from it, and `_metadata_value` is
    deliberately not a round-trip format (ledger 274). `edges_only` is the one that
    matters: it decides WHICH CELLS are kept (ledger 246), so reconstructing a recipe
    from metadata and guessing it would write a different package in silence.
    """

    from myflopy.grid_identity import MOVE_RECIPE

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    _sim, ghb = _eager_gpkg_simulation(tmp_path, a, edges_only=True, boundnames=False)
    recipe = ghb.metadata[MOVE_RECIPE]

    assert (recipe.nper, recipe.edges_only, recipe.boundnames) == (1, True, False)
    assert set(ghb.metadata).isdisjoint({"nper", "edges_only", "boundnames"}), (
        "none of the three is in metadata -- which is why the recipe is recorded, "
        "not re-derived"
    )
    assert type(recipe.values["conductance"]) is mf.Spread, "a LIVE object, not a dict"


def test_recording_the_recipe_moves_NO_digest(tmp_path):
    """It rides in `metadata`, which `package_digest` excludes by contract.

    If it did not, adding this would report every stored `.gpkg` package as drifted on
    the day it landed -- the reason ledger 250 kept the grid pin out of the digest too.
    """

    import dataclasses

    from myflopy.grid_identity import MOVE_RECIPE
    from myflopy.specs import package_digest

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    _sim, ghb = _eager_gpkg_simulation(tmp_path, a)
    without = dataclasses.replace(
        ghb, metadata={k: v for k, v in ghb.metadata.items() if k != MOVE_RECIPE}
    )

    assert package_digest(ghb) == package_digest(without)


def test_a_LAZY_package_is_not_relazified_twice(tmp_path):
    """A lazy declaration already rebuilds itself; re-deriving it would be a second
    source of truth for the same package."""

    from myflopy.grid_identity import MOVE_RECIPE

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    xmin, ymin, xmax, ymax = a.gdf_vorPolys.total_bounds
    path = tmp_path / "lazy.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [95.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.3, ymax)], crs=a.crs,
    ).to_file(path, layer="ghb", driver="GPKG")

    lazy = mf.ghb.gpkg(path, layer="ghb", nper=1, head="head")   # NO context=

    assert MOVE_RECIPE not in (lazy.metadata or {}), "the lazy form records no second copy"
    assert lazy.inputs, "it carries its recipe in `inputs`, where the builder reads it"


def test_the_OLD_pin_does_not_survive_the_relazify(tmp_path):
    """A pin describes the resolution that was just replaced.

    `package_move` reads `GRID_PIN_OPTIONS <= inputs` before anything else, so a stale
    pin sitting beside a fresh recipe is the exact collision ledger 275 voided -- and
    the coverage counts describe records the package no longer holds.
    """

    from myflopy.grid_identity import GRID_PIN, MOVE_RECIPE
    from myflopy.specs import _relazified

    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    _sim, ghb = _eager_gpkg_simulation(tmp_path, a)

    assert GRID_PIN in ghb.metadata and "coverage" in ghb.metadata

    moved = _relazified(ghb)

    assert GRID_PIN not in moved.metadata
    assert "coverage" not in moved.metadata
    assert MOVE_RECIPE in moved.metadata, "the recipe itself stays, for the next move"
