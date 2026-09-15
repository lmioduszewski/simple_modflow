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
no pin -- you cannot pin cell numbers the caller typed), `unchanged` (names no cells,
or asserted mesh-independent with `grid_free=`).
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.grid_identity import MOVED_FROM
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
