"""Everything that serializes a grid, tested on a grid flopy actually triangulated.

The suite's other grids come from `rectangular_voronoi`, which is also a real SciPy
Voronoi with ragged, unequal cells. The ONE thing it cannot exhibit is the dtype of
its `cell2d` vertex indices: python `int` throughout, where a `GridSpec.voronoi` grid
mixes `int` with `numpy.int64`. Ledger 263 first blamed raggedness; 265 measured it and
corrected the record. Raggedness proves nothing here -- only the dtype does.

`_json_value` passed `np.float64` (a `float` subclass) and refused `np.int64` (not an
`int` subclass), so on a real grid `to_dict` raised, `spec_recipe` was None,
`freeze_simulation` died, `Run._freeze` swallowed it as a warning, and the run could
not be rebuilt -- while the run itself executed with Normal termination. 3,315 tests
passed throughout.

This module is the one that would have failed.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)
from myflopy.project.freeze import thaw_simulation
from myflopy.specs import ModelContext, _json_value, disv_pin, grid_pin
from myflopy.workspace import Project


def _simulation(grid, name="real"):
    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    return mf.simulation(
        model, name=name, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )


def _index_kinds(grid):
    """The type names of every vertex index in ``grid``'s cell2d rows."""

    # `row` is (icell2d, xc, yc, ncvert, *iverts) -- the indices start at 4. Slicing
    # from 3 folds in `ncvert`, a plain `int` on both grids; measured, that does not
    # change either assertion's verdict, so this is precision rather than a fix.
    return {type(v).__name__ for row in grid.get_disv_gridprops()["cell2d"] for v in row[4:]}


def test_the_fixture_differs_from_the_suite_s_grid_in_the_way_that_matters(triangulated_grid):
    """The fixture earns its place only while it still carries numpy ints.

    Both halves are asserted, because only the CONTRAST is the point: if a future flopy
    normalised its output, or `rectangular_voronoi` started emitting int64 of its own,
    this module would quietly become a duplicate of the canonical-grid tests. The first
    version of this test also asserted "cells of differing widths", which reads as a
    second guard and is not one -- it measured `ncvert`, and both grids are ragged
    (ledger 265).
    """

    assert "int64" in _index_kinds(triangulated_grid), "the fixture no longer produces numpy ints"
    plain = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=4, nlay=1, nper=1))
    assert _index_kinds(plain) == {"int"}, (
        "`rectangular_voronoi` now emits numpy ints too, so it no longer hides the bug "
        "this module exists for -- re-derive what the fixture is still for, or drop it"
    )


def test_a_numpy_scalar_serializes():
    """It is a number, and whether one got through was an accident of subclassing."""

    # `type(...) is`, not `isinstance`: `np.float64` SUBCLASSES `float`, so an
    # isinstance check passes on the unconverted value and certifies the bug. That
    # is how `to_yaml` stayed broken after ledger 263 -- PyYAML dispatches on the
    # exact type where json.dumps does not (ledger 265).
    for raw, expected in [
        (np.int64(7), 7),
        (np.float64(1.5), 1.5),
        (np.float32(1.5), 1.5),
        (np.array(3), 3),
    ]:
        converted = _json_value(raw)
        assert converted == pytest.approx(expected)
        assert type(converted) is type(expected), f"{raw!r} came back as {type(converted)}"
    assert _json_value(np.bool_(True)) is True


def test_a_real_grid_model_serializes(triangulated_grid):
    """`to_dict` is the gate every other one of these stands behind."""

    assert _simulation(triangulated_grid).to_dict()["kind"] == "SimulationSpec"


def test_a_real_grid_run_is_frozen_and_rebuildable(tmp_path, triangulated_grid):
    """The claim ledger 251 makes, on the grid a user actually has."""

    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_simulation(triangulated_grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = project.prepare_run("r1", "real")
        run.built.simulation.write_simulation(silent=True)

    assert run.frozen is not None
    assert run.frozen["readable"] is True
    thawed = thaw_simulation(run.workspace, run_name="r1")
    assert thawed is not None
    assert thawed.models[0].context.grid.ncpl == triangulated_grid.ncpl


def test_a_real_grid_run_reopens_and_rebuilds(tmp_path, triangulated_grid):
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_simulation(triangulated_grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = project.prepare_run("r1", "real")
        run.built.simulation.write_simulation(silent=True)

    reopened = Project(tmp_path / "proj", name="p").reopen_run("r1")
    reopened.build()
    assert reopened.spec is not None


def test_a_real_grid_fingerprints_consistently(tmp_path, triangulated_grid):
    """`grid_pin` and `disv_pin` must agree on a mesh whose indices are numpy ints."""

    built = _simulation(triangulated_grid).build_flopy(tmp_path / "ws")
    disv = built.built_model("m").model.get_package("disv")
    assert disv_pin(disv) == grid_pin(triangulated_grid)


def test_a_real_grid_project_saves_and_reloads(tmp_path, triangulated_grid):
    """The whole chain: save, reload in a fresh Project, rebuild."""

    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_simulation(triangulated_grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = project.save()
    assert report.skipped == ()

    reloaded = Project.load(tmp_path / "proj")
    assert reloaded.prepare_run("r", "real") is not None


# --- the refusing direction, on two real ragged meshes ---------------------- #
#
# Everything above runs one real grid through the passing direction. The guards'
# refusals were only ever measured on hand-built meshes, where every ring is the
# same length -- and `_canonical_ring` and `_mesh_residual` are exactly the code a
# variable ring length could plausibly break.


def _chd_source(tmp_path, grid):
    """A one-feature CHD GeoPackage over the western fifth of ``grid``."""

    import geopandas as gpd
    from shapely.geometry import box

    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [9.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)],
        crs=grid.crs,
    ).to_file(path, layer="chd", driver="GPKG")
    return path


def test_two_real_meshes_pin_differently(triangulated_grid, triangulated_grid_finer):
    """The fixtures are only useful as a pair if they are actually different."""

    assert triangulated_grid.ncpl < triangulated_grid_finer.ncpl
    assert grid_pin(triangulated_grid) != grid_pin(triangulated_grid_finer)


def test_a_package_resolved_on_the_coarse_mesh_is_refused_on_the_fine_one(
    tmp_path, triangulated_grid, triangulated_grid_finer
):
    """Every cellid is IN RANGE on the finer grid, and every one names another cell.

    That is the range check's documented gap (ledger 249) and the whole reason the
    connectivity pin exists (ledger 250) -- tested here on meshes with ragged rings,
    where it had only been tested on uniform ones.
    """

    chd = mf.chd.gpkg(
        _chd_source(tmp_path, triangulated_grid), layer="chd", nper=1, head="head",
        context=ModelContext(
            grid=triangulated_grid,
            domain=np.ones((1, triangulated_grid.ncpl), dtype=int),
        ),
    )
    resolved = chd.options["stress_period_data"][0]
    assert max(record[0][1] for record in resolved) < triangulated_grid_finer.ncpl, (
        "the fine grid must be strictly larger, or the range check would catch this "
        "and the pin would not be what is under test"
    )

    model = mf.gwf(
        "m",
        context=ModelContext(
            grid=triangulated_grid_finer,
            domain=np.ones((1, triangulated_grid_finer.ncpl), dtype=int),
        ),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **triangulated_grid_finer.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), chd, mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="x", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])

    with pytest.raises(ValueError, match="resolved against a different mesh"):
        sim.build_flopy(workspace=tmp_path / "x")


def test_a_renumbered_real_mesh_is_refused(tmp_path, triangulated_grid):
    """The same renumbering `test_build_grid_consistency` pins, on a RAGGED mesh.

    That module proves it on cells that all have the same ring length. `_mesh_residual`
    is index-aligned and `_canonical_ring` rotates each ring to a canonical start, so a
    variable ring length is the one input that could plausibly break either -- and
    neither had ever seen one.
    """

    props = dict(triangulated_grid.get_disv_gridprops())
    ring_lengths = {len(row) - 4 for row in props["cell2d"]}
    assert len(ring_lengths) > 1, "the point of doing this here is the ragged rings"

    order = list(range(triangulated_grid.ncpl))
    order = order[7:] + order[:7]                       # a rotation of the numbering
    props["cell2d"] = [[new, *props["cell2d"][old][1:]] for new, old in enumerate(order)]

    model = mf.gwf(
        "m",
        context=ModelContext(
            grid=triangulated_grid,
            domain=np.ones((1, triangulated_grid.ncpl), dtype=int),
        ),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **props),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="r", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])

    with pytest.raises(ValueError, match=r"A renumbering counts"):
        sim.build_flopy(workspace=tmp_path / "r")


def test_a_thawed_real_grid_writes_the_same_mf6_input(tmp_path, triangulated_grid, mf6_body):
    """The written input comes from the PACKAGE options, so that is what this proves.

    Stated carefully because the first version of this docstring did not: a DISV spec
    carries its own `vertices`/`cell2d`, and those ride through the snapshot as package
    options. The frozen GRID block populates `context.grid` and reaches the written
    files not at all -- measured, every vertex in the snapshot can be shifted a full
    foot and these files are byte-identical.

    So the two halves are asserted separately: the files, which cover the package
    options and the `$npz` round trip, and the thawed grid's own vertices, which
    nothing else compares. `grid_pin` cannot stand in for the second -- it is
    connectivity-only by design, so it is blind to a vertex that moved (ledger 250).
    """

    project = Project(tmp_path / "p", name="p")
    project.add_simulation(_simulation(triangulated_grid, name="t"))
    first = project.prepare_run("one", "t")
    first.build()
    first.write()

    thawed = thaw_simulation(first.workspace)
    second = thawed.build_flopy(workspace=tmp_path / "again")
    second.sim.write_simulation()

    for name in ("m.disv", "m.npf", "m.ic"):
        # `mf6_body`, not `read_text`: FloPy stamps the wall clock into line 1 of every
        # file, so a raw comparison fails ~1% of the time -- accusing the serialization,
        # which is a false lead (ledger 266).
        original = mf6_body(first.workspace / name)
        rebuilt = mf6_body(tmp_path / "again" / name)
        assert original == rebuilt, f"{name} differs after a freeze/thaw round trip"

    thawed_grid = thawed.models[0].context.grid
    assert np.array_equal(np.asarray(thawed_grid.verts), np.asarray(triangulated_grid.verts)), (
        "the snapshot's own grid block came back with different vertices"
    )
    assert [list(ring) for ring in thawed_grid.iverts] == [
        list(ring) for ring in triangulated_grid.iverts
    ], "the snapshot's own grid block came back with different connectivity"
