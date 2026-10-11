"""A model whose packages were resolved on another mesh is refused at build.

`replace_grid` swapped a grid and left every package's baked cellids and dimension
counts on the old mesh; nothing looked, so MF6 ran the old cell numbers against the
new geometry and reported Normal termination. Measured on that path: a lazy CHD
raised IndexError on a finer grid and resolved to ZERO cells on a coarser one.

The guard is a range-and-identity check, not a correctness proof: old cellids that
all happen to fall inside a larger new grid still pass, which is what the grid pin
is for. `test_a_refined_swap_is_the_documented_gap` pins that limit so closing it
later is a deliberate act.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext, _cellids_in, _discretization_complaints


def _grid(nrow: int, ncol: int):
    """A rectangular Voronoi grid of ``nrow x ncol`` cells."""

    return rectangular_voronoi(CanonicalModelConfig(nrow=nrow, ncol=ncol, nlay=1, nper=1))


def _build(tmp_path, tag, disv_grid, ctx_grid, *, extra=(), domain=None, surfaces=None):
    """Build a minimal model whose DISV comes from one grid and whose context from another."""

    if domain is None:
        domain = np.ones((1, ctx_grid.ncpl), dtype=int)
    context = ModelContext(grid=ctx_grid, domain=domain, surfaces=surfaces)
    model = mf.gwf(
        "m",
        context=context,
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **disv_grid.get_disv_gridprops()),
            mf.npf(k=1.0),
            mf.ic(strt=5.0),
            *extra,
            mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(
        model, name=tag, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )
    return sim.build_flopy(workspace=tmp_path / tag)


# ---------------------------------------------------------------- the reproductions


def test_a_coarser_grid_than_the_discretization_is_refused(tmp_path):
    """100 cells of DISV against a 300-cell context grid."""

    with pytest.raises(ValueError, match=r"DISV package has ncpl=100 .* has 300 cells"):
        _build(tmp_path, "a", _grid(10, 10), _grid(20, 15))


def test_a_finer_grid_than_the_discretization_is_refused(tmp_path):
    """And the other direction, which used to resolve to zero cells."""

    with pytest.raises(ValueError, match=r"DISV package has ncpl=300 .* has 100 cells"):
        _build(tmp_path, "b", _grid(20, 15), _grid(10, 10))


def test_the_same_cell_count_is_not_the_same_mesh(tmp_path):
    """A 3x2 and a 2x3 grid both have six cells, and ncpl alone cannot tell them apart."""

    assert _grid(3, 2).ncpl == _grid(2, 3).ncpl == 6
    with pytest.raises(ValueError, match=r"same cell count but different geometry"):
        _build(tmp_path, "c", _grid(3, 2), _grid(2, 3))


def test_the_same_centres_with_a_different_vertex_count_is_refused(tmp_path):
    """Same cell count AND same centres is still not the same mesh.

    Two meshes can share their cell centres and differ in their boundaries -- the
    vertex count is what separates them, and the cell-centre check cannot see it. A
    mutant that dropped the ``nvert`` comparison survived every other test here.
    """

    grid = _grid(6, 6)
    props = dict(grid.get_disv_gridprops())
    # One extra, unreferenced vertex: every cell centre is untouched.
    props["vertices"] = [*props["vertices"], [len(props["vertices"]), 1.0, 1.0]]
    props["nvert"] = props["nvert"] + 1
    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **props),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="nv", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with pytest.raises(ValueError, match=r"nvert=\d+ but the context grid has \d+ vertices"):
        sim.build_flopy(workspace=tmp_path / "nv")


def test_a_stale_cellid_names_its_package(tmp_path):
    """The half-lazy case: a lazy disv beside a package holding the old grid's cells."""

    stale = mf.chd(stress_period_data={0: [[(0, 290), 9.0]]})
    with pytest.raises(ValueError, match=r"package 'chd' holds 1 cellid\(s\) this model's"):
        _build(tmp_path, "d", _grid(10, 10), _grid(10, 10), extra=(stale,))


def test_the_first_cell_past_the_end_is_stale(tmp_path):
    """The boundary: cellid == ncpl is the first invalid one, and ncpl - 1 is valid.

    Without this an off-by-one (``>`` for ``>=``) still catches a wildly stale cellid
    and passes every other test.
    """

    ok = mf.chd(stress_period_data={0: [[(0, 99), 9.0]]})
    assert _build(tmp_path, "b0", _grid(10, 10), _grid(10, 10), extra=(ok,)) is not None
    with pytest.raises(ValueError, match=r"\(0, 100\)"):
        _build(tmp_path, "b1", _grid(10, 10), _grid(10, 10),
               extra=(mf.chd(stress_period_data={0: [[(0, 100), 9.0]]}),))


def test_the_first_layer_past_the_end_is_stale(tmp_path):
    """Same boundary on the layer axis: layer == nlay is the first invalid one."""

    with pytest.raises(ValueError, match=r"\(1, 4\)"):
        _build(tmp_path, "b2", _grid(10, 10), _grid(10, 10),
               extra=(mf.chd(stress_period_data={0: [[(1, 4), 9.0]]}),))


@pytest.mark.parametrize("cellid", [(0, -1), (-1, 4)])
def test_a_negative_cellid_is_refused(tmp_path, cellid):
    """The LOWER bound: `(0, -1)` passed `>= ncpl`, and FloPy wrote it as cell 0 --
    one-based, so a cell that does not exist and not the one typed (ledger 382)."""

    bad = mf.chd(stress_period_data={0: [[cellid, 9.0]]})
    with pytest.raises(ValueError, match=r"package 'chd' holds 1 cellid\(s\) with a negative part"):
        _build(tmp_path, "neg", _grid(10, 10), _grid(10, 10), extra=(bad,))


def test_a_negative_cellid_is_refused_on_a_structured_model(tmp_path):
    """No DISV, so no `ncpl` range check -- the lower bound holds on any discretization."""

    model = mf.gwf(
        "m",
        packages=[
            mf.dis(nlay=1, nrow=1, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=[0.0]),
            mf.npf(k=1.0),
            mf.ic(strt=5.0),
            mf.chd(stress_period_data={0: [[(0, 0, -1), 9.0]]}),
        ],
    )
    with pytest.raises(ValueError, match=r"negative part: \(0, 0, -1\)"):
        mf.simulation(model, name="s").build_flopy(workspace=tmp_path / "dis")


def test_a_negative_node_is_refused_on_a_disu_model(tmp_path):
    """DISU's cellid is `(node,)`, one integer -- which the range checks' cellid reader
    skips, so the lower bound reads it on its own (ledger 382)."""

    disu = mf.disu(
        nodes=2, nja=4, top=10.0, bot=0.0, area=100.0, iac=[2, 2], ja=[0, 1, 1, 0],
        ihc=[1, 1, 1, 1], cl12=[0.0, 5.0, 0.0, 5.0], hwva=[0.0, 10.0, 0.0, 10.0],
    )
    model = mf.gwf(
        "m",
        packages=[
            disu, mf.npf(k=1.0), mf.ic(strt=5.0),
            mf.chd(stress_period_data={0: [[(-1,), 9.0]]}),
        ],
    )
    with pytest.raises(ValueError, match=r"negative part: \(-1,\)"):
        mf.simulation(model, name="s").build_flopy(workspace=tmp_path / "disu")


def test_a_stale_layer_is_caught_too(tmp_path):
    """A cellid can be stale in its LAYER without being stale in its cell."""

    stale = mf.chd(stress_period_data={0: [[(3, 4), 9.0]]})
    with pytest.raises(ValueError, match=r"nlay=1, ncpl=100"):
        _build(tmp_path, "e", _grid(10, 10), _grid(10, 10), extra=(stale,))


def test_a_refined_swap_is_the_documented_gap(tmp_path):
    """A stale cellid that still falls inside a LARGER new grid cannot be range-checked.

    This builds, deliberately. Closing it needs the grid pin -- a hash of the mesh a
    package's cell numbers refer to -- not a wider range check. Pinned so that making
    it raise later is a deliberate change with a ledger entry, not a silent one.
    """

    stale = mf.chd(stress_period_data={0: [[(0, 95), 9.0]]})
    built = _build(tmp_path, "f", _grid(20, 15), _grid(20, 15), extra=(stale,))
    assert built.built_model("m").model.get_package("chd") is not None


# ---------------------------------------------------------------- what must still build


def test_a_consistent_model_builds(tmp_path):
    """The control. The guard must be invisible to every correct model."""

    built = _build(tmp_path, "g", _grid(10, 10), _grid(10, 10))
    assert int(built.built_model("m").model.get_package("disv").ncpl.get_data()) == 100


def test_the_same_mesh_in_another_frame_builds(tmp_path):
    """A converted model's own frame is deliberately not the world's (ledgers 159, 218).

    ``to_mf6(local_origin=True)`` shifts to a local origin because MF6 loses DISV
    conductance precision on State Plane coordinates. The DISV is then written in the
    model's frame while the grid stays in the world's -- on the synthetic USG fixture
    that is 1.34e+06 length units apart, and comparing positions condemned it.
    """

    grid = _grid(6, 6)
    props = grid.get_disv_gridprops()
    shifted = dict(props)
    shifted["vertices"] = [[i, x + 1_340_000.0, y + 400_000.0] for i, x, y in props["vertices"]]
    shifted["cell2d"] = [
        [row[0], row[1] + 1_340_000.0, row[2] + 400_000.0, *row[3:]] for row in props["cell2d"]
    ]
    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **shifted),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="h", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    assert sim.build_flopy(workspace=tmp_path / "h") is not None


# ---------------------------------------------------------------- context shapes


def test_a_domain_shaped_for_another_grid_is_refused(tmp_path):
    """The stale piece ``replace_grid`` leaves behind."""

    with pytest.raises(ValueError, match=r"context domain is shaped .* against .* ncpl=100"):
        _build(tmp_path, "i", _grid(10, 10), _grid(10, 10), domain=np.ones((1, 37), dtype=int))


def test_a_domain_with_the_wrong_layer_count_is_refused(tmp_path):
    """Right number of cells, wrong number of layers."""

    with pytest.raises(ValueError, match=r"which is 4 layers, against .* nlay=1"):
        _build(tmp_path, "j", _grid(10, 10), _grid(10, 10), domain=np.ones((4, 100), dtype=int))


def test_a_surfaces_value_is_not_second_guessed(tmp_path):
    """``context.surfaces`` is deliberately NOT checked, and must not block a build.

    It is documented loosely -- "vor.gdf_topbtm or equivalent" -- and legitimately
    arrives as a frame, a ``LayerBuildResult``, or an array laid out (nlay, ncpl) or
    (nlay + 1, ncpl). There is no reliable "rows are cells" invariant to test, so an
    earlier version keyed on ``.shape[0]`` and refused a documented value with a
    message about a table. ``specs._layer_elevations`` is where surfaces are
    interpreted, and where a disagreement belongs.
    """

    import geopandas as gpd
    from shapely.geometry import Point

    odd = gpd.GeoDataFrame(
        {0: [10.0] * 12, 1: [0.0] * 12}, geometry=[Point(0, 0)] * 12, crs=_grid(10, 10).crs
    )
    assert _build(tmp_path, "sf", _grid(10, 10), _grid(10, 10), surfaces=odd) is not None
    assert _build(tmp_path, "sa", _grid(10, 10), _grid(10, 10),
                  surfaces=np.zeros((2, 100))) is not None


def test_every_disagreement_is_named_at_once(tmp_path):
    """One build, one error, listing all of it -- not one fix per rebuild."""

    stale = mf.chd(stress_period_data={0: [[(0, 290), 9.0]]})
    with pytest.raises(ValueError) as caught:
        _build(
            tmp_path, "l", _grid(10, 10), _grid(10, 10),
            extra=(stale,), domain=np.ones((1, 37), dtype=int),
        )
    message = str(caught.value)
    assert "context domain is shaped" in message
    assert "holds 1 cellid(s)" in message


# ---------------------------------------------------------------- the narrow reader


def test_only_a_leading_tuple_of_ints_reads_as_a_cellid():
    """``_cellids_in`` must not guess: a float or a bare int in position 0 is not a cellid."""

    assert list(_cellids_in({"stress_period_data": {0: [[(0, 4), 9.0]]}})) == [(0, 4)]
    assert list(_cellids_in({"stress_period_data": {0: [[0, 4, 9.0]]}})) == []
    assert list(_cellids_in({"stress_period_data": {0: [[(0, 4.5), 9.0]]}})) == []
    assert list(_cellids_in({"stress_period_data": {0: [[True, 9.0]]}})) == []
    # The bool guard is on the cellid's PARTS: True == 1, so an unguarded reader would
    # report (True, 4) as layer 1 and give a stale-cellid error for a type error.
    assert list(_cellids_in({"stress_period_data": {0: [[(True, 4), 9.0]]}})) == []
    assert list(_cellids_in({"stress_period_data": {0: []}})) == []
    assert list(_cellids_in({"packagedata": [[0, (0, 4), 1]]})) == []
    assert list(_cellids_in({})) == []


def test_a_three_part_structured_cellid_is_read():
    """A DIS cellid is (layer, row, col); the reader takes its last part as the cell."""

    assert list(_cellids_in({"stress_period_data": {0: [[(0, 1, 2), 9.0]]}})) == [(0, 1, 2)]


def test_a_grid_that_cannot_describe_itself_weakens_the_guard_rather_than_replacing_it():
    """A guard that cannot compare must not become the error the user sees."""

    class Hostile:
        ncpl = 4

        def get_disv_gridprops(self):
            raise RuntimeError("no triangulation on this grid")

    disv = SimpleNamespace(
        ncpl=SimpleNamespace(get_data=lambda: 4),
        nlay=SimpleNamespace(get_data=lambda: 1),
        nvert=SimpleNamespace(get_data=lambda: 9),
        cell2d=SimpleNamespace(get_data=lambda: [[0, 0.0, 0.0]]),
    )
    model = SimpleNamespace(get_package=lambda name: disv if name == "disv" else None)
    context = SimpleNamespace(grid=Hostile(), domain=None, surfaces=None)
    assert _discretization_complaints(model, context, []) == []


def test_a_model_without_a_disv_is_left_alone():
    """DIS and DISU carry no ncpl, so there is no range here to check (only the lower
    bound, which needs none -- `test_a_negative_cellid_is_refused_on_a_structured_model`)."""

    model = SimpleNamespace(get_package=lambda name: None)
    context = SimpleNamespace(grid=_grid(3, 3), domain=np.ones((1, 9), dtype=int), surfaces=None)
    assert _discretization_complaints(model, context, []) == []


# ------------------------------------------------- what the sorted-radii version missed


def _disv_from(props, grid, tag, tmp_path, extra=()):
    """Build a model whose DISV comes from hand-edited gridprops."""

    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **props),
            mf.npf(k=1.0), mf.ic(strt=5.0), *extra, mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name=tag, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    return sim.build_flopy(workspace=tmp_path / tag)


def test_a_renumbered_mesh_is_refused(tmp_path):
    """Cell ORDER is the whole meaning of a cellid, so a renumbering is staleness.

    The first version of this guard compared each point's distance from its own
    centroid, SORTED, so that a rotated frame would still pass. Sorting threw away
    exactly the thing a cellid means: a renumbered mesh passed while every package
    cellid pointed somewhere else -- measured on a 10x10 grid, cell 5 moved 865.7
    length units and the build was clean.
    """

    grid = _grid(10, 10)
    props = dict(grid.get_disv_gridprops())
    order = list(range(grid.ncpl))
    order = order[7:] + order[:7]                       # a rotation of the numbering
    props["cell2d"] = [[new, *props["cell2d"][old][1:]] for new, old in enumerate(order)]
    with pytest.raises(ValueError, match=r"A renumbering counts"):
        _disv_from(props, grid, "renum", tmp_path)


def test_a_mirrored_mesh_is_refused(tmp_path):
    """Mirroring preserves every distance from the centroid, and moves every cell."""

    grid = _grid(8, 8)
    props = dict(grid.get_disv_gridprops())
    props["vertices"] = [[i, -x, y] for i, x, y in props["vertices"]]
    props["cell2d"] = [[r[0], -r[1], r[2], *r[3:]] for r in props["cell2d"]]
    with pytest.raises(ValueError, match=r"same cell count but different geometry"):
        _disv_from(props, grid, "mirror", tmp_path)


def test_coordinates_written_to_a_hundredth_of_a_foot_still_build(tmp_path):
    """The same mesh, less precisely written, is the same mesh.

    An absolute tolerance refused this: the centres move by up to 0.005 length units,
    which is more than any fixed epsilon small enough to be meaningful elsewhere. The
    tolerance is a thousandth of a CELL now, taken from the data.
    """

    grid = _grid(8, 8)
    props = dict(grid.get_disv_gridprops())
    props["vertices"] = [[i, round(x, 2), round(y, 2)] for i, x, y in props["vertices"]]
    props["cell2d"] = [[r[0], round(r[1], 2), round(r[2], 2), *r[3:]] for r in props["cell2d"]]
    assert _disv_from(props, grid, "round", tmp_path) is not None


def test_a_float32_rendering_at_state_plane_still_builds(tmp_path):
    """float32 steps ~0.06 length units at 1.34e+06, which no fixed epsilon survives."""

    grid = _grid(8, 8)
    props = dict(grid.get_disv_gridprops())
    shift = 1_340_000.0
    props["vertices"] = [
        [i, float(np.float32(x + shift)), float(np.float32(y + shift))]
        for i, x, y in props["vertices"]
    ]
    props["cell2d"] = [
        [r[0], float(np.float32(r[1] + shift)), float(np.float32(r[2] + shift)), *r[3:]]
        for r in props["cell2d"]
    ]
    assert _disv_from(props, grid, "f32", tmp_path) is not None


# ------------------------------------------------- a guard must not become the error


def test_a_disv_that_leaves_nvert_to_flopy_still_builds(tmp_path):
    """FloPy computes NVERT at write time, so ``nvert.get_data()`` is None until then."""

    from myflopy.specs import _as_int

    assert _as_int(None) is None
    assert _as_int(np.array([1, 2])) is None      # an UnstructuredGrid's per-layer ncpl
    assert _as_int(7) == 7
    assert _as_int(np.int64(7)) == 7


def test_a_grid_whose_gridprops_omit_nvert_still_builds(tmp_path):
    """FloPy's own ``get_gridprops_vertexgrid`` deletes that key, so it is not exotic."""

    grid = _grid(6, 6)
    full = grid.get_disv_gridprops()

    class NoNvert:
        ncpl = grid.ncpl
        crs = grid.crs

        def get_disv_gridprops(self):
            trimmed = dict(full)
            trimmed.pop("nvert")
            return trimmed

    model = mf.gwf(
        "m",
        context=ModelContext(grid=NoNvert(), domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **full),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="nn", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    assert sim.build_flopy(workspace=tmp_path / "nn") is not None


def test_a_disabled_package_is_not_scanned(tmp_path):
    """A package that never reaches the model cannot make the model stale."""

    import dataclasses

    stale = dataclasses.replace(
        mf.chd(stress_period_data={0: [[(0, 290), 9.0]]}), enabled=False
    )
    assert _build(tmp_path, "off", _grid(10, 10), _grid(10, 10), extra=(stale,)) is not None


# ------------------------------------------------- the reader's reach


def test_a_bare_list_of_records_is_read_too(tmp_path):
    """``stress_period_data`` as a list, not a dict, used to bypass the check entirely."""

    assert list(_cellids_in({"stress_period_data": [[(0, 4), 9.0]]})) == [(0, 4)]
    with pytest.raises(ValueError, match=r"holds 1 cellid"):
        _build(tmp_path, "bare", _grid(10, 10), _grid(10, 10),
               extra=(mf.chd(stress_period_data=[[(0, 290), 9.0]]),))


def test_both_halves_of_a_barrier_pair_are_read():
    """An HFB record is cellid, cellid, hydchr -- and either half can go stale."""

    found = list(_cellids_in({"stress_period_data": {0: [[(0, 1), (0, 2), 1e-3]]}}))
    assert found == [(0, 1), (0, 2)]


def test_the_reader_stops_before_the_value_fields():
    """It walks LEADING cellids only, so a value that happens to be a pair is not one."""

    assert list(_cellids_in({"stress_period_data": {0: [[(0, 1), 9.0, (3, 4)]]}})) == [(0, 1)]
    assert list(_cellids_in({"stress_period_data": {0: [[9.0, (0, 1)]]}})) == []


def test_advanced_package_data_is_deliberately_not_read():
    """SFR's packagedata[1] is a cellid and LAK's is a stage, so neither is guessed at."""

    assert list(_cellids_in({"packagedata": [[0, (0, 4), 1]]})) == []
    assert list(_cellids_in({"connectiondata": [[0, (0, 4)]]})) == []


def test_a_flat_domain_is_accepted_for_any_layer_count(tmp_path):
    """``(ncpl,)`` means the same activity in every layer, which is a valid statement."""

    flat = np.ones(100, dtype=int)
    assert _build(tmp_path, "flat", _grid(10, 10), _grid(10, 10), domain=flat) is not None


def test_a_one_row_domain_on_a_multi_layer_model_points_at_the_flat_form(tmp_path):
    """``(1, ncpl)`` on a 3-layer model cannot answer for layers 1 and 2.

    The flat ``(ncpl,)`` spelling of the same intent IS accepted, which looks
    inconsistent until you read ``domain_active``: a two-dimensional domain is
    indexed per layer and raises for any layer past its rows, so a one-row domain on
    a three-layer model is broken for every consumer. Refusing it here reports it
    once, by name, instead of deferring it to whichever package asks first.
    """

    grid = _grid(10, 10)
    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=3, top=30.0, botm=[20.0, 10.0, 0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="1row", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with pytest.raises(ValueError, match=r"give it as a flat \(ncpl,\) array"):
        sim.build_flopy(workspace=tmp_path / "1row")


def test_many_stale_cellids_are_capped_and_counted(tmp_path):
    """Five are listed, the rest counted, and duplicates collapse to one."""

    records = [[(0, 200 + i), 9.0] for i in range(8)] + [[(0, 200), 9.0]] * 3
    with pytest.raises(ValueError, match=r"holds 8 cellid\(s\).*and 3 more"):
        _build(tmp_path, "many", _grid(10, 10), _grid(10, 10),
               extra=(mf.chd(stress_period_data={0: records}),))


# ------------------------------------------- what each refusal says it is about


def _dis_with_a_negative_cellid():
    return mf.gwf(
        "m",
        packages=[
            mf.dis(nlay=1, nrow=1, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=[0.0]),
            mf.npf(k=1.0),
            mf.ic(strt=5.0),
            mf.chd(stress_period_data={0: [[(0, 0, -1), 9.0]]}),
        ],
    )


@pytest.mark.parametrize(
    ("case", "heading"),
    [
        ("mesh", "does not agree with the grid it was built against"),
        ("domain", "has a context domain that does not fit the DISV it was built with"),
        ("past-the-grid", "names cells its grid does not have"),
        ("negative-on-dis", "names cells its grid does not have"),
    ],
)
def test_each_kind_of_build_refusal_says_what_it_is_about(tmp_path, case, heading):
    """One preamble, about cell numbers resolved on another mesh, used to head every
    refusal -- a negative cellid on a DIS model, with no mesh and no move, was told to
    re-declare its packages "against this grid" (ledger 427)."""

    with pytest.raises(ValueError) as raised:
        if case == "mesh":
            # A domain shaped for the DISV, so the mesh is the only complaint.
            _build(tmp_path, "a", _grid(10, 10), _grid(20, 15), domain=np.ones((1, 100), dtype=int))
        elif case == "domain":
            _build(tmp_path, "d", _grid(10, 10), _grid(10, 10), domain=np.ones((1, 37), dtype=int))
        elif case == "past-the-grid":
            _build(tmp_path, "p", _grid(10, 10), _grid(10, 10),
                   extra=(mf.chd(stress_period_data={0: [[(0, 100), 9.0]]}),))
        else:
            mf.simulation(_dis_with_a_negative_cellid(), name="s").build_flopy(
                workspace=tmp_path / "dis"
            )
    text = str(raised.value)
    assert text.startswith(f"Model 'm' {heading}:\n  - "), text
    assert text.count("Model 'm'") == 1
    assert ("resolved on another mesh" in text) == (case == "mesh"), text


def test_a_refusal_of_several_kinds_heads_each_group_in_turn(tmp_path):
    """Two kinds of complaint, two paragraphs, each with its own heading and closing."""

    with pytest.raises(ValueError) as raised:
        _build(tmp_path, "two", _grid(10, 10), _grid(10, 10),
               domain=np.ones((1, 37), dtype=int),
               extra=(mf.chd(stress_period_data={0: [[(0, -1), 9.0]]}),))
    text = str(raised.value)
    domain = text.index("Model 'm' has a context domain that does not fit")
    cellid = text.index("Model 'm' names cells its grid does not have:\n  - package 'chd' holds 1")
    assert domain == 0 and cellid > domain
    assert "Correct the context's `domain=`." in text[:cellid]
    assert text.rstrip().endswith("re-declare the package against this grid.")
