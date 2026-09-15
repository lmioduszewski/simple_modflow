"""A resolved package records the mesh its cell numbers refer to, and build checks it.

`_discretization_complaints` (ledger 249) range-checks cellids, which cannot see a
package resolved on a SMALL grid and rebuilt on a LARGER one: every number is still
in range, and every one names a different cell. The pin is what closes that -- it is
the mesh's identity, not a bound on its size.

The fingerprint is CONNECTIVITY only -- cell count, vertex count, and each cell's
vertex-index list in cell order. Nothing continuous, so no tolerance: the same mesh
written to fewer decimals, rendered float32, or shifted to a local origin pins the
same, while a renumbering does not.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import GRID_PIN, ModelContext, grid_pin


def _grid(nrow: int, ncol: int):
    """A rectangular Voronoi grid of ``nrow x ncol`` cells."""

    return rectangular_voronoi(CanonicalModelConfig(nrow=nrow, ncol=ncol, nlay=1, nper=1))


@pytest.fixture
def source(tmp_path):
    """A one-feature CHD GeoPackage covering the western edge of a 6x6 grid."""

    grid = _grid(6, 6)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [9.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)],
        crs=grid.crs,
    ).to_file(path, layer="chd", driver="GPKG")
    return path


def _model(grid, tag, tmp_path, extra=()):
    """Build a minimal model on ``grid``."""

    spec = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), *extra, mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(spec, name=tag, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    return sim.build_flopy(workspace=tmp_path / tag)


# ---------------------------------------------------------------- the fingerprint


def test_a_mesh_pins_the_same_every_time():
    """A fingerprint that moved between calls would refuse every model."""

    grid = _grid(5, 5)
    assert grid_pin(grid) == grid_pin(grid)


def test_two_different_meshes_pin_differently():
    assert grid_pin(_grid(5, 5)) != grid_pin(_grid(6, 6))


def test_the_same_mesh_built_twice_pins_the_same():
    """Two instances of one mesh are one mesh; the pin is content, not identity."""

    assert grid_pin(_grid(5, 5)) == grid_pin(_grid(5, 5))


def test_a_local_origin_shift_does_not_change_the_pin():
    """The USG importer shifts to a local origin (ledgers 159, 218) and renumbers nothing."""

    grid = _grid(5, 5)
    before = grid_pin(grid)
    grid.verts = np.asarray(grid.verts, dtype=float) + np.array([1_340_000.0, 400_000.0])
    assert grid_pin(grid) == before


def test_a_renumbered_mesh_pins_differently():
    """Cell order is the whole meaning of a cellid, so it is the whole content of a pin."""

    grid = _grid(5, 5)
    before = grid_pin(grid)
    order = list(range(grid.ncpl))
    grid.iverts = [grid.iverts[i] for i in (order[3:] + order[:3])]
    assert grid_pin(grid) != before


def test_an_ungrindable_grid_is_unpinned_rather_than_an_error():
    """A pin that cannot be taken means unpinned -- never an error at the caller."""

    class Hostile:
        def get_disv_gridprops(self):
            raise RuntimeError("no triangulation here")

    assert grid_pin(Hostile()) is None
    assert grid_pin(None) is None


# ---------------------------------------------------------------- the pin in use


def test_a_gpkg_package_records_its_grid(source):
    """Every list BC reaches one provenance method, so pinning there pins all seven."""

    grid = _grid(6, 6)
    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
    )
    assert chd.metadata[GRID_PIN] == grid_pin(grid)


def test_a_package_resolved_on_a_smaller_grid_is_refused(tmp_path, source):
    """The gap ledger 249 left open: every cellid in range, every one the wrong cell."""

    small = _grid(6, 6)
    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=small, domain=np.ones((1, small.ncpl), dtype=int)),
    )
    with pytest.raises(ValueError, match=r"was resolved against a different mesh"):
        _model(_grid(12, 12), "bigger", tmp_path, extra=(chd,))


def test_the_refusal_names_both_meshes(tmp_path, source):
    """Two hashes in the message, so the reader can tell which package moved."""

    small = _grid(6, 6)
    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=small, domain=np.ones((1, small.ncpl), dtype=int)),
    )
    with pytest.raises(ValueError) as caught:
        _model(_grid(12, 12), "named", tmp_path, extra=(chd,))
    message = str(caught.value)
    assert grid_pin(small) in message
    assert grid_pin(_grid(12, 12)) in message


def test_a_package_built_on_its_own_grid_is_untouched(tmp_path, source):
    """The control: pinning must be invisible to every model that is consistent."""

    grid = _grid(6, 6)
    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
    )
    assert _model(grid, "same", tmp_path, extra=(chd,)) is not None


def test_a_package_built_on_an_equal_mesh_is_untouched(tmp_path, source):
    """A rebuilt grid object of the same mesh is the same mesh."""

    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=_grid(6, 6), domain=np.ones((1, 36), dtype=int)),
    )
    assert _model(_grid(6, 6), "equal", tmp_path, extra=(chd,)) is not None


def test_a_hand_written_package_stays_unpinned(tmp_path):
    """You cannot pin what you were not told: hand-typed cellids name no source mesh.

    This is why the range check in ledger 249 is still needed, and why its documented
    gap stays documented -- a pin can only speak for records the library resolved.
    """

    chd = mf.chd(stress_period_data={0: [[(0, 5), 9.0]]})
    assert GRID_PIN not in (chd.metadata or {})
    assert _model(_grid(12, 12), "hand", tmp_path, extra=(chd,)) is not None


# ---------------------------------------------------------------- what the pin is not


def test_the_pin_is_not_in_the_digest(source):
    """Digesting it would report every stored run as drifted for no change in output."""

    from myflopy.specs import package_digest

    grid = _grid(6, 6)
    ctx = ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int))
    pinned = mf.chd.gpkg(source, layer="chd", nper=1, head="head", context=ctx)
    unpinned = pinned.with_metadata(**{GRID_PIN: "something else entirely"})
    assert package_digest(pinned) == package_digest(unpinned)


def test_the_pin_survives_a_round_trip(source):
    """It rides in metadata, which ``to_dict``/``from_dict`` carry."""

    from myflopy.specs import PackageSpec

    grid = _grid(6, 6)
    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
    )
    assert PackageSpec.from_dict(chd.to_dict()).metadata[GRID_PIN] == chd.metadata[GRID_PIN]


# ------------------------------------------------- the advanced builders (ledger 264)


def test_the_areal_builders_pin_their_grid(tmp_path):
    """`mf.rch(context=)` resolves per-cell records through RCHBuilder, not a .gpkg.

    Measured before this: built on a 36-cell grid, recharge landed on the wrong 16
    cells and MF6 ran to completion.
    """

    grid = _grid(4, 4)
    rch = mf.rch(
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        nper=1, recharge=1e-3,
    )
    assert rch.metadata[GRID_PIN] == grid_pin(grid)


def test_a_package_from_an_advanced_builder_is_refused_on_another_mesh(tmp_path):
    small, big = _grid(4, 4), _grid(6, 6)
    rch = mf.rch(
        context=ModelContext(grid=small, domain=np.ones((1, small.ncpl), dtype=int)),
        nper=1, recharge=1e-3,
    )
    with pytest.raises(ValueError, match=r"package 'rch' was resolved against a different mesh"):
        _model(big, "wrong", tmp_path, extra=(rch,))


def test_pin_metadata_is_empty_without_a_pinnable_grid():
    """A builder with no context, or a context with no grid, stamps nothing."""

    from myflopy.grid_identity import pin_metadata

    assert pin_metadata(None) == {}
    assert pin_metadata(ModelContext()) == {}
    grid = _grid(3, 3)
    assert pin_metadata(ModelContext(grid=grid)) == {GRID_PIN: grid_pin(grid)}


def test_a_model_with_no_grid_at_all_still_compares_package_pins(tmp_path, source):
    """The guard's LIVE hole: no grid object and no recorded pin.

    `mf.gwf(name, packages=[...])` with no context is used by this repo's own
    examples. The built DISV can answer for the mesh, which is what gives the package
    pins something to compare against -- measured before: BUILT CLEAN.
    """

    small = _grid(6, 6)
    chd = mf.chd.gpkg(
        source, layer="chd", nper=1, head="head",
        context=ModelContext(grid=small, domain=np.ones((1, small.ncpl), dtype=int)),
    )
    big = _grid(12, 12)
    model = mf.gwf(
        "m",
        context=ModelContext(domain=np.ones((1, big.ncpl), dtype=int)),   # no grid
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **big.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), chd, mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="nogrid",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with pytest.raises(ValueError, match=r"was resolved against a different mesh"):
        sim.build_flopy(tmp_path / "ws")


def test_the_recorded_pin_is_consulted_before_the_disv_s_own(tmp_path):
    """Order matters: taking the DISV's fingerprint first makes it agree with itself.

    That mistake was made and caught here -- it silently disabled the reload check
    added in ledger 257 while every other test still passed.
    """

    from myflopy.grid_identity import _discretization_complaints

    built = _model(_grid(6, 6), "ord", tmp_path).built_model("m").model
    stale = ModelContext(metadata={GRID_PIN: "0000deadbeef0000"})
    complaints = _discretization_complaints(built, stale, [])
    assert any("is not the mesh this model was declared against" in c for c in complaints)


# --- every grid-bound helper, on two REAL meshes ---------------------------- #
#
# Ledger 250 stamped the seven `.gpkg` forms and ledger 264 the areal LIST forms,
# and both reported the gap closed. It was not: `lak`, `hfb` and the ARRAY forms
# carried no pin, and the first two build clean on the wrong mesh because the range
# check cannot see their cellids at all. The coverage test below is the part meant
# to keep this from happening a third time -- it BUILDS each helper and reads its
# metadata, which is how the `rch.array` gap was finally found, rather than reading
# the source and believing it (ledger 267).


def _two_layer_context(grid):
    """A context with the surfaces the advanced builders need."""

    stack = mf.LayerStack(grid, top=mf.Surface.flat(10.0))
    stack.add("upper", thickness=5.0)
    stack.add("lower", thickness=5.0)
    result = stack.build(grid, attach=True)
    return ModelContext(grid=grid, domain=result.idomain, surfaces=grid.gdf_topbtm)


def _lake_source(tmp_path, grid):
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "lake.gpkg"
    mid_x, mid_y = (xmin + xmax) / 2, (ymin + ymax) / 2
    span = (xmax - xmin) * 0.15
    gpd.GeoDataFrame(
        {"name": ["L"], "strt": [5.0]},
        geometry=[box(mid_x - span, mid_y - span, mid_x + span, mid_y + span)],
        crs=grid.crs,
    ).to_file(path, layer="lake", driver="GPKG")
    return path


def _grid_bound_helpers(tmp_path, grid):
    """Every helper that resolves cell numbers from a context, built on ``grid``."""

    from shapely.geometry import LineString

    context = _two_layer_context(grid)
    xmin, _ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    mid = (xmin + xmax) / 2
    return {
        "lak": mf.lak(
            lakes=_lake_source(tmp_path, grid), context=context,
            starting_stage=5.0, lake_bottom=0.0, nper=1,
        ),
        "hfb.line": mf.hfb.line(
            LineString([(mid, _ymin), (mid, ymax)]), hydchr=1e-5, context=context
        ),
        "rch.array": mf.rch.array(recharge=0.001, context=context),
        "evt.array": mf.evt.array(
            surface=10.0, rate=0.001, depth=2.0, context=context
        ),
    }


def test_every_grid_bound_helper_pins_the_mesh_it_resolved_on(tmp_path, triangulated_grid):
    """Both directions: a new unpinned helper fails, and so does a stale entry here.

    Reading the source is what missed `rch.array` twice -- it reaches
    `advanced.py::_array_bc_spec`, not the `areal.py` funnel the ledger named. This
    builds each one and asks the spec.
    """

    built = _grid_bound_helpers(tmp_path, triangulated_grid)
    pinned = {name for name, spec in built.items() if GRID_PIN in (spec.metadata or {})}

    assert pinned == set(built), f"grid-bound helpers with no pin: {sorted(set(built) - pinned)}"
    for name, spec in built.items():
        assert spec.metadata[GRID_PIN] == grid_pin(triangulated_grid), name


def test_a_grid_free_array_package_is_not_pinned(tmp_path, triangulated_grid):
    """`irch=0` means the same layer on every mesh, so pinning it would refuse a move.

    The pin is decided on the resolved PAYLOAD -- does it name cells -- not on
    whether a `context=` was passed. Pinning on the latter is the easy mistake and
    it makes a correct move impossible.
    """

    context = _two_layer_context(triangulated_grid)

    grid_free = mf.rch.array(recharge=0.001, irch=0, context=context)
    per_cell = mf.rch.array(recharge=np.full(triangulated_grid.ncpl, 0.001), irch=None,
                            context=context)

    assert GRID_PIN not in (grid_free.metadata or {})
    assert GRID_PIN in (per_cell.metadata or {}), "a per-column array names cells"


def _refuses_on(tmp_path, spec, grid, tag):
    """Build ``spec`` against ``grid`` and return the refusal, if any."""

    model = mf.gwf(
        "m",
        context=_two_layer_context(grid),
        packages=[
            mf.disv(nlay=2, top=10.0, botm=[5.0, 0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), spec, mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(
        model, name="moved", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )
    sim.build_flopy(tmp_path / f"ws_{tag}")


@pytest.mark.parametrize("helper", ["lak", "hfb.line"])
def test_a_helper_the_range_check_cannot_see_is_refused_on_another_mesh(
    tmp_path, triangulated_grid, triangulated_grid_finer, helper
):
    """These two are the ones that BUILT CLEAN, and only the pin catches them.

    `_cellids_in` refuses to walk `packagedata`/`connectiondata` -- a cellid sits at a
    different index per package, so guessing would refuse correct models -- and a
    barrier's records are cell PAIRS with no cellid column at all. Measured before the
    pin: both built clean across 28 -> 92, on ragged real meshes.
    """

    spec = _grid_bound_helpers(tmp_path, triangulated_grid)[helper]

    with pytest.raises(ValueError, match=r"was resolved against a different mesh"):
        _refuses_on(tmp_path, spec, triangulated_grid_finer, helper.replace(".", "_"))


@pytest.mark.parametrize("helper", ["rch.array", "evt.array"])
def test_an_array_helper_is_refused_on_a_different_mesh_of_the_SAME_size(tmp_path, helper):
    """Same `ncpl`, different mesh -- the case that was silent, and the pin's whole job.

    A DIFFERENT `ncpl` already raised before the pin, but as FloPy's
    `MFDataException ... Expected data dimensions`, which fires while the package is
    still being built and so beats the guard to it. That is loud and unhelpful rather
    than silent, and improving it is not this commit's subject. The SAME-`ncpl` swap
    below said nothing at all and applied the wrong per-cell layer indices.

    `rectangular_voronoi` is used deliberately: these two must differ in GEOMETRY while
    agreeing in size, which two `GridSpec.voronoi` grids cannot be made to do. The pin
    is connectivity-only, so the vertex dtype that matters elsewhere (ledger 263) does
    not matter here.
    """

    grid_a = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=4, nlay=2, nper=1))
    grid_b = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=3, nlay=2, nper=1))
    assert grid_a.ncpl == grid_b.ncpl and grid_pin(grid_a) != grid_pin(grid_b)

    spec = _grid_bound_helpers(tmp_path, grid_a)[helper]

    with pytest.raises(ValueError, match=r"was resolved against a different mesh"):
        _refuses_on(tmp_path, spec, grid_b, helper.replace(".", "_"))


# --- two models, one mesh --------------------------------------------------- #
#
# `_discretization_complaints` validates each model against its OWN context and never
# compares siblings. MF6 does compare them -- `exg-gwfgwt.f90:224` terminates on a
# node-count mismatch and `:239-255` on an idomain mismatch -- but it compares COUNTS
# and IDOMAIN only, never geometry, so a mesh of the same size passes both of its
# checks and it maps node-for-node onto the wrong one (ledger 268).


def _flow(name, grid):
    return mf.gwf(
        name,
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord=f"{name}.hds"),
        ],
    )


def _transport(name, grid):
    return mf.gwt(
        name,
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.mst(porosity=0.3), mf.ic(strt=0.0), mf.adv(),
            mf.oc(concentration_filerecord=f"{name}.ucn"),
        ],
    )


def _coupled(flow_grid, transport_grid):
    flow, transport = _flow("flow", flow_grid), _transport("tran", transport_grid)
    return mf.simulation(
        flow, transport, name="c",
        packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=["flow"]), mf.ims(models=["tran"], name="ims_t"),
        ],
        exchanges=[mf.exchange(flow, transport)],
    )


def test_coupled_models_on_different_meshes_of_the_same_size_are_refused(tmp_path):
    """The measured failure: `ncpl=6` both sides, built clean, ran, zero warnings."""

    flow_grid = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=2, nlay=1, nper=1))
    transport_grid = rectangular_voronoi(CanonicalModelConfig(nrow=2, ncol=3, nlay=1, nper=1))
    assert flow_grid.ncpl == transport_grid.ncpl, "the same size is the whole point"

    with pytest.raises(ValueError, match=r"joins models on different meshes"):
        _coupled(flow_grid, transport_grid).build_flopy(tmp_path / "ws")


def test_the_refusal_names_the_exchange_and_both_meshes(tmp_path):
    """You cannot act on "something is wrong" -- it has to say which model moved."""

    flow_grid = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=2, nlay=1, nper=1))
    transport_grid = rectangular_voronoi(CanonicalModelConfig(nrow=2, ncol=3, nlay=1, nper=1))

    with pytest.raises(ValueError) as raised:
        _coupled(flow_grid, transport_grid).build_flopy(tmp_path / "ws")

    message = str(raised.value)
    assert "flow" in message and "tran" in message
    assert grid_pin(flow_grid) in message and grid_pin(transport_grid) in message


def test_coupled_models_on_one_mesh_build(tmp_path):
    """The ordinary case, which must stay ordinary."""

    grid = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=2, nlay=1, nper=1))

    assert _coupled(grid, grid).build_flopy(tmp_path / "ws") is not None


def test_a_gwfgwf_pair_on_two_meshes_still_builds(tmp_path):
    """Coupling a coarse model to a fine one is what that exchange is FOR.

    Excluding `gwfgwf` is the point of scoping the check to the three exchanges MF6
    itself requires to share a discretization.
    """

    coarse = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=3, nlay=1, nper=1))
    fine = rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1))
    sim = mf.simulation(
        _flow("coarse", coarse), _flow("fine", fine), name="two",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims()],
        exchanges=[
            mf.exchange.gwfgwf(
                "coarse", "fine", nexg=1,
                exchangedata=[[(0, 0), (0, 0), 1, 25.0, 25.0, 50.0]],
            )
        ],
    )

    assert sim.build_flopy(tmp_path / "ws") is not None
