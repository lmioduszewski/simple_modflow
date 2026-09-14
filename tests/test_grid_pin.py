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

import numpy as np
import pytest
import geopandas as gpd
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
