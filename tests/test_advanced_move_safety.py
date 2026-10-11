"""The advanced packages must not move silently (ledgers 302, 303).

Two holes, both silent-wrong-model rather than wrong-error:

* a GRIDLESS ``mf.uzf`` carried no pin and no cell marker, so ``package_move`` answered
  "unchanged -- names no cells" and a move BUILT CLEAN with the old mesh's UZF objects
  sitting on a prefix of the new mesh's columns;
* ``mf.mvr`` names its endpoints as bare integers and nothing compared them against the
  providers, so a mover beside a moved SFR routed water to a reach that no longer means
  what it meant.

Both were invisible in the same way: the verdict said the move was fine.
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.grid_identity import CELLWISE_OPTIONS, GRID_PIN, package_move
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext


def _grid(n, nlay=3):
    return rectangular_voronoi(
        CanonicalModelConfig(nrow=n, ncol=n, nlay=nlay, nper=1)
    )


@pytest.fixture(scope="module")
def meshes():
    return _grid(5), _grid(9)


def _uzf(context):
    return mf.uzf(
        context=context, nper=1, finf=1e-3, pet=1e-3,
        vks=1.0, thtr=0.05, thts=0.35, thti=0.08, eps=4.0,
    )


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    for name in ("a", "b", "c"):
        stack.add(name, thickness=20.0)
    return stack


# ------------------------------------------------------------------ C1: gridless uzf


def test_a_gridless_uzf_is_marked_as_addressing_cells(meshes):
    """It has no mesh to pin to, so it gets the weaker marker rather than nothing."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    spec = _uzf(ModelContext(domain=domain))

    assert GRID_PIN not in spec.metadata, "there was no grid to pin to"
    assert spec.metadata[CELLWISE_OPTIONS] == ["packagedata"]
    how, why = package_move(spec)
    assert how == "unverifiable", why


def test_a_grid_bearing_uzf_still_pins(meshes):
    """The eager, documented declaration must not move -- it was already safe."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    spec = _uzf(ModelContext(grid=small, domain=domain))

    assert GRID_PIN in spec.metadata
    assert CELLWISE_OPTIONS not in spec.metadata
    assert package_move(spec)[0] == "stale"


def test_the_marker_does_not_change_what_was_built(meshes):
    """A stamp is provenance; it must not touch the package's own records."""

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    gridless = _uzf(ModelContext(domain=domain))
    pinned = _uzf(ModelContext(grid=small, domain=domain))
    assert gridless.options["packagedata"] == pinned.options["packagedata"]
    assert gridless.options["perioddata"] == pinned.options["perioddata"]


def test_the_marker_moves_no_digest(meshes):
    """`package_digest` excludes metadata by contract -- assert it rather than trust it.

    Without this, closing the hole would report every stored run as edited on the day it
    landed, which is ledger 293's churn.
    """

    from dataclasses import replace

    from myflopy.specs import package_digest

    small, _ = meshes
    spec = _uzf(ModelContext(domain=np.ones((3, small.ncpl), dtype=int)))
    assert package_digest(spec) == package_digest(replace(spec, metadata={}))


def test_moving_a_gridless_uzf_is_now_REFUSED_at_build(meshes, tmp_path):
    """The defect itself: this used to build clean on a prefix of the new mesh.

    Measured before the fix: 75 UZF objects on columns 0..24 of an 81-cell mesh, zero
    warnings -- UZF over 31% of the domain.
    """

    small, large = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small, domain=domain),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                         _uzf(ModelContext(domain=domain))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    moved = sim.with_grid(large)
    assert next(r["how"] for r in moved.move_report() if r["package"] == "uzf") == "unverifiable"
    with pytest.raises(ValueError):
        moved.build(tmp_path / "moved")


def test_a_gridless_uzf_still_builds_on_its_OWN_mesh(meshes, tmp_path):
    """The fix must not narrow the legitimate case.

    A gridless declaration is correct for a model that never moves, and
    `_lazy_resolution_context`'s own docstring records that refusing one would make the
    lazy form narrower than the eager one (ledger 284). Two tests in
    `test_uzf_builder.py` declare exactly this shape.
    """

    small, _ = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small, domain=domain),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                         _uzf(ModelContext(domain=domain))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    built = sim.build(tmp_path / "own")
    uzf = built.simulation.get_model("m").get_package("uzf")
    assert len(uzf.packagedata.get_data()) == len(sim.model("m").packages[-1].options["packagedata"])


def test_the_limit_is_the_MOVE_not_every_wrong_grid(meshes, tmp_path):
    """What the marker does NOT close, pinned so the commit cannot overclaim.

    The stamp is read by `package_move`, which only `with_grid` and `move_report`
    consult. A model ASSEMBLED BY HAND against a different grid never goes through
    either, and still builds -- measured. Closing that is the build-time consistency
    guard's job (ledger 249), and it is a range check, not a proof.
    """

    small, large = meshes
    domain = np.ones((3, small.ncpl), dtype=int)
    hand_made = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=large, domain=np.ones((3, large.ncpl), int)),
               packages=[mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=50.0),
                         _uzf(ModelContext(domain=domain))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    built = hand_made.build(tmp_path / "byhand")
    rows = built.simulation.get_model("m").get_package("uzf").packagedata.get_data()
    columns = {int(r[1][1]) for r in rows}
    assert max(columns) < large.ncpl - 1, (
        "this test exists because the small mesh's cells are a PREFIX of the large "
        "one's, and nothing on this path objects"
    )


# ------------------------------------------------------------------------- C2: the mover


@pytest.fixture(scope="module")
def surface_water(tmp_path_factory):
    """A stream centreline and a lake polygon, both spanning the shared footprint."""

    import geopandas as gpd
    from shapely.geometry import LineString, box

    small, _ = _grid(6, nlay=1), None
    folder = tmp_path_factory.mktemp("sw")
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    streams = folder / "streams.gpkg"
    gpd.GeoDataFrame(
        {"name": ["main"], "rwid": [5.0], "rhk": [1.0], "rgrd": [0.001], "rbth": [1.0]},
        geometry=[LineString([(xmin + 1, (ymin + ymax) / 2), (xmax - 1, (ymin + ymax) / 2)])],
        crs=small.crs,
    ).to_file(streams, layer="s", driver="GPKG")
    lakes = folder / "lakes.gpkg"
    gpd.GeoDataFrame(
        {"name": ["pond"], "strt": [50.0], "bl": [0.1]},
        geometry=[box(xmin + 200, ymin + 200, xmin + 400, ymin + 400)],
        crs=small.crs,
    ).to_file(lakes, layer="l", driver="GPKG")
    return streams, lakes


def _flat_meshes():
    """Two single-layer meshes over the SAME footprint, so a centreline spans both."""

    return (
        rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1, cell_size=100.0)),
        rectangular_voronoi(CanonicalModelConfig(nrow=12, ncol=12, nlay=1, nper=1, cell_size=50.0)),
    )


def _surface_water_model(vor, surface_water):
    import pandas as pd

    streams, lakes = surface_water
    surfaces = pd.DataFrame({0: np.full(vor.ncpl, 60.0), 1: np.full(vor.ncpl, 0.0)})
    context = ModelContext(
        grid=vor, domain=np.ones((1, vor.ncpl), dtype=int), surfaces=surfaces
    )
    stack = mf.LayerStack(top=mf.Surface.flat(60.0))
    stack.add("a", thickness=60.0)
    sfr = mf.sfr(
        streams=streams, stream_id="name", context=context, nper=1, width="rwid",
        streambed_k="rhk", gradient="rgrd", streambed_thickness="rbth",
        reach_top=50.0, mover=True,
    )
    lak = mf.lak(
        lakes=lakes, lake_id_field="name", context=context, nper=1,
        starting_stage="strt", bed_leakance="bl", lake_bottom=40.0, mover=True,
    )
    mvr = mf.mvr(
        nper=1,
        moves=[mf.Move(mf.sfr_connection(sfr, "main"), mf.lak_connection(lak, "pond"))],
    )
    return mf.simulation(
        mf.gwf("m", context=context,
               packages=[mf.disv(stack), mf.npf(k=1.0), mf.ic(strt=50.0), sfr, lak, mvr]),
        name="s", packages=[mf.tdis(nper=1)],
    ), sfr, mvr


def test_a_mover_names_the_packages_it_indexes(surface_water):
    """The stamp both entry points must carry -- `mvr_spec` is the shared funnel."""

    from myflopy.grid_identity import INDEXED_PACKAGES

    small, _ = _flat_meshes()
    _sim, _sfr, mvr = _surface_water_model(small, surface_water)
    assert mvr.metadata[INDEXED_PACKAGES] == ["lak", "sfr"]


def test_the_raw_flopy_mover_carries_the_stamp_too(surface_water):
    """`mf.mvr.flopy(...)` reaches `mvr_spec` directly and would otherwise carry nothing."""

    from myflopy.grid_identity import INDEXED_PACKAGES

    spec = mf.mvr.flopy(
        packages=[["sfr"], ["lak"]],
        perioddata={0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]},
    )
    assert spec.metadata[INDEXED_PACKAGES] == ["lak", "sfr"]


def test_the_mover_index_really_does_move_with_the_mesh(surface_water):
    """The premise. Without this the verdict below is guarding nothing."""

    small, large = _flat_meshes()
    _s, sfr_a, mvr_a = _surface_water_model(small, surface_water)
    _l, sfr_b, mvr_b = _surface_water_model(large, surface_water)

    assert len(sfr_a.options["packagedata"]) != len(sfr_b.options["packagedata"])
    index_a = mvr_a.options["perioddata"][0][0][1]
    index_b = mvr_b.options["perioddata"][0][0][1]
    assert index_a != index_b, "the fixture must put the outlet at a mesh-dependent reach"


def test_a_mover_beside_a_moved_provider_is_no_longer_unchanged(surface_water):
    """The defect: it read `unchanged` -- "names no cells" -- while holding a stale reach.

    Measured before the fix: a mover pointing at reach 9 of 10, on a mesh whose same
    declaration resolves reach 20 of 21, reported as a clean move.
    """

    small, _ = _flat_meshes()
    sim, _sfr, _mvr = _surface_water_model(small, surface_water)

    rows = {row["package"]: row for row in sim.move_report()}
    assert rows["sfr"]["how"] == "stale"
    assert rows["mvr"]["how"] == "unverifiable"
    assert "by INDEX" in rows["mvr"]["why"]
    assert "sfr" in rows["mvr"]["why"]


def test_a_mover_whose_providers_all_stay_put_is_left_alone(surface_water):
    """The fix must not condemn every mover -- only one whose providers moved."""

    from myflopy.grid_identity import package_moves

    spec = mf.mvr.flopy(
        packages=[["sfr"]], perioddata={0: [["sfr", 0, "sfr", 1, "FACTOR", 1.0]]}
    )
    quiet = mf.npf(k=10.0)
    verdicts = package_moves([mf.npf(k=1.0, name="sfr"), spec, quiet])
    assert verdicts["sfr"][0] == "unchanged"
    assert verdicts["mvr"][0] == "unchanged", verdicts["mvr"][1]


def test_grid_free_exempts_a_mover(surface_water):
    """`grid_free=` is the caller asserting mesh-independence; it must reach here too."""

    from myflopy.grid_identity import package_moves

    small, _ = _flat_meshes()
    sim, _sfr, mvr = _surface_water_model(small, surface_water)
    packages = tuple(sim.model("m").packages)
    assert package_moves(packages, grid_free=("mvr",))["mvr"][0] == "unchanged"


def test_all_three_callers_agree_about_the_mover(surface_water, tmp_path):
    """Two of three would leave the run table and the report contradicting each other."""

    small, large = _flat_meshes()
    sim, _sfr, _mvr = _surface_water_model(small, surface_water)
    moved = sim.with_grid(large)

    report = {row["package"]: row["how"] for row in moved.move_report()}
    assert report["mvr"] == "unverifiable"

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_simulation(sim)
    run = project.prepare_run("r", sim)
    resolution = {row["package"]: row["how"] for row in run.resolution}
    assert resolution["mvr"] == "unverifiable", resolution


def test_the_stamp_moves_no_digest(surface_water):
    from dataclasses import replace

    from myflopy.specs import package_digest

    spec = mf.mvr.flopy(
        packages=[["sfr"]], perioddata={0: [["sfr", 0, "sfr", 1, "FACTOR", 1.0]]}
    )
    assert package_digest(spec) == package_digest(replace(spec, metadata={}))


def test_an_hfb_is_not_condemned_by_its_requires(surface_water):
    """`requires=("disv",)` means ORDERING, not an index reference -- hence a new key.

    Reading `requires` instead would refuse every barrier whenever the discretization
    rebuilt, which is both wrong and very common (a lazy disv is always `rebuilt`).
    """

    import pandas as pd
    from shapely.geometry import LineString

    from myflopy.grid_identity import INDEXED_PACKAGES, package_moves

    small, _ = _flat_meshes()
    surfaces = pd.DataFrame({0: np.full(small.ncpl, 60.0), 1: np.full(small.ncpl, 0.0)})
    context = ModelContext(
        grid=small, domain=np.ones((1, small.ncpl), dtype=int), surfaces=surfaces
    )
    xmin, ymin, _xmax, ymax = small.gdf_vorPolys.total_bounds
    hfb = mf.hfb.line(
        LineString([(xmin + 250, ymin), (xmin + 250, ymax)]), context=context, hydchr=1e-6
    )
    assert INDEXED_PACKAGES not in (hfb.metadata or {})

    # SYNTHETIC, and it has to be: the distinction only shows on a package that reads
    # `unchanged` AND carries `requires=`. A real `mf.hfb` is PINNED, so it reads `stale`
    # and the second pass skips it either way -- which is why the first version of this
    # test could not tell the two keys apart and a mutation reading `requires=` passed
    # all sixteen. This spec is what an ordering dependency looks like without the pin
    # (`mf.hfb` carried `requires=("disv",)` until ledger 378).
    import flopy

    from myflopy.specs import PackageSpec

    ordering_only = PackageSpec(
        "ordered", flopy.mf6.ModflowGwfnpf, {"k": 10.0}, requires=("disv",)
    )
    stack = mf.LayerStack(top=mf.Surface.flat(60.0))
    stack.add("a", thickness=60.0)
    verdicts = package_moves([mf.disv(stack), hfb, ordering_only])
    assert verdicts["disv"][0] == "rebuilt", "the fixture needs a provider that MOVED"
    assert verdicts["ordered"][0] == "unchanged", (
        "a `requires=` dependency is an ORDERING requirement, not an index reference -- "
        f"reading it as one refuses this package: {verdicts['ordered'][1]}"
    )
