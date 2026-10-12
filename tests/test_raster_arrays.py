"""A raster fills a per-cell array parameter, lazily and saveably (ledger 541).

``mf.RasterValues`` is the raster sibling of ``mf.Polygons`` (ledger 299): a RECIPE
resolved at BUILD against whatever grid the model supplies, through the one raster
sampler. Before it, a K, storage or recharge raster reached a model only as an
``ndarray`` sampled by hand -- ``Surface.raster(path).values(vor)`` -- which names
cells by POSITION, so the model could not move to another mesh (``unverifiable``) and
``Project.save()`` persisted no simulation at all (``to_dict`` refuses an ndarray,
ledger 193). The arrays here are asserted EQUAL to that hand-sampled array, so moving
onto recipes changes no number.

A product (``recipe * factor``) multiplies factors cell by cell, which is how a polygon
multiplier is declared without a hand-edited raster: ``RasterValues(...) *
Polygons(path, value=0.1, default=1.0)``.
"""

from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import myflopy as mf
from myflopy import geopackage
from myflopy.geopackage import array_from_raster, resolve_array_recipe
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.modflow.mf6.grid import surfaces as grid_surfaces
from myflopy.sources import (
    RASTER_MEANS,
    ArrayProduct,
    Polygons,
    RasterValues,
    source_from_dict,
)
from myflopy.specs import (
    ModelContext,
    PackageSpec,
    SimulationSpec,
    _json_value,
    _spec_value,
    package_digest,
)


def _grid(n):
    return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1))


@pytest.fixture(scope="module")
def meshes():
    """Two DIFFERENT meshes -- the pair is what a recipe has to survive."""

    return _grid(6), _grid(10)


def _tif(path, vor, fn=None, *, pad=50.0, px=7.0, nodata=-9999.0, bands=1):
    """A GeoTIFF over ``vor``'s footprint (plus ``pad``), pixel value ``fn(x, y)``."""

    fn = fn or (lambda x, y: 1.0 + (x - x.min()) / 100.0)
    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    xmin, ymin, xmax, ymax = xmin - pad, ymin - pad, xmax + pad, ymax + pad
    width, height = int(np.ceil((xmax - xmin) / px)), int(np.ceil((ymax - ymin) / px))
    x, y = np.meshgrid(xmin + (np.arange(width) + 0.5) * px, ymax - (np.arange(height) + 0.5) * px)
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=bands,
        dtype="float32", crs=vor.crs, transform=from_origin(xmin, ymax, px, px),
        nodata=nodata,
    ) as dst:
        for band in range(1, bands + 1):
            dst.write((fn(x, y) * band).astype("float32"), band)
    return path


@pytest.fixture(scope="module")
def k_tif(tmp_path_factory, meshes):
    """Written over the LARGER mesh, so it covers both."""

    _, large = meshes
    return _tif(tmp_path_factory.mktemp("rv") / "k.tif", large)


def _stack(nlay=1):
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    for index in range(nlay):
        stack.add(f"u{index}", thickness=50.0)
    return stack


def _sim(grid, package, *, nlay=1, extra=()):
    return mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid),
               packages=[mf.disv(_stack(nlay)), package, *extra]),
        name="s", packages=[mf.tdis(nper=1)],
    )


def _built_array(grid, package, tmp_path, name, *, attribute, pname="npf", nlay=1, extra=()):
    extra = extra or ((mf.ic(strt=5.0),) if pname != "ic" else ())
    built = _sim(grid, package, nlay=nlay, extra=extra).build(tmp_path / name)
    return getattr(built.simulation.get_model("m").get_package(pname), attribute).get_data()


# --------------------------------------------------------------------------- the point


def test_one_declaration_fills_k_on_two_meshes_exactly_as_the_hand_sample(meshes, k_tif, tmp_path):
    """The reason this exists, and the promise that it changes no number.

    The arithmetic default reproduces ``Surface.raster(path).values(vor)`` -- the array
    the user samples by hand today -- BIT FOR BIT, on each of two meshes, from ONE
    declaration.
    """

    recipe = mf.RasterValues(k_tif)
    for index, grid in enumerate(meshes):
        k = _built_array(grid, mf.npf(k=recipe), tmp_path, f"m{index}", attribute="k")
        assert k.shape == (1, grid.ncpl)
        assert np.array_equal(k[0], mf.Surface.raster(k_tif).values(grid))


def test_a_per_layer_list_puts_each_raster_in_its_own_layer(meshes, tmp_path):
    small, _ = meshes
    paths = [_tif(tmp_path / f"k_{n}.tif", small, lambda x, y, n=n: n + 0.0 * x) for n in (1, 2, 3)]
    k = _built_array(
        small,
        mf.npf(k=[mf.RasterValues(paths[0]), 7.5, mf.RasterValues(paths[2])]),
        tmp_path, "layers", attribute="k", nlay=3,
    )
    assert k.shape == (3, small.ncpl)
    assert np.allclose(k[0], 1.0) and np.allclose(k[1], 7.5) and np.allclose(k[2], 3.0)


def test_a_polygon_multiplier_in_one_layer_is_the_hand_computed_product(meshes, k_tif, tmp_path):
    """``raster * Polygons(value=0.1, default=1.0)``: 0.1 where the pit has the MAJORITY.

    The expectation is computed independently, from each cell's own intersection with
    the pit: with ``default=`` the part of a cell the pit does not cover competes as the
    default (ledger 551), so under ``majority`` the raster is multiplied by
    0.1 exactly in the cells MORE than half inside, under ``mean`` by the area blend
    ``0.1 * share + 1.0 * (1 - share)``, and under ``centroid`` in the cells whose centre
    is inside. Before, ``majority`` multiplied every cell the pit OVERLAPS at all (11
    here where 6 are more than half inside, ledger 541) -- the ``touched`` count, which
    this asserts is larger, so a revert to the old rule fails here by name.
    """

    small, _ = meshes
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    pit_box = box(xmin, ymin, xmin + 0.4 * (xmax - xmin), ymin + 0.5 * (ymax - ymin))
    pit = tmp_path / "pit.gpkg"
    gpd.GeoDataFrame({"name": ["pit"]}, geometry=[pit_box], crs=small.crs).to_file(pit)

    sampled = mf.Surface.raster(k_tif).values(small)
    cells = small.gdf_vorPolys.geometry
    share = np.array([cell.intersection(pit_box).area / cell.area for cell in cells])
    touched = share > 1e-9
    majority = share > 0.5
    centred = np.array([pit_box.contains(cell.centroid) for cell in cells])
    assert 0 < majority.sum() < touched.sum() < small.ncpl
    assert (0 < share[touched & ~majority]).all() and (share[touched & ~majority] < 0.5).all()

    k33 = _built_array(
        small,
        mf.npf(k=1.0, k33=[2.0, mf.RasterValues(k_tif) * mf.Polygons(pit, value=0.1, default=1.0)]),
        tmp_path, "pit", attribute="k33", nlay=2,
    )
    assert np.allclose(k33[0], 2.0)
    assert np.array_equal(k33[1], sampled * np.where(majority, 0.1, 1.0))

    by_centre = _built_array(
        small,
        mf.npf(k=1.0, k33=mf.RasterValues(k_tif)
               * mf.Polygons(pit, value=0.1, default=1.0, rule="centroid")),
        tmp_path, "pit_centre", attribute="k33",
    )
    assert np.array_equal(by_centre[0], sampled * np.where(centred, 0.1, 1.0))

    blended = _built_array(
        small,
        mf.npf(k=1.0, k33=mf.RasterValues(k_tif)
               * mf.Polygons(pit, value=0.1, default=1.0, rule="mean")),
        tmp_path, "pit_mean", attribute="k33",
    )
    covered = np.where(share > 1 - 1e-9, 1.0, np.where(share < 1e-9, 0.0, share))
    assert np.allclose(blended[0], sampled * (0.1 * covered + 1.0 * (1 - covered)), rtol=1e-12)


def test_a_product_with_a_number_scales_every_cell(meshes, k_tif, tmp_path):
    small, _ = meshes
    strt = _built_array(
        small, mf.ic(strt=0.5 * mf.RasterValues(k_tif)), tmp_path, "ic",
        attribute="strt", pname="ic", extra=(mf.npf(k=1.0),),
    )
    assert np.array_equal(strt[0], 0.5 * mf.Surface.raster(k_tif).values(small))


def test_sto_takes_a_raster_and_icelltype_still_takes_a_plain_value(meshes, k_tif, tmp_path):
    small, _ = meshes
    ss = _built_array(
        small, mf.sto(ss=mf.RasterValues(k_tif) * 1e-6, sy=0.2, transient={0: True}),
        tmp_path, "sto", attribute="ss", pname="sto", extra=(mf.npf(k=1.0), mf.ic(strt=5.0)),
    )
    assert np.allclose(ss[0], mf.Surface.raster(k_tif).values(small) * 1e-6)


# ------------------------------------------------------------------ moves and saves


def test_a_raster_backed_array_moves_rebuilt_where_a_sampled_one_cannot(meshes, k_tif, tmp_path):
    """``unverifiable`` -> ``rebuilt`` (ledgers 295, 299)."""

    small, large = meshes
    plain = _sim(small, mf.npf(k=mf.Surface.raster(k_tif).values(small)))
    assert _verdict(plain.with_grid(large), "npf") == "unverifiable"

    moved = _sim(small, mf.npf(k=mf.RasterValues(k_tif) * mf.RasterValues(k_tif))).with_grid(large)
    assert _verdict(moved, "npf") == "rebuilt"
    k = moved.build(tmp_path / "moved").simulation.get_model("m").get_package("npf").k.get_data()
    assert np.allclose(k[0], mf.Surface.raster(k_tif).values(large) ** 2)


def _verdict(sim, package):
    return next(row["how"] for row in sim.move_report() if row["package"] == package)


def _project(root, grid, package, name="p"):
    project = mf.Project(root, name=name)
    project.add_grid("mesh", grid)
    sim = mf.simulation(
        mf.gwf("m", grid=mf.grid_ref("mesh"),
               packages=[mf.disv(_stack()), package, mf.ic(strt=5.0)]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    project.add_simulation(sim)
    return project, sim


def test_a_project_carrying_a_raster_recipe_saves_and_loads(meshes, k_tif, tmp_path):
    """The hand-sampled array saves NOTHING; the recipe saves and rebuilds the same K."""

    small, _ = meshes
    eager = mf.npf(k=mf.Surface.raster(k_tif).values(small))
    lazy = mf.npf(k=[mf.RasterValues(k_tif) * mf.Polygons(_whole(tmp_path, small), value=2.0)])
    for index, (package, expected) in enumerate(((eager, 0), (lazy, 1))):
        project, _ = _project(tmp_path / f"proj{index}", small, package)
        assert project.save().simulations == expected
        loaded = mf.Project.load(tmp_path / f"proj{index}")
        assert len(loaded.simulations) == expected
    k = _run_array(loaded, "npf", "k")
    assert np.array_equal(k[0], 2.0 * mf.Surface.raster(k_tif).values(small))


def _run_array(project, pname, attribute, run="r"):
    """Build a run of the project's one simulation and read one written array."""

    built = project.prepare_run(run, next(iter(project.simulations.values()))).build()
    return getattr(built.simulation.get_model("m").get_package(pname), attribute).get_data()


def _whole(tmp_path, grid):
    path = tmp_path / "whole.gpkg"
    if not path.exists():
        gpd.GeoDataFrame(
            {"name": ["all"]}, geometry=[box(*grid.gdf_vorPolys.total_bounds).buffer(1.0)],
            crs=grid.crs,
        ).to_file(path)
    return path


def test_a_saved_product_round_trips_through_json_and_yaml(meshes, k_tif, tmp_path):
    small, _ = meshes
    product = mf.RasterValues(k_tif, band=1, mean="geometric", default=2.0) * mf.Polygons(
        _whole(tmp_path, small), value=0.5, default=1.0
    ) * 3.0
    revived = _spec_value(_json_value(product))
    assert revived == product
    sim = _sim(small, mf.npf(k=[product]))
    assert SimulationSpec.from_yaml(sim.to_yaml()).to_dict() == sim.to_dict()


def test_a_frozen_run_thaws_with_its_raster_deleted(meshes, k_tif, tmp_path):
    """The run freezes the SAMPLED arrays, not the recipe (ledgers 251, 270)."""

    from myflopy.project.freeze import freeze_simulation, thaw_simulation

    small, _ = meshes
    copied = tmp_path / "throwaway.tif"
    copied.write_bytes(k_tif.read_bytes())
    sim = _sim(small, mf.npf(k=mf.RasterValues(copied)), extra=(mf.ic(strt=5.0),))
    built = sim.build(tmp_path / "freeze").built
    assert isinstance(built.models["m"].resolved["npf"].options["k"], np.ndarray)
    freeze_simulation(
        sim, tmp_path / "snap", run_name="s",
        contexts={name: model.context for name, model in built.models.items()},
        resolved={name: dict(model.resolved or {}) for name, model in built.models.items()},
    )
    copied.unlink()
    packages = {p.name: p for p in thaw_simulation(tmp_path / "snap").models[0].packages}
    assert np.array_equal(
        np.asarray(packages["npf"].options["k"])[0], mf.Surface.raster(k_tif).values(small)
    )


def test_bundling_archives_a_raster_inside_a_product(meshes, tmp_path):
    small, _ = meshes
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    raster = _tif(outside / "k33.tif", small)
    pit = _whole(outside, small)
    project, _ = _project(
        tmp_path / "proj", small,
        mf.npf(k=1.0, k33=mf.RasterValues(raster) * mf.Polygons(pit, value=0.1, default=1.0)),
    )
    report = project.save(bundle_sources=True)
    assert report.not_bundled == ()
    names = {getattr(destination, "name", str(destination).rsplit("/", 1)[-1])
             for _, destination in report.bundled}
    assert {"k33.tif", "whole.gpkg"} <= names
    assert (project.root / "inputs" / "k33.tif").exists()
    loaded = mf.Project.load(project.root)
    raster.unlink()
    pit.unlink()
    k33 = _run_array(loaded, "npf", "k33")
    assert k33.shape == (1, small.ncpl)


def test_drift_reads_the_raster_and_notices_an_edit(meshes, tmp_path):
    small, _ = meshes
    raster = _tif(tmp_path / "k.tif", small)
    project, sim = _project(tmp_path / "proj", small, mf.npf(k=mf.RasterValues(raster)))
    project.prepare_run("r1", sim).build()
    assert set(project.runs.drift()["source"]) == {"ok"}
    _tif(raster, small, lambda x, y: 9.0 + 0.0 * x)
    assert set(project.runs.drift()["source"]) == {"edited"}


def test_the_digest_reads_content_not_location(meshes, tmp_path):
    """Ledger 293's rule, for a bare recipe AND inside a product."""

    small, _ = meshes
    first = _tif(tmp_path / "a.tif", small)
    copied = tmp_path / "b.tif"
    copied.write_bytes(first.read_bytes())
    other = _tif(tmp_path / "c.tif", small, lambda x, y: 4.0 + 0.0 * x)
    for wrap in (lambda r: r, lambda r: r * 2.0):
        digests = [package_digest(mf.npf(k=wrap(mf.RasterValues(path)))) for path in (first, copied, other)]
        assert digests[0] == digests[1], "a byte-identical copy must not move the digest"
        assert digests[0] != digests[2], "a different raster must move it"


def test_a_plain_declaration_did_not_move():
    spec = mf.npf(k=[10.0, 1.0])
    assert spec.builder is flopy.mf6.ModflowGwfnpf
    assert spec.inputs == {}


# ------------------------------------------------------------ gaps, means and bands


def test_a_cell_off_the_raster_is_refused_filled_or_defaulted(meshes, tmp_path):
    """Ledger 247: a NaN reaching MF6 is a NaN budget under Normal termination."""

    small, _ = meshes
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    half = tmp_path / "half.tif"
    _tif(half, small, lambda x, y: np.where(x < (xmin + xmax) / 2, 3.0, -9999.0), pad=0.0)

    with pytest.raises(ValueError, match="gives no value for .* cells") as caught:
        array_from_raster(RasterValues(path=half), small, parameter="npf.k")
    assert "gaps='nearest'" in str(caught.value) and "default=" in str(caught.value)

    nearest = array_from_raster(RasterValues(path=half, gaps="nearest"), small, parameter="npf.k")
    assert np.allclose(nearest, 3.0)
    defaulted = array_from_raster(RasterValues(path=half, default=0.25), small, parameter="npf.k")
    assert set(np.unique(defaulted)) == {0.25, 3.0}


def test_a_raster_that_reaches_no_cell_is_refused_even_with_a_default(meshes, tmp_path):
    """``default=`` everywhere would paint a misplaced raster in silence."""

    small, _ = meshes
    far = tmp_path / "far.tif"
    shifted = rectangular_voronoi(CanonicalModelConfig(nrow=2, ncol=2, nlay=1, nper=1))
    _tif(far, shifted, pad=0.0)
    import shapely.affinity

    moved = small.gdf_vorPolys.copy()
    moved["geometry"] = moved.geometry.apply(lambda g: shapely.affinity.translate(g, 1e6, 1e6))

    class Far:
        gdf_vorPolys = moved
        crs = small.crs
        centroids = (np.asarray(small.centroids[0]) + 1e6, np.asarray(small.centroids[1]) + 1e6)

    with pytest.raises(ValueError, match="ANY of the") as caught:
        array_from_raster(RasterValues(path=far, default=1.0), Far(), parameter="npf.k")
    assert "the grid covers" in str(caught.value)


def test_geometric_and_harmonic_means_are_what_they_say(tmp_path):
    """Four pixels in ONE cell, so each mean is checkable by hand."""

    cells = gpd.GeoDataFrame(geometry=[box(0, 0, 2, 2)], crs="EPSG:2927")

    class One:
        gdf_vorPolys = cells
        crs = "EPSG:2927"
        centroids = (np.array([1.0]), np.array([1.0]))

    path = tmp_path / "four.tif"
    pixels = np.array([[1.0, 10.0], [100.0, 1000.0]], dtype="float32")
    with rasterio.open(path, "w", driver="GTiff", height=2, width=2, count=1, dtype="float32",
                       crs="EPSG:2927", transform=from_origin(0, 2, 1, 1)) as dst:
        dst.write(pixels, 1)
    values = pixels.ravel().astype(float)
    got = {mean: array_from_raster(RasterValues(path=path, mean=mean), One(), parameter="k")[0]
           for mean in RASTER_MEANS}
    assert got["arithmetic"] == pytest.approx(values.mean())
    assert got["geometric"] == pytest.approx(np.exp(np.log(values).mean()))
    assert got["harmonic"] == pytest.approx(1.0 / (1.0 / values).mean())


def test_a_geometric_mean_refuses_a_non_positive_pixel(meshes, tmp_path):
    small, _ = meshes
    path = _tif(tmp_path / "zero.tif", small, lambda x, y: np.where(x < 100, 0.0, 5.0))
    with pytest.raises(ValueError, match="at or below zero"):
        array_from_raster(RasterValues(path=path, mean="geometric"), small, parameter="npf.k")


def test_band_reads_that_band_and_a_missing_band_is_refused(meshes, tmp_path):
    small, _ = meshes
    path = _tif(tmp_path / "two.tif", small, bands=2)
    first = array_from_raster(RasterValues(path=path), small, parameter="k")
    second = array_from_raster(RasterValues(path=path, band=2), small, parameter="k")
    assert np.allclose(second, 2.0 * first)
    with pytest.raises(ValueError, match="has 2 band"):
        array_from_raster(RasterValues(path=path, band=3), small, parameter="k")


def test_the_two_mean_vocabularies_agree():
    """Re-stated a layer up because the sampler cannot import the recipe classes."""

    assert grid_surfaces.RASTER_MEANS == RASTER_MEANS


# -------------------------------------- the two accepted-then-MFDataException defects


def test_a_polygon_map_in_a_per_layer_list_builds_the_right_array(meshes, tmp_path):
    """``mf.npf(k=[mf.Polygons(...), 1.0])`` was accepted and died in FloPy at build."""

    small, _ = meshes
    k = _built_array(
        small, mf.npf(k=[mf.Polygons(_whole(tmp_path, small), value=3.0), 1.0]),
        tmp_path, "list", attribute="k", nlay=2,
    )
    assert np.allclose(k[0], 3.0) and np.allclose(k[1], 1.0)


def test_rch_array_takes_a_polygon_map_and_a_raster(meshes, k_tif, tmp_path):
    """``mf.rch.array(recharge=mf.Polygons(...))`` was accepted and died in FloPy."""

    small, _ = meshes
    for index, recipe in enumerate(
        (mf.Polygons(_whole(tmp_path, small), value=1e-3), mf.RasterValues(k_tif) * 1e-4)
    ):
        spec = mf.rch.array(recharge=recipe)
        assert "recharge" in spec.inputs
        recharge = _built_array(small, spec, tmp_path, f"rch{index}",
                                attribute="recharge", pname="rch", extra=(mf.npf(k=1.0),))
        expected = (np.full(small.ncpl, 1e-3) if index == 0
                    else mf.Surface.raster(k_tif).values(small) * 1e-4)
        assert np.allclose(np.asarray(recharge[0]).ravel(), expected)


def test_each_distinct_recipe_is_sampled_once_per_build(meshes, k_tif, tmp_path, monkeypatch):
    """Twelve monthly rasters named by many periods read twelve times, not once a period."""

    small, _ = meshes
    second = _tif(tmp_path / "r2.tif", small, lambda x, y: 2e-4 + 0.0 * x)
    calls = []
    real = geopackage.array_from_raster

    def counting(recipe, grid, *, parameter):
        calls.append(recipe.path.name)
        return real(recipe, grid, parameter=parameter)

    monkeypatch.setattr(geopackage, "array_from_raster", counting)
    first = mf.RasterValues(k_tif) * 1e-4
    periods = {period: (first if period % 2 == 0 else mf.RasterValues(second)) for period in range(6)}
    built = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small),
               packages=[mf.disv(_stack()), mf.npf(k=mf.RasterValues(k_tif)), mf.ic(strt=5.0),
                         mf.rch.array(recharge=periods)]),
        name="s", packages=[mf.tdis(nper=6)],
    ).build(tmp_path / "memo")
    assert sorted(calls) == ["k.tif", "r2.tif"]
    recharge = built.simulation.get_model("m").get_package("rch").recharge.get_data()
    assert np.allclose(np.asarray(recharge[5]).ravel(), 2e-4)


@pytest.mark.slow
def test_mf6_runs_on_raster_k_and_its_recharge_budget_is_rate_times_area(meshes, k_tif, tmp_path):
    """A real MF6 run: K from a raster, recharge from a raster, Normal termination.

    The RCH budget MF6 writes equals the sampled rate times each cell's plan area --
    every column active, so ``irch='top_active'`` is layer 0 everywhere -- over every
    cell but the two CHD cells, where MF6 applies no boundary flux (measured: 152.907
    against 160.725 with them counted).
    """

    small, _ = meshes
    rate = mf.RasterValues(k_tif) * 1e-4
    built = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small), packages=[
            mf.disv(_stack(2)),
            mf.npf(k=[mf.RasterValues(k_tif), mf.RasterValues(k_tif) * 0.1], icelltype=0),
            mf.ic(strt=90.0),
            mf.chd(stress_period_data={0: [((0, 0), 90.0), ((0, small.ncpl - 1), 90.0)]}),
            mf.rch.array(recharge=rate),
            mf.oc(budget_filerecord="m.cbc", head_filerecord="m.hds",
                  saverecord=[("BUDGET", "ALL"), ("HEAD", "ALL")]),
        ]),
        name="s", packages=[mf.tdis(nper=1)],
    ).build(tmp_path / "run")
    simulation = built.simulation
    simulation.write_simulation(silent=True)
    success, report = simulation.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])

    budget = flopy.utils.CellBudgetFile(str(tmp_path / "run" / "m.cbc"))
    record = budget.get_data(text="RCH")[0]
    budget.close()
    areas = small.gdf_vorPolys.geometry.area.to_numpy()
    free = np.ones(small.ncpl, dtype=bool)
    free[[0, small.ncpl - 1]] = False
    expected = (mf.Surface.raster(k_tif).values(small) * 1e-4 * areas)[free].sum()
    assert record["q"].sum() == pytest.approx(expected, rel=1e-6)
    head_file = flopy.utils.HeadFile(str(tmp_path / "run" / "m.hds"))
    heads = head_file.get_data()
    head_file.close()
    assert np.isfinite(heads).all()


# ------------------------------------------------------------------------- refusals


def test_the_recipe_round_trips_as_itself(k_tif):
    recipe = mf.RasterValues(k_tif, band=2, method="area", mean="harmonic", nodata=0.0,
                             gaps="nearest")
    assert source_from_dict(recipe.to_dict()) == recipe
    assert set(mf.RasterValues(k_tif).to_dict()) == {"kind", "path"}


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"gaps": "warn"}, "'warn' and 'ignore' would hand"),
        ({"gaps": "ignore"}, "is not one of"),
        ({"method": "nearest"}, "is not one of"),
        ({"mean": "median"}, "is not one of"),
        ({"method": "centroid", "mean": "geometric"}, "ONE pixel"),
        ({"band": 0}, "numbered from 1"),
        ({"band": 1.5}, "whole number"),
        ({"gaps": "nearest", "default": 1.0}, "two answers"),
        ({"default": float("nan")}, "finite"),
    ],
)
def test_a_recipe_that_cannot_resolve_is_refused_at_the_declaration(kwargs, match):
    with pytest.raises(ValueError, match=match):
        RasterValues(path="k.tif", **kwargs)


def test_polygons_takes_exactly_one_of_column_and_value():
    with pytest.raises(ValueError, match="not both"):
        Polygons(path="p.gpkg", column="k", value=1.0)
    with pytest.raises(ValueError, match="cannot be empty"):
        Polygons(path="p.gpkg")
    with pytest.raises(ValueError, match="must be a number"):
        Polygons(path="p.gpkg", value=True)


def test_an_array_or_surface_factor_is_refused(k_tif):
    recipe = mf.RasterValues(k_tif)
    with pytest.raises(TypeError, match="remap=True"):
        recipe * np.ones(4)
    with pytest.raises(TypeError, match="remap=True"):
        np.ones(4) * recipe
    with pytest.raises(TypeError, match="ELEVATION"):
        recipe * mf.Surface.flat(1.0)
    with pytest.raises(ValueError, match="at least one"):
        ArrayProduct(factors=(1.0, 2.0))


def test_a_missing_raster_is_refused_at_the_declaration(tmp_path):
    with pytest.raises(FileNotFoundError, match="Raster source not found"):
        mf.npf(k=mf.RasterValues(tmp_path / "absent.tif"))


@pytest.mark.parametrize(("helper", "parameter"), [("npf", "icelltype"), ("sto", "iconvert")])
def test_an_integer_flag_refuses_a_raster(helper, parameter, k_tif):
    kwargs = {parameter: mf.RasterValues(k_tif)} | ({"k": 1.0} if helper == "npf" else {})
    with pytest.raises(TypeError, match="integer flag"):
        getattr(mf, helper)(**kwargs)
    with pytest.raises(TypeError, match="integer flag"):
        getattr(mf, helper)(**{**kwargs, parameter: [mf.RasterValues(k_tif), 1]})


def test_irch_refuses_a_recipe(k_tif):
    with pytest.raises(TypeError, match="LAYER index"):
        mf.rch.array(recharge=1e-3, irch=mf.RasterValues(k_tif))


def test_a_per_layer_list_refuses_an_array_entry_and_a_layered_map(k_tif, tmp_path, meshes):
    small, _ = meshes
    with pytest.raises(TypeError, match="entry 1 of a per-layer list"):
        mf.npf(k=[mf.RasterValues(k_tif), np.ones(4)])
    layered = mf.Polygons(_whole(tmp_path, small), value=1.0, layer_field="layer")
    with pytest.raises(ValueError, match="paints layers of its own"):
        mf.npf(k=[mf.RasterValues(k_tif) * layered, 1.0])


def test_an_areal_slot_refuses_a_list_a_layered_map_and_a_context(k_tif, tmp_path, meshes):
    small, _ = meshes
    with pytest.raises(TypeError, match="per-layer list means nothing"):
        mf.rch.array(recharge=[mf.RasterValues(k_tif)])
    with pytest.raises(ValueError, match="Drop layer_field"):
        mf.rch.array(recharge=mf.Polygons(_whole(tmp_path, small), value=1.0, layer_field="layer"))
    with pytest.raises(TypeError, match="Drop `context=`"):
        mf.rch.array(recharge=mf.RasterValues(k_tif),
                     context=ModelContext(domain=np.ones((1, small.ncpl), dtype=int)))


def test_a_surface_in_a_per_layer_list_is_refused_naming_rastervalues():
    with pytest.raises(TypeError, match="mf.RasterValues"):
        mf.npf(k=[mf.Surface.flat(1.0), 2.0])


def test_a_per_layer_list_of_the_wrong_length_is_refused_at_build(meshes, k_tif):
    small, _ = meshes
    with pytest.raises(ValueError, match="2 entries, and the model has 3"):
        resolve_array_recipe([mf.RasterValues(k_tif), 1.0], small, nlay=3, parameter="npf.k")


def test_a_recipe_beside_a_same_named_option_is_refused(k_tif):
    spec = PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 10.0},
                       inputs={"k": mf.RasterValues(k_tif)})
    with pytest.raises(TypeError, match="supplied twice"):
        spec.build(object())


def test_numpy_scalars_reach_the_payload_as_python_values():
    """Ledger 265's class: PyYAML dispatches on the EXACT type and refuses a numpy one."""

    import yaml

    recipe = mf.RasterValues("k.tif", band=np.int64(2), nodata=np.float32(0.0),
                             default=np.float64(1.5))
    product = recipe * np.float64(2.0) * mf.Polygons("p.gpkg", value=np.float32(0.5))
    payload = _json_value(product)
    yaml.safe_dump(payload)
    factors = payload[ArrayProduct.TOKEN]["factors"]
    assert type(factors[0]["band"]) is int and type(factors[0]["default"]) is float
    assert type(factors[1]) is float and type(factors[2]["value"]) is float


def test_a_recharge_recipe_moves_rebuilt(meshes, k_tif):
    """The areal slot moves as the layered one does; a per-column ``irch`` array does not."""

    small, large = meshes
    for irch, verdict in (("top_active", "rebuilt"), (0, "rebuilt"),
                          (np.zeros(small.ncpl, dtype=int), "unverifiable")):
        sim = _sim(small, mf.rch.array(recharge=mf.RasterValues(k_tif), irch=irch),
                   extra=(mf.npf(k=1.0), mf.ic(strt=5.0)))
        assert _verdict(sim.with_grid(large), "rch") == verdict


# --------------------------------------------------- review of R2-rasters (542)


def test_an_edit_route_refuses_a_recipe_it_would_put_in_the_options(k_tif, tmp_path):
    """Ledger 542: a recipe resolves from ``inputs`` only.

    Every edit route writes ``options``, and a recipe there reached FloPy unresolved:
    ``packages.fork("npf", "npf2", k=mf.RasterValues(p))`` constructed fine and died at
    build with an ``MFDataException`` naming neither the recipe nor the fork.
    """

    recipe = mf.RasterValues(k_tif)
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf", mf.npf(k=1.0))
    attempts = {
        "fork": lambda: project.packages.fork("npf", "npf2", k=recipe),
        "edit": lambda: project.packages.edit("npf", k33=[recipe, 1.0]),
        "with_options": lambda: mf.npf(k=1.0).with_options(k=recipe * 2.0),
        "by hand": lambda: PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": recipe}),
        "period dict": lambda: PackageSpec(
            "rch", flopy.mf6.ModflowGwfrcha, {"recharge": {0: recipe}}
        ),
    }
    for route, attempt in attempts.items():
        with pytest.raises(TypeError, match="array recipe") as caught:
            attempt()
        assert "through its helper" in str(caught.value), route
    assert "npf2" not in project.packages
    assert project.packages["npf"].options["k"] == 1.0


def test_nearest_on_a_raster_blank_over_the_grid_names_the_raster(meshes, tmp_path):
    """``gaps='nearest'`` with no valued cell died in a k-d tree with ``IndexError``."""

    small, _ = meshes
    blank = _tif(tmp_path / "blank.tif", small, lambda x, y: -9999.0 + 0.0 * x)
    with pytest.raises(ValueError, match="no valued cell") as caught:
        grid_surfaces.get_raster_vals_at_centroids(small, [blank], ["v"], gaps="nearest")
    assert "blank.tif" in str(caught.value)


@pytest.mark.slow
def test_mf6_reads_a_polygon_multiplier_only_where_the_pit_has_the_majority(meshes, k_tif, tmp_path):
    """A real MF6 run of ``raster * Polygons(value=0.1, default=1.0)`` (ledger 551).

    MF6 echoes the K it READ through NPF's ``EXPORT_ARRAY_ASCII``; that array is the
    raster times 0.1 exactly in the cells the pit covers more than half of and the
    raster elsewhere -- not in every cell the pit touches -- and the run terminates
    normally with finite heads.
    """

    small, _ = meshes
    xmin, ymin, xmax, ymax = small.gdf_vorPolys.total_bounds
    pit_box = box(xmin, ymin, xmin + 0.4 * (xmax - xmin), ymin + 0.5 * (ymax - ymin))
    pit = tmp_path / "pit.gpkg"
    gpd.GeoDataFrame({"name": ["pit"]}, geometry=[pit_box], crs=small.crs).to_file(pit)
    share = np.array([cell.intersection(pit_box).area / cell.area
                      for cell in small.gdf_vorPolys.geometry])
    assert (share > 0).sum() > (share > 0.5).sum() > 0

    built = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=small), packages=[
            mf.disv(_stack()),
            mf.npf(k=mf.RasterValues(k_tif) * mf.Polygons(pit, value=0.1, default=1.0),
                   icelltype=0, export_array_ascii=True),
            mf.ic(strt=90.0),
            mf.chd(stress_period_data={0: [((0, 0), 90.0), ((0, small.ncpl - 1), 80.0)]}),
            mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
        ]),
        name="s", packages=[mf.tdis(nper=1)],
    ).build(tmp_path / "run")
    simulation = built.simulation
    simulation.write_simulation(silent=True)
    success, report = simulation.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])

    exported = tmp_path / "run" / "m-npf.k.l1.txt"   # MF6's own echo, not FloPy's input
    read_by_mf6 = np.loadtxt(exported).ravel()
    expected = mf.Surface.raster(k_tif).values(small) * np.where(share > 0.5, 0.1, 1.0)
    assert read_by_mf6 == pytest.approx(expected, rel=1e-6)
    head_file = flopy.utils.HeadFile(str(tmp_path / "run" / "m.hds"))
    heads = head_file.get_data()
    head_file.close()
    assert np.isfinite(heads).all()
