"""`mf.Raster(path, nodata=...)` -- declaring a sentinel the file forgot to.

A raster clipped to a boundary but exported WITHOUT a no-data value in its
header is a very common GIS product: the fill (usually 0) is then read as a real
elevation. Measured on a real project 2026-08-27 -- 4.3 million pixels of exact
0.0 outside the domain, read as ground at sea level, which dragged every contact
capped to it down and collapsed all three layers to `min_sep`.

`nodata=` declares it at read time so the file need not be rewritten.
"""

from __future__ import annotations

from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pytest
import shapely
from shapely.geometry import box

import myflopy as mf


@pytest.fixture
def dem(tmp_path):
    """8x8 raster, left half real ground (500), right half filled 0, NO nodata."""

    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    values = np.full((8, 8), 500.0, dtype="float32")
    values[:, 4:] = 0.0
    path = tmp_path / "clipped_but_unlabelled.tif"
    with rasterio.open(
        path, "w", driver="GTiff", height=8, width=8, count=1, dtype="float32",
        crs="EPSG:2927", transform=from_origin(0, 8, 1, 1), nodata=None,
    ) as handle:
        handle.write(values, 1)
    return path


def _vor(*boxes):
    gdf = gpd.GeoDataFrame({"geometry": list(boxes)}, crs="EPSG:2927")
    c = gdf.geometry.centroid
    return SimpleNamespace(
        gdf_vorPolys=gdf, centroids=(c.x.to_numpy(), c.y.to_numpy()),
        crs="EPSG:2927", ncpl=len(gdf),
    )


@pytest.mark.parametrize("method", ["area", "centroid"])
def test_the_fill_is_masked_when_declared(dem, method):
    """Both samplers, because they mask in different places."""

    vor = _vor(box(0, 0, 4, 8), box(4, 0, 8, 8))     # real half, filled half

    unlabelled = np.asarray(mf.Raster(dem).values(vor, method=method), float)
    declared = np.asarray(mf.Raster(dem, nodata=0).values(vor, method=method), float)

    assert unlabelled.tolist() == [500.0, 0.0], "0 is read as ground without nodata="
    assert declared[0] == 500.0
    assert np.isnan(declared[1]), "the filled half must come back NaN"


def test_the_centroid_fallback_honours_the_override(dem):
    """The bug this test exists for, found while writing the feature.

    A cell whose every pixel is masked counts as UNCOVERED, so it falls through
    to the centroid fallback at the end of `_area_weighted_sample` -- which
    initially did not forward `nodata` and handed the sentinel straight back,
    undoing the masking for precisely the cells that needed it.
    """

    fully_filled = _vor(box(5, 1, 7, 3))             # entirely in the 0 half
    got = np.asarray(mf.Raster(dem, nodata=0).values(fully_filled, method="area"), float)
    assert np.isnan(got[0])


def test_a_straddling_cell_averages_only_the_valid_pixels(dem):
    """Masking happens BEFORE the per-cell mean, so an edge cell is the mean of
    its real pixels -- not a blend that would read as plausible-but-wrong ground."""

    straddles = _vor(box(2, 0, 6, 8))                # half real, half filled
    blended = np.asarray(mf.Raster(dem).values(straddles, method="area"), float)
    masked = np.asarray(mf.Raster(dem, nodata=0).values(straddles, method="area"), float)

    assert blended[0] == pytest.approx(250.0), "unmasked: averaged with the fill"
    assert masked[0] == pytest.approx(500.0), "masked: only the real pixels count"


def test_a_declared_header_nodata_still_wins_when_no_override(tmp_path):
    """`nodata=` overrides; it must not be required to get the normal behaviour."""

    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path = tmp_path / "proper.tif"
    values = np.full((4, 4), 500.0, dtype="float32")
    values[:, 2:] = -9999.0
    with rasterio.open(
        path, "w", driver="GTiff", height=4, width=4, count=1, dtype="float32",
        crs="EPSG:2927", transform=from_origin(0, 4, 1, 1), nodata=-9999.0,
    ) as handle:
        handle.write(values, 1)

    got = np.asarray(mf.Raster(path).values(_vor(box(2, 0, 4, 4))), float)
    assert np.isnan(got[0])


def test_nodata_defaults_to_none():
    assert mf.Raster("x.tif").nodata is None


# --- NaN cells must not poison an interpolated surface ---------------------- #

INTERPOLATORS = ["griddata", "rbf", "linearND", "cloughTocher2D"]


def test_the_parametrized_tests_name_every_interpolator():
    """The first fix reached two interpolators of four, because a fix at
    `_finite_samples` reaches only what calls it. A fifth must join the list."""

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    assert sorted(InterpolatedSurface()._interpolators) == sorted(INTERPOLATORS)


def test_a_surface_interpolates_around_cells_with_no_value():
    """Reported as "surfaces are not plotting". They were built and were entirely
    NaN, so plotly drew nothing.

    `griddata` triangulates every point it is handed, so ONE NaN vertex makes
    every output cell of every triangle touching it NaN. Scattered NaNs -- a DEM
    masked to a boundary, inactive cells in a head array -- therefore wipe out
    the whole surface. Measured on a real stack: 300 valid cells of 582 produced
    0 finite pixels of 14,400.

    Pre-existing rather than new: any raster with a properly declared no-data
    value hit it too. `mf.Raster(nodata=)` merely made it easy to reach.
    """

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    rng = np.random.default_rng(0)
    xs = rng.uniform(0, 100, 300)
    ys = rng.uniform(0, 100, 300)
    zs = 50 + 0.1 * xs
    zs[rng.uniform(size=zs.size) < 0.4] = np.nan       # 40% of cells have no value

    grid = InterpolatedSurface(
        xs=xs, ys=ys, zs=zs, surf_type="lyr", resolution=40, crs="EPSG:2927"
    ).surface
    finite = np.isfinite(np.asarray(grid, float))

    assert finite.any(), "NaN samples wiped out the entire interpolated surface"
    assert finite.mean() > 0.5


@pytest.mark.parametrize("interpolator", INTERPOLATORS)
def test_every_interpolator_interpolates_around_cells_with_no_value(interpolator):
    """The test above names no interpolator, so it only ever reached `griddata`.
    `linearND` and `cloughTocher2D` still triangulated the raw arrays, and on
    these samples left 292 and 0 finite pixels of 1,600: a head surface with
    inactive cells, or a masked DEM, drawn mostly or entirely empty.

    Inside the hull of the cells that DO have a value no pixel has a reason to be
    NaN, so that part is exact.
    """

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    rng = np.random.default_rng(0)
    xs = rng.uniform(0, 100, 300)
    ys = rng.uniform(0, 100, 300)
    zs = 50 + 0.1 * xs
    zs[rng.uniform(size=zs.size) < 0.4] = np.nan       # 40% of cells have no value

    surf = InterpolatedSurface(
        xs=xs, ys=ys, zs=zs, surf_type="lyr", resolution=40, crs="EPSG:2927",
        interpolator=interpolator,
    )
    assert surf.interpolator == interpolator, "an unknown name falls back to griddata"
    finite = np.isfinite(np.asarray(surf.surface, float))

    assert finite.mean() > 0.5, f"{finite.sum()} finite pixels of {finite.size}"
    valued = np.isfinite(zs)
    hull = shapely.MultiPoint(np.column_stack((xs[valued], ys[valued]))).convex_hull
    inside = shapely.contains_xy(hull, *surf.xy_meshgrid)
    assert finite[inside].all(), "a cell with no value left a hole among cells that have one"


def test_an_all_nan_surface_says_so_rather_than_drawing_nothing():
    """The degenerate case still has to be an error, not an empty picture."""

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    xs = np.linspace(0, 10, 25)
    surf = InterpolatedSurface(
        xs=xs, ys=xs, zs=np.full(25, np.nan), surf_type="lyr",
        resolution=10, crs="EPSG:2927",
    )
    with pytest.raises(ValueError, match="every cell in this surface is NaN"):
        _ = surf.surface


@pytest.mark.parametrize("interpolator", INTERPOLATORS)
def test_every_interpolator_says_so_when_every_cell_is_nan(interpolator):
    """`linearND` and `cloughTocher2D` returned an all-NaN surface here without a
    word -- the empty picture the error exists to replace. The points are
    scattered, not in a line: collinear points fail to triangulate and raise with
    or without the fix, so they never reach the silent case."""

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    rng = np.random.default_rng(2)
    surf = InterpolatedSurface(
        xs=rng.uniform(0, 10, 25), ys=rng.uniform(0, 10, 25), zs=np.full(25, np.nan),
        surf_type="lyr", resolution=10, crs="EPSG:2927", interpolator=interpolator,
    )
    assert surf.interpolator == interpolator, "an unknown name falls back to griddata"
    with pytest.raises(ValueError, match="every cell in this surface is NaN"):
        _ = surf.surface


def test_a_surface_with_no_nans_is_unchanged():
    """The fix must not perturb the ordinary path."""

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    rng = np.random.default_rng(1)
    xs, ys = rng.uniform(0, 100, 200), rng.uniform(0, 100, 200)
    zs = 50 + 0.1 * xs + 0.05 * ys
    grid = np.asarray(InterpolatedSurface(
        xs=xs, ys=ys, zs=zs, surf_type="lyr", resolution=30, crs="EPSG:2927"
    ).surface, float)

    inside = grid[np.isfinite(grid)]
    assert inside.size > 0
    assert 50 <= inside.min() and inside.max() <= 50 + 0.1 * 100 + 0.05 * 100 + 1e-6


# --- MF6's no-value sentinels must not reach a head surface ------------------ #
# `sentinel_run` (tests/conftest.py) is one real MF6 run with a DRY patch, an
# idomain-0 block and heads on both sides of 10,000 ft.


def test_a_head_surface_blanks_dry_and_inactive_cells_and_keeps_every_real_head(
    sentinel_run,
):
    """``zs`` cut MF6's sentinels at ``> 10_000``, which is wrong in both directions.

    A DRY cell is -1e30, which is not above 10,000, so it stayed a sample; and a
    real head above 10,000 ft is, so it was erased. The rule is the one every
    other reader of the head file uses, ``|value| >= 1e29``, so the surface and
    ``hds.array()``/``hds.get()`` must agree cell for cell (ledger 315).
    """

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    raw = np.asarray(sentinel_run.hds.array(masked=False), dtype=float)
    real = raw[np.abs(raw) < 1e29]
    assert (raw == -1e30).any(), "the fixture must hold a DRY cell"
    assert (raw == 1e30).any(), "the fixture must hold an INACTIVE cell"
    assert (real > 10_000).any() and (real < 10_000).any(), "heads on both sides of 10,000"

    zs = np.asarray(
        InterpolatedSurface(model=sentinel_run, kstpkper=(0, 0), layer=0, resolution=40).zs,
        dtype=float,
    ).ravel()

    assert np.isnan(zs[raw == -1e30]).all(), "a dry cell's -1e30 was kept as a head"
    assert np.isnan(zs[raw == 1e30]).all(), "an inactive cell's +1e30 was kept as a head"
    np.testing.assert_array_equal(zs[np.abs(raw) < 1e29], real)   # every real head, kept
    np.testing.assert_array_equal(zs, sentinel_run.hds.array())
    np.testing.assert_array_equal(zs, sentinel_run.hds.get(per=0, layer=0)["head"].to_numpy())


@pytest.mark.parametrize("interpolator", INTERPOLATORS)
def test_every_interpolator_keeps_a_head_surface_inside_the_real_heads(
    sentinel_run, interpolator
):
    """What is DRAWN, under each interpolator: no crater, and no ceiling at 10,000.

    Within a foot of the real heads, not exactly inside them: ``rbf`` overshoots a
    little on any field. Measured with the old cutoff on a 207-cell run: one dry
    patch put every finite ``griddata`` pixel outside the real 30-50 ft range,
    cubic overshoot reaching +1.6e29 as well as -1.2e30; and heads of
    9,950-10,050 ft drew a surface that stopped at 9,998
    (``docs/reproducers/head_surface_sentinels.py``).
    """

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    raw = np.asarray(sentinel_run.hds.array(masked=False), dtype=float)
    real = raw[np.abs(raw) < 1e29]
    surf = InterpolatedSurface(
        model=sentinel_run, kstpkper=(0, 0), layer=0, resolution=40,
        interpolator=interpolator,
    )
    assert surf.interpolator == interpolator, "an unknown name falls back to griddata"
    drawn = np.asarray(surf.fig.data[0].z, dtype=float)
    drawn = drawn[np.isfinite(drawn)]

    assert drawn.size > 0.5 * 40 * 40, f"{drawn.size} finite pixels of {40 * 40}"
    assert drawn.min() >= real.min() - 1.0, f"a crater to {drawn.min():.4g}"
    assert drawn.max() <= real.max() + 1.0, f"a spike to {drawn.max():.4g}"
    assert drawn.max() >= real.max() - 1.0, (
        f"the surface stops at {drawn.max():,.2f}; the highest head is {real.max():,.2f}"
    )


def test_a_sentinel_is_either_sign_of_1e30_and_never_a_real_value():
    """The rule itself. One sign only is the trap: MF6 writes +1e30 for an
    inactive cell and -1e30 for a dry one, and each reader that tested one
    side kept the other."""

    from myflopy.modflow.mf6.headsplus import _is_sentinel

    values = [1e30, -1e30, 1e29, -1e29, 12_000.0, -12_000.0, 9.9e28, 0.0, np.nan]
    assert _is_sentinel(values).tolist() == [
        True, True, True, True, False, False, False, False, False
    ]
