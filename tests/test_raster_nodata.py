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
