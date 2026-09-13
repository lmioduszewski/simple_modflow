"""A cell the raster gave no value for is reported, not handed back as a bare NaN.

`rasterio` yields the no-data value for a point outside coverage, which the
sampler turns into NaN, and the area sampler already falls back to the centroid
for a cell holding no pixel CENTRE. What was left was a cell whose centroid is off
the raster's data too: it came back NaN and the frame was returned with nothing
said. Sampled into `strt` or `gdf_topbtm`, that NaN reaches MF6, which reports
success and returns a NaN budget.

The motivating case is remapping between grids: `InterpolatedSurface` writes a
raster with `griddata`, which fills only the CONVEX HULL of its input points, so a
raster built from one grid's CELL CENTROIDS is blank in a rim half a cell wide
around the whole model -- and a boundary cell of a different grid lands in it.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from myflopy.modflow.mf6.grid.surfaces import get_raster_vals_at_centroids

CRS = "EPSG:2927"


def _raster(path, data, *, nodata=-9999.0, west=0.0, north=4.0, px=1.0, crs=CRS):
    """A small north-up GeoTIFF; pixel centres sit at west+0.5, north-0.5, ..."""

    data = np.asarray(data, dtype="float64")
    h, w = data.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=h, width=w, count=1, dtype="float64",
        crs=crs, transform=from_origin(west, north, px, px), nodata=nodata,
    ) as dst:
        dst.write(data, 1)
    return path


def _vor(polys, crs=CRS):
    """A grid stand-in backed by real cell polygons."""

    gdf = gpd.GeoDataFrame({"geometry": list(polys)}, crs=crs)
    cents = gdf.geometry.centroid
    return SimpleNamespace(
        gdf_vorPolys=gdf,
        centroids=(cents.x.to_numpy(), cents.y.to_numpy()),
        crs=crs,
        ncpl=len(gdf),
    )


def _row_of_cells(n=4):
    """``n`` 1 ft cells along y in [3, 4], so their centroids sit on the top pixel row."""

    return _vor([box(i, 3, i + 1, 4) for i in range(n)])


@pytest.fixture
def gappy(tmp_path):
    """A 4x4 raster whose top row is valued except its LAST pixel, which is no-data."""

    data = np.array(
        [
            [10.0, 20.0, 30.0, -9999.0],
            [40.0, 50.0, 60.0, 70.0],
            [80.0, 90.0, 100.0, 110.0],
            [120.0, 130.0, 140.0, 150.0],
        ]
    )
    return _raster(tmp_path / "gappy.tif", data)


def test_a_gap_is_warned_about_by_default(gappy):
    """The warning the user never got, naming the count and the cell."""

    with pytest.warns(UserWarning, match=r"gave no value for 1 of 4 cells"):
        got = get_raster_vals_at_centroids(_row_of_cells(), [gappy], ["z"], method="centroid")
    assert np.isnan(got["z"].to_numpy()[3])


def test_the_warning_names_the_cell_and_the_consequence(gappy):
    """It says which cell, and why a NaN matters, because MF6 will not say."""

    with pytest.warns(UserWarning, match=r"column 'z', which are NaN: 3"):
        get_raster_vals_at_centroids(_row_of_cells(), [gappy], ["z"], method="centroid")
    with pytest.warns(UserWarning, match=r"NaN budget"):
        get_raster_vals_at_centroids(_row_of_cells(), [gappy], ["z"], method="centroid")


def test_ignore_is_silent_and_keeps_the_nan(gappy, recwarn):
    """For a caller whose next stage fills gaps -- Surface(fill='propagate')."""

    got = get_raster_vals_at_centroids(
        _row_of_cells(), [gappy], ["z"], method="centroid", gaps="ignore"
    )
    assert np.isnan(got["z"].to_numpy()[3])
    assert [w for w in recwarn if "gave no value" in str(w.message)] == []


def test_raise_refuses_the_sample(gappy):
    """For a caller that would rather stop than carry a hole."""

    with pytest.raises(ValueError, match=r"gave no value for 1 of 4 cells"):
        get_raster_vals_at_centroids(
            _row_of_cells(), [gappy], ["z"], method="centroid", gaps="raise"
        )


def test_nearest_takes_the_gap_from_its_nearest_valued_cell(gappy):
    """Cell 3's nearest valued neighbour is cell 2, which sampled 30.0."""

    got = get_raster_vals_at_centroids(
        _row_of_cells(), [gappy], ["z"], method="centroid", gaps="nearest"
    )["z"].to_numpy()
    assert got.tolist() == [10.0, 20.0, 30.0, 30.0]
    assert np.isfinite(got).all()


def test_nearest_leaves_the_valued_cells_alone(gappy):
    """Filling must not perturb a cell that had an answer."""

    plain = get_raster_vals_at_centroids(
        _row_of_cells(), [gappy], ["z"], method="centroid", gaps="ignore"
    )["z"].to_numpy()
    filled = get_raster_vals_at_centroids(
        _row_of_cells(), [gappy], ["z"], method="centroid", gaps="nearest"
    )["z"].to_numpy()
    assert filled[:3].tolist() == plain[:3].tolist()


def test_nearest_logs_what_it_substituted(gappy, caplog):
    """A silent substitution is the defect being fixed, so the fill says so too."""

    with caplog.at_level(logging.INFO, logger="myflopy.modflow.mf6.grid.surfaces"):
        get_raster_vals_at_centroids(
            _row_of_cells(), [gappy], ["z"], method="centroid", gaps="nearest"
        )
    assert "filled 1 of 4 cells" in caplog.text


def test_an_unknown_mode_is_refused(gappy):
    """A misspelled mode must not read as 'ignore'."""

    with pytest.raises(ValueError, match=r"Unknown gaps mode 'fill'"):
        get_raster_vals_at_centroids(
            _row_of_cells(), [gappy], ["z"], method="centroid", gaps="fill"
        )


def test_a_raster_somewhere_else_entirely_names_both_extents(tmp_path):
    """No overlap at all is a location or CRS mismatch, reported as two boxes."""

    far = _raster(tmp_path / "far.tif", np.ones((4, 4)), west=10_000.0, north=10_004.0)
    with pytest.raises(ValueError, match=r"gave no value for ANY of the 4 cells"):
        get_raster_vals_at_centroids(_row_of_cells(), [far], ["z"], method="centroid")
    with pytest.raises(ValueError, match=r"covers .* the grid covers"):
        get_raster_vals_at_centroids(_row_of_cells(), [far], ["z"], method="centroid")


def test_an_all_nodata_raster_that_DOES_overlap_is_not_called_a_mismatch(tmp_path):
    """Blank data where the extents agree is a different bug; do not misdirect."""

    blank = _raster(tmp_path / "blank.tif", np.full((4, 4), -9999.0))
    with pytest.warns(UserWarning, match=r"gave no value for 4 of 4 cells"):
        got = get_raster_vals_at_centroids(
            _row_of_cells(), [blank], ["z"], method="centroid"
        )
    assert np.isnan(got["z"].to_numpy()).all()


def test_each_column_is_resolved_on_its_own(tmp_path):
    """Two rasters, gaps in different cells: neither answer may borrow the other's."""

    a = _raster(tmp_path / "a.tif", np.array([[1.0, 2.0, 3.0, -9999.0]] + [[9.0] * 4] * 3))
    b = _raster(tmp_path / "b.tif", np.array([[-9999.0, 5.0, 6.0, 7.0]] + [[9.0] * 4] * 3))
    got = get_raster_vals_at_centroids(
        _row_of_cells(), [a, b], ["a", "b"], method="centroid", gaps="nearest"
    )
    assert got["a"].to_numpy().tolist() == [1.0, 2.0, 3.0, 3.0]
    assert got["b"].to_numpy().tolist() == [5.0, 5.0, 6.0, 7.0]


def test_a_constant_column_is_never_a_gap(gappy, recwarn):
    """A number in raster_files is a constant, not a file; it cannot miss."""

    got = get_raster_vals_at_centroids(
        _row_of_cells(), [12.5], ["k"], method="centroid"
    )
    assert got["k"].to_numpy().tolist() == [12.5] * 4
    assert [w for w in recwarn if "gave no value" in str(w.message)] == []


def test_the_interpolated_surface_rim_is_the_motivating_case(tmp_path):
    """A raster valued only where an inner grid's centroids were.

    This is the shape of the reported failure: heads interpolated from one model's
    cell centroids leave a blank rim, and the new grid's boundary cell lands in it.
    One cell of five comes back NaN, and ``nearest`` gives it its neighbour's head
    -- which is what the reporter did by hand.
    """

    # Valued in the middle three pixels of the top row; blank at both ends.
    data = np.array([[-9999.0, 100.0, 101.0, 102.0, -9999.0]] + [[-9999.0] * 5] * 3)
    hull = _raster(tmp_path / "hds.tif", data, north=4.0)
    grid = _row_of_cells(5)

    with pytest.warns(UserWarning, match=r"gave no value for 2 of 5 cells"):
        raw = get_raster_vals_at_centroids(grid, [hull], ["strt"], method="centroid")
    assert np.isnan(raw["strt"].to_numpy()).sum() == 2

    filled = get_raster_vals_at_centroids(
        grid, [hull], ["strt"], method="centroid", gaps="nearest"
    )["strt"].to_numpy()
    assert filled.tolist() == [100.0, 100.0, 101.0, 102.0, 102.0]


def _real_grid(n=4):
    """A real VoronoiGridPlus: ``n`` 1 ft square cells in a row, y in [3, 4]."""

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    verts = np.array(
        [[float(i), 3.0] for i in range(n + 1)] + [[float(i), 4.0] for i in range(n + 1)]
    )
    iverts = [[i, n + 1 + i, n + 2 + i, i + 1] for i in range(n)]
    xcyc = np.array([[i + 0.5, 3.5] for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc, crs=CRS)


def test_the_grid_method_forwards_gaps(gappy):
    """``vor.get_raster_vals_at_centroids`` is the scope a user samples from.

    It forwarded neither ``method`` nor ``nodata``, which also put ``gaps`` out of
    reach of the only place it matters.
    """

    vor = _real_grid()
    with pytest.warns(UserWarning, match=r"gave no value for 1 of 4 cells"):
        vor.get_raster_vals_at_centroids([gappy], ["z"], method="centroid")
    filled = vor.get_raster_vals_at_centroids(
        [gappy], ["z"], method="centroid", gaps="nearest"
    )
    assert filled["z"].to_numpy().tolist() == [10.0, 20.0, 30.0, 30.0]


def test_the_grid_method_forwards_nodata(tmp_path):
    """The sibling parameter that was also unreachable: a sentinel not in the header."""

    sentinel = _raster(tmp_path / "s.tif", np.array([[1.0, 2.0, 3.0, 0.0]] + [[9.0] * 4] * 3), nodata=None)
    vor = _real_grid()
    assert vor.get_raster_vals_at_centroids([sentinel], ["z"], method="centroid")[
        "z"
    ].to_numpy().tolist() == [1.0, 2.0, 3.0, 0.0]
    with pytest.warns(UserWarning, match=r"gave no value for 1 of 4 cells"):
        vor.get_raster_vals_at_centroids([sentinel], ["z"], method="centroid", nodata=0.0)


def test_a_layer_surface_does_not_warn_about_its_own_gaps(gappy, recwarn):
    """``Surface`` fills gaps itself, so warning per cell there would be noise.

    A raster covering part of the model is normal in the layers pipeline:
    ``fill="propagate"`` inherits the rest from the surface above and the pinch-out
    rules turn what is left into idomain.
    """

    from myflopy.surfaces import Surface

    values = Surface.raster(gappy).values(_real_grid(), method="centroid")
    assert np.isnan(values[3])
    assert [w for w in recwarn if "gave no value" in str(w.message)] == []
