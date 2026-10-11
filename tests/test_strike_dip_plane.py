"""A strike/dip surface is the geological plane, sampled at pixel centres (ledger 362).

Convention: azimuths in degrees clockwise from north (+y); the plane dips toward
``strike + 90`` (the right-hand rule), descends at ``tan(dip)`` per unit distance in
that direction, is level along the strike, and passes through the known point.
``get_normal_from_strike_and_dip`` is its UPWARD UNIT normal.

Before: the raster fell toward azimuth ``strike`` (level along the dip direction)
by ``sin(2 * dip) / 2`` per unit distance -- ``dip=45`` gave 0.5 for ``tan`` 1.0 and
``dip=80`` gave the same slope as ``dip=10`` -- the "normal" of a level plane was
horizontal, every pixel held the plane at its upper-left corner, and the default
output file was ``raster.tif`` in the directory current when myflopy was IMPORTED.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import xy

from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid import get_normal_from_strike_and_dip
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

STRIKES = (0.0, 37.0, 90.0, 200.0, 315.0)
DIPS = (0.0, 10.0, 45.0, 80.0)
KNOWN = (130.0, 270.0, 500.0)


@pytest.fixture(scope="module")
def grid():
    return irregular_voronoi_grid(nrow=4, ncol=4, cell_size=100.0, refined=False)


def _plane(strike, dip, x, y, known=KNOWN):
    """The geological plane, written independently of the code under test."""
    s, d = np.radians(strike), np.radians(dip)
    dip_dir = np.array([np.sin(s + np.pi / 2), np.cos(s + np.pi / 2)])  # strike + 90
    along_dip = (np.asarray(x) - known[0]) * dip_dir[0] + (np.asarray(y) - known[1]) * dip_dir[1]
    return known[2] - np.tan(d) * along_dip


def _read_with_centres(path):
    with rasterio.open(path) as src:
        values = src.read(1).astype(float)
        rows, cols = np.indices(values.shape)
        xs, ys = xy(src.transform, rows.ravel(), cols.ravel(), offset="center")
    return values.ravel(), np.asarray(xs), np.asarray(ys)


@pytest.mark.parametrize(("strike", "dip"), list(itertools.product(STRIKES, DIPS)))
def test_the_raster_is_the_strike_dip_plane_at_each_pixel_centre(grid, tmp_path, strike, dip):
    out = tmp_path / "plane.tif"
    grid.get_raster_from_strike_dip(strike, dip, KNOWN, pixel_size=10, output_filename=out)
    values, xs, ys = _read_with_centres(out)
    expected = _plane(strike, dip, xs, ys)
    # float32 band: relative 1e-6 of the largest magnitude in the raster.
    tol = 1e-6 * max(1.0, np.abs(expected).max())
    np.testing.assert_allclose(values, expected, atol=tol, rtol=0)

    # A least-squares fit of the raster recovers tan(dip) toward strike + 90 and a
    # level line along the strike.
    coef = np.linalg.lstsq(np.c_[np.ones_like(xs), xs, ys], values, rcond=None)[0]
    s = np.radians(strike)
    gradient = np.array(coef[1:])
    strike_dir = np.array([np.sin(s), np.cos(s)])
    dip_dir = np.array([np.cos(s), -np.sin(s)])
    assert gradient @ strike_dir == pytest.approx(0.0, abs=1e-5)
    assert -(gradient @ dip_dir) == pytest.approx(np.tan(np.radians(dip)), rel=1e-5, abs=1e-6)


def test_a_cell_sample_carries_no_half_pixel_offset(grid, tmp_path):
    # Area-weighted means of a linear plane, compared with the plane at each cell
    # centroid: individual cells differ by how pixels tile them, but the BIAS is
    # zero. Corner registration put it at about half a pixel of slope (2.5 ft here).
    strike, dip, pixel = 0.0, 45.0, 5.0
    frame = grid.get_raster_from_strike_dip(
        strike, dip, KNOWN, pixel_size=pixel, output_filename=tmp_path / "p.tif"
    )
    cx, cy = np.asarray(grid.centroids[0]), np.asarray(grid.centroids[1])
    bias = float(np.mean(frame["elev"].to_numpy() - _plane(strike, dip, cx, cy)))
    half_pixel_of_slope = np.tan(np.radians(dip)) * pixel / 2
    assert abs(bias) < 0.1 * half_pixel_of_slope


@pytest.mark.parametrize(("strike", "dip"), list(itertools.product(STRIKES, DIPS)))
def test_the_normal_is_the_planes_upward_unit_normal(strike, dip):
    n = get_normal_from_strike_and_dip(strike, dip)
    s, d = np.radians(strike), np.radians(dip)
    strike_line = np.array([np.sin(s), np.cos(s), 0.0])
    dip_line = np.array([np.cos(s) * np.cos(d), -np.sin(s) * np.cos(d), -np.sin(d)])
    assert np.linalg.norm(n) == pytest.approx(1.0, abs=1e-12)
    assert n @ strike_line == pytest.approx(0.0, abs=1e-12)
    assert n @ dip_line == pytest.approx(0.0, abs=1e-12)
    assert n[2] > 0
    # The bound method is the same function.
    np.testing.assert_array_equal(VoronoiGridPlus.get_normal_from_strike_and_dip(strike, dip), n)


def test_a_level_plane_has_a_vertical_normal():
    np.testing.assert_allclose(get_normal_from_strike_and_dip(0, 0), [0.0, 0.0, 1.0], atol=1e-15)
    np.testing.assert_allclose(get_normal_from_strike_and_dip(123, 0), [0.0, 0.0, 1.0], atol=1e-15)


@pytest.mark.parametrize("dip", [90, 90.0, 95.0, 180.0])
def test_a_vertical_or_overturned_dip_is_refused(grid, tmp_path, dip):
    with pytest.raises(ValueError, match=r"not in \[0, 90\)"):
        get_normal_from_strike_and_dip(10.0, dip)
    with pytest.raises(ValueError, match=r"not in \[0, 90\)"):
        grid.get_raster_from_strike_dip(10.0, dip, KNOWN, output_filename=tmp_path / "x.tif")
    assert not (tmp_path / "x.tif").exists()


def test_a_negative_dip_names_the_strike_that_spells_it():
    with pytest.raises(ValueError, match=r"strike=190\.0, dip=20\.0"):
        get_normal_from_strike_and_dip(10.0, -20.0)


@pytest.mark.parametrize(("strike", "dip", "error"), [
    (float("nan"), 10.0, ValueError),
    (10.0, float("inf"), ValueError),
    ("10", 10.0, TypeError),
    (10.0, True, TypeError),
])
def test_a_strike_or_dip_that_is_not_a_finite_number_is_refused(strike, dip, error):
    with pytest.raises(error):
        get_normal_from_strike_and_dip(strike, dip)


@pytest.mark.parametrize("bound", [True, False], ids=["method", "function"])
def test_the_default_output_file_is_in_the_working_directory_at_the_call(
    grid, tmp_path, monkeypatch, bound
):
    from myflopy.modflow.mf6.grid import get_raster_from_strike_dip

    def write():
        if bound:
            grid.get_raster_from_strike_dip(0.0, 10.0, KNOWN, pixel_size=20)
        else:
            get_raster_from_strike_dip(grid, 0.0, 10.0, KNOWN, pixel_size=20)

    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()
    monkeypatch.chdir(first)
    write()
    monkeypatch.chdir(second)
    write()
    assert (first / "raster.tif").exists()
    assert (second / "raster.tif").exists()


def test_a_str_output_filename_is_written_and_sampled(grid, tmp_path):
    out = str(tmp_path / "plane.tif")
    frame = grid.get_raster_from_strike_dip(0.0, 10.0, KNOWN, pixel_size=20, output_filename=out)
    assert Path(out).exists()
    assert np.isfinite(frame["elev"].to_numpy()).all()
