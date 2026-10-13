"""No NaN layer elevation reaches a DISV built from a layer stack (ledger 587).

A raster covering only part of the grid, used as a unit bottom or as the stack's top,
leaves its cells ``NaN`` -- ``Surface.values`` samples with ``gaps="ignore"``, since a
``fill="propagate"`` unit fills them from the surface above. Nothing refused the rest:
``stack.to_disv``, ``mf.disv(stack)`` (eager and lazy), ``mf.disv(result)`` and
``LayerSurfaces.to_disv`` wrote the ``NaN`` into the DISV, and MF6 ran it to normal
termination with NaN heads (``docs/reproducers/nan_layer_elevations.py``; ledger 555
measured 200 of 400 such bottoms and refused only a basin's). Every route goes through
``surfaces.refuse_unbounded_cells`` now, which names the elevations, the cells and the
source; a fill (``fill="propagate"``, ``mf.Max``) passes, and so does a ``NaN`` in an
INACTIVE cell, which MF6 never uses.
"""

from __future__ import annotations

import dataclasses
import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.surfaces import LayerSurfaces

rasterio = pytest.importorskip("rasterio")


@pytest.fixture(scope="module")
def vor():
    return irregular_voronoi_grid(nrow=4, ncol=4, cell_size=100.0)


def _west_half(path):
    """A raster of 60.0 covering x < 200 of the 400 ft square, no-data east of it."""

    from rasterio.transform import from_origin

    with rasterio.open(path, "w", driver="GTiff", height=42, width=20, count=1,
                       dtype="float64", crs="EPSG:2927",
                       transform=from_origin(0.0, 410.0, 10.0, 10.0), nodata=-9999.0) as dst:
        dst.write(np.full((1, 42, 20), 60.0))
    return path


@pytest.fixture
def half(tmp_path):
    return _west_half(tmp_path / "half.tif")


def _east(vor):
    return np.flatnonzero(vor.gdf_vorPolys.centroid.x.to_numpy() > 200.0)


def test_a_partial_raster_as_a_unit_bottom_is_refused_naming_it(vor, half):
    stack = mf.LayerStack(vor, top=mf.Flat(100.0)).add("sand", bottom=mf.Raster(half)).add(
        "clay", thickness=10.0)
    east = _east(vor)
    with pytest.raises(ValueError) as refused:
        stack.to_disv()
    message = str(refused.value)
    assert f"the bottom of layer 'sand' (raster 'half.tif'): {east.size} of {vor.ncpl} cells" in message
    assert str(east[:10].tolist()) in message
    assert "the bottom of layer 'clay' (measured from the surface above" in message
    assert "fill='propagate'" in message
    for route in (lambda: mf.disv(stack), lambda: mf.disv(stack.build()),
                  lambda: stack.build().to_disv()):
        with pytest.raises(ValueError, match=r"raster 'half\.tif'"):
            route()


def test_a_partial_raster_as_the_top_is_refused(vor, half):
    stack = mf.LayerStack(vor, top=mf.Raster(half)).add("sand", thickness=20.0)
    east = _east(vor).size
    with pytest.raises(ValueError, match=rf"the model top \(raster 'half\.tif'\): {east} of 16 cells"):
        mf.disv(stack)


def test_the_engine_refuses_too(vor, half):
    surfaces = LayerSurfaces([mf.Flat(100.0), mf.Raster(half)], labels=["top", "sand"])
    with pytest.raises(ValueError, match=r"the bottom of layer 'sand' \(raster 'half\.tif'\)"):
        surfaces.to_disv(vor, attach=False)


def test_a_fill_passes(vor, half):
    propagated = mf.LayerStack(vor, top=mf.Flat(100.0)).add(
        "sand", bottom=mf.Raster(half, fill="propagate")).add("clay", thickness=10.0)
    fallback = mf.LayerStack(vor, top=mf.Max(mf.Raster(half), mf.Flat(70.0))).add(
        "sand", thickness=20.0)
    for stack in (propagated, fallback):
        spec = mf.disv(stack)
        assert np.isfinite(spec.options["top"]).all() and np.isfinite(spec.options["botm"]).all()


def test_a_nan_in_an_inactive_cell_passes(vor, half):
    """MF6 reads it and never uses it: every head as with a finite value there."""

    result = mf.LayerStack(vor, top=mf.Flat(100.0)).add("sand", bottom=mf.Raster(half)).build()
    idomain = result.idomain.copy()
    idomain[:, _east(vor)] = 0
    spec = dataclasses.replace(result, idomain=idomain).to_disv()
    assert np.isnan(spec.options["botm"]).sum() == _east(vor).size


def test_a_lazy_disv_is_refused_at_build(tmp_path, half):
    path = tmp_path / "dom.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 400, 400)], crs=2927).to_file(path)
    model = mf.gwf("m", packages=[
        mf.disv(mf.LayerStack(top=mf.Flat(100.0)).add("sand", bottom=mf.Raster(half))),
        mf.npf(k=1.0), mf.ic(strt=90.0), mf.oc(head_filerecord="m.hds"),
    ]).with_grid(mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=path),
                                     max_area=20_000.0))
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match=r"the bottom of layer 'sand' \(raster 'half\.tif'\)"):
            mf.Project(tmp_path / "p", name="p").prepare_run("r", simulation)


def test_a_split_unit_names_its_raster(vor, half):
    """A split unit's slices are ``toward`` surfaces: the refusal names the raster under
    them, where it said only "a 'toward' surface" (review of ledger 587)."""

    stack = mf.LayerStack(vor, top=mf.Flat(100.0)).add("sand", bottom=mf.Raster(half), split=2)
    with pytest.raises(ValueError) as refused:
        stack.to_disv()
    message = str(refused.value)
    for name in ("sand_1", "sand_2"):
        assert f"the bottom of layer {name!r} (a 'toward' surface of raster 'half.tif')" in message


def test_qc_counts_what_to_disv_refuses_a_pass_through_cell_included(vor, half):
    """A pass-through cell bounded by NaN is refused by ``to_disv`` (MF6 disconnects the
    cell below it in silence), and ``qc()`` -- which the refusal points at -- counted
    only ``idomain == 1``, so it read ok (review of ledger 587)."""

    result = mf.LayerStack(vor, top=mf.Flat(100.0)).add("sand", bottom=mf.Raster(half)).add(
        "clay", thickness=10.0).build()
    east = _east(vor)
    idomain = result.idomain.copy()
    idomain[0, east] = -1   # sand: pass-through over the gap
    idomain[1, east] = 0    # clay: inactive, so only the pass-through cells count
    passing = dataclasses.replace(result, idomain=idomain)
    report = passing.qc()
    assert report.nan_active_cells == east.size and not report.ok
    with pytest.raises(ValueError, match=r"the bottom of layer 'sand'"):
        passing.to_disv()
