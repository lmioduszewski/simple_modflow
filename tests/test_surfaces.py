from __future__ import annotations

import importlib
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy.surfaces import LayerSurfaces, Surface


class _FakeVor:
    """Minimal grid stand-in exposing what surface sampling needs."""

    def __init__(self, n: int = 3):
        points = [Point(float(i), 0.0) for i in range(n)]
        self.gdf_vorPolys = gpd.GeoDataFrame({"geometry": points}, crs="EPSG:2927")
        self.centroids = ([p.x for p in points], [p.y for p in points])
        self.ncpl = n
        self.crs = "EPSG:2927"


def test_modules_import_without_grass():
    # Importing must never require GRASS (it is a system, not pip, dependency).
    importlib.import_module("myflopy.surfaces")
    importlib.import_module("myflopy.modflow.utils.contour_interp")


def test_surfaces_exported_at_top_level():
    assert mf.Surface is Surface
    assert mf.LayerSurfaces is LayerSurfaces


def test_flat_surface_and_reconcile():
    vor = _FakeVor(3)
    # Middle surface (60) sits ABOVE the top (50) -- reconcile must fix ordering.
    layers = LayerSurfaces([Surface.flat(50), Surface.flat(60), Surface.flat(40)])

    gdf = layers.sample(vor, reconcile=True, min_sep=1, trigger_sep=1)

    assert gdf[0].iloc[0] == 50.0
    assert gdf[1].iloc[0] == 49.0  # lowered to top - min_sep
    assert gdf[2].iloc[0] == 40.0  # already fits, untouched


def test_flat_surface_stays_flat_without_reconcile():
    vor = _FakeVor(3)
    gdf = LayerSurfaces([Surface.flat(50), Surface.flat(40)]).sample(
        vor, reconcile=False
    )
    assert list(gdf[0]) == [50.0, 50.0, 50.0]
    assert list(gdf[1]) == [40.0, 40.0, 40.0]


def test_attach_writes_gdf_topbtm():
    vor = _FakeVor(3)
    LayerSurfaces([Surface.flat(50), Surface.flat(40)]).attach(vor, reconcile=False)
    assert vor.gdf_topbtm is not None
    assert list(vor.gdf_topbtm[0]) == [50.0, 50.0, 50.0]


def test_top_botm_arrays_shape():
    vor = _FakeVor(4)
    top, botm = LayerSurfaces(
        [Surface.flat(50), Surface.flat(40), Surface.flat(30)]
    ).top_botm(vor, reconcile=False)

    assert list(top) == [50.0] * 4
    assert botm.shape == (2, 4)
    assert list(botm[0]) == [40.0] * 4


def test_raster_surface_samples_values(tmp_path):
    import rasterio
    from rasterio.transform import from_origin

    path = tmp_path / "r.tif"
    data = np.array([[10.0, 20.0, 30.0]], dtype="float64")  # one row, three cols
    transform = from_origin(-0.5, 0.5, 1.0, 1.0)  # pixel centers at x=0,1,2 y=0
    with rasterio.open(
        path, "w", driver="GTiff", height=1, width=3, count=1,
        dtype="float64", crs="EPSG:2927", transform=transform,
    ) as dst:
        dst.write(data, 1)

    vor = _FakeVor(3)  # centroids at (0,0), (1,0), (2,0)
    assert list(Surface.raster(path).values(vor)) == [10.0, 20.0, 30.0]


def test_points_surface_interpolates_a_plane():
    vor = _FakeVor(3)  # centroids at x = 0, 1, 2
    # Four corners of a plane z = 10*x.
    surface = Surface.from_points([0, 2, 0, 2], [-1, -1, 1, 1], [0, 20, 0, 20])
    assert np.allclose(surface.values(vor), [0.0, 10.0, 20.0])


def test_from_contours_is_lazy_and_does_not_run_grass():
    surface = Surface.from_contours(
        "contours.gpkg", z="elev", region_raster="region.tif", out="out.tif"
    )
    assert surface.kind == "contours"
    assert surface.out == Path("out.tif")
    assert surface.z == "elev"
