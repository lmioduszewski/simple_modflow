from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point, box

import myflopy as mf
from myflopy.modflow.mf6.grid.surfaces import get_raster_vals_at_centroids
from myflopy.surfaces import LayerSurfaces, Surface


def _write_raster(path, data, *, nodata=None, west=0.0, north=4.0, px=1.0, crs="EPSG:2927"):
    """Write a small north-up GeoTIFF for sampling tests."""
    import rasterio
    from rasterio.transform import from_origin

    data = np.asarray(data, dtype="float64")
    h, w = data.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=h, width=w, count=1, dtype="float64",
        crs=crs, transform=from_origin(west, north, px, px), nodata=nodata,
    ) as dst:
        dst.write(data, 1)


def _poly_vor(polys, crs="EPSG:2927"):
    """Minimal grid stand-in backed by real cell *polygons* (for area sampling)."""
    gdf = gpd.GeoDataFrame({"geometry": list(polys)}, crs=crs)
    cents = gdf.geometry.centroid
    return SimpleNamespace(
        gdf_vorPolys=gdf,
        centroids=(cents.x.to_numpy(), cents.y.to_numpy()),
        crs=crs,
        ncpl=len(gdf),
    )


class _FakeVor:
    """Minimal grid stand-in exposing what surface sampling needs."""

    def __init__(self, n: int = 3):
        points = [Point(float(i), 0.0) for i in range(n)]
        self.gdf_vorPolys = gpd.GeoDataFrame({"geometry": points}, crs="EPSG:2927")
        self.centroids = ([p.x for p in points], [p.y for p in points])
        self.ncpl = n
        self.crs = "EPSG:2927"
        self.gdf_topbtm = None

    def get_disv_gridprops(self) -> dict:
        """Minimal DISV geometry stand-in (enough for Surface.to_disv)."""

        return {
            "ncpl": self.ncpl,
            "nvert": self.ncpl + 1,
            "vertices": [[i, float(i), 0.0] for i in range(self.ncpl + 1)],
            "cell2d": [[i, float(i), 0.0, 1, i] for i in range(self.ncpl)],
        }


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


def test_to_disv_builds_package_spec():
    vor = _FakeVor(3)
    spec = LayerSurfaces(
        [Surface.flat(50), Surface.flat(40), Surface.flat(30)]
    ).to_disv(vor, reconcile=False)

    assert spec.name == "disv"
    assert spec.options["nlay"] == 2          # 3 surfaces -> 2 layers
    assert spec.options["ncpl"] == 3
    assert spec.options["nvert"] == 4
    assert list(spec.options["top"]) == [50.0, 50.0, 50.0]
    assert spec.options["botm"].shape == (2, 3)
    assert vor.gdf_topbtm is not None         # attach=True by default


def test_to_disv_without_attach_leaves_grid_untouched():
    vor = _FakeVor(3)
    LayerSurfaces([Surface.flat(50), Surface.flat(40)]).to_disv(
        vor, attach=False, reconcile=False
    )
    assert vor.gdf_topbtm is None


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


def test_relative_surfaces_offset_from_previous():
    vor = _FakeVor(3)
    top, botm = LayerSurfaces(
        [Surface.flat(100), Surface.constant_thickness(10), Surface.offset_below(20)]
    ).top_botm(vor, reconcile=False)
    assert list(top) == [100.0] * 3
    assert list(botm[0]) == [90.0] * 3   # constant_thickness(10) below top
    assert list(botm[1]) == [70.0] * 3   # offset_below(20) below that


def test_relative_surface_cannot_be_model_top():
    vor = _FakeVor(3)
    with pytest.raises(ValueError, match="model top"):
        LayerSurfaces([Surface.offset_below(10), Surface.flat(0)]).sample(
            vor, reconcile=False
        )


def test_pinch_out_marks_thin_cells_pass_through():
    vor = _FakeVor(3)  # centroids x = 0, 1, 2
    # Bottom is the plane z = 10 - 5x -> thickness below flat(10) is [0, 5, 10].
    layers = LayerSurfaces(
        [Surface.flat(10), Surface.from_points([-1, 3, -1, 3], [-1, -1, 1, 1], [15, -5, 15, -5])]
    )
    top, botm, idomain = layers.top_botm_idomain(
        vor, minimum_thickness=1.0, reconcile=False
    )
    assert idomain.shape == (1, 3)
    # Thin cell (thickness 0) -> -1 (pass-through), not 0 (which would block flow).
    assert idomain[0].tolist() == [-1, 1, 1]


def test_pinch_out_requires_min_sep_below_minimum_thickness():
    vor = _FakeVor(3)
    layers = LayerSurfaces([Surface.flat(10), Surface.flat(5)])
    with pytest.raises(ValueError, match="min_sep"):
        # reconcile on by default; min_sep >= minimum_thickness must be rejected.
        layers.top_botm_idomain(vor, minimum_thickness=1.0, min_sep=2.0)


def test_to_disv_pinch_out_sets_idomain():
    vor = _FakeVor(3)
    layers = LayerSurfaces(
        [Surface.flat(10), Surface.from_points([-1, 3, -1, 3], [-1, -1, 1, 1], [15, -5, 15, -5])]
    )
    spec = layers.to_disv(vor, pinch_out=True, minimum_thickness=1.0, reconcile=False)
    idom = spec.options["idomain"]
    assert idom is not None
    assert idom[0].tolist() == [-1, 1, 1]


def test_to_disv_rejects_idomain_and_pinch_out_together():
    vor = _FakeVor(3)
    layers = LayerSurfaces([Surface.flat(10), Surface.flat(5)])
    with pytest.raises(ValueError, match="both"):
        layers.to_disv(vor, pinch_out=True, idomain=[1], reconcile=False)


def test_raster_fill_propagate_inherits_surface_above(tmp_path):
    import rasterio
    from rasterio.transform import from_origin

    path = tmp_path / "r.tif"
    nodata = -9999.0
    data = np.array([[10.0, nodata, 30.0]], dtype="float64")  # middle cell is nodata
    transform = from_origin(-0.5, 0.5, 1.0, 1.0)  # pixel centers x=0,1,2 y=0
    with rasterio.open(
        path, "w", driver="GTiff", height=1, width=3, count=1,
        dtype="float64", crs="EPSG:2927", transform=transform, nodata=nodata,
    ) as dst:
        dst.write(data, 1)

    vor = _FakeVor(3)
    top, botm = LayerSurfaces(
        [Surface.flat(50), Surface.raster(path, fill="propagate")]
    ).top_botm(vor, reconcile=False)
    # Nodata cell inherits the surface above (50); the rest sample normally.
    assert list(botm[0]) == [10.0, 50.0, 30.0]


def test_fill_propagate_on_model_top_raises(tmp_path):
    import rasterio
    from rasterio.transform import from_origin

    path = tmp_path / "r.tif"
    nodata = -9999.0
    with rasterio.open(
        path, "w", driver="GTiff", height=1, width=3, count=1, dtype="float64",
        crs="EPSG:2927", transform=from_origin(-0.5, 0.5, 1.0, 1.0), nodata=nodata,
    ) as dst:
        dst.write(np.array([[10.0, nodata, 30.0]], dtype="float64"), 1)

    vor = _FakeVor(3)
    with pytest.raises(ValueError, match="model top"):
        LayerSurfaces([Surface.raster(path, fill="propagate"), Surface.flat(0)]).sample(
            vor, reconcile=False
        )


def test_thickness_report_summarizes_layers():
    vor = _FakeVor(3)
    report = LayerSurfaces(
        [Surface.flat(10), Surface.flat(4), Surface.flat(0)]
    ).thickness_report(vor, reconcile=False)
    assert "2 layers" in report
    assert "layer 0" in report and "layer 1" in report


def test_from_contours_is_lazy_and_does_not_run_grass():
    surface = Surface.from_contours(
        "contours.gpkg", z="elev", region_raster="region.tif", out="out.tif"
    )
    assert surface.kind == "contours"
    assert surface.out == Path("out.tif")
    assert surface.z == "elev"


def test_grass_launcher_discovery_picks_latest_main_launcher(tmp_path):
    from myflopy.modflow.utils.contour_interp import _find_grass_launcher

    (tmp_path / "grass83.bat").write_text("")
    (tmp_path / "grass84.bat").write_text("")
    (tmp_path / "python-grass84.bat").write_text("")  # the wrapper, must be skipped

    found = _find_grass_launcher([tmp_path / "missing", tmp_path])
    assert found is not None
    assert found.name == "grass84.bat"


def test_grass_bin_env_var_takes_precedence(tmp_path, monkeypatch):
    from myflopy.modflow.utils.contour_interp import _default_grass_bin

    monkeypatch.setenv("GRASS_BIN", str(tmp_path / "my_grass.bat"))
    assert _default_grass_bin() == tmp_path / "my_grass.bat"


# --- area-weighted raster sampling -------------------------------------------
def test_area_weighted_vs_centroid_sampling(tmp_path):
    path = tmp_path / "r.tif"
    data = np.zeros((4, 4))
    data[1, 1] = 90.0  # single hot pixel at center of the sampled block
    _write_raster(path, data)
    # box(0,1,3,4) covers the 3x3 top-left block; centroid (1.5,2.5) = hot pixel.
    vor = _poly_vor([box(0, 1, 3, 4)])

    area = get_raster_vals_at_centroids(vor, [path], ["z"], method="area")["z"].to_numpy()
    cent = get_raster_vals_at_centroids(vor, [path], ["z"], method="centroid")["z"].to_numpy()

    assert area[0] == pytest.approx(10.0)   # mean of 9 pixels = 90/9
    assert cent[0] == pytest.approx(90.0)   # centroid lands on the hot pixel


def test_area_weighted_default_is_area(tmp_path):
    path = tmp_path / "r.tif"
    data = np.zeros((4, 4))
    data[1, 1] = 90.0
    _write_raster(path, data)
    vor = _poly_vor([box(0, 1, 3, 4)])
    # No method= -> default must be area-weighted, not centroid.
    default = get_raster_vals_at_centroids(vor, [path], ["z"])["z"].to_numpy()
    assert default[0] == pytest.approx(10.0)


def test_area_weighted_excludes_nodata(tmp_path):
    path = tmp_path / "r.tif"
    data = [
        [10, 20, 30, 0],
        [40, -9999, 60, 0],
        [70, 80, 90, 0],
        [0, 0, 0, 0],
    ]
    _write_raster(path, data, nodata=-9999)
    vor = _poly_vor([box(0, 1, 3, 4)])  # 3x3 block; center pixel is nodata

    area = get_raster_vals_at_centroids(vor, [path], ["z"], method="area")["z"].to_numpy()
    cent = get_raster_vals_at_centroids(vor, [path], ["z"], method="centroid")["z"].to_numpy()

    assert area[0] == pytest.approx(50.0)   # mean of the 8 valid pixels
    assert np.isnan(cent[0])                # centroid sits on the nodata pixel


def test_area_sampling_falls_back_to_centroid_for_subpixel_cells(tmp_path):
    path = tmp_path / "r.tif"
    data = np.arange(1, 17, dtype="float64").reshape(4, 4)  # data[0,0] == 1
    _write_raster(path, data)
    # A cell entirely inside one pixel contains no pixel center -> 0 coverage.
    vor = _poly_vor([box(0.05, 3.05, 0.45, 3.45)])  # centroid (0.25,3.25) -> pixel [0,0]

    area = get_raster_vals_at_centroids(vor, [path], ["z"], method="area")["z"].to_numpy()
    assert area[0] == pytest.approx(1.0)  # fell back to centroid value, not NaN


def test_area_sampling_warns_on_large_grid(tmp_path, monkeypatch):
    import myflopy.modflow.mf6.grid.surfaces as gs

    monkeypatch.setattr(gs, "_AREA_SAMPLE_WARN_NCPL", 0)  # force the warning
    path = tmp_path / "r.tif"
    _write_raster(path, np.ones((4, 4)))
    vor = _poly_vor([box(0, 1, 3, 4)])
    with pytest.warns(UserWarning, match="large grid"):
        gs.get_raster_vals_at_centroids(vor, [path], ["z"], method="area")


def test_layersurfaces_threads_sampling_method(tmp_path):
    path = tmp_path / "r.tif"
    data = np.zeros((4, 4))
    data[1, 1] = 90.0
    _write_raster(path, data)
    vor = _poly_vor([box(0, 1, 3, 4)])

    top_area, _ = LayerSurfaces([Surface.raster(path), Surface.flat(-100)]).top_botm(
        vor, reconcile=False, method="area"
    )
    top_cent, _ = LayerSurfaces([Surface.raster(path), Surface.flat(-100)]).top_botm(
        vor, reconcile=False, method="centroid"
    )
    assert top_area[0] == pytest.approx(10.0)
    assert top_cent[0] == pytest.approx(90.0)


# --- per-layer pinch thresholds & policies -----------------------------------
def test_idomain_per_layer_thresholds():
    thickness = np.array([[0.5, 5.0, 10.0], [10.0, 0.5, 10.0]])
    idom = LayerSurfaces._idomain_from_thickness(thickness, 1.0)
    assert idom[0].tolist() == [-1, 1, 1]
    assert idom[1].tolist() == [1, -1, 1]
    # per-layer thresholds: layer 0 strict (2.0), layer 1 lax (0.1)
    idom2 = LayerSurfaces._idomain_from_thickness(thickness, [2.0, 0.1])
    assert idom2[0].tolist() == [-1, 1, 1]
    assert idom2[1].tolist() == [1, 1, 1]


def test_idomain_pinch_policies():
    thickness = np.array([[0.5, 5.0]])
    assert LayerSurfaces._idomain_from_thickness(thickness, 1.0, "passthrough")[0].tolist() == [-1, 1]
    assert LayerSurfaces._idomain_from_thickness(thickness, 1.0, "inactive")[0].tolist() == [0, 1]
    assert LayerSurfaces._idomain_from_thickness(thickness, 1.0, "floor")[0].tolist() == [1, 1]


def test_idomain_per_layer_policies():
    thickness = np.array([[0.5, 5.0], [0.5, 5.0]])
    idom = LayerSurfaces._idomain_from_thickness(thickness, 1.0, ["passthrough", "inactive"])
    assert idom[0].tolist() == [-1, 1]
    assert idom[1].tolist() == [0, 1]


def test_idomain_unknown_policy_raises():
    with pytest.raises(ValueError, match="pinch policy"):
        LayerSurfaces._idomain_from_thickness(np.array([[0.5]]), 1.0, "bogus")


def test_idomain_per_layer_length_mismatch_raises():
    thickness = np.array([[0.5, 1.0], [2.0, 3.0]])  # 2 layers
    with pytest.raises(ValueError, match="2 layers"):
        LayerSurfaces._idomain_from_thickness(thickness, [1.0, 2.0, 3.0])


def test_validate_pinch_invariant_exempts_floor_layers():
    # A "floor" layer with a tiny threshold must NOT raise...
    LayerSurfaces._validate_pinch_invariant(
        [0.05, 1.0], ["floor", "passthrough"], 2,
        {"reconcile": True, "min_sep": 0.1},
    )
    # ...but a pinching layer with threshold <= min_sep must raise.
    with pytest.raises(ValueError, match="layer 1"):
        LayerSurfaces._validate_pinch_invariant(
            [1.0, 0.05], ["passthrough", "passthrough"], 2,
            {"reconcile": True, "min_sep": 0.1},
        )


def test_to_disv_pinch_inactive_policy():
    vor = _FakeVor(3)
    layers = LayerSurfaces(
        [Surface.flat(10), Surface.from_points([-1, 3, -1, 3], [-1, -1, 1, 1], [15, -5, 15, -5])]
    )
    spec = layers.to_disv(
        vor, pinch_out=True, minimum_thickness=1.0, pinch="inactive", reconcile=False
    )
    assert spec.options["idomain"][0].tolist() == [0, 1, 1]


# --- surface algebra & isopach -----------------------------------------------
def _sloping(vor=None):
    # z = 10 - 5x -> [10, 5, 0] at cells x = 0, 1, 2
    return Surface.from_points([-1, 3, -1, 3], [-1, -1, 1, 1], [15, -5, 15, -5])


def test_surface_minimum_and_maximum():
    vor = _FakeVor(3)
    a = _sloping()        # [10, 5, 0]
    b = Surface.flat(4)   # [4, 4, 4]
    assert list(Surface.minimum(a, b).values(vor)) == [4, 4, 0]
    assert list(Surface.maximum(a, b).values(vor)) == [10, 5, 4]


def test_surface_clamp_constants():
    vor = _FakeVor(3)
    clamped = Surface.clamp(_sloping(), lower=2, upper=8).values(vor)  # clamp [10,5,0]
    assert list(clamped) == [8, 5, 2]


def test_surface_clamp_with_relative_bound_uses_previous():
    vor = _FakeVor(3)
    prev = np.array([102.0, 102.0, 102.0])
    # base flat 100, but kept at least 5 below the surface above (102 -> 97).
    s = Surface.clamp(Surface.flat(100), upper=Surface.offset_below(5))
    assert list(s.values(vor, previous=prev)) == [97, 97, 97]


def test_surface_where_selects_by_zone():
    vor = _poly_vor([box(0, 0, 1, 1), box(1, 0, 2, 1)])  # centroids (0.5,.5),(1.5,.5)
    zone = box(-0.1, -0.1, 1.0, 1.1)  # contains the first centroid only
    out = Surface.where(zone, inside=Surface.flat(100), outside=Surface.flat(0)).values(vor)
    assert list(out) == [100, 0]


def test_surface_isopach_subtracts_thickness_from_previous():
    vor = _FakeVor(3)
    prev = np.array([100.0, 100.0, 100.0])
    assert list(Surface.isopach(Surface.flat(30)).values(vor, previous=prev)) == [70, 70, 70]


def test_isopach_cannot_be_model_top():
    vor = _FakeVor(3)
    with pytest.raises(ValueError, match="isopach cannot be the model top"):
        Surface.isopach(Surface.flat(30)).values(vor, previous=None)


# --- fluent surface algebra (capped_at / floored_at / between / offsets) ------
def test_fluent_capped_and_floored_match_min_max():
    vor = _FakeVor(3)
    a = _sloping()        # [10, 5, 0]
    b = Surface.flat(4)   # [4, 4, 4]
    assert list(a.capped_at(b).values(vor)) == [4, 4, 0]     # never above b == minimum
    assert list(a.floored_at(b).values(vor)) == [10, 5, 4]   # never below b == maximum


def test_fluent_between_matches_clamp():
    vor = _FakeVor(3)
    assert list(_sloping().between(lower=2, upper=8).values(vor)) == [8, 5, 2]


def test_fluent_offsets_and_arithmetic():
    vor = _FakeVor(3)
    assert list(Surface.flat(100).below(5).values(vor)) == [95, 95, 95]
    assert list(Surface.flat(100).above(5).values(vor)) == [105, 105, 105]
    assert list((Surface.flat(100) - 5).values(vor)) == [95, 95, 95]
    assert list((Surface.flat(100) + 5).values(vor)) == [105, 105, 105]
    assert list((5 + Surface.flat(100)).values(vor)) == [105, 105, 105]   # __radd__


def test_fluent_within_matches_where():
    vor = _poly_vor([box(0, 0, 1, 1), box(1, 0, 2, 1)])
    zone = box(-0.1, -0.1, 1.0, 1.1)
    out = Surface.flat(100).within(zone, outside=Surface.flat(0)).values(vor)
    assert list(out) == [100, 0]


def test_fluent_capped_at_shifted_surface():
    # "bedrock capped at 5 ft below ground" -- the readable form of the notebook's
    # Clamp(bedrock, upper=ground - 5).
    vor = _FakeVor(3)
    ground = Surface.flat(102)
    s = Surface.flat(100).capped_at(ground - 5)              # min(100, 97) -> 97
    assert list(s.values(vor)) == [97, 97, 97]


def test_subtracting_a_surface_from_a_surface_is_unsupported():
    with pytest.raises(TypeError):
        Surface.flat(100) - Surface.flat(5)
