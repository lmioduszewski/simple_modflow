"""`mf.Contours` interpolates where its contours ARE (ledger 602).

GRASS works in ONE frame, its location's CRS, and myflopy imports with ``-o``, so
GRASS never compared the two: contours in EPSG:4326 were interpolated at their raw
degree coordinates inside an EPSG:2927 location -- 0 of 6,000 pixels, no warning --
and a region polygon in another CRS the same. The contours and a region polygon are
now reprojected into the location's CRS and staged as the copies GRASS imports; a
file carrying NO CRS is read as already in that frame (the user's rule, 2026-10-16);
and what cannot be reconciled is refused before GRASS starts.

Every test but the last runs a FAKE GRASS (`_grass_modules` patched) that records
what it was asked to import, so it reads the staged file itself; the last runs the
real one and is skipped without it. `grassdata=` is under `tmp_path` everywhere.
"""

from __future__ import annotations

import os
import shutil
import types
import warnings

import geopandas as gpd
import numpy as np
import pytest
import rasterio
import shapely
import shapely.affinity
from pyproj import CRS
from rasterio.transform import from_origin

from myflopy.modflow.utils import contour_interp
from myflopy.modflow.utils.derived_raster import DerivedRaster
from myflopy.surfaces import Surface

#: Horizontal-ish contours over a 1000 x 600 ft box at the origin (EPSG:2927).
_LINES = [shapely.LineString([(0, y), (1000, y + 20)]) for y in range(50, 600, 50)]
_ELEV = [100.0 - y / 10 for y in range(50, 600, 50)]


class _FakeGrass:
    """Stands in for GRASS: records what each `v.in_ogr`/`r.in_gdal` was given.

    A staged file lives in a temporary directory removed when `run()` returns, so
    the import is READ here, at the call, not merely remembered by path.
    """

    def __init__(self):
        self.started = False
        self.imported: dict[str, gpd.GeoDataFrame] = {}
        self.rasters: dict[str, str] = {}

    def modules(self, grass_bin=None):
        self.started = True
        session = types.SimpleNamespace(finish=lambda: None)

        def in_ogr(*, input, output, **_):
            self.imported[output] = gpd.read_file(input)

        def in_gdal(*, input, output, **_):
            self.rasters[output] = input

        def out_gdal(*, output, **_):
            with open(output, "w") as handle:
                handle.write("fake raster")

        noop = lambda **_: None  # noqa: E731
        r = types.SimpleNamespace(in_gdal=in_gdal, surf_contour=noop, mapcalc=noop, out_gdal=out_gdal)
        g = types.SimpleNamespace(region=noop)
        v = types.SimpleNamespace(in_ogr=in_ogr, to_rast=noop)
        gsetup = types.SimpleNamespace(init=lambda *a, **k: session)
        return r, g, v, gsetup


@pytest.fixture
def grass(monkeypatch):
    fake = _FakeGrass()
    monkeypatch.setattr(contour_interp, "_grass_modules", fake.modules)
    return fake


@pytest.fixture
def files(tmp_path):
    """The contours in 2927, a 4326 copy, a no-CRS copy, and the box region in both."""

    contours = gpd.GeoDataFrame({"elev": _ELEV}, geometry=_LINES, crs=2927)
    region = gpd.GeoDataFrame(geometry=[shapely.box(0, 0, 1000, 600)], crs=2927)
    paths = types.SimpleNamespace(
        c2927=tmp_path / "c2927.gpkg", c4326=tmp_path / "c4326.gpkg",
        cnone=tmp_path / "cnone.gpkg", r2927=tmp_path / "r2927.gpkg",
        r4326=tmp_path / "r4326.gpkg", grassdata=tmp_path / "grassdata",
        tmp=tmp_path,
    )
    contours.to_file(paths.c2927)
    contours.to_crs(4326).to_file(paths.c4326)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)   # pyogrio: "'crs' was not provided"
        gpd.GeoDataFrame({"elev": _ELEV}, geometry=_LINES).to_file(paths.cnone)
    region.to_file(paths.r2927)
    region.to_crs(4326).to_file(paths.r4326)
    return paths


def _surface(files, contours, *, region_vector=None, region_raster=None, z="elev", **kwargs):
    if region_vector is None and region_raster is None:
        region_vector = files.r2927
    return Surface.from_contours(
        contours, z=z, region_vector=region_vector, region_raster=region_raster,
        resolution=10, out=files.tmp / f"{contours.stem}.interp.tif",
        grassdata=files.grassdata, **kwargs,
    )


def _location(files, name="myflopy_interp", **projection_files):
    permanent = files.grassdata / name / "PERMANENT"
    permanent.mkdir(parents=True)
    for file_name, text in projection_files.items():
        (permanent / file_name).write_text(text)


def _raster(path, *, crs, west=0.0, north=600.0, size=10.0):
    with rasterio.open(
        path, "w", driver="GTiff", height=60, width=100, count=1, dtype="float32",
        crs=crs, transform=from_origin(west, north, size, size),
    ) as handle:
        handle.write(np.ones((1, 60, 100), dtype="float32"))


# --------------------------------------------------------------- reprojection


def test_contours_in_another_crs_reach_grass_in_the_location_frame(grass, files):
    _surface(files, files.c4326).resolve_source()

    imported = grass.imported["vectContours"]
    assert imported.crs.to_epsg() == 2927
    np.testing.assert_allclose(imported.total_bounds, (0, 50, 1000, 570), atol=1e-3)
    assert imported["elev"].tolist() == _ELEV


def test_a_region_vector_in_another_crs_reaches_grass_in_the_location_frame(grass, files):
    _surface(files, files.c2927, region_vector=files.r4326).resolve_source()

    region = grass.imported["VRegion"]
    assert region.crs.to_epsg() == 2927
    np.testing.assert_allclose(region.total_bounds, (0, 0, 1000, 600), atol=1e-3)


def test_a_contour_file_without_a_crs_is_read_in_the_location_frame(grass, files):
    """The user's rule (2026-10-16): no CRS is the location's frame, never refused."""

    surface = _surface(files, files.cnone)
    surface.resolve_source()

    imported = grass.imported["vectContours"]
    assert imported.crs.to_epsg() == 2927
    np.testing.assert_allclose(imported.total_bounds, (0, 50, 1000, 570))
    assert "reprojected" not in surface._derived_raster().params


def test_the_staged_copy_of_a_file_in_the_frame_holds_the_same_features(grass, files):
    _surface(files, files.c2927).resolve_source()

    imported = grass.imported["vectContours"]
    original = gpd.read_file(files.c2927)
    assert imported.geometry.geom_equals_exact(original.geometry, tolerance=0).all()


def test_a_shapefile_whose_fid_column_is_not_a_key_is_staged(grass, files):
    """GeoPackage takes a column named `fid` as its key (ledger 604).

    QGIS's "Merge vector layers" writes a shapefile whose `fid` repeats; GRASS
    imported it as it was, and the staged copy refused it ("UNIQUE constraint
    failed"). The copy carries the geometry and `z` alone.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        gpd.GeoDataFrame({"fid": [1] * len(_ELEV), "elev": _ELEV, "note": "x"},
                         geometry=_LINES, crs=2927).to_crs(4326).to_file(files.tmp / "merged.shp")
        gpd.GeoDataFrame({"fid": [1.5]}, geometry=[shapely.box(0, 0, 1000, 600)],
                         crs=2927).to_file(files.tmp / "region.shp")

    _surface(files, files.tmp / "merged.shp",
             region_vector=files.tmp / "region.shp").resolve_source()

    contours = grass.imported["vectContours"]
    assert list(contours.columns) == ["elev", "geometry"]
    assert contours["elev"].tolist() == _ELEV
    assert list(grass.imported["VRegion"].columns) == ["geometry"]


def test_z_is_matched_ignoring_case_as_grass_matches_it(grass, files):
    """GRASS's SQLite table reads `ELEV` for a column `elev` (measured on GRASS 8.4)."""

    _surface(files, files.c2927, z="ELEV").resolve_source()

    assert list(grass.imported["vectContours"].columns) == ["elev", "geometry"]


def test_contours_without_the_z_attribute_are_refused_before_grass(grass, files):
    with pytest.raises(ValueError, match=r"carry no attribute 'height'.*They carry: 'elev'"):
        _surface(files, files.c2927, z="height").resolve_source()
    assert not grass.started


def test_a_crs_no_transformation_reaches_is_refused_naming_set_crs(grass, files):
    """An engineering CRS (a CAD export's) cannot be carried into 2927 at all."""

    local = ('LOCAL_CS["Unknown",LOCAL_DATUM["Unknown",32767],UNIT["foot",0.3048],'
             'AXIS["Easting",EAST],AXIS["Northing",NORTH]]')
    gpd.GeoDataFrame({"elev": _ELEV}, geometry=_LINES, crs=local).to_file(files.tmp / "cad.gpkg")

    with pytest.raises(ValueError, match=r"contours 'cad.gpkg' is in Unknown, which cannot "
                                         r"be carried into the GRASS location's EPSG:2927.*set_crs"):
        _surface(files, files.tmp / "cad.gpkg").resolve_source()
    assert not grass.started


# ------------------------------------------------------------------ refusals


@pytest.mark.parametrize("projection", [
    {"PROJ_SRID": "EPSG:2927\n"},
    {"PROJ_WKT": CRS(2927).to_wkt()},
    {"PROJ_EPSG": "epsg: 2927\n"},                       # GRASS 7's spelling
], ids=["srid", "wkt", "grass7"])
def test_epsg_contradicting_an_existing_location_is_refused_before_grass(grass, files, projection):
    _location(files, **projection)
    surface = _surface(files, files.c2927, epsg="32610")

    with pytest.raises(ValueError, match=r"location='myflopy_interp_32610'") as caught:
        surface.resolve_source()

    assert "already exists in EPSG:2927" in str(caught.value)
    assert "epsg='2927'" in str(caught.value)
    assert not grass.started


def test_a_contradicted_location_refuses_even_a_fresh_cache(grass, files):
    """The frame is checked BEFORE the cache: a raster GRASS wrote in the location's CRS
    while `epsg=` said another is never handed back as fresh."""

    surface = _surface(files, files.c2927, epsg="32610")
    surface.resolve_source()                       # no location yet: created in 32610
    assert surface.cache_status() == "fresh"
    grass.started = False
    _location(files, PROJ_SRID="EPSG:2927\n")      # ...but the location says otherwise

    with pytest.raises(ValueError, match="location="):
        surface.resolve_source()
    assert not grass.started


def test_an_existing_location_in_the_asked_crs_is_used(grass, files):
    _location(files, PROJ_SRID="EPSG:2927\n")
    _surface(files, files.c2927).resolve_source()
    assert grass.started


def test_contours_reaching_no_part_of_the_region_are_refused_before_grass(grass, files):
    far = gpd.GeoDataFrame({"elev": _ELEV}, geometry=[shapely.affinity.translate(g, 1e6, 0)
                                                      for g in _LINES], crs=2927)
    far.to_file(files.tmp / "far.gpkg")

    with pytest.raises(ValueError, match="None of the 11 contours in 'far.gpkg'") as caught:
        _surface(files, files.tmp / "far.gpkg").resolve_source()

    message = str(caught.value)
    assert "(1e+06, 50, 1.001e+06, 570)" in message and "(0, 0, 1,000, 600)" in message
    assert "EPSG:2927" in message and "carries no CRS" not in message
    assert not grass.started


def test_a_contour_file_without_a_crs_far_from_the_region_is_refused(grass, files):
    """No-CRS contours in another frame are caught by the zero-overlap guard (decision 1)."""

    degrees = gpd.read_file(files.c4326).set_crs(None, allow_override=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        degrees.to_file(files.tmp / "degrees_nocrs.gpkg")

    with pytest.raises(ValueError, match="carries no CRS and was read as EPSG:2927"):
        _surface(files, files.tmp / "degrees_nocrs.gpkg").resolve_source()
    assert not grass.started


def test_the_zero_overlap_refusal_names_a_moved_files_own_crs(grass, files):
    far = gpd.GeoDataFrame({"elev": _ELEV}, geometry=[shapely.affinity.translate(g, 1e6, 0)
                                                      for g in _LINES], crs=2927).to_crs(4326)
    far.to_file(files.tmp / "far4326.gpkg")

    with pytest.raises(ValueError, match="The contours' file is in EPSG:4326"):
        _surface(files, files.tmp / "far4326.gpkg").resolve_source()


def test_a_region_raster_in_another_crs_is_refused(grass, files):
    _raster(files.tmp / "dem4326.tif", crs="EPSG:4326", west=-123.0, north=46.0, size=1e-4)

    with pytest.raises(ValueError, match=r"Pass region_vector=") as caught:
        _surface(files, files.c2927, region_raster=files.tmp / "dem4326.tif").resolve_source()

    assert "'dem4326.tif' is in EPSG:4326" in str(caught.value)
    assert not grass.started


@pytest.mark.parametrize("crs", ["EPSG:2927", None], ids=["in-frame", "no-crs"])
def test_a_region_raster_in_the_frame_or_without_a_crs_is_read_as_it_is(grass, files, crs):
    _raster(files.tmp / "dem.tif", crs=crs)
    _surface(files, files.c2927, region_raster=files.tmp / "dem.tif").resolve_source()

    assert grass.rasters["Region"] == str(files.tmp / "dem.tif")
    assert "VRegion" not in grass.imported


def test_missing_files_are_not_refused_by_the_frame(files):
    """`cache_status()` is asked about files not yet written: nothing is read."""

    surface = _surface(files, files.tmp / "nope.gpkg", region_vector=files.tmp / "nope2.gpkg")
    assert surface.cache_status() == "missing"
    surface._contour_frame().refuse()


# ------------------------------------------------- one frame, spelled twice

#: EPSG:2927 as the user's lidar DEM and ground-surface contours carry it: no
#: authority code and a false easting of 1640416.66666667 where EPSG's WKT reads
#: 1640416.667. pyproj's `equals` says it is not EPSG:2927; reprojecting moves a
#: point 0.00033 ft (ledger 603).
_SPELLED_2927 = CRS(2927).to_wkt("WKT1_GDAL").replace("1640416.667", "1640416.66666667").replace(
    ',AUTHORITY["EPSG","2927"]]', "]")


def test_the_spelled_2927_is_not_equal_but_is_one_frame():
    """The premise of the three tests below, so none of them passes for the wrong reason."""

    from myflopy.areas import _same_crs

    assert not _same_crs(_SPELLED_2927, "EPSG:2927")
    assert contour_interp._same_frame(_SPELLED_2927, "EPSG:2927")
    # A real datum change still moves feet, and degrees cannot be carried at all.
    assert not contour_interp._same_frame("EPSG:2286", "EPSG:2927")
    assert not contour_interp._same_frame("EPSG:4326", "EPSG:2927")


def test_a_contour_file_in_the_locations_crs_spelled_otherwise_is_not_moved(grass, files):
    spelled = files.tmp / "spelled.gpkg"
    gpd.GeoDataFrame({"elev": _ELEV}, geometry=_LINES, crs=_SPELLED_2927).to_file(spelled)
    surface = _surface(files, spelled)

    assert surface._derived_raster().params == _pre_fix_params(files)    # no cache moves
    surface.resolve_source()
    imported = grass.imported["vectContours"]
    assert imported.crs.to_epsg() == 2927
    assert imported.geometry.geom_equals_exact(gpd.GeoSeries(_LINES, crs=2927), tolerance=0).all()


def test_a_region_raster_in_the_locations_crs_spelled_otherwise_is_read(grass, files):
    _raster(files.tmp / "lidar.tif", crs=_SPELLED_2927)
    _surface(files, files.c2927, region_raster=files.tmp / "lidar.tif").resolve_source()

    assert grass.rasters["Region"] == str(files.tmp / "lidar.tif")


def test_a_location_with_no_epsg_code_is_refused_naming_location_only(grass, files):
    _location(files, PROJ_WKT=CRS("+proj=tmerc +lat_0=0 +lon_0=-121 +k=1 +x_0=0 +y_0=0 "
                                  "+ellps=GRS80 +units=m").to_wkt())
    with pytest.raises(ValueError, match="location='myflopy_interp_2927'") as caught:
        _surface(files, files.c2927).resolve_source()
    assert "or epsg=" not in str(caught.value)


def test_given_both_regions_the_interpolator_checks_the_raster_it_uses(grass, files):
    """`_set_region` reads the raster first; the overlap check and staging follow it."""

    _raster(files.tmp / "dem.tif", crs="EPSG:2927")
    far = gpd.GeoDataFrame(geometry=[shapely.box(1e6, 0, 1.001e6, 600)], crs=2927)
    far.to_file(files.tmp / "far_region.gpkg")
    contour_interp.ContourSurfaceInterpolator(
        files.c2927, out=files.tmp / "both.tif", region_raster=files.tmp / "dem.tif",
        region_vector=files.tmp / "far_region.gpkg", z_field="elev", resolution=10,
        grassdata=files.grassdata,
    ).run()

    assert grass.rasters["Region"] == str(files.tmp / "dem.tif")
    assert "VRegion" not in grass.imported


# --------------------------------------------------------------------- cache


def _pre_fix_params(files):
    """What a surface's cache key held before reprojection, verbatim."""

    return {"z": "elev", "resolution": 10, "epsg": "2927", "clip": True,
            "grass_kwargs": {"grassdata": str(files.grassdata)}}


def test_a_same_crs_surface_keeps_its_cache_key(files):
    """No correct cache moves: the key of a surface whose files never move is the old one."""

    for contours in (files.c2927, files.cnone):
        assert _surface(files, contours)._derived_raster().params == _pre_fix_params(files)


def test_a_moved_input_is_recorded_in_the_cache_key(files):
    params = _surface(files, files.c4326, region_vector=files.r4326)._derived_raster().params
    assert params["reprojected"] == {"contours": "EPSG:4326", "region_vector": "EPSG:4326"}


def test_a_cache_interpolated_from_an_unreprojected_file_is_refused_not_reused(grass, files):
    surface = _surface(files, files.c4326)
    old = DerivedRaster(surface.out, [files.c4326, files.r2927], _pre_fix_params(files),
                        lambda: surface.out.write_text("blank, made at degree coordinates"))
    old.ensure()                                   # the sidecar the old code wrote

    with pytest.raises(ValueError, match=r"resolve_source\(refresh=True\)") as caught:
        surface.resolve_source()
    assert "contours in EPSG:4326 reached GRASS read as EPSG:2927" in str(caught.value)
    with pytest.raises(ValueError, match="not reused"):
        surface.resolve_source(warn=False)         # `warn=False` does not silence it
    assert not grass.started

    surface.resolve_source(refresh=True)
    assert grass.started
    assert surface.out.read_text() == "fake raster"
    assert surface.cache_status() == "fresh"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        surface.resolve_source()                   # the rebuilt cache is fresh


# ---------------------------------------------------------------- real GRASS


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("grass") is None and not os.environ.get("GRASS_BIN"),
                    reason="mf.Contours needs GRASS GIS")
def test_a_4326_contour_file_interpolates_as_if_reprojected_first(files):
    """Measured before the fix: 0 of 6,000 pixels valid, where the 2927 file gives 6,000."""

    def interpolate(contours, region):
        surface = _surface(files, contours, region_vector=region)
        with rasterio.open(surface.resolve_source(refresh=True)) as handle:
            return handle.read(1, masked=True), handle.crs.to_epsg()

    original, crs = interpolate(files.c2927, files.r2927)
    moved, moved_crs = interpolate(files.c4326, files.r4326)

    assert crs == moved_crs == 2927
    assert original.count() == moved.count() == 6000
    assert (original.mask == moved.mask).all()
    assert float(np.abs(original - moved).max()) < 1e-6
