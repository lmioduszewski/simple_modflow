"""The shared footprint reader, cell rule and bed sampler (the excavation design, C3).

``mf.lak(lakes=)`` and ``LayerStack.excavate(within=)`` read their polygons with ONE
reader (``geometry.read_footprints``), choose cells with ONE rule
(``geometry.footprint_cells``) and sample a bed with ONE sampler
(``surfaces.sample_bed``), so a lake and the pit cut for it cannot disagree about a
file's layer, a cell, or the float at a cell. These tests pin each helper against
the code it was extracted from (``cell_shares`` per feature, ``Surface.values`` with
the model top as ``previous``) on four footprint shapes and on both grid kinds --
``triangulated_grid`` carries ``numpy.int64`` vertex indices, the hand-built one
does not (ledger 265) -- and pin ``mf.lak(layer=)``, which the shared reader gave
the lake.
"""

from __future__ import annotations

import warnings

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
import shapely as shp
from rasterio.transform import from_origin

import myflopy as mf
from myflopy import LAKBuilder, ModelContext
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import (
    cell_shares,
    footprint_cells,
    read_footprints,
    refuse_footprint_misses,
)
from myflopy.surfaces import Surface, sample_bed


@pytest.fixture(scope="module")
def rect_grid():
    """A hand-built 8 x 8 Voronoi grid (python-int vertex indices)."""

    return irregular_voronoi_grid(nrow=8, ncol=8, cell_size=50.0)


@pytest.fixture(params=["rect", "triangulated"])
def grid(request, rect_grid, triangulated_grid):
    """Both grid kinds: the hand-built one and a real ``GridSpec.voronoi`` mesh."""

    return rect_grid if request.param == "rect" else triangulated_grid


def _shapes(grid):
    """Four footprint shapes over ``grid``, in its own frame."""

    x0, y0, x1, y1 = grid.gdf_vorPolys.total_bounds
    mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    span = min(x1 - x0, y1 - y0)
    one = shp.Point(mx, my).buffer(0.22 * span)
    west = shp.Point(x0 + 0.2 * span, y0 + 0.2 * span).buffer(0.11 * span)
    east = shp.Point(x1 - 0.2 * span, y1 - 0.2 * span).buffer(0.11 * span)
    detached = shp.MultiPolygon(
        [shp.Point(mx, my).buffer(0.17 * span),
         shp.box(x1 - 0.15 * span, y0 + 0.03 * span, x1 - 0.04 * span, y0 + 0.14 * span)]
    )
    return {"one": [one], "two_disjoint": [west, east], "detached_part": [detached]}


def test_footprint_cells_is_cell_shares_per_feature(grid):
    """The rule equals the code it was extracted from, feature by feature, on one
    lake, two disjoint lakes and a lake with a detached part."""

    polys = grid.gdf_vorPolys
    for name, geometries in _shapes(grid).items():
        cells = footprint_cells(polys, geometries)
        assert cells == [[cell for cell, _ in cell_shares(polys, g)] for g in geometries], name
        assert all(cells), name
        assert all(type(cell) is int for found in cells for cell in found), name
        assert all(found == sorted(found) for found in cells), name


def test_a_4326_round_trip_covers_the_same_cells(grid, tmp_path):
    """The fourth shape: a footprint written in EPSG:4326 and read back by path
    reaches exactly the cells it reaches in the grid's own frame -- the slivers the
    round trip leaves are dropped by the shared rule (ledger 241)."""

    polys = grid.gdf_vorPolys
    one = _shapes(grid)["one"]
    path = tmp_path / "pond.gpkg"
    gpd.GeoDataFrame({"name": ["pond"]}, geometry=one, crs=grid.crs).to_crs(4326).to_file(path)

    frame = read_footprints(path, crs=grid.crs, id_field="name", noun="lake")
    assert frame.crs == grid.crs
    assert str(frame.attrs["source_crs"].to_epsg()) == "4326"
    assert footprint_cells(polys, frame.geometry) == footprint_cells(polys, one)


def test_lake_cells_are_the_footprint_rule(grid):
    """``LAKBuilder.lake_cells`` IS the footprint rule -- what C4's excavation record
    will be compared against, cell for cell."""

    shapes = _shapes(grid)
    lakes = gpd.GeoDataFrame(
        {"name": ["a", "b"]}, geometry=shapes["two_disjoint"], crs=grid.crs
    )
    builder = LAKBuilder(
        context=ModelContext(grid=grid), nper=1, lakes=lakes, lake_id_field="name",
        starting_stage=1.0, lake_bottom=0.0,
    )
    assert builder.lake_cells == dict(
        zip(["a", "b"], footprint_cells(grid.gdf_vorPolys, shapes["two_disjoint"]), strict=True)
    )


# -- the reader -------------------------------------------------------------------


@pytest.fixture(scope="module")
def two_layers(tmp_path_factory, rect_grid):
    """A GeoPackage with a west layer ``ponds`` and an east layer ``basins``."""

    x0, y0, x1, y1 = rect_grid.gdf_vorPolys.total_bounds
    path = tmp_path_factory.mktemp("footprints") / "site.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"]}, geometry=[shp.box(x0 + 20, y0 + 20, x0 + 120, y0 + 120)],
        crs=rect_grid.crs,
    ).to_file(path, layer="ponds")
    gpd.GeoDataFrame(
        {"name": ["east"]}, geometry=[shp.box(x1 - 120, y1 - 120, x1 - 20, y1 - 20)],
        crs=rect_grid.crs,
    ).to_file(path, layer="basins")
    return path


def test_the_reader_reads_the_layer_it_is_asked_for(two_layers):
    assert read_footprints(two_layers, id_field="name").index.tolist() == ["west"]
    assert read_footprints(two_layers, layer="basins", id_field="name").index.tolist() == ["east"]
    assert read_footprints(two_layers, layer=1, id_field="name").index.tolist() == ["east"]
    with pytest.raises(ValueError, match=r"has no layer 'lakes'\. Layers: 'ponds', 'basins'"):
        read_footprints(two_layers, layer="lakes")
    with pytest.raises(ValueError, match="index 2 is out of range"):
        read_footprints(two_layers, layer=2)


def test_mf_lak_reads_the_named_layer(rect_grid, two_layers):
    """``mf.lak(layer=)`` reaches the reader: the second layer's lake lands on the
    EAST cells. Before it existed only the first layer could be read."""

    ctx = mf.ModelContext(grid=rect_grid)
    rect_grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: np.full(rect_grid.ncpl, 10.0), 1: np.full(rect_grid.ncpl, 0.0)},
        geometry=rect_grid.gdf_vorPolys.geometry, crs=rect_grid.crs,
    )
    kwargs = dict(context=ctx, nper=1, lakes=two_layers, lake_id_field="name",
                  starting_stage=9.0, lake_bottom=5.0, connection_modes="rectangular",
                  lake_top=10.0)
    first = mf.lak(**kwargs)
    east = mf.lak(**kwargs, layer="basins")
    assert first.metadata["lake_ids"] == ["west"]
    assert east.metadata["lake_ids"] == ["east"]
    xc = np.asarray(rect_grid.centroids[0])
    east_cells = {row[2][1] for row in east.options["connectiondata"]}
    assert min(xc[list(east_cells)]) > np.mean(xc)


def test_the_reader_takes_any_integer_index_and_refuses_a_bool(two_layers):
    """``numpy.int64`` subclasses nothing, so an ``isinstance(layer, int)`` test sent it
    down the NAME branch: "has no layer np.int64(1)" (ledger 265's class)."""

    assert read_footprints(two_layers, layer=np.int64(1), id_field="name").index.tolist() == ["east"]
    with pytest.raises(ValueError, match="index 2 is out of range"):
        read_footprints(two_layers, layer=np.int32(2))
    with pytest.raises(TypeError, match="not a bool"):
        read_footprints(two_layers, layer=True)


def test_mf_lak_names_its_vector_layer_and_refuses_an_integer(rect_grid, two_layers):
    """Beside ``only_layer`` an integer reads as a MODEL layer, and read as the file's
    first layer ``layer=0`` was the default -- accepted, changing nothing. Refused,
    pointing at ``only_layer``; a numpy string is stored as a plain ``str``."""

    layers = mf.LayerStack(top=mf.Flat(10.0)).add("a", thickness=10.0).build(rect_grid)
    ctx = mf.ModelContext(grid=rect_grid, domain=layers.idomain, surfaces=layers)
    kwargs = dict(context=ctx, nper=1, lakes=two_layers, lake_id_field="name",
                  starting_stage=9.0, lake_bottom=5.0, connection_modes="rectangular",
                  lake_top=10.0)
    assert mf.lak(**kwargs).metadata["lake_ids"] == ["west"]   # the premise: it builds
    for bad in (0, 1, np.int64(0), True):
        with pytest.raises(TypeError, match=r"VECTOR layer .* a model layer is only_layer="):
            mf.lak(**kwargs, layer=bad)
        with pytest.raises(TypeError, match="only_layer="):
            LAKBuilder(**kwargs, layer=bad)
    builder = LAKBuilder(**kwargs, layer=np.str_("basins"))
    assert type(builder.layer) is str and list(builder.lake_ids) == ["east"]
    with pytest.raises(TypeError, match="only_layer="):
        builder.with_updates(layer=0)


def test_a_layer_with_a_frame_raises_rather_than_being_ignored(rect_grid):
    frame = gpd.GeoDataFrame(geometry=[shp.box(0, 0, 10, 10)], crs=rect_grid.crs)
    with pytest.raises(ValueError, match="has none -- drop layer="):
        read_footprints(frame, layer="ponds")


def test_labels_follow_the_lake_rules(tmp_path, rect_grid):
    """Stems for one feature per file, the id field when given, else ``<noun>_<i>``."""

    a, b = tmp_path / "west.gpkg", tmp_path / "east.gpkg"
    gpd.GeoDataFrame(geometry=[shp.box(0, 0, 10, 10)], crs=rect_grid.crs).to_file(a)
    gpd.GeoDataFrame(geometry=[shp.box(20, 0, 30, 10)], crs=rect_grid.crs).to_file(b)
    assert read_footprints([a, b]).index.tolist() == ["west", "east"]
    two = gpd.GeoDataFrame(
        {"n": [7, 8]}, geometry=[shp.box(0, 0, 1, 1), shp.box(2, 0, 3, 1)], crs=rect_grid.crs
    )
    assert read_footprints(two, noun="pond").index.tolist() == ["pond_0", "pond_1"]
    assert read_footprints(two, id_field="n").index.tolist() == ["7", "8"]
    assert all(type(label) is str for label in read_footprints(two, id_field="n").index)
    with pytest.raises(ValueError, match="Pond ID field is missing: name"):
        read_footprints(two, id_field="name", noun="pond")
    with pytest.raises(ValueError, match="At least one lake geometry is required"):
        read_footprints([], noun="lake")


def test_the_reader_refuses_what_a_footprint_cannot_be(rect_grid):
    line = gpd.GeoDataFrame(geometry=[shp.LineString([(0, 0), (5, 5)])], crs=rect_grid.crs)
    with pytest.raises(ValueError, match="Excavation geometries must be Polygon or MultiPolygon"):
        read_footprints(line, noun="excavation")
    flat = gpd.GeoDataFrame(
        {"name": ["sliver"]}, geometry=[shp.Polygon([(0, 0), (5, 0), (10, 0)])], crs=rect_grid.crs
    )
    with pytest.raises(ValueError, match="Lake geometry has no area: sliver"):
        read_footprints(flat, id_field="name", noun="lake")
    twice = gpd.GeoDataFrame(
        {"name": ["a", "a"]}, geometry=[shp.box(0, 0, 1, 1)] * 2, crs=rect_grid.crs
    )
    with pytest.raises(ValueError, match="Lake IDs must be unique: a"):
        read_footprints(twice, id_field="name", noun="lake")


def test_the_reader_copies_a_frame_and_never_edits_it(rect_grid):
    frame = gpd.GeoDataFrame({"name": ["p"]}, geometry=[shp.box(0, 0, 9, 9)], crs=rect_grid.crs)
    before = frame.copy()
    read_footprints(frame, crs=4326, id_field="name")
    pd.testing.assert_frame_equal(frame, before)


# -- the miss report ---------------------------------------------------------------


def _pair(rect_grid, *, both_off=False):
    x0, y0, _, _ = rect_grid.gdf_vorPolys.total_bounds
    off = shp.box(x0 - 900, y0 - 900, x0 - 800, y0 - 800)
    inside = rect_grid.gdf_vorPolys.geometry.iloc[10]
    frame = gpd.GeoDataFrame(
        {"name": ["in", "out"]}, geometry=[off if both_off else inside, off], crs=rect_grid.crs
    )
    footprints = read_footprints(frame, id_field="name")
    return footprints, footprint_cells(rect_grid.gdf_vorPolys, footprints.geometry)


def test_raise_refuses_any_miss_with_both_extents(rect_grid):
    footprints, cells = _pair(rect_grid)
    with pytest.raises(ValueError) as error:
        refuse_footprint_misses(rect_grid.gdf_vorPolys, footprints, cells, noun="lakes",
                                consequence=", so they would carry no LAK connection")
    message = str(error.value)
    assert message.startswith("1 of 2 lakes cover no grid cell, so they would carry no LAK")
    assert "'out' (outside the grid)" in message and "The lakes span (" in message


def test_warn_reports_a_partial_miss_and_refuses_a_total_one(rect_grid):
    footprints, cells = _pair(rect_grid)
    with pytest.warns(UserWarning, match=r"1 of 2 ponds cover no grid cell and are left out: 'out'"):
        left_out = refuse_footprint_misses(
            rect_grid.gdf_vorPolys, footprints, cells, noun="ponds", partial="warn"
        )
    assert left_out == ["out"]
    footprints, cells = _pair(rect_grid, both_off=True)
    with pytest.raises(ValueError, match="2 of 2 ponds cover no grid cell"):
        refuse_footprint_misses(rect_grid.gdf_vorPolys, footprints, cells, noun="ponds",
                                partial="warn")
    with pytest.raises(ValueError, match="partial must be 'raise' or 'warn'"):
        refuse_footprint_misses(rect_grid.gdf_vorPolys, footprints, cells, partial="ignore")


def test_nothing_missing_is_silent(rect_grid):
    footprints, cells = _pair(rect_grid)
    footprints, cells = footprints.iloc[:1], cells[:1]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert refuse_footprint_misses(rect_grid.gdf_vorPolys, footprints, cells) == []


# -- the bed sampler --------------------------------------------------------------


@pytest.fixture(scope="module")
def bowl(tmp_path_factory, rect_grid):
    """A bowl-shaped bed raster over ``rect_grid``, in its CRS."""

    x0, y0, x1, y1 = rect_grid.gdf_vorPolys.total_bounds
    res = 5.0
    nx, ny = int((x1 - x0) / res) + 2, int((y1 - y0) / res) + 2
    xs = x0 + res * (np.arange(nx) + 0.5)
    ys = y1 + res - res * (np.arange(ny) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    z = 90.0 + ((xx - mx) ** 2 + (yy - my) ** 2) / 900.0
    path = tmp_path_factory.mktemp("bed") / "bowl.tif"
    with rasterio.open(path, "w", driver="GTiff", height=ny, width=nx, count=1,
                       dtype="float64", crs=rect_grid.crs,
                       transform=from_origin(x0, y1 + res, res, res), nodata=-9999.0) as dst:
        dst.write(z, 1)
    return path


def test_sample_bed_is_an_area_sample_against_the_model_top(rect_grid, bowl):
    ground = np.full(rect_grid.ncpl, 100.0)
    area = np.asarray(Surface.raster(bowl).values(rect_grid, previous=ground, method="area"))
    for source in (bowl, str(bowl), Surface.raster(bowl)):
        np.testing.assert_array_equal(sample_bed(source, rect_grid, ground), area)
    centroid = np.asarray(Surface.raster(bowl).values(rect_grid, method="centroid"))
    assert not np.array_equal(area, centroid)  # the bed ignores a build's method


def test_sample_bed_converts_units_and_refuses_them_on_a_relative_surface(rect_grid):
    ground = np.full(rect_grid.ncpl, 100.0)
    feet = sample_bed(Surface.flat(10.0, units="meters"), rect_grid, ground, length_units="feet")
    np.testing.assert_allclose(feet, 32.808398950131235)
    np.testing.assert_array_equal(sample_bed(Surface.offset_below(4.0), rect_grid, ground), 96.0)
    relative = Surface.offset_below(4.0)
    object.__setattr__(relative, "units", "meters")
    with pytest.raises(ValueError, match="only supported on absolute surfaces"):
        sample_bed(relative, rect_grid, ground, length_units="feet")


def test_the_lake_samples_with_the_shared_sampler(rect_grid, bowl, monkeypatch):
    """``LAKBuilder`` samples through ``sample_bed`` -- patching it moves the lake."""

    from myflopy.modflow.mf6 import lakes as lakes_module

    rect_grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: np.full(rect_grid.ncpl, 120.0), 1: np.full(rect_grid.ncpl, 0.0)},
        geometry=rect_grid.gdf_vorPolys.geometry, crs=rect_grid.crs,
    )
    seen = {}

    def spy(surface, grid, ground, *, length_units=None):
        seen["ground"] = np.asarray(ground)
        return sample_bed(surface, grid, ground, length_units=length_units)

    monkeypatch.setattr(lakes_module, "sample_bed", spy)
    lakes = gpd.GeoDataFrame({"name": ["p"]}, geometry=[_shapes(rect_grid)["one"][0]],
                             crs=rect_grid.crs)
    builder = LAKBuilder(context=ModelContext(grid=rect_grid), nper=1, lakes=lakes,
                         lake_id_field="name", starting_stage=110.0, lake_bottom=bowl)
    cell = builder.lake_cells["p"][0]
    assert builder._bottom("p", cell) == sample_bed(bowl, rect_grid, np.full(rect_grid.ncpl, 120.0))[cell]
    np.testing.assert_array_equal(seen["ground"], 120.0)  # the model top is `previous`
