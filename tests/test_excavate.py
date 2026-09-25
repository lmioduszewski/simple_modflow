"""``LayerStack.excavate`` / ``LayerBuildResult.excavate``: a void cut into the stack (ledger 320).

The excavation design (sections 2.1-2.2, 3.1-3.8), C4. What these tests pin:

- after the cut, the first active cell under every cut column has its top EXACTLY at the
  bed, everything above it is idomain 0 (never -1) and zero thick, the ground outside is
  untouched, and a column's thickness below the bed is conserved;
- the sliver rule ABSORBS (the bed never moves): a thin or inactive cell directly under
  the bed is folded into the cell below; a thin cell with nothing active below is kept;
- the footprint's cells are ``mf.lak``'s cells and the bed is ``mf.lak``'s float, through
  the shared reader, cell rule and sampler (C3), on a hand-built grid AND on
  ``GridSpec.voronoi`` grids, whose ``cell2d`` mixes ``numpy.int64`` in (ledger 265);
- the refusals at declaration and at build, the overlap refusal, order independence;
- ``to_disv == build``, ``LayerBuildResult.excavate``, ``attach_to_grid`` publishing the
  records beside ``gdf_topbtm``, the ``prune_isolated`` re-check, QC/report counting;
- and MF6 reading it (slow): the budget closes -- the PERCENT DISCREPANCY is asserted,
  never Normal termination alone -- and MF6's own grid file holds the cut.
"""

from __future__ import annotations

import contextlib
import io
import logging
import pathlib
import types
import warnings

import geopandas as gpd
import numpy as np
import pytest
import rasterio
import shapely as shp
from rasterio.transform import from_origin

import myflopy as mf
from myflopy import LAKBuilder, ModelContext
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import top_active_layer
from myflopy.surfaces import sample_bed

SIDE = 600.0


@pytest.fixture(scope="module")
def rect_grid():
    """A hand-built 10 x 10 Voronoi grid over the same 600 ft square (python-int indices)."""

    return irregular_voronoi_grid(nrow=10, ncol=10, cell_size=60.0)


@pytest.fixture(params=["rect", "triangulated", "triangulated_finer"])
def grid(request, rect_grid, triangulated_grid, triangulated_grid_finer):
    """Three meshes: one hand-built, two real ``GridSpec.voronoi`` ones (28 and 92 cells)."""

    return {
        "rect": rect_grid,
        "triangulated": triangulated_grid,
        "triangulated_finer": triangulated_grid_finer,
    }[request.param]


def _raster(path, z, *, res=5.0, nodata=-9999.0):
    """Write ``z(x, y)`` over the 600 ft square as a GeoTIFF in EPSG:2927."""

    n = int(SIDE / res)
    xs = res * (np.arange(n) + 0.5)
    ys = SIDE - res * (np.arange(n) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    values = np.asarray(z(xx, yy), dtype=float)
    values = np.where(np.isfinite(values), values, nodata)
    with rasterio.open(path, "w", driver="GTiff", height=n, width=n, count=1, dtype="float64",
                       crs="EPSG:2927", transform=from_origin(0.0, SIDE, res, res),
                       nodata=nodata) as dst:
        dst.write(values, 1)
    return path


@pytest.fixture(scope="module")
def gis(tmp_path_factory):
    """A two-layer footprint file, a bowl bed raster, and a raster covering half of it."""

    root = tmp_path_factory.mktemp("excavate_gis")
    site = root / "site.gpkg"
    gpd.GeoDataFrame(
        {"name": ["pond"], "floor": [93.0]}, geometry=[shp.Point(300, 300).buffer(150)], crs=2927
    ).to_file(site, layer="ponds")
    gpd.GeoDataFrame(
        {"name": ["basin"], "floor": [96.0]}, geometry=[shp.box(470, 40, 590, 160)], crs=2927
    ).to_file(site, layer="basins")
    # 92 ft at the centre, 94.05 (a 0.05 ft silt sliver) at r = 101 ft, 96.5 at the
    # pond's 150 ft edge: every pond cell is cut, through fill, silt and into sand.
    bowl = _raster(root / "bowl.tif", lambda x, y: 92.0 + ((x - 300) ** 2 + (y - 300) ** 2) / 5000.0)
    # ...and a steeper one reaching the ground (100) at r = 126 ft, so the pond has a
    # RIM whose bed is above the ground.
    steep = _raster(root / "steep.tif", lambda x, y: 92.0 + ((x - 300) ** 2 + (y - 300) ** 2) / 2000.0)
    half = _raster(root / "half.tif",
                   lambda x, y: np.where(x < 300, 90.0 + 0 * y, np.nan))
    return types.SimpleNamespace(root=root, site=site, bowl=bowl, steep=steep, half=half)


def _stack(**pinch):
    """Top 100; fill 100-95, silt 95-94 (a 1 ft unit), sand 94-64 in two slices, clay 64-44."""

    return (mf.LayerStack(top=mf.Flat(100.0))
            .add("fill", thickness=5.0)
            .add("silt", thickness=1.0, **pinch)
            .add("sand", thickness=30.0, split=2)
            .add("clay", thickness=20.0))


def _layer_tops(result, cell):
    return np.concatenate([[result.top[cell]], result.botm[:-1, cell]])


def _connected(result, cell):
    """The first ACTIVE layer of a column (what a lake's vertical connection reaches)."""

    return int(np.flatnonzero(result.idomain[:, cell] > 0)[0])


def _pond(gis, **kw):
    """The bowl pond, named ``pond`` (the default would be the file stem, ``site``)."""

    kw.setdefault("name", "pond")
    return _stack().excavate(bottom=mf.Raster(gis.bowl), within=gis.site, layer="ponds", **kw)


# -- the cut ---------------------------------------------------------------------------


def test_the_connected_top_is_the_bed_and_nothing_active_is_left_inside(grid, gis):
    baseline = _stack().build(grid)
    result = _pond(gis).build(grid)
    record = result.excavations["pond"]
    cut = 0
    for cell, bed in zip(record.cells, record.bed, strict=True):
        if not np.isfinite(bed) or bed >= baseline.top[cell]:
            assert result.top[cell] == baseline.top[cell]  # a rim cell is not cut
            continue
        cut += 1
        k = _connected(result, cell)
        # EXACT: the connected cell's top is the declared bed, not a nudged copy
        assert _layer_tops(result, cell)[k] == bed
        assert result.top[cell] == bed
        assert (result.idomain[:k, cell] == 0).all()          # 0, never -1
        assert (result.thickness[:k, cell] == 0).all()
        assert result.excavated[:k, cell].all() and not result.excavated[k:, cell].any()
        assert result.thickness[k, cell] >= record.min_remnant
        # the thickness under the bed is conserved and the stack's base did not move
        assert result.botm[-1, cell] == baseline.botm[-1, cell]
        assert result.thickness[:, cell].sum() == pytest.approx(bed - baseline.botm[-1, cell])
        # nothing active sits above the bed anywhere in the pond
        tops = _layer_tops(result, cell)
        assert not ((result.idomain[:, cell] > 0) & (tops > bed)).any()
    assert cut >= 3
    outside = np.setdiff1d(np.arange(grid.ncpl), record.cells)
    np.testing.assert_array_equal(result.top[outside], baseline.top[outside])
    np.testing.assert_array_equal(result.botm[:, outside], baseline.botm[:, outside])
    np.testing.assert_array_equal(result.idomain[:, outside], baseline.idomain[:, outside])
    np.testing.assert_array_equal(result.ground, baseline.top)
    assert not result.excavated[:, outside].any()


def test_the_record_holds_python_types(grid, gis):
    """Ledger 265: nothing numpy-typed leaks into the record, on a grid whose cell2d
    carries numpy.int64."""

    record = _pond(gis).build(grid).excavations["pond"]
    assert type(record.name) is str
    assert all(type(cell) is int for cell in record.cells)
    assert list(record.cells) == sorted(record.cells)
    assert record.bed.dtype == np.float64 and not record.bed.flags.writeable
    assert type(record.min_remnant) is float


def test_numpy_scalars_are_accepted_and_normalized(rect_grid, tmp_path):
    area = mf.Area(rect_grid.gdf_vorPolys.geometry.loc[[44, 45]].union_all())
    result = _stack().excavate(bottom=np.float64(93.0), within=area, name=np.str_("pond"),
                               min_remnant=np.float32(0.25)).build(rect_grid)
    record = result.excavations["pond"]
    assert type(record.name) is str and type(record.min_remnant) is float
    assert record.min_remnant == pytest.approx(0.25)
    assert list(result.excavations) == ["pond"] and type(next(iter(result.excavations))) is str


def test_the_cells_are_the_lakes_cells(grid, gis, tmp_path):
    """Same polygons, same cells as ``mf.lak`` -- by path (with its layer), through an
    ``Area``, through a GeoDataFrame, and from an EPSG:4326 copy of the file."""

    lake = LAKBuilder(context=ModelContext(grid=grid), nper=1, lakes=gis.site, layer="ponds",
                      lake_id_field="name", starting_stage=100.0, lake_bottom=90.0)
    lake_cells = tuple(sorted({c for cells in lake.lake_cells.values() for c in cells}))
    frame = gpd.read_file(gis.site, layer="ponds")
    as_4326 = tmp_path / "pond4326.gpkg"
    frame.to_crs(4326).to_file(as_4326)
    routes = {
        "path": dict(within=gis.site, layer="ponds"),
        "area": dict(within=mf.Area.coerce(frame)),
        "frame": dict(within=frame),
        "4326": dict(within=as_4326),
    }
    for label, kw in routes.items():
        record = _stack().excavate(bottom=90.0, name="pond", **kw).build(grid).excavations["pond"]
        assert record.cells == lake_cells, label


def test_two_disjoint_lakes_and_a_detached_part(rect_grid, tmp_path):
    path = tmp_path / "two.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west", "east"]},
        geometry=[shp.Point(150, 300).buffer(80),
                  shp.MultiPolygon([shp.Point(450, 300).buffer(80), shp.box(500, 520, 580, 590)])],
        crs=2927,
    ).to_file(path)
    lake = LAKBuilder(context=ModelContext(grid=rect_grid), nper=1, lakes=path,
                      lake_id_field="name", starting_stage=100.0, lake_bottom=90.0)
    union = tuple(sorted({c for cells in lake.lake_cells.values() for c in cells}))
    record = _stack().excavate(bottom=90.0, within=path).build(rect_grid).excavations["two"]
    assert record.cells == union


def test_the_bed_is_the_lakes_float(grid, gis):
    """The record's bed is ``sample_bed`` at its cells -- area-weighted whatever
    ``build(method=)`` says -- and the lake, declared on the EXCAVATED grid, samples the
    identical float (what C7's exact bed check rests on)."""

    baseline = _stack().build(grid)
    for method in ("area", "centroid"):
        result = _pond(gis).build(grid, method=method, attach=True)
        record = result.excavations["pond"]
        expected = sample_bed(mf.Raster(gis.bowl), grid, baseline.top)[list(record.cells)]
        np.testing.assert_array_equal(record.bed, expected)
    lake = LAKBuilder(context=ModelContext(grid=grid), nper=1, lakes=gis.site, layer="ponds",
                      lake_id_field="name", starting_stage=100.0, lake_bottom=mf.Raster(gis.bowl))
    np.testing.assert_array_equal(lake._sample(lake.lake_bottom)[list(record.cells)], record.bed)


# -- the sliver rule ---------------------------------------------------------------------


def _one_cell(grid, cell=44):
    return mf.Area(grid.gdf_vorPolys.geometry.loc[cell])


def test_a_sliver_is_absorbed_and_the_bed_never_moves(rect_grid):
    """Bed 94.05: 0.05 ft of silt would be left. Absorbed -- its bottom rises to the bed,
    it goes idomain 0 -- and sand_1 connects with its top EXACTLY at 94.05, 15.05 thick.
    A deepen rule would have moved the bed to 94.0."""

    result = _stack().excavate(94.05, _one_cell(rect_grid), name="p").build(rect_grid)
    col = result.idomain[:, 44].tolist()
    assert col == [0, 0, 1, 1, 1]
    assert _layer_tops(result, 44)[2] == 94.05 and result.top[44] == 94.05
    assert result.thickness[2, 44] == pytest.approx(15.05)
    assert result.absorbed[1, 44] == pytest.approx(0.05)
    assert result.excavated[:2, 44].all()
    report = result.excavation_report().loc["p"]
    assert report["absorbed"] == 1 and report["absorbed_units"] == {"silt": 1}
    assert report["connected_layers"] == {"sand_1": 1}


def test_min_remnant_zero_keeps_the_sliver_and_a_larger_one_absorbs_more(rect_grid):
    kept = _stack().excavate(94.05, _one_cell(rect_grid), min_remnant=0).build(rect_grid)
    assert kept.idomain[:, 44].tolist() == [0, 1, 1, 1, 1]
    assert kept.thickness[1, 44] == pytest.approx(0.05)
    half_foot = _stack().excavate(94.5, _one_cell(rect_grid)).build(rect_grid)
    assert half_foot.idomain[:, 44].tolist() == [0, 1, 1, 1, 1]      # 0.5 >= 0.1: kept
    one_foot = _stack().excavate(94.5, _one_cell(rect_grid), min_remnant=1.0).build(rect_grid)
    assert one_foot.idomain[:, 44].tolist() == [0, 0, 1, 1, 1]       # 0.5 < 1.0: absorbed
    assert _layer_tops(one_foot, 44)[2] == 94.5


@pytest.mark.parametrize("policy, value", [("inactive", 0), ("passthrough", -1)])
def test_an_inactive_cell_under_the_bed_is_absorbed_whatever_min_remnant(rect_grid, policy, value):
    """A pinched silt (1 ft < its 2 ft minimum) under a bed at 94.5: a lake's vertical
    connection cannot reach through it, so it is absorbed even at min_remnant=0, and the
    pond reads 0 there, never -1."""

    stack = _stack(pinch=policy, min_thickness=2.0)
    baseline = stack.build(rect_grid, min_sep=0.1)
    assert baseline.idomain[1, 44] == value
    for remnant in (0.1, 0):
        result = (_stack(pinch=policy, min_thickness=2.0)
                  .excavate(94.5, _one_cell(rect_grid), min_remnant=remnant)
                  .build(rect_grid, min_sep=0.1))
        assert result.idomain[:, 44].tolist() == [0, 0, 1, 1, 1]
        assert _layer_tops(result, 44)[2] == 94.5
        assert result.excavation_report().iloc[0]["absorbed"] == 1


def test_a_chain_of_absorptions_counts_each_cells_own_material(rect_grid):
    """A 0.05 ft silt remnant, then a 0.01 ft pass-through cell under it: both absorbed.
    Each records only its OWN material under the bed, so the column sums to what lay
    between the bed and the connected cell's ORIGINAL top -- 0.06 ft, not the 0.11 a
    cumulative ``bed - botm`` gave by counting the silt's remnant a second time."""

    def stack():
        return (mf.LayerStack(top=mf.Flat(100.0)).add("fill", thickness=5.0)
                .add("silt", thickness=1.0)
                .add("mid", thickness=0.5, pinch="passthrough", min_thickness=1.0)
                .add("sand", thickness=30.0))

    base = stack().build(rect_grid, min_sep=0.01)
    assert base.idomain[:, 44].tolist() == [1, 1, -1, 1]      # the premise: mid passes through
    assert base.thickness[2, 44] == pytest.approx(0.01)
    result = stack().excavate(94.05, _one_cell(rect_grid), name="p").build(rect_grid, min_sep=0.01)
    assert result.idomain[:, 44].tolist() == [0, 0, 0, 1]
    assert _layer_tops(result, 44)[3] == 94.05                   # the bed never moves
    np.testing.assert_allclose(result.absorbed[:, 44], [0.0, 0.05, 0.01, 0.0], atol=1e-9)
    assert result.absorbed[:, 44].sum() == pytest.approx(94.05 - base.botm[2, 44])
    area = float(rect_grid.gdf_vorPolys.geometry.area.iloc[44])
    row = result.excavation_report().loc["p"]
    assert row["absorbed"] == 2 and row["absorbed_units"] == {"silt": 1, "mid": 1}
    assert row["absorbed_volume"] == pytest.approx((94.05 - base.botm[2, 44]) * area)
    assert result.qc().absorbed == 2


def test_a_zero_thickness_cell_absorbed_below_a_sliver_still_counts(rect_grid):
    """Its own share is 0.0 (``reconcile=False`` leaves a zero-thickness unit), so the
    count cannot read ``absorbed > 0``: every excavated cell from a column's first
    absorbed one down was absorbed, not cut."""

    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("fill", thickness=5.0)
             .add("gone", bottom=mf.Flat(95.0), pinch="inactive", min_thickness=1.0)
             .add("sand", thickness=30.0))
    result = stack.excavate(95.05, _one_cell(rect_grid), name="p").build(rect_grid, reconcile=False)
    assert result.idomain[:, 44].tolist() == [0, 0, 1]
    np.testing.assert_allclose(result.absorbed[:, 44], [0.05, 0.0, 0.0], atol=1e-9)
    row = result.excavation_report().loc["p"]
    assert row["absorbed"] == 2 and row["absorbed_units"] == {"fill": 1, "gone": 1}
    assert result.qc().absorbed == 2 and result.qc().excavated == [1, 1, 0]


def test_a_thin_cell_with_nothing_active_below_is_kept(rect_grid):
    result = _stack().excavate(44.05, _one_cell(rect_grid), name="deep").build(rect_grid)
    assert result.idomain[:, 44].tolist() == [0, 0, 0, 0, 1]
    assert result.thickness[4, 44] == pytest.approx(0.05)
    assert result.excavation_report().loc["deep", "kept_slivers"] == 1


def test_a_column_with_no_active_cell_under_the_bed_raises(rect_grid):
    with pytest.raises(ValueError, match=r"no active cell under its bed.*'deep': \[44\]"):
        _stack().excavate(43.0, _one_cell(rect_grid), name="deep").build(rect_grid)
    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=36.0)
             .add("clay", thickness=20.0, pinch="inactive", min_thickness=30.0)
             .excavate(63.0, _one_cell(rect_grid), name="deep"))
    with pytest.raises(ValueError, match="no active cell under its bed"):
        stack.build(rect_grid)


def test_rim_cells_are_not_cut_and_are_logged(rect_grid, gis, caplog):
    with caplog.at_level(logging.INFO, logger="myflopy.surfaces"):
        result = (_stack().excavate(mf.Raster(gis.steep), gis.site, layer="ponds", name="pond")
                  .build(rect_grid))
    record = result.excavations["pond"]
    rim = [c for c, b in zip(record.cells, record.bed, strict=True) if b >= 100.0]
    assert rim, "the premise: the bowl rises above the ground inside the pond"
    assert (result.top[rim] == 100.0).all() and (result.idomain[:, rim] == 1).all()
    assert any("at or above the ground" in r.message for r in caplog.records)
    assert result.excavation_report().loc["pond", "rim"] == len(rim)


def test_cells_with_no_bed_are_not_cut_and_warned(rect_grid, gis):
    with pytest.warns(UserWarning, match="the bottom has no data at"):
        result = (_stack().excavate(mf.Raster(gis.half), gis.site, layer="ponds")
                  .build(rect_grid))
    record = result.excavations["site"]
    missing = [c for c, b in zip(record.cells, record.bed, strict=True) if not np.isfinite(b)]
    assert missing and (result.top[missing] == 100.0).all()
    empty = _raster(gis.root / "empty.tif", lambda x, y: np.where(x > 590, 90.0 + 0 * y, np.nan))
    with pytest.raises(ValueError, match="no data at any of its"):
        _stack().excavate(mf.Raster(empty), gis.site, layer="ponds").build(rect_grid)


def test_fill_propagate_leaves_a_no_data_cell_uncut_without_a_warning(rect_grid, gis):
    """``fill="propagate"`` takes the ground where the bed has no data, so those cells
    read as rim (bed == ground) and are not cut -- the same top ``mf.lak`` falls back
    to -- and nothing warns, because the caller asked for it."""

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = (_stack().excavate(mf.Surface.raster(gis.half, fill="propagate"), gis.site,
                                    layer="ponds", name="pond").build(rect_grid))
    record = result.excavations["pond"]
    east = [c for c, b in zip(record.cells, record.bed, strict=True) if b == 100.0]
    west = [c for c, b in zip(record.cells, record.bed, strict=True) if b == 90.0]
    assert east and west and np.isfinite(record.bed).all()
    assert (result.top[east] == 100.0).all() and (result.top[west] == 90.0).all()
    # ...and the lake, sampling the same bed against the EXCAVATED top, reads the
    # record's float at every cell: a bare propagate source is idempotent here.
    bed = mf.Surface.raster(gis.half, fill="propagate")
    lake = LAKBuilder(context=ModelContext(grid=rect_grid, domain=result.idomain, surfaces=result),
                      nper=1, lakes=gis.site, layer="ponds", lake_id_field="name",
                      starting_stage=100.0, lake_bottom=bed)
    np.testing.assert_array_equal(lake._sample(lake.lake_bottom)[list(record.cells)], record.bed)


@pytest.mark.parametrize(
    "build",
    [
        lambda half: mf.Surface.raster(half, fill="propagate") - 2.0,
        lambda half: mf.Min(mf.Surface.raster(half, fill="propagate"), mf.Flat(95.0)),
        lambda half: mf.Surface.clamp(mf.Surface.raster(half, fill="propagate") + 1.0, upper=99.0),
    ],
    ids=["shift", "min", "nested"],
)
def test_algebra_over_a_propagate_bed_is_refused(rect_grid, gis, build):
    """Where the source has no data it reads the surface ABOVE, so algebra over it is
    not idempotent: ``- 2`` would be cut to ``ground - 2`` by the stack and read as
    ``ground - 4`` by the lake off the excavated top (14 of 27 cells 2.0 ft apart)."""

    with pytest.raises(ValueError, match="built on a fill='propagate' source") as error:
        _stack().excavate(build(gis.half), gis.site, layer="ponds")
    assert "pass the propagate source BARE" in str(error.value)


# -- declaration and build refusals ------------------------------------------------------


@pytest.mark.parametrize(
    "bottom, match",
    [
        (mf.Surface.offset_below(5.0), "measured from the surface ABOVE"),
        (mf.Surface.constant_thickness(5.0), "measured from the surface ABOVE"),
        (mf.Isopach(mf.Flat(5.0)), "measured from the surface ABOVE"),
        (mf.Toward(mf.Flat(50.0), 0.5), "measured from the surface ABOVE"),
        (mf.Min(mf.Flat(95.0), mf.Surface.offset_below(2.0)), "measured from the surface ABOVE"),
    ],
)
def test_a_relative_bottom_is_refused_naming_the_absolute_spelling(rect_grid, bottom, match):
    with pytest.raises(ValueError, match=match) as error:
        _stack().excavate(bottom, _one_cell(rect_grid))
    assert "ground - 9" in str(error.value)


def test_an_absolute_shift_is_accepted(rect_grid):
    ground = mf.Flat(100.0)
    result = _stack().excavate(ground - 9, _one_cell(rect_grid)).build(rect_grid)
    assert result.top[44] == 91.0


def test_declaration_refusals(rect_grid, gis, tmp_path):
    cell = _one_cell(rect_grid)
    with pytest.raises(ValueError, match="needs within= to be a FILE"):
        _stack().excavate("floor", cell)
    with pytest.raises(ValueError, match="neither a column of site.gpkg"):
        _stack().excavate("nope", gis.site, layer="ponds")
    with pytest.raises(ValueError, match="layer='ponds' selects a layer of a FILE"):
        _stack().excavate(90.0, cell, layer="ponds")
    with pytest.raises(ValueError, match="has no layer 'lakes'"):
        _stack().excavate(90.0, gis.site, layer="lakes")
    with pytest.raises(FileNotFoundError):
        _stack().excavate(90.0, tmp_path / "missing.gpkg")
    with pytest.raises(FileNotFoundError):
        _stack().excavate(pathlib.Path(tmp_path / "missing.tif"), cell)
    with pytest.raises(ValueError, match="must be polygonal"):
        _stack().excavate(90.0, shp.LineString([(0, 0), (300, 300)]))
    with pytest.raises(ValueError, match="must be polygonal"):
        _stack().excavate(90.0, mf.Area(shp.Point(300, 300)))
    line = tmp_path / "line.gpkg"
    gpd.GeoDataFrame(geometry=[shp.LineString([(0, 0), (300, 300)])], crs=2927).to_file(line)
    with pytest.raises(ValueError, match="Excavation geometries must be Polygon"):
        _stack().excavate(90.0, line)
    for positional in ([44, 45], np.array([44, 45]), np.zeros(rect_grid.ncpl, bool), range(3)):
        with pytest.raises(TypeError, match="not cell numbers or a boolean mask"):
            _stack().excavate(90.0, positional)
    with pytest.raises(TypeError, match="got a dict"):
        _stack().excavate(90.0, {"x": 1})
    with pytest.raises(TypeError, match="not a bool"):
        _stack().excavate(True, cell)
    with pytest.raises(TypeError, match="got a list"):
        _stack().excavate([90.0], cell)
    for bad in (-0.1, float("nan")):
        with pytest.raises(ValueError, match="finite number >= 0"):
            _stack().excavate(90.0, cell, min_remnant=bad)
    with pytest.raises(TypeError, match="min_remnant must be a number"):
        _stack().excavate(90.0, cell, min_remnant="0.1")
    stack = _stack().excavate(90.0, cell, name="pond")
    with pytest.raises(ValueError, match="already exists"):
        stack.excavate(90.0, _one_cell(rect_grid, 12), name="pond")
    for taken in ("sand", "sand_1", "top"):
        with pytest.raises(ValueError, match="names a unit or layer"):
            _stack().excavate(90.0, cell, name=taken)
    with pytest.raises(ValueError, match="already names an excavation"):
        stack.add("pond", thickness=5.0)
    # ...and no layer a unit BUILDS later may take it: a generated split name, an
    # explicit names= entry, through add, insert_below and replace alike.
    with pytest.raises(ValueError, match=r"\['deep_1'\] already names an excavation"):
        _stack().excavate(90.0, cell, name="deep_1").add("deep", thickness=5.0, split=2)
    with pytest.raises(ValueError, match=r"\['deep_2'\] already names an excavation"):
        _stack().excavate(90.0, cell, name="deep_2").insert_below("silt", "deep", thickness=1.0,
                                                                   split=2)
    with pytest.raises(ValueError, match=r"\['pond'\] already names an excavation"):
        stack.replace("clay", split=2, names=["pond", "deep"])
    with pytest.raises(ValueError, match=r"\['clay_2'\] already names an excavation"):
        _stack().excavate(90.0, cell, name="clay_2").replace("clay", split=2)
    assert _stack().excavate(90.0, cell, name="clay_3").replace("clay", split=2).names


def test_a_null_column_floor_is_refused_at_the_call(rect_grid, tmp_path):
    path = tmp_path / "nulls.gpkg"
    gpd.GeoDataFrame({"name": ["a", "b"], "floor": [93.0, None]},
                     geometry=[shp.box(100, 100, 200, 200), shp.box(300, 300, 400, 400)],
                     crs=2927).to_file(path)
    with pytest.raises(ValueError, match="a NULL floor is not a floor"):
        _stack().excavate("floor", path)


def test_a_footprint_on_no_cell_raises_with_both_extents(rect_grid, tmp_path):
    away = tmp_path / "away.gpkg"
    gpd.GeoDataFrame(geometry=[shp.box(5000, 5000, 5100, 5100)], crs=2927).to_file(away)
    with pytest.raises(ValueError, match=r"1 of 1 features of excavation 'away' cover no grid cell"
                       r".*span \(5,000\.0"):
        _stack().excavate(90.0, away).build(rect_grid)


def test_a_partial_miss_warns_by_name_and_digs_the_rest(rect_grid, tmp_path):
    two = tmp_path / "two.gpkg"
    gpd.GeoDataFrame(geometry=[shp.Point(300, 300).buffer(60), shp.box(5000, 5000, 5100, 5100)],
                     crs=2927).to_file(two)
    # named by ROW: the reader's own fallback, 'excavation_1', read as another
    # excavation's default name
    with pytest.warns(UserWarning, match=r"1 of 2 features of excavation 'two' cover no grid cell"
                      r".*: 'row 1' \(outside the grid\)"):
        result = _stack().excavate(90.0, two).build(rect_grid)
    assert len(result.excavations["two"].cells) >= 1


def test_overlapping_excavations_are_refused_naming_both(rect_grid):
    a = mf.Area(rect_grid.gdf_vorPolys.geometry.loc[[44, 45]].union_all())
    b = mf.Area(rect_grid.gdf_vorPolys.geometry.loc[[45, 46]].union_all())
    stack = _stack().excavate(93.0, a, name="a").excavate(92.0, b, name="b")
    with pytest.raises(ValueError, match=r"'a' and 'b' share 1 cell\(s\) \[45\]"):
        stack.build(rect_grid)
    ok = _stack().excavate(93.0, a, name="a").excavate(
        92.0, mf.Area(rect_grid.gdf_vorPolys.geometry.loc[[47]].union_all()), name="b")
    assert set(ok.build(rect_grid).excavations) == {"a", "b"}


def test_several_excavations_do_not_depend_on_declaration_order(rect_grid, gis):
    east = mf.Area(shp.box(470, 40, 590, 160), name="basin")
    ab = _pond(gis).excavate(94.05, east).build(rect_grid)
    ba = (_stack().excavate(94.05, east)
          .excavate(mf.Raster(gis.bowl), gis.site, layer="ponds", name="pond").build(rect_grid))
    for field in ("top", "botm", "idomain", "thickness", "excavated", "absorbed"):
        np.testing.assert_array_equal(getattr(ab, field), getattr(ba, field))
    assert ab.excavations == ba.excavations and set(ab.excavations) == {"pond", "basin"}


def test_a_column_floor_per_feature(rect_grid, gis, tmp_path):
    result = _stack().excavate("floor", gis.site, layer="ponds").build(rect_grid)
    record = result.excavations["site"]
    assert set(record.bed.tolist()) == {93.0}
    clash = tmp_path / "clash.gpkg"
    gpd.GeoDataFrame({"floor": [93.0, 95.0]},
                     geometry=[shp.box(200, 200, 330, 330), shp.box(310, 310, 400, 400)],
                     crs=2927).to_file(clash)
    with pytest.raises(ValueError, match="give one cell different floor floors"):
        _stack().excavate("floor", clash).build(rect_grid)


def test_a_bed_in_other_units_is_converted(rect_grid):
    result = (_stack().excavate(mf.Surface.flat(28.5, units="meters"), _one_cell(rect_grid))
              .build(rect_grid))
    assert result.top[44] == pytest.approx(28.5 * 3.280839895013123)


# -- the result ------------------------------------------------------------------------


def test_to_disv_is_build(grid, gis):
    result = _pond(gis).build(grid)
    spec = _pond(gis).to_disv(grid)
    np.testing.assert_array_equal(spec.options["top"], result.top)
    np.testing.assert_array_equal(spec.options["botm"], result.botm)
    np.testing.assert_array_equal(spec.options["idomain"], result.idomain)
    eager = mf.disv(result)
    np.testing.assert_array_equal(eager.options["top"], result.top)
    assert list(grid.excavations) == ["pond"]  # to_disv attaches by default
    bound = mf.disv(_pond(gis).for_grid(grid))  # a stack that carries its grid: eager
    np.testing.assert_array_equal(bound.options["top"], result.top)
    np.testing.assert_array_equal(bound.options["idomain"], result.idomain)


def test_a_built_result_excavates_through_the_same_engine(rect_grid, gis):
    stack_cut = _pond(gis).build(rect_grid)
    plain = _stack().build(rect_grid)
    before = {f: np.copy(getattr(plain, f)) for f in ("top", "botm", "idomain")}
    result_cut = plain.excavate(mf.Raster(gis.bowl), gis.site, layer="ponds", name="pond")
    for field in ("top", "botm", "idomain", "thickness", "excavated", "absorbed", "ground"):
        np.testing.assert_array_equal(getattr(result_cut, field), getattr(stack_cut, field))
    assert result_cut.excavations == stack_cut.excavations
    for field, value in before.items():
        np.testing.assert_array_equal(getattr(plain, field), value)  # self unchanged
    assert plain.excavations == {} and plain.ground is None
    with pytest.raises(ValueError, match="share"):
        result_cut.excavate(90.0, _one_cell(rect_grid, 44), name="again")
    corner = result_cut.excavate(95.5, _one_cell(rect_grid, 0), name="corner")
    assert set(corner.excavations) == {"pond", "corner"}
    np.testing.assert_array_equal(corner.ground, plain.top)
    assert corner.top[0] == 95.5


def test_attach_publishes_the_records_and_any_other_table_clears_them(rect_grid, gis):
    result = _pond(gis).build(rect_grid, attach=True)
    assert rect_grid.excavations == result.excavations
    record = result.excavations["pond"]
    cut = [c for c, b in zip(record.cells, record.bed, strict=True) if b < 100.0]
    np.testing.assert_array_equal(rect_grid.gdf_topbtm[0].to_numpy()[cut], record.bed[
        [record.cells.index(c) for c in cut]])
    rect_grid.gdf_topbtm = rect_grid.gdf_topbtm.copy()  # any other writer of the table...
    assert rect_grid.excavations == {}                   # ...drops the records
    _stack().build(rect_grid, attach=True)               # another stack's attach
    assert rect_grid.excavations == {}


def test_prune_isolated_rechecks_the_excavated_columns(rect_grid):
    """A connected cell whose neighbours are all inactive is isolated; pruning it would
    leave the pond's column with nothing under the bed -- refused, naming the cell."""

    base = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", thickness=10.0)
            .add("lower", thickness=10.0)).build(rect_grid)
    neighbours = [int(c) for c in rect_grid.find_adjacent_cells(44)]
    idomain = base.idomain.copy()
    idomain[1, neighbours] = 0
    from dataclasses import replace
    walled = replace(base, idomain=idomain)
    dug = walled.excavate(90.0, _one_cell(rect_grid), name="pond")
    assert dug.idomain[:, 44].tolist() == [0, 1]
    with pytest.raises(ValueError, match=r"prune_isolated\(\) left an excavated column.*'pond': \[44\]"):
        dug.prune_isolated()


def test_qc_report_and_n_pinched_count_the_cut_on_its_own(rect_grid):
    plain = _stack().build(rect_grid)
    dug = _stack().excavate(94.05, _one_cell(rect_grid)).build(rect_grid)
    assert dug.n_pinched == plain.n_pinched == 0
    qc = dug.qc()
    assert qc.ok and qc.pinched == [0] * 5 and qc.excavated == [1, 1, 0, 0, 0] and qc.absorbed == 1
    assert qc.thin == plain.qc().thin
    assert "excavated: 2 cell(s) removed by excavation, 1 of them absorbed" in str(qc)
    assert "2 excavated in 1 excavation(s)" in dug.report()
    assert plain.qc().excavated is None and "excavated" not in str(plain.qc())
    # LayerStack.qc measures reconcile against the PRE-cut surfaces
    stack_qc = _stack().excavate(94.05, _one_cell(rect_grid)).qc(rect_grid)
    assert stack_qc.reconcile_adjusted == _stack().qc(rect_grid).reconcile_adjusted
    assert stack_qc.excavated == [1, 1, 0, 0, 0]


def test_the_report_counts_what_was_dug(rect_grid, gis):
    result = _pond(gis).build(rect_grid)
    report = result.excavation_report()
    row = report.loc["pond"]
    record = result.excavations["pond"]
    assert row["cells"] == len(record.cells) == row["cut"] + row["rim"] + row["no_data"]
    assert row["removed"] == int(result.excavated.sum())
    area = rect_grid.gdf_vorPolys.geometry.area.to_numpy()
    cells = list(record.cells)
    assert row["removed_volume"] == pytest.approx(((result.ground - result.top)[cells] * area[cells]).sum())
    assert sum(row["connected_layers"].values()) == row["cut"]
    assert row["connected_thin"] == 0 or row["connected_thin"] <= row["cut"]
    assert result.excavation_report(thinner_than=100.0).loc["pond", "connected_thin"] == row["cut"]
    assert list(_stack().build(rect_grid).excavation_report().columns) == list(report.columns)


def test_connected_thin_judges_a_split_unit_whole_as_qc_does(rect_grid):
    """sand is 30 ft, min_thickness 20, split in three 10 ft slices. A bed on sand_1's
    base connects a FULL sand_2 slice -- thinner than its unit's minimum by construction,
    and NOT thin: 20 ft of sand is left. A bed 5 ft into sand_2 leaves 15 ft: thin, in
    the report and in ``qc()`` alike."""

    stack = lambda: (mf.LayerStack(top=mf.Flat(100.0)).add("fill", thickness=5.0)  # noqa: E731
                     .add("sand", thickness=30.0, min_thickness=20.0, split=3)
                     .add("clay", thickness=20.0))
    two = mf.Area(rect_grid.gdf_vorPolys.geometry.loc[[44, 45]].union_all())
    whole = stack().excavate(85.0, two, name="p").build(rect_grid)
    row = whole.excavation_report().loc["p"]
    assert row["connected_layers"] == {"sand_2": 2} and row["connected_min"] == pytest.approx(10.0)
    assert row["connected_thin"] == 0 == sum(whole.qc().thin)
    assert whole.excavation_report(thinner_than=11.0).loc["p", "connected_thin"] == 2
    short = stack().excavate(80.0, two, name="p").build(rect_grid)
    row = short.excavation_report().loc["p"]
    assert row["connected_layers"] == {"sand_2": 2}
    assert row["connected_thin"] == 2 == short.qc().thin[2]


def test_cache_status_prefixes_an_excavation_bed(gis, tmp_path):
    contours = tmp_path / "bed_contours.gpkg"
    gpd.GeoDataFrame({"z": [92.0]}, geometry=[shp.LineString([(250, 300), (350, 300)])],
                     crs=2927).to_file(contours)
    bed = mf.Contours(contours, z="z", region_vector=gis.site, out=tmp_path / "bed.tif")
    stack = _stack().excavate(bed, gis.site, layer="ponds", name="pond")
    assert stack.cache_status() == {"excavation:pond": "missing"}


def test_a_stack_with_an_excavation_is_not_written_down_yet(rect_grid, gis):
    """Interim, until the lifecycle commit (ledger 320): refused loudly at the call, never
    dropped -- and ``for_grid`` keeps the pond."""

    stack = _pond(gis)
    with pytest.raises(NotImplementedError, match="cannot be written down yet"):
        stack.to_dict()
    with pytest.raises(NotImplementedError, match=r"mf\.disv\(stack\.build\(vor\)\)"):
        mf.disv(stack)
    bound = stack.for_grid(rect_grid)
    np.testing.assert_array_equal(bound.build().top, stack.build(rect_grid).top)
    assert _stack().to_dict()  # an un-excavated stack still serializes


def test_top_active_recharge_lands_under_the_bed(rect_grid, gis):
    """Section 3.8: with the excavated domain, ``irch="top_active"`` names the connected
    cell of every cut column (MF6 zeroes it while a lake is wet there -- C7 warns)."""

    result = _pond(gis).build(rect_grid)
    irch = top_active_layer(result.idomain)
    record = result.excavations["pond"]
    for cell, bed in zip(record.cells, record.bed, strict=True):
        if np.isfinite(bed) and bed < 100.0:
            assert irch[cell] == _connected(result, cell)
            assert _layer_tops(result, cell)[irch[cell]] == bed


def test_the_pictures_draw_an_excavated_stack(rect_grid, gis, tmp_path):
    """No default reconciled build produced ``top == botm`` before; every picture of one."""

    result = _stack().excavate(94.05, gis.site, layer="ponds").build(rect_grid)
    assert (result.thickness[:2][result.excavated[:2]] == 0).all()
    result.plot.map().save(tmp_path / "map.png")
    result.plot.map("silt").save(tmp_path / "silt.png")
    fig = result.plot.map("silt", basemap=True, logscale=True).fig
    assert fig.data
    result.plot.section(y=300.0).save(tmp_path / "section.png")
    result.plot.section(y=300.0, color_by="thickness").save(tmp_path / "thk.png")
    assert result.vertex_grid().top[44] == pytest.approx(result.top[44])
    pytest.importorskip("pyvista")
    scene = result.plot.grid(backend="vtk")
    assert scene.meshes


# -- MF6 ------------------------------------------------------------------------------


def _mf6_run(root, grid, result, tag):
    """A steady NEWTON model on ``result``: CHD at 98 ft on the east edge (the water table
    stands above the pond's bed), array recharge on the top active cell, IMS closed tight."""

    east = [int(c) for c in np.flatnonzero(grid.gdf_vorPolys.geometry.bounds.maxx >= SIDE - 1e-6)]
    context = ModelContext(grid=grid, domain=result.idomain)
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result),
        mf.npf(k=10.0, k33=1.0, icelltype=1),
        mf.ic(strt=99.0),
        mf.chd(stress_period_data={0: [[(4, c), 98.0] for c in east]}),
        mf.rch.array(context=context, recharge=1e-3),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-6, inner_dvclose=1e-8,
               rcloserecord=[1e-4, "STRICT"]),
    ])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root / f"p_{tag}", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    listing = (pathlib.Path(run.workspace) / "m.lst").read_text()
    discrepancy = float(
        [line for line in listing.splitlines() if "PERCENT DISCREPANCY" in line][-1].split()[-1]
    )
    return run, discrepancy


@pytest.mark.slow
def test_mf6_reads_the_cut_and_closes_its_budget(tmp_path, triangulated_grid_finer, gis, caplog):
    """Run MF6 on an excavated DISV (NEWTON, a bowl bed cut through fill and a 1 ft silt,
    an absorbed sliver). Asserted: the PERCENT DISCREPANCY -- the design measured a
    steady toy at Normal termination with a 53.89% error, so Normal termination alone
    proves nothing -- MF6's own grid file holding the cut, a head in every connected
    cell and none in a dug one. The same model un-excavated is run beside it: with IMS's
    DEFAULT closure they read -0.21% and -0.41% (a closure artefact on 360 ft3/d under
    Normal termination; ``docs/reproducers/excavation_build_checks.py``), which is why
    the solver is closed tight here, where both read 0.00. A PRT release point in layer
    0 over a dug cell lands on no active cell and is reported, not moved down (ledger
    242)."""

    import flopy

    from myflopy.modflow.mf6.prt import PRTReleasePoints

    grid = triangulated_grid_finer
    result = _pond(gis).build(grid, attach=True)
    record = result.excavations["pond"]
    cut = [c for c, b in zip(record.cells, record.bed, strict=True) if b < 100.0]
    assert len(cut) >= 5 and (result.absorbed > 0).any(), "the premise: a cut with a sliver"

    _, baseline = _mf6_run(tmp_path, grid, _stack().build(grid), "baseline")
    run, discrepancy = _mf6_run(tmp_path, grid, result, "excavated")
    assert abs(baseline) < 0.1 and abs(discrepancy) < 0.1

    workspace = pathlib.Path(run.workspace)
    grb = flopy.mf6.utils.MfGrdFile(str(workspace / "m.disv.grb"))
    np.testing.assert_allclose(np.asarray(grb.top), result.top, rtol=0, atol=1e-6)
    np.testing.assert_allclose(np.asarray(grb.bot).reshape(result.botm.shape), result.botm,
                               rtol=0, atol=1e-6)
    np.testing.assert_array_equal(np.asarray(grb.idomain).reshape(result.idomain.shape),
                                  result.idomain)
    head = np.asarray(flopy.utils.HeadFile(str(workspace / "m.hds")).get_data()).reshape(
        result.idomain.shape)
    for cell in cut:
        k = _connected(result, cell)
        assert abs(head[k, cell]) < 1e29                     # solved
        assert (np.abs(head[:k, cell]) >= 1e29).all()        # the dug cells are out

    # the transmissivity map of the model built on the cut draws (design section 3.7)
    live = run.model("m")
    assert live.plot.map(type="trans").fig.data

    # PRT reports a partial miss through its LOGGER only (not warnings.warn -- ledger 320
    # records that ledger 288's rule has not reached it), so it is read from there.
    dug_top = [c for c in cut if result.idomain[0, c] == 0]   # layer 0 removed there
    assert dug_top
    inside = grid.gdf_vorPolys.geometry.loc[dug_top[0]].representative_point()
    with caplog.at_level(logging.WARNING, logger="myflopy.modflow.mf6.prt"):
        points = PRTReleasePoints.from_points(live, [(inside.x, inside.y), (30.0, 30.0)], layer=0)
    assert len(points.packagedata) == 1 and points.packagedata[0][1][0] == 0
    assert any("1 of 2 release points land on no active cell in layer 0" in r.getMessage()
               for r in caplog.records)
