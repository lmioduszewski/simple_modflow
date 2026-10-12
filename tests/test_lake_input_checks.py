"""What ``mf.lak``, ``LakeTableBuilder`` and ``LakeTable`` accept, and refuse, since ledger 318.

Ledger 318 listed what the lake family accepted without checking -- a footprint path
as a ``str``, a footprint never reprojected, a NaN pixel, tables MF6 refuses, integer
keys, arrays written whole into every period, a storage coefficient of 0, counts that
disagreed with their blocks -- and documented each one as behaviour. Each is fixed; this
file pins the fixed behaviour, and the MF6 runs pin the claims about MF6 that the fixes
rest on: ``LakeTable``'s checks ARE MF6's, the DEM mode's one kept zero-volume row gives
MF6 the requested curve, and an array rainfall now runs where it killed MF6.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString, box

import myflopy as mf
from myflopy import (
    LAKBuilder,
    LakeConnection,
    LakeOutlet,
    LakeTable,
    LakeTableBuilder,
    ModelContext,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

# --------------------------------------------------------------------------- fixtures


def _toy_grid() -> VoronoiGridPlus:
    """4 x 4 unit cells, two layers (tops 12, 8, 0), written clockwise so MF6 can run it."""

    n = 4
    verts = np.array([[x, y] for y in range(n + 1) for x in range(n + 1)], float)
    iverts, centers = [], []
    for r in range(n):
        for c in range(n):
            v0 = r * (n + 1) + c
            iverts.append([v0, v0 + n + 1, v0 + n + 2, v0 + 1])  # clockwise, as MF6 wants
            centers.append([c + 0.5, r + 0.5])
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centers))
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [12.0] * 16, 1: [8.0] * 16, 2: [0.0] * 16}, geometry=grid.gdf_vorPolys.geometry, crs=grid.crs
    )
    grid.get_disu_connectivity(validate=True)
    return grid


@pytest.fixture(scope="module")
def toy():
    grid = _toy_grid()
    lakes = gpd.GeoDataFrame({"name": ["pond"], "num": [101]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)
    return {"grid": grid, "context": ModelContext(grid=grid, domain=np.ones((2, 16), dtype=int)),
            "lakes": lakes}


def _lak(toy, *, key="name", nper=3, **kw):
    """``mf.lak`` on the toy pond: a rectangular basin, bottom 9, top 11."""

    base = dict(context=toy["context"], nper=nper, lakes=toy["lakes"], lake_id_field=key,
                starting_stage=11.5, lake_bottom=9.0, connection_modes="rectangular",
                lake_top=11.0, bed_leakance=0.1, unexcavated="ignore")
    return mf.lak(**{**base, **kw})


@pytest.fixture(scope="module")
def bed(tmp_path_factory) -> Path:
    """A 2 x 2 raster of 1 ft pixels, EPSG:2927: bed 8, 9 / 9, 10."""

    path = tmp_path_factory.mktemp("lakebed") / "bed.tif"
    _raster(path, [[8.0, 9.0], [9.0, 10.0]])
    return path


def _raster(path, values, crs="EPSG:2927", nodata=None):
    values = np.asarray(values, dtype="float64")
    with rasterio.open(path, "w", driver="GTiff", height=values.shape[0], width=values.shape[1],
                       count=1, dtype="float64", crs=crs, transform=from_origin(0, 2, 1, 1),
                       nodata=nodata) as dst:
        dst.write(values, 1)
    return path


STAGES = [9, 10, 11]
WHOLE = ((9.0, 1.0, 1.0), (10.0, 4.0, 3.0), (11.0, 8.0, 4.0))
LEFT = ((9.0, 1.0, 1.0), (10.0, 3.0, 2.0), (11.0, 5.0, 2.0))
RIGHT = ((9.0, 0.0, 0.0), (10.0, 1.0, 1.0), (11.0, 3.0, 2.0))


def _dem(bed, **kw):
    return LakeTableBuilder(dem=bed, **{"stages": STAGES, **kw}).build().rows


# ------------------------------------------------------------- (a) footprint types


def test_every_footprint_form_reads_the_same_outline(bed, tmp_path):
    path = tmp_path / "fp.gpkg"
    frame = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(0, 0, 2, 2)], crs=2927)
    frame.to_file(path)
    forms = {
        "Path": path,
        "str": str(path),
        "list of paths": [path],
        "GeoDataFrame": frame,
        "GeoSeries": frame.geometry,
        "Area": mf.Area(box(0, 0, 2, 2), crs="EPSG:2927"),
        "geometry": box(0, 0, 2, 2),
    }
    for label, footprint in forms.items():
        assert _dem(bed, footprint=footprint) == WHOLE, label


def test_footprint_layer_reads_the_layer_named(bed, tmp_path):
    two = tmp_path / "two.gpkg"
    gpd.GeoDataFrame(geometry=[box(0, 0, 1, 2)], crs=2927).to_file(two, layer="a")
    gpd.GeoDataFrame(geometry=[box(1, 0, 2, 2)], crs=2927).to_file(two, layer="b")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # pyogrio's "More than one layer found"
        assert _dem(bed, footprint=two) == LEFT, "None reads the first layer, GDAL's default"
    assert _dem(bed, footprint=two, footprint_layer="b") == RIGHT
    assert _dem(bed, footprint=str(two), footprint_layer=np.int64(1)) == RIGHT, "an index"
    with pytest.raises(ValueError, match=r"has no layer 'c'\. Layers: 'a', 'b'"):
        _dem(bed, footprint=two, footprint_layer="c")
    with pytest.raises(ValueError, match="selects a layer of a footprint FILE"):
        _dem(bed, footprint=box(0, 0, 2, 2), footprint_layer="b")
    with pytest.raises(TypeError, match="not a bool"):
        _dem(bed, footprint=two, footprint_layer=True)


def test_a_footprint_that_is_not_an_area_is_refused(bed, tmp_path):
    with pytest.raises(ValueError, match="must be Polygon or MultiPolygon"):
        _dem(bed, footprint=gpd.GeoSeries([LineString([(0, 0), (2, 2)])], crs=2927))
    with pytest.raises(TypeError, match="footprint got a dict"):
        _dem(bed, footprint={"pond": box(0, 0, 2, 2)})
    with pytest.raises(Exception, match="No such file") as missing:
        _dem(bed, footprint=tmp_path / "nope.gpkg")
    assert type(missing.value).__name__ == "DataSourceError"
    with pytest.raises(ValueError, match="shrinks the footprint to nothing"):
        _dem(bed, footprint=box(0, 0, 2, 2), footprint_buffer=-5.0)


# ------------------------------------------------------ (b) reprojection and misses


def test_a_footprint_in_another_crs_is_reprojected_to_the_rasters(bed, tmp_path):
    path = tmp_path / "fp4326.gpkg"
    gpd.GeoDataFrame(geometry=[box(0, 0, 2, 2)], crs=2927).to_crs(4326).to_file(path)
    assert _dem(bed, footprint=path) == WHOLE, "it built a table of zeros before ledger 318"


@pytest.mark.parametrize(
    "footprint, reason",
    [
        (box(-122.5, 47.0, -122.4, 47.1), "it lies off the raster"),
        (box(100, 100, 101, 101), "it lies off the raster"),
        (box(0, 0, 0.4, 2), r"no pixel CENTRE falls inside it"),
    ],
    ids=["degrees-over-state-plane", "off-the-raster", "narrower-than-a-pixel"],
)
def test_a_footprint_reaching_no_pixel_is_refused_with_both_boxes(bed, footprint, reason):
    with pytest.raises(ValueError, match=reason) as refused:
        _dem(bed, footprint=footprint)
    message = str(refused.value)
    assert "which MF6 refuses" in message
    assert "The footprint spans (" in message and "and the raster (0.0, 0.0, 2.0, 2.0)" in message
    assert "carries no CRS, so it was read as the raster's (EPSG:2927)" in message


def test_the_refusal_names_a_reprojection(bed, tmp_path):
    far = tmp_path / "far.gpkg"
    gpd.GeoDataFrame(geometry=[box(-122.5, 47.0, -122.4, 47.1)], crs=4326).to_file(far)
    with pytest.raises(ValueError, match="reprojected from EPSG:4326 to the raster's EPSG:2927"):
        _dem(bed, footprint=far)


def test_a_footprint_partly_off_the_raster_warns_with_the_area(bed):
    with pytest.warns(UserWarning, match=r"4 of the lake-table footprint's 8 .* lies off bed\.tif"):
        assert _dem(bed, footprint=box(0, 0, 4, 2)) == WHOLE
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # less than one pixel off: within the table's resolution
        assert _dem(bed, footprint=box(0, 0, 2.5, 2)) == WHOLE


def test_a_raster_in_degrees_is_refused(tmp_path):
    degrees = _raster(tmp_path / "deg.tif", [[8.0, 9.0], [9.0, 10.0]], crs="EPSG:4326")
    with pytest.raises(ValueError, match=r"geographic CRS \(EPSG:4326\): its pixels measure degrees"):
        _dem(degrees, footprint=box(0, 0, 2, 2))


# ------------------------------------------------------------- (c) no-data pixels


def test_a_nan_pixel_is_no_bed_and_is_reported(tmp_path):
    nan = _raster(tmp_path / "nan.tif", [[8.0, np.nan], [9.0, 10.0]])
    with pytest.warns(UserWarning, match=r"1 of the 4 pixels of nan\.tif .* have no bed"):
        rows = _dem(nan, footprint=box(0, 0, 2, 2))
    assert rows == ((9.0, 1.0, 1.0), (10.0, 3.0, 2.0), (11.0, 6.0, 3.0)), "NaN volumes before"
    declared = _raster(tmp_path / "nodata.tif", [[8.0, -9999.0], [9.0, 10.0]], nodata=-9999.0)
    with pytest.warns(UserWarning, match="have no bed"):
        assert _dem(declared, footprint=box(0, 0, 2, 2)) == rows, "NaN is read as declared no-data is"
    empty = _raster(tmp_path / "empty.tif", np.full((2, 2), np.nan))
    with pytest.raises(ValueError, match=r"all 4 pixels whose centres fall inside it have no bed"):
        _dem(empty, footprint=box(0, 0, 2, 2))


def test_each_lake_table_is_built_once_per_mf_lak_call(tmp_path, toy, monkeypatch):
    """``prepared_tables`` rebuilt every ``LakeTableBuilder`` on each read, and ``validate()``
    read it inside its per-lake loop: two DEM-table lakes built 12 tables, and each table's
    no-data warning came 6 times, growing with the square of the lake count."""

    nan = _raster(tmp_path / "nan.tif", [[8.0, np.nan], [9.0, 10.0]])
    built = []
    real = LakeTableBuilder.build
    monkeypatch.setattr(LakeTableBuilder, "build", lambda self: built.append(self) or real(self))
    lakes = gpd.GeoDataFrame({"name": ["pond", "west"]}, geometry=[box(1, 1, 3, 3), box(0, 3, 1, 4)],
                             crs=toy["grid"].crs)
    tables = {name: LakeTableBuilder(dem=nan, footprint=box(0, 0, 2, 2), stages=STAGES)
              for name in ("pond", "west")}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        spec = _lak({**toy, "lakes": lakes}, starting_stage=10.5, tables=tables, unexcavated="warn")
    assert len(built) == 2, "one build per table"
    assert sum("have no bed" in str(w.message) for w in caught) == 2, "one warning per table"
    assert len(spec.options["lake_tables"]) == 2

    builder = LAKBuilder(context=toy["context"], nper=1, lakes=toy["lakes"], lake_id_field="name",
                         starting_stage=10.5, lake_bottom=9.0,
                         tables={"pond": LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)})
    first = builder.prepared_tables
    table = first["pond"]
    first.clear()
    assert builder.prepared_tables["pond"] is table, "built once, and a copy handed out, not the cache"
    builder.tables["pond"] = LakeTableBuilder(area=8.0, bottom=9.0, top=12.0)
    assert builder.prepared_tables["pond"].rows != table.rows, "a replaced table is rebuilt"


# ------------------------------------------------- (d) LakeTable is MF6's table check


def test_lake_table_refuses_what_mf6_refuses_and_passes_what_it_passes():
    with pytest.raises(ValueError, match=r"volumes must STRICTLY increase.*keep only the last"):
        LakeTable([(8.0, 0.0, 0.0), (9.0, 0.0, 0.0), (10.0, 5.0, 4.0)])
    with pytest.raises(ValueError, match="volumes must STRICTLY increase"):
        LakeTable([(8.0, 0.0, 0.0), (9.0, 5.0, 4.0), (10.0, 5.0, 4.0)])
    with pytest.raises(ValueError, match="surface areas must not fall"):
        LakeTable([(8.0, 0.0, 5.0), (9.0, 1.0, 4.0)])
    with pytest.raises(ValueError, match="stages must be strictly increasing"):
        LakeTable([(8.0, 0.0, 0.0), (8.0, 1.0, 4.0)])
    assert LakeTable([(8.0, 0.0, 4.0), (9.0, 1.0, 4.0)]).rows[1] == (9.0, 1.0, 4.0), "equal area passes"
    assert LakeTable([(8.0, 3.0, 0.0)]).rows == ((8.0, 3.0, 0.0),), "one row, nothing to compare"
    # Beyond MF6's neighbouring-row checks: finite values and no negative volume or
    # area (MF6 runs a negative one, and stops on a NaN only because FloPy cuts the row).
    for rows, match in (
        ([(1.0, 0.0, 0.0), (2.0, float("nan"), 1.0)], "volume of nan"),
        ([(1.0, 0.0, 0.0), (float("inf"), 1.0, 1.0)], "stage of inf"),
        ([(1.0, -5.0, 0.0), (2.0, 1.0, 1.0)], "negative volume"),
        ([(1.0, 0.0, -1.0), (2.0, 1.0, 0.0)], "negative surface area"),
    ):
        with pytest.raises(ValueError, match=match):
            LakeTable(rows)


def test_stages_at_or_below_the_bed_keep_only_the_last_zero_row(bed, caplog):
    with caplog.at_level(logging.INFO, logger="myflopy.modflow.mf6.lakes"):
        rows = _dem(bed, footprint=box(0, 0, 2, 2), stages=[6, 7, 7.5, 8, 9, 10])
    assert rows == ((8.0, 0.0, 0.0), (9.0, 1.0, 1.0), (10.0, 4.0, 3.0))
    assert "4 stages at or below the lowest bed" in caplog.text and "kept only the last, 8" in caplog.text
    assert _dem(bed, footprint=box(0, 0, 2, 2), stages=[8, 9]) == ((8.0, 0.0, 0.0), (9.0, 1.0, 1.0))
    with pytest.raises(ValueError, match=r"Every stage \(6 to 8\) is at or below the lowest bed .*, 8:"):
        _dem(bed, footprint=box(0, 0, 2, 2), stages=[6, 7, 8])


@pytest.mark.parametrize(
    "stages, match",
    [([], "stages is empty"), ([9, float("nan")], "must be finite"), ([9, 11, 10], "strictly increasing")],
    ids=["empty", "nan", "unsorted"],
)
def test_dem_stages_are_checked_before_the_raster_is_read(tmp_path, stages, match):
    with pytest.raises(ValueError, match=match):
        LakeTableBuilder(dem=tmp_path / "missing.tif", footprint=box(0, 0, 2, 2), stages=stages).build()


# ------------------------------------------------------------ (f) rows=DataFrame


def test_a_dataframe_is_read_by_its_column_names_or_by_position():
    want = ((96.0, 0.0, 0.0), (98.0, 1.0, 2.0))
    named = pd.DataFrame({"Area": [0.0, 2.0], "note": ["a", "b"], "STAGE": [96.0, 98.0],
                          "volume": [0.0, 1.0]}, index=[10, 20])
    assert LakeTableBuilder(rows=named).build().rows == want, "by name, any case, extras ignored"
    assert LakeTable(named).rows == want
    assert LakeTable(pd.DataFrame([[96.0, 0.0, 0.0], [98.0, 1.0, 2.0]])).rows == want, "by position"
    assert LakeTable(pd.DataFrame({"s": [96.0, 98.0], "v": [0.0, 1.0], "a": [0.0, 2.0]})).rows == want
    with pytest.raises(ValueError, match=r"has a stage and a volume column but no area"):
        LakeTable(pd.DataFrame({"stage": [96.0], "volume": [0.0], "sarea": [0.0]}))
    with pytest.raises(ValueError, match=r"has 4 columns \(0, 1, 2, 3\), none named"):
        LakeTable(pd.DataFrame([[96.0, 0.0, 0.0, 1.0]]))
    with pytest.raises(TypeError, match="not a dict"):
        LakeTable({"stage": [96.0]})


# ----------------------------------------------------------- (g) storage_coefficient


@pytest.mark.parametrize("coefficient", [0, 0.0, -1.0, 1.5, float("nan")])
def test_a_storage_coefficient_outside_zero_to_one_is_refused(bed, coefficient):
    with pytest.raises(ValueError, match=r"storage_coefficient must be greater than 0 and at most 1"):
        LakeTableBuilder(area=1.0, bottom=9.0, top=11.0, storage_coefficient=coefficient).build()
    with pytest.raises(ValueError, match=r"storage_coefficient must be greater than 0 and at most 1"):
        _dem(bed, footprint=box(0, 0, 2, 2), storage_coefficient=coefficient)


def test_a_storage_coefficient_in_range_builds_and_rows_mode_ignores_it():
    assert LakeTableBuilder(area=1.0, bottom=9.0, top=11.0, storage_coefficient=1).build().rows[-1] == (11.0, 2.0, 1.0)
    assert LakeTableBuilder(area=1.0, bottom=9.0, top=11.0, storage_coefficient=0.25).build().rows[-1][1] == 0.5
    with pytest.raises(TypeError, match="must be a number"):
        LakeTableBuilder(area=1.0, bottom=9.0, top=11.0, storage_coefficient=True).build()
    rows = [(9.0, 0.0, 1.0), (10.0, 1.0, 1.0)]
    assert LakeTableBuilder(rows=rows, storage_coefficient=0.0).build().rows == tuple(rows)


# ------------------------------------------------------------- (e) keys by lake id


def test_an_integer_key_names_the_lake_of_its_string(toy):
    table = LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)
    for key in ("101", 101, np.int64(101)):
        assert _lak(toy, key="num", tables={key: table}).options["tables"] == [[0, "lak_101.lak.tab"]], key
    with pytest.raises(ValueError, match="tables names lake '101' twice"):
        _lak(toy, key="num", tables={101: table, "101": table})
    with pytest.raises(ValueError, match=r"unknown lake 'pond2' \(lakes: 'pond'\)"):
        _lak(toy, tables={"pond2": table})
    with pytest.raises(TypeError, match=r"tables\['pond'\] is a tuple"):
        _lak(toy, tables={"pond": ((9.0, 0.0, 4.0),)})
    rows = _lak(toy, key="num", only_vertical={101: True}).options["connectiondata"]
    assert {row[3] for row in rows} == {"VERTICAL"}, "was ignored silently"
    with pytest.raises(ValueError, match="names lake '101' twice"):
        _lak(toy, key="num", starting_stage={101: 11.5, "101": 11.5})
    # Every other per-lake mapping reads its key the same way; each of these was reverted
    # alone with every lake test still passing, two of them silent (review of ledger 318).
    rows = _lak(toy, key="num", lake_bottom=6.0, only_layer={101: 1}).options["connectiondata"]
    assert {row[2][0] for row in rows} == {1}, "only_layer={101: 1} was read as no restriction"

    def horizontal(spec):
        return sum(row[3] == "HORIZONTAL" for row in spec.options["connectiondata"])

    assert horizontal(_lak(toy, key="num")) == 8, "edges only"
    assert horizontal(_lak(toy, key="num", rectangular_interior={101: True})) == 16, "was read as False"
    rows = _lak(toy, key="num", connection_modes={101: [LakeConnection((0, 5))]}).options["connectiondata"]
    assert rows == [[0, 0, (0, 5), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]]
    leakance = {(101, "vertical"): 0.3, (np.int64(101), "horizontal"): 0.7}
    rows = _lak(toy, key="num", bed_leakance=leakance).options["connectiondata"]
    assert {(row[3], row[4]) for row in rows} == {("VERTICAL", 0.3), ("HORIZONTAL", 0.7)}
    assert _amounts(toy, key="num", rainfall={0: {101: 0.2}}) == [[0.2], [], []]


def test_an_integer_indexed_series_is_read_by_lake_not_by_position(toy):
    """The Series branch tested ``lake in index`` with the STRING id, so a Series keyed by the
    integer ids was read by position: out of lake order, each lake took another's value."""

    lakes = gpd.GeoDataFrame({"num": [101, 102]}, geometry=[box(0, 0, 2, 2), box(2, 2, 4, 4)],
                             crs=toy["grid"].crs)
    two = {**toy, "lakes": lakes}

    def stages(**kw):
        return [row[1] for row in _lak(two, key="num", **kw).options["packagedata"]]

    def per_lake(field, **kw):
        rows = _lak(two, key="num", **kw).options["connectiondata"]
        return sorted({(row[0], row[field]) for row in rows})

    assert stages(starting_stage={102: 12.0, 101: 11.5}) == [11.5, 12.0], "the dict"
    assert stages(starting_stage=pd.Series({102: 12.0, 101: 11.5})) == [11.5, 12.0], "was [12.0, 11.5]"
    index = pd.Index(np.array([102, 101], dtype=np.int64))
    assert stages(starting_stage=pd.Series([12.0, 11.5], index=index)) == [11.5, 12.0]
    assert stages(starting_stage=pd.Series([11.5, 12.0])) == [11.5, 12.0], "no lake label: by position"
    assert per_lake(4, bed_leakance=pd.Series({102: 0.5, 101: 0.1})) == [(0, 0.1), (1, 0.5)]
    assert per_lake(5, lake_bottom=pd.Series({102: 9.5, 101: 9.0})) == [(0, 9.0), (1, 9.5)], (
        "one bed per lake; it was read by cell and refused"
    )
    with pytest.raises(ValueError, match=r"starting_stage names lake '101' twice \(101, '101'\)"):
        stages(starting_stage=pd.Series([11.5, 11.5, 12.0], index=pd.Index([101, "101", 102], dtype=object)))


def test_a_whole_grid_bed_over_integer_lake_ids_is_still_read_by_cell(toy):
    """An integer is a CELL number to ``lake_bottom``/``lake_top`` as well as a lake id. A
    Series holding every lake cell stays the per-cell map it always was, even when lakes
    numbered 1 and 2 make its index name every lake too -- read by label, each lake would
    get ONE bed, the value at the cell numbered like it."""

    lakes = gpd.GeoDataFrame({"num": [1, 2]}, geometry=[box(0, 0, 2, 2), box(2, 2, 4, 4)],
                             crs=toy["grid"].crs)
    bed = pd.Series({cell: 9.0 + 0.1 * cell for cell in range(16)})
    rows = _lak({**toy, "lakes": lakes}, key="num", lake_bottom=bed,
                connection_modes="bathy").options["connectiondata"]
    floors = {lake: sorted(row[5] for row in rows if row[3] == "VERTICAL" and row[0] == lake)
              for lake in (0, 1)}
    assert floors == {0: [bed[c] for c in (0, 1, 4, 5)], 1: [bed[c] for c in (10, 11, 14, 15)]}


# ------------------------------------------------------------- (h) per-period values


def _amounts(toy, key="name", **kw):
    data = _lak(toy, key=key, **kw).options["perioddata"]
    return [[row[2] for row in data.get(period, [])] for period in range(3)]


def test_an_array_or_series_is_read_by_period(toy):
    assert _amounts(toy, rainfall=np.array([0.1, 0.2, 0.3])) == [[0.1], [0.2], [0.3]]
    assert _amounts(toy, rainfall=pd.Series([0.1, 0.2])) == [[0.1], [0.2], []]
    assert _amounts(toy, rainfall=pd.Series({"pond": 0.4})) == [[0.4], [0.4], [0.4]], "a lake-id index"
    assert _amounts(toy, rainfall={"pond": np.array([0.1, None, 0.3], dtype=object)}) == [[0.1], [], [0.3]]
    amounts = _amounts(toy, rainfall=np.float64(0.5))
    assert amounts == [[0.5]] * 3 and all(type(a[0]) is float for a in amounts), "no numpy scalar written"
    rate = _lak(toy, outlets=[LakeOutlet("pond", rate=np.array([-1.0, -2.0]))]).options["perioddata"]
    assert rate == {0: [[0, "RATE", -1.0]], 1: [[0, "RATE", -2.0]]}  # period 2 sets nothing
    assert LakeOutlet("pond", rate=np.array([-1.0, -2.0])).rate == [-1.0, -2.0]
    assert LakeOutlet("pond", rate=pd.Series({1: -2.0})).rate == {1: -2.0}
    with pytest.raises(ValueError, match=r"2-D array of shape \(2, 3\)"):
        _lak(toy, rainfall=np.zeros((2, 3)))
    with pytest.raises(ValueError, match="2-D array"):
        LakeOutlet("pond", rate=np.zeros((2, 2)))


def test_a_value_that_would_never_be_written_is_refused(toy):
    with pytest.raises(ValueError, match=r"rainfall: 101 is not a stress period .* key a lake by its string id \('101'\)"):
        _lak(toy, key="num", rainfall={101: 0.5})
    assert _amounts(toy, key="num", rainfall={"101": 0.5}) == [[0.5]] * 3
    with pytest.raises(ValueError, match="rainfall has 4 values for 3 stress period"):
        _lak(toy, rainfall=[0.1, 0.2, 0.3, 0.4])
    with pytest.raises(ValueError, match=r"rainfall\['pond'\]: 5 is not a stress period"):
        _lak(toy, rainfall={"pond": {5: 0.1}})
    with pytest.raises(ValueError, match=r"outlets\[0\]\.rate: 3 is not a stress period"):
        _lak(toy, outlets=[LakeOutlet("pond", rate={3: -1.0})])
    with pytest.raises(ValueError, match=r"outlets\[0\]\.rate has 4 values"):
        _lak(toy, outlets=[LakeOutlet("pond", rate=[-1.0] * 4)])
    with pytest.raises(ValueError, match=r"rainfall names lake 'lake9', which no lake has"):
        _lak(toy, rainfall={"lake9": 0.1})


def test_integer_keys_that_are_both_periods_and_lake_ids_are_refused(toy):
    lakes = gpd.GeoDataFrame({"num": [1]}, geometry=[box(1, 1, 3, 3)], crs=toy["grid"].crs)
    small = {**toy, "lakes": lakes}
    with pytest.raises(ValueError, match=r"every key \(1\) is both a stress period and a lake id"):
        _lak(small, key="num", rainfall={1: 0.5})
    assert _amounts(small, key="num", rainfall={0: 0.5, 1: 0.2}) == [[0.5], [0.2], []], "0 is no lake"
    assert _amounts(small, key="num", rainfall={1: {"1": 0.2}}) == [[], [0.2], []]
    assert _amounts(small, key="num", rainfall={"1": 0.3}) == [[0.3]] * 3


def test_values_mf6_refuses_are_refused_here(toy):
    for name in ("rainfall", "evaporation", "runoff", "inflow", "withdrawals"):
        with pytest.raises(ValueError, match=f"{name} is -1.0: MF6 refuses a negative"):
            _lak(toy, **{name: -1.0})
    with pytest.raises(ValueError, match="status 'BOGUS' is not a lake status"):
        _lak(toy, status="BOGUS")
    with pytest.raises(ValueError, match="stage value nan is not finite"):
        _lak(toy, stage=float("nan"))
    with pytest.raises(ValueError, match=r"rainfall value \[0\.1\] is not a number"):
        _lak(toy, rainfall={0: [0.1]})
    data = _lak(toy, status=["constant", "Active"], stage=11.0, rainfall="rain_ts").options["perioddata"]
    assert [row[2] for row in data[0]] == [11.0, "rain_ts", "constant"], "any case; a string is a time series"


# ---------------------------------------------------------- (h) options, walls, flopy


def test_nlakes_and_pname_are_refused_in_options(toy):
    with pytest.raises(TypeError, match=r"nlakes=5 is not for \*\*options"):
        _lak(toy, nlakes=5)
    with pytest.raises(TypeError, match=r"pname='zz' is not for \*\*options: rename the package with name="):
        _lak(toy, pname="zz")
    spec = _lak(toy, filename="pond.lak", tables={"pond": LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)})
    assert spec.options["filename"] == "pond.lak", "names the input file only"
    assert (spec.name, spec.metadata["lak_index"]["package"]) == ("lak", "lak")
    assert spec.options["lake_tables"][0]["pname"] == "lak_pond_table"
    renamed = _lak(toy, name="pond", tables={"pond": LakeTableBuilder(area=4.0, bottom=9.0, top=12.0)})
    assert (renamed.options["pname"], renamed.metadata["lak_index"]["package"]) == ("pond", "pond")


def test_a_rectangular_lake_refuses_a_bed_no_layer_holds_as_bathy_does(toy):
    base = dict(context=toy["context"], nper=1, lakes=toy["lakes"], lake_id_field="name",
                starting_stage=11.5, bed_leakance=0.1, unexcavated="ignore")
    for mode in ({"connection_modes": "rectangular", "lake_top": 11.0}, {}):
        with pytest.raises(ValueError, match=r"bottom 6.0 does not intersect active cell 5: .*only_layer=0"):
            LAKBuilder(**base, **mode, lake_bottom=6.0, only_layer=0).connectiondata


def test_a_horizontal_lake_connection_needs_a_length_and_a_width():
    with pytest.raises(ValueError, match=r"needs connection_length > 0 \(got 0.0\)"):
        LakeConnection((0, 1), "HORIZONTAL")
    with pytest.raises(ValueError, match=r"needs connection_width > 0 \(got 0.0\)"):
        LakeConnection((0, 1), "horizontal", connection_length=1.0)
    with pytest.raises(ValueError, match="connection_length > 0"):
        LakeConnection((0, 1), "HORIZONTAL", connection_length=float("nan"), connection_width=1.0)
    assert LakeConnection((0, 1)).connection_length == 0.0, "VERTICAL reads neither"
    assert LakeConnection((0, 1), "HORIZONTAL", connection_length=0.5, connection_width=1.0).connection_width == 1.0


def test_mf_lak_flopy_refuses_counts_that_disagree_with_their_blocks():
    records = dict(packagedata=[[0, 101.0, 1]], perioddata={0: []},
                   connectiondata=[[0, 0, (0, 42), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]])
    outlet = [[0, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]]
    assert mf.lak.flopy(**records, nlakes=1).options["nlakes"] == 1
    assert mf.lak.flopy(**records, outlets=outlet).options["noutlets"] == 0, "FloPy raises it to the rows"
    assert mf.lak.flopy(**records, noutlets=1, outlets=outlet).options["noutlets"] == 1
    assert mf.lak.flopy(**{**records, "packagedata": [[0, 101.0, 1, "pond"]]}, boundnames=True)
    assert mf.lak.flopy(**{**records, "packagedata": [[0, 101.0, 1, 5.0]]}, auxiliary="temp")
    cases = {
        "nlakes=3 but packagedata has 1 row": dict(nlakes=3),
        "noutlets=2 but outlets has 1 row": dict(noutlets=2, outlets=outlet),
        "noutlets=1 but outlets has 0 row": dict(noutlets=1),
        "ntables=1 but tables has 0 row": dict(ntables=1),
        "carries a value past .* a boundname": dict(packagedata=[[0, 101.0, 1, "pond"]]),
        r"perioddata\[0\] record \[0, 'RAINFALL'\] holds an array": dict(
            perioddata={0: [[0, "RAINFALL", np.array([0.1])]]}),
    }
    for match, overrides in cases.items():
        with pytest.raises(ValueError, match=match):
            mf.lak.flopy(**{**records, **overrides})
    with pytest.raises(ValueError, match=r"(?s)nlakes=3 .*\n  - ntables=1"):
        mf.lak.flopy(**records, nlakes=3, ntables=1)


def test_mf_lak_flopy_counts_auxiliary_names_in_every_shape_flopy_writes(tmp_path):
    """``len(auxiliary)`` counted ``[("a", "b")]`` and FloPy's own ``get_data()`` recarray as
    ONE auxiliary, so a row carrying both values was refused as a boundname -- records FloPy
    writes ``AUXILIARY a b`` for and MF6 runs (``docs/reproducers/lak_docstring_claims.py``)."""

    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(tmp_path))
    flopy.mf6.ModflowTdis(sim, nper=1)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, top=8.0, botm=-10.0)
    records = dict(packagedata=[[0, 9.0, 1, 5.0, 6.0]], perioddata={0: []},
                   connectiondata=[[0, 0, (0, 1, 1), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]])
    flopys_own = flopy.mf6.ModflowGwflak(gwf, nlakes=1, auxiliary=["a", "b"], **records).auxiliary.get_data()
    assert len(flopys_own) == 1, "the shape that was counted as one name"
    shapes = (["a", "b"], [("a", "b")], ("a", "b"), ["auxiliary", "a", "b"], [("AUXILIARY", "a", "b")],
              np.array(["a", "b"]), flopys_own)
    for auxiliary in shapes:
        assert mf.lak.flopy(**records, auxiliary=auxiliary).options["packagedata"] == records["packagedata"]
    boundname = {**records, "packagedata": [[0, 9.0, 1, 5.0, 6.0, "pond"]]}
    for auxiliary in ([("a", "b")], flopys_own):
        with pytest.raises(ValueError, match=r"nlakeconn, aux, aux\] -- a boundname"):
            mf.lak.flopy(**boundname, auxiliary=auxiliary)
    assert mf.lak.flopy(**boundname, auxiliary=[("a", 2)]), "not all names: left to FloPy"


# ------------------------------------------------------------------------------ MF6


def _table_run(ws: Path, table, stages) -> tuple[bool, pd.DataFrame | None, str]:
    """DIS 3 x 3 of 10 ft cells, one lake on the centre cell with ``table``, held at each stage.

    One steady period per stage in ``stages``; returns (normal termination, the lake's
    STAGE/VOLUME/SURFACE-AREA observations, MF6's first error line).
    """

    shutil.rmtree(ws, ignore_errors=True)
    ws.mkdir(parents=True)
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=len(stages), perioddata=[(1.0, 1, 1.0)] * len(stages))
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10, delc=10, top=8.0, botm=-10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 5.0]])
    flopy.mf6.ModflowGwflak(
        gwf, nlakes=1, ntables=1, tables=[[0, "t.lak.tab"]],
        packagedata=[[0, float(stages[0]), 1]],
        connectiondata=[[0, 0, (0, 1, 1), "VERTICAL", 0.001, 0.0, 0.0, 0.0, 0.0]],
        perioddata={i: [[0, "STATUS", "CONSTANT"], [0, "STAGE", float(s)]] for i, s in enumerate(stages)},
        observations={"lak.obs.csv": [(n.lower().replace("-", ""), n, 1)
                                      for n in ("STAGE", "VOLUME", "SURFACE-AREA")]},
    )
    flopy.mf6.ModflowUtllaktab(gwf, nrow=len(table), ncol=3, table=[list(r) for r in table],
                               filename="t.lak.tab", pname="tab")
    sim.write_simulation(silent=True)
    proc = subprocess.run(["mf6"], cwd=ws, capture_output=True, text=True)
    lst = (ws / "mfsim.lst").read_text() if (ws / "mfsim.lst").exists() else ""
    error = lst.split("ERROR REPORT:")[1].strip().splitlines()[0].strip() if "ERROR REPORT:" in lst else ""
    ok = proc.returncode == 0 and "Normal termination" in lst
    return ok, (pd.read_csv(ws / "lak.obs.csv") if ok else None), error


@pytest.mark.slow
def test_lake_table_checks_are_mf6s_own(tmp_path):
    """Each of MF6's three rules, both ways: LakeTable refuses exactly what MF6 refuses."""

    cases = {
        "equal volume": ([(8.0, 0.0, 0.0), (9.0, 0.0, 0.0), (10.0, 5.0, 4.0)], "VOLUME ENTRY"),
        "falling area": ([(8.0, 0.0, 5.0), (9.0, 1.0, 4.0), (10.0, 5.0, 4.0)], "SURFACE AREA ENTRY"),
        "equal stage": ([(8.0, 0.0, 0.0), (8.0, 1.0, 4.0), (10.0, 5.0, 4.0)], "STAGE ENTRY"),
        "equal area": ([(8.0, 0.0, 4.0), (9.0, 1.0, 4.0), (10.0, 5.0, 4.0)], None),
    }
    for label, (rows, refused) in cases.items():
        ok, _obs, error = _table_run(tmp_path / label.replace(" ", "_"), rows, [9.0])
        assert ok == (refused is None) and (refused or "") in error, (label, error)
        if refused is None:
            assert LakeTable(rows).rows == tuple(rows)
        else:
            with pytest.raises(ValueError):
                LakeTable(rows)
    # What LakeTable refuses BEYOND MF6's checks, and what MF6 does with it instead.
    ok, obs, _ = _table_run(tmp_path / "negative", [(8.0, -5.0, 0.0), (9.0, 1.0, 4.0)], [8.5])
    assert ok and obs["VOLUME"].iloc[-1] == -2.0, "MF6 runs a negative volume as given"
    ok, _obs, error = _table_run(tmp_path / "nan", [(8.0, 0.0, 0.0), (9.0, np.nan, 4.0)], [8.5])
    assert not ok and "Could not read variable of type DOUBLE PRECISION" in error, "FloPy cut the row"


@pytest.mark.slow
def test_the_kept_zero_row_gives_mf6_the_requested_curve(bed, tmp_path):
    """Stages 6..11 over a bed of 8..10: MF6 cannot read the full list (four zero rows),
    and reads the kept table as the full list's curve at every stage -- below the kept
    first row too, because MF6 reads a lower stage as the first row."""

    requested = [6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0]
    full = [(s, float(v), float(a)) for s, v, a in zip(
        requested, [0, 0, 0, 0, 1, 4, 8], [0, 0, 0, 0, 1, 3, 4], strict=True)]
    kept = LakeTableBuilder(dem=bed, footprint=box(0, 0, 2, 2), stages=requested).build().rows
    assert kept == tuple(full[3:])
    ok, _obs, error = _table_run(tmp_path / "full", full, [9.0])
    assert not ok and "VOLUME ENTRY" in error, "the premise: MF6 refuses the full list"
    held = [6.5, 7.25, 8.0, 8.5, 9.5, 10.75]
    ok, obs, error = _table_run(tmp_path / "kept", kept, held)
    assert ok, error
    stages, volumes, areas = zip(*full, strict=True)
    assert np.allclose(obs["VOLUME"], np.interp(held, stages, volumes), atol=1e-9)
    assert np.allclose(obs["SURFACEAREA"], np.interp(held, stages, areas), atol=1e-9)


@pytest.mark.slow
def test_an_array_forcing_now_runs_in_mf6(toy, tmp_path):
    """A rainfall and an outlet rate given as arrays: written per period, and MF6 runs.

    Before ledger 318 each was written whole into every period and MF6 died reading it
    (a segmentation fault)."""

    grid = toy["grid"]
    ws = tmp_path / "array"
    ws.mkdir()
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdisv(gwf, nlay=2, top=12.0, botm=[8.0, 0.0], **grid.get_disv_gridprops())
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1e-5, iconvert=1, transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 10.0]])
    lak = _lak(toy, nper=2, rainfall=np.array([0.1, 0.2]),
               outlets=[LakeOutlet("pond", rate=np.array([-0.5, -1.0]))],
               tables={"pond": LakeTableBuilder(area=4.0, bottom=9.0, top=13.0, stage_step=0.5)},
               observations={"lak.obs.csv": [("rain", "RAINFALL", "pond"), ("out", "EXT-OUTFLOW", "pond")]})
    lak.build(gwf)
    sim.write_simulation(silent=True)
    rates = [line.split() for line in (ws / "m.lak").read_text().splitlines()
             if "RAINFALL" in line or "RATE" in line]
    assert [(row[1], float(row[2])) for row in rates] == [
        ("RAINFALL", 0.1), ("RATE", -0.5), ("RAINFALL", 0.2), ("RATE", -1.0)
    ], rates
    proc = subprocess.run(["mf6"], cwd=ws, capture_output=True, text=True)
    assert proc.returncode == 0 and "Normal termination" in (ws / "mfsim.lst").read_text()
    obs = pd.read_csv(ws / "lak.obs.csv")
    assert np.allclose(obs["RAIN"], [0.4, 0.8]), "rate x the table's 4.0 ft2"
    assert np.allclose(obs["OUT"].abs(), [0.5, 1.0])
