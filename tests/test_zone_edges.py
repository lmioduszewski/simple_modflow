"""A zone boundary drawn through a row of cell centroids left that row in no zone.

Both zone assignments -- ``ZoneSpec.from_polygons`` for PEST and
``surfaces._cells_in_zone`` for a layer surface -- asked whether a cell's centroid
was ``within`` the zone polygon, which is FALSE for a centroid lying exactly on the
polygon's edge. A grid's centroids are where its data is, so that edge is exactly
where a zone boundary gets drawn. Measured on a 4x4 grid of 25 ft cells:

* two zones split at ``y=37.5``, the line through the row-1 centroids, gave those
  four cells the ``nodata_zone`` label 0, while a split at ``y=37.4`` gave them
  zone 2. pyEMU gives a zone ``<= 0`` no multiplier, so that row of the model went
  uncalibrated under a definition that names it;
* a lone zone ``box(0, 0, 12.5, 100)``, whose east edge runs along the column-0
  centroids, selected NO cells at all -- an all-nodata ``ZoneSpec``, in silence;
* ``surfaces._cells_in_zone(grid, box(0, 0, 100, 37.5))`` selected ``[0, 1, 2, 3]``
  and not ``[4, 5, 6, 7]``, which feeds ``Surface.where``, hence botm, and through
  pinch-out idomain.

The centroid rule itself is deliberate (ledgers 118 and 213) and is kept: a zone
labels a WHOLE cell, so its centre decides, not the area rule ``cell_shares``
applies to a boundary condition. Only the boundary case moves -- a zone is CLOSED,
and a centroid on its edge is in it. Where two zones share that edge the FIRST ROW
wins, which is what both docstrings always claimed; ``covered_by`` alone does not
deliver it, because the winner then comes out of the spatial index's own order.
"""

from __future__ import annotations

import inspect
import os
import warnings
from types import SimpleNamespace

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Point, Polygon, box

from myflopy import surfaces as surfaces_module
from myflopy.modflow.mf6.grid import geometry as geometry_module
from myflopy.modflow.mf6.grid.geometry import centroid_zones
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.pest import zones as zones_module
from myflopy.modflow.mf6.pest.zones import ZoneSpec
from myflopy.surfaces import Surface, _cells_in_zone

CRS = "EPSG:2927"
# Row-major from the bottom left: cell 5 spans 25..50 in both axes, and every
# centroid sits on a multiple of 12.5.
LOWER, UPPER = box(0.0, 0.0, 100.0, 37.5), box(0.0, 37.5, 100.0, 100.0)
OFF_GRID = box(500.0, 500.0, 600.0, 600.0)
NOTHING_LANDED = "No cell centroid falls in any zone"


@pytest.fixture(scope="module")
def grid() -> VoronoiGridPlus:
    """4x4 square cells, 25 ft on a side, origin at (0, 0)."""

    n, d = 4, 25.0
    verts = np.array(
        [[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float
    )
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


@pytest.fixture(scope="module")
def model(grid) -> SimpleNamespace:
    """A model-like object: `ZoneSpec.from_polygons` reads only `model.vor`."""

    return SimpleNamespace(vor=grid)


def _frame(geometries, labels, crs=CRS) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame({"z": list(labels)}, geometry=list(geometries), crs=crs)


def _labels(model, geometries, labels, **kwargs) -> list[int]:
    frame = _frame(geometries, labels)
    return ZoneSpec.from_polygons(frame, model, column="z", **kwargs).values.tolist()


# ------------------------------------------------------------------ the two measured cases
def test_a_split_through_a_row_of_centroids_leaves_no_cell_unzoned(model):
    """The measured case: zones split at y=37.5 gave cells 4-7 the label 0, which
    pyEMU reads as 'no multiplier', so the row was silently uncalibrated."""

    values = _labels(model, [LOWER, UPPER], [1, 2])

    assert values == [1] * 8 + [2] * 8
    assert 0 not in values, "a zone boundary on a row of centroids left that row out"
    spec = ZoneSpec.from_polygons(_frame([LOWER, UPPER], [1, 2]), model, column="z")
    assert spec.ids == [1, 2], "the nodata label is no longer one of the zones"


def test_a_split_just_off_the_centroids_is_unchanged(model):
    """The control the defect was measured against: at y=37.4 nothing was ever
    wrong, and the fix must not move it."""

    assert _labels(model, [box(0.0, 0.0, 100.0, 37.4), box(0.0, 37.4, 100.0, 100.0)],
                   [1, 2]) == [1] * 4 + [2] * 12


def test_a_lone_zone_along_a_column_of_centroids_had_selected_nothing(model):
    """`box(0, 0, 12.5, 100)`'s east edge IS the column-0 centroid line, so the
    whole definition came back as nodata -- one zone, nothing parameterized."""

    values = _labels(model, [box(0.0, 0.0, 12.5, 100.0)], [5])

    assert values == [5, 0, 0, 0, 5, 0, 0, 0, 5, 0, 0, 0, 5, 0, 0, 0]
    assert set(values) == {0, 5}


def test_a_centroid_on_the_outer_edge_of_a_lone_zone_is_inside_it(model):
    """No second zone to share the edge with: the rule is that a zone is CLOSED,
    not that ties are resolved, so the boundary belongs to the only zone there is."""

    assert _labels(model, [LOWER], [5]) == [5] * 8 + [0] * 8
    # ...including a centroid on a CORNER of the zone rather than on a side.
    corner = box(12.5, 12.5, 100.0, 100.0)
    assert Point(12.5, 12.5).covered_by(corner) and not Point(12.5, 12.5).within(corner)
    assert _labels(model, [corner], [5])[0] == 5


def test_two_zones_sharing_an_edge_put_every_centroid_in_exactly_one(model):
    """Exactly one, and deterministically the FIRST ROW -- so swapping the rows
    swaps the shared row's label, and nothing else."""

    lower_first = _labels(model, [LOWER, UPPER], [1, 2])
    upper_first = _labels(model, [UPPER, LOWER], [1, 2])

    assert lower_first == [1] * 8 + [2] * 8
    assert upper_first == [2] * 4 + [1] * 12
    # Same partition either way: 16 cells, every one labelled, none twice.
    for values in (lower_first, upper_first):
        assert len(values) == 16 and 0 not in values


def test_the_shared_edge_goes_to_the_first_row_not_to_the_index_s_order(model):
    """`covered_by` alone would not have been deterministic. Measured on this very
    frame: the join's first match gave cells 4-7 the SECOND row's label, because
    the spatial index answers in tree order."""

    rows, labels = [UPPER, LOWER], [1, 2]
    positions = centroid_zones(model.vor.gdf_vorPolys, _frame(rows, labels))

    assert positions.tolist() == [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert _labels(model, rows, labels)[4:8] == [1, 1, 1, 1]


def test_a_centroid_in_no_zone_keeps_the_nodata_label(model):
    """Unchanged, and the reason the defect was silent rather than loud."""

    with pytest.warns(UserWarning, match=NOTHING_LANDED):
        assert _labels(model, [OFF_GRID], [9]) == [0] * 16
    with pytest.warns(UserWarning, match=NOTHING_LANDED):
        assert _labels(model, [OFF_GRID], [9], nodata_zone=-1) == [-1] * 16
    # A zone reaching part of the grid leaves the rest at nodata, and says NOTHING:
    # most zone layers cover part of a model on purpose, so only a total miss warns.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert _labels(model, [box(0.0, 0.0, 100.0, 20.0)], [4]) == [4] * 4 + [0] * 12
    assert not [w for w in caught if NOTHING_LANDED in str(w.message)]


def test_the_area_rule_is_still_not_the_zone_rule(model):
    """The centroid rule is deliberate (ledgers 118, 213), and this is the test that
    fails if someone later reaches for `cell_shares`: `box(0, 0, 100, 37.4)` covers
    49.6% of cells 4-7 by AREA and none of their centres, so it does not label
    them."""

    zone = box(0.0, 0.0, 100.0, 37.4)
    polygons = model.vor.gdf_vorPolys
    share = zone.intersection(polygons.geometry.loc[5]).area / polygons.geometry.loc[5].area
    assert 0.49 < share < 0.50

    assert _labels(model, [zone], [1])[4:8] == [0, 0, 0, 0]
    assert [cell for cell, _ in geometry_module.cell_shares(polygons, zone)][4:8] == [4, 5, 6, 7]


# ------------------------------------------------------------------ a zone that is not one
def test_a_row_whose_label_is_null_is_not_a_zone(model, tmp_path):
    """A NULL zone id in the file. `int(nan)` is -9223372036854775808, which pyEMU
    would calibrate like any other label, so the row is taken out of the running
    instead and its cells fall to whatever else covers them."""

    path = tmp_path / "zones.gpkg"
    gpd.GeoDataFrame({"z": [None, 2]}, geometry=[LOWER, UPPER], crs=CRS).to_file(
        path, layer="z", driver="GPKG"
    )
    # What `gpd.read_file` hands back for a NULL integer attribute.
    assert np.isnan(gpd.read_file(path)["z"].to_numpy()[0])

    spec = ZoneSpec.from_polygons(path, model, column="z")

    # Cells 0-3 are covered only by the unlabelled row; 4-7 by both, and the
    # labelled one is the only zone in the running.
    assert spec.values.tolist() == [0] * 4 + [2] * 12
    assert spec.ids == [0, 2]
    assert min(spec.values) >= 0, "a missing label must not become a real zone id"


def test_a_null_label_in_any_dtype_is_treated_the_same(model):
    """float64 NaN is what a file gives; a nullable `Int64` column carries `pd.NA`,
    whose `!= itself` cannot even be evaluated as a bool."""

    for column in (pd.array([pd.NA, 2], dtype="Int64"), [np.nan, 2.0], [None, 2]):
        frame = gpd.GeoDataFrame({"z": column}, geometry=[LOWER, UPPER], crs=CRS)
        values = ZoneSpec.from_polygons(frame, model, column="z").values
        assert values.tolist() == [0] * 4 + [2] * 12, column


def test_a_zone_layer_that_reaches_no_cell_says_so(model):
    """A frame in the wrong place, or silently in another CRS, produced an
    all-nodata map -- and array targets skip nodata, so nothing was parameterized.
    The extents-and-CRS sentence is the shared one (ledger 211's wording)."""

    with pytest.warns(UserWarning, match=NOTHING_LANDED) as caught:
        ZoneSpec.from_polygons(_frame([OFF_GRID], [9]), model, column="z")

    message = str(caught[0].message)
    assert "500.0" in message and "100.0" in message, message
    assert "zone polygons" in message and CRS in message, message


def test_a_zone_layer_with_no_rows_at_all_warns_rather_than_raising(model):
    """A degenerate layer still has to answer; `total_bounds` of nothing is NaN."""

    empty = gpd.GeoDataFrame({"z": []}, geometry=[], crs=CRS)
    with pytest.warns(UserWarning, match=NOTHING_LANDED):
        assert ZoneSpec.from_polygons(empty, model, column="z").values.tolist() == [0] * 16


# ------------------------------------------------------------------ the other routes
def test_the_array_route_is_untouched(model):
    assert ZoneSpec.from_array(np.arange(16) % 3).values.tolist() == [
        0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0
    ]


def test_the_raster_route_is_untouched(model, tmp_path):
    """`from_raster` majority-votes the pixels in a cell and falls back to the label
    UNDER the centroid (ledger 213). Pinned to the values HEAD produced.

    A raster label boundary cannot be made to run along a centroid line without
    halving the cell, which is a MAJORITY tie rather than a predicate question --
    which is the point: no polygon predicate is asked anywhere on this route, so the
    closed-zone rule cannot reach it. The companion assertion below is the probe."""

    band = np.repeat(
        np.where(np.arange(4)[:, None] < 2, 3, 1), 4, axis=1
    ).astype("int32")  # top half (y > 50) zone 3, bottom half zone 1
    path = tmp_path / "hsu.tif"
    with rasterio.open(
        path, "w", driver="GTiff", height=4, width=4, count=1, dtype="int32",
        crs=CRS, transform=from_origin(0.0, 100.0, 25.0, 25.0), nodata=-9999,
    ) as destination:
        destination.write(band, 1)

    assert ZoneSpec.from_raster(path, model).values.tolist() == [1] * 8 + [3] * 8

    source = inspect.getsource(zones_module._majority_sample) + inspect.getsource(
        zones_module._centroid_labels
    )
    assert "centroid_zones" not in source and "covered_by" not in source


def test_a_polygon_layer_in_another_crs_still_lands(model):
    """Reprojection moves the boundary OFF the centroids by ~1e-9, so this route was
    never the defect; it still has to label every cell."""

    frame = _frame([LOWER, UPPER], [1, 2]).to_crs("EPSG:4326")
    values = ZoneSpec.from_polygons(frame, model, column="z").values.tolist()

    assert 0 not in values
    assert values[:4] == [1] * 4 and values[8:] == [2] * 8


def test_a_null_geometry_row_is_skipped_rather_than_raising(model):
    """A zone layer read from a file can carry one; the join tolerated it."""

    assert _labels(model, [None, box(0.0, 0.0, 100.0, 20.0)], [7, 8]) == [8] * 4 + [0] * 12


def test_an_unnamed_column_still_numbers_the_rows(model):
    assert _labels(model, [LOWER, UPPER], [1, 2]) == ZoneSpec.from_polygons(
        _frame([LOWER, UPPER], [11, 22]), model
    ).values.tolist()


# ------------------------------------------------------------------ surfaces
def test_a_surface_zone_includes_a_centroid_on_its_edge(grid, tmp_path):
    """`_cells_in_zone` fed `Surface.where`, hence botm, and through pinch-out
    idomain. It selected `[0, 1, 2, 3]` for a zone whose edge ran along the row-1
    centroids."""

    assert np.flatnonzero(_cells_in_zone(grid, LOWER)).tolist() == [0, 1, 2, 3, 4, 5, 6, 7]
    # The control, unchanged.
    assert np.flatnonzero(
        _cells_in_zone(grid, box(0.0, 0.0, 100.0, 37.4))
    ).tolist() == [0, 1, 2, 3]

    path = tmp_path / "zone.gpkg"
    gpd.GeoDataFrame({"n": ["a"]}, geometry=[LOWER], crs=CRS).to_file(
        path, layer="z", driver="GPKG"
    )
    assert np.flatnonzero(_cells_in_zone(grid, path)).tolist() == [0, 1, 2, 3, 4, 5, 6, 7]


def test_a_multi_row_surface_zone_file_is_the_union_of_its_rows(grid, tmp_path):
    """One zone, no labels: a vector file's rows are unioned, so two disjoint
    polygons both count -- there is no first row to prefer, which is why this side
    keeps `union_all()` rather than the row order the PEST side needs."""

    path = tmp_path / "two_rows.gpkg"
    gpd.GeoDataFrame(
        {"n": ["west", "east"]},
        geometry=[box(0.0, 0.0, 12.5, 37.5), box(62.5, 0.0, 100.0, 37.5)],
        crs=CRS,
    ).to_file(path, layer="z", driver="GPKG")

    # Column 0 from the west row (its east edge is the centroid line) and columns
    # 2-3 from the east row, rows 0 and 1 (whose top edge is the centroid line).
    assert np.flatnonzero(_cells_in_zone(grid, path)).tolist() == [0, 2, 3, 4, 6, 7]


def test_the_elevation_a_zoned_surface_hands_a_cell_moves_with_it(grid):
    """The consequence: those four cells took the OUTSIDE elevation, 2.0 where the
    zone says 15.0 -- a 13 ft error in a layer bottom."""

    surface = Surface.flat(15.0).within(LOWER, outside=Surface.flat(2.0))

    assert surface.values(grid, None).tolist() == [15.0] * 8 + [2.0] * 8


# ------------------------------------------------------------------ one rule, one copy
def test_both_zone_callers_share_the_one_centroid_rule(grid):
    """Two copies of the predicate were two chances to fix one of them."""

    assert zones_module.centroid_zones is geometry_module.centroid_zones
    assert surfaces_module.centroid_zones is geometry_module.centroid_zones


def test_the_shared_rule_answers_row_positions(grid):
    """Positions into `zones`, not labels and not cell ids -- both callers map them
    themselves, one to a label array and one to a boolean mask."""

    polygons = grid.gdf_vorPolys

    assert centroid_zones(polygons, []).tolist() == [-1] * 16
    assert centroid_zones(polygons, [OFF_GRID]).tolist() == [-1] * 16
    # A GeoSeries, a GeoDataFrame and a plain list are all accepted.
    for zones in ([LOWER, UPPER], gpd.GeoSeries([LOWER, UPPER], crs=CRS),
                  _frame([LOWER, UPPER], [1, 2])):
        assert centroid_zones(polygons, zones).tolist() == [0] * 8 + [1] * 8
    # An empty or absent geometry occupies its row without ever winning it -- the
    # spatial index matches neither, so the answer is still by row position, which
    # is what lets a caller blank a row out (`from_polygons` does, for a null label).
    assert centroid_zones(polygons, [Polygon(), None, LOWER]).tolist() == [2] * 8 + [-1] * 8
