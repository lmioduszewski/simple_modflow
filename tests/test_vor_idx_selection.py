"""The two public ``get_vor_idx_*`` grid selectors pick cells by coverage, in the grid's frame.

Both used a bare spatial-index ``intersects`` and answered about cells a feature
only borders: on this 4x4 grid of 25 ft cells a polygon that was exactly cell 5
came back as 9 cells and a point on the vertex shared by 5, 6, 9 and 10 as all
four. ``get_vor_idx_from_geometry`` also forced ``crs=EPSG:2927`` and then compared
in EPSG:4326, so on a UTM grid the same polygon came back as cells 4, 5, 8 and 9,
a query in any other CRS came back EMPTY from the ``_idx`` twin and as six cells
from the dict one, and a naive geometry raised ``Cannot transform naive
geometries`` from one while the other answered (with its own 9 cells -- it read no
CRS at all).

Under an explicit ``predicate`` the two asked the spatial index its two sides
opposite ways round, so one word meant two questions: ``'contains'`` over a box
inside cell 5 gave ``[5]`` from the dict one, which is what
``get_vor_cells_as_series`` answers, and ``[]`` from the ``_idx`` one. Both ask it
the selector's way now, of the cell first, and that parity is pinned here in both
directions by two predicates that are each other's inverse.

Every expectation here is worked out from row/column arithmetic on a regular grid,
not from the code under test; the tests that ask the shared selector do so only to
pin that the three spellings agree.
"""

from __future__ import annotations

import inspect
import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, box

from myflopy.modflow.mf6.grid.geometry import (
    get_vor_idx_from_geometry,
    get_vor_idx_from_geometry_idx,
    selected_cells,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

LOGGER = "myflopy.modflow.mf6.grid.geometry"

N = 4                 # 4 x 4 cells
D = 25.0              # each 25 ft square, so the grid spans 0..100

CELL_5 = box(25.0, 25.0, 50.0, 50.0)          # exactly cell 5
ON_VERTEX = Point(50.0, 50.0)                 # the vertex shared by 5, 6, 9 and 10
ON_FACE = Point(25.0, 37.5)                   # the face between 4 and 5
THROUGH_ROW_0 = LineString([(10.0, 10.0), (90.0, 10.0)])   # inside row 0, no face touched
DIAGONAL = LineString([(0.0, 0.0), (100.0, 100.0)])        # corner to corner
ALONG_FACE = LineString([(0.0, 25.0), (100.0, 25.0)])      # the whole row-0/row-1 face
OFF_GRID = box(500.0, 500.0, 600.0, 600.0)
INSIDE_5 = box(30.0, 30.0, 40.0, 40.0)        # strictly inside cell 5, touching no face
HOLDS_FOUR = box(25.0, 25.0, 75.0, 75.0)      # exactly cells 5, 6, 9 and 10
CORNER_ONLY = LineString([(-10.0, 10.0), (10.0, -10.0)])   # meets the grid at (0, 0) alone


def _grid(crs: str = "EPSG:2927", x0: float = 0.0, y0: float = 0.0) -> VoronoiGridPlus:
    """N x N square cells of side D, row-major from the bottom left of (x0, y0)."""

    verts = np.array(
        [[x0 + i * D, y0 + j * D] for j in range(N + 1) for i in range(N + 1)], dtype=float
    )
    iverts = [
        [j * (N + 1) + i, (j + 1) * (N + 1) + i, (j + 1) * (N + 1) + i + 1, j * (N + 1) + i + 1]
        for j in range(N)
        for i in range(N)
    ]
    xcyc = np.array(
        [[x0 + (i + 0.5) * D, y0 + (j + 0.5) * D] for j in range(N) for i in range(N)], dtype=float
    )
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc, crs=crs)


@pytest.fixture(scope="module")
def grid() -> VoronoiGridPlus:
    return _grid()


def _frame(geometry, crs: str = "EPSG:2927") -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame({"name": ["feature"]}, geometry=[geometry], crs=crs)


def _four_answers(vor, frame, **kwargs) -> dict[str, list[int]]:
    """What each of the four public entry points says, flattened the same way.

    Two free functions and the two ``VoronoiGridPlus`` methods that forward to
    them: the defaults were restated on the methods, so fixing the free functions
    alone would have reached nothing a user calls.
    """

    return {
        "free dict": selected_cells(
            list(get_vor_idx_from_geometry(vor, gdf_to_query=frame, **kwargs).values())
        ),
        "free idx": get_vor_idx_from_geometry_idx(vor, frame, 0, **kwargs),
        "method dict": selected_cells(
            list(vor.get_vor_idx_from_geometry(gdf_to_query=frame, **kwargs).values())
        ),
        "method idx": vor.get_vor_idx_from_geometry_idx(frame, 0, **kwargs),
    }


def _agreed(vor, geometry, frame_crs: str = "EPSG:2927", **kwargs) -> list[int]:
    """The one answer all four entry points give, or a failure naming who differs."""

    answers = _four_answers(vor, _frame(geometry, frame_crs), **kwargs)
    assert len(set(map(tuple, answers.values()))) == 1, answers
    return answers["free idx"]


# ---------------------------------------------------------------- coverage


def test_a_cell_sized_polygon_reaches_that_one_cell(grid):
    """9 cells before, from both: the eight neighbours it only borders came too."""

    assert _agreed(grid, CELL_5) == [5]


def test_a_point_on_a_vertex_goes_to_the_lowest_cell_sharing_it(grid):
    """[5, 6, 9, 10] before: one well became four."""

    assert _agreed(grid, ON_VERTEX) == [5]


def test_a_point_on_a_face_goes_to_the_lower_of_the_two(grid):
    """The two functions DISAGREED here: [4, 5] from the index one, taking both
    cells, and [5] from the dict one, which reached that through a reprojection to
    EPSG:4326. One cell holds a point, and it is the lower of the pair."""

    assert _agreed(grid, ON_FACE) == [4]


def test_a_line_reaches_the_cells_it_runs_through(grid):
    assert _agreed(grid, THROUGH_ROW_0) == [0, 1, 2, 3]


def test_a_diagonal_line_crosses_four_cells_not_ten(grid):
    """It grazes the corner of six more, touching them at a single point."""

    assert _agreed(grid, DIAGONAL) == [0, 5, 10, 15]


def test_a_line_along_a_face_reaches_both_sides(grid):
    """Ledger 212: an unsplit field reaches both cells of the pair."""

    assert _agreed(grid, ALONG_FACE) == [0, 1, 2, 3, 4, 5, 6, 7]


def test_the_any_touch_answer_is_still_reachable(grid):
    """The escape hatch, for a caller that wants topology rather than coverage --
    and one ORIENTATION for it now, where the pair asked the index opposite ways."""

    assert _agreed(grid, CELL_5, predicate="intersects") == [0, 1, 2, 4, 5, 6, 8, 9, 10]
    assert _agreed(grid, ON_VERTEX, predicate="intersects") == [5, 6, 9, 10]


@pytest.mark.parametrize(
    ("geometry", "predicate", "expected"),
    [
        (INSIDE_5, "contains", [5]),
        (INSIDE_5, "covers", [5]),
        (HOLDS_FOUR, "within", [5, 6, 9, 10]),
    ],
    ids=["contains", "covers", "within"],
)
def test_an_asymmetric_predicate_is_asked_of_the_cell_first(grid, geometry, predicate, expected):
    """The orientation, pinned in BOTH directions by two predicates that are each
    other's inverse. ``get_vor_cells_as_series`` joins the CELLS to the query, so
    ``'contains'`` means "cells containing the feature"; the index selector asked
    ``contains(feature, cell)`` and answered the opposite question under the same
    name. Only a symmetric predicate -- the old default -- hid it, and one entry
    point agreeing with the selector while its twin disagreed is the defect this
    file exists to close."""

    frame = _frame(geometry)
    expected_from_selector = selected_cells(
        grid.get_vor_cells_as_series(frame, predicate=predicate, name_field="name")
    )

    assert expected_from_selector == expected
    for label, answer in _four_answers(grid, frame, predicate=predicate).items():
        assert answer == expected, label


def test_the_predicate_answer_is_ascending_whatever_order_the_grid_is_held_in():
    """A grid frame need not be in cell order -- and the spatial index answers in
    the frame's order, so the cells came back descending from a reversed one. Both
    branches sort; ``polys.index`` is also what turns a row POSITION into a cell."""

    reversed_grid = _grid()
    reversed_grid.gdf_vorPolys = reversed_grid.gdf_vorPolys.iloc[::-1]
    strip = box(10.0, 10.0, 90.0, 40.0)

    assert reversed_grid.gdf_vorPolys.index.tolist()[:3] == [15, 14, 13]
    assert _agreed(reversed_grid, strip) == [0, 1, 2, 3, 4, 5, 6, 7]
    assert _agreed(reversed_grid, strip, predicate="intersects") == [0, 1, 2, 3, 4, 5, 6, 7]
    # The dict function's own lists, unflattened, so the ORDER is what is compared
    # and not what `selected_cells` sorted on the way out.
    frame = _frame(strip)
    assert get_vor_idx_from_geometry(reversed_grid, gdf_to_query=frame) == {
        0: [0, 1, 2, 3, 4, 5, 6, 7],
    }
    assert get_vor_idx_from_geometry(
        reversed_grid, gdf_to_query=frame, predicate="intersects"
    ) == {0: [0, 1, 2, 3, 4, 5, 6, 7]}


# ---------------------------------------------------------------- the grid's own frame


def test_a_geometry_in_another_crs_answers_about_the_right_cells(grid, caplog):
    """Before: EMPTY from the index function (degrees measured against feet) and
    six cells from the dict one. The reconcile says what it did -- that wording is
    shared with ``get_vor_cells_as_series``, which is the only place it was."""

    in_degrees = gpd.GeoSeries([CELL_5], crs=grid.crs).to_crs("EPSG:4326").iloc[0]

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        assert _agreed(grid, in_degrees, frame_crs="EPSG:4326") == [5]

    messages = [record.getMessage() for record in caplog.records]
    assert any(
        "geometries are in EPSG:4326 but the grid is in EPSG:2927" in message
        and "reprojecting them to match the grid" in message
        for message in messages
    ), messages


def test_the_shared_reconcile_says_the_same_thing_for_the_selector(grid, caplog):
    """One copy of that warning: ``get_vor_cells_as_series`` now delegates to it."""

    in_degrees = _frame(
        gpd.GeoSeries([CELL_5], crs=grid.crs).to_crs("EPSG:4326").iloc[0], crs="EPSG:4326"
    )

    with caplog.at_level(logging.WARNING):
        grid.get_vor_cells_as_series(in_degrees, name_field="name")

    assert any(
        "geometries are in EPSG:4326 but the grid is in EPSG:2927" in record.getMessage()
        for record in caplog.records
    ), [record.getMessage() for record in caplog.records]


def test_crs_cannot_act_over_a_query_that_carries_its_own_frame(grid, caplog):
    """It names the frame of a NAIVE query, so over a georeferenced one it is
    ignored -- out loud. Silently ignored, it reads as the frame the query will be
    measured in, which is the defect this whole fix is about."""

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        answer = get_vor_idx_from_geometry(
            grid, gdf_to_query=_frame(CELL_5), crs="EPSG:4326", name_col="name"
        )

    assert answer == {"feature": [5]}
    assert any(
        "EPSG:4326" in record.getMessage() and "ignored" in record.getMessage()
        for record in caplog.records
    ), [record.getMessage() for record in caplog.records]


def test_crs_naming_the_query_s_own_frame_is_not_complained_about(grid, caplog):
    """Redundant is not wrong: a caller passing the frame the query already carries
    gets no warning."""

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        answer = get_vor_idx_from_geometry(grid, gdf_to_query=_frame(CELL_5), crs=grid.crs)

    assert answer == {0: [5]}
    assert [record.getMessage() for record in caplog.records] == []


def test_a_grid_outside_the_default_crs_is_not_answered_about_in_it():
    """A UTM grid: cells 4, 5, 8 and 9 before, for a polygon that was cell 5."""

    utm = _grid(crs="EPSG:32610", x0=5.0e5, y0=5.0e6)
    cell_5 = box(5.0e5 + 25.0, 5.0e6 + 25.0, 5.0e5 + 50.0, 5.0e6 + 50.0)

    assert _agreed(utm, cell_5, frame_crs="EPSG:32610") == [5]


def test_a_query_carrying_no_crs_is_read_as_the_grid_s(grid):
    """The dict function raised ``Cannot transform naive geometries`` on every one
    of these, while its twin answered."""

    naive = gpd.GeoDataFrame({"name": ["feature"]}, geometry=[CELL_5])
    answers = _four_answers(grid, naive)

    assert naive.crs is None
    assert answers == dict.fromkeys(answers, [5])


def test_crs_names_the_frame_of_a_naive_query(grid):
    """Its one surviving job: a file of degrees with no CRS written into it."""

    in_degrees = gpd.GeoSeries([CELL_5], crs=grid.crs).to_crs("EPSG:4326").iloc[0]
    naive = gpd.GeoDataFrame({"name": ["feature"]}, geometry=[in_degrees])

    assert get_vor_idx_from_geometry(
        grid, gdf_to_query=naive, crs="EPSG:4326", name_col="name"
    ) == {"feature": [5]}
    # The `_idx` twin has no `crs=` of its own -- signature kept -- so the same
    # naive degrees are read as the grid's feet and reach nothing.
    assert get_vor_idx_from_geometry_idx(grid, naive, 0) == []


def test_a_vector_file_is_read_in_its_own_crs(grid, tmp_path):
    """Not reprojected into the default one and then compared in lat/lon."""

    path = tmp_path / "query.gpkg"
    gpd.GeoDataFrame(
        {"name": ["cell_5"]}, geometry=[CELL_5], crs=grid.crs
    ).to_crs("EPSG:4326").to_file(path, driver="GPKG")

    assert get_vor_idx_from_geometry(grid, shp_to_query=path, name_col="name") == {"cell_5": [5]}


@pytest.mark.parametrize("spelling", [4326, "epsg:4326", "EPSG:4326", None])
def test_crs_latlon_is_read_as_a_crs_not_as_a_string(grid, spelling):
    """``4326`` is an int geopandas accepts and worked before; so does lower case.
    A string comparison would refuse the default value spelled another way."""

    assert get_vor_idx_from_geometry(
        grid, gdf_to_query=_frame(CELL_5), crs_latlon=spelling
    ) == {0: [5]}


def test_crs_latlon_cannot_act_and_says_so(grid):
    """It was the comparison frame -- the defect itself. Refused rather than
    accepted and ignored."""

    with pytest.raises(ValueError, match="crs_latlon"):
        get_vor_idx_from_geometry(grid, gdf_to_query=_frame(CELL_5), crs_latlon="EPSG:3857")


def test_the_default_crs_latlon_the_method_passes_is_accepted(grid):
    """The method repeats the literal, so a drift between the two would raise on
    every call through it."""

    default = inspect.signature(VoronoiGridPlus.get_vor_idx_from_geometry).parameters["crs_latlon"]

    assert grid.get_vor_idx_from_geometry(gdf_to_query=_frame(CELL_5)) == {0: [5]}
    assert get_vor_idx_from_geometry(
        grid, gdf_to_query=_frame(CELL_5), crs_latlon=default.default
    ) == {0: [5]}


# ---------------------------------------------------------------- keys, misses, defaults


def test_a_geometry_off_the_grid_comes_back_empty_and_is_named(grid, caplog):
    """Silently empty before. The key stays -- that shape is what callers index --
    and when NOTHING landed, the warning carries both extents and both CRSs, which
    is what ledger 211 says an empty result is acted on with. The failure this site
    actually produced was a CRS mismatch, so the diagnosis belongs in the one case
    where every feature missed."""

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        answer = get_vor_idx_from_geometry(
            grid, gdf_to_query=_frame(OFF_GRID), name_col="name"
        )

    assert answer == {"feature": []}
    assert get_vor_idx_from_geometry_idx(grid, _frame(OFF_GRID), 0) == []
    message = " ".join(record.getMessage() for record in caplog.records)
    assert "'feature'" in message
    assert "outside the grid" in message
    assert "500.0, 500.0, 600.0, 600.0" in message      # where the feature is
    assert "0.0, 0.0, 100.0, 100.0" in message          # where the grid is
    assert "both are in EPSG:2927" in message


def test_a_total_miss_names_the_frame_the_features_came_in(grid, caplog):
    """The sentence's job: the features' OWN bounds, in the CRS it names, because a
    mismatch is two boxes and reprojecting one of them first hides it."""

    away = gpd.GeoSeries([OFF_GRID], crs=grid.crs).to_crs("EPSG:4326")

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        get_vor_idx_from_geometry(
            grid, gdf_to_query=_frame(away.iloc[0], crs="EPSG:4326"), name_col="name"
        )

    message = " ".join(record.getMessage() for record in caplog.records)
    assert "reprojected from EPSG:4326 to the grid's EPSG:2927" in message
    # Degrees -- the frame the features were handed over in, not the grid's feet.
    assert "-126.9, 45.2" in message


def test_a_line_meeting_the_grid_at_one_corner_says_it_is_on_no_cell(grid, caplog):
    """Not "outside the grid" -- it reaches the grid, and covers nothing of it. The
    wording is ``miss_reason``'s, asked for a cell."""

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        answer = get_vor_idx_from_geometry(
            grid, gdf_to_query=_frame(CORNER_ONLY), name_col="name"
        )

    assert answer == {"feature": []}
    assert any("on no cell" in record.getMessage() for record in caplog.records), [
        record.getMessage() for record in caplog.records
    ]


def test_a_long_list_of_misses_is_truncated_and_counted(grid, caplog):
    """Ten named, the rest counted: the convention ``_measured_cells`` reports by."""

    away = gpd.GeoDataFrame(
        {"name": [f"m{i}" for i in range(12)]},
        geometry=[box(500.0 + 10 * i, 500.0, 505.0 + 10 * i, 505.0) for i in range(12)],
        crs=grid.crs,
    )

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        answer = get_vor_idx_from_geometry(grid, gdf_to_query=away, name_col="name")

    message = " ".join(record.getMessage() for record in caplog.records)
    assert list(answer) == [f"m{i}" for i in range(12)]
    assert "'m9'" in message
    assert "'m10'" not in message
    assert "and 2 more" in message


def test_a_partial_miss_keeps_the_rest(grid, caplog):
    features = gpd.GeoDataFrame(
        {"name": ["inside", "away"]}, geometry=[CELL_5, OFF_GRID], crs=grid.crs
    )

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        answer = get_vor_idx_from_geometry(grid, gdf_to_query=features, name_col="name")

    assert answer == {"inside": [5], "away": []}
    assert any("'away'" in record.getMessage() for record in caplog.records)


def test_a_filtered_frame_keeps_its_own_labels(grid):
    """``gdf.iloc[[2]]`` carries the label 2, which was read as row 2 of a one-row
    frame: ``IndexError: single positional indexer is out-of-bounds``."""

    features = gpd.GeoDataFrame(
        {"name": ["a", "b", "c"]},
        geometry=[box(26.0, 26.0, 49.0, 49.0), box(51.0, 51.0, 74.0, 74.0), box(76.0, 76.0, 99.0, 99.0)],
        crs=grid.crs,
    )
    third = features.iloc[[2]]

    assert third.index.tolist() == [2]
    assert get_vor_idx_from_geometry(grid, gdf_to_query=third) == {2: [15]}
    assert get_vor_idx_from_geometry(grid, gdf_to_query=third, name_col="name") == {"c": [15]}
    assert get_vor_idx_from_geometry_idx(grid, third, 0) == [15]
    assert get_vor_idx_from_geometry_idx(grid, features, 2) == [15]


def test_every_feature_keeps_its_own_key(grid):
    features = gpd.GeoDataFrame(
        {"name": ["five", "vertex"]}, geometry=[CELL_5, ON_VERTEX], crs=grid.crs
    )

    assert get_vor_idx_from_geometry(grid, gdf_to_query=features, name_col="name") == {
        "five": [5], "vertex": [5],
    }
    assert get_vor_idx_from_geometry(grid, gdf_to_query=features) == {0: [5], 1: [5]}


def test_neither_source_still_raises(grid):
    with pytest.raises(ValueError, match="shp_to_query or gdf_to_query"):
        get_vor_idx_from_geometry(grid)


def test_a_plain_dataframe_with_a_geometry_column_still_works(grid):
    """It reached the index function and answered (9 cells, any touch), so it is
    coerced rather than broken: a bare pandas ``Series`` has no CRS to reconcile,
    and asking it for one raised ``AttributeError: 'Series' object has no attribute
    'set_crs'`` a couple of frames away from the call."""

    plain = pd.DataFrame({"name": ["five"], "geometry": [CELL_5]})

    assert not isinstance(plain, gpd.GeoDataFrame)
    assert get_vor_idx_from_geometry(grid, gdf_to_query=plain, name_col="name") == {"five": [5]}
    assert get_vor_idx_from_geometry_idx(grid, plain, 0) == [5]


def test_a_query_carrying_no_geometry_is_named(grid):
    """Rather than an AttributeError from inside the reconcile."""

    with pytest.raises(ValueError, match="must carry geometry"):
        get_vor_idx_from_geometry(grid, gdf_to_query=pd.DataFrame({"name": ["five"]}))


@pytest.mark.parametrize(
    "func",
    [
        get_vor_idx_from_geometry,
        get_vor_idx_from_geometry_idx,
        VoronoiGridPlus.get_vor_idx_from_geometry,
        VoronoiGridPlus.get_vor_idx_from_geometry_idx,
    ],
)
def test_every_entry_point_defaults_to_the_measured_rule(func):
    assert inspect.signature(func).parameters["predicate"].default is None


def test_no_entry_point_forces_a_crs():
    for func in (get_vor_idx_from_geometry, VoronoiGridPlus.get_vor_idx_from_geometry):
        assert inspect.signature(func).parameters["crs"].default is None


# ---------------------------------------------------------------- one rule, three spellings


@pytest.mark.parametrize(
    "geometry",
    [CELL_5, ON_VERTEX, ON_FACE, THROUGH_ROW_0, DIAGONAL, ALONG_FACE, OFF_GRID],
    ids=["cell_5", "on_vertex", "on_face", "through_row_0", "diagonal", "along_face", "off_grid"],
)
def test_both_functions_and_both_methods_agree_with_the_shared_selector(grid, geometry):
    """``get_vor_cells_as_series`` is the rule everything else in the library goes
    through; these two answered differently from it, and from each other."""

    frame = _frame(geometry)
    expected = selected_cells(grid.get_vor_cells_as_series(frame, name_field="name"))

    for label, answer in _four_answers(grid, frame).items():
        assert answer == expected, label


def test_the_selector_agrees_on_a_query_in_another_crs(grid):
    """Where the dict function's forced reprojection put it six cells out."""

    in_degrees = _frame(
        gpd.GeoSeries([CELL_5], crs=grid.crs).to_crs("EPSG:4326").iloc[0], crs="EPSG:4326"
    )

    assert selected_cells(grid.get_vor_cells_as_series(in_degrees, name_field="name")) == [5]
    for label, answer in _four_answers(grid, in_degrees).items():
        assert answer == [5], label
