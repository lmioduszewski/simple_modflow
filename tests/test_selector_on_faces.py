"""The shared selector picks the cells a geometry covers, not every cell it touches.

``get_vor_cells_as_series`` joined with ``intersects`` and no measure, so a
polygon that was exactly cell 5 selected 9 cells, a face-aligned block of 4
selected all 16, and a point on a vertex selected 4. Everything that asks "which
cells does this shape pick?" goes through it: regions and ``map(select=...)``,
particle release, cross-sections, top/bottom reconciliation, and the legacy
boundary classes. The default is now the measured rule; an explicit predicate
still asks for the raw topological join.
"""

from __future__ import annotations

import inspect
import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, box

from myflopy.modflow.mf6.grid.selection import get_vor_cells_as_dict, get_vor_cells_as_series
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.regions import RegionRegistry

LOGGER = "myflopy.modflow.mf6.grid.selection"
CELL_5 = box(25.0, 25.0, 50.0, 50.0)


def _grid(n: int = 4, length: float = 100.0) -> VoronoiGridPlus:
    """n x n square cells, row-major from the bottom left."""

    d = length / n
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


@pytest.fixture
def grid():
    return _grid()


def _cells(grid, geometry, **kwargs) -> list[int]:
    series = get_vor_cells_as_series(grid.gdf_vorPolys, geometry, **kwargs)
    return sorted(cell for group in series.tolist() for cell in group)


@pytest.mark.parametrize(
    "geometry, expected",
    [
        (CELL_5, [5]),                                   # was 9: every neighbour it touched
        (box(25.0, 25.0, 75.0, 75.0), [5, 6, 9, 10]),    # was all 16
        (box(25.0, 25.0, 50.0, 75.0), [5, 9]),
        (Point(37.5, 37.5), [5]),                        # inside
        (Point(25.0, 37.5), [4]),                        # on a face: was [4, 5]
        (Point(25.0, 25.0), [0]),                        # on a vertex: was [0, 1, 4, 5]
        (LineString([(0.0, 37.5), (100.0, 37.5)]), [4, 5, 6, 7]),
        (LineString([(0.0, 0.0), (100.0, 100.0)]), [0, 5, 10, 15]),   # was 10 cells
    ],
)
def test_a_geometry_reaches_the_cells_it_covers(grid, geometry, expected):
    assert _cells(grid, geometry) == expected


def test_a_line_along_a_face_reaches_both_rows(grid):
    """Ledger 212: an unsplit field reaches both sides, and the length is split."""

    assert _cells(grid, LineString([(0.0, 25.0), (100.0, 25.0)])) == [0, 1, 2, 3, 4, 5, 6, 7]


def test_a_cell_polygon_survives_a_crs_round_trip(grid):
    """The near miss is commoner than exact coincidence: 4 cells, three of them
    slivers, before the rule."""

    wandered = (
        gpd.GeoSeries([CELL_5], crs=grid.crs).to_crs("EPSG:4326").to_crs(grid.crs).iloc[0]
    )

    assert _cells(grid, wandered) == [5]


def test_an_explicit_predicate_still_asks_for_the_raw_join(grid):
    """The escape hatch: topology, for a caller that wants it."""

    assert len(_cells(grid, CELL_5, predicate="intersects")) == 9


def test_a_feature_that_reaches_nothing_is_an_empty_series_and_a_warning(grid, caplog):
    """It used to come back as an empty GeoDataFrame, which callers indexed as a
    Series one line later."""

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        series = get_vor_cells_as_series(grid.gdf_vorPolys, box(500.0, 500.0, 600.0, 600.0))

    assert isinstance(series, pd.Series)
    assert series.empty and series.name == "cells"
    assert any("reach no grid cell" in record.getMessage() for record in caplog.records)


def test_a_partial_miss_keeps_the_rest_and_names_the_missed(grid, caplog):
    features = gpd.GeoDataFrame(
        {"ExploName": ["inside", "away"]},
        geometry=[CELL_5, box(500.0, 500.0, 600.0, 600.0)],
        crs=grid.crs,
    )

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        series = get_vor_cells_as_series(grid.gdf_vorPolys, features)

    assert series.index.tolist() == ["inside"]
    assert any("'away'" in record.getMessage() for record in caplog.records)


@pytest.mark.parametrize(
    "func",
    [
        get_vor_cells_as_series,
        get_vor_cells_as_dict,
        VoronoiGridPlus.get_vor_cells_as_series,
        VoronoiGridPlus.get_vor_cells_as_dict,
        RegionRegistry.add_from_geometry,
        SimulationBase.add_region_from_geometry,
    ],
)
def test_every_entry_point_defaults_to_the_measured_rule(func):
    """The default was restated at every forwarder, so fixing the free function
    alone reached almost nothing."""

    assert inspect.signature(func).parameters["predicate"].default is None


def test_the_grid_wrapper_gives_what_the_free_function_gives(grid):
    assert grid.get_vor_cells_as_series(CELL_5).tolist() == [[5]]


def test_the_dict_variant_finds_a_point_on_a_face(grid, tmp_path):
    """It asked for `contains`, which a point on a face satisfies for NO cell, and
    then said the point did not intersect the grid at all."""

    path = tmp_path / "locs.gpkg"
    gpd.GeoDataFrame(
        {"name": ["middle", "on_face"]},
        geometry=[Point(37.5, 37.5), Point(25.0, 37.5)],
        crs=grid.crs,
    ).to_file(path, driver="GPKG")

    found = get_vor_cells_as_dict(grid, path, loc_name_field="name")

    assert found == {"middle": [[5]], "on_face": [[4]]}


def test_the_layer_guard_still_refuses_an_in_memory_geometry(grid):
    with pytest.raises(ValueError, match="selects a layer inside a vector FILE"):
        get_vor_cells_as_series(grid.gdf_vorPolys, CELL_5, layer=0)


def test_no_caller_asks_the_selector_for_contains():
    """``contains`` is false for EVERY cell when a point sits on a shared face, so
    the two callers that asked for it dropped such a point and then reported it as
    outside the grid altogether. A ratchet, because the argument reads harmless."""

    src = Path(__file__).resolve().parents[1] / "src"
    offenders = [
        f"{path.relative_to(src)}:{number}"
        for path in sorted(src.rglob("*.py"))
        for number, line in enumerate(path.read_text().splitlines(), start=1)
        if "predicate='contains'" in line or 'predicate="contains"' in line
    ]

    assert offenders == []
