"""Regression tests for the cached CRS row pointer ``VoronoiGridPlus.ia``.

``find_adjacent_cells`` and ``LAKBuilder._adjacent_with_metrics`` used to recompute
``sum(iac[:cell])`` on every call -- an O(N) scan per cell, i.e. O(N^2) when sweeping
a whole lake footprint. On a large lake (e.g. a mine pit spanning thousands of cells)
that turned an automatic LAK build into a multi-minute hang; ``lakes.py`` made it worse
by using the Python builtin ``sum()`` over a numpy slice.

The per-cell offset is now the cached MODFLOW row pointer ``ia = [0, cumsum(iac)]``
indexed in O(1). These tests pin ``ia``'s correctness and prove the connectivity
lookups are byte-for-byte unchanged from the old running-sum form.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

from myflopy import LAKBuilder, ModelContext
from myflopy.modflow.mf6.grid.geometry import find_adjacent_cells
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _grid() -> VoronoiGridPlus:
    """A deterministic 3x2 rectangular Voronoi grid (6 cells, 2 layers)."""
    vertices = np.array(
        [
            [0, 0], [1, 0], [2, 0], [3, 0],
            [0, 1], [1, 1], [2, 1], [3, 1],
            [0, 2], [1, 2], [2, 2], [3, 2],
        ],
        dtype=float,
    )
    iverts = [
        [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6],
        [4, 5, 9, 8], [5, 6, 10, 9], [6, 7, 11, 10],
    ]
    centers = np.array(
        [[0.5, 0.5], [1.5, 0.5], [2.5, 0.5], [0.5, 1.5], [1.5, 1.5], [2.5, 1.5]]
    )
    grid = VoronoiGridPlus(verts=vertices, iverts=iverts, xcyc=centers)
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [12.0] * 6, 1: [8.0] * 6, 2: [0.0] * 6},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    grid.get_disu_connectivity(validate=True)
    return grid


def test_ia_equals_cumulative_iac():
    vor = _grid()
    iac = np.asarray(vor.iac)
    ia = vor.ia

    assert len(ia) == len(iac) + 1
    assert int(ia[0]) == 0
    assert int(ia[-1]) == int(vor.nja) == int(iac.sum())
    # O(1) ia[c] must equal the old O(N) running sum for every cell.
    for c in range(len(iac)):
        assert int(ia[c]) == int(np.sum(iac[:c]))


def test_find_adjacent_cells_offsets_unchanged():
    """The ia-based start reproduces exactly the old sum(iac[:c]) slicing."""
    vor = _grid()
    iac = np.asarray(vor.iac)
    ja = np.asarray(vor.ja)

    for c in range(len(iac)):
        start = int(np.sum(iac[:c]))  # the retired O(N) formula
        expected = [int(x) for x in ja[start:start + int(iac[c])] if int(x) != c]
        assert find_adjacent_cells(vor, c) == expected

    # Adjacency must be symmetric (c is a neighbor of each of its neighbors).
    for c in range(len(iac)):
        for neighbor in find_adjacent_cells(vor, c):
            assert c in find_adjacent_cells(vor, neighbor)


def test_ia_cache_invalidates_when_connectivity_recomputes():
    vor = _grid()
    first = vor.ia.copy()
    assert vor._ia is not None

    vor.get_disu_connectivity(validate=True)  # rebuilds iac -> ia must be dropped
    assert vor._ia is None

    np.testing.assert_array_equal(vor.ia, first)  # lazily recomputes identically


def test_lak_adjacent_metrics_read_cl12_hwva_at_row_pointer():
    """LAKBuilder connection metrics index cl12/hwva at ia[cell] + idx + 1."""
    vor = _grid()
    # A polygon wholly inside the centre cell (centre (1.5, 1.5) -> cell 4).
    lakes = gpd.GeoDataFrame(
        {"name": ["a"]},
        geometry=[Polygon([(1.05, 1.05), (1.95, 1.05), (1.95, 1.95), (1.05, 1.95)])],
        crs=vor.crs,
    )
    builder = LAKBuilder(
        context=ModelContext(grid=vor, domain=np.ones((2, 6), dtype=int)),
        nper=1,
        lakes=lakes,
        lake_id_field="name",
        starting_stage={"a": 11.0},
        lake_bottom={"a": 8.0},
        bed_leakance=0.1,
    )

    cell = builder.lake_cells["a"][0]
    metrics = builder._adjacent_with_metrics(cell)
    neighbors = find_adjacent_cells(vor, cell)
    assert [m[0] for m in metrics] == neighbors

    start = int(vor.ia[cell])
    cl12, hwva = np.asarray(vor.cl12), np.asarray(vor.hwva)
    for idx, (_other, length, width) in enumerate(metrics):
        assert length == float(cl12[start + idx + 1])
        assert width == float(hwva[start + idx + 1])
