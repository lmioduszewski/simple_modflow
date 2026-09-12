"""Edge cells are the ACTIVE cells bordering the outside, or an inactive cell.

Given an idomain, this outlined the active area and asked which cells its ring
touched. The inactive set came from a touch join, so a face-aligned strip marked
a whole extra column inactive and the answer lost the true edge column while
gaining four interior cells (14 against 16); a ring runs through vertices, so
corner-only neighbours were kept; and outer and hole rings were concatenated
without de-duplication, so the cells around a one-cell hole came back 36 times
for 27 cells -- which the simple model's drain mode wrote as 36 DRN records.

Every expectation here is worked out from row/column arithmetic on a regular
grid, not from the code under test.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

from myflopy.modflow.mf6.grid.selection import get_grid_edge_cells
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

N = 6                 # 6 x 6 cells
D = 25.0              # each 25 ft square, so the grid spans 0..150


def _grid(order: list[int] | None = None, idomain_path=None) -> VoronoiGridPlus:
    """N x N square cells, row-major from the bottom left; ``order`` permutes the ids."""

    verts = np.array([[i * D, j * D] for j in range(N + 1) for i in range(N + 1)], dtype=float)
    iverts = [
        [j * (N + 1) + i, (j + 1) * (N + 1) + i, (j + 1) * (N + 1) + i + 1, j * (N + 1) + i + 1]
        for j in range(N)
        for i in range(N)
    ]
    xcyc = [[(i + 0.5) * D, (j + 0.5) * D] for j in range(N) for i in range(N)]
    if order is not None:
        iverts = [iverts[i] for i in order]
        xcyc = [xcyc[i] for i in order]
    kwargs = {} if idomain_path is None else {"idomain_path": idomain_path}
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(xcyc), **kwargs)


def _id(row: int, col: int) -> int:
    return row * N + col


def _perimeter() -> set[int]:
    """The cells on the mesh's own boundary: 20 of the 36."""

    return {_id(r, c) for r in range(N) for c in range(N) if r in (0, N - 1) or c in (0, N - 1)}


def _column_block_edge(last_col: int) -> set[int]:
    """The edge of the active block of columns 0..last_col."""

    return {
        _id(r, c)
        for r in range(N)
        for c in range(last_col + 1)
        if r in (0, N - 1) or c in (0, last_col)
    }


def _neighbours(cell: int) -> set[int]:
    row, col = divmod(cell, N)
    return {
        _id(r, c)
        for r, c in ((row - 1, col), (row + 1, col), (row, col - 1), (row, col + 1))
        if 0 <= r < N and 0 <= c < N
    }


@pytest.fixture
def grid():
    return _grid()


def _strip_file(tmp_path, geometry, crs="EPSG:2927", name="idomain.gpkg"):
    path = tmp_path / name
    gpd.GeoDataFrame({"zone": ["inactive"]}, geometry=[geometry], crs=crs).to_file(path, driver="GPKG")
    return path


def test_with_no_idomain_the_edge_is_the_mesh_boundary(grid):
    assert get_grid_edge_cells(grid) == sorted(_perimeter())


def test_a_face_aligned_inactive_strip_keeps_the_true_edge_column(grid, tmp_path):
    """The defect: 14 cells, missing column 3 entirely and adding four interior ones."""

    inactive = [_id(r, c) for r in range(N) for c in (4, 5)]
    expected = sorted(_column_block_edge(3))

    assert get_grid_edge_cells(grid, idomain=inactive) == expected
    assert len(expected) == 16
    path = _strip_file(tmp_path, box(100.0, 0.0, 150.0, 150.0))
    assert get_grid_edge_cells(grid, idomain_path=path) == expected


def test_a_corner_cell_does_not_drag_in_its_diagonal(grid):
    """A ring passes through the corner vertex, so cell 28 came along."""

    edge = get_grid_edge_cells(grid, idomain=[35])

    assert edge == sorted(_perimeter() - {35})
    assert _id(4, 4) not in edge


def test_a_one_cell_hole_rings_it_once(grid):
    """36 entries for 27 cells before, of which four were only vertex contacts."""

    hole = _id(2, 2)

    edge = get_grid_edge_cells(grid, idomain=[hole])

    assert edge == sorted(_perimeter() | _neighbours(hole))
    assert len(edge) == len(set(edge)) == 24
    assert not {_id(1, 1), _id(1, 3), _id(3, 1), _id(3, 3)} & set(edge)


def test_a_two_by_two_hole_rings_it_once(grid):
    hole = {_id(2, 2), _id(2, 3), _id(3, 2), _id(3, 3)}
    ring = {cell for member in hole for cell in _neighbours(member)} - hole

    edge = get_grid_edge_cells(grid, idomain=sorted(hole))

    assert edge == sorted(_perimeter() | ring)
    assert len(edge) == 28


def test_include_interiors_false_drops_the_cells_around_a_hole(grid):
    hole = _id(2, 2)

    assert get_grid_edge_cells(grid, idomain=[hole], include_interiors=False) == sorted(_perimeter())


def test_an_island_inside_a_ring_keeps_its_own_edge(grid):
    """include_interiors=False drops a hole's ring, not a separate component's."""

    ring = {cell for cell in range(N * N) if cell not in _perimeter() and cell != _id(2, 2)}

    edge = get_grid_edge_cells(grid, idomain=sorted(ring), include_interiors=False)

    assert _id(2, 2) in edge


def test_an_explicit_idomain_wins_over_the_grids_own_file(tmp_path):
    """The file was consulted first, so idomain= was ignored and idomain=[] fell
    through to it."""

    path = _strip_file(tmp_path, box(100.0, 0.0, 150.0, 150.0))
    grid = _grid(idomain_path=path)

    assert get_grid_edge_cells(grid, idomain=[35]) == sorted(_perimeter() - {35})
    assert get_grid_edge_cells(grid, idomain=[]) == sorted(_perimeter())
    assert get_grid_edge_cells(grid) == sorted(_column_block_edge(3))


def test_both_arguments_are_refused(grid, tmp_path):
    path = _strip_file(tmp_path, box(100.0, 0.0, 150.0, 150.0))

    with pytest.raises(ValueError, match="not both"):
        get_grid_edge_cells(grid, idomain=[35], idomain_path=path)


def test_an_array_of_ids_is_accepted(grid):
    """np.array([34, 35]) raised 'truth value of an array is ambiguous'."""

    # 28 sits under the inactive 34, so it joins the edge; 29 and 33 already were on it.
    expected = sorted((_perimeter() - {34, 35}) | {_id(4, 4)})

    assert get_grid_edge_cells(grid, idomain=np.array([34, 35])) == expected


@pytest.mark.parametrize(
    "mask",
    [
        [True] * (N * N),
        [1] * (N * N - 2) + [0, 0],
    ],
)
def test_a_mask_is_refused_naming_the_conversion(grid, mask):
    """idomain= means INACTIVE ids; a 0/1 array reads as ids 0 and 1."""

    with pytest.raises(ValueError, match="flatnonzero"):
        get_grid_edge_cells(grid, idomain=mask)


def test_an_id_the_grid_does_not_have_is_refused(grid):
    with pytest.raises(ValueError, match=r"\[71\]"):
        get_grid_edge_cells(grid, idomain=[35, 71])


def test_a_polygon_that_only_grazes_a_face_marks_nothing(grid, tmp_path):
    """0.005 ft of overlap is noise, not an inactive column (ledger 210)."""

    path = _strip_file(tmp_path, box(100.0 - 0.005, 0.0, 150.0, 150.0))

    assert get_grid_edge_cells(grid, idomain_path=path) == sorted(_column_block_edge(3))


def test_a_polygon_covering_part_of_a_column_marks_it(grid, tmp_path):
    """40% of a cell is coverage, so column 3 goes inactive and the edge moves in."""

    path = _strip_file(tmp_path, box(90.0, 0.0, 150.0, 150.0))

    assert get_grid_edge_cells(grid, idomain_path=path) == sorted(_column_block_edge(2))


def test_an_idomain_file_in_another_crs_is_reprojected(grid, tmp_path):
    """union_all() dropped the file's CRS, and a 4326 file then raised KeyError: 0."""

    wgs84 = gpd.GeoSeries([box(100.0, 0.0, 150.0, 150.0)], crs=grid.crs).to_crs("EPSG:4326")
    path = _strip_file(tmp_path, wgs84.iloc[0], crs="EPSG:4326")

    assert get_grid_edge_cells(grid, idomain_path=path) == sorted(_column_block_edge(3))


def test_an_idomain_file_that_covers_nothing_raises(grid, tmp_path):
    path = _strip_file(tmp_path, box(5_000.0, 5_000.0, 6_000.0, 6_000.0), name="away.gpkg")

    with pytest.raises(ValueError, match="away.gpkg"):
        get_grid_edge_cells(grid, idomain_path=path)


def test_the_answer_is_ascending_whatever_order_the_grid_is_numbered_in():
    """The ring branch returned per-ring order, so a multi-ring grid came back
    unsorted -- and the simple model zips this list against per-cell heads."""

    order = list(reversed(range(N * N)))
    scrambled = _grid(order=order)
    hole = order.index(_id(2, 2))
    expected = sorted(order.index(cell) for cell in _perimeter() | _neighbours(_id(2, 2)))

    edge = get_grid_edge_cells(scrambled, idomain=[hole])

    assert edge == expected == sorted(edge)


def test_the_simple_models_perimeter_cells_are_each_listed_once(grid):
    """Its drain mode zips this with conductances, so a duplicate was a second DRN
    record on the same cell: 36 records over 27 cells, 3600 against 2400."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.simplemodel import _resolve_boundary_cells

    cells = _resolve_boundary_cells(SimpleNamespace(boundary_cells=None, vor=grid))

    assert cells == sorted(set(cells)) == sorted(_perimeter())


def test_two_active_components_still_come_back_ascending(grid):
    """Column 2 inactive splits the grid in two. The components are found lowest
    id first, but their members interleave (0, 1, 6, 7, ... against 3, 4, 5, ...),
    so concatenating them is not sorted."""

    split = [_id(r, 2) for r in range(N)]

    edge = get_grid_edge_cells(grid, idomain=split, include_interiors=False)

    left = {_id(r, c) for r in range(N) for c in (0, 1)}            # 2 wide: all of it
    right = {
        _id(r, c) for r in range(N) for c in (3, 4, 5) if c in (3, 5) or r in (0, N - 1)
    }

    assert edge == sorted(edge)
    assert edge == sorted(left | right)


def test_a_short_boolean_list_is_refused_too(grid):
    """A full-length mask is caught as an MF6 array; a shorter one would otherwise
    read as the ids 0 and 1."""

    with pytest.raises(ValueError, match="mask"):
        get_grid_edge_cells(grid, idomain=[True, False, True])
