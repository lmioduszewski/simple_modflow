"""A barrier cuts CONNECTIONS: the ones whose centre-to-centre segment crosses it.

``barrier_faces`` used to take the faces a trace crosses. That is right only for a
trace drawn exactly along faces. A trace through cells crosses the faces BETWEEN
consecutive cells along the fault, which block nothing across it: MF6 measured
68-100% of the cross-fault flow still crossing on every grid tried -- square,
irregular Voronoi, refined along the fault, the canonical grid -- while the HFB
summary read as sealed, because it sums |q| on the barrier's own faces.

These tests need no MF6 run: a barrier seals exactly when no cell on one side can
reach a cell on the other without crossing a cut connection.
"""

from __future__ import annotations

from collections import deque

import numpy as np
import pytest
import shapely as shp
from shapely.geometry import LineString, MultiLineString, box

from myflopy.modflow.mf6.grid.barriers import barrier_faces
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _square(nx: int = 40, ny: int = 20, d: float = 50.0) -> VoronoiGridPlus:
    """nx x ny square cells of side d, numbered row-major from the bottom left."""

    verts = np.array([[i * d, j * d] for j in range(ny + 1) for i in range(nx + 1)], dtype=float)
    iverts = [
        [j * (nx + 1) + i, (j + 1) * (nx + 1) + i, (j + 1) * (nx + 1) + i + 1, j * (nx + 1) + i + 1]
        for j in range(ny)
        for i in range(nx)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(ny) for i in range(nx)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _irregular(on_the_trace: float | None = 1000.0) -> VoronoiGridPlus:
    """A Voronoi grid of 400 random generators over 2000 x 1000, some ON ``x=on_the_trace``."""

    rng = np.random.default_rng(7)
    points = rng.uniform([10.0, 10.0], [1990.0, 990.0], size=(400, 2))
    if on_the_trace is not None:
        points[:12, 0] = on_the_trace
        points[:12, 1] = np.linspace(40.0, 960.0, 12)
    frame = box(0.0, 0.0, 2000.0, 1000.0)
    regions = shp.voronoi_polygons(shp.MultiPoint(points), extend_to=frame, ordered=True)
    cells = [shp.orient_polygons(region.intersection(frame), exterior_cw=True) for region in regions.geoms]
    index: dict[tuple[float, float], int] = {}
    iverts = []
    for cell in cells:
        ring = []
        for x, y in list(cell.exterior.coords)[:-1]:
            key = (round(x, 6), round(y, 6))
            ring.append(index.setdefault(key, len(index)))
        iverts.append(ring)
    verts = np.array(list(index), dtype=float)
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=points)


def _centres(vor):
    return np.asarray(vor.points, dtype=float)


def _sealed(vor, faces, left, right) -> bool:
    """Whether no cell in ``left`` reaches a cell in ``right`` without crossing a cut."""

    blocked = {tuple(pair) for pair in faces}
    adjacency = vor.adjacent_cells_idx
    seen = set(left)
    queue = deque(left)
    while queue:
        cell = queue.popleft()
        for neighbour in adjacency[cell]:
            neighbour = int(neighbour)
            if neighbour in seen or (min(cell, neighbour), max(cell, neighbour)) in blocked:
                continue
            seen.add(neighbour)
            queue.append(neighbour)
    return not (seen & set(right))


def _sides(vor, trace, *, margin: float = 1.0):
    """Cells clearly on either side of a straight trace (off it by more than ``margin``)."""

    (x0, y0), (x1, y1) = trace.coords[0], trace.coords[-1]
    c = _centres(vor)
    cross = (x1 - x0) * (c[:, 1] - y0) - (y1 - y0) * (c[:, 0] - x0)
    scale = np.hypot(x1 - x0, y1 - y0)
    return (
        [int(i) for i in np.flatnonzero(cross > margin * scale)],
        [int(i) for i in np.flatnonzero(cross < -margin * scale)],
    )


TRACES = {
    "along-faces": LineString([(1000.0, -10.0), (1000.0, 1010.0)]),
    "through-centres": LineString([(1025.0, -10.0), (1025.0, 1010.0)]),
    "off-centre": LineString([(1010.0, -10.0), (1010.0, 1010.0)]),
    "45-degrees": LineString([(450.0, -50.0), (1550.0, 1050.0)]),
    "68-degrees": LineString([(800.0, -10.0), (1212.0, 1010.0)]),
}


@pytest.mark.parametrize("name", list(TRACES))
def test_a_trace_across_a_square_grid_seals_whatever_it_passes_through(name):
    vor, trace = _square(), TRACES[name]
    left, right = _sides(vor, trace)
    assert _sealed(vor, barrier_faces(vor, trace), left, right)


def test_a_trace_along_faces_still_resolves_to_exactly_those_faces():
    """The one case the old rule got right is unchanged: 20 faces at x=1000."""

    faces = barrier_faces(_square(), TRACES["along-faces"])
    assert len(faces) == 20
    assert all(b - a == 1 and a % 40 == 19 for a, b in faces)


@pytest.mark.parametrize("on_the_trace", [1000.0, None], ids=["generators-on-the-trace", "none-on-it"])
def test_a_trace_across_an_irregular_voronoi_grid_seals(on_the_trace):
    vor = _irregular(on_the_trace)
    trace = LineString([(1000.0, -10.0), (1000.0, 1010.0)])
    left, right = _sides(vor, trace, margin=0.0)
    assert _sealed(vor, barrier_faces(vor, trace), left, right)


@pytest.mark.parametrize(
    ("trace", "cut_between_columns"),
    [
        (LineString([(37.5, -1.0), (37.5, 101.0)]), (1, 2)),  # drawn north: left is west
        (LineString([(37.5, 101.0), (37.5, -1.0)]), (0, 1)),  # drawn south: left is east
    ],
    ids=["drawn-north", "drawn-south"],
)
def test_a_centre_on_the_trace_belongs_to_its_left(trace, cut_between_columns):
    """Column 1's centres lie on x=37.5. They go LEFT, looking along the trace."""

    west, east = cut_between_columns
    faces = barrier_faces(_square(nx=4, ny=4, d=25.0), trace)
    assert faces == [(4 * row + west, 4 * row + east) for row in range(4)]


def test_a_closed_ring_is_watertight():
    """It used to leave every cell connected: a ring passes THROUGH cells."""

    vor = _square()
    ring = box(512.0, 237.0, 1488.0, 763.0).exterior
    c = _centres(vor)
    inside = shp.contains_xy(box(512.0, 237.0, 1488.0, 763.0), c[:, 0], c[:, 1])
    assert _sealed(
        vor, barrier_faces(vor, ring),
        [int(i) for i in np.flatnonzero(inside)], [int(i) for i in np.flatnonzero(~inside)],
    )


def test_an_open_ended_trace_does_not_seal_and_cuts_nothing_past_its_end():
    vor = _square()
    trace = LineString([(1010.0, -10.0), (1010.0, 600.0)])
    faces = barrier_faces(vor, trace)
    left, right = _sides(vor, LineString([(1010.0, -10.0), (1010.0, 1010.0)]))
    assert not _sealed(vor, faces, left, right)
    c = _centres(vor)
    assert all(max(c[a, 1], c[b, 1]) < 650.0 for a, b in faces)


def test_every_cut_connection_crosses_the_trace():
    vor = _irregular()
    trace = MultiLineString([[(300.0, -10.0), (700.0, 1010.0)], [(1500.0, -10.0), (1300.0, 1010.0)]])
    c = _centres(vor)
    faces = barrier_faces(vor, trace)
    assert faces
    for a, b in faces:
        assert LineString([c[a], c[b]]).distance(trace) <= 1e-9


def test_a_trace_that_starts_at_a_cell_centre_cuts_nothing_there():
    """A 0.01 stub from a centre separates nothing. Given a side like a centre ON a
    trace, its start cut the connection to the neighbour on its right."""

    vor = _square(nx=4, ny=4, d=25.0)
    x, y = _centres(vor)[5]
    assert barrier_faces(vor, LineString([(x, y), (x + 0.01, y)])) == []
