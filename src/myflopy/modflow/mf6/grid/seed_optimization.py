from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from shapely.geometry import Point, Polygon


@dataclass(slots=True)
class RelaxationResult:
    points: list[tuple[float, float]]
    iterations_run: int
    max_move: float
    free_vertex_count: int
    moved_vertex_count: int


def _protected_vertex_ids(
    node,
    protected_geometries: list[Polygon] | None,
    tolerance: float,
) -> set[int]:
    protected: set[int] = set()
    if not protected_geometries:
        return protected

    buffered = [geom.boundary.buffer(tolerance) for geom in protected_geometries if not geom.is_empty]
    if not buffered:
        return protected

    for row in node:
        point = Point(float(row["x"]), float(row["y"]))
        if any(geom.covers(point) for geom in buffered):
            protected.add(int(row["ivert"]))
    return protected


def get_fixed_and_free_vertex_ids(
    *,
    node,
    protected_geometries: list[Polygon] | None = None,
    tolerance: float = 1e-6,
) -> tuple[set[int], list[int]]:
    fixed = {int(row["ivert"]) for row in node if int(row["boundary_marker"]) != 0}
    fixed.update(_protected_vertex_ids(node, protected_geometries, tolerance))
    free_ids = sorted(int(row["ivert"]) for row in node if int(row["ivert"]) not in fixed)
    return fixed, free_ids


def _coord_key(x: float, y: float, digits: int = 8) -> tuple[float, float]:
    return (round(float(x), digits), round(float(y), digits))


def map_voronoi_cells_to_vertex_ids(
    *,
    node,
    vor_points,
    digits: int = 8,
) -> dict[int, int]:
    point_lookup = {
        _coord_key(point[0], point[1], digits=digits): idx
        for idx, point in enumerate(np.asarray(vor_points, dtype=float))
    }
    mapping: dict[int, int] = {}
    for row in node:
        key = _coord_key(row["x"], row["y"], digits=digits)
        if key in point_lookup:
            mapping[int(row["ivert"])] = point_lookup[key]
    return mapping


def cvt_relaxed_full_seed_points(
    *,
    node,
    vor_points,
    vor_polygons,
    domain: Polygon,
    damping: float = 0.35,
    min_move: float = 1e-3,
    protected_geometries: list[Polygon] | None = None,
    tolerance: float = 1e-6,
) -> RelaxationResult:
    damping = float(np.clip(damping, 0.0, 1.0))
    coords = {
        int(row["ivert"]): np.array([float(row["x"]), float(row["y"])], dtype=float)
        for row in node
    }
    _, free_ids = get_fixed_and_free_vertex_ids(
        node=node,
        protected_geometries=protected_geometries,
        tolerance=tolerance,
    )
    if not free_ids:
        return RelaxationResult(
            points=[],
            iterations_run=0,
            max_move=0.0,
            free_vertex_count=0,
            moved_vertex_count=0,
        )

    mapping = map_voronoi_cells_to_vertex_ids(node=node, vor_points=vor_points)
    proposed: dict[int, np.ndarray] = {}
    max_move = 0.0

    for vertex in free_ids:
        polygon_idx = mapping.get(vertex)
        if polygon_idx is None or polygon_idx >= len(vor_polygons):
            continue

        polygon = vor_polygons[polygon_idx]
        if polygon.is_empty or polygon.area <= 0:
            continue

        current = coords[vertex]
        centroid = polygon.centroid
        target = np.array([float(centroid.x), float(centroid.y)], dtype=float)
        candidate = current + damping * (target - current)

        if not domain.covers(Point(float(candidate[0]), float(candidate[1]))):
            continue

        move = float(np.linalg.norm(candidate - current))
        if move < min_move:
            continue

        proposed[vertex] = candidate
        max_move = max(max_move, move)

    coords.update(proposed)
    full_points = [tuple(float(value) for value in coords[vertex]) for vertex in free_ids]

    return RelaxationResult(
        points=full_points,
        iterations_run=1 if proposed else 0,
        max_move=max_move,
        free_vertex_count=len(free_ids),
        moved_vertex_count=len(proposed),
    )
