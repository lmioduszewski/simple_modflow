from __future__ import annotations

import numpy as np


def _triangle_coords(triangle) -> np.ndarray:
    vertices = triangle.node
    vert_lookup = {
        int(row["ivert"]): np.array([float(row["x"]), float(row["y"])], dtype=float)
        for row in vertices
    }

    coords = []
    for cell in triangle.ele:
        coords.append(
            np.vstack(
                [
                    vert_lookup[int(cell["iv1"])],
                    vert_lookup[int(cell["iv2"])],
                    vert_lookup[int(cell["iv3"])],
                ]
            )
        )
    return np.asarray(coords, dtype=float)


def _triangle_edge_lengths(coords: np.ndarray) -> np.ndarray:
    a = np.linalg.norm(coords[:, 1] - coords[:, 0], axis=1)
    b = np.linalg.norm(coords[:, 2] - coords[:, 1], axis=1)
    c = np.linalg.norm(coords[:, 0] - coords[:, 2], axis=1)
    return np.vstack([a, b, c]).T


def _triangle_areas(coords: np.ndarray) -> np.ndarray:
    return 0.5 * np.abs(
        (coords[:, 1, 0] - coords[:, 0, 0]) * (coords[:, 2, 1] - coords[:, 0, 1])
        - (coords[:, 2, 0] - coords[:, 0, 0]) * (coords[:, 1, 1] - coords[:, 0, 1])
    )


def _triangle_angles(coords: np.ndarray) -> np.ndarray:
    lengths = _triangle_edge_lengths(coords)
    a = lengths[:, 0]
    b = lengths[:, 1]
    c = lengths[:, 2]

    def angle(opposite: np.ndarray, side1: np.ndarray, side2: np.ndarray) -> np.ndarray:
        denom = 2.0 * side1 * side2
        cos_theta = np.ones_like(denom)
        valid = denom > 0
        cos_theta[valid] = (side1[valid] ** 2 + side2[valid] ** 2 - opposite[valid] ** 2) / denom[valid]
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        return np.degrees(np.arccos(cos_theta))

    return np.vstack(
        [
            angle(c, a, b),
            angle(a, b, c),
            angle(b, c, a),
        ]
    ).T


def _summary_stats(values: np.ndarray, prefix: str) -> dict[str, float]:
    if values.size == 0:
        return {
            f"{prefix}_min": 0.0,
            f"{prefix}_mean": 0.0,
            f"{prefix}_max": 0.0,
            f"{prefix}_std": 0.0,
        }
    return {
        f"{prefix}_min": float(np.min(values)),
        f"{prefix}_mean": float(np.mean(values)),
        f"{prefix}_max": float(np.max(values)),
        f"{prefix}_std": float(np.std(values)),
    }


def _triangle_edge_ratio(lengths: np.ndarray) -> np.ndarray:
    min_lengths = np.min(lengths, axis=1)
    max_lengths = np.max(lengths, axis=1)
    ratios = np.full(len(lengths), np.inf, dtype=float)
    valid = min_lengths > 0
    ratios[valid] = max_lengths[valid] / min_lengths[valid]
    return ratios


def _triangle_neighbor_area_ratios(triangle, areas: np.ndarray) -> np.ndarray:
    edge_to_triangles: dict[tuple[int, int], list[int]] = {}
    for tri_idx, cell in enumerate(triangle.ele):
        verts = [int(cell["iv1"]), int(cell["iv2"]), int(cell["iv3"])]
        for edge in ((verts[0], verts[1]), (verts[1], verts[2]), (verts[2], verts[0])):
            edge_key = tuple(sorted(edge))
            edge_to_triangles.setdefault(edge_key, []).append(tri_idx)

    ratios: list[float] = []
    for tri_indices in edge_to_triangles.values():
        if len(tri_indices) != 2:
            continue
        area_a = float(areas[tri_indices[0]])
        area_b = float(areas[tri_indices[1]])
        smaller = min(area_a, area_b)
        if smaller <= 0:
            continue
        ratios.append(max(area_a, area_b) / smaller)
    return np.asarray(ratios, dtype=float)


def triangle_quality_report(
    triangle,
    *,
    include_voronoi: bool = False,
    validate_voronoi: bool = False,
) -> dict[str, float | int | str]:
    if not hasattr(triangle, "ele") or len(triangle.ele) == 0:
        raise ValueError("triangle mesh must be built before quality_report() can be computed")

    coords = _triangle_coords(triangle)
    areas = _triangle_areas(coords)
    lengths = _triangle_edge_lengths(coords)
    angles = _triangle_angles(coords)
    edge_ratios = _triangle_edge_ratio(lengths)
    neighbor_area_ratios = _triangle_neighbor_area_ratios(triangle, areas)

    node_markers = triangle.node["boundary_marker"]
    boundary_vertices = int(np.count_nonzero(node_markers))
    interior_vertices = int(len(node_markers) - boundary_vertices)
    tiny_area_threshold = max(float(np.mean(areas)) * 1e-6, 1e-12)

    report: dict[str, float | int | str] = {
        "num_triangles": int(len(coords)),
        "num_vertices": int(len(triangle.node)),
        "num_boundary_vertices": boundary_vertices,
        "num_interior_vertices": interior_vertices,
        "duplicate_vertex_count": int(
            len(triangle.node) - len(np.unique(np.c_[triangle.node["x"], triangle.node["y"]], axis=0))
        ),
        "zero_area_triangle_count": int(np.count_nonzero(areas <= 0.0)),
        "tiny_triangle_count": int(np.count_nonzero(areas <= tiny_area_threshold)),
        "sliver_triangle_count": int(np.count_nonzero(np.min(angles, axis=1) < 20.0)),
    }
    report.update(_summary_stats(areas, "triangle_area"))
    report.update(_summary_stats(lengths.reshape(-1), "edge_length"))
    report.update(_summary_stats(angles.reshape(-1), "triangle_angle"))
    report.update(_summary_stats(edge_ratios, "triangle_edge_ratio"))
    report.update(_summary_stats(neighbor_area_ratios, "neighbor_area_ratio"))
    report["triangle_angle_min_overall"] = float(np.min(angles))
    report["triangle_angle_max_overall"] = float(np.max(angles))

    if include_voronoi or validate_voronoi:
        from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus

        try:
            vor = VoronoiGridPlus(triangle)
        except Exception as exc:
            report["voronoi_status"] = "error"
            report["voronoi_error"] = str(exc)
        else:
            report["voronoi_status"] = "built"
            cell_areas = np.asarray(vor.get_cell_areas(), dtype=float)
            report["voronoi_cell_count"] = int(len(cell_areas))
            report["voronoi_invalid_cell_count"] = int(np.count_nonzero(cell_areas <= 0.0))
            if include_voronoi:
                report.update(_summary_stats(cell_areas, "voronoi_area"))

    return report
