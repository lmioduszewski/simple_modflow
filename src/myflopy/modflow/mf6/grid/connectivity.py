from __future__ import annotations

import numpy as np
from shapely.geometry import LineString, MultiLineString
from shapely.prepared import prep


def build_disu_connectivity(
    gdf_vor_polys,
    adjacent_cells_idx,
    *,
    use_representative_point: bool = False,
    tol: float = 0.0,
    validate: bool = True,
):
    """
    Build DISU connectivity vectors (iac, ja, cl12, hwva, nja) from planar polygons.
    """
    geoms = gdf_vor_polys.geometry.values

    if use_representative_point:
        centers_seq = gdf_vor_polys.representative_point()
    else:
        centers_seq = gdf_vor_polys.centroid
    centers = np.array([(point.x, point.y) for point in centers_seq], dtype=np.float64)

    iac, ja, cl12, hwva = [], [], [], []

    print('getting connectivity properties (iac, ja, cl12, hwva, nja)')
    for i, nbrs in enumerate(adjacent_cells_idx):
        poly_i = geoms[i]
        prep_i = prep(poly_i)

        edge_neighbors = []
        faces = []
        for j in nbrs:
            poly_j = geoms[j]
            if not prep_i.intersects(poly_j):
                continue

            ij_int = poly_i.boundary.intersection(poly_j.boundary)

            if isinstance(ij_int, LineString):
                length = ij_int.length
                segment = ij_int if length > 0.0 else None
            elif isinstance(ij_int, MultiLineString) and len(ij_int.geoms) > 0:
                length = sum(seg.length for seg in ij_int.geoms)
                segment = max(ij_int.geoms, key=lambda seg: seg.length)
                if segment.length == 0.0:
                    segment = None
            else:
                length, segment = 0.0, None

            if segment is not None and length > tol:
                edge_neighbors.append(j)
                faces.append((float(length), segment))

        if edge_neighbors:
            order = np.argsort(edge_neighbors)
            edge_neighbors = [edge_neighbors[k] for k in order]
            faces = [faces[k] for k in order]

        iac.append(1 + len(edge_neighbors))
        ja.append(i)
        cl12.append(0.0)
        hwva.append(0.0)

        ci = centers[i]
        for j, (length, segment) in zip(edge_neighbors, faces, strict=False):
            hwva.append(length)

            (x1, y1) = segment.coords[0]
            (x2, y2) = segment.coords[-1]
            ex, ey = (x2 - x1), (y2 - y1)
            edge_length = np.hypot(ex, ey)

            dij = centers[j] - ci
            if edge_length == 0.0:
                spacing = float(np.hypot(dij[0], dij[1]))
            else:
                nx, ny = -ey / edge_length, ex / edge_length
                spacing = dij[0] * nx + dij[1] * ny
                if spacing < 0.0:
                    spacing = -spacing
            cl12.append(float(spacing))
            ja.append(j)

    iac = np.asarray(iac, dtype=np.int64)
    ja = np.asarray(ja, dtype=np.int64)
    cl12 = np.asarray(cl12, dtype=np.float64)
    hwva = np.asarray(hwva, dtype=np.float64)
    nja = int(iac.sum())

    assert nja == ja.size, f"nja ({nja}) != len(ja) ({ja.size})"

    if validate:
        pos = 0
        n = len(iac)
        for i in range(n):
            k = iac[i]
            assert k >= 1, f"iac[{i}] < 1"
            assert ja[pos] == i, f"Row {i} header in ja is {ja[pos]}, expected {i}"
            assert cl12[pos] == 0.0, f"cl12 self-entry at row {i} must be 0"
            assert hwva[pos] == 0.0, f"hwva self-entry at row {i} must be 0"
            if k > 1 and (hwva[pos + 1:pos + k] <= tol).any():
                bad_js = ja[pos + 1:pos + k][hwva[pos + 1:pos + k] <= tol]
                raise ValueError(
                    f"Zero/short faces found from cell {i} to neighbors {bad_js.tolist()} (tol={tol})"
                )
            pos += k
        assert pos == ja.size == cl12.size == hwva.size

        neigh = [set() for _ in range(n)]
        pos = 0
        for i in range(n):
            k = iac[i]
            row = ja[pos:pos + k]
            for j in row[1:]:
                neigh[i].add(int(j))
            pos += k

        asym = []
        for i in range(n):
            for j in neigh[i]:
                if i not in neigh[j]:
                    asym.append((i, j))
        if asym:
            raise ValueError(f"Asymmetric neighbor pairs detected (count={len(asym)}), e.g. {asym[:10]}")

        pos = 0
        edge_map = {}
        for i in range(n):
            k = iac[i]
            row_js = ja[pos:pos + k]
            row_c = cl12[pos:pos + k]
            row_h = hwva[pos:pos + k]
            for j, c, h in zip(row_js[1:], row_c[1:], row_h[1:], strict=False):
                key = (min(i, int(j)), max(i, int(j)))
                if key in edge_map:
                    c0, h0 = edge_map[key]
                    if not (abs(c - c0) <= 1e-9 * max(1.0, c0) and abs(h - h0) <= 1e-9 * max(1.0, h0)):
                        raise ValueError(f"Non-symmetric cl12/hwva for edge {key}: ({c0},{h0}) vs ({c},{h})")
                else:
                    edge_map[key] = (float(c), float(h))
            pos += k

    return iac, ja, cl12, hwva, nja
