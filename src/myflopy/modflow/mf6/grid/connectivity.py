from __future__ import annotations

import numpy as np
from shapely.geometry import LineString, MultiLineString
from shapely.prepared import prep

from myflopy._logging import get_logger

logger = get_logger(__name__)


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

    logger.info('building connectivity (iac, ja, cl12, hwva, nja) for %d cells',
                len(adjacent_cells_idx))
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


class DisvFaces:
    """MODFLOW 6's own horizontal-connection geometry for a DISV grid, per cell pair.

    ``grid.cl12`` is the full centre-to-centre spacing (the loop above measures it
    between polygon CENTROIDS), which is not what MF6 means by a connection length.
    A DISV connection's length on each side is the PERPENDICULAR distance from that
    cell's DISV centre (``cell2d``'s ``xc``/``yc``) to the LINE through the two
    vertices the cells share, and its width is the distance between those vertices
    (``DisvGeom.f90``: ``shared_edge``, ``distance_normal``). A LAK HORIZONTAL
    connection written into a cell is one side of such a face, so its ``CONNLEN`` is
    that cell's own distance -- measured 1.47-3.44x shorter than ``cl12`` on a real
    Voronoi pond (median 1.975), and the choice that reproduced a pit lake's
    through-flow to 0.5% (ledger 325).

    Built from :meth:`VoronoiGridPlus.get_disv_gridprops` -- the arrays the DISV is
    written from -- so the numbers are the ones MF6 will compute from the file.

    Parameters
    ----------
    gridprops : dict
        ``get_disv_gridprops()``'s result: ``"vertices"`` rows ``[iv, x, y]`` and
        ``"cell2d"`` rows ``[icell, xc, yc, ncvert, iv1, iv2, ...]``.

    Examples
    --------
    >>> import numpy as np
    >>> from myflopy.modflow.mf6.grid.connectivity import DisvFaces
    >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 2.0, 0.0],
    ...                       [3, 0.0, 1.0], [4, 1.0, 1.0], [5, 2.0, 1.0]],
    ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 3, 4, 1], [1, 1.25, 0.5, 4, 1, 4, 5, 2]]}
    >>> faces = DisvFaces(props)
    >>> faces.length_and_width(0, 1, toward=1)   # cell 1's centre is 0.25 from x = 1
    (0.25, 1.0)
    """

    def __init__(self, gridprops) -> None:
        self._xy = {int(row[0]): (float(row[1]), float(row[2])) for row in gridprops["vertices"]}
        self._rings: dict[int, list[int]] = {}
        self._centres: dict[int, tuple[float, float]] = {}
        for row in gridprops["cell2d"]:
            row = list(row)
            cell, count = int(row[0]), int(row[3])
            self._centres[cell] = (float(row[1]), float(row[2]))
            ring = [int(vertex) for vertex in row[4 : 4 + count]]
            if len(ring) > 1 and ring[0] == ring[-1]:
                ring = ring[:-1]
            self._rings[cell] = ring

    def shared_edge(self, a: int, b: int) -> tuple[int, int] | None:
        """The two vertices cells ``a`` and ``b`` share as an EDGE, or ``None``.

        MF6's rule: walk ``a``'s ring in order and take the first consecutive pair
        that appears reversed in ``b``'s. A pair that appears in the same order (a
        mesh whose rings are not consistently oriented) is accepted too -- MF6 would
        find no face there and connect nothing, which is a mesh defect this cannot
        repair and should not hide by returning ``None`` for a face that exists.

        Parameters
        ----------
        a : int
            One cell (a zero-based ``cell2d`` number); its ring order decides which
            of several shared edges is taken, as in MF6.
        b : int
            The other cell.

        Returns
        -------
        tuple of (int, int), or None
            The two vertex numbers, in ``a``'s ring order; ``None`` when the cells
            share no edge (a corner or nothing).

        Examples
        --------
        >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
        ...                       [4, 2.0, 0.0], [5, 2.0, 1.0]],
        ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 2, 3, 1], [1, 1.5, 0.5, 4, 1, 3, 5, 4]]}
        >>> DisvFaces(props).shared_edge(0, 1)
        (3, 1)
        """

        ring_a, ring_b = self._rings[int(a)], self._rings[int(b)]
        edges_b = {
            (ring_b[i], ring_b[(i + 1) % len(ring_b)]) for i in range(len(ring_b))
        }
        for i in range(len(ring_a)):
            first, second = ring_a[i], ring_a[(i + 1) % len(ring_a)]
            if (second, first) in edges_b or (first, second) in edges_b:
                return first, second
        return None

    def length_and_width(self, a: int, b: int, *, toward: int) -> tuple[float, float] | None:
        """``(CONNLEN, CONNWIDTH)`` of the face ``a|b`` measured from cell ``toward``.

        Parameters
        ----------
        a : int
            One of the two cells sharing the face (a zero-based ``cell2d`` number).
        b : int
            The other one.
        toward : int
            Whose DISV centre the length is measured from -- for a LAK wall, the
            aquifer cell the connection is written into.

        Returns
        -------
        tuple of float, or None
            The perpendicular distance from ``toward``'s centre to the line of the
            shared edge, and the edge's length; ``None`` when the cells share no edge
            (only a corner) or the edge has no length.

        Examples
        --------
        >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
        ...                       [4, 2.0, 0.0], [5, 2.0, 1.0]],
        ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 2, 3, 1], [1, 1.25, 0.5, 4, 1, 3, 5, 4]]}
        >>> DisvFaces(props).length_and_width(0, 1, toward=0)
        (0.5, 1.0)
        """

        edge = self.shared_edge(a, b)
        if edge is None:
            return None
        (x1, y1), (x2, y2) = self._xy[edge[0]], self._xy[edge[1]]
        width = float(np.hypot(x2 - x1, y2 - y1))
        if width == 0.0:
            return None
        x0, y0 = self._centres[int(toward)]
        length = abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) / width
        return float(length), width
