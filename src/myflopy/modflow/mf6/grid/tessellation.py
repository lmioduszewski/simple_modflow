"""Voronoi cells that tile their domain: FloPy's ``tri2vor`` checked and repaired.

FloPy's ``tri2vor`` builds each cell from a Qhull Voronoi diagram of EVERY Triangle node,
drops the Voronoi vertices a ray-cast test puts outside the domain, and orders what is
left by angle about the vertices' MEAN. Ledger 316 measured what that costs:

- a cell whose node IS a reflex corner is not convex, its mean lands beside the node, and
  on the notch side the corner vertex is sorted to the wrong place -- a hole inside the
  boundary, cell area outside it, and a face MF6 never connects (327 of 504 corners on
  ``GridSpec.voronoi`` defaults, every one at 270 deg);
- a Voronoi vertex exactly ON a right-hand or top boundary edge is dropped, and the
  interior cell that needed it loses a triangle (the L at max_area 10000: seven faces,
  MF6 flow 5.12% low);
- where that vertex is KEPT it duplicates the boundary midpoint under another index, so
  two correct polygons share no index pair and MF6 drops their face;
- across a thin notch, Qhull's Delaunay (of all nodes) is not Triangle's, and cells take
  vertices from triangles the mesh does not have;
- a boundary that is not ONE simple polygon keeps cells in its first part only
  (ledger 317).

The TRUTH used here is the dual of Triangle's own conforming Delaunay mesh
(:func:`dual_cells`): a node's cell is the circumcentres of its triangles walked in order
around it, plus -- for a boundary node -- the midpoints of its two boundary edges and the
node itself. A certificate proves the result tiles the meshed domain exactly.

:func:`tessellate` keeps FloPy's output wherever it is right, index for index, and
rewrites only the cells that are not, in FloPy's own vertex ids. A grid with nothing
wrong therefore comes back byte-identical -- the same ``grid_pin``, the same ``.disv``
bytes -- and only a grid that had a broken cell re-pins. Where FloPy's cell SET is wrong
(a part of the domain lost, a stray node kept as a cell) nothing can be repaired in place,
and the dual's cells are taken whole: renumbered, and logged at WARNING.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from myflopy._logging import get_logger

logger = get_logger(__name__)

TWO_PI = 2.0 * math.pi

#: Two dual vertices closer than this fraction of their mesh edge are ONE vertex:
#: cocircular triangles, or a right triangle's circumcentre on its boundary hypotenuse.
#: Measured over 349 grids: merged pairs were at most 8.05e-12 of the edge, the shortest
#: kept dual segment 2.63e-6 -- this sits orders of magnitude from both.
MERGE_RTOL = 1e-9

#: How far (as a fraction of the circumradius) FloPy's vertex may sit from a triangle's
#: exact circumcentre and still be taken as that circumcentre. Qhull merges
#: near-cocircular facets, and its roundoff grows with the coordinates: at state-plane
#: magnitudes, measured 1e-6 to 1.4e-2 of the circumradius. A mesh that is not globally
#: Delaunay across a thin notch puts the nearest FloPy vertex 0.24-0.34 away -- a
#: different point, which is appended instead.
QHULL_MERGE_RTOL = 0.05

#: Coincidence tolerance, in units of the coordinates' own spacing. Circumcentre against
#: FloPy's vertex, and a kept on-boundary vertex against its midpoint twin, measured 1-25
#: ulp apart; genuinely distinct vertices never closer than 3.5e-5 ft.
COINCIDENT_ULPS = 1000.0


class NotConformingDelaunay(ValueError):
    """A triangulation whose dual is not a tiling: a circumcentre on the wrong side.

    Triangle's ``-D`` switch (which FloPy always passes) guarantees the property this
    tests; a mesh that fails it came from other switches or was edited by hand, and its
    Voronoi cells overlap or leave gaps whoever builds them.
    """


@dataclass(frozen=True)
class TessellationReport:
    """What :func:`tessellate` did to FloPy's cells.

    Attributes
    ----------
    rewritten
        ``{cell: reason}`` for every cell whose ring was rewritten in place, on FloPy's
        (unchanged) numbering. Reasons: ``"corner-order"`` (a boundary cell's vertices
        out of order -- the reflex-corner defect), ``"missing-vertex"`` (a vertex FloPy
        dropped), ``"twin-vertex"`` (the right polygon, but on a duplicate vertex id MF6
        cannot match), ``"foreign-vertex"`` (a vertex from a triangle the mesh does not
        have), ``"interior-order"``.
    appended
        Vertices added because FloPy had none at a location a cell needs.
    adopted
        ``True`` when FloPy's cell SET was wrong and the dual's cells were taken whole --
        the grid is renumbered relative to FloPy's. ``reason`` says why.
    legacy_pin
        ``grid_pin`` of FloPy's own cells when the repair changed any ring, else ``None``.
        A package pinned to this was resolved on the unrepaired grid.
    """

    rewritten: dict = field(default_factory=dict)
    appended: int = 0
    adopted: bool = False
    reason: str | None = None
    flopy_ncpl: int | None = None
    ncpl: int | None = None
    legacy_pin: str | None = None

    @property
    def changed(self) -> bool:
        """Whether the grid differs from FloPy's at all."""

        return bool(self.rewritten) or self.adopted or self.appended > 0

    def summary(self) -> str:
        """One sentence: what was repaired, or that nothing needed it."""

        if self.adopted:
            return (f"FloPy's cells were replaced whole ({self.reason}): {self.flopy_ncpl} "
                    f"cells became {self.ncpl}, renumbered")
        if not self.rewritten:
            return "FloPy's cells tile the domain; nothing was rewritten"
        why = Counter(self.rewritten.values())
        return (f"{len(self.rewritten)} of {self.ncpl} cells rewritten in place ("
                + ", ".join(f"{n} {r}" for r, n in sorted(why.items()))
                + f"); {self.appended} vertices appended; no cell renumbered")


# --------------------------------------------------------------------------------------
# the dual of Triangle's mesh


def _circumcentres(P, T):
    """Circumcentres of triangles ``T`` over points ``P``, computed relative to each corner."""

    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    ba, ca = b - a, c - a
    d = 2.0 * (ba[:, 0] * ca[:, 1] - ba[:, 1] * ca[:, 0])
    b2, c2 = (ba**2).sum(1), (ca**2).sum(1)
    ux = (ca[:, 1] * b2 - ba[:, 1] * c2) / d
    uy = (ba[:, 0] * c2 - ca[:, 0] * b2) / d
    return a + np.column_stack([ux, uy])


class _UnionFind:
    def __init__(self, n):
        self.parent = np.arange(n)

    def find(self, i):
        parent = self.parent
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root

    def union(self, i, j):
        ri, rj = self.find(i), self.find(j)
        if ri != rj:
            self.parent[max(ri, rj)] = min(ri, rj)


def dual_cells(tri, *, straight_sin: float = 1e-5):
    """``(verts, iverts, points)``: the Voronoi cells of Triangle's own mesh.

    Reads only ``tri.verts`` and ``tri.iverts`` -- no Qhull, no angle sort, no
    point-in-polygon test, no boundary markers -- so holes and a domain in several parts
    need no special handling: a boundary edge is a mesh edge with one triangle.

    One cell per Triangle node, in node order; a node with no triangle (a point outside
    the domain) gets none and its point is dropped. Rings are clockwise, FloPy's
    convention, and vertex ids are shared between neighbours, so MF6 DISV finds each face
    as a shared vertex pair. A boundary node's own vertex is left out where the boundary
    is straight through it, by FloPy's own test.

    Raises :class:`NotConformingDelaunay` unless every interior edge is locally Delaunay
    and every boundary triangle has its circumcentre on the domain side -- together they
    make every ring star-shaped from its node, and so the cells an exact tiling. Raises
    it too where the domain touches itself at a node (two parts meeting at a corner, a hole
    touching the edge): that node's cell is two pieces, which one DISV ring cannot hold.
    """

    P_abs = np.asarray(tri.verts, dtype=float)
    origin = P_abs.min(axis=0) if len(P_abs) else np.zeros(2)
    # Shifted to the mesh's own origin, so the relative tolerances below mean the same at
    # state-plane coordinates as at 0,0.
    P = P_abs - origin
    T = np.asarray(tri.iverts, dtype=np.int64).reshape(-1, 3).copy()
    N, M = len(P), len(T)

    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    if np.any(cross == 0):
        raise NotConformingDelaunay(f"{int(np.sum(cross == 0))} zero-area triangle(s) in the mesh")
    flip = cross < 0
    T[flip, 1], T[flip, 2] = T[flip, 2].copy(), T[flip, 1].copy()
    cc = _circumcentres(P, T)

    # Directed edges: e = 3t + i runs T[t,i] -> T[t,i+1], and is also the CORNER of node
    # T[t,i] in triangle t. Its twin runs the other way in the neighbouring triangle.
    start = T.ravel()
    end = T[:, [1, 2, 0]].ravel()
    key = start * N + end
    if len(np.unique(key)) != len(key):
        raise NotConformingDelaunay("a directed edge appears twice: the mesh is not a manifold triangulation")
    order = np.argsort(key, kind="stable")
    rkey = end * N + start
    pos = np.minimum(np.searchsorted(key[order], rkey), len(key) - 1)
    twin = np.where(key[order][pos] == rkey, order[pos], -1)
    tri_of = np.arange(3 * M) // 3
    prev_edge = 3 * tri_of + (np.arange(3 * M) % 3 + 2) % 3
    nxt = twin[prev_edge]  # the node's next corner counter-clockwise; -1 at the boundary
    boundary = np.flatnonzero(twin < 0)  # domain on the left
    interior = np.flatnonzero(twin > np.arange(3 * M))

    d = P[end] - P[start]
    L = np.hypot(d[:, 0], d[:, 1])
    rot = np.column_stack([-d[:, 1], d[:, 0]]) / L[:, None]
    s_int = ((cc[tri_of[interior]] - cc[tri_of[twin[interior]]]) * rot[interior]).sum(1)
    s_bnd = ((cc[tri_of[boundary]] - P[start[boundary]]) * rot[boundary]).sum(1)
    bad_int = s_int < -MERGE_RTOL * L[interior]
    bad_bnd = s_bnd < -MERGE_RTOL * L[boundary]
    if bad_int.any() or bad_bnd.any():
        raise NotConformingDelaunay(
            f"{int(bad_int.sum())} interior edge(s) not locally Delaunay and {int(bad_bnd.sum())} "
            "boundary edge(s) with their triangle's circumcentre outside the domain. Triangle's "
            "-D switch (conforming Delaunay, which FloPy always passes) guarantees neither "
            "happens; check any additional_args given to Triangle."
        )

    # Vertex identities: circumcentres, merged where they coincide; a right triangle's
    # circumcentre on its boundary hypotenuse IS that edge's midpoint.
    B = len(boundary)
    uf = _UnionFind(M + B)
    for e in interior[np.abs(s_int) <= MERGE_RTOL * L[interior]]:
        uf.union(int(tri_of[e]), int(tri_of[twin[e]]))
    for k in np.flatnonzero(np.abs(s_bnd) <= MERGE_RTOL * L[boundary]):
        uf.union(int(tri_of[boundary[k]]), M + int(k))
    roots = np.array([uf.find(i) for i in range(M + B)])
    mids = 0.5 * (P[start[boundary]] + P[end[boundary]])
    _, group = np.unique(roots, return_inverse=True)
    G = int(group.max()) + 1 if len(group) else 0
    xy = np.zeros((G, 2))
    np.add.at(xy, group, np.vstack([cc, mids]))
    xy /= np.bincount(group, minlength=G)[:, None]
    mcount = np.bincount(group[M:], minlength=G)
    if mcount.any():  # a group holding a midpoint sits exactly on the boundary
        mid_xy = np.zeros((G, 2))
        np.add.at(mid_xy, group[M:], mids)
        on_edge = mcount > 0
        xy[on_edge] = mid_xy[on_edge] / mcount[on_edge, None]
    tri_vid = group[:M]
    mid_vid = np.full(3 * M, -1)
    mid_vid[boundary] = group[M:]
    bnodes = np.unique(start[boundary])
    node_vid = np.full(N, -1)
    node_vid[bnodes] = G + np.arange(len(bnodes))
    all_xy = np.vstack([xy, P[bnodes]])

    corners_of = [[] for _ in range(N)]
    for cid, n in enumerate(start):
        corners_of[n].append(cid)
    rings, kept = [], []
    for n in range(N):
        cs = corners_of[n]
        if not cs:
            continue
        starts = [x for x in cs if twin[x] < 0]
        if len(starts) > 1:
            x0, y0 = P_abs[n]
            raise NotConformingDelaunay(
                f"the domain touches itself at ({x0:.6g}, {y0:.6g}): two parts meet at a "
                "corner, or a hole touches the edge. The cell there would be two pieces, "
                "which a DISV cell cannot be -- separate them or merge them."
            )
        if starts:
            walk = [starts[0]]
            while nxt[walk[-1]] >= 0:
                walk.append(int(nxt[walk[-1]]))
            if len(walk) != len(cs):
                raise NotConformingDelaunay(f"node {n}: its triangles do not form one fan")
            m0, m1 = int(mid_vid[walk[0]]), int(mid_vid[prev_edge[walk[-1]]])
            ring = [m0] + [int(tri_vid[tri_of[x]]) for x in walk] + [m1]
            p0, p2 = all_xy[m1] - P[n], all_xy[m0] - P[n]
            sinang = (p0[0] * p2[1] - p0[1] * p2[0]) / (np.hypot(*p0) * np.hypot(*p2))
            if not (abs(sinang) < straight_sin and (p0 @ p2) < 0):
                ring = [int(node_vid[n])] + ring
        else:
            walk = [cs[0]]
            while True:
                x = int(nxt[walk[-1]])
                if x == walk[0]:
                    break
                if x < 0:
                    raise NotConformingDelaunay(f"interior node {n}: its fan does not close")
                walk.append(x)
            if len(walk) != len(cs):
                raise NotConformingDelaunay(f"node {n}: its triangles do not form one fan")
            ring = [int(tri_vid[tri_of[x]]) for x in walk]
        ring = [x for i, x in enumerate(ring) if x != ring[i - 1]] if len(ring) > 1 else ring
        if len(ring) < 3:
            continue
        rings.append(ring[::-1])  # clockwise
        kept.append(n)

    new = np.full(len(all_xy), -1)
    count = 0
    for r in rings:
        for x in r:
            if new[x] < 0:
                new[x] = count
                count += 1
    used = new >= 0
    verts = np.zeros((count, 2))
    verts[new[used]] = all_xy[used] + origin
    iverts = [[int(new[x]) for x in r] for r in rings]
    return verts, iverts, P_abs[np.asarray(kept, dtype=np.int64)]


# --------------------------------------------------------------------------------------
# FloPy's cells, repaired in place


class _NotRepairable(Exception):
    """FloPy's output does not have the shape the in-place repair relies on."""


def _clusters(V, tol):
    """A location id per vertex: the smallest index among vertices within ``tol`` of it."""

    parent = np.arange(len(V))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in cKDTree(V).query_pairs(tol):
        ri, rj = root(i), root(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    return np.array([root(i) for i in range(len(V))], dtype=np.int64)


def _near180_filter(ring, V):
    """FloPy ``get_sorted_vertices``' straight-angle filter, so a rewritten ring drops what FloPy drops."""

    out, n = [], len(ring)
    for i in range(n):
        x0, y0 = V[ring[i - 1]]
        x1, y1 = V[ring[i]]
        x2, y2 = V[ring[(i + 1) % n]]
        s0x, s0y, s2x, s2y = x0 - x1, y0 - y1, x2 - x1, y2 - y1
        m = math.hypot(s0x, s0y) * math.hypot(s2x, s2y)
        if m > 0 and abs((s0x * s2y - s0y * s2x) / m) < 1e-5 and s0x * s2x + s0y * s2y < 0:
            continue
        out.append(ring[i])
    return out


def _signed_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _repair_in_place(verts, iverts, points, tri):
    """``(iverts, verts, rewritten, appended)``: FloPy's cells, wrong ones rewritten.

    Every cell's ring is compared against what the dual of Triangle's mesh says it holds
    -- the circumcentres of its node's triangles and, on the boundary, its two edge
    midpoints and the node -- with sets and angles, no polygon clipping. A cell that
    passes keeps its ring exactly as FloPy wrote it; one that does not is rebuilt by
    walking the triangle fan around its node, on ONE FloPy vertex id per location (the
    one most rings already use there), so a face is the same index pair in both cells.
    Raises :class:`_NotRepairable` when FloPy's output is not laid out as expected.
    """

    V = np.asarray(verts, dtype=float)
    P = np.asarray(tri.verts, dtype=float)
    T = np.asarray(tri.iverts, dtype=np.int64)
    ncell = len(iverts)
    tol = COINCIDENT_ULPS * float(np.spacing(max(float(np.abs(P).max()), 1.0)))

    pts = np.asarray(points, dtype=float)
    if ncell == len(P):
        node_of = np.arange(ncell)
        if not np.array_equal(pts, P):
            raise _NotRepairable("cell centres are not Triangle's nodes")
    else:
        key = {(float(x), float(y)): i for i, (x, y) in enumerate(P)}
        try:
            node_of = np.array([key[(float(x), float(y))] for x, y in pts], dtype=np.int64)
        except KeyError as exc:
            raise _NotRepairable("a cell centre is not a Triangle node") from exc
    cell_of = np.full(len(P), -1, dtype=np.int64)
    cell_of[node_of] = np.arange(ncell)

    # The boundary edges tri2vor closes cells on, and its vertex layout for them:
    # [kept Voronoi vertices..., (midpoint, endpoint1, endpoint2) per boundary edge].
    marker = tri.edge["boundary_marker"]
    try:
        nexterior = len(tri._polygons[0]) + sum(len(tri._polygons[i + 1]) for i in range(len(tri._holes)))
    except (IndexError, TypeError) as exc:
        raise _NotRepairable("the triangulation's holes are not polygons") from exc
    be = tri.edge[(marker > 0) & (marker <= nexterior)]
    nbe = len(be)
    base = len(V) - 3 * nbe
    e1, e2 = be["endpoint1"].astype(np.int64), be["endpoint2"].astype(np.int64)
    mid = base + 3 * np.arange(nbe)
    if base < 0 or not (np.array_equal(V[mid], (P[e1] + P[e2]) / 2) and np.array_equal(V[mid + 1], P[e1])
                        and np.array_equal(V[mid + 2], P[e2])):
        raise _NotRepairable("tri2vor's vertex layout is not the one this repair reads")

    cc = _circumcentres(P, T)
    loc = _clusters(V, tol)
    tree = cKDTree(V)
    dist, nearest = tree.query(cc)
    cc_loc = np.where(dist <= tol, loc[nearest], -1)
    appended_rows = []
    unmatched = np.flatnonzero(cc_loc < 0)
    if len(unmatched):
        # Qhull merged near-cocircular facets (keep its vertex: the three cells already
        # share it), or the mesh is not globally Delaunay across a thin notch (FloPy has
        # no vertex there at all: append the circumcentre).
        radius = np.hypot(*(P[T[unmatched, 0]] - cc[unmatched]).T)
        ring_locs = {}
        for t, R in zip(unmatched, radius, strict=True):
            common = None
            for q in T[t]:
                j = int(cell_of[q])
                if j < 0:
                    continue
                if j not in ring_locs:
                    ring_locs[j] = {int(loc[v]) for v in iverts[j]}
                common = set(ring_locs[j]) if common is None else common & ring_locs[j]
            best = min(common, key=lambda x: float(np.hypot(*(V[x] - cc[t])))) if common else None
            if best is not None and float(np.hypot(*(V[best] - cc[t]))) < QHULL_MERGE_RTOL * R:
                cc_loc[t] = best
            else:
                appended_rows.append(int(t))
        if appended_rows:
            rows = np.asarray(appended_rows)
            cc_loc[rows] = len(V) + np.arange(len(rows))
            loc = np.concatenate([loc, len(V) + np.arange(len(rows))])
            V = np.vstack([V, cc[rows]])
    NL = len(V)

    lens = np.fromiter((len(r) for r in iverts), dtype=np.int64, count=ncell)
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)
    flat = np.fromiter((v for r in iverts for v in r), dtype=np.int64, count=int(lens.sum()))
    cellid = np.repeat(np.arange(ncell), lens)
    floc = loc[flat]

    # the canonical id at each location: the most used, ties to the smallest
    uses = np.bincount(flat, minlength=NL)
    order = np.lexsort((np.arange(NL), -uses, loc))
    first = np.ones(NL, dtype=bool)
    first[1:] = loc[order][1:] != loc[order][:-1]
    canon = np.full(NL, -1, dtype=np.int64)
    canon[loc[order][first]] = order[first]

    own_loc = np.full(ncell, -1, dtype=np.int64)  # a boundary cell's own node location
    ok1, ok2 = cell_of[e1] >= 0, cell_of[e2] >= 0
    own_loc[cell_of[e1[ok1]]] = loc[mid[ok1] + 1]
    own_loc[cell_of[e2[ok2]]] = loc[mid[ok2] + 2]
    is_own = floc == own_loc[cellid]

    expected_cell = np.concatenate([cell_of[T[:, 0]], cell_of[T[:, 1]], cell_of[T[:, 2]], cell_of[e1], cell_of[e2]])
    expected_loc = np.concatenate([cc_loc, cc_loc, cc_loc, loc[mid], loc[mid]])
    keep = expected_cell >= 0
    expected = np.unique(expected_cell[keep] * NL + expected_loc[keep])
    held, held_count = np.unique(cellid[~is_own] * NL + floc[~is_own], return_counts=True)
    flags = {}
    for k in np.unique(np.setdiff1d(held, expected) // NL):
        flags[int(k)] = "foreign-vertex"
    for k in np.unique(np.setdiff1d(expected, held) // NL):
        flags[int(k)] = "missing-vertex"
    twin_id = (~is_own) & (canon[floc] != flat)
    for k in np.unique(np.concatenate([cellid[twin_id], held[held_count > 1] // NL])):
        flags.setdefault(int(k), "twin-vertex")

    # interior cells: clockwise about the node is exactly one descent of atan2(-dy, dx)
    node_xy = P[node_of][cellid]
    ang = np.arctan2(node_xy[:, 1] - V[flat, 1], V[flat, 0] - node_xy[:, 0])
    after = np.arange(len(flat)) + 1
    after[starts + lens - 1] = starts
    descents = np.add.reduceat((ang[after] < ang).astype(np.int64), starts)
    for k in np.flatnonzero((own_loc < 0) & (descents != 1)):
        flags.setdefault(int(k), "interior-order")

    tv = T.ravel()
    torder = np.argsort(tv, kind="stable")
    tptr = np.searchsorted(tv[torder], np.arange(len(P) + 1))

    def incident(p):
        return (torder[tptr[p]:tptr[p + 1]] // 3).tolist()

    bnode = defaultdict(list)  # node -> [(other endpoint, midpoint id, own-vertex id)]
    for j in range(nbe):
        bnode[int(e1[j])].append((int(e2[j]), int(mid[j]), int(mid[j] + 1)))
        bnode[int(e2[j])].append((int(e1[j]), int(mid[j]), int(mid[j] + 2)))

    def first_last(p):
        """(first, last) midpoint of boundary node p: the domain lies clockwise from first."""
        entries = bnode[p]
        if len(entries) != 2:
            return None
        (qa, ma, _), (_qb, mb, _) = entries
        shared = [t for t in incident(p) if qa in T[t]]
        if len(shared) != 1:
            return None
        qc = next(int(q) for q in T[shared[0]] if q != p and q != qa)
        u, w = P[qa] - P[p], P[qc] - P[p]
        return (ma, mb) if u[0] * w[1] - u[1] * w[0] < 0 else (mb, ma)

    # boundary cells: [first midpoint, ...clockwise about the node..., last midpoint, (node)]
    fl_of = {}
    for p in bnode:
        k = int(cell_of[p])
        if k < 0:
            continue
        fl = first_last(p)
        if fl is None:
            raise _NotRepairable(f"boundary node {p} does not have exactly two boundary edges")
        fl_of[p] = fl
        if k in flags:
            continue
        ring = list(iverts[k])
        rl = [int(loc[v]) for v in ring]
        lf, ll, own = int(loc[fl[0]]), int(loc[fl[1]]), int(own_loc[k])
        if lf not in rl:
            flags[k] = "missing-vertex"
            continue
        s = rl.index(lf)
        rr, rrl = ring[s:] + ring[:s], rl[s:] + rl[:s]
        body = rr[:-1] if rrl[-1] == own else rr
        good = int(loc[body[-1]]) == ll and own not in rrl[:-1]
        if good:
            px, py = P[p]
            a0 = math.atan2(py - V[fl[0], 1], V[fl[0], 0] - px)
            phis = [(math.atan2(py - V[v, 1], V[v, 0] - px) - a0) % TWO_PI for v in body[1:]]
            good = all(0.0 < x <= y for x, y in zip(phis, phis[1:], strict=False)) and (not phis or phis[0] > 0.0)
        if not good:
            flags[k] = "corner-order"

    new_iverts = list(iverts)
    for k, _reason in flags.items():
        p = int(node_of[k])
        tris = incident(p)
        fl = fl_of.get(p)
        across = defaultdict(list)
        for t in tris:
            for q in T[t]:
                if q != p:
                    across[int(q)].append(t)

        def walk(t, q, p=p, across=across):
            seq, seen = [], set()
            while t is not None and t not in seen:
                seen.add(t)
                seq.append(int(cc_loc[t]))
                q = next(int(r) for r in T[t] if r != p and r != q)
                t = next((s for s in across[q] if s != t), None)
            return seq

        if fl is not None:
            q0 = next(qa for qa, m, _ in bnode[p] if m == fl[0])
            seq = [int(loc[fl[0]])] + walk(across[q0][0], q0) + [int(loc[fl[1]])]
        else:
            t0 = tris[0]
            seq = walk(t0, next(int(r) for r in T[t0] if r != p))
        locs = []
        for x in seq:
            if not locs or locs[-1] != x:
                locs.append(x)
        if fl is None and len(locs) > 1 and locs[0] == locs[-1]:
            locs.pop()
        ring = [int(canon[x]) for x in locs]
        if fl is not None:
            own = int(own_loc[k])
            ring.append(next((int(v) for v in iverts[k] if int(loc[v]) == own), bnode[p][0][2]))
        ring = _near180_filter(ring, V)
        if _signed_area(V[ring]) > 0:
            ring = ring[::-1]
        new_iverts[k] = ring
    rewritten = {k: r for k, r in sorted(flags.items()) if list(new_iverts[k]) != list(iverts[k])}
    return new_iverts, V, rewritten, len(appended_rows)


def tessellate(tri, verts, iverts, points) -> tuple[np.ndarray, list, np.ndarray, TessellationReport]:
    """FloPy's ``tri2vor`` output made to tile the domain: ``(verts, iverts, points, report)``.

    ``verts``/``iverts``/``points`` are what ``tri2vor(tri)`` returned. Cells FloPy got
    right come back exactly as it wrote them; the rest are rewritten on its vertex ids,
    appending a vertex only where FloPy has none. When FloPy's cell SET differs from the
    dual's -- a domain part lost, a node outside the domain kept as a cell, cells Qhull
    dropped -- the dual's cells are returned whole, renumbered, and a WARNING is logged.
    Raises :class:`NotConformingDelaunay` for a mesh whose dual is not a tiling.
    """

    d_verts, d_iverts, d_points = dual_cells(tri)
    flopy_points = np.asarray(points, dtype=float)
    ncell = len(iverts)

    def adopt(reason):
        logger.warning(
            "the Voronoi grid's cells were rebuilt from Triangle's mesh (%s): FloPy's %d "
            "cells became %d, renumbered (ledgers 316, 317)", reason, ncell, len(d_iverts))
        return d_verts, d_iverts, d_points, TessellationReport(
            adopted=True, reason=reason, flopy_ncpl=ncell, ncpl=len(d_iverts))

    if flopy_points.shape != d_points.shape or not np.array_equal(flopy_points, d_points):
        lost = len({tuple(p) for p in d_points.tolist()} - {tuple(p) for p in flopy_points.tolist()})
        extra = len({tuple(p) for p in flopy_points.tolist()} - {tuple(p) for p in d_points.tolist()})
        return adopt(f"FloPy dropped {lost} cell(s) the mesh has and kept {extra} it has not")
    try:
        new_iverts, new_verts, rewritten, appended = _repair_in_place(verts, iverts, flopy_points, tri)
    except _NotRepairable as exc:
        return adopt(str(exc))
    report = TessellationReport(rewritten=rewritten, appended=appended, flopy_ncpl=ncell, ncpl=ncell)
    if report.changed:
        logger.info("Voronoi grid repaired: %s", report.summary())
    return new_verts, new_iverts, flopy_points, report
