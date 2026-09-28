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
- a vertex lying straight between its two neighbours in ONE cell is dropped from that cell
  by FloPy's 180-degree filter and kept by the cell beside it, which loses their face;
- a boundary that is not ONE simple polygon keeps cells in its first part only
  (ledger 317).

The TRUTH used here is the dual of Triangle's own conforming Delaunay mesh
(:func:`dual_cells`): a node's cell is the circumcentres of its triangles walked in order
around it, plus -- for a boundary node -- the midpoints of its two boundary edges and the
node itself. A certificate proves the result tiles the meshed domain exactly.

:func:`tessellate` gives every one of the dual's vertices a FloPy vertex id -- the one
FloPy has at that location, or an appended one where it has none -- so every ring it
returns IS the dual's ring, and no face a DISV file can hold is lost between two cells. A
FloPy ring that already is that ring is returned exactly as FloPy wrote it. A grid with
nothing wrong therefore comes back byte-identical -- the same ``grid_pin``, the same
``.disv`` bytes -- and only a grid that had a broken cell re-pins. Where FloPy's cell SET
is wrong (a part of the domain lost, a stray node kept as a cell) the dual's cells are
taken whole: renumbered, and logged at WARNING.

A corner is never split into two vertices the written file cannot tell apart
(:func:`~myflopy.grid_identity.written_resolution`). Four nodes almost on one circle give
two circumcentres a hair apart -- 7.4e-5 ft on a real 20,323-cell model at state plane,
where FloPy writes 0.01 ft -- and FloPy drew ONE corner there. Kept as two, the face
between them was written with zero length and MF6 returned NaN for every head under
normal termination. Such a vertex takes FloPy's corner, and two the repair would append
are one.
"""

from __future__ import annotations

import dataclasses
import types
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
from flopy.utils.cvfdutil import get_disv_gridprops
from scipy.spatial import cKDTree

from myflopy._logging import get_logger
from myflopy.grid_identity import grid_pin, written_resolution

logger = get_logger(__name__)

#: Two dual vertices closer than this fraction of their mesh edge are ONE vertex:
#: cocircular triangles, or a right triangle's circumcentre on its boundary hypotenuse.
#: Measured over 349 grids: merged pairs were at most 8.05e-12 of the edge, the shortest
#: kept dual segment 2.63e-6 -- this sits orders of magnitude from both.
MERGE_RTOL = 1e-9

#: Coincidence tolerance, in units of the coordinates' own spacing, used twice: FloPy's
#: vertex IS a dual vertex within it, and two circumcentres within it are one vertex.
#: Circumcentre against FloPy's vertex, and a kept on-boundary vertex against its midpoint
#: twin, measured 1-25 ulp apart on grids at the origin, where distinct vertices were never
#: closer than 3.5e-5 ft. At state-plane magnitudes the match still holds -- no grid FloPy
#: draws right at (1.2e6, 6.5e5) or (5e6, 5e6) had a ring rewritten
#: (``voronoi_tiling_sweep.py`` part F), and on a spiky star at 1.2e6, where FloPy's cells
#: ARE wrong, 96% of circumcentres sat on FloPy's vertex and the nearest of the rest was
#: 4.6e-4 ft away, 2,000 times the tolerance (``voronoi_tiling_gaps.py`` section 11) -- but
#: DISTINCT circumcentres do come closer than it there, and are merged (part E counts them).
#: The first repair took a vertex up to 5% of the circumradius away as the same one, and
#: lost faces that way.
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
        dropped: outside by its ray test, or straight between two others), ``"twin-vertex"``
        (the right polygon, but on a duplicate vertex id MF6 cannot match, or holding one
        location twice -- a zero-width connection), ``"foreign-vertex"`` (a vertex that is
        not one of this cell's: a circumcentre of a triangle the mesh does not have, a
        neighbour's across a thin notch, or Qhull's merge of two), ``"interior-order"``.
    appended
        Vertices added because FloPy had none at a location a cell needs.
    unsplit
        The dual's vertices joined to a neighbour across a face shorter than a DISV file
        FloPy writes can record (:func:`~myflopy.grid_identity.written_resolution`), so a
        corner FloPy drew once is not split in two. Kept apart, such a face can be written
        with zero length, and MF6 returns NaN for every head.
    split_pin
        ``grid_pin`` of the grid with those vertices kept apart -- what the repair built
        from bb56668 (2026-09-26) until this rule -- when that differs, else ``None``. A
        package pinned to it was resolved on that grid: same cells, same numbers.
    adopted
        ``True`` when FloPy's cell SET was wrong and the dual's cells were taken whole --
        the grid is renumbered relative to FloPy's. ``reason`` says why.
    legacy_pin
        ``grid_pin`` of FloPy's own cells when the repair changed any ring, else ``None``.
        A package pinned to this was resolved on the unrepaired grid.
    """

    rewritten: dict = field(default_factory=dict)
    appended: int = 0
    unsplit: int = 0
    split_pin: str | None = None
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

        kept = (f"; {self.unsplit} corners closer than the written file records kept as one"
                if self.unsplit else "")
        if self.adopted:
            return (f"FloPy's cells were replaced whole ({self.reason}): {self.flopy_ncpl} "
                    f"cells became {self.ncpl}, renumbered{kept}")
        if not self.rewritten:
            return "FloPy's cells tile the domain; nothing was rewritten" + kept
        why = Counter(self.rewritten.values())
        return (f"{len(self.rewritten)} of {self.ncpl} cells rewritten in place ("
                + ", ".join(f"{n} {r}" for r, n in sorted(why.items()))
                + f"); {self.appended} vertices appended; no cell renumbered{kept}")


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
    make every ring star-shaped from its node, and so the cells an exact tiling. (The two
    circumcentres of an edge both lie on its bisector, and a fan triangle's turn is the
    sign of their separation along it, so the certificate IS the orientation of every fan
    triangle whose vertices the merge below leaves where they were.) The certificate's
    tolerance is 1e-9 of the edge and no wider: widened to 1000 ulp, it passed a boundary
    circumcentre 2.1e-7 ft outside the domain and returned a cell crossing itself.
    Circumcentres closer than 1000 ulp of the coordinates (``COINCIDENT_ULPS``) are then ONE
    vertex, so no edge a few ulps long reaches a written grid; a merged vertex moves by up
    to that much, which can turn a fan triangle that small over -- a mesh measured doing so
    still gave valid cells. Because of that merge, a node whose triangles' circumcentres
    all fall that close together has no cell left: that, and two nodes closer than 1000
    ulp, are refused rather than dropped. So is a domain that touches itself at a node
    (two parts meeting at a corner, a hole touching the edge): that node's cell is two
    pieces, which one DISV ring cannot hold.
    """

    verts, iverts, points, _ = _dual(tri, straight_sin=straight_sin)
    return verts, iverts, points


def _dual(tri, *, straight_sin: float = 1e-5):
    """:func:`dual_cells`, plus whether each cell's node is on the boundary."""

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
    resolvable = COINCIDENT_ULPS * float(np.spacing(max(float(np.abs(P_abs).max()), 1.0)))
    if L.min() < resolvable:
        e = int(np.argmin(L))
        x0, y0 = P_abs[start[e]]
        raise NotConformingDelaunay(
            f"two of the mesh's nodes are {L[e]:.3g} apart at ({x0:.6g}, {y0:.6g}), closer than "
            f"1000 ulp of its coordinates ({resolvable:.3g}), inside which circumcentres are "
            "merged into one vertex, so the cells between them would have no area. (FloPy "
            "writes a vertex to nine significant digits, so cells far larger than that reach "
            "MF6 with zero area too.) A point constraint on a boundary edge, or within "
            "rounding of one, can make Triangle pack nodes that close; move it clear of the edge."
        )
    rot = np.column_stack([-d[:, 1], d[:, 0]]) / L[:, None]
    s_int = ((cc[tri_of[interior]] - cc[tri_of[twin[interior]]]) * rot[interior]).sum(1)
    s_bnd = ((cc[tri_of[boundary]] - P[start[boundary]]) * rot[boundary]).sum(1)
    # The certificate stays at 1e-9 of the edge; two circumcentres closer than 1000 ulp are
    # MERGED whatever the edge: kept apart, they write an edge a few ulps long that a
    # reprojection can turn over.
    bad_int = s_int < -MERGE_RTOL * L[interior]
    bad_bnd = s_bnd < -MERGE_RTOL * L[boundary]
    near_int = np.maximum(MERGE_RTOL * L[interior], resolvable)
    near_bnd = np.maximum(MERGE_RTOL * L[boundary], resolvable)
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
    for e in interior[np.abs(s_int) <= near_int]:
        uf.union(int(tri_of[e]), int(tri_of[twin[e]]))
    for k in np.flatnonzero(np.abs(s_bnd) <= near_bnd):
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
    rings, kept, on_boundary = [], [], []
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
            x0, y0 = P_abs[n]
            raise NotConformingDelaunay(
                f"the node at ({x0:.6g}, {y0:.6g}) has no cell left: the circumcentres of its "
                f"triangles lie within 1000 ulp of each other ({resolvable:.3g}), where they are "
                "one vertex. A point constraint on a boundary edge, or within rounding of one, "
                "can make Triangle pack nodes that close; move it clear of the edge."
            )
        rings.append(ring[::-1])  # clockwise
        kept.append(n)
        on_boundary.append(bool(starts))

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
    points = P_abs[np.asarray(kept, dtype=np.int64)]
    return verts, iverts, points, np.asarray(on_boundary, dtype=bool)


# --------------------------------------------------------------------------------------
# FloPy's cells, checked against the dual


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


def _short_faces(D, d_iverts, unit, apart=None):
    """A root per dual vertex: the vertices joined along every face shorter than one
    written unit on BOTH axes (``unit``, :func:`~myflopy.grid_identity.written_resolution`),
    which a DISV file can write as a single point -- the lowest index of each group.

    Two groups that each hold an ``apart`` vertex are not joined: FloPy drew both corners,
    and a face between them is FloPy's own, not a corner the repair split. A group whose
    joining would leave a cell with fewer than three corners, or one twice, stays split.
    """

    parent = np.arange(len(D))
    held = np.zeros(len(D), dtype=bool) if apart is None else np.array(apart, dtype=bool)

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if len(d_iverts):
        u = np.concatenate([np.asarray(r, dtype=np.int64) for r in d_iverts])
        v = np.concatenate([np.roll(np.asarray(r, dtype=np.int64), -1) for r in d_iverts])
        for e in np.flatnonzero((np.abs(D[u] - D[v]) < unit).all(axis=1) & (u != v)):
            ru, rv = root(int(u[e])), root(int(v[e]))
            if ru != rv and not (held[ru] and held[rv]):
                lo, hi = min(ru, rv), max(ru, rv)
                parent[hi] = lo
                held[lo] = held[lo] or held[hi]
    roots = np.array([root(i) for i in range(len(D))], dtype=np.int64)
    # A group that would leave a cell with fewer than three corners, or one corner twice,
    # stays split: that cell is smaller than the file can write whichever way it is drawn,
    # and the build-time check (``zero_length_faces``) says so for the file actually written
    # -- where more digits (``float_precision``) can still hold it. Joined, the GRID refused.
    undo = set()
    grouped = np.bincount(roots, minlength=len(D))[roots] > 1
    for ring in d_iverts if grouped.any() else ():
        if not grouped[list(ring)].any():
            continue
        rr = [int(roots[x]) for x in ring]
        rr = [x for i, x in enumerate(rr) if x != rr[i - 1]]
        if len(rr) < 3 or len(set(rr)) != len(rr):
            undo.update(int(roots[x]) for x in ring)
    if undo:
        split = np.isin(roots, list(undo))
        roots[split] = np.flatnonzero(split)
    return roots


def _rings_on(label, d_iverts, where):
    """The dual's rings under ``label``, a vertex repeated in a row dropped; refused when a
    cell is left with fewer than three corners or one corner twice."""

    rings = []
    for k, ring in enumerate(d_iverts):
        want = [int(label[x]) for x in ring]
        want = [x for i, x in enumerate(want) if x != want[i - 1]]
        if len(want) < 3 or len(set(want)) != len(want):
            x0, y0 = where[int(ring[0])]
            raise NotConformingDelaunay(
                f"cell {k} near ({x0:.6g}, {y0:.6g}) has corners closer together than a DISV "
                "file FloPy writes can record, so it would reach MF6 with no area or "
                "crossing itself. Nodes packed that close usually come from a point "
                "constraint on or near a boundary edge, or two points almost on top of each "
                "other; move them apart."
            )
        rings.append(want)
    return rings


def _onto_flopy(d_verts, d_iverts, verts, iverts, on_boundary):
    """``(verts, iverts, rewritten, appended, unsplit, split)``: the dual's cells on FloPy's
    vertex ids.

    Each of the dual's vertices is given the FloPy vertex at its location -- the one most
    of FloPy's rings use there, since FloPy can hold one location under several ids -- or,
    where FloPy has none, is appended. Every ring is then the dual's ring under that one
    relabelling, so two neighbours share their face as the same index pair by
    construction. A FloPy ring that already IS that ring is kept exactly as FloPy wrote it,
    starting vertex included; any other is replaced, and its reason recorded.

    One exception, the file's: the dual's vertices joined by a face shorter than one
    written unit on both axes (:func:`~myflopy.grid_identity.written_resolution`) are ONE
    corner, not split -- the FloPy vertex one of them holds, else a FloPy vertex no dual
    vertex holds that lies within a unit of all of them (FloPy's single corner, drawn a
    hair off), else one appended vertex. A face that short cannot be written: kept, it
    reached MF6 with zero length and NaN heads. ``unsplit`` counts the vertices joined, and
    ``split`` is ``(verts, iverts)`` with them kept apart -- the grid the repair built from
    bb56668 until this rule -- when there are any; else ``None``.
    """

    V = np.asarray(verts, dtype=float)
    D = np.asarray(d_verts, dtype=float)
    scale = max(float(np.abs(V).max()) if len(V) else 0.0, float(np.abs(D).max()) if len(D) else 0.0, 1.0)
    tol = COINCIDENT_ULPS * float(np.spacing(scale))
    ncell = len(iverts)

    lens = np.fromiter((len(r) for r in iverts), dtype=np.int64, count=ncell)
    flat = np.fromiter((v for r in iverts for v in r), dtype=np.int64, count=int(lens.sum()))
    loc = _clusters(V, tol)
    uses = np.bincount(flat, minlength=len(V))
    order = np.lexsort((np.arange(len(V)), -uses, loc))
    first = np.ones(len(V), dtype=bool)
    first[1:] = loc[order][1:] != loc[order][:-1]
    canon = np.full(len(V), -1, dtype=np.int64)
    canon[loc[order][first]] = order[first]

    # One FloPy id per dual vertex and one dual vertex per FloPy id, nearest pairs first,
    # so two dual vertices closer than the tolerance never collapse onto one id.
    ids = np.unique(canon[loc])
    near = cKDTree(V[ids]).query_ball_point(D, tol)
    pair_d = np.fromiter((i for i, js in enumerate(near) for _ in js), dtype=np.int64)
    pair_f = ids[np.fromiter((j for js in near for j in js), dtype=np.int64, count=len(pair_d))]
    dist = np.hypot(*(D[pair_d] - V[pair_f]).T)
    label = np.full(len(D), -1, dtype=np.int64)
    dual_at = {}  # FloPy canonical id -> the dual vertex holding it
    for n in np.argsort(dist, kind="stable"):
        i, c = int(pair_d[n]), int(pair_f[n])
        if label[i] < 0 and c not in dual_at:
            dual_at[c] = i
            label[i] = c
    # Not split: the dual's vertices joined by a face the written file cannot hold are one
    # corner -- FloPy's where it drew one, else one appended vertex.
    unit = written_resolution(np.vstack([V, D]))
    roots = _short_faces(D, d_iverts, unit, apart=label >= 0)
    extra = np.flatnonzero(label < 0)
    split, unsplit = None, len(D) - len(np.unique(roots))
    if unsplit:
        split_label = label.copy()
        split_label[extra] = len(V) + np.arange(len(extra))
        split = (np.vstack([V, D[extra]]), [[int(split_label[x]) for x in r] for r in d_iverts])
        claimed = set(dual_at)
        free = np.array([c for c in ids if int(c) not in claimed], dtype=np.int64)
        tree = cKDTree(V[free] / unit) if len(free) else None
        joined = np.flatnonzero(np.bincount(roots, minlength=len(D)) > 1)
        for r in joined:
            members = np.flatnonzero(roots == r)
            held = label[members][label[members] >= 0]
            if len(held):
                label[members] = held[0]
                continue
            if tree is None:
                continue
            reach = [int(free[j]) for j in tree.query_ball_point(D[members[0]] / unit, 1.0, p=np.inf)]
            reach = [c for c in reach if c not in claimed and (np.abs(V[c] - D[members]) < unit).all()]
            if reach:
                c = min(reach, key=lambda c: float(np.hypot(*(V[c] - D[members[0]]))))
                label[members] = c
                claimed.add(c)
    extra = np.flatnonzero(label < 0)
    fresh = {}
    for j in extra:  # ascending, so each group is placed at its lowest-index member
        fresh.setdefault(int(roots[j]), (len(V) + len(fresh), int(j)))
        label[j] = fresh[int(roots[j])][0]
    if fresh:
        V = np.vstack([V, D[[j for _, j in fresh.values()]]])
    rings = _rings_on(label, d_iverts, D)

    new_iverts, rewritten = [], {}
    location = canon[loc]
    for k in range(ncell):
        want = rings[k]
        have = [int(v) for v in iverts[k]]
        if len(have) == len(want) and have[0] in want:
            i = want.index(have[0])
            if want[i:] + want[:i] == have:
                new_iverts.append(iverts[k])
                continue
        if have[0] in want:
            i = want.index(have[0])
            want = want[i:] + want[:i]
        new_iverts.append(want)
        at = [int(location[v]) for v in have]
        if set(at) - set(want):
            reason = "foreign-vertex"
        elif set(want) - set(at):
            reason = "missing-vertex"
        elif len(at) != len(set(at)) or sorted(have) != sorted(want):
            reason = "twin-vertex"  # one location under another id, or twice in the ring
        else:
            reason = "corner-order" if on_boundary[k] else "interior-order"
        rewritten[k] = reason
    return V, new_iverts, rewritten, len(V) - len(verts), unsplit, split


def tessellate(tri, verts, iverts, points) -> tuple[np.ndarray, list, np.ndarray, TessellationReport]:
    """FloPy's ``tri2vor`` output made to tile the domain: ``(verts, iverts, points, report)``.

    ``verts``/``iverts``/``points`` are what ``tri2vor(tri)`` returned. The cells are the
    dual's (:func:`dual_cells`) on FloPy's vertex ids: a cell FloPy got right comes back
    exactly as it wrote it, the rest are rewritten, and a vertex is appended only where
    FloPy has none. When FloPy's cell SET differs from the dual's -- a domain part lost, a
    node outside the domain kept as a cell -- the dual's cells are returned whole,
    renumbered, and a WARNING is logged. Raises :class:`NotConformingDelaunay` for a mesh
    whose dual is not a tiling.
    """

    d_verts, d_iverts, d_points, on_boundary = _dual(tri)
    flopy_points = np.asarray(points, dtype=float)
    ncell = len(iverts)

    def pin(v, rings, xy):
        props = types.SimpleNamespace(get_disv_gridprops=lambda: get_disv_gridprops(v, rings, xcyc=xy))
        return grid_pin(props)

    if flopy_points.shape != d_points.shape or not np.array_equal(flopy_points, d_points):
        lost = len({tuple(p) for p in d_points.tolist()} - {tuple(p) for p in flopy_points.tolist()})
        extra = len({tuple(p) for p in flopy_points.tolist()} - {tuple(p) for p in d_points.tolist()})
        reason = f"FloPy dropped {lost} cell(s) the mesh has and kept {extra} it has not"
        logger.warning(
            "the Voronoi grid's cells were rebuilt from Triangle's mesh (%s): FloPy's %d "
            "cells became %d, renumbered (ledgers 316, 317)", reason, ncell, len(d_iverts))
        # The dual's own vertices, a corner joined across any face the file cannot hold.
        roots = _short_faces(d_verts, d_iverts, written_resolution(d_verts))
        _, first, group = np.unique(roots, return_index=True, return_inverse=True)
        unsplit = len(d_verts) - len(first)
        rings = _rings_on(group, d_iverts, d_verts) if unsplit else d_iverts
        return d_verts[first], rings, d_points, TessellationReport(
            adopted=True, reason=reason, flopy_ncpl=ncell, ncpl=len(d_iverts), unsplit=unsplit,
            split_pin=pin(d_verts, d_iverts, d_points) if unsplit else None)
    new_verts, new_iverts, rewritten, appended, unsplit, split = _onto_flopy(
        d_verts, d_iverts, verts, iverts, on_boundary)
    report = TessellationReport(rewritten=rewritten, appended=appended, unsplit=unsplit,
                                flopy_ncpl=ncell, ncpl=ncell)
    if split is not None:
        split_pin = pin(split[0], split[1], flopy_points)
        if split_pin != pin(new_verts, new_iverts, flopy_points):
            report = dataclasses.replace(report, split_pin=split_pin)
    if report.changed:
        logger.info("Voronoi grid repaired: %s", report.summary())
    return new_verts, new_iverts, flopy_points, report
