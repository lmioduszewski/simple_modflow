"""Turning barrier geometry into the cell-pair faces MODFLOW 6 wants.

A horizontal flow barrier is a thin, low-permeability feature -- a fault gouge, a
slurry cutoff wall, a till finger -- that sits **on the face between two cells**
rather than inside either of them. MODFLOW 6 therefore addresses it as a *pair*
of cells, and FloPy will take that pair and no more: it validates each cellid
against ``idomain`` but never checks that the two are actually neighbours, so a
wrong pair is discovered by MODFLOW 6 aborting mid-run:

    HFB no. 1 is between two unconnected cells: (1,1) and (1,441)

What a modeller actually has is a LINE -- a fault trace digitised in a GeoPackage.
This module is the bridge, and it is the whole reason ``mf.hfb`` is worth more
than the FloPy call it wraps.

**The adjacency used here is exactly MODFLOW 6's.** ``vor.adjacent_cells_idx`` was
checked against the ``ia``/``ja`` arrays read back from a model's own ``.disv.grb``:
1247 pairs against 1247, both difference sets empty. So a pair this module accepts
is a pair MODFLOW 6 will accept.

**A closed wall is not the faces its ring crosses.** A ring passes *through* cells,
so the crossed-face set has gaps at every cell the ring enters and leaves -- on the
canonical 441-cell grid a 540x540 ring crossed 23 faces and left every cell still
hydraulically connected. Sealing a region is a different question: it is the *cut*
between the enclosed cell set and everything else (38 faces there, isolating 28
cells). That is why :func:`enclosed_faces` exists beside :func:`barrier_faces`
rather than as a flag on it.
"""

from __future__ import annotations

from collections import deque
from typing import TYPE_CHECKING, Any

import numpy as np
import shapely as shp
from shapely.errors import GEOSException

from myflopy._logging import get_logger

if TYPE_CHECKING:
    from collections.abc import Sequence

logger = get_logger(__name__)

__all__ = [
    "barrier_faces",
    "enclosed_faces",
    "is_watertight",
    "shared_face_segment",
    "validate_barrier_pairs",
]

#: Faces shorter than this are treated as a corner touch rather than a shared
#: edge. Matches the tolerance ``grid/connectivity.py`` uses when it builds the
#: DISU face list, so the two agree about what "adjacent" means.
FACE_TOLERANCE = 1e-9


def shared_face_segment(vor, cell_a: int, cell_b: int):
    """Return the ``LineString`` two cells share, or ``None`` if they only touch.

    ``grid/connectivity.py`` already computes this while building the DISU face
    list and then keeps only its *length*; the segment itself is what a barrier
    needs, so it is recomputed here rather than left discarded.

    A pair meeting at a single corner returns ``None``: a zero-length face is a
    touch, not a connection, and MODFLOW 6 does not connect such cells.

    Parameters
    ----------
    vor
        A grid exposing ``gdf_vorPolys``.
    cell_a, cell_b
        Zero-based plan cell indices.

    Returns
    -------
    shapely.LineString or None
        The longest shared segment. ``None`` when the cells share no edge.
    """

    geoms = vor.gdf_vorPolys.geometry
    try:
        shared = geoms.iloc[cell_a].boundary.intersection(geoms.iloc[cell_b].boundary)
    except GEOSException as error:
        # An invalid ring from a degenerate mesh; the pair simply has no usable
        # face, which is the same answer as "not adjacent".
        logger.debug("no shared face for cells %d/%d; treated as unconnected: %s", cell_a, cell_b, error)
        return None

    if isinstance(shared, shp.LineString):
        return shared if shared.length > FACE_TOLERANCE else None
    if isinstance(shared, shp.MultiLineString) and not shared.is_empty:
        longest = max(shared.geoms, key=lambda seg: seg.length)
        return longest if longest.length > FACE_TOLERANCE else None
    return None


def _as_geometry(barrier, *, crs=None):
    """Coerce any accepted barrier input to one shapely geometry in the grid's CRS."""

    import geopandas as gpd

    if isinstance(barrier, (gpd.GeoDataFrame, gpd.GeoSeries)):
        if crs is not None and barrier.crs is not None and barrier.crs != crs:
            barrier = barrier.to_crs(crs)
        return barrier.geometry.union_all() if hasattr(barrier, "geometry") else barrier.union_all()
    return barrier


def barrier_faces(
    vor,
    barrier: Any,
    *,
    strict: bool = True,
) -> list[tuple[int, int]]:
    """Return the cell-pair faces a barrier line crosses.

    Parameters
    ----------
    vor
        A grid exposing ``gdf_vorPolys`` and ``adjacent_cells_idx``.
    barrier
        A ``LineString``, ``MultiLineString``, ``GeoSeries`` or ``GeoDataFrame``.
        Reprojected to the grid's CRS when it declares one.
    strict
        Ignore a crossing that only touches a face's END POINT. A trace digitised
        *along* a Voronoi edge otherwise picks up both faces meeting that edge --
        three faces where the modeller drew one.

    Returns
    -------
    list of (int, int)
        Canonically ordered ``(low, high)`` plan-cell pairs, sorted and unique.
    """

    geometry = _as_geometry(barrier, crs=getattr(vor, "crs", None))
    if geometry is None or geometry.is_empty:
        logger.warning("barrier geometry is empty; no faces crossed")
        return []

    # Only cells the line actually touches can own a crossed face, so the
    # neighbour scan starts from them rather than from all ncpl cells.
    touched = vor.gdf_vorPolys.geometry.sindex.query(geometry, predicate="intersects")
    adjacency = vor.adjacent_cells_idx

    faces: set[tuple[int, int]] = set()
    for cell in np.atleast_1d(touched):
        cell = int(cell)
        for neighbour in adjacency[cell]:
            neighbour = int(neighbour)
            pair = (min(cell, neighbour), max(cell, neighbour))
            if pair in faces:
                continue
            segment = shared_face_segment(vor, *pair)
            if segment is None or not geometry.intersects(segment):
                continue
            if strict and not _crosses_interior(geometry, segment):
                continue
            faces.add(pair)

    logger.debug("barrier crosses %d face(s)", len(faces))
    return sorted(faces)


def _crosses_interior(geometry, segment) -> bool:
    """True when ``geometry`` meets ``segment`` somewhere other than its end points.

    A barrier drawn along a cell edge touches the two faces meeting that edge at
    a single shared vertex. Those are not crossings, and counting them turns one
    drawn barrier into three.
    """

    try:
        meeting = geometry.intersection(segment)
    except GEOSException as error:
        logger.debug("could not intersect barrier with a face; face skipped: %s", error)
        return False
    if meeting.is_empty:
        return False
    if meeting.length > FACE_TOLERANCE:
        return True
    ends = shp.MultiPoint([segment.coords[0], segment.coords[-1]])
    return not meeting.within(ends.buffer(FACE_TOLERANCE))


def enclosed_faces(vor, polygon: Any) -> tuple[list[tuple[int, int]], list[int]]:
    """Return the cut that seals the cells inside ``polygon`` from the rest.

    This is deliberately NOT the set of faces the polygon's ring crosses. A ring
    passes *through* cells, so the crossed set has a gap wherever the ring enters
    and leaves one; on the canonical grid a 540x540 ring crossed 23 faces and the
    interior stayed connected to every other cell. The cut -- every face with one
    cell inside and one outside -- is watertight by construction.

    Parameters
    ----------
    vor
        A grid exposing ``gdf_vorPolys`` and ``adjacent_cells_idx``.
    polygon
        The region to seal. A cell is "inside" when its centroid is.

    Returns
    -------
    (faces, interior)
        The cut as canonical ``(low, high)`` pairs, and the enclosed cell indices.
    """

    geometry = _as_geometry(polygon, crs=getattr(vor, "crs", None))
    centroids = vor.gdf_vorPolys.geometry.centroid
    inside = set(np.flatnonzero(shp.contains(geometry, shp.points(np.c_[centroids.x, centroids.y]))).tolist())
    if not inside:
        logger.warning("no cell centroid falls inside the polygon; nothing to enclose")
        return [], []

    adjacency = vor.adjacent_cells_idx
    faces = {
        (min(cell, int(neighbour)), max(cell, int(neighbour)))
        for cell in inside
        for neighbour in adjacency[cell]
        if int(neighbour) not in inside
    }
    logger.debug("enclosing %d cell(s) behind %d face(s)", len(inside), len(faces))
    return sorted(faces), sorted(inside)


def is_watertight(vor, faces: Sequence[tuple[int, int]], seed: int) -> bool:
    """True when ``faces`` disconnect ``seed`` from the rest of the grid.

    A flood fill from ``seed`` that may not cross any listed face. An open wall
    legitimately fails this; a cutoff wall that fails it is leaking, which is
    always a bug.
    """

    blocked = {(min(a, b), max(a, b)) for a, b in faces}
    adjacency = vor.adjacent_cells_idx
    ncpl = len(vor.gdf_vorPolys)

    seen = {seed}
    queue = deque([seed])
    while queue:
        cell = queue.popleft()
        for neighbour in adjacency[cell]:
            neighbour = int(neighbour)
            if (min(cell, neighbour), max(cell, neighbour)) in blocked or neighbour in seen:
                continue
            seen.add(neighbour)
            queue.append(neighbour)
    return len(seen) < ncpl


def validate_barrier_pairs(vor, pairs: Sequence[tuple[Any, Any]]) -> list[tuple[Any, Any]]:
    """Check barrier pairs against the grid, and return them deduplicated.

    Three checks, all of which MODFLOW 6 either rejects loudly or -- worse --
    accepts and gets wrong:

    1. **Duplicates are collapsed.** MODFLOW 6 accepts a repeated face without a
       word and applies its series formula *twice* (measured conductance 14.4727
       -> 0.8418 once -> 0.4335 twice). Worse, ``condsat_reset`` then restores the
       already-modified value, so the face stays wrong for the rest of the run
       even after the barrier is removed.
    2. **A self-pair is rejected.** MODFLOW 6 scans a cell's connections skipping
       the diagonal, so ``(c, c)`` reports as "unconnected".
    3. **Adjacency is required.** Laterally, the two plan cells must share a face;
       vertically, they must be the same plan cell in adjacent layers. Vertical
       barriers on DISV are legal as of MODFLOW 6 **6.7.0** -- FloPy 3.10's
       embedded definition still says otherwise, and is wrong.

    Parameters
    ----------
    vor
        A grid exposing ``adjacent_cells_idx``.
    pairs
        ``((lay, cell), (lay, cell))`` tuples.

    Returns
    -------
    list
        The same pairs, duplicates collapsed, in input order.

    Raises
    ------
    ValueError
        Naming every self-pair and every non-adjacent pair, echoing MODFLOW 6's
        own wording so the message reads the same as the one a run would give.
    """

    adjacency = vor.adjacent_cells_idx
    seen: dict[tuple, tuple] = {}
    self_pairs, unconnected = [], []

    for first, second in pairs:
        (layer_a, cell_a), (layer_b, cell_b) = tuple(first), tuple(second)
        key = tuple(sorted([(int(layer_a), int(cell_a)), (int(layer_b), int(cell_b))]))
        if key in seen:
            continue

        if (layer_a, cell_a) == (layer_b, cell_b):
            self_pairs.append(key)
        elif layer_a == layer_b:
            if int(cell_b) not in {int(n) for n in adjacency[int(cell_a)]}:
                unconnected.append(key)
        elif cell_a != cell_b or abs(int(layer_a) - int(layer_b)) != 1:
            unconnected.append(key)

        seen[key] = (first, second)

    duplicates = len(pairs) - len(seen)
    if duplicates:
        logger.warning(
            "collapsed %d duplicate barrier face(s); MODFLOW 6 would have applied each twice "
            "and left the face permanently wrong",
            duplicates,
        )

    problems = []
    if self_pairs:
        problems.append(
            "a barrier cannot sit between a cell and itself: "
            + ", ".join(f"({a[0] + 1},{a[1] + 1})" for a, _ in self_pairs[:5])
        )
    if unconnected:
        shown = ", ".join(
            f"({a[0] + 1},{a[1] + 1}) and ({b[0] + 1},{b[1] + 1})" for a, b in unconnected[:5]
        )
        more = f" and {len(unconnected) - 5} more" if len(unconnected) > 5 else ""
        problems.append(
            f"{len(unconnected)} barrier(s) are between two unconnected cells: {shown}{more}. "
            "MODFLOW 6 aborts on these; laterally the cells must share a face, vertically they "
            "must be the same cell in adjacent layers"
        )
    if problems:
        raise ValueError("; ".join(problems))

    return list(seen.values())
