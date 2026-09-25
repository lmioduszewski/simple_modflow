"""MODFLOW 6's horizontal-connection geometry for a DISV grid, in ONE place.

A DISV file holds no connection lengths or face widths: MF6 computes them itself,
per cell pair, from the vertices and the cell centres (``DisvGeom.f90``:
``shared_edge``, ``distance``, ``distance_normal``). Everything myflopy writes that
MF6 reads back as such a number -- a DISU's ``CL12``/``HWVA``, a LAK wall's
``CONNLEN``/``CONNWIDTH`` -- has to be that same computation, or it describes a
different aquifer from the one MF6 builds. :class:`DisvFaces` is it, and
:func:`build_disu_connectivity` is a view of it (ledger 328).
"""

from __future__ import annotations

from collections.abc import Mapping
from itertools import chain

import numpy as np

from myflopy._logging import get_logger

logger = get_logger(__name__)


def _distance(x1, y1, x2, y2):
    """MF6's ``distance``: the length of the segment, written as MF6 writes it."""

    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def _distance_normal(x0, y0, x1, y1, x2, y2):
    """MF6's ``distance_normal``: from ``(x0, y0)`` perpendicular to the LINE through the edge.

    The expression is MF6's, term for term, so the numbers agree to the last bit
    with what MF6 computes from the same coordinates. Arrays or scalars.
    """

    return np.abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)) / _distance(x1, y1, x2, y2)


class DisvFaces:
    """MODFLOW 6's own horizontal-connection geometry for a DISV grid, per cell pair.

    What MF6 means by a connection length is the distance from ONE cell's centre to
    the shared face, per side -- DISU's ``CL12`` is "the center of cell n [to] the
    shared face with each adjacent m cell", LAK's ``CONNLEN`` "the connected GWF
    CELLID node [to] the lake". For a DISV grid MF6 computes it itself: the
    PERPENDICULAR distance from the cell's DISV centre (``cell2d``'s ``xc``/``yc``,
    which for a Voronoi grid is its generator point, not its polygon centroid) to the
    LINE through the two vertices the cells share, and the face width is the
    distance between those vertices (``DisvGeom.f90``: ``shared_edge``,
    ``distance_normal``). This class is that computation, and the ONE place myflopy
    makes it: ``grid.cl12``/``grid.hwva`` (:func:`build_disu_connectivity`) and every
    LAK wall's ``CONNLEN``/``CONNWIDTH`` read it, so a DISU written from these arrays
    reproduces the DISV it came from (ledger 328). (``grid.calculate_distance``
    measures something else -- a polygon centroid to the nearest point of the shared
    boundary -- and nothing that writes a model reads it.)

    Which cells connect is MF6's rule too: two cells connect when their vertex rings
    share an EDGE -- two consecutive vertex NUMBERS -- never a corner alone, and never
    by coordinates, so two vertices at one point with different numbers do not join
    their cells (in MF6 or here). On every mesh measured the pairs equal the pairs
    whose polygons share a boundary segment: two real models (18,107 and 30,444 cells)
    and five generated meshes, a reflex-cornered one included.

    Built from :meth:`VoronoiGridPlus.get_disv_gridprops` -- the arrays the DISV is
    written from -- so the numbers are the ones MF6 will compute from the file. The
    table of every pair is built once, on first use.

    Parameters
    ----------
    gridprops : dict
        ``get_disv_gridprops()``'s result: ``"vertices"`` rows ``[iv, x, y]`` and
        ``"cell2d"`` rows ``[icell, xc, yc, ncvert, iv1, iv2, ...]``, where ``icell``
        counts 0, 1, 2, ... in row order (MF6's own rule: it stops on anything
        else). A ring may repeat its first vertex at the end; it is read once.
        Other keys are ignored.

    Raises
    ------
    KeyError
        ``gridprops`` has no ``"vertices"`` or no ``"cell2d"``.
    ValueError
        A ``cell2d`` row's ``icell`` is not its row number (the message names the
        first one), or a number in a row cannot be read by ``int()``/``float()``.

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
    >>> faces.length_and_width(0, 1, toward=0)   # and cell 0's is 0.5 from it
    (0.5, 1.0)
    """

    def __init__(self, gridprops) -> None:
        vertices = gridprops["vertices"]
        cell2d = gridprops["cell2d"]
        self._vertex_xy = {int(row[0]): (float(row[1]), float(row[2])) for row in vertices}
        centres, rings = [], []
        for position, row in enumerate(cell2d):
            row = list(row)
            if int(row[0]) != position:
                raise ValueError(
                    f"DisvFaces: cell2d row {position} is numbered {int(row[0])}; MF6 reads "
                    "cells in order (0, 1, 2, ...) and stops on anything else"
                )
            count = int(row[3])
            centres.append((float(row[1]), float(row[2])))
            ring = [int(vertex) for vertex in row[4 : 4 + count]]
            if len(ring) > 1 and ring[0] == ring[-1]:
                ring = ring[:-1]
            rings.append(ring)
        self._centres = np.asarray(centres, dtype=np.float64).reshape(-1, 2)
        self._rings = rings
        self._table: dict[str, np.ndarray] | None = None
        self._rows: dict[tuple[int, int], int] | None = None

    @property
    def ncpl(self) -> int:
        """The number of cells in one layer (``cell2d`` rows).

        Examples
        --------
        >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 1.0, 1.0]],
        ...          "cell2d": [[0, 0.6, 0.3, 3, 0, 2, 1]]}
        >>> DisvFaces(props).ncpl
        1
        """

        return len(self._rings)

    # -- the table: one row per connected pair, built once ------------------------------

    def _edges(self):
        """Every ring edge as flat arrays: owning cell, position in its ring, and its two vertices."""

        lengths = np.fromiter((len(ring) for ring in self._rings), dtype=np.int64, count=self.ncpl)
        total = int(lengths.sum())
        first = np.fromiter(chain.from_iterable(self._rings), dtype=np.int64, count=total)
        cell = np.repeat(np.arange(self.ncpl, dtype=np.int64), lengths)
        start = np.repeat(np.cumsum(lengths) - lengths, lengths)
        position = np.arange(total, dtype=np.int64) - start
        following = np.arange(total, dtype=np.int64) + 1
        ends = position == np.repeat(lengths - 1, lengths)
        following[ends] = start[ends]
        second = first[following]
        keep = first != second  # a ring repeating a vertex has no edge there
        return cell[keep], position[keep], first[keep], second[keep]

    def _build(self) -> dict[str, np.ndarray]:
        """The pair table: MF6's ``vertexconnect`` + ``cprops``, vectorized."""

        if self._table is not None:
            return self._table
        cell, position, first, second = self._edges()
        known = np.fromiter(self._vertex_xy, dtype=np.int64, count=len(self._vertex_xy))
        size = int(max(known.max(initial=-1), first.max(initial=-1)) + 1)
        held = np.zeros(size, dtype=bool)
        held[known] = True
        xy = np.full((size, 2), np.nan)
        for number, point in self._vertex_xy.items():
            xy[number] = point
        missing = np.flatnonzero(~held[first])
        if missing.size:
            raise KeyError(int(first[missing[0]]))

        low = np.minimum(first, second)
        key = low * size + np.maximum(first, second)
        order = np.lexsort((position, cell, key))
        cell, position, first, second, key = (
            cell[order], position[order], first[order], second[order], key[order]
        )
        # One candidate per two records of the same edge in DIFFERENT cells. Almost every
        # edge has one or two owners; three or more (a degenerate mesh) pair them all,
        # as MF6's vertexconnect would.
        bounds = np.flatnonzero(np.diff(key)) + 1
        starts = np.concatenate([[0], bounds])
        stops = np.concatenate([bounds, [key.size]])
        sizes = stops - starts
        a_index = list(starts[sizes == 2])
        b_index = list(starts[sizes == 2] + 1)
        for begin, end in zip(starts[sizes > 2], stops[sizes > 2], strict=True):
            for i in range(begin, end):
                for j in range(i + 1, end):
                    a_index.append(i)
                    b_index.append(j)
        a_index = np.asarray(a_index, dtype=np.int64)
        b_index = np.asarray(b_index, dtype=np.int64)
        distinct = cell[a_index] != cell[b_index]
        a_index, b_index = a_index[distinct], b_index[distinct]
        swap = cell[a_index] > cell[b_index]
        lo_index = np.where(swap, b_index, a_index)
        hi_index = np.where(swap, a_index, b_index)
        # MF6 builds each connection ONCE, from the lower node (disvconnections skips
        # m < n), and takes the first shared edge in THAT cell's ring order.
        choose = np.lexsort((position[lo_index], cell[hi_index], cell[lo_index]))
        lo_index, hi_index = lo_index[choose], hi_index[choose]
        pair = cell[lo_index] * self.ncpl + cell[hi_index]
        first_of_pair = np.concatenate([[True], pair[1:] != pair[:-1]]) if pair.size else pair.astype(bool)
        lo_index, hi_index = lo_index[first_of_pair], hi_index[first_of_pair]

        lo, hi = cell[lo_index], cell[hi_index]
        v1, v2 = first[lo_index], second[lo_index]  # the edge in the LOWER cell's ring order
        x1, y1 = xy[v1, 0], xy[v1, 1]
        x2, y2 = xy[v2, 0], xy[v2, 1]
        width = _distance(x1, y1, x2, y2)
        with np.errstate(divide="ignore", invalid="ignore"):
            length_lo = _distance_normal(self._centres[lo, 0], self._centres[lo, 1], x1, y1, x2, y2)
            length_hi = _distance_normal(self._centres[hi, 0], self._centres[hi, 1], x1, y1, x2, y2)
        self._table = {
            "lo": lo, "hi": hi,
            "lo_edge": np.column_stack([v1, v2]),
            "hi_edge": np.column_stack([first[hi_index], second[hi_index]]),
            "width": width, "length_lo": length_lo, "length_hi": length_hi,
            "xy": xy,
        }
        return self._table

    def _row(self, a, b) -> int | None:
        """The table row of the pair ``{a, b}``, or ``None`` when they share no edge."""

        a, b = int(a), int(b)
        for cell in (a, b):
            if not 0 <= cell < self.ncpl:
                raise KeyError(cell)
        if self._rows is None:
            table = self._build()
            self._rows = {
                (int(lo), int(hi)): row
                for row, (lo, hi) in enumerate(zip(table["lo"], table["hi"], strict=True))
            }
        return self._rows.get((min(a, b), max(a, b)))

    # -- the public reads -------------------------------------------------------------

    def shared_edge(self, a: int, b: int) -> tuple[int, int] | None:
        """The two vertices cells ``a`` and ``b`` share as an EDGE, or ``None``.

        MF6's rule: of the edges the two cells share, the FIRST in the LOWER-numbered
        cell's ring order -- MF6 builds each connection once, from the lower node --
        so which edge is taken does not depend on the argument order. Only a mesh
        where two cells share more than one edge (a vertex in the middle of a straight
        face) has a choice to make; none of the meshes measured for ledger 328 did. A
        pair that appears in the SAME order in both rings (a mesh whose rings are not
        consistently oriented) is accepted too -- MF6 would find no face there and
        connect nothing, which is a mesh defect this cannot repair and should not
        hide by returning ``None`` for a face that exists.

        Parameters
        ----------
        a : int or numpy.integer
            One cell (a zero-based ``cell2d`` number). The edge is returned in ITS
            ring order.
        b : int or numpy.integer
            The other cell.

        Returns
        -------
        tuple of (int, int), or None
            The two vertex numbers, in ``a``'s ring order; ``None`` when the cells
            share no edge (a corner or nothing).

        Raises
        ------
        KeyError
            ``a`` or ``b`` is not a cell of ``cell2d``, or a ring names a vertex
            ``vertices`` does not hold (the message is the number).
        ValueError
            ``a`` or ``b`` cannot be read by ``int()``.

        Examples
        --------
        >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
        ...                       [4, 2.0, 0.0], [5, 2.0, 1.0]],
        ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 2, 3, 1], [1, 1.5, 0.5, 4, 1, 3, 5, 4]]}
        >>> DisvFaces(props).shared_edge(0, 1)
        (3, 1)
        >>> DisvFaces(props).shared_edge(1, 0)
        (1, 3)
        """

        row = self._row(a, b)
        if row is None:
            return None
        table = self._build()
        edge = table["lo_edge"][row] if int(a) == int(table["lo"][row]) else table["hi_edge"][row]
        return int(edge[0]), int(edge[1])

    def length_and_width(self, a: int, b: int, *, toward: int) -> tuple[float, float] | None:
        """``(CONNLEN, CONNWIDTH)`` of the face ``a|b`` measured from cell ``toward``.

        Parameters
        ----------
        a : int or numpy.integer
            One of the two cells sharing the face (a zero-based ``cell2d`` number).
        b : int or numpy.integer
            The other one.
        toward : int or numpy.integer
            Whose DISV centre the length is measured from -- for a LAK wall, the
            aquifer cell the connection is written into; for a DISU row ``n``, cell
            ``n``. Normally ``a`` or ``b``; any cell is accepted and not checked
            against the pair.

        Returns
        -------
        tuple of float, or None
            The perpendicular distance from ``toward``'s centre to the line of the
            shared edge (model length units), and the edge's length; ``None`` when
            the cells share no edge (only a corner) or the edge has no length.

        Raises
        ------
        KeyError
            ``a``, ``b`` or ``toward`` is not a cell of ``cell2d``, or a ring names a
            vertex ``vertices`` does not hold (the message is the number).
        ValueError
            A cell number ``int()`` cannot read.

        Examples
        --------
        >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
        ...                       [4, 2.0, 0.0], [5, 2.0, 1.0]],
        ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 2, 3, 1], [1, 1.25, 0.5, 4, 1, 3, 5, 4]]}
        >>> DisvFaces(props).length_and_width(0, 1, toward=0)
        (0.5, 1.0)
        """

        toward = int(toward)
        if not 0 <= toward < self.ncpl:
            raise KeyError(toward)
        row = self._row(a, b)
        if row is None:
            return None
        table = self._build()
        width = float(table["width"][row])
        if width == 0.0:
            return None
        if toward == int(table["lo"][row]):
            length = table["length_lo"][row]
        elif toward == int(table["hi"][row]):
            length = table["length_hi"][row]
        else:
            (v1, v2), xy = table["lo_edge"][row], table["xy"]
            x0, y0 = self._centres[toward]
            length = _distance_normal(x0, y0, xy[v1, 0], xy[v1, 1], xy[v2, 0], xy[v2, 1])
        return float(length), width

    def disu_arrays(self, *, tol: float = 0.0):
        """The DISU connection table MF6 would build from this DISV: ``iac, ja, cl12, hwva, nja``.

        Row ``n`` holds cell ``n`` itself first (``cl12 = hwva = 0``), then every cell
        it shares an edge with, in increasing order; ``cl12`` at ``(n, m)`` is the
        distance from ``n``'s centre to the face (so ``cl12`` is NOT symmetric) and
        ``hwva`` the face width (it is).

        Parameters
        ----------
        tol : float, default 0.0
            A face whose width is ``<= tol`` (model length units) is not connected.
            ``0.0`` drops only zero-width faces, which MF6 would divide by.

        Returns
        -------
        tuple
            ``iac`` (int64, ``(ncpl,)``), ``ja`` (int64, ``(nja,)``, zero-based),
            ``cl12`` (float64, ``(nja,)``), ``hwva`` (float64, ``(nja,)``), and ``nja``
            (int).

        Raises
        ------
        KeyError
            A ring names a vertex ``vertices`` does not hold.
        TypeError
            ``tol`` is not a number.

        Examples
        --------
        >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
        ...                       [4, 2.0, 0.0], [5, 2.0, 1.0]],
        ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 2, 3, 1], [1, 1.25, 0.5, 4, 1, 3, 5, 4]]}
        >>> iac, ja, cl12, hwva, nja = DisvFaces(props).disu_arrays()
        >>> iac.tolist(), ja.tolist(), cl12.tolist(), hwva.tolist(), nja
        ([2, 2], [0, 1, 1, 0], [0.0, 0.5, 0.0, 0.25], [0.0, 1.0, 0.0, 1.0], 4)
        """

        if isinstance(tol, bool) or not isinstance(tol, (int, float, np.integer, np.floating)):
            raise TypeError(f"tol must be a number, got {type(tol).__name__}")
        table = self._build()
        keep = table["width"] > float(tol)
        lo, hi = table["lo"][keep], table["hi"][keep]
        width = table["width"][keep]
        row = np.concatenate([lo, hi])
        other = np.concatenate([hi, lo])
        length = np.concatenate([table["length_lo"][keep], table["length_hi"][keep]])
        widths = np.concatenate([width, width])
        order = np.lexsort((other, row))
        counts = np.bincount(row, minlength=self.ncpl).astype(np.int64)
        iac = counts + 1
        nja = int(iac.sum())
        ia = np.concatenate([[0], np.cumsum(iac)])
        diagonal = np.zeros(nja, dtype=bool)
        diagonal[ia[:-1]] = True
        ja = np.empty(nja, dtype=np.int64)
        cl12 = np.zeros(nja, dtype=np.float64)
        hwva = np.zeros(nja, dtype=np.float64)
        ja[diagonal] = np.arange(self.ncpl, dtype=np.int64)
        ja[~diagonal] = other[order]
        cl12[~diagonal] = length[order]
        hwva[~diagonal] = widths[order]
        return iac, ja, cl12, hwva, nja


def build_disu_connectivity(gridprops, *, tol: float = 0.0, validate: bool = True):
    """DISU connectivity (``iac, ja, cl12, hwva, nja``) with MF6's own DISV geometry.

    The numbers MF6 computes for itself when it reads the same grid as a DISV
    (:class:`DisvFaces`), so a DISU written from them reproduces that DISV -- heads
    and ``FLOW-JA-FACE`` to rounding. ``cl12`` is per SIDE, as MF6 defines it: at
    ``(n, m)``, the distance from cell ``n``'s DISV centre to the shared face, NOT the
    centre-to-centre spacing it held from 2025-07-02 until ledger 328 (which made a
    DISU's horizontal conductance about half the DISV's). The polygon centroid is not
    the centre: MF6 measures from ``cell2d``'s ``xc``/``yc``.

    Parameters
    ----------
    gridprops : dict
        ``get_disv_gridprops()``'s result (``"vertices"`` and ``"cell2d"``), as for
        :class:`DisvFaces`. Any other mapping with those two keys works; anything that
        is not a mapping -- a ``GeoDataFrame`` of the cell polygons, say -- is refused
        by name (``TypeError``): polygons carry no DISV centres.
    tol : float, default 0.0
        A face whose width is ``<= tol`` (model length units) is not connected.
    validate : bool, default True
        Check the result before returning it: each row starts with its own cell
        and zero ``cl12``/``hwva``, every connection appears in both rows, ``hwva``
        agrees between them, and every ``cl12`` is finite and non-negative with at
        least one side positive (MF6 divides by ``cl1 + cl2``). ``False`` skips the
        checks; the arrays are the same.

    Returns
    -------
    tuple
        ``iac`` (int64, ``(ncpl,)``), ``ja`` (int64, ``(nja,)``, zero-based), ``cl12``
        (float64), ``hwva`` (float64), ``nja`` (int) -- see
        :meth:`DisvFaces.disu_arrays`.

    Raises
    ------
    TypeError
        ``gridprops`` is not a mapping, or ``tol`` is not a number. The two-argument
        ``(polygons, adjacency)`` call this function took before ledger 328 fails too,
        but with Python's own arity message ("takes 1 positional argument but 2 were
        given"), before the body runs.

    See Also
    --------
    DisvFaces : The rule, one face at a time.
    myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus.get_disu_connectivity : The same, cached on the grid.
    KeyError
        ``gridprops`` has no ``"vertices"``/``"cell2d"``, or a ring names a vertex
        ``vertices`` does not hold.
    ValueError
        ``cell2d`` is not numbered 0, 1, 2, ... in order, or ``validate`` found a
        defect (the message names the cells).

    Examples
    --------
    >>> from myflopy.modflow.mf6.grid.connectivity import build_disu_connectivity
    >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
    ...                       [4, 2.0, 0.0], [5, 2.0, 1.0]],
    ...          "cell2d": [[0, 0.5, 0.5, 4, 0, 2, 3, 1], [1, 1.5, 0.5, 4, 1, 3, 5, 4]]}
    >>> iac, ja, cl12, hwva, nja = build_disu_connectivity(props)
    >>> cl12.tolist()   # half a cell each side, not the 1.0 between the centres
    [0.0, 0.5, 0.0, 0.5]
    """

    if not isinstance(gridprops, Mapping):
        raise TypeError(
            "build_disu_connectivity takes get_disv_gridprops()'s dict (ledger 328): "
            f"got {type(gridprops).__name__}. The (polygons, adjacency) form measured "
            "from polygon centroids, which are not the centres MF6 uses."
        )
    faces = DisvFaces(gridprops)
    logger.info('building connectivity (iac, ja, cl12, hwva, nja) for %d cells', faces.ncpl)
    iac, ja, cl12, hwva, nja = faces.disu_arrays(tol=tol)
    if validate:
        _validate(iac, ja, cl12, hwva, nja)
    return iac, ja, cl12, hwva, nja


def _validate(iac, ja, cl12, hwva, nja) -> None:
    """Refuse a connection table MF6 could not use, naming the first cells at fault."""

    ncpl = iac.size
    ia = np.concatenate([[0], np.cumsum(iac)])
    if nja != ja.size or not (ja.size == cl12.size == hwva.size):
        raise ValueError(f"nja ({nja}) and the array sizes ({ja.size}, {cl12.size}, {hwva.size}) differ")
    if (iac < 1).any():
        raise ValueError(f"iac < 1 for cells {np.flatnonzero(iac < 1)[:10].tolist()}")
    heads = ia[:-1]
    bad = np.flatnonzero(ja[heads] != np.arange(ncpl))
    if bad.size:
        raise ValueError(f"rows {bad[:10].tolist()} do not start with their own cell")
    if (cl12[heads] != 0.0).any() or (hwva[heads] != 0.0).any():
        raise ValueError("a row's own entry must carry cl12 = hwva = 0")
    row = np.repeat(np.arange(ncpl), iac)
    off = np.ones(nja, dtype=bool)
    off[heads] = False
    n, m = row[off], ja[off]
    forward = n * ncpl + m
    backward = m * ncpl + n
    order = np.argsort(forward)
    where = np.searchsorted(forward[order], backward)
    where = np.minimum(where, forward.size - 1)
    partner = order[where]
    unmatched = forward[partner] != backward
    if unmatched.any():
        first = np.flatnonzero(unmatched)[:10]
        raise ValueError(
            f"Asymmetric neighbor pairs detected (count={int(unmatched.sum())}), "
            f"e.g. {list(zip(n[first].tolist(), m[first].tolist(), strict=True))}"
        )
    width, other_width = hwva[off], hwva[off][partner]
    if (np.abs(width - other_width) > 1e-9 * np.maximum(1.0, width)).any():
        raise ValueError("hwva differs between the two rows of a connection")
    length, other_length = cl12[off], cl12[off][partner]
    bad = ~np.isfinite(length) | (length < 0.0) | (length + other_length <= 0.0)
    if bad.any():
        first = np.flatnonzero(bad)[:10]
        raise ValueError(
            "cl12 must be finite, non-negative and positive on one side at least: "
            f"cells {list(zip(n[first].tolist(), m[first].tolist(), strict=True))}"
        )
