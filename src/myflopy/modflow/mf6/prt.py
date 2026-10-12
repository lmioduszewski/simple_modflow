"""Canonical MODFLOW 6 PRT workflow.

PRT is a separate MF6 simulation that consumes a completed groundwater-flow
model's head and budget outputs through FMI. The resulting track CSV is exposed
as view nouns on :class:`PRTRunResults` (``pathlines``/``travel_time``/
``endpoints``/``capture``) and feeds the same plotting and 3D visualization layer
as other particle-tracking engines.
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any

import flopy
import geopandas as gpd
import pandas as pd
import shapely as shp

from myflopy._length_units import mf6_length_mirror
from myflopy._logging import get_logger
from myflopy._time_units import check_tdis_time_units
from myflopy.modflow.mf6.grid.domain import model_domain
from myflopy.modflow.mf6.grid.geometry import (
    cell_shares,
    domain_active,
    miss_reason,
    where_they_are,
)
from myflopy.modflow.mf6.hydraulics import declared_time_unit

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

logger = get_logger(__name__)


def _mfdata_value(value, default=None):
    """Unwrap a FloPy MFData value (calling ``get_data()`` if present), or return ``default``."""

    if value is None:
        return default
    getter = getattr(value, "get_data", None)
    return getter() if callable(getter) else value


def _default_prt_name(flow_name: str) -> str:
    """A default PRT model name derived from the flow model name (sanitized, ``_prt`` suffix)."""

    stem = re.sub(r"[^A-Za-z0-9_]", "_", str(flow_name))
    return f"{stem[:12]}_prt"


def _output_filename(model: SimulationBase, kind: str) -> str:
    """The GWF head/budget output filename from its OC record (falling back to ``<name>.hds/.cbc``)."""

    oc = getattr(model.gwf, "oc", None)
    record = getattr(oc, f"{kind}_filerecord", None) if oc is not None else None
    data = _mfdata_value(record)
    if data is not None and len(data):
        value = data[0][0]
        return value.decode() if isinstance(value, bytes) else str(value)
    suffix = "hds" if kind == "head" else "cbc"
    return f"{model.name}.{suffix}"


def _grid_filename(model: SimulationBase) -> str:
    """Return the GWF binary-grid filename required by PRT local-z releases."""

    grid = getattr(model.gwf, "dis", None) or getattr(model.gwf, "disv", None) or getattr(model.gwf, "disu", None)
    record = getattr(grid, "grb_filerecord", None) if grid is not None else None
    data = _mfdata_value(record)
    if data is not None and len(data):
        value = data[0][0]
        return value.decode() if isinstance(value, bytes) else str(value)
    package_type = getattr(grid, "package_type", "disv")
    return f"{model.name}.{package_type}.grb"


def _tdis_perioddata(model: SimulationBase) -> list[tuple[float, int, float]]:
    """The flow model's TDIS ``(perlen, nstp, tsmult)`` rows (defaulting to one steady period)."""

    tdis = getattr(model.sim, "tdis", None)
    data = _mfdata_value(getattr(tdis, "perioddata", None))
    if data is None:
        return [(1.0, 1, 1.0)]
    return [
        (float(row["perlen"]), int(row["nstp"]), float(row["tsmult"]))
        for row in data
    ]


def _copy_grid_to_prt(flow_model: SimulationBase, prt_model):
    """Replicate the GWF DIS/DISV discretization onto a PRT model (raises for unsupported grids).

    The ``length_units`` word goes through :func:`~myflopy._length_units.mf6_length_mirror`:
    a short spelling a loaded run's DIS holds (``"ft"``) becomes MF6's word, and a
    word MF6 read as undefined stays as it was (ledger 456).
    """

    gwf = flow_model.gwf
    if getattr(gwf, "disv", None) is not None:
        disv = gwf.disv
        return flopy.mf6.ModflowPrtdisv(
            prt_model,
            length_units=mf6_length_mirror(_mfdata_value(disv.length_units)),
            xorigin=_mfdata_value(disv.xorigin),
            yorigin=_mfdata_value(disv.yorigin),
            angrot=_mfdata_value(disv.angrot),
            nlay=int(_mfdata_value(disv.nlay)),
            ncpl=int(_mfdata_value(disv.ncpl)),
            nvert=int(_mfdata_value(disv.nvert)),
            top=_mfdata_value(disv.top),
            botm=_mfdata_value(disv.botm),
            idomain=_mfdata_value(disv.idomain),
            vertices=_mfdata_value(disv.vertices),
            cell2d=_mfdata_value(disv.cell2d),
        )
    if getattr(gwf, "dis", None) is not None:
        dis = gwf.dis
        return flopy.mf6.ModflowPrtdis(
            prt_model,
            length_units=mf6_length_mirror(_mfdata_value(dis.length_units)),
            xorigin=_mfdata_value(dis.xorigin),
            yorigin=_mfdata_value(dis.yorigin),
            angrot=_mfdata_value(dis.angrot),
            nlay=int(_mfdata_value(dis.nlay)),
            nrow=int(_mfdata_value(dis.nrow)),
            ncol=int(_mfdata_value(dis.ncol)),
            delr=_mfdata_value(dis.delr),
            delc=_mfdata_value(dis.delc),
            top=_mfdata_value(dis.top),
            botm=_mfdata_value(dis.botm),
            idomain=_mfdata_value(dis.idomain),
        )
    raise NotImplementedError("Canonical PRT currently supports GWF DIS and DISV grids.")


def _cell_frame(model: SimulationBase):
    """The model's cell polygons indexed by cell id, or ``None`` for another grid kind.

    ``None`` wherever the shared coverage rule cannot speak the model's cell ids:
    a structured (DIS) grid, whose cell is ``(row, col)``, and a model carrying no
    ``vor``, which has no cell frame at all. Those keep FloPy's own ``intersect``.

    The frame is the ``vor``'s, and it is the WORLD frame in every path myflopy
    builds one: a USG import's comes from the .gsf in world coordinates while its
    DISV declares the origin the mesh was shifted by, and a reopened run's is
    rebuilt from ``modelgrid.verts``, which has the offset applied. The model's
    OWN frame is reached from there by ``get_local_coords`` (ledger 218) -- a
    different conversion, and still the one a PRP record is written with.
    """

    if getattr(model.gwf.modelgrid, "grid_type", None) != "vertex":
        return None
    return getattr(getattr(model, "vor", None), "gdf_vorPolys", None)


def _release_frame(points) -> gpd.GeoSeries:
    """The release locations as a GeoSeries, whatever form they were given in.

    Point geometries, a GeoDataFrame or GeoSeries of them, or bare ``(x, y)``
    pairs -- which carry no CRS, and are read as the grid's exactly as they always
    have been. One frame, so the miss message can report where the points are.
    """

    geometries = list(getattr(points, "geometry", points))
    coords = [
        (float(point.x), float(point.y))
        if hasattr(point, "x")
        else (float(point[0]), float(point[1]))
        for point in geometries
    ]
    return gpd.GeoSeries(
        [shp.Point(x, y) for x, y in coords], crs=getattr(points, "crs", None)
    )


def _in_grid_crs(given: gpd.GeoSeries, polys) -> gpd.GeoSeries:
    """The release points in the grid's own CRS, which is what cells are resolved in.

    Every other myflopy source reprojects to it (the GeoPackage source, SFR, the
    observation targets); without it a frame read in degrees over a State Plane
    grid lands on no cell at all, and a reader of the refusal would be told the
    points are somewhere they are not. This is the WORLD frame either way -- the
    model's own frame is a separate conversion, applied to every record below.
    """

    grid_crs = getattr(polys, "crs", None)
    if given.crs is None or grid_crs is None or given.crs == grid_crs:
        return given
    return given.to_crs(grid_crs)


def _flopy_cellid(grid, point, layer: int) -> tuple[int, ...]:
    """The cellid FloPy's own ``intersect`` puts a point in (structured grids only).

    Knows nothing about idomain and answers with the lowest CELL2D on a shared
    face, which is the defect ledger 242 fixed for DISV. It is kept for the grids
    the shared rule cannot index, where it is the only answer available.
    """

    cell = grid.intersect(float(point.x), float(point.y))
    if isinstance(cell, tuple):
        return (
            (int(layer), *(int(value) for value in cell))
            if len(cell) == 2
            else tuple(int(value) for value in cell)
        )
    return (int(layer), int(cell))


def _refuse_or_report_missed(polys, given: gpd.GeoSeries, missed, *, layer: int) -> None:
    """Ledger 211 for release points: raise when none land, warn by name when some do.

    Neither miss used to be reported HERE, and both were reported badly later. A
    point resolved into an inactive cell was written, and FloPy's PRP refused the
    package from ``_check_valid_cellids`` -- ``Cellid (0, 4) is outside of the
    active model grid``, raised while ``PRTProject`` is built, naming a cellid the
    caller never wrote and neither the point that produced it nor the active
    neighbour the point also touches (measured: that neighbour releases and
    tracks). A point off the grid raised FloPy's ``point given is outside of the
    model area``, which said neither where the grid is nor which CRS either side is
    in, and took the other ninety-nine points of a well file with it. A point that
    lands nowhere now leaves no record: a PARTIAL miss can be deliberate -- one
    well file serving several grids -- so it warns, naming the points by their
    position in the input and the coordinates they were GIVEN as. Nothing landing
    means the points and the model are in different places, which is what the two
    extents and the two CRSs in the message are read for.

    ``missed`` pairs each point's index with its LOCATED geometry, since that is
    the frame ``miss_reason`` measures against the grid in.
    """

    target = f"active cell in layer {layer}"
    reasons = [
        f"#{index} ({given.iloc[index].x:,.1f}, {given.iloc[index].y:,.1f}) "
        f"({miss_reason(polys, point, target=target)})"
        for index, point in missed
    ]
    listing = ", ".join(reasons[:10])
    if len(reasons) > 10:
        listing += f", and {len(reasons) - 10} more"
    if len(missed) == len(given):
        raise ValueError(
            f"No release point lands on an {target}, so the PRP package would "
            f"release no particles: {listing}. "
            + where_they_are(given, polys, noun="release points", source_crs=given.crs)
        )
    logger.warning(
        "%d of %d release points land on no %s and are left out: %s",
        len(missed),
        len(given),
        target,
        listing,
    )


@dataclass(frozen=True)
class PRTReleasePoints:
    """The set of particle starting locations for a PRT (particle-tracking) run.

    Holds the PRP (particle release point) ``packagedata`` -- one record per
    released particle, giving its id and starting cell/coordinates -- in the form
    MF6's PRT model expects. Build it with the :meth:`from_cells` constructor
    (release particles at given cells) or :meth:`from_points` (release them at
    given coordinates) rather than assembling records by hand, combine sets with
    :meth:`merge`, then hand it to :class:`PRTProject` / ``ParticleTracking.prt(...)``.

    Attributes
    ----------
    packagedata : tuple of tuple
        Normalized PRP release records, one tuple per particle. Records are
        5-tuples ``(irpt, cellid, x, y, local_z)``, or 6-tuples ending in a
        release-group boundname when ``group=`` was given. ``irpt`` is ZERO-based
        (FloPy adds one when it writes the PRP file); ``cellid`` is ``(layer, cell)``
        on a DISV grid or ``(layer, row, col)`` on a DIS grid, zero-based; ``x``/``y``
        are in the model's OWN frame (what a PRP record is written in, which on an
        offset grid differs from the world frame); the fifth value is MF6's ``zrpt``,
        read as a local z within the cell when :class:`PRTProject` is built with
        ``local_z=True`` (its default) and as an elevation otherwise.
    """

    packagedata: tuple[tuple[Any, ...], ...]

    @property
    def has_groups(self) -> bool:
        """Whether any record carries a release-group boundname (a 6-tuple)."""

        return any(len(row) >= 6 for row in self.packagedata)

    @property
    def groups(self) -> tuple[str, ...]:
        """The distinct release-group labels, in first-seen order (empty when ungrouped).

        These are the labels MF6 echoes -- **uppercased** -- into the track CSV's
        ``name`` column, which is what ``results.capture`` groups by.
        """

        seen: list[str] = []
        for row in self.packagedata:
            if len(row) >= 6 and row[5] not in seen:
                seen.append(str(row[5]))
        return tuple(seen)

    @staticmethod
    def _resolve_groups(group, count: int) -> list[str] | None:
        """Broadcast ``group`` to one label per release point (``None`` stays ungrouped)."""

        if group is None:
            return None
        if isinstance(group, str):
            return [group] * count
        labels = [str(value) for value in group]
        if len(labels) != count:
            raise ValueError(
                f"group must be one label or one per release point; got {len(labels)} "
                f"labels for {count} points."
            )
        return labels

    @classmethod
    def merge(cls, *sets: PRTReleasePoints) -> PRTReleasePoints:
        """Combine release-point sets into one PRP, renumbering ``irpt`` sequentially.

        The practical way to build a multi-group release: make one grouped set per
        capture zone with :meth:`from_cells`, then merge them into the single PRP
        package MF6 wants. Group labels ride along untouched.

        Parameters
        ----------
        *sets : PRTReleasePoints
            The release-point sets to combine, in order. Every record keeps its
            cellid, coordinates, ``local_z`` value and group label; only ``irpt`` is
            rewritten, numbering the combined records ``0, 1, 2, ...`` in the order
            the sets (and their records) were given. Called with no sets it returns
            an empty set, which :class:`PRTProject` refuses.

        Returns
        -------
        PRTReleasePoints
            One set holding every record of ``sets``.

        Examples
        --------
        >>> west = mf.PRTReleasePoints.from_cells(model, [10, 11], group="west_wells")
        >>> east = mf.PRTReleasePoints.from_cells(model, [80, 81], group="east_wells")
        >>> rp = mf.PRTReleasePoints.merge(west, east)
        """

        rows: list[tuple[Any, ...]] = []
        for release_set in sets:
            for row in release_set.packagedata:
                rows.append((len(rows), *tuple(row)[1:]))
        return cls(tuple(rows))

    @classmethod
    def from_cells(
        cls,
        model: SimulationBase,
        cells: Sequence[int | tuple[int, ...]],
        *,
        layer: int = 0,
        local_z: float = 0.5,
        group: str | Sequence[str] | None = None,
    ) -> PRTReleasePoints:
        """Create particle release points at the centers of given model cells.

        Parameters
        ----------
        model : SimulationBase
            The flow model whose grid (``model.gwf.modelgrid``) supplies the cell
            centers. Each particle starts at its cell's FloPy cell center, converted
            to the model's own frame.
        cells : Sequence[int or tuple]
            Cell ids to release from, one particle each, in order -- zero-based DISV
            cell numbers, or structured ``(row, col)`` / ``(layer, row, col)`` tuples.
            A tuple of any other length raises ``ValueError``.
        layer : int, default 0
            Zero-based layer to place the particles in when ``cells`` are bare cell
            numbers or ``(row, col)`` tuples; ignored for a ``(layer, row, col)``
            tuple, which names its own layer.
        local_z : float, default 0.5
            The vertical position written as every record's ``zrpt``. Under
            :class:`PRTProject`'s default ``local_z=True`` it is a local coordinate
            within the cell: 0 = bottom, 1 = effective top (the water table in a
            partially saturated convertible cell). Built with ``local_z=False``, MF6
            reads it as an elevation in model length units instead.
        group : str or Sequence[str], optional
            Release-group label(s) written as PRP boundnames -- one label for every
            point, or one per cell (a sequence of another length raises
            ``ValueError``). MF6 echoes them (uppercased) into the track CSV, which
            is what ``results.capture.map()`` groups particles by. ``None`` (the
            default) writes ungrouped 5-tuple records.

        Returns
        -------
        PRTReleasePoints
            Pass it to ``model.particle_tracking.prt(release_points=...)``.

        Notes
        -----
        The cell ids are TAKEN, not resolved: this neither range-checks them nor
        reads idomain, so an inactive cell is written and FloPy's PRP refuses the
        package later, naming the cellid. :meth:`from_points` resolves a POINT,
        and that one applies the shared coverage rule -- lowest cell active in
        ``layer`` (ledger 242).

        Examples
        --------
        >>> rp = mf.PRTReleasePoints.from_cells(model, cells=[100, 120, 140],
        ...                                     layer=0, local_z=0.5)
        >>> west = mf.PRTReleasePoints.from_cells(model, cells=[100, 120],
        ...                                       group="west_wells")
        """

        cells = list(cells)
        groups = cls._resolve_groups(group, len(cells))
        grid = model.gwf.modelgrid
        rows = []
        for irpt, cell in enumerate(cells):
            if isinstance(cell, tuple):
                if len(cell) == 2:
                    cellid = (int(layer), int(cell[0]), int(cell[1]))
                elif len(cell) == 3:
                    cellid = tuple(int(value) for value in cell)
                else:
                    raise ValueError("Structured PRT cell tuples must be (row, col) or (layer, row, col).")
                row, col = cellid[-2:]
                x = grid.xcellcenters[row, col]
                y = grid.ycellcenters[row, col]
            else:
                cell = int(cell)
                cellid = (int(layer), cell)
                x = grid.xcellcenters[cell]
                y = grid.ycellcenters[cell]
            # FloPy's cell centres are WORLD coordinates; a PRP record is in the
            # model's own frame. On an offset grid -- every USG import -- those are
            # different places, and MF6 released the particles off the grid.
            x, y = grid.get_local_coords(x, y)
            record: tuple[Any, ...] = (
                irpt,
                cellid,
                float(x),
                float(y),
                float(local_z),
            )
            rows.append(record if groups is None else (*record, groups[irpt]))
        return cls(tuple(rows))

    @classmethod
    def from_points(
        cls,
        model: SimulationBase,
        points,
        *,
        layer: int = 0,
        local_z: float = 0.5,
        group: str | Sequence[str] | None = None,
    ) -> PRTReleasePoints:
        """Create release points from point geometries or an iterable of ``(x, y)``.

        ``group`` labels the points as one release group (or one label per point) --
        see :meth:`from_cells`. Points are WORLD coordinates, in the grid's CRS (a
        frame in another one is reprojected to it); they are written in the model's
        own frame, which on an offset grid differs.

        Parameters
        ----------
        model : SimulationBase
            The flow model whose grid the points are located on. On a DISV model
            carrying a ``vor``, its cell polygons and ``idomain`` decide the cell
            (see Notes); on a structured model, or a vertex model with no ``vor``,
            FloPy's ``model.gwf.modelgrid.intersect`` does.
        points : iterable of shapely.Point or (x, y) pairs, GeoSeries, or GeoDataFrame
            The release locations, one particle each, in WORLD coordinates. Anything
            with ``.x``/``.y`` is read as a point, anything else as an ``(x, y)``
            pair; a GeoDataFrame is read through its ``geometry`` column. A
            GeoSeries or GeoDataFrame whose CRS differs from the grid polygons' is
            reprojected to it first (DISV models with a ``vor`` only, and only when
            both carry a CRS); bare pairs carry no CRS and are taken as the grid's.
        layer : int, default 0
            Zero-based layer every point is released in. On a DISV model with a
            ``vor``, only cells active (``idomain > 0``) in this layer can take a
            point.
        local_z : float, default 0.5
            The vertical position written as every record's ``zrpt`` -- a local z
            within the cell (0 = bottom, 1 = effective top) under :class:`PRTProject`'s
            default ``local_z=True``, an elevation otherwise. See :meth:`from_cells`.
        group : str or Sequence[str], optional
            Release-group label(s) written as PRP boundnames -- one label for every
            point, or one per point in input order (a sequence of another length
            raises ``ValueError``). A point left out as a miss drops its label with
            it. ``None`` (the default) writes ungrouped records.

        Returns
        -------
        PRTReleasePoints
            One record per point that landed, ``irpt`` numbered over those only.

        Raises
        ------
        ValueError
            On a DISV model with a ``vor``, when no point lands on a cell active in
            ``layer``; the message names both extents and both CRSs. Also when
            ``group`` has the wrong length, and (same models) when ``layer`` is past
            the last layer of a layered ``(nlay, ncpl)`` idomain.

        Notes
        -----
        A point resolves to the ONE cell that holds it: on a shared face or
        vertex, the lowest-numbered cell ACTIVE in ``layer`` (ledger 242). FloPy's
        ``modelgrid.intersect`` answers with the lowest CELL2D and knows nothing
        about idomain, so a point on the face between an inactive cell and an
        active one resolved to the inactive cell, and one on a vertex to whichever
        of its four cells is numbered lowest. That cost the run, late and
        obscurely: FloPy's PRP checks every cellid, so ``PRTProject`` raised
        ``Cellid (0, 4) is outside of the active model grid`` -- an id the caller
        never wrote -- where the active neighbour the point also touches releases
        and tracks.

        A point landing on no active cell is left out and named in a warning, and
        ``irpt`` is renumbered over the survivors as :meth:`merge` does, because
        MF6 wants release-point numbers 1..NRELEASEPTS with none missing. If NONE
        of them land that raises, naming both extents and both CRSs (ledger 211),
        where FloPy's own refusal named only the point. One point off the grid used
        to end the whole call.

        A structured (DIS) grid keeps FloPy's ``intersect``, and with it FloPy's
        own refusal: the shared rule is indexed by cell id and a DIS cell is
        ``(row, col)``, so neither the activity rule nor the miss handling above
        reaches a structured model.
        """

        given = _release_frame(points)
        groups = cls._resolve_groups(group, len(given))
        grid = model.gwf.modelgrid
        polys = _cell_frame(model)
        located = given if polys is None else _in_grid_crs(given, polys)
        domain = None if polys is None else model_domain(model, len(polys))
        rows: list[tuple[Any, ...]] = []
        missed: list[tuple[int, Any]] = []
        for index, point in enumerate(located):
            if polys is None:
                cellid: tuple[int, ...] = _flopy_cellid(grid, point, layer)
            else:
                hits = cell_shares(
                    polys,
                    point,
                    # Eligibility only. ``cell_shares`` sorts whatever this
                    # returns, so the lowest-numbered active cell wins a tie
                    # without a second copy of that rule here.
                    usable=lambda cells: [
                        cell for cell in cells if domain_active(domain, layer, cell)
                    ],
                )
                if not hits:
                    missed.append((index, point))
                    continue
                cellid = (int(layer), int(hits[0][0]))
            local_x, local_y = grid.get_local_coords(float(point.x), float(point.y))
            record: tuple[Any, ...] = (
                len(rows), cellid, float(local_x), float(local_y), float(local_z)
            )
            rows.append(record if groups is None else (*record, groups[index]))
        if missed:
            _refuse_or_report_missed(polys, given, missed, layer=layer)
        return cls(tuple(rows))


@dataclass
class PRTRunResults:
    """Read access to the pathlines/output of a finished PRT simulation.

    A file-backed handle to one completed PRT run: it points at the track CSV (and
    optional budget) in ``workspace`` and exposes the simulated particle pathlines
    for plotting and analysis, joined back to the originating ``flow_model``'s grid.
    Created for you by :class:`PRTProject` after a run, or reopened from disk with
    :func:`open_prt_run` without rebuilding the simulation.

    Results are read through view nouns, each answering the usual
    ``get``/``summary``/``plot``/``map``/``mosaic`` verbs
    (``myflopy.modflow.mf6.prt_maps``)::

        results.pathlines.map()                  # the tracks over the head map
        results.pathlines.get()                  # normalized track records
        results.travel_time.map(stat="median")   # time-of-travel choropleth
        results.endpoints.get()                  # counts per terminating cell
        results.capture.map()                    # one panel per release group

    The three per-cell maps are time-integrated over the whole run, so they take
    no ``per``. :attr:`pathlines` is a **view**, not a frame -- call ``.get()``
    for the records (or :attr:`track_records` for the untouched CSV).

    Attributes
    ----------
    flow_model : SimulationBase
        The flow model the particles were tracked through; its grid is what the
        pathlines are drawn and aggregated on.
    workspace : pathlib.Path
        Directory holding the PRT output files.
    name : str
        PRT model name (the stem of its output files).
    track_csv_path : pathlib.Path
        The particle-track CSV MF6 wrote (``<name>.trk.csv`` in ``workspace``).
        Nothing checks it exists until :attr:`track_records` reads it, which raises
        ``FileNotFoundError`` when it does not.
    budget_path : pathlib.Path or None, optional
        The PRT budget file. :class:`PRTProject` and :func:`open_prt_run` both set
        ``<name>.bud`` in ``workspace`` without checking it exists; ``None`` by
        default when the dataclass is built directly.
    success : bool or None, optional
        MF6's run verdict, set by :meth:`PRTProject.run` (always ``True`` there,
        because a failed run raises instead of returning). ``None`` (the default)
        when the handle was reopened with :func:`open_prt_run` or built by
        :meth:`PRTProject.results` without one.
    report : tuple of str, optional
        The MF6 screen output captured by :meth:`PRTProject.run` with
        ``report=True``; empty (the default) otherwise.
    """

    flow_model: SimulationBase
    workspace: Path
    name: str
    track_csv_path: Path
    budget_path: Path | None = None
    success: bool | None = None
    report: tuple[str, ...] = ()

    @property
    def engine(self) -> str:
        """The particle-tracking engine that produced these results (``"mf6-prt"``)."""

        return "mf6-prt"

    @cached_property
    def track_records(self) -> pd.DataFrame:
        """Load and cache the PRT track CSV exactly as MF6 wrote it.

        The raw table: one-based ``icell``, no derived columns. FloPy and PyVista
        consume it directly; for analysis prefer ``pathlines.get()``, which adds
        zero-based ``cell``/``layer``, ``travel_time``, and the release group.
        """

        if not self.track_csv_path.exists():
            raise FileNotFoundError(f"PRT track CSV not found: {self.track_csv_path}")
        return pd.read_csv(self.track_csv_path)

    def refresh(self) -> PRTRunResults:
        """Clear cached file-backed results so subsequent access rereads disk."""

        self.__dict__.pop("track_records", None)
        return self

    @property
    def terminal_points(self) -> pd.DataFrame:
        """The terminating pathline records (``ireason == 3``), one per particle endpoint."""

        data = self.track_records
        if "ireason" not in data.columns:
            return data.iloc[0:0].copy()
        return data.loc[data["ireason"] == 3].copy()

    @staticmethod
    def _maps():
        """The derived-map views module (imported lazily -- it pulls in the plot layer)."""

        from myflopy.modflow.mf6 import prt_maps

        return prt_maps

    @property
    def pathlines(self):
        """Pathline view: the tracks themselves (``get``/``summary``/``plot``/``map``/``mosaic``).

        A **view**, not a DataFrame -- ``pathlines.get()`` returns the normalized
        records and ``pathlines.map()`` draws them over the model map.
        """

        return self._maps().PRTPathlineView(self)

    @property
    def travel_time(self):
        """Travel-time view: ``get``/``summary``/``plot``/``map``/``mosaic`` (time of travel per cell)."""

        return self._maps().PRTTravelTimeView(self)

    @property
    def endpoints(self):
        """Endpoint view: particle-termination counts per cell, with the same verbs."""

        return self._maps().PRTEndpointsView(self)

    @property
    def capture(self):
        """Capture-zone view: endpoints split by release group (``map()`` gives one panel per group)."""

        return self._maps().PRTCaptureView(self)

    def scene(self, **kwargs):
        """The 3-D pathline scene, as a :class:`~myflopy.viz.VtkScene` Picture.

        The same picture as ``flow_model.plot.grid(pathlines=..., backend="vtk")``
        -- this is the spelling you reach for when you already have the run
        results in hand. Display it, ``.show()`` it, or ``.html(path)`` it.
        Requires ``pyvista``.

        Parameters
        ----------
        **kwargs
            Forwarded unchanged to
            :func:`~myflopy.modflow.mf6.interactive_plotting.build_particle_tracking_scene`,
            which accepts ``vertical_exaggeration`` (float, default 1.0),
            ``model_style`` (PyVista mesh style for the grid, default
            ``"wireframe"``), ``model_opacity`` (default 0.25), ``show_edges``
            (default ``True``), ``pathline_cmap`` (colormap for the time-colored
            tubes, default ``"viridis"``), ``pathline_width`` (default 4.0) and
            ``off_screen`` (default ``True``). Any other keyword raises
            ``TypeError``. The flow model and the raw :attr:`track_records` are
            passed for you.

        Returns
        -------
        VtkScene
            The grid and pathline tubes, as a :class:`~myflopy.viz.Picture`.
        """

        from myflopy.modflow.mf6.interactive_plotting import build_particle_tracking_scene

        return build_particle_tracking_scene(self.flow_model, self.track_records, **kwargs)

    def plot_map(self, **kwargs):
        """Plot the particle pathlines on a plan-view map of the flow model.

        A static Matplotlib map drawn by FloPy's ``PlotMapView``. For the
        interactive picture over the head map, use ``pathlines.map()``.

        Parameters
        ----------
        **kwargs
            Forwarded unchanged to
            :func:`~myflopy.modflow.mf6.interactive_plotting.plot_particle_pathlines`,
            which accepts ``layer`` (``"all"``, the default, draws every track
            point; an int keeps the points whose layer column equals it, compared as
            given with no base conversion), ``ax`` (an existing Matplotlib Axes; a
            new figure when ``None``), ``head_data`` / ``head_layer`` (an optional
            head array
            shaded beneath the tracks, and the layer of it to draw, default 0),
            ``pathline_color``, ``pathline_alpha`` (default 0.7),
            ``pathline_linewidth`` (default 1.0), ``show_grid`` (default ``True``),
            ``grid_alpha`` (default 0.25) and ``title``. Any other keyword raises
            ``TypeError``. The flow model and the raw :attr:`track_records` are
            passed for you.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """

        from myflopy.modflow.mf6.interactive_plotting import plot_particle_pathlines

        return plot_particle_pathlines(self.flow_model, self.track_records, **kwargs)

    def export_3d_html(self, output_path: str | Path, **kwargs) -> Path:
        """Write the 3-D pathline scene to a standalone interactive HTML file.

        Kept as a named convenience because "give me a file I can send" is a
        distinct intent; it is exactly ``self.scene(**kwargs).html(output_path)``,
        and closes the plotter afterwards.

        Parameters
        ----------
        output_path : str or pathlib.Path
            The HTML file to write. Missing parent directories are created, and an
            existing file is overwritten.
        **kwargs
            Forwarded to :meth:`scene`, and from there to
            ``build_particle_tracking_scene`` -- the same keywords :meth:`scene`
            lists; any other raises ``TypeError``.

        Returns
        -------
        pathlib.Path
            The path written.
        """

        scene = self.scene(**kwargs)
        try:
            return scene.html(output_path)
        finally:
            scene.scene.close()


class PRTProject:
    """Build, run, and collect an MF6 particle-tracking run off an existing flow model.

    The high-level driver for native MF6 PRT: given a completed/loaded flow
    ``model`` and a set of ``release_points``, it assembles the PRT model (MIP/PRP/OC
    with the supplied porosity, tracking-time, and termination options), couples it
    to the flow model, writes and runs it in ``workspace``, and returns a
    :class:`PRTRunResults` for the pathlines. Construct it through
    ``model.particle_tracking.prt(...)`` (the :class:`ParticleTracking` front door)
    or directly. For the declarative spec API, see ``mf.prt(...)``.

    Building one only assembles the FloPy objects (and creates ``workspace``); call
    :meth:`write` and :meth:`run` to produce output. The PRT simulation reads the
    flow model's head, budget and binary-grid files through FMI, by paths relative
    to ``workspace``, so the flow model must have been run first.

    Parameters
    ----------
    model : SimulationBase
        The flow model to track particles through. Its DIS or DISV discretization
        is copied onto the PRT model (any other grid raises
        ``NotImplementedError``), its TDIS supplies the defaults for
        ``perioddata``/``time_units``, and its OC records name the head and budget
        files FMI reads (``<name>.hds``/``<name>.cbc`` when OC names none).
    workspace : str or pathlib.Path
        Directory the PRT simulation is written to and run in; created (with its
        parents) at construction.
    name : str, optional
        PRT model and simulation name, also the stem of the output files
        (``<name>.trk.csv``, ``<name>.bud``). ``None`` (the default) derives it from
        the flow model's name: non-alphanumerics become ``_``, the first 12
        characters are kept and ``_prt`` is appended. More than 16 characters
        raises ``ValueError`` (MF6's limit).
    release_points : PRTReleasePoints or Sequence of tuple
        Where particles start: a :class:`PRTReleasePoints`, or raw PRP
        ``packagedata`` records in the same layout (each converted to a tuple).
        Required, and an empty set raises ``ValueError``. Records carrying a sixth
        (group) field switch the PRP ``BOUNDNAMES`` option on.
    porosity : float or Sequence[float] or array_like, default 0.2
        Aquifer porosity used to convert flow to seepage velocity, passed unchanged
        to FloPy's ``ModflowPrtmip`` ``porosity``: one value for every cell, a
        per-layer list (one value per layer), or a per-cell array.
    perioddata : Sequence of (float, int, float), optional
        The PRT simulation's own TDIS ``(perlen, nstp, tsmult)`` rows. ``None`` (or
        an empty sequence) copies the flow model's TDIS periods, or one period
        ``(1.0, 1, 1.0)`` when it has no TDIS period data. PRT tracks through the flow
        model's saved flows, so these should match the flow run's periods.
    time_units : str, optional
        TDIS time units for the PRT simulation: one of MF6's words ``"unknown"``,
        ``"seconds"``, ``"minutes"``, ``"hours"``, ``"days"``, ``"years"``, any case
        (e.g. ``"DAYS"``); any other word raises ``ValueError`` before anything is
        written, naming MF6's word for a known short spelling -- MF6 stops on it
        (ledger 456). ``None`` (the
        default) copies the flow simulation's TDIS ``TIME_UNITS``, read through
        :func:`~myflopy.modflow.mf6.hydraulics.declared_time_unit`; when the flow
        simulation declares none (a plain ``mf.tdis(...)``, ``unknown``, or no TDIS
        at all) the PRT TDIS declares none either, which MF6 reads as unknown
        units -- FloPy answers an unset option with ``None``, and writing that
        through ``str()`` gave ``TIME_UNITS none``, which MF6 refuses (ledger 397). Travel
        times are in the flow model's time unit whatever this says -- it only
        labels the PRT listing.
    release_perioddata : dict, optional
        PRP ``perioddata``: zero-based stress period -> that period's MF6 release
        settings (``ALL``, ``FIRST``, ``LAST``, ``FREQUENCY``, ``STEPS``), in the
        form FloPy's ``ModflowPrtprp`` takes and passed to it unchanged. ``None``
        (the default) passes ``{0: ["FIRST"]}``: release at the start of the first
        time step of stress period 0.
    extend_tracking : bool, default False
        PRP ``EXTEND_TRACKING``: keep tracking past the end of the simulation's
        final time step, on that step's flows, until particles terminate or reach
        ``stoptime``/``stoptraveltime``. By default MF6 terminates particles at the
        end of the final time step.
    stoptime : float, optional
        PRP ``STOPTIME``: the simulation time (model time units) at which every
        particle still moving terminates. ``None`` (the default) sets none.
    stoptraveltime : float, optional
        PRP ``STOPTRAVELTIME``: the travel time (model time units) after which a
        particle terminates. ``None`` (the default) sets none. When both stop
        options are given, whichever is reached first stops the particle.
    drape : bool, default True
        PRP ``DRAPE``: a particle whose release cell is inactive at release time is
        moved to the topmost active cell below it, if any. With ``False`` MF6 does
        not release such a particle. It does not apply to a cell with
        ``idomain <= 0``, which MF6 treats as an error.
    local_z : bool, default True
        PRP ``LOCAL_Z``: read each record's z (``zrpt``) as a local coordinate in
        the cell, 0 at the bottom and 1 at the effective top (the water table in a
        partially saturated convertible cell). ``False`` reads it as an elevation
        in model length units. :meth:`PRTReleasePoints.from_cells` and
        :meth:`~PRTReleasePoints.from_points` write a local value (0.5 by default),
        so they expect the default.
    stop_at_weak_sink : bool, default False
        PRP ``STOP_AT_WEAK_SINK``: terminate a particle when it enters a weak-sink
        cell. By default particles pass through weak sinks.
    track_times : Sequence[float], optional
        Times from the simulation start (model time units) written to the PRT OC
        ``TRACKTIMES`` block, at which MF6 adds a track record (``ireason`` 5) for
        each particle still being tracked. ``None`` (the default) writes no block.
    exe_name : str, default "mf6"
        The MODFLOW 6 executable FloPy runs (a name on ``PATH`` or a path).
    """

    def __init__(
        self,
        model: SimulationBase,
        *,
        workspace: str | Path,
        name: str | None = None,
        release_points: PRTReleasePoints | Sequence[tuple[Any, ...]],
        porosity: float | Sequence[float] = 0.2,
        perioddata: Sequence[tuple[float, int, float]] | None = None,
        time_units: str | None = None,
        release_perioddata: dict | None = None,
        extend_tracking: bool = False,
        stoptime: float | None = None,
        stoptraveltime: float | None = None,
        drape: bool = True,
        local_z: bool = True,
        stop_at_weak_sink: bool = False,
        track_times: Sequence[float] | None = None,
        exe_name: str = "mf6",
    ):
        """Configure a PRT run off ``model``: validate the name/release points and set output paths.

        See the class docstring for the tracking-time, placement, and termination
        parameters. Raises if the name exceeds 16 characters or no release points
        are given.
        """

        self.flow_model = model
        self.workspace = Path(workspace)
        self.name = _default_prt_name(model.name) if name is None else str(name)
        if len(self.name) > 16:
            raise ValueError("MF6 PRT model names cannot exceed 16 characters.")
        self.release_points = (
            release_points
            if isinstance(release_points, PRTReleasePoints)
            else PRTReleasePoints(tuple(tuple(row) for row in release_points))
        )
        if not self.release_points.packagedata:
            raise ValueError("At least one PRT release point is required.")
        # The declaration rule of mf.tdis (ledger 456): MF6 stops on any other word.
        check_tdis_time_units(time_units, helper="PRTProject")
        self.workspace.mkdir(parents=True, exist_ok=True)
        self.track_csv_path = self.workspace / f"{self.name}.trk.csv"
        self.budget_path = self.workspace / f"{self.name}.bud"

        self.sim = flopy.mf6.MFSimulation(
            sim_name=self.name,
            version="mf6",
            exe_name=exe_name,
            sim_ws=self.workspace,
        )
        self.tdis = flopy.mf6.ModflowTdis(
            self.sim,
            time_units=time_units or declared_time_unit(model),
            nper=len(perioddata or _tdis_perioddata(model)),
            perioddata=list(perioddata or _tdis_perioddata(model)),
        )
        self.prt = flopy.mf6.ModflowPrt(
            self.sim,
            modelname=self.name,
            model_nam_file=f"{self.name}.nam",
        )
        self.discretization = _copy_grid_to_prt(model, self.prt)
        self.mip = flopy.mf6.ModflowPrtmip(self.prt, pname="mip", porosity=porosity)
        self.prp = flopy.mf6.ModflowPrtprp(
            self.prt,
            pname="prp",
            # Release-group labels are PRP boundnames; without BOUNDNAMES the sixth
            # field is rejected, so the flag follows the record width.
            boundnames=self.release_points.has_groups or None,
            nreleasepts=len(self.release_points.packagedata),
            packagedata=list(self.release_points.packagedata),
            perioddata={0: ["FIRST"]} if release_perioddata is None else release_perioddata,
            extend_tracking=extend_tracking,
            stoptime=stoptime,
            stoptraveltime=stoptraveltime,
            drape=drape,
            local_z=local_z,
            stop_at_weak_sink=stop_at_weak_sink,
            coordinate_check_method=None,
        )
        oc_kwargs = {}
        if track_times is not None:
            oc_kwargs.update(
                ntracktimes=len(track_times),
                tracktimes=[(float(value),) for value in track_times],
            )
        self.oc = flopy.mf6.ModflowPrtoc(
            self.prt,
            pname="oc",
            budget_filerecord=[self.budget_path.name],
            trackcsv_filerecord=[self.track_csv_path.name],
            saverecord=[("BUDGET", "ALL")],
            **oc_kwargs,
        )
        flow_workspace = Path(model.model_output_folder_path)
        self.fmi = flopy.mf6.ModflowPrtfmi(
            self.prt,
            pname="fmi",
            packagedata=[
                ("GWFHEAD", Path(os.path.relpath(flow_workspace / _output_filename(model, "head"), self.workspace)).as_posix()),
                (
                    "GWFBUDGET",
                    Path(os.path.relpath(flow_workspace / _output_filename(model, "budget"), self.workspace)).as_posix(),
                ),
                (
                    "GWFGRID",
                    Path(os.path.relpath(flow_workspace / _grid_filename(model), self.workspace)).as_posix(),
                ),
            ],
        )
        self.ems = flopy.mf6.ModflowEms(self.sim, pname="ems", filename=f"{self.name}.ems")

    def write(self, *, silent: bool = True) -> Path:
        """Write the PRT simulation files to the workspace and return its path.

        Parameters
        ----------
        silent : bool, default True
            Passed to FloPy's ``MFSimulation.write_simulation``; ``False`` lets
            FloPy print its progress as it writes each file.

        Returns
        -------
        pathlib.Path
            :attr:`workspace`.
        """

        self.sim.write_simulation(silent=silent)
        return self.workspace

    def results(self, *, success: bool | None = None, report: Sequence[str] = ()) -> PRTRunResults:
        """Build a :class:`PRTRunResults` handle over this project's output files.

        Nothing is read or checked here: the handle points at ``<name>.trk.csv`` and
        ``<name>.bud`` in the workspace, and the track CSV is read on first access.

        Parameters
        ----------
        success : bool or None, optional
            Stored as :attr:`PRTRunResults.success`. :meth:`run` passes MF6's
            verdict; ``None`` (the default) means no run was observed.
        report : Sequence[str], default ()
            MF6 screen-output lines, stored as :attr:`PRTRunResults.report` (as a
            tuple). :meth:`run` passes the lines it captured.

        Returns
        -------
        PRTRunResults
        """

        return PRTRunResults(
            flow_model=self.flow_model,
            workspace=self.workspace,
            name=self.name,
            track_csv_path=self.track_csv_path,
            budget_path=self.budget_path,
            success=success,
            report=tuple(report),
        )

    def run(self, *, write: bool = True, silent: bool = True, report: bool = True) -> PRTRunResults:
        """Write (optionally) and run the PRT simulation, returning its :class:`PRTRunResults`.

        Requires the GWF head/budget/grid outputs to already exist (raises with the
        missing paths otherwise); raises ``RuntimeError`` if the MF6 run fails.

        Parameters
        ----------
        write : bool, default True
            Write the simulation files first (:meth:`write`, with the same
            ``silent``). ``False`` runs whatever is already in the workspace.
        silent : bool, default True
            Passed to FloPy's ``write_simulation`` and ``run_simulation``; ``False``
            echoes FloPy's and MF6's screen output as they run.
        report : bool, default True
            Passed to FloPy's ``run_simulation``: capture MF6's screen output. The
            captured lines become :attr:`PRTRunResults.report`, and on a failed run
            their last 20 lines are quoted in the ``RuntimeError``. ``False``
            captures nothing, so both are empty.

        Returns
        -------
        PRTRunResults
            A handle over this run's output, with ``success=True``.

        Raises
        ------
        FileNotFoundError
            When the flow model's head, budget or binary-grid file is missing from
            its output folder (checked before anything is written).
        RuntimeError
            When MF6 does not report normal termination.
        """

        missing = [
            path
            for path in (
                Path(self.flow_model.model_output_folder_path) / _output_filename(self.flow_model, "head"),
                Path(self.flow_model.model_output_folder_path) / _output_filename(self.flow_model, "budget"),
                Path(self.flow_model.model_output_folder_path) / _grid_filename(self.flow_model),
            )
            if not path.exists()
        ]
        if missing:
            formatted = "\n".join(f"- {path}" for path in missing)
            raise FileNotFoundError(
                "PRT requires completed GWF head and budget outputs before it can run:\n"
                f"{formatted}"
            )
        if write:
            self.write(silent=silent)
        success, output = self.sim.run_simulation(silent=silent, report=report)
        if not success:
            tail = "\n".join(output[-20:]) if output else "No MF6 report was returned."
            raise RuntimeError(f"PRT simulation failed in {self.workspace}\n{tail}")
        return self.results(success=success, report=output)


def open_prt_run(
    model: SimulationBase,
    workspace: str | Path,
    *,
    name: str | None = None,
) -> PRTRunResults:
    """Reopen a finished PRT run from disk as :class:`PRTRunResults`.

    The review counterpart to :class:`PRTProject`: point it at a ``workspace`` that
    already contains PRT output and it locates the track/budget files and returns a
    results handle -- without re-assembling or re-running the PRT simulation. Use it
    to revisit pathlines from an earlier run.

    Parameters
    ----------
    model : SimulationBase
        The flow model the particles were tracked through; stored as
        :attr:`PRTRunResults.flow_model` and used for every map and scene.
    workspace : str or pathlib.Path
        Directory containing the completed PRT output.
    name : str, optional
        PRT model name: the handle points at ``<name>.trk.csv`` and ``<name>.bud``
        in ``workspace``, neither checked here. ``None`` (the default) infers it
        from the one ``*.trk.csv`` file in ``workspace``.

    Returns
    -------
    PRTRunResults
        A handle to the existing run's pathline output (``success`` is ``None``
        and ``report`` empty -- the run was not observed).

    Raises
    ------
    FileNotFoundError
        With ``name=None``, when ``workspace`` holds no ``*.trk.csv`` file or more
        than one.
    """

    workspace = Path(workspace)
    if name is None:
        csv_candidates = sorted(workspace.glob("*.trk.csv"))
        if len(csv_candidates) != 1:
            raise FileNotFoundError(
                f"Expected one *.trk.csv file in {workspace}, found {len(csv_candidates)}."
            )
        track_path = csv_candidates[0]
        name = track_path.name.removesuffix(".trk.csv")
    else:
        track_path = workspace / f"{name}.trk.csv"
    return PRTRunResults(
        flow_model=model,
        workspace=workspace,
        name=name,
        track_csv_path=track_path,
        budget_path=workspace / f"{name}.bud",
    )


class ParticleTracking:
    """A model's particle-tracking front door (native MF6 PRT and legacy MP3DU).

    The small accessor exposed as ``model.particle_tracking``: it bundles the
    particle-tracking workflows available for a flow ``model`` so you do not have to
    import the project classes yourself. :meth:`prt` builds a native MF6 PRT
    simulation (returning a :class:`PRTProject`, which you then ``.run()``),
    :meth:`open_prt` reopens a finished one, and :meth:`mp3du` covers the older
    mod-PATH3DU workflow.

    Parameters
    ----------
    model : SimulationBase
        The flow model these particle-tracking workflows operate on; every method
        passes it as the flow model.
    """

    def __init__(self, model: SimulationBase):
        """Bind the particle-tracking front door to a flow ``model``."""

        self.model = model

    def prt(self, *, workspace: str | Path, release_points, **kwargs) -> PRTProject:
        """Build a native MF6 PRT project for this model (see :class:`PRTProject`).

        Nothing is written or run: call ``.run()`` on the returned project.

        Parameters
        ----------
        workspace : str or pathlib.Path
            Directory the PRT simulation is written to and run in (created now).
            See :class:`PRTProject`.
        release_points : PRTReleasePoints or Sequence of tuple
            Where particles start -- a :class:`PRTReleasePoints` (from
            :meth:`PRTReleasePoints.from_cells` / ``from_points`` / ``merge``) or raw
            PRP records. See :class:`PRTProject`.
        **kwargs
            Forwarded to :class:`PRTProject`: ``name``, ``porosity``,
            ``perioddata``, ``time_units``, ``release_perioddata``,
            ``extend_tracking``, ``stoptime``, ``stoptraveltime``, ``drape``,
            ``local_z``, ``stop_at_weak_sink``, ``track_times`` and ``exe_name``,
            with the defaults and meanings documented there. Any other keyword
            raises ``TypeError``.

        Returns
        -------
        PRTProject
        """

        return PRTProject(
            self.model,
            workspace=workspace,
            release_points=release_points,
            **kwargs,
        )

    def open_prt(self, workspace: str | Path, *, name: str | None = None) -> PRTRunResults:
        """Reopen a completed PRT run in ``workspace`` (see :func:`open_prt_run`).

        Parameters
        ----------
        workspace : str or pathlib.Path
            Directory containing the completed PRT output.
        name : str, optional
            PRT model name, the stem of ``<name>.trk.csv``. ``None`` (the default)
            infers it from the one ``*.trk.csv`` file in ``workspace`` and raises
            ``FileNotFoundError`` when there is none or more than one.

        Returns
        -------
        PRTRunResults
        """

        return open_prt_run(self.model, workspace, name=name)

    def mp3du(self, particles, **kwargs):
        """Prepare a legacy MP3DU (MODPATH-style) particle-tracking run for this model.

        Exactly :func:`~myflopy.modflow.mp3du.prepare_particle_tracking` with
        ``model`` bound, so by default it writes the mod-PATH3DU inputs and does NOT
        execute them (``execute=False``).

        Parameters
        ----------
        particles : str, pathlib.Path, or Sequence[int]
            The release source: a shapefile/GeoPackage path, or zero-based cell ids
            (one particle at each cell's centroid). See
            :func:`~myflopy.modflow.mp3du.prepare_particle_tracking`.
        **kwargs
            Forwarded to ``prepare_particle_tracking`` -- its ``porosity``,
            ``output_path``, ``execute``, ``convert_output`` and
            ``write_diagnostics``, and through it any
            :class:`~myflopy.modflow.mp3du.ParticleTrackingInput` keyword
            (``direction``, ``simulation_end_time``, ``iface_overrides``, ...). A
            keyword neither accepts raises ``TypeError``.

        Returns
        -------
        tuple of (ParticleTrackingInput, ParticleTrackingResult)
            The configured tracker and the result of its ``run()``.
        """

        from myflopy.modflow.mp3du import prepare_particle_tracking

        return prepare_particle_tracking(model=self.model, particles=particles, **kwargs)
