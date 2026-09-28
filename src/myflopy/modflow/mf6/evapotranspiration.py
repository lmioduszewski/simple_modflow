"""File-less, domain-aware EVT (evapotranspiration) package builder.

:class:`EVTBuilder` is the ET analogue of
:class:`~myflopy.modflow.mf6.recharge.RCHBuilder`: it self-selects top-active
cells from the model domain and broadcasts the ET ``rate`` and extinction
``depth`` across them, exactly like RCH broadcasts ``recharge``. The one piece
with no RCH counterpart is the per-cell ET ``surface`` elevation, which is
resolved through the shared :class:`~myflopy.geopackage.SurfaceResolver` /
:class:`~myflopy.geopackage.CellSurfaceOffset` engine rather than by reading the
grid's ``gdf_topbtm`` columns directly.

This is the engine under ``mf.evt(context=..., nper=..., rate=..., depth=...)``;
use that facade, not this class, to build models.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Real
from typing import Any

from myflopy.advanced import evt_spec
from myflopy.geopackage import CellSurfaceOffset, SurfaceResolver
from myflopy.modflow.mf6.areal import CellId, _ArealBuilder
from myflopy.specs import PackageSpec


@dataclass(frozen=True, slots=True, kw_only=True)
class EVTBuilder(_ArealBuilder):
    """Prepare a list-based EVT package from explicit ET configuration.

    ``rate`` (maximum ET flux) and ``depth`` (extinction depth) accept the same
    shapes as :class:`~myflopy.modflow.mf6.recharge.RCHBuilder`'s ``recharge``:
    a scalar / time-series name, a per-cell sequence, a ``{cellid: value}``
    mapping, or a ``{period: ...}`` mapping.

    ``surface`` is the ET surface *elevation* per cell and defaults to the model
    top / land surface (``CellSurfaceOffset("model_top")``, resolved through the
    surface engine). ``model_top`` resolves for every column regardless of which
    layer is top-active, which is the physically sensible ET reference; use
    ``CellSurfaceOffset("cell_top")`` only when the ET surface should follow the
    top of a deeper top-active cell (which then needs per-layer top surfaces).
    ``surface`` also accepts a numeric constant, an explicit per-cell sequence /
    ``{cellid: value}`` mapping, or a ``CellSurfaceOffset`` with a numeric
    offset/minimum (e.g. ``CellSurfaceOffset("model_top", offset=-2.0)`` for two
    length-units below ground). Field-name (string) offsets require the
    ``mf.evt.gpkg(...)`` path -- there are no feature rows here to name a field.

    Only the single-segment record ``(cellid, surface, rate, depth)`` is
    produced (``nseg=1``); segmented ET (``pxdp``/``petm``) is a direct-form
    feature (see the compromise ledger).

    All parameters are keyword-only. Only ``nper``, a string ``cells`` and the
    ``surface`` form are checked when the builder is created; everything else waits
    for :meth:`build`, which validates, writes one record per selected cell for
    every period, and returns the ``PackageSpec``, stamped with the pin of
    ``context``'s grid when it has one.

    Parameters
    ----------
    context : ModelContext
        Where the cells come from: ``context.domain`` (the idomain, ``(ncpl,)`` or
        ``(nlay, ncpl)``) selects and checks them, and with no domain every cell of
        ``context.grid.ncpl`` is used. A ``CellSurfaceOffset`` surface also reads
        its elevations here -- ``context.surfaces``, else ``context.grid.gdf_topbtm``
        -- and raises ``ValueError`` at :meth:`build` when there are none. A context
        with neither a domain nor a grid raises ``ValueError`` at :meth:`build` for
        a generated selection.
    nper : int
        Number of stress periods; every period ``0 .. nper - 1`` gets a record for
        every selected cell. A value below 1 raises ``ValueError``.
    cells : {"top_active", "all_active", "surface_only"} or sequence, default "top_active"
        Which cells receive ET. ``"top_active"`` takes, per column, the highest
        layer with ``idomain > 0`` and skips a column with none; ``"all_active"``
        takes every active ``(layer, cell)``; ``"surface_only"`` is accepted and
        currently behaves exactly like ``"top_active"``. With no domain all three
        take every column of the grid, and with a one-dimensional domain its active
        columns, all on ``layer``. A sequence names the cells itself:
        ``(layer, cell)`` pairs, or bare zero-based cell numbers placed on
        ``layer``; at :meth:`build` each is checked against the domain (when there
        is one), and an inactive cell or one past the domain's extent raises
        ``ValueError``. Any other string raises ``ValueError``.
    layer : int, default 0
        Zero-based layer for bare cell numbers in ``cells``, and for every
        generated cell when the domain is missing or one-dimensional. Ignored by
        the string selectors on an ``(nlay, ncpl)`` domain.
    name_by_cell : Mapping, optional
        Boundary names, keyed by ``(layer, cell)`` or by bare cell number (the
        pair is tried first). Read only with ``boundnames=True``; a cell it does not
        name falls back to ``"<name>_<layer>_<cell>"``.
    boundnames : bool, default False
        Append a boundary name to every record and write MF6's ``BOUNDNAMES``
        option.
    name : str, default "evt"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>``, the
        built spec's name, and the prefix of generated boundary names.
    options : Mapping[str, Any], optional
        Extra keywords for :class:`flopy.mf6.ModflowGwfevt` (``fixed_cell``,
        ``print_input``, ``timeseries``, ...), passed through
        :func:`~myflopy.advanced.evt_spec` and checked against that class's
        signature at :meth:`build`, where an unknown name raises ``TypeError``. They
        may override ``save_flows=True``. ``nseg`` is always passed as 1, so an
        ``nseg`` here collides with it and raises ``TypeError``, as does a ``name``
        or ``boundnames``. Defaults to none.
    rate : float, str, sequence, or Mapping
        The maximum ET flux, in model length/time. A number is applied to every
        selected cell; a ``str`` is passed through unchanged, as an MF6 time-series
        name; a sequence must hold exactly one value per selected cell, in
        selection order; a cell mapping is looked up by ``(layer, cell)`` and then
        by bare cell number. A mapping whose keys are ALL integers in
        ``range(nper)`` is read as per-period (each value any of the other shapes),
        and a period it leaves out raises ``ValueError``. A missing cell, a
        wrong-length sequence, or a negative number raises ``ValueError`` at
        :meth:`build`.
    depth : float, str, sequence, or Mapping
        The extinction depth below ``surface``, in model length units, at which ET
        falls to zero. The same shapes and checks as ``rate``, negative numbers
        refused.
    surface : CellSurfaceOffset, float, str, sequence, or Mapping, optional
        The ET surface elevation per cell, the same in every period (no per-period
        form). Defaults to ``CellSurfaceOffset("model_top")``, the land surface. A
        :class:`~myflopy.geopackage.CellSurfaceOffset` is resolved per
        selected cell from the context's elevations, and its ``offset``/``minimum``
        must be numbers: a column name raises ``ValueError`` here (use
        ``mf.evt.gpkg(...)`` for that). A number is used everywhere, a ``str`` is
        passed through as an MF6 time-series name, a sequence needs one value per
        selected cell, and a mapping is looked up by ``(layer, cell)`` then bare
        cell number -- a missing cell raises ``ValueError`` at :meth:`build`.
    """

    rate: Any
    depth: Any
    surface: Any = CellSurfaceOffset("model_top")
    name: str = "evt"

    def __post_init__(self) -> None:
        """Validate the base contract, plus that a ``CellSurfaceOffset`` surface is numeric."""

        super().__post_init__()
        surface = self.surface
        if isinstance(surface, CellSurfaceOffset) and (
            isinstance(surface.offset, str) or isinstance(surface.minimum, str)
        ):
            raise ValueError(
                "file-less mf.evt(context=...) requires a numeric CellSurfaceOffset "
                "offset/minimum (there are no feature rows to name a field); use "
                "mf.evt.gpkg(...) for field-driven surfaces."
            )

    def _resolve_surface(self, cells: tuple[CellId, ...]) -> dict[CellId, Any]:
        """Resolve the ET surface elevation for each selected cell.

        Delegates ``CellSurfaceOffset`` to the shared :class:`SurfaceResolver`
        (the same engine the ``mf.evt.gpkg`` path uses), so the ``gdf_topbtm``
        column layout is never re-encoded here.
        """

        surface = self.surface
        if isinstance(surface, CellSurfaceOffset):
            resolver = SurfaceResolver(self.context)
            return {
                (layer, cell): surface.resolve(resolver, {}, layer=layer, cell=cell)
                for (layer, cell) in cells
            }
        if isinstance(surface, Real):
            return {cellid: float(surface) for cellid in cells}
        return {
            cellid: self._value_for_cell(surface, cellid, cells, label="surface")
            for cellid in cells
        }

    def _prepared(self, cells: tuple[CellId, ...]) -> dict[CellId, Any]:
        """Resolve the (period-independent) per-cell ET surface once per build."""

        return self._resolve_surface(cells)

    def _row_for(self, period: int, cellid: CellId, cells: tuple[CellId, ...], prepared: Any) -> list[Any]:
        """One EVT record ``[cellid, surface, rate, depth]`` (rate and depth nonnegative)."""

        rate = self._value_for_period_cell(self.rate, period, cellid, cells, label="rate")
        depth = self._value_for_period_cell(self.depth, period, cellid, cells, label="depth")
        if isinstance(rate, Real) and float(rate) < 0.0:
            raise ValueError("rate must be nonnegative.")
        if isinstance(depth, Real) and float(depth) < 0.0:
            raise ValueError("depth must be nonnegative.")
        return [cellid, prepared[cellid], rate, depth]

    def _make_spec(self, stress_period_data: dict[int, list[list[Any]]]) -> PackageSpec:
        """Build the list-based single-segment (``nseg=1``) EVT spec."""

        return evt_spec(
            stress_period_data,
            name=self.name,
            nseg=1,
            boundnames=self.boundnames,
            **dict(self.options),
        )

    def _build_metadata(self) -> dict[str, Any]:
        """Manifest metadata for the built EVT spec."""

        return {
            "builder": "EVTBuilder",
            "cells": list(self.selected_cells),
            "surface_form": type(self.surface).__name__,
        }
