"""Wrappers for common MF6 discretization and temporal-discretization setup."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import flopy

from myflopy._logging import get_logger

logger = get_logger(__name__)

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def _first_timestep_length(perlen: float, nstp, tsmult: float) -> float:
    """Length of the *first* sub-step MF6's TDIS would take for one period.

    Mirrors MF6's geometric time-step series so an ATS period can start where
    the equivalent fixed-step period would have started.
    """
    nstp = int(nstp) if nstp else 1
    if tsmult in (0, 1) or nstp <= 1:
        return perlen / nstp
    return perlen * (tsmult - 1.0) / (tsmult ** nstp - 1.0)


def _resolve_ats_periods(ats, nper: int) -> dict[int, dict]:
    """Normalise the ``ats`` argument to ``{period: overrides_dict}``.

    Stress periods are **zero-based** throughout myflopy (like FloPy's own
    ``stress_period_data`` keys). ``ats`` may be:
      * ``None``/``False`` -> no ATS (``{}``)
      * ``True``           -> every period
      * an iterable of zero-based period indices
      * a mapping ``{period: {overrides}}``
    """
    if ats is None or ats is False:
        return {}
    if ats is True:
        return {p: {} for p in range(nper)}
    if isinstance(ats, Mapping):
        periods = {int(p): (dict(ov) if ov else {}) for p, ov in ats.items()}
    else:
        periods = {int(p): {} for p in ats}
    for p in periods:
        if not (0 <= p < nper):
            raise ValueError(
                f"ats stress period {p} is out of range 0..{nper - 1} "
                f"(stress periods are zero-based)"
            )
    return periods


def _build_ats_records(
    period_data: list,
    ats_periods: dict[int, dict],
    dt0, dtmin, dtmax, dtadj, dtfailadj,
) -> list:
    """Build ATS period records with zero-based ``iperats``.

    Stress periods are zero-based throughout myflopy. FloPy's ``ats_perioddata``
    ``iperats`` column is likewise zero-based (FloPy writes ``iperats + 1`` to the
    1-based MF6 input file), so the zero-based period index passes straight
    through. Any of ``dt0/dtmin/dtmax`` left ``None`` is auto-derived from that
    period's ``[perlen, nstp, tsmult]`` record.
    """
    records = []
    for p in sorted(ats_periods):
        rec = period_data[p]
        perlen = rec[0]
        nstp = rec[1] if len(rec) > 1 else 1
        tsmult = rec[2] if len(rec) > 2 else 1.0
        ov = ats_periods[p]

        _dt0 = ov.get("dt0", dt0)
        if _dt0 is None:
            _dt0 = _first_timestep_length(perlen, nstp, tsmult)
        _dtmin = ov.get("dtmin", dtmin)
        if _dtmin is None:
            _dtmin = perlen * 1e-5
        _dtmax = ov.get("dtmax", dtmax)
        if _dtmax is None:
            _dtmax = perlen
        _dtadj = ov.get("dtadj", dtadj)
        _dtfailadj = ov.get("dtfailadj", dtfailadj)

        records.append([p, _dt0, _dtmin, _dtmax, _dtadj, _dtfailadj])
    return records


def _set_maxats(sim, n: int):
    """Set ``maxats`` on the ATS package FloPy auto-creates from ``ats_perioddata``.

    FloPy does not size the ``MAXATS`` dimension from the record list, so the
    written file would otherwise declare ``MAXATS 1`` regardless of how many
    ATS periods were supplied. Returns the ATS package (or ``None``).
    """
    for pkg in sim.sim_package_list:
        if "ats" in (getattr(pkg, "package_type", "") or "").lower():
            pkg.maxats.set_data(n)
            return pkg
    return None


class DisuGrid:
    """Create a one-layer DISU grid from a ``VoronoiGridPlus`` helper."""

    def __init__(
            self,
            vor: Vor,
            model: SimulationBase,
            top=None,
            bottom=None
    ):
        """Parameters
        ----------
        vor
            Voronoi/discretization helper providing connectivity and geometry.
        model
            Target model receiving the DISU package.
        top, bottom
            Optional top and bottom elevations. If omitted, the wrapper tries to
            read them from ``vor.gdf_topbtm``.
        """
        self.nlay = 1
        if top is None:
            try:
                top = vor.gdf_topbtm["top"].to_list()
            except ValueError:
                logger.warning('the grid carries no top surface; DISU top is unset')
        if bottom is None:
            try:
                bottom = vor.gdf_topbtm["bottom"].to_list()
            except ValueError:
                logger.warning('the grid carries no bottom surface; DISU botm is unset')
        grid_props = vor.get_disv_gridprops()
        self.disu = flopy.mf6.ModflowGwfdisu(
            model.gwf,
            vertices=grid_props["vertices"],
            cell2d=grid_props["cell2d"],
            length_units="FEET",
            top=top,
            bot=bottom,
            filename=f"{model.name}.disu",
            nvert=len(grid_props["vertices"]),
            nodes=len(grid_props["cell2d"]),
            nja=vor.nja,
            iac=vor.iac,
            ja=vor.ja,
            area=vor.get_cell_areas(),
            ihc=1,
            cl12=vor.cl12,
            hwva=vor.hwva,
            idomain=[1 for _ in range(vor.ncpl)],
        )


class DisvGrid:
    """Create a DISV grid from a ``VoronoiGridPlus`` helper."""

    def __init__(
            self,
            vor: Vor = None,
            model: SimulationBase = None,
            top=None,
            bottom=None,
            nlay=1,
            idomain=None
    ):
        """Parameters
        ----------
        vor
            Voronoi/discretization helper providing cell2d and vertex geometry.
        model
            Target model receiving the DISV package.
        top, bottom
            Top and bottom elevation inputs passed to FloPy.
        nlay
            Number of model layers.
        idomain
            Optional idomain array passed through to the DISV package.
        """
        self.nlay = nlay
        model.nlay = nlay
        vor = model.vor if vor is None else vor
        grid_props = vor.get_disv_gridprops()
        self.disv = flopy.mf6.ModflowGwfdisv(
            model.gwf,
            length_units="FEET",
            nlay=nlay,
            ncpl=grid_props['ncpl'],
            nvert=len(grid_props["vertices"]),
            vertices=grid_props['vertices'],
            cell2d=grid_props['cell2d'],
            pname='disv',
            filename=f'{model.name}.disv',
            top=top,
            botm=bottom,
            idomain=idomain
        )


class TemporalDiscretization:
    """Create the MF6 TDIS package for a model."""

    def __init__(
            self,
            model: SimulationBase,
            time_units: str = 'DAYS',
            per_len: int = 1,
            period_data: list = None,
            num_steps=10,
            multiplier=1.1,
            ats=None,
            ats_dt0: float | None = None,
            ats_dtmin: float | None = None,
            ats_dtmax: float | None = None,
            ats_dtadj: float = 2.0,
            ats_dtfailadj: float = 5.0,
    ):
        """Parameters
        ----------
        model
            Target model receiving the TDIS package.
        time_units
            MF6 time-units label.
        per_len
            Default period length used when ``period_data`` is omitted.
        period_data
            Optional explicit MF6 period-data records.
        num_steps
            Default number of timesteps per period when ``period_data`` is omitted.
        multiplier
            Default timestep multiplier when ``period_data`` is omitted.
        ats
            Enable MF6 Adaptive Time Stepping. ATS lets MF6 shrink the time step
            when the solver struggles (and retry a *failed* step at a smaller dt
            instead of giving up) and grow it back on easy stretches -- ideal when
            only a few stress periods are numerically hard (e.g. flashy inflows into
            a lake). Accepts:

              * ``None``/``False`` -- no ATS (default; unchanged behaviour).
              * ``True`` -- ATS on every period (quiet periods take big steps, hard
                ones auto-subdivide). The simplest "set it and forget it" choice.
              * an iterable of **zero-based** period indices -- ATS on only those
                periods, e.g. ``ats=[7, 9]`` for the 8th and 10th stress periods.
              * a mapping ``{period: {overrides}}`` -- per-period control,
                e.g. ``ats={7: {"dtmin": 1e-4}}``. Overrides may set any of
                ``dt0/dtmin/dtmax/dtadj/dtfailadj``.

            Stress periods are zero-based here (as everywhere in myflopy, and as
            in FloPy's own ``stress_period_data`` keys); MF6's listing prints them
            1-based, so myflopy period ``p`` is "PERIOD ``p + 1``" in the .lst.
        ats_dt0
            Initial step for ATS periods. ``None`` -> the first sub-step the
            equivalent fixed-step period would have taken.
        ats_dtmin
            Minimum allowed step for ATS periods. ``None`` -> ``perlen * 1e-5``.
        ats_dtmax
            Maximum allowed step for ATS periods. ``None`` -> ``perlen`` (a full
            period in one step when nothing is straining the solver).
        ats_dtadj
            Factor to grow/shrink the step by based on solver effort (must be 0, 1,
            or > 1; 0/1 disables growth). Default ``2.0``.
        ats_dtfailadj
            Divisor applied to retry a *failed* step at a smaller dt (must be 0 or
            > 1; 0 means a failed step stops the run). Default ``5.0``.
        """
        nper = model.nper
        if period_data is None:
            period_data = [[per_len, num_steps, multiplier] for _ in range(nper)]
        model.num_steps = num_steps
        model.per_len = per_len

        ats_periods = _resolve_ats_periods(ats, nper)
        ats_records = None
        self.ats = None
        self.ats_perioddata = None
        if ats_periods:
            ats_records = _build_ats_records(
                period_data, ats_periods,
                ats_dt0, ats_dtmin, ats_dtmax, ats_dtadj, ats_dtfailadj,
            )
            self.ats_perioddata = ats_records

        self.tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
            model.sim,
            pname="tdis",
            time_units=time_units,
            nper=nper,
            perioddata=period_data,
            ats_perioddata=ats_records,
            filename=f"{model.name}.tdis"
        )

        if ats_records:
            # FloPy auto-creates the ATS package but leaves MAXATS at 1; size it.
            self.ats = _set_maxats(model.sim, len(ats_records))
            model.ats_perioddata = self.ats_perioddata
