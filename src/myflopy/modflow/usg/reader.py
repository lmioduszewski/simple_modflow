"""The front door: read a MODFLOW-USG model into a :class:`UsgModel`."""

from __future__ import annotations

import re
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.usg._io import NameFile
from myflopy.modflow.usg.cln import read_cln
from myflopy.modflow.usg.model import UsgModel
from myflopy.modflow.usg.packages import (
    layered_report,
    read_bas,
    read_disu,
    read_ets,
    read_hfb,
    read_list_bc,
    read_lpf,
    read_oc_words,
    read_rch,
    read_sms,
)

logger = get_logger(__name__)

__all__ = ["read_usg"]

#: List boundary conditions this reader understands, and the numeric field count
#: of one record (including the leading node number).
_LIST_PACKAGES = {"CHD": 3, "DRN": 3, "GHB": 3, "RIV": 4, "WEL": 2}

#: A date written in a DISU stress-period trailing comment, e.g. ``7/31/2023``.
_DATE = re.compile(r"\b(\d{1,2})/(\d{1,2})/(\d{4})\b")


def read_usg(
    nam: str | Path,
    *,
    gsf: str | Path | None = None,
    crs: str | None = None,
    nper: int | None = None,
    require_layered: bool = True,
    read_boundaries: bool = True,
) -> UsgModel:
    """Read a MODFLOW-USG model.

    Parameters
    ----------
    nam
        Path to the MODFLOW-USG name file. Note that a model directory often
        holds several, differing in their DISU and therefore in how many stress
        periods they run -- read the one you actually intend to convert.
    gsf
        Path to the ``.gsf`` grid specification file. A USG ``DISU`` carries no
        coordinates, so without this the model has connectivity but no geometry:
        it can still be inspected, but it cannot be converted, plotted, or
        re-gridded. Defaults to the only ``.gsf`` beside the name file, if there
        is exactly one.
    crs
        Coordinate reference system for the grid. Defaults to
        :meth:`VoronoiGridPlus.from_gsf`'s own default.
    nper
        Override the number of stress periods to read. Defaults to the DISU's
        ``NPER``, which is the right answer: boundary-condition files routinely
        carry more period blocks than the model runs, and MODFLOW reads the
        first ``NPER`` and ignores the rest.
    require_layered
        Raise when the DISU is not layered. A non-layered (nested or
        ghost-node-refined) USG grid has no MODFLOW 6 DISV equivalent, so a
        conversion would be wrong rather than merely lossy. Set ``False`` to
        read such a model anyway for inspection.
    read_boundaries
        Read the boundary-condition packages. ``False`` reads only the grid,
        properties and timing, which is much faster on a model whose stress
        files run to tens of megabytes.

    Returns
    -------
    UsgModel

    Examples
    --------
    >>> usg = read_usg("flow-tt01_USE_5yr.nam", gsf="flow-tt01.gsf")
    >>> print(usg.report())
    >>> sim = usg.to_mf6("tentrails", crs="EPSG:2927")
    """

    nam = Path(nam)
    if not nam.is_file():
        raise FileNotFoundError(f"name file not found: {nam}")

    name_file = NameFile.read(nam)
    logger.info("reading MODFLOW-USG model %s (%s)", nam.name, ", ".join(name_file.package_types))

    disu = read_disu(name_file)
    layered, why = layered_report(disu)
    if not layered:
        message = (
            f"{nam.name}: this DISU is not layered -- {why}. MODFLOW 6's DISV requires "
            "one plan mesh repeated per layer, so converting this grid would misplace "
            "cells. Pass require_layered=False to read it for inspection anyway."
        )
        if require_layered:
            raise ValueError(message)
        logger.warning("%s", message)
    else:
        logger.debug("DISU is layered: %s", why)

    periods = int(nper if nper is not None else disu.nper)
    grid = _read_grid(name_file, gsf=gsf, crs=crs)

    bas = read_bas(name_file, disu)
    lpf = read_lpf(name_file, disu)
    sms = read_sms(name_file) if name_file.package("SMS") else None

    boundaries = {}
    rch = ets = None
    hfb = np.empty((0, 3))
    if read_boundaries:
        present = set(name_file.package_types)
        for ftype, n_values in _LIST_PACKAGES.items():
            if ftype in present:
                boundaries[ftype] = read_list_bc(
                    name_file,
                    ftype,
                    nper=periods,
                    n_gwf_nodes=disu.nodes,
                    n_values=n_values,
                )
        if "RCH" in present:
            rch = read_rch(name_file, nper=periods, ncpl=disu.ncpl)
        if "ETS" in present:
            ets = read_ets(name_file, nper=periods, ncpl=disu.ncpl)
        if "HFB6" in present or "HFB" in present:
            hfb = read_hfb(name_file)

    cln = read_cln(name_file, ncpl=disu.ncpl) if "CLN" in set(name_file.package_types) else None

    model = UsgModel(
        name_file=name_file,
        disu=disu,
        bas=bas,
        lpf=lpf,
        sms=sms,
        grid=grid,
        boundaries=boundaries,
        rch=rch,
        ets=ets,
        hfb=hfb,
        cln=cln,
        oc_words=read_oc_words(name_file),
        length_units=disu.length_units,
        time_units=disu.time_units,
        start_date_time=_infer_start(name_file, disu),
    )
    if grid is not None:
        # The grid came from a .gsf, which carries geometry only. The model has
        # the layer elevations, so publish them: the layer hover, the mounding
        # colorscale and the surface-aware builders all read `grid.gdf_topbtm`.
        model.attach_layers_to_grid()

    logger.info(
        "read %d layers x %d cells, %d periods, %d boundary package(s)%s",
        model.nlay,
        model.ncpl,
        periods,
        len(boundaries),
        f", CLN with {cln.nclnnds} nodes" if cln else "",
    )
    return model


def _read_grid(name_file: NameFile, *, gsf: str | Path | None, crs: str | None):
    """Build the plan grid from a ``.gsf``, or return ``None`` when there is none."""

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    if gsf is None:
        candidates = sorted(name_file.workspace.glob("*.gsf"))
        if len(candidates) == 1:
            gsf = candidates[0]
            logger.debug("using the only .gsf beside the name file: %s", gsf.name)
        elif not candidates:
            logger.warning(
                "no .gsf beside %s -- the model will have no geometry, so it cannot be "
                "converted, plotted or re-gridded. Pass gsf= to supply one.",
                name_file.path.name,
            )
            return None
        else:
            logger.warning(
                "%d .gsf files beside %s; pass gsf= to choose one",
                len(candidates),
                name_file.path.name,
            )
            return None

    kwargs = {"crs": crs} if crs else {}
    return VoronoiGridPlus.from_gsf(Path(gsf), **kwargs)


def _infer_start(name_file: NameFile, disu) -> str | None:
    """Infer a TDIS start datetime from the DISU period lines' trailing comments.

    MODFLOW reads ``PERLEN NSTP TSMULT Ss/Tr`` and ignores the rest of each
    stress-period line, so anything further along is free text -- and being free
    text, it is not necessarily *right*. The Ten Trails DISU parks a date there
    whose year is the literal constant 2023 on all 612 lines, so the sequence
    runs ``10/31/2023, 11/30/2023, 12/31/2023, 1/31/2023`` and jumps backwards
    every January.

    A date column is therefore only trusted when it is strictly increasing over
    the periods the model runs. Otherwise this returns ``None`` and says why:
    a start date invented from a broken column would mislabel every result and
    every time series drawn from it, silently.
    """

    entry = name_file.package("DISU")
    if entry is None:
        return None
    try:
        text = entry.path.read_text(errors="replace")
    except OSError as error:
        logger.debug("could not read DISU for a start date; TDIS will have none: %s", error)
        return None

    dates: list[datetime] = []
    for match in _DATE.finditer(text):
        month, day, year = (int(g) for g in match.groups())
        try:
            dates.append(datetime(year, month, day))
        except ValueError:
            continue
        if len(dates) >= disu.nper:
            break

    if len(dates) < 2:
        return None
    if any(b <= a for a, b in zip(dates, dates[1:], strict=False)):
        logger.warning(
            "%s carries a date on each stress-period line, but the dates do not increase "
            "(%s then %s) -- so they cannot give a start date. Pass start_date_time= to "
            "to_mf6() to set one.",
            entry.path.name,
            dates[0].date(),
            next(b.date() for a, b in zip(dates, dates[1:], strict=False) if b <= a),
        )
        return None

    start = dates[0] - timedelta(days=float(disu.perlen[0]))
    logger.debug("inferred start date %s from the DISU period labels", start.date())
    return start.strftime("%Y-%m-%d %H:%M:%S")
