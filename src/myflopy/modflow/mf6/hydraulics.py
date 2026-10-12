"""Transmissivity of a model's layers -- the one computation every spelling shares.

``model.plot.map(type="trans")``, ``type="trans_max"`` and the nouns
``model.npf.transmissivity`` / ``.transmissivity_max`` all read their numbers
from :func:`model_transmissivity`, so the free map, the noun, its mosaic and its
table cannot disagree about what transmissivity is (ledger 308).

**Saturated is the default, and it follows MODFLOW's own rule rather than
"head minus bottom" everywhere.** A layer's ICELLTYPE decides what thickness MF6
conducts through:

* ``icelltype != 0`` (convertible): ``min(head, top) - bottom``, floored at 0 and
  capped at the full thickness -- what MF6 writes as ``DATA-SAT`` (measured equal
  to it, to the digit, on the canonical model);
* ``icelltype == 0`` (confined): the FULL thickness, whatever the head;
* ``icelltype < 0`` under ``THICKSTRT``: confined, with the thickness taken from
  the STARTING head, ``strt - bottom``, fixed for the run.

The maximum -- ``K * (top - bottom)`` -- needs no heads, so it is the one a model
that has not been run can show.

Pure numpy and duck-typed on the model, importing only the unit tables
(``myflopy._length_units``, ``myflopy._time_units``), so it sits in import layer 1 and
every plotting module can reach it.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from myflopy._length_units import mf6_length_abbreviation
from myflopy._logging import get_logger
from myflopy._time_units import MF6_TIME_UNITS

__all__ = [
    "MF6_TIME_UNITS",
    "Transmissivity",
    "declared_time_unit",
    "flow_model",
    "length_label",
    "transmissivity",
    "model_transmissivity",
    "unit_labels",
]

logger = get_logger(__name__)

class Transmissivity(NamedTuple):
    """Per-layer transmissivity and the two numbers it is the product of.

    Every array is ``(nlay, ncpl)``. ``thickness`` is the SATURATED thickness for a
    saturated result and the full layer thickness for a maximum one, so
    ``transmissivity == k * thickness`` holds cell for cell. Inactive cells
    (``idomain <= 0``) are NaN in all three.
    """

    transmissivity: np.ndarray
    k: np.ndarray
    thickness: np.ndarray
    saturated: bool


def _layer_cells(values, nlay: int, ncpl: int, label: str) -> np.ndarray:
    """``values`` broadcast to ``(nlay, ncpl)`` -- a scalar, one per layer, or one per cell."""

    array = np.asarray(values, dtype=float)
    if array.ndim == 0:
        return np.full((nlay, ncpl), float(array))
    flat = array.reshape(-1)
    if flat.size == nlay * ncpl:
        return flat.reshape(nlay, ncpl)
    if flat.size == nlay:
        return np.repeat(flat.reshape(nlay, 1), ncpl, axis=1)
    if flat.size == ncpl:
        return np.tile(flat.reshape(1, ncpl), (nlay, 1))
    raise ValueError(
        f"{label} must hold one value, one per layer ({nlay}) or one per cell "
        f"({nlay}x{ncpl}); got {array.size} values."
    )


def transmissivity(
    k,
    top,
    botm,
    *,
    idomain=None,
    heads=None,
    icelltype=None,
    strt=None,
    thickstrt: bool = False,
) -> Transmissivity:
    """Transmissivity per layer and cell, saturated when ``heads`` are given.

    Parameters
    ----------
    k
        Horizontal hydraulic conductivity (MF6's ``K``, the major axis), a scalar,
        one per layer, or ``(nlay, ncpl)``.
    top, botm
        The model top ``(ncpl,)`` and the layer bottoms ``(nlay, ncpl)``.
    idomain
        Cells with ``idomain <= 0`` come back NaN: MF6 removes a ``-1``
        (vertical pass-through) cell from horizontal flow as well as a ``0``.
    heads
        ``(nlay, ncpl)`` heads at one time. Omit for the MAXIMUM transmissivity,
        ``K * (top - bottom)``.
    icelltype, strt, thickstrt
        MF6's rule for which thickness conducts -- see the module docstring. With
        ``heads`` given and ``icelltype`` omitted, every layer is treated as
        convertible.
    """

    botm = np.asarray(botm, dtype=float)
    if botm.ndim == 1:
        botm = botm.reshape(1, -1)
    nlay, ncpl = botm.shape
    top = _layer_cells(top, 1, ncpl, "top")[0]
    layer_top = np.vstack([top.reshape(1, ncpl), botm[:-1]])
    full = np.clip(layer_top - botm, 0.0, None)
    k = _layer_cells(k, nlay, ncpl, "k")

    if heads is None:
        thickness = full
    else:
        heads = _layer_cells(heads, nlay, ncpl, "heads")
        celltype = (
            np.ones((nlay, ncpl), dtype=int)
            if icelltype is None
            else _layer_cells(icelltype, nlay, ncpl, "icelltype").astype(int)
        )
        # A DRY cell's head is MF6's DHDRY, -1e30: below every bottom, so the clip
        # already gives it zero thickness. (A sentinel mask here was deleted as
        # unfalsifiable -- no head MF6 writes for an ACTIVE cell reaches it; an
        # inactive cell's +1e30 is masked by idomain below.)
        saturated = np.clip(np.minimum(heads, layer_top) - botm, 0.0, full)
        thickness = np.where(celltype != 0, saturated, full)
        if thickstrt:
            if strt is None:
                raise ValueError("thickstrt=True needs strt, the starting heads it reads.")
            start = _layer_cells(strt, nlay, ncpl, "strt")
            from_start = np.clip(np.minimum(start, layer_top) - botm, 0.0, full)
            thickness = np.where(celltype < 0, from_start, thickness)

    values = k * thickness
    if idomain is not None:
        inactive = _layer_cells(idomain, nlay, ncpl, "idomain") <= 0
        values = np.where(inactive, np.nan, values)
        k = np.where(inactive, np.nan, k)
        thickness = np.where(inactive, np.nan, thickness)
    return Transmissivity(values, k, thickness, heads is not None)


def flow_model(model):
    """The FloPy GWF model behind ``model``, or raise naming what is missing."""

    if model is None:
        raise ValueError(
            "Transmissivity needs a MODEL -- its NPF package's K and cell type -- "
            "and a bare grid has none. Map it from the model: model.plot.map(type=...)."
        )
    gwf = getattr(model, "gwf", None)
    npf = getattr(gwf, "npf", None) if gwf is not None else None
    if npf is None:
        kind = getattr(model, "model_type", None) or type(gwf).__name__
        raise ValueError(
            "Transmissivity needs a FLOW model's NPF package (its K and cell type); "
            f"this model ({kind}) has none. Map the GWF model it is coupled to."
        )
    return gwf


def _option(package, name):
    """A FloPy package option's value, or None when it is unset."""

    data = getattr(package, name, None)
    getter = getattr(data, "get_data", None)
    return getter() if callable(getter) else data


def model_transmissivity(model, *, heads=None) -> Transmissivity:
    """:func:`transmissivity` for a model: its K, layers, idomain and cell types.

    ``heads`` -- ``(nlay, ncpl)`` at one time -- makes the result saturated; the
    caller reads them, because which time and how is the caller's question.
    """

    gwf = flow_model(model)
    grid = gwf.modelgrid
    npf = gwf.npf
    icelltype = getattr(getattr(npf, "icelltype", None), "array", None)
    strt = None
    thickstrt = bool(_option(npf, "thickstrt"))
    if thickstrt:
        ic = getattr(gwf, "ic", None)
        strt = getattr(getattr(ic, "strt", None), "array", None)
    return transmissivity(
        npf.k.array,
        grid.top,
        grid.botm,
        idomain=grid.idomain,
        heads=heads,
        icelltype=icelltype,
        strt=strt,
        thickstrt=thickstrt,
    )


def _simulation_tdis(model):
    """The TDIS of the simulation ``model`` belongs to, or None.

    A myflopy model names its simulation ``sim`` (``simulation`` on some
    wrappers); a FloPy model reaches it through ``gwf.simulation``.
    """

    for simulation in (
        getattr(model, "sim", None),
        getattr(model, "simulation", None),
        getattr(getattr(model, "gwf", None), "simulation", None),
    ):
        tdis = getattr(simulation, "tdis", None) if simulation is not None else None
        if tdis is not None:
            return tdis
    return None


def declared_time_unit(model) -> str | None:
    """The simulation's TDIS ``TIME_UNITS`` as MF6's lower-case word, or None.

    The one reader of a declared time unit: the PRT simulation's TDIS mirrors it
    and every label -- a budget's ``ft³/d``, the transmissivity hover's ``ft²/d``,
    a pathline's travel time -- abbreviates it through :data:`MF6_TIME_UNITS`
    (ledger 442).

    ``None`` means the simulation declares no unit MF6 would use: no TDIS, a TDIS
    with no ``TIME_UNITS`` (what a plain ``mf.tdis(...)`` writes -- FloPy answers
    the unset option with ``None``), a blank one, or ``unknown``, which MF6 reads
    as undefined. A word outside :data:`MF6_TIME_UNITS` -- FloPy writes any
    word it is given, ``"d"`` and ``"day"`` included -- is also ``None``, logged
    as a warning: MF6 stops on such a TDIS, so it names no unit to label with,
    and a PRT TDIS that mirrored it would stop too.

    Parameters
    ----------
    model : object
        A myflopy model (its ``sim`` or ``simulation``) or anything whose ``gwf``
        is a FloPy model (its ``gwf.simulation``); the first of those three that
        has a ``tdis`` is read. On a file-backed model reading ``sim`` performs
        the simulation load, and whatever that raises propagates -- a caller that
        must never fail (a hover) catches it.

    Returns
    -------
    str or None
        ``"seconds"``, ``"minutes"``, ``"hours"``, ``"days"`` or ``"years"``;
        ``None`` as above.
    """

    tdis = _simulation_tdis(model)
    if tdis is None:
        return None
    value = _option(tdis, "time_units")
    word = "" if value is None else str(value).strip().lower()
    if not word:
        return None
    if word not in MF6_TIME_UNITS:
        logger.warning(
            "TDIS TIME_UNITS %r is not a word MF6 reads (%s), and MF6 stops on it; "
            "read as no time unit.", word, ", ".join(MF6_TIME_UNITS),
        )
        return None
    return None if word == "unknown" else word


def unit_labels(model) -> tuple[str | None, str | None]:
    """``(length, time)`` abbreviations from the model's DIS and TDIS, or None each.

    Read, never assumed: the house head hover hardcoded ``"ft"`` (until ledger
    403), and a transmissivity label is a product of TWO units, so a wrong
    guess would be wrong twice. The map hovers read it for the head, the
    transmissivity and the footer's cell area (ledger 401).
    """

    length = length_label(getattr(model, "gwf", None))
    unit = declared_time_unit(model)
    time = MF6_TIME_UNITS[unit] if unit is not None else None
    return length, time


def length_label(flopy_model) -> str | None:
    """The abbreviation of a FloPy model's DIS ``LENGTH_UNITS`` -- ``ft``, ``m``, ``cm`` -- or None.

    The one reader of a declared LENGTH unit: :func:`unit_labels` takes its
    length half from here, and so do the cross-section axis labels, which wrote
    ``(ft)`` whatever the model said (ledger 464). ``None`` is a model with
    no discretization, a DIS naming no unit, or ``unknown`` -- MF6's undefined --
    and a label then carries no unit rather than a guessed one (ledger 434's rule).
    The word abbreviates through the one table every label reads
    (:func:`myflopy._length_units.mf6_length_abbreviation`, ledger 456).

    Parameters
    ----------
    flopy_model : flopy.mf6.MFModel or None
        The FloPy model whose ``disv``/``dis``/``disu`` package is read, in that
        order; a myflopy model's is its ``.gwf``. Anything without
        ``get_package`` (``None`` included) reads as no unit.

    Returns
    -------
    str or None
    """

    if flopy_model is None or not hasattr(flopy_model, "get_package"):
        return None
    for name in ("disv", "dis", "disu"):
        package = flopy_model.get_package(name)
        if package is not None:
            return mf6_length_abbreviation(_option(package, "length_units"))
    return None
