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

The same rule, over a well screen's open interval instead of the whole layer,
weights a screened head target (:func:`screened_head`, ledger 556):
both read their thickness from ONE function, :func:`saturated_thickness`.

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
    "FlowProperties",
    "Transmissivity",
    "declared_time_unit",
    "flow_model",
    "length_label",
    "model_flow_properties",
    "saturated_thickness",
    "screened_head",
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


def saturated_thickness(heads, top, bottom, *, icelltype=None, start=None) -> np.ndarray:
    """The thickness MF6 conducts through between ``bottom`` and ``top``, by its ICELLTYPE rule.

    The ONE statement of the rule (see the module docstring), read by
    :func:`transmissivity` over whole layers and by :func:`screened_head` over a
    screen's open interval in each layer (ledger 556).

    Parameters
    ----------
    heads : float or array_like
        The head in each cell, broadcast against ``top`` and ``bottom``. A DRY
        cell's head -- MF6's ``-1e30``, or under NEWTON a finite head below the
        cell's bottom -- gives zero thickness, through the clip.
    top, bottom : array_like
        The interval's top and bottom elevations, same shape as each other: a
        layer's top and bottom for a transmissivity, or the part of a screen
        inside each layer (``open_top``/``open_bottom``) for a screened head. An
        interval with ``top <= bottom`` is zero thickness.
    icelltype : array_like or None, optional
        MF6's NPF ``ICELLTYPE`` per cell, broadcast against ``top``: ``0`` is
        confined (the FULL interval, whatever the head); anything else is
        convertible (``min(head, top) - bottom``, floored at 0, capped at the full
        interval). Default ``None``: every cell convertible.
    start : array_like or None, optional
        The STARTING heads (``IC`` ``strt``) when the NPF says ``THICKSTRT``: a cell
        with ``icelltype < 0`` then conducts through ``min(strt, top) - bottom``
        (floored and capped as above), fixed for the run. Default ``None``: no
        THICKSTRT, so ``icelltype < 0`` is convertible like any non-zero type.

    Returns
    -------
    numpy.ndarray
        The thickness, in the elevations' length units, broadcast shape of the
        inputs.
    """

    top = np.asarray(top, dtype=float)
    bottom = np.asarray(bottom, dtype=float)
    heads = np.asarray(heads, dtype=float)
    full = np.clip(top - bottom, 0.0, None)
    celltype = np.ones(np.broadcast(top, bottom).shape, dtype=int) if icelltype is None \
        else np.asarray(icelltype).astype(int)
    # A DRY cell's head is MF6's DHDRY, -1e30: below every bottom, so the clip
    # already gives it zero thickness. (A sentinel mask here was deleted as
    # unfalsifiable -- no head MF6 writes for an ACTIVE cell reaches it; an
    # inactive cell's +1e30 is masked by idomain by the callers.)
    saturated = np.clip(np.minimum(heads, top) - bottom, 0.0, full)
    thickness = np.where(celltype != 0, saturated, full)
    if start is not None:
        from_start = np.clip(np.minimum(np.asarray(start, dtype=float), top) - bottom, 0.0, full)
        thickness = np.where(celltype < 0, from_start, thickness)
    return thickness


def screened_head(heads, k, open_bottom, open_top, *, icelltype=None, fixed_thickness=None):
    """The head a well screened across several layers reads: their heads weighted by K x saturated open length.

    A screened head target is ONE observation (ledger 556): the heads
    of the layers its screen opens into, each weighted by that layer's
    horizontal K times the length of the screen inside the layer that is
    SATURATED -- :func:`saturated_thickness` over the open interval, MF6's
    ICELLTYPE rule as :func:`transmissivity` states it. The weighted mean is
    ``sum(K_i * b_i * h_i) / sum(K_i * b_i)``. When every weight is zero -- every
    open interval dry -- it is the head of the DEEPEST open layer (the last), which
    under NEWTON is a finite head below the screen that still moves with the
    parameters (on Ten Trails 12,290 of 90,446 active cells hold such a head), so
    a calibration still sees a gradient.

    Pure numpy; the layers run along the LAST axis of every argument.

    Parameters
    ----------
    heads : array_like
        The head in each open layer, ``(..., n)``: one well's ``(n,)``, or
        ``(ntimes, n)`` for one per time. Ordered shallowest open layer FIRST,
        so the last is the deepest, which the all-dry case reads. Heads are in
        the model's length units; a dry head (``-1e30``, or below the interval's
        bottom) weighs zero.
    k : array_like
        Horizontal hydraulic conductivity (NPF ``K``) of the cell in each open
        layer, ``(n,)`` or broadcast against ``heads``, in length per time; only
        ratios matter. NaN is not allowed: inactive layers are not open layers.
    open_bottom, open_top : array_like
        The screen's interval inside each open layer, ``(n,)``: the screen
        clipped to the cell (``max(screen bottom, cell bottom)``,
        ``min(screen top, cell top)``), in the model's length units.
    icelltype : array_like or None, optional
        NPF ``ICELLTYPE`` of each open layer's cell, ``(n,)``. ``0`` conducts
        through the whole open length whatever the head; non-zero through its
        saturated part. Default ``None``: every layer convertible.
    fixed_thickness : array_like or None, optional
        ``(n,)``; where finite, the open length that conducts, used in place of
        the head-dependent one -- a THICKSTRT layer's, fixed at build from the
        starting heads. Default ``None``: none fixed.

    Returns
    -------
    float or numpy.ndarray
        The screened head, ``heads.shape[:-1]`` (a float for one well at one
        time).

    Raises
    ------
    ValueError
        When the arguments do not share their last axis, or that axis is empty.

    Examples
    --------
    Two layers, K 10 and 1, each open 5 ft and saturated:

    >>> screened_head([100.0, 90.0], [10.0, 1.0], [80.0, 70.0], [85.0, 75.0],
    ...               icelltype=[0, 0])
    99.0909...
    """

    heads = np.asarray(heads, dtype=float)
    k = np.asarray(k, dtype=float)
    open_bottom = np.asarray(open_bottom, dtype=float)
    open_top = np.asarray(open_top, dtype=float)
    n = heads.shape[-1] if heads.ndim else 0
    if n == 0 or open_bottom.shape[-1:] != (n,) or open_top.shape[-1:] != (n,):
        raise ValueError(
            "screened_head(...): heads, open_bottom and open_top must share their last "
            f"axis, one entry per open layer (got {heads.shape}, {open_bottom.shape}, "
            f"{open_top.shape})."
        )
    thickness = saturated_thickness(
        heads, open_top, open_bottom,
        icelltype=None if icelltype is None else np.broadcast_to(np.asarray(icelltype), (n,)),
    )
    if fixed_thickness is not None:
        fixed = np.broadcast_to(np.asarray(fixed_thickness, dtype=float), (n,))
        thickness = np.where(np.isfinite(fixed), fixed, thickness)
    weight = k * thickness
    total = weight.sum(axis=-1)
    # A zero weight times a dry -1e30 is 0, never NaN; masked anyway so a stray
    # +inf head in a weightless layer cannot reach the sum.
    weighted = np.where(weight > 0, weight * heads, 0.0).sum(axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = weighted / total
    result = np.where(total > 0, mean, heads[..., -1])
    return float(result) if np.ndim(result) == 0 else result


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
        The model top ``(ncpl,)`` and the layer bottoms ``(nlay, ncpl)``; a DIS
        grid's ``(nrow, ncol)`` top and ``(nlay, nrow, ncol)`` bottoms are read
        row-major, and a DISU's per-node ``botm`` is one layer of every node.
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
    elif botm.ndim > 2:
        # A DIS grid's `(nlay, nrow, ncol)`, row-major -- the order the views number
        # a DIS cell `row * ncol + col` by (ledger 419). It raised "too many values
        # to unpack" (ledger 496).
        botm = botm.reshape(botm.shape[0], -1)
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
        if thickstrt and strt is None:
            raise ValueError("thickstrt=True needs strt, the starting heads it reads.")
        thickness = saturated_thickness(
            heads,
            layer_top,
            botm,
            icelltype=celltype,
            start=None if not thickstrt else _layer_cells(strt, nlay, ncpl, "strt"),
        )

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
    if npf is None and callable(getattr(model, "package", None)):
        # A run opened with `mf.load_mf6_run` loads its discretization and OC
        # alone until a package is asked for BY NAME, so `gwf.npf` was None and
        # transmissivity refused every such run, DIS included (ledger 496).
        # `model.package` loads the whole simulation once, and `gwf` is the
        # reloaded model after it.
        try:
            model.package("npf")
        except AttributeError:
            logger.debug("no NPF package on %r; transmissivity is refused", model)
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


class FlowProperties(NamedTuple):
    """What MF6's conductance rule reads off a flow model: K, elevations, cell types.

    Read once by :func:`model_flow_properties` for :func:`model_transmissivity` and
    the screened head targets (ledger 556), so the two cannot
    read a different K or a different set of elevations.

    Attributes
    ----------
    k : numpy.ndarray
        NPF ``K`` as FloPy holds it (``npf.k.array``, factors applied).
    top : numpy.ndarray
        The model top as the grid holds it (``modelgrid.top``).
    botm : numpy.ndarray
        The layer bottoms (``modelgrid.botm``).
    idomain : numpy.ndarray or None
        ``modelgrid.idomain``.
    icelltype : numpy.ndarray or None
        NPF ``ICELLTYPE``, or ``None`` when the package holds none.
    strt : numpy.ndarray or None
        IC ``strt`` when the NPF says ``THICKSTRT``, else ``None``.
    thickstrt : bool
        Whether the NPF says ``THICKSTRT``.
    """

    k: np.ndarray
    top: np.ndarray
    botm: np.ndarray
    idomain: np.ndarray | None
    icelltype: np.ndarray | None
    strt: np.ndarray | None
    thickstrt: bool


def model_flow_properties(model) -> FlowProperties:
    """K, elevations and cell types of ``model``'s flow model, as MF6's conductance reads them.

    Parameters
    ----------
    model : SimulationBase or object
        A myflopy model, or anything :func:`flow_model` accepts; its NPF, IC and
        discretization are read (a run opened with ``mf.load_mf6_run`` loads its
        simulation once for the NPF).

    Returns
    -------
    FlowProperties

    Raises
    ------
    ValueError
        When the model has no flow model or no NPF (:func:`flow_model`).
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
    return FlowProperties(
        npf.k.array, grid.top, grid.botm, grid.idomain, icelltype, strt, thickstrt
    )


def model_transmissivity(model, *, heads=None) -> Transmissivity:
    """:func:`transmissivity` for a model: its K, layers, idomain and cell types.

    ``heads`` -- ``(nlay, ncpl)`` at one time -- makes the result saturated; the
    caller reads them, because which time and how is the caller's question.
    """

    flow = model_flow_properties(model)
    return transmissivity(
        flow.k,
        flow.top,
        flow.botm,
        idomain=flow.idomain,
        heads=heads,
        icelltype=flow.icelltype,
        strt=flow.strt,
        thickstrt=flow.thickstrt,
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
