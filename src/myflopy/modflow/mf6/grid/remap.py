"""Carry a package's per-cell arrays from one mesh onto another.

B3's C6. The overlay engine: it intersects the OLD grid's own cell polygons with the
NEW grid's and uses the resulting areas to combine old values into new cells. That is
what lets a model whose arrays came from nowhere reproducible -- a PEST-calibrated K
field, an array read out of an imported model -- move to another mesh at all, where
``with_grid`` could previously only stamp it ``unverifiable`` and let the build refuse it.

**A remap is the WORSE route whenever a source exists, and that is measured.**
Re-declaring the same polygon map on the new mesh (``mf.Polygons``, ledger 299) held to a
maximum 7.1% and a median 2.1% boundary-flux error against the fine mesh, where the best
reducer's own worst case was 12.5% and the worst was four orders of magnitude. So this is
the route for an array with NO recipe, not a general-purpose grid changer, and the
warning says so.

**What it deliberately does not serve** (ledger 300, each with its measurement):
an EXTENSIVE quantity -- a rate or a total -- because the right denominator is the OLD
cell's area rather than the new cell's covered area, and because ``FieldSpec.kind``
cannot tell a rate from a total; and an INTEGER FLAG, because averaging a cell type or a
zero-based layer index produces a number MF6 has no meaning for.
"""

from __future__ import annotations

import warnings

import numpy as np

from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package
from myflopy.modflow.mf6.grid.geometry import (
    REMAP_RULES,
    cell_overlay,
    refuse_a_different_frame,
)

logger = get_logger(__name__)

#: Per-cell arrays remapped with the plain area-weighted mean. Each is a quantity whose
#: value does not scale with the cell holding it and which may legitimately be zero or
#: negative -- a head especially -- so the geometric and harmonic rules, which take a
#: logarithm and a reciprocal, are not offered for them.
INTENSIVE_ARRAYS: frozenset[str] = frozenset({"strt", "ss", "sy"})

#: Conductivities. These are the reason ``rule=`` exists: the correct average depends on
#: the FLOW DIRECTION, which is a property of the solution and not of the field, so no
#: default can be right everywhere. Measured in MF6 on one 2:1 coarsening -- with the
#: contrast ACROSS the flow the arithmetic mean gave 264x the fine-mesh boundary flux and
#: the harmonic mean was nearly exact; with the SAME contrast rotated 90 degrees the
#: arithmetic mean was exact and the harmonic wrong by about 250x. The ordering REVERSES.
CONDUCTIVITY_ARRAYS: frozenset[str] = frozenset({"k", "k22", "k33"})

#: Integer flags, refused by name. Every rule here returns a float and two of the four
#: can return one that is in no source cell. Measured on an 8x8 -> 5x5 coarsening of an
#: ``icelltype`` contact: the mean gave 5 of 25 cells values like 0.3851 and 0.5459, and
#: even ``majority``, which invents nothing, moved the convertible AREA fraction from
#: 0.4949 to 0.4569. ``irch``/``ievt`` are worse again -- they are zero-based LAYER
#: INDICES, so an average of them is meaningless in a way nothing downstream catches.
FLAG_ARRAYS: frozenset[str] = frozenset({"icelltype", "iconvert", "irch", "ievt"})

#: What a new cell with NO old cell under it means.
UNCOVERED_MODES: tuple[str, ...] = ("nearest", "raise")

#: Ratio between the arithmetic and harmonic means above which the choice of rule is
#: reported as material. 2.0 is a factor of two in a conductivity, which is the point at
#: which a modeller would want to know rather than a threshold with physics behind it.
_SPREAD_REPORTED_ABOVE = 2.0


class RemapRefused(ValueError):
    """A package's arrays cannot be carried onto the new mesh, and why."""


def _nearest_old_value(old_polys, new_polys, values, targets, *, active=None):
    """Fill ``targets`` from the nearest old cell; returns ``(filled, worst_distance)``.

    Centroid-to-centroid through a KD-tree. Polygon distance is the more defensible
    measure and is what ledger 288 chose for snapping a FEATURE to a cell -- the two
    measured 3.07x apart on one real miss -- but a feature is one geometry against a
    grid, where this is thousands against thousands, and the answer is only ever used to
    fill a cell that has NO data at all. The distance is RETURNED rather than logged,
    because it is the number the caller's warning has to name: the existing raster fill
    computes it and throws it away.
    """

    from scipy.spatial import cKDTree

    usable = np.ones(len(old_polys), dtype=bool) if active is None else np.asarray(active, bool)
    usable &= np.isfinite(np.asarray(values, dtype=float))
    if not usable.any():
        return None, None
    old_centres = np.column_stack(
        (old_polys.geometry.centroid.x.to_numpy()[usable],
         old_polys.geometry.centroid.y.to_numpy()[usable])
    )
    new_centres = np.column_stack(
        (new_polys.geometry.centroid.x.to_numpy()[targets],
         new_polys.geometry.centroid.y.to_numpy()[targets])
    )
    distance, index = cKDTree(old_centres).query(new_centres)
    return np.asarray(values, dtype=float)[usable][index], float(np.max(distance))


def remap_arrays(
    options: dict,
    *,
    old_polys,
    new_polys,
    cellwise: tuple[str, ...],
    package: str,
    rule: str = "geometric",
    uncovered: str = "nearest",
    domain=None,
    verb: str = "with_grid(remap=True)",
) -> tuple[dict, dict]:
    """Return ``(new options, report)`` with every per-cell array on the new mesh.

    ``cellwise`` names the options that hold one value per cell, as
    :func:`myflopy.grid_identity.cellwise_options` decides it -- this does not re-derive
    that, because a per-LAYER list such as ``k=[10.0, 1.0]`` means the same thing on any
    mesh and must not be touched.

    ``domain`` is the OLD model's idomain, and it masks the WEIGHTS per layer rather than
    the values. Measured on a pinching wedge, handing an inactive old cell's value
    through as ``NaN`` left ten reached new cells ``NaN``, three of them ACTIVE on the new
    mesh -- and ``reached`` says nothing about it.
    """

    if rule not in REMAP_RULES:
        raise ValueError(
            f"rule={rule!r} is not a rule. Use one of "
            + ", ".join(repr(name) for name in REMAP_RULES) + "."
        )
    if uncovered not in UNCOVERED_MODES:
        raise ValueError(
            f"uncovered={uncovered!r} is not a mode. Use one of "
            + ", ".join(repr(name) for name in UNCOVERED_MODES) + "."
        )

    refused = sorted(set(cellwise) & FLAG_ARRAYS)
    if refused:
        raise RemapRefused(
            f"{package}: "
            + ", ".join(repr(name) for name in refused)
            + (" is an integer flag" if len(refused) == 1 else " are integer flags")
            + " held one value per cell, and combining several old cells into one new "
            "cell averages them -- there is no cell type 0.4, and irch/ievt are "
            "zero-based LAYER INDICES. Re-declare "
            + ("it" if len(refused) == 1 else "them")
            + " on the new mesh, or name this package in `grid_free=` if the value is "
            "uniform."
        )
    unknown = sorted(set(cellwise) - INTENSIVE_ARRAYS - CONDUCTIVITY_ARRAYS)
    if unknown:
        raise RemapRefused(
            f"{package}: "
            + ", ".join(repr(name) for name in unknown)
            + " holds one value per cell and this engine has no measured rule for "
            + ("it" if len(unknown) == 1 else "them")
            + ". It carries intensive arrays ("
            + ", ".join(sorted(INTENSIVE_ARRAYS))
            + ") and conductivities ("
            + ", ".join(sorted(CONDUCTIVITY_ARRAYS))
            + "). An EXTENSIVE quantity -- a rate or a total -- is deliberately not "
            "carried: its denominator is the OLD cell's area, and a rate and a total "
            "need opposite treatment while declaring the same kind (ledger 300)."
        )

    # The frame guard, before any value is touched and before `uncovered=` is consulted.
    refuse_a_different_frame(old_polys, new_polys, verb=verb)
    overlay = cell_overlay(old_polys, new_polys)

    idomain = None if domain is None else np.asarray(domain)
    moved: dict = dict(options)
    report: dict = {"filled": 0, "filled_area": 0.0, "distance": None, "spread": {}}
    new_area = new_polys.geometry.area.to_numpy()

    for name in cellwise:
        source = np.asarray(options[name], dtype=float)
        layers = source if source.ndim == 2 else source.reshape(1, -1)
        chosen = rule if name in CONDUCTIVITY_ARRAYS else "mean"
        carried = []
        for index, layer_values in enumerate(layers):
            active = _active_mask(idomain, index, overlay.n_old)
            painted, reached = overlay.reduce(layer_values, rule=chosen, active=active)
            if name in CONDUCTIVITY_ARRAYS:
                _note_spread(overlay, layer_values, active, report, name, index)
            if not reached.all():
                painted = _fill(
                    painted, reached, layer_values, overlay, old_polys, new_polys,
                    active=active, uncovered=uncovered, package=package, name=name,
                    new_area=new_area, report=report, verb=verb,
                )
            carried.append(painted)
        stacked = np.array(carried)
        moved[name] = stacked if source.ndim == 2 else stacked[0]

    _warn(report, package=package, cellwise=cellwise, rule=rule)
    return moved, report


def _active_mask(idomain, layer: int, n_old: int):
    """The old mesh's active cells in ``layer``, or ``None`` when it declares none."""

    if idomain is None:
        return None
    if idomain.ndim == 1:
        return idomain > 0
    if layer < idomain.shape[0]:
        return idomain[layer] > 0
    # Fewer idomain layers than array layers: the array is authoritative about its own
    # depth, and guessing a mask for a layer the domain does not describe would mask
    # cells that are in fact active. Left unmasked rather than assumed inactive.
    return None


def _note_spread(overlay, values, active, report, name, layer) -> None:
    """Record how far the candidate rules are apart, for the warning to name."""

    try:
        ratio, worst = overlay.rule_spread(values, active=active)
    except ValueError:
        # A non-positive contribution: the spread is undefined and the chosen rule has
        # already raised if it needed positivity. Nothing to report rather than a crash
        # inside a diagnostic.
        return
    if worst is None or ratio.size == 0:
        return
    material = int(np.count_nonzero(ratio > _SPREAD_REPORTED_ABOVE))
    if material:
        report["spread"][f"{name}[{layer}]"] = (material, float(np.nanmax(ratio)))


def _fill(painted, reached, layer_values, overlay, old_polys, new_polys, *,
          active, uncovered, package, name, new_area, report, verb):
    """Fill the new cells with NO old cell under them, per ``uncovered``."""

    targets = np.flatnonzero(~reached)
    short_area = float(new_area[targets].sum())
    if uncovered == "raise":
        raise RemapRefused(
            f"{verb}: {package}.{name} has {targets.size} new cells with no old cell "
            f"under them at all ({short_area:,.0f} in model units of area). Pass "
            "`uncovered='nearest'` to fill them from the nearest old cell, or move to a "
            "mesh the old one covers."
        )
    filled, distance = _nearest_old_value(
        old_polys, new_polys, layer_values, targets, active=active
    )
    if filled is None:
        raise RemapRefused(
            f"{verb}: {package}.{name} has {targets.size} uncovered new cells and NO "
            "usable old cell to fill them from -- every old cell is inactive in this "
            "layer or holds a non-finite value."
        )
    painted = painted.copy()
    painted[targets] = filled
    report["filled"] += int(targets.size)
    report["filled_area"] += short_area
    report["distance"] = max(report["distance"] or 0.0, distance or 0.0)
    return painted


def _warn(report, *, package: str, cellwise, rule: str) -> None:
    """Say what was interpolated and what was guessed -- in a channel that reaches a human.

    ``warnings.warn`` as well as the logger, because ``import myflopy`` installs a
    ``NullHandler`` on the package root and nothing configures logging, so
    ``logging.lastResort`` never fires and a log line reaches NOBODY -- measured at zero
    bytes on stderr (ledger 288). The existing nearest-neighbour fill for rasters reports
    through ``logger.info`` and is the anti-precedent.

    The PACKAGE and the parameters are named in the text because ``warnings.warn``
    de-duplicates on ``(text, category, module, lineno)``, and one move remaps several
    arrays from one line of code -- without them the second package's warning would be
    swallowed as a repeat of the first.
    """

    pieces = []
    if report["filled"]:
        pieces.append(
            f"{report['filled']} new cells had no old cell under them and were filled "
            f"from the nearest one, up to {report['distance']:,.0f} away, covering "
            f"{report['filled_area']:,.0f} in model units of area"
        )
    for key, (count, worst) in sorted(report["spread"].items()):
        pieces.append(
            f"{key}: the arithmetic and harmonic means differ by more than "
            f"{_SPREAD_REPORTED_ABOVE:g}x on {count} cells (worst {worst:,.1f}x), so "
            f"rule={rule!r} is a real choice there and not a formality"
        )
    if not pieces:
        return
    message = (
        f"{package}: remapped "
        + ", ".join(sorted(cellwise))
        + " onto the new mesh by area overlay. "
        + "; ".join(pieces)
        + ". Values were INTERPOLATED, not re-derived: where this package has a source "
        "(a polygon map, a GeoPackage, a raster) re-declaring it on the new mesh is "
        "measurably more faithful -- 7.1% worst-case flux error against this engine's "
        "12.5% at best (ledger 300)."
    )
    logger.warning("%s", message)
    warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())
