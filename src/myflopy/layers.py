"""Layer-centric authoring facade for model discretization.

`LayerStack` is a thin, friendly front door over :class:`~myflopy.surfaces.LayerSurfaces`.
You declare a model top, then ``add`` named layers -- each carrying its own bottom
definition and (optionally) its own ``min_thickness`` / ``pinch`` policy::

    from myflopy.layers import LayerStack, Raster, Flat, Contours

    stack = LayerStack(vor, top=Raster("ground.tif"), length_units="feet")
    stack.add("alluvium", bottom=Raster("allu_bot.tif"))
    stack.add("clay",     thickness=30, min_thickness=1, pinch="passthrough")
    stack.add("bedrock",  bottom=Contours("bedrock.gpkg", z="elev"), fill="propagate")

    result = stack.build()          # result.top / .botm / .idomain / .report()
    disv   = stack.to_disv(vor)     # ready-to-use mf.disv spec

This module adds **no new geometry logic**: it translates "top + N named layers"
into the "N+1 surfaces" list and delegates sampling, reconcile, and pinch-out to
the existing, tested `LayerSurfaces` engine.
"""

from __future__ import annotations

import numbers
import tempfile
from dataclasses import dataclass, field
from dataclasses import replace as _dc_replace
from pathlib import Path
from typing import Any

import numpy as np
import shapely.geometry as shp

from myflopy._logging import get_logger
from myflopy._optional import require
from myflopy.areas import Area
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.plotting import GridPlots, _as_linestring
from myflopy.modflow.mf6.grid.triangle import TriangleGrid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.utils.datatypes.extent import resolve_extent
from myflopy.modflow.utils.datatypes.readers import _vector_layer_name, read_shp_gpkg
from myflopy.package_api import disv as disv_spec
from myflopy.surfaces import (
    _MF6_LENGTH_UNITS,
    ExcavationRecord,
    LayerSurfaces,
    Surface,
    _content_keyed,
    _depends_on_previous,
    _excavate_layers,
    _read_excavation_file,
    _refuse_overlapping_excavations,
    _resolve_excavation,
)
from myflopy.viz import Fig, MplPicture, Picture, VtkScene, mpl_axes

logger = get_logger(__name__)

# Convenience source aliases -- the existing Surface constructors under friendlier
# names for layer authoring (NOT new implementations).
Raster = Surface.raster
Flat = Surface.flat
Contours = Surface.from_contours
Points = Surface.from_points
Array = Surface.from_array
Isopach = Surface.isopach
# Fractional subdivision: the contact a share of the way down to a target.
Toward = Surface.toward
# Surface algebra (lower/upper envelopes, clamping, zone selection).
Min = Surface.minimum
Max = Surface.maximum
Clamp = Surface.clamp
Where = Surface.where

_UNSET = object()

#: What ``section(color_by=...)`` will actually key the fill to. Unknown values
#: used to fall through to the layer-coloured branch and draw a picture that
#: silently was not the one asked for.
_SECTION_COLOR_BY = ("layer", "thickness")

#: Per-cell scalars ``grid(backend="vtk")`` can colour by, each with the
#: colormap that suits it. ``layer`` is discrete -- one flat colour per layer,
#: named in the bar -- and the rest are continuous fields. ``elevation`` is the
#: odd one out and stays for compatibility: it is a POINT array of vertex z, so
#: it ramps smoothly within a cell, where ``top`` and ``botm`` are the flat
#: per-cell contact elevations.
_VTK_GRID_CMAPS = {
    "layer": "tab10",
    "thickness": "viridis",
    "top": "terrain",
    "botm": "terrain",
    "cellid": "tab20",
    "elevation": "terrain",
}


def _section_line_points(line) -> list[tuple[float, float]]:
    """Coerce a section line to a list of ``(x, y)`` points.

    Accepts everything the ``section`` docstrings advertise: a shapely
    ``LineString`` / ``MultiLineString``, a path to a vector file (``.shp`` /
    ``.gpkg``, read with the same reader ``vor.plot.section`` uses), or a plain
    sequence of points.

    Raises here, at the call, rather than deep inside FloPy: a section is a lazy
    Picture, so an unusable line used to survive ``stack.plot.section(...)``
    untouched and fail much later on ``.axes`` -- ``TypeError: 'LineString'
    object is not iterable`` -- pointing nowhere near the argument that caused
    it.
    """

    if isinstance(line, (str, Path)):
        # Not wrapped in a try: a broad handler here would eat the argument
        # validation below it, which is the exact bug already found twice in
        # `read_gpkg` and `contour_line_segments`.
        line = _as_linestring(read_shp_gpkg(Path(line)).union_all())
    if isinstance(line, (shp.LineString, shp.MultiLineString)):
        line = list(_as_linestring(line).coords)

    try:
        points = [tuple(float(v) for v in pt)[:2] for pt in line]
    except TypeError:
        raise TypeError(
            f"a section line must be a LineString, a path to a vector file, or "
            f"a sequence of (x, y) points; got {type(line).__name__}."
        ) from None
    if len(points) < 2 or any(len(pt) != 2 for pt in points):
        raise ValueError(
            f"a section line needs at least two (x, y) points; got {points!r}."
        )
    return points


def _check_layer_name(name, layers) -> None:
    """Reject a layer name that would collide with another surface.

    Names have to be unique because they IDENTIFY things: a contact in
    ``surface_names``, and -- since the 3-D volume draws one actor per layer --
    an actor in the PyVista scene, where a repeat silently REPLACES its
    predecessor rather than adding one. ``"top"`` is reserved for the model top
    for the same reason; without this it failed much later, inside pandas, as
    ``cannot reindex on an axis with duplicate labels``.
    """

    if name == "top":
        raise ValueError(
            "'top' is reserved for the model top surface; name the layer for "
            "the material it is, e.g. 'topsoil' or 'upper_sand'."
        )
    if any(layer.name == name for layer in layers):
        raise ValueError(f"layer {name!r} already exists.")


def _coerce_surface(x) -> Surface:
    """Accept a Surface, or coerce a path -> Raster and a number -> Flat."""

    if isinstance(x, Surface):
        return x
    if isinstance(x, (str, Path)):
        return Surface.raster(x)
    if isinstance(x, (int, float)):
        return Surface.flat(float(x))
    raise TypeError(f"Expected a Surface, path, or number; got {type(x).__name__}.")


def _make_surface(bottom, thickness, fill) -> Surface:
    """Build a layer's bottom surface from ``bottom=`` or ``thickness=``."""

    if (bottom is None) == (thickness is None):
        raise ValueError("Provide exactly one of bottom= or thickness=.")
    if thickness is not None:
        if isinstance(thickness, Surface):
            if thickness.kind not in _THICKNESS_KINDS:
                # Any OTHER Surface was passed through verbatim, and a Surface
                # is an ELEVATION -- so `thickness=Flat(20)` under a top of 100
                # set the bottom to 20 (an 80-thick layer), silently, while the
                # docstring promised a thickness. A raw ndarray at least failed
                # loudly. The three kinds above really do measure from the
                # surface above, so they mean here exactly what they say.
                raise TypeError(
                    f"thickness= got a {thickness.kind!r} surface, which is an "
                    f"ELEVATION -- it would set the layer BOTTOM to those values "
                    f"rather than cut that thickness below the surface above. "
                    f"For a per-cell thickness map write "
                    f"`bottom=Isopach(<map>)`; for a constant, `thickness=20.0`; "
                    f"to place the contact itself, `bottom=<surface>`."
                )
            surface = thickness
        else:
            surface = Surface.constant_thickness(float(thickness))
    else:
        surface = _coerce_surface(bottom)
    if fill is not None:
        surface = _dc_replace(surface, fill=fill)
    return surface


#: Surface kinds that genuinely express "this far below the surface above", and
#: so mean what they say when handed to ``thickness=``. Every other kind is an
#: elevation, and was silently used as the layer BOTTOM.
_THICKNESS_KINDS = ("isopach", "constant_thickness", "offset_below")


def _split_shares(split, name: str) -> list[float]:
    """Normalize a ``split=`` argument to shares of the unit, summing to 1.

    ``None`` -> one layer; an int ``N`` -> N equal shares; a sequence -> those
    shares verbatim.
    """

    if split is None:
        return [1.0]
    if isinstance(split, bool):  # bool is an int; `split=True` means nothing
        raise TypeError(
            f"layer {name!r}: split= takes a count or a list of shares, not a bool."
        )
    if isinstance(split, (int, np.integer)):
        n = int(split)
        if n < 1:
            raise ValueError(f"layer {name!r}: split= must be at least 1, got {n}.")
        return [1.0 / n] * n
    # `str` before the loop: it is iterable, so it would otherwise get as far as
    # float('t') and report a confusing conversion error rather than the type.
    if isinstance(split, str):
        raise TypeError(
            f"layer {name!r}: split= takes a count (split=3) or a sequence of "
            f"shares (split=[0.3, 0.7]); got the string {split!r}."
        )
    try:
        shares = [float(f) for f in split]
    except TypeError:
        raise TypeError(
            f"layer {name!r}: split= takes a count (split=3) or a sequence of "
            f"shares (split=[0.3, 0.7]); got {type(split).__name__}."
        ) from None
    if not shares:
        raise ValueError(f"layer {name!r}: split= got an empty sequence.")
    if any(s <= 0 for s in shares):
        raise ValueError(
            f"layer {name!r}: every split share must be positive, got {shares}."
        )
    total = sum(shares)
    if abs(total - 100.0) < 1e-6:
        raise ValueError(
            f"layer {name!r}: split shares must sum to 1, not 100 -- these look "
            f"like percentages. Write {[s / 100 for s in shares]} instead of "
            f"{shares}."
        )
    if abs(total - 1.0) > 1e-6:
        raise ValueError(
            f"layer {name!r}: split shares must sum to 1, got {total:g} from "
            f"{shares}."
        )
    return shares


def _split_names(names, shares: list[float], unit: str) -> list[str] | None:
    """Validate an explicit ``names=`` list for a split unit, else ``None``.

    ``None`` keeps the generated ``<unit>_1 .. <unit>_N``. An explicit list has
    to name EVERY slice: a partial list would leave the rest generated, so which
    MODFLOW layer a given name refers to would depend on where the caller
    stopped counting. Validated at declaration rather than at build, so a bad
    list is reported at the ``add()`` that wrote it.
    """

    if names is None:
        return None
    # `str` before the coercion: it is iterable, so a bare name would otherwise
    # be silently accepted as a list of one-character names.
    if isinstance(names, str):
        raise TypeError(
            f"layer {unit!r}: names= takes a sequence of names, one per split "
            f"layer; got the string {names!r}."
        )
    try:
        given = list(names)
    except TypeError:
        raise TypeError(
            f"layer {unit!r}: names= takes a sequence of names, one per split "
            f"layer; got {type(names).__name__}."
        ) from None
    if len(shares) == 1:
        raise ValueError(
            f"layer {unit!r}: names= needs a split of 2 or more to name -- an "
            f"unsplit unit is already called {unit!r}. Pass split= as well, or "
            f"rename the unit itself."
        )
    if len(given) != len(shares):
        raise ValueError(
            f"layer {unit!r}: names= must have one name per split layer, got "
            f"{len(given)} name(s) for {len(shares)} layers."
        )
    for candidate in given:
        if not isinstance(candidate, str) or not candidate.strip():
            raise ValueError(
                f"layer {unit!r}: every split name must be a non-empty string, "
                f"got {given!r}."
            )
        if candidate == "top":
            raise ValueError(
                f"layer {unit!r}: 'top' is reserved for the model top surface; "
                f"name the layer for the material it is."
            )
    repeats = sorted({n for n in given if given.count(n) > 1})
    if repeats:
        raise ValueError(
            f"layer {unit!r}: split names must be distinct -- a name identifies "
            f"a contact, a DataFrame column and a 3-D scene actor -- got "
            f"repeats {repeats}."
        )
    return given


def _sub_surface(surface: Surface, share: float, step: float, unit: str) -> Surface:
    """One sub-layer's bottom, for a unit being split.

    Two cases, because a unit is declared in one of two ways. A unit given a
    ``thickness=`` is RELATIVE, and splitting it just divides that thickness --
    no interpolation, nothing to interpolate toward. A unit given an absolute
    ``bottom=`` needs :meth:`~myflopy.surfaces.Surface.toward`, whose ``step``
    is the share of the REMAINING interval (see :func:`_split_steps`).
    """

    if surface.is_relative:
        return _dc_replace(surface, value=float(surface.value) * share)
    if _depends_on_previous(surface):
        raise TypeError(
            f"layer {unit!r} cannot be split: its bottom is a {surface.kind!r} "
            f"surface, which is measured from the layer above, so each cut "
            f"would measure from the previous cut and the contacts would drift "
            f"downward rather than divide the unit. Give the unit an absolute "
            f"bottom (a raster, contours, or points), or declare it with "
            f"thickness= and split that."
        )
    return Surface.toward(surface, step)


#: A pinching slice's reconcile trigger sits this fraction UNDER its share of the
#: unit's minimum thickness. The pinch rule is STRICTLY below the minimum, so a unit
#: exactly at it is present and the engine's "at or below the trigger" must not reach
#: it; 1e-9 is far above the ~1e-15 a split's ``toward`` arithmetic leaves in a slice
#: and far below any thickness a modeller declares (ledger 377).
_UNDER_THE_MINIMUM = 1.0 - 1e-9


def _unit_thickness(thickness: np.ndarray, units: dict[str, list[int]]) -> np.ndarray:
    """Each layer's thickness replaced by its UNIT's total.

    Pinch policy is a statement about a geologic unit -- "where this clay is
    thinner than a foot, take it out" -- and a split is a numerical choice that
    should not change the answer. Evaluated per slice it does: a 2.5 ft unit
    split three ways comes back ``idomain [0, 1, 0]``, holes inside a unit that
    is fully present, and the zeros BLOCK vertical flow. Feeding the unit total
    in makes the verdict all-or-nothing.

    Note dividing ``min_thickness`` by N instead does NOT work -- reconcile has
    already crushed the outer slices by then, so the verdict is the same ragged
    one, and it breaks the ``min_sep < min_thickness`` invariant past N ~ 9.
    """

    if all(len(idx) == 1 for idx in units.values()):
        return thickness
    totals = np.array(thickness, dtype=float)
    for indices in units.values():
        if len(indices) > 1:
            totals[indices] = thickness[indices].sum(axis=0)
    return totals


def _unit_names(name: str, split, names) -> set[str]:
    """The unit's own name and every model-layer name it builds (``<name>_i`` or ``names=``).

    What an excavation's name must not collide with, whichever side comes first: the
    excavation is declared after the unit, or the unit is added, inserted or
    re-split after the excavation.
    """

    shares = _split_shares(split, name)
    given = _split_names(names, shares, name)
    if given:
        return {str(name), *(str(n) for n in given)}
    if len(shares) > 1:
        return {str(name), *(f"{name}_{i}" for i in range(1, len(shares) + 1))}
    return {str(name)}


def _absorbed_cells(excavated, absorbed) -> np.ndarray:
    """Which excavated cells the sliver rule ABSORBED, as opposed to cut above the bed.

    An absorbed and a cut cell both end at ``botm == bed``, so only ``absorbed`` tells
    them apart. It holds each absorbed cell's OWN thickness under the bed, and the
    first absorbed cell of a column always has some (it lay directly under the bed,
    ``bed > botm``). Absorption then runs contiguously down from it, so every excavated
    cell at or below that one was absorbed -- a zero-thickness cell among them too
    (``reconcile=False`` leaves those), whose own share is 0.0 and which a plain
    ``absorbed > 0`` test would count as cut.
    """

    started = np.cumsum(np.asarray(absorbed) > 0, axis=0) > 0
    return np.asarray(excavated, dtype=bool) & started


def _split_steps(shares: list[float]) -> list[float]:
    """Turn shares of the WHOLE unit into shares of what is left at each cut.

    ``Surface.toward`` resolves against the contact immediately above it, not
    the unit top -- so feeding it cumulative fractions 1/3, 2/3, 1 yields
    thicknesses [20, 26.67, 13.33], not equal thirds. The conversion is
    ``g_i = share_i / (1 - sum(shares before i))``. The last step is forced to
    exactly 1.0: it must land ON the unit's base, and floating-point residue in
    the running remainder would otherwise leave it fractionally short.
    """

    steps, remaining = [], 1.0
    for share in shares:
        steps.append(share / remaining)
        remaining -= share
    steps[-1] = 1.0
    return steps


def _reconcile_args(reconcile) -> tuple[bool, str]:
    """Map the facade ``reconcile`` argument to ``LayerSurfaces.sample`` kwargs."""

    if reconcile in (False, None):
        return False, "bottom"
    if reconcile is True:
        return True, "bottom"
    if reconcile in ("bottom", "top"):
        return True, reconcile
    raise ValueError("reconcile must be 'bottom', 'top', or False.")


def modflow_surfaces(source, *, resample: bool = True):
    """Read an existing MODFLOW model's surfaces as :class:`Surface` objects.

    ``source`` is a flopy model (anything exposing ``.modelgrid``) or a flopy
    modelgrid directly. Returns ``(top_surface, [bottom_surface, ...])``. With
    ``resample=True`` (default) each surface is interpolated from the source cell
    centres, so it transfers onto a *different* grid; with ``resample=False`` the
    arrays are used verbatim (the target grid must match cell-for-cell).

    Parameters
    ----------
    source : flopy model or flopy modelgrid
        Anything with a ``.modelgrid`` attribute (a FloPy model), or a FloPy
        modelgrid itself. Its ``top`` is flattened to one value per column and its
        ``botm`` to ``(nlay, ncpl)``, which fits a structured (DIS) or vertex (DISV)
        grid; its ``xcellcenters``/``ycellcenters`` are read when ``resample``.
    resample : bool, optional
        Default True: each surface is a :meth:`Surface.from_points` over the source
        cell centres, interpolated linearly (``scipy.interpolate.griddata``) at the
        target grid's centroids when sampled -- a target cell outside the convex
        hull of those centres samples NaN. False: each is a
        :meth:`Surface.from_array` of the raw per-cell values, valid only on a grid
        with exactly the source's cells. Both kinds hold per-cell data, so a stack
        built from them cannot be written by :meth:`LayerStack.to_dict`.
    """
    mg = getattr(source, "modelgrid", source)
    top = np.asarray(mg.top, dtype=float).ravel()
    botm = np.asarray(mg.botm, dtype=float)
    botm = botm.reshape(botm.shape[0], -1)
    if resample:
        xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
        yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
        top_s = Surface.from_points(xc, yc, top)
        botm_s = [Surface.from_points(xc, yc, botm[k]) for k in range(botm.shape[0])]
    else:
        top_s = Surface.from_array(top)
        botm_s = [Surface.from_array(botm[k]) for k in range(botm.shape[0])]
    return top_s, botm_s


@dataclass
class _Layer:
    name: str
    surface: Surface
    min_thickness: float | None = None
    pinch: str | None = None
    #: How many model layers this geologic UNIT becomes. `None`/1 = itself; an
    #: int = that many equal shares; a sequence = those shares. The unit stays
    #: ONE record until `_expanded_layers` compiles it, so `min_thickness` and
    #: `pinch` keep meaning "this unit", not "each slice of it".
    split: Any = None
    #: Explicit names for a split unit's layers, one per share. `None` keeps the
    #: generated `<name>_1 .. <name>_N`. Validated where it is declared, so a
    #: stack never reaches `_expanded_layers` with a names/split length mismatch.
    names: Any = None


@dataclass
class _Excavation:
    """One declared excavation of a stack: what ``LayerStack.excavate`` records.

    Nothing here is resolved: the footprint is a path or an inline ``Area`` and the
    bottom a Surface or a column name, so the declaration builds on any grid.
    """

    name: str
    #: The floor as an absolute Surface, or None when ``column`` names it.
    bottom: Surface | None
    #: A column of the footprint file giving each feature a flat floor.
    column: str | None
    #: The footprint file, resolved at declaration (ledger 293's two-bases trap).
    zone: Path | None
    #: The footprint inline, when not a file.
    area: Area | None
    #: The vector layer of ``zone``, by name (a plain ``str``).
    layer: str | None
    min_remnant: float
    #: ``{"zone": <16-char content hash>}`` once the stack becomes part of a LAZY package
    #: (``stamp_source_hashes``, ledger 293) -- the digest's stand-in for the footprint
    #: file's path. A STORED field, carried by ``to_dict``/``from_dict``: a ``from_dict``
    #: refusing it would make the stamping raise rather than hand back a copy.
    source_shas: dict = field(default_factory=dict)


#: The keys ``LayerStack.to_dict`` writes at the top level; ``from_dict`` refuses others,
#: so a payload from a later version is never read as a stack without its pond.
_STACK_KEYS = frozenset({"top", "layers", "length_units", "time_units", "excavations"})

#: The keys an excavation entry may carry. ``zone`` is the footprint FILE -- under the
#: same key ``Surface.where`` uses, so the bundler, the build-time anchor, ``drift()``'s
#: ``built_from`` and the hash stamping find it with no per-type wiring
#: (``specs.SURFACE_INPUT_KEYS``) -- and ``area`` is an inline footprint.
_EXCAVATION_KEYS = frozenset(
    {"name", "bottom", "zone", "area", "layer", "min_remnant", "source_shas"}
)


def _excavation_payload(excavation: _Excavation) -> dict[str, Any]:
    """One excavation as plain JSON types: the DECLARATION, resolved on no grid."""

    entry: dict[str, Any] = {"name": excavation.name}
    entry["bottom"] = (
        excavation.bottom.to_dict() if excavation.bottom is not None
        else {"column": excavation.column}
    )
    if excavation.zone is not None:
        entry["zone"] = excavation.zone.as_posix()
    else:
        entry["area"] = excavation.area.to_dict()
    if excavation.layer is not None:
        entry["layer"] = excavation.layer
    entry["min_remnant"] = excavation.min_remnant
    if excavation.source_shas:
        entry["source_shas"] = dict(excavation.source_shas)
    return entry


def _excavation_from_payload(entry: dict[str, Any], *, taken, reserved) -> _Excavation:
    """The inverse of :func:`_excavation_payload`, built DIRECTLY -- and refused as ``excavate`` refuses.

    Never through :meth:`LayerStack.excavate`, whose declaration reads the footprint
    file: a bundled payload names a project-RELATIVE path that the build anchors to the
    project root, and checking it here would check it against the working directory --
    ledger 293's two-bases trap. Every refusal that reads no file still runs, through
    the same :func:`_checked_scalars` / :func:`_refuse_declaration` pair the verb uses:
    a hand-written payload (``SimulationSpec.from_yaml``) reaches this path, and
    building it unchecked let a repeated name drop one excavation's RECORD while its cut
    stayed in the arrays, and a relative bed be dug from the ground.
    """

    unknown = sorted(set(entry) - _EXCAVATION_KEYS)
    if unknown:
        raise ValueError(
            f"an excavation entry carries keys this version does not know: {unknown}. "
            f"Known: {sorted(_EXCAVATION_KEYS)}."
        )
    if ("zone" in entry) == ("area" in entry):
        raise ValueError(
            f"excavation {entry.get('name')!r} must name its footprint once -- a file "
            "under 'zone' or an inline 'area' -- not both, and not neither."
        )
    min_remnant, layer = _checked_scalars(entry["min_remnant"], entry.get("layer"))
    bottom = entry["bottom"]
    column = surface = None
    if isinstance(bottom, dict) and set(bottom) == {"column"}:
        column = str(bottom["column"])
    else:
        surface = Surface.from_dict(bottom)
    excavation = _Excavation(
        name=str(entry["name"]),
        bottom=surface,
        column=column,
        zone=Path(entry["zone"]) if "zone" in entry else None,
        area=Area.from_dict(entry["area"]) if "area" in entry else None,
        layer=layer,
        min_remnant=min_remnant,
        source_shas=dict(entry.get("source_shas") or {}),
    )
    _refuse_declaration(excavation, taken=taken, reserved=reserved)
    return excavation


def _checked_scalars(min_remnant, layer) -> tuple[float, str | None]:
    """``min_remnant`` and ``layer`` as the declaration stores them, or the refusal.

    The half of :func:`_refuse_declaration` that must run BEFORE the footprint is read:
    the reader takes ``layer``, and an integer handed to it would read the file's first
    layer rather than be refused. Shared by ``excavate`` and ``from_dict``.
    """

    if isinstance(min_remnant, (bool, np.bool_)) or not isinstance(min_remnant, numbers.Real):
        raise TypeError(f"min_remnant must be a number (length units), not {min_remnant!r}.")
    min_remnant = float(min_remnant)
    if not np.isfinite(min_remnant) or min_remnant < 0:
        raise ValueError(f"min_remnant must be a finite number >= 0, not {min_remnant!r}.")
    if layer is not None:
        # By NAME only: in a layer stack an integer reads as a MODEL layer.
        layer = _vector_layer_name(
            layer, of="within", model_layer="not chosen by an excavation, which cuts every "
            "layer above its bed"
        )
    return min_remnant, layer


def _refuse_declaration(excavation: _Excavation, *, taken, reserved) -> None:
    """Every refusal of one excavation that needs no file read, in one place.

    Called by ``LayerStack.excavate`` / ``LayerBuildResult.excavate`` (through
    :func:`_declare_excavation`) and by ``LayerStack.from_dict`` (through
    :func:`_excavation_from_payload`), so a payload cannot build what the verb refuses.
    ``taken`` is the names already declared, ``reserved`` every unit and layer name.
    """

    from_file = excavation.zone is not None
    if excavation.layer is not None and not from_file:
        raise ValueError(
            f"layer={excavation.layer!r} selects a layer of a FILE; within= is not one. "
            "Drop layer=, or pass the path."
        )
    if excavation.column is not None and not from_file:
        raise ValueError(
            f"bottom={excavation.column!r} names a column of the footprint's attribute "
            "table, and a column floor needs within= to be a FILE -- an mf.Area or a "
            "geometry has none."
        )
    surface = excavation.bottom
    if surface is not None and _depends_on_previous(surface):
        raise ValueError(
            f"excavate(bottom=) got a {surface.kind!r} surface, which is measured from "
            "the surface ABOVE it. The stack would measure it from the ground and "
            "mf.lak from the excavated top -- a different bed (measured +5.0 ft in 70 "
            "of 70 cells). Give the floor as an absolute elevation, e.g. "
            "`ground - 9` for 9 ft below the ground surface `ground`."
        )
    if surface is not None and _propagates_inside_algebra(surface):
        raise ValueError(
            f"excavate(bottom=) got a {surface.kind!r} surface built on a "
            "fill='propagate' source. Where that source has no data it reads the "
            "surface ABOVE, so the stack would take it from the ground and mf.lak from "
            "the excavated top -- a different bed (measured: 14 of 27 cells 2.0 ft "
            "apart for `raster(p, fill='propagate') - 2`). Fill the raster's gaps, or "
            "pass the propagate source BARE: a bare one reads the ground where it has "
            "no data, so those cells are simply not cut."
        )
    name = excavation.name
    if not name.strip():
        raise ValueError(
            "an excavation's name may not be empty: it keys `result.excavations` and "
            "labels every report and refusal. Pass name=, or leave it out for the default."
        )
    if name in taken:
        raise ValueError(f"an excavation named {name!r} already exists; pass another name=.")
    if name in reserved or name == "top":
        raise ValueError(
            f"{name!r} names a unit or layer of the stack; an excavation is not a unit, "
            "so give it a name of its own (name=)."
        )


def _looks_positional(value) -> bool:
    """Whether ``value`` reads as cell numbers or a boolean mask rather than a place."""

    if isinstance(value, (range, np.ndarray)):
        return True
    if isinstance(value, (list, tuple)) and value:
        return all(isinstance(item, (numbers.Integral, bool, np.bool_)) for item in value)
    return False


def _propagates_inside_algebra(surface) -> bool:
    """Whether a ``fill="propagate"`` source sits INSIDE a composite, at any depth.

    ``fill="propagate"`` makes a source read the surface above wherever it has no data.
    A bare one is idempotent as an excavation bed -- the stack reads the ground there,
    so the cell is a rim and is not cut, and ``mf.lak`` reads the same unchanged top --
    but algebra over it is not: ``raster(p, fill="propagate") - 2`` gives the stack
    ``ground - 2`` there, which IS cut, and then gives the lake ``(ground - 2) - 2``
    off the excavated top (measured: 14 of 27 cells 2.0 ft apart).
    """

    def carries(node) -> bool:
        return node.fill == "propagate" or any(
            carries(op) for op in node.operands if op is not None
        )

    return any(carries(op) for op in surface.operands if op is not None)


def _declare_excavation(
    bottom, within, *, name, layer, min_remnant, taken, reserved, index
) -> _Excavation:
    """Validate one ``excavate(...)`` call and record it, resolving nothing on a grid.

    Every declaration-time refusal of ``LayerStack.excavate`` lives here, so the
    stack verb and ``LayerBuildResult.excavate`` refuse the same things the same way.
    The ones that read no file are :func:`_checked_scalars` and
    :func:`_refuse_declaration`, which ``LayerStack.from_dict`` runs too.
    """

    import geopandas as gpd

    min_remnant, layer = _checked_scalars(min_remnant, layer)

    zone = area = frame = None
    if isinstance(within, (str, Path)):
        path = Path(within)
        if not path.exists():
            raise FileNotFoundError(f"excavate(within=) names a file that does not exist: {path}")
        zone = path.resolve()
        # Read NOW, so a line, a zero-area feature or a missing layer refuses at the
        # call rather than at build; the read is repeated at build, which is what a
        # file edited in between gets.
        frame = _read_excavation_file(zone, layer=layer)
    else:
        if isinstance(within, Area):
            area = within
        elif isinstance(within, (gpd.GeoDataFrame, gpd.GeoSeries)) or hasattr(within, "geom_type"):
            area = Area.coerce(within)
        elif _looks_positional(within):
            raise TypeError(
                "excavate(within=) takes the footprint's POLYGONS -- a vector path, an "
                "mf.Area, a GeoDataFrame/GeoSeries or a shapely polygon -- not cell "
                "numbers or a boolean mask: those name cells by position on one mesh and "
                "would dig different cells after `with_grid`."
            )
        else:
            raise TypeError(
                f"excavate(within=) got a {type(within).__name__}; pass a vector path, "
                "an mf.Area, a GeoDataFrame/GeoSeries or a shapely polygon."
            )
        geometry = area.geometry
        if geometry.geom_type not in ("Polygon", "MultiPolygon") or not geometry.area > 0:
            raise ValueError(
                f"excavate(within=) must be polygonal with area; got a "
                f"{geometry.geom_type} of area {geometry.area:g}. A pond covers cells by "
                "AREA -- a line or point would dig the cells along it by length, which "
                "mf.lak refuses."
            )

    column = surface = None
    if isinstance(bottom, Surface):
        surface = bottom
    elif isinstance(bottom, (bool, np.bool_)):
        raise TypeError("excavate(bottom=) is an elevation, not a bool.")
    elif isinstance(bottom, numbers.Real):
        surface = Surface.flat(float(bottom))
    elif isinstance(bottom, str):
        bottom = str(bottom)
        if frame is not None and bottom in frame.columns and bottom != "geometry":
            column = bottom
        elif Path(bottom).is_file():
            surface = Surface.raster(Path(bottom))
        elif frame is None:
            raise ValueError(
                f"bottom={bottom!r} is not an existing raster file, and a column floor "
                "needs within= to be a FILE with an attribute table -- an mf.Area or a "
                "geometry has none."
            )
        else:
            names = ", ".join(str(c) for c in frame.columns if c not in ("geometry", "_source_name"))
            raise ValueError(
                f"bottom={bottom!r} is neither a column of {zone.name} ({names}) nor an "
                "existing raster file."
            )
    elif isinstance(bottom, Path):
        if not bottom.is_file():
            raise FileNotFoundError(f"excavate(bottom=) names a raster that does not exist: {bottom}")
        surface = Surface.raster(bottom)
    else:
        raise TypeError(
            f"excavate(bottom=) got a {type(bottom).__name__}; pass a Surface "
            "(mf.Raster, mf.Contours, mf.Flat, ...), a number, a raster path, or the "
            "name of a column of the footprint file."
        )
    if column is not None:
        for label, value in zip(frame.index, frame[column], strict=True):
            try:
                number = float(value)
            except (TypeError, ValueError):
                number = float("nan")
            if not np.isfinite(number):
                raise ValueError(
                    f"excavate(bottom={column!r}): feature {label!r} of {zone.name} has "
                    f"{column}={value!r}, not a number -- a NULL floor is not a floor."
                )

    if name is None:
        if zone is not None:
            name = zone.stem
        elif area is not None and area.name:
            name = area.name
        else:
            name = f"excavation_{index}"
    excavation = _Excavation(str(name), surface, column, zone, area, layer, min_remnant)
    _refuse_declaration(excavation, taken=taken, reserved=reserved)
    return excavation


@dataclass
class LayerQCReport:
    """Diagnostics for a built layer stack -- the problems to fix before MF6 runs.

    The structured result of QC-ing a :class:`LayerBuildResult` (via
    ``result.qc()`` or :meth:`LayerStack.qc`). It counts the geometry pathologies
    that make a DISV grid fail or behave oddly -- cells with no source coverage
    (NaN top/bottom), active cells bounded by a NaN surface or with non-positive
    thickness (fatal), overly thin cells, pinched-out cells, isolated active cells
    with no connection, and how many disconnected active components exist. When
    produced by ``LayerStack.qc`` it also reports how much top-down reconciliation
    moved each surface. Check :attr:`ok` for a fatal/clean verdict, or ``str(report)``
    for a human-readable per-layer breakdown.

    Cells an excavation removed (:meth:`LayerStack.excavate`) are counted on their
    own line, never as pinched or thin: MF6's idomain cannot tell an excavated 0 from
    a pinched 0, and reading a dug pond as "139 pinched cells" would send you looking
    for a geology problem that is not there.

    Parameters
    ----------
    nlay : int
        Model layers.
    ncpl : int
        Cells per layer.
    names : list of str
        The layer names, top to bottom.
    nan_top : int
        Cells with no top elevation (no source coverage).
    nan_botm : list of int
        Per layer: cells with no bottom elevation.
    nan_active_cells : int
        Active cells bounded by a NaN surface -- fatal.
    nonpositive_active : list of int
        Per layer: ACTIVE cells whose thickness is <= 0 -- fatal (MF6 refuses them).
    thin : list of int
        Per layer: cells whose UNIT is thinner than its ``min_thickness`` (a split
        unit is judged whole), excavated cells excluded.
    pinched : list of int
        Per layer: cells removed by the pinch policy (idomain != 1), excavated cells
        excluded.
    isolated_active : list of tuple of int
        ``(layer, cell)`` active cells with no active neighbour -- fatal;
        :meth:`LayerBuildResult.prune_isolated` deactivates them.
    n_active_components : int
        Connected components of active cells (pass-through bridges vertically).
    component_sizes : list of int
        Their sizes, largest first.
    reconcile_adjusted : list of int, optional
        Per layer: cells reconcile moved. Filled only by :meth:`LayerStack.qc`,
        which compares against the PRE-cut surfaces so an excavation does not read
        as reconcile moving the ground.
    reconcile_max_shift : list of float, optional
        Per layer: reconcile's largest move (length units), likewise.
    excavated : list of int, optional
        Per layer: cells an excavation removed (cut above its bed, or absorbed under
        it). ``None`` when the stack has no excavation.
    absorbed : int, optional
        How many of those were ABSORBED -- a sliver or an inactive cell directly
        under a bed, folded into the cell below (:meth:`LayerStack.excavate`,
        ``min_remnant``). ``None`` when the stack has no excavation.

    See Also
    --------
    LayerBuildResult.qc : Produces one.
    LayerBuildResult.excavation_report : Per-excavation detail.

    Examples
    --------
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> import myflopy as mf
    >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
    >>> report = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).build(vor).qc()
    >>> report.ok, report.pinched, report.excavated
    (True, [0], None)
    """

    nlay: int
    ncpl: int
    names: list[str]
    nan_top: int                       # cells with no top elevation (no coverage)
    nan_botm: list[int]                # per layer: cells with no bottom elevation
    nan_active_cells: int              # active cells bounded by a NaN surface (fatal)
    nonpositive_active: list[int]      # per layer: active cells with thickness <= 0
    thin: list[int]                    # per layer: cells thinner than min_thickness
    pinched: list[int]                 # per layer: cells removed (idomain != 1)
    isolated_active: list[tuple]       # (layer, cell) active cells with no connection
    n_active_components: int           # connected components of active cells
    component_sizes: list[int]         # sizes, largest first
    # reconcile diagnostics -- only filled by LayerStack.qc (needs the surfaces)
    reconcile_adjusted: list[int] | None = None    # per layer: cells reconcile moved
    reconcile_max_shift: list[float] | None = None  # per layer: largest move
    # excavation diagnostics -- only when the stack was excavated
    excavated: list[int] | None = None  # per layer: cells an excavation removed
    absorbed: int | None = None          # of those, cells absorbed under a bed

    @property
    def ok(self) -> bool:
        """True when nothing fatal was found (no NaN-bounded or isolated active cells,
        no non-positive active thickness). Thin/pinched cells are expected, not fatal."""
        return (
            self.nan_active_cells == 0
            and sum(self.nonpositive_active) == 0
            and len(self.isolated_active) == 0
        )

    def __str__(self) -> str:
        """Render a human-readable multi-line QC summary (status, warnings, per-layer stats)."""

        head = "OK" if self.ok else "PROBLEMS FOUND"
        lines = [
            f"LayerStack QC [{head}]: {self.nlay} layers, {self.ncpl} cells, "
            f"{self.n_active_components} active component(s)"
        ]
        if self.nan_active_cells:
            lines.append(f"  ** {self.nan_active_cells} active cell(s) bounded by a NaN "
                         "surface (no source coverage) -- fix the source or fill='propagate'")
        if self.isolated_active:
            lines.append(f"  ** {len(self.isolated_active)} isolated active cell(s) with no "
                         "connection -- prune_isolated() deactivates them")
        if sum(self.nonpositive_active):
            lines.append(f"  ** {sum(self.nonpositive_active)} active cell(s) with thickness <= 0")
        if self.excavated is not None:
            lines.append(f"  excavated: {sum(self.excavated)} cell(s) removed by excavation, "
                         f"{self.absorbed} of them absorbed under a bed (not pinched, not thin)")
        if self.n_active_components > 1:
            shown = ", ".join(str(s) for s in self.component_sizes[:5])
            lines.append(f"  note: active cells split into {self.n_active_components} components "
                         f"(sizes: {shown}{'...' if self.n_active_components > 5 else ''})")
        for i, name in enumerate(self.names):
            extra = ""
            if self.reconcile_adjusted is not None:
                extra = (f"  reconcile_moved={self.reconcile_adjusted[i]} "
                         f"(max {self.reconcile_max_shift[i]:.2f})")
            dug = "" if self.excavated is None else f" excavated={self.excavated[i]}"
            lines.append(
                f"  [{i}] {name:<14} nan_bottom={self.nan_botm[i]} "
                f"thin={self.thin[i]} pinched={self.pinched[i]}{dug} "
                f"thickness<=0(active)={self.nonpositive_active[i]}{extra}"
            )
        return "\n".join(lines)


@dataclass
class LayerBuildResult:
    """The disv-ready arrays produced by :meth:`LayerStack.build`.

    Bundles the discretization arrays a layer stack resolves to -- ``top``
    ``(ncpl,)``, ``botm`` ``(nlay, ncpl)``, ``idomain`` ``(nlay, ncpl)`` (1 active,
    -1 pass-through, 0 inactive), and the derived ``thickness`` -- alongside the
    per-layer metadata (``names``, ``min_thickness``, ``pinch`` policy, units) and a
    back-reference to the ``vor`` grid. Feed it straight to ``mf.disv(result)`` (or
    ``.top``/``.botm``/``.idomain`` to a hand-written one). Inspect quality with
    :meth:`report` (a per-layer thickness/pinch summary) or ``.qc()`` (a
    :class:`LayerQCReport`), and publish the surfaces onto the grid for the GIS-aware
    builders with :meth:`attach_to_grid` (or ``LayerStack.build(attach=True)``).

    An EXCAVATED result (:meth:`LayerStack.excavate`, :meth:`excavate`) carries the
    cut in its arrays -- ``top`` is the excavated top, the cells above each bed have
    zero thickness and idomain 0 -- plus the four fields that say what was dug:
    ``excavations``, ``ground``, ``excavated`` and ``absorbed``. An excavation is not
    a unit: ``nlay``, ``names``, ``units`` and :meth:`per_layer` are unchanged by it,
    and the thinner cell under a pond keeps its unit's K.

    Parameters
    ----------
    top : numpy.ndarray of float, shape (ncpl,)
        The model top -- AFTER any excavation (see ``ground``).
    botm : numpy.ndarray of float, shape (nlay, ncpl)
        Layer bottoms.
    idomain : numpy.ndarray of int, shape (nlay, ncpl)
        1 active, -1 vertical pass-through, 0 inactive (a pinch-out, a pruned cell,
        or an excavated one -- ``excavated`` tells them apart).
    thickness : numpy.ndarray of float, shape (nlay, ncpl)
        ``top``-of-layer minus bottom, per cell; 0.0 in an excavated cell.
    names : list of str
        The MODEL layer names, top to bottom (a split unit's slices named
        ``<unit>_1..`` unless ``names=`` gave them).
    min_thickness : list of float
        Per model layer: its unit's minimum thickness.
    pinch : list of str
        Per model layer: its unit's pinch policy (``"floor"``, ``"inactive"``,
        ``"passthrough"``).
    length_units : str
        The stack's length unit, written to DISV (``"feet"`` -> ``FEET``).
    time_units : str
        The stack's time unit.
    vor : VoronoiGridPlus, optional
        The grid the result was built on, for the views and :meth:`to_disv`.
    units : dict of str to list of int, optional
        Declared unit name -> the model layers it became, e.g.
        ``{"fill": [0], "sand": [1, 2, 3], "clay": [4]}``. Every unit appears,
        split or not, so it is also the map from what you DECLARED to what MODFLOW
        got. See :meth:`per_layer`.
    excavations : dict of str to ExcavationRecord, optional
        Every excavation cut into this result, by name: its cells and the SAMPLED
        bed at each (:class:`~myflopy.surfaces.ExcavationRecord`). What ``mf.lak``
        reads to tell an excavated lake from one nobody dug; :meth:`attach_to_grid`
        publishes it on the grid as ``grid.excavations``. Empty when nothing was dug.
    ground : numpy.ndarray of float, shape (ncpl,), optional
        The model top BEFORE the cut -- unrecoverable from the arrays after it, and
        what a depth-below-ground map, a section's ground line and a lake table's
        rim want. ``None`` when nothing was dug (``top`` is the ground then).
    excavated : numpy.ndarray of bool, shape (nlay, ncpl), optional
        True where an excavation removed the cell -- cut above a bed, or absorbed
        under one. A fact of the build: :meth:`prune_isolated` only deactivates
        active cells, so it never changes this. ``None`` when nothing was dug.
    absorbed : numpy.ndarray of float, shape (nlay, ncpl), optional
        The thickness of its OWN material each ABSORBED cell had under the bed when
        the sliver rule folded it into the cell below (``min_remnant``); 0.0
        elsewhere. Its own: in a chain -- a sliver, then a pass-through cell under
        it -- each cell counts only what lay between its own top and bottom, so a
        column's sum is the material between the bed and the connected cell's
        ORIGINAL top. The absorbed cells are the excavated ones from a column's
        first positive entry down (a zero-thickness cell absorbed below it holds
        0.0). ``None`` when nothing was dug.

    See Also
    --------
    LayerStack.build : Produces one.
    excavate : Cut a pond, pit or trench into this result.
    excavation_report : What each excavation removed and connected.

    Examples
    --------
    >>> import myflopy as mf
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
    >>> result = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).build(vor)
    >>> result.nlay, result.top[:2].tolist(), result.excavations, result.ground
    (1, [100.0, 100.0], {}, None)
    """

    top: np.ndarray          # (ncpl,)
    botm: np.ndarray         # (nlay, ncpl)
    idomain: np.ndarray      # (nlay, ncpl): 1 active, -1 pass-through, 0 inactive
    thickness: np.ndarray    # (nlay, ncpl)
    names: list[str]
    min_thickness: list[float]
    pinch: list[str]
    length_units: str
    time_units: str
    vor: Any = None          # grid the result was built on (for views)
    #: Declared unit name -> the model layers it became, e.g.
    #: ``{"fill": [0], "sand": [1, 2, 3], "clay": [4]}``. Every unit appears,
    #: split or not, so it is also the map from what you DECLARED to what MODFLOW
    #: got. See :meth:`per_layer`.
    units: dict[str, list[int]] = field(default_factory=dict)
    #: name -> ExcavationRecord for every excavation cut into this result.
    excavations: dict[str, ExcavationRecord] = field(default_factory=dict)
    #: The model top before the cut (None when nothing was dug).
    ground: np.ndarray | None = None
    #: (nlay, ncpl) True where an excavation removed the cell (None when nothing was dug).
    excavated: np.ndarray | None = None
    #: (nlay, ncpl) each absorbed cell's OWN thickness under the bed (None when nothing was dug).
    absorbed: np.ndarray | None = None

    def per_layer(self, mapping, *, default=_UNSET) -> list:
        """Expand a ``{unit: value}`` mapping to one value per MODEL layer.

        The bridge between geology and discretization. Splitting a unit changes
        ``nlay``, and every per-layer argument downstream -- ``mf.npf(k=...)``,
        ``sto``, ``ic`` -- is a POSITIONAL list, so a hand-written one silently
        means something new the moment a split is added or removed::

            mf.npf(k=layers.per_layer({"fill": 30.0, "sand": 25.0, "clay": 0.05}))

        A list of the wrong length is at least rejected by FloPy; a stale list of
        the RIGHT length is accepted with new meaning, which is the failure this
        exists to prevent. Layer names work as keys too, so one slice of a unit
        can differ from its siblings -- a later key wins.

        Parameters
        ----------
        mapping : dict
            Unit names (or individual layer names) to values.
        default : optional
            Value for units the mapping does not name. Without it, an unnamed
            unit raises.
        """

        by_layer: list = [_UNSET] * self.nlay
        for unit, indices in self.units.items():
            if unit in mapping:
                for i in indices:
                    by_layer[i] = mapping[unit]
        for i, name in enumerate(self.names):  # a layer name overrides its unit
            if name in mapping:
                by_layer[i] = mapping[name]

        unknown = set(mapping) - set(self.units) - set(self.names)
        if unknown:
            raise KeyError(
                f"{sorted(unknown)} name neither a unit nor a layer; have units "
                f"{sorted(self.units)} and layers {self.names}."
            )
        missing = [self.names[i] for i, v in enumerate(by_layer) if v is _UNSET]
        if missing:
            if default is _UNSET:
                raise KeyError(
                    f"no value for {missing}; give one per unit, or pass "
                    f"default= to fill the rest."
                )
            by_layer = [default if v is _UNSET else v for v in by_layer]
        return by_layer

    @property
    def nlay(self) -> int:
        """Number of layers (rows of ``botm``)."""

        return self.botm.shape[0]

    @property
    def n_pinched(self) -> int:
        """Cells the pinch policy (or a prune) removed: ``idomain != 1``, excavated cells excluded.

        An excavated cell is idomain 0 too, and MF6 cannot tell the two apart, so
        :attr:`excavated` is what separates them: a dug pond is not a pinch-out.

        Returns
        -------
        int
            The count over every layer.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> stack = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).add("clay", thickness=20.0)
        >>> stack.excavate(75.0, mf.Area(vor.gdf_vorPolys.geometry.loc[4])).build(vor).n_pinched
        0
        """

        removed = self.idomain != 1
        if self.excavated is not None:
            removed = removed & ~self.excavated
        return int(removed.sum())

    def _dug(self, k: int) -> np.ndarray:
        """Layer ``k``'s excavated mask (all False when nothing was dug)."""

        if self.excavated is None:
            return np.zeros(self.idomain.shape[1], dtype=bool)
        return self.excavated[k]

    def report(self) -> str:
        """Per-layer thickness, thin, pinched and excavated counts, as text.

        Thickness statistics and the ``thin`` count leave excavated cells out -- they
        are zero thick because they are lake, not because the layer thins there --
        and ``pinched`` counts only what the pinch policy removed. An excavated
        result adds an ``excavated=`` count per layer and one header clause.

        Returns
        -------
        str
            One header line, then one line per model layer.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> print(mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).build(vor).report())
        LayerStack: 1 layers, 9 cells [feet/days], 0 pinched cells
          [0] sand           thk min=20.00 mean=20.00 max=20.00  thin(<1.0)=0 pinched=0 (floor)
        """
        n_units = len(self.units) or self.nlay
        split = "" if n_units == self.nlay else f" from {n_units} units"
        dug = ""
        if self.excavated is not None:
            dug = (f", {int(self.excavated.sum())} excavated in "
                   f"{len(self.excavations)} excavation(s)")
        lines = [
            f"LayerStack: {self.nlay} layers{split}, {self.top.size} cells "
            f"[{self.length_units}/{self.time_units}], {self.n_pinched} pinched cells{dug}"
        ]
        # `thin` counts against the UNIT total, matching the pinch verdict: a
        # slice of a split unit is thinner than its unit by construction, and
        # reporting each slice as thin would flag every split stack.
        unit_thk = _unit_thickness(self.thickness, self.units)
        for i, name in enumerate(self.names):
            kept = ~self._dug(i)
            t = self.thickness[i][kept]
            thin = int(((unit_thk[i] < float(self.min_thickness[i])) & kept).sum())
            pinched = int(((self.idomain[i] != 1) & kept).sum())
            stats = (f"thk min={t.min():.2f} mean={t.mean():.2f} max={t.max():.2f}"
                     if t.size else "thk (every cell excavated)")
            extra = "" if self.excavated is None else f" excavated={int((~kept).sum())}"
            lines.append(
                f"  [{i}] {name:<14} {stats}  thin(<{self.min_thickness[i]})={thin} "
                f"pinched={pinched} ({self.pinch[i]}){extra}"
            )
        return "\n".join(lines)

    def to_disv(self, vor=None, *, name: str = "disv"):
        """Return the ``mf.disv`` spec for **this** build -- no second resolution.

        ``LayerStack.to_disv`` samples the stack again, so it is a different
        resolution from the one you are holding, and any argument the two calls do
        not share writes a model that is not the one you inspected. That drifted
        8.5 ft on a real stack, and because :meth:`attach_to_grid` publishes *this*
        build's surfaces, ``mf.CellSurfaceOffset("cell_bottom", ...)`` then placed
        drains against bottoms MODFLOW never received.

        Building the spec from the arrays already on this result removes the
        possibility by construction: there is one geometry, and it is the one you
        checked with :meth:`qc`. ``mf.disv(result)`` is this call.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            The grid to take the MESH from -- ``ncpl``, the vertices and ``cell2d``,
            through its ``get_disv_gridprops()``. Omitted, the grid this result was
            built on (:attr:`vor`). The elevations and ``idomain`` always come from
            this result, so the grid must have the result's cell count: one that does
            not RAISES (it wrote a DISV whose ``top`` FloPy then refused at build
            with an ``MFDataException`` naming neither count). A different mesh with
            the SAME count cannot be told apart and is accepted.
        name : str, default "disv"
            The spec's name -- the key ``with_packages(slot=)``, ``without_packages``
            and a run's frozen snapshot use. As on ``mf.disv``, it does not reach
            FloPy: the built package is always ``disv``, written to
            ``<model>.disv``.

        Returns
        -------
        PackageSpec
            An EAGER DISV spec: ``options`` hold ``nlay``, ``ncpl``, ``nvert``,
            ``vertices``, ``cell2d``, ``top``, ``botm``, ``idomain`` and, when the
            result declares them, ``length_units`` (``"feet"`` -> ``"FEET"``). The
            ``top``/``botm``/``idomain`` arrays are COPIES of this result's: an
            in-place edit of the result after the call -- ``result.idomain[0, :5] =
            0`` -- does not reach a package already declared, a library copy of it,
            or its digest (ledger 270, amended 2026-09-25; it did, measured, because
            ``np.asarray`` does not copy). Declare it again to take the edit.

        Raises
        ------
        ValueError
            No grid: ``vor`` omitted and the result carries none (``result.vor`` is
            ``None``). Or ``vor``'s cell count differs from the result's.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6 import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> result = mf.LayerStack(top=mf.Surface.flat(100.0)).add(
        ...     "sand", thickness=20.0).build(vor)
        >>> spec = result.to_disv()
        >>> spec.options["ncpl"], spec.options["length_units"]
        (9, 'FEET')
        >>> result.idomain[0, :3] = 0                 # an edit AFTER the declaration...
        >>> int(spec.options["idomain"].min())        # ...does not reach the package
        1
        """

        grid = self.vor if vor is None else vor
        if grid is None:
            raise ValueError("No grid; pass vor= or build the stack with a grid.")
        props = grid.get_disv_gridprops()
        ncpl = int(props["ncpl"])
        if ncpl != int(np.asarray(self.top).size):
            # FloPy took the mismatch and refused it only at BUILD, as "Unable to set
            # data 'top' layer 0. Data is not in a valid format" -- naming neither count.
            raise ValueError(
                f"This layer result has {int(np.asarray(self.top).size)} cells per layer "
                f"and the grid given has {ncpl}: its elevations and idomain belong to "
                "another mesh. Build the stack on this grid (`stack.build(vor)`), or "
                "hand `mf.disv` the STACK so it resolves on whichever grid it is given."
            )
        options = {}
        if self.length_units is not None:
            options["length_units"] = _MF6_LENGTH_UNITS.get(
                self.length_units.lower(), self.length_units.upper()
            )
        # COPIES (`np.array` copies; `np.asarray` did not). The spec is a declaration,
        # and it held this result's own arrays, so an in-place edit of the result --
        # the shape an excavation of a built result takes -- reached a package already
        # declared, a library copy and its digest (ledger 270, amended).
        return disv_spec(
            nlay=self.nlay,
            ncpl=props["ncpl"],
            nvert=len(props["vertices"]),
            vertices=props["vertices"],
            cell2d=props["cell2d"],
            top=np.array(self.top, dtype=float),
            botm=np.array(self.botm, dtype=float),
            idomain=np.array(self.idomain),
            name=name,
            **options,
        )

    def attach_to_grid(self, vor=None):
        """Publish ``top``/``botm`` onto ``vor.gdf_topbtm``, and the excavations beside it.

        Writes a GeoDataFrame with integer columns (``0`` = model top, ``1..nlay`` =
        layer bottoms), geometry first -- the format the surface-aware builders read
        (SFR reach tops, LAK lake-cell layering, ``mf.CellSurfaceOffset``). On an
        excavated result ``0`` is the EXCAVATED top, so ``CellSurfaceOffset("model_top",
        -2)`` inside a pond is the bed less 2.

        It also sets ``grid.excavations`` to this result's
        :attr:`excavations` (``{}`` when nothing was dug) -- an attribute, not a
        column, because consumers read ``gdf_topbtm`` positionally (ledger 196). The
        ``gdf_topbtm`` setter CLEARS ``grid.excavations`` and this method sets it
        after, so any other writer of the table (a direct assignment, another
        stack's attach) leaves no stale record behind: the records and the
        elevations always come from one result.

        Usually you do not call this directly: pass ``attach=True`` to
        :meth:`LayerStack.build`.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            The grid to publish on; omitted, the grid this result was built on.

        Returns
        -------
        VoronoiGridPlus
            The grid, with ``gdf_topbtm`` and ``excavations`` set.

        Raises
        ------
        ValueError
            No grid: ``vor`` omitted and the result carries none.

        Examples
        --------
        >>> import myflopy as mf, shapely
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> pond = mf.Area(vor.gdf_vorPolys.geometry.loc[4])
        >>> result = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        ...           .add("clay", thickness=20.0)
        ...           .excavate(bottom=90.0, within=pond, name="pond").build(vor))
        >>> grid = result.attach_to_grid()
        >>> float(grid.gdf_topbtm[0].loc[4]), list(grid.excavations)
        (90.0, ['pond'])
        """
        grid = self.vor if vor is None else vor
        if grid is None:
            raise ValueError("No grid to attach to; pass vor= or build the stack with a grid.")
        grid.gdf_topbtm = self._elevation_table(grid)
        # AFTER the table: its setter clears the records, so they are always the
        # ones that came with these elevations.
        grid.excavations = dict(self.excavations)
        return grid

    def _elevation_table(self, grid):
        """This result's elevations as ``gdf_topbtm``'s geometry-first table, on ``grid``.

        What :meth:`attach_to_grid` publishes on the grid, and what a lazy ``mf.disv``
        publishes on the MODEL CONTEXT -- which must be built from the result it just
        cut, not read back off the grid: with ``attach=False`` the grid still held
        whatever was there before (measured, a top of 500.0 beside a DISV top of 100.0).
        """

        import geopandas as gpd

        columns = {0: np.asarray(self.top, dtype=float)}
        for i in range(self.nlay):
            columns[i + 1] = np.asarray(self.botm[i], dtype=float)
        return gpd.GeoDataFrame(
            {"geometry": grid.gdf_vorPolys.geometry, **columns},
            geometry="geometry", crs=grid.crs,
        )

    # -- QC / validation -------------------------------------------------- #
    def _active_components(self):
        """Union-find over active cells (idomain == 1). Returns ``(labels, sizes)``:
        ``labels`` is an ``(nlay, ncpl)`` array of component roots (-1 where not
        active) and ``sizes`` maps a root to its cell count.

        Connections follow MF6: horizontal between active plan-neighbors, and
        vertical down a column where ``idomain == -1`` (pass-through) bridges
        active cells while ``idomain == 0`` (inactive) blocks them."""
        from collections import Counter

        nlay, ncpl = self.idomain.shape
        parent = list(range(nlay * ncpl))

        def find(a):
            """Union-find root of flat cell index ``a``, with path compression."""

            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            """Merge the union-find sets containing flat cell indices ``a`` and ``b``."""

            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        active = self.idomain == 1
        adj = self.vor.adjacent_cells_idx
        for k in range(nlay):
            base = k * ncpl
            for i in range(ncpl):
                if not active[k, i]:
                    continue
                for j in adj[i]:
                    if 0 <= j < ncpl and active[k, j]:
                        union(base + i, k * ncpl + j)
        for i in range(ncpl):
            last = None
            for k in range(nlay):
                d = self.idomain[k, i]
                if d == 1:
                    if last is not None:
                        union(last * ncpl + i, k * ncpl + i)
                    last = k
                elif d == 0:
                    last = None      # inactive blocks vertical flow
                # d == -1: pass-through keeps ``last`` reachable
        labels = np.full((nlay, ncpl), -1, dtype=int)
        sizes = Counter()
        for k in range(nlay):
            for i in range(ncpl):
                if active[k, i]:
                    r = find(k * ncpl + i)
                    labels[k, i] = r
                    sizes[r] += 1
        return labels, sizes

    def qc(self) -> LayerQCReport:
        """Check the built stack for problems MODFLOW 6 would choke on.

        Catches **NaN-bounded active cells** (a surface had no source coverage),
        **non-positive active thickness**, and **isolated active cells** (no
        connection to any neighbour) -- plus thin/pinched counts and the number
        of connected active components. ``LayerStack.qc`` adds reconcile
        diagnostics on top of this. Excavated cells are counted on their own
        (``excavated`` per layer, ``absorbed``), never as pinched or thin.

        Returns
        -------
        LayerQCReport
            Its ``ok`` is False if anything fatal was found.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> pond = mf.Area(vor.gdf_vorPolys.geometry.loc[4])
        >>> report = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        ...           .add("clay", thickness=20.0)
        ...           .excavate(bottom=75.0, within=pond).build(vor).qc())
        >>> report.ok, report.pinched, report.excavated, report.absorbed
        (True, [0, 0], [1, 0], 0)
        """
        nlay, ncpl = self.idomain.shape
        active = self.idomain == 1
        dug = (np.zeros((nlay, ncpl), dtype=bool) if self.excavated is None
               else np.asarray(self.excavated, dtype=bool))
        surf_top = np.vstack([self.top[None, :], self.botm[:-1]])  # top of each layer
        nan_bounded = np.isnan(surf_top) | np.isnan(self.botm)
        labels, sizes = self._active_components()
        isolated = [
            (k, i)
            for k in range(nlay)
            for i in range(ncpl)
            if active[k, i] and sizes[labels[k, i]] == 1
        ]
        unit_thk = _unit_thickness(self.thickness, self.units)
        excavated = absorbed = None
        if self.excavated is not None:
            excavated = [int(dug[k].sum()) for k in range(nlay)]
            taken = np.zeros(dug.shape) if self.absorbed is None else self.absorbed
            absorbed = int(_absorbed_cells(dug, taken).sum())
        return LayerQCReport(
            nlay=nlay,
            ncpl=ncpl,
            names=list(self.names),
            nan_top=int(np.isnan(self.top).sum()),
            nan_botm=[int(np.isnan(self.botm[k]).sum()) for k in range(nlay)],
            nan_active_cells=int((nan_bounded & active).sum()),
            nonpositive_active=[int(((self.thickness[k] <= 0) & active[k]).sum())
                                for k in range(nlay)],
            # Against the UNIT total, matching the pinch verdict: a slice of a
            # split unit is thinner than its unit by construction, so counting
            # per slice would flag every split stack as thin.
            thin=[int(((unit_thk[k] < float(self.min_thickness[k])) & ~dug[k]).sum())
                  for k in range(nlay)],
            pinched=[int(((self.idomain[k] != 1) & ~dug[k]).sum()) for k in range(nlay)],
            isolated_active=isolated,
            n_active_components=len(sizes),
            component_sizes=sorted(sizes.values(), reverse=True),
            excavated=excavated,
            absorbed=absorbed,
        )

    def validate(self) -> LayerBuildResult:
        """Raise ``ValueError`` if :meth:`qc` finds anything fatal; else return self."""
        report = self.qc()
        if not report.ok:
            raise ValueError("Layer stack failed QC:\n" + str(report))
        return self

    def prune_isolated(self) -> LayerBuildResult:
        """Deactivate (``idomain -> 0``) active cells that have no connection.

        Isolated cells have no edges, so removing them never disconnects anything
        else -- a single pass is enough. On an excavated result it then checks
        again that every excavated column still has an active cell under its bed:
        a pruned connected cell would leave a lake there connecting to nothing.

        Returns
        -------
        LayerBuildResult
            A new result (``self`` when nothing is isolated). ``excavated`` and
            ``absorbed`` are carried over unchanged: pruning only deactivates
            ACTIVE cells, so they stay facts of the build.

        Raises
        ------
        ValueError
            Pruning left an excavated column with no active cell under its bed,
            naming the excavation and the cells.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> result = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).build(vor)
        >>> result.prune_isolated() is result      # nothing isolated
        True
        """
        isolated = self.qc().isolated_active
        if not isolated:
            return self
        idomain = self.idomain.copy()
        for k, i in isolated:
            idomain[k, i] = 0
        pruned = _dc_replace(self, idomain=idomain)
        pruned._refuse_stranded_excavations("prune_isolated()")
        return pruned

    def _refuse_stranded_excavations(self, verb: str) -> None:
        """Raise if an excavated column has no active cell left under its bed."""

        stranded = {}
        ground = self.top if self.ground is None else self.ground
        for record in self.excavations.values():
            cut = [
                cell for cell, bed in zip(record.cells, record.bed, strict=True)
                if np.isfinite(bed) and bed < ground[cell]
            ]
            lost = [cell for cell in cut if not (self.idomain[:, cell] > 0).any()]
            if lost:
                stranded[record.name] = lost
        if stranded:
            detail = "; ".join(f"{name!r}: {cells[:10]}" for name, cells in stranded.items())
            raise ValueError(
                f"{verb} left an excavated column with no active cell under its bed "
                f"({detail}); a lake there would connect to nothing."
            )

    # -- excavation --------------------------------------------------------- #
    def excavate(self, bottom, within, *, name=None, layer=None, min_remnant=0.1) -> LayerBuildResult:
        """Cut a pond, basin, pit or trench into this BUILT result -- a new result.

        The same cut as :meth:`LayerStack.excavate`, through the same engine, for a
        model with no stack to declare it on: arrays loaded from a run and handed to
        ``LayerBuildResult``, or ``LayerStack.from_modflow`` (which cannot recover an
        earlier excavation). Inside ``within`` every surface above ``bottom`` is
        lowered to it, the cells above the bed go idomain 0, and a sliver or inactive
        cell directly under the bed is absorbed into the cell below -- see
        :meth:`LayerStack.excavate` for every rule. An already-excavated result can be
        cut again; the new footprint must not share a cell with an earlier one.

        Parameters
        ----------
        bottom : Surface, float, str or pathlib.Path
            The floor, an ABSOLUTE elevation in this result's ``length_units``: a
            :class:`~myflopy.surfaces.Surface` (``mf.Raster``, ``mf.Contours``,
            ``mf.Flat``, algebra on an absolute one), any real number (numpy scalars
            included), a raster path, or -- with a file ``within`` -- the name of a
            column of its attribute table giving each feature a flat floor. Sampled
            area-weighted onto the cells, with this result's ``ground`` (or ``top``)
            as the surface above. A surface measured from the one above it is refused,
            and so is algebra over a ``fill="propagate"`` source (a bare one is
            accepted: it takes the ground where it has no data, which is not cut).
        within : str, pathlib.Path, Area, geopandas.GeoDataFrame, geopandas.GeoSeries or shapely geometry
            The footprint: the polygons you pass to ``mf.lak(lakes=)``. A vector
            path (read with ``layer``), an ``mf.Area``, or a frame/series/geometry
            (coerced to an inline ``Area``; one with no CRS is in the grid's frame).
            Polygons only; cell numbers and boolean masks are refused.
        name : str, optional
            The record's key; defaults to the file stem, then ``Area.name``, then
            ``"excavation_<n>"``. Must be non-empty, and differ from every excavation
            already cut and from every unit and layer name.
        layer : str, optional
            The vector layer of a ``within`` FILE, by name; refused with any other
            ``within``, and an integer (a MODEL layer, in a layer result) is refused.
        min_remnant : float, default 0.1
            The thinnest cell allowed directly under the bed (length units); a
            thinner one is absorbed. ``0`` turns sliver absorption off.

        Returns
        -------
        LayerBuildResult
            A new result: excavated arrays, ``excavations`` gaining the new record,
            ``ground`` (kept from the first cut), ``excavated`` and ``absorbed``.
            ``self`` is unchanged.

        Raises
        ------
        ValueError
            No grid on the result; any declaration refusal of
            :meth:`LayerStack.excavate`; a footprint covering no cell; overlap with
            an existing excavation; a bed with no data anywhere in the footprint; a
            column left with no active cell under the bed.
        TypeError
            A ``bottom`` or ``within`` of an unsupported type, or a ``layer`` that
            is not a ``str``.

        See Also
        --------
        LayerStack.excavate : The declaration on a stack (and every rule, in full).
        excavation_report : What the cut removed and connected.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> result = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        ...           .add("clay", thickness=20.0).build(vor))
        >>> dug = result.excavate(90.0, mf.Area(vor.gdf_vorPolys.geometry.loc[4]), name="pond")
        >>> float(dug.top[4]), dug.idomain[:, 4].tolist(), float(result.top[4])
        (90.0, [1, 1], 100.0)
        """

        declaration = _declare_excavation(
            bottom, within, name=name, layer=layer, min_remnant=min_remnant,
            taken=set(self.excavations),
            reserved=set(self.units) | set(self.names),
            index=len(self.excavations),
        )
        return self._excavated_by([declaration])

    def _excavated_by(self, declarations, *, refresh: bool = False) -> LayerBuildResult:
        """Resolve ``declarations`` on this result's grid and cut them in: the one engine.

        ``LayerStack.build`` hands it the stack's excavations on the PRE-cut result;
        :meth:`excavate` hands it one more on an already-cut one. Every declaration is
        resolved against the same arrays and cut together, so the answer does not
        depend on declaration order (design section 3.6).
        """

        if not declarations:
            return self
        if self.vor is None:
            raise ValueError(
                "Excavating a layer result needs its grid; build it with one "
                "(`stack.build(vor)`) or pass the arrays through a result that carries it."
            )
        ground = np.array(self.top if self.ground is None else self.ground, dtype=float)
        records = [
            _resolve_excavation(
                name=declaration.name,
                bottom=declaration.bottom,
                column=declaration.column,
                within=declaration.zone if declaration.zone is not None else declaration.area,
                layer=declaration.layer,
                min_remnant=declaration.min_remnant,
                grid=self.vor,
                ground=ground,
                length_units=self.length_units,
                refresh=refresh,
            )
            for declaration in declarations
        ]
        _refuse_overlapping_excavations([*self.excavations.values(), *records])
        cut = _excavate_layers(self.top, self.botm, self.idomain, records)
        thickness = np.vstack([cut.top[None, :], cut.botm[:-1]]) - cut.botm
        excavated = cut.excavated if self.excavated is None else (self.excavated | cut.excavated)
        absorbed = cut.absorbed if self.absorbed is None else (self.absorbed + cut.absorbed)
        return _dc_replace(
            self,
            top=cut.top,
            botm=cut.botm,
            idomain=cut.idomain,
            thickness=thickness,
            excavations={**self.excavations, **{record.name: record for record in records}},
            ground=ground,
            excavated=excavated,
            absorbed=absorbed,
        )

    def excavation_report(self, *, thinner_than: float | None = None):
        """What each excavation removed, absorbed and connected -- computed now, from the arrays.

        Nothing here is stored on the :class:`~myflopy.surfaces.ExcavationRecord`
        (which holds only the declaration resolved on the grid), so it cannot go stale
        through :meth:`prune_isolated` or a copy of the result.

        **Why the connected cell matters.** MF6's vertical lake conductance is
        ``area / (1/bedleak + 0.5 * thickness / K33)``, so with a transparent bed (a
        large ``bed_leakance``) the connected cell's thickness sets the infiltration
        rate: measured on a real pond at ``bed_leakance=3000``, the half-cell term was
        98.6% of the resistance at the median connection, over cells 0.10-48.99 ft
        thick. Read ``connected_min``/``connected_max`` and ``connected_layers``
        before trusting a transparent bed; ``min_remnant`` and the layering under the
        pond are hydraulic choices there.

        Parameters
        ----------
        thinner_than : float, optional
            The threshold ``connected_thin`` counts against, in length units: a
            connected CELL thinner than it counts -- pass your build's ``min_sep`` to
            ask how many are thinner than the separation you built with. Omitted, a
            connected cell counts when its UNIT -- what is left of the unit in that
            column, a split unit judged whole -- is thinner than the unit's
            ``min_thickness``: the rule :meth:`qc` and :meth:`report` count ``thin``
            by, so the three agree. (Judged per slice, every connected slice of a
            split unit would count, being thinner than its unit by construction.)

        Returns
        -------
        pandas.DataFrame
            One row per excavation, indexed by name, columns:

            * ``cells`` -- footprint cells; ``cut`` -- of those, where the bed is
              below the ground; ``rim`` -- bed at or above the ground (not cut);
              ``no_data`` -- no bed there (not cut);
            * ``removed`` -- cells made inactive (cut above the bed, or absorbed);
              ``absorbed`` -- of those, cells absorbed under the bed;
              ``absorbed_units`` -- ``{unit: count}``; ``absorbed_volume`` -- the
              material absorbed (thickness x cell area, length units cubed);
            * ``removed_volume`` -- the void: ``(ground - top) x area`` over the
              footprint;
            * ``connected_layers`` -- ``{layer name: count}`` of the cell directly
              under each cut column; ``connected_min``/``connected_median``/
              ``connected_max`` -- its thickness; ``connected_thin`` -- how many are
              thin, by ``thinner_than``'s rule; ``kept_slivers`` -- how many are
              thinner than the excavation's ``min_remnant`` and were kept because
              nothing active lies below them.

            Empty (with those columns) when nothing was dug.

        Raises
        ------
        ValueError
            No grid on the result, when it has excavations (the volumes need cell
            areas).

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> result = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        ...           .add("clay", thickness=20.0)
        ...           .excavate(bottom=80.05, within=mf.Area(vor.gdf_vorPolys.geometry.loc[4]),
        ...                     name="pond").build(vor))
        >>> row = result.excavation_report().loc["pond"]
        >>> int(row["cut"]), int(row["absorbed"]), row["connected_layers"]
        (1, 1, {'clay': 1})
        """

        import pandas as pd

        columns = [
            "cells", "cut", "rim", "no_data", "removed", "absorbed", "absorbed_units",
            "absorbed_volume", "removed_volume", "connected_layers", "connected_min",
            "connected_median", "connected_max", "connected_thin", "kept_slivers",
        ]
        if not self.excavations:
            return pd.DataFrame(columns=columns)
        if self.vor is None:
            raise ValueError("excavation_report() needs the result's grid for cell areas.")
        area = np.asarray(self.vor.gdf_vorPolys.geometry.area, dtype=float)
        ground = self.top if self.ground is None else self.ground
        dug = np.zeros(self.idomain.shape, bool) if self.excavated is None else self.excavated
        absorbed = np.zeros(self.idomain.shape) if self.absorbed is None else self.absorbed
        unit_of = {k: unit for unit, indices in self.units.items() for k in indices}
        unit_thk = _unit_thickness(self.thickness, self.units)
        rows = {}
        for record in self.excavations.values():
            cells = np.asarray(record.cells, dtype=int)
            bed = record.bed
            finite = np.isfinite(bed)
            cut_mask = finite & (bed < ground[cells])
            cut_cells = cells[cut_mask]
            taken = _absorbed_cells(dug[:, cells], absorbed[:, cells])
            absorbed_units: dict[str, int] = {}
            for k in np.flatnonzero(taken.any(axis=1)):
                unit = unit_of.get(int(k), self.names[int(k)])
                absorbed_units[unit] = absorbed_units.get(unit, 0) + int(taken[k].sum())
            connected, thickness, thin = {}, [], 0
            kept = 0
            for cell in cut_cells:
                active = np.flatnonzero(self.idomain[:, cell] > 0)
                if not active.size:
                    continue
                k = int(active[0])
                connected[self.names[k]] = connected.get(self.names[k], 0) + 1
                t = float(self.thickness[k, cell])
                thickness.append(t)
                if thinner_than is None:   # the unit left in the column, as qc() judges it
                    thin += int(float(unit_thk[k, cell]) < float(self.min_thickness[k]))
                else:
                    thin += int(t < float(thinner_than))
                kept += int(t < record.min_remnant)
            rows[record.name] = {
                "cells": len(cells),
                "cut": int(cut_mask.sum()),
                "rim": int((finite & ~cut_mask).sum()),
                "no_data": int((~finite).sum()),
                "removed": int(dug[:, cells].sum()),
                "absorbed": int(taken.sum()),
                "absorbed_units": absorbed_units,
                "absorbed_volume": float((absorbed[:, cells] * area[cells]).sum()),
                "removed_volume": float(((ground[cells] - self.top[cells]) * area[cells]).sum()),
                "connected_layers": connected,
                "connected_min": float(np.min(thickness)) if thickness else float("nan"),
                "connected_median": float(np.median(thickness)) if thickness else float("nan"),
                "connected_max": float(np.max(thickness)) if thickness else float("nan"),
                "connected_thin": thin,
                "kept_slivers": kept,
            }
        return pd.DataFrame.from_dict(rows, orient="index")[columns]

    # -- views (thin wrappers over flopy / the surface API) --------------- #
    def _resolve_line(self, line, x, y) -> dict:
        """Build a flopy cross-section line from an explicit line, x=, y=, or a
        default West-East line through the grid centre.

        ``line`` is anything :func:`_section_line_points` accepts -- geometry,
        vector file, or points."""
        if line is not None:
            return {"line": _section_line_points(line)}
        minx, miny, maxx, maxy = (float(v) for v in self.vor.gdf_vorPolys.total_bounds)
        if x is not None:
            return {"line": [(float(x), miny), (float(x), maxy)]}
        if y is not None:
            return {"line": [(minx, float(y)), (maxx, float(y))]}
        ymid = (miny + maxy) / 2.0  # default: West-East section through the centre
        return {"line": [(minx, ymid), (maxx, ymid)]}

    def vertex_grid(self):
        """A flopy ``VertexGrid`` carrying this result's geometry."""
        import flopy

        p = self.vor.get_disv_gridprops()
        return flopy.discretization.VertexGrid(
            vertices=p["vertices"], cell2d=p["cell2d"],
            top=self.top, botm=self.botm, idomain=self.idomain,
            nlay=self.nlay, ncpl=p["ncpl"], crs=str(getattr(self.vor, "crs", None)),
        )

    def _require_finite_geometry(self, verb: str) -> None:
        """Refuse to draw a stack whose elevations contain ``NaN``.

        A NaN top or bottom propagates into the renderer's axis limits, where
        Matplotlib rejects it as ``Axis limits cannot be NaN or Inf`` -- eight
        frames deep, naming neither the layer nor the cause. The cause is always
        the same: a source surface had no data over some cells. :meth:`qc`
        already counts them per layer, so point at it.
        """

        nan_top = int(np.isnan(self.top).sum())
        nan_botm = [int(np.isnan(self.botm[k]).sum()) for k in range(self.nlay)]
        if not nan_top and not any(nan_botm):
            return
        worst = [
            f"{self.names[k]!r} ({n} cell{'s' if n != 1 else ''})"
            for k, n in enumerate(nan_botm) if n
        ]
        detail = ""
        if nan_top:
            detail += f"the model top ({nan_top} cells)"
        if worst:
            detail += (", " if detail else "") + ", ".join(worst)
        raise ValueError(
            f"cannot draw {verb}: the stack has NaN elevations in {detail}. A "
            f"surface had no source data over those cells -- a raster that does "
            f"not cover the whole grid, contours interpolated inside a smaller "
            f"hull, or a nodata value read as elevation. Run .qc() for the "
            f"per-layer counts, then either extend the source, pass "
            f"fill='propagate' on that layer to inherit the surface above, or "
            f"deactivate the cells with idomain."
        )

    def _draw_cross_section(
        self, line=None, *, x=None, y=None, color_by="layer", ax=None,
        cmap="tab10", show_grid=True, legend=True, title=None,
    ):
        """Render the layer cross-section into ``ax`` and return it.

        Private: this is the RENDERER. ``stack.plot.section(y=300)`` is the verb,
        and it returns a Picture that answers ``.show()``/``.save()``/``.html()``
        like every other. Matplotlib-native by necessity -- see
        :class:`LayerSection`.

        Layer coloring delegates to the shared
        :func:`~myflopy.modflow.mf6.cross_section_plotting.plot_layered_cross_section`
        renderer (the same core behind the model-aware ``plot_model_cross_section``)."""
        import matplotlib.pyplot as plt

        from myflopy.modflow.mf6.cross_section_plotting import (
            ModelCrossSectionStyle,
            plot_layered_cross_section,
        )

        self._require_finite_geometry("a cross-section")
        vg = self.vertex_grid()
        line_spec = self._resolve_line(line, x, y)

        if color_by == "thickness":
            import flopy

            if ax is None:
                _, ax = mpl_axes(figsize=(9, 4))
            xsec = flopy.plot.PlotCrossSection(modelgrid=vg, line=line_spec, ax=ax)
            xsec.plot_array(self.thickness, cmap="viridis")
            if show_grid:
                xsec.plot_grid(lw=0.25, color="0.3")
            ax.set_title(title or "Layer cross-section")
            ax.set_xlabel(f"distance along section [{self.length_units}]")
            ax.set_ylabel(f"elevation [{self.length_units}]")
            return ax

        cm = plt.get_cmap(cmap)
        colors = [cm(i % 10) for i in range(self.nlay)]  # matplotlib RGBA tuples
        style = ModelCrossSectionStyle(
            figsize=(9, 4), grid_linewidth=0.25, grid_color="0.3",
            layer_alpha=0.75, title_fontsize=12, label_fontsize=10,
            legend_loc="upper right", legend_frameon=True, legend_fontsize=8,
        )
        _, ax = plot_layered_cross_section(
            vg, line_spec, ax=ax, style=style,
            layer_colors=colors, layer_labels=self.names,
            show_grid=show_grid, show_layers=True, show_head=False,
            show_legend=legend, title=title or "Layer cross-section",
            xlabel=f"distance along section [{self.length_units}]",
            ylabel=f"elevation [{self.length_units}]",
        )
        return ax

    @property
    def surface_names(self) -> list[str]:
        """Plottable surface names: the model top plus each layer's bottom contact."""
        return ["top"] + list(self.names)

    def _surface_z(self, name) -> np.ndarray:
        """Elevation array for a surface name (``"top"`` or a layer's bottom)."""
        if name == "top":
            return self.top
        return self.botm[self.names.index(name)]

    def _resolve_surface_names(self, layer) -> list[str]:
        """Normalize the ``layer`` selector to a list of valid surface names."""
        valid = self.surface_names
        if isinstance(layer, str) and layer == "all":
            return valid
        names = list(layer) if isinstance(layer, (list, tuple)) else [layer]
        bad = [n for n in names if n not in valid]
        if bad:
            raise KeyError(f"no surface named {bad!r}; choose from {valid} or 'all'.")
        # Deduplicate, order preserved: drawing one contact twice means nothing,
        # and in the VTK scene each name is an ACTOR name, where a repeat
        # silently replaces the first rather than adding a second.
        return list(dict.fromkeys(names))

    @staticmethod
    def _rgb(rgba) -> str:
        """Format a 0-1 RGBA tuple as a Plotly ``rgb(r,g,b)`` string (0-255)."""

        r, g, b = (int(round(255 * c)) for c in rgba[:3])
        return f"rgb({r},{g},{b})"

    def _surface_fig(
        self, layer="top", *, resolution=120, colorscale="Earth_r",
        color_by=None, opacity=None, height=None,
    ):
        """The 3-D surface figure. Reached as ``stack.plot.surface(...)``.

        Private because it is Layer 1 only: it BUILDS the figure and nothing
        else. Writing it out and opening it are `.html(path)` / `.show()` on the
        Picture that wraps it -- this used to take `html_path=`/`browser=` and do
        both itself.

        Interactive 3D surface(s) of one or more layers (plotly).

        ``layer`` selects which surface(s) to draw and may be:

        * ``"top"`` or any layer name -- that layer's *bottom* contact (default ``"top"``),
        * a list of names, e.g. ``["top", "clay", "bedrock"]``, overlaid in one scene,
        * ``"all"`` -- the model top plus every layer bottom.

        Discover the choices with ``result.surface_names``. A single surface with
        relief is shaded by elevation (``colorscale``); a flat surface or several
        surfaces get distinct solid colors. Override with
        ``color_by="elevation"`` / ``"surface"`` and tune ``opacity``.

        ``height`` sets the figure height in pixels (default ``None`` = fill the
        container / browser window)."""

        names = self._resolve_surface_names(layer)
        crs = str(getattr(self.vor, "crs", None))
        xs = np.asarray(self.vor.centroids[0])
        ys = np.asarray(self.vor.centroids[1])
        # Every contact is masked to the grid's footprint, taken ONCE here: the
        # interpolator alone fills the convex hull of the cell centres, which
        # spans every concavity of a domain that is not convex (ledger 314).
        footprint = self.vor.gdf_vorPolys.union_all()

        # Each surface is built by the shared InterpolatedSurface.surface_trace
        # so there is one go.Surface builder across the toolkit; we keep the
        # interpolated z-array to set a common elevation scale below.
        grids = []
        for name in names:
            isurf = InterpolatedSurface(
                xs=xs, ys=ys, zs=np.asarray(self._surface_z(name)),
                surf_type="lyr", resolution=resolution, crs=crs, clip=footprint,
            )
            grids.append((name, isurf, isurf.clipped_surface))

        # Overall elevation span -> colour mode, shared colour range and z-axis.
        zmin = min(float(np.nanmin(zz)) for *_, zz in grids)
        zmax = max(float(np.nanmax(zz)) for *_, zz in grids)
        flat = zmax - zmin < 1e-6
        if color_by is None:
            color_by = "elevation" if (len(names) == 1 and not flat) else "surface"
        if opacity is None:
            opacity = 1.0 if len(names) == 1 else 0.85
        zaxis = dict(title=f"elev [{self.length_units}]")
        aspectmode = "manual"
        if flat:
            # A perfectly flat surface has zero vertical extent; a lone flat
            # go.Surface then fails to render in hardware-WebGL viewers (its faces
            # get a 0/0 colour and the trace has no z-depth), while multi-surface
            # scenes dodge this because their combined z-range is non-zero. Give
            # the sheet a faint, non-planar relief and a real z-axis so it always
            # draws; ``cube`` aspect keeps the near-flat sheet prominent.
            mid = 0.5 * (zmin + zmax)
            pad = max(1.0, abs(mid) * 0.01)
            gx0, gy0 = grids[0][1].xy_meshgrid
            xspan = float(np.nanmax(gx0) - np.nanmin(gx0)) or 1.0
            yspan = float(np.nanmax(gy0) - np.nanmin(gy0)) or 1.0
            ripple = (
                ((gx0 - np.nanmin(gx0)) / xspan - 0.5)
                + ((gy0 - np.nanmin(gy0)) / yspan - 0.5)
            ) * (0.5 * pad)                      # ~+/-pad/2 of imperceptible relief
            grids = [(n, s, zz + ripple) for (n, s, zz) in grids]
            zmin, zmax = mid - pad, mid + pad
            zaxis["range"] = [zmin, zmax]
            aspectmode = "cube"

        fig = Fig()
        if color_by == "elevation":
            for i, (name, isurf, zz) in enumerate(grids):
                fig.add_trace(isurf.surface_trace(
                    surface=zz, colorscale=colorscale, name=name,
                    cmin=zmin, cmax=zmax, opacity=opacity, showscale=(i == 0),
                    colorbar=dict(title=f"elev [{self.length_units}]"),
                ))
        else:  # one solid color per surface, distinguished by a legend
            import matplotlib.pyplot as plt

            cmap = plt.get_cmap("tab10")
            for i, (name, isurf, zz) in enumerate(grids):
                c = self._rgb(cmap(i % 10))
                fig.add_trace(isurf.surface_trace(
                    surface=zz, colorscale=[[0, c], [1, c]], name=name,
                    cmin=zmin, cmax=zmax, opacity=opacity,
                    showscale=False, showlegend=True,
                ))

        multi = len(names) > 1
        scene = dict(aspectmode=aspectmode, zaxis=zaxis)
        if aspectmode == "manual":
            scene["aspectratio"] = dict(x=1, y=0.6, z=0.45)
        fig.update_layout(
            title=("layer surfaces" if multi else f"{names[0]} surface"),
            autosize=True, height=height, margin=dict(l=0, r=0, t=40, b=0),
            showlegend=(color_by == "surface" and multi), scene=scene,
        )
        return fig

    def _resolve_layer_indices(self, layers) -> list[int]:
        """Normalize a layer selector to sorted, unique layer indices.

        Accepts ``None``/``"all"`` (every layer), a single layer name or integer
        index, a UNIT name (which expands to every layer it was split into), or a
        list mixing them."""
        if layers is None or (isinstance(layers, str) and layers == "all"):
            return list(range(self.nlay))
        items = list(layers) if isinstance(layers, (list, tuple)) else [layers]
        idx = []
        for it in items:
            if isinstance(it, (int, np.integer)) and not isinstance(it, bool):
                i = int(it)
                if not 0 <= i < self.nlay:
                    raise IndexError(f"layer index {i} out of range [0, {self.nlay}).")
                idx.append(i)
            elif it in self.names:
                idx.append(self.names.index(it))
            elif it in self.units:
                # A unit name selects every layer it became, so `grid("sand")`
                # still works once `sand` is split into `sand_1..sand_3`.
                idx.extend(self.units[it])
            else:
                choices = sorted(set(self.names) | set(self.units))
                raise KeyError(
                    f"no layer or unit named {it!r}; choose from {choices} or an index."
                )
        return sorted(set(idx))

    def _vtk_surface_plotter(
        self, layer="all", *, resolution=120, scale=8, cmap="tab10",
        opacity=1.0, width=900, height=580, show_edges=False,
    ):
        """Build the PyVista scene for one or more CONTACT SURFACES.

        Private, Layer 1 only: it builds the scene, and the
        :class:`~myflopy.viz.VtkScene` that wraps it does the showing and
        writing. Reached as ``stack.plot.surface(..., backend="vtk")``.

        The sibling of :meth:`_vtk_plotter`, and a different subject: that one
        renders the layered cell VOLUME, this renders the contacts as separate
        sheets. Sheets are what you want when the question is "where does this
        contact go" -- each is its own actor, so a viewer can hide the ones above
        to look underneath, which a single fused volume will not let you do.

        Each surface is interpolated by the shared
        :class:`~myflopy.modflow.mf6.grid.interpolated_surface.InterpolatedSurface`,
        so a VTK sheet and its Plotly counterpart are the same numbers.

        Returns the plotter and the sheet meshes, which the wrapping
        :class:`~myflopy.viz.VtkScene` carries as ``.meshes`` -- the handle for
        writing one out (``.save("contact.vtu")``) and opening it elsewhere.
        """

        pv = require("pyvista", feature="interactive 3-D scenes")
        import matplotlib as mpl

        names = self._resolve_surface_names(layer)
        xs = np.asarray(self.vor.centroids[0])
        ys = np.asarray(self.vor.centroids[1])
        colors = mpl.colormaps[cmap]

        # `off_screen`, matching `_vtk_plotter`: without it there is no render
        # window outside Jupyter, and `.save("sheets.png")` fails with "Nothing
        # to screenshot" in a plain script. A notebook forces it on anyway.
        plotter = pv.Plotter(off_screen=True, window_size=(width, height))
        # Masked to the grid's footprint like the Plotly sheets, so the two stay
        # the same numbers (ledger 314).
        footprint = self.vor.gdf_vorPolys.union_all()
        sheets = []
        for i, name in enumerate(names):
            interp = InterpolatedSurface(
                xs=xs, ys=ys, zs=np.asarray(self._surface_z(name)),
                surf_type="lyr", resolution=resolution,
                crs=str(getattr(self.vor, "crs", None)), clip=footprint,
            )
            gx, gy = interp.xy_meshgrid
            gz = np.asarray(interp.clipped_surface, dtype=float)
            grid = pv.StructuredGrid(
                np.asarray(gx, float), np.asarray(gy, float), gz * scale
            )
            # Cells with no data would otherwise render as a sheet pinned to z=0,
            # which reads as a real contact at sea level.
            grid["elevation"] = gz.ravel(order="F")
            sheet = grid.threshold(scalars="elevation")
            plotter.add_mesh(
                sheet,
                # PyVista wants floats or hex, not plotly's `rgb(r,g,b)` string,
                # so `_rgb` (which the Plotly path uses) is deliberately skipped.
                color=tuple(colors(i % colors.N)[:3]),
                opacity=opacity,
                show_edges=show_edges,
                label=str(name),
                # `name=` as well as `label=`: the label is what the LEGEND
                # says, the name is what `plotter.actors` is KEYED by. Without
                # it the keys are addresses -- `UnstructuredGrid(Addr=0x...)` --
                # so "hide the sheets above", which is the whole reason these
                # are separate actors, has no way to say which sheet it means.
                name=str(name),
                smooth_shading=True,
            )
            sheets.append(sheet)
        if len(names) > 1:
            plotter.add_legend(bcolor=None)
        plotter.add_axes()
        return plotter, sheets

    def _vtk_plotter(
        self, layers=None, *, color_by="layer", scale=8, cmap=None,
        width=900, height=580,
    ):
        """Build the PyVista plotter for the layered grid volume.

        Private: this is Layer 1 only -- it BUILDS the scene. Displaying it and
        writing it out belong to the :class:`~myflopy.viz.VtkScene` that wraps it.
        Reached as ``stack.plot.grid(backend="vtk")``.

        ``layers`` selects which layers to show: ``None``/``"all"`` (default) for
        every layer, a single layer name or index, or a list mixing them
        (e.g. ``["sand", "clay"]`` or ``[0, 2]``). Colors stay keyed to each
        layer's position, so a subset keeps the same colors it has in the full
        stack.

        Every cell carries ``layer``, ``thickness``, ``top``, ``botm`` and
        ``cellid`` whichever one ``color_by`` draws, so a picked cell can report
        all of them and a written ``.vtu`` keeps them all. ``cmap=None`` takes
        the default that suits the chosen scalar (:data:`_VTK_GRID_CMAPS`).

        Returns the plotter and, as the one mesh, the ASSEMBLED volume -- not
        the per-layer pieces the actors draw. Those are views of it, and it is
        the whole selection with every array on it, which is what makes
        ``scene.meshes[0].save("stack.vtu")`` the useful thing to hand to
        ParaView."""
        import tempfile

        import flopy
        from flopy.export.vtk import Vtk

        pv = require("pyvista", feature="interactive 3-D scenes")

        if color_by not in _VTK_GRID_CMAPS:
            raise ValueError(
                f"color_by must be one of {list(_VTK_GRID_CMAPS)}, not "
                f"{color_by!r}."
            )
        sel = self._resolve_layer_indices(layers)
        p = self.vor.get_disv_gridprops()
        # FloPy's VTK exporter wants a MODEL, so one is built for the geometry and
        # thrown away. Nothing is written into its workspace, and the mesh is in
        # memory by the end of the block -- which is what lets the directory go:
        # a `mkdtemp` here left one behind on every 3-D stack view (ledger 335).
        with tempfile.TemporaryDirectory(prefix="layer_view_") as ws:
            sim = flopy.mf6.MFSimulation(sim_name="layerview", sim_ws=ws)
            flopy.mf6.ModflowTdis(sim)
            flopy.mf6.ModflowIms(sim)
            gwf = flopy.mf6.ModflowGwf(sim, modelname="layers")
            flopy.mf6.ModflowGwfdisv(
                gwf, nlay=self.nlay, ncpl=p["ncpl"], nvert=len(p["vertices"]),
                vertices=p["vertices"], cell2d=p["cell2d"],
                top=self.top, botm=self.botm, idomain=self.idomain,
            )
            vtk = Vtk(model=gwf, vertical_exageration=1, binary=True, smooth=False)
            vtk.add_model(gwf)
            # Always attach a per-cell layer index so a subset can be selected.
            # Integer, so `add_array` leaves it alone -- see the float note below.
            vtk.add_array(np.repeat(np.arange(self.nlay)[:, None], p["ncpl"], axis=1), "layer")
            mesh = vtk.to_pyvista()
        if isinstance(mesh, pv.MultiBlock):
            mesh = mesh.combine()

        # Per-cell fields, attached HERE rather than through `Vtk.add_array`,
        # which NaN-masks FLOAT arrays wherever idomain == 0 -- i.e. on exactly
        # the cells a pinched-out layer creates, which are the ones worth
        # inspecting. (Measured: 138 of 843 cells on a stack with one pinching
        # layer; integer arrays come through intact.) `top` and `botm` are
        # already here from `add_model` and are OVERWRITTEN, not added: flopy's
        # own `top` is NaN for every layer below 0, so a per-layer top has to be
        # built from the contacts.
        mesh.cell_data["thickness"] = self.thickness.ravel()
        mesh.cell_data["top"] = np.vstack([self.top[None, :], self.botm[:-1]]).ravel()
        mesh.cell_data["botm"] = self.botm.ravel()
        mesh.cell_data["cellid"] = np.tile(np.arange(p["ncpl"]), self.nlay)

        if len(sel) != self.nlay:  # keep only the requested layers' cells
            layer_cell = np.asarray(mesh.cell_data["layer"]).astype(int)
            mesh = mesh.extract_cells(np.isin(layer_cell, sel))
        if color_by == "elevation":
            # The one POINT array, and the odd one out: vertex z, so it ramps
            # within a cell where `top`/`botm` are flat per-cell contacts.
            mesh["elevation"] = mesh.points[:, 2]

        mesh_kwargs = dict(
            show_edges=True,
            cmap=_VTK_GRID_CMAPS[color_by] if cmap is None else cmap,
        )
        if color_by == "layer":
            # Discrete colors keyed to each layer's global index; label only the
            # layers actually shown.
            mesh_kwargs.update(
                n_colors=self.nlay,
                clim=[-0.5, self.nlay - 0.5],
                annotations={float(i): self.names[i] for i in sel},
                scalar_bar_args=dict(title="layer", n_labels=0),
            )
        else:
            # An explicit clim over the whole selection, because each layer is
            # its own actor below and PyVista would otherwise scale each one to
            # its OWN range -- making a thin layer and a thick one look alike.
            finite = np.asarray(mesh[color_by], dtype=float)
            finite = finite[np.isfinite(finite)]
            mesh_kwargs.update(
                clim=[float(finite.min()), float(finite.max())] if finite.size else None,
                scalar_bar_args=dict(title=color_by),
            )

        plotter = pv.Plotter(off_screen=True, window_size=(width - 40, height - 40))
        # One actor per layer, named for it, rather than one fused volume. The
        # fused mesh drew the same picture, but a viewer could not take it
        # apart: `scene.scene.actors["clay"].visibility = False` needs a `clay`
        # actor to exist. Costs nlay draw calls instead of one.
        #
        # Each actor asks for the colour bar and the scene still gets exactly
        # one: PyVista keys bars by TITLE and reuses the existing one. That is
        # only correct because every actor shares the clim and cmap set above --
        # give them separate ranges and the single bar would describe one layer
        # while colouring all of them.
        layer_cell = np.asarray(mesh.cell_data["layer"]).astype(int)
        for k in sel:
            sheet = mesh.extract_cells(layer_cell == k)
            if sheet.n_cells == 0:  # selected but wholly absent from the mesh
                continue
            plotter.add_mesh(
                sheet, scalars=color_by, name=self.names[k], **mesh_kwargs,
            )
        plotter.set_scale(zscale=scale)
        plotter.add_axes()
        plotter.camera_position = "yz"
        return plotter, [mesh]

    def _thickness_values(self, layer=None):
        """Per-cell thickness (total, or a single named layer) and its label."""

        if layer is None:
            return self.thickness.sum(axis=0), "total thickness"
        return self.thickness[self.names.index(layer)], f"{layer!r} thickness"

    def _draw_thickness_map(self, *, layer=None, ax=None):
        """Render the thickness choropleth into ``ax`` (geopandas, no basemap)."""

        values, title = self._thickness_values(layer)
        gdf = self.vor.gdf_vorPolys.copy().assign(_thickness=values)
        if ax is None:
            _, ax = mpl_axes()
        gdf.plot(column="_thickness", ax=ax, legend=True)
        ax.set_title(f"{title} [{self.length_units}]")
        ax.set_aspect("equal")
        return ax

    @property
    def plot(self) -> StackPlots:
        """The plotting verbs for this layer geometry: ``map``, ``section``,
        ``surface``, ``grid``.

        Same verbs as ``model.plot`` and ``vor.plot``, and everything returned is
        a :class:`~myflopy.viz.Picture`. A stack has no results and no time, so
        there is no ``animate``; ``qc()`` stays a report rather than becoming a
        picture, because it is text you read.

        ``surface`` is a height field -- one contact as ``z(x, y)``. The layered
        cell VOLUME is ``grid(backend="vtk")``, a different shape entirely.

        Replaces ``thickness_map()``/``preview()`` (now ``plot.map()``),
        ``cross_section()`` (``plot.section()``), ``surface_3d()``
        (``plot.surface()``), ``vtk_3d()`` (``plot.grid()``) and ``views()``
        (compose what you want with ``myflopy.plot.mosaic``).

        Examples
        --------
        >>> stack.plot.map()                          # total thickness
        >>> stack.plot.map("sand", basemap=True)      # one layer, on a basemap
        >>> stack.plot.section(y=300)
        >>> stack.plot.surface("all")
        >>> stack.plot.grid(["sand", "clay"], scale=12)
        """

        return StackPlots(self)


class LayerSection(MplPicture):
    """A filled, layer-coloured cross-section through a built stack.

    Matplotlib by necessity, not by preference: the renderer is FloPy's
    ``PlotCrossSection``, and the Plotly ``GridSection`` draws cell outlines
    without layer fills. Rather than exempt it from the picture grammar, it is an
    :class:`~myflopy.viz.MplPicture` -- same verbs, Axes underneath.
    """

    def __init__(self, result, line=None, *, x=None, y=None, **kwargs):
        """Bind a section of ``result`` along ``line`` (or ``x=``/``y=``).

        The line and ``color_by`` are checked HERE, not at draw time. A Picture
        is lazy by design, but an argument that can never work is a caller
        mistake, and reporting it from ``.axes`` puts the traceback in the wrong
        place entirely."""

        self._result = result
        self._line = None if line is None else _section_line_points(line)
        self._x, self._y = x, y
        color_by = kwargs.get("color_by", "layer")
        if color_by not in _SECTION_COLOR_BY:
            raise ValueError(
                f"color_by must be one of {list(_SECTION_COLOR_BY)}, not "
                f"{color_by!r}."
            )
        self._kwargs = kwargs
        self.title = kwargs.get("title") or "Layer cross-section"

    def draw(self, ax=None, **kwargs):
        """Render the section into ``ax`` (or a new one) and return the Axes."""

        return self._result._draw_cross_section(
            self._line, x=self._x, y=self._y, ax=ax, **{**self._kwargs, **kwargs}
        )


class LayerThicknessMap(MplPicture):
    """Per-cell layer thickness, drawn on the grid without a basemap.

    Deliberately NOT a :class:`Choro`: a stack under construction is often on
    synthetic or local coordinates, and a web basemap would put it in the ocean.
    ``stack.plot.map(basemap=True)`` gives the georeferenced choropleth when the
    grid really is where it says it is.
    """

    def __init__(self, result, layer=None, *, extent=None, padding: float = 0.05):
        """Bind a thickness map of ``result`` (total, or one named layer).

        ``extent`` is the window it opens on (see :func:`myflopy.plot.map`),
        resolved here so a box off the grid raises at the call.
        """

        self._result = result
        self._layer = layer
        _, self.title = result._thickness_values(layer)
        vor = result.vor
        self._extent = resolve_extent(
            extent,
            crs=getattr(vor, "crs", None),
            bounds=vor.gdf_vorPolys.total_bounds,
        )
        self._padding = float(padding)

    def draw(self, ax=None, **kwargs):
        """Render the thickness map into ``ax`` (or a new one)."""

        ax = self._result._draw_thickness_map(layer=self._layer, ax=ax, **kwargs)
        if self._extent is not None:
            xmin, ymin, xmax, ymax = self._extent.padded(self._padding)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
        return ax


class LayerSurface(Picture):
    """One or more layer contacts as an interactive 3-D Plotly surface."""

    def __init__(self, result, layer="top", **kwargs):
        """Bind a 3-D surface of ``result`` for the chosen contact(s)."""

        self._result = result
        self._layer = layer
        self._kwargs = kwargs
        self._built = None

    @property
    def fig(self) -> Fig:
        """The assembled 3-D figure (built once, then cached)."""

        if self._built is None:
            self._built = self._result._surface_fig(self._layer, **self._kwargs)
        return self._built


#: `_vtk_plotter`'s own defaults, so `grid(backend="plotly")` can tell a value
#: the caller chose from one it merely inherited. Mirrored, not imported, because
#: the signature is the public contract and a test pins the two equal.
_VTK_GRID_DEFAULTS = {
    "color_by": "layer", "scale": 8, "cmap": None, "width": 900, "height": 580,
}


class StackPlots:
    """The plotting verbs for layer geometry -- ``stack.plot.map()`` (plan 8.5a).

    Three verbs, because a layer stack can answer three questions: how thick is
    it (:meth:`map`), what does it look like in section (:meth:`section`), and
    what shape is a given contact (:meth:`surface`). No ``animate`` -- a stack has
    no stress periods. ``qc()`` stays a method on the stack, because a QC report
    is text you read, not a picture you look at.

    Everything returned is a :class:`~myflopy.viz.Picture`, so it renders inline
    and answers ``.show()`` / ``.save(path)`` / ``.html(path)``. Two of the three
    are Matplotlib underneath and answer those over ``.axes``; only
    :meth:`surface` has a Plotly ``.fig``.
    """

    def __init__(self, result):
        """Bind the plotting verbs to a built :class:`LayerBuildResult`."""

        self.result = result

    def __repr__(self):
        """Name the verbs, since tab-completion is how this gets found."""

        return f"StackPlots({self.result.nlay} layers: map, section, surface, grid)"

    def map(self, layer=None, *, basemap: bool = False, **kwargs):
        """Per-cell thickness -- total, or one named layer.

        Draws on the grid with no basemap by default, because a stack is often
        still on synthetic coordinates. ``basemap=True`` routes through the
        shared choropleth instead, for a grid that really is georeferenced.

        .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.map`.

        The single map verb. Contours, observation markers, a hillshade and particle
        pathlines are all **options** here rather than verbs of their own, because a
        map is a plan view whatever is drawn on it.

        Every parameter below is spelled out in the signature rather than swept into
        ``**kwargs``, so an editor can complete and type-check them. The defaults are
        not restated here by hand -- each mirrors the default of whichever link in the
        ``map -> _choropleth_factory -> Choro`` chain owns that argument, and
        ``test_plot_vocabulary`` fails if the two ever drift apart.

        .. end generated

        Parameters
        ----------
        layer : str or int, optional
            One layer by name or index. With none, total thickness.
        basemap : bool, default False
            Route through the shared choropleth (``vor.plot.map``) so the cells
            sit on a web basemap. Requires a real CRS.
        **kwargs
            Choropleth styling, forwarded to :meth:`~myflopy.modflow.mf6.grid
            .plotting.GridPlots.map`. **Only meaningful with ``basemap=True``**
            -- the default renderer is Matplotlib and takes none of them.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        values : sequence of float, optional
            One value per cell -- heads, drawdown, K, a zone id, a residual, any
            per-cell array. Overrides whatever ``type`` would have read. Length must
            equal ``vor.ncpl``.
        per : int, optional
            Stress period to read (0-based). Mutually exclusive with ``kstpkper``;
            with neither, the model's first output time is used.
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)`` to read, as MODFLOW reports it. Use
            ``model.kstpkper`` to list what is available.
        per_timestep : {'last', 'first'} or int, default 'last'
            Which timestep WITHIN ``per`` to read, when a period has several. Ignored
            when ``kstpkper`` names the timestep outright.
        type : {'hds', 'conc', 'temp', 'rch', 'ks', 'trans', 'trans_max', 'custom'}, default 'hds'
            Which field to read when ``values`` is not given. ``'hds'`` heads,
            ``'conc'`` GWT concentration, ``'temp'`` GWE temperature -- these three
            also select the default sectioned hover -- ``'rch'`` the simulated
            recharge rate, and ``'ks'`` the horizontal conductivity ``npf.k``.
            ``'trans'`` is TRANSMISSIVITY, K times the SATURATED thickness at this
            map's ``per``/``kstpkper`` -- ``min(head, top) - bottom`` in a
            convertible layer, the full thickness in a confined one (MF6's own
            ICELLTYPE rule) -- and ``'trans_max'`` is K times the FULL thickness;
            both hover T with the K and thickness it multiplies, in the model's
            units. ``'ks'`` and ``'trans_max'`` are INPUTS, so they draw on a model
            that has not been run; the others are results and raise there. Case does
            not matter; any other value raises, naming this set, rather than drawing
            node numbers. The same two fields are nouns too:
            ``model.npf.transmissivity`` / ``.transmissivity_max``, whose
            ``.mosaic()`` draws one panel per layer.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        bgs : bool, default False
            Map depth below ground: each cell's head minus the model top, negative
            below it. Blank where the layer has no head, dry or inactive (ledger 315).
            ``show_mounding`` takes precedence.
        zoom : int, default 13
            Initial map zoom. Ignored when ``fit_bounds`` is True.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_layers : {'active', 'active+strip', 'all', 'none'}, optional
            How the per-layer profile renders in the hover.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        hover_fields : sequence of str, optional
            Extra columns to append to the hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        show_mounding : bool, default False
            Map mounding instead of head: each cell's head above its layer's bottom
            (above the model top with ``layer=-1``), clipped at 0, and add it to the
            hover. A DRY cell reads 0 and an inactive one is blank (ledger 315).
        hover_heads : bool, default True
            Include heads in the legacy flat hover.
        hover_ks : bool, default False
            Include hydraulic conductivity in the legacy flat hover.
        custom_hover : dict, optional
            Legacy flat hover: ``{label: per-cell sequence}``. Supplying it
            suppresses the default sectioned hover.
        rch_scale : float, optional
            Multiplier applied to recharge values when ``type='rch'``.
        animation_kstpkpers : sequence of tuple, optional
            The output times ``.ani`` steps through. Defaults to every time the model
            wrote.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.
        **trace_kwargs
            Anything else rides through to the ``go.Choroplethmap`` trace --
            ``zmid``, ``colorbar``, ``reversescale``, ``showscale``. These are
            genuinely open-ended and Plotly owns their names, so they are validated
            LATE, at render time, not here.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default), a
            :class:`~myflopy.viz.Picture`: it renders itself in Jupyter, and answers
            ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. It also
            carries ``.plot_mpl()`` for a static rendering and ``.ani`` for the
            animation over periods.

            With ``backend='mpl'``, a bare Matplotlib ``Figure`` -- not a Picture, so
            use ``.savefig(path)`` and ``.axes[0]`` rather than the picture verbs.

            .. end generated

        Raises
        ------
        TypeError
            If styling arguments are given without ``basemap=True``. They used
            to be accepted and silently discarded, which quietly produced an
            unstyled picture.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        See Also
        --------
        section : the same data as a vertical slice.
        grid : the mesh with no values and no basemap (and no CRS needed).
        myflopy.plot.animate : flip a sequence of these through time.

        Examples
        --------
        >>> model.plot.map(layer=0)                          # this model's heads
        >>> model.plot.map(values=drawdown, layer=0)         # any per-cell array
        >>> model.plot.map(layer=0, contours=True, contour_levels=8)
        >>> model.plot.map(layer=0, features="wells.gpkg", hillshade_path="hs.tif")
        >>> model.plot.map(values=k, logscale=True, colorscale="Viridis")  # K over decades
        >>> model.plot.map(layer=0, features={"Wells": "wells.gpkg",
        ...                                   "Ponds": plot.Feature("ponds.gpkg", fill=True)})
        >>> model.plot.map(layer=0, zmin=100, zmax=125).save("heads.png")
        >>> vor.plot.map(values=node_ids)                    # a bare grid
        """

        if basemap:
            values, _ = self.result._thickness_values(layer)
            # `GridPlots`, not `myflopy.plot`: that module is a layer ABOVE this
            # one in the import graph, so reaching it would need a deferred
            # import and the exact-match ratchet only moves down. Same function
            # either way -- `vor.plot.map` is what the front door calls too.
            return GridPlots(self.result.vor).map(values=list(values), **kwargs)
        # The one option the static thickness map honours -- framing, which a
        # Matplotlib picture cannot do by zooming (ledger 309).
        extent = kwargs.pop("extent", None)
        padding = kwargs.pop("bounds_padding", 0.05) if extent is not None else 0.05
        if kwargs:
            raise TypeError(
                f"{', '.join(sorted(kwargs))} style the choropleth, which is "
                f"only drawn with basemap=True; the default thickness map is "
                f"Matplotlib and honours only extent=/bounds_padding=."
            )
        return LayerThicknessMap(self.result, layer=layer, extent=extent, padding=padding)

    def section(
        self,
        line=None,
        *,
        x=None,
        y=None,
        color_by: str = "layer",
        cmap: str = "tab10",
        show_grid: bool = True,
        legend: bool = True,
        title: str | None = None,
        **kwargs,
    ) -> LayerSection:
        """A filled, layer-coloured cross-section: ``stack.plot.section(y=300)``.

        .. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.section`.

        ``section`` draws one of two pictures, and ``fill`` chooses which:

        ===============================  =============================================
        ``fill="layer"`` (the default)   The CROSS-SECTION: every cell the line cuts,
                                         drawn from its bottom to its top and coloured
                                         by layer, with the simulated water level on
                                         top.
        ``fill="results"``               The same cells, coloured by the model's own
                                         result -- head on GWF, concentration on GWT,
                                         temperature on GWE -- with a colorbar.
        ``fill=<(nlay, ncpl) array>``    The same cells, coloured by any field of your
                                         own: K, a zone id, a residual.
        ``fill=None``                    The LINE PROFILE: head against distance along
                                         the line, one line per ``layer=``, plus the
                                         model top. No cells.
        ===============================  =============================================

        Both pictures draw on both backends, and each backend draws the same cells.
        ``backend="mpl"`` -- the default -- returns a static Matplotlib ``Figure``;
        ``backend="plotly"`` returns an interactive picture, which for the line
        profile can also step through the saved times (``.ani``).

        Through a bare GRID (``vor.plot.section(line)``) there are no results: the
        cross-section is the layer geometry alone, ``fill="results"`` raises, and
        ``fill=None`` draws the cell outlines.

        **Which arguments go with which picture.** An argument given to the other
        picture RAISES and says which picture it belongs to, rather than being
        quietly ignored:

        - *either picture:* ``line`` or ``cells``, ``per`` or ``kstpkper``,
          ``reverse``, ``section_name``, ``backend``;
        - *the cross-section:* ``layers``, ``head_layers``, ``layer_labels``,
          ``show_grid``, ``legend``, ``fill_cmap``, ``fill_label``;
        - *the line profile (``fill=None``):* ``layer``, ``interpolate``,
          ``interpolator``, ``use_rbf``, ``num_points``, ``clip``, ``x_or_y``,
          ``show_model_top``, ``show_model_btm``, ``surf_type``,
          ``animation_kstpkpers``.

        Mind ``layers=`` against ``layer=``: the cross-section's chooses whose CELLS
        are drawn, the profile's whose head LINES.

        .. end generated

        Parameters
        ----------
        line : LineString, MultiLineString, path, or sequence of points, optional
            The section line: a shapely geometry, a path to a ``.shp``/``.gpkg``
            to read it from, or the points themselves as ``[(x, y), ...]``.
            Alternatively give ``x=`` or ``y=`` for an axis-aligned slice. With
            none of them, a West-East line through the grid centre.
        x, y : float, optional
            Draw the section along a constant x or constant y.
        color_by : {'layer', 'thickness'}, default 'layer'
            Cell scalar the fill is keyed to: discrete layer colours, or a
            continuous viridis thickness field.
        cmap : str, default 'tab10'
            Colormap for that scalar.
        show_grid : bool, default True
            Draw cell edges over the fill.
        legend : bool, default True
            Include the layer legend.
        title : str, optional
            Plot title. Defaults to "Layer cross-section".
        **kwargs
            Forwarded to :class:`LayerSection`.

            .. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
        cells : int or list of int, optional
            The section path as cell numbers, in order; the line runs through their
            centroids.
        per : int, optional
            Stress period, zero-based: the LAST step saved in it is drawn, so
            ``per=0`` is the end of the first period as the output saved it, on
            both pictures (and on ``map`` and ``surface``). With neither ``per`` nor
            ``kstpkper``, both pictures draw the model's FIRST saved time. Giving
            both raises. A model only.
        kstpkper : tuple of (int, int), optional
            The exact ``(time step, period)`` to draw; one the output did not save
            raises. A model only.
        fill : {'layer', 'results'}, array-like, or None, default 'layer'
            Which picture, and what colours the cells -- see the table above.
            ``'layer'``: by layer, water level on top. ``'results'``: by the
            model's own result, with a colorbar; a model only. An ``(nlay, ncpl)``
            array: by that field; a flat ``(ncpl,)`` array raises rather than
            painting every layer the same. ``None``: the line profile, no cells.
        layers : int or sequence of int, optional
            *(cross-section)* Which layers' CELLS to draw, zero-based; every layer
            by default. The rest are masked out and the vertical extent is cropped
            to what remains. A layer outside the model raises.
        head_layers : int, sequence of int, or None, default 0
            *(cross-section, ``fill='layer'``)* Whose water levels go on top:
            layer 0 by default, a list draws one labelled line each, ``None`` none.
            Ignored by ``'results'`` and an array, which paint the field onto the
            cells instead. A model only.
        layer_labels : sequence of str, optional
            *(cross-section)* Names for the layers, one per layer in model order,
            for the legend and the water levels. Defaults to the names the model was
            built with (``ModelContext(surfaces=stack.build(vor))`` keeps them),
            else ``Layer 1..N``.
        fill_cmap : str, default 'viridis'
            *(cross-section, ``fill='results'`` or an array)* A Matplotlib colormap
            name. Both backends colour from it.
        fill_label : str, optional
            *(cross-section, ``fill='results'`` or an array)* The colorbar title.
        layer : int or list of int, default 0
            *(line profile)* Whose head lines to draw; a list overlays several. An
            interpolated profile draws one layer, so a list of more raises there.
        interpolate : bool, default False
            *(line profile)* Sample an interpolated surface along the line, rather
            than stepping from cell centre to cell centre.
        interpolator : str, optional
            *(line profile)* The interpolation method by name -- ``'griddata'``,
            ``'rbf'``, ``'linearND'`` or ``'cloughTocher2D'``. Giving one turns
            ``interpolate`` on.
        use_rbf : bool, default True
            *(line profile, interpolating)* With no ``interpolator``: radial basis
            functions (``True``) or ``griddata`` (``False``). ``False`` without
            interpolation raises.
        num_points : int, default 100
            *(line profile, interpolating)* How many points to sample along the line.
            Changed without interpolation, it raises: nothing is sampled then.
        clip : Path or geometry, optional
            *(line profile, interpolating)* Clip the interpolated surface to this
            region; raises without interpolation.
        x_or_y : {'x', 'y'}, default 'x'
            *(line profile)* With a single cell in ``cells``, run the line east-west
            (``'x'``) or north-south (``'y'``) through it.
        show_model_top : bool, default True
            *(line profile)* Draw the model top.
        show_model_btm : bool, default False
            *(line profile)* Draw the bottoms of the layers in ``layer``.
        surf_type : {'hds', 'lyr'}, default 'hds'
            *(line profile, interpolating)* Interpolate the head field, or with
            ``'lyr'`` the TOP of ``layer`` -- the model top for ``layer=0``, the bottom
            of layer ``k - 1`` for ``layer=k`` -- from the model's own
            discretization. ``'lyr'`` without interpolation raises.
        animation_kstpkpers : sequence of tuple, optional
            *(line profile, ``backend='plotly'``)* The times ``.ani`` steps through;
            every saved time by default. A Matplotlib figure has no ``.ani``, so
            this raises with ``backend='mpl'``.
        reverse : bool, default False
            Draw the line from its other end: a section drawn west to east comes back
            east to west, with distance measured from the new start.
        section_name : str, optional
            The cross-section's title, or the line profile's legend name.
        backend : {'mpl', 'plotly'}, default 'mpl'
            ``'mpl'`` returns a static :class:`matplotlib.figure.Figure`;
            ``'plotly'`` an interactive :class:`~myflopy.viz.Picture`. Accepts
            ``'matplotlib'``/``'static'`` and ``'interactive'`` as aliases; anything
            else raises.

            .. end generated

        Returns
        -------
        LayerSection
            A Matplotlib :class:`~myflopy.viz.Picture`; ``.axes`` rather than
            ``.fig``.

        Raises
        ------
        TypeError
            If ``line`` is not a geometry, a path, or a sequence of points.
        ValueError
            If ``line`` resolves to fewer than two points, or ``color_by`` is
            not one of the values above. Both are raised HERE, not later from
            ``.axes``.

            .. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
        See Also
        --------
        map : the same model in plan view.
        myflopy.plot.grid : the mesh itself.

        Examples
        --------
        >>> model.plot.section(line=line)                    # the cross-section, by layer
        >>> model.plot.section(line=line, fill="results")    # cells coloured by head
        >>> model.plot.section(line=line, layers=[0, 1], head_layers=[0, 1],
        ...                    layer_labels=["fill", "sand", "clay", "till"])
        >>> model.plot.section(line=line, legend="outside right").savefig("section.png")
        >>> model.plot.section(line=line, backend="plotly")  # the same, interactive
        >>> model.plot.section(line=line, fill=None)         # the head line profile
        >>> model.plot.section(cells=[1653, 651, 1241], fill=None, layer=[0, 1])
        >>> model.plot.section(line=line, fill=None, interpolate=True, num_points=200)
        >>> model.plot.section(cells=cells, fill=None, backend="plotly").ani
        >>> vor.plot.section(line)                           # a bare grid's layers
        """

        return LayerSection(
            self.result, line, x=x, y=y, color_by=color_by, cmap=cmap,
            show_grid=show_grid, legend=legend, title=title, **kwargs,
        )

    def surface(
        self,
        layer="top",
        *,
        backend: str = "plotly",
        resolution: int = 120,
        colorscale: str = "Earth_r",
        color_by: str | None = None,
        opacity: float | None = None,
        height: int | None = None,
        scale: float = 8,
        cmap: str = "tab10",
        width: int = 900,
        show_edges: bool = False,
        **kwargs,
    ):
        """One or more contacts as an interactive 3-D surface.

        ``layer`` is a name, a list of names, or ``"all"``; discover them with
        ``result.surface_names``. That is a height field ``z(x, y)``; the layered
        grid VOLUME is :meth:`grid`, a different shape entirely.

        .. generated from myflopy.plot.surface -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.surface`.

        Always Plotly, and always a height field. The 3-D grid VOLUME and 3-D
        particle pathlines are different shapes and live on :func:`grid` with
        ``backend="vtk"``: a ``backend=`` switch should change the renderer, not
        what is being drawn.

        .. end generated

        Parameters
        ----------
        layer : str or list of str, default 'top'
            Which contact(s) to draw. ``"all"`` draws every one.
        resolution : int, default 120
            Interpolation grid size per axis.
        colorscale : str, default 'Earth_r'
            Plotly colorscale for the height field.
        color_by : str, optional
            Colour by a scalar other than elevation.
        opacity : float, optional
            Surface opacity, useful when stacking several contacts.
        height : int, optional
            Figure height in pixels.
        backend : {'plotly', 'vtk'}, default 'plotly'
            ``'plotly'`` overlays the surfaces in one figure. ``'vtk'`` renders
            each contact as its own interactive sheet in a PyVista scene, which
            is the one to reach for when you need to look at surfaces
            INDIVIDUALLY -- every sheet is a separate actor, so you can hide the
            ones above and see underneath. Needs the ``viz3d`` extra.
        scale : float, default 8
            *(vtk only)* Vertical exaggeration.
        cmap : str, default 'tab10'
            *(vtk only)* Colormap the sheets are coloured from, one per surface.
        width : int, default 900
            *(vtk only)* Scene width in pixels.
        show_edges : bool, default False
            *(vtk only)* Draw the interpolation mesh on each sheet.
        **kwargs
            Forwarded to :class:`LayerSurface` (plotly) or the scene builder.

            .. generated from myflopy.plot.surface -- edit there, then run scripts/derive_docstrings.py
        per : int, optional
            Stress period (0-based), drawn at its LAST saved step, as on ``map`` and
            ``section`` (ledger 337). Mutually exclusive with ``kstpkper``: both
            raise.
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)``. Defaults to the model's FIRST output time.
        surf_type : {'hds', 'lyr'}, default 'hds'
            Interpolate the head field, or with ``'lyr'`` the TOP of ``layer`` -- the
            model top for ``layer=0``, the bottom of layer ``k - 1`` for ``layer=k``
            -- read from the model's own discretization (a bare grid's
            ``gdf_topbtm``). Anything else raises.
        use_rbf : bool, default False
            Radial-basis interpolation instead of linear.
        interpolator : str, optional
            Override the interpolation method by name.
        clip : Path, str or (Multi)Polygon, optional
            Restrict the surface to this region: only the cells it covers are
            interpolated, and the drawn surface is masked to it. Without one the
            surface is masked to the grid's own footprint, because the interpolator
            alone fills the convex hull of the cell centres -- which spans every
            concavity of a domain that is not convex. ``.clipped_fig(clip=False)`` on
            the result draws that unmasked interpolation.
        crs : str, optional
            Override the grid's CRS.
        xs, ys, zs : array-like, optional
            Supply the point cloud directly instead of reading it from a model.

            .. end generated

        Returns
        -------
        LayerSurface or VtkScene
            A :class:`~myflopy.viz.Picture` either way. The VTK scene exposes
            ``.scene`` (the PyVista ``Plotter``) rather than ``.fig``.

        Raises
        ------
        ValueError
            If ``backend`` is neither ``'plotly'`` nor ``'vtk'``.

            .. generated from myflopy.plot.surface -- edit there, then run scripts/derive_docstrings.py
        See Also
        --------
        grid : the mesh itself, including the 3-D volume via ``backend="vtk"``.

        Examples
        --------
        >>> model.plot.surface(layer=0)
        >>> model.plot.surface(layer=0, surf_type="lyr")     # layer elevations
        >>> model.plot.surface(layer=0, resolution=400).html("surface.html")
        """

        if backend == "vtk":
            plotter, sheets = self.result._vtk_surface_plotter(
                layer, resolution=resolution, scale=scale, cmap=cmap,
                opacity=1.0 if opacity is None else opacity,
                width=width, height=height or 580, show_edges=show_edges,
                **kwargs,
            )
            return VtkScene(plotter, meshes=sheets, title="layer surfaces")
        if backend != "plotly":
            raise ValueError(f"backend must be 'plotly' or 'vtk', not {backend!r}.")
        return LayerSurface(
            self.result, layer, resolution=resolution, colorscale=colorscale,
            color_by=color_by, opacity=opacity, height=height, **kwargs,
        )

    def grid(
        self,
        layers=None,
        *,
        backend: str = "vtk",
        color_by: str = "layer",
        scale: float = 8,
        cmap: str | None = None,
        width: int = 900,
        height: int = 580,
    ):
        """The layered grid mesh itself.

        ``backend="vtk"`` (the default here) renders the cell VOLUME in 3-D,
        coloured by layer -- the picture the old ``vtk_3d()`` drew, minus its
        habit of writing an HTML file into the working directory on every call.
        ``backend="plotly"`` gives the flat 2-D mesh instead, the same picture as
        ``vor.plot.grid()``.

        A ``backend`` switch is honest here because both branches draw the SAME
        subject -- this grid -- and differ only in renderer. That is why the 3-D
        volume is `grid`, not `surface`: `surface` means a height field.

        .. generated from myflopy.plot.grid -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.grid`.

        The one picture :func:`map` cannot give you. A choropleth colours cells
        against a web basemap and so **requires a CRS**; this draws in the grid's own
        coordinates and needs none, which makes it the view for a grid you are still
        refining, before there is a model or a projection.

        Both backends draw the same subject -- this grid -- so ``backend=`` switches
        only the renderer. That is why the 3-D volume is ``grid`` and not
        ``surface``: ``surface`` means a height field.

        .. end generated

        Parameters
        ----------
        layers : str or int or list, optional
            *(vtk only)* Which layers to show: a name, an index, or a list mixing
            them. Colours stay keyed to each layer's position, so a subset looks
            the same as it does in the full stack.
        backend : {'vtk', 'plotly'}, default 'vtk'
            ``'vtk'`` renders the 3-D volume and needs the ``viz3d`` extra;
            ``'plotly'`` draws the flat 2-D mesh.
        color_by : {'layer', 'thickness', 'top', 'botm', 'cellid', 'elevation'}, default 'layer'
            *(vtk only)* Cell scalar to colour by. ``'layer'`` is discrete, with
            each layer's name in the colour bar; the rest are continuous fields.
            ``'top'`` and ``'botm'`` are the flat per-cell contacts, where
            ``'elevation'`` is vertex z and so ramps within a cell.

            Every one of these is attached to the mesh regardless of which is
            drawn, so ``scene.meshes[0]`` carries them all into a ``.vtu``.
        scale : float, default 8
            *(vtk only)* Vertical exaggeration.
        cmap : str, optional
            *(vtk only)* Colormap for ``color_by``. Defaults to the one that
            suits the chosen scalar -- ``tab10`` for layers, ``viridis`` for
            thickness, ``terrain`` for elevations.
        width, height : int, default 900, 580
            *(vtk only)* Scene size in pixels.

            .. generated from myflopy.plot.grid -- edit there, then run scripts/derive_docstrings.py
        pathlines : DataFrame, optional
            Particle track records, drawn as time-coloured tubes over the 3-D mesh.
            Requires ``backend="vtk"``; in plan view the equivalent is
            ``map(pathlines=...)``. Passing it with ``backend="plotly"`` raises.

        vertical_exaggeration : float, default 1.0
            *(vtk only)* Multiplier on z, to make a thin model legible.
        model_style : {'wireframe', 'surface', 'points'}, default 'wireframe'
            *(vtk only)* How the grid itself is drawn beneath the tracks.
        model_opacity : float, default 0.25
            *(vtk only)* Opacity of the grid, so tracks inside it stay visible.
        pathline_cmap : str, default 'viridis'
            *(vtk only)* Colormap for the time-coloured tubes.
        pathline_width : float, default 4.0
            *(vtk only)* Tube width.
        show_edges : bool, default True
            *(vtk only)* Draw cell edges on the mesh.
        off_screen : bool, default True
            *(vtk only)* Render without opening a window -- the right default in a
            notebook or on a headless machine.
        **kwargs
            Forwarded to the backend's builder.

            .. end generated

        Returns
        -------
        VtkScene or GridMesh
            A :class:`~myflopy.viz.Picture` either way. ``VtkScene`` exposes
            ``.scene`` instead of ``.fig``.

        Raises
        ------
        ValueError
            If ``backend`` is neither value, if ``color_by`` is not one of the
            values above, or if a vtk-only argument is given with
            ``backend="plotly"``.

            .. generated from myflopy.plot.grid -- edit there, then run scripts/derive_docstrings.py
        See Also
        --------
        map : values over the cells, on a basemap.
        surface : a 3-D height field, which is a different shape.
        myflopy.layers.StackPlots.grid : the LAYER-stack 3-D volume, which takes
            ``layers``/``scale``/``color_by``/``cmap`` -- a different builder, and
            not reachable through this function.

        Examples
        --------
        >>> vor.plot.grid()                         # the 2-D mesh
        >>> vor.plot()                              # the same thing, shorthand
        >>> vor.plot.grid().plot_mpl()              # FloPy's matplotlib renderer
        >>> stack.plot.grid(["sand", "clay"])       # 3-D volume, a subset of layers
        >>> model.plot.grid(pathlines=run.track_records, backend="vtk")
        """

        scene_args = {"color_by": color_by, "scale": scale, "cmap": cmap,
                      "width": width, "height": height}
        if backend == "plotly":
            stray = [n for n, v in scene_args.items() if v != _VTK_GRID_DEFAULTS[n]]
            if layers is not None:
                stray.insert(0, "layers")
            if stray:
                raise ValueError(
                    f"{', '.join(stray)} configure the 3-D scene; the flat mesh "
                    f"has no use for them. Drop backend='plotly', or drop these."
                )
            return GridPlots(self.result.vor).grid()
        if backend != "vtk":
            raise ValueError(
                f"backend must be 'vtk' or 'plotly', not {backend!r}."
            )
        plotter, meshes = self.result._vtk_plotter(layers, **scene_args)
        return VtkScene(plotter, meshes=meshes, title="layered grid")


class LayerStack:
    """Author a model's vertical layering from a top surface + named layers.

    ``LayerStack`` is the **friendly facade** for building MODFLOW layer geometry.
    You declare the model top, then ``.add(...)`` one named layer at a time (by its
    *bottom* surface or its *thickness*), and ``.build()`` resolves the stack into
    DISV-ready ``top`` / ``botm`` / ``idomain`` arrays. It is a thin layer over the
    :class:`~myflopy.surfaces.LayerSurfaces` engine -- it compiles to it (see
    :meth:`_layer_surfaces`) and never re-implements the sampling / reconcile /
    pinch-out logic. Each surface is an atomic :class:`~myflopy.surfaces.Surface`
    (use the aliases :class:`Raster`, :class:`Contours`, :class:`Points`,
    :class:`Flat`, :class:`Array`).

    What :meth:`build` does for you: samples every surface onto the grid cells
    (area-weighted by default), resolves them **top-down** so a flat/relative
    surface stays flat where it fits and is lowered only where it would intrude on
    the surface above ("flat where possible, fit between"), and turns sub-minimum
    or inverted thicknesses into **pinch-outs** (idomain) per your per-layer policy.

    Parameters
    ----------
    vor : VoronoiGridPlus or None, optional
        The grid whose cells the surfaces are sampled onto, used by :meth:`build`,
        :meth:`qc`, :meth:`to_disv` and :attr:`plot` whenever they are not given
        one. Default None: the DEFERRED form -- declare the layering now and pass
        the grid to those methods later, or bind it with :meth:`for_grid`
        (:attr:`plot` then draws on a coarse :meth:`draft_grid`). A
        :class:`~myflopy.surfaces.Surface` here raises ``TypeError``: it is the
        first positional argument, so ``LayerStack(ground)`` would bind the top as
        the grid.
    top : Surface, str, pathlib.Path, int or float
        The model-top surface. Required (the ``None`` default raises
        ``TypeError``). A :class:`~myflopy.surfaces.Surface` is used as given; a
        ``str``/``Path`` is read as a :class:`Raster`; an ``int``/``float`` is a
        :class:`Flat` elevation. Anything else -- a bare ``numpy`` array included
        -- raises ``TypeError``; wrap per-cell values in :class:`Array`.
    length_units : str, optional
        The model length unit. Default ``"feet"``. At :meth:`build` every surface
        declaring its own ``units=`` is converted into it (known spellings:
        ``feet``/``ft``/``foot``/``us-ft``/``usft``, ``meters``/``meter``/
        ``metres``/``m``, ``centimeters``/``cm``; any other raises ``ValueError``
        there). :meth:`to_disv` writes it as DISV ``LENGTH_UNITS``: ``"feet"``,
        ``"meters"`` and ``"centimeters"`` become ``FEET``/``METERS``/
        ``CENTIMETERS``, and any other spelling is written upper-cased as given.
        Also labels the plot axes.
    time_units : str, optional
        The model time unit. Default ``"days"``. Carried onto the built
        :class:`LayerBuildResult` (its report prints it) and saved by
        :meth:`to_dict`; nothing in the stack converts with it.

    Examples
    --------
    >>> from myflopy import LayerStack
    >>> from myflopy.layers import Raster, Contours
    >>> stack = (
    ...     LayerStack(vor, top=Raster("ground.tif"))   # land surface from a DEM
    ...     .add("alluvium",   thickness=25)            # 25-ft upper aquifer
    ...     .add("aquitard",   thickness=8, pinch="inactive")   # pinches out where thin
    ...     .add("bedrock",    bottom=Contours("bedrock_top.shp"))
    ... )
    >>> layers = stack.build()          # -> layers.top, layers.botm, layers.idomain, layers.nlay
    >>> print(stack.qc())               # NaN/thickness/connectivity report
    >>> layers.plot.section(x=500)      # quick section to eyeball it

    Feed the result straight into ``mf.disv`` and the model context::

        gp = vor.get_disv_gridprops()
        ctx = mf.ModelContext(grid=vor, domain=layers.idomain)
        flow = mf.gwf("flow", context=ctx, packages=[
            mf.disv(nlay=layers.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                    vertices=gp["vertices"], cell2d=gp["cell2d"],
                    top=layers.top, botm=layers.botm, idomain=layers.idomain),
            ...,
        ])

    See Also
    --------
    from_modflow : Seed an editable stack from an existing model's top/botm.
    build : Resolve the stack into DISV arrays (the per-layer pinch options).
    qc : Geometry quality-control report (NaN, thin/pinched, connectivity).
    """

    def __init__(
        self,
        vor=None,
        top=None,
        *,
        length_units: str = "feet",
        time_units: str = "days",
    ):
        """Start an editable layer stack with the model-top surface ``top``.

        Layers are added below the top with :meth:`add`; ``length_units`` /
        ``time_units`` are carried onto the built :class:`LayerBuildResult`.

        **``vor`` is optional.** Surfaces are lazy and ``add`` resolves nothing, so
        the whole layering -- names, sources, per-layer ``min_thickness`` and
        ``pinch`` -- can be declared before a grid exists::

            stack = (mf.LayerStack(top=ground)          # no grid yet
                     .add("sand", thickness=40.0, pinch="inactive")
                     .add("clay", bottom=clay_base))
            layers = stack.build(vor)                   # grid arrives here

        Supply it later to :meth:`build`, :meth:`qc` or :meth:`to_disv`, or bind
        it once with :meth:`for_grid`. Passing it here still works and is right
        whenever the grid already exists -- the deferred form exists so that
        declaring the layering does not force the grid to be built first.

        Raises
        ------
        TypeError
            If ``top`` is missing, or if a :class:`~myflopy.surfaces.Surface` is
            passed as ``vor`` -- ``LayerStack(ground)`` reads as "the deferred
            form" but binds ``ground`` to the grid, which would fail much later
            and somewhere else.
        """

        if top is None:
            hint = (
                " It looks like you passed the top surface positionally: the "
                "first positional argument is the GRID, so write "
                "`LayerStack(top=...)` for the deferred form, or "
                "`LayerStack(vor, top)`."
            ) if isinstance(vor, Surface) else ""
            raise TypeError(f"LayerStack needs a `top` surface.{hint}")
        if isinstance(vor, Surface):
            raise TypeError(
                "LayerStack's first argument is the grid, not a surface. Write "
                "`LayerStack(top=...)` to defer the grid, or `LayerStack(vor, top)`."
            )

        self.vor = vor
        self._top = _coerce_surface(top)
        self._layers: list[_Layer] = []
        #: Declared excavations (:meth:`excavate`), cut in after reconcile and pinch.
        self._excavations: list[_Excavation] = []
        #: Cached draft grid for `.plot` on a gridless stack; see `draft_grid`.
        self._draft = None
        self.length_units = length_units
        self.time_units = time_units

    #: Token :func:`myflopy.specs._json_value` writes this object under. See
    #: :attr:`myflopy.surfaces.Surface._SPEC_TOKEN` for why it is not a ``kind`` tag.
    _SPEC_TOKEN = "$layer_stack"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready payload: the DECLARATION, not the built arrays.

        What rides: the top surface, each unit's bottom/thickness surface and its
        ``min_thickness``/``pinch``/``split``/``names``, the units, and -- only when
        there are any -- the excavations (ledger 323). What does NOT: ``vor`` (a stack is
        declarable without one, and a grid cannot be inlined -- ledger 193) and
        ``_draft`` (a cached preview grid).

        **An excavation entry** is ``{"name", "bottom", "zone" | "area", "layer"?,
        "min_remnant", "source_shas"?}``: ``bottom`` is the floor's Surface payload, or
        ``{"column": <name>}`` for a per-feature floor read from the footprint file;
        ``zone`` is the footprint FILE as a posix path (resolved to an absolute one at
        declaration) and ``area`` an inline ``mf.Area`` document (WKT + CRS). The file
        rides under ``"zone"`` -- the key ``Surface.where`` already uses -- so the
        bundler, the build-time anchor, ``runs.drift()``'s ``built_from`` and the hash
        stamping all find it with no wiring of their own. Nothing resolved on a grid is
        written (no cells, no sampled bed): the stack re-reads the file and re-covers
        the cells on whatever mesh it is built on.

        **A stack does not determine its own arrays**, which is why this is a recipe
        and not a result: ``reconcile``, ``min_sep``, ``trigger_sep`` and ``method``
        live on the ``build()``/``to_disv()`` call and move the answer -- same stack,
        same grid, a different ``reconcile`` measured a different ``botm``. Whoever
        rebuilds from this payload has to supply them again (ledger 269).

        Returns
        -------
        dict
            Plain Python types only (``str``, ``float``, ``int``, ``bool``, ``None``,
            lists, dicts) -- tuples inside a surface payload are tagged by the spec
            serializer, not here.

        Raises
        ------
        ValueError
            A surface holding per-cell arrays (``Surface.from_array``/``from_points``)
            cannot be written down (ledger 269).

        Examples
        --------
        >>> import myflopy as mf
        >>> import shapely
        >>> pond = mf.Area(shapely.box(0.0, 0.0, 50.0, 50.0), crs="EPSG:2927")
        >>> stack = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        ...          .excavate(bottom=93.0, within=pond, name="pond"))
        >>> entry = stack.to_dict()["excavations"][0]
        >>> sorted(entry), entry["bottom"], entry["min_remnant"]
        (['area', 'bottom', 'min_remnant', 'name'], {'kind': 'flat', 'value': 93.0}, 0.1)
        >>> mf.LayerStack.from_dict(stack.to_dict()).to_dict() == stack.to_dict()
        True
        """

        payload: dict[str, Any] = {
            "top": self._top.to_dict(),
            "layers": [
                {
                    "name": layer.name,
                    "surface": layer.surface.to_dict(),
                    **{
                        key: getattr(layer, key)
                        for key in ("min_thickness", "pinch", "split", "names")
                        if getattr(layer, key) is not None
                    },
                }
                for layer in self._layers
            ],
            "length_units": self.length_units,
            "time_units": self.time_units,
        }
        # Only when there are any: an un-excavated stack's payload -- and so its digest
        # -- is byte-identical to what it was before excavations existed.
        if self._excavations:
            payload["excavations"] = [
                _excavation_payload(excavation) for excavation in self._excavations
            ]
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """:meth:`to_dict`, with each file fingerprinted by CONTENT where it can be.

        ``to_dict`` nests each surface as a plain dict, so ``specs._canonical`` never
        meets the ``Surface`` objects themselves and their own ``digest_payload`` would
        never be consulted. This is what reaches them (ledger 293's second pass) -- and
        an excavation's footprint file, whose hash replaces its path through the same
        helper (``surfaces._content_keyed``). An inline ``mf.Area``'s ``source`` and
        ``note`` are dropped: provenance and a comment, not content.

        Returns
        -------
        dict
            The payload a package's digest fingerprints: moving a stamped file keeps
            it, editing the file moves it.

        Raises
        ------
        ValueError
            As :meth:`to_dict`.

        Examples
        --------
        >>> import myflopy as mf
        >>> import shapely
        >>> pond = mf.Area(shapely.box(0.0, 0.0, 50.0, 50.0), note="drawn by hand")
        >>> stack = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        ...          .excavate(bottom=93.0, within=pond, name="pond"))
        >>> "note" in stack.to_dict()["excavations"][0]["area"]
        True
        >>> "note" in stack.digest_payload()["excavations"][0]["area"]
        False
        """

        payload = self.to_dict()
        payload["top"] = self._top.digest_payload()
        for entry, layer in zip(payload["layers"], self._layers, strict=True):
            entry["surface"] = layer.surface.digest_payload()
        for entry, excavation in zip(
            payload.get("excavations", ()), self._excavations, strict=True
        ):
            if excavation.bottom is not None:
                entry["bottom"] = excavation.bottom.digest_payload()
            if "area" in entry:
                entry["area"] = {
                    key: value for key, value in entry["area"].items()
                    if key not in ("source", "note")
                }
            _content_keyed(entry, excavation.source_shas)
        return payload

    @classmethod
    def from_dict(cls, data: dict[str, Any], *, vor=None) -> LayerStack:
        """Rebuild a stack from :meth:`to_dict`'s payload, optionally bound to ``vor``.

        The excavations are rebuilt DIRECTLY from their entries, never through
        :meth:`excavate`: that verb reads the footprint file, and a bundled payload
        names a project-RELATIVE path that only the build anchors to the project root
        -- checking it here would check it against the working directory (ledger 293's
        two-bases trap). A missing file is found at build, where it is read. Every
        OTHER refusal of :meth:`excavate` runs here, through the same code: a
        hand-written payload (``SimulationSpec.from_yaml``) reaches this path, and
        without them a repeated name dropped one excavation's record while its cut
        stayed in the arrays (ledger 323).

        Parameters
        ----------
        data : dict
            A payload :meth:`to_dict` wrote, after a JSON or YAML round trip if you
            like (tuples tagged by the spec serializer must already be untagged, as
            ``specs._spec_value`` does). Keys it does not know -- at the top level or
            in an excavation entry -- are REFUSED, so a payload written by a later
            version is never read as a stack missing whatever it carried.
        vor : VoronoiGridPlus, optional
            The grid to bind the rebuilt stack to; omitted, it is gridless (the
            payload never carries one).

        Returns
        -------
        LayerStack
            A new stack; its surfaces, excavations and units are new objects.

        Raises
        ------
        ValueError
            An unknown key; an excavation entry naming both or neither of ``zone`` and
            ``area``; an unknown ``Surface`` field in any surface payload. And each
            refusal :meth:`excavate` makes without reading a file: an excavation name
            that is empty, repeated, or names a unit or layer; a ``bottom`` measured
            from the surface above (``offset_below`` and the like) or algebra over a
            ``fill="propagate"`` source; ``layer`` or a column floor
            (``{"column": ...}``) with an inline ``area`` rather than a ``zone`` file;
            ``min_remnant`` negative or not finite.
        TypeError
            A ``min_remnant`` that is not a number, or a ``layer`` that is not a
            ``str`` (an integer reads as a MODEL layer, as in :meth:`excavate`).
        KeyError
            A payload missing ``top``, or an excavation entry missing ``name``,
            ``bottom`` or ``min_remnant``.

        Examples
        --------
        >>> import myflopy as mf
        >>> stack = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        >>> mf.LayerStack.from_dict(stack.to_dict()).names
        ['sand']
        >>> mf.LayerStack.from_dict({**stack.to_dict(), "ponds": []})
        Traceback (most recent call last):
        ...
        ValueError: A LayerStack payload carries keys this version does not know: ['ponds']. Known: ['excavations', 'layers', 'length_units', 'time_units', 'top'].
        >>> import shapely
        >>> pond = {"name": "pond", "bottom": {"kind": "flat", "value": 93.0},
        ...         "area": mf.Area(shapely.box(0, 0, 10, 10)).to_dict(), "min_remnant": 0.1}
        >>> mf.LayerStack.from_dict({**stack.to_dict(), "excavations": [pond, pond]})
        Traceback (most recent call last):
        ...
        ValueError: an excavation named 'pond' already exists; pass another name=.
        """

        unknown = sorted(set(data) - _STACK_KEYS)
        if unknown:
            raise ValueError(
                f"A LayerStack payload carries keys this version does not know: {unknown}. "
                f"Known: {sorted(_STACK_KEYS)}."
            )
        stack = cls(
            vor,
            top=Surface.from_dict(data["top"]),
            length_units=data.get("length_units"),
            time_units=data.get("time_units"),
        )
        for entry in data.get("layers", ()):
            stack._layers.append(
                _Layer(
                    name=entry["name"],
                    surface=Surface.from_dict(entry["surface"]),
                    min_thickness=entry.get("min_thickness"),
                    pinch=entry.get("pinch"),
                    split=entry.get("split"),
                    names=entry.get("names"),
                )
            )
        reserved = stack._reserved_names()
        for entry in data.get("excavations", ()):
            stack._excavations.append(
                _excavation_from_payload(
                    entry,
                    taken={excavation.name for excavation in stack._excavations},
                    reserved=reserved,
                )
            )
        return stack

    def _require_grid(self, vor, caller: str):
        """Resolve the grid for ``caller``, preferring an explicit one.

        One message for every consumer, because "NoneType has no attribute
        ncpl" three frames down is the failure a deferred stack would otherwise
        produce.
        """

        resolved = self.vor if vor is None else vor
        if resolved is None:
            raise ValueError(
                f"LayerStack.{caller}() needs a grid. This stack was declared "
                f"without one, so pass it now -- `stack.{caller}(vor)` -- or bind "
                f"it once with `stack.for_grid(vor)`."
            )
        return resolved

    def _georeferenced_surfaces(self) -> list[Surface]:
        """Every surface in the stack that carries its own extent, composites included.

        A raster or contour surface knows where it is; ``flat``/``array`` and the
        algebra kinds do not, but they are *built from* ones that might, so the
        walk descends through ``operands``.
        """

        found, seen = [], set()

        def walk(surface):
            if id(surface) in seen:
                return
            seen.add(id(surface))
            if surface.kind in ("raster", "contours"):
                found.append(surface)
            for operand in surface.operands or ():
                if isinstance(operand, Surface):
                    walk(operand)

        walk(self._top)
        for layer in self._layers:
            if isinstance(layer.surface, Surface):
                walk(layer.surface)
        return found

    def draft_grid(self, *, cells: int = 400, extent=None, crs: str | None = None):
        """Build a throwaway Voronoi grid covering this stack's own extent.

        For LOOKING at a stack before its real grid exists. The extent comes from
        the first raster- or contour-backed surface in the stack -- those are the
        only kinds that know where they are -- and the mesh is deliberately coarse.

        **Not a modelling grid.** It has no boundary polygon, no refinement, and
        an extent that is whatever your DEM happens to cover rather than your
        model domain. :meth:`build`, :meth:`qc` and :meth:`to_disv` will not use
        it and still demand a real grid; only :attr:`plot` falls back to it.

        Parameters
        ----------
        cells : int, default 400
            Roughly how many cells to aim for. Coarse on purpose -- a few hundred
            renders in well under a second and is plenty to see shape, ordering
            and where a layer pinches.
        extent : tuple of float, optional
            ``(xmin, ymin, xmax, ymax)`` to use instead of the surfaces' own --
            for previewing over your real domain, or a corner of it.
        crs : str, optional
            Overrides the CRS read from the source raster.

        Returns
        -------
        VoronoiGridPlus

        Raises
        ------
        ValueError
            If no surface carries an extent and none was given. A stack of
            ``Flat``/``Array`` surfaces describes thicknesses with no notion of
            where they are, so there is nothing to infer.
        """

        if extent is None:
            rasterio = require("rasterio", feature="reading a surface extent for draft_grid()")
            for surface in self._georeferenced_surfaces():
                source = surface.resolve_source()
                if isinstance(source, Path):
                    with rasterio.open(source) as handle:
                        extent = tuple(handle.bounds)
                        crs = crs or str(handle.crs)
                    break
            else:
                raise ValueError(
                    "draft_grid() needs an extent: no surface in this stack carries "
                    "one. `Raster` and `Contours` know where they are; `Flat`, "
                    "`Array` and the algebra built on them do not. Pass "
                    "extent=(xmin, ymin, xmax, ymax), or use a real grid."
                )

        xmin, ymin, xmax, ymax = (float(v) for v in extent)
        width, height = xmax - xmin, ymax - ymin
        if width <= 0 or height <= 0:
            raise ValueError(
                f"draft_grid() got an empty extent: {(xmin, ymin, xmax, ymax)}."
            )

        # Triangle writes its input and output files here; `build()` reads the mesh
        # back into memory and the Voronoi grid is made from that, so the files are
        # spent once the grid exists. A `mkdtemp` left them behind on every draft
        # grid (ledger 335).
        with tempfile.TemporaryDirectory(prefix="draft_grid_") as workspace:
            tri = TriangleGrid(model_ws=workspace, angle=30)
            tri.set_domain_rectangle(
                x_dist=width, y_dist=height, origin=(xmin, ymin),
                max_area=(width * height) / max(int(cells), 1),
            )
            tri.build()
            return VoronoiGridPlus(tri, crs=crs)

    def for_grid(self, vor) -> LayerStack:
        """Return a copy of this stack bound to ``vor``.

        The layering is the expensive thing to write and the grid is the thing
        you change, so binding returns a NEW stack rather than mutating this one:
        one declaration can serve a coarse test grid and a fine production grid
        without being restated.

            coarse = stack.for_grid(vor_coarse).build()
            fine   = stack.for_grid(vor_fine).build()

        Also what makes ``.plot`` reachable on a deferred stack, since a property
        cannot take a grid argument.

        Parameters
        ----------
        vor : VoronoiGridPlus or None
            The grid the copy is bound to -- what its :meth:`build`, :meth:`qc`,
            :meth:`to_disv` and :attr:`plot` use when not given one. It REPLACES any
            grid this stack holds; ``None`` gives an unbound (deferred) copy. The
            copy shares this stack's top, unit declarations, excavations and units;
            this stack is not changed. A :class:`~myflopy.surfaces.Surface` raises
            ``TypeError``, as in the constructor.
        """

        bound = LayerStack(
            vor,
            self._top,
            length_units=self.length_units,
            time_units=self.time_units,
        )
        bound._layers = list(self._layers)
        # The excavations ride along: a bound copy that dropped them would build the
        # stack with no pond, silently -- and `mf.disv` copies an array-holding stack
        # through here.
        bound._excavations = list(self._excavations)
        return bound

    @classmethod
    def from_modflow(
        cls,
        vor,
        source,
        *,
        names: list[str] | None = None,
        resample: bool = True,
        length_units: str = "feet",
        time_units: str = "days",
    ) -> LayerStack:
        """Build a stack on ``vor`` from an existing MODFLOW model's top/botm.

        ``source`` is a flopy model or modelgrid; its top and per-layer bottoms
        become this stack's surfaces (interpolated onto ``vor`` when
        ``resample=True``, used verbatim when ``False``). Name the layers with
        ``names`` (defaults to ``layer1..N``). Edit the returned stack like any
        other -- e.g. ``.replace(...)`` a bottom, then ``.build()`` / ``.to_disv()``.

        Parameters
        ----------
        vor : VoronoiGridPlus or None
            The grid the new stack is bound to (the constructor's ``vor``); ``None``
            gives a deferred stack whose grid is supplied at :meth:`build`. With
            ``resample=False`` it must be the SAME grid, cell for cell, as
            ``source``.
        source : flopy model or flopy modelgrid
            Where the surfaces come from: anything with a ``.modelgrid`` (a FloPy
            model), or a FloPy modelgrid itself. Its ``top`` becomes the stack's top
            and each layer of its ``botm`` one unit's bottom, read through
            :func:`modflow_surfaces` -- so a structured or vertex (DISV) grid, one
            ``top`` value per column.
        names : list of str, optional
            One unit name per source layer, top to bottom. Default None:
            ``"layer1"``, ``"layer2"``, ... A list of the wrong length raises
            ``ValueError``.
        resample : bool, optional
            Default True: each surface is interpolated from the source's CELL
            CENTRES onto ``vor``'s centroids (``scipy.interpolate.griddata``,
            linear), so it transfers onto a DIFFERENT grid; a target cell outside
            the convex hull of the source centres samples NaN. False: the arrays are
            used verbatim as per-cell values, so ``vor`` must have exactly the
            source's cells. Either way the surfaces hold per-cell data, which
            :meth:`to_dict` refuses -- a stack from here cannot be saved as a recipe.
        length_units : str, optional
            The stack's length unit; see the class ``length_units``. Default
            ``"feet"``. The source's elevations are taken as already in it -- no
            conversion is applied to them.
        time_units : str, optional
            The stack's time unit; see the class ``time_units``. Default ``"days"``.

        Raises
        ------
        ValueError
            ``names`` does not have one entry per source layer.
        """
        top_s, botm_s = modflow_surfaces(source, resample=resample)
        if names is None:
            names = [f"layer{i + 1}" for i in range(len(botm_s))]
        if len(names) != len(botm_s):
            raise ValueError(f"{len(names)} names given for {len(botm_s)} layers.")
        stack = cls(vor, top=top_s, length_units=length_units, time_units=time_units)
        for name, surface in zip(names, botm_s, strict=False):
            stack.add(name, bottom=surface)
        return stack

    # -- authoring -------------------------------------------------------- #
    def add(
        self,
        name: str,
        *,
        bottom=None,
        thickness=None,
        min_thickness: float | None = None,
        pinch: str | None = None,
        fill: str | None = None,
        split=None,
        names=None,
    ) -> LayerStack:
        """Append a named geologic unit beneath the current bottom. Returns ``self``.

        Define the unit either by its **bottom** surface or its **thickness**
        (exactly one). Thickness is measured down from the surface above, so layers
        compose naturally as you stack them.

        ``split=`` discretizes the unit into several MODEL layers without
        changing its geometry -- see below.

        Parameters
        ----------
        name : str
            Unit name (used in QC, plots, and as the surface label). ``"top"`` is
            reserved, a name already used by a unit raises, and so does one naming an
            excavation of this stack (:meth:`excavate`) -- as does any layer name the
            unit builds, a generated ``<name>_i`` of ``split`` or one of ``names``.
        bottom : Surface, str, pathlib.Path or float, optional
            The unit's bottom as a :class:`Surface` / alias (e.g.
            ``Contours("base.shp")``, ``Array(values)``, ``Flat(90.0)``), a raster
            path (read as ``Raster``) or a number (``Flat``). Mutually exclusive with
            ``thickness``.
        thickness : float or Surface, optional
            Constant thickness below the surface above. Mutually exclusive with
            ``bottom``. For a per-cell thickness MAP write ``bottom=Isopach(map)``
            -- a bare elevation Surface here is refused, because it would set the
            bottom rather than cut a thickness (the relative kinds ``Isopach``,
            ``constant_thickness`` and ``offset_below`` are accepted).
        min_thickness : float, optional
            Minimum thickness for the UNIT; thinner cells are handled per
            ``pinch``. Defaults to the stack-wide value passed to :meth:`build`.
            Not divided among a split unit's layers -- see ``split``.
        pinch : str, optional
            What to do where the unit is thinner than ``min_thickness``:
            ``"floor"`` (hold the cell open at the minimum; :meth:`build`'s
            ``default_pinch``, used when this is omitted), ``"passthrough"``
            (idomain -1: removed, flow passes vertically through), or
            ``"inactive"`` (idomain 0 -- a true pinch-out). This entry used to say
            ``"passthrough"`` was the default; :meth:`build`'s is ``"floor"``.
        fill : str, optional
            For raster/derived bottoms, how to fill cells with no source data
            (``"propagate"`` inherits the surface above -> the unit pinches there).
        split : int or sequence of float, optional
            Discretize this unit into several model layers. ``split=3`` gives
            three equal shares; ``split=[0.3, 0.7]`` gives those shares, which
            must sum to 1. The layers are named ``<name>_1 .. <name>_N``, and
            ``<name>`` remains addressable through
            :attr:`LayerBuildResult.units`.

            The unit's GEOMETRY is unchanged: the last cut lands exactly on the
            declared bottom. ``min_thickness`` and ``pinch`` stay unit-scoped, so
            a unit either pinches out whole or not at all -- evaluating them per
            slice produces holes inside a unit that physically exists.

            A unit whose bottom is measured from the layer above (an
            ``Isopach``) cannot be split; give it an absolute bottom, or declare
            it with ``thickness=`` and split that.
        names : sequence of str, optional
            Names for the split layers, top to bottom, replacing the generated
            ``<name>_1 .. <name>_N``. Requires ``split`` of 2 or more, and must
            give one name per layer -- naming only some of them would leave
            which layer a name refers to depending on where you stopped
            counting. The UNIT keeps ``name``, so
            :attr:`LayerBuildResult.units` and
            :meth:`LayerBuildResult.per_layer` are unaffected by the choice.

        Returns
        -------
        LayerStack
            ``self``.

        Raises
        ------
        ValueError
            ``name`` reserved or taken by a unit; ``name`` or a layer name it builds
            (``<name>_i``, ``names``) taken by an excavation; both or neither of
            ``bottom``/``thickness``; a bad ``split`` or ``names``.
        TypeError
            An elevation Surface given as ``thickness``; a ``bottom`` of an
            unsupported type; ``names`` that is not a sequence of names.

        Examples
        --------
        >>> import myflopy as mf
        >>> stack = mf.LayerStack(top=mf.Flat(100.0))
        >>> stack.add("sand", thickness=20).names                   # 20-ft unit
        ['sand']
        >>> _ = stack.add("clay", thickness=5, pinch="inactive")    # pinches out where thin
        >>> _ = stack.add("till", thickness=60, split=[0.25, 0.75]) # -> 15 ft, then 45 ft
        >>> _ = stack.add("rock", bottom=mf.Flat(0.0), split=3,
        ...               names=["upper rock", "mid rock", "lower rock"])  # unit stays "rock"
        >>> stack.names
        ['sand', 'clay', 'till', 'rock']
        """
        _check_layer_name(name, self._layers)
        surface = _make_surface(bottom, thickness, fill)
        shares = _split_shares(split, name)  # reject a bad split here, not at build
        _split_names(names, shares, name)    # ...and a bad names= list with it
        self._refuse_excavation_name(name, split, names)
        self._layers.append(_Layer(name, surface, min_thickness, pinch, split, names))
        return self

    def insert_below(
        self, name: str, new_name: str, *, bottom=None, thickness=None,
        min_thickness: float | None = None, pinch: str | None = None,
        fill: str | None = None, split=None, names=None,
    ) -> LayerStack:
        """Insert a new unit directly below the existing unit ``name``. Returns ``self``.

        Parameters
        ----------
        name : str
            The existing unit to insert below; an unknown name raises ``KeyError``.
        new_name : str
            The new unit's name -- the rules of :meth:`add`'s ``name``.
        bottom : Surface, str, pathlib.Path or float, optional
            As on :meth:`add`.
        thickness : float or Surface, optional
            As on :meth:`add`.
        min_thickness : float, optional
            As on :meth:`add`.
        pinch : str, optional
            As on :meth:`add` (``"floor"``, ``"passthrough"`` or ``"inactive"``;
            omitted, :meth:`build`'s ``default_pinch``).
        fill : str, optional
            As on :meth:`add`.
        split : int or sequence of float, optional
            As on :meth:`add`.
        names : sequence of str, optional
            As on :meth:`add`.

        Returns
        -------
        LayerStack
            ``self``.

        Raises
        ------
        KeyError
            ``name`` names no unit.
        ValueError
            Any refusal of :meth:`add` (``new_name`` taken by a unit, it or a layer
            name it builds taken by an excavation, a bad ``split``/``names``, both or
            neither of ``bottom``/``thickness``).
        TypeError
            As on :meth:`add`.

        Examples
        --------
        >>> import myflopy as mf
        >>> stack = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20).add("rock", thickness=50)
        >>> stack.insert_below("sand", "clay", thickness=5).names
        ['sand', 'clay', 'rock']
        """
        _check_layer_name(new_name, self._layers)
        surface = _make_surface(bottom, thickness, fill)
        _split_names(names, _split_shares(split, new_name), new_name)
        self._refuse_excavation_name(new_name, split, names)
        self._layers.insert(
            self._index(name) + 1,
            _Layer(new_name, surface, min_thickness, pinch, split, names),
        )
        return self

    def replace(
        self, name: str, *, bottom=None, thickness=None,
        min_thickness=_UNSET, pinch=_UNSET, fill=None, split=_UNSET, names=_UNSET,
    ) -> LayerStack:
        """Update an existing unit in place; unspecified fields are kept.

        This is also how a unit is split after the fact --
        ``stack.replace("sand", split=3)``. There is deliberately no ``.split()``
        verb: this method already edits a declared unit by name, already returns
        ``self``, and already keeps what you do not mention.

        ``split`` and ``names`` are kept independently but validated TOGETHER,
        so changing one to a length the other cannot match raises here rather
        than at build. Dropping a split from a named unit therefore reads
        ``replace(name, split=None, names=None)`` -- the alternative, silently
        discarding names the caller wrote, is the kind of quiet renaming this
        module exists to prevent.

        Parameters
        ----------
        name : str
            The unit to change; an unknown name raises ``KeyError``. The unit keeps
            it (a rename is ``remove`` + :meth:`insert_below`).
        bottom : Surface, str, pathlib.Path or float, optional
            A new bottom, as on :meth:`add`. Mutually exclusive with ``thickness``;
            omitted with ``thickness`` too, the unit keeps its current bottom.
        thickness : float or Surface, optional
            A new constant thickness, as on :meth:`add`. Mutually exclusive with
            ``bottom``.
        min_thickness : float or None, optional
            As on :meth:`add`; ``None`` falls back to :meth:`build`'s
            ``default_min_thickness``. Omitted, the unit keeps its current value.
        pinch : str or None, optional
            As on :meth:`add` (``"floor"``, ``"passthrough"`` or ``"inactive"``);
            ``None`` falls back to :meth:`build`'s ``default_pinch``. Omitted, kept.
        fill : str, optional
            As on :meth:`add`. With a new ``bottom``/``thickness`` it applies to
            that; given alone it changes the fill of the unit's CURRENT bottom.
            Omitted (``None``), the current bottom's fill is kept.
        split : int, sequence of float or None, optional
            As on :meth:`add`; ``None`` makes the unit one layer again. Omitted,
            kept. Validated with ``names`` as a pair.
        names : sequence of str or None, optional
            As on :meth:`add`; ``None`` restores the generated
            ``<name>_1 .. <name>_N``. Omitted, kept. Validated with ``split`` as a
            pair.

        Returns
        -------
        LayerStack
            ``self``.

        Raises
        ------
        KeyError
            ``name`` names no unit.
        ValueError
            Both ``bottom`` and ``thickness``; a ``split``/``names`` pair that
            does not match (checked only when either is given); a layer the new
            ``split``/``names`` would build -- a generated ``<name>_i`` or one of
            ``names`` -- already naming an excavation of this stack
            (:meth:`excavate`).
        TypeError
            As on :meth:`add`.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)
        >>> stack = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20).add("clay", thickness=5)
        >>> stack.replace("sand", split=3).build(vor).names
        ['sand_1', 'sand_2', 'sand_3', 'clay']
        >>> stack.replace("sand", thickness=10.0, split=None).build(vor).botm[:, 0].tolist()
        [90.0, 85.0]
        """
        idx = self._index(name)
        layer = self._layers[idx]
        if bottom is not None or thickness is not None:
            surface = _make_surface(bottom, thickness, fill)
        else:
            surface = layer.surface if fill is None else _dc_replace(layer.surface, fill=fill)
        resolved_split = layer.split if split is _UNSET else split
        resolved_names = layer.names if names is _UNSET else names
        if split is not _UNSET or names is not _UNSET:
            _split_names(
                resolved_names, _split_shares(resolved_split, name), name
            )
            self._refuse_excavation_name(name, resolved_split, resolved_names)
        self._layers[idx] = _Layer(
            name,
            surface,
            layer.min_thickness if min_thickness is _UNSET else min_thickness,
            layer.pinch if pinch is _UNSET else pinch,
            resolved_split,
            resolved_names,
        )
        return self

    def remove(self, name: str) -> LayerStack:
        """Remove a layer by name.

        Edits this stack IN PLACE and returns it, so calls chain. Excavations are
        not touched.

        Parameters
        ----------
        name : str
            The UNIT name as declared by :meth:`add` (one of :attr:`names`) -- a
            split unit is removed whole, by its unit name; a slice name such as
            ``"sand_2"`` is not accepted.

        Raises
        ------
        KeyError
            No unit is named ``name``; the message lists the names there are.
        """
        del self._layers[self._index(name)]
        return self

    def _refuse_excavation_name(self, name, split=None, names=None) -> None:
        """A unit or layer may not take an excavation's name (they share reports and keys).

        Checks every name the unit BUILDS -- the generated ``<name>_1..N`` of a split
        as well as ``names=`` -- because those are what ``LayerBuildResult.names``
        holds beside ``excavations``.
        """

        dug = {excavation.name for excavation in self._excavations}
        clash = sorted(dug & _unit_names(name, split, names))
        if clash:
            raise ValueError(
                f"{clash} already names an excavation of this stack; a unit and its "
                "layers need names of their own."
            )

    def excavate(self, bottom, within, *, name=None, layer=None, min_remnant=0.1) -> LayerStack:
        """Declare a void cut into the stack -- a pond, basin, pit, trench, or a lake. Returns ``self``.

        Inside ``within`` every layer surface above ``bottom`` is lowered to it: the
        first active cell under the bed then has its top EXACTLY at the bed, the
        cells above it go idomain 0 (never -1 -- the space is lake water, not a
        pass-through layer), and the ground outside is untouched. For a natural lake
        whose DEM top is its water surface, and for a dug basin alike. Nothing is
        sampled until :meth:`build`, and ``build`` cuts AFTER reconcile and after the
        pinch policy -- cutting before them was measured wrong under every pinch
        policy and reconcile mode (a dug layer came back active at its minimum
        thickness under ``"floor"``, the connected top sat ``min_sep`` below the bed
        under ``"inactive"``, a split unit was squeezed rather than removed).

        It mutates and returns ``self``, like :meth:`add`, :meth:`insert_below`,
        :meth:`replace` and :meth:`remove`. So build the un-excavated baseline FIRST
        if you want one: after ``pond = stack.excavate(...)``, ``pond is stack``, and
        a baseline built from ``stack`` later is excavated too. A package already
        declared from the stack holds its own copy and is not affected.

        **The footprint is the lake's.** Cells are chosen by the rule
        ``mf.lak(lakes=)`` uses -- a polygon covers a cell by AREA, a sliver of a
        cell is noise (ledger 241) -- through the same reader and the same function,
        so pass the SAME polygons you pass to ``mf.lak(lakes=)``. (A centroid rule,
        ``Surface.where``'s, missed 50 of a real pond's 612 lake cells.) The bed is
        sampled the lake's way too -- area-weighted, whatever ``build(method=)``
        says -- so a lake declared on the excavated stack reads the same float at
        every cell.

        **The sliver rule: absorb, not deepen.** Where the cell left directly under
        the bed is thinner than ``min_remnant`` and an active cell lies below it, it
        is ABSORBED: its bottom rises to the bed (zero thickness, idomain 0) and the
        cell below becomes the connected one, its top still exactly at the bed -- the
        column's thickness is conserved and the bed is never moved. An inactive or
        pass-through cell directly under the bed is absorbed the same way whatever
        ``min_remnant`` says, because a lake's vertical connection cannot reach
        through it. A thin cell with nothing active below it is kept (and reported).

        **With a transparent bed the connected cell sets the infiltration rate**
        (MF6's vertical lake conductance is ``area / (1/bedleak + 0.5 * thk / K33)``):
        a larger ``min_remnant`` is then a GEOLOGY change, not a solver setting --
        on a real pond ``1.0`` moved 172 of the lake's vertical connections from a
        1-ft K33 = 11.7 unit onto a K33 = 0.41 one and lowered the mound under the
        pond 17%. Read :meth:`LayerBuildResult.excavation_report` before choosing it.

        **Several excavations** are cut together from the same pre-cut arrays, so
        declaration order does not matter; two footprints sharing a cell are
        REFUSED at build (a berm narrower than a cell does that -- measured, ponds 2
        or 10 ft apart on 17-ft cells shared 2 cells).

        **What "model top" means afterwards.** The built ``top`` is the excavated
        top and ``LayerBuildResult.ground`` keeps the ground. With ``attach=True``
        the grid's ``gdf_topbtm`` carries the excavated top and ``grid.excavations``
        the records ``mf.lak`` reads. So ``mf.rch.array(irch="top_active")`` lands on
        the cell under the bed (which MF6 then zeroes while a lake is wet there --
        give the lake ``rainfall=``), a ``.gpkg`` boundary on a dug cell is dropped
        with the partial-miss warning, ``CellSurfaceOffset("model_top", -2)`` is the
        bed less 2, and a PRT release point in layer 0 over a dug cell lands on no
        active cell and is reported, not moved down. A package declared against a
        context built BEFORE the cut still names the removed cells, and the model's
        build REFUSES it, naming this excavation: a list boundary or an observation
        on a dug cell, an ``irch``/``ievt`` array pointing at one (ledger 321).
        ``irch=None`` -- one layer everywhere -- is warned instead, with the rate
        MF6 will drop (ledger 322).

        Parameters
        ----------
        bottom : Surface, float, str or pathlib.Path
            The floor, an ABSOLUTE elevation in the stack's ``length_units``.

            * A :class:`~myflopy.surfaces.Surface` or alias -- ``mf.Contours``,
              ``mf.Raster``, ``mf.Flat``, algebra on an absolute surface
              (``ground - 9``). ``units=`` on it is converted to the stack's
              ``length_units``, as a unit bottom is. A BARE ``fill="propagate"``
              source takes the ground where it has no data, so those cells are
              not cut (algebra over one is refused, below).
            * Any real number (``int``, ``float``, numpy scalars), a flat floor.
            * A ``str`` naming a COLUMN of ``within``'s attribute table -- one flat
              floor per feature (``"flat_bottom_elev"`` over several ponds) -- when
              ``within`` is a file that has it; otherwise a ``str`` or ``Path`` naming
              an existing raster file.

            Always sampled area-weighted over each cell (the lake's sampler,
            :func:`myflopy.surfaces.sample_bed`), with the built ground as the
            surface above. REFUSED: a surface measured from the one above it
            (``offset_below``, ``constant_thickness``, ``Isopach``, ``toward``, and
            anything built on them) -- the stack would measure it from the ground and
            ``mf.lak`` from the excavated top (+5.0 ft in 70 of 70 cells, measured);
            write ``ground - 9`` instead. REFUSED too: algebra over a
            ``fill="propagate"`` source (``Surface.raster(p, fill="propagate") - 2``),
            which reads the surface above wherever the source has no data -- the
            stack would cut those cells to ``ground - 2`` and ``mf.lak`` would read
            ``ground - 4`` there (14 of 27 cells 2.0 ft apart, measured). Fill the
            raster's gaps, or pass the propagate source bare.
        within : str, pathlib.Path, Area, geopandas.GeoDataFrame, geopandas.GeoSeries or shapely geometry
            The footprint -- the polygons you pass to ``mf.lak(lakes=)``.

            * A vector path (``str`` or ``Path``; a bare string is always a path),
              read at ``layer`` by the reader ``mf.lak`` uses. Resolved to an
              absolute path at the call and read then, so a missing file, a missing
              layer, a line or a zero-area feature refuses here; read again at
              build. Saved as that PATH (ledger 323): ``save(bundle_sources=True)``
              carries the file, a lazy ``mf.disv`` digests its content, and
              ``runs.drift()`` reports it edited.
            * An ``mf.Area`` -- kept inline, and saved inline as WKT + CRS;
              ``crs=None`` means the grid's frame.
            * A GeoDataFrame, GeoSeries or shapely geometry, coerced to an inline
              ``Area`` (its features unioned, so a multi-part footprint selects per
              PART; a bare geometry with no CRS is in the grid's frame).

            Polygons only: a line or point would dig the cells along it by length.
            REFUSED: cell numbers and boolean masks -- they name cells by position
            on one mesh and would dig different cells after ``with_grid``.
        name : str, optional
            The excavation's key in ``LayerBuildResult.excavations`` and its label
            in reports and refusals (converted with ``str``). Defaults to the file
            stem, then ``Area.name``, then ``"excavation_<n>"``. Must be non-empty,
            unique among the stack's excavations, and must not name a unit or layer -- nor
            may a unit added, inserted or re-split later build a layer of that name
            (a generated ``<unit>_i`` included; :meth:`add`, :meth:`insert_below` and
            :meth:`replace` refuse it).
        layer : str, optional
            The vector layer of a ``within`` FILE, by NAME (a numpy string is
            converted with ``str``); omitted, the first. Refused with any other
            ``within`` -- it could not act -- and an integer is refused
            (``TypeError``): in a layer stack it reads as a MODEL layer, and an
            excavation cuts every layer above its bed. Pass ``mf.lak`` the same
            ``layer=``.
        min_remnant : float, default 0.1
            The thinnest cell allowed directly under the bed, in length units,
            converted with ``float``. A thinner one is a sliver and is absorbed
            (above). ``0`` turns sliver absorption off; inactive cells under the bed
            are still absorbed. Deliberately not the build's ``min_sep`` nor a
            unit's ``min_thickness``: both would change the geology the stack
            declares (see the note on transparent beds).

        Returns
        -------
        LayerStack
            ``self``.

        Raises
        ------
        ValueError
            At the call: a ``bottom`` measured from the surface above, or algebra
            over a ``fill="propagate"`` source; a ``str``
            bottom that is neither a column of ``within`` nor a raster file (a
            column needs a FILE ``within``); a NULL or non-numeric column value;
            ``layer`` without a file; a line, point or zero-area footprint; an
            empty or duplicate ``name``, or one naming a unit or layer; ``min_remnant``
            negative or not finite. At :meth:`build`: a footprint covering no cell
            (both extents and both CRSs, ledger 211); two excavations sharing a
            cell, or two features of one file giving one cell different column
            floors; a bed with no data at every footprint cell; a column left with
            no active cell under the bed (a bed below the stack's bottom included).
        TypeError
            ``bottom`` or ``within`` of an unsupported type (cell numbers and
            masks named as such), a non-numeric ``min_remnant``, or a ``layer``
            that is not a ``str``.
        FileNotFoundError
            A ``within`` or ``bottom`` path that does not exist.

        Warns
        -----
        UserWarning
            At build: some (not all) features cover no cell, named -- by the file's
            stem when it holds one feature, else by zero-based row (``'row 1'``), or
            ``'<name> part <i>'`` for a part of an ``Area``; footprint cells where
            the bed has no data, which are not cut.

        See Also
        --------
        LayerBuildResult.excavate : The same cut on a built result.
        LayerBuildResult.excavation_report : What each excavation removed and connected.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=50.0)
        >>> pond = mf.Area(vor.gdf_vorPolys.geometry.loc[[5, 6]].union_all())
        >>> stack = (mf.LayerStack(top=mf.Flat(100.0))
        ...          .add("fill", thickness=5.0)
        ...          .add("sand", thickness=30.0)
        ...          .add("clay", thickness=20.0))
        >>> baseline = stack.build(vor)                      # BEFORE excavating
        >>> dug = stack.excavate(bottom=93.0, within=pond, name="pond").build(vor, attach=True)
        >>> float(dug.top[5]), dug.idomain[:, 5].tolist(), float(dug.ground[5])
        (93.0, [0, 1, 1], 100.0)
        >>> dug.excavations["pond"].cells
        (5, 6)
        >>> vor.excavations["pond"] is dug.excavations["pond"]
        True
        """

        self._excavations.append(
            _declare_excavation(
                bottom, within, name=name, layer=layer, min_remnant=min_remnant,
                taken={excavation.name for excavation in self._excavations},
                reserved=self._reserved_names(),
                index=len(self._excavations),
            )
        )
        return self

    def _reserved_names(self) -> set[str]:
        """Every unit name and every model-layer name the stack would build."""

        names: set[str] = set()
        for layer in self._layers:
            names |= _unit_names(layer.name, layer.split, layer.names)
        return names

    @property
    def names(self) -> list[str]:
        """The layer names, top-to-bottom (excluding the model top)."""

        return [layer.name for layer in self._layers]

    def _index(self, name: str) -> int:
        """The position of the layer named ``name`` (raises ``KeyError`` if absent)."""

        for i, layer in enumerate(self._layers):
            if layer.name == name:
                return i
        raise KeyError(f"no layer named {name!r} (have {self.names}).")

    # -- compilation ------------------------------------------------------ #
    def _expanded_layers(self) -> tuple[list[_Layer], dict[str, list[int]]]:
        """Expand each declared UNIT into the model layers it becomes.

        Returns the flat layer list and the ``unit name -> layer indices`` map
        that :attr:`LayerBuildResult.units` carries forward. A unit with no
        ``split`` passes through unchanged and keeps its own name, so declaring
        ``split=1`` -- or adding a split later -- never renames anything. A unit
        carrying ``names`` uses those instead of the generated ones; either way
        the UNIT key in the returned map is the declared ``name``.
        """

        out: list[_Layer] = []
        units: dict[str, list[int]] = {}
        for layer in self._layers:
            shares = _split_shares(layer.split, layer.name)
            sub_names = _split_names(layer.names, shares, layer.name)
            first = len(out)
            if len(shares) == 1:
                out.append(_dc_replace(layer, split=None, names=None))
            else:
                steps = _split_steps(shares)
                for i, (share, step) in enumerate(zip(shares, steps, strict=True), 1):
                    out.append(_dc_replace(
                        layer,
                        name=sub_names[i - 1] if sub_names else f"{layer.name}_{i}",
                        surface=_sub_surface(layer.surface, share, step, layer.name),
                        split=None,
                        names=None,
                    ))
            units[layer.name] = list(range(first, len(out)))

        seen: dict[str, str] = {}
        for unit, indices in units.items():
            for i in indices:
                clash = seen.get(out[i].name)
                if clash is not None:
                    raise ValueError(
                        f"splitting {unit!r} produces the layer name "
                        f"{out[i].name!r}, which unit {clash!r} already uses. "
                        f"Layer names identify a contact, a DataFrame column and "
                        f"a 3-D scene actor, so they have to be distinct -- "
                        f"rename one of the two units."
                    )
                seen[out[i].name] = unit
        return out, units

    def _layer_surfaces(self) -> LayerSurfaces:
        """Compile to the :class:`LayerSurfaces` engine: the top plus each layer bottom, labeled."""

        layers, _ = self._expanded_layers()
        surfaces = [self._top] + [layer.surface for layer in layers]
        labels = ["top"] + [layer.name for layer in layers]
        return LayerSurfaces(surfaces, labels=labels)

    def _per_layer_config(self, default_min_thickness, default_pinch):
        """Resolve each layer's ``(min_thickness, pinch)``, filling unset ones with the defaults.

        Per MODEL layer, but the values come from the UNIT: every slice of a
        split unit inherits the unit's policy undivided, which is what makes the
        pinch verdict a property of the unit rather than of an arbitrary
        discretization choice.
        """

        layers, _ = self._expanded_layers()
        min_thk = [
            default_min_thickness if layer.min_thickness is None else layer.min_thickness
            for layer in layers
        ]
        pinch = [
            default_pinch if layer.pinch is None else layer.pinch
            for layer in layers
        ]
        return min_thk, pinch

    def _layer_shares(self) -> list[float]:
        """Each MODEL layer's share of its unit, in :meth:`_expanded_layers` order."""

        return [
            share
            for layer in self._layers
            for share in _split_shares(layer.split, layer.name)
        ]

    def _reconcile_separations(self, min_thk, pinch, min_sep: float, trigger_sep):
        """Per-layer ``(min_sep, trigger_sep)`` for reconcile, from each layer's policy.

        One rule: **reconcile never reaches a unit at or above its own
        ``min_thickness``** (ledger 377). The engine SETS every contact at most
        ``trigger_sep`` apart to exactly ``min_sep``, so a trigger above a unit's
        minimum thins a unit that is present.

        A ``"floor"`` layer declares a real minimum thickness, so reconcile enforces
        **its own** ``min_thickness`` -- that is what makes ``min_thickness`` set
        geometry rather than only decide an idomain value. Its separation and its
        trigger are both that minimum, so reconcile only ever OPENS it. An explicit
        ``trigger_sep`` replaces a floor layer's trigger, as documented on
        :meth:`build`; it then holds a thicker unit at its minimum too.

        A ``"passthrough"`` or ``"inactive"`` layer must be allowed to come out too
        thin, because being too thin is the signal that it pinches out. Reconcile
        gives those only the fallback ``min_sep``, just enough to keep the stack
        ordered, and ``min_thickness`` stays a pure threshold for them. Its trigger
        is the stack trigger (``trigger_sep``, else ``1.0 / <largest split>``) but
        never more than :data:`_UNDER_THE_MINIMUM` of its own minimum: a 0.8 ft clay
        with ``min_thickness=0.5`` was crushed to 0.1 ft under the 1.0 trigger and
        switched off. A unit below its minimum is pinched exactly as before, and
        crushed to ``min_sep`` only where it is also within the stack trigger: a
        1.5 ft clay with ``min_thickness=2`` keeps its 1.5 ft and reads idomain 0.

        ``min_thickness`` is declared for the UNIT, so a split unit's minimum is
        divided among its slices BY SHARE: three equal slices of a 3 ft minimum get
        1 ft each and the unit still comes out 3 ft. Handing every slice the unit's
        own figure would inflate the unit by its split factor -- the same failure
        ledger 158 records for the pinch threshold -- and dividing it equally held
        a 0.2 ft slice of a ``split=[0.9, 0.1]`` unit at 1.0 ft, pushing the unit's
        base 0.8 ft down.

        The stack trigger is scaled by the largest split because the engine's
        1.0 length unit is sized for geologic units, and applied unchanged to their
        SLICES it was destructive: a 2.4 ft unit split three ways came back
        ``[0.1, 1.5, 0.1]`` (ledger 158). An unsplit stack keeps the historical 1.0.
        """

        stack_trigger = (
            float(trigger_sep) if trigger_sep is not None else 1.0 / self._max_split()
        )
        seps, triggers = [], []
        for thickness, policy, share in zip(min_thk, pinch, self._layer_shares(), strict=True):
            own = float(thickness) * share
            if policy == "floor":
                seps.append(own)
                triggers.append(own if trigger_sep is None else stack_trigger)
            else:
                seps.append(float(min_sep))
                triggers.append(min(stack_trigger, own * _UNDER_THE_MINIMUM))
        return seps, triggers

    def _max_split(self) -> int:
        """The largest number of model layers any one unit becomes."""

        return max(
            (len(_split_shares(layer.split, layer.name)) for layer in self._layers),
            default=1,
        )

    def _named_surfaces(self) -> list[tuple[str, Surface]]:
        """``(key, surface)`` for the top, each unit's bottom and each excavation's bed.

        An excavation's key is ``"excavation:<name>"``, so it cannot collide with a
        unit's.
        """

        named = [("top", self._top)] + [(layer.name, layer.surface) for layer in self._layers]
        named += [
            (f"excavation:{excavation.name}", excavation.bottom)
            for excavation in self._excavations
            if excavation.bottom is not None
        ]
        return named

    def refresh(self) -> LayerStack:
        """Rebuild the cached raster of every derived (contour) surface now.

        Excavation beds included: a pond's bed is usually a ``Contours``.

        Returns
        -------
        LayerStack
            ``self``.

        Examples
        --------
        >>> import myflopy as mf
        >>> mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).refresh().names
        ['sand']
        """
        for _name, surface in self._named_surfaces():
            if getattr(surface, "is_derived", False):
                surface.resolve_source(refresh=True, warn=False)
        return self

    def cache_status(self) -> dict[str, str]:
        """Map ``name -> "missing"/"fresh"/"stale"`` for each derived surface.

        Returns
        -------
        dict of str to str
            Keyed by ``"top"``, a unit name, or ``"excavation:<name>"`` for an
            excavation's bed -- prefixed, so a unit and an excavation can never
            share a key. Only derived (contour) surfaces appear.

        Examples
        --------
        >>> import myflopy as mf
        >>> mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).cache_status()
        {}
        """
        return {
            name: surface.cache_status()
            for name, surface in self._named_surfaces()
            if getattr(surface, "is_derived", False)
        }

    def build(
        self,
        vor=None,
        *,
        default_min_thickness: float = 1.0,
        default_pinch: str = "floor",
        reconcile="bottom",
        min_sep: float = 0.1,
        trigger_sep: float | None = None,
        method: str = "area",
        refresh: bool = False,
        attach: bool = False,
    ) -> LayerBuildResult:
        """Resolve, reconcile, and pinch out the stack into DISV-ready arrays.

        Samples every surface onto the grid, reconciles crossing surfaces, applies
        the per-layer pinch policy, and returns the ``top``/``botm``/``idomain``
        arrays ready for ``mf.disv(...)``. A stale derived (contour) surface is
        reused with a warning; pass ``refresh=True`` to rebuild its cache first.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            The grid to sample onto. Required only if this stack was declared
            without one (``LayerStack(top=...)``); an explicit grid here wins
            over the stack's own. See :meth:`for_grid` to bind one instead.
        default_min_thickness : float, default 1.0
            Minimum layer thickness for layers that do not declare their own with
            ``.add(..., min_thickness=)``.
        default_pinch : str, default "floor"
            Thin-layer policy for layers that do not declare their own: ``"floor"``
            (hold the cell open at its minimum thickness), ``"inactive"`` (idomain
            0, a true pinch-out), or ``"passthrough"`` (idomain -1, flow passes
            vertically through and the cell can carry NO boundary condition).
        reconcile : {'bottom', 'top', True, False, None}, default 'bottom'
            What to do where two surfaces cross or come closer than 1 length unit
            of each other -- interpolated contacts routinely do, and a crossing
            means a negative layer thickness.

            ``'bottom'``
                Trust the surface ABOVE: push the lower one down to
                ``upper - min_sep``. **The usual choice**, and the only one that
                leaves your model top where you put it -- which matters when the
                top is a measured DEM. Corrections cascade downward, so one
                intruding contact pushes every contact below it clear.
            ``'top'``
                Trust the surface BELOW: raise the upper one to
                ``lower + min_sep``. Corrections propagate upward and CAN MOVE
                THE MODEL TOP -- with ``top=100, a=105`` it returns ``top=105.1``.
                Reach for it only when the deeper contact is the reliable one
                (say a well-picked bedrock surface under a coarse interpolated
                top), and check ``qc()`` afterwards.
            ``True``
                Same as ``'bottom'``.
            ``False`` or ``None``
                Do not reconcile. Crossing surfaces stay crossed, so layers get
                zero or negative thickness and the pinch policy is what saves
                you -- useful for seeing the raw sampled surfaces, not for
                building a model.

            Either way ``qc()`` reports, per layer, how many cells reconcile had
            to move and the largest move, which is how you tell "tidied two
            cells" from "rebuilt the geometry".
        min_sep : float, default 0.1
            **Fallback** separation, used only for layers that pinch
            (``"passthrough"`` / ``"inactive"``) -- enough to keep the stack ordered
            while still letting them come out thin enough to pinch. A ``"floor"``
            layer ignores it and uses its own ``min_thickness`` instead (a split
            unit's divided among its slices by share), which is how
            ``min_thickness`` comes to set geometry.
        trigger_sep : float, optional
            The separation at or below which reconcile ACTS -- surfaces this close
            are treated as conflicting even if they never actually cross, and the
            layer between them is SET to its separation. Defaults to
            ``1.0 / <largest split>``: the historical 1 length unit for an
            unsplit stack, and proportionally smaller once a unit is cut into
            slices, so that splitting a unit does not shrink it. Pass a number to
            override. It never reaches a pinching unit at or above its own
            ``min_thickness``: such a layer's trigger is capped just under its
            (share of its) minimum, so a 0.8 ft clay with ``min_thickness=0.5``
            stays 0.8 ft and active (ledger 377). A unit below its minimum
            is pinched, and crushed to ``min_sep`` where it is also within this
            trigger -- a 1.5 ft clay with ``min_thickness=2`` keeps its 1.5 ft and
            is pinched without being moved. A ``"floor"`` layer's
            trigger is its own minimum when this is omitted, and this number
            when given -- so a large one holds a thicker floor unit at its
            minimum too.
        method : str, default "area"
            Raster sampling method: ``"area"`` (area-weighted) or ``"centroid"``.
        refresh : bool, default False
            Rebuild any stale derived (contour) surface caches before sampling.
        attach : bool, default False
            Also publish the result onto ``vor.gdf_topbtm`` (see
            :meth:`LayerBuildResult.attach_to_grid`) so surface-aware builders --
            ``mf.sfr`` reach tops, ``mf.lak`` lake-cell layering -- can read the
            elevations, and the excavation records onto ``vor.excavations``. Use it
            when your model has SFR/LAK on this grid.

        Returns
        -------
        LayerBuildResult
            Bundles ``top`` ``(ncpl,)``, ``botm`` / ``idomain`` ``(nlay, ncpl)``,
            derived ``thickness``, and per-layer metadata -- with every declared
            :meth:`excavate` cut in, AFTER reconcile and the pinch policy, and
            ``excavations``/``ground``/``excavated``/``absorbed`` saying what was dug.
            ``refresh=True`` rebuilds an excavation bed's contour cache too.

        Raises
        ------
        ValueError
            If no layers have been added with :meth:`add`; any build-time refusal
            of :meth:`excavate` (a footprint on no cell, overlapping excavations,
            a bed with no data, a column with no active cell under its bed).

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> layers = (mf.LayerStack(top=mf.Flat(100.0))
        ...           .add("sand", thickness=20)
        ...           .add("clay", bottom=mf.Flat(79.5), pinch="inactive")
        ...           .build(vor, attach=True))
        >>> layers.nlay, layers.idomain[:, 0].tolist()    # clay is 0.5 ft < 1.0: pinched
        (2, [1, 0])
        >>> spec = mf.disv(layers)                         # the DISV from these arrays

        With a GIS top it reads ``mf.LayerStack(top=mf.Raster("ground.tif"))`` and a
        ``bottom=mf.Contours("base.gpkg", z="elev")``; nothing else changes.
        """
        vor = self._require_grid(vor, "build")
        result = self._build(
            vor,
            default_min_thickness=default_min_thickness, default_pinch=default_pinch,
            reconcile=reconcile, min_sep=min_sep, trigger_sep=trigger_sep,
            method=method, refresh=refresh,
        )
        if attach:
            result.attach_to_grid(vor)
        return result

    def _build(
        self, vor, *, default_min_thickness, default_pinch, reconcile, min_sep,
        trigger_sep, method, refresh, excavate: bool = True,
    ) -> LayerBuildResult:
        """:meth:`build` without the attach; ``excavate=False`` stops before the cut.

        The cut is the LAST step (design section 3.1): sampling, reconcile, the pinch
        invariant and the idomain all run on the un-excavated stack, then
        :meth:`LayerBuildResult._excavated_by` cuts every declared excavation in.
        ``qc`` asks for the pre-cut result to measure reconcile against.
        """
        if not self._layers:
            raise ValueError("Add at least one layer with .add(...) before build().")
        layers, units = self._expanded_layers()
        ls = self._layer_surfaces()
        rec_on, which = _reconcile_args(reconcile)
        min_thk, pinch = self._per_layer_config(default_min_thickness, default_pinch)
        seps, triggers = self._reconcile_separations(min_thk, pinch, min_sep, trigger_sep)
        gdf = ls.sample(
            vor, reconcile=rec_on, which=which, min_sep=seps,
            trigger_sep=triggers,
            method=method, length_units=self.length_units, refresh=refresh,
        )
        top, botm = ls._split_top_botm(gdf)
        thickness = ls._thickness(gdf)
        ls._validate_pinch_invariant(
            min_thk, pinch, thickness.shape[0],
            {"reconcile": rec_on, "min_sep": seps},
        )
        idomain = ls._idomain_from_thickness(
            _unit_thickness(thickness, units), min_thk, pinch
        )
        result = LayerBuildResult(
            top=top, botm=botm, idomain=idomain, thickness=thickness,
            names=[layer.name for layer in layers],
            min_thickness=min_thk, pinch=pinch,
            length_units=self.length_units, time_units=self.time_units,
            vor=vor, units=units,
        )
        if excavate and self._excavations:
            result = result._excavated_by(self._excavations, refresh=refresh)
        return result

    def qc(
        self,
        vor=None,
        *,
        default_min_thickness: float = 1.0,
        default_pinch: str = "floor",
        reconcile="bottom",
        min_sep: float = 0.1,
        trigger_sep: float | None = None,
        method: str = "area",
        refresh: bool = False,
    ) -> LayerQCReport:
        """Build the stack and run QC, including reconcile diagnostics.

        On top of :meth:`LayerBuildResult.qc` (NaN coverage, isolated active
        cells, thickness), this samples the surfaces *without* reconciling and
        reports, per layer, how many cells reconcile had to move and the largest
        move -- showing where surfaces were crossing before reconcile fixed them.

        Read this before trusting a build. A large ``reconcile_moved`` count
        means the geometry you got is substantially not the geometry you
        described, which is worth knowing before it becomes an idomain.

        Every build option is passed straight through to :meth:`build`, so QC
        reports on the geometry you are actually going to build -- including every
        declared :meth:`excavate`.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            Required only if the stack was declared without a grid.
        default_min_thickness : float, default 1.0
            As on :meth:`build`.
        default_pinch : str, default "floor"
            As on :meth:`build`.
        reconcile : {'bottom', 'top', True, False, None}, default 'bottom'
            As on :meth:`build`, where each option is spelled out -- they are what
            this report is diagnosing.
        min_sep : float, default 0.1
            As on :meth:`build`.
        trigger_sep : float, optional
            As on :meth:`build`.
        method : str, default "area"
            As on :meth:`build`.
        refresh : bool, default False
            As on :meth:`build` (excavation beds included).

        Returns
        -------
        LayerQCReport
            The excavated result's report (excavated cells counted on their own
            line), with ``reconcile_adjusted``/``reconcile_max_shift`` measured
            against the PRE-cut surfaces.

        Raises
        ------
        ValueError
            No grid; any refusal of :meth:`build`.

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> stack = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0).add("clay", thickness=20.0)
        >>> report = stack.excavate(75.0, mf.Area(vor.gdf_vorPolys.geometry.loc[4])).qc(vor)
        >>> report.reconcile_adjusted, report.excavated
        ([0, 0], [1, 0])
        """
        vor = self._require_grid(vor, "qc")
        precut = self._build(
            vor,
            default_min_thickness=default_min_thickness, default_pinch=default_pinch,
            reconcile=reconcile, min_sep=min_sep, trigger_sep=trigger_sep,
            method=method, refresh=refresh, excavate=False,
        )
        # Reconcile is measured against the PRE-cut surfaces; the dug result would
        # read the excavation as "reconcile moved the contact" in every pond cell.
        result = precut._excavated_by(self._excavations, refresh=refresh)
        report = result.qc()
        ls = self._layer_surfaces()
        raw = ls.sample(
            vor, reconcile=False, method=method,
            length_units=self.length_units, refresh=refresh,
        )
        _, raw_botm = ls._split_top_botm(raw)
        adjusted, max_shift = [], []
        for k in range(precut.nlay):
            shift = np.abs(np.asarray(raw_botm[k], float) - np.asarray(precut.botm[k], float))
            finite = shift[np.isfinite(shift)]
            adjusted.append(int((finite > 1e-6).sum()))
            max_shift.append(float(finite.max()) if finite.size else 0.0)
        report.reconcile_adjusted = adjusted
        report.reconcile_max_shift = max_shift
        return report

    # -- views: build with defaults, then view the result ----------------- #
    @property
    def plot(self) -> StackPlots:
        """The plotting verbs for this stack: ``map``, ``section``, ``surface``.

        Builds the stack with default options and returns the result's namespace,
        so ``stack.plot.map()`` is ``stack.build().plot.map()``. Build with
        non-default options first if you need them.

        **Works without a grid.** A stack declared with none draws on a coarse
        :meth:`draft_grid` covering the surfaces' own extent, so you can look at
        the layering before committing to a mesh. One INFO line records the
        extent and cell count, because a blocky preview should not be mistaken
        for your model. The draft is cached, so repeated ``.plot`` calls reuse it.

        Only pictures do this. :meth:`build`, :meth:`qc` and :meth:`to_disv`
        still require a real grid -- an approximate picture is useful, invented
        model geometry is not. Control the draft with
        ``stack.for_grid(stack.draft_grid(cells=2000)).plot``.
        """

        if self.vor is not None:
            return self.build().plot

        if self._draft is None:
            self._draft = self.draft_grid()
            xmin, ymin, xmax, ymax = self._draft.gdf_vorPolys.total_bounds
            logger.info(
                "LayerStack.plot: no grid on this stack, drawing on a draft grid "
                "of %d cells over (%.0f, %.0f)-(%.0f, %.0f). Pictures only -- "
                "build()/qc()/to_disv() still need a real grid.",
                int(self._draft.ncpl), xmin, ymin, xmax, ymax,
            )
        return self.build(self._draft).plot


    def to_disv(
        self,
        vor=None,
        *,
        name: str = "disv",
        attach: bool = True,
        default_min_thickness: float = 1.0,
        default_pinch: str = "floor",
        reconcile="bottom",
        min_sep: float = 0.1,
        trigger_sep: float | None = None,
        method: str = "area",
        refresh: bool = False,
    ):
        """Build a ready-to-use ``mf.disv`` spec (with pinch-out idomain) -- ``build(...).to_disv()``.

        The one-call alternative to ``build()`` + ``mf.disv(result)``. It is exactly
        that: the stack is resolved ONCE, by :meth:`build`, and the spec is
        :meth:`LayerBuildResult.to_disv` of that result, so every declared
        :meth:`excavate` is in it. It used to sample the stack a second time, which
        matched ``build()`` only by passing each per-layer separation along and which
        would have written the un-excavated top beside the excavated idomain;
        rerouting it moved no digest and no option (12 raster stacks, both sampling
        methods, three pinch policies, compared before and after).

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            The grid; required only if the stack was declared without one.
        name : str, default "disv"
            The spec's name -- the key ``with_packages(slot=)`` and a run's frozen
            snapshot use. It does not reach FloPy: the built package is always
            ``disv``, written to ``<model>.disv`` (as on ``mf.disv``).
        attach : bool, default True
            Publish the built elevations onto ``vor.gdf_topbtm`` and the excavation
            records onto ``vor.excavations`` (:meth:`LayerBuildResult.attach_to_grid`)
            -- note this defaults to True here and to False on :meth:`build`.
            Surface-aware builders (``mf.sfr`` reach tops, ``mf.lak``) and the
            layer-elevation hover rows read them.
        default_min_thickness : float, default 1.0
            As on :meth:`build`: the minimum thickness of a unit that declares none.
        default_pinch : str, default "floor"
            As on :meth:`build`: the pinch policy of a unit that declares none.
        reconcile : {'bottom', 'top', True, False, None}, default 'bottom'
            As on :meth:`build`, where each option is spelled out.
        min_sep : float, default 0.1
            As on :meth:`build`: the fallback separation for pinching layers.
        trigger_sep : float, optional
            As on :meth:`build`: the separation below which reconcile acts.
        method : str, default "area"
            As on :meth:`build`: raster sampling of the unit surfaces. An
            excavation's bed is always area-weighted.
        refresh : bool, default False
            As on :meth:`build`: rebuild stale contour caches, excavation beds
            included.

        Returns
        -------
        PackageSpec
            An EAGER DISV spec holding copies of the result's ``top``/``botm``/
            ``idomain``.

        Raises
        ------
        ValueError
            No layers added; no grid; any refusal of :meth:`build` (the excavation
            ones included).

        Examples
        --------
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
        >>> stack = mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
        >>> spec = stack.excavate(90.0, mf.Area(vor.gdf_vorPolys.geometry.loc[4])).to_disv(vor)
        >>> float(spec.options["top"][4]), list(vor.excavations)
        (90.0, ['excavation_0'])
        """
        if not self._layers:
            raise ValueError("Add at least one layer before to_disv().")
        vor = self._require_grid(vor, "to_disv")
        # ONE resolution, the result's own (design section 3.1). This used to sample
        # the stack a second time through `LayerSurfaces.to_disv` -- which could only
        # match `build()` by passing every per-layer separation along, and which, once
        # the stack is excavated, would write the UN-excavated top beside the
        # excavated idomain. Measured neutral for a stack with no excavation: the
        # spec's options and digest are unchanged on raster stacks under both
        # sampling methods.
        result = self.build(
            vor,
            default_min_thickness=default_min_thickness, default_pinch=default_pinch,
            reconcile=reconcile, min_sep=min_sep, trigger_sep=trigger_sep,
            method=method, refresh=refresh,
        )
        spec = result.to_disv(vor, name=name)
        if attach:
            # The arrays this spec was built from, with the INTEGER columns the
            # surface-aware builders read -- and the excavation records beside them.
            result.attach_to_grid(vor)
        return spec


__all__ = [
    "LayerStack",
    "LayerBuildResult",
    "LayerQCReport",
    "modflow_surfaces",
    "Raster",
    "Flat",
    "Contours",
    "Points",
    "Array",
    "Isopach",
    "Toward",
    "Min",
    "Max",
    "Clamp",
    "Where",
]
