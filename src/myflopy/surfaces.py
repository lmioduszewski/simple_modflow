"""Layer top/bottom surfaces for MODFLOW 6 discretization.

A :class:`Surface` describes one elevation surface -- a raster, a flat constant,
a surface offset below the one above it, an interpolation from drawn elevation
contours (via GRASS), or interpolation from scattered points. :class:`LayerSurfaces`
stacks them into ordered model layers, samples them onto a grid's cells, and
reconciles overlaps into DISV-ready top/bottom elevations. Thin layers can
optionally be pinched out (``idomain = -1``, vertical pass-through).

The grid does not own this logic. ``LayerSurfaces.attach(grid)`` writes the
result back to ``grid.gdf_topbtm`` so existing consumers (choropleth mapping,
layer-elevation lookups, ``mf.disv``) keep working unchanged.
"""

from __future__ import annotations

import numbers
import types
import warnings
from collections.abc import Sequence
from dataclasses import MISSING, dataclass, field, fields
from pathlib import Path
from typing import Any, NamedTuple

import geopandas as gpd
import numpy as np
import shapely as shp

from myflopy._defaults import DEFAULT_EPSG
from myflopy._length_units import (
    _declares_no_unit,
    _length_factor,
    _mf6_length_units,
)
from myflopy._logging import get_logger
from myflopy.areas import Area
from myflopy.modflow.mf6.grid.geometry import (
    BASIN_SHAPES,
    _per_contact,
    basin_depths,
    basin_radius,
    centroid_zones,
    footprint_cells,
    geometry_first,
    read_footprints,
    reconcile_surfaces,
    refuse_footprint_misses,
)
from myflopy.modflow.mf6.grid.surfaces import _area_weighted_sample, get_raster_vals_at_centroids
from myflopy.modflow.utils.datatypes.readers import _vector_layer_name

_logger = get_logger(__name__)

# Surfaces defined relative to the surface immediately above them.
_RELATIVE_KINDS = ("offset_below", "constant_thickness")

# Kinds that need `previous` but are NOT a fixed offset from it, so `values`
# handles them before the `is_relative` branch.
_PREVIOUS_KINDS = ("isopach", "toward")

# Kinds whose numbers are the surface's OWN elevations, and so the only ones a
# `units=` can describe: every other kind is relative (a distance from the surface
# above, in model units) or algebra over operands that each declare their own
# (ledger 405).
_UNIT_KINDS = ("flat", "raster", "contours", "points", "array")


def _depends_on_previous(surface) -> bool:
    """True if sampling ``surface`` consults the surface above it, at any depth.

    Composites carry it upward: ``Flat(150).capped_at(Isopach(map))`` is a `min`
    node, which looks absolute, over an operand that is not. `toward` needs this
    because a target measured from `previous` re-resolves against each cut it
    makes, drifting instead of dividing.
    """

    if surface.kind in _RELATIVE_KINDS or surface.kind in _PREVIOUS_KINDS:
        return True
    return any(
        _depends_on_previous(op) for op in surface.operands if op is not None
    )

# Pinch-out policies for thin cells (thickness < minimum_thickness):
#   passthrough -> idomain = -1 (cell removed, vertical flow passes through)
#   inactive    -> idomain = 0  (cell removed, blocks vertical flow)
#   floor       -> stay active  (reconcile enforces minimum spacing)
_PINCH_POLICIES = ("passthrough", "inactive", "floor")


def _as_per_layer(value, nlay: int, name: str) -> list:
    """Broadcast a scalar to one value per layer, or validate a per-layer sequence."""

    if np.isscalar(value):
        return [value] * nlay
    seq = list(value)
    if len(seq) != nlay:
        raise ValueError(
            f"{name} has {len(seq)} entries but there are {nlay} layers."
        )
    return seq


def _cells_in_zone(vor, zone) -> np.ndarray:
    """Boolean per-cell mask: True where the cell centroid lies inside ``zone``.

    ``zone`` is a shapely geometry (in the grid CRS) or a path to a vector file,
    whose rows are unioned into one zone. The zone is CLOSED -- a centroid exactly
    on its edge is inside it, through the shared
    :func:`~myflopy.modflow.mf6.grid.geometry.centroid_zones` rule. It used to be
    ``within``'s strict inside, so a zone whose edge ran along a row of centroids
    reached none of that row: ``box(0, 0, 100, 37.5)`` on a 4x4 grid of 25 ft cells
    selected cells 0-3 and not 4-7 (ledger 243).
    """

    if isinstance(zone, (str, Path)):
        gdf = gpd.read_file(zone)
        vor_crs = getattr(vor, "crs", None)
        if vor_crs is not None and gdf.crs is not None:
            gdf = gdf.to_crs(vor_crs)
        geom = gdf.union_all()
    else:
        geom = zone
    return centroid_zones(vor.gdf_vorPolys, [geom]) >= 0


def _content_keyed(payload: dict[str, Any], source_shas: dict) -> dict[str, Any]:
    """``payload`` with every path ``source_shas`` holds a hash for REMOVED, in place.

    The one statement of "the hash replaces the path" (ledger 293) that both carriers
    of a file inside a layer stack use -- a surface's raster or contours, and an
    excavation's footprint file. The hash itself rides beside the path under
    ``source_shas``, which the payload keeps, so the digest sees the file's CONTENT
    and not where it sits: moving a file keeps the digest, editing it moves the
    digest. A path the stamping could not read has no hash and stays, so two
    unreadable files never digest alike.
    """

    for name in source_shas:
        payload.pop(name, None)
    return payload


def _without_area_names(payload: dict[str, Any]) -> dict[str, Any]:
    """``payload`` with the ``name`` of every inline ``mf.Area`` footprint removed, in place.

    A name is a label, not content: an area taken from ``model.aoi`` carries one and
    the same polygon drawn inline does not, and the two must digest alike. Walks the
    ``operands`` too, which :meth:`Surface.to_dict` nests as plain dicts.
    """

    zone = payload.get("zone")
    if isinstance(zone, dict) and zone.get("kind") == "Area":
        zone.pop("name", None)
    for operand in payload.get("operands") or ():
        if isinstance(operand, dict):
            _without_area_names(operand)
    return payload


def _carries_propagate(node) -> bool:
    """Whether ``node``, or any surface under it, reads ``fill="propagate"``.

    The one walk behind two refusals: ``LayerStack.excavate``'s (a propagate source
    INSIDE algebra makes the stack and ``mf.lak`` measure different beds) and
    :meth:`Surface.basin`'s (its water surface may carry none at all).
    """

    return node.fill == "propagate" or any(
        _carries_propagate(op) for op in node.operands if op is not None
    )


def _carries_basin(node) -> bool:
    """Whether ``node``, or any surface under it, is a ``kind="basin"`` surface."""

    return node.kind == "basin" or any(
        _carries_basin(op) for op in node.operands if op is not None
    )


def _refuse_a_basin_gap(surface, values, *, label) -> None:
    """Refuse a LAYER surface that reads a basin's ``NaN`` off its lake (ledger 555).

    A basin has a bed only on its footprint's cells (:meth:`Surface.basin`); as a model
    top or a unit's bottom every other cell is ``NaN``, which nothing in the stack fills,
    and MF6 ran such a DISV to normal termination with every head ``NaN``. Composition
    that gives every cell a value -- ``mf.Min(<surface>, basin)``, whose ``nanmin`` takes
    the other surface off the lake -- passes, because only what reaches a cell is read.
    """

    missing = np.isnan(values)
    if _carries_basin(surface) and missing.any():
        cells = np.flatnonzero(missing)
        raise ValueError(
            f"layer surface {label!r} reads an mf.Basin, which has no bed off its lake's "
            f"cells -- {cells.size} of {values.size} cells here (first "
            f"{cells[:10].tolist()}). A layer elevation needs a value at every cell, "
            "and a NaN one runs in MF6 to normal termination with every head NaN. To cut "
            "a lake into the stack use stack.excavate(bottom=<basin>, within=<footprint>); "
            "to bake its bed into this surface give it one to fall back on, "
            "mf.Min(<surface>, <basin>)."
        )


def _describe_source(surface) -> str:
    """A few words naming where a layer surface's values come from, for a refusal."""

    if surface is None:
        return ""
    kind = getattr(surface, "kind", None)
    if kind == "raster" and surface.path is not None:
        return f"raster '{Path(surface.path).name}'"
    if kind == "contours" and surface.contours is not None:
        return f"contours '{Path(surface.contours).name}'"
    if surface.is_relative:
        return "measured from the surface above, which has no value there"
    return f"a {kind!r} surface"


def refuse_unbounded_cells(top, botm, idomain=None, *, labels, surfaces=None) -> None:
    """Refuse layer elevations that leave a cell of the DISV without a top or a bottom.

    The one rule every layer-stack route to a DISV goes through
    (``LayerBuildResult.to_disv`` -- so ``stack.to_disv``, ``mf.disv(stack)`` and
    ``mf.disv(result)`` -- and :meth:`LayerSurfaces.to_disv`): a cell with ``idomain != 0``
    whose top (the model top, or the bottom of the layer above) or bottom is not finite
    is refused, naming the elevations, the cells and their source. A raster covering
    only part of the grid as a unit bottom or a stack top wrote ``NaN`` there in silence
    (ledgers 555, 587). MF6 runs such a DISV to normal termination
    (``docs/reproducers/nan_layer_elevations.py``): every head ``NaN`` for an active
    cell, and a pass-through (``idomain -1``) cell silently disconnecting the cell below
    it, which kept its starting head. An INACTIVE cell (``idomain 0``) may hold one: MF6
    reads it and never uses it, and every head comes out as with a finite value there.

    Parameters
    ----------
    top : array_like of float, shape (ncpl,)
        The model top.
    botm : array_like of float, shape (nlay, ncpl)
        The layer bottoms, top layer first.
    idomain : array_like of int, shape (nlay, ncpl), optional
        ``None`` (the default) counts every cell.
    labels : sequence of str
        One name per elevation: the top's, then each layer bottom's.
    surfaces : sequence of Surface or None, optional
        The surface each elevation was sampled from, in ``labels`` order, to name the
        source; ``None`` (the default) names only the labels.

    Raises
    ------
    ValueError
        Naming each elevation with a missing value at a counted cell, how many cells
        and the first ten, its source, and the spellings that fill it.
    """

    top = np.asarray(top, dtype=float).reshape(-1)
    botm = np.asarray(botm, dtype=float).reshape(-1, top.size)
    elevations = np.vstack([top[None, :], botm])
    counted = (np.ones(botm.shape, dtype=bool) if idomain is None
               else np.asarray(idomain).reshape(botm.shape) != 0)
    # An elevation is needed where a counted cell has it as its top or its bottom.
    needed = np.zeros(elevations.shape, dtype=bool)
    needed[:-1] |= counted
    needed[1:] |= counted
    missing = needed & ~np.isfinite(elevations)
    if not missing.any():
        return
    sources = list(surfaces) if surfaces is not None else [None] * len(elevations)
    lines = []
    for index in np.flatnonzero(missing.any(axis=1)):
        cells = np.flatnonzero(missing[index])
        label = labels[index] if index < len(labels) else f"surface {index}"
        what = "the model top" if index == 0 else f"the bottom of layer {label!r}"
        source = _describe_source(sources[index] if index < len(sources) else None)
        lines.append(
            f"  {what}{f' ({source})' if source else ''}: {cells.size} of {top.size} cells, "
            f"first {cells[:10].tolist()}"
        )
    raise ValueError(
        "The layer elevations have no value at cells the DISV would make active or "
        "pass-through:\n" + "\n".join(lines) + "\n"
        "A raster covering only part of the grid leaves its cells NaN, and MF6 runs a NaN "
        "elevation to normal termination with NaN heads (or, under a pass-through cell, "
        "a cell silently disconnected). Fill them: a unit bottom from a raster takes the "
        "surface above where it has none with mf.Raster(path, fill='propagate'); any "
        "surface takes another's value where it is blank through mf.Max(<surface>, <other>) "
        "or mf.Min(...) (nanmax/nanmin) -- the model top has no surface above, so it needs "
        "a source covering the grid or such a fallback. Or make those cells inactive "
        "(idomain 0). result.qc() counts them."
    )


def _looks_positional(value) -> bool:
    """Whether ``value`` reads as cell numbers or a boolean mask rather than a place."""

    if isinstance(value, (range, np.ndarray)):
        return True
    if isinstance(value, (list, tuple)) and value:
        return all(isinstance(item, (numbers.Integral, bool, np.bool_)) for item in value)
    return False


def _read_footprint_file(path, *, layer=None, crs=None, noun: str = "excavation"):
    """Read a footprint FILE through the lake's reader, labelled for messages.

    One feature is labelled with the file's stem, as ``mf.lak`` labels it. Several are
    labelled ``"row <i>"`` -- the zero-based row of the layer read -- because the
    reader's own fallback, ``"<noun>_<i>"``, is the shape of a default excavation
    NAME, so "'excavation_1' (outside the grid)" read as another excavation missing.
    Shared by an excavation's footprint and a basin's (:meth:`Surface.basin`).
    """

    frame = read_footprints(path, layer=layer, crs=crs, noun=noun)
    if len(frame) > 1:
        frame.index = [f"row {i}" for i in range(len(frame))]
    return frame


def _coerce_footprint(within, *, layer, verb: str, noun: str):
    """``within`` as ``(zone, area, frame)``: a resolved file with its features, or an inline Area.

    The one reading of a footprint argument behind ``LayerStack.excavate(within=)`` and
    :meth:`Surface.basin`, so the two take the same spellings and refuse the same ones.
    A path must exist, is RESOLVED (ledger 293's two-bases trap) and is read NOW, so a
    line, a zero-area feature or a missing layer refuses at the call; ``frame`` is that
    read (``None`` for an inline footprint). An ``mf.Area`` is kept; a GeoDataFrame,
    GeoSeries or geometry is coerced to one; cell numbers or a mask are refused, since
    they name cells by position on one mesh.
    """

    zone = area = frame = None
    if isinstance(within, (str, Path)):
        path = Path(within)
        if not path.exists():
            raise FileNotFoundError(f"{verb}(within=) names a file that does not exist: {path}")
        zone = path.resolve()
        frame = _read_footprint_file(zone, layer=layer, noun=noun)
        return zone, area, frame
    if isinstance(within, Area):
        area = within
    elif isinstance(within, (gpd.GeoDataFrame, gpd.GeoSeries)) or hasattr(within, "geom_type"):
        area = Area.coerce(within)
    elif _looks_positional(within):
        raise TypeError(
            f"{verb}(within=) takes the footprint's POLYGONS -- a vector path, an "
            "mf.Area, a GeoDataFrame/GeoSeries or a shapely polygon -- not cell "
            "numbers or a boolean mask: those name cells by position on one mesh and "
            "would dig different cells after `with_grid`."
        )
    else:
        raise TypeError(
            f"{verb}(within=) got a {type(within).__name__}; pass a vector path, "
            "an mf.Area, a GeoDataFrame/GeoSeries or a shapely polygon."
        )
    geometry = area.geometry
    if geometry.geom_type not in ("Polygon", "MultiPolygon") or not geometry.area > 0:
        raise ValueError(
            f"{verb}(within=) must be polygonal with area; got a "
            f"{geometry.geom_type} of area {geometry.area:g}. A pond covers cells by "
            "AREA -- a line or point would dig the cells along it by length, which "
            "mf.lak refuses."
        )
    return zone, area, frame


def _footprint_frame(zone, area, *, layer, crs, name: str, noun: str):
    """The footprint as a frame in ``crs``: the file's features, or the Area's PARTS.

    An inline :class:`~myflopy.areas.Area` gives one feature per part, labelled
    ``name`` (one part) or ``"<name> part <i>"``; a file is read through
    :func:`_read_footprint_file`. Shared by an excavation and a basin at sampling time.
    """

    if area is not None:
        parts = [part for part in shp.get_parts(area.geometry) if not part.is_empty]
        labels = [name] if len(parts) == 1 else [f"{name} part {i}" for i in range(len(parts))]
        frame = gpd.GeoDataFrame({"_part": labels}, geometry=parts, crs=area.crs)
        return read_footprints(frame, crs=crs, id_field="_part", noun=noun)
    return _read_footprint_file(zone, layer=layer, crs=crs, noun=noun)


def _basin_number(value, *, label: str, column: str, what: str) -> float:
    """A basin's per-feature column value as a float, or a refusal naming the feature (ledger 247)."""

    try:
        number = float(value)
    except (TypeError, ValueError):
        number = float("nan")
    if not np.isfinite(number):
        raise ValueError(
            f"basin: feature {label!r} has {column}={value!r}, which is not a number -- a "
            f"NULL {what} is not a {what}; fill it or drop the feature."
        )
    return number


def _refuse_basin(surface) -> None:
    """Every refusal of a ``kind="basin"`` surface that reads no file.

    Run by ``Surface.__post_init__``, so :meth:`Surface.basin`, :meth:`Surface.from_dict`
    (a hand-written ``SimulationSpec.from_yaml`` payload) and a raw ``Surface(...)``
    refuse the same things.
    """

    if surface.shape not in BASIN_SHAPES:
        raise ValueError(
            f"Surface.basin(shape=) must be one of {list(BASIN_SHAPES)}, not {surface.shape!r}: "
            "'bowl' is a paraboloid (mean depth = depth / 2), 'cone' is linear (mean "
            "depth = depth / 3)."
        )
    if surface.units is not None:
        raise ValueError(
            f"units={surface.units!r} on a basin: its depth is a distance in model length "
            "units, and its water surface converts itself -- declare units= on the "
            "water_surface (e.g. mf.Raster(dem, units='meters'))."
        )
    if (surface.value is None) == (surface.depth_column is None):
        raise ValueError("a basin takes its depth once: a number (value) or a column (depth_column).")
    if surface.value is not None:
        depth = float(surface.value)
        if not (np.isfinite(depth) and depth > 0.0):
            raise ValueError(
                f"Surface.basin(depth=) must be a positive, finite length, not {surface.value!r}."
            )
    if (len(surface.operands) == 1) == (surface.surface_column is not None) or len(surface.operands) > 1:
        raise ValueError(
            "a basin takes its water surface once: one surface (operands) or a column "
            "(surface_column)."
        )
    if isinstance(surface.zone, Area):
        from_file = False
    elif isinstance(surface.zone, (str, Path)):
        from_file = True
    else:
        raise TypeError(
            f"a basin's footprint (zone) is a vector path or an mf.Area, not a "
            f"{type(surface.zone).__name__}."
        )
    if not from_file:
        for name, value in (
            ("layer", surface.zone_layer),
            ("depth", surface.depth_column),
            ("water_surface", surface.surface_column),
        ):
            if value is not None:
                raise ValueError(
                    f"Surface.basin({name}={value!r}) names "
                    + ("a layer of" if name == "layer" else "a column of")
                    + " the footprint FILE, and within= is not one -- an mf.Area or a "
                    "geometry has no layers and no attribute table. Pass the path."
                )
    if surface.operands:
        water = surface.operands[0]
        if _depends_on_previous(water):
            raise ValueError(
                f"Surface.basin(water_surface=) got a {water.kind!r} surface, which is "
                "measured from the surface ABOVE it. A basin's bed is measured from it, "
                "and the stack would take it from the ground while mf.lak takes it from "
                "the excavated top -- a different bed. Give the water surface as an "
                "absolute elevation: a number, a raster, or a column of the footprint."
            )
        if _carries_propagate(water):
            raise ValueError(
                "Surface.basin(water_surface=) carries a fill='propagate' source, which "
                "reads the surface ABOVE where it has no data -- the ground in the stack, "
                "the excavated bed in mf.lak. Fill the raster's gaps over the lake, or "
                "give the water surface as a number or a column."
            )


def _basin_column_refusal(argument: str, value, zone, frame) -> str:
    """The refusal text for a ``str`` argument of :meth:`Surface.basin` that names no column."""

    if frame is None:
        return (
            f"Surface.basin({argument}={value!r}) names no column: a column needs within= "
            "to be a FILE with an attribute table -- an mf.Area or a geometry has none."
        )
    names = ", ".join(str(c) for c in frame.columns if c not in ("geometry", "_source_name"))
    return f"Surface.basin({argument}={value!r}) is not a column of {zone.name} ({names})."


def _refuse_overlapping_basins(frame, *, name: str) -> None:
    """Refuse two features of one basin footprint that overlap by area, naming them.

    Each feature is profiled from its own shore, so where two overlap a point would
    sit at two depths -- the deeper is no more right than the sum.
    """

    geometries = list(frame.geometry)
    tree = shp.STRtree(geometries)
    left, right = tree.query(geometries, predicate="intersects")
    pairs = [
        (frame.index[i], frame.index[j])
        for i, j in zip(left, right, strict=True)
        if i < j and geometries[i].intersection(geometries[j]).area > 1e-9 * min(
            geometries[i].area, geometries[j].area
        )
    ]
    if pairs:
        shown = ", ".join(f"{a!r} and {b!r}" for a, b in pairs[:10])
        raise ValueError(
            f"basin {name!r}: features of its footprint overlap ({shown}). Each lake is "
            "profiled from its own shoreline, so an overlap has two depths; dissolve them "
            "into one feature or clip them apart."
        )


#: The most lattice points one basin feature is sampled on (see `_basin_cell_depths`).
_BASIN_MAX_PIXELS = 2_000_000


def _basin_cell_depths(cells, parts, *, shape: str, crs, label: str) -> np.ndarray:
    """The area-mean depth of each of ``cells`` (a slice of ``gdf_vorPolys``) under one feature.

    The profile is evaluated on a lattice over the cells' bounds -- the pixel the
    smaller of a hundredth of the feature's largest inscribed radius and a third of its
    smallest cell's width (each cell then holds at least ~9 pixels), measured to keep
    the lake's volume within 0.06% (ledger 554) -- and averaged per cell by the
    ONE area sampler, :func:`~myflopy.modflow.mf6.grid.surfaces._area_weighted_sample`,
    with depth 0 outside the lake. A lattice past ``_BASIN_MAX_PIXELS`` is coarsened to
    fit and says so at INFO.
    """

    from rasterio.io import MemoryFile
    from rasterio.transform import from_origin

    radius = max(basin_radius(part) for part, _depth in parts)
    width = float(np.sqrt(cells.area.to_numpy()).min())
    pixel = min(radius / 100.0, width / 3.0)
    xmin, ymin, xmax, ymax = (float(v) for v in cells.total_bounds)
    count = np.ceil((xmax - xmin) / pixel) * np.ceil((ymax - ymin) / pixel)
    if count > _BASIN_MAX_PIXELS:
        coarser = pixel * float(np.sqrt(count / _BASIN_MAX_PIXELS))
        _logger.info(
            "basin %s: sampling its depth on %.3g-unit pixels rather than %.3g, to keep "
            "the lattice under %d points", label, coarser, pixel, _BASIN_MAX_PIXELS,
        )
        pixel = coarser
    nx = max(int(np.ceil((xmax - xmin) / pixel)), 1)
    ny = max(int(np.ceil((ymax - ymin) / pixel)), 1)
    xs = xmin + pixel * (np.arange(nx) + 0.5)
    ys = ymax - pixel * (np.arange(ny) + 0.5)
    grid_x, grid_y = np.meshgrid(xs, ys)
    depth = basin_depths(grid_x.ravel(), grid_y.ravel(), parts, shape=shape).reshape(ny, nx)
    centres = np.array([geometry.centroid.coords[0] for geometry in cells.geometry])
    profile = dict(
        driver="GTiff", height=ny, width=nx, count=1, dtype="float64",
        transform=from_origin(xmin, ymax, pixel, pixel), crs=crs,
    )
    with MemoryFile() as memory:
        with memory.open(**profile) as sink:
            sink.write(depth, 1)
        with memory.open() as source, warnings.catch_warnings():
            # The lattice is drawn in the grid's own frame, so there is nothing to
            # reproject; on a grid with no CRS the sampler would say it is missing one.
            warnings.filterwarnings("ignore", message="CRS missing")
            return _area_weighted_sample(
                types.SimpleNamespace(gdf_vorPolys=cells), source, crs,
                centres[:, 0], centres[:, 1], f"basin {label}",
            )


#: Surface kinds whose content is per-cell arrays, so they cannot ride in a spec
#: (ledger 193 refused inlining arrays; ledger 269 names them here).
_UNSERIALIZABLE_SURFACES = frozenset({"array", "points"})


@dataclass(frozen=True)
class Surface:
    """One elevation surface, sampled onto a grid's cells when needed.

    Build one with a constructor rather than the raw fields: :meth:`raster`,
    :meth:`flat`, :meth:`offset_below`, :meth:`constant_thickness`,
    :meth:`from_contours`, :meth:`from_points`, or :meth:`basin` (an assumed lakebed)
    -- or combine surfaces with
    :meth:`minimum`, :meth:`maximum`, :meth:`clamp`, :meth:`where`, :meth:`isopach`,
    :meth:`shift`, :meth:`toward` and their fluent forms. Frozen: a constructor or an
    operation returns a NEW surface. Declaring one samples nothing; :meth:`values`
    does, normally called by a ``LayerStack``.

    Attributes
    ----------
    kind : str
        What the surface is, and so which other fields it reads: ``"flat"``,
        ``"raster"``, ``"contours"``, ``"points"``, ``"array"``, ``"offset_below"``,
        ``"constant_thickness"``, ``"basin"``, or one of the algebra kinds ``"min"``,
        ``"max"``, ``"clamp"``, ``"where"``, ``"isopach"``, ``"shift"``, ``"toward"``.
        Not validated when set, bar ``"basin"`` (checked as :meth:`basin` checks its
        arguments): an unknown kind raises ``TypeError`` when the surface is sampled.
    value : float, optional
        The number the kind needs: the elevation of a ``"flat"`` surface, the distance
        below the surface above for ``"offset_below"``/``"constant_thickness"``, the
        vertical displacement of a ``"shift"``, the fraction of a ``"toward"``, the
        maximum depth of a ``"basin"`` (``None`` when ``depth_column`` gives it).
        Unused by the other kinds.
    path : pathlib.Path, optional
        The raster file of a ``"raster"`` surface (set by :meth:`raster`).
    fill : {'propagate'}, optional
        For a raster or contour surface: ``"propagate"`` copies the surface above into
        cells the raster has no data for. ``None`` (the default) leaves them ``NaN``.
        Another value raises ``ValueError`` when sampled.
    nodata : float, optional
        For a raster or contour surface: a pixel value to treat as no-data,
        overriding the raster's own (see :meth:`raster`).
    units : str, optional
        The length unit the surface's elevations are in (``"meters"``, ``"feet"``,
        ``"ft"``, ``"m"``, ``"cm"``...). Converted by a single factor when sampled
        through something that declares the model's length unit (a ``LayerStack``
        with ``length_units=``, :meth:`LayerSurfaces.sample`, :func:`sample_bed`,
        ``mf.lak``'s bed). Sampled with NO model length unit (``length_units=None``,
        or a context's ``"unknown"``) it raises ``ValueError`` naming
        ``length_units=`` -- it was read as model units in silence (ledger
        433). Read wherever the surface sits: an
        operand inside an algebra surface (``minimum``, ``shift``, ``where``, ...) is
        converted before the operands combine (ledger 405). Only a ``flat``,
        ``raster``, ``contours``, ``points`` or ``array`` surface takes it -- on a
        relative or an algebra surface it raises ``ValueError`` there -- and an
        unknown name raises ``ValueError`` listing the known ones.
    contours : pathlib.Path, optional
        The contour-line vector file of a ``"contours"`` surface. This and the next
        nine fields are set by :meth:`from_contours`, which documents them.
    z : str, default "Elev"
        The contours' elevation attribute.
    region_raster : pathlib.Path, optional
        Raster whose extent is the GRASS region.
    region_vector : pathlib.Path, optional
        Polygon whose extent (and shape, with ``clip``) is the GRASS region.
    clip : bool, default True
        Null the interpolation outside the region.
    resolution : float, default 4
        Cell size of the interpolated raster, in the CRS's units.
    out : pathlib.Path, optional
        Where the interpolated GeoTIFF (the cache) is written.
    epsg : str, default "2927"
        EPSG code of the GRASS location.
    grass_kwargs : dict, optional
        Extra keywords for the GRASS interpolator; part of the cache signature.
    progress : bool, default False
        Show GRASS's own progress while interpolating; NOT part of the cache
        signature.
    points : tuple of three array_like, optional
        ``(xs, ys, zs)`` of a ``"points"`` surface (set by :meth:`from_points`).
    method : str, default "linear"
        The :func:`scipy.interpolate.griddata` method of a ``"points"`` surface. Not
        the raster sampling method, which is :meth:`values`' own ``method=``.
    array : numpy.ndarray, optional
        The per-cell elevations of an ``"array"`` surface (set by :meth:`from_array`).
    operands : tuple of Surface, optional
        The surfaces an algebra kind combines: two or more for ``"min"``/``"max"``;
        ``(surface, lower, upper)`` with ``None`` for an absent bound for ``"clamp"``;
        ``(inside, outside)`` for ``"where"``; one for ``"isopach"``, ``"shift"`` and
        ``"toward"``, and a ``"basin"``'s water surface (empty when
        ``surface_column`` gives it). Empty for every other kind.
    zone : shapely geometry, str, pathlib.Path or mf.Area, optional
        The zone of a ``"where"`` surface (see :meth:`where`), or a ``"basin"``'s
        footprint: its resolved vector file, or an inline :class:`~myflopy.areas.Area`.
    source_shas : dict, optional
        ``{field name: 16-character content hash}`` for the files this surface reads,
        stamped when it becomes part of a LAZY package (``mf.disv`` on a gridless
        stack, a lazy ``mf.evt``); empty otherwise. The digest reads it in place of the
        paths (:meth:`digest_payload`).
    depth_column : str, optional
        A ``"basin"``'s footprint column giving each feature its maximum depth (model
        length units), in place of ``value``.
    surface_column : str, optional
        A ``"basin"``'s footprint column giving each feature a flat water surface, in
        place of an operand.
    shape : {'bowl', 'cone'}, optional
        A ``"basin"``'s profile (see :meth:`basin`); ``None`` for every other kind.
    zone_layer : str, optional
        The layer of a ``"basin"``'s footprint file, by name; ``None`` reads the first.
    """

    kind: str
    #: Token :func:`myflopy.specs._json_value` writes this object under, and
    #: :func:`myflopy.specs._spec_value` revives it from. A `$`-prefixed key rather
    #: than a ``kind`` tag because ``Surface`` HAS a ``kind`` field of its own
    #: ("flat", "raster", "shift", ...) and the source-spec dispatcher keys on
    #: ``kind`` -- tagging one as the other makes a surface come back a plain dict.
    _SPEC_TOKEN = "$surface"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready payload, or raise for a surface that holds arrays.

        ``from_array`` and ``from_points`` are REFUSED by name. Their content is
        per-cell data, and inlining arrays into a spec is what ledger 193 refused --
        ~3.3 MB per model at 9,405 cells, and it floods `Project.validate`. They stay
        usable; they just cannot be written down, which is what a run's frozen
        snapshot is for (ledger 269).
        """

        if self.kind in _UNSERIALIZABLE_SURFACES:
            raise ValueError(
                f"A '{self.kind}' surface holds per-cell arrays and cannot be written "
                "into a spec. Declare the layer from a raster, contours or a constant "
                "so the stack can be rebuilt on another mesh, or keep the built arrays "
                "and let the run's frozen snapshot carry them."
            )
        payload: dict[str, Any] = {}
        for item in fields(self):
            value = getattr(self, item.name)
            # `default_factory` fields have `default is MISSING`, so comparing against
            # `item.default` alone wrote an empty `grass_kwargs: {}` into every payload.
            default = (
                item.default_factory()
                if item.default is MISSING and item.default_factory is not MISSING
                else item.default
            )
            if item.name != "kind" and type(value) is type(default) and value == default:
                continue
            if item.name == "operands":
                payload[item.name] = [operand.to_dict() for operand in value]
            elif isinstance(value, Area):
                # A basin's inline footprint: WKT + CRS, as an excavation's `area` is.
                payload[item.name] = value.to_dict()
            elif isinstance(value, Path):
                payload[item.name] = value.as_posix()
            else:
                payload[item.name] = value
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """The declaration with file CONTENT standing in for file LOCATION.

        Same rule as :meth:`myflopy.sources.BoundarySource.digest_payload`, and it
        exists for the same measured reason one carrier over (ledger 293). A surface
        reaches ``PackageSpec.inputs`` through a lazy ``mf.disv`` or ``mf.evt``, where
        the file it names is read at BUILD time -- so digesting the path meant
        ``save(bundle_sources=True)`` moved the package's digest on a byte-identical
        archive, and a different raster at the same path digested identically.

        A path is replaced only where :attr:`source_shas` has a hash for it, so a
        surface whose file could not be read still digests by path rather than
        colliding with every other unreadable one.

        ``out`` is dropped outright rather than hashed. It is the GRASS interpolation
        CACHE -- a file this surface WRITES -- and its content is fully determined by
        ``contours`` plus the build parameters, both of which are digested here. So it
        carries nothing the digest needs, while bundling repoints it like any other
        path and would move the digest for a file nobody edited. An inline footprint's
        ``name`` is dropped too: a label, not content.
        """

        payload = _without_area_names(self.to_dict())
        payload.pop("out", None)
        return _content_keyed(payload, self.source_shas)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Surface:
        """Rebuild a surface from :meth:`to_dict`'s payload.

        Parameters
        ----------
        data : dict
            A payload as :meth:`to_dict` writes it: field name to value, with
            ``"kind"`` required and every omitted field taking its default. The path
            fields (``path``, ``contours``, ``region_raster``, ``region_vector``,
            ``out``, and a basin's ``zone``) are revived as :class:`~pathlib.Path`, an
            inline footprint written by :class:`~myflopy.areas.Area`'s ``to_dict`` as an
            ``Area``, and each ``operands`` entry is rebuilt recursively. A ``"basin"``
            payload is checked as :meth:`basin` checks its arguments, short of reading
            the footprint file (whose path may be project-relative until the build
            anchors it).

        Returns
        -------
        Surface

        Raises
        ------
        ValueError
            ``data`` holds a key that is not a :class:`Surface` field, or a
            ``"basin"`` payload :meth:`basin` would refuse (an unknown shape, a depth
            that is not positive, a water surface measured from the surface above or
            carrying ``fill="propagate"``, a layer or column with an inline footprint).
        TypeError
            A ``"basin"`` payload whose ``zone`` is neither a path nor an Area.
        """

        known = {item.name for item in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"Unknown Surface fields: {sorted(unknown)}")
        values = dict(data)
        for name in ("path", "contours", "region_raster", "region_vector", "out"):
            if values.get(name) is not None:
                values[name] = Path(values[name])
        if "operands" in values:
            values["operands"] = tuple(cls.from_dict(item) for item in values["operands"])
        if isinstance(values.get("zone"), dict):
            values["zone"] = Area.from_dict(values["zone"])
        elif values.get("kind") == "basin" and isinstance(values.get("zone"), str):
            values["zone"] = Path(values["zone"])
        return cls(**values)

    value: float | None = None
    path: Path | None = None
    fill: str | None = None
    #: Overrides the raster's declared no-data value; see `Surface.raster`.
    nodata: float | None = None
    units: str | None = None  # source units, if different from the model units
    # contour-interpolation parameters
    contours: Path | None = None
    z: str = "Elev"
    region_raster: Path | None = None
    region_vector: Path | None = None
    clip: bool = True
    resolution: float = 4
    out: Path | None = None
    epsg: str = str(DEFAULT_EPSG)
    grass_kwargs: dict = field(default_factory=dict)
    #: Show GRASS's own per-step progress while interpolating. Deliberately NOT
    #: part of `grass_kwargs`: those feed the derived-raster cache signature, so
    #: routing a display flag through them would invalidate the cache and force a
    #: silent re-interpolation the first time anyone asked to watch.
    progress: bool = False
    # point-interpolation parameters
    points: tuple[Any, Any, Any] | None = None
    method: str = "linear"
    # direct per-cell elevation array (e.g. an existing model's top/botm)
    array: Any = None
    # surface-algebra / isopach parameters
    operands: tuple = ()
    zone: Any = None
    #: ``{field name: 16-char content hash}`` for the files this surface reads, filled
    #: in when the surface becomes part of a LAZY package (``mf.disv`` on a gridless
    #: stack, ``mf.evt(surface=...)``). Empty everywhere else: a bare
    #: ``Surface.raster(path)`` deliberately opens nothing, and hashing a multi-hundred-
    #: megabyte DEM at construction would make the cheapest call in the layer API the
    #: most expensive. It is the DIGEST's stand-in for the path (:meth:`digest_payload`)
    #: -- see ledger 293's second pass.
    source_shas: dict = field(default_factory=dict)
    # assumed-bathymetry (``kind="basin"``) parameters, set by :meth:`basin`
    #: The footprint column giving each feature its maximum depth.
    depth_column: str | None = None
    #: The footprint column giving each feature its water-surface elevation.
    surface_column: str | None = None
    #: The bed's profile, one of ``BASIN_SHAPES``.
    shape: str | None = None
    #: The vector layer of a basin's footprint file, by name.
    zone_layer: str | None = None

    def __post_init__(self) -> None:
        """Refuse a basin no grid could sample (:func:`_refuse_basin`); nothing else is checked."""

        if self.kind == "basin":
            _refuse_basin(self)

    @classmethod
    def raster(
        cls,
        path: Path | str,
        *,
        fill: str | None = None,
        units: str | None = None,
        nodata: float | None = None,
    ) -> Surface:
        """A surface sampled from an existing raster.

        Parameters
        ----------
        path : Path or str
            The raster to sample.
        fill : {'propagate'}, optional
            ``"propagate"`` carries the surface above down into cells where the
            raster has no data (nodata or outside coverage), so the layer
            pinches against the surface above instead of leaving a ``NaN``
            elevation.
        units : str, optional
            e.g. ``"meters"``. Converts the sampled elevations to the model
            length units when sampled through a stack that declares them.
        nodata : float, optional
            Treat this value as no-data, overriding whatever the file declares.

            For rasters clipped to a boundary but written WITHOUT a no-data
            value in the header -- a very common GIS export -- where the fill
            (usually ``0``) is then read as a real elevation. Symptoms: a
            surface that is rectangular rather than domain-shaped, ground at
            elevation zero outside the boundary, and every layer beneath it
            collapsing to ``min_sep`` once the contacts get capped to it.

            Masking happens BEFORE the per-cell average, so a cell straddling
            the edge is the mean of its valid pixels rather than a blend with
            the sentinel.

            ``mf.Raster(dem, nodata=0)`` is the fix without rewriting the file;
            declaring it properly in the raster (``gdalwarp -dstnodata``) is the
            fix that helps every other tool too.

        Returns
        -------
        Surface

        Examples
        --------
        >>> mf.Raster("ground.tif")                       # header nodata is honoured
        >>> mf.Raster("clipped.tif", nodata=0)            # 0 is fill, not ground
        >>> mf.Raster("bedrock.tif", fill="propagate")    # gaps inherit from above
        """

        return cls(
            kind="raster", path=Path(path), fill=fill, units=units, nodata=nodata
        )

    @classmethod
    def flat(cls, value: float, *, units: str | None = None) -> Surface:
        """A flat, constant-elevation surface (replaces the magic integer).

        Parameters
        ----------
        value : float
            The elevation, the same at every cell. Anything ``float()`` accepts,
            numpy scalars included.
        units : str, optional
            Length unit of ``value`` (``"meters"``, ``"feet"``, ...), converted to the
            model's when sampled through a stack that declares ``length_units``.
            ``None`` (the default) takes ``value`` as already in model units.

        Returns
        -------
        Surface
            A ``kind="flat"`` surface.
        """

        return cls(kind="flat", value=float(value), units=units)

    @classmethod
    def offset_below(cls, distance: float) -> Surface:
        """A surface a fixed vertical distance below the surface above it.

        Parameters
        ----------
        distance : float
            How far below the surface immediately above, in model length units, per
            cell: the result is ``previous - distance``. Relative, so it cannot be the
            model top (``ValueError`` when sampled) and takes no ``units=``.

        Returns
        -------
        Surface
            A ``kind="offset_below"`` surface.
        """

        return cls(kind="offset_below", value=float(distance))

    @classmethod
    def constant_thickness(cls, thickness: float) -> Surface:
        """A surface giving the layer above it a fixed thickness.

        Numerically identical to :meth:`offset_below`, but named for the common
        "this layer is N units thick" intent.

        Parameters
        ----------
        thickness : float
            The layer's thickness in model length units: the result is
            ``previous - thickness`` per cell. Relative, so it cannot be the model top
            (``ValueError`` when sampled).

        Returns
        -------
        Surface
            A ``kind="constant_thickness"`` surface.
        """

        return cls(kind="constant_thickness", value=float(thickness))

    @classmethod
    def from_contours(
        cls,
        contours: Path | str,
        *,
        z: str = "Elev",
        region_raster: Path | str | None = None,
        region_vector: Path | str | None = None,
        resolution: float = 4,
        out: Path | str | None = None,
        epsg: str = str(DEFAULT_EPSG),
        units: str | None = None,
        clip: bool = True,
        progress: bool = False,
        **grass_kwargs: Any,
    ) -> Surface:
        """A surface interpolated from drawn elevation contours via GRASS.

        The interpolated raster is written to ``out`` (default: alongside the
        contours with a ``.interp.tif`` suffix) the first time the surface is
        resolved, then reused -- interpolate once, sample many times.

        **When does the interpolation run?** Not here. Declaring the surface is
        free; the first call that needs its values -- or merely its extent, which
        includes ``LayerStack.draft_grid`` and so a gridless ``stack.plot`` --
        resolves it. After that the GeoTIFF is reused across builds, sessions and
        processes. :meth:`cache_status` reports ``fresh``/stale/absent, and
        ``refresh=True`` on ``build``/``qc``/``to_disv`` forces a rebuild.

        Parameters
        ----------
        contours : Path or str
            Vector file of labelled elevation contour lines. Stored as a ``Path``.
        z : str, default "Elev"
            The attribute holding each contour's elevation.
        region_raster : Path or str, optional
            Raster whose extent becomes the GRASS region. **Exactly one of this
            or ``region_vector`` is required** -- without a region GRASS has no
            extent to work in. Stored as a ``Path``.

            Note it uses the WHOLE raster. A DEM covering more ground than your
            model makes every interpolation proportionally slower, and the extra
            area is the worst case for ``r.surf.contour`` because it has no
            contours to work from out there. Prefer ``region_vector``.
        region_vector : Path or str, optional
            Polygon whose extent becomes the GRASS region -- usually your model
            domain, which is normally what you want interpolated and no more.
            Stored as a ``Path``: a ``str`` was kept as given while
            :meth:`from_dict` revives a ``Path``, so a saved-and-reloaded contour
            surface compared ``!=`` the original with an identical payload (found
            2026-09-25).
        resolution : float, default 4
            Output cell size, in the CRS's units.
        out : Path or str, optional
            Where the interpolated GeoTIFF lands. Defaults to the contours path
            with a ``.interp.tif`` suffix. Stored as a ``Path``.
        epsg : str, default "2927"
            EPSG code for the GRASS location.
        units : str, optional
            Length unit of the elevations, converted to the model's if they differ.
        clip : bool, default True
            Null out everything beyond the region. With ``region_vector`` that
            means the POLYGON, not merely its bounding box -- a domain-shaped
            result, with the bbox corners left as nodata. Already what you want
            in nearly every case; ``False`` keeps the full rectangle, which is
            occasionally useful when a neighbouring model needs the overlap.
        progress : bool, default False
            Print GRASS's own per-step percentages while it works.

            GRASS modules are subprocesses that write to the inherited file
            descriptors, so **in a terminal you already see them** and this flag
            changes nothing. In a **Jupyter** notebook they go to the kernel's
            console instead of the cell, and `contextlib.redirect_stdout` cannot
            catch them because the writing happens below Python. This captures
            them at the fd level and re-emits them through Python, so they land
            in the cell.

            Off by default: library code logs rather than prints, and the flag is
            the human asking. It is NOT part of the cache signature, so turning it
            on does not re-interpolate anything.
        **grass_kwargs
            Forwarded to :class:`~myflopy.modflow.utils.contour_interp
            .ContourSurfaceInterpolator`, whose remaining keywords are
            ``grass_bin`` (the GRASS launcher), ``grassdata`` (the GRASS database
            directory, default ``~/grassdata``) and ``location`` (default
            ``"myflopy_interp"``). These DO feed the cache signature, so changing
            one invalidates the cached raster. A name it does not take is accepted
            HERE and refused only when the surface is first resolved, with a
            ``TypeError`` from that class (measured with ``extra=(1, 2)``; found
            2026-09-25, not changed).

        Returns
        -------
        Surface
            A ``kind="contours"`` surface. Nothing is read or interpolated yet; see
            "When does the interpolation run?" above.

        Raises
        ------
        ValueError
            Neither or both of ``region_raster`` and ``region_vector``.

        Examples
        --------
        >>> import myflopy as mf
        >>> ground = mf.Contours("ground.gpkg", z="Elev", region_vector="domain.gpkg")
        >>> ground.kind, ground.region_vector.as_posix(), ground.out.as_posix()
        ('contours', 'domain.gpkg', 'ground.interp.tif')
        >>> ground == mf.Surface.from_dict(ground.to_dict())   # what a reload builds
        True
        >>> watch = mf.Contours("clay.gpkg", region_vector="domain.gpkg", progress=True)
        """

        # Checked here, not in `_set_region`: that runs after GRASS has started a
        # session, so the message arrives late and buried. This is a typo you can
        # catch while writing the line.
        if (region_raster is None) == (region_vector is None):
            both = region_raster is not None
            raise ValueError(
                "Contours need exactly one region, and got "
                + ("both" if both else "neither")
                + ". Pass region_vector=<your domain polygon> (usually what you "
                "want) or region_raster=<a raster whose extent to use>. GRASS has "
                "no extent to interpolate in without one."
            )

        contours = Path(contours)
        out = Path(out) if out is not None else contours.with_suffix(".interp.tif")
        return cls(
            kind="contours",
            contours=contours,
            z=z,
            # `Path`, as `contours` and `out` are and as `from_dict` revives them: a
            # `str` kept as given made a round-tripped surface compare unequal.
            region_raster=Path(region_raster) if region_raster is not None else None,
            region_vector=Path(region_vector) if region_vector is not None else None,
            clip=clip,
            resolution=resolution,
            out=out,
            epsg=epsg,
            units=units,
            grass_kwargs=grass_kwargs,
            progress=progress,
        )

    @classmethod
    def from_points(
        cls, xs, ys, zs, *, method: str = "linear", units: str | None = None
    ) -> Surface:
        """A surface interpolated from scattered (x, y, z) points (SciPy).

        Interpolated at each cell CENTROID with :func:`scipy.interpolate.griddata`
        when sampled. Holds arrays, so :meth:`to_dict` refuses it: a stack using one
        cannot be saved in a project spec.

        Parameters
        ----------
        xs, ys : array_like of float
            The points' coordinates, in the grid's own coordinate frame (nothing is
            reprojected). Converted with ``numpy.asarray``.
        zs : array_like of float
            The elevation at each point, one per ``(xs, ys)`` pair.
        method : {'linear', 'nearest', 'cubic'}, default 'linear'
            The ``griddata`` method. ``'linear'`` and ``'cubic'`` give ``NaN`` for a
            cell centroid outside the points' convex hull; ``'nearest'`` fills every
            cell.
        units : str, optional
            Length unit of ``zs``, converted to the model's when sampled through a
            stack that declares ``length_units``.

        Returns
        -------
        Surface
            A ``kind="points"`` surface.
        """

        return cls(
            kind="points",
            points=(np.asarray(xs), np.asarray(ys), np.asarray(zs)),
            method=method,
            units=units,
        )

    @classmethod
    def from_array(cls, values, *, units: str | None = None) -> Surface:
        """A surface whose per-cell elevations are supplied directly.

        ``values`` must have one entry per grid cell (``len == ncpl``); use this
        to reuse an existing model's top/botm arrays on the *same* grid. To move
        an existing model's surfaces onto a *different* grid, interpolate instead
        (:meth:`from_points`, or :meth:`LayerStack.from_modflow`).

        Parameters
        ----------
        values : array_like of float
            One elevation per cell, in cell order. Stored as a float array and
            flattened when sampled, so a ``(1, ncpl)`` array works too; any other
            total size raises ``ValueError`` when sampled. Holds an array, so
            :meth:`to_dict` refuses it.
        units : str, optional
            Length unit of ``values``, converted to the model's when sampled through a
            stack that declares ``length_units``.

        Returns
        -------
        Surface
            A ``kind="array"`` surface.
        """

        return cls(kind="array", array=np.asarray(values, dtype=float), units=units)

    @staticmethod
    def _coerce(x) -> Surface:
        """Coerce a Surface / path / number operand into a Surface."""
        if isinstance(x, Surface):
            return x
        if isinstance(x, (int, float)):
            return Surface.flat(float(x))
        if isinstance(x, (str, Path)):
            return Surface.raster(x)
        raise TypeError(f"Expected a Surface, path, or number; got {type(x).__name__}.")

    @staticmethod
    def _reject_lone_envelope(name: str, fluent: str, surfaces: tuple) -> None:
        """Refuse a one-surface envelope -- nearly always a misbound instance call.

        ``minimum``/``maximum`` are classmethods taking ``*surfaces``, so
        ``a.maximum(b)`` binds ``a`` to ``cls`` and arrives here as the single
        operand ``b``: a plausible-looking Surface that quietly ignores ``a``,
        and downstream a wrong contact elevation nobody sees until the heads
        look odd. The variadic signature is what hides it -- ``a.shift(-25)``
        raises for the missing ``distance``, but ``*surfaces`` swallows the
        arity error. An envelope of one surface is a no-op even when it is
        deliberate, so it costs nothing to make it the tripwire.
        """

        if len(surfaces) >= 2:
            return
        raise TypeError(
            f"Surface.{name}() takes two or more surfaces, got {len(surfaces)}. "
            f"If you wrote `a.{name}(b)`: {name} is a classmethod, so `a` binds to "
            f"`cls` and is dropped -- the result would be `b` alone. Write "
            f"`Surface.{name}(a, b)`, or the instance-bound `a.{fluent}(b)`."
        )

    @classmethod
    def minimum(cls, *surfaces) -> Surface:
        """Per-cell minimum (lower envelope) of two or more surfaces.

        A classmethod: call it as ``Surface.minimum(a, b)``, or reach for the
        instance-bound :meth:`capped_at`. One surface raises :class:`TypeError`,
        because that is what a misbound ``a.minimum(b)`` looks like from here.

        Parameters
        ----------
        *surfaces : Surface, str, pathlib.Path, int or float
            Two or more operands, each a :class:`Surface`, a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`), or a python ``int``/``float`` (read as
            :meth:`flat`). Anything else -- a ``numpy.int64`` included -- raises ``TypeError``.
            Sampled with ``numpy.nanmin``, so a cell where one operand is ``NaN`` takes the others'
            minimum.

        Returns
        -------
        Surface
            A ``kind="min"`` surface.

        Raises
        ------
        TypeError
            Fewer than two operands, or an operand of another type.
        """
        cls._reject_lone_envelope("minimum", "capped_at", surfaces)
        return cls(kind="min", operands=tuple(cls._coerce(s) for s in surfaces))

    @classmethod
    def maximum(cls, *surfaces) -> Surface:
        """Per-cell maximum (upper envelope) of two or more surfaces.

        A classmethod: call it as ``Surface.maximum(a, b)``, or reach for the
        instance-bound :meth:`floored_at`. One surface raises :class:`TypeError`,
        because that is what a misbound ``a.maximum(b)`` looks like from here.

        Parameters
        ----------
        *surfaces : Surface, str, pathlib.Path, int or float
            Two or more operands, each a :class:`Surface`, a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`), or a python ``int``/``float`` (read as
            :meth:`flat`). Anything else -- a ``numpy.int64`` included -- raises ``TypeError``.
            Sampled with ``numpy.nanmax``, so a cell where one operand is ``NaN`` takes the others'
            maximum.

        Returns
        -------
        Surface
            A ``kind="max"`` surface.

        Raises
        ------
        TypeError
            Fewer than two operands, or an operand of another type.
        """
        cls._reject_lone_envelope("maximum", "floored_at", surfaces)
        return cls(kind="max", operands=tuple(cls._coerce(s) for s in surfaces))

    @classmethod
    def clamp(cls, surface, *, lower=None, upper=None) -> Surface:
        """Constrain ``surface`` to stay within ``[lower, upper]`` (per cell).

        ``lower``/``upper`` may be surfaces, paths, or constants; either may be
        a relative surface (e.g. ``upper=Surface.offset_below(5)`` to keep this
        surface at least 5 below the one above). At least one is required: with
        neither, the result is ``surface`` unchanged, which is also exactly what
        a misbound ``a.clamp(b)`` produces, so both raise :class:`TypeError`.

        Parameters
        ----------
        surface : Surface, str, pathlib.Path, int or float
            The surface to constrain: a :class:`Surface`, a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`), or a python ``int``/``float`` (read as
            :meth:`flat`).
        lower : Surface, str, pathlib.Path, int or float, optional
            The floor, per cell (``numpy.maximum`` with it). ``None`` (the default)
            sets no floor.
        upper : Surface, str, pathlib.Path, int or float, optional
            The ceiling, per cell (``numpy.minimum`` with it), applied after
            ``lower`` -- so where the two cross, ``upper`` wins. ``None`` (the
            default) sets no ceiling.

        Returns
        -------
        Surface
            A ``kind="clamp"`` surface.

        Raises
        ------
        TypeError
            Neither ``lower`` nor ``upper`` given, or an operand of another type.
        """
        if lower is None and upper is None:
            raise TypeError(
                "Surface.clamp() needs lower= or upper=; with neither it would "
                "return `surface` unchanged. If you wrote `a.clamp(b)`: clamp is a "
                "classmethod, so `a` binds to `cls` and is dropped. Write "
                "`Surface.clamp(a, lower=..., upper=...)`, or the instance-bound "
                "`a.between(lower=..., upper=...)`."
            )
        return cls(
            kind="clamp",
            operands=(
                cls._coerce(surface),
                None if lower is None else cls._coerce(lower),
                None if upper is None else cls._coerce(upper),
            ),
        )

    @classmethod
    def where(cls, zone, inside, outside) -> Surface:
        """Use ``inside`` for cells whose centroid is in ``zone``, else ``outside``.

        ``zone`` is a shapely geometry or a path to a vector file. A centroid
        exactly ON the zone's edge counts as inside it (ledger 243).

        Parameters
        ----------
        zone : shapely geometry, str or pathlib.Path
            The area that selects ``inside``. A geometry is taken to be in the grid's
            CRS. A path is read when the surface is sampled, reprojected to the grid's
            CRS when both declare one, and all its features are unioned into one zone.
        inside : Surface, str, pathlib.Path, int or float
            The surface for cells whose centroid is in ``zone``: a :class:`Surface`, a raster path
            (``str`` or :class:`~pathlib.Path`, read as :meth:`raster`), or a python
            ``int``/``float`` (read as :meth:`flat`).
        outside : Surface, str, pathlib.Path, int or float
            The surface everywhere else, coerced the same way.

        Returns
        -------
        Surface
            A ``kind="where"`` surface.

        Raises
        ------
        TypeError
            ``inside`` or ``outside`` is of another type.
        """
        return cls(
            kind="where", zone=zone, operands=(cls._coerce(inside), cls._coerce(outside))
        )

    @classmethod
    def isopach(cls, source) -> Surface:
        """A layer bottom defined by a *thickness* map subtracted from above.

        Parameters
        ----------
        source : Surface, str, pathlib.Path, int or float
            The THICKNESS map, in model length units: a :class:`Surface`, a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`), or a python ``int``/``float`` (read as
            :meth:`flat`). The result is ``previous - thickness`` per cell, so it cannot be the
            model top (``ValueError`` when sampled).

        Returns
        -------
        Surface
            A ``kind="isopach"`` surface.
        """
        return cls(kind="isopach", operands=(cls._coerce(source),))

    @classmethod
    def shift(cls, source, distance: float) -> Surface:
        """``source`` displaced vertically by ``distance`` (negative = downward).

        Parameters
        ----------
        source : Surface, str, pathlib.Path, int or float
            The surface to move: a :class:`Surface`, a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`), or a python ``int``/``float`` (read as
            :meth:`flat`).
        distance : float
            The vertical displacement in model length units, added to every cell;
            negative moves it down.

        Returns
        -------
        Surface
            A ``kind="shift"`` surface.
        """
        return cls(kind="shift", value=float(distance), operands=(cls._coerce(source),))

    @classmethod
    def toward(cls, target, fraction: float) -> Surface:
        """A surface ``fraction`` of the way from the one above down to ``target``.

        ``previous + fraction * (target - previous)``, per cell. This is the one
        primitive vertical subdivision needs: ``fraction=1/3`` of the way from a
        unit's top to its base is a contact that no amount of ``shift``/envelope
        composition can express, because the interval is per-cell and unknown
        until the surfaces are sampled.

        It exists INSTEAD of general ``Surface - Surface`` and ``Surface *
        scalar``. Those are deliberately absent (see
        ``test_subtracting_a_surface_from_a_surface_is_unsupported``) and adding
        them would widen the elevation-vs-thickness ambiguity that already makes
        a thickness map easy to confuse with a contact. The recurrence form also
        dodges the fact that :meth:`values` carries only ONE ``previous`` slot:
        each cut resolves against the contact immediately above it, so a caller
        subdividing a unit must renormalise cumulative fractions into per-step
        ones. :meth:`myflopy.layers.LayerStack.add`'s ``split=`` does that for
        you, and is the spelling to reach for.

        ``target`` must be ABSOLUTE. A relative target (``offset_below`` /
        ``constant_thickness``) re-resolves against each sub-contact rather than
        the unit top, which yields plausible, wrong, monotonically drifting
        elevations -- ``toward(constant_thickness(60), f)`` at thirds gives
        ``[80, 50, -10]`` where the answer is ``[80, 60, 40]``.

        Parameters
        ----------
        target : Surface, str, pathlib.Path, int or float
            The ABSOLUTE surface to descend toward: a :class:`Surface`, a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`), or a python ``int``/``float`` (read as
            :meth:`flat`). A target containing a relative, ``isopach`` or ``toward``
            surface, at any depth of an algebra tree, is refused.
        fraction : float
            The share of the remaining interval to descend, in ``(0, 1]``: ``1.0``
            lands exactly on ``target``. Needs a surface above it, so it cannot be the
            model top (``ValueError`` when sampled).

        Returns
        -------
        Surface
            A ``kind="toward"`` surface.

        Raises
        ------
        TypeError
            If ``target`` is a relative surface, or ``fraction`` is outside
            ``(0, 1]``.
        """

        target = cls._coerce(target)
        if _depends_on_previous(target):
            raise TypeError(
                f"Surface.toward() needs a target that does not depend on the "
                f"surface above it; {target.kind!r} does. Such a target measures "
                f"from each cut rather than from the unit top, so the contacts "
                f"drift downward instead of dividing the unit. Give the unit's "
                f"actual base surface -- e.g. Raster('sand_base.tif') -- or, for "
                f"a unit declared by thickness, split it with "
                f"`.add(name, thickness=T, split=N)`, which needs no target."
            )
        fraction = float(fraction)
        if not 0.0 < fraction <= 1.0:
            raise TypeError(
                f"Surface.toward() fraction must be in (0, 1], got {fraction}. "
                f"It is the share of the REMAINING interval to descend, so 1.0 "
                f"lands exactly on the target and 0 would be a zero-thickness "
                f"layer."
            )
        return cls(kind="toward", value=fraction, operands=(target,))

    @classmethod
    def basin(
        cls,
        within,
        depth,
        *,
        water_surface,
        shape: str = "bowl",
        layer: str | None = None,
    ) -> Surface:
        """An assumed lakebed: a lake's water surface less a depth that grows from its shore.

        For a lake with no bathymetry survey. Each footprint part is given a profile
        that is 0 at its shoreline and the full ``depth`` at its pole of
        inaccessibility -- the centre of the largest circle inside it, at radius
        ``R``: at distance ``d`` inside the shore the lake is
        ``depth * f(min(d / R, 1))`` deep, ``f(s) = s * (2 - s)`` for a ``"bowl"`` and
        ``f(s) = s`` for a ``"cone"``. The bed is ``water_surface`` less that depth.
        On a circle the bowl is a paraboloid of revolution whose mean depth is
        ``depth / 2`` and the cone's is ``depth / 3``, so at one maximum depth the
        bowl holds 1.5x the water (``docs/reproducers/basin_volume.py``). The profile
        is a property of the polygon alone, so the same declaration gives the same bed
        on any mesh.

        A cell the footprint COVERS (the lake's own cell rule,
        :func:`~myflopy.modflow.mf6.grid.geometry.footprint_cells`) reads the AREA
        MEAN of the depth over the whole cell -- 0 where the cell is outside the
        lake -- so the lake's volume, the sum of cell area times depth, is the
        profile's on every mesh (within 0.06% measured, where a centroid sample
        missed it by 10% at 200 ft cells on a 300 ft lake); a consequence is that no
        cell reaches ``depth`` itself. Every other cell is ``NaN``: ``mf.lak`` reads a lake neighbour with
        no bed as the model top, as it does for a clipped bathymetry raster. So it is
        not a LAYER surface on its own: as a stack's top or a unit's bottom (or inside
        algebra that keeps the ``NaN``) sampling raises naming the cells -- compose
        ``mf.Min(<surface>, basin)`` to bake the bed into a surface that has a value
        everywhere (ledger 555).

        It is an ABSOLUTE surface: pass it as ``stack.excavate(bottom=...)`` (with the
        same footprint as ``within=``) and declare the lake with ``mf.lak(lakes=...)``,
        which takes its bed from the excavation's record, or pass it to
        ``mf.lak(lake_bottom=...)`` directly. It saves in a spec like a raster
        surface: a footprint FILE rides under ``zone`` (so a bundled project carries
        it and ``runs.drift()`` sees it edited), an inline footprint as WKT + CRS.

        Parameters
        ----------
        within : str, pathlib.Path, mf.Area, geopandas.GeoDataFrame, geopandas.GeoSeries or shapely geometry
            The lake outline(s), polygonal with area. A path is RESOLVED and read now
            (with ``layer``) -- a missing file or layer, a line, a point or a
            zero-area feature raises here -- and read again when the surface is
            sampled, reprojected to the grid's CRS; each FEATURE is one lake, and each
            polygon part of it is profiled on its own radius. An ``mf.Area`` is kept
            inline; a frame, series or geometry is coerced to one (one with no CRS is
            in the grid's frame). Cell numbers or a mask raise ``TypeError``: they
            name cells on one mesh.
        depth : float or str
            The MAXIMUM depth, in model length units, reached at each part's pole of
            inaccessibility: a positive number for every lake, or the name of a column
            of the footprint FILE giving each feature its own (a NULL, non-numeric,
            zero or negative value raises naming the feature). A column needs
            ``within`` to be a file.
        water_surface : float, str, pathlib.Path or Surface
            The elevation the depth is measured down from -- over a natural lake, the
            lidar ground surface IS the water surface. A number (model units); a
            :class:`Surface` (``mf.Raster(dem)``, algebra; declare ``units=`` on it
            when it is not in model units); a raster path (``str`` or
            :class:`~pathlib.Path`, read as :meth:`raster`); or the name of a column
            of the footprint file giving each lake one flat stage (a NULL or
            non-numeric value raises naming the feature). A ``str`` is a column
            first, then a raster file; neither raises ``ValueError``. Sampled with the
            stack's own method and ``length_units``. A surface measured from the one
            above it (``offset_below``, ``isopach``, ...) or carrying a
            ``fill="propagate"`` source anywhere raises ``ValueError``.
        shape : {'bowl', 'cone'}, default 'bowl'
            The profile: ``"bowl"`` (paraboloid, mean depth ``depth / 2`` on a
            circle) or ``"cone"`` (linear, mean depth ``depth / 3``). Anything else
            raises ``ValueError`` listing both.
        layer : str, optional
            The layer of the footprint FILE to read, by name; ``None`` (the default)
            reads the first. An integer raises ``TypeError`` (in a layer stack it would
            read as a MODEL layer), and ``layer`` with an inline footprint raises
            ``ValueError``.

        Returns
        -------
        Surface
            A ``kind="basin"`` surface.

        Raises
        ------
        ValueError
            An unknown ``shape``; a ``depth`` that is not a positive finite number or
            a column of the file; a column (``depth`` or ``water_surface``) or
            ``layer`` with an inline footprint; a NULL, non-numeric or non-positive
            column value, naming the feature; a relative or ``fill="propagate"``
            water surface; a non-polygonal or zero-area footprint. When the surface
            is SAMPLED: features covering no cell (all of them -- some only warn, by
            name, ledger 211), two features overlapping, a column no longer in the
            file, or two features giving one cell different water surfaces; and when a
            layer stack samples it as a top or unit bottom that is ``NaN`` off the lake.
        TypeError
            ``within``, ``depth`` or ``water_surface`` of another type (a ``bool``
            included), or an integer ``layer``.
        FileNotFoundError
            ``within`` or a ``water_surface`` path that does not exist.

        See Also
        --------
        myflopy.layers.LayerStack.excavate : Cut the stack to the bed.
        myflopy.modflow.mf6.grid.geometry.basin_depths : The profile itself.

        Examples
        --------
        >>> import numpy as np, shapely
        >>> import myflopy as mf
        >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
        >>> vor = irregular_voronoi_grid(nrow=10, ncol=10, cell_size=20.0)
        >>> x0, y0, x1, y1 = vor.gdf_vorPolys.total_bounds
        >>> lake = mf.Area(shapely.Point((x0 + x1) / 2, (y0 + y1) / 2).buffer(60.0, 64))
        >>> bed = mf.Basin(lake, 8.0, water_surface=100.0)
        >>> values = bed.values(vor)
        >>> bool(np.isnan(values).any()), bool((values[np.isfinite(values)] < 100.0).all())
        (True, True)
        >>> bed.to_dict()["shape"]
        'bowl'
        """

        if layer is not None:
            layer = _vector_layer_name(
                layer, of="within", model_layer="not chosen by a basin, which is one surface"
            )
        zone, area, frame = _coerce_footprint(
            within, layer=layer, verb="Surface.basin", noun="basin"
        )
        if layer is not None and frame is None:
            raise ValueError(
                f"Surface.basin(layer={layer!r}) names a layer of the footprint FILE, and "
                "within= is not one -- an mf.Area or a geometry has no layers. Pass the path."
            )

        def column_of(value) -> str | None:
            """``value`` when it names a column of the footprint file, else None."""

            if frame is not None and value in frame.columns and value not in ("geometry", "_source_name"):
                return value
            return None

        depth_value = depth_column = None
        if isinstance(depth, (bool, np.bool_)):
            raise TypeError("Surface.basin(depth=) is a length, not a bool.")
        if isinstance(depth, numbers.Real):
            depth_value = float(depth)
        elif isinstance(depth, str):
            depth_column = column_of(str(depth))
            if depth_column is None:
                raise ValueError(_basin_column_refusal("depth", depth, zone, frame))
            for label, value in zip(frame.index, frame[depth_column], strict=True):
                number = _basin_number(value, label=str(label), column=depth_column, what="depth")
                if not number > 0.0:
                    raise ValueError(
                        f"basin: feature {label!r} has {depth_column}={value!r}; a depth "
                        "must be positive."
                    )
        else:
            raise TypeError(
                f"Surface.basin(depth=) got a {type(depth).__name__}; pass a positive "
                "number (model length units) or the name of a column of the footprint file."
            )

        operands: tuple = ()
        surface_column = None
        if isinstance(water_surface, Surface):
            operands = (water_surface,)
        elif isinstance(water_surface, (bool, np.bool_)):
            raise TypeError("Surface.basin(water_surface=) is an elevation, not a bool.")
        elif isinstance(water_surface, numbers.Real):
            operands = (Surface.flat(float(water_surface)),)
        elif isinstance(water_surface, str):
            surface_column = column_of(str(water_surface))
            if surface_column is not None:
                for label, value in zip(frame.index, frame[surface_column], strict=True):
                    _basin_number(value, label=str(label), column=surface_column, what="water surface")
            elif Path(water_surface).is_file():
                operands = (Surface.raster(Path(water_surface)),)
            else:
                raise ValueError(
                    _basin_column_refusal("water_surface", water_surface, zone, frame)
                    + " Nor is it an existing raster file."
                )
        elif isinstance(water_surface, Path):
            if not water_surface.is_file():
                raise FileNotFoundError(
                    f"Surface.basin(water_surface=) names a raster that does not exist: {water_surface}"
                )
            operands = (Surface.raster(water_surface),)
        else:
            raise TypeError(
                f"Surface.basin(water_surface=) got a {type(water_surface).__name__}; pass "
                "a number, a Surface (mf.Raster, ...), a raster path, or the name of a "
                "column of the footprint file."
            )
        return cls(
            kind="basin",
            value=depth_value,
            depth_column=depth_column,
            operands=operands,
            surface_column=surface_column,
            zone=zone if zone is not None else area,
            zone_layer=layer,
            shape=shape,
        )

    # -- fluent surface algebra ------------------------------------------- #
    # Read base-surface-first, left to right: ``Raster(bedrock).capped_at(ground - 5)``.
    # These are thin wrappers over the constructors above (same engine), so
    # ``Min``/``Max``/``Clamp``/``Where`` keep working unchanged.
    def capped_at(self, other) -> Surface:
        """This surface, but never *above* ``other`` (per cell) -- a ceiling.

        The fluent form of :meth:`minimum` (lower envelope); e.g. a layer top
        ``capped_at`` an eroding valley floor.

        Parameters
        ----------
        other : Surface, str, pathlib.Path, int or float
            The ceiling: a :class:`Surface`, a raster path (read as :meth:`raster`),
            or a python ``int``/``float`` (read as :meth:`flat`). Another type raises
            ``TypeError``.

        Returns
        -------
        Surface
            ``Surface.minimum(self, other)``.
        """
        return Surface.minimum(self, other)

    def floored_at(self, other) -> Surface:
        """This surface, but never *below* ``other`` (per cell) -- a floor.

        The fluent form of :meth:`maximum` (upper envelope).

        Parameters
        ----------
        other : Surface, str, pathlib.Path, int or float
            The floor: a :class:`Surface`, a raster path (read as :meth:`raster`), or
            a python ``int``/``float`` (read as :meth:`flat`). Another type raises
            ``TypeError``.

        Returns
        -------
        Surface
            ``Surface.maximum(self, other)``.
        """
        return Surface.maximum(self, other)

    def between(self, *, lower=None, upper=None) -> Surface:
        """This surface constrained to ``[lower, upper]`` per cell (fluent :meth:`clamp`).

        Parameters
        ----------
        lower : Surface, str, pathlib.Path, int or float, optional
            The floor, coerced as in :meth:`clamp`; may be relative, e.g.
            ``Surface.offset_below(5)``. ``None`` (the default) sets none.
        upper : Surface, str, pathlib.Path, int or float, optional
            The ceiling, applied after ``lower`` (it wins where they cross). ``None``
            (the default) sets none. Giving neither raises ``TypeError``.

        Returns
        -------
        Surface
            ``Surface.clamp(self, lower=lower, upper=upper)``.
        """
        return Surface.clamp(self, lower=lower, upper=upper)

    def below(self, distance: float) -> Surface:
        """A surface ``distance`` units *below* this one (same as ``self - distance``).

        Parameters
        ----------
        distance : float
            How far down, in model length units; a negative value moves it up.

        Returns
        -------
        Surface
            A ``kind="shift"`` surface over this one.
        """
        return Surface.shift(self, -float(distance))

    def above(self, distance: float) -> Surface:
        """A surface ``distance`` units *above* this one (same as ``self + distance``).

        Parameters
        ----------
        distance : float
            How far up, in model length units; a negative value moves it down.

        Returns
        -------
        Surface
            A ``kind="shift"`` surface over this one.
        """
        return Surface.shift(self, float(distance))

    def within(self, zone, *, outside) -> Surface:
        """Use this surface for cells inside ``zone``, ``outside`` elsewhere.

        The fluent form of :meth:`where`; e.g.
        ``Flat(15).within(channel, outside=Flat(2))``.

        Parameters
        ----------
        zone : shapely geometry, str or pathlib.Path
            The area where this surface applies, judged by cell centroid (an edge
            counts as inside); a geometry is taken to be in the grid's CRS, a vector
            file is reprojected to it (when both declare one) and its features
            unioned. See :meth:`where`.
        outside : Surface, str, pathlib.Path, int or float
            The surface for every other cell: a :class:`Surface`, a raster path, or a
            python ``int``/``float``. Required, keyword-only.

        Returns
        -------
        Surface
            ``Surface.where(zone, self, outside)``.
        """
        return Surface.where(zone, self, outside)

    def thickness(self) -> Surface:
        """A layer bottom set by this *thickness* map below the surface above
        (fluent :meth:`isopach`)."""
        return Surface.isopach(self)

    def __sub__(self, distance) -> Surface:
        """``surface - distance``: a new surface shifted down by a constant elevation."""

        if not isinstance(distance, (int, float, np.integer, np.floating)):
            return NotImplemented
        return Surface.shift(self, -float(distance))

    def __add__(self, distance) -> Surface:
        """``surface + distance``: a new surface shifted up by a constant elevation."""

        if not isinstance(distance, (int, float, np.integer, np.floating)):
            return NotImplemented
        return Surface.shift(self, float(distance))

    __radd__ = __add__

    @property
    def is_relative(self) -> bool:
        """True if this surface is defined relative to the surface above it."""

        return self.kind in _RELATIVE_KINDS

    @property
    def is_derived(self) -> bool:
        """True if this surface resolves to a cached, regenerable raster."""

        return self.kind == "contours"

    def _derived_raster(self):
        """Return the :class:`DerivedRaster` cache handle for a derived surface."""

        from myflopy.modflow.utils.derived_raster import DerivedRaster

        def produce() -> None:
            """Regenerate the derived raster by interpolating this surface's contours."""

            from myflopy.modflow.utils.contour_interp import (
                interpolate_contours_to_raster,
            )

            interpolate_contours_to_raster(
                self.contours,
                out=self.out,
                z_field=self.z,
                region_raster=self.region_raster,
                region_vector=self.region_vector,
                clip=self.clip,
                resolution=self.resolution,
                epsg=self.epsg,
                progress=self.progress,
                **self.grass_kwargs,
            )

        sources = [self.contours]
        # BOTH region kinds belong here, not just the raster. `sources` is
        # content-tracked while `params` holds only values, so while
        # `region_vector` rode `grass_kwargs` its path was recorded and its
        # CONTENT was not -- editing the domain polygon left the cache "fresh"
        # and returned a raster interpolated over the old extent.
        for region in (self.region_raster, self.region_vector):
            if region is not None:
                sources.append(region)
        params = {
            "z": self.z,
            "resolution": self.resolution,
            "epsg": self.epsg,
            # Geometry-affecting, so it belongs here. `progress` deliberately
            # does NOT: it changes what you see, not what gets written.
            "clip": self.clip,
            "grass_kwargs": {k: str(v) for k, v in self.grass_kwargs.items()},
        }
        return DerivedRaster(self.out, sources, params, produce)

    def cache_status(self) -> str | None:
        """``"missing"``/``"fresh"``/``"stale"`` for derived surfaces, else ``None``."""

        return self._derived_raster().status() if self.is_derived else None

    def resolve_source(self, *, refresh: bool = False, warn: bool = True) -> Path | float:
        """Return a raster path or flat constant for centroid sampling.

        Derived surfaces (contours) interpolate to a cached raster on first use;
        a stale cache is reused with a warning unless ``refresh=True``.

        Parameters
        ----------
        refresh : bool, default False
            For a contour surface, re-run the GRASS interpolation even when the cached
            raster is fresh. A missing cache is always built. Ignored by the other
            kinds.
        warn : bool, default True
            For a contour surface, warn when a STALE cache (its contours, region or
            parameters changed) is reused without ``refresh``. ``False`` reuses it
            silently. Ignored by the other kinds.

        Returns
        -------
        pathlib.Path or float
            ``value`` for a ``"flat"`` surface, ``path`` for a ``"raster"`` one, the
            cached GeoTIFF's path for a ``"contours"`` one.

        Raises
        ------
        TypeError
            Any other kind -- relative, points, array or algebra surfaces have no
            single raster or constant.
        """

        if self.kind == "flat":
            return self.value
        if self.kind == "raster":
            return self.path
        if self.kind == "contours":
            return self._derived_raster().ensure(refresh=refresh, warn=warn)
        raise TypeError(f"Surface kind '{self.kind}' has no raster/constant source.")

    def _basin_label(self) -> str:
        """What a basin is called in messages: its file's stem, its Area's name, or ``"basin"``."""

        if isinstance(self.zone, Area):
            return self.zone.name or "basin"
        return Path(self.zone).stem

    def _basin_values(self, vor, sample) -> np.ndarray:
        """Sample a ``kind="basin"`` surface: water surface less the cell-mean depth.

        ``sample`` samples an operand with the caller's arguments (:meth:`values`'
        ``_op``). See :meth:`basin` for the rule; NaN off the footprint's cells.
        """

        polys = vor.gdf_vorPolys
        grid_crs = getattr(vor, "crs", None)
        name = self._basin_label()
        frame = _footprint_frame(
            None if isinstance(self.zone, Area) else Path(self.zone),
            self.zone if isinstance(self.zone, Area) else None,
            layer=self.zone_layer, crs=grid_crs, name=name, noun="basin",
        )
        for column in (self.depth_column, self.surface_column):
            if column is not None and column not in frame:
                raise ValueError(
                    f"basin {name!r}: {column!r} is no longer a column of its footprint file."
                )
        per_feature = footprint_cells(polys, frame.geometry)
        refuse_footprint_misses(
            polys, frame, per_feature,
            noun=f"features of basin '{name}'",
            consequence=", so they carry no bed",
            partial="warn",
            source_crs=frame.attrs.get("source_crs"),
            grid_crs=grid_crs,
        )
        _refuse_overlapping_basins(frame, name=name)
        ncpl = len(polys)
        depth_mean = np.zeros(ncpl, dtype=float)
        water = np.full(ncpl, np.nan, dtype=float)
        claimed: dict[int, float] = {}
        conflicts = []
        for label, geometry, found, row in zip(
            frame.index, frame.geometry, per_feature, range(len(frame)), strict=True
        ):
            if not found:
                continue
            if self.depth_column is None:
                depth = float(self.value)
            else:
                value = frame[self.depth_column].iloc[row]
                depth = _basin_number(value, label=str(label), column=self.depth_column, what="depth")
                if not depth > 0.0:
                    raise ValueError(
                        f"basin {name!r}: feature {label!r} has {self.depth_column}={value!r}; "
                        "a depth must be positive."
                    )
            parts = [(part, depth) for part in shp.get_parts(geometry) if part.area > 0.0]
            depth_mean[found] += _basin_cell_depths(
                polys.loc[found], parts, shape=self.shape, crs=grid_crs, label=f"{name}/{label}"
            )
            if self.surface_column is not None:
                value = frame[self.surface_column].iloc[row]
                stage = _basin_number(value, label=str(label), column=self.surface_column, what="water surface")
                for cell in found:
                    if cell in claimed and claimed[cell] != stage:
                        conflicts.append((cell, claimed[cell], stage))
                    claimed[cell] = stage
                water[found] = stage
        if conflicts:
            shown = ", ".join(f"cell {c} ({a:g} vs {b:g})" for c, a, b in conflicts[:10])
            raise ValueError(
                f"basin {name!r}: two features give one cell different "
                f"{self.surface_column} water surfaces: {shown}. A cell has one bed; move "
                "the lakes apart or give them one water surface."
            )
        cells = sorted({cell for found in per_feature for cell in found})
        if self.surface_column is None:
            water = sample(self.operands[0])
        bed = np.full(ncpl, np.nan, dtype=float)
        bed[cells] = water[cells] - depth_mean[cells]
        return bed

    def values(
        self,
        vor,
        previous: np.ndarray | None = None,
        *,
        method: str = "area",
        refresh: bool = False,
        length_units: str | None = None,
    ) -> np.ndarray:
        """Return per-cell elevation values sampled on ``vor``.

        ``previous`` is the surface immediately above (or ``None`` for the model
        top). Relative surfaces (:meth:`offset_below`, :meth:`constant_thickness`)
        and raster ``fill="propagate"`` require it. ``method`` controls raster
        sampling: ``"area"`` (area-weighted, default) or ``"centroid"``.
        ``refresh`` rebuilds a derived (contour) surface's cached raster.
        ``length_units`` converts every surface in the tree that declares ``units=``
        -- an algebra surface's operands included -- before anything combines.

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto: its cell polygons (``gdf_vorPolys``), centroids
            and CRS are read.
        previous : numpy.ndarray of float, shape (ncpl,), optional
            The surface immediately above, already sampled on ``vor``. Required by the
            kinds that consult it -- ``offset_below``, ``constant_thickness``,
            ``isopach``, ``toward``, and a raster or contour surface with
            ``fill="propagate"`` that has gaps -- which raise ``ValueError`` without
            it. ``None`` (the default) means this is the model top. Algebra kinds pass
            it to every operand.
        method : {'area', 'centroid'}, default 'area'
            How a raster or contour surface is sampled: the area-weighted mean of the
            pixels in each cell, or the single pixel under its centroid. Another value
            raises ``ValueError`` there; the other kinds ignore it.
        refresh : bool, default False
            Rebuild a contour surface's cached raster before sampling (see
            :meth:`resolve_source`). Passed down to every operand.
        length_units : str, optional
            The model's length unit (any spelling :class:`LayerStack` takes:
            ``"feet"``, ``"ft"``, ``"meters"``, ``"m"``, ...). Each surface of the
            tree that declares ``units=`` -- a ``flat``, ``raster``, ``contours``,
            ``points`` or ``array`` one, at any depth of an algebra surface -- is
            converted to it by one factor BEFORE it is combined, compared or used to
            fill: ``Surface.maximum(Surface.flat(30, units="meters"), 0.0)`` is
            ``98.43`` in feet, and a ``fill="propagate"`` raster in metres takes the
            surface above as it is (already in model units). A relative surface
            (``offset_below``/``constant_thickness``) or an algebra surface carrying
            ``units`` itself raises ``ValueError`` -- the first is a distance in model
            units, and on the second write ``units=`` on the operands. ``None`` (the
            default) or ``"unknown"`` declares no model unit: a tree with no
            ``units=`` anywhere is sampled as it is, and a surface declaring
            ``units=`` raises ``ValueError`` naming ``length_units=`` (ledger
            433). Passed down to every operand.

        Returns
        -------
        numpy.ndarray of float, shape (ncpl,)
            One elevation per cell, in ``length_units`` when given, else as the
            surfaces hold them (none of which declares ``units=``). ``NaN`` where a raster has no data and nothing
            filled it, and where a ``"points"`` interpolation has none.

        Raises
        ------
        ValueError
            A relative kind with no ``previous``; a ``from_array`` surface of the
            wrong size; an unknown ``fill``; an unknown ``method`` on a raster; with
            ``units=`` on a relative or algebra surface, an unknown length unit on
            either side, or ``units=`` on a surface sampled with no
            ``length_units``.
        TypeError
            An unknown ``kind``.
        """

        # The conversion is per LEAF, before anything combines (ledger 405). The
        # callers converted only the surface they were handed, so a metre DEM inside
        # `capped_at`/`shift`/`where`/... reached the stack as feet -- 3.28x off,
        # silently -- and a `fill="propagate"` raster in metres multiplied the
        # surface it inherited, already in feet, by the same factor.
        factor = 1.0
        if self.units is not None:
            if self.kind in _RELATIVE_KINDS:
                raise ValueError(
                    "units= is only supported on absolute surfaces, not relative "
                    "(offset_below / constant_thickness): their value is a distance "
                    "in model length units."
                )
            if self.kind not in _UNIT_KINDS:
                raise ValueError(
                    f"units={self.units!r} on a {self.kind!r} surface: an algebra "
                    "surface takes its units from its operands. Declare units= on "
                    "each operand that is not in model units -- each is converted "
                    "before they combine."
                )
            if _declares_no_unit(length_units):
                # With no model unit a declared `units=` cannot be converted, and it was
                # read as model units in silence: a metre DEM sampled by a stack declared
                # `length_units=None` came back 3.28x off (ledger 433).
                raise ValueError(
                    f"A {self.kind!r} surface declares units={self.units!r}, but it was "
                    f"sampled with length_units={length_units!r} -- no model length unit "
                    "to convert it into. Declare the model's unit: length_units= on the "
                    "LayerStack, LayerSurfaces.sample/to_disv or sample_bed call, or "
                    "ModelContext(length_units=) for mf.lak's lake_bottom. Or drop "
                    "units= if the surface is already in model units."
                )
            factor = _length_factor(self.units, length_units)

        def _op(operand):
            """Sample one operand surface on ``vor`` as a float array (same args as the parent)."""

            return np.asarray(
                operand.values(
                    vor, previous, method=method, refresh=refresh,
                    length_units=length_units,
                ),
                dtype=float,
            )

        if self.kind in ("min", "max"):
            stacked = np.vstack([_op(op) for op in self.operands])
            return np.nanmin(stacked, axis=0) if self.kind == "min" else np.nanmax(stacked, axis=0)

        if self.kind == "shift":
            return _op(self.operands[0]) + float(self.value)

        if self.kind == "clamp":
            base, lower, upper = self.operands
            vals = _op(base)
            if lower is not None:
                vals = np.maximum(vals, _op(lower))
            if upper is not None:
                vals = np.minimum(vals, _op(upper))
            return vals

        if self.kind == "where":
            inside, outside = self.operands
            return np.where(_cells_in_zone(vor, self.zone), _op(inside), _op(outside))

        if self.kind == "basin":
            return self._basin_values(vor, _op)

        if self.kind == "isopach":
            if previous is None:
                raise ValueError("Surface.isopach cannot be the model top.")
            return previous - _op(self.operands[0])

        if self.kind == "toward":
            # Handled here rather than via `is_relative`, matching `isopach`
            # above: both NEED `previous`, but neither is a fixed offset from
            # it, so `_RELATIVE_KINDS` (which `values` reads as
            # `previous - value`) would be exactly wrong for them.
            if previous is None:
                raise ValueError(
                    "Surface.toward cannot be the model top -- there is no "
                    "surface above to descend from."
                )
            return previous + float(self.value) * (_op(self.operands[0]) - previous)

        if self.is_relative:
            if previous is None:
                raise ValueError(
                    f"Surface.{self.kind} cannot be the model top -- there is no "
                    "surface above to offset from."
                )
            return previous - float(self.value)

        if self.kind == "flat":
            return np.full(len(np.asarray(vor.centroids[0])), float(self.value) * factor)

        if self.kind == "array":
            vals = np.asarray(self.array, dtype=float).ravel()
            n = len(np.asarray(vor.centroids[0]))
            if vals.size != n:
                raise ValueError(
                    f"from_array surface has {vals.size} values but the grid has {n} "
                    "cells; use Surface.from_points to interpolate onto another grid."
                )
            return vals * factor

        if self.kind == "points":
            from scipy.interpolate import griddata

            cx = np.asarray(vor.centroids[0])
            cy = np.asarray(vor.centroids[1])
            xs, ys, zs = self.points
            return griddata((xs, ys), zs, (cx, cy), method=self.method) * factor

        gdf = get_raster_vals_at_centroids(
            vor, [self.resolve_source(refresh=refresh)], [0], method=method,
            nodata=self.nodata,
            # A raster covering part of the model is normal here: `fill="propagate"`
            # inherits the rest from the surface above, and `.qc()` counts what is
            # left. Warning per cell would be noise; a NaN that would reach a DISV is
            # refused there instead (`refuse_unbounded_cells`, ledger 587).
            gaps="ignore",
        )
        # Converted BEFORE the fill: `previous` is in model units already.
        vals = gdf[0].to_numpy(dtype=float) * factor
        if self.fill == "propagate":
            mask = np.isnan(vals)
            if mask.any():
                if previous is None:
                    raise ValueError(
                        "fill='propagate' cannot be used on the model top -- there "
                        "is no surface above to inherit from."
                    )
                vals = vals.copy()
                vals[mask] = previous[mask]
        elif self.fill is not None:
            raise ValueError(
                f"Unknown fill mode {self.fill!r}; expected None or 'propagate'."
            )
        return vals


def sample_bed(
    surface, grid, ground, *, length_units: str | None = None, refresh: bool = False
) -> np.ndarray:
    """Sample a lakebed or excavation floor onto every cell of ``grid``.

    The one bed sampler behind ``mf.lak(lake_bottom=)`` and
    ``LayerStack.excavate(bottom=)``, so a stack cut to a bed and a lake declared on
    it read the SAME float at every cell (the excavation design, section 3.5) --
    which is what lets ``mf.lak`` compare its bed with an excavation's record
    exactly. It is always AREA-weighted, whatever a stack's ``build(method=)`` says:
    that is how the lake samples, and a centroid sample of the same raster was
    measured up to 1.027 ft away from it.

    Parameters
    ----------
    surface : Surface, str or pathlib.Path
        The bed: any :class:`Surface` (``mf.Raster``, ``mf.Contours``, ``mf.Flat``,
        algebra), or a raster file path, read as ``Surface.raster(path)``.
    grid : VoronoiGridPlus
        The grid to sample onto.
    ground : array_like of float, shape (ncpl,)
        The model top -- the surface ABOVE the bed. A surface that consults the one
        above it (``offset_below``, ``fill="propagate"``) is resolved against it; an
        absolute one ignores it.
    length_units : str, optional
        The model's length unit. Every surface declaring ``units=`` is converted to
        it (``"meters"`` -> ``"feet"`` multiplies by 3.2808...), an operand of an
        algebra surface included -- ``mf.Raster(dem_m, units="meters").below(3)`` is
        the DEM in feet less 3 ft (see :meth:`Surface.values`). ``None`` (or
        ``"unknown"``) declares no unit: a surface declaring ``units=`` then raises
        ``ValueError`` naming ``length_units=`` (ledger 433), and one declaring
        none is sampled as it is.
    refresh : bool, default False
        Rebuild a derived (contour) surface's cached raster before sampling, as
        ``LayerStack.build(refresh=True)`` does for every other surface of the stack.
        A stale cache is otherwise reused with a warning.

    Returns
    -------
    numpy.ndarray of float, shape (ncpl,)
        One elevation per cell; ``NaN`` where the source has no data.

    Raises
    ------
    ValueError
        ``units=`` on a relative surface (``offset_below`` / ``constant_thickness``)
        or on an algebra surface itself, which a unit factor cannot convert; an
        unknown length unit; ``units=`` with no ``length_units``; any refusal :meth:`Surface.values` makes (a
        ``from_array`` of the wrong length, ...).

    Examples
    --------
    >>> import numpy as np
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> from myflopy.surfaces import Surface, sample_bed
    >>> vor = irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)
    >>> sample_bed(Surface.flat(30.0, units="meters"), vor, np.full(vor.ncpl, 200.0),
    ...            length_units="feet").round(2).tolist()
    [98.43, 98.43, 98.43, 98.43]
    >>> sample_bed(Surface.flat(30.0, units="meters").floored_at(0.0), vor,
    ...            np.full(vor.ncpl, 200.0), length_units="feet").round(2).tolist()
    [98.43, 98.43, 98.43, 98.43]
    >>> sample_bed(Surface.offset_below(5.0), vor, np.full(vor.ncpl, 200.0)).tolist()
    [195.0, 195.0, 195.0, 195.0]
    """

    if isinstance(surface, (str, Path)):
        surface = Surface.raster(surface)
    previous = np.asarray(ground, dtype=float)
    return np.asarray(
        surface.values(grid, previous=previous, refresh=refresh, length_units=length_units),
        dtype=float,
    )


@dataclass(frozen=True, eq=False)
class ExcavationRecord:
    """One excavation of a layer stack, resolved on one grid: which cells, and the bed at each.

    What ``LayerStack.excavate`` / ``LayerBuildResult.excavate`` leave behind on the
    built result (``result.excavations[name]``) and publish on the grid beside
    ``gdf_topbtm`` (``grid.excavations``) when the result is attached. It is the
    DECLARATION resolved on this grid and nothing derived from the cut: the connected
    layer, the absorbed slivers, the connected cell's thickness and the removed volume
    are computed on demand by ``LayerBuildResult.excavation_report()``, so nothing here
    can go stale through ``prune_isolated`` or a copy of the result.

    Parameters
    ----------
    name : str
        The excavation's key (converted with ``str``): the label in reports and
        refusals, and the key of ``result.excavations``.
    cells : sequence of int
        The footprint's cells, by the lake's area rule
        (:func:`~myflopy.modflow.mf6.grid.geometry.footprint_cells`) -- the same
        cells ``mf.lak`` resolves from the same polygons. Stored sorted, as a tuple
        of python ``int`` (a ``numpy.int64`` is converted), duplicates refused.
    bed : array_like of float, shape (len(cells),)
        The floor at each cell, in the stack's length units, as SAMPLED --
        area-weighted, before any capping at the ground. ``NaN`` where the source has
        no data (that cell was not cut). A cell whose bed is at or above the ground
        (a rim cell) keeps its sampled value here; it was not cut either. Stored as
        a read-only ``float64`` array, in ``cells`` order.
    min_remnant : float, default 0.1
        The declaration's thinnest cell allowed directly under the bed (length
        units), kept so the report can say which connected cells were left thinner
        than it because nothing active lay below them.

    Raises
    ------
    ValueError
        ``bed`` of another length than ``cells``, a repeated cell, a negative cell,
        or a negative or non-finite ``min_remnant``.

    See Also
    --------
    myflopy.layers.LayerStack.excavate : Declares one.
    myflopy.layers.LayerBuildResult.excavation_report : What the cut did, per excavation.

    Examples
    --------
    >>> import numpy as np
    >>> from myflopy.surfaces import ExcavationRecord
    >>> record = ExcavationRecord("pond", [5, 4], [95.0, 96.5])
    >>> record.cells, record.bed.tolist()
    ((4, 5), [96.5, 95.0])
    >>> record.bed_at(5)
    95.0
    """

    name: str
    cells: tuple[int, ...]
    bed: np.ndarray
    min_remnant: float = 0.1

    def __post_init__(self) -> None:
        """Normalize to plain types: a sorted tuple of ``int`` and a read-only float bed."""

        cells = [int(cell) for cell in self.cells]
        bed = np.asarray(self.bed, dtype=np.float64).ravel()
        if bed.size != len(cells):
            raise ValueError(
                f"excavation {self.name!r}: {bed.size} bed values for {len(cells)} cells."
            )
        if len(set(cells)) != len(cells):
            raise ValueError(f"excavation {self.name!r}: a cell is listed twice.")
        if any(cell < 0 for cell in cells):
            raise ValueError(f"excavation {self.name!r}: a cell number is negative.")
        order = np.argsort(cells, kind="stable")
        bed = bed[order].copy()
        bed.setflags(write=False)
        min_remnant = float(self.min_remnant)
        if not np.isfinite(min_remnant) or min_remnant < 0:
            raise ValueError(
                f"excavation {self.name!r}: min_remnant must be a finite number >= 0, "
                f"not {self.min_remnant!r}."
            )
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(self, "cells", tuple(cells[i] for i in order))
        object.__setattr__(self, "bed", bed)
        object.__setattr__(self, "min_remnant", min_remnant)

    def __eq__(self, other) -> bool:
        """Equal when the name, cells, ``min_remnant`` and every bed value agree (NaN == NaN)."""

        if not isinstance(other, ExcavationRecord):
            return NotImplemented
        return (
            self.name == other.name
            and self.cells == other.cells
            and self.min_remnant == other.min_remnant
            and np.array_equal(self.bed, other.bed, equal_nan=True)
        )

    __hash__ = None  # a record holds an array; equality is by value, so it is unhashable

    def __copy__(self) -> ExcavationRecord:
        """The record itself: it is immutable (frozen, a tuple of cells, a read-only bed)."""

        return self

    def __deepcopy__(self, memo) -> ExcavationRecord:
        """The record itself, for the same reason -- and pandas deep-copies a frame's
        ``attrs`` into EVERY frame derived from it. A lazy DISV publishes the records in
        its elevation table's ``attrs`` (ledger 323), where copying a 612-cell record
        cost ~35 us per derived frame (329 against 225 us for three column reads on a
        9,405-cell table); sharing it costs ~3 us."""

        return self

    def bed_at(self, cell: int) -> float:
        """The sampled bed at one of the record's cells.

        Parameters
        ----------
        cell : int
            A cell of :attr:`cells` (any integer type).

        Returns
        -------
        float
            The bed there, ``nan`` where the source had no data.

        Raises
        ------
        KeyError
            A cell the excavation does not cover.

        Examples
        --------
        >>> ExcavationRecord("pond", [3], [90.0]).bed_at(3)
        90.0
        """

        try:
            return float(self.bed[self.cells.index(int(cell))])
        except ValueError:
            raise KeyError(f"cell {cell} is not in excavation {self.name!r}.") from None

    def __repr__(self) -> str:
        finite = self.bed[np.isfinite(self.bed)]
        span = f"{finite.min():.6g}..{finite.max():.6g}" if finite.size else "no data"
        return (
            f"ExcavationRecord({self.name!r}, {len(self.cells)} cells, bed {span}, "
            f"min_remnant={self.min_remnant:g})"
        )


def _feature_floor(value, *, name: str, label: str, column: str) -> float:
    """A column floor as a finite float, or a refusal naming the feature (ledger 247)."""

    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(
            f"excavation {name!r}: feature {label!r} has {column}={value!r}, which is "
            "not a number."
        ) from None
    if not np.isfinite(number):
        raise ValueError(
            f"excavation {name!r}: feature {label!r} has no {column} ({value!r}) -- a "
            "NULL floor is not a floor; fill it or drop the feature."
        )
    return number


def _resolve_excavation(
    *,
    name: str,
    bottom,
    column: str | None,
    within,
    layer,
    min_remnant: float,
    grid,
    ground,
    length_units: str | None,
    refresh: bool = False,
) -> ExcavationRecord:
    """Resolve one declared excavation onto ``grid``: its cells and the bed at each.

    ``within`` is a vector path (read at ``layer``) or an :class:`~myflopy.areas.Area`
    (one feature per PART). Cells come from the shared footprint rule and the bed from
    the shared sampler, so a lake declared on the same polygons and floor resolves the
    same cells and the same float (design sections 3.4-3.5). A feature missing the grid
    warns by name and a total miss raises (ledger 211); a column floor is read per
    feature, and two features giving one cell different floors raise. Cells where the
    bed has no data are kept (NaN) and warned about; all of them raises.
    """

    polys = grid.gdf_vorPolys
    grid_crs = getattr(grid, "crs", None)
    footprints = _footprint_frame(
        None if isinstance(within, Area) else within,
        within if isinstance(within, Area) else None,
        layer=layer, crs=grid_crs, name=name, noun="excavation",
    )
    per_feature = footprint_cells(polys, footprints.geometry)
    refuse_footprint_misses(
        polys,
        footprints,
        per_feature,
        noun=f"features of excavation '{name}'",
        consequence=", so nothing is dug for them",
        partial="warn",
        source_crs=footprints.attrs.get("source_crs"),
        grid_crs=grid_crs,
    )
    cells = sorted({cell for found in per_feature for cell in found})
    if column is not None:
        if column not in footprints:
            raise ValueError(
                f"excavation {name!r}: bottom={column!r} is no longer a column of its "
                "footprint file."
            )
        floors: dict[int, float] = {}
        conflicts = []
        for label, found, value in zip(footprints.index, per_feature, footprints[column], strict=True):
            floor = _feature_floor(value, name=name, label=str(label), column=column)
            for cell in found:
                if cell in floors and floors[cell] != floor:
                    conflicts.append((cell, floors[cell], floor))
                floors[cell] = floor
        if conflicts:
            shown = ", ".join(f"cell {c} ({a:g} vs {b:g})" for c, a, b in conflicts[:10])
            raise ValueError(
                f"excavation {name!r}: two features of its footprint give one cell "
                f"different {column} floors: {shown}. A cell has one floor; move the "
                "features apart or give them one value."
            )
        bed = np.array([floors[cell] for cell in cells], dtype=float)
    else:
        sampled = sample_bed(bottom, grid, ground, length_units=length_units, refresh=refresh)
        bed = np.asarray(sampled, dtype=float)[cells]
    missing = [cell for cell, value in zip(cells, bed, strict=True) if not np.isfinite(value)]
    if cells and len(missing) == len(cells):
        raise ValueError(
            f"excavation {name!r}: its bottom has no data at any of its {len(cells)} "
            "footprint cells -- does the bed raster cover the footprint?"
        )
    if missing:
        message = (
            f"excavation {name!r}: the bottom has no data at {len(missing)} of its "
            f"{len(cells)} cells {missing[:10]}; they are NOT cut (mf.lak falls back to "
            "the same model top there). Usually cells the footprint only clips at its edge."
        )
        _logger.warning(message)
        warnings.warn(message, UserWarning, stacklevel=4)
    return ExcavationRecord(name, tuple(cells), bed, min_remnant)


def _refuse_overlapping_excavations(records) -> None:
    """Refuse two excavations sharing a cell, naming the cells and both names (F2).

    A shared cell would carry two floors and, once lakes are declared on them,
    connections from two lakes. Measured with 17 ft cells: ponds 2 or 10 ft apart
    shared 2 cells, 20 ft or more shared none (the design's ``q5_overlap.py``).
    """

    owner: dict[int, str] = {}
    shared: dict[tuple[str, str], list[int]] = {}
    for record in records:
        for cell in record.cells:
            first = owner.setdefault(cell, record.name)
            if first != record.name:
                shared.setdefault((first, record.name), []).append(cell)
    if shared:
        detail = "; ".join(
            f"{a!r} and {b!r} share {len(found)} cell(s) {sorted(found)[:10]}"
            for (a, b), found in shared.items()
        )
        raise ValueError(
            f"Excavations may not overlap: {detail}. A cell has one floor and, with a "
            "lake on each, would take connections from two lakes. Move the footprints "
            "apart (with the grid's cells ~a cell width between them), or merge them "
            "into one excavation."
        )


class _Cut(NamedTuple):
    """What :func:`_excavate_layers` returns: the cut arrays and the facts of the cut."""

    top: np.ndarray
    botm: np.ndarray
    idomain: np.ndarray
    excavated: np.ndarray
    absorbed: np.ndarray


def _excavate_layers(top, botm, idomain, records) -> _Cut:
    """Cut ``records`` into built layer arrays: the post-reconcile, post-pinch step.

    Runs AFTER reconcile and the pinch policy (design section 3.1): cutting before
    them was wrong under every pinch policy and reconcile mode measured. For each
    record cell with a finite bed ``b`` below its top: ``top = b``, every
    ``botm[k] = min(botm[k], b)``, so the layers above the bed have zero thickness and
    go idomain 0 -- never -1: the space is lake water, not a pass-through layer. Pinch
    verdicts are the pre-cut ones and are not re-run.

    Then the sliver rule (section 3.3), per column: while the first cell under the bed
    is inactive (idomain <= 0, a pass-through -1 included) -- a LAK vertical
    connection cannot reach through it -- or is thinner than the record's
    ``min_remnant`` with an active cell somewhere below it, it is ABSORBED: its bottom
    rises to the bed (zero thickness, idomain 0) and the cell below becomes the
    connected one, its top still exactly at the bed. Column thickness is conserved.
    A thin cell with nothing active below it is kept. A column left with no active
    cell under the bed -- a bed below ``botm[-1]`` included -- raises, naming the
    cells. A bed at or above the top (a rim cell) is not cut.

    Records must not overlap (:func:`_refuse_overlapping_excavations`).
    """

    top = np.array(top, dtype=float)
    botm = np.array(botm, dtype=float)
    idomain = np.array(idomain)
    nlay, ncpl = botm.shape
    excavated = np.zeros((nlay, ncpl), dtype=bool)
    absorbed = np.zeros((nlay, ncpl), dtype=float)
    stranded: dict[str, list[int]] = {}
    for record in records:
        rim, kept = [], []
        for cell, bed in zip(record.cells, record.bed, strict=True):
            if not np.isfinite(bed):
                continue
            if bed >= top[cell]:
                rim.append(cell)
                continue
            top[cell] = bed
            np.minimum(botm[:, cell], bed, out=botm[:, cell])
            below = np.flatnonzero(botm[:, cell] < bed)
            first = int(below[0]) if below.size else nlay
            excavated[:first, cell] = True
            idomain[:first, cell] = 0
            k = first
            # The top cell ``k`` had BEFORE any absorption above it: the bed for the
            # first cell under it, then the previous cell's own (unraised) bottom.
            own_top = bed
            while k < nlay:
                active_below = bool((idomain[k + 1:, cell] > 0).any())
                # What the cell would be if it connected: everything absorbed above it
                # is now part of it, so the min_remnant test is on this, from the bed.
                thickness = bed - botm[k, cell]
                if idomain[k, cell] > 0 and (thickness >= record.min_remnant or not active_below):
                    if thickness < record.min_remnant:
                        kept.append(cell)
                    break
                if not active_below:
                    k = nlay  # an inactive cell with nothing active below it
                    break
                # ...but what it CONTRIBUTES is only its own material: in a chain of
                # absorptions the cumulative figure counts every earlier remnant again.
                absorbed[k, cell] = own_top - botm[k, cell]
                own_top = botm[k, cell]
                botm[k, cell] = bed
                idomain[k, cell] = 0
                excavated[k, cell] = True
                k += 1
            if k >= nlay:
                stranded.setdefault(record.name, []).append(int(cell))
        if rim:
            _logger.info(
                "excavation %r: %d footprint cell(s) have a bed at or above the ground "
                "and are not cut (a rim; mf.lak caps its bed at the top there): %s",
                record.name, len(rim), rim[:10],
            )
        if kept:
            _logger.info(
                "excavation %r: %d connected cell(s) are thinner than min_remnant=%g and "
                "kept, because no active cell lies below them: %s",
                record.name, len(kept), record.min_remnant, kept[:10],
            )
    if stranded:
        detail = "; ".join(f"{name!r}: {cells[:10]}" for name, cells in stranded.items())
        raise ValueError(
            "An excavation leaves a column with no active cell under its bed (the bed "
            f"is below the stack's bottom, or everything below it is inactive): {detail}. "
            "A lake there would connect to nothing. Raise the bottom there, deepen the "
            "stack, or shrink the footprint."
        )
    return _Cut(top, botm, idomain, excavated, absorbed)


class LayerSurfaces:
    """An ordered stack of layer surfaces sampled onto a grid's cells.

    The first surface is the model top; each subsequent surface is the bottom of
    a layer. ``sample`` returns a centroid GeoDataFrame with one column per
    surface (matching ``gdf_topbtm``); ``top_botm`` returns DISV-ready arrays.
    """

    def __init__(self, surfaces: Sequence[Surface], labels: Sequence | None = None):
        """Wrap an ordered surface stack (top first) with matching column ``labels``.

        Nothing is sampled here; :meth:`sample` and the methods built on it do that.
        ``LayerStack`` compiles to this class (its labels are ``"top"`` plus each
        layer's name) and is the facade to build with.

        Parameters
        ----------
        surfaces : sequence of Surface
            The stack, top first: the first surface is the model top and every later
            one is the bottom of one layer, so ``n`` surfaces make ``n - 1`` layers.
            Copied into a list; the surfaces themselves are not copied or checked.
        labels : sequence of hashable, optional
            One column name per surface, in the same order, for the frame
            :meth:`sample` returns. ``None`` (the default) uses the integer positions
            ``0 .. n-1``, the integer columns ``gdf_topbtm`` has always carried. Each
            label must be distinct: a repeated one raises ``ValueError`` naming it
            (it kept only the LAST of its surfaces' values under every position
            naming it -- one layer bottom silently replaced by another).

        Raises
        ------
        ValueError
            ``labels`` has a different length from ``surfaces``, or names one column
            twice.
        """

        self.surfaces = list(surfaces)
        self.labels = (
            list(labels) if labels is not None else list(range(len(self.surfaces)))
        )
        if len(self.labels) != len(self.surfaces):
            raise ValueError("labels length must match surfaces length")
        # One column per surface (ledger 406). `sample` builds a dict keyed by
        # label, so a repeat kept the LAST surface's values under every position
        # naming it: labels t, a, a over 100, 90, 80 sampled 100, 80, 80.
        positions: dict = {}
        for index, label in enumerate(self.labels):
            positions.setdefault(label, []).append(index)
        repeated = {label: at for label, at in positions.items() if len(at) > 1}
        if repeated:
            named = "; ".join(
                f"{label!r} at positions {at}" for label, at in repeated.items()
            )
            raise ValueError(
                f"LayerSurfaces labels must be distinct, one column per surface: "
                f"{named}. A repeated label keeps only its last surface's values, "
                "replacing a layer bottom with another. Give each surface its own "
                "label (LayerStack names them 'top' plus each layer's name)."
            )

    def sample(
        self,
        vor,
        *,
        reconcile: bool = True,
        min_sep: float = 0.1,
        trigger_sep: float = 1,
        which: str = "bottom",
        method: str = "area",
        length_units: str | None = None,
        refresh: bool = False,
    ) -> gpd.GeoDataFrame:
        """Sample every surface onto ``vor`` cells and reconcile overlaps.

        Surfaces are resolved top-down so that each one can depend on the surface
        above it (relative surfaces and ``fill="propagate"``). With
        ``reconcile=True`` (default), a flat surface stays flat where it fits and
        is lowered where it would intrude on the surface above -- the "flat where
        possible, fit between" behavior.

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto.
        reconcile : bool, default True
            Enforce layer ordering where surfaces cross or crowd. ``False``
            returns the raw sampled surfaces, crossings and all.
        min_sep : float or sequence of float, default 0.1
            The separation reconcile sets where it acts: a scalar for every contact,
            or one value per layer (one fewer than the surfaces), else
            ``ValueError``.
        trigger_sep : float or sequence of float, default 1
            How close counts as a conflict: a contact at most this far apart is SET
            to exactly ``min_sep``, even if the surfaces never cross -- which THINS
            a layer between ``min_sep`` and ``trigger_sep`` thick (with the
            defaults a 0.8 ft layer comes back 0.1 ft, and under
            :meth:`top_botm_idomain` it then pinches whatever its minimum
            thickness). A scalar, or one value per layer as ``min_sep``.
            ``LayerStack`` caps it per layer so that a unit at or above its own
            ``min_thickness`` is never reached (ledger 377).
        which : {'bottom', 'top'}, default 'bottom'
            Which side gives way. ``'bottom'`` lowers the deeper surface to
            ``upper - min_sep`` and cascades downward, leaving the model top
            untouched. ``'top'`` raises the shallower one to ``lower + min_sep``
            and propagates upward, so it CAN MOVE THE MODEL TOP -- with
            ``top=100, a=105`` it returns ``top=105.1``. Prefer ``'bottom'``
            unless the deeper contact is the one you trust.
        method : {'area', 'centroid'}, default 'area'
            Raster sampling: area-weighted over the cell, or a single centroid
            lookup.
        length_units : str, optional
            The model length unit. Every surface declaring different ``units`` is
            converted to it, an operand inside an algebra surface included, before
            it combines with anything (:meth:`Surface.values`); ``units=`` on a
            relative or an algebra surface itself raises ``ValueError``, as does an
            unknown unit name. ``None`` (the default) declares no unit: a stack
            whose surfaces declare no ``units=`` is sampled as it is, and one that
            does raises ``ValueError`` naming ``length_units=`` (ledger
            433).
        refresh : bool, default False
            Rebuild any derived (contour) caches before sampling.

        Returns
        -------
        geopandas.GeoDataFrame
            Cell centroids with one column per surface, in stack order.

        Raises
        ------
        ValueError
            A surface that reads an ``mf.Basin`` (:meth:`Surface.basin`) is ``NaN`` at
            some cell -- off the lake, where a basin has no bed -- naming the surface and
            the cells (ledger 555); cut a lake with ``LayerStack.excavate``
            instead, or compose ``mf.Min(<surface>, <basin>)``. A unit-conversion or
            sampling error from :meth:`Surface.values`.
        """

        columns: dict = {}
        previous: np.ndarray | None = None
        for label, surface in zip(self.labels, self.surfaces, strict=False):
            vals = np.asarray(
                surface.values(
                    vor, previous=previous, method=method, refresh=refresh,
                    length_units=length_units,
                ),
                dtype=float,
            )
            _refuse_a_basin_gap(surface, vals, label=label)
            columns[label] = vals
            previous = vals

        gdf = gpd.GeoDataFrame(
            columns, geometry=vor.gdf_vorPolys.centroid.values, crs=vor.crs
        )
        gdf = gdf[["geometry"] + self.labels]
        if reconcile:
            df = reconcile_surfaces(
                vor, df=gdf, min_sep=min_sep, trigger_sep=trigger_sep, which=which
            )
            # Re-apply the ordering: `reconcile_surfaces` returns a plain frame
            # and `gpd.GeoDataFrame(df, geometry=...)` APPENDS geometry at the
            # end, undoing the `[["geometry"] + labels]` three lines above. That
            # made a keyword about physics move the column every positional
            # reader depends on -- see `geometry_first`.
            gdf = geometry_first(
                gpd.GeoDataFrame(df, geometry=gdf.geometry.values, crs=gdf.crs)
            )
        return gdf

    def attach(self, vor, **kwargs) -> gpd.GeoDataFrame:
        """Sample and store the result on ``vor.gdf_topbtm`` for downstream use.

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto (see :meth:`sample`) and to store the result on.
            Writing ``gdf_topbtm`` also CLEARS ``vor.excavations``: the grid's setter
            drops excavation records with the elevations they came with.
        **kwargs
            Forwarded to :meth:`sample` -- ``reconcile``, ``min_sep``,
            ``trigger_sep``, ``which``, ``method``, ``length_units``, ``refresh``. Any
            other keyword raises ``TypeError`` there, before ``vor`` is written.

        Returns
        -------
        geopandas.GeoDataFrame
            The sampled frame, as :meth:`sample` returns it.
        """

        gdf = self.sample(vor, **kwargs)
        vor.gdf_topbtm = gdf
        return gdf

    @staticmethod
    def _split_top_botm(gdf: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
        """Split a sampled centroid frame into ``(top, botm)`` arrays."""

        columns = [c for c in gdf.columns if c != "geometry"]
        elevations = gdf[columns].to_numpy().T  # (n_surfaces, n_cells)
        return elevations[0], elevations[1:]

    @staticmethod
    def _thickness(gdf: gpd.GeoDataFrame) -> np.ndarray:
        """Per-layer thickness ``(nlay, n_cells)`` from a sampled centroid frame."""

        columns = [c for c in gdf.columns if c != "geometry"]
        elev = gdf[columns].to_numpy().T  # (n_surfaces, n_cells)
        return elev[:-1] - elev[1:]

    @staticmethod
    def _idomain_from_thickness(
        thickness: np.ndarray, minimum_thickness=1.0, pinch="passthrough"
    ) -> np.ndarray:
        """Per-layer idomain from thickness, applying each layer's pinch policy.

        ``minimum_thickness`` and ``pinch`` may each be a scalar (applied to
        every layer) or a per-layer sequence. Policies: ``"passthrough"`` ->
        ``-1`` (vertical pass-through), ``"inactive"`` -> ``0`` (blocks flow),
        ``"floor"`` -> stay active (reconcile enforces minimum spacing).
        """

        nlay = thickness.shape[0]
        min_thk = _as_per_layer(minimum_thickness, nlay, "minimum_thickness")
        policy = _as_per_layer(pinch, nlay, "pinch")
        idomain = np.ones(thickness.shape, dtype=int)
        for i in range(nlay):
            if policy[i] not in _PINCH_POLICIES:
                raise ValueError(
                    f"Unknown pinch policy {policy[i]!r}; expected one of "
                    f"{_PINCH_POLICIES}."
                )
            thin = thickness[i] < float(min_thk[i])
            if policy[i] == "passthrough":
                idomain[i, thin] = -1
            elif policy[i] == "inactive":
                idomain[i, thin] = 0
            # "floor": leave active (idomain stays 1)
        return idomain

    @staticmethod
    def _validate_pinch_invariant(
        minimum_thickness, pinch, nlay: int, sample_kwargs: dict
    ) -> None:
        """Guard the fallback separation against each *pinching* layer's threshold.

        A ``"floor"`` layer has reconcile enforce its own ``min_thickness``, so it is
        exempt -- that is the whole mechanism by which ``min_thickness`` sets geometry.

        A ``"passthrough"`` or ``"inactive"`` layer is different: it must be allowed to
        come out thinner than its threshold, because that is the signal it pinches out.
        Reconcile gives it only the fallback ``min_sep``, so if the fallback is not
        below the threshold nothing can ever be thin enough to pinch and the policy is
        a silent no-op.

        ``min_sep`` is read as reconcile reads it -- a scalar or one value per layer,
        through the same ``_per_contact`` -- and each layer is checked against its
        OWN separation (ledger 377; a list raised ``TypeError`` here).
        """

        if not sample_kwargs.get("reconcile", True):
            return
        seps = _per_contact(sample_kwargs.get("min_sep", 0.1), nlay, "min_sep")
        min_thk = _as_per_layer(minimum_thickness, nlay, "minimum_thickness")
        policy = _as_per_layer(pinch, nlay, "pinch")
        for i in range(nlay):
            if policy[i] in ("passthrough", "inactive") and not (
                seps[i] < float(min_thk[i])
            ):
                raise ValueError(
                    f"layer {i}: the fallback min_sep ({seps[i]}) must be below its "
                    f"min_thickness ({min_thk[i]}) for pinch={policy[i]!r} to ever "
                    "fire -- reconcile would otherwise inflate every thin cell past "
                    "the threshold. Either lower min_sep, raise this layer's "
                    "min_thickness, or use pinch='floor' if you want the layer kept "
                    "at min_thickness rather than pinched out."
                )

    def top_botm(self, vor, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(top, botm)`` elevation arrays ready for ``mf.disv``.

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto (see :meth:`sample`). Its ``gdf_topbtm`` is not
            written; :meth:`attach` does that.
        **kwargs
            Forwarded to :meth:`sample` -- ``reconcile``, ``min_sep``,
            ``trigger_sep``, ``which``, ``method``, ``length_units``, ``refresh``.

        Returns
        -------
        top : numpy.ndarray of float, shape (ncpl,)
            The first surface, sampled (and reconciled unless ``reconcile=False``).
        botm : numpy.ndarray of float, shape (nlay, ncpl)
            Every later surface, one row per layer, top layer first.
        """

        return self._split_top_botm(self.sample(vor, **kwargs))

    def top_botm_idomain(
        self,
        vor,
        *,
        minimum_thickness=1.0,
        pinch="passthrough",
        **sample_kwargs,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(top, botm, idomain)`` with thin cells pinched out.

        ``minimum_thickness`` and ``pinch`` may be scalars (every layer) or
        per-layer sequences. Default policy ``"passthrough"`` gives thin cells
        ``idomain = -1`` (MODFLOW 6 vertical pass-through).

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto (see :meth:`sample`). Its ``gdf_topbtm`` is not
            written.
        minimum_thickness : float or sequence of float, default 1.0
            The thickness, in the sampled (converted) length units, below which a
            cell is thin -- strictly below, so a cell exactly this thick stays
            active. A scalar applies to every layer; a sequence gives one value per
            layer and must have one entry per bottom surface, or ``ValueError``.
        pinch : {'passthrough', 'inactive', 'floor'} or sequence of those, default 'passthrough'
            What a thin cell becomes: ``'passthrough'`` -> ``idomain = -1`` (removed,
            vertical flow passes through it), ``'inactive'`` -> ``0`` (removed, blocks
            flow), ``'floor'`` -> stays ``1`` at its sampled thickness. A scalar
            applies to every layer, a sequence gives one per layer (same length rule
            as ``minimum_thickness``). Any other value raises ``ValueError``. A
            ``'passthrough'`` or ``'inactive'`` layer is also checked against the
            reconcile separation (see ``**sample_kwargs``).
        **sample_kwargs
            Forwarded to :meth:`sample` -- ``reconcile``, ``min_sep``,
            ``trigger_sep``, ``which``, ``method``, ``length_units``, ``refresh``.
            ``reconcile`` and ``min_sep`` are read here too: while reconcile is on
            (the default), every ``'passthrough'``/``'inactive'`` layer needs
            ``min_sep`` (default 0.1) strictly below its ``minimum_thickness``,
            because reconcile would otherwise open every thin contact past the
            threshold and the policy could never fire -- that raises ``ValueError``.
            A per-layer ``min_sep`` is checked layer by layer, each pinching layer
            against its own separation.

        Returns
        -------
        top : numpy.ndarray of float, shape (ncpl,)
            The sampled model top.
        botm : numpy.ndarray of float, shape (nlay, ncpl)
            The sampled layer bottoms, top layer first.
        idomain : numpy.ndarray of int, shape (nlay, ncpl)
            ``1`` active, ``0`` or ``-1`` for a thin cell as ``pinch`` says.

        Raises
        ------
        ValueError
            A per-layer ``minimum_thickness`` or ``pinch`` of the wrong length, an
            unknown ``pinch``, or ``min_sep`` not below a pinching layer's
            ``minimum_thickness``.
        """

        gdf = self.sample(vor, **sample_kwargs)
        top, botm = self._split_top_botm(gdf)
        thickness = self._thickness(gdf)
        self._validate_pinch_invariant(
            minimum_thickness, pinch, thickness.shape[0], sample_kwargs
        )
        idomain = self._idomain_from_thickness(thickness, minimum_thickness, pinch)
        return top, botm, idomain

    def thickness_report(
        self, vor, *, minimum_thickness: float = 1.0, **sample_kwargs
    ) -> str:
        """Per-layer thickness summary and thin-cell counts (for diagnostics).

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto (see :meth:`sample`). Its ``gdf_topbtm`` is not
            written.
        minimum_thickness : float or sequence of float, default 1.0
            The threshold a cell is counted as thin below (strictly), in the sampled
            (converted) length units. Only counted and printed -- nothing is pinched.
            A scalar applies to every layer; a sequence gives one value per layer
            and must have one entry per bottom surface, or ``ValueError``.
        **sample_kwargs
            Forwarded to :meth:`sample` -- ``reconcile``, ``min_sep``,
            ``trigger_sep``, ``which``, ``method``, ``length_units``, ``refresh`` --
            so the report describes the reconciled thicknesses unless
            ``reconcile=False`` is passed.

        Returns
        -------
        str
            A header line (layer and cell counts) and one line per layer with its
            min/mean/max thickness and its thin-cell count. Returned, not printed.
        """

        thickness = self._thickness(self.sample(vor, **sample_kwargs))
        nlay, ncell = thickness.shape
        min_thk = _as_per_layer(minimum_thickness, nlay, "minimum_thickness")
        lines = [f"LayerSurfaces: {nlay} layers, {ncell} cells"]
        for i in range(nlay):
            t = thickness[i]
            thin = int((t < float(min_thk[i])).sum())
            lines.append(
                f"  layer {i}: thickness min={t.min():.2f} mean={t.mean():.2f} "
                f"max={t.max():.2f}  thin={thin} (<{min_thk[i]})"
            )
        return "\n".join(lines)

    def to_disv(
        self,
        vor,
        *,
        nlay: int | None = None,
        idomain: Any | None = None,
        pinch_out: bool = False,
        minimum_thickness=1.0,
        pinch="passthrough",
        length_units: str | None = None,
        name: str = "disv",
        attach: bool = True,
        **sample_kwargs,
    ):
        """Build a ready-to-use ``mf.disv`` package spec from ``vor`` + surfaces.

        Samples the surfaces onto ``vor`` (reconciling overlaps by default),
        takes the cell geometry from ``vor.get_disv_gridprops()``, and returns a
        :class:`~myflopy.specs.PackageSpec` for DISV. With ``pinch_out=True``,
        cells thinner than ``minimum_thickness`` get ``idomain = -1`` (vertical
        pass-through). With ``attach=True`` (default) the sampled elevations are
        also written to ``vor.gdf_topbtm`` so choropleth/mapping code keeps
        working.

        Parameters
        ----------
        vor : VoronoiGridPlus
            The grid to sample onto (see :meth:`sample`); its
            ``get_disv_gridprops()`` supplies ``ncpl``, ``vertices`` and ``cell2d``
            (``nvert`` is the vertex count). Written to only with ``attach=True``.
        nlay : int, optional
            The layer count written. ``None`` (the default) writes
            ``len(surfaces) - 1``, one layer per bottom surface. A given value is
            passed to ``mf.disv(nlay=)`` as it is: nothing is resampled, so ``botm``
            (and a derived ``idomain``) still has one row per bottom surface.
        idomain : array_like of int, shape (nlay, ncpl), optional
            An explicit IDOMAIN, forwarded to ``mf.disv(idomain=)``: ``1`` active,
            ``0`` inactive, ``-1`` vertical pass-through. ``None`` (the default)
            with ``pinch_out=False`` writes none, and MF6 treats every cell as
            active. Refused together with ``pinch_out=True`` (``ValueError``).
        pinch_out : bool, default False
            Derive IDOMAIN from each layer's sampled thickness, with
            ``minimum_thickness`` and ``pinch``, exactly as
            :meth:`top_botm_idomain` does. ``False`` ignores both of them.
        minimum_thickness : float or sequence of float, default 1.0
            With ``pinch_out=True`` only: the thickness, in the sampled (converted)
            length units, strictly below which a cell is thin. A scalar applies to
            every layer; a sequence gives one per layer and must have one entry per
            bottom surface, or ``ValueError``.
        pinch : {'passthrough', 'inactive', 'floor'} or sequence of those, default 'passthrough'
            With ``pinch_out=True`` only: what a thin cell becomes.
            ``'passthrough'`` writes ``idomain = -1``; ``'inactive'`` writes ``0``;
            ``'floor'`` leaves it ``1`` at its sampled thickness. Scalar or per
            layer, as ``minimum_thickness``; any other value raises ``ValueError``.
        length_units : str, optional
            The model length unit, which does two things. It is passed to
            :meth:`sample`, converting every surface of the stack that declares its
            own ``units``, an operand inside an algebra surface included (a relative
            or algebra surface declaring them raises ``ValueError`` there). And it is
            written as DISV ``LENGTH_UNITS`` in MF6's own word: every spelling the
            conversion takes -- ``'feet'``/``'ft'``/``'foot'``/``'us-ft'``/
            ``'usft'``, ``'meters'``/``'meter'``/``'metres'``/``'m'``,
            ``'centimeters'``/``'centimeter'``/``'cm'``, any case -- becomes ``FEET``, ``METERS`` or
            ``CENTIMETERS``. Any other name raises ``ValueError`` naming those,
            before anything is sampled or attached: MF6 reads a word it does not
            know (``FT`` was written for ``'ft'``) as length unit 0, undefined,
            without a message. ``None`` (the default) declares no unit: it writes no
            ``LENGTH_UNITS``, and a surface declaring ``units=`` raises
            ``ValueError`` in :meth:`sample` (ledger 433).
        name : str, default "disv"
            Forwarded to ``mf.disv(name=)``: the returned spec's name, the key
            ``with_packages(slot=)`` and ``without_packages`` use.
        attach : bool, default True
            Also write the sampled frame to ``vor.gdf_topbtm`` (which clears
            ``vor.excavations``). Written LAST, once the spec is built: a call that
            raises -- an unknown ``length_units``, ``pinch_out`` beside an
            ``idomain``, a pinch check, a refusal of ``mf.disv`` -- leaves ``vor``
            as it was. ``False`` leaves ``vor`` untouched.
        **sample_kwargs
            Forwarded to :meth:`sample` -- ``reconcile``, ``min_sep``,
            ``trigger_sep``, ``which``, ``method``, ``refresh`` (``length_units`` is
            this method's own). With ``pinch_out=True``, ``reconcile`` and
            ``min_sep`` are also read by the separation check that
            :meth:`top_botm_idomain` describes. Any other keyword raises
            ``TypeError`` in :meth:`sample`.

        Returns
        -------
        PackageSpec
            The hand-written ``mf.disv`` form: the mesh, ``top``, ``botm``, the
            ``idomain`` when there is one, and ``length_units`` when given.

        Raises
        ------
        ValueError
            An unknown ``length_units``; ``pinch_out=True`` together with an
            ``idomain``; with ``pinch_out=True``, a ``minimum_thickness``/``pinch``
            of the wrong length, an unknown ``pinch``, or ``min_sep`` not below a
            pinching layer's threshold; a cell with ``idomain != 0`` (every cell,
            with no ``idomain``) left with no finite top or bottom -- a raster
            covering only part of the grid -- naming the surfaces and the cells
            (:func:`refuse_unbounded_cells`; ledger 587). None of them has
            written ``vor``.
        """

        from myflopy.package_api import disv  # lazy: avoid an import cycle

        # Every refusal runs before `attach` writes the grid (ledger 407): the
        # write came straight after sampling, so a call refused below had already
        # replaced `vor.gdf_topbtm` -- and its setter cleared `vor.excavations`.
        if pinch_out and idomain is not None:
            raise ValueError(
                "Pass either pinch_out=True or an explicit idomain, not both."
            )
        disv_options = {}
        if length_units is not None:
            disv_options["length_units"] = _mf6_length_units(length_units)
        gdf = self.sample(vor, length_units=length_units, **sample_kwargs)
        top, botm = self._split_top_botm(gdf)
        if pinch_out:
            thickness = self._thickness(gdf)
            self._validate_pinch_invariant(
                minimum_thickness, pinch, thickness.shape[0], sample_kwargs
            )
            idomain = self._idomain_from_thickness(
                thickness, minimum_thickness, pinch
            )
        refuse_unbounded_cells(top, botm, idomain, labels=self.labels, surfaces=self.surfaces)
        props = vor.get_disv_gridprops()
        spec = disv(
            nlay=len(self.surfaces) - 1 if nlay is None else nlay,
            ncpl=props["ncpl"],
            nvert=len(props["vertices"]),
            vertices=props["vertices"],
            cell2d=props["cell2d"],
            top=top,
            botm=botm,
            idomain=idomain,
            name=name,
            **disv_options,
        )
        if attach:
            vor.gdf_topbtm = gdf
        return spec


__all__ = ["ExcavationRecord", "Surface", "LayerSurfaces", "sample_bed"]
