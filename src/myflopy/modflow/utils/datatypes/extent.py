"""``map(extent=...)`` -- the window a map opens on, resolved once (ledger 309).

One resolver for every map, scope and backend, so what a box MEANS cannot differ
between them. An extent sets the VIEW and nothing else (the user's decision): no
cell is cropped and the colour scale still spans the whole field -- ``zmin`` /
``zmax`` is how to tighten it for a close-up.

Accepted, each resolving to one box in the GRID's CRS and the same box in WGS84:

* four numbers ``(xmin, ymin, xmax, ymax)`` -- shapely/GeoPandas order, as
  ``LayerStack.draft_grid(extent=)`` uses; NOT FloPy's or matplotlib's
  ``(xmin, xmax, ymin, ymax)``;
* a run of ``(x, y)`` points -- their bounding box;
* any shapely geometry -- its bounds (a polygon, a line along a transect ...);
* a ``GeoSeries``/``GeoDataFrame`` -- reprojected from its own CRS;
* an area object answering ``to_geoseries()`` (``mf.Area``);
* the NAME of an area of interest, looked up in the map's own ``names``.

Numbers and bare geometries are in the grid's CRS, the rule ``features=`` uses; a
box that misses the grid RAISES, printing both boxes and both CRSs, which is what
a lon/lat box passed to a feet-based model looks like.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Real

import geopandas as gpd
import numpy as np
import shapely as shp

from myflopy._defaults import DEFAULT_CRS
from myflopy.viz import is_extent_list

__all__ = ["Extent", "is_extent_list", "resolve_extent"]


@dataclass(frozen=True)
class Extent:
    """A resolved view window: one box in the grid's CRS and the same in WGS84."""

    native: tuple[float, float, float, float]
    wgs84: tuple[float, float, float, float]
    label: str | None = None

    def padded(self, padding: float) -> tuple[float, float, float, float]:
        """The native box grown by ``padding`` of its own size on every side."""

        xmin, ymin, xmax, ymax = self.native
        dx = (xmax - xmin) * float(padding)
        dy = (ymax - ymin) * float(padding)
        return xmin - dx, ymin - dy, xmax + dx, ymax + dy


def _is_number(value) -> bool:
    return isinstance(value, Real) and not isinstance(value, bool)


def _numbers(value) -> list[float] | None:
    """``value`` as a flat list of numbers, or None if it is not a sequence of them."""

    if isinstance(value, (str, bytes)) or not isinstance(value, (Sequence, np.ndarray)):
        return None
    flat = np.asarray(value, dtype=object).reshape(-1).tolist()
    if not flat or not all(_is_number(item) for item in flat):
        return None
    return [float(item) for item in flat]


def _box_of(value, *, crs, names):
    """``(geometry, its crs, label)`` for one accepted spelling."""

    if isinstance(value, str):
        if names is None or value not in names:
            known = sorted(names) if names is not None else []
            hint = (
                f" Known here: {known}." if known
                else " This map has no named areas; define one with model.aoi / "
                "project.aoi, or pass a box or a geometry."
            )
            raise KeyError(f"extent={value!r} names no area of interest.{hint}")
        area, label = names[value], value
        return (*_box_of(area, crs=crs, names=None)[:2], label)
    if isinstance(value, (gpd.GeoSeries, gpd.GeoDataFrame)):
        frame = value if value.crs is not None else value.set_crs(crs)
        return shp.box(*frame.to_crs(crs).total_bounds), crs, None
    if hasattr(value, "to_geoseries"):
        series = value.to_geoseries()
        label = getattr(value, "name", None)
        return (*_box_of(series, crs=crs, names=None)[:2], label)
    if hasattr(value, "geom_type"):
        return shp.box(*value.bounds), crs, None
    numbers = _numbers(value)
    if numbers is not None:
        shape = np.asarray(value, dtype=float).shape
        if len(shape) == 1 and shape[0] == 4:
            xmin, ymin, xmax, ymax = numbers
            if not (xmin < xmax and ymin < ymax):
                raise ValueError(
                    f"extent={tuple(numbers)} is not (xmin, ymin, xmax, ymax) with "
                    "xmin < xmax and ymin < ymax. The order is shapely's -- NOT "
                    "FloPy's or matplotlib's (xmin, xmax, ymin, ymax)."
                )
            return shp.box(xmin, ymin, xmax, ymax), crs, None
        if len(shape) == 2 and shape[1] == 2:
            return shp.box(*shp.MultiPoint(np.asarray(value, dtype=float)).bounds), crs, None
    raise TypeError(
        f"extent={value!r} is not an extent: pass (xmin, ymin, xmax, ymax), a list of "
        "(x, y) points, a shapely geometry, a GeoDataFrame, an area, or an area's name."
    )


def resolve_extent(value, *, crs=None, bounds=None, names: Mapping | None = None) -> Extent | None:
    """One accepted ``extent=`` spelling as an :class:`Extent`, or None for None.

    ``crs`` is the grid's (the frame numbers and bare geometries are read in);
    ``bounds`` its ``(xmin, ymin, xmax, ymax)``, which the extent must overlap;
    ``names`` maps area-of-interest names to areas.
    """

    if value is None:
        return None
    if is_extent_list(value):
        raise TypeError(
            "extent= got a LIST of extents -- one per panel, which only a mosaic can "
            "use. A single map takes one: a box, a geometry or a name."
        )
    crs = crs if crs is not None else DEFAULT_CRS
    box, _, label = _box_of(value, crs=crs, names=names)
    xmin, ymin, xmax, ymax = (float(v) for v in box.bounds)
    if not (xmin < xmax and ymin < ymax):
        raise ValueError(
            f"extent={value!r} has no area -- bounds ({xmin}, {ymin}, {xmax}, {ymax}). "
            "A point, or a line running straight along one axis, frames nothing; "
            "pass a box around it."
        )
    if bounds is not None:
        gxmin, gymin, gxmax, gymax = (float(v) for v in bounds)
        if xmax < gxmin or xmin > gxmax or ymax < gymin or ymin > gymax:
            raise ValueError(
                "extent= lies entirely off the grid, so the map would open on empty "
                f"space. extent ({xmin:.6g}, {ymin:.6g}, {xmax:.6g}, {ymax:.6g}) vs "
                f"grid ({gxmin:.6g}, {gymin:.6g}, {gxmax:.6g}, {gymax:.6g}), both read "
                f"in {crs}. Numbers and bare geometries are taken to be in the "
                "GRID's CRS; a box in another one belongs in a GeoDataFrame that "
                "declares it."
            )
    # Densified before reprojecting: a box's edges are not straight in WGS84.
    edge = max(xmax - xmin, ymax - ymin) / 32.0
    outline = gpd.GeoSeries([shp.box(xmin, ymin, xmax, ymax)], crs=crs).segmentize(edge)
    west, south, east, north = (float(v) for v in outline.to_crs(epsg=4326).total_bounds)
    return Extent((xmin, ymin, xmax, ymax), (west, south, east, north), label)
