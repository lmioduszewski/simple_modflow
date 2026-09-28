from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import shapely as shp
from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPolygon, Polygon

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.helpers import densify_poly

logger = get_logger(__name__)

#: A part the profile's simplification changes by more than this fraction of its area is
#: narrower than the tolerance and is kept as drawn -- boundary parts only (ledger 316).
#: The real boundaries measured moved 0-10%; a collapse moves ~50% (a quad becomes a triangle).
COLLAPSE_FRACTION = 0.25


@dataclass(slots=True)
class CleanupStats:
    polygons_in: int
    polygons_out: int
    vertices_in: int
    vertices_out: int
    vanished: tuple = ()  # outline length of each part resampling left with no area


def iter_polygons(geometry: Polygon | MultiPolygon | GeometryCollection):
    """Yield each constituent ``Polygon`` of a polygon/multipolygon/collection (raises otherwise)."""

    if isinstance(geometry, Polygon):
        yield geometry
        return
    if isinstance(geometry, MultiPolygon):
        yield from geometry.geoms
        return
    if isinstance(geometry, GeometryCollection):
        for geom in geometry.geoms:
            if isinstance(geom, Polygon):
                yield geom
            elif isinstance(geom, MultiPolygon):
                yield from geom.geoms
        return
    raise TypeError("geometry must be a Polygon, MultiPolygon, or GeometryCollection")


def polygon_vertex_count(geometry: Polygon | MultiPolygon | GeometryCollection) -> int:
    """Total vertex count across all exterior and interior rings of a polygonal geometry."""

    total = 0
    for polygon in iter_polygons(geometry):
        total += len(polygon.exterior.coords)
        total += sum(len(interior.coords) for interior in polygon.interiors)
    return total


def _resample_linestring(line: LineString, target_segment_length: float) -> LineString:
    """Re-sample a line to roughly even ``target_segment_length`` segments (no-op if short/disabled)."""

    if target_segment_length <= 0 or line.length <= target_segment_length:
        return line

    distances = np.arange(0.0, line.length, target_segment_length, dtype=float)
    coords = [line.interpolate(distance).coords[0] for distance in distances]
    coords.append(line.coords[-1])

    deduped: list[tuple[float, float]] = []
    for x, y in coords:
        point = (float(x), float(y))
        if not deduped or point != deduped[-1]:
            deduped.append(point)

    if len(deduped) < 2:
        return line
    return LineString(deduped)


def resample_polygon_boundaries(
    geometry: Polygon | MultiPolygon,
    *,
    target_segment_length: float,
) -> Polygon | MultiPolygon:
    """Densify each polygon's boundary to ~``target_segment_length`` spacing (no-op if disabled)."""

    if target_segment_length <= 0:
        return geometry

    polygons: list[Polygon] = []
    for polygon in iter_polygons(geometry):
        resampled = densify_poly(polygon, target_segment_length)
        polygons.append(resampled.buffer(0))

    if len(polygons) == 1:
        return polygons[0]
    return MultiPolygon(polygons)


def resample_linear_geometry(
    geometry: LineString | MultiLineString,
    *,
    target_segment_length: float,
) -> LineString | MultiLineString:
    """Re-sample a line/multiline to ~``target_segment_length`` segments (no-op if disabled)."""

    if target_segment_length <= 0:
        return geometry
    if isinstance(geometry, LineString):
        return _resample_linestring(geometry, target_segment_length)
    if isinstance(geometry, MultiLineString):
        return MultiLineString(
            [_resample_linestring(line, target_segment_length) for line in geometry.geoms]
        )
    raise TypeError("geometry must be a LineString or MultiLineString")


def cleanup_polygonal_geometry(
    geometry: Polygon | MultiPolygon,
    *,
    snap_tolerance: float = 0.0,
    simplify_tolerance: float | None = None,
    min_feature_area: float = 0.0,
    target_segment_length: float | None = None,
    keep_collapsing_parts: bool = False,
) -> tuple[Polygon | MultiPolygon, CleanupStats]:
    """Clean a polygonal geometry (snap, simplify, drop tiny features, re-sample) with before/after stats.

    ``keep_collapsing_parts=True`` (the model boundary) keeps as drawn, with a warning, any
    part the simplification would change by more than ``COLLAPSE_FRACTION`` of its area: a
    part narrower than ``simplify_tolerance`` loses a corner and half its area otherwise --
    a 1 ft square became a triangle meshed into 20,681 cells (ledger 316).
    """

    polygons_in = sum(1 for _ in iter_polygons(geometry))
    vertices_in = polygon_vertex_count(geometry)

    cleaned = geometry.buffer(0)
    if snap_tolerance and snap_tolerance > 0:
        cleaned = shp.snap(cleaned, cleaned, snap_tolerance)

    if simplify_tolerance and simplify_tolerance > 0:
        simplified = cleaned.simplify(simplify_tolerance, preserve_topology=True)
        before, after = list(iter_polygons(cleaned)), list(iter_polygons(simplified))
        if keep_collapsing_parts and len(before) == len(after):
            kept = [i for i, (b, a) in enumerate(zip(before, after, strict=True))
                    if b.area > 0 and abs(a.area - b.area) > COLLAPSE_FRACTION * b.area]
            if kept:
                where = ", ".join(
                    f"the part at ({before[i].bounds[0]:.10g}, {before[i].bounds[1]:.10g})-"
                    f"({before[i].bounds[2]:.10g}, {before[i].bounds[3]:.10g}), "
                    f"{before[i].area:.4g} -> {after[i].area:.4g}" for i in kept
                )
                message = (
                    f"simplify_tolerance={simplify_tolerance:g} would change the area of "
                    f"{len(kept)} boundary part(s) by more than {COLLAPSE_FRACTION:.0%} ({where}): "
                    "narrower than the tolerance, they are kept as drawn."
                )
                logger.warning(message)
                import warnings

                warnings.warn(message, UserWarning, stacklevel=2)
                parts = [before[i] if i in kept else a for i, a in enumerate(after)]
                simplified = parts[0] if len(parts) == 1 else MultiPolygon(parts)
        cleaned = simplified

    cleaned = cleaned.buffer(0)

    polygons = [polygon.buffer(0) for polygon in iter_polygons(cleaned)]
    polygons = [polygon for polygon in polygons if not polygon.is_empty and polygon.area > 0]

    if min_feature_area > 0:
        filtered = [polygon for polygon in polygons if polygon.area >= min_feature_area]
        polygons = filtered or polygons

    if not polygons:
        raise ValueError("geometry cleanup removed all polygonal area")

    cleaned = polygons[0] if len(polygons) == 1 else MultiPolygon(polygons)

    vanished = ()
    if target_segment_length and target_segment_length > 0:
        vanished = tuple(
            part.exterior.length for part in iter_polygons(cleaned)
            if densify_poly(part, target_segment_length).buffer(0).area <= 0
        )
        cleaned = resample_polygon_boundaries(
            cleaned,
            target_segment_length=target_segment_length,
        ).buffer(0)

    polygons_out = sum(1 for _ in iter_polygons(cleaned))
    vertices_out = polygon_vertex_count(cleaned)

    return cleaned, CleanupStats(
        polygons_in=polygons_in,
        polygons_out=polygons_out,
        vertices_in=vertices_in,
        vertices_out=vertices_out,
        vanished=vanished,
    )
