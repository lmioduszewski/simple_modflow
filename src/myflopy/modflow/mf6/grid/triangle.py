"""Triangle-based mesh setup utilities for Voronoi MODFLOW workflows.

The public surface in this module is centered on :class:`TriangleGrid`, which
wraps FloPy's Triangle helper with a more intuitive workflow for defining a
domain, adding refinement regions from code or GIS data, cleaning the geometry,
and optionally running a conservative CVT/Lloyd-style smoothing pass.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Any
import warnings

import geopandas as gpd
import numpy as np
import shapely as shp
from flopy.utils.geospatial_utils import GeoSpatialUtil
from flopy.utils.triangle import Triangle
from shapely.geometry import GeometryCollection, MultiPoint, MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry

from myflopy._logging import get_logger
from myflopy.grid_identity import written_resolution
from myflopy.modflow.mf6.grid.geometry_cleanup import cleanup_polygonal_geometry
from myflopy.modflow.mf6.grid.helpers import densify_poly
from myflopy.modflow.mf6.grid.mesh_quality import triangle_quality_report
from myflopy.modflow.mf6.grid.seed_optimization import (
    cvt_relaxed_full_seed_points,
    get_fixed_and_free_vertex_ids,
)
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg


logger = get_logger(__name__)


class TriangleFailed(subprocess.CalledProcessError):
    """The Triangle program crashed, ran out of memory, hung or wrote without end.

    A ``CalledProcessError``, as FloPy's own run raised, so the smoothing loop's fallback
    still catches it. A point constraint a fraction of an ulp inside an edge made Triangle
    SIGSEGV, hang, or grow by 0.3-0.5 GB a second until the machine ran out (ledger 316),
    so every run is now bounded: a memory cap, a timeout and a cap on its output.
    """

    def __str__(self) -> str:
        return str(self.output or super().__str__())


class NoAreaInside(ValueError):
    """A boundary or region with no polygonal area left to mesh.

    ``cause`` says which step removed it: ``"outside"`` (a region with no area inside the
    model boundary -- disjoint from it, or touching it only along a line or at a point),
    ``"clipped"`` (a buffered point or line outside the boundary but within its buffer of
    it: it is clipped to the boundary BEFORE it is buffered, so nothing is left to
    buffer), ``"buffer"`` (a negative buffer at least half its width) or ``"empty"``
    (zero area as given). ``GridSpec.voronoi`` skips an ``"outside"`` feature with a warning when others
    in its source land, and names the rest (ledger 316).
    """

    def __init__(self, message: str = "", *, cause: str = "empty"):
        super().__init__(message)
        self.cause = cause

    def __reduce__(self):
        return (type(self), (str(self),), {"cause": self.cause})


@dataclass(slots=True)
class DomainSpec:
    """Describes the model domain and its optional background cell sizing."""

    geometry: Polygon | MultiPolygon
    label: str = "domain"
    max_area: float | int | None = None
    attribute: int | float = 0


@dataclass(slots=True)
class RegionSpec:
    """Describes a requested refinement region before Triangle inputs are built.

    ``segments=False`` hands Triangle the region's POINT only, never its ring: Triangle then
    floods ``max_area`` through everything no other ring encloses. Set on the background
    sizing region when ``max_area`` is too large for the domain to hold one (ledger 316).
    """

    geometry: Polygon | MultiPolygon
    max_area: float | int
    label: str | None = None
    priority: int = 0
    attribute: int | float = 0
    point_hint: tuple[float, float] | None = None
    source: str = "geometry"
    segments: bool = True


@dataclass(slots=True)
class PreparedRegion:
    """Stores a region after overlaps and interior points are resolved."""

    geometry: Polygon
    claim_geometry: Polygon
    point: tuple[float, float]
    max_area: float | int
    label: str
    priority: int = 0
    attribute: int | float = 0
    source: str = "geometry"
    segments: bool = True


@dataclass(slots=True)
class MeshBuildProfile:
    """A named bundle of mesh cleanup + smoothing settings for triangulation.

    Groups the knobs that control how :meth:`TriangleGrid.build_mesh` cleans input
    geometry and relaxes the triangulation -- snapping/simplifying tolerances,
    minimum feature area, target segment length, and the constrained CVT/Lloyd
    relaxation of :meth:`TriangleGrid.optimize_seeds` (iterations, damping, minimum
    move). Rather than tuning each by hand, pick a built-in preset with
    :meth:`from_name` (``"fast"``, ``"balanced"``, ``"smooth"``, ``"clean"``,
    ``"baseline"``) and override individual fields as needed -- either by building
    your own instance or by passing the field as a keyword to
    :meth:`TriangleGrid.build_mesh`. Higher quality (more iterations, shorter
    segments) costs build time.

    The six cleanup fields are the keywords of :meth:`TriangleGrid.clean_geometry`
    (see :meth:`cleanup_kwargs`); lengths are in model length units and areas in
    model length units squared.

    Parameters
    ----------
    name : str
        Profile label. :meth:`TriangleGrid.build_mesh` returns it as the
        ``"profile"`` entry of its report; it does not change the mesh.
    cleanup : bool, optional
        Run :meth:`TriangleGrid.clean_geometry` before triangulating. Default
        ``True``. ``False`` skips cleanup, so the six cleanup fields below are
        not used.
    snap_tolerance : float, optional
        Distance within which the vertices and segments of the domain and of each
        region are snapped to that geometry's own vertices (``shapely.snap`` of
        the geometry to itself). Default ``0.0``; ``0`` or less disables snapping.
    simplify_tolerance : float or None, optional
        ``shapely`` simplification tolerance for the domain and every region,
        with topology preserved. Default ``2.0``; ``None`` or ``0`` disables it.
    min_feature_area : float, optional
        Polygon parts smaller than this area are dropped from the domain and
        from each region -- unless that would drop every part, in which case all
        are kept. Default ``0.0`` (keep everything).
    target_segment_length : float or None, optional
        Spacing to resample boundaries to. The ring is rebuilt from points spaced
        this far apart along it, so vertices between them (corners included)
        are not kept, and interior rings (holes) are dropped. Applied only where
        ``resample_domain_boundary`` / ``resample_region_sources`` select it.
        Default ``120.0``; ``None`` or ``0`` disables resampling.
    resample_domain_boundary : bool, optional
        Resample the domain boundary to ``target_segment_length`` as well.
        Default ``False``. Resampling drops the domain's holes.
    resample_region_sources : tuple of str, optional
        Region ``source`` tags whose boundaries are resampled to
        ``target_segment_length`` (``"line"`` is the tag
        :meth:`TriangleGrid.add_line_feature` gives). Default ``("line",)``.
    optimization_iterations : int, optional
        Maximum number of relax-and-rebuild rounds of
        :meth:`TriangleGrid.optimize_seeds`. Default ``3``. With
        ``build_mesh(optimize=None)``, ``0`` means no optimization pass.
    damping : float, optional
        Fraction of the way each free mesh vertex moves toward its Voronoi cell
        centroid per round; values outside ``[0, 1]`` are clipped. Default
        ``0.25``.
    min_move : float, optional
        A vertex whose proposed move is shorter than this distance stays put.
        Default ``1e-3``.
    """

    name: str
    cleanup: bool = True
    snap_tolerance: float = 0.0
    simplify_tolerance: float | None = 2.0
    min_feature_area: float = 0.0
    target_segment_length: float | None = 120.0
    resample_domain_boundary: bool = False
    resample_region_sources: tuple[str, ...] = ("line",)
    optimization_iterations: int = 3
    damping: float = 0.25
    min_move: float = 1e-3

    @classmethod
    def from_name(cls, profile: str):
        """Return one of the built-in mesh build presets by name.

        Every preset keeps the field defaults (``simplify_tolerance=2.0``, no
        snapping, no minimum feature area, only ``"line"`` regions resampled,
        ``min_move=1e-3``) except:

        ============  ======================  ===========  =======  =======
        name          target_segment_length   iterations   damping  cleanup
        ============  ======================  ===========  =======  =======
        ``fast``      150                     1            0.2      yes
        ``balanced``  120                     3            0.25     yes
        ``smooth``    100                     5            0.3      yes
        ``clean``     120                     0            0.0      yes
        ``baseline``  None                    0            0.0      no
        ============  ======================  ===========  =======  =======

        Parameters
        ----------
        profile : str
            Preset name, one of ``"fast"``, ``"balanced"``, ``"smooth"``,
            ``"clean"`` or ``"baseline"``. Surrounding whitespace and case are
            ignored.

        Returns
        -------
        MeshBuildProfile
            A new instance holding the preset's settings.

        Raises
        ------
        ValueError
            If ``profile`` names no preset.
        """
        normalized = profile.strip().lower()
        presets = {
            "fast": cls(
                name="fast",
                target_segment_length=150.0,
                optimization_iterations=1,
                damping=0.2,
            ),
            "balanced": cls(
                name="balanced",
                target_segment_length=120.0,
                optimization_iterations=3,
                damping=0.25,
            ),
            "smooth": cls(
                name="smooth",
                target_segment_length=100.0,
                optimization_iterations=5,
                damping=0.3,
            ),
            "clean": cls(
                name="clean",
                target_segment_length=120.0,
                optimization_iterations=0,
                damping=0.0,
            ),
            "baseline": cls(
                name="baseline",
                cleanup=False,
                target_segment_length=None,
                optimization_iterations=0,
                damping=0.0,
            ),
        }
        if normalized not in presets:
            valid = ", ".join(sorted(presets))
            raise ValueError(f"Unknown mesh profile {profile!r}. Expected one of: {valid}")
        return presets[normalized]

    def cleanup_kwargs(self) -> dict[str, object]:
        """Return cleanup-related settings as keyword arguments."""
        return {
            "snap_tolerance": self.snap_tolerance,
            "simplify_tolerance": self.simplify_tolerance,
            "min_feature_area": self.min_feature_area,
            "target_segment_length": self.target_segment_length,
            "resample_domain_boundary": self.resample_domain_boundary,
            "resample_region_sources": self.resample_region_sources,
        }


class TriangleGrid(Triangle):
    """High-level wrapper over FloPy's Triangle for building Voronoi-ready meshes.

    Extends ``flopy.utils.triangle.Triangle`` with a declarative,
    geometry-driven workflow: declare the model domain and refinement regions
    (from shapes/lines/points), build a quality-controlled triangular mesh (via
    a :class:`MeshBuildProfile`), and produce the triangulation that
    :class:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus` dualizes into the
    Voronoi cells used for DISV/DISU grids. The triangulation is the intermediate
    step; you usually consume its Voronoi dual rather than the triangles directly.

    Lengths are in the units of the input coordinates (model length units) and
    areas in those units squared.

    Parameters
    ----------
    angle : float or None, optional
        Minimum interior angle, in degrees, that Triangle enforces by adding
        vertices (its ``-q`` switch): higher gives better-shaped triangles and
        more cells. Default ``32`` (FloPy's own default is 20). ``None`` omits
        ``-q``, so no angle constraint is imposed.
    region_point_tolerance : float or None, optional
        Minimum distance between the interior marker points chosen for
        refinement regions in :meth:`prepare`. The same number sets the gap left
        between a higher-priority region's claimed area and a lower-priority
        region's (a quarter of it), and, in :meth:`optimize_seeds`, how close a
        mesh vertex must be to a protected region's boundary to be held fixed.
        Default ``None`` scales it from the geometry: the largest of ``1e-6``,
        ``1e-6`` x the longer side of the domain's bounding box and ``1e-4`` x the
        longer side of the region's (in :meth:`optimize_seeds`, of the domain's),
        each side counted as at least ``1``; for the marker points only, also
        ``0.05 * sqrt(max_area)`` of the region.
    domain_clip_tolerance : float or None, optional
        Distance the domain is inset by before each refinement region is clipped
        to it in :meth:`prepare`, so no region segment coincides with the domain
        boundary (which makes Triangle abort). Default ``None`` uses ``1e-4`` x
        the diagonal of the domain's bounding box; ``0`` or less clips to the
        domain itself. An inset that would empty the domain is not applied.
    *args
        Positional values after ``domain_clip_tolerance``, forwarded to
        ``flopy.utils.triangle.Triangle`` as its first three parameters,
        ``model_ws``, ``exe_name`` and ``maximum_area``, in that order (so
        ``TriangleGrid("ws")`` sets ``angle``, not ``model_ws``). A fourth value
        collides with ``angle`` and raises ``TypeError``; pass these by keyword
        instead.
    **kwargs
        Forwarded to ``flopy.utils.triangle.Triangle``:

        - ``model_ws`` (str or path, default ``"."``): directory the Triangle
          input and output files are written to; :meth:`build`,
          :meth:`build_mesh` and :meth:`optimize_seeds` create it.
        - ``exe_name`` (str, default ``"triangle"``): the Triangle executable,
          resolved at construction -- ``FileNotFoundError`` if it cannot be found.
        - ``maximum_area`` (float): a global maximum triangle area. FloPy then
          passes ``-a<value>`` and never a bare ``-a``, so every per-region
          ``max_area`` (including the domain's background ``max_area``) is NOT
          applied -- measured, a 200 x 200 region at ``max_area=100`` inside a
          1000 x 1000 domain meshed to 1100 triangles without it and to 174
          with ``maximum_area=1e4``.
        - ``nodes`` (array of shape ``(n, 2)``): accepted but discarded --
          :meth:`prepare` rebuilds the node list from :meth:`add_points`, so use
          that instead.
        - ``additional_args`` (list of str): extra Triangle command-line switches.

        ``angle`` cannot be passed here; it is a named parameter.
    """

    def __init__(
        self,
        angle=32,
        region_point_tolerance: float | None = None,
        domain_clip_tolerance: float | None = None,
        *args,
        **kwargs,
    ):
        """Initialize a TriangleGrid with a minimum ``angle`` and region/domain tolerances.

        Starts with no domain or regions declared; ``None`` tolerances auto-scale
        from the geometry when regions are prepared. The parameters are
        documented on the class.
        """

        super().__init__(angle=angle, *args, **kwargs)
        self.domain_spec: DomainSpec | None = None
        self.region_specs: list[RegionSpec] = []
        self._prepared_regions: list[PreparedRegion] = []
        self._manual_regions: list[list] = []
        self._point_constraints: list[tuple[float, float]] = []
        self._optimization_points: list[tuple[float, float]] = []
        self.region_point_tolerance = region_point_tolerance
        # Refinement regions are clipped to the domain inset by this tolerance so
        # their boundary segments never coincide with the domain boundary (which
        # makes Triangle abort). None auto-scales to the domain extent.
        self.domain_clip_tolerance = domain_clip_tolerance
        self._last_cleanup_report: dict[str, float | int] | None = None
        self._last_quality_report: dict[str, float | int] | None = None
        self._last_optimization_report: dict[str, float | int | str] | None = None

    @property
    def domain_geometry(self) -> Polygon | MultiPolygon | None:
        """Return the current domain geometry, if one has been defined."""
        return None if self.domain_spec is None else self.domain_spec.geometry

    @staticmethod
    def _circle_polygon(
        radius: float = 100.0,
        center_coords: tuple = (0, 0),
        radians_step: float = 0.1,
    ) -> Polygon:
        """A circular polygon of ``radius`` about ``center_coords``, sampled every ``radians_step``."""

        theta = np.arange(0.0, 2 * np.pi, radians_step)
        x = radius * np.cos(theta) + center_coords[0]
        y = radius * np.sin(theta) + center_coords[1]
        circle_poly = [(x_coord, y_coord) for x_coord, y_coord in zip(x, y)]
        return shp.Polygon(circle_poly)

    @staticmethod
    def _rectangle_polygon(
        x_dist=100,
        y_dist=100,
        origin=(0, 0),
    ) -> Polygon:
        """An axis-aligned rectangle of ``x_dist`` by ``y_dist`` with its lower-left at ``origin``."""

        x_min, y_min = origin[0], origin[1]
        x_max, y_max = x_min + x_dist, y_min + y_dist
        polygon_coords = ((x_min, y_min), (x_min, y_max), (x_max, y_max), (x_max, y_min))
        return shp.Polygon(polygon_coords)

    @staticmethod
    def _coerce_geometry(geometry: Any):
        """Coerce a path / GeoDataFrame / shapely geometry / coord list to one shapely geometry."""

        if isinstance(geometry, Path):
            return read_shp_gpkg(geometry).union_all()
        if isinstance(geometry, gpd.GeoDataFrame | gpd.GeoSeries):
            return geometry.union_all()
        # Any shapely geometry passes through as-is. Points and lines are valid
        # inputs here -- they are turned into an area by the buffer step in
        # _polygonize_geometry (e.g. add_polygon(uic_point, buffer=4.5)).
        if isinstance(geometry, BaseGeometry):
            return geometry
        return shp.Polygon(geometry)

    @staticmethod
    def _iter_polygons(geometry: Polygon | MultiPolygon | GeometryCollection):
        """Yield each constituent ``Polygon`` of a polygon/multipolygon/collection geometry."""

        if isinstance(geometry, Polygon):
            yield geometry
            return
        if isinstance(geometry, MultiPolygon):
            for polygon in geometry.geoms:
                yield polygon
            return
        if isinstance(geometry, GeometryCollection):
            for geom in geometry.geoms:
                if isinstance(geom, Polygon):
                    yield geom
                elif isinstance(geom, MultiPolygon):
                    yield from geom.geoms
            return
        raise TypeError("geometry must be a Polygon, MultiPolygon, or GeometryCollection")

    @classmethod
    def _polygonize_geometry(
        cls,
        geometry: Any,
        *,
        domain: Polygon | MultiPolygon | Path | None = None,
        buffer: int | float = 0,
        simplify_tolerance=None,
        densify_dist: int | None = None,
        what: str = "The geometry",
    ) -> Polygon | MultiPolygon:
        """Turn arbitrary region input into a clean polygon/multipolygon for meshing.

        Coerces the geometry, optionally clips it to ``domain``, buffers points/lines
        into areas, simplifies/densifies, and drops empty pieces; raises
        :class:`NoAreaInside`, naming ``what`` and the cause, if nothing with positive area
        remains.
        """

        geometry = cls._coerce_geometry(geometry)
        drawn = geometry

        clipped = None
        if domain is not None:
            domain = cls._coerce_geometry(domain)
            geometry = clipped = geometry.intersection(domain)
            # Nothing of it inside: disjoint, or (with no buffer to widen it) touching only
            # along a line or at a point. A buffered edge-touching polygon builds a strip
            # today, so that one is left alone; an EMPTY feature (a line buffered by 0) is not
            # "outside" but empty, below. A feature outside but within ``buffer`` of the
            # boundary would reach inside once buffered, and the clip comes first: refused as
            # it was, never skipped.
            if not drawn.is_empty and (geometry.is_empty or (buffer <= 0 and drawn.area > 0 and geometry.area == 0)):
                if buffer > 0 and drawn.distance(domain) < buffer:
                    raise NoAreaInside(
                        f"{what} lies outside the model boundary, {drawn.distance(domain):.4g} from "
                        f"it -- within its buffer={buffer:g}, so its buffered area would reach inside, "
                        "but it is clipped to the boundary BEFORE it is buffered. Move it inside the "
                        "boundary, or draw it as a polygon.",
                        cause="clipped",
                    )
                raise NoAreaInside(cls._no_area_inside(what, drawn, domain, geometry), cause="outside")

        if buffer != 0:
            geometry = geometry.buffer(buffer)

        if not isinstance(geometry, (Polygon, MultiPolygon, GeometryCollection)):
            raise ValueError(
                f"A {geometry.geom_type} region has no area; pass buffer>0 to give "
                "it width, e.g. add_polygon(point_or_line, buffer=...)."
            )

        if simplify_tolerance:
            geometry = geometry.simplify(simplify_tolerance)

        polygons = [polygon.buffer(0) for polygon in cls._iter_polygons(geometry)]
        polygons = [polygon for polygon in polygons if not polygon.is_empty and polygon.area > 0]
        if densify_dist:
            polygons = [densify_poly(polygon, densify_dist) for polygon in polygons]

        if not polygons:
            source = clipped if clipped is not None else drawn
            area = [polygon for polygon in cls._iter_polygons(source) if polygon.area > 0]
            if buffer < 0 and area:
                width = 2.0 * float(shp.maximum_inscribed_circle(MultiPolygon(area)).length)
                raise NoAreaInside(
                    f"{what} has no area left after buffer={buffer:g}: its part inside the boundary "
                    f"is at most {width:.4g} wide, so the inward buffer removes all of it. Use a "
                    "smaller buffer.",
                    cause="buffer",
                )
            raise NoAreaInside(
                f"{what} has no polygonal area (area {drawn.area:.4g}, a {drawn.geom_type}); "
                "geometry with zero area cannot be meshed.",
                cause="empty",
            )

        if len(polygons) == 1:
            return polygons[0]
        return MultiPolygon(polygons)

    @staticmethod
    def _no_area_inside(what, drawn, domain, clipped) -> str:
        """Why a region clipped to the model boundary has no area: outside it, or touching it."""

        bounds = ", ".join(f"{v:.6g}" for v in drawn.bounds)
        boundary = ", ".join(f"{v:.6g}" for v in domain.bounds)
        if clipped.is_empty:
            where = f"lies entirely outside the model boundary ({drawn.distance(domain):.4g} away)"
        else:
            where = (f"touches the model boundary only along a line or at a point (the overlap "
                     f"is a {clipped.geom_type}), so no part of it is inside")
        return (f"{what} {where}. Its bounds are ({bounds}); the boundary's are ({boundary}). "
                "Check its location and CRS.")

    def _region_point_spacing(
        self,
        geometry: Polygon | None = None,
        max_area: float | int | None = None,
    ) -> float:
        """The minimum spacing between region marker points (explicit tolerance, else auto-scaled).

        Auto-scales from the domain/region extent and ``max_area`` so markers are
        distinct without crowding.
        """

        if self.region_point_tolerance is not None:
            return float(self.region_point_tolerance)

        spacing = 1e-6

        if self.domain_geometry is not None:
            xmin, ymin, xmax, ymax = self.domain_geometry.bounds
            span = max(xmax - xmin, ymax - ymin, 1.0)
            spacing = max(spacing, span * 1e-6)

        if max_area is not None and max_area > 0:
            spacing = max(spacing, float(np.sqrt(max_area)) * 0.05)

        if geometry is not None:
            xmin, ymin, xmax, ymax = geometry.bounds
            span = max(xmax - xmin, ymax - ymin, 1.0)
            spacing = max(spacing, span * 1e-4)

        return spacing

    @staticmethod
    def _point_tuple(point: Point | tuple[float, float]) -> tuple[float, float]:
        """Normalize a shapely ``Point`` or coordinate pair to an ``(x, y)`` float tuple."""

        if isinstance(point, Point):
            return (point.x, point.y)
        return (float(point[0]), float(point[1]))

    @staticmethod
    def _sample_candidate_grid(geometry: Polygon, sample_count: int = 5) -> list[tuple[float, float]]:
        """Interior points from a ``sample_count`` x ``sample_count`` grid over the geometry's bounds."""

        xmin, ymin, xmax, ymax = geometry.bounds
        xs = np.linspace(xmin, xmax, sample_count)
        ys = np.linspace(ymin, ymax, sample_count)
        candidates: list[tuple[float, float]] = []
        for x in xs:
            for y in ys:
                pt = Point(x, y)
                if geometry.covers(pt):
                    candidates.append((pt.x, pt.y))
        return candidates

    @staticmethod
    def _inward_shells(geometry: Polygon) -> list[Polygon]:
        """The geometry plus progressively inward-buffered copies, for finding a deep interior point."""

        shells = [geometry]
        if geometry.is_empty:
            return shells

        min_clearance = geometry.minimum_clearance
        if not np.isfinite(min_clearance) or min_clearance <= 0:
            return shells

        for fraction in (0.1, 0.25, 0.4):
            buffered = geometry.buffer(-(min_clearance * fraction))
            if buffered.is_empty:
                continue
            for piece in TriangleGrid._iter_polygons(buffered):
                if piece.area > 0:
                    shells.append(piece)
        return shells

    def _candidate_points(
        self,
        geometry: Polygon,
        point_hint: tuple[float, float] | None = None,
    ) -> list[tuple[float, float]]:
        """De-duplicated candidate interior points for a region marker (hint, reps, centroids, grid)."""

        candidates: list[tuple[float, float]] = []

        if point_hint is not None:
            hint_point = Point(point_hint)
            if geometry.covers(hint_point):
                candidates.append(self._point_tuple(hint_point))

        for shell in self._inward_shells(geometry):
            rep = shell.representative_point()
            candidates.append((rep.x, rep.y))

            centroid = shell.centroid
            if shell.covers(centroid):
                candidates.append((centroid.x, centroid.y))

            candidates.extend(self._sample_candidate_grid(shell, sample_count=5))

        deduped: list[tuple[float, float]] = []
        for candidate in candidates:
            if candidate not in deduped:
                deduped.append(candidate)
        return deduped

    def _choose_region_point(
        self,
        geometry: Polygon,
        *,
        existing_points: list[tuple[float, float]],
        point_hint: tuple[float, float] | None = None,
        max_area: float | int | None = None,
    ) -> tuple[float, float]:
        """A unique interior marker point for a region, at least ``spacing`` from existing markers.

        Raises if no candidate lies inside the geometry and clear of every existing point.
        """

        spacing = self._region_point_spacing(geometry=geometry, max_area=max_area)
        for candidate in self._candidate_points(geometry, point_hint=point_hint):
            point = Point(candidate)
            if not geometry.covers(point):
                continue
            if all(point.distance(Point(existing)) > spacing for existing in existing_points):
                return candidate
        raise ValueError("Could not find a unique interior point for a refinement region")

    def _claimed_region_geometry(
        self,
        polygon: Polygon,
        occupied_geometries: list[Polygon],
        *,
        label: str,
        priority: int,
    ) -> Polygon:
        """The largest part of ``polygon`` not already claimed by higher-priority regions.

        Subtracts the union of ``occupied_geometries`` (with a small buffer) and
        raises if no unique area remains.
        """

        if not occupied_geometries:
            return polygon

        spacing = self._region_point_spacing(geometry=polygon)
        excluded = shp.union_all(occupied_geometries)
        candidate = polygon.difference(excluded.buffer(spacing * 0.25))
        pieces = [piece for piece in self._iter_polygons(candidate) if not piece.is_empty and piece.area > 0]
        if not pieces:
            raise ValueError(
                f"Region {label!r} does not have any unique interior area left after resolving overlaps "
                f"at priority {priority}."
            )
        return max(pieces, key=lambda piece: piece.area)

    def _coerce_domain_region(
        self,
        geometry: Polygon | MultiPolygon,
        *,
        max_area: float | int | None,
        label: str,
    ):
        """Record the domain geometry and, when ``max_area`` is set, add a background sizing region."""

        self._refuse_a_self_touching_domain(geometry)
        self.domain_spec = DomainSpec(geometry=geometry, label=label, max_area=max_area)
        if max_area is not None:
            domain_region_geometry = self._domain_background_geometry(geometry, max_area=max_area)
            self.add_region_polygon(
                domain_region_geometry,
                max_area=max_area,
                label=f"{label}_size" if label else "domain_size",
                priority=-10,
                source="domain",
            )
            if domain_region_geometry is geometry:
                # No point of the domain is 0.5*sqrt(max_area) from its boundary, so the
                # inset is empty and the sizing region is the domain itself. Its ring, clipped
                # to the domain inset by 1e-4 of the diagonal, was a copy of the boundary
                # 0.14 ft inside it on the L, and Triangle filled that band with ~65,000-
                # 285,000 cells whatever max_area was (the L at 1e6: 149,385, smallest
                # 0.00019 ft2). Its point alone floods max_area through the whole domain.
                self.region_specs[-1].segments = False
                message = (
                    f"max_area={max_area:g} is larger than the domain '{label}' can hold: no "
                    f"point of it lies more than 0.5*sqrt(max_area) = {0.5 * np.sqrt(max_area):g} "
                    "from its boundary, so the mesh is as coarse as the boundary allows -- a "
                    "handful of cells. Lower max_area for a finer mesh."
                )
                logger.warning(message)
                warnings.warn(message, UserWarning, stacklevel=3)
        return geometry

    @classmethod
    def _refuse_a_self_touching_domain(cls, geometry: Polygon | MultiPolygon) -> None:
        """Refuse a domain that touches itself: two parts meeting, a hole meeting an edge.

        Shapely calls both valid, and Triangle meshes them, but the node where the domain
        touches itself has two separate fans of triangles, so its Voronoi cell would be two
        pieces -- which one DISV cell cannot be. Refused here, where the place can be named,
        rather than as FloPy's "duplicate points" or a fan that does not close.
        """

        parts = list(cls._iter_polygons(geometry))
        rings = []
        for index, polygon in enumerate(parts):
            rings.append((f"part {index}'s edge", polygon.exterior))
            rings.extend((f"part {index}'s hole {j}", ring) for j, ring in enumerate(polygon.interiors))
        for i, (name_a, ring_a) in enumerate(rings):
            for name_b, ring_b in rings[i + 1:]:
                if ring_a.intersects(ring_b):
                    x, y = ring_a.intersection(ring_b).representative_point().coords[0]
                    raise ValueError(
                        f"The model boundary touches itself at ({x:.6g}, {y:.6g}): {name_a} "
                        f"meets {name_b}. A Voronoi cell there would be two pieces, which a "
                        "DISV cell cannot be. Separate them, or merge them into one polygon."
                    )

    def _domain_background_geometry(
        self,
        geometry: Polygon | MultiPolygon,
        *,
        max_area: float | int,
    ) -> Polygon | MultiPolygon:
        """An inset copy of the domain used as the background sizing region.

        Buffered inward so its marker point stays clear of the domain boundary
        (where a coincident region point would make Triangle abort).
        """

        polygons = list(self._iter_polygons(geometry))
        if not polygons:
            return geometry

        xmin, ymin, xmax, ymax = geometry.bounds
        span = max(xmax - xmin, ymax - ymin, 1.0)
        inward = max(span * 1e-4, float(np.sqrt(max_area)) * 0.5)

        inset_polygons: list[Polygon] = []
        for polygon in polygons:
            inset = polygon.buffer(-inward)
            if inset.is_empty:
                continue
            inset_polygons.extend(
                piece.buffer(0)
                for piece in self._iter_polygons(inset)
                if not piece.is_empty and piece.area > 0
            )

        if not inset_polygons:
            return geometry
        if len(inset_polygons) == 1:
            return inset_polygons[0]
        return MultiPolygon(inset_polygons)

    def set_domain_polygon(
        self,
        polygon: Polygon | MultiPolygon | Path,
        *,
        label: str = "domain",
        max_area: float | int | None = None,
        buffer: int | float = 0,
        simplify_tolerance=None,
        densify_dist: int | None = None,
    ):
        """Define the model domain from polygonal geometry.

        The geometry is buffered, simplified and densified, in that order, and
        stored as :attr:`domain_spec`, replacing any earlier domain. A
        multi-part geometry gives a multi-part domain; interior rings become
        holes in the mesh.

        Parameters
        ----------
        polygon : shapely geometry, pathlib.Path, GeoDataFrame, GeoSeries or sequence
            The domain. A ``Path`` names a ``.shp`` or ``.gpkg`` file (a
            GeoPackage's first layer) whose features are unioned; a
            ``GeoDataFrame``/``GeoSeries`` is unioned; a sequence of ``(x, y)``
            pairs is read as one polygon's exterior ring. Coordinates are used as
            they are, with no reprojection. A ``str`` is not read as a path --
            wrap it in ``pathlib.Path``. A point or line needs ``buffer > 0``.
        label : str, optional
            Name of the domain, kept as ``domain_spec.label``. When ``max_area``
            is given it also names the background sizing region
            ``"<label>_size"`` (``"domain_size"`` for an empty label), the label
            :meth:`preview_regions` shows. Default ``"domain"``.
        max_area : float, optional
            Background maximum triangle area for the whole domain, in model
            length units squared. When given, a copy of the domain inset by
            ``max(0.5 * sqrt(max_area), 1e-4 x its longer bounding-box side)``
            (the domain itself if that inset would empty it) is added as a
            refinement region with priority ``-10`` and source ``"domain"``, so
            every region of higher priority (the default ``0`` included) claims
            its area first. Default ``None`` adds no background region: triangles
            outside every refinement region then have no area limit, only the
            ``angle`` constraint.
        buffer : float, optional
            Distance to buffer the geometry by (``shapely`` ``buffer``; negative
            shrinks it). Default ``0``, no buffer.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance, topology preserved, applied
            after the buffer. Default ``None``; ``None`` or ``0`` skips it.
        densify_dist : float, optional
            Resample each polygon's exterior to points this far apart, applied
            last. Vertices between the samples (corners included) are not kept,
            and interior rings (holes) are dropped. Default ``None``; ``None`` or
            ``0`` skips it.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The stored domain geometry.

        Raises
        ------
        ValueError
            If no polygonal area remains after preprocessing (a point or line
            with no ``buffer``, or a buffer that removes everything), or if the
            boundary touches itself (two parts meeting, or a hole meeting an
            edge), which a Voronoi cell cannot represent.
        """
        geometry = self._polygonize_geometry(
            polygon,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
            what=f"The model boundary '{label}'",
        )
        return self._coerce_domain_region(geometry, max_area=max_area, label=label)

    def set_domain_rectangle(
        self,
        x_dist=100,
        y_dist=100,
        origin=(0, 0),
        *,
        label: str = "domain",
        max_area: float | int | None = None,
    ):
        """Define the domain as an axis-aligned rectangle.

        Parameters
        ----------
        x_dist : float, optional
            Width of the rectangle along x, in model length units. Default ``100``.
        y_dist : float, optional
            Height of the rectangle along y, in model length units. Default ``100``.
        origin : tuple of float, optional
            ``(x, y)`` of the rectangle's lower-left corner. Default ``(0, 0)``.
        label : str, optional
            Name of the domain; as in :meth:`set_domain_polygon`. Default
            ``"domain"``.
        max_area : float, optional
            Background maximum triangle area for the whole domain; as in
            :meth:`set_domain_polygon`. Default ``None``, no background region.

        Returns
        -------
        shapely.geometry.Polygon
            The stored domain rectangle.
        """
        polygon = self._rectangle_polygon(x_dist=x_dist, y_dist=y_dist, origin=origin)
        return self._coerce_domain_region(polygon, max_area=max_area, label=label)

    def set_domain_circle(
        self,
        radius: float = 100.0,
        center_coords: tuple = (0, 0),
        *,
        radians_step: float = 0.1,
        label: str = "domain",
        max_area: float | int | None = None,
    ):
        """Define the domain as an approximate circular polygon.

        Parameters
        ----------
        radius : float, optional
            Circle radius, in model length units. Default ``100.0``.
        center_coords : tuple of float, optional
            ``(x, y)`` of the centre. Default ``(0, 0)``.
        radians_step : float, optional
            Angle, in radians, between successive vertices: the polygon has a
            vertex at every multiple of it below ``2 * pi`` (63 at the default).
            Default ``0.1``.
        label : str, optional
            Name of the domain; as in :meth:`set_domain_polygon`. Default
            ``"domain"``.
        max_area : float, optional
            Background maximum triangle area for the whole domain; as in
            :meth:`set_domain_polygon`. Default ``None``, no background region.

        Returns
        -------
        shapely.geometry.Polygon
            The stored domain polygon.
        """
        polygon = self._circle_polygon(radius=radius, center_coords=center_coords, radians_step=radians_step)
        return self._coerce_domain_region(polygon, max_area=max_area, label=label)

    def set_domain_file(
        self,
        shp_gpkg: Path,
        *,
        label: str = "domain",
        max_area: float | int | None = None,
        buffer: int | float = 0,
        simplify_tolerance=None,
        densify_dist: int | None = None,
    ):
        """Load the model domain from GIS data and register it.

        The same as :meth:`set_domain_polygon`, which it calls with every
        argument unchanged.

        Parameters
        ----------
        shp_gpkg : pathlib.Path
            A ``.shp`` or ``.gpkg`` file (a GeoPackage's first layer) whose
            features are unioned into the domain. Any other input
            :meth:`set_domain_polygon` takes works too; a ``str`` path does not.
        label : str, optional
            Name of the domain; as in :meth:`set_domain_polygon`. Default
            ``"domain"``.
        max_area : float, optional
            Background maximum triangle area for the whole domain; as in
            :meth:`set_domain_polygon`. Default ``None``, no background region.
        buffer : float, optional
            Distance to buffer the geometry by; negative shrinks it. Default ``0``.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance, topology preserved. Default
            ``None``, no simplification.
        densify_dist : float, optional
            Resample each exterior to points this far apart; corners between the
            samples and holes are lost. Default ``None``, no resampling.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The stored domain geometry.

        Raises
        ------
        ValueError
            As :meth:`set_domain_polygon`.
        """
        return self.set_domain_polygon(
            shp_gpkg,
            label=label,
            max_area=max_area,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
        )

    def add_region_polygon(
        self,
        polygon: Polygon | MultiPolygon | Path,
        *,
        max_area: float | int,
        label: str | None = None,
        priority: int = 0,
        point_hint: tuple[float, float] | None = None,
        source: str = "geometry",
        domain: Polygon | MultiPolygon | Path | None = None,
        buffer: int | float = 0,
        simplify_tolerance=None,
        densify_dist: int | None = None,
    ):
        """Add a polygonal refinement region.

        The geometry is clipped to ``domain`` (or the stored domain), then
        buffered, simplified and densified, and recorded in
        :attr:`region_specs`. Overlaps with other regions are resolved later, in
        :meth:`prepare`, which also clips every region to the domain inset by
        ``domain_clip_tolerance``.

        Parameters
        ----------
        polygon : shapely geometry, pathlib.Path, GeoDataFrame, GeoSeries or sequence
            The region, in any form :meth:`set_domain_polygon` takes (a ``Path``
            to a ``.shp``/``.gpkg`` whose features are unioned, a shapely
            geometry, a GeoDataFrame/GeoSeries, or a sequence of ``(x, y)`` read
            as a polygon's exterior ring; not a ``str``). A point or line needs
            ``buffer > 0``.
        max_area : float
            Maximum triangle area requested inside the region, in model length
            units squared. Required.
        label : str, optional
            Region name, shown by :meth:`preview_regions` and matched by
            ``protected_labels`` in :meth:`optimize_seeds`/:meth:`build_mesh`.
            Default ``None``: :meth:`prepare` names it ``"region_<i>"``, ``i``
            being its position in priority order. A region that ends up in
            several parts gets ``"_<j>"`` appended per part.
        priority : int, optional
            Higher-priority regions claim overlapping area before lower
            priority regions during :meth:`prepare`; at equal priority the
            smaller region claims first. A region left with no area of its own
            raises ``ValueError`` there. Default ``0``.
        point_hint : tuple of float, optional
            Optional interior point to try first when selecting Triangle's
            region point. Used only if it lies inside the area the region claims
            and clear of the other regions' points; otherwise a point is chosen
            automatically, without a warning. Default ``None``.
        source : str, optional
            Short source label used by cleanup and feature-protection logic:
            :meth:`clean_geometry`'s ``resample_region_sources`` and
            :meth:`build_mesh`'s ``protect_sources`` select regions by it, and
            :meth:`preview_regions` shows it. The helpers set ``"line"``,
            ``"rectangle"``, ``"circle"``, ``"file"``, ``"polygon"`` and others.
            Default ``"geometry"``.
        domain : shapely geometry, pathlib.Path, GeoDataFrame, GeoSeries or sequence, optional
            Geometry to clip the region to instead of the stored domain, in the
            same forms as ``polygon``. Default ``None`` clips to the stored
            domain, or not at all when none is set yet.
        buffer : float, optional
            Distance to buffer the region by AFTER it is clipped (``shapely``
            ``buffer``); a negative value pulls it back from the clip edge.
            Default ``0``, no buffer.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance, topology preserved, applied
            after the buffer. Default ``None``; ``None`` or ``0`` skips it.
        densify_dist : float, optional
            Resample each polygon's exterior to points this far apart, applied
            last. Vertices between the samples (corners included) are not kept,
            and interior rings are dropped. Default ``None``; ``None`` or ``0``
            skips it.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The region geometry as recorded.

        Raises
        ------
        ValueError
            If no polygonal area remains after clipping and preprocessing.
        """
        clip_domain = self.domain_geometry if domain is None else domain
        geometry = self._polygonize_geometry(
            polygon,
            domain=clip_domain,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
            what=f"{'Breakline' if source == 'line' else 'Refinement region'} '{label or 'unnamed'}'",
        )
        spec = RegionSpec(
            geometry=geometry,
            max_area=max_area,
            label=label,
            priority=priority,
            point_hint=point_hint,
            source=source,
        )
        self.region_specs.append(spec)
        return geometry

    def add_region_rectangle(
        self,
        *,
        x_dist=100,
        y_dist=100,
        origin=(0, 0),
        max_area: float | int,
        label: str | None = None,
        priority: int = 0,
        point_hint: tuple[float, float] | None = None,
    ):
        """Add a rectangular refinement region.

        Calls :meth:`add_region_polygon` with source ``"rectangle"``, so the
        rectangle is clipped to the stored domain.

        Parameters
        ----------
        x_dist : float, optional
            Width of the rectangle along x, in model length units. Default ``100``.
        y_dist : float, optional
            Height of the rectangle along y, in model length units. Default ``100``.
        origin : tuple of float, optional
            ``(x, y)`` of the rectangle's lower-left corner. Default ``(0, 0)``.
        max_area : float
            Maximum triangle area inside the region, in model length units
            squared. Required.
        label : str, optional
            Region name; as in :meth:`add_region_polygon`. Default ``None``,
            named ``"region_<i>"`` by :meth:`prepare`.
        priority : int, optional
            Overlap priority; as in :meth:`add_region_polygon`. Default ``0``.
        point_hint : tuple of float, optional
            ``(x, y)`` to try first as Triangle's region point; as in
            :meth:`add_region_polygon`. Default ``None``.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The region geometry as recorded (the rectangle clipped to the domain).
        """
        polygon = self._rectangle_polygon(x_dist=x_dist, y_dist=y_dist, origin=origin)
        return self.add_region_polygon(
            polygon,
            max_area=max_area,
            label=label,
            priority=priority,
            point_hint=point_hint,
            source="rectangle",
        )

    def add_region_circle(
        self,
        *,
        radius: float = 100.0,
        center_coords: tuple = (0, 0),
        max_area: float | int,
        label: str | None = None,
        priority: int = 0,
        point_hint: tuple[float, float] | None = None,
        radians_step: float = 0.1,
    ):
        """Add a circular refinement region.

        Calls :meth:`add_region_polygon` with source ``"circle"``, so the circle
        is clipped to the stored domain.

        Parameters
        ----------
        radius : float, optional
            Circle radius, in model length units. Default ``100.0``.
        center_coords : tuple of float, optional
            ``(x, y)`` of the centre. Default ``(0, 0)``.
        max_area : float
            Maximum triangle area inside the region, in model length units
            squared. Required.
        label : str, optional
            Region name; as in :meth:`add_region_polygon`. Default ``None``,
            named ``"region_<i>"`` by :meth:`prepare`.
        priority : int, optional
            Overlap priority; as in :meth:`add_region_polygon`. Default ``0``.
        point_hint : tuple of float, optional
            ``(x, y)`` to try first as Triangle's region point; as in
            :meth:`add_region_polygon`. Default ``None``.
        radians_step : float, optional
            Angle, in radians, between successive vertices of the circle
            polygon (a vertex at every multiple below ``2 * pi``). Default ``0.1``.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The region geometry as recorded (the circle clipped to the domain).
        """
        polygon = self._circle_polygon(radius=radius, center_coords=center_coords, radians_step=radians_step)
        return self.add_region_polygon(
            polygon,
            max_area=max_area,
            label=label,
            priority=priority,
            point_hint=point_hint,
            source="circle",
        )

    def add_region_file(
        self,
        shp_gpkg: Path,
        *,
        max_area: float | int,
        label: str | None = None,
        priority: int = 0,
        point_hint: tuple[float, float] | None = None,
        explode: bool = False,
        buffer: int | float = 0,
        simplify_tolerance=None,
        densify_dist: int | None = None,
    ):
        """Load one or more refinement regions from GIS data.

        The features are unioned, clipped to the stored domain, buffered,
        simplified and densified (as in :meth:`add_region_polygon`), and
        recorded with source ``"file"``.

        Parameters
        ----------
        shp_gpkg : pathlib.Path
            A ``.shp`` or ``.gpkg`` file (a GeoPackage's first layer) whose
            features are unioned. Any other input :meth:`add_region_polygon`
            takes works too (a GeoDataFrame, a shapely geometry, ...); a ``str``
            path and a list of paths do not.
        max_area : float
            Maximum triangle area inside the region(s), in model length units
            squared. Required.
        label : str, optional
            Region name; as in :meth:`add_region_polygon`. With ``explode=True``
            each part is named ``"<label>_<i>"``. Default ``None``, named
            ``"region_<i>"`` by :meth:`prepare`.
        priority : int, optional
            Overlap priority; as in :meth:`add_region_polygon`. Default ``0``.
        point_hint : tuple of float, optional
            ``(x, y)`` to try first as Triangle's region point; as in
            :meth:`add_region_polygon`. With ``explode=True`` every part is given
            the same hint. Default ``None``.
        explode : bool, optional
            ``True`` records each polygon part as its own region. Default
            ``False`` records the whole geometry as one region, which
            :meth:`prepare` still splits into one Triangle region per part.
        buffer : float, optional
            Distance to buffer by after clipping; negative shrinks. Default ``0``.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance, topology preserved. Default
            ``None``, no simplification.
        densify_dist : float, optional
            Resample each exterior to points this far apart; corners between the
            samples and holes are lost. Default ``None``, no resampling.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The region geometry after preprocessing.

        Raises
        ------
        ValueError
            If no polygonal area remains after clipping and preprocessing.
        """
        geometry = self._polygonize_geometry(
            shp_gpkg,
            domain=self.domain_geometry,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
            what=f"Refinement file {shp_gpkg}",
        )
        if explode:
            for idx, polygon in enumerate(self._iter_polygons(geometry)):
                region_label = None if label is None else f"{label}_{idx}"
                self.region_specs.append(
                    RegionSpec(
                        geometry=polygon,
                        max_area=max_area,
                        label=region_label,
                        priority=priority,
                        point_hint=point_hint,
                        source="file",
                    )
                )
            return geometry

        return self.add_region_polygon(
            geometry,
            max_area=max_area,
            label=label,
            priority=priority,
            point_hint=point_hint,
            source="file",
        )

    def add_line_feature(
        self,
        line: Path | Any,
        *,
        buffer: int | float = 10,
        max_area: float | int,
        label: str | None = None,
        priority: int = 0,
        point_hint: tuple[float, float] | None = None,
        simplify_tolerance: int | float | None = 10,
        densify_dist: int | None = None,
        negative_buffer_after_clipping: float | int = 0,
    ):
        """Add refinement around a linear feature by buffering it into polygons.

        The line is buffered by ``buffer``, then handed to
        :meth:`add_region_polygon` with source ``"line"``: clipped to the stored
        domain, buffered again by ``negative_buffer_after_clipping``, simplified
        and densified. Source ``"line"`` is what :meth:`build_mesh` protects and
        :meth:`clean_geometry` resamples by default.

        Parameters
        ----------
        line : shapely geometry, pathlib.Path, GeoDataFrame or GeoSeries
            The feature: a ``LineString``/``MultiLineString`` (any shapely
            geometry is buffered as it is), a ``Path`` to a ``.shp``/``.gpkg``
            whose features are unioned, or a GeoDataFrame/GeoSeries, unioned. A
            bare sequence of ``(x, y)`` pairs is read as a POLYGON ring, not a
            line, so the buffer surrounds the filled shape -- pass
            ``shapely.LineString(coords)`` instead. Not a ``str`` path.
        buffer : float, optional
            Half-width of the refinement band around the line, in model length
            units. Default ``10``.
        max_area : float
            Maximum triangle area inside the band, in model length units squared.
            Required.
        label : str, optional
            Region name; as in :meth:`add_region_polygon`. Default ``None``,
            named ``"region_<i>"`` by :meth:`prepare`.
        priority : int, optional
            Overlap priority; as in :meth:`add_region_polygon`. Default ``0``.
        point_hint : tuple of float, optional
            ``(x, y)`` to try first as Triangle's region point; as in
            :meth:`add_region_polygon`. Default ``None``.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance for the band's outline,
            topology preserved, applied after clipping. Default ``10``; ``None``
            or ``0`` skips it.
        densify_dist : float, optional
            Resample the band's exterior to points this far apart; corners
            between the samples and holes are lost. Default ``None``, no
            resampling.
        negative_buffer_after_clipping : float, optional
            A buffer of zero or less applied after the band is clipped to the
            domain, pulling it back from the domain boundary. Default ``0``, none.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The band as recorded.

        Raises
        ------
        ValueError
            If ``negative_buffer_after_clipping`` is greater than zero, or if no
            area remains after clipping and preprocessing.
        """
        if negative_buffer_after_clipping > 0:
            raise ValueError("negative_buffer_after_clipping must be less than or equal to zero")
        line_geom = self._coerce_geometry(line)
        region_geom = line_geom.buffer(buffer)
        return self.add_region_polygon(
            region_geom,
            max_area=max_area,
            label=label,
            priority=priority,
            point_hint=point_hint,
            source="line",
            buffer=negative_buffer_after_clipping,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
        )

    def add_circle(
        self,
        radius: float = 100.0,
        center_coords: tuple = (0, 0),
        polygon_to_add=None,
        point_to_add=None,
        point_region_size_max=1,
        return_only: bool = False,
        radians_step: int = 0.1,
    ):
        """Backward-compatible helper for adding a circular polygon feature.

        Prefer :meth:`set_domain_circle` / :meth:`add_region_circle`. With no
        domain set, the circle becomes the domain (no background ``max_area``).
        Otherwise it goes through :meth:`add_polygon`, which records it as a
        refinement region with ``max_area=1`` -- one square model length unit,
        which can produce a very large mesh.

        Parameters
        ----------
        radius : float, optional
            Circle radius, in model length units. Default ``100.0``.
        center_coords : tuple of float, optional
            ``(x, y)`` of the centre. Default ``(0, 0)``.
        polygon_to_add : sequence of (x, y), optional
            Exterior ring of an extra polygon, added after the circle through
            :meth:`add_polygon` -- so as a refinement region with ``max_area=1``.
            Default ``None``, none.
        point_to_add : tuple of float, optional
            ``(x, y)`` of a raw Triangle region point, added through
            :meth:`add_region` with ``point_region_size_max`` as its area.
            Default ``None``, none.
        point_region_size_max : float, optional
            Maximum triangle area for the ``point_to_add`` region, in model
            length units squared. Ignored without ``point_to_add``. Default ``1``.
        return_only : bool, optional
            ``True`` returns the circle polygon without adding anything, and
            ignores ``polygon_to_add``/``point_to_add``. Default ``False``.
        radians_step : float, optional
            Angle, in radians, between successive vertices of the circle
            polygon. Default ``0.1``.

        Returns
        -------
        shapely.geometry.Polygon
            The circle polygon, unclipped.
        """
        circle_poly = self._circle_polygon(
            radius=radius,
            center_coords=center_coords,
            radians_step=radians_step,
        )
        if return_only:
            return circle_poly

        if self.domain_spec is None and polygon_to_add is None and point_to_add is None:
            self.set_domain_polygon(circle_poly)
        else:
            self.add_polygon(circle_poly)

        if polygon_to_add:
            self.add_polygon(polygon=Polygon(shell=polygon_to_add))
        if point_to_add:
            self.add_region(point=point_to_add, maximum_area=point_region_size_max)

        return circle_poly

    def add_rectangle(
        self,
        x_dist=100,
        y_dist=100,
        origin=(0, 0),
        return_only=False,
        max_area=None,
    ):
        """Backward-compatible helper for adding a rectangular polygon feature.

        Prefer :meth:`set_domain_rectangle` / :meth:`add_region_rectangle`. With
        no domain set, the rectangle becomes the domain; otherwise it is recorded
        as a refinement region with source ``"rectangle"``.

        Parameters
        ----------
        x_dist : float, optional
            Width of the rectangle along x, in model length units. Default ``100``.
        y_dist : float, optional
            Height of the rectangle along y, in model length units. Default ``100``.
        origin : tuple of float, optional
            ``(x, y)`` of the rectangle's lower-left corner. Default ``(0, 0)``.
        return_only : bool, optional
            ``True`` returns the rectangle without adding anything. Default
            ``False``.
        max_area : float, optional
            As the domain, its background maximum triangle area (``None``: no
            background region, see :meth:`set_domain_polygon`). As a region, its
            maximum triangle area, with ``None`` meaning ``1`` square model
            length unit. Default ``None``.

        Returns
        -------
        shapely.geometry.Polygon
            The rectangle, unclipped.
        """
        polygon = self._rectangle_polygon(x_dist=x_dist, y_dist=y_dist, origin=origin)
        if return_only:
            return polygon

        if self.domain_spec is None:
            self.set_domain_polygon(polygon, max_area=max_area)
        else:
            self.add_region_polygon(polygon, max_area=max_area if max_area is not None else 1, source="rectangle")

        return polygon

    def add_region(self, point, attribute=0, maximum_area=None):
        """Add a raw Triangle region point and optional area constraint.

        Unlike :meth:`add_region_polygon`, this adds no boundary segments and
        takes no part in :meth:`prepare`'s overlap resolution: the point is
        passed to Triangle as it is, after the prepared regions, and Triangle
        applies it to whichever segment-bounded area contains it.

        Parameters
        ----------
        point : tuple of float or shapely.geometry.Point
            ``(x, y)`` of the region point (anything FloPy's ``GeoSpatialUtil``
            reads as a point).
        attribute : int or float, optional
            Regional attribute Triangle assigns to the triangles of that area.
            Default ``0``.
        maximum_area : float, optional
            Maximum triangle area for that area, in model length units squared.
            Default ``None``, no area constraint from this point.
        """
        point = GeoSpatialUtil(point, shapetype="point").points
        self._manual_regions.append([point, attribute, maximum_area])

    def add_regions(self, points, attributes=None, maximum_areas=None):
        """Add multiple raw Triangle region points in one call.

        Each point is recorded as :meth:`add_region` records one. The three
        sequences are paired by position and stop at the shortest.

        Parameters
        ----------
        points : sequence of tuple of float or shapely.geometry.Point
            ``(x, y)`` of each region point.
        attributes : sequence of int or float, optional
            Regional attribute per point. Default ``None``, ``0`` for every point.
        maximum_areas : sequence of float or None, optional
            Maximum triangle area per point, in model length units squared.
            Default ``None``, no area constraint for any point.
        """
        if attributes is None:
            attributes = [0 for _ in points]
        if maximum_areas is None:
            maximum_areas = [None for _ in points]

        for region in zip(points, attributes, maximum_areas):
            self._manual_regions.append(list(region))
        return

    def add_poly_regions(
        self,
        shp_gpkg: list | Path,
        points: list | tuple = None,
        use_representative_point: bool = True,
        maximum_areas: list | int | float = None,
        *args,
        **kwargs,
    ):
        """Legacy helper that adds one or more polygon regions from GIS data.

        Prefer :meth:`add_region_file`, which this calls once per file.

        Parameters
        ----------
        shp_gpkg : pathlib.Path or list of pathlib.Path
            One ``.shp``/``.gpkg`` file, or a list of them, each becoming its
            own region (see :meth:`add_region_file`). Not a ``str``.
        points : tuple of float or list of tuple of float, optional
            Region-point hints: one ``(x, y)`` for a single file, or one per file
            for a list. Used only when ``use_representative_point=False``.
            Default ``None``.
        use_representative_point : bool, optional
            ``True`` ignores ``points`` and lets the region point be chosen
            automatically; ``False`` passes ``points`` as ``point_hint``.
            Default ``True``.
        maximum_areas : float or list of float, optional
            Maximum triangle area, in model length units squared: one value for
            every file, or a list with one per file. Effectively required: the
            default ``None`` records regions with no area, which
            :meth:`prepare` refuses with ``ValueError``.
        *args
            Forwarded positionally to :meth:`add_region_file` after the path.
            That method takes no further positional arguments, so any value here
            raises ``TypeError``.
        **kwargs
            Forwarded to :meth:`add_region_file` for every file: ``label``,
            ``priority``, ``explode``, ``buffer``, ``simplify_tolerance``,
            ``densify_dist``. ``max_area`` and ``point_hint`` are set from
            ``maximum_areas``/``points`` and raise ``TypeError`` if given here.
        """
        if isinstance(shp_gpkg, list):
            for idx, path in enumerate(shp_gpkg):
                area = maximum_areas[idx] if isinstance(maximum_areas, list) else maximum_areas
                point_hint = points[idx] if points is not None else None
                self.add_region_file(
                    path,
                    max_area=area,
                    point_hint=point_hint if not use_representative_point else None,
                    *args,
                    **kwargs,
                )
            return

        self.add_region_file(
            shp_gpkg,
            max_area=maximum_areas,
            point_hint=None if use_representative_point else points,
            *args,
            **kwargs,
        )

    def add_polygon(
        self,
        polygon: Polygon | Path,
        domain: Polygon | Path = None,
        buffer: int | float = 0,
        simplify_tolerance=None,
        ignore_holes=True,
        max_area=None,
        densify_dist: int = None,
    ):
        """Backward-compatible polygon entry point.

        If no domain exists yet, the polygon becomes the domain. Otherwise it
        is added as a refinement region.

        Prefer :meth:`set_domain_polygon` / :meth:`add_region_polygon`, which
        this calls.

        Parameters
        ----------
        polygon : shapely geometry, pathlib.Path, GeoDataFrame, GeoSeries or sequence
            The geometry, in any form :meth:`set_domain_polygon` takes.
        domain : shapely geometry, pathlib.Path, GeoDataFrame, GeoSeries or sequence, optional
            Geometry to clip the region to. Giving it always adds a refinement
            region, even when no domain is set (this does not set one). Default
            ``None``: the polygon becomes the domain if none is set, otherwise a
            region clipped to the stored domain.
        buffer : float, optional
            Distance to buffer the geometry by (after clipping, for a region);
            negative shrinks. Default ``0``.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance, topology preserved. Default
            ``None``, no simplification.
        ignore_holes : bool, optional
            Ignored; kept for backward compatibility. Holes in a domain polygon
            are always meshed as holes, and a region's holes are kept as
            segments. Default ``True``.
        max_area : float, optional
            As the domain, its background maximum triangle area (``None``: no
            background region, see :meth:`set_domain_polygon`). As a region, its
            maximum triangle area, with ``None`` meaning ``1`` square model
            length unit. Default ``None``.
        densify_dist : float, optional
            Resample each exterior to points this far apart; corners between the
            samples and holes are lost. Default ``None``, no resampling.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The stored domain, or the region geometry as recorded (source
            ``"polygon"``).

        Raises
        ------
        ValueError
            As :meth:`set_domain_polygon` or :meth:`add_region_polygon`.
        """
        del ignore_holes
        if self.domain_spec is None and domain is None:
            return self.set_domain_polygon(
                polygon,
                max_area=max_area,
                buffer=buffer,
                simplify_tolerance=simplify_tolerance,
                densify_dist=densify_dist,
            )

        if max_area is None:
            max_area = 1

        return self.add_region_polygon(
            polygon,
            domain=domain,
            max_area=max_area,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
            source="polygon",
        )

    def add_points(self, points: shp.MultiPoint):
        """Add fixed point constraints to the Triangle node list.

        Each point becomes a mesh vertex, and so the seed of a Voronoi cell.
        Calls accumulate. At :meth:`prepare`, points that repeat one another or a
        domain/region vertex (compared at 8 decimal places) are dropped.

        Parameters
        ----------
        points : shapely.geometry.MultiPoint or sequence
            The points: a ``MultiPoint``, or anything ``shapely.MultiPoint``
            accepts -- a sequence of ``(x, y)`` pairs or of ``Point`` objects.

        Raises
        ------
        TypeError
            If ``points`` cannot be converted to a ``shapely.MultiPoint``
            (shapely's own ``ValueError`` for malformed coordinates passes
            through unchanged).
        """
        if not isinstance(points, shp.MultiPoint):
            try:
                points = shp.MultiPoint(points)
            except TypeError as exc:
                raise TypeError("points must be convertible to shapely.MultiPoint") from exc

        self._point_constraints.extend(
            (float(point.x), float(point.y))
            for point in points.geoms
        )

    def clean_geometry(
        self,
        *,
        snap_tolerance: float = 0.0,
        simplify_tolerance: float | None = None,
        min_feature_area: float = 0.0,
        target_segment_length: float | None = None,
        resample_domain_boundary: bool = False,
        resample_region_sources: tuple[str, ...] = ("line",),
    ):
        """Clean stored domain and region geometry before meshing.

        Rewrites :attr:`domain_spec` and every entry of :attr:`region_specs` in
        place. Each geometry is repaired (``buffer(0)``), snapped, simplified,
        stripped of small parts and, where selected, resampled, in that order.
        Lengths are in model length units, areas in their square.

        Parameters
        ----------
        snap_tolerance : float, optional
            Distance within which each geometry's vertices and segments are
            snapped to its own vertices (``shapely.snap`` of the geometry to
            itself). Default ``0.0``; ``0`` or less disables snapping.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance, topology preserved, for the
            domain and every region. Default ``None``; ``None`` or ``0`` skips it.
        min_feature_area : float, optional
            Polygon parts smaller than this area are dropped from the domain and
            from each region -- unless every part is smaller, in which case all
            are kept. Default ``0.0`` (keep everything).
        target_segment_length : float, optional
            Spacing to resample boundaries to: the ring is rebuilt from points
            this far apart along it, so vertices between them (corners included)
            are not kept, and interior rings (holes) are dropped. Applies to the
            domain only with ``resample_domain_boundary=True`` and to the regions
            whose source is in ``resample_region_sources``. Default ``None``;
            ``None`` or ``0`` disables resampling.
        resample_domain_boundary : bool, optional
            Resample the domain boundary too. Default ``False``. Resampling drops
            the domain's holes.
        resample_region_sources : tuple of str, optional
            Region ``source`` tags (see :meth:`add_region_polygon`) whose
            boundaries are resampled. Default ``("line",)``, the bands
            :meth:`add_line_feature` adds.

        Returns
        -------
        dict
            Small report describing the cleanup settings plus the vertex and
            polygon counts before and after cleanup.

        Raises
        ------
        ValueError
            If no domain is set, a region has no ``max_area``, or cleanup removes
            all polygonal area from a geometry.
        """
        self.validate_setup()

        domain_geometry, domain_stats = cleanup_polygonal_geometry(
            self.domain_spec.geometry,
            snap_tolerance=snap_tolerance,
            simplify_tolerance=simplify_tolerance,
            min_feature_area=min_feature_area,
            target_segment_length=target_segment_length if resample_domain_boundary else None,
            keep_collapsing_parts=True,
        )
        self.domain_spec = DomainSpec(
            geometry=domain_geometry,
            label=self.domain_spec.label,
            max_area=self.domain_spec.max_area,
            attribute=self.domain_spec.attribute,
        )

        cleaned_specs: list[RegionSpec] = []
        vertices_before = domain_stats.vertices_in
        vertices_after = domain_stats.vertices_out
        polygons_before = domain_stats.polygons_in
        polygons_after = domain_stats.polygons_out

        for region in self.region_specs:
            resample = (
                target_segment_length
                if target_segment_length and region.source in set(resample_region_sources)
                else None
            )
            geometry, stats = cleanup_polygonal_geometry(
                region.geometry,
                snap_tolerance=snap_tolerance,
                simplify_tolerance=simplify_tolerance,
                min_feature_area=min_feature_area,
                target_segment_length=resample,
            )
            if resample and stats.vanished:
                # Resampling rebuilds an outline from points `resample` apart, so one shorter
                # than twice that comes back with no area -- a breakline of ~35-95 ft at the
                # default 10 ft buffer and 120 ft spacing, or a piece of a longer one the
                # boundary cut off. prepare() skipped a whole such region saying it was "empty
                # after clipping to the domain interior", and a piece went without a word.
                gone = not self._polygonal_parts(geometry)
                lengths = ", ".join(f"{length:.3g}" for length in stats.vanished)
                message = (
                    f"Region '{region.label}' ({region.source}): "
                    + ("it was removed" if gone else f"{len(stats.vanished)} of its piece(s) were removed")
                    + f" when resampled at target_segment_length={resample:g} -- an outline "
                    f"shorter than twice that spacing ({lengths} long) has no area left. Pass a "
                    "smaller target_segment_length, resample_region_sources=(), or "
                    "profile='baseline' to keep it."
                )
                logger.warning(message)
                warnings.warn(message, UserWarning, stacklevel=2)
                if gone:
                    continue
            cleaned_specs.append(
                RegionSpec(
                    geometry=geometry,
                    max_area=region.max_area,
                    label=region.label,
                    priority=region.priority,
                    attribute=region.attribute,
                    point_hint=region.point_hint,
                    source=region.source,
                    segments=region.segments,
                )
            )
            vertices_before += stats.vertices_in
            vertices_after += stats.vertices_out
            polygons_before += stats.polygons_in
            polygons_after += stats.polygons_out

        self.region_specs = cleaned_specs
        self._last_cleanup_report = {
            "snap_tolerance": float(snap_tolerance),
            "simplify_tolerance": None if simplify_tolerance is None else float(simplify_tolerance),
            "min_feature_area": float(min_feature_area),
            "target_segment_length": None if target_segment_length is None else float(target_segment_length),
            "resample_domain_boundary": bool(resample_domain_boundary),
            "resample_region_sources": ",".join(resample_region_sources),
            "vertices_before": int(vertices_before),
            "vertices_after": int(vertices_after),
            "polygons_before": int(polygons_before),
            "polygons_after": int(polygons_after),
        }
        return self._last_cleanup_report

    def add_line_buffer(
        self,
        line: Path,
        buffer: int = 10,
        simplify_tolerance: int = 10,
        densify_dist: int = None,
        domain: shp.Polygon | Path = None,
        max_area: int = None,
        negative_buffer_after_clipping: float | int = 0,
    ):
        """Legacy alias for :meth:`add_line_feature`.

        Parameters
        ----------
        line : shapely geometry, pathlib.Path, GeoDataFrame or GeoSeries
            The linear feature; as in :meth:`add_line_feature` (a bare sequence
            of ``(x, y)`` pairs is read as a polygon ring, not a line).
        buffer : float, optional
            Half-width of the refinement band around the line, in model length
            units. Default ``10``.
        simplify_tolerance : float, optional
            ``shapely`` simplification tolerance for the band, topology
            preserved. Default ``10``; ``None`` or ``0`` skips it.
        densify_dist : float, optional
            Resample the band's exterior to points this far apart; corners
            between the samples and holes are lost. Default ``None``, no
            resampling.
        domain : shapely.geometry.Polygon or pathlib.Path, optional
            Ignored; kept for backward compatibility. The band is clipped to the
            stored domain when one is set. Default ``None``.
        max_area : float, optional
            Maximum triangle area inside the band, in model length units squared.
            Default ``None``, which means ``1``.
        negative_buffer_after_clipping : float, optional
            A buffer of zero or less applied after clipping, pulling the band
            back from the domain boundary. Default ``0``.

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            The band as recorded.

        Raises
        ------
        ValueError
            As :meth:`add_line_feature`.
        """
        del domain
        return self.add_line_feature(
            line,
            buffer=buffer,
            max_area=max_area if max_area is not None else 1,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
            negative_buffer_after_clipping=negative_buffer_after_clipping,
        )

    @staticmethod
    def generate_dissipating_point_cloud(
        polygon: shp.Polygon | Path = None,
        buffer_dist: int | float = 1000,
        num_buffers: int = 5,
        min_area: int | float = None,
        max_area: int | float = None,
        min_spacing: int | float = 200,
        max_spacing: int | float = 1000,
        method: str = 'power',
        exponent: int | float = 2,
    ):
        """Generate outward polygons/point clouds with increasing spacing.

        Ring ``i`` (of ``num_buffers``) is a polygon drawn through points spaced
        ``spacing[i]`` apart along the outline of ``polygon`` buffered outward by
        ``min_spacing + curve[i] * (buffer_dist - min_spacing)``, where
        ``curve`` runs from 0 to 1 by ``method`` and ``spacing`` runs linearly
        from ``min_spacing`` to ``max_spacing``. So the first ring sits
        ``min_spacing`` out with the densest vertices and the last sits
        ``buffer_dist`` out with the sparsest. Lengths are in model length units.

        Parameters
        ----------
        polygon : shapely.geometry.Polygon or pathlib.Path
            The core polygon, or a ``.shp``/``.gpkg`` whose features union to a
            single ``Polygon``. Required despite the ``None`` default: anything
            that is not a ``Polygon`` (a ``MultiPolygon`` included) fails an
            ``assert``.
        buffer_dist : float, optional
            Outward distance of the outermost ring. Default ``1000``.
        num_buffers : int, optional
            Number of rings. Default ``5``.
        min_area : float, optional
            Used only when ``min_spacing`` is ``None``, which is then
            ``sqrt(min_area)``. Default ``None``.
        max_area : float, optional
            Used only when ``max_spacing`` is ``None``, which is then
            ``sqrt(max_area)``. Default ``None``.
        min_spacing : float, optional
            Vertex spacing along the innermost ring, and its outward offset.
            Default ``200``.
        max_spacing : float, optional
            Vertex spacing along the outermost ring. Default ``1000``.
        method : {'power', 'exponential'}, optional
            Shape of ``curve``: ``"power"`` is ``t ** exponent``,
            ``"exponential"`` is ``(exp(t) - 1) / (e - 1)``, ``t`` running
            linearly from 0 to 1. Default ``"power"``.
        exponent : float, optional
            Exponent of the ``"power"`` curve; ignored by ``"exponential"``.
            Default ``2``.

        Returns
        -------
        list of shapely.geometry.Polygon
            ``polygon`` itself followed by the ``num_buffers`` rings, innermost
            first.

        Raises
        ------
        ValueError
            If ``method`` is neither ``"power"`` nor ``"exponential"``.
        """
        if isinstance(polygon, Path):
            polygon = read_shp_gpkg(polygon).union_all()
        assert isinstance(polygon, shp.Polygon), 'provided argument cannot be converted to a polygon'

        min_spacing = np.sqrt(min_area) if min_spacing is None else min_spacing
        max_spacing = np.sqrt(max_area) if max_spacing is None else max_spacing

        if method == 'exponential':
            linear_space = np.linspace(0, 1, num_buffers)
            curve = np.exp(linear_space) - 1
            curve /= curve[-1]
        elif method == 'power':
            linear_space = np.linspace(0, 1, num_buffers)
            curve = linear_space ** exponent
        else:
            raise ValueError("Method must be 'exponential' or 'power'.")

        buffer_distances = min_spacing + curve * (buffer_dist - min_spacing)
        spacings = np.linspace(min_spacing, max_spacing, num_buffers)

        pols = [polygon]
        for dist, spac in zip(buffer_distances, spacings):
            pnts = []
            for interp_len in np.arange(0, polygon.buffer(dist).exterior.length, spac):
                pnts.append(polygon.buffer(dist).exterior.line_interpolate_point(interp_len))
            pols.append(shp.Polygon(pnts))

        return pols

    def add_dissipating_buffer_zones(
        self,
        max_area: int | float = 100,
        **kwargs,
    ):
        """Add several outward buffer zones with progressively lower priority.

        Builds the polygons of :meth:`generate_dissipating_point_cloud` and adds
        each through :meth:`add_region_polygon` with source
        ``"dissipating_buffer"``: the core polygon at priority ``0`` and ring
        ``i`` at priority ``-i``, so each ring claims only the band outside the
        one before it. Every zone gets the SAME ``max_area``; the grading comes
        from the rings' vertex spacing, which Triangle keeps as mesh vertices --
        unless cleanup simplifies them away first: every :meth:`build_mesh`
        preset but ``"baseline"`` simplifies at ``2.0``, which drops the evenly
        spaced vertices along a ring's straight runs (measured, a ring around a
        square went from 48 vertices to 11), so pass ``simplify_tolerance=None``
        or ``cleanup=False`` to keep them.

        Parameters
        ----------
        max_area : float, optional
            Maximum triangle area for every zone, in model length units squared.
            Default ``100``.
        **kwargs
            Forwarded to :meth:`generate_dissipating_point_cloud`: ``polygon``
            (required), ``buffer_dist``, ``num_buffers``, ``min_area``,
            ``min_spacing``, ``max_spacing``, ``method``, ``exponent``. Its own
            ``max_area`` cannot be passed (the name is taken here), so
            ``max_spacing=None`` raises ``TypeError``.
        """
        polys = self.generate_dissipating_point_cloud(**kwargs)
        for i, geom in enumerate(polys):
            if i == 0:
                self.add_region_polygon(geom, max_area=max_area, source="dissipating_buffer")
            else:
                self.add_region_polygon(geom, max_area=max_area, priority=-i, source="dissipating_buffer")

    def add_polygons_with_multiregions(
        self,
        shp_gpkg,
        buff_num=5,
        buff_sep=200,
        min_area=50,
        max_area=5_000,
    ):
        """Legacy helper that adds concentric buffered regions around polygons.

        Adds ``buff_num + 1`` regions through :meth:`add_region_polygon`, with
        source ``"multiregion"``: the unioned polygons themselves at priority
        ``0``, then copies buffered outward by ``buff_sep``, ``2 * buff_sep``,
        ... at priorities ``-1``, ``-2``, ... Their maximum triangle areas run
        linearly from ``min_area`` (innermost) to ``max_area`` (outermost). Each
        is simplified at a tolerance of 10 and resampled every 50 model length
        units (``densify_dist``, so holes are dropped), both fixed.

        Parameters
        ----------
        shp_gpkg : pathlib.Path or list of pathlib.Path
            ``.shp``/``.gpkg`` file(s) whose features are unioned into the core
            geometry. Not a ``str``.
        buff_num : int, optional
            Number of buffered copies around the core. Default ``5``.
        buff_sep : float, optional
            Distance between successive buffers, in model length units. Default
            ``200``.
        min_area : float, optional
            Maximum triangle area of the core region, in model length units
            squared. Default ``50``.
        max_area : float, optional
            Maximum triangle area of the outermost buffer. Default ``5_000``.
        """
        poly = read_shp_gpkg(shp_gpkg).union_all()
        pols = [poly] + [poly.buffer(buff_sep * mult) for mult in range(1, buff_num + 1)]
        buff_areas = np.linspace(min_area, max_area, buff_num + 1)

        assert len(buff_areas) == len(pols)
        for i, pol in enumerate(pols):
            self.add_region_polygon(
                pol,
                max_area=float(buff_areas[i]),
                simplify_tolerance=10,
                densify_dist=50,
                priority=-i,
                source="multiregion",
            )

    def validate_setup(self):
        """Validate that a buildable domain and region configuration exists."""
        if self.domain_spec is None:
            raise ValueError("A domain must be defined before building the mesh")
        for idx, region in enumerate(self.region_specs):
            if region.max_area is None:
                raise ValueError(f"Region {idx} is missing max_area")
        return True

    def _domain_interior(self) -> Polygon | MultiPolygon | None:
        """Return the domain inset by the clip tolerance (for trimming regions).

        A refinement region whose boundary coincides with the domain boundary
        makes Triangle abort with a topological inconsistency. Clipping each
        region to a slightly inset domain keeps its segments strictly interior.
        ``None`` when no domain is set; the original domain when the tolerance is
        non-positive or would empty the domain.
        """
        if self.domain_spec is None:
            return None
        geometry = self.domain_spec.geometry
        tol = self.domain_clip_tolerance
        if tol is None:
            minx, miny, maxx, maxy = geometry.bounds
            tol = ((maxx - minx) ** 2 + (maxy - miny) ** 2) ** 0.5 * 1e-4
        if tol <= 0:
            return geometry
        inset = geometry.buffer(-tol)
        return inset if not inset.is_empty else geometry

    @staticmethod
    def _polygonal_parts(geometry: Any) -> list[Polygon]:
        """Return the non-empty polygonal pieces of any geometry (ignore lines/points)."""
        if isinstance(geometry, Polygon):
            parts = [geometry]
        elif isinstance(geometry, MultiPolygon):
            parts = list(geometry.geoms)
        elif isinstance(geometry, GeometryCollection):
            parts = []
            for geom in geometry.geoms:
                if isinstance(geom, Polygon):
                    parts.append(geom)
                elif isinstance(geom, MultiPolygon):
                    parts.extend(geom.geoms)
        else:
            parts = []
        return [poly for poly in parts if not poly.is_empty and poly.area > 0]

    def prepare(self) -> list[PreparedRegion]:
        """Resolve regions into Triangle-ready polygons and interior points."""
        self.validate_setup()

        prepared: list[PreparedRegion] = []
        occupied_geometries: list[Polygon] = []
        occupied_points: list[tuple[float, float]] = []
        interior = self._domain_interior()

        sorted_specs = sorted(
            self.region_specs,
            key=lambda region: (-region.priority, region.geometry.area),
        )

        for idx, region in enumerate(sorted_specs):
            geometry = region.geometry
            if interior is not None:
                geometry = geometry.intersection(interior)
            region_polygons = self._polygonal_parts(geometry)
            if not region_polygons:
                warnings.warn(
                    f"Region '{region.label or f'region_{idx}'}' is empty after "
                    "clipping to the domain interior; skipping it.",
                    stacklevel=2,
                )
                continue
            for poly_idx, polygon in enumerate(region_polygons):
                label = region.label or f"region_{idx}"
                if len(region_polygons) > 1:
                    label = f"{label}_{poly_idx}"
                claim_geometry = self._claimed_region_geometry(
                    polygon,
                    occupied_geometries,
                    label=label,
                    priority=region.priority,
                )

                point = self._choose_region_point(
                    claim_geometry,
                    existing_points=occupied_points,
                    point_hint=region.point_hint,
                    max_area=region.max_area,
                )

                prepared_region = PreparedRegion(
                    geometry=polygon,
                    claim_geometry=claim_geometry,
                    point=point,
                    max_area=region.max_area,
                    label=label,
                    priority=region.priority,
                    attribute=region.attribute,
                    source=region.source,
                    segments=region.segments,
                )
                prepared.append(prepared_region)
                occupied_geometries.append(claim_geometry)
                occupied_points.append(point)

        self._prepared_regions = prepared
        self._sync_triangle_inputs()
        return prepared

    def _sync_triangle_inputs(self):
        """Rebuild the underlying Triangle inputs (polygons, regions, nodes) from the current specs."""

        self._polygons = []
        self._holes = []
        self._regions = []
        self._nodes = None

        if self.domain_spec is None:
            return

        # Every part's exterior, then every hole as a ring polygon PLUS a point inside it:
        # the order FloPy's tri2vor reads (domain, holes, then regions). FloPy's own
        # add_polygon hands an interior ring to add_hole, which wants a point -- the ring
        # was never a segment and the "point" was its first two vertices, so a boundary
        # with a hole did not build (ledger 317).
        parts = list(self._iter_polygons(self.domain_spec.geometry))
        for polygon in parts:
            super().add_polygon(polygon, ignore_holes=True)
        holes = []
        for polygon in parts:
            for ring in polygon.interiors:
                hole = Polygon(ring)
                super().add_polygon(hole, ignore_holes=True)
                # The seed must be in the hole and in no other part: an island in the hole
                # took the seed at its centre, and Triangle ate the island and meshed the water.
                islands = [other for other in parts if other is not polygon and other.intersects(hole)]
                seedable = hole.difference(shp.union_all(islands)) if islands else hole
                Triangle.add_hole(self, seedable.representative_point().coords[0])
                holes.append(hole)

        for region in self._prepared_regions:
            if not region.segments:
                Triangle.add_region(self, point=region.point, attribute=region.attribute,
                                    maximum_area=region.max_area)
                continue
            super().add_polygon(region.geometry, ignore_holes=True)
            # A region's own ring is a SEGMENT, never a hole: handed to add_hole it carved
            # the region's middle out of the mesh -- an annular refinement lost its whole
            # disc. A ring that IS a domain hole is a segment already, and adding it a
            # second time, from an intersection's coordinates, is refused by FloPy as
            # duplicate points (measured with domain_clip_tolerance=0).
            for ring in region.geometry.interiors:
                inner = Polygon(ring)
                if not any(inner.symmetric_difference(h).area <= 1e-9 * max(h.area, 1.0) for h in holes):
                    super().add_polygon(inner, ignore_holes=True)
            Triangle.add_region(
                self,
                point=region.point,
                attribute=region.attribute,
                maximum_area=region.max_area,
            )

        for manual_region in self._manual_regions:
            point, attribute, maximum_area = manual_region
            Triangle.add_region(self, point=point, attribute=attribute, maximum_area=maximum_area)

        excluded = self._geometry_coordinate_tuples(self.domain_spec.geometry)
        for region in self._prepared_regions:
            excluded.extend(self._geometry_coordinate_tuples(region.geometry))

        node_points = self._dedupe_node_points(
            [*self._point_constraints, *self._optimization_points],
            excluded=excluded,
        )
        node_points = self._place_point_constraints(node_points, n_parts=len(parts), n_holes=len(holes))
        if node_points:
            self._nodes = np.asarray(node_points, dtype=float)

    def _place_point_constraints(self, node_points, *, n_parts: int, n_holes: int):
        """The user's point constraints placed where Triangle can mesh them; the rest unchanged.

        A free node strictly INSIDE the domain and nearer a segment than the file writes makes
        Triangle bisect that segment toward it: below a few ulp Triangle crashed, hung or ran
        out of memory; below about one written unit the grid could not be written (ledger
        316, measured over 394 offsets). Wells snapped onto a boundary in GIS, or reprojected
        through latitude and longitude, land there. Within ``tol`` -- two written units
        (:func:`~myflopy.grid_identity.written_resolution`: 0.02 ft at state plane), floored
        at 4000 ulp -- of a segment Triangle is given (the boundary, a hole, a region's ring,
        the background sizing ring), a point is inserted INTO that ring as a vertex; within
        ``tol`` of a ring vertex, of another point, or of both segments at a corner, it is
        merged into it. A point exactly ON a segment is left to Triangle, which splits the
        segment there, unless that leaves a piece shorter than ``tol``. A point OUTSIDE the
        domain is left as it was -- Triangle meshes around it and it gets no cell, as before --
        and a warning says so. Anything farther than ``tol`` is untouched, so every grid whose
        points are clear of the edges keeps its pin.
        """

        from collections import defaultdict

        from scipy.spatial import cKDTree
        from shapely import STRtree
        from shapely.geometry import LineString

        user = {(float(x), float(y)) for x, y in self._point_constraints}
        rings = [[(float(v[0]), float(v[1])) for v in ring] for ring in self._polygons]
        if not user or not node_points or not rings:
            return node_points
        coords = np.array([v for ring in rings for v in ring], dtype=float)
        tol = max(2.0 * float(written_resolution(coords).max()), 4000.0 * float(np.spacing(np.abs(coords).max())))
        segments, owner = [], []
        for ri, ring in enumerate(rings):
            for j in range(len(ring)):
                segments.append(LineString([ring[j], ring[(j + 1) % len(ring)]]))
                owner.append((ri, j))
        tree, corners = STRtree(segments), cKDTree(coords)
        domain_rings = rings[: n_parts + n_holes]

        def collinear(k, a, b):
            c, d = rings[owner[k][0]][owner[k][1]], rings[owner[k][0]][(owner[k][1] + 1) % len(rings[owner[k][0]])]
            return all(self._segment_side([a, b], 0, q)[0] == 0 for q in (c, d))

        kept, placed = [], []  # placed: points that will be nodes INSIDE or on the domain
        inserts, moved, merged, outside = defaultdict(list), [], [], []
        on_segment = defaultdict(list)  # a segment -> the user points kept exactly on it
        for p in node_points:
            if p not in user:
                kept.append(p)
                continue
            near = [int(i) for i in tree.query(Point(p), predicate="dwithin", distance=tol)]
            side = {i: self._segment_side(rings[owner[i][0]], owner[i][1], p) for i in near}
            to_vertex, vertex = corners.query(p)
            crowd = min([float(to_vertex)] + [float(np.hypot(p[0] - q[0], p[1] - q[1])) for q in placed])
            exact = [i for i in near if side[i][0] == 0 and 0 <= side[i][1] <= 1]
            if exact:
                if 0 < crowd < tol:
                    merged.append((p, crowd))
                    continue
                kept.append(p)
                placed.append(p)
                for i in exact:
                    on_segment[owner[i]].append((side[i][1], p))
                continue
            if not self._inside_rings(p, domain_rings, exact=bool(near) or crowd < tol):
                kept.append(p)  # Triangle meshes around it, as before; never merged into
                outside.append(p)
                continue
            if 0 < crowd < tol:
                merged.append((p, crowd))
                continue
            if not near:
                kept.append(p)
                placed.append(p)
                continue
            first = min(near, key=lambda i: side[i][2])
            ri, j = owner[first]
            a, b = rings[ri][j], rings[ri][(j + 1) % len(rings[ri])]
            # Every copy of this segment -- two zones along one contact, a region on the sizing
            # ring -- takes the point, or the copy left straight is the pile all over again.
            copies = [k for k in near if collinear(k, a, b)]
            ends = [q for k in near for q in (rings[owner[k][0]][owner[k][1]],
                                              rings[owner[k][0]][(owner[k][1] + 1) % len(rings[owner[k][0]])])]
            at_an_end = any(float(side[k][1]) * side[k][3] < tol or (1 - float(side[k][1])) * side[k][3] < tol
                            for k in copies)
            if len(copies) < len(near) or at_an_end:
                # near two lines (a corner, a wedge's tip) or a segment's end: that vertex
                nearest = min(ends, key=lambda q: np.hypot(p[0] - q[0], p[1] - q[1]))
                merged.append((p, float(np.hypot(p[0] - nearest[0], p[1] - nearest[1]))))
                continue
            for k in copies:
                inserts[owner[k]].append((side[k][1], p))
            placed.append(p)
            moved.append((p, side[first][2]))
        # A point kept exactly ON a segment that now bends to take another would sit off the
        # new chord: it goes into the ring too, in order.
        for key in list(inserts):
            for t, q in on_segment.get(key, []):
                inserts[key].append((t, q))
                if q in kept:
                    kept.remove(q)
        for ri, j in sorted(inserts, key=lambda key: (key[0], -key[1])):
            ring = list(self._polygons[ri])
            self._polygons[ri] = ring[: j + 1] + [q for _, q in sorted(inserts[(ri, j)])] + ring[j + 1:]

        messages = []
        if moved or merged:
            worst = max(d for _, d in moved + merged)
            messages.append(
                f"{len(moved)} point constraint(s) within {tol:.3g} of an edge were placed ON it, and "
                f"{len(merged)} within {tol:.3g} of a ring vertex, a corner or another point were merged "
                f"into it (moved at most {worst:.3g}): nearer than the DISV file writes, Triangle "
                "crashes, hangs or builds a grid that cannot be written."
            )
        if outside:
            messages.append(
                f"{len(outside)} point constraint(s) lie outside the meshed domain and get no cell, "
                f"e.g. ({outside[0][0]:.10g}, {outside[0][1]:.10g}). A point drawn on the boundary can "
                "round outside it, or sit outside after the profile simplifies the boundary."
            )
        issued = self.__dict__.setdefault("_point_messages_issued", set())
        for message in messages:
            if message not in issued:
                issued.add(message)
                logger.warning(message)
                warnings.warn(message, UserWarning, stacklevel=4)
        return kept

    @staticmethod
    def _segment_side(ring, j, p):
        """``(cross, t, distance, length)`` of ``p`` against segment ``j`` of ``ring``: the
        cross product and the foot's parameter EXACTLY (``Fraction``), since a float distance
        reads 0 for a point 1e-15 inside; distance and length as floats."""

        from fractions import Fraction

        a, b = ring[j], ring[(j + 1) % len(ring)]
        ax, ay, bx, by, px, py = (Fraction(v) for v in (*a, *b, *p))
        dx, dy = bx - ax, by - ay
        cross = dx * (py - ay) - dy * (px - ax)
        squared = dx * dx + dy * dy
        t = (dx * (px - ax) + dy * (py - ay)) / squared if squared else Fraction(0)
        length = float(squared) ** 0.5
        return cross, t, abs(float(cross)) / length if length else float(np.hypot(*(np.subtract(p, a)))), length

    @staticmethod
    def _inside_rings(p, rings, *, exact):
        """Whether ``p`` is inside the domain bounded by ``rings`` (every part's exterior and
        every hole): an odd number of them around it, so an island in another part's hole is
        inside. Exactly, by crossing number on ``Fraction``s, when ``exact``."""

        if not exact:
            point = Point(p)
            return sum(Polygon(ring).contains(point) for ring in rings) % 2 == 1

        from fractions import Fraction

        px, py = Fraction(p[0]), Fraction(p[1])

        def crossing(ring):
            inside = False
            for k in range(len(ring)):
                (ax, ay), (bx, by) = ((Fraction(c) for c in ring[k]), (Fraction(c) for c in ring[(k + 1) % len(ring)]))
                if (ay > py) != (by > py):
                    left = (bx - ax) * (py - ay) - (px - ax) * (by - ay)
                    if (left > 0) == (by > ay):
                        inside = not inside
            return inside

        return sum(crossing(r) for r in rings) % 2 == 1

    @staticmethod
    def _geometry_coordinate_tuples(geometry: Polygon | MultiPolygon | GeometryCollection) -> list[tuple[float, float]]:
        """Every exterior and interior-ring vertex of a geometry as ``(x, y)`` tuples."""

        coords: list[tuple[float, float]] = []
        for polygon in TriangleGrid._iter_polygons(geometry):
            coords.extend((float(x), float(y)) for x, y in polygon.exterior.coords)
            for interior in polygon.interiors:
                coords.extend((float(x), float(y)) for x, y in interior.coords)
        return coords

    @staticmethod
    def _dedupe_node_points(
        points: list[tuple[float, float]],
        *,
        excluded: list[tuple[float, float]] | None = None,
        digits: int = 8,
    ) -> list[tuple[float, float]]:
        """Drop duplicate node points (and any coinciding with ``excluded``), rounded to ``digits``."""

        seen = set()
        if excluded:
            seen.update((round(float(x), digits), round(float(y), digits)) for x, y in excluded)

        deduped: list[tuple[float, float]] = []
        for x, y in points:
            key = (round(float(x), digits), round(float(y), digits))
            if key in seen:
                continue
            seen.add(key)
            deduped.append((float(x), float(y)))
        return deduped

    @staticmethod
    def _thin_points(
        points: list[tuple[float, float]],
        *,
        min_spacing: float,
        excluded: list[tuple[float, float]] | None = None,
    ) -> list[tuple[float, float]]:
        """Greedily thin points so no two accepted points are closer than ``min_spacing``."""

        if min_spacing <= 0:
            return points

        accepted: list[tuple[float, float]] = []
        all_points = [] if excluded is None else list(excluded)
        for point in points:
            candidate = Point(point)
            if all(candidate.distance(Point(existing)) > min_spacing for existing in all_points):
                accepted.append(point)
                all_points.append(point)
        return accepted

    @staticmethod
    def _limit_points(
        points: list[tuple[float, float]],
        *,
        max_points: int,
    ) -> list[tuple[float, float]]:
        """Downsample to at most ``max_points`` points by even index selection (no-op if under)."""

        if max_points <= 0 or len(points) <= max_points:
            return points

        indices = np.linspace(0, len(points) - 1, max_points, dtype=int)
        return [points[int(index)] for index in indices]

    def preview_regions(self) -> gpd.GeoDataFrame:
        """Return a GeoDataFrame preview of prepared region polygons."""
        prepared = self.prepare()
        return gpd.GeoDataFrame(
            {
                "label": [region.label for region in prepared],
                "max_area": [region.max_area for region in prepared],
                "claim_area": [region.claim_geometry.area for region in prepared],
                "priority": [region.priority for region in prepared],
                "source": [region.source for region in prepared],
                "point_x": [region.point[0] for region in prepared],
                "point_y": [region.point[1] for region in prepared],
            },
            geometry=[region.geometry for region in prepared],
        )

    def get_region_points(self) -> gpd.GeoDataFrame:
        """Return the prepared region points as a GeoDataFrame."""
        prepared = self.prepare()
        return gpd.GeoDataFrame(
            {
                "label": [region.label for region in prepared],
                "max_area": [region.max_area for region in prepared],
                "priority": [region.priority for region in prepared],
                "source": [region.source for region in prepared],
            },
            geometry=[Point(region.point) for region in prepared],
        )

    def _build_triangle_mesh(self, verbose=False):
        """Prepare inputs, ensure the workspace exists, and run the underlying Triangle build."""

        self.prepare()
        Path(self.model_ws).mkdir(parents=True, exist_ok=True)
        return self._bounded_build(verbose=verbose)

    #: Bounds on one run of the Triangle program: its address space, its wall time and how
    #: much output is kept. A mesh of ~1e6 cells takes seconds and ~0.3 GB.
    TRIANGLE_MEMORY_LIMIT = 4 * 2**30
    TRIANGLE_TIMEOUT = 300.0
    TRIANGLE_OUTPUT_LIMIT = 64 * 2**20

    def _bounded_build(self, verbose=False):
        """FloPy's ``Triangle.build`` -- the same files, switches and results -- under bounds.

        FloPy runs Triangle through ``subprocess.check_output`` with no timeout, holding all
        of its output in memory; a runaway run grew by 0.3-0.5 GB a second there. This caps
        Triangle's address space (POSIX), kills it after :attr:`TRIANGLE_TIMEOUT` seconds or
        :attr:`TRIANGLE_OUTPUT_LIMIT` bytes of output, and raises :class:`TriangleFailed`.
        """

        import os
        import threading

        self.clean()
        self._write_nodefile(os.path.join(self.model_ws, f"{self.file_prefix}.0.node"))
        self._write_polyfile(os.path.join(self.model_ws, f"{self.file_prefix}.0.poly"))
        cmds = [self.exe_name, f"-a{self.maximum_area}" if self.maximum_area is not None else "-a"]
        if self.angle is not None:
            cmds.append(f"-q{self.angle}")
        if self.additional_args is not None:
            cmds += list(self.additional_args)
        cmds += ["-A", "-p", "-V", "-D", "-e", "-n", f"{self.file_prefix}.0"]

        limit = self.TRIANGLE_MEMORY_LIMIT

        def cap_memory():  # runs in the child, before exec: lower the SOFT limit only
            import resource

            try:
                soft, hard = resource.getrlimit(resource.RLIMIT_AS)
                cap = limit if hard == resource.RLIM_INFINITY else min(limit, hard)
                if soft == resource.RLIM_INFINITY or soft > cap:
                    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
            except (ValueError, OSError):
                pass  # a cap that cannot be set never blocks the run

        try:
            process = subprocess.Popen(cmds, cwd=self.model_ws, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                       preexec_fn=cap_memory if os.name == "posix" else None)
        except (OSError, subprocess.SubprocessError) as exc:
            raise TriangleFailed(-1, cmds, output=f"Triangle could not be started in {self.model_ws}: {exc}") from exc
        timed_out = threading.Event()

        def kill_on_timeout():
            timed_out.set()
            process.kill()

        timer = threading.Timer(self.TRIANGLE_TIMEOUT, kill_on_timeout)
        timer.start()
        output, runaway = bytearray(), False
        try:
            for chunk in iter(lambda: process.stdout.read(1 << 16), b""):
                if len(output) + len(chunk) > self.TRIANGLE_OUTPUT_LIMIT:
                    runaway = True
                    process.kill()
                    break
                output += chunk
            returncode = process.wait()
        finally:
            timer.cancel()
            if process.poll() is None:  # an interrupt or error while reading: never leave it running
                process.kill()
                process.wait()
            process.stdout.close()
        text = output.decode(errors="replace")
        if returncode != 0 or runaway:
            if timed_out.is_set():
                why = f"was still running after {self.TRIANGLE_TIMEOUT:g} s"
            elif runaway:
                why = f"wrote more than {self.TRIANGLE_OUTPUT_LIMIT // 2**20} MB of output"
            elif "Out of memory" in text or returncode == -9:
                why = f"ran out of memory (capped at {limit / 2**30:g} GB)"
            else:
                why = f"exited with code {returncode}"
            last = next((line for line in reversed(text.splitlines()) if line.strip()), "")
            raise TriangleFailed(
                returncode or -9, cmds,
                output=(f"Triangle {why} meshing {self.model_ws}; its last line: {last[:160]!r}. A point "
                        "constraint within rounding of an edge, or two edges drawn a hair apart, "
                        "can do this (ledger 316)."),
            )
        if verbose:
            print(text)
        self._load_results()
        self.ncpl = self.ele.shape[0]
        self.nvert = self.node.shape[0]
        self.verts = np.array(self.node[["x", "y"]].tolist(), float)
        self.iverts = [[row[1], row[2], row[3]] for row in self.ele]
        return None

    def _shared_nodes(self):
        """``(nodes, rings)`` with ONE Triangle node per distinct ``(x, y)``, or ``None``.

        FloPy writes every ring's vertices as nodes of their own, so a point two rings share
        -- two refinement zones along one contact, a region corner on the background sizing
        ring's corner -- became two input nodes. Triangle marks one unused and still copies
        both to its output, and FloPy's ``tri2vor`` then refused the mesh with a bare
        "duplicate points" ``Exception`` (ledger 316; Triangle's triangles were measured the
        same with and without the duplicate). ``None`` when nothing coincides, and under
        Triangle's ``-j`` (which already jettisons it and numbers by its own sort): FloPy's
        writers then run unchanged, so every grid that builds today gets the same bytes.
        """

        index, nodes, rings, shared = {}, [], [], False
        for polygon in self._polygons:
            ids = []
            for vertex in polygon:
                key = (float(vertex[0]), float(vertex[1]))
                if key in index:
                    shared = True
                else:
                    index[key] = len(nodes)
                    nodes.append((vertex[0], vertex[1]))
                ids.append(index[key])
            rings.append(ids)
        extra = [] if self._nodes is None else [(self._nodes[i, 0], self._nodes[i, 1]) for i in range(self._nodes.shape[0])]
        for x, y in extra:
            key = (float(x), float(y))
            if key in index:
                shared = True
            else:
                index[key] = len(nodes)
                nodes.append((x, y))
        jettison = any(isinstance(a, str) and a.startswith("-") and "j" in a[1:] for a in (self.additional_args or ()))
        if not shared or jettison:
            return None
        return nodes, rings

    def _write_nodefile(self, fname):
        """FloPy's ``.node`` writer, with a point several rings share written once."""

        shared = self._shared_nodes()
        if shared is None:
            return super()._write_nodefile(fname)
        nodes, _ = shared
        with open(fname, "w") as f:
            f.write(f"{len(nodes)} 2 0 0\n")
            for ip, (x, y) in enumerate(nodes):
                f.write(f"{ip} {x} {y}\n")
        return None

    def _write_polyfile(self, fname):
        """FloPy's ``.poly`` writer, each ring's segments on the shared node ids.

        Segment count, order and markers (``iseg + 1``) are FloPy's, and ``_polygons``,
        ``_holes`` and ``_regions`` are untouched, so ``tri2vor``'s boundary-marker reading
        is unchanged.
        """

        shared = self._shared_nodes()
        if shared is None:
            return super()._write_polyfile(fname)
        _, rings = shared
        with open(fname, "w") as f:
            f.write("0 0 0 0\n")
            f.write(f"{sum(len(ring) for ring in rings)} 1\n")
            iseg = 0
            for ids in rings:
                for i in range(len(ids)):
                    f.write(f"{iseg} {ids[i]} {ids[(i + 1) % len(ids)]} {iseg + 1}\n")
                    iseg += 1
            f.write(f"{len(self._holes)}\n")
            for i, hole in enumerate(self._holes):
                f.write(f"{i} {hole[0]} {hole[1]}\n")
            f.write(f"{len(self._regions)}\n")
            for i, (pt, attribute, maxarea) in enumerate((r[0], r[1], r[2]) for r in self._regions):
                f.write(f"{i} {pt[0]} {pt[1]} {attribute} {-1.0 if maxarea is None else maxarea}\n")
        return None

    @staticmethod
    def _build_voronoi_for_mesh(mesh):
        """Wrap a triangulation in a :class:`VoronoiGridPlus`, silencing benign numpy warnings."""

        from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Mean of empty slice.",
                category=RuntimeWarning,
            )
            warnings.filterwarnings(
                "ignore",
                message="invalid value encountered in divide",
                category=RuntimeWarning,
            )
            return VoronoiGridPlus(mesh)

    def build(self, verbose=False):
        """Prepare Triangle inputs and run a single direct Triangle build.

        Runs :meth:`prepare`, creates ``model_ws`` if needed, and calls FloPy's
        ``Triangle.build``, which writes the Triangle input files there, runs the
        executable and loads the mesh (``node``, ``ele``, ``ncpl``, ``verts``,
        ...). No cleanup and no optimization: :meth:`build_mesh` is the
        recommended workflow.

        Parameters
        ----------
        verbose : bool, optional
            Print the Triangle program's own output. Default ``False``.

        Returns
        -------
        None
            The mesh is stored on the instance.

        Raises
        ------
        ValueError
            If no domain is set, a region has no ``max_area``, or overlap
            resolution leaves a region no area of its own (see :meth:`prepare`).
        """
        return self._build_triangle_mesh(verbose=verbose)

    def quality_report(self, *, include_voronoi: bool = False):
        """Return triangle quality metrics and optional Voronoi build checks.

        The mesh must already be built (:meth:`build` or :meth:`build_mesh`).
        The report holds triangle counts, degeneracies (duplicate vertices,
        zero-area, tiny and sliver triangles) and summary statistics of triangle
        areas, edge lengths, angles, edge ratios and neighbour-area ratios, plus
        ``cleanup_vertices_before``/``_after`` once :meth:`clean_geometry` has
        run. It is also kept on the instance for later reads.

        Parameters
        ----------
        include_voronoi : bool, optional
            Also try to build the Voronoi grid from this mesh
            (:class:`~myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus`) and add
            ``voronoi_status`` (``"built"`` or ``"error"``, with
            ``voronoi_error``), ``voronoi_cell_count``,
            ``voronoi_invalid_cell_count`` (cells of area <= 0) and Voronoi
            cell-area statistics. Default ``False``, triangle metrics only.

        Returns
        -------
        dict
            Metric name to value.

        Raises
        ------
        ValueError
            If the mesh has not been built.
        """
        report = triangle_quality_report(
            self,
            include_voronoi=include_voronoi,
            validate_voronoi=include_voronoi,
        )
        if self._last_cleanup_report:
            report["cleanup_vertices_before"] = self._last_cleanup_report["vertices_before"]
            report["cleanup_vertices_after"] = self._last_cleanup_report["vertices_after"]
        self._last_quality_report = report
        return report

    @staticmethod
    def _mesh_report_is_safe_for_voronoi(report: dict[str, float | int | str]) -> bool:
        """True if a mesh report has no duplicate/zero-area triangles and the Voronoi built."""

        return (
            int(report.get("duplicate_vertex_count", 0)) == 0
            and int(report.get("zero_area_triangle_count", 0)) == 0
            and report.get("voronoi_status", "built") == "built"
        )

    @staticmethod
    def _quality_regression_reasons(
        baseline: dict[str, float | int | str],
        candidate: dict[str, float | int | str],
    ) -> list[str]:
        """Reasons a ``candidate`` mesh is worse than ``baseline`` (empty if it is acceptable).

        Flags Voronoi-unsafe meshes and regressions in tiny/sliver counts, minimum
        angle, and edge/neighbor-area ratios beyond tolerance.
        """

        reasons: list[str] = []
        if not TriangleGrid._mesh_report_is_safe_for_voronoi(candidate):
            reasons.append("voronoi_unsafe")
            return reasons

        tiny_before = int(baseline.get("tiny_triangle_count", 0))
        tiny_after = int(candidate.get("tiny_triangle_count", 0))
        if tiny_after > tiny_before:
            reasons.append("tiny_triangle_increase")

        sliver_before = int(baseline.get("sliver_triangle_count", 0))
        sliver_after = int(candidate.get("sliver_triangle_count", 0))
        sliver_limit = sliver_before + max(10, int(sliver_before * 0.5))
        if sliver_after > sliver_limit:
            reasons.append("sliver_triangle_increase")

        min_angle_before = float(baseline.get("triangle_angle_min_overall", 0.0))
        min_angle_after = float(candidate.get("triangle_angle_min_overall", 0.0))
        min_angle_limit = max(1.0, min_angle_before * 0.85)
        if min_angle_after < min_angle_limit:
            reasons.append("min_angle_regression")

        edge_ratio_before = float(baseline.get("triangle_edge_ratio_mean", 0.0))
        edge_ratio_after = float(candidate.get("triangle_edge_ratio_mean", 0.0))
        if edge_ratio_before > 0 and edge_ratio_after > (edge_ratio_before * 1.1):
            reasons.append("edge_ratio_regression")

        neighbor_ratio_before = float(baseline.get("neighbor_area_ratio_mean", 0.0))
        neighbor_ratio_after = float(candidate.get("neighbor_area_ratio_mean", 0.0))
        if neighbor_ratio_before > 0 and neighbor_ratio_after > (neighbor_ratio_before * 1.1):
            reasons.append("neighbor_area_ratio_regression")

        return reasons

    @staticmethod
    def _free_vertex_points(
        *,
        node,
        protected_geometries: list[Polygon] | None = None,
        tolerance: float = 1e-6,
    ) -> list[tuple[float, float]]:
        """The ``(x, y)`` coordinates of mesh vertices free to move (not on a protected geometry)."""

        _, free_ids = get_fixed_and_free_vertex_ids(
            node=node,
            protected_geometries=protected_geometries,
            tolerance=tolerance,
        )
        row_lookup = {int(row["ivert"]): row for row in node}
        return [
            (float(row_lookup[vertex]["x"]), float(row_lookup[vertex]["y"]))
            for vertex in free_ids
        ]

    def _resolve_protected_labels(
        self,
        protected_labels: list[str] | None = None,
        protect_sources: tuple[str, ...] | None = None,
    ) -> list[str]:
        """The set of region labels to protect from optimization: explicit labels plus any from ``protect_sources``."""

        labels = [] if protected_labels is None else list(dict.fromkeys(protected_labels))
        if protect_sources:
            self.prepare()
            for region in self._prepared_regions:
                if region.source in protect_sources and region.label not in labels:
                    labels.append(region.label)
        return labels

    def optimize_seeds(
        self,
        *,
        iterations: int = 3,
        damping: float = 0.35,
        min_move: float = 1e-3,
        protected_labels: list[str] | None = None,
        max_optimization_points: int = 0,
        verbose: bool = False,
    ):
        """Run the experimental constrained CVT/Lloyd optimization loop.

        The optimizer keeps protected feature vertices fixed, rebuilds the mesh
        after each move, and rejects candidate meshes that would materially
        degrade quality or fail Voronoi conversion.

        Builds the mesh first if it has not been built. Each round moves every
        free vertex part of the way toward the centroid of its Voronoi cell and
        re-runs Triangle with the free vertices as fixed input nodes. Every vertex
        on an input segment is always fixed -- the domain boundary, its holes and
        every refinement region's ring, the background ``max_area`` region's
        included -- because FloPy gives each segment a nonzero boundary marker,
        which Triangle passes on to the vertices along it. A move that would
        leave the domain's largest part is not made. A rejected round restores
        the last accepted mesh and ends the loop. Distances are in model length
        units.

        Parameters
        ----------
        iterations : int, optional
            Maximum number of relax-and-rebuild rounds. The loop also stops early
            when no vertex would move by at least ``min_move``, or when a round
            is rejected. Default ``3``; ``0`` builds the Voronoi grid once and
            moves nothing (status ``"not_attempted"``).
        damping : float, optional
            Fraction of the way from a vertex to its cell centroid that it moves
            per round. Values outside ``[0, 1]`` are clipped when applied (the
            report records the value passed). Default ``0.35``.
        min_move : float, optional
            A vertex whose proposed move is shorter than this stays put. Default
            ``1e-3``.
        protected_labels : list of str, optional
            Labels of prepared regions (as :meth:`preview_regions` shows them,
            ``"_<j>"`` suffixes included) to protect: every mesh vertex within
            ``region_point_tolerance`` (or its automatic value) of such a region's
            boundary is held fixed. The vertices ON every region's ring are fixed
            whether or not it is protected (see above), so this adds only those
            near a ring but not on it. A label matching no region protects
            nothing, silently. Default ``None``, no region protected.
        max_optimization_points : int, optional
            Recorded in the returned report only; it does not limit the number of
            relaxed points. Default ``0``.
        verbose : bool, optional
            Print progress for each round, and the Triangle program's output for
            each rebuild. Default ``False``.

        Returns
        -------
        dict
            The optimization report: ``status`` (for example ``"applied"``,
            ``"converged_no_move"``, ``"fallback_original_mesh"`` or
            ``"not_attempted"``), ``iterations_run``, ``max_move``,
            ``moved_vertex_count``, ``rejection_reasons`` and the settings used.
        """
        if not hasattr(self, "node") or not hasattr(self, "ele"):
            self._build_triangle_mesh(verbose=verbose)

        if verbose:
            print("[TriangleGrid] Starting constrained CVT/Lloyd optimization.")

        protected_labels = protected_labels or []
        protected_geometries = [
            region.geometry
            for region in self._prepared_regions
            if region.label in protected_labels
        ]

        tolerance = self._region_point_spacing(geometry=self.domain_geometry)
        baseline_free_point_count = len(
            self._free_vertex_points(
                node=self.node,
                protected_geometries=protected_geometries,
                tolerance=tolerance,
            )
        )
        accepted_points: list[tuple[float, float]] = []
        accepted_report = {
            "status": "not_attempted",
            "iterations_run": 0,
            "max_move": 0.0,
            "moved_vertex_count": 0,
            "free_vertex_count": baseline_free_point_count,
        }
        regression_reasons: list[str] = []
        attempts = 0
        domain_polygon = max(self._iter_polygons(self.domain_geometry), key=lambda polygon: polygon.area)
        try:
            current_vor = self._build_voronoi_for_mesh(self)
        except Exception as exc:  # noqa: BLE001 - flopy raises a BARE Exception
            # Not narrowable. `flopy.utils.voronoi` contains a literal
            # `raise Exception(...)` for a degenerate triangulation (verified
            # against flopy 3.10), so the only clause that catches the case this
            # guard exists for is `Exception` itself. The failure is REPORTED,
            # not swallowed: it becomes a `voronoi_build_failed` reason the
            # caller acts on.
            regression_reasons = [f"voronoi_build_failed:{exc}"]
            self._last_optimization_report = {
                "method": "constrained_cvt_lloyd",
                "iterations_requested": int(iterations),
                "iterations_run": 0,
                "damping": float(damping),
                "min_move": float(min_move),
                "free_vertex_count": int(accepted_report["free_vertex_count"]),
                "moved_vertex_count": 0,
                "max_move": 0.0,
                "candidate_point_count": 0,
                "optimization_point_count": 0,
                "max_optimization_points": int(max_optimization_points),
                "attempts": 0,
                "status": "fallback_original_mesh",
                "baseline_voronoi_status": "error",
                "rejection_reasons": ",".join(regression_reasons),
            }
            if verbose:
                print(
                    "[TriangleGrid] "
                    f"Optimization finished with status={self._last_optimization_report['status']}."
                )
            return self._last_optimization_report

        quality_before = triangle_quality_report(self, validate_voronoi=False)
        quality_before["voronoi_status"] = "built"
        quality_before["voronoi_cell_count"] = int(len(current_vor.vor_list))
        accepted_quality = quality_before

        for iteration in range(iterations):
            attempts += 1
            if verbose:
                print(f"[TriangleGrid] Lloyd iteration {iteration + 1}/{iterations}: building Voronoi centroids.")

            result = cvt_relaxed_full_seed_points(
                node=self.node,
                vor_points=current_vor.points,
                vor_polygons=current_vor.vor_list,
                domain=domain_polygon,
                damping=damping,
                min_move=min_move,
                protected_geometries=protected_geometries,
                tolerance=tolerance,
            )

            if verbose:
                print(
                    "[TriangleGrid] "
                    f"free seeds={result.free_vertex_count}, moved={result.moved_vertex_count}, "
                    f"max move={result.max_move:.3f}"
                )

            if not result.points or result.moved_vertex_count == 0:
                accepted_report["status"] = "converged_no_move"
                break

            self._optimization_points = result.points
            try:
                self._build_triangle_mesh(verbose=verbose)
            except subprocess.CalledProcessError:
                regression_reasons = ["triangle_build_failed"]
                self._optimization_points = accepted_points
                self._build_triangle_mesh(verbose=verbose)
                accepted_report["status"] = "fallback_original_mesh"
                break

            try:
                candidate_vor = self._build_voronoi_for_mesh(self)
            except Exception as exc:  # noqa: BLE001 - flopy raises a BARE Exception
                # Not narrowable. `flopy.utils.voronoi` contains a literal
                # `raise Exception(...)` for a degenerate triangulation (verified
                # against flopy 3.10), so the only clause that catches the case this
                # guard exists for is `Exception` itself. The failure is REPORTED,
                # not swallowed: it becomes a `voronoi_build_failed` reason the
                # caller acts on.
                regression_reasons = [f"voronoi_build_failed:{exc}"]
                self._optimization_points = accepted_points
                self._build_triangle_mesh(verbose=verbose)
                accepted_report["status"] = "fallback_original_mesh"
                break

            quality_after = triangle_quality_report(self, validate_voronoi=False)
            quality_after["voronoi_status"] = "built"
            quality_after["voronoi_cell_count"] = int(len(candidate_vor.vor_list))
            regression_reasons = self._quality_regression_reasons(accepted_quality, quality_after)
            if regression_reasons:
                if verbose:
                    print(
                        "[TriangleGrid] "
                        f"Rejected iteration {iteration + 1}: {', '.join(regression_reasons)}"
                    )
                self._optimization_points = accepted_points
                self._build_triangle_mesh(verbose=verbose)
                accepted_report["status"] = "fallback_original_mesh"
                break

            accepted_points = self._free_vertex_points(
                node=self.node,
                protected_geometries=protected_geometries,
                tolerance=tolerance,
            )
            accepted_quality = quality_after
            current_vor = candidate_vor
            accepted_report = {
                "status": "applied",
                "iterations_run": iteration + 1,
                "max_move": float(result.max_move),
                "moved_vertex_count": int(result.moved_vertex_count),
                "free_vertex_count": int(result.free_vertex_count),
            }

            if verbose:
                print(
                    "[TriangleGrid] "
                    f"Accepted iteration {iteration + 1}: "
                    f"min angle={quality_after['triangle_angle_min_overall']:.2f}, "
                    f"slivers={quality_after['sliver_triangle_count']}, "
                    f"tiny={quality_after['tiny_triangle_count']}"
                )

            if result.max_move < min_move:
                accepted_report["status"] = "converged_small_move"
                break

        self._last_optimization_report = {
            "method": "constrained_cvt_lloyd",
            "iterations_requested": int(iterations),
            "iterations_run": int(accepted_report["iterations_run"]),
            "damping": float(damping),
            "min_move": float(min_move),
            "free_vertex_count": int(accepted_report["free_vertex_count"]),
            "moved_vertex_count": int(accepted_report["moved_vertex_count"]),
            "max_move": float(accepted_report["max_move"]),
            "candidate_point_count": int(0 if self._nodes is None else self._nodes.shape[0]),
            "optimization_point_count": int(len(self._optimization_points)),
            "max_optimization_points": int(max_optimization_points),
            "attempts": int(attempts),
            "status": str(accepted_report["status"]),
            "baseline_voronoi_status": str(quality_before.get("voronoi_status", "unknown")),
            "protected_label_count": int(len(protected_labels)),
            "protected_labels": ",".join(protected_labels),
            "rejection_reasons": ",".join(regression_reasons),
        }
        if verbose:
            print(
                "[TriangleGrid] "
                f"Optimization finished with status={self._last_optimization_report['status']}."
            )
        return self._last_optimization_report

    def build_optimized(
        self,
        *,
        cleanup: bool = True,
        snap_tolerance: float = 0.0,
        simplify_tolerance: float | None = None,
        min_feature_area: float = 0.0,
        target_segment_length: float | None = None,
        resample_domain_boundary: bool = False,
        resample_region_sources: tuple[str, ...] = ("line",),
        optimization_iterations: int = 3,
        damping: float = 0.35,
        min_move: float = 1e-3,
        protected_labels: list[str] | None = None,
        max_optimization_points: int = 250,
        verbose: bool = False,
    ):
        """Build the mesh and then run the experimental optimization pass.

        Runs :meth:`clean_geometry` (when ``cleanup``), a fresh Triangle build,
        :meth:`optimize_seeds`, and a final :meth:`quality_report`-style check
        that includes the Voronoi build. :meth:`build_mesh` calls this with a
        profile's settings; the defaults here differ from the profiles' (no
        simplification, no resampling and ``damping=0.35`` here, against
        ``simplify_tolerance=2.0``, ``target_segment_length=120.0`` and
        ``damping=0.25`` in ``"balanced"``).

        Parameters
        ----------
        cleanup : bool, optional
            Run :meth:`clean_geometry` first. Default ``True``. ``False`` makes
            the six cleanup parameters below inert.
        snap_tolerance : float, optional
            Cleanup snapping distance; see :meth:`clean_geometry`. Default
            ``0.0``, no snapping.
        simplify_tolerance : float, optional
            Cleanup simplification tolerance; see :meth:`clean_geometry`.
            Default ``None``, no simplification.
        min_feature_area : float, optional
            Cleanup: drop polygon parts smaller than this area; see
            :meth:`clean_geometry`. Default ``0.0``.
        target_segment_length : float, optional
            Cleanup resampling spacing; see :meth:`clean_geometry`. Default
            ``None``, no resampling.
        resample_domain_boundary : bool, optional
            Cleanup: resample the domain boundary too; see
            :meth:`clean_geometry`. Default ``False``.
        resample_region_sources : tuple of str, optional
            Cleanup: region source tags to resample; see :meth:`clean_geometry`.
            Default ``("line",)``.
        optimization_iterations : int, optional
            Maximum relaxation rounds; :meth:`optimize_seeds`' ``iterations``.
            Default ``3``.
        damping : float, optional
            Fraction of the move toward each cell centroid per round; see
            :meth:`optimize_seeds`. Default ``0.35``.
        min_move : float, optional
            Shortest vertex move made, in model length units; see
            :meth:`optimize_seeds`. Default ``1e-3``.
        protected_labels : list of str, optional
            Prepared-region labels to protect during optimization; see
            :meth:`optimize_seeds`. Default ``None``.
        max_optimization_points : int, optional
            Recorded in the optimization report only; it limits nothing (see
            :meth:`optimize_seeds`). Default ``250``.
        verbose : bool, optional
            Print optimization progress and the Triangle program's output.
            Default ``False``.

        Returns
        -------
        dict
            ``{"cleanup": ..., "optimization": ..., "quality": ...}``: the
            :meth:`clean_geometry` report (the previous one, or ``None``, when
            ``cleanup=False``), the :meth:`optimize_seeds` report and the
            quality report with Voronoi checks.
        """
        if cleanup:
            self.clean_geometry(
                snap_tolerance=snap_tolerance,
                simplify_tolerance=simplify_tolerance,
                min_feature_area=min_feature_area,
                target_segment_length=target_segment_length,
                resample_domain_boundary=resample_domain_boundary,
                resample_region_sources=resample_region_sources,
            )

        self._optimization_points = []
        self._build_triangle_mesh(verbose=verbose)
        optimization_report = self.optimize_seeds(
            iterations=optimization_iterations,
            damping=damping,
            min_move=min_move,
            protected_labels=protected_labels,
            max_optimization_points=max_optimization_points,
            verbose=verbose,
        )
        quality_report = triangle_quality_report(self, include_voronoi=True, validate_voronoi=True)
        self._last_quality_report = quality_report
        return {
            "cleanup": self._last_cleanup_report,
            "optimization": optimization_report,
            "quality": quality_report,
        }

    def build_mesh(
        self,
        *,
        profile: str | MeshBuildProfile = "balanced",
        optimize: bool | None = None,
        protected_labels: list[str] | None = None,
        protect_sources: tuple[str, ...] | None = ("line",),
        verbose: bool = False,
        **overrides,
    ):
        """Run the recommended high-level mesh workflow.

        Three modes, picked from the profile and ``optimize`` and returned as
        ``"mode"``: ``"optimized"`` (:meth:`build_optimized`, which cleans first
        when ``cleanup``), ``"cleanup_only"`` (:meth:`clean_geometry`, then one
        Triangle build) and ``"baseline"`` (no cleanup: :meth:`prepare`, then one
        Triangle build). Every mode ends with a quality report that includes the
        Voronoi build.

        Parameters
        ----------
        profile : str or MeshBuildProfile, optional
            Either a preset profile name or a :class:`MeshBuildProfile`
            instance. A name goes through :meth:`MeshBuildProfile.from_name`,
            which raises ``ValueError`` for an unknown one. Default
            ``"balanced"``.
        optimize : bool, optional
            ``True`` always runs the optimization pass (with ``0`` iterations it
            only reports ``"not_attempted"``); ``False`` never does. Default
            ``None`` optimizes when the profile's (or the override's)
            ``optimization_iterations`` is above zero -- so not for ``"clean"``
            or ``"baseline"``.
        protected_labels : list of str, optional
            Prepared-region labels (see :meth:`preview_regions`) to protect
            during optimization, added to those ``protect_sources`` selects; see
            :meth:`optimize_seeds` for what protection fixes. Without
            optimization they are only echoed in the report. Default ``None``.
        protect_sources : tuple of str, optional
            Region ``source`` tags (see :meth:`add_region_polygon`) whose
            regions are protected during optimization, as ``protected_labels``.
            Default ``("line",)``, the bands :meth:`add_line_feature` adds.
            ``None`` protects none. Ignored without optimization.
        verbose : bool, optional
            When ``True``, print the Triangle program's output for each build
            and progress for each optimization round. Default ``False``.
        **overrides
            Replace individual settings of the profile for this call:
            ``cleanup``, ``optimization_iterations``, the cleanup fields
            (``snap_tolerance``, ``simplify_tolerance``, ``min_feature_area``,
            ``target_segment_length``, ``resample_domain_boundary``,
            ``resample_region_sources``), ``damping``, ``min_move``, and
            :meth:`build_optimized`'s ``max_optimization_points``. In the
            optimized mode they are forwarded to :meth:`build_optimized`, which
            raises ``TypeError`` on any other name. Without optimization only the
            cleanup fields and ``cleanup`` act, and every other name -- a typo
            included -- is ignored without a warning.

        Returns
        -------
        dict
            Combined report containing the selected profile plus any cleanup,
            optimization, and mesh quality diagnostics collected along the way.
        """
        profile_obj = profile if isinstance(profile, MeshBuildProfile) else MeshBuildProfile.from_name(profile)

        cleanup = bool(overrides.pop("cleanup", profile_obj.cleanup))
        optimization_iterations = int(
            overrides.pop("optimization_iterations", profile_obj.optimization_iterations)
        )
        do_optimize = optimization_iterations > 0 if optimize is None else bool(optimize)

        cleanup_kwargs = profile_obj.cleanup_kwargs()
        build_kwargs = {
            "cleanup": cleanup,
            "optimization_iterations": optimization_iterations,
            "damping": profile_obj.damping,
            "min_move": profile_obj.min_move,
            "verbose": verbose,
        }
        build_kwargs.update(cleanup_kwargs)
        build_kwargs.update(overrides)

        resolved_protected_labels = self._resolve_protected_labels(
            protected_labels=protected_labels,
            protect_sources=protect_sources if do_optimize else None,
        )

        if do_optimize:
            result = self.build_optimized(
                protected_labels=resolved_protected_labels,
                **build_kwargs,
            )
            result["profile"] = profile_obj.name
            result["mode"] = "optimized"
            return result

        if cleanup:
            self.clean_geometry(**cleanup_kwargs | {k: build_kwargs[k] for k in cleanup_kwargs if k in build_kwargs})
        else:
            self.prepare()

        self._optimization_points = []
        self._build_triangle_mesh(verbose=verbose)
        quality_report = triangle_quality_report(self, include_voronoi=True, validate_voronoi=True)
        self._last_quality_report = quality_report
        return {
            "cleanup": self._last_cleanup_report if cleanup else None,
            "optimization": None,
            "quality": quality_report,
            "profile": profile_obj.name,
            "mode": "cleanup_only" if cleanup else "baseline",
            "protected_labels": resolved_protected_labels,
        }
