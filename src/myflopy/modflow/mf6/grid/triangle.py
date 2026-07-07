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

from myflopy.modflow.mf6.grid.geometry_cleanup import cleanup_polygonal_geometry
from myflopy.modflow.mf6.grid.helpers import densify_poly
from myflopy.modflow.mf6.grid.mesh_quality import triangle_quality_report
from myflopy.modflow.mf6.grid.seed_optimization import (
    cvt_relaxed_full_seed_points,
    get_fixed_and_free_vertex_ids,
)
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg


@dataclass(slots=True)
class DomainSpec:
    """Describes the model domain and its optional background cell sizing."""

    geometry: Polygon | MultiPolygon
    label: str = "domain"
    max_area: float | int | None = None
    attribute: int | float = 0


@dataclass(slots=True)
class RegionSpec:
    """Describes a requested refinement region before Triangle inputs are built."""

    geometry: Polygon | MultiPolygon
    max_area: float | int
    label: str | None = None
    priority: int = 0
    attribute: int | float = 0
    point_hint: tuple[float, float] | None = None
    source: str = "geometry"


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


@dataclass(slots=True)
class MeshBuildProfile:
    """A named bundle of mesh cleanup + smoothing settings for triangulation.

    Groups the knobs that control how :meth:`TriangleGrid.build_mesh` cleans input
    geometry and relaxes the triangulation -- snapping/simplifying tolerances,
    minimum feature area, target segment length, and Laplacian-style optimization
    (iterations, damping, minimum move). Rather than tuning each by hand, pick a
    built-in preset with :meth:`from_name` (``"fast"``, ``"balanced"``,
    ``"smooth"``, ``"clean"``, ``"baseline"``) and override individual fields as
    needed. Higher quality (more iterations, shorter segments) costs build time.

    Attributes
    ----------
    name
        Profile label.
    cleanup
        Run geometry cleanup before triangulating.
    snap_tolerance, simplify_tolerance, min_feature_area, target_segment_length
        Geometry-cleanup tolerances applied to the input regions/boundary.
    resample_domain_boundary, resample_region_sources
        Whether/which sources to resample to ``target_segment_length``.
    optimization_iterations, damping, min_move
        Mesh-smoothing controls applied after triangulation.
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
        """Return one of the built-in mesh build presets by name."""
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

    Parameters
    ----------
    angle
        Minimum interior angle constraint passed to Triangle (higher = better
        shaped, more cells).
    region_point_tolerance
        Tolerance for matching region marker points to regions.
    domain_clip_tolerance
        Tolerance for clipping geometry to the domain boundary.
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
        to the domain extent when the mesh is built.
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
    ) -> Polygon | MultiPolygon:
        """Turn arbitrary region input into a clean polygon/multipolygon for meshing.

        Coerces the geometry, optionally clips it to ``domain``, buffers points/lines
        into areas, simplifies/densifies, and drops empty pieces; raises if nothing
        with positive area remains.
        """

        geometry = cls._coerce_geometry(geometry)

        if domain is not None:
            domain = cls._coerce_geometry(domain)
            geometry = geometry.intersection(domain)

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
            raise ValueError("geometry does not contain any polygonal area after preprocessing")

        if len(polygons) == 1:
            return polygons[0]
        return MultiPolygon(polygons)

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
        return geometry

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

        Parameters
        ----------
        polygon
            Domain geometry as a shapely polygon, multipolygon, or GIS file
            path.
        label
            Name used in previews and diagnostics.
        max_area
            Optional background cell area for the whole domain.
        buffer, simplify_tolerance, densify_dist
            Preprocessing controls applied before the domain is stored.
        """
        geometry = self._polygonize_geometry(
            polygon,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
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
        """Define the domain as an axis-aligned rectangle."""
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
        """Define the domain as an approximate circular polygon."""
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
        """Load the model domain from GIS data and register it."""
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

        Parameters
        ----------
        polygon
            Region geometry or path to polygonal GIS data.
        max_area
            Maximum triangle area requested inside the region.
        label
            Optional region name used in diagnostics and optimization
            protection.
        priority
            Higher-priority regions claim overlapping area before lower
            priority regions during :meth:`prepare`.
        point_hint
            Optional interior point to try first when selecting Triangle's
            region point.
        source
            Short source label used by cleanup and feature-protection logic.
        """
        clip_domain = self.domain_geometry if domain is None else domain
        geometry = self._polygonize_geometry(
            polygon,
            domain=clip_domain,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
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
        """Add a rectangular refinement region."""
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
        """Add a circular refinement region."""
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
        """Load one or more refinement regions from GIS data."""
        geometry = self._polygonize_geometry(
            shp_gpkg,
            domain=self.domain_geometry,
            buffer=buffer,
            simplify_tolerance=simplify_tolerance,
            densify_dist=densify_dist,
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
        """Add refinement around a linear feature by buffering it into polygons."""
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
        """Backward-compatible helper for adding a circular polygon feature."""
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
        """Backward-compatible helper for adding a rectangular polygon feature."""
        polygon = self._rectangle_polygon(x_dist=x_dist, y_dist=y_dist, origin=origin)
        if return_only:
            return polygon

        if self.domain_spec is None:
            self.set_domain_polygon(polygon, max_area=max_area)
        else:
            self.add_region_polygon(polygon, max_area=max_area if max_area is not None else 1, source="rectangle")

        return polygon

    def add_region(self, point, attribute=0, maximum_area=None):
        """Add a raw Triangle region point and optional area constraint."""
        point = GeoSpatialUtil(point, shapetype="point").points
        self._manual_regions.append([point, attribute, maximum_area])

    def add_regions(self, points, attributes=None, maximum_areas=None):
        """Add multiple raw Triangle region points in one call."""
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
        """Legacy helper that adds one or more polygon regions from GIS data."""
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
        """Add fixed point constraints to the Triangle node list."""
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

        Returns
        -------
        dict
            Small report describing the cleanup settings plus the vertex and
            polygon counts before and after cleanup.
        """
        self.validate_setup()

        domain_geometry, domain_stats = cleanup_polygonal_geometry(
            self.domain_spec.geometry,
            snap_tolerance=snap_tolerance,
            simplify_tolerance=simplify_tolerance,
            min_feature_area=min_feature_area,
            target_segment_length=target_segment_length if resample_domain_boundary else None,
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
            geometry, stats = cleanup_polygonal_geometry(
                region.geometry,
                snap_tolerance=snap_tolerance,
                simplify_tolerance=simplify_tolerance,
                min_feature_area=min_feature_area,
                target_segment_length=(
                    target_segment_length
                    if target_segment_length and region.source in set(resample_region_sources)
                    else None
                ),
            )
            cleaned_specs.append(
                RegionSpec(
                    geometry=geometry,
                    max_area=region.max_area,
                    label=region.label,
                    priority=region.priority,
                    attribute=region.attribute,
                    point_hint=region.point_hint,
                    source=region.source,
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
        """Legacy alias for :meth:`add_line_feature`."""
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
        """Generate outward polygons/point clouds with increasing spacing."""
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
        """Add several outward buffer zones with progressively lower priority."""
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
        """Legacy helper that adds concentric buffered regions around polygons."""
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

        for polygon in self._iter_polygons(self.domain_spec.geometry):
            super().add_polygon(polygon)

        for region in self._prepared_regions:
            super().add_polygon(region.geometry)
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
        if node_points:
            self._nodes = np.asarray(node_points, dtype=float)

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
        return super().build(verbose=verbose)

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
        """Prepare Triangle inputs and run a single direct Triangle build."""
        return self._build_triangle_mesh(verbose=verbose)

    def quality_report(self, *, include_voronoi: bool = False):
        """Return triangle quality metrics and optional Voronoi build checks."""
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
        except Exception as exc:
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
            except Exception as exc:
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
        """Build the mesh and then run the experimental optimization pass."""
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

        Parameters
        ----------
        profile
            Either a preset profile name or a :class:`MeshBuildProfile`
            instance.
        protected_labels, protect_sources
            Labels or source tags that should remain fixed during optimization.
        verbose
            When ``True``, print progress updates during cleanup and smoothing.

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
