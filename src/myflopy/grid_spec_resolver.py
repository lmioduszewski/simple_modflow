"""Resolve durable grid specs into existing Triangle/Voronoi grid objects."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import shapely as shp

from myflopy.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.sources import (
    DataSourceSpec,
    GeoPackageSourceSpec,
    ShapeSource,
    TableSource,
)
from myflopy.specs import GridSpec


_BUILD_OPTION_KEYS = {
    "cleanup",
    "damping",
    "max_optimization_points",
    "min_feature_area",
    "min_move",
    "optimization_iterations",
    "optimize",
    "profile",
    "protect_sources",
    "protected_labels",
    "resample_domain_boundary",
    "resample_region_sources",
    "simplify_tolerance",
    "snap_tolerance",
    "target_segment_length",
    "verbose",
}


def _source_path(
    source: DataSourceSpec, project_root: Path | str | None = None
) -> Path:
    path = Path(source.path)
    if path.is_absolute() or source.external or project_root is None:
        return path
    return Path(project_root) / path


def _read_source(
    source: DataSourceSpec,
    project_root: Path | str | None = None,
    *,
    target_crs: str | None = None,
) -> gpd.GeoDataFrame:
    path = _source_path(source, project_root)
    if isinstance(source, GeoPackageSourceSpec):
        kwargs = {} if source.layer is None else {"layer": source.layer}
        gdf = gpd.read_file(path, **kwargs)
        if source.query:
            gdf = gdf.query(source.query)
    elif isinstance(source, ShapeSource):
        gdf = gpd.read_file(path)
    elif isinstance(source, TableSource):
        frame = pd.read_csv(path)
        if "geometry" in frame:
            gdf = gpd.GeoDataFrame(
                frame, geometry=gpd.GeoSeries.from_wkt(frame["geometry"])
            )
        else:
            x_field = source.metadata.get("x_field", "x")
            y_field = source.metadata.get("y_field", "y")
            if x_field not in frame or y_field not in frame:
                raise ValueError(
                    f"Table source {path} needs geometry WKT or x/y columns "
                    f"({x_field!r}, {y_field!r}) for grid resolution."
                )
            geometry = gpd.points_from_xy(frame[x_field], frame[y_field])
            gdf = gpd.GeoDataFrame(
                frame, geometry=geometry, crs=source.metadata.get("crs")
            )
    else:
        raise TypeError(f"Unsupported grid source type: {type(source).__name__}")

    source_crs = getattr(source, "crs", None) or source.metadata.get("crs")
    if source_crs is not None and gdf.crs is None:
        gdf = gdf.set_crs(source_crs)
    if (
        target_crs is not None
        and gdf.crs is not None
        and str(gdf.crs) != str(target_crs)
    ):
        gdf = gdf.to_crs(target_crs)
    return gdf


def _source_value(source: DataSourceSpec, key: str, default: Any = None) -> Any:
    fields = getattr(source, "fields", {})
    if isinstance(fields, dict) and key in fields:
        return fields[key]
    return default


def _row_value(row: Any, source: DataSourceSpec, key: str, default: Any = None) -> Any:
    field = _source_value(source, key)
    if field is not None and field in row:
        value = row[field]
        if value is not None and pd.notna(value):
            return value
    return default


def _coalesce(value: Any, default: Any) -> Any:
    return default if value is None else value


def _option(options: dict[str, Any], *names: str, default: Any = None) -> Any:
    for name in names:
        if name in options:
            return options[name]
    return default


def _union_geometry(
    source: DataSourceSpec,
    project_root: Path | str | None,
    *,
    target_crs: str | None,
) -> shp.Geometry:
    gdf = _read_source(source, project_root, target_crs=target_crs)
    if gdf.empty:
        raise ValueError(
            f"Grid source has no features: {_source_path(source, project_root)}"
        )
    return gdf.geometry.union_all()


def _add_refinement_source(
    tri: TriangleGrid,
    source: DataSourceSpec,
    *,
    project_root: Path | str | None,
    options: dict[str, Any],
    target_crs: str | None,
) -> None:
    gdf = _read_source(source, project_root, target_crs=target_crs)
    if gdf.empty:
        return
    default_area = _option(options, "refinement_max_area", "default_refinement_area")
    if default_area is None and _source_value(source, "area") is None:
        raise ValueError(
            "Refinement sources require fields={'area': '<column>'} or "
            "options['refinement_max_area']."
        )
    for index, row in gdf.iterrows():
        geometry = row.geometry
        if geometry is None or geometry.is_empty:
            continue
        label = _row_value(row, source, "label", f"refinement_{index}")
        max_area = _row_value(row, source, "area", default_area)
        priority = int(
            _coalesce(
                _row_value(row, source, "priority"),
                _option(options, "refinement_priority", default=0),
            )
        )
        tri.add_region_polygon(
            geometry,
            max_area=float(max_area),
            label=str(label),
            priority=priority,
            source="refinement",
            buffer=float(_option(options, "refinement_buffer", default=0)),
            simplify_tolerance=_option(options, "refinement_simplify_tolerance"),
            densify_dist=_option(options, "refinement_densify_dist"),
        )


def _add_breakline_source(
    tri: TriangleGrid,
    source: DataSourceSpec,
    *,
    project_root: Path | str | None,
    options: dict[str, Any],
    target_crs: str | None,
) -> None:
    gdf = _read_source(source, project_root, target_crs=target_crs)
    if gdf.empty:
        return
    default_area = _option(options, "breakline_max_area", "default_breakline_area")
    if default_area is None and _source_value(source, "area") is None:
        raise ValueError(
            "Breakline sources require fields={'area': '<column>'} or "
            "options['breakline_max_area']."
        )
    for index, row in gdf.iterrows():
        geometry = row.geometry
        if geometry is None or geometry.is_empty:
            continue
        label = _row_value(row, source, "label", f"breakline_{index}")
        max_area = _row_value(row, source, "area", default_area)
        priority = int(
            _coalesce(
                _row_value(row, source, "priority"),
                _option(options, "breakline_priority", default=1),
            )
        )
        buffer = float(_option(options, "breakline_buffer", "line_buffer", default=10))
        negative_buffer = float(
            _option(options, "breakline_negative_buffer_after_clipping", default=0)
        )
        if negative_buffer > 0:
            raise ValueError(
                "breakline_negative_buffer_after_clipping must be less than or equal to zero"
            )
        tri.add_region_polygon(
            geometry.buffer(buffer),
            max_area=float(max_area),
            label=str(label),
            priority=priority,
            source="line",
            buffer=negative_buffer,
            simplify_tolerance=_option(
                options, "breakline_simplify_tolerance", default=10
            ),
            densify_dist=_option(options, "breakline_densify_dist"),
        )


def _add_point_source(
    tri: TriangleGrid,
    source: DataSourceSpec,
    *,
    project_root: Path | str | None,
    target_crs: str | None,
) -> None:
    gdf = _read_source(source, project_root, target_crs=target_crs)
    points = []
    for geometry in gdf.geometry:
        if geometry is None or geometry.is_empty:
            continue
        if isinstance(geometry, shp.Point):
            points.append((geometry.x, geometry.y))
        elif isinstance(geometry, shp.MultiPoint):
            points.extend((point.x, point.y) for point in geometry.geoms)
        else:
            points.append(
                (geometry.representative_point().x, geometry.representative_point().y)
            )
    if points:
        tri.add_points(points)


def _build_options(options: dict[str, Any]) -> dict[str, Any]:
    build_options = {
        key: value for key, value in options.items() if key in _BUILD_OPTION_KEYS
    }
    if "profile" not in build_options:
        build_options["profile"] = "balanced"
    if isinstance(build_options.get("profile"), str):
        build_options["profile"] = MeshBuildProfile.from_name(build_options["profile"])
    return build_options


def _resolve_voronoi(
    spec: GridSpec,
    *,
    project_root: Path | str | None,
    workspace: Path | str | None,
    build: bool,
    return_triangle: bool,
) -> TriangleGrid | VoronoiGridPlus | tuple[VoronoiGridPlus, TriangleGrid]:
    options = dict(spec.options)
    model_ws = workspace or _option(
        options, "model_ws", "workspace", default=Path.cwd() / "_triangle"
    )
    tri = TriangleGrid(
        angle=float(_option(options, "angle", "min_angle", default=32)),
        region_point_tolerance=_option(options, "region_point_tolerance"),
        model_ws=str(model_ws),
    )
    tri.myflopy_grid_spec = spec

    if spec.boundary is None:
        raise ValueError("Voronoi GridSpec requires a boundary source.")
    tri.set_domain_polygon(
        _union_geometry(spec.boundary, project_root, target_crs=spec.crs),
        label=str(_option(options, "boundary_label", default="domain")),
        max_area=_option(options, "boundary_max_area", "default_cell_area", "max_area"),
        buffer=float(_option(options, "boundary_buffer", default=0)),
        simplify_tolerance=_option(options, "boundary_simplify_tolerance"),
        densify_dist=_option(options, "boundary_densify_dist"),
    )

    if spec.refinement is not None:
        _add_refinement_source(
            tri,
            spec.refinement,
            project_root=project_root,
            options=options,
            target_crs=spec.crs,
        )
    for source in spec.breaklines:
        _add_breakline_source(
            tri,
            source,
            project_root=project_root,
            options=options,
            target_crs=spec.crs,
        )
    for source in spec.points:
        _add_point_source(tri, source, project_root=project_root, target_crs=spec.crs)

    if not build:
        tri.prepare()
        return tri

    build_options = _build_options(options)
    verbose = bool(build_options.pop("verbose", False))
    tri.build_mesh(verbose=verbose, **build_options)
    vor = VoronoiGridPlus(
        tri,
        crs=spec.crs or _option(options, "crs", default="EPSG:2927"),
        name=spec.name,
        qhull_options=_option(options, "qhull_options"),
    )
    vor.myflopy_grid_spec = spec
    return (vor, tri) if return_triangle else vor


def resolve_grid_spec(
    spec: GridSpec,
    *,
    project_root: Path | str | None = None,
    workspace: Path | str | None = None,
    build: bool = True,
    return_triangle: bool = False,
) -> Any:
    """Resolve a :class:`GridSpec` into an existing grid implementation."""

    if spec.method != "voronoi":
        raise NotImplementedError(
            "GridSpec.resolve currently wires generated Voronoi specs. "
            f"Received method={spec.method!r}."
        )
    if spec.engine != "triangle_voronoi_plus":
        raise NotImplementedError(
            "GridSpec.resolve currently supports engine='triangle_voronoi_plus'. "
            f"Received engine={spec.engine!r}."
        )
    return _resolve_voronoi(
        spec,
        project_root=project_root,
        workspace=workspace,
        build=build,
        return_triangle=return_triangle,
    )


__all__ = [
    "resolve_grid_spec",
]
