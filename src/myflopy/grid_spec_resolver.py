"""Resolve durable grid specs into existing Triangle/Voronoi grid objects."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
import shapely as shp

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.sources import (
    DataSourceSpec,
    GeoPackageSourceSpec,
    ShapeSource,
    TableSource,
)
from myflopy.specs import GridSpec

logger = get_logger(__name__)

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

_TRIANGLE_OPTION_KEYS = {
    "additional_args",
    "angle",
    "exe_name",
    "maximum_area",
    "nodes",
    "region_point_tolerance",
}


def _source_path(
    source: DataSourceSpec, project_root: Path | str | None = None
) -> Path:
    """Resolve a source's path to absolute, anchoring relative non-external paths to the project root."""

    path = Path(source.path)
    if path.is_absolute() or source.external or project_root is None:
        return path
    return Path(project_root) / path


def _project_path(value: Path | str, project_root: Path | str | None = None) -> Path:
    """Resolve ``value`` to an absolute path, anchoring a relative one to the project root."""

    path = Path(value)
    if path.is_absolute() or project_root is None:
        return path
    return Path(project_root) / path


def _resolved_project_root(
    project_root: Path | str | None,
    *,
    fallback: Path | str | None = None,
) -> Path:
    """The project root to use: the given value, else ``fallback``, else the CWD."""

    if project_root is not None:
        return Path(project_root)
    if fallback is not None:
        return Path(fallback)
    return Path.cwd()


def _path_option(value: Any, project_root: Path | str | None = None) -> Any:
    """Recursively resolve any path-like option value(s) against the project root.

    Sources, paths, and strings become absolute paths; lists/tuples are mapped
    element-wise; anything else is returned unchanged.
    """

    if isinstance(value, DataSourceSpec):
        return _source_path(value, project_root)
    if isinstance(value, Path):
        return _project_path(value, project_root)
    if isinstance(value, str):
        return _project_path(value, project_root)
    if isinstance(value, list):
        return [_path_option(item, project_root) for item in value]
    if isinstance(value, tuple):
        return tuple(_path_option(item, project_root) for item in value)
    return value


def _input_path(value: Any, project_root: Path | str | None = None) -> Path | None:
    """The local path an input value should exist at, or ``None`` for external sources."""

    if isinstance(value, DataSourceSpec):
        if value.external:
            return None
        return _source_path(value, project_root)
    if isinstance(value, (Path, str)):
        return _project_path(value, project_root)
    return None


def _validate_inputs(inputs: tuple[Any, ...], project_root: Path | str | None) -> None:
    """Raise ``FileNotFoundError`` if any local (non-external) grid input is missing."""

    for value in inputs:
        path = _input_path(value, project_root)
        if path is not None and not path.exists():
            raise FileNotFoundError(f"GridSpec input does not exist: {path}")


def _load_builder(script: Path, function: str) -> Any:
    """Import ``function`` from a standalone grid-builder ``script`` and return it.

    Loads the script under a hash-derived module name (its own directory on
    ``sys.path`` during import) and returns the named callable; raises if the
    module cannot load, the function is absent, or it is not callable.
    """

    module_hash = hashlib.sha1(str(script.resolve()).encode()).hexdigest()[:12]
    module_name = f"_myflopy_grid_builder_{module_hash}"
    module_spec = importlib.util.spec_from_file_location(module_name, script)
    if module_spec is None or module_spec.loader is None:
        raise ImportError(f"Could not load grid builder script: {script}")

    module = importlib.util.module_from_spec(module_spec)
    search_paths = [str(script.parent)]
    prior_sys_path = list(sys.path)
    try:
        for path in reversed(search_paths):
            if path not in sys.path:
                sys.path.insert(0, path)
        module_spec.loader.exec_module(module)
    finally:
        sys.path[:] = prior_sys_path

    try:
        builder = getattr(module, function)
    except AttributeError as error:
        raise AttributeError(
            f"Grid builder function {function!r} was not found in {script}."
        ) from error
    if not callable(builder):
        raise TypeError(f"Grid builder {function!r} in {script} is not callable.")
    return builder


def _read_source(
    source: DataSourceSpec,
    project_root: Path | str | None = None,
    *,
    target_crs: str | None = None,
) -> gpd.GeoDataFrame:
    """Read a grid data source into a GeoDataFrame, reprojected to ``target_crs``.

    Handles GeoPackage (with optional layer/query), shapefile, and CSV table
    sources (geometry from a WKT column or x/y columns); applies the source's
    declared CRS when the file has none.
    """

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
    """The column name a source maps to logical ``key`` (via its ``fields`` map), else ``default``."""

    fields = getattr(source, "fields", {})
    if isinstance(fields, dict) and key in fields:
        return fields[key]
    return default


def _row_value(row: Any, source: DataSourceSpec, key: str, default: Any = None) -> Any:
    """One row's value for logical ``key`` (resolved through the source's field map), else ``default``."""

    field = _source_value(source, key)
    if field is not None and field in row:
        value = row[field]
        if value is not None and pd.notna(value):
            return value
    return default


def _coalesce(value: Any, default: Any) -> Any:
    """``value`` unless it is ``None``, in which case ``default``."""

    return default if value is None else value


def _option(options: dict[str, Any], *names: str, default: Any = None) -> Any:
    """The first of ``names`` present in ``options`` (an alias lookup), else ``default``."""

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
    """The unioned geometry of all features in ``source`` (raises if the source is empty)."""

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
    """Add each polygon in a refinement source to the triangle grid as a max-area region.

    Cell-size targets come from the source's ``area`` field or
    ``options['refinement_max_area']`` (one is required); per-feature label and
    priority fall back to option defaults.
    """

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
    """Add each line in a breakline source to the triangle grid as a buffered region.

    Lines are buffered (``breakline_buffer``/``line_buffer``) into refinement
    polygons with an optional negative post-clip buffer; cell-size targets come
    from the source's ``area`` field or ``options['breakline_max_area']``.
    """

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
    """Add each feature of a point source to the triangle grid as a fixed mesh vertex.

    Points and multipoints contribute their coordinates directly; other
    geometries contribute a representative point.
    """

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
    """Filter ``options`` to mesh-build keys, defaulting/resolving the ``profile``.

    A missing profile defaults to ``"balanced"`` and a string profile is
    resolved to a :class:`MeshBuildProfile`.
    """

    build_options = {
        key: value for key, value in options.items() if key in _BUILD_OPTION_KEYS
    }
    if "profile" not in build_options:
        build_options["profile"] = "balanced"
    if isinstance(build_options.get("profile"), str):
        build_options["profile"] = MeshBuildProfile.from_name(build_options["profile"])
    return build_options


def _triangle_options(spec: GridSpec, options: dict[str, Any]) -> dict[str, Any]:
    """Merge the spec's Triangle options with matching top-level options (spec wins)."""

    triangle_options = dict(spec.triangle_options)
    if "min_angle" in options and "angle" not in triangle_options:
        triangle_options["angle"] = options["min_angle"]
    for key in _TRIANGLE_OPTION_KEYS:
        if key in options and key not in triangle_options:
            triangle_options[key] = options[key]
    return triangle_options


def _mesh_options(spec: GridSpec, options: dict[str, Any]) -> dict[str, Any]:
    """Merge the spec's mesh options with matching top-level options, then build-filter them."""

    mesh_options = dict(spec.mesh_options)
    for key in _BUILD_OPTION_KEYS:
        if key in options and key not in mesh_options:
            mesh_options[key] = options[key]
    return _build_options(mesh_options)


def _voronoi_options(
    spec: GridSpec,
    options: dict[str, Any],
    *,
    project_root: Path | str | None,
) -> dict[str, Any]:
    """Merge the spec's Voronoi options with top-level options; default CRS/name, resolve paths."""

    voronoi_options = dict(spec.voronoi_options)
    for key in ("idomain", "idomain_path", "name", "qhull_options", "rasters"):
        if key in options and key not in voronoi_options:
            voronoi_options[key] = options[key]
    voronoi_options.setdefault("crs", spec.crs or options.get("crs", "EPSG:2927"))
    voronoi_options.setdefault("name", spec.name)
    for key in ("idomain_path", "rasters"):
        if key in voronoi_options:
            voronoi_options[key] = _path_option(voronoi_options[key], project_root)
    return voronoi_options


def _resolve_voronoi(
    spec: GridSpec,
    *,
    project_root: Path | str | None,
    workspace: Path | str | None,
    build: bool,
    return_triangle: bool,
) -> TriangleGrid | VoronoiGridPlus | tuple[VoronoiGridPlus, TriangleGrid]:
    """Resolve a Voronoi :class:`GridSpec` into the built grid object(s).

    Builds a :class:`TriangleGrid` from the boundary + refinement/breakline/point
    sources, then (when ``build``) meshes it and wraps it in a
    :class:`VoronoiGridPlus`. With ``build=False`` the prepared (unmeshed)
    triangle grid is returned; ``return_triangle`` additionally returns the
    triangle grid alongside the Voronoi grid.
    """

    options = dict(spec.options)
    triangle_options = _triangle_options(spec, options)
    model_ws = workspace or _option(
        triangle_options,
        "model_ws",
        "workspace",
        default=_option(
            options, "model_ws", "workspace", default=Path.cwd() / "_triangle"
        ),
    )
    triangle_options.pop("model_ws", None)
    triangle_options.pop("workspace", None)
    region_point_tolerance = triangle_options.pop(
        "region_point_tolerance",
        _option(options, "region_point_tolerance"),
    )
    tri = TriangleGrid(
        model_ws=str(model_ws),
        region_point_tolerance=region_point_tolerance,
        **triangle_options,
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

    build_options = _mesh_options(spec, options)
    verbose = bool(build_options.pop("verbose", False))
    tri.build_mesh(verbose=verbose, **build_options)
    voronoi_options = _voronoi_options(spec, options, project_root=project_root)
    vor = VoronoiGridPlus(tri, **voronoi_options)
    vor.myflopy_grid_spec = spec
    return (vor, tri) if return_triangle else vor


def _resolve_python(
    spec: GridSpec,
    *,
    project_root: Path | str | None,
    workspace: Path | str | None,
) -> Any:
    """Resolve a Python :class:`GridSpec` by importing and calling its builder function.

    Validates the script/function and declared inputs exist, then invokes the
    builder to produce the grid object.
    """

    if spec.script is None:
        raise ValueError("Python GridSpec requires a script.")
    if spec.function is None:
        raise ValueError("Python GridSpec requires a function.")

    root = _resolved_project_root(project_root)
    script = _project_path(spec.script, root)
    if not script.exists():
        raise FileNotFoundError(f"Grid builder script does not exist: {script}")
    if not script.is_file():
        raise ValueError(f"Grid builder script is not a file: {script}")

    _validate_inputs(spec.inputs, root)
    grid_workspace = (
        Path(workspace) if workspace is not None else root / "_grid" / spec.name
    )
    grid_workspace.mkdir(parents=True, exist_ok=True)

    builder = _load_builder(script, spec.function)
    grid = builder(project_root=root, workspace=grid_workspace, spec=spec)
    if grid is None:
        raise ValueError(f"Grid builder {spec.function!r} in {script} returned None.")
    try:
        grid.myflopy_grid_spec = spec
    except Exception:  # noqa: BLE001 - see below; `grid` came from user code
        # Deliberately broad, and the asymmetry is the whole argument. `grid` is
        # whatever a USER'S builder script returned, so the set of types a
        # `setattr` on it can raise is genuinely open: AttributeError for
        # __slots__ / a read-only property / a frozen dataclass, TypeError for a
        # class rather than an instance, and -- measured against the pydantic
        # installed here -- a plain ValueError for an undeclared field. All that
        # is lost is a provenance tag; narrowing would trade that for failing
        # the user's whole grid resolution.
        logger.debug(
            "could not tag the grid from %r with its GridSpec", spec.name,
            exc_info=True,
        )
    return grid


def resolve_grid_spec(
    spec: GridSpec,
    *,
    project_root: Path | str | None = None,
    workspace: Path | str | None = None,
    build: bool = True,
    return_triangle: bool = False,
) -> Any:
    """Resolve a :class:`GridSpec` into an existing grid implementation."""

    if spec.method == "python":
        if return_triangle:
            raise ValueError("return_triangle is only supported for Voronoi GridSpec.")
        return _resolve_python(spec, project_root=project_root, workspace=workspace)

    if spec.method != "voronoi":
        raise NotImplementedError(
            "GridSpec.resolve currently wires Python (GridSpec.python) and "
            "generated Voronoi (GridSpec.voronoi) specs only. To use an existing "
            "or structured grid, wrap a built grid with GridSpec.from_object(grid) "
            f"or load one with GridSpec.from_pickle(path). Received method={spec.method!r}."
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
