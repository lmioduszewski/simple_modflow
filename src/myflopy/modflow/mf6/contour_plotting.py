"""Quick contour plotting helpers for cell-centered MF6 data."""

from __future__ import annotations
from myflopy.viz import mpl_axes

from collections.abc import Sequence

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.figure import Figure
from scipy.interpolate import CloughTocher2DInterpolator
from shapely.geometry import GeometryCollection, LineString, MultiLineString, Point


def _matplotlib_cmap(name):
    """Return a Matplotlib colormap, accepting Plotly-style capitalization."""

    if name is None:
        return plt.get_cmap("viridis")
    try:
        return plt.get_cmap(name)
    except ValueError:
        return plt.get_cmap(str(name).lower())


def _cell_center_xy(vor) -> tuple[np.ndarray, np.ndarray]:
    """Return one x/y centroid coordinate pair per Voronoi cell."""

    if hasattr(vor, "centroids_x") and hasattr(vor, "centroids_y"):
        return np.asarray(vor.centroids_x, dtype=float), np.asarray(vor.centroids_y, dtype=float)
    centroids = vor.gdf_vorPolys.geometry.centroid
    return centroids.x.to_numpy(dtype=float), centroids.y.to_numpy(dtype=float)


def _as_cell_values(values, *, ncpl: int, label: str) -> np.ndarray:
    """Return one numeric value per cell."""

    arr = np.asarray(values, dtype=float).squeeze()
    if arr.size != int(ncpl):
        raise ValueError(f"{label} must provide one value per cell; got shape={arr.shape}, ncpl={ncpl}.")
    return arr.reshape(int(ncpl))


def _can_contour(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> bool:
    """Whether contouring is possible: at least 3 finite points and 2 distinct z values."""

    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if int(finite.sum()) < 3:
        return False
    return bool(np.unique(z[finite]).size >= 2)


def _resolve_contour_levels(levels: int | float | Sequence[float], values: np.ndarray) -> np.ndarray:
    """Resolve contour levels, treating scalar values as contour intervals."""

    if isinstance(levels, bool):
        raise ValueError("contour levels must be an interval or a sequence of explicit levels.")
    if isinstance(levels, (int, float, np.integer, np.floating)):
        interval = float(levels)
        if interval <= 0:
            raise ValueError("contour interval must be greater than zero.")
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return np.asarray([], dtype=float)
        start = np.floor(float(np.nanmin(finite)) / interval) * interval
        stop = np.ceil(float(np.nanmax(finite)) / interval) * interval
        count = int(round((stop - start) / interval)) + 1
        resolved = start + (np.arange(max(count, 1), dtype=float) * interval)
        return np.round(resolved, decimals=10)
    return np.asarray(list(levels), dtype=float)


def _segments_from_contour_set(contour_set):
    """Flatten a matplotlib contour set into ``{level, x, y}`` polyline dicts (2+ points each)."""

    segments = []
    for level, level_segments in zip(contour_set.levels, contour_set.allsegs, strict=False):
        for segment in level_segments:
            if len(segment) < 2:
                continue
            segments.append(
                {
                    "level": float(level),
                    "x": np.asarray(segment[:, 0], dtype=float),
                    "y": np.asarray(segment[:, 1], dtype=float),
                }
            )
    return segments


def _iter_lines(geometry):
    """Yield each ``LineString`` in a line/multiline/collection geometry (skipping empties)."""

    if geometry is None or geometry.is_empty:
        return
    if isinstance(geometry, LineString):
        yield geometry
        return
    if isinstance(geometry, MultiLineString | GeometryCollection):
        for part in geometry.geoms:
            yield from _iter_lines(part)


def _clip_contour_segments(segments, clip_geometry):
    """Clip each contour polyline to ``clip_geometry``, dropping degenerate pieces (no-op if ``None``)."""

    if clip_geometry is None:
        return segments
    clipped = []
    for item in segments:
        line = LineString(zip(item["x"], item["y"], strict=False))
        for clipped_line in _iter_lines(line.intersection(clip_geometry)):
            if clipped_line.length <= 0:
                continue
            xs, ys = clipped_line.xy
            if len(xs) < 2:
                continue
            clipped.append(
                {
                    "level": item["level"],
                    "x": np.asarray(xs, dtype=float),
                    "y": np.asarray(ys, dtype=float),
                }
            )
    return clipped


def _filter_points_to_clip(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    clip_geometry,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Drop contour input points outside the active clipping geometry."""

    if clip_geometry is None:
        return x, y, z
    keep = np.asarray([clip_geometry.covers(Point(px, py)) for px, py in zip(x, y, strict=False)], dtype=bool)
    return x[keep], y[keep], z[keep]


def _linear_contour_segments(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    levels: np.ndarray,
):
    """Contour segments via linear triangulation of the scattered ``(x, y, z)`` points."""

    fig = Figure()
    ax = fig.subplots()
    triangulation = mtri.Triangulation(x, y)
    contour_set = ax.tricontour(triangulation, z, levels=levels)
    return _segments_from_contour_set(contour_set)


def _cubic_contour_segments(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    levels: np.ndarray,
    resolution: int,
):
    """Contour segments via cubic interpolation onto a ``resolution``-square grid, then contouring."""

    xmin, xmax = float(np.nanmin(x)), float(np.nanmax(x))
    ymin, ymax = float(np.nanmin(y)), float(np.nanmax(y))
    if xmin == xmax or ymin == ymax:
        return []
    grid_x, grid_y = np.meshgrid(
        np.linspace(xmin, xmax, int(resolution)),
        np.linspace(ymin, ymax, int(resolution)),
    )
    interpolator = CloughTocher2DInterpolator(np.column_stack((x, y)), z, fill_value=np.nan)
    grid_z = interpolator(grid_x, grid_y)
    if not np.isfinite(grid_z).any():
        return []
    fig = Figure()
    ax = fig.subplots()
    contour_set = ax.contour(grid_x, grid_y, grid_z, levels=levels)
    return _segments_from_contour_set(contour_set)


def contour_line_segments(
    vor,
    values,
    *,
    levels: int | float | Sequence[float] = 10,
    label: str = "value",
    clip_geometry=None,
    resolution: int = 150,
    method: str = "linear",
):
    """Return contour line segments in the grid CRS.

    Each returned dictionary has ``level``, ``x``, and ``y`` keys. Empty output
    means the values were too sparse, flat, or geometrically unsuitable for
    contouring. Scalar ``levels`` values are interpreted as contour intervals.
    Sequence ``levels`` values are interpreted as exact contour values. When
    ``clip_geometry`` is provided, returned lines are clipped before projection
    or display. ``method="linear"`` uses fast triangulated contours;
    ``method="cubic"`` contours a Clough-Tocher interpolated surface.
    """

    x, y = _cell_center_xy(vor)
    z = _as_cell_values(values, ncpl=len(x), label=label)
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if not _can_contour(x, y, z):
        return []
    x = x[finite]
    y = y[finite]
    z = z[finite]
    x, y, z = _filter_points_to_clip(x, y, z, clip_geometry)
    if not _can_contour(x, y, z):
        return []
    try:
        resolved_levels = _resolve_contour_levels(levels, z)
        if resolved_levels.size == 0:
            return []
        normalized_method = str(method).lower()
        if normalized_method in {"linear", "tri", "tricontour"}:
            segments = _linear_contour_segments(x, y, z, levels=resolved_levels)
        elif normalized_method in {"cubic", "clough", "clough_tocher", "cloughtocher"}:
            segments = _cubic_contour_segments(x, y, z, levels=resolved_levels, resolution=int(resolution))
        else:
            raise ValueError("contour method must be 'linear' or 'cubic'.")
        return _clip_contour_segments(segments, clip_geometry)
    except Exception:
        return []


def contour_line_segments_latlon(
    vor,
    values,
    *,
    levels: int | float | Sequence[float] = 10,
    label: str = "value",
    clip_geometry=None,
    resolution: int = 150,
    method: str = "linear",
):
    """Return contour line segments as longitude/latitude arrays."""

    segments = contour_line_segments(
        vor,
        values,
        levels=levels,
        label=label,
        clip_geometry=clip_geometry,
        resolution=resolution,
        method=method,
    )
    if not segments:
        return []
    crs = getattr(vor, "crs", None)
    if crs is None:
        return [{"level": item["level"], "lon": item["x"], "lat": item["y"]} for item in segments]

    lines = [LineString(zip(item["x"], item["y"], strict=False)) for item in segments]
    geo = gpd.GeoSeries(lines, crs=crs).to_crs(epsg=4326)
    latlon_segments = []
    for item, line in zip(segments, geo, strict=False):
        xs, ys = line.xy
        latlon_segments.append(
            {
                "level": item["level"],
                "lon": np.asarray(xs, dtype=float),
                "lat": np.asarray(ys, dtype=float),
            }
        )
    return latlon_segments


def plot_cell_contours(
    vor,
    values,
    *,
    label: str = "value",
    levels: int | float | Sequence[float] = 10,
    filled: bool = False,
    ax=None,
    cmap: str = "viridis",
    colors=None,
    linewidths: float = 1.0,
    show_grid: bool = True,
    show_points: bool = False,
    label_contours: bool = True,
    colorbar: bool = True,
    title: str | None = None,
    grid_color: str = "#666666",
    grid_linewidth: float = 0.4,
    point_size: float = 16.0,
    **kwargs,
):
    """Plot contours from cell-centered values and return a Matplotlib figure.

    The helper uses cell centroids as interpolation points. If the values are too
    sparse or flat for contour lines, it falls back to a colored centroid scatter
    while still returning a figure with attached values for inspection.
    """

    x, y = _cell_center_xy(vor)
    z = _as_cell_values(values, ncpl=len(x), label=label)
    if ax is None:
        fig, ax = mpl_axes(figsize=(7, 6))
    else:
        fig = ax.figure

    if show_grid:
        vor.gdf_vorPolys.boundary.plot(ax=ax, color=grid_color, linewidth=grid_linewidth)

    contour_set = None
    scatter = None
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    cmap_obj = _matplotlib_cmap(cmap)
    if _can_contour(x, y, z):
        triangulation = mtri.Triangulation(x[finite], y[finite])
        resolved_levels = _resolve_contour_levels(levels, z[finite])
        if resolved_levels.size == 0:
            resolved_levels = levels
        if filled:
            contour_set = ax.tricontourf(triangulation, z[finite], levels=resolved_levels, cmap=cmap_obj, **kwargs)
        else:
            line_kwargs = {"linewidths": linewidths, **kwargs}
            if colors is None:
                line_kwargs["cmap"] = cmap_obj
            else:
                line_kwargs["colors"] = colors
            contour_set = ax.tricontour(triangulation, z[finite], levels=resolved_levels, **line_kwargs)
            if label_contours:
                ax.clabel(contour_set, inline=True, fontsize=8)
        if colorbar:
            colorbar_obj = fig.colorbar(contour_set, ax=ax, shrink=0.9)
            colorbar_obj.set_label(label)
    else:
        scatter = ax.scatter(x[finite], y[finite], c=z[finite], cmap=cmap_obj, s=point_size)
        if colorbar:
            colorbar_obj = fig.colorbar(scatter, ax=ax, shrink=0.9)
            colorbar_obj.set_label(label)

    if show_points and finite.any():
        ax.scatter(x[finite], y[finite], s=point_size, facecolor="none", edgecolor="#222222", linewidth=0.6)
    ax.set_aspect("equal")
    ax.set_axis_off()
    if title:
        ax.set_title(title)
    fig.tight_layout()
    fig._myflopy_contour_values = np.asarray(z, dtype=float)
    fig._myflopy_contour_set = contour_set
    fig._myflopy_contour_scatter = scatter
    return fig
