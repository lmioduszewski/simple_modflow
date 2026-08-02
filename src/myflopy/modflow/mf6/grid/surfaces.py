from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.errors import WindowError
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window, from_bounds

from myflopy._logging import get_logger

logger = get_logger(__name__)

# Above this cell count, area-weighted raster sampling is warned as potentially
# slow (the user can switch to method="centroid"). Not an auto-fallback.
_AREA_SAMPLE_WARN_NCPL = 50_000


def _centroid_sample(src, vor_crs, xs, ys, raster_path) -> np.ndarray:
    """Sample the raster at each cell centroid (point sampling)."""

    # Reproject centroids from the grid CRS to the raster CRS so we sample the
    # correct pixels when the two differ -- otherwise the coordinates are read
    # in the raster's CRS and sampling returns silent garbage.
    if vor_crs is not None and src.crs is not None:
        sample_xs, sample_ys = warp_transform(vor_crs, src.crs, list(xs), list(ys))
    else:
        warnings.warn(
            f"CRS missing (grid={vor_crs}, raster={src.crs}); sampling "
            f"'{raster_path}' without reprojection.",
            stacklevel=2,
        )
        sample_xs, sample_ys = list(xs), list(ys)
    sampled = np.array(
        [val[0] for val in src.sample(zip(sample_xs, sample_ys, strict=False))], dtype=float
    )
    # rasterio.sample() yields the nodata value (e.g. -9999) for nodata and
    # out-of-bounds points, never None, so convert those to NaN explicitly.
    if src.nodata is not None:
        sampled[sampled == src.nodata] = np.nan
    return sampled


def _area_weighted_sample(vor, src, vor_crs, xs, ys, raster_path) -> np.ndarray:
    """Per-cell mean of raster pixels falling within each Voronoi cell polygon.

    Cells that cover no raster pixel (too small relative to the raster, or
    outside coverage) fall back to point-at-centroid sampling so every cell
    still gets a value.
    """

    polys = vor.gdf_vorPolys
    if vor_crs is not None and src.crs is not None:
        polys_r = polys.to_crs(src.crs)
    else:
        warnings.warn(
            f"CRS missing (grid={vor_crs}, raster={src.crs}); area-sampling "
            f"'{raster_path}' without reprojection.",
            stacklevel=2,
        )
        polys_r = polys
    geoms = list(polys_r.geometry.values)
    ncpl = len(geoms)
    mean = np.full(ncpl, np.nan, dtype=float)
    covered = np.zeros(ncpl, dtype=bool)

    # Read only the window covering the cells, clamped to the raster extent.
    # Degenerate (zero-area, e.g. point geometries) or non-overlapping bounds
    # leave every cell uncovered -> centroid fallback below.
    minx, miny, maxx, maxy = polys_r.total_bounds
    win = None
    if np.all(np.isfinite([minx, miny, maxx, maxy])) and maxx > minx and maxy > miny:
        try:
            win = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
            win = win.round_offsets().round_lengths()
            win = win.intersection(Window(0, 0, src.width, src.height))
        except WindowError:
            win = None
    if win is not None and win.width >= 1 and win.height >= 1:
        band = src.read(1, window=win).astype(float)
        if src.nodata is not None:
            band[band == src.nodata] = np.nan
        # Burn cell ids (1..ncpl) onto the raster grid, then average the valid
        # pixels of each cell in one vectorized pass.
        shapes = (
            (geom, i + 1)
            for i, geom in enumerate(geoms)
            if geom is not None and not geom.is_empty
        )
        labels_arr = rasterize(
            shapes,
            out_shape=band.shape,
            transform=src.window_transform(win),
            fill=0,
            dtype="int32",
        )
        flat_labels = labels_arr.ravel()
        flat_vals = band.ravel()
        valid = (flat_labels > 0) & np.isfinite(flat_vals)
        sums = np.bincount(
            flat_labels[valid], weights=flat_vals[valid], minlength=ncpl + 1
        )[1:]
        counts = np.bincount(flat_labels[valid], minlength=ncpl + 1)[1:]
        covered = counts > 0
        mean[covered] = sums[covered] / counts[covered]

    missing = ~covered
    if missing.any():
        mean[missing] = _centroid_sample(
            src, vor_crs, np.asarray(xs)[missing], np.asarray(ys)[missing], raster_path
        )
    return mean


def get_raster_vals_at_centroids(
    vor,
    raster_files: Path | list[Path | int | float] | None = None,
    labels: str | list[str] | list[int] | None = None,
    *,
    method: str = "area",
) -> gpd.GeoDataFrame | None:
    """Sample raster values or constants onto Voronoi cells.

    ``method="area"`` (default) returns the area-weighted mean of the raster
    pixels within each cell polygon; ``method="centroid"`` samples a single
    pixel at each cell centroid (faster, less faithful on coarse layers).
    """
    if raster_files is None:
        logger.warning('no raster files provided; there are no surfaces to sample')
        return None

    if isinstance(raster_files, Path):
        raster_files = [raster_files]
    if labels is None:
        labels = list(range(len(raster_files)))
    if isinstance(labels, str):
        labels = [labels]

    if len(labels) != len(raster_files):
        raise ValueError("Labels length must match raster_files length")
    if method not in ("area", "centroid"):
        raise ValueError(f"Unknown method {method!r}; expected 'area' or 'centroid'.")

    xs, ys = np.array(vor.centroids[0]), np.array(vor.centroids[1])
    centroid_vals = vor.gdf_vorPolys.copy()
    vor_crs = getattr(vor, 'crs', None)

    numeric_vals = {
        labels[i]: val
        for i, val in enumerate(raster_files)
        if isinstance(val, (int, float))
    }
    raster_items = [
        (labels[i], path)
        for i, path in enumerate(raster_files)
        if isinstance(path, Path)
    ]

    if method == "area" and raster_items and len(xs) > _AREA_SAMPLE_WARN_NCPL:
        warnings.warn(
            f"Area-weighted raster sampling on a large grid ({len(xs)} cells) "
            "may be slow; pass method='centroid' for faster point sampling.",
            stacklevel=2,
        )

    for label, raster_path in raster_items:
        logger.info('sampling raster %s', raster_path)
        with rasterio.open(raster_path) as src:
            if method == "area":
                sampled = _area_weighted_sample(
                    vor, src, vor_crs, xs, ys, raster_path
                )
            else:
                sampled = _centroid_sample(src, vor_crs, xs, ys, raster_path)
            centroid_vals[label] = sampled

    for label, val in numeric_vals.items():
        centroid_vals[label] = float(val)

    ordered_cols = ['geometry'] + list(labels)
    centroid_vals = centroid_vals[ordered_cols]
    centroid_vals.geometry = centroid_vals.centroid
    return centroid_vals


def get_gdf_topbtm_multilyr(
    vor,
    rasters: list,
    labels: list[str] | list[int] | None = None,
) -> gpd.GeoDataFrame:
    """
    Build centroid-sampled top and bottom elevations for each model layer.
    """
    if labels is None:
        labels = list(range(len(rasters)))
    return get_raster_vals_at_centroids(vor, raster_files=rasters, labels=labels)


def get_cell_areas(vor) -> list[float]:
    """
    Return polygon areas for each Voronoi cell.
    """
    logger.debug('computing cell areas')
    return [poly.area for poly in vor.gdf_vorPolys['geometry']]


def get_origin_xy(vor) -> tuple[float, float]:
    """
    Get the lower-left origin of the Voronoi vertex set.
    """
    verts = np.asarray(vor.verts)
    xmin = float(np.min(verts[:, 0]))
    ymin = float(np.min(verts[verts[:, 0] == xmin][:, 1]))
    origin_xy = (xmin, ymin)
    vor.origin_xy = origin_xy
    return origin_xy


def get_domain(vor):
    """
    Return the unioned model domain polygon.
    """
    return vor.gdf_vorPolys.union_all()


def get_normal_from_strike_and_dip(strike: int, dip: int) -> np.ndarray:
    """
    Return the normal vector of a plane defined by strike and dip.
    """
    strike_rad = np.radians(strike)
    dip_rad = np.radians(90 - dip)
    return np.array(
        [
            np.sin(dip_rad) * np.sin(strike_rad),
            np.sin(dip_rad) * np.cos(strike_rad),
            np.cos(dip_rad),
        ]
    )


def get_raster_from_strike_dip(
    vor,
    strike: int,
    dip: int,
    known_point: tuple,
    pixel_size: int = 1,
    output_filename: Path = Path.cwd().joinpath('raster.tif'),
) -> gpd.GeoDataFrame:
    """
    Generate a raster for a sloping plane and sample it at Voronoi centroids.
    """
    domain = get_domain(vor)
    min_x, min_y, max_x, max_y = domain.bounds

    min_x = np.floor(min_x / pixel_size) * pixel_size
    max_x = np.ceil(max_x / pixel_size) * pixel_size
    min_y = np.floor(min_y / pixel_size) * pixel_size
    max_y = np.ceil(max_y / pixel_size) * pixel_size

    width = int((max_x - min_x) / pixel_size)
    height = int((max_y - min_y) / pixel_size)

    dip_rad = np.radians(90 - dip)
    normal = get_normal_from_strike_and_dip(strike, dip)
    transform = from_origin(min_x, max_y, pixel_size, pixel_size)

    x_coords = min_x + np.arange(width) * pixel_size
    y_coords = max_y - np.arange(height) * pixel_size
    x_grid, y_grid = np.meshgrid(x_coords, y_coords)

    known_x, known_y, known_elevation = known_point
    point_vectors = np.stack([x_grid - known_x, y_grid - known_y, np.zeros_like(x_grid)], axis=-1)

    norm_normal = np.linalg.norm(normal)
    distance_along_normal = np.dot(point_vectors, normal) / norm_normal
    elevation_data = known_elevation - distance_along_normal * np.cos(dip_rad)

    with rasterio.open(
        output_filename,
        'w',
        driver='GTiff',
        height=height,
        width=width,
        count=1,
        dtype=rasterio.float32,
        crs=vor.crs,
        transform=transform,
    ) as dst:
        dst.write(elevation_data.astype(np.float32), 1)

    return get_raster_vals_at_centroids(vor, [output_filename], ['elev'])
