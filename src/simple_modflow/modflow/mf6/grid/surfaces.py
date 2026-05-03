from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin


def get_raster_vals_at_centroids(
    vor,
    raster_files: Path | list[Path | int | float] | None = None,
    labels: str | list[str] | list[int] | None = None,
) -> gpd.GeoDataFrame | None:
    """
    Sample raster values or constants at Voronoi cell centroids.
    """
    if raster_files is None:
        print('No raster files provided')
        return None

    if isinstance(raster_files, Path):
        raster_files = [raster_files]
    if labels is None:
        labels = list(range(len(raster_files)))
    if isinstance(labels, str):
        labels = [labels]

    if len(labels) != len(raster_files):
        raise ValueError("Labels length must match raster_files length")

    xs, ys = np.array(vor.centroids[0]), np.array(vor.centroids[1])
    centroid_vals = vor.gdf_vorPolys.copy()

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

    for label, raster_path in raster_items:
        print(f'reading raster file {raster_path}')
        with rasterio.open(raster_path) as src:
            sampled = np.array(
                [val[0] if val is not None else np.nan for val in src.sample(zip(xs, ys))]
            )
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
    print('getting cell areas')
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
