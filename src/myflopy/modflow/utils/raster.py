from __future__ import annotations
from typing import TYPE_CHECKING, Dict, Sequence, Union, List, Optional, Any

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import rasterio
import shapely as shp
from pathlib import Path
import geopandas as gpd
import pickle
import pandas as pd
from osgeo import gdal, ogr, osr
import numpy as np
import os
from rasterio.features import geometry_mask
from shapely.geometry import LineString, mapping, Point
from dataclasses import dataclass


def _gaussian_kernel_1d(sigma: float, radius: int) -> np.ndarray:
    """A normalized 1-D Gaussian kernel of half-width ``radius`` and spread ``sigma``."""

    if sigma <= 0:
        raise ValueError("sigma must be > 0")
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    k = np.exp(-(x * x) / (2.0 * sigma * sigma))
    k /= k.sum()
    return k


def _convolve1d_reflect(arr: np.ndarray, kernel: np.ndarray, axis: int) -> np.ndarray:
    """Convolve ``arr`` along ``axis`` with ``kernel`` using reflect-padded edges."""

    radius = (kernel.size - 1) // 2
    pad_width = [(0, 0)] * arr.ndim
    pad_width[axis] = (radius, radius)
    padded = np.pad(arr, pad_width, mode="reflect")

    padded = np.moveaxis(padded, axis, -1)
    out = np.empty(padded.shape[:-1] + (padded.shape[-1] - 2 * radius,), dtype=np.float64)

    for i in range(out.shape[-1]):
        window = padded[..., i:i + kernel.size]
        out[..., i] = np.tensordot(window, kernel, axes=([-1], [0]))

    out = np.moveaxis(out, -1, axis)
    return out


def _gaussian_smooth_nodata_aware(
    data: np.ndarray,
    nodata: Optional[float],
    sigma: float,
    radius: Optional[int] = None,
) -> np.ndarray:
    """A separable 2-D Gaussian smooth that ignores nodata cells (normalized by valid weight)."""

    data = data.astype(np.float64, copy=False)

    if radius is None:
        radius = int(np.ceil(3.0 * sigma))
    if radius < 1:
        return data.copy()

    if nodata is None or (isinstance(nodata, float) and np.isnan(nodata)):
        mask = np.ones_like(data, dtype=np.float64)
        data0 = np.where(np.isfinite(data), data, 0.0)
    else:
        valid = (data != nodata) & np.isfinite(data)
        mask = valid.astype(np.float64)
        data0 = np.where(valid, data, 0.0)

    k = _gaussian_kernel_1d(sigma=sigma, radius=radius)

    num = _convolve1d_reflect(data0, k, axis=0)
    num = _convolve1d_reflect(num,  k, axis=1)

    den = _convolve1d_reflect(mask, k, axis=0)
    den = _convolve1d_reflect(den,  k, axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0, num / den, np.nan)

    if nodata is not None and np.isfinite(nodata):
        out = np.where(np.isfinite(out), out, nodata)

    return out


def geotiff_to_contours(
    tif_path: str,
    out_gpkg: str,
    layer_name: str = "contours",
    band_index: int = 1,
    interval: float = 1.0,
    base: float = 0.0,
    nodata: float | None = None,
    attr_name: str = "elev",
    ignore_nodata: bool = True,
    smoothing: Optional[Dict[str, Any]] = None,
):
    r"""
            # Example usage:
        # geotiff_to_contours(
        #     tif_path=r"C:\path\to\dem.tif",
        #     out_gpkg=r"C:\path\to\dem_contours.gpkg",
        #     layer_name="dem_contours",
        #     interval=5.0,
        #     base=0.0,
        #     attr_name="z",
        #     smoothing={"method": "gaussian", "sigma": 1.0},  # sigma in pixels
        # )
    :param tif_path:
    :param out_gpkg:
    :param layer_name:
    :param band_index:
    :param interval:
    :param base:
    :param nodata:
    :param attr_name:
    :param ignore_nodata:
    :param smoothing:
    :return:
    """
    ds = gdal.Open(tif_path, gdal.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f"Could not open: {tif_path}")

    band = ds.GetRasterBand(band_index)

    if nodata is None:
        nodata = band.GetNoDataValue()
    if nodata is not None and ignore_nodata:
        band.SetNoDataValue(nodata)

    contour_band = band
    mem_ds = None

    if smoothing is not None:
        method = str(smoothing.get("method", "gaussian")).lower()
        if method != "gaussian":
            raise ValueError(f"Unsupported smoothing method: {method!r} (only 'gaussian' supported)")

        sigma = float(smoothing.get("sigma", 0.0))
        if sigma <= 0:
            raise ValueError("For gaussian smoothing, provide smoothing['sigma'] > 0 (pixels).")

        radius = smoothing.get("radius", None)
        if radius is not None:
            radius = int(radius)

        arr = band.ReadAsArray()
        if arr is None:
            raise RuntimeError("Failed to read raster band as array.")

        smoothed = _gaussian_smooth_nodata_aware(
            arr,
            nodata if ignore_nodata else None,
            sigma=sigma,
            radius=radius,
        )

        mem_ds = gdal.GetDriverByName("MEM").Create(
            "", ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32
        )
        mem_ds.SetGeoTransform(ds.GetGeoTransform())
        mem_ds.SetProjection(ds.GetProjection())

        mem_band = mem_ds.GetRasterBand(1)
        if nodata is not None and ignore_nodata:
            mem_band.SetNoDataValue(nodata)

        mem_band.WriteArray(smoothed.astype(np.float32))
        mem_band.FlushCache()
        contour_band = mem_band

    if os.path.exists(out_gpkg):
        os.remove(out_gpkg)

    drv = ogr.GetDriverByName("GPKG")
    out_ds = drv.CreateDataSource(out_gpkg)
    if out_ds is None:
        raise RuntimeError(f"Could not create: {out_gpkg}")

    wkt = ds.GetProjection()
    srs = osr.SpatialReference()
    if wkt:
        srs.ImportFromWkt(wkt)

    layer = out_ds.CreateLayer(layer_name, srs=srs, geom_type=ogr.wkbLineString)
    if layer is None:
        raise RuntimeError("Could not create output layer")

    layer.CreateField(ogr.FieldDefn(attr_name, ogr.OFTReal))

    # Correct way for your binding: pass everything through options
    options = [
        "ELEV_FIELD=0",
        f"IGNORE_NODATA={'YES' if ignore_nodata else 'NO'}",
        f"LEVEL_INTERVAL={float(interval)}",
        f"LEVEL_BASE={float(base)}",
    ]
    if nodata is not None:
        options.append(f"NODATA={nodata}")

    gdal.ContourGenerateEx(contour_band, layer, options=options)

    layer = None
    out_ds = None
    mem_ds = None
    ds = None

    return out_gpkg


def sample_raster_cells_crossed_by_lines(
    lines: Sequence[LineString],
    raster_path: str,
    *,
    band: int = 1,
    nodata_to_nan: bool = True,
    keep_na: bool = True,
) -> Dict[int, pd.Series]:

    """
    Samples raster cell values crossed by given lines and maps them to their projected distance
    along the lines. This function intersects the lines with the raster and retrieves the pixel
    values, associating them with the distance along the line. Optionally handles nodata values
    and allows filtering of missing data.

    :param lines: Input lines to be intersected with the raster.
    :type lines: Sequence[LineString]
    :param raster_path: Path to the raster file to be sampled.
    :type raster_path: str
    :param band: Band number of the raster to read data from. Defaults to 1.
    :type band: int
    :param nodata_to_nan: Flag to convert nodata values in the raster to NaN. Defaults to True.
    :type nodata_to_nan: bool
    :param keep_na: Flag to decide whether to keep NaN values in the output. Defaults to True.
    :type keep_na: bool
    :return: A dictionary where keys are indices of the lines, and values are pandas Series
        representing the raster values along the line distance.
    :rtype: Dict[int, pd.Series]
    """
    out: Dict[int, pd.Series] = {}

    with rasterio.open(raster_path) as ds:
        arr = ds.read(band)
        nod = ds.nodata

        for i, line in enumerate(lines):
            # mask True where pixel is intersected by the line
            m = geometry_mask([mapping(line)], out_shape=arr.shape, transform=ds.transform, invert=True)
            rows, cols = np.where(m)

            if len(rows) == 0:
                out[i] = pd.Series(
                    dtype=float,
                    index=pd.Index([], name="dist"),
                    name=raster_path,
                )
                continue

            # pixel center coordinates for intersected pixels
            xs, ys = rasterio.transform.xy(ds.transform, rows, cols, offset="center")
            vals = arr[rows, cols].astype(float)

            if nodata_to_nan and nod is not None:
                vals[vals == nod] = np.nan

            # distance along line: project pixel centers onto line
            dists = np.array(
                [line.project(Point(x, y)) for x, y in zip(xs, ys)],
                dtype=float,
            )

            # sort along line
            order = np.argsort(dists)
            dists = dists[order]
            vals = vals[order]

            s = pd.Series(vals, index=pd.Index(dists, name="dist"), name=raster_path)

            if not keep_na:
                s = s.dropna()

            out[i] = s

    return out


class RasterData:

    def __init__(self, raster_path: Path, point: shp.Point = None, vor: Vor = None):
        """Wrap a raster for sampling at a point or across a Voronoi grid's cells."""

        self.raster_path = raster_path
        self.point = point
        self.vor = vor

    @property
    def raster_elevs(self):
        """The raster's cell values (elevations) read from the file."""

        with rasterio.open(self.raster_path) as src:
            elevations = src.read(1)
        return elevations

    def sample_raster(self):
        """
        Returns the raster values for the given iterable of Shapely points.
        """

        if self.point:
            print(f'sampling raster from single point {self.point}')
            points = [self.point]
        elif self.vor:
            print(f'sampling raster at all centroids in voronoi grid')
            points = self.vor.gdf_vorPolys.centroid.to_list()
        else:
            raise ValueError('No valid sample points provided')

        with rasterio.open(self.raster_path) as src:
            # Convert Shapely points to (x, y) coordinate pairs
            coords = [(point.x, point.y) for point in points]

            # rasterio.sample returns an iterator of arrays (one per band)
            # For a single-band raster, each returned array has 1 value: e.g. array([value])
            # Extract that single value with val[0].
            values = [val[0] for val in src.sample(coords)]

        return pd.Series(values)


if __name__ == '__main__':

    prism_raster = Path(r"C:\Users\lukem\mf6\Cumberland general\PRISM_ppt_30yr_normal_800mM4_annual_bil\PRISM_ppt_30yr_normal_annual_inches_EPSG_2926.tif")
    weather_stn = Path(r"C:\Users\lukem\mf6\Cumberland general\landsburg.gpkg")
    vor_path_v2c = Path(r"C:\Users\lukem\mf6\Cumberland general\cumberland_v2c.vor")

    with open(vor_path_v2c, 'rb') as file:
        vor: Vor = pickle.load(file)
    point = gpd.read_file(weather_stn).geometry.iloc[0]

    rs = RasterData(prism_raster, vor=vor)
    print(rs.sample_raster())

