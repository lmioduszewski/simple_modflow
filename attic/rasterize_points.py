"""rasterize_points — attic'd from modflow/utils/surfaces.py (plan 4.3).

Generic x/y/z-points -> GeoTIFF interpolation helper with zero importers in
src/, tests/, or examples/ at retirement time (2026-07-17).
"""

import numpy as np
import rasterio
from rasterio.transform import from_origin
from scipy.interpolate import RBFInterpolator, griddata
from scipy.spatial import cKDTree

def rasterize_points(
    df: pd.DataFrame,
    xcol="x",
    ycol="y",
    zcol="z",
    cellsize=20.0,
    crs="EPSG:4326",            # set to your projected CRS
    bounds=None,                 # (xmin, ymin, xmax, ymax); if None uses data extent
    method="griddata",                # "idw" | "griddata" | "rbf"
    idw_power=2.0,
    idw_k=16,                    # neighbors for IDW
    idw_radius=None,             # optional search radius (same units as coords)
    rbf_kernel="linear",         # e.g., "linear", "thin_plate_spline"
    rbf_smoothing=0.0,           # >0 smooths
    nodata=np.nan,
    out_path=None
):
    """
    Interpolates a raster from a DataFrame of x, y, z - points.
    Returns (Z, transform) if out_path is None, else saves GeoTIFF and returns out_path.
    df must have columns x, y, z. Coordinates should be in a projected CRS (meters/feet).

    Example usage:

        rasterize_points(
            df,
            cellsize=20.0,
            crs="EPSG:2927",
            method="griddata",   # or "rbf"
            out_path="surface_z.tif"
        )
    """
    # --- prep points ---
    x = df[xcol].to_numpy(dtype=float)
    y = df[ycol].to_numpy(dtype=float)
    z = df[zcol].to_numpy(dtype=float)

    if bounds is None:
        xmin, xmax = np.nanmin(x), np.nanmax(x)
        ymin, ymax = np.nanmin(y), np.nanmax(y)
        # tiny padding to fully include edges
        pad = 1e-9
        xmin -= pad
        ymin -= pad
        xmax += pad
        ymax += pad
    else:
        xmin, ymin, xmax, ymax = bounds

    # build grid (top-left origin)
    xi = np.arange(xmin, xmax + cellsize, cellsize)
    yi = np.arange(ymax, ymin - cellsize, -cellsize)  # descending so row 0 is top
    cols = xi.size
    rows = yi.size
    XI, YI = np.meshgrid(xi, yi)
    grid_pts = np.column_stack([XI.ravel(), YI.ravel()])
    transform = from_origin(xmin, ymax, cellsize, cellsize)

    # --- interpolation methods ---
    if method.lower() == "idw":
        tree = cKDTree(np.column_stack([x, y]))
        # query k neighbors (fast, vectorized). SciPy >=1.6 supports workers=-1
        dist, idx = tree.query(grid_pts, k=min(idw_k, len(x)), workers=-1)

        # handle 1 neighbor vs many neighbors consistently
        if dist.ndim == 1:
            dist = dist[:, None]
            idx = idx[:, None]

        if idw_radius is not None:
            mask = dist > idw_radius
        else:
            mask = np.zeros_like(dist, dtype=bool)

        with np.errstate(divide="ignore"):
            w = 1.0 / np.power(dist, idw_power)
        w[mask] = 0.0
        # exact hits (dist=0) -> set weight=inf so they dominate
        w[np.isinf(w)] = 1e12

        vals = z[idx]
        wsum = w.sum(axis=1)
        out = np.full(grid_pts.shape[0], np.nan, dtype=float)
        good = wsum > 0
        out[good] = (w[good] * vals[good]).sum(axis=1) / wsum[good]

        Z = out.reshape(rows, cols)

    elif method.lower() == "griddata":
        pts = np.column_stack([x, y])
        # linear inside hull
        Z = griddata(pts, z, (XI, YI), method="linear")
        # nearest fill for holes / outside hull
        Zn = griddata(pts, z, (XI, YI), method="nearest")
        if Z is None:
            Z = Zn
        else:
            fill = np.isnan(Z)
            Z[fill] = Zn[fill]

    elif method.lower() == "rbf":
        pts = np.column_stack([x, y])
        rbf = RBFInterpolator(pts, z, kernel=rbf_kernel, smoothing=rbf_smoothing)
        Z = rbf(grid_pts).reshape(rows, cols)
    else:
        raise ValueError("method must be 'idw', 'griddata', or 'rbf'")

    # --- write GeoTIFF or return ---
    if out_path:
        profile = {
            "driver": "GTiff",
            "height": rows,
            "width": cols,
            "count": 1,
            "dtype": "float32",
            "crs": crs,
            "transform": transform,
            "nodata": nodata if np.isnan(nodata) else float(nodata),
            "compress": "deflate",
            "tiled": True,
            "predictor": 3,
        }
        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(Z.astype("float32"), 1)
        return out_path
    else:
        return Z, transform

