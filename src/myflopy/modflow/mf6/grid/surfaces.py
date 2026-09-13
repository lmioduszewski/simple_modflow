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
from rasterio.warp import transform_bounds
from rasterio.windows import Window, from_bounds

from myflopy._logging import get_logger

logger = get_logger(__name__)

# Above this cell count, area-weighted raster sampling is warned as potentially
# slow (the user can switch to method="centroid"). Not an auto-fallback.
_AREA_SAMPLE_WARN_NCPL = 50_000


def _centroid_sample(src, vor_crs, xs, ys, raster_path, nodata=None) -> np.ndarray:
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
    # `nodata` overrides the file's own, for rasters that use a sentinel (often
    # 0) without declaring it in the header.
    fill = src.nodata if nodata is None else nodata
    if fill is not None:
        sampled[sampled == fill] = np.nan
    return sampled


def _area_weighted_sample(vor, src, vor_crs, xs, ys, raster_path, nodata=None) -> np.ndarray:
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
        # Masked BEFORE the per-cell mean: averaging first would blend the
        # sentinel into every cell that straddles the edge, and no later
        # comparison could recover it (half of 0 and 500 is a plausible 250).
        fill = src.nodata if nodata is None else nodata
        if fill is not None:
            band[band == fill] = np.nan
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
        # `nodata` MUST be forwarded here. A cell whose every pixel was masked
        # counts as uncovered, so it lands in this fallback -- and without the
        # override the centroid read hands the sentinel straight back, undoing
        # the masking for exactly the cells that needed it most.
        mean[missing] = _centroid_sample(
            src, vor_crs, np.asarray(xs)[missing], np.asarray(ys)[missing],
            raster_path, nodata=nodata,
        )
    return mean


GAP_MODES = ("warn", "ignore", "nearest", "raise")


def _fill_nearest(values: np.ndarray, xs, ys, missing: np.ndarray) -> np.ndarray:
    """Copy each gap's value from the nearest cell that HAS one, by centroid distance."""

    from scipy.spatial import cKDTree           # heavy; only this mode needs it

    have = np.flatnonzero(np.isfinite(values))
    tree = cKDTree(np.column_stack([np.asarray(xs)[have], np.asarray(ys)[have]]))
    _distance, nearest = tree.query(
        np.column_stack([np.asarray(xs)[missing], np.asarray(ys)[missing]])
    )
    filled = values.copy()
    filled[missing] = values[have[nearest]]
    return filled


def _disjoint(vor, src) -> bool:
    """Whether the raster and the grid do not overlap at all, in the raster's CRS.

    Only asked to decide whether an all-no-data column should be reported as a
    misplaced file. Unknowable without both CRSs, and a transform can fail on a
    projection a point set survives, so either way the answer is "not disjoint" --
    the caller then gives the ordinary per-cell report, which is true regardless.
    """

    grid_crs = getattr(vor, "crs", None)
    if grid_crs is None or src.crs is None:
        return False
    left, bottom, right, top = vor.gdf_vorPolys.total_bounds
    try:
        gl, gb, gr, gt = transform_bounds(grid_crs, src.crs, left, bottom, right, top)
    except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
        # Deliberately broad, and it only ever downgrades an error message. This
        # reaches pyproj (CRSError, ProjError), rasterio (CRSError), and a plain
        # ValueError for a projection that cannot represent the grid's corners --
        # and a bare EPSG string that geopandas accepts is a caller input, not a bug
        # (see `normalize_crs`). Failing here must not hide the real complaint.
        logger.debug("could not compare raster and grid extents", exc_info=True)
        return False
    rl, rb, rr, rt = src.bounds
    return gl >= rr or gr <= rl or gb >= rt or gt <= rb


def _resolve_gaps(values, vor, src, xs, ys, *, label, raster_path, gaps: str) -> np.ndarray:
    """Say something about the cells a raster gave no value for -- or fill them.

    ``rasterio`` yields the no-data value for a point outside coverage, which this
    module turns into NaN, and the area sampler already falls back to the centroid
    for a cell holding no pixel CENTRE. What was left was a cell whose centroid is
    off the raster's data as well: it came back NaN, and the frame was handed over
    with nothing said. An interpolated surface makes this the normal case rather
    than the odd one -- ``griddata`` fills only the CONVEX HULL of its input points,
    so a raster written from one grid's CELL CENTROIDS has no data in a rim half a
    cell wide around the whole model, and every boundary cell of a DIFFERENT grid
    can land in it. Sampled into ``strt`` or ``gdf_topbtm``, that NaN reaches MF6.

    Gaps are legitimate mid-pipeline -- ``Surface(fill="propagate")`` inherits them
    from the surface above -- so the caller says which it is. ``"ignore"`` is for
    those; the default ``"warn"`` is for a caller reading values off a grid, who
    has no later stage to catch them.
    """

    if gaps not in GAP_MODES:
        raise ValueError(f"Unknown gaps mode {gaps!r}; expected one of {GAP_MODES}.")
    values = np.asarray(values, dtype=float)
    missing = np.flatnonzero(~np.isfinite(values))
    if missing.size == 0 or gaps == "ignore":
        return values

    if missing.size == values.size and _disjoint(vor, src):
        # Nothing landed AND the two do not even overlap: a location or CRS mismatch,
        # reported as two boxes the way an off-grid vector source is (ledger 211).
        # Every cell being no-data where the extents DO overlap is a different thing
        # -- a raster that is blank here -- and falls through to the ordinary report,
        # because claiming a mismatch would send the reader looking for the wrong bug.
        raise ValueError(
            f"Raster '{raster_path}' gave no value for ANY of the {values.size} cells "
            f"of this grid, so column {label!r} would be entirely no-data. The raster "
            f"covers {tuple(round(b, 2) for b in src.bounds)} in {src.crs}; the grid "
            f"covers {tuple(round(b, 2) for b in vor.gdf_vorPolys.total_bounds)} in "
            f"{getattr(vor, 'crs', None)}. Check they are the same place in the same CRS."
        )

    listing = ", ".join(str(cell) for cell in missing[:10])
    if missing.size > 10:
        listing += f", and {missing.size - 10} more"
    if gaps == "raise":
        raise ValueError(
            f"Raster '{raster_path}' gave no value for {missing.size} of "
            f"{values.size} cells in column {label!r}: {listing}. Pass "
            "gaps='nearest' to take each one from its nearest valued cell, or "
            "gaps='warn' to receive them as NaN."
        )
    if gaps == "nearest":
        filled = _fill_nearest(values, xs, ys, missing)
        logger.info(
            "column %r: filled %d of %d cells from the nearest valued cell (%s)",
            label, missing.size, values.size, listing,
        )
        return filled

    warnings.warn(
        f"Raster '{raster_path}' gave no value for {missing.size} of {values.size} "
        f"cells in column {label!r}, which are NaN: {listing}. A NaN reaching an MF6 "
        "input makes MF6 report success and return a NaN budget. Pass gaps='nearest' "
        "to take each from its nearest valued cell, gaps='raise' to refuse, or "
        "gaps='ignore' if a later stage fills them.",
        stacklevel=3,
    )
    return values


def get_raster_vals_at_centroids(
    vor,
    raster_files: Path | list[Path | int | float] | None = None,
    labels: str | list[str] | list[int] | None = None,
    *,
    method: str = "area",
    nodata: float | None = None,
    gaps: str = "warn",
) -> gpd.GeoDataFrame | None:
    """Sample raster values or constants onto Voronoi cells.

    ``method="area"`` (default) returns the area-weighted mean of the raster
    pixels within each cell polygon; ``method="centroid"`` samples a single
    pixel at each cell centroid (faster, less faithful on coarse layers).

    ``nodata`` overrides each raster's declared no-data value, for files that use
    a sentinel without recording it in the header. It applies to every raster in
    the call, so pass one raster at a time when they differ.

    ``gaps`` decides what happens to a cell the raster gave no value for, which
    used to be handed back as NaN with nothing said -- see :func:`_resolve_gaps`:

    ``"warn"`` (default)
        Keep the NaN and warn, naming the cells.
    ``"nearest"``
        Take each gap from its nearest valued cell, by centroid distance, and log
        which. This is the fix for sampling a raster interpolated from another
        grid's centroids, whose data stops half a cell inside the model edge.
    ``"raise"``
        Refuse the sample.
    ``"ignore"``
        Keep the NaN silently -- for a caller whose next stage fills them, which is
        what ``Surface(fill="propagate")`` and the pinch-out rules do.

    A raster that covers NO cell of the grid raises whatever the mode, except
    ``"ignore"``: that is a location or CRS mismatch, and the error prints both
    extents and both CRSs.
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
                    vor, src, vor_crs, xs, ys, raster_path, nodata=nodata
                )
            else:
                sampled = _centroid_sample(
                    src, vor_crs, xs, ys, raster_path, nodata=nodata
                )
            centroid_vals[label] = _resolve_gaps(
                sampled, vor, src, xs, ys,
                label=label, raster_path=raster_path, gaps=gaps,
            )

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
