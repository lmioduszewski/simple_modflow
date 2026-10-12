"""A 3-D surface drawn across its domain's concavities, and the three clip defects under it.

Ledger 314. Reported as: `InterpolatedSurface(xs, ys, zs, vor).show()` "looks clipped
to a convex hull but not the actual domain", while `plot_heatmap()` was clipped.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/surface_domain_mask.py

THE SETUP. A real triangulated grid -- `GridSpec.voronoi`, `max_area=5000` -- on an L,
built with `tessellation="flopy"`: FloPy's own cells, which these numbers were measured on
before ledger 316 repaired them (the default grid has no hole at the corner now, and
section 4 then prints none). The L:
1000 x 1000 ft with the north-east 600 x 600 removed, so the domain is NOT convex. On
a convex domain the hull of the cell centres sits inside the footprint and none of
this is visible, which is why no test on the usual rectangular fixtures saw it. The
field is a plane, so every interpolator reproduces it and only coverage differs.
Resolution 300 per axis.

WHAT IT MEASURES
  1. `.show()` drew `.fig` = `clipped_fig(clip=False)`: never masked. The edge it
     showed is `griddata`'s -- NaN outside the convex hull of the points.
     `clipped_fig(clip=False)` still draws exactly that, so it is the "before".
  2. The raster `transform` put its ORIGIN (a pixel corner, to `from_origin`) on the
     first cell centre, where `xy_meshgrid` puts the first pixel's CENTRE: every
     pixel sat half a pixel east and south of where the surface drew it. The count
     is of drawn pixels a mask tested there gets wrong.
  3. The clip filled with 0 and `clipped_fig`/`save_raster` then turned every 0 to
     NaN. Linear interpolation reproduces an exactly-zero region exactly (the default
     cubic `griddata` produced no exact zero on this field), so the loss is counted
     on `interpolator="linearND"`.
  4. The grid's own footprint has a hole at the L's inner corner -- the cells leave
     a gap there -- and a masked surface shows it. That is the mask being faithful,
     not a defect of it. The hole is FloPy's ``tri2vor`` mis-sorting the corner cell:
     ledger 316, ``voronoi_tiling_gaps.py``, which the default grid now repairs.
"""

from __future__ import annotations

import atexit
import pathlib
import shutil
import tempfile
import warnings

import geopandas as gpd
import numpy as np
import shapely
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from shapely.geometry import Polygon, mapping

import myflopy as mf
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

RESOLUTION = 300


def _l_grid(root: pathlib.Path):
    shape = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"name": ["domain"]}, geometry=[shape], crs=2927).to_file(boundary)
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=5_000.0, tessellation="flopy"
    ).resolve(workspace=root / "grid")


def main() -> None:
    warnings.simplefilter("ignore")
    root = pathlib.Path(tempfile.mkdtemp())
    atexit.register(shutil.rmtree, root, ignore_errors=True)  # /tmp is a quota'd tmpfs (ledger 460)
    grid = _l_grid(root)
    xs = np.asarray(grid.centroids_x, dtype=float)
    ys = np.asarray(grid.centroids_y, dtype=float)
    zs = 100.0 + 0.01 * xs + 0.02 * ys
    domain = grid.gdf_vorPolys.union_all()
    hull = shapely.MultiPoint(np.column_stack([xs, ys])).convex_hull
    print(f"grid: {len(xs)} cells on the L, resolution {RESOLUTION}")
    print(f"domain {domain.area:,.0f} ft2; hull of the cell centres {hull.area:,.0f} ft2, "
          f"{hull.difference(domain).area:,.0f} ft2 of it outside the domain")

    surface = InterpolatedSurface(xs, ys, zs, grid, resolution=RESOLUTION)
    grid_x, grid_y = surface.xy_meshgrid
    inside = shapely.contains_xy(domain, grid_x, grid_y)

    # 1. What .show() drew, and what it draws now.
    before = np.asarray(surface.clipped_fig(clip=False).data[0].z, dtype=float)
    after = np.asarray(surface.fig.data[0].z, dtype=float)
    drawn = np.isfinite(before)
    print(f"1. unmasked (what .show() drew): {(drawn & ~inside).sum():,} of {drawn.sum():,} "
          f"drawn pixels outside the domain ({(drawn & ~inside).sum() / drawn.sum():.1%})")
    print(f"   .fig now: {(np.isfinite(after) & ~inside).sum()} outside; "
          f"{(~np.isfinite(after) & inside).sum()} pixels inside the domain have no value "
          "(outside the hull -- the rim ledger 248 describes; masking does not fill it)")

    # 2. The half-pixel registration.
    xres = (xs.max() - xs.min()) / (RESOLUTION - 1)
    yres = (ys.max() - ys.min()) / (RESOLUTION - 1)
    old = from_origin(xs.min(), ys.max(), xres, yres)
    shape = (RESOLUTION, RESOLUTION)
    old_mask = geometry_mask([mapping(domain)], out_shape=shape, transform=old, invert=True)
    new_mask = geometry_mask([mapping(domain)], out_shape=shape, transform=surface.transform, invert=True)
    x0, y0 = old * (0.5, 0.5)
    print(f"2. old transform: pixel (0, 0) at ({x0 - grid_x[0, 0]:+.2f}, {y0 - grid_y[0, 0]:+.2f}) "
          f"ft from where it is drawn, a {xres:.2f} ft pixel")
    print(f"   drawn pixels a mask tested there gets wrong: {(drawn & (old_mask != inside)).sum()} "
          f"(dropped inside {(drawn & inside & ~old_mask).sum()}, kept outside "
          f"{(drawn & ~inside & old_mask).sum()}); on the registered transform: "
          f"{(drawn & (new_mask != inside)).sum()}")

    # 3. The 0 sentinel.
    zero_west = np.where(xs < 500.0, 0.0, (xs - 500.0) * 0.01)
    linear = InterpolatedSurface(xs, ys, zero_west, grid, resolution=RESOLUTION, interpolator="linearND")
    raw = np.asarray(linear.surface, dtype=float)
    zeros = (raw == 0.0) & inside
    kept = np.asarray(linear.clipped_fig(clip=True).data[0].z, dtype=float)[zeros]
    print(f"3. linearND, zero west of x=500: {zeros.sum():,} in-domain pixels are exactly 0; "
          f"the old `== 0 -> NaN` erased every one, now {(kept == 0.0).sum():,} kept")

    # 4. The footprint's own hole.
    holes = [Polygon(ring) for ring in getattr(domain, "interiors", [])]
    print(f"4. the grid's footprint has {len(holes)} hole(s), "
          f"{sum(h.area for h in holes):,.1f} ft2 in all, at "
          f"{[tuple(round(v) for v in h.centroid.coords[0]) for h in holes]}")


if __name__ == "__main__":
    main()
