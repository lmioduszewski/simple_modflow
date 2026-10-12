"""Raw raster sampling vs a reconciled LayerStack: where `gdf_topbtm` is wrong.

Ledger 279's deferral, checked in because the figure that first prompted it was
quoted without its setup and could not be reproduced.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/raw_vs_reconciled_surfaces.py

THE MECHANISM. Two routes sample the same rasters and only one reconciles:

  * `grid.rasters = [top, bot]` -> `VoronoiGridPlus.gdf_topbtm` calls
    `get_raster_vals_at_centroids` -- RAW per-cell sampling, no reconciliation.
  * `LayerStack.build()` -- pushes each bottom below the one above, enforces
    `min_thickness`, applies pinch rules.

Where the rasters CROSS -- which real interpolated surfaces do, since
`InterpolatedSurface` fills only the convex hull of its control points (ledger 248)
and a USG export's conditioned tops rise in 17-38% of nodes (ledger 169) -- the two
disagree, and `_layer_elevations` consults the RAW one first.

ALL THREE CONDITIONS ARE REQUIRED, which is why the default path is safe:
  1. the grid carries `rasters=`,
  2. the stack is built with `attach=False`  (`attach=True` is the DEFAULT), and
  3. the rasters actually cross.
"""

from __future__ import annotations

import atexit
import pathlib
import shutil
import tempfile

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext


def _grid():
    return rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1))


def main() -> None:
    tmp = pathlib.Path(tempfile.mkdtemp())
    atexit.register(shutil.rmtree, tmp, ignore_errors=True)  # /tmp is a quota'd tmpfs (ledger 460)
    grid = _grid()
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds

    def raster(name, value_at, n=60):
        path = tmp / name
        xs, ys = np.linspace(xmin, xmax, n), np.linspace(ymax, ymin, n)
        data = np.array([[value_at(x, y) for x in xs] for y in ys], dtype="float32")
        with rasterio.open(
            path, "w", driver="GTiff", height=n, width=n, count=1, dtype="float32",
            crs=grid.crs,
            transform=from_origin(xmin, ymax, (xmax - xmin) / n, (ymax - ymin) / n),
        ) as dst:
            dst.write(data, 1)
        return path

    # A flat top, and a bottom that RISES THROUGH it across the eastern third.
    top = raster("top.tif", lambda x, y: 100.0)
    bot = raster("bot.tif", lambda x, y: 60.0 + 60.0 * (x - xmin) / (xmax - xmin))

    print("top 100.0 ft flat;  bottom ramps 60.0 -> 120.0 ft, so it crosses the top\n")

    print("=== 1. the two tables ===")
    for label, attach in (("attach=False", False), ("attach=True (DEFAULT)", True)):
        g = _grid()
        g.rasters = [top, bot]
        stack = mf.LayerStack(top=mf.Raster(top))
        stack.add("unit", bottom=mf.Raster(bot))
        result = stack.build(g, attach=attach)
        table = np.asarray(g.gdf_topbtm.iloc[:, 2], dtype=float)
        reconciled = np.asarray(result.botm[0], dtype=float)
        differing = np.abs(table - reconciled) > 1e-6
        print(
            f"  {label:22} differ {differing.sum():2d}/{g.ncpl}"
            f"  max |diff| {np.abs(table - reconciled).max():7.3f} ft"
            f"  botm above top: {(table > 100.0).sum():2d} cells"
        )

    print("\n=== 2. does it reach a boundary condition? (attach=False) ===")
    g = _grid()
    stack = mf.LayerStack(top=mf.Raster(top))
    stack.add("unit", bottom=mf.Raster(bot))
    built = stack.build(g, attach=False)
    g.rasters = [top, bot]

    drn_path = tmp / "drn.gpkg"
    gpd.GeoDataFrame(
        {"name": ["d"], "layer": [1], "cond": [10.0]},
        geometry=[LineString([(xmin, (ymin + ymax) / 2), (xmax, (ymin + ymax) / 2)])],
        crs=g.crs,
    ).to_file(drn_path, layer="drn", driver="GPKG")

    drn = mf.drn.gpkg(
        drn_path, layer="drn", nper=1,
        elevation=mf.CellSurfaceOffset("cell_bottom", offset=0.0),
        conductance=mf.Spread("cond"),
        context=ModelContext(grid=g, domain=built.idomain),
    )
    records = drn.options["stress_period_data"][0]
    cells = [int(record[0][1]) for record in records]
    resolved = np.array([float(record[1]) for record in records])
    from_disv = np.asarray(built.botm[0], dtype=float)[cells]
    from_table = np.asarray(g.gdf_topbtm.iloc[:, 2], dtype=float)[cells]

    print(f"  {len(records)} DRN records")
    print(f"    vs the DISV's own botm : max |diff| {np.abs(resolved - from_disv).max():7.3f} ft")
    print(f"    vs the RAW gdf_topbtm  : max |diff| {np.abs(resolved - from_table).max():7.3f} ft")
    print(f"    -> the BC followed the RAW table, and placed "
          f"{(resolved > 100.0).sum()} of {len(resolved)} inverts ABOVE the model top")


if __name__ == "__main__":
    main()
