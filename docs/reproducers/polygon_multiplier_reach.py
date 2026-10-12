"""How many cells a polygon MULTIPLIER reaches, per rule (ledgers 541, 551).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/polygon_multiplier_reach.py

No MF6. Under a second.

The model: ``rectangular_voronoi`` meshes of 6x6, 10x10 and 20x20 cells over a square
footprint, and one "pit" box covering the south-west 40% x 50% of it -- the shape of
``mf.RasterValues(...) * mf.Polygons("pit.gpkg", value=0.1, default=1.0)``.

With ``default=`` given (the multiplier's case) the part of a cell no polygon covers
COMPETES as the default (ledger 551): under ``rule="majority"`` (the
default) the pit takes a cell only when it covers MORE than half of it, and under
``rule="mean"`` a partly covered cell blends 0.1 and 1.0 by area. Without ``default=``
the remainder is no competitor and ``majority`` reaches every cell the pit OVERLAPS at
all -- what ledger 541 measured, printed here as ``no default``. Printed per mesh: the
cells the pit touches, the cells more than half inside, the cells each rule multiplies
by 0.1 (to 1e-12) with ``default=1.0`` (``mean`` also: the cells it blends), the cells
``majority`` reaches without a default, and the area multiplied by 0.1 under
``majority`` over the pit's own area.

Measured 2026-10-12 (``mean``'s count is the cells wholly inside, the only ones it
multiplies by 0.1 unblended; the rest of the cells it touches it blends)::

    6x6    touched 11  >half 6   majority 6   centroid 6   mean 3 (+8 blended)   no default 11  area 0.992x
    10x10  touched 23  >half 15  majority 15  centroid 15  mean 13 (+10 blended) no default 23  area 0.867x
    20x20  touched 76  >half 64  majority 64  centroid 64  mean 56 (+20 blended) no default 76  area 0.991x

Before ledger 551 ``majority`` with ``default=1.0`` multiplied the
``no default`` count -- 11, 23 and 76 cells, 1.865x, 1.495x and 1.215x the pit's area.
"""

from __future__ import annotations

import warnings

import numpy as np
from shapely.geometry import box

from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.modflow.mf6.grid.geometry import polygon_cell_values

warnings.simplefilter("ignore")

for n in (6, 10, 20):
    grid = rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1))
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pit = box(xmin, ymin, xmin + 0.4 * (xmax - xmin), ymin + 0.5 * (ymax - ymin))
    cells = grid.gdf_vorPolys.geometry
    share = np.array([cell.intersection(pit).area / cell.area for cell in cells])
    painted = {
        rule: polygon_cell_values(grid.gdf_vorPolys, [pit], [0.1], rule=rule, default=1.0)[0]
        for rule in ("majority", "centroid", "mean")
    }
    exact = {rule: np.isclose(values, 0.1, rtol=1e-12) for rule, values in painted.items()}
    blended = int((~exact["mean"] & (painted["mean"] < 1.0)).sum())
    no_default = int(polygon_cell_values(grid.gdf_vorPolys, [pit], [0.1])[1].sum())
    multiplied = float(
        (exact["majority"] * cells.area.to_numpy()).sum() / pit.area
    )
    print(
        f"{n}x{n}  touched {int((share > 0).sum())}  >half {int((share > 0.5).sum())}  "
        + "  ".join(f"{rule} {int(hit.sum())}" for rule, hit in exact.items())
        + f" (+{blended} blended)  no default {no_default}  area {multiplied:.3f}x"
    )
