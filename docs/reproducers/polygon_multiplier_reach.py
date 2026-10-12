"""How many cells a polygon MULTIPLIER reaches, per rule (ledger 541).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/polygon_multiplier_reach.py

No MF6. Under a second.

The model: ``rectangular_voronoi`` meshes of 6x6, 10x10 and 20x20 cells over a square
footprint, and one "pit" box covering the south-west 40% x 50% of it -- the shape of
``mf.RasterValues(...) * mf.Polygons("pit.gpkg", value=0.1, default=1.0)``.

Under ``rule="majority"`` (the default) a cell takes the value of whichever POLYGON
covers most of it; the uncovered rest of the cell is not a competitor, and ``default=``
fills only cells no polygon reaches (ledger 299's contract). So the multiplier reaches
every cell the pit OVERLAPS at all. Printed: the cells each rule multiplies, the cells
more than half inside, and the multiplied area over the pit's own area.

Measured 2026-10-12::

    6x6    touched 11  >half 6   majority 11  centroid 6   mean 11  area 1.865x
    10x10  touched 23  >half 15  majority 23  centroid 15  mean 23  area 1.495x
    20x20  touched 76  >half 64  majority 76  centroid 64  mean 76  area 1.215x
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
    reached = {
        rule: int(polygon_cell_values(grid.gdf_vorPolys, [pit], [0.1], rule=rule)[1].sum())
        for rule in ("majority", "centroid", "mean")
    }
    multiplied = float(((share > 0) * cells.area.to_numpy()).sum() / pit.area)
    print(
        f"{n}x{n}  touched {int((share > 0).sum())}  >half {int((share > 0.5).sum())}  "
        + "  ".join(f"{rule} {count}" for rule, count in reached.items())
        + f"  area {multiplied:.3f}x"
    )
