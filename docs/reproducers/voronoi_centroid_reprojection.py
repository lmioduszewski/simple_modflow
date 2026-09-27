"""A cell valid where it was drawn and crossing itself once reprojected: an edge of two ulps.

Ledger 316 (the fuzz pass). A one-hole star at (5e6, 5e6) ft in EPSG:2927, meshed on the
"smooth" profile, did not build on the first version of the tiling repair: GEOS refused
the union ``get_grid_centroid`` takes of the cells REPROJECTED to EPSG:4326 ("unable to
assign free hole to a shell"). 17 of its 156,582 cells were invalid once reprojected and
valid before, and every one of the 17 had an edge 9.3e-10 to 2.95e-9 ft long -- one to
three ulps at 5e6 ft, two circumcentres the dual kept apart because its merge tolerance
(1e-9 of the mesh edge) was finer than the coordinates. That grid took about two minutes
to build and is not rebuilt here; the first of the 17 cells is below, verbatim.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/voronoi_centroid_reprojection.py

A second.

WHAT IT MEASURES
  1. The cell as it was written: valid in EPSG:2927, and after ``to_crs(4326)``.
  2. The same cell with its two closest vertices made one, as ``dual_cells`` now does for
     circumcentres closer than 1000 ulp of the coordinates: after ``to_crs(4326)``.
The repair's own test of the merge is
``tests/test_voronoi_tiling.py::test_circumcentres_closer_than_the_coordinates_resolve_are_one_vertex``;
``get_grid_centroid`` also falls back to the model's own coordinates when the reprojected
union is refused, which no grid measured since reaches.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Polygon

CELL = [(5000397.5528519675, 5000563.556042625), (5000397.559567844, 5000563.544307067),
        (5000397.559567842, 5000563.544307066), (5000397.54014068, 5000563.5423113825),
        (5000397.5351805575, 5000563.552451226)]


def _report(label: str, ring) -> None:
    cell = Polygon(ring)
    edges = np.hypot(*np.diff(np.vstack([ring, ring[:1]]), axis=0).T)
    lonlat = gpd.GeoSeries([cell], crs=2927).to_crs(4326).iloc[0]
    print(f"{label}: {len(ring)} vertices, shortest edge {edges.min():.3g} ft, {cell.area:.3g} ft2; valid in "
          f"EPSG:2927 {cell.is_valid}; reprojected {shapely.is_valid_reason(lonlat)}")


def main() -> None:
    ring = np.asarray(CELL)
    _report("1. as written", ring)
    _report("2. merged", np.delete(ring, 2, axis=0))


if __name__ == "__main__":
    main()
