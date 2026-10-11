"""A strike/dip raster is the geological plane, at pixel centres (ledger 362).

Writes ``vor.get_raster_from_strike_dip`` on a 4x4 irregular Voronoi grid (10 ft
pixels, plane through (0, 0) at 100 ft), fits a plane to the raster at each pixel's
CENTRE, and prints the slope down the dip direction (strike + 90, the right-hand
rule) and along the strike, against ``tan(dip)``; then the normal of a level plane
and pixel (0, 0).

Before the fix (f4817ce and earlier), measured:

    strike=90 dip=45: fell toward azimuth 90 (EAST) by 0.5000, level north-south
    strike=0  dip=10: fell toward azimuth 0 (NORTH) by 0.1710 (tan 0.1763)
    strike=90 dip=80: 0.1710 -- the same as dip=10 (tan 5.6713)
    get_normal_from_strike_and_dip(0, 0) = [0, 1, 6e-17] (horizontal)
    pixel (0, 0) at strike=90 dip=45 = 100.0 (the plane at its upper-left corner;
    at its centre the old plane read 97.5)

i.e. it fell toward ``strike`` (along the strike line) by ``sin(2 * dip) / 2``.
After: down-dip slope equals tan(dip) to float32, along-strike slope 0, the level
normal is [0, 0, 1], and pixel (0, 0) holds the plane at its centre.

Run: ``PYTHONPATH=src python docs/reproducers/strike_dip_plane.py``
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import xy

from myflopy.modflow.mf6.canonical import irregular_voronoi_grid


def main() -> None:
    vor = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=100.0, refined=False)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "plane.tif"
        for strike, dip in [(90, 45), (0, 10), (90, 80), (37, 45)]:
            vor.get_raster_from_strike_dip(strike, dip, (0.0, 0.0, 100.0), pixel_size=10, output_filename=out)
            with rasterio.open(out) as src:
                values = src.read(1).astype(float)
                rows, cols = np.indices(values.shape)
                xs, ys = xy(src.transform, rows.ravel(), cols.ravel(), offset="center")
            xs, ys = np.asarray(xs), np.asarray(ys)
            coef = np.linalg.lstsq(np.c_[np.ones_like(xs), xs, ys], values.ravel(), rcond=None)[0]
            s = np.radians(strike)
            down_dip = -(coef[1] * np.cos(s) - coef[2] * np.sin(s))
            along_strike = coef[1] * np.sin(s) + coef[2] * np.cos(s)
            print(
                f"strike={strike:>3} dip={dip:>2}: down-dip slope {down_dip:.4f} "
                f"(tan {np.tan(np.radians(dip)):.4f}), along strike {along_strike:.1e}; "
                f"pixel (0, 0) centre ({xs[0]:.0f}, {ys[0]:.0f}) holds {values[0, 0]:.4f}"
            )
    print("normal of a level plane:", vor.get_normal_from_strike_and_dip(0, 0))


if __name__ == "__main__":
    main()
