"""`mf.Contours` interpolated its contours at their RAW coordinates (ledger 602).

GRASS works in one frame, its LOCATION's CRS, and myflopy imported with ``v.in.ogr -o``
(override the projection check), so a contour file in another CRS was read as if its
numbers were already in the location's: EPSG:4326 degrees inside an EPSG:2927 location
land near the origin, reach none of the region, and the interpolation comes back blank --
with no warning, and the blank cache then reads ``fresh`` forever.

This runs both routes on the same synthetic data, through the same GRASS pipeline
(`ContourSurfaceInterpolator._run_grass`): RAW hands GRASS the files as given -- what
``run()`` did before -- and STAGED is ``run()`` now, which reprojects the contours and a
region polygon into the location's CRS first and refuses what it cannot reconcile.
Measured here (GRASS 8.4.2, 10 ft cells over a 1000 x 600 ft box at the origin):

- contours in 4326: RAW 0 of 6,000 pixels valid; STAGED 6,000, equal to the 2927 file
  (max |diff| 0.0, same mask);
- a region polygon in 4326: RAW a one-pixel raster with nothing valid (the region is a
  box of degrees); STAGED equal to the 2927 region, 6,000 valid;
- ``epsg="32610"`` with 32610 inputs on an EXISTING 2927 location: RAW writes a raster
  LABELLED EPSG:2927 holding metre coordinates (553 valid pixels; GRASS reads ``epsg=`` only when it
  creates a location); STAGED refuses, naming ``location=``;
- a file already in 2927: the staged copy interpolates exactly as the file itself
  (max |diff| 0.0), so no surface whose files are in the frame changes.

Needs GRASS GIS (exits with a message without it). Writes only in a temporary directory.

Run: `python docs/reproducers/contours_in_another_crs.py`
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
import shapely

from myflopy.modflow.utils import contour_interp
from myflopy.modflow.utils.contour_interp import ContourSurfaceInterpolator

if shutil.which("grass") is None and not os.environ.get("GRASS_BIN"):
    sys.exit("GRASS GIS is not installed (no `grass` on PATH, no GRASS_BIN): nothing to measure.")


def _raw(interp: ContourSurfaceInterpolator) -> Path:
    """The pipeline as it ran before: the files handed to GRASS as given."""

    r, g, v, gsetup = contour_interp._grass_modules(interp.grass_bin)
    return interp._run_grass(r, g, v, gsetup, contours=interp.contours,
                             region_vector=interp.region_vector)


def _read(path: Path):
    with rasterio.open(path) as handle:
        return handle.read(1, masked=True), handle.crs.to_string()


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="contours-crs-") as tmp:
        work = Path(tmp)
        lines = [shapely.LineString([(0, y), (1000, y + 20)]) for y in range(50, 600, 50)]
        contours = gpd.GeoDataFrame({"elev": [100.0 - y / 10 for y in range(50, 600, 50)]},
                                    geometry=lines, crs=2927)
        region = gpd.GeoDataFrame(geometry=[shapely.box(0, 0, 1000, 600)], crs=2927)
        for name, frame in (("c2927", contours), ("c4326", contours.to_crs(4326)),
                            ("c32610", contours.to_crs(32610)), ("r2927", region),
                            ("r4326", region.to_crs(4326)), ("r32610", region.to_crs(32610))):
            frame.to_file(work / f"{name}.gpkg")

        def interp(c, rg, out, **kwargs):
            kwargs.setdefault("grassdata", work / "grassdata")
            return ContourSurfaceInterpolator(work / f"{c}.gpkg", region_vector=work / f"{rg}.gpkg",
                                              out=work / f"{out}.tif", z_field="elev",
                                              resolution=10, **kwargs)

        truth, _ = _read(interp("c2927", "r2927", "truth").run())
        print(f"2927 contours, 2927 region: {truth.count()} of {truth.size} pixels valid")

        for c, rg, label in (("c4326", "r2927", "contours in 4326"),
                             ("c2927", "r4326", "region polygon in 4326")):
            raw, _ = _read(_raw(interp(c, rg, f"raw_{c}_{rg}")))
            staged, _ = _read(interp(c, rg, f"staged_{c}_{rg}").run())
            same = bool((staged.mask == truth.mask).all())
            print(f"{label}: RAW {raw.count()} of {raw.size} valid; STAGED {staged.count()} "
                  f"valid, max |diff| vs 2927 {float(np.abs(staged - truth).max()):.3g}, "
                  f"same mask {same}")

        raw, label = _read(_raw(interp("c32610", "r32610", "raw_32610", epsg="32610")))
        print(f"epsg='32610' on the existing 2927 location: RAW writes a raster labelled "
              f"{label} with {raw.count()} valid pixels at metre coordinates")
        try:
            interp("c32610", "r32610", "staged_32610", epsg="32610").run()
        except ValueError as refusal:
            print(f"  STAGED refuses: {str(refusal)[:110]}...")

        raw, _ = _read(_raw(interp("c2927", "r2927", "raw_2927")))
        print(f"a file already in 2927: staged vs direct import, max |diff| "
              f"{float(np.abs(raw - truth).max()):.3g}, same mask {bool((raw.mask == truth.mask).all())}")


if __name__ == "__main__":
    main()
