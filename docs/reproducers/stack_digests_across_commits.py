"""Lazy-DISV digests of a spread of UN-excavated stacks, to diff across two commits.

Behind ledger 323's "no stored digest moves": writing excavations into the
``$layer_stack`` payload must leave a stack WITHOUT one fingerprinting exactly as
before. A digest is only comparable between two source trees, so run this under each
and diff (the files it reads are written into the directory you name, the same bytes
both times)::

    git archive <older-commit> src | tar -x -C /tmp/older
    PYTHONPATH=/tmp/older/src python docs/reproducers/stack_digests_across_commits.py /tmp/dg > before.txt
    PYTHONPATH=src            python docs/reproducers/stack_digests_across_commits.py /tmp/dg > after.txt
    diff before.txt after.txt

Measured 2026-09-25, 51a59ad (C5) against the C6 tree: 16 of 16 identical -- eight
stacks (flat; a raster top; a tuple split with a pinch; per-unit names; a raster in
metres; a contour bottom with a region vector; a ``Surface.where`` zone file; algebra
over a raster) under two build-kwarg sets.
"""

from __future__ import annotations

import pathlib
import sys
import warnings

import geopandas as gpd
import numpy as np
import rasterio
import shapely as shp
from rasterio.transform import from_origin

import myflopy as mf
from myflopy.specs import package_digest

warnings.simplefilter("ignore")


def _inputs(root: pathlib.Path):
    root.mkdir(parents=True, exist_ok=True)
    tif, zone, contours = root / "top.tif", root / "zone.gpkg", root / "c.gpkg"
    if not tif.exists():
        with rasterio.open(tif, "w", driver="GTiff", height=12, width=12, count=1,
                           dtype="float32", crs="EPSG:2927",
                           transform=from_origin(0, 600, 50, 50)) as handle:
            handle.write(np.full((12, 12), 100.0, dtype="float32"), 1)
    if not zone.exists():
        gpd.GeoDataFrame({"n": [1]}, geometry=[shp.box(0, 0, 300, 300)], crs=2927).to_file(zone)
    if not contours.exists():
        gpd.GeoDataFrame({"Elev": [90.0]}, geometry=[shp.LineString([(0, 0), (600, 600)])],
                         crs=2927).to_file(contours)
    return tif, zone, contours


def main(root: pathlib.Path) -> None:
    tif, zone, contours = _inputs(root)
    flat = mf.Surface.flat
    stacks = {
        "flat2": lambda: mf.LayerStack(top=flat(100.0)).add("sand", thickness=30.0)
        .add("clay", thickness=20.0),
        "raster": lambda: mf.LayerStack(top=mf.Surface.raster(tif)).add("sand", thickness=50),
        "split": lambda: mf.LayerStack(top=flat(100.0)).add("a", thickness=20.0, split=(0.3, 0.7))
        .add("b", thickness=5.0, pinch="inactive", min_thickness=2.0),
        "names": lambda: mf.LayerStack(top=flat(100.0)).add("g", thickness=5.0, split=2,
                                                             names=["u", "l"]),
        "units": lambda: mf.LayerStack(top=mf.Surface.raster(tif, units="meters"))
        .add("s", bottom=flat(40.0)),
        "contours": lambda: mf.LayerStack(top=flat(100.0)).add(
            "s", bottom=mf.Contours(contours, region_vector=zone, out=root / "c.tif")),
        "where": lambda: mf.LayerStack(top=flat(100.0)).add(
            "s", bottom=mf.Surface.where(zone, flat(80.0), flat(70.0))),
        "algebra": lambda: mf.LayerStack(top=mf.Surface.raster(tif) - 2.0)
        .add("s", thickness=10.0),
    }
    for name, make in stacks.items():
        for kwargs in ({}, {"default_pinch": "inactive", "min_sep": 0.5}):
            print(name, kwargs, package_digest(mf.disv(make(), **kwargs)))


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
