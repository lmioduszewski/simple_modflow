"""Five ``GridSpec.voronoi`` inputs a fuzz pass found failing, and what each fix changes.

Ledger 316 ("The five, diagnosed and resolved 2026-09-28").

    /home/lukem/python/envs/gw/bin/python docs/reproducers/gridspec_build_failures.py

Every build here runs in a ``tempfile.TemporaryDirectory``. Each section builds the failing
input and a working neighbour, and prints the working neighbour's ``grid_pin`` -- the pins
printed are the ones the grids had before the fixes, which is the claim: nothing that built
moved. Run it once against a checkout from before 2026-09-28 (``PYTHONPATH=<old>/src``) to
see the failures.

  1. Two rings sharing an exact point: the L with box(50, 50, 300, 300) at max_area 10000
     (its corner on the background sizing ring's corner) and two adjacent refinement zones
     raised FloPy's bare "duplicate points" Exception on every profile. The working
     neighbour, box(60, 60, 300, 300), keeps pin b8a48a446b0910ad on the default profile.
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Polygon, box

import myflopy as mf
from myflopy.grid_identity import grid_pin

L_SHAPE = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])


def build(root, boundary, *, refinement=(), **options):
    """``GridSpec.voronoi`` on ``boundary``: ``(ncpl, pin)`` or the error it raised."""
    root.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame({"n": [0]}, geometry=[boundary], crs=2927).to_file(root / "b.gpkg")
    kwargs = {}
    if refinement:
        gpd.GeoDataFrame({"area": [a for _, a in refinement]}, geometry=[g for g, _ in refinement],
                         crs=2927).to_file(root / "r.gpkg")
        kwargs["refinement"] = mf.GeoPackageSourceSpec(path=root / "r.gpkg", fields={"area": "area"})
    try:
        vor = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "b.gpkg"), **kwargs,
                                  **options).resolve(workspace=root / "grid")
    except Exception as exc:  # noqa: BLE001 - a failing input is what this prints
        return f"{type(exc).__name__}: {str(exc)[:90]}"
    return f"{vor.ncpl} cells, pin {grid_pin(vor)}"


def section_1(root):
    cases = {
        "box(50,50,300,300), corner on the sizing ring": [(box(50, 50, 300, 300), 2500.0)],
        "two adjacent zones": [(box(100, 100, 200, 200), 2500.0), (box(200, 100, 300, 200), 2500.0)],
        "box(60,60,300,300), working neighbour": [(box(60, 60, 300, 300), 2500.0)],
    }
    for profile in ("balanced", "clean", "baseline"):
        for name, regions in cases.items():
            out = build(root / f"1-{profile}-{len(name)}", L_SHAPE, refinement=regions, max_area=10000.0,
                        min_angle=28, profile=profile)
            print(f"1. {profile:8s} {name}: {out}")


def main():
    warnings.simplefilter("ignore")
    with tempfile.TemporaryDirectory() as root:
        section_1(Path(root))


if __name__ == "__main__":
    main()
