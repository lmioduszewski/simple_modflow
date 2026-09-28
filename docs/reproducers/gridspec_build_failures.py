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
  2. A max_area the domain cannot hold (no point 0.5*sqrt(max_area) from the boundary; for a
     W x H rectangle, max_area >= H**2): the sizing region was the domain itself, clipped to
     a ring 1e-4 of the diagonal inside it, and Triangle filled that band -- the L at 1e6
     gave 149,385 cells, smallest 0.00019 ft2, the same at 1e7. Now 8 cells and a warning.
     The L at 1e5 (207 cells, pin 4394bb5bf7eefb22) and just below the threshold at 1.58e5
     (2,851 cells, pin 84a3746283f05e8d -- a thin inset, fine-grained, unchanged on
     purpose) keep their pins.
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


def section_2(root):
    for max_area in (1e5, 1.58e5, 1e6):
        print(f"2. L at max_area {max_area:g}: {build(root / f'2-L-{max_area:g}', L_SHAPE, max_area=max_area)}")
    rectangle = box(0, 0, 2000, 1000)
    for max_area in (0.9e6, 1e6):
        out = build(root / f"2-rect-{max_area:g}", rectangle, max_area=max_area)
        print(f"2. 2000 x 1000 rectangle at max_area {max_area:g}: {out}")


def main():
    warnings.simplefilter("ignore")
    with tempfile.TemporaryDirectory() as root:
        section_1(Path(root))
        section_2(Path(root))


if __name__ == "__main__":
    main()
