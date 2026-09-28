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
  3. Breaklines on the L at (1.2e6, 6.5e5), profile "clean", min_angle 23, max_area 3000,
     breakline_max_area 60 (the fuzz input). Line B alone, a 20 ft line and the fuzz's five
     rows (A, B, B, C, C) raised "A linearring requires at least 4 coordinates"; they build,
     the duplicated rows dropped with a warning. A 40 ft line still vanishes, now named. Line
     A alone keeps pin c4c813f7b898e39b and a 100 ft line pin 824e8318f6c2acf2 -- A's pin on
     the repair from bb56668 was 4714a10d3ef12cf9, a grid that wrote a zero-length face.
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


def section_3(root):
    from shapely import affinity
    from shapely.geometry import LineString

    origin = (1.2e6, 6.5e5)
    domain = affinity.translate(L_SHAPE, *origin)
    a = LineString([(1199808.0861394587, 649924.4522525535), (1201191.9138605413, 651075.5477474465)])
    b = LineString([(1200400, 650581.5845702824), (1200509.9330967965, 650172.0865051051)])
    c = LineString([(1200464.0038896762, 650400), (1200428.6498325951, 650165.1040641002)])

    def line(length):
        return LineString([(origin[0] + 100, origin[1] + 200), (origin[0] + 100 + length, origin[1] + 200)])

    cases = {"A alone": [a], "B alone": [b], "A, B, B, C, C": [a, b, b, c, c], "20 ft line": [line(20)],
             "40 ft line": [line(40)], "100 ft line": [line(100)]}
    for i, (name, lines) in enumerate(cases.items()):
        folder = root / f"3-{i}"
        folder.mkdir(parents=True)
        gpd.GeoDataFrame({"n": list(range(len(lines)))}, geometry=lines, crs=2927).to_file(folder / "l.gpkg")
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            out = build(folder, domain, profile="clean", min_angle=23, max_area=3000.0, breakline_max_area=60.0,
                        breaklines=[mf.GeoPackageSourceSpec(path=folder / "l.gpkg")])
        said = sorted({str(w.message)[:60] for w in record if "duplicate" in str(w.message) or "resampled" in str(w.message)})
        print(f"3. {name}: {out}" + (f"  {said}" if said else ""))


def main():
    warnings.simplefilter("ignore")
    with tempfile.TemporaryDirectory() as root:
        section_1(Path(root))
        section_2(Path(root))
        section_3(Path(root))


if __name__ == "__main__":
    main()
