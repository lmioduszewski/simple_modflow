"""How general ledger 316's two ``tri2vor`` defects are: the sweep behind its rates.

Ledger 316; the mechanism, the MF6 consequence and the consumers are in
``voronoi_tiling_gaps.py`` beside this, whose helpers (the reference cell above all) this
imports so the two cannot disagree about what a bad cell is.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/voronoi_tiling_sweep.py

1,266 grids on 10 processes: several minutes.

THE SETUP. A 1000 ft square less a wedge at (400, 400), so the domain has ONE reflex
corner, of angle theta, rotated about (500, 500) so no answer rests on an axis-aligned
edge. A corner FAILS when its node's cell differs from its reference (mode 1). A MODE 2
EVENT is a Voronoi vertex exactly on the meshed boundary that ``point_in_polygon`` drops.
The hole AT the corner is the gap touching its node, reported against the corner's own
reference cell, against the mean cell and against ``max_area``.

WHAT IT MEASURES
  A. ``GridSpec.voronoi`` defaults (profile "balanced", min angle 32): theta x 6 rotations
     x max_area 500-20,000. The failure table, the hole against cell size, mode 2 events.
  B. The ``min_angle=`` quality setting, 20-34, at max_area 5000 (3 rotations).
  C. The build path: flopy's own Triangle + VoronoiGrid at angle 20 (its default) and 32,
     and myflopy's "baseline" profile (no cleanup, no smoothing), at max_area 1000 and 5000.
  D. Convex domains, where only mode 2 can occur: a square and a 1500 x 600 rectangle.
"""

from __future__ import annotations

import collections
import pathlib
import sys
import tempfile
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import shapely
from flopy.utils.geometry import point_in_polygon
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid
from scipy.spatial import Voronoi
from shapely import affinity
from shapely.geometry import Point, Polygon, box

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from voronoi_tiling_gaps import _bad, _cells, _reference, _wedge_notch  # noqa: E402

THETAS = (190, 200, 210, 225, 240, 250, 260, 270, 280, 290, 300, 315, 330, 345)
ROTATIONS = (0, 7, 17, 29, 41, 53)
AREAS = (500, 1000, 2000, 5000, 10000, 20000)
ANGLES = (20, 25, 28, 30, 32, 33, 34)
CONVEX = {"square": box(0, 0, 1000, 1000), "rectangle": box(0, 0, 1500, 600)}


def _build(path: str, domain: Polygon, max_area: float, angle: float, profile: str | None):
    root = pathlib.Path(tempfile.mkdtemp())
    if path == "flopy":
        tri = Triangle(model_ws=str(root), angle=angle, maximum_area=max_area)
        tri.add_polygon(domain)
        tri.add_region(domain.representative_point().coords[0], 0, maximum_area=max_area)
        tri.build()
        return VoronoiGrid(tri), tri
    options = {} if profile is None else {"profile": profile}
    if angle != 32:
        options["min_angle"] = angle
    import geopandas as gpd

    import myflopy as mf

    gpd.GeoDataFrame({"n": [0]}, geometry=[domain], crs=2927).to_file(root / "d.gpkg")
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "d.gpkg"),
                               max_area=float(max_area), **options)
    return spec.resolve(workspace=root / "g", return_triangle=True)


def _mode2_events(tri, meshed: Polygon) -> int:
    vertices = Voronoi(tri.verts).vertices
    on_edge = [v for v in vertices if meshed.exterior.distance(Point(v)) < 1e-6]
    if not on_edge:
        return 0
    kept = point_in_polygon(np.array([[v[0]] for v in on_edge]), np.array([[v[1]] for v in on_edge]),
                            list(tri._polygons[0])).ravel()
    return int(np.count_nonzero(~kept))


def _one(task):
    part, path, name, theta, rot, max_area, angle, profile = task
    warnings.simplefilter("ignore")
    domain = _wedge_notch(theta, rot) if name == "wedge" else affinity.rotate(CONVEX[name], rot, origin="centroid")
    try:
        vor, tri = _build(path, domain, max_area, angle, profile)
    except Exception as exc:  # noqa: BLE001 - a failed build is a row, not the end of the sweep
        return dict(part=part, error=f"{type(exc).__name__}: {exc}")
    reference, meshed = _reference(vor, tri)
    bad = _bad(vor, reference)
    gap = meshed.difference(shapely.union_all(_cells(vor)))
    row = dict(part=part, path=path, name=name, theta=theta, rot=rot, max_area=max_area, angle=angle,
               profile=profile, ncell=vor.ncpl, mean_cell=meshed.area / vor.ncpl, gap=gap.area,
               mode2=_mode2_events(tri, meshed))
    if name == "wedge":
        apex = affinity.rotate(Point(400, 400), rot, origin=(500, 500))
        node = [k for k in range(vor.ncpl) if Point(vor.points[k]).distance(apex) < 1e-6]
        row["corner_ok"] = bool(node)
        if node:
            k = node[0]
            row["failed"] = k in bad
            row["corner_cell"] = reference[k].area
            row["hole"] = sum(g.area for g in shapely.get_parts(gap) if g.distance(apex) < 1e-6)
            outside = shapely.union_all(_cells(vor)).difference(meshed)
            row["outside"] = sum(g.area for g in shapely.get_parts(outside) if g.distance(apex) < 1e-6)
    return row


def _tasks():
    for theta in THETAS:
        for rot in ROTATIONS:
            for area in AREAS:
                yield ("A", "myflopy", "wedge", theta, rot, area, 32, None)
        for rot in (0, 17, 41):
            for angle in ANGLES:
                yield ("B", "myflopy", "wedge", theta, rot, 5000, angle, None)
            for area in (1000, 5000):
                yield ("C", "flopy", "wedge", theta, rot, area, 20, None)
                yield ("C", "flopy", "wedge", theta, rot, area, 32, None)
                yield ("C", "myflopy", "wedge", theta, rot, area, 32, "baseline")
    for name in CONVEX:
        for rot in ROTATIONS:
            for area in AREAS:
                yield ("D", "flopy", name, 0, rot, area, 20, None)
                yield ("D", "flopy", name, 0, rot, area, 32, None)
                yield ("D", "myflopy", name, 0, rot, area, 32, None)


def _failure_table(rows, column, label):
    print(f"\n{label}: corners that failed / built  [median hole, % of the corner's own cell]")
    keys = sorted({column(r) for r in rows})
    print("  theta " + "".join(f"{str(k):>20s}" for k in keys))
    for theta in THETAS:
        cells = []
        for key in keys:
            group = [r for r in rows if r["theta"] == theta and column(r) == key and r.get("corner_ok")]
            failed = [r for r in group if r["failed"]]
            pct = f"{np.median([100 * r['hole'] / r['corner_cell'] for r in failed]):4.0f}%" if failed else "   -"
            cells.append(f"{len(failed):2d}/{len(group):2d} [{pct}]")
        print(f"  {theta:5d} " + "".join(f"{c:>20s}" for c in cells))


def main() -> None:
    tasks = list(_tasks())
    start = time.time()
    with ProcessPoolExecutor(max_workers=10) as pool:
        rows = list(pool.map(_one, tasks, chunksize=4))
    errors = [r for r in rows if "error" in r]
    rows = [r for r in rows if "error" not in r]
    print(f"{len(tasks)} grids in {time.time() - start:.0f} s; {len(errors)} failed to build "
          f"{dict(collections.Counter(e['error'][:60] for e in errors))}")
    missing = [r for r in rows if r["name"] == "wedge" and not r["corner_ok"]]
    print(f"wedges whose corner is not a node (not scored): {len(missing)}")

    part = [r for r in rows if r["part"] == "A"]
    _failure_table(part, lambda r: r["max_area"], "A. GridSpec.voronoi defaults, 6 rotations per cell")
    scored = [r for r in part if r.get("corner_ok")]
    print(f"  all: {sum(r['failed'] for r in scored)} of {len(scored)} corners failed; of those, "
          f"{sum(r['failed'] and r['hole'] < 1e-6 and r['outside'] < 1e-6 for r in scored)} have neither a "
          "hole nor cell area outside at the corner (a reference artifact would look like that)")
    print("\n  the hole against cell size, over failed corners of 250 deg and more:")
    print("  max_area    n  median ft2  max ft2   / max_area   / mean cell   / corner cell")
    for area in AREAS:
        f = [r for r in scored if r["max_area"] == area and r["failed"] and r["theta"] >= 250]
        if not f:
            continue
        hole = np.array([r["hole"] for r in f])
        print(f"  {area:8d} {len(f):4d} {np.median(hole):10.1f} {hole.max():8.1f} "
              f"{np.median(hole / area):11.4f} {np.median([r['hole'] / r['mean_cell'] for r in f]):12.3f} "
              f"{np.median([r['hole'] / r['corner_cell'] for r in f]):13.3f}")
    events = [r for r in part if r["mode2"]]
    print(f"  mode 2: {len(events)} of {len(part)} grids drop an on-boundary vertex "
          f"({sum(r['mode2'] for r in events)} vertices); "
          + ", ".join(f"theta {r['theta']} rot {r['rot']} max_area {r['max_area']}" for r in events[:8]))

    _failure_table([r for r in rows if r["part"] == "B"], lambda r: r["angle"],
                   "B. min_angle=, max_area 5000, 3 rotations per cell")
    _failure_table([r for r in rows if r["part"] == "C"],
                   lambda r: f"{r['path']} {r['profile'] or ''} {r['angle']}".replace("  ", " "),
                   "C. build path, max_area 1000 and 5000, 3 rotations each")

    print("\nD. convex domains, 6 rotations x 6 max_area each:")
    for name in CONVEX:
        for path, angle in (("flopy", 20), ("flopy", 32), ("myflopy", 32)):
            group = [r for r in rows if r["part"] == "D" and r["name"] == name and r["path"] == path and r["angle"] == angle]
            hit = [r for r in group if r["gap"] > 1e-6]
            worst = max(hit, key=lambda r: r["gap"] / (r["mean_cell"] * r["ncell"]), default=None)
            tail = (f"; worst {worst['gap']:,.0f} ft2 ({100 * worst['gap'] / (worst['mean_cell'] * worst['ncell']):.1f}% "
                    f"of the domain) at max_area {worst['max_area']}") if worst else ""
            print(f"  {name:9s} {path:7s} angle {angle}: {len(hit)} of {len(group)} grids have a gap{tail}")


if __name__ == "__main__":
    main()
