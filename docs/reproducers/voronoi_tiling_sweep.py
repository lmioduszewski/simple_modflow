"""How general ledger 316's two ``tri2vor`` defects are: the sweep behind its rates.

Ledger 316; the mechanism, the MF6 consequence and the consumers are in
``voronoi_tiling_gaps.py`` beside this, whose helpers (the reference cell above all) this
imports so the two cannot disagree about what a bad cell is.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/voronoi_tiling_sweep.py

1,510 grids on 10 processes: about 3-8 minutes (half an hour before the sizing fix of 2026-09-28,
which removed the ~200,000-cell sliver meshes).

THE SETUP. A 1000 ft square less a wedge at (400, 400), so the domain has ONE reflex
corner, of angle theta, rotated about (500, 500) so no answer rests on an axis-aligned
edge. Every rate is of FloPy's cells (``tessellation="flopy"``); each grid is also
repaired -- ``VoronoiGridPlus(tri)``, the default -- and the gap, overshoot and overlap it
leaves against the union of Triangle's triangles are summed at the end. A corner FAILS when its node's cell differs from its reference (mode 1). A MODE 2
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
  E. Spiky stars -- 8 to 40 points, alternating long and short radii, seeded -- at the
     origin, at (1.2e6, 6.5e5) and at (5e6, 5e6) ft, on the "baseline" and "clean"
     profiles at max_area 1000 and 10000. Across a thin spike Qhull's Delaunay (of all
     nodes) is not Triangle's, and state-plane coordinates are where a fuzz pass found the
     first repair losing faces.
  F. The grids FloPy gets right, moved: the square and rectangle of D (4 rotations x
     max_area 500-20,000) and wedges of 190-210 deg (3 rotations x max_area 1000 and 5000),
     at (1.2e6, 6.5e5) and (5e6, 5e6) ft. Qhull's roundoff grows with the coordinates, so
     this is where a repair matching FloPy's vertices to 1000 ulp would re-pin a grid with
     nothing wrong, if it were going to.

Every grid, every part, also counts LOST FACES: mesh edges whose dual segment is longer
than 1e-7 of the edge (an interior edge between two distinct circumcentres, a boundary
edge whose circumcentre is not its midpoint) AND than what a DISV file FloPy writes can
hold -- one unit in its last digit, on at least one axis (``written_resolution``: 0.01 ft
at x = 1.2e6, 1e-8 ft below 1e5) -- and that the two cells do not share as a vertex pair,
which is how MF6 DISV finds a face -- for FloPy's cells and for the repaired grid. A face
shorter than that is written as one point: kept, it reaches MF6 with zero length and MF6
returns NaN for every head, which a real 20,323-cell model at state plane did after the
repair split a corner 7.4e-5 ft wide (ledger 316, 2026-09-28). So each grid also counts
the faces WRITTEN with zero length -- on FloPy's cells, on the repaired grid, and on the
grid the repair built from bb56668 until corners closer than the file holds were joined
(``_onto_flopy``'s ``split``) -- and how many dual vertices the repair joined; and the faces
written with BOTH cell centres on their line, the build check's second refusal (CONNLEN 0,
a division by zero in MF6). And for each grid FloPy got right (no lost face; no gap, overshoot or overlap; no
ring holding one location twice, which is a zero-width MF6 connection), it says whether
the repair changed it anyway: that would re-pin a grid with nothing wrong.

A gap, overshoot or overlap is measured against the union of Triangle's triangles with
GEOS, and on a grid of 200,000 cells a tenth of a foot across at state-plane coordinates
the union itself is off (measured: 0.00065 ft2 of "overshoot" to which no cell contributes
1e-9). So a misfit of the repaired grid above 1e-6 ft2 is re-measured cell by cell, where
it is: each cell's own area outside the triangles near it, and the gap no cell covers.
"""

from __future__ import annotations

import collections
import math
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
STAR_SEEDS = range(12)
OFFSETS = ((0.0, 0.0), (1.2e6, 6.5e5), (5e6, 5e6))


def _build(path: str, domain: Polygon, max_area: float, angle: float, profile: str | None,
           root: pathlib.Path | None = None):
    root = pathlib.Path(tempfile.mkdtemp()) if root is None else root
    if path == "flopy":
        tri = Triangle(model_ws=str(root), angle=angle, maximum_area=max_area)
        tri.add_polygon(domain)
        tri.add_region(domain.representative_point().coords[0], 0, maximum_area=max_area)
        tri.build()
        return VoronoiGrid(tri), tri
    options = {"tessellation": "flopy"} if profile is None else {"profile": profile, "tessellation": "flopy"}
    if angle != 32:
        options["min_angle"] = angle
    import geopandas as gpd

    import myflopy as mf

    gpd.GeoDataFrame({"n": [0]}, geometry=[domain], crs=2927).to_file(root / "d.gpkg")
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "d.gpkg"),
                               max_area=float(max_area), **options)
    return spec.resolve(workspace=root / "g", return_triangle=True)


def _lost_faces(vor, tri, short: bool = False) -> int:
    """Mesh faces the cells do not share as a vertex pair (see the module docstring).

    ``short=True`` counts instead the faces longer than 1e-7 of their edge but shorter than
    the written file holds: the ones the repair may make one vertex. Pass the repaired grid,
    which has a cell for every node.
    """
    from myflopy.grid_identity import written_resolution

    P = np.asarray(tri.verts, dtype=float)
    T = np.asarray(tri.iverts, dtype=np.int64)
    cell_of = {(float(x), float(y)): k for k, (x, y) in enumerate(np.asarray(vor.points))}
    node_cell = [cell_of.get((float(x), float(y))) for x, y in P]
    P = P - P.min(axis=0)  # the circumcentre below squares coordinates: 1e12 at state plane
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    ba, ca = b - a, c - a
    d = 2.0 * (ba[:, 0] * ca[:, 1] - ba[:, 1] * ca[:, 0])
    b2, c2 = (ba**2).sum(1), (ca**2).sum(1)
    cc = a + np.column_stack([(ca[:, 1] * b2 - ba[:, 1] * c2) / d, (ba[:, 0] * c2 - ca[:, 0] * b2) / d])
    edges = collections.defaultdict(list)
    for t, (i, j, k) in enumerate(T):
        for u, v in ((i, j), (j, k), (k, i)):
            edges[(min(u, v), max(u, v))].append(t)
    unit = written_resolution(np.asarray(tri.verts, dtype=float))
    want = set()
    for (u, v), ts in edges.items():
        if node_cell[u] is None or node_cell[v] is None:
            continue
        end = cc[ts[1]] if len(ts) == 2 else (P[u] + P[v]) / 2
        dual, floor = np.hypot(*(cc[ts[0]] - end)), 1e-7 * np.hypot(*(P[u] - P[v]))
        unwritable = bool((np.abs(cc[ts[0]] - end) < unit).all())
        if (floor < dual and unwritable) if short else (dual > floor and not unwritable):
            want.add((min(node_cell[u], node_cell[v]), max(node_cell[u], node_cell[v])))
    if short:
        return len(want)
    owner, have = {}, set()
    for k, ring in enumerate(vor.iverts):
        for i in range(len(ring)):
            edge = frozenset((int(ring[i]), int(ring[(i + 1) % len(ring)])))
            if edge in owner and owner[edge] != k:
                have.add((min(owner[edge], k), max(owner[edge], k)))
            owner[edge] = k
    return len(want - have)


def _unwritable(verts, iverts, points) -> collections.Counter:
    """Faces a DISV file FloPy writes (its default number format) holds so MF6 cannot compute
    them, by why: ``"zero length"`` (NaN heads) or ``"centres on it"`` (a division by zero).
    The centres are the generators, as ``VoronoiGrid.get_disv_gridprops`` writes them."""
    from types import SimpleNamespace

    from flopy.utils.cvfdutil import get_disv_gridprops

    from myflopy.grid_identity import unwritable_faces

    props = get_disv_gridprops(np.asarray(verts, dtype=float), [list(r) for r in iverts], xcyc=np.asarray(points))
    vertices = np.rec.fromrecords([tuple(v) for v in props["vertices"]], names="iv,xv,yv")
    disv = SimpleNamespace(vertices=SimpleNamespace(get_data=lambda: vertices),
                           cell2d=SimpleNamespace(get_data=lambda: props["cell2d"]))
    return collections.Counter(face[-1] for face in unwritable_faces(disv))


def _zero_written(verts, iverts, points) -> int:
    """Faces a DISV file FloPy writes (its default number format) holds with zero length."""
    return _unwritable(verts, iverts, points)["zero length"]


def _split_grid(vor, tri):
    """The repaired grid as bb56668 built it: every dual vertex kept, however close."""
    from myflopy.modflow.mf6.grid.tessellation import _dual, _onto_flopy

    d_verts, d_iverts, d_points, on_boundary = _dual(tri)
    if not np.array_equal(np.asarray(vor.points, dtype=float), d_points):
        return d_verts, d_iverts, d_points  # adopted whole; bb56668 joined nothing there either
    verts, iverts, *_, split = _onto_flopy(d_verts, d_iverts, vor.verts, vor.iverts, on_boundary)
    return (split[0], split[1], d_points) if split is not None else (verts, iverts, d_points)


def _star(seed: int) -> Polygon | None:
    """8 to 40 points around (500, 500), alternating long and short radii, as a fuzz pass drew them."""
    rng = np.random.default_rng(seed)
    angles = np.unique(np.round(np.sort(rng.uniform(0, 2 * math.pi, int(rng.integers(8, 40)))), 4))
    radius = np.where(np.arange(len(angles)) % 2 == 0, rng.uniform(380, 450, len(angles)),
                      rng.uniform(25, 135, len(angles)))
    star = Polygon([(500 + r * math.cos(t), 500 + r * math.sin(t)) for t, r in zip(angles, radius, strict=True)])
    return star if star.is_valid else None


def _mode2_events(tri, meshed: Polygon) -> int:
    vertices = Voronoi(tri.verts).vertices
    on_edge = [v for v in vertices if meshed.exterior.distance(Point(v)) < 1e-6]
    if not on_edge:
        return 0
    kept = point_in_polygon(np.array([[v[0]] for v in on_edge]), np.array([[v[1]] for v in on_edge]),
                            list(tri._polygons[0])).ravel()
    return int(np.count_nonzero(~kept))


def _one(task):
    """One grid, measured in a temporary directory that is removed after it."""
    with tempfile.TemporaryDirectory() as tmp:
        return _measure(task, pathlib.Path(tmp))


def _measure(task, root: pathlib.Path):
    part, path, name, theta, rot, max_area, angle, profile = task
    warnings.simplefilter("ignore")
    import logging

    logging.disable(logging.WARNING)  # the repair's adopted-grid WARNING, once per star
    rotation, offset = rot if part == "F" else (rot, None)  # F carries (rotation, offset)
    if name == "wedge":
        domain = _wedge_notch(theta, rotation)
    elif name == "star":
        domain = affinity.translate(_star(theta), *rot)  # theta is the seed, rot the offset
    else:
        domain = affinity.rotate(CONVEX[name], rotation, origin="centroid")
    if offset is not None:
        domain = affinity.translate(domain, *offset)
    try:
        vor, tri = _build(path, domain, max_area, angle, profile, root)
    except Exception as exc:  # noqa: BLE001 - a failed build is a row, not the end of the sweep
        return dict(part=part, error=f"{type(exc).__name__}: {exc}")
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    try:
        repaired = VoronoiGridPlus(tri)
    except Exception as exc:  # noqa: BLE001 - a refused mesh is a row, not the end of the sweep
        return dict(part=part, error=f"repair: {type(exc).__name__}: {exc}", task=task, ncell=vor.ncpl,
                    zero_before=_zero_written(vor.verts, vor.iverts, vor.points),
                    zero_split=_zero_written(*_split_grid(vor, tri)))
    pieces = [Polygon(np.asarray(tri.verts)[t]) for t in np.asarray(tri.iverts)]
    triangles = shapely.union_all(pieces)

    def misfit(grid, local=False):
        cells = _cells(grid)
        union = shapely.union_all(cells)
        gap, over = triangles.difference(union), union.difference(triangles)
        total = gap.area + over.area + abs(sum(c.area for c in cells) - union.area)
        if not local or total <= 1e-6:
            return total
        cell_tree, piece_tree = shapely.STRtree(cells), shapely.STRtree(pieces)
        out = 0.0
        for k in cell_tree.query(over):
            near = shapely.union_all([pieces[j] for j in piece_tree.query(cells[k])])
            out += cells[k].difference(near).area
        uncovered = gap.difference(shapely.union_all([cells[k] for k in cell_tree.query(gap)])).area
        overlap = 0.0
        region = gap.union(over)
        suspects = sorted(set(cell_tree.query(region.buffer(1.0))))
        for k in suspects:
            for j in cell_tree.query(cells[k]):
                if j > k:
                    overlap += cells[k].intersection(cells[j]).area
        return out + uncovered + overlap

    flopy_rings = [np.asarray(vor.verts, dtype=float)[list(ring)] for ring in vor.iverts]
    ulps = 1000 * float(np.spacing(np.abs(np.asarray(tri.verts, dtype=float)).max()))
    twice = sum(bool(np.any(np.hypot(*(ring - np.roll(ring, -1, axis=0)).T) < ulps)) for ring in flopy_rings)

    row = dict(part=part, path=path, name=name, theta=theta, rot=rot, max_area=max_area, angle=angle,
               profile=profile, ncell=vor.ncpl, after=misfit(repaired, local=True), before=misfit(vor), twice=twice,
               lost=_lost_faces(vor, tri), lost_after=_lost_faces(repaired, tri), short=_lost_faces(repaired, tri, True),
               changed=repaired.tessellation_report.changed, unsplit=repaired.tessellation_report.unsplit,
               zero_before=_zero_written(vor.verts, vor.iverts, vor.points),
               zero_split=_zero_written(*_split_grid(vor, tri)))
    faces = _unwritable(repaired.verts, repaired.iverts, repaired.points)
    row.update(zero_after=faces["zero length"], centres_after=faces["centres on it"],
               centres_before=_unwritable(vor.verts, vor.iverts, vor.points)["centres on it"])
    if name == "star" or part == "F":
        return row
    reference, meshed = _reference(vor, tri)
    bad = _bad(vor, reference)
    gap = meshed.difference(shapely.union_all(_cells(vor)))
    row.update(mean_cell=meshed.area / vor.ncpl, gap=gap.area, mode2=_mode2_events(tri, meshed))
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
    for seed in STAR_SEEDS:
        if _star(seed) is None:
            continue
        for offset in OFFSETS:
            for profile in ("baseline", "clean"):
                for area in (1000, 10000):
                    yield ("E", "myflopy", "star", seed, offset, area, 21, profile)
    for offset in OFFSETS[1:]:
        for name in CONVEX:
            for rot in ROTATIONS[:4]:
                for area in (500, 2000, 5000, 20000):
                    yield ("F", "myflopy", name, 0, (rot, offset), area, 32, None)
        for theta in (190, 200, 210):
            for rot in (0, 17, 41):
                for area in (1000, 5000):
                    yield ("F", "myflopy", "wedge", theta, (rot, offset), area, 32, None)


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
    missing = [r for r in rows if r["name"] == "wedge" and r["part"] != "F" and not r.get("corner_ok")]
    print(f"wedges whose corner is not a node (not scored): {len(missing)}")

    part = [r for r in rows if r["part"] == "A"]
    for r in part:
        r.setdefault("corner_ok", False)
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
    print("\nE. spiky stars, 'baseline' and 'clean' profiles, max_area 1000 and 10000:")
    for offset in OFFSETS:
        group = [r for r in rows if r["part"] == "E" and tuple(r["rot"]) == offset]
        print(f"  at {offset}: {sum(r['lost'] > 0 for r in group)} of {len(group)} grids lose faces on FloPy's "
              f"cells ({sum(r['lost'] for r in group)} faces), {sum(r['lost_after'] > 0 for r in group)} on the "
              f"repaired grid; {sum(r['before'] > 1e-6 for r in group)} leave a gap, overshoot or overlap on "
              f"FloPy's cells, {sum(r['after'] > 1e-6 for r in group)} repaired; {sum(r['short'] for r in group)} "
              f"faces shorter than the written file holds in {sum(r['short'] > 0 for r in group)} grids; "
              f"{sum(r['unsplit'] for r in group)} dual vertices joined in {sum(r['unsplit'] > 0 for r in group)}")
    print("\nF. grids FloPy gets right, at state-plane offsets:")
    for offset in OFFSETS[1:]:
        group = [r for r in rows if r["part"] == "F" and tuple(r["rot"][1]) == offset]
        good = [r for r in group if r["lost"] == 0 and r["before"] <= 1e-6 and not r["twice"]]
        print(f"  at {offset}: FloPy's cells right in {len(good)} of {len(group)} grids; the repair changed "
              f"{sum(r['changed'] for r in good)} of those; repaired, {sum(r['lost_after'] > 0 for r in group)} lose "
              f"a face and {sum(r['after'] > 1e-6 for r in group)} leave a gap, overshoot or overlap")
    worst = max(rows, key=lambda r: r["after"])
    right = [r for r in rows if r["lost"] == 0 and r["before"] <= 1e-6 and not r["twice"]]
    print(f"\nREPAIRED (the default): {sum(r['after'] > 1e-6 for r in rows)} of {len(rows)} grids leave a gap, "
          f"overshoot or overlap above 1e-6 ft2 (the largest {worst['after']:.2g} ft2), and "
          f"{sum(r['lost_after'] > 0 for r in rows)} lose a face, where FloPy's cells lost "
          f"{sum(r['lost'] for r in rows)} faces in {sum(r['lost'] > 0 for r in rows)} grids. Of the "
          f"{len(right)} grids FloPy got right, the repair changed {sum(r['changed'] for r in right)}; "
          f"{sum(r['twice'] > 0 and r['lost'] == 0 and r['before'] <= 1e-6 for r in rows)} more tiled but "
          "held a location twice in a ring, and the repair took the zero-width connection out.")
    for label, key in (("FloPy's cells", "zero_before"), ("the repair as bb56668 built it", "zero_split"),
                       ("the repaired grid", "zero_after")):
        print(f"WRITTEN WITH ZERO LENGTH (NaN in MF6), {label}: {sum(r[key] for r in rows)} faces in "
              f"{sum(r[key] > 0 for r in rows)} grids")
    def where(r):
        return f"{r['part']} {r['name']} {r['theta']} {r['rot']} max_area {r['max_area']} {r['profile']} ({r['ncell']} cells)"

    for label, key in (("FloPy's cells", "centres_before"), ("the repaired grid", "centres_after")):
        print(f"WRITTEN WITH BOTH CELL CENTRES ON THE FACE (CONNLEN 0, MF6 divides by zero), {label}: "
              f"{sum(r[key] for r in rows)} faces in {sum(r[key] > 0 for r in rows)} grids"
              + "".join(f"\n  {where(r)}: {r[key]}" for r in rows if r[key]))

    for r in errors:
        if "task" in r:
            print(f"  REFUSED {r['task']} ({r['ncell']} cells): FloPy's cells wrote {r['zero_before']} "
                  f"zero-length faces, bb56668's {r['zero_split']}")
    for r in rows:
        if r["zero_after"] or r["after"] > 1e-6:
            print(f"  {where(r)}: joined {r['unsplit']}, misfit {r['after']:.3g} ft2 (FloPy's {r['before']:.3g}); "
                  f"zero-length faces FloPy {r['zero_before']}, bb56668 {r['zero_split']}, now {r['zero_after']}")
    bad_right = [r for r in right if r["zero_before"]]
    print(f"  of the {len(right)} grids FloPy got right, {len(bad_right)} still write a zero-length face on "
          f"FloPy's own cells: " + "; ".join(where(r) for r in bad_right[:10]))
    joined = [r for r in rows if r["unsplit"]]
    print(f"JOINED: {sum(r['unsplit'] for r in joined)} dual vertices in {len(joined)} grids "
          f"({dict(collections.Counter((r['part'], str(r['rot'])) for r in joined))}); of the {len(right)} grids "
          f"FloPy got right, {sum(r['unsplit'] > 0 for r in right)} had any joined.")


if __name__ == "__main__":
    main()
