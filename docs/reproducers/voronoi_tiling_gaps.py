"""Voronoi cells that did not tile their boundary: FloPy's ``tri2vor``, and the repair.

Ledgers 316 and 317 (section 7). Ledger 314 found a 77.4 ft2 blank at the L's inner corner
of a masked surface and put it down to the grid's own footprint. This is where that hole
came from, how general it was, what read it -- and, on each, what the repaired default
grid does instead.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/voronoi_tiling_gaps.py

About a minute; it runs MF6 four times.

THE SETUP. The L of surface_domain_mask.py -- 1000 x 1000 ft with the north-east 600 x 600
removed, EPSG:2927 -- through ``GridSpec.voronoi(max_area=...)`` on its defaults (profile
"balanced", min angle 32). "FloPy's cells" are that grid built with
``tessellation="flopy"``, which keeps ``tri2vor``'s output as it is; "repaired" is the
default, ``VoronoiGridPlus(tri)`` on the same triangulation. Every cell is compared with a
REFERENCE: the node's restricted
Voronoi cell, i.e. the bisector half-planes of its Triangle neighbours intersected, clipped
to the polygon Triangle meshed, keeping the piece that touches the node. Not
``shapely.voronoi_polygons``: on cocircular nodes (an arc of boundary) GEOS returned a
region that is not a Voronoi cell, and good cells read as bad. A cell is BAD when it
differs from its reference by more than 1e-7 of the reference's area. The reference has a
limit of its own: at reflex corners of about 330 deg and sharper it is not an exact
tiling (a node across the notch owns a sliver beside the corner), so there it can flag a
correct cell. The hole and outside areas never use it, and the sweep counts how many of
its failures come with neither.

WHAT IT MEASURES
  1. MODE 1, a cell whose node IS a reflex corner. ``get_sorted_vertices`` orders a
     cell's vertices by angle about their MEAN. The corner cell is not convex: an angular
     sort is right only from inside its kernel, the wedge of 360 - theta at the corner.
     The two boundary midpoints pull the mean into the notch and the circumcentres pull
     it out, so it lands near the corner node (0.06 ft from it on the L) and the side it
     lands on decides the cell. Outside the wedge the corner vertex is sorted to the wrong
     place, so the polygon cuts across the notch (area OUTSIDE the boundary) and dips to the corner
     (a HOLE inside it), and it no longer shares a face with the cell across the hole.
  2. The same with no myflopy code at all: flopy's own Triangle + VoronoiGrid.
  3. How general, on the default path: reflex angle x rotation, at max_area 5000.
  4. MODE 2, a Voronoi vertex exactly ON the boundary -- the circumcentre of a right
     triangle whose hypotenuse is a boundary segment. ``point_in_polygon`` is a half-open
     ray test, so such a vertex on a right-hand or top edge counts as outside and is
     dropped, and the INTERIOR cell that needed it loses a triangle. The L at max_area
     10000, same defaults.
  5. What MF6 makes of it: connections from MF6's own .grb, and the steady flow around the
     L, on FloPy's cells and on the repaired grid (the same nodes; only broken cells
     re-outlined, on FloPy's own vertex ids).
  6. What myflopy's geometry makes of it: a feature inside the hole, the grid-edge cells,
     cell_coverage between two grids of the one boundary, and total area.
  7. A boundary that is not ONE simple polygon: with an interior ring it does not build,
     and in two parts FloPy's cells cover one part only -- tri2vor's domain is
     ``_polygons[0]`` alone. The repaired grid takes the dual's cells there.
  8. Not tri2vor: the default profile ("balanced", simplify_tolerance=2.0) simplifies the
     DOMAIN before meshing, so on a curved boundary the cells tile the simplified polygon,
     not the one supplied. A corridor 300 ft wide along a sine, at max_area 5000.
  9. What the repair costs: its own CPU time inside ``GridSpec.resolve`` of the L, against
     the rest of the resolve. It runs on every grid the smoothing loop builds (six per
     resolve on the "balanced" profile), not only the last.
"""

from __future__ import annotations

import math
import pathlib
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
import shapely
from flopy.mf6.utils import MfGrdFile
from flopy.utils.geometry import point_in_polygon
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid
from scipy.spatial import Voronoi
from shapely import affinity
from shapely.geometry import LineString, Point, Polygon, box

import myflopy as mf
from myflopy.modflow.mf6.grid.geometry import (
    cell_coverage,
    cell_overlay,
    cell_shares,
    find_adjacent_polygons,
    miss_reason,
)
from myflopy.modflow.mf6.grid.selection import get_grid_edge_cells
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

L_SHAPE = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])


def _grid(root: pathlib.Path, domain: Polygon, max_area: float, name: str, tessellation: str = "flopy"):
    """``GridSpec.voronoi`` on its defaults, returning the grid and its triangulation.

    ``tessellation="flopy"`` by default: this measures what FloPy's cells did.
    """
    path = root / f"{name}.gpkg"
    gpd.GeoDataFrame({"name": ["domain"]}, geometry=[domain], crs=2927).to_file(path)
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=path), max_area=float(max_area),
                               tessellation=tessellation)
    return spec.resolve(workspace=root / name, return_triangle=True)


def _wedge_notch(theta: float, rot: float) -> Polygon:
    """A 1000 ft square less a wedge at (400, 400): one reflex corner of angle theta."""
    alpha, bisector, far = math.radians(360.0 - theta), math.radians(45.0), 1.0e4
    ray = lambda a: (400 + far * math.cos(a), 400 + far * math.sin(a))  # noqa: E731
    wedge = Polygon([(400, 400), ray(bisector - alpha / 2), ray(bisector), ray(bisector + alpha / 2)])
    return affinity.rotate(box(0, 0, 1000, 1000).difference(wedge), rot, origin=(500, 500))


def _reflex_corners(domain: Polygon) -> dict:
    """{(x, y): interior angle} for every reflex vertex of the exterior ring."""
    coords = np.asarray(domain.exterior.coords)[:-1]
    out = {}
    for i in range(len(coords)):
        v1, v2 = coords[i - 1] - coords[i], coords[(i + 1) % len(coords)] - coords[i]
        turn = math.degrees(math.atan2(v1[0] * v2[1] - v1[1] * v2[0], v1 @ v2)) % 360
        angle = 360 - turn if domain.exterior.is_ccw else turn
        if angle > 180 + 1e-6:
            out[(round(coords[i][0], 6), round(coords[i][1], 6))] = angle
    return out


def _reference(vor, tri):
    """The reference cell of every grid cell, the triangulation's neighbour sets, and the meshed polygon."""
    meshed = Polygon(tri._polygons[0])
    tverts = np.asarray(tri.verts, dtype=float)
    index = {(float(x), float(y)): i for i, (x, y) in enumerate(tverts)}
    neighbours = [set() for _ in tverts]
    for a, b, c in tri.iverts:
        neighbours[a] |= {b, c}
        neighbours[b] |= {a, c}
        neighbours[c] |= {a, b}
    x0, y0, x1, y1 = meshed.bounds
    reach = 4 * max(x1 - x0, y1 - y0)
    cells = []
    for p in np.asarray(vor.points, dtype=float):
        cell = box(x0 - reach, y0 - reach, x1 + reach, y1 + reach)
        for j in neighbours[index[(float(p[0]), float(p[1]))]]:
            q = tverts[j]
            mid, d = (p + q) / 2, (q - p) / np.hypot(*(q - p))
            n = np.array([-d[1], d[0]]) * 3 * reach
            cell = cell.intersection(Polygon([mid + n, mid + n - d * 3 * reach, mid - n - d * 3 * reach, mid - n]))
        pieces = [g for g in shapely.get_parts(cell.intersection(meshed)) if g.geom_type == "Polygon"]
        cells.append(min(pieces, key=lambda g: g.distance(Point(p))))
    return cells, meshed


def _xy(point) -> str:
    return f"({float(point[0]):.1f}, {float(point[1]):.1f})"


def _cells(vor, iverts=None):
    verts = np.asarray(vor.verts, dtype=float)
    return [Polygon(verts[iv]) for iv in (vor.iverts if iverts is None else iverts)]


def _bad(vor, reference):
    return [k for k, (c, r) in enumerate(zip(_cells(vor), reference, strict=True)) if c.symmetric_difference(r).area > 1e-7 * r.area]


def _faces(iverts):
    return [{frozenset((iv[i], iv[(i + 1) % len(iv)])) for i in range(len(iv))} for iv in iverts]


def _mf6(vor, iverts, ws: pathlib.Path):
    """Steady confined flow from the east arm's end (10 ft) to the north arm's end (0 ft)."""
    points = np.asarray(vor.points)
    props = flopy.utils.cvfdutil.get_disv_gridprops(np.asarray(vor.verts), iverts, xcyc=points)
    sim = flopy.mf6.MFSimulation(sim_name="l", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-9, inner_dvclose=1e-10)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="l", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, nlay=1, top=10.0, botm=[0.0], **props)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    east = [[(0, k), 10.0] for k, (x, _) in enumerate(points) if x > 1000 - 1e-6]
    north = [[(0, k), 0.0] for k, (_, y) in enumerate(points) if y > 1000 - 1e-6]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=east + north)
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="l.hds", budget_filerecord="l.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    if not ok:
        raise RuntimeError(f"MF6 failed in {ws}")
    chd = gwf.output.budget().get_data(text="CHD")[0]
    grb = MfGrdFile(str(ws / "l.disv.grb"))
    return dict(ia=grb.ia, ja=grb.ja, flow=float(chd["q"][chd["q"] > 0].sum()),
                head=gwf.output.head().get_data().ravel())


def main() -> None:
    warnings.simplefilter("ignore")
    root = pathlib.Path(tempfile.mkdtemp())

    # 1. Mode 1 on the L.
    vor, tri = _grid(root, L_SHAPE, 5000, "l5000")
    reference, meshed = _reference(vor, tri)
    bad = _bad(vor, reference)
    union = shapely.union_all(_cells(vor))
    print(f"1. L, max_area 5000: {vor.ncpl} cells, and the ones that differ from the reference: "
          + ", ".join(f"cell {k} at {_xy(vor.points[k])}" for k in bad))
    print(f"   hole inside the boundary {meshed.difference(union).area:.1f} ft2, cell area outside it "
          f"{union.difference(meshed).area:.1f} ft2")
    corner = bad[0]
    ring = np.asarray(vor.verts)[vor.iverts[corner]]
    offset = ring.mean(axis=0) - vor.points[corner]
    print(f"   the corner cell's {len(ring)} vertices average ({offset[0]:+.3f}, {offset[1]:+.3f}) ft from "
          f"its node: north-east of it, in the notch")
    verts = np.asarray(vor.verts)
    for label, centre in (("their mean, as tri2vor does", ring.mean(axis=0)),
                          ("a point 0.01 ft south-west of the node", vor.points[corner] - 0.01)):
        order = sorted(vor.iverts[corner], key=lambda v: math.atan2(centre[1] - verts[v][1], verts[v][0] - centre[0]))
        miss = Polygon(verts[order]).symmetric_difference(reference[corner]).area
        print(f"   the same {len(ring)} vertices sorted by angle about {label}: {miss:.1f} ft2 from the reference")
    faces = _faces(vor.iverts)
    lost = [j for j in range(vor.ncpl) if j != corner
            and reference[corner].boundary.intersection(reference[j].buffer(1e-6)).length > 1e-4
            and not faces[corner] & faces[j]]
    print(f"   cells whose reference shares a face with it but whose polygon shares no edge: {lost}")
    fixed = VoronoiGridPlus(tri)
    funion = shapely.union_all(_cells(fixed))
    print(f"   repaired: hole {meshed.difference(funion).area:.2g} ft2, outside {funion.difference(meshed).area:.2g}; "
          f"{fixed.tessellation_report.summary()}")

    # 2. No myflopy code.
    failed = total = 0
    for theta in (240, 270, 300, 330):
        for rot in (0, 17, 41):
            domain = _wedge_notch(theta, rot)
            ws = root / f"flopy_{theta}_{rot}"
            ws.mkdir()
            ftri = Triangle(model_ws=str(ws), angle=20, maximum_area=5000.0)
            ftri.add_polygon(domain)
            ftri.add_region(domain.representative_point().coords[0], 0, maximum_area=5000.0)
            ftri.build()
            fvor = VoronoiGrid(ftri)
            fref, _ = _reference(fvor, ftri)
            apex = affinity.rotate(Point(400, 400), rot, origin=(500, 500))
            at_apex = [k for k in _bad(fvor, fref) if Point(fvor.points[k]).distance(apex) < 1e-6]
            failed += bool(at_apex)
            total += 1
    print(f"2. flopy's own Triangle (angle 20) + VoronoiGrid, reflex corners of 240-330 deg: "
          f"{failed} of {total} corners come out wrong")

    # 3. How general, on the default path.
    print("3. GridSpec.voronoi defaults, max_area 5000, one reflex corner, 3 rotations each:")
    for theta in (210, 240, 270, 300, 330):
        wrong, holes, after = 0, [], 0.0
        for rot in (0, 17, 41):
            gvor, gtri = _grid(root, _wedge_notch(theta, rot), 5000, f"w{theta}_{rot}")
            gref, gmeshed = _reference(gvor, gtri)
            apex = affinity.rotate(Point(400, 400), rot, origin=(500, 500))
            hit = [k for k in _bad(gvor, gref) if Point(gvor.points[k]).distance(apex) < 1e-6]
            if hit:
                wrong += 1
                gap = gmeshed.difference(shapely.union_all(_cells(gvor)))
                near = sum(g.area for g in shapely.get_parts(gap) if g.distance(apex) < 1e-6)
                holes.append(100 * near / gref[hit[0]].area)
            gfixed = shapely.union_all(_cells(VoronoiGridPlus(gtri)))
            after += gmeshed.difference(gfixed).area + gfixed.difference(gmeshed).area
        tail = f"; the hole is {', '.join(f'{h:.0f}%' for h in holes)} of the corner cell" if holes else ""
        print(f"   {theta} deg: {wrong} of 3 wrong{tail}; repaired, {after:.2g} ft2 of gap and overshoot in all")

    # 4. Mode 2 on the L at max_area 10000.
    vor10, tri10 = _grid(root, L_SHAPE, 10000, "l10000")
    reference10, meshed10 = _reference(vor10, tri10)
    bad10 = _bad(vor10, reference10)
    union10 = shapely.union_all(_cells(vor10))
    reflex = _reflex_corners(meshed10)
    kinds = ["reflex corner" if (round(vor10.points[k][0], 6), round(vor10.points[k][1], 6)) in reflex
             else ("interior" if meshed10.exterior.distance(Point(vor10.points[k])) > 1e-6 else "boundary")
             for k in bad10]
    print("4. L, max_area 10000: bad cells " + ", ".join(f"{k} ({kind} at {_xy(vor10.points[k])})" for k, kind in zip(bad10, kinds, strict=True)))
    print(f"   hole {meshed10.difference(union10).area:,.1f} ft2 of {meshed10.area:,.0f}; cell area outside the "
          f"boundary {union10.difference(meshed10).area:.1f} ft2")
    qhull = Voronoi(tri10.verts)
    on_edge = [v for v in qhull.vertices if meshed10.exterior.distance(Point(v)) < 1e-6]
    kept = point_in_polygon(np.array([[v[0]] for v in on_edge]), np.array([[v[1]] for v in on_edge]),
                            list(tri10._polygons[0])).ravel() if on_edge else []
    print(f"   Voronoi vertices exactly on the boundary: {len(on_edge)}, kept by point_in_polygon: "
          f"{int(np.sum(kept))} -- dropped at " + ", ".join(_xy(v) for v, k in zip(on_edge, kept, strict=True) if not k))

    # 5. What MF6 makes of it.
    print("5. MF6, steady flow around the L (CHD 10 ft on the east arm's end, 0 ft on the north arm's):")
    for label, g, t in (("5000", vor, tri), ("10000", vor10, tri10)):
        built = _mf6(g, g.iverts, root / f"mf6_built_{label}")
        repaired = VoronoiGridPlus(t)
        fixed = _mf6(repaired, repaired.iverts, root / f"mf6_fixed_{label}")
        restored = (len(fixed["ja"]) - len(built["ja"])) // 2
        dh = np.abs(fixed["head"] - built["head"])
        print(f"   max_area {label}: {restored} face(s) MF6 does not connect on FloPy's cells; flow "
              f"{built['flow']:.4f} on them vs {fixed['flow']:.4f} ft3/d repaired "
              f"({100 * (built['flow'] - fixed['flow']) / fixed['flow']:+.2f}%), largest head change {dh.max():.3f} ft")
        if label == "5000":
            neighbours = sorted(int(j) for j in built["ja"][built["ia"][corner] + 1: built["ia"][corner + 1]])
            print(f"   the corner cell on FloPy's cells connects to {neighbours}")

    # 6. What myflopy's geometry makes of it.
    polys = vor.gdf_vorPolys
    print(f"6. myflopy's own polygon adjacency (the DISU path) gives the corner cell "
          f"{sorted(int(j) for j in find_adjacent_polygons(polys)[corner])}")
    hole = meshed.difference(union)
    inside = hole.representative_point()
    print(f"   a point inside the hole: cell_shares {cell_shares(polys, inside)}, miss_reason "
          f"{miss_reason(polys, inside)!r}")
    line, square = LineString([(380, 391), (420, 391)]), box(390, 380, 410, 400)
    for name, feature in (("a 40 ft line through the corner", line), ("a 20 x 20 ft polygon at it", square)):
        total_share = sum(share for _, share in cell_shares(polys, feature))
        print(f"   {name}: its shares sum to {total_share:.3f} -- mf.Spread's default 'clip' writes that much")
    for label, g, t in (("5000", vor, tri), ("10000", vor10, tri10)):
        fixed = VoronoiGridPlus(t)
        extra = sorted(set(get_grid_edge_cells(g)) - set(get_grid_edge_cells(fixed)))
        boundary = Polygon(L_SHAPE).exterior
        distance = min((boundary.distance(c) for c in _cells(g, [g.iverts[k] for k in extra])), default=float("nan"))
        print(f"   max_area {label}: get_grid_edge_cells adds cells {extra} -- the nearest of them "
              f"{distance:.1f} ft from the boundary")
    vor2000, _ = _grid(root, L_SHAPE, 2000, "l2000")
    for old, new, label in ((vor, vor2000, "5000 -> 2000"), (vor10, vor2000, "10000 -> 2000")):
        print(f"   cell_coverage {label}, two grids of the ONE boundary: {cell_coverage(old.gdf_vorPolys, new.gdf_vorPolys).summary()}")
        overlay = cell_overlay(old.gdf_vorPolys, new.gdf_vorPolys)
        carried, reached = overlay.reduce(np.full(overlay.n_old, 7.0))
        print(f"   ...and a constant 7.0 remapped across it (CellOverlay.reduce, what with_grid(remap=True) "
              f"uses): {int((~reached).sum())} new cells unreached, largest departure "
              f"{np.nanmax(np.abs(carried - 7.0)):.1e}")
    for label, g in (("2000", vor2000), ("5000", vor), ("10000", vor10)):
        print(f"   max_area {label}: total cell area {g.gdf_vorPolys.area.sum() - L_SHAPE.area:+,.1f} ft2 against the boundary")

    # 7. A boundary with a hole.
    ringed = Polygon(box(0, 0, 1000, 1000).exterior.coords, [box(300, 300, 600, 600).exterior.coords])
    try:
        _grid(root, ringed, 5000, "ringed")
        print("7. a boundary with an interior ring: built")
    except Exception as exc:  # noqa: BLE001 - reported, whatever it is
        print(f"7. a boundary with an interior ring: {type(exc).__name__}: {exc}")
    parts = [box(0, 0, 1000, 1000), box(1500, 0, 2500, 1000)]
    path = root / "two_parts.gpkg"
    gpd.GeoDataFrame({"name": ["west", "east"]}, geometry=parts, crs=2927).to_file(path)
    pvor, ptri = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=path), max_area=5000.0,
                                     tessellation="flopy").resolve(workspace=root / "two_parts", return_triangle=True)

    def per_part(points):
        return " and ".join(f"{int(np.sum([Point(p).within(part.buffer(1e-6)) for p in points]))}" for part in parts)

    print(f"   a boundary of two features, 1,000,000 ft2 each: Triangle meshes {per_part(ptri.verts)} nodes in "
          f"them; FloPy's cells cover {per_part(pvor.points)} of them, the repaired grid "
          f"{per_part(VoronoiGridPlus(ptri).points)}")

    # 8. The boundary the default profile meshes.
    xs = np.linspace(0, 2000, 81)
    corridor = LineString(np.column_stack([xs, 300 * np.sin(xs / 2000 * 3 * np.pi)])).buffer(150, quad_segs=8)
    cvor, ctri = _grid(root, corridor, 5000, "corridor")
    cmeshed = Polygon(ctri._polygons[0])
    cunion = shapely.union_all(_cells(cvor))
    print(f"8. a {corridor.area:,.0f} ft2 corridor: the meshed polygon differs from it by "
          f"{cmeshed.symmetric_difference(corridor).area:,.0f} ft2 and keeps {len(_reflex_corners(cmeshed))} of "
          f"its {len(_reflex_corners(corridor))} reflex vertices; the cells leave {corridor.difference(cunion).area:,.0f} "
          f"ft2 of it uncovered and put {cunion.difference(corridor).area:,.0f} ft2 outside it")

    # 9. The cost: the repair's own CPU time inside a resolve, as a share of the resolve.
    #    A ratio of two timings taken together is what a busy machine distorts least; two
    #    separate resolves compared wall clock to wall clock measured -12% at 2,136 cells.
    import os

    import myflopy.modflow.mf6.grid.voronoi as voronoi_module

    def cpu():
        t = os.times()
        return t.user + t.system + t.children_user + t.children_system

    spent = {"cpu": 0.0, "calls": 0}
    original = voronoi_module.tessellate

    def timed(*args, **kwargs):
        begin = cpu()
        try:
            return original(*args, **kwargs)
        finally:
            spent["cpu"] += cpu() - begin
            spent["calls"] += 1

    voronoi_module.tessellate = timed
    try:
        for area in (5000, 500, 100):
            spent.update(cpu=0.0, calls=0)
            begin = cpu()
            g, _ = _grid(root, L_SHAPE, area, f"cost_{area}", "dual")
            total = cpu() - begin
            print(f"9. GridSpec.resolve, the L at max_area {area} ({g.ncpl} cells): {total:.2f} CPU s, of which the "
                  f"repair {spent['cpu']:.2f} ({100 * spent['cpu'] / (total - spent['cpu']):.0f}% on top of the rest) "
                  f"over {spent['calls']} grid builds")
    finally:
        voronoi_module.tessellate = original

if __name__ == "__main__":
    main()
