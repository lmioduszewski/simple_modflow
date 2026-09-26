"""Voronoi cells tile their domain: FloPy's ``tri2vor`` repaired (ledger 316).

The TRUTH here is computed from Triangle's mesh directly, not from the repair's own
builder: the meshed region is the union of Triangle's triangles, and the faces are the
mesh edges whose dual segment has length -- an interior edge between two distinct
circumcentres, a boundary edge whose triangle's circumcentre is not its midpoint. A test
comparing the repair against ``dual_cells`` alone would compare the module with itself.
"""

from __future__ import annotations

import math
import pickle

import numpy as np
import pytest
import shapely
from shapely.geometry import Polygon, box

L_SHAPE = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])


def _grid(tmp_path_factory, domain, max_area, **options):
    import geopandas as gpd

    import myflopy as mf

    root = tmp_path_factory.mktemp("tiling")
    gpd.GeoDataFrame({"name": ["domain"]}, geometry=[domain], crs=2927).to_file(root / "b.gpkg")
    spec = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=root / "b.gpkg"), max_area=float(max_area), **options
    )
    return spec.resolve(workspace=root / "grid", return_triangle=True)


@pytest.fixture(scope="module")
def l5000(tmp_path_factory):
    return _grid(tmp_path_factory, L_SHAPE, 5000)


@pytest.fixture(scope="module")
def l5000_flopy(tmp_path_factory):
    return _grid(tmp_path_factory, L_SHAPE, 5000, tessellation="flopy")


@pytest.fixture(scope="module")
def l10000(tmp_path_factory):
    return _grid(tmp_path_factory, L_SHAPE, 10000)


@pytest.fixture(scope="module")
def l10000_flopy(tmp_path_factory):
    return _grid(tmp_path_factory, L_SHAPE, 10000, tessellation="flopy")


def _cells(vor):
    verts = np.asarray(vor.verts, dtype=float)
    return [Polygon(verts[iv]) for iv in vor.iverts]


def _tiling(vor, tri):
    """(gap, outside, overlap) of the cells against the union of Triangle's triangles."""

    P, T = np.asarray(tri.verts), np.asarray(tri.iverts)
    meshed = shapely.union_all([Polygon(P[t]) for t in T])
    cells = _cells(vor)
    union = shapely.union_all(cells)
    return (meshed.difference(union).area, union.difference(meshed).area,
            sum(c.area for c in cells) - union.area)


def _disv_faces(iverts):
    """Cell pairs sharing a vertex PAIR -- how MF6 DISV finds a face."""

    owner, faces = {}, set()
    for k, ring in enumerate(iverts):
        for i in range(len(ring)):
            edge = frozenset((int(ring[i]), int(ring[(i + 1) % len(ring)])))
            if edge in owner and owner[edge] != k:
                faces.add((min(owner[edge], k), max(owner[edge], k)))
            owner[edge] = k
    return faces


def _mesh_faces(vor, tri):
    """The faces the mesh says exist, as cell pairs, from Triangle's triangles alone."""

    P, T = np.asarray(tri.verts, dtype=float), np.asarray(tri.iverts, dtype=int)
    cell_of = {(float(x), float(y)): k for k, (x, y) in enumerate(np.asarray(vor.points))}
    node_cell = [cell_of.get((float(x), float(y))) for x, y in P]

    def circumcentre(t):
        (ax, ay), (bx, by), (cx, cy) = P[t]
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
        uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
        return np.array([ux, uy])

    edges = {}
    for t in T:
        for u, v in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edges.setdefault((min(u, v), max(u, v)), []).append(t)
    faces = set()
    for (u, v), ts in edges.items():
        a, b = node_cell[u], node_cell[v]
        if a is None or b is None:
            continue
        length = math.dist(P[u], P[v])
        if len(ts) == 2:
            dual = np.linalg.norm(circumcentre(ts[0]) - circumcentre(ts[1]))
        else:
            dual = np.linalg.norm(circumcentre(ts[0]) - (P[u] + P[v]) / 2)
        if dual > 1e-7 * length:
            faces.add((min(a, b), max(a, b)))
    return faces


def test_a_reflex_corner_cell_tiles_its_corner(l5000, l5000_flopy):
    """Mode 1: the L's corner cell was drawn out of order -- a 77.4 ft2 hole, 43.9 outside."""

    vor, tri = l5000
    gap, outside, overlap = _tiling(vor, tri)
    assert gap < 1e-6 and outside < 1e-6 and abs(overlap) < 1e-6, (gap, outside, overlap)
    before = _tiling(*l5000_flopy)
    assert before[0] > 70 and before[1] > 40  # what FloPy's cells leave


def test_a_vertex_on_the_boundary_keeps_its_faces(l10000, l10000_flopy):
    """Mode 2: a Voronoi vertex ON a right-hand edge was dropped; seven faces went with it."""

    vor, tri = l10000
    assert _disv_faces(vor.iverts) == _mesh_faces(vor, tri)
    assert _tiling(vor, tri)[0] < 1e-6
    old_vor, old_tri = l10000_flopy
    assert len(_mesh_faces(old_vor, old_tri) - _disv_faces(old_vor.iverts)) == 7
    assert "missing-vertex" in set(vor.tessellation_report.rewritten.values())


def test_the_report_names_the_cells_and_the_old_pin(l5000, l10000, l5000_flopy):
    from myflopy.grid_identity import grid_pin

    vor, _ = l5000
    report = vor.tessellation_report
    assert report.rewritten == {3: "corner-order"}
    assert not report.adopted and report.appended == 0
    assert report.legacy_pin == grid_pin(l5000_flopy[0]) != grid_pin(vor)
    assert sorted(l10000[0].tessellation_report.rewritten) == [3, 6, 12, 15]


def test_every_cell_flopy_got_right_keeps_its_ring(l5000, l5000_flopy, l10000, l10000_flopy):
    """Only broken cells are rewritten, index for index, so an unbroken grid keeps its pin."""

    for (new, _), (old, _) in ((l5000, l5000_flopy), (l10000, l10000_flopy)):
        changed = [k for k in range(old.ncpl) if list(new.iverts[k]) != list(old.iverts[k])]
        assert changed == sorted(new.tessellation_report.rewritten)
        assert np.array_equal(new.points, old.points)


def test_a_grid_with_nothing_wrong_is_byte_identical(triangulated_grid, triangulated_grid_finer):
    from myflopy.grid_identity import grid_pin

    for vor in (triangulated_grid, triangulated_grid_finer):
        assert not vor.tessellation_report.changed
        assert vor.tessellation_report.legacy_pin is None
        assert grid_pin(vor) is not None


def test_the_repair_does_not_move_the_mesh(l10000, l10000_flopy):
    """The smoothing loop builds grids too; the final Triangle nodes are bit-identical."""

    assert np.array_equal(np.asarray(l10000[1].verts), np.asarray(l10000_flopy[1].verts))


def test_a_sharp_slit_tiles(tmp_path):
    """A 355 deg corner: cells beside a slit took vertices across it; FloPy left gaps and overlaps."""

    from myflopy.modflow.mf6.grid.triangle import TriangleGrid
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    alpha, far = math.radians(5.0), 1.0e4
    ray = lambda a: (400 + far * math.cos(a), 400 + far * math.sin(a))  # noqa: E731
    wedge = Polygon([(400, 400), ray(math.radians(45) - alpha / 2), ray(math.radians(45) + alpha / 2)])
    domain = box(0, 0, 1000, 1000).difference(wedge)
    tri = TriangleGrid(model_ws=str(tmp_path), angle=32)
    tri.set_domain_polygon(domain, max_area=10000.0)
    tri.build_mesh(profile="balanced")
    vor = VoronoiGridPlus(tri)
    gap, outside, overlap = _tiling(vor, tri)
    assert gap < 1e-6 and outside < 1e-6 and abs(overlap) < 1e-6
    assert _disv_faces(vor.iverts) == _mesh_faces(vor, tri)
    before = _tiling(VoronoiGridPlus(tri, tessellation="flopy"), tri)
    assert before[0] > 1.0


def test_twin_vertex_ids_share_one_face(tmp_path):
    """Mode 2b: FloPy kept a boundary vertex AND its midpoint twin; MF6 lost 20 of 97 faces."""

    from flopy.utils.triangle import Triangle

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    tri = Triangle(model_ws=str(tmp_path), angle=20, maximum_area=20000.0)
    tri.add_polygon(box(0, 0, 1000, 1000))
    tri.add_region((500.0, 500.0), 0, maximum_area=20000.0)
    tri.build()
    vor = VoronoiGridPlus(tri)
    assert _disv_faces(vor.iverts) == _mesh_faces(vor, tri)
    assert "twin-vertex" in set(vor.tessellation_report.rewritten.values())
    old = VoronoiGridPlus(tri, tessellation="flopy")
    assert len(_disv_faces(old.iverts)) < len(_disv_faces(vor.iverts))


def test_a_point_outside_the_boundary_gets_no_cell(tmp_path):
    """A constraint point 2 ft into the L's notch: FloPy kept a cell centred outside the domain."""

    from myflopy.modflow.mf6.grid.triangle import TriangleGrid
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    tri = TriangleGrid(model_ws=str(tmp_path), angle=32)
    tri.set_domain_polygon(L_SHAPE, max_area=10000.0)
    tri.add_points([(401.0, 402.0)])
    tri.build_mesh(profile="clean")
    vor = VoronoiGridPlus(tri)
    assert all(L_SHAPE.covers(shapely.Point(p)) for p in np.asarray(vor.points))
    gap, outside, overlap = _tiling(vor, tri)
    assert outside < 1e-6 and abs(overlap) < 1e-6
    assert vor.tessellation_report.adopted
    old = VoronoiGridPlus(tri, tessellation="flopy")
    assert not all(L_SHAPE.covers(shapely.Point(p)) for p in np.asarray(old.points))


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_tessellation_flopy_is_tri2vor_exactly(l5000, l5000_flopy):
    from flopy.utils.voronoi import tri2vor

    _, tri = l5000
    _, iverts, _ = tri2vor(tri)
    assert [list(map(int, r)) for r in l5000_flopy[0].iverts] == [list(map(int, r)) for r in iverts]
    assert l5000_flopy[0].tessellation_report is None


def test_an_unknown_tessellation_is_refused(l5000):
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    with pytest.raises(ValueError, match="not a tessellation"):
        VoronoiGridPlus(l5000[1], tessellation="exact")


def test_a_grid_rebuilt_from_arrays_is_never_repaired(l5000_flopy):
    """A snapshot, a pickle or a reloaded run is the grid that ran -- defects and all."""

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    old, _ = l5000_flopy
    rebuilt = VoronoiGridPlus(verts=np.asarray(old.verts), iverts=old.iverts, xcyc=np.asarray(old.points))
    assert [list(r) for r in rebuilt.iverts] == [list(r) for r in old.iverts]
    assert rebuilt.tessellation_report is None
    thawed = pickle.loads(pickle.dumps(old))
    assert [list(r) for r in thawed.iverts] == [list(r) for r in old.iverts]


def test_a_grid_pickled_before_the_repair_reads_no_report(l5000):
    vor, _ = l5000
    state = vor.__getstate__()
    state.pop("tessellation_report")
    fresh = type(vor).__new__(type(vor))
    fresh.__setstate__(state)
    assert fresh.tessellation_report is None


def test_a_mesh_that_is_not_delaunay_is_refused():
    from types import SimpleNamespace

    from myflopy.modflow.mf6.grid.tessellation import NotConformingDelaunay, dual_cells

    obtuse = SimpleNamespace(verts=np.array([[0.0, 0.0], [2.0, 0.0], [1.0, 0.2]]), iverts=[[0, 1, 2]])
    with pytest.raises(NotConformingDelaunay, match="-D"):
        dual_cells(obtuse)


def test_the_refusal_says_the_package_predates_the_repair(l5000, l5000_flopy):
    from types import SimpleNamespace

    from myflopy.grid_identity import GRID_PIN, _discretization_complaints, grid_pin

    vor, _ = l5000
    package = SimpleNamespace(name="chd", metadata={GRID_PIN: grid_pin(l5000_flopy[0])})
    complaints = _discretization_complaints(None, SimpleNamespace(grid=vor, metadata={}), [package])
    assert any("BEFORE ledger 316's repair" in c and "tessellation='flopy'" in c for c in complaints), complaints


@pytest.mark.slow
def test_mf6_connects_every_face_and_the_flow_comes_back(l10000, l10000_flopy, tmp_path):
    """The L at max_area 10000: MF6's own .grb, and the flow the lost faces cost (5.12%)."""

    import flopy
    from flopy.mf6.utils import MfGrdFile

    def run(vor, name):
        points = np.asarray(vor.points)
        props = flopy.utils.cvfdutil.get_disv_gridprops(np.asarray(vor.verts), vor.iverts, xcyc=points)
        sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(tmp_path / name), exe_name="mf6")
        flopy.mf6.ModflowTdis(sim)
        flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-9, inner_dvclose=1e-10)
        gwf = flopy.mf6.ModflowGwf(sim, modelname="l", save_flows=True)
        flopy.mf6.ModflowGwfdisv(gwf, nlay=1, top=10.0, botm=[0.0], **props)
        flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=0)
        flopy.mf6.ModflowGwfic(gwf, strt=5.0)
        east = [[(0, k), 10.0] for k, (x, _) in enumerate(points) if x > 1000 - 1e-6]
        north = [[(0, k), 0.0] for k, (_, y) in enumerate(points) if y > 1000 - 1e-6]
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=east + north)
        flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="l.cbc", saverecord=[("BUDGET", "ALL")])
        sim.write_simulation(silent=True)
        assert sim.run_simulation(silent=True)[0]
        chd = gwf.output.budget().get_data(text="CHD")[0]
        return MfGrdFile(str(tmp_path / name / "l.disv.grb")), float(chd["q"][chd["q"] > 0].sum())

    vor, tri = l10000
    grb, flow = run(vor, "repaired")
    assert (grb.nja - vor.ncpl) // 2 == len(_mesh_faces(vor, tri))
    _, old_flow = run(l10000_flopy[0], "flopy")
    assert old_flow < 0.96 * flow
