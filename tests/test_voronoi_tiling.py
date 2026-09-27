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
    P = P - P.min(axis=0)  # the formula below squares coordinates: at state plane, 1e12

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


def test_a_package_from_the_first_repair_is_told_its_cells_kept_their_numbers(l5000, triangulated_grid):
    """The first version of the repair pinned no record of itself, so a package resolved on
    its grid meets the generic refusal -- which alone says the cells are another mesh's."""

    from types import SimpleNamespace

    from myflopy.grid_identity import GRID_PIN, _discretization_complaints

    package = SimpleNamespace(name="chd", metadata={GRID_PIN: "0123456789abcdef"})
    repaired = _discretization_complaints(None, SimpleNamespace(grid=l5000[0], metadata={}), [package])
    assert any("first version of ledger 316's repair" in c for c in repaired), repaired
    untouched = _discretization_complaints(None, SimpleNamespace(grid=triangulated_grid, metadata={}), [package])
    assert untouched and not any("first version" in c for c in untouched), untouched


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


STAR = [(832.7, 272.9), (612.4, 462.2), (921.4, 419.7), (615.6, 483.0), (786.3, 790.1), (556.0, 600.6),
        (678.2, 912.3), (548.8, 619.9), (619.0, 873.6), (483.2, 612.2), (424.1, 884.4), (447.4, 618.0),
        (314.3, 889.7), (429.1, 601.4), (223.9, 828.8), (422.2, 586.8), (201.7, 831.2), (421.2, 587.3),
        (214.9, 762.5), (415.1, 575.3), (148.3, 653.1), (374.1, 537.9), (110.1, 468.9), (398.8, 425.0),
        (228.7, 147.6), (422.9, 396.5), (324.1, 127.7), (485.8, 376.2), (626.8, 108.6), (576.1, 415.4)]


@pytest.fixture(scope="module")
def star(tmp_path_factory):
    return _grid(tmp_path_factory, Polygon([(1.2e6 + x, 6.5e5 + y) for x, y in STAR]), 10000,
                 profile="baseline", min_angle=21)


def test_a_spiky_boundary_at_state_plane_keeps_every_face(star, tmp_path_factory):
    """Found by a fuzz pass: across 30 thin spikes near x = 1.2e6, Qhull's Delaunay is not
    Triangle's, and cells FloPy drew from both lost faces the first repair did not give back."""

    vor, tri = star
    assert _disv_faces(vor.iverts) == _mesh_faces(vor, tri)
    gap, outside, overlap = _tiling(vor, tri)
    assert gap < 1e-6 and outside < 1e-6 and abs(overlap) < 1e-6
    polygon = Polygon([(1.2e6 + x, 6.5e5 + y) for x, y in STAR])
    old, _ = _grid(tmp_path_factory, polygon, 10000, profile="baseline", min_angle=21, tessellation="flopy")
    assert len(_mesh_faces(old, tri) - _disv_faces(old.iverts)) > 10


def test_every_vertex_is_where_the_mesh_puts_one(star):
    """Faces and coverage cannot see a vertex moved in every cell at once. The first repair
    took FloPy's vertex up to 5% of the circumradius from the circumcentre it stood for;
    put back, every test above passed with cells 0.3 ft off their own shape."""

    vor, tri = star
    P, T = np.asarray(tri.verts, dtype=float), np.asarray(tri.iverts, dtype=int)
    origin = P.min(axis=0)
    local = P - origin
    a, b, c = local[T[:, 0]], local[T[:, 1]], local[T[:, 2]]
    ba, ca = b - a, c - a
    d = 2 * (ba[:, 0] * ca[:, 1] - ba[:, 1] * ca[:, 0])
    b2, c2 = (ba**2).sum(1), (ca**2).sum(1)
    centres = a + np.column_stack([(ca[:, 1] * b2 - ba[:, 1] * c2) / d, (ba[:, 0] * c2 - ca[:, 0] * b2) / d])
    edges = {}
    for t in T:
        for u, v in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edges[(min(u, v), max(u, v))] = edges.get((min(u, v), max(u, v)), 0) + 1
    rim = np.array([k for k, n in edges.items() if n == 1])
    truth = np.vstack([centres, (local[rim[:, 0]] + local[rim[:, 1]]) / 2, local[np.unique(rim)]]) + origin
    used = np.unique(np.concatenate([np.asarray(r) for r in vor.iverts]))
    from scipy.spatial import cKDTree

    off = cKDTree(truth).query(np.asarray(vor.verts, dtype=float)[used])[0]
    assert off.max() <= 2 * 1000 * np.spacing(np.abs(P).max()), off.max()


def _triangle_at_state_plane(tmp_path, nodes):
    """Triangle's own mesh of a 1000 ft square at (1.2e6, 6.5e5), with ``nodes`` as constraints."""

    from flopy.utils.triangle import Triangle

    corner = np.array([1.2e6, 6.5e5])
    tri = Triangle(model_ws=str(tmp_path), maximum_area=2000.0, angle=28.0, nodes=corner + np.asarray(nodes))
    tri.add_polygon(corner + np.array([[0.0, 0.0], [1000.0, 0.0], [1000.0, 1000.0], [0.0, 1000.0]]))
    tri.build()
    return tri


def _resolvable():
    from myflopy.modflow.mf6.grid.tessellation import COINCIDENT_ULPS

    return COINCIDENT_ULPS * float(np.spacing(1.2e6 + 1000.0))


def test_two_nodes_closer_than_1000_ulp_are_refused(tmp_path):
    """Two wells 50 ulp apart: a mesh the certificate ACCEPTS, whose cells between the two
    the merge would collapse. Refused by name, naming how far apart and where."""

    from myflopy.modflow.mf6.grid.tessellation import NotConformingDelaunay, dual_cells

    tri = _triangle_at_state_plane(tmp_path, [[400.0, 600.0], [400.0 + 50 * np.spacing(1.2e6), 600.0]])
    with pytest.raises(NotConformingDelaunay, match="closer than 1000 ulp"):
        dual_cells(tri)


def test_a_node_whose_circumcentres_merge_is_refused_not_dropped(tmp_path):
    """A point 1.24 x 1000 ulp off the south edge passes the node refusal, and its
    triangles' circumcentres all merge onto the edge: its cell was dropped without a word."""

    from myflopy.modflow.mf6.grid.tessellation import NotConformingDelaunay, dual_cells

    tri = _triangle_at_state_plane(tmp_path, [[706.7, 1.24 * _resolvable()]])
    with pytest.raises(NotConformingDelaunay, match="has no cell left"):
        dual_cells(tri)


def test_the_merge_does_not_widen_the_certificate():
    """A boundary triangle whose circumcentre lies 2.1e-7 ft OUTSIDE the domain -- inside
    the merge's 1000 ulp. Certified at that width, the cell beside it crossed itself."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.grid.tessellation import NotConformingDelaunay, dual_cells

    half_disc = SimpleNamespace(
        verts=np.array([(1199980.0, 650000.0), (1200020.0, 650000.0), (1200019.1794925705, 650005.6698380732),
                        (1200015.737745389, 650012.3419357049), (1200014.7466624975, 650013.5105866779),
                        (1200014.43531111, 650013.8427525797), (1200013.3826761066, 650014.8628386959),
                        (1200002.9509658685, 650019.7810971513), (1199986.4341472546, 650014.6958374687)]),
        iverts=[[1, 6, 0], [1, 2, 6], [8, 0, 6], [2, 3, 6], [7, 8, 6], [5, 6, 4], [3, 4, 6]])
    with pytest.raises(NotConformingDelaunay, match="1 boundary edge"):
        dual_cells(half_disc)
    # Exactly Delaunay and Gabriel, so certified: the merge snaps one circumcentre onto a
    # boundary midpoint, which turns a fan triangle of 1e-6 ft2 over -- and every cell is
    # still a valid polygon.
    certified = SimpleNamespace(
        verts=np.array([(1199980.0, 650000.0), (1200020.0, 650000.0), (1200019.0202360365, 650006.1830917775),
                        (1199985.5275980216, 650013.8039699646), (1199983.2061290708, 650010.8612109976)]),
        iverts=[[2, 3, 1], [1, 3, 0], [4, 0, 3]])
    verts, iverts, _ = dual_cells(certified)
    assert all(Polygon(verts[ring]).is_valid for ring in iverts)


def test_a_grid_whose_reprojected_cells_do_not_union_still_has_a_centroid():
    """17 of 156,582 cells at (5e6, 5e6) ft crossed themselves once reprojected, every one
    valid where drawn, and GEOS refused the union: the grid did not build."""

    from types import SimpleNamespace

    import geopandas as gpd

    from myflopy.modflow.mf6.grid.geometry import get_grid_centroid

    class Refusing:
        def union_all(self):
            raise shapely.errors.GEOSException("TopologyException: unable to assign free hole to a shell")

    cells = gpd.GeoDataFrame(geometry=[box(1.2e6, 6.5e5, 1.2e6 + 100, 6.5e5 + 100),
                                       box(1.2e6 + 100, 6.5e5, 1.2e6 + 300, 6.5e5 + 100)], crs=2927)
    vor = SimpleNamespace(gdf_latlon=Refusing(), gdf_vorPolys=cells, crs_latlon="EPSG:4326")
    expected = gpd.GeoSeries([shapely.Point(1.2e6 + 150, 6.5e5 + 50)], crs=2927).to_crs(4326).iloc[0]
    assert get_grid_centroid(vor).distance(expected) < 1e-9



def test_circumcentres_closer_than_the_coordinates_resolve_are_one_vertex():
    """At (5e6, 5e6) two circumcentres 2e-9 ft apart were kept as two vertices -- an edge of
    two ulps, which a reprojection turned over. Merging only within 1e-9 of the edge kept them."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.grid.tessellation import dual_cells

    # two triangles on edge AB, cocircular but for D sitting 5e-8 ft outside the circle
    local = np.array([[0.0, 0.0], [10.0, 0.0], [5.0, 5.0], [5.0, -5.0 - 5e-8]])
    mesh = SimpleNamespace(verts=local + [1.2e6, 6.5e5], iverts=[[0, 1, 2], [0, 3, 1]])
    verts, _, _ = dual_cells(mesh)
    near_centre = np.hypot(*(verts - [1.2e6 + 5.0, 6.5e5]).T) < 1e-6
    assert near_centre.sum() == 1


def test_one_flopy_vertex_is_never_given_to_two_dual_vertices():
    from myflopy.modflow.mf6.grid.tessellation import _onto_flopy

    dual = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [1e-14, 0.0], [0.0, -1.0]])
    rings = [[0, 1, 2, 4], [3, 2, 4]]  # vertices 0 and 3 sit 1e-14 apart, inside the tolerance
    flopy = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, -1.0], [0.0, -1.0]])
    verts, iverts, _, appended = _onto_flopy(dual, rings, flopy, [[0, 1, 2, 3], [0, 2, 3]], np.array([False, False]))
    assert appended == 1 and iverts[1][0] == 4 and len(verts) == 5
