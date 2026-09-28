"""Five ``GridSpec.voronoi`` inputs a fuzz pass found failing (ledger 316, resolved 2026-09-28).

Each fix changes only inputs that fail, or silently misbehave, today: a grid that builds
today keeps its Triangle input and its ``grid_pin``. The tests below build the failing input
AND a working neighbour, so reverting a fix fails one and a fix that reached too far fails
the other.
"""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Polygon, box

L_SHAPE = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])


def _voronoi(root, boundary, *, refinement=(), breaklines=(), points=(), **options):
    """``(vor, tri)`` for ``boundary`` in EPSG:2927, each source written to ``root``."""

    import geopandas as gpd

    import myflopy as mf

    root.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame({"n": [0]}, geometry=[boundary], crs=2927).to_file(root / "b.gpkg")
    kwargs = {}
    if refinement:
        gpd.GeoDataFrame({"area": [a for _, a in refinement]}, geometry=[g for g, _ in refinement],
                         crs=2927).to_file(root / "r.gpkg")
        kwargs["refinement"] = mf.GeoPackageSourceSpec(path=root / "r.gpkg", fields={"area": "area"})
    if breaklines:
        gpd.GeoDataFrame({"n": list(range(len(breaklines)))}, geometry=list(breaklines), crs=2927).to_file(root / "l.gpkg")
        kwargs["breaklines"] = [mf.GeoPackageSourceSpec(path=root / "l.gpkg")]
    if points:
        from shapely.geometry import Point

        gpd.GeoDataFrame({"n": list(range(len(points)))}, geometry=[Point(p) for p in points], crs=2927).to_file(root / "p.gpkg")
        kwargs["points"] = [mf.GeoPackageSourceSpec(path=root / "p.gpkg")]
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "b.gpkg"), **kwargs, **options)
    return spec.resolve(workspace=root / "grid", return_triangle=True)


def _triangle_areas(tri):
    P, T = np.asarray(tri.verts, dtype=float), np.asarray(tri.iverts, dtype=np.int64)
    a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
    cross = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    return (a + b + c) / 3.0, 0.5 * np.abs(cross)


# --------------------------------------------------------------------------------------------
# 1. Two rings sharing an exact point: one Triangle node, not two.


@pytest.mark.parametrize("profile", ["balanced", "clean", "baseline"])
def test_adjacent_refinement_zones_build(tmp_path, profile):
    """Two zones along one contact: FloPy wrote the shared corners twice and tri2vor raised
    its bare "duplicate points" Exception on every profile. Each zone keeps its max_area."""

    tiles = [box(100, 100, 200, 200), box(200, 100, 300, 200)]
    vor, tri = _voronoi(tmp_path, L_SHAPE, refinement=[(t, 2500.0) for t in tiles], max_area=10000.0,
                        min_angle=28, profile=profile)
    centroid, area = _triangle_areas(tri)
    for tile in tiles:
        x0, y0, x1, y1 = tile.bounds
        inside = (centroid[:, 0] > x0) & (centroid[:, 0] < x1) & (centroid[:, 1] > y0) & (centroid[:, 1] < y1)
        assert inside.any() and area[inside].max() <= 2500.0 * (1 + 1e-9)
    assert len({(float(x), float(y)) for x, y in np.asarray(tri.verts)}) == len(tri.verts)


def test_a_region_corner_on_the_sizing_ring_builds(tmp_path):
    """The fuzz case: box(50, 50, 300, 300) at max_area 10000, whose corner is exactly
    0.5 * sqrt(max_area) from two edges -- ON the corner of the background sizing ring."""

    vor, tri = _voronoi(tmp_path, L_SHAPE, refinement=[(box(50, 50, 300, 300), 2500.0)], max_area=10000.0,
                        min_angle=28)
    assert vor.ncpl > 0 and not vor.tessellation_report.adopted


def test_input_with_nothing_shared_is_flopys_own_bytes(tmp_path):
    """FloPy's writers run unchanged when no point coincides: every grid that builds today
    gets the same Triangle input. With a shared point, it is written once."""

    from flopy.utils.triangle import Triangle

    _, plain = _voronoi(tmp_path / "plain", L_SHAPE, refinement=[(box(60, 60, 300, 300), 2500.0)], max_area=10000.0)
    _, shared = _voronoi(tmp_path / "shared", L_SHAPE, refinement=[(box(50, 50, 300, 300), 2500.0)], max_area=10000.0)
    for name in ("_write_nodefile", "_write_polyfile"):
        getattr(plain, name)(tmp_path / "ours")
        getattr(Triangle, name)(plain, tmp_path / "flopys")
        assert (tmp_path / "ours").read_bytes() == (tmp_path / "flopys").read_bytes()
    shared._write_nodefile(tmp_path / "ours")
    Triangle._write_nodefile(shared, tmp_path / "flopys")
    ours = int((tmp_path / "ours").read_text().split()[0])
    assert ours == int((tmp_path / "flopys").read_text().split()[0]) - 1
    shared.additional_args = ["-j"]  # Triangle already drops the copy under -j: keep its numbering
    assert shared._shared_nodes() is None


# --------------------------------------------------------------------------------------------
# 2. A max_area the domain cannot hold: the sizing region is a point, not a 0.14 ft band.


def test_a_max_area_too_large_for_the_domain_gives_a_coarse_mesh_and_warns(tmp_path):
    """No point of the L is 500 ft from its boundary, so its sizing region was the whole L,
    clipped to a ring 0.14 ft inside it -- and Triangle filled that band with 149,385 cells
    (smallest 0.00019 ft2) whatever max_area was. Its point alone gives what was asked."""

    with pytest.warns(UserWarning, match="larger than the domain 'domain' can hold"):
        vor, tri = _voronoi(tmp_path, L_SHAPE, max_area=1e6)
    assert vor.ncpl < 50
    assert [r.segments for r in tri._prepared_regions if r.source == "domain"] == [False]


def test_a_max_area_the_domain_holds_keeps_its_sizing_ring(tmp_path):
    """Just below the threshold the inset is a thin sliver and the mesh is fine-grained
    (2,851 cells) -- unchanged on purpose: it builds today, and this pin is its pin."""

    from myflopy.grid_identity import grid_pin

    vor, tri = _voronoi(tmp_path, L_SHAPE, max_area=1.58e5)
    assert all(r.segments for r in tri._prepared_regions)
    assert (vor.ncpl, grid_pin(vor)) == (2851, "84a3746283f05e8d")


# --------------------------------------------------------------------------------------------
# 3. Breaklines: a band piece too short to resample, duplicated rows, and a band resampled
#    to nothing.

LINE_B = [(400.0, 581.5845702824), (509.9330967965, 172.0865051051)]  # starts on the notch edge, heads out


def _breaklines(root, lines, areas=None, **options):
    """``GridSpec.voronoi`` of the L with ``lines`` as breaklines on the "clean" profile."""

    import geopandas as gpd
    from shapely.geometry import LineString

    import myflopy as mf

    root.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame({"n": [0]}, geometry=[L_SHAPE], crs=2927).to_file(root / "b.gpkg")
    columns = {"n": list(range(len(lines)))} if areas is None else {"area": list(areas)}
    gpd.GeoDataFrame(columns, geometry=[LineString(c) for c in lines], crs=2927).to_file(root / "l.gpkg")
    source = mf.GeoPackageSourceSpec(path=root / "l.gpkg", fields={} if areas is None else {"area": "area"})
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "b.gpkg"), breaklines=[source],
                               profile="clean", min_angle=23, max_area=3000.0, breakline_max_area=60.0, **options)
    return spec.resolve(workspace=root / "grid", return_triangle=True)


@pytest.mark.parametrize("line", [LINE_B, [(100.0, 200.0), (120.0, 200.0)]], ids=["off-the-edge", "20-ft"])
def test_a_band_piece_too_short_to_resample_builds(tmp_path, line):
    """A piece of the 10 ft band whose outline is no longer than the 120 ft resampling
    spacing -- line B's sliver inside the notch, a 20 ft line -- became a 2-point ring and
    shapely raised "A linearring requires at least 4 coordinates"."""

    vor, tri = _breaklines(tmp_path, [line])
    assert vor.ncpl > 0 and any(r.source == "line" for r in tri._prepared_regions)


def test_an_exact_duplicate_breakline_row_is_dropped_with_a_warning(tmp_path):
    """Two identical rows claimed no unique area and the grid was refused; they mean one
    line. The same line at a DIFFERENT max_area is still refused, as documented."""

    line = [(464.0, 400.0), (428.6, 165.1)]
    with pytest.warns(UserWarning, match="exact duplicate of 'breakline_0'"):
        vor, tri = _breaklines(tmp_path / "same", [line, line[::-1]])
    assert sum(r.source == "line" for r in tri._prepared_regions) == 1
    with pytest.raises(ValueError, match="unique interior area"):
        _breaklines(tmp_path / "different", [line, line], areas=[60.0, 30.0])


def test_a_band_resampled_to_nothing_is_named_not_blamed_on_the_clip(tmp_path):
    """A 40 ft line: its band's outline is under twice the spacing, so resampling leaves no
    area. It was skipped as "empty after clipping to the domain interior", which it is not;
    the grid is unchanged (the no-breakline grid, as before) and the warning names why."""

    from myflopy.grid_identity import grid_pin

    with pytest.warns(UserWarning, match=r"it was removed when resampled at target_segment_length=120") as record:
        vor, _ = _breaklines(tmp_path / "line", [[(100.0, 200.0), (140.0, 200.0)]])
    assert not any("empty after clipping" in str(w.message) for w in record)
    plain, _ = _voronoi(tmp_path / "plain", L_SHAPE, profile="clean", min_angle=23, max_area=3000.0)
    assert grid_pin(vor) == grid_pin(plain)


# --------------------------------------------------------------------------------------------
# 4. "No polygonal area": a feature that misses the boundary is skipped and named, and a
#    boundary part narrower than the 2 ft simplification is kept as drawn.


def test_a_refinement_feature_outside_the_boundary_is_skipped_by_name(tmp_path):
    """One stray feature refused the whole grid, naming neither it nor its file. Now the
    grid is the one without it, and the warning says which feature and where."""

    from myflopy.grid_identity import grid_pin

    inside = (box(100, 100, 300, 300), 2500.0)
    with pytest.warns(UserWarning, match=r"1 of 2 refinement features in .*r\.gpkg .* skipped: 'refinement_1'"):
        vor, _ = _voronoi(tmp_path / "both", L_SHAPE, refinement=[inside, (box(2000, 2000, 2100, 2100), 2500.0)],
                          max_area=10000.0)
    alone, _ = _voronoi(tmp_path / "alone", L_SHAPE, refinement=[inside], max_area=10000.0)
    assert grid_pin(vor) == grid_pin(alone)


def test_a_point_feature_given_width_is_skipped_when_it_is_outside(tmp_path):
    """A point refinement widened by refinement_buffer: outside, it has no area inside."""

    import geopandas as gpd
    from shapely.geometry import Point

    import myflopy as mf

    gpd.GeoDataFrame({"n": [0]}, geometry=[L_SHAPE], crs=2927).to_file(tmp_path / "b.gpkg")
    gpd.GeoDataFrame({"area": [400.0, 400.0]}, geometry=[Point(200, 200), Point(3000, 3000)],
                     crs=2927).to_file(tmp_path / "p.gpkg")
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=tmp_path / "b.gpkg"), max_area=10000.0,
                               refinement=mf.GeoPackageSourceSpec(path=tmp_path / "p.gpkg", fields={"area": "area"}),
                               refinement_buffer=50.0)
    with pytest.warns(UserWarning, match="1 of 2 refinement features"):
        assert spec.resolve(workspace=tmp_path / "grid").ncpl > 0


@pytest.mark.parametrize(
    ("region", "why"),
    [(box(2000, 2000, 2100, 2100), "lies entirely outside the model boundary"),
     (box(1000, 100, 1100, 200), "touches the model boundary only along a line")],
    ids=["disjoint", "edge-touching"],
)
def test_a_source_with_no_feature_inside_raises_naming_the_file_and_the_cause(tmp_path, region, why):
    """Every feature missing still refuses (ledger 211's rule): an edge-touching polygon was
    told to "pass buffer>0", advice for a point or a line."""

    from myflopy.modflow.mf6.grid.triangle import NoAreaInside

    with pytest.raises(NoAreaInside, match=rf"No refinement feature in .*r\.gpkg .*{why}"):
        _voronoi(tmp_path, L_SHAPE, refinement=[(region, 2500.0)], max_area=10000.0)


def test_a_negative_buffer_that_erases_a_breakline_names_the_buffer_and_the_width(tmp_path):
    from myflopy.modflow.mf6.grid.triangle import NoAreaInside

    with pytest.raises(NoAreaInside, match=r"Breakline 'breakline_0' has no area left after buffer=-10: its part inside "
                                           r"the boundary is at most 20 wide"):
        _breaklines(tmp_path, [[(100.0, 200.0), (300.0, 200.0)]], breakline_negative_buffer_after_clipping=-10)


def test_a_boundary_part_narrower_than_the_simplification_is_kept_as_drawn(tmp_path):
    """The 2 ft profile simplification took a corner off a 1000 x 2 ft strip and meshed half
    its area, saying nothing. A 10 ft square, which it does not change, is untouched."""

    from myflopy.grid_identity import grid_pin

    with pytest.warns(UserWarning) as record:  # item 2's warning too: a 2 ft strip holds no 10000
        vor, _ = _voronoi(tmp_path / "strip", box(0, 0, 1000, 2), max_area=10000.0)
    assert any("the part at (0, 0)-(1000, 2), 2000 -> 1000" in str(w.message) and "kept as drawn" in str(w.message)
               for w in record)
    assert sum(Polygon(np.asarray(vor.verts)[list(r)]).area for r in vor.iverts) == pytest.approx(2000.0)
    with pytest.warns(UserWarning, match="larger than the domain") as record:  # item 2's, not this one
        square, _ = _voronoi(tmp_path / "square", box(0, 0, 10, 10), max_area=10000.0)
    assert not any("kept as drawn" in str(w.message) for w in record)
    assert grid_pin(square) == "399b78376886f595"


# --------------------------------------------------------------------------------------------
# 5. A point constraint on an edge, or within rounding of one: placed ON the edge.

STAR = Polygon([(661.8362256750904, 883.1029801823399), (511.91496377356594, 626.1164105850135),
                (397.21997094758257, 885.1137352756316), (386.40892476116454, 557.6707655713188),
                (109.1024319089716, 485.9970101374019), (370.24990640694506, 465.76037745899043),
                (347.44602898891804, 135.9209639441857), (475.2659797603329, 367.603618660851),
                (619.0565569568803, 109.47346703950711), (597.078140534828, 424.58165921244836)])
SQUARE_AT_STATE_PLANE = box(1.2e6, 6.5e5, 1.2e6 + 1000, 6.5e5 + 1000)


def _zero_length_faces(vor):
    import flopy

    from myflopy.grid_identity import zero_length_faces

    props = vor.get_disv_gridprops()
    with __import__("tempfile").TemporaryDirectory() as ws:
        sim = flopy.mf6.MFSimulation(sim_ws=ws)
        gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
        disv = flopy.mf6.ModflowGwfdisv(gwf, nlay=1, ncpl=props["ncpl"], nvert=props["nvert"], top=1.0, botm=[0.0],
                                        vertices=props["vertices"], cell2d=props["cell2d"])
        return zero_length_faces(disv, sim.simulation_data)


def test_a_point_a_fraction_of_an_ulp_inside_an_edge_builds_on_it(tmp_path):
    """The fuzz's point, 1.4e-14 ft inside a star's edge: Triangle packed 448 nodes within a
    foot and the repair refused two of them 5.7e-14 ft apart. Now it is a vertex of that edge."""

    point = (393.48212885521264, 771.9026644960962)
    with pytest.warns(UserWarning, match=r"1 point constraint\(s\) within .* were placed ON"):
        vor, tri = _voronoi(tmp_path, STAR, points=[point], max_area=1000.0, profile="clean", min_angle=29)
    assert any(tuple(v) == point for v in tri._polygons[0])
    assert (np.hypot(*(np.asarray(tri.verts) - point).T) < 1.0).sum() < 20


def test_a_point_a_millionth_of_a_foot_inside_an_edge_at_state_plane_is_writable(tmp_path):
    """0.01 ft is what FloPy writes there: a well 1e-6 ft inside built 303 zero-length faces."""

    point = (1.2e6 + 1e-6, 6.5e5 + 437.123)
    with pytest.warns(UserWarning, match="placed ON"):
        vor, _ = _voronoi(tmp_path, SQUARE_AT_STATE_PLANE, points=[point], max_area=2000.0)
    assert _zero_length_faces(vor) == []


def test_a_point_clear_of_the_edges_is_untouched_and_one_outside_is_named(tmp_path):
    """0.5 ft inside is 25 written units from the edge: Triangle keeps it as it was. Outside,
    it still gets no cell, and the warning now says so."""

    import warnings

    inside, outside = (1.2e6 + 0.5, 6.5e5 + 437.123), (1.2e6 - 50.0, 6.5e5 + 400.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _, tri = _voronoi(tmp_path / "in", SQUARE_AT_STATE_PLANE, points=[inside], max_area=2000.0)
    assert any(tuple(v) == inside for v in np.asarray(tri.verts))
    with pytest.warns(UserWarning, match=r"1 point constraint\(s\) lie outside the meshed domain"):
        _voronoi(tmp_path / "out", SQUARE_AT_STATE_PLANE, points=[outside], max_area=2000.0)


def test_triangle_runs_bounded_and_otherwise_as_flopy_ran_it(tmp_path):
    """FloPy held all of Triangle's output in memory with no timeout; a runaway run grew by
    0.3-0.5 GB a second. A bounded run writes the same mesh, and a run past a bound raises
    a CalledProcessError, which the smoothing loop's fallback catches."""

    import subprocess

    from flopy.utils.triangle import Triangle

    from myflopy.modflow.mf6.grid.triangle import TriangleFailed, TriangleGrid

    def grid(ws):
        tri = TriangleGrid(model_ws=str(ws), angle=28)
        tri.set_domain_polygon(L_SHAPE, max_area=5000.0)
        tri.prepare()
        return tri

    ours, flopys = grid(tmp_path / "ours"), grid(tmp_path / "flopys")
    (tmp_path / "ours").mkdir()
    (tmp_path / "flopys").mkdir()
    ours._bounded_build()
    Triangle.build(flopys)
    for suffix in ("1.node", "1.ele", "1.edge", "1.neigh"):
        assert (tmp_path / "ours" / f"_triangle.{suffix}").read_bytes() == (tmp_path / "flopys" / f"_triangle.{suffix}").read_bytes()
    capped = grid(tmp_path / "capped")
    (tmp_path / "capped").mkdir()
    capped.TRIANGLE_OUTPUT_LIMIT = 16
    with pytest.raises(TriangleFailed, match="wrote more than"):
        capped.build()  # the public route runs bounded too
    assert issubclass(TriangleFailed, subprocess.CalledProcessError)


# --------------------------------------------------------------------------------------------
# What a review of the five found, each reproduced before it was fixed.


@pytest.mark.skipif(__import__("os").name != "posix", reason="RLIMIT_AS is POSIX")
def test_triangle_runs_under_a_hard_memory_limit_below_its_own_cap(tmp_path):
    """The cap RAISED the hard limit, which an unprivileged process cannot: under
    `ulimit -v 4000000` every grid build failed with a bare SubprocessError."""

    import subprocess
    import sys
    import textwrap

    code = textwrap.dedent(f"""
        import resource
        limit = int(3.5 * 2**30)
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
        from shapely.geometry import box
        from myflopy.modflow.mf6.grid.triangle import TriangleGrid
        tri = TriangleGrid(model_ws={str(tmp_path)!r}, angle=28)
        tri.set_domain_polygon(box(0, 0, 100, 100), max_area=100.0)
        tri.build()
        print("built", tri.ncpl)
    """)
    run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert run.returncode == 0 and "built" in run.stdout, run.stderr[-2000:]


@pytest.mark.skipif(__import__("os").name != "posix", reason="a shell script and SIGALRM")
def test_an_interrupt_while_triangle_runs_kills_it(tmp_path, monkeypatch):
    """The finally block cancelled the timeout and left a hung Triangle running."""

    import os
    import signal
    import subprocess

    from myflopy.modflow.mf6.grid.triangle import TriangleGrid

    fake = tmp_path / "triangle.sh"
    fake.write_text("#!/bin/sh\nsleep 60\n")
    fake.chmod(0o755)
    spawned = []

    class Spy(subprocess.Popen):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            spawned.append(self)

    monkeypatch.setattr(subprocess, "Popen", Spy)

    def interrupt(signum, frame):
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGALRM, interrupt)
    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=28)
    tri.set_domain_polygon(box(0, 0, 100, 100), max_area=100.0)
    tri.exe_name = str(fake)
    try:
        signal.setitimer(signal.ITIMER_REAL, 0.5)
        with pytest.raises(KeyboardInterrupt):
            tri.build()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
    assert spawned and all(p.poll() is not None for p in spawned)
    assert not any(os.path.exists(f"/proc/{p.pid}") and p.returncode is None for p in spawned)


def test_a_point_exactly_on_an_edge_that_another_point_bends_goes_into_the_ring_too(tmp_path):
    """p0 exactly on the bottom edge was left to Triangle; q, 0.005 ft above the same edge,
    was spliced into it, and p0 then sat off the bent chord: outside it, with no cell."""

    p0, q = (1.2e6 + 400.0, 6.5e5 + 0.0), (1.2e6 + 500.0, 6.5e5 + 0.005)
    with pytest.warns(UserWarning, match="placed ON"):
        vor, tri = _voronoi(tmp_path, SQUARE_AT_STATE_PLANE, points=[p0, q], max_area=2000.0)
    ring = [tuple(v) for v in tri._polygons[0]]
    assert p0 in ring and q in ring and abs(ring.index(p0) - ring.index(q)) == 1  # in the ring's own direction
    assert _zero_length_faces(vor) == []


def test_a_point_near_a_contact_two_zones_share_goes_into_both_rings(tmp_path):
    """The contact is a segment of BOTH zones' rings; spliced into one, the other kept its
    straight segment 1e-6 ft from the new vertex and Triangle ran out of memory."""

    origin = np.array([1.2e6, 6.5e5])
    tiles = [box(*(origin + (100, 100)), *(origin + (200, 200))), box(*(origin + (200, 100)), *(origin + (300, 200)))]
    domain = Polygon([tuple(origin + c) for c in L_SHAPE.exterior.coords])
    well = tuple(origin + (200 + 1e-6, 150.123))
    with pytest.warns(UserWarning, match="placed ON"):
        vor, tri = _voronoi(tmp_path, domain, refinement=[(t, 2500.0) for t in tiles], points=[well],
                            max_area=10000.0, min_angle=28)
    assert sum(well in [tuple(v) for v in ring] for ring in tri._polygons) >= 2
    assert _zero_length_faces(vor) == []


def test_a_point_on_an_island_in_a_hole_is_inside(tmp_path):
    """A part inside another part's hole was called outside: a false 'no cell' warning,
    and a well a hair inside the island's edge was not placed (388 zero-length faces)."""

    import warnings

    origin = np.array([1.2e6, 6.5e5])
    donut = box(*origin, *(origin + 1000)).difference(box(*(origin + 300), *(origin + 700)))
    island = box(*(origin + (320, 320)), *(origin + (400, 400)))
    domain = donut.union(island)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        vor, _ = _voronoi(tmp_path / "mid", domain, points=[tuple(origin + (360, 361))], max_area=2000.0)
    assert not any("outside the meshed domain" in str(w.message) for w in record)
    with pytest.warns(UserWarning, match="placed ON"):
        vor, _ = _voronoi(tmp_path / "edge", domain, points=[tuple(origin + (320 + 1e-6, 361))], max_area=2000.0)
    assert _zero_length_faces(vor) == []


def test_an_inside_copy_is_never_merged_into_an_outside_copy(tmp_path):
    """One well recorded twice either side of the edge, outside first: the inside record
    was merged into the outside one, which gets no cell -- the well vanished."""

    inside = (1.2e6 + 0.004, 6.5e5 + 437.123)
    with pytest.warns(UserWarning):
        _, tri = _voronoi(tmp_path, SQUARE_AT_STATE_PLANE, points=[(1.2e6 - 0.004, 6.5e5 + 437.123), inside],
                          max_area=2000.0)
    assert inside in [tuple(v) for v in tri._polygons[0]]


def test_a_hole_seed_is_never_on_an_island_in_it(tmp_path):
    """The hole's representative point was the island's centre: Triangle ate the island and
    meshed the water (975,888 ft2 of cells on an 880,000 ft2 domain)."""

    origin = np.array([1.2e6, 6.5e5])
    donut = box(*origin, *(origin + 1000)).difference(box(*(origin + 300), *(origin + 700)))
    domain = donut.union(box(*(origin + 400), *(origin + 600)))
    vor, _ = _voronoi(tmp_path, domain, max_area=2000.0)
    meshed = sum(Polygon(np.asarray(vor.verts)[list(r)]).area for r in vor.iverts)
    assert meshed == pytest.approx(domain.area, rel=1e-9)


def test_a_line_with_no_buffer_is_refused_naming_the_buffer(tmp_path):
    """breakline_buffer=0 made a line an EMPTY region, which read as 'outside' (at nan
    ft) and was skipped with a warning; before, it was refused."""

    with pytest.raises(ValueError, match=r"breakline_buffer=0: a line needs a positive buffer"):
        _breaklines(tmp_path, [[(100.0, 200.0), (300.0, 200.0)]], breakline_buffer=0)


def test_a_buffered_point_just_outside_the_boundary_is_refused_not_skipped(tmp_path):
    """Its 50 ft disc is half inside, but the clip comes before the buffer: skipping it
    silently dropped a boundary well's refinement. Refused, as before, now saying why."""

    import geopandas as gpd
    from shapely.geometry import Point

    import myflopy as mf
    from myflopy.modflow.mf6.grid.triangle import NoAreaInside

    gpd.GeoDataFrame({"n": [0]}, geometry=[L_SHAPE], crs=2927).to_file(tmp_path / "b.gpkg")
    gpd.GeoDataFrame({"area": [400.0, 400.0]}, geometry=[Point(200, 200), Point(200, 1000 + 1e-6)],
                     crs=2927).to_file(tmp_path / "p.gpkg")
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=tmp_path / "b.gpkg"), max_area=10000.0,
                               refinement=mf.GeoPackageSourceSpec(path=tmp_path / "p.gpkg", fields={"area": "area"}),
                               refinement_buffer=50.0)
    with pytest.raises(NoAreaInside, match=r"within its buffer=50.*clipped to the boundary BEFORE it is buffered"):
        spec.resolve(workspace=tmp_path / "grid")


def test_the_face_check_reads_the_numbers_mf6_reads_and_skips_inactive_cells(tmp_path):
    """Across FloPy's 1e5 format switch '100000.00000000' and '1.00000000E+05' are one
    number to MF6; and a face beside an inactive cell makes no connection."""

    import flopy

    from myflopy.grid_identity import unwritable_faces

    rings = [[0, 7, 8, 9, 1], [1, 9, 3, 2], [8, 5, 4, 3, 9], [7, 6, 5, 8]]  # corner 8 split: face 8-9

    def quads(x0, gap, idomain=None):
        y0 = 5e4
        verts = np.array([(x0 - 100, y0 - 100), (x0, y0 - 100), (x0 + 100, y0 - 100), (x0 + 100, y0),
                          (x0 + 100, y0 + 100), (x0, y0 + 100), (x0 - 100, y0 + 100), (x0 - 100, y0),
                          (x0 - 3e-9, y0), (x0 + gap, y0)])
        sim = flopy.mf6.MFSimulation(sim_ws=str(tmp_path))
        gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
        centres = [(-50, -50), (50, -50), (50, 50), (-50, 50)]
        disv = flopy.mf6.ModflowGwfdisv(
            gwf, nlay=1, ncpl=4, nvert=len(verts), top=10.0, botm=[0.0], idomain=idomain,
            vertices=[[i, x, y] for i, (x, y) in enumerate(verts)],
            cell2d=[[k, x0 + cx, 5e4 + cy, len(r), *r] for k, ((cx, cy), r) in enumerate(zip(centres, rings, strict=True))])
        return unwritable_faces(disv, sim.simulation_data)

    assert [f[4] for f in quads(1e5, 3e-4)] == ["zero length"]
    assert quads(1e5, 3e-4, idomain=[[0, 1, 1, 1]]) == []


def test_a_face_whose_cell_centres_are_both_written_on_it_is_refused(tmp_path):
    """A 0.006 ft hairline notch at state plane, 'baseline' profile: no face is written with
    zero length, but three have BOTH cell centres written on their line, and MF6 stopped
    on a division by zero -- a model the zero-length check passed."""

    import myflopy as mf

    x, y = 1.2e6, 6.5e5
    notch = Polygon([(x, y), (x + 500, y), (x + 500, y + 400), (x + 500.006, y + 400), (x + 500.006, y),
                     (x + 1000, y), (x + 1000, y + 1000), (x, y + 1000)])
    vor, _ = _voronoi(tmp_path, notch, max_area=5000.0, profile="baseline")
    props = vor.get_disv_gridprops()
    model = mf.gwf("m", packages=[
        mf.disv(nlay=1, ncpl=props["ncpl"], nvert=props["nvert"], top=10.0, botm=[0.0],
                vertices=props["vertices"], cell2d=props["cell2d"]),
        mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds")])
    sim = mf.simulation(model, name="n", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with pytest.raises(ValueError, match="centres on it"):
        sim.build_flopy(workspace=tmp_path / "sim")


def test_small_things_the_review_found():
    """NoAreaInside survives a process pool; a jettison switch combined with others counts."""

    import pickle

    from myflopy.modflow.mf6.grid.triangle import NoAreaInside, TriangleGrid

    back = pickle.loads(pickle.dumps(NoAreaInside("gone", cause="outside")))
    assert (str(back), back.cause) == ("gone", "outside")
    tri = TriangleGrid(angle=28, additional_args=["-jz"])
    tri._polygons = [[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)], [(0.0, 0.0), (1.0, 1.0), (0.0, 1.0)]]
    assert tri._shared_nodes() is None
    tri.additional_args = ["-z"]
    assert tri._shared_nodes() is not None


def test_a_short_boundary_edge_flopy_drew_is_not_joined():
    """Edge 9-0 belongs to cell 0 alone: no connection, so MF6 runs it as written. Joining it
    re-pinned a grid that runs, against 'no grid that runs today moves'."""

    from myflopy.modflow.mf6.grid.tessellation import _onto_flopy

    x0, y0 = 1.2e6 + 100.0, 6.5e5 + 100.0
    verts = np.array([(x0 - 100, y0 - 100), (x0, y0 - 100), (x0 + 100, y0 - 100), (x0 + 100, y0),
                      (x0 + 100, y0 + 100), (x0, y0 + 100), (x0 - 100, y0 + 100), (x0 - 100, y0), (x0, y0),
                      (x0 - 100 + 7.4e-5, y0 - 100)])
    rings = [[0, 7, 8, 1, 9], [1, 8, 3, 2], [8, 5, 4, 3], [7, 6, 5, 8]]
    _, iverts, rewritten, appended, unsplit, _ = _onto_flopy(verts, rings, verts, rings, np.ones(4, dtype=bool))
    assert (rewritten, appended, unsplit) == ({}, 0, 0) and [list(r) for r in iverts] == rings


def test_a_package_on_flopys_own_cells_is_told_so_not_blamed_on_bb56668(tmp_path):
    """When the only join is a face FloPy drew itself, the grid bb56668 built IS FloPy's:
    its pin is legacy_pin, and the refusal must say the package predates the repair."""

    from types import SimpleNamespace

    from myflopy.grid_identity import GRID_PIN, _discretization_complaints

    corners = [(0, 0), (1100.3, 130.7), (1250.9, 940.2), (310.4, 1180.6)]
    quad = Polygon([(1.2e6 + x, 6.5e5 + y) for x, y in corners])
    vor, _ = _voronoi(tmp_path, quad, max_area=2000.0)
    report = vor.tessellation_report
    assert report.legacy_pin is not None and report.split_pin in (None, report.legacy_pin)
    package = SimpleNamespace(name="chd", metadata={GRID_PIN: report.legacy_pin})
    [complaint] = _discretization_complaints(None, SimpleNamespace(grid=vor, metadata={}), [package])
    assert "BEFORE ledger 316's repair" in complaint and "bb56668" not in complaint


def test_a_band_piece_resampled_to_nothing_is_named_even_when_the_rest_survives(tmp_path):
    """The notch cuts this band into a 901 ft piece and a 192 ft one; the short piece
    resampled to nothing without a word, because only a WHOLE region's loss was reported."""

    with pytest.warns(UserWarning, match=r"1 of its piece\(s\) were removed when resampled"):
        _breaklines(tmp_path, [[(100.0, 850.0), (600.0, 350.0)]])


def test_a_duplicate_row_counts_as_landed_and_names_both_files(tmp_path):
    """Across two sources: the warning named neither file ('breakline_0' is a duplicate of
    'breakline_0'), and the second source, its inside row deduplicated, was refused as if
    every feature missed."""

    import geopandas as gpd
    from shapely.geometry import LineString

    import myflopy as mf

    inside, outside = LineString([(100, 200), (300, 200)]), LineString([(2000, 2000), (2100, 2000)])
    gpd.GeoDataFrame({"n": [0]}, geometry=[L_SHAPE], crs=2927).to_file(tmp_path / "b.gpkg")
    gpd.GeoDataFrame({"n": [0]}, geometry=[inside], crs=2927).to_file(tmp_path / "one.gpkg")
    gpd.GeoDataFrame({"n": [0, 1]}, geometry=[inside, outside], crs=2927).to_file(tmp_path / "two.gpkg")
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=tmp_path / "b.gpkg"), max_area=3000.0,
                               breakline_max_area=60.0, breaklines=[mf.GeoPackageSourceSpec(path=tmp_path / "one.gpkg"),
                                                                    mf.GeoPackageSourceSpec(path=tmp_path / "two.gpkg")])
    with pytest.warns(UserWarning) as record:
        spec.resolve(workspace=tmp_path / "grid")
    said = " ".join(str(w.message) for w in record)
    assert "two.gpkg) is an exact duplicate of 'breakline_0' (feature 0 of" in said and "one.gpkg)" in said
    assert "1 of 2 breakline features in" in said
