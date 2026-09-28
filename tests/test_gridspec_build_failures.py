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

    with pytest.warns(UserWarning, match="removed when its outline was resampled") as record:
        vor, _ = _breaklines(tmp_path / "line", [[(100.0, 200.0), (140.0, 200.0)]])
    assert not any("empty after clipping" in str(w.message) for w in record)
    plain, _ = _voronoi(tmp_path / "plain", L_SHAPE, profile="clean", min_angle=23, max_area=3000.0)
    assert grid_pin(vor) == grid_pin(plain)
