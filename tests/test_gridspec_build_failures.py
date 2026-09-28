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
