"""A boundary that is not ONE simple polygon: parts and holes (ledger 317).

FloPy's ``tri2vor`` reads ``_polygons[0]`` as the domain and ``_polygons[1..nholes]`` as
its holes, and FloPy's ``Triangle.add_polygon`` hands an interior ring to ``add_hole``,
which wants a POINT. So a boundary drawn as two features kept cells in one part only, and
one with a hole did not build at all. The truth each test holds the grid to is geometric:
the cells cover the boundary exactly, and no cell sits in a hole.
"""

from __future__ import annotations

import numpy as np
import pytest
import shapely
from shapely.geometry import MultiPolygon, Point, Polygon, box

L_SHAPE = Polygon([(0, 0), (1000, 0), (1000, 400), (400, 400), (400, 1000), (0, 1000)])


def _resolve(tmp_path, geometries, max_area=5000.0, **options):
    import geopandas as gpd

    import myflopy as mf

    path = tmp_path / "boundary.gpkg"
    gpd.GeoDataFrame({"name": [f"f{i}" for i in range(len(geometries))]}, geometry=geometries, crs=2927).to_file(path)
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=path), max_area=max_area, **options)
    return spec.resolve(workspace=tmp_path / "grid", return_triangle=True)


def _cover(vor, boundary):
    cells = [Polygon(np.asarray(vor.verts)[iv]) for iv in vor.iverts]
    union = shapely.union_all(cells)
    return (boundary.difference(union).area, union.difference(boundary).area,
            sum(c.area for c in cells) - union.area)


def test_a_boundary_of_two_features_gets_cells_in_both(tmp_path):
    parts = [box(0, 0, 1000, 1000), box(1500, 0, 2500, 1000)]
    vor, _ = _resolve(tmp_path, parts)
    for part in parts:
        assert sum(part.covers(Point(p)) for p in np.asarray(vor.points)) > 300
    gap, outside, overlap = _cover(vor, MultiPolygon(parts))
    assert gap < 1e-6 and outside < 1e-6 and abs(overlap) < 1e-6
    assert vor.tessellation_report.adopted  # FloPy's cells covered the first part only


def test_a_boundary_with_a_hole_builds_and_tiles(tmp_path):
    holed = Polygon(box(0, 0, 1000, 1000).exterior.coords, [box(300, 300, 600, 600).exterior.coords])
    vor, _ = _resolve(tmp_path, [holed])
    gap, outside, overlap = _cover(vor, holed)
    assert gap < 1e-6 and outside < 1e-6 and abs(overlap) < 1e-6
    assert not any(box(300, 300, 600, 600).buffer(-1e-6).contains(Point(p)) for p in np.asarray(vor.points))
    # the hole's four corners are reflex corners of the domain: ledger 316's repair
    assert "corner-order" in set(vor.tessellation_report.rewritten.values())


def test_an_island_in_an_l_and_two_parts_with_a_hole(tmp_path):
    (tmp_path / "island").mkdir()
    (tmp_path / "pair").mkdir()
    island = Polygon(L_SHAPE.exterior.coords, [Point(200, 200).buffer(60).exterior.coords])
    vor, tri = _resolve(tmp_path / "island", [island])
    # a round island is simplified by 2 ft before meshing (profile "balanced"); the cells
    # tile the boundary that was MESHED
    assert max(abs(x) for x in _cover(vor, tri.domain_spec.geometry)) < 1e-6
    pair = [box(0, 0, 1000, 1000).difference(box(400, 400, 600, 600)), box(1500, 0, 2500, 1000)]
    vor, _ = _resolve(tmp_path / "pair", pair)
    assert max(abs(x) for x in _cover(vor, MultiPolygon(pair))) < 1e-6


def test_an_annular_refinement_keeps_its_middle_meshed(tmp_path):
    """A region's own ring went to add_hole too, and Triangle carved its disc out of the mesh."""

    from myflopy.modflow.mf6.grid.triangle import TriangleGrid
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    tri = TriangleGrid(model_ws=str(tmp_path), angle=32)
    tri.set_domain_polygon(L_SHAPE, max_area=5000.0)
    tri.add_region_polygon(Point(200, 200).buffer(120).difference(Point(200, 200).buffer(60)), max_area=500.0)
    tri.build_mesh(profile="clean")
    P, T = np.asarray(tri.verts), np.asarray(tri.iverts)
    meshed = shapely.union_all([Polygon(P[t]) for t in T])
    assert abs(meshed.area - L_SHAPE.area) < 1e-6
    assert max(abs(x) for x in _cover(VoronoiGridPlus(tri), L_SHAPE)) < 1e-6


def test_a_refinement_ring_on_a_hole_is_not_added_twice(tmp_path):
    """With no clip inset, a refinement around a hole has the hole as its ring: FloPy refused the duplicate."""

    from myflopy.modflow.mf6.grid.triangle import TriangleGrid
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    holed = Polygon(box(0, 0, 1000, 1000).exterior.coords, [box(300, 300, 600, 600).exterior.coords])
    tri = TriangleGrid(model_ws=str(tmp_path), angle=32, domain_clip_tolerance=0.0)
    tri.set_domain_polygon(holed, max_area=5000.0)
    tri.add_region_polygon(box(250, 250, 650, 650), max_area=1000.0)
    tri.build_mesh(profile="clean")
    assert max(abs(x) for x in _cover(VoronoiGridPlus(tri), holed)) < 1e-6


@pytest.mark.parametrize(
    "geometry, where",
    [
        (Polygon(box(0, 0, 1000, 1000).exterior.coords, [[(0, 0), (300, 200), (200, 300)]]), "hole 0"),
        (MultiPolygon([box(0, 0, 500, 500), box(500, 500, 1000, 1000)]), "part 1"),
    ],
    ids=["hole touching the edge", "parts meeting at a corner"],
)
def test_a_boundary_that_touches_itself_is_refused_by_name(tmp_path, geometry, where):
    from myflopy.modflow.mf6.grid.triangle import TriangleGrid

    tri = TriangleGrid(model_ws=str(tmp_path), angle=32)
    with pytest.raises(ValueError, match=rf"touches itself at .*{where}"):
        tri.set_domain_polygon(geometry, max_area=5000.0)
