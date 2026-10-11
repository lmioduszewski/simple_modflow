"""What sizes Triangle meshes to: the regions, the domain's background and a global cap.

Ledger 408: a global ``maximum_area`` used to REPLACE every regional limit (FloPy
passes ``-a<value>`` instead of the bare ``-a`` that reads the .poly file's areas).
Ledger 409: re-declaring the domain on one TriangleGrid replaces its background
sizing region instead of keeping the previous one beside it.
"""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Point, box

from myflopy.modflow.mf6.grid.triangle import TriangleGrid


def _triangles(tri):
    """(centroids, areas) of the triangles Triangle wrote."""

    verts = np.asarray(tri.verts, dtype=float)
    rings = np.asarray(tri.iverts, dtype=int)
    p = verts[rings]  # (ntri, 3, 2)
    areas = 0.5 * np.abs(
        (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
        - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    )
    return p.mean(axis=1), areas


def test_a_global_maximum_area_keeps_every_refinement_region(tmp_path):
    region = box(400, 400, 600, 600)
    tri = TriangleGrid(model_ws=str(tmp_path / "capped"), angle=32, maximum_area=1.0e4)
    tri.set_domain_rectangle(1000, 1000)
    tri.add_region_polygon(region, max_area=100)
    tri.build()
    centroids, areas = _triangles(tri)
    inside = np.array([region.contains(Point(c)) for c in centroids])

    # The region at max_area=100 needs at least 400 triangles to cover 40,000 ft2;
    # FloPy's switch gave the whole domain 174.
    assert inside.sum() >= 400
    assert areas[inside].max() <= 100.0 * (1 + 1e-9)
    assert areas.max() <= 1.0e4 * (1 + 1e-9)

    plain = TriangleGrid(model_ws=str(tmp_path / "plain"), angle=32)
    plain.set_domain_rectangle(1000, 1000)
    plain.add_region_polygon(region, max_area=100)
    plain.build()
    # The cap only adds triangles where the background had none smaller than it.
    assert len(tri.ele) >= len(plain.ele)


def test_a_global_maximum_area_keeps_the_domain_background_size(tmp_path):
    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=32, maximum_area=1.0e5)
    tri.set_domain_rectangle(1000, 1000, max_area=5000)
    tri.build()
    _, areas = _triangles(tri)
    assert areas.max() <= 5000.0 * (1 + 1e-9)


@pytest.mark.parametrize("cap", [500.0, 5000.0])
def test_the_smaller_of_the_cap_and_the_region_limit_wins(tmp_path, cap):
    region = box(300, 300, 700, 700)
    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=32, maximum_area=cap)
    tri.set_domain_rectangle(1000, 1000, max_area=20_000)
    tri.add_region_polygon(region, max_area=1000)
    tri.build()
    centroids, areas = _triangles(tri)
    inside = np.array([region.contains(Point(c)) for c in centroids])
    assert areas[inside].max() <= min(cap, 1000.0) * (1 + 1e-9)
    assert areas[~inside].max() <= min(cap, 20_000.0) * (1 + 1e-9)


def test_re_declaring_the_domain_replaces_its_background_sizing_region(tmp_path):
    tri = TriangleGrid(model_ws=str(tmp_path / "twice"), angle=32)
    tri.set_domain_rectangle(1000, 1000, max_area=5000)
    tri.add_region_polygon(box(100, 100, 200, 200), max_area=50, label="kept")
    tri.set_domain_rectangle(500, 500, max_area=100)
    assert [(r.label, r.source, r.max_area) for r in tri.region_specs] == [
        ("kept", "geometry", 50),
        ("domain_size", "domain", 100),
    ]
    tri.build()

    fresh = TriangleGrid(model_ws=str(tmp_path / "once"), angle=32)
    fresh.set_domain_rectangle(500, 500, max_area=100)
    fresh.add_region_polygon(box(100, 100, 200, 200), max_area=50, label="kept")
    fresh.build()
    # The stale 5000 ft2 region's ring sat 0.07 ft inside the new edge: 85,839
    # triangles (of 0.0014 ft2) against 4,880 for the domain declared once.
    assert len(tri.ele) == len(fresh.ele)
    assert np.array_equal(np.asarray(tri.verts), np.asarray(fresh.verts))


def test_re_declaring_the_domain_without_max_area_drops_the_old_background(tmp_path):
    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=32)
    tri.set_domain_rectangle(1000, 1000, max_area=5000)
    tri.set_domain_circle(300.0, (500.0, 500.0))
    assert tri.region_specs == []


def test_the_domain_tag_is_refused_from_a_caller(tmp_path):
    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=32)
    tri.set_domain_rectangle(1000, 1000)
    with pytest.raises(ValueError, match=r"set_domain_\*\(max_area"):
        tri.add_region_polygon(box(100, 100, 200, 200), max_area=50, source="domain")
