"""Which mesh vertices the optimizer holds fixed, and what freeing more of them costs.

Ledgers 411 and 445. `seed_optimization.get_fixed_and_free_vertex_ids` holds every
vertex with a nonzero Triangle boundary marker fixed, and FloPy's ``.poly`` writer gives
EVERY segment a marker (``iseg + 1``) -- refinement rings included, the domain's
background sizing ring too -- so every ring vertex is fixed. `TriangleGrid.optimize_seeds`
also holds every vertex within ``region_point_tolerance`` of a ``"line"`` band's ring (the
bands ``add_line_feature`` adds). That second rule was the default of the removed
``build_mesh(protect_sources=("line",))``; it is fixed now, and the ``line hold`` switch
below turns it off by patching the rule, which is what ``protect_sources=None`` did.

Part 1 builds ledger 411's 1000 x 1000 ft domain (background ``max_area=5000``) with a
river band, a lake circle and a wellfield polygon, then runs ``build_mesh`` under each
optimizing profile with three marker rules:

* ``as shipped`` -- the rule above;
* ``rings free`` -- only vertices ON the domain boundary keep their marker fixed, so a
  refinement ring may relax;
* ``bg fixed`` -- as ``rings free``, but the domain's background sizing ring stays fixed
  too, so only the rings of refinement regions a caller added may relax;

each with the line hold on and off. A last pair compares ``as shipped`` and ``bg fixed``
under ``"balanced"`` on a domain holding the background size and the river band only.

Part 2 builds the refined workflow domain of ``tests/test_mf6_refactor_smoke.py`` and
``refined_end_to_end_preferred_api_workflow.ipynb`` (5000 x 4000 ft, two stream bands, two
lakes, four refinement boxes) as shipped, with the line hold on and off, under the
settings that test and notebook pass (``"balanced"``, ``target_segment_length=100``,
``optimization_iterations=1``) and under ``build_mesh()``'s plain defaults.

Measured 2026-10-11 (triangle 1.6, FloPy 3.10):

* Part 1, as shipped, the line hold changes nothing: ``balanced`` 842 free vertices and the
  same 2,346 triangles either way (``fast`` 653 / 2,011, ``smooth`` 2,041 / 5,380).
* rings free, ``balanced`` frees 184 more vertices without the hold (1,026) and 47 with it
  (889). The FIRST round is then REJECTED under every profile -- ``tiny_triangle_increase``
  and ``neighbor_area_ratio_regression`` -- because a Steiner point moved off its segment
  lands beside it and Triangle splits the segment again. So the optimizer falls back to
  the UNOPTIMIZED mesh: ``balanced`` 2,124 triangles (mean neighbour-area ratio 1.316
  against 1.267 optimized), ``smooth`` 1,809 against 5,380.
* bg fixed falls back the same way on this domain (its lake and wellfield are not line
  bands) -- ``balanced`` 812 free with the hold, 949 without, and under every profile the
  same mesh as rings free. With the background size and the river band only it is
  byte-identical to as shipped (603 free, 1,718 triangles, the same vertex hash).
* Part 2 is where the line hold acts (445). Under the test's settings Triangle packs
  66 unmarked vertices within 0.05 ft of the ``stream_branch`` band's ring near
  (1850, 1513.5), where it meets the ``refine_zone_1`` box; the automatic tolerance is
  0.5 ft on this domain. Held (as shipped): 2,292 free, the round is applied, 9,818
  triangles at a mean neighbour-area ratio of 1.403. Not held: 2,358 free, the round is
  REJECTED (``tiny_triangle_increase``) and the optimizer falls back to the unoptimized
  7,614 triangles at 1.449. Under the plain defaults the hold changes nothing (2,972 free,
  9,397 triangles either way): the cluster comes with ``target_segment_length=100``.

So no ring vertex may be freed (Part 1: every rule that frees one falls back on a domain
holding a refinement region that is not a line band), and the hold near line bands is not
inert (Part 2): removing it falls back a build that left ``protect_sources`` at its
default. A ring vertex would have to slide
ALONG its segment for freeing one to mean something.

Run: ``PYTHONPATH=src python docs/reproducers/mesh_region_rings_fixed.py``
"""

from __future__ import annotations

import hashlib
import tempfile
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point, Polygon, box

import myflopy.modflow.mf6.grid.seed_optimization as seed_optimization
import myflopy.modflow.mf6.grid.triangle as triangle
from myflopy.modflow.mf6.grid.triangle import TriangleGrid

AS_SHIPPED = seed_optimization.get_fixed_and_free_vertex_ids


def _markers_on(keep):
    """Fix a marked vertex only when it lies on the boundary of a geometry in ``keep``."""

    def fixed_and_free(*, node, protected_geometries=None, tolerance=1e-6):
        bands = [geometry.boundary.buffer(tolerance) for geometry in keep]
        fixed = {
            int(row["ivert"])
            for row in node
            if int(row["boundary_marker"]) != 0
            and any(band.covers(Point(float(row["x"]), float(row["y"]))) for band in bands)
        }
        fixed.update(seed_optimization._protected_vertex_ids(node, protected_geometries, tolerance))
        return fixed, sorted(int(row["ivert"]) for row in node if int(row["ivert"]) not in fixed)

    return fixed_and_free


def _without_line_hold(rule):
    """``rule`` with the line-band hold switched off (what ``protect_sources=None`` did)."""

    def fixed_and_free(*, node, protected_geometries=None, tolerance=1e-6):
        return rule(node=node, protected_geometries=None, tolerance=tolerance)

    return fixed_and_free


def _run(tri, rule, line_hold, **build_mesh):
    chosen = rule if line_hold else _without_line_hold(rule)
    # Both modules hold the name: the relaxation step reads seed_optimization's, the
    # free-point list triangle's.
    for module in (seed_optimization, triangle):
        module.get_fixed_and_free_vertex_ids = chosen
    try:
        report = tri.build_mesh(**build_mesh)
    finally:
        for module in (seed_optimization, triangle):
            module.get_fixed_and_free_vertex_ids = AS_SHIPPED
    opt, quality = report["optimization"], report["quality"]
    digest = hashlib.sha1(np.asarray(tri.verts, dtype=float).tobytes()).hexdigest()[:10]
    return (
        f"line_hold={line_hold!s:5s} free={opt['free_vertex_count']:5d} "
        f"status={opt['status']:22s} triangles={len(tri.ele):5d} "
        f"neighbour_area_ratio={quality['neighbor_area_ratio_mean']:.3f} "
        f"verts={digest} rejected={opt['rejection_reasons'] or '-'}"
    )


def part_1(profile, rule, line_hold, workspace, *, features=True):
    tri = TriangleGrid(model_ws=workspace, angle=32)
    tri.set_domain_rectangle(1000, 1000, max_area=5000)
    tri.add_line_feature(LineString([(100, 100), (900, 800)]), buffer=15, max_area=200, label="river")
    if features:
        tri.add_region_circle(center_coords=(300, 700), radius=80, max_area=300, label="lake")
        tri.add_region_polygon(Polygon([(600, 100), (900, 100), (900, 400), (600, 400)]), max_area=800, label="wf")
    chosen = {
        "as shipped": lambda: AS_SHIPPED,
        "rings free": lambda: _markers_on([tri.domain_geometry]),
        "bg fixed": lambda: _markers_on(
            [tri.domain_geometry] + [r.geometry for r in tri.region_specs if r.source == "domain"]
        ),
    }[rule]()
    return f"{profile:8s} {rule:10s} " + _run(tri, chosen, line_hold, profile=profile)


def part_2(line_hold, workspace, **build_mesh):
    """The refined workflow domain, built as `tests/test_mf6_refactor_smoke.py` builds it."""

    root = Path(workspace)
    root.mkdir(parents=True)

    def gpkg(name, geometry):
        path = root / f"{name}.gpkg"
        gpd.GeoDataFrame({"name": [name]}, geometry=[geometry], crs="EPSG:2927").to_file(path, driver="GPKG")
        return path

    tri = TriangleGrid(model_ws=str(root / "mesh"), angle=30)
    tri.set_domain_file(
        gpkg("domain", box(0.0, 0.0, 5000.0, 4000.0)),
        simplify_tolerance=5, densify_dist=125, max_area=20000, label="domain",
    )
    streams = {
        "stream_main": [(250.0, 3500.0), (1200.0, 3050.0), (2300.0, 2400.0), (3550.0, 1650.0), (4750.0, 850.0)],
        "stream_branch": [(650.0, 950.0), (1500.0, 1400.0), (2500.0, 1900.0), (3650.0, 2250.0), (4550.0, 2500.0)],
    }
    for label, coords in streams.items():
        tri.add_line_feature(
            gpkg(label, LineString(coords)),
            buffer=55, simplify_tolerance=5, densify_dist=75, max_area=3000, label=label, priority=5,
        )
    lakes = {
        "lake_west": Point(1650.0, 2550.0).buffer(290.0, quad_segs=24),
        "lake_east": Point(3350.0, 1500.0).buffer(250.0, quad_segs=24),
    }
    for label, geometry in lakes.items():
        tri.add_region_file(
            gpkg(label, geometry), simplify_tolerance=5, densify_dist=50, max_area=2500, label=label, priority=6
        )
    zones = [box(350.0, 2200.0, 2200.0, 3800.0), box(1850.0, 900.0, 3400.0, 2500.0), box(3100.0, 1300.0, 4700.0, 3400.0)]
    for index, geometry in enumerate(zones):
        tri.add_region_file(
            gpkg(f"refine_zone_{index}", geometry),
            simplify_tolerance=5, densify_dist=80, max_area=5500, label=f"refine_zone_{index}", priority=4,
        )
    tri.add_region_polygon(
        box(900.0, 900.0, 4100.0, 3200.0), densify_dist=100, max_area=10000, label="interior_refinement", priority=3
    )
    settings = "test   " if build_mesh else "default"
    return f"refined {settings} " + _run(tri, AS_SHIPPED, line_hold, **build_mesh)


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    with tempfile.TemporaryDirectory() as root:
        print("Part 1:")
        for profile in ("balanced", "fast", "smooth"):
            for rule in ("as shipped", "rings free", "bg fixed"):
                for line_hold in (True, False):
                    print(part_1(profile, rule, line_hold, f"{root}/{profile}_{rule[0]}_{line_hold}"))
        print("background size and river band only:")
        for rule in ("as shipped", "bg fixed"):
            print(part_1("balanced", rule, True, f"{root}/only_{rule[0]}", features=False))
        print("Part 2:")
        test_settings = {"target_segment_length": 100, "optimization_iterations": 1}
        for line_hold in (True, False):
            print(part_2(line_hold, f"{root}/refined_test_{line_hold}", **test_settings))
            print(part_2(line_hold, f"{root}/refined_default_{line_hold}"))
