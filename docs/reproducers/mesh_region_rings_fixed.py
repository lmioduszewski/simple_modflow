"""Why `protect_sources` / `protected_labels` cannot change a mesh, and what freeing the rings costs.

Ledger 411. `seed_optimization.get_fixed_and_free_vertex_ids` holds every vertex with
a nonzero Triangle boundary marker fixed, and FloPy's ``.poly`` writer gives EVERY segment
a marker (``iseg + 1``) -- refinement rings included, the domain's background sizing ring
too -- so every ring vertex is already fixed and protecting one adds nothing.

This script builds one 1000 x 1000 ft domain (background ``max_area=5000``) with a river
band (``add_line_feature``), a lake circle and a wellfield polygon, then runs
``build_mesh`` under each optimizing profile two ways:

* ``as shipped`` -- the rule above;
* ``rings free`` -- only vertices ON the domain boundary are fixed (plus protected ones),
  so an unprotected refinement ring may relax.

and with ``protect_sources=("line",)`` (the default) and ``None``.

Measured 2026-10-11 (triangle 1.6, FloPy 3.10):

* as shipped, protection changes nothing: ``balanced`` 842 free vertices and the same
  2,346 triangles either way (``fast`` 653 / 2,011, ``smooth`` 2,041 / 5,380).
* rings free, ``balanced`` frees 184 more vertices unprotected (1,026) and 47 with the
  river protected (889). The FIRST round is then REJECTED under every profile --
  ``tiny_triangle_increase`` and ``neighbor_area_ratio_regression`` -- because a Steiner
  point moved off its segment lands beside it and Triangle splits the segment again. So
  the optimizer falls back to the UNOPTIMIZED mesh: ``balanced`` 2,124 triangles (mean
  neighbour-area ratio 1.316 against 1.267 optimized), ``smooth`` 1,809 against 5,380.
  Protection changes the free count but not the mesh, since both fall back.

So a naive "free the unprotected rings" moves every default optimized mesh (``build_mesh``
defaults to ``"balanced"``, and so does ``GridSpec.voronoi``) and makes it worse. A ring
vertex would have to slide ALONG its segment for protection to mean something.

Run: ``PYTHONPATH=src python docs/reproducers/mesh_region_rings_fixed.py``
"""

from __future__ import annotations

import tempfile
import warnings

from shapely.geometry import LineString, Point, Polygon

import myflopy.modflow.mf6.grid.seed_optimization as seed_optimization
import myflopy.modflow.mf6.grid.triangle as triangle
from myflopy.modflow.mf6.grid.triangle import TriangleGrid

AS_SHIPPED = seed_optimization.get_fixed_and_free_vertex_ids


def _rings_free(domain):
    def fixed_and_free(*, node, protected_geometries=None, tolerance=1e-6):
        edge = domain.boundary.buffer(tolerance)
        fixed = {
            int(row["ivert"])
            for row in node
            if int(row["boundary_marker"]) != 0 and edge.covers(Point(float(row["x"]), float(row["y"])))
        }
        fixed.update(seed_optimization._protected_vertex_ids(node, protected_geometries, tolerance))
        return fixed, sorted(int(row["ivert"]) for row in node if int(row["ivert"]) not in fixed)

    return fixed_and_free


def build(profile, rule, protect_sources, workspace):
    tri = TriangleGrid(model_ws=workspace, angle=32)
    tri.set_domain_rectangle(1000, 1000, max_area=5000)
    tri.add_line_feature(LineString([(100, 100), (900, 800)]), buffer=15, max_area=200, label="river")
    tri.add_region_circle(center_coords=(300, 700), radius=80, max_area=300, label="lake")
    tri.add_region_polygon(Polygon([(600, 100), (900, 100), (900, 400), (600, 400)]), max_area=800, label="wf")
    chosen = AS_SHIPPED if rule == "as shipped" else _rings_free(tri.domain_geometry)
    # Both modules hold the name: the relaxation step reads seed_optimization's, the
    # free-point list triangle's.
    for module in (seed_optimization, triangle):
        module.get_fixed_and_free_vertex_ids = chosen
    try:
        report = tri.build_mesh(profile=profile, protect_sources=protect_sources)
    finally:
        for module in (seed_optimization, triangle):
            module.get_fixed_and_free_vertex_ids = AS_SHIPPED
    opt, quality = report["optimization"], report["quality"]
    return (
        f"{profile:8s} {rule:10s} protect={protect_sources!s:9s} free={opt['free_vertex_count']:5d} "
        f"status={opt['status']:22s} triangles={len(tri.ele):5d} "
        f"neighbour_area_ratio={quality['neighbor_area_ratio_mean']:.3f} "
        f"rejected={opt['rejection_reasons'] or '-'}"
    )


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    with tempfile.TemporaryDirectory() as root:
        for profile in ("balanced", "fast", "smooth"):
            for rule in ("as shipped", "rings free"):
                for protect_sources in (("line",), None):
                    print(build(profile, rule, protect_sources, f"{root}/{profile}_{rule[0]}_{protect_sources is None}"))
