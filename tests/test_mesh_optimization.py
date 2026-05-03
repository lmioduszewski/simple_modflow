from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Polygon

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from simple_modflow.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid  # noqa: E402
from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from simple_modflow.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from simple_modflow.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization  # noqa: E402
from simple_modflow.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    InitialConditions,
    KFlow,
    OutputControl,
    Storage,
)


def _project_temp_dir(name: str) -> Path:
    root = ROOT / ".pytest-work" / name
    if root.exists():
        shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    return root


def _max_segment_length(polygon: Polygon) -> float:
    coords = list(polygon.exterior.coords)
    lengths = []
    for start, end in zip(coords[:-1], coords[1:]):
        dx = float(end[0]) - float(start[0])
        dy = float(end[1]) - float(start[1])
        lengths.append((dx * dx + dy * dy) ** 0.5)
    return max(lengths)


def _build_complex_triangle_grid(workspace: Path) -> tuple[TriangleGrid, dict[str, object]]:
    tri = TriangleGrid(model_ws=str(workspace / "triangle"))
    tri.set_domain_rectangle(x_dist=550, y_dist=350, origin=(0, 0))
    tri.add_region_polygon(
        Polygon([(35, 25), (515, 25), (515, 325), (35, 325)]),
        max_area=16000,
        label="background_refine",
    )
    tri.add_line_feature(
        LineString([(50, 290), (190, 240), (360, 170), (490, 80)]),
        buffer=35,
        max_area=2600,
        label="stream_refine",
        priority=3,
        simplify_tolerance=1,
    )
    tri.add_region_circle(
        center_coords=(190, 120),
        radius=50,
        max_area=1200,
        label="lake_refine",
        priority=4,
    )
    tri.add_region_polygon(
        Polygon([(340, 210), (490, 210), (490, 300), (340, 300)]),
        max_area=3500,
        label="wetland_refine",
        priority=2,
    )

    report = tri.build_optimized(
        cleanup=True,
        target_segment_length=75,
        optimization_iterations=2,
        damping=0.3,
        protected_labels=["stream_refine", "lake_refine"],
        max_optimization_points=50,
        verbose=False,
    )
    return tri, report


def test_clean_geometry_selectively_resamples_feature_regions():
    workspace = _project_temp_dir("mesh_optimization_cleanup")
    try:
        tri = TriangleGrid(model_ws=str(workspace / "triangle"))
        tri.set_domain_rectangle(x_dist=400, y_dist=240, origin=(0, 0))
        tri.add_region_polygon(
            Polygon([(30, 20), (370, 20), (370, 220), (30, 220)]),
            max_area=9000,
            label="background_refine",
        )
        tri.add_line_feature(
            LineString([(30, 200), (160, 150), (330, 60)]),
            buffer=30,
            max_area=1800,
            label="stream_refine",
            priority=3,
            simplify_tolerance=1,
        )

        domain_vertices_before = len(tri.domain_spec.geometry.exterior.coords)
        region_vertices_before = {
            region.label: len(max(tri._iter_polygons(region.geometry), key=lambda poly: poly.area).exterior.coords)
            for region in tri.region_specs
        }
        stream_before = max(
            tri._iter_polygons(next(region.geometry for region in tri.region_specs if region.label == "stream_refine")),
            key=lambda poly: poly.area,
        )

        cleanup_report = tri.clean_geometry(
            target_segment_length=50,
        )

        domain_vertices_after = len(tri.domain_spec.geometry.exterior.coords)
        region_vertices_after = {
            region.label: len(max(tri._iter_polygons(region.geometry), key=lambda poly: poly.area).exterior.coords)
            for region in tri.region_specs
        }
        stream_after = max(
            tri._iter_polygons(next(region.geometry for region in tri.region_specs if region.label == "stream_refine")),
            key=lambda poly: poly.area,
        )

        assert cleanup_report["resample_domain_boundary"] is False
        assert cleanup_report["resample_region_sources"] == "line"
        assert domain_vertices_after == domain_vertices_before
        assert region_vertices_after["background_refine"] == region_vertices_before["background_refine"]
        assert _max_segment_length(stream_after) < _max_segment_length(stream_before)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_sync_triangle_inputs_stores_point_constraints_and_cvt_points_as_nodes():
    workspace = _project_temp_dir("mesh_optimization_nodes")
    try:
        tri = TriangleGrid(model_ws=str(workspace / "triangle"))
        tri.set_domain_rectangle(x_dist=100, y_dist=80, origin=(0, 0))
        tri.add_region_polygon(
            Polygon([(10, 10), (90, 10), (90, 70), (10, 70)]),
            max_area=1000,
            label="background_refine",
        )
        tri.add_points([(15, 15), (25, 25)])
        tri._optimization_points = [(35, 35), (45, 45)]

        tri.prepare()

        assert tri._nodes is not None
        assert tri._nodes.shape == (4, 2)
        assert np.allclose(tri._nodes[0], [15.0, 15.0])
        assert np.allclose(tri._nodes[-1], [45.0, 45.0])
        assert all(len(polygon) > 2 for polygon in tri._polygons)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_build_mesh_profile_api_runs_cleanup_only_mode():
    workspace = _project_temp_dir("mesh_profile_clean")
    try:
        tri = TriangleGrid(model_ws=str(workspace / "triangle"))
        tri.set_domain_rectangle(x_dist=220, y_dist=160, origin=(0, 0))
        tri.add_region_rectangle(origin=(20, 20), x_dist=180, y_dist=120, max_area=5000, label="background")
        tri.add_line_feature(
            LineString([(20, 140), (90, 100), (200, 45)]),
            buffer=18,
            max_area=900,
            label="stream_refine",
            priority=3,
        )

        report = tri.build_mesh(profile="clean", verbose=False)

        assert report["profile"] == "clean"
        assert report["mode"] == "cleanup_only"
        assert report["optimization"] is None
        assert report["quality"]["voronoi_status"] == "built"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_build_mesh_profile_api_auto_protects_feature_sources():
    workspace = _project_temp_dir("mesh_profile_balanced")
    try:
        tri = TriangleGrid(model_ws=str(workspace / "triangle"))
        tri.set_domain_rectangle(x_dist=320, y_dist=240, origin=(0, 0))
        tri.add_region_rectangle(origin=(25, 25), x_dist=270, y_dist=190, max_area=7000, label="background")
        tri.add_line_feature(
            LineString([(20, 200), (120, 150), (280, 50)]),
            buffer=20,
            max_area=1200,
            label="stream_refine",
            priority=4,
        )
        tri.add_region_circle(
            center_coords=(110, 80),
            radius=30,
            max_area=800,
            label="lake_refine",
            priority=5,
        )

        report = tri.build_mesh(profile=MeshBuildProfile.from_name("fast"), verbose=False)

        assert report["profile"] == "fast"
        assert report["mode"] == "optimized"
        assert report["optimization"]["method"] == "constrained_cvt_lloyd"
        assert "stream_refine" in report["optimization"]["protected_labels"]
        assert "lake_refine" not in report["optimization"]["protected_labels"]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_optimized_triangle_grid_builds_voronoi_safe_mesh_and_runs_mf6():
    workspace = _project_temp_dir("mesh_optimization_run")
    try:
        tri, report = _build_complex_triangle_grid(workspace)
        quality = report["quality"]
        optimization = report["optimization"]

        assert optimization["status"] in {"applied", "fallback_original_mesh"}
        assert optimization["method"] == "constrained_cvt_lloyd"
        assert optimization["iterations_run"] <= optimization["iterations_requested"]
        if optimization["status"] == "applied":
            assert optimization["rejection_reasons"] == ""
        else:
            assert optimization["rejection_reasons"] != ""
        assert quality["duplicate_vertex_count"] == 0
        assert quality["zero_area_triangle_count"] == 0
        assert quality["voronoi_status"] == "built"
        assert quality["triangle_angle_min_overall"] > 1.0
        assert quality["num_vertices"] > 1000
        assert quality["neighbor_area_ratio_max"] > 1.0

        vor = VoronoiGridPlus(tri)
        top = 120.0 - (0.018 * np.asarray(vor.centroids_x)) + (0.01 * np.asarray(vor.centroids_y))
        bottom = top - 35.0
        vor.gdf_topbtm = gpd.GeoDataFrame(
            {
                "geometry": vor.gdf_vorPolys.geometry,
                0: top,
                1: bottom,
            },
            geometry="geometry",
            crs=vor.crs,
        )

        model = SimulationBase(name="meshopt", mf_folder_path=workspace / "run", vor=vor, nper=1)
        DisvGrid(vor=vor, model=model, top=top.tolist(), bottom=[bottom.tolist()], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=(top - 2.0).tolist())
        KFlow(model=model, k=[12.0] * vor.ncpl, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)

        west = Polygon([(0, 0), (40, 0), (40, 350), (0, 350)])
        east = Polygon([(510, 0), (550, 0), (550, 350), (510, 350)])
        west_cells = sorted(set(vor.get_vor_cells_as_series(west).iloc[0]))
        east_cells = sorted(set(vor.get_vor_cells_as_series(east).iloc[0]))
        CHD(
            model=model,
            stress_period_data={
                0: [[(0, cell), 121.0] for cell in west_cells] + [[(0, cell), 98.0] for cell in east_cells]
            },
        )

        success, _ = model.run_simulation()

        assert success is True
        assert vor.ncpl == quality["voronoi_cell_count"]
        assert (workspace / "run" / "meshopt" / "meshopt.hds").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_quality_regression_guard_rejects_obviously_worse_mesh():
    baseline = {
        "duplicate_vertex_count": 0,
        "zero_area_triangle_count": 0,
        "voronoi_status": "built",
        "tiny_triangle_count": 0,
        "sliver_triangle_count": 2,
        "triangle_angle_min_overall": 16.0,
        "triangle_edge_ratio_mean": 1.35,
        "neighbor_area_ratio_mean": 1.25,
    }
    candidate = {
        "duplicate_vertex_count": 0,
        "zero_area_triangle_count": 0,
        "voronoi_status": "built",
        "tiny_triangle_count": 25,
        "sliver_triangle_count": 450,
        "triangle_angle_min_overall": 0.3,
        "triangle_edge_ratio_mean": 1.8,
        "neighbor_area_ratio_mean": 1.9,
    }

    reasons = TriangleGrid._quality_regression_reasons(baseline, candidate)

    assert "tiny_triangle_increase" in reasons
    assert "sliver_triangle_increase" in reasons
    assert "min_angle_regression" in reasons
    assert "edge_ratio_regression" in reasons
    assert "neighbor_area_ratio_regression" in reasons
