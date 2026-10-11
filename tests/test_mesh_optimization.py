from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Polygon

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from myflopy.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid  # noqa: E402
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from myflopy.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from myflopy.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.packages import (  # noqa: E402
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
    for start, end in zip(coords[:-1], coords[1:], strict=False):
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


def test_build_mesh_profile_api_runs_optimized_mode():
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
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_optimized_triangle_grid_builds_voronoi_safe_mesh_and_runs_mf6():
    workspace = _project_temp_dir("mesh_optimization_run")
    try:
        tri, report = _build_complex_triangle_grid(workspace)
        quality = report["quality"]
        optimization = report["optimization"]

        assert optimization["status"] in {
            "applied", "fallback_original_mesh", "fallback_last_accepted_mesh",
        }
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


def test_add_polygon_buffers_point_and_line_features(tmp_path):
    # F1: point/line refinement features are polygonized by buffer (the old API
    # buffered them; the rewrite must too).
    from shapely.geometry import Point

    tri = TriangleGrid(model_ws=str(tmp_path / "pt"), angle=30)
    tri.set_domain_rectangle(x_dist=1000, y_dist=1000, origin=(0, 0), max_area=40000)
    tri.add_polygon(Point(500, 500), buffer=60, max_area=400)
    tri.build(verbose=False)
    assert VoronoiGridPlus(tri, crs="EPSG:2927").ncpl > 0


def test_add_polygon_point_without_buffer_errors_clearly(tmp_path):
    import pytest
    from shapely.geometry import Point

    tri = TriangleGrid(model_ws=str(tmp_path / "pt2"), angle=30)
    tri.set_domain_rectangle(x_dist=1000, y_dist=1000, origin=(0, 0))
    with pytest.raises(ValueError, match="no area"):
        tri.add_polygon(Point(200, 200))


def test_region_touching_domain_boundary_builds(tmp_path):
    # F2: a refinement region whose boundary coincides with the domain boundary
    # used to abort Triangle ("topological inconsistency"). It is auto-clipped to
    # the domain interior and now builds.
    tri = TriangleGrid(model_ws=str(tmp_path / "edge"), angle=30)
    tri.set_domain_rectangle(x_dist=1000, y_dist=1000, origin=(0, 0), max_area=40000)
    # Region shares the domain's left, top, and bottom edges.
    tri.add_region_polygon(
        Polygon([(0, 0), (300, 0), (300, 1000), (0, 1000)]),
        max_area=2000,
        label="edge_region",
    )
    tri.build(verbose=False)
    assert VoronoiGridPlus(tri, crs="EPSG:2927").ncpl > 0


def test_max_optimization_points_is_gone_everywhere(tmp_path):
    # Ledger 410: it was recorded in the report and limited nothing -- the helper
    # that would have applied it was never called -- so it is removed, not wired.
    import inspect

    import pytest

    import myflopy as mf

    for method in (TriangleGrid.optimize_seeds, TriangleGrid.build_optimized):
        assert "max_optimization_points" not in inspect.signature(method).parameters
    with pytest.raises(TypeError, match="max_optimization_points"):
        mf.GridSpec.voronoi(
            boundary=mf.ShapeSource(tmp_path / "domain.gpkg"), max_optimization_points=10
        )

    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=30)
    tri.set_domain_rectangle(x_dist=300, y_dist=200, max_area=2000)
    report = tri.build_mesh(profile="fast")
    assert "max_optimization_points" not in report["optimization"]
    with pytest.raises(TypeError, match="max_optimization_points"):
        tri.build_mesh(profile="fast", max_optimization_points=10)


def test_protection_options_are_gone_everywhere(tmp_path):
    # Ledger 445 (the user's call, on ledger 411): `protect_sources` /
    # `protected_labels` are removed. Every name raises at the call, the
    # non-optimizing build_mesh modes included -- they used to ignore any override.
    import inspect

    import pytest

    import myflopy as mf

    for method in (TriangleGrid.optimize_seeds, TriangleGrid.build_optimized, TriangleGrid.build_mesh):
        parameters = inspect.signature(method).parameters
        assert "protect_sources" not in parameters
        assert "protected_labels" not in parameters

    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=30)
    tri.set_domain_rectangle(x_dist=300, y_dist=200, max_area=2000)
    for profile in ("fast", "clean", "baseline"):
        with pytest.raises(TypeError, match="protect_sources"):
            tri.build_mesh(profile=profile, protect_sources=("line",))
        with pytest.raises(TypeError, match="protected_labels"):
            tri.build_mesh(profile=profile, protected_labels=["river"])
    with pytest.raises(TypeError, match="target_segment_lenght"):
        tri.build_mesh(profile="clean", target_segment_lenght=50)
    with pytest.raises(TypeError, match="protected_labels"):
        tri.build_optimized(protected_labels=["river"])
    with pytest.raises(TypeError, match="protected_labels"):
        tri.optimize_seeds(protected_labels=["river"])

    report = tri.build_mesh(profile="fast")
    assert not any("protect" in key for key in report["optimization"])
    assert not any("protect" in key for key in tri.build_mesh(profile="clean"))

    boundary = mf.ShapeSource(tmp_path / "domain.gpkg")
    with pytest.raises(TypeError, match="protect_sources"):
        mf.GridSpec.voronoi(boundary=boundary, protect_sources=("line",))
    for where in ("options", "mesh_options"):
        for key in ("protect_sources", "protected_labels"):
            with pytest.raises(TypeError, match=rf"{where}=\{{'{key}': \.\.\.\}}, which was removed"):
                mf.GridSpec.voronoi(boundary=boundary, **{where: {key: ("line",)}})


@pytest.mark.parametrize("where", ["options", "mesh_options"])
@pytest.mark.parametrize("key", ["protect_sources", "protected_labels"])
def test_a_saved_gridspec_carrying_a_removed_mesh_option_is_refused_at_resolve(tmp_path, where, key):
    # Ledger 445: a spec saved while the options were routed still carries them.
    # Filtering them at resolve, as any other stray key, would drop what the caller
    # wrote in silence, so resolve refuses naming the key -- before reading any file.
    import myflopy as mf

    payload = mf.GridSpec.voronoi(
        name="saved", boundary=mf.ShapeSource(tmp_path / "missing.gpkg"), boundary_max_area=100.0
    ).to_dict()
    payload[where] = {**payload.get(where, {}), key: ["river"]}
    restored = mf.GridSpec.from_dict(payload)
    with pytest.raises(TypeError, match=rf"GridSpec 'saved' carries {where}=\{{'{key}'"):
        restored.resolve(workspace=tmp_path / "ws")


def test_optimize_seeds_holds_vertices_near_a_line_band_ring(tmp_path):
    # Ledger 445: what the removed `protect_sources=("line",)` default did is the
    # fixed rule. Marked vertices -- every ring's -- are fixed anyway; this holds the
    # UNMARKED ones within the tolerance of a line band's ring, which Triangle can pack
    # there and which, freed, fell the refined workflow domain back to its unoptimized
    # mesh (docs/reproducers/mesh_region_rings_fixed.py, Part 2).
    from shapely.geometry import Point

    from myflopy.modflow.mf6.grid.seed_optimization import get_fixed_and_free_vertex_ids

    tri = TriangleGrid(model_ws=str(tmp_path / "ws"), angle=30)
    tri.set_domain_rectangle(400, 300, max_area=2000)
    tri.add_line_feature(LineString([(50, 50), (350, 250)]), buffer=20, max_area=300, label="river")
    tri.prepare()
    ring = next(r for r in tri._prepared_regions if r.source == "line").geometry.exterior
    on_ring = ring.interpolate(0.37, normalized=True)
    near = (on_ring.x + 0.007, on_ring.y + 0.007)
    tolerance = tri._region_point_spacing(geometry=tri.domain_geometry)
    assert 0 < ring.distance(Point(near)) < tolerance
    tri.add_points([near])

    report = tri.optimize_seeds(iterations=0)

    _, unmarked = get_fixed_and_free_vertex_ids(node=tri.node)
    rows = {int(row["ivert"]): Point(float(row["x"]), float(row["y"])) for row in tri.node}
    held = {vertex for vertex in unmarked if ring.distance(rows[vertex]) <= tolerance}
    near_ids = [vertex for vertex, point in rows.items() if (point.x, point.y) == near]
    # Triangle packs a cluster of unmarked vertices into the hair-wide gap, as on the
    # refined workflow domain; the point itself is one of them.
    assert len(near_ids) == 1 and near_ids[0] in held and len(held) > 1
    assert report["free_vertex_count"] == len(unmarked) - len(held)


def _optimizer_statuses_assigned():
    """Every ``status`` string literal ``optimize_seeds`` can put in its report."""

    import ast
    import inspect
    import textwrap

    found = set()
    nodes = [
        node
        for method in (TriangleGrid.optimize_seeds, TriangleGrid._mark_fallback)
        for node in ast.walk(ast.parse(textwrap.dedent(inspect.getsource(method))))
    ]
    for node in nodes:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant):
            for target in node.targets:
                if (
                    isinstance(target, ast.Subscript)
                    and isinstance(target.slice, ast.Constant)
                    and target.slice.value == "status"
                ):
                    found.add(node.value.value)
        if isinstance(node, ast.Dict):
            for key, value in zip(node.keys, node.values, strict=True):
                if isinstance(key, ast.Constant) and key.value == "status" and isinstance(value, ast.Constant):
                    found.add(value.value)
    return found


def test_every_status_the_optimizer_can_report_is_one_it_documents():
    # Ledger 410: "converged_small_move" was assigned behind a test no run could
    # pass (a move under min_move is never made, so max_move >= min_move whenever it
    # ran). The report's status is the reason the loop stopped; these five are all
    # (ledger 418 split the fallback in two).
    documented = {
        "applied", "converged_no_move", "fallback_original_mesh",
        "fallback_last_accepted_mesh", "not_attempted",
    }
    assert _optimizer_statuses_assigned() == documented
    doc = TriangleGrid.optimize_seeds.__doc__
    assert all(f'``"{status}"``' in doc for status in documented)


def test_the_optimizer_reports_why_it_stopped(tmp_path):
    def grid(name):
        tri = TriangleGrid(model_ws=str(tmp_path / name), angle=30)
        tri.set_domain_rectangle(x_dist=300, y_dist=200, max_area=2000)
        tri.build()
        return tri

    report = grid("still").optimize_seeds(iterations=3, min_move=1.0e6)
    assert (report["status"], report["iterations_run"]) == ("converged_no_move", 0)

    report = grid("ran").optimize_seeds(iterations=1)
    assert (report["status"], report["iterations_run"]) == ("applied", 1)

    report = grid("none").optimize_seeds(iterations=0)
    assert (report["status"], report["iterations_run"]) == ("not_attempted", 0)


@pytest.mark.parametrize(
    "reject_round, status, accepted", [(1, "fallback_original_mesh", 0), (2, "fallback_last_accepted_mesh", 1)]
)
def test_a_rejected_round_names_the_mesh_it_restores(tmp_path, monkeypatch, reject_round, status, accepted):
    # Ledger 418: one status, "fallback_original_mesh", covered both cases, while
    # a rejection after an accepted round restores THAT round's mesh. The rejection
    # is forced here (no default domain measured rejects a later round) -- what is
    # under test is the loop's bookkeeping, not the quality rules.
    def grid(name):
        tri = TriangleGrid(model_ws=str(tmp_path / name), angle=30)
        tri.set_domain_rectangle(x_dist=300, y_dist=200, max_area=2000)
        tri.build()
        return tri

    def vertices(tri):
        return np.unique(np.round(np.column_stack([tri.node["x"], tri.node["y"]]), 6), axis=0)

    expected = grid("expected")
    if accepted:
        expected.optimize_seeds(iterations=accepted)
    calls = {"n": 0}

    def forced(before, after):
        calls["n"] += 1
        return ["forced"] if calls["n"] == reject_round else []

    monkeypatch.setattr(TriangleGrid, "_quality_regression_reasons", staticmethod(forced))
    tri = grid("rejected")
    report = tri.optimize_seeds(iterations=3)
    assert (report["status"], report["iterations_run"]) == (status, accepted)
    assert report["rejection_reasons"] == "forced"
    np.testing.assert_array_equal(vertices(tri), vertices(expected))
