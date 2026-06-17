from __future__ import annotations

import ast
import json
import os
import shutil
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest
from shapely.geometry import LineString, Point, Polygon, box

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import figs  # noqa: E402
import myflopy  # noqa: E402
from myflopy import (  # noqa: E402
    ModelContext,
    Project,
    SimpleModelConfig as PackageSimpleModelConfig,
    SimulationBase as PackageSimulationBase,
    TriangleGrid as PackageTriangleGrid,
    UZFBuilder,
    VoronoiGridPlus as PackageVoronoiGridPlus,
    read_gpkg,
    read_shp_gpkg,
    simple_model_spec as package_simple_model_spec,
)
from myflopy.modflow import geotiff_to_contours, get_iheads  # noqa: E402
from myflopy.modflow.mf6 import (  # noqa: E402
    DRN as PackageDRN,
    DisvGrid as PackageDisvGrid,
    GHB as PackageGHB,
    ModelRegion as PackageModelRegion,
    OutputControl as PackageOutputControl,
    RegionGroup as PackageRegionGroup,
    RegionRegistry as PackageRegionRegistry,
    SimulationBase as PackageMf6SimulationBase,
    SimpleModelConfig as PackageMf6SimpleModelConfig,
    TemporalDiscretization as PackageTemporalDiscretization,
    VoronoiGridPlus as PackageMf6Voronoi,
    simple_model_spec as mf6_simple_model_spec,
)
from myflopy.modflow.mf6.grid.connectivity import build_disu_connectivity  # noqa: E402
from myflopy.modflow.mf6.grid.helpers import densify_poly, signed_area  # noqa: E402
from myflopy.modflow.mf6.grid.plotting import GridSection  # noqa: E402
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as GridVoronoiGridPlus  # noqa: E402
from myflopy.modflow.mf6.cross_section_plotting import (  # noqa: E402
    ModelCrossSectionStyle,
    plot_model_cross_section,
)
from myflopy.modflow.mf6.grid.selection import (  # noqa: E402
    get_grid_edge_cells,
    get_vor_cells_as_series,
)
from myflopy.modflow.mf6.boundaries import Boundaries  # noqa: E402
from myflopy.modflow.mf6.chd import CHDFromVector  # noqa: E402
from myflopy.modflow.mf6.drn import DRN, DRNFromVector  # noqa: E402
from myflopy.modflow.mf6.ghb import GHB, GHBFromVector  # noqa: E402
from myflopy.modflow.mf6.kflow import KFromVector  # noqa: E402
from myflopy.modflow.gwt.gwt import GWT  # noqa: E402
from myflopy.modflow.mf6.lakes import (  # noqa: E402
    LAKBuilder,
    LakeTableBuilder,
)
from myflopy.modflow.mf6.mvr import MVRBuilder, Move, MoverConnection  # noqa: E402
from myflopy.modflow.mf6.mfsimbase import SimulationBase  # noqa: E402
from myflopy.modflow.mf6.recharge import RechargeFromShp, RCHFromVector  # noqa: E402
from myflopy.modflow.mf6.sfr import SFRBuilder  # noqa: E402
from myflopy.modflow.mf6.simplemodel import (  # noqa: E402
    SimpleModelConfig,
    simple_model_spec,
)
from myflopy.modflow.mf6.simulation import SimulationBase as SimulationPkgBase  # noqa: E402
from myflopy.modflow.mf6.simulation.base import SimulationBase as SimulationBaseDirect  # noqa: E402
from myflopy.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisuGrid,
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.regions import (  # noqa: E402
    ModelRegion,
    RegionGroup,
    RegionRegistry,
)
from myflopy.modflow.utils.datatypes.xsections import XSection  # noqa: E402
from myflopy.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)
from myflopy.modflow.mp3du.particles import ParticleTrackingInput  # noqa: E402
from myflopy.modflow.mf6.voronoiplus import TriangleGrid, VoronoiGridPlus  # noqa: E402


def _model_context(model, vor):
    return ModelContext(grid=vor, domain=np.asarray(model.gwf.modelgrid.idomain))


def _surface_cells(context):
    domain = np.asarray(context.domain)
    if domain.ndim == 1:
        domain = domain[np.newaxis, :]
    return [
        (int(np.flatnonzero(domain[:, cell] > 0)[0]), cell)
        for cell in range(domain.shape[1])
        if np.any(domain[:, cell] > 0)
    ]


def _uzf_finf_from_rch(rch_dict, surface_cells, nper):
    values = {}
    for period in range(nper):
        by_cell = {tuple(row[0]): row[1] for row in rch_dict.get(period, [])}
        values[period] = [by_cell.get(cellid, 0.0) for cellid in surface_cells]
    return values


def _attach_sfr(
    model,
    vor,
    stream_paths,
    *,
    inflows=None,
    widths=10.0,
    gradients=0.001,
    mannings=0.03,
    streambed_k=1.0,
    streambed_thickness=1.0,
    mover=False,
    region_name_prefix=None,
    combined_region_name=None,
    region_tags=None,
    overwrite_regions=False,
):
    builder = SFRBuilder(
        context=_model_context(model, vor),
        nper=model.nper,
        streams=stream_paths,
        inflow=inflows,
        width=widths,
        gradient=gradients,
        roughness=mannings,
        streambed_k=streambed_k,
        streambed_thickness=streambed_thickness,
        mover=mover,
    )
    builder.build().build(model.gwf)
    for stream_id, cells in builder.stream_cells.items():
        name = f"{region_name_prefix}_{stream_id}" if region_name_prefix else stream_id
        model.add_region_from_cells(
            name,
            [(0, cell) for cell in cells],
            category="boundary",
            package="sfr",
            tags=region_tags or ["sfr"],
            geometry=builder.stream_table.loc[stream_id].geometry,
            overwrite=overwrite_regions,
        )
    if combined_region_name:
        model.add_region_from_cells(
            combined_region_name,
            [(0, cell) for cells in builder.stream_cells.values() for cell in cells],
            category="boundary",
            package="sfr",
            tags=region_tags or ["sfr"],
            overwrite=overwrite_regions,
        )
    return builder


def _attach_lak(
    model,
    vor,
    lake_paths,
    *,
    starting_stage,
    lake_bottom,
    bed_leakance=1.0,
    status="ACTIVE",
    mover=False,
    region_name_prefix=None,
    combined_region_name=None,
    region_tags=None,
    overwrite_regions=False,
):
    builder = LAKBuilder(
        context=_model_context(model, vor),
        nper=model.nper,
        lakes=lake_paths,
        lake_id_field="name",
        starting_stage=starting_stage,
        lake_bottom=lake_bottom,
        bed_leakance=bed_leakance,
        connection_modes="rectangular",
        status=status,
        mover=mover,
    )
    builder.build().build(model.gwf)
    for lake_id, cells in builder.lake_cells.items():
        name = f"{region_name_prefix}_{lake_id}" if region_name_prefix else lake_id
        model.add_region_from_cells(
            name,
            [(0, cell) for cell in cells],
            category="boundary",
            package="lak",
            tags=region_tags or ["lak"],
            geometry=builder.lake_table.loc[lake_id].geometry,
            overwrite=overwrite_regions,
        )
    if combined_region_name:
        model.add_region_from_cells(
            combined_region_name,
            [(0, cell) for cells in builder.lake_cells.values() for cell in cells],
            category="boundary",
            package="lak",
            tags=region_tags or ["lak"],
            overwrite=overwrite_regions,
        )
    return builder


def _two_cell_grid() -> gpd.GeoDataFrame:
    polys = [
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        Polygon([(1, 0), (2, 0), (2, 1), (1, 1)]),
    ]
    return gpd.GeoDataFrame({"cell": [0, 1]}, geometry=polys, crs="EPSG:2927").set_index("cell")


def _two_cell_vor():
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ],
        dtype=float,
    )
    iverts = [[0, 1, 2, 3], [1, 4, 5, 2]]
    xcyc = np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float)
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _two_cell_vor_clockwise():
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ],
        dtype=float,
    )
    iverts = [[0, 3, 2, 1], [1, 2, 5, 4]]
    xcyc = np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float)
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _project_temp_dir(name: str) -> Path:
    root = ROOT / ".pytest-work" / f"{name}_{time.time_ns()}"
    root.mkdir(parents=True, exist_ok=True)
    return root


def _write_gpkg(path: Path, gdf: gpd.GeoDataFrame) -> Path:
    gdf.to_file(path, driver="GPKG")
    return path


def _write_features(path: Path, rows: list[dict], crs: str) -> Path:
    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs=crs)
    gdf.to_file(path, driver="GPKG")
    return path


def _two_cell_model(name: str, workspace: Path, *, nper: int = 2):
    vor = _two_cell_vor()
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: [10.0, 12.0],
            1: [2.0, 4.0],
        },
        geometry="geometry",
        crs=vor.crs,
    )
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=nper)
    DisvGrid(vor=vor, model=model, top=[10.0, 12.0], bottom=[[2.0, 4.0]], nlay=1, idomain=[[1, 1]])
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    KFlow(model=model, k=[1.0, 1.0], k33_vert=[0.25, 0.4])
    return model, vor


def _four_cell_vor_clockwise():
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ],
        dtype=float,
    )
    iverts = [
        [0, 3, 4, 1],
        [1, 4, 5, 2],
        [3, 6, 7, 4],
        [4, 7, 8, 5],
    ]
    xcyc = np.array(
        [
            [0.5, 0.5],
            [1.5, 0.5],
            [0.5, 1.5],
            [1.5, 1.5],
        ],
        dtype=float,
    )
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _write_refined_workflow_vectors(workspace: Path) -> dict[str, object]:
    """Create synthetic-but-realistic geopackage inputs for the end-to-end workflow test."""

    crs = "EPSG:2927"
    domain_geom = box(0.0, 0.0, 5000.0, 4000.0)
    west_band = box(0.0, 0.0, 225.0, 4000.0)
    east_band = box(4775.0, 0.0, 5000.0, 4000.0)
    north_band = box(0.0, 3650.0, 5000.0, 4000.0)

    stream_main = LineString(
        [
            (250.0, 3500.0),
            (1200.0, 3050.0),
            (2300.0, 2400.0),
            (3550.0, 1650.0),
            (4750.0, 850.0),
        ]
    )
    stream_branch = LineString(
        [
            (650.0, 950.0),
            (1500.0, 1400.0),
            (2500.0, 1900.0),
            (3650.0, 2250.0),
            (4550.0, 2500.0),
        ]
    )
    lake_west = Point(1650.0, 2550.0).buffer(290.0, quad_segs=24)
    lake_east = Point(3350.0, 1500.0).buffer(250.0, quad_segs=24)
    recharge_west = box(350.0, 2200.0, 2200.0, 3800.0)
    recharge_central = box(1850.0, 900.0, 3400.0, 2500.0)
    recharge_east = box(3100.0, 1300.0, 4700.0, 3400.0)
    k_west = box(0.0, 0.0, 2100.0, 4000.0)
    k_central = box(1650.0, 700.0, 3450.0, 3200.0)
    k_east = box(2900.0, 0.0, 5000.0, 4000.0)

    domain_path = _write_gpkg(
        workspace / "refined_domain.gpkg",
        gpd.GeoDataFrame({"name": ["domain"]}, geometry=[domain_geom], crs=crs),
    )
    stream_paths = [
        _write_gpkg(
            workspace / "stream_main.gpkg",
            gpd.GeoDataFrame({"name": ["stream_main"]}, geometry=[stream_main], crs=crs),
        ),
        _write_gpkg(
            workspace / "stream_branch.gpkg",
            gpd.GeoDataFrame({"name": ["stream_branch"]}, geometry=[stream_branch], crs=crs),
        ),
    ]
    lake_paths = [
        _write_gpkg(
            workspace / "lake_west.gpkg",
            gpd.GeoDataFrame({"name": ["lake_west"]}, geometry=[lake_west], crs=crs),
        ),
        _write_gpkg(
            workspace / "lake_east.gpkg",
            gpd.GeoDataFrame({"name": ["lake_east"]}, geometry=[lake_east], crs=crs),
        ),
    ]
    chd_path = _write_gpkg(
        workspace / "boundary_chd.gpkg",
        gpd.GeoDataFrame({"name": ["west_chd"], "elev": [0.0], "layer": [1]}, geometry=[west_band], crs=crs),
    )
    ghb_path = _write_gpkg(
        workspace / "boundary_ghb.gpkg",
        gpd.GeoDataFrame(
            {"name": ["east_ghb"], "elev": [0.0], "height": [0.0], "cond": [20.0], "layer": [1], "min_elev": [0.0]},
            geometry=[east_band],
            crs=crs,
        ),
    )
    drn_path = _write_gpkg(
        workspace / "boundary_drn.gpkg",
        gpd.GeoDataFrame(
            {"name": ["north_drn"], "height": [-1.0], "cond": [5.0], "layer": [1], "min_elev": [0.0]},
            geometry=[north_band],
            crs=crs,
        ),
    )
    rch_path = _write_gpkg(
        workspace / "recharge_zones.gpkg",
        gpd.GeoDataFrame(
            {
                "zone": ["west_uplands", "central_infiltration", "east_uplands"],
                "rch_0": [0.00005, 0.00008, 0.00004],
                "rch_1": [0.00008, 0.00012, 0.00006],
            },
            geometry=[recharge_west, recharge_central, recharge_east],
            crs=crs,
        ),
    )
    k_path = _write_gpkg(
        workspace / "k_zones.gpkg",
        gpd.GeoDataFrame(
            {
                "name": ["west_k", "central_k", "east_k"],
                "k": [22.0, 12.0, 16.0],
                "layer": [1, 1, 1],
            },
            geometry=[k_west, k_central, k_east],
            crs=crs,
        ),
    )
    refinement_paths = [
        _write_gpkg(
            workspace / "refine_recharge_west.gpkg",
            gpd.GeoDataFrame({"name": ["west_uplands"]}, geometry=[recharge_west], crs=crs),
        ),
        _write_gpkg(
            workspace / "refine_recharge_central.gpkg",
            gpd.GeoDataFrame({"name": ["central_infiltration"]}, geometry=[recharge_central], crs=crs),
        ),
        _write_gpkg(
            workspace / "refine_recharge_east.gpkg",
            gpd.GeoDataFrame({"name": ["east_uplands"]}, geometry=[recharge_east], crs=crs),
        ),
    ]

    return {
        "crs": crs,
        "domain": domain_path,
        "streams": stream_paths,
        "lakes": lake_paths,
        "chd": chd_path,
        "ghb": ghb_path,
        "drn": drn_path,
        "rch": rch_path,
        "k": k_path,
        "refinement_regions": refinement_paths,
    }


def _build_refined_vector_test_triangle(workspace: Path):
    """Build a refined optimized triangle/voronoi grid from synthetic vector inputs."""

    inputs = _write_refined_workflow_vectors(workspace)
    tri = TriangleGrid(model_ws=str(workspace / "refined_mesh"), angle=30)
    tri.set_domain_file(
        inputs["domain"],
        simplify_tolerance=5,
        densify_dist=125,
        max_area=20000,
        label="domain",
    )
    for stream_path, label in zip(inputs["streams"], ["stream_main", "stream_branch"], strict=True):
        tri.add_line_feature(
            stream_path,
            buffer=55,
            simplify_tolerance=5,
            densify_dist=75,
            max_area=3000,
            label=label,
            priority=5,
        )
    for lake_path, label in zip(inputs["lakes"], ["lake_west", "lake_east"], strict=True):
        tri.add_region_file(
            lake_path,
            simplify_tolerance=5,
            densify_dist=50,
            max_area=2500,
            label=label,
            priority=6,
        )
    for idx, refinement_path in enumerate(inputs["refinement_regions"]):
        tri.add_region_file(
            refinement_path,
            simplify_tolerance=5,
            densify_dist=80,
            max_area=5500,
            label=f"refine_zone_{idx}",
            priority=4,
        )
    tri.add_region_polygon(
        box(900.0, 900.0, 4100.0, 3200.0),
        densify_dist=100,
        max_area=10000,
        label="interior_refinement",
        priority=3,
    )

    report = tri.build_mesh(
        profile="balanced",
        protect_sources=("line",),
        protected_labels=["stream_main", "stream_branch", "lake_west", "lake_east"],
        target_segment_length=100,
        optimization_iterations=1,
        verbose=False,
    )
    vor = VoronoiGridPlus(tri, crs=inputs["crs"], name="refined_vector_test")
    return tri, vor, report, inputs


def _refined_vor_top_bottom(vor: VoronoiGridPlus) -> tuple[list[float], list[float]]:
    """Create smooth synthetic top and bottom surfaces for a refined workflow test."""

    centroids = vor.gdf_vorPolys.geometry.centroid
    top = 142.0 - (centroids.x.to_numpy() / 150.0) - (centroids.y.to_numpy() / 500.0)
    botm = top - 75.0 - 4.0 * np.sin(centroids.x.to_numpy() / 850.0)
    botm = np.minimum(botm, top - 5.0)
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: top,
            1: botm,
        },
        geometry="geometry",
        crs=vor.crs,
    )
    return top.tolist(), botm.tolist()


class DummyVor:
    def __init__(self, gdf: gpd.GeoDataFrame):
        self.gdf_vorPolys = gdf
        self.adjacent_cells_idx = [[1], [0]]
        self.idomain_path = None
        self.crs = "EPSG:2927"


def test_imports_and_custom_figs_are_available():
    assert figs is not None
    assert myflopy.__version__
    assert SimulationBase is not None
    assert TriangleGrid is not None
    assert DRN is not None
    assert DRNFromVector is not None
    assert GHB is not None
    assert GHBFromVector is not None
    assert RechargeFromShp is not None
    assert RCHFromVector is not None
    assert LakeTableBuilder is not None
    assert SFRBuilder is not None
    assert GWT is not None
    assert ParticleTrackingInput is not None
    assert PackageSimulationBase is SimulationBase
    assert PackageSimpleModelConfig is SimpleModelConfig
    assert PackageTriangleGrid is TriangleGrid
    assert PackageVoronoiGridPlus is VoronoiGridPlus
    assert GridVoronoiGridPlus is VoronoiGridPlus
    assert PackageMf6Voronoi is VoronoiGridPlus
    assert PackageMf6SimulationBase is SimulationBase
    assert PackageMf6SimpleModelConfig is SimpleModelConfig
    assert PackageModelRegion is ModelRegion
    assert PackageRegionGroup is RegionGroup
    assert PackageRegionRegistry is RegionRegistry
    assert SimulationPkgBase is SimulationBase
    assert SimulationBaseDirect is SimulationBase
    assert package_simple_model_spec is simple_model_spec
    assert mf6_simple_model_spec is simple_model_spec
    assert PackageDisvGrid is DisvGrid
    assert PackageTemporalDiscretization is TemporalDiscretization
    assert PackageOutputControl is OutputControl
    assert PackageDRN is DRN
    assert PackageGHB is GHB
    assert read_gpkg is not None
    assert read_shp_gpkg is not None
    assert get_iheads is not None
    assert geotiff_to_contours is not None
    assert myflopy.validate_surface_water_configuration is not None


def test_examples_layout_smoke():
    examples_dir = ROOT / "examples" / "mf6"
    sample_dir = examples_dir / "sample_model_output"
    notebooks_dir = examples_dir / "notebooks"

    assert examples_dir.exists()
    assert (examples_dir / "artifacts").exists()
    assert notebooks_dir.exists()
    assert sample_dir.exists()
    assert (sample_dir / "top_raster.tif").exists()
    assert (notebooks_dir / "feature_rich_model_workflow.ipynb").exists()
    assert (notebooks_dir / "refined_feature_rich_model_workflow.ipynb").exists()
    assert (notebooks_dir / "code_geometry_refined_model_workflow.ipynb").exists()
    assert (notebooks_dir / "refined_end_to_end_preferred_api_workflow.ipynb").exists()


def test_example_notebooks_are_valid_json():
    notebook_names = [
        "simple_model_workflow.ipynb",
        "triangle_voronoi_simplemodel_workflow.ipynb",
        "feature_rich_model_workflow.ipynb",
        "refined_feature_rich_model_workflow.ipynb",
        "code_geometry_refined_model_workflow.ipynb",
        "cumberland_mesh_optimization_workflow.ipynb",
        "cumberland_cvt_diagnostics_workflow.ipynb",
        "refined_end_to_end_preferred_api_workflow.ipynb",
    ]
    notebooks_dir = ROOT / "examples" / "mf6" / "notebooks"

    for notebook_name in notebook_names:
        notebook_path = notebooks_dir / notebook_name
        payload = json.loads(notebook_path.read_text(encoding="utf-8"))
        assert payload["nbformat"] == 4
        assert len(payload["cells"]) >= 3


def test_internal_modules_do_not_import_legacy_facades():
    source_root = ROOT / "src" / "myflopy"
    allowed_legacy_modules = {
        source_root / "modflow" / "mf6" / "mfsimbase.py",
        source_root / "modflow" / "mf6" / "voronoiplus.py",
    }
    violations: list[str] = []

    for path in source_root.rglob("*.py"):
        if path in allowed_legacy_modules or "archive" in path.parts:
            continue

        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            imported_modules: list[str] = []
            if isinstance(node, ast.ImportFrom) and node.module:
                imported_modules.append(node.module)
            elif isinstance(node, ast.Import):
                imported_modules.extend(alias.name for alias in node.names)

            for module_name in imported_modules:
                if module_name.endswith(("mfsimbase", "voronoiplus")):
                    rel_path = path.relative_to(ROOT)
                    violations.append(f"{rel_path}:{node.lineno}:{module_name}")

    assert not violations, "Legacy façade imports remain in internal modules:\n" + "\n".join(violations)


def test_simulation_base_minimal_smoke():
    workspace = _project_temp_dir("simulation_base_smoke")
    try:
        model = SimulationBase(name="smoke_model", mf_folder_path=workspace, nper=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        OutputControl(model=model)

        assert model.sim.name == "smoke_model"
        assert model.gwf.name == "smoke_model"
        assert model.model_output_folder_path == workspace / "smoke_model"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_grid_helper_functions_smoke():
    coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    assert signed_area(coords) == 1.0

    poly = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    dense = densify_poly(poly, 1)
    assert len(dense.exterior.coords) > len(poly.exterior.coords)


def test_build_disu_connectivity_for_two_adjacent_cells():
    gdf = _two_cell_grid()
    iac, ja, cl12, hwva, nja = build_disu_connectivity(
        gdf,
        [[1], [0]],
        validate=True,
    )

    assert iac.tolist() == [2, 2]
    assert ja.tolist() == [0, 1, 1, 0]
    assert nja == 4
    assert np.allclose(hwva[[1, 3]], [1.0, 1.0])
    assert np.allclose(cl12[[1, 3]], [1.0, 1.0])


def test_selection_helpers_find_cells_and_edges():
    gdf = _two_cell_grid()
    vor = DummyVor(gdf)

    hit = Polygon([(0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75)])
    cells = get_vor_cells_as_series(gdf, hit)

    assert isinstance(cells, pd.Series)
    assert cells.iloc[0] == [0]

    edges = get_grid_edge_cells(vor)
    assert edges == [0, 1]


def test_direct_voronoi_grid_smoke():
    vor = _two_cell_vor()

    assert vor.ncpl == 2
    assert len(vor.get_voronoi_polygons()) == 2
    assert vor.adjacent_cells_idx == [[1], [0]]
    assert vor.get_grid_edge() == [0, 1]
    assert vor.get_cell_areas() == [1.0, 1.0]
    assert vor.get_origin_xy() == (0.0, 0.0)
    assert np.isclose(vor.get_domain().area, 2.0)
    assert vor.generate_grid_coordinates(1) == ([0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1])
    assert np.allclose(vor.centroids[0], [0.5, 1.5])
    assert np.allclose(vor.centroids[1], [0.5, 0.5])
    assert np.isclose(vor.calculate_distance(vor.gdf_vorPolys, 0, 1), 0.5)
    assert np.isclose(
        vor.shared_face_length(vor.gdf_vorPolys.geometry.iloc[0], vor.gdf_vorPolys.geometry.iloc[1]),
        1.0,
    )
    assert np.isclose(vor.get_overlapping_area(cell_list=[0]), 1.0)

    topbtm = vor.get_gdf_topbtm_multilyr([10.0, 0.0])
    assert topbtm[0].tolist() == [10.0, 10.0]
    assert topbtm[1].tolist() == [0.0, 0.0]
    vor.gdf_topbtm = topbtm

    query = gpd.GeoDataFrame(
        {"name": ["left", "both"]},
        geometry=[Point(0.5, 0.5), Polygon([(0.25, 0.25), (1.75, 0.25), (1.75, 0.75), (0.25, 0.75)])],
        crs="EPSG:2927",
    )
    idx_map = vor.get_vor_idx_from_geometry(gdf_to_query=query, name_col="name")
    assert idx_map["left"] == [0]
    assert idx_map["both"] == [0, 1]
    assert vor.get_vor_idx_from_geometry_idx(query, idx=1) == [0, 1]

    assert vor.set_k_vor({5.0: [1]}, k_default=1.0) == [1.0, 5.0]

    fig2d = vor.plot2d()
    fig3d = vor.plot3d()
    assert isinstance(fig2d, go.Figure)
    assert isinstance(fig3d, go.Figure)
    assert len(fig2d.data) == 2
    assert len(fig3d.data) == 2

    section = vor.cross_section(LineString([(0.5, -0.5), (0.5, 1.5)]))
    assert isinstance(section, GridSection)
    assert len(section.poly_coords) >= 1
    assert len(section.figure.data) >= 1
    section_df = section.to_frame()
    assert {"distance", "elevation", "series", "polygon_id"}.issubset(section_df.columns)
    mpl_fig, mpl_ax = section.plot_mpl()
    assert mpl_fig is not None
    assert len(mpl_ax.lines) >= 1

    iac, ja, cl12, hwva, nja = vor.get_disu_connectivity(validate=True)
    assert iac.tolist() == [2, 2]
    assert ja.tolist() == [0, 1, 1, 0]
    assert nja == 4
    assert np.allclose(cl12[[1, 3]], [1.0, 1.0])
    assert np.allclose(hwva[[1, 3]], [1.0, 1.0])


def test_xsection_mpl_adapter_smoke():
    class DummyXSection(XSection):
        def __init__(self):
            self.interpolate = False
            self.section_name = "dummy"
            self._layer = [0]
            self.show_model_top = False
            self.show_model_btm = False

        @property
        def xsect(self):
            return [0.0, 10.0], [[100.0, 95.0]]

    section = DummyXSection()
    section_df = section.to_frame()

    assert {"distance", "elevation", "series", "kind", "layer"}.issubset(section_df.columns)

    mpl_fig, mpl_ax = section.plot_mpl(show_legend=False)
    assert mpl_fig is not None
    assert len(mpl_ax.lines) == 1


def test_model_cross_section_plotting_smoke():
    workspace = _project_temp_dir("cross_section_plotting")
    try:
        model, _ = _two_cell_model("cross_section_plotting", workspace, nper=1)
        head_data = np.array([[8.0, 9.5]])
        style = ModelCrossSectionStyle(
            title="Smoke Cross Section",
            use_figs_theme=True,
        )
        fig, ax = plot_model_cross_section(
            model,
            LineString([(0.5, -0.5), (0.5, 1.5)]),
            head_data=head_data,
            head_layer=0,
            style=style,
            ylim=(0.0, 12.0),
            layer_colors=["#fff6cc"],
            layer_labels=["Layer 1"],
        )

        assert fig is not None
        assert ax.get_title() == "Smoke Cross Section"
        assert len(ax.lines) >= 1
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_minimal_model_can_write_inputs():
    workspace = _project_temp_dir("myflopy_write")
    try:
        vor = _two_cell_vor()
        model = SimulationBase(name="write_smoke", mf_folder_path=workspace, vor=vor, nper=1)

        top = [10.0, 10.0]
        botm = [[0.0, 0.0]]
        DisvGrid(vor=vor, model=model, top=top, bottom=botm, nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[5.0, 5.0])
        KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=model)
        OutputControl(model=model)

        model.sim.write_simulation()
        outdir = workspace / "write_smoke"

        assert outdir.exists()
        expected = {
            "mfsim.nam",
            "write_smoke.disv",
            "write_smoke.ic",
            "write_smoke.ims",
            "write_smoke.nam",
            "write_smoke.npf",
            "write_smoke.oc",
            "write_smoke.sto",
            "write_smoke.tdis",
        }
        for _ in range(10):
            actual = {path.name for path in outdir.iterdir()}
            if expected.issubset(actual):
                break
            time.sleep(0.1)
        assert expected.issubset(actual)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_minimal_model_can_run_and_read_heads():
    workspace = _project_temp_dir("myflopy_run")
    try:
        vor = _two_cell_vor_clockwise()
        model = SimulationBase(name="run_smoke", mf_folder_path=workspace, vor=vor, nper=1)

        top = [10.0, 10.0]
        botm = [[0.0, 0.0]]
        DisvGrid(vor=vor, model=model, top=top, bottom=botm, nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [
                    [(0, 0), 10.0],
                    [(0, 1), 9.0],
                ]
            },
        )

        success, _ = model.run_simulation()
        assert success is True

        outdir = workspace / "run_smoke"
        assert (outdir / "run_smoke.hds").exists()
        assert (outdir / "run_smoke.model").exists()

        heads = model.hds.get_data(kstpkper=(0, 0)).squeeze()
        assert np.allclose(heads, [10.0, 9.0])

        all_heads = model.all_heads
        heads_long = model.hds.long()
        heads_wide = model.hds.wide()
        assert len(all_heads) == 2
        assert np.allclose(all_heads["elev"].astype(float).to_numpy(), [10.0, 9.0])
        assert heads_long.index.names == ["kstpkper", "layer", "cell"]
        assert heads_long.name == "elev"
        assert np.allclose(heads_long.astype(float).to_numpy(), [10.0, 9.0])
        assert {"layer", "cell", "kstpkper_0_0"}.issubset(heads_wide.columns)
        assert np.allclose(heads_wide.sort_values("cell")["kstpkper_0_0"].astype(float).to_numpy(), [10.0, 9.0])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_region_registry_supports_cells_geometry_and_head_queries():
    workspace = _project_temp_dir("myflopy_region_registry")
    try:
        vor = _two_cell_vor_clockwise()
        model = SimulationBase(name="region_smoke", mf_folder_path=workspace, vor=vor, nper=1)

        top = [10.0, 10.0]
        botm = [[0.0, 0.0]]
        DisvGrid(vor=vor, model=model, top=top, bottom=botm, nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [
                    [(0, 0), 10.0],
                    [(0, 1), 9.0],
                ]
            },
        )

        success, _ = model.run_simulation()
        assert success is True

        left_region = model.add_region_from_cells(
            "left_cell",
            cellids=[(0, 0)],
            category="boundary",
            package="chd",
            tags=["left", "boundary"],
            metadata={"description": "left boundary cell"},
        )
        geometry_region = model.add_region_from_geometry(
            "right_geom",
            Polygon([(1.1, 0.1), (1.9, 0.1), (1.9, 0.9), (1.1, 0.9)]),
            category="selection",
            tags=["right"],
        )

        assert isinstance(model.regions, RegionRegistry)
        assert isinstance(left_region, ModelRegion)
        assert left_region.cells == [0]
        assert left_region.layer == 0
        assert geometry_region.cells == [1]
        assert model.get_region("left_cell").metadata["description"] == "left boundary cell"
        assert model.get_region_cells("right_geom") == [1]

        summary = model.list_regions()
        assert set(summary["name"]) == {"left_cell", "right_geom"}
        assert summary.set_index("name").loc["left_cell", "package"] == "chd"

        left_heads = model.region_heads("left_cell", per=0)
        right_heads = model.region_heads("right_geom", per=0)

        assert isinstance(left_heads, gpd.GeoDataFrame)
        assert left_heads["region"].unique().tolist() == ["left_cell"]
        assert right_heads["region"].unique().tolist() == ["right_geom"]
        assert left_heads["cell"].tolist() == [0]
        assert right_heads["cell"].tolist() == [1]
        assert np.allclose(left_heads["elev"].astype(float).to_numpy(), [10.0])
        assert np.allclose(right_heads["elev"].astype(float).to_numpy(), [9.0])
        assert left_heads.geometry.iloc[0].equals(vor.gdf_vorPolys.geometry.loc[0])

        removed = model.remove_region("right_geom")
        assert removed.name == "right_geom"
        assert model.list_regions()["name"].tolist() == ["left_cell"]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_region_groups_resolve_nested_membership_and_deduplicate_cells():
    workspace = _project_temp_dir("myflopy_region_groups")
    try:
        vor = _two_cell_vor_clockwise()
        model = SimulationBase(name="group_smoke", mf_folder_path=workspace, vor=vor, nper=1)

        top = [10.0, 10.0]
        botm = [[0.0, 0.0]]
        DisvGrid(vor=vor, model=model, top=top, bottom=botm, nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [
                    [(0, 0), 10.0],
                    [(0, 1), 9.0],
                ]
            },
        )

        success, _ = model.run_simulation()
        assert success is True

        region_a = model.add_region_from_cells("region_a", cellids=[(0, 0)], tags=["a"])
        region_b = model.add_region_from_cells("region_b", cellids=[(0, 0), (0, 1)], tags=["b"])
        region_c = model.add_region_from_geometry(
            "region_c",
            Polygon([(1.1, 0.1), (1.9, 0.1), (1.9, 0.9), (1.1, 0.9)]),
            tags=["c"],
        )

        group_1 = model.add_group("group_1", members=["region_a", "region_b"], tags=["group"])
        group_2 = model.add_group("group_2", members=["group_1", "region_c"], tags=["group"])

        assert isinstance(region_a, ModelRegion)
        assert isinstance(group_1, RegionGroup)
        assert model.list_regions()["name"].tolist() == ["region_a", "region_b", "region_c"]
        assert model.list_groups()["name"].tolist() == ["group_1", "group_2"]

        assert model.get_region_cells("region_a") == [0]
        assert model.get_region_cells("region_b") == [0, 1]
        assert model.get_region_cells("group_1") == [0, 1]
        assert model.get_region_cells("group_2") == [0, 1]
        assert model.resolve_region_cells("group_2") == [0, 1]

        resolved_cells, trace = model.resolve_region_cells_with_trace("group_2")
        assert resolved_cells == [0, 1]
        assert trace[0] == ["region_a", "region_b"]
        assert trace[1] == ["region_b", "region_c"]

        group_heads = model.region_heads("group_2", per=0)
        assert group_heads["region"].unique().tolist() == ["group_2"]
        assert sorted(group_heads["cell"].tolist()) == [0, 1]
        assert np.allclose(sorted(group_heads["elev"].astype(float).to_numpy()), [9.0, 10.0])

        model.remove_from_group("group_1", "region_b")
        assert model.get_group("group_1").members == ["region_a"]
        assert model.get_region_cells("group_1") == [0]

        model.add_to_group("group_1", "region_b")
        assert model.get_region_cells("group_1") == [0, 1]

        with pytest.raises(ValueError, match="cyclic group relationship"):
            model.add_to_group("group_1", "group_2")
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simple_model_spec_can_build_and_run_through_project():
    workspace = _project_temp_dir("simple_model_factory_run")
    try:
        vor = _two_cell_vor_clockwise()
        config = SimpleModelConfig(
            vor=vor,
            name="factory_smoke",
            nper=1,
            nlay=1,
            grid_type="disv",
            top=[10.0, 9.0],
            bottom=[[0.0, 0.0]],
            initial_heads=[10.0, 9.0],
            k=[1.0, 1.0],
            save_specific_discharge=False,
            boundary_mode="chd",
            boundary_cells=[0, 1],
            boundary_head=[10.0, 9.0],
            sto_steady={0: True},
            sto_transient={},
        )

        spec = simple_model_spec(config)
        project = Project(workspace / "project")
        run = project.run("baseline", spec)

        assert run.success is True
        assert spec.models[0].name == "factory_smoke"
        assert spec.models[0].context.grid is vor
        assert spec.models[0].context.metadata["grid_type"] == "disv"
        assert [package.name for package in spec.models[0].packages] == [
            "disv",
            "ic",
            "npf",
            "sto",
            "oc",
            "chd",
        ]
        assert (run.workspace / "factory_smoke.chd").exists()

        heads = run.simulation.get_model("factory_smoke").output.head().get_data(
            kstpkper=(0, 0)
        ).squeeze()
        assert np.allclose(heads, [10.0, 9.0])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simple_model_spec_exposes_replaceable_package_specs():
    workspace = _project_temp_dir("simple_model_class_smoke")
    try:
        vor = _two_cell_vor_clockwise()
        config = SimpleModelConfig(
            vor=vor,
            name="class_smoke",
            top=[10.0, 10.0],
            bottom=[[0.0, 0.0]],
            k=[1.0, 1.0],
            boundary_mode="drain",
            boundary_conductance=5.0,
        )
        baseline = simple_model_spec(config)
        flow = baseline.model("class_smoke")
        npf = flow.package("npf")
        high_k = baseline.with_model(flow.with_package(npf.with_options(k=[25.0, 25.0])))

        assert [package.name for package in flow.packages][-1] == "drn"
        assert flow.package("drn").options["stress_period_data"] == [
            [(0, 0), 0.1, 5.0],
            [(0, 1), 0.1, 5.0],
        ]
        assert high_k.model("class_smoke").package("npf").options["k"] == [25.0, 25.0]
        assert npf.options["k"] == [1.0, 1.0]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simple_model_spec_can_run_disu_grid():
    workspace = _project_temp_dir("simple_model_disu_spec")
    try:
        config = SimpleModelConfig(
            vor=_two_cell_vor_clockwise(),
            name="disu_spec",
            grid_type="disu",
            top=[10.0, 9.0],
            bottom=[0.0, 0.0],
            idomain=np.array([1, 1]),
            initial_heads=[10.0, 9.0],
            k=[1.0, 1.0],
            save_specific_discharge=False,
            boundary_mode="chd",
            boundary_cells=[0, 1],
            boundary_head=[10.0, 9.0],
            sto_transient={},
        )

        spec = simple_model_spec(config)
        run = Project(workspace / "project").run("baseline", spec)

        assert run.success is True
        assert spec.model("disu_spec").packages[0].name == "disu"
        assert (run.workspace / "disu_spec.disu").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simple_model_config_validates_model_name_length():
    with pytest.raises(ValueError, match="16 characters or fewer"):
        SimpleModelConfig(vor=_two_cell_vor(), name="name_too_long_for_mf6")


def test_rch_and_uzf_can_run_together():
    workspace = _project_temp_dir("myflopy_rch_uzf_run")
    try:
        vor = _two_cell_vor_clockwise()
        model = SimulationBase(name="rch_uzf_smoke", mf_folder_path=workspace, vor=vor, nper=2)

        top = [10.0, 10.0]
        botm = [[0.0, 0.0]]
        DisvGrid(vor=vor, model=model, top=top, bottom=botm, nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 10.0])
        KFlow(model=model, k=[1.0, 1.0], k33_vert=[1.0, 1.0])
        Storage(model=model, sto_steady={0: True}, sto_transient={1: True})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [[(0, 0), 10.0], [(0, 1), 10.0]],
                1: [[(0, 0), 10.0], [(0, 1), 10.0]],
            },
        )

        recharge_path = _write_gpkg(
            workspace / "recharge.gpkg",
            gpd.GeoDataFrame(
                {
                    "zone": ["left"],
                    "rch_0": [0.02],
                    "rch_1": [0.03],
                },
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
                crs=vor.crs,
            ),
        )

        recharge_builder = RCHFromVector(
            model=model,
            vor=vor,
            shp_gpkg=recharge_path,
            uid="zone",
            rch_fields=["rch_0", "rch_1"],
            rch_fields_to_pers=[0, 1],
            background_rch=0.0,
            grid_type="disv",
            limit_to_k33=False,
        )
        rch_dict = recharge_builder.from_vector(
            register_regions=True,
            region_name_prefix="rch_zone",
            combined_region_name="all_rch_zones",
            region_tags=["rch"],
            overwrite_regions=True,
        )
        Recharge(model=model, vor=vor, rch_dict=rch_dict)

        uzf_context = _model_context(model, vor)
        uzf = UZFBuilder(
            context=uzf_context,
            nper=model.nper,
            vks=1.0,
            thtr=0.1,
            thts=0.3,
            thti=0.2,
            finf=_uzf_finf_from_rch(rch_dict, _surface_cells(uzf_context), model.nper),
        )
        uzf.build().build(model.gwf)
        model.add_region_from_cells(
            "uzf_all", uzf.uzf_cells, category="boundary", package="uzf",
            tags=["uzf"], overwrite=True,
        )

        success, _ = model.run_simulation()
        assert success is True
        assert np.allclose(uzf.finf[0], [0.019, 0.0])
        assert np.allclose(uzf.finf[1], [0.0285, 0.0])
        assert model.get_region_cells("rch_zone_left") == [0]
        assert set(model.get_region_cells("all_rch_zones")) == {0}
        assert set(model.get_region_cells("uzf_all")) == {0, 1}

        outdir = workspace / "rch_uzf_smoke"
        assert (outdir / "rch_uzf_smoke.uzf").exists()
        assert (outdir / "rch_uzf_smoke.rch").exists()
        assert (outdir / "rch_uzf_smoke_budget.uzf").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_minimal_disu_model_can_run():
    workspace = _project_temp_dir("myflopy_disu_run")
    try:
        vor = _two_cell_vor_clockwise()
        vor.get_disu_connectivity(validate=True)
        model = SimulationBase(name="disu_smoke", mf_folder_path=workspace, vor=vor, nper=1)
        model.nlay = 1

        DisuGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[0.0, 0.0])
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [
                    [0, 10.0],
                    [1, 9.0],
                ]
            },
        )

        success, _ = model.run_simulation()
        assert success is True

        outdir = workspace / "disu_smoke"
        assert (outdir / "disu_smoke.disu").exists()
        heads = model.gwf.output.head().get_data(kstpkper=(0, 0)).squeeze()
        assert np.allclose(heads, [10.0, 9.0])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_lak_and_sfr_can_run_together():
    workspace = _project_temp_dir("myflopy_lak_sfr_run")
    try:
        vor = _four_cell_vor_clockwise()
        vor.gdf_topbtm = gpd.GeoDataFrame(
            {
                "geometry": vor.gdf_vorPolys.geometry,
                0: [12.0, 12.0, 12.0, 12.0],
                1: [0.0, 0.0, 0.0, 0.0],
            },
            geometry="geometry",
            crs=vor.crs,
        )
        model = SimulationBase(name="lak_sfr_smoke", mf_folder_path=workspace, vor=vor, nper=1)

        DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[11.0] * 4)
        KFlow(model=model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [[(0, 3), 11.0]]
            },
        )
        model.vor.get_disu_connectivity(validate=True)

        lake_path = _write_gpkg(
            workspace / "lake.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=vor.crs,
            ),
        )
        stream_path = _write_gpkg(
            workspace / "stream.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream_0"]},
                geometry=[LineString([(0.1, 1.5), (1.9, 1.5)])],
                crs=vor.crs,
            ),
        )

        lak = _attach_lak(
            model,
            vor,
            [lake_path],
            starting_stage=11.0,
            lake_bottom=9.0,
            bed_leakance=0.1,
            status="ACTIVE",
            region_name_prefix="lake_zone",
            combined_region_name="all_lakes",
            region_tags=["lak"],
            overwrite_regions=True,
        )
        sfr = _attach_sfr(
            model,
            vor,
            [stream_path],
            widths=1.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            region_name_prefix="sfr_group",
            combined_region_name="all_streams",
            region_tags=["sfr"],
            overwrite_regions=True,
        )
        sfr_package_spec = sfr.build()

        success, _ = model.run_simulation()
        assert success is True

        outdir = workspace / "lak_sfr_smoke"
        assert (outdir / "lak_sfr_smoke.lak").exists()
        assert (outdir / "lak_sfr_smoke.sfr").exists()
        assert sfr_package_spec.name == "sfr"
        assert sfr_package_spec.options["nreaches"] == sfr.total_nreaches
        assert len(lak.connectiondata) >= 1
        assert sfr.total_nreaches >= 1
        assert model.get_region_cells("lake_zone_lake_0") == [0]
        assert set(model.get_region_cells("all_lakes")) == {0}
        assert set(model.get_region_cells("sfr_group_stream")) == {2, 3}
        assert set(model.get_region_cells("all_streams")) == {2, 3}
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_small_vector_sfr_lak_mvr_workflow_can_run():
    workspace = _project_temp_dir("myflopy_sfr_lak_mvr_run")
    try:
        vor = _four_cell_vor_clockwise()
        vor.gdf_topbtm = gpd.GeoDataFrame(
            {
                "geometry": vor.gdf_vorPolys.geometry,
                0: [12.0, 12.0, 12.0, 12.0],
                1: [0.0, 0.0, 0.0, 0.0],
            },
            geometry="geometry",
            crs=vor.crs,
        )
        model = SimulationBase(name="sfr_lak_mvr", mf_folder_path=workspace, vor=vor, nper=1)

        DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[11.2] * 4)
        KFlow(model=model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(model=model, stress_period_data={0: [[(0, 3), 11.0]]})

        lake_path = _write_gpkg(
            workspace / "lake_mvr.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=vor.crs,
            ),
        )
        stream_path = _write_gpkg(
            workspace / "stream_mvr.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream_0"]},
                geometry=[LineString([(0.2, 1.8), (1.8, 1.2)])],
                crs=vor.crs,
            ),
        )

        lak = _attach_lak(
            model,
            vor,
            [lake_path],
            starting_stage=11.0,
            lake_bottom=9.0,
            bed_leakance=0.1,
            status="ACTIVE",
            mover=True,
            region_name_prefix="lake_zone",
            combined_region_name="all_lakes",
            overwrite_regions=True,
        )

        sfr = _attach_sfr(
            model,
            vor,
            [stream_path],
            inflows={0: [(0, 0.5)]},
            widths=5.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            mover=True,
            region_name_prefix="sfr_group",
            combined_region_name="all_streams",
            overwrite_regions=True,
        )

        mvr = MVRBuilder(
            nper=1,
            moves={
                0: [
                    Move(
                        source=sfr.connection(sfr.stream_ids[0]),
                        receiver=lak.connection(lak.lake_ids[0]),
                    )
                ]
            },
        )
        validation_report = model.validate_surface_water(
            nlakes=1,
            lak_packagedata=lak.packagedata,
            lak_connectiondata=lak.connectiondata,
            lak_perioddata=lak.perioddata,
            sfr=sfr,
            maxmvr=1,
            maxpackages=2,
            mvr_packages=mvr.packages,
            mvr_perioddata=mvr.perioddata,
            raise_on_error=True,
        )
        assert validation_report.ok is True
        mvr.build().build(model.gwf)

        success, _ = model.run_simulation()
        assert success is True

        outdir = workspace / "sfr_lak_mvr"
        assert (outdir / "sfr_lak_mvr.lak").exists()
        assert (outdir / "sfr_lak_mvr.sfr").exists()
        assert (outdir / "sfr_lak_mvr.mvr").exists()
        assert (outdir / "sfr_lak_mvr.mvr.bud").exists()
        assert sfr.total_nreaches >= 1
        assert model.get_region_cells("all_lakes")
        assert model.get_region_cells("all_streams")
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_feature_rich_small_model_workflow_can_run():
    workspace = _project_temp_dir("myflopy_feature_rich_workflow")
    try:
        vor = _four_cell_vor_clockwise()
        vor.gdf_topbtm = gpd.GeoDataFrame(
            {
                "geometry": vor.gdf_vorPolys.geometry,
                0: [12.0, 12.0, 12.0, 12.0],
                1: [0.0, 0.0, 0.0, 0.0],
            },
            geometry="geometry",
            crs=vor.crs,
        )
        model = SimulationBase(name="feat_rich_demo", mf_folder_path=workspace, vor=vor, nper=2)

        DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[11.0] * 4)
        KFlow(model=model, k=[1.0] * 4, k33_vert=[1.0] * 4, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={1: True})
        OutputControl(model=model)
        CHD(
            model=model,
            stress_period_data={
                0: [[(0, 1), 10.5]],
                1: [[(0, 1), 10.25]],
            },
        )

        drain_path = _write_gpkg(
            workspace / "drain.gpkg",
            gpd.GeoDataFrame(
                {
                    "name": ["northwest_drain"],
                    "height": [1.0],
                    "cond": [5.0],
                    "layer": [1],
                    "min_elev": [2.0],
                },
                geometry=[Polygon([(0.0, 1.0), (0.95, 1.0), (0.95, 2.0), (0.0, 2.0)])],
                crs=vor.crs,
            ),
        )
        drn_builder = DRNFromVector(model=model, vor=vor, shp_gpkg=drain_path, uid="name", idomain=[1, 1, 1, 1])
        drn_dict = drn_builder.from_vector(
            edges_only=True,
            register_regions=True,
            region_name_prefix="drn_group",
            combined_region_name="all_drains",
            region_tags=["drn"],
            overwrite_regions=True,
        )
        myflopy.modflow.mf6.simulation.packages.Drains(model=model, stress_period_data=drn_dict)

        ghb_path = _write_gpkg(
            workspace / "ghb.gpkg",
            gpd.GeoDataFrame(
                {
                    "name": ["northeast_ghb"],
                    "elev": [11.25],
                    "height": [0.0],
                    "cond": [7.0],
                    "layer": [1],
                    "min_elev": [10.0],
                },
                geometry=[Polygon([(1.05, 1.0), (2.0, 1.0), (2.0, 2.0), (1.05, 2.0)])],
                crs=vor.crs,
            ),
        )
        ghb_builder = GHBFromVector(model=model, vor=vor, shp_gpkg=ghb_path, uid="name", idomain=[1, 1, 1, 1])
        ghb_dict = ghb_builder.from_vector(
            register_regions=True,
            region_name_prefix="ghb_group",
            combined_region_name="all_ghb",
            region_tags=["ghb"],
            overwrite_regions=True,
        )
        myflopy.modflow.mf6.simulation.packages.GHB(model=model, stress_period_data=ghb_dict)

        recharge_path = _write_gpkg(
            workspace / "recharge.gpkg",
            gpd.GeoDataFrame(
                {
                    "zone": ["left_recharge", "right_recharge"],
                    "rch_0": [0.015, 0.01],
                    "rch_1": [0.02, 0.012],
                },
                geometry=[
                    Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)]),
                    Polygon([(1.05, 0.0), (2.0, 0.0), (2.0, 1.0), (1.05, 1.0)]),
                ],
                crs=vor.crs,
            ),
        )
        recharge_builder = RCHFromVector(
            model=model,
            vor=vor,
            shp_gpkg=recharge_path,
            uid="zone",
            rch_fields=["rch_0", "rch_1"],
            rch_fields_to_pers=[0, 1],
            background_rch=0.0,
            grid_type="disv",
            limit_to_k33=False,
        )
        rch_dict = recharge_builder.from_vector(
            register_regions=True,
            region_name_prefix="rch_zone",
            combined_region_name="all_rch",
            region_tags=["rch"],
            overwrite_regions=True,
        )
        Recharge(model=model, vor=vor, rch_dict=rch_dict)
        uzf_context = _model_context(model, vor)
        uzf = UZFBuilder(
            context=uzf_context,
            nper=model.nper,
            vks=1.0,
            thtr=0.1,
            thts=0.3,
            thti=0.2,
            finf=_uzf_finf_from_rch(rch_dict, _surface_cells(uzf_context), model.nper),
        )
        uzf.build().build(model.gwf)
        model.add_region_from_cells(
            "uzf_all", uzf.uzf_cells, category="boundary", package="uzf",
            tags=["uzf"], overwrite=True,
        )

        lake_path = _write_gpkg(
            workspace / "lake.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=vor.crs,
            ),
        )
        lak = _attach_lak(
            model,
            vor,
            [lake_path],
            starting_stage=11.0,
            lake_bottom=9.0,
            bed_leakance=0.1,
            status=["ACTIVE", "ACTIVE"],
            region_name_prefix="lake_zone",
            combined_region_name="all_lakes",
            region_tags=["lak"],
            overwrite_regions=True,
        )

        stream_path = _write_gpkg(
            workspace / "stream.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream_0"]},
                geometry=[LineString([(0.1, 1.5), (1.9, 1.5)])],
                crs=vor.crs,
            ),
        )
        _attach_sfr(
            model,
            vor,
            [stream_path],
            widths=1.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            region_name_prefix="sfr_group",
            combined_region_name="all_streams",
            region_tags=["sfr"],
            overwrite_regions=True,
        )

        model.add_group(
            "boundary_features",
            members=["all_drains", "all_ghb", "all_rch", "uzf_all", "all_lakes", "all_streams"],
        )

        success, _ = model.run_simulation()
        assert success is True

        resolved_cells, trace = model.resolve_region_cells_with_trace("boundary_features")
        heads = model.region_heads("boundary_features", per=1)

        assert resolved_cells == [0, 1, 2, 3]
        assert trace[0] == ["all_drains", "all_lakes", "all_rch", "uzf_all"]
        assert trace[1] == ["all_ghb", "all_rch", "uzf_all"]
        assert trace[2] == ["all_drains", "all_rch", "all_streams", "uzf_all"]
        assert trace[3] == ["all_ghb", "all_rch", "all_streams", "uzf_all"]
        assert heads["region"].unique().tolist() == ["boundary_features"]
        assert sorted(heads["cell"].tolist()) == [0, 1, 2, 3]

        outdir = workspace / "feat_rich_demo"
        assert (outdir / "feat_rich_demo.drn").exists()
        assert (outdir / "feat_rich_demo.ghb").exists()
        assert (outdir / "feat_rich_demo.rch").exists()
        assert (outdir / "feat_rich_demo.uzf").exists()
        assert (outdir / "feat_rich_demo.lak").exists()
        assert (outdir / "feat_rich_demo.sfr").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_refined_end_to_end_model_can_run_with_preferred_builder_api():
    workspace = _project_temp_dir("refined_end_to_end")
    try:
        _, vor, _, inputs = _build_refined_vector_test_triangle(workspace)
        top, botm = _refined_vor_top_bottom(vor)
        top_arr = np.asarray(top, dtype=float)
        botm_arr = np.asarray(botm, dtype=float)
        initial_heads = np.maximum(botm_arr + 10.0, top_arr - 8.0).tolist()

        model = SimulationBase(name="refined_e2e", mf_folder_path=workspace, vor=vor, nper=2)
        DisvGrid(vor=vor, model=model, top=top, bottom=[botm], nlay=1, idomain=[[1] * vor.ncpl])
        TemporalDiscretization(
            model=model,
            period_data=[
                [1.0, 1, 1.0],
                [1.0, 1, 1.0],
            ],
        )
        InitialConditions(model=model, vor=vor, nlay=1, strt=initial_heads)

        k_builder = KFromVector(model=model, vor=vor, shp_gpkg=inputs["k"], uid="name")
        k_array = k_builder.from_vector(defaults=[18.0])
        KFlow(
            model=model,
            k=k_array[0].tolist(),
            k33_vert=(k_array[0] * 0.2).tolist(),
            save_specific_discharge=False,
        )
        Storage(model=model, sto_steady={0: True}, sto_transient={1: True})
        OutputControl(model=model)

        west_stage = float(np.nanpercentile(top_arr, 92) - 5.0)
        east_stage = west_stage - 20.0

        chd_builder = CHDFromVector(model=model, vor=vor, shp_gpkg=inputs["chd"], uid="name")
        chd_dict = chd_builder.from_vector(
            head_reference={"west_chd": [west_stage, west_stage]},
            register_regions=True,
            region_name_prefix="chd_group",
            combined_region_name="all_chd",
            overwrite_regions=True,
        )
        CHD(model=model, stress_period_data=chd_dict)

        ghb_builder = GHBFromVector(model=model, vor=vor, shp_gpkg=inputs["ghb"], uid="name")
        ghb_dict = ghb_builder.from_vector(
            elev_reference={"east_ghb": [east_stage, east_stage - 0.5]},
            register_regions=True,
            region_name_prefix="ghb_group",
            combined_region_name="all_ghb",
            overwrite_regions=True,
        )
        myflopy.modflow.mf6.simulation.packages.GHB(model=model, stress_period_data=ghb_dict)

        drn_builder = DRNFromVector(model=model, vor=vor, shp_gpkg=inputs["drn"], uid="name")
        drn_dict = drn_builder.from_vector(
            edges_only=False,
            top_drain=True,
            register_regions=True,
            region_name_prefix="drn_group",
            combined_region_name="all_drains",
            overwrite_regions=True,
        )
        myflopy.modflow.mf6.simulation.packages.Drains(model=model, stress_period_data=drn_dict)

        recharge_builder = RCHFromVector(
            model=model,
            vor=vor,
            shp_gpkg=inputs["rch"],
            uid="zone",
            rch_fields=["rch_0", "rch_1"],
            rch_fields_to_pers=[0, 1],
            background_rch=0.0,
            grid_type="disv",
            limit_to_k33=False,
        )
        rch_dict = recharge_builder.from_vector(
            register_regions=True,
            region_name_prefix="rch_zone",
            combined_region_name="all_rch",
            overwrite_regions=True,
        )
        Recharge(model=model, vor=vor, rch_dict=rch_dict)
        uzf_context = _model_context(model, vor)
        uzf = UZFBuilder(
            context=uzf_context,
            nper=model.nper,
            vks=0.05,
            thtr=0.08,
            thts=0.28,
            thti=0.18,
            finf=_uzf_finf_from_rch(rch_dict, _surface_cells(uzf_context), model.nper),
        )
        uzf.build().build(model.gwf)
        model.add_region_from_cells(
            "uzf_all", uzf.uzf_cells, category="boundary", package="uzf", overwrite=True,
        )

        lak = _attach_lak(
            model,
            vor,
            [inputs["lakes"][0]],
            starting_stage=129.0,
            lake_bottom=118.0,
            bed_leakance=0.001,
            status=["ACTIVE", "ACTIVE"],
            mover=True,
            region_name_prefix="lake_zone",
            combined_region_name="all_lakes",
            region_tags=["lak"],
            overwrite_regions=True,
        )

        sfr = _attach_sfr(
            model,
            vor,
            [inputs["streams"][0]],
            inflows={
                0: [(0, 0.005)],
                1: [(0, 0.005)],
            },
            widths=5.0,
            gradients=0.0008,
            mannings=0.03,
            streambed_k=0.05,
            streambed_thickness=2.0,
            mover=True,
            region_name_prefix="sfr_group",
            combined_region_name="all_streams",
            overwrite_regions=True,
        )
        mvr = MVRBuilder(
            nper=2,
            moves={
                period: [
                    Move(
                        source=sfr.connection(sfr.stream_ids[0]),
                        receiver=lak.connection(lak.lake_ids[0]),
                        value=0.25,
                    )
                ]
                for period in range(2)
            },
        )
        validation_report = model.validate_surface_water(
            nlakes=1,
            lak_packagedata=lak.packagedata,
            lak_connectiondata=lak.connectiondata,
            lak_perioddata=lak.perioddata,
            sfr=sfr,
            maxmvr=1,
            maxpackages=2,
            mvr_packages=mvr.packages,
            mvr_perioddata=mvr.perioddata,
            raise_on_error=True,
        )
        assert validation_report.ok is True
        assert validation_report.summary()["num_errors"] == 0
        mvr.build().build(model.gwf)

        success, _ = model.run_simulation()
        assert success is True

        outdir = workspace / "refined_e2e"
        assert vor.ncpl >= 700
        assert len(chd_dict[0]) > 0
        assert len(ghb_dict[0]) > 0
        assert len(drn_dict[0]) > 0
        assert len(rch_dict[0]) > 0
        assert len(lak.connectiondata) >= 1
        assert sfr.total_nreaches > 5
        assert len(uzf.packagedata) == vor.ncpl
        assert list(outdir.glob("*.disv"))
        assert list(outdir.glob("*.chd"))
        assert list(outdir.glob("*.ghb"))
        assert list(outdir.glob("*.drn"))
        assert list(outdir.glob("*.rch"))
        assert list(outdir.glob("*.uzf"))
        assert list(outdir.glob("*.lak"))
        assert list(outdir.glob("*.sfr"))
        assert list(outdir.glob("*.mvr"))
        assert model.get_region_cells("all_chd")
        assert model.get_region_cells("all_ghb")
        assert model.get_region_cells("all_drains")
        assert model.get_region_cells("all_rch")
        assert model.get_region_cells("all_lakes")
        assert model.get_region_cells("all_streams")
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_surface_water_validation_catches_invalid_mvr_source_index():
    workspace = _project_temp_dir("surface_water_validation_invalid_mvr")
    try:
        vor = _four_cell_vor_clockwise()
        vor.gdf_topbtm = gpd.GeoDataFrame(
            {
                "geometry": vor.gdf_vorPolys.geometry,
                0: [12.0, 12.0, 12.0, 12.0],
                1: [0.0, 0.0, 0.0, 0.0],
            },
            geometry="geometry",
            crs=vor.crs,
        )
        model = SimulationBase(name="surface_water_validation_invalid", mf_folder_path=workspace, vor=vor, nper=1)
        DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[11.2] * 4)
        KFlow(model=model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(model=model, stress_period_data={0: [[(0, 3), 11.0]]})

        lake_path = _write_gpkg(
            workspace / "lake_invalid.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=vor.crs,
            ),
        )
        stream_path = _write_gpkg(
            workspace / "stream_invalid.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream_0"]},
                geometry=[LineString([(0.2, 1.8), (1.8, 1.2)])],
                crs=vor.crs,
            ),
        )

        lak = _attach_lak(
            model,
            vor,
            [lake_path],
            starting_stage=11.0,
            lake_bottom=9.0,
            bed_leakance=0.1,
            status="ACTIVE",
            mover=True,
        )

        sfr = _attach_sfr(
            model,
            vor,
            [stream_path],
            inflows={0: [(0, 0.5)]},
            widths=5.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            mover=True,
        )

        mvr = MVRBuilder(
            nper=1,
            moves={
                0: [
                    Move(
                        source=MoverConnection("sfr", sfr.total_nreaches),
                        receiver=lak.connection(lak.lake_ids[0]),
                    )
                ]
            },
        )
        report = model.validate_surface_water(
            nlakes=1,
            lak_packagedata=lak.packagedata,
            lak_connectiondata=lak.connectiondata,
            lak_perioddata=lak.perioddata,
            sfr=sfr,
            maxmvr=1,
            maxpackages=2,
            mvr_packages=mvr.packages,
            mvr_perioddata=mvr.perioddata,
        )

        assert report.ok is False
        error_codes = {issue.code for issue in report.errors}
        assert "mvr_source_range" in error_codes
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_code_defined_geometries_can_drive_refinement_regions_and_packages():
    workspace = _project_temp_dir("code_geometry_workflow")
    try:
        domain_geom = Polygon([(0, 0), (600, 0), (600, 450), (0, 450)])
        lake_geom = Point(220, 170).buffer(60)
        stream_line = LineString([(50, 360), (260, 270), (540, 140)])
        drain_geom = Polygon([(0, 330), (180, 330), (180, 450), (0, 450)])
        ghb_geom = Polygon([(0, 0), (180, 0), (180, 120), (0, 120)])
        uplands_geom = Polygon([(0, 220), (600, 220), (600, 450), (0, 450)])
        lowlands_geom = Polygon([(0, 0), (600, 0), (600, 220), (0, 220)])

        tri = TriangleGrid(model_ws=str(workspace / "triangle_build"))
        tri.set_domain_polygon(domain_geom)
        tri.add_region_polygon(
            Polygon([(40, 40), (560, 40), (560, 410), (40, 410)]),
            max_area=18000,
            label="mid_refine",
        )
        tri.add_region_polygon(
            stream_line.buffer(30),
            max_area=4000,
            label="stream_refine",
            priority=2,
            source="line",
        )
        tri.add_region_polygon(
            lake_geom,
            max_area=1200,
            label="lake_refine",
            priority=3,
            source="circle",
        )
        tri.build(verbose=False)

        vor = VoronoiGridPlus(tri)
        top = 115.0 - (0.02 * np.asarray(vor.centroids_x)) + (0.01 * np.asarray(vor.centroids_y))
        bottom = top - 30.0
        vor.gdf_topbtm = gpd.GeoDataFrame(
            {
                "geometry": vor.gdf_vorPolys.geometry,
                0: top,
                1: bottom,
            },
            geometry="geometry",
            crs=vor.crs,
        )

        model = SimulationBase(name="geom_shapes", mf_folder_path=workspace, vor=vor, nper=1)
        DisvGrid(vor=vor, model=model, top=top.tolist(), bottom=[bottom.tolist()], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=(top - 5.0).tolist())
        KFlow(model=model, k=[15.0] * vor.ncpl, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)

        east_strip = Polygon([(560, 0), (600, 0), (600, 450), (560, 450)])
        east_cells = sorted(set(vor.get_vor_cells_as_series(east_strip).iloc[0]))
        CHD(model=model, stress_period_data={0: [[(0, cell), 100.0] for cell in east_cells]})

        lake_region = model.add_region_from_geometry(
            "lake_circle_region",
            lake_geom,
            category="custom",
            overwrite=True,
        )
        stream_region = model.add_region_from_geometry(
            "stream_corridor_region",
            stream_line.buffer(30),
            category="custom",
            overwrite=True,
        )

        drain_path = _write_features(
            workspace / "drain.gpkg",
            [{"name": "northwest_drain", "height": 2.0, "cond": 500.0, "layer": 1, "min_elev": 80.0, "geometry": drain_geom}],
            vor.crs,
        )
        ghb_path = _write_features(
            workspace / "ghb.gpkg",
            [{"name": "southwest_ghb", "elev": 90.0, "height": 0.0, "cond": 650.0, "layer": 1, "min_elev": 78.0, "geometry": ghb_geom}],
            vor.crs,
        )
        recharge_path = _write_features(
            workspace / "recharge.gpkg",
            [
                {"zone": "uplands", "rch_0": 0.002, "geometry": uplands_geom},
                {"zone": "lowlands", "rch_0": 0.0015, "geometry": lowlands_geom},
            ],
            vor.crs,
        )
        lake_path = _write_features(
            workspace / "lake.gpkg",
            [{"name": "lake_0", "geometry": lake_geom}],
            vor.crs,
        )
        stream_path = _write_features(
            workspace / "stream.gpkg",
            [{"name": "stream_0", "geometry": stream_line}],
            vor.crs,
        )

        drn_builder = DRNFromVector(model=model, vor=vor, shp_gpkg=drain_path, uid="name", idomain=[1] * vor.ncpl)
        drn_dict = drn_builder.from_vector(
            edges_only=True,
            register_regions=True,
            region_name_prefix="drn_group",
            combined_region_name="all_drains",
            overwrite_regions=True,
        )
        myflopy.modflow.mf6.simulation.packages.Drains(model=model, stress_period_data=drn_dict)

        ghb_builder = GHBFromVector(model=model, vor=vor, shp_gpkg=ghb_path, uid="name", idomain=[1] * vor.ncpl)
        ghb_dict = ghb_builder.from_vector(
            register_regions=True,
            region_name_prefix="ghb_group",
            combined_region_name="all_ghb",
            overwrite_regions=True,
        )
        myflopy.modflow.mf6.simulation.packages.GHB(model=model, stress_period_data=ghb_dict)

        recharge_builder = RCHFromVector(
            model=model,
            vor=vor,
            shp_gpkg=recharge_path,
            uid="zone",
            rch_fields=["rch_0"],
            rch_fields_to_pers=[0],
            background_rch=0.0,
            grid_type="disv",
            limit_to_k33=False,
        )
        rch_dict = recharge_builder.from_vector(
            register_regions=True,
            region_name_prefix="rch_zone",
            combined_region_name="all_rch",
            overwrite_regions=True,
        )
        Recharge(model=model, vor=vor, rch_dict=rch_dict)
        uzf_context = _model_context(model, vor)
        uzf = UZFBuilder(
            context=uzf_context,
            nper=model.nper,
            vks=0.5,
            thtr=0.1,
            thts=0.3,
            thti=0.2,
            finf=_uzf_finf_from_rch(rch_dict, _surface_cells(uzf_context), model.nper),
        )
        uzf.build().build(model.gwf)
        model.add_region_from_cells(
            "uzf_all", uzf.uzf_cells, category="boundary", package="uzf", overwrite=True,
        )

        lak = _attach_lak(
            model,
            vor,
            [lake_path],
            starting_stage=100.0,
            lake_bottom=92.0,
            bed_leakance=0.05,
            status="ACTIVE",
            region_name_prefix="lake_zone",
            combined_region_name="all_lakes",
            overwrite_regions=True,
        )

        sfr = _attach_sfr(
            model,
            vor,
            [stream_path],
            widths=12.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=2.0,
            streambed_thickness=1.5,
            region_name_prefix="sfr_group",
            combined_region_name="all_streams",
            overwrite_regions=True,
        )

        assert len(lake_region.cells) > 0
        assert len(stream_region.cells) > 0
        assert len(drn_dict[0]) > 0
        assert len(ghb_dict[0]) > 0
        assert len(rch_dict[0]) > 0
        assert len(uzf.packagedata) > 0
        assert len(lak.connectiondata) > 0
        assert len(lak.packagedata) == 1
        assert len(sfr.stream_cells) == 1
        assert model.get_region_cells("all_lakes") == sorted(set(model.get_region_cells("all_lakes")))
        assert model.get_region_cells("lake_circle_region") == sorted(set(model.get_region_cells("lake_circle_region")))
        assert set(model.get_region_cells("all_streams")).issubset(set(stream_region.cells))

        success, _ = model.run_simulation()
        assert success is True
        assert (model.model_output_folder_path / f"{model.name}.hds").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_surface_reconciliation_helpers_smoke():
    vor = _two_cell_vor()
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: [10.0, 10.0],
            1: [10.5, 9.0],
            2: [7.0, 6.0],
        },
        geometry="geometry",
        crs=vor.crs,
    )

    reconciled = vor.reconcile_surfaces(min_sep=0.5, trigger_sep=1.0)
    assert reconciled.loc[0, 1] == 9.5
    assert reconciled.loc[0, 2] == 7.0

    adjusted = vor.adjust_cells_by_id([1], adjustment=-2.0, layer=0, reconcile=False)
    assert adjusted.loc[1, 0] == 8.0


def test_triangle_grid_helpers_smoke():
    workspace = _project_temp_dir("triangle_grid_smoke")
    try:
        tri = TriangleGrid(model_ws=str(workspace))

        rect = tri.add_rectangle(x_dist=2, y_dist=1, return_only=True)
        circle = tri.add_circle(radius=1, return_only=True, radians_step=0.5)
        clouds = tri.generate_dissipating_point_cloud(
            polygon=rect,
            buffer_dist=3,
            num_buffers=3,
            min_spacing=1,
            max_spacing=2,
            method="power",
        )

        assert np.isclose(rect.area, 2.0)
        assert circle.area > 2.5
        assert len(clouds) == 4

        tri.add_polygon(rect, max_area=1.0)
        prepared = tri.prepare()
        assert tri.domain_geometry.equals(rect)
        assert len(tri.region_specs) >= 1
        assert len(prepared) >= 1
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_triangle_grid_domain_and_region_api_smoke():
    workspace = _project_temp_dir("triangle_grid_domain_region_smoke")
    try:
        tri = TriangleGrid(model_ws=str(workspace))
        tri.set_domain_rectangle(x_dist=10, y_dist=8, origin=(0, 0))
        tri.add_region_rectangle(
            origin=(2, 2),
            x_dist=6,
            y_dist=4,
            max_area=4.0,
            label="outer",
        )
        tri.add_region_rectangle(
            origin=(4, 3),
            x_dist=2,
            y_dist=1,
            max_area=1.0,
            label="inner",
            priority=1,
        )

        preview = tri.preview_regions()
        points = tri.get_region_points()

        assert np.isclose(tri.domain_geometry.area, 80.0)
        assert set(preview["label"]) == {"outer", "inner"}
        assert len(points) == 2

        outer_geom = preview.set_index("label").loc["outer", "geometry"]
        inner_geom = preview.set_index("label").loc["inner", "geometry"]
        outer_point = points.set_index("label").loc["outer", "geometry"]
        inner_point = points.set_index("label").loc["inner", "geometry"]

        assert tri.domain_geometry.covers(outer_point)
        assert tri.domain_geometry.covers(inner_point)
        assert outer_geom.covers(outer_point)
        assert inner_geom.covers(inner_point)
        assert not inner_geom.covers(outer_point)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_triangle_grid_file_driven_region_setup_smoke():
    workspace = _project_temp_dir("triangle_grid_file_smoke")
    try:
        domain_path = _write_gpkg(
            workspace / "domain.gpkg",
            gpd.GeoDataFrame(
                {"name": ["domain"]},
                geometry=[Polygon([(0, 0), (6, 0), (6, 6), (0, 6)])],
                crs="EPSG:2927",
            ),
        )
        region_path = _write_gpkg(
            workspace / "region.gpkg",
            gpd.GeoDataFrame(
                {"name": ["left", "right"]},
                geometry=[
                    Polygon([(0.5, 0.5), (2.5, 0.5), (2.5, 2.5), (0.5, 2.5)]),
                    Polygon([(3.5, 3.5), (5.5, 3.5), (5.5, 5.5), (3.5, 5.5)]),
                ],
                crs="EPSG:2927",
            ),
        )

        tri = TriangleGrid(model_ws=str(workspace))
        tri.set_domain_file(domain_path)
        tri.add_region_file(region_path, max_area=2.0, label="zones", explode=True)

        preview = tri.preview_regions()
        points = tri.get_region_points()

        assert len(preview) == 2
        assert set(preview["label"]) == {"zones_0", "zones_1"}
        assert len(points) == 2
        assert all(tri.domain_geometry.covers(point) for point in points.geometry)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_triangle_grid_nested_regions_claim_unique_point_areas():
    workspace = _project_temp_dir("triangle_grid_nested_regions")
    try:
        tri = TriangleGrid(model_ws=str(workspace))
        tri.set_domain_rectangle(x_dist=12, y_dist=12, origin=(0, 0))
        tri.add_region_rectangle(
            origin=(2, 2),
            x_dist=8,
            y_dist=8,
            max_area=9.0,
            label="outer",
            priority=0,
        )
        tri.add_region_rectangle(
            origin=(4, 4),
            x_dist=4,
            y_dist=4,
            max_area=1.0,
            label="inner",
            priority=2,
        )

        preview = tri.preview_regions().set_index("label")
        points = tri.get_region_points().set_index("label")

        outer_geom = preview.loc["outer", "geometry"]
        inner_geom = preview.loc["inner", "geometry"]
        outer_point = points.loc["outer", "geometry"]
        inner_point = points.loc["inner", "geometry"]

        assert preview.loc["outer", "claim_area"] < outer_geom.area
        assert preview.loc["inner", "claim_area"] == inner_geom.area
        assert outer_geom.covers(outer_point)
        assert inner_geom.covers(inner_point)
        assert not inner_geom.covers(outer_point)
        assert outer_point.distance(inner_point) > 0
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_refined_triangle_voronoi_workflow_builds_substantially_detailed_grid():
    workspace = _project_temp_dir("triangle_refined_workflow_smoke")
    try:
        tri = TriangleGrid(model_ws=str(workspace))
        tri.set_domain_rectangle(x_dist=2000, y_dist=1500, origin=(0, 0))
        tri.add_region_rectangle(
            origin=(150, 150),
            x_dist=1700,
            y_dist=1200,
            max_area=60000,
            label="mid_refine",
        )
        tri.add_region_polygon(
            LineString([(120, 1100), (1880, 500)]).buffer(75),
            max_area=12000,
            label="stream_refine",
            priority=2,
            source="line",
        )
        tri.add_region_circle(
            center_coords=(700, 500),
            radius=180,
            max_area=5000,
            label="lake_refine",
            priority=3,
        )
        tri.build(verbose=False)

        import warnings

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            vor = VoronoiGridPlus(tri)
        preview = tri.preview_regions().set_index("label")
        flopy_geometry_warnings = [
            warning
            for warning in caught
            if "flopy\\utils\\geometry.py" in str(warning.filename).lower()
        ]

        assert vor.ncpl > 500
        assert np.isclose(vor.get_domain().area, 3_000_000.0)
        assert set(preview.index) == {"mid_refine", "stream_refine", "lake_refine"}
        assert preview.loc["lake_refine", "priority"] == 3
        assert preview.loc["stream_refine", "claim_area"] > 0
        assert flopy_geometry_warnings == []
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_triangle_grid_build_creates_missing_model_workspace():
    workspace = _project_temp_dir("triangle_missing_workspace")
    try:
        missing_ws = workspace / "nested" / "triangle_build"
        tri = TriangleGrid(model_ws=str(missing_ws))
        tri.set_domain_rectangle(x_dist=200, y_dist=150, origin=(0, 0))
        tri.add_region_circle(
            center_coords=(70, 55),
            radius=18,
            max_area=500,
            label="refine",
            priority=1,
        )

        tri.build(verbose=False)

        assert missing_ws.exists()
        assert any(path.name.startswith("_triangle.") for path in missing_ws.iterdir())
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_triangle_grid_rejects_fully_overlapping_regions_without_unique_claim_area():
    workspace = _project_temp_dir("triangle_grid_overlap_error")
    try:
        tri = TriangleGrid(model_ws=str(workspace))
        tri.set_domain_rectangle(x_dist=10, y_dist=10, origin=(0, 0))
        tri.add_region_rectangle(
            origin=(2, 2),
            x_dist=6,
            y_dist=6,
            max_area=4.0,
            label="first",
            priority=0,
        )
        tri.add_region_rectangle(
            origin=(2, 2),
            x_dist=6,
            y_dist=6,
            max_area=1.0,
            label="second",
            priority=0,
        )

        with pytest.raises(ValueError, match="unique interior area"):
            tri.prepare()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_boundary_package_builders_smoke():
    workspace = _project_temp_dir("boundary_packages_smoke")
    try:
        model, vor = _two_cell_model("boundary_smoke", workspace)

        boundary_gdf = gpd.GeoDataFrame(
            {
                "name": ["left"],
                "height": [1.5],
                "cond": [7.0],
                "layer": [1],
                "min_elev": [3.0],
                "elev": [6.0],
            },
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("name")

        drn = DRNFromVector(model=model, vor=vor, idomain=[1, 0])
        drn._gdf = boundary_gdf
        drn._intersections = pd.Series({"left": [0, 1]}, name="intersect")
        drn._edge_intersections = pd.Series({"left": [0, 1]}, name="intersect")

        drn_rows = drn.get_drn_stress_period_data(
            cells=[0, 1],
            bottom_addition=1.0,
            conductance=[5.0, 9.0],
            disMf="disv",
            layer=0,
            region_name="manual_drn_group",
            region_tags=["drn", "manual"],
            region_metadata={"source": "unit_test"},
        )
        assert drn_rows == [[(0, 0), 3.0, 5.0]]
        assert model.get_region_cells("manual_drn_group") == [0]
        assert model.get_region("manual_drn_group").metadata["source"] == "unit_test"

        drn_poly = drn.from_vector(
            edges_only=True,
            register_regions=True,
            region_name_prefix="drn_group",
            combined_region_name="all_drn_groups",
            region_tags=["drn", "polygon"],
            overwrite_regions=True,
        )
        assert drn_poly[0] == [[(0, 0), 3.5, 7.0]]
        assert model.get_region_cells("drn_group_left") == [0]
        assert set(model.get_region_cells("all_drn_groups")) == {0}

        ghb = GHBFromVector(model=model, vor=vor, idomain=[1, 0])
        ghb._gdf = boundary_gdf
        ghb._edge_intersections = pd.Series({"left": [0, 1]}, name="intersect")

        ghb_poly = ghb.from_vector(
            elev_reference={"left": [7.0, 8.0]},
            reference_offset=0.5,
            register_regions=True,
            region_name_prefix="ghb_group",
            combined_region_name="all_ghb_groups",
            region_tags=["ghb", "polygon"],
            overwrite_regions=True,
        )
        assert ghb_poly[0] == [[(0, 0), 7.5, 7.0]]
        assert ghb_poly[1] == [[(0, 0), 8.5, 7.0]]
        assert model.get_region_cells("ghb_group_left") == [0]
        assert set(model.get_region_cells("all_ghb_groups")) == {0}
        assert model.get_region("ghb_group_left").metadata["uses_elev_reference"] is True

        ghb_updated = ghb.add_to_dict(existing_dict={0: [[(0, 0), 1.0, 2.0]]}, dict_to_add={0: [[(0, 0), 9.0, 3.0]]})
        assert ghb_updated[0] == [[(0, 0), 9.0, 3.0]]

        recharge_gdf = gpd.GeoDataFrame(
            {"zone": ["zone_a"], "rch_a": [1.0], "rch_b": [0.6]},
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("zone")

        recharge = RCHFromVector(
            model=model,
            vor=vor,
            uid="zone",
            rch_fields=["rch_a", "rch_b"],
            rch_fields_to_pers=[0, 1],
            background_rch=0.05,
            grid_type="disv",
            limit_to_k33=True,
            limit_to_k33_by=0.5,
        )
        recharge._gdf = recharge_gdf
        recharge._rch_scale = pd.Series({"zone_a": 1.0})
        recharge._intersections_no_duplicates = pd.DataFrame(
            {"no_dup": [[0]]},
            index=pd.Index(["zone_a"], name="zone"),
        )

        rch_dict = recharge.from_vector(
            register_regions=True,
            region_name_prefix="rch_zone",
            combined_region_name="all_rch",
            region_tags=["rch", "unit"],
            overwrite_regions=True,
        )
        assert rch_dict[0] == [[(0, 0), 0.125], [(0, 1), 0.05]]
        assert rch_dict[1] == [[(0, 0), 0.125], [(0, 1), 0.05]]
        assert model.get_region_cells("rch_zone_zone_a") == [0]
        assert set(model.get_region_cells("all_rch")) == {0}

        merged_rch = RCHFromVector.add_to_rch_dict(
            rch_dict,
            {0: [[(0, 0), 0.25]], 1: [[(0, 1), 0.02]]},
            replace=False,
        )
        assert merged_rch[0] == [[(0, 0), 0.375], [(0, 1), 0.05]]
        assert merged_rch[1] == [[(0, 0), 0.125], [(0, 1), 0.07]]

        uzf = UZFBuilder(
            context=_model_context(model, vor),
            nper=model.nper,
            cells=[(0, 0), (0, 1)],
            vks=1.0,
            thtr=0.1,
            thts=0.3,
            thti=0.2,
            pet={1: [0.01, 0.02]},
            finf=_uzf_finf_from_rch(rch_dict, [(0, 0), (0, 1)], model.nper),
        )

        assert uzf.finf[0] == [0.125, 0.05]
        assert uzf.finf[1] == [0.125, 0.05]
        assert len(uzf.packagedata) == 2
        assert uzf.packagedata[0][1] == (0, 0)
        assert uzf.perioddata[0][0][1] == 0.125
        assert uzf.perioddata[1][1][2] == 0.02
        uzf.build().build(model.gwf)
        model.add_region_from_cells(
            "uzf_zone", uzf.uzf_cells, category="boundary", package="uzf",
            tags=["uzf", "unit"], overwrite=True,
        )
        assert set(model.get_region_cells("uzf_zone")) == {0, 1}
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_boundaries_shared_base_exposes_polygon_iteration_only():
    assert hasattr(Boundaries, "iter_polygon_boundary_features")
    assert not hasattr(Boundaries, "get_rch_dict")
    assert not hasattr(Boundaries, "get_ghb_from_shp")
    assert not hasattr(Boundaries, "get_chd_from_shp")
    assert not hasattr(Boundaries, "get_k_from_shp")


def test_boundaries_polygon_feature_iteration_filters_inactive_cells():
    workspace = _project_temp_dir("boundary_feature_iteration")
    try:
        model, vor = _two_cell_model("boundary_feature_iteration", workspace)
        boundary_gdf = gpd.GeoDataFrame(
            {
                "name": ["left"],
                "height": [1.5],
                "cond": [7.0],
                "layer": [1],
                "min_elev": [3.0],
            },
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("name")

        boundaries = Boundaries(model=model, vor=vor, uid="name", idomain=[1, 0])
        boundaries._gdf = boundary_gdf
        boundaries._intersections = pd.Series({"left": [0, 1]}, name="intersect")

        features = list(boundaries.iter_polygon_boundary_features(name_field="name"))

        assert len(features) == 1
        name, row, active_cells = features[0]
        assert name == "left"
        assert row["height"] == 1.5
        assert active_cells == [0]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_chd_from_vector_builds_stress_period_data_and_regions():
    workspace = _project_temp_dir("chd_from_vector")
    try:
        model, vor = _two_cell_model("chd_from_vector", workspace)
        chd_gdf = gpd.GeoDataFrame(
            {
                "name": ["left"],
                "elev": [10.0],
                "layer": [1],
            },
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("name")

        chd = CHDFromVector(model=model, vor=vor, uid="name", idomain=[1, 0])
        chd._gdf = chd_gdf
        chd._intersections = pd.Series({"left": [0, 1]}, name="intersect")

        chd_dict = chd.from_vector(
            head_reference={"left": [11.0, 12.0]},
            reference_offset=0.5,
            register_regions=True,
            region_name_prefix="chd_group",
            combined_region_name="all_chd",
            region_tags=["chd", "polygon"],
            overwrite_regions=True,
        )

        assert chd_dict[0] == [[(0, 0), 11.5]]
        assert chd_dict[1] == [[(0, 0), 12.5]]
        assert model.get_region_cells("chd_group_left") == [0]
        assert set(model.get_region_cells("all_chd")) == {0}
        assert model.get_region("chd_group_left").metadata["uses_head_reference"] is True
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_k_from_vector_assigns_all_matching_cells_and_layer_defaults():
    workspace = _project_temp_dir("k_from_vector")
    try:
        model, vor = _two_cell_model("k_from_vector", workspace)
        k_gdf = gpd.GeoDataFrame(
            {
                "name": ["both_cells"],
                "k": [3.25],
                "layer": [1],
            },
            geometry=[Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("name")

        k_builder = KFromVector(model=model, vor=vor, uid="name")
        k_builder._gdf = k_gdf
        k_builder._intersections = pd.Series({"both_cells": [0, 1]}, name="intersect")

        k_df = k_builder.to_frame(defaults=[0.5])
        k_array = k_builder.to_array(defaults=[0.5])

        assert k_df.loc[(0, 0), "k"] == pytest.approx(3.25)
        assert k_df.loc[(0, 1), "k"] == pytest.approx(3.25)
        assert k_array.shape == (1, 2)
        assert k_array[0, 0] == pytest.approx(3.25)
        assert k_array[0, 1] == pytest.approx(3.25)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_vector_builder_aliases_are_consistent_across_packages():
    workspace = _project_temp_dir("vector_builder_aliases")
    try:
        model, vor = _two_cell_model("vector_builder_aliases", workspace)

        boundary_gdf = gpd.GeoDataFrame(
            {
                "name": ["left"],
                "height": [1.5],
                "cond": [7.0],
                "layer": [1],
                "min_elev": [3.0],
                "elev": [6.0],
            },
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("name")

        drn = DRN(model=model, vor=vor, idomain=[1, 0])
        drn._gdf = boundary_gdf
        drn._intersections = pd.Series({"left": [0, 1]}, name="intersect")
        assert drn.from_vector() == drn.from_polygons()
        assert drn.get_drn_from_poly() == drn.from_polygons()

        ghb = GHB(model=model, vor=vor, idomain=[1, 0])
        ghb._gdf = boundary_gdf
        ghb._edge_intersections = pd.Series({"left": [0, 1]}, name="intersect")
        assert ghb.from_vector() == ghb.from_polygons()
        assert ghb.get_from_poly() == ghb.from_polygons()

        chd = CHDFromVector(model=model, vor=vor, idomain=[1, 0])
        chd._gdf = boundary_gdf.loc[:, ["elev", "layer", "geometry"]].copy()
        chd._intersections = pd.Series({"left": [0, 1]}, name="intersect")
        assert chd.from_vector() == chd.from_polygons()
        assert chd.get_from_poly() == chd.from_polygons()

        recharge_gdf = gpd.GeoDataFrame(
            {"zone": ["zone_a"], "rch_a": [1.0], "rch_b": [0.6]},
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("zone")
        recharge = RechargeFromShp(
            model=model,
            vor=vor,
            uid="zone",
            rch_fields=["rch_a", "rch_b"],
            rch_fields_to_pers=[0, 1],
            background_rch=0.05,
            grid_type="disv",
            limit_to_k33=False,
        )
        recharge._gdf = recharge_gdf
        recharge._rch_scale = pd.Series({"zone_a": 1.0})
        recharge._intersections_no_duplicates = pd.DataFrame(
            {"no_dup": [[0]]},
            index=pd.Index(["zone_a"], name="zone"),
        )
        assert recharge.from_vector() == recharge.from_polygons()
        assert recharge.get_rch() == recharge.from_polygons()

        k_gdf = gpd.GeoDataFrame(
            {"name": ["both_cells"], "k": [3.25], "layer": [1]},
            geometry=[Polygon([(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0)])],
            crs=vor.crs,
        ).set_index("name")
        k_builder = KFromVector(model=model, vor=vor, uid="name")
        k_builder._gdf = k_gdf
        k_builder._intersections = pd.Series({"both_cells": [0, 1]}, name="intersect")
        np.testing.assert_allclose(
            k_builder.from_vector(),
            k_builder.from_polygons(),
        )
        pd.testing.assert_frame_equal(
            k_builder.from_vector(return_array=False),
            k_builder.from_polygons(return_array=False),
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_geopackage_source_is_the_public_vector_builder_api():
    assert myflopy.GeoPackageSource is not None
    for legacy_name in (
        "DRNFromVector",
        "GHBFromVector",
        "CHDFromVector",
        "RCHFromVector",
        "KFromVector",
    ):
        assert not hasattr(myflopy, legacy_name)
