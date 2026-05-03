from __future__ import annotations

import json
import os
import subprocess
import shutil
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import flopy
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, Polygon

plt.rcParams["figure.max_open_warning"] = 0

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from simple_modflow import (  # noqa: E402
    DiscoveredRun,
    LoadedMf6Run,
    ModelGroup,
    ModelSpec,
    PackageArtifact,
    PackageCompatibilityError,
    ProjectCatalog,
    RunComparison,
    RunExplorer,
    RunLoader,
    RunRecord,
    RunSpec,
    discover_existing_runs,
    explore_runs,
    import_run_archive,
    load_mf6_run,
    patch_simulation_plot,
    summarize_discovered_runs,
)
from simple_modflow.modflow.mf6.grid.triangle import TriangleGrid  # noqa: E402
from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from simple_modflow.modflow.mf6.lakes import (  # noqa: E402
    LakeConnectionData,
    LakePackageData,
    LakePeriodData,
)
from simple_modflow.modflow.mf6.sfr import SFR  # noqa: E402
from simple_modflow.modflow.mf6.budget import DRNBudget  # noqa: E402
from simple_modflow.modflow.mf6.budget import LakStage, SFRStage  # noqa: E402
from simple_modflow.modflow.mf6 import budget as budget_module  # noqa: E402
from simple_modflow.modflow.mf6.package_explorer import build_lak_q_map_payload  # noqa: E402
from simple_modflow.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from simple_modflow.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisuGrid,
    DisvGrid,
    TemporalDiscretization,
)
from simple_modflow.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    Drains,
    GHB,
    InitialConditions,
    KFlow,
    LAK as LAKPackage,
    MVR as MVRPackage,
    OutputControl,
    Recharge,
    Storage,
    UZF as UZFPackage,
)
from simple_modflow.project.manifest_io import load_run_record  # noqa: E402
from simple_modflow.modflow.mf6.uzf import UZFPackageData  # noqa: E402


def _project_temp_dir(name: str) -> Path:
    root = ROOT / ".pytest-work" / name
    if root.exists():
        shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    return root


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


def _three_cell_vor_clockwise():
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [3.0, 0.0],
            [3.0, 1.0],
        ],
        dtype=float,
    )
    iverts = [[0, 3, 2, 1], [1, 2, 5, 4], [4, 5, 7, 6]]
    xcyc = np.array([[0.5, 0.5], [1.5, 0.5], [2.5, 0.5]], dtype=float)
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _four_cell_vor_clockwise():
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.0, 2.0],
            [2.0, 2.0],
        ],
        dtype=float,
    )
    iverts = [[0, 3, 2, 1], [1, 2, 5, 4], [3, 6, 7, 2], [2, 7, 8, 5]]
    xcyc = np.array(
        [
            [0.5, 0.5],
            [1.5, 0.5],
            [0.5, 1.5],
            [1.5, 1.5],
        ],
        dtype=float,
    )
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: [12.0, 12.0, 12.0, 12.0],
            1: [0.0, 0.0, 0.0, 0.0],
        },
        geometry="geometry",
        crs=vor.crs,
    )
    return vor


def _write_gpkg(path: Path, gdf: gpd.GeoDataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _build_small_optimized_voronoi_model(workspace: Path):
    """Build a modest refined Voronoi grid plus geometry inputs for artifact tests."""

    domain_geom = Polygon([(0.0, 0.0), (520.0, 0.0), (520.0, 420.0), (0.0, 420.0)])
    lake_geom = Point(180.0, 160.0).buffer(52.0)
    stream_line = LineString([(35.0, 355.0), (210.0, 265.0), (485.0, 120.0)])

    tri = TriangleGrid(model_ws=str(workspace / "triangle_build"))
    tri.set_domain_polygon(domain_geom)
    tri.add_region_polygon(
        Polygon([(45.0, 45.0), (475.0, 45.0), (475.0, 375.0), (45.0, 375.0)]),
        max_area=18000,
        label="mid_refine",
    )
    tri.add_region_polygon(
        stream_line.buffer(28.0),
        max_area=3500,
        label="stream_refine",
        priority=2,
        source="line",
    )
    tri.add_region_polygon(
        lake_geom,
        max_area=900,
        label="lake_refine",
        priority=3,
        source="circle",
    )
    mesh_report = tri.build_mesh(profile="balanced", verbose=False)

    vor = VoronoiGridPlus(tri)
    top = 118.0 - (0.02 * np.asarray(vor.centroids_x)) + (0.009 * np.asarray(vor.centroids_y))
    bottom = top - 34.0
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: top,
            1: bottom,
        },
        geometry="geometry",
        crs=vor.crs,
    )

    geometries = {
        "domain": domain_geom,
        "lake": lake_geom,
        "stream": stream_line,
        "east_boundary": Polygon([(470.0, 0.0), (520.0, 0.0), (520.0, 420.0), (470.0, 420.0)]),
    }
    return vor, mesh_report, geometries


def _build_and_run_two_cell_model(record: RunRecord, heads: tuple[float, float] = (10.0, 9.0)) -> SimulationBase:
    vor = _two_cell_vor_clockwise()
    model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
    DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=list(heads))
    KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[(0, 0), heads[0]], [(0, 1), heads[1]]]})
    success, _ = model.run_simulation()
    assert success is True
    return model


def _build_and_run_two_cell_rch_uzf_model(
    record: RunRecord,
    *,
    heads: tuple[float, float] = (10.0, 9.0),
    recharge: tuple[float, float] = (0.001, 0.00075),
    finf: tuple[float, float] = (0.0005, 0.00025),
) -> SimulationBase:
    vor = _two_cell_vor_clockwise()
    model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
    DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=list(heads))
    KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[(0, 0), heads[0]], [(0, 1), heads[1]]]} )
    Recharge(
        model=model,
        rch_dict={0: [[(0, 0), recharge[0]], [(0, 1), recharge[1]]]},
    )
    UZFPackageData(
        model=model,
        vor=vor,
        uzf_cells=[(0, 0), (0, 1)],
        vks=1.0,
        thtr=0.1,
        thts=0.3,
        thti=0.2,
        finf={0: list(finf)},
        add_uzf=True,
    )
    success, _ = model.run_simulation()
    assert success is True
    return model


def _build_and_run_four_cell_lake_model(
    record: RunRecord,
    workspace: Path,
    *,
    starting_stage: float = 10.5,
) -> SimulationBase:
    vor = _four_cell_vor_clockwise()
    model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
    DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=[11.9, 10.1, 11.7, 9.9])
    KFlow(model=model, k=[1.0] * 4, save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[(0, 1), 12.0], [(0, 3), 10.0]]})

    lake_path = _write_gpkg(
        workspace / f"{record.run_id}_lake.gpkg",
        gpd.GeoDataFrame(
            {"name": ["lake_0"]},
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.95), (0.0, 1.95)])],
            crs=vor.crs,
        ),
    )
    lake_connections = LakeConnectionData(
        model=model,
        vor=vor,
        paths=[lake_path],
        bed_leakance=1.0,
        horizontal_connections={0: [11.5, 9.0]},
        use_reconciled_surfaces=False,
        only_vertical=False,
    )
    connectiondata = lake_connections.connection_data
    lake_packagedata = LakePackageData(nlakes=1, starting_stage=[starting_stage], connectiondata=connectiondata)
    lake_perioddata = LakePeriodData(
        model=model,
        lake_ids=[0],
        lake_stages=[starting_stage],
        status=["ACTIVE"],
    )
    LAKPackage(
        model=model,
        nlakes=1,
        noutlets=0,
        ntables=0,
        packagedata=lake_packagedata.packagedata,
        connectiondata=connectiondata,
        perioddata=lake_perioddata.perioddata,
        mover=False,
    )
    success, _ = model.run_simulation()
    assert success is True
    return model


def _build_and_run_four_cell_lak_sfr_model(
    record: RunRecord,
    workspace: Path,
    *,
    lake_stage: float = 10.5,
    sfr_inflow: float = 0.5,
) -> SimulationBase:
    vor = _four_cell_vor_clockwise()
    model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
    DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=[11.9, 10.1, 11.7, 9.9])
    KFlow(model=model, k=[1.0] * 4, save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[(0, 1), 12.0], [(0, 3), 10.0]]})

    lake_path = _write_gpkg(
        workspace / f"{record.run_id}_lake.gpkg",
        gpd.GeoDataFrame(
            {"name": ["lake_0"]},
            geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 1.95), (0.0, 1.95)])],
            crs=vor.crs,
        ),
    )
    stream_path = _write_gpkg(
        workspace / f"{record.run_id}_stream.gpkg",
        gpd.GeoDataFrame(
            {"name": ["stream_0"]},
            geometry=[LineString([(1.2, 1.8), (1.8, 1.2)])],
            crs=vor.crs,
        ),
    )

    lake_connections = LakeConnectionData(
        model=model,
        vor=vor,
        paths=[lake_path],
        bed_leakance=1.0,
        horizontal_connections={0: [11.5, 9.0]},
        use_reconciled_surfaces=False,
        only_vertical=False,
    )
    connectiondata = lake_connections.connection_data
    lake_packagedata = LakePackageData(nlakes=1, starting_stage=[lake_stage], connectiondata=connectiondata)
    lake_perioddata = LakePeriodData(model=model, lake_ids=[0], lake_stages=[lake_stage], status=["ACTIVE"])
    LAKPackage(
        model=model,
        nlakes=1,
        noutlets=0,
        ntables=0,
        packagedata=lake_packagedata.packagedata,
        connectiondata=connectiondata,
        perioddata=lake_perioddata.perioddata,
        mover=False,
    )

    SFR(
        model=model,
        vor=vor,
        stream_paths=[stream_path],
        inflows={0: [(0, sfr_inflow)]},
        widths=5.0,
        gradients=0.001,
        mannings=0.03,
        streambed_k=1.0,
        streambed_thickness=1.0,
        mover=False,
        add_sfr=True,
    )

    success, _ = model.run_simulation()
    assert success is True
    return model


def _build_and_run_four_cell_standard_budget_model(
    record: RunRecord,
    *,
    chd_heads: tuple[float, float] = (11.0, 10.5),
    drn_elev: float = 10.8,
    drn_cond: float = 0.4,
    ghb_head: float = 9.8,
    ghb_cond: float = 0.35,
    recharge: tuple[float, float, float, float] = (0.0008, 0.0008, 0.0006, 0.0006),
) -> SimulationBase:
    vor = _four_cell_vor_clockwise()
    model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
    DisvGrid(vor=vor, model=model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=[11.4, 11.0, 10.2, 10.0])
    KFlow(model=model, k=[1.0] * 4, save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[(0, 1), chd_heads[0]], [(0, 3), chd_heads[1]]]})
    Drains(model=model, stress_period_data={0: [[(0, 0), drn_elev, drn_cond]]})
    GHB(model=model, stress_period_data={0: [[(0, 2), ghb_head, ghb_cond]]})
    Recharge(
        model=model,
        vor=vor,
        rch_dict={0: [[(0, cell), recharge[cell]] for cell in range(vor.ncpl)]},
    )
    success, _ = model.run_simulation()
    assert success is True
    return model


def _build_legacy_two_cell_run(workspace: Path, run_name: str, heads: tuple[float, float] = (10.0, 9.0)) -> RunRecord:
    record = RunRecord(
        run_id=run_name,
        model_spec="legacy",
        workspace=workspace,
        status="completed",
    )
    _build_and_run_two_cell_model(record, heads=heads)
    return record


def _build_and_run_two_cell_disu_model(
    record: RunRecord,
    heads: tuple[float, float] = (10.0, 9.0),
) -> SimulationBase:
    vor = _two_cell_vor_clockwise()
    model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
    DisuGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[0.0, 0.0])
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=list(heads))
    KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)
    CHD(model=model, stress_period_data={0: [[0, heads[0]], [1, heads[1]]]})
    success, _ = model.run_simulation()
    assert success is True
    return model


def test_project_catalog_public_imports_available():
    assert DiscoveredRun is not None
    assert LoadedMf6Run is not None
    assert ProjectCatalog is not None
    assert ModelSpec is not None
    assert RunLoader is not None
    assert RunComparison is not None
    assert RunExplorer is not None
    assert RunSpec is not None
    assert RunRecord is not None
    assert discover_existing_runs is not None
    assert explore_runs is not None
    assert import_run_archive is not None
    assert load_mf6_run is not None
    assert patch_simulation_plot is not None
    assert summarize_discovered_runs is not None


def test_project_catalog_registers_specs_and_runs_roundtrip():
    workspace = _project_temp_dir("project_catalog_roundtrip")
    try:
        catalog = ProjectCatalog(workspace / "demo_project", name="demo_project")
        spec = ModelSpec(
            name="baseline_model",
            description="baseline setup",
            grid_ref="grid_v1",
            tags=["baseline"],
            default_packages={"drn": "drn_v1", "rch": "rch_v1"},
            default_regions=["all_drains"],
        )
        catalog.register_model_spec(spec)

        run_spec = RunSpec(
            run_id="baseline_run",
            model_spec="baseline_model",
            scenario="baseline",
            package_versions={"drn": "drn_v2", "rch": "rch_v3"},
            tags=["calibration"],
        )
        record = catalog.create_run(run_spec)

        assert record.workspace == catalog.runs_dir / "baseline_run"
        assert record.manifest_path.exists()
        assert record.paths["workspace"] == "."
        assert record.paths["simulation_name_file"] == "mfsim.nam"
        assert record.paths["model_name_file"] == "baseline_run.nam"

        loaded_spec = catalog.load_model_spec("baseline_model")
        loaded_run = catalog.load_run("baseline_run")

        assert loaded_spec.default_packages["drn"] == "drn_v1"
        assert loaded_run.package_versions["rch"] == "rch_v3"
        assert loaded_run.workspace == record.workspace
        assert loaded_run.simulation_kwargs()["name"] == "baseline_run"
        assert loaded_run.simulation_kwargs()["mf_folder_path"] == record.workspace.parent

        specs = catalog.list_model_specs()
        runs = catalog.list_runs()

        assert specs["name"].tolist() == ["baseline_model"]
        assert runs["run_id"].tolist() == ["baseline_run"]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_run_manifest_stays_portable_when_workspace_is_copied():
    workspace = _project_temp_dir("project_catalog_copyable_run")
    try:
        catalog = ProjectCatalog(workspace / "portable_project", name="portable_project")
        catalog.register_model_spec(ModelSpec(name="portable_model"))
        record = catalog.create_run(RunSpec(run_id="portable_run", model_spec="portable_model"))

        copied_workspace = workspace / "copied_run"
        shutil.copytree(record.workspace, copied_workspace)
        copied_record = load_run_record(copied_workspace / "run.toml")

        assert copied_record.workspace == copied_workspace
        assert copied_record.get_path("simulation_name_file") == copied_workspace / "mfsim.nam"
        assert copied_record.get_path("run_manifest") == copied_workspace / "run.toml"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_run_workspace_can_hold_complete_mf6_run():
    workspace = _project_temp_dir("project_catalog_mf6_run")
    try:
        catalog = ProjectCatalog(workspace / "run_project", name="run_project")
        catalog.register_model_spec(
            ModelSpec(
                name="tiny_model",
                grid_ref="two_cell_grid",
                default_packages={"oc": "oc_v1", "chd": "chd_v1"},
            )
        )
        record = catalog.create_run(
            RunSpec(
                run_id="tiny_run",
                model_spec="tiny_model",
                scenario="smoke",
                package_versions={"oc": "oc_v1", "chd": "chd_v1"},
            )
        )

        model = _build_and_run_two_cell_model(record, heads=(10.0, 9.0))

        assert record.workspace.exists()
        expected_files = {
            "mfsim.nam",
            "mfsim.lst",
            "run.toml",
            "tiny_run.disv",
            "tiny_run.hds",
            "tiny_run.nam",
            "tiny_run.oc",
            "tiny_run.tdis",
        }
        assert expected_files.issubset({path.name for path in record.workspace.iterdir()})

        reloaded = catalog.load_run("tiny_run")
        assert reloaded.workspace == record.workspace
        assert reloaded.get_path("heads_file") == record.workspace / "tiny_run.hds"
        assert np.allclose(model.hds.get_data(kstpkper=(0, 0)).squeeze(), [10.0, 9.0])
        assert model.package_names == ["CHD", "DISV", "IC", "NPF", "OC", "STO"]
        assert model.grid_type == "disv"
        assert model.summary().set_index("name").loc["tiny_run", "source"] == "simulation"
        assert model.package_summary().set_index("package").loc["CHD", "package_type"] == "chd"
        assert bool(model.output_summary().iloc[0]["has_heads"]) is True
        assert model.grid_summary().iloc[0]["ncpl"] == 2
        assert model.result_summary().iloc[0]["head_mean"] == 9.5
        assert set(model.file_summary()["category"]) == {"input", "output"}
        assert {path.name for path in model.list_input_files()} >= {"mfsim.nam", "tiny_run.disv", "tiny_run.nam"}
        assert {path.name for path in model.list_output_files()} >= {"mfsim.lst", "tiny_run.hds", "tiny_run.cbc"}
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_supports_multiple_realistic_runs_and_reopen():
    workspace = _project_temp_dir("project_catalog_multiple_runs")
    try:
        catalog_root = workspace / "multi_run_project"
        catalog = ProjectCatalog(catalog_root, name="multi_run_project")
        catalog.register_model_spec(
            ModelSpec(
                name="tiny_model",
                description="two-cell steady-state model",
                grid_ref="two_cell_grid",
                default_packages={"oc": "oc_v1", "chd": "chd_v1"},
            )
        )

        baseline = catalog.create_run(
            RunSpec(
                run_id="baseline_run",
                model_spec="tiny_model",
                scenario="baseline",
                package_versions={"oc": "oc_v1", "chd": "chd_v1"},
                tags=["baseline"],
            )
        )
        variant = catalog.create_run(
            RunSpec(
                run_id="variant_run",
                model_spec="tiny_model",
                scenario="lower_right_head",
                package_versions={"oc": "oc_v1", "chd": "chd_v2"},
                tags=["scenario"],
            )
        )

        baseline_model = _build_and_run_two_cell_model(baseline, heads=(10.0, 9.0))
        variant_model = _build_and_run_two_cell_model(variant, heads=(10.0, 8.5))

        reopened = ProjectCatalog(catalog_root, create=False)
        runs = reopened.list_runs().sort_values("run_id").reset_index(drop=True)
        discovered = reopened.discover_runs()

        assert runs["run_id"].tolist() == ["baseline_run", "variant_run"]
        assert {run.run_id for run in discovered} == {"baseline_run", "variant_run"}
        assert reopened.load_run("baseline_run").package_versions["chd"] == "chd_v1"
        assert reopened.load_run("variant_run").package_versions["chd"] == "chd_v2"

        baseline_heads = flopy.utils.HeadFile(baseline.get_path("heads_file")).get_data(kstpkper=(0, 0)).squeeze()
        variant_heads = flopy.utils.HeadFile(variant.get_path("heads_file")).get_data(kstpkper=(0, 0)).squeeze()

        assert np.allclose(baseline_heads, [10.0, 9.0])
        assert np.allclose(variant_heads, [10.0, 8.5])
        assert not np.allclose(baseline_heads, variant_heads)
        assert np.allclose(baseline_model.hds.get_data(kstpkper=(0, 0)).squeeze(), baseline_heads)
        assert np.allclose(variant_model.hds.get_data(kstpkper=(0, 0)).squeeze(), variant_heads)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_copied_run_directory_can_be_rerun_with_plain_mf6():
    workspace = _project_temp_dir("project_catalog_plain_mf6_rerun")
    try:
        catalog = ProjectCatalog(workspace / "portable_project", name="portable_project")
        catalog.register_model_spec(ModelSpec(name="portable_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(
                run_id="portable_run",
                model_spec="portable_model",
                package_versions={"chd": "chd_v1"},
            )
        )

        _build_and_run_two_cell_model(record, heads=(10.0, 9.0))

        copied_workspace = workspace / "portable_run_copy"
        shutil.copytree(record.workspace, copied_workspace)

        for file_name in ["mfsim.lst", "portable_run.hds", "portable_run.cbc"]:
            copied_file = copied_workspace / file_name
            if copied_file.exists():
                copied_file.unlink()

        completed = subprocess.run(
            ["mf6"],
            cwd=copied_workspace,
            capture_output=True,
            text=True,
            check=False,
        )

        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert (copied_workspace / "mfsim.lst").exists()
        assert (copied_workspace / "portable_run.hds").exists()

        copied_heads = flopy.utils.HeadFile(copied_workspace / "portable_run.hds").get_data(kstpkper=(0, 0)).squeeze()
        assert np.allclose(copied_heads, [10.0, 9.0])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_copied_project_catalog_rediscovers_runs_in_new_location():
    workspace = _project_temp_dir("project_catalog_copy_project")
    try:
        original_root = workspace / "original_project"
        catalog = ProjectCatalog(original_root, name="original_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        baseline = catalog.create_run(RunSpec(run_id="run_a", model_spec="tiny_model", package_versions={"chd": "v1"}))
        variant = catalog.create_run(RunSpec(run_id="run_b", model_spec="tiny_model", package_versions={"chd": "v2"}))

        _build_and_run_two_cell_model(baseline, heads=(10.0, 9.0))
        _build_and_run_two_cell_model(variant, heads=(10.0, 8.5))

        copied_root = workspace / "copied_project"
        shutil.copytree(original_root, copied_root)

        copied_catalog = ProjectCatalog(copied_root, create=False)
        copied_runs = {run.run_id: run for run in copied_catalog.discover_runs()}

        assert set(copied_runs) == {"run_a", "run_b"}
        assert copied_runs["run_a"].workspace == copied_root / "runs" / "run_a"
        assert copied_runs["run_b"].workspace == copied_root / "runs" / "run_b"

        copied_run_a_heads = flopy.utils.HeadFile(copied_runs["run_a"].get_path("heads_file")).get_data(kstpkper=(0, 0)).squeeze()
        copied_run_b_heads = flopy.utils.HeadFile(copied_runs["run_b"].get_path("heads_file")).get_data(kstpkper=(0, 0)).squeeze()

        assert np.allclose(copied_run_a_heads, [10.0, 9.0])
        assert np.allclose(copied_run_b_heads, [10.0, 8.5])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_discover_existing_runs_supports_plain_mf6_workspaces_without_model_pickle():
    workspace = _project_temp_dir("project_discovery_plain_mf6")
    try:
        run_workspace = workspace / "plain_run"
        record = RunRecord(
            run_id="plainrun",
            model_spec="plain",
            workspace=run_workspace,
            status="completed",
        )
        _build_and_run_two_cell_model(record, heads=(10.0, 9.0))

        for model_object in workspace.rglob("*.model"):
            model_object.unlink()

        discovered = discover_existing_runs(workspace)

        assert len(discovered) == 1
        assert discovered[0].run_id == "plainrun"
        assert discovered[0].model_file is None
        assert discovered[0].metadata["has_heads"] is True
        assert discovered[0].metadata["grid_type"] == "disv"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_run_explorer_browses_filters_and_opens_runs_naturally():
    workspace = _project_temp_dir("project_run_explorer")
    try:
        first = _build_legacy_two_cell_run(workspace / "cumb_alpha", "cumb_alpha", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(workspace / "cumb_beta", "cumb_beta", heads=(10.0, 8.5))
        third = _build_legacy_two_cell_run(workspace / "lkpt_gamma", "lkpt_gamma", heads=(10.0, 8.0))

        explorer = explore_runs(workspace, model_spec="archive_runs")

        assert len(explorer) == 3
        assert set(explorer.run_ids) == {"cumb_alpha", "cumb_beta", "lkpt_gamma"}

        summary = explorer.summary().sort_values("run_id").reset_index(drop=True)
        assert summary["run_id"].tolist() == ["cumb_alpha", "cumb_beta", "lkpt_gamma"]
        assert summary["has_heads"].tolist() == [True, True, True]

        families = explorer.families()
        family_counts = dict(zip(families["family"], families["count"], strict=False))
        assert family_counts["cumb"] == 2
        assert family_counts["lkpt"] == 1

        cumb_only = explorer.filter(family="cumb")
        assert set(cumb_only.run_ids) == {"cumb_alpha", "cumb_beta"}

        searched = explorer.filter(text="gamma")
        assert searched.run_ids == ["lkpt_gamma"]

        loaded = explorer.open("cumb_beta")
        heads = loaded.hds.get_data(kstpkper=(0, 0)).squeeze()
        assert np.allclose(heads, [10.0, 8.5])

        comparison = explorer.compare.compare_heads("cumb_alpha", "cumb_beta")
        assert "head_diff" in comparison.columns
        assert np.isclose(comparison["head_diff"].iloc[-1], -0.5)

        assert explorer["lkpt_gamma"].workspace == third.workspace
        assert explorer.load_run("cumb_alpha").workspace == first.workspace
        assert explorer.load_run("cumb_beta").workspace == second.workspace
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_duplicate_registration_and_overwrite_rules():
    workspace = _project_temp_dir("project_catalog_overwrite_rules")
    try:
        catalog = ProjectCatalog(workspace / "overwrite_project", name="overwrite_project")
        spec = ModelSpec(name="tiny_model", description="first version", grid_ref="grid_a")
        catalog.register_model_spec(spec)

        with pytest.raises(FileExistsError, match="Model spec already exists"):
            catalog.register_model_spec(ModelSpec(name="tiny_model", description="second version"))

        catalog.register_model_spec(
            ModelSpec(name="tiny_model", description="second version", grid_ref="grid_b"),
            overwrite=True,
        )
        assert catalog.load_model_spec("tiny_model").grid_ref == "grid_b"

        first_record = catalog.create_run(RunSpec(run_id="run_a", model_spec="tiny_model"))
        assert first_record.manifest_path.exists()

        with pytest.raises(FileExistsError, match="Run workspace already exists"):
            catalog.create_run(RunSpec(run_id="run_a", model_spec="tiny_model"))

        replacement = catalog.create_run(
            RunSpec(run_id="run_a", model_spec="tiny_model", package_versions={"chd": "v2"}),
            overwrite=True,
        )
        assert replacement.package_versions["chd"] == "v2"
        assert catalog.load_run("run_a").package_versions["chd"] == "v2"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_supports_custom_project_layout_and_relative_run_workspace():
    workspace = _project_temp_dir("project_catalog_custom_layout")
    try:
        catalog = ProjectCatalog(
            workspace / "custom_project",
            name="custom_project",
            runs_dir_name="model_runs",
            model_specs_dir_name="spec_library",
        )
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        record = catalog.create_run(
            RunSpec(
                run_id="custom_run",
                model_spec="tiny_model",
                workspace=Path("manual_runs") / "custom_run",
                package_versions={"chd": "v1"},
            )
        )

        reopened = ProjectCatalog(workspace / "custom_project", create=False)

        assert catalog.runs_dir == workspace / "custom_project" / "model_runs"
        assert catalog.model_specs_dir == workspace / "custom_project" / "spec_library"
        assert record.workspace == workspace / "custom_project" / "manual_runs" / "custom_run"
        assert reopened.runs_dir == workspace / "custom_project" / "model_runs"
        assert reopened.model_specs_dir == workspace / "custom_project" / "spec_library"

        manual_manifest = record.workspace / "run.toml"
        loaded = load_run_record(manual_manifest)
        assert loaded.workspace == record.workspace
        assert loaded.package_versions["chd"] == "v1"
        assert reopened.discover_runs() == []
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_register_run_record_roundtrip_with_completed_status_and_metadata():
    workspace = _project_temp_dir("project_catalog_register_run_record")
    try:
        catalog = ProjectCatalog(workspace / "record_project", name="record_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model"))

        record = RunRecord(
            run_id="completed_run",
            model_spec="tiny_model",
            workspace=catalog.runs_dir / "completed_run",
            status="completed",
            completed_at="2026-04-24T12:00:00+00:00",
            tags=["baseline", "shared"],
            package_versions={"rch": "rch_v1", "drn": "drn_v2"},
            metadata={"author": "luke", "purpose": "sharing"},
            parameter_overrides={"rch_factor": 0.8},
        )
        manifest_path = catalog.register_run_record(record)
        loaded = catalog.load_run("completed_run")

        assert manifest_path == record.manifest_path
        assert loaded.status == "completed"
        assert loaded.completed_at == "2026-04-24T12:00:00+00:00"
        assert loaded.metadata["purpose"] == "sharing"
        assert loaded.parameter_overrides["rch_factor"] == 0.8
        assert loaded.package_versions["drn"] == "drn_v2"

        with pytest.raises(FileExistsError, match="Run manifest already exists"):
            catalog.register_run_record(record)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_loader_and_compare_support_realistic_run_analysis():
    workspace = _project_temp_dir("project_catalog_loader_compare")
    try:
        catalog = ProjectCatalog(workspace / "analysis_project", name="analysis_project")
        catalog.register_model_spec(
            ModelSpec(
                name="tiny_model",
                description="two-cell steady-state model",
                grid_ref="two_cell_grid",
                default_packages={"oc": "oc_v1", "chd": "chd_v1"},
            )
        )

        baseline = catalog.create_run(
            RunSpec(
                run_id="baseline_run",
                model_spec="tiny_model",
                scenario="baseline",
                package_versions={"oc": "oc_v1", "chd": "chd_v1", "rch": "rch_v1"},
            )
        )
        variant = catalog.create_run(
            RunSpec(
                run_id="variant_run",
                model_spec="tiny_model",
                scenario="variant",
                package_versions={"oc": "oc_v1", "chd": "chd_v2", "rch": "rch_v1"},
            )
        )

        _build_and_run_two_cell_model(baseline, heads=(10.0, 9.0))
        _build_and_run_two_cell_model(variant, heads=(10.0, 8.5))

        assert isinstance(catalog.loader, RunLoader)
        assert isinstance(catalog.compare, RunComparison)

        baseline_heads = catalog.loader.load_heads_array("baseline_run")
        variant_heads = catalog.loader.load_heads_array("variant_run")
        baseline_frame = catalog.loader.load_heads_frame("baseline_run")
        budget_file = catalog.loader.load_budget_file("baseline_run")

        assert np.allclose(baseline_heads, [10.0, 9.0])
        assert np.allclose(variant_heads, [10.0, 8.5])
        assert baseline_frame["cell"].tolist() == [0, 1]
        assert np.allclose(baseline_frame["head"].to_numpy(), [10.0, 9.0])
        assert budget_file is not None

        package_compare = catalog.compare.compare_package_versions("baseline_run", "variant_run")
        head_compare = catalog.compare.compare_heads("baseline_run", "variant_run")
        stat_compare = catalog.compare.compare_head_stats("baseline_run", "variant_run")
        region_compare = catalog.compare.compare_region_heads(
            "baseline_run",
            "variant_run",
            region_cells=[1],
            region_name="right_cell",
        )

        package_compare = package_compare.set_index("package")
        assert bool(package_compare.loc["oc", "changed"]) is False
        assert bool(package_compare.loc["chd", "changed"]) is True
        assert package_compare.loc["chd", "run_a"] == "chd_v1"
        assert package_compare.loc["chd", "run_b"] == "chd_v2"

        assert head_compare["cell"].tolist() == [0, 1]
        assert np.allclose(head_compare["baseline_run"].to_numpy(), [10.0, 9.0])
        assert np.allclose(head_compare["variant_run"].to_numpy(), [10.0, 8.5])
        assert np.allclose(head_compare["head_diff"].to_numpy(), [0.0, -0.5])

        stat_compare = stat_compare.set_index("metric")
        assert stat_compare.loc["count", "value"] == 2.0
        assert stat_compare.loc["mean_a", "value"] == 9.5
        assert stat_compare.loc["mean_b", "value"] == 9.25
        assert stat_compare.loc["mean_diff", "value"] == -0.25
        assert stat_compare.loc["max_abs_diff", "value"] == 0.5

        assert region_compare["region"].unique().tolist() == ["right_cell"]
        assert region_compare["cell"].tolist() == [1]
        assert np.allclose(region_compare["head_diff"].to_numpy(), [-0.5])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_package_artifact_can_be_reused_on_related_run():
    workspace = _project_temp_dir("project_catalog_package_artifact_reuse")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        source_run = catalog.create_run(
            RunSpec(
                run_id="source_run",
                model_spec="tiny_model",
                package_versions={},
            )
        )
        source_model = _build_and_run_two_cell_model(source_run, heads=(10.0, 8.75))
        catalog.attach_run(source_model, source_run)

        artifact = catalog.create_package_artifact(
            "baseline_chd",
            model=source_model,
            package_name="chd",
            description="Two-cell baseline constant heads",
            tags=["baseline", "boundary"],
        )

        assert isinstance(artifact, PackageArtifact)
        listed = catalog.list_package_artifacts().set_index("artifact_id")
        assert listed.loc["baseline_chd", "package_type"] == "chd"

        related_run = catalog.create_run(
            RunSpec(
                run_id="related_run",
                model_spec="tiny_model",
                scenario="reused_chd",
                package_versions={"chd": "baseline_chd"},
            )
        )

        vor = _two_cell_vor_clockwise()
        related_model = SimulationBase(vor=vor, nper=1, **related_run.simulation_kwargs())
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        attached_run = catalog.attach_run(related_model, related_run)
        attached_packages = catalog.apply_run_package_artifacts(related_model, related_run)
        success, _ = related_model.run_simulation()

        assert attached_run.run_id == "related_run"
        assert "chd" in attached_packages
        assert success is True
        assert related_model.run_record.run_id == "related_run"
        assert related_model.project_catalog is catalog
        assert related_model.attached_package_artifacts["chd"] == "baseline_chd"
        assert np.allclose(related_model.hds.get_data(kstpkper=(0, 0)).squeeze(), [10.0, 8.75])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simulationbase_can_register_run_and_package_artifact_during_build():
    workspace = _project_temp_dir("project_catalog_simulationbase_run_build")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        chd = CHD(
            model=source_model,
            stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.75]]},
            artifact_id="baseline_chd",
            artifact_description="Created during source model build",
        )
        success, _ = source_model.run_simulation()

        stored_run = catalog.load_run("source_run")
        stored_artifact = catalog.load_package_artifact("baseline_chd")

        assert success is True
        assert source_model.run_record.run_id == "source_run"
        assert source_model.model_spec == "tiny_model"
        assert source_model.workspace == catalog.runs_dir / "source_run"
        assert source_model.name == "source_run"
        assert source_model.run_record.package_versions["chd"] == "baseline_chd"
        assert stored_run.package_versions["chd"] == "baseline_chd"
        assert chd.package_artifact.artifact_id == "baseline_chd"
        assert stored_artifact.source_run_id == "source_run"

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            run_package_versions={"chd": "baseline_chd"},
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)
        related_model.apply_registered_package_artifacts()
        success, _ = related_model.run_simulation()

        assert success is True
        assert related_model.run_record.package_versions["chd"] == "baseline_chd"
        assert np.allclose(related_model.hds.get_data(kstpkper=(0, 0)).squeeze(), [10.0, 8.75])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simulationbase_auto_registers_run_when_catalog_metadata_are_provided():
    workspace = _project_temp_dir("project_catalog_simulationbase_auto_register")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="auto_run",
            model_spec="tiny_model",
            package_versions={"chd": "baseline_chd"},
        )

        assert model.run_record is not None
        assert model.run_record.run_id == "auto_run"
        assert model.run_record.package_versions["chd"] == "baseline_chd"
        assert model.workspace == catalog.runs_dir / "auto_run"

        loaded_run = catalog.load_run("auto_run")
        assert loaded_run.run_id == "auto_run"
        assert loaded_run.package_versions["chd"] == "baseline_chd"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_recharge_wrapper_stays_quiet_and_headsplus_reuses_one_reader(capsys):
    workspace = _project_temp_dir("project_catalog_wrapper_cleanup")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(run_id="quiet_run", model_spec="tiny_model"),
            overwrite=True,
        )

        vor = _two_cell_vor_clockwise()
        model = SimulationBase(vor=vor, nper=1, **record.simulation_kwargs())
        DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(model=model, stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 9.0]]})
        _ = capsys.readouterr()
        Recharge(
            model=model,
            rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.00075]]},
        )

        captured = capsys.readouterr()
        assert captured.out == ""

        success, _ = model.run_simulation()
        assert success is True
        heads = model.hds
        assert heads.hds is heads
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_headsplus_observation_and_plot_api_still_work_after_split():
    workspace = _project_temp_dir("headsplus_split_helpers")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(run_id="heads_run", model_spec="tiny_model"),
            overwrite=True,
        )

        model = _build_and_run_two_cell_rch_uzf_model(record)
        obs_path = _write_gpkg(
            workspace / "obs_points.gpkg",
            gpd.GeoDataFrame(
                {"ExploName": ["obs_a", "obs_b"]},
                geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
                crs=model.vor.crs,
            ),
        )

        obs_cells = model.hds.get_obs_cells(obs_path)
        obs_heads = model.hds.get_obs_heads(obs_path)
        obs_heads_long = model.hds.get_obs_heads(obs_path, long_format=True)
        fig = model.hds.plot_heads(locs=[0, 1], plot_fig=False, return_fig=True)
        choro = model.hds.choropleth(kstpkper=(0, 0), layer=0)

        assert obs_cells == {"obs_a": 0, "obs_b": 1}
        assert list(obs_heads.columns) == ["obs_a", "obs_b"]
        assert obs_heads_long.index.names == ["locs", "layer", "kstpkper"]
        assert len(fig.data) == 2
        assert len(choro.data) >= 1
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_budget_obs_respects_zero_based_budget_indexing():
    workspace = _project_temp_dir("budget_obs_zero_based")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="four_cell_grid"))
        record = catalog.create_run(
            RunSpec(run_id="budget_run", model_spec="tiny_model"),
            overwrite=True,
        )

        model = _build_and_run_four_cell_standard_budget_model(record)
        zone_path = _write_gpkg(
            workspace / "drn_zone.gpkg",
            gpd.GeoDataFrame(
                {"name": ["drn_zone"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=model.vor.crs,
            ),
        )

        frame = model.bud("drn").plot_budget_obs(shp_gpkg=zone_path, plot_fig=False)

        assert "drn_zone" in frame.columns
        assert frame["drn_zone"].iloc[0] > 0
        assert frame["total"].iloc[0] >= frame["drn_zone"].iloc[0]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_drn_budget_df_is_zero_based_exactly_once():
    workspace = _project_temp_dir("budget_df_zero_based_once")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="four_cell_grid"))
        record = catalog.create_run(
            RunSpec(run_id="drn_budget_run", model_spec="tiny_model"),
            overwrite=True,
        )

        model = _build_and_run_four_cell_standard_budget_model(record)
        drn_df = model.bud("drn").df.reset_index()

        nodes = drn_df["node"].astype(int).tolist()
        assert nodes == [0]
        assert min(nodes) == 0
        assert all(node >= 0 for node in nodes)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_drn_budget_plot_choro_uses_zero_based_cell_positions():
    workspace = _project_temp_dir("budget_plot_choro_zero_based")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="four_cell_grid"))
        record = catalog.create_run(
            RunSpec(run_id="drn_plot_run", model_spec="tiny_model"),
            overwrite=True,
        )

        model = _build_and_run_four_cell_standard_budget_model(record)
        captured: dict[str, object] = {}

        class _DummyPlot:
            def plot(self):
                captured["plotted"] = True

        def fake_choro(*, per, custom_zs, zmin, zmax):
            captured["per"] = per
            captured["custom_zs"] = list(custom_zs)
            captured["zmin"] = zmin
            captured["zmax"] = zmax
            return _DummyPlot()

        model.choro = fake_choro
        DRNBudget(model).plot_choro(per=0)

        custom_zs = captured["custom_zs"]
        assert captured["per"] == 0
        assert captured["plotted"] is True
        assert custom_zs[0] > 0
        assert custom_zs[1:] == [0, 0, 0]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_outputs_namespace_is_the_canonical_home_for_lak_and_sfr_stage():
    workspace = _project_temp_dir("outputs_stage_namespace")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="four_cell_grid"))

        source_run = catalog.create_run(
            RunSpec(run_id="source_run", model_spec="tiny_model"),
            overwrite=True,
        )
        vor = _four_cell_vor_clockwise()
        source_model = SimulationBase(vor=vor, nper=1, **source_run.simulation_kwargs())
        DisvGrid(vor=vor, model=source_model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[11.5] * 4)
        KFlow(model=source_model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(model=source_model, stress_period_data={0: [[(0, 1), 11.0], [(0, 3), 10.5]]})

        lake_path = _write_gpkg(
            workspace / "lake_stage.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=vor.crs,
            ),
        )
        stream_path = _write_gpkg(
            workspace / "stream_stage.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream"]},
                geometry=[LineString([(0.2, 1.8), (1.8, 1.2)])],
                crs=vor.crs,
            ),
        )

        lake_connections = LakeConnectionData(
            model=source_model,
            vor=vor,
            paths=[lake_path],
            bed_leakance=0.1,
            horizontal_connections={0: [11.0, 9.0]},
            use_reconciled_surfaces=False,
        )
        connectiondata = lake_connections.connection_data
        lake_packagedata = LakePackageData(nlakes=1, starting_stage=[11.0], connectiondata=connectiondata)
        lake_perioddata = LakePeriodData(model=source_model, lake_ids=[0], lake_stages=[11.0], status=["ACTIVE"])
        LAKPackage(
            model=source_model,
            nlakes=1,
            noutlets=0,
            ntables=0,
            packagedata=lake_packagedata.packagedata,
            connectiondata=connectiondata,
            perioddata=lake_perioddata.perioddata,
            mover=False,
        )
        SFR(
            model=source_model,
            vor=vor,
            stream_paths=[stream_path],
            inflows={0: [(0, 0.5)]},
            widths=5.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            mover=False,
            add_sfr=True,
        )
        success, _ = source_model.run_simulation()
        assert success is True

        lak_from_outputs = source_model.outputs.lak.stage.get()
        lak_from_wrapper = LakStage(source_model).get()
        sfr_from_outputs = source_model.outputs.sfr.stage.get()
        sfr_from_wrapper = SFRStage(source_model).get()

        assert np.allclose(lak_from_outputs, lak_from_wrapper)
        pd.testing.assert_frame_equal(sfr_from_outputs, sfr_from_wrapper)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_budget_requires_package_before_df():
    workspace = _project_temp_dir("budget_requires_package")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(run_id="budget_base_run", model_spec="tiny_model"),
            overwrite=True,
        )
        model = _build_and_run_two_cell_model(record)

        with pytest.raises(ValueError, match="gwf_package is required"):
            _ = model.bud().df
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_lak_and_sfr_budget_get_require_bud_type():
    workspace = _project_temp_dir("package_budget_requires_type")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="four_cell_grid"))

        source_run = catalog.create_run(
            RunSpec(run_id="source_run", model_spec="tiny_model"),
            overwrite=True,
        )
        vor = _four_cell_vor_clockwise()
        source_model = SimulationBase(vor=vor, nper=1, **source_run.simulation_kwargs())
        DisvGrid(vor=vor, model=source_model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[11.5] * 4)
        KFlow(model=source_model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(model=source_model, stress_period_data={0: [[(0, 1), 11.0], [(0, 3), 10.5]]})

        lake_path = _write_gpkg(
            workspace / "lake_budget_type.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[Polygon([(0.0, 0.0), (0.95, 0.0), (0.95, 0.95), (0.0, 0.95)])],
                crs=vor.crs,
            ),
        )
        stream_path = _write_gpkg(
            workspace / "stream_budget_type.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream"]},
                geometry=[LineString([(0.2, 1.8), (1.8, 1.2)])],
                crs=vor.crs,
            ),
        )

        lake_connections = LakeConnectionData(
            model=source_model,
            vor=vor,
            paths=[lake_path],
            bed_leakance=0.1,
            horizontal_connections={0: [11.0, 9.0]},
            use_reconciled_surfaces=False,
        )
        connectiondata = lake_connections.connection_data
        lake_packagedata = LakePackageData(nlakes=1, starting_stage=[11.0], connectiondata=connectiondata)
        lake_perioddata = LakePeriodData(model=source_model, lake_ids=[0], lake_stages=[11.0], status=["ACTIVE"])
        LAKPackage(
            model=source_model,
            nlakes=1,
            noutlets=0,
            ntables=0,
            packagedata=lake_packagedata.packagedata,
            connectiondata=connectiondata,
            perioddata=lake_perioddata.perioddata,
            mover=False,
        )
        SFR(
            model=source_model,
            vor=vor,
            stream_paths=[stream_path],
            inflows={0: [(0, 0.5)]},
            widths=5.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            mover=False,
            add_sfr=True,
        )
        success, _ = source_model.run_simulation()
        assert success is True

        with pytest.raises(ValueError, match="bud_type is required when requesting LAK budget data"):
            source_model.outputs.lak.bud.get()
        with pytest.raises(ValueError, match="bud_type is required when requesting SFR budget data"):
            source_model.outputs.sfr.bud.get()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_lak_and_sfr_budget_warn_and_return_raw_when_period_lengths_do_not_match(monkeypatch):
    class _DummyPackageBudget:
        def get_data(self, text=None):
            return [{"dummy": 1}]

    class _DummyPackage:
        class _Output:
            def budget(self):
                return _DummyPackageBudget()

        class _Count:
            data = 1

        output = _Output()
        nreaches = _Count()

    class _DummyModel:
        kstpkper = [(0, 0), (0, 1)]
        lak = _DummyPackage()
        sfr = _DummyPackage()

    def always_fail(*args, **kwargs):
        raise ValueError("mismatch")

    monkeypatch.setattr(budget_module, "_package_output_df", always_fail)

    lak_accessor = budget_module.LakBudget(_DummyModel())
    with pytest.warns(UserWarning, match="Length of LAK budget periods and record arrays do not match"):
        lak_raw = lak_accessor.get("FLOW")
    assert lak_raw == [{"dummy": 1}]

    sfr_accessor = budget_module.SFRBudget(_DummyModel())
    with pytest.warns(UserWarning, match="Length of SFR budget periods and record arrays do not match"):
        sfr_raw = sfr_accessor.get("FLOW")
    assert sfr_raw == [{"dummy": 1}]


def test_simulationbase_can_reuse_ic_npf_rch_artifacts_on_related_run():
    workspace = _project_temp_dir("project_catalog_reuse_ic_npf_rch")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(
            model=source_model,
            vor=vor,
            nlay=1,
            strt=[10.0, 8.75],
            artifact_id="baseline_ic",
        )
        KFlow(
            model=source_model,
            k=[1.0, 1.0],
            save_specific_discharge=False,
            artifact_id="baseline_npf",
        )
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(
            model=source_model,
            stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.75]]},
            artifact_id="baseline_chd",
        )
        Recharge(
            model=source_model,
            rch_dict={0: [[(0, 0), 0.001], [(0, 1), 0.001]]},
            artifact_id="baseline_rch",
        )
        success, _ = source_model.run_simulation()
        assert success is True

        stored_run = catalog.load_run("source_run")
        assert stored_run.package_versions["ic"] == "baseline_ic"
        assert stored_run.package_versions["npf"] == "baseline_npf"
        assert stored_run.package_versions["chd"] == "baseline_chd"
        assert stored_run.package_versions["rch"] == "baseline_rch"

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            package_versions={
                "ic": "baseline_ic",
                "npf": "baseline_npf",
                "chd": "baseline_chd",
                "rch": "baseline_rch",
            },
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        attached = related_model.apply_registered_package_artifacts()
        success, _ = related_model.run_simulation()

        assert success is True
        assert set(attached) == {"ic", "npf", "chd", "rch"}
        assert np.allclose(related_model.hds.get_data(kstpkper=(0, 0)).squeeze(), [10.0, 8.75])
        listed = catalog.list_package_artifacts().set_index("artifact_id")
        assert listed.loc["baseline_ic", "package_type"] == "ic"
        assert listed.loc["baseline_npf", "package_type"] == "npf"
        assert listed.loc["baseline_rch", "package_type"] == "rch"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_related_run_artifact_workflow_can_reuse_and_modify_packages_on_small_optimized_voronoi_model():
    workspace = _project_temp_dir("project_catalog_related_run_feature_workflow")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="small_feature_model", grid_ref="optimized_small_voronoi"))

        vor, mesh_report, geometries = _build_small_optimized_voronoi_model(workspace)
        assert mesh_report["profile"] == "balanced"
        assert mesh_report["mode"] == "optimized"
        assert mesh_report["optimization"]["iterations_requested"] == 3
        assert mesh_report["quality"]["voronoi_status"] == "built"
        assert mesh_report["quality"]["duplicate_vertex_count"] == 0
        assert mesh_report["quality"]["zero_area_triangle_count"] == 0

        top = vor.gdf_topbtm[0].tolist()
        bottom = vor.gdf_topbtm[1].tolist()
        centroids_x = np.asarray(vor.centroids_x, dtype=float)
        left_cells = [cell for cell, x in enumerate(centroids_x) if x < 260.0]
        right_cells = [cell for cell, x in enumerate(centroids_x) if x >= 260.0]
        east_cells = sorted(set(vor.get_vor_cells_as_series(geometries["east_boundary"]).iloc[0]))
        assert left_cells and right_cells and east_cells

        lake_path = _write_gpkg(
            workspace / "lake.gpkg",
            gpd.GeoDataFrame(
                {"name": ["lake_0"]},
                geometry=[geometries["lake"]],
                crs=vor.crs,
            ),
        )
        stream_path = _write_gpkg(
            workspace / "stream.gpkg",
            gpd.GeoDataFrame(
                {"name": ["stream_0"]},
                geometry=[geometries["stream"]],
                crs=vor.crs,
            ),
        )

        baseline_east_head = 98.75
        baseline_rch_rows = []
        baseline_finf = []
        for cell in range(vor.ncpl):
            if cell in left_cells:
                baseline_rch_rows.append([(0, cell), 0.0018])
                baseline_finf.append(0.00045)
            else:
                baseline_rch_rows.append([(0, cell), 0.0011])
                baseline_finf.append(0.00025)

        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="small_feature_model",
        )
        DisvGrid(vor=vor, model=source_model, top=top, bottom=[bottom], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(
            model=source_model,
            vor=vor,
            nlay=1,
            strt=(np.asarray(top) - 4.0).tolist(),
            artifact_id="baseline_ic",
        )
        KFlow(
            model=source_model,
            k=[12.0] * vor.ncpl,
            save_specific_discharge=False,
            artifact_id="baseline_npf",
        )
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(
            model=source_model,
            stress_period_data={0: [[(0, cell), baseline_east_head] for cell in east_cells]},
            artifact_id="baseline_chd",
        )
        Recharge(
            model=source_model,
            rch_dict={0: baseline_rch_rows},
            artifact_id="baseline_rch",
        )
        UZFPackageData(
            model=source_model,
            vor=vor,
            uzf_cells=[(0, cell) for cell in range(vor.ncpl)],
            vks=0.45,
            thtr=0.1,
            thts=0.3,
            thti=0.2,
            finf={0: baseline_finf},
            add_uzf=True,
            artifact_id="baseline_uzf",
        )

        lake_connections = LakeConnectionData(
            model=source_model,
            vor=vor,
            paths=[lake_path],
            bed_leakance=0.05,
            horizontal_connections={0: [100.0, 95.0]},
            use_reconciled_surfaces=False,
            only_vertical=True,
        )
        lake_packagedata = LakePackageData(
            nlakes=1,
            starting_stage=[100.0],
            connectiondata=lake_connections.connection_data,
        )
        lake_perioddata = LakePeriodData(
            model=source_model,
            lake_ids=[0],
            lake_stages=[100.0],
            status=["ACTIVE"],
        )
        LAKPackage(
            model=source_model,
            nlakes=1,
            noutlets=0,
            ntables=0,
            packagedata=lake_packagedata.packagedata,
            connectiondata=lake_connections.connection_data,
            perioddata=lake_perioddata.perioddata,
            mover=False,
            artifact_id="baseline_lak",
        )
        SFR(
            model=source_model,
            vor=vor,
            stream_paths=[stream_path],
            inflows={0: [(0, 0.2)]},
            widths=10.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            mover=False,
            add_sfr=True,
            artifact_id="baseline_sfr",
        )

        success, _ = source_model.run_simulation()
        assert success is True
        source_heads = source_model.hds.get_data(kstpkper=(0, 0)).squeeze()

        catalog.derive_package_artifact(
            "baseline_chd",
            artifact_id="variant_chd",
            description="Lower east boundary stages",
            package_data_updates={
                "stress_period_data": {
                    "0": [[[0, int(cell)], baseline_east_head - 1.0] for cell in east_cells],
                }
            },
        )
        catalog.derive_package_artifact(
            "baseline_rch",
            artifact_id="variant_rch",
            description="Reduced recharge for related scenario",
            package_data_updates={
                "stress_period_data": {
                    "0": [
                        [[0, int(cell)], 0.00135 if cell in left_cells else 0.0007]
                        for cell in range(vor.ncpl)
                    ],
                }
            },
        )

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="small_feature_model",
            package_versions={
                "ic": "baseline_ic",
                "npf": "baseline_npf",
                "chd": "variant_chd",
                "rch": "variant_rch",
                "uzf": "baseline_uzf",
                "lak": "baseline_lak",
                "sfr": "baseline_sfr",
            },
        )
        DisvGrid(vor=vor, model=related_model, top=top, bottom=[bottom], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        validation = related_model.validate_referenced_package_artifacts()
        assert not validation.empty
        assert set(validation["package_key"]) == {"ic", "npf", "chd", "rch", "uzf", "lak", "sfr"}
        assert bool(validation["valid"].all()) is True

        success, _ = related_model.run_simulation()
        assert success is True
        related_heads = related_model.hds.get_data(kstpkper=(0, 0)).squeeze()

        assert related_model.pending_package_artifact_keys() == []
        assert set(related_model.attached_package_artifacts) == {"ic", "npf", "chd", "rch", "uzf", "lak", "sfr"}
        assert related_model.run_record.package_versions["chd"] == "variant_chd"
        assert related_model.run_record.package_versions["rch"] == "variant_rch"
        assert not np.allclose(source_heads, related_heads)
        assert float(np.mean(related_heads)) < float(np.mean(source_heads))

        related_outdir = related_model.workspace
        assert (related_outdir / "related_run.chd").exists()
        assert (related_outdir / "related_run.rch").exists()
        assert (related_outdir / "related_run.uzf").exists()
        assert (related_outdir / "related_run.lak").exists()
        assert (related_outdir / "related_run.sfr").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simulationbase_auto_applies_registered_package_artifacts_on_run():
    workspace = _project_temp_dir("project_catalog_auto_apply_package_artifacts")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(
            model=source_model,
            stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.75]]},
            artifact_id="baseline_chd",
        )
        success, _ = source_model.run_simulation()
        assert success is True

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            package_versions={"chd": "baseline_chd"},
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        assert related_model.auto_apply_package_artifacts is True
        assert related_model.pending_package_artifact_keys() == ["chd"]

        success, _ = related_model.run_simulation()

        assert success is True
        assert related_model.pending_package_artifact_keys() == []
        assert related_model.attached_package_artifacts["chd"] == "baseline_chd"
        assert np.allclose(related_model.hds.get_data(kstpkper=(0, 0)).squeeze(), [10.0, 8.75])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simulationbase_summary_is_backward_compatible_with_missing_grid_type_override():
    workspace = _project_temp_dir("simulationbase_summary_backward_compatibility")
    try:
        record = RunRecord(
            run_id="legacy_like_run",
            model_spec="tiny_model",
            workspace=workspace / "legacy_like_run",
            status="created",
        )
        model = _build_and_run_two_cell_model(record)

        if hasattr(model, "_grid_type_override"):
            delattr(model, "_grid_type_override")

        summary = model.summary()

        assert model.grid_type == "disv"
        assert summary.set_index("name").loc["legacy_like_run", "grid_type"] == "disv"
        assert summary.set_index("name").loc["legacy_like_run", "source"] == "simulation"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_validate_referenced_package_artifacts_reports_ready_related_run():
    workspace = _project_temp_dir("project_catalog_validate_ready_artifacts")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(
            model=source_model,
            stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.75]]},
            artifact_id="baseline_chd",
        )

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            package_versions={"chd": "baseline_chd"},
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        report = related_model.validate_referenced_package_artifacts()

        assert list(report["package_key"]) == ["chd"]
        assert bool(report.loc[0, "valid"]) is True
        assert report.loc[0, "message"] == "ok"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_uzf_package_artifact_can_be_created_and_reused_on_same_grid():
    workspace = _project_temp_dir("project_catalog_reuse_uzf")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(model=source_model, stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 9.0]]})
        uzf_source = UZFPackageData(
            model=source_model,
            vor=vor,
            uzf_cells=[(0, 0), (0, 1)],
            vks=1.0,
            thtr=0.1,
            thts=0.3,
            thti=0.2,
            finf={0: [0.125, 0.05]},
            add_uzf=True,
            artifact_id="baseline_uzf",
        )
        success, _ = source_model.run_simulation()
        assert success is True
        assert uzf_source.package_artifact.artifact_id == "baseline_uzf"

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            run_package_versions={"uzf": "baseline_uzf"},
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)
        CHD(model=related_model, stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 9.0]]})
        attached = related_model.apply_registered_package_artifacts()
        success, _ = related_model.run_simulation()

        assert success is True
        assert "uzf" in attached
        assert related_model.run_record.package_versions["uzf"] == "baseline_uzf"
        assert (related_model.workspace / "related_run.uzf").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_lak_and_sfr_package_artifacts_can_be_reused_on_same_grid():
    workspace = _project_temp_dir("project_catalog_reuse_lak_sfr")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="four_cell_grid"))

        vor = _four_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=source_model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[11.5] * 4)
        KFlow(model=source_model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(model=source_model, stress_period_data={0: [[(0, 1), 11.0], [(0, 3), 10.5]]})

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
                {"name": ["stream"]},
                geometry=[LineString([(0.2, 1.8), (1.8, 1.2)])],
                crs=vor.crs,
            ),
        )

        lake_connections = LakeConnectionData(
            model=source_model,
            vor=vor,
            paths=[lake_path],
            bed_leakance=0.1,
            horizontal_connections={0: [11.0, 9.0]},
            use_reconciled_surfaces=False,
        )
        connectiondata = lake_connections.connection_data
        lake_packagedata = LakePackageData(nlakes=1, starting_stage=[11.0], connectiondata=connectiondata)
        lake_perioddata = LakePeriodData(model=source_model, lake_ids=[0], lake_stages=[11.0], status=["ACTIVE"])
        lake_package = LAKPackage(
            model=source_model,
            nlakes=1,
            noutlets=0,
            ntables=0,
            packagedata=lake_packagedata.packagedata,
            connectiondata=connectiondata,
            perioddata=lake_perioddata.perioddata,
            mover=False,
            artifact_id="baseline_lak",
        )
        sfr_source = SFR(
            model=source_model,
            vor=vor,
            stream_paths=[stream_path],
            inflows={0: [(0, 0.5)]},
            widths=5.0,
            gradients=0.001,
            mannings=0.03,
            streambed_k=1.0,
            streambed_thickness=1.0,
            mover=False,
            add_sfr=True,
            artifact_id="baseline_sfr",
        )
        success, _ = source_model.run_simulation()
        assert success is True
        assert lake_package.package_artifact.artifact_id == "baseline_lak"
        assert sfr_source.package_artifact.artifact_id == "baseline_sfr"

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            run_package_versions={"lak": "baseline_lak", "sfr": "baseline_sfr"},
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=related_model, top=[12.0] * 4, bottom=[[0.0] * 4], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[11.5] * 4)
        KFlow(model=related_model, k=[1.0] * 4, save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)
        CHD(model=related_model, stress_period_data={0: [[(0, 1), 11.0], [(0, 3), 10.5]]})
        attached = related_model.apply_registered_package_artifacts()
        success, _ = related_model.run_simulation()

        assert success is True
        assert set(attached) == {"lak", "sfr"}
        assert (related_model.workspace / "related_run.lak").exists()
        assert (related_model.workspace / "related_run.sfr").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_mvr_package_artifact_can_be_reused_when_dependencies_exist():
    workspace = _project_temp_dir("project_catalog_reuse_mvr")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        vor = _two_cell_vor_clockwise()
        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        UZFPackage(
            model=source_model,
            nuzfcells=2,
            packagedata=[
                [0, (0, 0), 1, -1, 0.001, 1.0, 0.1, 0.3, 0.2, 4.0],
                [1, (0, 1), 1, -1, 0.001, 1.0, 0.1, 0.3, 0.2, 4.0],
            ],
            perioddata={0: [[0, 0.1, 0, 0, 0, 0, 0, 0], [1, 0.1, 0, 0, 0, 0, 0, 0]]},
            mover=True,
            artifact_id="baseline_uzf_pkg",
        )
        LAKPackage(
            model=source_model,
            nlakes=1,
            noutlets=0,
            ntables=0,
            packagedata=[[0, 9.5, 1]],
            connectiondata=[[0, 0, (0, 0), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
            perioddata={0: [[0, "STATUS", "ACTIVE"]]},
            mover=True,
            artifact_id="baseline_lak_pkg",
        )
        mvr = MVRPackage(
            model=source_model,
            maxmvr=1,
            maxpackages=2,
            packages=[["lak"], ["uzf"]],
            perioddata={0: [["lak", 0, "uzf", 1, "FACTOR", 1.0]]},
            artifact_id="baseline_mvr",
        )
        assert mvr.package_artifact.artifact_id == "baseline_mvr"

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            run_package_versions={
                "uzf": "baseline_uzf_pkg",
                "lak": "baseline_lak_pkg",
                "mvr": "baseline_mvr",
            },
            register_run=True,
            overwrite_run=True,
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 9.0])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)
        attached = related_model.apply_registered_package_artifacts()
        related_model.sim.write_simulation()

        assert set(attached) == {"uzf", "lak", "mvr"}
        assert (related_model.workspace / "related_run.mvr").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_simulationbase_package_versions_alias_updates_registered_run():
    workspace = _project_temp_dir("project_catalog_package_versions_alias")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        vor = _two_cell_vor_clockwise()

        model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="alias_run",
            model_spec="tiny_model",
            package_versions={"chd": "baseline_chd"},
            parameter_overrides={"rch_factor": 1.2},
            register_run=True,
            overwrite_run=True,
        )

        loaded = catalog.load_run("alias_run")
        assert model.run_record.package_versions["chd"] == "baseline_chd"
        assert loaded.package_versions["chd"] == "baseline_chd"
        assert loaded.parameter_overrides["rch_factor"] == 1.2
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_validate_referenced_package_artifacts_reports_missing_dependencies():
    workspace = _project_temp_dir("project_catalog_validate_missing_dependencies")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        vor = _two_cell_vor_clockwise()

        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        UZFPackage(
            model=source_model,
            nuzfcells=1,
            packagedata=[[0, (0, 0), 1, -1, 0.001, 1.0, 0.1, 0.3, 0.2, 4.0]],
            perioddata={0: [[0, 0.1, 0, 0, 0, 0, 0, 0]]},
            mover=True,
            artifact_id="baseline_uzf_pkg",
        )
        LAKPackage(
            model=source_model,
            nlakes=1,
            noutlets=0,
            ntables=0,
            packagedata=[[0, 10.0, 1]],
            connectiondata=[[0, 0, (0, 1), "vertical", 1.0, 1.0, 1.0, 0.0, 0.0]],
            perioddata={0: [["rainfall", 0, 0.0]]},
            mover=True,
            artifact_id="baseline_lak_pkg",
        )
        MVRPackage(
            model=source_model,
            maxmvr=1,
            maxpackages=2,
            packages=[["lak"], ["uzf"]],
            perioddata={0: [["lak", 0, "uzf", 0, "factor", 1.0]]},
            artifact_id="baseline_mvr",
        )

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            package_versions={"mvr": "baseline_mvr"},
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        report = related_model.validate_referenced_package_artifacts()
        assert bool(report.loc[0, "valid"]) is False
        assert report.loc[0, "missing_dependencies"] == ["lak", "uzf"]
        assert "missing dependent packages: lak, uzf" in report.loc[0, "message"]

        with pytest.raises(PackageCompatibilityError, match="missing dependent packages: lak, uzf"):
            related_model.validate_referenced_package_artifacts(raise_on_error=True)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_can_derive_package_artifact_with_lineage_and_modified_data():
    workspace = _project_temp_dir("project_catalog_derive_artifact")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        vor = _two_cell_vor_clockwise()

        source_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="source_run",
            model_spec="tiny_model",
        )
        DisvGrid(vor=vor, model=source_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=source_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=source_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=source_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=source_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=source_model)
        CHD(
            model=source_model,
            stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.75]]},
            artifact_id="baseline_chd",
        )

        derived = catalog.derive_package_artifact(
            "baseline_chd",
            artifact_id="raised_chd",
            description="Raised constant heads",
            package_data_updates={
                "stress_period_data": {
                    "0": [[[0, 0], 11.0], [[0, 1], 9.25]],
                }
            },
        )

        assert derived.derived_from_artifact_id == "baseline_chd"
        assert derived.description == "Raised constant heads"

        related_model = SimulationBase(
            vor=vor,
            nper=1,
            project_catalog=catalog,
            run_id="related_run",
            model_spec="tiny_model",
            package_versions={"chd": "raised_chd"},
        )
        DisvGrid(vor=vor, model=related_model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=related_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=related_model, vor=vor, nlay=1, strt=[10.0, 8.75])
        KFlow(model=related_model, k=[1.0, 1.0], save_specific_discharge=False)
        Storage(model=related_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=related_model)

        success, _ = related_model.run_simulation()
        assert success is True
        assert np.allclose(related_model.hds.get_data(kstpkper=(0, 0)).squeeze(), [11.0, 9.25])

        listed = catalog.list_package_artifacts().set_index("artifact_id")
        assert listed.loc["raised_chd", "derived_from_artifact_id"] == "baseline_chd"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_package_artifact_rejects_incompatible_grid():
    workspace = _project_temp_dir("project_catalog_package_artifact_compatibility")
    try:
        catalog = ProjectCatalog(workspace / "artifact_project", name="artifact_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="mixed_grids"))

        source_run = catalog.create_run(RunSpec(run_id="source_run", model_spec="tiny_model"))
        source_model = _build_and_run_two_cell_model(source_run, heads=(10.0, 9.0))
        catalog.attach_run(source_model, source_run)
        catalog.create_package_artifact("baseline_chd", model=source_model, package_name="chd")

        incompatible_run = catalog.create_run(
            RunSpec(
                run_id="incompatible_run",
                model_spec="tiny_model",
                package_versions={"chd": "baseline_chd"},
            )
        )

        vor = _three_cell_vor_clockwise()
        incompatible_model = SimulationBase(vor=vor, nper=1, **incompatible_run.simulation_kwargs())
        DisvGrid(vor=vor, model=incompatible_model, top=[10.0, 10.0, 10.0], bottom=[[0.0, 0.0, 0.0]], nlay=1)
        TemporalDiscretization(model=incompatible_model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=incompatible_model, vor=vor, nlay=1, strt=[10.0, 9.5, 9.0])
        KFlow(model=incompatible_model, k=[1.0, 1.0, 1.0], save_specific_discharge=False)
        Storage(model=incompatible_model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=incompatible_model)

        catalog.attach_run(incompatible_model, incompatible_run)

        with pytest.raises(PackageCompatibilityError, match="different grid layout|node_count"):
            catalog.apply_run_package_artifacts(incompatible_model, incompatible_run)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_run_loader_can_open_file_backed_disv_run_without_model_object():
    workspace = _project_temp_dir("project_catalog_file_backed_disv")
    try:
        catalog = ProjectCatalog(workspace / "file_backed_project", name="file_backed_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(
                run_id="file_run",
                model_spec="tiny_model",
                package_versions={"chd": "chd_v1"},
            )
        )

        _build_and_run_two_cell_model(record, heads=(10.0, 8.75))
        model_object_path = record.get_path("model_object_file")
        if model_object_path.exists():
            model_object_path.unlink()

        loaded = catalog.load_run_model("file_run")
        assert loaded._vor is None
        summary = loaded.summary()
        assert loaded._vor is None
        package_summary = loaded.package_summary()
        output_summary = loaded.output_summary()
        grid_summary = loaded.grid_summary()
        result_summary = loaded.result_summary()
        file_summary = loaded.file_summary()
        assert loaded._vor is None
        heads = loaded.hds.get_data(kstpkper=(0, 0)).squeeze()
        all_heads = loaded.all_heads
        fig = loaded.hds.plot_heads(locs=[0, 1], plot_fig=False, return_fig=True)
        sim_axes = loaded.sim.plot()
        run_axes = loaded.plot()
        gwf_axes = loaded.gwf.plot()

        loaded.add_region_from_cells("left_cell", cellids=[(0, 0)], overwrite=True)
        loaded.add_region_from_geometry(
            "right_geom",
            loaded.vor.gdf_vorPolys.geometry.iloc[1],
            overwrite=True,
        )
        loaded.add_group("all_cells", members=["left_cell", "right_geom"], overwrite=True)

        assert isinstance(loaded, LoadedMf6Run)
        assert loaded.source == "mf6_files"
        assert loaded.grid_type == "disv"
        assert loaded.package_names == ["CHD", "DISV", "IC", "NPF", "OC", "STO"]
        assert loaded.record.run_id == "file_run"
        assert loaded.workspace == record.workspace
        assert summary.set_index("name").loc["file_run", "grid_type"] == "disv"
        assert package_summary.set_index("package").loc["CHD", "package_type"] == "chd"
        assert bool(output_summary.iloc[0]["has_heads"]) is True
        assert grid_summary.iloc[0]["ncpl"] == 2
        assert result_summary.iloc[0]["head_mean"] == 9.375
        assert {"input", "output"} == set(file_summary["category"])
        assert loaded.vor.ncpl == 2
        assert loaded._vor is not None
        assert np.allclose(heads, [10.0, 8.75])
        assert np.allclose(all_heads["elev"].astype(float).to_numpy(), [10.0, 8.75])
        assert fig is not None
        assert sim_axes is not None
        assert run_axes is not None
        assert gwf_axes is not None
        plt.close("all")
        assert {path.name for path in loaded.list_input_files()} >= {"mfsim.nam", "file_run.disv", "file_run.nam"}
        assert {path.name for path in loaded.list_output_files()} >= {"mfsim.lst", "file_run.hds", "file_run.cbc"}
        assert loaded.resolve_region_cells("all_cells") == [0, 1]
        assert loaded.region_heads("all_cells", per=0)["cell"].tolist() == [0, 1]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_loaded_run_progress_messages_explain_open_and_lazy_voronoi(capsys):
    workspace = _project_temp_dir("project_catalog_progress_messages")
    try:
        catalog = ProjectCatalog(workspace / "progress_project", name="progress_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(
                run_id="progress_run",
                model_spec="tiny_model",
                package_versions={"chd": "chd_v1"},
            )
        )

        _build_and_run_two_cell_model(record, heads=(10.0, 8.75))
        model_object_path = record.get_path("model_object_file")
        if model_object_path.exists():
            model_object_path.unlink()

        loaded = catalog.load_run_model("progress_run", verbosity_level=1)
        open_output = capsys.readouterr().out

        assert "[RunLoader] Opening run 'progress_run'" in open_output
        assert "[LoadedMf6Run] Prepared lazy file-backed model 'progress_run'" in open_output
        assert "Grid, outputs, and package definitions will load only when needed." in open_output
        assert loaded._gwf is None
        assert loaded._vor is None

        _ = loaded.vor
        vor_output = capsys.readouterr().out
        assert "[LoadedMf6Run] Loading MF6 packages for core access:" in vor_output
        assert "[LoadedMf6Run] Locating GWF model 'progress_run'" in vor_output
        assert "[LoadedMf6Run] Building Voronoi grid view from the MF6 modelgrid" in vor_output
        assert "[LoadedMf6Run] Voronoi grid ready (ncpl=2)" in vor_output
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_load_mf6_run_supports_direct_workspace_loading_for_disu_runs():
    workspace = _project_temp_dir("project_catalog_file_backed_disu")
    try:
        record = RunRecord(
            run_id="disu_file_run",
            model_spec="tiny_model",
            workspace=workspace / "disu_file_run",
            status="completed",
        )
        _build_and_run_two_cell_disu_model(record, heads=(10.0, 8.5))
        model_object_path = record.get_path("model_object_file")
        if model_object_path.exists():
            model_object_path.unlink()

        loaded = load_mf6_run(record.workspace)
        heads = loaded.hds.get_data(kstpkper=(0, 0)).squeeze()

        assert isinstance(loaded, LoadedMf6Run)
        assert loaded.grid_type == "disu"
        assert loaded.name == "disu_file_run"
        assert loaded.vor.ncpl == 2
        assert np.allclose(heads, [10.0, 8.5])
        assert loaded.gwf.modelgrid.__class__.__name__ == "UnstructuredGrid"
        assert {path.name for path in loaded.list_input_files()} >= {"disu_file_run.disu", "disu_file_run.nam", "mfsim.nam"}
        assert {path.name for path in loaded.list_output_files()} >= {"disu_file_run.hds", "disu_file_run.cbc", "mfsim.lst"}
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_loaded_run_lazily_loads_grid_outputs_and_selected_packages():
    workspace = _project_temp_dir("project_catalog_lazy_load_selected_packages")
    try:
        catalog = ProjectCatalog(workspace / "lazy_project", name="lazy_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(
            RunSpec(
                run_id="lazy_run",
                model_spec="tiny_model",
                package_versions={"rch": "rch_v1", "uzf": "uzf_v1"},
            )
        )

        _build_and_run_two_cell_rch_uzf_model(record)
        model_object_path = record.get_path("model_object_file")
        if model_object_path.exists():
            model_object_path.unlink()

        loaded = catalog.load_run_model("lazy_run")

        assert loaded._gwf is None
        assert loaded._vor is None
        assert loaded.package_names == ["CHD", "DISV", "IC", "NPF", "OC", "RCH", "STO", "UZF"]

        summary = loaded.summary().iloc[0]
        package_summary = loaded.package_summary()
        output_summary = loaded.output_summary().iloc[0]

        assert loaded._gwf is None
        assert loaded._vor is None
        assert summary["ncpl"] is None or np.isnan(summary["ncpl"])
        assert summary["grid_type"] == "disv"
        assert summary["source"] == "mf6_files"
        assert set(package_summary["package"]) >= {"DISV", "RCH", "UZF"}
        assert bool(output_summary["has_heads"]) is True
        assert bool(output_summary["has_budget"]) is True

        vor = loaded.vor
        assert loaded._gwf is not None
        assert vor.ncpl == 2
        assert loaded._loaded_package_types >= {"DISV", "OC"}

        packages = loaded.package(["rch", "uzf"])
        assert set(packages) == {"rch", "uzf"}
        assert packages["rch"] is loaded.rch
        assert packages["uzf"] is loaded.uzf
        assert loaded._loaded_package_types >= {"DISV", "OC", "RCH", "UZF"}

        ifno_to_cellid = loaded.outputs.uzf.ifno_to_cellid
        assert ifno_to_cellid.index.name == "ifno"
        assert ifno_to_cellid["cellid"].tolist() == [0, 1]

        loaded.load_all()
        assert loaded._fully_loaded is True
        assert loaded._loaded_package_types >= {"CHD", "DISV", "IC", "NPF", "OC", "RCH", "STO", "UZF"}
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_loaded_run_budget_access_does_not_force_core_grid_loading():
    workspace = _project_temp_dir("project_catalog_lazy_budget_only")
    try:
        catalog = ProjectCatalog(workspace / "lazy_project", name="lazy_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        record = catalog.create_run(RunSpec(run_id="lazy_budget", model_spec="tiny_model"))

        _build_and_run_two_cell_model(record, heads=(10.0, 8.75))
        model_object_path = record.get_path("model_object_file")
        if model_object_path.exists():
            model_object_path.unlink()

        loaded = catalog.load_run_model("lazy_budget")
        assert loaded._gwf is None
        assert loaded._vor is None

        budget_df = loaded.bud("chd").df

        assert budget_df.empty is False
        assert loaded._gwf is None
        assert loaded._vor is None
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_heads_and_budget_compare_runs_against_reference():
    workspace = _project_temp_dir("project_catalog_model_group_heads_budget")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        pre = catalog.create_run(RunSpec(run_id="pre", model_spec="tiny_model"))
        post = catalog.create_run(RunSpec(run_id="post", model_spec="tiny_model"))

        _build_and_run_two_cell_rch_uzf_model(
            pre,
            heads=(10.0, 9.0),
            recharge=(0.0010, 0.0007),
            finf=(0.0005, 0.0002),
        )
        _build_and_run_two_cell_rch_uzf_model(
            post,
            heads=(10.0, 8.5),
            recharge=(0.0014, 0.0009),
            finf=(0.0006, 0.0003),
        )

        group = ModelGroup(
            {
                "pre": catalog.load_run_model("pre"),
                "post": catalog.load_run_model("post"),
            },
            reference="pre",
        )

        heads = group.hds.get(per=0)
        head_compare = group.hds.compare(per=0)
        chd_budget = group.bud("chd").get()
        chd_compare = group.bud("chd").compare()

        assert set(heads["model"]) == {"pre", "post"}
        assert set(heads["cell"]) == {0, 1}
        assert set(head_compare["model"]) == {"post"}
        assert set(head_compare["reference_model"]) == {"pre"}
        assert set(head_compare["cell"]) == {0, 1}
        assert np.allclose(
            head_compare.sort_values("cell")["diff"].to_numpy(),
            [0.0, -0.5],
        )

        assert set(chd_budget["model"]) == {"pre", "post"}
        assert "node" in chd_budget.columns
        assert chd_budget["node"].min() == 0
        assert set(chd_compare["model"]) == {"post"}
        assert set(chd_compare["reference_model"]) == {"pre"}
        assert "node" in chd_compare.columns
        assert chd_compare["node"].min() == 0
        assert np.any(np.abs(chd_compare["diff"].astype(float).to_numpy()) > 0.0)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_outputs_lak_stage_returns_all_models():
    workspace = _project_temp_dir("project_catalog_model_group_lak_outputs")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="lake_model", grid_ref="four_cell_grid"))

        lake_a = catalog.create_run(RunSpec(run_id="lake_a", model_spec="lake_model"))
        lake_b = catalog.create_run(RunSpec(run_id="lake_b", model_spec="lake_model"))

        _build_and_run_four_cell_lake_model(lake_a, workspace, starting_stage=11.0)
        _build_and_run_four_cell_lake_model(lake_b, workspace, starting_stage=10.7)

        group = ModelGroup(
            {
                "lake_a": catalog.load_run_model("lake_a"),
                "lake_b": catalog.load_run_model("lake_b"),
            },
            reference="lake_a",
        )

        stages = group.outputs.lak.stage()

        assert set(stages["model"]) == {"lake_a", "lake_b"}
        assert set(stages["lake"]) == {0}
        assert stages["kstpkper"].tolist() == [(0, 0), (0, 0)]
        pivot = stages.pivot(index="lake", columns="model", values="stage")
        expected_a = float(group.models["lake_a"].outputs.lak.stage.get()[0, 0])
        expected_b = float(group.models["lake_b"].outputs.lak.stage.get()[0, 0])
        assert float(pivot.loc[0, "lake_a"]) == pytest.approx(expected_a)
        assert float(pivot.loc[0, "lake_b"]) == pytest.approx(expected_b)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_tiny_lake_fixture_has_gaining_and_losing_exchange():
    workspace = _project_temp_dir("project_catalog_tiny_lake_exchange_signs")
    try:
        record = RunRecord(run_id="tiny_lake", model_spec="lake_model", workspace=workspace / "tiny_lake", status="completed")
        model = _build_and_run_four_cell_lake_model(record, workspace)

        budget = model.outputs.lak.bud.get("GWF")
        node_exchange = budget.groupby("node2", as_index=False)["q"].sum()
        cell_exchange = (
            model.packages.lak.results.q.get(per=0)
            .groupby("cell", as_index=False)["q"]
            .sum()
        )

        assert node_exchange.empty is False
        assert budget["FLOW-AREA"].gt(0.0).any()
        assert (node_exchange["q"] > 0.0).any()
        assert (node_exchange["q"] < 0.0).any()
        assert set(cell_exchange["cell"]) == {0, 2}
        assert (cell_exchange["q"] > 0.0).any()
        assert (cell_exchange["q"] < 0.0).any()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_run_explorer_group_opens_runs_and_supports_group_heads():
    workspace = _project_temp_dir("project_catalog_explorer_group")
    try:
        legacy_root = workspace / "legacy_runs"
        _build_legacy_two_cell_run(legacy_root / "pre", "pre", heads=(10.0, 9.0))
        _build_legacy_two_cell_run(legacy_root / "post", "post", heads=(10.0, 8.5))

        explorer = explore_runs(legacy_root, model_spec="legacy_archive")
        group = explorer.group(["pre", "post"], reference="pre")
        comparison = group.hds.compare(per=0)

        assert group.reference == "pre"
        assert group.run_ids == ["pre", "post"]
        assert set(comparison["model"]) == {"post"}
        assert np.allclose(
            comparison.sort_values("cell")["diff"].to_numpy(),
            [0.0, -0.5],
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_can_be_created_from_workspace_directories():
    workspace = _project_temp_dir("project_catalog_model_group_from_dirs")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "pre", "pre", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(legacy_root / "post", "post", heads=(10.0, 8.5))

        group = ModelGroup([first.workspace, second.workspace], reference="pre")
        comparison = group.hds.compare(per=0)

        assert group.run_ids == ["pre", "post"]
        assert set(comparison["model"]) == {"post"}
        assert np.allclose(
            comparison.sort_values("cell")["diff"].to_numpy(),
            [0.0, -0.5],
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_can_be_created_from_named_workspace_mapping():
    workspace = _project_temp_dir("project_catalog_model_group_named_dirs")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "pre", "pre", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(legacy_root / "post", "post", heads=(10.0, 8.5))

        group = ModelGroup(
            {
                "predev": first.workspace,
                "postdev": second.workspace,
            },
            reference="predev",
        )
        heads = group.hds.get(per=0)
        comparison = group.hds.compare(per=0)

        assert group.run_ids == ["predev", "postdev"]
        assert set(heads["model"]) == {"predev", "postdev"}
        assert set(comparison["model"]) == {"postdev"}
        assert set(comparison["reference_model"]) == {"predev"}
        assert np.allclose(
            comparison.sort_values("cell")["diff"].to_numpy(),
            [0.0, -0.5],
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_shared_grid_reuses_one_voronoi_for_identical_workspaces():
    workspace = _project_temp_dir("project_catalog_model_group_shared_grid")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "pre", "pre", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(legacy_root / "post", "post", heads=(10.0, 8.5))

        group = ModelGroup([first.workspace, second.workspace], reference="pre", shared_grid=True)
        pre = group.models["pre"]
        post = group.models["post"]

        assert getattr(pre, "_shared_vor_source", None) is None
        assert getattr(post, "_shared_vor_source", None) is pre
        assert pre._vor is None
        assert post._vor is None

        pre_vor = pre.vor
        post_vor = post.vor

        assert pre_vor is post_vor
        assert pre._vor is pre_vor
        assert post._vor is None
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_shared_grid_rejects_mismatched_discretizations():
    workspace = _project_temp_dir("project_catalog_model_group_shared_grid_mismatch")
    try:
        two_cell = RunRecord(run_id="two", model_spec="tiny", workspace=workspace / "two", status="completed")
        four_cell = RunRecord(run_id="four", model_spec="standard", workspace=workspace / "four", status="completed")

        _build_and_run_two_cell_model(two_cell, heads=(10.0, 9.0))
        _build_and_run_four_cell_standard_budget_model(four_cell)

        with pytest.raises(ValueError, match="same discretization"):
            ModelGroup([two_cell.workspace, four_cell.workspace], reference="two", shared_grid=True)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_budget_supports_standard_boundary_packages():
    workspace = _project_temp_dir("project_catalog_model_group_standard_packages")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="standard_model", grid_ref="four_cell_grid"))

        base = catalog.create_run(RunSpec(run_id="base", model_spec="standard_model"))
        variant = catalog.create_run(RunSpec(run_id="variant", model_spec="standard_model"))

        _build_and_run_four_cell_standard_budget_model(
            base,
            chd_heads=(11.0, 10.5),
            drn_elev=10.8,
            drn_cond=0.4,
            ghb_head=9.8,
            ghb_cond=0.35,
            recharge=(0.0008, 0.0008, 0.0006, 0.0006),
        )
        _build_and_run_four_cell_standard_budget_model(
            variant,
            chd_heads=(10.8, 10.2),
            drn_elev=10.6,
            drn_cond=0.55,
            ghb_head=9.6,
            ghb_cond=0.5,
            recharge=(0.0011, 0.0010, 0.0008, 0.00075),
        )

        group = ModelGroup(
            {
                "base": catalog.load_run_model("base"),
                "variant": catalog.load_run_model("variant"),
            },
            reference="base",
        )

        for package in ["drn", "ghb", "chd", "rch"]:
            budget = group.bud(package).get()
            comparison = group.bud(package).compare()

            assert budget.empty is False
            assert comparison.empty is False
            assert set(budget["model"]) == {"base", "variant"}
            assert set(comparison["model"]) == {"variant"}
            assert set(comparison["reference_model"]) == {"base"}
            if "node" in budget.columns:
                assert budget["node"].min() >= 0
            if "node2" in budget.columns:
                assert budget["node2"].min() >= 0
            if "node" in comparison.columns:
                assert comparison["node"].min() >= 0
            if "node2" in comparison.columns:
                assert comparison["node2"].min() >= 0
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_budget_access_does_not_force_grid_build_for_file_backed_runs():
    workspace = _project_temp_dir("project_catalog_group_budget_only")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "pre", "pre", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(legacy_root / "post", "post", heads=(10.0, 8.5))

        group = ModelGroup([first.workspace, second.workspace], reference="pre")
        pre = group.models["pre"]
        post = group.models["post"]

        assert pre._gwf is None
        assert post._gwf is None
        assert pre._vor is None
        assert post._vor is None

        budget = group.bud("chd").get()

        assert budget.empty is False
        assert pre._gwf is None
        assert post._gwf is None
        assert pre._vor is None
        assert post._vor is None
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_package_inputs_support_standard_boundary_packages():
    workspace = _project_temp_dir("project_catalog_model_group_package_inputs")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="standard_model", grid_ref="four_cell_grid"))

        base = catalog.create_run(RunSpec(run_id="base", model_spec="standard_model"))
        variant = catalog.create_run(RunSpec(run_id="variant", model_spec="standard_model"))

        _build_and_run_four_cell_standard_budget_model(
            base,
            chd_heads=(11.0, 10.5),
            drn_elev=10.8,
            drn_cond=0.4,
            ghb_head=9.8,
            ghb_cond=0.35,
            recharge=(0.0008, 0.0008, 0.0006, 0.0006),
        )
        _build_and_run_four_cell_standard_budget_model(
            variant,
            chd_heads=(10.8, 10.2),
            drn_elev=10.6,
            drn_cond=0.55,
            ghb_head=9.6,
            ghb_cond=0.5,
            recharge=(0.0011, 0.0010, 0.0008, 0.00075),
        )

        group = ModelGroup(
            {
                "base": catalog.load_run_model("base"),
                "variant": catalog.load_run_model("variant"),
            },
            reference="base",
        )

        package_cases = {
            "rch": ("recharge", "recharge_diff"),
            "chd": ("head", "head_diff"),
            "drn": ("cond", "cond_diff"),
            "ghb": ("cond", "cond_diff"),
        }

        for accessor_name, (value_column, diff_column) in package_cases.items():
            accessor = getattr(group, accessor_name)
            data = accessor.get()
            comparison = accessor.compare()

            assert data.empty is False
            assert comparison.empty is False
            assert set(data["model"]) == {"base", "variant"}
            assert set(comparison["model"]) == {"variant"}
            assert set(comparison["reference_model"]) == {"base"}
            assert set(data["per"]) == {0}
            assert data["layer"].min() >= 0
            assert data["cell"].min() >= 0
            assert value_column in data.columns
            assert diff_column in comparison.columns
            assert np.any(np.abs(comparison[diff_column].astype(float).to_numpy()) > 0.0)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_packages_cell_inputs_support_get_summary_and_map():
    workspace = _project_temp_dir("project_catalog_model_packages_inputs")
    try:
        record = RunRecord(run_id="cell_inputs", model_spec="standard", workspace=workspace / "cell_inputs", status="completed")
        model = _build_and_run_four_cell_standard_budget_model(
            record,
            chd_heads=(11.0, 10.5),
            drn_elev=10.8,
            drn_cond=0.4,
            ghb_head=9.8,
            ghb_cond=0.35,
            recharge=(0.0008, 0.0008, 0.0006, 0.0006),
        )

        package_cases = {
            "rch": ("recharge", [0.0008, 0.0008, 0.0006, 0.0006]),
            "chd": ("head", [0.0, 11.0, 0.0, 10.5]),
            "drn": ("elev", [10.8, 0.0, 0.0, 0.0]),
            "ghb": ("bhead", [0.0, 0.0, 9.8, 0.0]),
        }

        for package_name, (value_column, expected_zs) in package_cases.items():
            accessor = getattr(model.packages, package_name).inputs
            data = accessor.get(per=0)
            summary = accessor.summary()
            choro = accessor.map(per=0, value_column=value_column)

            assert data.empty is False
            assert {"model", "package", "per", "layer", "cell", value_column}.issubset(data.columns)
            assert summary.iloc[0]["label"] == f"{package_name}.inputs"
            assert summary.iloc[0]["records"] == len(data)
            assert choro.type == "custom"
            assert choro.hover_heads is False
            assert choro.hover_ks is False
            np.testing.assert_allclose(np.asarray(choro.zs, dtype=float), np.asarray(expected_zs, dtype=float))
            assert choro.hover_dict["Cell"] == [0, 1, 2, 3]
            assert choro.hover_dict["Period"] == [0, 0, 0, 0]
            assert value_column in choro.hover_dict
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_packages_uzf_finf_support_get_summary_and_map():
    workspace = _project_temp_dir("project_catalog_model_packages_uzf")
    try:
        record = RunRecord(run_id="uzf_inputs", model_spec="tiny", workspace=workspace / "uzf_inputs", status="completed")
        model = _build_and_run_two_cell_rch_uzf_model(
            record,
            finf=(0.0005, 0.00025),
            recharge=(0.0010, 0.0007),
        )

        accessor = model.packages.uzf.inputs.finf
        data = accessor.get(per=0)
        summary = accessor.summary()
        choro = accessor.map(per=0)

        assert data.empty is False
        assert {"model", "package", "per", "ifno", "layer", "cell", "finf"}.issubset(data.columns)
        assert summary.iloc[0]["label"] == "uzf.inputs.finf"
        assert summary.iloc[0]["records"] == 2
        np.testing.assert_allclose(np.asarray(choro.zs, dtype=float), np.asarray([0.0005, 0.00025], dtype=float))
        assert choro.hover_dict["Cell"] == [0, 1]
        assert choro.hover_dict["finf"] == [0.0005, 0.00025]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_loaded_mf6_run_packages_namespace_supports_input_exploration():
    workspace = _project_temp_dir("project_catalog_loaded_run_packages")
    try:
        record = RunRecord(run_id="loaded_inputs", model_spec="tiny", workspace=workspace / "loaded_inputs", status="completed")
        _build_and_run_two_cell_rch_uzf_model(
            record,
            finf=(0.0005, 0.00025),
            recharge=(0.0010, 0.0007),
        )

        loaded = load_mf6_run(record.workspace)
        assert loaded._gwf is None
        assert loaded._vor is None

        rch = loaded.packages.rch.inputs.get(per=0)
        uzf_choro = loaded.packages.uzf.inputs.finf.map(per=0)

        assert rch.empty is False
        assert {"recharge", "layer", "cell"}.issubset(rch.columns)
        assert loaded._gwf is not None
        assert loaded._vor is not None
        np.testing.assert_allclose(np.asarray(uzf_choro.zs, dtype=float), np.asarray([0.0005, 0.00025], dtype=float))
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_packages_results_support_budget_and_stage_exploration():
    workspace = _project_temp_dir("project_catalog_model_package_results")
    try:
        cell_record = RunRecord(run_id="cell_results", model_spec="standard", workspace=workspace / "cell_results", status="completed")
        cell_model = _build_and_run_four_cell_standard_budget_model(cell_record)

        chd_results = cell_model.packages.chd.results.q.get(per=0)
        chd_summary = cell_model.packages.chd.results.q.summary()
        chd_map = cell_model.packages.chd.results.q.map(per=0)
        uzf_results_model = _build_and_run_two_cell_rch_uzf_model(
            RunRecord(run_id="uzf_results", model_spec="tiny", workspace=workspace / "uzf_results", status="completed"),
            finf=(0.0005, 0.00025),
            recharge=(0.0010, 0.0007),
        )
        uzf_gwrch = uzf_results_model.packages.uzf.results.gwrch.get(per=0)
        uzf_sat = uzf_results_model.packages.uzf.results.sat.get(per=0)
        uzf_sat_map = uzf_results_model.packages.uzf.results.sat.map(per=0)

        assert chd_results.empty is False
        assert {"model", "package", "kstpkper", "per", "layer", "cell", "q"}.issubset(chd_results.columns)
        assert chd_summary.iloc[0]["label"] == "chd.results.q"
        assert len(chd_map.zs) == cell_model.vor.ncpl
        assert {"gwrch", "cell", "layer"}.issubset(uzf_gwrch.columns)
        assert {"sat", "cell", "layer"}.issubset(uzf_sat.columns)
        assert len(uzf_sat_map.zs) == uzf_results_model.vor.ncpl

        routed_record = RunRecord(run_id="stage_results", model_spec="sfr_lak", workspace=workspace / "stage_results", status="completed")
        routed_model = _build_and_run_four_cell_lak_sfr_model(routed_record, workspace)

        lak_connections = routed_model.packages.lak.connections.get()
        lak_connections_summary = routed_model.packages.lak.connections.summary()
        lak_connections_map = routed_model.packages.lak.connections.map()
        lak_stage = routed_model.packages.lak.results.stage.get(per=0)
        lak_stage_map = routed_model.packages.lak.results.stage.map(per=0)
        lak_stage_fig = routed_model.packages.lak.results.stage.plot_timeseries()
        lak_q = routed_model.packages.lak.results.q.get(per=0)
        lak_q_summary = routed_model.packages.lak.results.q.budget_summary(per=0)
        lak_q_map = routed_model.packages.lak.results.q.map(per=0)
        lak_q_fig = routed_model.packages.lak.results.q.plot_budget(per=0)
        sfr_stage = routed_model.packages.sfr.results.stage.get(per=0)
        sfr_stage_map = routed_model.packages.sfr.results.stage.map(per=0)
        sfr_q = routed_model.packages.sfr.results.q.get(per=0)
        sfr_q_map = routed_model.packages.sfr.results.q.map(per=0)
        sfr_q_profile = routed_model.packages.sfr.results.q.profile(per=0)
        sfr_q_profile_fig = routed_model.packages.sfr.results.q.plot_profile(per=0)
        sfr_q_profile_by_reach_fig = routed_model.packages.sfr.results.q.plot_profile(per=0, x="reach")
        sfr_stage_profile = routed_model.packages.sfr.results.stage.profile(per=0)
        sfr_stage_profile_fig = routed_model.packages.sfr.results.stage.plot_profile(per=0)
        sfr_long_profile = routed_model.packages.sfr.results.long_profile(per=0)
        sfr_long_profile_fig = routed_model.packages.sfr.results.plot_long_profile(per=0)
        sfr_long_profile_reach_fig = routed_model.packages.sfr.results.plot_long_profile(
            per=0,
            x="reach",
            include_exchange=False,
        )

        assert lak_connections.empty is False
        assert {"model", "package", "lake", "iconn", "layer", "cell", "claktype", "connection_area"}.issubset(
            lak_connections.columns
        )
        assert float(lak_connections_summary.loc[0, "total_connection_area"]) > 0.0
        assert lak_stage.empty is False
        assert {"model", "package", "per", "lake", "layer", "cell", "stage"}.issubset(lak_stage.columns)
        assert lak_q.empty is False
        assert {"q", "q_per_area", "flow_area", "claktype", "lake", "cell"}.issubset(lak_q.columns)
        assert lak_q_summary.empty is False
        assert {"per", "lake", "claktype", "record_count", "q", "flow_area", "q_per_area"}.issubset(
            lak_q_summary.columns
        )
        lak_q_by_cell = lak_q.groupby("cell", as_index=False)["q"].sum()
        assert (lak_q_by_cell["q"] > 0.0).any()
        assert (lak_q_by_cell["q"] < 0.0).any()
        assert sfr_stage.empty is False
        assert {"model", "package", "per", "reach", "layer", "cell", "stage"}.issubset(sfr_stage.columns)
        assert {"q", "q_per_length", "reach", "cell", "rlen"}.issubset(sfr_q.columns)
        assert len(lak_connections_map.zs) == routed_model.vor.ncpl
        assert len(lak_stage_map.zs) == routed_model.vor.ncpl
        assert len(lak_q_map.zs) == routed_model.vor.ncpl
        assert len(sfr_stage_map.zs) == routed_model.vor.ncpl
        assert len(sfr_q_map.zs) == routed_model.vor.ncpl
        assert {"distance_start", "distance_mid", "distance_end", "rlen"}.issubset(sfr_q_profile.columns)
        assert {"distance_start", "distance_mid", "distance_end", "rlen"}.issubset(sfr_stage_profile.columns)
        assert {
            "streambed_top",
            "streambed_bottom",
            "stage",
            "q",
            "distance_mid",
            "rgrd",
            "rwid",
        }.issubset(sfr_long_profile.columns)
        assert sfr_q_profile["reach"].is_monotonic_increasing
        assert sfr_q_profile["distance_mid"].is_monotonic_increasing
        assert sfr_stage_profile["reach"].is_monotonic_increasing
        assert sfr_stage_profile["distance_mid"].is_monotonic_increasing
        assert sfr_long_profile["reach"].is_monotonic_increasing
        assert sfr_long_profile["distance_mid"].is_monotonic_increasing
        expected_sfr_q_map_frame = sfr_q.groupby("cell", as_index=False).agg({"q": "sum", "rlen": "sum"})
        expected_sfr_q_map = dict(
            zip(
                expected_sfr_q_map_frame["cell"].astype(int),
                expected_sfr_q_map_frame["q"].astype(float) / expected_sfr_q_map_frame["rlen"].astype(float),
                strict=False,
            )
        )
        for cell, normalized_q in expected_sfr_q_map.items():
            assert sfr_q_map.zs[int(cell)] == pytest.approx(normalized_q)
        expected_lak_q_map_frame = lak_q.groupby("cell", as_index=False).agg({"q": "sum", "flow_area": "sum"})
        expected_lak_q_map = dict(
            zip(
                expected_lak_q_map_frame["cell"].astype(int),
                expected_lak_q_map_frame["q"].astype(float) / expected_lak_q_map_frame["flow_area"].astype(float),
                strict=False,
            )
        )
        expected_lak_connections_map = dict(
            zip(
                lak_connections.groupby("cell", as_index=False)["connection_area"].sum()["cell"].astype(int),
                lak_connections.groupby("cell", as_index=False)["connection_area"].sum()["connection_area"].astype(float),
                strict=False,
            )
        )
        for cell, connection_area in expected_lak_connections_map.items():
            assert lak_connections_map.zs[int(cell)] == pytest.approx(connection_area)
        for cell, normalized_q in expected_lak_q_map.items():
            assert lak_q_map.zs[int(cell)] == pytest.approx(normalized_q)
        assert len(lak_stage_fig.axes) == 1
        assert len(lak_stage_fig.axes[0].lines) >= 1
        assert len(lak_q_fig.axes) == 1
        assert len(lak_q_fig.axes[0].patches) >= 1
        assert len(sfr_q_profile_fig.data) == 1
        assert len(sfr_q_profile_by_reach_fig.data) == 1
        assert len(sfr_stage_profile_fig.data) == 1
        assert [trace.name for trace in sfr_long_profile_fig.data] == [
            "Streambed Top",
            "Streambed Bottom",
            "Stage",
            "Exchange q",
        ]
        assert [trace.name for trace in sfr_long_profile_reach_fig.data] == [
            "Streambed Top",
            "Streambed Bottom",
            "Stage",
        ]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_uzf_finf_supports_get_and_compare():
    workspace = _project_temp_dir("project_catalog_model_group_uzf_finf")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        pre = catalog.create_run(RunSpec(run_id="pre", model_spec="tiny_model"))
        post = catalog.create_run(RunSpec(run_id="post", model_spec="tiny_model"))

        _build_and_run_two_cell_rch_uzf_model(
            pre,
            finf=(0.0005, 0.0002),
            recharge=(0.0010, 0.0007),
        )
        _build_and_run_two_cell_rch_uzf_model(
            post,
            finf=(0.00065, 0.00035),
            recharge=(0.0012, 0.00085),
        )

        group = ModelGroup(
            {
                "pre": catalog.load_run_model("pre"),
                "post": catalog.load_run_model("post"),
            },
            reference="pre",
        )

        data = group.uzf.finf.get()
        comparison = group.uzf.finf.compare()

        assert data.empty is False
        assert comparison.empty is False
        assert set(data["model"]) == {"pre", "post"}
        assert set(comparison["model"]) == {"post"}
        assert set(comparison["reference_model"]) == {"pre"}
        assert set(data["per"]) == {0}
        assert set(data["ifno"]) == {0, 1}
        assert set(data["cell"]) == {0, 1}
        assert set(data["layer"]) == {0}
        assert np.allclose(
            comparison.sort_values("ifno")["finf_diff"].astype(float).to_numpy(),
            [0.00015, 0.00015],
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_package_results_support_get_compare_and_maps():
    workspace = _project_temp_dir("project_catalog_group_package_results")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="standard_model", grid_ref="four_cell_grid"))

        base = catalog.create_run(RunSpec(run_id="base", model_spec="standard_model"))
        variant = catalog.create_run(RunSpec(run_id="variant", model_spec="standard_model"))

        _build_and_run_four_cell_standard_budget_model(
            base,
            chd_heads=(11.0, 10.5),
            drn_elev=10.8,
            drn_cond=0.4,
            ghb_head=9.8,
            ghb_cond=0.35,
            recharge=(0.0008, 0.0008, 0.0006, 0.0006),
        )
        _build_and_run_four_cell_standard_budget_model(
            variant,
            chd_heads=(10.8, 10.2),
            drn_elev=10.6,
            drn_cond=0.55,
            ghb_head=9.6,
            ghb_cond=0.5,
            recharge=(0.0011, 0.0010, 0.0008, 0.00075),
        )

        group = ModelGroup(
            {
                "base": catalog.load_run_model("base"),
                "variant": catalog.load_run_model("variant"),
            },
            reference="base",
        )

        data = group.packages.chd.results.q.get(per=0)
        comparison = group.packages.chd.results.q.compare(model_name="variant", per=0)
        raw_map = group.packages.chd.results.q.map(model_name="variant", per=0)
        diff_map = group.packages.chd.results.q.compare_map(model_name="variant", per=0)

        assert data.empty is False
        assert comparison.empty is False
        assert {"q", "reference_q", "q_diff"}.issubset(comparison.columns)
        assert len(raw_map.zs) == group.models["variant"].vor.ncpl
        assert len(diff_map.zs) == group.models["base"].vor.ncpl

        uzf_catalog = ProjectCatalog(workspace / "uzf_group_project", name="uzf_group_project")
        uzf_catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))
        pre = uzf_catalog.create_run(RunSpec(run_id="pre", model_spec="tiny_model"))
        post = uzf_catalog.create_run(RunSpec(run_id="post", model_spec="tiny_model"))
        _build_and_run_two_cell_rch_uzf_model(pre, finf=(0.0005, 0.0002), recharge=(0.0010, 0.0007))
        _build_and_run_two_cell_rch_uzf_model(post, finf=(0.00065, 0.00035), recharge=(0.0012, 0.00085))
        uzf_group = ModelGroup(
            {
                "pre": uzf_catalog.load_run_model("pre"),
                "post": uzf_catalog.load_run_model("post"),
            },
            reference="pre",
        )
        sat = uzf_group.packages.uzf.results.sat.get()
        gwrch_diff = uzf_group.packages.uzf.results.gwrch.compare(model_name="post")
        assert sat.empty is False
        assert gwrch_diff.empty is False

        sfr_catalog = ProjectCatalog(workspace / "sfr_group_project", name="sfr_group_project")
        sfr_catalog.register_model_spec(ModelSpec(name="sfr_model", grid_ref="four_cell_grid"))
        sfr_a = sfr_catalog.create_run(RunSpec(run_id="sfr_a", model_spec="sfr_model"))
        sfr_b = sfr_catalog.create_run(RunSpec(run_id="sfr_b", model_spec="sfr_model"))
        _build_and_run_four_cell_lak_sfr_model(sfr_a, workspace, sfr_inflow=0.5)
        _build_and_run_four_cell_lak_sfr_model(sfr_b, workspace, sfr_inflow=0.75)
        sfr_group = ModelGroup(
            {
                "sfr_a": sfr_catalog.load_run_model("sfr_a"),
                "sfr_b": sfr_catalog.load_run_model("sfr_b"),
            },
            reference="sfr_a",
        )
        sfr_q = sfr_group.packages.sfr.results.q.get()
        sfr_q_diff = sfr_group.packages.sfr.results.q.compare(model_name="sfr_b")
        sfr_q_map = sfr_group.packages.sfr.results.q.compare_map(model_name="sfr_b", per=0)
        assert sfr_q.empty is False
        assert sfr_q_diff.empty is False
        assert {"q_per_length", "q_per_length_diff"}.issubset(sfr_q_diff.columns)
        assert len(sfr_q_map.zs) == sfr_group.models["sfr_a"].vor.ncpl
        base_frame = (
            sfr_q.loc[sfr_q["model"] == "sfr_a"]
            .groupby("cell", as_index=False)
            .agg({"q": "sum", "rlen": "sum"})
        )
        variant_frame = (
            sfr_q.loc[sfr_q["model"] == "sfr_b"]
            .groupby("cell", as_index=False)
            .agg({"q": "sum", "rlen": "sum"})
        )
        base_norm = dict(
            zip(
                base_frame["cell"].astype(int),
                base_frame["q"].astype(float) / base_frame["rlen"].astype(float),
                strict=False,
            )
        )
        variant_norm = dict(
            zip(
                variant_frame["cell"].astype(int),
                variant_frame["q"].astype(float) / variant_frame["rlen"].astype(float),
                strict=False,
            )
        )
        for cell, base_value in base_norm.items():
            diff = variant_norm.get(cell, 0.0) - base_value
            assert sfr_q_map.zs[int(cell)] == pytest.approx(diff)

        lak_catalog = ProjectCatalog(workspace / "lak_group_project", name="lak_group_project")
        lak_catalog.register_model_spec(ModelSpec(name="lak_model", grid_ref="four_cell_grid"))
        lak_a = lak_catalog.create_run(RunSpec(run_id="lak_a", model_spec="lak_model"))
        lak_b = lak_catalog.create_run(RunSpec(run_id="lak_b", model_spec="lak_model"))
        _build_and_run_four_cell_lake_model(lak_a, workspace, starting_stage=10.5)
        _build_and_run_four_cell_lake_model(lak_b, workspace, starting_stage=11.0)
        lak_group = ModelGroup(
            {
                "lak_a": lak_catalog.load_run_model("lak_a"),
                "lak_b": lak_catalog.load_run_model("lak_b"),
            },
            reference="lak_a",
        )
        lak_q = lak_group.packages.lak.results.q.get()
        lak_q_diff = lak_group.packages.lak.results.q.compare(model_name="lak_b")
        lak_q_map = lak_group.packages.lak.results.q.compare_map(model_name="lak_b", per=0)
        assert lak_q.empty is False
        assert lak_q_diff.empty is False
        assert {"q_per_area", "q_per_area_diff"}.issubset(lak_q_diff.columns)
        assert len(lak_q_map.zs) == lak_group.models["lak_a"].vor.ncpl
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_build_lak_q_map_payload_normalizes_exchange_by_area():
    frame = pd.DataFrame(
        {
            "package": ["lak", "lak", "lak"],
            "per": [0, 0, 0],
            "layer": [0, 0, 0],
            "cell": [0, 0, 1],
            "lake": [0, 0, 0],
            "claktype": ["VERTICAL", "HORIZONTAL", "VERTICAL"],
            "q": [10.0, -4.0, -3.0],
            "flow_area": [2.0, 2.0, 1.5],
        }
    )

    values, hover = build_lak_q_map_payload(frame, ncpl=3, per=0, layer=0)

    assert values[0] == pytest.approx((10.0 - 4.0) / (2.0 + 2.0))
    assert values[1] == pytest.approx(-3.0 / 1.5)
    assert values[2] == pytest.approx(0.0)
    assert hover["q"][0] == pytest.approx(6.0)
    assert hover["flow_area"][0] == pytest.approx(4.0)
    assert hover["q_per_area"][0] == pytest.approx(1.5)


def test_model_group_packages_namespace_matches_existing_group_accessors():
    workspace = _project_temp_dir("project_catalog_group_packages_namespace")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        pre = catalog.create_run(RunSpec(run_id="pre", model_spec="tiny_model"))
        post = catalog.create_run(RunSpec(run_id="post", model_spec="tiny_model"))

        _build_and_run_two_cell_rch_uzf_model(
            pre,
            finf=(0.0005, 0.0002),
            recharge=(0.0010, 0.0007),
        )
        _build_and_run_two_cell_rch_uzf_model(
            post,
            finf=(0.00065, 0.00035),
            recharge=(0.0012, 0.00085),
        )

        group = ModelGroup(
            {
                "pre": catalog.load_run_model("pre"),
                "post": catalog.load_run_model("post"),
            },
            reference="pre",
        )

        pd.testing.assert_frame_equal(
            group.packages.rch.inputs.get().sort_index(axis=1),
            group.rch.get().sort_index(axis=1),
        )
        pd.testing.assert_frame_equal(
            group.packages.rch.inputs.compare().sort_index(axis=1),
            group.rch.compare().sort_index(axis=1),
        )
        pd.testing.assert_frame_equal(
            group.packages.uzf.inputs.finf.get().sort_index(axis=1),
            group.uzf.finf.get().sort_index(axis=1),
        )
        pd.testing.assert_frame_equal(
            group.packages.uzf.inputs.finf.compare().sort_index(axis=1),
            group.uzf.finf.compare().sort_index(axis=1),
        )

        uzf_results = group.packages.uzf.results.gwrch.get()
        assert {"model", "package", "kstpkper", "per", "layer", "cell", "gwrch"}.issubset(uzf_results.columns)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_model_group_package_maps_support_raw_and_difference_choropleths():
    workspace = _project_temp_dir("project_catalog_group_package_maps")
    try:
        catalog = ProjectCatalog(workspace / "group_project", name="group_project")
        catalog.register_model_spec(ModelSpec(name="tiny_model", grid_ref="two_cell_grid"))

        pre = catalog.create_run(RunSpec(run_id="pre", model_spec="tiny_model"))
        post = catalog.create_run(RunSpec(run_id="post", model_spec="tiny_model"))

        _build_and_run_two_cell_rch_uzf_model(
            pre,
            finf=(0.0005, 0.0002),
            recharge=(0.0010, 0.0007),
        )
        _build_and_run_two_cell_rch_uzf_model(
            post,
            finf=(0.00065, 0.00035),
            recharge=(0.0012, 0.00085),
        )

        group = ModelGroup(
            {
                "pre": catalog.load_run_model("pre"),
                "post": catalog.load_run_model("post"),
            },
            reference="pre",
        )

        rch_map = group.packages.rch.inputs.map(model_name="post", per=0)
        rch_diff_map = group.packages.rch.inputs.compare_map(model_name="post", per=0)
        uzf_map = group.packages.uzf.inputs.finf.map(model_name="post", per=0)
        uzf_diff_map = group.packages.uzf.inputs.finf.compare_map(model_name="post", per=0)

        np.testing.assert_allclose(np.asarray(rch_map.zs, dtype=float), np.asarray([0.0012, 0.00085], dtype=float))
        np.testing.assert_allclose(
            np.asarray(rch_diff_map.zs, dtype=float),
            np.asarray([0.0002, 0.00015], dtype=float),
            atol=1e-12,
        )
        np.testing.assert_allclose(np.asarray(uzf_map.zs, dtype=float), np.asarray([0.00065, 0.00035], dtype=float))
        np.testing.assert_allclose(
            np.asarray(uzf_diff_map.zs, dtype=float),
            np.asarray([0.00015, 0.00015], dtype=float),
            atol=1e-12,
        )

        assert rch_map.hover_dict["Cell"] == [0, 1]
        assert rch_map.hover_dict["Package"] == ["rch", "rch"]
        assert rch_diff_map.hover_dict["Model"] == ["post", "post"]
        assert rch_diff_map.hover_dict["Reference Model"] == ["pre", "pre"]
        np.testing.assert_allclose(
            np.asarray(rch_diff_map.hover_dict["recharge_diff"], dtype=float),
            np.asarray([0.0002, 0.00015], dtype=float),
            atol=1e-12,
        )
        np.testing.assert_allclose(
            np.asarray(uzf_diff_map.hover_dict["finf_diff"], dtype=float),
            np.asarray([0.00015, 0.00015], dtype=float),
            atol=1e-12,
        )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_discover_existing_runs_detects_realistic_legacy_run_folders_and_duplicate_names():
    workspace = _project_temp_dir("project_catalog_discovery")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "group_a" / "legacy_dup", "legacy_dup", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(
            legacy_root / "group_b" / "pest" / "legacy_dup",
            "legacy_dup",
            heads=(10.0, 8.5),
        )
        third = _build_legacy_two_cell_run(legacy_root / "silcrk_test", "silcrk_test", heads=(11.0, 10.5))

        discovered = discover_existing_runs(legacy_root)
        by_workspace = {run.workspace: run for run in discovered}

        assert len(discovered) == 3
        assert {run.workspace for run in discovered} == {first.workspace, second.workspace, third.workspace}
        assert isinstance(discovered[0], DiscoveredRun)
        assert by_workspace[first.workspace].run_id == "legacy_dup"
        assert by_workspace[second.workspace].run_id.startswith("legacy_dup__")
        assert by_workspace[second.workspace].run_id != by_workspace[first.workspace].run_id
        assert by_workspace[third.workspace].tags == ["silcrk"]
        assert by_workspace[first.workspace].metadata["grid_type"] == "disv"
        assert by_workspace[first.workspace].metadata["has_heads"] is True
        assert by_workspace[first.workspace].metadata["has_budget"] is True
        assert by_workspace[first.workspace].paths["model_object_file"] == "legacy_dup.model"
        assert by_workspace[first.workspace].paths["heads_file"] == "legacy_dup.hds"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_catalog_can_import_existing_runs_without_modifying_legacy_workspaces():
    workspace = _project_temp_dir("project_catalog_import_existing")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "group_a" / "legacy_dup", "legacy_dup", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(
            legacy_root / "group_b" / "pest" / "legacy_dup",
            "legacy_dup",
            heads=(10.0, 8.5),
        )

        catalog = ProjectCatalog(workspace / "import_project", name="import_project")
        imported = catalog.import_existing_runs(legacy_root, model_spec="legacy_archive")

        assert len(imported) == 2
        assert not (first.workspace / "run.toml").exists()
        assert not (second.workspace / "run.toml").exists()
        assert len(list(catalog.run_records_dir.glob("*.toml"))) == 2

        imported_runs = {run.workspace: run for run in catalog.discover_runs()}
        assert imported_runs[first.workspace].model_spec == "legacy_archive"
        assert imported_runs[second.workspace].model_spec == "legacy_archive"
        assert imported_runs[first.workspace].status == "imported"
        assert imported_runs[first.workspace].metadata["relative_workspace"] == "group_a/legacy_dup"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_imported_legacy_runs_support_lazy_loading_model_opening_and_comparison():
    workspace = _project_temp_dir("project_catalog_imported_compare")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "cumb_case_a", "cumb_case_a", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(legacy_root / "cumb_case_b", "cumb_case_b", heads=(10.0, 8.5))

        catalog = ProjectCatalog(workspace / "import_compare_project", name="import_compare_project")
        imported = catalog.import_existing_runs(legacy_root, model_spec="legacy_archive")
        imported_ids = {record.workspace: record.run_id for record in imported}

        loaded_model = catalog.loader.load_model_object(imported_ids[first.workspace])
        head_compare = catalog.compare.compare_heads(imported_ids[first.workspace], imported_ids[second.workspace])
        head_stats = catalog.compare.compare_head_stats(imported_ids[first.workspace], imported_ids[second.workspace])

        assert isinstance(loaded_model, SimulationBase)
        assert loaded_model.name == "cumb_case_a"
        assert np.allclose(catalog.loader.load_heads_array(imported_ids[first.workspace]), [10.0, 9.0])
        assert np.allclose(catalog.loader.load_heads_array(imported_ids[second.workspace]), [10.0, 8.5])
        assert np.allclose(head_compare["head_diff"].to_numpy(), [0.0, -0.5])
        assert head_stats.set_index("metric").loc["max_abs_diff", "value"] == 0.5
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_imported_legacy_runs_can_be_opened_from_mf6_files_without_pickle():
    workspace = _project_temp_dir("project_catalog_imported_file_open")
    try:
        legacy_root = workspace / "legacy_runs"
        _build_legacy_two_cell_run(legacy_root / "cumb_case_a", "cumb_case_a", heads=(10.0, 9.0))

        catalog = ProjectCatalog(workspace / "import_file_project", name="import_file_project")
        imported = catalog.import_existing_runs(legacy_root, model_spec="legacy_archive")
        model_object_path = imported[0].get_path("model_object_file")
        if model_object_path.exists():
            model_object_path.unlink()

        loaded = catalog.load_run_model(imported[0].run_id)
        assert isinstance(loaded, LoadedMf6Run)
        assert loaded.record.run_id == imported[0].run_id
        assert np.allclose(loaded.hds.get_data(kstpkper=(0, 0)).squeeze(), [10.0, 9.0])
        assert loaded.summary().set_index("name").loc["cumb_case_a", "source"] == "mf6_files"
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_project_helpers_support_archive_summary_and_import_workflow():
    workspace = _project_temp_dir("project_catalog_helpers")
    try:
        legacy_root = workspace / "legacy_runs"
        first = _build_legacy_two_cell_run(legacy_root / "cumb_case_a", "cumb_case_a", heads=(10.0, 9.0))
        second = _build_legacy_two_cell_run(legacy_root / "ssb_case_b", "ssb_case_b", heads=(10.0, 8.5))

        summary = summarize_discovered_runs(legacy_root).sort_values("run_id").reset_index(drop=True)
        assert summary["run_id"].tolist() == ["cumb_case_a", "ssb_case_b"]
        assert summary.set_index("run_id").loc["cumb_case_a", "family"] == "cumb"
        assert bool(summary.set_index("run_id").loc["ssb_case_b", "has_heads"]) is True

        catalog = import_run_archive(
            legacy_root,
            workspace / "imported_catalog",
            project_name="imported_catalog",
            model_spec="legacy_archive",
        )
        run_summary = catalog.run_summary().sort_values("run_id").reset_index(drop=True)

        assert catalog.load_model_spec("legacy_archive").grid_ref == "imported"
        assert run_summary["run_id"].tolist() == ["cumb_case_a", "ssb_case_b"]
        assert run_summary.set_index("run_id").loc["cumb_case_a", "family"] == "cumb"
        assert run_summary.set_index("run_id").loc["ssb_case_b", "grid_type"] == "disv"

        imported_ids = {record.workspace: record.run_id for record in catalog.discover_runs()}
        assert np.allclose(catalog.loader.load_heads_array(imported_ids[first.workspace]), [10.0, 9.0])
        assert np.allclose(catalog.loader.load_heads_array(imported_ids[second.workspace]), [10.0, 8.5])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_import_run_archive_is_idempotent_by_default():
    workspace = _project_temp_dir("project_catalog_import_idempotent")
    try:
        legacy_root = workspace / "legacy_runs"
        _build_legacy_two_cell_run(legacy_root / "cumb_case_a", "cumb_case_a", heads=(10.0, 9.0))
        _build_legacy_two_cell_run(legacy_root / "ssb_case_b", "ssb_case_b", heads=(10.0, 8.5))

        catalog_root = workspace / "imported_catalog"
        first = import_run_archive(
            legacy_root,
            catalog_root,
            project_name="imported_catalog",
            model_spec="legacy_archive",
        )
        second = import_run_archive(
            legacy_root,
            catalog_root,
            project_name="imported_catalog",
            model_spec="legacy_archive",
        )

        first_summary = first.run_summary().sort_values("run_id").reset_index(drop=True)
        second_summary = second.run_summary().sort_values("run_id").reset_index(drop=True)

        assert first_summary["run_id"].tolist() == ["cumb_case_a", "ssb_case_b"]
        assert second_summary["run_id"].tolist() == ["cumb_case_a", "ssb_case_b"]
        assert len(list((catalog_root / "run_records").glob("*.toml"))) == 2
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_import_existing_runs_example_script_and_notebook_are_valid():
    workspace = _project_temp_dir("project_catalog_example_script")
    try:
        legacy_root = workspace / "legacy_runs"
        _build_legacy_two_cell_run(legacy_root / "cumb_case_a", "cumb_case_a", heads=(10.0, 9.0))
        _build_legacy_two_cell_run(legacy_root / "ssb_case_b", "ssb_case_b", heads=(10.0, 8.5))

        script_path = ROOT / "examples" / "mf6" / "import_existing_runs_workflow.py"
        notebook_path = ROOT / "examples" / "mf6" / "notebooks" / "import_existing_runs_workflow.ipynb"
        catalog_root = workspace / "imported_catalog"

        notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
        assert notebook["nbformat"] == 4
        assert len(notebook["cells"]) >= 5

        completed = subprocess.run(
            [
                sys.executable,
                str(script_path),
                "--archive-root",
                str(legacy_root),
                "--catalog-root",
                str(catalog_root),
                "--project-name",
                "example_catalog",
                "--compare",
                "cumb_case_a",
                "ssb_case_b",
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert "Discovered Runs" in completed.stdout
        assert "Imported Catalog Summary" in completed.stdout
        assert "Head Differences: cumb_case_a vs ssb_case_b" in completed.stdout
        assert (catalog_root / "run_records" / "cumb_case_a.toml").exists()
        assert (catalog_root / "run_records" / "ssb_case_b.toml").exists()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)
