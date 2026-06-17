"""Compact synthetic PEST workflow for end-to-end calibration demos.

This module builds a small but expressive transient model with:

* a direct DISV/Voronoi grid under 1,000 cells
* several common boundary-condition types (CHD, GHB, recharge, drains)
* a deliberately wrong starting K field and drain conductances
* pseudo-observations generated from a separate "truth" model

The goal is to make K-field and drain-parameter changes obvious when PEST++
optimizes the model, while keeping runtimes manageable enough for demos and
tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace
from typing import Sequence

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

from myflopy.modflow.mf6.drn import DRNFromVector
from myflopy.modflow.mf6.kflow import KFromVector
from myflopy.modflow.mf6.observations import HeadTargets
from myflopy.modflow.mf6.pest.project import PestProject
from myflopy.modflow.mf6.pest.results import open_pest_run
from myflopy.modflow.mf6.pest.specs import (
    DrainConductanceParameter,
    ExpGeoStruct,
    HeadTargetObservationSpec,
    KPilotPointParameter,
    VectorParameterSource,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization
from myflopy.modflow.mf6.simulation.packages import (
    CHD,
    GHB,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)


@dataclass(slots=True)
class SyntheticPestDemoConfig:
    """Configuration for the compact synthetic PEST demonstration."""

    artifact_root: Path
    run_family: str = "synthetic_compact_pest_demo"
    model_name: str = "syn_demo"
    nx: int = 20
    ny: int = 20
    cell_size: float = 100.0
    nper: int = 6
    start_date: str = "2020-01-01"
    west_head: float = 100.0
    east_ghb_cond: float = 250.0
    start_k: float = 15.0
    pp_spacing: float = 300.0
    truth_drain_cond: tuple[float, float, float] = (80.0, 360.0, 950.0)
    start_drain_cond: tuple[float, float, float] = (220.0, 140.0, 180.0)
    recharge_by_period: tuple[float, ...] = (
        6.0e-4,
        1.25e-3,
        2.0e-4,
        1.85e-3,
        4.5e-4,
        1.4e-3,
    )
    east_stage_by_period: tuple[float, ...] = (
        92.0,
        91.4,
        92.8,
        90.7,
        93.1,
        91.9,
    )
    observation_fractions: tuple[tuple[float, float], ...] = (
        (0.10, 0.20),
        (0.25, 0.20),
        (0.40, 0.18),
        (0.55, 0.20),
        (0.70, 0.18),
        (0.85, 0.20),
        (0.12, 0.45),
        (0.28, 0.50),
        (0.42, 0.48),
        (0.58, 0.52),
        (0.72, 0.48),
        (0.88, 0.50),
        (0.15, 0.72),
        (0.34, 0.78),
        (0.48, 0.70),
        (0.62, 0.76),
        (0.78, 0.72),
        (0.90, 0.78),
    )
    transient_period_length_days: float = 30.0
    transient_num_steps: int = 3
    crs: str = "EPSG:2927"
    top_base: float = 102.0
    bottom_offset: float = 42.0

    def __post_init__(self):
        self.artifact_root = Path(self.artifact_root)
        if self.nper < 2:
            raise ValueError("Synthetic demo must include at least one steady and one transient period.")
        if len(self.recharge_by_period) != self.nper:
            raise ValueError("recharge_by_period length must equal nper.")
        if len(self.east_stage_by_period) != self.nper:
            raise ValueError("east_stage_by_period length must equal nper.")
        if self.nx * self.ny > 1000:
            raise ValueError("Synthetic demo grid must stay at or below 1,000 cells.")


@dataclass(slots=True)
class SyntheticPestDemoRun:
    """Artifacts returned after building or running the synthetic demo."""

    config: SyntheticPestDemoConfig
    workspace_root: Path
    inputs_dir: Path
    truth_workspace: Path
    model_workspace: Path
    pest_workspace: Path
    run_info_path: Path
    control_file: str
    targets: HeadTargets
    truth_target_values: pd.DataFrame
    observation_locations: gpd.GeoDataFrame
    start_k_source: Path
    truth_drn_source: Path
    start_drn_source: Path
    review_dir: Path
    result_workspace: Path | None = None
    review_summary: dict | None = None


def build_workspace_root(artifact_root: Path, run_family: str) -> Path:
    """Create a fresh timestamped workspace root and update the latest-pointer file."""

    artifact_root = Path(artifact_root)
    artifact_root.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    workspace_root = artifact_root / f"{run_family}_{timestamp}"
    suffix = 1
    while workspace_root.exists():
        workspace_root = artifact_root / f"{run_family}_{timestamp}_{suffix:02d}"
        suffix += 1
    workspace_root.mkdir(parents=True, exist_ok=True)
    (artifact_root / f"{run_family}_latest.txt").write_text(str(workspace_root), encoding="utf-8")
    return workspace_root


def _write_gpkg(path: Path, gdf: gpd.GeoDataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _build_rectangular_voronoi(
    *,
    nx: int,
    ny: int,
    cell_size: float,
    crs: str,
    name: str,
) -> VoronoiGridPlus:
    """Build a direct rectangular DISV/Voronoi grid without Triangle."""

    verts: list[list[float]] = []
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            verts.append([ix * cell_size, iy * cell_size])
    verts_arr = np.asarray(verts, dtype=float)

    def vertex_id(ix: int, iy: int) -> int:
        return iy * (nx + 1) + ix

    iverts: list[list[int]] = []
    centers: list[list[float]] = []
    for iy in range(ny):
        for ix in range(nx):
            ll = vertex_id(ix, iy)
            lr = vertex_id(ix + 1, iy)
            ur = vertex_id(ix + 1, iy + 1)
            ul = vertex_id(ix, iy + 1)
            iverts.append([ll, ul, ur, lr])
            centers.append([(ix + 0.5) * cell_size, (iy + 0.5) * cell_size])

    vor = VoronoiGridPlus(
        verts=verts_arr,
        iverts=iverts,
        xcyc=np.asarray(centers, dtype=float),
        crs=crs,
        name=name,
    )
    return vor


def _build_surfaces(vor: VoronoiGridPlus, config: SyntheticPestDemoConfig) -> tuple[np.ndarray, np.ndarray]:
    """Create gently sloped top and bottom surfaces for the synthetic grid."""

    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    y = centers.y.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    ymax = max(float(y.max()), 1.0)
    x_norm = x / xmax
    y_norm = y / ymax

    top = (
        config.top_base
        - 6.0 * x_norm
        - 1.2 * np.exp(-((y_norm - 0.55) / 0.14) ** 2)
        - 0.8 * np.exp(-((x_norm - 0.75) / 0.12) ** 2 - ((y_norm - 0.25) / 0.16) ** 2)
    )
    bottom = top - config.bottom_offset
    return np.asarray(top, dtype=float), np.asarray(bottom, dtype=float)


def _truth_k_field(vor: VoronoiGridPlus, config: SyntheticPestDemoConfig) -> np.ndarray:
    """Return the synthetic "truth" K field used to generate pseudo-observations."""

    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    y = centers.y.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    ymax = max(float(y.max()), 1.0)
    x_norm = x / xmax
    y_norm = y / ymax

    k = np.full(vor.ncpl, float(config.start_k), dtype=float)

    # A broad high-K corridor that should become obvious after calibration.
    corridor = np.exp(-((y_norm - 0.52) / 0.08) ** 2)
    k *= 1.0 + 4.5 * corridor

    # A low-K barrier that interrupts the corridor and creates a clear localized mismatch.
    barrier = (
        (x_norm > 0.42)
        & (x_norm < 0.56)
        & (y_norm > 0.22)
        & (y_norm < 0.84)
    )
    k[barrier] *= 0.12

    # A southeast transmissive lens to create a second nonuniform feature.
    lens = (((x_norm - 0.78) / 0.10) ** 2 + ((y_norm - 0.22) / 0.08) ** 2) < 1.0
    k[lens] *= 3.5

    return np.clip(k, 0.5, None)


def _domain_polygon(vor: VoronoiGridPlus) -> Polygon:
    return Polygon(vor.gdf_vorPolys.total_bounds.reshape(2, 2)[[0, 0, 1, 1], :])


def _build_start_k_source(vor: VoronoiGridPlus, config: SyntheticPestDemoConfig) -> gpd.GeoDataFrame:
    """Create a single-zone K source for pilot-point multiplier calibration."""

    minx, miny, maxx, maxy = vor.gdf_vorPolys.total_bounds
    polygon = Polygon([(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)])
    return gpd.GeoDataFrame(
        {
            "name": ["zone_all"],
            "k": [float(config.start_k)],
            "layer": [1],
        },
        geometry=[polygon],
        crs=config.crs,
    )


def _build_drain_sources(vor: VoronoiGridPlus, config: SyntheticPestDemoConfig) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Create truth and starting drain-feature geopackages."""

    minx, miny, maxx, maxy = vor.gdf_vorPolys.total_bounds
    width = maxx - minx
    height = maxy - miny

    def rect(x0, x1, y0, y1) -> Polygon:
        return Polygon(
            [
                (minx + width * x0, miny + height * y0),
                (minx + width * x0, miny + height * y1),
                (minx + width * x1, miny + height * y1),
                (minx + width * x1, miny + height * y0),
            ]
        )

    names = ["north_relief", "central_valley", "south_capture"]
    polygons = [
        rect(0.18, 0.84, 0.68, 0.82),
        rect(0.12, 0.90, 0.42, 0.56),
        rect(0.22, 0.78, 0.16, 0.30),
    ]
    heights = [-4.0, -6.5, -3.5]
    min_elev = [70.0, 68.0, 69.0]
    layer = [1, 1, 1]

    base = pd.DataFrame(
        {
            "name": names,
            "par_name": names,
            "height": heights,
            "layer": layer,
            "min_elev": min_elev,
        }
    )
    truth = gpd.GeoDataFrame(
        base.assign(cond=list(config.truth_drain_cond)),
        geometry=polygons,
        crs=config.crs,
    )
    start = gpd.GeoDataFrame(
        base.assign(cond=list(config.start_drain_cond)),
        geometry=polygons,
        crs=config.crs,
    )
    return truth, start


def _build_recharge_dict(vor: VoronoiGridPlus, recharge_by_period: Sequence[float]) -> dict[int, list[list[object]]]:
    """Create a uniform recharge dictionary for all periods."""

    cell_ids = [(0, int(cell)) for cell in range(vor.ncpl)]
    return {
        int(per): [[cellid, float(recharge)] for cellid in cell_ids]
        for per, recharge in enumerate(recharge_by_period)
    }


def _build_chd_spd(config: SyntheticPestDemoConfig) -> dict[int, list[list[object]]]:
    """Create west-edge CHD stress-period data."""

    spd: dict[int, list[list[object]]] = {}
    for per in range(config.nper):
        rows = []
        for iy in range(config.ny):
            cell = iy * config.nx
            rows.append([(0, cell), float(config.west_head)])
        spd[per] = rows
    return spd


def _build_ghb_spd(config: SyntheticPestDemoConfig) -> dict[int, list[list[object]]]:
    """Create east-edge GHB stress-period data."""

    spd: dict[int, list[list[object]]] = {}
    for per, stage in enumerate(config.east_stage_by_period):
        rows = []
        for iy in range(config.ny):
            cell = iy * config.nx + (config.nx - 1)
            rows.append([(0, cell), float(stage), float(config.east_ghb_cond)])
        spd[int(per)] = rows
    return spd


def _initial_heads(vor: VoronoiGridPlus, config: SyntheticPestDemoConfig) -> np.ndarray:
    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    east_start = float(config.east_stage_by_period[0])
    return np.asarray(
        config.west_head - (config.west_head - east_start) * (x / xmax),
        dtype=float,
    )


def _period_data(config: SyntheticPestDemoConfig) -> list[list[float]]:
    data = [[1.0, 1, 1.0]]
    for _ in range(config.nper - 1):
        data.append([float(config.transient_period_length_days), int(config.transient_num_steps), 1.1])
    return data


def _build_model(
    *,
    workspace: Path,
    name: str,
    config: SyntheticPestDemoConfig,
    k_array: np.ndarray,
    drn_source_path: Path,
) -> SimulationBase:
    """Build one synthetic MF6 model workspace."""

    vor = _build_rectangular_voronoi(
        nx=config.nx,
        ny=config.ny,
        cell_size=config.cell_size,
        crs=config.crs,
        name=f"{name}_grid",
    )
    top, bottom = _build_surfaces(vor, config)
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: top.astype(float),
            1: bottom.astype(float),
        },
        geometry="geometry",
        crs=config.crs,
    )

    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=config.nper)
    DisvGrid(
        vor=vor,
        model=model,
        top=top.tolist(),
        bottom=[bottom.tolist()],
        nlay=1,
        idomain=[[1] * int(vor.ncpl)],
    )
    TemporalDiscretization(model=model, period_data=_period_data(config))
    model.per_dates = pd.date_range(config.start_date, periods=config.nper, freq="MS")
    OutputControl(model=model, save_record=(("HEAD", "ALL"), ("BUDGET", "ALL")))
    InitialConditions(model=model, vor=vor, nlay=1, strt=_initial_heads(vor, config).tolist())
    KFlow(model=model, k=np.asarray(k_array, dtype=float).tolist(), k33_vert=(np.asarray(k_array, dtype=float) / 8.0).tolist())
    Storage(
        model=model,
        specific_storage=1.0e-4,
        specific_yield=0.18,
        sto_steady={0: True},
        sto_transient={per: True for per in range(1, config.nper)},
    )
    Recharge(model=model, vor=vor, rch_dict=_build_recharge_dict(vor, config.recharge_by_period))
    CHD(model=model, stress_period_data=_build_chd_spd(config))
    GHB(model=model, stress_period_data=_build_ghb_spd(config))

    drn_fields = {
        "name": "name",
        "height_over_btm": "height",
        "conductance": "cond",
        "layer": "layer",
        "min_elev": "min_elev",
    }
    drn_spd = DRNFromVector(
        model=model,
        vor=vor,
        shp_gpkg=drn_source_path,
        uid="name",
        crs=2927,
    ).from_vector(fields=drn_fields, edges_only=False, top_drain=True)
    Drains(model=model, stress_period_data=drn_spd)
    return model


def _build_start_k_array(vor: VoronoiGridPlus, start_k_source: Path, config: SyntheticPestDemoConfig) -> np.ndarray:
    """Map the start K source polygon to the model cells."""

    pseudo_model = SimpleNamespace(vor=vor, nper=1)
    return KFromVector(
        model=pseudo_model,
        shp_gpkg=start_k_source,
        uid="name",
        crs=2927,
    ).from_vector(nlay=1, defaults=[float(config.start_k)])[0]


def _select_observation_cells(config: SyntheticPestDemoConfig) -> list[int]:
    """Choose target cells from normalized x/y fractions across the grid."""

    cells: list[int] = []
    seen: set[int] = set()
    for fx, fy in config.observation_fractions:
        ix = min(config.nx - 1, max(0, int(round(fx * (config.nx - 1)))))
        iy = min(config.ny - 1, max(0, int(round(fy * (config.ny - 1)))))
        cell = iy * config.nx + ix
        if cell not in seen:
            seen.add(cell)
            cells.append(cell)
    return cells


def _build_observation_locations(
    vor: VoronoiGridPlus,
    cells: Sequence[int],
    *,
    crs: str,
) -> gpd.GeoDataFrame:
    centers = vor.gdf_vorPolys.geometry.centroid
    rows = []
    for idx, cell in enumerate(cells, start=1):
        point = Point(float(centers.iloc[int(cell)].x), float(centers.iloc[int(cell)].y))
        rows.append(
            {
                "name": f"OBS_{idx:02d}",
                "layer": 0,
                "group": "synthetic_heads",
                "weight": 1.0,
                "cell": int(cell),
                "geometry": point,
            }
        )
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=crs)


def _extract_head_values(model: SimulationBase, observation_locations: gpd.GeoDataFrame) -> pd.DataFrame:
    """Extract a wide head-target table from a completed MF6 run."""

    head_path = model.model_output_folder_path / f"{model.name}.hds"
    hds = flopy.utils.HeadFile(str(head_path))
    try:
        rows = []
        by_period = {int(kper): (int(kstp), int(kper)) for kstp, kper in hds.get_kstpkper()}
        for per in sorted(by_period):
            layer_data = np.asarray(hds.get_data(kstpkper=by_period[per])[0], dtype=float).reshape(-1)
            row = {"per": int(per)}
            for target in observation_locations.itertuples(index=False):
                row[str(target.name)] = float(layer_data[int(target.cell)])
            rows.append(row)
        return pd.DataFrame(rows)
    finally:
        hds.close()


def _run_model(model: SimulationBase, *, label: str):
    """Run one MF6 model and raise a helpful error on failure."""

    success, buff = model.run_simulation()
    if not success:
        raise RuntimeError(f"{label} failed to run: {buff}")


def _build_targets(truth_model: SimulationBase, config: SyntheticPestDemoConfig) -> tuple[HeadTargets, gpd.GeoDataFrame, pd.DataFrame]:
    """Build pseudo-observation targets from the truth model outputs."""

    observation_cells = _select_observation_cells(config)
    locations = _build_observation_locations(truth_model.vor, observation_cells, crs=config.crs)
    values = _extract_head_values(truth_model, locations)
    targets = HeadTargets(
        locations=locations.loc[:, ["name", "layer", "group", "weight", "geometry"]],
        values=values,
        time_column="per",
    )
    return targets, locations, values


def _build_pest_project(
    *,
    model: SimulationBase,
    pest_workspace: Path,
    start_k_source: Path,
    start_drn_source: Path,
    targets: HeadTargets,
    config: SyntheticPestDemoConfig,
) -> PestProject:
    """Create the synthetic demo PEST project."""

    pest = PestProject(
        model=model,
        name=config.run_family,
        workspace=pest_workspace,
        start_datetime=config.start_date,
    )
    pest.add_parameter(
        KPilotPointParameter(
            name="hk",
            source=VectorParameterSource(
                path=start_k_source,
                value_column="k",
                zone_column="name",
                layer_column="layer",
                crs=2927,
            ),
            parameter_space="multiplier",
            bounds=(0.05, 12.0),
            bounds_mode="multiplier",
            transform="log",
            pp_spacing=float(config.pp_spacing),
            geostruct=ExpGeoStruct(range=max(2.0 * config.cell_size, config.pp_spacing), transform="log"),
        )
    )
    pest.add_parameter(
        DrainConductanceParameter(
            name="drn_cond",
            source=VectorParameterSource(
                path=start_drn_source,
                value_column="cond",
                feature_id_column="par_name",
                layer_column="layer",
                crs=2927,
            ),
            bounds=(0.05, 12.0),
            bounds_mode="multiplier",
            transform="log",
        )
    )
    pest.add_observation(HeadTargetObservationSpec(targets=targets))
    return pest


def run_pestpp_live(pest_workspace: Path, control_file: str, exe_name: str):
    """Run one serial PEST++ calibration and stream output live."""

    print(f"Running PEST++ in: {pest_workspace}")
    print(f"Control file: {control_file}")
    proc = subprocess.Popen(
        [exe_name, control_file],
        cwd=pest_workspace,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="")
    proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, proc.args)


def run_pestpp_parallel(
    *,
    workspace_root: Path,
    pest_workspace: Path,
    control_file: str,
    exe_name: str,
    n_workers: int,
    keep_workers: bool,
) -> Path:
    """Run PEST++ with local workers and return the master directory."""

    import pyemu

    worker_root = workspace_root / "workers"
    worker_root.mkdir(parents=True, exist_ok=True)
    master_dir = workspace_root / "pest_master"
    print(f"Running parallel PEST++ with {n_workers} workers")
    print(f"Template workspace: {pest_workspace}")
    print(f"Master workspace: {master_dir}")
    print(f"Worker root: {worker_root}")
    pyemu.os_utils.start_workers(
        str(pest_workspace),
        exe_name,
        control_file,
        num_workers=int(n_workers),
        worker_root=str(worker_root),
        master_dir=str(master_dir),
        cleanup=not keep_workers,
        verbose=True,
    )
    return master_dir


def _write_review_outputs(
    *,
    workspace_root: Path,
    run_dir: Path,
    truth_k_gdf: gpd.GeoDataFrame,
    truth_drn_gdf: gpd.GeoDataFrame,
    start_drn_gdf: gpd.GeoDataFrame,
) -> dict:
    """Reopen a completed run and write summary review artifacts."""

    review_dir = workspace_root / "review"
    review_dir.mkdir(parents=True, exist_ok=True)

    run = open_pest_run(run_dir)
    review = run.review()
    review.residual_compare.to_csv(review_dir / "residual_compare.csv", index=False)
    review.stats.to_csv(review_dir / "residual_stats.csv", index=False)
    k_review = review.k_geodata
    if k_review.crs is None and truth_k_gdf.crs is not None:
        k_review = k_review.set_crs(truth_k_gdf.crs, allow_override=True)
    k_review.to_file(review_dir / "k_review.gpkg", driver="GPKG")
    truth_k_gdf.to_file(review_dir / "truth_k.gpkg", driver="GPKG")

    final_drain_values = pd.read_csv(run.pest_workspace / "drn_cond_drain_conductance.csv")
    drain_summary = start_drn_gdf.loc[:, ["name", "cond"]].rename(columns={"cond": "start_cond"})
    drain_summary["truth_cond"] = truth_drn_gdf["cond"].to_numpy(dtype=float)
    drain_summary["final_multiplier"] = final_drain_values["value"].to_numpy(dtype=float)
    drain_summary["final_cond"] = drain_summary["start_cond"] * drain_summary["final_multiplier"]
    drain_summary.to_csv(review_dir / "drain_conductance_summary.csv", index=False)

    summary = {
        "result_workspace": str(run_dir),
        "review_dir": str(review_dir),
        "n_targets": int(len(review.targets.to_long())),
        "mae_baseline": float(review.stats["mae_baseline"].iloc[0]),
        "mae_calibrated": float(review.stats["mae_calibrated"].iloc[0]),
        "rmse_baseline": float(review.stats["rmse_baseline"].iloc[0]),
        "rmse_calibrated": float(review.stats["rmse_calibrated"].iloc[0]),
        "k_ratio_min": float(k_review["k_ratio"].min()),
        "k_ratio_max": float(k_review["k_ratio"].max()),
        "drain_final_min": float(drain_summary["final_cond"].min()),
        "drain_final_max": float(drain_summary["final_cond"].max()),
    }
    (review_dir / "review_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def build_and_optionally_run_synthetic_demo(
    config: SyntheticPestDemoConfig,
    *,
    run_pestpp: bool = False,
    noptmax: int = 10,
    n_workers: int = 1,
    keep_workers: bool = False,
    pestpp_exe: str = "pestpp-glm",
) -> SyntheticPestDemoRun:
    """Build the full compact synthetic workflow and optionally run PEST++."""

    workspace_root = build_workspace_root(config.artifact_root, config.run_family)
    inputs_dir = workspace_root / "inputs"
    truth_workspace = workspace_root / "truth"
    model_workspace = workspace_root / "model"
    pest_workspace = workspace_root / "pest"
    review_dir = workspace_root / "review"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    review_dir.mkdir(parents=True, exist_ok=True)

    # Create shared geometry, parameter, and truth-reference inputs.
    shared_vor = _build_rectangular_voronoi(
        nx=config.nx,
        ny=config.ny,
        cell_size=config.cell_size,
        crs=config.crs,
        name=f"{config.model_name}_shared_grid",
    )
    truth_k = _truth_k_field(shared_vor, config)
    truth_k_gdf = shared_vor.gdf_vorPolys.copy()
    truth_k_gdf["cell"] = truth_k_gdf.index.astype(int)
    truth_k_gdf["k_truth"] = truth_k.astype(float)
    truth_k_path = _write_gpkg(inputs_dir / "truth_k.gpkg", truth_k_gdf.loc[:, ["cell", "k_truth", "geometry"]])

    start_k_source_gdf = _build_start_k_source(shared_vor, config)
    start_k_source = _write_gpkg(inputs_dir / "start_k_source.gpkg", start_k_source_gdf)

    truth_drn_gdf, start_drn_gdf = _build_drain_sources(shared_vor, config)
    truth_drn_source = _write_gpkg(inputs_dir / "truth_drains.gpkg", truth_drn_gdf)
    start_drn_source = _write_gpkg(inputs_dir / "start_drains.gpkg", start_drn_gdf)

    # Build and run the truth model.
    truth_model = _build_model(
        workspace=truth_workspace,
        name=f"{config.model_name}_truth",
        config=config,
        k_array=truth_k,
        drn_source_path=truth_drn_source,
    )
    _run_model(truth_model, label="Synthetic truth model")

    targets, observation_locations, truth_target_values = _build_targets(truth_model, config)
    _write_gpkg(inputs_dir / "observation_locations.gpkg", observation_locations)
    truth_target_values.to_csv(inputs_dir / "truth_targets.csv", index=False)

    # Build and run the wrong starting model.
    start_k_array = _build_start_k_array(truth_model.vor, start_k_source, config)
    start_model = _build_model(
        workspace=model_workspace,
        name=config.model_name,
        config=config,
        k_array=start_k_array,
        drn_source_path=start_drn_source,
    )
    start_model.targets.heads = targets
    _run_model(start_model, label="Synthetic starting model")

    baseline_stats = targets.stats(start_model)
    baseline_stats.to_csv(review_dir / "baseline_stats.csv", index=False)

    pest = _build_pest_project(
        model=start_model,
        pest_workspace=pest_workspace,
        start_k_source=start_k_source,
        start_drn_source=start_drn_source,
        targets=targets,
        config=config,
    )
    pst = pest.build_pst("synthetic_compact.pst")
    pst.control_data.noptmax = int(noptmax)
    pst.write(pest_workspace / "synthetic_compact.pst")

    run_info = {
        "run_family": config.run_family,
        "workspace_root": str(workspace_root),
        "inputs_dir": str(inputs_dir),
        "truth_workspace": str(truth_workspace),
        "model_workspace": str(model_workspace),
        "pest_workspace": str(pest_workspace),
        "truth_k_path": str(truth_k_path),
        "start_k_source": str(start_k_source),
        "truth_drn_source": str(truth_drn_source),
        "start_drn_source": str(start_drn_source),
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "noptmax": int(noptmax),
        "n_workers": int(n_workers),
        "run_pestpp": bool(run_pestpp),
    }
    run_info_path = workspace_root / "run_info.json"
    run_info_path.write_text(json.dumps(run_info, indent=2), encoding="utf-8")

    result_workspace: Path | None = None
    review_summary: dict | None = None
    if run_pestpp:
        if int(n_workers) > 1:
            result_workspace = run_pestpp_parallel(
                workspace_root=workspace_root,
                pest_workspace=pest_workspace,
                control_file="synthetic_compact.pst",
                exe_name=pestpp_exe,
                n_workers=int(n_workers),
                keep_workers=keep_workers,
            )
        else:
            run_pestpp_live(pest_workspace, "synthetic_compact.pst", pestpp_exe)
            result_workspace = pest_workspace

        final_par = result_workspace / "synthetic_compact.par"
        if final_par.exists():
            review_summary = _write_review_outputs(
                workspace_root=workspace_root,
                run_dir=result_workspace,
                truth_k_gdf=truth_k_gdf.loc[:, ["cell", "k_truth", "geometry"]],
                truth_drn_gdf=truth_drn_gdf,
                start_drn_gdf=start_drn_gdf,
            )

    return SyntheticPestDemoRun(
        config=config,
        workspace_root=workspace_root,
        inputs_dir=inputs_dir,
        truth_workspace=truth_workspace,
        model_workspace=model_workspace,
        pest_workspace=pest_workspace,
        run_info_path=run_info_path,
        control_file="synthetic_compact.pst",
        targets=targets,
        truth_target_values=truth_target_values,
        observation_locations=observation_locations,
        start_k_source=start_k_source,
        truth_drn_source=truth_drn_source,
        start_drn_source=start_drn_source,
        review_dir=review_dir,
        result_workspace=result_workspace,
        review_summary=review_summary,
    )
