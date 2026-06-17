"""Gold-standard synthetic MF6 + PEST workflow for the preferred API.

This module builds a compact but feature-rich synthetic model that exercises the
preferred ``myflopy`` workflow from end to end:

* DISV/Voronoi groundwater model under 1,000 cells
* UZF infiltration
* two lakes
* one diverted stream network feeding both lakes
* drain seepage zones
* canonical targets for heads, lake stage, SFR stage, SFR flow, and DRN flow
* PEST parameterization using the currently implemented parameter families:
  pilot-point ``K``, drain elevation offsets, and drain-conductance multipliers

The synthetic "truth" model uses a distinct K field plus biased drain inputs so
that calibration movement is visually and numerically obvious.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import subprocess
from typing import Any, Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, Polygon

from myflopy.modflow.mf6.drn import DRN as DRNFromVector
from myflopy.modflow.mf6.kflow import KFromVector
from myflopy.modflow.mf6.lakes import LAKBuilder
from myflopy.modflow.mf6.mvr import MVRBuilder, Move
from myflopy.modflow.mf6.observations import (
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.pest.project import PestProject
from myflopy.modflow.mf6.pest.results import open_pest_run
from myflopy.modflow.mf6.pest.specs import (
    DrainConductanceParameter,
    DrainElevationParameter,
    DrnFlowObservationSpec,
    ExpGeoStruct,
    HeadTargetObservationSpec,
    KPilotPointParameter,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
    VectorParameterSource,
)
from myflopy.modflow.mf6.pest.synthetic_demo import (
    _build_rectangular_voronoi,
    _period_data,
    _run_model,
    _write_gpkg,
    build_workspace_root,
)
from myflopy.modflow.mf6.sfr import SFRBuilder
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization
from myflopy.modflow.mf6.simulation.packages import (
    CHD,
    GHB,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Storage,
)
from myflopy.modflow.mf6.uzf import UZFBuilder
from myflopy.specs import ModelContext


@dataclass(slots=True)
class GoldStandardPestDemoConfig:
    """Configuration for the gold-standard synthetic workflow."""

    artifact_root: Path
    run_family: str = "gold_standard_pest_demo"
    model_name: str = "gold_demo"
    nx: int = 18
    ny: int = 16
    cell_size: float = 120.0
    nper: int = 6
    start_date: str = "2020-01-01"
    west_head: float = 106.0
    east_ghb_cond: float = 450.0
    start_k: float = 12.0
    pp_spacing: float = 360.0
    truth_drain_cond: tuple[float, float, float] = (35.0, 125.0, 260.0)
    start_drain_cond: tuple[float, float, float] = (90.0, 45.0, 70.0)
    truth_drain_height: tuple[float, float, float] = (6.0, 8.0, 5.0)
    start_drain_height: tuple[float, float, float] = (3.0, 4.5, 2.5)
    uzf_infiltration_by_period: tuple[float, ...] = (
        8.0e-5,
        1.7e-4,
        5.0e-5,
        2.2e-4,
        7.5e-5,
        1.9e-4,
    )
    uzf_pet_by_period: tuple[float, ...] = (
        1.2e-4,
        1.5e-4,
        2.3e-4,
        2.1e-4,
        1.6e-4,
        1.3e-4,
    )
    east_stage_by_period: tuple[float, ...] = (
        95.5,
        94.8,
        96.2,
        94.0,
        96.6,
        95.3,
    )
    stream_inflow_by_period: tuple[float, ...] = (
        2.4,
        3.2,
        1.8,
        3.8,
        2.1,
        3.0,
    )
    north_lake_stage_start: float = 103.8
    south_lake_stage_start: float = 101.2
    transient_period_length_days: float = 30.0
    transient_num_steps: int = 1
    crs: str = "EPSG:2927"
    top_base: float = 108.0
    bottom_offset: float = 44.0
    observation_fractions: tuple[tuple[float, float], ...] = (
        (0.12, 0.22),
        (0.25, 0.18),
        (0.39, 0.24),
        (0.56, 0.22),
        (0.72, 0.18),
        (0.88, 0.24),
        (0.15, 0.48),
        (0.30, 0.54),
        (0.44, 0.46),
        (0.58, 0.52),
        (0.74, 0.48),
        (0.90, 0.52),
        (0.18, 0.78),
        (0.35, 0.72),
        (0.52, 0.80),
        (0.68, 0.74),
        (0.84, 0.80),
    )

    def __post_init__(self):
        self.artifact_root = Path(self.artifact_root)
        if self.nx * self.ny > 1000:
            raise ValueError("Gold-standard demo grid must stay at or below 1,000 cells.")
        if self.nper < 2:
            raise ValueError("Gold-standard demo must include at least one steady and one transient period.")
        expected = self.nper
        for name, values in (
            ("uzf_infiltration_by_period", self.uzf_infiltration_by_period),
            ("uzf_pet_by_period", self.uzf_pet_by_period),
            ("east_stage_by_period", self.east_stage_by_period),
            ("stream_inflow_by_period", self.stream_inflow_by_period),
        ):
            if len(values) != expected:
                raise ValueError(f"{name} length must equal nper.")


@dataclass(slots=True)
class GoldStandardModelContext:
    """Package and target-location metadata derived while building one model."""

    lake_polygons: gpd.GeoDataFrame
    stream_paths: list[Path]
    drn_zone_locations: gpd.GeoDataFrame
    lake_targets_empty: LakeStageTargets
    sfr_stage_targets_empty: SfrStageTargets
    sfr_flow_targets_empty: SfrFlowTargets
    drn_targets_empty: DrnFlowTargets
    sfr_builder: SFRBuilder | None = None


@dataclass(slots=True)
class GoldStandardPestDemoRun:
    """Artifacts returned after building or running the gold-standard demo."""

    config: GoldStandardPestDemoConfig
    workspace_root: Path
    inputs_dir: Path
    truth_workspace: Path
    model_workspace: Path
    pest_workspace: Path
    review_dir: Path
    run_info_path: Path
    control_file: str
    targets: dict[str, Any]
    truth_k_path: Path
    start_k_source: Path
    truth_drn_source: Path
    start_drn_source: Path
    result_workspace: Path | None = None
    review_summary: dict[str, Any] | None = None


def _bounds(vor) -> tuple[float, float, float, float]:
    minx, miny, maxx, maxy = vor.gdf_vorPolys.total_bounds
    return float(minx), float(miny), float(maxx), float(maxy)


def _initial_heads(vor, config: GoldStandardPestDemoConfig) -> np.ndarray:
    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    east_start = float(config.east_stage_by_period[0])
    return np.asarray(
        config.west_head - (config.west_head - east_start) * (x / xmax),
        dtype=float,
    )


def _build_surfaces(vor, config: GoldStandardPestDemoConfig) -> tuple[np.ndarray, np.ndarray]:
    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    y = centers.y.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    ymax = max(float(y.max()), 1.0)
    x_norm = x / xmax
    y_norm = y / ymax

    top = (
        config.top_base
        - 7.0 * x_norm
        - 1.6 * np.exp(-((y_norm - 0.72) / 0.10) ** 2)
        - 2.8 * np.exp(-((x_norm - 0.78) / 0.10) ** 2 - ((y_norm - 0.68) / 0.12) ** 2)
        - 1.9 * np.exp(-((x_norm - 0.74) / 0.12) ** 2 - ((y_norm - 0.30) / 0.12) ** 2)
    )
    bottom = top - config.bottom_offset
    return np.asarray(top, dtype=float), np.asarray(bottom, dtype=float)


def _truth_k_field(vor, config: GoldStandardPestDemoConfig) -> np.ndarray:
    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    y = centers.y.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    ymax = max(float(y.max()), 1.0)
    x_norm = x / xmax
    y_norm = y / ymax

    k = np.full(vor.ncpl, float(config.start_k), dtype=float)

    main_corridor = np.exp(-((y_norm - 0.72) / 0.08) ** 2)
    branch_corridor = np.exp(-((y_norm - (0.72 - 0.52 * np.clip((x_norm - 0.45) / 0.40, 0.0, 1.0))) / 0.09) ** 2)
    k *= 1.0 + 4.2 * main_corridor + 2.8 * branch_corridor

    ridge = (
        (x_norm > 0.42)
        & (x_norm < 0.62)
        & (y_norm > 0.18)
        & (y_norm < 0.90)
    )
    k[ridge] *= 0.18

    north_lens = (((x_norm - 0.78) / 0.11) ** 2 + ((y_norm - 0.70) / 0.09) ** 2) < 1.0
    south_lens = (((x_norm - 0.78) / 0.12) ** 2 + ((y_norm - 0.28) / 0.10) ** 2) < 1.0
    k[north_lens] *= 2.6
    k[south_lens] *= 3.2

    return np.clip(k, 0.35, None)


def _build_start_k_source(vor, config: GoldStandardPestDemoConfig) -> gpd.GeoDataFrame:
    minx, miny, maxx, maxy = _bounds(vor)
    polygon = Polygon([(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)])
    return gpd.GeoDataFrame(
        {"name": ["zone_all"], "k": [float(config.start_k)], "layer": [1]},
        geometry=[polygon],
        crs=config.crs,
    )


def _build_lake_polygons(vor, config: GoldStandardPestDemoConfig) -> gpd.GeoDataFrame:
    minx, miny, maxx, maxy = _bounds(vor)
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

    return gpd.GeoDataFrame(
        {"name": ["north_lake", "south_lake"]},
        geometry=[rect(0.72, 0.90, 0.60, 0.82), rect(0.68, 0.88, 0.14, 0.34)],
        crs=config.crs,
    )


def _build_stream_geometries(vor, config: GoldStandardPestDemoConfig) -> list[gpd.GeoDataFrame]:
    minx, miny, maxx, maxy = _bounds(vor)
    width = maxx - minx
    height = maxy - miny

    def xy(fx: float, fy: float) -> tuple[float, float]:
        return (minx + width * fx, miny + height * fy)

    main_line = LineString([xy(0.08, 0.74), xy(0.30, 0.71), xy(0.58, 0.72), xy(0.88, 0.75)])
    return [gpd.GeoDataFrame({"name": ["main_stream"]}, geometry=[main_line], crs=config.crs)]


def _build_drain_sources(vor, config: GoldStandardPestDemoConfig) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    minx, miny, maxx, maxy = _bounds(vor)
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

    names = ["north_relief", "central_valley", "south_relief"]
    polygons = [
        rect(0.58, 0.96, 0.56, 0.78),
        rect(0.38, 0.82, 0.34, 0.56),
        rect(0.56, 0.94, 0.10, 0.30),
    ]
    base = pd.DataFrame(
        {
            "name": names,
            "par_name": names,
            "layer": [1, 1, 1],
            "min_elev": [74.0, 72.0, 71.0],
        }
    )
    truth = gpd.GeoDataFrame(
        base.assign(
            height=list(config.truth_drain_height),
            cond=list(config.truth_drain_cond),
        ),
        geometry=polygons,
        crs=config.crs,
    )
    start = gpd.GeoDataFrame(
        base.assign(
            height=list(config.start_drain_height),
            cond=list(config.start_drain_cond),
        ),
        geometry=polygons,
        crs=config.crs,
    )
    return truth, start


def _build_drn_zone_locations(vor, config: GoldStandardPestDemoConfig) -> gpd.GeoDataFrame:
    minx, miny, maxx, maxy = _bounds(vor)
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

    return gpd.GeoDataFrame(
        {
            "name": ["north_seepage", "central_seepage", "south_seepage"],
            "group": ["drn", "drn", "drn"],
            "weight": [0.75, 0.75, 0.75],
        },
        geometry=[
            rect(0.56, 0.98, 0.54, 0.82),
            rect(0.34, 0.86, 0.32, 0.60),
            rect(0.54, 0.98, 0.08, 0.34),
        ],
        crs=config.crs,
    )


def _build_chd_spd(config: GoldStandardPestDemoConfig) -> dict[int, list[list[object]]]:
    spd: dict[int, list[list[object]]] = {}
    for per in range(config.nper):
        rows = []
        for iy in range(config.ny):
            rows.append([(0, iy * config.nx), float(config.west_head)])
        spd[per] = rows
    return spd


def _build_ghb_spd(config: GoldStandardPestDemoConfig) -> dict[int, list[list[object]]]:
    spd: dict[int, list[list[object]]] = {}
    for per, stage in enumerate(config.east_stage_by_period):
        rows = []
        for iy in range(config.ny):
            rows.append([(0, iy * config.nx + (config.nx - 1)), float(stage), float(config.east_ghb_cond)])
        spd[int(per)] = rows
    return spd


def _build_uzf_finf(vor, uzf_cells: Sequence[tuple[int, int]], config: GoldStandardPestDemoConfig) -> dict[int, list[float]]:
    centers = vor.gdf_vorPolys.geometry.centroid
    x = centers.x.to_numpy(dtype=float)
    y = centers.y.to_numpy(dtype=float)
    xmax = max(float(x.max()), 1.0)
    ymax = max(float(y.max()), 1.0)
    x_norm = x / xmax
    y_norm = y / ymax
    pattern = (
        0.82
        + 0.35 * np.exp(-((y_norm - 0.78) / 0.14) ** 2)
        + 0.22 * np.exp(-((x_norm - 0.24) / 0.18) ** 2)
        + 0.18 * np.exp(-((x_norm - 0.70) / 0.12) ** 2 - ((y_norm - 0.30) / 0.14) ** 2)
    )
    cell_indices = np.asarray([int(cell) for _, cell in uzf_cells], dtype=int)
    return {
        int(per): (base * pattern[cell_indices]).astype(float).tolist()
        for per, base in enumerate(config.uzf_infiltration_by_period)
    }


def _build_uzf_pet(uzf_cells: Sequence[tuple[int, int]], config: GoldStandardPestDemoConfig) -> dict[int, list[float]]:
    count = len(uzf_cells)
    return {int(per): [float(base)] * count for per, base in enumerate(config.uzf_pet_by_period)}


def _select_head_cells(config: GoldStandardPestDemoConfig) -> list[int]:
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


def _build_head_locations(vor, config: GoldStandardPestDemoConfig) -> gpd.GeoDataFrame:
    centers = vor.gdf_vorPolys.geometry.centroid
    rows = []
    for idx, cell in enumerate(_select_head_cells(config), start=1):
        point = Point(float(centers.iloc[int(cell)].x), float(centers.iloc[int(cell)].y))
        rows.append(
            {
                "name": f"GW_{idx:02d}",
                "layer": 0,
                "group": "heads",
                "weight": 1.0,
                "cell": int(cell),
                "geometry": point,
            }
        )
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=config.crs)


def _build_empty_targets(
    *,
    model: SimulationBase,
    sfr: SFRBuilder,
    drn_zone_locations: gpd.GeoDataFrame,
    config: GoldStandardPestDemoConfig,
) -> dict[str, Any]:
    head_locations = _build_head_locations(model.vor, config)
    main_reaches = list(sfr.stream_reaches[sfr.stream_ids[0]])
    sfr_stage_locations = {
        "main_upstream_stage": int(main_reaches[0]),
        "main_split_stage": int(main_reaches[len(main_reaches) // 2]),
        "main_outlet_stage": int(main_reaches[-1]),
    }
    sfr_flow_locations = {
        "main_upstream_flow": int(main_reaches[0]),
        "main_mid_flow": int(main_reaches[len(main_reaches) // 2]),
        "main_outlet_flow": int(main_reaches[-1]),
    }
    return {
        "heads": HeadTargets(
            locations=head_locations.loc[:, ["name", "layer", "group", "weight", "geometry"]],
            values=pd.DataFrame({"time": []}),
            time_column="time",
        ),
        "lake_stage": LakeStageTargets(locations={"north_lake": 0, "south_lake": 1}, values=None),
        "sfr_stage": SfrStageTargets(locations=sfr_stage_locations, values=None),
        "sfr_flow": SfrFlowTargets(locations=sfr_flow_locations, values=None),
        "drn_flow": DrnFlowTargets(locations=drn_zone_locations, values=None),
    }


def _materialize_truth_targets(truth_model: SimulationBase, empty_targets: dict[str, Any]) -> dict[str, Any]:
    heads_empty: HeadTargets = empty_targets["heads"]
    head_values = heads_empty.simulated_heads(truth_model).rename(columns={"per": "time"})
    heads = HeadTargets(
        locations=heads_empty.locations_gdf.copy(),
        values=head_values,
        time_column="time",
    )

    lake_empty: LakeStageTargets = empty_targets["lake_stage"]
    lake_targets = LakeStageTargets(
        locations=lake_empty.get().loc[:, ["name", "lake"]],
        values=lake_empty.simulated_series(truth_model),
        time_column="time",
        value_column="stage",
    )

    sfr_stage_empty: SfrStageTargets = empty_targets["sfr_stage"]
    sfr_stage_targets = SfrStageTargets(
        locations=sfr_stage_empty.get().loc[:, ["name", "reach"]],
        values=sfr_stage_empty.simulated_series(truth_model),
        time_column="time",
        value_column="stage",
    )

    sfr_flow_empty: SfrFlowTargets = empty_targets["sfr_flow"]
    sfr_flow_targets = SfrFlowTargets(
        locations=sfr_flow_empty.get().loc[:, ["name", "reach"]],
        values=sfr_flow_empty.simulated_series(truth_model),
        time_column="time",
        value_column="flow",
    )

    drn_empty: DrnFlowTargets = empty_targets["drn_flow"]
    drn_flow_targets = DrnFlowTargets(
        locations=drn_empty.get(model=truth_model),
        values=drn_empty.simulated_series(truth_model),
        time_column="time",
        value_column="flow",
    )
    return {
        "heads": heads,
        "lake_stage": lake_targets,
        "sfr_stage": sfr_stage_targets,
        "sfr_flow": sfr_flow_targets,
        "drn_flow": drn_flow_targets,
    }


def _attach_targets_and_obs(model: SimulationBase, targets: dict[str, Any]):
    import flopy

    model.targets.heads = targets["heads"]
    model.targets.lake_stage = targets["lake_stage"]
    model.targets.sfr_stage = targets["sfr_stage"]
    model.targets.sfr_flow = targets["sfr_flow"]
    model.targets.drn_flow = targets["drn_flow"]

    model.targets.heads.attach_flopy_obs(filename="gold_heads.obs", csv_name="gold_heads.csv")
    model.targets.lake_stage.attach_flopy_obs(filename="gold_lakes.obs", csv_name="gold_lakes.csv")
    sfr_continuous = {}
    sfr_continuous.update(model.targets.sfr_stage.targets.to_flopy_obs(csv_name="gold_sfr_stage.csv"))
    flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
        model.gwf.sfr,
        pname="sfr_obs",
        continuous=sfr_continuous,
        filename="gold_sfr.obs",
    )
    model.targets.drn_flow.attach_flopy_obs(filename="gold_drn.obs", csv_name="gold_drn.csv")


def _build_model(
    *,
    workspace: Path,
    inputs_dir: Path,
    name: str,
    config: GoldStandardPestDemoConfig,
    k_array: np.ndarray,
    drn_source_path: Path,
    lake_paths: list[Path],
    stream_paths: list[Path],
) -> tuple[SimulationBase, GoldStandardModelContext]:
    vor = _build_rectangular_voronoi(
        nx=config.nx,
        ny=config.ny,
        cell_size=config.cell_size,
        crs=config.crs,
        name=f"{name}_grid",
    )
    top, bottom = _build_surfaces(vor, config)
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {"geometry": vor.gdf_vorPolys.geometry, 0: top.astype(float), 1: bottom.astype(float)},
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
    KFlow(
        model=model,
        k=np.asarray(k_array, dtype=float).tolist(),
        k33_vert=(np.asarray(k_array, dtype=float) / 10.0).tolist(),
    )
    Storage(
        model=model,
        specific_storage=1.5e-4,
        specific_yield=0.18,
        sto_steady={0: True},
        sto_transient={per: True for per in range(1, config.nper)},
    )
    CHD(model=model, stress_period_data=_build_chd_spd(config))
    GHB(model=model, stress_period_data=_build_ghb_spd(config))

    lak = LAKBuilder(
        context=ModelContext(grid=vor, domain=np.ones((1, int(vor.ncpl)), dtype=int)),
        nper=config.nper,
        lakes=lake_paths,
        lake_id_field="name",
        starting_stage={
            "north_lake": config.north_lake_stage_start,
            "south_lake": config.south_lake_stage_start,
        },
        lake_bottom={
            "north_lake": config.north_lake_stage_start - 7.5,
            "south_lake": config.south_lake_stage_start - 7.0,
        },
        bed_leakance=0.12,
        connection_modes="rectangular",
        status={"north_lake": ["ACTIVE"] * config.nper, "south_lake": ["ACTIVE"] * config.nper},
        mover=True,
        length_conversion=3.28081,
        time_conversion=86_400.0,
    )
    lak.build().build(model.gwf)
    for lake_id, cells in lak.lake_cells.items():
        model.add_region_from_cells(
            f"lake_zone_{lake_id}",
            [(0, cell) for cell in cells],
            category="boundary",
            package="lak",
            tags=["lak"],
            geometry=lak.lake_table.loc[lake_id, "geometry"],
            overwrite=True,
        )
    model.add_region_from_cells(
        "all_lakes",
        [(0, cell) for cells in lak.lake_cells.values() for cell in cells],
        category="boundary",
        package="lak",
        tags=["lak"],
        overwrite=True,
    )

    inflows = {int(per): [(0, float(flow))] for per, flow in enumerate(config.stream_inflow_by_period)}
    sfr = SFRBuilder(
        context=ModelContext(grid=vor, domain=np.ones((1, int(vor.ncpl)), dtype=int)),
        nper=config.nper,
        streams=stream_paths,
        inflow=inflows,
        width=18.0,
        gradient=0.0008,
        roughness=0.030,
        streambed_k=1.0e-5,
        streambed_thickness=2.0,
        mover=True,
    )
    sfr.build().build(model.gwf)
    model.add_region_from_cells(
        "all_streams",
        [(0, cell) for cells in sfr.stream_cells.values() for cell in cells],
        category="boundary",
        package="sfr",
        tags=["sfr"],
        overwrite=True,
    )

    stream_id = sfr.stream_ids[0]
    mvr = MVRBuilder(
        nper=config.nper,
        moves={
            period: [
                Move(sfr.connection(stream_id), lak.connection("north_lake"), value=0.60),
                Move(sfr.connection(stream_id), lak.connection("south_lake"), value=0.40),
            ]
            for period in range(config.nper)
        },
    )
    mvr.build().build(model.gwf)

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

    lake_cells = {
        int(cell)
        for values in lak.lake_cells.values()
        for cell in values
    }
    stream_cells = {int(cell) for cells in sfr.stream_cells.values() for cell in cells}
    uzf_cells = [(0, int(cell)) for cell in range(vor.ncpl) if cell not in lake_cells and cell not in stream_cells]
    uzf_builder = UZFBuilder(
        context=ModelContext(grid=vor, domain=np.ones((1, int(vor.ncpl)), dtype=int)),
        nper=config.nper,
        cells=uzf_cells,
        vks=0.35,
        thtr=0.08,
        thts=0.33,
        thti=0.15,
        eps=4.2,
        finf=_build_uzf_finf(vor, uzf_cells, config),
        pet=_build_uzf_pet(uzf_cells, config),
        extdp=6.0,
    )
    uzf = uzf_builder.build().build(model.gwf)
    model.add_region_from_cells(
        "uzf_active",
        uzf_builder.uzf_cells,
        category="boundary",
        package="uzf",
        tags=["uzf"],
        metadata={"nuzfcells": len(uzf_builder.uzf_cells)},
        overwrite=True,
    )

    drn_zone_locations = _build_drn_zone_locations(vor, config)
    context = GoldStandardModelContext(
        lake_polygons=pd.concat([gpd.read_file(path) for path in lake_paths], ignore_index=True),
        stream_paths=stream_paths,
        drn_zone_locations=drn_zone_locations,
        lake_targets_empty=LakeStageTargets(locations={"north_lake": 0, "south_lake": 1}),
        sfr_stage_targets_empty=SfrStageTargets(locations={}),
        sfr_flow_targets_empty=SfrFlowTargets(locations={}),
        drn_targets_empty=DrnFlowTargets(locations=drn_zone_locations),
        sfr_builder=sfr,
    )
    empty_targets = _build_empty_targets(model=model, sfr=sfr, drn_zone_locations=drn_zone_locations, config=config)
    context.lake_targets_empty = empty_targets["lake_stage"]
    context.sfr_stage_targets_empty = empty_targets["sfr_stage"]
    context.sfr_flow_targets_empty = empty_targets["sfr_flow"]
    context.drn_targets_empty = empty_targets["drn_flow"]
    return model, context


def _build_pest_project(
    *,
    model: SimulationBase,
    pest_workspace: Path,
    start_k_source: Path,
    start_drn_source: Path,
    targets: dict[str, Any],
    config: GoldStandardPestDemoConfig,
) -> PestProject:
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
            bounds=(0.05, 15.0),
            bounds_mode="multiplier",
            transform="log",
            pp_spacing=float(config.pp_spacing),
            geostruct=ExpGeoStruct(range=max(2.0 * config.cell_size, config.pp_spacing), transform="log"),
        )
    )
    pest.add_parameter(
        DrainElevationParameter(
            name="drn_elev",
            source=VectorParameterSource(
                path=start_drn_source,
                value_column="height",
                feature_id_column="par_name",
                layer_column="layer",
                crs=2927,
            ),
            bounds=(-8.0, 8.0),
            bounds_mode="absolute",
            transform="none",
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

    pest.add_observation(HeadTargetObservationSpec(targets=targets["heads"]))
    pest.add_observation(LakeStageObservationSpec(targets=targets["lake_stage"], prefix="lak_stage"))
    pest.add_observation(SfrStageObservationSpec(targets=targets["sfr_stage"], prefix="sfr_stage"))
    pest.add_observation(SfrFlowObservationSpec(targets=targets["sfr_flow"], prefix="sfr_flow"))
    pest.add_observation(DrnFlowObservationSpec(targets=targets["drn_flow"], prefix="drn_flow"))
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


def _merge_baseline_calibrated(
    *,
    baseline: pd.DataFrame,
    calibrated: pd.DataFrame,
    join_columns: list[str],
    sim_column: str,
    target_column: str,
) -> pd.DataFrame:
    base = baseline.copy().rename(
        columns={
            sim_column: f"{sim_column}_baseline",
            "residual": "residual_baseline",
            "abs_residual": "abs_residual_baseline",
        }
    )
    calib = calibrated.copy().rename(
        columns={
            sim_column: f"{sim_column}_calibrated",
            "residual": "residual_calibrated",
            "abs_residual": "abs_residual_calibrated",
        }
    )
    keep_cols = join_columns + [
        target_column,
        f"{sim_column}_calibrated",
        "residual_calibrated",
        "abs_residual_calibrated",
    ]
    merged = base.merge(calib.loc[:, keep_cols], on=join_columns + [target_column], how="inner")
    merged["abs_residual_improvement"] = merged["abs_residual_baseline"] - merged["abs_residual_calibrated"]
    return merged


def _dataset_review(
    *,
    label: str,
    targets,
    baseline_model,
    calibrated_model,
    compare_path: Path,
    sim_column: str,
    target_column: str,
    join_columns: list[str],
) -> pd.DataFrame:
    baseline = targets.compare(baseline_model)
    calibrated = targets.compare(calibrated_model)
    merged = _merge_baseline_calibrated(
        baseline=baseline,
        calibrated=calibrated,
        join_columns=join_columns,
        sim_column=sim_column,
        target_column=target_column,
    )
    merged.to_csv(compare_path, index=False)
    stats = pd.DataFrame(
        {
            "dataset": [label],
            "mae_baseline": [float(baseline["abs_residual"].mean())],
            "mae_calibrated": [float(calibrated["abs_residual"].mean())],
            "rmse_baseline": [float(np.sqrt((baseline["residual"] ** 2).mean()))],
            "rmse_calibrated": [float(np.sqrt((calibrated["residual"] ** 2).mean()))],
        }
    )
    return stats


def _write_review_outputs(
    *,
    workspace_root: Path,
    run_dir: Path,
    truth_k_gdf: gpd.GeoDataFrame,
    truth_drn_gdf: gpd.GeoDataFrame,
    start_drn_gdf: gpd.GeoDataFrame,
) -> dict[str, Any]:
    review_dir = workspace_root / "review"
    review_dir.mkdir(parents=True, exist_ok=True)

    run = open_pest_run(run_dir)
    run.export_review(review_dir / "head_review", timeseries_names=["GW_01", "GW_08", "GW_14"])

    baseline_model = run.load_baseline_model()
    calibrated_model = run.load_calibrated_model()

    targets = {
        "heads": run.load_head_targets(),
        "lake_stage": run.load_lake_stage_targets(),
        "sfr_stage": run.load_sfr_stage_targets(),
        "sfr_flow": run.load_sfr_flow_targets(),
        "drn_flow": run.load_drn_flow_targets(),
    }

    stats_frames = [
        _dataset_review(
            label="heads",
            targets=targets["heads"],
            baseline_model=baseline_model,
            calibrated_model=calibrated_model,
            compare_path=review_dir / "heads_compare.csv",
            sim_column="sim_head",
            target_column="head_target",
            join_columns=["name", "layer", "time", "per"],
        ),
        _dataset_review(
            label="lake_stage",
            targets=targets["lake_stage"],
            baseline_model=baseline_model,
            calibrated_model=calibrated_model,
            compare_path=review_dir / "lake_stage_compare.csv",
            sim_column="sim_stage",
            target_column="stage_target",
            join_columns=["name", "lake", "time", "per"],
        ),
        _dataset_review(
            label="sfr_stage",
            targets=targets["sfr_stage"],
            baseline_model=baseline_model,
            calibrated_model=calibrated_model,
            compare_path=review_dir / "sfr_stage_compare.csv",
            sim_column="sim_stage",
            target_column="stage_target",
            join_columns=["name", "reach", "time", "per"],
        ),
        _dataset_review(
            label="sfr_flow",
            targets=targets["sfr_flow"],
            baseline_model=baseline_model,
            calibrated_model=calibrated_model,
            compare_path=review_dir / "sfr_flow_compare.csv",
            sim_column="sim_flow",
            target_column="flow_target",
            join_columns=["name", "reach", "time", "per"],
        ),
        _dataset_review(
            label="drn_flow",
            targets=targets["drn_flow"],
            baseline_model=baseline_model,
            calibrated_model=calibrated_model,
            compare_path=review_dir / "drn_flow_compare.csv",
            sim_column="sim_flow",
            target_column="flow_target",
            join_columns=["name", "time", "per"],
        ),
    ]
    stats = pd.concat(stats_frames, ignore_index=True)
    stats["mae_improvement"] = stats["mae_baseline"] - stats["mae_calibrated"]
    stats["rmse_improvement"] = stats["rmse_baseline"] - stats["rmse_calibrated"]
    stats.to_csv(review_dir / "target_stats.csv", index=False)

    truth_k_out = review_dir / "truth_k.gpkg"
    if truth_k_out.exists():
        truth_k_out.unlink()
    truth_k_gdf.to_file(truth_k_out, driver="GPKG")

    final_cond = pd.read_csv(run.pest_workspace / "drn_cond_drain_conductance.csv")
    final_elev = pd.read_csv(run.pest_workspace / "drn_elev_drain_elevation.csv")
    drain_summary = start_drn_gdf.loc[:, ["name", "height", "cond"]].rename(
        columns={"height": "start_height", "cond": "start_cond"}
    )
    drain_summary["truth_height"] = truth_drn_gdf["height"].to_numpy(dtype=float)
    drain_summary["truth_cond"] = truth_drn_gdf["cond"].to_numpy(dtype=float)
    drain_summary["final_height_offset"] = final_elev["value"].to_numpy(dtype=float)
    drain_summary["final_cond_multiplier"] = final_cond["value"].to_numpy(dtype=float)
    drain_summary["final_height"] = drain_summary["start_height"] + drain_summary["final_height_offset"]
    drain_summary["final_cond"] = drain_summary["start_cond"] * drain_summary["final_cond_multiplier"]
    drain_summary.to_csv(review_dir / "drain_parameter_summary.csv", index=False)

    summary = {
        "result_workspace": str(run_dir),
        "review_dir": str(review_dir),
        "datasets": stats["dataset"].tolist(),
        "head_rmse_baseline": float(stats.loc[stats["dataset"] == "heads", "rmse_baseline"].iloc[0]),
        "head_rmse_calibrated": float(stats.loc[stats["dataset"] == "heads", "rmse_calibrated"].iloc[0]),
        "lake_rmse_baseline": float(stats.loc[stats["dataset"] == "lake_stage", "rmse_baseline"].iloc[0]),
        "lake_rmse_calibrated": float(stats.loc[stats["dataset"] == "lake_stage", "rmse_calibrated"].iloc[0]),
        "sfr_flow_rmse_baseline": float(stats.loc[stats["dataset"] == "sfr_flow", "rmse_baseline"].iloc[0]),
        "sfr_flow_rmse_calibrated": float(stats.loc[stats["dataset"] == "sfr_flow", "rmse_calibrated"].iloc[0]),
        "drn_flow_rmse_baseline": float(stats.loc[stats["dataset"] == "drn_flow", "rmse_baseline"].iloc[0]),
        "drn_flow_rmse_calibrated": float(stats.loc[stats["dataset"] == "drn_flow", "rmse_calibrated"].iloc[0]),
        "k_ratio_min": float(run.k_geodata()["k_ratio"].min()),
        "k_ratio_max": float(run.k_geodata()["k_ratio"].max()),
    }
    (review_dir / "review_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def build_and_optionally_run_gold_standard_demo(
    config: GoldStandardPestDemoConfig,
    *,
    run_pestpp: bool = False,
    noptmax: int = 10,
    n_workers: int = 1,
    keep_workers: bool = False,
    pestpp_exe: str = "pestpp-glm",
) -> GoldStandardPestDemoRun:
    """Build the full feature-rich synthetic workflow and optionally run PEST++."""

    workspace_root = build_workspace_root(config.artifact_root, config.run_family)
    inputs_dir = workspace_root / "inputs"
    truth_workspace = workspace_root / "truth"
    model_workspace = workspace_root / "model"
    pest_workspace = workspace_root / "pest"
    review_dir = workspace_root / "review"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    review_dir.mkdir(parents=True, exist_ok=True)

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

    start_k_source = _write_gpkg(inputs_dir / "start_k_source.gpkg", _build_start_k_source(shared_vor, config))
    truth_drn_gdf, start_drn_gdf = _build_drain_sources(shared_vor, config)
    truth_drn_source = _write_gpkg(inputs_dir / "truth_drains.gpkg", truth_drn_gdf)
    start_drn_source = _write_gpkg(inputs_dir / "start_drains.gpkg", start_drn_gdf)

    lake_polygons = _build_lake_polygons(shared_vor, config)
    lake_path = _write_gpkg(inputs_dir / "lakes.gpkg", lake_polygons)
    lake_paths: list[Path] = []
    for row in lake_polygons.itertuples(index=False):
        one = gpd.GeoDataFrame({"name": [row.name]}, geometry=[row.geometry], crs=lake_polygons.crs)
        one_path = inputs_dir / f"{row.name}.gpkg"
        _write_gpkg(one_path, one)
        lake_paths.append(one_path)
    stream_geometries = _build_stream_geometries(shared_vor, config)
    stream_paths: list[Path] = []
    for idx, gdf in enumerate(stream_geometries, start=1):
        path = inputs_dir / f"stream_{idx}.gpkg"
        _write_gpkg(path, gdf)
        stream_paths.append(path)

    truth_model, truth_context = _build_model(
        workspace=truth_workspace,
        inputs_dir=inputs_dir,
        name=f"{config.model_name}_truth",
        config=config,
        k_array=truth_k,
        drn_source_path=truth_drn_source,
        lake_paths=lake_paths,
        stream_paths=stream_paths,
    )
    empty_targets = _build_empty_targets(
        model=truth_model,
        sfr=truth_context.sfr_builder,
        drn_zone_locations=truth_context.drn_zone_locations,
        config=config,
    )
    _attach_targets_and_obs(truth_model, empty_targets)
    _run_model(truth_model, label="Gold-standard truth model")

    targets = _materialize_truth_targets(truth_model, empty_targets)
    for name, target in targets.items():
        if hasattr(target, "to_wide"):
            target.to_wide().to_csv(inputs_dir / f"{name}_truth_targets.csv", index=False)
    if isinstance(targets["heads"], HeadTargets):
        _write_gpkg(inputs_dir / "head_locations.gpkg", targets["heads"].locations_gdf.copy())
    if isinstance(truth_context.drn_zone_locations, gpd.GeoDataFrame):
        _write_gpkg(inputs_dir / "drn_target_zones.gpkg", truth_context.drn_zone_locations.copy())

    start_k_array = KFromVector(
        model=truth_model,
        shp_gpkg=start_k_source,
        uid="name",
        crs=2927,
    ).from_vector(nlay=1, defaults=[float(config.start_k)])[0]
    start_model, _ = _build_model(
        workspace=model_workspace,
        inputs_dir=inputs_dir,
        name=config.model_name,
        config=config,
        k_array=start_k_array,
        drn_source_path=start_drn_source,
        lake_paths=lake_paths,
        stream_paths=stream_paths,
    )
    _attach_targets_and_obs(start_model, targets)
    _run_model(start_model, label="Gold-standard starting model")

    baseline_stats = pd.concat(
        [
            targets["heads"].stats(start_model).assign(dataset="heads"),
            targets["lake_stage"].stats(start_model).assign(dataset="lake_stage"),
            targets["sfr_stage"].stats(start_model).assign(dataset="sfr_stage"),
            targets["sfr_flow"].stats(start_model).assign(dataset="sfr_flow"),
            targets["drn_flow"].stats(start_model).assign(dataset="drn_flow"),
        ],
        ignore_index=True,
    )
    baseline_stats.to_csv(review_dir / "baseline_target_stats.csv", index=False)

    pest = _build_pest_project(
        model=start_model,
        pest_workspace=pest_workspace,
        start_k_source=start_k_source,
        start_drn_source=start_drn_source,
        targets=targets,
        config=config,
    )
    pst = pest.build_pst("gold_standard.pst")
    pst.control_data.noptmax = int(noptmax)
    pst.write(pest_workspace / "gold_standard.pst")

    run_info = {
        "run_family": config.run_family,
        "workspace_root": str(workspace_root),
        "truth_workspace": str(truth_workspace),
        "model_workspace": str(model_workspace),
        "pest_workspace": str(pest_workspace),
        "truth_k_path": str(truth_k_path),
        "start_k_source": str(start_k_source),
        "truth_drn_source": str(truth_drn_source),
        "start_drn_source": str(start_drn_source),
        "lake_path": str(lake_path),
        "lake_paths": [str(path) for path in lake_paths],
        "stream_paths": [str(path) for path in stream_paths],
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "noptmax": int(noptmax),
        "n_workers": int(n_workers),
        "run_pestpp": bool(run_pestpp),
        "implemented_parameter_families": ["k_pilotpoints", "drn_elev", "drn_cond"],
        "implemented_target_families": ["heads", "lake_stage", "sfr_stage", "sfr_flow", "drn_flow"],
    }
    run_info_path = workspace_root / "run_info.json"
    run_info_path.write_text(json.dumps(run_info, indent=2), encoding="utf-8")

    result_workspace: Path | None = None
    review_summary: dict[str, Any] | None = None
    if run_pestpp:
        if int(n_workers) > 1:
            result_workspace = run_pestpp_parallel(
                workspace_root=workspace_root,
                pest_workspace=pest_workspace,
                control_file="gold_standard.pst",
                exe_name=pestpp_exe,
                n_workers=int(n_workers),
                keep_workers=keep_workers,
            )
        else:
            run_pestpp_live(pest_workspace, "gold_standard.pst", pestpp_exe)
            result_workspace = pest_workspace

        if (result_workspace / "gold_standard.par").exists():
            review_summary = _write_review_outputs(
                workspace_root=workspace_root,
                run_dir=result_workspace,
                truth_k_gdf=truth_k_gdf.loc[:, ["cell", "k_truth", "geometry"]],
                truth_drn_gdf=truth_drn_gdf,
                start_drn_gdf=start_drn_gdf,
            )

    return GoldStandardPestDemoRun(
        config=config,
        workspace_root=workspace_root,
        inputs_dir=inputs_dir,
        truth_workspace=truth_workspace,
        model_workspace=model_workspace,
        pest_workspace=pest_workspace,
        review_dir=review_dir,
        run_info_path=run_info_path,
        control_file="gold_standard.pst",
        targets=targets,
        truth_k_path=truth_k_path,
        start_k_source=start_k_source,
        truth_drn_source=truth_drn_source,
        start_drn_source=start_drn_source,
        result_workspace=result_workspace,
        review_summary=review_summary,
    )
