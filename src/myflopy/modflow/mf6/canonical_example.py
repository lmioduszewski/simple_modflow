"""Build the authoritative package-rich canonical myflopy model."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Polygon

from myflopy.modflow.mf6.lakes import LAKBuilder
from myflopy.modflow.mf6.mvr import MVRBuilder, Move
from myflopy.modflow.mf6.canonical import (
    CANONICAL_MODEL_CONTRACT,
    irregular_voronoi_grid,
)
from myflopy.modflow.mf6.observations import (
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.sfr import SFRBuilder
from myflopy.modflow.mf6.uzf import UZFBuilder
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.discretization import (
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.packages import (
    CHD,
    GHB,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
    Wells,
)
from myflopy.specs import ModelContext


@dataclass(frozen=True)
class CanonicalModelConfig:
    """Size and timing controls for the master example model."""

    nrow: int = 100
    ncol: int = 100
    nlay: int = 4
    nper: int = 6
    cell_size: float = 100.0
    period_length: float = 30.0
    steps_per_period: int = 2

    @classmethod
    def validation(cls) -> "CanonicalModelConfig":
        """Return a smaller profile for automated API validation."""

        return cls(nrow=50, ncol=50, nlay=4, nper=6, cell_size=100.0)

    @property
    def ncpl(self) -> int:
        return self.nrow * self.ncol


def rectangular_voronoi(config: CanonicalModelConfig) -> VoronoiGridPlus:
    """Build the deterministic irregular Voronoi grid shared by all profiles."""

    return irregular_voronoi_grid(
        nrow=config.nrow,
        ncol=config.ncol,
        cell_size=config.cell_size,
    )


def build_transient_model(
    workspace: str | Path,
    *,
    config: CanonicalModelConfig | None = None,
    name: str = "viz_prt_master",
) -> SimulationBase:
    """Build a package-rich transient multilayer GWF model for result review."""

    config = CanonicalModelConfig() if config is None else config
    workspace = Path(workspace)
    vor = rectangular_voronoi(config)
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=config.nper)

    centers = np.asarray(vor.points, dtype=float)
    width = config.ncol * config.cell_size
    height = config.nrow * config.cell_size
    x = centers[:, 0] / width
    y = centers[:, 1] / height
    top = 137.0 - 24.0 * x + 3.5 * np.sin(np.pi * y)
    bottom = np.vstack(
        [
            top - 18.0,
            top - 40.0,
            top - 49.0,
            top - 82.0,
        ]
    )
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {
            "geometry": vor.gdf_vorPolys.geometry,
            0: top,
            **{layer + 1: bottom[layer] for layer in range(config.nlay)},
        },
        geometry="geometry",
        crs=vor.crs,
    )
    idomain = np.ones((config.nlay, config.ncpl), dtype=int)
    DisvGrid(vor=vor, model=model, top=top, bottom=bottom, nlay=config.nlay, idomain=idomain)
    TemporalDiscretization(
        model=model,
        period_data=[
            (config.period_length, config.steps_per_period, 1.0)
            for _ in range(config.nper)
        ],
    )

    initial = np.vstack(
        [
            128.0 - 23.0 * x - 0.55 * layer + 1.2 * np.sin(np.pi * y)
            for layer in range(config.nlay)
        ]
    )
    paleochannel = 115.0 * np.exp(-((y - (0.70 - 0.12 * x)) / 0.075) ** 2)
    low_k_barrier = np.where((x > 0.44) & (x < 0.51) & (y < 0.62), 0.18, 1.0)
    aquifer_pattern = (38.0 + 34.0 * x + paleochannel) * low_k_barrier
    conductivity = np.vstack(
        [
            aquifer_pattern,
            aquifer_pattern * 0.75,
            np.maximum(aquifer_pattern * 0.0015, 0.03),
            aquifer_pattern * 0.50,
        ]
    )
    InitialConditions(model=model, vor=vor, nlay=config.nlay, strt=initial)
    KFlow(
        model=model,
        k=conductivity,
        k33_vert=np.maximum(conductivity * 0.05, 0.002),
        icelltype=list(CANONICAL_MODEL_CONTRACT.icelltype),
        save_specific_discharge=True,
    )
    Storage(
        model=model,
        specific_storage=2.0e-5,
        specific_yield=0.18,
        sto_steady={0: True},
        sto_transient={period: True for period in range(1, config.nper)},
        iconvert=list(CANONICAL_MODEL_CONTRACT.icelltype),
    )

    left_cells = [row * config.ncol for row in range(config.nrow)]
    right_cells = [row * config.ncol + config.ncol - 1 for row in range(config.nrow)]
    stress_period_data = {}
    for period in range(config.nper):
        seasonal = 2.2 * np.sin(period * 2.0 * np.pi / max(config.nper - 1, 1))
        records = []
        for layer in range(config.nlay):
            records.extend(
                [
                    (
                        (layer, cell),
                        max(
                            float(bottom[layer, cell] + 0.5),
                            float(129.0 + seasonal - 0.45 * layer + 1.5 * y[cell]),
                        ),
                    )
                    for cell in left_cells
                ]
            )
        stress_period_data[period] = records
    CHD(model=model, stress_period_data=stress_period_data)

    ghb_data = {}
    for period in range(config.nper):
        ghb_data[period] = [
            (
                (layer, cell),
                max(
                    float(bottom[layer, cell] + 1.0),
                    float(104.0 - 0.7 * layer + 0.8 * np.cos(period * 2.0 * np.pi / config.nper)),
                ),
                750.0,
            )
            for layer in range(config.nlay)
            for cell in right_cells
        ]
    GHB(model=model, stress_period_data=ghb_data)

    mountain_front_cells = [
        row * config.ncol + col
        for row in range(config.nrow)
        for col in range(max(1, config.ncol // 12))
    ]
    mountain_front_recharge = {
        period: [
            ((0, cell), float(2.2e-4 + 1.2e-4 * np.sin(period * 2.0 * np.pi / config.nper)))
            for cell in mountain_front_cells
        ]
        for period in range(config.nper)
    }
    Recharge(model=model, rch_dict=mountain_front_recharge)

    unconfined_well = (config.nrow // 3) * config.ncol + (2 * config.ncol // 3)
    confined_well = (2 * config.nrow // 3) * config.ncol + (3 * config.ncol // 4)
    well_data = {
        period: [
            (
                (1, unconfined_well),
                -2200.0 - (4200.0 if period in {2, 3} else 0.0),
                "unconfined_supply",
            ),
            (
                (3, confined_well),
                -3200.0 - (5600.0 if period in {3, 4} else 0.0),
                "confined_supply",
            ),
        ]
        for period in range(config.nper)
    }
    Wells(model=model, stress_period_data=well_data)

    unconfined_seepage_cells = np.flatnonzero(
        (x > 0.48) & (x < 0.66) & (y > 0.66) & (y < 0.88)
    )[::2].tolist()
    confined_seepage_cells = np.flatnonzero(
        (x > 0.72) & (x < 0.91) & (y > 0.08) & (y < 0.30)
    )[::2].tolist()
    drain_cells = [*unconfined_seepage_cells, *confined_seepage_cells]
    drain_data = {
        period: (
            [
                ((0, cell), float(top[cell] - 12.0), 70.0 + 5.0 * (cell % 7))
                for cell in unconfined_seepage_cells
            ]
            + [
                ((3, cell), float(top[cell] - 21.0), 12.0 + 1.5 * (cell % 5))
                for cell in confined_seepage_cells
            ]
        )
        for period in range(config.nper)
    }
    Drains(model=model, stress_period_data=drain_data)
    model.add_region_from_cells(
        "unconfined_seepage_slope",
        [(0, cell) for cell in unconfined_seepage_cells],
        category="seepage",
        package="drn",
        overwrite=True,
    )
    model.add_region_from_cells(
        "confined_seepage_slope",
        [(3, cell) for cell in confined_seepage_cells],
        category="seepage",
        package="drn",
        overwrite=True,
    )

    inputs = _write_surface_water_inputs(workspace, config)
    lak = LAKBuilder(
        context=ModelContext(grid=vor, domain=idomain),
        nper=config.nper,
        lakes=inputs["lakes"],
        lake_id_field="name",
        starting_stage={"north_lake": 120.0, "south_lake": 118.0},
        lake_bottom={"north_lake": 108.0, "south_lake": 106.0},
        bed_leakance=0.003,
        connection_modes="automatic",
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

    inflows = {
        period: [(0, 90_000.0 + 15_000.0 * np.sin(period * 2.0 * np.pi / max(config.nper, 1)))]
        for period in range(config.nper)
    }
    sfr = SFRBuilder(
        context=ModelContext(grid=vor, domain=idomain),
        nper=config.nper,
        streams=[inputs["stream"]],
        inflow=inflows,
        width=28.0,
        gradient=0.0015,
        roughness=0.028,
        streambed_k=0.01,
        streambed_thickness=1.5,
        length_conversion=3.28081,
        time_conversion=86_400.0,
        mover=True,
    )
    stream_id = sfr.stream_ids[0]
    reach_cells = sfr.stream_cells[stream_id]
    reach_fraction = np.linspace(0.0, 1.0, len(reach_cells))
    gw_reference = initial[0, reach_cells]
    exchange_offset = 1.3 * np.sin(5.0 * np.pi * reach_fraction)
    reach_elevations = np.minimum.accumulate(gw_reference + exchange_offset)
    sfr = sfr.with_updates(
        reach_top=reach_elevations.tolist(),
        width=[24.0 + 10.0 * value for value in reach_fraction],
        streambed_k=[0.004 + 0.014 * abs(np.sin(3.0 * np.pi * value)) for value in reach_fraction],
    )
    sfr.build().build(model.gwf)
    model.add_region_from_cells(
        "all_streams",
        [cellid for cellids in sfr.stream_cells.values() for cellid in [(0, cell) for cell in cellids]],
        category="boundary",
        package="sfr",
        tags=["sfr"],
        overwrite=True,
    )
    mvr = MVRBuilder(
        nper=config.nper,
        moves={
            period: [
                Move(sfr.connection(stream_id), lak.connection("north_lake"), value=0.20),
                Move(sfr.connection(stream_id), lak.connection("south_lake"), value=0.15),
            ]
            for period in range(config.nper)
        },
    )
    mvr.build().build(model.gwf)

    lake_cells = {
        int(cell)
        for values in lak.lake_cells.values()
        for cell in values
    }
    stream_cells = {int(cell) for cells in sfr.stream_cells.values() for cell in cells}
    uzf_cells = [
        (0, cell)
        for cell in range(config.ncpl)
        if cell not in lake_cells and cell not in stream_cells
    ]
    infiltration_pond_cells = [
        cell
        for _layer, cell in uzf_cells
        if ((x[cell] - 0.31) / 0.075) ** 2 + ((y[cell] - 0.28) / 0.065) ** 2 <= 1.0
    ]
    pond_cell_set = set(infiltration_pond_cells)
    finf = {
        period: [
            float(
                (
                    9.5e-4
                    if cell in pond_cell_set and period in {1, 2, 3}
                    else 1.5e-4
                    + 1.25e-4 * np.sin(period * 2.0 * np.pi / max(config.nper, 1))
                )
            )
            for _layer, cell in uzf_cells
        ]
        for period in range(config.nper)
    }
    pet = {
        period: [float(1.0e-4 + 5.0e-5 * np.cos(period * 2.0 * np.pi / max(config.nper, 1)))]
        * len(uzf_cells)
        for period in range(config.nper)
    }
    uzf_builder = UZFBuilder(
        context=ModelContext(grid=vor, domain=idomain),
        nper=config.nper,
        cells=uzf_cells,
        vks=0.25,
        thtr=0.08,
        thts=0.34,
        thti=0.17,
        eps=4.0,
        finf=finf,
        pet=pet,
        extdp=7.0,
    )
    uzf_builder.build().build(model.gwf)
    model.add_region_from_cells(
        "uzf_active",
        uzf_builder.uzf_cells,
        category="boundary",
        package="uzf",
        tags=["uzf"],
        metadata={"nuzfcells": len(uzf_builder.uzf_cells)},
        overwrite=True,
    )
    model.add_region_from_cells(
        "infiltration_pond",
        [(0, cell) for cell in infiltration_pond_cells],
        category="infiltration",
        package="uzf",
        tags=["mounding", "visualization"],
        overwrite=True,
    )

    _attach_canonical_targets(
        model,
        config,
        sfr,
        unconfined_seepage_cells=unconfined_seepage_cells,
        confined_seepage_cells=confined_seepage_cells,
        infiltration_pond_cells=infiltration_pond_cells,
        unconfined_well=unconfined_well,
        confined_well=confined_well,
    )
    OutputControl(
        model=model,
        save_record=(("HEAD", "ALL"), ("BUDGET", "ALL")),
    )
    model.canonical_contract = CANONICAL_MODEL_CONTRACT
    CANONICAL_MODEL_CONTRACT.validate(model, full_profile=config.ncpl >= 10_000)
    return model


def _write_gpkg(path: Path, gdf: gpd.GeoDataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _write_surface_water_inputs(workspace: Path, config: CanonicalModelConfig) -> dict[str, object]:
    width = config.ncol * config.cell_size
    height = config.nrow * config.cell_size
    crs = "EPSG:2927"

    def rectangle(x0, x1, y0, y1):
        return Polygon(
            [
                (width * x0, height * y0),
                (width * x0, height * y1),
                (width * x1, height * y1),
                (width * x1, height * y0),
            ]
        )

    north_lake = _write_gpkg(
        workspace / "inputs" / "north_lake.gpkg",
        gpd.GeoDataFrame({"name": ["north_lake"]}, geometry=[rectangle(0.72, 0.88, 0.62, 0.80)], crs=crs),
    )
    south_lake = _write_gpkg(
        workspace / "inputs" / "south_lake.gpkg",
        gpd.GeoDataFrame({"name": ["south_lake"]}, geometry=[rectangle(0.68, 0.86, 0.20, 0.38)], crs=crs),
    )
    stream = _write_gpkg(
        workspace / "inputs" / "stream.gpkg",
        gpd.GeoDataFrame(
            {"name": ["main_stream"]},
            geometry=[
                LineString(
                    [
                        (width * 0.03, height * 0.78),
                        (width * 0.13, height * 0.69),
                        (width * 0.23, height * 0.77),
                        (width * 0.33, height * 0.65),
                        (width * 0.43, height * 0.73),
                        (width * 0.53, height * 0.59),
                        (width * 0.63, height * 0.68),
                        (width * 0.73, height * 0.53),
                        (width * 0.83, height * 0.61),
                        (width * 0.93, height * 0.48),
                        (width * 0.98, height * 0.52),
                    ]
                )
            ],
            crs=crs,
        ),
    )
    return {"lakes": [north_lake, south_lake], "stream": stream}


def _attach_canonical_targets(
    model: SimulationBase,
    config: CanonicalModelConfig,
    sfr: SFRBuilder,
    *,
    unconfined_seepage_cells: list[int],
    confined_seepage_cells: list[int],
    infiltration_pond_cells: list[int],
    unconfined_well: int,
    confined_well: int,
) -> None:
    import flopy

    selected = representative_cells(config)
    pond_cell = infiltration_pond_cells[len(infiltration_pond_cells) // 2]
    head_locations = {
        "regional_center": (0, int(selected["center"])),
        "pond_mound": (0, pond_cell),
        "unconfined_pumping": (1, unconfined_well),
        "confined_pumping": (3, confined_well),
    }
    model.targets.heads = HeadTargets(
        locations=[
            {"name": name, "layer": layer, "cell": cell}
            for name, (layer, cell) in head_locations.items()
        ],
        values=pd.DataFrame(
            {
                "time": [period for period in range(config.nper) for _ in head_locations],
                "name": [
                    name
                    for _period in range(config.nper)
                    for name in head_locations
                ],
                "head": np.nan,
            }
        ),
    )
    model.targets.lake_stage = LakeStageTargets(locations={"north_lake": 0, "south_lake": 1})
    reaches = [int(value) for value in sfr.stream_reaches[sfr.stream_ids[0]]]
    model.targets.sfr_stage = SfrStageTargets(
        locations={"sfr_upstream": reaches[0], "sfr_midpoint": reaches[len(reaches) // 2]}
    )
    model.targets.sfr_flow = SfrFlowTargets(
        locations={"sfr_midflow": reaches[len(reaches) // 2], "sfr_outflow": reaches[-1]}
    )
    model.targets.drn_flow = DrnFlowTargets(
        locations=[
            {
                "name": "unconfined_seepage",
                "group": "seepage",
                "layer": 0,
                "cells": unconfined_seepage_cells,
            },
            {
                "name": "confined_seepage",
                "group": "seepage",
                "layer": 3,
                "cells": confined_seepage_cells,
            },
        ]
    )
    model.targets.heads.attach_flopy_obs(filename="master_heads.obs", csv_name="master_heads.csv")
    model.targets.lake_stage.attach_flopy_obs(filename="master_lakes.obs", csv_name="master_lakes.csv")
    sfr_continuous = {}
    sfr_continuous.update(model.targets.sfr_stage.to_flopy_obs(csv_name="master_sfr_stage.csv"))
    sfr_continuous.update(model.targets.sfr_flow.to_flopy_obs(csv_name="master_sfr_flow.csv"))
    flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
        model.gwf.sfr,
        pname="sfr_obs",
        continuous=sfr_continuous,
        filename="master_sfr.obs",
    )
    model.targets.drn_flow.attach_flopy_obs(filename="master_drn.obs", csv_name="master_drn.csv")


def representative_cells(config: CanonicalModelConfig) -> dict[str, int | list[int]]:
    """Return stable cells used by map, cross-section, and PRT examples."""

    center_row = config.nrow // 2
    center = center_row * config.ncol + config.ncol // 2
    cross_section = [center_row * config.ncol + col for col in range(config.ncol)]
    releases = [
        row * config.ncol + 2
        for row in (config.nrow // 4, config.nrow // 2, 3 * config.nrow // 4)
    ]
    return {"center": center, "cross_section": cross_section, "releases": releases}


# Preferred names. The older master-example names remain aliases for notebooks
# and scripts created before the canonical model became a package-level API.
MasterExampleConfig = CanonicalModelConfig
build_canonical_model = build_transient_model
