"""Build the authoritative package-rich canonical myflopy model.

The canonical model is a small but realistic **intermontane alluvial valley**.
Its conceptual story drives where every boundary condition goes, so each one
produces a *visibly distinct* signal that the diagnostic notebooks can read:

* The valley runs **west -> east**. Ground and the regional water table both
  decline down-valley and keep declining *past* the terminal lake toward a low
  general-head boundary at the valley mouth.
* Two tributary streams enter at the up-valley flanks, **converge** into one
  main stem, and the main stem **discharges into a terminal lake** (SFR -> LAK
  via MVR).
* Because the regional water table is high up-valley (mountain-front recharge)
  and low down-valley (the GHB underflow boundary plus deep pumping), the
  streams **gain** from the aquifer in the headwaters and **lose** to it near
  the lake -- the lake is perched over a declining regional water table.

Layering (icelltype ``(1, 1, 0, 0)``):

==== ======================================== ===========
lay  role                                      converts?
==== ======================================== ===========
L1   upper unconfined alluvium                 yes
L2   lower unconfined alluvium (main aquifer)  yes
L3   lacustrine-clay aquitard                  no
L4   confined basin-fill aquifer               no
==== ======================================== ===========
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Polygon

from myflopy.modflow.mf6.canonical import (
    CANONICAL_MODEL_CONTRACT,
    irregular_voronoi_grid,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.lakes import LAKBuilder
from myflopy.modflow.mf6.mvr import Move, MVRBuilder
from myflopy.modflow.mf6.observations import (
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)
from myflopy.modflow.mf6.sfr import SFRBuilder, StreamConnection
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
from myflopy.modflow.mf6.uzf import UZFBuilder
from myflopy.specs import ModelContext


@dataclass(frozen=True)
class CanonicalModelConfig:
    """Grid size and time-stepping knobs for the canonical valley example model.

    Selects the resolution and temporal extent of the canonical alluvial-valley
    model that :func:`build_canonical_model` assembles -- horizontal extent
    (``nrow``/``ncol`` of the irregular Voronoi grid, ``cell_size``), ``nlay``, and
    the stress-period schedule (``nper``/``period_length``/``steps_per_period``).
    The defaults are the full demonstration profile; use the smaller
    :meth:`validation` profile (50x50) for fast automated tests.

    Attributes
    ----------
    nrow, ncol, cell_size
        Grid resolution and cell size used to generate the Voronoi mesh.
    nlay
        Number of layers (the canonical contract expects 4).
    nper, period_length, steps_per_period
        Stress-period count, length, and sub-stepping.
    """

    nrow: int = 100
    ncol: int = 100
    nlay: int = 4
    nper: int = 6
    cell_size: float = 100.0
    period_length: float = 30.0
    steps_per_period: int = 2

    @classmethod
    def validation(cls) -> CanonicalModelConfig:
        """Return a smaller profile for automated API validation."""

        return cls(nrow=50, ncol=50, nlay=4, nper=6, cell_size=100.0)

    @classmethod
    def testing(cls) -> CanonicalModelConfig:
        """Return the smallest contract-complete profile, for the test suite.

        Feature-complete but small: every canonical package (13) and
        observation family (5) is still built and the full
        ``CANONICAL_MODEL_CONTRACT`` still validates — the grid is simply as
        coarse as the contract floors allow (SFR reaches >= 40, pond/spring
        cells present, three non-empty regions), and single time steps keep
        the transient schedule at the ``nper >= 6`` floor. The notebooks and
        the weekly slow CI lane keep using :meth:`validation` /
        the full default profile; see ``docs/phase_baselines.md``.
        """

        return cls(nrow=21, ncol=21, nlay=4, nper=6, cell_size=100.0, steps_per_period=1)

    @property
    def ncpl(self) -> int:
        """Cells per layer (``nrow * ncol``)."""

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
    """Assemble the canonical transient alluvial-valley GWF model in ``workspace``.

    Builds the project's reference model end to end -- the irregular Voronoi grid,
    layered valley surfaces, and the full package set (NPF/IC/STO/OC plus the
    surface-water network: GHB/DRN/RCH/UZF/SFR/LAK/MVR) sized and timed by
    ``config`` -- and returns the assembled :class:`SimulationBase` (not yet run).
    This single model is reused across the examples, the integration tests
    (validated against :class:`CanonicalModelContract`), and the PEST notebooks.
    Exported under the alias ``build_canonical_model``.

    Built in the legacy imperative builder style for its computed-cell geometry;
    that is historical, not a recommendation -- new model assembly should use the
    package-first API.

    Parameters
    ----------
    workspace
        Directory to build the model in.
    config
        Size/timing controls; defaults to the full :class:`CanonicalModelConfig`.
    name
        Model name (must satisfy MF6 length limits).

    Returns
    -------
    SimulationBase
        The assembled (unrun) canonical model.
    """

    config = CanonicalModelConfig() if config is None else config
    workspace = Path(workspace)
    vor = rectangular_voronoi(config)
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=config.nper)

    centers = np.asarray(vor.points, dtype=float)
    width = config.ncol * config.cell_size
    height = config.nrow * config.cell_size
    # xn: 0 at the up-valley (west) head, 1 at the (east) valley mouth.
    # yn: 0..1 across the valley; yn == 0.5 is the valley axis.
    xn = centers[:, 0] / width
    yn = centers[:, 1] / height
    across = 2.0 * np.abs(yn - 0.5)  # 0 on the axis, 1 at either wall

    def nearest_cell(x_fraction: float, y_fraction: float) -> int:
        """The cell nearest a fractional ``(x, y)`` position (0-1) within the domain."""

        return int(np.argmin((xn - x_fraction) ** 2 + (yn - y_fraction) ** 2))

    # --- Valley topography ---------------------------------------------------
    # A strong down-valley slope (the dominant gradient) plus a gentle
    # cross-valley "U" so the upper layers stay wet across the whole domain.
    axis_top = 150.0 - 52.0 * xn
    cross_relief = 14.0 * across ** 1.3 * (1.0 - 0.25 * xn)
    top = axis_top + cross_relief
    bottom = np.vstack(
        [
            top - 28.0,  # base of L1, upper unconfined alluvium
            top - 52.0,  # base of L2, lower unconfined (main) aquifer
            top - 64.0,  # base of L3, 12 ft lacustrine-clay aquitard
            top - 100.0,  # base of L4, confined basin-fill aquifer
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

    # --- Regional water table ------------------------------------------------
    # Declines from a high up-valley head to a low head at the mouth and keeps
    # declining past the lake. This is what makes the streams gain in the
    # headwaters and lose near the (perched) lake.
    regional = 146.0 - 60.0 * xn + 0.45 * cross_relief
    initial = np.vstack(
        [np.maximum(bottom[layer] + 1.0, regional - 0.6 * layer) for layer in range(config.nlay)]
    )
    InitialConditions(model=model, vor=vor, nlay=config.nlay, strt=initial)

    # --- Hydraulic conductivity ----------------------------------------------
    # A high-K paleochannel ribbon down the valley axis channels flow toward the
    # streams and lake; L3 is a distinct low-K lacustrine aquitard. A cross-
    # valley bedrock constriction near mid-valley forces the regional water
    # table to step down, so the lower valley + lake drain below the streambed
    # (the streams gain above the narrows and lose below it).
    paleochannel = 30.0 * np.exp(-((yn - 0.5) / 0.13) ** 2)
    # Coarser, more transmissive basin fill toward the mouth lets the low GHB
    # head propagate up-valley; a low-K constriction makes the step.
    constriction = np.where((xn > 0.52) & (xn < 0.60), 0.06, 1.0)
    aquifer = (12.0 + 18.0 * xn + paleochannel) * constriction
    conductivity = np.vstack(
        [
            aquifer,
            aquifer * 0.85,
            np.full_like(aquifer, 0.02),  # lacustrine clay aquitard
            np.maximum(aquifer * 0.55, 0.5),
        ]
    )
    KFlow(
        model=model,
        k=conductivity,
        k33_vert=np.vstack(
            [
                aquifer * 0.1,
                aquifer * 0.1,
                np.full_like(aquifer, 0.0015),
                np.maximum(aquifer * 0.05, 0.05),
            ]
        ),
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

    # --- CHD: up-valley regional inflow (west edge) --------------------------
    west_cells = [row * config.ncol for row in range(config.nrow)]
    chd_data = {
        period: [
            ((layer, cell), float(max(bottom[layer, cell] + 0.5, top[cell] - 6.0)))
            for layer in range(config.nlay)
            for cell in west_cells
        ]
        for period in range(config.nper)
    }
    CHD(model=model, stress_period_data=chd_data)

    # --- GHB: underflow leaving at the valley mouth (east edge) --------------
    # Heads are well below the lake stage so the regional gradient continues
    # down-valley past the lake.
    east_cells = [row * config.ncol + config.ncol - 1 for row in range(config.nrow)]
    ghb_data = {
        period: [
            (
                (layer, cell),
                float(max(bottom[layer, cell] + 1.0, 80.0 + 0.3 * cross_relief[cell] - 0.4 * layer)),
                4000.0,
            )
            for layer in range(config.nlay)
            for cell in east_cells
        ]
        for period in range(config.nper)
    }
    GHB(model=model, stress_period_data=ghb_data)

    # --- RCH: mountain-front recharge on the two valley walls ----------------
    wall = np.flatnonzero(across >= 0.62)
    mountain_front_recharge = {
        period: [
            (
                (0, int(cell)),
                float((3.2e-4 - 1.6e-4 * xn[cell]) * (1.0 + 0.35 * np.sin(period * 2.0 * np.pi / config.nper))),
            )
            for cell in wall
        ]
        for period in range(config.nper)
    }
    Recharge(model=model, rch_dict=mountain_front_recharge)

    # --- WEL: shallow (L2) and deep (L4) valley-floor pumping ----------------
    # Period 0 is a no-pump steady-state baseline so drawdown/capture is the
    # difference between any pumping period and period 0.
    shallow_well = nearest_cell(0.58, 0.44)
    deep_well = nearest_cell(0.52, 0.55)
    well_data = {
        period: [
            (
                (1, shallow_well),
                0.0 if period == 0 else -2400.0 - (3800.0 if period in {2, 3} else 0.0),
                "shallow_supply",
            ),
            (
                (3, deep_well),
                0.0 if period == 0 else -3000.0 - (5200.0 if period in {3, 4} else 0.0),
                "deep_supply",
            ),
        ]
        for period in range(config.nper)
    }
    Wells(model=model, stress_period_data=well_data)

    # --- DRN: toe-of-slope springs in the gaining upper valley ---------------
    north_spring_cells = np.flatnonzero(
        (xn > 0.12) & (xn < 0.48) & (yn > 0.76) & (yn < 0.86)
    )[::2].tolist()
    south_spring_cells = np.flatnonzero(
        (xn > 0.12) & (xn < 0.48) & (yn > 0.14) & (yn < 0.24)
    )[::2].tolist()
    drain_data = {
        period: (
            [
                ((0, cell), float(regional[cell] - 1.0), 60.0 + 4.0 * (cell % 7))
                for cell in north_spring_cells
            ]
            + [
                ((0, cell), float(regional[cell] - 1.0), 60.0 + 4.0 * (cell % 5))
                for cell in south_spring_cells
            ]
        )
        for period in range(config.nper)
    }
    Drains(model=model, stress_period_data=drain_data)
    model.add_region_from_cells(
        "north_seepage_springs",
        [(0, cell) for cell in north_spring_cells],
        category="seepage",
        package="drn",
        overwrite=True,
    )
    model.add_region_from_cells(
        "south_seepage_springs",
        [(0, cell) for cell in south_spring_cells],
        category="seepage",
        package="drn",
        overwrite=True,
    )

    # --- LAK / SFR / MVR surface-water network -------------------------------
    inputs = _write_surface_water_inputs(workspace, config)
    lake_id = "valley_lake"
    lak = LAKBuilder(
        context=ModelContext(grid=vor, domain=idomain),
        nper=config.nper,
        lakes=inputs["lakes"],
        lake_id_field="name",
        starting_stage={lake_id: 101.0},
        lake_bottom={lake_id: 96.0},
        lake_top={lake_id: 101.0},  # flat-bottom lake -> rectangular; rim at the stage
        bed_leakance=0.11,
        connection_modes="rectangular",
        status={lake_id: ["ACTIVE"] * config.nper},
        mover=True,
        length_conversion=3.28081,
        time_conversion=86_400.0,
    )
    lak.build().build(model.gwf)
    for built_lake_id, cells in lak.lake_cells.items():
        model.add_region_from_cells(
            f"lake_zone_{built_lake_id}",
            [(0, cell) for cell in cells],
            category="boundary",
            package="lak",
            tags=["lak"],
            geometry=lak.lake_table.loc[built_lake_id, "geometry"],
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
        period: {
            "north_trib": 18_000.0 + 4_000.0 * np.sin(period * 2.0 * np.pi / max(config.nper, 1)),
            "south_trib": 12_000.0 + 3_000.0 * np.sin(period * 2.0 * np.pi / max(config.nper, 1)),
        }
        for period in range(config.nper)
    }
    sfr = SFRBuilder(
        context=ModelContext(grid=vor, domain=idomain),
        nper=config.nper,
        streams=inputs["streams"],
        connection_mode="automatic",
        connections=(
            StreamConnection("north_trib", "main_stem"),
            StreamConnection("south_trib", "main_stem"),
        ),
        inflow=inflows,
        width=18.0,
        gradient=0.0012,
        roughness=0.030,
        streambed_k=0.05,
        streambed_thickness=1.5,
        length_conversion=3.28081,
        time_conversion=86_400.0,
        mover=True,
    )
    # Streambed tops: a few feet below the regional water table in the
    # headwaters (gaining) rising above it toward the lake (losing).
    reach_xn = np.asarray(
        [xn[int(cellid[1])] for cellid in sfr.reaches["cellid"]], dtype=float
    )
    streambed_offset = -4.0 + 9.0 * np.clip((reach_xn - 0.30) / 0.55, 0.0, 1.0)
    reach_regional = 146.0 - 60.0 * reach_xn
    reach_tops = (reach_regional + streambed_offset).tolist()
    sfr = sfr.with_updates(
        reach_top=reach_tops,
        width=[14.0 + 12.0 * value for value in reach_xn],
        streambed_k=[0.02 + 0.10 * value for value in reach_xn],
    )
    sfr.build().build(model.gwf)
    main_stem = "main_stem"
    model.add_region_from_cells(
        "all_streams",
        [
            (0, cell)
            for cells in sfr.stream_cells.values()
            for cell in cells
        ],
        category="boundary",
        package="sfr",
        tags=["sfr"],
        overwrite=True,
    )
    mvr = MVRBuilder(
        nper=config.nper,
        moves={
            period: [
                Move(sfr.connection(main_stem, "downstream"), lak.connection(lake_id), value=1.0)
            ]
            for period in range(config.nper)
        },
    )
    mvr.build().build(model.gwf)

    # --- UZF: unsaturated zone + ET across the valley floor ------------------
    lake_cells = {int(cell) for cells in lak.lake_cells.values() for cell in cells}
    stream_cells = {int(cell) for cells in sfr.stream_cells.values() for cell in cells}
    floor = (across < 0.62)
    uzf_cells = [
        (0, cell)
        for cell in range(config.ncpl)
        if floor[cell] and cell not in lake_cells and cell not in stream_cells
    ]
    infiltration_pond_cells = [
        cell
        for _layer, cell in uzf_cells
        if ((xn[cell] - 0.27) / 0.07) ** 2 + ((yn[cell] - 0.5) / 0.06) ** 2 <= 1.0
    ]
    pond_cell_set = set(infiltration_pond_cells)
    finf = {
        period: [
            float(
                5.0e-4
                if cell in pond_cell_set and period in {1, 2, 3}
                else 3.0e-5 + 2.0e-5 * np.sin(period * 2.0 * np.pi / max(config.nper, 1))
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
        main_stem=main_stem,
        lake_id=lake_id,
        north_spring_cells=north_spring_cells,
        south_spring_cells=south_spring_cells,
        infiltration_pond_cells=infiltration_pond_cells,
        shallow_well=shallow_well,
        deep_well=deep_well,
    )
    OutputControl(
        model=model,
        save_record=(("HEAD", "ALL"), ("BUDGET", "ALL")),
    )
    model.canonical_contract = CANONICAL_MODEL_CONTRACT
    CANONICAL_MODEL_CONTRACT.validate(model, full_profile=config.ncpl >= 10_000)
    return model


def _write_gpkg(path: Path, gdf: gpd.GeoDataFrame) -> Path:
    """Write a GeoDataFrame to a GeoPackage at ``path`` (overwriting any existing file)."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _write_surface_water_inputs(workspace: Path, config: CanonicalModelConfig) -> dict[str, object]:
    """Write the terminal lake and the two converging tributaries to GeoPackages."""

    width = config.ncol * config.cell_size
    height = config.nrow * config.cell_size
    crs = "EPSG:2927"

    def point(x_fraction: float, y_fraction: float) -> tuple[float, float]:
        """Model coordinates for a fractional ``(x, y)`` position (0-1) within the domain."""

        return (width * x_fraction, height * y_fraction)

    def rectangle(x0, x1, y0, y1):
        """A rectangle polygon spanning the fractional bounds ``[x0, x1] x [y0, y1]``."""

        return Polygon([point(x0, y0), point(x0, y1), point(x1, y1), point(x1, y0)])

    lake = _write_gpkg(
        workspace / "inputs" / "valley_lake.gpkg",
        gpd.GeoDataFrame(
            {"name": ["valley_lake"]},
            geometry=[rectangle(0.80, 0.92, 0.40, 0.60)],
            crs=crs,
        ),
    )
    # Two tributaries meet at the confluence (0.40, 0.50); the main stem carries
    # the combined flow east to the lake's west shore (0.80, 0.50).
    confluence = point(0.40, 0.50)
    north_trib = _write_gpkg(
        workspace / "inputs" / "north_trib.gpkg",
        gpd.GeoDataFrame(
            {"name": ["north_trib"]},
            geometry=[
                LineString([point(0.03, 0.85), point(0.14, 0.78), point(0.26, 0.66), confluence])
            ],
            crs=crs,
        ),
    )
    south_trib = _write_gpkg(
        workspace / "inputs" / "south_trib.gpkg",
        gpd.GeoDataFrame(
            {"name": ["south_trib"]},
            geometry=[
                LineString([point(0.03, 0.15), point(0.14, 0.22), point(0.26, 0.34), confluence])
            ],
            crs=crs,
        ),
    )
    main_stem = _write_gpkg(
        workspace / "inputs" / "main_stem.gpkg",
        gpd.GeoDataFrame(
            {"name": ["main_stem"]},
            geometry=[
                LineString([confluence, point(0.55, 0.50), point(0.68, 0.50), point(0.80, 0.50)])
            ],
            crs=crs,
        ),
    )
    return {"lakes": [lake], "streams": [north_trib, south_trib, main_stem]}


def _attach_canonical_targets(
    model: SimulationBase,
    config: CanonicalModelConfig,
    sfr: SFRBuilder,
    *,
    main_stem: str,
    lake_id: str,
    north_spring_cells: list[int],
    south_spring_cells: list[int],
    infiltration_pond_cells: list[int],
    shallow_well: int,
    deep_well: int,
) -> None:
    """Register the canonical model's head/stage/flow observation targets on ``model``."""

    import flopy

    selected = representative_cells(config)
    pond_cell = infiltration_pond_cells[len(infiltration_pond_cells) // 2]
    head_locations = {
        "regional_center": (0, int(selected["center"])),
        "pond_mound": (0, pond_cell),
        "shallow_pumping": (1, shallow_well),
        "deep_pumping": (3, deep_well),
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
    model.targets.lake_stage = LakeStageTargets(locations={lake_id: 0})
    reaches = [int(value) for value in sfr.stream_reaches[main_stem]]
    model.targets.sfr_stage = SfrStageTargets(
        locations={"sfr_upstream": reaches[0], "sfr_midpoint": reaches[len(reaches) // 2]}
    )
    model.targets.sfr_flow = SfrFlowTargets(
        locations={"sfr_midflow": reaches[len(reaches) // 2], "sfr_outflow": reaches[-1]}
    )
    model.targets.drn_flow = DrnFlowTargets(
        locations=[
            {
                "name": "north_springs",
                "group": "seepage",
                "layer": 0,
                "cells": north_spring_cells,
            },
            {
                "name": "south_springs",
                "group": "seepage",
                "layer": 0,
                "cells": south_spring_cells,
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
    # PRT release points along the up-valley mountain front and floor.
    releases = [
        row * config.ncol + 2
        for row in (config.nrow // 5, config.nrow // 2, 4 * config.nrow // 5)
    ]
    return {"center": center, "cross_section": cross_section, "releases": releases}


# Preferred names. The older master-example names remain aliases for notebooks
# and scripts created before the canonical model became a package-level API.
MasterExampleConfig = CanonicalModelConfig
build_canonical_model = build_transient_model
