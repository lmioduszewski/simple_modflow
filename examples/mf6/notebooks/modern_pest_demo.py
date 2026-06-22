"""Compact calibration demo model for the modern PEST notebooks.

Builds a small single-layer Voronoi groundwater model with a known "truth"
hydraulic conductivity, samples synthetic head observations (and one forecast)
from it, then resets the model to a deliberately wrong starting K. The returned
model is the *starting* model to calibrate; the targets carry the truth-derived
"measured" values.

Kept deliberately small (a few hundred cells, one layer, steady state) so the
whole calibrate-and-assess loop -- including a real PESTPP-IES run -- finishes
in a notebook-friendly amount of time.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations import HeadTargets
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization
from myflopy.modflow.mf6.simulation.packages import (
    CHD,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)


@dataclass
class CalibrationDemo:
    """Everything the PEST notebooks need to calibrate the compact demo."""

    model: SimulationBase
    head_targets: HeadTargets
    forecast_targets: HeadTargets
    truth_k: float
    start_k: float
    workspace: Path


def _cell_centroids(vor) -> np.ndarray:
    """Return an (ncpl, 2) array of cell-centroid coordinates."""

    centroids = vor.gdf_vorPolys.geometry.centroid
    return np.column_stack([centroids.x.to_numpy(), centroids.y.to_numpy()])


def build_calibration_demo(
    workspace: str | Path,
    *,
    nrow: int = 12,
    ncol: int = 12,
    cell_size: float = 100.0,
    truth_k: float = 2.0,
    start_k: float = 10.0,
    recharge: float = 5.0e-4,
    n_wells: int = 8,
    seed: int = 12,
) -> CalibrationDemo:
    """Build the compact demo: truth-derived observations + a wrong starting model.

    Parameters
    ----------
    workspace
        Directory for the model files.
    nrow, ncol, cell_size
        Controls the irregular Voronoi grid extent/resolution.
    truth_k
        Uniform "true" hydraulic conductivity used to generate observations.
    start_k
        Uniform (wrong) starting conductivity the calibration must correct.
    recharge
        Uniform recharge rate (L/T).
    n_wells
        Number of synthetic head-observation wells.
    seed
        RNG seed for reproducible well placement.

    Returns
    -------
    CalibrationDemo
        ``model`` is the starting model (wrong K); ``head_targets`` carry the
        truth-derived measured heads; ``forecast_targets`` is one head
        prediction with no calibration weight.
    """

    workspace = Path(workspace)
    vor = irregular_voronoi_grid(nrow=nrow, ncol=ncol, cell_size=cell_size)
    ncpl = int(vor.ncpl)
    width = ncol * cell_size

    top = np.full(ncpl, 50.0)
    bottom = np.full((1, ncpl), 0.0)
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {"geometry": vor.gdf_vorPolys.geometry, 0: top, 1: bottom[0]},
        geometry="geometry",
        crs=vor.crs,
    )

    model = SimulationBase(name="pest_demo", mf_folder_path=workspace, vor=vor, nper=1)
    DisvGrid(vor=vor, model=model, top=top, bottom=bottom, nlay=1, idomain=np.ones((1, ncpl), dtype=int))
    TemporalDiscretization(model=model, per_len=1.0, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=np.full(ncpl, 30.0))
    KFlow(model=model, k=np.full(ncpl, truth_k), k33_vert=np.full(ncpl, truth_k * 0.1),
          save_specific_discharge=False)
    Storage(model=model, sto_steady={0: True}, sto_transient={})
    OutputControl(model=model)

    centroids = _cell_centroids(vor)
    x = centroids[:, 0]
    # Constant-head boundaries on the left (high) and right (low) edges drive flow.
    left_cells = np.where(x <= cell_size)[0]
    right_cells = np.where(x >= width - cell_size)[0]
    chd_rows = [[(0, int(cell)), 40.0] for cell in left_cells]
    chd_rows += [[(0, int(cell)), 20.0] for cell in right_cells]
    CHD(model=model, stress_period_data={0: chd_rows})
    Recharge(model=model, vor=vor, rch_dict={0: [[(0, int(c)), recharge] for c in range(ncpl)]})

    # Run the truth model and sample synthetic observations from it.
    model.run_simulation()

    rng = np.random.default_rng(seed)
    boundary = set(left_cells.tolist()) | set(right_cells.tolist())
    interior = [c for c in range(ncpl) if c not in boundary]
    well_cells = rng.choice(interior, size=min(n_wells, len(interior)), replace=False)
    forecast_cell = int(rng.choice([c for c in interior if c not in set(well_cells.tolist())]))

    def _targets_for(cells, prefix):
        names = [f"{prefix}_{i:02d}" for i in range(len(cells))]
        points = [Point(*centroids[int(c)]) for c in cells]
        wells = gpd.GeoDataFrame(
            {"name": names, "layer": [0] * len(cells), "group": ["heads"] * len(cells),
             "weight": [1.0] * len(cells)},
            geometry=points, crs=vor.crs,
        )
        placeholder = pd.DataFrame({"per": [0], **{name: [0.0] for name in names}})
        sampled = HeadTargets(locations=wells, values=placeholder, time_column="per").compare(model)
        wide = sampled.pivot(index="time", columns="name", values="sim_head").reset_index()
        wide = wide.rename(columns={"time": "per"})
        return HeadTargets(locations=wells, values=wide, time_column="per")

    head_targets = _targets_for(list(well_cells), "obs")
    forecast_targets = _targets_for([forecast_cell], "pred")

    # Reset K to the wrong, uniform starting value: this is the model to calibrate.
    model.gwf.npf.k.set_data(np.full(ncpl, start_k))
    model.run_simulation()

    return CalibrationDemo(
        model=model,
        head_targets=head_targets,
        forecast_targets=forecast_targets,
        truth_k=truth_k,
        start_k=start_k,
        workspace=workspace,
    )
