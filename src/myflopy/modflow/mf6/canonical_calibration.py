"""Synthetic calibration setup built on the canonical valley model.

The PEST notebooks calibrate the **canonical valley model itself** rather than a
throwaway demo. :func:`build_canonical_calibration_demo` builds the canonical
model as the synthetic *truth*, samples truth-derived head observations (plus a
downgradient head forecast), then resets hydraulic conductivity to a
deliberately wrong starting value. The returned model is the *starting* model to
calibrate; the targets carry the truth-derived "measured" values.

This replaces the retired standalone demo models (``gold_standard_demo``,
``synthetic_demo``, ``modern_pest_demo.build_calibration_demo``) so every
notebook uses one model.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    build_canonical_model,
)
from myflopy.modflow.mf6.observations import HeadTargets
from myflopy.modflow.mf6.simulation.base import SimulationBase


@dataclass
class CanonicalCalibrationDemo:
    """Everything the PEST notebooks need to calibrate the canonical model.

    Attributes
    ----------
    model
        The *starting* canonical model (K perturbed to a wrong value); this is
        what PEST calibrates.
    head_targets
        Truth-derived head observations spread across the valley floor -- the
        weighted history-matching target.
    forecast_targets
        A single downgradient head prediction (near the lake), registered with
        zero weight as the forecast of interest for uncertainty analysis.
    start_k_factor
        Factor the truth K field was multiplied by to make the wrong start.
    workspace
        Directory the starting model was written to.
    """

    model: SimulationBase
    head_targets: HeadTargets
    forecast_targets: HeadTargets
    start_k_factor: float
    workspace: Path


def _truth_head_targets(model, cells, prefix: str, *, layer: int = 0) -> HeadTargets:
    """Sample simulated heads at ``cells`` and return them as measured targets."""

    cells = [int(c) for c in cells]
    names = (
        [f"{prefix}_{i:02d}" for i in range(len(cells))] if len(cells) > 1 else [prefix]
    )
    locations = [
        {"name": name, "layer": layer, "cell": cell}
        for name, cell in zip(names, cells)
    ]
    placeholder = pd.DataFrame(
        {"per": list(range(model.nper)), **{name: [np.nan] * model.nper for name in names}}
    )
    sampler = HeadTargets(locations=locations, values=placeholder, time_column="per")
    measured = sampler.simulated_heads(model)
    return HeadTargets(locations=locations, values=measured, time_column="per")


def build_canonical_calibration_demo(
    workspace: str | Path,
    *,
    config: CanonicalModelConfig | None = None,
    n_head_wells: int = 16,
    start_k_factor: float = 3.0,
    seed: int = 2026,
) -> CanonicalCalibrationDemo:
    """Build the canonical model as truth, sample observations, perturb the start.

    Parameters
    ----------
    workspace
        Directory for the starting model files.
    config
        Canonical model profile. Defaults to
        :meth:`CanonicalModelConfig.validation` (50x50). Use the full profile
        for a faithful but slower calibration.
    n_head_wells
        Number of truth-sampled head observation wells placed across the valley
        floor (the ``uzf_active`` region, away from lake/stream/boundary cells).
    start_k_factor
        Factor applied to the truth K field to create the wrong starting model
        the calibration must correct (default ``3.0`` = 3x too transmissive).
    seed
        RNG seed for reproducible well placement.

    Returns
    -------
    CanonicalCalibrationDemo
        ``model`` is the starting model (wrong K); ``head_targets`` carry the
        truth-derived measured heads; ``forecast_targets`` is one downgradient
        head prediction.
    """

    config = config or CanonicalModelConfig.validation()
    workspace = Path(workspace)

    truth = build_canonical_model(workspace, config=config)
    success, report = truth.run_simulation()
    if not success:
        raise RuntimeError("Truth canonical model failed to run:\n" + "\n".join(report[-20:]))

    # Observation wells: valley-floor cells (the UZF region excludes lake, stream
    # and boundary cells), sampled reproducibly.
    floor = sorted(int(c) for c in truth.get_region_cells("uzf_active"))
    rng = np.random.default_rng(seed)
    well_cells = sorted(
        int(c) for c in rng.choice(floor, size=min(n_head_wells, len(floor)), replace=False)
    )

    # Forecast point: a head down-valley near the lake (not an observation well).
    centers = np.asarray(truth.vor.points, dtype=float)
    width = config.ncol * config.cell_size
    height = config.nrow * config.cell_size
    xn = centers[:, 0] / width
    yn = centers[:, 1] / height
    forecast_cell = int(np.argmin((xn - 0.72) ** 2 + (yn - 0.5) ** 2))

    head_targets = _truth_head_targets(truth, well_cells, "obs", layer=0)
    forecast_targets = _truth_head_targets(truth, [forecast_cell], "fore_lakehead", layer=0)

    # Reset K to the wrong, uniformly-too-high starting value -- the model to
    # calibrate. set_all_data_external() later splits this into per-layer files.
    truth_k = np.asarray(truth.gwf.npf.k.get_data(), dtype=float)
    truth.gwf.npf.k.set_data(truth_k * float(start_k_factor))
    truth.run_simulation()

    return CanonicalCalibrationDemo(
        model=truth,
        head_targets=head_targets,
        forecast_targets=forecast_targets,
        start_k_factor=float(start_k_factor),
        workspace=workspace,
    )
