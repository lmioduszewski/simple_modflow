"""Dataclass specifications for the first ``myflopy`` PEST slice."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from myflopy.modflow.mf6.observations import (
    DrnFlowTargets,
    HeadTargets,
    LakeStageTargets,
    SfrFlowTargets,
    SfrStageTargets,
)


@dataclass
class VectorParameterSource:
    """Describe a GIS vector source used to define calibration parameters.

    Parameters
    ----------
    path
        Path to a shapefile/geopackage layer.
    value_column
        Base-value column such as ``"k"``, ``"elev"``, or ``"cond"``.
    feature_id_column
        Unique feature identifier column used for feature-based parameters.
    zone_column
        Optional spatial zone column.
    group_column
        Optional logical group column for grouping parameters or observations.
    layer_column
        Optional layer column.
    lower_bound_column, upper_bound_column
        Optional GIS-defined bound columns.
    """

    path: str | Path
    value_column: str
    feature_id_column: str | None = "name"
    zone_column: str | None = None
    group_column: str | None = None
    layer_column: str | None = None
    crs: int | None = 2927
    lower_bound_column: str | None = None
    upper_bound_column: str | None = None


@dataclass
class ExpGeoStruct:
    """Simple exponential geostatistical structure definition.

    This wrapper stores values in a ``myflopy``-friendly form and can be
    converted to a real ``pyemu.geostats.GeoStruct`` lazily when pyEMU is
    available.
    """

    range: float
    contribution: float = 1.0
    anisotropy: float = 1.0
    bearing: float = 0.0
    nugget: float = 0.0
    transform: str = "none"


@dataclass
class HeadTargetObservationSpec:
    """Use :class:`HeadTargets` as a pyEMU-compatible observation source."""

    targets: HeadTargets
    simulated_values: str | Path | pd.DataFrame | None = None
    prefix: str = "hds"


@dataclass
class LakeStageObservationSpec:
    """Lake-stage observation specification.

    Parameters
    ----------
    lake_names
        Names of lake-stage series to expose as observations.
    values
        Optional target-value table in the same wide/long conventions used by
        :class:`HeadTargets`.
    time_column
        Time/per identifier column when ``values`` is long-form.
    value_column
        Stage-value column when ``values`` is long-form.
    weight
        Optional uniform observation weight to apply when target values are
        supplied.
    prefix
        pyEMU observation-name prefix.
    """

    lake_names: list[str] | None = None
    values: str | Path | pd.DataFrame | None = None
    time_column: str = "time"
    value_column: str = "stage"
    weight: float | None = None
    prefix: str = "stage"
    targets: LakeStageTargets | None = None


@dataclass
class SfrStageObservationSpec:
    """Use :class:`SfrStageTargets` as a pyEMU-compatible observation source."""

    targets: SfrStageTargets
    prefix: str = "sfr_stage"


@dataclass
class SfrFlowObservationSpec:
    """Use :class:`SfrFlowTargets` as a pyEMU-compatible observation source."""

    targets: SfrFlowTargets
    prefix: str = "sfr_flow"


@dataclass
class DrnFlowObservationSpec:
    """Use :class:`DrnFlowTargets` as a pyEMU-compatible observation source."""

    targets: DrnFlowTargets
    prefix: str = "drn_flow"
