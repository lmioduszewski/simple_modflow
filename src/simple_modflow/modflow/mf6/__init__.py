"""Public MF6-facing exports for simple_modflow."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.boundaries import Boundaries
    from simple_modflow.modflow.mf6.canonical import (
        CANONICAL_MODEL_CONTRACT,
        CanonicalModelContract,
        canonical_feature_signals,
        canonical_head_signals,
        canonical_partition_mask,
        canonical_sfr_signals,
        irregular_voronoi_grid,
    )
    from simple_modflow.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )
    from simple_modflow.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from simple_modflow.modflow.mf6.interactive_plotting import (
        FrameExportProgress,
        ModelMapStyle,
        ModelVisualization,
        ParticleTrackingScene,
        StandaloneHtmlSlider,
    )
    from simple_modflow.modflow.mf6.observations import (
        DrnFlowTargets,
        HeadTargets,
        LakeStageTargets,
        SfrFlowTargets,
        SfrStageTargets,
    )
    from simple_modflow.modflow.mf6.parallel import (
        ParallelCompatibilityError,
        ParallelEnvironment,
        ParallelModelWorkflow,
        ParallelSplitResults,
        ParallelSplitRun,
    )
    from simple_modflow.modflow.mf6.pest import PestProject, PestRunResults, PestRunReview, open_pest_run
    from simple_modflow.modflow.mf6.prt import (
        ParticleTracking,
        PRTProject,
        PRTReleasePoints,
        PRTRunResults,
        open_prt_run,
    )
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase
    from simple_modflow.modflow.mf6.simulation.discretization import (
        DisuGrid,
        DisvGrid,
        TemporalDiscretization,
    )
    from simple_modflow.modflow.mf6.simulation.packages import (
        CHD,
        LAK,
        MVR,
        UZF,
        Drains,
        InitialConditions,
        KFlow,
        OutputControl,
        Recharge,
        Storage,
        Wells,
    )


_SUBMODULES = {
    "grid": "simple_modflow.modflow.mf6.grid",
    "pest": "simple_modflow.modflow.mf6.pest",
    "simulation": "simple_modflow.modflow.mf6.simulation",
}

_EXPORTS = {
    "Boundaries": ("simple_modflow.modflow.mf6.boundaries", "Boundaries"),
    "CHD": ("simple_modflow.modflow.mf6.simulation.packages", "CHD"),
    "CANONICAL_MODEL_CONTRACT": ("simple_modflow.modflow.mf6.canonical", "CANONICAL_MODEL_CONTRACT"),
    "CanonicalModelContract": ("simple_modflow.modflow.mf6.canonical", "CanonicalModelContract"),
    "CanonicalModelConfig": ("simple_modflow.modflow.mf6.canonical_example", "CanonicalModelConfig"),
    "CHDFromVector": ("simple_modflow.modflow.mf6.chd", "CHDFromVector"),
    "ModelCrossSectionStyle": ("simple_modflow.modflow.mf6.cross_section_plotting", "ModelCrossSectionStyle"),
    "ModelMapStyle": ("simple_modflow.modflow.mf6.interactive_plotting", "ModelMapStyle"),
    "FrameExportProgress": ("simple_modflow.modflow.mf6.interactive_plotting", "FrameExportProgress"),
    "ModelVisualization": ("simple_modflow.modflow.mf6.interactive_plotting", "ModelVisualization"),
    "ParticleTrackingScene": ("simple_modflow.modflow.mf6.interactive_plotting", "ParticleTrackingScene"),
    "ParticleTracking": ("simple_modflow.modflow.mf6.prt", "ParticleTracking"),
    "ParallelCompatibilityError": ("simple_modflow.modflow.mf6.parallel", "ParallelCompatibilityError"),
    "ParallelEnvironment": ("simple_modflow.modflow.mf6.parallel", "ParallelEnvironment"),
    "ParallelModelWorkflow": ("simple_modflow.modflow.mf6.parallel", "ParallelModelWorkflow"),
    "ParallelSplitResults": ("simple_modflow.modflow.mf6.parallel", "ParallelSplitResults"),
    "ParallelSplitRun": ("simple_modflow.modflow.mf6.parallel", "ParallelSplitRun"),
    "StandaloneHtmlSlider": ("simple_modflow.modflow.mf6.interactive_plotting", "StandaloneHtmlSlider"),
    "DRN": ("simple_modflow.modflow.mf6.drn", "DRN"),
    "DRNFromVector": ("simple_modflow.modflow.mf6.drn", "DRNFromVector"),
    "DisuGrid": ("simple_modflow.modflow.mf6.simulation.discretization", "DisuGrid"),
    "DisvGrid": ("simple_modflow.modflow.mf6.simulation.discretization", "DisvGrid"),
    "Drains": ("simple_modflow.modflow.mf6.simulation.packages", "Drains"),
    "GHB": ("simple_modflow.modflow.mf6.ghb", "GHB"),
    "GHBFromVector": ("simple_modflow.modflow.mf6.ghb", "GHBFromVector"),
    "HeadTargets": ("simple_modflow.modflow.mf6.observations", "HeadTargets"),
    "DrnFlowTargets": ("simple_modflow.modflow.mf6.observations", "DrnFlowTargets"),
    "LakeStageTargets": ("simple_modflow.modflow.mf6.observations", "LakeStageTargets"),
    "InitialConditions": ("simple_modflow.modflow.mf6.simulation.packages", "InitialConditions"),
    "KFlow": ("simple_modflow.modflow.mf6.simulation.packages", "KFlow"),
    "KFromVector": ("simple_modflow.modflow.mf6.kflow", "KFromVector"),
    "KPilotPointParameter": ("simple_modflow.modflow.mf6.pest", "KPilotPointParameter"),
    "LAK": ("simple_modflow.modflow.mf6.simulation.packages", "LAK"),
    "LakeAreaVolumeRelationship": ("simple_modflow.modflow.mf6.lakes", "LakeAreaVolumeRelationship"),
    "LakeStageObservationSpec": ("simple_modflow.modflow.mf6.pest", "LakeStageObservationSpec"),
    "SfrStageTargets": ("simple_modflow.modflow.mf6.observations", "SfrStageTargets"),
    "SfrFlowTargets": ("simple_modflow.modflow.mf6.observations", "SfrFlowTargets"),
    "SfrStageObservationSpec": ("simple_modflow.modflow.mf6.pest", "SfrStageObservationSpec"),
    "SfrFlowObservationSpec": ("simple_modflow.modflow.mf6.pest", "SfrFlowObservationSpec"),
    "DrnFlowObservationSpec": ("simple_modflow.modflow.mf6.pest", "DrnFlowObservationSpec"),
    "MeshBuildProfile": ("simple_modflow.modflow.mf6.grid.triangle", "MeshBuildProfile"),
    "ModelRegion": ("simple_modflow.modflow.mf6.simulation.regions", "ModelRegion"),
    "OutputControl": ("simple_modflow.modflow.mf6.simulation.packages", "OutputControl"),
    "PestProject": ("simple_modflow.modflow.mf6.pest", "PestProject"),
    "PestRunReview": ("simple_modflow.modflow.mf6.pest", "PestRunReview"),
    "PestRunResults": ("simple_modflow.modflow.mf6.pest", "PestRunResults"),
    "PRTProject": ("simple_modflow.modflow.mf6.prt", "PRTProject"),
    "PRTReleasePoints": ("simple_modflow.modflow.mf6.prt", "PRTReleasePoints"),
    "PRTRunResults": ("simple_modflow.modflow.mf6.prt", "PRTRunResults"),
    "Recharge": ("simple_modflow.modflow.mf6.simulation.packages", "Recharge"),
    "RechargeFromShp": ("simple_modflow.modflow.mf6.recharge", "RechargeFromShp"),
    "RCHFromVector": ("simple_modflow.modflow.mf6.recharge", "RCHFromVector"),
    "RegionGroup": ("simple_modflow.modflow.mf6.simulation.regions", "RegionGroup"),
    "RegionRegistry": ("simple_modflow.modflow.mf6.simulation.regions", "RegionRegistry"),
    "SFR": ("simple_modflow.modflow.mf6.sfr", "SFR"),
    "SimpleModel": ("simple_modflow.modflow.mf6.simplemodel", "SimpleModel"),
    "SimpleModelConfig": ("simple_modflow.modflow.mf6.simplemodel", "SimpleModelConfig"),
    "SimulationBase": ("simple_modflow.modflow.mf6.simulation.base", "SimulationBase"),
    "Storage": ("simple_modflow.modflow.mf6.simulation.packages", "Storage"),
    "Wells": ("simple_modflow.modflow.mf6.simulation.packages", "Wells"),
    "irregular_voronoi_grid": ("simple_modflow.modflow.mf6.canonical", "irregular_voronoi_grid"),
    "SurfaceWaterValidationIssue": (
        "simple_modflow.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationIssue",
    ),
    "SurfaceWaterValidationReport": (
        "simple_modflow.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationReport",
    ),
    "TemporalDiscretization": ("simple_modflow.modflow.mf6.simulation.discretization", "TemporalDiscretization"),
    "TriangleGrid": ("simple_modflow.modflow.mf6.grid.triangle", "TriangleGrid"),
    "UZF": ("simple_modflow.modflow.mf6.simulation.packages", "UZF"),
    "UZFPackageData": ("simple_modflow.modflow.mf6.uzf", "UZFPackageData"),
    "VectorParameterSource": ("simple_modflow.modflow.mf6.pest", "VectorParameterSource"),
    "VoronoiGridPlus": ("simple_modflow.modflow.mf6.grid.voronoi", "VoronoiGridPlus"),
    "open_pest_run": ("simple_modflow.modflow.mf6.pest", "open_pest_run"),
    "open_prt_run": ("simple_modflow.modflow.mf6.prt", "open_prt_run"),
    "DrainConductanceParameter": ("simple_modflow.modflow.mf6.pest", "DrainConductanceParameter"),
    "DrainElevationParameter": ("simple_modflow.modflow.mf6.pest", "DrainElevationParameter"),
    "ExpGeoStruct": ("simple_modflow.modflow.mf6.pest", "ExpGeoStruct"),
    "HeadTargetObservationSpec": ("simple_modflow.modflow.mf6.pest", "HeadTargetObservationSpec"),
    "build_simple_model": ("simple_modflow.modflow.mf6.simplemodel", "build_simple_model"),
    "build_canonical_model": ("simple_modflow.modflow.mf6.canonical_example", "build_canonical_model"),
    "canonical_head_signals": ("simple_modflow.modflow.mf6.canonical", "canonical_head_signals"),
    "canonical_feature_signals": ("simple_modflow.modflow.mf6.canonical", "canonical_feature_signals"),
    "canonical_sfr_signals": ("simple_modflow.modflow.mf6.canonical", "canonical_sfr_signals"),
    "canonical_partition_mask": ("simple_modflow.modflow.mf6.canonical", "canonical_partition_mask"),
    "plot_model_cross_section": ("simple_modflow.modflow.mf6.cross_section_plotting", "plot_model_cross_section"),
    "plot_model_head_map": ("simple_modflow.modflow.mf6.interactive_plotting", "plot_model_head_map"),
    "plot_particle_pathlines": ("simple_modflow.modflow.mf6.interactive_plotting", "plot_particle_pathlines"),
    "export_matplotlib_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_matplotlib_slider_html",
    ),
    "export_cross_section_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_cross_section_slider_html",
    ),
    "export_head_map_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_head_map_slider_html",
    ),
    "export_head_layer_mosaic_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_head_layer_mosaic_slider_html",
    ),
    "build_particle_tracking_scene": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "build_particle_tracking_scene",
    ),
    "export_particle_tracking_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_particle_tracking_html",
    ),
    "validate_surface_water_configuration": (
        "simple_modflow.modflow.mf6.surface_water_validation",
        "validate_surface_water_configuration",
    ),
}

__all__ = sorted([*_SUBMODULES.keys(), *_EXPORTS.keys()])


def __getattr__(name: str):
    if name in _SUBMODULES:
        return import_module(_SUBMODULES[name])

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    return sorted(list(globals().keys()) + __all__)
