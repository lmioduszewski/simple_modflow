"""Public MF6-facing exports for myflopy."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.boundaries import Boundaries
    from myflopy.modflow.mf6.canonical import (
        CANONICAL_MODEL_CONTRACT,
        CanonicalModelContract,
        canonical_feature_signals,
        canonical_head_signals,
        canonical_partition_mask,
        canonical_sfr_signals,
        irregular_voronoi_grid,
    )
    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )
    from myflopy.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from myflopy.modflow.mf6.interactive_plotting import (
        FrameExportProgress,
        ModelMapStyle,
        ModelVisualization,
        ParticleTrackingScene,
        StandaloneHtmlSlider,
    )
    from myflopy.modflow.mf6.observations import (
        DrnFlowTargets,
        HeadTargets,
        LakeStageTargets,
        SfrFlowTargets,
        SfrStageTargets,
    )
    from myflopy.modflow.mf6.parallel import (
        ParallelCompatibilityError,
        ParallelEnvironment,
        ParallelModelWorkflow,
        ParallelSplitResults,
        ParallelSplitRun,
    )
    from myflopy.modflow.mf6.pest import PestProject, PestRunResults, PestRunReview, open_pest_run
    from myflopy.modflow.mf6.prt import (
        ParticleTracking,
        PRTProject,
        PRTReleasePoints,
        PRTRunResults,
        open_prt_run,
    )
    from myflopy.modflow.mf6.simulation.base import SimulationBase
    from myflopy.modflow.mf6.simulation.discretization import (
        DisuGrid,
        DisvGrid,
        TemporalDiscretization,
    )
    from myflopy.modflow.mf6.simulation.packages import (
        CHD,
        UZF,
        Drains,
        InitialConditions,
        KFlow,
        OutputControl,
        Recharge,
        Storage,
        Wells,
    )
    from myflopy.modflow.mf6.recharge import RCHBuilder


_SUBMODULES = {
    "grid": "myflopy.modflow.mf6.grid",
    "pest": "myflopy.modflow.mf6.pest",
    "simulation": "myflopy.modflow.mf6.simulation",
}

_EXPORTS = {
    "Boundaries": ("myflopy.modflow.mf6.boundaries", "Boundaries"),
    "CHD": ("myflopy.modflow.mf6.simulation.packages", "CHD"),
    "CANONICAL_MODEL_CONTRACT": ("myflopy.modflow.mf6.canonical", "CANONICAL_MODEL_CONTRACT"),
    "CanonicalModelContract": ("myflopy.modflow.mf6.canonical", "CanonicalModelContract"),
    "CanonicalModelConfig": ("myflopy.modflow.mf6.canonical_example", "CanonicalModelConfig"),
    "ModelCrossSectionStyle": ("myflopy.modflow.mf6.cross_section_plotting", "ModelCrossSectionStyle"),
    "ModelMapStyle": ("myflopy.modflow.mf6.interactive_plotting", "ModelMapStyle"),
    "FrameExportProgress": ("myflopy.modflow.mf6.interactive_plotting", "FrameExportProgress"),
    "ModelVisualization": ("myflopy.modflow.mf6.interactive_plotting", "ModelVisualization"),
    "ParticleTrackingScene": ("myflopy.modflow.mf6.interactive_plotting", "ParticleTrackingScene"),
    "ParticleTracking": ("myflopy.modflow.mf6.prt", "ParticleTracking"),
    "ParallelCompatibilityError": ("myflopy.modflow.mf6.parallel", "ParallelCompatibilityError"),
    "ParallelEnvironment": ("myflopy.modflow.mf6.parallel", "ParallelEnvironment"),
    "ParallelModelWorkflow": ("myflopy.modflow.mf6.parallel", "ParallelModelWorkflow"),
    "ParallelSplitResults": ("myflopy.modflow.mf6.parallel", "ParallelSplitResults"),
    "ParallelSplitRun": ("myflopy.modflow.mf6.parallel", "ParallelSplitRun"),
    "StandaloneHtmlSlider": ("myflopy.modflow.mf6.interactive_plotting", "StandaloneHtmlSlider"),
    "DRN": ("myflopy.modflow.mf6.drn", "DRN"),
    "DisuGrid": ("myflopy.modflow.mf6.simulation.discretization", "DisuGrid"),
    "DisvGrid": ("myflopy.modflow.mf6.simulation.discretization", "DisvGrid"),
    "Drains": ("myflopy.modflow.mf6.simulation.packages", "Drains"),
    "GHB": ("myflopy.modflow.mf6.ghb", "GHB"),
    "HeadTargets": ("myflopy.modflow.mf6.observations", "HeadTargets"),
    "DrnFlowTargets": ("myflopy.modflow.mf6.observations", "DrnFlowTargets"),
    "LakeStageTargets": ("myflopy.modflow.mf6.observations", "LakeStageTargets"),
    "InitialConditions": ("myflopy.modflow.mf6.simulation.packages", "InitialConditions"),
    "KFlow": ("myflopy.modflow.mf6.simulation.packages", "KFlow"),
    "KPilotPointParameter": ("myflopy.modflow.mf6.pest", "KPilotPointParameter"),
    "LAKBuilder": ("myflopy.modflow.mf6.lakes", "LAKBuilder"),
    "LakeConnection": ("myflopy.modflow.mf6.lakes", "LakeConnection"),
    "LakeOutlet": ("myflopy.modflow.mf6.lakes", "LakeOutlet"),
    "LakeTable": ("myflopy.modflow.mf6.lakes", "LakeTable"),
    "LakeTableBuilder": ("myflopy.modflow.mf6.lakes", "LakeTableBuilder"),
    "LakeStageObservationSpec": ("myflopy.modflow.mf6.pest", "LakeStageObservationSpec"),
    "MVRBuilder": ("myflopy.modflow.mf6.mvr", "MVRBuilder"),
    "Move": ("myflopy.modflow.mf6.mvr", "Move"),
    "MoverConnection": ("myflopy.modflow.mf6.mvr", "MoverConnection"),
    "SfrStageTargets": ("myflopy.modflow.mf6.observations", "SfrStageTargets"),
    "SfrFlowTargets": ("myflopy.modflow.mf6.observations", "SfrFlowTargets"),
    "SfrStageObservationSpec": ("myflopy.modflow.mf6.pest", "SfrStageObservationSpec"),
    "SfrFlowObservationSpec": ("myflopy.modflow.mf6.pest", "SfrFlowObservationSpec"),
    "DrnFlowObservationSpec": ("myflopy.modflow.mf6.pest", "DrnFlowObservationSpec"),
    "MeshBuildProfile": ("myflopy.modflow.mf6.grid.triangle", "MeshBuildProfile"),
    "ModelRegion": ("myflopy.modflow.mf6.simulation.regions", "ModelRegion"),
    "OutputControl": ("myflopy.modflow.mf6.simulation.packages", "OutputControl"),
    "PestProject": ("myflopy.modflow.mf6.pest", "PestProject"),
    "PestRunReview": ("myflopy.modflow.mf6.pest", "PestRunReview"),
    "PestRunResults": ("myflopy.modflow.mf6.pest", "PestRunResults"),
    "PRTProject": ("myflopy.modflow.mf6.prt", "PRTProject"),
    "PRTReleasePoints": ("myflopy.modflow.mf6.prt", "PRTReleasePoints"),
    "PRTRunResults": ("myflopy.modflow.mf6.prt", "PRTRunResults"),
    "Recharge": ("myflopy.modflow.mf6.simulation.packages", "Recharge"),
    "RCHBuilder": ("myflopy.modflow.mf6.recharge", "RCHBuilder"),
    "RechargeFromShp": ("myflopy.modflow.mf6.recharge", "RechargeFromShp"),
    "RegionGroup": ("myflopy.modflow.mf6.simulation.regions", "RegionGroup"),
    "RegionRegistry": ("myflopy.modflow.mf6.simulation.regions", "RegionRegistry"),
    "SFRBuilder": ("myflopy.modflow.mf6.sfr", "SFRBuilder"),
    "StreamConnection": ("myflopy.modflow.mf6.sfr", "StreamConnection"),
    "StreamDiversion": ("myflopy.modflow.mf6.sfr", "StreamDiversion"),
    "StreamNetwork": ("myflopy.modflow.mf6.sfr", "StreamNetwork"),
    "SimpleModelConfig": ("myflopy.modflow.mf6.simplemodel", "SimpleModelConfig"),
    "SimulationBase": ("myflopy.modflow.mf6.simulation.base", "SimulationBase"),
    "Storage": ("myflopy.modflow.mf6.simulation.packages", "Storage"),
    "Wells": ("myflopy.modflow.mf6.simulation.packages", "Wells"),
    "irregular_voronoi_grid": ("myflopy.modflow.mf6.canonical", "irregular_voronoi_grid"),
    "SurfaceWaterValidationIssue": (
        "myflopy.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationIssue",
    ),
    "SurfaceWaterValidationReport": (
        "myflopy.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationReport",
    ),
    "TemporalDiscretization": ("myflopy.modflow.mf6.simulation.discretization", "TemporalDiscretization"),
    "TriangleGrid": ("myflopy.modflow.mf6.grid.triangle", "TriangleGrid"),
    "UZF": ("myflopy.modflow.mf6.simulation.packages", "UZF"),
    "UZFBuilder": ("myflopy.modflow.mf6.uzf", "UZFBuilder"),
    "VectorParameterSource": ("myflopy.modflow.mf6.pest", "VectorParameterSource"),
    "VoronoiGridPlus": ("myflopy.modflow.mf6.grid.voronoi", "VoronoiGridPlus"),
    "open_pest_run": ("myflopy.modflow.mf6.pest", "open_pest_run"),
    "open_prt_run": ("myflopy.modflow.mf6.prt", "open_prt_run"),
    "DrainConductanceParameter": ("myflopy.modflow.mf6.pest", "DrainConductanceParameter"),
    "DrainElevationParameter": ("myflopy.modflow.mf6.pest", "DrainElevationParameter"),
    "ExpGeoStruct": ("myflopy.modflow.mf6.pest", "ExpGeoStruct"),
    "HeadTargetObservationSpec": ("myflopy.modflow.mf6.pest", "HeadTargetObservationSpec"),
    "simple_model_spec": ("myflopy.modflow.mf6.simplemodel", "simple_model_spec"),
    "build_canonical_model": ("myflopy.modflow.mf6.canonical_example", "build_canonical_model"),
    "canonical_head_signals": ("myflopy.modflow.mf6.canonical", "canonical_head_signals"),
    "canonical_feature_signals": ("myflopy.modflow.mf6.canonical", "canonical_feature_signals"),
    "canonical_sfr_signals": ("myflopy.modflow.mf6.canonical", "canonical_sfr_signals"),
    "canonical_partition_mask": ("myflopy.modflow.mf6.canonical", "canonical_partition_mask"),
    "plot_model_cross_section": ("myflopy.modflow.mf6.cross_section_plotting", "plot_model_cross_section"),
    "plot_model_head_map": ("myflopy.modflow.mf6.interactive_plotting", "plot_model_head_map"),
    "plot_particle_pathlines": ("myflopy.modflow.mf6.interactive_plotting", "plot_particle_pathlines"),
    "export_matplotlib_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_matplotlib_slider_html",
    ),
    "export_cross_section_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_cross_section_slider_html",
    ),
    "export_head_map_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_head_map_slider_html",
    ),
    "export_head_layer_mosaic_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_head_layer_mosaic_slider_html",
    ),
    "build_particle_tracking_scene": (
        "myflopy.modflow.mf6.interactive_plotting",
        "build_particle_tracking_scene",
    ),
    "export_particle_tracking_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_particle_tracking_html",
    ),
    "validate_surface_water_configuration": (
        "myflopy.modflow.mf6.surface_water_validation",
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
