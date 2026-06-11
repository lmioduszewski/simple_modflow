"""Public package surface for simple_modflow."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    # These imports are for static analysis and IDE completion only. Runtime
    # access still goes through the lazy ``_EXPORTS`` / ``__getattr__`` path.
    from simple_modflow.modflow.mf6.drn import DRNFromVector
    from simple_modflow.modflow.mf6.ghb import GHBFromVector
    from simple_modflow.modflow.mf6.grid.triangle import MeshBuildProfile, TriangleGrid
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from simple_modflow.modflow.mf6.chd import CHDFromVector
    from simple_modflow.modflow.mf6.kflow import KFromVector
    from simple_modflow.modflow.mf6.observations import (
        DrnFlowTargets,
        HeadTargets,
        LakeStageTargets,
        SfrFlowTargets,
        SfrStageTargets,
    )
    from simple_modflow.modflow.mf6.interactive_plotting import (
        FrameExportProgress,
        ModelMapStyle,
        ModelVisualization,
        ParticleTrackingScene,
        StandaloneHtmlSlider,
        build_particle_tracking_scene,
        export_cross_section_slider_html,
        export_head_layer_mosaic_slider_html,
        export_head_map_slider_html,
        export_matplotlib_slider_html,
        export_particle_tracking_html,
        plot_model_head_map,
        plot_particle_pathlines,
    )
    from simple_modflow.modflow.mf6.pest import (
        DrainConductanceParameter,
        DrainElevationParameter,
        DrnFlowObservationSpec,
        ExpGeoStruct,
        HeadTargetObservationSpec,
        KPilotPointParameter,
        LakeStageObservationSpec,
        PestProject,
        PestRunResults,
        PestRunReview,
        SfrFlowObservationSpec,
        SfrStageObservationSpec,
        VectorParameterSource,
        open_pest_run,
    )
    from simple_modflow.modflow.mf6.prt import (
        ParticleTracking,
        PRTProject,
        PRTReleasePoints,
        PRTRunResults,
        open_prt_run,
    )
    from simple_modflow.modflow.mf6.parallel import (
        ParallelCompatibilityError,
        ParallelEnvironment,
        ParallelModelWorkflow,
        ParallelSplitResults,
        ParallelSplitRun,
    )
    from simple_modflow.modflow.mf6.canonical import (
        CANONICAL_MODEL_CONTRACT,
        CanonicalModelContract,
        canonical_feature_signals,
        canonical_head_signals,
        canonical_partition_mask,
        canonical_sfr_signals,
    )
    from simple_modflow.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )
    from simple_modflow.modflow.mf6.recharge import RCHFromVector
    from simple_modflow.modflow.mf6.surface_water_validation import (
        SurfaceWaterValidationIssue,
        SurfaceWaterValidationReport,
        validate_surface_water_configuration,
    )
    from simple_modflow.modflow.mf6.simplemodel import SimpleModel, SimpleModelConfig, build_simple_model
    from simple_modflow.modflow.mp3du import ParticleTrackingInput, prepare_particle_tracking, run_particle_tracking
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase
    from simple_modflow.modflow.mf6.simulation.packages import Wells
    from simple_modflow.modflow.utils.datatypes.readers import read_gpkg, read_shp_gpkg
    from simple_modflow.project import (
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

try:
    from importlib.metadata import version

    __version__ = version("simple_modflow")
except Exception:
    __version__ = "0.1.0"


_EXPORTS = {
    "DiscoveredRun": ("simple_modflow.project", "DiscoveredRun"),
    "LoadedMf6Run": ("simple_modflow.project", "LoadedMf6Run"),
    "ModelSpec": ("simple_modflow.project", "ModelSpec"),
    "ModelGroup": ("simple_modflow.project", "ModelGroup"),
    "PackageArtifact": ("simple_modflow.project", "PackageArtifact"),
    "PackageCompatibilityError": ("simple_modflow.project", "PackageCompatibilityError"),
    "ParallelCompatibilityError": ("simple_modflow.modflow.mf6.parallel", "ParallelCompatibilityError"),
    "ParallelEnvironment": ("simple_modflow.modflow.mf6.parallel", "ParallelEnvironment"),
    "ParallelModelWorkflow": ("simple_modflow.modflow.mf6.parallel", "ParallelModelWorkflow"),
    "ParallelSplitResults": ("simple_modflow.modflow.mf6.parallel", "ParallelSplitResults"),
    "ParallelSplitRun": ("simple_modflow.modflow.mf6.parallel", "ParallelSplitRun"),
    "CANONICAL_MODEL_CONTRACT": ("simple_modflow.modflow.mf6.canonical", "CANONICAL_MODEL_CONTRACT"),
    "CanonicalModelContract": ("simple_modflow.modflow.mf6.canonical", "CanonicalModelContract"),
    "CanonicalModelConfig": ("simple_modflow.modflow.mf6.canonical_example", "CanonicalModelConfig"),
    "ProjectCatalog": ("simple_modflow.project", "ProjectCatalog"),
    "RunComparison": ("simple_modflow.project", "RunComparison"),
    "RunExplorer": ("simple_modflow.project", "RunExplorer"),
    "RunRecord": ("simple_modflow.project", "RunRecord"),
    "RunLoader": ("simple_modflow.project", "RunLoader"),
    "RunSpec": ("simple_modflow.project", "RunSpec"),
    "discover_existing_runs": ("simple_modflow.project", "discover_existing_runs"),
    "explore_runs": ("simple_modflow.project", "explore_runs"),
    "import_run_archive": ("simple_modflow.project", "import_run_archive"),
    "load_mf6_run": ("simple_modflow.project", "load_mf6_run"),
    "MeshBuildProfile": ("simple_modflow.modflow.mf6.grid.triangle", "MeshBuildProfile"),
    "CHDFromVector": ("simple_modflow.modflow.mf6.chd", "CHDFromVector"),
    "DRNFromVector": ("simple_modflow.modflow.mf6.drn", "DRNFromVector"),
    "DrnFlowTargets": ("simple_modflow.modflow.mf6.observations", "DrnFlowTargets"),
    "GHBFromVector": ("simple_modflow.modflow.mf6.ghb", "GHBFromVector"),
    "HeadTargets": ("simple_modflow.modflow.mf6.observations", "HeadTargets"),
    "ModelMapStyle": ("simple_modflow.modflow.mf6.interactive_plotting", "ModelMapStyle"),
    "FrameExportProgress": ("simple_modflow.modflow.mf6.interactive_plotting", "FrameExportProgress"),
    "ModelVisualization": ("simple_modflow.modflow.mf6.interactive_plotting", "ModelVisualization"),
    "ParticleTrackingScene": ("simple_modflow.modflow.mf6.interactive_plotting", "ParticleTrackingScene"),
    "ParticleTracking": ("simple_modflow.modflow.mf6.prt", "ParticleTracking"),
    "StandaloneHtmlSlider": ("simple_modflow.modflow.mf6.interactive_plotting", "StandaloneHtmlSlider"),
    "LakeStageTargets": ("simple_modflow.modflow.mf6.observations", "LakeStageTargets"),
    "SfrStageTargets": ("simple_modflow.modflow.mf6.observations", "SfrStageTargets"),
    "SfrFlowTargets": ("simple_modflow.modflow.mf6.observations", "SfrFlowTargets"),
    "KFromVector": ("simple_modflow.modflow.mf6.kflow", "KFromVector"),
    "KPilotPointParameter": ("simple_modflow.modflow.mf6.pest", "KPilotPointParameter"),
    "patch_simulation_plot": ("simple_modflow.project", "patch_simulation_plot"),
    "DrainConductanceParameter": ("simple_modflow.modflow.mf6.pest", "DrainConductanceParameter"),
    "DrainElevationParameter": ("simple_modflow.modflow.mf6.pest", "DrainElevationParameter"),
    "DrnFlowObservationSpec": ("simple_modflow.modflow.mf6.pest", "DrnFlowObservationSpec"),
    "ExpGeoStruct": ("simple_modflow.modflow.mf6.pest", "ExpGeoStruct"),
    "HeadTargetObservationSpec": ("simple_modflow.modflow.mf6.pest", "HeadTargetObservationSpec"),
    "LakeStageObservationSpec": ("simple_modflow.modflow.mf6.pest", "LakeStageObservationSpec"),
    "SfrStageObservationSpec": ("simple_modflow.modflow.mf6.pest", "SfrStageObservationSpec"),
    "SfrFlowObservationSpec": ("simple_modflow.modflow.mf6.pest", "SfrFlowObservationSpec"),
    "PestProject": ("simple_modflow.modflow.mf6.pest", "PestProject"),
    "PestRunReview": ("simple_modflow.modflow.mf6.pest", "PestRunReview"),
    "PestRunResults": ("simple_modflow.modflow.mf6.pest", "PestRunResults"),
    "PRTProject": ("simple_modflow.modflow.mf6.prt", "PRTProject"),
    "PRTReleasePoints": ("simple_modflow.modflow.mf6.prt", "PRTReleasePoints"),
    "PRTRunResults": ("simple_modflow.modflow.mf6.prt", "PRTRunResults"),
    "RCHFromVector": ("simple_modflow.modflow.mf6.recharge", "RCHFromVector"),
    "SimpleModel": ("simple_modflow.modflow.mf6.simplemodel", "SimpleModel"),
    "SimpleModelConfig": ("simple_modflow.modflow.mf6.simplemodel", "SimpleModelConfig"),
    "SimulationBase": ("simple_modflow.modflow.mf6.simulation.base", "SimulationBase"),
    "SurfaceWaterValidationIssue": (
        "simple_modflow.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationIssue",
    ),
    "SurfaceWaterValidationReport": (
        "simple_modflow.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationReport",
    ),
    "TriangleGrid": ("simple_modflow.modflow.mf6.grid.triangle", "TriangleGrid"),
    "VectorParameterSource": ("simple_modflow.modflow.mf6.pest", "VectorParameterSource"),
    "open_pest_run": ("simple_modflow.modflow.mf6.pest", "open_pest_run"),
    "open_prt_run": ("simple_modflow.modflow.mf6.prt", "open_prt_run"),
    "VoronoiGridPlus": ("simple_modflow.modflow.mf6.grid.voronoi", "VoronoiGridPlus"),
    "build_simple_model": ("simple_modflow.modflow.mf6.simplemodel", "build_simple_model"),
    "build_canonical_model": ("simple_modflow.modflow.mf6.canonical_example", "build_canonical_model"),
    "canonical_head_signals": ("simple_modflow.modflow.mf6.canonical", "canonical_head_signals"),
    "canonical_feature_signals": ("simple_modflow.modflow.mf6.canonical", "canonical_feature_signals"),
    "canonical_sfr_signals": ("simple_modflow.modflow.mf6.canonical", "canonical_sfr_signals"),
    "canonical_partition_mask": ("simple_modflow.modflow.mf6.canonical", "canonical_partition_mask"),
    "build_particle_tracking_scene": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "build_particle_tracking_scene",
    ),
    "export_cross_section_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_cross_section_slider_html",
    ),
    "export_head_layer_mosaic_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_head_layer_mosaic_slider_html",
    ),
    "export_head_map_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_head_map_slider_html",
    ),
    "export_matplotlib_slider_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_matplotlib_slider_html",
    ),
    "export_particle_tracking_html": (
        "simple_modflow.modflow.mf6.interactive_plotting",
        "export_particle_tracking_html",
    ),
    "plot_model_head_map": ("simple_modflow.modflow.mf6.interactive_plotting", "plot_model_head_map"),
    "plot_particle_pathlines": ("simple_modflow.modflow.mf6.interactive_plotting", "plot_particle_pathlines"),
    "ParticleTrackingInput": ("simple_modflow.modflow.mp3du", "ParticleTrackingInput"),
    "prepare_particle_tracking": ("simple_modflow.modflow.mp3du", "prepare_particle_tracking"),
    "run_particle_tracking": ("simple_modflow.modflow.mp3du", "run_particle_tracking"),
    "modflow": ("simple_modflow", "modflow"),
    "project": ("simple_modflow", "project"),
    "read_gpkg": ("simple_modflow.modflow.utils.datatypes.readers", "read_gpkg"),
    "read_shp_gpkg": ("simple_modflow.modflow.utils.datatypes.readers", "read_shp_gpkg"),
    "summarize_discovered_runs": ("simple_modflow.project", "summarize_discovered_runs"),
    "validate_surface_water_configuration": (
        "simple_modflow.modflow.mf6.surface_water_validation",
        "validate_surface_water_configuration",
    ),
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    if name == "modflow":
        return import_module("simple_modflow.modflow")
    if name == "project":
        return import_module("simple_modflow.project")

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    return sorted(list(globals().keys()) + __all__)
