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
    from simple_modflow.modflow.mf6.recharge import RCHFromVector
    from simple_modflow.modflow.mf6.surface_water_validation import (
        SurfaceWaterValidationIssue,
        SurfaceWaterValidationReport,
        validate_surface_water_configuration,
    )
    from simple_modflow.modflow.mf6.simplemodel import SimpleModel, SimpleModelConfig, build_simple_model
    from simple_modflow.modflow.mp3du import ParticleTrackingInput, prepare_particle_tracking, run_particle_tracking
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase
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
    "VoronoiGridPlus": ("simple_modflow.modflow.mf6.grid.voronoi", "VoronoiGridPlus"),
    "build_simple_model": ("simple_modflow.modflow.mf6.simplemodel", "build_simple_model"),
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
