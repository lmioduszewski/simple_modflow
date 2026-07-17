"""Public package surface for myflopy.

Package-first helpers such as ``mf.lak(...)`` and ``mf.rch(...)`` are the
preferred top-level API. Lower-level builders and ``*_spec`` factories remain
available as engine imports (``__engine__``), but they are intentionally
omitted from ``__all__`` and ``dir(myflopy)`` so discovery points at the
package-first path. Deprecated names are tracked in ``__compatibility__``
(see ``docs/deprecation_policy.md``).
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    # These imports are for static analysis and IDE completion only. Runtime
    # access still goes through the lazy ``_EXPORTS`` / ``__getattr__`` path.
    from myflopy.advanced import (
        chd_spec,
        drn_spec,
        ghb_spec,
        lak_spec,
        mvr_spec,
        rch_spec,
        sfr_spec,
        uzf_spec,
        wel_spec,
    )
    from myflopy.builders import (
        PackageBuilder,
        build_gwf_gwe_exchange,
        build_gwf_gwf_exchange,
        build_gwf_gwt_exchange,
        build_gwf_prt_exchange,
        build_ims,
    )
    from myflopy.geopackage import CellSurfaceOffset, GeoPackageSource
    from myflopy.layers import (
        LayerBuildResult,
        LayerQCReport,
        LayerStack,
        modflow_surfaces,
    )
    from myflopy.modflow.mf6.canonical import (
        CANONICAL_MODEL_CONTRACT,
        CanonicalModelContract,
        canonical_feature_signals,
        canonical_head_signals,
        canonical_partition_mask,
        canonical_sfr_signals,
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
        build_particle_tracking_scene,
        export_cross_section_slider_html,
        export_head_layer_mosaic_slider_html,
        export_head_map_slider_html,
        export_matplotlib_slider_html,
        export_particle_tracking_html,
        plot_model_head_map,
        plot_particle_pathlines,
    )
    from myflopy.modflow.mf6.lakes import (
        LAKBuilder,
        LakeConnection,
        LakeOutlet,
        LakeTable,
        LakeTableBuilder,
    )
    from myflopy.modflow.mf6.mvr import Move, MoverConnection, MVRBuilder
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
    from myflopy.modflow.mf6.pest import (
        DrnFlowObservationSpec,
        ExpGeoStruct,
        HeadTargetObservationSpec,
        LakeStageObservationSpec,
        PestProject,
        PestRunHandle,
        SfrFlowObservationSpec,
        SfrStageObservationSpec,
        VectorParameterSource,
        find_pest_runs,
    )
    from myflopy.modflow.mf6.prt import (
        ParticleTracking,
        PRTProject,
        PRTReleasePoints,
        PRTRunResults,
        open_prt_run,
    )
    from myflopy.modflow.mf6.recharge import RCHBuilder
    from myflopy.modflow.mf6.sfr import (
        SFRBuilder,
        StreamConnection,
        StreamDiversion,
        StreamNetwork,
    )
    from myflopy.modflow.mf6.simplemodel import SimpleModelConfig, simple_model_spec
    from myflopy.modflow.mf6.simulation.base import SimulationBase
    from myflopy.modflow.mf6.simulation.packages import Wells
    from myflopy.modflow.mf6.surface_water_validation import (
        SurfaceWaterValidationIssue,
        SurfaceWaterValidationReport,
        validate_surface_water_configuration,
    )
    from myflopy.modflow.mf6.uzf import UZFBuilder
    from myflopy.modflow.mp3du import (
        ParticleTrackingInput,
        prepare_particle_tracking,
        run_particle_tracking,
    )
    from myflopy.modflow.utils.datatypes.hover import HoverSpec, HoverStyle
    from myflopy.modflow.utils.datatypes.readers import read_gpkg, read_shp_gpkg
    from myflopy.package_api import (
        chd,
        disv,
        drn,
        ghb,
        gwe,
        gwf,
        gwt,
        ic,
        ims,
        lak,
        mvr,
        npf,
        oc,
        prt,
        rch,
        sfr,
        simulation,
        sto,
        tdis,
        uzf,
        wel,
    )
    from myflopy.project import (
        LoadedMf6Run,
        ModelGroup,
        PackageArtifact,
        PackageCompatibilityError,
        load_mf6_run,
        patch_simulation_plot,
    )
    from myflopy.sources import (
        DataSourceSpec,
        GeoPackageSourceSpec,
        LiteralSource,
        RasterSource,
        ShapeSource,
        TableSource,
    )
    from myflopy.specs import (
        BuiltModel,
        BuiltSimulation,
        ExchangeSpec,
        GridRef,
        GridSpec,
        GweModel,
        GwfModel,
        GwtModel,
        Mf6Model,
        Mf6Simulation,
        ModelContext,
        ModelSpec,
        ModelType,
        PackageRef,
        PackageSpec,
        PostBuildHook,
        PrtModel,
        SimulationSpec,
        SpecBuildContext,
        grid_ref,
        ref,
    )
    from myflopy.surfaces import LayerSurfaces, Surface
    from myflopy.workspace import ModelView, Project, ProjectLayout, Run, load_run

try:
    from importlib.metadata import version

    __version__ = version("myflopy")
except Exception:
    __version__ = "0.1.0"


_EXPORTS = {
    "LoadedMf6Run": ("myflopy.project", "LoadedMf6Run"),
    "BuiltModel": ("myflopy.specs", "BuiltModel"),
    "BuiltSimulation": ("myflopy.specs", "BuiltSimulation"),
    "GweModel": ("myflopy.specs", "GweModel"),
    "GwfModel": ("myflopy.specs", "GwfModel"),
    "GwtModel": ("myflopy.specs", "GwtModel"),
    "Mf6Model": ("myflopy.specs", "Mf6Model"),
    "Mf6Simulation": ("myflopy.specs", "Mf6Simulation"),
    "PrtModel": ("myflopy.specs", "PrtModel"),
    "build_gwf_gwe_exchange": ("myflopy.builders", "build_gwf_gwe_exchange"),
    "build_gwf_gwf_exchange": ("myflopy.builders", "build_gwf_gwf_exchange"),
    "build_gwf_gwt_exchange": ("myflopy.builders", "build_gwf_gwt_exchange"),
    "build_gwf_prt_exchange": ("myflopy.builders", "build_gwf_prt_exchange"),
    "build_ims": ("myflopy.builders", "build_ims"),
    "PackageBuilder": ("myflopy.builders", "PackageBuilder"),
    "ExchangeSpec": ("myflopy.specs", "ExchangeSpec"),
    "ModelContext": ("myflopy.specs", "ModelContext"),
    "ModelSpec": ("myflopy.specs", "ModelSpec"),
    "ModelType": ("myflopy.specs", "ModelType"),
    "GridSpec": ("myflopy.specs", "GridSpec"),
    "GridRef": ("myflopy.specs", "GridRef"),
    "PackageRef": ("myflopy.specs", "PackageRef"),
    "PackageSpec": ("myflopy.specs", "PackageSpec"),
    "ref": ("myflopy.specs", "ref"),
    "grid_ref": ("myflopy.specs", "grid_ref"),
    "PostBuildHook": ("myflopy.specs", "PostBuildHook"),
    "SpecBuildContext": ("myflopy.specs", "SpecBuildContext"),
    "DataSourceSpec": ("myflopy.sources", "DataSourceSpec"),
    "GeoPackageSourceSpec": ("myflopy.sources", "GeoPackageSourceSpec"),
    "LiteralSource": ("myflopy.sources", "LiteralSource"),
    "RasterSource": ("myflopy.sources", "RasterSource"),
    "ShapeSource": ("myflopy.sources", "ShapeSource"),
    "TableSource": ("myflopy.sources", "TableSource"),
    "CellSurfaceOffset": ("myflopy.geopackage", "CellSurfaceOffset"),
    "GeoPackageSource": ("myflopy.geopackage", "GeoPackageSource"),
    "Surface": ("myflopy.surfaces", "Surface"),
    "LayerSurfaces": ("myflopy.surfaces", "LayerSurfaces"),
    "LayerStack": ("myflopy.layers", "LayerStack"),
    "LayerBuildResult": ("myflopy.layers", "LayerBuildResult"),
    "LayerQCReport": ("myflopy.layers", "LayerQCReport"),
    "modflow_surfaces": ("myflopy.layers", "modflow_surfaces"),
    "ProjectLayout": ("myflopy.workspace", "ProjectLayout"),
    "ModelView": ("myflopy.workspace", "ModelView"),
    "Project": ("myflopy.workspace", "Project"),
    "Run": ("myflopy.workspace", "Run"),
    "load_run": ("myflopy.workspace", "load_run"),
    "SimulationSpec": ("myflopy.specs", "SimulationSpec"),
    "SFRBuilder": ("myflopy.modflow.mf6.sfr", "SFRBuilder"),
    "LAKBuilder": ("myflopy.modflow.mf6.lakes", "LAKBuilder"),
    "LakeConnection": ("myflopy.modflow.mf6.lakes", "LakeConnection"),
    "LakeOutlet": ("myflopy.modflow.mf6.lakes", "LakeOutlet"),
    "LakeTable": ("myflopy.modflow.mf6.lakes", "LakeTable"),
    "LakeTableBuilder": ("myflopy.modflow.mf6.lakes", "LakeTableBuilder"),
    "MVRBuilder": ("myflopy.modflow.mf6.mvr", "MVRBuilder"),
    "Move": ("myflopy.modflow.mf6.mvr", "Move"),
    "MoverConnection": ("myflopy.modflow.mf6.mvr", "MoverConnection"),
    "StreamConnection": ("myflopy.modflow.mf6.sfr", "StreamConnection"),
    "StreamDiversion": ("myflopy.modflow.mf6.sfr", "StreamDiversion"),
    "StreamNetwork": ("myflopy.modflow.mf6.sfr", "StreamNetwork"),
    "RCHBuilder": ("myflopy.modflow.mf6.recharge", "RCHBuilder"),
    "UZFBuilder": ("myflopy.modflow.mf6.uzf", "UZFBuilder"),
    "Wells": ("myflopy.modflow.mf6.simulation.packages", "Wells"),
    "simulation": ("myflopy.package_api", "simulation"),
    "tdis": ("myflopy.package_api", "tdis"),
    "ims": ("myflopy.package_api", "ims"),
    "disv": ("myflopy.package_api", "disv"),
    "ic": ("myflopy.package_api", "ic"),
    "npf": ("myflopy.package_api", "npf"),
    "sto": ("myflopy.package_api", "sto"),
    "oc": ("myflopy.package_api", "oc"),
    "chd": ("myflopy.package_api", "chd"),
    "ghb": ("myflopy.package_api", "ghb"),
    "gwe": ("myflopy.package_api", "gwe"),
    "gwf": ("myflopy.package_api", "gwf"),
    "gwt": ("myflopy.package_api", "gwt"),
    "drn": ("myflopy.package_api", "drn"),
    "wel": ("myflopy.package_api", "wel"),
    "prt": ("myflopy.package_api", "prt"),
    "rch": ("myflopy.package_api", "rch"),
    "uzf": ("myflopy.package_api", "uzf"),
    "sfr": ("myflopy.package_api", "sfr"),
    "lak": ("myflopy.package_api", "lak"),
    "mvr": ("myflopy.package_api", "mvr"),
    "sfr_connection": ("myflopy.package_api", "sfr_connection"),
    "lak_connection": ("myflopy.package_api", "lak_connection"),
    "chd_spec": ("myflopy.advanced", "chd_spec"),
    "drn_spec": ("myflopy.advanced", "drn_spec"),
    "ghb_spec": ("myflopy.advanced", "ghb_spec"),
    "lak_spec": ("myflopy.advanced", "lak_spec"),
    "mvr_spec": ("myflopy.advanced", "mvr_spec"),
    "rch_spec": ("myflopy.advanced", "rch_spec"),
    "sfr_spec": ("myflopy.advanced", "sfr_spec"),
    "uzf_spec": ("myflopy.advanced", "uzf_spec"),
    "wel_spec": ("myflopy.advanced", "wel_spec"),
    "ModelGroup": ("myflopy.project", "ModelGroup"),
    "PackageArtifact": ("myflopy.project", "PackageArtifact"),
    "PackageCompatibilityError": ("myflopy.project", "PackageCompatibilityError"),
    "ParallelCompatibilityError": (
        "myflopy.modflow.mf6.parallel",
        "ParallelCompatibilityError",
    ),
    "ParallelEnvironment": ("myflopy.modflow.mf6.parallel", "ParallelEnvironment"),
    "ParallelModelWorkflow": ("myflopy.modflow.mf6.parallel", "ParallelModelWorkflow"),
    "ParallelSplitResults": ("myflopy.modflow.mf6.parallel", "ParallelSplitResults"),
    "ParallelSplitRun": ("myflopy.modflow.mf6.parallel", "ParallelSplitRun"),
    "CANONICAL_MODEL_CONTRACT": (
        "myflopy.modflow.mf6.canonical",
        "CANONICAL_MODEL_CONTRACT",
    ),
    "CanonicalModelContract": (
        "myflopy.modflow.mf6.canonical",
        "CanonicalModelContract",
    ),
    "CanonicalModelConfig": (
        "myflopy.modflow.mf6.canonical_example",
        "CanonicalModelConfig",
    ),
    "load_mf6_run": ("myflopy.project", "load_mf6_run"),
    "MeshBuildProfile": ("myflopy.modflow.mf6.grid.triangle", "MeshBuildProfile"),
    "DrnFlowTargets": ("myflopy.modflow.mf6.observations", "DrnFlowTargets"),
    "HeadTargets": ("myflopy.modflow.mf6.observations", "HeadTargets"),
    "ModelMapStyle": ("myflopy.modflow.mf6.interactive_plotting", "ModelMapStyle"),
    "FrameExportProgress": (
        "myflopy.modflow.mf6.interactive_plotting",
        "FrameExportProgress",
    ),
    "ModelVisualization": (
        "myflopy.modflow.mf6.interactive_plotting",
        "ModelVisualization",
    ),
    "ParticleTrackingScene": (
        "myflopy.modflow.mf6.interactive_plotting",
        "ParticleTrackingScene",
    ),
    "ParticleTracking": ("myflopy.modflow.mf6.prt", "ParticleTracking"),
    "StandaloneHtmlSlider": (
        "myflopy.modflow.mf6.interactive_plotting",
        "StandaloneHtmlSlider",
    ),
    "LakeStageTargets": ("myflopy.modflow.mf6.observations", "LakeStageTargets"),
    "SfrStageTargets": ("myflopy.modflow.mf6.observations", "SfrStageTargets"),
    "SfrFlowTargets": ("myflopy.modflow.mf6.observations", "SfrFlowTargets"),
    "patch_simulation_plot": ("myflopy.project", "patch_simulation_plot"),
    "DrnFlowObservationSpec": ("myflopy.modflow.mf6.pest", "DrnFlowObservationSpec"),
    "ExpGeoStruct": ("myflopy.modflow.mf6.pest", "ExpGeoStruct"),
    "HeadTargetObservationSpec": (
        "myflopy.modflow.mf6.pest",
        "HeadTargetObservationSpec",
    ),
    "LakeStageObservationSpec": (
        "myflopy.modflow.mf6.pest",
        "LakeStageObservationSpec",
    ),
    "SfrStageObservationSpec": ("myflopy.modflow.mf6.pest", "SfrStageObservationSpec"),
    "SfrFlowObservationSpec": ("myflopy.modflow.mf6.pest", "SfrFlowObservationSpec"),
    "PestProject": ("myflopy.modflow.mf6.pest", "PestProject"),
    "PestRunHandle": ("myflopy.modflow.mf6.pest", "PestRunHandle"),
    "PRTProject": ("myflopy.modflow.mf6.prt", "PRTProject"),
    "PRTReleasePoints": ("myflopy.modflow.mf6.prt", "PRTReleasePoints"),
    "PRTRunResults": ("myflopy.modflow.mf6.prt", "PRTRunResults"),
    "SimpleModelConfig": ("myflopy.modflow.mf6.simplemodel", "SimpleModelConfig"),
    "SimulationBase": ("myflopy.modflow.mf6.simulation.base", "SimulationBase"),
    "SurfaceWaterValidationIssue": (
        "myflopy.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationIssue",
    ),
    "SurfaceWaterValidationReport": (
        "myflopy.modflow.mf6.surface_water_validation",
        "SurfaceWaterValidationReport",
    ),
    "TriangleGrid": ("myflopy.modflow.mf6.grid.triangle", "TriangleGrid"),
    "VectorParameterSource": ("myflopy.modflow.mf6.pest", "VectorParameterSource"),
    "find_pest_runs": ("myflopy.modflow.mf6.pest", "find_pest_runs"),
    "open_prt_run": ("myflopy.modflow.mf6.prt", "open_prt_run"),
    "VoronoiGridPlus": ("myflopy.modflow.mf6.grid.voronoi", "VoronoiGridPlus"),
    "simple_model_spec": ("myflopy.modflow.mf6.simplemodel", "simple_model_spec"),
    "build_canonical_model": (
        "myflopy.modflow.mf6.canonical_example",
        "build_canonical_model",
    ),
    "canonical_head_signals": (
        "myflopy.modflow.mf6.canonical",
        "canonical_head_signals",
    ),
    "canonical_feature_signals": (
        "myflopy.modflow.mf6.canonical",
        "canonical_feature_signals",
    ),
    "canonical_sfr_signals": ("myflopy.modflow.mf6.canonical", "canonical_sfr_signals"),
    "canonical_partition_mask": (
        "myflopy.modflow.mf6.canonical",
        "canonical_partition_mask",
    ),
    "build_particle_tracking_scene": (
        "myflopy.modflow.mf6.interactive_plotting",
        "build_particle_tracking_scene",
    ),
    "export_cross_section_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_cross_section_slider_html",
    ),
    "export_head_layer_mosaic_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_head_layer_mosaic_slider_html",
    ),
    "export_head_map_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_head_map_slider_html",
    ),
    "export_matplotlib_slider_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_matplotlib_slider_html",
    ),
    "export_particle_tracking_html": (
        "myflopy.modflow.mf6.interactive_plotting",
        "export_particle_tracking_html",
    ),
    "plot_model_head_map": (
        "myflopy.modflow.mf6.interactive_plotting",
        "plot_model_head_map",
    ),
    "plot_particle_pathlines": (
        "myflopy.modflow.mf6.interactive_plotting",
        "plot_particle_pathlines",
    ),
    "ParticleTrackingInput": ("myflopy.modflow.mp3du", "ParticleTrackingInput"),
    "prepare_particle_tracking": ("myflopy.modflow.mp3du", "prepare_particle_tracking"),
    "run_particle_tracking": ("myflopy.modflow.mp3du", "run_particle_tracking"),
    "modflow": ("myflopy", "modflow"),
    "project": ("myflopy", "project"),
    "HoverSpec": ("myflopy.modflow.utils.datatypes.hover", "HoverSpec"),
    "HoverStyle": ("myflopy.modflow.utils.datatypes.hover", "HoverStyle"),
    "read_gpkg": ("myflopy.modflow.utils.datatypes.readers", "read_gpkg"),
    "read_shp_gpkg": ("myflopy.modflow.utils.datatypes.readers", "read_shp_gpkg"),
    "validate_surface_water_configuration": (
        "myflopy.modflow.mf6.surface_water_validation",
        "validate_surface_water_configuration",
    ),
}

_SECOND_TIER_EXPORTS = {
    "LAKBuilder",
    "MVRBuilder",
    "PackageBuilder",
    "RCHBuilder",
    "SFRBuilder",
    "UZFBuilder",
    "Wells",
    "build_gwf_gwe_exchange",
    "build_gwf_gwf_exchange",
    "build_gwf_gwt_exchange",
    "build_gwf_prt_exchange",
    "build_ims",
    "chd_spec",
    "drn_spec",
    "ghb_spec",
    "lak_spec",
    "mvr_spec",
    "rch_spec",
    "sfr_spec",
    "uzf_spec",
    "wel_spec",
}

__preferred__ = tuple(
    sorted(name for name in _EXPORTS if name not in _SECOND_TIER_EXPORTS)
)
__engine__ = tuple(sorted(_SECOND_TIER_EXPORTS))
__all__ = list(__preferred__)

# Warned compatibility aliases — the authoritative registry required by
# docs/deprecation_policy.md. Every name here resolves with a
# DeprecationWarning via __getattr__ only (hidden from __all__/dir()/
# TYPE_CHECKING per D12) and may be removed two tagged releases after the
# release that deprecated it. tests/test_deprecation.py cross-checks this
# tuple against the registry in myflopy/_deprecation.py.
__compatibility__ = (
    "myflopy.modflow.mf6.DRN",
    "myflopy.modflow.mf6.GHB",
    "myflopy.modflow.mf6.drn.DRN",
    "myflopy.modflow.mf6.ghb.GHB",
    "myflopy.modflow.mp3du.particles.PRT",
    "myflopy.modflow.mp3du.particles.PrtDisv",
    "myflopy.modflow.mp3du.particles.PrtFmi",
    "myflopy.modflow.mp3du.particles.PrtMip",
    "myflopy.modflow.mp3du.particles.PrtOc",
    "myflopy.modflow.mp3du.particles.PrtPrp",
    "myflopy.project.model_group.ModelGroup.chd",
    "myflopy.project.model_group.ModelGroup.drn",
    "myflopy.project.model_group.ModelGroup.ghb",
    "myflopy.project.model_group.ModelGroup.rch",
    "myflopy.project.model_group.ModelGroup.uzf",
    "myflopy.project.model_group.ModelGroup.wel",
)


def __getattr__(name: str) -> Any:
    """Lazily import a public ``myflopy`` export (or subpackage) on first attribute access."""

    if name == "modflow":
        return import_module("myflopy.modflow")
    if name == "project":
        return import_module("myflopy.project")

    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = _EXPORTS[name]
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__():
    """Advertise the lazily-exported names for tab-completion."""

    return sorted(list(globals().keys()) + __all__)
