from __future__ import annotations

from pathlib import Path
from typing import get_type_hints

import simple_modflow as mf
from simple_modflow.modflow.mf6.observations import (
    BoundDrnFlowTargets,
    BoundHeadTargets,
    BoundLakeStageTargets,
    BoundSfrFlowTargets,
    BoundSfrStageTargets,
    TargetRegistry,
)
from simple_modflow.modflow.mf6.simulation.base import SimulationBase


def test_package_declares_inline_typing_support():
    package_root = Path(mf.__file__).resolve().parent
    assert (package_root / "py.typed").exists()


def test_primary_model_namespaces_have_ide_visible_return_annotations():
    expected = {
        "targets": "'TargetRegistry'",
        "visualize": "'ModelVisualization'",
        "particle_tracking": "'ParticleTracking'",
        "parallel": "'ParallelModelWorkflow'",
        "outputs": "'ModelOutputs'",
        "packages": "'ModelPackages'",
    }
    for name, return_annotation in expected.items():
        assert getattr(SimulationBase, name).fget.__annotations__["return"] == return_annotation


def test_canonical_target_names_are_explicit_typed_properties():
    expected = {
        "heads": BoundHeadTargets,
        "lake_stage": BoundLakeStageTargets,
        "sfr_stage": BoundSfrStageTargets,
        "sfr_flow": BoundSfrFlowTargets,
        "drn_flow": BoundDrnFlowTargets,
    }
    for name, return_annotation in expected.items():
        prop = getattr(TargetRegistry, name)
        assert isinstance(prop, property)
        assert get_type_hints(prop.fget)["return"] is return_annotation


def test_explicit_target_properties_preserve_missing_attribute_semantics():
    registry = TargetRegistry(object())
    assert not hasattr(registry, "heads")
