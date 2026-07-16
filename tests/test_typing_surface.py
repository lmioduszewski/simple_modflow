from __future__ import annotations

from pathlib import Path
from typing import get_type_hints

import myflopy as mf
from myflopy.modflow.mf6.observations import (
    BoundDrnFlowTargets,
    BoundHeadTargets,
    BoundLakeStageTargets,
    BoundSfrFlowTargets,
    BoundSfrStageTargets,
    TargetRegistry,
)
from myflopy.modflow.mf6.simulation.base import SimulationBase


def test_package_declares_inline_typing_support():
    package_root = Path(mf.__file__).resolve().parent
    assert (package_root / "py.typed").exists()


def test_top_level_exports_signal_package_first_api():
    assert "lak" in mf.__all__
    assert "sfr" in mf.__all__
    assert "rch" in mf.__all__
    assert "LAKBuilder" not in mf.__all__
    assert "lak_spec" not in mf.__all__
    assert "build_ims" not in mf.__all__
    assert "Wells" not in mf.__all__

    assert "lak" in dir(mf)
    assert "LAKBuilder" not in dir(mf)
    assert "lak_spec" not in dir(mf)
    assert "LAKBuilder" in mf.__compatibility__
    assert "lak_spec" in mf.__compatibility__
    assert "Wells" in mf.__compatibility__


def test_second_tier_top_level_exports_remain_explicitly_importable():
    from myflopy import LAKBuilder, Wells, lak_spec

    assert LAKBuilder.__name__ == "LAKBuilder"
    assert Wells.__name__ == "Wells"
    assert callable(lak_spec)


def test_primary_model_namespaces_have_ide_visible_return_annotations():
    # Bare names (not quote-nested): the module uses `from __future__ import
    # annotations`, so every annotation is already a string and ruff UP037
    # strips the redundant inner quotes.
    expected = {
        "targets": "TargetRegistry",
        "visualize": "ModelVisualization",
        "particle_tracking": "ParticleTracking",
        "parallel": "ParallelModelWorkflow",
        "outputs": "ModelOutputs",
        "packages": "ModelPackages",
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
