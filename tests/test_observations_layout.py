"""The observations/ package layout (implementation plan 4.1).

Both import paths must resolve to the same objects: the historical module
path (``myflopy.modflow.mf6.observations``) and the split submodules. The
pest package's private cross-import (``_normalize_row_labels``) is part of
the contract.
"""

from __future__ import annotations

import myflopy.modflow.mf6.observations as obs_pkg

EXPECTED_HOMES = {
    "HeadTargets": "heads",
    "BoundHeadTargets": "heads",
    "BoundHeadTargetPlots": "heads",
    "LakeStageTargets": "lake",
    "BoundLakeStageTargets": "lake",
    "SfrStageTargets": "sfr",
    "SfrFlowTargets": "sfr",
    "BoundSfrStageTargets": "sfr",
    "BoundSfrFlowTargets": "sfr",
    "DrnFlowTargets": "drn",
    "BoundDrnFlowTargets": "drn",
    "BoundNamedSeriesTargetPlots": "plots",
    "TargetRegistry": "registry",
}


def test_package_root_reexports_everything():
    assert sorted(obs_pkg.__all__) == sorted(EXPECTED_HOMES)
    for name in EXPECTED_HOMES:
        assert getattr(obs_pkg, name) is not None


def test_submodule_and_root_paths_are_the_same_objects():
    import importlib

    for name, stem in EXPECTED_HOMES.items():
        submodule = importlib.import_module(f"myflopy.modflow.mf6.observations.{stem}")
        assert getattr(submodule, name) is getattr(obs_pkg, name), name
        # __module__ points at the real home — pickles written before the
        # split still resolve through the package-root re-export.
        assert getattr(obs_pkg, name).__module__ == submodule.__name__, name


def test_pest_private_cross_import_stays_alive():
    from myflopy.modflow.mf6.observations import _normalize_row_labels
    from myflopy.modflow.mf6.observations._shared import (
        _normalize_row_labels as shared_fn,
    )

    assert _normalize_row_labels is shared_fn
