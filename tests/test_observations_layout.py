"""The observations/ package layout (implementation plan 4.1).

Both import paths must resolve to the same objects: the historical module
path (``myflopy.modflow.mf6.observations``) and the split submodules.

The pest package used to import ``_normalize_row_labels`` from here to rebuild
pyEMU's observation names by SPELLING. It matches a target to the row PEST reads
by number and position now (ledger 340), so that private cross-import is gone,
and the test below holds it gone.
"""

from __future__ import annotations

import myflopy.modflow.mf6.observations as obs_pkg

EXPECTED_HOMES = {
    "HeadTargets": "heads",
    "ConcTargets": "conc",
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
    "RivFlowTargets": "riv",
    "BoundRivFlowTargets": "riv",
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


def test_no_copy_of_pyemus_row_label_rule_survives():
    """Matching a target by the spelling of its pst name dropped every target
    whose key pyEMU spelled differently (``per:1`` against ``per:1.0``)."""

    import myflopy.modflow.mf6.observations._shared as shared

    assert not hasattr(obs_pkg, "_normalize_row_labels")
    assert not hasattr(shared, "_normalize_row_labels")
