from __future__ import annotations

import os
import shutil
from pathlib import Path
from uuid import uuid4

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

from myflopy.modflow.mf6.canonical import CANONICAL_MODEL_CONTRACT
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    build_canonical_model,
)


def _pytest_temp_root() -> Path:
    configured = os.environ.get("SIMPLE_MODFLOW_PYTEST_TMP_ROOT")
    if configured:
        return Path(configured)
    return _REPO_ROOT / ".pytest-work" / "custom_tmp"


@pytest.fixture
def tmp_path():
    """Repo-local replacement for pytest's builtin tmp_path fixture.

    On this Windows setup, pytest's temp-path cleanup can intermittently leave
    behind folders with broken permissions. This fixture uses a normal temp
    directory we control and ignores cleanup failures instead of failing the
    whole test session.
    """

    root = _pytest_temp_root()
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"case_{uuid4().hex[:10]}"
    path.mkdir(parents=True, exist_ok=False)
    yield path
    try:
        shutil.rmtree(path)
    except OSError:
        pass


@pytest.fixture
def canonical_config():
    """Return the scaled validation profile of the authoritative model."""

    return CanonicalModelConfig.validation()


@pytest.fixture
def canonical_model(tmp_path, canonical_config):
    """Build and contract-check the canonical model for integration tests."""

    model = build_canonical_model(tmp_path / "canonical", config=canonical_config)
    CANONICAL_MODEL_CONTRACT.validate(model)
    return model


@pytest.fixture
def canonical_run(canonical_model):
    """Run and return the contract-checked canonical integration model."""

    success, report = canonical_model.run_simulation()
    assert success, "\n".join(report[-30:])
    CANONICAL_MODEL_CONTRACT.validate(canonical_model)
    return canonical_model


# --- Slow-test marking -------------------------------------------------------
# A minority of tests shell out to MODFLOW 6 / MODPATH binaries, run PEST
# builds, or exercise the parallel-model workflow. They dominate the suite's
# wall time (see ``pytest --durations``); ``test_parallel_model`` alone is ~70%.
# They are auto-marked ``slow`` here so the fast inner loop is:
#
#     pytest -m "not slow"
#
# A bare ``pytest`` still runs everything (CI / pre-push), so coverage is
# unchanged. Update the lists below if heavy tests are renamed or added.

# Whole modules whose every test runs a PEST build (all >1s).
_SLOW_MODULES = {
    "test_synthetic_pest_demo",
    "test_gold_standard_pest_demo",
}

# Individual heavy tests living in otherwise-fast modules. Any test using the
# ``canonical_run`` fixture is also marked automatically (see below), so this
# list only needs the heavy tests that build/run a model *without* that fixture.
_SLOW_TESTS = {
    # test_parallel_model.py -- partition/split runs (the rest are <1s units)
    "test_canonical_model_prepares_contiguous_partitions_across_representative_part_counts",
    # test_master_visualization_prt_example.py -- model build/run + contract
    # (the remaining tests in that module are sub-second doc checks)
    "test_master_example_full_profile_boundary_heads_stay_above_cell_bottoms",
    "test_master_example_validation_profile_runs_complex_package_topology",
    "test_canonical_fixture_enforces_complete_model_contract",
    # test_mf6_refactor_smoke.py -- end-to-end models that run MF6
    "test_refined_end_to_end_model_can_run_with_preferred_builder_api",
    "test_code_defined_geometries_can_drive_refinement_regions_and_packages",
    "test_feature_rich_small_model_workflow_can_run",
    "test_small_vector_sfr_lak_mvr_workflow_can_run",
    # test_mf6_prt.py
    "test_canonical_prt_run_stops_at_flow_end_and_writes_pathlines",
    "test_real_gwf_to_prt_run_and_shared_scene",
    # test_interactive_plotting.py -- large multilayer external-frame rendering
    "test_large_multilayer_results_use_external_frames_and_mask_dry_values",
    "test_large_multilayer_shared_color_scale_ignores_dry_and_nan_values",
    # test_mesh_optimization.py
    "test_region_touching_domain_boundary_builds",
    "test_optimized_triangle_grid_builds_voronoi_safe_mesh_and_runs_mf6",
    # test_mf6_pest.py -- forward run / reopen a completed run
    "test_pest_forward_run_applies_k_and_drn_parameters_end_to_end",
    "test_pest_run_results_reopen_completed_artifact_and_compare_heads",
    "test_native_pstfrom_parameterize_build_and_forward_run_end_to_end",
    "test_run_ies_end_to_end_and_assess_with_ies_results",
    # test_workspace.py -- executes MF6
    "test_run_builds_writes_executes_and_reopens_with_flopy_310",
    "test_project_multi_model_run_uses_per_model_subdirs",
    # test_layers.py -- flopy VTK -> pyvista export (heavy)
    "test_vtk_3d_writes_html_and_returns_iframe",
    "test_vtk_3d_shows_a_subset_of_layers",
    "test_views_runs_all_four",
}


def pytest_collection_modifyitems(config, items):
    """Auto-mark binary-running / parallel tests as ``slow``."""

    slow = pytest.mark.slow
    for item in items:
        if (
            item.path.stem in _SLOW_MODULES
            or item.originalname in _SLOW_TESTS
            or "canonical_run" in item.fixturenames
        ):
            item.add_marker(slow)
