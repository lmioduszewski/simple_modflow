from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from simple_modflow.modflow.mf6.pest.results import open_pest_run  # noqa: E402
from simple_modflow.modflow.mf6.pest.synthetic_demo import (  # noqa: E402
    SyntheticPestDemoConfig,
    build_and_optionally_run_synthetic_demo,
)


def _temp_artifact_root(name: str) -> Path:
    root = ROOT / ".pytest-work" / name
    if root.exists():
        shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    return root


def test_synthetic_demo_build_only_writes_reopenable_workspace():
    artifact_root = _temp_artifact_root("synthetic_demo_build")
    config = SyntheticPestDemoConfig(
        artifact_root=artifact_root,
        run_family="synthetic_demo_build",
        nx=8,
        ny=6,
        pp_spacing=250.0,
    )

    run = build_and_optionally_run_synthetic_demo(
        config,
        run_pestpp=False,
        noptmax=0,
    )

    assert run.workspace_root.exists()
    assert (run.pest_workspace / run.control_file).exists()
    assert (run.review_dir / "baseline_stats.csv").exists()
    assert (run.inputs_dir / "truth_targets.csv").exists()
    assert (run.pest_workspace / "simple_modflow_pest_metadata.json").exists()
    assert len(run.targets.to_long()) == len(run.observation_locations) * config.nper

    pest_run = open_pest_run(run.pest_workspace)
    loaded_targets = pest_run.load_head_targets()
    assert len(loaded_targets.to_long()) == len(run.targets.to_long())
    baseline_model = pest_run.load_baseline_model()
    stats = loaded_targets.stats(baseline_model)
    assert float(stats["rmse"].iloc[0]) > 0.0


def test_open_pest_run_prefers_completed_master_workspace_when_present():
    artifact_root = _temp_artifact_root("synthetic_demo_master_discovery")
    config = SyntheticPestDemoConfig(
        artifact_root=artifact_root,
        run_family="synthetic_demo_master_discovery",
        nx=8,
        ny=6,
        pp_spacing=250.0,
    )

    run = build_and_optionally_run_synthetic_demo(
        config,
        run_pestpp=False,
        noptmax=0,
    )

    master_dir = run.workspace_root / "pest_master"
    shutil.copytree(run.pest_workspace, master_dir)
    (master_dir / "synthetic_compact.par").write_text(
        "single point\nhk_pp_0000 1.0\n",
        encoding="utf-8",
    )

    pest_run = open_pest_run(run.workspace_root)
    assert pest_run.pest_workspace == master_dir
    assert pest_run.has_final_parameters
    assert pest_run.par_path is not None
    assert pest_run.par_path.name == "synthetic_compact.par"

    loaded_targets = pest_run.load_head_targets()
    assert len(loaded_targets.to_long()) == len(run.targets.to_long())


@pytest.mark.skipif(shutil.which("pestpp-glm") is None, reason="pestpp-glm is not available on PATH")
def test_synthetic_demo_short_pest_run_improves_fit_and_reopens():
    pytest.importorskip("pyemu")

    artifact_root = _temp_artifact_root("synthetic_demo_pest")
    config = SyntheticPestDemoConfig(
        artifact_root=artifact_root,
        run_family="synthetic_demo_pest",
        nx=8,
        ny=6,
        pp_spacing=250.0,
    )

    run = build_and_optionally_run_synthetic_demo(
        config,
        run_pestpp=True,
        noptmax=1,
        n_workers=1,
        pestpp_exe="pestpp-glm",
    )

    assert run.result_workspace is not None
    assert run.review_summary is not None
    assert run.review_summary["mae_calibrated"] < run.review_summary["mae_baseline"]
    assert run.review_summary["rmse_calibrated"] < run.review_summary["rmse_baseline"]

    review = open_pest_run(run.result_workspace).review()
    assert float(review.stats["mae_calibrated"].iloc[0]) < float(review.stats["mae_baseline"].iloc[0])
    assert float(review.k_geodata["k_ratio"].max()) > 1.2
    assert float(review.k_geodata["k_ratio"].min()) < 0.9

    drain_summary = pd.read_csv(run.review_dir / "drain_conductance_summary.csv")
    relative_change = (drain_summary["final_cond"] - drain_summary["start_cond"]).abs() / drain_summary["start_cond"]
    assert (relative_change > 0.10).any()
