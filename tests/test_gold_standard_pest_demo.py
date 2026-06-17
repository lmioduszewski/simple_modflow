from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from myflopy.modflow.mf6.pest.gold_standard_demo import (  # noqa: E402
    GoldStandardPestDemoConfig,
    build_and_optionally_run_gold_standard_demo,
)
from myflopy.modflow.mf6.pest.results import open_pest_run  # noqa: E402


def _temp_artifact_root(name: str) -> Path:
    root = ROOT / ".pytest-work" / name
    if root.exists():
        shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True, exist_ok=True)
    return root


def test_gold_standard_demo_build_only_writes_reopenable_workspace():
    artifact_root = _temp_artifact_root("gold_standard_demo_build")
    config = GoldStandardPestDemoConfig(
        artifact_root=artifact_root,
        run_family="gold_standard_demo_build",
        nx=8,
        ny=6,
        pp_spacing=360.0,
    )

    run = build_and_optionally_run_gold_standard_demo(
        config,
        run_pestpp=False,
        noptmax=0,
    )

    assert run.workspace_root.exists()
    assert (run.pest_workspace / run.control_file).exists()
    assert (run.pest_workspace / "myflopy_pest_metadata.json").exists()
    assert (run.model_workspace / config.model_name / "gold_sfr_stage.csv").exists()
    assert (run.model_workspace / config.model_name / f"{config.model_name}_sfr_budget.sfr").exists()

    pest_run = open_pest_run(run.pest_workspace)
    assert not pest_run.load_head_targets().to_long().empty
    assert not pest_run.load_lake_stage_targets().to_long().empty
    assert not pest_run.load_sfr_stage_targets().to_long().empty
    assert not pest_run.load_sfr_flow_targets().to_long().empty
    assert not pest_run.load_drn_flow_targets().to_long().empty

    mvr = pd.read_csv(run.model_workspace / config.model_name / f"{config.model_name}.mvr.csv")
    assert (mvr["TOTAL_IN"].abs() > 0.0).any()
    assert (mvr["LAK_OUT"].abs() > 0.0).any()

    lake_budget = pd.read_csv(run.model_workspace / config.model_name / f"{config.model_name}_lake_budget.csv")
    assert (lake_budget["FROM-MVR_IN"].abs() > 0.0).any()
    assert ((lake_budget["GWF_IN"].abs() > 0.0) | (lake_budget["GWF_OUT"].abs() > 0.0)).any()


def test_gold_standard_demo_forward_run_regenerates_named_series_outputs():
    artifact_root = _temp_artifact_root("gold_standard_demo_forward")
    config = GoldStandardPestDemoConfig(
        artifact_root=artifact_root,
        run_family="gold_standard_demo_forward",
        nx=8,
        ny=6,
        pp_spacing=360.0,
    )

    run = build_and_optionally_run_gold_standard_demo(
        config,
        run_pestpp=False,
        noptmax=0,
    )

    expected_outputs = [
        run.pest_workspace / "hds_simulated_heads.csv",
        run.pest_workspace / "lak_stage_simulated_lake_stage.csv",
        run.pest_workspace / "sfr_stage_simulated_sfr_stage.csv",
        run.pest_workspace / "sfr_flow_simulated_sfr_flow.csv",
        run.pest_workspace / "drn_flow_simulated_drn_flow.csv",
    ]
    for path in expected_outputs:
        if path.exists():
            path.unlink()

    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=run.pest_workspace,
        check=True,
        capture_output=True,
        text=True,
    )

    for path in expected_outputs:
        assert path.exists(), f"Expected regenerated output {path.name}"
        frame = pd.read_csv(path)
        assert not frame.empty
