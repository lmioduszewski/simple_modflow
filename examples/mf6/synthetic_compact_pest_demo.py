"""Run a compact synthetic PEST workflow with obvious K and drain changes."""

from __future__ import annotations

import argparse
from pathlib import Path

from simple_modflow.modflow.mf6.pest.synthetic_demo import (
    SyntheticPestDemoConfig,
    build_and_optionally_run_synthetic_demo,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a small synthetic transient model, generate pseudo-observations "
            "from a truth model, write a PEST workspace, and optionally run pestpp-glm."
        )
    )
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=Path(__file__).resolve().parent / "artifacts",
        help="Directory where timestamped synthetic demo workspaces will be created.",
    )
    parser.add_argument(
        "--run-family",
        default="synthetic_compact_pest_demo",
        help="Prefix used for the timestamped workspace directory.",
    )
    parser.add_argument("--nx", type=int, default=20, help="Number of columns in the compact DISV grid.")
    parser.add_argument("--ny", type=int, default=20, help="Number of rows in the compact DISV grid.")
    parser.add_argument("--cell-size", type=float, default=100.0, help="Square-cell width/height in model units.")
    parser.add_argument(
        "--pp-spacing",
        type=float,
        default=300.0,
        help="Pilot-point spacing used for K multipliers in map space.",
    )
    parser.add_argument(
        "--noptmax",
        type=int,
        default=10,
        help="PEST++ noptmax value. Use 0 for setup/check only. Defaults to 10 for a real optimization run.",
    )
    parser.add_argument(
        "--run-pestpp",
        action="store_true",
        help="Launch pestpp-glm after writing the control file.",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=1,
        help="Number of local workers to launch when --run-pestpp is used. Use 8 or more for parallel runs.",
    )
    parser.add_argument(
        "--keep-workers",
        action="store_true",
        help="Keep worker folders on disk instead of cleaning them up after a parallel run.",
    )
    parser.add_argument(
        "--pestpp-exe",
        default="pestpp-glm",
        help="PEST++ executable name or full path.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    config = SyntheticPestDemoConfig(
        artifact_root=args.artifact_root,
        run_family=args.run_family,
        nx=args.nx,
        ny=args.ny,
        cell_size=args.cell_size,
        pp_spacing=args.pp_spacing,
    )

    run = build_and_optionally_run_synthetic_demo(
        config,
        run_pestpp=args.run_pestpp,
        noptmax=args.noptmax,
        n_workers=args.n_workers,
        keep_workers=args.keep_workers,
        pestpp_exe=args.pestpp_exe,
    )

    print(f"Workspace root: {run.workspace_root}")
    print(f"Truth workspace: {run.truth_workspace}")
    print(f"Model workspace: {run.model_workspace}")
    print(f"PEST workspace: {run.pest_workspace}")
    print(f"Observation targets: {len(run.targets.to_long())}")
    print(f"Control file: {run.pest_workspace / run.control_file}")
    if run.result_workspace is not None:
        print(f"Result workspace: {run.result_workspace}")
    if run.review_summary:
        print("Review summary:")
        for key, value in run.review_summary.items():
            print(f"  {key}: {value}")
        print(f"Review artifacts: {run.review_dir}")
    else:
        print("Review artifacts were not written yet. Run with --run-pestpp to generate calibrated outputs.")


if __name__ == "__main__":
    main()
