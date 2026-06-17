"""Run the gold-standard synthetic myflopy + PEST workflow."""

from __future__ import annotations

import argparse
from pathlib import Path

from myflopy.modflow.mf6.pest.gold_standard_demo import (
    GoldStandardPestDemoConfig,
    build_and_optionally_run_gold_standard_demo,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a compact but feature-rich synthetic MF6 model with UZF, lakes, "
            "a diverted stream network, drain seepage zones, canonical targets, and "
            "a PEST workspace using the currently implemented calibration parameters."
        )
    )
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=Path(__file__).resolve().parent / "artifacts",
        help="Directory where timestamped gold-standard demo workspaces will be created.",
    )
    parser.add_argument(
        "--run-family",
        default="gold_standard_pest_demo",
        help="Prefix used for the timestamped workspace directory.",
    )
    parser.add_argument("--nx", type=int, default=18, help="Number of columns in the DISV grid.")
    parser.add_argument("--ny", type=int, default=16, help="Number of rows in the DISV grid.")
    parser.add_argument("--cell-size", type=float, default=120.0, help="Square-cell width/height in model units.")
    parser.add_argument(
        "--pp-spacing",
        type=float,
        default=360.0,
        help="Pilot-point spacing used for K multipliers in map space.",
    )
    parser.add_argument(
        "--noptmax",
        type=int,
        default=10,
        help="PEST++ noptmax value. Defaults to 10 for a real optimization run.",
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
    config = GoldStandardPestDemoConfig(
        artifact_root=args.artifact_root,
        run_family=args.run_family,
        nx=args.nx,
        ny=args.ny,
        cell_size=args.cell_size,
        pp_spacing=args.pp_spacing,
    )

    run = build_and_optionally_run_gold_standard_demo(
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
    print(f"Control file: {run.pest_workspace / run.control_file}")
    print("Implemented calibration parameter families: k_pilotpoints, drn_elev, drn_cond")
    print("Saved target families: heads, lake_stage, sfr_stage, sfr_flow, drn_flow")
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
