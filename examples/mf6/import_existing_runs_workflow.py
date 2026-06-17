from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from myflopy import Project


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Discover native MF6 workspaces through the clean Project API.",
    )
    parser.add_argument("--archive-root", type=Path, required=True, help="Root directory that contains existing MF6 run folders.")
    parser.add_argument("--project-root", type=Path, required=True, help="Project workspace to create or reopen.")
    parser.add_argument("--project-name", default=None, help="Optional project name.")
    return parser.parse_args()


def main():
    args = parse_args()
    project = Project(args.project_root, name=args.project_name)
    discovered = project.discover_native_runs(args.archive_root)

    print("\nDiscovered Native MF6 Runs")
    print("==========================")
    for run in discovered:
        print(f"{run.name}: {run.workspace}")


if __name__ == "__main__":
    main()
