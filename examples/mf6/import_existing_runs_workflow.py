from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from simple_modflow import import_run_archive, summarize_discovered_runs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Discover legacy MF6 runs, import them into a project catalog, and print simple summaries.",
    )
    parser.add_argument("--archive-root", type=Path, required=True, help="Root directory that contains existing MF6 run folders.")
    parser.add_argument("--catalog-root", type=Path, required=True, help="Catalog directory to create or update.")
    parser.add_argument("--project-name", default=None, help="Optional project catalog name.")
    parser.add_argument("--model-spec", default="legacy_archive", help="Model spec name to use for imported runs.")
    parser.add_argument("--compare", nargs=2, metavar=("RUN_A", "RUN_B"), help="Optional pair of run ids to compare.")
    return parser.parse_args()


def print_frame(title: str, frame: pd.DataFrame, rows: int = 10):
    print(f"\n{title}")
    print("=" * len(title))
    if frame.empty:
        print("<empty>")
        return
    print(frame.head(rows).to_string(index=False))


def main():
    args = parse_args()

    discovered = summarize_discovered_runs(args.archive_root).sort_values("run_id").reset_index(drop=True)
    print_frame("Discovered Runs", discovered, rows=20)

    catalog = import_run_archive(
        args.archive_root,
        args.catalog_root,
        project_name=args.project_name,
        model_spec=args.model_spec,
    )
    imported = catalog.run_summary().sort_values("run_id").reset_index(drop=True)
    print_frame("Imported Catalog Summary", imported, rows=20)

    if args.compare:
        run_a, run_b = args.compare
        package_diff = catalog.compare.compare_package_versions(run_a, run_b)
        head_diff = catalog.compare.compare_heads(run_a, run_b)
        head_stats = catalog.compare.compare_head_stats(run_a, run_b)
        print_frame(f"Package Differences: {run_a} vs {run_b}", package_diff)
        print_frame(f"Head Differences: {run_a} vs {run_b}", head_diff)
        print_frame(f"Head Stats: {run_a} vs {run_b}", head_stats)


if __name__ == "__main__":
    main()
