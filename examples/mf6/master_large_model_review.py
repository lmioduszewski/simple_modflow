"""Run bounded review sections for the large visualization and PRT example."""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime
from pathlib import Path

import numpy as np

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import representative_cells


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=Path(__file__).resolve().parent / "artifacts",
    )
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d_%H%M%S"))
    parser.add_argument("--frame-stride", type=int, default=4)
    parser.add_argument("--max-frames", type=int, default=3)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--export-plotly", action="store_true")
    parser.add_argument("--run-prt", action="store_true")
    parser.add_argument("--export-pyvista", action="store_true")
    return parser.parse_args()


def timed(summary: dict, name: str, callback):
    start = time.perf_counter()
    result = callback()
    summary[f"{name}_seconds"] = round(time.perf_counter() - start, 3)
    return result


def main():
    args = parse_args()
    config = mf.CanonicalModelConfig()
    root = args.artifact_root / f"master_large_review_{args.run_id}"
    html_dir = root / "html"
    html_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "workspace": str(root),
        "ncpl": config.ncpl,
        "nlay": config.nlay,
        "total_cells": config.ncpl * config.nlay,
        "nper": config.nper,
    }

    model = timed(
        summary,
        "build",
        lambda: mf.build_canonical_model(root / "gwf", config=config, name="master_review"),
    )
    success, report = timed(summary, "gwf_run", model.run_simulation)
    if not success:
        raise RuntimeError("\n".join(report[-30:]))

    head_reader = model.gwf.output.head()
    kstpkpers = list(head_reader.get_kstpkper())
    final_heads = np.asarray(head_reader.get_data(kstpkper=kstpkpers[-1]), dtype=float)
    finite = final_heads[np.isfinite(final_heads) & (np.abs(final_heads) < 1.0e29)]
    summary.update(
        {
            "head_frames": len(kstpkpers),
            "head_min": float(finite.min()),
            "head_max": float(finite.max()),
            "head_mean": float(finite.mean()),
            "packages": sorted(model.gwf.package_names),
            "lake_stage_rows": len(model.packages.lak.results.stage.get()),
            "sfr_exchange_rows": len(model.packages.sfr.results.q.get()),
            "uzf_cells": len(model.outputs.uzf.ifno_to_cellid),
        }
    )

    style = mf.ModelMapStyle(show_contours=False, dpi=55, figsize=(6, 5))
    head_slider = timed(
        summary,
        "head_map_export",
        lambda: model.visualize.head_map_slider_html(
            html_dir / "head_map_review.html",
            layer=0,
            style=style,
            embed_frames=False,
            frame_stride=args.frame_stride,
            max_frames=args.max_frames,
            resume=args.resume,
            progress=True,
        ),
    )
    mosaic_slider = timed(
        summary,
        "mosaic_export",
        lambda: model.visualize.head_layer_mosaic_slider_html(
            html_dir / "head_mosaic_review.html",
            layers=list(range(config.nlay)),
            ncols=2,
            style=style,
            embed_frames=False,
            frame_stride=args.frame_stride,
            max_frames=args.max_frames,
            resume=args.resume,
            progress=True,
        ),
    )
    summary["exports"] = {
        "head_map": str(head_slider.path),
        "head_mosaic": str(mosaic_slider.path),
        "head_map_frames": head_slider.frame_count,
        "mosaic_frames": mosaic_slider.frame_count,
    }

    if args.export_plotly:
        plotly_path = html_dir / "plotly_head_map_review.html"
        plotly_map = timed(
            summary,
            "plotly_head_map_export",
            lambda: model.visualize.plotly_head_map_animation(
                output_path=plotly_path,
                layer=0,
                show_layer_elevs=False,
                frame_stride=args.frame_stride,
                max_frames=args.max_frames,
                include_plotlyjs="cdn",
            ),
        )
        summary["exports"].update(
            {
                "plotly_head_map": str(plotly_path),
                "plotly_head_map_frames": len(plotly_map.frames),
                "plotly_head_map_bytes": plotly_path.stat().st_size,
            }
        )

    if args.run_prt:
        releases = mf.PRTReleasePoints.from_cells(model, representative_cells(config)["releases"])
        project = model.particle_tracking.prt(
            workspace=root / "prt",
            release_points=releases,
            porosity=0.25,
            stop_at_weak_sink=False,
        )
        result = timed(summary, "prt_run", lambda: project.run(silent=True))
        summary["prt_rows"] = len(result.pathlines)
        summary["prt_terminal_rows"] = len(result.terminal_points)
        if args.export_pyvista:
            path = timed(
                summary,
                "pyvista_export",
                lambda: result.export_3d_html(html_dir / "prt_scene.html", vertical_exaggeration=5.0),
            )
            summary["prt_scene"] = str(path)

    summary_path = root / "review_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
