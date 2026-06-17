from __future__ import annotations

import argparse
import json
import pickle
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from myflopy.modflow.mp3du.particles import ParticleTrackingInput


@dataclass(frozen=True)
class CumberlandForwardPaths:
    model_pickle: Path = Path(r"C:\Users\lukem\mf6\cum9a_algo4\cum9a_algo4.model")
    edge_particles: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\particles\edge_particles.shp")
    deep_lake_particles: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\particles\deep_lake_particles.gpkg")
    backward_particles: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\particles\cumb_backward_particles.gpkg")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare and run a Cumberland forward MP3DU particle-tracking job using a "
            "geometry-normalized particle release file."
        )
    )
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=ROOT / "examples" / "mf6" / "artifacts",
        help="Directory where timestamped MP3DU run workspaces will be created.",
    )
    parser.add_argument(
        "--run-family",
        default="cumberland_forward_mp3du",
        help="Prefix used for the timestamped artifact directory.",
    )
    parser.add_argument(
        "--model-pickle",
        type=Path,
        default=CumberlandForwardPaths.model_pickle,
        help="Pickled Cumberland model object to use as the flow-model source.",
    )
    parser.add_argument(
        "--particle-mode",
        choices=["explicit-cells", "interior", "existing"],
        default="explicit-cells",
        help=(
            "How to build the MP3DU start locations. 'explicit-cells' uses exact zero-based "
            "cell IDs, 'interior' generates new starts from eligible interior cells, and "
            "'existing' rebuilds starts from an existing particle file."
        ),
    )
    parser.add_argument(
        "--particle-file",
        type=Path,
        default=CumberlandForwardPaths.edge_particles,
        help="Source particle point file used when --particle-mode existing.",
    )
    parser.add_argument(
        "--cell-ids",
        default=None,
        help="Comma-separated zero-based cell IDs used when --particle-mode explicit-cells.",
    )
    parser.add_argument(
        "--cell-ids-file",
        type=Path,
        default=None,
        help=(
            "Text, CSV, or JSON file containing zero-based cell IDs used when "
            "--particle-mode explicit-cells."
        ),
    )
    parser.add_argument(
        "--interior-edge-rings",
        type=int,
        default=5,
        help="Number of outer active-cell edge rings to exclude for --particle-mode interior.",
    )
    parser.add_argument(
        "--interior-package-buffer-rings",
        type=int,
        default=2,
        help="Number of adjacency rings to exclude around selected boundary packages for interior starts.",
    )
    parser.add_argument(
        "--zloc",
        type=float,
        default=0.95,
        help="Normalized local-z release elevation to write into the prepared forward particle file.",
    )
    parser.add_argument(
        "--release-time",
        type=float,
        default=None,
        help="Override release time for every particle. Defaults to the source TimeRel values.",
    )
    parser.add_argument(
        "--simulation-end-time",
        type=float,
        default=None,
        help="Optional MP3DU simulation end time.",
    )
    parser.add_argument(
        "--max-particles",
        type=int,
        default=None,
        help="Optional cap on the number of prepared particles, applied after filtering.",
    )
    parser.add_argument(
        "--keep-inactive",
        action="store_true",
        help="Keep particles whose geometry maps to inactive cells.",
    )
    parser.add_argument(
        "--exclude-package",
        action="append",
        default=None,
        help=(
            "Boundary package to exclude from the prepared start set based on geometry-mapped "
            "start cells. Repeat for multiple packages."
        ),
    )
    parser.add_argument(
        "--iface",
        action="append",
        default=None,
        help="Override MP3DU IFACE as PACKAGE=VALUE. Repeat for multiple packages.",
    )
    parser.add_argument(
        "--porosity",
        type=float,
        default=0.2,
        help="Layer porosity used in the MP3DU path file.",
    )
    parser.add_argument(
        "--flow-thread-count",
        type=int,
        default=10,
        help="MP3DU flow-model thread count.",
    )
    parser.add_argument(
        "--pathline-thread-count",
        type=int,
        default=4,
        help="MP3DU pathline thread count.",
    )
    parser.add_argument(
        "--write-only",
        action="store_true",
        help="Only prepare inputs and diagnostics. Do not execute MP3DU.",
    )
    parser.add_argument(
        "--skip-output-conversion",
        action="store_true",
        help="Skip writeP3DOutput conversion after the MP3DU run.",
    )
    return parser.parse_args()


def build_workspace_root(artifact_root: Path, run_family: str) -> Path:
    artifact_root.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    workspace_root = artifact_root / f"{run_family}_{timestamp}"
    suffix = 1
    while workspace_root.exists():
        workspace_root = artifact_root / f"{run_family}_{timestamp}_{suffix:02d}"
        suffix += 1
    workspace_root.mkdir(parents=True)
    (artifact_root / f"{run_family}_latest.txt").write_text(str(workspace_root), encoding="utf-8")
    return workspace_root


def load_model(model_pickle: Path):
    with open(model_pickle, "rb") as handle:
        model = pickle.load(handle)
    idomain_path = getattr(model, "idomain_path", None)
    if isinstance(idomain_path, str):
        model.idomain_path = Path(idomain_path)
    return model


def _flatten_first_cell(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return _flatten_first_cell(value.tolist())
    if isinstance(value, list):
        return _flatten_first_cell(value[0]) if value else None
    if isinstance(value, tuple):
        if not value:
            return None
        if len(value) == 2 and all(isinstance(item, (int, np.integer)) for item in value):
            return int(value[1])
        return _flatten_first_cell(value[0])
    if isinstance(value, np.generic):
        value = value.item()
    return int(value)


def _default_iface_overrides() -> dict[str, int]:
    return {
        "CHD": 2,
        "DRN": 7,
        "EVT": 6,
        "GHB": 2,
        "RCH": 6,
        "RIV": 6,
        "SFR": 6,
        "WEL": 0,
    }


def parse_iface_overrides(values: list[str] | None) -> dict[str, int]:
    overrides = _default_iface_overrides()
    if not values:
        return overrides
    for entry in values:
        if "=" not in entry:
            raise ValueError(f"Invalid --iface value '{entry}'. Expected PACKAGE=VALUE.")
        package, raw_value = entry.split("=", 1)
        overrides[package.strip().upper()] = int(raw_value.strip())
    return overrides


def summarize_overlap(start_cells: list[int], boundary_sets: dict[str, set[int]]) -> dict[str, dict[str, float]]:
    total = len(start_cells)
    summary: dict[str, dict[str, float]] = {}
    for package_name, package_cells in sorted(boundary_sets.items()):
        count = sum(cell in package_cells for cell in start_cells)
        summary[package_name] = {
            "count": count,
            "pct": round((count / total) * 100.0, 2) if total else 0.0,
        }
    return summary


def parse_zero_based_cell_ids(*, cell_ids: str | None, cell_ids_file: Path | None) -> list[int]:
    values: list[int] = []
    if cell_ids:
        raw = cell_ids.replace("\n", ",").replace("\r", ",")
        values.extend(part.strip() for part in raw.split(",") if part.strip())

    if cell_ids_file is not None:
        text = cell_ids_file.read_text(encoding="utf-8")
        stripped = text.strip()
        if stripped.startswith("["):
            payload = json.loads(stripped)
            values.extend(str(item).strip() for item in payload if str(item).strip())
        else:
            normalized = stripped.replace("\n", ",").replace("\r", ",").replace(";", ",")
            values.extend(part.strip() for part in normalized.split(",") if part.strip())

    parsed = [int(value) for value in values]
    if not parsed:
        raise ValueError("No cell IDs were provided for --particle-mode explicit-cells.")
    return parsed


def sanitize_boundary_sets(model, boundary_sets: dict[str, set[int]]) -> dict[str, set[int]]:
    ncpl = int(model.vor.ncpl)
    return {
        package_name: {int(cell) for cell in cells if 0 <= int(cell) < ncpl}
        for package_name, cells in boundary_sets.items()
    }


def expand_cell_rings(adjacency: list[list[int]], seed_cells: set[int], rings: int) -> set[int]:
    seen = {int(cell) for cell in seed_cells if 0 <= int(cell) < len(adjacency)}
    frontier = set(seen)
    for _ in range(max(int(rings), 0)):
        next_frontier: set[int] = set()
        for cell in frontier:
            next_frontier.update(int(neighbor) for neighbor in adjacency[cell] if 0 <= int(neighbor) < len(adjacency))
        next_frontier -= seen
        seen |= next_frontier
        frontier = next_frontier
        if not frontier:
            break
    return seen


def prepare_existing_particle_file(
    *,
    model,
    particle_file: Path,
    output_dir: Path,
    zloc: float,
    release_time: float | None,
    max_particles: int | None,
    keep_inactive: bool,
    exclude_packages: list[str],
) -> tuple[Path, dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = gpd.read_file(particle_file).copy()
    geometry_cells = model.vor.get_vor_cells_as_series(frame)
    mapped_cells = [_flatten_first_cell(value) for value in geometry_cells.reset_index(drop=True).tolist()]

    frame["ModelCellZeroBased"] = mapped_cells
    frame["_source_order"] = np.arange(len(frame), dtype=int)
    frame["SourceFile"] = particle_file.name
    if "Node" in frame.columns:
        frame["SourceNode"] = frame["Node"]
    if "LocName" not in frame.columns:
        frame["LocName"] = [f"{particle_file.stem}_{idx + 1}" for idx in range(len(frame))]
    if release_time is not None or "TimeRel" not in frame.columns:
        frame["TimeRel"] = 0.0 if release_time is None else float(release_time)
    else:
        frame["TimeRel"] = frame["TimeRel"].astype(float)
    frame["ZLoc"] = float(zloc)

    unmapped_mask = frame["ModelCellZeroBased"].isna()
    kept = frame.loc[~unmapped_mask].copy()
    kept["ModelCellZeroBased"] = kept["ModelCellZeroBased"].astype(int)

    helper = ParticleTrackingInput(
        model=model,
        output_path=output_dir / "_prep_helper",
        porosities_by_layer=[0.2 for _ in range(model.gwf.modelgrid.nlay)],
        particle_shp=particle_file,
    )
    boundary_sets = sanitize_boundary_sets(model, helper.collect_boundary_cell_sets())
    inactive_cells = helper.collect_inactive_cells()
    excluded_cell_sets = {name.upper(): boundary_sets.get(name.upper(), set()) for name in exclude_packages}
    excluded_cells = set().union(*excluded_cell_sets.values()) if excluded_cell_sets else set()

    prep_summary: dict[str, Any] = {
        "source_particle_file": str(particle_file),
        "source_particle_count": int(len(frame)),
        "unmapped_particles": int(unmapped_mask.sum()),
        "source_overlap_by_geometry_cell": summarize_overlap(
            [int(cell) for cell in kept["ModelCellZeroBased"].tolist()],
            boundary_sets,
        ),
    }

    bad_mask = kept["ModelCellZeroBased"].isin(excluded_cells)
    if not keep_inactive:
        bad_mask = bad_mask | kept["ModelCellZeroBased"].isin(inactive_cells)

    good = kept.loc[~bad_mask].copy()
    bad = kept.loc[bad_mask].copy()

    prep_summary["reassigned_from_inactive"] = (
        0 if keep_inactive else int(bad["ModelCellZeroBased"].isin(inactive_cells).sum())
    )
    prep_summary["reassigned_from_package"] = {
        package_name: int(bad["ModelCellZeroBased"].isin(package_cells).sum())
        for package_name, package_cells in excluded_cell_sets.items()
    }

    relocated = gpd.GeoDataFrame(columns=kept.columns, geometry="geometry", crs=kept.crs)
    if len(bad):
        all_cells = set(range(model.vor.ncpl))
        eligible_cells = sorted(all_cells - excluded_cells - (set() if keep_inactive else inactive_cells))
        if not eligible_cells:
            raise RuntimeError("No eligible Cumberland start cells remain after applying the exclusion filters.")

        eligible = model.vor.gdf_vorPolys.loc[eligible_cells, ["geometry"]].copy()
        eligible["geometry"] = eligible.geometry.centroid
        eligible = eligible.reset_index().rename(columns={"index": "AssignedCellZeroBased"})

        bad_points = bad.drop(columns=["ModelCellZeroBased"]).copy().reset_index(drop=True)
        relocated = gpd.sjoin_nearest(
            bad_points,
            eligible[["AssignedCellZeroBased", "geometry"]],
            how="left",
            distance_col="RelocateDistance",
        )
        relocated = relocated.rename(columns={"AssignedCellZeroBased": "ModelCellZeroBased"})
        relocated["ModelCellZeroBased"] = relocated["ModelCellZeroBased"].astype(int)
        relocated["Relocated"] = True

    good["RelocateDistance"] = 0.0
    good["Relocated"] = False
    assigned = gpd.GeoDataFrame(
        pd.concat([good, relocated], ignore_index=True),
        geometry="geometry",
        crs=kept.crs,
    )
    assigned = assigned.sort_values("_source_order").reset_index(drop=True)
    centroid_lookup = model.vor.gdf_vorPolys.geometry.centroid
    assigned["geometry"] = assigned["ModelCellZeroBased"].map(centroid_lookup)

    if max_particles is not None:
        assigned = assigned.iloc[: int(max_particles)].copy()

    assigned["P3D_CellID"] = assigned["ModelCellZeroBased"] + 1
    prepared_columns = [
        "P3D_CellID",
        "TimeRel",
        "ZLoc",
        "LocName",
        "SourceFile",
        "ModelCellZeroBased",
        "Relocated",
        "RelocateDistance",
    ]
    if "SourceNode" in assigned.columns:
        prepared_columns.append("SourceNode")
    prepared = assigned.loc[:, prepared_columns + ["geometry"]].copy()
    prepared = prepared.rename(
        columns={
            "SourceFile": "SrcFile",
            "ModelCellZeroBased": "Cell0",
            "Relocated": "Moved",
            "RelocateDistance": "MoveDist",
            "SourceNode": "SrcNode",
        }
    )

    prepared_path = output_dir / f"{particle_file.stem}_forward_prepared.shp"
    for sidecar in output_dir.glob(f"{particle_file.stem}_forward_prepared.*"):
        sidecar.unlink()
    prepared.to_file(prepared_path, driver="ESRI Shapefile")

    prep_summary.update(
        {
            "prepared_particle_file": str(prepared_path),
            "prepared_particle_count": int(len(prepared)),
            "relocated_particles": int(prepared["Moved"].sum()),
            "max_relocation_distance": float(prepared["MoveDist"].max()) if len(prepared) else 0.0,
            "mean_relocation_distance": float(prepared["MoveDist"].mean()) if len(prepared) else 0.0,
            "prepared_overlap_by_geometry_cell": summarize_overlap(
                prepared["Cell0"].astype(int).tolist(),
                boundary_sets,
            ),
            "zloc": float(zloc),
            "release_time_override": None if release_time is None else float(release_time),
            "excluded_packages": [name.upper() for name in exclude_packages],
        }
    )
    return prepared_path, prep_summary


def prepare_interior_particle_file(
    *,
    model,
    output_dir: Path,
    zloc: float,
    release_time: float | None,
    max_particles: int | None,
    exclude_packages: list[str],
    edge_rings: int,
    package_buffer_rings: int,
) -> tuple[Path, dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)

    helper = ParticleTrackingInput(
        model=model,
        output_path=output_dir / "_prep_helper",
        porosities_by_layer=[0.2 for _ in range(model.gwf.modelgrid.nlay)],
        particle_shp=CumberlandForwardPaths.edge_particles,
    )
    boundary_sets = sanitize_boundary_sets(model, helper.collect_boundary_cell_sets())
    inactive_cells = helper.collect_inactive_cells()
    active_cells = set(range(model.vor.ncpl)) - inactive_cells
    edge_cells = set(model.vor.get_grid_edge(idomain=sorted(inactive_cells)))
    near_edge = expand_cell_rings(model.vor.adjacent_cells_idx, edge_cells, max(int(edge_rings) - 1, 0))

    excluded_package_sets = {name.upper(): boundary_sets.get(name.upper(), set()) for name in exclude_packages}
    excluded_package_cells = set().union(*excluded_package_sets.values()) if excluded_package_sets else set()
    near_package_cells = expand_cell_rings(
        model.vor.adjacent_cells_idx,
        excluded_package_cells,
        int(package_buffer_rings),
    )

    eligible_cells = sorted(active_cells - near_edge - near_package_cells)
    if not eligible_cells:
        raise RuntimeError("No eligible interior Cumberland cells remain after the requested filters.")

    if max_particles is not None and len(eligible_cells) > int(max_particles):
        sample_idx = np.linspace(0, len(eligible_cells) - 1, int(max_particles), dtype=int)
        eligible_cells = [eligible_cells[idx] for idx in sample_idx.tolist()]

    centroids = model.vor.gdf_vorPolys.geometry.centroid
    release_value = 0.0 if release_time is None else float(release_time)
    prepared = gpd.GeoDataFrame(
        {
            "P3D_CellID": [cell + 1 for cell in eligible_cells],
            "TimeRel": [release_value for _ in eligible_cells],
            "ZLoc": [float(zloc) for _ in eligible_cells],
            "LocName": [f"interior_{cell}" for cell in eligible_cells],
            "SrcFile": ["generated_interior" for _ in eligible_cells],
            "Cell0": eligible_cells,
            "Moved": [False for _ in eligible_cells],
            "MoveDist": [0.0 for _ in eligible_cells],
        },
        geometry=[centroids.iloc[cell] for cell in eligible_cells],
        crs=model.vor.gdf_vorPolys.crs,
    )

    prepared_path = output_dir / "cumberland_interior_forward_prepared.shp"
    for sidecar in output_dir.glob("cumberland_interior_forward_prepared.*"):
        sidecar.unlink()
    prepared.to_file(prepared_path, driver="ESRI Shapefile")

    prep_summary = {
        "source_particle_file": None,
        "source_particle_count": 0,
        "unmapped_particles": 0,
        "generated_mode": "interior",
        "active_cell_count": int(len(active_cells)),
        "edge_cell_count": int(len(edge_cells)),
        "excluded_edge_rings": int(edge_rings),
        "excluded_package_buffer_rings": int(package_buffer_rings),
        "excluded_packages": [name.upper() for name in exclude_packages],
        "excluded_package_seed_counts": {
            package_name: int(len(package_cells))
            for package_name, package_cells in excluded_package_sets.items()
        },
        "eligible_cell_count_before_sampling": int(len(active_cells - near_edge - near_package_cells)),
        "prepared_particle_file": str(prepared_path),
        "prepared_particle_count": int(len(prepared)),
        "relocated_particles": 0,
        "max_relocation_distance": 0.0,
        "mean_relocation_distance": 0.0,
        "prepared_overlap_by_geometry_cell": summarize_overlap(eligible_cells, boundary_sets),
        "zloc": float(zloc),
        "release_time_override": None if release_time is None else float(release_time),
    }
    return prepared_path, prep_summary


def prepare_explicit_cell_particle_file(
    *,
    model,
    output_dir: Path,
    zloc: float,
    release_time: float | None,
    cell_ids: list[int],
) -> tuple[Path, dict[str, Any]]:
    output_dir.mkdir(parents=True, exist_ok=True)

    helper = ParticleTrackingInput(
        model=model,
        output_path=output_dir / "_prep_helper",
        porosities_by_layer=[0.2 for _ in range(model.gwf.modelgrid.nlay)],
        particle_shp=CumberlandForwardPaths.edge_particles,
    )
    boundary_sets = sanitize_boundary_sets(model, helper.collect_boundary_cell_sets())
    inactive_cells = helper.collect_inactive_cells()

    invalid_cells = sorted({int(cell) for cell in cell_ids if int(cell) < 0 or int(cell) >= model.vor.ncpl})
    if invalid_cells:
        raise ValueError(f"Explicit cell list includes out-of-range IDs: {invalid_cells[:10]}")

    explicit_cells = [int(cell) for cell in cell_ids]
    centroids = model.vor.gdf_vorPolys.geometry.centroid
    release_value = 0.0 if release_time is None else float(release_time)
    prepared = gpd.GeoDataFrame(
        {
            "P3D_CellID": [cell + 1 for cell in explicit_cells],
            "TimeRel": [release_value for _ in explicit_cells],
            "ZLoc": [float(zloc) for _ in explicit_cells],
            "LocName": [f"cell_{cell}" for cell in explicit_cells],
            "SrcFile": ["explicit_cells" for _ in explicit_cells],
            "Cell0": explicit_cells,
            "Moved": [False for _ in explicit_cells],
            "MoveDist": [0.0 for _ in explicit_cells],
        },
        geometry=[centroids.iloc[cell] for cell in explicit_cells],
        crs=model.vor.gdf_vorPolys.crs,
    )

    prepared_path = output_dir / "cumberland_explicit_forward_prepared.shp"
    for sidecar in output_dir.glob("cumberland_explicit_forward_prepared.*"):
        sidecar.unlink()
    prepared.to_file(prepared_path, driver="ESRI Shapefile")

    prep_summary = {
        "source_particle_file": None,
        "source_particle_count": int(len(explicit_cells)),
        "unmapped_particles": 0,
        "generated_mode": "explicit-cells",
        "inactive_particle_count": int(sum(cell in inactive_cells for cell in explicit_cells)),
        "prepared_particle_file": str(prepared_path),
        "prepared_particle_count": int(len(prepared)),
        "relocated_particles": 0,
        "max_relocation_distance": 0.0,
        "mean_relocation_distance": 0.0,
        "prepared_overlap_by_geometry_cell": summarize_overlap(explicit_cells, boundary_sets),
        "zloc": float(zloc),
        "release_time_override": None if release_time is None else float(release_time),
        "cell_ids_zero_based": explicit_cells,
    }
    return prepared_path, prep_summary


def main():
    args = parse_args()
    default_excludes = ["DRN", "LAK", "SFR"] if args.particle_mode == "interior" else ["DRN", "LAK"]
    exclude_packages = [name.upper() for name in (args.exclude_package or default_excludes)]
    iface_overrides = parse_iface_overrides(args.iface)

    workspace_root = build_workspace_root(args.artifact_root, args.run_family)
    mp3du_workspace = workspace_root / "mp3du"
    particle_workspace = workspace_root / "particles"
    run_info_path = workspace_root / "run_info.json"

    print(f"Workspace root: {workspace_root}")
    print(f"Model pickle: {args.model_pickle}")
    if args.particle_mode == "existing":
        print(f"Source particle file: {args.particle_file}")
    else:
        print("Source particle file: generated interior starts")

    model = load_model(args.model_pickle)
    if args.particle_mode == "existing":
        prepared_particle_path, preparation_summary = prepare_existing_particle_file(
            model=model,
            particle_file=args.particle_file,
            output_dir=particle_workspace,
            zloc=args.zloc,
            release_time=args.release_time,
            max_particles=args.max_particles,
            keep_inactive=bool(args.keep_inactive),
            exclude_packages=exclude_packages,
        )
    else:
        prepared_particle_path, preparation_summary = prepare_interior_particle_file(
            model=model,
            output_dir=particle_workspace,
            zloc=args.zloc,
            release_time=args.release_time,
            max_particles=args.max_particles,
            exclude_packages=exclude_packages,
            edge_rings=args.interior_edge_rings,
            package_buffer_rings=args.interior_package_buffer_rings,
        )

    porosities = [float(args.porosity) for _ in range(model.gwf.modelgrid.nlay)]
    tracker = ParticleTrackingInput(
        model=model,
        output_path=mp3du_workspace,
        porosities_by_layer=porosities,
        particle_shp=prepared_particle_path,
        particle_field_map={
            "CELLID_ATTR": "P3D_CellID",
            "TIME_ATTR": "TimeRel",
            "ZLOC_ATTR": "ZLoc",
            "ADDTL_ATTR": ["LocName", "SrcFile", "Cell0"],
        },
        cellid_index_base=1,
        iface_overrides=iface_overrides,
        direction="FORWARD",
        simulation_end_time=args.simulation_end_time,
        flow_thread_count=args.flow_thread_count,
        pathline_thread_count=args.pathline_thread_count,
    )
    result = tracker.run(
        execute=not args.write_only,
        convert_output=(not args.write_only) and (not args.skip_output_conversion),
        write_diagnostics=True,
    )

    run_info = {
        "run_family": args.run_family,
        "workspace_root": str(workspace_root),
        "mp3du_workspace": str(mp3du_workspace),
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "model_pickle": str(args.model_pickle),
        "model_name": model.name,
        "particle_mode": args.particle_mode,
        "source_particle_file": None if args.particle_mode == "interior" else str(args.particle_file),
        "prepared_particle_file": str(prepared_particle_path),
        "iface_overrides": iface_overrides,
        "particle_preparation": preparation_summary,
        "mp3du_result": {
            "json_file": str(result["json_file"]),
            "path_file": str(result["path_file"]),
            "output_json": None if result["output_json"] is None else str(result["output_json"]),
            "diagnostics_file": None if result["diagnostics_file"] is None else str(result["diagnostics_file"]),
            "start_cell_diagnostics": result["start_cell_diagnostics"],
            "endpoint_summary": result["endpoint_summary"],
        },
    }
    run_info_path.write_text(json.dumps(run_info, indent=2), encoding="utf-8")

    print("\nPrepared forward particle summary")
    print(json.dumps(preparation_summary, indent=2))

    print("\nMP3DU start-cell diagnostics")
    print(json.dumps(result["start_cell_diagnostics"], indent=2))

    if result["endpoint_summary"] is not None:
        print("\nEndpoint termination summary")
        print(json.dumps(result["endpoint_summary"], indent=2))
    else:
        print("\nEndpoint termination summary not available yet.")

    print(f"\nRun info written to: {run_info_path}")


if __name__ == "__main__":
    main()
