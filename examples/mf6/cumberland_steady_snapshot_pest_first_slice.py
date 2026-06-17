from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import myflopy as mf
from myflopy.modflow.mf6.simulation import DisvGrid, TemporalDiscretization
from myflopy.modflow.mf6.simulation.packages import (
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)

from cumberland_transient_observed_pest_first_slice import (
    CumberlandPaths,
    build_single_observation_rows,
    build_recharge_dict,
    build_workspace_root,
    configure_solver,
    load_and_adjust_voronoi,
    prepare_pest_project,
    run_pestpp_live,
    run_pestpp_parallel,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a steady-state Cumberland snapshot calibration workflow using "
            "the current GIS inputs, one steady-state stress period, a representative "
            "observed head snapshot, and the first-slice PESTProject API."
        )
    )
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=ROOT / "examples" / "mf6" / "artifacts",
        help="Directory where timestamped model and PEST workspaces will be created.",
    )
    parser.add_argument(
        "--run-family",
        default="cumberland_steady_snapshot_pest_first_slice",
        help="Prefix used for the timestamped artifact directory.",
    )
    parser.add_argument(
        "--snapshot-per",
        type=int,
        default=None,
        help=(
            "Specific Cumberland observed period to use as the steady-state target snapshot. "
            "If omitted, the script picks the period with the broadest well coverage."
        ),
    )
    parser.add_argument(
        "--pp-spacing",
        type=float,
        default=2500.0,
        help="Pilot-point spacing for the K parameter field.",
    )
    parser.add_argument(
        "--fast-mode",
        action="store_true",
        help=(
            "Use a lighter first calibration setup: wider pilot-point spacing and "
            "K-only parameters."
        ),
    )
    parser.add_argument(
        "--k-only",
        action="store_true",
        help="Only expose K pilot points to PEST++. Skip drain-conductance parameters.",
    )
    parser.add_argument(
        "--mfr-factor",
        type=float,
        default=1.0,
        help="Fallback multiplier for MFR polygons that do not define their own factor field.",
    )
    parser.add_argument(
        "--steady-recharge",
        type=float,
        default=0.007,
        help="Uniform recharge applied during the single steady-state period.",
    )
    parser.add_argument(
        "--supplemental-weight",
        type=float,
        default=0.25,
        help="Observation weight assigned to one-time supplemental head targets.",
    )
    parser.add_argument(
        "--noptmax",
        type=int,
        default=0,
        help="PEST++ noptmax value. Use 0 for setup/check only.",
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
        help="Number of local PEST++ workers to launch when --run-pestpp is used.",
    )
    parser.add_argument(
        "--keep-workers",
        action="store_true",
        help="Keep parallel worker folders on disk for debugging instead of cleaning them up.",
    )
    parser.add_argument(
        "--pestpp-exe",
        default="pestpp-glm",
        help="PEST++ executable name or full path.",
    )
    return parser.parse_args()


def select_snapshot_row(
    obs_frame: pd.DataFrame,
    *,
    snapshot_per: int | None,
) -> pd.Series:
    obs_frame = obs_frame.sort_values("per").reset_index(drop=True)
    well_columns = [column for column in obs_frame.columns if column not in {"per", "Deep Lake"}]
    coverage = obs_frame[well_columns].notna().sum(axis=1)

    if snapshot_per is not None:
        selected = obs_frame.loc[obs_frame["per"].astype(int) == int(snapshot_per)].copy()
        if selected.empty:
            raise ValueError(f"Could not find snapshot period {snapshot_per} in calib_observations.xlsx.")
        return selected.iloc[0]

    best_idx = int(coverage.idxmax())
    return obs_frame.iloc[best_idx]


def build_snapshot_targets(
    paths: CumberlandPaths,
    obs_frame: pd.DataFrame,
    snapshot_row: pd.Series,
    *,
    supplemental_weight: float,
) -> tuple[mf.HeadTargets, pd.DataFrame]:
    locs = gpd.read_file(paths.calib_locs).copy()
    locs["name"] = locs["ExploName"].astype(str)
    locs["layer_num"] = 0
    locs["weight"] = 1.0
    allowed_names = set(locs["name"].tolist())

    rows: list[dict[str, float | int | str]] = []
    for name, value in snapshot_row.items():
        if name == "per" or name == "Deep Lake" or name not in allowed_names or pd.isna(value):
            continue
        rows.append({"per": 0, "name": str(name), "head": float(value)})
    selected_names = {row["name"] for row in rows}
    rows.extend(
        build_single_observation_rows(
            obs_frame,
            allowed_names=allowed_names,
            selected_names=selected_names,
            weight=float(supplemental_weight),
        )
    )

    supplemental_alias = {
        "MW-1 (HWA)": "MW1_HWA",
        "MW-2 (HWA)": "MW2_HWA",
        "MW-3 (HWA)": "MW3_HWA",
        "MW-4 (HWA)": "MW4_HWA",
        "B-1 (GD)": "B1_GD",
        "B-2 (GD)": "B2_GD",
    }
    supplemental_names = list(supplemental_alias.keys())

    merged = gpd.read_file(paths.merged_wells).copy()
    merged["name"] = merged["ExploName"].astype(str)
    merged = merged.loc[merged["name"].isin(supplemental_names), ["name", "geometry"]].copy()
    merged["name"] = merged["name"].map(supplemental_alias)
    merged["layer_num"] = 0
    merged["weight"] = float(supplemental_weight)

    boring = pd.read_csv(paths.boring_summary).copy()
    boring["name"] = boring["Boring"].astype(str)
    boring["gw_elev"] = pd.to_numeric(boring["GW Elev. (ft)"], errors="coerce")
    boring = boring.loc[boring["name"].isin(supplemental_names), ["name", "gw_elev"]].copy()
    overrides = {"B-1 (GD)": 836.0, "B-2 (GD)": 800.0}
    for name, value in overrides.items():
        if name in boring["name"].tolist():
            boring.loc[boring["name"] == name, "gw_elev"] = value
        else:
            boring = pd.concat([boring, pd.DataFrame([{"name": name, "gw_elev": value}])], ignore_index=True)
    boring["name"] = boring["name"].map(supplemental_alias)
    rows.extend(
        {"per": 0, "name": str(row.name), "head": float(row.gw_elev)}
        for row in boring.itertuples(index=False)
        if pd.notna(row.gw_elev)
    )

    target_values = pd.DataFrame(rows)
    active_names = sorted(target_values["name"].unique().tolist())
    locs = locs.loc[locs["name"].isin(active_names), ["name", "layer_num", "weight", "geometry"]].copy()
    merged = merged.loc[merged["name"].isin(active_names), ["name", "layer_num", "weight", "geometry"]].copy()
    combined_locs = pd.concat([locs, merged], ignore_index=True).drop_duplicates(subset=["name"], keep="first")
    if "weight_override" in target_values.columns:
        override_weights = (
            target_values.loc[target_values["weight_override"].notna(), ["name", "weight_override"]]
            .drop_duplicates(subset=["name"], keep="last")
            .rename(columns={"weight_override": "weight"})
        )
        combined_locs = combined_locs.merge(override_weights, on="name", how="left", suffixes=("", "_override"))
        mask = combined_locs["weight_override"].notna()
        combined_locs.loc[mask, "weight"] = combined_locs.loc[mask, "weight_override"]
        combined_locs = combined_locs.drop(columns=["weight_override"])
        target_values = target_values.drop(columns=["weight_override"])
    targets = mf.HeadTargets(
        locations=combined_locs,
        values=target_values,
        name_column="name",
        layer_column="layer_num",
        time_column="per",
        value_column="head",
    )
    return targets, target_values


def build_model(
    *,
    paths: CumberlandPaths,
    model_workspace: Path,
    model_name: str,
    mfr_factor: float,
    steady_recharge: float,
):
    vor, iheads, idomain, inactive_cells, top, botm = load_and_adjust_voronoi(paths)
    model = mf.SimulationBase(
        name=model_name,
        nper=1,
        vor=vor,
        mf_folder_path=model_workspace,
    )
    DisvGrid(vor=vor, model=model, nlay=1, top=top, bottom=botm, idomain=idomain)
    TemporalDiscretization(model=model, period_data=[[1.0, 15, 1.1]])
    model.per_dates = [pd.Timestamp("2020-01-01")]
    OutputControl(model=model)

    k_array = mf.KFromVector(
        model=model,
        vor=vor,
        shp_gpkg=paths.ks_v5,
        uid="name",
        crs=2926,
        idomain_path=paths.idomain_path,
    ).from_vector(nlay=1, defaults=[10.0])
    KFlow(model=model, k=k_array[0], k33_vert=(k_array[0] / 10.0), save_specific_discharge=False)

    InitialConditions(model=model, vor=vor, strt=[iheads])
    Storage(
        model=model,
        specific_yield=0.2,
        specific_storage=1e-4,
        sto_steady={0: True},
        sto_transient={},
    )

    rch_dict = build_recharge_dict(
        vor,
        idomain,
        [],
        mfr_path=paths.mfr_polys,
        default_mfr_factor=mfr_factor,
        steady_recharge=steady_recharge,
    )
    Recharge(model=model, vor=vor, rch_dict={0: rch_dict[0]})

    drn_fields = {
        "name": "name",
        "height_over_btm": "height",
        "conductance": "cond",
        "layer": "layer",
        "min_elev": "min_elev",
    }
    drn_spd = mf.DRNFromVector(
        model=model,
        vor=vor,
        shp_gpkg=paths.drn_path,
        uid="name",
        crs=2926,
        idomain_path=paths.idomain_path,
    ).from_vector(fields=drn_fields, edges_only=False, top_drain=True)
    Drains(model=model, stress_period_data=drn_spd)

    configure_solver(model)
    return model


def main():
    args = parse_args()
    if args.fast_mode:
        args.k_only = True
        if args.pp_spacing < 4000.0:
            args.pp_spacing = 4000.0

    paths = CumberlandPaths()
    workspace_root = build_workspace_root(args.artifact_root, args.run_family)
    model_workspace = workspace_root / "model"
    pest_workspace = workspace_root / "pest"
    run_info = {
        "run_family": args.run_family,
        "workspace_root": str(workspace_root),
        "model_workspace": str(model_workspace),
        "pest_workspace": str(pest_workspace),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    run_info_path = workspace_root / "run_info.json"
    run_info_path.write_text(json.dumps(run_info, indent=2), encoding="utf-8")

    print(f"Workspace root: {workspace_root}")
    print(f"Model workspace parent: {model_workspace}")
    print(f"PEST workspace: {pest_workspace}")

    obs_frame = pd.read_excel(paths.calib_obs)
    snapshot_row = select_snapshot_row(obs_frame, snapshot_per=args.snapshot_per)
    snapshot_per = int(snapshot_row["per"])
    forcing_frame = pd.read_excel(paths.forcing_table)
    forcing_row = forcing_frame.loc[forcing_frame["per"].astype(int) == snapshot_per].copy()
    snapshot_month = (
        pd.to_datetime(forcing_row["Month"].iloc[0]).date().isoformat()
        if not forcing_row.empty and "Month" in forcing_row.columns
        else None
    )
    print(f"Using steady-state snapshot period: {snapshot_per}")
    if snapshot_month is not None:
        print(f"Snapshot month: {snapshot_month}")

    targets, target_values = build_snapshot_targets(
        paths,
        obs_frame,
        snapshot_row,
        supplemental_weight=args.supplemental_weight,
    )
    target_values_path = workspace_root / "head_targets_values.csv"
    target_values.to_csv(target_values_path, index=False)
    target_locations_path = workspace_root / "head_target_locations.gpkg"
    if target_locations_path.exists():
        target_locations_path.unlink()
    targets.locations_gdf.to_file(target_locations_path, driver="GPKG")

    run_info.update(
        {
            "snapshot_per": snapshot_per,
            "snapshot_month": snapshot_month,
            "steady_recharge": float(args.steady_recharge),
            "supplemental_weight": float(args.supplemental_weight),
            "pp_spacing": float(args.pp_spacing),
            "k_only": bool(args.k_only),
            "fast_mode": bool(args.fast_mode),
            "head_target_locations": str(target_locations_path),
            "head_target_values": str(target_values_path),
        }
    )
    run_info_path.write_text(json.dumps(run_info, indent=2), encoding="utf-8")
    print(
        f"Built {len(target_values)} steady-state head targets across "
        f"{target_values['name'].nunique()} wells."
    )

    model = build_model(
        paths=paths,
        model_workspace=model_workspace,
        model_name="cpre_ss_obs",
        mfr_factor=args.mfr_factor,
        steady_recharge=args.steady_recharge,
    )
    success, _ = model.run_simulation()
    if not success:
        raise RuntimeError("Baseline steady-state Cumberland model did not solve successfully.")

    baseline_stats = targets.stats(model)
    print("\nBaseline residual stats")
    print(baseline_stats.to_string(index=False))

    start_datetime = snapshot_month or "2020-01-01"
    pest = prepare_pest_project(
        model=model,
        workspace_root=workspace_root,
        pest_workspace=pest_workspace,
        ks_v5=paths.ks_v5,
        drn_path=paths.drn_path,
        targets=targets,
        start_datetime=start_datetime,
        pp_spacing=args.pp_spacing,
        include_drn_cond=not args.k_only,
    )
    pst = pest.build_pst(f"{model.name}.pst")
    pst.control_data.noptmax = int(args.noptmax)
    pst_path = pest_workspace / f"{model.name}.pst"
    pst.write(pst_path)

    print("\nPEST workspace summary")
    print(f"npar_adj: {pst.npar_adj}")
    print(f"nobs: {pst.nobs}")
    print(f"Control file: {pst_path}")
    print(f"K-only mode: {args.k_only}")
    print(f"Pilot-point spacing: {args.pp_spacing}")

    result_workspace = workspace_root
    if args.run_pestpp:
        if int(args.n_workers) > 1:
            result_workspace = run_pestpp_parallel(
                workspace_root=workspace_root,
                pest_workspace=pest_workspace,
                control_file=pst_path.name,
                exe_name=args.pestpp_exe,
                n_workers=args.n_workers,
                keep_workers=args.keep_workers,
            )
        else:
            run_pestpp_live(pest_workspace, pst_path.name, args.pestpp_exe)
        pest_run = mf.open_pest_run(result_workspace)
        print("\nCompleted run summary")
        print(pest_run.summary().to_string(index=False))
        if pest_run.has_final_parameters:
            print("\nResidual comparison stats")
            print(pest_run.compare_head_target_stats(targets).to_string(index=False))

    print("\nNext steps")
    print(f"- Reopen results with: mf.open_pest_run(r'{result_workspace}')")
    print(f"- Latest pointer file: {args.artifact_root / f'{args.run_family}_latest.txt'}")


if __name__ == "__main__":
    main()
