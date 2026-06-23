from __future__ import annotations

import argparse
import json
import pickle
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import myflopy as mf
from myflopy import read_shp_gpkg
from myflopy.modflow.mf6.simulation import DisvGrid, TemporalDiscretization
from myflopy.modflow.mf6.simulation.packages import (
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)


@dataclass(frozen=True)
class CumberlandPaths:
    ks_v5: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\ks_v5.gpkg")
    drn_path: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\drn.gpkg")
    idomain_path: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\idomain.gpkg")
    vor_algo: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\cumb_v14c_existing_algomesh_v2.1.vor")
    iheads_path: Path = Path(r"C:\Users\lukem\mf6\cumb_v6z_botm14V\iheads.hds")
    pit2: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\pit_2_footprint.gpkg")
    pit4: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Herrera\pit4.gpkg")
    rej: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Herrera\min_elevs_to_fix_inf_rej.gpkg")
    mfr_polys: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\mfr.gpkg")
    calib_locs: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\calibration_locs.gpkg")
    calib_obs: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Tables\calib_observations.xlsx")
    forcing_table: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Tables\all_data_summary.xlsx")
    merged_wells: Path = Path(r"C:\Users\lukem\QGIS\SHP\Cumberland\explo cumberland merged.gpkg")
    boring_summary: Path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boring Summary.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a transient Cumberland first-slice PEST workflow with one initial "
            "steady-state stress period followed by 12 monthly transient periods."
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
        default="cumberland_transient_observed_pest_first_slice",
        help="Prefix used for the timestamped artifact directory.",
    )
    parser.add_argument(
        "--source-start-per",
        type=int,
        default=57,
        help=(
            "First monthly Cumberland source stress period to map onto model period 1. "
            "Model period 0 remains a separate steady-state spinup period."
        ),
    )
    parser.add_argument(
        "--n-transient",
        type=int,
        default=12,
        help="Number of monthly transient stress periods after the initial steady period.",
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
        help="Uniform recharge applied during the initial steady-state spinup period.",
    )
    parser.add_argument(
        "--supplemental-weight",
        type=float,
        default=0.25,
        help="Observation weight assigned to one-time supplemental head targets.",
    )
    parser.add_argument(
        "--min-transient-recharge",
        type=float,
        default=1e-4,
        help="Minimum monthly transient recharge floor used to improve dry-season stability.",
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


def load_and_adjust_voronoi(paths: CumberlandPaths):
    with open(paths.iheads_path, "rb") as handle:
        iheads = pickle.load(handle)

    with open(paths.vor_algo, "rb") as handle:
        vor = pickle.load(handle)

    vor.gdf_topbtm[0] = vor.gdf_topbtm[0].astype(float)
    vor.gdf_topbtm[1] = vor.gdf_topbtm[1].astype(float)
    reconciled = vor.reconcile_surfaces(min_sep=5, trigger_sep=5, which="top")
    vor.gdf_topbtm[0] = reconciled[0].astype(float)
    vor.gdf_topbtm[1] = reconciled[1].astype(float)

    inactive_cells = vor.get_vor_cells_as_series(paths.idomain_path).explode().tolist()
    idomain = pd.Series(1 for _ in range(vor.ncpl))
    idomain.iloc[inactive_cells] = 0

    pit2_cells = vor.get_vor_cells_as_series(paths.pit2)[0]
    tb2 = vor.gdf_topbtm.loc[pit2_cells]
    pit2_thickness = tb2.loc[:, 0] - tb2.loc[:, 1]
    thin_cells = pit2_thickness[pit2_thickness < 25].index
    vor.gdf_topbtm.loc[thin_cells, 0] = vor.gdf_topbtm.loc[thin_cells, 1] + 25

    pit4_cells = vor.get_vor_cells_as_series(paths.pit4)[0]
    tb4 = vor.gdf_topbtm.loc[pit4_cells]
    pit4_thickness = tb4.loc[:, 0] - tb4.loc[:, 1]
    thin_cells = pit4_thickness[pit4_thickness < 20].index
    vor.gdf_topbtm.loc[thin_cells, 0] = vor.gdf_topbtm.loc[thin_cells, 1] + 20

    rej_cells = vor.get_vor_cells_as_series(paths.rej, name_field="name")
    rej_df = pd.concat([read_shp_gpkg(paths.rej).set_index("name"), rej_cells], axis=1)
    for name, row in rej_df.iterrows():
        cells_to_check = vor.gdf_topbtm.loc[row.cells, 0]
        low_cells = cells_to_check.loc[cells_to_check < row.thickness].index.to_list()
        vor.gdf_topbtm.loc[low_cells, 0] = row.thickness

    liz_cells = vor.get_vor_cells_as_series(read_shp_gpkg(paths.mfr_polys).loc[8, "geometry"])[0]
    lizmfr_cells = [cell for cell in liz_cells if cell not in inactive_cells]
    vor.gdf_topbtm.loc[lizmfr_cells, 0] += 10

    reconciled = vor.reconcile_surfaces(min_sep=5, trigger_sep=5, which="top")
    vor.gdf_topbtm[0] = reconciled[0].astype(float)
    vor.gdf_topbtm[1] = reconciled[1].astype(float)

    top = vor.gdf_topbtm[0].to_numpy(dtype=float)
    botm = vor.gdf_topbtm[1].to_numpy(dtype=float)
    return vor, iheads, idomain, inactive_cells, top, botm


def select_source_window(
    obs_frame: pd.DataFrame,
    forcing_frame: pd.DataFrame,
    *,
    total_periods: int,
    source_start_per: int | None,
) -> pd.DataFrame:
    obs_frame = obs_frame.sort_values("per").reset_index(drop=True)
    forcing_frame = forcing_frame.sort_values("per").reset_index(drop=True)
    forcing_pers = set(pd.to_numeric(forcing_frame["per"], errors="coerce").dropna().astype(int))

    if source_start_per is not None:
        source_pers = list(range(source_start_per, source_start_per + total_periods))
        selected = obs_frame.loc[obs_frame["per"].isin(source_pers)].copy()
        if sorted(selected["per"].astype(int).tolist()) != source_pers:
            raise ValueError(
                "Requested source-start period does not produce a complete contiguous window "
                f"of {total_periods} periods in calib_observations.xlsx."
            )
        if not set(source_pers).issubset(forcing_pers):
            raise ValueError(
                "Requested source-start period does not align with all_data_summary.xlsx "
                f"for all periods: {source_pers}"
            )
        return selected.sort_values("per").reset_index(drop=True)

    well_columns = [column for column in obs_frame.columns if column not in {"per", "Deep Lake"}]
    obs_frame["coverage"] = obs_frame[well_columns].notna().sum(axis=1)
    periods = obs_frame["per"].astype(int).tolist()

    best_indices: tuple[int, int] | None = None
    best_score: tuple[int, int] | None = None
    for start_idx in range(0, len(obs_frame) - total_periods + 1):
        window = obs_frame.iloc[start_idx : start_idx + total_periods]
        window_pers = window["per"].astype(int).tolist()
        if window_pers != list(range(window_pers[0], window_pers[0] + total_periods)):
            continue
        if not set(window_pers).issubset(forcing_pers):
            continue
        score = (int(window["coverage"].sum()), int(window_pers[-1]))
        if best_score is None or score > best_score:
            best_score = score
            best_indices = (start_idx, start_idx + total_periods)

    if best_indices is None:
        raise ValueError("Could not find a contiguous observed/forcing window for the requested period count.")

    start_idx, end_idx = best_indices
    return obs_frame.iloc[start_idx:end_idx].copy().drop(columns="coverage").reset_index(drop=True)


def build_single_observation_rows(
    obs_frame: pd.DataFrame,
    *,
    allowed_names: set[str],
    selected_names: set[str],
    weight: float,
) -> list[dict[str, float | int | str]]:
    """Return one-time target rows for wells with exactly one observed head."""

    rows: list[dict[str, float | int | str]] = []
    well_columns = [column for column in obs_frame.columns if column not in {"per", "Deep Lake"}]
    for name in well_columns:
        if name not in allowed_names or name in selected_names:
            continue
        series = pd.to_numeric(obs_frame[name], errors="coerce").dropna()
        if len(series) != 1:
            continue
        rows.append(
            {
                "per": 0,
                "name": str(name),
                "head": float(series.iloc[0]),
                "weight_override": float(weight),
            }
        )
    return rows


def build_head_targets(
    paths: CumberlandPaths,
    obs_frame: pd.DataFrame,
    selected_obs: pd.DataFrame,
    *,
    supplemental_weight: float,
) -> tuple[mf.HeadTargets, pd.DataFrame]:
    locs = gpd.read_file(paths.calib_locs).copy()
    locs["name"] = locs["ExploName"].astype(str)
    locs["layer_num"] = 0
    locs["weight"] = 1.0

    allowed_names = set(locs["name"].tolist())
    rows: list[dict[str, float | int | str]] = []
    well_columns = [column for column in selected_obs.columns if column != "per"]
    for transient_offset, (_, row) in enumerate(selected_obs.iterrows(), start=1):
        for name in well_columns:
            if name not in allowed_names:
                continue
            value = row[name]
            if pd.isna(value):
                continue
            rows.append({"per": transient_offset, "name": str(name), "head": float(value)})
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

    supplemental_rows = [
        {"per": 0, "name": str(row.name), "head": float(row.gw_elev)}
        for row in boring.itertuples(index=False)
        if pd.notna(row.gw_elev)
    ]
    rows.extend(supplemental_rows)

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


def _cells_from_geometry(vor, geometry) -> list[int]:
    cells = vor.get_vor_cells_as_series(geometry)
    if isinstance(cells, pd.Series):
        exploded = cells.explode()
        return sorted({int(cell) for cell in exploded.dropna().tolist()})
    return sorted({int(cell) for cell in cells})


def build_recharge_dict(
    vor,
    idomain: pd.Series,
    recharge_values: list[float],
    *,
    mfr_path: Path,
    default_mfr_factor: float,
    steady_recharge: float,
) -> dict[int, list[list[tuple[int, int] | float]]]:
    active_cells = sorted(idomain.loc[idomain == 1].index.astype(int).tolist())
    rch_dict: dict[int, list[list[tuple[int, int] | float]]] = {}

    mfr = gpd.read_file(mfr_path)
    mfr_groups: list[tuple[list[int], float]] = []
    for _, row in mfr.iterrows():
        factor = row.get("factor", default_mfr_factor)
        factor = default_mfr_factor if pd.isna(factor) or float(factor) <= 0 else float(factor)
        cells = [cell for cell in _cells_from_geometry(vor, row.geometry) if cell in active_cells]
        if cells:
            mfr_groups.append((cells, factor))

    recharge_by_per = [max(float(steady_recharge), 0.0)] + [max(float(value), 0.0) for value in recharge_values]
    for per, per_value in enumerate(recharge_by_per):
        cell_values = {cell: per_value for cell in active_cells}
        for cells, factor in mfr_groups:
            scaled = per_value * factor
            for cell in cells:
                cell_values[cell] = scaled
        rch_dict[per] = [[(0, int(cell)), float(cell_values[cell])] for cell in active_cells]
    return rch_dict


def configure_solver(model):
    ims = model.sim.ims
    ims.under_relaxation_momentum = 0.0
    ims.backtracking_tolerance = 1.1
    ims.backtracking_number = 20
    ims.preconditioner_levels = 5
    ims.preconditioner_drop_tolerance = 1e-4
    ims.outer_dvclose = 1e-5
    ims.inner_dvclose = 1e-6
    ims.under_relaxation_theta = 0.7
    ims.under_relaxation_kappa = 0.1
    ims.under_relaxation_gamma = 0.2
    ims.backtracking_reduction_factor = 0.2
    ims.backtracking_residual_limit = 10
    ims.rcloserecord = 1e-3
    ims.complexity = "MODERATE"
    ims.under_relaxation = "dbd"
    ims.linear_acceleration = "BICGSTAB"
    ims.outer_maximum = 600
    ims.inner_maximum = 300
    ims.relaxation_factor = 1.0


def build_model(
    *,
    paths: CumberlandPaths,
    model_workspace: Path,
    selected_forcing: pd.DataFrame,
    pp_model_name: str,
    mfr_factor: float,
    steady_recharge: float,
):
    vor, iheads, idomain, inactive_cells, top, botm = load_and_adjust_voronoi(paths)

    nper = len(selected_forcing) + 1
    model = mf.SimulationBase(
        name=pp_model_name,
        nper=nper,
        vor=vor,
        mf_folder_path=model_workspace,
    )
    transient_dates = pd.to_datetime(selected_forcing["Month"]).tolist()
    steady_date = transient_dates[0] - pd.offsets.MonthBegin(1)
    model.per_dates = [steady_date] + transient_dates

    DisvGrid(vor=vor, model=model, nlay=1, top=top, bottom=botm, idomain=idomain)

    month_days = pd.to_datetime(selected_forcing["Month"]).dt.days_in_month.astype(float).tolist()
    period_data = [[1.0, 15, 1.1]] + [[days, 10, 1.05] for days in month_days]
    TemporalDiscretization(model=model, period_data=period_data)
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
        sto_transient={1: True},
    )

    recharge_values = selected_forcing["recharge_fpd"].astype(float).tolist()
    rch_dict = build_recharge_dict(
        vor,
        idomain,
        recharge_values,
        mfr_path=paths.mfr_polys,
        default_mfr_factor=mfr_factor,
        steady_recharge=steady_recharge,
    )
    Recharge(model=model, vor=vor, rch_dict=rch_dict)

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


def prepare_pest_project(
    *,
    model,
    workspace_root: Path,
    pest_workspace: Path,
    ks_v5: Path,
    drn_path: Path,
    targets: mf.HeadTargets,
    start_datetime: str,
    pp_spacing: float,
    include_drn_cond: bool,
):
    raw_drn = gpd.read_file(drn_path).copy()
    raw_drn["par_name"] = raw_drn["name"].astype(str).str.lower().map(
        lambda s: re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    )
    param_drn = workspace_root / "drn_param_source.gpkg"
    if param_drn.exists():
        param_drn.unlink()
    raw_drn.to_file(param_drn, driver="GPKG")

    pest = mf.PestProject(
        model=model,
        name=f"{model.name}_pest",
        workspace=pest_workspace,
        start_datetime=start_datetime,
    )
    pest.add_parameter(
        mf.KPilotPointParameter(
            name="hk",
            source=mf.VectorParameterSource(
                path=ks_v5,
                value_column="k",
                feature_id_column="name",
                zone_column="name",
                layer_column="layer",
                crs=2926,
            ),
            bounds=(0.25, 4.0),
            bounds_mode="multiplier",
            transform="log",
            pp_spacing=pp_spacing,
            geostruct=mf.ExpGeoStruct(range=4000.0, transform="log"),
        )
    )
    if include_drn_cond:
        pest.add_parameter(
            mf.DrainConductanceParameter(
                name="drn_cond",
                source=mf.VectorParameterSource(
                    path=param_drn,
                    value_column="cond",
                    feature_id_column="par_name",
                    layer_column="layer",
                    crs=2926,
                ),
                bounds=(0.25, 4.0),
                bounds_mode="multiplier",
                transform="log",
            )
        )
    pest.add_observation(mf.HeadTargetObservationSpec(targets=targets))
    return pest


def run_pestpp_live(pest_workspace: Path, control_file: str, exe_name: str):
    print(f"Running PEST++ in: {pest_workspace}")
    print(f"Control file: {control_file}")
    proc = subprocess.Popen(
        [exe_name, control_file],
        cwd=pest_workspace,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="")
    proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, proc.args)


def run_pestpp_parallel(
    *,
    workspace_root: Path,
    pest_workspace: Path,
    control_file: str,
    exe_name: str,
    n_workers: int,
    keep_workers: bool,
) -> Path:
    import pyemu

    worker_root = workspace_root / "workers"
    worker_root.mkdir(parents=True, exist_ok=True)
    master_dir = workspace_root / "pest_master"
    print(f"Running parallel PEST++ with {n_workers} workers")
    print(f"Template workspace: {pest_workspace}")
    print(f"Master workspace: {master_dir}")
    print(f"Worker root: {worker_root}")

    pyemu.os_utils.start_workers(
        str(pest_workspace),
        exe_name,
        control_file,
        num_workers=int(n_workers),
        worker_root=str(worker_root),
        master_dir=str(master_dir),
        cleanup=not keep_workers,
        verbose=True,
    )
    return master_dir


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
    forcing_frame = pd.read_excel(paths.forcing_table)
    total_periods = args.n_transient
    selected_obs = select_source_window(
        obs_frame,
        forcing_frame,
        total_periods=total_periods,
        source_start_per=args.source_start_per,
    )
    source_pers = selected_obs["per"].astype(int).tolist()
    selected_forcing = forcing_frame.loc[forcing_frame["per"].isin(source_pers), ["Month", "per", "precip", "ET"]].copy()
    selected_forcing = selected_forcing.sort_values("per").reset_index(drop=True)
    selected_forcing["Month"] = pd.to_datetime(selected_forcing["Month"])
    selected_forcing["recharge_fpd"] = (selected_forcing["precip"] - selected_forcing["ET"]).clip(
        lower=float(args.min_transient_recharge)
    )
    selected_forcing_path = workspace_root / "selected_forcing.csv"
    selected_forcing.to_csv(selected_forcing_path, index=False)

    print(f"Selected Cumberland source periods: {source_pers}")
    print(
        "Source window dates: "
        f"{selected_forcing['Month'].iloc[0].date()} to {selected_forcing['Month'].iloc[-1].date()}"
    )

    targets, target_values = build_head_targets(
        paths,
        obs_frame,
        selected_obs,
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
            "source_periods": source_pers,
            "source_start_per": int(source_pers[0]),
            "source_end_per": int(source_pers[-1]),
            "n_transient": int(args.n_transient),
            "steady_recharge": float(args.steady_recharge),
            "min_transient_recharge": float(args.min_transient_recharge),
            "supplemental_weight": float(args.supplemental_weight),
            "pp_spacing": float(args.pp_spacing),
            "k_only": bool(args.k_only),
            "fast_mode": bool(args.fast_mode),
            "head_target_locations": str(target_locations_path),
            "head_target_values": str(target_values_path),
            "selected_forcing": str(selected_forcing_path),
        }
    )
    run_info_path.write_text(json.dumps(run_info, indent=2), encoding="utf-8")
    print(
        f"Built {len(target_values)} transient head targets across "
        f"{target_values['name'].nunique()} wells and {target_values['per'].nunique()} model periods."
    )

    model = build_model(
        paths=paths,
        model_workspace=model_workspace,
        selected_forcing=selected_forcing,
        pp_model_name="cpre_obs_13p",
        mfr_factor=args.mfr_factor,
        steady_recharge=args.steady_recharge,
    )
    success, _ = model.run_simulation()
    if not success:
        raise RuntimeError("Baseline transient Cumberland model did not solve successfully.")

    baseline_stats = targets.stats(model)
    print("\nBaseline residual stats")
    print(baseline_stats.to_string(index=False))

    pest = prepare_pest_project(
        model=model,
        workspace_root=workspace_root,
        pest_workspace=pest_workspace,
        ks_v5=paths.ks_v5,
        drn_path=paths.drn_path,
        targets=targets,
        start_datetime=str(selected_forcing["Month"].iloc[0].date()),
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
