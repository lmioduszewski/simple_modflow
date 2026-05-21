"""Forward-run helpers for ``simple_modflow`` PEST/pyEMU workflows."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _expand_cellid_columns(frame):
    """Expand FloPy cellid tuples into zero-based ``layer`` and ``cell`` columns."""

    layers = []
    cells = []
    for cellid in frame["cellid"].tolist():
        if isinstance(cellid, tuple):
            layers.append(int(cellid[0]))
            cells.append(int(cellid[1]))
        else:
            layers.append(0)
            cells.append(int(cellid))
    frame = frame.copy()
    frame["layer"] = layers
    frame["cell"] = cells
    return frame


def _build_drn_frame(package):
    """Normalize the MF6 DRN package data into one DataFrame."""

    frames = []
    for per, data in package.stress_period_data.data.items():
        frame = pd.DataFrame(data).copy()
        frame = _expand_cellid_columns(frame)
        frame["per"] = int(per)
        frame["row_index"] = np.arange(len(frame), dtype=int)
        frames.append(frame.loc[:, ["per", "row_index", "cellid", "layer", "cell", "elev", "cond"]])
    if not frames:
        return pd.DataFrame(columns=["per", "row_index", "cellid", "layer", "cell", "elev", "cond"])
    return pd.concat(frames, ignore_index=True)


def _rebuild_drn_stress_period_data(frame):
    """Convert a normalized DRN frame back into MF6 stress-period data."""

    result = {}
    for per, per_frame in frame.sort_values(["per", "row_index"]).groupby("per"):
        rows = []
        for row in per_frame.itertuples(index=False):
            rows.append([tuple(row.cellid), float(row.elev), float(row.cond)])
        result[int(per)] = rows
    return result


def _load_parameter_values(path):
    """Load template-populated parameter values from a CSV file."""

    frame = pd.read_csv(path)
    if not {"parnme", "value"}.issubset(frame.columns):
        raise ValueError(f"Parameter CSV {path!s} must contain 'parnme' and 'value' columns.")
    return pd.Series(pd.to_numeric(frame["value"], errors="coerce").to_numpy(), index=frame["parnme"].astype(str))


def _apply_drain_specs(gwf, specs):
    """Apply drain parameter specs to the loaded MF6 model."""

    if not specs:
        return
    current = _build_drn_frame(gwf.drn)
    for spec in specs:
        support = pd.read_csv(spec["support_csv"])
        parameter_values = _load_parameter_values(spec["parameter_csv"])
        support["parameter_value"] = support["parnme"].map(parameter_values)
        support = support.dropna(subset=["parameter_value"])
        updated = support.loc[:, ["per", "row_index"]].copy()
        if spec["kind"] == "drn_elev":
            updated["elev_new"] = pd.to_numeric(support["base_elev"], errors="coerce") + pd.to_numeric(
                support["parameter_value"], errors="coerce"
            )
            current = current.merge(updated, on=["per", "row_index"], how="left")
            mask = current["elev_new"].notna()
            current.loc[mask, "elev"] = current.loc[mask, "elev_new"]
            current = current.drop(columns=["elev_new"])
        elif spec["kind"] == "drn_cond":
            updated["cond_new"] = pd.to_numeric(support["base_cond"], errors="coerce") * pd.to_numeric(
                support["parameter_value"], errors="coerce"
            )
            current = current.merge(updated, on=["per", "row_index"], how="left")
            mask = current["cond_new"].notna()
            current.loc[mask, "cond"] = current.loc[mask, "cond_new"]
            current = current.drop(columns=["cond_new"])
        else:
            raise ValueError(f"Unsupported drain parameter kind {spec['kind']!r}.")
    gwf.drn.stress_period_data.set_data(_rebuild_drn_stress_period_data(current))


def _interpolate_idw(x, y, px, py, values):
    """Inverse-distance weighting for one target point."""

    dx = px - x
    dy = py - y
    dist2 = dx * dx + dy * dy
    exact = dist2 == 0
    if exact.any():
        return float(values[np.where(exact)[0][0]])
    weights = 1.0 / np.maximum(dist2, 1.0e-12)
    return float(np.sum(weights * values) / np.sum(weights))


def _apply_k_specs(gwf, specs):
    """Apply pilot-point K multiplier specs to the loaded MF6 model."""

    if not specs:
        return
    k_array = np.asarray(gwf.npf.k.array, dtype=float).copy()
    for spec in specs:
        point_meta = pd.read_csv(spec["points_meta_csv"])
        cell_meta = pd.read_csv(spec["cells_meta_csv"])
        parameter_values = _load_parameter_values(spec["parameter_csv"])
        point_meta["parameter_value"] = point_meta["parnme"].map(parameter_values)
        if point_meta["parameter_value"].isna().any():
            raise ValueError(f"Missing pilot-point values in {spec['parameter_csv']!r}.")
        for (layer, zone), group in cell_meta.groupby(["layer", "zone"], dropna=False):
            point_group = point_meta[(point_meta["layer"] == layer) & (point_meta["zone"] == zone)]
            if point_group.empty:
                point_group = point_meta[point_meta["layer"] == layer]
            px = point_group["x"].to_numpy(dtype=float)
            py = point_group["y"].to_numpy(dtype=float)
            pv = point_group["parameter_value"].to_numpy(dtype=float)
            for row in group.itertuples(index=False):
                factor = _interpolate_idw(float(row.x), float(row.y), px, py, pv)
                if spec.get("parameter_space", "multiplier") == "multiplier":
                    new_value = float(row.base_k) * factor
                else:
                    new_value = factor
                k_array[int(layer), int(row.cell)] = new_value
    gwf.npf.k.set_data(k_array)


def _write_head_target_csv(model_name, mapping_csv, output_csv):
    """Write simulated heads by stress period for named target locations."""

    import flopy

    mapping = pd.read_csv(mapping_csv)
    hds = flopy.utils.HeadFile(f"{model_name}.hds")
    try:
        by_period = {}
        for kstp, kper in hds.get_kstpkper():
            by_period[int(kper)] = (int(kstp), int(kper))
        rows = []
        for per in sorted(by_period):
            data = np.asarray(hds.get_data(kstpkper=by_period[per]), dtype=float)
            row = {"per": int(per)}
            for target in mapping.itertuples(index=False):
                layer_values = np.asarray(data[int(target.layer)], dtype=float).reshape(-1)
                row[str(target.name)] = float(layer_values[int(target.cell)])
            rows.append(row)
        pd.DataFrame(rows).to_csv(output_csv, index=False)
    finally:
        hds.close()


def _write_simulation_with_retry(sim, *, attempts=4, delay_seconds=1.0):
    """Write the MF6 simulation with a small retry loop for transient Windows errors."""

    last_error = None
    for attempt in range(1, attempts + 1):
        try:
            sim.write_simulation(silent=True)
            return
        except OSError as exc:
            last_error = exc
            if attempt == attempts:
                break
            time.sleep(delay_seconds)
    raise last_error


def apply_pest_forward_run(config_path="pest_forward_config.json"):
    """Apply current parameter values, run MF6, and regenerate observation CSVs."""

    import flopy
    import json
    from pathlib import Path

    config = json.loads(Path(config_path).read_text(encoding="utf-8"))
    sim = flopy.mf6.MFSimulation.load(sim_ws=".", verbosity_level=0)
    if config.get("exe_name"):
        sim.exe_name = config["exe_name"]
    model_name = config.get("model_name")
    gwf = sim.get_model(model_name) if model_name else sim.get_model()

    _apply_k_specs(gwf, config.get("k_specs", []))
    _apply_drain_specs(gwf, config.get("drain_specs", []))

    _write_simulation_with_retry(sim)
    success, buff = sim.run_simulation(silent=True)
    if not success:
        raise RuntimeError(f"MF6 forward run failed: {buff}")

    for target in config.get("head_target_outputs", []):
        _write_head_target_csv(
            model_name=gwf.name,
            mapping_csv=target["mapping_csv"],
            output_csv=target["output_csv"],
        )

    return True
