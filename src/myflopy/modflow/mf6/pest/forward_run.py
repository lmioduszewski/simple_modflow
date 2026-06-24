"""Forward-run helpers for ``myflopy`` PEST/pyEMU workflows."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def apply_pilotpoints_to_array(parameter_csv, points_meta_csv, cells_meta_csv, k_file,
                               lower_limit=-1.0e30, upper_limit=1.0e30):
    """Interpolate pilot-point multipliers onto a flattened K layer file (IDW).

    Runs as a pre-model command in the native forward run: reads the
    template-populated pilot-point values, inverse-distance-weights them onto
    each cell, multiplies the captured base K, clamps to the physical limits,
    and rewrites the external K layer file (one value per line) that MF6 reads.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so it must not depend on other module-level helpers.
    """

    import numpy as np
    import pandas as pd

    pv_frame = pd.read_csv(parameter_csv)
    values = pd.Series(
        pd.to_numeric(pv_frame["value"], errors="coerce").to_numpy(),
        index=pv_frame["parnme"].astype(str),
    )
    points = pd.read_csv(points_meta_csv)
    cells = pd.read_csv(cells_meta_csv)
    points["pv"] = points["parnme"].astype(str).map(values)
    if points["pv"].isna().any():
        raise ValueError("Missing pilot-point values in %r." % parameter_csv)
    px = points["x"].to_numpy(dtype=float)
    py = points["y"].to_numpy(dtype=float)
    pv = points["pv"].to_numpy(dtype=float)
    cx = cells["x"].to_numpy(dtype=float)
    cy = cells["y"].to_numpy(dtype=float)
    base = cells["base_k"].to_numpy(dtype=float)
    cell_ids = cells["cell"].to_numpy(dtype=int)
    # The captured K file is treated as a model *output* by pyEMU, so the input
    # .txt no longer exists -- rebuild the whole layer array from the captured
    # base K and write the file MF6 reads.
    arr = np.zeros(int(cell_ids.max()) + 1, dtype=float)
    for index in range(len(cell_ids)):
        d2 = (px - cx[index]) ** 2 + (py - cy[index]) ** 2
        zero = np.where(d2 == 0)[0]
        if zero.size:
            factor = float(pv[zero[0]])
        else:
            weights = 1.0 / np.maximum(d2, 1.0e-12)
            factor = float(np.sum(weights * pv) / np.sum(weights))
        value = base[index] * factor
        arr[cell_ids[index]] = min(max(value, float(lower_limit)), float(upper_limit))
    np.savetxt(k_file, arr.reshape(-1, 1), fmt="%.10E")


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


def write_named_series_targets(kind, locations_file, output_csv, sim_ws="."):
    """Regenerate one lake/SFR/DRN simulated-target CSV in the native forward run.

    Post-model command for named-series observations (lake stage, SFR stage/flow,
    DRN seepage): it reloads the just-run MF6 workspace, rebuilds the target set
    from the saved definition file (the ``*_target_locations`` CSV/GPKG written at
    build time), and writes the simulated series pyEMU reads back as observations.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so it loads the model and reads the saved locations inline rather
    than calling other module-level helpers.
    """

    from pathlib import Path

    import myflopy as mf

    locations_path = Path(locations_file)
    if locations_path.suffix.lower() == ".gpkg":
        import geopandas as gpd

        locations = gpd.read_file(locations_path)
    else:
        locations = pd.read_csv(locations_path)
        if "cells" in locations.columns:
            locations["cells"] = locations["cells"].fillna("").apply(
                lambda text: [int(v) for v in str(text).split(",") if str(v).strip() != ""]
            )

    model = mf.load_mf6_run(Path(sim_ws), verbosity_level=0)
    model.load_all()

    builders = {
        "lake_stage": mf.LakeStageTargets,
        "sfr_stage": mf.SfrStageTargets,
        "sfr_flow": mf.SfrFlowTargets,
        "drn_flow": mf.DrnFlowTargets,
    }
    kind_key = str(kind).strip().lower()
    if kind_key not in builders:
        raise ValueError("Unsupported named-series observation kind %r." % (kind,))
    targets = builders[kind_key](locations=locations, values=None)
    targets.simulated_series(model).to_csv(output_csv, index=False)

