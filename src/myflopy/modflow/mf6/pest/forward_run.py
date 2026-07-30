"""Forward-run helpers for ``myflopy`` PEST/pyEMU workflows."""

from __future__ import annotations

import numpy as np
import pandas as pd


def apply_pilotpoints_to_array(parameter_csv, points_meta_csv, cells_meta_csv, array_file,
                               lower_limit=-1.0e30, upper_limit=1.0e30):
    """Interpolate pilot-point multipliers onto a flattened array layer file (IDW).

    Runs as a pre-model command in the native forward run: reads the
    template-populated pilot-point values, inverse-distance-weights them onto
    each cell, multiplies the captured base values, clamps to the physical
    limits, and rewrites the external layer file (one value per line) that MF6
    reads.

    Nothing here is K-specific -- the base values are whichever array the target
    names (``base_value`` in the cells CSV). The parameters were called ``k_file``
    and ``base_k`` until 2026-07-30, which is how the build side came to hardwire
    NPF K for every target and write it into K33.

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
        raise ValueError(f"Missing pilot-point values in {parameter_csv!r}.")
    px = points["x"].to_numpy(dtype=float)
    py = points["y"].to_numpy(dtype=float)
    pv = points["pv"].to_numpy(dtype=float)
    cx = cells["x"].to_numpy(dtype=float)
    cy = cells["y"].to_numpy(dtype=float)
    base = cells["base_value"].to_numpy(dtype=float)
    cell_ids = cells["cell"].to_numpy(dtype=int)
    # The captured array file is treated as a model *output* by pyEMU, so the
    # input .txt no longer exists -- rebuild the whole layer array from the
    # captured base values and write the file MF6 reads.
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
    np.savetxt(array_file, arr.reshape(-1, 1), fmt="%.10E")


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


def _write_conc_target_csv(model_name, mapping_csv, output_csv):
    """Write simulated concentrations by stress period for named target locations.

    The concentration twin of ``_write_head_target_csv``, and deliberately the
    same shape: open the binary directly and index it. The named-series helper
    would instead reload the whole simulation on every forward run, and its
    ``all_conc`` read is kind-gated, so it would have to resolve to the GWT model
    rather than the flow model it loads by default.

    ``model_name`` is the TRANSPORT model's name -- the ``.ucn`` sits beside the
    flow model's ``.hds`` at the workspace root because the simulation is flat.
    """

    import flopy

    mapping = pd.read_csv(mapping_csv)
    ucn = flopy.utils.HeadFile(f"{model_name}.ucn", text="concentration")
    try:
        by_period = {}
        for kstp, kper in ucn.get_kstpkper():
            by_period[int(kper)] = (int(kstp), int(kper))
        rows = []
        for per in sorted(by_period):
            data = np.asarray(ucn.get_data(kstpkper=by_period[per]), dtype=float)
            row = {"per": int(per)}
            for target in mapping.itertuples(index=False):
                layer_values = np.asarray(data[int(target.layer)], dtype=float).reshape(-1)
                row[str(target.name)] = float(layer_values[int(target.cell)])
            rows.append(row)
        pd.DataFrame(rows).to_csv(output_csv, index=False)
    finally:
        ucn.close()


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

    # Deliberately duplicated: pyEMU extracts this function's SOURCE into the
    # generated forward_run.py, whose header does not import Path — the import
    # must live inside the function to survive extraction.
    from pathlib import Path  # noqa: F811

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
        raise ValueError(f"Unsupported named-series observation kind {kind!r}.")
    targets = builders[kind_key](locations=locations, values=None)
    targets.simulated_series(model).to_csv(output_csv, index=False)

