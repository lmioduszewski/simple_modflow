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


def _write_head_target_csv(model_name, mapping_csv, output_csv, perlen=None, sim_ws="."):
    """Write simulated heads at each stress period's END for named target locations.

    Post-model command in the native forward run, and the ONE writer of this file:
    the build writes the base run's copy by calling this same function (ledger
    340). PEST's instruction file reads the CSV by POSITION -- the n-th comma of
    the m-th line -- and names what it reads after the BASE copy's header, so the
    two copies must hold the same columns in the same order. They did not: the
    base copy was a pivot, whose columns pandas SORTS, and this wrote mapping
    order, so targets named W2, W1, W10 handed PEST W1's head as W10's in every
    realization while ``compare()`` paired them correctly.

    A period's head is the one that ENDS it, picked by the one rule every target
    reads (``myflopy.modflow.mf6.period_end``, ledger 336) on the head file's own
    clock.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so its imports live inside it.

    Parameters
    ----------
    model_name : str
        The flow model's name; its head file is ``<model_name>.hds`` in
        ``sim_ws``.
    mapping_csv : str or pathlib.Path
        The target map the build wrote: ``obsname``/``name``, ``layer``, ``cell``.
        Its row order IS the column order written, and its labels are read as
        text, verbatim (a well named ``007`` is not the number 7, ledger 341).
    output_csv : str or pathlib.Path
        The simulated CSV to write: ``per``, then one column per mapping row
        (its ``obsname``, else its ``name``) in mapping order; one row per stress
        period the head file saved.
    perlen : list of float, optional
        Each stress period's length, recorded from the TDIS at build time. With
        it a period whose final step the head file does not hold is REFUSED --
        the run fails rather than hand PEST a value from before the period ended;
        ATS can move a ``STEPS`` save off the end in one realization and not
        another. ``None`` (a forward run generated before ledger 336) takes each
        period's last saved step unchecked.
    sim_ws : str or pathlib.Path, default "."
        The directory holding the head file: the forward run's working
        directory, or the model workspace when the build writes the base copy.

    Raises
    ------
    ValueError
        With ``perlen``, when a period the head file saved lacks its final step.
    """

    # Imported here, not hoisted: pyEMU copies only this function's source.
    from pathlib import Path

    import flopy

    from myflopy.modflow.mf6.period_end import period_end_steps, saved_steps

    # Labels read VERBATIM: pandas reads a well named 007 as the integer 7 and one
    # named NA as a missing value, so the columns written would be 7 and nan while
    # the target map, the targets and the base table's header say 007 and NA
    # (ledger 341). A converter sees the raw text before either reading.
    mapping = pd.read_csv(mapping_csv, converters={"obsname": str, "name": str})
    hds = flopy.utils.HeadFile(str(Path(sim_ws) / f"{model_name}.hds"))
    try:
        ends = period_end_steps(
            saved_steps(hds), perlen=perlen, file=f"{model_name}.hds", saves="HEAD"
        )
        rows = []
        for per, kstpkper in ends.items():
            data = np.asarray(hds.get_data(kstpkper=kstpkper), dtype=float)
            row = {"per": int(per)}
            for target in mapping.itertuples(index=False):
                layer_values = np.asarray(data[int(target.layer)], dtype=float).reshape(-1)
                # The mapping's own label: one column per (name, layer). Keyed on
                # the bare name, two screens of one well overwrote each other here
                # while the build read the other one. A mapping written before the
                # label existed has no such column, and its name still serves.
                label = str(getattr(target, "obsname", "") or target.name)
                row[label] = float(layer_values[int(target.cell)])
            rows.append(row)
        pd.DataFrame(rows).to_csv(output_csv, index=False)
    finally:
        hds.close()


def _write_conc_target_csv(model_name, mapping_csv, output_csv, perlen=None, sim_ws="."):
    """Write simulated concentrations at each stress period's END for named target locations.

    The concentration twin of ``_write_head_target_csv``, and deliberately the
    same shape: open the binary directly and index it, at the record the one
    period-end rule picks (ledger 336). The named-series helper would instead
    reload the whole simulation on every forward run, and its ``all_conc`` read
    is kind-gated, so it would have to resolve to the GWT model rather than the
    flow model it loads by default. Like it, this is the ONE writer of its file:
    the build writes the base run's copy through it (ledger 340), so the columns
    PEST's positional instruction file names and the ones every realization
    writes are one list in one order.

    Parameters
    ----------
    model_name : str
        The TRANSPORT model's name -- the ``.ucn`` sits beside the flow model's
        ``.hds`` in ``sim_ws`` because the simulation is flat.
    mapping_csv : str or pathlib.Path
        The target map the build wrote: ``obsname``/``name``, ``layer``, ``cell``.
        Its row order IS the column order written, and its labels are read as
        text, verbatim (a well named ``007`` is not the number 7, ledger 341).
    output_csv : str or pathlib.Path
        The simulated CSV to write: ``per``, then one column per mapping row in
        mapping order; one row per stress period the file saved.
    perlen : list of float, optional
        Each stress period's length, recorded from the TDIS at build time; with
        it a period whose final step the file does not hold is refused (see
        ``_write_head_target_csv``). ``None`` takes each period's last saved step
        unchecked.
    sim_ws : str or pathlib.Path, default "."
        The directory holding the ``.ucn``: the forward run's working directory,
        or the model workspace when the build writes the base copy.

    Raises
    ------
    ValueError
        With ``perlen``, when a period the concentration file saved lacks its
        final step.
    """

    # Imported here, not hoisted: pyEMU copies only this function's source.
    from pathlib import Path

    import flopy

    from myflopy.modflow.mf6.period_end import period_end_steps, saved_steps

    # Labels read VERBATIM: pandas reads a well named 007 as the integer 7 and one
    # named NA as a missing value, so the columns written would be 7 and nan while
    # the target map, the targets and the base table's header say 007 and NA
    # (ledger 341). A converter sees the raw text before either reading.
    mapping = pd.read_csv(mapping_csv, converters={"obsname": str, "name": str})
    ucn = flopy.utils.HeadFile(
        str(Path(sim_ws) / f"{model_name}.ucn"), text="concentration"
    )
    try:
        ends = period_end_steps(
            saved_steps(ucn), perlen=perlen, file=f"{model_name}.ucn", saves="CONCENTRATION"
        )
        rows = []
        for per, kstpkper in ends.items():
            data = np.asarray(ucn.get_data(kstpkper=kstpkper), dtype=float)
            row = {"per": int(per)}
            for target in mapping.itertuples(index=False):
                layer_values = np.asarray(data[int(target.layer)], dtype=float).reshape(-1)
                # The mapping's own label: one column per (name, layer). Keyed on
                # the bare name, two screens of one well overwrote each other here
                # while the build read the other one. A mapping written before the
                # label existed has no such column, and its name still serves.
                label = str(getattr(target, "obsname", "") or target.name)
                row[label] = float(layer_values[int(target.cell)])
            rows.append(row)
        pd.DataFrame(rows).to_csv(output_csv, index=False)
    finally:
        ucn.close()


def write_named_series_targets(kind, locations_file, output_csv, sim_ws=".", model_name=None, package=None):
    """Regenerate one lake/SFR/DRN simulated-target CSV in the native forward run.

    Post-model command for named-series observations (lake stage, SFR stage/flow,
    DRN seepage): it reloads the just-run MF6 workspace, rebuilds the target set
    from the saved definition file (the ``*_target_locations`` CSV/GPKG written at
    build time), and writes the simulated series pyEMU reads back as observations.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so it loads the model and reads the saved locations inline rather
    than calling other module-level helpers.

    A lake with no stage is ``NaN`` in ``LakeStageTargets.simulated_series`` and is
    WRITTEN here as 1e30, the value the base run was built with (ledger 329).
    Measured with PEST++ (``docs/reproducers/pest_no_value_fields.py``): an empty
    field and ``nan`` each fail the run, and so does a file one row short of the
    base run, at end of file. A dry lake is otherwise legitimate -- zero-weighted,
    1e30 changes nothing; weighted, it makes that realization's phi ~1e60.

    Each period's value is the one at its END (ledger 336) -- an SFR or DRN flow
    included, which used to be the sum of the steps the period saved -- because
    this calls the same ``simulated_series`` the base run's build did, and that
    reads the one period-end rule on each file's own clock, checked against the
    TDIS it reloads. A realization whose output lacks a period's end fails the
    run rather than hand PEST a value from before the period ended.

    **On the observation route that is one row per stress period too** (ledger
    346): MF6 writes a continuous observation at every time step, and the route
    handed every step on, so under ATS -- which decides the number of steps at run
    time -- a realization wrote another number of rows than the base run (12, 13
    or 14 at K = 10, 1, 2), and PEST, reading each realization by the BASE run's
    line positions, took the value at t = 23 as the end of the last period without
    a word. Every realization writes ``nper`` rows now, whatever its steps. A
    template built on the observation route before that commit expects one line
    per TIME STEP: with more than one step per period (or ATS) it fails at end of
    file on its next run and must be rebuilt; with one step per period the rows
    it reads are the same values under new ``time`` labels, which it skips.

    Parameters
    ----------
    kind : {"lake_stage", "sfr_stage", "sfr_flow", "drn_flow"}
        Which target family the saved definition file describes.
    locations_file : str or pathlib.Path
        The ``*_target_locations`` snapshot the build wrote: a CSV (lake and
        reach ids, or DRN zones with their ``layer`` and comma-joined resolved
        ``cells``) or a GeoPackage (DRN zones drawn as geometry, with the same
        columns). A CSV's names are read as text, verbatim, as the base run's
        table has them (ledger 341); a DRN zone is read in its snapshot's
        ``layer`` (ledger 342 -- the snapshot carried none, so it was read in
        layer 0).
    output_csv : str or pathlib.Path
        The simulated CSV to write, the file pyEMU's instruction file reads.
    sim_ws : str or pathlib.Path, default "."
        The just-run MF6 workspace to reload.
    model_name : str or None, optional
        The CALIBRATED model's name, which the build records. ``None`` (a
        forward run generated before ledger 340) reloads the simulation's first
        GWF model -- which, in a simulation of several, handed PEST another
        model's lake, reach or zone every realization while the base run's
        values came from the calibrated one.
    package : str or None, optional
        The targets' ``package`` -- the name of the LAK, SFR or DRN package they
        belong to -- which the build records. ``None`` (the default, and what a
        forward run generated before ledger 344 passes) reads the model's one
        package of the type; with several, the family refuses as it did at build.

    Raises
    ------
    ValueError
        For an unknown ``kind``; or, as the family's ``simulated_series`` does,
        when a period the output file saved lacks its final step, or its
        observation route refuses what the package declares (ledger 344).
    FileNotFoundError
        As the family's ``simulated_series``: a declared observation file is not
        in ``sim_ws`` -- only the files the package's OBS names are read, so a
        stale copy left in the template is never taken for this run's output.
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
        # Names read VERBATIM, as the base run's table has them (ledger 341): a
        # zone named 007 read as the integer 7, one named NA as a missing value.
        locations = pd.read_csv(locations_path, converters={"name": str})
        if "cells" in locations.columns:
            locations["cells"] = locations["cells"].fillna("").apply(
                lambda text: [int(v) for v in str(text).split(",") if str(v).strip() != ""]
            )

    model = mf.load_mf6_run(Path(sim_ws), model_name=model_name, verbosity_level=0)
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
    targets = builders[kind_key](locations=locations, values=None, package=package)
    simulated = targets.simulated_series(model)
    if kind_key == "lake_stage":
        simulated = simulated.fillna(targets._PEST_NO_STAGE)
    simulated.to_csv(output_csv, index=False)

