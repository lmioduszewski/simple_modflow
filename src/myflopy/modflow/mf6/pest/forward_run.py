"""Forward-run helpers for ``myflopy`` PEST/pyEMU workflows."""

from __future__ import annotations

import numpy as np
import pandas as pd

#: The columns of a build's snapshot files that hold a target's NAME, or a label or
#: column derived from it: the target map's ``name``/``obsname``, a locations or
#: values snapshot's ``name``, the values snapshot's ``column`` (ledger 358).
TARGET_NAME_COLUMNS = ("name", "obsname", "column")

#: The columns of a build's snapshot files read back as TEXT besides the names: a zone
#: snapshot's comma-joined ``cells``, which pandas reads as FLOATS (``27.0``) when
#: every row holds one cell or none (ledger 550).
TARGET_TEXT_COLUMNS = (*TARGET_NAME_COLUMNS, "cells")


def read_target_table(path):
    """Read one of a PEST build's snapshot CSVs with every target name as its TEXT.

    The ONE rule for reading a target name back from a file the build wrote (ledger
    341): a converter hands each field of :data:`TARGET_NAME_COLUMNS` over as the
    text written, before pandas reads it. A plain ``pandas.read_csv`` reads a well
    named ``007`` as the integer 7 and one named ``NA`` as a missing value, so the
    forward-run writers wrote columns ``7`` and ``nan`` (fixed in 341 with this
    converter, inline), and the IES review's residual map -- which read the same
    files plainly -- lost the cell of every such target and dropped such a DRN zone
    entirely, which 341 recorded as not fixed. ``dtype=str`` is not this rule: it keeps ``007`` but
    still reads ``NA`` as missing. A blank name reads as ``""``.

    A zone snapshot's ``cells`` column is read as its text too (ledger
    550): a zone given no ``layer`` is snapshotted as one row per
    layer, each holding its columns ACTIVE there, so a one-column zone over a layer
    where that column is inactive writes ``""`` beside ``"27"``. Pandas read that
    column as floats, ``27.0``, which neither the forward run nor the IES review
    could parse as a cell: every realization failed after a clean build. A blank
    ``cells`` reads as ``""``, a zone of no cells.

    The forward-run writers call this through an import inside their own bodies,
    because pyEMU copies only a writer's source into the generated forward run.

    Parameters
    ----------
    path : str or pathlib.Path
        A CSV the build wrote beside the control file: a target map
        (``<prefix>_head_target_map.csv``), a locations snapshot
        (``<prefix>_target_locations.csv``) or a values snapshot
        (``<prefix>_target_values.csv``).

    Returns
    -------
    pandas.DataFrame
        The table, each of its :data:`TARGET_TEXT_COLUMNS` (the names and a zone's
        ``cells``) as ``str``; every other column as ``pandas.read_csv`` reads it.
    """

    header = pd.read_csv(path, nrows=0).columns
    converters = {column: str for column in TARGET_TEXT_COLUMNS if column in header}
    return pd.read_csv(path, converters=converters)


def apply_pilotpoints_to_array(array_files, manifest="pilot_point_sets.json"):
    """Multiply every pilot-point set on some array files into them (a forward-run pre-command).

    The build adds ONE call per group of array files covered by the same sets --
    one per K layer file, one for every period file of an array RCH (ledger
    547) -- after pyEMU's ``apply_list_and_array_pars``. Each set's
    factor is interpolated ONCE per call: the set's template-populated pilot values,
    combined through the weights the build stored (``factors``: per covered cell, its
    nearest points and their inverse-distance-squared weights) --
    ``10 ** sum(weight * log10(value))`` for a log-transformed parameter,
    ``sum(weight * value)`` otherwise. Each file then starts from itself as pyEMU has
    just rewritten it from its ``org/`` copy when pyEMU parameterizes that file too
    (``base_from_file``), else from the base the build captured for it
    (``base_file``) -- so a worker directory reused across realizations writes the
    same file every run -- and each cell is multiplied by the PRODUCT of the factors
    of every set covering it. A covered cell is then clamped to the tightest of its
    sets' ``physical=`` limits and the file entry's (the pyEMU multipliers' on it,
    ledger 528); a cell no set covers keeps the value it started from. Each file is
    rewritten in the layout it was read in (ledger 525).

    Until 525 this rebuilt the whole layer from the build-time base after
    pyEMU's apply, which erased every constant/zone/grid multiplier on the file and
    every other pilot-point set on it, interpolated over EVERY point of the layer in
    linear space (a field 13-60% above its points' geometric mean), and recomputed
    each cell's distances on every run. A template built before it keeps its own
    copy of that function.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so its imports live inside it.

    Parameters
    ----------
    array_files : str or list of str
        The external array files MF6 reads (``<model>.npf_k_layer1.txt``, or every
        ``<model>.rch_recharge_<n>.txt`` of one period group), relative to the
        working directory, which is the template or a copy of it; one name is one
        file. Every file must be covered by the same sets.
    manifest : str, optional
        The build's pilot-point manifest, ``pilot_point_sets.json`` (the default),
        naming each set's ``array_files``, ``parameter_csv``, ``points_csv``,
        ``factors``, ``interpolation`` and limits and, per file,
        ``base_from_file``, ``base_file`` and the pyEMU multipliers'
        ``lower_limit``/``upper_limit`` (``None`` for none).

    Raises
    ------
    ValueError
        When the manifest has no entry for a file, or the files are not covered by
        the same sets; when a set's value file lacks one of its points' values, or
        holds one that is not a finite number; or when a log-interpolated value is
        not above zero.
    """

    import json

    import pandas as pd

    if isinstance(array_files, str):
        array_files = [array_files]
    with open(manifest, encoding="utf-8") as handle:
        info = json.load(handle)
    entries = {}
    for array_file in array_files:
        found = [entry for entry in info["files"] if entry["array_file"] == array_file]
        if len(found) != 1:
            raise ValueError(f"{manifest} names no pilot-point sets for {array_file!r}.")
        entries[array_file] = found[0]
    names = {tuple(entry["sets"]) for entry in entries.values()}
    if len(names) != 1:
        raise ValueError(
            f"{manifest}: the files {list(array_files)} are not covered by the same "
            "pilot-point sets, so they cannot share one interpolation."
        )
    wanted = set(next(iter(names)))
    starts = {}
    for array_file, entry in entries.items():
        source = array_file if entry["base_from_file"] else entry["base_file"]
        starts[array_file] = np.loadtxt(source, ndmin=2)
    sizes = {start.size for start in starts.values()}
    if len(sizes) != 1:
        raise ValueError(
            f"{list(array_files)} do not hold one value per cell each ({sorted(sizes)} values), "
            "so one interpolation cannot serve them."
        )
    size = sizes.pop()
    factor = np.ones(size)
    lower = np.full(size, -np.inf)
    upper = np.full(size, np.inf)
    covered = np.zeros(size, dtype=bool)
    for record in info["sets"]:
        if record["name"] not in wanted:
            continue
        frame = pd.read_csv(record["parameter_csv"])
        given = pd.Series(
            pd.to_numeric(frame["value"], errors="coerce").to_numpy(),
            index=frame["parnme"].astype(str).str.lower(),
        )
        points = pd.read_csv(record["points_csv"])
        pv = points["parnme"].astype(str).str.lower().map(given).to_numpy(dtype=float)
        if not np.isfinite(pv).all():
            raise ValueError(
                f"{record['parameter_csv']} holds no finite value for some of "
                f"{record['name']!r}'s pilot points."
            )
        stored = np.load(record["factors"])
        cells, index, weight = stored["cells"], stored["index"], stored["weight"]
        if record["interpolation"] == "log":
            if (pv <= 0).any():
                raise ValueError(
                    f"{record['name']!r} is interpolated in log space and a pilot value in "
                    f"{record['parameter_csv']} is not above zero."
                )
            set_factor = 10.0 ** np.sum(weight * np.log10(pv)[index], axis=1)
        else:
            set_factor = np.sum(weight * pv[index], axis=1)
        factor[cells] *= set_factor
        covered[cells] = True
        if record["lower_limit"] is not None:
            lower[cells] = np.maximum(lower[cells], record["lower_limit"])
        if record["upper_limit"] is not None:
            upper[cells] = np.minimum(upper[cells], record["upper_limit"])
    for array_file, entry in entries.items():
        start = starts[array_file]
        # The pyEMU multipliers' physical= on this file: pyEMU clamped its own
        # product, and a pilot factor must not carry a cell back past it (ledger 528).
        file_lower, file_upper = lower, upper
        if entry.get("lower_limit") is not None:
            file_lower = np.maximum(lower, entry["lower_limit"])
        if entry.get("upper_limit") is not None:
            file_upper = np.minimum(upper, entry["upper_limit"])
        result = start.reshape(-1).astype(float) * factor
        result[covered] = np.clip(result[covered], file_lower[covered], file_upper[covered])
        np.savetxt(array_file, result.reshape(start.shape), fmt="%.10E")


def _write_head_target_csv(
    model_name, mapping_csv, output_csv, perlen=None, sim_ws=".", head_file=None, k_sources=None
):
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
    clock -- through the ONE period-end reader ``HeadTargets.compare()`` reads too,
    :func:`myflopy.modflow.mf6.target_values.read_target_values` (ledger
    557). A SCREENED well's rows (several under one ``obsname``) are
    ONE column, their heads weighted by K x saturated open length (ledger
    556), with the K THIS realization's model read: the files
    ``k_sources`` names, which pyEMU's multipliers and the pilot points have just
    rewritten.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so its imports live inside it.

    Parameters
    ----------
    model_name : str
        The flow model's name; its head file is ``<model_name>.hds`` in
        ``sim_ws`` unless ``head_file`` names it.
    mapping_csv : str or pathlib.Path
        The target map the build wrote: ``obsname``/``name``, ``layer``, ``cell``,
        and for a screened well one row per open layer under one ``obsname`` with
        ``open_bottom``, ``open_top``, ``icelltype`` and ``fixed_b``
        (``target_values.SCREEN_COLUMNS``). The order in which labels first appear IS the
        column order written, and its labels are read as text, verbatim (a well
        named ``007`` is not the number 7, ledger 341).
    output_csv : str or pathlib.Path
        The simulated CSV to write: ``per``, then one column per label (its
        ``obsname``, else its ``name``) in mapping order; one row per stress
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
    head_file : str or None, optional
        The head file, relative to ``sim_ws``: the one the model's OC declares
        (``HEAD FILEOUT``), which the build resolves and records (ledger
        505). Default ``None`` -- what a forward run generated before
        then passes -- reads ``<model_name>.hds``, which is wrong whenever the
        OC named the file otherwise: every realization failed, and the build
        could read a STALE ``<model_name>.hds`` left in the run directory.
    k_sources : list of dict or None, optional
        Where this run's horizontal K is, per layer, when the map holds a screened
        well: ``{"layer", "file", "factor"}`` (the NPF's external file for the
        layer, relative to the working directory) or ``{"layer", "constant"}``
        (``target_values.k_from_sources``). The build records it (ledger
        556). Default ``None``: no screened well, so no K read.

    Raises
    ------
    ValueError
        With ``perlen``, when a period the head file saved lacks its final step;
        when the map holds a screened well and ``k_sources`` is ``None``.
    FileNotFoundError
        When the head file, or a K file ``k_sources`` names, is not there.
    """

    # Imported here, not hoisted: pyEMU copies only this function's source.
    from pathlib import Path

    from myflopy.modflow.mf6.pest.forward_run import read_target_table
    from myflopy.modflow.mf6.target_values import k_from_sources, read_target_values

    # Labels read VERBATIM: pandas reads a well named 007 as the integer 7 and one
    # named NA as a missing value, so the columns written would be 7 and nan while
    # the target map, the targets and the base table's header say 007 and NA
    # (ledger 341). ``read_target_table`` is the one rule for that read.
    mapping = read_target_table(mapping_csv)
    # The file the OC declares, which the build recorded (ledger 505).
    head_name = str(head_file) if head_file else f"{model_name}.hds"
    # The K this realization's model READ, after pyEMU and the pilot points
    # rewrote its files (ledger 556); never the build's K.
    k = None if not k_sources else k_from_sources(k_sources, base=sim_ws)
    read_target_values(
        Path(sim_ws) / head_name, mapping, perlen=perlen, file=head_name, saves="HEAD", k=k,
    ).to_csv(output_csv, index=False)


def _write_conc_target_csv(model_name, mapping_csv, output_csv, perlen=None, sim_ws=".", conc_file=None):
    """Write simulated concentrations at each stress period's END for named target locations.

    The concentration twin of ``_write_head_target_csv``, and deliberately the
    same shape: the one period-end reader both use
    (``target_values.read_target_values``, ledger 557) opens the binary
    directly and indexes it, at the record the one period-end rule picks (ledger 336). The named-series helper would instead
    reload the whole simulation on every forward run, and its ``all_conc`` read
    is kind-gated, so it would have to resolve to the GWT model rather than the
    flow model it loads by default. Like it, this is the ONE writer of its file:
    the build writes the base run's copy through it (ledger 340), so the columns
    PEST's positional instruction file names and the ones every realization
    writes are one list in one order.

    Parameters
    ----------
    model_name : str
        The TRANSPORT model's name: its concentration file is
        ``<model_name>.ucn`` in ``sim_ws`` unless ``conc_file`` names it.
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
        The simulation directory: the forward run's working directory, or the
        model workspace when the build writes the base copy.
    conc_file : str or None, optional
        The concentration file, relative to ``sim_ws``: the one the transport
        model's OC declares (``CONCENTRATION FILEOUT``), which the build
        resolves and records (ledger 505) -- in the model's own
        subdirectory when the simulation keeps one per model. Default ``None``
        -- a forward run generated before then -- reads ``<model_name>.ucn``.

    Raises
    ------
    ValueError
        With ``perlen``, when a period the concentration file saved lacks its
        final step.
    FileNotFoundError
        When the concentration file is not in ``sim_ws`` (from FloPy).
    """

    # Imported here, not hoisted: pyEMU copies only this function's source.
    from pathlib import Path

    from myflopy.modflow.mf6.pest.forward_run import read_target_table
    from myflopy.modflow.mf6.target_values import read_target_values

    # Labels read VERBATIM: pandas reads a well named 007 as the integer 7 and one
    # named NA as a missing value, so the columns written would be 7 and nan while
    # the target map, the targets and the base table's header say 007 and NA
    # (ledger 341). ``read_target_table`` is the one rule for that read.
    mapping = read_target_table(mapping_csv)
    # The file the OC declares, which the build recorded (ledger 505).
    conc_name = str(conc_file) if conc_file else f"{model_name}.ucn"
    # The one period-end reader (ledger 557); a concentration
    # map never holds a screened well (refused, ledger 556).
    read_target_values(
        Path(sim_ws) / conc_name, mapping, perlen=perlen, file=conc_name,
        saves="CONCENTRATION", text="concentration",
    ).to_csv(output_csv, index=False)


def write_named_series_targets(sets, sim_ws=".", model_name=None):
    """Regenerate every lake/SFR/DRN/RIV simulated-target CSV in the native forward run.

    Post-model command for named-series observations (lake stage, SFR stage/flow,
    DRN seepage, RIV exchange): it opens the just-run MF6 workspace ONCE, rebuilds
    each target set from the saved definition file (the ``*_target_locations``
    CSV/GPKG written at build time), and writes the simulated series pyEMU reads
    back as observations. It loaded the whole simulation once PER SET, every
    realization -- about 9 s a set on an externalized Ten Trails template beside a
    20 s MF6 run (ledger 540). Now the sets share one file-backed model,
    which loads what the first set needs once.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so it loads the model and reads the saved locations inline rather
    than calling other module-level helpers.

    A DRY lake is its BOTTOM in ``LakeStageTargets.simulated_series`` on both of
    its routes (ledger 351) -- MF6's own STAGE observation of it -- and is written
    as that number. An INACTIVE lake is ``NaN`` there and is WRITTEN here as 1e30,
    the value the base run was built with (ledger 329). Measured with PEST++
    (``docs/reproducers/pest_no_value_fields.py``): an empty field and ``nan`` each
    fail the run, and so does a file one row short of the base run, at end of
    file. Zero-weighted, 1e30 changes nothing; weighted, it makes that
    realization's phi ~1e60 -- which is why a DRY lake is no longer written so:
    PESTPP-IES never recovered its lambda from it
    (``docs/reproducers/dry_lake_bottom.py``).

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
    sets : list of dict
        One entry per named-series set, in the order to write them, each with
        ``kind`` (``"lake_stage"``, ``"sfr_stage"``, ``"sfr_flow"``,
        ``"drn_flow"`` or ``"riv_flow"``: which target family the saved
        definition file describes), ``locations_file`` (the
        ``*_target_locations`` snapshot the build wrote: a CSV -- lake and reach
        ids, or zones with their ``layer`` and comma-joined resolved ``cells`` --
        or a GeoPackage for zones drawn as geometry, with the same columns; a
        CSV's names are read as text, verbatim, as the base run's table has them,
        ledger 341, and a zone is read in its snapshot's ``layer``, ledger 342),
        ``output_csv`` (the simulated CSV to write, the file pyEMU's instruction
        file reads) and ``package`` (the targets' ``package`` -- the name of the
        LAK, SFR, DRN or RIV package they belong to, which the build records;
        ``None`` reads the model's one package of the type, or for a zone family
        every one on the budget route). A missing ``package`` key is ``None``.
    sim_ws : str or pathlib.Path, default "."
        The just-run MF6 workspace to open.
    model_name : str or None, optional
        The CALIBRATED model's name, which the build records. ``None`` opens the
        simulation's first GWF model -- in a simulation of several, possibly
        another model's lake, reach or zone (ledger 340).

    Raises
    ------
    ValueError
        For an unknown ``kind``; or, as the family's ``simulated_series`` does,
        when a period the output file saved lacks its final step, its observation
        route refuses what the package declares (ledger 344), or a zone's package
        wrote no budget record (ledger 536).
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
    from myflopy.modflow.mf6.pest.forward_run import read_target_table

    builders = {
        "lake_stage": mf.LakeStageTargets,
        "sfr_stage": mf.SfrStageTargets,
        "sfr_flow": mf.SfrFlowTargets,
        "drn_flow": mf.DrnFlowTargets,
        "riv_flow": mf.RivFlowTargets,
    }
    for entry in sets:
        if str(entry["kind"]).strip().lower() not in builders:
            raise ValueError(f"Unsupported named-series observation kind {entry['kind']!r}.")

    # ONE file-backed model for every set: it loads the simulation when a set
    # first needs it, and every later set reads what is already loaded.
    model = mf.load_mf6_run(Path(sim_ws), model_name=model_name, verbosity_level=0)
    for entry in sets:
        kind_key = str(entry["kind"]).strip().lower()
        locations_path = Path(entry["locations_file"])
        if locations_path.suffix.lower() == ".gpkg":
            import geopandas as gpd

            locations = gpd.read_file(locations_path)
        else:
            # Names read VERBATIM, as the base run's table has them (ledger 341): a
            # zone named 007 read as the integer 7, one named NA as a missing value.
            locations = read_target_table(locations_path)
            if "cells" in locations.columns:
                locations["cells"] = locations["cells"].fillna("").apply(
                    lambda text: [int(v) for v in str(text).split(",") if str(v).strip() != ""]
                )
        targets = builders[kind_key](locations=locations, values=None, package=entry.get("package"))
        simulated = targets.simulated_series(model)
        if kind_key == "lake_stage":
            simulated = simulated.fillna(targets._PEST_NO_STAGE)
        simulated.to_csv(entry["output_csv"], index=False)


def carry_river_bottom(model_file, pair_csv, stage_col, rbot_col):
    """Move each river record's bottom by its stage's offset, keeping the declared depth.

    Runs as a pre-model command of the forward run AFTER pyEMU's
    ``apply_list_and_array_pars`` (which pyEMU puts first), on a RIV list file a
    ``riv.stage`` parameter offsets. pyEMU moves the STAGE column only; this sets
    ``rbot = rbot_start + (stage - stage_start)`` row by row, from the values the
    build recorded in ``pair_csv``, so a record keeps its depth: lowering a stage
    alone pushed it below its bottom, which MF6 refuses at any icelltype, and changed
    the depth, which moves the river's flux (ledgers 290/291's pair rule; ledger
    519). Idempotent: run twice, it writes the same file.

    Self-contained: pyEMU copies only this function's source into the generated
    forward run, so it must not depend on other module-level helpers.

    Parameters
    ----------
    model_file : str
        The RIV external list file MF6 reads, as pyEMU just wrote it, relative to the
        forward run's directory.
    pair_csv : str
        The build's record of the file's starting ``stage`` and ``rbot``, one row per
        record in file order.
    stage_col : int
        Zero-based column of the stage in ``model_file``, on the model's grid (two
        further right of the cellid's first column on DISV, ledger 420).
    rbot_col : int
        Zero-based column of the river bottom in ``model_file``.

    Raises
    ------
    ValueError
        When ``model_file`` holds a different number of records than ``pair_csv``.
    """

    import pandas as pd

    rows = pd.read_csv(model_file, sep=r"\s+", header=None, dtype=str, keep_default_na=False)
    start = pd.read_csv(pair_csv)
    if len(rows) != len(start):
        raise ValueError(
            f"{model_file} holds {len(rows)} records and {pair_csv} {len(start)}: the "
            "river bottom cannot be carried with its stage record by record."
        )
    stage = pd.to_numeric(rows[stage_col]).to_numpy()
    rbot = start["rbot"].to_numpy(dtype=float) + (stage - start["stage"].to_numpy(dtype=float))
    rows[rbot_col] = [repr(float(value)) for value in rbot]
    with open(model_file, "w", encoding="utf-8") as handle:
        for values in rows.itertuples(index=False):
            handle.write(" ".join(str(value) for value in values) + "\n")
