"""Budget and output-derived table builders for package explorers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from flopy.mf6.mfbase import FlopyException, MFDataException

from myflopy._logging import get_logger
from myflopy.modflow.mf6.hydraulics import unit_labels

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_explorer_utils import (
    _aggregate_hover_strings,
    _coerce_numeric_like_columns,
    _filter_normalized_table,
    # Re-exported: both live one layer down so the LEGACY budget path
    # (``budget_tables``) can share them without dragging four modules up a
    # layer. They are still part of this module's public surface.
    _model_budget_record_frame,
    _normalize_budget_nodes,
    _normalize_connection_type_filter,
    _normalize_surface_water_include,
    _normalize_term_filter,
    budget_term_attribute,
    cell_index_column,
    cell_numbering,
    split_cellid_columns,
)
from myflopy.modflow.mf6.package_registry import (
    _package_type_name,
    model_packages_of_type,
    one_package,
    package_name_of,
    package_outputs,
)
from myflopy.modflow.mf6.package_tables import (
    build_sfr_reach_table,
)
from myflopy.modflow.mf6.period_end import at_last_saved_step, at_last_saved_steps

logger = get_logger(__name__)

#: What "ask a model for its declared units" can raise.
#:
#: The expensive part is not the attribute walk in ``hydraulics.unit_labels`` --
#: it is that on a file-backed model (``LoadedMf6Run``) both ``model.gwf`` and
#: ``model.sim`` are LAZY properties, so reading the units pays for the whole
#: ``MFSimulation.load``. Measured failure modes of that load: flopy's own
#: ``MFDataException``/``FlopyException`` (neither subclasses a builtin), a
#: ``UnicodeDecodeError`` (a ``ValueError``) on a binary/truncated ``mfsim.nam``,
#: ``OSError`` on an unreadable workspace, and ``TypeError``/``ZeroDivisionError``
#: from myflopy's own post-load ``int(sim.tdis.nper.data)`` and cell-count
#: division when the loaded simulation is incomplete.
_UNITS_UNAVAILABLE = (
    AttributeError, TypeError, ValueError, ZeroDivisionError, OSError,
    MFDataException, FlopyException,
)


def budget_value_units(model: SimulationBase) -> str:
    """Return the unit label for a model-budget flow value, by model kind.

    A GWF budget is VOLUMETRIC (L3/T); a GWT budget carries MASS per time and a
    GWE budget ENERGY per time. Every call site used to hardcode ``"ft3/d"``,
    which was wrong twice over: the wrong *dimension* on any transport model,
    and the wrong *units* on any GWF model not declared in feet and days.

    Both halves are read from where MF6 records them (``dis.length_units`` and
    ``tdis.time_units``), through :func:`~myflopy.modflow.mf6.hydraulics.unit_labels`
    -- the reader the map hovers use, so a budget and a hover cannot disagree
    about a unit (ledger 442) -- and fall back to the dimensional
    placeholders ``L`` and ``T`` when a model declares neither (a TDIS saying
    ``unknown`` declares none) -- which transport models built by
    ``mf.gwt``/``mf.gwe`` typically do not. The canonical model declares feet
    and days, so it still reads ``ft3/d``.
    """

    try:
        length_label, time_label = unit_labels(model)
    except _UNITS_UNAVAILABLE:
        # A budget label must never break a table or a map, so the placeholders
        # stay -- but not silently: an entire failed simulation load used to show
        # up only as a label reading "L³/T".
        logger.debug("could not read the declared units", exc_info=True)
        length_label = time_label = None
    time_label = time_label or "T"
    kind = str(getattr(model, "model_type", "") or "").lower()
    if kind.startswith("gwt"):
        return f"M/{time_label}"
    if kind.startswith("gwe"):
        return f"E/{time_label}"
    volume_label = f"{length_label}³" if length_label else "L³"
    return f"{volume_label}/{time_label}"


def build_budget_result_table(
    model: SimulationBase,
    *,
    budget_text: str,
    package_name: str,
    value_name: str = "q",
    packages: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Return a normalized cell-result table from one budget record type.

    EVERY record of the term is read -- one per package per saved time -- each row
    tagged with the package that wrote it (``pname``, MF6's ``paknam2``). It read
    the FIRST record of each time alone, so a second DRN package's flows were
    absent from ``model.drn.results.q`` in silence (ledger 352, measured: ``drn_2``'s
    -199.79 ft3/d on cell 21 missing beside ``drn``'s -207.02 on cell 20).

    Records are matched on the term EXACTLY -- MF6's array-form ``RCHA``/``EVTA``
    records counting as ``RCH``/``EVT`` (``package_registry._package_type_name``).
    FloPy's ``get_data(text=...)`` takes the first record NAME containing the text
    ("RCH" is inside "UZF-GWRCH"), and the ``paknam=`` meant to pin it never matched
    -- MF6 writes the MODEL's name there and the package's in ``paknam2`` -- so every
    read fell back to the text alone, and the answer rested on the file's record
    order. Measured, that order did not bite: MF6 wrote ``RCH`` before ``UZF-GWRCH``
    on a model declaring UZF first. Matching exactly no longer depends on it.

    Parameters
    ----------
    model : SimulationBase
        Live or file-backed model object.
    budget_text : str
        Budget record name such as ``"RCH"`` or ``"UZF-GWRCH"``, any case; matched
        whole against the file's record names, after MF6's array forms are read
        as their type.
    package_name : str
        Logical package name reported in the table's ``package`` column (the
        type: ``"drn"``); it selects nothing.
    value_name : str, default "q"
        Public column name to use for the normalized budget quantity: the
        record's ``q`` renamed -- unless the record carries an AUXILIARY of that
        name (NPF's ``DATA-SAT`` writes the saturation as aux ``sat`` beside a
        placeholder flow of 0.0), which is then the value, and ``q`` is dropped.
    packages : iterable of str or None, optional
        Package NAMES whose records to read, compared with each record's
        ``paknam2`` stripped and upper-cased -- the one package a LAK/SFR/UZF view
        reads. Default ``None`` reads every package's.

    Returns
    -------
    pandas.DataFrame
        Tidy DataFrame with ``model``, ``package``, ``pname`` (the writing package's
        name, upper-cased as MF6 writes it; ``None`` for a record no package wrote,
        such as ``STO-SS``), ``kstpkper``, zero-based ``per``, ``layer``, ``cell``
        and the normalized result field, one row per record entry, in file order.
    """

    reader = model._get_budget_reader()
    headers = reader.recordarray
    wanted = _package_type_name(str(budget_text).strip())
    texts = {raw: _package_type_name(raw.decode().strip()) for raw in set(headers["text"].tolist())}
    selected = np.isin(headers["text"], [raw for raw, text in texts.items() if text == wanted])
    if packages is not None:
        names = {str(name).strip().upper() for name in packages}
        owners = {raw: raw.decode().strip().upper() for raw in set(headers["paknam2"].tolist())}
        selected &= np.isin(headers["paknam2"], [raw for raw, name in owners.items() if name in names])

    frames: list[pd.DataFrame] = []
    for index in np.flatnonzero(selected):
        header = headers[index]
        kstpkper = (int(header["kstp"]) - 1, int(header["kper"]) - 1)
        frame = _model_budget_record_frame(
            model, reader.get_record(int(index)), budget_text=budget_text
        )
        if frame.empty:
            continue
        # ``node`` is already zero-based; only the package-specific second id is
        # left to normalize.
        frame = _normalize_budget_nodes(frame, columns=("node2",))
        nodes = frame["node"].astype(int).tolist()
        layers: list[int] = []
        cells: list[int] = []
        for node in nodes:
            layer, cell = model.node_to_lni[int(node)]
            layers.append(int(layer))
            cells.append(int(cell))
        frame["layer"] = layers
        frame["cell"] = cells
        frame["model"] = model.name
        frame["package"] = str(package_name).lower()
        frame["pname"] = header["paknam2"].decode().strip().upper() or None
        frame["kstpkper"] = [kstpkper] * len(frame)
        frame["per"] = kstpkper[1]
        if value_name != "q" and value_name in frame.columns and "q" in frame.columns:
            # The record carries the value as an AUXILIARY of that name -- NPF's
            # DATA-SAT writes the saturation as aux ``sat`` beside a placeholder
            # flow of 0.0 -- so the aux column is the value and the placeholder
            # goes. Renaming it as below made two ``sat`` columns, and every
            # ``uzf.results.sat`` map and comparison raised (ledger 471).
            frame = frame.drop(columns=["q"])
        elif value_name != "q" and "q" in frame.columns:
            # TRUE rename (drop FloPy's raw "q"): the public column carries the
            # reference frame in its name (q_gwf / q_lake). Renaming rather than
            # duplicating means any consumer that still reads a bare "q" fails
            # loudly instead of silently reading a stale copy.
            frame = frame.rename(columns={"q": value_name})
        frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=["model", "package", "pname", "kstpkper", "per", "layer", "cell", value_name]
        )

    combined = pd.concat(frames, ignore_index=True)
    for column in ("per", "layer", "cell"):
        combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(
        combined,
        exclude={"model", "package", "pname", "kstpkper"},
    )
    return combined


def build_sfr_stage_result_table(
    model: SimulationBase, *, masked: bool = True, package: str | None = None
) -> pd.DataFrame:
    """Return a normalized table of SFR stage results by reach and saved step.

    One row per reach per saved stage record, labelled from the stage file's OWN
    record headers: ``kstpkper`` is the ``(kstp, kper)`` a record names and ``per``
    its stress period. A view that wants one value per period reads each period's
    last saved step (``period_end.at_last_saved_steps``).

    ``stage`` is ``NaN`` wherever MF6 wrote no stage -- a DRY reach (-1e30 in the
    stage file: its depth is zero) or an INACTIVE one (+1e30, ``STATUS
    INACTIVE``) -- because it is read through ``model.outputs.sfr.stage.get()``,
    which masks both by the head file's rule (ledger 353, LAK's ledger 329 for
    SFR). The row stays, so a dry period is a gap in every view built on this
    table, never a missing record. ``masked=False`` keeps the two numbers, which is
    what ``SfrStageTargets`` reads: the user kept a reach TARGET as it was, so a
    dry reach still reaches it as -1e30.

    Parameters
    ----------
    model : SimulationBase
        A run model with an SFR package declaring a ``stage_filerecord``.
    masked : bool, default True
        Passed to ``model.outputs.sfr.stage.get``: ``True`` makes MF6's two
        no-stage sentinels ``NaN``; ``False`` keeps -1e30 (dry) and +1e30
        (inactive) as MF6 wrote them. Every view reads the default.
    package : str or None, optional
        The SFR package's NAME, any case. Default ``None``: the model's only SFR
        package -- several raise ``SeveralPackagesError`` naming
        ``model.sfr(package=...)``, because each package numbers its own reaches
        (ledger 352).

    Returns
    -------
    pandas.DataFrame
        ``model``, ``package`` (``"sfr"``), ``kstpkper`` (zero-based tuple),
        ``per`` (zero-based stress period), ``reach`` (zero-based), ``layer``,
        ``cell``, ``rlen``, ``distance_start``, ``distance_mid``,
        ``distance_end`` (model length units along the network) and ``stage``
        (model length units; ``NaN`` where MF6 wrote none, or its -1e30/+1e30
        under ``masked=False``).

    Raises
    ------
    AttributeError
        If the model has no SFR package or it declares no ``stage_filerecord``.
    SeveralPackagesError
        If the model carries several SFR packages and ``package`` names none.
    ValueError
        When ``package`` names no SFR package.
    """

    one_package(model, "sfr", package)  # several and none named: refused by the view's spelling
    stage_frame = package_outputs(model, "sfr", package).stage.get(masked=masked).copy()
    stage_frame.index.name = "reach"
    stage_long = stage_frame.reset_index().melt(
        id_vars="reach", var_name="kstpkper", value_name="stage"
    )
    stage_long["reach"] = stage_long["reach"].astype(int)
    stage_long["kstpkper"] = stage_long["kstpkper"].apply(
        lambda value: (int(value[0]), int(value[1]))
    )
    stage_long["per"] = stage_long["kstpkper"].apply(lambda value: int(value[1]))

    reach_table = build_sfr_reach_table(model, package=package)
    stage_long = stage_long.merge(reach_table, on="reach", how="left")
    stage_long["model"] = model.name
    stage_long["package"] = "sfr"
    stage_long["layer"] = cell_index_column(stage_long["layer"])
    stage_long["cell"] = cell_index_column(stage_long["cell"])
    ordered = [
        "model",
        "package",
        "kstpkper",
        "per",
        "reach",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        "stage",
    ]
    return stage_long[ordered]


def build_sfr_long_profile_table(
    model: SimulationBase, *, per: int = 0, value_name: str = "q_gwf", package: str | None = None
) -> pd.DataFrame:
    """Return one merged SFR long-profile table for a selected stress period.

    Stage and exchange are each read at the period's LAST saved step on their own
    file's clock -- the stage file follows ``SAVE HEAD``, the SFR budget file ``SAVE
    BUDGET`` (ledgers 301, 337). The exchange ADDED every step the period saved
    until ledger 348 (about 2x at two steps per period under ``SAVE BUDGET ALL``)
    while the stage beside it was the last step's.

    Parameters
    ----------
    model : SimulationBase
        Live or file-backed model object with an SFR package.
    per : int, default 0
        Zero-based stress period to extract, read at its last saved step.
    value_name : str, default "q_gwf"
        Name given to the exchange column (MF6's raw aquifer-referenced ``q``);
        ``"q"`` keeps the raw name.
    package : str or None, optional
        The SFR package's NAME, any case. Default ``None``: the model's only SFR
        package -- several raise ``SeveralPackagesError`` (ledger 352).

    Returns
    -------
    pandas.DataFrame
        Reach-ordered table combining geometry metadata, packagedata fields,
        stage results, and stream-groundwater exchange results. The table is
        intended as the canonical input to higher-level SFR profile plots.

    Raises
    ------
    ValueError
        When the stage file, or an SFR budget file that holds records, saved
        nothing in ``per`` -- naming the periods it has.
    """

    sfr = one_package(model, "sfr", package)
    packagedata = pd.DataFrame(sfr.packagedata.get_data()).copy()
    packagedata = split_cellid_columns(packagedata, cell_numbering(model))
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "reach"})
    packagedata["reach"] = packagedata["reach"].astype(int)
    packagedata["layer"] = cell_index_column(packagedata["layer"])
    packagedata["cell"] = cell_index_column(packagedata["cell"])
    for column in ("rlen", "rwid", "rgrd", "rtp", "rbth", "rhk", "man"):
        if column in packagedata.columns:
            packagedata[column] = pd.to_numeric(packagedata[column], errors="coerce")
    packagedata["streambed_top"] = pd.to_numeric(
        packagedata.get("rtp"), errors="coerce"
    )
    packagedata["streambed_bottom"] = pd.to_numeric(
        packagedata.get("rtp"), errors="coerce"
    ) - pd.to_numeric(packagedata.get("rbth"), errors="coerce")

    reach_table = build_sfr_reach_table(model, package=package)
    packagedata = packagedata.merge(
        reach_table,
        on=["reach", "layer", "cell", "rlen"],
        how="left",
    )

    # Each file's OWN last saved step of the period (ledger 337's rule, ledger 348):
    # the stage file follows SAVE HEAD, the budget file SAVE BUDGET. The exchange
    # ADDED every step the period saved until ledger 348 -- about 2x at the canonical
    # profile's two steps per period -- beside a stage taken from the last step.
    stage = build_sfr_stage_result_table(model, package=package)
    stage = (
        at_last_saved_step(stage, per)
        .loc[:, ["reach", "stage"]]
        .drop_duplicates(subset=["reach"])
        .reset_index(drop=True)
    )

    q_frame = build_sfr_budget_result_table(model, budget_text="SFR", value_name="q", package=package)
    if not q_frame.empty:
        q_frame = (
            at_last_saved_step(q_frame, per)
            .loc[:, ["reach", "q"]]
            .dropna(subset=["reach"])
            .groupby("reach", as_index=False)["q"]
            .sum()
        )
        q_frame["reach"] = q_frame["reach"].astype(int)
    else:
        q_frame = pd.DataFrame(columns=["reach", "q"])

    profile = packagedata.merge(stage, on="reach", how="left").merge(
        q_frame, on="reach", how="left"
    )
    profile["q_per_length"] = np.where(
        pd.to_numeric(profile["rlen"], errors="coerce") > 0.0,
        pd.to_numeric(profile["q"], errors="coerce")
        / pd.to_numeric(profile["rlen"], errors="coerce"),
        np.nan,
    )
    # The exchange column is aquifer-referenced (the "gwf" frame); present it
    # under its frame-explicit public name (q_gwf) while the derived
    # q_per_length keeps its own name.
    if value_name != "q":
        profile = profile.rename(columns={"q": value_name})
    profile["model"] = model.name
    profile["package"] = "sfr"
    profile["per"] = int(per)
    ordered = [
        "model",
        "package",
        "per",
        "reach",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        "rwid",
        "rgrd",
        "rtp",
        "rbth",
        "rhk",
        "man",
        "streambed_top",
        "streambed_bottom",
        "stage",
        value_name,
        "q_per_length",
    ]
    remaining = [column for column in profile.columns if column not in ordered]
    return (
        profile[ordered + remaining]
        .sort_values(["reach", "cell"])
        .reset_index(drop=True)
    )



def build_sfr_budget_result_table(
    model: SimulationBase,
    *,
    budget_text: str = "SFR",
    value_name: str = "q",
    package: str | None = None,
) -> pd.DataFrame:
    """Return the SFR exchange-result table with reach ids attached.

    The ``q`` column carries MF6's RAW sign from the SFR cell-by-cell record:
    flow FROM the reach TO the GWF cell (aquifer-referenced, the "gwf" frame),
    so a NEGATIVE value means the reach gains water from the aquifer. It is not
    normalized -- ``ResultSpec.reference_frame`` declares the frame so plots can
    orient their colours, and the column is presented as ``q_gwf`` downstream.

    Each record is joined to the reach it names (``node2``), one row per record; a
    cell holding several reaches carries one row per reach. The join was on the
    CELL until ledger 470, which gave a cell of ``k`` reaches ``k * k`` rows.

    ONE package's records: ``package`` names it (any case), ``None`` -- the default
    -- is the model's only SFR package, and several raise ``SeveralPackagesError``,
    because each numbers its own reaches (ledger 352). It read the first ``SFR``
    record of each time, whichever package wrote it, and joined it to the reaches
    of the package FloPy's name-first lookup found.
    """

    if package is None and not model_packages_of_type(model, "sfr"):
        names = []  # no SFR package, so no SFR record: the empty table, as before
    else:
        names = [package_name_of(one_package(model, "sfr", package))]
    frame = build_budget_result_table(
        model,
        budget_text=budget_text,
        package_name="sfr",
        value_name=value_name,
        packages=names,
    )
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "kstpkper",
                "per",
                "reach",
                "layer",
                "cell",
                value_name,
            ]
        )

    # Each record names its OWN reach: MF6 writes the GWF cell as ``node`` and the
    # reach as ``node2`` (one-based on file, zero-based here). Joining on the CELL
    # instead gave a cell of k reaches k*k rows, each reach carrying every flow in
    # its cell -- canonical cell 195: reaches 22, 37 and 38 each read the cell's
    # -493.56, where MF6's own SFR budget gives them -20.91, -471.96 and -0.69.
    reach_map = build_sfr_reach_table(model, package=package)
    frame["reach"] = pd.to_numeric(frame["node2"], errors="coerce").astype("Int64")
    frame = frame.merge(
        reach_map.drop(columns=["layer", "cell"]).drop_duplicates(subset=["reach"]),
        on="reach",
        how="left",
    )
    has_per_length = value_name in frame.columns and "rlen" in frame.columns
    if has_per_length:
        # Keep MF6's RAW sign -- do NOT normalize. This is the SFR cell-by-cell
        # record, aquifer-referenced (the "gwf" frame):
        #   negative = the reach GAINS water from the aquifer
        #   positive = the reach LOSES water to the aquifer
        # q_per_length inherits that sign so the table and the map agree. Which
        # sign means "gaining" is declared by ResultSpec.reference_frame, not
        # baked into the data -- see that field for why this is safer than the
        # negation that briefly inverted the SFR map on 2026-07-19. The exchange
        # column is named by ``value_name`` (``q_gwf`` on the public path, ``q``
        # when an internal caller asks for the raw name).
        q_series = pd.to_numeric(frame[value_name], errors="coerce")
        rlen_series = pd.to_numeric(frame["rlen"], errors="coerce")
        frame["q_per_length"] = np.where(
            rlen_series > 0.0, q_series / rlen_series, np.nan
        )
    ordered = [
        "model",
        "package",
        "kstpkper",
        "per",
        "reach",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        value_name,
    ]
    if has_per_length:
        ordered.append("q_per_length")
    remaining = [column for column in frame.columns if column not in ordered]
    return frame[ordered + remaining]


def build_lak_stage_result_table(
    model: SimulationBase, *, masked: bool = True, package: str | None = None
) -> pd.DataFrame:
    """Return a normalized table of LAK stage results mapped to connected cells.

    One row per saved stage record per connected cell, each carrying the lake's one
    stage. **Every row is labelled from the stage file's OWN clock** (ledger 348):
    ``kstpkper`` is the ``(kstp, kper)`` the record's header names
    (``model.outputs.lak.stage.kstpkper()``) and ``per`` is its stress period, the
    columns and order of :func:`build_sfr_stage_result_table`. ``per`` was the
    record's INDEX until then (ledger 301's stage clock), so at ``nstp=3`` under
    ``SAVE HEAD ALL`` the nine records read as periods 0..8 and ``per=1`` selected
    step 2 of period 0. A view that wants one value per period reads each period's
    last saved step (``period_end.at_last_saved_steps``); ``get(per=N)`` keeps every
    step of period ``N``.

    ``stage`` is ``NaN`` wherever MF6 wrote no stage -- a DRY lake (-1e30 in the
    stage file) or an INACTIVE one (+1e30) -- because it is read through
    ``model.outputs.lak.stage.get()``, which masks both by the head file's rule
    (ledger 329). The row stays, so a dry period is a gap in every view built on
    this table, never a missing record. ``masked=False`` keeps the two numbers,
    which is how ``LakeStageTargets`` tells a dry lake (calibrated against its
    bottom, ledger 351) from an inactive one.

    Parameters
    ----------
    model : SimulationBase
        A run model with a LAK package declaring a ``stage_filerecord``.
    masked : bool, default True
        Passed to ``model.outputs.lak.stage.get``: ``True`` makes MF6's two
        no-stage sentinels ``NaN``; ``False`` keeps -1e30 (dry) and +1e30
        (inactive) as MF6 wrote them. Every view reads the default.
    package : str or None, optional
        The LAK package's NAME, any case. Default ``None``: the model's only LAK
        package -- several raise ``SeveralPackagesError`` naming
        ``model.lak(package=...)``, because each package numbers its own lakes
        (ledger 352: FloPy's name-first lookup read the one NAMED ``lak`` beside
        ``lak_2`` as the model's, and the other lake vanished from the table).

    Returns
    -------
    pandas.DataFrame
        ``model``, ``package`` (``"lak"``), ``kstpkper`` (zero-based tuple),
        ``per`` (zero-based stress period), ``lake`` (zero-based), ``layer``,
        ``cell`` and ``stage`` (model length units; ``NaN`` where MF6 wrote none,
        or its -1e30/+1e30 under ``masked=False``).

    Raises
    ------
    AttributeError
        If the model has no LAK package or it declares no ``stage_filerecord``.
    SeveralPackagesError
        If the model carries several LAK packages and ``package`` names none.
    ValueError
        When ``package`` names no LAK package; when the stage file's saved times
        do not line up with its records -- a
        count that differs from the stage array's rows -- naming both: the records
        cannot be labelled, and zipping them would put a stage at the wrong time
        (ledger 301's rule, "a reader whose records do not line up must RAISE").
    """

    lak = one_package(model, "lak", package)
    stage_output = package_outputs(model, "lak", package).stage
    stage_array = np.asarray(stage_output.get(masked=masked), dtype=float)
    saved = stage_output.kstpkper()
    nlakes = int(lak.nlakes.data)
    records = stage_array.reshape(-1, nlakes)
    if len(saved) != records.shape[0]:
        raise ValueError(
            f"The LAK stage file of model {model.name!r} holds {records.shape[0]} "
            f"record(s) of {nlakes} lake(s) but names {len(saved)} saved time(s) in its "
            "record headers, so its records cannot be labelled with the time each was "
            "saved at. Is the file complete -- was it read while MF6 was still writing it?"
        )
    periods = pd.DataFrame(records)
    periods["kstpkper"] = [(int(kstp), int(kper)) for kstp, kper in saved]
    stage_long = periods.melt(id_vars="kstpkper", var_name="lake", value_name="stage")
    stage_long["per"] = stage_long["kstpkper"].map(lambda value: int(value[1]))
    stage_long["lake"] = stage_long["lake"].astype(int)

    connectiondata = pd.DataFrame(lak.connectiondata.get_data()).copy()
    connectiondata = split_cellid_columns(connectiondata, cell_numbering(model))
    if "ifno" in connectiondata.columns:
        connectiondata = connectiondata.rename(columns={"ifno": "lake"})
    cell_map = connectiondata.loc[:, ["lake", "layer", "cell"]].drop_duplicates().copy()
    stage_long = stage_long.merge(cell_map, on="lake", how="left")
    stage_long["model"] = model.name
    stage_long["package"] = "lak"
    stage_long["layer"] = cell_index_column(stage_long["layer"])
    stage_long["cell"] = cell_index_column(stage_long["cell"])
    return stage_long[
        ["model", "package", "kstpkper", "per", "lake", "layer", "cell", "stage"]
    ]


def build_lak_budget_result_table(
    model: SimulationBase,
    *,
    budget_text: str = "GWF",
    value_name: str = "q",
    package: str | None = None,
) -> pd.DataFrame:
    """Return a normalized LAK exchange-result table with connection geometry.

    ONE package's budget file and connections: ``package`` names it (any case),
    ``None`` -- the default -- is the model's only LAK package, and several raise
    ``SeveralPackagesError``, because each numbers its own lakes (ledger 352).

    Notes
    -----
    The returned table carries one row per LAK budget record together with the
    matching LAK connection metadata. Per the MF6 I/O specification, the LAK
    package-output ``GWF`` budget rows do not include the input-file ``iconn``
    identifier. They only carry the lake id, the connected GWF cell/node id, the
    simulated flow ``q``, and the auxiliary ``FLOW-AREA`` value. The helper
    therefore reconstructs connection-specific metadata by matching repeated
    LAK budget rows back to ``lak.connectiondata`` in the written within-cell
    connection order for each stress period. Where the record carries
    ``FLOW-AREA`` the table also has:

    - ``flow_area``: MF6's own ``FLOW-AREA`` for the connection at that step --
      the WETTED area it computed the exchange through, not a geometric constant:
      a vertical connection reads its cell's whole area while the lake is wet over
      it and 0.0 while it is not (13 of the testing canonical lake's 33 vertical
      connections read 0.0 in period 0), and a horizontal one the wetted part of
      its face. The static interface area is
      :func:`~myflopy.modflow.mf6.package_tables.build_lak_connection_table`'s
      ``connection_area``. This said ``flow_area`` was that geometry until ledger
      474; the column itself always held MF6's number.
    - ``q_per_area``: ``q / flow_area`` (length per time), ``NaN`` where
      ``flow_area`` is 0.
    """

    del budget_text
    lak = one_package(model, "lak", package)
    frame = package_outputs(model, "lak", package).bud.get("GWF").copy()
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "kstpkper",
                "per",
                "lake",
                "layer",
                "cell",
                value_name,
            ]
        )
    frame = _normalize_budget_nodes(frame)
    frame["lake"] = pd.to_numeric(frame["node"], errors="coerce").astype("Int64")
    node2 = pd.to_numeric(frame["node2"], errors="coerce")
    layers: list[int] = []
    cells: list[int] = []
    for node in node2.astype(int).tolist():
        layer_id, cell_id = model.node_to_lni[int(node)]
        layers.append(int(layer_id))
        cells.append(int(cell_id))
    frame["layer"] = layers
    frame["cell"] = cells
    frame["model"] = model.name
    frame["package"] = "lak"
    if "kstpkper" in frame.columns:
        frame["kstpkper"] = frame["kstpkper"].apply(
            lambda values: tuple(int(value) for value in values)
        )
        frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
    else:
        frame["kstpkper"] = [(0, 0)] * len(frame)
        frame["per"] = 0
    if value_name != "q" and "q" in frame.columns:
        # TRUE rename: LAK's GWF record is feature-referenced, so the public
        # column names its frame (q_lake). Renaming rather than duplicating means
        # a missed bare-"q" read fails loudly. See ResultSpec.reference_frame.
        frame = frame.rename(columns={"q": value_name})

    connectiondata = pd.DataFrame(lak.connectiondata.get_data()).copy()
    connectiondata = split_cellid_columns(connectiondata, cell_numbering(model))
    if "ifno" in connectiondata.columns:
        connectiondata = connectiondata.rename(columns={"ifno": "lake"})
    if "iconn" in connectiondata.columns:
        connectiondata["iconn"] = pd.to_numeric(
            connectiondata["iconn"], errors="coerce"
        ).astype("Int64")
    connectiondata["lake"] = pd.to_numeric(
        connectiondata["lake"], errors="coerce"
    ).astype("Int64")
    connectiondata["layer"] = pd.to_numeric(
        connectiondata["layer"], errors="coerce"
    ).astype(int)
    connectiondata["cell"] = pd.to_numeric(
        connectiondata["cell"], errors="coerce"
    ).astype(int)
    if "claktype" in connectiondata.columns:
        connectiondata["claktype"] = (
            connectiondata["claktype"].astype("string").str.upper()
        )
    raw_iconn_column = next(
        (column for column in frame.columns if str(column).lower() == "iconn"), None
    )
    use_raw_iconn = False
    if raw_iconn_column is not None and "iconn" in connectiondata.columns:
        if raw_iconn_column != "iconn":
            frame = frame.rename(columns={raw_iconn_column: "iconn"})
        frame["iconn"] = pd.to_numeric(frame["iconn"], errors="coerce").astype("Int64")
        grouped_iconn = frame.groupby(["lake", "layer", "cell"], dropna=False)["iconn"]
        ambiguous_groups = (grouped_iconn.size() > 1) & (
            grouped_iconn.nunique(dropna=False) <= 1
        )
        use_raw_iconn = not bool(ambiguous_groups.any())

    if use_raw_iconn:
        frame = frame.merge(
            connectiondata,
            on=["lake", "iconn"],
            how="left",
            suffixes=("", "_conn"),
        )
        for column in ("layer", "cell"):
            connection_column = f"{column}_conn"
            if connection_column in frame.columns:
                frame[column] = (
                    pd.to_numeric(frame[column], errors="coerce")
                    .astype("Int64")
                    .combine_first(
                        pd.to_numeric(frame[connection_column], errors="coerce").astype(
                            "Int64"
                        )
                    )
                )
                frame = frame.drop(columns=[connection_column])
        if "lake_conn" in frame.columns:
            frame = frame.drop(columns=["lake_conn"])
    else:
        connectiondata["_connection_row_order"] = connectiondata.groupby(
            ["lake", "layer", "cell"]
        ).cumcount()
        frame_group_keys = ["lake", "layer", "cell"]
        if "kstpkper" in frame.columns:
            # LAK package-output budgets repeat the same lake-cell connection rows
            # every time step, so the within-cell row order must reset each period.
            frame_group_keys = ["kstpkper", *frame_group_keys]
        frame["_connection_row_order"] = frame.groupby(frame_group_keys).cumcount()
        frame = frame.merge(
            connectiondata,
            on=["lake", "layer", "cell", "_connection_row_order"],
            how="left",
            suffixes=("", "_conn"),
        )
        frame = frame.drop(columns=["_connection_row_order"], errors="ignore")

    fallback_columns = [
        column
        for column in (
            "iconn",
            "claktype",
            "bedleak",
            "belev",
            "telev",
            "connlen",
            "connwidth",
        )
        if column in connectiondata.columns
    ]
    if fallback_columns:
        fallback_source = connectiondata.loc[
            :, ["lake", "layer", "cell", *fallback_columns]
        ].copy()

        def _first_non_null(series: pd.Series):
            """The first non-null value of a series (NaN if all are null)."""

            values = series.dropna()
            if values.empty:
                return np.nan
            return values.iloc[0]

        fallback_agg: dict[str, object] = {}
        for column in fallback_columns:
            if column == "claktype":
                fallback_agg[column] = _aggregate_hover_strings
            else:
                fallback_agg[column] = _first_non_null
        fallback = (
            fallback_source.groupby(
                ["lake", "layer", "cell"], dropna=False, as_index=False
            )
            .agg(fallback_agg)
            .rename(
                columns={column: f"{column}_fallback" for column in fallback_columns}
            )
        )
        frame = frame.merge(fallback, on=["lake", "layer", "cell"], how="left")
        for column in fallback_columns:
            fallback_column = f"{column}_fallback"
            if column in frame.columns and fallback_column in frame.columns:
                missing_mask = frame[column].isna()
                frame.loc[missing_mask, column] = frame.loc[
                    missing_mask, fallback_column
                ]
                frame = frame.drop(columns=[fallback_column])

    has_per_area = value_name in frame.columns and "FLOW-AREA" in frame.columns
    if has_per_area:
        flow_area = pd.to_numeric(frame.get("FLOW-AREA"), errors="coerce")
        q_series = pd.to_numeric(frame[value_name], errors="coerce")
        frame["flow_area"] = flow_area
        frame["q_per_area"] = np.where(flow_area > 0.0, q_series / flow_area, np.nan)

    ordered = [
        "model",
        "package",
        "kstpkper",
        "per",
        "lake",
        "layer",
        "cell",
        "claktype",
        "belev",
        "telev",
        "connlen",
        "connwidth",
        value_name,
    ]
    if has_per_area:
        ordered.extend(["flow_area", "q_per_area"])
    remaining = [column for column in frame.columns if column not in ordered]
    return frame[ordered + remaining]


#: Every term MODFLOW 6 writes to a LAK or SFR package budget file (``lak.f90`` /
#: ``sfr.f90``'s budget objects). A requested term outside it AND outside the file is a
#: misspelling and is refused; one inside it the file lacks (no mover, no auxiliary) is
#: simply absent, and reads an empty table.
_PACKAGE_BUDGET_TERMS = {
    "lak": frozenset({
        "FLOW-JA-FACE", "GWF", "RAINFALL", "EVAPORATION", "RUNOFF", "EXT-INFLOW",
        "WITHDRAWAL", "EXT-OUTFLOW", "STORAGE", "CONSTANT", "FROM-MVR", "TO-MVR",
        "AUXILIARY",
    }),
    "sfr": frozenset({
        "FLOW-JA-FACE", "GWF", "RAINFALL", "EVAPORATION", "RUNOFF", "EXT-INFLOW",
        "EXT-OUTFLOW", "STORAGE", "FROM-MVR", "TO-MVR", "AUXILIARY",
    }),
}


def _selected_budget_terms(
    package_type: str, available: list[str], requested: list[str] | None
) -> list[str]:
    """The file's terms a ``term=`` filter keeps, refusing a term no budget file has.

    ``budget.get(term="ext_inflow")`` matched nothing and returned an empty table
    without a word (ledger 473): the filter compares MF6's spelling only. A term
    that is neither in the file nor one MF6 writes for the package raises, naming the
    file's terms and, where the attribute spelling was given, the MF6 one.
    """

    if requested is None:
        return list(available)
    known = _PACKAGE_BUDGET_TERMS[package_type] | set(available)
    unknown = [value for value in requested if value not in known]
    if unknown:
        by_attribute = {budget_term_attribute(name): name for name in sorted(known)}
        hints = [
            f"{value.lower()!r} is spelled {by_attribute[budget_term_attribute(value)]!r}"
            for value in unknown
            if budget_term_attribute(value) in by_attribute
        ]
        raise ValueError(
            f"{package_type.upper()} budget term(s) {unknown} are not terms MODFLOW 6 "
            f"writes for this package; the file holds {available}."
            + (f" {'; '.join(hints)}." if hints else "")
        )
    return [value for value in available if value in requested]


def build_lak_budget_term_table(
    model: SimulationBase,
    *,
    term: str | Iterable[str] | None = None,
    package: str | None = None,
) -> pd.DataFrame:
    """Return one canonical long dataframe of LAK package-output budget terms.

    Parameters
    ----------
    model
        Live or file-backed model object exposing ``model.outputs.lak.bud``.
    term
        Optional single term or iterable of terms to include. Terms follow the
        MF6 package-output identifiers such as ``"GWF"``, ``"STORAGE"``, or
        ``"TO-MVR"``.
    package
        The LAK package's NAME, any case. Default ``None``: the model's only LAK
        package -- several raise ``SeveralPackagesError`` (ledger 352).

    Returns
    -------
    pandas.DataFrame
        One normalized long table containing a ``term`` column and common
        metadata columns such as ``model``, ``package``, ``kstpkper``, ``per``,
        ``lake``, and ``q``. Term-specific columns, including ``layer``,
        ``cell``, ``iconn``, ``FLOW-AREA``, ``flow_area``, ``q_per_area``, and
        ``VOLUME``, are preserved when applicable.
    """

    requested_terms = _normalize_term_filter(term)
    one_package(model, "lak", package)  # several and none named: refused by the view's spelling
    budget = package_outputs(model, "lak", package).bud
    available_terms = [str(value).strip().upper() for value in budget.types]
    selected_terms = _selected_budget_terms("lak", available_terms, requested_terms)

    frames: list[pd.DataFrame] = []
    for current_term in selected_terms:
        if current_term == "GWF":
            frame = build_lak_budget_result_table(
                model, budget_text=current_term, value_name="q", package=package
            ).copy()
            if frame.empty:
                continue
            frame["term"] = current_term
            frame["node"] = pd.to_numeric(frame.get("lake"), errors="coerce")
            frames.append(frame)
            continue

        # No isinstance guard: ``bud.get`` returns a DataFrame or raises. It
        # used to warn and hand back a raw record list when the clocks
        # disagreed, and the guard that caught that here CONTINUED -- dropping
        # the term with no signal (measured: STORAGE 38 rows -> 0). Ledger 301.
        frame = budget.get(current_term).copy()
        if frame.empty:
            continue
        frame = _normalize_budget_nodes(frame)
        frame["model"] = model.name
        frame["package"] = "lak"
        frame["term"] = current_term
        if "kstpkper" in frame.columns:
            frame["kstpkper"] = frame["kstpkper"].apply(
                lambda values: tuple(int(value) for value in values)
            )
            frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        else:
            frame["kstpkper"] = [(0, 0)] * len(frame)
            frame["per"] = 0
        if "node" in frame.columns:
            frame["lake"] = pd.to_numeric(frame["node"], errors="coerce").astype(
                "Int64"
            )
        if current_term == "FLOW-JA-FACE" and "node2" in frame.columns:
            frame["lake_to"] = pd.to_numeric(frame["node2"], errors="coerce").astype(
                "Int64"
            )
        frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "term",
                "kstpkper",
                "per",
                "lake",
                "q",
            ]
        )

    combined = pd.concat(frames, ignore_index=True, sort=False)
    ordered = [
        "model",
        "package",
        "term",
        "kstpkper",
        "per",
        "lake",
        "lake_to",
        "layer",
        "cell",
        "iconn",
        "claktype",
        "q",
        "FLOW-AREA",
        "flow_area",
        "q_per_area",
        "VOLUME",
        "node",
        "node2",
    ]
    remaining = [column for column in combined.columns if column not in ordered]
    return combined[
        [column for column in ordered if column in combined.columns] + remaining
    ]


def build_sfr_budget_term_table(
    model: SimulationBase,
    *,
    term: str | Iterable[str] | None = None,
    package: str | None = None,
) -> pd.DataFrame:
    """Return one canonical long dataframe of SFR package-output budget terms.

    ``package`` names the SFR package (any case); ``None``, the default, is the
    model's only one -- several raise ``SeveralPackagesError`` (ledger 352).
    """

    requested_terms = _normalize_term_filter(term)
    one_package(model, "sfr", package)  # several and none named: refused by the view's spelling
    budget = package_outputs(model, "sfr", package).bud
    available_terms = [str(value).strip().upper() for value in budget.types]
    selected_terms = _selected_budget_terms("sfr", available_terms, requested_terms)
    reach_table = build_sfr_reach_table(model, package=package)

    frames: list[pd.DataFrame] = []
    for current_term in selected_terms:
        # No isinstance guard: ``bud.get`` returns a DataFrame or raises. It
        # used to warn and hand back a raw record list when the clocks
        # disagreed, and the guard that caught that here CONTINUED -- dropping
        # the term with no signal (measured: STORAGE 38 rows -> 0). Ledger 301.
        frame = budget.get(current_term).copy()
        if frame.empty:
            continue
        frame = _normalize_budget_nodes(frame)
        frame["model"] = model.name
        frame["package"] = "sfr"
        frame["term"] = current_term
        if "kstpkper" in frame.columns:
            frame["kstpkper"] = frame["kstpkper"].apply(
                lambda values: tuple(int(value) for value in values)
            )
            frame["per"] = frame["kstpkper"].apply(lambda values: int(values[1]))
        else:
            frame["kstpkper"] = [(0, 0)] * len(frame)
            frame["per"] = 0
        if "node" in frame.columns:
            frame["reach"] = pd.to_numeric(frame["node"], errors="coerce").astype(
                "Int64"
            )
        if "reach" in frame.columns:
            frame = frame.merge(
                reach_table, on="reach", how="left", suffixes=("", "_pkg")
            )
        if current_term == "GWF" and "node2" in frame.columns:
            node2_numeric = pd.to_numeric(frame["node2"], errors="coerce")
            layers: list[int | None] = []
            cells: list[int | None] = []
            for node_value in node2_numeric.tolist():
                if pd.isna(node_value):
                    layers.append(None)
                    cells.append(None)
                    continue
                layer_id, cell_id = model.node_to_lni[int(node_value)]
                layers.append(int(layer_id))
                cells.append(int(cell_id))
            frame["layer"] = layers
            frame["cell"] = cells
            q_series = pd.to_numeric(frame.get("q"), errors="coerce")
            rlen_series = pd.to_numeric(frame.get("rlen"), errors="coerce")
            frame["q_per_length"] = np.where(
                rlen_series > 0.0, q_series / rlen_series, np.nan
            )
        elif current_term == "FLOW-JA-FACE" and "node2" in frame.columns:
            frame["reach_to"] = pd.to_numeric(frame["node2"], errors="coerce").astype(
                "Int64"
            )
        else:
            if "q" in frame.columns and "rlen" in frame.columns:
                q_series = pd.to_numeric(frame.get("q"), errors="coerce")
                rlen_series = pd.to_numeric(frame.get("rlen"), errors="coerce")
                frame["q_per_length"] = np.where(
                    rlen_series > 0.0, q_series / rlen_series, np.nan
                )
        frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "term",
                "kstpkper",
                "per",
                "reach",
                "q",
            ]
        )

    combined = pd.concat(frames, ignore_index=True, sort=False)
    ordered = [
        "model",
        "package",
        "term",
        "kstpkper",
        "per",
        "reach",
        "reach_to",
        "layer",
        "cell",
        "rlen",
        "distance_start",
        "distance_mid",
        "distance_end",
        "q",
        "FLOW-AREA",
        "q_per_length",
        "VOLUME",
        "node",
        "node2",
    ]
    remaining = [column for column in combined.columns if column not in ordered]
    return combined[
        [column for column in ordered if column in combined.columns] + remaining
    ]


def build_surface_water_exchange_cell_table(
    model: SimulationBase,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    include: str | Iterable[str] | None = None,
    lak_connection_type: str | Iterable[str] | None = None,
) -> pd.DataFrame:
    """Return a normalized cell table combining SFR and LAK exchange intensity.

    One row per period per cell per package, each period read at its LAST saved
    step on the budget file's own clock (ledgers 337, 348) -- the rows of that one
    step are what add up in a cell. Every step the period saved was ADDED until
    ledger 348.

    Notes
    -----
    The returned ``exchange_intensity`` field uses one unified physical sign
    convention across packages:

    - positive = groundwater gaining into the surface-water feature
    - negative = surface-water losing to groundwater

    LAK already follows that convention in its raw ``q`` sign, so
    ``exchange_intensity = q_per_area``. SFR reports ``GWF`` exchange from the
    stream to groundwater, so its sign is inverted here and
    ``exchange_intensity = -q_per_length``.
    """

    include_packages = _normalize_surface_water_include(include)
    frames: list[pd.DataFrame] = []

    # EVERY LAK and SFR package (ledger 352): the rows are CELLS, so several packages
    # combine as MF6 superposes them -- each package's reaches or lakes joined within
    # it, then the flows and lengths or areas of one cell add. The reach and lake
    # numbers are dropped by the groupby below, so no package's numbering leaks. A
    # type the model does not carry adds no rows: an SFR-only model raised
    # AttributeError on the LAK half, where a LAK-only one read (the SFR half read
    # the budget file, which simply held no SFR record).
    if "sfr" in include_packages:
        sfr_packages = model_packages_of_type(model, "sfr")
        sfr_frame = _concat_packages(
            build_sfr_budget_result_table(
                model, budget_text="SFR", value_name="q",
                package=package_name_of(sfr) if len(sfr_packages) > 1 else None,
            )
            for sfr in sfr_packages
        )
        # One step per period, its last saved (ledger 348): the groupby below is by
        # PERIOD, and it added every step the period saved.
        sfr_frame = at_last_saved_steps(
            _filter_normalized_table(sfr_frame, per=per, layer=layer, cells=None)
        )
        if not sfr_frame.empty:
            grouped = sfr_frame.groupby(
                ["model", "package", "per", "layer", "cell"], as_index=False
            ).agg(q=("q", "sum"), rlen=("rlen", "sum"))
            grouped["q_per_length"] = np.where(
                pd.to_numeric(grouped["rlen"], errors="coerce") > 0.0,
                pd.to_numeric(grouped["q"], errors="coerce")
                / pd.to_numeric(grouped["rlen"], errors="coerce"),
                np.nan,
            )
            # exchange_intensity is myflopy's OWN unified field (positive =
            # feature gains; see this function's docstring), NOT a raw q. SFR's
            # q_per_length is aquifer-referenced (gaining is negative), so negate
            # it here to match LAK, whose feature-referenced q is already
            # positive-means-gaining. The negation lives on this derived field,
            # never on the raw q column, which keeps MF6's sign.
            grouped["exchange_intensity"] = -pd.to_numeric(
                grouped["q_per_length"], errors="coerce"
            )
            grouped["source"] = "sfr"
            frames.append(grouped)

    if "lak" in include_packages:
        lak_packages = model_packages_of_type(model, "lak")
        lak_frame = _concat_packages(
            build_lak_budget_result_table(
                model, budget_text="GWF", value_name="q",
                package=package_name_of(lak) if len(lak_packages) > 1 else None,
            )
            for lak in lak_packages
        )
        lak_frame = at_last_saved_steps(
            _filter_normalized_table(lak_frame, per=per, layer=layer, cells=None)
        )
        connection_types = _normalize_connection_type_filter(lak_connection_type)
        if connection_types is not None and "claktype" in lak_frame.columns:
            lak_frame = lak_frame.loc[
                lak_frame["claktype"]
                .astype("string")
                .str.upper()
                .isin(connection_types)
            ].copy()
        if not lak_frame.empty:
            grouped = lak_frame.groupby(
                ["model", "package", "per", "layer", "cell"], as_index=False
            ).agg(q=("q", "sum"), flow_area=("flow_area", "sum"))
            grouped["q_per_area"] = np.where(
                pd.to_numeric(grouped["flow_area"], errors="coerce") > 0.0,
                pd.to_numeric(grouped["q"], errors="coerce")
                / pd.to_numeric(grouped["flow_area"], errors="coerce"),
                np.nan,
            )
            grouped["exchange_intensity"] = pd.to_numeric(
                grouped["q_per_area"], errors="coerce"
            )
            grouped["source"] = "lak"
            frames.append(grouped)

    if not frames:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "per",
                "layer",
                "cell",
                "source",
                "exchange_intensity",
            ]
        )

    combined = pd.concat(frames, ignore_index=True, sort=False)
    for column in ("per", "layer", "cell"):
        if column in combined.columns:
            combined[column] = pd.to_numeric(combined[column], errors="coerce").astype(
                int
            )
    return combined.reset_index(drop=True)


def _concat_packages(frames) -> pd.DataFrame:
    """One table of several packages' tables (empty when there are none)."""

    frames = [frame for frame in frames if not frame.empty]
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(
        columns=["model", "package", "kstpkper", "per", "layer", "cell"]
    )


def build_lak_stage_change_table(model: SimulationBase, *, package: str | None = None) -> pd.DataFrame:
    """Return per-transition lake-stage changes for one model.

    The resulting table carries one row per lake and stress-period transition,
    with ``stage_change = stage_1 - stage_0``. **Each period is its LAST saved
    step** on the stage file's own clock (``period_end.at_last_saved_steps``, ledger
    337's rule; ledger 348): the table's ``per`` was a saved-record index until
    ledger 348, so at ``nstp=3`` under ``SAVE HEAD ALL`` it reported eight
    "period" transitions between consecutive TIME STEPS of a three-period run. A
    transition into or out of a period with no stage (a dry or inactive lake,
    ``NaN`` in the stage table) is KEPT, with ``stage_change`` ``NaN``: the change
    is undefined, not zero, and dropping the row would make a lake that went dry
    look as if it had no transition at all (ledger 329).

    Parameters
    ----------
    model : SimulationBase
        A run model with a LAK package declaring a ``stage_filerecord``.
    package : str or None, optional
        The LAK package's NAME, any case. Default ``None``: the model's only LAK
        package -- several raise ``SeveralPackagesError`` (ledger 352).

    Returns
    -------
    pandas.DataFrame
        ``model``, ``package``, ``lake``, ``per0``/``per1`` (consecutive stress
        periods the stage file saved, zero-based), ``stage0``/``stage1`` (each
        period's last saved stage) and ``stage_change`` (``stage1 - stage0``, model
        length units).

    Raises
    ------
    ValueError
        As :func:`build_lak_stage_result_table`: a stage file whose saved times do
        not line up with its records.
    """

    stage_table = at_last_saved_steps(build_lak_stage_result_table(model, package=package))
    if stage_table.empty:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "lake",
                "per0",
                "per1",
                "stage0",
                "stage1",
                "stage_change",
            ]
        )

    base = (
        stage_table.sort_values(["lake", "per", "cell"])
        .drop_duplicates(["lake", "per"])
        .loc[:, ["model", "package", "lake", "per", "stage"]]
        .reset_index(drop=True)
    )
    shifted = base.copy()
    shifted["per0"] = shifted["per"]
    shifted["stage0"] = shifted["stage"]
    shifted["per1"] = shifted.groupby("lake")["per0"].shift(-1)
    shifted["stage1"] = shifted.groupby("lake")["stage0"].shift(-1)
    # Only a lake's LAST period has no successor. ``stage1`` is NaN also for a
    # transition INTO a dry period, and dropping on it deleted exactly those rows
    # while keeping the way back out (ledger 329).
    shifted = shifted.dropna(subset=["per1"]).copy()
    shifted["per0"] = shifted["per0"].astype(int)
    shifted["per1"] = shifted["per1"].astype(int)
    shifted["stage_change"] = shifted["stage1"].astype(float) - shifted[
        "stage0"
    ].astype(float)
    return shifted[
        ["model", "package", "lake", "per0", "per1", "stage0", "stage1", "stage_change"]
    ]


__all__ = [
    "_model_budget_record_frame",
    "_normalize_budget_nodes",
    "budget_value_units",
    "build_budget_result_table",
    "build_sfr_stage_result_table",
    "build_sfr_long_profile_table",
    "build_sfr_budget_result_table",
    "build_lak_stage_result_table",
    "build_lak_budget_result_table",
    "build_lak_budget_term_table",
    "build_sfr_budget_term_table",
    "build_surface_water_exchange_cell_table",
    "build_lak_stage_change_table",
]
