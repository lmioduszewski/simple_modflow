"""Observation builders for the ``myflopy`` PEST build: what PEST is handed matches the targets.

Each target set becomes a simulated CSV that pyEMU turns into one observation per
cell of the table, read back on every forward run by an instruction file that finds
each value by POSITION -- the n-th comma of the m-th line. The build holds three
rules about that table (ledger 340):

* **A target is matched to its row by NUMBER and to its column by name**, and its
  observation is whichever one the instruction file reads at that row and column.
  The key is read as a stress period by the one rule every target family uses
  (:func:`~myflopy.modflow.mf6.observations._shared._numeric_periods`), so ``1.0``
  is period 1 in the build as in ``compare()``. The build used to rebuild pyEMU's
  observation NAME from the key's spelling -- ``per:1.0`` where pyEMU wrote
  ``per:1``, ``time:10`` where it wrote ``time:10.0`` -- and skip any name it did
  not find, so those targets reached nothing, silently.
* **A target no row holds is REFUSED**, naming it and the keys the table has.
* **Every row no target names carries weight 0.** pyEMU gives each row weight 1
  and the BASE run's own value, so every unmeasured period -- and every period of a
  forecast set but the ones forecast -- history-matched the calibration back to
  the uncalibrated model (ledger 81's build-side note). Its value stays the base
  run's, so ensembles still report it.

Head and concentration tables are written, base run and every forward run alike,
by ONE writer (``forward_run._write_head_target_csv`` / ``_write_conc_target_csv``),
and the build asserts the base header is the order the forward run writes: the base
copy used to be a pivot, whose columns pandas sorts, so a well named W1 read W10's
head in every realization when names were not already in sorted order.
"""

from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.observations._shared import (
    _measured_periods,
    _period_key_hints,
    _shown,
)
from myflopy.modflow.mf6.period_end import period_lengths
from myflopy.modflow.mf6.pest.forward_run import (
    _write_conc_target_csv,
    _write_head_target_csv,
)
from myflopy.modflow.mf6.pest.model_lookup import resolve_transport_model_name
from myflopy.modflow.mf6.pest.specs import (
    ConcObservationSpec,
    DrnFlowObservationSpec,
    HeadTargetObservationSpec,
    LakeStageObservationSpec,
    SfrFlowObservationSpec,
    SfrStageObservationSpec,
)

#: The columns of a prepared set's ``target_frame``: who the target is, the column
#: of the simulated table it reads, its key as given and as a period, its value
#: (``NaN`` = blank: unmeasured, or a forecast without a placeholder) and weight.
_TARGET_COLUMNS = ["name", "column", "time", "per", "value", "weight"]

#: An observation name in a PEST instruction line: ``!name!``.
_INSTRUCTION_NAME = re.compile(r"!([^!]+)!")


def _write_frame(frame: pd.DataFrame, path: Path):
    """Write a DataFrame to CSV with a predictable location."""

    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    return path


def _write_locations_snapshot(gdf, path: Path):
    """Write a normalized observation-location snapshot to a vector file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    gdf.to_file(path, driver="GPKG")
    return path


def _point_locations_snapshot(project, targets) -> gpd.GeoDataFrame:
    """One target set's locations as a point GeoDataFrame, for the residual map.

    Shared by every POINT-geometry observation family (heads, concentration) so
    the two cannot drift: a family that snapshots its points differently would
    still map, just wrong. Targets built from coordinates already carry a
    geometry; targets built from cell ids do not, so their coordinates come from
    ``match_to_model``, which resolves cell centroids.
    """

    locations_frame = targets.locations_gdf.copy()
    if isinstance(locations_frame, gpd.GeoDataFrame):
        # Only the targets that MATCHED: one left out of the mapping must not be
        # written to the snapshot, or counted in `n_locations`, as if it were used.
        matched = set(targets.match_to_model(project.model)["name"].astype(str))
        return locations_frame.loc[
            locations_frame["name"].astype(str).isin(matched),
            ["name", "layer", "group", "weight", locations_frame.geometry.name],
        ].copy()
    snapshot = targets.match_to_model(project.model).loc[
        :, ["name", "layer", "group", "weight", "x", "y"]
    ].copy()
    crs = getattr(project.model.vor, "crs", None)
    return gpd.GeoDataFrame(
        snapshot.drop(columns=["x", "y"]),
        geometry=gpd.points_from_xy(snapshot["x"], snapshot["y"], crs=crs),
        crs=crs,
    )


def _column_key(label) -> str:
    """How a target's name meets a simulated column's: stripped and lower-cased, as pyEMU names it."""

    return str(label).strip().lower()


def _format_keys(values) -> str:
    """A table's row keys, briefly: whole numbers as integers, the rest to 10 figures."""

    keys = []
    for value in pd.unique(pd.Series(values, dtype=object)):
        number = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if pd.notna(number) and float(number).is_integer():
            keys.append(str(int(number)))
        elif pd.notna(number):
            keys.append(f"{float(number):.10g}")
        else:
            keys.append(repr(value))
    if len(keys) > 12:
        keys = [*keys[:6], "...", *keys[-3:]]
    return ", ".join(keys) if keys else "none"


def _caller(owner: str, prefix: str) -> str:
    """How a refusal names the set being built."""

    return f"{owner} (PEST build, prefix {prefix!r})"


def _named_series_targets(targets, *, value_column: str, prefix: str) -> pd.DataFrame:
    """A lake, reach or zone target set's values, one row per (target, time), in the build's shape.

    Read from the VALUES the set was given -- ``to_long()`` joins them to the
    location rows, so a DRN zone given as two rows repeated each value, and a set
    with no values returned its location table, which has no ``time`` (the
    ``KeyError: 'time'`` a value-less forecast or capture template died on).
    """

    values = targets._values
    if values.empty:
        return pd.DataFrame(columns=_TARGET_COLUMNS)
    per = _measured_periods(
        values, value_column=value_column, caller=_caller(type(targets).__name__, prefix)
    )
    names = values["name"].astype(str)
    return pd.DataFrame(
        {
            "name": names,
            "column": names,
            "time": values["time"],
            "per": per,
            "value": values[value_column].astype(float),
            "weight": np.nan,
        }
    ).reset_index(drop=True)


def _zone_weights(targets, frame: pd.DataFrame, *, prefix: str) -> pd.Series:
    """Each DRN zone's observation weight, from its location rows; refuse a zone given two."""

    locations = targets.locations_gdf
    if "weight" not in locations.columns or frame.empty:
        return frame["weight"]
    table = pd.DataFrame(
        {
            "name": locations["name"].astype(str),
            "weight": pd.to_numeric(locations["weight"], errors="coerce"),
        }
    ).dropna(subset=["weight"])
    spread = table.groupby("name")["weight"].nunique()
    if (spread > 1).any():
        zones = ", ".join(repr(name) for name in spread[spread > 1].index[:5])
        raise ValueError(
            f"{_caller(type(targets).__name__, prefix)}: zone(s) {zones} are given as "
            "several rows carrying DIFFERENT weights, and one zone is one observation "
            "per period with one weight. Give each zone's rows the same weight."
        )
    return frame["name"].map(table.groupby("name")["weight"].first())


def _locate_targets(
    frame: pd.DataFrame,
    simulated: pd.DataFrame,
    *,
    owner: str,
    prefix: str,
    table: str,
    perlen=None,
) -> pd.DataFrame:
    """Each target's row and column in the table PEST reads, found by NUMBER and name.

    The row is the one whose index value (``per``, or ``time``) EQUALS the
    target's period; the column the one whose header is the target's column,
    compared stripped and lower-cased as pyEMU names it.

    Parameters
    ----------
    frame : pandas.DataFrame
        The set's targets (``_TARGET_COLUMNS``).
    simulated : pandas.DataFrame
        The base-run table exactly as written: its first column the row index,
        every other column one observation per row.
    owner : str
        The target class, for the error.
    prefix : str
        The set's observation prefix, for the error.
    table : str
        The table's file name, for the error.
    perlen : sequence of float or None, optional
        The model's period lengths, for the error: a key that is a period's END in
        simulated time -- the observation route's spelling before ledger 346 --
        is named with the period to key it by. ``None``: no such hint.

    Returns
    -------
    pandas.DataFrame
        ``frame`` plus ``row`` (the data row, from 0) and ``col`` (the position
        among the observation columns, from 0).

    Raises
    ------
    ValueError
        When a target matches no row or no column -- naming each such target,
        its key and the table's rows and columns -- or when the table holds two
        rows or two columns a target could mean.
    """

    index_column = simulated.columns[0]
    use_cols = list(simulated.columns[1:])
    keys = pd.to_numeric(simulated[index_column], errors="coerce").to_numpy(dtype=float)
    by_column: dict[str, list[int]] = {}
    for position, label in enumerate(use_cols):
        by_column.setdefault(_column_key(label), []).append(position)

    rows: list[int] = []
    cols: list[int] = []
    missing: list[str] = []
    ambiguous: list[str] = []
    for target in frame.itertuples(index=False):
        found_cols = by_column.get(_column_key(target.column), [])
        found_rows = np.flatnonzero(keys == float(target.per))
        where = f"{_shown(target.name)} at time {_shown(target.time)} (period {int(target.per)})"
        if len(found_cols) > 1 or len(found_rows) > 1:
            ambiguous.append(where)
        elif not found_cols:
            missing.append(f"{where}: no column {target.column!r}")
        elif not len(found_rows):
            missing.append(f"{where}: no row with {index_column} = {int(target.per)}")
        rows.append(int(found_rows[0]) if len(found_rows) == 1 else -1)
        cols.append(found_cols[0] if len(found_cols) == 1 else -1)

    if ambiguous:
        raise ValueError(
            f"{_caller(owner, prefix)}: the simulated table {table} holds two rows with "
            f"one {index_column}, or two columns under one name, so these targets could "
            "each be read from either: " + "; ".join(ambiguous[:5]) + "."
        )
    if missing:
        more = len(missing) - 5
        unplaced = [
            target.time
            for target, row in zip(frame.itertuples(index=False), rows, strict=True)
            if row < 0
        ]
        hints = _period_key_hints(pd.unique(pd.Series(unplaced, dtype=object)), perlen)
        raise ValueError(
            f"{_caller(owner, prefix)}: {len(missing)} target(s) match nothing in the "
            f"table PEST reads, {table}, so PEST would never see them: "
            + "; ".join(missing[:5])
            + (f"; and {more} more" if more > 0 else "")
            + f". Its rows are {index_column} = {_format_keys(simulated[index_column])}; "
            f"its columns are {', '.join(repr(str(c)) for c in use_cols[:12])}"
            + (" ..." if len(use_cols) > 12 else "")
            + ". A target reads the row whose index EQUALS its time read as a zero-based "
            "stress period -- the observation route's too since ledger 346, where MF6's "
            "simulated time was the key. "
            + "".join(f"{hint}. " for hint in hints[:5])
            + "The build used to skip such a target silently and leave the base run's "
            "own value in its place (ledger 340)."
        )
    return frame.assign(row=rows, col=cols).reset_index(drop=True)


def _register_table(project, simulated: pd.DataFrame, simulated_path: Path, *, prefix: str):
    """Hand pyEMU the base copy of a simulated table: one observation per cell of it."""

    project.pf.add_observations(
        simulated_path.name,
        insfile=f"{simulated_path.name}.ins",
        index_cols=simulated.columns[0],
        use_cols=[column for column in simulated.columns[1:]],
        prefix=prefix,
    )


def _prepared_item(
    *,
    prefix: str,
    owner: str,
    located: pd.DataFrame,
    valueless: bool,
    simulated: pd.DataFrame,
    simulated_path: Path,
    values_path: Path | None,
    value_label: str,
) -> dict:
    """The fields every prepared set carries into :func:`finalize_observations`."""

    return {
        "prefix": prefix,
        "owner": owner,
        "target_frame": located,
        "valueless": bool(valueless),
        "instruction_file": f"{simulated_path.name}.ins",
        "shape": (int(len(simulated)), int(len(simulated.columns) - 1)),
        "values_path": values_path,
        "value_label": value_label,
    }


def _write_values_snapshot(located: pd.DataFrame, path: Path, *, value_label: str) -> Path:
    """The ``<prefix>_target_values.csv`` review reads: each target's key, value and period."""

    snapshot = located.loc[:, ["time", "name", "value", "per"]].rename(
        columns={"value": value_label}
    )
    if "obsnme" in located.columns:
        snapshot["obsnme"] = located["obsnme"].to_numpy()
    return _write_frame(snapshot, path)


def _prepare_named_series_observations(
    project,
    *,
    targets,
    prefix: str,
    value_column: str,
    simulated_wide: pd.DataFrame,
    metadata_kind: str,
    locations_snapshot: pd.DataFrame | gpd.GeoDataFrame | None = None,
    geometry: str | None = None,
    weight: float | None = None,
):
    """Create one pyEMU list-style observation source from a named-series target set."""

    owner = type(targets).__name__
    simulated_path = project.template_workspace / f"{prefix}_simulated_{metadata_kind}.csv"
    _write_frame(simulated_wide, simulated_path)
    _register_table(project, simulated_wide, simulated_path, prefix=prefix)

    frame = _named_series_targets(targets, value_column=value_column, prefix=prefix)
    if metadata_kind == "drn_flow":
        frame["weight"] = _zone_weights(targets, frame, prefix=prefix)
    if weight is not None and not frame.empty:
        frame["weight"] = float(weight)
    located = _locate_targets(
        frame, simulated_wide, owner=owner, prefix=prefix, table=simulated_path.name,
        perlen=period_lengths(project.model),
    )
    values_path = None
    if not located.empty:
        values_path = project.template_workspace / f"{prefix}_target_values.csv"
        _write_values_snapshot(located, values_path, value_label=value_column)

    use_cols = list(simulated_wide.columns[1:])
    metadata = {
        "kind": metadata_kind,
        "prefix": prefix,
        "time_column": "time",
        "value_column": value_column,
        "n_locations": int(len(use_cols)),
        "n_rows": int(len(located)),
    }
    if values_path is not None:
        metadata["values_file"] = values_path.name
    locations_path = None
    if locations_snapshot is not None:
        if isinstance(locations_snapshot, gpd.GeoDataFrame):
            locations_path = project.template_workspace / f"{prefix}_target_locations.gpkg"
            _write_locations_snapshot(locations_snapshot, locations_path)
        else:
            locations_path = project.template_workspace / f"{prefix}_target_locations.csv"
            _write_frame(pd.DataFrame(locations_snapshot), locations_path)
        metadata["locations_file"] = locations_path.name
        # Absent for lake/SFR: they snapshot a lake or reach NUMBER, which is
        # not a place on the grid (ledger 76). Only a family that declares its
        # geometry reaches the residual map.
        if geometry is not None:
            metadata["geometry"] = geometry

    return _prepared_item(
        prefix=prefix,
        owner=owner,
        located=located,
        valueless=targets._values.empty,
        simulated=simulated_wide,
        simulated_path=simulated_path,
        values_path=values_path,
        value_label=value_column,
    ) | {
        "named_series_forward_run_config": {
            "kind": metadata_kind,
            "locations_file": locations_path.name if locations_path is not None else None,
            "output_csv": simulated_path.name,
            # The package the targets name (ledger 344), so every realization
            # reads the observations the base run did.
            "package": getattr(targets, "package", None),
        },
        "metadata": metadata,
    }


def _one_writer_table(
    project,
    spec,
    *,
    writer,
    model_name: str,
    mapping: pd.DataFrame,
    mapping_path: Path,
    simulated_path: Path,
) -> pd.DataFrame:
    """The base copy of a head or concentration table, written by the forward run's own writer.

    PEST's instruction file reads each realization's copy by POSITION and names what
    it reads after the base copy's header, so the two must list the same columns in
    the same order. Writing both with one function makes that true by construction;
    a ``simulated_values`` table the caller supplies cannot be, so the header check
    after it refuses one in any other order.
    """

    if spec.simulated_values is None:
        writer(
            model_name,
            mapping_path,
            simulated_path,
            perlen=period_lengths(project.model),
            sim_ws=Path(project.model.workspace),
        )
        simulated = pd.read_csv(simulated_path)
    else:
        simulated = (
            spec.simulated_values.copy()
            if isinstance(spec.simulated_values, pd.DataFrame)
            else pd.read_csv(spec.simulated_values)
        )
        _write_frame(simulated, simulated_path)

    labels = mapping["obsname"] if "obsname" in mapping.columns else mapping["name"]
    expected = ["per", *[str(label) for label in labels]]
    header = [str(column) for column in simulated.columns]
    if header != expected:
        raise ValueError(
            f"{_caller(type(spec.targets).__name__, spec.prefix)}: the base simulated "
            f"table {simulated_path.name} has columns {header}, and every forward run "
            f"writes {expected}. PEST reads that file by POSITION and names each value "
            "after the base table's header, so a different order hands every "
            "realization one location's value under another's name. Give "
            "`simulated_values` the forward run's columns in its order, or leave it "
            "None to have the build write it."
        )
    return simulated


def _point_target_frame(targets, labels: pd.DataFrame, *, prefix: str) -> pd.DataFrame:
    """A head or concentration target set's values, one row per (screen, time), in the build's shape.

    Joined from the VALUES table to the matched locations' labels, so a name
    screened in two layers is two targets -- one per label.

    Raises
    ------
    ValueError
        When a values row names a location that lands on no active cell of the
        model (ledger 341). ``match_to_model`` only WARNS about such a partial miss
        -- one target file may serve several grids (ledger 211) -- and that warning
        goes to a logger nobody reads, so the inner join here dropped the target
        from the control file in silence: the one route left by which a target
        with a value never reached PEST.
    """

    values = targets._values
    caller = _caller(type(targets).__name__, prefix)
    per = _measured_periods(
        values,
        value_column="head_target",
        name_column="name_key",
        caller=caller,
    )
    left_out = values.loc[~values["name_key"].isin(set(labels["name_key"]))]
    if not left_out.empty:
        listing = [
            f"{_shown(name)} at time(s) {_format_keys(group['time'])}"
            for name, group in left_out.groupby("name", sort=True)
        ]
        more = len(listing) - 5
        raise ValueError(
            f"{caller}: {len(listing)} target(s) name a location that lands on no "
            "active cell of this model, so PEST would never see them: "
            + "; ".join(listing[:5])
            + (f"; and {more} more" if more > 0 else "")
            + ". Comparing a target file against a model leaves such a location "
            "out with a warning, because one file may serve several grids (ledger "
            "211); a PEST build is one model, and each of these is a target the "
            "calibration would lose without a word (ledger 341). Build the targets "
            "from the locations inside the model's active cells."
        )
    joined = (
        values.drop(columns=["name"], errors="ignore")
        .assign(per=per.to_numpy())
        .merge(labels.loc[:, ["name_key", "name", "obsname", "weight"]], on="name_key", how="inner")
    )
    return pd.DataFrame(
        {
            "name": joined["name"].astype(str),
            "column": joined["obsname"].astype(str),
            "time": joined["time"],
            "per": joined["per"].astype(int),
            "value": joined["head_target"].astype(float),
            "weight": pd.to_numeric(joined["weight"]).astype(float),
        }
    ).reset_index(drop=True)


def _prepare_point_observations(
    project,
    spec,
    *,
    writer,
    model_name: str,
    kind: str,
    value_label: str,
    mapping_name: str,
    simulated_name: str,
):
    """Heads or concentrations: map, base table (one writer), locations, targets."""

    owner = type(spec.targets).__name__
    labels = spec.targets.observation_labels(project.model)
    mapping = labels.loc[:, ["obsname", "name", "layer", "cell"]]
    mapping_path = project.template_workspace / f"{spec.prefix}_{mapping_name}"
    _write_frame(mapping, mapping_path)
    simulated_path = project.template_workspace / f"{spec.prefix}_{simulated_name}"
    simulated = _one_writer_table(
        project,
        spec,
        writer=writer,
        model_name=model_name,
        mapping=mapping,
        mapping_path=mapping_path,
        simulated_path=simulated_path,
    )
    _register_table(project, simulated, simulated_path, prefix=spec.prefix)

    locations_snapshot = _point_locations_snapshot(project, spec.targets)
    locations_path = project.template_workspace / f"{spec.prefix}_target_locations.gpkg"
    _write_locations_snapshot(locations_snapshot, locations_path)

    frame = _point_target_frame(spec.targets, labels, prefix=spec.prefix)
    located = _locate_targets(
        frame, simulated, owner=owner, prefix=spec.prefix, table=simulated_path.name,
        perlen=period_lengths(project.model),
    )
    values_path = project.template_workspace / f"{spec.prefix}_target_values.csv"
    _write_values_snapshot(located, values_path, value_label=value_label)

    item = _prepared_item(
        prefix=spec.prefix,
        owner=owner,
        located=located,
        valueless=spec.targets._values.empty,
        simulated=simulated,
        simulated_path=simulated_path,
        values_path=values_path,
        value_label=value_label,
    )
    item["metadata"] = {
        "kind": kind,
        "prefix": spec.prefix,
        "geometry": "points",
        "locations_file": locations_path.name,
        "mapping_file": mapping_path.name,
        "values_file": values_path.name,
        "name_column": "name",
        "layer_column": "layer",
        "group_column": "group",
        "weight_column": "weight",
        "time_column": "time",
        "value_column": value_label,
        "n_locations": int(mapping["name"].nunique()),
        "n_rows": int(len(located)),
    }
    item["mapping_csv"] = mapping_path.name
    item["output_csv"] = simulated_path.name
    return item


def prepare_head_target_observations(project, spec: HeadTargetObservationSpec):
    """Create simulated-head CSVs and register pyEMU observations.

    The base table is written by the forward run's own writer
    (``forward_run._write_head_target_csv``) from the model's head file, so its
    columns are the forward run's in the forward run's order (ledger 340).

    Parameters
    ----------
    project : PestProject
        The calibration being built (its ``model``, ``pf`` and
        ``template_workspace``).
    spec : HeadTargetObservationSpec
        The head targets and their prefix.

    Returns
    -------
    dict
        The prepared set :func:`finalize_observations` completes.

    Raises
    ------
    ValueError
        When a target's key is not a whole number, two values name one target
        in one period, a target's location lands on no active cell (ledger 341),
        or a target matches no row or column of the table; or when a
        ``simulated_values`` table's columns are not the forward run's, in its
        order.
    """

    item = _prepare_point_observations(
        project,
        spec,
        writer=_write_head_target_csv,
        model_name=project.model.name,
        kind="head_targets",
        value_label="head",
        mapping_name="head_target_map.csv",
        simulated_name="simulated_heads.csv",
    )
    item["forward_run_config"] = {
        "mapping_csv": item.pop("mapping_csv"),
        "output_csv": item.pop("output_csv"),
    }
    return item


def prepare_conc_observations(project, spec: ConcObservationSpec):
    """Create simulated-concentration CSVs and register pyEMU observations.

    The concentration twin of :func:`prepare_head_target_observations`: the base
    table is written by ``forward_run._write_conc_target_csv`` from the TRANSPORT
    model's ``.ucn`` beside the flow model's head file.

    Parameters
    ----------
    project : PestProject
        The calibration being built; its ``model`` is the FLOW model.
    spec : ConcObservationSpec
        The concentration targets, their prefix and, optionally, the GWT model's
        name (found from the simulation when ``None``).

    Returns
    -------
    dict
        The prepared set :func:`finalize_observations` completes.

    Raises
    ------
    ValueError
        As :func:`prepare_head_target_observations`; or when the simulation holds
        no GWT model, or several and none is named.
    """

    transport_name = spec.transport_model_name or resolve_transport_model_name(project)
    item = _prepare_point_observations(
        project,
        spec,
        writer=_write_conc_target_csv,
        model_name=transport_name,
        kind="conc_targets",
        value_label="conc",
        mapping_name="conc_target_map.csv",
        simulated_name="simulated_conc.csv",
    )
    item["metadata"]["transport_model"] = transport_name
    item["conc_forward_run_config"] = {
        "model_name": transport_name,
        "mapping_csv": item.pop("mapping_csv"),
        "output_csv": item.pop("output_csv"),
    }
    return item


def prepare_lake_stage_observations(project, spec: LakeStageObservationSpec):
    """Create simulated lake-stage CSVs and register pyEMU observations."""

    targets = spec.targets
    if targets is None:
        raise ValueError(
            "LakeStageObservationSpec now expects targets=LakeStageTargets(...) "
            "as the canonical workflow."
        )

    # A lake with no stage is NaN in the series and 1e30 in the base run, exactly
    # as the forward run writes it: PEST++ reads neither an empty field nor "nan"
    # (ledger 329), and the two files must agree row for row.
    pivot = targets.simulated_series(project.model).fillna(targets._PEST_NO_STAGE)
    locations_snapshot = targets.get().loc[:, ["name", "lake"]].drop_duplicates()
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="stage",
        simulated_wide=pivot,
        metadata_kind="lake_stage",
        locations_snapshot=locations_snapshot,
        weight=spec.weight,
    )


def prepare_sfr_stage_observations(project, spec: SfrStageObservationSpec):
    """Create simulated SFR-stage CSVs and register pyEMU observations."""

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.get().loc[:, ["name", "reach"]].drop_duplicates()
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="stage_target",
        simulated_wide=pivot,
        metadata_kind="sfr_stage",
        locations_snapshot=locations_snapshot,
    )


def prepare_sfr_flow_observations(project, spec: SfrFlowObservationSpec):
    """Create simulated SFR-flow CSVs and register pyEMU observations."""

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.get().loc[:, ["name", "reach"]].drop_duplicates()
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="flow_target",
        simulated_wide=pivot,
        metadata_kind="sfr_flow",
        locations_snapshot=locations_snapshot,
    )


def prepare_drn_flow_observations(project, spec: DrnFlowObservationSpec):
    """Create simulated DRN-zone seepage CSVs and register pyEMU observations.

    The zone snapshot the forward run rebuilds the targets from carries each row's
    RESOLVED ``cells`` and its ``layer``, so every realization sums the drains the
    base run did. It dropped ``layer`` (ledger 342), and the forward run then read
    each zone in layer 0.
    """

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    locations_snapshot = targets.zone_definitions(project.model)
    columns = ["name", "group", "weight", "layer", "cells"]
    if isinstance(locations_snapshot, gpd.GeoDataFrame):
        snapshot = locations_snapshot.loc[:, [*columns, locations_snapshot.geometry.name]].copy()
    else:
        snapshot = locations_snapshot.loc[:, [column for column in columns if column in locations_snapshot.columns]].copy()
    if "cells" in snapshot.columns:
        snapshot["cells"] = snapshot["cells"].apply(lambda values: ",".join(str(value) for value in values))
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="flow_target",
        simulated_wide=pivot,
        metadata_kind="drn_flow",
        locations_snapshot=snapshot,
        geometry="zones",
    )


def _instruction_names(path: Path, *, shape: tuple[int, int]) -> list[list[str]]:
    """The observation each cell of a simulated table is read into, by POSITION.

    pyEMU writes one instruction line per line of the CSV -- a bare ``l1`` for the
    header, then one line per data row naming, in column order, the observation
    read at each comma -- and PEST++ reads the file exactly so. That is the
    ground truth of which observation holds which number; the observation's NAME
    is pyEMU's business and is never rebuilt here.

    Parameters
    ----------
    path : pathlib.Path
        The instruction file pyEMU wrote beside the table.
    shape : tuple of (int, int)
        The table's data rows and observation columns.

    Returns
    -------
    list of list of str
        ``names[row][col]``, lower-cased as the control file holds them.

    Raises
    ------
    ValueError
        When the file's lines do not line up with the table: another number of
        data lines, a header line reading an observation, or a line reading
        another number of columns.
    """

    lines = path.read_text().splitlines()
    rows = [
        [name.strip().lower() for name in _INSTRUCTION_NAME.findall(line)]
        for line in lines[1:]
    ]
    nrows, ncols = shape
    if (
        not lines
        or not lines[0].lower().startswith("pif")
        or len(rows) != nrows + 1
        or rows[0]
        or any(len(row) != ncols for row in rows[1:])
    ):
        read = [len(row) for row in rows]
        raise ValueError(
            f"The instruction file {path.name} does not line up with the table it "
            f"reads ({nrows} data rows of {ncols} observations after one header "
            f"line); it reads {read[:8]}{' ...' if len(read) > 8 else ''} observations "
            "per line. Which observation holds which value cannot be told, so none "
            "is assigned."
        )
    return rows[1:]


def finalize_observations(project, prepared: list[dict]):
    """Give each set's observations their targets' values and weights, and every other row weight 0.

    Runs after the control file is built. For each prepared set it reads which
    observation the instruction file takes from each cell of the set's table
    (:func:`_instruction_names`), then:

    * every observation of the set gets weight 0 -- a row no target names keeps
      the base run's value as its ``obsval`` but is never history-matched;
    * a history-matching target with a value gets that value and its weight (the
      target's own, else 1.0); a blank one is unmeasured and stays at weight 0;
    * a forecast target stays at weight 0, takes its placeholder value when it has
      one, and is recorded in ``item["forecast_names"]`` for
      ``PestProject._apply_forecasts``. A forecast set given NO values names every
      row of its table -- its locations at every period, the only reading under
      which declaring it forecasts anything.

    The ``<prefix>_target_values.csv`` snapshot is rewritten with each target's
    observation name, which is what the review joins on.

    Parameters
    ----------
    project : PestProject
        The calibration being built, its ``pst`` built.
    prepared : list of dict
        The prepared sets (``_prepare_observation_specs``), each carrying
        ``is_forecast`` when it came from ``forecast()``.

    Raises
    ------
    ValueError
        When an instruction file does not line up with its table, or names an
        observation the control file does not hold.
    """

    obs = project.pst.observation_data
    for item in prepared:
        names = _instruction_names(
            project.template_workspace / item["instruction_file"], shape=item["shape"]
        )
        every = [name for row in names for name in row]
        absent = [name for name in every if name not in obs.index]
        if absent:
            raise ValueError(
                f"Observation set {item['prefix']!r}: its instruction file reads "
                f"{len(absent)} observation(s) the control file does not hold "
                f"({', '.join(absent[:3])}), so the two were not built together."
            )
        obs.loc[every, "weight"] = 0.0

        frame = item["target_frame"]
        frame["obsnme"] = [names[row][col] for row, col in zip(frame["row"], frame["col"], strict=True)]
        valued = frame[frame["value"].notna()]
        if item.get("is_forecast"):
            obs.loc[valued["obsnme"], "obsval"] = valued["value"].to_numpy(dtype=float)
            named = every if item.get("valueless") else frame["obsnme"].tolist()
            item["forecast_names"] = list(dict.fromkeys(named))
        else:
            if valued["obsnme"].duplicated().any():
                twice = sorted(set(valued.loc[valued["obsnme"].duplicated(), "obsnme"]))
                raise ValueError(
                    f"Observation set {item['prefix']!r}: two targets carry a value for "
                    f"one observation ({', '.join(twice[:3])})."
                )
            obs.loc[valued["obsnme"], "obsval"] = valued["value"].to_numpy(dtype=float)
            obs.loc[valued["obsnme"], "weight"] = valued["weight"].fillna(1.0).to_numpy(dtype=float)
        if item.get("values_path") is not None:
            _write_values_snapshot(frame, item["values_path"], value_label=item["value_label"])
