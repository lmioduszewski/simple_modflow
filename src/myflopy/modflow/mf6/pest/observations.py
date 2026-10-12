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
* **A target's group is its PEST observation group** (ledger 358): a measured
  observation's ``obgnme`` is its target's ``group`` label made a PEST group name
  (:func:`pest_group_name`), else its set's prefix; every other row of a set --
  unmeasured, or a period nobody forecast -- is ``<prefix>_unmeasured``, and a
  forecast is ``<prefix>``, so neither ever joins a user's group. pyEMU gave each
  simulated column a group of its own, so every well was a group of one and
  PEST++'s phi-by-group could not see the user's grouping.

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

from myflopy.modflow.mf6.grid.geometry import model_cell_polygons
from myflopy.modflow.mf6.observations._shared import (
    _is_blank,
    _measured_periods,
    _period_key_hints,
    _shown,
)
from myflopy.modflow.mf6.output_files import _declared_output_path
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
#: (``NaN`` = blank: unmeasured, or a forecast without a placeholder), weight and
#: group label as given (``<NA>``: none; lake and SFR targets carry none).
_TARGET_COLUMNS = ["name", "column", "time", "per", "value", "weight", "group"]

#: An observation name in a PEST instruction line: ``!name!``.
_INSTRUCTION_NAME = re.compile(r"!([^!]+)!")

#: Every run of characters a PEST group name does not keep -- anything but a
#: lower-case letter, a digit, ``_``, ``.``, ``:`` and ``-`` -- becomes one ``_``.
#: PEST++ 5.2.16 refuses the version-2 control file over a comma or a double
#: quote in a group name, and whitespace splits a classic (version-1) one
#: (``docs/reproducers/pest_observation_groups.py``).
_GROUP_NAME_DROPS = re.compile(r"[^a-z0-9_.:\-]+")

#: A group whose name starts with one of these is an INEQUALITY to PESTPP-IES --
#: ``l_``/``less`` "simulated at most", ``g_``/``greater`` "at least" -- matched on
#: the bare start (``lessard`` and ``greaterville`` too): a measured group read so
#: scores 0 wherever the inequality holds (measured, PEST++ 5.2.16, ledger 358).
_INEQUALITY_GROUP_STARTS = ("l_", "less", "g_", "greater")

#: The group of a set's rows no target measures -- unmeasured periods, blank
#: targets, the periods a forecast set does not forecast: ``<prefix>`` + this.
UNMEASURED_GROUP_SUFFIX = "_unmeasured"

#: What ``Series.astype(str)`` writes for a MISSING label -- ``NaN``, ``None``,
#: ``pd.NA`` -- compared lower-cased. Read as a label, it made a group named
#: ``nan`` of every target left unlabelled (ledger 359), so it is refused.
_STRINGIFIED_MISSING = frozenset({"nan", "none", "<na>"})


def pest_group_name(label: str) -> str:
    """The PEST observation group a group label or a set prefix names.

    Lower-cased (pyEMU lower-cases every name it reads), stripped, and every run
    of characters other than letters, digits, ``_``, ``.``, ``:`` and ``-``
    turned into one ``_`` -- so ``"North Wells"`` is ``north_wells``. No length
    limit is applied: PEST++ 5.2.16 read a 5000-character group name from the
    version-2 control file myflopy writes, and pyEMU's own per-column groups
    (``oname:hds_otype:lst_usecol:w1``) are already past classic PEST's 12.

    Parameters
    ----------
    label : str
        The text: a target's ``group`` label, or an observation set's prefix.

    Returns
    -------
    str
        The group name. Two labels may give one name (``"North"``, ``"north"``);
        :func:`assign_observation_groups` refuses that.
    """

    return _GROUP_NAME_DROPS.sub("_", str(label).strip().lower())


def _label_text(label) -> str | None:
    """A target's group label as stripped text, ``None`` when it carries none.

    Missing (``None``, ``NaN``, ``pd.NA``) and blank text are no label -- the
    one rule a blank target VALUE follows too (``_shared._is_blank``).

    Raises
    ------
    TypeError
        When the label is anything else that is not text -- a number, a bool: a
        PEST group is a name, and reading ``1`` as ``"1"`` and ``1.0`` as
        ``"1.0"`` would make one zone two groups (pandas reads an integer column
        holding a blank as floats).
    """

    if _is_blank(label):
        return None
    if isinstance(label, str):
        return label.strip()
    raise TypeError(type(label).__name__)


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
    crs = model_cell_polygons(project.model).crs
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
            # A reading's own weight, when the values table gave one (ledger
            # 504); the zone's and the set's fill in below.
            "weight": (
                values["weight"].astype(float) if "weight" in values.columns else np.nan
            ),
            "group": pd.Series(pd.NA, index=values.index, dtype=object),
        }
    ).reset_index(drop=True)


def _zone_groups(targets, frame: pd.DataFrame, *, prefix: str) -> pd.Series:
    """Each DRN zone's group label, from its location rows; refuse a zone given two.

    A zone is ONE observation per period, so it is in one group: a zone given as
    several rows carrying different labels is refused, as different weights are
    (:func:`_zone_weights`). Rows with no label defer to the ones that have one --
    as ``DrnFlowTargets.compare()`` reads them (``_zone_details``, ledger 359).
    """

    locations = targets.locations_gdf
    if "group" not in locations.columns or frame.empty:
        return frame["group"]
    table = pd.DataFrame(
        {"name": locations["name"].astype(str), "group": locations["group"].astype(object)}
    )
    table = table.loc[~table["group"].map(_is_blank).astype(bool)]
    if table.empty:
        return frame["group"]
    spread = table.groupby("name")["group"].agg(
        lambda labels: len({str(label).strip() for label in labels})
    )
    if (spread > 1).any():
        zones = ", ".join(
            f"{name!r} ({', '.join(sorted(map(repr, set(table.loc[table['name'] == name, 'group']))))})"
            for name in spread[spread > 1].index[:5]
        )
        raise ValueError(
            f"{_caller(type(targets).__name__, prefix)}: zone(s) {zones} are given as "
            "several rows carrying DIFFERENT group labels, and one zone is one "
            "observation per period in one PEST observation group (ledger 358). Give "
            "each zone's rows the same group."
        )
    return frame["name"].map(table.groupby("name")["group"].first()).astype(object)


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
        # The base run's value in each cell of the table, for the check that a
        # WEIGHTED observation starts from a number (ledger 506).
        "base_values": simulated.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float),
    }


def _write_values_snapshot(located: pd.DataFrame, path: Path, *, value_label: str) -> Path:
    """The ``<prefix>_target_values.csv`` review reads: each target's key, value and period.

    It names the simulated column each target reads (``column``: a well screened
    in two layers is two), and once the control file is built each target's
    observation (``obsnme``, ledger 340) and PEST observation group (``obgnme``,
    ledger 358) -- what ``IesResults`` draws one time series per location from.
    """

    snapshot = located.loc[:, ["time", "name", "value", "per"]].rename(
        columns={"value": value_label}
    )
    for column in ("obsnme", "column", "obgnme"):
        if column in located.columns:
            snapshot[column] = located[column].to_numpy()
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
    # A reading's own weight, else its zone's, else the set's, else 1.0 (ledger
    # 504). The set's weight used to overwrite every row, and a values
    # table's weight column was dropped before it got here.
    if metadata_kind == "drn_flow":
        frame["weight"] = frame["weight"].fillna(_zone_weights(targets, frame, prefix=prefix))
        frame["group"] = _zone_groups(targets, frame, prefix=prefix)
    if weight is not None and not frame.empty:
        frame["weight"] = frame["weight"].fillna(float(weight))
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


def _declared_field_file(project, *, model_name: str, kind: str) -> str:
    """The head or concentration file the model's OC declares, relative to the run directory.

    Read by the one rule the run's own readers use (``output_files._declared_output_path``,
    ledger 301). The writers hard-coded ``<model>.hds``/``.ucn``: an OC naming its
    file otherwise failed every realization, and the build read a STALE
    ``<model>.hds`` a previous run had left in the run directory (ledger 505).
    """

    workspace = Path(project.model.workspace)
    return _declared_output_path(workspace, model_name, kind).relative_to(workspace).as_posix()


def _absent_field_file(error, *, owner: str, prefix: str, model_name: str, field_file: dict) -> ValueError:
    """The refusal for a build whose OC-declared head or concentration file is not in the run."""

    (keyword, name), = field_file.items()
    kind, word = {"head_file": ("head", "HEAD"), "conc_file": ("concentration", "CONCENTRATION")}[keyword]
    return ValueError(
        f"{_caller(owner, prefix)}: model {model_name!r} writes its {kind} to {name!r} -- what "
        f"its OC declares ({word} FILEOUT), or <model>.{Path(name).suffix.lstrip('.')} when it "
        f"names none -- and the run directory holds no such file ({error}). Run the model "
        "with that output control before building; the build reads only the file the OC "
        "declares, never another left in the directory (ledger 505)."
    )


def _one_writer_table(
    project,
    spec,
    *,
    writer,
    model_name: str,
    mapping: pd.DataFrame,
    mapping_path: Path,
    simulated_path: Path,
    field_file: dict[str, str],
) -> pd.DataFrame:
    """The base copy of a head or concentration table, written by the forward run's own writer.

    PEST's instruction file reads each realization's copy by POSITION and names what
    it reads after the base copy's header, so the two must list the same columns in
    the same order. Writing both with one function makes that true by construction;
    a ``simulated_values`` table the caller supplies cannot be, so the header check
    after it refuses one in any other order. ``field_file`` is the writer's keyword
    for the OC-declared file (``{"head_file": "renamed.hds"}``), which the forward
    run is given too.
    """

    if spec.simulated_values is None:
        try:
            writer(
                model_name,
                mapping_path,
                simulated_path,
                perlen=period_lengths(project.model),
                sim_ws=Path(project.model.workspace),
                **field_file,
            )
        except FileNotFoundError as error:
            raise _absent_field_file(
                error, owner=type(spec.targets).__name__, prefix=spec.prefix,
                model_name=model_name, field_file=field_file,
            ) from error
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
        .rename(columns={"weight": "row_weight"})
        .assign(per=per.to_numpy())
        .merge(
            labels.loc[:, ["name_key", "name", "obsname", "weight", "group"]],
            on="name_key",
            how="inner",
        )
    )
    # A reading's own weight beats its location's (ledger 504).
    weight = pd.to_numeric(joined["weight"]).astype(float)
    if "row_weight" in joined.columns:
        weight = joined["row_weight"].astype(float).fillna(weight)
    return pd.DataFrame(
        {
            "name": joined["name"].astype(str),
            "column": joined["obsname"].astype(str),
            "time": joined["time"],
            "per": joined["per"].astype(int),
            "value": joined["head_target"].astype(float),
            "weight": weight,
            "group": joined["group"].astype(object),
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
    field_kind = {"head_targets": "head", "conc_targets": "concentration"}[kind]
    field_file = {
        {"head": "head_file", "concentration": "conc_file"}[field_kind]: _declared_field_file(
            project, model_name=model_name, kind=field_kind
        )
    }
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
        field_file=field_file,
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
    item["field_file"] = field_file
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
        # The OC's HEAD FILEOUT, relative to the run (ledger 505).
        "head_file": item.pop("field_file").get("head_file"),
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
        # The transport OC's CONCENTRATION FILEOUT, relative to the run (ledger 505).
        "conc_file": item.pop("field_file").get("conc_file"),
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

    # A DRY lake is its bottom in the series (ledger 351); an INACTIVE one is NaN
    # there and 1e30 in the base run, exactly as the forward run writes it: PEST++
    # reads neither an empty field nor "nan" (ledger 329), and the two files must
    # agree row for row.
    pivot = targets.simulated_series(project.model).fillna(targets._PEST_NO_STAGE)
    # Every DEFINED lake, its number checked and a boundname resolved: the forward
    # run rebuilds the targets from this, and a snapshot of the value rows lost a
    # lake given no value -- every realization wrote a column fewer than the base
    # run, and PEST's positional read failed (flows-3).
    locations_snapshot = targets._definitions_on(project.model)
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
    # Every DEFINED reach (see ``prepare_lake_stage_observations``).
    locations_snapshot = targets._locations_on(project.model)
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="stage_target",
        simulated_wide=pivot,
        metadata_kind="sfr_stage",
        locations_snapshot=locations_snapshot,
        weight=spec.weight,
    )


def prepare_sfr_flow_observations(project, spec: SfrFlowObservationSpec):
    """Create simulated SFR-flow CSVs and register pyEMU observations."""

    targets = spec.targets
    pivot = targets.simulated_series(project.model)
    # Every DEFINED reach (see ``prepare_lake_stage_observations``).
    locations_snapshot = targets._locations_on(project.model)
    return _prepare_named_series_observations(
        project,
        targets=targets,
        prefix=spec.prefix,
        value_column="flow_target",
        simulated_wide=pivot,
        metadata_kind="sfr_flow",
        locations_snapshot=locations_snapshot,
        weight=spec.weight,
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


def _listed(entries, limit: int = 5) -> str:
    """Up to ``limit`` entries joined, and how many more."""

    entries = list(entries)
    more = len(entries) - limit
    return "; ".join(entries[:limit]) + (f"; and {more} more" if more > 0 else "")


def _set_groups(prepared: list[dict], others: set[str]) -> dict[str, str]:
    """Each set's own two groups, ``<prefix>`` and ``<prefix>_unmeasured``; refuse a clash.

    Stamps ``prefix_group`` and ``unmeasured_group`` on each item and returns what
    each such name is (for a refusal naming a label that takes one).
    """

    owned: dict[str, str] = {}
    for item in prepared:
        prefix = str(item["prefix"])
        own = pest_group_name(prefix)
        where = f"observation set {prefix!r}"
        if not re.search(r"[a-z0-9]", own):
            raise ValueError(
                f"The {where}: a set's prefix is the PEST observation group of its "
                "targets that carry no group label (ledger 358), and this one holds no "
                "letter or digit. Pass observe(targets, prefix=...) a name."
            )
        if own.startswith(_INEQUALITY_GROUP_STARTS):
            raise ValueError(
                f"The {where}: a set's prefix is the PEST observation group of its "
                f"targets that carry no group label (ledger 358), and PESTPP-IES reads a "
                f"group whose name starts with {', '.join(map(repr, _INEQUALITY_GROUP_STARTS))} "
                f"as an INEQUALITY, so {own!r} would score 0 wherever the model is on the "
                "allowed side of a target, not when it meets it. Pass a prefix that does "
                "not start so."
            )
        item["prefix_group"] = own
        item["unmeasured_group"] = own + UNMEASURED_GROUP_SUFFIX
        for name, what in (
            (own, f"the group of {where}'s targets with no group label"),
            (own + UNMEASURED_GROUP_SUFFIX, f"the group of {where}'s rows no target measures"),
        ):
            clash = owned.get(name) or (
                f"a group the control file already holds ({name!r})" if name in others else None
            )
            if clash:
                raise ValueError(
                    f"Two observation groups would share the name {name!r}: {what}, and "
                    f"{clash}. Give the sets prefixes that do not become one group name "
                    "(ledger 358)."
                )
            owned[name] = what
    return owned


def _label_groups(prepared: list[dict], owned: dict[str, str], others: set[str]) -> None:
    """Each measured target's group, in ``target_frame["obgnme"]``; refuse a label that cannot be one."""

    labels: dict[str, dict[str, list[str]]] = {}
    not_text: list[str] = []
    stringified: list[str] = []
    for item in prepared:
        frame = item["target_frame"]
        prefix = str(item["prefix"])
        if item.get("is_forecast"):
            frame["obgnme"] = item["prefix_group"]
            continue
        groups = []
        for name, value, label in zip(frame["name"], frame["value"], frame["group"], strict=True):
            if pd.isna(value):
                groups.append(item["unmeasured_group"])
                continue
            try:
                text = _label_text(label)
            except TypeError:
                not_text.append(
                    f"{_shown(name)} in set {prefix!r}: {label!r} ({type(label).__name__})"
                )
                groups.append(None)
                continue
            if text is None:
                groups.append(item["prefix_group"])
                continue
            if text.lower() in _STRINGIFIED_MISSING:
                stringified.append(f"{_shown(name)} in set {prefix!r}: {text!r}")
                groups.append(None)
                continue
            group = pest_group_name(text)
            labels.setdefault(group, {}).setdefault(text, []).append(
                f"{_shown(name)} in set {prefix!r}"
            )
            groups.append(group)
        frame["obgnme"] = groups

    # A text column that keeps a blank BLANK: `astype(str)` writes it 'nan'.
    as_text = (
        "`frame['group'].astype('string')` -- `.astype('Int64').astype('string')` for "
        "whole-number ids, which writes 1.0 as '1'"
    )
    if not_text:
        not_text = list(dict.fromkeys(not_text))
        raise ValueError(
            f"{len(not_text)} target(s) carry a group label that is not text: "
            + _listed(not_text)
            + ". A target's group is its PEST observation group (ledger 358), and a "
            f"group is a NAME: give the column as text that keeps a blank blank ({as_text}; "
            "`astype(str)` writes a blank as 'nan', ledger 359), or labels such as "
            "'zone1', or leave it blank for the set's own group."
        )
    if stringified:
        stringified = list(dict.fromkeys(stringified))
        raise ValueError(
            f"{len(stringified)} target(s) carry the group label 'nan', 'None' or "
            "'<NA>': " + _listed(stringified) + ". That is what `astype(str)` writes for a "
            "MISSING label, so it would make one group of every target left unlabelled "
            f"(ledger 359). Keep a blank blank -- {as_text} -- so those targets are in "
            "their set's own group, or give the group another name."
        )
    nameless = [group_texts for group, group_texts in labels.items() if not re.search(r"[a-z0-9]", group)]
    if nameless:
        raise ValueError(
            "Group label(s) "
            + ", ".join(repr(text) for texts in nameless for text in texts)
            + " hold no letter or digit, so they name no PEST observation group "
            "(ledger 358). Give each a name, or leave it blank for the set's own group."
        )
    inequality = sorted(group for group in labels if group.startswith(_INEQUALITY_GROUP_STARTS))
    if inequality:
        raise ValueError(
            "Group label(s) "
            + ", ".join(repr(text) for group in inequality for text in labels[group])
            + f" become PEST group(s) {', '.join(map(repr, inequality))}, and PESTPP-IES "
            f"reads a group whose name starts with "
            f"{', '.join(map(repr, _INEQUALITY_GROUP_STARTS))} as an INEQUALITY: its "
            "targets would score 0 wherever the model is on the allowed side of them, "
            "not when it meets them (ledger 358). Name the groups otherwise."
        )
    taken = [
        f"{', '.join(map(repr, labels[group]))} becomes {group!r}, which is "
        + (owned.get(group) or "a group the control file already holds")
        for group in sorted(labels)
        if group in owned or group in others
    ]
    if taken:
        raise ValueError(
            "Group label(s) take a name the build gives another group: "
            + _listed(taken)
            + ". A set's targets with no label are grouped under its prefix and its "
            "unmeasured rows under '<prefix>_unmeasured' (ledger 358); name the group "
            "otherwise, or leave those targets unlabeled."
        )
    merged = [
        f"{' and '.join(repr(text) for text in texts)} both become {group!r} "
        f"({'; '.join(where[0] for where in texts.values())})"
        for group, texts in sorted(labels.items())
        if len(texts) > 1
    ]
    if merged:
        raise ValueError(
            "Two different group labels become one PEST observation group: "
            + _listed(merged)
            + ". A group name is lower-cased, and every run of characters other than "
            "letters, digits, '_', '.', ':' and '-' becomes one '_' (ledger 358), so "
            "one of these is very likely a typo of the other. Spell each group one way."
        )


def assign_observation_groups(obs: pd.DataFrame, prepared: list[dict]) -> None:
    """Give every observation of every target set its PEST observation group (ledger 358).

    * A history-matching target WITH a value -- a measured observation -- is in
      its target's ``group`` label, made a PEST group name by
      :func:`pest_group_name`; with no label (missing, or blank text) it is in
      its set's prefix (``hds``, ``lak_stage``, ``sfr_flow``, ...), so phi splits by
      data type out of the box.
    * A forecast target is in its set's prefix (``fore1``, ...); its label is not
      read.
    * Every other row of a set -- a period no target measures, a blank target, a
      period a forecast set does not forecast -- is in ``<prefix>_unmeasured``.

    So a user's group holds exactly the observations measured under its label:
    zero-weight rows and forecasts never join one. pyEMU named each simulated
    column's observations a group of their own (``oname:hds_otype:lst_usecol:w1``),
    so a group was one well and PEST++'s phi by group could not see the user's
    grouping. Observations no set owns (captured parameter fields) keep theirs.

    Parameters
    ----------
    obs : pandas.DataFrame
        The control file's observation data, indexed by name; its ``obgnme`` is
        rewritten for every observation the sets read.
    prepared : list of dict
        The prepared sets, as :func:`finalize_observations` completes them: each
        with ``observation_names`` (every observation its instruction file reads),
        ``forecast_names`` for a forecast, and ``target_frame`` with ``obsnme``.
        Each item gains ``prefix_group`` and ``unmeasured_group``, and its
        ``target_frame`` an ``obgnme`` column.

    Raises
    ------
    ValueError
        When a label is not text (a number, a bool); is ``nan``, ``None`` or
        ``<NA>`` in any case -- what ``astype(str)`` writes for a missing label
        (ledger 359); holds no letter or digit;
        becomes a group PESTPP-IES reads as an inequality (``l_``, ``less``,
        ``g_``, ``greater``); becomes a group the build gives a set's unlabeled
        targets or unmeasured rows, or the control file already holds; or when
        two different labels become one group name -- naming them. Also when a
        set's prefix holds no letter or digit, starts as an inequality, or
        becomes a group another set's does.
    """

    in_sets = {name for item in prepared for name in item["observation_names"]}
    others = set(obs.loc[~obs.index.isin(in_sets), "obgnme"].astype(str).str.lower())
    owned = _set_groups(prepared, others)
    _label_groups(prepared, owned, others)
    for item in prepared:
        obs.loc[item["observation_names"], "obgnme"] = item["unmeasured_group"]
        frame = item["target_frame"]
        if item.get("is_forecast"):
            obs.loc[item["forecast_names"], "obgnme"] = item["prefix_group"]
            continue
        valued = frame.loc[frame["value"].notna()]
        obs.loc[valued["obsnme"], "obgnme"] = valued["obgnme"].to_numpy(dtype=object)


def _refuse_weighted_rows_with_no_base_value(item: dict, valued: pd.DataFrame) -> None:
    """Refuse a positively weighted observation whose BASE run value is not a finite number.

    The base table is where PEST's instruction file reads every realization, and a
    cell the base run left empty -- a DRN zone with no drain in a period on the
    budget route, a supplied ``simulated_values`` table with a gap -- is empty in
    every realization written the same way, so PEST++ fails the run on reading it
    (``docs/reproducers/pest_no_value_fields.py``). The build accepted it and the
    first IES run failed (ledger 506). A weight-0 row is left alone: it is
    never history-matched. MF6's no-value numbers (``|value| >= 1e29``) ARE finite
    and are not refused here.
    """

    base = item.get("base_values")
    if base is None or valued.empty:
        return
    weights = valued["weight"].fillna(1.0).to_numpy(dtype=float)
    rows = valued["row"].to_numpy(dtype=int)
    cols = valued["col"].to_numpy(dtype=int)
    empty = (weights > 0) & ~np.isfinite(base[rows, cols])
    if not empty.any():
        return
    listing = [
        f"{_shown(target.name)} at time {_shown(target.time)} (period {int(target.per)}, "
        f"observation {target.obsnme})"
        for target in valued.loc[empty].itertuples(index=False)
    ]
    raise ValueError(
        f"{_caller(item['owner'], item['prefix'])}: {len(listing)} weighted target(s) read a "
        "cell the base run's simulated table holds no number in, so every realization "
        f"would hand PEST++ an empty field there and fail: {_listed(listing)}. Find why "
        "the model writes nothing there (a zone with no drain in that period, a gap in a "
        "supplied simulated table), or give the target weight 0."
    )


def finalize_observations(project, prepared: list[dict]):
    """Give each set's observations their targets' values, weights and groups; every other row weight 0.

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

    Then every observation of every set gets its group
    (:func:`assign_observation_groups`, ledger 358), and the
    ``<prefix>_target_values.csv`` snapshot is rewritten with each target's
    observation name and group, which is what the review joins on.

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
        observation the control file does not hold; when a positively weighted
        target reads a cell the base run left without a finite number (ledger
        506); or when a group label cannot be a PEST group
        (:func:`assign_observation_groups`).
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
        item["observation_names"] = every
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
            _refuse_weighted_rows_with_no_base_value(item, valued)
    assign_observation_groups(obs, prepared)
    for item in prepared:
        if item.get("values_path") is not None:
            _write_values_snapshot(
                item["target_frame"], item["values_path"], value_label=item["value_label"]
            )
