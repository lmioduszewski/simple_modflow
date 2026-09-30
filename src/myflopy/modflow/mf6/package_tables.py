"""Normalized package input and static table builders."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from flopy.mf6.mfbase import MFDataException

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_explorer_utils import (
    _coerce_numeric_like_columns,
    _extract_structured_column,
    _get_package_explorer_cache,
    _normalize_iterable_filter,
    split_cellid_columns,
)
from myflopy.modflow.mf6.package_registry import (
    model_packages_of_type,
    one_package,
    package_name_of,
    package_outputs,
)

#: What "build a normalized table for this package on this model" can raise
#: when the package is absent, empty, or in a form this tier cannot read.
#:
#: Callers that tolerate a missing package -- the diff tiers, the group
#: accessors -- catch this rather than `Exception`, so a genuine bug in the
#: builders still surfaces. Each member is a measured path:
#:
#: * ``AttributeError`` -- ``model.package(name)`` raises it for a package that
#:   is not attached, and an ARRAY-form package (``ModflowGwfrcha``) has no
#:   ``stress_period_data`` attribute at all.
#: * ``KeyError`` -- flopy's ``get_data()`` returns None for empty storage, and
#:   ``pd.DataFrame(None)`` is a zero-COLUMN frame, so the next column access
#:   raises rather than yielding an empty table.
#: * ``TypeError`` -- on DISU flopy returns ``cellid`` as a 1-TUPLE, which
#:   ``split_cellid_columns`` passes to ``int()``.
#: * ``ValueError`` -- pandas integer casting on a column carrying NaN.
#: * ``OSError`` / ``MFDataException`` -- a file-backed model resolving its
#:   packages performs a lazy ``MFSimulation.load``.
PACKAGE_TABLE_UNAVAILABLE = (
    AttributeError, KeyError, TypeError, ValueError, OSError, MFDataException,
)


def build_cell_package_input_table(
    model: SimulationBase,
    package_name: str,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
) -> pd.DataFrame:
    """Return a normalized table for a cell-based MF6 stress-period package.

    Parameters
    ----------
    model
        Live or file-backed model object.
    package_name
        MF6 package name such as ``"rch"``, ``"chd"``, ``"drn"``, or ``"ghb"``.

    Returns
    -------
    pandas.DataFrame
        Tidy DataFrame with at least ``package``, ``per``, ``layer``, and
        ``cell`` columns plus the original package data fields.
    """

    layer_values = _normalize_iterable_filter(layer)
    cell_values = _normalize_iterable_filter(cells)
    if per is None and layer_values is None and cell_values is None:
        cache_key = ("cell_package_input_table", str(package_name).lower())
        cached = _get_package_explorer_cache(model).get(cache_key)
        if cached is not None:
            return cached.copy()
    # EVERY package of the type -- two CHDs, a list RCH beside an array one -- each
    # row tagged with its package's own name. Two same-type packages used to give 0
    # rows (FloPy handed back a list, read as an array package with no fields) or
    # the first package only, and an array RCH was never asked for (ledger 229).
    frames = [
        _one_package_table(
            model, package, package_name, per=per, layer_values=layer_values, cell_values=cell_values
        )
        for package in _packages_of(model, package_name)
    ]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame(columns=["model", "package", "pname", "per", "layer", "cell"])

    combined = pd.concat(frames, ignore_index=True)
    for column in ("layer", "cell", "per"):
        if column in combined.columns:
            combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(combined, exclude={"model", "package", "pname"})
    if per is None and layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[
            ("cell_package_input_table", str(package_name).lower())
        ] = combined.copy()
    return combined


def _packages_of(model, package_name) -> list:
    """The packages a table for ``package_name`` covers.

    A TYPE -- ``"rch"`` -- is every package of it, array form included: MF6 has one
    recharge package, and ``RCHA`` is FloPy's class for its array form. Anything
    else is one package's own NAME -- ``"rch_array"``. Both through the model's own
    lookup (``model_packages_of_type``, then ``model.package``, ledger 352), which a
    duck-typed model answers too.
    """

    typed = model_packages_of_type(model, package_name)
    if typed:
        return typed
    found = model.package(package_name)  # a NAME, or raises naming what is not there
    return list(found) if isinstance(found, list) else [found]


def _one_package_table(model, package, package_name, *, per, layer_values, cell_values) -> pd.DataFrame:
    """One package's rows, each tagged with the package's own name."""

    if not hasattr(package, "stress_period_data"):
        # An array-form (READASARRAYS) package -- `mf.rch.array` / `mf.evt.array`,
        # or one loaded from a model that used them. It has no
        # `stress_period_data` at all; its fields are transient ARRAYS, one value
        # per column per period. Expanding them into the same tidy frame is what
        # lets the inputs tier, the diff tier and the group tier read an array
        # package at all -- without it they raise AttributeError.
        frame = _array_package_input_table(
            model,
            package,
            package_name,
            per=per,
            layer_values=layer_values,
            cell_values=cell_values,
        )
    else:
        if per is not None:
            period_items = [(int(per), package.stress_period_data.get_data(key=int(per)))]
        else:
            period_items = package.stress_period_data.get_data().items()
        frames: list[pd.DataFrame] = []
        for period, data in period_items:
            frame = pd.DataFrame(data).copy()
            if frame.empty:
                continue
            frame = split_cellid_columns(frame)
            if layer_values is not None and "layer" in frame.columns:
                frame = frame[frame["layer"].isin(layer_values)]
            if cell_values is not None and "cell" in frame.columns:
                frame = frame[frame["cell"].isin(cell_values)]
            if frame.empty:
                continue
            frame["package"] = str(package_name).lower()
            frame["per"] = int(period)
            frame["model"] = model.name
            frames.append(frame)
        frame = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not frame.empty:
        frame["pname"] = str(package.package_name).upper()
    return frame



def _array_package_input_table(
    model,
    package,
    package_name,
    *,
    per,
    layer_values,
    cell_values,
) -> pd.DataFrame:
    """Tidy an array-form (``READASARRAYS``) package into the cell-wise frame.

    Each transient array yields one value per column per period. The layer a
    column's flux lands in comes from the package's own ``irch``/``ievt`` index
    where it has one, and is layer 0 otherwise -- which is exactly MODFLOW 6's
    own rule for an array package without an index array.
    """

    index_array = None
    for candidate in ("irch", "ievt"):
        entry = getattr(package, candidate, None)
        if entry is not None and hasattr(entry, "get_data"):
            index_array = entry
            break

    fields = {
        field: getattr(package, field)
        for field in ("recharge", "rate", "depth", "surface")
        if getattr(package, field, None) is not None
        and hasattr(getattr(package, field), "get_data")
    }
    if not fields:
        return pd.DataFrame(columns=["model", "package", "per", "layer", "cell"])

    periods: dict[int, dict[str, np.ndarray]] = {}
    for field, entry in fields.items():
        data = entry.get_data()
        if not isinstance(data, dict):
            data = {0: data}
        for period, values in data.items():
            if values is None:
                continue
            periods.setdefault(int(period), {})[field] = np.asarray(values).ravel()

    frames: list[pd.DataFrame] = []
    for period in sorted(periods):
        if per is not None and int(period) != int(per):
            continue
        columns = periods[period]
        length = max(len(v) for v in columns.values())
        frame = pd.DataFrame({name: np.resize(values, length) for name, values in columns.items()})
        frame["cell"] = np.arange(length)

        layers = np.zeros(length, dtype=int)
        if index_array is not None:
            raw = index_array.get_data()
            entry = raw.get(period) if isinstance(raw, dict) else raw
            if entry is None and isinstance(raw, dict) and raw:
                entry = raw[max(k for k in raw if k <= period)] if any(k <= period for k in raw) else None
            if entry is not None:
                layers = np.resize(np.asarray(entry).ravel().astype(int), length)
        frame["layer"] = layers

        if layer_values is not None:
            frame = frame[frame["layer"].isin(layer_values)]
        if cell_values is not None:
            frame = frame[frame["cell"].isin(cell_values)]
        if frame.empty:
            continue
        frame["package"] = str(package_name).lower()
        frame["per"] = int(period)
        frame["model"] = model.name
        frames.append(frame)

    if not frames:
        return pd.DataFrame(columns=["model", "package", "per", "layer", "cell"])
    combined = pd.concat(frames, ignore_index=True)
    for column in ("layer", "cell", "per"):
        combined[column] = combined[column].astype(int)
    return combined


def build_uzf_field_input_table(
    model: SimulationBase,
    field_name: str,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
    package: str | None = None,
) -> pd.DataFrame:
    """Return a normalized table for one UZF perioddata field.

    Parameters
    ----------
    model
        Live or file-backed model object.
    field_name
        Perioddata field name such as ``"finf"``.
    package : str or None, optional
        The UZF package's NAME, any case. Default ``None``: the model's only UZF
        package -- several raise ``SeveralPackagesError`` naming them, because
        ``ifno`` numbers the UZF cells within ONE package (ledger 352).

    Returns
    -------
    pandas.DataFrame
        Tidy DataFrame with ``per``, ``ifno``, zero-based ``layer``/``cell``,
        and the requested field.
    """

    uzf = one_package(model, "uzf", package)
    name = package_name_of(uzf)
    layer_values = _normalize_iterable_filter(layer)
    cell_values = _normalize_iterable_filter(cells)
    # Keyed by the PACKAGE too: two UZF packages each have a `finf` (ledger 352).
    cache_key = ("uzf_field_input_table", str(field_name).lower(), name.lower())
    if per is None and layer_values is None and cell_values is None:
        cached = _get_package_explorer_cache(model).get(cache_key)
        if cached is not None:
            return cached.copy()
    ifno_to_cellid = package_outputs(model, "uzf", package).ifno_to_cellid
    layer_lookup = ifno_to_cellid["layer"].to_numpy()
    cell_lookup = ifno_to_cellid["cellid"].to_numpy()
    frames: list[pd.DataFrame] = []
    if per is not None:
        period_items = [(int(per), uzf.perioddata.get_data(key=int(per)))]
    else:
        period_items = uzf.perioddata.get_data().items()
    for period, data in period_items:
        if data is None or len(data) == 0:
            continue
        try:
            ifno_raw = _extract_structured_column(data, "ifno")
            value_raw = _extract_structured_column(data, field_name)
        except KeyError as exc:
            raise KeyError(f"UZF field {field_name!r} was not found.") from exc

        ifno_series = pd.to_numeric(pd.Series(ifno_raw), errors="coerce")
        valid_mask = ifno_series.notna().to_numpy()
        if not valid_mask.any():
            continue
        ifno_values = ifno_series.loc[valid_mask].astype(int).to_numpy()
        layer_values_for_rows = np.take(layer_lookup, ifno_values)
        cell_values_for_rows = np.take(cell_lookup, ifno_values)
        field_values = np.asarray(value_raw)[valid_mask]

        row_mask = np.ones(len(ifno_values), dtype=bool)
        if layer_values is not None:
            row_mask &= np.isin(layer_values_for_rows, layer_values)
        if cell_values is not None:
            row_mask &= np.isin(cell_values_for_rows, cell_values)
        if not row_mask.any():
            continue

        ifno_values = ifno_values[row_mask]
        layer_values_for_rows = layer_values_for_rows[row_mask]
        cell_values_for_rows = cell_values_for_rows[row_mask]
        field_values = field_values[row_mask]
        frames.append(
            pd.DataFrame(
                {
                    "model": model.name,
                    "package": "uzf",
                    "per": int(period),
                    "ifno": ifno_values,
                    "layer": layer_values_for_rows,
                    "cell": cell_values_for_rows,
                    field_name: field_values,
                }
            )
        )

    if not frames:
        return pd.DataFrame(
            columns=["model", "package", "per", "ifno", "layer", "cell", field_name]
        )

    combined = pd.concat(frames, ignore_index=True)
    for column in ("per", "ifno", "layer", "cell"):
        combined[column] = combined[column].astype(int)
    combined = _coerce_numeric_like_columns(combined, exclude={"model", "package"})
    if per is None and layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[cache_key] = combined.copy()
    return combined


def build_uzf_field_input_wide_table(
    model: SimulationBase,
    field_name: str,
    *,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
    package: str | None = None,
) -> pd.DataFrame:
    """Return a compact wide UZF table with one row per UZF record/cell.

    ``package`` names the UZF package (any case); ``None``, the default, is the
    model's only one -- several raise ``SeveralPackagesError`` (ledger 352).
    """

    uzf = one_package(model, "uzf", package)
    name = package_name_of(uzf)
    cache_key = ("uzf_field_input_wide_table", str(field_name).lower(), name.lower())
    layer_values = _normalize_iterable_filter(layer)
    cell_values = _normalize_iterable_filter(cells)
    if layer_values is None and cell_values is None:
        cached = _get_package_explorer_cache(model).get(cache_key)
        if cached is not None:
            return cached.copy()

    ifno_to_cellid = package_outputs(model, "uzf", package).ifno_to_cellid.copy()
    ifno_to_cellid.index.name = "ifno"
    base = ifno_to_cellid.reset_index().rename(columns={"cellid": "cell"})
    if layer_values is not None:
        base = base[base["layer"].isin(layer_values)]
    if cell_values is not None:
        base = base[base["cell"].isin(cell_values)]
    if base.empty:
        return pd.DataFrame(columns=["ifno", "layer", "cell"])

    wide = base.set_index("ifno").copy()
    period_items = uzf.perioddata.get_data().items()
    for period, data in period_items:
        if data is None or len(data) == 0:
            continue
        try:
            ifno_raw = _extract_structured_column(data, "ifno")
            value_raw = _extract_structured_column(data, field_name)
        except KeyError as exc:
            raise KeyError(f"UZF field {field_name!r} was not found.") from exc
        ifno_series = pd.to_numeric(pd.Series(ifno_raw), errors="coerce")
        value_series = pd.Series(value_raw)
        valid_mask = ifno_series.notna().to_numpy()
        if not valid_mask.any():
            continue
        period_series = pd.Series(
            value_series.loc[valid_mask].to_numpy(),
            index=ifno_series.loc[valid_mask].astype(int).to_numpy(),
            name=f"per_{int(period)}",
        )
        wide = wide.join(period_series, how="left")

    wide = wide.reset_index()
    wide = _coerce_numeric_like_columns(wide)
    if layer_values is None and cell_values is None:
        _get_package_explorer_cache(model)[cache_key] = wide.copy()
    return wide


def build_sfr_reach_table(model: SimulationBase, *, package: str | None = None) -> pd.DataFrame:
    """Return per-reach geometry metadata used by SFR tables and profiles.

    ``package`` names the SFR package (any case); ``None``, the default, is the
    model's only one -- several raise ``SeveralPackagesError``, because each
    package numbers its own reaches (ledger 352).
    """

    packagedata = pd.DataFrame(one_package(model, "sfr", package).packagedata.get_data()).copy()
    packagedata = split_cellid_columns(packagedata)
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "reach"})
    packagedata["reach"] = packagedata["reach"].astype(int)
    packagedata["layer"] = packagedata["layer"].astype(int)
    packagedata["cell"] = packagedata["cell"].astype(int)
    packagedata["rlen"] = pd.to_numeric(packagedata["rlen"], errors="coerce")
    packagedata = packagedata.sort_values("reach").reset_index(drop=True)
    packagedata["distance_start"] = packagedata["rlen"].fillna(
        0.0
    ).cumsum() - packagedata["rlen"].fillna(0.0)
    packagedata["distance_end"] = packagedata["distance_start"] + packagedata[
        "rlen"
    ].fillna(0.0)
    packagedata["distance_mid"] = packagedata["distance_start"] + (
        packagedata["rlen"].fillna(0.0) / 2.0
    )
    return packagedata[
        [
            "reach",
            "layer",
            "cell",
            "rlen",
            "distance_start",
            "distance_mid",
            "distance_end",
        ]
    ].drop_duplicates()


def build_lak_connection_table(model: SimulationBase, *, package: str | None = None) -> pd.DataFrame:
    """Return a normalized table of LAK connection geometry by cell.

    ``package`` names the LAK package (any case); ``None``, the default, is the
    model's only one -- several raise ``SeveralPackagesError``, because each
    package numbers its own lakes (ledger 352).

    Notes
    -----
    The returned table keeps one row per LAK connection record and computes a
    ``connection_area`` field that can be mapped directly as a cell
    choropleth:

    - vertical connections use the plan-view Voronoi cell area
    - horizontal/embedded horizontal connections use ``connwidth * (telev - belev)``

    The values therefore represent the physical interface area available for
    exchange between the lake and groundwater.
    """

    connectiondata = pd.DataFrame(one_package(model, "lak", package).connectiondata.get_data()).copy()
    if connectiondata.empty:
        return pd.DataFrame(
            columns=[
                "model",
                "package",
                "lake",
                "iconn",
                "layer",
                "cell",
                "claktype",
                "belev",
                "telev",
                "connlen",
                "connwidth",
                "connection_area",
            ]
        )

    connectiondata = split_cellid_columns(connectiondata)
    if "ifno" in connectiondata.columns:
        connectiondata = connectiondata.rename(columns={"ifno": "lake"})
    connectiondata["lake"] = pd.to_numeric(
        connectiondata["lake"], errors="coerce"
    ).astype("Int64")
    if "iconn" in connectiondata.columns:
        connectiondata["iconn"] = pd.to_numeric(
            connectiondata["iconn"], errors="coerce"
        ).astype("Int64")
    connectiondata["layer"] = pd.to_numeric(
        connectiondata["layer"], errors="coerce"
    ).astype(int)
    connectiondata["cell"] = pd.to_numeric(
        connectiondata["cell"], errors="coerce"
    ).astype(int)
    connectiondata["package"] = "lak"
    connectiondata["model"] = model.name
    if "claktype" in connectiondata.columns:
        connectiondata["claktype"] = (
            connectiondata["claktype"].astype("string").str.upper()
        )
    else:
        connectiondata["claktype"] = "UNKNOWN"

    cell_area_lookup = {
        int(cell): float(area)
        for cell, area in enumerate(
            np.asarray(model.vor.area_list, dtype=float).reshape(-1)
        )
    }
    top = pd.to_numeric(connectiondata.get("telev"), errors="coerce")
    bottom = pd.to_numeric(connectiondata.get("belev"), errors="coerce")
    width = pd.to_numeric(connectiondata.get("connwidth"), errors="coerce")
    vertical_mask = connectiondata["claktype"].astype(str).str.upper().eq("VERTICAL")
    thickness = (top - bottom).clip(lower=0.0)
    horizontal_area = width * thickness
    vertical_area = connectiondata["cell"].map(cell_area_lookup).astype(float)
    connectiondata["connection_area"] = np.where(
        vertical_mask, vertical_area, horizontal_area
    )

    ordered = [
        "model",
        "package",
        "lake",
        "iconn",
        "layer",
        "cell",
        "claktype",
        "belev",
        "telev",
        "connlen",
        "connwidth",
        "connection_area",
    ]
    remaining = [column for column in connectiondata.columns if column not in ordered]
    return connectiondata[ordered + remaining]


#: MF6's LAK connection types whose lake bottom comes from a lake TABLE, which
#: :func:`lake_bottoms` does not read (``gwf-lak.f90``, ``lak_read_table``).
_EMBEDDED_LAK_CONNECTIONS = ("EMBEDDEDH", "EMBEDDEDV")


def _cell_top_and_bottom(model: SimulationBase):
    """MF6's ``dis%top`` / ``dis%bot`` of a cell, looked up by its FloPy ``cellid``.

    The model's discretization package, found by its exact TYPE (FloPy's
    ``get_package("dis")`` prefix-matches a DISV). A cell's top is the model top
    in layer 0 and the bottom of the layer above it below that, as MF6's
    ``Disv.f90`` and ``Dis.f90`` fill ``top``; a DISU cell carries its own.
    """

    gwf = model.gwf
    dis = next(
        (package for package in gwf.packagelist
         if getattr(package, "package_type", None) in ("dis", "disv", "disu")),
        None,
    )
    if dis is None:
        raise AttributeError(
            f"Model '{model.name}' has no DIS, DISV or DISU package, so its cells' "
            "elevations -- which MF6 derives a lake's bottom from -- are unknown."
        )
    if dis.package_type == "disu":
        top = np.asarray(dis.top.array, dtype=float).reshape(-1)
        bottom = np.asarray(dis.bot.array, dtype=float).reshape(-1)
        return lambda cellid: (float(top[cellid[0]]), float(bottom[cellid[0]]))
    top = np.asarray(dis.top.array, dtype=float)
    botm = np.asarray(dis.botm.array, dtype=float)
    surfaces = np.concatenate([top[np.newaxis, ...], botm], axis=0)

    def lookup(cellid):
        layer, rest = int(cellid[0]), tuple(int(value) for value in cellid[1:])
        return float(surfaces[(layer, *rest)]), float(surfaces[(layer + 1, *rest)])

    return lookup


def lake_bottoms(
    model: SimulationBase, *, lakes: Iterable[int] | None = None, package: str | None = None
) -> pd.Series:
    """Each lake's BOTTOM as MF6 computes it: the stage MF6 gives a DRY lake.

    MF6 derives a lake's bottom (``lakebot``) from its connections and writes it
    nowhere, yet it is the number a dry lake reads: the solver never lets a
    lake's stage below it, the STAGE observation reports that clamped stage,
    and the stage file writes -1e30 once the stage is at or below it
    (``gwf-lak.f90``, ``lak_ot_dv``). So it is what a dry lake is CALIBRATED
    against (``LakeStageTargets``, ledger 351). The rule, from
    ``gwf-lak.f90`` (``lak_read_lake_connections``):

    - a VERTICAL connection gives its cell's TOP -- MF6 overwrites the BELEV
      written for a vertical connection, so a bed declared inside the cell does
      not lower it (ledger 327);
    - a HORIZONTAL connection gives its BELEV as written, or its cell's BOTTOM
      when BELEV equals TELEV exactly (MF6 then takes the whole cell face);
    - the lake's bottom is the LOWEST of these over its connections.

    Read from the model's own LAK and discretization packages -- exactly what
    MF6 read when the model was loaded from its files (``mf.load_mf6_run``,
    which PEST's forward run uses; measured equal to MF6's observation to the
    bit). A model still held in memory from its build carries its elevations
    unrounded, where FloPy WROTE them rounded -- a number between 0.001 and
    100,000 at 8 decimals -- so there the two differ by up to 5e-9 (measured
    8.8e-10 and 1.2e-9, ledger 351).

    Parameters
    ----------
    model : SimulationBase
        A model with a LAK package and a DIS, DISV or DISU package; run or not --
        nothing here reads an output file.
    lakes : iterable of int or None, optional
        The lakes to compute, as zero-based MF6 lake numbers (each passed through
        ``int``). Default ``None``: every lake of the package, ``0 ..
        nlakes - 1``. Only the lakes asked for are checked, so an EMBEDDED lake
        (below) refuses only when it is one of them.
    package : str or None, optional
        The LAK package's NAME, any case. Default ``None``: the model's only LAK
        package (ledger 352).

    Returns
    -------
    pandas.Series
        Named ``"lake_bottom"``, indexed by zero-based lake number in the order
        asked (ascending for ``None``), in the model's length units.

    Raises
    ------
    AttributeError
        When the model has no LAK package, or no discretization package.
    SeveralPackagesError
        When the model carries several LAK packages and ``package`` names none
        (``package_registry.one_package``).
    ValueError
        When ``package`` names no LAK package of the model; when a lake asked for
        is not a lake of the package, has no connection,
        or connects ``EMBEDDEDH``/``EMBEDDEDV``: MF6 takes an embedded lake's
        bottom from its cell's bottom and then from the leading all-zero rows of
        its lake TABLE (``lak_read_table``), which this does not read. MF6's
        own STAGE observation of that lake reports it
        (``LakeStageTargets.attach_flopy_obs``).

    Examples
    --------
    >>> lake_bottoms(model)            # one VERTICAL connection under a 100 ft top
    0    100.0
    Name: lake_bottom, dtype: float64
    """

    lak = one_package(model, "lak", package)
    nlakes = int(lak.nlakes.get_data())
    wanted = list(range(nlakes)) if lakes is None else [int(lake) for lake in lakes]
    unknown = [lake for lake in wanted if not 0 <= lake < nlakes]
    if unknown:
        raise ValueError(
            f"Lake(s) {unknown} are not lakes of model '{model.name}'s LAK package, "
            f"which numbers its {nlakes} lake(s) 0..{nlakes - 1} (zero-based)."
        )
    connections = pd.DataFrame(lak.connectiondata.get_data())
    connections["ifno"] = connections["ifno"].astype(int)
    connections["claktype"] = connections["claktype"].astype(str).str.strip().str.upper()
    elevations = _cell_top_and_bottom(model)

    bottoms: dict[int, float] = {}
    for lake in wanted:
        rows = connections.loc[connections["ifno"] == lake]
        if rows.empty:
            raise ValueError(
                f"Lake {lake} of model '{model.name}' has no connection in its LAK "
                "CONNECTIONDATA, so MF6 gives it no bottom."
            )
        embedded = sorted(set(rows["claktype"]) & set(_EMBEDDED_LAK_CONNECTIONS))
        if embedded:
            raise ValueError(
                f"Lake {lake} of model '{model.name}' connects {embedded[0]}: MF6 takes "
                "an embedded lake's bottom from its cell's bottom and then from the "
                "leading all-zero rows of its lake TABLE (gwf-lak.f90, lak_read_table), "
                "which this does not read. MF6's own STAGE observation of the lake "
                "reports it: declare one (LakeStageTargets(...).attach_flopy_obs)."
            )
        candidates = []
        for row in rows.itertuples(index=False):
            top, bottom = elevations(tuple(row.cellid))
            if row.claktype == "VERTICAL":
                candidates.append(top)
            elif float(row.belev) == float(row.telev):
                candidates.append(bottom)
            else:
                candidates.append(float(row.belev))
        bottoms[lake] = min(candidates)
    return pd.Series(bottoms, name="lake_bottom", dtype=float).rename_axis("lake")


def build_sfr_input_table(model: SimulationBase, *, package: str | None = None) -> pd.DataFrame:
    """Return normalized static and stress-period SFR inputs mapped to cells.

    ``package`` names the SFR package (any case); ``None``, the default, is the
    model's only one -- several raise ``SeveralPackagesError`` (ledger 352).
    """

    sfr = one_package(model, "sfr", package)
    base = pd.DataFrame(sfr.packagedata.get_data()).copy()
    base = split_cellid_columns(base)
    if "ifno" in base.columns:
        base = base.rename(columns={"ifno": "reach"})
    base["reach"] = pd.to_numeric(base["reach"], errors="coerce").astype(int)
    base["layer"] = pd.to_numeric(base["layer"], errors="coerce").astype(int)
    base["cell"] = pd.to_numeric(base["cell"], errors="coerce").astype(int)
    settings = _numeric_period_settings(sfr, id_column="ifno", nper=model.nper)
    frames = []
    for period in range(int(model.nper)):
        frame = base.copy()
        for reach, values in settings[period].items():
            mask = frame["reach"] == int(reach)
            for name, value in values.items():
                frame.loc[mask, name] = value
        frame["per"] = period
        frame["package"] = "sfr"
        frame["model"] = model.name
        frames.append(frame)
    return _coerce_numeric_like_columns(
        pd.concat(frames, ignore_index=True), exclude={"model", "package"}
    )


def build_lak_input_table(model: SimulationBase, *, package: str | None = None) -> pd.DataFrame:
    """Return normalized static and stress-period LAK inputs mapped to connection cells.

    ``package`` names the LAK package (any case); ``None``, the default, is the
    model's only one -- several raise ``SeveralPackagesError`` (ledger 352).
    """

    lak = one_package(model, "lak", package)
    connections = build_lak_connection_table(model, package=package)
    packagedata = pd.DataFrame(lak.packagedata.get_data()).copy()
    if "ifno" in packagedata.columns:
        packagedata = packagedata.rename(columns={"ifno": "lake"})
    packagedata["lake"] = pd.to_numeric(packagedata["lake"], errors="coerce").astype(
        int
    )
    base = connections.merge(
        packagedata, on="lake", how="left", suffixes=("", "_package")
    )
    settings = _numeric_period_settings(lak, id_column="number", nper=model.nper)
    frames = []
    for period in range(int(model.nper)):
        frame = base.copy()
        for lake, values in settings[period].items():
            mask = frame["lake"] == int(lake)
            for name, value in values.items():
                frame.loc[mask, name] = value
        frame["per"] = period
        frame["package"] = "lak"
        frame["model"] = model.name
        frames.append(frame)
    return _coerce_numeric_like_columns(
        pd.concat(frames, ignore_index=True), exclude={"model", "package"}
    )


def summarize_input_table(
    frame: pd.DataFrame, *, label: str, value_columns: list[str]
) -> pd.DataFrame:
    """Build a compact one-row summary for a normalized input table."""

    if frame.empty:
        return pd.DataFrame(
            [
                {
                    "label": label,
                    "records": 0,
                    "periods": 0,
                    "layers": 0,
                    "cells": 0,
                    "value_columns": value_columns,
                }
            ]
        )

    return pd.DataFrame(
        [
            {
                "label": label,
                "records": int(len(frame)),
                "periods": int(frame["per"].nunique()) if "per" in frame.columns else 0,
                "layers": int(frame["layer"].nunique())
                if "layer" in frame.columns
                else 0,
                "cells": int(frame["cell"].nunique()) if "cell" in frame.columns else 0,
                "value_columns": value_columns,
            }
        ]
    )


def _numeric_period_settings(
    package, *, id_column: str, nper: int
) -> dict[int, dict[int, dict[str, float]]]:
    """Return carried-forward numeric package settings keyed by period and feature id."""

    current: dict[int, dict[str, float]] = {}
    by_period: dict[int, dict[int, dict[str, float]]] = {}
    perioddata = package.perioddata.get_data()
    for period in range(int(nper)):
        records = perioddata.get(period)
        if records is not None:
            frame = pd.DataFrame(records).copy()
            if not frame.empty:
                identifier = (
                    id_column if id_column in frame.columns else frame.columns[0]
                )
                setting_column = next(
                    (
                        column
                        for column in frame.columns
                        if str(column).lower().endswith("setting")
                    ),
                    None,
                )
                data_column = next(
                    (
                        column
                        for column in frame.columns
                        if str(column).lower().endswith("setting_data")
                    ),
                    None,
                )
                if setting_column is not None and data_column is not None:
                    for row in frame.itertuples(index=False):
                        feature_id = int(getattr(row, identifier))
                        setting = (
                            str(getattr(row, setting_column))
                            .strip()
                            .lower()
                            .replace("-", "_")
                        )
                        raw_value = getattr(row, data_column)
                        numeric = pd.to_numeric(
                            pd.Series([raw_value]), errors="coerce"
                        ).iloc[0]
                        if pd.notna(numeric):
                            current.setdefault(feature_id, {})[setting] = float(numeric)
        by_period[period] = {
            feature_id: dict(settings) for feature_id, settings in current.items()
        }
    return by_period


__all__ = [
    "build_cell_package_input_table",
    "build_uzf_field_input_table",
    "build_uzf_field_input_wide_table",
    "build_sfr_reach_table",
    "build_lak_connection_table",
    "lake_bottoms",
    "build_sfr_input_table",
    "build_lak_input_table",
    "summarize_input_table",
    "_numeric_period_settings",
]
