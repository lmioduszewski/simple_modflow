"""Cross-cutting helpers for package explorer modules."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from flopy.mf6.mfbase import FlopyException, MFDataException

from myflopy._logging import get_logger

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

logger = get_logger(__name__)

#: What "open this model's budget file and ask it something" can raise.
#:
#: ``_get_budget_reader`` is ``self.gwf.output.budget()``, and flopy RETURNS
#: None from that when the ``.cbc`` is missing (it catches ``OSError``
#: internally) -- so the dominant failure is an ``AttributeError`` on None, not
#: a file error. The rest were measured against flopy 3.10: ``ValueError`` on an
#: empty or truncated budget file, ``NotImplementedError`` from the base
#: ``Grid.shape``/``.nnodes`` that ``CellBudgetFile.__init__`` touches
#: unconditionally when a model has no discretization package, ``OSError`` on an
#: unreadable path, and flopy's own two exception classes from a lazy load.
BUDGET_READER_UNAVAILABLE = (
    AttributeError, TypeError, ValueError, NotImplementedError, OSError,
    MFDataException, FlopyException,
)


def _get_package_explorer_cache(model: SimulationBase) -> dict[tuple, pd.DataFrame]:
    """Return a per-model cache for expensive normalized explorer tables."""

    cache = getattr(model, "_package_explorer_cache", None)
    if cache is None:
        cache = {}
        model._package_explorer_cache = cache
    return cache


def _extract_structured_column(data, column_name: str):
    """Return one column from structured MF6 record data with a safe fallback."""

    dtype = getattr(data, "dtype", None)
    names = getattr(dtype, "names", None)
    if names and column_name in names:
        return np.asarray(data[column_name])
    frame = pd.DataFrame(data).copy()
    if column_name not in frame.columns:
        raise KeyError(column_name)
    return frame[column_name].to_numpy()


def _coerce_numeric_like_columns(
    frame: pd.DataFrame, *, exclude: Iterable[str] = ()
) -> pd.DataFrame:
    """Convert object/string columns to numeric when every non-null value is numeric-like."""

    excluded = {str(column) for column in exclude}
    coerced = frame.copy()
    for column in coerced.columns:
        if str(column) in excluded:
            continue
        series = coerced[column]
        if not (
            pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series)
        ):
            continue
        non_null = series.notna()
        if not non_null.any():
            continue
        converted = pd.to_numeric(series, errors="coerce")
        if converted.loc[non_null].notna().all():
            coerced[column] = converted
    return coerced


def split_cellid_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Expand a FloPy ``cellid`` column into zero-based ``layer`` and ``cell``.

    Parameters
    ----------
    frame
        Input DataFrame that may contain a FloPy-style ``cellid`` object
        column. The function accepts structured tuples such as ``(layer, cell)``
        as well as simpler scalar cell identifiers.

    Returns
    -------
    pandas.DataFrame
        A copy of ``frame`` where ``cellid`` has been replaced with zero-based
        ``layer`` and ``cell`` columns when possible.
    """

    if "cellid" not in frame.columns:
        return frame.copy()

    expanded = frame.copy()
    layers: list[int | None] = []
    cells: list[int | None] = []
    for cellid in expanded["cellid"].tolist():
        if isinstance(cellid, tuple) and len(cellid) >= 2:
            layers.append(int(cellid[0]))
            cells.append(int(cellid[1]))
        elif cellid is None:
            layers.append(None)
            cells.append(None)
        else:
            layers.append(0)
            cells.append(int(cellid))
    expanded["layer"] = layers
    expanded["cell"] = cells
    return expanded.drop(columns=["cellid"])


def _normalize_term_filter(term: str | Iterable[str] | None) -> list[str] | None:
    """Normalize an optional package-budget term filter to uppercase strings."""

    if term is None:
        return None
    if isinstance(term, str):
        return [str(term).strip().upper()]
    return [str(value).strip().upper() for value in term]


def _normalize_connection_type_filter(
    connection_type: str | Iterable[str] | None,
) -> list[str] | None:
    """Normalize optional connection-type filters to uppercase labels."""

    if connection_type is None:
        return None
    if isinstance(connection_type, str):
        return [connection_type.upper()]
    normalized: list[str] = []
    for value in connection_type:
        normalized.append(str(value).upper())
    return normalized


def _normalize_surface_water_include(include: str | Iterable[str] | None) -> list[str]:
    """Normalize selected surface-water package names."""

    if include is None:
        normalized = ["sfr", "lak"]
    elif isinstance(include, str):
        normalized = [include.lower()]
    else:
        normalized = [str(value).lower() for value in include]
    allowed = {"sfr", "lak"}
    invalid = sorted(set(normalized) - allowed)
    if invalid:
        raise ValueError(
            f"Unsupported surface-water packages: {invalid!r}. Allowed values are 'sfr' and 'lak'."
        )
    ordered: list[str] = []
    for package_name in ("sfr", "lak"):
        if package_name in normalized:
            ordered.append(package_name)
    return ordered


def _infer_default_value_column(
    frame: pd.DataFrame, *, fallback: str | None = None
) -> str:
    """Infer the primary numeric value column for a normalized input table."""

    if fallback is not None and fallback in frame.columns:
        return fallback

    excluded = {"model", "package", "per", "layer", "cell", "ifno"}
    numeric = [
        column
        for column in frame.columns
        if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
    ]
    if not numeric:
        raise ValueError(
            "Could not infer a numeric value column for this package table."
        )
    return numeric[0]


def _normalize_iterable_filter(values: int | Iterable[int] | None) -> list[int] | None:
    """Normalize optional scalar-or-iterable selectors to integer lists."""

    if values is None:
        return None
    if isinstance(values, Iterable) and not isinstance(values, (str, bytes)):
        return [int(value) for value in values]
    return [int(values)]


def _filter_normalized_table(
    frame: pd.DataFrame,
    *,
    per: int | None = None,
    layer: int | Iterable[int] | None = None,
    cells: Iterable[int] | None = None,
) -> pd.DataFrame:
    """Apply standard ``per/layer/cell`` filters to a normalized input table."""

    selected = frame.copy()
    if per is not None and "per" in selected.columns:
        selected = selected[selected["per"] == int(per)]
    layer_values = _normalize_iterable_filter(layer)
    if layer_values is not None and "layer" in selected.columns:
        selected = selected[selected["layer"].isin(layer_values)]
    cell_values = _normalize_iterable_filter(cells)
    if cell_values is not None and "cell" in selected.columns:
        selected = selected[selected["cell"].isin(cell_values)]
    return selected.reset_index(drop=True)


def _aggregate_hover_strings(series: pd.Series) -> str:
    """Aggregate a cell's non-numeric values into a concise hover string."""

    values = [str(value) for value in series.dropna().tolist()]
    if not values:
        return ""
    unique_values = list(dict.fromkeys(values))
    return ", ".join(unique_values)


def _default_show_layer_elevs(model) -> bool:
    """Return whether choropleths should show layer elevations by default."""

    return getattr(model.vor, "gdf_topbtm", None) is not None


def _normalize_budget_nodes(
    frame: pd.DataFrame, *, columns: tuple[str, ...] = ("node", "node2")
) -> pd.DataFrame:
    """Normalize budget node columns to zero-based indexing.

    ``columns`` selects which id columns to shift; the default covers both and
    is what every package-output caller wants. The model-budget path splits the
    two -- :func:`_model_budget_record_frame` owns ``node`` (which is always a
    model cell) and re-enters here with ``columns=("node2",)`` for the second
    id, whose meaning is package-specific.
    """

    normalized = frame.copy()
    for column in columns:
        if column not in normalized.columns:
            continue
        numeric = pd.to_numeric(normalized[column], errors="coerce")
        mask = numeric.notna()
        if not mask.any():
            continue
        ints = numeric.loc[mask].astype(int)
        normalized[column] = numeric
        if int(ints.min()) >= 1:
            normalized.loc[mask, column] = (ints - 1).astype(float)
        else:
            normalized.loc[mask, column] = ints.astype(float)
    return normalized


#: MF6 model-budget records written as imeth=1 full arrays that are indexed by
#: cell CONNECTION rather than by cell. This is an MF6 *file-format* fact, not a
#: per-package one, which is why it is a literal here rather than a registry
#: field. If MF6 ever adds a second connection-indexed full-array record, it
#: belongs in this set -- see ledger 96 for why size cannot be used instead.
_CONNECTION_INDEXED_BUDGET_RECORDS = frozenset({"FLOW-JA-FACE"})


def _resolve_budget_record_names(model, budget_text: str) -> list[str]:
    """Return the real record name(s) a caller's ``budget_text`` selects.

    Callers may pass a substring -- ``model.bud("flow")`` reaches
    ``FLOW-JA-FACE`` -- so a guard that tests the caller's string directly would
    miss exactly the alias that needs catching. The budget reader is cached on
    the model (``_bud``), so this costs one lookup, not a file reopen.
    """

    try:
        reader = model._get_budget_reader()
        names = [
            str(name).strip()
            for name in reader.get_unique_record_names(decode=True)
        ]
    except BUDGET_READER_UNAVAILABLE:
        logger.debug(
            "could not list budget record names; matching %r literally instead",
            budget_text, exc_info=True,
        )
        return [str(budget_text).strip()]

    wanted = str(budget_text).strip().upper()
    exact = [name for name in names if name.upper() == wanted]
    if exact:
        return exact
    return [name for name in names if wanted in name.upper()] or [wanted]


def _model_budget_record_frame(
    model: SimulationBase,
    record,
    *,
    budget_text: str,
) -> pd.DataFrame:
    """Return one MF6 *model*-budget record as a frame with a zero-based ``node``.

    MF6 writes model-budget records in two shapes, and only one of them is the
    ``(node, node2, q)`` recarray this layer used to assume:

    * **imeth 2/5 -- a recarray** (``DRN``, ``CHD``, ``SOURCE-SINK MIX``, ...):
      one row per boundary, whose ``node`` is a 1-based model cell.
    * **imeth 1 -- a plain full array** (``STO-SS``, ``STORAGE-AQUEOUS``,
      ``STORAGE-CELLBLK``, ``FLOW-JA-FACE``): one value per *position*, with no
      ``node`` column at all.

    The full-array shape used to be dropped on the floor. ``from_records`` turns
    a ``(nlay, 1, ncpl)`` float array into a nonsense ``(1, 1)`` frame with an
    integer column name, which the caller's ``"node" not in frame.columns``
    guard then skipped -- so **two of a transport model's three budget terms
    came back as an empty table with no error at all** (measured 2026-07-27 on
    real GWT and GWE runs).

    A full array is only cell-mappable when it is indexed by cell. Records that
    are not -- ``FLOW-JA-FACE``, indexed by cell *connection* -- raise here
    rather than inventing a mapping. Raising is the point: the bug being fixed
    was silence, and a term that cannot be a cell table must say so.

    That test is made on the record's NAME, not on its length, because **length
    cannot tell the two apart**. Under an idomain reduction MF6 expands cell
    arrays back to ``nodesuser`` while ``FLOW-JA-FACE`` stays at the reduced
    ``nja``, so the two counts are independent and do collide -- measured
    2026-07-27 on a 1-layer DISV with ``ncpl=4`` and two active adjacent cells,
    where ``nja == nodesuser == 4`` and *both* records arrive with the identical
    shape ``(1, 1, 4)``. A size-only guard accepted that FLOW-JA-FACE as four
    cells and reported flow through two ``idomain=0`` cells, with no error.

    ``node2`` is deliberately left as MF6 wrote it; its meaning is
    package-specific (a boundary index for list packages, a feature id for
    SFR/LAK/UZF), so each caller applies its own rule.
    """

    array = np.asarray(record)
    if array.dtype.names:
        frame = pd.DataFrame.from_records(record).copy()
        if "node" not in frame.columns:
            raise ValueError(
                f"Budget record {budget_text!r} has no 'node' column "
                f"(columns: {list(frame.columns)}), so its rows cannot be "
                "mapped to model cells."
            )
        return _normalize_budget_nodes(frame, columns=("node",))

    values = array.ravel()
    node_count = len(model.node_to_lni)
    connection_indexed = [
        name
        for name in _resolve_budget_record_names(model, budget_text)
        if name.strip().upper() in _CONNECTION_INDEXED_BUDGET_RECORDS
    ]
    if connection_indexed:
        raise ValueError(
            f"Budget record {budget_text!r} resolves to "
            f"{connection_indexed[0]!r}, which is indexed by cell CONNECTION "
            f"rather than by cell, so it cannot be mapped to cells. (Its "
            f"length, {values.size}, can coincidentally equal the model's node "
            f"count on an idomain-reduced model, so size alone cannot tell "
            f"them apart.)"
        )
    if values.size != node_count:
        raise ValueError(
            f"Budget record {budget_text!r} is a full array of {values.size} "
            f"values, which is not one per model node ({node_count}), so it "
            "cannot be mapped to cells."
        )
    # Full arrays are POSITIONAL: index 0 is node 0. There is no 1-based shift
    # to undo here, unlike the recarray branch above.
    return pd.DataFrame(
        {
            "node": np.arange(node_count, dtype=float),
            "q": values.astype(float),
        }
    )


__all__ = [
    "_get_package_explorer_cache",
    "_extract_structured_column",
    "_coerce_numeric_like_columns",
    "split_cellid_columns",
    "_normalize_term_filter",
    "_normalize_connection_type_filter",
    "_normalize_surface_water_include",
    "_infer_default_value_column",
    "_normalize_iterable_filter",
    "_filter_normalized_table",
    "_aggregate_hover_strings",
    "_default_show_layer_elevs",
    "_normalize_budget_nodes",
    "_model_budget_record_frame",
]
