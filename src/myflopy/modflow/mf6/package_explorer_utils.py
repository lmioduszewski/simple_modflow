"""Cross-cutting helpers for package explorer modules."""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
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


#: The cellid each discretization writes, by FloPy's ``grid_type``: the three grid
#: types :func:`cell_numbering` reads.
CELLID_FORMS = {
    "structured": "(layer, row, col)",
    "vertex": "(layer, cell)",
    "unstructured": "(node,)",
}


@dataclass(frozen=True, slots=True)
class CellNumbering:
    """How a model's cellids name the ``(layer, cell)`` its views are keyed by.

    ONE rule for every reader that turns a FloPy cellid into a cell (ledger
    419): the package input tables, the SFR/LAK tables, the diff's value tier
    and the region registry. Every view of a model numbers a cell ``(layer, cell)``
    with ``cell`` the in-layer index the model's own grid, its heads table and its
    cell polygons use, so a cellid maps as

    * DIS ``(layer, row, col)`` -> ``(layer, row * ncol + col)``, FloPy's row-major
      order (``split_cellid_columns`` read the ROW as the cell, so two identical
      DIS runs diffed as changed);
    * DISV ``(layer, cell)`` -> itself;
    * DISU ``(node,)`` -> ``(0, node)``: MF6 numbers a DISU model's cells by node
      alone and writes its heads as one record, so its heads table and its cell
      polygons are one layer of ``nodes`` cells -- whatever layering FloPy infers
      from ``IHC``.

    A cellid of another length, or off the grid, raises ``ValueError`` naming the
    form: off the grid, a DIS column past ``ncol`` would otherwise land on the next
    row's cell.

    Attributes
    ----------
    kind : str
        FloPy's ``grid_type``: ``"structured"``, ``"vertex"`` or ``"unstructured"``.
    nlay : int
        Layers the views number: the grid's layers, and 1 on DISU.
    ncpl : int
        Cells per layer the views number: ``nrow * ncol`` on DIS, the DISV's
        ``ncpl``, and the node count on DISU.
    ncol : int or None
        Columns of a DIS grid; ``None`` on DISV and DISU.
    """

    kind: str
    nlay: int
    ncpl: int
    ncol: int | None = None

    @property
    def form(self) -> str:
        """The cellid this discretization writes, e.g. ``"(layer, row, col)"``."""

        return CELLID_FORMS[self.kind]

    @property
    def length(self) -> int:
        """How many items that cellid has: 3 on DIS, 2 on DISV, 1 on DISU."""

        return {"structured": 3, "vertex": 2, "unstructured": 1}[self.kind]

    def layer_cells(self, cellids: Iterable) -> tuple[np.ndarray, np.ndarray]:
        """Each cellid as ``(layer, cell)``, as two integer arrays.

        Parameters
        ----------
        cellids : iterable of tuple
            FloPy cellids, each a tuple of this grid's :attr:`length` integers
            (a bare integer is read as a one-item tuple, so it is a DISU node and
            refused elsewhere).

        Returns
        -------
        tuple of numpy.ndarray
            ``(layers, cells)``, one entry per cellid, in order.

        Raises
        ------
        ValueError
            A cellid has another number of items, or names a layer, row, column,
            cell or node off the grid. The message names the offending cellids and
            this grid's form.
        """

        items = [tuple(value) if isinstance(value, (tuple, list, np.ndarray, np.void))
                 else (value,) for value in cellids]
        wrong = [item for item in items if len(item) != self.length]
        if wrong:
            raise ValueError(
                f"Cellid {wrong[0]!r} has {len(wrong[0])} item(s), but on this "
                f"{self.kind} grid a cellid is {self.form}"
                + (f" ({len(wrong)} such cellids)." if len(wrong) > 1 else ".")
            )
        if not items:
            empty = np.zeros(0, dtype=int)
            return empty, empty.copy()
        table = np.asarray(items, dtype=int).reshape(len(items), self.length)
        if self.kind == "structured":
            nrow = self.ncpl // self.ncol
            layers, rows, cols = table[:, 0], table[:, 1], table[:, 2]
            off = (
                (layers < 0) | (layers >= self.nlay) | (rows < 0) | (rows >= nrow)
                | (cols < 0) | (cols >= self.ncol)
            )
            bounds = (
                f"layers 0..{self.nlay - 1}, rows 0..{nrow - 1} and columns "
                f"0..{self.ncol - 1}"
            )
            cells = rows * self.ncol + cols
        elif self.kind == "vertex":
            layers, cells = table[:, 0], table[:, 1]
            off = (layers < 0) | (layers >= self.nlay) | (cells < 0) | (cells >= self.ncpl)
            bounds = f"layers 0..{self.nlay - 1} and cells 0..{self.ncpl - 1}"
        else:
            cells = table[:, 0]
            layers = np.zeros_like(cells)
            off = (cells < 0) | (cells >= self.ncpl)
            bounds = f"nodes 0..{self.ncpl - 1}"
        if off.any():
            bad = [items[index] for index in np.flatnonzero(off)[:5]]
            raise ValueError(
                f"Cellid(s) {', '.join(repr(item) for item in bad)} are off the grid: "
                f"a zero-based {self.form} here has {bounds}"
                + (f" ({int(off.sum())} such cellids)." if off.sum() > 1 else ".")
            )
        return layers.astype(int), cells.astype(int)

    def layer_cell(self, cellid) -> tuple[int, int]:
        """One cellid as ``(layer, cell)``; see :meth:`layer_cells`.

        Parameters
        ----------
        cellid : tuple
            One FloPy cellid of this grid's :attr:`length` integers.

        Returns
        -------
        tuple of (int, int)
            The ``(layer, cell)`` the model's views number it.
        """

        layers, cells = self.layer_cells([cellid])
        return int(layers[0]), int(cells[0])

    def cellid(self, layer: int, cell: int) -> tuple[int, ...]:
        """The FloPy cellid of the ``(layer, cell)`` the model's views number: the inverse
        of :meth:`layer_cell`.

        ``(layer, cell // ncol, cell % ncol)`` on DIS, ``(layer, cell)`` on DISV and
        ``(cell,)`` on DISU, whose views number one layer of every node. What a
        writer of a cellid the model's own package will read uses -- a DRN zone's
        OBS records wrote ``(layer, cell)`` on every grid, which MF6 reads as no
        cell on DIS and DISU (ledger 455).

        Parameters
        ----------
        layer : int
            The zero-based layer, ``0..nlay - 1`` (``0`` alone on DISU).
        cell : int
            The zero-based in-layer cell number (the node on DISU), ``0..ncpl - 1``.

        Returns
        -------
        tuple of int
            The cellid, zero-based as FloPy holds it, :attr:`length` items long.

        Raises
        ------
        ValueError
            When the layer or the cell is off the grid (:meth:`check`).
        """

        self.check(layers=[layer], cells=[cell])
        layer, cell = int(layer), int(cell)
        if self.kind == "structured":
            return (layer, cell // self.ncol, cell % self.ncol)
        if self.kind == "vertex":
            return (layer, cell)
        return (cell,)

    def check(self, *, layers: Iterable[int] = (), cells: Iterable[int] = ()) -> None:
        """Refuse a layer or an in-layer cell number off the grid.

        Parameters
        ----------
        layers : iterable of int, optional
            Zero-based layers to check against ``0..nlay - 1``. Default none.
        cells : iterable of int, optional
            Zero-based in-layer cell numbers (a node on DISU) to check against
            ``0..ncpl - 1``. Default none.

        Raises
        ------
        ValueError
            Naming the first values off the grid and the grid's range.
        """

        bad_layers = sorted({int(v) for v in layers if not 0 <= int(v) < self.nlay})
        if bad_layers:
            span = "one layer, 0" if self.nlay == 1 else f"layers 0..{self.nlay - 1}"
            raise ValueError(
                f"Layer(s) {bad_layers[:5]} are off the grid: this {self.kind} grid has "
                f"{span} (zero-based)."
            )
        bad_cells = sorted({int(v) for v in cells if not 0 <= int(v) < self.ncpl})
        if bad_cells:
            noun = "node" if self.kind == "unstructured" else "in-layer cell"
            raise ValueError(
                f"Cell(s) {bad_cells[:5]} are off the grid: this {self.kind} grid has "
                f"{noun} numbers 0..{self.ncpl - 1} (zero-based)."
            )


def _model_grid(source):
    """The FloPy grid of ``source``: its ``modelgrid``, or ``source`` itself when it has none.

    Parameters
    ----------
    source : object
        A model (anything with a ``modelgrid``) or a FloPy grid.

    Returns
    -------
    object
        ``source.modelgrid`` when it exists, else ``source``.
    """

    return getattr(source, "modelgrid", source)


def has_cell_numbering(source) -> bool:
    """Whether :func:`cell_numbering` can read ``source``'s cellids.

    Parameters
    ----------
    source : object or None
        A model, a FloPy grid, or ``None``.

    Returns
    -------
    bool
        ``True`` when ``source``'s grid is of one of the three types in
        :data:`CELLID_FORMS`; ``False`` for ``None``, a model with no
        discretization (FloPy's untyped base grid, ``grid_type`` ``None``) or a
        double without a grid.
    """

    return source is not None and getattr(_model_grid(source), "grid_type", None) in CELLID_FORMS


def views_nlay(source) -> int:
    """How many layers ``source``'s views number: the grid's, and ONE on a DISU.

    The one rule every per-layer view reads its layer count by (ledger
    448): :attr:`CellNumbering.nlay` when the grid is
    typed, else FloPy's ``modelgrid.nlay``. FloPy INFERS a DISU's layering from
    ``IHC`` (two on a 2-layer fixture), while MF6 writes its heads as one record
    of nodes and every table keys a node as ``(0, node)`` -- so a view counting
    FloPy's layers invents a layer 1 that duplicates layer 0.

    Parameters
    ----------
    source : object
        A model (anything with a ``modelgrid``, such as a FloPy GWF model or a
        myflopy model's ``gwf``) or a FloPy grid. A source whose grid has no
        ``grid_type`` (a test double, a model with no discretization) falls back
        to its ``nlay``.

    Returns
    -------
    int
        The layer count: ``nlay`` on DIS and DISV, ``1`` on DISU.
    """

    if has_cell_numbering(source):
        return int(cell_numbering(source).nlay)
    return int(_model_grid(source).nlay)


def cell_numbering(source) -> CellNumbering:
    """The :class:`CellNumbering` of a model or of a FloPy model grid.

    Parameters
    ----------
    source : model or flopy.discretization.Grid
        Anything with a ``modelgrid`` (a myflopy model, a FloPy model), or a FloPy
        grid itself. Its ``grid_type`` says how its cellids read.

    Returns
    -------
    CellNumbering
        How that model's cellids map to ``(layer, cell)``.

    Raises
    ------
    ValueError
        The grid's ``grid_type`` is none of ``"structured"``, ``"vertex"`` and
        ``"unstructured"`` (a model with no discretization package), so no cellid
        can be read.
    """

    grid = _model_grid(source)
    kind = getattr(grid, "grid_type", None)
    if kind == "structured":
        return CellNumbering(
            kind, nlay=int(grid.nlay), ncpl=int(grid.nrow) * int(grid.ncol), ncol=int(grid.ncol)
        )
    if kind == "vertex":
        return CellNumbering(kind, nlay=int(grid.nlay), ncpl=int(grid.ncpl))
    if kind == "unstructured":
        return CellNumbering(kind, nlay=1, ncpl=int(np.sum(np.asarray(grid.ncpl, dtype=int))))
    raise ValueError(
        f"Cannot read a cellid on a grid of type {kind!r}: a cellid is "
        "(layer, row, col) on DIS, (layer, cell) on DISV and (node,) on DISU, and this "
        "model's grid is none of them (does it have a discretization package?)."
    )


def cell_index_column(values) -> pd.Series:
    """A ``layer`` or ``cell`` column as integers, nullable where a cell is missing.

    ``int`` when every value is present; pandas' nullable ``Int64`` when one is not --
    an SFR reach connected to NO cell (MF6 ``NONE``), which
    :func:`split_cellid_columns` gives a missing ``layer`` and ``cell``. A plain
    ``astype(int)`` raised on it in every SFR table.

    Parameters
    ----------
    values : pandas.Series
        The column, numeric or numeric-like; a value that is not a number reads as
        missing.

    Returns
    -------
    pandas.Series
        Same index, ``int64`` or ``Int64``.
    """

    numbers = pd.to_numeric(pd.Series(values), errors="coerce")
    if numbers.isna().any():
        return numbers.astype("Int64")
    return numbers.astype(int)


def split_cellid_columns(frame: pd.DataFrame, numbering: CellNumbering) -> pd.DataFrame:
    """Expand a FloPy ``cellid`` column into zero-based ``layer`` and ``cell``.

    Parameters
    ----------
    frame : pandas.DataFrame
        Records that may carry a FloPy ``cellid`` column; one without it is
        returned as a copy, unchanged.
    numbering : CellNumbering
        The model's numbering (:func:`cell_numbering`), which reads each cellid:
        ``(layer, row, col)`` on DIS is ``(layer, row * ncol + col)``, ``(layer,
        cell)`` on DISV is itself, ``(node,)`` on DISU is ``(0, node)``. A cellid of
        ``None`` or MF6's ``NONE`` (an SFR reach connected to no cell, which FloPy
        reads back as the string ``'none'``) gives a missing value in both columns, which
        are then pandas' nullable ``Int64`` (:func:`cell_index_column`).

    Returns
    -------
    pandas.DataFrame
        A copy of ``frame`` with ``cellid`` replaced by ``layer`` and ``cell``.

    Raises
    ------
    ValueError
        A cellid is not this grid's form, or is off the grid
        (:meth:`CellNumbering.layer_cells`).
    """

    if "cellid" not in frame.columns:
        return frame.copy()

    expanded = frame.copy()
    raw = expanded["cellid"].tolist()
    # MF6 writes an SFR reach connected to no cell as NONE, which FloPy reads back as
    # the string 'none' (a reach built in memory may hold None).
    present = [
        index for index, cellid in enumerate(raw)
        if cellid is not None
        and not (isinstance(cellid, str) and cellid.strip().lower() == "none")
    ]
    layers: list[int | None] = [None] * len(raw)
    cells: list[int | None] = [None] * len(raw)
    if present:
        found_layers, found_cells = numbering.layer_cells([raw[index] for index in present])
        for index, layer, cell in zip(present, found_layers, found_cells, strict=True):
            layers[index] = int(layer)
            cells[index] = int(cell)
    expanded["layer"] = cell_index_column(pd.Series(layers, index=expanded.index, dtype=object))
    expanded["cell"] = cell_index_column(pd.Series(cells, index=expanded.index, dtype=object))
    return expanded.drop(columns=["cellid"])


def _normalize_term_filter(term: str | Iterable[str] | None) -> list[str] | None:
    """Normalize an optional package-budget term filter to uppercase strings."""

    if term is None:
        return None
    if isinstance(term, str):
        return [str(term).strip().upper()]
    return [str(value).strip().upper() for value in term]


def budget_term_attribute(term: str) -> str:
    """Return the Python attribute name that reaches one MF6 budget record.

    ``"SOURCE-SINK MIX"`` -> ``source_sink_mix``, ``"FLOW-JA-FACE"`` ->
    ``flow_ja_face``, ``"STO-SS"`` -> ``sto_ss``, ``"UZF-GWRCH"`` ->
    ``uzf_gwrch``. MF6 separates words with hyphens and, in one case, a space;
    both collapse to the underscore.

    The one normalizer to run this direction. It is what lets
    ``ModelBudgetNamespace.__getitem__`` accept either spelling, and what the LAK/SFR
    ``budget.get(term=...)`` refusal names a term by (``"ext_inflow" is spelled
    "EXT-INFLOW"``; ledger 473) -- that filter compares MODFLOW 6's spelling
    only, and returned an empty table in silence for the attribute one. Lives here,
    in import layer 0, so the budget-table builders below the views can ask it.

    Parameters
    ----------
    term : str
        A budget record name as MODFLOW 6 writes it, any case and padding.

    Returns
    -------
    str
        The lower-case attribute spelling, words joined by ``_``.
    """

    return re.sub(r"[^0-9a-z]+", "_", str(term).strip().lower()).strip("_")


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


def _stress_period_count(dataset) -> int | None:
    """The simulation's ``nper``, as the FloPy dataset's own TDIS states it.

    ``None`` when the dataset cannot say (a test double, a package detached from a
    simulation): the period is then not range-checked.
    """

    try:
        model_dim = dataset.data_dimensions.package_dim.model_dim[0]
        return int(model_dim.simulation_time.get_num_stress_periods())
    except (AttributeError, IndexError, TypeError, ValueError, MFDataException, FlopyException):
        logger.debug("the dataset states no stress-period count; per is not range-checked")
        return None


def period_number(per) -> int:
    """``per`` as a Python ``int`` stress period, or ``TypeError`` naming it.

    The type half of :func:`checked_period`, for a caller that must refuse a
    mistyped period BEFORE a ``PACKAGE_TABLE_UNAVAILABLE`` handler (which catches
    ``TypeError``) could read it as a missing package -- the group tier.

    Parameters
    ----------
    per : int
        Zero-based stress period: a Python or NumPy integer, or a float holding a
        whole number. ``bool``, any other float, a string and ``None`` raise
        ``TypeError``. The range is not checked here.

    Returns
    -------
    int
        ``per`` as a Python ``int``.
    """

    if isinstance(per, bool) or not isinstance(per, (int, np.integer, float, np.floating)):
        raise TypeError(f"per must be a zero-based stress period (an int), not {per!r}")
    if not float(per).is_integer():
        raise TypeError(f"per must be a whole stress period number, not {per!r}")
    return int(per)


def checked_period(dataset, per) -> int:
    """``per`` as a zero-based stress period of ``dataset``'s simulation, or refused.

    Parameters
    ----------
    dataset : flopy.mf6.data.mfdata.MFTransientData
        The FloPy transient dataset ``per`` indexes (``stress_period_data``,
        ``perioddata`` or one array of an array-form package); its simulation's
        TDIS supplies ``nper``.
    per : int
        Zero-based stress period. A Python or NumPy integer, or a float that is a
        whole number; ``bool``, any other float and a string raise ``TypeError``.

    Returns
    -------
    int
        ``per`` as a Python ``int``.

    Raises
    ------
    TypeError
        ``per`` is not an integer.
    IndexError
        ``per`` is negative, or not below the simulation's ``nper`` (when the
        dataset states one). ``IndexError`` rather than ``ValueError`` on purpose:
        the group and diff tiers read a ``ValueError`` from a table builder as "this
        member does not carry the package" (``PACKAGE_TABLE_UNAVAILABLE``), which
        would turn a mistyped period into a silently empty table.
    """

    period = period_number(per)
    nper = _stress_period_count(dataset)
    if period < 0 or (nper is not None and period >= nper):
        bound = "" if nper is None else f" (nper={nper}: periods 0 to {nper - 1})"
        raise IndexError(f"per={per!r} is not a stress period of this simulation{bound}")
    return period


def written_period_blocks(dataset) -> dict[int, object | None]:
    """Every PERIOD block a FloPy transient dataset writes, read without changing it.

    The one reader of a package's period blocks for the read views (ledger
    449). FloPy's ``get_data(key=N)`` -- and ``has_data(key=N)`` -- CREATE the
    storage for a period they are asked about, and wherever FloPy still holds a block
    header for that period (a header outlives ``remove_transient_key``) it writes
    that storage as an EMPTY ``PERIOD`` block, which MF6 reads as "no records from
    here on": a view that read GHB at such a period turned 84 records into 0 from
    that period on, under Normal termination. This asks only for keys already stored.

    Parameters
    ----------
    dataset : flopy.mf6.data.mfdata.MFTransientData or None
        A transient list or array dataset -- a list package's
        ``stress_period_data``, an advanced package's ``perioddata``, one array of
        an array-form package. ``None`` gives ``{}``.

    Returns
    -------
    dict of int to object or None
        ``{period: data}`` in period order, one entry per block the package file
        writes: the records (a recarray or DataFrame) or array FloPy stores there,
        or ``None`` for an EMPTY block -- one FloPy marks empty (``empty_keys``,
        what ``{N: []}`` declares and a file's ``BEGIN PERIOD N`` / ``END PERIOD``
        loads as), or empty storage under a block header FloPy already holds. Empty
        storage with NO header -- what any ``get_data(key=N)`` leaves behind, a
        user's own included -- is not a block: FloPy writes nothing for it, and MF6
        runs the block before it there (ledger 449).
    """

    if dataset is None:
        return {}
    return dict(_written_blocks(dataset))


def period_block_in_force(dataset, per) -> tuple[int | None, object | None]:
    """The block MODFLOW 6 runs a list package's period ``per`` on, without changing it.

    MF6's rule for a list boundary package (CHD, DRN, GHB, RIV, WEL, list RCH/EVT)
    and for HFB: a ``PERIOD`` block replaces every record before it, and a period
    with no block of its own reuses the last block before it -- so a package written
    once, in period 0, runs those records in every period, and an EMPTY block
    switches the package off until the next block (ledger 449). Read through
    :func:`written_period_blocks`, so asking never adds a block.

    Parameters
    ----------
    dataset : flopy.mf6.data.mfdata.MFTransientData
        The package's ``stress_period_data`` (or any transient dataset whose blocks
        replace one another).
    per : int
        Zero-based stress period, checked by :func:`checked_period`: a non-integer
        raises ``TypeError``, a period outside ``0 .. nper - 1`` ``IndexError``.

    Returns
    -------
    tuple of (int or None, object or None)
        The period of the block in force and its data; ``(None, None)`` before the
        first block, and ``(block, None)`` when the block in force is empty -- MF6
        applies no records in either case.
    """

    period = checked_period(dataset, per)
    blocks = list(_written_blocks(dataset, upto=period))
    return blocks[-1] if blocks else (None, None)


def feature_records_in_force(dataset, per, *, id_field: str) -> pd.DataFrame | None:
    """The settings MODFLOW 6 holds for each feature of an advanced package at ``per``.

    UZF's rule (``gwf-uzf.f90``'s ``uzf_rp``), unlike a list package's: a ``PERIOD``
    block sets the values of the features it LISTS and leaves every other feature at
    what an earlier block set, so the settings in force are, feature by feature, its
    row in the last block at or before ``per`` that lists it. An empty block changes
    nothing. Read through :func:`written_period_blocks`, so asking never adds a block
    (ledger 449).

    Parameters
    ----------
    dataset : flopy.mf6.data.mfdata.MFTransientData
        The package's ``perioddata``.
    per : int
        Zero-based stress period, checked by :func:`checked_period` (``TypeError``
        for a non-integer, ``IndexError`` outside ``0 .. nper - 1``).
    id_field : str
        The record field naming the feature -- ``"ifno"`` for UZF.

    Returns
    -------
    pandas.DataFrame or None
        One row per feature some block at or before ``per`` lists -- the record's
        own fields as FloPy stores them -- ordered by ``id_field``; ``None`` when
        no such block lists any.
    """

    period = checked_period(dataset, per)
    frames = []
    for _key, data in _written_blocks(dataset, upto=period):
        if data is None:
            continue
        frame = pd.DataFrame(data)
        if id_field not in frame.columns:
            raise KeyError(f"the period records carry no {id_field!r} field")
        frames.append(frame)
    if not frames:
        return None
    combined = pd.concat(frames, ignore_index=True)
    # A later block's row for a feature replaces an earlier one; the rest stand.
    latest = combined.drop_duplicates(subset=id_field, keep="last")
    return latest.sort_values(id_field, kind="stable").reset_index(drop=True)


def _written_blocks(dataset, *, upto: int | None = None):
    """``(period, data)`` for every block ``dataset`` WRITES, in period order.

    FloPy's own write rule (``MFBlock.write``): a key whose storage holds records is
    written, so is one it marks EMPTY (``empty_keys``), and so is empty storage under
    a block header FloPy already holds (a header outlives ``remove_transient_key``,
    and a read re-creating the storage then writes as an empty block). Empty storage
    with NO header -- which ``get_data(key=N)`` leaves, whoever called it -- is
    written as nothing, so it is no block and MF6 runs the block before it (ledger
    449). Read off the dataset's own bookkeeping, never through
    ``get_data(key=)`` for a key it does not store. ``upto`` stops after that period.
    ``data`` is ``None`` for an empty block.
    """

    stored = {
        int(key) for key in (dataset.get_active_key_dict() if hasattr(dataset, "get_active_key_dict") else {})
        if isinstance(key, (int, np.integer))
    }
    empty = {
        int(key) for key, flag in (getattr(dataset, "empty_keys", None) or {}).items()
        if flag and isinstance(key, (int, np.integer))
    }
    for key in sorted(stored | empty):
        if upto is not None and key > upto:
            return
        data = _block_data(dataset, key, key in stored)
        if data is None and key not in empty and not _block_header_exists(dataset, key):
            continue
        yield key, data


def _block_header_exists(dataset, key: int) -> bool:
    """Whether FloPy holds a ``PERIOD`` block header for ``key`` (what it writes).

    ``True`` when the dataset names no block to ask (a test double), so its empty
    storage is read as an empty block.
    """

    block = getattr(dataset, "_block", None)
    header_exists = getattr(block, "header_exists", None)
    if header_exists is None:
        return True
    return bool(header_exists(key))


def _block_data(dataset, key: int, stored: bool):
    """The data FloPy stores for an EXISTING block ``key``; ``None`` when it holds none.

    ``get_data(key=)`` adds no storage for a key already stored, which is the only
    key this is asked about.
    """

    if not stored:
        return None
    data = dataset.get_data(key=key)
    if data is None:
        return None
    if isinstance(data, pd.DataFrame):
        return None if data.empty else data
    try:
        return None if len(data) == 0 else data
    except TypeError:  # a 0-d array or a scalar constant: a value, not a record list
        return data


def _aggregate_hover_strings(series: pd.Series) -> str:
    """Aggregate a cell's non-numeric values into a concise hover string."""

    values = [str(value) for value in series.dropna().tolist()]
    if not values:
        return ""
    unique_values = list(dict.fromkeys(values))
    return ", ".join(unique_values)


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
    "CELLID_FORMS",
    "CellNumbering",
    "cell_index_column",
    "cell_numbering",
    "has_cell_numbering",
    "views_nlay",
    "split_cellid_columns",
    "_normalize_term_filter",
    "budget_term_attribute",
    "_normalize_connection_type_filter",
    "_normalize_surface_water_include",
    "_infer_default_value_column",
    "_normalize_iterable_filter",
    "_filter_normalized_table",
    "_aggregate_hover_strings",
    "_normalize_budget_nodes",
    "_model_budget_record_frame",
]
