"""Small shared helpers for boundary/workflow modules."""

from __future__ import annotations

from collections.abc import Iterable
from numbers import Real

import numpy as np


def normalize_grid_type(grid_type: str) -> str:
    """Validate and normalize a MODFLOW grid type string."""

    normalized = grid_type.lower()
    if normalized not in {"disu", "disv"}:
        raise ValueError(f"grid_type must be 'disu' or 'disv', got {grid_type!r}")
    return normalized


def build_cell_id(cell: int, grid_type: str = "disv", layer: int | None = 0):
    """Build the MF6-style cell id for a DISV or DISU cell reference."""

    grid_type = normalize_grid_type(grid_type)
    if grid_type == "disu":
        return int(cell)
    if layer is None:
        raise ValueError("layer is required when grid_type='disv'")
    return int(layer), int(cell)


def filter_inactive_cells(cells: Iterable[int], inactive_cells: Iterable[int] | None = None) -> list[int]:
    """Remove inactive cells from an iterable of candidate cell ids."""

    if inactive_cells is None:
        return [int(cell) for cell in cells]
    inactive = {int(cell) for cell in inactive_cells}
    return [int(cell) for cell in cells if int(cell) not in inactive]


def coerce_values_by_cell(
    cells: Iterable[int],
    values,
    *,
    name: str,
) -> dict[int, float | int | object]:
    """Coerce scalar/list/dict inputs into a ``cell -> value`` mapping."""

    cell_list = [int(cell) for cell in cells]
    if isinstance(values, dict):
        missing = [cell for cell in cell_list if cell not in values]
        if missing:
            raise ValueError(f"{name} is missing values for cells: {missing}")
        return {cell: values[cell] for cell in cell_list}
    if isinstance(values, Real):
        return {cell: values for cell in cell_list}
    if isinstance(values, list):
        if len(values) != len(cell_list):
            raise ValueError(f"{name} list length must equal number of cells")
        return dict(zip(cell_list, values, strict=False))
    raise TypeError(f"{name} must be a scalar, list, or dict")


def merge_stress_period_data(
    existing: dict[int, list[list]],
    updates: dict[int, list[list]],
    *,
    replace: bool = True,
) -> dict[int, list[list]]:
    """Merge stress-period data dictionaries keyed by stress period."""

    merged: dict[int, list[list]] = {per: [row.copy() for row in rows] for per, rows in existing.items()}
    for per, rows in updates.items():
        current = {row[0]: row[1:] for row in merged.get(per, [])}
        for row in rows:
            cell_id = row[0]
            if cell_id in current and not replace:
                current[cell_id][0] = current[cell_id][0] + row[1]
            else:
                current[cell_id] = row[1:]
        merged[per] = [[cell_id, *values] for cell_id, values in current.items()]
    return merged


def expand_periodic_cell_input(
    value,
    *,
    nper: int,
    nitems: int,
    default=0.0,
) -> dict[int, list]:
    """Expand scalar or sparse per-period input into a full MF6-style mapping."""

    if value is None:
        return {per: [default] * nitems for per in range(nper)}
    if isinstance(value, Real):
        return {per: [value] * nitems for per in range(nper)}
    if isinstance(value, dict):
        expanded: dict[int, list] = {}
        for per in range(nper):
            period_values = value.get(per, [default] * nitems)
            if isinstance(period_values, Real):
                expanded[per] = [period_values] * nitems
                continue
            if len(period_values) != nitems:
                raise ValueError(
                    f"period {per} must contain {nitems} values, got {len(period_values)}"
                )
            expanded[per] = list(period_values)
        return expanded
    raise TypeError("period input must be None, a scalar, or a dict keyed by stress period")


def bottom_clearance(bottom):
    """The smallest gap above a cell bottom that survives the trip through MF6's files.

    A boundary elevation set *exactly* equal to its cell bottom is not safe. FloPy
    writes ``botm`` fixed-point and a list-BC value in E format, so the two numbers
    MODFLOW 6 compares are not the two numbers Python held: the discarded digits can
    fall either way, and MF6's comparison is strict. Ledger 290 measured 13 of 24 real
    runs FAILING that way at bottoms drawn from [40, 600] -- round numbers being the
    accidental survivors -- so the gap is relative, not a flat ``1e-6`` (which is itself
    too small at a bottom of 3456.789012346).

    Accepts a scalar or an array and answers in kind. A scalar comes back as a plain
    ``float`` on purpose: a ``numpy.float64`` returned here would ride into a written
    boundary record, and a numpy scalar at a serialization boundary is this project's
    recurring bug class (ledger 265).
    """

    clearance = np.maximum(np.abs(bottom) * 1e-6, 1e-6)
    return float(clearance) if np.ndim(bottom) == 0 else clearance


def river_pair_lift(stage, rbot, floor):
    """Lift a river bottom to ``floor``, and carry its stage only if the stage needs it.

    **A RIV record is a pair, not two independent elevations.** MODFLOW 6 compares
    ``rbot`` against the cell bottom and enforces ``stage >= rbot`` separately, so
    raising ``rbot`` alone turns a model that ran into a hard failure whenever the stage
    was below the bottom too -- and clamping both to the floor collapses the declared
    river depth to zero.

    But a stage that is ALREADY inside its cell is a legal, physical number: it is the
    water surface the river imposes, and moving it moves the answer. So it is carried up
    by the same delta as ``rbot`` **only when it is itself below the floor**, which is
    exactly when MODFLOW 6 would refuse it.

    Measured on a one-cell bench, cell bottom 48.0 and cell top 60.0, at both
    ``icelltype`` settings (``docs/reproducers/riv_pair_and_icelltype.py``)::

        declared          rule            written        q @ ict0   q @ ict1
        stage 40/rbot 35  lift both       53.0 / 48.0     -109.09     -77.93
                          rbot only       40.0 / 48.0     REFUSED     REFUSED
                          clamp both      48.0 / 48.0     -381.82    -196.46  depth 0
        stage 55/rbot 40  as declared     55.0 / 40.0        0.00     REFUSED
                          lift both       63.0 / 48.0     +436.37    +378.82  ^ above top
                          rbot only       55.0 / 48.0        0.00       0.00

    Lifting the stage unconditionally put the water surface 3 ft ABOVE land surface on
    the second row and injected 436 ft3/d where the river as drawn moved nothing; with a
    deeper bed (stage 50 / rbot 20) it reversed the sign, -272.73 to +1254.55. Ledger 290
    measured only the first row and generalised from it; ledger 291 measured the second.
    This rule reproduces the best reading of both.

    Works elementwise, so a caller may pass scalars or arrays already masked down to the
    records that offend. A scalar stage comes back a plain ``float``: a ``numpy.float64``
    here would ride into a written boundary record (ledger 265). Returns ``(stage, rbot)``
    -- the order they appear in both a USG and a MODFLOW 6 river record.
    """

    lifted = np.where(stage >= floor, stage, stage + (floor - rbot))
    return (float(lifted) if np.ndim(stage) == 0 else lifted), floor
