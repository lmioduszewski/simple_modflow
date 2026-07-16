"""Small shared helpers for boundary/workflow modules."""

from __future__ import annotations

from collections.abc import Iterable
from numbers import Real


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
        return dict(zip(cell_list, values))
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
