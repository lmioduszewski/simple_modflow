"""Shared engine for file-less, domain-aware areal packages (RCH, EVT).

:class:`_ArealBuilder` owns everything the recharge and evapotranspiration
builders have in common: selecting cells from the model domain (top-active /
all-active / explicit), broadcasting scalar / per-cell / per-period values across
those cells, naming boundaries, and validating the result. Concrete builders
subclass it and supply only what differs -- the per-cell record shape
(:meth:`_row_for`), the spec factory (:meth:`_make_spec`), and the build metadata
(:meth:`_build_metadata`); a subclass that precomputes per-build data shared
across periods (e.g. EVT's per-cell surface) overrides :meth:`_prepared`.

This is the engine under the ``mf.rch(context=...)`` / ``mf.evt(context=...)``
package-first facade; use those, not this class, to build models.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Real
from typing import Any

import numpy as np

from myflopy.specs import ModelContext, PackageSpec

CellId = tuple[int, int]


@dataclass(frozen=True, slots=True, kw_only=True)
class _ArealBuilder:
    """Base for domain-aware areal package builders (see the module docstring).

    Cell ids are DISV-style ``(layer, cell)`` tuples. Integer cell ids are also
    accepted and are placed on ``layer``. Use tuple keys when a mapping might
    otherwise look like period keys, such as ``{0: ...}``.
    """

    context: ModelContext
    nper: int
    cells: str | Sequence[int | CellId] = "top_active"
    layer: int = 0
    name_by_cell: Mapping[int | CellId, str] | None = None
    boundnames: bool = False
    name: str = "areal"
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate ``nper >= 1`` and that a string ``cells`` selector is a known keyword."""

        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        if isinstance(self.cells, str) and self.cells not in {"top_active", "all_active", "surface_only"}:
            raise ValueError("cells must be 'top_active', 'surface_only', 'all_active', or explicit cells.")

    def with_updates(self, **updates: Any) -> _ArealBuilder:
        """Return a changed builder without modifying the original."""

        return replace(self, **updates)

    @property
    def domain(self):
        """Return the active-domain array carried by the model context."""

        return None if self.context.domain is None else np.asarray(self.context.domain)

    @property
    def selected_cells(self) -> tuple[CellId, ...]:
        """Return target cells as DISV-style ``(layer, cell)`` IDs."""

        if not isinstance(self.cells, str):
            return tuple(self._normalize_cell(cell) for cell in self.cells)

        domain = self.domain
        if domain is None:
            grid = self.context.grid
            if grid is None or not hasattr(grid, "ncpl"):
                raise ValueError("context.domain or context.grid.ncpl is required for generated cells.")
            # NOTE: with no idomain, 'top_active' and 'all_active' are identical
            # (every column has one nominal layer). Recorded as a deferred quirk.
            return tuple((self.layer, cell) for cell in range(int(grid.ncpl)))

        if domain.ndim == 1:
            active = np.flatnonzero(domain > 0)
            return tuple((self.layer, int(cell)) for cell in active)

        selected: list[CellId] = []
        for cell in range(domain.shape[1]):
            active_layers = np.flatnonzero(domain[:, cell] > 0)
            if len(active_layers) == 0:
                continue
            if self.cells == "all_active":
                selected.extend((int(layer), int(cell)) for layer in active_layers)
            else:
                # 'top_active' and (the currently unhandled) 'surface_only' both
                # take the first active layer per column. Recorded as a deferred quirk.
                selected.append((int(active_layers[0]), int(cell)))
        return tuple(selected)

    def _normalize_cell(self, cell: int | CellId) -> CellId:
        """Coerce a cell spec to a ``(layer, cell)`` pair (a bare int lands on ``self.layer``)."""

        if isinstance(cell, Real):
            return (int(self.layer), int(cell))
        if not isinstance(cell, Sequence) or isinstance(cell, str) or len(cell) != 2:
            raise ValueError("Explicit cells must be integer cells or (layer, cell) pairs.")
        layer, node = cell
        return (int(layer), int(node))

    def _cell_lookup_keys(self, cellid: CellId) -> tuple[Any, ...]:
        """The keys to try when looking up a cell in a user mapping: the pair, then the bare cell."""

        return (cellid, cellid[1])

    def _name_for_cell(self, cellid: CellId) -> str:
        """The boundname for a cell -- from ``name_by_cell`` if present, else ``<name>_<layer>_<cell>``."""

        if not self.name_by_cell:
            return f"{self.name}_{cellid[0]}_{cellid[1]}"
        for key in self._cell_lookup_keys(cellid):
            if key in self.name_by_cell:
                return str(self.name_by_cell[key])
        return f"{self.name}_{cellid[0]}_{cellid[1]}"

    @staticmethod
    def _period_keys(value: Mapping[Any, Any], nper: int) -> bool:
        """Whether a mapping is keyed by stress-period indices (all int keys in ``range(nper)``)."""

        return bool(value) and all(isinstance(key, int) and key in range(nper) for key in value)

    def _value_for_cell(
        self, value: Any, cellid: CellId, selected_cells: tuple[CellId, ...], *, label: str
    ) -> Any:
        """Resolve one cell's value from a cell mapping, a scalar/field, or a per-cell sequence."""

        if isinstance(value, Mapping):
            for key in self._cell_lookup_keys(cellid):
                if key in value:
                    return value[key]
            raise ValueError(f"{label} is missing cell {cellid}.")
        if isinstance(value, Real) or isinstance(value, str):
            return value
        values = list(value)
        if len(values) != len(selected_cells):
            raise ValueError(f"{label} sequence must contain one value per selected cell.")
        return values[selected_cells.index(cellid)]

    def _value_for_period_cell(
        self, value: Any, period: int, cellid: CellId, selected_cells: tuple[CellId, ...], *, label: str
    ) -> Any:
        """Resolve one cell's value for a given period (unwrapping a per-period mapping first)."""

        if isinstance(value, Mapping) and self._period_keys(value, self.nper):
            if period not in value:
                raise ValueError(f"{label} is missing period {period}.")
            value = value[period]
        return self._value_for_cell(value, cellid, selected_cells, label=label)

    # --- subclass hooks -------------------------------------------------------

    def _prepared(self, cells: tuple[CellId, ...]) -> Any:
        """Precompute per-build data shared across periods (e.g. EVT surface). Default: none."""

        return None

    def _row_for(self, period: int, cellid: CellId, cells: tuple[CellId, ...], prepared: Any) -> list[Any]:
        """Return one MF6 record (without the trailing boundname). Subclasses implement this."""

        raise NotImplementedError

    def _make_spec(self, stress_period_data: dict[int, list[list[Any]]]) -> PackageSpec:
        """Build the package spec from stress-period data. Subclasses implement this."""

        raise NotImplementedError

    def _build_metadata(self) -> dict[str, Any]:
        """Manifest metadata stamped onto the built spec. Subclasses implement this."""

        raise NotImplementedError

    # --- shared assembly ------------------------------------------------------

    @property
    def stress_period_data(self) -> dict[int, list[list[Any]]]:
        """Return normalized MF6 stress-period data for the selected cells."""

        cells = self.selected_cells
        prepared = self._prepared(cells)
        data: dict[int, list[list[Any]]] = {period: [] for period in range(self.nper)}
        for period in range(self.nper):
            for cellid in cells:
                row = self._row_for(period, cellid, cells, prepared)
                if self.boundnames:
                    row = [*row, self._name_for_cell(cellid)]
                data[period].append(row)
        return data

    def validate(self) -> None:
        """Raise clear errors for empty selection or inactive/out-of-domain cells."""

        cells = self.selected_cells
        if not cells:
            raise ValueError(f"{type(self).__name__} selected no active cells.")
        domain = self.domain
        if domain is not None and domain.ndim == 1:
            for _layer, cell in cells:
                if cell >= domain.shape[0] or domain[cell] <= 0:
                    raise ValueError(f"cell {cell} is inactive or outside the domain.")
        if domain is not None and domain.ndim > 1:
            for layer, cell in cells:
                if layer >= domain.shape[0] or cell >= domain.shape[1] or domain[layer, cell] <= 0:
                    raise ValueError(f"cell {(layer, cell)} is inactive or outside the domain.")
        self.stress_period_data

    def build(self) -> PackageSpec:
        """Validate stored configuration and return its package spec."""

        self.validate()
        return self._make_spec(self.stress_period_data).with_metadata(**self._build_metadata())
