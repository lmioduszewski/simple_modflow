"""Build MODFLOW 6 UZF package specifications."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Real
from typing import Any

import numpy as np

from myflopy.advanced import uzf_spec
from myflopy.specs import ModelContext, PackageSpec


CellId = tuple[int, int]
CellSelection = str | Sequence[CellId]


def _is_scalar(value: Any) -> bool:
    """True for a single scalar UZF input (a number or field-name string, not a per-cell sequence)."""

    return isinstance(value, Real) or isinstance(value, str)


@dataclass(frozen=True, slots=True)
class UZFBuilder:
    """Prepare one UZF package from a model context and stored configuration.

    By default, every vertically connected active cell beneath land surface is
    represented by a UZF cell. Only the uppermost UZF cell in each column is a
    land-surface cell and receives period inputs such as infiltration and PET.
    """

    context: ModelContext
    nper: int
    vks: Any
    thtr: Any
    thts: Any
    thti: Any
    cells: CellSelection = "all_active"
    eps: Any = 4.0
    surfdep: Any = 0.001
    finf: Any = 0.0
    pet: Any = None
    extdp: Any = None
    extwc: Any = None
    ha: Any = None
    hroot: Any = None
    rootact: Any = None
    name: str = "uzf"
    boundnames: bool = False
    mover: bool = False
    simulate_et: bool | None = None
    linear_gwet: bool = False
    square_gwet: bool = False
    simulate_gwseep: bool = False
    unsat_etwc: bool = False
    unsat_etae: bool = False
    ntrailwaves: int = 7
    nwavesets: int = 40
    options: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate ``nper``, the ``cells`` keyword, and that a domain exists for auto cell selection."""

        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        if isinstance(self.cells, str) and self.cells not in {"all_active", "surface_only"}:
            raise ValueError("cells must be 'all_active', 'surface_only', or explicit cell IDs.")
        if isinstance(self.cells, str) and self.context.domain is None:
            raise ValueError("ModelContext.domain is required for automatic UZF cell selection.")

    def with_updates(self, **updates: Any) -> UZFBuilder:
        """Return a changed builder without modifying the original."""

        return replace(self, **updates)

    @property
    def domain(self) -> np.ndarray | None:
        """Return the active-domain array normalized to ``(nlay, ncpl)``."""

        if self.context.domain is None:
            return None
        domain = np.asarray(self.context.domain)
        if domain.ndim == 1:
            domain = domain[np.newaxis, :]
        if domain.ndim != 2:
            raise ValueError("ModelContext.domain must be a 1D or 2D array.")
        return domain

    @property
    def uzf_cells(self) -> tuple[CellId, ...]:
        """Return the configured UZF cells in feature-number order."""

        if not isinstance(self.cells, str):
            selected = tuple((int(layer), int(cell)) for layer, cell in self.cells)
            if len(set(selected)) != len(selected):
                raise ValueError("Explicit UZF cells must be unique.")
            self._validate_explicit_cells(selected)
            return selected

        domain = self.domain
        assert domain is not None
        selected: list[CellId] = []
        for cell in range(domain.shape[1]):
            active_layers = np.flatnonzero(domain[:, cell] > 0).tolist()
            if not active_layers:
                continue

            top = active_layers[0]
            selected.append((top, cell))
            if self.cells == "surface_only":
                continue

            for layer in range(top + 1, domain.shape[0]):
                if domain[layer, cell] <= 0:
                    break
                selected.append((layer, cell))
        return tuple(selected)

    @property
    def surface_cells(self) -> tuple[CellId, ...]:
        """Return the uppermost UZF cell in each represented grid column."""

        top_by_column: dict[int, CellId] = {}
        for cellid in self.uzf_cells:
            layer, column = cellid
            current = top_by_column.get(column)
            if current is None or layer < current[0]:
                top_by_column[column] = cellid
        return tuple(cellid for cellid in self.uzf_cells if top_by_column[cellid[1]] == cellid)

    def _validate_explicit_cells(self, cells: tuple[CellId, ...]) -> None:
        """Assert each explicit ``(layer, cell)`` is in range and active in the domain (no-op if none)."""

        domain = self.domain
        if domain is None:
            return
        for layer, cell in cells:
            if layer < 0 or cell < 0 or layer >= domain.shape[0] or cell >= domain.shape[1]:
                raise ValueError(f"UZF cell {(layer, cell)} is outside ModelContext.domain.")
            if domain[layer, cell] <= 0:
                raise ValueError(f"UZF cell {(layer, cell)} is inactive.")

    def _vertical_connections(self) -> dict[CellId, int]:
        """Map each UZF cell to the feature number directly beneath it."""

        feature_number = {cellid: index for index, cellid in enumerate(self.uzf_cells)}
        connections: dict[CellId, int] = {}
        for layer, cell in self.uzf_cells:
            connections[(layer, cell)] = feature_number.get((layer + 1, cell), -1)
        return connections

    def _cell_values(self, value: Any, *, name: str) -> list[Any]:
        """Normalize one static property to all configured UZF cells."""

        cells = self.uzf_cells
        if _is_scalar(value):
            return [value] * len(cells)
        if isinstance(value, Mapping):
            missing = [cellid for cellid in cells if cellid not in value]
            if missing:
                raise ValueError(f"{name} is missing UZF cells: {missing}")
            return [value[cellid] for cellid in cells]

        values = np.asarray(value)
        if self.domain is not None and values.shape == self.domain.shape:
            return [values[cellid] for cellid in cells]
        if values.ndim == 1 and len(values) == len(cells):
            return values.tolist()
        raise ValueError(
            f"{name} must be a scalar, a cell mapping, an array shaped like domain, "
            f"or contain {len(cells)} UZF-cell values."
        )

    def _surface_values(self, value: Any, *, name: str, default: Any = 0.0) -> list[Any]:
        """Normalize one period property to land-surface UZF cells."""

        cells = self.surface_cells
        if value is None:
            return [default] * len(cells)
        if _is_scalar(value):
            return [value] * len(cells)
        if isinstance(value, Mapping):
            missing = [cellid for cellid in cells if cellid not in value]
            if missing:
                raise ValueError(f"{name} is missing surface UZF cells: {missing}")
            return [value[cellid] for cellid in cells]

        values = np.asarray(value)
        if self.domain is not None and values.shape == self.domain.shape:
            return [values[cellid] for cellid in cells]
        if values.ndim == 1 and len(values) == len(cells):
            return values.tolist()
        raise ValueError(
            f"{name} period values must be a scalar, a surface-cell mapping, "
            f"an array shaped like domain, or contain {len(cells)} surface values."
        )

    def _period_values(self, value: Any, *, name: str, default: Any = 0.0) -> dict[int, list[Any]]:
        """Normalize transient input to every stress period and surface cell."""

        if isinstance(value, Mapping) and all(isinstance(key, int) for key in value):
            bad_periods = sorted(set(value) - set(range(self.nper)))
            if bad_periods:
                raise ValueError(f"{name} contains periods outside the simulation: {bad_periods}")
            return {
                period: self._surface_values(value.get(period), name=name, default=default)
                for period in range(self.nper)
            }
        return {
            period: self._surface_values(value, name=name, default=default)
            for period in range(self.nper)
        }

    @property
    def packagedata(self) -> list[list[Any]]:
        """Return MF6 UZF packagedata records."""

        connections = self._vertical_connections()
        surface = set(self.surface_cells)
        static = {
            "surfdep": self._cell_values(self.surfdep, name="surfdep"),
            "vks": self._cell_values(self.vks, name="vks"),
            "thtr": self._cell_values(self.thtr, name="thtr"),
            "thts": self._cell_values(self.thts, name="thts"),
            "thti": self._cell_values(self.thti, name="thti"),
            "eps": self._cell_values(self.eps, name="eps"),
        }
        records = []
        for ifno, cellid in enumerate(self.uzf_cells):
            record = [
                ifno,
                cellid,
                int(cellid in surface),
                connections[cellid],
                static["surfdep"][ifno],
                static["vks"][ifno],
                static["thtr"][ifno],
                static["thts"][ifno],
                static["thti"][ifno],
                static["eps"][ifno],
            ]
            if self.boundnames:
                record.append(f"uzf_{cellid[0]}_{cellid[1]}")
            records.append(record)
        return records

    @property
    def perioddata(self) -> dict[int, list[list[Any]]]:
        """Return MF6 UZF perioddata for land-surface UZF cells."""

        inputs = {
            "finf": self._period_values(self.finf, name="finf"),
            "pet": self._period_values(self.pet, name="pet"),
            "extdp": self._period_values(self.extdp, name="extdp"),
            "extwc": self._period_values(self.extwc, name="extwc"),
            "ha": self._period_values(self.ha, name="ha"),
            "hroot": self._period_values(self.hroot, name="hroot"),
            "rootact": self._period_values(self.rootact, name="rootact"),
        }
        feature_number = {cellid: index for index, cellid in enumerate(self.uzf_cells)}
        return {
            period: [
                [
                    feature_number[cellid],
                    inputs["finf"][period][index],
                    inputs["pet"][period][index],
                    inputs["extdp"][period][index],
                    inputs["extwc"][period][index],
                    inputs["ha"][period][index],
                    inputs["hroot"][period][index],
                    inputs["rootact"][period][index],
                ]
                for index, cellid in enumerate(self.surface_cells)
            ]
            for period in range(self.nper)
        }

    def validate(self) -> None:
        """Raise a clear error when the configured UZF data are invalid."""

        if not self.uzf_cells:
            raise ValueError("UZFBuilder selected no cells.")

        packagedata = self.packagedata
        for record in packagedata:
            ifno, cellid = record[0], record[1]
            surfdep, vks, thtr, thts, thti, eps = record[4:10]
            if float(surfdep) < 0:
                raise ValueError(f"surfdep must be nonnegative for UZF cell {cellid}.")
            if float(vks) <= 0:
                raise ValueError(f"vks must be positive for UZF cell {cellid}.")
            if not float(thtr) < float(thts):
                raise ValueError(f"thtr must be less than thts for UZF cell {cellid}.")
            if not float(thtr) <= float(thti) <= float(thts):
                raise ValueError(f"thti must be between thtr and thts for UZF cell {cellid}.")
            if float(eps) <= 0:
                raise ValueError(f"eps must be positive for UZF cell {cellid}.")
            if record[3] >= len(packagedata):
                raise ValueError(f"UZF cell feature {ifno} has an invalid vertical connection.")

        for period, records in self.perioddata.items():
            for record in records:
                if isinstance(record[1], Real) and float(record[1]) < 0:
                    raise ValueError(f"finf must be nonnegative in period {period}.")
                if isinstance(record[2], Real) and float(record[2]) < 0:
                    raise ValueError(f"pet must be nonnegative in period {period}.")
                if isinstance(record[3], Real) and float(record[3]) < 0:
                    raise ValueError(f"extdp must be nonnegative in period {period}.")

    def build(self) -> PackageSpec:
        """Validate the stored configuration and return its package spec."""

        self.validate()
        simulate_et = (
            self.pet is not None or self.extdp is not None or self.extwc is not None
            if self.simulate_et is None
            else self.simulate_et
        )
        package = uzf_spec(
            self.packagedata,
            self.perioddata,
            name=self.name,
            mover=self.mover,
            simulate_et=simulate_et,
            boundnames=self.boundnames,
            linear_gwet=self.linear_gwet,
            square_gwet=self.square_gwet,
            simulate_gwseep=self.simulate_gwseep,
            unsat_etwc=self.unsat_etwc,
            unsat_etae=self.unsat_etae,
            ntrailwaves=self.ntrailwaves,
            nwavesets=self.nwavesets,
            **dict(self.options),
        )
        return package.with_metadata(
            builder="UZFBuilder",
            cell_selection=self.cells if isinstance(self.cells, str) else "explicit",
            uzf_cells=list(self.uzf_cells),
            surface_cells=list(self.surface_cells),
        )


__all__ = ["UZFBuilder"]
