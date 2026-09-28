"""Build MODFLOW 6 UZF package specifications."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Real
from typing import Any

import numpy as np

from myflopy.advanced import uzf_spec
from myflopy.grid_identity import CELLWISE_OPTIONS, pin_metadata
from myflopy.modflow.mf6.grid.geometry import top_active_layer
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

    This is the engine under ``mf.uzf(...)``, declared with or without
    ``context=``. The builder is frozen and resolves lazily: cells and records
    are computed from the stored fields each time ``.packagedata``,
    ``.perioddata`` or :meth:`build` asks for them (nothing is cached), and
    :meth:`build` runs :meth:`validate` first. Lengths are model length units and times model time
    units throughout. A UZF cell's number (MF6 ``IFNO``) is its position in
    :attr:`uzf_cells`.

    Parameters
    ----------
    context : ModelContext
        Supplies ``context.domain`` (the IDOMAIN, ``(nlay, ncpl)`` or a
        one-dimensional single layer), which picks the cells when ``cells`` is a
        keyword -- then it is required, ``ValueError`` at construction -- and
        checks explicit cells when present. The grid is not read except to pin
        the built package to it; without one the package is marked as
        addressing cells instead.
    nper : int
        Number of stress periods; ``perioddata`` gets every period
        ``0 .. nper - 1``. Less than 1 raises ``ValueError``.
    vks : float or array_like or Mapping
        Vertical saturated hydraulic conductivity of the unsaturated zone, MF6
        ``VKS`` (model length per model time), per UZF cell -- any of the
        per-cell shapes in the Notes. :meth:`validate` requires it positive.
    thtr : float or array_like or Mapping
        Residual water content, MF6 ``THTR`` (dimensionless), per UZF cell.
        :meth:`validate` requires ``thtr < thts``.
    thts : float or array_like or Mapping
        Saturated water content, MF6 ``THTS`` (dimensionless), per UZF cell.
    thti : float or array_like or Mapping
        Initial water content, MF6 ``THTI`` (dimensionless), per UZF cell.
        :meth:`validate` requires ``thtr <= thti <= thts``.
    cells : str or sequence of (int, int), optional
        Which cells get a UZF object. ``"all_active"`` (default): in each column,
        the top cell with ``idomain > 0`` and each cell below it, stopping
        before the first with ``idomain <= 0``. ``"surface_only"``: just that top
        cell. Or explicit zero-based ``(layer, cell)`` ids, in the order that
        numbers them; they must be unique, and in range and active when the
        context has a domain (``ValueError`` otherwise). Any other string raises
        ``ValueError`` at construction. In every form the uppermost UZF cell of a
        column is its land-surface cell, and each cell connects vertically to
        the UZF cell directly below it when there is one.
    eps : float or array_like or Mapping, optional
        Brooks-Corey epsilon exponent, MF6 ``EPS``, per UZF cell; default 4.0.
        :meth:`validate` requires it positive.
    surfdep : float or array_like or Mapping, optional
        Surface depression depth, MF6 ``SURFDEP`` (model length units), per UZF
        cell; default 0.001. :meth:`validate` refuses a negative value.
    finf : float, str, array_like or Mapping, optional
        Applied infiltration rate, MF6 ``FINF`` (model length per model time),
        on the land-surface cells -- any of the period shapes in the Notes.
        Default 0.0. :meth:`validate` refuses a negative number.
    pet : float, str, array_like, Mapping or None, optional
        Potential evapotranspiration rate, MF6 ``PET`` (model length per model
        time), on the land-surface cells; used by MF6 only with ET simulated.
        Default ``None``, written as 0.0. :meth:`validate` refuses a negative
        number.
    extdp : float, str, array_like, Mapping or None, optional
        Evapotranspiration extinction depth, MF6 ``EXTDP`` (model length units),
        on the land-surface cells. Default ``None``, written as 0.0.
        :meth:`validate` refuses a negative number.
    extwc : float, str, array_like, Mapping or None, optional
        Extinction water content, MF6 ``EXTWC`` (dimensionless), on the
        land-surface cells; MF6 uses it only with ``unsat_etwc``. Default
        ``None``, written as 0.0.
    ha : float, str, array_like, Mapping or None, optional
        Air-entry potential (head), MF6 ``HA`` (model length units), on the
        land-surface cells; MF6 uses it only with ``unsat_etae``. Default
        ``None``, written as 0.0.
    hroot : float, str, array_like, Mapping or None, optional
        Root potential (head), MF6 ``HROOT`` (model length units), on the
        land-surface cells; MF6 uses it only with ``unsat_etae``. Default
        ``None``, written as 0.0.
    rootact : float, str, array_like, Mapping or None, optional
        Root activity function, MF6 ``ROOTACT``, on the land-surface cells; MF6
        uses it only with ``unsat_etae``. Default ``None``, written as 0.0.
    name : str, optional
        Package name, default ``"uzf"``: FloPy's ``pname``, the input file
        ``<model>.<name>`` and the output files ``<model>_<name>_*``.
    boundnames : bool, optional
        Write MF6's ``BOUNDNAMES`` option and name every UZF cell
        ``"uzf_<layer>_<cell>"`` (zero-based). Default ``False``.
    mover : bool, optional
        Write MF6's ``MOVER`` option, which the package needs before an MVR
        package (``mf.mvr``) can move water to or from it. Default ``False``.
    simulate_et : bool or None, optional
        MF6 ``SIMULATE_ET``. Default ``None``: on when any of ``pet``, ``extdp``
        or ``extwc`` is given (``ha``/``hroot``/``rootact`` do not turn it on).
    linear_gwet : bool, optional
        MF6 ``LINEAR_GWET``: groundwater ET by MODFLOW-2005's original ET
        formulation. Default ``False``.
    square_gwet : bool, optional
        MF6 ``SQUARE_GWET``: groundwater ET at a constant rate for heads between
        land surface and ``extdp`` below it, reduced smoothly to zero there.
        Default ``False``.
    simulate_gwseep : bool, optional
        MF6 ``SIMULATE_GWSEEP``: simulate groundwater discharge to land surface.
        Default ``False``.
    unsat_etwc : bool, optional
        MF6 ``UNSAT_ETWC``: unsaturated-zone ET limited by ``extwc``. Default
        ``False``.
    unsat_etae : bool, optional
        MF6 ``UNSAT_ETAE``: unsaturated-zone ET by the capillary-pressure
        (``ha``/``hroot``/``rootact``) formulation. Default ``False``.
    ntrailwaves : int, optional
        MF6 ``NTRAILWAVES``, the number of trailing waves. Default 7.
    nwavesets : int, optional
        MF6 ``NWAVESETS``, the number of wave sets. Default 40.
    options : Mapping[str, Any], optional
        Extra keyword arguments for :class:`flopy.mf6.ModflowGwfuzf`, passed
        through :func:`myflopy.advanced.uzf_spec`; they override the defaults
        that function sets (``save_flows=True``, the output file records). A key
        :meth:`build` already passes -- ``packagedata``, ``perioddata``,
        ``name``, ``mover``, ``simulate_et``, ``boundnames``, the five ET/seepage
        flags from ``linear_gwet`` to ``unsat_etae``, ``ntrailwaves`` or
        ``nwavesets`` -- raises ``TypeError`` at build. Default empty.

    Notes
    -----
    **Per-cell inputs** (``vks``, ``thtr``, ``thts``, ``thti``, ``eps``,
    ``surfdep``) take a number (every UZF cell); a mapping ``{(layer, cell):
    value}`` covering every UZF cell (a missing one raises ``ValueError``); an
    array shaped like the domain (``(nlay, ncpl)``; a one-dimensional domain
    counts as ``(1, ncpl)``), read at each cell; or a flat sequence of one
    value per UZF cell, in :attr:`uzf_cells` order. Any other shape raises
    ``ValueError``. A string is written
    unchanged, and :meth:`validate` then fails unless it reads as a number.

    **Period inputs** (``finf``, ``pet``, ``extdp``, ``extwc``, ``ha``,
    ``hroot``, ``rootact``) are written for the land-surface cells only. A
    mapping whose keys are ALL integers is ``{period: value}``: a key outside
    ``0 .. nper - 1`` raises ``ValueError``, and a period left out is written as
    0.0 -- it does NOT carry the previous period forward. Any other value is
    used for every period. Each period's value takes the per-cell shapes above
    against the land-surface cells (a flat sequence is one value per
    :attr:`surface_cells` entry), and a string, such as a time-series name, is
    written unchanged.
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
        """Return a changed builder without modifying the original.

        Parameters
        ----------
        **updates
            Constructor fields to replace, by name (``finf=1e-4``,
            ``cells="surface_only"``), passed to :func:`dataclasses.replace`.
            Every other field is carried over. The construction checks run
            again, and a name that is not a constructor field raises
            ``TypeError``.

        Returns
        -------
        UZFBuilder
        """

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
        # The column's land surface is its top ACTIVE cell -- the shared rule,
        # `idomain > 0` (ledger 319).
        for cell, top in enumerate(top_active_layer(domain).tolist()):
            if top < 0:
                continue

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
        # A PIN when there was a mesh to pin to, and the weaker cell-addressing MARKER
        # otherwise -- the two-branch shape `_array_bc_spec` already uses (ledger 282).
        #
        # `pin_metadata` returns {} for a context carrying no grid, and every UZF record
        # holds a cellid at index 1 by construction, so without the fallback a GRIDLESS
        # declaration reached `package_move` with nothing at all to classify on and
        # answered "unchanged -- names no cells". Measured: a 25 -> 81 cell `with_grid`
        # then BUILT CLEAN with 75 UZF objects sitting on columns 0..24 of 81, no
        # warning -- UZF over 31% of the domain. Three separate misses had to line up for
        # that: no pin, `_cellids_in` deliberately walks only `stress_period_data`, and
        # the ncpl-exact fallback reads `None` from a record list whose second element is
        # a cellid tuple. Ledger 295's defect, one package over (ledger 302).
        #
        # UNCONDITIONAL, unlike `_array_bc_spec`'s `if not cellwise: return spec` -- that
        # one has a genuinely grid-free payload to protect (`mf.rch.array(irch=0)` means
        # the same layer on every mesh). A UZF package with no cell-addressing records
        # does not exist.
        stamped = pin_metadata(self.context, certifies=("packagedata",))
        return package.with_metadata(
            **(stamped or {CELLWISE_OPTIONS: ["packagedata"]}),
            builder="UZFBuilder",
            cell_selection=self.cells if isinstance(self.cells, str) else "explicit",
            uzf_cells=list(self.uzf_cells),
            surface_cells=list(self.surface_cells),
        )


__all__ = ["UZFBuilder"]
