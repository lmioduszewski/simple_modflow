"""Build a small, reusable MODFLOW 6 simulation specification."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import flopy

from myflopy.builders import build_ims
from myflopy.specs import ModelContext, ModelSpec, PackageSpec, SimulationSpec

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor


def _is_sequence(value: Any) -> bool:
    """True if ``value`` is a list or tuple (a per-cell/per-layer array, not a scalar)."""

    return isinstance(value, (list, tuple))


def _coerce_cell_values(values: Any, expected_len: int, name: str) -> list:
    """Broadcast a scalar to ``expected_len`` cells, or validate a sequence's length."""

    if _is_sequence(values):
        result = list(values)
        if len(result) != expected_len:
            raise ValueError(f"{name} must have {expected_len} values, got {len(result)}")
        return result
    return [values for _ in range(expected_len)]


def _coerce_layer_offsets(values: Any, nlay: int) -> list[float]:
    """Broadcast a scalar drain-bottom offset to ``nlay``, or validate a per-layer list."""

    if _is_sequence(values):
        result = list(values)
        if len(result) != nlay:
            raise ValueError(f"drain_bottom_addition must have {nlay} values, got {len(result)}")
        return result
    return [values for _ in range(nlay)]


def _surface_columns(vor: Vor) -> list:
    """The non-geometry surface columns of the grid's ``gdf_topbtm`` (empty if absent)."""

    if getattr(vor, "gdf_topbtm", None) is None:
        return []
    return [column for column in vor.gdf_topbtm.columns if column != "geometry"]


def _resolve_top_and_bottom(config: SimpleModelConfig) -> tuple[list, list | list[list]]:
    """Resolve model top and per-layer bottoms from the config or the grid surfaces.

    Falls back to the grid's ``gdf_topbtm`` columns when ``top``/``bottom`` are
    unset, and returns cell arrays sized to the grid (a flat bottom for DISU, one
    array per layer for DISV).
    """

    top = config.top
    bottom = config.bottom
    surface_columns = _surface_columns(config.vor)

    if top is None:
        if not surface_columns:
            raise ValueError("top must be provided when vor.gdf_topbtm is unavailable")
        top = config.vor.gdf_topbtm.loc[:, surface_columns[0]].to_list()

    if bottom is None:
        required_columns = config.nlay + 1
        if len(surface_columns) < required_columns:
            raise ValueError(
                "bottom must be provided when vor.gdf_topbtm does not include enough layer surfaces"
            )
        if config.grid_type == "disu":
            if config.nlay != 1:
                raise ValueError("DISU simple models currently support nlay=1 only")
            bottom = config.vor.gdf_topbtm.loc[:, surface_columns[1]].to_list()
        else:
            bottom = [
                config.vor.gdf_topbtm.loc[:, surface_columns[layer + 1]].to_list()
                for layer in range(config.nlay)
            ]

    top = _coerce_cell_values(top, config.vor.ncpl, "top")

    if config.grid_type == "disu":
        if config.nlay != 1:
            raise ValueError("DISU simple models currently support nlay=1 only")
        bottom = _coerce_cell_values(bottom, config.vor.ncpl, "bottom")
    else:
        if not _is_sequence(bottom):
            raise ValueError("bottom must be a list of layer arrays for DISV models")
        bottom_layers = list(bottom)
        if len(bottom_layers) != config.nlay:
            raise ValueError(f"bottom must include {config.nlay} layers, got {len(bottom_layers)}")
        bottom = [
            _coerce_cell_values(layer_values, config.vor.ncpl, f"bottom[{layer_idx}]")
            for layer_idx, layer_values in enumerate(bottom_layers)
        ]

    return top, bottom


def _flatten_bottom_cells(bottom: list | list[list], *, grid_type: str) -> list:
    """Flatten the bottom arrays to one list of cell elevations (per-layer arrays concatenated)."""

    if grid_type == "disu":
        return list(bottom)
    return [cell for layer_values in bottom for cell in layer_values]


def _resolve_initial_heads(config: SimpleModelConfig, bottom: list | list[list]):
    """The initial-head array: the config's if given, else each cell bottom + saturated thickness."""

    if config.initial_heads is not None:
        return config.initial_heads

    botm_cells = _flatten_bottom_cells(bottom, grid_type=config.grid_type)
    return [cell_elev + config.initial_sat_thickness for cell_elev in botm_cells]


def _resolve_boundary_cells(config: SimpleModelConfig) -> list[int]:
    """The perimeter boundary cells: the config's if given, else the grid's edge cells."""

    if config.boundary_cells is not None:
        return list(config.boundary_cells)
    return list(config.vor.get_grid_edge())


def _resolve_boundary_heads(
    config: SimpleModelConfig,
    cells: list[int],
    top: list[float],
    layer_idx: int,
) -> list[float]:
    """Boundary heads for one layer's ``cells``: cell tops by default, else the config value(s)."""

    boundary_head = config.boundary_head

    if boundary_head is None:
        return [top[cell] for cell in cells]

    if _is_sequence(boundary_head) and len(boundary_head) == config.nlay:
        layer_values = list(boundary_head)[layer_idx]
        if _is_sequence(layer_values):
            return _coerce_cell_values(layer_values, len(cells), f"boundary_head[{layer_idx}]")
        return [layer_values for _ in cells]

    return _coerce_cell_values(boundary_head, len(cells), "boundary_head")


def build_edge_drain_stress_period_data(
    config: SimpleModelConfig,
    bottom: list | list[list],
) -> list[list]:
    """Build edge-drain stress-period rows from a :class:`SimpleModelConfig`.

    One row per boundary cell per layer, ``[cellid, elevation, conductance]``:
    the cellid is ``(layer, cell)`` on DISV and the bare cell on DISU, and the
    elevation is that cell's bottom plus ``config.drain_bottom_addition``.

    Parameters
    ----------
    config : SimpleModelConfig
        Supplies ``boundary_cells`` (default the grid's perimeter),
        ``boundary_conductance``, ``drain_bottom_addition``, ``nlay`` and
        ``grid_type``. ``boundary_mode`` is not read.
    bottom : list
        Resolved bottom elevations: for DISV a list of ``nlay`` per-cell lists,
        for DISU one per-cell list (the form :func:`simple_model_spec` resolves).

    Returns
    -------
    list of list
        The rows, one period's worth, layer by layer.
    """

    boundary_cells = _resolve_boundary_cells(config)
    layer_offsets = _coerce_layer_offsets(config.drain_bottom_addition, config.nlay)
    conductances = _coerce_cell_values(
        config.boundary_conductance,
        len(boundary_cells),
        "boundary_conductance",
    )
    stress_period_data = []

    for layer_idx, bottom_addition in enumerate(layer_offsets):
        layer_bottom = bottom if config.grid_type == "disu" else bottom[layer_idx]
        for cell, conductance in zip(boundary_cells, conductances, strict=False):
            cell_id = cell if config.grid_type == "disu" else (layer_idx, cell)
            stress_period_data.append(
                [cell_id, layer_bottom[cell] + bottom_addition, conductance]
            )

    return stress_period_data


def build_constant_head_stress_period_data(
    config: SimpleModelConfig,
    top: list[float],
) -> dict[int, list[list]]:
    """Build constant-head stress-period data from a :class:`SimpleModelConfig`.

    Parameters
    ----------
    config : SimpleModelConfig
        Supplies ``boundary_cells`` (default the grid's perimeter),
        ``boundary_head``, ``nlay``, ``nper`` and ``grid_type``.
        ``boundary_mode`` is not read.
    top : list of float
        Resolved model top, one value per cell of a layer; it is the head used
        when ``config.boundary_head`` is None.

    Returns
    -------
    dict of {int: list of list}
        ``{period: [[cellid, head], ...]}`` for every zero-based period, one row
        per boundary cell per layer; the cellid is ``(layer, cell)`` on DISV and
        the bare cell on DISU.
    """

    boundary_cells = _resolve_boundary_cells(config)
    stress_period_data: dict[int, list[list]] = {}

    for per in range(config.nper):
        period_rows = []
        for layer_idx in range(config.nlay):
            layer_heads = _resolve_boundary_heads(config, boundary_cells, top, layer_idx)
            for cell, head in zip(boundary_cells, layer_heads, strict=False):
                cell_id = cell if config.grid_type == "disu" else (layer_idx, cell)
                period_rows.append([cell_id, head])
        stress_period_data[per] = period_rows

    return stress_period_data


def _default_sto_steady(config: SimpleModelConfig) -> dict[int, bool]:
    """Default STO steady-state flags: only the first stress period is steady."""

    return {0: True}


def _default_sto_transient(config: SimpleModelConfig) -> dict[int, bool]:
    """Default STO transient flags: every period after the first (empty for a single period)."""

    return {} if config.nper <= 1 else {per: True for per in range(1, config.nper)}


@dataclass(slots=True)
class SimpleModelConfig:
    """Compact, all-in-one description of a basic single-GWF MF6 model.

    Gathers everything :func:`simple_model_spec` needs to assemble a runnable
    flow model in one flat dataclass: the Voronoi grid (``vor``), discretization
    (``nlay``/``grid_type``/``top``/``bottom``/``idomain``), time stepping
    (``nper``/``per_len``/``num_steps``/...), aquifer properties (``k``/``k33_vert``/
    storage), initial heads, and an optional perimeter boundary (CHD or DRN via
    ``boundary_mode``/``boundary_cells``/...). Most fields have sensible defaults,
    so a minimal config is just a grid plus a name. This is a convenience
    front-end; for full multi-package or multi-physics models use the
    package-first API (``mf.gwf(...)`` with explicit package specs).

    The model's length unit is always feet (``length_units="FEET"`` on the
    discretization) and its time unit is ``time_units``. Every cell is
    convertible (NPF ``icelltype=1``, STO ``iconvert=1``). See
    :func:`simple_model_spec` for how the fields become a
    :class:`~myflopy.specs.SimulationSpec`.

    ``top``, ``bottom``, ``boundary_conductance``, ``boundary_head`` and
    ``drain_bottom_addition`` must be a scalar or a Python ``list`` / ``tuple``:
    there a NumPy array is not recognised as a sequence, so it is copied whole
    into every cell (measured: an ``ndarray`` ``top`` fails inside FloPy's DISV,
    and an ``ndarray`` ``bottom`` raises ``ValueError``). Call ``.tolist()``
    first. The fields passed to FloPy unchanged (``k``, ``k33_vert``,
    ``initial_heads``, ``idomain``, ``specific_storage``, ``specific_yield``)
    take an array as it is (measured for ``k`` and ``initial_heads``).

    Parameters
    ----------
    vor : VoronoiGridPlus
        The grid. Supplies the cells and vertices, the perimeter cells
        (``get_grid_edge()``) and, when ``top``/``bottom`` are omitted, the
        surfaces in ``vor.gdf_topbtm``; for ``grid_type="disu"`` also the
        connection arrays (``iac``, ``ja``, ``cl12``, ``hwva``, cell areas).
    name : str, optional
        Simulation and model name, and the stem of every file written
        (``<name>.disv``, ``<name>.hds``, ...). At most 16 characters, MODFLOW
        6's limit; longer raises ``ValueError``. Default ``"simplemodel"``.
    nper : int, optional
        Number of stress periods, at least 1 (``ValueError`` otherwise).
        Default 1.
    nlay : int, optional
        Number of layers, at least 1; must be 1 for ``grid_type="disu"``
        (``ValueError`` otherwise). Default 1.
    grid_type : str, optional
        ``"disv"`` (default) or ``"disu"``, case-insensitive (stored lower-case).
        Anything else raises ``ValueError``.
    top : float, list or tuple, optional
        Model top elevation in feet: one value for every cell, or a list/tuple
        of ``ncpl`` values (another length raises ``ValueError``). Default None:
        the first surface column of ``vor.gdf_topbtm``; ``ValueError`` if the
        grid has none.
    bottom : float, list or tuple, optional
        Layer bottom elevations in feet. For DISV, a list/tuple of ``nlay``
        entries, top layer first, each a scalar or a list/tuple of ``ncpl``
        values (a bare scalar raises ``ValueError``). For DISU, a scalar or a
        list/tuple of ``ncpl`` values. Default None: the ``nlay`` surface
        columns after the top in ``vor.gdf_topbtm``; ``ValueError`` if it has
        fewer than ``nlay + 1``.
    idomain : array_like, optional
        Passed to the discretization's ``idomain`` and kept as the model
        context's ``domain``. Default None: for DISV none is written (MF6 then
        treats every cell as active); for DISU a list of ``ncpl`` ones.
    initial_heads : float or array_like, optional
        Starting heads, passed to IC ``strt`` unchanged. Default None: each
        cell's bottom plus ``initial_sat_thickness``, listed layer by layer.
    initial_sat_thickness : float, optional
        Feet of saturation above each cell bottom for the default starting
        heads; ignored when ``initial_heads`` is given. Default 1.0.
    per_len : float, optional
        Length of EVERY stress period, in ``time_units``. Default 30.
    num_steps : int, optional
        Time steps in every stress period. Default 1.
    multiplier : float, optional
        Time-step multiplier (TDIS ``tsmult``) in every stress period.
        Default 1.0.
    time_units : str, optional
        TDIS ``time_units`` (``"DAYS"``, ``"SECONDS"``, ...). Default ``"DAYS"``.
    k : float or array_like, optional
        Horizontal hydraulic conductivity in feet per ``time_units``, passed to
        NPF ``k`` unchanged. Default 100.
    k33_vert : float or array_like, optional
        Vertical hydraulic conductivity, passed to NPF ``k33`` unchanged.
        Default None: not written, so MODFLOW 6 sets it equal to ``k``.
    perched : bool, optional
        NPF ``PERCHED`` option. Default False.
    save_specific_discharge : bool, optional
        NPF ``SAVE_SPECIFIC_DISCHARGE`` option. Default True. (``SAVE_FLOWS`` and
        ``SAVE_SATURATION`` are always on.)
    specific_storage : float or array_like, optional
        STO ``ss``, per foot. Default 0.0001.
    specific_yield : float or array_like, optional
        STO ``sy``, dimensionless. Default 0.2.
    sto_steady : dict of {int: bool}, optional
        STO ``steady_state``, keyed by zero-based stress period. Default None --
        and also an EMPTY dict -- gives ``{0: True}``.
    sto_transient : dict of {int: bool}, optional
        STO ``transient``, keyed by zero-based stress period. Default None:
        every period after the first (``{}`` for a single period). Filled
        independently of ``sto_steady``.
    boundary_mode : str, optional
        Perimeter boundary package: ``"drain"`` (DRN) or ``"chd"`` (CHD), placed
        on ``boundary_cells`` in every layer. Default None: none. The match is
        exact and unchecked -- ``"ghb"``, ``"drn"`` or ``"CHD"`` build NO
        boundary package and raise nothing (measured).
    boundary_cells : list of int, optional
        Zero-based cell indices within a layer for the boundary. Default None:
        the grid's perimeter cells, ``vor.get_grid_edge()``.
    boundary_conductance : float, list or tuple, optional
        DRN conductance (square feet per ``time_units``): one value, or one per
        boundary cell, the same in every layer. Used only when
        ``boundary_mode="drain"``. Default 1000.
    boundary_head : float, list or tuple, optional
        CHD head in feet, used only when ``boundary_mode="chd"``, and applied in
        every stress period. A list/tuple of exactly ``nlay`` entries is read
        per layer (each a scalar or one head per boundary cell); anything else
        is one value, or one head per boundary cell, for every layer. Default
        None: each boundary cell's model top.
    drain_bottom_addition : float, list or tuple, optional
        DRN elevation above each boundary cell's bottom, in feet: one value, or
        one per layer (``nlay`` entries). Used only when
        ``boundary_mode="drain"``. The drain records are given for period 0
        only, so MF6 repeats them in later periods. Default 0.1.
    rch_dict : dict, optional
        RCH ``stress_period_data`` (``{period: [(cellid, rate), ...]}``, rate in
        feet per ``time_units``), passed unchanged with ``maxbound`` set to the
        cell count of one layer. Default None: no RCH package.
    newton : bool, optional
        True writes the model's ``NEWTON UNDER_RELAXATION`` option; False writes
        no Newton option. Default True.
    complexity : str, optional
        IMS ``complexity`` (``"SIMPLE"``, ``"MODERATE"`` or ``"COMPLEX"``), passed
        to :func:`~myflopy.builders.build_ims`. Default ``"MODERATE"``.
    output_save_record : tuple, optional
        OC ``saverecord`` (heads go to ``<name>.hds``, budgets to
        ``<name>.cbc``). Default ``(("HEAD", "LAST"), ("BUDGET", "LAST"))``.
    output_print_record : tuple, optional
        OC ``printrecord``. Default None: the OC file carries no ``PRINT``
        records.
    """
    vor: Vor
    name: str = "simplemodel"
    nper: int = 1
    nlay: int = 1
    grid_type: str = "disv"
    top: Any = None
    bottom: Any = None
    idomain: Any = None
    initial_heads: Any = None
    initial_sat_thickness: float = 1.0
    per_len: int = 30
    num_steps: int = 1
    multiplier: float = 1.0
    time_units: str = "DAYS"
    k: Any = 100
    k33_vert: Any = None
    perched: bool = False
    save_specific_discharge: bool = True
    specific_storage: float = 0.0001
    specific_yield: float = 0.2
    sto_steady: dict[int, bool] | None = None
    sto_transient: dict[int, bool] | None = None
    boundary_mode: str | None = None
    boundary_cells: list[int] | None = None
    boundary_conductance: float | int | list = 1000
    boundary_head: Any = None
    drain_bottom_addition: float | int | list = 0.1
    rch_dict: dict | None = None
    newton: bool = True
    complexity: str = "MODERATE"
    output_save_record: tuple = (("HEAD", "LAST"), ("BUDGET", "LAST"))
    output_print_record: tuple | None = None

    def __post_init__(self):
        """Normalize ``grid_type`` and validate name length, ``nlay``/``nper``, and DISU limits."""

        self.grid_type = str(self.grid_type).lower()
        if self.grid_type not in {"disu", "disv"}:
            raise ValueError("grid_type must be either 'disu' or 'disv'")
        if len(self.name) > 16:
            raise ValueError("name must be 16 characters or fewer for MODFLOW 6")
        if self.nlay < 1:
            raise ValueError("nlay must be at least 1")
        if self.nper < 1:
            raise ValueError("nper must be at least 1")
        if self.grid_type == "disu" and self.nlay != 1:
            raise ValueError("DISU simple models currently support nlay=1 only")


def _discretization_spec(
    config: SimpleModelConfig,
    top: list,
    bottom: list | list[list],
) -> PackageSpec:
    """Return the configured DISV or DISU package specification."""

    grid_props = config.vor.get_disv_gridprops()
    if config.grid_type == "disv":
        return PackageSpec(
            "disv",
            flopy.mf6.ModflowGwfdisv,
            {
                "length_units": "FEET",
                "nlay": config.nlay,
                "ncpl": grid_props["ncpl"],
                "nvert": len(grid_props["vertices"]),
                "vertices": grid_props["vertices"],
                "cell2d": grid_props["cell2d"],
                "top": top,
                "botm": bottom,
                "idomain": config.idomain,
                "pname": "disv",
                "filename": f"{config.name}.disv",
            },
        )

    return PackageSpec(
        "disu",
        flopy.mf6.ModflowGwfdisu,
        {
            "length_units": "FEET",
            "vertices": grid_props["vertices"],
            "cell2d": grid_props["cell2d"],
            "top": top,
            "bot": bottom,
            "nvert": len(grid_props["vertices"]),
            "nodes": len(grid_props["cell2d"]),
            "nja": config.vor.nja,
            "iac": config.vor.iac,
            "ja": config.vor.ja,
            "area": config.vor.get_cell_areas(),
            "ihc": 1,
            "cl12": config.vor.cl12,
            "hwva": config.vor.hwva,
            "idomain": (
                config.idomain
                if config.idomain is not None
                else [1 for _ in range(config.vor.ncpl)]
            ),
            "pname": "disu",
            "filename": f"{config.name}.disu",
        },
    )


def _flow_package_specs(
    config: SimpleModelConfig,
    top: list,
    bottom: list | list[list],
    initial_heads: Any,
) -> tuple[PackageSpec, ...]:
    """Return the ordered GWF packages selected by the simple-model config."""

    packages = [
        _discretization_spec(config, top, bottom),
        PackageSpec(
            "ic",
            flopy.mf6.ModflowGwfic,
            {
                "strt": initial_heads,
                "pname": "ic",
                "filename": f"{config.name}.ic",
            },
        ),
        PackageSpec(
            "npf",
            flopy.mf6.ModflowGwfnpf,
            {
                "k": config.k,
                "k33": config.k33_vert,
                "icelltype": 1,
                "perched": config.perched,
                "save_flows": True,
                "save_saturation": True,
                "save_specific_discharge": config.save_specific_discharge,
                "pname": "npf",
                "filename": f"{config.name}.npf",
            },
        ),
        PackageSpec(
            "sto",
            flopy.mf6.ModflowGwfsto,
            {
                "ss": config.specific_storage,
                "sy": config.specific_yield,
                "iconvert": 1,
                "steady_state": config.sto_steady or _default_sto_steady(config),
                "transient": (
                    config.sto_transient
                    if config.sto_transient is not None
                    else _default_sto_transient(config)
                ),
                "save_flows": True,
                "pname": "sto",
                "filename": f"{config.name}.sto",
            },
        ),
        PackageSpec(
            "oc",
            flopy.mf6.ModflowGwfoc,
            {
                "head_filerecord": f"{config.name}.hds",
                "budget_filerecord": f"{config.name}.cbc",
                "saverecord": config.output_save_record,
                "printrecord": config.output_print_record,
                "pname": "oc",
                "filename": f"{config.name}.oc",
            },
        ),
    ]

    if config.boundary_mode == "drain":
        packages.append(
            PackageSpec(
                "drn",
                flopy.mf6.ModflowGwfdrn,
                {
                    "stress_period_data": build_edge_drain_stress_period_data(config, bottom),
                    "save_flows": True,
                    "pname": "drn",
                    "filename": f"{config.name}.drn",
                },
            )
        )
    elif config.boundary_mode == "chd":
        packages.append(
            PackageSpec(
                "chd",
                flopy.mf6.ModflowGwfchd,
                {
                    "stress_period_data": build_constant_head_stress_period_data(config, top),
                    "save_flows": True,
                    "pname": "chd",
                    "filename": f"{config.name}.chd",
                },
            )
        )

    if config.rch_dict is not None:
        packages.append(
            PackageSpec(
                "rch",
                flopy.mf6.ModflowGwfrch,
                {
                    "stress_period_data": config.rch_dict,
                    "maxbound": len(config.vor.iverts),
                    "save_flows": True,
                    "pname": "rch",
                    "filename": f"{config.name}.rch",
                },
            )
        )
    return tuple(packages)


def simple_model_spec(config: SimpleModelConfig) -> SimulationSpec:
    """Build a runnable single-GWF :class:`SimulationSpec` from a config.

    The convenience bridge from a flat :class:`SimpleModelConfig` to the
    declarative spec API: it resolves top/botm and initial heads, assembles the
    DIS(V)/NPF/IC/STO/OC packages plus any perimeter boundary, wraps them in a
    GWF :class:`~myflopy.specs.ModelSpec` with the grid attached via
    :class:`~myflopy.specs.ModelContext`, and returns a
    :class:`~myflopy.specs.SimulationSpec` (TDIS + IMS + the model) ready to
    build/run. Use it for quick single-model setups; reach for ``mf.gwf(...)``
    with explicit packages when you need full control or multiple packages/models.

    Parameters
    ----------
    config : SimpleModelConfig
        The compact model description. Its grid is attached as the model
        context's ``grid``, the resolved top and bottoms as its ``surfaces``,
        and ``idomain`` as its ``domain``. Nothing is written or run here.

    Returns
    -------
    SimulationSpec
        A composable, runnable simulation spec named ``config.name``: TDIS and
        IMS plus one GWF model with DISV (or DISU), IC, NPF, STO and OC, then
        DRN or CHD when ``boundary_mode`` asks for one and RCH when
        ``rch_dict`` is given.

    Raises
    ------
    ValueError
        ``top`` or ``bottom`` is omitted and ``vor.gdf_topbtm`` cannot supply
        it, or a per-cell or per-layer list has the wrong length.
    """

    top, bottom = _resolve_top_and_bottom(config)
    initial_heads = _resolve_initial_heads(config, bottom)
    flow = ModelSpec(
        config.name,
        "gwf",
        packages=_flow_package_specs(config, top, bottom, initial_heads),
        context=ModelContext(
            grid=config.vor,
            surfaces={"top": top, "bottom": bottom},
            domain=config.idomain,
            metadata={
                "grid_type": config.grid_type,
                "nlay": config.nlay,
                "nper": config.nper,
            },
        ),
        options={
            "model_nam_file": f"{config.name}.nam",
            "newtonoptions": "under_relaxation" if config.newton else None,
            "save_flows": True,
        },
    )
    perioddata = [
        (config.per_len, config.num_steps, config.multiplier)
        for _ in range(config.nper)
    ]
    return SimulationSpec(
        config.name,
        models=(flow,),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {
                    "time_units": config.time_units,
                    "nper": config.nper,
                    "perioddata": perioddata,
                    "pname": "tdis",
                    "filename": f"{config.name}.tdis",
                },
            ),
            PackageSpec(
                "ims",
                build_ims,
                {
                    "models": (config.name,),
                    "complexity": config.complexity,
                    "pname": "ims_gwf",
                    "print_option": "SUMMARY",
                },
            ),
        ),
        options={"exe_name": "mf6", "version": "mf6"},
    )


__all__ = [
    "SimpleModelConfig",
    "build_constant_head_stress_period_data",
    "build_edge_drain_stress_period_data",
    "simple_model_spec",
]
