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
    return isinstance(value, (list, tuple))


def _coerce_cell_values(values: Any, expected_len: int, name: str) -> list:
    if _is_sequence(values):
        result = list(values)
        if len(result) != expected_len:
            raise ValueError(f"{name} must have {expected_len} values, got {len(result)}")
        return result
    return [values for _ in range(expected_len)]


def _coerce_layer_offsets(values: Any, nlay: int) -> list[float]:
    if _is_sequence(values):
        result = list(values)
        if len(result) != nlay:
            raise ValueError(f"drain_bottom_addition must have {nlay} values, got {len(result)}")
        return result
    return [values for _ in range(nlay)]


def _surface_columns(vor: "Vor") -> list:
    if getattr(vor, "gdf_topbtm", None) is None:
        return []
    return [column for column in vor.gdf_topbtm.columns if column != "geometry"]


def _resolve_top_and_bottom(config: "SimpleModelConfig") -> tuple[list, list | list[list]]:
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
    if grid_type == "disu":
        return list(bottom)
    return [cell for layer_values in bottom for cell in layer_values]


def _resolve_initial_heads(config: "SimpleModelConfig", bottom: list | list[list]):
    if config.initial_heads is not None:
        return config.initial_heads

    botm_cells = _flatten_bottom_cells(bottom, grid_type=config.grid_type)
    return [cell_elev + config.initial_sat_thickness for cell_elev in botm_cells]


def _resolve_boundary_cells(config: "SimpleModelConfig") -> list[int]:
    if config.boundary_cells is not None:
        return list(config.boundary_cells)
    return list(config.vor.get_grid_edge())


def _resolve_boundary_heads(
    config: "SimpleModelConfig",
    cells: list[int],
    top: list[float],
    layer_idx: int,
) -> list[float]:
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
    config: "SimpleModelConfig",
    bottom: list | list[list],
) -> list[list]:
    """Build edge-drain stress-period rows from a :class:`SimpleModelConfig`."""

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
        for cell, conductance in zip(boundary_cells, conductances):
            cell_id = cell if config.grid_type == "disu" else (layer_idx, cell)
            stress_period_data.append(
                [cell_id, layer_bottom[cell] + bottom_addition, conductance]
            )

    return stress_period_data


def build_constant_head_stress_period_data(
    config: "SimpleModelConfig",
    top: list[float],
) -> dict[int, list[list]]:
    """Build constant-head stress-period data from a :class:`SimpleModelConfig`."""

    boundary_cells = _resolve_boundary_cells(config)
    stress_period_data: dict[int, list[list]] = {}

    for per in range(config.nper):
        period_rows = []
        for layer_idx in range(config.nlay):
            layer_heads = _resolve_boundary_heads(config, boundary_cells, top, layer_idx)
            for cell, head in zip(boundary_cells, layer_heads):
                cell_id = cell if config.grid_type == "disu" else (layer_idx, cell)
                period_rows.append([cell_id, head])
        stress_period_data[per] = period_rows

    return stress_period_data


def _default_sto_steady(config: "SimpleModelConfig") -> dict[int, bool]:
    return {0: True}


def _default_sto_transient(config: "SimpleModelConfig") -> dict[int, bool]:
    return {} if config.nper <= 1 else {per: True for per in range(1, config.nper)}


@dataclass(slots=True)
class SimpleModelConfig:
    """Compact inputs for :func:`simple_model_spec`."""
    vor: "Vor"
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
    """Translate a compact config into a composable, runnable simulation spec."""

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
