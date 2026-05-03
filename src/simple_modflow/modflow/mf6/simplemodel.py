"""Convenience helpers for assembling small MF6 models from a compact config."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from simple_modflow.modflow.mf6.drn import DRN
from simple_modflow.modflow.mf6.simulation.base import SimulationBase
from simple_modflow.modflow.mf6.simulation.discretization import (
    DisuGrid,
    DisvGrid,
    TemporalDiscretization,
)
from simple_modflow.modflow.mf6.simulation.packages import (
    CHD,
    Drains,
    InitialConditions,
    KFlow,
    OutputControl,
    Recharge,
    Storage,
)

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor


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
    model: SimulationBase,
    config: "SimpleModelConfig",
) -> list[list]:
    """Build edge-drain stress-period rows from a :class:`SimpleModelConfig`."""

    boundary_cells = _resolve_boundary_cells(config)
    layer_offsets = _coerce_layer_offsets(config.drain_bottom_addition, config.nlay)
    drain_builder = DRN(model=model, vor=model.vor)
    stress_period_data = []

    for layer_idx, bottom_addition in enumerate(layer_offsets):
        stress_period_data.extend(
            drain_builder.get_drn_stress_period_data(
                cells=boundary_cells,
                bottom_addition=bottom_addition,
                conductance=config.boundary_conductance,
                disMf=config.grid_type,
                layer=layer_idx,
            )
        )

    return stress_period_data


def build_constant_head_stress_period_data(
    model: SimulationBase,
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
    """Compact configuration for the legacy quick-start simple-model workflow."""
    vor: "Vor"
    name: str = "simplemodel"
    mf_folder_path: Path = Path().home() / "mf6"
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
    per_dates: list | None = None
    idomain_path: Path | None = None
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
        self.mf_folder_path = Path(self.mf_folder_path)
        if self.grid_type == "disu" and self.nlay != 1:
            raise ValueError("DISU simple models currently support nlay=1 only")


def _configure_simple_model(
    model: SimulationBase,
    config: SimpleModelConfig,
) -> SimulationBase:
    top, bottom = _resolve_top_and_bottom(config)
    initial_heads = _resolve_initial_heads(config, bottom)

    model.simple_model_config = config
    model.vor = config.vor
    model.name = config.name
    model.nper = config.nper
    model.nlay = config.nlay
    model.top = top
    model.bottom = bottom
    model.grid_type = config.grid_type
    model.boundary_mode = config.boundary_mode
    model.boundary_cells = _resolve_boundary_cells(config)
    model.rch_dict = config.rch_dict
    model.initial_sat_thickness = config.initial_sat_thickness
    model.initial_heads = initial_heads
    model.iheads = initial_heads
    model.k_input = config.k
    model.k33_vert_input = config.k33_vert

    model.tdis = TemporalDiscretization(
        model=model,
        time_units=config.time_units,
        per_len=config.per_len,
        num_steps=config.num_steps,
        multiplier=config.multiplier,
    )

    if config.grid_type == "disu":
        model.disu = DisuGrid(vor=config.vor, model=model, top=top, bottom=bottom)
    else:
        model.disv = DisvGrid(
            vor=config.vor,
            model=model,
            top=top,
            bottom=bottom,
            nlay=config.nlay,
            idomain=config.idomain,
        )

    model.ic = InitialConditions(
        model=model,
        vor=config.vor,
        botm_cells=_flatten_bottom_cells(bottom, grid_type=config.grid_type),
        initial_sat_thickness=config.initial_sat_thickness,
        nlay=config.nlay,
        strt=initial_heads,
    )
    model.k = KFlow(
        model=model,
        k=config.k,
        k33_vert=config.k33_vert,
        perched=config.perched,
        save_specific_discharge=config.save_specific_discharge,
    )
    model.oc = OutputControl(
        model=model,
        save_record=config.output_save_record,
        print_record=config.output_print_record,
    )
    model.sto = Storage(
        model=model,
        specific_yield=config.specific_yield,
        specific_storage=config.specific_storage,
        sto_steady=config.sto_steady or _default_sto_steady(config),
        sto_transient=(
            config.sto_transient
            if config.sto_transient is not None
            else _default_sto_transient(config)
        ),
    )

    if config.boundary_mode == "drain":
        model.drain_stress_period_data = build_edge_drain_stress_period_data(model, config)
        model.drn = Drains(model=model, stress_period_data=model.drain_stress_period_data)
    elif config.boundary_mode == "chd":
        model.chd_stress_period_data = build_constant_head_stress_period_data(model, config, top)
        model.chd = CHD(model=model, stress_period_data=model.chd_stress_period_data)

    if config.rch_dict is not None:
        model.rch = Recharge(model=model, vor=config.vor, rch_dict=config.rch_dict)

    return model


def build_simple_model(config: SimpleModelConfig) -> SimulationBase:
    """Construct and return a runnable ``SimulationBase`` from a simple config."""

    model = SimulationBase(
        name=config.name,
        mf_folder_path=config.mf_folder_path,
        nper=config.nper,
        vor=config.vor,
        per_dates=config.per_dates,
        idomain_path=config.idomain_path,
        newton=config.newton,
        complexity=config.complexity,
    )
    return _configure_simple_model(model, config)


class SimpleModel(SimulationBase):
    """Config-driven convenience class for building a small MF6 model."""

    def __init__(
        self,
        config: SimpleModelConfig | None = None,
        **config_kwargs,
    ):
        if config is None:
            if not config_kwargs:
                raise ValueError("Provide either a SimpleModelConfig or keyword arguments for SimpleModelConfig")
            config = SimpleModelConfig(**config_kwargs)
        elif config_kwargs:
            raise ValueError("Pass either config or keyword arguments, not both")

        super().__init__(
            name=config.name,
            mf_folder_path=config.mf_folder_path,
            nper=config.nper,
            vor=config.vor,
            per_dates=config.per_dates,
            idomain_path=config.idomain_path,
            newton=config.newton,
            complexity=config.complexity,
        )
        _configure_simple_model(self, config)


__all__ = [
    "SimpleModel",
    "SimpleModelConfig",
    "build_constant_head_stress_period_data",
    "build_edge_drain_stress_period_data",
    "build_simple_model",
]
