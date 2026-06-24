"""Build MODFLOW 6 LAK package specifications from lake geometries."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Real
from pathlib import Path
from typing import Any

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry.base import BaseGeometry

from myflopy.advanced import lak_spec
from myflopy.specs import ModelContext, PackageSpec


@dataclass(frozen=True, slots=True)
class LakeConnection:
    """One explicit connection between a lake and a groundwater cell."""

    cellid: tuple[int, int]
    connection_type: str = "VERTICAL"
    bed_leakance: float | None = None
    bottom_elevation: float = 0.0
    top_elevation: float = 0.0
    connection_length: float = 0.0
    connection_width: float = 0.0

    def __post_init__(self) -> None:
        connection_type = self.connection_type.upper()
        if connection_type not in {"VERTICAL", "HORIZONTAL"}:
            raise ValueError("connection_type must be 'VERTICAL' or 'HORIZONTAL'.")
        object.__setattr__(self, "connection_type", connection_type)


@dataclass(frozen=True, slots=True)
class LakeOutlet:
    """A semantic LAK outlet definition keyed by stable lake IDs."""

    source: str
    receiver: str | None = None
    outlet_type: str = "SPECIFIED"
    invert: float = 0.0
    width: float = 0.0
    roughness: float = 0.0
    slope: float = 0.0
    rate: Any = None

    def __post_init__(self) -> None:
        outlet_type = self.outlet_type.upper()
        if outlet_type not in {"SPECIFIED", "MANNING", "WEIR"}:
            raise ValueError("outlet_type must be 'SPECIFIED', 'MANNING', or 'WEIR'.")
        object.__setattr__(self, "outlet_type", outlet_type)


@dataclass(frozen=True, slots=True)
class LakeTable:
    """Prepared stage-volume-area data for one lake."""

    rows: tuple[tuple[float, float, float], ...]
    filename: str | None = None

    def __post_init__(self) -> None:
        rows = tuple(tuple(float(value) for value in row) for row in self.rows)
        if not rows:
            raise ValueError("A lake table must contain at least one row.")
        if any(len(row) != 3 for row in rows):
            raise ValueError("Lake table rows must contain stage, volume, and surface area.")
        stages = [row[0] for row in rows]
        volumes = [row[1] for row in rows]
        areas = [row[2] for row in rows]
        if any(current <= previous for previous, current in zip(stages, stages[1:])):
            raise ValueError("Lake table stages must be strictly increasing.")
        if any(current < previous for previous, current in zip(volumes, volumes[1:])):
            raise ValueError("Lake table volumes must be nondecreasing.")
        if any(value < 0.0 for value in areas):
            raise ValueError("Lake table surface areas cannot be negative.")
        object.__setattr__(self, "rows", rows)

    @property
    def minimum_stage(self) -> float:
        return self.rows[0][0]

    @property
    def maximum_stage(self) -> float:
        return self.rows[-1][0]


@dataclass(frozen=True, slots=True)
class LakeTableBuilder:
    """Build explicit, rectangular, or DEM-derived lake tables."""

    rows: Sequence[Sequence[float]] | None = None
    area: float | None = None
    top: float | None = None
    bottom: float | None = None
    stage_step: float = 1.0
    storage_coefficient: float = 1.0
    dem: Path | None = None
    footprint: Path | BaseGeometry | None = None
    stages: Sequence[float] | None = None
    footprint_buffer: float = 0.0
    filename: str | None = None

    def _footprint_geometry(self) -> BaseGeometry:
        if isinstance(self.footprint, BaseGeometry):
            geometry = self.footprint
        elif isinstance(self.footprint, Path):
            geometries = gpd.read_file(self.footprint).geometry
            geometry = geometries.union_all()
        else:
            raise ValueError("DEM-derived lake tables require a footprint geometry or path.")
        return geometry.buffer(self.footprint_buffer) if self.footprint_buffer else geometry

    def _build_from_dem(self) -> LakeTable:
        if self.stages is None:
            raise ValueError("DEM-derived lake tables require stages.")
        geometry = self._footprint_geometry()
        with rasterio.open(self.dem) as source:
            dem = source.read(1, masked=True)
            inside = geometry_mask(
                [geometry],
                out_shape=dem.shape,
                transform=source.transform,
                invert=True,
            )
            valid = inside & ~np.ma.getmaskarray(dem)
            elevations = np.asarray(dem.filled(np.nan), dtype=float)
            cell_area = abs(float(source.transform.a * source.transform.e))

        rows = []
        for stage in self.stages:
            depth = np.where(valid, np.maximum(float(stage) - elevations, 0.0), 0.0)
            wet = depth > 0.0
            rows.append(
                (
                    float(stage),
                    float(depth.sum() * cell_area * self.storage_coefficient),
                    float(wet.sum() * cell_area),
                )
            )
        return LakeTable(tuple(rows), filename=self.filename)

    def build(self) -> LakeTable:
        """Return a validated immutable lake table."""

        if self.rows is not None:
            return LakeTable(tuple(tuple(row) for row in self.rows), filename=self.filename)
        if self.dem is not None:
            return self._build_from_dem()
        if self.area is None or self.top is None or self.bottom is None:
            raise ValueError(
                "LakeTableBuilder requires rows; DEM, footprint, and stages; "
                "or area, top, and bottom."
            )
        if self.area <= 0.0:
            raise ValueError("area must be positive.")
        if self.top <= self.bottom:
            raise ValueError("top must be above bottom.")
        if self.stage_step <= 0.0:
            raise ValueError("stage_step must be positive.")

        stages = list(np.arange(self.bottom, self.top, self.stage_step))
        if not stages or not np.isclose(stages[-1], self.top):
            stages.append(self.top)
        rows = tuple(
            (
                float(stage),
                float(self.area * (stage - self.bottom) * self.storage_coefficient),
                float(self.area),
            )
            for stage in stages
        )
        return LakeTable(rows, filename=self.filename)


def _build_lak_with_tables(model, *, lake_tables=(), **options):
    """Build LAK and its child table packages."""

    options = {
        key: value.format(model_name=model.name) if isinstance(value, str) else value
        for key, value in options.items()
    }
    package = flopy.mf6.ModflowGwflak(model, **options)
    for table in lake_tables:
        flopy.mf6.ModflowUtllaktab(
            model=model,
            nrow=len(table["rows"]),
            ncol=3,
            table=table["rows"],
            filename=table["filename"],
            pname=table["pname"],
        )
    return package


@dataclass(frozen=True, slots=True)
class LAKBuilder:
    """Engine that turns lake polygon geometry into an MF6 LAK package.

    Finds the lake cells from each polygon footprint, builds the horizontal +
    vertical lake-aquifer bed connections, assigns per-lake inputs (stage, bottom,
    bed leakance, outlets, tables, forcings), and emits a
    :class:`~myflopy.specs.PackageSpec`. Resolves a mover endpoint via
    :meth:`connection` and lake ids via :meth:`lake_number`.

    This is the engine under the package-first ``mf.lak(...)`` facade -- prefer that
    front door for new work; use ``LAKBuilder`` directly only when you need the
    builder object (e.g. to read ``.lake_cells`` mid-build). Construct it with a
    :class:`~myflopy.specs.ModelContext` carrying the grid + surfaces, then call
    :meth:`build`.
    """

    context: ModelContext
    nper: int
    lakes: Path | str | Sequence[Path | str] | gpd.GeoDataFrame
    lake_id_field: str | None = None
    starting_stage: Any = None
    lake_bottom: Any = None
    bed_leakance: Any = 1.0
    connection_modes: str | Mapping[str, str | Sequence[LakeConnection]] = "automatic"
    tables: Mapping[str, LakeTable | LakeTableBuilder] = field(default_factory=dict)
    outlets: tuple[LakeOutlet, ...] = ()
    stage: Any = None
    rainfall: Any = None
    evaporation: Any = None
    runoff: Any = None
    withdrawals: Any = None
    inflow: Any = None
    status: Any = None
    name: str = "lak"
    mover: bool = False
    boundnames: bool = True
    length_conversion: float = 1.0
    time_conversion: float = 1.0
    maximum_iterations: int = 100
    maximum_stage_change: float = 1.0e-5
    options: Mapping[str, Any] = field(default_factory=dict)
    _lakes: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False, compare=False)
    _connections: tuple[tuple[str, LakeConnection], ...] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _bottom_cache: dict[tuple[str, int], float] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if self.context.grid is None:
            raise ValueError("ModelContext.grid is required for LAKBuilder.")
        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        object.__setattr__(self, "tables", dict(self.tables))
        object.__setattr__(self, "outlets", tuple(self.outlets))

    def with_updates(self, **updates: Any) -> LAKBuilder:
        """Return a changed builder without modifying the original."""

        return replace(self, **updates)

    @property
    def grid(self):
        return self.context.grid

    @property
    def lake_table(self) -> gpd.GeoDataFrame:
        """Return normalized lake records with stable lake IDs."""

        if self._lakes is not None:
            return self._lakes
        if isinstance(self.lakes, gpd.GeoDataFrame):
            records = self.lakes.copy()
        else:
            paths = [self.lakes] if isinstance(self.lakes, (str, Path)) else list(self.lakes)
            frames = []
            for path_value in paths:
                path = Path(path_value)
                frame = gpd.read_file(path)
                frame["_source_name"] = path.stem
                frames.append(frame)
            if not frames:
                raise ValueError("At least one lake geometry is required.")
            records = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=frames[0].crs)

        grid_crs = getattr(self.grid, "crs", None)
        if grid_crs is not None and records.crs is not None and records.crs != grid_crs:
            records = records.to_crs(grid_crs)
        if self.lake_id_field is not None:
            if self.lake_id_field not in records:
                raise ValueError(f"Lake ID field is missing: {self.lake_id_field}")
            ids = records[self.lake_id_field].astype(str)
        elif "_source_name" in records and not records["_source_name"].duplicated().any():
            ids = records["_source_name"].astype(str)
        else:
            ids = pd.Series([f"lake_{index}" for index in range(len(records))], index=records.index)
        if ids.duplicated().any():
            duplicates = sorted(ids[ids.duplicated(keep=False)].unique())
            raise ValueError(f"Lake IDs must be unique: {', '.join(duplicates)}")
        if any(geometry.geom_type not in {"Polygon", "MultiPolygon"} for geometry in records.geometry):
            raise ValueError("Lake geometries must be Polygon or MultiPolygon.")

        records = records.copy()
        records["lake_id"] = ids
        records["lake_number"] = range(len(records))
        records = records.set_index("lake_id", drop=False)
        object.__setattr__(self, "_lakes", records)
        return records

    @property
    def lake_ids(self) -> tuple[str, ...]:
        return tuple(self.lake_table.index)

    @property
    def lake_numbers(self) -> dict[str, int]:
        return self.lake_table["lake_number"].astype(int).to_dict()

    def lake_number(self, lake_id: str) -> int:
        try:
            return self.lake_numbers[str(lake_id)]
        except KeyError as error:
            raise ValueError(f"Unknown lake ID: {lake_id}") from error

    def connection(self, lake_id: str) -> "MoverConnection":
        """Return an MVR endpoint for a lake by stable lake ID."""

        from myflopy.modflow.mf6.mvr import MoverConnection

        return MoverConnection(self.name, self.lake_number(lake_id))

    @property
    def lake_cells(self) -> dict[str, list[int]]:
        polygons = self.grid.gdf_vorPolys.geometry
        result = {}
        for lake_id, row in self.lake_table.iterrows():
            intersections = polygons.intersection(row.geometry)
            result[lake_id] = [int(cell) for cell in intersections.index[intersections.area > 0.0]]
        return result

    def _lake_value(self, value: Any, lake_id: str, *, name: str, default: Any = None) -> Any:
        if value is None:
            if default is None:
                raise ValueError(f"{name} is required.")
            return default
        if isinstance(value, str):
            return self.lake_table.loc[lake_id, value] if value in self.lake_table.columns else value
        if isinstance(value, Mapping):
            if lake_id not in value:
                raise ValueError(f"{name} is missing lake '{lake_id}'.")
            return value[lake_id]
        if isinstance(value, Path):
            return value
        if isinstance(value, Real):
            return value
        values = list(value)
        if len(values) != len(self.lake_ids):
            raise ValueError(f"{name} must contain one value per lake.")
        return values[self.lake_number(lake_id)]

    def _connection_mode(self, lake_id: str) -> str | Sequence[LakeConnection]:
        mode = self.connection_modes
        if isinstance(mode, Mapping):
            if lake_id not in mode:
                raise ValueError(f"connection_modes is missing lake '{lake_id}'.")
            mode = mode[lake_id]
        if isinstance(mode, str):
            mode = mode.lower()
            if mode not in {"automatic", "rectangular"}:
                raise ValueError("Generated connection modes must be 'automatic' or 'rectangular'.")
            return mode
        if not all(isinstance(item, LakeConnection) for item in mode):
            raise ValueError("Explicit connection modes must contain LakeConnection objects.")
        return tuple(mode)

    def _surfaces(self) -> pd.DataFrame:
        surfaces = self.context.surfaces
        if surfaces is None:
            surfaces = getattr(self.grid, "gdf_topbtm", None)
        if surfaces is None:
            raise ValueError("ModelContext.surfaces or grid.gdf_topbtm is required.")
        return surfaces.drop(columns="geometry", errors="ignore")

    def _active(self, layer: int, cell: int) -> bool:
        domain = self.context.domain
        if domain is None:
            return True
        values = np.asarray(domain)
        return values.ndim < 2 or bool(values[layer, cell] > 0)

    def _layer_containing(self, cell: int, elevation: float) -> int | None:
        surfaces = self._surfaces().loc[cell].to_numpy(dtype=float)
        for layer, (top, bottom) in enumerate(zip(surfaces, surfaces[1:])):
            if top >= elevation >= bottom and self._active(layer, cell):
                return layer
        return None

    def _bottom(self, lake_id: str, cell: int) -> float:
        key = (lake_id, cell)
        if key in self._bottom_cache:
            return self._bottom_cache[key]
        value = self._lake_value(self.lake_bottom, lake_id, name="lake_bottom")
        if isinstance(value, Path):
            sampled = self.grid.get_raster_vals_at_centroids([value], [lake_id])
            value = sampled.loc[cell, lake_id]
        elif isinstance(value, (Mapping, pd.Series)):
            if cell not in value:
                raise ValueError(f"lake_bottom for '{lake_id}' is missing cell {cell}.")
            value = value[cell]
        resolved = float(value)
        self._bottom_cache[key] = resolved
        return resolved

    def _leakance(self, lake_id: str, connection_type: str, explicit: float | None = None) -> float:
        if explicit is not None:
            return float(explicit)
        if isinstance(self.bed_leakance, Mapping) and (lake_id, connection_type.lower()) in self.bed_leakance:
            value = self.bed_leakance[(lake_id, connection_type.lower())]
        else:
            value = self._lake_value(self.bed_leakance, lake_id, name="bed_leakance")
        if isinstance(value, Mapping):
            value = value.get(connection_type.lower(), value.get(connection_type.upper()))
        if value is None or float(value) <= 0.0:
            raise ValueError(f"bed_leakance must be positive for lake '{lake_id}'.")
        return float(value)

    def _adjacent_with_metrics(self, cell: int) -> list[tuple[int, float, float]]:
        adjacent = list(self.grid.find_adjacent_cells(cell))
        start = int(sum(self.grid.iac[:cell]))
        return [
            (int(other), float(self.grid.cl12[start + index + 1]), float(self.grid.hwva[start + index + 1]))
            for index, other in enumerate(adjacent)
        ]

    def _vertical_connection(self, lake_id: str, cell: int) -> LakeConnection:
        bottom = self._bottom(lake_id, cell)
        layer = self._layer_containing(cell, bottom)
        if layer is None:
            raise ValueError(f"Lake '{lake_id}' bottom {bottom} does not intersect active cell {cell}.")
        return LakeConnection(
            cellid=(layer, cell),
            connection_type="VERTICAL",
            bed_leakance=self._leakance(lake_id, "vertical"),
            bottom_elevation=bottom,
            top_elevation=bottom,
        )

    def _rectangular_vertical_connection(self, lake_id: str, cell: int) -> LakeConnection | None:
        bottom = self._bottom(lake_id, cell)
        top = float(self._lake_value(self.starting_stage, lake_id, name="starting_stage"))
        layer = self._layer_containing(cell, 0.5 * (top + bottom))
        if layer is None:
            return None
        return LakeConnection(
            cellid=(layer, cell),
            connection_type="VERTICAL",
            bed_leakance=self._leakance(lake_id, "vertical"),
            bottom_elevation=bottom,
            top_elevation=top,
        )

    def _sidewall_connections(self, lake_id: str, cell: int, *, rectangular: bool) -> list[LakeConnection]:
        lake_cells = set(self.lake_cells[lake_id])
        bottom = self._bottom(lake_id, cell)
        stage = float(self._lake_value(self.starting_stage, lake_id, name="starting_stage"))
        surfaces = self._surfaces().loc[cell].to_numpy(dtype=float)
        result = []
        for adjacent, length, width in self._adjacent_with_metrics(cell):
            if not rectangular and adjacent in lake_cells:
                continue
            for layer, (cell_top, cell_bottom) in enumerate(zip(surfaces, surfaces[1:])):
                if not self._active(layer, cell):
                    continue
                top = min(cell_top, stage)
                connection_bottom = max(cell_bottom, bottom)
                if top <= connection_bottom:
                    continue
                result.append(
                    LakeConnection(
                        cellid=(layer, cell),
                        connection_type="HORIZONTAL",
                        bed_leakance=self._leakance(lake_id, "horizontal"),
                        bottom_elevation=float(connection_bottom),
                        top_elevation=float(top),
                        connection_length=length,
                        connection_width=width,
                    )
                )
                if not rectangular:
                    break
        return result

    @property
    def connections(self) -> tuple[tuple[str, LakeConnection], ...]:
        """Return generated or explicit connections keyed by stable lake ID."""

        if self._connections is not None:
            return self._connections
        result = []
        for lake_id in self.lake_ids:
            mode = self._connection_mode(lake_id)
            if not isinstance(mode, str):
                result.extend((lake_id, connection) for connection in mode)
                continue
            for cell in self.lake_cells[lake_id]:
                vertical = (
                    self._rectangular_vertical_connection(lake_id, cell)
                    if mode == "rectangular"
                    else self._vertical_connection(lake_id, cell)
                )
                if vertical is not None:
                    result.append((lake_id, vertical))
                result.extend(
                    (lake_id, connection)
                    for connection in self._sidewall_connections(
                        lake_id, cell, rectangular=mode == "rectangular"
                    )
                )
        object.__setattr__(self, "_connections", tuple(result))
        return self._connections

    @property
    def connectiondata(self) -> list[list[Any]]:
        """Return MF6 LAK connection records."""

        counts = {lake_id: 0 for lake_id in self.lake_ids}
        rows = []
        for lake_id, connection in self.connections:
            number = self.lake_number(lake_id)
            index = counts[lake_id]
            counts[lake_id] += 1
            rows.append(
                [
                    number,
                    index,
                    connection.cellid,
                    connection.connection_type,
                    self._leakance(lake_id, connection.connection_type, connection.bed_leakance),
                    connection.bottom_elevation,
                    connection.top_elevation,
                    connection.connection_length,
                    connection.connection_width,
                ]
            )
        return rows

    @property
    def packagedata(self) -> list[list[Any]]:
        """Return MF6 LAK packagedata records."""

        counts = pd.Series([row[0] for row in self.connectiondata]).value_counts().to_dict()
        return [
            [
                self.lake_number(lake_id),
                float(self._lake_value(self.starting_stage, lake_id, name="starting_stage")),
                int(counts.get(self.lake_number(lake_id), 0)),
                *([lake_id] if self.boundnames else []),
            ]
            for lake_id in self.lake_ids
        ]

    def _period_setting(self, value: Any, keyword: str) -> dict[int, list[list[Any]]]:
        result = {period: [] for period in range(self.nper)}
        if value is None:
            return result

        def at_period(current: Any, period: int) -> Any:
            if isinstance(current, Mapping) and current and all(isinstance(key, int) for key in current):
                return current.get(period)
            if isinstance(current, Sequence) and not isinstance(current, (str, bytes)):
                return current[period] if period < len(current) else None
            return current

        for period in range(self.nper):
            if isinstance(value, Mapping) and value and all(isinstance(key, int) for key in value):
                current = value.get(period)
                if current is None:
                    continue
                if isinstance(current, Mapping):
                    for lake_id, amount in current.items():
                        result[period].append([self.lake_number(str(lake_id)), keyword, amount])
                    continue
                for lake_id in self.lake_ids:
                    result[period].append([self.lake_number(lake_id), keyword, current])
                continue

            if isinstance(value, Mapping):
                for lake_id, series in value.items():
                    amount = at_period(series, period)
                    if amount is not None:
                        result[period].append([self.lake_number(str(lake_id)), keyword, amount])
                continue

            for lake_id in self.lake_ids:
                amount = at_period(value, period)
                if amount is not None:
                    result[period].append([self.lake_number(lake_id), keyword, amount])
        return result

    @property
    def outletdata(self) -> list[list[Any]]:
        rows = []
        for outlet_number, outlet in enumerate(self.outlets):
            rows.append(
                [
                    outlet_number,
                    self.lake_number(outlet.source),
                    -1 if outlet.receiver is None else self.lake_number(outlet.receiver),
                    outlet.outlet_type,
                    outlet.invert,
                    outlet.width,
                    outlet.roughness,
                    outlet.slope,
                ]
            )
        return rows

    @property
    def perioddata(self) -> dict[int, list[list[Any]]]:
        """Return combined lake and outlet period settings."""

        result = {period: [] for period in range(self.nper)}
        for value, keyword in (
            (self.stage, "STAGE"),
            (self.rainfall, "RAINFALL"),
            (self.evaporation, "EVAPORATION"),
            (self.runoff, "RUNOFF"),
            (self.withdrawals, "WITHDRAWAL"),
            (self.inflow, "INFLOW"),
            (self.status, "STATUS"),
        ):
            settings = self._period_setting(value, keyword)
            for period, rows in settings.items():
                result[period].extend(rows)
        for outlet_number, outlet in enumerate(self.outlets):
            rate = outlet.rate
            if rate is None:
                continue
            for period in range(self.nper):
                if isinstance(rate, Mapping):
                    amount = rate.get(period)
                elif isinstance(rate, Sequence) and not isinstance(rate, (str, bytes)):
                    amount = rate[period] if period < len(rate) else None
                else:
                    amount = rate
                if amount is not None:
                    result[period].append([outlet_number, "RATE", amount])
        return result

    @property
    def prepared_tables(self) -> dict[str, LakeTable]:
        result = {}
        for lake_id, table in self.tables.items():
            if lake_id not in self.lake_numbers:
                raise ValueError(f"Lake table references unknown lake '{lake_id}'.")
            result[lake_id] = table.build() if isinstance(table, LakeTableBuilder) else table
        return result

    def validate(self) -> None:
        """Validate generated LAK inputs before building the package spec."""

        if any(not cells for cells in self.lake_cells.values()):
            missing = [lake_id for lake_id, cells in self.lake_cells.items() if not cells]
            raise ValueError(f"Lakes do not intersect the grid: {', '.join(missing)}")
        if any(not isinstance(outlet, LakeOutlet) for outlet in self.outlets):
            raise ValueError("outlets must contain LakeOutlet objects.")
        for outlet in self.outlets:
            self.lake_number(outlet.source)
            if outlet.receiver is not None:
                self.lake_number(outlet.receiver)
        for lake_id in self.lake_ids:
            stage = float(self._lake_value(self.starting_stage, lake_id, name="starting_stage"))
            if isinstance(self._connection_mode(lake_id), str):
                bottoms = [self._bottom(lake_id, cell) for cell in self.lake_cells[lake_id]]
                if stage < min(bottoms):
                    raise ValueError(f"starting_stage is below lake_bottom for lake '{lake_id}'.")
            table = self.prepared_tables.get(lake_id)
            if table is not None and not table.minimum_stage <= stage <= table.maximum_stage:
                raise ValueError(f"starting_stage is outside the lake table range for lake '{lake_id}'.")
        self.connectiondata
        self.prepared_tables

    def build(self) -> PackageSpec:
        """Validate stored configuration and return its package spec."""

        self.validate()
        tables = self.prepared_tables
        table_records = []
        child_tables = []
        for lake_id, table in tables.items():
            filename = table.filename or f"{self.name}_{lake_id}.lak.tab"
            table_records.append([self.lake_number(lake_id), filename])
            child_tables.append(
                {
                    "rows": [list(row) for row in table.rows],
                    "filename": filename,
                    "pname": f"{self.name}_{lake_id}_table",
                }
            )
        base = lak_spec(
            self.packagedata,
            self.connectiondata,
            self.perioddata,
            name=self.name,
            noutlets=len(self.outlets),
            ntables=len(tables),
            mover=self.mover,
            boundnames=self.boundnames,
            tables=table_records or None,
            outlets=self.outletdata or None,
            length_conversion=self.length_conversion,
            time_conversion=self.time_conversion,
            maximum_iterations=self.maximum_iterations,
            maximum_stage_change=self.maximum_stage_change,
            **dict(self.options),
        )
        return replace(
            base,
            builder=_build_lak_with_tables,
            options={**base.options, "lake_tables": tuple(child_tables)},
            metadata={
                **base.metadata,
                "builder": "LAKBuilder",
                "lake_ids": list(self.lake_ids),
                "connection_modes": {
                    lake_id: (
                        self._connection_mode(lake_id)
                        if isinstance(self._connection_mode(lake_id), str)
                        else "explicit"
                    )
                    for lake_id in self.lake_ids
                },
            },
        )


__all__ = [
    "LAKBuilder",
    "LakeConnection",
    "LakeOutlet",
    "LakeTable",
    "LakeTableBuilder",
]
