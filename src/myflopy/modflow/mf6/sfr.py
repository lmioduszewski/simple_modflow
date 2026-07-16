"""Build MODFLOW 6 SFR package specifications from stream geometries."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Real
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from myflopy.modflow.mf6.mvr import MoverConnection

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from myflopy.advanced import sfr_spec
from myflopy.specs import (
    ModelContext,
    PackageSpec,
    mf6_length_conversion,
    mf6_time_conversion,
)


StreamLocation = str | Point


@dataclass(frozen=True, slots=True)
class StreamConnection:
    """An explicit reach-to-reach link between two named streams.

    Pass these to ``mf.sfr(connections=(...))`` (or ``SFRBuilder``) to disambiguate
    confluences that automatic geometry-based topology cannot resolve -- e.g. two
    tributaries meeting a main stem. Water leaves ``source`` (by default at its
    downstream end) and enters ``receiver`` (by default at its nearest reach).

    Attributes
    ----------
    source, receiver
        Stream ids (the ``stream_id`` attribute values). ``receiver=None`` marks
        ``source`` as a terminal outflow.
    source_location, receiver_location
        Where on each stream to connect: ``"downstream"``/``"upstream"``,
        ``"nearest"``, ``"intersection"``, or a shapely ``Point``.

    Examples
    --------
    >>> mf.sfr(..., connections=(mf.StreamConnection("north_trib", "main_stem"),
    ...                          mf.StreamConnection("south_trib", "main_stem")))
    """

    source: str
    receiver: str | None
    source_location: StreamLocation = "downstream"
    receiver_location: StreamLocation = "nearest"


@dataclass(frozen=True, slots=True)
class StreamDiversion:
    """An intentional diversion that splits flow from one stream into another.

    Like :class:`StreamConnection` but moves only a portion of the flow (a
    fraction, an excess, or an up-to amount) -- a canal headgate, a bypass.

    Attributes
    ----------
    source, receiver
        Stream ids for the donor and the diversion channel.
    amount
        How much to divert, interpreted per ``priority``.
    source_location, receiver_location
        Connection points (see :class:`StreamConnection`).
    priority
        Diversion rule: ``"FRACTION"`` (default), ``"EXCESS"``, ``"THRESHOLD"`` or
        ``"UPTO"`` -- the MF6 SFR CPRIOR options.
    """

    source: str
    receiver: str
    amount: Any
    source_location: StreamLocation = "downstream"
    receiver_location: StreamLocation = "nearest"
    priority: str = "FRACTION"


@dataclass(frozen=True, slots=True)
class StreamNetwork:
    """The resolved stream-to-stream topology of an SFR network (grid-independent).

    A snapshot of how the named streams connect, produced by
    :class:`SFRBuilder` while assembling reaches. It captures the logical wiring --
    which streams flow into which (``connections``), which split off as diversions
    (``diversions``), and which streams could not be linked automatically
    (``unresolved``) -- separately from the per-reach, grid-dependent records. Use
    it to inspect or validate connectivity (e.g. confirm tributaries reach the
    main stem) before the SFR package is built.

    Attributes
    ----------
    stream_ids
        Names of all streams in the network.
    connections
        Resolved confluences as :class:`StreamConnection` records.
    diversions
        Resolved splits as :class:`StreamDiversion` records.
    unresolved
        Stream ids whose downstream link could not be inferred and need an
        explicit connection.
    """

    stream_ids: tuple[str, ...]
    connections: tuple[StreamConnection, ...]
    diversions: tuple[StreamDiversion, ...]
    unresolved: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class SFRBuilder:
    """Engine that turns stream centerline geometry into an MF6 SFR package.

    Discretizes each stream line into reaches on the grid, resolves the reach-to-
    reach topology (geometry-based, with explicit :class:`StreamConnection` /
    :class:`StreamDiversion` overrides), assigns per-reach properties, and emits a
    :class:`~myflopy.specs.PackageSpec`. Also resolves semantic mover endpoints via
    :meth:`connection` / :meth:`outlet_reach`.

    This is the engine under the package-first ``mf.sfr(...)`` facade -- prefer that
    front door for new work; use ``SFRBuilder`` directly only when you need the
    builder object (e.g. to read ``.reaches`` / ``.stream_cells`` mid-build).
    Construct it with a :class:`~myflopy.specs.ModelContext` carrying the grid, then
    call :meth:`build`.
    """

    context: ModelContext
    nper: int
    streams: Path | str | Sequence[Path | str] | gpd.GeoDataFrame
    stream_id: str | None = None
    from_node: str | None = None
    to_node: str | None = None
    connection_mode: str = "automatic"
    connections: tuple[StreamConnection, ...] = ()
    diversions: tuple[StreamDiversion, ...] = ()
    connection_tolerance: float | None = None
    reverse_streams: tuple[str, ...] = ()
    reach_layer: int | str | Mapping[str, int] = "top_active"
    width: Any = 10.0
    gradient: Any = 0.001
    reach_top: Any = None
    roughness: Any = 0.03
    streambed_k: Any = 1.0
    streambed_thickness: Any = 1.0
    inflow: Any = None
    rainfall: Any = None
    evaporation: Any = None
    runoff: Any = None
    status: Any = None
    name: str = "sfr"
    mover: bool = False
    length_conversion: float | None = None  # None -> derive from context.length_units
    time_conversion: float | None = None     # None -> derive from context.time_units
    maximum_picard_iterations: int = 1
    maximum_iterations: int = 1000
    maximum_depth_change: float = 0.01
    options: Mapping[str, Any] = field(default_factory=dict)
    _streams: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False, compare=False)
    _network: StreamNetwork | None = field(default=None, init=False, repr=False, compare=False)
    _reaches: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate the grid, ``nper``, and ``connection_mode``, and freeze the sequence fields."""

        if self.context.grid is None:
            raise ValueError("ModelContext.grid is required for SFRBuilder.")
        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        if self.connection_mode not in {"automatic", "nodes", "explicit"}:
            raise ValueError("connection_mode must be 'automatic', 'nodes', or 'explicit'.")
        object.__setattr__(self, "connections", tuple(self.connections))
        object.__setattr__(self, "diversions", tuple(self.diversions))
        object.__setattr__(self, "reverse_streams", tuple(str(value) for value in self.reverse_streams))

    def with_updates(self, **updates: Any) -> SFRBuilder:
        """Return a changed builder without modifying the original."""

        return replace(self, **updates)

    @property
    def grid(self):
        """The Voronoi grid carried by the model context."""

        return self.context.grid

    @property
    def stream_table(self) -> gpd.GeoDataFrame:
        """Return normalized stream records with stable stream IDs."""

        if self._streams is not None:
            return self._streams

        if isinstance(self.streams, gpd.GeoDataFrame):
            records = self.streams.copy()
        else:
            paths = [self.streams] if isinstance(self.streams, (str, Path)) else list(self.streams)
            frames = []
            for path_value in paths:
                path = Path(path_value)
                frame = gpd.read_file(path)
                frame["_source_name"] = path.stem
                frames.append(frame)
            records = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=frames[0].crs)

        grid_crs = getattr(self.grid, "crs", None)
        if grid_crs is not None and records.crs is not None and records.crs != grid_crs:
            records = records.to_crs(grid_crs)

        if self.stream_id is not None:
            if self.stream_id not in records:
                raise ValueError(f"Stream ID field is missing: {self.stream_id}")
            ids = records[self.stream_id].astype(str)
        elif "_source_name" in records and not records["_source_name"].duplicated().any():
            ids = records["_source_name"].astype(str)
        else:
            ids = pd.Series([f"stream_{index}" for index in range(len(records))], index=records.index)

        if ids.duplicated().any():
            duplicates = sorted(ids[ids.duplicated(keep=False)].unique())
            raise ValueError(f"Stream IDs must be unique: {', '.join(duplicates)}")

        records = records.copy()
        records["stream_id"] = ids
        for index, row in records.iterrows():
            if row.geometry.geom_type != "LineString":
                raise ValueError(
                    f"Stream '{row['stream_id']}' must have LineString geometry, "
                    f"got {row.geometry.geom_type}."
                )
            if row["stream_id"] in self.reverse_streams:
                records.at[index, "geometry"] = row.geometry.reverse()
        records = records.set_index("stream_id", drop=False)
        object.__setattr__(self, "_streams", records)
        return records

    @property
    def tolerance(self) -> float:
        """Return the geometry tolerance used for automatic connections."""

        if self.connection_tolerance is not None:
            return float(self.connection_tolerance)
        areas = np.asarray(self.grid.gdf_vorPolys.geometry.area, dtype=float)
        return max(float(np.sqrt(np.median(areas)) * 0.05), 1.0e-9)

    def _automatic_connections(self) -> tuple[list[StreamConnection], list[str]]:
        """Discover downstream connections by snapping each stream's end to the nearest stream.

        Returns ``(connections, unresolved)``; streams with no neighbor within
        ``tolerance`` are unresolved, and an ambiguous nearest tie raises.
        """

        explicit_sources = {connection.source for connection in self.connections}
        connections: list[StreamConnection] = []
        unresolved: list[str] = []
        records = self.stream_table

        for stream_id, row in records.iterrows():
            if stream_id in explicit_sources:
                continue
            endpoint = Point(row.geometry.coords[-1])
            candidates = [
                (str(other_id), endpoint.distance(other.geometry))
                for other_id, other in records.iterrows()
                if other_id != stream_id
            ]
            matches = sorted((candidate for candidate in candidates if candidate[1] <= self.tolerance), key=lambda x: x[1])
            if not matches:
                unresolved.append(str(stream_id))
                continue
            if len(matches) > 1 and np.isclose(matches[0][1], matches[1][1]):
                names = ", ".join(match[0] for match in matches if np.isclose(match[1], matches[0][1]))
                raise ValueError(f"Automatic connection for '{stream_id}' is ambiguous: {names}")
            connections.append(StreamConnection(str(stream_id), matches[0][0]))
        return connections, unresolved

    def _node_connections(self) -> tuple[list[StreamConnection], list[str]]:
        """Build connections from the ``from_node``/``to_node`` fields (topology table).

        Returns ``(connections, unresolved)``; a stream whose to-node matches
        several from-nodes raises as ambiguous.
        """

        if self.from_node is None or self.to_node is None:
            raise ValueError("from_node and to_node fields are required for connection_mode='nodes'.")
        for field_name in (self.from_node, self.to_node):
            if field_name not in self.stream_table:
                raise ValueError(f"Stream node field is missing: {field_name}")

        by_from_node: dict[Any, list[str]] = {}
        for stream_id, row in self.stream_table.iterrows():
            by_from_node.setdefault(row[self.from_node], []).append(str(stream_id))

        connections: list[StreamConnection] = []
        unresolved: list[str] = []
        for stream_id, row in self.stream_table.iterrows():
            receivers = by_from_node.get(row[self.to_node], [])
            if not receivers:
                unresolved.append(str(stream_id))
            elif len(receivers) > 1:
                raise ValueError(f"Node connection for '{stream_id}' has multiple receivers: {receivers}")
            else:
                connections.append(StreamConnection(str(stream_id), receivers[0]))
        return connections, unresolved

    @property
    def network(self) -> StreamNetwork:
        """Return resolved semantic stream connections before grid mapping."""

        if self._network is not None:
            return self._network

        explicit = {connection.source: connection for connection in self.connections}
        if self.connection_mode == "automatic":
            discovered, unresolved = self._automatic_connections()
        elif self.connection_mode == "nodes":
            discovered, unresolved = self._node_connections()
        else:
            discovered, unresolved = [], []

        resolved = {connection.source: connection for connection in discovered}
        resolved.update(explicit)
        stream_ids = tuple(str(value) for value in self.stream_table.index)
        unknown = sorted(
            {
                value
                for connection in resolved.values()
                for value in (connection.source, connection.receiver)
                if value is not None and value not in stream_ids
            }
        )
        if unknown:
            raise ValueError(f"Stream connections reference unknown streams: {', '.join(unknown)}")
        unknown_diversions = sorted(
            {
                value
                for diversion in self.diversions
                for value in (diversion.source, diversion.receiver)
                if value not in stream_ids
            }
        )
        if unknown_diversions:
            raise ValueError(f"Stream diversions reference unknown streams: {', '.join(unknown_diversions)}")
        invalid_priorities = sorted(
            {
                diversion.priority
                for diversion in self.diversions
                if diversion.priority.upper() not in {"FRACTION", "EXCESS", "THRESHOLD", "UPTO"}
            }
        )
        if invalid_priorities:
            raise ValueError(f"Unsupported diversion priorities: {', '.join(invalid_priorities)}")

        network = StreamNetwork(
            stream_ids=stream_ids,
            connections=tuple(resolved.values()),
            diversions=self.diversions,
            unresolved=tuple(stream for stream in unresolved if stream not in resolved),
        )
        self._validate_no_cycles(network.connections)
        object.__setattr__(self, "_network", network)
        return network

    @staticmethod
    def _validate_no_cycles(connections: Sequence[StreamConnection]) -> None:
        """Raise if the downstream routing graph contains a cycle."""

        downstream = {item.source: item.receiver for item in connections if item.receiver is not None}
        for start in downstream:
            seen = set()
            current = start
            while current in downstream:
                if current in seen:
                    raise ValueError(f"Stream connections contain a routing cycle involving '{current}'.")
                seen.add(current)
                current = downstream[current]

    def _ordered_stream_ids(self) -> list[str]:
        """Order source streams before their receivers for valid signed routing."""

        ids = list(self.network.stream_ids)
        edges = [(item.source, item.receiver) for item in self.network.connections if item.receiver is not None]
        ordered: list[str] = []
        remaining = set(ids)
        while remaining:
            ready = [item for item in ids if item in remaining and not any(receiver == item and source in remaining for source, receiver in edges)]
            if not ready:
                raise ValueError("Stream network could not be ordered.")
            ordered.extend(ready)
            remaining.difference_update(ready)
        return ordered

    def _layer_for_cell(self, stream_id: str, cell: int) -> int:
        """The model layer a reach sits in: a constant, a per-stream map, a field, or ``top_active``."""

        if isinstance(self.reach_layer, int):
            return self.reach_layer
        if isinstance(self.reach_layer, Mapping):
            return int(self.reach_layer[stream_id])
        if self.reach_layer != "top_active":
            row = self.stream_table.loc[stream_id]
            if self.reach_layer not in row:
                raise ValueError(f"Reach-layer field is missing: {self.reach_layer}")
            return int(row[self.reach_layer])

        domain = self.context.domain
        if domain is None:
            return 0
        values = np.asarray(domain)
        if values.ndim == 1:
            return 0
        active = np.flatnonzero(values[:, cell] > 0)
        if len(active) == 0:
            raise ValueError(f"Stream '{stream_id}' intersects inactive grid column {cell}.")
        return int(active[0])

    @property
    def reaches(self) -> gpd.GeoDataFrame:
        """Return the inspectable grid-dependent SFR reach table."""

        if self._reaches is not None:
            return self._reaches

        rows = []
        rno = 0
        for stream_id in self._ordered_stream_ids():
            line = self.stream_table.loc[stream_id].geometry
            intersecting = self.grid.gdf_vorPolys[self.grid.gdf_vorPolys.intersects(line)]
            segments = []
            for cell, grid_row in intersecting.iterrows():
                segment = grid_row.geometry.intersection(line)
                if segment.length <= 0:
                    continue
                midpoint = segment.interpolate(0.5, normalized=True)
                segments.append((line.project(midpoint), int(cell), segment))
            segments.sort(key=lambda item: item[0])
            if not segments:
                raise ValueError(f"Stream '{stream_id}' does not intersect any grid cells.")
            for sequence, (_, cell, segment) in enumerate(segments):
                rows.append(
                    {
                        "rno": rno,
                        "stream_id": stream_id,
                        "sequence": sequence,
                        "cellid": (self._layer_for_cell(stream_id, cell), cell),
                        "rlen": float(segment.length),
                        "geometry": segment,
                    }
                )
                rno += 1
        reaches = gpd.GeoDataFrame(rows, geometry="geometry", crs=self.stream_table.crs).set_index("rno", drop=False)
        object.__setattr__(self, "_reaches", reaches)
        return reaches

    @property
    def total_nreaches(self) -> int:
        """The total number of generated reaches across all streams."""

        return len(self.reaches)

    @property
    def stream_ids(self) -> tuple[str, ...]:
        """The distinct stream IDs, in reach order."""

        return tuple(self.reaches["stream_id"].drop_duplicates())

    @property
    def stream_reaches(self) -> dict[str, list[int]]:
        """Map each stream ID to its ordered list of reach numbers."""

        return {
            stream_id: group["rno"].astype(int).tolist()
            for stream_id, group in self.reaches.groupby("stream_id", sort=False)
        }

    @property
    def stream_cells(self) -> dict[str, list[int]]:
        """Map each stream ID to the grid cells its reaches occupy, in order."""

        return {
            stream_id: [int(cellid[1]) for cellid in group["cellid"]]
            for stream_id, group in self.reaches.groupby("stream_id", sort=False)
        }

    def outlet_reach(self, stream_id: str) -> int:
        """Return the final generated reach for one stream."""

        return self.stream_reaches[str(stream_id)][-1]

    def connection(self, stream_id: str, location: StreamLocation = "downstream") -> "MoverConnection":
        """Return an MVR endpoint for a stream by stable stream ID and location."""

        from myflopy.modflow.mf6.mvr import MoverConnection

        return MoverConnection(self.name, self._reach_at(str(stream_id), location))

    def _point_for_location(self, stream_id: str, location: StreamLocation, *, other: str | None = None) -> Point:
        """Resolve a stream location keyword/Point to a coordinate on the stream line.

        Supports an explicit ``Point``, ``upstream``/``downstream`` ends, and
        ``nearest``/``intersection`` against ``other`` (required for those).
        """

        line = self.stream_table.loc[stream_id].geometry
        if isinstance(location, Point):
            return location
        if location in {"downstream", "end"}:
            return Point(line.coords[-1])
        if location in {"upstream", "start"}:
            return Point(line.coords[0])
        if location in {"nearest", "automatic", "intersection"}:
            if other is None:
                raise ValueError(f"Location '{location}' requires another stream.")
            other_line = self.stream_table.loc[other].geometry
            if location == "intersection":
                intersection = line.intersection(other_line)
                if intersection.is_empty:
                    raise ValueError(f"Streams '{stream_id}' and '{other}' do not intersect.")
                return intersection.representative_point()
            return line.interpolate(line.project(other_line.representative_point()))
        raise ValueError(f"Unsupported stream location: {location!r}")

    def _reach_at(self, stream_id: str, location: StreamLocation, *, other: str | None = None) -> int:
        """The reach number on ``stream_id`` closest to the resolved ``location`` point."""

        point = self._point_for_location(stream_id, location, other=other)
        reaches = self.reaches[self.reaches["stream_id"] == stream_id]
        distances = reaches.geometry.distance(point)
        return int(distances.idxmin())

    @property
    def resolved_connections(self) -> tuple[tuple[int, int], ...]:
        """Return source/receiver reach pairs for ordinary stream connections."""

        pairs = []
        for connection in self.network.connections:
            if connection.receiver is None:
                continue
            source = self._reach_at(
                connection.source,
                connection.source_location,
                other=connection.receiver,
            )
            receiver = self._reach_at(
                connection.receiver,
                connection.receiver_location,
                other=connection.source,
            )
            pairs.append((source, receiver))
        return tuple(pairs)

    @property
    def resolved_diversions(self) -> tuple[tuple[StreamDiversion, int, int, int], ...]:
        """Return diversion definitions with generated reach and diversion IDs."""

        counts: dict[int, int] = {}
        resolved = []
        for diversion in self.diversions:
            source = self._reach_at(diversion.source, diversion.source_location, other=diversion.receiver)
            receiver = self._reach_at(diversion.receiver, diversion.receiver_location, other=diversion.source)
            idiv = counts.get(source, 0)
            counts[source] = idiv + 1
            resolved.append((diversion, source, receiver, idiv))
        return tuple(resolved)

    @property
    def connectiondata(self) -> list[list[int]]:
        """Return signed MF6 SFR reach connections."""

        upstream: dict[int, list[int]] = {rno: [] for rno in self.reaches.index}
        downstream: dict[int, list[int]] = {rno: [] for rno in self.reaches.index}
        for reaches in self.stream_reaches.values():
            for source, receiver in zip(reaches, reaches[1:]):
                downstream[source].append(receiver)
                upstream[receiver].append(source)
        for source, receiver in self.resolved_connections:
            downstream[source].append(receiver)
            upstream[receiver].append(source)

        return [
            [int(rno), *sorted(set(upstream[rno])), *[-value for value in sorted(set(downstream[rno]))]]
            for rno in self.reaches.index
        ]

    def _reach_values(self, value: Any, *, name: str, default: Any = None) -> list[Any]:
        """Normalize a static property to generated reaches."""

        if value is None:
            if default is None:
                raise ValueError(f"{name} is required.")
            value = default
        if isinstance(value, Real) or isinstance(value, str) and value not in self.stream_table:
            if isinstance(value, str) and value in self.stream_table.columns:
                pass
            else:
                return [value] * self.total_nreaches
        if isinstance(value, str):
            if value not in self.stream_table.columns:
                raise ValueError(f"Stream field is missing for {name}: {value}")
            by_stream = self.stream_table[value].to_dict()
            return [by_stream[stream_id] for stream_id in self.reaches["stream_id"]]
        if isinstance(value, Mapping):
            return [
                value.get(int(rno), value.get(stream_id, value.get(cellid)))
                for rno, stream_id, cellid in self.reaches[["rno", "stream_id", "cellid"]].itertuples(index=False, name=None)
            ]
        values = list(value)
        if len(values) == self.total_nreaches:
            return values
        if len(values) == len(self.stream_ids):
            by_stream = dict(zip(self.stream_ids, values))
            return [by_stream[stream_id] for stream_id in self.reaches["stream_id"]]
        raise ValueError(f"{name} must contain one value per reach or stream.")

    def _default_reach_tops(self) -> list[float]:
        """Per-reach streambed-top elevations from the grid surface, enforced monotone downstream.

        Samples each reach cell's layer top, then clamps so a reach is never higher
        than the reach above it in the same stream.
        """

        if getattr(self.grid, "gdf_topbtm", None) is None:
            raise ValueError("reach_top is required when grid surface elevations are unavailable.")
        values = []
        for layer, cell in self.reaches["cellid"]:
            top_column = 0 if layer == 0 else layer
            values.append(float(self.grid.gdf_topbtm.loc[cell, top_column]))
        for reaches in self.stream_reaches.values():
            previous = None
            for rno in reaches:
                if previous is not None and values[rno] > previous:
                    values[rno] = previous
                previous = values[rno]
        return values

    @property
    def packagedata(self) -> list[list[Any]]:
        """Return MF6 SFR packagedata records."""

        widths = self._reach_values(self.width, name="width")
        gradients = self._reach_values(self.gradient, name="gradient")
        tops = self._default_reach_tops() if self.reach_top is None else self._reach_values(self.reach_top, name="reach_top")
        thickness = self._reach_values(self.streambed_thickness, name="streambed_thickness")
        conductivity = self._reach_values(self.streambed_k, name="streambed_k")
        roughness = self._reach_values(self.roughness, name="roughness")
        ncon = [len(row) - 1 for row in self.connectiondata]
        diversion_counts: dict[int, int] = {}
        receiving_diversions = set()
        for _, source, receiver, _ in self.resolved_diversions:
            diversion_counts[source] = diversion_counts.get(source, 0) + 1
            receiving_diversions.add(receiver)

        return [
            [
                int(rno),
                self.reaches.loc[rno, "cellid"],
                float(self.reaches.loc[rno, "rlen"]),
                widths[rno],
                gradients[rno],
                tops[rno],
                thickness[rno],
                conductivity[rno],
                roughness[rno],
                ncon[rno],
                0.0 if rno in receiving_diversions else 1.0,
                diversion_counts.get(rno, 0),
            ]
            for rno in self.reaches.index
        ]

    def _period_setting(self, value: Any, keyword: str) -> dict[int, list[list[Any]]]:
        """Expand an inflow/rainfall/etc setting into per-period MF6 ``[reach, keyword, amount]`` rows.

        Accepts a per-period mapping, a per-reach/location mapping, a scalar
        (applied to the first reach), or explicit row tuples.
        """

        if value is None:
            return {period: [] for period in range(self.nper)}
        period_values = value if isinstance(value, Mapping) and all(isinstance(key, int) for key in value) else {period: value for period in range(self.nper)}
        result = {period: [] for period in range(self.nper)}
        for period in range(self.nper):
            current = period_values.get(period)
            if current is None:
                continue
            if isinstance(current, Mapping):
                for location, amount in current.items():
                    rno = int(location) if isinstance(location, int) else self._reach_at(str(location), "upstream")
                    result[period].append([rno, keyword, amount])
            elif isinstance(current, Real) or isinstance(current, str):
                result[period].append([int(self.reaches.index[0]), keyword, current])
            else:
                for row in current:
                    result[period].append([int(row[0]), keyword, *row[1:]])
        return result

    @property
    def perioddata(self) -> dict[int, list[list[Any]]]:
        """Return normalized MF6 SFR perioddata."""

        result = {period: [] for period in range(self.nper)}
        for keyword, value in (
            ("INFLOW", self.inflow),
            ("RAINFALL", self.rainfall),
            ("EVAPORATION", self.evaporation),
            ("RUNOFF", self.runoff),
            ("STATUS", self.status),
        ):
            settings = self._period_setting(value, keyword)
            for period in result:
                result[period].extend(settings[period])

        for diversion, source, _, idiv in self.resolved_diversions:
            amounts = diversion.amount if isinstance(diversion.amount, Mapping) else {period: diversion.amount for period in range(self.nper)}
            for period in result:
                if period in amounts:
                    result[period].append([source, "DIVERSION", idiv, amounts[period]])
        return result

    @property
    def diversiondata(self) -> list[list[Any]] | None:
        """MF6 SFR ``diversions`` records (source reach, index, receiver, priority), or ``None``."""

        rows = [
            [source, idiv, receiver, diversion.priority.upper()]
            for diversion, source, receiver, idiv in self.resolved_diversions
        ]
        return rows or None

    def validate(self) -> None:
        """Raise clear errors for invalid reach, routing, and package data."""

        if self.total_nreaches == 0:
            raise ValueError("SFRBuilder generated no reaches.")
        for row in self.packagedata:
            if row[2] <= 0:
                raise ValueError(f"SFR reach {row[0]} has nonpositive length.")
            if row[3] <= 0:
                raise ValueError(f"SFR reach {row[0]} has nonpositive width.")
            if row[4] < 0:
                raise ValueError(f"SFR reach {row[0]} has negative gradient.")
            if row[6] <= 0:
                raise ValueError(f"SFR reach {row[0]} has nonpositive streambed thickness.")
            if row[7] < 0:
                raise ValueError(f"SFR reach {row[0]} has negative streambed conductivity.")
            if row[8] <= 0:
                raise ValueError(f"SFR reach {row[0]} has nonpositive roughness.")
        valid_reaches = set(self.reaches.index)
        for period, settings in self.perioddata.items():
            for setting in settings:
                if setting[0] not in valid_reaches:
                    raise ValueError(f"SFR period {period} references unknown reach {setting[0]}.")

    def resolved_length_conversion(self) -> float:
        """SFR ``LENGTH_CONVERSION`` -- explicit value, else from context units.

        Manning's streamflow equation is defined in SI; this factor converts it to the
        model's length units. Left at ``None`` (the default) it is taken from
        ``context.length_units`` (default ``"feet"`` -> 3.28081), so a feet/days model
        gets the right constant without the caller remembering to set it.
        """
        if self.length_conversion is not None:
            return self.length_conversion
        return mf6_length_conversion(getattr(self.context, "length_units", "feet"))

    def resolved_time_conversion(self) -> float:
        """SFR ``TIME_CONVERSION`` -- explicit value, else from context units."""
        if self.time_conversion is not None:
            return self.time_conversion
        return mf6_time_conversion(getattr(self.context, "time_units", "days"))

    def build(self) -> PackageSpec:
        """Validate stored configuration and return its package spec."""

        self.validate()
        package = sfr_spec(
            self.packagedata,
            self.connectiondata,
            self.perioddata,
            name=self.name,
            mover=self.mover,
            diversions=self.diversiondata,
            length_conversion=self.resolved_length_conversion(),
            time_conversion=self.resolved_time_conversion(),
            maximum_picard_iterations=self.maximum_picard_iterations,
            maximum_iterations=self.maximum_iterations,
            maximum_depth_change=self.maximum_depth_change,
            **dict(self.options),
        )
        return package.with_metadata(
            builder="SFRBuilder",
            stream_ids=list(self.stream_ids),
            connections=[
                {"source": item.source, "receiver": item.receiver}
                for item in self.network.connections
            ],
            unresolved=list(self.network.unresolved),
        )


__all__ = ["SFRBuilder", "StreamConnection", "StreamDiversion", "StreamNetwork"]
