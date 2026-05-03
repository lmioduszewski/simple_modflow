"""Named model-region and region-group helpers for simulation analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import geopandas as gpd
import pandas as pd

from simple_modflow.modflow.mf6.grid.selection import get_vor_cells_as_series


@dataclass(slots=True)
class ModelRegion:
    """Represents one named set of model cells plus optional metadata."""

    name: str
    cells: list[int]
    layer: int | list[int] | None = None
    geometry: Any = None
    category: str = "selection"
    package: str | None = None
    tags: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    cellids: list[Any] = field(default_factory=list)


@dataclass(slots=True)
class RegionGroup:
    """Represents a named group of regions and/or other groups."""

    name: str
    members: list[str] = field(default_factory=list)
    tags: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    category: str = "group"


def normalize_region_inputs(
    cellids: Iterable[Any] | None,
    *,
    layer: int | list[int] | None = None,
) -> tuple[list[Any], list[int], int | list[int] | None]:
    """Normalize region cell inputs into raw IDs, cell numbers, and layers."""
    if cellids is None:
        raw_cellids: list[Any] = []
    else:
        raw_cellids = list(cellids)

    cells: list[int] = []
    derived_layers: list[int] = []

    for cellid in raw_cellids:
        if isinstance(cellid, tuple):
            if len(cellid) >= 2:
                derived_layers.append(int(cellid[0]))
                cells.append(int(cellid[1]))
            elif len(cellid) == 1:
                cells.append(int(cellid[0]))
        else:
            cells.append(int(cellid))

    unique_cells = sorted(set(cells))

    if layer is None and derived_layers:
        unique_layers = sorted(set(derived_layers))
        if len(unique_layers) == 1:
            layer = unique_layers[0]
        else:
            layer = unique_layers

    return raw_cellids, unique_cells, layer


class RegionRegistry:
    """Registry of named regions and groups attached to one model."""

    def __init__(self, model):
        self.model = model
        self._regions: dict[str, ModelRegion] = {}
        self._groups: dict[str, RegionGroup] = {}

    def __contains__(self, name: str) -> bool:
        return name in self._regions or name in self._groups

    def __getitem__(self, name: str) -> ModelRegion | RegionGroup:
        return self.get(name)

    def add(self, region: ModelRegion, *, overwrite: bool = False) -> ModelRegion:
        """Register a concrete region object."""
        if region.name in self._groups:
            raise ValueError(f"Cannot add region {region.name!r}; a group with that name already exists")
        if not overwrite and region.name in self._regions:
            raise ValueError(f"Region {region.name!r} already exists")
        self._regions[region.name] = region
        return region

    def add_group(
        self,
        name: str,
        *,
        members: Iterable[str] | None = None,
        tags: list[str] | None = None,
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
    ) -> RegionGroup:
        """Create or replace a named group of regions/groups."""
        if name in self._regions:
            raise ValueError(f"Cannot add group {name!r}; a region with that name already exists")
        if not overwrite and name in self._groups:
            raise ValueError(f"Group {name!r} already exists")

        group = RegionGroup(
            name=name,
            members=[] if members is None else [],
            tags=[] if tags is None else list(tags),
            metadata={} if metadata is None else dict(metadata),
        )
        self._groups[name] = group

        if members is not None:
            for member_name in members:
                self.add_to_group(name, member_name)

        return group

    def add_to_group(self, group_name: str, member_name: str) -> RegionGroup:
        """Add an existing region or group to another group."""
        group = self.get_group(group_name)
        if member_name not in self:
            raise KeyError(f"Cannot add {member_name!r} to {group_name!r}; member does not exist")
        if member_name == group_name:
            raise ValueError(f"Group {group_name!r} cannot contain itself")
        self._assert_no_cycle(group_name, member_name)
        if member_name not in group.members:
            group.members.append(member_name)
        return group

    def remove_from_group(self, group_name: str, member_name: str) -> RegionGroup:
        """Remove one member from a named group."""
        group = self.get_group(group_name)
        group.members = [member for member in group.members if member != member_name]
        return group

    def remove(self, name: str) -> ModelRegion | RegionGroup:
        """Remove a region or group and detach it from all groups."""
        if name in self._regions:
            removed = self._regions.pop(name)
        elif name in self._groups:
            removed = self._groups.pop(name)
        else:
            raise KeyError(f"{name!r} does not exist")

        for group in self._groups.values():
            group.members = [member for member in group.members if member != name]

        return removed

    def get(self, name: str) -> ModelRegion | RegionGroup:
        """Return either a region or a group by name."""
        if name in self._regions:
            return self._regions[name]
        if name in self._groups:
            return self._groups[name]
        raise KeyError(f"{name!r} does not exist")

    def get_region(self, name: str) -> ModelRegion:
        """Return one concrete region by name."""
        if name not in self._regions:
            raise KeyError(f"Region {name!r} does not exist")
        return self._regions[name]

    def get_group(self, name: str) -> RegionGroup:
        """Return one region group by name."""
        if name not in self._groups:
            raise KeyError(f"Group {name!r} does not exist")
        return self._groups[name]

    def list(self) -> list[str]:
        """List all region and group names."""
        return sorted([*self._regions.keys(), *self._groups.keys()])

    def list_regions(self) -> list[str]:
        """List only concrete region names."""
        return sorted(self._regions)

    def list_groups(self) -> list[str]:
        """List only region-group names."""
        return sorted(self._groups)

    def region_summary(self) -> pd.DataFrame:
        """Return a tabular summary of concrete regions."""
        rows = []
        for region in self._regions.values():
            rows.append(
                {
                    "name": region.name,
                    "kind": "region",
                    "category": region.category,
                    "package": region.package,
                    "layer": region.layer,
                    "num_cells": len(region.cells),
                    "tags": list(region.tags),
                }
            )
        return pd.DataFrame(rows).sort_values("name") if rows else pd.DataFrame(
            columns=["name", "kind", "category", "package", "layer", "num_cells", "tags"]
        )

    def group_summary(self) -> pd.DataFrame:
        """Return a tabular summary of region groups."""
        rows = []
        for group in self._groups.values():
            rows.append(
                {
                    "name": group.name,
                    "kind": "group",
                    "category": group.category,
                    "num_members": len(group.members),
                    "members": list(group.members),
                    "tags": list(group.tags),
                }
            )
        return pd.DataFrame(rows).sort_values("name") if rows else pd.DataFrame(
            columns=["name", "kind", "category", "num_members", "members", "tags"]
        )

    def summary(self) -> pd.DataFrame:
        """Return a combined summary of regions and groups."""
        region_summary = self.region_summary()
        group_summary = self.group_summary()
        if region_summary.empty and group_summary.empty:
            return pd.DataFrame(columns=["name", "kind", "category"])
        return pd.concat([region_summary, group_summary], ignore_index=True, sort=False).sort_values("name")

    def add_from_cells(
        self,
        name: str,
        cellids: Iterable[Any],
        *,
        layer: int | list[int] | None = None,
        category: str = "selection",
        package: str | None = None,
        tags: list[str] | None = None,
        geometry: Any = None,
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
    ) -> ModelRegion:
        """Create and register a region from explicit cell IDs."""
        raw_cellids, cells, resolved_layer = normalize_region_inputs(cellids, layer=layer)
        region = ModelRegion(
            name=name,
            cells=cells,
            layer=resolved_layer,
            geometry=geometry,
            category=category,
            package=package,
            tags=[] if tags is None else list(tags),
            metadata={} if metadata is None else dict(metadata),
            cellids=raw_cellids if raw_cellids else list(cells),
        )
        return self.add(region, overwrite=overwrite)

    def add_from_geometry(
        self,
        name: str,
        geometry: Any,
        *,
        predicate: str = "intersects",
        layer: int | list[int] | None = None,
        category: str = "selection",
        package: str | None = None,
        tags: list[str] | None = None,
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
    ) -> ModelRegion:
        """Create and register a region from geometry intersection."""
        cell_series = get_vor_cells_as_series(
            self.model.vor.gdf_vorPolys,
            geometry,
            predicate=predicate,
        )
        if not isinstance(cell_series, pd.Series) or cell_series.empty:
            raise ValueError(f"Geometry for region {name!r} does not intersect any model cells")

        cells = sorted({int(cell) for group in cell_series.tolist() for cell in group})
        return self.add_from_cells(
            name,
            cells,
            layer=layer,
            category=category,
            package=package,
            tags=tags,
            geometry=geometry,
            metadata=metadata,
            overwrite=overwrite,
        )

    def resolve_cells(self, name: str) -> list[int]:
        """Resolve a region or group name into unique model cell IDs."""
        cells, _ = self.resolve_cells_with_trace(name)
        return cells

    def resolve_cells_with_trace(self, name: str) -> tuple[list[int], dict[int, list[str]]]:
        """Resolve a region/group and record which names contributed each cell."""
        resolved_cells: set[int] = set()
        trace: dict[int, set[str]] = {}
        self._resolve_into(name, resolved_cells, trace, stack=[])
        trace_out = {cell: sorted(names) for cell, names in trace.items()}
        return sorted(resolved_cells), trace_out

    def _resolve_into(
        self,
        name: str,
        resolved_cells: set[int],
        trace: dict[int, set[str]],
        *,
        stack: list[str],
    ):
        if name in self._regions:
            region = self._regions[name]
            for cell in region.cells:
                resolved_cells.add(cell)
                trace.setdefault(cell, set()).add(region.name)
            return

        group = self.get_group(name)
        if group.name in stack:
            cycle = " -> ".join([*stack, group.name])
            raise ValueError(f"Cycle detected while resolving groups: {cycle}")

        next_stack = [*stack, group.name]
        for member_name in group.members:
            self._resolve_into(member_name, resolved_cells, trace, stack=next_stack)

    def _group_contains(self, group_name: str, target_name: str) -> bool:
        group = self.get_group(group_name)
        for member_name in group.members:
            if member_name == target_name:
                return True
            if member_name in self._groups and self._group_contains(member_name, target_name):
                return True
        return False

    def _assert_no_cycle(self, group_name: str, member_name: str):
        if member_name in self._groups and self._group_contains(member_name, group_name):
            raise ValueError(
                f"Cannot add {member_name!r} to {group_name!r}; that would create a cyclic group relationship"
            )


def list_model_regions(model) -> pd.DataFrame:
    """Return the model's concrete region summary."""
    return model.regions.region_summary()


def list_model_groups(model) -> pd.DataFrame:
    """Return the model's region-group summary."""
    return model.regions.group_summary()


def get_region_cells(model, name: str) -> list[int]:
    """Resolve one region or group name to cell IDs."""
    return model.regions.resolve_cells(name)


def filter_region_heads(
    model,
    name: str,
    *,
    per: int | None = None,
    kstpkper: tuple | None = None,
    layer: int | list[int] | None = None,
) -> gpd.GeoDataFrame:
    """Filter ``all_heads`` down to one region or group.

    Parameters
    ----------
    model
        Model whose heads table should be filtered.
    name
        Region or group name to resolve through the model's registry.
    per, kstpkper
        Optional stress-period selectors.
    layer
        Optional layer filter applied after the region cells are resolved.
    """
    node = model.regions.get(name)
    resolved_cells = model.regions.resolve_cells(name)
    heads = model.all_heads.reset_index()
    heads = heads.loc[heads["cell"].isin(resolved_cells)].copy()

    target_layer = layer if layer is not None else getattr(node, "layer", None)
    if target_layer is not None:
        if isinstance(target_layer, list):
            heads = heads.loc[heads["layer"].isin(target_layer)]
        else:
            heads = heads.loc[heads["layer"] == int(target_layer)]

    if kstpkper is not None:
        heads = heads.loc[heads["kstpkper"] == tuple(kstpkper)]
    elif per is not None:
        heads = heads.loc[heads["kstpkper"].map(lambda value: value[1]) == int(per)]

    heads["geometry"] = heads["cell"].map(model.vor.gdf_vorPolys.geometry.to_dict())
    heads["region"] = name
    return gpd.GeoDataFrame(heads, geometry="geometry", crs=model.vor.crs)
