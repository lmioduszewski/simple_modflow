"""Named model-region and region-group helpers for simulation analysis."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy.modflow.mf6.grid.geometry import selected_cells
from myflopy.modflow.mf6.grid.selection import get_vor_cells_as_series
from myflopy.modflow.mf6.headsplus import _is_sentinel
from myflopy.modflow.mf6.package_explorer_utils import (
    CellNumbering,
    cell_numbering,
    check_view_layers,
    has_cell_numbering,
)
from myflopy.modflow.mf6.period_end import last_saved_step
from myflopy.modflow.utils.validators import whole_number


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
    numbering: CellNumbering | None = None,
) -> tuple[list[Any], list[int], int | list[int] | None]:
    """Normalize region cell inputs into raw IDs, cell numbers, and layers.

    Parameters
    ----------
    cellids : iterable of int or tuple, or None
        The cells. An integer is an in-layer cell number as given, and a one-item
        tuple ``(cell,)`` the same -- on a DISU grid that is MF6's own cellid
        ``(node,)``, which a DISU model's heads table and cell polygons number as
        cell ``node`` of layer 0 (measured, ledger 419). A longer tuple is the
        grid's own cellid, read through ``numbering`` (``CellNumbering``, the one
        rule the package tables use): ``(layer, row, col)`` on DIS, whose cell
        number is the row-major ``row * ncol + col``, and ``(layer, cell)`` on DISV;
        on DISU a longer tuple is refused naming ``(node,)``. With no grid known
        (``numbering=None``) a longer tuple reads as ``(layer, cell)``. A tuple of
        another length raises ``ValueError`` naming the form (ledger 376:
        ``(layer, row, col)`` used to read the ROW as the cell). An empty tuple is
        dropped. ``None`` is no cells.
    layer : int or list of int or None, optional
        The region's layer(s). ``None`` (the default) takes the layers the tuple
        cellids name -- one as an ``int``, several as a sorted list -- else none.
    numbering : CellNumbering or None, optional
        The model grid's numbering (``package_explorer_utils.cell_numbering``).
        Given, every layer (a tuple's and ``layer=``'s) is checked against its
        ``nlay`` and every cell against its ``ncpl`` -- a DIS row and column against
        ``nrow``/``ncol`` -- and one off the grid raises ``ValueError`` (ledger
        419: ``(5, 2, 3)`` on a 2-layer grid registered layer 5, and a bare
        cell past the layer was kept). ``None`` (the default) checks nothing.

    Returns
    -------
    tuple
        ``(raw_cellids, cells, layer)``: the cellids as given, the sorted unique
        in-layer cell numbers, and the resolved layer(s).

    Raises
    ------
    ValueError
        A tuple of the wrong length for the grid, or a layer or cell off it.
    """
    if cellids is None:
        raw_cellids: list[Any] = []
    else:
        raw_cellids = list(cellids)

    cells: list[int] = []
    derived_layers: list[int] = []
    if numbering is None:
        expected, form, grid = 2, "(layer, cell)", "a grid of unknown type"
    else:
        expected, form = numbering.length, numbering.form
        grid = {
            "structured": "a structured (DIS)",
            "vertex": "a vertex (DISV)",
            "unstructured": "an unstructured (DISU)",
        }[numbering.kind]

    for cellid in raw_cellids:
        if isinstance(cellid, tuple):
            if len(cellid) == 0:
                continue
            if len(cellid) == 1:
                cells.append(whole_number(cellid[0], what="cell"))
                continue
            if len(cellid) != expected:
                bare = "an in-layer cell number" if expected > 1 else "a bare node number"
                raise ValueError(
                    f"Region cellid {cellid!r} has {len(cellid)} items; on {grid} grid a "
                    f"cellid is {form}, or {bare}."
                )
            if numbering is None:
                derived_layers.append(int(cellid[0]))
                cells.append(int(cellid[1]))
                continue
            found_layer, found_cell = numbering.layer_cell(cellid)
            derived_layers.append(found_layer)
            cells.append(found_cell)
        elif isinstance(cellid, (list, np.ndarray)):
            raise TypeError(
                f"Region cellid {cellid!r} is a {type(cellid).__name__}; a cellid with "
                "a layer is a TUPLE, and a bare cell number an int."
            )
        else:
            # A whole number only: int() truncated 5.7 to cell 5 and read "5"
            # (ledger 478).
            cells.append(whole_number(cellid, what="cell"))

    if layer is not None:
        # A layer is a whole number too: 1.7 was stored as given (ledger 478).
        if isinstance(layer, (list, tuple, np.ndarray)):
            layer = [whole_number(value, what="layer") for value in layer]
        else:
            layer = whole_number(layer, what="layer")
    if numbering is not None:
        numbering.check(cells=cells)
        given = [] if layer is None else ([layer] if isinstance(layer, int) else layer)
        numbering.check(layers=given)

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
        """Create an empty region/group registry bound to ``model``."""

        self.model = model
        self._regions: dict[str, ModelRegion] = {}
        self._groups: dict[str, RegionGroup] = {}

    def __contains__(self, name: str) -> bool:
        """Whether ``name`` is a registered region or group."""

        return name in self._regions or name in self._groups

    def __getitem__(self, name: str) -> ModelRegion | RegionGroup:
        """The registered region or group named ``name``: ``model.regions["site"]``.

        The same lookup as :meth:`get`.

        Parameters
        ----------
        name : str
            A region's or a group's name -- the two share one namespace. An unknown
            name raises ``KeyError``.

        Returns
        -------
        ModelRegion or RegionGroup
            The registered object itself, not a copy.
        """

        return self.get(name)

    def add(self, region: ModelRegion, *, overwrite: bool = False) -> ModelRegion:
        """Register a region object built by hand.

        The object is stored as given, not copied, so a later edit to its ``cells``
        shows through every lookup and every group naming it. Nothing here checks
        it beyond its type -- not its cells against the model's grid, not its
        ``layer`` -- which is what :meth:`add_from_cells` and
        :meth:`add_from_geometry` do on the way to calling this.

        Parameters
        ----------
        region : ModelRegion
            The region, registered under ``region.name``. Anything else -- a
            :class:`RegionGroup` included, which is :meth:`add_group`'s -- raises
            ``TypeError`` (a group registered here broke :meth:`region_summary` and
            :meth:`summary` until ledger 478).
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a region of that name is
            registered; ``True`` replaces it, and every group naming it resolves to
            the new one. A name a GROUP holds raises ``ValueError`` either way.

        Returns
        -------
        ModelRegion
            ``region`` itself.

        Raises
        ------
        ValueError
            If ``region.name`` is a group's name, or an existing region's without
            ``overwrite=True``.
        """
        if not isinstance(region, ModelRegion):
            hint = " A group is add_group's." if isinstance(region, RegionGroup) else ""
            raise TypeError(
                f"RegionRegistry.add registers a ModelRegion, not a {type(region).__name__} "
                f"(ledger 478).{hint}"
            )
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
        """Create or replace a named group of regions/groups.

        A group holds region and group NAMES, resolved to cells only when read
        (:meth:`resolve_cells`), so a later change to a member shows through.

        Parameters
        ----------
        name : str
            Group name. A name a REGION holds raises ``ValueError`` whatever
            ``overwrite`` says.
        members : iterable of str or None, optional
            Region or group names to add, in order, each through
            :meth:`add_to_group`: an unknown name raises ``KeyError``, one that
            would make a cycle (``name`` itself, or a group already containing
            ``name``) ``ValueError``, and a repeated name is kept once. A refused
            member leaves the registry exactly as it was before the call -- no
            group ``name`` if there was none, the previous one under
            ``overwrite=True`` (ledger 375). Default ``None`` creates an empty
            group.
        tags : iterable of str or None, optional
            Labels, copied into a list on the group and shown by
            :meth:`group_summary`; nothing reads them. Any iterable is copied item
            by item, so a bare string is split into its characters
            (``tags="wells"`` stores five one-letter tags) -- pass a list. Default
            ``None`` stores none.
        metadata : dict or None, optional
            Free key/values, copied into a dict on the group; nothing reads them.
            Default ``None`` stores an empty dict.
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a group called ``name``
            exists; ``True`` replaces it with the new group once every member has
            been accepted. Groups that named the old one now resolve to the new.

        Returns
        -------
        RegionGroup
            The registered group.

        Raises
        ------
        ValueError
            If ``name`` is a region's name, an existing group's without
            ``overwrite=True``, or a member would make a cycle.
        KeyError
            If a ``members`` name is neither a region nor a group.
        """
        if name in self._regions:
            raise ValueError(f"Cannot add group {name!r}; a region with that name already exists")
        if not overwrite and name in self._groups:
            raise ValueError(f"Group {name!r} already exists")

        group = RegionGroup(
            name=name,
            members=[],
            tags=[] if tags is None else list(tags),
            metadata={} if metadata is None else dict(metadata),
        )
        # Every member is checked against the registry WITH the new group in it
        # (a cycle runs through it), and a refused one puts the registry back as
        # it was: the group used to stay registered half-filled, and under
        # overwrite=True the group it replaced was already gone (ledger 375).
        member_names = [] if members is None else list(members)
        previous = self._groups.get(name)
        self._groups[name] = group
        filled = False
        try:
            for member_name in member_names:
                self.add_to_group(name, member_name)
            filled = True
        finally:
            if not filled:
                if previous is None:
                    del self._groups[name]
                else:
                    self._groups[name] = previous

        return group

    def add_to_group(self, group_name: str, member_name: str) -> RegionGroup:
        """Add an existing region or group to another group.

        Parameters
        ----------
        group_name : str
            The group to add to. An unknown name, or a region's, raises
            ``KeyError``.
        member_name : str
            A registered region or group name; an unknown one raises ``KeyError``.
            ``group_name`` itself, or a group that already contains ``group_name``
            directly or through its subgroups (a cycle), raises ``ValueError``.
            Adding a member already present changes nothing.

        Returns
        -------
        RegionGroup
            The updated group.
        """
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
        """Remove one member from a named group.

        Only the membership goes; the member itself stays registered (:meth:`remove`
        deletes it).

        Parameters
        ----------
        group_name : str
            The group to remove from. An unknown name, or a region's, raises
            ``KeyError``.
        member_name : str
            The member to drop. A name that is not a member -- registered or not --
            is ignored without an error.

        Returns
        -------
        RegionGroup
            The updated group.
        """
        group = self.get_group(group_name)
        group.members = [member for member in group.members if member != member_name]
        return group

    def remove(self, name: str) -> ModelRegion | RegionGroup:
        """Remove a region or group and detach it from all groups.

        Parameters
        ----------
        name : str
            A region or group name; an unknown one raises ``KeyError``. The name is
            also dropped from the members of every group. Removing a group does not
            remove its members.

        Returns
        -------
        ModelRegion or RegionGroup
            The removed object.
        """
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
        """Return either a region or a group by name.

        Parameters
        ----------
        name : str
            A region's or a group's name -- the two share one namespace. An unknown
            name raises ``KeyError``.

        Returns
        -------
        ModelRegion or RegionGroup
            The registered object itself, not a copy.
        """
        if name in self._regions:
            return self._regions[name]
        if name in self._groups:
            return self._groups[name]
        raise KeyError(f"{name!r} does not exist")

    def get_region(self, name: str) -> ModelRegion:
        """Return one concrete region by name.

        Parameters
        ----------
        name : str
            A region's name. A group's name, or an unknown one, raises ``KeyError``
            (:meth:`get_group` is a group's, :meth:`get` either).

        Returns
        -------
        ModelRegion
            The registered region itself, not a copy.
        """
        if name not in self._regions:
            raise KeyError(f"Region {name!r} does not exist")
        return self._regions[name]

    def get_group(self, name: str) -> RegionGroup:
        """Return one region group by name.

        Parameters
        ----------
        name : str
            A group's name. A region's name, or an unknown one, raises ``KeyError``
            (:meth:`get_region` is a region's, :meth:`get` either).

        Returns
        -------
        RegionGroup
            The registered group itself, not a copy.
        """
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
        """Create and register a region from explicit cell IDs.

        The region is kept in this in-memory registry only (not saved with a
        project). ``model.add_region`` is this call.

        Parameters
        ----------
        name : str
            Region name. A name a GROUP holds raises ``ValueError`` whatever
            ``overwrite`` says.
        cellids : iterable of int or tuple, or None
            The cells. An ``int`` (a numpy integer too) is an in-layer cell number,
            and a one-item tuple ``(cell,)`` the same -- on DISU that is MF6's own
            cellid ``(node,)``, cell ``node`` of layer 0; an empty tuple is dropped
            and ``None`` is no cells. A longer tuple is the model grid's own
            cellid, read by the rule every package table uses (ledger 419):
            ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS, whose cell
            number is FloPy's row-major ``row * ncol + col``; on DISU a longer
            tuple is refused naming ``(node,)``. A tuple of another length raises
            ``ValueError`` naming the expected form (ledger 376), and so does any
            layer, row, column or cell off the grid -- a bare cell number
            included (ledger 419). The model's grid is read whenever a cellid or
            ``layer`` is given; a registry bound to no model, or a model with no
            discretization package, knows no grid, reads a longer tuple as
            ``(layer, cell)`` and checks nothing. Anything that is not a tuple
            must be a WHOLE NUMBER -- an integer, or an integral float such as
            ``5.0`` -- and anything else raises ``TypeError``: ``int()`` truncated
            ``5.7`` to cell 5 and read ``"5"`` in silence until ledger 478. A
            list ``[0, 5]`` or an array row raises ``TypeError`` -- a cellid with a
            layer must be a TUPLE.
            The region stores the sorted unique in-layer cell
            numbers, which is what :meth:`resolve_cells` returns and
            ``model.region_heads`` matches against the heads table's ``cell``
            column, and ``cellids`` as given.
        layer : int or list of int or None, optional
            Zero-based layer (or layers) stored on the region; ``model.region_heads``
            keeps only those layers when called without its own ``layer``. Nothing
            else reads it -- ``map(select=name)`` draws the region's cells on
            whatever layer the map shows. A layer off the model's grid raises
            ``ValueError`` (on DISU the only layer is 0; ledger 419), and one that
            is not a whole number ``TypeError`` (ledger 478). Default
            ``None`` takes the layers the
            tuple ``cellids`` name (one as an ``int``, several as a sorted list),
            else leaves the region layer-free. A value given here WINS over those
            layers without a word: ``[(0, 5), (1, 6)]`` with ``layer=2`` stores
            layer 2, the tuples' 0 and 1 surviving only in ``cellids``.
        category : str, optional
            Free label stored on the region and shown by :meth:`region_summary`;
            nothing reads it. Default ``"selection"``.
        package : str or None, optional
            A package name stored on the region as a label and shown by
            :meth:`region_summary`; nothing is looked up or checked. Default
            ``None``.
        tags : iterable of str or None, optional
            Labels, copied into a list on the region; nothing reads them. Any
            iterable is copied item by item, so a bare string is split into its
            characters (``tags="wells"`` stores five one-letter tags) -- pass a
            list. Default ``None`` stores none.
        geometry : object, optional
            Anything, kept on the region for reference; it selects no cells
            (:meth:`add_from_geometry` does that). Default ``None``.
        metadata : dict or None, optional
            Free key/values, copied into a dict on the region; nothing reads them.
            Default ``None`` stores an empty dict.
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a region called ``name``
            exists; ``True`` replaces it, and every group naming it resolves to the
            new cells.

        Returns
        -------
        ModelRegion
            The registered region.

        Raises
        ------
        ValueError
            If ``name`` is a group's, or an existing region's without
            ``overwrite=True``, or a tuple cellid is not of the grid's form, or a
            layer, row, column or cell is off the model's grid.
        TypeError
            If a cellid is a list or an array row (neither a tuple nor a number).
        """
        # ``None`` is no cells, as ``normalize_region_inputs`` reads it.
        cellids = [] if cellids is None else list(cellids)
        # The model's grid reads the tuples and bounds every layer and cell
        # (ledger 419). An unbound registry (``model=None``) or a model with
        # no discretization package knows no grid.
        needs_grid = bool(cellids) or layer is not None
        known = needs_grid and has_cell_numbering(self.model)
        raw_cellids, cells, resolved_layer = normalize_region_inputs(
            cellids,
            layer=layer,
            numbering=cell_numbering(self.model) if known else None,
        )
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
        predicate: str | None = None,
        layer: int | list[int] | None = None,
        category: str = "selection",
        package: str | None = None,
        tags: list[str] | None = None,
        metadata: dict[str, Any] | None = None,
        overwrite: bool = False,
    ) -> ModelRegion:
        """Create and register a region from geometry intersection.

        The cells are selected on the model's grid (``model.vor.gdf_vorPolys``,
        so in-layer cell numbers), flattened into one list and registered through
        :meth:`add_from_cells`, so every rule there applies.
        ``model.add_region_from_geometry`` is this call.

        Parameters
        ----------
        name : str
            Region name, as in :meth:`add_from_cells`.
        geometry : shapely geometry, GeoDataFrame, GeoSeries, str or pathlib.Path
            What selects the cells: any shapely geometry, a GeoDataFrame or
            GeoSeries (every feature's cells, pooled), or a path to a vector file
            (its first layer is read). One carrying a CRS is reprojected to the
            grid's; one carrying none is read in the grid's frame. Any other type
            raises ``ValueError``. Also stored on the region as given.
        predicate : str or None, optional
            Default ``None`` selects by MEASURE: a polygon reaches the cells it
            covers by area, a line those it runs through, a point the one cell
            holding it (ledgers 210-212 and 237); a feature reaching no cell is
            left out, named in a logged warning. A geopandas predicate
            (``"intersects"``, ``"contains"``, ...) takes the raw spatial join
            instead, which counts any touch -- a polygon that is exactly one cell
            selects its neighbours across shared faces and vertices too.
        layer : int or list of int or None, optional
            Zero-based layer (or layers) stored on the region, as in
            :meth:`add_from_cells`; the geometry selects in plan view and names no
            layer, so ``None`` (the default) leaves the region layer-free.
        category : str, optional
            Free label stored on the region. Default ``"selection"``.
        package : str or None, optional
            A package name stored on the region as a label; nothing is looked up.
            Default ``None``.
        tags : iterable of str or None, optional
            Labels, copied into a list, as in :meth:`add_from_cells` (a bare string
            is split into characters). Default ``None`` stores none.
        metadata : dict or None, optional
            Free key/values, copied into a dict on the region. Default ``None``
            stores an empty dict.
        overwrite : bool, optional
            Default ``False`` raises ``ValueError`` when a region called ``name``
            exists; ``True`` replaces it.

        Returns
        -------
        ModelRegion
            The registered region.

        Raises
        ------
        ValueError
            If ``geometry`` is of an unsupported type or reaches no cell, or for a
            name collision as in :meth:`add_from_cells`.
        """
        cell_series = get_vor_cells_as_series(
            self.model.vor.gdf_vorPolys,
            geometry,
            predicate=predicate,
        )
        if not isinstance(cell_series, pd.Series) or cell_series.empty:
            raise ValueError(f"Geometry for region {name!r} does not intersect any model cells")

        # One entry per feature, each a LIST of cells: one flattener.
        cells = selected_cells(cell_series)
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
        """Resolve a region or group name into unique model cell IDs.

        Parameters
        ----------
        name : str
            A region or group name; a group is resolved through its members,
            recursively. An unknown name raises ``KeyError`` saying it is neither
            (it said "Group" until ledger 478); a cycle raises ``ValueError``.

        Returns
        -------
        list of int
            The sorted union of the in-layer cell numbers the reached regions
            store. Layers are not applied: two regions on different layers sharing
            a cell number give it once.
        """
        cells, _ = self.resolve_cells_with_trace(name)
        return cells

    def resolve_cells_with_trace(self, name: str) -> tuple[list[int], dict[int, list[str]]]:
        """Resolve a region/group and record which names contributed each cell.

        Parameters
        ----------
        name : str
            A region or group name, resolved as in :meth:`resolve_cells`; an
            unknown name raises ``KeyError``, a cycle ``ValueError``.

        Returns
        -------
        tuple of (list of int, dict of int to list of str)
            The sorted unique cells, and for each cell the sorted names of the
            REGIONS that contributed it (never a group's name; a cell in several
            regions lists each).
        """
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
        """Recursively add a region/group's cells to ``resolved_cells`` (with per-cell provenance).

        ``stack`` tracks the resolution path to detect and reject cyclic group
        references.
        """

        if name in self._regions:
            region = self._regions[name]
            for cell in region.cells:
                resolved_cells.add(cell)
                trace.setdefault(cell, set()).add(region.name)
            return
        if name not in self._groups:
            # Not "Group ... does not exist": the name may have meant a region.
            raise KeyError(f"{name!r} is neither a region nor a group of this model")

        group = self.get_group(name)
        if group.name in stack:
            cycle = " -> ".join([*stack, group.name])
            raise ValueError(f"Cycle detected while resolving groups: {cycle}")

        next_stack = [*stack, group.name]
        for member_name in group.members:
            self._resolve_into(member_name, resolved_cells, trace, stack=next_stack)

    def _group_contains(self, group_name: str, target_name: str) -> bool:
        """Whether ``group_name`` contains ``target_name`` directly or transitively via subgroups."""

        group = self.get_group(group_name)
        for member_name in group.members:
            if member_name == target_name:
                return True
            if member_name in self._groups and self._group_contains(member_name, target_name):
                return True
        return False

    def _assert_no_cycle(self, group_name: str, member_name: str):
        """Raise if adding ``member_name`` to ``group_name`` would create a cyclic group relationship."""

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
    model : SimulationBase
        Model whose heads table should be filtered.
    name : str
        Region or group name to resolve through the model's registry.
    per : int or None, optional
        Zero-based stress period, read at its LAST saved step on the head file's
        clock (``period_end.last_saved_step``, ledger 337's rule; ledger 348): the
        same time ``kstpkper`` names, one per call. It kept EVERY saved step of the
        period until ledger 348. A period the head file saved nothing in raises
        ``ValueError``. Ignored when ``kstpkper`` is given.
    kstpkper : tuple of (int, int) or None, optional
        One zero-based saved time; takes precedence over ``per``. A time not in the
        file gives an empty frame.
    layer : int or list of int or None, optional
        Layer filter applied after the region cells are resolved; ``None`` uses the
        region's own layer, if any. A layer the model's views do not number -- ONE
        on a DISU -- raises ``ValueError`` naming the count (ledger 496); it
        filtered to an empty frame.

    Returns
    -------
    geopandas.GeoDataFrame
        The kept rows of ``all_heads`` with each cell's polygon and a ``region``
        column. A cell MF6 wrote as dry (-1e30) or inactive (+1e30) reads
        ``NaN``, its row kept (ledger 495): those numbers are not heads, and
        a mean over the region skips them.
    """
    node = model.regions.get(name)
    resolved_cells = model.regions.resolve_cells(name)
    heads = model.all_heads.reset_index()
    heads = heads.loc[heads["cell"].isin(resolved_cells)].copy()
    # MF6's no-value numbers -- +1e30 inactive, -1e30 DRY -- are NaN, the row kept,
    # by the head file's one rule (`headsplus._is_sentinel`, ledger 315): a region
    # mean averaged the 1e30 in (ledger 495).
    for column in heads.columns.difference(["kstpkper", "layer", "cell"]):
        values = pd.to_numeric(heads[column], errors="coerce")
        heads[column] = values.mask(_is_sentinel(values))

    if layer is not None:
        # A layer the model's views do not number is refused, not filtered to an
        # empty frame -- ONE layer on a DISU (`check_view_layers`, ledger 496).
        target_layer = check_view_layers(model.gwf, layer)
    else:
        stored = getattr(node, "layer", None)
        target_layer = (
            None if stored is None
            else [int(value) for value in stored] if isinstance(stored, (list, tuple))
            else [int(stored)]
        )
    if target_layer is not None:
        heads = heads.loc[heads["layer"].isin(target_layer)]

    if kstpkper is not None:
        heads = heads.loc[heads["kstpkper"] == tuple(kstpkper)]
    elif per is not None:
        # The head file's own clock (`model.kstpkper`), the one `all_heads` is read on.
        step = last_saved_step(model.kstpkper, per)
        heads = heads.loc[heads["kstpkper"].map(lambda value: tuple(map(int, value))) == step]

    heads["geometry"] = heads["cell"].map(model.vor.gdf_vorPolys.geometry.to_dict())
    heads["region"] = name
    return gpd.GeoDataFrame(heads, geometry="geometry", crs=model.vor.crs)
