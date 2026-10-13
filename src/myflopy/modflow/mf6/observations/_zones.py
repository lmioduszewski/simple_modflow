"""Cell-zone flow targets: the ONE implementation behind ``DrnFlowTargets`` and ``RivFlowTargets``.

A zone is a named set of ``(layer, cell)`` pairs whose boundary flow -- MF6's
``DRN`` or ``RIV`` term -- is summed into one series. Both families share every
rule here (ledger 534): how a zone is drawn and resolved, the
observation route's declared identity (ledger 344), the budget route's reader, and
what each refuses. A family is a thin subclass naming its MF6 quantity, its sign
rule and its defaults.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry.base import BaseGeometry

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.domain import model_domain
from myflopy.modflow.mf6.grid.geometry import (
    cell_shares,
    domain_active,
    miss_reason,
    model_cell_polygons,
)
from myflopy.modflow.mf6.observations._shared import (
    _aggregate_named_series,
    _check_inequality,
    _compare_without_values,
    _copy_frame,
    _declared_route,
    _default_compare_stats,
    _default_obs_csv_name,
    _in_grid_crs,
    _is_blank,
    _measured_periods,
    _named_package,
    _normalize_compare_time,
    _normalize_named_series_values,
    _observation_label,
    _ObservedQuantity,
    _owning_package,
    _packages_of_type,
    _placement_note,
    _read_declared_observations,
    _refuse_flows_leaving_the_frame,
    _refuse_keys_off_the_periods,
    _refuse_observation_kind,
    _refuse_reserved_target_names,
    _report_missed_features,
    _score_bounds,
    _wide_by_period,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.package_budget import build_budget_result_table
from myflopy.modflow.mf6.package_explorer_utils import (
    cell_numbering,
    has_cell_numbering,
    views_ncpl,
    views_nlay,
)
from myflopy.modflow.mf6.package_registry import package_name_of
from myflopy.modflow.mf6.period_end import period_ends

logger = get_logger(__name__)

#: What the budget route names as the file it reads.
_BUDGET_FILE = "the model's cell-budget file"


def _zone_layer(value, name, *, owner: str):
    """One zone row's layer: a zero-based int, or ``None`` when blank; anything else raises."""

    if value is None or (not isinstance(value, str) and pd.isna(value)) or str(value).strip() == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = float("nan")
    if not np.isfinite(number) or number < 0 or number != int(number):
        raise ValueError(
            f"{owner} zone {str(name)!r}: layer {value!r} is not a zero-based layer "
            "number (a whole number, 0 or more)."
        )
    return int(number)


def _with_geometry(frame: pd.DataFrame, *, owner: str, name_column: str = "name") -> pd.DataFrame:
    """A plain frame carrying a ``geometry`` column, as the ``GeoDataFrame`` it describes.

    A dict with a ``geometry`` key and a list of record dicts carrying one were
    read into a plain ``DataFrame``, whose geometry the column filter then dropped,
    so the zones were refused as having none ("polygon zones require a
    GeoDataFrame"). The geometry is taken in the GRID's coordinates: a plain dict
    carries no CRS.
    """

    if isinstance(frame, gpd.GeoDataFrame) or "geometry" not in frame.columns:
        return frame
    values = [
        None if value is None or (isinstance(value, float) and np.isnan(value)) else value
        for value in frame["geometry"]
    ]
    wrong = [
        (str(name), type(value).__name__)
        for name, value in zip(frame[name_column], values, strict=True)
        if value is not None and not isinstance(value, BaseGeometry)
    ]
    if wrong:
        listing = ", ".join(f"{name!r} ({kind})" for name, kind in wrong[:5])
        raise ValueError(
            f"{owner} locations: a 'geometry' entry must be a shapely geometry "
            f"(a Point, LineString or Polygon), not {listing}. Parse text first -- "
            "shapely.wkt.loads(...) -- or give the zones as a GeoDataFrame or vector file."
        )
    located = frame.copy()
    located["geometry"] = values
    return gpd.GeoDataFrame(located, geometry="geometry")


def _coerce_cells(values):
    """One zone row's cells as a list of ints; ``None`` when the row gives none.

    ``None`` (or ``NaN``) means "resolve this row from its geometry"; an empty
    list or empty text is an explicit zone of no cells.
    """

    if isinstance(values, str):
        return [int(value) for value in values.split(",") if str(value).strip() != ""]
    if isinstance(values, (list, tuple, np.ndarray, pd.Series)):
        return [int(value) for value in values]
    if values is None or pd.isna(values):
        return None
    return [int(values)]


def _normalize_zone_locations(
    locations,
    *,
    owner: str,
    name_column: str = "name",
    group_column: str | None = "group",
    weight_column: str | None = "weight",
) -> pd.DataFrame:
    """Normalize zone definitions from cells or geometry.

    Every row comes back with ``name``, ``group``, ``weight``, ``layer`` and
    ``cells`` when the input has that column: a list of ints, or ``None`` where the
    row gives none -- a row :func:`_resolve_zone_cells` resolves from its geometry.
    An EMPTY list (``[]``, or ``""`` in a CSV) is an explicit zone of no cells.

    ``layer`` is an int where the row gives one, and ``None`` where it is blank:
    the row then reads EVERY layer of its columns -- a row resolved from geometry
    the columns it covers (ledger 535, ledger 482's area rule), and a row giving
    ``cells`` those columns (ledger 550), both expanded by
    :func:`_resolve_zone_cells`. Either read layer 0 alone: a drain four layers
    deep under a seepage polygon read 1.8% of its flow, and the Ten Trails ``ss``
    south drains listed by column read 245 of their 935 records.
    """

    if isinstance(locations, pd.Series):
        frame = locations.rename("cells").rename_axis(name_column).reset_index()
        if frame.columns[0] != name_column:
            frame = frame.rename(columns={frame.columns[0]: name_column})
    elif isinstance(locations, dict):
        lower_keys = {str(key).strip().lower() for key in locations.keys()}
        if "name" in lower_keys or "cells" in lower_keys or "geometry" in lower_keys:
            frame = pd.DataFrame(locations)
        else:
            frame = pd.DataFrame({name_column: list(locations.keys()), "cells": list(locations.values())})
    elif isinstance(locations, (list, tuple)):
        frame = pd.DataFrame(locations)
    elif isinstance(locations, (pd.DataFrame, gpd.GeoDataFrame)):
        frame = locations.copy()
    else:
        frame = gpd.read_file(locations)

    if name_column not in frame.columns:
        raise ValueError(f"{owner} locations must include a name column.")
    frame = _with_geometry(frame, owner=owner, name_column=name_column)
    keep = [name_column]
    for column in ("cells", group_column, weight_column, "layer"):
        if column is not None and column in frame.columns:
            keep.append(column)
    if isinstance(frame, gpd.GeoDataFrame):
        keep.append(frame.geometry.name)
    frame = frame.loc[:, list(dict.fromkeys(keep))].copy()
    frame = frame.rename(columns={name_column: "name"})
    if group_column and group_column in frame.columns:
        frame = frame.rename(columns={group_column: "group"})
    if weight_column and weight_column in frame.columns:
        frame = frame.rename(columns={weight_column: "weight"})
    if "group" not in frame.columns:
        frame["group"] = pd.NA
    if "weight" not in frame.columns:
        frame["weight"] = np.nan
    if "cells" in frame.columns:
        frame["cells"] = pd.Series(
            [_coerce_cells(values) for values in frame["cells"]], index=frame.index, dtype=object
        )
    # A blank layer is ``None`` -- every layer -- whichever way the row names its cells.
    layers = frame["layer"] if "layer" in frame.columns else [None] * len(frame)
    resolved_layers = [
        _zone_layer(value, name, owner=owner)
        for value, name in zip(layers, frame["name"], strict=True)
    ]
    # Whole numbers as int64 once every row has one, so a snapshot written to a
    # GeoPackage reads back as numbers; ``None`` (every layer) keeps it object.
    every = any(layer is None for layer in resolved_layers)
    frame["layer"] = pd.Series(resolved_layers, index=frame.index, dtype=object if every else "int64")
    return frame


def _zone_listing(names) -> str:
    """Zone names for a message: the first ten, quoted, and a count of the rest."""

    names = list(dict.fromkeys(str(name) for name in names))
    listing = ", ".join(repr(name) for name in names[:10])
    if len(names) > 10:
        listing += f", and {len(names) - 10} more"
    return listing


def _layer_count(model, domain, ncpl: int) -> int:
    """How many layers the model's grid has: the domain's, else its node count over ``ncpl``."""

    values = None if domain is None else np.asarray(domain)
    if values is not None and values.ndim == 2:
        return int(values.shape[0])
    lookup = getattr(model, "node_to_lni", None)
    if lookup:
        return max(1, len(lookup) // max(int(ncpl), 1))
    nlay = getattr(model, "nlay", None)
    return int(nlay) if nlay else 1


def _is_every_layer(layer) -> bool:
    """Whether a normalized zone row's ``layer`` is blank: the row reads every layer."""

    return layer is None or (not isinstance(layer, (int, np.integer)) and pd.isna(layer))


def _grid_layers(model, *, owner_name: str) -> tuple[int, int]:
    """The model's layers and cells per layer, as its views number them.

    :func:`~myflopy.modflow.mf6.package_explorer_utils.views_nlay` and
    ``views_ncpl`` -- the rule the observation route writes cellids by (one layer
    of nodes on a DISU) -- else, for a model whose grid has no type, the
    ``(layer, cell)`` pairs of its ``node_to_lni``, the cells the budget route
    reads.

    Raises
    ------
    ValueError
        When the model offers neither, naming ``layer``.
    """

    if has_cell_numbering(model):
        return views_nlay(model), views_ncpl(model)
    lookup = getattr(model, "node_to_lni", None)
    if lookup:
        pairs = list(lookup.values())
        return 1 + max(int(layer) for layer, _cell in pairs), 1 + max(int(cell) for _layer, cell in pairs)
    raise ValueError(
        f"{owner_name}: model {getattr(model, 'name', '?')!r} has no grid whose layers "
        "can be read, so a zone given no `layer` cannot be read in every layer; give each "
        "zone row its `layer`."
    )


def _expand_listed_columns(frame: pd.DataFrame, model, *, owner_name: str) -> pd.DataFrame:
    """Each row LISTING ``cells`` with no ``layer`` as one row per layer: every layer of its columns.

    The rule a zone drawn as geometry with no layer follows (ledger 535,
    ledger 482's area rule), for a zone whose columns are LISTED (ledger
    550): one row per layer of the model, each holding the listed
    columns ACTIVE in that layer, since MF6 writes neither a budget entry nor an
    observation where there is no cell. It read layer 0 alone: the Ten Trails
    ``ss`` south drains listed by column read 245 of their 935 records
    (-33,300.5 ft3/d where every layer holds -2,613,784.4). A row with a
    ``layer`` keeps its pairs as given, an inactive cell included (ledger 343).

    Refused, naming the zones: a listed column the grid does not have (it would
    match nothing in any layer), and one ACTIVE IN NO LAYER -- a column every
    layer drops is a column no boundary can be in, so the zone would lose a cell
    its caller typed, in silence.

    Raises
    ------
    ValueError
        When ``model`` is ``None`` (the layers are the model's), when the model's
        layers cannot be read, or naming a listed column the grid does not have
        or one active in no layer.
    """

    listed = [
        cells is not None and _is_every_layer(layer)
        for cells, layer in zip(frame["cells"], frame["layer"], strict=True)
    ] if "cells" in frame.columns else [False] * len(frame)
    if not any(listed):
        return frame
    names = [name for name, every in zip(frame["name"], listed, strict=True) if every]
    if model is None:
        raise ValueError(
            f"{owner_name} zone(s) {_zone_listing(names)} list `cells` with no `layer`, "
            "which reads EVERY layer of those columns -- and the layers are the model's. "
            "Pass the model (it is None here), or give each row its `layer`: "
            "[{'name': ..., 'layer': 0, 'cells': [...]}]."
        )
    nlay, ncpl = _grid_layers(model, owner_name=owner_name)
    domain = model_domain(model, ncpl)
    off: dict[str, list[int]] = {}
    dead: dict[str, list[int]] = {}
    rows = []
    for every, (_index, row) in zip(listed, frame.iterrows(), strict=True):
        if not every:
            rows.append(row)
            continue
        name = str(row["name"])
        columns = [int(cell) for cell in row["cells"]]
        stray = [cell for cell in columns if not 0 <= cell < ncpl]
        if stray:
            off.setdefault(name, []).extend(stray)
            continue
        per_layer = [
            [cell for cell in columns if domain_active(domain, k, cell)] for k in range(nlay)
        ]
        held = {cell for cells in per_layer for cell in cells}
        unheld = [cell for cell in columns if cell not in held]
        if unheld:
            dead.setdefault(name, []).extend(unheld)
            continue
        for k, cells in enumerate(per_layer):
            entry = row.copy()
            entry["layer"] = k
            entry["cells"] = cells
            rows.append(entry)
    if off:
        raise ValueError(
            f"{owner_name}: zone(s) "
            + "; ".join(f"{name!r} {sorted(set(cells))[:10]}" for name, cells in list(off.items())[:10])
            + f" list cells this model's grid does not have: {nlay} layer(s) of {ncpl} "
            "cells, numbered from 0. A zone's `cells` are ZERO-based numbers within a layer."
        )
    if dead:
        raise ValueError(
            f"{owner_name}: zone(s) "
            + "; ".join(f"{name!r} {sorted(set(cells))[:10]}" for name, cells in list(dead.items())[:10])
            + " list columns active in NO layer of the model, so no boundary can be "
            "there; a zone given no `layer` reads each listed column in every layer it "
            "is active in. Remove those cells, or check the cell numbers."
        )
    if isinstance(frame, gpd.GeoDataFrame):
        expanded = gpd.GeoDataFrame(rows, geometry=frame.geometry.name, crs=frame.crs)
    else:
        expanded = pd.DataFrame(rows)
    expanded = expanded.reset_index(drop=True)
    expanded["cells"] = pd.Series(list(expanded["cells"]), index=expanded.index, dtype=object)
    layers = list(expanded["layer"])
    every = any(_is_every_layer(layer) for layer in layers)
    expanded["layer"] = pd.Series(
        [None if _is_every_layer(layer) else int(layer) for layer in layers],
        index=expanded.index, dtype=object if every else "int64",
    )
    return expanded


def _resolve_zone_cells(
    locations: pd.DataFrame, model, *, owner, owner_name: str, subject: str
) -> pd.DataFrame:
    """Attach explicit cell lists to each zone, resolving geometry against the grid.

    A row that gives ``cells`` keeps them, as given; a row that gives none -- no
    ``cells`` column, or a missing entry in one -- is resolved from its geometry,
    so one ``GeoDataFrame`` may mix the two (ledger 342: its geometry rows were
    skipped whenever a ``cells`` column existed, and came back with no cells).

    A zone COVERS the cells it is drawn over, by the shared rule
    (:func:`myflopy.modflow.mf6.grid.geometry.cell_shares`): a seepage polygon
    reaches the cells whose AREA it overlaps, a trace the cells it runs through,
    and a point the ONE cell holding it. The join here was
    ``predicate="intersects"``, which is true of a shared edge and of a shared
    corner, so two zones abutting on a column of faces BOTH claimed the cells on
    it, and every boundary there was calibrated twice.

    **A geometry row with a layer reads that layer; one with none reads EVERY
    layer** of the columns it covers (ledger 535, ledger 482's rule
    for an observation area): it comes back as one row per layer, each holding the
    covered columns active in that layer. It read layer 0 alone, in silence -- a
    seepage polygon over drains stacked in four layers read 1.8% of their flow.

    Only cells ACTIVE in a row's layer are kept, since MF6 writes neither a budget
    entry nor an observation where there is no cell. Membership is deliberately NOT
    restricted to cells that carry a boundary: a zone is geometry, resolved with no
    reference to the package, and a covered cell without one adds zero.

    A zone given as SEVERAL ROWS is ONE zone: each geometry row carries the union
    of the cells its zone's geometry rows cover IN ITS LAYER, and the zone has
    missed only when it holds no cell at all, which is also how the set is counted.

    Misses follow ledger 211. A zone set that lands nowhere raises, naming both
    extents and both CRSs; a partial miss warns, by name, once per set rather than
    once per call, because every consumer resolves the zones again.

    **A row LISTING ``cells`` with no layer reads every layer of those columns**
    too, by the same rule (:func:`_expand_listed_columns`, ledger
    550); a row listing them with a layer keeps its pairs as given.
    """

    def drawn(frame):
        """Which rows are to be resolved from geometry: those giving no ``cells``."""

        if "cells" in frame.columns:
            return frame["cells"].isna().to_numpy()
        return np.ones(len(frame), dtype=bool)

    frame = locations.copy()
    pending = drawn(frame)
    if pending.any():
        waiting = frame.loc[pending, "name"]
        if not isinstance(frame, gpd.GeoDataFrame):
            raise ValueError(
                f"{owner_name} zone(s) {_zone_listing(waiting)} give no cells and no geometry. "
                "Give each a 'cells' list, or give the zones as a GeoDataFrame, a vector file, "
                "or record dicts carrying a 'geometry'."
            )
        empty = frame.geometry.isna().to_numpy() | frame.geometry.is_empty.to_numpy()
        if (pending & empty).any():
            raise ValueError(
                f"{owner_name} zone(s) {_zone_listing(frame.loc[pending & empty, 'name'])} give "
                "neither cells nor a geometry, so there is nothing to resolve them from."
            )
        if model is None:
            raise ValueError(
                f"{owner_name} zone(s) {_zone_listing(waiting)} are drawn as geometry, which is "
                "resolved against a model's grid; pass the model (it is None here), or give "
                "their 'cells'."
            )
    frame = _expand_listed_columns(frame, model, owner_name=owner_name)
    pending = drawn(frame)
    if not pending.any():
        return frame
    polys = model_cell_polygons(model)
    domain = model_domain(model, len(polys))
    nlay = _layer_count(model, domain, len(polys))
    located = _in_grid_crs(frame, model)
    # (name, layer) -> covered columns; layer None is "every layer".
    covered: dict[tuple[str, object], list[int]] = {}
    reasons: dict[str, str] = {}
    for name, layer, geometry in zip(
        located.loc[pending, "name"], located.loc[pending, "layer"],
        located.geometry[pending], strict=True,
    ):
        every = _is_every_layer(layer)
        key = (str(name), None if every else int(layer))
        if every:
            def usable(cells):
                """Columns active in ANY layer: the area reads each layer it has."""

                return [
                    cell for cell in cells
                    if any(domain_active(domain, k, cell) for k in range(nlay))
                ]
            target = "column active in any layer"
        else:
            def usable(cells, layer=int(layer)):
                """Cells active in the row's own layer."""

                return [cell for cell in cells if domain_active(domain, layer, cell)]
            target = f"active cell in layer {int(layer)}"
        # Eligibility only. ``cell_shares`` sorts whatever ``usable`` returns, so
        # the lowest-numbered usable cell wins a tie without a second copy of that
        # rule here.
        hits = cell_shares(polys, geometry, usable=usable)
        if not hits:
            reasons.setdefault(key[0], miss_reason(polys, geometry, target=target))
        covered[key] = sorted(set(covered.get(key, ())) | {cell for cell, _share in hits})

    rows = []
    for position, (_index, row) in enumerate(located.iterrows()):
        if not pending[position]:
            entry = row.copy()
            entry["cells"] = list(row["cells"])
            rows.append(entry)
            continue
        layer = row["layer"]
        every = _is_every_layer(layer)
        if not every:
            entry = row.copy()
            entry["cells"] = list(covered[(str(row["name"]), int(layer))])
            rows.append(entry)
            continue
        columns = covered[(str(row["name"]), None)]
        for k in range(nlay):
            entry = row.copy()
            entry["layer"] = k
            # A list per row, not the shared one: two rows of one zone must not alias.
            entry["cells"] = [cell for cell in columns if domain_active(domain, k, cell)]
            rows.append(entry)
    resolved = gpd.GeoDataFrame(rows, geometry=located.geometry.name, crs=located.crs)
    resolved = resolved.reset_index(drop=True)
    resolved["cells"] = pd.Series(list(resolved["cells"]), index=resolved.index, dtype=object)
    resolved["layer"] = pd.Series([int(value) for value in resolved["layer"]], index=resolved.index, dtype="int64")
    # Counted by NAME, not by row: a multipart zone has one row per part (and an
    # every-layer zone one per layer), and counting rows said "1 of 2 zones land on
    # no active cell" of a single zone given as two rows.
    members = _zone_members(resolved)
    missed = [f"{name!r} ({reasons[name]})" for name, cells in members.items()
              if not cells and name in reasons]
    if missed:
        _report_missed_features(
            owner,
            missed,
            total=len(members),
            subject=subject,
            nothing_left="sum",
            where=lambda: _placement_note(frame, model, label="zones"),
            log=logger,
        )
    return resolved


def _zone_members(zones: pd.DataFrame) -> dict[str, set[tuple[int, int]]]:
    """Each zone's member cells by NAME, as ``(layer, cell)``: several rows are ONE zone.

    A cell is a ``(layer, cell)`` pair because that is what MF6 names: an
    observation's cellid, and the model-wide node a cell-budget record carries
    (``layer * ncpl + cell`` on a DISV). A zone keyed by its in-layer numbers
    alone read, below layer 0, whichever layer-0 boundary shared a number
    (ledger 342). A zone given as several rows -- explicit cells, or geometry, in
    one layer or several -- is united here, so its flow is summed once over all
    of them rather than once per row.

    Parameters
    ----------
    zones : pandas.DataFrame
        Zones with ``name``, ``layer`` and ``cells`` (``zone_definitions``); a row
        whose ``cells`` is ``None`` (unresolved geometry) or whose ``layer`` is
        blank (an unresolved every-layer row) adds nothing.

    Returns
    -------
    dict of str to set of tuple of (int, int)
        Zone name to its ``(layer, cell)`` pairs, zero-based, in the order the
        zones first appear.
    """

    members: dict[str, set[tuple[int, int]]] = {}
    for row in zones.itertuples(index=False):
        cells = getattr(row, "cells", None)
        entry = members.setdefault(str(row.name), set())
        layer = getattr(row, "layer", None)
        if not cells or _is_every_layer(layer):
            continue
        entry.update((int(layer), int(value)) for value in cells)
    return members


def _refuse_cells_off_the_grid(
    members: dict[str, set[tuple[int, int]]], node_to_lni, *, owner: str
) -> None:
    """Refuse a zone ``(layer, cell)`` the model's grid does not have, naming the zones.

    An explicit ``cells`` entry is taken as given, and on the budget-file route a
    pair the grid does not have -- a cell number past its layer's count, a
    negative one, a layer below the bottom -- matched no budget row, so the zone
    summed the cells that DID exist and said nothing: a zone of cells ``[20, 99]``
    on 36-cell layers read cell 20's flow as the zone's (ledger 343). MF6 refuses
    such a cellid on the observation route ("outside of the grid"), so the two
    routes agree. A cell the grid has but whose ``idomain`` removes is NOT refused:
    it holds no boundary and adds zero, as any covered cell without one does.

    Parameters
    ----------
    members : dict of str to set of tuple of (int, int)
        :func:`_zone_members`.
    node_to_lni : dict of int to tuple of (int, int)
        The model's ``node_to_lni``: every cell of its grid as ``(layer, cell)``.
    owner : str
        The target class, named in the error.

    Raises
    ------
    ValueError
        Naming each zone (the first five) with the pairs (its first five) the
        grid does not have, and the grid's layers and cells.
    """

    known = set(node_to_lni.values())
    stray = {
        name: sorted(pair for pair in pairs if pair not in known)
        for name, pairs in members.items()
    }
    stray = {name: pairs for name, pairs in stray.items() if pairs}
    if not stray:
        return
    per_layer = Counter(int(layer) for layer, _cell in known)
    counts = sorted(set(per_layer.values()))
    shape = (
        f"{len(per_layer)} layer(s) of {counts[0]} cells, numbered from 0"
        if len(counts) == 1
        else "layers " + ", ".join(f"{layer}: {per_layer[layer]} cells" for layer in sorted(per_layer))
    )
    listing = "; ".join(
        f"{name!r} {pairs[:5]}" + (f" and {len(pairs) - 5} more" if len(pairs) > 5 else "")
        for name, pairs in list(stray.items())[:5]
    )
    raise ValueError(
        f"{owner}: zone(s) name (layer, cell) pairs this model's grid does not "
        f"have, zero-based: {listing}. The grid has {shape}; correct the zone's "
        "`cells` or its `layer`."
    )


def _refuse_zones_holding_no_boundary(
    members: dict[str, set[tuple[int, int]]],
    names: list[str],
    budget: pd.DataFrame,
    *,
    model,
    quantity: _ObservedQuantity,
    package: str | None,
) -> None:
    """Refuse a zone whose cells hold no boundary of the packages read in ANY saved record.

    A zone empty in only SOME periods reads 0.0 there -- an absent boundary took no
    water (ledger 537). One empty in EVERY period the budget file saved is not a
    zone of zero flow but a zone in the wrong place: its cells in the wrong layer
    (a row given a ``layer`` reads that layer alone), or the boundaries in a
    package ``package=`` leaves out. It read 0.0 in every period,
    and a PEST build calibrated the target to a flow the model could never move
    (ledger 537, amended).

    Parameters
    ----------
    members : dict of str to set of tuple of (int, int)
        :func:`_zone_members`.
    names : list of str
        The zones to check: those with at least one member cell (a zone with none
        has no column and is refused by name where a target names it).
    budget : pandas.DataFrame
        :func:`_zone_budget_frame`: every saved record of the packages in scope,
        one row per boundary, with ``layer`` and ``cell``.
    model : SimulationBase
        The model read, whose packages in scope (:func:`_packages_in_scope`) the
        message names; asked only when a zone is refused.
    quantity : _ObservedQuantity
        The family's MF6 quantity, for the message.
    package : str or None
        The targets' ``package``, for the message.

    Raises
    ------
    ValueError
        Naming each empty zone (the first ten) with its layer(s), the packages
        read and ``package=``, and that a zone given no layer reads every layer.
    """

    held = set(zip(budget["layer"].astype(int), budget["cell"].astype(int), strict=True))
    empty = [name for name in names if not members[name] & held]
    if not empty:
        return
    scope = _packages_in_scope(model, quantity, package)
    text = quantity.package_type.upper()
    zones = ", ".join(
        f"{name!r} (layer{'s' if len(layers) > 1 else ''} {', '.join(str(layer) for layer in layers)})"
        for name in empty[:10]
        for layers in [sorted({int(layer) for layer, _cell in members[name]})]
    ) + (f", and {len(empty) - 10} more" if len(empty) > 10 else "")
    packages = ", ".join(repr(name.lower()) for name in scope)
    read = (
        f"package={package!r} ({text} package {packages})"
        if package is not None
        else f"package=None (every {text} package: {packages})"
    )
    other = (
        f"; package=None reads every {text} package"
        if package is not None
        else ""
    )
    raise ValueError(
        f"{quantity.owner}: zone(s) {zones} hold no {quantity.boundary} of {read} in any "
        f"period {_BUDGET_FILE} saved, so their flow would read 0.0 in every period -- "
        "a zone in the wrong layer or package, not a zone of zero flow. A zone given a "
        "`layer` reads that layer alone: check its `layer` -- or give none, to read every "
        "layer of its columns -- or check which "
        f"package holds its boundaries{other}."
    )


def _zone_labels(
    members: dict[str, set[tuple[int, int]]], *, owner: str
) -> dict[str, list[tuple[str, tuple[int, int]]]]:
    """Each zone's MF6 observation labels, one per member cell, and the cellid each names.

    ``<name>_c<cell>`` for a zone in one layer -- the label ``to_flopy_obs`` has
    always written -- and ``<name>_lay<layer>_c<cell>`` for a zone reaching
    several, whose cell numbers can repeat between layers (the head targets'
    ``_lay`` convention). Every character that is not a letter or digit becomes
    ``_``. Two zones whose labels coincide once written -- ``"a b"`` and ``"a_b"``
    -- are refused: MF6 accepts the duplicate and the reader could then only take
    one of the two columns for both.

    Parameters
    ----------
    members : dict of str to set of tuple of (int, int)
        :func:`_zone_members`.
    owner : str
        The target class, named in the error.

    Returns
    -------
    dict of str to list of tuple of (str, tuple of (int, int))
        Zone name to ``(label, (layer, cell))``, in ascending ``(layer, cell)``.

    Raises
    ------
    ValueError
        When two cells of the set would carry one label, naming the label and the
        zones.
    """

    labels: dict[str, list[tuple[str, tuple[int, int]]]] = {}
    owners: dict[str, list[str]] = {}
    for name, pairs in members.items():
        layered = len({layer for layer, _cell in pairs}) > 1
        entries = []
        for layer, cell in sorted(pairs):
            label = _observation_label(
                f"{name}_lay{layer}_c{cell}" if layered else f"{name}_c{cell}"
            )
            owners.setdefault(label.lower(), []).append(name)
            entries.append((label, (layer, cell)))
        labels[name] = entries
    clashes = {label: names for label, names in owners.items() if len(names) > 1}
    if clashes:
        listing = "; ".join(
            f"{label!r} from zones {_zone_listing(names)}" for label, names in list(clashes.items())[:5]
        )
        raise ValueError(
            f"{owner}: two zone cells would carry one MF6 observation label once "
            f"written (letters and digits only): {listing}. Rename one of the zones."
        )
    return labels


#: Where a reading's own weight waits while its values meet the zone table, which
#: carries a ``weight`` column of its own (ledger 504).
_READING_WEIGHT = "_reading_weight"


def _readings(values: pd.DataFrame) -> pd.DataFrame:
    """``values`` with a reading's own ``weight`` set aside before it meets the zone table."""

    return values.rename(columns={"weight": _READING_WEIGHT})


def _reading_weight_first(merged: pd.DataFrame) -> pd.DataFrame:
    """One ``weight``: a reading's own, else its zone's (ledger 504).

    The values table may carry a per-reading ``weight`` and the zone table always
    carries the zone's, so a plain merge returned ``weight_x`` and ``weight_y``
    from ``get`` and ``compare``, and ``residuals`` raised ``KeyError`` naming
    ``weight``. This is the order the PEST build applies
    (``pest.observations._zone_weights``).
    """

    if _READING_WEIGHT not in merged.columns:
        return merged
    reading = pd.to_numeric(merged.pop(_READING_WEIGHT), errors="coerce").astype(float)
    if "weight" in merged.columns:
        reading = reading.fillna(pd.to_numeric(merged["weight"], errors="coerce").astype(float))
    merged["weight"] = reading
    return merged


def _packages_in_scope(model, quantity: _ObservedQuantity, package: str | None) -> list[str]:
    """The NAMES of the packages a zone set reads on the budget route, upper-cased as MF6 writes them.

    The one named (``package``, any case; one the model does not have is refused),
    else every package of the family's type, found by TYPE (ledger 352).

    Raises
    ------
    ValueError
        When ``package`` names no package of the type, or the model has none.
    """

    if package is not None:
        return [package_name_of(_named_package(model, quantity, package)).strip().upper()]
    names = [package_name_of(found).strip().upper() for found in _packages_of_type(model, quantity.package_type)]
    if not names:
        raise ValueError(
            f"{quantity.owner}: model {getattr(model, 'name', '?')!r} has no "
            f"{quantity.package_type.upper()} package, so there is no "
            f"{quantity.boundary} flow to read."
        )
    return names


def _zone_budget_frame(model, quantity: _ObservedQuantity, package: str | None) -> pd.DataFrame:
    """The family's cell-budget records as a tidy frame: ``kstpkper``, ``per``, ``pname``, ``layer``, ``cell``, ``q``.

    Read through :func:`~myflopy.modflow.mf6.package_budget.build_budget_result_table`
    -- the reader ``model.<pkg>.results.q`` uses (ledger 352) -- so each record is
    matched to its term EXACTLY and to the package that wrote it by ``paknam2``, on
    the budget file's own clock (ledger 301), and each node is read as
    ``(layer, cell)`` through ``model.node_to_lni`` (ledger 342). It read
    ``model.bud(<type>).df``, which carries no package, so ``package=`` could not
    select one (ledger 534).

    **Every package in scope must have written a record** (ledger 536):
    a package declared without ``SAVE_FLOWS`` writes none, and its boundaries were
    left out of every zone in silence -- a zone over two drains, one in a package
    without it, read -339.52 where both removed about 679.

    Raises
    ------
    ValueError
        When ``package`` names no package of the type, the model has none, a
        package in scope wrote no record to the cell-budget file (naming it and
        ``save_flows=True``), or the file names a node the grid does not have.
    """

    scope = _packages_in_scope(model, quantity, package)
    reader = model._get_budget_reader()
    headers = reader.recordarray
    text = quantity.package_type.upper()
    written = {
        bytes(owner).decode().strip().upper()
        for name, owner in zip(headers["text"], headers["paknam2"], strict=True)
        if bytes(name).decode().strip().upper() == text
    }
    unsaved = [name for name in scope if name not in written]
    if unsaved:
        raise ValueError(
            f"{quantity.owner}: {text} package(s) {', '.join(repr(name.lower()) for name in unsaved)} "
            f"wrote no record to {_BUDGET_FILE}, so the {quantity.boundary} flow in "
            "their cells cannot be read from it and the zones would leave it out in "
            "silence. Declare them with save_flows=True"
            + (
                f", or name the package to read: {quantity.owner}(..., package=<its name>)."
                if package is None and len(scope) > len(unsaved)
                else "."
            )
        )
    return build_budget_result_table(
        model, budget_text=text, package_name=text.lower(), packages=scope
    )


@dataclass
class _CellZoneFlowTargets:
    """The shared body of the cell-zone flow families (``DrnFlowTargets``, ``RivFlowTargets``).

    A subclass names, as class variables, the MF6 quantity it reads
    (``_QUANTITY``), the sign rule it refuses by (``_SIGN``: keyword arguments of
    ``_refuse_flows_leaving_the_frame``, or ``None`` for a two-signed flow), the
    zone noun used in messages (``_SUBJECT``), its default observation file names
    and its plot title. Its class docstring documents the fields, which are
    declared here once.
    """

    _QUANTITY: ClassVar[_ObservedQuantity]
    _SIGN: ClassVar[dict | None] = None
    _SUBJECT: ClassVar[str] = "zone"
    _DEFAULT_CSV: ClassVar[str] = "zone_flow_targets.csv"

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    name_column: str = "name"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "flow"
    times: list | tuple | pd.Index | pd.Series | None = None
    package: str | None = None
    inequality: str | None = None

    def __post_init__(self):
        """Normalize the zone definitions and observed flows into internal tables."""

        owner = type(self).__name__
        _check_inequality(self.inequality, owner=owner)
        self._locations = _normalize_zone_locations(
            self.locations,
            owner=owner,
            name_column=self.name_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
        )
        _refuse_reserved_target_names(
            self._locations["name"], owner=owner, columns=("flow_target", "sim_flow"),
        )
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="flow_target",
            times=self.times,
            owner=owner,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations["name"].unique()) != 1:
                raise ValueError(f"{owner} values without explicit names require exactly one zone.")
            self._values["name"] = str(self._locations["name"].iloc[0])
        if self._SIGN is not None:
            _refuse_flows_leaving_the_frame(self._values, "flow_target", owner=owner, **self._SIGN)

    @property
    def locations_gdf(self):
        """A copy of the normalized zone-definition table (geometry preserved)."""

        return self._locations.copy()

    def zone_definitions(self, model) -> pd.DataFrame:
        """Resolve each zone's member cells against ``model`` (covering geometry if needed).

        A row given ``cells`` AND a ``layer`` keeps its pairs as given. A row given
        ``cells`` with no layer reads EVERY layer of those columns, one row per
        layer, each holding the columns active in that layer (ledger
        550). A row given no cells is resolved from its geometry: it
        reaches the cells it COVERS -- the rule ``.gpkg`` boundaries and head
        targets use -- active in its own layer, or, when the row gives no layer,
        in EVERY layer of the columns it covers, one row per layer (ledger 535);
        it comes back in the GRID's CRS, whatever the zone file's own was.

        Parameters
        ----------
        model : SimulationBase or None
            The model whose grid (``model.vor``), layers and active domain the
            rows without a layer, and the geometry rows, are resolved against.
            Not read when every row gives its ``cells`` AND its ``layer``, which
            are returned as given; ``None`` then works.

        Returns
        -------
        pandas.DataFrame or geopandas.GeoDataFrame
            One row per ``locations`` row -- and per LAYER for a row that gives
            none: ``name``, ``cells`` (zero-based cell numbers in the row's layer:
            a row with a layer, its own; a listed row with none, its columns
            active in that layer; a geometry row's, the union over every geometry
            row of that name IN THAT LAYER), ``group``, ``weight``, ``layer`` (an
            int on every row), and, where the zones carry geometry, the geometry,
            converted to the grid's CRS when both CRSs are known.

        Raises
        ------
        ValueError
            When a row gives no cells and has no geometry to resolve them from,
            or ``model`` is ``None`` while a row needs it (a geometry row, or
            ``cells`` with no layer); when a row listing ``cells`` with no layer
            names a column the grid does not have, or one active in no layer;
            or when no zone lands on an active cell, naming the zones (the first
            ten) and both extents and CRSs (ledger 211).

        Warns
        -----
        UserWarning
            When SOME zones land on no active cell, naming them, once per missed
            set on this object -- through ``warnings.warn`` as well as the
            logger, because a log line alone reaches nobody by default (ledgers
            288 and 342).
        """

        return _resolve_zone_cells(
            self._locations, model, owner=self, owner_name=type(self).__name__,
            subject=self._SUBJECT,
        )

    def get(self, model=None) -> pd.DataFrame:
        """Return the zone definitions (cells resolved when ``model`` is given), merged with values.

        Parameters
        ----------
        model : SimulationBase or None, optional
            When given, zones are resolved through :meth:`zone_definitions`,
            whose ``Raises`` apply. Default ``None``: the zones as normalized at
            construction -- a geometry row's ``cells`` is ``None``, unresolved,
            any row's ``layer`` is ``None`` when it gives none (every layer, not
            yet expanded), and geometry is in the zone file's own CRS.

        Returns
        -------
        pandas.DataFrame
            The zone table alone when there are no ``values``; otherwise the
            observed values -- ``time``, ``name``, ``flow_target`` -- joined on
            ``name`` to the zone table's columns other than geometry: one row
            per value and zone ROW, so a zone given as several rows repeats each
            of its values, and a value naming no zone gets one row of ``NaN``.
            ONE ``weight`` column: a reading's own weight where the values table
            gave one, else its zone's (ledger 504) -- the weight the PEST
            build gives the observation, ``NaN`` meaning 1.0 there.
        """

        definitions = self._locations.copy()
        if model is not None:
            definitions = self.zone_definitions(model)
        if self._values.empty:
            return definitions
        merged = _readings(self._values).merge(
            definitions.drop(columns="geometry", errors="ignore"), on="name", how="left"
        )
        return _copy_frame(_reading_weight_first(merged))

    def to_long(self, model=None) -> pd.DataFrame:
        """Return the long-format zone-flow target table (the same table as :meth:`get`).

        Parameters
        ----------
        model : SimulationBase or None, optional
            Passed to :meth:`get`: resolves geometry zones when given. Default
            ``None``.

        Returns
        -------
        pandas.DataFrame
            :meth:`get`'s table.
        """

        return self.get(model=model)

    def to_wide(self) -> pd.DataFrame:
        """Return the observed flows in wide format with one column per zone.

        One row per stress period a ``time`` names -- ``time`` is the zero-based
        period as ``int``, read as :meth:`compare` and the PEST build read it, so
        ``"1"`` and ``1`` are one row (ledger 422) -- and one column per target
        name holding its observed flow, in model flow units and MF6's sign (ledger
        350); a set given no values gives ``time`` and one empty column per
        target. Two rows HOLDING a value for one target in one stress period are
        REFUSED -- names compared stripped and lower-cased -- by the check
        :meth:`compare` and the PEST build read (ledgers 340, 371). A blank row is
        unmeasured and never collides; a target whose every row is blank keeps its
        column, all blank (ledger 422).

        Returns
        -------
        pandas.DataFrame
            ``time`` (``int64``, the period) plus one float column per target
            name, sorted.

        Raises
        ------
        ValueError
            When a ``time`` label is not a zero-based stress period -- not a
            number, or one with a fraction (``1.0`` is period 1) -- or two rows
            hold a value for one target in one period, naming the target, the
            period and the count.
        """

        caller = f"{type(self).__name__}.to_wide()"
        if self._values.empty:
            return _wide_by_period(
                self._values, value_column="flow_target", caller=caller,
                names=self._locations["name"],
            )
        _measured_periods(self._values, value_column="flow_target", caller=caller)
        # After the check a (target, time) cell holds at most ONE value, so
        # "first" only prefers it over blank rows.
        return _wide_by_period(self._values, value_column="flow_target", caller=caller)

    def summary(self, model=None) -> pd.DataFrame:
        """Return a compact summary (zone count, resolved cell count, row count).

        Parameters
        ----------
        model : SimulationBase or None, optional
            When given, zones are resolved through :meth:`zone_definitions`
            first, whose ``Raises`` apply. Default ``None``: geometry rows and
            rows with no ``layer`` are unresolved and count no cells.

        Returns
        -------
        pandas.DataFrame
            One row: ``n_zones`` (distinct names), ``n_cells`` (distinct
            ``(layer, cell)`` pairs over the whole set -- a cell is counted once
            however many rows or zones name it) and ``n_rows`` (observed values).
        """

        definitions = self.zone_definitions(model) if model is not None else self._locations.copy()
        zone_count = int(definitions["name"].nunique())
        cells = set().union(*_zone_members(definitions).values()) if len(definitions) else set()
        return pd.DataFrame([{"n_zones": zone_count, "n_cells": len(cells), "n_rows": int(len(self._values))}])

    def _zone_flow_table(self, model, zones: pd.DataFrame) -> tuple[pd.DataFrame, list[int]]:
        """Each zone's flow at each stress period's END, from the model's budget file.

        A zone's flow is the sum of its cells' budget rows in ONE time step
        (``kind="extensive"`` -- a zone IS that sum), read at the record that ends
        each period on the budget file's own clock (``SAVE BUDGET``, ledger 301),
        through the one period-end rule every target reads (ledger 336). A budget
        row belongs to a zone by its ``(layer, cell)`` (ledger 342), from the
        packages in scope -- the one ``package`` names, else every package of the
        type (:func:`_zone_budget_frame`).

        **A period in which a zone's cells hold no boundary reads 0.0** (ledger
        537), as the observation route reads it: the flow is extensive,
        and a boundary that is absent took no water. It read ``NaN``, which the
        PEST build wrote as an empty field that failed every run's instruction
        read. **A zone whose cells hold no boundary in ANY saved record is
        refused** (:func:`_refuse_zones_holding_no_boundary`, ledger 537 as
        amended): it came back with no rows, and the build died on a bare
        ``IndexError``, then read 0.0 in every period -- a zone in the wrong layer
        or package, calibrated against a flow no parameter moves.

        Returns
        -------
        tuple of (pandas.DataFrame, list of int)
            ``per``, ``name`` and ``sim_flow`` (MF6's sign), one row per zone
            WITH cells per period the file saved; and those periods.

        Raises
        ------
        ValueError
            As :func:`_zone_budget_frame`, :func:`_refuse_cells_off_the_grid` and
            the period-end rule; and when a zone with cells holds no boundary of
            the packages read in any record the budget file saved.
        """

        budget = _zone_budget_frame(model, self._QUANTITY, self.package)
        members = _zone_members(zones)
        _refuse_cells_off_the_grid(members, model.node_to_lni, owner=type(self).__name__)
        ends = period_ends(
            model,
            reader=lambda: model._get_budget_reader(),
            file=_BUDGET_FILE,
            saves="BUDGET",
            fallback=budget["kstpkper"],
        )
        periods = sorted(int(per) for per in ends)
        names = [name for name, pairs in members.items() if pairs]
        if not names or not periods:
            return pd.DataFrame(columns=["per", "name", "sim_flow"]), periods
        _refuse_zones_holding_no_boundary(
            members, names, budget, model=model, quantity=self._QUANTITY, package=self.package,
        )
        located = pd.MultiIndex.from_arrays([budget["layer"], budget["cell"]])
        pieces = [
            budget.loc[located.isin(sorted(members[name]))].assign(zone=index)
            for index, name in enumerate(names)
        ]
        rows = pd.concat(pieces, ignore_index=True)
        table = _aggregate_named_series(
            rows, id_column="zone", value_column="q", kind="extensive", ends=ends
        )
        full = pd.MultiIndex.from_product([periods, range(len(names))], names=["per", "zone"])
        table = (
            table.set_index(["per", "zone"])["q"].reindex(full).fillna(0.0).reset_index()
        )
        table["name"] = [names[int(zone)] for zone in table["zone"]]
        return table.rename(columns={"q": "sim_flow"}).loc[:, ["per", "name", "sim_flow"]], periods

    def _zone_details(self, zones: pd.DataFrame) -> pd.DataFrame:
        """One row per zone NAME: its ``group``, ``weight`` and ``cells`` as ``(layer, cell)``.

        What :meth:`compare` joins to each simulated row, on both routes. A zone
        given as several rows takes its ``group`` and ``weight`` from the first
        row that carries one, as the PEST build does
        (``pest.observations._zone_groups`` / ``_zone_weights``, which refuse rows
        that DISAGREE; ledger 359).
        """

        rows = [
            {
                "name": str(row.name),
                "group": getattr(row, "group", pd.NA),
                "weight": getattr(row, "weight", np.nan),
            }
            for row in zones.itertuples(index=False)
        ]
        group: dict[str, object] = {}
        weight: dict[str, float] = {}
        for row in rows:
            if not _is_blank(row["group"]):
                group.setdefault(row["name"], row["group"])
            number = pd.to_numeric(pd.Series([row["weight"]]), errors="coerce").iloc[0]
            if pd.notna(number):
                weight.setdefault(row["name"], float(number))
        details = pd.DataFrame(rows, columns=["name", "group", "weight"]).drop_duplicates(subset="name")
        details["group"] = pd.Series(
            [group.get(name, pd.NA) for name in details["name"]], index=details.index, dtype=object
        )
        details["weight"] = [weight.get(name, np.nan) for name in details["name"]]
        members = _zone_members(zones)
        details["cells"] = pd.Series(
            [sorted(members[name]) for name in details["name"]], index=details.index, dtype=object
        )
        return details.reset_index(drop=True)

    def _observed(self, model) -> tuple[pd.DataFrame, pd.DataFrame | None]:
        """The resolved zones, and the observation route's wide table (``None`` for the budget file).

        Each zone cell's label is the one ``to_flopy_obs`` writes, from the one
        :func:`_zone_labels`, and it must be declared in the package's own OBS as
        the family's observation OF THAT CELL (ledger 344); a zone reads the sum
        of its cells' columns in each time step, one row per stress period, at its
        end (ledger 346).

        Raises
        ------
        ValueError
            As :meth:`zone_definitions` and :func:`_zone_labels`; and as the
            observation route refuses (:meth:`simulated_series`).
        FileNotFoundError
            When a declared observation file is not in the workspace.
        """

        zones = self.zone_definitions(model)
        quantity = self._QUANTITY
        labels = _zone_labels(_zone_members(zones), owner=quantity.owner)
        wanted = {
            f"{zone_name}__{index}": (zone_name, label, cellid)
            for zone_name, entries in labels.items()
            for index, (label, cellid) in enumerate(entries)
        }
        chosen = _declared_route(model, quantity, wanted, package=self.package)
        if chosen is None:
            return zones, None
        observed = _read_declared_observations(model, chosen, owner=quantity.owner)
        result = pd.DataFrame({"time": observed["time"]})
        for zone_name, entries in labels.items():
            keys = [f"{zone_name}__{index}" for index in range(len(entries))]
            result[zone_name] = observed.loc[:, keys].sum(axis=1) if keys else 0.0
        return zones, result

    def _targets_by_period(self, caller: str) -> pd.DataFrame:
        """The values with ``per`` read from ``time`` by the one period rule (ledger 340)."""

        targets = self._values.copy()
        targets["per"] = _measured_periods(targets, value_column="flow_target", caller=caller)
        return targets

    def _compare_from_obs(self, model) -> pd.DataFrame | None:
        """Build the target-vs-simulated compare table from observation output, or ``None``."""

        zones, wide = self._observed(model)
        if wide is None:
            return None

        details = self._zone_details(zones)
        simulated = wide.melt(id_vars=["time"], var_name="name", value_name="sim_flow")
        simulated["per"] = simulated["time"]
        simulated = details.merge(simulated, on="name", how="left")
        if self._values.empty:
            compare = _compare_without_values(simulated, target_column="flow_target")
        else:
            caller = f"{type(self).__name__}.compare(...)"
            targets = self._targets_by_period(caller)
            _refuse_keys_off_the_periods(targets, wide["time"], model=model, caller=caller)
            compare = _reading_weight_first(
                _readings(targets).merge(simulated, on=["name", "per"], how="left")
            )
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        # A bound scores 0 where it holds, on both routes (ledger 564).
        return _copy_frame(_score_bounds(compare, self.inequality))

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-flow table: a ``time`` column and one column per zone.

        Read from the package's own observations (the family's MF6 type, ``DRN``
        or ``RIV``) when EVERY zone cell's label (as ``to_flopy_obs`` writes it) is
        declared there, for that cell -- only the files that OBS names are read,
        each record from its own column, one row per stress period at the time
        step that ENDS it (ledger 346), the zone the sum of its cells in that step
        -- else, when NO zone cell is declared, from the model's cell-budget file
        (ledger 344).

        **On the budget-file route a period holds ONE flow per zone: the flow at
        the period's END** (ledger 336), summed over the zone's cells in that step,
        from the package ``package`` names or, when it names none, from every
        package of the family's type -- each matched by the name MF6 writes on its
        records (ledger 534). A period in which the zone's cells hold no
        boundary reads 0.0, as the observation route reads it (ledger
        537); a zone whose cells hold none in ANY record the file saved is
        refused (ledger 537, amended), naming its layer(s) and ``package=``: it is
        a zone in the wrong layer or package -- a zone given a ``layer`` reads that
        layer alone -- not a zone of zero flow. A period whose
        final step was not saved is refused, because its end-of-period flow does
        not exist in the file; so is a package in scope that wrote no record
        (ledger 536).

        Parameters
        ----------
        model : SimulationBase
            A run model with a package of the family's type, as
            ``mf.load_mf6_run(ws)`` or ``run.model()`` returns; zones drawn as
            geometry are resolved against its grid. The observation route reads
            the package's declared OBS and the files it names in the model's
            workspace.

        Returns
        -------
        pandas.DataFrame
            ``time`` (the zero-based stress period) plus one float column per
            zone, in MF6's sign: water leaving the aquifer is NEGATIVE. On the
            budget-file route a zone with no member cell gets no column; on the
            observation route it reads 0.0.

        Raises
        ------
        FileNotFoundError
            On the observation route, when a file the OBS declares is not in the
            workspace.
        ValueError
            When a zone cell's label is declared in the package's OBS as another
            quantity only, for another cell or twice; when some zone cells of the
            set are declared and others are not; when the model carries several
            packages of the type, a zone cell is declared in them and ``package``
            does not name one; when ``package`` names no package of the type; or
            when a declared file's columns do not line up with its records (ledger
            344), or it does not end every stress period (ledger 346). On the
            budget-file route, when a package in scope wrote no record (``save_flows``
            off), a period the budget file saved lacks its final step -- naming the
            period, the steps saved and ``SAVE BUDGET LAST`` -- the file names a
            node the model's grid does not have, a zone names a
            ``(layer, cell)`` it does not have (ledger 343), or a zone's cells
            hold no boundary of the packages read in any record the file saved,
            naming the zone, its layer(s) and ``package=`` (ledger 537, amended);
            or as :meth:`zone_definitions` (no zone landing on an active cell,
            ledger 211) and ``to_flopy_obs`` (two zones sharing an observation
            label).
        """

        zones, observed = self._observed(model)
        if observed is not None:
            return observed

        frame, periods = self._zone_flow_table(model, zones)
        if frame.empty:
            return pd.DataFrame({"time": pd.Series(periods, dtype="int64")})
        # Columns in pivot order (sorted by name), as every template built so far
        # wrote its base table: a PEST instruction file reads by POSITION.
        wide = frame.pivot_table(index="per", columns="name", values="sim_flow", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare observed against simulated per-zone flow (with residuals).

        ``sim_flow`` comes from the same two routes as :meth:`simulated_series`:
        one flow per zone per stress period, the flow at the period's END (ledger
        336). Targets are matched on ``per``, read from their ``time`` as a
        zero-based integer -- the stress period, on both routes (ledger 346). A key
        naming no period the route holds is refused on both routes, naming the
        period each simulated-time key ends (it read ``NaN`` on the budget route and
        ``stats()`` dropped it in silence, ledger 538).

        Parameters
        ----------
        model : SimulationBase
            A run model with a package of the family's type, as
            ``mf.load_mf6_run(ws)`` or ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row -- one per zone and
            period, on both routes -- when the targets hold no values): ``name``,
            ``time``, ``per``, ``flow_target``, ``sim_flow``, ``residual``
            (simulated minus target), ``abs_residual``, and the zone's ``group``,
            ``weight`` and ``cells`` -- ONE set per zone however many rows define
            it, its cells the sorted ``(layer, cell)`` pairs MF6 reads; ``weight``
            is a reading's own where the values table gave one, else the zone's
            (ledger 504). ``flow_target`` and ``sim_flow`` are both in MF6's sign
            (ledger 350). For a BOUND set (``inequality``) ``residual`` is
            PEST++'s, 0 where the bound holds, and ``raw_residual`` after it is
            simulated minus target (ledger 564). With no values ``flow_target`` and the residuals are
            ``NaN``, and ``time`` is the period, equal to ``per``.

        Raises
        ------
        FileNotFoundError
            As :meth:`simulated_series`: a declared observation file is missing.
        ValueError
            As :meth:`simulated_series`; when a target's key is no stress period
            the route holds -- on the observation route a period of the model
            (ledger 346), on the budget route a period the budget file saved; or
            when a target's ``time`` is not a whole number, or two rows hold a
            value for one target in one period (ledger 340).
        """

        observed = self._compare_from_obs(model)
        if observed is not None:
            return observed

        zones = self.zone_definitions(model)
        simulated, periods = self._zone_flow_table(model, zones)
        simulated = simulated.merge(self._zone_details(zones), on="name", how="left")
        if self._values.empty:
            compare = _compare_without_values(simulated, target_column="flow_target")
        else:
            caller = f"{type(self).__name__}.compare(...)"
            targets = self._targets_by_period(caller)
            _refuse_keys_off_the_periods(
                targets, periods, model=model, caller=caller, route=_BUDGET_FILE,
            )
            compare = _reading_weight_first(
                _readings(targets).merge(simulated, on=["name", "per"], how="left")
            )
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        # A bound scores 0 where it holds, on both routes (ledger 564).
        return _copy_frame(_score_bounds(compare, self.inequality))

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table (target/sim/residual, with group/weight) for the zones.

        Parameters
        ----------
        model : SimulationBase
            A run model with a package of the family's type; passed to
            :meth:`compare`, whose routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``group``, ``time``,
            ``per``, ``flow_target``, ``sim_flow``, ``residual``,
            ``abs_residual``, ``weight`` and ``cells`` (each zone's
            ``(layer, cell)`` pairs; absent only when :meth:`compare` has none),
            and ``raw_residual`` for a bound set (ledger 564).
            A set given no ``values`` -- a capture template -- gets one row per
            simulated row, its targets and residuals ``NaN`` (ledger 342).
        """

        frame = self.compare(model)
        columns = ["name", "group", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual", "weight"]
        if "cells" in frame.columns:
            columns.append("cells")
        if "raw_residual" in frame.columns:
            columns.append("raw_residual")
        return frame.loc[:, columns]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the zone-flow target set.

        Parameters
        ----------
        model : SimulationBase
            A run model with a package of the family's type; passed to
            :meth:`compare`, whose routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual), ``mean_error``, ``mae`` and
            ``rmse`` of simulated minus target, in model flow units (length
            cubed per time). ``n`` is 0 and the rest ``NaN`` when no row has a
            residual, as with no ``values``.
        """

        return _default_compare_stats(self.compare(model))

    def _to_flopy_obs(self, model, *, csv_name: str, kind: str) -> dict:
        """The ``continuous`` dict ``to_flopy_obs`` documents, for the family's quantity."""

        kind = _refuse_observation_kind(kind, self._QUANTITY)
        members = _zone_members(self.zone_definitions(model))
        # The cellid the model's own package reads (ledger 455): `(layer, cell)`
        # was written on every grid, which on DIS and DISU names no such cell.
        numbering = cell_numbering(model) if has_cell_numbering(model) else None
        records = [
            (label, kind, cellid if numbering is None else numbering.cellid(*cellid))
            for entries in _zone_labels(members, owner=self._QUANTITY.owner).values()
            for label, cellid in entries
        ]
        return {str(csv_name): records}

    def _attach_flopy_obs(
        self, model, *, pname: str, filename, csv_name: str | None, kind: str, package
    ):
        """The OBS package ``attach_flopy_obs`` documents, for the family's quantity."""

        import flopy

        owner = package if package is not None else _owning_package(model, self._QUANTITY, self.package)
        csv_name = (
            _default_obs_csv_name(filename, self._DEFAULT_CSV) if csv_name is None else str(csv_name)
        )
        continuous = self._to_flopy_obs(model, csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )


class _BoundCellZoneFlowTargets:
    """The shared body of the model-bound zone-flow helpers (``model.targets.drn_flow`` / ``riv_flow``).

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from, and ``get``,
    ``to_long`` and ``summary`` resolve geometry zones against it. Flows are in
    MF6's sign (ledger 350).
    """

    _PLOT_TITLE: ClassVar[str] = "Observed vs simulated zone flow"

    def __init__(self, model, targets: _CellZoneFlowTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The zone definitions with member cells resolved against the bound model."""

        return self.targets.get(model=self.model)

    def to_long(self) -> pd.DataFrame:
        """The long-format zone-flow target table (resolved against the bound model)."""

        return self.targets.to_long(model=self.model)

    def to_wide(self) -> pd.DataFrame:
        """The observed flows in wide format (one column per zone)."""

        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        """A compact summary of the zone target set (resolved against the bound model)."""

        return self.targets.summary(model=self.model)

    def simulated_series(self, model=None) -> pd.DataFrame:
        """The wide simulated-flow series per zone (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to read; zones drawn as geometry are resolved against
            ITS grid. Default ``None`` uses the bound model. Read by the wrapped
            set's ``simulated_series`` -- the package's declared observations when
            every zone cell is declared there, else the model's cell-budget file
            -- one row per stress period at its END, each zone the sum of its
            cells, whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            ``time`` (the zero-based stress period) plus one flow column per
            zone, in model length cubed per time and MF6's sign.
        """

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated zone flow (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against; zones drawn as geometry are
            resolved against ITS grid. Default ``None`` uses the bound model;
            another run compares the same targets against it -- what the plots'
            ``baseline=`` does. Read as the wrapped set's ``compare`` reads it
            (each target's ``time`` a zero-based stress period), whose ``Raises``
            apply.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row -- one per zone and
            period, ``flow_target`` and the residuals ``NaN`` -- when the set
            holds no values): ``name``, ``time``, ``per``, ``flow_target``,
            ``sim_flow``, ``residual`` (simulated minus target), ``abs_residual``,
            and the zone's ``group``, ``weight`` and ``cells``.
        """

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the zone targets (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``group``, ``time``,
            ``per``, ``flow_target``, ``sim_flow``, ``residual``,
            ``abs_residual``, ``weight`` and ``cells`` (absent only when
            :meth:`compare` has none), as the wrapped set's ``residuals``
            reduces them.
        """

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            One row: ``n``, ``mean_error``, ``mae`` and ``rmse`` of simulated
            minus target flow, in model length cubed per time.
        """

        return self.targets.stats(self.model if model is None else model)

    @property
    def plot(self):
        """Plotting helpers (obs-vs-sim, time series) for these zone-flow targets."""

        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="flow_target",
                simulated_column="sim_flow",
                title=self._PLOT_TITLE,
                yaxis_title="Flow",
            )
        return self._plot
