"""DRN flow observation targets (`DrnFlowTargets` + bound view)."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

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
    _compare_without_values,
    _copy_frame,
    _declared_route,
    _default_compare_stats,
    _default_obs_csv_name,
    _in_grid_crs,
    _is_blank,
    _measured_periods,
    _normalize_compare_time,
    _normalize_named_series_values,
    _observation_label,
    _ObservedQuantity,
    _owning_package,
    _placement_note,
    _read_declared_observations,
    _refuse_flows_leaving_the_frame,
    _refuse_keys_off_the_periods,
    _refuse_observation_kind,
    _refuse_reserved_target_names,
    _report_missed_features,
    _wide_by_period,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.period_end import period_ends

logger = get_logger(__name__)

#: What a seepage zone reads on the observation route: a DRN package's DRN
#: observations, one per zone cell (ledger 344). The budget route sums every
#: drain in a zone's cells whichever DRN package holds it.
_DRN_FLOW = _ObservedQuantity(
    owner="DrnFlowTargets",
    package_type="drn",
    obstype="drn",
    feature="cell",
    file_route="the model's cell-budget file",
    file_route_reads_every_package=True,
)


def _zone_layer(value, name) -> int:
    """One zone row's layer as a zero-based int; blank is 0, anything else not a layer raises."""

    if value is None or (not isinstance(value, str) and pd.isna(value)) or str(value).strip() == "":
        return 0
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = float("nan")
    if not np.isfinite(number) or number < 0 or number != int(number):
        raise ValueError(
            f"DrnFlowTargets zone {str(name)!r}: layer {value!r} is not a zero-based layer "
            "number (a whole number, 0 or more)."
        )
    return int(number)


def _with_geometry(frame: pd.DataFrame, *, name_column: str = "name") -> pd.DataFrame:
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
            "DrnFlowTargets locations: a 'geometry' entry must be a shapely geometry "
            f"(a Point, LineString or Polygon), not {listing}. Parse text first -- "
            "shapely.wkt.loads(...) -- or give the zones as a GeoDataFrame or vector file."
        )
    located = frame.copy()
    located["geometry"] = values
    return gpd.GeoDataFrame(located, geometry="geometry")


def _normalize_drn_zone_locations(
    locations,
    *,
    name_column: str = "name",
    group_column: str | None = "group",
    weight_column: str | None = "weight",
) -> pd.DataFrame:
    """Normalize drain-zone definitions from cells or geometry.

    Every row comes back with ``name``, ``group``, ``weight``, ``layer`` (an int,
    0 when blank) and ``cells`` when the input has that column: a list of ints, or
    ``None`` where the row gives none -- a row :func:`_resolve_drn_zone_cells`
    resolves from its geometry. An EMPTY list (``[]``, or ``""`` in a CSV) is an
    explicit zone of no cells, taken as given.
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
        raise ValueError("DrnFlowTargets locations must include a name column.")
    frame = _with_geometry(frame, name_column=name_column)
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
    if "layer" not in frame.columns:
        frame["layer"] = 0
    frame["layer"] = [
        _zone_layer(value, name) for value, name in zip(frame["layer"], frame["name"], strict=True)
    ]
    if "cells" in frame.columns:
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

        frame["cells"] = pd.Series(
            [_coerce_cells(values) for values in frame["cells"]], index=frame.index, dtype=object
        )
    return frame


def _zone_listing(names) -> str:
    """Zone names for a message: the first ten, quoted, and a count of the rest."""

    names = list(dict.fromkeys(str(name) for name in names))
    listing = ", ".join(repr(name) for name in names[:10])
    if len(names) > 10:
        listing += f", and {len(names) - 10} more"
    return listing


def _resolve_drn_zone_cells(locations: pd.DataFrame, model, *, owner=None) -> pd.DataFrame:
    """Attach explicit cell lists to each DRN zone, resolving geometry against the grid.

    A row that gives ``cells`` keeps them, as given; a row that gives none -- no
    ``cells`` column, or a missing entry in one -- is resolved from its geometry,
    so one ``GeoDataFrame`` may mix the two (ledger 342: its geometry rows were
    skipped whenever a ``cells`` column existed, and came back with no cells).

    A zone COVERS the cells it is drawn over, by the shared rule
    (:func:`myflopy.modflow.mf6.grid.geometry.cell_shares`): a seepage polygon
    reaches the cells whose AREA it overlaps, a drain trace the cells it runs
    through, and a spring point the ONE cell holding it. The join here was
    ``predicate="intersects"``, which is true of a shared edge and of a shared
    corner, so on a 4x4 grid of 25 ft cells a zone drawn exactly on cell 5 claimed
    ``[0, 1, 2, 4, 5, 6, 8, 9, 10]``, a point on a face claimed both its cells and
    one on a vertex four. Worst of all, two zones abutting on the column of faces
    at x=50 BOTH claimed ``[1, 2, 5, 6, 9, 10, 13, 14]``, so every drain in that
    column was summed into two different observations and calibrated twice.

    Only cells ACTIVE in the zone's own layer are kept, since MF6 writes neither a
    DRN budget entry nor a DRN observation where there is no cell. Membership is
    deliberately NOT restricted to cells that carry a drain: a zone is geometry,
    resolved with no reference to the DRN package, and a covered cell with no
    drain adds zero to the budget sum.

    A zone given as SEVERAL ROWS is ONE zone, as the grouped join made it: each
    geometry row carries the union of the cells its zone's geometry rows cover IN
    ITS LAYER, and the zone has missed only when it holds no cell at all, which is
    also how the set is counted.

    Misses follow ledger 211. A zone set that lands nowhere raises, naming both
    extents and both CRSs, where an off-grid zone used to come back as ``[]`` and
    its flow as a silent zero; a partial miss warns, by name, once per set rather
    than once per call, because every consumer resolves the zones again.
    """

    frame = locations.copy()
    if "cells" in frame.columns:
        pending = frame["cells"].isna().to_numpy()
    else:
        pending = np.ones(len(frame), dtype=bool)
    if not pending.any():
        return frame
    waiting = frame.loc[pending, "name"]
    if not isinstance(frame, gpd.GeoDataFrame):
        raise ValueError(
            f"DrnFlowTargets zone(s) {_zone_listing(waiting)} give no cells and no geometry. "
            "Give each a 'cells' list, or give the zones as a GeoDataFrame, a vector file, "
            "or record dicts carrying a 'geometry'."
        )
    empty = frame.geometry.isna().to_numpy() | frame.geometry.is_empty.to_numpy()
    if (pending & empty).any():
        raise ValueError(
            f"DrnFlowTargets zone(s) {_zone_listing(frame.loc[pending & empty, 'name'])} give "
            "neither cells nor a geometry, so there is nothing to resolve them from."
        )
    if model is None:
        raise ValueError(
            f"DrnFlowTargets zone(s) {_zone_listing(waiting)} are drawn as geometry, which is "
            "resolved against a model's grid; pass the model (it is None here), or give "
            "their 'cells'."
        )
    polys = model_cell_polygons(model)
    domain = model_domain(model, len(polys))
    located = _in_grid_crs(frame, model)
    covered: dict[tuple[str, int], list[int]] = {}
    reasons: dict[str, str] = {}
    for name, layer, geometry in zip(
        located.loc[pending, "name"], located.loc[pending, "layer"],
        located.geometry[pending], strict=True,
    ):
        key = (str(name), int(layer))
        hits = cell_shares(
            polys,
            geometry,
            # Eligibility only. ``cell_shares`` sorts whatever this returns, so
            # the lowest-numbered usable cell wins a tie without a second copy of
            # that rule here.
            usable=lambda cells, layer=int(layer): [
                cell for cell in cells if domain_active(domain, layer, cell)
            ],
        )
        if not hits:
            reason = miss_reason(polys, geometry, target=f"active cell in layer {int(layer)}")
            reasons.setdefault(key[0], reason)
        covered[key] = sorted(set(covered.get(key, ())) | {cell for cell, _share in hits})
    resolved = located.copy()
    # A list per row, not the shared one: two rows of one zone must not alias.
    resolved["cells"] = pd.Series(
        [
            list(covered[(str(name), int(layer))]) if waits else list(cells)
            for name, layer, cells, waits in zip(
                resolved["name"], resolved["layer"],
                resolved["cells"] if "cells" in resolved.columns else [None] * len(resolved),
                pending, strict=True,
            )
        ],
        index=resolved.index,
        dtype=object,
    )
    # Counted by NAME, not by row. A multipart zone has one row per part, and
    # accumulating the misses inside the loop said "1 of 2 DRN zones land on no
    # active cell" of a single zone given as two rows -- naming the zone that had
    # in fact landed, and which the frame still carried.
    members = _zone_members(resolved)
    missed = [f"{name!r} ({reasons[name]})" for name, cells in members.items()
              if not cells and name in reasons]
    if missed:
        _report_missed_features(
            owner,
            missed,
            total=len(members),
            subject="DRN zone",
            nothing_left="sum",
            where=lambda: _placement_note(frame, model, label="zones"),
            log=logger,
        )
    return resolved


def _zone_members(zones: pd.DataFrame) -> dict[str, set[tuple[int, int]]]:
    """Each zone's member cells by NAME, as ``(layer, cell)``: several rows are ONE zone.

    A cell is a ``(layer, cell)`` pair because that is what MF6 names: a DRN
    observation's cellid, and the model-wide node a cell-budget record carries
    (``layer * ncpl + cell`` on a DISV). A zone keyed by its in-layer numbers
    alone read, below layer 0, whichever layer-0 drain shared a number
    (ledger 342). A zone given as several rows -- explicit cells, or geometry, in
    one layer or several -- is united here, so its flow is summed once over all
    of them rather than once per row.

    Parameters
    ----------
    zones : pandas.DataFrame
        Zones with ``name``, ``layer`` and ``cells``
        (``DrnFlowTargets.zone_definitions``); a row whose ``cells`` is ``None``
        (unresolved geometry) adds nothing.

    Returns
    -------
    dict of str to set of tuple of (int, int)
        Zone name to its ``(layer, cell)`` pairs, zero-based, in the order the
        zones first appear.
    """

    members: dict[str, set[tuple[int, int]]] = {}
    for row in zones.itertuples(index=False):
        layer = int(getattr(row, "layer", 0))
        members.setdefault(str(row.name), set()).update(
            (layer, int(value)) for value in (getattr(row, "cells", None) or ())
        )
    return members


def _refuse_cells_off_the_grid(members: dict[str, set[tuple[int, int]]], node_to_lni) -> None:
    """Refuse a zone ``(layer, cell)`` the model's grid does not have, naming the zones.

    An explicit ``cells`` entry is taken as given, and on the budget-file route a
    pair the grid does not have -- a cell number past its layer's count, a
    negative one, a layer below the bottom -- matched no budget row, so the zone
    summed the cells that DID exist and said nothing: a zone of cells ``[20, 99]``
    on 36-cell layers read cell 20's flow as the zone's (ledger 343). MF6 refuses
    such a cellid on the observation route ("outside of the grid"), so the two
    routes now agree. A cell the grid has but whose ``idomain`` removes is NOT
    refused: it holds no drain and adds zero, as any covered cell without one does.

    Parameters
    ----------
    members : dict of str to set of tuple of (int, int)
        :func:`_zone_members`.
    node_to_lni : dict of int to tuple of (int, int)
        The model's ``node_to_lni``: every cell of its grid as ``(layer, cell)``.

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
        "DrnFlowTargets: zone(s) name (layer, cell) pairs this model's grid does not "
        f"have, zero-based: {listing}. The grid has {shape}; correct the zone's "
        "`cells` or its `layer`."
    )


def _zone_labels(members: dict[str, set[tuple[int, int]]]) -> dict[str, list[tuple[str, tuple[int, int]]]]:
    """Each zone's MF6 observation labels, one per member cell, and the cellid each names.

    ``<name>_c<cell>`` for a zone in one layer -- the label :meth:`to_flopy_obs`
    has always written -- and ``<name>_lay<layer>_c<cell>`` for a zone reaching
    several, whose cell numbers can repeat between layers (the head targets'
    ``_lay`` convention). Every character that is not a letter or digit becomes
    ``_``. Two zones whose labels coincide once written -- ``"a b"`` and ``"a_b"``
    -- are refused: MF6 accepts the duplicate and the reader could then only take
    one of the two columns for both.

    Parameters
    ----------
    members : dict of str to set of tuple of (int, int)
        :func:`_zone_members`.

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
            "DrnFlowTargets: two zone cells would carry one MF6 observation label once "
            f"written (letters and digits only): {listing}. Rename one of the zones."
        )
    return labels


@dataclass
class DrnFlowTargets:
    """Observed seepage/discharge targets aggregated from DRN budget cells.

    Sums the simulated DRN flow over a *group* of drain cells to form one flow
    series per named seepage zone, then aligns it with observed discharge -- the
    pattern for calibrating spring/seep or drain-discharge measurements where the
    field datum is a lumped flux over many cells. Unlike the SFR/LAK targets,
    ``locations`` maps each ``name`` to a set of drain cells (via geometry or an
    explicit cell list), and the comparison aggregates the DRN cell-by-cell
    budget. Pass an instance to ``cal.observe(...)`` / ``cal.forecast(...)``, or
    call ``.compare(model)`` / ``.stats(model)`` for standalone review.

    **A period's flow is its END-of-period flow** (ledger 336): the zone's cells
    ADD within the time step that ends the period -- that is what a lumped seepage
    measurement is -- and the steps a period saved do not; a period whose final
    step output control did not save is refused rather than read early. On the
    observation route the zone is its cells' observations at the time step that
    ends each period, keyed by the period (ledger 346).

    **A flow is in MF6's own sign, on both routes**: MF6 writes water a drain
    takes out of the aquifer as NEGATIVE, so ``sim_flow`` is zero or negative and
    a measured spring or drain discharge ``Q`` is entered as ``-Q``. A target
    value above zero can never be matched -- a drain cannot add water -- and is
    REFUSED when the targets are built, naming ``-Q`` (ledger 350); zero (a dry
    spring) is accepted. A positive discharge was accepted before and scored a
    residual of ``-2Q`` at a perfect fit, on either route.

    **A zone's cell is a ``(layer, cell)`` pair, on every route** (ledger 342): the
    budget file names a drain by its model-wide node, and matching that against a
    zone's in-layer numbers read a zone in layer 1 as whichever layer-0 drain
    shared a number. A zone is summed ONCE over the union of its rows' cells, and
    written to MF6 as one observation per cell.

    **The observation route reads what the DRN package DECLARES** (ledger 344):
    a zone reads MF6's ``DRN`` observations only when every one of its cells'
    labels (as :meth:`to_flopy_obs` writes them) is declared in the DRN
    package's own OBS as a ``DRN`` observation OF THAT CELL, from the file that
    OBS names; a label declared as another quantity only, for another cell or
    twice raises, and so does a set of which some zone cells are declared and
    some are not. It used to read the first workspace CSV whose header held the
    labels -- a stale copy sorting first froze every PEST realization at it.

    Parameters
    ----------
    locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, list or tuple
        The seepage zones, one name each; rows sharing a name are ONE zone, its
        cells united across the rows and their layers. Each row names its cells
        one of two ways. EXPLICIT: a ``cells`` entry holding a list of ints, one
        int or a comma-separated string (``"317, 342"``) -- ZERO-based cell
        numbers within the row's ``layer``, taken as given (an empty list, or
        ``""``, is a zone of no cells), except that reading the budget file
        refuses a ``(layer, cell)`` the model's grid does not have (ledger 343;
        MF6 refuses it on the observation route); the ``{name: cells}`` dict and a
        ``Series`` of cells indexed by name are shorthands for it. GEOMETRY: a
        row with no ``cells`` entry (no such column, or a missing value in it)
        and a geometry -- :meth:`zone_definitions` resolves it to the cells it
        COVERS that are active in its layer. One table may mix the two. Geometry
        comes from a ``GeoDataFrame``, from any other value (a vector-file path)
        read with ``geopandas.read_file``, or from a ``geometry`` column of
        shapely geometries in a dict, a list of record dicts or a ``DataFrame``,
        which becomes a ``GeoDataFrame`` in the GRID's coordinates (it carries no
        CRS; anything but a geometry there raises ``ValueError``). A dict with a
        ``name``, ``cells`` or ``geometry`` key is read as columns, and a list or
        tuple as record dicts. Kept: ``name_column``, ``cells``,
        ``group_column``, ``weight_column``, ``layer`` and the geometry; every
        other column is DROPPED, a ``node`` or ``cell`` column included.
        ``layer`` is zero-based, default 0, a blank read as 0: the layer a row's
        cells are in, which every route reads -- MF6's cell-budget node is
        ``layer * ncpl + cell`` on a DISV, and the budget route matches the pair
        (ledger 342). ``ValueError`` when ``name_column`` is missing, or a
        ``layer`` is not a whole number 0 or more.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed discharge through time, in model length cubed per model time
        unit, in MF6's own sign: water leaving the aquifer through a drain is
        NEGATIVE, so enter a measured discharge ``Q`` as ``-Q``. :meth:`compare`
        subtracts it, as given, from the simulated flow, which carries MF6's sign
        on both routes. A value above zero is refused (``ValueError`` naming the
        zone, the time and ``-Q``, ledger 350); zero is accepted. A path is read with
        ``pandas.read_csv``. LONG form has ``name`` and ``value_column`` columns
        (headers matched ignoring case); WIDE form has one column per zone name
        beside ``time_column``; a table of only ``time_column`` and
        ``value_column``, or a ``Series`` indexed by time, is one series and
        needs exactly one zone (``ValueError`` otherwise). Names are NOT checked
        against ``locations`` and are joined to them exactly (case counts): a
        row naming no zone gets no cells and no simulated flow, and the PEST
        build refuses it. A blank value (missing, or empty text) is unmeasured; any other
        entry that is not a finite number raises ``ValueError`` naming the target and
        its column (ledger 340). Default ``None`` makes a
        capture template: observed, its rows register at weight 0; forecast,
        every row is a forecast.
    name_column : str, optional
        Column of ``locations`` holding each zone's name, matched exactly and
        renamed ``name``. Default ``"name"``. For the ``{name: cells}`` dict and
        ``Series`` shorthands it is only the column the names are put in.
    group_column : str or None, optional
        Column of ``locations`` holding a group label, renamed ``group``.
        Default ``"group"``; ``None`` or an absent column leaves ``group``
        ``pandas.NA``. Carried into :meth:`get`, :meth:`compare`,
        :meth:`residuals` and the zone table a PEST build writes. It is the
        PEST observation group of the zone's MEASURED observations (ledger
        358), made a PEST name and refused as
        :class:`~myflopy.modflow.mf6.observations.heads.HeadTargets`'
        ``group_column`` says; a missing or blank label puts them in the set's
        prefix (``drn_flow``). A zone given as several rows carrying DIFFERENT
        labels is refused by the build -- one zone is one observation per
        period -- and rows with no label take the others', in the build and in
        :meth:`compare` alike (ledger 359).
    weight_column : str or None, optional
        Column of ``locations`` holding each zone's observation weight, renamed
        ``weight``. Default ``"weight"``; ``None`` or an absent column leaves it
        ``NaN``. A PEST build sets a zone's non-``NaN`` weight on each
        observation of that zone that has a target value, and weight 1 where it
        is ``NaN``; a zone given as rows with DIFFERENT weights is refused, and
        every period without a target value gets weight 0 (ledger 340).
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` and the PEST build read it as a
        zero-based stress period, by number (``1.0`` is period 1; a fraction
        raises), on both routes (the observation route since ledger 346, where a
        key naming no period is refused -- MF6's simulated time was its key).
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed flow, matched
        ignoring case. Default ``"flow"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.
    package : str or None, optional
        The NAME of the DRN package whose observations the zones read (FloPy's
        ``package_name``, any case), needed when the model carries several DRN
        packages and the zones' cells are declared in their observations: a DRN
        observation measures ONE package's drain in a cell. Default ``None``:
        the model's one DRN package; with several, reading raises
        ``ValueError`` naming them when a zone cell's label is declared in any
        of them, and otherwise the budget route sums every package's drains in
        each zone, as before. A named package that declares none of the zone
        cells raises too -- the budget route cannot read one package's drains
        apart from the others' -- as does a ``package`` the model does not
        have. :meth:`attach_flopy_obs` attaches to it.

    Examples
    --------
    >>> DrnFlowTargets(locations=seep_zones_gdf,
    ...                values=observed_seepage_df)       # -> cal.observe(...)

    See Also
    --------
    SfrFlowTargets, HeadTargets, LakeStageTargets
    """

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None = None
    name_column: str = "name"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "flow"
    times: list | tuple | pd.Index | pd.Series | None = None
    package: str | None = None

    def __post_init__(self):
        """Normalize the zone definitions and observed discharge into internal tables."""

        self._locations = _normalize_drn_zone_locations(
            self.locations,
            name_column=self.name_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
        )
        _refuse_reserved_target_names(
            self._locations["name"], owner=type(self).__name__,
            columns=("flow_target", "sim_flow"),
        )
        self._values = _normalize_named_series_values(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            output_column="flow_target",
            times=self.times,
            owner=type(self).__name__,
        )
        if not self._values.empty and "name" not in self._values.columns:
            if len(self._locations["name"].unique()) != 1:
                raise ValueError("DrnFlowTargets values without explicit names require exactly one zone.")
            self._values["name"] = str(self._locations["name"].iloc[0])
        _refuse_flows_leaving_the_frame(
            self._values,
            "flow_target",
            owner=type(self).__name__,
            quantity="MF6's DRN flow",
            leaving="the aquifer through a drain",
            measured="a measured discharge",
        )

    @property
    def locations_gdf(self):
        """A copy of the normalized zone-definition table (geometry preserved)."""

        return self._locations.copy()

    def zone_definitions(self, model) -> pd.DataFrame:
        """Resolve each zone's member cells against ``model`` (covering geometry if needed).

        A row given ``cells`` keeps them as given; a row given none is resolved from
        its geometry: it reaches the cells it COVERS, active in its own layer -- the
        rule ``.gpkg`` boundaries and head targets use -- and comes back in the
        GRID's CRS, whatever the zone file's own was. See
        :func:`_resolve_drn_zone_cells` for what that rule replaced and for what a
        zone landing on nothing does.

        Parameters
        ----------
        model : SimulationBase or None
            The model whose grid (``model.vor``) and active domain geometry rows
            are resolved against. Not read when every row gives its ``cells``,
            which are returned as given; ``None`` then works.

        Returns
        -------
        pandas.DataFrame or geopandas.GeoDataFrame
            One row per ``locations`` row: ``name``, ``cells`` (zero-based cell
            numbers in the row's layer: an explicit row's own; a geometry row's,
            the union over every geometry row of that name IN THAT LAYER),
            ``group``, ``weight``, ``layer``, and, where the zones carry
            geometry, the geometry, converted to the grid's CRS when both CRSs
            are known.

        Raises
        ------
        ValueError
            When a row gives no cells and has no geometry to resolve them from,
            or ``model`` is ``None`` while one does; or when no zone lands on a
            cell active in its layer, naming the zones (the first ten) and both
            extents and CRSs (ledger 211).

        Warns
        -----
        UserWarning
            When SOME zones land on no active cell, naming them, once per missed
            set on this object -- through ``warnings.warn`` as well as the
            logger, because a log line alone reaches nobody by default (ledgers
            288 and 342).
        """

        return _resolve_drn_zone_cells(self._locations, model, owner=self)

    def get(self, model=None) -> pd.DataFrame:
        """Return the zone definitions (cells resolved when ``model`` is given), merged with values.

        Parameters
        ----------
        model : SimulationBase or None, optional
            When given, zones are resolved through :meth:`zone_definitions`,
            whose ``Raises`` apply. Default ``None``: the zones as normalized at
            construction -- a geometry row's ``cells`` is ``None``, unresolved,
            and its geometry in the zone file's own CRS.

        Returns
        -------
        pandas.DataFrame
            The zone table alone when there are no ``values``; otherwise the
            observed values -- ``time``, ``name``, ``flow_target`` -- joined on
            ``name`` to the zone table's columns other than geometry: one row
            per value and zone ROW, so a zone given as several rows repeats each
            of its values, and a value naming no zone gets one row of ``NaN``.
        """

        definitions = self._locations.copy()
        if model is not None:
            definitions = self.zone_definitions(model)
        if self._values.empty:
            return definitions
        return _copy_frame(self._values.merge(definitions.drop(columns="geometry", errors="ignore"), on="name", how="left"))

    def to_long(self, model=None) -> pd.DataFrame:
        """Return the long-format DRN-flow target table (the same table as :meth:`get`).

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
        """Return the observed discharge in wide format with one column per zone.

        One row per stress period a ``time`` names -- ``time`` is the zero-based
        period as ``int``, read as :meth:`compare` and the PEST build read it, so ``"1"``
        and ``1`` are one row (they were two, ledger 422) -- and one column per target
        name holding its observed flow, in model flow units and MF6's sign (ledger 350); a set given no values gives ``time`` and one
        empty column per target. Two rows HOLDING a value for one target in one
        stress period are REFUSED -- names compared stripped and lower-cased --
        by the check :meth:`compare` and the PEST build read (ledger 340): this
        used to keep the first in silence (``pivot_table(aggfunc="first")``), so
        the wide table and the build disagreed about the same set (ledger
        371). A blank row is unmeasured and never collides; a target whose
        every row is blank keeps its column, all blank, as :meth:`compare` keeps
        its rows (ledger 422).

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

        if self._values.empty:
            return _wide_by_period(
                self._values, value_column="flow_target", caller="DrnFlowTargets.to_wide()",
                names=self._locations["name"],
            )
        _measured_periods(self._values, value_column="flow_target", caller="DrnFlowTargets.to_wide()")
        # After the check a (target, time) cell holds at most ONE value, so
        # "first" only prefers it over blank rows.
        return _wide_by_period(self._values, value_column="flow_target", caller="DrnFlowTargets.to_wide()")

    def summary(self, model=None) -> pd.DataFrame:
        """Return a compact summary (zone count, resolved cell count, row count).

        Parameters
        ----------
        model : SimulationBase or None, optional
            When given, zones are resolved through :meth:`zone_definitions`
            first, whose ``Raises`` apply. Default ``None``: geometry rows are
            unresolved and count no cells.

        Returns
        -------
        pandas.DataFrame
            One row: ``n_zones`` (distinct names), ``n_cells`` (distinct
            ``(layer, cell)`` pairs over the whole set -- a cell is counted once
            however many rows or zones name it; it counted every row's list, so a
            geometry zone given as several rows, each carrying the zone's union,
            was counted once per row) and ``n_rows`` (observed values).
        """

        definitions = self.zone_definitions(model) if model is not None else self._locations.copy()
        zone_count = int(definitions["name"].nunique())
        cells = set().union(*_zone_members(definitions).values()) if len(definitions) else set()
        return pd.DataFrame([{"n_zones": zone_count, "n_cells": len(cells), "n_rows": int(len(self._values))}])

    def _drn_budget_frame(self, model) -> pd.DataFrame:
        """The DRN cell-by-cell budget as a tidy frame: ``per``, ``node``, ``layer``, ``cell``, ``q``.

        ``node`` is MF6's model-wide USER node, zero-based -- ``layer * ncpl +
        cell`` on a DISV, whatever the ``idomain`` removes (MF6 writes user, not
        reduced, node numbers; pinned by ``tests/test_drn_zone_layers.py``) --
        and ``layer``/``cell`` are that node through ``model.node_to_lni``, the
        one conversion the package budget tables use.

        Raises
        ------
        ValueError
            When the budget carries neither ``per`` nor ``kstpkper``, no ``node``,
            or a node that is not a cell of the model's grid -- the file and the
            model do not line up.
        """

        budget = model.bud("drn").df.reset_index().copy()
        if "per" not in budget.columns:
            if "kstpkper" not in budget.columns:
                raise ValueError("DRN budget dataframe requires 'per' or 'kstpkper' columns.")
            budget["per"] = budget["kstpkper"].apply(lambda item: int(item[1]))
        if "node" not in budget.columns:
            raise ValueError("DRN budget dataframe requires a zero-based 'node' column.")
        budget["node"] = pd.to_numeric(budget["node"], errors="coerce").astype(int)
        budget["q"] = pd.to_numeric(budget["q"], errors="coerce")
        lookup = model.node_to_lni
        nodes = pd.unique(budget["node"])
        stray = sorted(int(node) for node in nodes if int(node) not in lookup)
        if stray:
            raise ValueError(
                f"The DRN cell-budget record names node(s) {stray[:5]} (zero-based), which "
                f"are not cells of this model's grid ({len(lookup)} nodes): the budget "
                "file and the model do not line up."
            )
        pairs = {int(node): lookup[int(node)] for node in nodes}
        budget["layer"] = [int(pairs[int(node)][0]) for node in budget["node"]]
        budget["cell"] = [int(pairs[int(node)][1]) for node in budget["node"]]
        return budget

    def _zone_flow_table(self, model, zones: pd.DataFrame) -> pd.DataFrame:
        """Each zone's DRN flow at each stress period's END, from the model's budget file.

        A zone's flow is the sum of its cells' DRN budget rows in ONE time step
        (``kind="extensive"`` -- a seepage zone IS that sum), read at the record
        that ends each period on the budget file's own clock (``SAVE BUDGET``,
        ledger 301), through the one period-end rule every target reads (ledger
        336). It used to add every step a period saved as well: a one-cell zone
        on the canonical model at ``nstp=3`` under ``SAVE BUDGET ALL`` read
        -292.51 where each step held -97.50.

        A budget row belongs to a zone by its ``(layer, cell)``, not its in-layer
        cell number: the record's node is model-wide, so matching it against a
        zone's cell numbers read, for a zone in layer 1, the LAYER-0 drain
        sharing a number (ledger 342: -6.77 for a zone whose drains held
        -5975.73).

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package and a cell-budget file.
        zones : pandas.DataFrame
            The resolved zones (:meth:`zone_definitions`): ``name``, ``layer``
            and ``cells``.

        Returns
        -------
        pandas.DataFrame
            ``per``, ``name`` and ``sim_flow`` (MF6's sign: a drain's outflow is
            negative), one row per zone with cells per period the file saved.

        Raises
        ------
        ValueError
            When a period the budget file saved lacks its final step -- naming
            the period, the steps saved and ``SAVE BUDGET LAST``; when the
            budget names a node the grid does not have; or when a zone names a
            ``(layer, cell)`` the grid does not have (ledger 343: it matched no
            row, and the zone summed the rest silently).
        """

        budget = self._drn_budget_frame(model)
        members = _zone_members(zones)
        _refuse_cells_off_the_grid(members, model.node_to_lni)
        names = list(members)
        located = pd.MultiIndex.from_arrays([budget["layer"], budget["cell"]])
        pieces = [
            budget.loc[located.isin(sorted(pairs))].assign(zone=index)
            for index, pairs in enumerate(members.values())
            if pairs
        ]
        if not pieces:
            return pd.DataFrame(columns=["per", "name", "sim_flow"])
        rows = pd.concat(pieces, ignore_index=True)
        ends = None
        if "kstpkper" in rows.columns:
            ends = period_ends(
                model,
                reader=lambda: model._get_budget_reader(),
                file="the model's cell-budget file",
                saves="BUDGET",
                fallback=budget["kstpkper"],
            )
        table = _aggregate_named_series(
            rows, id_column="zone", value_column="q", kind="extensive", ends=ends
        )
        table["name"] = [names[int(zone)] for zone in table["zone"]]
        return table.rename(columns={"q": "sim_flow"}).loc[:, ["per", "name", "sim_flow"]]

    def _zone_details(self, zones: pd.DataFrame) -> pd.DataFrame:
        """One row per zone NAME: its ``group``, ``weight`` and ``cells`` as ``(layer, cell)``.

        What :meth:`compare` joins to each simulated row, on both routes. The
        observation route joined the zone table's ROWS, so a zone given as two
        rows doubled every comparison row and ``stats`` counted each residual
        twice. A zone given as several rows takes its ``group`` and ``weight``
        from the first row that carries one, as the PEST build does
        (``pest.observations._zone_groups`` / ``_zone_weights``, which refuse
        rows that DISAGREE): it took the first row's, so a zone whose first row
        was unlabelled read ``<NA>`` here and its label in PEST (ledger 359).
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

        Each zone cell's label is the one :meth:`to_flopy_obs` writes, from the one
        :func:`_zone_labels`, and it must be declared in the DRN package's own OBS
        as a ``DRN`` observation OF THAT CELL (ledger 344); a zone reads the sum of
        its cells' columns in each time step. A zone given as several rows is ONE
        zone (:func:`_zone_members`): its labels are those of the union of its
        rows' cells, each once -- it was keyed per row by assignment, so the LAST
        row won and a zone of rows ``[20]`` and ``[21, 27]`` summed cells 21 and 27
        only on this route (ledger 336). One row per stress period, at its end
        (ledger 346).

        Parameters
        ----------
        model : SimulationBase
            The run model the zones are resolved against.

        Returns
        -------
        tuple of (pandas.DataFrame, pandas.DataFrame or None)
            The resolved zones (:meth:`zone_definitions`), and ``time`` (the
            zero-based stress period) plus one column per zone -- a zone of no
            cells reads 0.0 -- or ``None`` when no zone cell is declared.

        Raises
        ------
        ValueError
            As :meth:`zone_definitions` and :func:`_zone_labels`; and as the
            observation route refuses (:meth:`simulated_series`).
        FileNotFoundError
            When a declared observation file is not in the workspace.
        """

        zones = self.zone_definitions(model)
        labels = _zone_labels(_zone_members(zones))
        wanted = {
            f"{zone_name}__{index}": (zone_name, label, cellid)
            for zone_name, entries in labels.items()
            for index, (label, cellid) in enumerate(entries)
        }
        chosen = _declared_route(model, _DRN_FLOW, wanted, package=self.package)
        if chosen is None:
            return zones, None
        observed = _read_declared_observations(model, chosen, owner=_DRN_FLOW.owner)
        result = pd.DataFrame({"time": observed["time"]})
        for zone_name, entries in labels.items():
            keys = [f"{zone_name}__{index}" for index in range(len(entries))]
            result[zone_name] = observed.loc[:, keys].sum(axis=1) if keys else 0.0
        return zones, result

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
            targets = self._values.copy()
            targets["per"] = _measured_periods(
                targets, value_column="flow_target", caller="DrnFlowTargets.compare(...)"
            )
            _refuse_keys_off_the_periods(
                targets, wide["time"], model=model, caller="DrnFlowTargets.compare(...)"
            )
            compare = targets.merge(simulated, on=["name", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-discharge table: a ``time`` column and one column per zone.

        Read from the DRN package's own ``DRN`` observations when EVERY zone
        cell's label (as :meth:`to_flopy_obs` writes it) is declared there, for
        that cell -- only the files that OBS names are read, each record from
        its own column; each row is one stress period, keyed by the period and
        read at the time step that ENDS it (ledger 346 -- it was every time step,
        keyed by MF6's simulated time), the zone the sum of its cells in that
        step -- else, when NO zone cell is declared, from the model's cell-budget
        file (ledger 344; see the class notes for what is refused).

        **On the budget-file route a period holds ONE flow per zone: the flow at
        the period's END** (ledger 336), summed over the zone's cells in that
        step -- the record that ends the period on the budget file's own clock,
        the rule head and stage targets read too. The route used to ADD every
        time step a period saved as well (3.0x at ``nstp=3`` under ``SAVE BUDGET
        ALL``). A period whose final step was not saved (``FIRST``, ``STEPS``
        short of it) is refused, because its end-of-period flow does not exist in
        the file.

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns; zones drawn as geometry are resolved against
            its grid. The observation route reads the DRN package's declared OBS
            and the files it names in the model's workspace.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per zone, in MF6's sign: zero or
            NEGATIVE, the water the zone's drains take out of the aquifer. On the
            budget-file route a zone with no member cell gets no column.

        Raises
        ------
        FileNotFoundError
            On the observation route, when a file the OBS declares is not in the
            workspace.
        ValueError
            When a zone cell's label is declared in the DRN's OBS as another
            quantity only, for another cell or twice as a ``DRN``; when some zone
            cells of the set are declared and others are not; when the model
            carries several DRN packages, a zone cell is declared in them and
            ``package`` does not name one, or ``package`` names one that declares
            none of them; or when a declared file's columns do not line up with
            its records (ledger 344), or it does not end every stress period (a
            run that stopped, a truncated copy -- ledger 346). On the budget-file
            route, when a period
            the budget file saved lacks its final step -- naming the period, the
            steps saved and the ``SAVE BUDGET LAST`` setting that fixes it, or
            when the file names a node
            the model's grid does not have, or a zone names a ``(layer, cell)``
            it does not have (ledger 343); or as :meth:`zone_definitions` (no
            zone landing on an active cell, ledger 211) and :meth:`to_flopy_obs`
            (two zones sharing an observation label).

        Examples
        --------
        >>> seeps = mf.DrnFlowTargets(locations={"north": [317, 342]})
        >>> series = seeps.simulated_series(model)   # one row per stress period
        """

        zones, observed = self._observed(model)
        if observed is not None:
            return observed

        frame = self._zone_flow_table(model, zones)
        if frame.empty:
            return pd.DataFrame(columns=["time", *zones["name"].astype(str).tolist()])
        wide = frame.pivot_table(index="per", columns="name", values="sim_flow", aggfunc="first").reset_index()
        wide.columns.name = None
        return wide.rename(columns={"per": "time"})

    def compare(self, model) -> pd.DataFrame:
        """Compare observed discharge against simulated per-zone DRN flow (with residuals).

        ``sim_flow`` comes from the same two routes as :meth:`simulated_series`:
        on the budget-file route, one flow per zone per stress period, the flow at
        the period's END (ledger 336). Targets are matched on ``per``, read from
        their ``time`` as a zero-based integer -- the stress period, on both routes
        (ledger 346; the observation route was keyed by MF6's simulated time,
        where a target keyed by period matched no row). On the observation route
        a key naming no period of the model is refused, naming the period each
        simulated-time key ends.

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row -- one per zone and
            period, on both routes -- when the targets hold no values): ``name``, ``time``, ``per``,
            ``flow_target``, ``sim_flow``, ``residual`` (simulated minus target),
            ``abs_residual``, and the zone's ``group``, ``weight`` and ``cells``
            -- ONE set per zone however many rows define it, its cells the
            sorted ``(layer, cell)`` pairs MF6 reads. ``flow_target`` and
            ``sim_flow`` are both in MF6's sign -- zero or negative, water
            leaving the aquifer (ledger 350) -- so a residual above zero means
            the model drains LESS than was measured. With no values
            ``flow_target`` and the residuals are ``NaN``, and ``time`` is the
            period, equal to ``per`` -- on the budget-file route since ledger 342
            (it had no ``per`` column, so :meth:`residuals` raised ``KeyError``),
            on the observation route since ledger 346.

        Raises
        ------
        FileNotFoundError
            As :meth:`simulated_series`: a declared observation file is missing.
        ValueError
            As :meth:`simulated_series`: a zone cell's label declared as
            something else, a set declared in part, several DRN packages and no
            ``package``, or a file that does not line up or does not end every
            period; on the observation route, when a target's key is no stress
            period of the model (ledger 346); on the budget-file
            route, when a period the budget file saved lacks
            its final step or a zone names a ``(layer, cell)`` the model's grid
            does not have; or when a target's
            ``time`` is not a whole number, or two rows hold a value for one
            target in one period (the period rule PEST's build reads too, ledger
            340).
        """

        observed = self._compare_from_obs(model)
        if observed is not None:
            return observed

        zones = self.zone_definitions(model)
        simulated = self._zone_flow_table(model, zones)
        simulated = simulated.merge(self._zone_details(zones), on="name", how="left")
        if self._values.empty:
            compare = _compare_without_values(simulated, target_column="flow_target")
        else:
            targets = self._values.copy()
            targets["per"] = _measured_periods(
                targets, value_column="flow_target", caller="DrnFlowTargets.compare(...)"
            )
            compare = targets.merge(simulated, on=["name", "per"], how="left")
        compare = _normalize_compare_time(compare)
        compare["residual"] = compare["sim_flow"] - compare["flow_target"]
        compare["abs_residual"] = compare["residual"].abs()
        return _copy_frame(compare)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table (target/sim/residual, with group/weight) for DRN zones.

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``group``, ``time``,
            ``per``, ``flow_target``, ``sim_flow``, ``residual``,
            ``abs_residual``, ``weight`` and ``cells`` (each zone's
            ``(layer, cell)`` pairs; absent only when :meth:`compare` has none).
            A set given no ``values`` -- a capture template -- gets one row per
            simulated row, its targets and residuals ``NaN``: it raised
            ``KeyError`` naming ``per`` on the budget-file route (ledger 342).
        """

        frame = self.compare(model)
        columns = ["name", "group", "time", "per", "flow_target", "sim_flow", "residual", "abs_residual", "weight"]
        if "cells" in frame.columns:
            columns.append("cells")
        return frame.loc[:, columns]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the DRN-flow target set.

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package; passed to :meth:`compare`, whose
            routes and ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual), ``mean_error``, ``mae`` and
            ``rmse`` of simulated minus target, in model flow units (length
            cubed per time). ``n`` is 0 and the rest ``NaN`` when no row has a
            residual, as with no ``values``.
        """

        return _default_compare_stats(self.compare(model))

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "drn_flow_targets.csv",
        kind: str = "DRN",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        """Return a FloPy ``continuous`` dict with one MF6 DRN observation per member cell.

        Parameters
        ----------
        model : SimulationBase or None
            Passed to :meth:`zone_definitions`, whose ``Raises`` apply: it
            resolves geometry zones. ``None`` works when every zone gives explicit
            ``cells`` -- as when the records go into the ``observations=`` of the
            DRN being declared, before any model exists.
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"drn_flow_targets.csv"``. Any name: :meth:`simulated_series`
            reads the file the DRN package's OBS declares, whatever it is called
            (ledger 344).
        kind : str, optional
            The MF6 observation type written on every record. Must be ``"DRN"``
            (any case), the default -- MF6's drain-flow observation, which
            carries MF6's sign (outflow negative) and is the only type the
            observation route sums as a zone's flow -- and any other raises
            ``ValueError`` (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, kind, (layer, cell)), ...]}``: ONE record per
            cell of each zone -- the ``(layer, cell)`` pairs its rows name,
            united, whether the zone is given as one row or several -- labelled
            ``<name>_c<cell>``, or ``<name>_lay<layer>_c<cell>`` for a zone
            reaching several layers, with every character that is not a letter
            or digit made ``_`` and leading ``_`` stripped; the cellid ZERO-based
            as FloPy takes it. Records followed the ROWS before ledger 342, so
            every cell of a geometry zone given as several rows was written once
            per row under one label. The shape ``ModflowUtlobs(continuous=)`` and
            a DRN package's ``observations=`` take.

        Raises
        ------
        ValueError
            When ``kind`` is not ``"DRN"``; as :meth:`zone_definitions`; or when
            two zones' cells would carry one label once written (``"a b"`` and
            ``"a_b"``), naming them.
        """

        kind = _refuse_observation_kind(kind, _DRN_FLOW)
        members = _zone_members(self.zone_definitions(model))
        records = [
            (label, kind, cellid)
            for entries in _zone_labels(members).values()
            for label, cellid in entries
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "drn_flow_obs",
        filename: str = "drn_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "DRN",
        package=None,
    ):
        """Attach an MF6 utility observation package to DRN for these seepage-zone targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase or None
            The model whose FloPy DRN package receives the OBS package unless
            ``package`` is given: the one the targets name (their ``package``
            field, ``ValueError`` when the model has no DRN package of that
            name), else the model's only DRN package, found by type (several
            raise ``ValueError`` naming them and ``package=``, none
            ``AttributeError``; ledger 352); also passed to :meth:`to_flopy_obs` to
            resolve the zones. With ``package`` given and every zone's ``cells``
            explicit, ``None`` works.
        pname : str, optional
            FloPy package name. Default ``"drn_flow_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory. Default ``"drn_flow_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"drn_flow_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`: ``"DRN"``, the
            default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            the one ``model`` resolves (see ``model``). Default ``None``. Reading
            the observations back needs the targets' own ``package`` field to
            name it when the model carries several DRN packages.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.

        Raises
        ------
        ValueError
            As :meth:`to_flopy_obs`; or when the targets name a package the model
            does not have.
        """

        import flopy

        if package is not None:
            owner = package
        else:
            owner = _owning_package(model, _DRN_FLOW, self.package)
        csv_name = _default_obs_csv_name(filename, "drn_flow_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(model, csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

class BoundDrnFlowTargets:
    """Model-bound DRN seepage-zone helper returned from ``model.targets.drn_flow``."""

    def __init__(self, model, targets: DrnFlowTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The zone definitions with member cells resolved against the bound model."""

        return self.targets.get(model=self.model)

    def to_long(self) -> pd.DataFrame:
        """The long-format DRN-flow target table (resolved against the bound model)."""

        return self.targets.to_long(model=self.model)

    def to_wide(self) -> pd.DataFrame:
        """The observed discharge in wide format (one column per zone)."""

        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        """A compact summary of the seepage-zone target set (resolved against the bound model)."""

        return self.targets.summary(model=self.model)

    def simulated_series(self, model=None) -> pd.DataFrame:
        """The wide simulated-discharge series per zone (defaults to the bound model)."""

        return self.targets.simulated_series(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated DRN discharge (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the seepage-zone targets (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(self, *, csv_name: str = "drn_flow_targets.csv", kind: str = "DRN"):
        """A FloPy ``continuous`` observation dict for the DRN seepage-zone targets."""

        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "drn_flow_obs",
        filename: str = "drn_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "DRN",
        package=None,
    ):
        """Attach an MF6 DRN observation package to the bound model for these zones."""

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    @property
    def plot(self):
        """Plotting helpers (obs-vs-sim, time series) for these DRN seepage targets."""

        if self._plot is None:
            self._plot = BoundNamedSeriesTargetPlots(
                self,
                target_column="flow_target",
                simulated_column="sim_flow",
                title="Observed vs simulated DRN seepage",
                yaxis_title="Flow",
            )
        return self._plot


