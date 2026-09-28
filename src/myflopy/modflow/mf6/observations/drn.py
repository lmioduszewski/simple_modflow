"""DRN flow observation targets (`DrnFlowTargets` + bound view)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.domain import model_domain
from myflopy.modflow.mf6.grid.geometry import cell_shares, domain_active, miss_reason
from myflopy.modflow.mf6.observations._shared import (
    _aggregate_named_series,
    _clean_observation_output_values,
    _copy_frame,
    _default_compare_stats,
    _default_obs_csv_name,
    _find_named_observation_output,
    _in_grid_crs,
    _measured_periods,
    _normalize_compare_time,
    _normalize_named_series_values,
    _observation_label,
    _placement_note,
    _report_missed_features,
)
from myflopy.modflow.mf6.observations.plots import BoundNamedSeriesTargetPlots
from myflopy.modflow.mf6.period_end import period_ends

logger = get_logger(__name__)


def _normalize_drn_zone_locations(
    locations,
    *,
    name_column: str = "name",
    group_column: str | None = "group",
    weight_column: str | None = "weight",
) -> pd.DataFrame:
    """Normalize drain-zone definitions from cells or polygons."""

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
    if "cells" in frame.columns:
        def _coerce_cells(values):
            """Coerce one zone's cell spec (CSV string / sequence / scalar / NaN) to a list of ints."""

            if isinstance(values, str):
                return [int(value) for value in values.split(",") if str(value).strip() != ""]
            if isinstance(values, (list, tuple, np.ndarray, pd.Series)):
                return [int(value) for value in values]
            if pd.isna(values):
                return []
            return [int(values)]

        frame["cells"] = frame["cells"].apply(_coerce_cells)
    return frame


def _resolve_drn_zone_cells(locations: pd.DataFrame, model, *, owner=None) -> pd.DataFrame:
    """Attach explicit cell lists to each DRN zone, resolving geometry against the grid.

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

    A zone given as SEVERAL ROWS is ONE zone, as the grouped join made it: its
    cells are the union over its rows and it has missed only when EVERY row of it
    has, which is also how the set is counted.

    Misses follow ledger 211. A zone set that lands nowhere raises, naming both
    extents and both CRSs, where an off-grid zone used to come back as ``[]`` and
    its flow as a silent zero; a partial miss warns, by name, once per set rather
    than once per call, because every consumer resolves the zones again.
    """

    frame = locations.copy()
    if "cells" in frame.columns and frame["cells"].notna().all():
        return frame
    if not isinstance(frame, gpd.GeoDataFrame):
        raise ValueError("DrnFlowTargets polygon zones require a GeoDataFrame or GIS layer input.")
    polys = model.vor.gdf_vorPolys
    domain = model_domain(model, len(polys))
    located = _in_grid_crs(frame, model)
    layers = located["layer"] if "layer" in located.columns else pd.Series(0, index=located.index)
    covered: dict[str, list[int]] = {}
    reasons: dict[str, str] = {}
    for name, layer, geometry in zip(located["name"], layers, located.geometry, strict=True):
        key = str(name)
        layer = 0 if pd.isna(layer) else int(layer)
        hits = cell_shares(
            polys,
            geometry,
            # Eligibility only. ``cell_shares`` sorts whatever this returns, so
            # the lowest-numbered usable cell wins a tie without a second copy of
            # that rule here.
            usable=lambda cells, layer=layer: [
                cell for cell in cells if domain_active(domain, layer, cell)
            ],
        )
        if not hits:
            reason = miss_reason(polys, geometry, target=f"active cell in layer {layer}")
            reasons.setdefault(key, reason)
        covered[key] = sorted(set(covered.get(key, ())) | {cell for cell, _share in hits})
    # Counted by NAME, not by row. A multipart zone has one row per part, and
    # accumulating the misses inside the loop said "1 of 2 DRN zones land on no
    # active cell" of a single zone given as two rows -- naming the zone that had
    # in fact landed, and which the frame still carried.
    missed = [f"{key!r} ({reasons[key]})" for key, cells in covered.items() if not cells]
    if missed:
        _report_missed_features(
            owner,
            missed,
            total=len(covered),
            subject="DRN zone",
            nothing_left="sum",
            where=lambda: _placement_note(frame, model, label="zones"),
            log=logger,
        )
    resolved = located.drop(columns=["cells"], errors="ignore").copy()
    # A list per row, not the shared one: two rows of one zone must not alias.
    resolved["cells"] = [list(covered[str(name)]) for name in resolved["name"]]
    return resolved


def _zone_members(zones: pd.DataFrame) -> dict[str, set[int]]:
    """Each zone's member cells by NAME: a zone given as several rows is ONE zone.

    ``_resolve_drn_zone_cells`` gives every row of a geometry zone the union of
    its rows' cells already; a zone given as several rows of explicit cells is
    united here, so its flow is summed once over all of them rather than once
    per row.

    Parameters
    ----------
    zones : pandas.DataFrame
        Resolved zones (``DrnFlowTargets.zone_definitions``): ``name`` and
        ``cells``.

    Returns
    -------
    dict of str to set of int
        Zone name to its cells, in the order the zones first appear.
    """

    members: dict[str, set[int]] = {}
    for row in zones.itertuples(index=False):
        members.setdefault(str(row.name), set()).update(
            int(value) for value in getattr(row, "cells", [])
        )
    return members


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
    step output control did not save is refused rather than read early.

    Parameters
    ----------
    locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, list or tuple
        The seepage zones, one name each; rows sharing a name are ONE zone, its
        cells united. A zone's cells come one of two ways. EXPLICIT: a ``cells``
        column holding a list of ints, one int or a comma-separated string
        (``"317, 342"``) -- ZERO-based cell numbers within the zone's layer,
        taken as given; the ``{name: cells}`` dict and a ``Series`` of cells
        indexed by name are shorthands for it. GEOMETRY: a ``GeoDataFrame``, or
        any other value (a vector-file path) read with ``geopandas.read_file``,
        with no ``cells`` column -- :meth:`zone_definitions` resolves each zone
        to the cells it COVERS that are active in its layer. A dict with a
        ``name``, ``cells`` or ``geometry`` key is read as columns, and a list or
        tuple as record dicts. Kept: ``name_column``, ``cells``,
        ``group_column``, ``weight_column``, ``layer`` (zero-based, default 0:
        the layer a geometry zone is resolved in and :meth:`to_flopy_obs`
        writes each record for) and a ``GeoDataFrame``'s
        geometry; every other column is DROPPED, a ``node`` or ``cell`` column
        and a geometry in a plain dict or ``DataFrame`` included. When a
        ``GeoDataFrame`` carries a ``cells`` column, every row is taken as
        explicit and a row whose cells are ``NaN`` gets none. ``ValueError``
        when ``name_column`` is missing.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed discharge through time, in model length cubed per model time
        unit. :meth:`compare` subtracts it, as given, from the simulated flow,
        which carries MF6's sign on both routes -- a drain's outflow is NEGATIVE
        -- so enter a measured discharge negative. A path is read with
        ``pandas.read_csv``. LONG form has ``name`` and ``value_column`` columns
        (headers matched ignoring case); WIDE form has one column per zone name
        beside ``time_column``; a table of only ``time_column`` and
        ``value_column``, or a ``Series`` indexed by time, is one series and
        needs exactly one zone (``ValueError`` otherwise). Names are NOT checked
        against ``locations`` and are joined to them exactly (case counts): a
        row naming no zone gets no cells and no simulated flow, and the PEST
        build refuses it. A blank value (missing, or empty text) is unmeasured; any other
        entry that is not a number raises ``ValueError`` naming the target and
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
        :meth:`residuals` and the zone table a PEST build writes.
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
        raises), matched on the budget-file route against the period.
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed flow, matched
        ignoring case. Default ``"flow"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.

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

    def __post_init__(self):
        """Normalize the zone definitions and observed discharge into internal tables."""

        self._locations = _normalize_drn_zone_locations(
            self.locations,
            name_column=self.name_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
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

    @property
    def locations_gdf(self):
        """A copy of the normalized zone-definition table (geometry preserved)."""

        return self._locations.copy()

    def zone_definitions(self, model) -> pd.DataFrame:
        """Resolve each zone's member cells against ``model`` (covering geometry if needed).

        A zone reaches the cells it COVERS, active in its own layer -- the rule
        ``.gpkg`` boundaries and head targets use -- and comes back in the GRID's
        CRS, whatever the zone file's own was. See
        :func:`_resolve_drn_zone_cells` for what that rule replaced and for what a
        zone landing on nothing does.

        Parameters
        ----------
        model : SimulationBase or None
            The model whose grid (``model.vor``) and active domain geometry zones
            are resolved against. Not read when the zones carry a ``cells``
            column, which is returned as given; ``None`` then works.

        Returns
        -------
        pandas.DataFrame or geopandas.GeoDataFrame
            One row per ``locations`` row: ``name``, ``cells`` (zero-based cell
            numbers in the row's layer; for a geometry zone, the union over
            every row of that name), ``group``, ``weight``, ``layer``, and for
            geometry zones their geometry, converted to the grid's CRS when both
            CRSs are known.

        Raises
        ------
        ValueError
            When zones without a ``cells`` column were not given as a
            ``GeoDataFrame`` or vector file; or when no zone lands on a cell
            active in its layer, naming the zones (the first ten) and both
            extents and CRSs (ledger 211). A PARTIAL miss is logged at WARNING instead, once per
            missed set on this object.
        """

        return _resolve_drn_zone_cells(self._locations, model, owner=self)

    def get(self, model=None) -> pd.DataFrame:
        """Return the zone definitions (cells resolved when ``model`` is given), merged with values.

        Parameters
        ----------
        model : SimulationBase or None, optional
            When given, zones are resolved through :meth:`zone_definitions`,
            whose ``Raises`` apply. Default ``None``: the zones as normalized at
            construction -- explicit cells only, geometry unresolved and in the
            zone file's own CRS.

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
        """Return the observed discharge in wide format with one column per zone."""

        if self._values.empty:
            return pd.DataFrame(columns=["time", *self._locations["name"].astype(str).tolist()])
        frame = self._values.pivot_table(index="time", columns="name", values="flow_target", aggfunc="first").reset_index()
        frame.columns.name = None
        return frame

    def summary(self, model=None) -> pd.DataFrame:
        """Return a compact summary (zone count, resolved cell count, row count).

        Parameters
        ----------
        model : SimulationBase or None, optional
            When given, zones are resolved through :meth:`zone_definitions`
            first, whose ``Raises`` apply. Default ``None``: geometry zones are
            unresolved and count no cells.

        Returns
        -------
        pandas.DataFrame
            One row: ``n_zones`` (distinct names), ``n_cells`` (the length of
            every row's cell list, added over rows -- so a geometry zone given as
            several rows, each carrying the zone's union, is counted once per
            row) and ``n_rows`` (observed values).
        """

        definitions = self.zone_definitions(model) if model is not None else self._locations.copy()
        zone_count = int(definitions["name"].nunique())
        cell_count = int(sum(len(value) for value in definitions.get("cells", pd.Series(dtype=object))))
        return pd.DataFrame([{"n_zones": zone_count, "n_cells": cell_count, "n_rows": int(len(self._values))}])

    def _drn_budget_frame(self, model) -> pd.DataFrame:
        """The DRN cell-by-cell budget as a tidy per/node/q frame (adds ``per`` if only kstpkper)."""

        budget = model.bud("drn").df.reset_index().copy()
        if "per" not in budget.columns:
            if "kstpkper" not in budget.columns:
                raise ValueError("DRN budget dataframe requires 'per' or 'kstpkper' columns.")
            budget["per"] = budget["kstpkper"].apply(lambda item: int(item[1]))
        if "node" not in budget.columns:
            raise ValueError("DRN budget dataframe requires a zero-based 'node' column.")
        budget["node"] = pd.to_numeric(budget["node"], errors="coerce").astype(int)
        budget["q"] = pd.to_numeric(budget["q"], errors="coerce")
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

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package and a cell-budget file.
        zones : pandas.DataFrame
            The resolved zones (:meth:`zone_definitions`): ``name`` and ``cells``.

        Returns
        -------
        pandas.DataFrame
            ``per``, ``name`` and ``sim_flow`` (MF6's sign: a drain's outflow is
            negative), one row per zone with cells per period the file saved.

        Raises
        ------
        ValueError
            When a period the budget file saved lacks its final step -- naming
            the period, the steps saved and ``SAVE BUDGET LAST``.
        """

        budget = self._drn_budget_frame(model)
        members = _zone_members(zones)
        names = list(members)
        pieces = [
            budget.loc[budget["node"].isin(cells)].assign(zone=index)
            for index, cells in enumerate(members.values())
            if cells
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

    def _obs_column_map(self, model) -> tuple[pd.DataFrame, dict[str, list[str]]]:
        """The resolved zones and, per zone NAME, the MF6 observation label of each member cell.

        A zone given as several rows is ONE zone here as on the budget-file route
        (:func:`_zone_members`): its labels are those of the union of its rows'
        cells, each once. It was keyed per row by assignment, so the LAST row won
        and a zone of rows ``[20]`` and ``[21, 27]`` summed cells 21 and 27 only on
        this route while the budget route summed all three (ledger 336).

        Parameters
        ----------
        model : SimulationBase
            The run model the zones are resolved against.

        Returns
        -------
        tuple of (pandas.DataFrame, dict of str to list of str)
            The resolved zones (:meth:`zone_definitions`), and each zone's labels
            in ascending cell order.
        """

        zones = self.zone_definitions(model)
        zone_columns = {
            name: [_observation_label(f"{name}_c{cell}") for cell in sorted(cells)]
            for name, cells in _zone_members(zones).items()
        }
        return zones, zone_columns

    def _simulated_series_from_obs(self, model) -> pd.DataFrame | None:
        """Sum per-cell observation-CSV outputs into one wide series per zone, or ``None`` if absent."""

        zones, zone_columns = self._obs_column_map(model)
        expected = {
            f"{zone_name}__{index}": column
            for zone_name, columns in zone_columns.items()
            for index, column in enumerate(columns)
        }
        observed = _find_named_observation_output(model, expected_columns=expected)
        if observed is None:
            return None

        result = pd.DataFrame({"time": observed["time"]})
        for zone_name, columns in zone_columns.items():
            if not columns:
                result[zone_name] = 0.0
                continue
            zone_frame = pd.DataFrame(
                {
                    column: _clean_observation_output_values(
                        observed[f"{zone_name}__{index}"]
                    )
                    for index, column in enumerate(columns)
                }
            )
            result[zone_name] = zone_frame.sum(axis=1)
        return result

    def _compare_from_obs(self, model) -> pd.DataFrame | None:
        """Build the target-vs-simulated compare table from observation output, or ``None``."""

        wide = self._simulated_series_from_obs(model)
        if wide is None:
            return None

        zones = self.zone_definitions(model)
        simulated = wide.melt(id_vars=["time"], var_name="name", value_name="sim_flow")
        simulated["per"] = simulated["time"]
        simulated = zones.merge(simulated, on="name", how="left")
        if self._values.empty:
            compare = simulated.copy()
            compare["flow_target"] = np.nan
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

    def simulated_series(self, model) -> pd.DataFrame:
        """Return one wide simulated-discharge table: a ``time`` column and one column per zone.

        Read from the model's DRN observations when one workspace CSV holds a
        column for every member cell (MF6's labels, as :meth:`to_flopy_obs` writes
        them) -- each row is one time step at MF6's simulated time, the zone the
        sum of its cells in that step -- else from the model's cell-budget file.

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
            its grid.

        Returns
        -------
        pandas.DataFrame
            ``time`` plus one float column per zone, in MF6's sign (a drain's
            outflow is negative). On the budget-file route a zone with no member
            cell gets no column.

        Raises
        ------
        ValueError
            On the budget-file route, when a period the budget file saved lacks
            its final step -- naming the period, the steps saved and the ``SAVE
            BUDGET LAST`` setting that fixes it; or when no zone lands on an
            active cell (ledger 211).

        Examples
        --------
        >>> seeps = mf.DrnFlowTargets(locations={"north": [317, 342]})
        >>> series = seeps.simulated_series(model)   # one row per stress period
        """

        observed = self._simulated_series_from_obs(model)
        if observed is not None:
            return observed

        zones = self.zone_definitions(model)
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
        their ``time`` as a zero-based integer -- compared with the stress period
        on the budget-file route, and with MF6's simulated time on the
        observation route, where a target keyed by period matches no row.

        Parameters
        ----------
        model : SimulationBase
            A run model with a DRN package, as ``mf.load_mf6_run(ws)`` or
            ``run.model()`` returns.

        Returns
        -------
        pandas.DataFrame
            One row per target row (per simulated row when the targets hold no
            values): ``name``, ``time``, ``per``, ``flow_target``, ``sim_flow``,
            ``residual`` (simulated minus target), ``abs_residual``, and the
            zone's ``group``, ``weight`` and ``cells``. On the budget-file route
            with no values there is no ``per`` column: the period is ``time``,
            and ``flow_target`` and the residuals are ``NaN``.

        Raises
        ------
        ValueError
            On the budget-file route, when a period the budget file saved lacks
            its final step (as :meth:`simulated_series`); or when a target's
            ``time`` is not a whole number, or two rows hold a value for one
            target in one period (the period rule PEST's build reads too, ledger
            340).
        """

        observed = self._compare_from_obs(model)
        if observed is not None:
            return observed

        zones = self.zone_definitions(model)
        simulated = self._zone_flow_table(model, zones)
        details = pd.DataFrame(
            [
                {
                    "name": str(row.name),
                    "group": getattr(row, "group", pd.NA),
                    "weight": getattr(row, "weight", np.nan),
                }
                for row in zones.itertuples(index=False)
            ]
        ).drop_duplicates(subset="name")
        members = _zone_members(zones)
        details["cells"] = [sorted(members[name]) for name in details["name"]]
        simulated = simulated.merge(details, on="name", how="left")
        if self._values.empty:
            compare = simulated.rename(columns={"per": "time"}).copy()
            compare["flow_target"] = np.nan
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
            ``abs_residual``, ``weight`` and, when present, ``cells``.

        Raises
        ------
        KeyError
            Naming ``per``, when the targets hold no ``values`` and the
            budget-file route is taken: :meth:`compare` then returns no ``per``
            column.
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
            Default ``"drn_flow_targets.csv"``. :meth:`simulated_series` reads
            observations back from a workspace CSV holding a column for every
            member cell, and skips any file whose name contains ``simulated_`` or
            ``_target_``.
        kind : str, optional
            The MF6 observation type written on every record, upper-cased.
            Default ``"DRN"``, MF6's drain-flow observation, which carries MF6's
            sign (outflow negative); the observation route of
            :meth:`simulated_series` sums whatever that CSV holds as the flow.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, kind, (layer, cell)), ...]}``: one record per
            cell of every zone row, labelled ``<name>_c<cell>`` with every
            character that is not a letter or digit made ``_`` and leading
            ``_`` stripped, the cellid
            ZERO-based as FloPy takes it. Records follow the ROWS, so a cell on
            two rows of one zone -- every cell of a geometry zone given as
            several rows -- is written once per row under one label. The shape
            ``ModflowUtlobs(continuous=)`` and a DRN package's
            ``observations=`` take.
        """

        zones = self.zone_definitions(model)
        records: list[tuple[str, str, tuple[int, int]]] = []
        for row in zones.itertuples(index=False):
            for cell in getattr(row, "cells", []):
                obsname = _observation_label(f"{row.name}_c{int(cell)}")
                records.append((obsname, str(kind).upper(), (int(getattr(row, "layer", 0)), int(cell))))
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
            The model whose FloPy DRN package, ``model.gwf.drn``, receives the
            OBS package unless ``package`` is given; also passed to
            :meth:`to_flopy_obs` to resolve the zones. With ``package`` given and
            every zone's ``cells`` explicit, ``None`` works.
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
            The MF6 observation type, as in :meth:`to_flopy_obs`. Default
            ``"DRN"``.
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            ``model.gwf.drn``. Default ``None``.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

        import flopy

        owner = model.gwf.drn if package is None else package
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


