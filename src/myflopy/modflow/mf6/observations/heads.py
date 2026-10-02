"""Head observation targets (`HeadTargets` + bound views/plots)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from myflopy._logging import get_logger
from myflopy.grid_identity import excavation_cells
from myflopy.modflow.mf6.grid.domain import model_domain
from myflopy.modflow.mf6.grid.geometry import cell_shares, domain_active, miss_reason
from myflopy.modflow.mf6.observations._shared import (
    _build_values_from_cells,
    _coerce_label_list,
    _coerce_scalar_or_list,
    _copy_frame,
    _default_obs_csv_name,
    _ensure_unique_observation_names,
    _in_grid_crs,
    _load_locations,
    _measured_periods,
    _normalize_identifier_series,
    _normalize_values_frame,
    _numeric_periods,
    _placement_note,
    _report_missed_features,
    _simulated_field_by_period,
)
from myflopy.viz import mpl_axes

logger = get_logger(__name__)


def _dug_cells(model) -> dict[int, str]:
    """``{cell: excavation name}`` for whatever records this model's context or grid carry.

    Only used to NAME why a target's cell is not active (ledger 321): a cell inside an
    excavation is lake water, and the head there belongs to the cell under the bed.
    """

    import types

    found: dict[int, str] = {}
    for context in (
        getattr(getattr(model, "gwf", None), "myflopy_context", None),
        getattr(model, "myflopy_context", None),
        types.SimpleNamespace(grid=getattr(model, "vor", None)),
    ):
        for cell, name in excavation_cells(context).items():
            found.setdefault(cell, name)
    return found


_DUG_ADVICE = (
    " A cell inside an excavation is lake water, not aquifer -- observe the head in the "
    "cell under its bed, the column's top active layer."
)


@dataclass
class HeadTargets:
    """Observation locations plus target heads through time or stress period.

    Parameters
    ----------
    locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, or list or tuple of dict
        Where each target is, one row per target, read at construction:

        - a path (``str`` or ``Path``) is read with ``geopandas.read_file``;
        - a ``GeoDataFrame`` locates each target by its POINT geometry, matched
          to the one cell holding it by :meth:`match_to_model` (reprojected to
          the grid's CRS first when both carry one). An empty geometry raises
          here; a line or polygon raises at matching;
        - a ``DataFrame`` or a list or tuple of record dicts must carry a ``cell`` (or
          ``node``) column: the zero-based cell index WITHIN a layer (the DISV
          ``cell2d`` number, ``0 .. ncpl-1``). ``node`` is read as that same
          index, not as MF6's layered node number;
        - a ``Series`` maps each name (its index) to a cell index;
        - a ``dict`` is read as columns when one of its keys is ``name``,
          ``cell`` or ``node`` (any case), and otherwise as ``{name: cell}``.

        A ``cell``/``node`` column wins over a geometry when both are present.
        Only the ``name_column``, ``cell``/``node``, ``layer_column``,
        ``group_column`` and ``weight_column`` columns (and the geometry) are
        kept. Raises ``ValueError`` when ``name_column`` is missing, when a
        table with no geometry has no ``cell``/``node`` column, or when a cell
        value is not numeric.
    values : str, pathlib.Path, pandas.DataFrame, pandas.Series, dict, list or tuple
        Target values through time. A path is read with ``pandas.read_csv``; a
        ``dict`` is read as columns and a list/tuple as rows (``times`` can
        supply the time column for either). Three table shapes are recognised:

        - LONG: a ``name`` column plus a ``value_column`` column, one row per
          target per time (both names matched ignoring case);
        - one series: only ``time_column`` and ``value_column`` and no ``name``,
          which is allowed only when there is exactly one location;
        - WIDE: every column other than ``time_column`` is one target, named by
          its header.

        A ``Series`` is one target's values indexed by time (so it also needs
        exactly one location). Every name must match a location name, compared
        stripped and lower-cased, or ``ValueError`` is raised. Values are read
        as floats: a BLANK entry (missing, or empty text) is ``NaN`` and means
        unmeasured -- the PEST build leaves that row out of the targets, at
        weight 0 -- and any other entry that is not a finite number raises
        ``ValueError`` naming the target and its column (ledger 340; it used to
        become ``NaN`` silently).
    name_column : str, optional
        Column of ``locations`` holding each target's name. Default ``"name"``.
        Required: ``ValueError`` when absent. For a ``Series`` or a
        ``{name: cell}`` dict it is the name given to the index or the keys.
        It applies to ``locations`` only -- a LONG ``values`` table always names
        its targets in a column called ``name``.
    layer_column : str or None, optional
        Column of ``locations`` holding each target's zero-based model layer.
        Default ``"layer"``. ``None``, or a column that is not there, puts
        every target in layer 0, and a blank or non-numeric entry is read as 0.
        A target only matches a cell that is ACTIVE in its own layer
        (:meth:`match_to_model`).
    group_column : str or None, optional
        Column of ``locations`` holding each location's observation-group label,
        carried as ``group`` into :meth:`get`, :meth:`residuals`,
        :meth:`summary`'s ``n_groups`` and the PEST locations snapshot. Default
        ``"group"``; ``None`` or an absent column leaves ``group`` as ``<NA>``.
        **It is the PEST observation group** (ledger 358): when
        ``cal.observe(...)`` builds the control file, each MEASURED observation
        at the location is in the group the label names -- lower-cased, and
        every run of characters other than letters, digits, ``_``, ``.``,
        ``:`` and ``-`` turned into one ``_`` (``"North Wells"`` is
        ``north_wells``), so phi is reported, and PEST++'s phi-factor file can
        weight, by your groups. A missing or blank label puts it in the set's
        prefix (``hds``, ``conc``); a period with no value is in
        ``<prefix>_unmeasured`` whatever the label, and a forecast is in its set's
        prefix. The build refuses a label that is not text (a number), holds no
        letter or digit, starts as PESTPP-IES reads an inequality (``l_``,
        ``less``, ``g_``, ``greater``), takes a group the build gives a set
        (``hds``, ``stage_unmeasured``, ...), or becomes the same name as a
        different label (``"North"`` and ``"north"``) -- naming each.
    weight_column : str or None, optional
        Column of ``locations`` holding a per-location observation weight,
        carried as ``weight``. When ``cal.observe(...)`` builds the PEST control
        file, a weight that is not ``NaN`` becomes the weight of each
        observation at that location that has a target value, and ``NaN``
        gives those observations weight 1; a period with no target value
        carries weight 0 (ledger 340). Default ``"weight"``; ``None`` or an
        absent column leaves it ``NaN``.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``. Required once ``times`` has been applied
        (``ValueError`` otherwise); :meth:`get` calls it ``time``.
        :meth:`compare` and the PEST build read it as a ZERO-BASED STRESS
        PERIOD, by number (``1.0`` is period 1), and both raise on a label that
        is not a whole number, and on two values for one target in one period.
    value_column : str, optional
        Column of a LONG ``values`` table holding the target value, matched
        ignoring case. Default ``"head"`` (``"conc"`` on
        :class:`~myflopy.modflow.mf6.observations.conc.ConcTargets`). Unused
        by a WIDE table; for a one-column list given with ``times`` it names
        that column.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels, one per row, for ``values`` given as a ``dict`` or a
        list/tuple with no ``time_column`` column; they are inserted as that
        column. A list/tuple of another length raises ``ValueError``, and so
        does a dict (from pandas). Ignored for a ``DataFrame``, a ``Series`` or
        a CSV path, and whenever ``values`` already has the column. Default
        ``None``.

    Notes
    -----
    A target keyed by stress period reads the value at that period's END (ledger
    336) -- the record that ends the period on the head (or concentration) file's
    own clock, the one rule the stage and flow targets read too. A period whose
    final step output control did not save (``SAVE HEAD FIRST``, ``STEPS`` short
    of it) is refused: its end-of-period value does not exist in the file, and it
    used to be read at the last step saved, wherever that fell.
    """

    #: model attribute holding the full value table for this field
    _TABLE_ATTRIBUTE: ClassVar[str] = "all_heads"
    #: the value column inside that table
    _STORE_COLUMN: ClassVar[str] = "elev"
    #: human-readable field name, used in errors and plot labels
    _FIELD_LABEL: ClassVar[str] = "head"

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple
    name_column: str = "name"
    layer_column: str | None = "layer"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "head"
    times: list | tuple | pd.Index | pd.Series | None = None

    def __post_init__(self):
        """Load and normalize locations and target tables."""

        self._locations = _load_locations(
            self.locations,
            name_column=self.name_column,
            layer_column=self.layer_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
        )
        self._values = _normalize_values_frame(
            self.values,
            time_column=self.time_column,
            value_column=self.value_column,
            times=self.times,
            owner=type(self).__name__,
        )
        if "name" not in self._values.columns:
            if len(self._locations["name_key"].unique()) != 1:
                raise ValueError(
                    "HeadTargets values without explicit names require exactly one target location."
                )
            only_name = str(self._locations["name"].iloc[0])
            self._values["name"] = only_name
            self._values["name_key"] = _normalize_identifier_series(self._values["name"])
        missing = sorted(
            set(self._values["name_key"].unique()) - set(self._locations["name_key"].unique())
        )
        if missing:
            raise ValueError(
                "Head target values reference locations not found in the locations layer: "
                + ", ".join(missing[:10])
            )

    @classmethod
    def from_cells(
        cls,
        *,
        cells,
        values,
        times=None,
        names=None,
        layers=0,
        groups=None,
        weights=None,
        time_column: str = "time",
        value_column: str = "head",
    ):
        """Create head targets directly from explicit model cell ids.

        Parameters
        ----------
        cells : sequence of int
            One zero-based cell index per target, WITHIN its layer (the DISV
            ``cell2d`` number, ``0 .. ncpl-1``); each is passed through ``int``.
            At least one is required (``ValueError``). The ids are checked
            against the grid when the targets are matched (:meth:`match_to_model`):
            an id the grid does not have, or one inactive in its layer, raises.
        values : array-like, dict, pandas.DataFrame or pandas.Series
            The target values:

            - a 2-D array-like of shape ``(n_times, len(cells))``: column ``j``
              is target ``j`` (another column count raises ``ValueError``);
            - a 1-D array-like: one target's values through time, allowed only
              with ONE cell;
            - a ``Series``: one target (ONE cell only), indexed by time;
            - a ``dict`` or ``DataFrame``: a table in any shape the constructor's
              ``values`` accepts -- columns named by target (``names``), or LONG
              with ``name`` and ``value_column`` columns. A ``DataFrame`` whose
              non-time columns are labelled ``0 .. n-1`` has them renamed to
              ``names`` in order.

            Anything with more than two dimensions raises ``ValueError``.
        times : sequence or None, optional
            Time labels, one per row of ``values``, becoming ``time_column``.
            Default ``None``: ``0, 1, 2, ...`` for array input, so each row is
            read as that zero-based stress period. Required for a ``dict`` or
            ``DataFrame`` that has no ``time_column`` column (``ValueError``);
            ignored for a ``Series`` (its index is the time) and for a table
            that already has the column.
        names : str, sequence of str or None, optional
            One name per cell, in the order of ``cells``. Default ``None`` names
            them ``OBS_001``, ``OBS_002``, ... A single string is allowed only
            for one cell; a sequence of another length raises ``ValueError``.
        layers : int, sequence of int or None, optional
            Zero-based model layer of each target. Default ``0``. A scalar
            applies to every cell and ``None`` means 0; a list, tuple, array or
            ``Series`` must have one entry per cell (``ValueError``).
        groups : object, sequence or None, optional
            Observation-group label per target -- its PEST observation group
            (see the constructor's ``group_column``, which says how a label
            becomes a group name and which labels a PEST build refuses: anything
            but text among them). A scalar, a string included, applies to every
            cell; a sequence must have one entry per cell. Default ``None``
            leaves ``group`` as ``<NA>``: the set's prefix in PEST.
        weights : float, sequence of float or None, optional
            PEST observation weight per target (see the constructor's
            ``weight_column``). A scalar applies to every cell; a sequence must
            have one entry per cell. Default ``None`` leaves ``weight`` as
            ``NaN``.
        time_column : str, optional
            Name of the time column in the table built here, and stored on the
            result. Default ``"time"``.
        value_column : str, optional
            Value column of a LONG ``dict``/``DataFrame`` given as ``values``,
            stored on the result. Default ``"head"`` -- also when called on
            :class:`~myflopy.modflow.mf6.observations.conc.ConcTargets`, whose
            class default is ``"conc"``. Unused by array and wide input.

        Returns
        -------
        HeadTargets
            An instance of the class it is called on (``ConcTargets`` included).
        """

        cell_list = [int(value) for value in list(cells)]
        count = len(cell_list)
        if count == 0:
            raise ValueError("from_cells(...) requires at least one cell.")
        names_list = _coerce_label_list(names, count=count, default_prefix="OBS")
        layer_list = [int(value) for value in _coerce_scalar_or_list(layers, count=count, default=0)]
        group_list = _coerce_scalar_or_list(groups, count=count, default=pd.NA)
        weight_list = _coerce_scalar_or_list(weights, count=count, default=np.nan)
        locations = pd.DataFrame(
            {
                "name": names_list,
                "layer": layer_list,
                "cell": cell_list,
                "group": group_list,
                "weight": weight_list,
            }
        )
        values_frame = _build_values_from_cells(
            values,
            names=names_list,
            times=times,
            time_column=time_column,
            value_column=value_column,
        )
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
        )

    @classmethod
    def from_records(
        cls,
        records,
        *,
        time_column: str = "time",
        value_column: str = "head",
    ):
        """Create head targets from long-format record dictionaries/dataframes.

        Parameters
        ----------
        records : list of dict, pandas.DataFrame or GeoDataFrame
            One row per target per time, as ``pandas.DataFrame(records)`` reads
            it. Required columns, matched exactly: ``name``, ``time_column``,
            ``value_column``, and a location given as ``cell`` or ``node`` (the
            zero-based cell index within a layer) or as a point ``geometry``.
            Optional: ``layer`` (zero-based, default 0), ``group``, ``weight``.
            A target's location is taken from its FIRST row. A ``geometry``
            column of plain shapely points (record dicts, a ``DataFrame``)
            carries no CRS, so its coordinates are read in the grid's own
            frame; a ``GeoDataFrame``'s CRS is kept, and reprojected to the
            grid's at matching when the grid has one. ``cell``/``node`` wins
            over a geometry when
            both are given. ``ValueError`` when there are no records or a
            required column is missing.
        time_column : str, optional
            Column holding each record's time label, read by :meth:`compare` as
            a zero-based stress period. Default ``"time"``.
        value_column : str, optional
            Column holding each record's target value. Default ``"head"``
            (also on ``ConcTargets``).

        Returns
        -------
        HeadTargets
            An instance of the class it is called on.
        """

        frame = pd.DataFrame(records).copy()
        if frame.empty:
            raise ValueError("from_records(...) requires at least one record.")
        if "name" not in frame.columns:
            raise ValueError("from_records(...) requires a 'name' column.")
        if time_column not in frame.columns:
            raise ValueError(f"from_records(...) requires the time column {time_column!r}.")
        if value_column not in frame.columns:
            raise ValueError(f"from_records(...) requires the value column {value_column!r}.")
        if "cell" not in frame.columns and "geometry" not in frame.columns and "node" not in frame.columns:
            raise ValueError("from_records(...) requires 'cell', 'node', or 'geometry' columns.")
        locations_columns = [
            column
            for column in ["name", "layer", "cell", "node", "group", "weight", "geometry", "x", "y"]
            if column in frame.columns
        ]
        locations = frame.loc[:, locations_columns].drop_duplicates(subset=["name"]).copy()
        if "geometry" in locations.columns:
            locations = gpd.GeoDataFrame(locations, geometry="geometry")
        values_frame = frame.loc[:, [time_column, "name", value_column]].copy()
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
        )

    @property
    def locations_gdf(self) -> gpd.GeoDataFrame:
        """Normalized observation locations table."""

        return self._locations.copy()

    def get(self) -> pd.DataFrame:
        """Return the merged long-format target table."""

        if isinstance(self._locations, gpd.GeoDataFrame):
            location_frame = self._locations.drop(columns=self._locations.geometry.name)
        else:
            location_frame = self._locations.copy()
        merged = self._values.merge(
            location_frame,
            on="name_key",
            how="left",
            suffixes=("", "_loc"),
        )
        if "name_loc" in merged.columns:
            merged = merged.drop(columns=["name"]).rename(columns={"name_loc": "name"})
        return _copy_frame(merged)

    def to_long(self) -> pd.DataFrame:
        """Return the long-format target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per location."""

        frame = self.get().pivot_table(
            index="time",
            columns="name",
            values="head_target",
            aggfunc="first",
        ).reset_index()
        frame.columns.name = None
        return frame

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the target dataset."""

        frame = self.get()
        summary = {
            "n_locations": int(frame["name_key"].nunique()),
            "n_rows": int(len(frame)),
            "n_groups": int(frame["group"].dropna().nunique()) if "group" in frame else 0,
            "n_weighted": int(frame["weight"].notna().sum()) if "weight" in frame else 0,
        }
        return pd.DataFrame([summary])

    def plot_locations(self, *, ax=None, **kwargs):
        """Plot target locations and return the matplotlib axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or None, optional
            Axes to draw on. Default ``None`` makes a new figure and axes.
        **kwargs
            Forwarded unchanged to ``geopandas.GeoDataFrame.plot`` (for example
            ``color=``, ``markersize=``).

        Returns
        -------
        matplotlib.axes.Axes
            The axes, titled "Head Target Locations" (on ``ConcTargets`` too).

        Raises
        ------
        ValueError
            When the locations carry no geometry -- targets given as cell ids.
        """

        if not isinstance(self._locations, gpd.GeoDataFrame):
            raise ValueError("plot_locations() requires geospatial target locations.")
        if ax is None:
            _, ax = mpl_axes()
        self._locations.plot(ax=ax, **kwargs)
        ax.set_title("Head Target Locations")
        return ax

    def match_to_model(self, model) -> pd.DataFrame:
        """Map target locations to zero-based model cell ids.

        Parameters
        ----------
        model : SimulationBase
            The model to match against -- a run (``run.model()``,
            ``mf.load_mf6_run(ws)``) or any model-like object exposing
            ``vor.gdf_vorPolys``, the grid's cell polygons. Which cells are
            active comes from its discretization's ``idomain``, else its
            ``myflopy_context.domain``; with neither, every cell is active.

        Returns
        -------
        pandas.DataFrame
            One row per MATCHED target location, with its ``cell`` and ``layer``.

        Notes
        -----
        A point resolves to the ONE cell that holds it: on a shared face or
        vertex, the lowest-numbered cell active in the target's own layer
        (ledger 233). It used to reach every cell it touched -- two on a face,
        four on a vertex -- so one well became four observations, and PEST's
        build and its forward run read different cells for it.

        Targets are reprojected to the grid's CRS first. One that lands on no
        active cell is left out and named in a warning -- a Python
        ``UserWarning`` as well as the log (ledger 342) -- once per missed set on
        this object; if NONE of them land,
        that raises, rather than handing back an observation set that is quietly
        about somewhere else. A PEST build refuses the one left out when the
        values name it (ledger 341): there it would be a target PEST never sees. A target given as a cell id is checked instead: an
        id this grid does not have, or one inactive in its layer, raises.
        """

        if "cell" in self._locations.columns:
            matches = self._validated_cells(model)
            if matches["x"].isna().any() or matches["y"].isna().any():
                vor = model.vor.gdf_vorPolys.reset_index().rename(columns={"index": "cell"})
                centroids = pd.DataFrame(
                    {
                        "cell": vor["cell"].astype(int),
                        "x": vor.geometry.centroid.x,
                        "y": vor.geometry.centroid.y,
                    }
                )
                matches = matches.drop(columns=["x", "y"], errors="ignore").merge(centroids, on="cell", how="left")
            return _copy_frame(
                matches.loc[:, ["name", "name_key", "layer", "group", "weight", "cell", "x", "y"]]
            )

        return self._matched_from_geometry(model)

    def _validated_cells(self, model) -> pd.DataFrame:
        """The locations, with their explicit cell ids checked against the grid.

        An id is an assertion about the model, so a wrong one raises rather than
        being warned about the way a missed point is. Before this, an id the grid
        did not have reached the value table as NaN and dropped out of ``stats``
        silently, and an inactive id was read where MF6 writes no value.
        """

        matches = self._locations.copy()
        ncpl = len(model.vor.gdf_vorPolys)
        domain = model_domain(model, ncpl)
        cells = matches["cell"].astype(int)
        layers = matches["layer"].astype(int)
        outside = [
            f"{name!r} (cell {cell})"
            for name, cell in zip(matches["name"], cells, strict=True)
            if not 0 <= cell < ncpl
        ]
        if outside:
            raise ValueError(
                f"These {self._FIELD_LABEL} targets name cells this grid does not "
                f"have (it has {ncpl}, numbered from 0): {', '.join(outside)}."
            )
        dug = _dug_cells(model)
        inactive = [
            f"{name!r} (cell {cell}, layer {layer}"
            + (f", inside excavation {dug[cell]!r}" if cell in dug else "")
            + ")"
            for name, cell, layer in zip(matches["name"], cells, layers, strict=True)
            if not domain_active(domain, layer, cell)
        ]
        if inactive:
            raise ValueError(
                f"These {self._FIELD_LABEL} targets name cells that are not active "
                f"in their own layer, where MF6 writes no {self._FIELD_LABEL}: "
                f"{', '.join(inactive)}."
                + (_DUG_ADVICE if any("inside excavation" in entry for entry in inactive) else "")
            )
        return matches

    def _matched_from_geometry(self, model) -> pd.DataFrame:
        """Resolve each target's point to the one cell it sits in."""

        polys = model.vor.gdf_vorPolys
        locations = self._locations_in_grid_crs(model)
        domain = model_domain(model, len(polys))
        columns = ["name", "name_key", "layer", "group", "weight", "cell", "x", "y"]
        rows: list[dict[str, Any]] = []
        missed: list[str] = []
        dug: dict[int, str] | None = None
        for _, row in locations.iterrows():
            geometry = row.geometry
            name = str(row["name"])
            if geometry is not None and not geometry.is_empty and not isinstance(geometry, Point):
                raise ValueError(
                    f"{self._FIELD_LABEL.capitalize()} target {name!r} is a "
                    f"{geometry.geom_type}. A {self._FIELD_LABEL} is measured at a "
                    f"POINT: give its coordinates, or the cell it sits in."
                )
            layer = int(row["layer"])
            hits = cell_shares(
                polys,
                geometry,
                # Eligibility only. ``cell_shares`` sorts whatever this returns,
                # so the lowest-numbered usable cell wins a tie without a second
                # copy of that rule here.
                usable=lambda cells, layer=layer: [
                    cell for cell in cells if domain_active(domain, layer, cell)
                ],
            )
            if not hits:
                reason = miss_reason(polys, geometry, target=f"active cell in layer {layer}")
                if dug is None:
                    dug = _dug_cells(model)
                inside = [cell for cell, _share in cell_shares(polys, geometry)] if dug else []
                if inside and inside[0] in dug:
                    reason += (
                        f", inside excavation {dug[inside[0]]!r} -- the head there is the "
                        "cell under its bed, in the column's top active layer"
                    )
                missed.append(f"{name!r} ({reason})")
                continue
            rows.append({column: row.get(column) for column in columns} | {"cell": hits[0][0]})
        if missed:
            self._report_missed(model, missed, total=len(locations))
        frame = pd.DataFrame(rows, columns=columns)
        frame["cell"] = frame["cell"].astype(int)
        frame["layer"] = frame["layer"].astype(int)
        return _copy_frame(frame)

    def _locations_in_grid_crs(self, model) -> gpd.GeoDataFrame:
        """The locations in the GRID's CRS, with x/y recomputed from them.

        Targets are read in whatever CRS their file carries, and the grid decides
        where a cell is. Before this, targets in degrees over a State Plane grid
        were located by their raw degrees: two wells whose cells were 5 and 10
        both resolved to cell 4.
        """

        locations = self._locations
        reprojected = _in_grid_crs(locations, model)
        if reprojected is locations:
            return locations
        reprojected["x"] = reprojected.geometry.x
        reprojected["y"] = reprojected.geometry.y
        return reprojected

    def _report_missed(self, model, missed: list[str], *, total: int) -> None:
        """Raise when no target lands on an active cell; warn, by name, when only some do.

        A partial miss can be deliberate -- one target file serving several grids
        -- so it warns. None of them landing means the targets and the model are
        in different places, which is what the extents and CRSs in the message
        are read for. Every consumer (``simulated_heads``, ``compare``, ``stats``,
        ``to_flopy_obs``) calls ``match_to_model`` again, so the warning is said
        once per set of missed targets rather than once per call.
        """

        _report_missed_features(
            self,
            missed,
            total=total,
            subject=f"{self._FIELD_LABEL} target",
            nothing_left="compare against",
            where=lambda: self._where(model),
            log=logger,
        )

    def _where(self, model) -> str:
        """Where the targets are and where the grid is -- what an empty match is acted on with."""

        return _placement_note(self._locations, model)

    def observation_labels(self, model) -> pd.DataFrame:
        """The matched targets, each with the label its observation is known by.

        ONE definition of that label, shared by the MF6 ``continuous`` records,
        the simulated CSV's columns, the mapping the forward run reads, and the
        pst rows. A name used once IS its label; a name used in two layers becomes
        ``<name>_lay<layer>`` (then ``_c<cell>``), which is what the FloPy records
        have always done.

        The PEST side keyed on the bare name instead, so a well screened in two
        layers was one column: the build pivoted ``first`` and read layer 0, the
        forward run wrote each mapping row in turn and left layer 1, and the two
        disagreed by whatever the layers differ by (5.0 against 1005.0 measured).

        Parameters
        ----------
        model : SimulationBase
            The model the targets are matched against, as in
            :meth:`match_to_model` (a target on no active cell is left out and
            warned about; none landing raises).

        Returns
        -------
        pandas.DataFrame
            :meth:`match_to_model`'s table plus ``obsname``, the label.
        """

        matches = self.match_to_model(model).copy()
        matches["obsname"] = _ensure_unique_observation_names(matches)
        return matches

    def simulated_heads(self, model) -> pd.DataFrame:
        """Return simulated heads for the target locations by stress period.

        The result is a wide table suitable for pyEMU list-style observation
        setup: one row per stress period and one column per observation LABEL
        (see :meth:`observation_labels`, which a repeated name disambiguates).
        Each period's value is the one at its END (ledger 336), read on the
        head file's own clock.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or, for ``ConcTargets``, the
            concentration) file, as ``mf.load_mf6_run(ws)`` or ``run.model()``
            returns.

        Returns
        -------
        pandas.DataFrame
            ``per`` plus one column per observation label.

        Raises
        ------
        ValueError
            When a period the file saved lacks its final step -- naming the
            period, the steps saved and the ``SAVE HEAD LAST`` setting that
            fixes it; or when no target lands on an active cell (ledger 211).
        """

        matches = self.observation_labels(model)
        sim = _simulated_field_by_period(
            model,
            table_attribute=self._TABLE_ATTRIBUTE,
            store_column=self._STORE_COLUMN,
        )
        merged = matches.merge(sim, on=["layer", "cell"], how="left")
        wide = merged.pivot_table(
            index="per",
            # The label, not the bare name: one column per (name, layer), so a
            # two-screen well is two observations rather than one column whose
            # value depended on which of its rows got there first.
            columns="obsname",
            values="sim_head",
            aggfunc="first",
        ).reset_index()
        wide.columns.name = None
        return wide

    def compare(self, model) -> pd.DataFrame:
        """Compare target heads against simulated heads by stress period.

        A target's ``time`` is read as a ZERO-BASED STRESS PERIOD -- the one
        reading every target family and PEST's build share (ledger 340): ``1.0``
        is period 1, and a key with a fraction is refused rather than truncated
        (``10.9`` used to read period 10). ``sim_head`` is the head at each
        period's END (ledger 336), as :meth:`simulated_heads` reads it. A blank
        target value is a row with ``NaN`` ``head_target`` and ``residual``.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or concentration) file.

        Returns
        -------
        pandas.DataFrame
            One row per target row: the target's columns plus ``per``,
            ``cell``, ``sim_head``, ``residual`` (simulated minus target) and
            ``abs_residual``.

        Raises
        ------
        ValueError
            When a target's ``time`` is not a number or not a whole number; when
            two rows hold a value for one target in one period (PEST keeps one
            observation per target per period); or when a period the file saved
            lacks its final step (as :meth:`simulated_heads`).
        """

        caller = f"{type(self).__name__}.compare(...)"
        _measured_periods(
            self._values, value_column="head_target", name_column="name_key", caller=caller
        )
        targets = self.get().copy()
        targets["per"] = _numeric_periods(targets, caller=caller)
        sim = _simulated_field_by_period(
            model,
            table_attribute=self._TABLE_ATTRIBUTE,
            store_column=self._STORE_COLUMN,
        )
        if "cell" in targets.columns and targets["cell"].notna().all():
            merged = targets.merge(sim, on=["per", "layer", "cell"], how="left")
            if (targets["x"].isna().any() or targets["y"].isna().any()) and hasattr(model, "vor"):
                matches = self.match_to_model(model).loc[:, ["name", "name_key", "layer", "cell", "x", "y"]]
                merged = merged.drop(columns=["x", "y"], errors="ignore").merge(
                    matches,
                    on=["name", "name_key", "layer", "cell"],
                    how="left",
                )
        else:
            matches = self.match_to_model(model).loc[:, ["name", "name_key", "layer", "cell", "x", "y"]]
            merged = (
                targets.merge(matches, on=["name", "name_key", "layer"], how="left")
                .merge(sim, on=["per", "layer", "cell"], how="left")
            )
        merged["residual"] = merged["sim_head"] - merged["head_target"]
        merged["abs_residual"] = merged["residual"].abs()
        return _copy_frame(merged)

    def residuals(self, model) -> pd.DataFrame:
        """Return a compact residual table.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or, for ``ConcTargets``, the
            concentration) file; passed to :meth:`compare`, whose ``Raises``
            apply.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced to ``name``, ``group``, ``time``,
            ``per``, ``layer``, ``cell``, ``head_target``, ``sim_head``,
            ``residual`` (simulated minus target), ``abs_residual`` and
            ``weight``.
        """

        frame = self.compare(model)
        return frame.loc[
            :,
            [
                "name",
                "group",
                "time",
                "per",
                "layer",
                "cell",
                "head_target",
                "sim_head",
                "residual",
                "abs_residual",
                "weight",
            ],
        ]

    def stats(self, model) -> pd.DataFrame:
        """Return aggregate residual statistics for the target set.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or concentration) file; passed to
            :meth:`compare`, whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row: ``n`` (rows with a residual), ``mean_error``, ``mae`` and
            ``rmse`` of simulated minus target, in the field's units. ``n`` is 0
            and the rest ``NaN`` when no row has a residual.
        """

        frame = self.compare(model).dropna(subset=["residual"])
        if frame.empty:
            return pd.DataFrame(
                [{"n": 0, "mean_error": np.nan, "mae": np.nan, "rmse": np.nan}]
            )
        residual = frame["residual"].to_numpy(dtype=float)
        stats = {
            "n": int(len(frame)),
            "mean_error": float(np.mean(residual)),
            "mae": float(np.mean(np.abs(residual))),
            "rmse": float(np.sqrt(np.mean(residual**2))),
        }
        return pd.DataFrame([stats])

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        """Return a FloPy ``continuous`` dict for MF6 head observations.

        Parameters
        ----------
        model : SimulationBase
            The model the targets are matched against (:meth:`observation_labels`):
            one record per matched target, a target on no active cell left out
            with a warning, and none landing raises ``ValueError``.
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"head_targets.csv"``.
        kind : str, optional
            The MF6 observation type written on every record, upper-cased.
            Default ``"HEAD"``. ``ConcTargets`` inherits this default, which a
            GWT model's OBS file refuses ("Observation type not found: HEAD",
            measured): pass ``kind="CONCENTRATION"`` there.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(obsname, kind, (layer, cell)), ...]}``, the cell id
            zero-based (FloPy writes it one-based), ``obsname`` the label from
            :meth:`observation_labels`. The shape ``ModflowUtlobs(continuous=)``
            and a package's ``observations=`` take.
        """

        matches = self.observation_labels(model)
        records = [
            (str(row.obsname), str(kind).upper(), (int(row.layer), int(row.cell)))
            for row in matches.itertuples(index=False)
        ]
        return {str(csv_name): records}

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str = "HEAD",
        package=None,
    ):
        """Attach an MF6 utility observation package for these head targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase
            The model whose targets are matched (:meth:`to_flopy_obs`), and whose
            FloPy model object, ``model.gwf``, receives the package unless
            ``package`` is given. On a transport model that object is the GWT
            model.
        pname : str, optional
            FloPy package name. Default ``"gwf_obs"``.
        filename : str, optional
            The OBS input file name FloPy writes, relative to the owning model's
            directory. Default ``"head_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"head_targets.csv"`` for the
            default ``filename``).
        kind : str, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`. Default
            ``"HEAD"`` -- on ``ConcTargets`` too, where MF6 refuses it: pass
            ``kind="CONCENTRATION"``.
        package : flopy object or None, optional
            The FloPy model or package to attach the OBS package to instead of
            ``model.gwf``. Default ``None``.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

        import flopy

        owner = model.gwf if package is None else package
        csv_name = _default_obs_csv_name(filename, "head_targets.csv") if csv_name is None else str(csv_name)
        continuous = self.to_flopy_obs(model, csv_name=csv_name, kind=kind)
        return flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            owner,
            pname=pname,
            continuous=continuous,
            filename=str(filename),
        )

    def calibration_plot(self, model, *, type: str = "calibration"):
        """Build a calibration plot directly from these targets and one model.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or concentration) file; the plot is
            drawn from :meth:`compare`, whose ``Raises`` apply.
        type : str, optional
            Which plot:

            - ``"calibration"`` (default): target against simulated, one point
              per row, with fit statistics;
            - ``"obs_vs_sim"``: the same cross plot without the statistics;
            - ``"heads"``: target and simulated values through time, one pair of
              lines per location and layer;
            - ``"residuals_by_period"``: mean absolute residual per ``time``.

            Any other string draws the ``"calibration"`` plot.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        return CalibrationPlot.from_targets(self, model=model, type=type)


class BoundHeadTargets:
    """Model-bound head-target helper returned from ``model.targets.heads``."""

    def __init__(self, model, targets: HeadTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        self.model = model
        self.targets = targets
        self._plot = None

    def get(self) -> pd.DataFrame:
        """The target definitions (delegates to the wrapped :class:`HeadTargets`)."""

        return self.targets.get()

    def to_long(self) -> pd.DataFrame:
        """The long-format target table (delegates to the wrapped targets)."""

        return self.targets.to_long()

    def to_wide(self) -> pd.DataFrame:
        """The wide-format target table (delegates to the wrapped targets)."""

        return self.targets.to_wide()

    def summary(self) -> pd.DataFrame:
        """A compact summary of the target set (delegates to the wrapped targets)."""

        return self.targets.summary()

    def plot_locations(self, *, ax=None, **kwargs):
        """Plot the target locations on a map axis (delegates to the wrapped targets)."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def match_to_model(self, model=None) -> pd.DataFrame:
        """Match each target location to its grid cell (defaults to the bound model)."""

        return self.targets.match_to_model(self.model if model is None else model)

    def simulated_heads(self, model=None) -> pd.DataFrame:
        """Simulated heads at the target locations (defaults to the bound model)."""

        return self.targets.simulated_heads(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated heads (defaults to the bound model)."""

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the head targets (defaults to the bound model)."""

        return self.targets.residuals(self.model if model is None else model)

    def stats(self, model=None) -> pd.DataFrame:
        """Aggregate residual statistics (defaults to the bound model)."""

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "head_targets.csv",
        kind: str = "HEAD",
    ):
        """A FloPy ``continuous`` observation dict for the head targets on the bound model."""

        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str = "HEAD",
        package=None,
    ):
        """Attach an MF6 head observation package to the bound model."""

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the head targets on the bound model."""

        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        """Plotting helpers for this model-bound target set."""

        if self._plot is None:
            self._plot = BoundHeadTargetPlots(self)
        return self._plot


class BoundHeadTargetPlots:
    """Plotting facade for model-bound head targets."""

    def __init__(self, bound_targets: BoundHeadTargets):
        """Wrap a :class:`BoundHeadTargets` to expose its plotting helpers."""

        self.bound_targets = bound_targets

    @property
    def model(self):
        """The model the underlying targets are bound to."""

        return self.bound_targets.model

    @property
    def targets(self):
        """The underlying :class:`HeadTargets` being plotted."""

        return self.bound_targets.targets

    def locations(self, *, ax=None, **kwargs):
        """Plot target locations on a map axis."""

        return self.targets.plot_locations(ax=ax, **kwargs)

    def calibration(self, *, type: str = "calibration"):
        """Return a calibration/heads plot for the bound targets."""

        return self.bound_targets.calibration_plot(type=type)

    def obs_vs_sim(self, *, baseline=None, ax=None, backend: str = "plotly"):
        """Return a target-vs-simulated cross plot for one or two models.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(current, baseline_compare=baseline_frame, backend=backend)

    def timeseries(self, name: str | None = None, *, baseline=None, ax=None, backend: str = "plotly"):
        """Return a time-series plot for one target location.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
            backend=backend,
        )

    def residuals_by_period(self, *, baseline=None, backend: str = "plotly"):
        """Return a by-period residual summary plot.

        ``backend`` selects ``"plotly"`` (default) or ``"matplotlib"``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(
            current,
            baseline_compare=baseline_frame,
            backend=backend,
        )


