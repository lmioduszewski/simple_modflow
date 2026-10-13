"""Head observation targets (`HeadTargets` + bound views/plots)."""

from __future__ import annotations

import warnings
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
from myflopy.modflow.mf6.grid.geometry import (
    cell_shares,
    domain_active,
    miss_reason,
    model_cell_polygons,
)
from myflopy.modflow.mf6.hydraulics import model_flow_properties, saturated_thickness
from myflopy.modflow.mf6.observations._shared import (
    _FIELD_FILES,
    _build_values_from_cells,
    _coerce_label_list,
    _coerce_scalar_or_list,
    _copy_frame,
    _default_obs_csv_name,
    _ensure_unique_observation_names,
    _in_grid_crs,
    _is_blank,
    _listing,
    _load_locations,
    _measured_periods,
    _normalize_identifier_series,
    _normalize_values_frame,
    _numeric_periods,
    _placement_note,
    _refuse_keys_off_the_periods,
    _refuse_reserved_target_names,
    _report_missed_features,
    _simulated_field_by_period,
    _wide_by_period,
)
from myflopy.modflow.mf6.observations.hydrographs import (
    bin_readings,
    model_clock,
    read_readings,
    refuse_another_clock,
)
from myflopy.modflow.mf6.output_files import _declared_output_path
from myflopy.modflow.mf6.package_explorer_utils import cell_numbering, has_cell_numbering
from myflopy.modflow.mf6.period_end import period_lengths
from myflopy.modflow.mf6.target_values import read_target_values, target_value_table
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


def _usable(domain, layer, cell: int) -> bool:
    """Whether a target may land on ``cell``: active in its layer, or -- a screened
    target with no layer -- active in ANY layer (its screen picks the layers later,
    ledger 556)."""

    if not _is_blank(layer):
        return domain_active(domain, int(layer), cell)
    if domain is None:
        return True
    values = np.asarray(domain)
    if values.ndim == 1:
        return bool(values[cell] > 0)
    return bool((values[:, cell] > 0).any())


def _layer_key(layers: pd.Series) -> np.ndarray:
    """A location's layer as a join key: ``-1`` for a screened target with none."""

    return np.array([-1 if _is_blank(value) else int(value) for value in layers], dtype=int)


def _records_of_table(frame: pd.DataFrame) -> list[tuple[int, np.ndarray]]:
    """``(per, (nlay, ncells))`` records from a ``per``/``layer``/``cell``/``sim_head`` table.

    How a test double holding only a value table reaches the one reader
    (:func:`~myflopy.modflow.mf6.target_values.target_value_table`); a cell the
    table lacks reads ``NaN``.
    """

    if frame.empty:
        return []
    nlay = int(frame["layer"].max()) + 1
    ncells = int(frame["cell"].max()) + 1
    records = []
    for per, part in frame.groupby("per", sort=True):
        values = np.full((nlay, ncells), np.nan)
        values[part["layer"].astype(int).to_numpy(), part["cell"].astype(int).to_numpy()] = (
            part["sim_head"].astype(float).to_numpy()
        )
        records.append((int(per), values))
    return records


def _layer_series(layers: pd.Series) -> pd.Series:
    """Layers as ``int``, or ``Int64`` with ``<NA>`` when a screened target gives none."""

    if any(_is_blank(value) for value in layers):
        return pd.array([pd.NA if _is_blank(value) else int(value) for value in layers], dtype="Int64")
    return layers.astype(int)


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
        ``group_column``, ``weight_column``, ``screen_top_column`` and
        ``screen_bottom_column`` columns (and the geometry) are kept. Raises ``ValueError`` when ``name_column`` is missing, when a
        table with no geometry has no ``cell``/``node`` column, when a cell
        value is not numeric, when a present ``layer_column`` holds anything
        but a whole number from 0, or when a name repeats within one layer
        (see ``layer_column``).
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
        A LONG or one-series table may carry a ``weight`` column (matched
        ignoring case): each READING's own PEST weight (ledger 504). It
        beats the location's weight (``weight_column``), which beats 1.0; a blank means "none of its own", and
        anything but a blank or a finite number from 0 raises ``ValueError``.
        PESTPP-IES draws a reading's noise with standard deviation 1/weight. A
        WIDE table cannot carry one (its ``weight`` column would be a target of
        that name) and raises, naming the LONG form. Any other extra column of
        a LONG table is not read.
    name_column : str, optional
        Column of ``locations`` holding each target's name. Default ``"name"``.
        Required: ``ValueError`` when absent. For a ``Series`` or a
        ``{name: cell}`` dict it is the name given to the index or the keys.
        It applies to ``locations`` only -- a LONG ``values`` table always names
        its targets in a column called ``name``.
    layer_column : str or None, optional
        Column of ``locations`` holding each target's zero-based model layer.
        Default ``"layer"``. ``None``, or a column that is not there, puts
        every target in layer 0. When the column IS there every row must give a
        whole number from 0 (``1.0`` is layer 1): a blank, a non-number
        (``"2,3"``, ``"L3"``), a fraction or a negative raises ``ValueError``
        naming the rows (ledger 501 -- a blank or non-number was read
        as 0, and a fraction truncated) -- except on a SCREENED row (see
        ``screen_top_column``), whose layer may be blank: its screen says which
        layers it reads, and a layer it does give must be one its screen opens
        into. A name may appear once per layer --
        twice in ONE layer raises, naming it (ledger 502); in two
        layers, with no screens, it is two observations of the one measured
        series, each at full weight (ledger 234) -- give the well a screen to
        make it one. A target only matches a cell that is ACTIVE in its own layer
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
        prefix. The build refuses a label that is not text (a number: give the
        column as ``astype("string")``, which keeps a blank blank), is ``nan``,
        ``None`` or ``<NA>`` (what ``astype(str)`` writes for a missing label,
        ledger 359), holds no letter or digit, starts as PESTPP-IES reads an
        inequality (``l_``, ``less``, ``g_``, ``greater``), takes a group the
        build gives a set (``hds``, ``lak_stage_unmeasured``, ...), or becomes the
        same name as a different label (``"North"`` and ``"north"``) -- naming
        each.
    weight_column : str or None, optional
        Column of ``locations`` holding a per-location observation weight,
        carried as ``weight``. When ``cal.observe(...)`` builds the PEST control
        file, a weight that is not ``NaN`` becomes the weight of each
        observation at that location that has a target value and no weight of
        its own in ``values`` (a reading's weight wins, ledger 504), and
        ``NaN`` gives those observations weight 1; a period with no target value
        carries weight 0 (ledger 340). Default ``"weight"``; ``None`` or an
        absent column leaves it ``NaN``.
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``. Required once ``times`` has been applied
        (``ValueError`` otherwise); :meth:`get` calls it ``time``.
        :meth:`compare` and the PEST build read it as a ZERO-BASED STRESS
        PERIOD, by number (``1.0`` is period 1), and both raise on a label that
        is not a whole number, is negative, or is a date or a duration (pandas
        reads a datetime as nanoseconds; ledger 499 -- dated readings are binned
        onto periods by :meth:`from_hydrographs`), on two values for
        one target in one period, and on a period the run did not save (ledger
        500).
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
    screen_top_column : str or None, optional
        Column of ``locations`` holding each well's SCREEN TOP elevation, in the
        model's length units. Default ``"screen_top"``; ``None``, or a table
        without the column (and without ``screen_bottom_column``), reads no
        screens -- every target is a layer target, as before. A row giving both
        a top and a bottom is a SCREENED target (ledger 556):
        ONE observation, under its bare name, of the heads of every ACTIVE
        layer its screen opens into at its cell, each weighted by the layer's
        horizontal K (NPF ``K``) times the length of the screen inside the layer
        that is saturated -- MF6's ICELLTYPE rule (a confined layer conducts
        through its whole open length, a convertible one through the part below
        the head; a THICKSTRT layer's length is fixed from the starting heads).
        When every open interval is dry it reads the DEEPEST open layer's head.
        :meth:`compare` weights with the model's K; the PEST forward run with
        each realization's K, read from the files MF6 reads. A row leaving both
        blank is a layer target. Refused (``ValueError``, naming the rows): one
        of the two blank, a value that is not a finite number, a top at or below
        its bottom, the column without its twin, and a screened name on more
        than one row. Refused when the model is read (:meth:`observation_labels`):
        a layer the screen does not open into, a screen opening into no ACTIVE
        layer at its cell, and a DISU model; a screen reaching above the model
        top, below its bottom or into an inactive layer is clipped to what it
        opens into, with a warning. On ``ConcTargets`` screens are refused.
    screen_bottom_column : str or None, optional
        Column of ``locations`` holding each well's SCREEN BOTTOM elevation, in
        the model's length units; the twin of ``screen_top_column``. Default
        ``"screen_bottom"``.

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
    #: The MF6 OBS type of the family's quantity -- the only one ``to_flopy_obs``
    #: writes (ledger 477, as ledger 344 made it for LAK/SFR/DRN).
    _OBS_KIND: ClassVar[str] = "HEAD"
    #: The :class:`~myflopy.modflow.mf6.observations.hydrographs.Binning` of a
    #: :meth:`from_hydrographs` set (its report and clock stamp); ``None`` otherwise.
    _binned: ClassVar[Any] = None
    #: Whether a location may carry a well screen (ledger 556). A
    #: flux-weighted concentration needs the flow model's heads and K beside the
    #: transport file, so ``ConcTargets`` refuses screens.
    _SCREENS: ClassVar[bool] = True

    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple
    name_column: str = "name"
    layer_column: str | None = "layer"
    group_column: str | None = "group"
    weight_column: str | None = "weight"
    time_column: str = "time"
    value_column: str = "head"
    times: list | tuple | pd.Index | pd.Series | None = None
    screen_top_column: str | None = "screen_top"
    screen_bottom_column: str | None = "screen_bottom"

    def __post_init__(self):
        """Load and normalize locations and target tables."""

        self._locations = _load_locations(
            self.locations,
            name_column=self.name_column,
            layer_column=self.layer_column,
            group_column=self.group_column,
            weight_column=self.weight_column,
            screen_top_column=self.screen_top_column,
            screen_bottom_column=self.screen_bottom_column,
        )
        if not self._SCREENS and self.screened:
            raise ValueError(
                f"{type(self).__name__}: the locations carry well screens "
                f"({self.screen_top_column!r}/{self.screen_bottom_column!r}), and a "
                f"screened {self._FIELD_LABEL} is not read: a flux-weighted value needs "
                "the flow model's heads and K beside the transport file, which is "
                "deferred (ledger 556). Give each sampling point its "
                "layer, or pass screen_top_column=None, screen_bottom_column=None to "
                "ignore the columns."
            )
        self._screen_warned: set[str] = set()
        _refuse_reserved_target_names(
            self._locations["name"], owner=type(self).__name__, columns=("head_target",)
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
        layers=None,
        groups=None,
        weights=None,
        time_column: str = "time",
        value_column: str = "head",
        screen_tops=None,
        screen_bottoms=None,
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
            Zero-based model layer of each target. Default ``None``: layer 0 for
            a target with no screen, and no layer for a screened one (its screen
            says which layers it reads). A scalar applies to every cell; a list,
            tuple, array or ``Series`` must have one entry per cell
            (``ValueError``). Each must be a whole number from 0: a blank, a
            fraction or a negative raises (ledger 501), except a blank on a
            screened target; a layer given for a screened target must be one its
            screen opens into.
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
        screen_tops : float, sequence of float or None, optional
            Each target's screen-top elevation, in the model's length units
            (the constructor's ``screen_top_column``, ledger 556): a
            scalar for every cell, or one entry per cell, ``None``/``NaN`` for a
            target with no screen. Default ``None``: no screens. Given without
            ``screen_bottoms`` it raises ``ValueError``.
        screen_bottoms : float, sequence of float or None, optional
            Each target's screen-bottom elevation, the twin of ``screen_tops``.
            Default ``None``.

        Returns
        -------
        HeadTargets
            An instance of the class it is called on (``ConcTargets`` included).
        """

        cell_list = [int(value) for value in list(cells)]
        count = len(cell_list)
        if count == 0:
            raise ValueError("from_cells(...) requires at least one cell.")
        if (screen_tops is None) != (screen_bottoms is None):
            raise ValueError(
                "from_cells(...): screen_tops and screen_bottoms come together -- a "
                "screen is an interval (ledger 556)."
            )
        names_list = _coerce_label_list(names, count=count, default_prefix="OBS")
        screens = {}
        if screen_tops is not None:
            screens = {
                "screen_top": list(_coerce_scalar_or_list(screen_tops, count=count, default=np.nan)),
                "screen_bottom": list(
                    _coerce_scalar_or_list(screen_bottoms, count=count, default=np.nan)
                ),
            }
        # Checked where every location table is (``_location_layers``): ``int()``
        # here truncated a fractional layer in silence (ledger 501). With no layer
        # given, a screened target has none (its screen decides) and any other is
        # in layer 0.
        if layers is None:
            layer_list = [
                pd.NA if screens and not _is_blank(top) else 0
                for top in screens.get("screen_top", [np.nan] * count)
            ]
        else:
            layer_list = list(_coerce_scalar_or_list(layers, count=count, default=0))
        group_list = _coerce_scalar_or_list(groups, count=count, default=pd.NA)
        weight_list = _coerce_scalar_or_list(weights, count=count, default=np.nan)
        locations = pd.DataFrame(
            {
                "name": names_list,
                "layer": layer_list,
                "cell": cell_list,
                "group": group_list,
                "weight": weight_list,
                **screens,
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
            screen_top_column="screen_top",
            screen_bottom_column="screen_bottom",
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
            Optional: ``layer`` (zero-based, default 0), ``group``,
            ``screen_top``/``screen_bottom`` (a screened well, read as the
            constructor's ``screen_top_column`` reads them; ledger
            556), and
            ``weight`` -- each RECORD's own PEST weight, a blank meaning 1.0
            (ledger 504; it was the name's first record's, given to
            every reading). A target's location, screen and group are taken from
            its FIRST row. A ``geometry``
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
            for column in [
                "name", "layer", "cell", "node", "group", "screen_top", "screen_bottom",
                "geometry", "x", "y",
            ]
            if column in frame.columns
        ]
        locations = frame.loc[:, locations_columns].drop_duplicates(subset=["name"]).copy()
        if "geometry" in locations.columns:
            locations = gpd.GeoDataFrame(locations, geometry="geometry")
        # A record's weight is ITS reading's (ledger 504): it was taken from
        # each name's first record and given to every reading of it.
        values_frame = frame.loc[
            :, [time_column, "name", value_column, *(["weight"] if "weight" in frame.columns else [])]
        ].copy()
        return cls(
            locations=locations,
            values=values_frame,
            name_column="name",
            layer_column="layer",
            group_column="group",
            weight_column="weight",
            time_column=time_column,
            value_column=value_column,
            screen_top_column="screen_top",
            screen_bottom_column="screen_bottom",
        )

    @classmethod
    def from_hydrographs(
        cls,
        locations,
        readings,
        *,
        model,
        mode: str = "history",
        rule: str | None = None,
        date_column: str = "date",
        value_column: str | None = None,
        name_column: str = "name",
        layer_column: str | None = "layer",
        group_column: str | None = "group",
        weight_column: str | None = "weight",
        max_offset=None,
        between=None,
        min_years: int | None = None,
        screen_top_column: str | None = "screen_top",
        screen_bottom_column: str | None = "screen_bottom",
    ):
        """Targets from DATED readings, binned onto ``model``'s stress periods by one stated rule.

        A target set is keyed by stress period (ledger 340) and a date is refused as
        a key (ledger 499); this is where dated readings become period keys --
        against ``model``'s own calendar (``model.calendar``, its TDIS
        ``START_DATE_TIME``, ``TIME_UNITS`` and period lengths; ledger 529) --
        by the rule :mod:`~myflopy.modflow.mf6.observations.hydrographs` states
        once for heads, concentrations and lake stages (ledger 543). The
        result is an ordinary period-keyed set: :meth:`compare`, the PEST build,
        groups and weights read it unchanged. It also carries :attr:`binning` --
        which reading or statistic each target is -- and :attr:`calendar`, the
        CLOCK STAMP: :meth:`compare`, :meth:`simulated_heads` and the PEST build
        refuse the set on a model whose periods are other dates (another cycle
        count, another start), naming this call.

        Parameters
        ----------
        locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, or list or tuple of dict
            Where each well is, read exactly as the constructor's ``locations`` (a
            vector file, points, or zero-based cells within a layer), with
            ``name_column``, ``layer_column``, ``group_column`` and
            ``weight_column`` as there.
        readings : str, pathlib.Path, pandas.DataFrame or pandas.Series
            The dated readings, in the model's length units (concentration units
            on ``ConcTargets``). A path is read with ``pandas.read_csv``. WIDE:
            ``date_column`` plus one column per well, headed by its name -- the
            usual hydrograph export; LONG: ``date_column``, ``name`` and
            ``value_column``, one row per reading; ONE SERIES (``date_column`` and
            ``value_column`` only, or a ``Series`` indexed by date) when there is one
            location. Dates may be the table's index instead of a column. A blank
            value is no reading; any other entry that is not a finite number raises
            ``ValueError``. A name not among ``locations`` (stripped, ignoring case)
            raises ``ValueError`` naming it. Two readings of one well at one instant
            collapse when equal and raise when they differ. A ``weight`` column
            raises: several readings become one target, so weight a well through
            ``weight_column``.
        model : SimulationBase
            The model the targets will be compared on -- ``run.model(name)`` or
            ``mf.load_mf6_run(ws)`` -- whose calendar the readings are binned
            against. Its TDIS must declare ``START_DATE_TIME`` and a
            ``TIME_UNITS`` (``ValueError`` naming ``mf.tdis(...)`` otherwise).
        mode : {"history", "average_year"}, default "history"
            ``"history"``: the model runs the dated record (Ten Trails' monthly
            2017-2023 periods); readings outside its span,
            ``(START_DATE_TIME, last period end]``, are left out with a warning
            naming each well and its count (``warnings.warn`` and the log, ledger
            288), and none in span raises, printing both spans. ``"average_year"``:
            the model is an ``mf.AverageYear`` cycle (ledger 530) -- checked by
            ``mf.AverageYear.from_model(model)``, and a model that is not one raises
            naming ``mf.AverageYear`` -- and each calendar month's target, from
            every year's readings of that month, lands on the month's period of
            the LAST cycle (``AverageYear.months_of(-1)``).
        rule : str or None, default None
            How a period's value is made. History: ``"nearest_end"`` (the default)
            -- every reading goes to the period END nearest it (a reading midway
            between two goes to the period holding it), and each period's target is
            the reading nearest its end (a tie goes to the EARLIER reading, the one
            inside the period). It matches what is simulated: the head at the
            period's end (ledger 336). With Ten Trails' TDIS starting 2017-09-30
            every period ends at 00:00 on its month's LAST day, so a reading on
            31 October at noon is period 0's, where a rule by containment would
            give it to period 1. ``"median"``/``"mean"``: the statistic of the
            readings inside the period, ``(start, end]``. Average year:
            ``"median"`` (the default) -- each year's monthly median, then the
            median of those over the years, so a year read daily counts once, as
            a year read monthly does (the user's rule, Q4); ``"mean"`` -- the same
            with means. ``"nearest_end"`` there raises ``ValueError`` naming
            ``rule='median'``.
        date_column : str, default "date"
            The readings' date column, matched exactly. Each entry must be a date
            with a four-digit year and no time zone (``"2018-01-15"``,
            ``"2018-01-15 08:30"``); anything else raises ``ValueError`` naming the
            rows. Text written ``d/m/y`` or ``m/d/y`` is read MONTH first, as
            pandas reads it, and a column written day first in some rows
            (``"13/02/2018"``) and read month first in others (``"01/02/2018"``)
            raises ``ValueError`` -- parse it with ``pandas.to_datetime(...,
            dayfirst=True)`` and pass the table (ledger 545).
        value_column : str or None, default None
            The value column of a LONG or one-series ``readings`` table, matched
            ignoring case, and the result's ``value_column``. ``None``: the class
            default, ``"head"`` (``"conc"`` on ``ConcTargets``). Unused by a WIDE
            table.
        name_column : str, default "name"
            Column of ``locations`` holding each well's name (the constructor's).
        layer_column : str or None, default "layer"
            Column of ``locations`` holding each well's zero-based layer (the
            constructor's).
        group_column : str or None, default "group"
            Column of ``locations`` holding each well's observation-group label --
            its PEST group (the constructor's).
        weight_column : str or None, default "weight"
            Column of ``locations`` holding each well's PEST weight (the
            constructor's); with no reading weight possible, it is every target's
            weight at that well.
        max_offset : str, pandas.Timedelta, datetime.timedelta, numpy.timedelta64 or None, default None
            History with ``rule="nearest_end"`` only: a reading farther than this
            from its period's end is not used, so a period with no reading that
            near is unmeasured (``"10D"``, ``"36h"``). A bare number raises
            (pandas reads it as nanoseconds), as does a duration that is not
            positive. Given in average-year mode, or with another rule, it raises.
            ``None``: no limit.
        between : tuple of two dates or None, default None
            Average year only: ``(first, last)``, the readings used, both included,
            either ``None`` for open (``("2005-10-01", None)``). A date is its
            midnight, so ``last="2005-12-31"`` leaves out a reading later that day.
            Given in history mode it raises. ``None``: every reading.
        min_years : int or None, default None
            Average year only: a month read in fewer distinct calendar years than
            this has no target (its :attr:`binning` row says ``measured`` False).
            ``None`` means 1; a whole number from 1. Given in history mode it
            raises.
        screen_top_column : str or None, default "screen_top"
            Column of ``locations`` holding each well's screen-top elevation (the
            constructor's): a screened well is ONE observation of the head its
            screen reads (ledger 556). ``None`` or an absent column:
            no screens.
        screen_bottom_column : str or None, default "screen_bottom"
            Column of ``locations`` holding each well's screen-bottom elevation (the
            constructor's).

        Returns
        -------
        HeadTargets
            An instance of the class it is called on, keyed by zero-based stress
            period, carrying :attr:`binning` and :attr:`calendar`.

        Raises
        ------
        ValueError
            As listed under each parameter; and as the constructor raises for the
            locations.

        Warns
        -----
        UserWarning
            History mode: readings outside the model's span, counted per well; any
            mode: a well whose readings all fell out, by name.

        Examples
        --------
        >>> flow = run.model("flow")
        >>> heads = mf.HeadTargets.from_hydrographs("wells.gpkg", "hydrographs.csv", model=flow)
        >>> heads.binning[["name", "per", "reading", "offset_days", "value"]].head()
        >>> year_heads = mf.HeadTargets.from_hydrographs(
        ...     "wells.gpkg", "hydrographs.csv", model=flow, mode="average_year",
        ...     between=("2005-10-01", None), min_years=3)
        """

        owner = cls.__name__
        if value_column is None:
            value_column = cls.__dataclass_fields__["value_column"].default
        if isinstance(locations, (str, Path)):
            locations = gpd.read_file(locations)
        known = _load_locations(
            locations,
            name_column=name_column,
            layer_column=layer_column,
            group_column=group_column,
            weight_column=weight_column,
            screen_top_column=screen_top_column,
            screen_bottom_column=screen_bottom_column,
        )
        calendar, months = model_clock(model, mode=mode, owner=owner)
        long = read_readings(readings, date_column=date_column, value_column=value_column, owner=owner)
        if "name" not in long.columns:
            if known["name_key"].nunique() != 1:
                raise ValueError(
                    f"{owner}.from_hydrographs(...): the readings name no well (one "
                    f"series), and there are {known['name_key'].nunique()} locations; give "
                    "a WIDE table, one column per well, or a LONG one with a 'name' column."
                )
            long["name"] = str(known["name"].iloc[0])
        unknown = sorted(
            set(_normalize_identifier_series(long["name"])) - set(known["name_key"])
        )
        if unknown:
            raise ValueError(
                f"{owner}.from_hydrographs(...): the readings name {len(unknown)} well(s) "
                f"not among the locations (compared stripped, ignoring case): "
                + ", ".join(unknown[:10]) + (" ..." if len(unknown) > 10 else "")
                + ". Add them to the locations, or drop their columns."
            )
        binning = bin_readings(
            long, calendar, mode=mode, rule=rule, max_offset=max_offset, between=between,
            min_years=min_years, months=months, owner=owner,
        )
        targets = cls(
            locations=locations,
            values=binning.values.rename(columns={"value": value_column}),
            name_column=name_column,
            layer_column=layer_column,
            group_column=group_column,
            weight_column=weight_column,
            time_column="time",
            value_column=value_column,
            screen_top_column=screen_top_column,
            screen_bottom_column=screen_bottom_column,
        )
        targets._binned = binning
        return targets

    @property
    def binning(self) -> pd.DataFrame | None:
        """How each target was made from dated readings, or ``None`` for a set keyed by period.

        Set by :meth:`from_hydrographs`: one row per (well, period) its readings
        reached -- ``name``, ``per``, ``month`` (the period's midpoint month),
        ``period_start``, ``period_end``, ``n_readings``, ``n_years``,
        ``first_reading``, ``last_reading``, ``reading`` (the reading used, under
        ``rule="nearest_end"``), ``offset_days`` (it minus the period end; negative
        before), ``value`` and ``measured`` (False for an average-year month under
        ``min_years``, which is no target).
        """

        binned = self._binned
        return None if binned is None else binned.report.copy()

    @property
    def calendar(self):
        """The :class:`~myflopy.modflow.mf6.period_calendar.PeriodCalendar` the targets were binned against, or ``None``.

        The CLOCK STAMP of a :meth:`from_hydrographs` set: comparing it with, or
        building PEST on, a model whose calendar differs raises (ledger
        543). ``None`` for a set keyed by period by hand.
        """

        binned = self._binned
        return None if binned is None else binned.calendar

    @property
    def screened(self) -> bool:
        """Whether any location carries a well screen (ledger 556).

        A screened location is ONE observation of the transmissivity-weighted head
        over the layers its screen opens into (see the constructor's
        ``screen_top_column``).
        """

        return "screen_top" in self._locations.columns and bool(
            self._locations["screen_top"].notna().any()
        )

    @property
    def locations_gdf(self) -> gpd.GeoDataFrame:
        """Normalized observation locations table."""

        return self._locations.copy()

    def get(self) -> pd.DataFrame:
        """Return the merged long-format target table.

        One row per values row and location row sharing its name (two for a well
        given in two layers; one for a SCREENED well, which is one location row,
        with its ``screen_top``/``screen_bottom``). ``weight`` is the reading's own weight where the
        values table gave one, else its location's (ledger 504) -- the
        weight the PEST build gives the observation, ``NaN`` meaning 1.0 there.

        Returns
        -------
        pandas.DataFrame
            The values (``time``, ``head_target``, ``name_key``) joined to the
            location columns (``name``, ``cell``, ``layer``, ``group``,
            ``weight``, ``x``, ``y``).
        """

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
        if "weight_loc" in merged.columns:
            # A reading's own weight beats its location's (ledger 504).
            merged["weight"] = merged["weight"].fillna(merged["weight_loc"])
            merged = merged.drop(columns=["weight_loc"])
        return _copy_frame(merged)

    def to_long(self) -> pd.DataFrame:
        """Return the long-format target table."""

        return self.get()

    def to_wide(self) -> pd.DataFrame:
        """Return the target table in wide format with one column per location.

        One row per stress period a ``time`` names, and one column per target
        name holding its target value (a head, or a concentration on
        ``ConcTargets``); a set given no values gives ``time`` and one empty column
        per target. ``time`` is the zero-based period as ``int``, read the way
        :meth:`compare` and the PEST build read it, so ``"1"`` and ``1`` are one
        row (they were two, ledger 422). Two rows HOLDING a value for one
        target in one stress period are REFUSED -- names compared stripped and
        lower-cased -- by the check :meth:`compare` and the PEST build read
        (ledger 340): this used to keep the first in silence
        (``pivot_table(aggfunc="first")``), so the wide table and the build
        disagreed about the same set (ledger 371). A blank row is unmeasured and
        never collides; a target whose every row is blank keeps its column, all
        blank, as :meth:`compare` keeps its rows (ledger 422).

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

        _measured_periods(
            self._values,
            value_column="head_target",
            name_column="name_key",
            caller=f"{type(self).__name__}.to_wide()",
        )
        # After the check a (target, period) cell holds at most ONE value, so
        # "first" only prefers it over blank rows (and, for a well given in two
        # layers, over the copy ``get()`` makes for its second layer row).
        return _wide_by_period(
            self.get(), value_column="head_target", caller=f"{type(self).__name__}.to_wide()",
            names=self._locations["name"],
        )

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
            ``vor.gdf_vorPolys``, the grid's cell polygons, or a FloPy ``modelgrid``
            (a DIS model, ledger 420). Which cells are
            active comes from its discretization's ``idomain``, else its
            ``myflopy_context.domain``; with neither, every cell is active.

        Returns
        -------
        pandas.DataFrame
            One row per MATCHED target location, with its ``cell`` and ``layer``
            (``<NA>`` for a screened target given no layer), and its
            ``screen_top``/``screen_bottom`` when the set carries screens.

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

        A target's ``layer`` must be one the model's views number (ledger 419), or
        ``ValueError`` names it: on DISU that is layer ``0`` alone, because MF6
        numbers a DISU model's cells by node and writes its heads as one record of
        every node. A point there resolves to the lowest-numbered node holding it;
        name a deeper node by its ``cell`` (ledger 452).

        A SCREENED target with no layer (ledger 556) matches the
        lowest-numbered cell holding it that is active in ANY layer; which of its
        layers it reads is :meth:`observation_labels`' question. A screened set on
        a DISU model is refused here: a screen reads a column's layers, and DISU has
        no columns.
        """

        self._refuse_layers_off_the_grid(model)
        if "cell" in self._locations.columns:
            matches = self._validated_cells(model)
            if matches["x"].isna().any() or matches["y"].isna().any():
                vor = model_cell_polygons(model).reset_index().rename(columns={"index": "cell"})
                centroids = pd.DataFrame(
                    {
                        "cell": vor["cell"].astype(int),
                        "x": vor.geometry.centroid.x,
                        "y": vor.geometry.centroid.y,
                    }
                )
                matches = matches.drop(columns=["x", "y"], errors="ignore").merge(centroids, on="cell", how="left")
            return _copy_frame(
                matches.loc[:, ["name", "name_key", "layer", "group", "weight", "cell", "x", "y",
                                *self._screen_columns()]]
            )

        return self._matched_from_geometry(model)

    def _screen_columns(self) -> list[str]:
        """``["screen_top", "screen_bottom"]`` when the locations carry screens, else none."""

        return [column for column in ("screen_top", "screen_bottom") if column in self._locations.columns]

    def _refuse_layers_off_the_grid(self, model) -> None:
        """Refuse a target layer the model's views do not number (ledger 419).

        A DISU target in layer 1 matched the node holding its point and was written
        as layer 1, and the PEST build died reading the head file's one record with
        ``IndexError: index 1 is out of bounds`` (ledger 452). A model whose
        grid has no cellid form (a test double, no discretization) is not checked.
        """

        if not has_cell_numbering(model):
            return
        numbering = cell_numbering(model)
        if numbering.kind == "unstructured" and self.screened:
            raise ValueError(
                f"{self._FIELD_LABEL.capitalize()} targets: a screened well reads the "
                "layers of its COLUMN, and a DISU model has no columns -- its cells are "
                "numbered by node, in one layer. Give each screened well's node by its "
                "`cell` with no screen (ledger 556)."
            )
        try:
            numbering.check(layers=self._locations["layer"].dropna().astype(int))
        except ValueError as error:
            hint = (
                " A DISU model's cells are numbered by node, in one layer: give a "
                "target layer 0, and a node below the top one by its `cell`, since a "
                "point resolves to the lowest-numbered node holding it."
                if numbering.kind == "unstructured" else ""
            )
            raise ValueError(f"{self._FIELD_LABEL.capitalize()} targets: {error}{hint}") from error

    def _validated_cells(self, model) -> pd.DataFrame:
        """The locations, with their explicit cell ids checked against the grid.

        An id is an assertion about the model, so a wrong one raises rather than
        being warned about the way a missed point is. Before this, an id the grid
        did not have reached the value table as NaN and dropped out of ``stats``
        silently, and an inactive id was read where MF6 writes no value.
        """

        matches = self._locations.copy()
        ncpl = len(model_cell_polygons(model))
        domain = model_domain(model, ncpl)
        cells = matches["cell"].astype(int)
        layers = matches["layer"]
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
            if not _usable(domain, layer, cell)
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

        polys = model_cell_polygons(model)
        locations = self._locations_in_grid_crs(model)
        domain = model_domain(model, len(polys))
        columns = ["name", "name_key", "layer", "group", "weight", "cell", "x", "y",
                   *self._screen_columns()]
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
            layer = row["layer"]
            hits = cell_shares(
                polys,
                geometry,
                # Eligibility only. ``cell_shares`` sorts whatever this returns,
                # so the lowest-numbered usable cell wins a tie without a second
                # copy of that rule here.
                usable=lambda cells, layer=layer: [
                    cell for cell in cells if _usable(domain, layer, cell)
                ],
            )
            if not hits:
                reason = miss_reason(
                    polys, geometry,
                    target="active cell" if _is_blank(layer) else f"active cell in layer {int(layer)}",
                )
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
        frame["layer"] = _layer_series(frame["layer"])
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
            :meth:`match_to_model`'s table plus ``obsname``, the label. A SCREENED
            location (ledger 556) is one label -- its bare name --
            on one row per ACTIVE layer its screen opens into at its cell, each
            with ``layer`` that layer, ``open_bottom``/``open_top`` the screen
            clipped to the cell, ``icelltype`` the layer's NPF cell type there and
            ``fixed_b`` a THICKSTRT layer's open length from the starting heads
            (``NaN`` otherwise); those four columns are present only for a set with
            screens, ``NaN`` on its unscreened rows.

        Raises
        ------
        ValueError
            As :meth:`match_to_model`; and for a screened location, when its given
            layer is not one its screen opens into, when it opens into no active
            layer at its cell, or when the model has no NPF to weight by.

        Warns
        -----
        UserWarning
            Once per screen reaching above the model top, below its bottom or
            into an inactive layer: it is clipped to what it opens into.
        """

        return self._expand_screens(self._labelled(model), model)

    def _labelled(self, model) -> pd.DataFrame:
        """:meth:`match_to_model` plus ``obsname``: one row per matched LOCATION row."""

        matches = self.match_to_model(model).copy()
        matches["obsname"] = _ensure_unique_observation_names(matches)
        return matches

    def _expand_screens(self, labelled: pd.DataFrame, model) -> pd.DataFrame:
        """A screened location as one row per ACTIVE layer its screen opens into (ledger 556).

        Each row carries the screen clipped to that layer's cell (``open_bottom``,
        ``open_top``), the layer's NPF ``icelltype`` there and, for a THICKSTRT
        layer, the open length fixed from the starting heads (``fixed_b``; ``NaN``
        otherwise) -- what :func:`~myflopy.modflow.mf6.hydraulics.screened_head`
        weights by. An unscreened location is passed through, its screen columns
        ``NaN``.
        """

        if not self.screened:
            return labelled
        screened = labelled["screen_top"].notna().to_numpy()
        try:
            flow = model_flow_properties(model)
        except ValueError as error:
            raise ValueError(
                f"{type(self).__name__}: a screened well's head is weighted by each "
                "layer's K and saturated open length, read off the model's NPF and "
                f"discretization, and this model cannot answer: {error}"
            ) from error
        botm = np.asarray(flow.botm, dtype=float)
        nlay = botm.shape[0]
        botm = botm.reshape(nlay, -1)
        ncells = botm.shape[1]
        top = np.asarray(flow.top, dtype=float).reshape(-1)
        layer_tops = np.vstack([top.reshape(1, -1), botm[:-1]])
        celltype = (
            np.ones((nlay, ncells), dtype=int) if flow.icelltype is None
            else np.broadcast_to(np.asarray(flow.icelltype).reshape(nlay, -1), (nlay, ncells)).astype(int)
        )
        start = (
            None if not flow.thickstrt or flow.strt is None
            else np.broadcast_to(np.asarray(flow.strt, dtype=float).reshape(nlay, -1), (nlay, ncells))
        )
        domain = model_domain(model, ncells)
        rows: list[dict] = []
        notes: list[str] = []
        for index, row in enumerate(labelled.to_dict("records")):
            if not screened[index]:
                rows.append(row | {"open_bottom": np.nan, "open_top": np.nan,
                                   "icelltype": np.nan, "fixed_b": np.nan})
                continue
            cell = int(row["cell"])
            name = str(row["name"])
            screen_top, screen_bottom = float(row["screen_top"]), float(row["screen_bottom"])
            open_top = np.minimum(screen_top, layer_tops[:, cell])
            open_bottom = np.maximum(screen_bottom, botm[:, cell])
            opens = np.flatnonzero(open_top > open_bottom)
            given = row["layer"]
            if not _is_blank(given) and int(given) not in set(opens.tolist()):
                raise ValueError(
                    f"{self._FIELD_LABEL.capitalize()} target {name!r} gives layer "
                    f"{int(given)} and a screen from {screen_bottom:g} to {screen_top:g}, "
                    f"which at its cell {cell} opens into "
                    + (f"layer(s) {opens.tolist()}" if len(opens) else "no layer")
                    + f" (the column runs from {top[cell]:g} down to {botm[-1, cell]:g}). "
                    "Correct the layer or the screen, or leave the layer blank to let the "
                    "screen decide (ledger 556)."
                )
            active = [layer for layer in opens if domain_active(domain, int(layer), cell)]
            if not active:
                where = (
                    "lies entirely above the model top" if screen_bottom >= top[cell]
                    else "lies entirely below the model bottom" if screen_top <= botm[-1, cell]
                    else f"opens only into layer(s) {opens.tolist()}, inactive there"
                )
                raise ValueError(
                    f"{self._FIELD_LABEL.capitalize()} target {name!r}: its screen, from "
                    f"{screen_bottom:g} to {screen_top:g}, {where} at its cell {cell} "
                    f"(the column runs from {top[cell]:g} down to {botm[-1, cell]:g}), so "
                    "MF6 writes no head it reads. Check the screen's elevations (the "
                    "model's length units and datum) and the well's location (ledger "
                    "341's rule: a target PEST could never see is refused)."
                )
            lost = []
            if screen_top > top[cell]:
                lost.append(f"{screen_top - top[cell]:g} above the model top")
            if screen_bottom < botm[-1, cell]:
                lost.append(f"{botm[-1, cell] - screen_bottom:g} below the model bottom")
            idle = [int(layer) for layer in opens if int(layer) not in set(int(a) for a in active)]
            if idle:
                lost.append(f"its part in inactive layer(s) {idle}")
            if lost:
                notes.append(f"{name!r} ({', '.join(lost)})")
            for layer in active:
                layer = int(layer)
                fixed = np.nan
                if start is not None and celltype[layer, cell] < 0:
                    fixed = float(saturated_thickness(
                        start[layer, cell], open_top[layer], open_bottom[layer],
                        icelltype=celltype[layer, cell], start=start[layer, cell],
                    ))
                rows.append(row | {
                    "layer": layer,
                    "open_bottom": float(open_bottom[layer]),
                    "open_top": float(open_top[layer]),
                    "icelltype": int(celltype[layer, cell]),
                    "fixed_b": fixed,
                })
        new_notes = [note for note in notes if note not in self._screen_warned]
        if new_notes:
            self._screen_warned.update(new_notes)
            message = (
                f"{self._FIELD_LABEL.capitalize()} targets: {len(new_notes)} screen(s) reach "
                f"past what they read, and are clipped to the active layers they open into: "
                f"{_listing(new_notes)} (ledger 556)."
            )
            logger.warning(message)
            warnings.warn(message, UserWarning, stacklevel=3)
        frame = pd.DataFrame(rows)
        frame["layer"] = frame["layer"].astype(int)
        frame["cell"] = frame["cell"].astype(int)
        return _copy_frame(frame)

    def simulated_heads(self, model) -> pd.DataFrame:
        """Return simulated heads for the target locations by stress period.

        The result is a wide table suitable for pyEMU list-style observation
        setup: one row per stress period and one column per observation LABEL
        (see :meth:`observation_labels`, which a repeated name disambiguates), in
        the order the forward run writes them. Each period's value is the one at
        its END (ledger 336), read on the head file's own clock -- the file the
        model's OC declares (``HEAD FILEOUT``, ledger 505), through the one
        period-end reader PEST's forward run uses (ledger 557). A
        SCREENED well's column is its layers' heads weighted by NPF ``K`` x
        saturated open length (ledger 556).

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

        refuse_another_clock(self, model, caller=f"{type(self).__name__}.simulated_heads(...)")
        # The label, not the bare name: one column per (name, layer), so a
        # two-layer well is two observations rather than one column whose value
        # depended on which of its rows got there first -- and one column per
        # SCREENED well, its layers weighted into one head (ledger 556).
        return self._simulated_table(model, self.observation_labels(model))

    def _simulated_table(self, model, mapping: pd.DataFrame) -> pd.DataFrame:
        """``per`` plus one column per label of ``mapping``: the ONE period-end reader (ledger 557).

        A run model's values come from the file its OC DECLARES (``HEAD FILEOUT``,
        ``CONCENTRATION FILEOUT``; ledger 505) through
        :func:`~myflopy.modflow.mf6.target_values.read_target_values`, the reader
        PEST's forward run uses, at each period's end on the file's own clock
        (ledger 336). It read ``<model>.hds`` through a ``HeadsPlus`` built twice
        per call, with the whole value table (7.0-9.5 s and 1.9 GB per
        ``compare()`` on Ten Trails at 38 periods). A duck-typed double with no
        workspace or TDIS is read from its value table, through the same
        :func:`~myflopy.modflow.mf6.target_values.target_value_table`.
        """

        k = None
        if "open_top" in mapping.columns and mapping["open_top"].notna().any():
            flow = model_flow_properties(model)
            k = np.asarray(flow.k, dtype=float).reshape(np.asarray(flow.botm).shape[0], -1)
        kind, saves, text = {
            "all_heads": ("head", "HEAD", "head"),
            "all_conc": ("concentration", "CONCENTRATION", "concentration"),
        }[self._TABLE_ATTRIBUTE]
        workspace = getattr(model, "workspace", None)
        perlen = period_lengths(model)
        if workspace is None or perlen is None:
            frame = _simulated_field_by_period(
                model, table_attribute=self._TABLE_ATTRIBUTE, store_column=self._STORE_COLUMN,
            )
            return target_value_table(_records_of_table(frame), mapping, k=k)
        gate = getattr(model, "_require_model_kind", None)
        if callable(gate) and self._TABLE_ATTRIBUTE == "all_conc":
            gate("gwt6", "all_conc", "concentration")
        path = _declared_output_path(Path(workspace), model.name, kind)
        try:
            shown = path.relative_to(Path(workspace)).as_posix()
        except ValueError:
            shown = str(path)
        return read_target_values(
            path, mapping, perlen=perlen, file=shown, saves=saves, text=text, k=k,
        )

    def compare(self, model) -> pd.DataFrame:
        """Compare target heads against simulated heads by stress period.

        A target's ``time`` is read as a ZERO-BASED STRESS PERIOD -- the one
        reading every target family and PEST's build share (ledger 340): ``1.0``
        is period 1, and a key with a fraction is refused rather than truncated
        (``10.9`` used to read period 10). ``sim_head`` is the head at each
        period's END (ledger 336), as :meth:`simulated_heads` reads it -- for a
        SCREENED well the K x saturated-open-length weighted head of its open
        layers (ledger 556), weighted with this model's NPF ``K``.
        A blank target value is a row with ``NaN`` ``head_target`` and
        ``residual``. Every target is matched to the model first
        (:meth:`match_to_model`), so a target given as a cell id the grid does not
        have, or one inactive in its layer, raises here too.

        Parameters
        ----------
        model : SimulationBase
            A run model holding the head (or concentration) file.

        Returns
        -------
        pandas.DataFrame
            One row per target row: the target's columns plus ``per``,
            ``cell``, ``sim_head``, ``residual`` (simulated minus target) and
            ``abs_residual``. ``x``/``y`` are the matched cell's or point's, in
            the GRID's frame, for a target given as a point and one given as a
            cell alike (a point target returned ``x_x``/``y_x``, its own CRS,
            beside ``x_y``/``y_y``, and no ``x``).

        Raises
        ------
        ValueError
            When a target's ``time`` is a date or a duration, not a number, not a
            whole number or negative; when two rows hold a value for one target in
            one period (PEST keeps one observation per target per period); when a
            target names a period the file did not save -- past the run's last
            period, say (ledger 500: it read ``NaN`` and ``stats()`` left
            it out in silence); or when a period the file saved lacks its final
            step (as :meth:`simulated_heads`).
        """

        caller = f"{type(self).__name__}.compare(...)"
        refuse_another_clock(self, model, caller=caller)
        _measured_periods(
            self._values, value_column="head_target", name_column="name_key", caller=caller
        )
        targets = self.get().copy()
        targets["per"] = _numeric_periods(targets, caller=caller)
        labelled = self._labelled(model)
        # The ONE period-end reader, shared with PEST's forward run (ledger
        # 557): one value per label, a screened well's layers weighted
        # into one (ledger 556).
        table = self._simulated_table(model, self._expand_screens(labelled, model))
        # A key past the run's periods matched no row, and stats() dropped it in
        # silence while the PEST build refused it (ledger 500).
        _refuse_keys_off_the_periods(
            targets, table["per"].unique(), model=model, caller=caller,
            route=_FIELD_FILES[self._TABLE_ATTRIBUTE][1],
        )
        sim = table.melt(id_vars="per", var_name="obsname", value_name="sim_head")
        # Each target row is its location row's label: by name and the location's
        # own layer (blank for a screened well with none, matched as one key).
        # The matched cell and x/y -- in the GRID's frame -- win over the targets'
        # own: both sides carried x/y, so pandas returned x_x/y_x beside x_y/y_y.
        matches = labelled.loc[:, ["name_key", "layer", "cell", "x", "y", "obsname"]].assign(
            _layer_key=_layer_key(labelled["layer"])
        ).drop(columns=["layer"])
        merged = (
            targets.drop(columns=["cell", "x", "y"], errors="ignore")
            .assign(_layer_key=_layer_key(targets["layer"]))
            .merge(matches, on=["name_key", "_layer_key"], how="left")
            .merge(sim, on=["per", "obsname"], how="left")
            .drop(columns=["_layer_key", "obsname"])
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
        kind: str | None = None,
    ) -> dict[str, list[tuple[str, str, tuple[int, ...]]]]:
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
        kind : str or None, optional
            The MF6 observation type written on every record. Default ``None``:
            the family's own -- ``"HEAD"``, or ``"CONCENTRATION"`` on
            ``ConcTargets``. A string must be that type (any case); any other
            raises ``ValueError`` -- it was written unchecked, and the default
            gave a ``ConcTargets`` ``"HEAD"``, which a GWT model's OBS file
            refuses, until ledger 477.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(obsname, kind, cellid), ...]}``, the cellid zero-based
            (FloPy writes it one-based) in the model's own form --
            ``(layer, row, col)`` on DIS, ``(layer, cell)`` on DISV, ``(node,)``
            on DISU (``CellNumbering.cellid``, ledger 455) -- and
            ``(layer, cell)`` when ``model`` has no discretization to read;
            ``obsname`` the label from :meth:`observation_labels`. The shape
            ``ModflowUtlobs(continuous=)`` and a package's ``observations=`` take.
        """

        kind = self._observation_kind(kind)
        if self.screened:
            raise ValueError(
                f"{type(self).__name__}.to_flopy_obs(...): these targets hold screened "
                "wells, and an MF6 HEAD observation reads ONE cell -- MF6's OBS has no "
                "multi-layer head. Compare them with compare(), or calibrate them through "
                "cal.observe(...), which weight each screen's layers by K x saturated "
                "open length (ledger 556)."
            )
        matches = self.observation_labels(model)
        # The cellid the model's own OBS reads (ledger 455): `(layer, cell)` was
        # written on every grid, which MF6 refuses on DIS and reads on DISU as node
        # `layer + 1` -- a target at node 6 observed node 1, in silence.
        numbering = cell_numbering(model) if has_cell_numbering(model) else None
        records = [
            (
                str(row.obsname),
                kind,
                (int(row.layer), int(row.cell)) if numbering is None
                else numbering.cellid(row.layer, row.cell),
            )
            for row in matches.itertuples(index=False)
        ]
        return {str(csv_name): records}

    def _observation_kind(self, kind: str | None) -> str:
        """The family's MF6 OBS type; ``kind`` must be it or ``None`` (ledger 477)."""

        if kind is None:
            return self._OBS_KIND
        written = str(kind).strip().upper()
        if written != self._OBS_KIND:
            raise ValueError(
                f"{type(self).__name__} writes {self._OBS_KIND} observations, the "
                f"quantity it compares; kind={kind!r} would be declared under the "
                "targets' names as another quantity. Drop kind=, or declare another "
                "observation type with FloPy directly, under a name of its own."
            )
        return written

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str | None = None,
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
        kind : str or None, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`: default
            ``None``, the family's own (``"HEAD"``, ``"CONCENTRATION"`` on
            ``ConcTargets``); any other raises ``ValueError``.
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
              lines per location and layer it is observed in (it raised
              ``KeyError`` unless every location was observed in every layer
              until ledger 477);
            - ``"residuals_by_period"``: mean absolute residual per ``time``.

            Any other string raises ``ValueError``.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        return CalibrationPlot.from_targets(self, model=model, type=type)


class BoundHeadTargets:
    """Model-bound head-target helper returned from ``model.targets.heads``.

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from. A ``ConcTargets`` set
    registered on ``model.targets`` is bound by this class too (it is a
    ``HeadTargets`` subclass), so "head" below reads "concentration" for one.
    """

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
        """Plot the target locations on a map axis (delegates to the wrapped targets).

        The bound model takes no part: the points are drawn in the CRS their
        locations carry, not reprojected to the grid's.

        Parameters
        ----------
        ax : matplotlib.axes.Axes or None, optional
            Axes to draw on. Default ``None`` makes a new figure and axes.
        **kwargs
            Forwarded unchanged to ``geopandas.GeoDataFrame.plot`` (for example
            ``color=``, ``markersize=``); one it does not take raises there.

        Returns
        -------
        matplotlib.axes.Axes
            The axes, titled "Head Target Locations".

        Raises
        ------
        ValueError
            When the locations carry no geometry -- targets given as cell ids,
            as the canonical model's are.
        """

        return self.targets.plot_locations(ax=ax, **kwargs)

    def match_to_model(self, model=None) -> pd.DataFrame:
        """Match each target location to its grid cell (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The model to match against. Default ``None`` uses the bound model;
            any other model-like object exposing ``vor.gdf_vorPolys``, or a FloPy
            ``modelgrid`` (a DIS model, ledger 420), is used instead. Matching is :meth:`HeadTargets.match_to_model`'s: one cell
            per point, active in the target's own layer, a miss warned about by
            name and none landing raising ``ValueError``. A target given as a
            cell id (the canonical model's are) is checked instead: an id the
            grid does not have, or one inactive in its layer, raises
            ``ValueError``.

        Returns
        -------
        pandas.DataFrame
            One row per MATCHED target location, with its ``cell`` and ``layer``
            (``<NA>`` for a screened target given no layer), and its
            ``screen_top``/``screen_bottom`` when the set carries screens.
        """

        return self.targets.match_to_model(self.model if model is None else model)

    def simulated_heads(self, model=None) -> pd.DataFrame:
        """Simulated heads at the target locations (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model whose head (or concentration) file is read. Default
            ``None`` uses the bound model; another run reads that run's file for
            the same targets. Each period's value is the one at its END (ledger
            336), as :meth:`HeadTargets.simulated_heads` reads it, whose
            ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            ``per`` plus one column per observation label, in the field's units
            (model length units for a head).
        """

        return self.targets.simulated_heads(self.model if model is None else model)

    def compare(self, model=None) -> pd.DataFrame:
        """Compare observed vs simulated heads (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against. Default ``None`` uses the bound
            model; another run compares the same targets against it -- what the
            plots' ``baseline=`` does. Read as :meth:`HeadTargets.compare` reads
            it (each target's ``time`` a zero-based stress period, the head at
            its END), whose ``Raises`` apply.

        Returns
        -------
        pandas.DataFrame
            One row per target row: the target's columns plus ``per``, ``cell``,
            ``sim_head``, ``residual`` (simulated minus target) and
            ``abs_residual``.
        """

        return self.targets.compare(self.model if model is None else model)

    def residuals(self, model=None) -> pd.DataFrame:
        """The residual table for the head targets (defaults to the bound model).

        Parameters
        ----------
        model : SimulationBase or None, optional
            The run model to compare against, as in :meth:`compare`. Default
            ``None`` uses the bound model.

        Returns
        -------
        pandas.DataFrame
            :meth:`compare`'s rows, reduced as :meth:`HeadTargets.residuals`
            reduces them (``name`` ... ``residual``, ``abs_residual``,
            ``weight``).
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
            minus target, in the field's units (:meth:`HeadTargets.stats`).
        """

        return self.targets.stats(self.model if model is None else model)

    def to_flopy_obs(
        self,
        *,
        csv_name: str = "head_targets.csv",
        kind: str | None = None,
    ):
        """A FloPy ``continuous`` observation dict for the head targets on the bound model.

        The targets are matched against the BOUND model
        (:meth:`HeadTargets.observation_labels`): one record per matched target,
        a target on no active cell left out with a warning, none landing
        ``ValueError``.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key,
            passed through ``str``. Default ``"head_targets.csv"``.
        kind : str or None, optional
            The MF6 observation type written on every record, as in
            :meth:`HeadTargets.to_flopy_obs`: default ``None``, the family's own
            (``"HEAD"``, ``"CONCENTRATION"`` on a bound ``ConcTargets``); any
            other raises ``ValueError``.

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(obsname, kind, cellid), ...]}``, the cellid zero-based
            in the bound model's own form, as :meth:`HeadTargets.to_flopy_obs`
            returns it.
        """

        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "gwf_obs",
        filename: str = "head_targets.obs",
        csv_name: str | None = None,
        kind: str | None = None,
        package=None,
    ):
        """Attach an MF6 head observation package to the bound model.

        The targets are matched against the bound model (:meth:`to_flopy_obs`)
        and the package is added to its FloPy objects in memory; writing and
        running the simulation afterwards is the caller's.

        Parameters
        ----------
        pname : str, optional
            FloPy package name of the OBS package. Default ``"gwf_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, relative to the owning model's
            directory; its stem names the default ``csv_name``. Default
            ``"head_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"head_targets.csv"`` for the
            default ``filename``).
        kind : str or None, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`: default
            ``None``, the family's own; any other raises ``ValueError``.
        package : flopy object or None, optional
            The FloPy model or package to attach the OBS package to. Default
            ``None`` attaches it to the bound model's FloPy model, ``model.gwf``
            (on a transport model, the GWT model). The records are matched
            against the bound model either way.

        Returns
        -------
        flopy.mf6.modflow.mfutlobs.ModflowUtlobs
            The observation package.
        """

        return self.targets.attach_flopy_obs(
            self.model,
            pname=pname,
            filename=filename,
            csv_name=csv_name,
            kind=kind,
            package=package,
        )

    def calibration_plot(self, *, type: str = "calibration"):
        """A calibration plot for the head targets on the bound model.

        Parameters
        ----------
        type : str, optional
            Which plot, drawn from :meth:`compare` on the bound model (whose
            ``Raises`` apply):

            - ``"calibration"`` (default): target against simulated, one point
              per row, with fit statistics;
            - ``"obs_vs_sim"``: the same cross plot without the statistics;
            - ``"heads"``: target and simulated values through time, one pair of
              lines per location and layer it is observed in (it raised
              ``KeyError`` unless every location was observed in every layer
              until ledger 477);
            - ``"residuals_by_period"``: mean absolute residual per ``time``.

            Any other string raises ``ValueError``.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

        return self.targets.calibration_plot(self.model, type=type)

    @property
    def plot(self):
        """Plotting helpers for this model-bound target set."""

        if self._plot is None:
            self._plot = BoundHeadTargetPlots(self)
        return self._plot


class BoundHeadTargetPlots:
    """Plotting facade for model-bound head targets (``model.targets.heads.plot``).

    Every picture is drawn from the bound targets' ``compare()`` on the bound
    model; a ``baseline`` model adds the same comparison on a second run.
    """

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
        """Plot target locations on a map axis (:meth:`BoundHeadTargets.plot_locations`).

        Parameters
        ----------
        ax : matplotlib.axes.Axes or None, optional
            Axes to draw on. Default ``None`` makes a new figure and axes.
        **kwargs
            Forwarded unchanged to ``geopandas.GeoDataFrame.plot`` (for example
            ``color=``, ``markersize=``); one it does not take raises there.

        Returns
        -------
        matplotlib.axes.Axes
            The axes, titled "Head Target Locations".

        Raises
        ------
        ValueError
            When the locations carry no geometry -- targets given as cell ids.
        """

        return self.targets.plot_locations(ax=ax, **kwargs)

    def calibration(self, *, type: str = "calibration"):
        """Return a calibration/heads plot for the bound targets.

        Parameters
        ----------
        type : str, optional
            Which plot, as :meth:`BoundHeadTargets.calibration_plot` reads it:
            ``"calibration"`` (default, the cross plot with fit statistics),
            ``"obs_vs_sim"`` (without them), ``"heads"`` (values through time,
            one pair of lines per location and layer it is observed in) or
            ``"residuals_by_period"``; any other string raises ``ValueError``.

        Returns
        -------
        myflopy.modflow.calcs.calibration.CalibrationPlot
            A Plotly figure (a ``viz.Fig``). There is no ``backend``: use
            :meth:`obs_vs_sim` or :meth:`residuals_by_period` for a static one.
        """

        return self.bound_targets.calibration_plot(type=type)

    def obs_vs_sim(self, *, baseline=None, backend: str = "plotly"):
        """Return a target-vs-simulated cross plot for one or two models.

        One point per comparison row holding both a target and a simulated
        value, in the field's units (model length units for a head), with a 1:1
        line spanning every value drawn; titled "Observed vs simulated heads",
        or "... concentrations" on a bound ``ConcTargets`` (it read "heads" there
        until ledger 477). A simulated value that is MODFLOW's
        no-value marker (``|value| >= 1e29``) is not drawn (ledger 355). When no
        row of the bound model's comparison holds both values, its series is
        empty, the figure carries a note naming the missing side, and a
        ``UserWarning`` says so.

        Parameters
        ----------
        baseline : SimulationBase or None, optional
            A second run model: the SAME targets compared against it are drawn
            as a "Baseline" series beside the bound model's "Model" series.
            Default ``None`` draws the bound model alone. Its comparison raises
            as the bound model's does.
        backend : str, optional
            ``"plotly"`` (default; also ``"interactive"``) returns a
            :class:`~myflopy.modflow.calcs.calibration.CalibrationPlot`;
            ``"matplotlib"`` (also ``"mpl"``, ``"static"``) a static
            ``matplotlib.figure.Figure``. Case and surrounding space are
            ignored; anything else raises ``ValueError``, after the comparisons
            have run.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            The cross plot, per ``backend``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(
            current,
            baseline_compare=baseline_frame,
            title=f"Observed vs simulated {self.targets._FIELD_LABEL}s",
            backend=backend,
        )

    def timeseries(self, name: str | None = None, *, baseline=None, backend: str = "plotly"):
        """Return a time-series plot for one target location.

        The target's values ("Target") and the bound model's ("Model") against
        the target's ``time`` -- the zero-based stress period -- in the field's
        units, titled with the target's name; the y axis is titled "Head", or
        "Concentration" on a bound ``ConcTargets`` (it read "Head" there until
        ledger 477). A simulated no-value marker
        (``|value| >= 1e29``) is a gap in the line (ledger 355).

        Parameters
        ----------
        name : str or None, optional
            The target to draw, matched against the comparison's ``name``
            column ignoring case (and the argument's surrounding space). Default
            ``None`` draws the set's only name, and raises ``ValueError`` when
            the set holds more than one; a name matching no row raises
            ``ValueError`` too. A well given in two layers (no screen) is drawn as
            one line per layer (one line zig-zagged between them until ledger
            477); a screened well is one observation, one line.
        baseline : SimulationBase or None, optional
            A second run model whose simulated values for the same target are
            drawn as a "Baseline" line. Default ``None`` draws the bound model
            alone.
        backend : str, optional
            ``"plotly"`` (default; also ``"interactive"``) or ``"matplotlib"`` (also
            ``"mpl"``, ``"static"``), as in :meth:`obs_vs_sim`;
            anything else raises ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            The time series, per ``backend``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
            yaxis_title=self.targets._FIELD_LABEL.capitalize(),
            backend=backend,
        )

    def residuals_by_period(self, *, baseline=None, backend: str = "plotly"):
        """Return a by-period residual summary plot.

        The mean absolute residual of the measured rows in each ``time`` (the
        zero-based stress period), "Model MAE", titled "Residual MAE by period";
        a period with no measured row is a gap, and a residual against MODFLOW's
        no-value marker is left out of its mean (ledger 355).

        Parameters
        ----------
        baseline : SimulationBase or None, optional
            A second run model whose per-period MAE for the same targets is drawn
            as "Baseline MAE". Default ``None`` draws the bound model alone.
        backend : str, optional
            ``"plotly"`` (default; also ``"interactive"``) or ``"matplotlib"`` (also
            ``"mpl"``, ``"static"``), as in :meth:`obs_vs_sim`;
            anything else raises ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            The period summary, per ``backend``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.targets.compare(self.model)
        baseline_frame = None if baseline is None else self.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(
            current,
            baseline_compare=baseline_frame,
            backend=backend,
        )
