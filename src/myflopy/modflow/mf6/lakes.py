"""Build MODFLOW 6 LAK package specifications from lake geometries."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Integral, Real
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from myflopy.modflow.mf6.mvr import MoverConnection

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import CRS
from rasterio.features import geometry_mask
from shapely.geometry import box
from shapely.geometry.base import BaseGeometry

from myflopy._length_units import mf6_length_conversion
from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package
from myflopy.advanced import lak_spec
from myflopy.areas import Area
from myflopy.grid_identity import (
    LAKE_CONNECTIONS,
    LAKE_RULES,
    _disv_arrays,
    _listing,
    bed_tolerance,
    lak_geometry_complaints,
    pin_metadata,
)
from myflopy.modflow.mf6.grid.connectivity import DisvFaces
from myflopy.modflow.mf6.grid.geometry import (
    domain_active,
    footprint_cells,
    read_footprints,
    refuse_footprint_misses,
)
from myflopy.modflow.mf6.lak_implicit import ModflowGwflakImplicit
from myflopy.modflow.utils.datatypes.readers import _vector_layer_name
from myflopy.specs import (
    ModelContext,
    PackageSpec,
    _layer_elevation_source,
    mf6_time_conversion,
)
from myflopy.surfaces import Surface, sample_bed

logger = get_logger(__name__)

#: A HORIZONTAL face thinner than this (model length units) is not written. MF6 reads
#: a face with ``belev == telev`` as the WHOLE cell (``gwf-lak.f90``), and FloPy writes
#: connection elevations to 8 fixed decimals, so a face of 1e-9 written as ``99.90000000
#: 99.90000000`` became a whole-cell face (measured). 1e-6 is 100x that precision; the
#: thinnest face on a real pond's legacy connections is 1.3e-4, so the tolerance drops
#: none of them (the excavation design's F9, ledger 325).
FACE_TOL = 1.0e-6

#: What ``unexcavated=`` accepts: the three words of the raster sampler's ``gaps=``
#: that make sense here -- there is no "nearest" way to excavate a lake.
_UNEXCAVATED_WORDS = ("warn", "raise", "ignore")

#: MF6's lake STATUS words (``gwf-lak.f90``). MF6 reads them in any case, and refuses
#: anything else by name, one run late (ledger 318, measured with ``BOGUS``).
_LAKE_STATUSES = ("ACTIVE", "INACTIVE", "CONSTANT")

#: The period settings MF6 refuses below zero, each by name (ledger 318, measured all
#: five). ``STAGE`` and an outlet's ``RATE`` take either sign.
_NONNEGATIVE_SETTINGS = ("RAINFALL", "EVAPORATION", "RUNOFF", "INFLOW", "WITHDRAWAL")

#: The ``mf.lak`` parameter each period keyword comes from, for the refusals.
_SETTING_NAMES = {
    "STAGE": "stage",
    "RAINFALL": "rainfall",
    "EVAPORATION": "evaporation",
    "RUNOFF": "runoff",
    "WITHDRAWAL": "withdrawals",
    "INFLOW": "inflow",
    "STATUS": "status",
}


def _plain_scalar(value: Any) -> Any:
    """A numpy scalar as its python value (ledger 265); anything else unchanged."""

    return value.item() if isinstance(value, np.generic) else value


def _per_period_form(value: Any, *, name: str) -> Any:
    """``value`` with every numpy/pandas container replaced by the python form it MEANS.

    A per-period value (a lake forcing, a ``LakeOutlet.rate``) is a scalar, a list by
    period, or a mapping. A ``numpy`` array or ``pandas`` Series is neither a python
    ``Sequence`` nor a ``Mapping``, so it was taken for a scalar and written WHOLE into
    every period -- ``1  RATE  [-1. -2.]`` -- and MF6 died reading it (a segmentation
    fault, ledger 318). A one-dimensional array is the list it holds; a Series is the
    mapping of its index, so an integer index is periods, exactly as a list's positions
    are, and a lake-id index is lakes; a 0-d array or numpy scalar is its python value.
    Applied through lists, tuples and mapping values, so ``{"pond": np.array(...)}``
    is read the same way. An array of two or more dimensions is refused: which axis is
    the period is a guess.
    """

    if isinstance(value, pd.Series):
        return {_plain_scalar(key): _per_period_form(item, name=name) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return _plain_scalar(value[()])
        if value.ndim > 1:
            raise ValueError(
                f"{name} is a {value.ndim}-D array of shape {value.shape}; a per-period "
                "value is ONE number per period -- a 1-D array or list (position 0 is "
                "the first stress period) or a mapping. Which axis is the period would "
                "be a guess, so it is refused."
            )
        return [_per_period_form(item, name=name) for item in value.tolist()]
    if isinstance(value, Mapping):
        return {_plain_scalar(key): _per_period_form(item, name=name) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        items = [_per_period_form(item, name=name) for item in value]
        return tuple(items) if isinstance(value, tuple) else items
    return _plain_scalar(value)


def _is_period_key(key: Any) -> bool:
    """Whether a mapping key reads as a stress period: any integer, numpy's included, never a bool."""

    return isinstance(key, Integral) and not isinstance(key, (bool, np.bool_))


def _lake_keyed(mapping: Mapping, lake_id: str, *, name: str) -> tuple[bool, Any]:
    """``(found, value)`` for ``lake_id`` in a ``{lake_id: ...}`` mapping whose keys may be ints.

    Lake ids are STRINGS -- the ``lake_id_field`` value through ``str()`` -- so an
    integer key ``101`` names lake ``"101"``, the way :meth:`LAKBuilder.lake_number`
    reads one. It used to be looked up as the int and missed: ``only_vertical={101:
    True}`` was ignored silently, and ``tables={101: ...}`` raised "unknown lake
    '101'", printing the int exactly as the id it missed (ledger 318). A tuple key
    (``bed_leakance``'s ``(lake_id, "vertical")``) is not a lake id and is skipped
    here. Two keys naming one lake (``101`` and ``"101"``) raise ``ValueError``.
    """

    hits = [
        key
        for key in mapping
        if not isinstance(key, tuple) and str(_plain_scalar(key)) == lake_id
    ]
    if len(hits) > 1:
        raise ValueError(
            f"{name} names lake {lake_id!r} twice ({', '.join(repr(key) for key in hits)}); "
            "lake ids are strings, so an integer key names the lake of that string."
        )
    if not hits:
        return False, None
    return True, mapping[hits[0]]


def _series_lake_labels(series: pd.Series, lake_ids: Sequence[str], *, name: str) -> dict[str, Any] | None:
    """``{lake_id: label}`` when ``series``' index names every lake, else ``None``.

    A label names a lake as a mapping key does (:func:`_lake_keyed`): through
    ``str()``, so ``101`` names lake ``"101"``. The Series branch of
    :meth:`LAKBuilder._lake_value` used to test ``lake in series.index`` with the
    STRING id, so a Series keyed by the integer ids failed it and fell through to the
    positional reading: ``starting_stage=pd.Series({102: 12.0, 101: 11.5})`` gave lake
    ``"101"`` 12.0 and ``"102"`` 11.5, where the same dict gives each its own (ledger
    318, review). Two labels naming one lake (``101`` and ``"101"``, or a repeated
    label) raise ``ValueError``.
    """

    wanted = set(lake_ids)
    found: dict[str, list[Any]] = {}
    for label in series.index:
        lake = str(_plain_scalar(label))
        if lake in wanted:
            found.setdefault(lake, []).append(label)
    for lake, labels in found.items():
        if len(labels) > 1:
            raise ValueError(
                f"{name} names lake {lake!r} twice ({', '.join(repr(label) for label in labels)}); "
                "lake ids are strings, so an integer label names the lake of that string."
            )
    if len(found) < len(wanted):
        return None
    return {lake: labels[0] for lake, labels in found.items()}


@dataclass(frozen=True, slots=True)
class LakeConnection:
    """One lake-aquifer connection you write yourself: one MF6 ``CONNECTIONDATA`` record.

    A LAK lake exchanges water with the aquifer only through its connections, one
    record per connected cell. ``mf.lak`` normally GENERATES them from the lake
    polygon (``connection_modes="bathy"`` or ``"rectangular"``). Build these by hand
    only to give one lake its connections explicitly:
    ``mf.lak(..., connection_modes={lake_id: [LakeConnection(...), ...]})``. That
    lake's connections are then written exactly as given, in the order given (MF6's
    connection number, from zero), and ``lake_bottom``, ``lake_top``,
    ``only_vertical``, ``only_layer`` and ``rectangular_interior`` do not touch them.

    Nothing here is checked against the grid: not the cell id, not the elevations,
    not whether the cell is active. FloPy checks the cell id against the idomain when
    the package is constructed. MF6 checks the rest when it reads the file, as each
    entry below says -- except the two a HORIZONTAL connection cannot do without,
    its ``connection_length`` and ``connection_width``, which are refused here when
    they are not greater than zero (ledger 318). All lengths are in the model's
    length unit.

    Parameters
    ----------
    cellid : tuple of (int, int)
        Required. The aquifer cell as ``(layer, cell)``, both zero-based, on the
        model's DISV grid (``mf.lak`` builds on DISV only). Written as given.
    connection_type : str, default "VERTICAL"
        ``"VERTICAL"`` or ``"HORIZONTAL"``, in any case (it is stored upper-case).
        Anything else raises ``ValueError``, including MF6's ``EMBEDDEDH`` and
        ``EMBEDDEDV``. What MF6 does with each (``gwf-lak.f90``):

        * ``"VERTICAL"`` -- exchange through the lakebed over the cell, using the
          cell's K33 (K11 when NPF gives no K33). Without a lake table, the cell's whole plan area is also lake
          surface area once the stage rises above the cell's top, and that is where
          the lake's volume over this cell starts. **MF6 IGNORES
          ``bottom_elevation`` and ``top_elevation`` for this type**: it replaces
          them with the cell's TOP and top + ``surfdep`` (the MF6 input guide:
          "Any value can be specified if CLAKTYPE is VERTICAL"). Measured: one
          VERTICAL connection under a cell whose top is 8.0, lake held at 12.0,
          elevations written as 3.0, 6.0 and 0.0; MF6 reported ``VOLUME`` 400.0
          (100 ft2 x (12 - 8)) and ``SURFACE-AREA`` 100.0 all three times.
        * ``"HORIZONTAL"`` -- exchange through a vertical face of the cell (a lake
          side wall), using the cell's K11, ``connection_length`` and
          ``connection_width`` and the two elevations. It adds no surface area and
          no volume.
    bed_leakance : float, optional
        Default ``None``. Lakebed leakance for this connection (1/T: the bed's
        hydraulic conductivity divided by its thickness). ``None`` takes the lake's
        value for this connection type from ``mf.lak(bed_leakance=...)``. A number
        here wins over that and is written as given, unchecked: ``0.0`` passes,
        and MF6 then moves no water through the connection (measured: flow 0.0),
        where ``mf.lak`` refuses a resolved ``bed_leakance`` of zero or less.
    bottom_elevation : float, default 0.0
        MF6's ``BELEV``, the bottom of the wetted face (length). Read for a
        HORIZONTAL connection only, where MF6 requires it to be at or above the
        cell's bottom. Ignored for VERTICAL (see ``connection_type``). When it
        EQUALS ``top_elevation``, as the two ``0.0`` defaults do, MF6 resets the
        face to the whole cell, bottom to top (measured: the same flow as a face
        written from the cell bottom to its top).
    top_elevation : float, default 0.0
        MF6's ``TELEV``, the top of the wetted face (length). HORIZONTAL only; MF6
        refuses one above the cell's top (measured: ``telev ... MUST BE <= cell
        top``). Ignored for VERTICAL.
    connection_length : float, default 0.0
        MF6's ``CONNLEN`` (length): the distance from the cell's node to the lake
        face, used in a HORIZONTAL connection's conductance. For HORIZONTAL it
        must be a finite number greater than zero, so the ``0.0`` default is
        refused here (``ValueError``) -- MF6 refuses it by name when it reads the
        file (measured: ``connection length (connlen) ... MUST BE > 0``), one run
        later. Ignored, and not checked, for VERTICAL.
    connection_width : float, default 0.0
        MF6's ``CONNWIDTH`` (length): the width of a HORIZONTAL face. For
        HORIZONTAL it must be a finite number greater than zero, so the ``0.0``
        default is refused here (``ValueError``): MF6 accepts a zero width, and
        the connection then carries no water at all (measured: flow 0.0), which
        is never what a wall written by hand is for. Ignored, and not checked,
        for VERTICAL.

    Raises
    ------
    ValueError
        If ``connection_type`` is not ``"VERTICAL"`` or ``"HORIZONTAL"``, in any
        letter case; or, for a HORIZONTAL connection, if ``connection_length``
        or ``connection_width`` is not a finite number greater than zero.

    See Also
    --------
    mf.lak : Takes these through ``connection_modes={lake_id: [...]}``, and
        generates them for the other lakes.

    Examples
    --------
    >>> floor = mf.LakeConnection((0, 42))
    >>> wall = mf.LakeConnection((1, 43), "horizontal", bottom_elevation=90.0,
    ...                          top_elevation=95.0, connection_length=12.5,
    ...                          connection_width=25.0)
    >>> wall.connection_type
    'HORIZONTAL'

    Then ``mf.lak(..., connection_modes={"vault": [floor, wall]})``. A wall
    missing its length is refused here rather than by MF6:

    >>> mf.LakeConnection((1, 43), "horizontal", connection_width=25.0)
    Traceback (most recent call last):
    ...
    ValueError: A HORIZONTAL LakeConnection at (1, 43) needs connection_length > 0 ...
    """

    cellid: tuple[int, int]
    connection_type: str = "VERTICAL"
    bed_leakance: float | None = None
    bottom_elevation: float = 0.0
    top_elevation: float = 0.0
    connection_length: float = 0.0
    connection_width: float = 0.0

    def __post_init__(self) -> None:
        """Uppercase and validate ``connection_type``; a HORIZONTAL one needs a length and a width."""

        connection_type = self.connection_type.upper()
        if connection_type not in {"VERTICAL", "HORIZONTAL"}:
            raise ValueError("connection_type must be 'VERTICAL' or 'HORIZONTAL'.")
        object.__setattr__(self, "connection_type", connection_type)
        if connection_type == "HORIZONTAL":
            # MF6 refuses a zero CONNLEN by name, one run late, and runs a zero
            # CONNWIDTH carrying nothing (ledger 318, measured). Both defaults are 0.0,
            # so a wall written without them is refused here.
            for field_name, consequence in (
                ("connection_length", "MF6 refuses a HORIZONTAL connection without one"),
                ("connection_width", "MF6 runs one of zero width carrying no water at all"),
            ):
                value = getattr(self, field_name)
                if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not (
                    np.isfinite(float(value)) and float(value) > 0.0
                ):
                    raise ValueError(
                        f"A HORIZONTAL LakeConnection at {self.cellid} needs {field_name} > 0 "
                        f"(got {value!r}): {consequence}. The default 0.0 is for a VERTICAL "
                        "connection, which reads neither."
                    )


@dataclass(frozen=True, slots=True)
class LakeOutlet:
    """One lake outlet: an MF6 ``OUTLETS`` record, naming its lakes by lake id.

    An outlet moves water out of a lake by a discharge law, to another lake or out
    of the model; a ``SPECIFIED`` outlet can also add water. ``source`` and
    ``receiver`` are LAKE IDS -- the strings ``mf.lak`` names lakes by (the
    ``lake_id_field`` values) -- never MF6's lake numbers. ``mf.lak`` resolves them
    when it builds and numbers the outlets by their position in ``outlets=``, from
    zero. Nothing here is checked until then: the ids, and every number, are stored
    as given, except that a ``rate`` held in a ``numpy`` array or ``pandas`` Series
    is stored as the list or mapping it means (see ``rate``).

    Parameters
    ----------
    source : str or int
        Required. The lake the outlet draws from (MF6 ``LAKEIN``), by lake id.
        ``mf.lak`` converts it with ``str()``, so an integer ``101`` finds lake
        ``"101"``. An id naming no lake raises ``ValueError: Unknown lake ID``
        when ``mf.lak`` builds.
    receiver : str, int or None, default None
        The lake the outflow goes to (MF6 ``LAKEOUT``), by lake id, resolved like
        ``source``. ``None`` sends it out of the model (written as ``LAKEOUT 0``,
        MF6's external boundary). A receiver equal to ``source`` is accepted.
    outlet_type : str, default "SPECIFIED"
        The discharge law (MF6 ``COUTTYPE``), in any case (stored upper-case). With
        ``d = stage - invert`` (``gwf-lak.f90``):

        * ``"SPECIFIED"`` -- the flow is ``rate``; nothing else is read.
        * ``"WEIR"`` -- a sharp-crested weir, ``2/3 x 0.61 x width x d x sqrt(2 g
          d)`` while ``d > 0``; reads ``invert`` and ``width``.
        * ``"MANNING"`` -- ``(1 / roughness) x width x d**(5/3) x sqrt(slope)``
          while ``d > 0``; reads ``invert``, ``width``, ``roughness`` and
          ``slope``.

        ``g`` and Manning's constant are SI values that MF6 converts with
        ``mf.lak``'s ``length_conversion`` and ``time_conversion``. Anything else
        raises ``ValueError``.
    invert : float, default 0.0
        WEIR and MANNING: the outlet's invert elevation (length). No water leaves
        while the stage is at or below it. Ignored for SPECIFIED.
    width : float, default 0.0
        WEIR and MANNING: the outlet width (length). With the ``0.0`` default the
        outlet never flows. Ignored for SPECIFIED.
    roughness : float, default 0.0
        MANNING only: Manning's n (MF6 ``ROUGH``). MF6 divides by it only when it
        is above zero, so the ``0.0`` default gives no flow.
    slope : float, default 0.0
        MANNING only: the outlet channel's bed slope (dimensionless). The ``0.0``
        default gives no flow.
    rate : float, str, list, tuple, dict of int to float, numpy.ndarray or pandas.Series, optional
        Default ``None``. The SPECIFIED flow rate (length cubed per time), written
        as a ``RATE`` period setting. **A NEGATIVE rate removes water from the
        lake and a POSITIVE rate ADDS water** (MF6's input guide: "A positive value
        indicates inflow and a negative value indicates outflow"). Measured on a
        100 ft2 lake: ``+5`` raised the stage 0.05 ft per day and ``-5`` lowered
        it 0.05. MF6 caps a negative rate at the water available to the lake in
        the time step. The forms:

        * a number: every stress period (a ``str`` is written as given, which is
          how MF6 names a time series);
        * a list or tuple, or a one-dimensional ``numpy`` array: position 0 is the
          first stress period, and ``None`` in a list skips that period;
        * ``{period: value}``: keyed by ZERO-based stress period, so
          ``{1: -2.0}`` is MF6's ``PERIOD 2``;
        * a ``pandas.Series``: its index is the period, as a dict's keys are.

        An array is stored here as the list it holds and a Series as the dict of
        its index; each used to be written WHOLE into every period, and MF6 died
        reading that (measured, a segmentation fault; ledger 318). An array of
        two or more dimensions raises ``ValueError`` here. A period the value
        does not reach -- past the end of the list, or with no key -- is written
        with no RATE, and MF6 keeps the previous period's rate (measured); before
        the first one, the rate is MF6's default, 0. ``None`` writes no RATE at
        all. ``mf.lak`` refuses (``ValueError``) a list with more values than
        stress periods, a period key outside ``range(nper)`` -- neither would
        ever be written -- and a number that is not finite. It writes the RATE
        for a WEIR or MANNING outlet too, and MF6 ignores it.

    Raises
    ------
    ValueError
        If ``outlet_type`` is not ``"SPECIFIED"``, ``"WEIR"`` or ``"MANNING"``, in
        any letter case, or ``rate`` is an array of two or more dimensions.

    See Also
    --------
    mf.lak : Takes these through ``outlets=[...]`` and resolves their lake ids.

    Examples
    --------
    A pumped withdrawal of 250 ft3/d out of the model, from the second stress
    period on, and an overflow weir from one lake into another:

    >>> pump = mf.LakeOutlet("north_pond", rate={1: -250.0})
    >>> spill = mf.LakeOutlet("north_pond", "south_pond", outlet_type="weir",
    ...                       invert=531.5, width=4.0)
    >>> spill.outlet_type
    'WEIR'

    Then ``mf.lak(..., outlets=[pump, spill])``.
    """

    source: str
    receiver: str | None = None
    outlet_type: str = "SPECIFIED"
    invert: float = 0.0
    width: float = 0.0
    roughness: float = 0.0
    slope: float = 0.0
    rate: Any = None

    def __post_init__(self) -> None:
        """Uppercase and validate ``outlet_type``; store an array or Series ``rate`` as what it means."""

        outlet_type = self.outlet_type.upper()
        if outlet_type not in {"SPECIFIED", "MANNING", "WEIR"}:
            raise ValueError("outlet_type must be 'SPECIFIED', 'MANNING', or 'WEIR'.")
        object.__setattr__(self, "outlet_type", outlet_type)
        object.__setattr__(self, "rate", _per_period_form(self.rate, name="LakeOutlet.rate"))


#: The column names a ``DataFrame`` of lake-table rows is read by, in the file's order.
_TABLE_COLUMNS = ("stage", "volume", "area")


def _table_rows(rows: Any) -> tuple[tuple[float, ...], ...]:
    """A lake table's rows as tuples of python ``float``; a DataFrame read by name or position.

    Iterating a bare ``pandas.DataFrame`` yields its column LABELS, so a frame used to
    fail ``float()`` on ``'s'`` (string labels) or on ``0`` (``'int' object is not
    iterable``) -- a message about neither the frame nor the table (ledger 318). It is
    read by its ``stage``/``volume``/``area`` columns when it names all three, by
    position when it has exactly three, and refused otherwise: which column is which
    would be a guess.
    """

    if isinstance(rows, (str, bytes, Mapping)):
        raise TypeError(
            "A lake table's rows are a sequence of (stage, volume, area) rows, a 2-D array "
            f"or a DataFrame, not a {type(rows).__name__}."
        )
    if isinstance(rows, pd.DataFrame):
        by_name = {str(column).strip().lower(): column for column in rows.columns}
        named = [word for word in _TABLE_COLUMNS if word in by_name]
        if len(named) == len(_TABLE_COLUMNS):
            rows = rows[[by_name[word] for word in _TABLE_COLUMNS]]
        elif named:
            missing = [word for word in _TABLE_COLUMNS if word not in by_name]
            raise ValueError(
                f"The lake table's DataFrame has a {' and a '.join(named)} column but no "
                f"{' or '.join(missing)}: name all three (stage, volume, area), or give "
                "exactly three columns, read in that order. Which column is which would be "
                "a guess."
            )
        elif rows.shape[1] != len(_TABLE_COLUMNS):
            raise ValueError(
                f"The lake table's DataFrame has {rows.shape[1]} columns "
                f"({', '.join(repr(column) for column in rows.columns)}), none named stage, "
                "volume or area: name those three, or give exactly three columns, read in "
                "that order."
            )
        rows = rows.to_numpy()
    return tuple(tuple(float(value) for value in row) for row in rows)


@dataclass(frozen=True, slots=True)
class LakeTable:
    """One lake's validated stage / volume / surface-area table, ready for MF6.

    What ``mf.lak`` writes as a LAK lake-table (``laktab``) child file: for a
    series of stages, the lake's stored volume and its surface area. It is
    immutable, and it is what :meth:`LakeTableBuilder.build` returns. Build one
    directly only when you already hold the final rows -- ``LakeTableBuilder(rows=...)``
    is the same thing, plus the DEM and prism modes. ``mf.lak(tables=...)``
    accepts either, keyed by lake id.

    Parameters
    ----------
    rows : sequence of (stage, volume, area) sequences, numpy.ndarray or pandas.DataFrame
        Required. One row per stage, each exactly three numbers: the stage (model
        length units), the lake's stored volume at that stage (length cubed) and
        its surface area (length squared). Every value goes through ``float()``,
        so a list of tuples, a 2-D ``numpy`` array and numeric strings all work,
        and the rows are stored as a tuple of 3-tuples of ``float``. A
        ``pandas.DataFrame`` is read by its column NAMES when it has columns
        named ``stage``, ``volume`` and ``area`` (in any letter case; other
        columns are ignored), and otherwise by POSITION when it has exactly three
        columns; its index is ignored. A frame naming only some of the three, or
        naming none with other than three columns, raises ``ValueError`` -- which
        column is which would be a guess. No unit conversion is done. Checked
        here, each failure a ``ValueError``; the first three are MF6's own table
        checks (``gwf-lak.f90``, ``lak_read_table``), so a table that passes here
        is one MF6 reads:

        * stages STRICTLY increasing -- MF6: ``TABLE STAGE ENTRY ... MUST BE
          GREATER THAN THE PREVIOUS STAGE ENTRY``;
        * volumes STRICTLY increasing, so two rows of equal volume (two
          zero-volume rows below a lake's bed, say) are refused -- MF6: ``TABLE
          VOLUME ENTRY ... MUST BE GREATER THAN THE PREVIOUS VOLUME ENTRY``;
        * surface areas NONDECREASING (an equal area passes) -- MF6: ``TABLE
          SURFACE AREA ENTRY ... MUST BE GREATER THAN OR EQUAL TO THE PREVIOUS
          SURFACE AREA ENTRY``;
        * at least one row (MF6: ``NROW MUST BE > 0``) and exactly three values in
          every row;
        * every value finite, and no NEGATIVE volume or area. MF6's own checks
          compare neighbouring rows only (measured, ledger 318): it runs a
          negative volume or area as given, which no lake holds (a first volume
          of -5.0 read ``VOLUME`` -2.0 at a stage inside the table, and a first
          area of -4.0 ``SURFACE-AREA`` -3.0), and a ``NaN``
          reaches it only as a row FloPy wrote cut short at the ``NaN``, which it
          stops on without naming the value, the column or the lake (``Could not
          read variable of type DOUBLE PRECISION from the following line:
          '9.00000000'``).

        These are the rules :class:`LakeTableBuilder` and ``mf.lak`` rely on; it
        used to accept an equal volume, a falling area, a ``NaN`` and a negative
        volume, and MF6 refused the first three one run later and ran the last
        (ledger 318).
    filename : str, optional
        Default ``None``. Name of the ``laktab`` file ``mf.lak`` writes for this
        table, beside the model's other input files. ``None`` or ``""`` lets
        ``mf.lak`` name it ``"<package name>_<lake_id>.lak.tab"``, e.g.
        ``lak_pond.lak.tab``.

    Attributes
    ----------
    minimum_stage : float
        The first row's stage, the lowest in the table.
    maximum_stage : float
        The last row's stage, the highest in the table.

    Raises
    ------
    ValueError
        On construction, if ``rows`` breaks any rule listed under ``rows``, a
        value is a string ``float()`` cannot convert, or a ``DataFrame``'s
        columns cannot be told apart.
    TypeError
        On construction, if a row or a value is not something ``float()`` or
        iteration accepts, or ``rows`` is a mapping or a string.

    See Also
    --------
    LakeTableBuilder : Builds one from explicit rows, a bathymetry raster
        (including ``mf.Contours``) or a constant-area prism.
    mf.lak : Writes it, via ``tables={lake_id: table}``.

    Examples
    --------
    >>> table = mf.LakeTable([(96.0, 0.0, 0.0), (98.0, 1.5e4, 1.2e4)], filename="pond.tab")
    >>> table.minimum_stage, table.maximum_stage
    (96.0, 98.0)

    A survey held in a DataFrame, read by its column names:

    >>> import pandas as pd
    >>> survey = pd.DataFrame({"Area": [0.0, 1.2e4], "Stage": [96.0, 98.0],
    ...                        "Volume": [0.0, 1.5e4]})
    >>> mf.LakeTable(survey).rows
    ((96.0, 0.0, 0.0), (98.0, 15000.0, 12000.0))

    Two rows of equal volume are refused here, as MF6 refuses them:

    >>> mf.LakeTable([(95.0, 0.0, 0.0), (96.0, 0.0, 0.0), (98.0, 1.5e4, 1.2e4)])
    Traceback (most recent call last):
    ...
    ValueError: Lake table volumes must STRICTLY increase, as MF6 requires: ...
    """

    rows: tuple[tuple[float, float, float], ...]
    filename: str | None = None

    def __post_init__(self) -> None:
        """Coerce, validate, and freeze the ``(stage, volume, area)`` rows by MF6's table rules."""

        rows = _table_rows(self.rows)
        if not rows:
            raise ValueError("A lake table must contain at least one row.")
        if any(len(row) != 3 for row in rows):
            raise ValueError("Lake table rows must contain stage, volume, and surface area.")
        for row in rows:
            for label, value in zip(("stage", "volume", "surface area"), row, strict=True):
                if not np.isfinite(value):
                    raise ValueError(
                        f"Lake table row {row} has a {label} of {value!r}: every value must be "
                        "finite. FloPy writes the row cut short at a NaN, and MF6 stops unable "
                        "to read it, naming neither the value nor the lake."
                    )
            if row[1] < 0.0 or row[2] < 0.0:
                raise ValueError(
                    f"Lake table row {row} has a negative "
                    f"{'volume' if row[1] < 0.0 else 'surface area'}: no lake holds one. MF6 "
                    "does not check the sign, and would run it as given."
                )
        for previous, current in zip(rows, rows[1:], strict=False):
            if not current[0] > previous[0]:
                raise ValueError(
                    "Lake table stages must be strictly increasing, as MF6 requires: stage "
                    f"{current[0]!r} follows {previous[0]!r}."
                )
            if not current[1] > previous[1]:
                hint = (
                    " Two zero-volume rows are stages at or below the lake's lowest bed: "
                    "keep only the last of them."
                    if current[1] == 0.0
                    else ""
                )
                raise ValueError(
                    "Lake table volumes must STRICTLY increase, as MF6 requires: the row at "
                    f"stage {current[0]!r} holds {current[1]!r}, not more than the "
                    f"{previous[1]!r} at stage {previous[0]!r}.{hint}"
                )
            if current[2] < previous[2]:
                raise ValueError(
                    "Lake table surface areas must not fall as the stage rises, as MF6 "
                    f"requires: {current[2]!r} at stage {current[0]!r} is below the "
                    f"{previous[2]!r} at stage {previous[0]!r}."
                )
        object.__setattr__(self, "rows", rows)

    @property
    def minimum_stage(self) -> float:
        """The first row's stage, the lowest in the table (stages strictly increase).

        ``mf.lak`` refuses a ``starting_stage`` below it. While the model runs,
        MF6 treats any stage below it as this row: the lake keeps the first row's
        volume and surface area.

        See Also
        --------
        maximum_stage : The other end of the table.

        Examples
        --------
        >>> mf.LakeTable([(96.0, 0.0, 0.0), (98.0, 1.5e4, 1.2e4)]).minimum_stage
        96.0
        """

        return self.rows[0][0]

    @property
    def maximum_stage(self) -> float:
        """The last row's stage, the highest in the table.

        ``mf.lak`` refuses a ``starting_stage`` above it. While the model runs,
        MF6 extends the table above it: the surface area stays at the last row's
        value, and the volume grows by that FULL area per unit of stage rise
        (``gwf-lak.f90``, ``lak_calculate_vol``). A
        :class:`LakeTableBuilder` ``storage_coefficient`` is not applied there,
        because MF6 never sees it.

        See Also
        --------
        minimum_stage : The other end of the table.

        Examples
        --------
        >>> mf.LakeTable([(96.0, 0.0, 0.0), (98.0, 1.5e4, 1.2e4)]).maximum_stage
        98.0
        """

        return self.rows[-1][0]


@dataclass(frozen=True, slots=True)
class LakeTableBuilder:
    """Recipe for one lake's stage / volume / surface-area table (a :class:`LakeTable`).

    MF6's LAK package can read one table per lake, a ``laktab`` file. It lists,
    for a series of stages, the lake's stored volume and its surface area.
    Without a table, MF6 computes both from the lake's VERTICAL connections alone:
    each adds its connected cell's full plan area once the stage rises above that
    CELL'S TOP. MF6 replaces the bottom and top elevations written for a VERTICAL
    connection with the cell top and top + ``surfdep`` (``gwf-lak.f90``), so the
    bed elevations ``mf.lak`` writes do not enter the volume -- only the layer
    each connection lands in does. HORIZONTAL connections add no area or volume.
    A table replaces only that volume/area calculation; the lake-aquifer exchange
    through VERTICAL and HORIZONTAL connections comes from the connections either
    way.

    This object only holds the inputs. :meth:`build` turns them into a validated
    :class:`LakeTable` -- validated by MF6's own table checks, so a table it
    returns is one MF6 reads. Hand either one to ``mf.lak(..., tables={lake_id:
    ...})``, which calls :meth:`build` for you (Notes, *Handing it to mf.lak*).

    **Three modes. The fields you set choose the mode.** :meth:`build` checks
    them in this order and uses the FIRST that applies. It does not read or
    validate the fields of the other two modes, so a leftover field is ignored
    silently:

    1. **Explicit rows** -- ``rows`` is set. Taken as given.
    2. **DEM** -- ``rows`` is ``None`` and ``dem`` is set. Needs ``footprint``
       and ``stages``. Optional: ``footprint_layer``, ``footprint_buffer``,
       ``storage_coefficient``.
    3. **Prism** -- neither of the above. Needs ``area``, ``top`` and
       ``bottom``. Optional: ``stage_step``, ``storage_coefficient``. The surface
       area is constant: a trench, vault or lined basin.

    ``filename`` applies to all three.

    **No unit conversion happens anywhere in this class.** Every length (stages,
    ``top``/``bottom``, the DEM's elevations AND its pixel size,
    ``footprint_buffer``) must already be in the MODEL's length unit. Areas come
    out in that unit squared and volumes in that unit cubed. The footprint IS
    reprojected to the raster's CRS, but that moves it, it does not convert
    units: a DEM in a metre CRS under a model in feet gives square metres and
    cubic metres, and nothing notices. A DEM in a GEOGRAPHIC CRS, whose pixels
    measure degrees, is refused.

    Parameters
    ----------
    rows : sequence of (stage, volume, area) sequences, numpy.ndarray or pandas.DataFrame, optional
        Default ``None``. **Rows mode.** When set, it wins: every field except
        ``filename`` is ignored. Each row is three numbers: stage (length),
        stored volume (length cubed) and surface area (length squared), in
        increasing stage -- a list of tuples, a 2-D ``numpy`` array, or a
        ``pandas.DataFrame``, which is read by its columns named ``stage``,
        ``volume`` and ``area`` (any letter case) when it has all three, else
        by position when it has exactly three columns, and refused
        (``ValueError``) otherwise. A bare DataFrame used to fail on its column
        LABELS (ledger 318). The rows go to :class:`LakeTable` as they are, so
        ``storage_coefficient`` is NOT applied, and :class:`LakeTable`'s
        validation -- MF6's own table checks, plus finite and non-negative
        values -- is the whole check.
    area : float, optional
        Default ``None``. **Prism mode.** The lake's constant surface area
        (length squared). Must be > 0. Every row gets this area, including the
        first row, at ``bottom``.
    top : float, optional
        Default ``None``. **Prism mode.** The highest stage (length). Must be >
        ``bottom``. The last row's stage is ``top`` (to ``numpy.isclose``).
    bottom : float, optional
        Default ``None``. **Prism mode.** The lowest stage (length), where the
        volume is zero. The first row is ``(bottom, 0.0, area)``.
    stage_step : float, default 1.0
        **Prism mode only.** The DEM mode takes its stages from ``stages`` and
        ignores this. It is the row spacing (length) and must be > 0. Rows sit at
        ``bottom, bottom + stage_step, ...`` below ``top``. ``top`` is added
        unless the last of those already equals it (``numpy.isclose``), so the
        last interval can be shorter than ``stage_step``: ``bottom=9, top=10,
        stage_step=0.4`` gives stages 9.0, 9.4, 9.8 and 10.0. A step larger than
        ``top - bottom`` gives just the two rows.
    storage_coefficient : float, default 1.0
        **Prism and DEM modes. The rows mode ignores it.** Multiplies every
        VOLUME, never the area. It is the fraction of the geometric volume that
        holds water, e.g. the porosity of a gravel-filled infiltration gallery,
        so :meth:`build` refuses anything but a number greater than 0 and at
        most 1 (``ValueError``, ``NaN`` included; ``TypeError`` for a ``bool``
        or a non-number). Each of those used to build (ledger 318): ``0`` a
        table of zero volumes, which MF6 refuses; a negative value falling
        volumes, refused by a message that did not name this field; and a value
        above 1 a basin holding more water than its own volume, which no fill
        does. MF6 knows nothing of it beyond the table: above the last row it
        grows the volume at the FULL last-row area (Notes, *Beyond the table*).
    dem : str or pathlib.Path, optional
        Default ``None``. **DEM mode.** Setting it selects the mode, unless
        ``rows`` is also set. It is a raster of the lake-BED elevation, anything
        :func:`rasterio.open` opens, and only band 1 is read. A
        :class:`~myflopy.surfaces.Surface` does NOT work here (``TypeError:
        invalid path or file``). Pass ``surface.resolve_source()`` instead, which
        returns:

        * for ``mf.Contours(...)``: the path of its interpolated GeoTIFF
          (``<contours>.interp.tif`` unless ``out=`` was given). GRASS runs on
          the first call and the cache is reused after that. A stale cache is
          reused with a warning; ``resolve_source(refresh=True)`` rebuilds it;
        * for ``mf.Raster(path)``: that path;
        * for ``mf.Flat(z)``: a float, which ``dem=`` refuses with ``TypeError``;
        * for any other Surface -- algebra such as ``.below()``,
          ``.capped_at()``, ``mf.Max`` or ``Surface.where``; a relative surface;
          ``Surface.from_points`` or ``Surface.from_array``; ``mf.Isopach`` --
          nothing, because it raises ``TypeError: Surface kind '...' has no
          raster/constant source``. No single raster sits behind those; write
          one yourself and pass its path.

        ``resolve_source()`` hands over the FILE, so the Surface's own ``units=``
        conversion and ``nodata=`` override do not apply here.

        Which pixels count: a pixel belongs to the lake when its CENTRE falls
        inside ``footprint`` (rasterio ``geometry_mask`` with
        ``all_touched=False``). A pixel holding the raster's DECLARED no-data
        value, or ``NaN``, has no bed and is left out; a ``NaN`` used to make
        every volume ``NaN`` (ledger 318). When SOME pixels inside the footprint
        are left out that way, :meth:`build` WARNS with their count and area,
        since the table is short by that area; when all of them are, it raises.
        A raster that declares no no-data value is read as it is: a ``0`` or
        ``-9999`` fill then counts as a real bed and inflates every volume.
        Each pixel's area is ``|a * e|`` from the raster's affine transform, in
        the raster CRS's units squared (rotation terms are ignored). A raster in
        a geographic CRS, whose pixels measure degrees, raises ``ValueError``.
    footprint : str, pathlib.Path, geopandas.GeoDataFrame, geopandas.GeoSeries, Area or shapely geometry, optional
        Default ``None``. **DEM mode, required** (``ValueError`` when missing).
        The lake outline that limits which DEM pixels count, read by the one
        footprint reader behind ``mf.lak``'s ``lakes`` and
        ``LayerStack.excavate``'s ``within``
        (``myflopy.modflow.mf6.grid.geometry.read_footprints``): a vector file
        as a ``str`` or ``pathlib.Path`` (or a list of them, stacked in order),
        read at ``footprint_layer``; a ``GeoDataFrame`` or ``GeoSeries``; an
        ``mf.Area``; or a shapely geometry. Anything else raises ``TypeError``.
        Every feature must be a ``Polygon`` or ``MultiPolygon`` with area
        (``ValueError``), and all of them are UNIONED into one outline, so a
        file holding several lakes gives one table for all of them together --
        pass the one lake's row, or put it in its own layer. It is REPROJECTED
        to the raster's CRS when both declare one; a footprint carrying no CRS
        (a bare shapely geometry, an ``mf.Area`` with ``crs=None``) is taken to
        be in the raster's frame already. A footprint that reaches NO pixel with
        a bed raises ``ValueError`` saying which of three things happened -- it
        lies off the raster, no pixel CENTRE falls inside it (it is narrower
        than a pixel), or every pixel inside is no-data -- with the extents and
        CRSs of both. More than one pixel's area of it lying off the raster
        WARNS, with that area. Until ledger 318 a ``str``, a frame and a series
        were refused, the footprint was never reprojected (a file's own CRS was
        ignored), only a file's first layer was read, and each of those misses
        built a table of zeros, which MF6 refuses.
    footprint_layer : str or int, optional
        Default ``None``. **DEM mode.** Which layer of a ``footprint`` FILE to
        read: its name, or its zero-based position in the file (any integer
        type, numpy's included; a ``bool`` raises ``TypeError``). ``None`` reads
        the first layer, GDAL's default -- for a GeoPackage holding several,
        behind pyogrio's generic "More than one layer found" warning. A name or
        position the file does not hold raises ``ValueError`` listing the
        layers present, and a ``footprint_layer`` given with a footprint that
        is not a file raises ``ValueError``, since it could not act. A position
        is accepted here, where ``mf.lak``'s ``layer=`` takes only a name,
        because that one sits beside a MODEL layer (``only_layer``) and this
        one beside nothing an integer could be mistaken for.
        ``geopandas.list_layers(path)`` gives the names.
    stages : sequence of float, optional
        Default ``None``. **DEM mode, required** (``ValueError`` when missing),
        and ignored by the other modes. These are the stages to tabulate
        (length), strictly increasing, one row each, e.g.
        ``np.arange(526.0, 536.25, 0.25)``. An empty sequence, a stage that is
        not finite, and unsorted or repeated stages raise ``ValueError``. At
        each stage, the volume is the sum over lake pixels of ``max(stage -
        bed, 0) * pixel area * storage_coefficient``. The area is the pixel area
        times the number of lake pixels whose bed is strictly below the stage
        (a pixel exactly AT the stage is dry). So every stage at or below the
        lowest bed gives ``(stage, 0.0, 0.0)``, and MF6 refuses two such rows
        (its volumes must strictly increase). **Only the LAST of them is
        kept** -- the highest stage at or below the lowest bed -- and the rest
        are dropped, with an INFO log line. That is not a guess: MF6 reads a
        stage below a table's first row as that row, and between two zero rows
        it reads zero, so the kept table gives MF6 the same volume and area at
        every stage as the full list would (measured with ``mf6``, ledger 318).
        The table's ``minimum_stage`` is then the kept stage, which is what
        ``mf.lak`` checks a ``starting_stage`` against. When EVERY stage is at
        or below the lowest bed the table would hold no water at any stage it
        lists, and ``ValueError`` names the lowest bed. Every stage above the
        highest bed gets the area of ALL lake pixels (their count x the pixel
        area, not the footprint polygon's own area), and the volume grows by
        that area x ``storage_coefficient`` per unit of stage. A first stage
        ABOVE the lowest bed makes MF6 give every lower stage the first row's
        volume (Notes), so start at or below it; the lowest bed is the minimum
        of the DEM inside the footprint, which need not be a contour label.
    footprint_buffer : float, default 0.0
        **DEM mode.** Buffers the footprint by this distance before masking, in
        the RASTER CRS's units -- after the footprint is reprojected there
        (shapely ``buffer``). A positive value grows it, e.g. to pick up pixel
        centres along a tightly drawn shoreline. A negative value shrinks it,
        and one that shrinks it to nothing raises ``ValueError``. ``0`` leaves
        it untouched.
    filename : str, optional
        Default ``None``. **All modes.** Carried onto the built
        :class:`LakeTable`. It names the ``laktab`` file ``mf.lak`` writes,
        beside the model's other input files. ``None`` or ``""`` lets ``mf.lak``
        name it ``"<package name>_<lake_id>.lak.tab"``, e.g.
        ``lak_pond.lak.tab`` for lake ``pond`` in the default package ``"lak"``.

    See Also
    --------
    LakeTable : The validated, immutable result of :meth:`build`.
    mf.lak : Takes it through ``tables={lake_id: ...}``.
    myflopy.surfaces.Surface.resolve_source : The raster path behind an
        ``mf.Contours`` or ``mf.Raster`` surface.

    Notes
    -----
    *Handing it to mf.lak.* Pass ``tables={lake_id: builder_or_table}``, with at
    most one table per lake. A key is read as a lake id the way the ids are
    made -- converted with ``str()`` -- so a lake whose ``lake_id_field`` value
    is the integer ``101`` may be keyed ``101`` or ``"101"``; both at once raise
    ``ValueError``. (An integer key used to raise ``unknown lake '101'``,
    printing the int exactly as the id it missed; ledger 318.) Without
    ``lake_id_field``, the ids are the file stems when every file holds one
    feature, and otherwise ``"lake_0"``, ``"lake_1"``, ... in row order. A key
    naming no lake raises ``ValueError: Lake table references unknown lake``,
    listing the ids there are. ``mf.lak`` calls :meth:`build`, then refuses a
    ``starting_stage`` outside ``[minimum_stage, maximum_stage]``. It writes
    each table as a child package named ``"<package name>_<lake_id>_table"``.
    Lakes left out of the mapping get no table.

    *What MF6 checks.* When it reads the file, MF6 needs the stages and the
    volumes STRICTLY increasing and the areas nondecreasing (``gwf-lak.f90``,
    ``lak_read_table``), and :class:`LakeTable` checks exactly that, so a table
    this builds is one MF6 reads. The DEM mode's areas never fall, and its
    volumes rise strictly once any pixel is wet; its zero-volume rows below the
    bed are reduced to one (``stages``), and a footprint reaching no pixel with
    a bed is refused (``footprint``). Before ledger 318 both of those built,
    and MF6 stopped with ``TABLE VOLUME ENTRY 2 ... MUST BE GREATER THAN THE
    PREVIOUS VOLUME ENTRY``.

    *Beyond the table.* Below the first stage, MF6 uses the first row's volume
    and area. Above the last stage, it keeps the last row's area and grows the
    volume at that FULL area; ``storage_coefficient`` is not applied there
    (``gwf-lak.f90``, ``lak_calculate_vol``). With ``storage_coefficient=1`` a
    DEM table need not go past the highest bed. With any other value, tabulate
    up to the highest stage the lake can reach, or MF6 treats the part above the
    table as open water. Measured: a 2 x 2 ft bed of 8, 9, 9 and 10 at
    ``storage_coefficient=0.5``, a lake held at 12.5 ft; the table tabulated to
    13 gave ``VOLUME`` 7.0 and the same table stopped at 11 gave 10.0, both
    under normal termination. Between rows MF6 interpolates volume and area
    linearly, so the stage spacing is the resolution of the curve.

    Examples
    --------
    A prism: an infiltration trench 40 x 10 ft and 6 ft deep, tabulated every
    0.5 ft, whose gravel fill holds 30% water:

    >>> trench = mf.LakeTableBuilder(area=400.0, bottom=96.0, top=102.0,
    ...                              stage_step=0.5, storage_coefficient=0.3)
    >>> trench.build().rows[:2]
    ((96.0, 0.0, 400.0), (96.5, 60.0, 400.0))

    Rows you already have, e.g. from a survey or a spreadsheet -- a DataFrame is
    read by its ``stage``/``volume``/``area`` columns:

    >>> import pandas as pd
    >>> survey = pd.DataFrame({"stage": [96.0, 98.0, 100.0],
    ...                        "volume": [0.0, 1.5e4, 4.2e4],
    ...                        "area": [0.0, 1.2e4, 1.6e4]})
    >>> mf.LakeTableBuilder(rows=survey).build().rows
    ((96.0, 0.0, 0.0), (98.0, 15000.0, 12000.0), (100.0, 42000.0, 16000.0))

    A bathymetry raster, masked to the lake outline -- the ``ponds`` layer of a
    GeoPackage, reprojected to the raster's CRS if it is in another:

    >>> mf.LakeTableBuilder(dem="bathy.tif", footprint="lakes.gpkg",
    ...                     footprint_layer="ponds", stages=np.arange(94.0, 102.25, 0.25))

    Bathymetry drawn as contours. The same ``mf.Contours`` can also be
    ``mf.lak``'s ``lake_bottom``, and ``resolve_source()`` hands this class its
    interpolated raster (GRASS runs the first time):

    >>> lak_bathy = mf.Contours("pond_contours.gpkg", z="Elev", epsg="2927",
    ...                         region_vector="pond.gpkg")
    >>> pond_table = mf.LakeTableBuilder(
    ...     dem=lak_bathy.resolve_source(),
    ...     footprint="pond.gpkg",
    ...     stages=np.arange(526.0, 536.25, 0.25),  # from the lowest bed up
    ... )

    Handing the tables to ``mf.lak``, keyed by lake id (the ``lake_id_field``
    values):

    >>> lak = mf.lak(context=ctx, nper=nper, lakes="pond.gpkg", lake_id_field="name",
    ...              starting_stage={"sw": 530.0}, lake_bottom=lak_bathy,
    ...              tables={"sw": pond_table})   # writes lak_sw.lak.tab
    """

    rows: Sequence[Sequence[float]] | pd.DataFrame | None = None
    area: float | None = None
    top: float | None = None
    bottom: float | None = None
    stage_step: float = 1.0
    storage_coefficient: float = 1.0
    dem: str | Path | None = None
    footprint: str | Path | Sequence[str | Path] | gpd.GeoDataFrame | gpd.GeoSeries | Area | BaseGeometry | None = None
    footprint_layer: str | int | None = None
    stages: Sequence[float] | None = None
    footprint_buffer: float = 0.0
    filename: str | None = None

    def _coefficient(self) -> float:
        """``storage_coefficient`` as a float in (0, 1], or a refusal naming the field (ledger 318)."""

        value = self.storage_coefficient
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
            raise TypeError(f"storage_coefficient must be a number, not {value!r}.")
        value = float(value)
        if not 0.0 < value <= 1.0:
            raise ValueError(
                f"storage_coefficient must be greater than 0 and at most 1 -- the fraction of "
                f"the geometric volume that holds water -- not {value!r}. 0 builds a table of "
                "zero volumes, which MF6 refuses; a negative value makes the volumes fall; and "
                "above 1 the basin would hold more water than its own volume."
            )
        return value

    def _dem_stages(self) -> list[float]:
        """``stages`` as python floats: present, not empty, finite and strictly increasing."""

        if self.stages is None:
            raise ValueError("DEM-derived lake tables require stages.")
        stages = [float(stage) for stage in self.stages]
        if not stages:
            raise ValueError("stages is empty: a DEM-derived lake table needs at least one stage.")
        for stage in stages:
            if not np.isfinite(stage):
                raise ValueError(f"stages must be finite; got {stage!r}.")
        for previous, current in zip(stages, stages[1:], strict=False):
            if not current > previous:
                raise ValueError(
                    "Lake table stages must be strictly increasing (one row each): "
                    f"{current!r} follows {previous!r}."
                )
        return stages

    def _footprint_frame(self, crs) -> gpd.GeoDataFrame:
        """``footprint`` read by the shared footprint reader, in ``crs`` (the raster's).

        :func:`~myflopy.modflow.mf6.grid.geometry.read_footprints` -- the reader
        ``mf.lak(lakes=)`` and ``LayerStack.excavate(within=)`` use -- so one file is
        read the same way, at the same layer and in the same frame, by all three. It
        refuses a feature that is not polygonal or has no area, and reprojects to
        ``crs``; a footprint carrying no CRS is taken to be in it already.
        """

        footprint = self.footprint
        is_file = isinstance(footprint, (str, Path)) or (
            isinstance(footprint, (list, tuple))
            and len(footprint) > 0
            and all(isinstance(item, (str, Path)) for item in footprint)
        )
        if not is_file and self.footprint_layer is not None:
            raise ValueError(
                f"footprint_layer={self.footprint_layer!r} selects a layer of a footprint "
                f"FILE, and the footprint is a {type(footprint).__name__}, which has none -- "
                "drop footprint_layer=, or pass the path."
            )
        if is_file or isinstance(footprint, gpd.GeoDataFrame):
            source = footprint
        elif isinstance(footprint, gpd.GeoSeries):
            source = gpd.GeoDataFrame(geometry=footprint)
        elif isinstance(footprint, Area):
            source = gpd.GeoDataFrame(geometry=[footprint.geometry], crs=footprint.crs)
        elif isinstance(footprint, BaseGeometry):
            source = gpd.GeoDataFrame(geometry=[footprint])
        else:
            raise TypeError(
                f"footprint got a {type(footprint).__name__}; pass a vector file path (str "
                "or Path), a GeoDataFrame or GeoSeries, an mf.Area, or a shapely polygon."
            )
        layer = self.footprint_layer
        if isinstance(layer, str):
            layer = str(layer)  # a numpy string is stored as a plain one (ledger 265)
        return read_footprints(source, layer=layer, crs=crs, noun="footprint")

    def _build_from_dem(self) -> LakeTable:
        """Integrate stage/volume/area from a DEM masked to the footprint at each stage.

        For every stage, sums inundation depth x PIXEL area over the valid pixels
        whose centres fall in the footprint (times the storage coefficient) for
        the volume, and the area of the pixels with depth > 0 for the surface
        area. Raster pixels, not model cells: the grid plays no part. The footprint
        is read and reprojected by the shared reader; a pixel with no bed (declared
        no-data or ``NaN``) is left out; a footprint reaching no pixel with a bed is
        refused, and a partial miss warns. The zero-volume rows at or below the
        lowest bed are reduced to the last of them (ledger 318).
        """

        stages = self._dem_stages()
        coefficient = self._coefficient()
        if self.footprint is None:
            raise ValueError(
                "DEM-derived lake tables require a footprint: a vector file path, a "
                "GeoDataFrame or GeoSeries, an mf.Area, or a shapely polygon."
            )
        with rasterio.open(self.dem) as source:
            raster_crs = CRS.from_user_input(source.crs) if source.crs else None
            name = Path(str(self.dem)).name
            if raster_crs is not None and raster_crs.is_geographic:
                raise ValueError(
                    f"{name} is in a geographic CRS ({raster_crs.to_string()}): its pixels "
                    "measure degrees, and a lake table's areas and volumes come from the "
                    "pixel size. Reproject the raster to the model's projected CRS first."
                )
            frame = self._footprint_frame(raster_crs)
            geometry = frame.geometry.union_all()
            if self.footprint_buffer:
                geometry = geometry.buffer(self.footprint_buffer)
            if geometry.is_empty:
                raise ValueError(
                    f"footprint_buffer={self.footprint_buffer!r} shrinks the footprint to "
                    "nothing: there are no pixels left to tabulate."
                )
            dem = source.read(1, masked=True)
            inside = geometry_mask(
                [geometry],
                out_shape=dem.shape,
                transform=source.transform,
                invert=True,
            )
            elevations = np.asarray(dem.astype(float).filled(np.nan), dtype=float)
            # A NaN pixel is no bed either: it made every volume NaN (ledger 318).
            valid = inside & ~np.ma.getmaskarray(dem) & np.isfinite(elevations)
            cell_area = abs(float(source.transform.a * source.transform.e))
            pixel = f"{abs(source.transform.a):g} x {abs(source.transform.e):g}"
            raster_box = box(*source.bounds)

        where = _footprint_and_raster(
            geometry, raster_box, source_crs=frame.attrs.get("source_crs"), raster_crs=raster_crs
        )
        inside_count = int(inside.sum())
        valid_count = int(valid.sum())
        if valid_count == 0:
            if not geometry.intersects(raster_box):
                why = "it lies off the raster"
            elif inside_count == 0:
                why = (
                    f"no pixel CENTRE falls inside it -- it is narrower there than the "
                    f"raster's {pixel} pixels; footprint_buffer= grows it"
                )
            else:
                why = (
                    f"all {inside_count} pixels whose centres fall inside it have no bed "
                    "(the raster's declared no-data value, or NaN)"
                )
            raise ValueError(
                f"The lake-table footprint reaches no pixel of {name} with a bed: {why}. "
                f"Every row would be (stage, 0.0, 0.0), which MF6 refuses. {where}"
            )
        if inside_count > valid_count:
            missing = inside_count - valid_count
            _warn_table(
                f"{missing} of the {inside_count} pixels of {name} whose centres fall inside "
                "the lake-table footprint have no bed (the raster's declared no-data value, "
                f"or NaN) and are left out, so the table's area and volume are short by up "
                f"to {missing * cell_area:g} (raster units squared) of lake."
            )
        off = geometry.difference(raster_box).area
        if off > cell_area:
            _warn_table(
                f"{off:g} of the lake-table footprint's {geometry.area:g} (raster units "
                f"squared) lies off {name}, more than one pixel: the table covers only the "
                f"part on the raster. {where}"
            )

        rows = []
        for stage in stages:
            depth = np.where(valid, np.maximum(stage - elevations, 0.0), 0.0)
            wet = depth > 0.0
            rows.append(
                (
                    float(stage),
                    float(depth.sum() * cell_area * coefficient),
                    float(wet.sum() * cell_area),
                )
            )
        # Every stage at or below the lowest bed gives (stage, 0.0, 0.0), and MF6 refuses
        # two equal volumes. They are a prefix (volume is 0 exactly until a pixel is wet),
        # and MF6 reads a stage below the first row AS the first row, so keeping only the
        # last of them gives MF6 the same curve at every stage (ledger 318, measured).
        dry = 0
        while dry < len(rows) and rows[dry][1] == 0.0:
            dry += 1
        if dry == len(rows):
            lowest = float(elevations[valid].min())
            raise ValueError(
                f"Every stage ({stages[0]:g} to {stages[-1]:g}) is at or below the lowest bed "
                f"inside the footprint, {lowest:g}: the table would hold no water at any stage "
                "it lists. Tabulate from the lowest bed up."
            )
        if dry > 1:
            logger.info(
                "LakeTableBuilder: %d stages at or below the lowest bed of %s are zero-volume "
                "rows, which MF6 refuses; kept only the last, %g -- MF6 reads a lower stage as "
                "that row, so the curve is unchanged.",
                dry,
                name,
                rows[dry - 1][0],
            )
            rows = rows[dry - 1:]
        return LakeTable(tuple(rows), filename=self.filename)

    def build(self) -> LakeTable:
        """Return the validated :class:`LakeTable` for the mode the fields select.

        The modes are checked in order, and the first that applies wins: ``rows``
        if set, else ``dem`` if set, else the prism (``area``/``top``/``bottom``).
        Only that mode's fields are read. ``mf.lak`` calls this for you. Call it
        yourself to look at the rows before building the package.

        Returns
        -------
        LakeTable
            The rows, carrying this builder's ``filename``. They pass MF6's own
            table checks (:class:`LakeTable`).

        Raises
        ------
        ValueError
            Rows mode: any :class:`LakeTable` rule, a string ``float()`` cannot
            convert, or a ``DataFrame`` whose columns cannot be told apart. DEM
            mode, in this order: ``stages`` missing, empty, not finite or not
            strictly increasing; ``storage_coefficient`` not in (0, 1];
            ``footprint`` missing; the raster in a geographic CRS;
            ``footprint_layer`` naming no layer of the file, or given with a
            footprint that is not a file; a footprint feature that is not a
            polygon or has no area; a ``footprint_buffer`` that shrinks it to
            nothing; a footprint reaching no pixel with a bed; every stage at or
            below the lowest bed. Prism mode: ``area``, ``top`` or ``bottom``
            missing (the message lists what each of the three modes needs),
            ``area <= 0``, ``top <= bottom``, ``stage_step <= 0``, or
            ``storage_coefficient`` not in (0, 1].
        TypeError
            Rows mode: a row that is not iterable, or ``rows`` a mapping or a
            string. DEM mode: ``dem`` is not something rasterio can open, e.g. a
            ``Surface`` itself, or the float that ``resolve_source()`` returns
            for ``mf.Flat``; a ``footprint`` of another type; a ``bool``
            ``footprint_layer``. Prism and DEM modes: a ``storage_coefficient``
            that is a ``bool`` or not a number.
        rasterio.errors.RasterioIOError
            DEM mode: ``dem`` names a file that is missing or not a raster.
        pyogrio.errors.DataSourceError
            DEM mode: ``footprint`` names a vector file that is missing or
            unreadable.

        Warns
        -----
        UserWarning
            DEM mode: some pixels inside the footprint have no bed (declared
            no-data or ``NaN``), or more than one pixel's area of the footprint
            lies off the raster -- the table is short by that area. Also logged.

        See Also
        --------
        LakeTable : What this returns.

        Examples
        --------
        >>> mf.LakeTableBuilder(area=100.0, bottom=9.0, top=10.0, stage_step=0.5).build().rows
        ((9.0, 0.0, 100.0), (9.5, 50.0, 100.0), (10.0, 100.0, 100.0))
        """

        if self.rows is not None:
            return LakeTable(self.rows, filename=self.filename)
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
        coefficient = self._coefficient()

        stages = list(np.arange(self.bottom, self.top, self.stage_step))
        if not stages or not np.isclose(stages[-1], self.top):
            stages.append(self.top)
        rows = tuple(
            (
                float(stage),
                float(self.area * (stage - self.bottom) * coefficient),
                float(self.area),
            )
            for stage in stages
        )
        return LakeTable(rows, filename=self.filename)


def _warn_table(message: str) -> None:
    """Say a lake table is short of lake, where a notebook sees it: logger AND ``warnings`` (ledger 288)."""

    logger.warning(message)
    warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())


def _footprint_and_raster(geometry, raster_box, *, source_crs, raster_crs) -> str:
    """Where a lake-table footprint is and where its raster is: two boxes and two frames (ledger 211).

    The raster-side twin of :func:`~myflopy.modflow.mf6.grid.geometry.where_they_are`,
    which words a GRID; the footprint's box is the one it was masked with, in the
    raster's frame after any reprojection.
    """

    def label(crs) -> str:
        return crs.to_string() if hasattr(crs, "to_string") else str(crs)

    if raster_crs is None:
        frame = "the raster declares no CRS, so the footprint's coordinates were used as they are"
    elif source_crs is None:
        frame = f"the footprint carries no CRS, so it was read as the raster's ({label(raster_crs)})"
    elif CRS.from_user_input(source_crs) != raster_crs:
        frame = f"the footprint was reprojected from {label(source_crs)} to the raster's {label(raster_crs)}"
    else:
        frame = f"both are in {label(raster_crs)}"
    span = ", ".join(f"{value:,.1f}" for value in geometry.bounds)
    raster = ", ".join(f"{value:,.1f}" for value in raster_box.bounds)
    return f"The footprint spans ({span}) and the raster ({raster}); {frame}."



def _record_count(value: Any) -> int | None:
    """How many rows a FloPy LAK block holds, or ``None`` when it is not rows in memory.

    A list, a tuple or an array is counted; anything else -- FloPy's external-file
    ``{"filename": ...}`` form, a DataFrame -- is not, and its count is left to FloPy
    and MF6.
    """

    return len(value) if isinstance(value, (list, tuple, np.ndarray)) else None


def _auxiliary_count(auxiliary: Any) -> int | None:
    """How many auxiliary NAMES ``auxiliary`` declares, or ``None`` when its shape is unknown.

    FloPy takes the names in several shapes and writes each one the same --
    ``["a", "b"]``, ``[("a", "b")]``, ``["auxiliary", "a", "b"]`` and the one-row
    recarray its own ``lak.auxiliary.get_data()`` returns, ``[("auxiliary", "a",
    "b")]``, all write ``AUXILIARY a b`` and MF6 runs them (measured). So the names
    are counted FLAT, through nested lists, tuples and record rows, and a leading
    ``auxiliary`` keyword is not a name. ``len(auxiliary)`` counted the nested and the
    recarray forms as ONE auxiliary, and :func:`_flopy_lak_refusals` then refused
    every packagedata row carrying both values as a boundname. An element that is not
    a string makes the shape unknown, and the caller then skips the check rather than
    guess.
    """

    if auxiliary is None:
        return 0
    names: list[str] = []

    def flatten(item: Any) -> bool:
        if isinstance(item, str):
            names.append(item)
            return True
        if isinstance(item, np.ndarray):
            return flatten(item.tolist())
        if isinstance(item, np.void):  # a record row: its fields, in order
            return flatten(item.tolist())
        if isinstance(item, (list, tuple)):
            return all(flatten(inner) for inner in item)
        return False

    if not flatten(auxiliary):
        return None
    if names and names[0].strip().lower() == "auxiliary":
        names = names[1:]
    return len(names)


def _flopy_lak_refusals(
    packagedata: Any,
    perioddata: Any,
    *,
    nlakes: int | None,
    noutlets: int,
    ntables: int,
    options: Mapping[str, Any],
) -> list[str]:
    """What ``mf.lak.flopy`` refuses before FloPy or MF6 would, one message each (ledger 318).

    The escape hatch keeps the caller's RECORDS as given; what it compares is each
    COUNT with the block it counts, and the shape of a period value -- every case here
    one FloPy fails at build or MF6 refuses by name when it reads the file, a run
    later (measured, ``docs/reproducers/lak_docstring_claims.py``, section 9). FloPy
    RAISES ``NOUTLETS``/``NTABLES`` to the rows given, so only a count ABOVE its rows
    is refused; ``nlakes`` it writes as given, so any disagreement is.
    """

    refusals = []
    lakes = _record_count(packagedata)
    if nlakes is not None and lakes is not None and int(nlakes) != lakes:
        refusals.append(
            f"nlakes={nlakes} but packagedata has {lakes} row(s); MF6 needs one row per lake "
            "and refuses either direction by name (NO DATA SPECIFIED FOR LAKE ..., lakeno MUST "
            "BE > 0 and <= ...). Leave nlakes out to count the rows."
        )
    for count_name, count, block in (("noutlets", noutlets, "outlets"), ("ntables", ntables, "tables")):
        given = options.get(block)
        rows = 0 if given is None else _record_count(given)
        if rows is not None and int(count) > rows:
            refusals.append(
                f"{count_name}={count} but {block} has {rows} row(s); FloPy raises the count to "
                "the rows given but never lowers it, and MF6 refuses a count above them by name."
            )
    naux = _auxiliary_count(options.get("auxiliary"))
    if lakes and naux is not None and not options.get("boundnames"):
        long_rows = [
            row for row in packagedata if isinstance(row, (list, tuple)) and len(row) > 3 + naux
        ]
        if long_rows:
            refusals.append(
                f"packagedata row {list(long_rows[0])} carries a value past [ifno, strt, "
                f"nlakeconn{', aux' * naux}] -- a boundname -- but boundnames=True is not "
                "set, and FloPy fails the build on it."
            )
    if isinstance(perioddata, Mapping):
        for period, records in perioddata.items():
            arrays = [
                record
                for record in (records if isinstance(records, (list, tuple)) else [])
                if isinstance(record, (list, tuple))
                and any(isinstance(item, (np.ndarray, pd.Series)) for item in record)
            ]
            if arrays:
                refusals.append(
                    f"perioddata[{period}] record {list(arrays[0][:2])} holds an array or "
                    "Series as its value: FloPy writes it whole, as text, and MF6 died reading "
                    "that (a segmentation fault). Give one number per record."
                )
    return refusals

def _build_lak_with_tables(model, *, _spec=None, lake_tables=(), implicit=False, **options):
    """Build LAK and its child table packages, after checking the lake against the DISV.

    When the model already carries its DISV, the lake's connections are compared with
    it BEFORE FloPy constructs the package (:func:`myflopy.grid_identity.
    lak_geometry_complaints`, ledger 325): FloPy validates each connection's cellid
    against the idomain at construction and would otherwise stop first with an
    ``MFDataException`` naming only the cellid. ``_spec`` -- the ``PackageSpec`` this
    builds, passed by ``PackageSpec.build`` -- is where the lake rules are read from.
    A DISV built AFTER the lake is covered by ``ModelSpec.build``'s own check.

    ``implicit`` (``mf.lak(implicit=True)``, carried in the spec's options) builds a
    :class:`~myflopy.modflow.mf6.lak_implicit.ModflowGwflakImplicit` instead of FloPy's
    ``ModflowGwflak``: FloPy 3.10 has no keyword for MF6's LAK ``IMPLICIT``, so that
    class's ``write`` inserts it into the file after FloPy writes it (ledger 326). The
    flag is consumed here, so FloPy never sees it.
    """

    options = {
        key: value.format(model_name=model.name) if isinstance(value, str) else value
        for key, value in options.items()
    }
    rules = (getattr(_spec, "metadata", None) or {}).get(LAKE_RULES)
    disv = model.get_package("disv") if rules and hasattr(model, "get_package") else None
    if disv is not None and not isinstance(disv, list):
        arrays = _disv_arrays(disv)
        name = getattr(_spec, "name", "lak")
        complaints = (
            lak_geometry_complaints(*arrays, options.get("connectiondata"), rules, name=name)
            if arrays is not None
            else []
        )
        if complaints:
            raise ValueError(
                f"Lake package {name!r} does not fit the DISV it is being built with:\n  - "
                + "\n  - ".join(complaints)
            )
    package_class = ModflowGwflakImplicit if implicit else flopy.mf6.ModflowGwflak
    package = package_class(model, **options)
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
    :meth:`build`. ``mf.lak(**kw, **options)`` is
    ``LAKBuilder(**kw, options=options).build()`` plus the ``lak_index`` that
    ``mf.lak_connection`` reads, so every parameter below means what it means
    there, with the same entry (a test holds the two texts equal). The one
    difference is TIMING: construction checks only ``context``, ``nper``,
    ``layer``, ``unexcavated`` and ``implicit``, and that ``tables`` reads as a dict
    whose keys name distinct lakes and ``outlets`` as a tuple; every other
    refusal an entry names is raised when :meth:`build` -- or the property that
    first reads the value (:attr:`lake_table`, :attr:`lake_cells`,
    :attr:`connectiondata`, ...) -- runs.

    **Per-lake values.** ``starting_stage``, ``bed_leakance``,
    ``lake_bottom``, ``lake_top``, ``only_vertical``, ``only_layer`` and
    ``rectangular_interior`` take one value per lake. Where an entry below
    says "any per-lake form", that is: one number for every lake; the name
    of a column of ``lakes``; a numeric string (one number for every lake);
    ``{lake_id: value}`` naming EVERY lake (a missing lake raises; a key is
    read through ``str()``, as the ids are made, so ``101`` names lake
    ``"101"`` -- it used to miss it, and ``only_vertical={101: True}`` was
    ignored silently, ledger 318); a ``pandas.Series`` whose index names
    every lake, read by label, a label naming a lake as a key does (an
    integer-indexed Series used to be read by POSITION, so one out of lake
    order gave each lake another's value; two labels naming one lake raise
    ``ValueError``), except that for ``lake_bottom`` and ``lake_top`` an
    integer-indexed Series that ALSO holds every cell of every lake is read
    by cell, its integers being cell numbers too; or any other sequence of
    one value per lake, in lake order (``LAKBuilder.lake_ids``).

    **Per-period values.** ``stage``, ``rainfall``, ``evaporation``,
    ``runoff``, ``withdrawals``, ``inflow`` and ``status`` each become one MF6
    period setting per lake per period, over ``range(nper)``. Periods are
    ZERO-based. The forms:

    * a scalar: every lake, every period;
    * a list or tuple, or a one-dimensional ``numpy`` array: every lake, by
      position (position 0 = the first stress period; ``None`` skips one);
    * ``{period: value}`` (every key an int): every lake, in the periods
      named;
    * ``{period: {lake_id: value}}``: the lakes named, in the periods named;
    * ``{lake_id: value}``: that lake, every period;
    * ``{lake_id: list}`` or ``{lake_id: {period: value}}``: that lake, per
      period;
    * a ``pandas.Series``: the mapping of its index, so an integer index is
      periods and a lake-id index is lakes.

    A period a value does not reach writes nothing for that lake, and MF6
    KEEPS the previous period's setting (measured: rainfall given for the
    first period only was still applied in the second). A mapping whose keys
    are ALL integers is read as periods and any other by lake id, so with
    integer-like lake ids key lakes by the STRING id. What would never be
    written, or written to a guess, raises ``ValueError``: a period outside
    ``range(nper)`` -- ``rainfall={101: 0.5}`` over lake ``"101"`` wrote
    nothing at all, silently -- a list longer than ``nper``, and an
    all-integer mapping whose every key is ALSO a lake id (give periods as
    a list, or as ``{period: {lake_id: value}}``). A ``numpy`` array or
    ``pandas`` Series was written WHOLE into every period, and MF6 died
    reading it (a segmentation fault); both are read as above now, and an
    array of two or more dimensions raises (ledger 318). Each value is
    checked here as MF6 would check it one run later: a ``status`` other
    than ``ACTIVE``, ``INACTIVE`` or ``CONSTANT`` (any case), a number that
    is not finite, and a negative rainfall, evaporation, runoff, inflow or
    withdrawal all raise ``ValueError``; a string elsewhere is written as
    given, which is how MF6 names a time series. MF6 refused the status and
    each negative forcing by name (measured, all five).

    Parameters
    ----------
    context : ModelContext
        Required (``None`` raises ``AttributeError``). The geometry the lakes
        resolve against. Its ``grid`` is required (``ValueError`` when
        ``None``); the lake cells, their neighbours and the connection
        lengths and widths come from it. The layer elevations come from
        ``surfaces`` when given (a table or a
        ``LayerBuildResult``), else from ``grid.gdf_topbtm``, and one of the
        two must exist (``ValueError``). The EXCAVATION RECORDS that decide
        whether a lake is written by the excavated rule (``connection_modes``)
        come from the same place, never from the geometry: a
        ``LayerBuildResult`` given as ``surfaces`` carries its
        ``excavations``; with no ``surfaces`` the grid carries the ones
        ``LayerStack.build(..., attach=True)`` published beside
        ``gdf_topbtm``; a bare table carries none, unless it IS the grid's
        own ``gdf_topbtm``. ``domain`` is the idomain; ``None``
        means every cell is active. ``length_units``/``time_units`` set the
        conversions when ``length_conversion``/``time_conversion`` are
        ``None``.
    nper : int
        Required. The number of stress periods the per-period values are
        expanded over, from zero. Must be at least 1 (``ValueError``). Give
        the simulation's TDIS count; the two are not compared here.
    lakes : str, pathlib.Path, sequence of str or Path, or geopandas.GeoDataFrame
        Required. The lake polygons, ``Polygon`` or ``MultiPolygon`` only
        (``ValueError`` otherwise, and for a polygon with no area). A path or
        a list of paths is read at ``layer`` (each file's first layer by
        default), and the files are stacked in order; files in different
        CRSs raise ``ValueError: Cannot determine common CRS``. The lakes are
        reprojected to the grid's CRS when both declare one; lakes carrying
        no CRS are read as the grid's. Lake order (MF6's lake numbers, from
        zero) is row order.
    layer : str, optional
        Default ``None``. Which layer of each ``lakes`` FILE to read, by NAME
        (a numpy string is converted with ``str``); ``None`` reads the first,
        which is all ``mf.lak`` could read before this parameter existed (a
        multi-layer GeoPackage reads its first layer behind pyogrio's generic
        warning). A name the file does not hold raises ``ValueError`` listing
        the layers present; so does a ``layer`` with a GeoDataFrame, which has
        none to choose. Read by the shared footprint reader
        (``myflopy.modflow.mf6.grid.geometry.read_footprints``). This is the
        VECTOR layer, not a model layer -- that is ``only_layer`` -- and so an
        integer is refused (``TypeError``) rather than read as a position in
        the file: ``layer=0`` reads as a model layer, and as the file's first
        layer it would have been the default, changing nothing.
        ``geopandas.list_layers(path)`` gives the names.
    lake_id_field : str, optional
        Default ``None``. The attribute naming each lake; its values are
        converted with ``str()`` to become the lake ids (``101`` becomes
        ``"101"``, a float ``101.0`` becomes ``"101.0"``). A name that is not
        a column raises ``ValueError``, and so do repeated ids. ``None`` names
        each lake by its FILE stem when ``lakes`` is paths and every stem is
        unique (one feature per file), and otherwise ``"lake_0"``,
        ``"lake_1"``, ... in row order. Every ``{lake_id: ...}`` mapping,
        ``tables`` key and :class:`~myflopy.modflow.mf6.lakes.LakeOutlet`
        uses these ids.
    starting_stage : float, str, dict of str to float, pandas.Series or sequence of float
        Required (``ValueError: starting_stage is required``). Each lake's
        stage at the start of the simulation (length), in any per-lake form.
        Refused below the lowest ``lake_bottom`` over a generated lake's
        cells, and outside a lake's table (``[minimum_stage,
        maximum_stage]``).
    lake_bottom : float, str, pathlib.Path, Surface, pandas.Series or dict, optional
        Default ``None``. The lakebed elevation (length). For an EXCAVATED
        lake (``connection_modes``) ``None`` takes the excavation's own bed --
        the bed is declared once, on ``LayerStack.excavate`` -- and a
        ``{lake_id: ...}`` mapping that leaves the lake out does the same. A
        value GIVEN for an excavated lake must be that bed: it is sampled as
        below and compared cell by cell, before any capping, and a difference
        beyond ``max(1e-6 * |bed|, 1e-6)`` at any cell raises ``ValueError``
        naming the worst (the stack was cut to another bottom). For every
        other lake it is required by the two generated modes (``ValueError``),
        and it is not read for a lake with explicit connections.
        One value per lake, in any per-lake form, is a flat bed -- what
        ``"rectangular"`` lakes take. ``"bathy"`` wants one per CELL:

        * a :class:`~myflopy.surfaces.Surface` -- ``mf.Raster("bathy.tif")``,
          ``mf.Contours(...)``, or algebra on one (``mf.Raster(dem).below(10)``);
        * a raster path, as a ``Path`` or a ``str`` (a ``str`` naming a column
          of ``lakes`` is read as the column);
        * a cell-indexed ``pd.Series`` (its index is cell numbers, not lake ids).

        Bare, it serves every lake; ``{lake_id: ...}`` gives each its own. A
        Surface or raster is sampled onto the grid once, area-weighted, and
        each lake reads its cells from that. Only the lake's own cells need a
        value: the neighbours ``"bathy"`` compares against are land, and sit at
        the model top wherever the source has none -- a raster clipped to the
        lake, or a Series of the lake's cells. A per-cell bed above its cell's
        model top is capped there, with a warning. A lake cell a raster has no
        data for (one the polygon clips at its edge) is set to the top, with a
        warning; a lake cell missing from a Series raises. What the bed DOES
        in MF6 is under ``connection_modes``: it picks each VERTICAL
        connection's layer, through which MF6 then counts from that cell's
        top, and it is the bottom of each HORIZONTAL face ``"bathy"`` writes on
        the lake cell, through which MF6 exchanges water from the bed up (out
        of the lake, in a transient period, no more than arrives).
    lake_top : float, str, pathlib.Path, Surface, pandas.Series or dict, optional
        Default ``None``. The facility top (length), **required for
        ``"rectangular"``** lakes (``ValueError``) and refused at or below
        ``lake_bottom`` there. It is the top of a rectangular lake's
        HORIZONTAL faces (and the ``telev`` written on its VERTICAL ones,
        which MF6 ignores). Takes the same forms as ``lake_bottom``, but is
        read at the lake's cells only and never filled or capped. Ignored by
        ``"bathy"``.
    only_vertical : bool or dict of str to bool, default False
        ``True``, for every lake, or ``{lake_id: True}`` (a lake left out is
        ``False``) drops every HORIZONTAL connection of a generated lake, in
        either mode: a sealed vault whose only exchange is through its floor.
    only_layer : int, dict of str to int, or None, default None
        A zero-based layer index, for every lake or ``{lake_id: layer}``
        (a lake left out is unrestricted), to which a generated lake's
        connections are confined. ``None`` means every layer a face spans. In
        either generated mode a cell whose bed is not in that layer raises
        ``ValueError: ... does not intersect active cell`` (a ``"rectangular"``
        lake used to get NO vertical connection there, silently; ledger 318).
        On an EXCAVATED lake it must name the top active layer under every
        lake cell (``ValueError`` otherwise), and its walls are confined to
        that layer.
    rectangular_interior : bool or dict of str to bool, default False
        ``"rectangular"`` lakes only. ``True`` (or ``{lake_id: True}``) also
        connects each lake cell toward its neighbours INSIDE the lake,
        modeling a permeable-fill basin (e.g. clean gravel) that is
        hydraulically continuous with the aquifer through every face.
        ``False`` = edges only (a lined trench or vault). It is a physical
        choice about the basin, not a solver setting, although all-faces also
        spreads the lake-stage coupling, which steadies a small basin that
        drains empty. Refused (``ValueError``) on an EXCAVATED lake: its
        floor is open water, with no interior face to conduct through -- a
        gravel-filled basin is aquifer, and is not excavated.
    bed_leakance : float, str, dict, pandas.Series or sequence of float, default 1.0
        The lakebed leakance (1/T: bed hydraulic conductivity divided by bed
        thickness) written on every generated connection, in any per-lake
        form. Two more forms reach one connection TYPE: ``{lake_id:
        {"vertical": v, "horizontal": h}}``, and a key ``(lake_id,
        "vertical")`` or ``(lake_id, "horizontal")``, which wins over the
        lake's own entry. The resolved value must be greater than zero
        (``ValueError``). An explicit
        :class:`~myflopy.modflow.mf6.lakes.LakeConnection` that sets its own
        ``bed_leakance`` wins over all of these.
    connection_modes : str or dict of str to (str or sequence of LakeConnection), default "bathy"
        How each lake's connections are made. A string, in any case, serves
        every lake; ``{lake_id: ...}`` gives each its own and must name every
        lake (``ValueError``). A lake cell can carry a VERTICAL connection and
        HORIZONTAL ones together. What each mode WRITES:

        * ``"bathy"`` (alias ``"automatic"``) -- a natural lake over a bed
          that varies by cell. Per lake cell, ONE VERTICAL connection in the
          highest active layer whose top and bottom bracket that cell's bed,
          written with ``belev = telev = bed`` (a bed in no active layer
          raises ``ValueError``). Then, toward each face neighbour whose bed
          is HIGHER -- another lake cell, or land, whose bed is the per-cell
          source's value there or else the model top -- one HORIZONTAL
          connection per active layer the step spans, written on the LAKE
          CELL ITSELF (the connection's cell id is the lake cell, not the
          neighbour), with ``belev = max(bed, layer bottom)``, ``telev =
          min(neighbour's bed, layer top)``, ``CONNLEN`` the distance from
          the LAKE cell's DISV centre to the line of the shared face
          (``grid.cl12``, MF6's own DISV geometry; until ledger 328 it was the
          full centre-to-centre spacing, about twice that) and ``CONNWIDTH``
          the shared face's length. A flat bed has no steps, so a lake of more
          than one cell gets NO horizontal connection, with a warning.
        * ``"rectangular"`` -- a box: trench, vault, lined basin. Needs a flat
          ``lake_bottom`` and a ``lake_top``. Per lake cell, ONE VERTICAL
          connection in the active layer holding ``lake_bottom``, written
          ``belev = lake_bottom``, ``telev = lake_top``; a cell where
          ``lake_bottom`` falls in no active layer raises ``ValueError``, as in
          ``"bathy"`` (it used to get no VERTICAL connection, silently; ledger
          318). Then,
          toward each neighbour OUTSIDE the lake (every neighbour with
          ``rectangular_interior=True``), one HORIZONTAL connection per active
          layer of the lake cell between ``lake_bottom`` and ``lake_top``,
          again on the lake cell itself, with the same ``CONNLEN`` and
          ``CONNWIDTH``.
        * a sequence of :class:`~myflopy.modflow.mf6.lakes.LakeConnection`
          (dict form only) -- written exactly as given; ``lake_bottom``,
          ``lake_top``, ``only_vertical``, ``only_layer`` and
          ``rectangular_interior`` do not touch that lake.

        **An EXCAVATED lake is written by its own rule, in either generated
        mode.** A lake is excavated when an excavation record reaches
        ``mf.lak`` (see ``context``) and the lake's cells are EXACTLY that
        excavation's cells -- ``LayerStack.excavate(within=<these lakes=
        polygons>)``, which covers cells by the same rule. Chosen from the
        record, never from the geometry (``LAKBuilder.lake_rules`` says which
        rule each lake got). Per lake cell: ONE VERTICAL connection to the top
        ACTIVE cell under it, whose top IS the bed, so MF6 counts from the
        bed; then, toward each face neighbour whose model top -- its bed, for
        a lake cell; the ground, for land; capped at ``lake_top`` for
        ``"rectangular"`` -- is above this cell's bed, one HORIZONTAL
        connection per active layer of the NEIGHBOUR the wall spans, written
        INTO the neighbour's cell (the aquifer the water meets), with ``belev
        = max(bed, layer bottom)``, ``telev = min(wall top, layer top)``,
        ``CONNLEN`` the distance from the NEIGHBOUR's DISV centre to the line
        of the shared face and ``CONNWIDTH`` the face's length -- MF6's own
        DISV geometry, the number ``grid.cl12`` holds in the neighbour's row. A
        cell whose declared bed did not reach below the ground (a rim cell)
        connects at the ground. Measured, a pit lake below the water table in
        uniform flow against the same pit as a K = 1e4 inclusion: through-flow
        1.005 of the reference, where ``"bathy"`` on the un-excavated stack
        gives 0.52 (0.39 before ledger 328). A lake covering only
        part of an excavation, or cells of two, raises ``ValueError``
        whatever ``unexcavated`` says. Declare the model with NEWTON: without
        it MF6 carries nothing through a wall whose cell is dry (the build
        warns). In EVERY mode a HORIZONTAL face thinner than 1e-6 is not
        written: MF6 would read it as the whole cell.

        **MF6 IGNORES the BELEV and TELEV written for a VERTICAL connection**
        and uses the connected cell's TOP (and top + ``surfdep``) instead (MF6
        input guide: "Any value can be specified if CLAKTYPE is VERTICAL";
        ``gwf-lak.f90``). So a bed that lies INSIDE a cell, below its top,
        does not lower the VERTICAL connection: through it MF6 exchanges
        water, and without a table stores it, only from that cell's top up.
        Measured: one VERTICAL connection under a cell whose top is 8.0, lake
        held at 12.0, bed written as 3.0, 6.0 and 0.0; MF6 reported
        ``VOLUME`` 400.0 = 100 ft2 x (12 - 8) all three times. For a
        HORIZONTAL connection MF6 DOES read both, requires them inside the
        cell, and treats ``belev == telev`` as the whole cell -- and
        ``"bathy"`` writes the lake cell's own HORIZONTAL faces from the bed
        up, so through THEM MF6 exchanges water from the bed, below the cell
        top, while the lake holds no volume there. What leaves the lake that
        way is still bounded by what it holds: in a TRANSIENT period MF6 caps
        the water leaving a lake at what is available to it (its storage at
        the start of the step plus what arrives), so a lake below its cell
        top takes water in through such a face freely but lets out only what
        arrives. Measured on the same cell, its VERTICAL connection plus a
        HORIZONTAL face from 3.0 (the bed) to 8.0, lake held at 6.0, below
        the top: ``VOLUME`` 0.0 and the VERTICAL connection 0.0 every time;
        the face carried 3.46 ft3/d out of the lake in a steady period, 0.0
        in a transient one, 2.0 there with an ``INFLOW`` of 2.0, and 5.15
        INTO the lake with the aquifer at 7.0. So the bed decides which LAYER
        a VERTICAL connection lands in and where each HORIZONTAL face starts.
        The elevations written never depend on the (transient) stage: MF6
        wets and dries each face itself, every time step.
    unexcavated : {"warn", "raise", "ignore"}, default "warn"
        What to do about a generated lake NO excavation record reached, when
        MF6 will not see its bed. ``"warn"`` issues a ``UserWarning`` (and a
        log line), ``"raise"`` raises ``ValueError`` with the same text, and
        ``"ignore"`` says nothing; any other value raises ``ValueError`` when
        the builder is constructed. Two findings trigger it. VERTICAL
        connections in a cell whose TOP is above the lake's bed: MF6 counts
        from that top, so each exchanges nothing until the stage reaches it
        and, without a table, the storage between the bed and the top is
        missing; the text gives how many, the worst and the median gap, and
        how many of those tops are above the starting stage and the table's
        highest stage (on a real pond, 577 of 612, worst 12.58). And a lake
        that LOOKS excavated -- more than half its cells have a zero-thickness
        layer above an active cell whose top is the bed -- but carries no
        record, so it was written without walls into its neighbours. A lake whose bed IS
        its cells' tops (a USG export's) triggers neither. Checked when the lakes are
        validated: at the ``mf.lak`` call, or at ``LAKBuilder.build``/``validate``. It
        never softens a partial match with an excavation, which always raises. The
        remedy both findings name: cut the pond into the layering and re-declare
        ``mf.lak`` -- ``stack.excavate(bottom=, within=)`` on a ``LayerStack``, rebuilt
        with ``attach=True`` or its ``LayerBuildResult`` passed as
        ``ModelContext(surfaces=...)``; or, with no stack (a DISV built from a loaded
        run's arrays), ``dug = result.excavate(bottom=, within=)`` on the
        ``LayerBuildResult`` the DISV was built from -- a NEW result, to build the DISV
        and the context from (``mf.disv(dug)``, ``ModelContext(surfaces=dug,
        domain=dug.idomain)``).
    tables : dict of str or int to (LakeTable or LakeTableBuilder), default {}
        Default ``{}``: no tables (``None`` raises ``TypeError`` here; ``mf.lak``
        reads it as ``{}``). Stage / volume /
        surface-area tables, at most one per lake. A value is a
        :class:`~myflopy.modflow.mf6.lakes.LakeTable` or a
        :class:`~myflopy.modflow.mf6.lakes.LakeTableBuilder`, which is built
        here; anything else raises ``TypeError``. A key is read as a lake id is
        made, through ``str()``, so a lake whose ``lake_id_field`` value is the
        integer ``101`` may be keyed ``101`` or ``"101"``; both at once raise
        ``ValueError``. (An int key used to raise ``unknown lake '101'``,
        printing the int exactly as the id it missed; ledger 318.) A key naming
        no lake raises ``ValueError`` listing the lake ids. Without a table,
        MF6 computes a lake's volume
        and surface area from the plan areas and TOPS of the cells its
        VERTICAL connections reach, not from the bed elevations. A table's
        stages must bracket that lake's ``starting_stage`` (``ValueError``
        otherwise). Each table is written to its own file: its ``filename``,
        or else ``"<name>_<lake_id>.lak.tab"``. To build one from a
        bathymetry raster, from ``mf.Contours`` (via ``.resolve_source()``),
        as a constant-area prism or from explicit rows, see
        ``mf.LakeTableBuilder``.
    outlets : tuple or list of LakeOutlet, default ()
        Lake outlets. Anything in it that is not a
        :class:`~myflopy.modflow.mf6.lakes.LakeOutlet` raises ``ValueError``.
        Each outlet names its ``source`` lake by lake id, and its
        ``receiver`` lake by lake id or ``None`` to discharge out of the
        model. An unknown id raises ``ValueError``. Outlets are numbered by
        position, from zero. ``outlet_type`` picks the discharge law:
        ``"SPECIFIED"`` (the default) uses ``rate``; ``"WEIR"`` uses
        ``invert`` and ``width``; ``"MANNING"`` uses those plus ``roughness``
        and ``slope``. A ``rate`` (a scalar for every period, a list with
        position 0 = the first stress period, or ``{period: value}`` keyed by
        ZERO-based stress period, or an array or Series read as
        ``LakeOutlet.rate`` says) is written as a ``RATE`` period setting, and
        refused (``ValueError``) where it would never be written or is not a
        finite number. MF6 applies it only to a SPECIFIED outlet, where a NEGATIVE rate
        removes water from the lake and a POSITIVE rate ADDS water. ``()``:
        no outlets.
    stage : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None`` (nothing written). The lake stage (length), a
        per-period value. MF6 applies it ONLY to a lake whose ``status`` is
        ``"CONSTANT"``, and holds the lake at it. A number must be finite.
    rainfall : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None``. Rainfall rate on the lake surface (length per
        time), a per-period value; MF6 multiplies it by the lake's surface
        area. Must be >= 0 (``ValueError``; MF6 refuses it by name).
    evaporation : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None``. The MAXIMUM evaporation rate from the lake surface
        (length per time), a per-period value; MF6 lowers it when the lake
        cannot supply it. Must be >= 0 (``ValueError``; MF6 refuses it by
        name).
    runoff : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None``. Runoff into the lake (length cubed per time), a
        per-period value. Must be >= 0 (``ValueError``; MF6 refuses it by
        name).
    withdrawals : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None``. The MAXIMUM withdrawal from the lake (length cubed
        per time), a per-period value, written as MF6's ``WITHDRAWAL``. Must
        be >= 0 (``ValueError``; MF6 refuses it by name): it is a positive
        number that removes water.
    inflow : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None``. External inflow to the lake (length cubed per time),
        a per-period value. Must be >= 0 (``ValueError``; MF6 refuses it by
        name).
    status : str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
        Default ``None`` (MF6's default, ``ACTIVE``). Each lake's status, a
        per-period value: ``"ACTIVE"``, ``"INACTIVE"`` or ``"CONSTANT"`` (the
        last holds the lake at ``stage``), in any case, written as given (MF6
        accepts any case, measured with ``"constant"``). Anything else raises
        ``ValueError`` here; MF6 refused it by name one run later (measured:
        ``Unknown LAK lak status keyword: BOGUS``).
    name : str, default "lak"
        The package name. It is the FloPy ``pname``, the input file
        ``<model>.<name>``, the stem of the four output files
        (``<model>_<name>_stage.lak``, ``_budget.lak``, ``_budget.csv``,
        ``_package_convergence.csv``), the prefix of each default table file
        (``<name>_<lake_id>.lak.tab``) and table child (``<name>_<lake_id>_table``),
        and the package a mover endpoint names (:meth:`connection`). Give
        two LAK packages in one model different names.
    mover : bool, default False
        Writes MF6's ``MOVER`` option, so the lakes can give water to, and
        take it from, other packages through ``mf.mvr``; address a lake with
        :meth:`connection` (``builder.connection(lake_id)``). The spec
        :meth:`build` returns carries no ``lak_index``, so
        ``mf.lak_connection`` refuses it (``ValueError``); ``mf.lak`` adds one.
    boundnames : bool, default True
        Writes each lake's id as its MF6 boundname (and the ``BOUNDNAMES``
        option), so observations and budgets can name lakes by id. ``False``
        leaves both out.
    length_conversion : float, optional
        Default ``None``. MF6's ``LENGTH_CONVERSION``, from metres to model
        length units; MF6 uses it only for WEIR and MANNING outlet flows.
        ``None`` derives it from ``context.length_units``, read ignoring
        case and surrounding spaces: 3.28081 for ``"feet"``/``"foot"``/
        ``"ft"``/``"us-ft"``/``"usft"``, 1.0 for ``"meters"``/``"meter"``/
        ``"metres"``/``"m"``, 100.0 for ``"centimeters"``/``"centimeter"``/
        ``"cm"``, 1.0 for ``None``, ``""`` or ``"unknown"``.
        ``ModelContext`` refuses any other name at construction (ledger
        433).
    time_conversion : float, optional
        Default ``None``. MF6's ``TIME_CONVERSION``, from seconds to model
        time units; also used only for WEIR and MANNING outlets. ``None``
        derives it from ``context.time_units``, read the same way: 86400.0
        for ``"days"``/``"day"``/``"d"``, 3600.0 for ``"hours"``/``"hour"``/
        ``"hr"``/``"h"``, 60.0 for ``"minutes"``/``"minute"``/``"min"``, 1.0
        for ``"seconds"``/``"second"``/``"sec"``/``"s"``, 31557600.0 for
        ``"years"``/``"year"``/``"yr"``, 1.0 for ``None``, ``""`` or
        ``"unknown"``, and ``ValueError`` for any other name.
    maximum_iterations : int, default 100
        MF6's ``MAXIMUM_ITERATIONS``: the most Newton-Raphson iterations
        MF6 spends on each lake's stage per outer iteration. Raise it only
        when a lake's budget error is large.
    maximum_stage_change : float, default 1e-05
        MF6's ``MAXIMUM_STAGE_CHANGE``: the lake-stage closure tolerance
        (length).
    implicit : bool or numpy.bool_, default False
        Writes MF6's LAK ``IMPLICIT`` option: each lake's stage is solved as one more
        unknown of the groundwater-flow matrix, WITH the heads, instead of by MF6's
        default substitution iteration, which solves the stage and the heads by turns.
        That can take many outer iterations, or fail to converge, when a lake and the
        aquifer are strongly coupled -- a bed that leaks freely, an excavated pond
        exchanging through its floor and walls: on a real one the default formulation
        took 81% more outer iterations and 16 failed time steps, and ``IMPLICIT`` 11%
        and none (the excavation design's runs). The converged answer is the same: a
        transient pit filled at 2,000 ft3/d reached the default formulation's stages to
        1e-8 ft (measured). MF6's requirements, each checked here: the solution's linear
        acceleration must be BICGSTAB -- MF6 stops on CG, so an IMS resolving to CG
        (``linear_acceleration="CG"``, or ``complexity="SIMPLE"`` or no complexity with
        no ``linear_acceleration``) is refused when the simulation is built, and the IMS
        ``mf.simulation`` makes for you picks BICGSTAB -- for a lake reached through
        ``mf.ref`` too, once the build resolves it; every lake must have a connection
        (``ValueError``); NEWTON is recommended for a lake that may become perched, and
        the build warns without it. ``maximum_iterations`` and ``maximum_stage_change``
        then do nothing (MF6 checks the stage in IMS), and the package-convergence file
        is not asked for, since MF6 writes none. **Availability:** MF6's develop branch
        (commit ``e2dfc84c``) and the 6.8 line, in no release before; the executable is
        probed once when the run executes, and one without the option is WARNED about,
        naming it and its version -- the run goes ahead, and MF6 stops with its own
        message. FloPy 3.10 has no keyword for it, so myflopy writes it into the LAK
        file after FloPy does (``ModflowGwflakImplicit``). It survives a run loaded by
        ``Run.open`` or ``mf.load_mf6_run`` and written again,
        ``ParallelModelWorkflow.split_model`` and a package artifact; a bare
        ``flopy.mf6.MFSimulation.load`` drops it silently, and ``model.config`` shows it
        as the lake's ``implicit`` setting. ``True`` rides in the spec's options (it is
        saved and moves the digest); ``False`` writes and stores nothing. A
        ``numpy.bool_`` is accepted and stored as ``True``/``False``; anything else
        raises ``TypeError``.
    options : dict of str to Any, default {}
        ``mf.lak``'s ``**options`` collected into one mapping: any other
        ``flopy.mf6.ModflowGwflak`` keyword, e.g. ``surfdep``,
        ``print_stage``, ``print_flows``, ``observations`` or
        ``auxiliary``. It is read by :meth:`build`, not at construction, so a
        name FloPy does not accept raises ``TypeError`` (with a did-you-mean)
        there. These override what the builder sets by default --
        ``save_flows=True`` and the four output file records -- and a string
        value may contain ``{model_name}``, filled in at build. Rename
        the package with ``name``: a ``pname`` here raises ``TypeError`` at
        :meth:`build`, because it renamed only what FloPy writes while the
        spec's name, :meth:`connection` and each table child kept ``name``
        (ledger 318). A ``filename`` here names the LAK input file only. What
        the builder computes itself is not for here: ``nlakes`` raises
        ``TypeError`` naming it (it used to replace the lake count silently,
        ledger 318), and ``noutlets``, ``ntables`` and the three data blocks
        raise ``TypeError: ... got multiple values``, all at :meth:`build`.

    Raises
    ------
    ValueError
        At construction, ``context.grid`` is ``None``, ``nper`` is below 1,
        ``unexcavated`` is not one of its three words, or two ``tables`` keys
        name one lake (``101`` and ``"101"``). Every other refusal
        named above (the full list is ``mf.lak``'s) is raised by :meth:`build`
        or the property that first needs the value -- the excavated rule's by
        :attr:`lake_rules` and :attr:`connections`, the ``unexcavated="raise"``
        refusal by :meth:`validate`.
    AttributeError
        ``context=None``.
    TypeError
        At construction, an ``implicit`` that is not a bool, a ``layer`` that is
        not a ``str`` (an integer reads as a model layer; that is
        ``only_layer``), a ``tables`` that ``dict()`` or an
        ``outlets`` that ``tuple()`` cannot read (``None`` included). At :meth:`build`, an
        ``options`` name ``ModflowGwflak`` does not accept, one the builder
        sets itself, ``nlakes`` or ``pname``; a ``tables`` value that is
        neither a ``LakeTable`` nor a ``LakeTableBuilder``.

    See Also
    --------
    myflopy.lak : The package-first front door, ``mf.lak(...)``.
    LakeTableBuilder : Stage / volume / area tables for ``tables=``.
    LakeOutlet : Outlets for ``outlets=``.
    LakeConnection : Hand-written connections for ``connection_modes=``.

    Examples
    --------
    >>> builder = mf.LAKBuilder(context=ctx, nper=1, lakes="pond.gpkg",
    ...                         lake_id_field="name", starting_stage=101.0,
    ...                         lake_bottom=mf.Raster("bathy.tif"), mover=True)
    >>> builder.lake_cells
    {'pond': [5, 6, 9, 10]}
    >>> spec = builder.build()
    >>> builder.connection("pond")
    MoverConnection(package='lak', index=0)
    """

    context: ModelContext
    nper: int
    lakes: Path | str | Sequence[Path | str] | gpd.GeoDataFrame
    layer: str | None = None
    lake_id_field: str | None = None
    starting_stage: Any = None
    lake_bottom: Any = None
    lake_top: Any = None
    only_vertical: Any = False
    only_layer: Any = None
    rectangular_interior: Any = False
    bed_leakance: Any = 1.0
    connection_modes: str | Mapping[str, str | Sequence[LakeConnection]] = "bathy"
    unexcavated: str = "warn"
    tables: Mapping[str | int, LakeTable | LakeTableBuilder] = field(default_factory=dict)
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
    length_conversion: float | None = None  # None -> derive from context.length_units
    time_conversion: float | None = None      # None -> derive from context.time_units
    maximum_iterations: int = 100
    maximum_stage_change: float = 1.0e-5
    implicit: bool = False
    options: Mapping[str, Any] = field(default_factory=dict)
    _lakes: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False, compare=False)
    _connections: tuple[tuple[str, LakeConnection], ...] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    # lake_bottom / lake_top resolved once per lake, and each raster or Surface
    # sampled onto the grid once, however many lakes and cells read it.
    _per_lake: dict[tuple[str, str], float | dict[int, float]] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )
    _sampled: dict[Any, np.ndarray] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )
    _lake_cells: dict[str, list[int]] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    # The CRS the lakes arrived in, kept because ``lake_table`` reprojects them to
    # the grid's: a miss is diagnosed from the two frames, and after that
    # reprojection the frame on hand can only report the grid's own.
    _source_crs: Any = field(default=None, init=False, repr=False, compare=False)
    # The elevation table (geometry dropped) and the excavation records beside it, read
    # ONCE through the one resolver (`specs._layer_elevation_source`), so the rule a lake
    # is given and the elevations it is written against cannot come from two places.
    _source: tuple[pd.DataFrame | None, dict] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _elevations: np.ndarray | None = field(default=None, init=False, repr=False, compare=False)
    # {lake_id: "excavated" | "legacy" | "explicit"} and each excavated lake's record.
    _rules: dict[str, str] | None = field(default=None, init=False, repr=False, compare=False)
    _records_of: dict[str, Any] | None = field(default=None, init=False, repr=False, compare=False)
    # (lake_id, lake cell or None, role) beside each entry of `connections`.
    _roles: tuple[tuple[str, int | None, str], ...] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _face_geometry: DisvFaces | None = field(default=None, init=False, repr=False, compare=False)
    # `prepared_tables`' built tables, beside the `tables` entries they were built from.
    _prepared: tuple[tuple[tuple[str, Any], ...], dict[str, LakeTable]] | None = field(
        default=None, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        """Validate a grid is present, ``nper >= 1`` and ``unexcavated``, and freeze tables/outlets."""

        if self.context.grid is None:
            raise ValueError("ModelContext.grid is required for LAKBuilder.")
        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        if self.unexcavated not in _UNEXCAVATED_WORDS:
            raise ValueError(
                f"unexcavated must be one of {', '.join(repr(w) for w in _UNEXCAVATED_WORDS)}, "
                f"not {self.unexcavated!r}."
            )
        # `np.bool_` is what a numpy predicate returns and is NOT a `bool`; stored as
        # one so the spec's options serialize and digest as a plain `True` (ledger 265).
        if not isinstance(self.implicit, (bool, np.bool_)):
            raise TypeError(f"implicit must be True or False, not {self.implicit!r}.")
        object.__setattr__(self, "implicit", bool(self.implicit))
        if self.layer is not None:
            object.__setattr__(
                self, "layer", _vector_layer_name(self.layer, of="lakes", model_layer="only_layer=")
            )
        # A table key is read as a lake id is made, through `str()`: an int key used to
        # miss its lake and raise "unknown lake '101'" (ledger 318).
        tables: dict[str, Any] = {}
        for key, table in dict(self.tables).items():
            label = str(_plain_scalar(key))
            if label in tables:
                raise ValueError(
                    f"tables names lake {label!r} twice; lake ids are strings, so an integer "
                    "key names the lake of that string."
                )
            tables[label] = table
        object.__setattr__(self, "tables", tables)
        object.__setattr__(self, "outlets", tuple(self.outlets))

    def with_updates(self, **updates: Any) -> LAKBuilder:
        """Return a changed builder without modifying the original.

        The way to vary one input -- a stage, a bed, a leakance -- across
        scenarios without re-reading the lake polygons by hand.

        Parameters
        ----------
        **updates : Any
            Field names of this builder (its Parameters) and their new values.
            The copy is made with :func:`dataclasses.replace`, so it is
            constructed again: its caches start empty and the construction
            checks run again.

        Returns
        -------
        LAKBuilder
            A new builder; ``self`` is unchanged.

        Raises
        ------
        TypeError
            A name that is not a field, or one of the private caches.
        ValueError
            The construction checks (``context.grid`` ``None``, ``nper`` below 1).

        See Also
        --------
        build : Turn either builder into its package spec.

        Examples
        --------
        >>> deeper = builder.with_updates(starting_stage=102.0)
        >>> deeper.starting_stage, builder.starting_stage
        (102.0, 101.0)
        """

        return replace(self, **updates)

    @property
    def grid(self):
        """The Voronoi grid carried by the model context."""

        return self.context.grid

    @property
    def lake_table(self) -> gpd.GeoDataFrame:
        """Return normalized lake records with stable lake IDs.

        Read by :func:`~myflopy.modflow.mf6.grid.geometry.read_footprints`, the one
        reader ``LayerStack.excavate(within=)`` uses too, so a lake and the pit cut
        for it can never read different features, layers or frames of one file.
        """

        if self._lakes is not None:
            return self._lakes
        records = read_footprints(
            self.lakes,
            layer=self.layer,
            crs=getattr(self.grid, "crs", None),
            id_field=self.lake_id_field,
            noun="lake",
        )
        object.__setattr__(self, "_source_crs", records.attrs.get("source_crs"))
        records = records.copy()
        records["lake_id"] = records.index.to_numpy(dtype=object)
        records["lake_number"] = range(len(records))
        records = records.set_index("lake_id", drop=False)
        object.__setattr__(self, "_lakes", records)
        return records

    @property
    def lake_ids(self) -> tuple[str, ...]:
        """The stable string lake IDs in build order."""

        return tuple(self.lake_table.index)

    @property
    def lake_numbers(self) -> dict[str, int]:
        """Map each stable lake ID to its zero-based MF6 lake number."""

        return self.lake_table["lake_number"].astype(int).to_dict()

    def lake_number(self, lake_id: str) -> int:
        """The zero-based MF6 lake number for a stable ``lake_id``.

        MF6 numbers lakes by row order in ``lakes``; this reads that number back
        from the id, which is what an outlet, a table and a mover name.

        Parameters
        ----------
        lake_id : str or object
            A lake id (see ``lake_id_field``). It is converted with ``str()``
            first, so ``101`` finds the lake whose id is ``"101"``.

        Returns
        -------
        int
            The lake's position in :attr:`lake_ids`, MF6's lake number from zero.

        Raises
        ------
        ValueError
            ``Unknown lake ID: ...`` for an id no lake has.

        See Also
        --------
        connection : The same number as a mover endpoint.

        Examples
        --------
        >>> builder.lake_number("pond")
        0
        """

        try:
            return self.lake_numbers[str(lake_id)]
        except KeyError as error:
            raise ValueError(f"Unknown lake ID: {lake_id}") from error

    def connection(self, lake_id: str) -> MoverConnection:
        """Return an MVR endpoint for a lake by stable lake ID.

        The builder's counterpart of ``mf.lak_connection``, which reads the
        ``lak_index`` only ``mf.lak`` stamps on its spec.

        Parameters
        ----------
        lake_id : str or object
            A lake id, read as :meth:`lake_number` reads it.

        Returns
        -------
        MoverConnection
            ``MoverConnection(name, lake_number)`` for ``mf.Move``. The package
            must be declared with ``mover=True`` for MF6 to move water through it.

        Raises
        ------
        ValueError
            An unknown lake id.

        See Also
        --------
        lake_number : The lake number alone.

        Examples
        --------
        >>> builder.connection("pond")
        MoverConnection(package='lak', index=0)
        """

        from myflopy.modflow.mf6.mvr import MoverConnection

        return MoverConnection(self.name, self.lake_number(lake_id))

    @property
    def lake_cells(self) -> dict[str, list[int]]:
        """Map each lake ID to the grid cells its polygon COVERS (cached; intersection is costly).

        A lake covers a cell by AREA, through the shared rule
        :func:`myflopy.modflow.mf6.grid.geometry.cell_shares` -- the same one a
        ``.gpkg`` boundary, a head target and a region are resolved with. Every
        cell whose overlap was above zero used to count, and each of those becomes
        a full LAK cell with its own vertical connection, so a lake drawn exactly
        on one cell of a State Plane grid and round-tripped through EPSG:4326 grew
        from 1 cell to 4 -- the extras overlapping by 7.41e-06, 9.07e-14 and
        7.41e-06 square feet -- and one grown by 1e-06 ft reached all nine cells
        around it. That changes the lake's sidewall topology and its budget. A
        piece narrower than ``_SLIVER_WIDTH`` and under ``_SLIVER_SHARE`` of its
        cell is noise and goes, UNLESS every piece of the lake is that thin, when
        they are the lake. That measure is per PIECE, so a multipart lake whose
        detached part is itself narrower than the floor loses that part's cell
        (measured: a 15 x 0.005 ft satellite of a one-cell lake, 0.075 square feet,
        goes; a 1 x 1 ft pond and a 15 x 0.1 ft strip stay) -- the price of one
        rule for every caller, and the same answer a ``.gpkg`` zone of that shape
        gets.

        Shares are discarded here: LAK connects to a cell or it does not, and
        nothing about a connection is apportioned by how much of the lake sits in
        one. Which LAYER each cell connects in is still
        :meth:`_layer_containing`'s, against that cell's own elevations.
        """

        # Cached: each access does a geometry intersection per candidate cell, and
        # connection building reads this once per lake cell -- recomputing it there
        # made an automatic build over a large lake O(cells x lakes x grid) and hang.
        if self._lake_cells is not None:
            return self._lake_cells
        table = self.lake_table
        cells = footprint_cells(self.grid.gdf_vorPolys, table.geometry)
        result = dict(zip(table.index, cells, strict=True))
        object.__setattr__(self, "_lake_cells", result)
        return result

    def _lake_value(
        self, value: Any, lake_id: str, *, name: str, default: Any = None, per_cell: bool = False
    ) -> Any:
        """Resolve a per-lake parameter for ``lake_id`` from many input shapes.

        Accepts ``None`` (uses ``default``), a lake-table column name, a numeric
        string, a per-lake mapping, a scalar, a Series keyed by lake id (read by
        label, never position, a label through ``str()`` as a key is --
        :func:`_series_lake_labels`), or a one-per-lake sequence; raises if a
        required value is missing. With ``per_cell`` (``lake_bottom``/``lake_top``)
        it also passes through what varies by CELL for :meth:`_per_cell`: a
        :class:`~myflopy.surfaces.Surface`, a raster ``Path`` or ``str`` naming an
        existing file, or a Series whose index does not name every lake -- or one
        whose integer labels name every lake AND hold every lake cell, since those
        integers are cell numbers too.
        """

        if value is None:
            if default is None:
                raise ValueError(f"{name} is required.")
            return default
        if isinstance(value, str):
            if value in self.lake_table.columns:
                return self.lake_table.loc[lake_id, value]
            if per_cell and Path(value).is_file():
                return Path(value)
            try:
                return float(value)
            except ValueError:
                columns = ", ".join(str(c) for c in self.lake_table.columns if c != "geometry")
                raise ValueError(
                    f"{name}={value!r} is not a column of the lake table ({columns})"
                    + (" or an existing raster file." if per_cell else ".")
                ) from None
        if isinstance(value, Mapping):
            found, inner = _lake_keyed(value, lake_id, name=name)
            if not found:
                hint = (
                    " A mapping here is keyed by lake id; for one value per cell pass a "
                    "cell-indexed pandas Series, a raster, or a Surface."
                    if per_cell
                    else ""
                )
                raise ValueError(f"{name} is missing lake '{lake_id}'.{hint}")
            if per_cell and isinstance(inner, str) and Path(inner).is_file():
                return Path(inner)
            return inner
        if isinstance(value, Surface):
            if not per_cell:
                raise ValueError(f"{name} is one value per lake; only lake_bottom and lake_top take a Surface.")
            return value
        if isinstance(value, (Path, Real)):
            return value
        if isinstance(value, pd.Series):
            labels = _series_lake_labels(value, self.lake_ids, name=name)
            if (
                labels is not None
                and per_cell
                and not all(isinstance(label, str) for label in labels.values())
            ):
                # An integer is a cell number here as well as a lake id: a Series that
                # also holds every cell of every lake is the per-cell map it was always
                # read as (a whole-grid bed over lakes numbered 1, 2, ... holds both).
                index = set(value.index)
                if all(cell in index for lake in self.lake_ids for cell in self.lake_cells[lake]):
                    labels = None
            if labels is not None:
                return value.loc[labels[lake_id]]
            if per_cell:
                return value  # indexed by cell: the same map for every lake
        values = list(value)
        if len(values) != len(self.lake_ids):
            hint = (
                " For one value per CELL pass a cell-indexed pandas Series or "
                "mf.Surface.from_array(values)."
                if per_cell
                else ""
            )
            raise ValueError(f"{name} must contain one value per lake.{hint}")
        return values[self.lake_number(lake_id)]

    def _connection_mode(self, lake_id: str) -> str | Sequence[LakeConnection]:
        """The connection strategy for a lake: ``"bathy"``/``"rectangular"`` or explicit connections.

        Resolves a per-lake mapping, normalizes the legacy ``"automatic"`` alias to
        ``"bathy"``, and validates explicit connection lists.
        """

        mode = self.connection_modes
        if isinstance(mode, Mapping):
            found, mode = _lake_keyed(mode, lake_id, name="connection_modes")
            if not found:
                raise ValueError(f"connection_modes is missing lake '{lake_id}'.")
        if isinstance(mode, str):
            mode = mode.lower()
            if mode == "automatic":  # retained alias for the natural-lake mode
                mode = "bathy"
            if mode not in {"bathy", "rectangular"}:
                raise ValueError(
                    "Generated connection modes must be 'bathy' or 'rectangular' "
                    "(legacy alias 'automatic' == 'bathy')."
                )
            return mode
        if not all(isinstance(item, LakeConnection) for item in mode):
            raise ValueError("Explicit connection modes must contain LakeConnection objects.")
        return tuple(mode)

    def _only_vertical(self, lake_id: str) -> bool:
        """Sealed-sides (concrete bottomless vault) case: only vertical (bottom)
        connections, no horizontal exchange at all -- not even at the perimeter.

        Distinct from a ``rectangular`` edges-only lake (``rectangular_interior``
        ``False``), which still connects through its PERIMETER sidewalls (a lined
        trench open at the edges). Use ``only_vertical`` when the walls are
        impermeable on every side and only the bottom exchanges with the aquifer.
        """
        value = self.only_vertical
        if isinstance(value, Mapping):
            found, value = _lake_keyed(value, lake_id, name="only_vertical")
            return bool(value) if found else False
        return bool(value)

    def _only_layer(self, lake_id: str) -> int | None:
        """Restrict a lake's connections to a single model layer (special cases only)."""
        value = self.only_layer
        if isinstance(value, Mapping):
            return _lake_keyed(value, lake_id, name="only_layer")[1]
        return value

    def _rectangular_interior(self, lake_id: str) -> bool:
        """Whether a rectangular lake also connects through its INTERIOR faces.

        This is a PHYSICAL choice about the basin's construction, not a numerical
        or solver tuning knob -- pick it from whether the interior faces actually
        conduct water, never to nudge convergence.

        * ``False`` (default) = edges only. The flat bottom is one infiltration
          surface and only the perimeter cells expose a vertical sidewall to the
          aquifer -- a lined trench, or an infiltration basin with an impermeable
          liner between adjacent lake cells.
        * ``True`` (or a per-lake mapping) connects every lake cell to every
          neighbor, representing a basin whose fill is permeable (e.g. clean
          gravel) so the aquifer stays hydraulically continuous through the
          interior faces too.

        For a lake sealed on EVERY side (a concrete bottomless vault) use
        ``only_vertical`` instead -- that removes the perimeter connections this
        mode still keeps.
        """
        value = self.rectangular_interior
        if isinstance(value, Mapping):
            found, value = _lake_keyed(value, lake_id, name="rectangular_interior")
            return bool(value) if found else False
        return bool(value)

    def _top(self, lake_id: str, cell: int) -> float:
        """Facility top elevation for rectangular lakes (the connection ``telev``)."""

        top = self._per_cell("lake_top", lake_id)
        return top if isinstance(top, float) else top[cell]

    def _elevation_source(self) -> tuple[pd.DataFrame | None, dict]:
        """The elevation table (geometry dropped) and the excavation records beside it.

        Read once, through ``specs._layer_elevation_source`` -- the one resolver that
        returns the records from the SAME place as the elevations (``context.surfaces``
        when given, else the grid's ``gdf_topbtm`` and ``excavations``), so a lake can
        never be given the excavated rule from one object and elevations from another.
        """

        if self._source is None:
            table, records = _layer_elevation_source(self.context)
            if table is not None:
                table = table.drop(columns="geometry", errors="ignore")
            object.__setattr__(self, "_source", (table, records))
        return self._source

    def _surfaces(self) -> pd.DataFrame:
        """The per-cell top/bottom elevation table (context surfaces, else the grid's), geometry dropped."""

        surfaces = self._elevation_source()[0]
        if surfaces is None:
            raise ValueError("ModelContext.surfaces or grid.gdf_topbtm is required.")
        return surfaces

    def _elevation_array(self) -> np.ndarray:
        """``(ncpl, nlay + 1)`` floats: each cell's model top, then its layer bottoms."""

        if self._elevations is None:
            numbers = self._surfaces().select_dtypes("number")
            array = numbers.reindex(range(self.grid.ncpl)).to_numpy(dtype=float)
            object.__setattr__(self, "_elevations", array)
        return self._elevations

    def _top_active(self, cell: int) -> int | None:
        """The highest layer at ``cell`` that is active AND has thickness, or ``None``.

        What an excavated lake's VERTICAL connection lands in (the design's section
        4.2). Activity alone is not enough: a layer the cut left at zero thickness has
        ``top == bottom == bed`` and would still "contain" the bed, and with no
        ``domain`` on the context every layer reads active.
        """

        column = self._elevation_array()[int(cell)]
        for layer in range(len(column) - 1):
            if column[layer] - column[layer + 1] > 0.0 and self._active(layer, cell):
                return layer
        return None

    def _active(self, layer: int, cell: int) -> bool:
        """Whether ``(layer, cell)`` is active in the context's idomain (True if none set).

        The three-valued idomain rule itself is
        :func:`myflopy.modflow.mf6.grid.geometry.domain_active`. This used to be a
        second copy of it that differed twice: a ONE-dimensional domain counted as
        active everywhere, so a lake connected through cells the model had switched
        off; and a layer beyond a ``(nlay, ncpl)`` domain raised ``IndexError``
        rather than naming the layer.
        """

        return domain_active(self.context.domain, layer, cell)

    def _layer_containing(self, cell: int, elevation: float, only_layer: int | None = None) -> int | None:
        """The active layer whose top/bottom bracket ``elevation`` at ``cell`` (or ``None``).

        ``only_layer`` restricts the search to a single layer.
        """

        surfaces = self._surfaces().loc[cell].to_numpy(dtype=float)
        for layer, (top, bottom) in enumerate(zip(surfaces, surfaces[1:], strict=False)):
            if only_layer is not None and layer != only_layer:
                continue
            if top >= elevation >= bottom and self._active(layer, cell):
                return layer
        return None

    def _bottom(self, lake_id: str, cell: int) -> float:
        """The lakebed elevation at ``cell`` -- one of the lake's cells, or a neighbour.

        For an EXCAVATED lake it is the bed each connection is written at, the top of
        the cell's top active layer (:meth:`_excavated_bed`); only lake cells are asked.
        """

        if self._rule(lake_id) == "excavated":
            return self._excavated_bed(lake_id)[int(cell)]
        bed = self._per_cell("lake_bottom", lake_id)
        return bed if isinstance(bed, float) else bed[cell]

    def _per_cell(self, name: str, lake_id: str) -> float | dict[int, float]:
        """``lake_bottom`` or ``lake_top`` for one lake: one flat value, or one per cell.

        Resolved once per lake. A Surface or raster path is sampled onto the grid
        by :meth:`_sample`; a cell-indexed Series or ``{cell: value}`` mapping is
        read as given. Anything else is one value for the whole lake.
        """

        key = (name, lake_id)
        if key not in self._per_lake:
            if name == "lake_bottom" and self.lake_bottom is None:
                raise ValueError(
                    f"lake_bottom is required for lake {lake_id!r}: no excavation record "
                    "reached mf.lak for it, and only an EXCAVATED lake takes its bed from "
                    "one (LayerStack.build(attach=True), or ModelContext(surfaces=<the "
                    "LayerBuildResult>)). Give lake_bottom=, or excavate the pond."
                )
            value = self._lake_value(getattr(self, name), lake_id, name=name, per_cell=True)
            if isinstance(value, (Path, Surface, Mapping, pd.Series)):
                self._per_lake[key] = self._cell_values(name, lake_id, value)
            else:
                self._per_lake[key] = float(value)
        return self._per_lake[key]

    def _cell_values(self, name: str, lake_id: str, source: Any) -> dict[int, float]:
        """Per-cell ``lake_bottom``/``lake_top`` over the lake's cells (and a bed's neighbours).

        ``bathy`` compares each lake cell's bed with its NEIGHBOURS', so a bed is
        read at the cells around the lake too. Those are land, and where the
        source has no value for one -- a raster clipped to the lake, a Series of
        the lake's own cells -- the bed there is the model top. That used to be an
        accident: a Series had to list them or raise, and a raster's NaN reached
        ``min(cell_top, nan)``, which returns the cell top only because it is the
        first argument.

        A bed is capped at its cell's model top, because a lakebed is not above the
        ground it is cut into. A lake cell capped that way WARNS, since more than
        rounding means the two surfaces disagree. A lake cell a raster has no data
        for (usually one the polygon clips at its edge) is set to the top, zero
        depth, and WARNS; one missing from a Series RAISES, because that is an
        omission rather than a sampling artefact. ``lake_top`` is read at the
        lake's cells only, and is neither filled nor capped.
        """

        bed = name == "lake_bottom"
        cells = self.lake_cells[lake_id]
        wanted = [*cells, *self._neighbours(lake_id)] if bed else list(cells)
        sampled = isinstance(source, (Path, Surface))
        if sampled:
            grid_values = self._sample(source)
            values = {cell: float(grid_values[cell]) for cell in wanted}
        else:
            values = {cell: float(source[cell]) if cell in source else np.nan for cell in wanted}
        missing = [cell for cell in cells if not np.isfinite(values[cell])]
        if missing and not (sampled and bed):
            if sampled:
                reason = "the raster has no data there"
            else:
                reason = "a per-cell value must cover every cell of the lake"
                if bed:
                    reason += " (only its neighbours default, to the model top)"
                if isinstance(source, pd.Series):
                    reason += "; a Series is read by CELL unless its index names every lake id"
            raise ValueError(f"{name} for '{lake_id}' has no value at lake cell(s) {missing}: {reason}.")
        if not bed:
            return values
        if len(missing) == len(cells):
            raise ValueError(
                f"lake_bottom has no data at any of the {len(cells)} cells of lake "
                f"'{lake_id}' -- does the raster cover the lake?"
            )

        model_top = self._surfaces().select_dtypes("number").iloc[:, 0]
        above = [cell for cell in cells if values[cell] > model_top.loc[cell]]
        excess = max((values[cell] - float(model_top.loc[cell]) for cell in above), default=0.0)
        for cell, value in values.items():
            top = float(model_top.loc[cell])
            if not np.isfinite(value) or value > top:
                values[cell] = top
        if missing:
            warnings.warn(
                f"lake_bottom has no data at {len(missing)} of the {len(cells)} cells of "
                f"lake '{lake_id}' {missing[:10]}; their bed is set to the model top "
                "(zero depth). Usually cells the lake polygon only clips at its edge.",
                UserWarning,
                stacklevel=2,
            )
        if above:
            warnings.warn(
                f"lake_bottom is above the model top at {len(above)} of the {len(cells)} "
                f"cells of lake '{lake_id}', by up to {excess:.3g}; capped at the top "
                "there. More than rounding means the bed and the model top disagree.",
                UserWarning,
                stacklevel=2,
            )
        return values

    def _neighbours(self, lake_id: str) -> list[int]:
        """The cells sharing a face with the lake that are not part of it."""

        cells = self.lake_cells[lake_id]
        around = {int(other) for cell in cells for other in self.grid.find_adjacent_cells(cell)}
        return sorted(around - set(cells))

    def _sample(self, source: Path | Surface) -> np.ndarray:
        """Sample a raster or Surface onto every cell of the grid, once per source.

        Each lake indexes the result. It used to be sampled once for every
        distinct cell read -- the whole grid, area-weighted, each time: 5 samples
        for a 2-cell pond with 3 neighbours, so the cost grew as (lake + shore
        cells) x grid size. A Surface defined relative to the one above it
        (``offset_below``, ``fill="propagate"``) takes the model top as that
        surface, and ``units=`` converts to the context's length units, as it
        does in a ``LayerStack``.
        """

        key = str(source) if isinstance(source, Path) else id(source)
        if key not in self._sampled:
            model_top = self._surfaces().select_dtypes("number").iloc[:, 0]
            previous = model_top.reindex(range(self.grid.ncpl)).to_numpy(dtype=float)
            self._sampled[key] = sample_bed(
                source,
                self.grid,
                previous,
                length_units=getattr(self.context, "length_units", None),
            )
        return self._sampled[key]

    def _leakance(self, lake_id: str, connection_type: str, explicit: float | None = None) -> float:
        """The (positive) lakebed leakance for one lake + connection type (explicit wins).

        Resolves per-lake and per-connection-type mappings; raises if the resolved
        value is missing or not positive.
        """

        if explicit is not None:
            return float(explicit)
        pair = (
            [
                key
                for key in self.bed_leakance
                if isinstance(key, tuple)
                and len(key) == 2
                and str(_plain_scalar(key[0])) == lake_id
                and key[1] == connection_type.lower()
            ]
            if isinstance(self.bed_leakance, Mapping)
            else []
        )
        if pair:
            value = self.bed_leakance[pair[0]]
        else:
            value = self._lake_value(self.bed_leakance, lake_id, name="bed_leakance")
        if isinstance(value, Mapping):
            value = value.get(connection_type.lower(), value.get(connection_type.upper()))
        if value is None or float(value) <= 0.0:
            raise ValueError(f"bed_leakance must be positive for lake '{lake_id}'.")
        return float(value)

    def _adjacent_with_metrics(self, cell: int) -> list[tuple[int, float, float]]:
        """Neighbors of ``cell`` as ``(other, connection_length, connection_width)`` triples.

        Reads the grid's CRS/connection arrays (``ia``/``cl12``/``hwva``) via O(1)
        row pointers so building a large lake stays linear. Row ``cell``, so the length
        is ``cell``'s OWN side: its DISV centre to the line of the shared face -- MF6's
        geometry, the rule :class:`~myflopy.modflow.mf6.grid.connectivity.DisvFaces`
        states once (ledger 328; until then ``cl12`` held the full centre-to-centre
        spacing and a legacy wall's ``CONNLEN`` was about twice MF6's).
        """

        adjacent = list(self.grid.find_adjacent_cells(cell))
        start = int(self.grid.ia[cell])  # O(1) row pointer; sum(iac[:cell]) was O(N^2) over a lake
        return [
            (int(other), float(self.grid.cl12[start + index + 1]), float(self.grid.hwva[start + index + 1]))
            for index, other in enumerate(adjacent)
        ]

    def _vertical_connection(self, lake_id: str, cell: int, mode: str) -> LakeConnection:
        """One vertical (lake-bottom) connection for a cell.

        ``bathy``: at the layer containing the cell's lake bottom.
        ``rectangular``: the flat-bottom infiltration face; the connection spans
        ``lake_bottom -> lake_top`` and lands in the layer holding the lake bottom.
        MF6 ignores ``belev``/``telev`` for vertical connections, but they are set to
        the physical extent for readability and parity with the legacy output. A bed
        in no active layer (of ``only_layer``, when given) raises in BOTH modes.
        """
        bottom = self._bottom(lake_id, cell)
        only_layer = self._only_layer(lake_id)
        layer = self._layer_containing(cell, bottom, only_layer)
        if layer is None:
            # Both modes refuse: a "rectangular" lake used to drop this cell's VERTICAL
            # connection silently -- 8 HORIZONTAL and 0 VERTICAL written, nothing raised
            # (ledger 318) -- where "bathy" raised in the same place.
            where = "any active layer" if only_layer is None else f"only_layer={only_layer} (or it is inactive there)"
            raise ValueError(
                f"Lake '{lake_id}' bottom {bottom} does not intersect active cell {cell}: "
                f"the bed is not inside {where}, so the cell's VERTICAL connection has no "
                "cell to land in. Move lake_bottom or only_layer, or give the lake its "
                "connections by hand (connection_modes={lake_id: [LakeConnection(...)]})."
            )
        top = self._top(lake_id, cell) if mode == "rectangular" else bottom
        return LakeConnection(
            cellid=(layer, cell),
            connection_type="VERTICAL",
            bed_leakance=self._leakance(lake_id, "vertical"),
            bottom_elevation=bottom,
            top_elevation=top,
        )

    def _sidewall_connections(self, lake_id: str, cell: int, mode: str) -> list[LakeConnection]:
        """Horizontal (sidewall) connections for one lake cell.

        For each adjacent cell we work out the vertical extent of face exposed to the
        lake, then clip it to every model layer it spans -- so a single neighbor can
        yield connections in more than one layer. The exposed face is:

        * ``rectangular`` (trench/vault): edges only -- skip neighbors that are also
          lake cells (a flat bottom has no interior face); the face runs from the flat
          ``lake_bottom`` up to the facility ``lake_top``.
        * ``bathy`` (natural lake): connect toward a neighbor only where this cell's
          lake bottom is *below* the neighbor's (an exposed step); the face runs from
          this cell's lake bottom up to the neighbor's lake bottom.

        ``belev``/``telev`` are the face clipped to each layer's [bottom, top]; the
        connection geometry never depends on the (transient) lake stage -- MF6 wets
        and dries the face itself each timestep.
        """
        lake_cells = set(self.lake_cells[lake_id])
        lake_bottom = self._bottom(lake_id, cell)
        surfaces = self._surfaces().loc[cell].to_numpy(dtype=float)
        only_layer = self._only_layer(lake_id)
        interior = self._rectangular_interior(lake_id)
        result = []
        for adjacent, length, width in self._adjacent_with_metrics(cell):
            if mode == "rectangular":
                if adjacent in lake_cells and not interior:
                    continue  # lined basin: no interior face (permeable fill -> interior=True)
                face_bottom = lake_bottom
                face_top = self._top(lake_id, cell)
            else:  # bathy: exposed step only where this cell is deeper than the neighbor
                neighbor_bottom = self._bottom(lake_id, adjacent)
                if lake_bottom >= neighbor_bottom:
                    continue
                face_bottom = lake_bottom
                face_top = neighbor_bottom
            for layer, (cell_top, cell_bottom) in enumerate(zip(surfaces, surfaces[1:], strict=False)):
                if only_layer is not None and layer != only_layer:
                    continue
                if not self._active(layer, cell):
                    continue
                belev = max(cell_bottom, face_bottom)
                telev = min(cell_top, face_top)
                # A face thinner than FACE_TOL is dropped: MF6 reads belev == telev as
                # the WHOLE cell, and 8-decimal output makes a 1e-9 face exactly that
                # (the design's F9). The rows are otherwise the rule before ledger 325's,
                # with CONNLEN the lake cell's own side of the face since ledger 328.
                if telev - belev < FACE_TOL:
                    continue
                result.append(
                    LakeConnection(
                        cellid=(layer, cell),
                        connection_type="HORIZONTAL",
                        bed_leakance=self._leakance(lake_id, "horizontal"),
                        bottom_elevation=float(belev),
                        top_elevation=float(telev),
                        connection_length=length,
                        connection_width=width,
                    )
                )
        return result

    # -- the excavated rule (ledger 325; the excavation design's section 4) ----------

    def _rule(self, lake_id: str) -> str:
        """``"excavated"``, ``"legacy"`` or ``"explicit"`` for one lake (see :attr:`lake_rules`)."""

        if self._rules is None:
            self._resolve_rules()
        return self._rules[lake_id]

    @property
    def lake_rules(self) -> dict[str, str]:
        """Which rule each lake's connections are written by, keyed by lake id.

        Chosen from the EXCAVATION RECORD, never from the geometry (the excavation
        design's section 4.1: a geometric test was measured wrong three ways -- it
        changed output for lakes nobody excavated, missed a stack cut deeper than the
        lake's bed, and refused a flat-floored basin that runs). The records come
        through the same resolver as the layer elevations, in its order:
        ``context.surfaces`` when given (a ``LayerBuildResult`` carries its
        ``excavations``; a bare table carries none unless it IS the grid's own
        ``gdf_topbtm``), else the grid's ``excavations``, which ``LayerStack.build(
        attach=True)`` publishes beside ``gdf_topbtm``.

        * ``"excavated"`` -- the lake's cells are exactly one excavation's cells. Its
          VERTICAL connections go to each cell's top active cell, whose top is the bed,
          and its HORIZONTAL ones into the NEIGHBOURING columns (``connection_modes``).
        * ``"legacy"`` -- no excavation touches the lake: the rule ``mf.lak`` always
          had, plus the ``unexcavated=`` guard.
        * ``"explicit"`` -- the lake's connections were given as ``LakeConnection``
          objects and are written as given; no excavation is consulted.

        Returns
        -------
        dict of str to str
            A copy, in lake order. Reading it resolves the lake cells and the records,
            and raises what that resolution refuses.

        Raises
        ------
        ValueError
            A lake covering cells of two excavations; or the lakes touching an
            excavation not covering exactly its cells (a lake cell the excavation did
            not dig, or a dug cell no lake of this package covers) -- whatever
            ``unexcavated=`` says.

        Examples
        --------
        >>> builder.lake_rules
        {'pond': 'excavated'}
        """

        if self._rules is None:
            self._resolve_rules()
        return dict(self._rules)

    def _resolve_rules(self) -> None:
        """Match every generated lake against the excavation records, refusing a partial match."""

        generated = [lake_id for lake_id in self.lake_ids if isinstance(self._connection_mode(lake_id), str)]
        records = self._elevation_source()[1] if generated else {}
        rules: dict[str, str] = {}
        matched: dict[str, Any] = {}
        touching: dict[str, list[str]] = {}
        for lake_id in self.lake_ids:
            if lake_id not in generated:
                rules[lake_id] = "explicit"
                continue
            cells = set(self.lake_cells[lake_id])
            hits = [record for record in records.values() if cells & set(record.cells)]
            if len(hits) > 1:
                names = ", ".join(repr(record.name) for record in hits)
                raise ValueError(
                    f"lake {lake_id!r} covers cells of {len(hits)} excavations ({names}). An "
                    "excavated lake is written at the bed of the ONE excavation cut for it; "
                    "a lake spanning several would take two beds' rules at once. Give each "
                    "pond its own lake (one feature per pond in lakes=), matching the "
                    "excavation cut for it."
                )
            if hits:
                rules[lake_id] = "excavated"
                matched[lake_id] = hits[0]
                touching.setdefault(hits[0].name, []).append(lake_id)
            else:
                rules[lake_id] = "legacy"
        for name, lakes in touching.items():
            record = matched[lakes[0]]
            covered = set().union(*(set(self.lake_cells[lake_id]) for lake_id in lakes))
            dug = set(record.cells)
            if covered == dug:
                continue
            undug, spare = sorted(covered - dug), sorted(dug - covered)
            parts = []
            if undug:
                parts.append(f"{len(undug)} lake cell(s) the excavation did not dig ({_listing(undug)})")
            if spare:
                parts.append(f"{len(spare)} dug cell(s) no lake of this package covers ({_listing(spare)})")
            label = ", ".join(repr(lake_id) for lake_id in lakes)
            raise ValueError(
                f"lake(s) {label} and excavation {name!r} cover different cells: "
                + " and ".join(parts)
                + ". mf.lak writes a lake by the excavated rule only when the lake IS the "
                "excavation -- the same polygons, e.g. LayerStack.excavate(within=<the "
                "lakes= polygons>, layer=<the same layer>). Part of a lake over a pond "
                "would exchange at the bed and the rest from the ground, which is neither "
                "rule. This is refused whatever unexcavated= says."
            )
        object.__setattr__(self, "_rules", rules)
        object.__setattr__(self, "_records_of", matched)

    def _check_declared_bed(self, lake_id: str, record: Any) -> None:
        """Refuse a ``lake_bottom`` that is not the bed ``record`` was cut to.

        Compared BEFORE any capping, cell by cell, within
        :func:`~myflopy.grid_identity.bed_tolerance` of the record's bed (the same
        Surface samples to the identical float; no data on both sides is equal). A
        ``lake_bottom`` left out -- ``None``, or a mapping without this lake -- IS the
        record's bed (the design's F13), so there is nothing to compare.
        """

        value = self.lake_bottom
        if value is None or (
            isinstance(value, Mapping) and not _lake_keyed(value, lake_id, name="lake_bottom")[0]
        ):
            return
        declared = self._lake_value(value, lake_id, name="lake_bottom", per_cell=True)
        cells = self.lake_cells[lake_id]
        if isinstance(declared, (Path, Surface)):
            sampled = self._sample(declared)
            raw = {cell: float(sampled[cell]) for cell in cells}
        elif isinstance(declared, (Mapping, pd.Series)):
            raw = {cell: float(declared[cell]) if cell in declared else np.nan for cell in cells}
        else:
            raw = dict.fromkeys(cells, float(declared))
        different = []
        for cell in cells:
            mine, theirs = raw[cell], record.bed_at(cell)
            if np.isnan(mine) and np.isnan(theirs):
                continue
            if np.isnan(mine) or np.isnan(theirs) or abs(mine - theirs) > bed_tolerance(theirs):
                different.append((cell, mine, theirs))
        if not different:
            return
        worst = max(
            different,
            key=lambda item: np.inf if np.isnan(item[1] - item[2]) else abs(item[1] - item[2]),
        )
        raise ValueError(
            f"lake_bottom for lake {lake_id!r} is not the bed excavation {record.name!r} was "
            f"cut to: {len(different)} of its {len(cells)} cells differ, the worst at cell "
            f"{worst[0]} (lake_bottom {worst[1]:.6g}, excavation {worst[2]:.6g}). The stack "
            "was cut to another bottom than this lake declares, so its connections would sit "
            "at a bed its own lake_bottom contradicts. Leave lake_bottom out -- an excavated "
            "lake takes its excavation's bed -- or excavate with the same surface (a surface "
            "defined relative to the one above it samples differently off the excavated top)."
        )

    def _excavated_bed(self, lake_id: str) -> dict[int, float]:
        """``{cell: bed}`` for an excavated lake: the top of each cell's top active layer.

        The excavated top, ``min(declared bed, ground)``: exactly the bed where the
        stack was cut, the ground at a rim cell the bed did not reach below. It is the
        float the DISV is written from, so a wall's ``belev`` equals a cell surface to
        the last digit (MF6 checks a face against its cell with no tolerance).
        """

        key = ("excavated_bed", lake_id)
        if key not in self._per_lake:
            record = self._records_of[lake_id]
            self._check_declared_bed(lake_id, record)
            beds, missing = {}, []
            elevations = self._elevation_array()
            for cell in self.lake_cells[lake_id]:
                layer = self._top_active(cell)
                if layer is None:
                    missing.append(cell)
                    continue
                beds[int(cell)] = float(elevations[int(cell), layer])
            if missing:
                raise ValueError(
                    f"excavated lake {lake_id!r}: {len(missing)} of its cells have no "
                    f"active cell with thickness in any layer ({_listing(missing)}), so "
                    "there is nothing under the bed to connect to. Check the context's "
                    "domain= is the excavated result's idomain."
                )
            self._per_lake[key] = beds
        return self._per_lake[key]

    def _faces(self) -> DisvFaces:
        """MF6's own DISV face geometry for the grid (built once per builder)."""

        if self._face_geometry is None:
            object.__setattr__(self, "_face_geometry", DisvFaces(self.grid.get_disv_gridprops()))
        return self._face_geometry

    def _excavated_connections(self, lake_id: str, mode: str) -> list[tuple[int, str, LakeConnection]]:
        """``(lake cell, role, connection)`` for an excavated lake, lake cell by lake cell.

        Per lake cell ``A`` at bed ``b``: one VERTICAL connection to the top active
        cell, whose top IS ``b``; then, toward each face neighbour ``N`` whose model top
        (its bed, for a lake cell; the ground, for land) is above ``b`` -- capped at
        ``lake_top`` for ``"rectangular"`` -- one HORIZONTAL connection per active layer
        of ``N`` the wall spans, written INTO ``N`` (the aquifer cell the water meets),
        ``belev = max(b, bottom)``, ``telev = min(wall top, top)``, ``CONNLEN`` the
        distance from ``N``'s DISV centre to the line of the shared face and
        ``CONNWIDTH`` the face's length -- MF6's own DISV geometry. Faces thinner than
        :data:`FACE_TOL` are dropped.
        """

        mode_label = "rectangular" if mode == "rectangular" else "bathy"
        if mode == "rectangular" and self._rectangular_interior(lake_id):
            raise ValueError(
                f"rectangular_interior=True on excavated lake {lake_id!r}: an excavated lake "
                "is open water over its whole floor, so there are no interior faces for a "
                "permeable fill to conduct through -- it would change nothing. A "
                "gravel-filled basin is aquifer, not lake: do not excavate it, and give it "
                "rectangular_interior=True on the un-excavated stack."
            )
        beds = self._excavated_bed(lake_id)
        record = self._records_of[lake_id]
        only_layer = self._only_layer(lake_id)
        only_vertical = self._only_vertical(lake_id)
        elevations = self._elevation_array()
        vertical_leakance = self._leakance(lake_id, "vertical")
        result: list[tuple[int, str, LakeConnection]] = []
        for cell in self.lake_cells[lake_id]:
            cell = int(cell)
            bed = beds[cell]
            layer = self._top_active(cell)
            if only_layer is not None and layer != only_layer:
                raise ValueError(
                    f"only_layer={only_layer} for excavated lake {lake_id!r}, but its bed at "
                    f"cell {cell} is the top of layer {layer}, the top active cell there. An "
                    "excavated lake's floor connects to the cell under the bed; drop "
                    "only_layer, or excavate so the bed sits on that layer."
                )
            declared = record.bed_at(cell)
            role = "rim" if (not np.isfinite(declared) or declared > bed + bed_tolerance(bed)) else "floor"
            top = self._top(lake_id, cell) if mode_label == "rectangular" else bed
            result.append(
                (cell, role, LakeConnection((layer, cell), "VERTICAL", vertical_leakance, bed, top))
            )
            if only_vertical:
                continue
            ceiling = self._top(lake_id, cell) if mode_label == "rectangular" else None
            for neighbour in self.grid.find_adjacent_cells(cell):
                neighbour = int(neighbour)
                wall_top = float(elevations[neighbour, 0])
                if ceiling is not None:
                    wall_top = min(ceiling, wall_top)
                if not wall_top > bed + FACE_TOL:
                    continue
                metrics = self._faces().length_and_width(cell, neighbour, toward=neighbour)
                if metrics is None:
                    continue  # they meet at a corner only: MF6 connects no face there
                length, width = metrics
                column = elevations[neighbour]
                for index in range(len(column) - 1):
                    if only_layer is not None and index != only_layer:
                        continue
                    if not self._active(index, neighbour):
                        continue
                    belev = max(bed, float(column[index + 1]))
                    telev = min(wall_top, float(column[index]))
                    if telev - belev < FACE_TOL:
                        continue
                    result.append(
                        (
                            cell,
                            "wall",
                            LakeConnection(
                                (index, neighbour),
                                "HORIZONTAL",
                                self._leakance(lake_id, "horizontal"),
                                belev,
                                telev,
                                length,
                                width,
                            ),
                        )
                    )
        return result

    def _unexcavated_report(self) -> str | None:
        """The ``unexcavated=`` guard's findings for every LEGACY lake, or ``None``.

        Two findings, both about a lake no excavation record reached:

        * VERTICAL connections in a cell whose TOP is above the lake's bed -- MF6
          counts a vertical connection from its cell's top, so the storage and the
          exchange between the bed and that top are missing (on a real pond, 577 of
          612 cells, worst 12.58 ft, median 3.35);
        * a lake that LOOKS excavated -- more than half its cells have a zero-thickness
          layer above an active cell whose top is the bed -- yet was written by the old
          rule, which puts no wall into the neighbours: the record never arrived.
        """

        lines = []
        rows = list(zip(self.connections, self._roles, strict=True))
        for lake_id in self.lake_ids:
            if self._rule(lake_id) != "legacy":
                continue
            cells = self.lake_cells[lake_id]
            elevations = self._elevation_array()
            gaps, tops = [], []
            for (owner, connection), _role in rows:
                if owner != lake_id or connection.connection_type != "VERTICAL":
                    continue
                layer, cell = connection.cellid
                cell_top = float(elevations[int(cell), int(layer)])
                bed = self._bottom(lake_id, cell)
                if cell_top - bed > bed_tolerance(bed):
                    gaps.append(cell_top - bed)
                    tops.append(cell_top)
            looks = 0
            for cell in cells:
                layer = self._top_active(cell)
                if not layer:
                    continue
                column = elevations[int(cell)]
                if not (column[:layer] - column[1 : layer + 1] <= 0.0).any():
                    continue
                bed = self._bottom(lake_id, cell)
                if abs(float(column[layer]) - bed) <= bed_tolerance(bed):
                    looks += 1
            parts = []
            if gaps:
                total = sum(
                    1
                    for (owner, connection), _role in rows
                    if owner == lake_id and connection.connection_type == "VERTICAL"
                )
                stage = float(self._lake_value(self.starting_stage, lake_id, name="starting_stage"))
                above = sum(1 for top in tops if top > stage)
                text = (
                    f"lake {lake_id!r}: {len(gaps)} of its {total} VERTICAL connections sit in a "
                    f"cell whose TOP is above the bed, by up to {max(gaps):.4g} (median "
                    f"{float(np.median(gaps)):.4g}); {above} of those tops are above its "
                    f"starting stage {stage:g}"
                )
                table = self.prepared_tables.get(lake_id)
                if table is not None:
                    over = sum(1 for top in tops if top > table.maximum_stage)
                    text += f", {over} above its table's highest stage {table.maximum_stage:g}"
                parts.append(
                    text + ". MF6 counts a VERTICAL connection from its cell's top, so each "
                    "exchanges nothing until the stage reaches that top, and the storage "
                    "between the bed and the top is missing without a table"
                )
            if looks > len(cells) / 2:
                parts.append(
                    f"lake {lake_id!r} LOOKS excavated -- {looks} of its {len(cells)} cells have "
                    "a zero-thickness layer above an active cell whose top is the bed -- but no "
                    f"excavation record reached mf.lak, so it was written by the "
                    f"{self._connection_mode(lake_id)!r} rule, with no walls into its "
                    "neighbours"
                )
            lines.extend(parts)
        if not lines:
            return None
        return (
            "mf.lak: lake(s) NOT excavated into the layer stack.\n  - "
            + "\n  - ".join(lines)
            + "\nCut the pond into the stack -- stack.excavate(bottom=<the same "
            "lake_bottom>, within=<the lakes= polygons>) -- rebuild it with attach=True (or "
            "pass the LayerBuildResult as ModelContext(surfaces=...)), and re-declare "
            "mf.lak: those two spellings are what carry an excavation's record to it. With "
            "no stack (a DISV built from a loaded run's arrays), cut the LayerBuildResult "
            "the DISV was built from instead: dug = result.excavate(bottom=<the same "
            "lake_bottom>, within=<the lakes= polygons>) returns a NEW result -- build the "
            "DISV from it (mf.disv(dug)) and pass ModelContext(surfaces=dug, "
            "domain=dug.idomain). unexcavated='ignore' silences this; 'raise' refuses it."
        )

    def _apply_unexcavated_guard(self) -> None:
        """Warn about, refuse, or ignore :meth:`_unexcavated_report`, per ``unexcavated=``."""

        if self.unexcavated == "ignore":
            return
        report = self._unexcavated_report()
        if report is None:
            return
        if self.unexcavated == "raise":
            raise ValueError(report)
        logger.warning(report)
        warnings.warn(report, UserWarning, stacklevel=stacklevel_outside_package())

    @property
    def connections(self) -> tuple[tuple[str, LakeConnection], ...]:
        """Return generated or explicit connections keyed by stable lake ID."""

        if self._connections is not None:
            return self._connections
        result = []
        roles: list[tuple[str, int | None, str]] = []
        for lake_id in self.lake_ids:
            mode = self._connection_mode(lake_id)
            if not isinstance(mode, str):
                for connection in mode:
                    result.append((lake_id, connection))
                    roles.append(
                        (lake_id, None, "floor" if connection.connection_type == "VERTICAL" else "wall")
                    )
                continue
            if self._rule(lake_id) == "excavated":
                for cell, role, connection in self._excavated_connections(lake_id, mode):
                    result.append((lake_id, connection))
                    roles.append((lake_id, cell, role))
                continue
            only_vertical = self._only_vertical(lake_id)
            lake_cells = self.lake_cells[lake_id]
            horizontal = 0
            for cell in lake_cells:
                vertical = self._vertical_connection(lake_id, cell, mode)
                if vertical is not None:
                    result.append((lake_id, vertical))
                    roles.append((lake_id, int(cell), "floor"))
                if only_vertical:  # vault: concrete sidewalls, no horizontal exchange
                    continue
                for connection in self._sidewall_connections(lake_id, cell, mode):
                    horizontal += 1
                    result.append((lake_id, connection))
                    roles.append((lake_id, int(cell), "wall"))
            # Guard the degenerate bathy+flat-bottom config: 'bathy' connects a
            # cell to a neighbor only where its lake bottom is BELOW the
            # neighbor's (an exposed step). A flat lake_bottom has no steps, so a
            # multi-cell bathy lake produces ZERO horizontal connections -- a
            # lake with no lateral aquifer exchange, which is almost never
            # intended (an earlier canonical config dropped from 39 horizontal
            # connections to 0 this way). The geometry is internally consistent,
            # so warn (rather than raise) and point at the three real choices.
            if (
                mode == "bathy"
                and not only_vertical
                and len(lake_cells) > 1
                and horizontal == 0
            ):
                warnings.warn(
                    f"Lake '{lake_id}' uses 'bathy' connections but produced no "
                    f"horizontal (sidewall) connections across its "
                    f"{len(lake_cells)} cells: its lake_bottom is flat, so there "
                    f"are no exposed steps to connect through, leaving the lake "
                    f"with no lateral aquifer exchange. If that is not intended, "
                    f"either give lake_bottom real bathymetry (a Surface such "
                    f"as mf.Raster or mf.Contours, a raster path, or a "
                    f"cell-indexed Series), or use "
                    f"connection_modes='rectangular' for a flat-bottomed "
                    f"infiltration basin (edges only), or only_vertical=True for "
                    f"a sealed concrete vault (bottom only).",
                    UserWarning,
                    stacklevel=2,
                )
        object.__setattr__(self, "_connections", tuple(result))
        object.__setattr__(self, "_roles", tuple(roles))
        return self._connections

    @property
    def connection_roles(self) -> dict[str, list[list[Any]]]:
        """What each written connection IS, per lake: ``[[iconn, lake_cell, role], ...]``.

        ``iconn`` is MF6's connection number within the lake (from zero), ``role`` is
        ``"floor"`` (a VERTICAL connection under a lake cell), ``"rim"`` (a VERTICAL one
        at a cell of an excavated lake whose declared bed did not reach below the ground,
        so it connects at the ground) or ``"wall"`` (a HORIZONTAL one), and
        ``lake_cell`` the lake cell the connection belongs to -- which a wall's own row
        no longer says for an excavated lake, whose walls are written into the
        NEIGHBOUR's cell. ``None`` for a lake whose connections were given by hand.
        Stamped on the spec's metadata as ``lake_connections``.

        Returns
        -------
        dict of str to list of list
            ``{lake_id: rows}``, in connection order, every value a plain python type.

        Examples
        --------
        >>> builder.connection_roles["pond"][:2]
        [[0, 5, 'floor'], [1, 5, 'wall']]
        """

        self.connections
        counts: dict[str, int] = {}
        result: dict[str, list[list[Any]]] = {lake_id: [] for lake_id in self.lake_ids}
        for lake_id, cell, role in self._roles:
            index = counts.get(lake_id, 0)
            counts[lake_id] = index + 1
            result[lake_id].append([index, None if cell is None else int(cell), role])
        return result

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
        """Expand a stage/rainfall/etc setting into per-period MF6 ``[lake, keyword, amount]`` rows.

        Accepts a scalar (applied to all lakes/periods), a per-period mapping, a
        per-lake mapping (scalar or per-period series), or a per-period sequence; a
        ``numpy`` array or ``pandas`` Series is first turned into the list or mapping it
        means (:func:`_per_period_form`). Refuses what would never be written or what
        MF6 refuses one run later (ledger 318): a period outside ``range(nper)``, a list
        longer than ``nper``, an all-integer mapping whose every key is also a lake id
        (periods or lakes would be a guess), an unknown status, a negative forcing, and
        a number that is not finite.
        """

        name = _SETTING_NAMES[keyword]
        result = {period: [] for period in range(self.nper)}
        if value is None:
            return result
        value = _per_period_form(value, name=name)

        if isinstance(value, Mapping) and value and all(_is_period_key(key) for key in value):
            self._check_period_keys(value, name, top_level=True)
            for period, current in value.items():
                if current is None:
                    continue
                if isinstance(current, Mapping):
                    for lake_id, amount in current.items():
                        number = self._setting_lake(lake_id, name)
                        result[int(period)].append(
                            [number, keyword, self._setting_amount(amount, keyword, name)]
                        )
                    continue
                for lake_id in self.lake_ids:
                    result[int(period)].append(
                        [self.lake_number(lake_id), keyword, self._setting_amount(current, keyword, name)]
                    )
            return result

        if isinstance(value, Mapping):
            for lake_id, series in value.items():
                number = self._setting_lake(lake_id, name)
                for period, amount in self._by_period(series, f"{name}[{lake_id!r}]"):
                    result[period].append([number, keyword, self._setting_amount(amount, keyword, name)])
            return result

        for period, amount in self._by_period(value, name):
            for lake_id in self.lake_ids:
                result[period].append(
                    [self.lake_number(lake_id), keyword, self._setting_amount(amount, keyword, name)]
                )
        return result

    def _setting_lake(self, lake_id: Any, name: str) -> int:
        """The lake number a per-period mapping's lake key names, or a refusal naming the value."""

        label = str(_plain_scalar(lake_id))
        if label not in self.lake_numbers:
            raise ValueError(
                f"{name} names lake {label!r}, which no lake has (lakes: "
                f"{', '.join(repr(lake) for lake in self.lake_ids)}). A mapping whose keys "
                "are ALL integers is read by stress period; any other is read by lake id."
            )
        return self.lake_number(label)

    def _check_period_keys(self, value: Mapping, name: str, *, top_level: bool = False) -> None:
        """Refuse a period key outside ``range(nper)``, and a mapping every key of which is also a lake id.

        A key past the last period is never written -- ``rainfall={101: 0.5}`` over
        integer lake ids wrote NOTHING, silently (ledger 318). And when every key of a
        TOP-LEVEL mapping is a period AND a lake id, which one was meant is a guess; a
        mapping inside a lake's entry can only be periods.
        """

        outside = [int(key) for key in value if not 0 <= int(key) < self.nper]
        lakes = (
            [str(int(key)) for key in value if str(int(key)) in self.lake_numbers] if top_level else []
        )
        hint = (
            f" Lake ids are strings, and every key here is an integer, so each is read as a "
            f"stress period; key a lake by its string id ({lakes[0]!r})."
            if lakes
            else ""
        )
        if outside:
            raise ValueError(
                f"{name}: {outside[0]} is not a stress period -- periods run from 0 to "
                f"{self.nper - 1} (zero-based), so a value there would never be written.{hint}"
            )
        scalar_values = not any(isinstance(item, Mapping) for item in value.values())
        if top_level and scalar_values and len(lakes) == len(value):
            raise ValueError(
                f"{name}: every key ({', '.join(str(key) for key in value)}) is both a stress "
                "period and a lake id, so which one was meant would be a guess. Key lakes by "
                "their string id, and give periods as a list (None skips one) or as "
                "{period: {lake_id: value}}."
            )

    def _by_period(self, series: Any, name: str) -> list[tuple[int, Any]]:
        """``(period, amount)`` for every period a per-period value reaches, ``None`` skipped.

        A scalar reaches every period; a list or tuple one per position; an all-integer
        mapping the periods it names. A list longer than ``nper`` and a period outside
        ``range(nper)`` are refused -- they would never be written (ledger 318).
        """

        if isinstance(series, Mapping):
            if not series:
                return []
            if not all(_is_period_key(key) for key in series):
                raise ValueError(
                    f"{name} is a mapping keyed by {next(iter(series))!r}; a per-period "
                    "mapping here is keyed by zero-based stress period (an integer)."
                )
            self._check_period_keys(series, name)
            return [(int(key), item) for key, item in series.items() if item is not None]
        if isinstance(series, (list, tuple)):
            if len(series) > self.nper:
                raise ValueError(
                    f"{name} has {len(series)} values for {self.nper} stress period(s): position "
                    "0 is the first period, so the last "
                    f"{len(series) - self.nper} would never be written."
                )
            return [(period, item) for period, item in enumerate(series) if item is not None]
        return [(period, series) for period in range(self.nper)]

    def _setting_amount(self, amount: Any, keyword: str, name: str) -> Any:
        """One period value, checked as MF6 will check it, as a plain python value.

        ``STATUS`` must be one of MF6's three words (any case, written as given); a
        numeric setting is a finite number -- non-negative for the five MF6 refuses
        below zero -- or a string, which MF6 reads as a time-series name. MF6 refused
        each of these by name, one run late (ledger 318).
        """

        amount = _plain_scalar(amount)
        if keyword == "STATUS":
            if not isinstance(amount, str) or amount.strip().upper() not in _LAKE_STATUSES:
                raise ValueError(
                    f"status {amount!r} is not a lake status: MF6 takes "
                    f"{', '.join(_LAKE_STATUSES)} (any case), and refuses anything else by name."
                )
            return amount
        if isinstance(amount, str):
            return amount
        if isinstance(amount, (bool, np.bool_)) or not isinstance(amount, Real):
            raise ValueError(
                f"{name} value {amount!r} is not a number (or a time-series name): one number "
                "per lake per period."
            )
        if not np.isfinite(amount):
            raise ValueError(f"{name} value {amount!r} is not finite.")
        if keyword in _NONNEGATIVE_SETTINGS and amount < 0.0:
            raise ValueError(
                f"{name} is {amount!r}: MF6 refuses a negative {keyword} by name. "
                + (
                    "It is a positive number that REMOVES water."
                    if keyword == "WITHDRAWAL"
                    else "Take water out with withdrawals= or an outlet."
                )
            )
        return amount

    @property
    def outletdata(self) -> list[list[Any]]:
        """MF6 LAK ``outlets`` records, resolving source/receiver lake names to numbers."""

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
            name = f"outlets[{outlet_number}].rate"
            for period, amount in self._by_period(rate, name):
                result[period].append([outlet_number, "RATE", self._setting_amount(amount, "RATE", name)])
        return result

    @property
    def prepared_tables(self) -> dict[str, LakeTable]:
        """Per-lake built :class:`LakeTable` objects (builders resolved; unknown lakes raise).

        Built ONCE per set of ``tables`` entries, and returned as a new dict each
        time. :meth:`validate` reads this inside its per-lake loop and :meth:`build`
        and the unexcavated check read it again, and each read rebuilt every
        ``LakeTableBuilder`` -- a DEM table re-read its raster and footprint every
        time, so one ``mf.lak`` call over two DEM-table lakes built 12 tables and
        issued each table's no-data warning 6 times, growing with the square of the
        lake count. The cache is keyed on the entries themselves (compared by
        identity, and held, so an id cannot be reused), so replacing a table in
        ``tables`` rebuilds.
        """

        entries = tuple(self.tables.items())
        cached = self._prepared
        if (
            cached is not None
            and len(cached[0]) == len(entries)
            and all(
                lake_id == seen_id and table is seen
                for (lake_id, table), (seen_id, seen) in zip(entries, cached[0], strict=True)
            )
        ):
            return dict(cached[1])
        result = {}
        for lake_id, table in entries:
            if lake_id not in self.lake_numbers:
                raise ValueError(
                    f"Lake table references unknown lake '{lake_id}' (lakes: "
                    f"{', '.join(repr(lake) for lake in self.lake_ids)})."
                )
            if not isinstance(table, (LakeTable, LakeTableBuilder)):
                raise TypeError(
                    f"tables[{lake_id!r}] is a {type(table).__name__}; give a LakeTable or a "
                    "LakeTableBuilder."
                )
            result[lake_id] = table.build() if isinstance(table, LakeTableBuilder) else table
        object.__setattr__(self, "_prepared", (entries, result))
        return dict(result)

    def _refuse_lakes_off_the_grid(self) -> None:
        """Raise when a lake covers no grid cell, saying why and where (ledger 211).

        A lake with no cells carries no LAK connection, which MF6 refuses, so ANY
        lake missing raises -- there is no partial-miss tier here as there is for a
        ``.gpkg`` source, where one layer can deliberately serve several grids. The
        message carries why each lake missed
        (:func:`~myflopy.modflow.mf6.grid.geometry.miss_reason`) and both extents
        with both CRSs
        (:func:`~myflopy.modflow.mf6.grid.geometry.where_they_are`), because a
        location mismatch is two boxes. It used to be ``Lakes do not intersect the
        grid: <names>``, which named the lake and nothing that would find it.

        Only TWO of ``miss_reason``'s answers can come back here -- "outside the
        grid" and "only borders cells" -- because a lake is areal and nothing
        filters its cells: if any candidate overlaps the lake by area,
        ``cell_shares`` keeps it (a piece that is a sliver is dropped only while
        another piece is not), so an empty answer means no overlap at all. Its
        ``target`` names what a RESTRICTED search was looking for, and this one is
        unrestricted, so the default stands and its wording is never printed.
        """

        table = self.lake_table
        refuse_footprint_misses(
            self.grid.gdf_vorPolys,
            table,
            [self.lake_cells[lake_id] for lake_id in table.index],
            noun="lakes",
            consequence=", so they would carry no LAK connection",
            partial="raise",
            source_crs=self._source_crs,
            grid_crs=getattr(self.grid, "crs", None),
        )

    def validate(self) -> None:
        """Validate generated LAK inputs before building the package spec.

        Runs every check :meth:`build` runs, without returning the spec: the
        lakes land on the grid, the outlets are ``LakeOutlet`` objects naming
        known lakes, each generated lake's ``starting_stage`` is at or above its
        lowest bed and inside its table, a rectangular lake's top is above its
        bottom, and the connections and tables resolve.

        Raises
        ------
        ValueError
            The first refusal met (the list is under the class's Raises).

        See Also
        --------
        build : Validates, then returns the spec.

        Examples
        --------
        >>> builder.validate()
        """

        self._refuse_lakes_off_the_grid()
        if any(not isinstance(outlet, LakeOutlet) for outlet in self.outlets):
            raise ValueError("outlets must contain LakeOutlet objects.")
        for outlet in self.outlets:
            self.lake_number(outlet.source)
            if outlet.receiver is not None:
                self.lake_number(outlet.receiver)
        for lake_id in self.lake_ids:
            stage = float(self._lake_value(self.starting_stage, lake_id, name="starting_stage"))
            mode = self._connection_mode(lake_id)
            if isinstance(mode, str):
                cells = self.lake_cells[lake_id]
                bottoms = [self._bottom(lake_id, cell) for cell in cells]
                if stage < min(bottoms):
                    raise ValueError(f"starting_stage is below lake_bottom for lake '{lake_id}'.")
                if mode == "rectangular":
                    # Rectangular facilities need an explicit flat top; the connection
                    # telev is lake_top, never the (transient) starting stage.
                    tops = [self._top(lake_id, cell) for cell in cells]
                    if any(top <= bottom for top, bottom in zip(tops, bottoms, strict=False)):
                        raise ValueError(
                            f"lake_top must be above lake_bottom for rectangular lake '{lake_id}'."
                        )
            table = self.prepared_tables.get(lake_id)
            if table is not None and not table.minimum_stage <= stage <= table.maximum_stage:
                raise ValueError(f"starting_stage is outside the lake table range for lake '{lake_id}'.")
        self.connectiondata
        self.prepared_tables
        if self.implicit:
            connected = {lake_id for lake_id, _connection in self.connections}
            bare = [lake_id for lake_id in self.lake_ids if lake_id not in connected]
            if bare:
                raise ValueError(
                    f"implicit=True, and lake {bare[0]!r} has no connection"
                    + (f" (nor do {', '.join(repr(i) for i in bare[1:])})" if bare[1:] else "")
                    + ": MF6's IMPLICIT formulation solves each stage through its lakebed "
                    "seepage and stops on a lake with none ('lake N has no connections; the "
                    "IMPLICIT option requires each lake to have at least one GWF "
                    "connection', measured). Give it a connection, or leave implicit off."
                )
        self._apply_unexcavated_guard()

    def build(self) -> PackageSpec:
        """Validate stored configuration and return its package spec.

        Runs :meth:`validate`, then writes the package data, connections,
        period settings, outlets and tables through ``lak_spec``.

        Returns
        -------
        PackageSpec
            The LAK package spec, with one child ``laktab`` per table. Its
            metadata carries the lake ids, each lake's connection mode and the
            grid pin, but NOT the ``lak_index`` ``mf.lak`` adds, so hand it to
            ``mf.mvr`` through :meth:`connection`, not ``mf.lak_connection``.

        Raises
        ------
        ValueError
            Any refusal :meth:`validate` makes.
        TypeError
            An ``options`` name ``ModflowGwflak`` does not accept, or one the
            builder sets itself (``noutlets``, ``ntables``, the data blocks).

        See Also
        --------
        validate : The same checks, without the spec.

        Examples
        --------
        >>> spec = builder.build()
        >>> spec.name, spec.metadata["lake_ids"]
        ('lak', ['pond'])
        """

        for key, why in (
            ("nlakes", "it is the number of lakes read from `lakes`"),
            ("pname", "rename the package with name=, which also names its outputs, its "
                      "table children and the endpoint a mover reads"),
        ):
            if key in self.options:
                raise TypeError(
                    f"LAK option {key}={self.options[key]!r} is not for **options: {why} "
                    "(ledger 318)."
                )
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
            length_conversion=(
                self.length_conversion if self.length_conversion is not None
                else mf6_length_conversion(getattr(self.context, "length_units", "feet"))
            ),
            time_conversion=(
                self.time_conversion if self.time_conversion is not None
                else mf6_time_conversion(getattr(self.context, "time_units", "days"))
            ),
            maximum_iterations=self.maximum_iterations,
            maximum_stage_change=self.maximum_stage_change,
            **dict(self.options),
        )
        options = {**base.options, "lake_tables": tuple(child_tables)}
        if self.implicit:
            # Only when True: `False` adds no key, so no existing digest moves. MF6 writes
            # no package-convergence CSV under IMPLICIT (it warns instead), so the record
            # for it is not asked for (ledger 326).
            options["implicit"] = True
            options.pop("package_convergence_filerecord", None)
        return replace(
            base,
            builder=_build_lak_with_tables,
            options=options,
            metadata={
                **base.metadata,
                # LAK is invisible to the range check: `_cellids_in` deliberately
                # refuses to walk `packagedata`/`connectiondata`, because a cellid
                # sits at a different index per package (SFR's `packagedata[1]` is a
                # cellid, LAK's is a stage). So the pin is the ONLY thing that can
                # catch a lake resolved on another mesh -- measured, a lake declared
                # on 28 cells BUILT CLEAN on 92 (ledger 267).
                **pin_metadata(self.context, certifies=("connectiondata",)),
                "builder": "LAKBuilder",
                "lake_ids": list(self.lake_ids),
                # Which rule wrote each lake, and what each connection is (ledger 325).
                # Read by the build-time lake check and warnings; voided with the pin
                # when an edit replaces `connectiondata`.
                LAKE_RULES: dict(self.lake_rules),
                LAKE_CONNECTIONS: self.connection_roles,
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
