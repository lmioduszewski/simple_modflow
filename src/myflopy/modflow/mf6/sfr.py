"""Build MODFLOW 6 SFR package specifications from stream geometries."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from numbers import Integral, Real
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from myflopy.modflow.mf6.mvr import MoverConnection

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
from shapely.ops import nearest_points

from myflopy._length_units import mf6_length_conversion
from myflopy._time_units import mf6_time_conversion
from myflopy.advanced import sfr_spec
from myflopy.geopackage import _layer_index
from myflopy.grid_identity import pin_metadata
from myflopy.modflow.mf6.grid.geometry import _SLIVER_WIDTH, top_active_layer
from myflopy.modflow.mf6.surface_water_validation import _REACH_VALUES, _reach_value_fault
from myflopy.specs import (
    ModelContext,
    PackageSpec,
    _layer_elevations,
)

StreamLocation = str | Point


def _is_integer(value: Any) -> bool:
    """Whether ``value`` is an integer -- a Python ``int`` or a numpy one.

    ``isinstance(value, int)`` is False for ``numpy.int64`` (it subclasses nothing),
    so a layer or a key read out of an array or a DataFrame was taken for a column
    name or a stream id (the numpy-at-the-boundary class, ledger 265). A ``bool`` is
    not a layer, a reach or a period.
    """

    return isinstance(value, Integral) and not isinstance(value, bool)


#: Period settings that are a rate per unit of streambed area (``RAINFALL``,
#: ``EVAPORATION``, length per time) or a reach's state (``STATUS``): a stream's
#: value holds on EVERY one of its reaches, as a scalar's holds on every reach.
_PER_REACH_SETTINGS = frozenset({"RAINFALL", "EVAPORATION", "STATUS"})


@dataclass(frozen=True, slots=True)
class StreamConnection:
    """An explicit reach-to-reach link between two named streams.

    Pass these to ``mf.sfr(connections=(...))`` (or ``SFRBuilder``) to disambiguate
    confluences that automatic geometry-based topology cannot resolve -- e.g. two
    tributaries meeting a main stem. Water leaves ``source`` (by default at its
    downstream end) and enters ``receiver`` (by default at its nearest reach).

    An explicit connection replaces whatever connection the builder would have
    inferred for the same ``source``, in every ``connection_mode``. Both
    stream ids must exist in the builder's stream table, or
    :attr:`SFRBuilder.network` raises ``ValueError``; so does a set of
    connections that routes water in a cycle.

    Parameters
    ----------
    source : str
        Stream id of the stream the water leaves -- a value of the builder's
        ``stream_id`` column, or the id it generated (the file stem, or
        ``"stream_<n>"``).
    receiver : str or None
        Stream id of the stream the water enters. ``None`` marks ``source`` as a
        terminal outflow: no reach connection is written for it, and it is not
        reported as ``unresolved``.
    source_location : str or shapely.geometry.Point, optional
        Where on ``source`` the water leaves. Default ``"downstream"`` -- the
        end of the line (``"end"`` is the same). ``"upstream"``/``"start"`` is
        the first vertex; ``"nearest"``/``"automatic"`` is the point of
        ``source`` nearest ``receiver``; ``"intersection"`` is a point where the
        two lines cross (``ValueError`` if they do not); a ``Point`` is used as
        given. The source reach is the reach of ``source`` nearest that point.
        Any other string raises ``ValueError`` when the connection is resolved.
        A source reach above the stream's last one already flows into the next
        reach of its stream; the second downstream reach makes MF6's upstream
        fractions sum to two, so :meth:`SFRBuilder.validate` refuses it -- a
        partial split is a :class:`StreamDiversion` (ledger 366).
    receiver_location : str or shapely.geometry.Point, optional
        Where on ``receiver`` the water enters, with the same vocabulary.
        Default ``"nearest"``: the point ``source_location`` resolved to,
        projected onto the receiver's line -- so a tributary joins the main stem
        where it actually meets it. The receiving reach is the reach of
        ``receiver`` nearest the resolved point; a point on the boundary between
        two reaches goes to the downstream one.

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

    Each diversion becomes one MF6 ``DIVERSIONS`` record on the source reach
    (numbered in the order the diversions are given, per source reach), and one
    ``DIVERSION`` period setting per stress period that has an amount. The
    receiving reach is written with ``USTRF`` 0.0 (MF6's fraction of upstream
    flow a reach takes). Both stream ids must exist in the builder's stream
    table, or :attr:`SFRBuilder.network` raises ``ValueError``.

    MF6 diverts only along a reach connection, so the diversion is also written
    as one, from the source reach to the receiving reach (ledger 364) --
    once, when a :class:`StreamConnection` names the same two reaches too. The
    streams need no other link. The water the diversion does not take goes on
    to the source reach's other downstream reaches by their ``USTRF``, so MF6
    needs the source reach to HAVE one: the default ``source_location``, the
    stream's last reach, works only when the stream flows on into another, and
    :meth:`SFRBuilder.validate` refuses it otherwise, naming a ``Point`` above
    the last reach or an onward :class:`StreamConnection`. And because the
    receiving reach's ``USTRF`` is 0, no other reach may flow into it alone --
    divert into a reach nothing else flows into (ledger 366). The two locations
    must resolve to two DIFFERENT reaches: a diversion leaving and entering one
    reach would connect it to itself, which MF6 cannot read, and
    :meth:`SFRBuilder.validate` refuses it naming the reach.

    Parameters
    ----------
    source : str
        Stream id of the donor stream.
    receiver : str
        Stream id of the diversion channel.
    amount : float, str or Mapping[int, float or str]
        MF6's ``DIVFLOW``: how much to divert, read according to ``priority`` --
        a fraction of the reach's outflow for ``"FRACTION"``, a volumetric rate
        (model length cubed per model time) for the other three. A mapping is
        ``{period: amount}`` with zero-based periods, and a period missing from
        it gets no ``DIVERSION`` setting; anything else is written unchanged for
        every period. Not checked here.
    source_location : str or shapely.geometry.Point, optional
        Where on ``source`` the diversion leaves; default ``"downstream"``, the
        stream's last reach -- refused by :meth:`SFRBuilder.validate` unless the
        stream flows on into another, since MF6 diverts only from a reach with a
        downstream reach of its own (above). Same vocabulary and resolution as
        :class:`StreamConnection`.
    receiver_location : str or shapely.geometry.Point, optional
        Where on ``receiver`` the diverted water enters; default ``"nearest"``.
        Same vocabulary and resolution as :class:`StreamConnection`. The reach
        it names must be one no other reach flows into alone -- ``"upstream"``
        on a stream nothing joins -- or :meth:`SFRBuilder.validate` refuses it.
    priority : str, optional
        MF6 SFR ``CPRIOR``: ``"FRACTION"`` (default), ``"EXCESS"``,
        ``"THRESHOLD"`` or ``"UPTO"``, case-insensitive and written upper-case.
        Any other value makes :attr:`SFRBuilder.network` raise ``ValueError``.
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
    stream_ids : tuple of str
        Ids of all streams in the network, in stream-table order.
    connections : tuple of StreamConnection
        At most one connection per source stream: the ones the builder
        inferred, with any explicit ``connections=`` entry replacing the
        inferred one for the same source.
    diversions : tuple of StreamDiversion
        The builder's ``diversions=``, unchanged.
    unresolved : tuple of str
        Stream ids whose downstream link could not be inferred (no stream end
        within tolerance, or no matching ``from_node``) and that no explicit
        connection covers. Such a stream ends the network; give it an explicit
        connection -- ``receiver=None`` for a true outlet -- to say so.
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

    The builder is frozen and resolves lazily: the stream table, network and
    reaches are computed on first access and cached (the MF6 records are
    recomputed from them on each access), so a construction error in a field
    (a missing column, an unknown stream) surfaces when ``.reaches``,
    ``.network``, ``.packagedata`` or :meth:`build` first needs it. Lengths are
    model length units and times model time units throughout.

    Parameters
    ----------
    context : ModelContext
        Carries the grid the reaches are cut on (``context.grid``, required --
        ``ValueError`` at construction without one), the ``domain`` read by
        ``reach_layer="top_active"``, the layer elevations read by the default
        ``reach_top`` (``context.surfaces``, else ``grid.gdf_topbtm``), and the
        ``length_units``/``time_units`` the default unit conversions come from.
        The grid's pin is stamped on the built package.
    nper : int
        Number of stress periods; ``perioddata`` gets an entry for each period
        ``0 .. nper - 1`` that sets something (a period with no setting has none, ledger
        468). Less than 1 raises ``ValueError``. Period keys at or
        beyond ``nper`` in a period input are ignored.
    streams : str, pathlib.Path, sequence of str or Path, or geopandas.GeoDataFrame
        Stream centerlines, one ``LineString`` per stream (anything else, a
        ``MultiLineString`` included, raises ``ValueError``). A path or a list of
        paths is read with :func:`geopandas.read_file` and concatenated; the
        files must share one CRS (geopandas refuses to concatenate otherwise).
        The lines are reprojected to the grid's CRS when both carry one. A line
        runs upstream to downstream: its first vertex is the headwater and its
        last the outlet (see ``reverse_streams``).
    stream_id : str or None, optional
        Column holding each stream's id; values are turned into strings and must
        be unique (``ValueError`` otherwise, and when the column is missing).
        Default ``None``: each file's stem when ``streams`` is paths that each
        hold one feature, else ``"stream_0"``, ``"stream_1"``, ... by row.
    from_node : str or None, optional
        Column naming each stream's upstream node, read only by
        ``connection_mode="nodes"`` (where it is required). Default ``None``.
    to_node : str or None, optional
        Column naming each stream's downstream node, read only by
        ``connection_mode="nodes"``: a stream flows into the stream whose
        ``from_node`` equals its ``to_node``. No match leaves it unresolved; more
        than one raises ``ValueError``. Default ``None``.
    connection_mode : str, optional
        How the stream-to-stream topology is inferred. ``"automatic"`` (default)
        connects each stream's downstream end to the nearest other stream within
        ``connection_tolerance`` of it (a tie for nearest raises ``ValueError``;
        none leaves it unresolved). ``"nodes"`` reads ``from_node``/``to_node``.
        ``"explicit"`` infers nothing, so only ``connections`` link streams. Any
        other value raises ``ValueError`` at construction.
    connections : sequence of StreamConnection, optional
        Explicit links, stored as a tuple. Each replaces the inferred connection
        for its ``source`` in every mode; for two with the same source the
        later wins. Default ``()``.
    diversions : sequence of StreamDiversion, optional
        Explicit partial splits, written as MF6 ``DIVERSIONS`` records,
        ``DIVERSION`` period settings and a reach connection from the source
        reach to the receiving reach (see :class:`StreamDiversion`, including
        the two placements MF6 cannot run and :meth:`validate` refuses).
        Default ``()``.
    connection_tolerance : float or None, optional
        Largest distance, in model length units, between a stream's downstream
        end and another stream for ``connection_mode="automatic"`` to connect
        them. Default ``None``: 5% of the square root of the median cell area,
        and never below 1e-9. Unused by the other modes.
    reverse_streams : str or sequence of str, optional
        Stream ids whose line is reversed when the stream table is built, for
        lines digitized mouth-first. Ids are compared as strings. A single string is
        ONE id: ``"main"`` reverses stream ``main`` (it was iterated character by
        character). Default ``()``.
    reach_layer : int, float, str or Mapping[str, int], optional
        Zero-based layer each reach sits in. A number is every reach's layer; a
        mapping is ``{stream_id: layer}`` (a stream missing from it raises
        ``KeyError``); a string other than ``"top_active"`` names a stream-table
        column holding each stream's layer. Every layer so given must be a WHOLE
        number -- a Python or numpy integer, or a float such as ``1.0`` -- at or
        above 0 and, when ``context.domain`` is two-dimensional, below its layer
        count. A whole number written as text (``"1"``, a text column) is read
        as that number, as the GeoPackage builders read a layer column; a
        fraction, a ``bool``, a NaN, text that is no number, or anything else
        that is not a number raises ``ValueError`` naming the value when the
        reaches are resolved. Default
        ``"top_active"``: the top layer with ``idomain > 0`` in the reach's column
        of ``context.domain`` -- a column with none raises ``ValueError`` -- or
        layer 0 when the context has no domain or a one-dimensional one.
    width : float, str, Mapping, or sequence, optional
        Reach width, MF6 ``RWID``, in model length units; default 10.0. Any of
        the per-reach shapes in the Notes. :meth:`validate` refuses a width
        that is not positive.
    gradient : float, str, Mapping, or sequence, optional
        Streambed slope, MF6 ``RGRD`` (dimensionless); default 0.001. Per-reach
        shapes as in the Notes. :meth:`validate` refuses a negative gradient.
    reach_top : float, str, Mapping, sequence, or None, optional
        Streambed-top elevation, MF6 ``RTP``, in model length units. Per-reach
        shapes as in the Notes. Default ``None``: the top of each reach's layer
        in its cell -- read from ``context.surfaces``, else
        ``grid.gdf_topbtm`` (``ValueError`` when neither is available) -- then
        lowered where needed so no reach is higher than the reach above it in
        the same stream. Any finite number; MF6 takes no time series here.
    roughness : float, str, Mapping, or sequence, optional
        Manning's roughness coefficient, MF6 ``MAN``; default 0.03. Per-reach
        shapes as in the Notes, or a time-series name. :meth:`validate` refuses
        a value that is not positive.
    streambed_k : float, str, Mapping, or sequence, optional
        Streambed hydraulic conductivity, MF6 ``RHK``, in model length per model
        time; default 1.0. Per-reach shapes as in the Notes, or a time-series
        name. :meth:`validate` refuses a negative value.
    streambed_thickness : float, str, Mapping, or sequence, optional
        Streambed thickness, MF6 ``RBTH``, in model length units; default 1.0.
        Per-reach shapes as in the Notes. :meth:`validate` refuses a value that
        is not positive.
    inflow : float, str, Mapping, sequence, or None, optional
        Volumetric inflow, MF6 ``INFLOW`` (model length cubed per model time).
        Period shapes as in the Notes; per period, a number or a time-series
        name enters at the FIRST reach of the reach table only -- the headwater
        of the first stream -- not at every stream's headwater, while
        ``{stream_id: q}`` enters at that stream's first reach. Default
        ``None``, no inflow.
    rainfall : float, str, Mapping, sequence, or None, optional
        Rainfall rate on the reach surface, MF6 ``RAINFALL`` (model length per
        model time -- a rate per unit of streambed area, so every reach takes
        it unchanged). Per period, a number or time-series name applies to
        EVERY reach, and ``{stream_id: rate}`` to every reach of that stream; a
        reach-number key names one reach and beats its stream's key. Other
        period shapes as in the Notes. Default ``None``.
    evaporation : float, str, Mapping, sequence, or None, optional
        Evaporation rate from the reach surface, MF6 ``EVAPORATION`` (model
        length per model time). Placed exactly like ``rainfall``: a scalar on
        every reach, ``{stream_id: rate}`` on every reach of that stream.
        Default ``None``.
    runoff : float, Mapping, sequence, or None, optional
        Overland runoff, MF6 ``RUNOFF`` (model length cubed per model time). Per
        period, a number is EACH stream's total and ``{stream_id: total}`` names
        one stream's (an unknown stream raises ``KeyError``); either is split
        over the stream's reaches in proportion to reach length, so refining the
        grid does not multiply it. An integer key names a reach and takes its
        amount as given, as do explicit rows. A time-series name cannot be split
        and raises ``ValueError`` unless given as rows. Default ``None``.
    status : str, Mapping, sequence, or None, optional
        MF6 reach ``STATUS``: ``"ACTIVE"``, ``"INACTIVE"`` or ``"SIMPLE"``. Per
        period, a string applies to every reach; ``{stream_id: status}`` to every
        reach of that stream, placed like ``rainfall``. Other period shapes as in
        the Notes. Default ``None``, no status setting (MF6's default, active).
    name : str, optional
        Package name, default ``"sfr"``: FloPy's ``pname``, the input file
        ``<model>.<name>``, the output files ``<model>_<name>_*``, and the
        package name :meth:`connection` puts in each mover endpoint.
    mover : bool, optional
        Write MF6's ``MOVER`` option so MVR can take water from and give water
        to these reaches. Default ``False``.
    length_conversion : float or None, optional
        MF6 ``LENGTH_CONVERSION`` for Manning's equation (meters to model length
        units). Default ``None``: from ``context.length_units``, read
        case-insensitively -- ``"feet"`` 3.28081, ``"meters"`` 1.0,
        ``"centimeters"`` 100.0 (and every other spelling a ``LayerStack``
        takes: ``"ft"``, ``"metres"``, ``"cm"``, ...), ``None``, ``""`` or
        ``"unknown"`` 1.0. ``ModelContext`` refuses any other unit at
        construction (ledger 433).
    time_conversion : float or None, optional
        MF6 ``TIME_CONVERSION`` (seconds to model time units). Default ``None``:
        from ``context.time_units``, read case-insensitively -- ``"days"``
        86400, ``"hours"`` 3600, ``"minutes"`` 60, ``"seconds"`` 1, ``"years"``
        31557600 (singulars and the short forms ``"d"``, ``"hr"``, ``"min"``,
        ``"s"``, ``"yr"``, ... too), ``None``, ``""`` or ``"unknown"`` 1.0.
        ``ModelContext`` refuses any other unit at construction (ledger
        456).
    min_reach_length : float, optional
        Shortest piece of a stream, in model length units, that becomes a
        reach. A piece in one cell no longer than ``max(0.01, min_reach_length)``
        is dropped -- MF6 divides by reach length, and a sliver reach can crash
        it while reading the package -- and the stream stays connected through
        its other reaches. A stream left with no reach raises ``ValueError``.
        Default 0.0, which drops only pieces of 0.01 or less.
    maximum_picard_iterations : int, optional
        MF6 SFR ``MAXIMUM_PICARD_ITERATIONS``, written as given. Default 1.
    maximum_iterations : int, optional
        MF6 SFR ``MAXIMUM_ITERATIONS`` (reach-solution iterations), written as
        given. Default 1000.
    maximum_depth_change : float, optional
        MF6 SFR ``MAXIMUM_DEPTH_CHANGE`` (model length units), written as given.
        Default 0.01.
    options : Mapping[str, Any], optional
        Extra keyword arguments for :class:`flopy.mf6.ModflowGwfsfr`, passed
        through :func:`myflopy.advanced.sfr_spec`. They override the defaults
        that function sets (``save_flows=True``, the stage/budget output file
        records). A key :meth:`build` already passes -- ``packagedata``,
        ``connectiondata``, ``perioddata``, ``name``, ``mover``, ``diversions``,
        the unit conversions or the three iteration options -- raises
        ``TypeError`` at build. Default empty.

    Notes
    -----
    **Per-reach inputs** (``width``, ``gradient``, ``reach_top``, ``roughness``,
    ``streambed_k``, ``streambed_thickness``) take:

    - a number, applied to every reach;
    - a string naming a stream-table column, read per stream. Any other string
      is written unchanged to every reach as an MF6 time-series name, which
      MF6's SFR takes for ``streambed_k`` (``RHK``) and ``roughness`` (``MAN``)
      only: for the other four :meth:`validate` raises ``ValueError`` naming the
      argument, the reach and its stream;
    - a mapping keyed by reach number, stream id or ``(layer, cell)`` cellid,
      looked up in that order for each reach -- a reach none of them matches
      has no value, and :meth:`validate` raises ``ValueError`` naming it;
    - a sequence of one value per reach (in reach order), or one per stream in
      :attr:`stream_ids` order, which puts source streams before the streams
      they flow into and can differ from the file's order. Any other length
      raises ``ValueError``.

    **Period inputs** (``inflow``, ``rainfall``, ``evaporation``, ``runoff``,
    ``status``) take ``None`` (no setting), or a mapping whose keys are ALL
    integers, read as ``{period: value}`` with zero-based periods -- a period
    left out gets no setting -- or any other value, used for every period. So
    ``{3: 100.0}`` means period 3, not reach 3: key a reach inside a period,
    ``{0: {3: 100.0}}``. Each period's value is a number or a string (placed as
    each parameter says), a mapping ``{stream_id or reach number: value}`` (a
    stream id means every reach of the stream for ``rainfall``, ``evaporation``
    and ``status``, its first reach for ``inflow``, and a share by reach length
    for ``runoff``; a stream not in the table raises ``KeyError`` naming the
    streams), or an iterable of explicit rows ``(reach, value, ...)`` written as
    given. Keys may be Python or numpy integers.

    Every per-reach value must end as a finite number, or a string where MF6
    takes a time series; anything else -- a ``NaN`` in a stream-table column
    included -- makes :meth:`validate` (and so :meth:`build`) raise
    ``ValueError`` naming the argument, the reach and its stream.
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
    reverse_streams: str | tuple[str, ...] = ()
    reach_layer: int | float | str | Mapping[str, int] = "top_active"
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
    min_reach_length: float = 0.0
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
        # A bare string is ONE id, as at `usg.to_mf6(include=)` (ledger 429): iterated, "main"
        # reversed streams "m", "a", "i" and "n" -- none of them, in silence.
        reverse = (
            (self.reverse_streams,) if isinstance(self.reverse_streams, str) else self.reverse_streams
        )
        object.__setattr__(self, "reverse_streams", tuple(str(value) for value in reverse))

    def with_updates(self, **updates: Any) -> SFRBuilder:
        """Return a changed builder without modifying the original.

        Parameters
        ----------
        **updates
            Constructor fields to replace, by name (``width=20.0``,
            ``streams=other``), passed to :func:`dataclasses.replace`. Every other
            field is carried over; the cached stream table, network and reaches
            are not, so the new builder resolves afresh. The construction checks
            run again, and a name that is not a constructor field raises
            ``TypeError``.

        Returns
        -------
        SFRBuilder
        """

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

        domain = self.context.domain
        if isinstance(self.reach_layer, Mapping):
            return self._layer_number(self.reach_layer[stream_id], f"reach_layer[{stream_id!r}]")
        if isinstance(self.reach_layer, str) and self.reach_layer != "top_active":
            row = self.stream_table.loc[stream_id]
            if self.reach_layer not in row:
                raise ValueError(f"Reach-layer field is missing: {self.reach_layer}")
            return self._layer_number(
                row[self.reach_layer],
                f"reach_layer column {self.reach_layer!r} for stream {stream_id!r}",
            )
        if not isinstance(self.reach_layer, str):
            return self._layer_number(self.reach_layer, "reach_layer")
        if domain is None:
            return 0
        values = np.asarray(domain)
        if values.ndim == 1:
            return 0
        # The shared rule, `idomain > 0` (ledger 319), asked of this one column.
        layer = int(top_active_layer(values[:, [cell]])[0])
        if layer < 0:
            raise ValueError(f"Stream '{stream_id}' intersects inactive grid column {cell}.")
        return layer

    @property
    def reaches(self) -> gpd.GeoDataFrame:
        """Return the inspectable grid-dependent SFR reach table."""

        if self._reaches is not None:
            return self._reaches

        rows = []
        rno = 0
        polys = self.grid.gdf_vorPolys
        for stream_id in self._ordered_stream_ids():
            line = self.stream_table.loc[stream_id].geometry
            # Lowest-numbered first: the cell that keeps a stretch lying on a face.
            intersecting = polys[polys.intersects(line)].sort_index()
            segments = []
            claimed = None
            for cell, grid_row in intersecting.iterrows():
                segment = grid_row.geometry.intersection(line)
                # A closed cell includes its boundary, so a stretch lying ON a face is
                # in both cells' pieces. Counted in both, it doubled rlen and zig-zagged
                # the chain between them. The lower-numbered cell keeps it (decided
                # 2026-09-11): halving it would make two parallel reaches.
                if claimed is not None:
                    segment = segment.difference(claimed)
                # A stream clipping the corner of a cell yields a sliver reach. MF6
                # divides by reach length, and a reach of 0.005 ft crashed it outright
                # (SIGFPE while reading the package) on a real 652-reach model, where
                # 36 reaches were under 5 ft. Dropping one leaves the stream connected
                # through its neighbours, because the topology is rebuilt from the
                # ordering below rather than from the discarded cell.
                # Under ``_SLIVER_WIDTH`` a piece is noise whatever ``min_reach_length``
                # says: a CRS round trip clipped corner cells by ~9e-7 ft.
                if segment.length <= max(_SLIVER_WIDTH, self.min_reach_length):
                    continue
                # Claimed only once KEPT, so a lower cell whose piece is dropped does
                # not take the stretch from the cell that would have kept it.
                along = segment.intersection(grid_row.geometry.boundary)
                if along.length > 0:
                    claimed = along if claimed is None else claimed.union(along)
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
        """Return the final generated reach for one stream.

        Parameters
        ----------
        stream_id : str
            The stream's id (compared as a string). An id with no reaches raises
            ``KeyError``.

        Returns
        -------
        int
            Zero-based reach number of the stream's last (most downstream) reach.
        """

        return self.stream_reaches[str(stream_id)][-1]

    def connection(self, stream_id: str, location: StreamLocation = "downstream") -> MoverConnection:
        """Return an MVR endpoint for a stream by stable stream ID and location.

        Parameters
        ----------
        stream_id : str
            The stream's id (compared as a string). An unknown id raises
            ``KeyError``.
        location : str or shapely.geometry.Point, optional
            Which reach of the stream. Default ``"downstream"`` (or ``"end"``): the
            reach nearest the line's last vertex. ``"upstream"``/``"start"``: the
            reach nearest its first vertex. A ``Point``: the reach nearest it, the
            downstream one when it sits on the boundary between two. The
            two-stream words ``"nearest"``, ``"automatic"`` and ``"intersection"``
            need a second stream and raise ``ValueError`` here, as does any other
            string.

        Returns
        -------
        MoverConnection
            ``MoverConnection(<builder name>, <reach number>)``, for ``mf.mvr``.
        """

        from myflopy.modflow.mf6.mvr import MoverConnection

        return MoverConnection(self.name, self._reach_at(str(stream_id), location))

    def _point_for_location(
        self,
        stream_id: str,
        location: StreamLocation,
        *,
        other: str | None = None,
        toward: Point | None = None,
    ) -> Point:
        """Resolve a stream location keyword/Point to a coordinate on the stream line.

        Supports an explicit ``Point``, ``upstream``/``downstream`` ends, and
        ``nearest``/``intersection`` against ``other`` (required for those).
        ``nearest`` is where ``other`` MEETS this stream: ``toward`` -- the point
        ``other`` itself resolved to, a tributary's outlet say -- projected onto
        this line, or without it the nearest point between the two lines.
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
            # This projected ``other_line.representative_point()`` -- for a two-vertex
            # tributary its upstream START -- so a tributary joining the main stem
            # mid-grid was connected to the main stem's LAST reach.
            anchor = toward if toward is not None else nearest_points(other_line, line)[0]
            return line.interpolate(line.project(anchor))
        raise ValueError(f"Unsupported stream location: {location!r}")

    def _reach_at(
        self,
        stream_id: str,
        location: StreamLocation,
        *,
        other: str | None = None,
        toward: Point | None = None,
    ) -> int:
        """The reach number on ``stream_id`` at the resolved ``location`` point.

        A point on the boundary between two reaches is as near one as the other.
        The DOWNSTREAM reach takes it: water joining there enters that reach at its
        upstream end, where the upstream reach would carry it along a stretch above
        the confluence.
        """

        point = self._point_for_location(stream_id, location, other=other, toward=toward)
        reaches = self.reaches[self.reaches["stream_id"] == stream_id]
        distances = reaches.geometry.distance(point)
        tied = distances[distances <= distances.min() + _SLIVER_WIDTH].index
        return int(reaches.loc[tied, "sequence"].idxmax())

    def _layer_number(self, value: Any, where: str) -> int:
        """One ``reach_layer`` value as a zero-based layer, or ``ValueError`` naming it.

        The layer rule the GeoPackage builders read a layer column by
        (``geopackage._layer_index``): a WHOLE number -- ``1``, ``numpy.int64(1)``,
        ``1.0`` -- at or above 0 and, where the context's two-dimensional ``domain``
        says how many layers there are, below that. A float fell through to the
        column-name branch ("Reach-layer field is missing: 1.0"), a fraction in a
        mapping or a column was truncated by ``int()`` into a layer nobody wrote, and
        a negative layer or one past the domain was written. A ``bool`` is not a
        layer (ledger 265's numeric-scalar rule; :func:`_is_integer`; ledger 412).
        """

        domain = self.context.domain
        nlay = int(np.shape(domain)[0]) if domain is not None and np.ndim(domain) == 2 else None
        if isinstance(value, (bool, np.bool_)):
            why = f"{value!r} is a bool, not a layer number"
        else:
            layer, why = _layer_index(value, base=0, nlay=nlay)
            if layer is not None:
                return layer
        raise ValueError(
            f"SFR {where} names no model layer: {why}. A layer is a zero-based whole number."
        )

    @property
    def resolved_connections(self) -> tuple[tuple[int, int], ...]:
        """Return source/receiver reach pairs for ordinary stream connections."""

        pairs = []
        for connection in self.network.connections:
            if connection.receiver is None:
                continue
            outlet = self._point_for_location(
                connection.source, connection.source_location, other=connection.receiver
            )
            source = self._reach_at(connection.source, outlet)
            receiver = self._reach_at(
                connection.receiver,
                connection.receiver_location,
                other=connection.source,
                toward=outlet,
            )
            pairs.append((source, receiver))
        return tuple(pairs)

    @property
    def resolved_diversions(self) -> tuple[tuple[StreamDiversion, int, int, int], ...]:
        """Return diversion definitions with generated reach and diversion IDs."""

        counts: dict[int, int] = {}
        resolved = []
        for diversion in self.diversions:
            outlet = self._point_for_location(
                diversion.source, diversion.source_location, other=diversion.receiver
            )
            source = self._reach_at(diversion.source, outlet)
            receiver = self._reach_at(
                diversion.receiver, diversion.receiver_location, other=diversion.source, toward=outlet
            )
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
            for source, receiver in zip(reaches, reaches[1:], strict=False):
                downstream[source].append(receiver)
                upstream[receiver].append(source)
        for source, receiver in self.resolved_connections:
            downstream[source].append(receiver)
            upstream[receiver].append(source)
        # MF6 diverts only along a reach connection ("Reach 3 is not a upstream reach
        # for reach 5 ... diversion 1 from reach 3 to reach 5 is not possible"), so a
        # diversion is one too (ledger 364). The sets below write it once when a
        # StreamConnection already names the same pair.
        for _, source, receiver, _ in self.resolved_diversions:
            downstream[source].append(receiver)
            upstream[receiver].append(source)

        # A downstream reach is written NEGATIVE, and reach 0 has no negative int:
        # FloPy reads the float -0.0 as "downstream to reach 0" (its
        # support_negative_index rule) and a plain 0 as "upstream from reach 0".
        # Nothing else connects INTO reach 0 -- streams are numbered sources first --
        # but a diversion's receiver can be reach 0 (ledger 364).
        return [
            [
                int(rno),
                *sorted(set(upstream[rno])),
                *[-0.0 if value == 0 else -int(value) for value in sorted(set(downstream[rno]))],
            ]
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
            by_stream = dict(zip(self.stream_ids, values, strict=False))
            return [by_stream[stream_id] for stream_id in self.reaches["stream_id"]]
        raise ValueError(f"{name} must contain one value per reach or stream.")

    def _default_reach_tops(self) -> list[float]:
        """Per-reach streambed-top elevations from the grid surface, enforced monotone downstream.

        Samples each reach cell's layer top, then clamps so a reach is never higher
        than the reach above it in the same stream.
        """

        # The shared resolver, not `self.grid.gdf_topbtm`: SFR used to be the one
        # builder that ignored `context.surfaces`, so it disagreed with LAK and
        # CellSurfaceOffset in the same model whenever the context was stale.
        surfaces = _layer_elevations(self.context)
        if surfaces is None:
            raise ValueError("reach_top is required when grid surface elevations are unavailable.")
        # Columns by POSITION among the numeric ones ([top, botm0, botm1, ...]),
        # so a unit-labelled table reads the same as an integer-labelled one.
        elevations = surfaces.select_dtypes("number")
        values = []
        for layer, cell in self.reaches["cellid"]:
            top_column = 0 if layer == 0 else layer
            values.append(float(elevations.loc[cell].iloc[top_column]))
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

        Accepts a per-period mapping, a per-reach/location mapping, a scalar, or
        explicit row tuples.

        A **scalar** is expanded by keyword: ``RAINFALL``/``EVAPORATION``/
        ``STATUS`` are per-reach quantities and broadcast to EVERY reach
        (matching :class:`LAKBuilder`'s per-lake broadcast), while ``INFLOW`` is
        a volumetric point source and enters at the single headwater (first)
        reach. ``RUNOFF`` is volumetric too, but spread: see :meth:`_runoff_rows`.
        A ``{stream_id: value}`` mapping follows the same rule per stream: see
        :meth:`_per_reach_rows`.
        """

        if value is None:
            return {period: [] for period in range(self.nper)}
        period_values = value if isinstance(value, Mapping) and all(_is_integer(key) for key in value) else {period: value for period in range(self.nper)}
        result = {period: [] for period in range(self.nper)}
        for period in range(self.nper):
            current = period_values.get(period)
            if current is None:
                continue
            if keyword == "RUNOFF" and (isinstance(current, (Mapping, Real, str))):
                result[period].extend(self._runoff_rows(current))
                continue
            if isinstance(current, Mapping) and keyword in _PER_REACH_SETTINGS:
                result[period].extend(self._per_reach_rows(current, keyword))
            elif isinstance(current, Mapping):
                for location, amount in current.items():
                    if _is_integer(location):
                        rno = int(location)
                    else:
                        self._reaches_of(str(location), keyword)
                        rno = self._reach_at(str(location), "upstream")
                    result[period].append([rno, keyword, amount])
            elif isinstance(current, Real) or isinstance(current, str):
                if keyword == "INFLOW":
                    # Volumetric point source: a scalar enters at the single
                    # headwater (first) reach. Broadcasting it would inject the
                    # inflow once per reach.
                    result[period].append(
                        [int(self.reaches.index[0]), keyword, current]
                    )
                else:
                    # Per-reach quantity (rate or status): a scalar applies to
                    # EVERY reach. Previously it landed on the first reach only,
                    # so e.g. rainfall=0.001 silently wetted one reach.
                    for rno in self.reaches.index:
                        result[period].append([int(rno), keyword, current])
            else:
                for row in current:
                    result[period].append([int(row[0]), keyword, *row[1:]])
        return result

    def _reaches_of(self, stream: str, keyword: str) -> list[int]:
        """``stream``'s reach numbers, in order; ``KeyError`` naming the streams when unknown."""

        if stream not in self.stream_reaches:
            raise KeyError(
                f"{keyword.lower()} names stream {stream!r}; the streams are "
                f"{sorted(self.stream_reaches)}."
            )
        return self.stream_reaches[stream]

    def _per_reach_rows(self, current: Mapping, keyword: str) -> list[list[Any]]:
        """One period's ``{stream_id or reach: value}`` for a per-reach setting.

        ``RAINFALL`` and ``EVAPORATION`` are rates per unit of streambed area (L/T)
        and ``STATUS`` a reach's state, so a stream's value is written to EVERY one
        of its reaches, unchanged -- not divided, which is RUNOFF's rule for a
        volume. It went to the stream's FIRST reach only (ledger 363), so
        ``rainfall={0: {"main": 0.001}}`` wetted one reach where ``rainfall=0.001``
        wets them all. A reach-number key names that reach and beats its stream's
        key whichever comes first in the mapping: the more specific value wins, as
        in the per-reach inputs' lookup. One row per reach, in reach order.
        """

        by_reach: dict[int, Any] = {}
        for location, amount in current.items():
            if not _is_integer(location):
                for rno in self._reaches_of(str(location), keyword):
                    by_reach[int(rno)] = amount
        for location, amount in current.items():
            if _is_integer(location):
                by_reach[int(location)] = amount
        return [[rno, keyword, by_reach[rno]] for rno in sorted(by_reach)]

    def _runoff_rows(self, current: Any) -> list[list[Any]]:
        """One period's RUNOFF rows: a TOTAL per stream, spread along it by reach length.

        RUNOFF is a volume per time on each reach, and a scalar was written to
        EVERY reach -- so a stream's runoff scaled with its reach count, and
        refining the grid under a stream multiplied it. A number is now each
        stream's total, ``{stream: total}`` gives each named stream its own, and
        either is split in proportion to ``rlen``. An integer key names a REACH and
        takes its amount as given. A time series cannot be split by length, so it
        is refused unless given per reach, as explicit rows.
        """

        def refuse():
            raise ValueError(
                "A RUNOFF time series cannot be spread along a stream by reach length; "
                "give it per reach -- runoff={period: [(reach, 'series'), ...]} -- or "
                "give each stream a number."
            )

        rows: list[list[Any]] = []
        totals: dict[str, Any] = {}
        if isinstance(current, Mapping):
            for location, amount in current.items():
                if _is_integer(location):
                    rows.append([int(location), "RUNOFF", amount])
                else:
                    totals[str(location)] = amount
        elif isinstance(current, str):
            refuse()
        else:
            totals = {stream: current for stream in self.stream_reaches}

        for stream, total in totals.items():
            if isinstance(total, str):
                refuse()
            reaches = self._reaches_of(stream, "RUNOFF")
            lengths = self.reaches.loc[reaches, "rlen"].astype(float)
            for rno, length in zip(reaches, lengths, strict=True):
                rows.append([int(rno), "RUNOFF", float(total) * float(length) / float(lengths.sum())])
        return rows

    @property
    def perioddata(self) -> dict[int, list[list[Any]]]:
        """MF6 SFR perioddata: ``{period: [record, ...]}`` for each period that sets something.

        A period with no setting has NO entry (ledger 468): MF6 keeps every reach
        setting until a period block changes it, so an empty block is the same model
        -- measured, budgets and stages byte-identical -- while FloPy 3.10 reads the
        package's FIRST period for its dtype, and an empty one made ``gwf.plot()`` and
        the package's ``.plot()`` raise. Every period empty gives ``{}``.
        """

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
        return {period: rows for period, rows in result.items() if rows}

    @property
    def diversiondata(self) -> list[list[Any]] | None:
        """MF6 SFR ``diversions`` records (source reach, index, receiver, priority), or ``None``."""

        rows = [
            [source, idiv, receiver, diversion.priority.upper()]
            for diversion, source, receiver, idiv in self.resolved_diversions
        ]
        return rows or None

    def validate(self) -> None:
        """Raise clear errors for invalid reach, routing, and package data.

        Every per-reach value must be a finite number (``ValueError`` naming the
        argument, the reach and its stream otherwise) -- ``streambed_k`` and
        ``roughness`` may instead be a string, read by MF6 as a time-series name,
        since those are the only two its SFR ``packagedata`` takes one for (ledger
        365). Then the bounds: positive length, width, streambed thickness and
        roughness; non-negative gradient and streambed conductivity; ``reach_top``
        any finite number. Every period setting must name a reach that exists.
        Last, the routing MF6 checks before it runs: every reach's downstream
        reaches must carry upstream fractions (``USTRF``) summing to one, which
        refuses a diversion from a reach with no downstream reach of its own, a
        diversion into a reach another reach flows into alone, and a
        :class:`StreamConnection` leaving a stream above its last reach (ledger
        366) -- after a diversion whose source and receiving locations resolve to
        ONE reach, which would connect that reach to itself.
        """

        if self.total_nreaches == 0:
            raise ValueError("SFRBuilder generated no reaches.")
        stream_of = self.reaches["stream_id"].to_dict()
        for row in self.packagedata:
            rno = int(row[0])
            where = f"reach {rno} (stream {stream_of[rno]!r})"
            if row[2] <= 0:
                raise ValueError(f"SFR {where} has nonpositive length.")
            for position, name, mf6_name, takes_series, bound in _REACH_VALUES:
                value = row[position]
                fault = _reach_value_fault(value, takes_series, bound)
                if fault == "number":
                    self._check_reach_number(value, name, mf6_name, takes_series, where)
                if fault == "positive":
                    raise ValueError(f"SFR {where} has nonpositive {name} ({mf6_name}): {value!r}.")
                if fault == "nonnegative":
                    raise ValueError(f"SFR {where} has negative {name} ({mf6_name}): {value!r}.")
        valid_reaches = set(self.reaches.index)
        for period, settings in self.perioddata.items():
            for setting in settings:
                if setting[0] not in valid_reaches:
                    raise ValueError(f"SFR period {period} references unknown reach {setting[0]}.")
        self._check_diversions_leave_their_reach(stream_of)
        self._check_upstream_fractions(stream_of)

    def _check_diversions_leave_their_reach(self, stream_of: Mapping[int, Any]) -> None:
        """Refuse a :class:`StreamDiversion` whose source and receiving reach are one reach.

        The diversion is written as a reach connection (ledger 364), so it connected
        the reach to itself, and MF6 stops reading the package on that: "Could not
        fill SFR package IA and JA connection data" (measured on a headwater reach
        and on a canal's first reach, ``docs/reproducers/sfr_stream_inputs.py`` E).
        On a reach something else flows into, the upstream-fraction check refused
        it first, blaming the reach above.
        """

        for diversion, source, receiver, _ in self.resolved_diversions:
            if source != receiver:
                continue
            raise ValueError(
                f"SFR StreamDiversion({diversion.source!r}, {diversion.receiver!r}) leaves and "
                f"enters the same reach, reach {source} (stream {stream_of[source]!r}): its "
                f"source_location={diversion.source_location!r} and "
                f"receiver_location={diversion.receiver_location!r} resolve to one reach, and MF6 "
                "cannot connect a reach to itself (it stops reading the package: 'Could not fill "
                "SFR package IA and JA connection data'). Give the two locations different reaches."
            )

    def _check_upstream_fractions(self, stream_of: Mapping[int, Any]) -> None:
        """Refuse a network MF6's ``sfr_check_ustrf`` stops on, naming the reaches.

        MF6 requires, for every reach with a downstream connection, that the
        ``USTRF`` of its downstream reaches sum to one -- diversion receivers
        included -- and a diversion's receiving reach to carry ``USTRF`` 0 (which
        this builder writes). So a diversion needs its source reach to have a
        downstream reach of its own, and its receiving reach to have no other
        reach flowing into it alone. Neither was checked, and MF6 stopped with
        "Upstream fractions for reach N is not equal to one (0.0)" -- for the
        DEFAULT ``source_location="downstream"`` on a stream with no outlet
        connection among them. The other way out, a receiving ``USTRF`` of 1,
        crashes MF6 6.8 in its own warning (``gwf-sfr.f90``, "End of record"),
        measured in ``docs/reproducers/sfr_stream_inputs.py`` D (ledger 366).
        """

        ustrf = {int(row[0]): float(row[10]) for row in self.packagedata}
        diverted_from: dict[int, list[int]] = {}
        for _, source, receiver, _ in self.resolved_diversions:
            diverted_from.setdefault(int(receiver), []).append(int(source))

        def named(rno: int) -> str:
            return f"reach {rno} (stream {stream_of[rno]!r})"

        for row in self.connectiondata:
            rno = int(row[0])
            downstream = [int(-value) for value in row[1:] if value < 0 or np.signbit(value)]
            if not downstream:
                continue
            total = sum(ustrf[reach] for reach in downstream)
            if abs(total - 1.0) <= 1.0e-6:
                continue
            receivers = [reach for reach in downstream if reach in diverted_from]
            if total == 0.0 and receivers and all(rno in diverted_from[reach] for reach in downstream):
                raise ValueError(
                    f"SFR diversion from {named(rno)} to "
                    f"{', '.join(named(reach) for reach in downstream)}: MF6 diverts only from a "
                    "reach that has a downstream reach of its own besides the receiving one, "
                    f"for the water the diversion does not take, and {named(rno)} has none -- it "
                    f"is the last reach of {stream_of[rno]!r}, which flows into no other stream. "
                    "Give the StreamDiversion a source_location above the stream's last reach "
                    f"(a Point on it), or connect {stream_of[rno]!r} onward with a "
                    "StreamConnection (ledger 366)."
                )
            if total == 0.0 and receivers:
                sources = sorted({source for reach in receivers for source in diverted_from[reach]})
                raise ValueError(
                    f"SFR {', '.join(named(reach) for reach in receivers)} receives a diversion from "
                    f"{', '.join(named(source) for source in sources)} and is also the only reach "
                    f"that {named(rno)} flows into: MF6 takes the upstream fraction (USTRF) of a "
                    f"diversion's receiving reach as 0, which leaves the water of {named(rno)} "
                    "nowhere to go. "
                    "Divert into a reach nothing else flows into -- receiver_location='upstream' "
                    "on a stream no other stream joins, or a Point on that reach "
                    "(ledger 366)."
                )
            raise ValueError(
                f"SFR {named(rno)} flows into "
                f"{', '.join(named(reach) for reach in downstream)}, whose upstream fractions "
                f"(USTRF) sum to {total:g}; MF6 requires one. A StreamConnection leaving a "
                "stream above its last reach gives that reach a second downstream reach -- "
                "use a StreamDiversion for a partial split (ledger 366)."
            )

    def _check_reach_number(
        self, value: Any, name: str, mf6_name: str, takes_series: bool, where: str
    ) -> None:
        """Refuse a per-reach value that is not a finite number, saying where it came from.

        A string that names no stream-table column was written unchanged to every
        reach and a mapping that matched none of a reach's keys left it ``None``:
        both reached ``<= 0`` as a bare ``TypeError``, and ``reach_top`` -- which
        nothing compared -- went to FloPy as written (ledger 365).
        """

        if isinstance(value, Real) and not isinstance(value, bool) and np.isfinite(value):
            return
        given = getattr(self, name)
        series = (
            ""
            if takes_series
            else (
                f"; MF6 reads {mf6_name} as a number -- of the per-reach values only "
                "streambed_k (RHK) and roughness (MAN) take a time-series name"
            )
        )
        if isinstance(given, str) and given in self.stream_table.columns:
            why = f"column {given!r} of the stream table holds {value!r} there{series}"
        elif given is None and name == "reach_top":
            why = f"the top of its layer in its cell, the default, is {value!r}"
        elif value is None and isinstance(given, Mapping):
            why = (
                f"the {name}= mapping names none of its reach number, its stream id "
                "or its (layer, cell) cellid"
            )
        elif isinstance(value, str):
            columns = [str(column) for column in self.stream_table.columns if column != "geometry"]
            why = f"{value!r} names no column of the stream table (its columns are {columns}){series}"
        elif value is None or isinstance(value, Real):
            why = f"{value!r} is not a finite number"
        else:
            why = f"{value!r} is a {type(value).__name__}, not a number"
        raise ValueError(f"SFR {name} ({mf6_name}) for {where}: {why}.")

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
            **pin_metadata(self.context, certifies=("packagedata",)),
            builder="SFRBuilder",
            stream_ids=list(self.stream_ids),
            connections=[
                {"source": item.source, "receiver": item.receiver}
                for item in self.network.connections
            ],
            unresolved=list(self.network.unresolved),
        )


__all__ = ["SFRBuilder", "StreamConnection", "StreamDiversion", "StreamNetwork"]
