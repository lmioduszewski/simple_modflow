"""RIV exchange-flow observation targets (`RivFlowTargets` + bound view).

A river cell's MF6 ``RIV`` flow is the exchange between the aquifer and the
river: NEGATIVE where the aquifer feeds the river (a gaining reach -- baseflow),
POSITIVE where the river leaks into the aquifer. The zone machinery is shared with
``DrnFlowTargets`` in :mod:`myflopy.modflow.mf6.observations._zones` (ledger
534); this module names what is RIV's own.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from myflopy.modflow.mf6.observations._shared import _ObservedQuantity
from myflopy.modflow.mf6.observations._zones import (
    _BoundCellZoneFlowTargets,
    _CellZoneFlowTargets,
)

#: What a river zone reads on the observation route: a RIV package's RIV
#: observations, one per zone cell (ledger 344). The budget route reads the
#: river cells of the package ``package=`` names, or of every RIV package.
_RIV_FLOW = _ObservedQuantity(
    owner="RivFlowTargets",
    package_type="riv",
    obstype="riv",
    feature="cell",
    file_route="the model's cell-budget file",
    file_route_by_package=True,
    boundary="river boundary",
)


@dataclass
class RivFlowTargets(_CellZoneFlowTargets):
    """Observed river exchange (baseflow) targets aggregated from RIV budget cells.

    Sums the simulated RIV flow -- the exchange between the aquifer and a river
    boundary -- over a *group* of river cells to form one flow series per named
    zone, then aligns it with observed flow: a creek's gauged baseflow gain
    between two points, or a whole creek's. It is the flux target that tells
    hydraulic conductivity from recharge, which heads alone cannot. ``locations``
    maps each ``name`` to a set of cells (via geometry or an explicit cell list),
    exactly as :class:`~myflopy.modflow.mf6.observations.drn.DrnFlowTargets`'
    does -- the two share every zone rule (ledger 534). Pass an
    instance to ``cal.observe(...)`` / ``cal.forecast(...)``, or call
    ``.compare(model)`` / ``.stats(model)`` for standalone review.

    **A flow is in MF6's own sign, on both routes, and either sign is a target**
    (ledger 350's rule): MF6 writes water the aquifer gives to the river as
    NEGATIVE and water the river leaks into the aquifer as POSITIVE, so a GAINING
    zone's baseflow ``Q`` is entered as ``-Q`` -- the spelling
    :class:`~myflopy.modflow.mf6.observations.sfr.SfrFlowTargets` uses for a
    gauge, so a target set keeps its values when the creek moves to SFR -- and a
    net-LOSING zone is a positive number. Unlike ``DrnFlowTargets`` no sign is
    refused: a losing reach is physical (89 of 1,021 river cells on a real
    two-creek model, ledger 534).

    **A period's flow is its END-of-period flow** (ledger 336), the zone's cells
    added within the time step that ends the period; a period whose final step
    output control did not save is refused rather than read early. On the
    observation route the zone is its cells' ``RIV`` observations at the time
    step that ends each period, keyed by the period (ledger 346). A period in
    which the zone's cells hold no river boundary reads 0.0 on both routes
    (ledger 537); a zone whose cells hold no river boundary of the packages read
    in ANY period the budget file saved is REFUSED on the budget route, naming
    its layer(s) and ``package=`` -- a zone in the wrong layer or package, not
    one of zero flow (ledger 537, amended; a zone given a ``layer`` reads that
    layer alone). A zone given NO ``layer`` reads every layer of its columns,
    listed or drawn (ledgers 535 and 550).

    **Several RIV packages** (a main stem and a tributary, say): on the
    budget-file route ``package=None`` sums every RIV package's records in the
    zone's cells and ``package=<name>`` reads that package's alone; every package
    read must have been declared with ``save_flows=True``, or reading raises
    naming it (ledger 536). On the observation route a ``RIV``
    observation measures ONE package's boundary, so with several packages
    declaring the zones' labels ``package=`` is required (ledger 352).

    Parameters
    ----------
    locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, list or tuple
        The river zones, one name each; rows sharing a name are ONE zone, its
        cells united across the rows and their layers. Each row names its cells
        one of two ways. EXPLICIT: a ``cells`` entry holding a list of ints, one
        int or a comma-separated string (``"317, 342"``) -- ZERO-based cell
        numbers within the row's ``layer``, taken as given -- with no layer,
        COLUMNS read in every layer they are active in -- (an empty list, or
        ``""``, is a zone of no cells), except that reading the budget file
        refuses a ``(layer, cell)`` the model's grid does not have (ledger 343;
        MF6 refuses it on the observation route); the ``{name: cells}`` dict and a
        ``Series`` of cells indexed by name are shorthands for it. GEOMETRY: a
        row with no ``cells`` entry (no such column, or a missing value in it)
        and a geometry -- a polygon around a gauged reach, a buffered centerline
        -- which :meth:`zone_definitions` resolves to the cells it COVERS. One
        table may mix the two. Geometry comes from a ``GeoDataFrame``, from any
        other value (a vector-file path) read with ``geopandas.read_file``, or
        from a ``geometry`` column of shapely geometries in a dict, a list of
        record dicts or a ``DataFrame``, which becomes a ``GeoDataFrame`` in the
        GRID's coordinates (it carries no CRS; anything but a geometry there
        raises ``ValueError``). A dict with a ``name``, ``cells`` or
        ``geometry`` key is read as columns, and a list or tuple as record
        dicts. Kept: ``name_column``, ``cells``, ``group_column``,
        ``weight_column``, ``layer`` and the geometry; every other column is
        DROPPED. ``layer`` is zero-based: the layer a row's cells are in, which
        every route reads (MF6's cell-budget node is ``layer * ncpl + cell`` on
        a DISV, ledger 342). A BLANK layer (no column, or a missing value) is
        EVERY layer, on an EXPLICIT row and on a GEOMETRY row alike -- one row
        per layer, each holding the listed or covered columns ACTIVE in that
        layer, ledger 482's area rule (ledgers 535 and 550; an
        explicit row read layer 0). The ``{name: cells}`` dict and the
        ``Series`` carry no layer, so they read every layer; give record dicts
        with a ``layer`` to read one. Such a row needs the model to be read
        (``to_flopy_obs(None)`` refuses it), and a listed column the grid does
        not have, or one active in no layer, raises ``ValueError``; so does a
        missing ``name_column``, or a ``layer`` that is not a whole number 0 or
        more. A zone by RIV
        BOUNDNAME is not offered (ledger 539): draw the creek, or
        list its cells.
    values : str, pathlib.Path, DataFrame, Series, dict, list, tuple or None, optional
        Observed exchange flow through time, in model length cubed per model
        time unit, in MF6's own sign: NEGATIVE for water the aquifer gives to the
        river, so a gauged baseflow gain ``Q`` is entered as ``-Q``; positive for
        a net-losing zone. Either sign is accepted. :meth:`compare` subtracts it,
        as given, from the simulated flow, which carries MF6's sign on both
        routes. A path is read with ``pandas.read_csv``. LONG form has ``name``
        and ``value_column`` columns (headers matched ignoring case); WIDE form
        has one column per zone name beside ``time_column``; a table of only
        ``time_column`` and ``value_column``, or a ``Series`` indexed by time, is
        one series and needs exactly one zone (``ValueError`` otherwise). Names
        are NOT checked against ``locations`` and are joined to them exactly
        (case counts): a row naming no zone gets no cells and no simulated flow,
        and the PEST build refuses it. A blank value (missing, or empty text) is
        unmeasured; any other entry that is not a finite number raises
        ``ValueError`` naming the target and its column (ledger 340). A LONG or
        one-series table may carry a ``weight`` column (matched ignoring case):
        each READING's own PEST weight, which beats the zone's (the locations'
        ``weight`` column), which beats 1.0 (ledger 504); a blank means none of
        its own, a weight that is not a finite number from 0 raises, and a WIDE
        table with a ``weight`` column raises, naming the LONG form. Default
        ``None`` makes a capture template: observed, its rows register at weight
        0; forecast, every row is a forecast.
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
        prefix (``riv_flow``). A zone given as several rows carrying DIFFERENT
        labels is refused by the build -- one zone is one observation per
        period -- and rows with no label take the others', in the build and in
        :meth:`compare` alike (ledger 359).
    weight_column : str or None, optional
        Column of ``locations`` holding each zone's observation weight, renamed
        ``weight``. Default ``"weight"``; ``None``, or the default absent, leaves
        it ``NaN``. Any OTHER name the table lacks raises ``ValueError`` naming
        the columns it has, as does any other name given with ``{name: cells}``
        dict or ``Series`` locations, which carry no columns (ledger 609).
        A PEST build sets a zone's non-``NaN`` weight on each
        observation of that zone that has a target value and no weight of its
        own, and weight 1 where both are ``NaN``; a zone given as rows with
        DIFFERENT weights is refused, and every period without a target value
        gets weight 0 (ledger 340). PESTPP-IES draws a reading's noise with
        standard deviation 1/weight, so a flow in ft3/d at weight 1 carries
        1 ft3/d of noise. A blank is "none of its own"; anything but a blank or a finite
        number from 0 raises ``ValueError`` naming the zone (ledger
        609: it was not checked, and ``'x'`` reached the control
        file as 1.0).
    time_column : str, optional
        Column of ``values`` holding each row's time label, matched exactly.
        Default ``"time"``; ``ValueError`` when ``values`` lacks it after
        ``times`` is applied. :meth:`compare` and the PEST build read it as a
        zero-based stress period, by number (``1.0`` is period 1; a fraction
        raises), on both routes; a key naming no period the route holds is
        refused on both (ledgers 346, 538).
    value_column : str, optional
        Column of a LONG ``values`` table holding the observed flow, matched
        ignoring case. Default ``"flow"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels inserted as ``time_column`` when ``values`` is a ``dict`` or
        list/tuple without that column; a list/tuple of another length raises
        ``ValueError``. Ignored for a ``DataFrame``, ``Series`` or CSV path.
        Default ``None``.
    package : str or None, optional
        The NAME of the RIV package the zones read (FloPy's ``package_name``,
        any case, e.g. ``"riv"`` or ``"crisp"``); one the model does not have
        raises ``ValueError``. On the budget-file route it reads that package's
        river cells alone; on the observation route its declared ``RIV``
        observations alone, which is needed when the model carries several RIV
        packages and the zones' cells are declared in their observations -- a
        ``RIV`` observation measures ONE package's boundary in a cell (ledger
        352). Default ``None``: on the budget route every RIV package's records
        in each zone, summed; on the observation route the model's one RIV
        package, and with several, reading raises ``ValueError`` naming them
        when a zone cell's label is declared in any of them. Every package the
        budget route reads must have been declared with ``save_flows=True``, or
        reading raises naming it (ledger 536).
        :meth:`attach_flopy_obs` attaches to the named package.
    inequality : str or None, optional
        Whether every value of the set is a one-sided BOUND (ledger
        564), on the number AS WRITTEN -- MF6's sign, water the
        aquifer gives the river negative (ledger 350). ``"at_most"``: the
        simulated exchange must not exceed the value; ``"at_least"``: must not
        fall below it. So "the reach gains at least 2 of baseflow" is ``-2.0``
        with ``inequality="at_most"`` (simulated at most -2), and "it loses at
        most 3" is ``3.0`` with ``"at_most"``. Default ``None``: an ordinary
        target; anything else raises ``ValueError``. In a PEST build each measured
        observation is in PESTPP-IES's ``l_``/``g_`` group before the zone's own
        (``l_riv_flow``), which scores 0 wherever the bound holds;
        :meth:`compare` scores it so too. ``cal.forecast(...)`` refuses a bound
        set.

    Examples
    --------
    >>> baseflow = mf.RivFlowTargets(
    ...     locations=gauged_reaches_gdf,                  # one polygon per gauged reach
    ...     values={"time": [0], "rock_creek": [-208_000.0]},  # gain Q entered as -Q
    ... )
    >>> cal.observe(baseflow)

    See Also
    --------
    DrnFlowTargets, SfrFlowTargets, HeadTargets, LakeStageTargets
    """

    _QUANTITY: ClassVar[_ObservedQuantity] = _RIV_FLOW
    _SIGN: ClassVar[dict | None] = None
    _SUBJECT: ClassVar[str] = "RIV zone"
    _DEFAULT_CSV: ClassVar[str] = "riv_flow_targets.csv"

    def to_flopy_obs(
        self,
        model,
        *,
        csv_name: str = "riv_flow_targets.csv",
        kind: str = "RIV",
    ) -> dict[str, list[tuple[str, str, tuple[int, int]]]]:
        """Return a FloPy ``continuous`` dict with one MF6 RIV observation per member cell.

        Parameters
        ----------
        model : SimulationBase or None
            Passed to :meth:`zone_definitions`, whose ``Raises`` apply: it
            resolves geometry zones and zones given no ``layer``. Its grid also says the form each cellid is
            written in. ``None`` works when every zone gives explicit ``cells``
            AND a ``layer`` (a row with none reads the model's every layer)
            -- as when the records go into the ``observations=`` of the RIV being
            declared, before any model exists -- and writes DISV's ``(layer,
            cell)``, having no grid to read: a DIS or DISU model's build refuses
            that form, naming the package (ledger 490), so there declare the
            RIV without ``observations=`` and call :meth:`attach_flopy_obs` with
            the built model (``run.model()``) before running it.
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"riv_flow_targets.csv"``. Any name: :meth:`simulated_series`
            reads the file the RIV package's OBS declares, whatever it is called
            (ledger 344).
        kind : str, optional
            The MF6 observation type written on every record. Must be ``"RIV"``
            (any case), the default -- MF6's river-exchange observation, in MF6's
            sign (water to the river negative), the only type the observation
            route sums as a zone's flow -- and any other raises ``ValueError``
            (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, kind, cellid), ...]}``: ONE record per cell of
            each zone -- the ``(layer, cell)`` pairs its rows name, united, and
            every layer's for a row given no layer -- labelled
            ``<name>_c<cell>``, or ``<name>_lay<layer>_c<cell>`` for a zone
            reaching several layers, with every character that is not a letter
            or digit made ``_`` and leading ``_`` stripped; the cellid ZERO-based
            as FloPy takes it, in the model's own form -- ``(layer, row, col)``
            on DIS, ``(layer, cell)`` on DISV, ``(node,)`` on DISU (ledger 455)
            -- and ``(layer, cell)`` when ``model`` is ``None`` or has no
            discretization to read. A cell holding no river boundary gets a
            record too: MF6 accepts it and writes its no-value number, which the
            reader takes as 0.0. The shape ``ModflowUtlobs(continuous=)`` and a
            RIV package's ``observations=`` take.

        Raises
        ------
        ValueError
            When ``kind`` is not ``"RIV"``; as :meth:`zone_definitions`; or when
            two zones' cells would carry one label once written (``"a b"`` and
            ``"a_b"``), naming them.
        """

        return self._to_flopy_obs(model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        model,
        *,
        pname: str = "riv_flow_obs",
        filename: str = "riv_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "RIV",
        package=None,
    ):
        """Attach an MF6 utility observation package to a RIV package for these zone targets.

        The package is added to the FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        model : SimulationBase or None
            The model whose FloPy RIV package receives the OBS package unless
            ``package`` is given: the one the targets name (their ``package``
            field, ``ValueError`` when the model has no RIV package of that
            name), else the model's only RIV package, found by type (several
            raise ``ValueError`` naming them and ``package=``, none
            ``AttributeError``; ledger 352); also passed to :meth:`to_flopy_obs`
            to resolve the zones. With ``package`` given and every zone's
            ``cells`` and ``layer`` explicit, ``None`` works.
        pname : str, optional
            FloPy package name. Default ``"riv_flow_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory. Default ``"riv_flow_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"riv_flow_targets.csv"`` for the
            default ``filename``).
        kind : str, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`: ``"RIV"``, the
            default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to instead of
            the one ``model`` resolves (see ``model``). Default ``None``. Reading
            the observations back needs the targets' own ``package`` field to
            name it when the model carries several RIV packages.

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

        return self._attach_flopy_obs(
            model, pname=pname, filename=filename, csv_name=csv_name, kind=kind, package=package,
        )


class BoundRivFlowTargets(_BoundCellZoneFlowTargets):
    """Model-bound RIV exchange-flow helper returned from ``model.targets.riv_flow``.

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from, and ``get``,
    ``to_long`` and ``summary`` resolve geometry zones against it. Flows are in
    MF6's sign -- negative where the aquifer feeds the river (ledger 350).
    """

    _PLOT_TITLE: ClassVar[str] = "Observed vs simulated RIV exchange"

    def __init__(self, model, targets: RivFlowTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        super().__init__(model, targets)

    def to_flopy_obs(self, *, csv_name: str = "riv_flow_targets.csv", kind: str = "RIV"):
        """A FloPy ``continuous`` observation dict for the RIV zone targets.

        Zones drawn as geometry are resolved against the BOUND model
        (:meth:`RivFlowTargets.zone_definitions`); one record per zone cell.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"riv_flow_targets.csv"``. Any name: the observation route
            reads whichever file the RIV package's OBS declares (ledger 344).
        kind : str, optional
            The MF6 observation type written on every record. Must be ``"RIV"``
            (any case), the default, and any other raises ``ValueError`` (ledger
            344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "RIV", cellid), ...]}``, one record per cell of
            each zone, the cellid ZERO-based in the bound model's own form, as
            :meth:`RivFlowTargets.to_flopy_obs` returns it (whose ``Raises``
            apply).
        """

        return self.targets.to_flopy_obs(self.model, csv_name=csv_name, kind=kind)

    def attach_flopy_obs(
        self,
        *,
        pname: str = "riv_flow_obs",
        filename: str = "riv_flow_targets.obs",
        csv_name: str | None = None,
        kind: str = "RIV",
        package=None,
    ):
        """Attach an MF6 RIV observation package to the bound model for these zones.

        Zones are resolved against the bound model (:meth:`to_flopy_obs`) and
        the package is added to its FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        pname : str, optional
            FloPy package name of the OBS package. Default ``"riv_flow_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory; its stem names the default ``csv_name``. Default
            ``"riv_flow_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"riv_flow_targets.csv"`` for the
            default ``filename``).
        kind : str, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`: ``"RIV"``, the
            default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to. Default
            ``None`` takes the bound model's RIV package: the one the targets
            name (their ``package`` field; ``ValueError`` when the model has no
            RIV package of that name), else its only RIV package, found by type
            (several raise ``ValueError`` naming them, none ``AttributeError``;
            ledger 352). Reading the observations back needs the targets' own
            ``package`` field to name it when the model carries several.

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
