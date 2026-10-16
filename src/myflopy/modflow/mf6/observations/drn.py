"""DRN flow observation targets (`DrnFlowTargets` + bound view).

The zone machinery is shared with ``RivFlowTargets`` in
:mod:`myflopy.modflow.mf6.observations._zones` (ledger 534); this
module names what is DRN's own -- MF6's ``DRN`` observation, a drain's one-signed
flow, the defaults.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from myflopy.modflow.mf6.observations._shared import _ObservedQuantity
from myflopy.modflow.mf6.observations._zones import (
    _BoundCellZoneFlowTargets,
    _CellZoneFlowTargets,
)

#: What a seepage zone reads on the observation route: a DRN package's DRN
#: observations, one per zone cell (ledger 344). The budget route reads the
#: drains of the package ``package=`` names, or of every DRN package
#: (ledger 534).
_DRN_FLOW = _ObservedQuantity(
    owner="DrnFlowTargets",
    package_type="drn",
    obstype="drn",
    feature="cell",
    file_route="the model's cell-budget file",
    file_route_by_package=True,
    boundary="drain",
)


@dataclass
class DrnFlowTargets(_CellZoneFlowTargets):
    """Observed seepage/discharge targets aggregated from DRN budget cells.

    Sums the simulated DRN flow over a *group* of drain cells to form one flow
    series per named seepage zone, then aligns it with observed discharge -- the
    pattern for calibrating spring/seep or drain-discharge measurements where the
    field datum is a lumped flux over many cells. Unlike the SFR/LAK targets,
    ``locations`` maps each ``name`` to a set of drain cells (via geometry or an
    explicit cell list), and the comparison aggregates the DRN cell-by-cell
    budget. Pass an instance to ``cal.observe(...)`` / ``cal.forecast(...)``, or
    call ``.compare(model)`` / ``.stats(model)`` for standalone review. Its river
    twin is :class:`~myflopy.modflow.mf6.observations.riv.RivFlowTargets`; the two
    share every zone rule (ledger 534).

    **A period's flow is its END-of-period flow** (ledger 336): the zone's cells
    ADD within the time step that ends the period -- that is what a lumped seepage
    measurement is -- and the steps a period saved do not; a period whose final
    step output control did not save is refused rather than read early. On the
    observation route the zone is its cells' observations at the time step that
    ends each period, keyed by the period (ledger 346). A period in which the
    zone's cells hold no drain reads 0.0 on both routes (ledger 537); a zone whose
    cells hold no drain of the packages read in ANY period the budget file saved
    is REFUSED on the budget route, naming its layer(s) and ``package=`` -- a zone
    in the wrong layer or package, not one of zero flow (ledger 537, amended; a
    zone given a ``layer`` reads that layer alone).

    **A flow is in MF6's own sign, on both routes**: MF6 writes water a drain
    takes out of the aquifer as NEGATIVE, so ``sim_flow`` is zero or negative and
    a measured spring or drain discharge ``Q`` is entered as ``-Q``. A target
    value above zero can never be matched -- a drain cannot add water -- and is
    REFUSED when the targets are built, naming ``-Q`` (ledger 350); zero (a dry
    spring) is accepted.

    **A zone's cell is a ``(layer, cell)`` pair, on every route** (ledger 342): the
    budget file names a drain by its model-wide node. A zone is summed ONCE over
    the union of its rows' cells, and written to MF6 as one observation per cell.
    A zone given no ``layer`` reads every layer of its columns, whether they are
    DRAWN as geometry (ledger 535) or LISTED as ``cells`` (ledger
    550); both read layer 0 alone, so a polygon over drains stacked
    in four layers read 1.8% of their flow, and the same drains listed by column
    read 245 of their 935 records.

    **The observation route reads what the DRN package DECLARES** (ledger 344):
    a zone reads MF6's ``DRN`` observations only when every one of its cells'
    labels (as :meth:`to_flopy_obs` writes them) is declared in the DRN
    package's own OBS as a ``DRN`` observation OF THAT CELL, from the file that
    OBS names; a label declared as another quantity only, for another cell or
    twice raises, and so does a set of which some zone cells are declared and
    some are not.

    Parameters
    ----------
    locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, list or tuple
        The seepage zones, one name each; rows sharing a name are ONE zone, its
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
        and a geometry -- :meth:`zone_definitions` resolves it to the cells it
        COVERS. One table may mix the two. Geometry comes from a
        ``GeoDataFrame``, from any other value (a vector-file path) read with
        ``geopandas.read_file``, or from a ``geometry`` column of shapely
        geometries in a dict, a list of record dicts or a ``DataFrame``, which
        becomes a ``GeoDataFrame`` in the GRID's coordinates (it carries no CRS;
        anything but a geometry there raises ``ValueError``). A dict with a
        ``name``, ``cells`` or ``geometry`` key is read as columns, and a list or
        tuple as record dicts. Kept: ``name_column``, ``cells``,
        ``group_column``, ``weight_column``, ``layer`` and the geometry; every
        other column is DROPPED, a ``node`` or ``cell`` column included.
        ``layer`` is zero-based: the layer a row's cells are in, which every
        route reads -- MF6's cell-budget node is ``layer * ncpl + cell`` on a
        DISV, and the budget route matches the pair (ledger 342). A BLANK layer
        (no column, or a missing value) is EVERY layer, on an EXPLICIT row and
        on a GEOMETRY row alike -- one row per layer, each holding the listed or
        covered columns ACTIVE in that layer, ledger 482's area rule (ledgers
        535 and 550; an explicit row read layer 0). The
        ``{name: cells}`` dict and the ``Series`` carry no layer, so they read every
        layer; give record dicts with a ``layer`` to read one. Such a row needs
        the model to be read (``to_flopy_obs(None)`` refuses it), and a listed
        column the grid does not have, or one active in no layer, raises
        ``ValueError``; so does a missing ``name_column``, or a ``layer`` that
        is not a whole number 0 or more.
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
        its column (ledger 340). A LONG or one-series table may carry a
        ``weight`` column (matched ignoring case): each READING's own PEST
        weight, which beats the zone's (the locations' ``weight`` column),
        which beats 1.0 (ledger 504); a blank means none of its own, a
        weight that is not a finite number from 0 raises, and a WIDE table
        with a ``weight`` column raises, naming the LONG form. Default ``None`` makes a
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
        ``weight``. Default ``"weight"``; ``None``, or the default absent, leaves
        it ``NaN``. Any OTHER name the table lacks raises ``ValueError`` naming
        the columns it has, as does any other name given with ``{name: cells}``
        dict or ``Series`` locations, which carry no columns (ledger 609).
        A PEST build sets a zone's non-``NaN`` weight on each
        observation of that zone that has a target value and no weight of its
        own, and weight 1 where both are ``NaN``; a zone given as rows with
        DIFFERENT weights is refused, and every period without a target value
        gets weight 0 (ledger 340). A blank is "none of its own"; anything but a
        blank or a finite number from 0 raises ``ValueError`` naming the zone
        (ledger 609: it was not checked, and ``'x'`` reached the
        control file as 1.0).
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
        The NAME of the DRN package the zones read (FloPy's ``package_name``,
        any case); one the model does not have raises ``ValueError``. On the
        budget-file route it reads that package's drains alone (ledger
        534; it was refused there); on the observation route its
        declared ``DRN`` observations alone, which is needed when the model
        carries several DRN packages and the zones' cells are declared in their
        observations -- a DRN observation measures ONE package's drain in a cell
        (ledger 352). Default ``None``: on the budget route every DRN package's
        drains in each zone, summed; on the observation route the model's one
        DRN package, and with several, reading raises ``ValueError`` naming them
        when a zone cell's label is declared in any of them. Every package the
        budget route reads must have been declared with ``save_flows=True``, or
        reading raises naming it (ledger 536).
        :meth:`attach_flopy_obs` attaches to the named package.
    inequality : str or None, optional
        Whether every value of the set is a one-sided BOUND (ledger
        564), on the number AS WRITTEN -- MF6's sign, a drain's
        discharge zero or negative (ledger 350). ``"at_most"``: the simulated flow
        must not exceed the value; ``"at_least"``: must not fall below it. So "a
        spring discharges at least 2" is ``-2.0`` with ``inequality="at_most"``
        (simulated at most -2), and "seepage of at most 5" is ``-5.0`` with
        ``"at_least"``. Default ``None``: an ordinary target; anything else raises
        ``ValueError``. A value above zero is refused either way. In a PEST build
        each measured observation is in PESTPP-IES's ``l_``/``g_`` group before the
        zone's own (``l_seeps``, ``g_drn_flow``), which scores 0 wherever the bound
        holds; :meth:`compare` scores it so too. ``cal.forecast(...)`` refuses a
        bound set.

    Examples
    --------
    >>> DrnFlowTargets(locations=seep_zones_gdf,
    ...                values=observed_seepage_df)       # -> cal.observe(...)

    See Also
    --------
    RivFlowTargets, SfrFlowTargets, HeadTargets, LakeStageTargets
    """

    _QUANTITY: ClassVar[_ObservedQuantity] = _DRN_FLOW
    _SIGN: ClassVar[dict | None] = {
        "quantity": "MF6's DRN flow",
        "leaving": "the aquifer through a drain",
        "measured": "a measured discharge",
    }
    _SUBJECT: ClassVar[str] = "DRN zone"
    _DEFAULT_CSV: ClassVar[str] = "drn_flow_targets.csv"

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
            resolves geometry zones and zones given no ``layer``. Its grid also says the form each cellid is
            written in. ``None`` works when every zone gives explicit ``cells``
            AND a ``layer`` (a row with none reads the model's every layer)
            -- as when the records go into the ``observations=`` of the DRN being
            declared, before any model exists -- and writes DISV's ``(layer,
            cell)``, having no grid to read: a DIS or DISU model's build refuses
            that form, naming the package (ledger 490), so there declare the
            DRN without ``observations=`` and call :meth:`attach_flopy_obs` with
            the built model (``run.model()``) before running it.
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
            ``{csv_name: [(label, kind, cellid), ...]}``: ONE record per
            cell of each zone -- the ``(layer, cell)`` pairs its rows name,
            united, whether the zone is given as one row or several, and every layer's for a row given no layer -- labelled
            ``<name>_c<cell>``, or ``<name>_lay<layer>_c<cell>`` for a zone
            reaching several layers, with every character that is not a letter
            or digit made ``_`` and leading ``_`` stripped; the cellid ZERO-based
            as FloPy takes it, in the model's own form -- ``(layer, row, col)``
            on DIS, ``(layer, cell)`` on DISV, ``(node,)`` on DISU
            (``CellNumbering.cellid``, ledger 455) -- and ``(layer, cell)``
            when ``model`` is ``None`` or has no discretization to read. The
            shape ``ModflowUtlobs(continuous=)`` and a DRN package's
            ``observations=`` take.

        Raises
        ------
        ValueError
            When ``kind`` is not ``"DRN"``; as :meth:`zone_definitions`; or when
            two zones' cells would carry one label once written (``"a b"`` and
            ``"a_b"``), naming them.
        """

        return self._to_flopy_obs(model, csv_name=csv_name, kind=kind)

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
            resolve the zones. With ``package`` given and every zone's
            ``cells`` and ``layer`` explicit, ``None`` works.
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

        return self._attach_flopy_obs(
            model, pname=pname, filename=filename, csv_name=csv_name, kind=kind, package=package,
        )


class BoundDrnFlowTargets(_BoundCellZoneFlowTargets):
    """Model-bound DRN seepage-zone helper returned from ``model.targets.drn_flow``.

    Every method reads the BOUND model unless given another: ``model=None`` (the
    default) means the model the registry was read from, and ``get``,
    ``to_long`` and ``summary`` resolve geometry zones against it. Flows are in
    MF6's sign -- zero or negative, water the drains take out (ledger 350).
    """

    _PLOT_TITLE: ClassVar[str] = "Observed vs simulated DRN seepage"

    def __init__(self, model, targets: DrnFlowTargets):
        """Bind ``targets`` to ``model`` so its methods default to that model."""

        super().__init__(model, targets)

    def to_flopy_obs(self, *, csv_name: str = "drn_flow_targets.csv", kind: str = "DRN"):
        """A FloPy ``continuous`` observation dict for the DRN seepage-zone targets.

        Zones drawn as geometry are resolved against the BOUND model
        (:meth:`DrnFlowTargets.zone_definitions`); one record per zone cell.

        Parameters
        ----------
        csv_name : str, optional
            The CSV file MF6 writes these observations to -- the dict's one key.
            Default ``"drn_flow_targets.csv"``. Any name: the observation route
            reads whichever file the DRN package's OBS declares (ledger 344).
        kind : str, optional
            The MF6 observation type written on every record. Must be ``"DRN"``
            (any case), the default -- MF6's drain-flow observation, outflow
            negative -- and any other raises ``ValueError`` (ledger 344).

        Returns
        -------
        dict of str to list of tuple
            ``{csv_name: [(label, "DRN", cellid), ...]}``, one record per
            cell of each zone, the cellid ZERO-based in the bound model's own
            form, as :meth:`DrnFlowTargets.to_flopy_obs` returns it (whose
            ``Raises`` apply).
        """

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
        """Attach an MF6 DRN observation package to the bound model for these zones.

        Zones are resolved against the bound model (:meth:`to_flopy_obs`) and
        the package is added to its FloPy objects in memory; writing and running
        the simulation afterwards is the caller's.

        Parameters
        ----------
        pname : str, optional
            FloPy package name of the OBS package. Default ``"drn_flow_obs"``.
        filename : str or pathlib.Path, optional
            The OBS input file name FloPy writes, in the owning model's
            directory; its stem names the default ``csv_name``. Default
            ``"drn_flow_targets.obs"``.
        csv_name : str or None, optional
            The CSV file MF6 writes the observations to. Default ``None`` takes
            ``filename``'s stem plus ``.csv`` (``"drn_flow_targets.csv"`` for
            the default ``filename``).
        kind : str, optional
            The MF6 observation type, as in :meth:`to_flopy_obs`: ``"DRN"``, the
            default, and nothing else (``ValueError``).
        package : flopy object or None, optional
            The FloPy package (or model) to attach the OBS package to. Default
            ``None`` takes the bound model's DRN package: the one the targets
            name (their ``package`` field; ``ValueError`` when the model has no
            DRN package of that name), else its only DRN package, found by type
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
