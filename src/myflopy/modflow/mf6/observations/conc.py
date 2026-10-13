"""Concentration observation targets for GWT calibration (``ConcTargets``)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from myflopy.modflow.mf6.observations.heads import HeadTargets


@dataclass
class ConcTargets(HeadTargets):
    """Measured vs simulated CONCENTRATION at monitoring points, for PEST.

    The transport twin of :class:`~myflopy.modflow.mf6.observations.heads.\
HeadTargets`, and deliberately a *configuration* of it rather than a copy:
    concentration is sampled exactly the way head is -- a point value at a
    ``(layer, cell)`` per stress period -- so the only differences are which
    model table is read and what the field is called.

    That is why this subclasses ``HeadTargets`` instead of the named-series
    targets (``DrnFlowTargets`` and friends): those aggregate a FLUX over a zone
    of cells, which is the wrong semantics for a state variable, and they carry
    no layer in their identity while concentration is inherently three
    dimensional.

    Concentration is the observation that makes transport calibration worth
    doing: heads are insensitive to porosity and only weakly constrain flow
    paths, whereas concentration responds strongly to both. Measured on the
    canonical model, a 2.5x porosity change moves concentration by 54% and a
    2.5x K change by 38%.

    Parameters
    ----------
    locations : str, pathlib.Path, GeoDataFrame, DataFrame, Series, dict, or list or tuple of dict
        Where each sampling point is, read exactly as
        :class:`~myflopy.modflow.mf6.observations.heads.HeadTargets` reads it:
        a vector file path or ``GeoDataFrame`` of POINTS (matched to the one
        cell holding each, active in its layer), or a table, ``Series`` or
        ``{name: cell}`` dict of zero-based cell indices within a layer
        (``cell`` or ``node`` column). ``ValueError`` when ``name_column`` is
        missing, or a table with no geometry has no ``cell``/``node`` column.
    values : str, pathlib.Path, pandas.DataFrame, pandas.Series, dict, list or tuple
        Measured concentrations through time: a CSV path or a table in LONG
        form (``name`` and ``value_column`` columns), WIDE form (one column per
        location) or as one series when there is a single location. Every name
        must match a location (stripped, ignoring case) or ``ValueError``.
        Units are the transport model's concentration units. A blank value (missing, or empty text) is unmeasured; any other
        entry that is not a finite number raises ``ValueError`` naming the target and
        its column (ledger 340).
    name_column : str, optional
        Column of ``locations`` holding each location's name. Default
        ``"name"``; required.
    layer_column : str or None, optional
        Column of ``locations`` holding each location's zero-based model layer.
        Default ``"layer"``; ``None`` or an absent column means layer 0.
    group_column : str or None, optional
        Column of ``locations`` holding each location's observation-group label,
        carried as ``group``. Default ``"group"``; ``None`` or absent leaves
        ``<NA>``. It is the PEST observation group of the location's MEASURED
        observations (ledger 358), made a PEST name and refused as
        :class:`~myflopy.modflow.mf6.observations.heads.HeadTargets`'
        ``group_column`` says; a missing or blank label puts them in the set's
        prefix (``conc``).
    weight_column : str or None, optional
        Column of ``locations`` holding a per-location PEST weight, carried as
        ``weight``; a value that is not ``NaN`` becomes the weight of each of
        that location's observations that has a target value. Default
        ``"weight"``; ``None`` or absent leaves ``NaN``.
    time_column : str, optional
        Column of ``values`` holding each row's time label, read as a zero-based
        stress period, by number, by :meth:`compare` and the PEST build: ``1.0``
        is period 1, and a fraction, or two values for one location in one
        period, raises ``ValueError``. Default ``"time"``; required once
        ``times`` is applied.
    value_column : str, optional
        Column of a LONG ``values`` table holding the concentration, matched
        ignoring case. Default ``"conc"`` -- the one default that differs from
        ``HeadTargets`` (``"head"``). ``from_cells``/``from_records`` still
        default it to ``"head"``.
    times : list, tuple, pandas.Index, pandas.Series or None, optional
        Time labels for ``values`` given as a ``dict`` or list/tuple with no
        ``time_column`` column; a length mismatch raises ``ValueError``.
        Ignored for a ``DataFrame``, ``Series`` or CSV. Default ``None``.
    screen_top_column : str or None, optional
        Column of ``locations`` that would hold a well's screen-top elevation.
        Default ``"screen_top"``. A concentration target is not read over a
        screen: a flux-weighted concentration needs the flow model's heads and K
        beside the transport file, which is deferred (ledger 556),
        so a locations table carrying this column and ``screen_bottom_column``
        raises ``ValueError`` -- give each sampling point its layer, or pass
        ``None`` for both to ignore the columns.
    screen_bottom_column : str or None, optional
        The screen-bottom twin of ``screen_top_column``, refused alike. Default
        ``"screen_bottom"``.
    inequality : str or None, optional
        ``"at_most"`` (the simulated concentration must not exceed each value -- a
        sample below detection: the detection limit) or ``"at_least"`` (must not
        fall below it) makes every value of the set a one-sided BOUND, read as
        :class:`~myflopy.modflow.mf6.observations.heads.HeadTargets`'
        ``inequality`` says (ledger 564): PESTPP-IES's ``l_``/``g_``
        group in a PEST build, PEST++'s residual in :meth:`compare`. Default
        ``None``: ordinary targets; anything else raises ``ValueError``.

    Examples
    --------
    ::

        wells = ConcTargets.from_cells(cells=[120, 344], values=measured)
        cal.observe(wells)

    Notes
    -----
    The normalized frames keep ``head_target``/``sim_head`` as their internal
    value columns for every field kind. They are frame-internal names meaning
    "the target value" and "the simulated value"; sharing them is what lets the
    calibration-plot and PEST observation machinery serve concentration without
    change (ledger 102).
    """

    _TABLE_ATTRIBUTE: ClassVar[str] = "all_conc"
    _STORE_COLUMN: ClassVar[str] = "conc"
    _FIELD_LABEL: ClassVar[str] = "concentration"
    _OBS_KIND: ClassVar[str] = "CONCENTRATION"
    #: Screens are refused (ledger 556).
    _SCREENS: ClassVar[bool] = False

    #: the long-form input column name callers supply measured values under
    value_column: str = "conc"

    def simulated_conc(self, model):
        """Return simulated concentration by stress period, wide by name.

        The concentration spelling of :meth:`HeadTargets.simulated_heads`; both
        read the same table, and either may be called.

        Parameters
        ----------
        model : SimulationBase
            A run GWT model holding the concentration file, as
            ``run.model("<gwt name>")`` returns; a GWF model raises
            ``AttributeError`` (``all_conc`` is gated to GWT models). Each
            period's value is the one at its END, and a period whose final step
            was not saved raises ``ValueError`` naming ``SAVE ... LAST``.

        Returns
        -------
        pandas.DataFrame
            ``per`` plus one column per observation label.
        """

        return self.simulated_heads(model)


__all__ = ["ConcTargets"]
