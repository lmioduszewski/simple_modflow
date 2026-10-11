"""LAK, SFR, and combined surface-water package explorers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from flopy.mf6.mfbase import FlopyException, MFDataException

from myflopy import viz as figs
from myflopy._deprecation import deprecated_instance_getattr
from myflopy._logging import get_logger
from myflopy.viz import mpl_axes

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.mf6.package_budget import (
    build_lak_budget_result_table,
    build_lak_budget_term_table,
    build_lak_stage_change_table,
    build_lak_stage_result_table,
    build_sfr_budget_result_table,
    build_sfr_budget_term_table,
    build_sfr_long_profile_table,
    build_sfr_stage_result_table,
    build_surface_water_exchange_cell_table,
)
from myflopy.modflow.mf6.package_explorer_utils import (
    _filter_normalized_table,
    _normalize_connection_type_filter,
)
from myflopy.modflow.mf6.package_plotting import (
    FieldMappable,
    SpatialView,
    _apply_backend,
    _blue_white_red_diverging_colorscale,
    _exchange_colorscale,
    _symmetric_color_limit,
    build_cell_input_map_payload,
    build_lak_q_map_payload,
    build_sfr_q_map_payload,
    build_surface_water_q_map_payload,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_registry import (
    checked_package_name,
    get_default_budget_term,
    get_package_result_spec,
    one_package,
    package_outputs,
)
from myflopy.modflow.mf6.package_results import (
    CellBudgetResultsExplorer,
    PackageBudgetTermExplorer,
    StageResultsExplorer,
)
from myflopy.modflow.mf6.package_tables import (
    build_lak_connection_table,
    build_lak_input_table,
    build_sfr_input_table,
    summarize_input_table,
)
from myflopy.modflow.mf6.period_end import (
    at_last_saved_step,
    at_last_saved_steps,
    last_saved_steps,
)
from myflopy.modflow.utils.datatypes.hover import (
    cell_input_hover,
    lak_hover,
    sfr_hover,
    surface_water_hover,
)

logger = get_logger(__name__)


def _model_units(model) -> dict[str, str | None]:
    """``{"length": ..., "time": ...}`` -- the model's own abbreviations, for a hover."""

    length, time = unit_labels(model)
    return {"length": length, "time": time}

#: Everything that can come back from "ask the model for a stage result".
#:
#: The stage lookups below wrap a whole pipeline -- resolve the package, open
#: its binary output, normalize the frame -- not a single call, so the tuple is
#: correspondingly wide. Each member was measured against flopy 3.10:
#: ``AttributeError`` when the package is not attached (or its ``output.stage()``
#: is None because it declares no stage fileout), ``EOFError`` on a TRUNCATED
#: stage file, ``ValueError`` on an empty or all-zero one, ``OSError`` on an
#: unreadable path, and flopy's own two exception classes -- which subclass
#: ``Exception`` directly, so no builtin covers them -- from the lazy simulation
#: load a file-backed model performs on first access.
_STAGE_UNAVAILABLE = (
    AttributeError, TypeError, ValueError, EOFError, OSError,
    MFDataException, FlopyException,
)


def _budget_file(model, package_type: str, package: str | None):
    """One LAK or SFR package's budget-file reader, refused by the VIEW's spelling.

    ``package`` is the package's NAME, ``None`` the model's only one; several and
    none named raise ``SeveralPackagesError`` naming ``model.<type>(package=...)``
    (ledger 352), before ``model.outputs`` could name its own spelling.
    """

    one_package(model, package_type, package)
    return package_outputs(model, package_type, package).bud


def _join_feature_stage(frame, stage_table, *, per=None):
    """Join per-cell feature stage onto an exchange frame for the hover.

    ``stage_table`` is a normalized stage result table (one row per feature-cell
    per saved step, labelled by the stage file's own ``kstpkper``). With ``per``,
    the stage joined is that period's LAST saved step on the STAGE file's clock
    (ledger 337's rule) -- the exchange beside it is drawn at its own budget file's
    last saved step, each file read on its own clock (ledger 301). It averaged every
    step of the period until ledger 348, and on LAK, whose ``per`` was a record
    index, joined a record of an EARLIER period. A stage file that saved nothing in
    ``per`` leaves the frame without a stage, as a missing stage file does. Cells
    touched by more than one feature take the mean stage, matching
    :class:`StageResultsExplorer`'s series aggregate.

    Parameters
    ----------
    frame : pandas.DataFrame
        The exchange rows being drawn, with a ``cell`` column.
    stage_table : pandas.DataFrame
        ``build_lak_stage_result_table`` / ``build_sfr_stage_result_table``.
    per : int or None, optional
        The stress period drawn, zero-based. ``None`` joins the mean over every
        record.

    Returns
    -------
    pandas.DataFrame
        ``frame`` with a ``stage`` column, or ``frame`` unchanged when it is empty,
        already has one, or there is no stage to join.
    """

    if frame.empty or "stage" in frame.columns or stage_table.empty:
        return frame
    selected = stage_table
    if per is not None and "kstpkper" in selected.columns:
        if int(per) not in last_saved_steps(selected["kstpkper"]):
            logger.debug("the stage file saved nothing in period %s; the hover omits it", per)
            return frame
        selected = at_last_saved_step(selected, per)
    elif per is not None and "per" in selected.columns:
        selected = selected[selected["per"] == int(per)]
    if selected.empty or "stage" not in selected.columns:
        return frame
    mapping = (
        pd.to_numeric(selected["stage"], errors="coerce")
        .groupby(selected["cell"])
        .mean()
    )
    joined = frame.copy()
    joined["stage"] = joined["cell"].map(mapping)
    return joined


def join_sfr_stage(model, frame, *, per=None, package=None):
    """Add SFR reach stage to an exchange frame; a no-op if stage is unavailable.

    ``package`` names the SFR package whose stage file is read (``None``: the
    model's only one); several packages and none named still raise
    ``SeveralPackagesError`` -- not a missing stage, and not swallowed (ledger 352).
    """

    try:
        stage_table = build_sfr_stage_result_table(model, package=package)
    except _STAGE_UNAVAILABLE:
        logger.debug(
            "no SFR stage for %s; the exchange hover will omit it",
            getattr(model, "name", model), exc_info=True,
        )
        return frame
    return _join_feature_stage(frame, stage_table, per=per)


def join_lak_stage(model, frame, *, per=None, package=None):
    """Add lake stage to an exchange frame; a no-op if stage is unavailable.

    ``package`` names the LAK package whose stage file is read (``None``: the
    model's only one); several packages and none named still raise
    ``SeveralPackagesError`` -- not a missing stage, and not swallowed (ledger 352).
    """

    try:
        stage_table = build_lak_stage_result_table(model, package=package)
    except _STAGE_UNAVAILABLE:
        logger.debug(
            "no lake stage for %s; the exchange hover will omit it",
            getattr(model, "name", model), exc_info=True,
        )
        return frame
    return _join_feature_stage(frame, stage_table, per=per)


class SfrReachProfileView:
    """A single-field, reach-ordered SFR profile: one table and its line figure.

    The field-level analogue of :class:`SfrProfileView` -- reached from a field
    explorer (``sfr.results.q.profile`` / ``sfr.results.stage.profile``), it
    plots one field (exchange or stage) along the reach / cumulative-distance
    axis. Follows the house view shape (``docs/view_layer_conventions.md``):
    ``get`` for the frame, ``plot`` for the figure, ``summary`` for the digest;
    calling the view rebinds the stress period, so these are the same figure::

        model.sfr.results.q.profile.plot(per=3)
        model.sfr.results.q.profile(per=3).plot()
    """

    def __init__(self, explorer, *, y_column: str, y_label: str, label: str, per: int = 0):
        """Bind the profile to ``explorer`` (a field explorer), field, and period."""

        self.explorer = explorer
        self.y_column = str(y_column)
        self.y_label = str(y_label)
        self.label = str(label)
        self.per = int(per)

    def __repr__(self) -> str:
        """Show the plotted field and bound period (the view's only state)."""

        return f"{type(self).__name__}(field={self.y_label!r}, per={self.per})"

    def __call__(self, *, per: int) -> SfrReachProfileView:
        """Return an equivalent view bound to stress period ``per``.

        Parameters
        ----------
        per : int
            Zero-based stress period the new view reads when its verbs are called
            without one; required, keyword-only, passed through ``int()`` (so
            ``None`` raises ``TypeError``). This view is left unchanged.

        Returns
        -------
        SfrReachProfileView
            A view of the same field bound to ``per``.
        """

        return type(self)(
            self.explorer, y_column=self.y_column, y_label=self.y_label, label=self.label, per=per
        )

    def _period(self, per: int | None) -> int:
        """Resolve an explicit ``per`` against the period bound to this view."""

        return self.per if per is None else int(per)

    def get(self, *, per: int | None = None) -> pd.DataFrame:
        """Return the reach-ordered profile table for the selected period.

        One row per reach-cell, read at the period's LAST saved step on the field's
        own file's clock -- the stage file for ``stage``, the budget file for the
        exchange (``period_end.at_last_saved_step``, ledger 337's rule), as
        ``sfr.results.profile`` beside it reads both. It held EVERY step the period
        saved until the review of ledger 348, so at ``nstp=3`` under ``SAVE ... ALL``
        each reach carried three values and :meth:`plot` zig-zagged through them.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period. ``None`` (default) is the period bound to the
            view (``profile(per=N)``, 0 unless rebound).

        Returns
        -------
        pandas.DataFrame
            The field explorer's rows at that step, sorted by reach and cell; empty
            when the period has none.
        """

        period = self._period(per)
        frame = self.explorer.get(per=period)
        if frame.empty:
            return frame
        if "kstpkper" in frame.columns:
            frame = at_last_saved_step(frame, period)
        return frame.sort_values(["reach", "cell"]).reset_index(drop=True)

    def summary(self, *, per: int | None = None) -> pd.DataFrame:
        """Return a compact digest of the profiled field.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period, read at its last saved step as :meth:`get`
            reads it. ``None`` (default) is the period bound to the view
            (``profile(per=N)``, 0 unless rebound).

        Returns
        -------
        pandas.DataFrame
            One row: the view's label and the record, period, layer and cell
            counts of :meth:`get`'s table, naming the profiled column.
        """

        frame = self.get(per=per)
        return summarize_input_table(
            frame,
            label=self.label,
            value_columns=[column for column in (self.y_column,) if column in frame.columns],
        )

    def plot(
        self,
        *,
        per: int | None = None,
        x: str = "distance",
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot the field along the reach axis, by cumulative distance or reach number.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period, read at its last saved step as :meth:`get`
            reads it. ``None`` (default) is the period bound to the view
            (``profile(per=N)``, 0 unless rebound). A period with no rows draws an
            empty figure rather than raising.
        x : str, default "distance"
            The horizontal axis. ``"distance"`` (any case) is ``distance_mid``, the
            reach mid-point's running total of ``rlen`` in reach-number order
            (model length units); ANY other value is the zero-based reach number --
            ``"reach"`` is the spelling, but a misspelling draws by reach too,
            without a word.
        plot_fig : bool, default False
            Call ``show()`` on the figure before returning.
        return_fig : bool, default True
            Return the figure. With ``False`` the method returns ``None``, so
            ``plot_fig=False, return_fig=False`` builds a figure nobody sees.

        Returns
        -------
        myflopy.viz.Fig or None
            One ``lines+markers`` trace of the field (``q_gwf`` or ``stage``) along
            the chosen axis, or ``None`` when ``return_fig=False``.

        Raises
        ------
        KeyError
            When the table lacks the chosen axis column.
        """

        period = self._period(per)
        frame = self.get(per=period)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Profile x-axis column {x_column!r} was not found.")
        fig.add_scattergl(
            x=frame[x_column],
            y=frame[self.y_column],
            mode="lines+markers",
            name=f"SFR {self.y_label} per {period}",
            customdata=np.column_stack([frame["reach"], frame["cell"]]),
            hovertemplate=(
                "reach=%{customdata[0]}<br>"
                "cell=%{customdata[1]}<br>"
                f"{self.y_label}=%{{y}}<extra></extra>"
            ),
        )
        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title=self.y_label,
            title=f"SFR {self.y_label} profile (per={period})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None


class SfrBudgetResultsExplorer(CellBudgetResultsExplorer):
    """SFR-specific result explorer with reach-profile helpers."""

    def __init__(
        self,
        model: SimulationBase,
        *,
        budget_text: str = "SFR",
        value_name: str = "q_gwf",
        package: str | None = None,
    ):
        """Bind an SFR budget-result explorer (defaults to the ``SFR`` term's ``q_gwf``).

        ``value_name`` is the emitted column (``q_gwf`` -- the SFR exchange is
        aquifer-referenced); the accessor stays ``results.q`` via ``result_name``.
        ``package`` is the SFR package's NAME (``None``: the model's only one).
        """

        super().__init__(model, "sfr", budget_text, value_name, result_name="q", package=package)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized SFR budget-result table for selected rows.

        The GWF model budget's ``SFR`` records -- the stream-aquifer exchange, in
        MODFLOW 6's own aquifer-referenced sign (``q_gwf`` NEGATIVE where the reach
        GAINS), model length³/time -- one row per record at every saved step,
        joined to the reach table for ``reach``, ``rlen`` and the ``distance_*``
        columns. **The join is on the CELL, not on the record's own reach**
        (``node2``): a cell holding several reaches gets every one of its records
        once per reach, each tagged with every reach -- measured on the testing
        canonical model, cell 195 holds reaches 22, 37 and 38 and carries 9 rows
        for its 3 records. Read a multi-reach cell's flows through ``node2``.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period to keep, every saved step of it -- a FILTER.
            One period (a list raises ``TypeError``); ``None`` (default) keeps
            every period, and a period the file saved nothing in gives an empty
            table.
        layer : int or iterable of int or None, optional
            Zero-based layer(s) of the groundwater cell to keep. ``None``
            (default) keeps every layer.
        cells : int or iterable of int or None, optional
            Zero-based cell number(s) within a layer (the ``cell`` column), ANDed
            with ``layer``. ``None`` (default) keeps every cell.

        Returns
        -------
        pandas.DataFrame
            :class:`CellBudgetResultsExplorer`'s columns (``q_gwf`` the value)
            plus ``reach``, ``rlen``, ``distance_start``/``distance_mid``/
            ``distance_end`` and ``q_per_length`` (``q_gwf / rlen``).
        """

        frame = build_sfr_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
            package=self.package,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    @property
    def profile(self) -> SfrReachProfileView:
        """Return the reach exchange-profile view: ``.get()`` table, ``.plot()`` line, ``.summary()``.

        The single-field (``q_gwf``) profile noun. For the merged multi-field
        profile (stage + streambed + exchange) use ``sfr.results.profile``.
        """

        return SfrReachProfileView(
            self, y_column=self.value_name, y_label=self.value_name, label="sfr.results.q.profile"
        )

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        per: int = 0,
        layer: int = 0,
        # -- table to cells --------------------------------------------------
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing --------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover -----------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """Stream-aquifer exchange as a choropleth, per unit of reach length.

        A RESULTS noun over the SFR ``GWF`` budget term. The mapped value is
        ``sum(q) / sum(rlen)`` within each cell, which is what keeps the picture
        about the physics: a raw volumetric exchange makes a cell look busy
        merely for containing a longer piece of stream.

        **The signs are MODFLOW 6's own.** SFR's cell record is flow FROM the
        reach TO the groundwater cell, so a POSITIVE value is a LOSING reach and
        a negative one is gaining. The diverging colorscale is derived from the
        package registry's frame rather than written as a literal, so it cannot
        drift from the data -- the same defect once shipped inverted on the LAK
        map because a frame literal had been copied across.

        Because the quantity is SIGNED, the scale is centred on zero:
        ``zmin``/``zmax`` default to ``-max|q|``/``+max|q|`` and the trace gets
        ``zmid=0``, which is what makes gaining and losing different colours
        rather than two shades of one. Pass ``zmin``/``zmax`` yourself to pin a
        scale across several periods -- that also disables the symmetric
        default, so pass both.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into the tooltip on every cell. Passed anyway, they ride
        ``**trace_kwargs`` to the free verb. ``show_layer_elevs`` and
        ``hover_surfaces`` ARE named (ledger 306): the layer frame is the grid's,
        so they act on any map.

        Parameters
        ----------
        per : int, default 0
            Stress period to draw, zero-based, at its LAST saved step on the SFR
            budget file's own clock, the hover's stage at the stage file's (ledger
            337's rule): at ``nstp=3`` under ``SAVE BUDGET ALL`` the map ADDED the
            period's three steps, and averaged its stages, until ledger 348. One
            period per picture; use :meth:`animate` to flip through them.
        layer : int, default 0
            Zero-based layer. Reaches connected in other layers are not drawn.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, or ``-1.0`` to
            read gaining as positive. That changes the PICTURE only: ``.get()``
            keeps MF6's convention, so the colorbar will then disagree with the
            table, and the caption should say so. Applied to every cell after the
            per-cell normalization, ``fill_value`` included -- except when the
            selection holds no record at all, which draws the fill unscaled.
        fill_value : float, default 0.0
            Value given to cells with no stream reach in them. Zero is honest for
            an exchange -- no reach, no exchange -- but ``float("nan")`` leaves
            them uncoloured, which reads better on a sparse network.
        agg : str, default 'sum'
            **Accepted and not used.** The normalization is computed explicitly
            from total exchange over total length per cell, so there is no
            free choice of reducer left to make. Kept in the signature because
            removing it would break callers that pass it; it is a no-op, not a
            silent alternative.
        **trace_kwargs
            Every other keyword, forwarded to ``model.plot.map`` after two checks:
            ``values``, ``type``, ``custom_hover``, ``hover_heads`` and ``hover_ks``
            raise ``TypeError`` (this noun supplies its own field and hover), and
            ``hover_spec`` is read as the old spelling of ``hover``, which wins when
            both are given. A free-verb parameter this noun does not name -- ``zoom``,
            ``bgs``, ``kstpkper``, ``per_timestep``, ``hover_layers``,
            ``show_mounding``, ... -- reaches the verb by that name, and ``custom_zs``,
            which this noun computes, raises ``TypeError`` as a second value. The rest
            are properties of the Plotly ``Choroplethmap`` trace -- ``zmid``,
            ``colorbar``, ``reversescale``, ``showscale`` -- checked only when the
            figure is built, so a misspelled one raises ``ValueError`` at ``.fig`` or
            ``.show()``, and is IGNORED under ``backend="mpl"``, which never builds
            that trace. ``zmid`` defaults to ``0.0`` here, the centre of the
            symmetric scale; pass one to move it.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        get : the per-reach exchanges behind the picture, as a DataFrame.
        profile : the same term along the stream rather than in plan view.
        myflopy.modflow.mf6.package_surface_water.SurfaceWaterExchangeResultsExplorer.map :
            SFR and LAK together on one shared scale.

        Examples
        --------
        >>> model.sfr.results.q.map()
        >>> model.sfr.results.q.map(per=5)
        >>> model.sfr.results.q.map(zmin=-0.5, zmax=0.5)
        >>> model.sfr.results.q.map(fill_value=float("nan"))
        >>> model.sfr.results.q.map(select="all_streams", select_style="both")
        >>> model.sfr.results.q.map(backend="mpl").savefig("sfr_q.png")
        """

        del agg
        refuse_noun_parameters(self.label, self.value_name, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = join_sfr_stage(
            self.model, self.get(per=per, layer=layer), per=per, package=self.package
        )
        values, cell_hover = build_sfr_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            value_column=self.value_name,
        )
        absmax = _symmetric_color_limit(values)
        if absmax > 0:
            zmin = -absmax if zmin is None else zmin
            zmax = absmax if zmax is None else zmax
        # `zmid` is a Plotly trace property, not a verb parameter, so it rides
        # the open tail the way `colorbar`/`reversescale` do.
        trace_kwargs.setdefault("zmid", 0.0)
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
            # In the model's own units, never assumed feet (ledger 434).
            hover_spec=sfr_hover(**_model_units(self.model)),
            zmin=zmin,
            zmax=zmax,
            # SFR's cell record is flow FROM reach TO cell (the "gwf" frame), so
            # gaining is negative. The orientation is derived from the registry,
            # not assumed, so it cannot drift from the data.
            colorscale=colorscale or _exchange_colorscale(_exchange_frame("sfr")),
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)

    # -- backing method for the pre-view spelling -----------------------------
    # ``plot_profile`` became ``profile.plot``. This body preserves the OLD return
    # (the exchange line figure); it resolves solely through __getattr__ so the
    # retired spelling stays out of dir()/completion (D12). Named ``_legacy_*``
    # rather than echoing the old name, which would resurface it in completion.
    def _legacy_profile_plot(self, **kwargs):
        """Back the retired ``plot_profile`` spelling; returns the same figure."""

        return self.profile.plot(**kwargs)

    __getattr__ = deprecated_instance_getattr(
        {
            "plot_profile": (
                "_legacy_profile_plot",
                "model.sfr.results.q.profile.plot",
                "0.1.0",
            ),
        },
        "myflopy.modflow.mf6.package_surface_water.SfrBudgetResultsExplorer",
    )


class LakBudgetView:
    """LAK groundwater-exchange budget by connection type: table + bar figure.

    Reached from the LAK exchange field (``lak.results.q.budget``); follows the
    house view shape (``docs/view_layer_conventions.md``) -- ``get`` for the
    per-lake / per-connection-type summary table, ``plot`` for the bar figure,
    ``summary`` for a compact digest. Calling the view rebinds the stress period.
    """

    def __init__(self, explorer, *, per: int = 0):
        """Bind the budget view to a LAK exchange explorer at period ``per``."""

        self.explorer = explorer
        self.per = int(per)

    def __repr__(self) -> str:
        """Show the bound period (the view's only state)."""

        return f"{type(self).__name__}(per={self.per})"

    def __call__(self, *, per: int) -> LakBudgetView:
        """Return an equivalent view bound to stress period ``per``.

        Parameters
        ----------
        per : int
            Zero-based stress period the new view reads when its verbs are called
            without one; required, keyword-only, passed through ``int()`` (so
            ``None`` raises ``TypeError``). This view is left unchanged.

        Returns
        -------
        LakBudgetView
            A view of the same LAK exchange bound to ``per``.
        """

        return type(self)(self.explorer, per=per)

    def _period(self, per: int | None) -> int:
        """Resolve an explicit ``per`` against the period bound to this view."""

        return self.per if per is None else int(per)

    def get(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return the exchange summary grouped by period, lake, and connection type.

        ONE period, at its LAST saved step on the LAK budget file's own clock: the
        connections of that step add within each (lake, connection type).

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period. ``None`` (default) is the period bound to the
            view (``budget(per=N)``, 0 unless rebound) -- never every period. A
            period the file saved nothing in gives an empty table.
        connection_type : str or iterable of str or None, optional
            LAK connection type(s) to keep -- ``"vertical"``, ``"horizontal"``,
            ``"embeddedh"``, ``"embeddedv"``, any case. ``None`` (default) keeps
            every type. A type no connection has gives an empty table, not an
            error.

        Returns
        -------
        pandas.DataFrame
            ``per``, ``lake``, ``claktype``, ``record_count`` (connections), the
            summed exchange ``q_lake`` (model length³/time, the LAKE's sign:
            positive is the lake gaining), the summed ``flow_area`` (length²) and
            ``q_per_area`` (``q_lake / flow_area``, length/time; ``NaN`` where the
            area is zero).
        """

        return self.explorer._build_budget_summary(
            per=self._period(per), connection_type=connection_type
        )

    def summary(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return a compact digest of the budget summary.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period, as :meth:`get` reads it. ``None`` (default) is
            the period bound to the view.
        connection_type : str or iterable of str or None, optional
            LAK connection type(s) to keep -- ``"vertical"``, ``"horizontal"``,
            ``"embeddedh"``, ``"embeddedv"``, any case. ``None`` (default) keeps
            every type. A type no connection has gives an empty table, not an
            error.

        Returns
        -------
        pandas.DataFrame
            One row: ``label``, ``records`` (rows of :meth:`get`), ``lakes`` and
            ``connection_types`` (distinct counts), zeros when nothing matches.
        """

        frame = self.get(per=per, connection_type=connection_type)
        if frame.empty:
            return pd.DataFrame(
                [{"label": "lak.results.q.budget", "records": 0, "lakes": 0, "connection_types": 0}]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.results.q.budget",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "connection_types": int(frame["claktype"].nunique()),
                }
            ]
        )

    def plot(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
        value: str = "q_per_area",
        ax=None,
        return_fig: bool = True,
    ):
        """Plot a compact LAK budget summary by connection type for one period.

        Positive (the lake gains / inflow) draws blue, negative (loses / outflow)
        red -- the house rule for the feature-referenced LAK exchange.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period, as :meth:`get` reads it. ``None`` (default) is
            the period bound to the view. A period with no rows draws labelled,
            empty axes rather than raising.
        connection_type : str or iterable of str or None, optional
            LAK connection type(s) to keep -- ``"vertical"``, ``"horizontal"``,
            ``"embeddedh"``, ``"embeddedv"``, any case. ``None`` (default) keeps
            every type. A type no connection has gives an empty table, not an
            error.
        value : {"q_per_area", "q", "q_lake", "flow_area"}, default "q_per_area"
            The :meth:`get` column drawn, one bar per (lake, connection type):
            ``"q_per_area"`` (length/time), ``"q"`` -- read as ``"q_lake"`` -- the
            summed exchange (length³/time), or ``"flow_area"`` (length², drawn in
            one colour). Anything else raises ``ValueError``.
        ax : matplotlib.axes.Axes or None, optional
            Axes to draw onto; the figure returned is then ``ax.figure``, and
            ``tight_layout()`` is called on it whenever bars are drawn (not on an
            empty period). ``None`` (default) makes a new 8 x 4 inch figure.
        return_fig : bool, default True
            Return the figure. With ``False`` the method returns ``None``, which
            suits drawing onto your own ``ax``.

        Returns
        -------
        matplotlib.figure.Figure or None
            The figure the bars are on, or ``None`` when ``return_fig=False``.

        Raises
        ------
        ValueError
            When ``value`` is not one of the four names.
        """

        period = self._period(per)
        summary = self.get(per=period, connection_type=connection_type)
        q = self.explorer.value_name
        # Accept the accessor-facing "q" as an alias for the frame-named column.
        if value == "q":
            value = q
        if value not in {q, "flow_area", "q_per_area"}:
            raise ValueError(
                f"value must be one of: 'q', {q!r}, 'flow_area', 'q_per_area'"
            )
        if ax is None:
            fig, ax = mpl_axes(figsize=(8, 4))
        else:
            fig = ax.figure
        if summary.empty:
            ax.set_title(f"LAK {value} summary (per={period})")
            ax.set_xlabel("Connection Type")
            ax.set_ylabel(value)
            if return_fig:
                return fig
            return None

        summary = summary.copy()
        labels = summary.apply(
            lambda row: f"Lake {int(row['lake'])}\n{row['claktype']}",
            axis=1,
        )
        colors = []
        if value == "flow_area":
            colors = ["#4c78a8" for _ in range(len(summary))]
        else:
            for current in pd.to_numeric(summary[value], errors="coerce").fillna(0.0):
                colors.append("#1f77b4" if current >= 0.0 else "#d62728")
        ax.bar(labels, summary[value].astype(float), color=colors)
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
        ax.set_title(f"LAK {value} summary (per={period})")
        ax.set_xlabel("Lake / Connection Type")
        ax.set_ylabel(value)
        ax.tick_params(axis="x", rotation=0)
        fig.tight_layout()
        if return_fig:
            return fig
        return None


class LakBudgetResultsExplorer(CellBudgetResultsExplorer):
    """LAK-specific result explorer with area-normalized exchange maps."""

    def __init__(
        self,
        model: SimulationBase,
        *,
        budget_text: str = "GWF",
        value_name: str = "q_lake",
        package: str | None = None,
    ):
        """Bind a LAK budget-result explorer (defaults to the ``GWF`` exchange ``q_lake``).

        ``value_name`` is the emitted column (``q_lake`` -- LAK's exchange is
        feature-referenced); the accessor stays ``results.q`` via ``result_name``.
        ``package`` is the LAK package's NAME (``None``: the model's only one).
        """

        super().__init__(model, "lak", budget_text, value_name, result_name="q", package=package)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK budget-result table for selected rows.

        The LAK package budget file's ``GWF`` records -- one row per lake-aquifer
        connection at every saved step -- matched back to the connection data for
        ``claktype``, ``belev``/``telev``, ``connlen``/``connwidth`` and
        ``bedleak``. ``q_lake`` is in the LAKE's sign (positive: the lake gains
        from the aquifer), model length³/time.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period to keep, every saved step of it -- a FILTER.
            One period (a list raises ``TypeError``); ``None`` (default) keeps
            every period, and a period the file saved nothing in gives an empty
            table.
        layer : int or iterable of int or None, optional
            Zero-based layer(s) of the connected groundwater cell to keep. ``None``
            (default) keeps every layer.
        cells : int or iterable of int or None, optional
            Zero-based cell number(s) within a layer of the connected cell (the
            ``cell`` column), ANDed with ``layer``. ``None`` (default) keeps every
            cell.
        connection_type : str or iterable of str or None, optional
            LAK connection type(s) to keep -- ``"vertical"``, ``"horizontal"``,
            ``"embeddedh"``, ``"embeddedv"``, any case -- ANDed with the other
            selectors. ``None`` (default) keeps every type; a type no connection
            has gives an empty table, not an error.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``kstpkper``, ``per``, ``lake``, ``layer``,
            ``cell``, ``claktype``, the connection geometry, ``q_lake``,
            ``flow_area`` and ``q_per_area`` (``q_lake / flow_area``,
            length/time; ``NaN`` where the area is zero). ``flow_area`` (length²)
            is MODFLOW 6's own ``FLOW-AREA`` for the connection at that step --
            the WETTED area, not the connection table's ``connection_area``: on
            the testing canonical model a vertical connection reads its cell's
            whole plan area or 0.0, 11 to 14 of 33 of them 0.0 depending on the
            period.
        """

        frame = build_lak_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
            package=self.package,
        )
        frame = _filter_normalized_table(frame, per=per, layer=layer, cells=cells)
        connection_types = _normalize_connection_type_filter(connection_type)
        if connection_types is not None and "claktype" in frame.columns:
            frame = frame.loc[
                frame["claktype"].astype("string").str.upper().isin(connection_types)
            ].copy()
        return frame.reset_index(drop=True)

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        per: int = 0,
        layer: int = 0,
        connection_type: str | Iterable[str] | None = None,
        # -- table to cells --------------------------------------------------
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing --------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover -----------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """Lake-aquifer exchange as a choropleth, per unit of connection area.

        A RESULTS noun over the LAK ``GWF`` budget term. The mapped value is
        ``sum(q) / sum(flow_area)`` within each cell, giving a signed exchange
        INTENSITY in length-per-time rather than a raw volumetric exchange that
        would scale with how much lakebed a cell happens to hold.

        **The signs are MODFLOW 6's own, and they are NOT the SFR frame.** LAK's
        ``GWF`` record comes from the LAK package budget, written from the LAKE's
        perspective, so a POSITIVE value means the lake GAINS from the aquifer --
        the opposite of SFR's cell-side record. This map once shipped inverted,
        losing lakes drawing blue, because a frame literal had been copied across
        from SFR; the orientation is now derived from the package registry, which
        removes that failure mode rather than fixing one instance of it.

        Because the quantity is SIGNED, the scale is centred on zero:
        ``zmin``/``zmax`` default to ``-max|q|``/``+max|q|`` and the trace gets
        ``zmid=0``, which is what makes gaining and losing different colours
        rather than two shades of one. Pass ``zmin``/``zmax`` yourself to pin a
        scale across several periods -- that also disables the symmetric
        default, so pass both.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into the tooltip on every cell. Passed anyway, they ride
        ``**trace_kwargs`` to the free verb. ``show_layer_elevs`` and
        ``hover_surfaces`` ARE named (ledger 306): the layer frame is the grid's,
        so they act on any map.

        Parameters
        ----------
        per : int, default 0
            Stress period to draw, zero-based, at its LAST saved step on the LAK
            budget file's own clock, the hover's stage at the stage file's (ledger
            337's rule): at ``nstp=3`` under ``SAVE BUDGET ALL`` the map ADDED the
            period's three flows and three copies of each area, and the hover
            joined a stage from an EARLIER period, until ledger 348. One period per
            picture; use :meth:`animate` to flip through them.
        layer : int, default 0
            Zero-based layer. A lake usually connects across several, so this is
            how you pick the one you mean -- ``.get()`` shows which exist.
        connection_type : str or iterable of str, optional
            Restrict to particular LAK connection types (``"vertical"``,
            ``"horizontal"``, ``"embeddedh"``, ``"embeddedv"``). With none, every
            connection in the layer contributes. Worth setting when a lake's
            vertical lakebed leakage and its horizontal shoreline exchange are
            physically different things you do not want summed into one colour.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, or ``-1.0`` to
            flip the sign convention. That changes the PICTURE only: ``.get()``
            keeps MF6's convention, so say so in the caption if you use it.
            Applied to every cell after the per-cell normalization,
            ``fill_value`` included -- except when the selection holds no record
            at all, which draws the fill unscaled.
        fill_value : float, default 0.0
            Value given to cells with no lake connection. Zero is honest for an
            exchange, but ``float("nan")`` leaves them uncoloured, which reads
            better when the lakes cover a small part of the grid.
        agg : str, default 'sum'
            **Accepted and not used.** The normalization is computed explicitly
            from total exchange over total flow area per cell, so no free choice
            of reducer remains. Kept in the signature because removing it would
            break callers that pass it; it is a no-op, not a silent alternative.
        **trace_kwargs
            Every other keyword, forwarded to ``model.plot.map`` after two checks:
            ``values``, ``type``, ``custom_hover``, ``hover_heads`` and ``hover_ks``
            raise ``TypeError`` (this noun supplies its own field and hover), and
            ``hover_spec`` is read as the old spelling of ``hover``, which wins when
            both are given. A free-verb parameter this noun does not name -- ``zoom``,
            ``bgs``, ``kstpkper``, ``per_timestep``, ``hover_layers``,
            ``show_mounding``, ... -- reaches the verb by that name, and ``custom_zs``,
            which this noun computes, raises ``TypeError`` as a second value. The rest
            are properties of the Plotly ``Choroplethmap`` trace -- ``zmid``,
            ``colorbar``, ``reversescale``, ``showscale`` -- checked only when the
            figure is built, so a misspelled one raises ``ValueError`` at ``.fig`` or
            ``.show()``, and is IGNORED under ``backend="mpl"``, which never builds
            that trace. ``zmid`` defaults to ``0.0`` here, the centre of the
            symmetric scale; pass one to move it.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        get : the per-connection exchanges behind the picture, as a DataFrame.
        budget : the lake's whole water balance, term by term.
        myflopy.modflow.mf6.package_surface_water.SurfaceWaterExchangeResultsExplorer.map :
            SFR and LAK together on one shared scale and one sign convention.

        Examples
        --------
        >>> model.lak.results.q.map()
        >>> model.lak.results.q.map(per=5, connection_type="vertical")
        >>> model.lak.results.q.map(layer=1, fill_value=float("nan"))
        >>> model.lak.results.q.map(zmin=-0.2, zmax=0.2)
        >>> model.lak.results.q.map(contours=True, contour_levels=6)
        >>> model.lak.results.q.map(backend="mpl").savefig("lak_q.png")
        """

        del agg
        refuse_noun_parameters(self.label, self.value_name, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = join_lak_stage(
            self.model,
            self.get(per=per, layer=layer, connection_type=connection_type),
            per=per,
            package=self.package,
        )
        values, cell_hover = build_lak_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            value_column=self.value_name,
        )
        absmax = _symmetric_color_limit(values)
        if absmax > 0:
            zmin = -absmax if zmin is None else zmin
            zmax = absmax if zmax is None else zmax
        # `zmid` is a Plotly trace property, not a verb parameter, so it rides
        # the open tail the way `colorbar`/`reversescale` do.
        trace_kwargs.setdefault("zmid", 0.0)
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
            # In the model's own units, never assumed feet (ledger 434).
            hover_spec=lak_hover(**_model_units(self.model)),
            zmin=zmin,
            zmax=zmax,
            # NOT the SFR frame: LAK's GWF record comes from the LAK package
            # budget, written from the LAKE's perspective (the "feature" frame),
            # so gaining is POSITIVE. This map once shipped inverted -- losing
            # lakes drew blue -- because a frame literal was copied from SFR;
            # deriving it from the registry removes that whole failure mode.
            colorscale=colorscale or _exchange_colorscale(_exchange_frame("lak")),
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)

    @property
    def budget(self) -> LakBudgetView:
        """Return the budget-summary view: ``.get()`` table, ``.plot()`` bar, ``.summary()``.

        Summarizes lake-groundwater exchange by period, lake, and connection type.
        """

        return LakBudgetView(self)

    def _build_budget_summary(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Summarize lake-groundwater exchange by period, lake, and connection type.

        Each period is its LAST saved step on the LAK budget file's own clock
        (ledgers 337, 348): the connections of that one step add. Every step the
        period saved was added until ledger 348 -- three steps' flows, and three
        copies of each connection's area and count, at ``nstp=3`` under ``SAVE
        BUDGET ALL``.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period filter. When omitted, all periods
            are summarized.
        connection_type : str or iterable of str or None, optional
            LAK connection types to keep (``"VERTICAL"``, ``"HORIZONTAL"``, ...,
            any case). ``None`` keeps every type.

        Returns
        -------
        pandas.DataFrame
            Summary table with signed volumetric exchange ``q`` and
            area-normalized exchange ``q_per_area`` grouped by period, lake,
            and connection type.
        """

        frame = at_last_saved_steps(self.get(per=per, connection_type=connection_type))
        q = self.value_name  # feature-referenced exchange column (q_lake)
        if frame.empty:
            return pd.DataFrame(
                columns=[
                    "per",
                    "lake",
                    "claktype",
                    "record_count",
                    q,
                    "flow_area",
                    "q_per_area",
                ]
            )
        summary = (
            frame.groupby(["per", "lake", "claktype"], dropna=False, as_index=False)
            .agg(
                record_count=("cell", "size"),
                **{q: (q, "sum")},
                flow_area=("flow_area", "sum"),
            )
            .sort_values(["per", "lake", "claktype"])
            .reset_index(drop=True)
        )
        summary["q_per_area"] = np.where(
            pd.to_numeric(summary["flow_area"], errors="coerce") > 0.0,
            pd.to_numeric(summary[q], errors="coerce")
            / pd.to_numeric(summary["flow_area"], errors="coerce"),
            np.nan,
        )
        return summary

    # -- backing methods for the pre-view spellings ---------------------------
    # ``budget_summary``/``plot_budget`` became ``budget.get``/``budget.plot``.
    # These bodies preserve the OLD returns exactly (a DataFrame and the bar
    # figure) and resolve solely through __getattr__ so the retired spellings
    # stay out of dir()/completion (D12). Named ``_legacy_*`` to keep the old
    # spellings from resurfacing in IDE completion.
    def _legacy_budget_summary(self, **kwargs) -> pd.DataFrame:
        """Back the retired ``budget_summary`` spelling; returns the same table."""

        return self.budget.get(**kwargs)

    def _legacy_budget_plot(self, **kwargs):
        """Back the retired ``plot_budget`` spelling; returns the same figure."""

        return self.budget.plot(**kwargs)

    __getattr__ = deprecated_instance_getattr(
        {
            "budget_summary": (
                "_legacy_budget_summary",
                "model.lak.results.q.budget.get",
                "0.1.0",
            ),
            "plot_budget": (
                "_legacy_budget_plot",
                "model.lak.results.q.budget.plot",
                "0.1.0",
            ),
        },
        "myflopy.modflow.mf6.package_surface_water.LakBudgetResultsExplorer",
    )


class LakStageResultsExplorer(StageResultsExplorer):
    """LAK stage explorer -- the unified grammar's ``plot()`` draws one line
    per lake by stress period (replaced the old ``plot_timeseries``)."""

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind a LAK stage-result explorer using the lake stage table builder.

        ``package`` is the LAK package's NAME (``None``: the model's only one).
        """

        super().__init__(model, "lak", build_lak_stage_result_table, package=package)


class LakStageChangeExplorer:
    """Explorer for lake-stage changes between stress periods."""

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the lake stage-change explorer to ``model``'s LAK package ``package``.

        ``package`` is the package's NAME (``None``: the model's only LAK).
        """

        self.model = model
        self.package = package

    def get(
        self,
        *,
        lake: int | None = None,
        per0: int | None = None,
        per1: int | None = None,
    ) -> pd.DataFrame:
        """Return stage-change rows for one lake and/or one period transition.

        One row per lake per pair of CONSECUTIVE periods the stage file saved,
        each period read at its last saved step (``stage_change = stage1 -
        stage0``, model length units; ``NaN`` into or out of a dry or inactive
        period, the row kept). The three selectors combine by AND; each matches
        exactly, so a value no row has gives an empty table, not an error.

        Parameters
        ----------
        lake : int or None, optional
            Zero-based lake number (within this LAK package). One lake: a list
            raises ``TypeError``. ``None`` (default) keeps every lake.
        per0 : int or None, optional
            Zero-based stress period the transition starts FROM. ``None``
            (default) keeps every start.
        per1 : int or None, optional
            Zero-based stress period the transition ends IN -- the next period the
            stage file saved after ``per0``. ``None`` (default) keeps every end.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``lake``, ``per0``, ``per1``, ``stage0``,
            ``stage1`` and ``stage_change``.
        """

        frame = build_lak_stage_change_table(self.model, package=self.package)
        if lake is not None:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        if per0 is not None:
            frame = frame.loc[frame["per0"] == int(per0)].copy()
        if per1 is not None:
            frame = frame.loc[frame["per1"] == int(per1)].copy()
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available lake-stage transitions."""

        frame = self.get()
        if frame.empty:
            return pd.DataFrame(
                [
                    {
                        "label": "lak.results.stage_change",
                        "records": 0,
                        "lakes": 0,
                        "transitions": 0,
                    }
                ]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.results.stage_change",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "transitions": int(
                        frame[["per0", "per1"]].drop_duplicates().shape[0]
                    ),
                }
            ]
        )

    def plot(
        self,
        *,
        lake: int | None = None,
        backend: str = "plotly",
        title: str | None = None,
    ):
        """Series panel: stage change by stress-period transition, per lake.

        The grammar's ``plot`` verb for this node (its x-axis is the period
        *transition*, not the period, so it does not use the generic series
        engine). ``backend="plotly"`` returns a ``viz.Fig``; ``backend="mpl"``
        a matplotlib figure.

        Parameters
        ----------
        lake : int or None, optional
            Zero-based lake number to draw, as :meth:`get` selects it. ``None``
            (default) draws one line per lake. A lake with no rows draws empty
            axes rather than raising.
        backend : {"plotly", "mpl"}, default "plotly"
            Any case. ``"plotly"`` (or ``"interactive"``) returns a ``viz.Fig``;
            ``"mpl"`` (or ``"matplotlib"``, ``"static"``) a Matplotlib figure, 8 x
            4 inches, with a zero line. Anything else raises ``ValueError``.
        title : str or None, optional
            The figure title. ``None`` (default) -- or an empty string -- is
            ``"LAK stage change by transition"``.

        Returns
        -------
        myflopy.viz.Fig or matplotlib.figure.Figure
            One line per lake, x the ``"per0->per1"`` transition label and y
            ``stage_change`` (model length units).

        Raises
        ------
        ValueError
            When ``backend`` is not one of the five spellings.
        """

        from myflopy.viz import Fig

        frame = self.get(lake=lake)
        heading = title or "LAK stage change by transition"
        normalized = str(backend).lower()
        if normalized in ("plotly", "interactive"):
            fig = Fig()
            for lake_id, group in frame.groupby("lake", dropna=False):
                labels = [
                    f"{int(start)}->{int(end)}"
                    for start, end in zip(group["per0"], group["per1"], strict=False)
                ]
                fig.add_scatter(
                    x=labels,
                    y=group["stage_change"].astype(float).to_numpy(),
                    mode="lines+markers",
                    name=f"Lake {int(lake_id)}",
                )
            fig.update_layout(
                title=heading,
                xaxis_title="Stress-Period Transition",
                yaxis_title="Stage Change",
            )
            return fig
        if normalized not in ("mpl", "matplotlib", "static"):
            raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")
        fig, ax = mpl_axes(figsize=(8, 4))
        for lake_id, group in frame.groupby("lake", dropna=False):
            labels = [
                f"{int(start)}->{int(end)}"
                for start, end in zip(group["per0"], group["per1"], strict=False)
            ]
            ax.plot(
                labels,
                group["stage_change"].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"Lake {int(lake_id)}",
            )
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
        ax.set_title(heading)
        ax.set_xlabel("Stress-Period Transition")
        ax.set_ylabel("Stage Change")
        if not frame.empty:
            ax.legend()
        fig.tight_layout()
        return fig


class LakConnectionsExplorer:
    """Explorer for LAK connection geometry and exchange interface area."""

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the LAK connection-geometry explorer to ``model``'s LAK package ``package``.

        ``package`` is the package's NAME (``None``: the model's only LAK).
        """

        self.model = model
        self.package = package

    def get(
        self,
        *,
        lake: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK connection table.

        One row per connection record of ``connectiondata`` -- input geometry, no
        time axis. The selectors combine by AND; a value no row has gives an empty
        table, not an error.

        Parameters
        ----------
        lake : int or None, optional
            Zero-based lake number (within this LAK package). One lake: a list
            raises ``TypeError``. ``None`` (default) keeps every lake.
        layer : int or iterable of int or None, optional
            Zero-based layer(s) of the connected groundwater cell. ``None``
            (default) keeps every layer.
        cells : int or iterable of int or None, optional
            Zero-based cell number(s) within a layer of the connected cell (the
            ``cell`` column). ``None`` (default) keeps every cell.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``lake``, ``iconn``, ``layer``, ``cell``,
            ``claktype``, ``belev``, ``telev``, ``connlen``, ``connwidth``,
            ``connection_area`` (length²: the cell's plan area for a vertical
            connection, ``connwidth * (telev - belev)`` -- the INPUT elevations,
            a negative thickness read as zero -- for every other type, embedded
            ones included, although MODFLOW 6 reads ``connwidth`` on a HORIZONTAL
            connection only and accepts any value on an embedded one) and
            ``bedleak``.
        """

        frame = build_lak_connection_table(self.model, package=self.package)
        frame = _filter_normalized_table(frame, per=None, layer=layer, cells=cells)
        if lake is not None and "lake" in frame.columns:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of LAK connection geometry."""

        frame = self.get()
        if frame.empty:
            return pd.DataFrame(
                [
                    {
                        "label": "lak.connections",
                        "records": 0,
                        "lakes": 0,
                        "layers": 0,
                        "cells": 0,
                        "total_connection_area": 0.0,
                    }
                ]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.connections",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "layers": int(frame["layer"].nunique()),
                    "cells": int(frame["cell"].nunique()),
                    "total_connection_area": float(
                        pd.to_numeric(frame["connection_area"], errors="coerce").sum()
                    ),
                }
            ]
        )

    def map(
        self,
        *,
        # -- which connections -----------------------------------------------
        lake: int | None = None,
        layer: int = 0,
        value_column: str = "connection_area",
        # -- table to cells --------------------------------------------------
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing --------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover -----------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """LAK connection GEOMETRY as a choropleth: one colour per cell.

        An INPUTS noun. It draws how a lake is wired to the aquifer -- the area
        or width of each connection -- not what flows through it. For the flow,
        see ``model.lak.results.q.map()``.

        Connection geometry has no time axis, so there is no ``per``: a lake's
        connections are declared once and do not change with the stress period.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into the tooltip on every cell. Passed anyway, they ride
        ``**trace_kwargs`` to the free verb. ``show_layer_elevs`` and
        ``hover_surfaces`` ARE named (ledger 306): the layer frame is the grid's,
        so they act on any map.

        Parameters
        ----------
        lake : int, optional
            Restrict to one lake, by zero-based lake id. With none, every lake's
            connections are drawn together -- fine for seeing the whole lake
            system, misleading if two lakes share a cell and you meant one.
        layer : int, default 0
            Zero-based layer to render. A lake usually connects across several;
            ``.get()`` shows which, in its ``layer`` column.
        value_column : str, default 'connection_area'
            Which connection field carries the colour. ``"connection_area"`` and
            ``"connwidth"`` are the usual choices; anything present in
            ``.get()`` works, and a name that is not raises ``KeyError`` naming
            it rather than drawing an empty map.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, default 'sum'
            How several connections landing in ONE cell are reduced to one
            number. ``sum`` is right for an area -- two connections in a cell
            wet the sum of their areas. Use ``mean`` for a width, where adding
            them would invent a wider connection than any that exists.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, most often.
            Applied AFTER ``agg`` (so under ``-1.0``, ``agg="max"`` draws
            minus the largest record) and to every cell, ``fill_value`` included:
            ``fill_value=1.0, multiplier=2.0`` draws 2.0 where there is no record
            -- but 1.0 when the selection holds no record at all.
        fill_value : float, default 0.0
            Value given to cells with no lake connection. Pass ``float("nan")``
            to leave them uncoloured, which is usually clearer here: a lake
            covers a small part of the grid, and a 0.0 backdrop takes over the
            colour scale.
        **trace_kwargs
            Every other keyword, forwarded to ``model.plot.map`` after two checks:
            ``values``, ``type``, ``custom_hover``, ``hover_heads`` and ``hover_ks``
            raise ``TypeError`` (this noun supplies its own field and hover), and
            ``hover_spec`` is read as the old spelling of ``hover``, which wins when
            both are given. A free-verb parameter this noun does not name -- ``zoom``,
            ``bgs``, ``kstpkper``, ``per_timestep``, ``hover_layers``,
            ``show_mounding``, ... -- reaches the verb by that name, and ``custom_zs``,
            which this noun computes, raises ``TypeError`` as a second value -- and so
            does ``per``: connection geometry has no time axis, this signature names
            no ``per``, and the noun passes ``per=0`` itself (the ``per`` entry below
            is spliced from the free verb and does not apply here). The rest are
            properties of the Plotly ``Choroplethmap`` trace -- ``zmid``,
            ``colorbar``, ``reversescale``, ``showscale`` -- checked only when the
            figure is built, so a misspelled one raises ``ValueError`` at ``.fig`` or
            ``.show()``, and is IGNORED under ``backend="mpl"``, which never builds
            that trace.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        per : int, optional
            Stress period to read (0-based). Mutually exclusive with ``kstpkper``;
            with neither, the model's first output time is used.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        get : the connection records behind the picture, as a DataFrame.
        summary : the same records reduced to one row per lake.
        myflopy.modflow.mf6.package_surface_water.LakBudgetResultsExplorer.map :
            what actually flows through these connections.

        Examples
        --------
        >>> model.lak.connections.map()
        >>> model.lak.connections.map(lake=0, layer=1)
        >>> model.lak.connections.map(value_column="connwidth", agg="mean")
        >>> model.lak.connections.map(fill_value=float("nan"))
        >>> model.lak.connections.map(logscale=True)
        >>> model.lak.connections.map(backend="mpl").savefig("conn.png")
        """

        refuse_noun_parameters("lak.connections", value_column, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.get(lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        choro = self.model.plot.map(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
            hover_spec=cell_input_hover(value_column),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale or "earth",
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class SfrStageResultsExplorer(StageResultsExplorer):
    """SFR stage explorer with reach-profile helpers."""

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind an SFR stage-result explorer using the reach stage table builder.

        ``package`` is the SFR package's NAME (``None``: the model's only one).
        """

        super().__init__(model, "sfr", build_sfr_stage_result_table, package=package)

    @property
    def profile(self) -> SfrReachProfileView:
        """Return the reach stage-profile view: ``.get()`` table, ``.plot()`` line, ``.summary()``."""

        return SfrReachProfileView(
            self, y_column="stage", y_label="stage", label="sfr.results.stage.profile"
        )

    def _legacy_profile_plot(self, **kwargs):
        """Back the retired ``plot_profile`` spelling; returns the same stage figure."""

        return self.profile.plot(**kwargs)

    __getattr__ = deprecated_instance_getattr(
        {
            "plot_profile": (
                "_legacy_profile_plot",
                "model.sfr.results.stage.profile.plot",
                "0.1.0",
            ),
        },
        "myflopy.modflow.mf6.package_surface_water.SfrStageResultsExplorer",
    )


class LakResultsNamespace(FieldMappable):
    """Namespace for LAK result explorers.

    Mappable fields: ``q`` (lake-groundwater exchange, the default) and
    ``stage`` -- use ``results.map(field="stage")`` or ``results.stage.map()``.
    ``stage_change`` stays a first-class accessor (``results.stage_change.plot()``)
    but is a per-transition series, not a spatial field, so it is not in
    ``field_names()``.
    """

    _default_field = "q"

    def _field_names(self):
        """The mappable LAK result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the LAK results namespace to ``model``'s LAK package ``package``.

        ``package`` is the package's NAME (``None``: the model's only LAK).
        """

        self.model = model
        self.package = package

    @property
    def stage(self) -> LakStageResultsExplorer:
        """Return the lake stage explorer mapped to connected cells."""

        return LakStageResultsExplorer(self.model, package=self.package)

    @property
    def stage_change(self) -> LakStageChangeExplorer:
        """Return the lake-stage change explorer."""

        return LakStageChangeExplorer(self.model, package=self.package)

    @property
    def q(self) -> LakBudgetResultsExplorer:
        """Return the lake-groundwater exchange result explorer.

        Exchange maps are normalized to lake connection area, so
        ``map(...)`` renders ``sum(q) / sum(flow_area)`` by cell.
        """

        budget_text, value_name = get_default_budget_term("lak") or ("GWF", "q_lake")
        return LakBudgetResultsExplorer(
            self.model, budget_text=budget_text, value_name=value_name, package=self.package
        )


class LakBudgetNamespace:
    """Namespace for all MF6-defined LAK package-output budget terms."""

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the LAK budget-term namespace to ``model``'s LAK package ``package``.

        ``package`` is the package's NAME (``None``: the model's only LAK).
        """

        self.model = model
        self.package = package

    @property
    def types(self) -> list[str]:
        """Return the LAK package-output budget term names available for the model."""

        return [
            str(value).strip().upper()
            for value in _budget_file(self.model, "lak", self.package).types
        ]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of LAK package-output budget terms.

        Every record of the LAK package's OWN budget file (not the GWF model
        budget), one row per entry per saved step -- a FILTER on ``per``, so every
        step a period saved stays; :meth:`summary` and :meth:`wide` read each
        period's last. The selectors combine by AND, and none raises on a value
        nothing matches: the table comes back empty.

        Parameters
        ----------
        term : str or iterable of str or None, optional
            Term(s) to keep, MODFLOW 6's spelling in any case -- ``"GWF"``,
            ``["GWF", "STORAGE"]``, ``"EXT-INFLOW"`` (:attr:`types` lists what the
            file holds). Compared after ``strip().upper()`` and nothing else, so
            the attribute spelling ``"ext_inflow"`` matches NOTHING and returns an
            empty table without a word (use ``model.lak.budget.ext_inflow``).
            ``None`` (default) keeps every term.
        per : int or iterable of int or None, optional
            Zero-based stress period(s) to keep, every saved step of each.
            ``None`` (default) keeps every period.
        lakes : int or iterable of int or None, optional
            Zero-based lake number(s) to keep (within this LAK package; the
            ``lake`` column, the record's ``node``). On ``FLOW-JA-FACE`` MODFLOW 6
            writes each lake-to-lake outlet TWICE, once from each end, with
            ``lake_to`` the other lake and ``q`` positive INTO ``lake`` -- so a
            lake selects both its outflows (negative) and its inflows
            (positive). ``None`` (default) keeps every lake.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``term``, ``kstpkper``, ``per``, ``lake``,
            ``q`` (model length³/time, in the sign the LAK budget file records:
            positive is an INFLOW to the lake's water budget -- on ``STORAGE``,
            water released as the lake's volume FALLS) and each term's own
            columns -- ``layer``/``cell``/``claktype``/``flow_area``/
            ``q_per_area`` on ``GWF``, ``lake_to`` on ``FLOW-JA-FACE``,
            ``VOLUME`` on ``STORAGE``.
        """

        frame = build_lak_budget_term_table(self.model, term=term, package=self.package)
        if per is not None:
            if isinstance(per, Iterable) and not isinstance(per, (str, bytes)):
                periods = {int(value) for value in per}
                frame = frame.loc[frame["per"].isin(periods)].copy()
            else:
                frame = frame.loc[frame["per"] == int(per)].copy()
        if lakes is not None and "lake" in frame.columns:
            if isinstance(lakes, Iterable) and not isinstance(lakes, (str, bytes)):
                lake_ids = {int(value) for value in lakes}
            else:
                lake_ids = {int(lakes)}
            frame = frame.loc[
                pd.to_numeric(frame["lake"], errors="coerce").isin(lake_ids)
            ].copy()
        return frame.reset_index(drop=True)

    def summary(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
    ) -> pd.DataFrame:
        """Summarize LAK budget terms by selected grouping columns.

        Each period is its LAST saved step on the LAK budget file's own clock
        (``period_end.at_last_saved_steps``, ledgers 337, 348): the records of that
        one step combine, never the steps of a period, as every picture of the
        period reads it. Every step the period saved was ADDED until the review of
        ledger 348 -- three steps' flows, and three copies of each connection's
        area, at ``nstp=3`` under ``SAVE BUDGET ALL``.

        Parameters
        ----------
        term : str or iterable of str or None, optional
            LAK budget term(s) to keep (``"GWF"``, ``["GWF", "STORAGE"]``, any
            case). ``None`` (default) keeps every term the file holds.
        per : int or iterable of int or None, optional
            Zero-based stress period(s) to keep. ``None`` (default) keeps every one.
        lakes : iterable of int or None, optional
            Zero-based lake number(s) to keep. ``None`` (default) keeps every one.
        by : list or tuple of str or None, optional
            Columns to group by. ``None`` (default) is ``["per", "lake", "term"]``.
            Naming ``"kstpkper"`` keeps every saved step, one group each, rather than
            each period's last.

        Returns
        -------
        pandas.DataFrame
            The ``by`` columns, ``record_count``, the summed ``q`` (model flow
            units, in the sign the package's budget file records) and the summed
            extra columns the term carries.
        """

        frame = self.get(term=term, per=per, lakes=lakes)
        group_columns = list(by) if by is not None else ["per", "lake", "term"]
        if "kstpkper" in frame.columns and "kstpkper" not in group_columns:
            frame = at_last_saved_steps(frame)
        if frame.empty:
            columns = [*group_columns, "record_count", "q"]
            return pd.DataFrame(columns=columns)

        agg_map: dict[str, tuple[str, str]] = {
            "record_count": ("q", "size"),
            "q": ("q", "sum"),
        }
        for column in ("FLOW-AREA", "flow_area", "VOLUME"):
            if column in frame.columns:
                agg_map[column] = (column, "sum")
        summary = (
            frame.groupby(group_columns, dropna=False, as_index=False)
            .agg(**agg_map)
            .sort_values(group_columns)
            .reset_index(drop=True)
        )
        if "flow_area" in summary.columns:
            summary["q_per_area"] = np.where(
                pd.to_numeric(summary["flow_area"], errors="coerce") > 0.0,
                pd.to_numeric(summary["q"], errors="coerce")
                / pd.to_numeric(summary["flow_area"], errors="coerce"),
                np.nan,
            )
        return summary

    def wide(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
    ) -> pd.DataFrame:
        """Pivot the long LAK budget table to one wide table by term.

        Each ``per`` row is the period's LAST saved step on the LAK budget file's own
        clock (``period_end.at_last_saved_steps``, ledgers 337, 348); every step the
        period saved was ADDED until the review of ledger 348.

        Parameters
        ----------
        term : str or iterable of str or None, optional
            LAK budget term(s) to keep as columns. ``None`` (default) keeps every
            term the file holds.
        per : int or iterable of int or None, optional
            Zero-based stress period(s) to keep. ``None`` (default) keeps every one.
        lakes : iterable of int or None, optional
            Zero-based lake number(s) to keep. ``None`` (default) keeps every one.
        index : list or tuple of str, default ("per", "lake")
            Columns labelling the rows. Naming ``"kstpkper"`` keeps every saved
            step, one row each, rather than each period's last.
        values : str, default "q"
            The column to pivot; records of one step on one row add.

        Returns
        -------
        pandas.DataFrame
            The ``index`` columns, then one column per term. Empty when nothing
            matches.

        Raises
        ------
        KeyError
            When ``values`` is not a column of the budget table.
        """

        frame = self.get(term=term, per=per, lakes=lakes)
        if "kstpkper" in frame.columns and "kstpkper" not in index:
            frame = at_last_saved_steps(frame)
        if frame.empty:
            return pd.DataFrame()
        if values not in frame.columns:
            raise KeyError(f"LAK budget column {values!r} was not found.")
        wide = frame.pivot_table(
            index=list(index),
            columns="term",
            values=values,
            aggfunc="sum",
        ).sort_index()
        if isinstance(wide.columns, pd.Index):
            wide.columns.name = None
        return wide.reset_index()

    @property
    def gwf(self) -> PackageBudgetTermExplorer:
        """Lake-groundwater exchange term helper."""

        return PackageBudgetTermExplorer(self, term="GWF", label="lak.budget.gwf")

    @property
    def storage(self) -> PackageBudgetTermExplorer:
        """Lake storage term helper."""

        return PackageBudgetTermExplorer(
            self, term="STORAGE", label="lak.budget.storage"
        )

    @property
    def runoff(self) -> PackageBudgetTermExplorer:
        """Lake runoff term helper."""

        return PackageBudgetTermExplorer(self, term="RUNOFF", label="lak.budget.runoff")

    @property
    def rainfall(self) -> PackageBudgetTermExplorer:
        """Lake rainfall term helper."""

        return PackageBudgetTermExplorer(
            self, term="RAINFALL", label="lak.budget.rainfall"
        )

    @property
    def evaporation(self) -> PackageBudgetTermExplorer:
        """Lake evaporation term helper."""

        return PackageBudgetTermExplorer(
            self, term="EVAPORATION", label="lak.budget.evaporation"
        )

    @property
    def withdrawal(self) -> PackageBudgetTermExplorer:
        """Lake withdrawal term helper."""

        return PackageBudgetTermExplorer(
            self, term="WITHDRAWAL", label="lak.budget.withdrawal"
        )

    @property
    def constant(self) -> PackageBudgetTermExplorer:
        """Lake constant-stage balancing flow term helper."""

        return PackageBudgetTermExplorer(
            self, term="CONSTANT", label="lak.budget.constant"
        )

    @property
    def ext_inflow(self) -> PackageBudgetTermExplorer:
        """External inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-INFLOW", label="lak.budget.ext_inflow"
        )

    @property
    def ext_outflow(self) -> PackageBudgetTermExplorer:
        """External outflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-OUTFLOW", label="lak.budget.ext_outflow"
        )

    @property
    def from_mvr(self) -> PackageBudgetTermExplorer:
        """Mover inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="FROM-MVR", label="lak.budget.from_mvr"
        )

    @property
    def to_mvr(self) -> PackageBudgetTermExplorer:
        """Mover outflow term helper."""

        return PackageBudgetTermExplorer(self, term="TO-MVR", label="lak.budget.to_mvr")

    @property
    def flow_ja_face(self) -> PackageBudgetTermExplorer:
        """Lake-to-lake outlet/routing connection term helper."""

        return PackageBudgetTermExplorer(
            self, term="FLOW-JA-FACE", label="lak.budget.flow_ja_face"
        )

    @property
    def auxiliary(self) -> PackageBudgetTermExplorer:
        """Auxiliary term helper."""

        return PackageBudgetTermExplorer(
            self, term="AUXILIARY", label="lak.budget.auxiliary"
        )

    @property
    def mvr(self) -> PackageBudgetTermExplorer:
        """Combined mover-related LAK budget term helper."""

        return PackageBudgetTermExplorer(
            self, term=["FROM-MVR", "TO-MVR"], label="lak.budget.mvr"
        )

    @property
    def lake_fluxes(self) -> PackageBudgetTermExplorer:
        """Combined lake-level flux term helper excluding connection-level GWF rows."""

        return PackageBudgetTermExplorer(
            self,
            term=[
                "EXT-INFLOW",
                "RUNOFF",
                "RAINFALL",
                "EVAPORATION",
                "WITHDRAWAL",
                "STORAGE",
                "CONSTANT",
                "EXT-OUTFLOW",
                "FROM-MVR",
                "TO-MVR",
            ],
            label="lak.budget.lake_fluxes",
        )


class SfrBudgetNamespace:
    """Namespace for all MF6-defined SFR package-output budget terms."""

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the SFR budget-term namespace to ``model``'s SFR package ``package``.

        ``package`` is the package's NAME (``None``: the model's only SFR).
        """

        self.model = model
        self.package = package

    @property
    def types(self) -> list[str]:
        """Return the SFR package-output budget term names available for the model."""

        return [
            str(value).strip().upper()
            for value in _budget_file(self.model, "sfr", self.package).types
        ]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of SFR package-output budget terms.

        Every record of the SFR package's OWN budget file (not the GWF model
        budget), one row per entry per saved step, each joined to the reach table
        -- a FILTER on ``per``, so every step a period saved stays; :meth:`summary`
        and :meth:`wide` read each period's last. The selectors combine by AND,
        and none raises on a value nothing matches: the table comes back empty.

        Parameters
        ----------
        term : str or iterable of str or None, optional
            Term(s) to keep, MODFLOW 6's spelling in any case -- ``"GWF"``,
            ``["GWF", "STORAGE"]``, ``"EXT-INFLOW"`` (:attr:`types` lists what the
            file holds). Compared after ``strip().upper()`` and nothing else, so
            the attribute spelling ``"ext_inflow"`` matches NOTHING and returns an
            empty table without a word (use ``model.sfr.budget.ext_inflow``).
            ``None`` (default) keeps every term.
        per : int or iterable of int or None, optional
            Zero-based stress period(s) to keep, every saved step of each.
            ``None`` (default) keeps every period.
        reaches : int or iterable of int or None, optional
            Zero-based reach number(s) to keep (the ``reach`` column, the record's
            ``node``). On ``FLOW-JA-FACE`` MODFLOW 6 writes each reach-to-reach
            connection TWICE, once from each end, with ``reach_to`` the other
            reach and ``q`` positive INTO ``reach`` -- measured at period 1 of the
            testing canonical model: reach 0 to its downstream reach 1 reads
            -21464.10 as ``(0, 1)`` and +21464.10 as ``(1, 0)`` -- so a reach
            selects both its outflows and its inflows.
            ``None`` (default) keeps every reach.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``term``, ``kstpkper``, ``per``, ``reach``,
            the reach's ``layer``/``cell`` (on ``GWF`` the connected
            groundwater cell's) and geometry (``rlen``, ``distance_*``), ``q``
            (model length³/time, in the sign the SFR budget file records:
            positive is an INFLOW to the reach's water budget -- the OPPOSITE of
            ``model.sfr.results.q``'s ``q_gwf`` on ``GWF``) and each term's own
            columns --
            ``q_per_length`` (``q / rlen``) on every term but ``FLOW-JA-FACE``,
            which carries ``reach_to``; ``FLOW-AREA`` and ``VOLUME`` where MODFLOW
            6 writes them.
        """

        frame = build_sfr_budget_term_table(self.model, term=term, package=self.package)
        if per is not None:
            if isinstance(per, Iterable) and not isinstance(per, (str, bytes)):
                periods = {int(value) for value in per}
                frame = frame.loc[frame["per"].isin(periods)].copy()
            else:
                frame = frame.loc[frame["per"] == int(per)].copy()
        if reaches is not None and "reach" in frame.columns:
            if isinstance(reaches, Iterable) and not isinstance(reaches, (str, bytes)):
                reach_ids = {int(value) for value in reaches}
            else:
                reach_ids = {int(reaches)}
            frame = frame.loc[
                pd.to_numeric(frame["reach"], errors="coerce").isin(reach_ids)
            ].copy()
        return frame.reset_index(drop=True)

    def summary(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
    ) -> pd.DataFrame:
        """Summarize SFR budget terms by selected grouping columns.

        Each period is its LAST saved step on the SFR budget file's own clock
        (``period_end.at_last_saved_steps``, ledgers 337, 348): the records of that
        one step combine, never the steps of a period, as every picture of the
        period reads it. Every step the period saved was ADDED until the review of
        ledger 348 -- three steps' flows at ``nstp=3`` under ``SAVE BUDGET ALL``.

        Parameters
        ----------
        term : str or iterable of str or None, optional
            SFR budget term(s) to keep (``"GWF"``, ``["GWF", "STORAGE"]``, any
            case). ``None`` (default) keeps every term the file holds.
        per : int or iterable of int or None, optional
            Zero-based stress period(s) to keep. ``None`` (default) keeps every one.
        reaches : iterable of int or None, optional
            Zero-based reach number(s) to keep. ``None`` (default) keeps every one.
        by : list or tuple of str or None, optional
            Columns to group by. ``None`` (default) is ``["per", "reach", "term"]``.
            Naming ``"kstpkper"`` keeps every saved step, one group each, rather than
            each period's last.

        Returns
        -------
        pandas.DataFrame
            The ``by`` columns, ``record_count``, the summed ``q`` (model flow
            units, in the sign the package's budget file records) and the summed
            extra columns the term carries.
        """

        frame = self.get(term=term, per=per, reaches=reaches)
        group_columns = list(by) if by is not None else ["per", "reach", "term"]
        if "kstpkper" in frame.columns and "kstpkper" not in group_columns:
            frame = at_last_saved_steps(frame)
        if frame.empty:
            columns = [*group_columns, "record_count", "q"]
            return pd.DataFrame(columns=columns)

        agg_map: dict[str, tuple[str, str]] = {
            "record_count": ("q", "size"),
            "q": ("q", "sum"),
        }
        if "rlen" in frame.columns and "reach" in group_columns:
            agg_map["rlen"] = ("rlen", "first")
        for column in ("FLOW-AREA", "VOLUME"):
            if column in frame.columns:
                agg_map[column] = (column, "sum")
        summary = (
            frame.groupby(group_columns, dropna=False, as_index=False)
            .agg(**agg_map)
            .sort_values(group_columns)
            .reset_index(drop=True)
        )
        if "rlen" in summary.columns:
            summary["q_per_length"] = np.where(
                pd.to_numeric(summary["rlen"], errors="coerce") > 0.0,
                pd.to_numeric(summary["q"], errors="coerce")
                / pd.to_numeric(summary["rlen"], errors="coerce"),
                np.nan,
            )
        return summary

    def wide(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "reach"),
        values: str = "q",
    ) -> pd.DataFrame:
        """Pivot the long SFR budget table to one wide table by term.

        Each ``per`` row is the period's LAST saved step on the SFR budget file's own
        clock (``period_end.at_last_saved_steps``, ledgers 337, 348); every step the
        period saved was ADDED until the review of ledger 348.

        Parameters
        ----------
        term : str or iterable of str or None, optional
            SFR budget term(s) to keep as columns. ``None`` (default) keeps every
            term the file holds.
        per : int or iterable of int or None, optional
            Zero-based stress period(s) to keep. ``None`` (default) keeps every one.
        reaches : iterable of int or None, optional
            Zero-based reach number(s) to keep. ``None`` (default) keeps every one.
        index : list or tuple of str, default ("per", "reach")
            Columns labelling the rows. Naming ``"kstpkper"`` keeps every saved
            step, one row each, rather than each period's last.
        values : str, default "q"
            The column to pivot; records of one step on one row add.

        Returns
        -------
        pandas.DataFrame
            The ``index`` columns, then one column per term. Empty when nothing
            matches.

        Raises
        ------
        KeyError
            When ``values`` is not a column of the budget table.
        """

        frame = self.get(term=term, per=per, reaches=reaches)
        if "kstpkper" in frame.columns and "kstpkper" not in index:
            frame = at_last_saved_steps(frame)
        if frame.empty:
            return pd.DataFrame()
        if values not in frame.columns:
            raise KeyError(f"SFR budget column {values!r} was not found.")
        wide = frame.pivot_table(
            index=list(index),
            columns="term",
            values=values,
            aggfunc="sum",
        ).sort_index()
        if isinstance(wide.columns, pd.Index):
            wide.columns.name = None
        return wide.reset_index()

    @property
    def gwf(self) -> PackageBudgetTermExplorer:
        """Stream-groundwater exchange term helper."""

        return PackageBudgetTermExplorer(self, term="GWF", label="sfr.budget.gwf")

    @property
    def flow_ja_face(self) -> PackageBudgetTermExplorer:
        """Reach-to-reach routing connection term helper."""

        return PackageBudgetTermExplorer(
            self, term="FLOW-JA-FACE", label="sfr.budget.flow_ja_face"
        )

    @property
    def ext_inflow(self) -> PackageBudgetTermExplorer:
        """External inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-INFLOW", label="sfr.budget.ext_inflow"
        )

    @property
    def runoff(self) -> PackageBudgetTermExplorer:
        """Runoff term helper."""

        return PackageBudgetTermExplorer(self, term="RUNOFF", label="sfr.budget.runoff")

    @property
    def rain(self) -> PackageBudgetTermExplorer:
        """Rainfall term helper."""

        return PackageBudgetTermExplorer(self, term="RAIN", label="sfr.budget.rain")

    @property
    def evaporation(self) -> PackageBudgetTermExplorer:
        """Evaporation term helper."""

        return PackageBudgetTermExplorer(
            self, term="EVAPORATION", label="sfr.budget.evaporation"
        )

    @property
    def ext_outflow(self) -> PackageBudgetTermExplorer:
        """External outflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-OUTFLOW", label="sfr.budget.ext_outflow"
        )

    @property
    def storage(self) -> PackageBudgetTermExplorer:
        """Storage term helper."""

        return PackageBudgetTermExplorer(
            self, term="STORAGE", label="sfr.budget.storage"
        )

    @property
    def from_mvr(self) -> PackageBudgetTermExplorer:
        """Mover inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="FROM-MVR", label="sfr.budget.from_mvr"
        )

    @property
    def to_mvr(self) -> PackageBudgetTermExplorer:
        """Mover outflow term helper."""

        return PackageBudgetTermExplorer(self, term="TO-MVR", label="sfr.budget.to_mvr")

    @property
    def auxiliary(self) -> PackageBudgetTermExplorer:
        """Auxiliary term helper."""

        return PackageBudgetTermExplorer(
            self, term="AUXILIARY", label="sfr.budget.auxiliary"
        )

    @property
    def mvr(self) -> PackageBudgetTermExplorer:
        """Combined mover-related SFR budget term helper."""

        return PackageBudgetTermExplorer(
            self, term=["FROM-MVR", "TO-MVR"], label="sfr.budget.mvr"
        )

    @property
    def stream_fluxes(self) -> PackageBudgetTermExplorer:
        """Combined reach-level flux term helper excluding GWF and routing rows."""

        return PackageBudgetTermExplorer(
            self,
            term=[
                "EXT-INFLOW",
                "RUNOFF",
                "RAIN",
                "EVAPORATION",
                "EXT-OUTFLOW",
                "STORAGE",
                "FROM-MVR",
                "TO-MVR",
            ],
            label="sfr.budget.stream_fluxes",
        )


def _exchange_frame(package: str) -> str:
    """The declared reference frame for a package's exchange ``q``.

    Reads ``ResultSpec.reference_frame`` from the registry so a map's colour
    orientation is derived from the same field that documents the sign, never a
    literal at the call site that could fall out of step with the data.
    """

    spec = get_package_result_spec(package, "q")
    return spec.reference_frame if spec is not None else "gwf"


def _signed_exchange_colors() -> tuple[str, str]:
    """Return the ``(gaining, losing)`` bar colors for signed exchange plots.

    Both are read off ``_blue_white_red_diverging_colorscale`` rather than
    written out, so discrete bars and continuous maps cannot drift apart. The
    colours are frame-independent (gaining is always blue, losing always red);
    which *sign* of ``q`` is gaining is a per-frame fact the caller applies. The
    signed SFR profile bars are in the "gwf" frame, where negative ``q`` is the
    reach gaining (BLUE) and positive ``q`` is the reach losing (RED).
    """

    scale = _blue_white_red_diverging_colorscale()
    return str(scale[0][1]), str(scale[-1][1])


class SfrProfileView:
    """The SFR long profile: one merged reach-ordered table and its figure.

    Follows the house view shape (``docs/view_layer_conventions.md``) -- a noun
    reached from the results namespace, carrying ``get`` for the frame, ``plot``
    for the figure, and ``summary`` for the compact digest. Calling the view
    rebinds the stress period, so these are the same figure::

        model.sfr.results.profile.plot(per=3)
        model.sfr.results.profile(per=3).plot()

    The table merges reach geometry, streambed elevations, simulated stage, and
    stream-groundwater exchange, so it feeds custom analysis as readily as it
    feeds ``plot``.
    """

    def __init__(self, model: SimulationBase, *, per: int = 0, package: str | None = None):
        """Bind the profile view to ``model`` at stress period ``per``.

        ``package`` is the SFR package's NAME (``None``: the model's only SFR).
        """

        self.model = model
        self.per = int(per)
        self.package = package

    def __repr__(self) -> str:
        """Show the bound period, since it is the view's only state."""

        return f"{type(self).__name__}(per={self.per})"

    def __call__(self, *, per: int) -> SfrProfileView:
        """Return an equivalent view bound to stress period ``per``.

        Parameters
        ----------
        per : int
            Zero-based stress period the new view reads when its verbs are called
            without one; required, keyword-only, passed through ``int()`` (so
            ``None`` raises ``TypeError``). This view is left unchanged.

        Returns
        -------
        SfrProfileView
            A view of the same SFR package bound to ``per``.
        """

        return type(self)(self.model, per=per, package=self.package)

    def _period(self, per: int | None) -> int:
        """Resolve an explicit ``per`` against the period bound to this view."""

        return self.per if per is None else int(per)

    def get(self, *, per: int | None = None) -> pd.DataFrame:
        """Return the merged reach-ordered profile table.

        One row per reach: geometry, ``packagedata`` fields, the stage and the
        stream-aquifer exchange, each read at the period's LAST saved step on its
        own file's clock (stage file, SFR budget file).

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period. ``None`` (default) is the period bound to the
            view (``profile(per=N)``, 0 unless rebound). A period the stage file,
            or an SFR budget file holding records, saved nothing in raises
            ``ValueError`` naming the periods it has.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``per``, ``reach``, ``layer``, ``cell``,
            ``rlen``, ``distance_start``/``distance_mid``/``distance_end``
            (running totals of ``rlen`` in REACH NUMBER order, model length units
            -- a downstream distance only where the reaches are numbered along one
            channel), the ``packagedata`` fields, ``streambed_top``/
            ``streambed_bottom``, ``stage`` and ``q_gwf`` (MF6's
            aquifer-referenced sign: negative where the reach gains) with
            ``q_per_length``. ``q_gwf`` is summed from
            ``model.sfr.results.q.get()``, which joins records to reaches on the
            CELL, so every reach of a cell holding several reads the CELL's total
            (testing canonical cell 195 at period 3: -493.56 on reaches 22, 37 and
            38, whose own flows are -20.91, -471.96 and -0.69; ledger 440).
        """

        return build_sfr_long_profile_table(self.model, per=self._period(per), package=self.package)

    def summary(self, *, per: int | None = None) -> pd.DataFrame:
        """Return a compact digest of the profile's stage and exchange fields.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period, read as :meth:`get` reads it. ``None``
            (default) is the period bound to the view.

        Returns
        -------
        pandas.DataFrame
            One row: ``label``, the record, period, layer and cell counts, and the
            value columns present among ``stage``, ``streambed_top``, ``q_gwf``,
            ``q_per_length``.
        """

        frame = self.get(per=per)
        return summarize_input_table(
            frame,
            label="sfr.results.profile",
            value_columns=[
                column
                for column in ("stage", "streambed_top", "q_gwf", "q_per_length")
                if column in frame.columns
            ],
        )

    def plot(
        self,
        *,
        per: int | None = None,
        x: str = "distance",
        include_stage: bool = True,
        include_streambed: bool = True,
        include_exchange: bool = True,
        signed_exchange: bool = True,
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot the long profile with the common hydrologic overlays.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period, read as :meth:`get` reads it (a period with
            no saved record raises ``ValueError``). ``None`` (default) is the
            period bound to the view.
        x : str, default "distance"
            The horizontal axis. ``"distance"`` (any case) is ``distance_mid``, the
            reach mid-point's running total of ``rlen`` in reach-number order
            (model length units); ANY other value is the zero-based reach number --
            ``"reach"`` is the spelling, but a misspelling draws by reach too,
            without a word.
        include_stage : bool, default True
            Draw the simulated stage (left axis, model length units). Drawn only
            when the table has a ``stage`` column.
        include_streambed : bool, default True
            Draw ``streambed_top`` and, dashed, ``streambed_bottom`` (left axis).
        include_exchange : bool, default True
            Draw ``q_gwf`` on a second, right-hand axis (model length³/time).
        signed_exchange : bool, default True
            How ``include_exchange`` draws it; no effect without it. ``True``: one
            bar per reach coloured by sign -- blue where the reach gains (NEGATIVE
            ``q_gwf``: SFR's cell record is aquifer-referenced, MF6's raw sign), red
            where it loses, matching the SFR map colorscale. ``False``: one line,
            uncoloured by sign.
        plot_fig : bool, default False
            Call ``show()`` on the figure before returning.
        return_fig : bool, default True
            Return the figure. With ``False`` the method returns ``None``, so
            ``plot_fig=False, return_fig=False`` builds a figure nobody sees.

        Returns
        -------
        myflopy.viz.Fig or None
            The long profile, or ``None`` when ``return_fig=False``.

        Raises
        ------
        ValueError
            When the stage or budget file saved nothing in the period.
        KeyError
            When the table lacks the chosen axis column.
        """

        period = self._period(per)
        frame = self.get(per=period)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Long-profile x-axis column {x_column!r} was not found.")

        customdata = np.column_stack([frame["reach"], frame["cell"]])
        if include_streambed and "streambed_top" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["streambed_top"],
                mode="lines",
                name="Streambed Top",
                line={"color": "#8c564b", "width": 2},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "streambed_top=%{y}<extra></extra>"
                ),
            )
        if include_streambed and "streambed_bottom" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["streambed_bottom"],
                mode="lines",
                name="Streambed Bottom",
                line={"color": "#c49c94", "width": 2, "dash": "dash"},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "streambed_bottom=%{y}<extra></extra>"
                ),
            )
        if include_stage and "stage" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["stage"],
                mode="lines+markers",
                name="Stage",
                line={"color": "#1f77b4", "width": 3},
                marker={"size": 7},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "stage=%{y}<extra></extra>"
                ),
            )
        if include_exchange and "q_gwf" in frame.columns:
            self._add_exchange_trace(
                fig,
                frame=frame,
                x_column=x_column,
                customdata=customdata,
                signed=signed_exchange,
            )
            fig.update_layout(
                yaxis2={
                    "title": "Exchange q",
                    "overlaying": "y",
                    "side": "right",
                    "showgrid": False,
                    "zeroline": True,
                }
            )

        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title="Elevation / Stage",
            title=f"SFR long profile (per={period})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None

    @staticmethod
    def _add_exchange_trace(fig, *, frame, x_column, customdata, signed: bool) -> None:
        """Add the exchange trace, as signed bars or as one unsigned line."""

        hovertemplate = (
            "reach=%{customdata[0]}<br>"
            "cell=%{customdata[1]}<br>"
            "q=%{y}<extra></extra>"
        )
        if not signed:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["q_gwf"],
                mode="lines+markers",
                name="Exchange q",
                line={"color": "#d62728", "width": 2},
                marker={"size": 6, "symbol": "diamond"},
                yaxis="y2",
                customdata=customdata,
                hovertemplate=hovertemplate,
            )
            return

        gaining_color, losing_color = _signed_exchange_colors()
        q = pd.to_numeric(frame["q_gwf"], errors="coerce").to_numpy(float)
        positions = pd.to_numeric(frame[x_column], errors="coerce").to_numpy(float)
        # One bar per reach, sized just under the reach spacing so adjacent
        # bars read as separate reaches. A single reach has no spacing to
        # measure, so fall back to plotly's own default width.
        spacing = np.diff(np.sort(positions[np.isfinite(positions)]))
        spacing = spacing[spacing > 0.0]
        width = float(np.median(spacing)) * 0.9 if spacing.size else None
        fig.add_bar(
            x=frame[x_column],
            y=frame["q_gwf"],
            name="Exchange q (blue gains, red loses)",
            marker={
                # SFR's q keeps MF6's raw sign (the "gwf" frame): NEGATIVE q is
                # the reach gaining -> blue; positive q is losing -> red.
                "color": np.where(q < 0.0, gaining_color, losing_color).tolist(),
                "line": {"width": 0},
            },
            opacity=0.55,
            width=width,
            yaxis="y2",
            customdata=customdata,
            hovertemplate=hovertemplate,
        )


class SfrResultsNamespace(FieldMappable):
    """Namespace for SFR result explorers.

    Fields: ``q`` (stream-groundwater exchange, the default) and ``stage``.
    Use ``results.map(field="stage")`` or ``results.stage.map()``.
    """

    _default_field = "q"

    def _field_names(self):
        """The mappable SFR result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the SFR results namespace to ``model``'s SFR package ``package``.

        ``package`` is the package's NAME (``None``: the model's only SFR).
        """

        self.model = model
        self.package = package

    @property
    def stage(self) -> StageResultsExplorer:
        """Return the stream stage explorer mapped to reach cells."""

        return SfrStageResultsExplorer(self.model, package=self.package)

    @property
    def q(self) -> SfrBudgetResultsExplorer:
        """Return the stream-groundwater exchange result explorer."""

        budget_text, value_name = get_default_budget_term("sfr") or ("SFR", "q_gwf")
        return SfrBudgetResultsExplorer(
            self.model, budget_text=budget_text, value_name=value_name, package=self.package
        )

    @property
    def profile(self) -> SfrProfileView:
        """Return the long-profile view: ``.get()`` for the table, ``.plot()``.

        Named ``profile`` to match ``results.stage.profile`` one level down --
        inside an SFR namespace a profile can only be longitudinal, so "long"
        carried no information. This namespace-level view merges every field;
        the field-level ones cover a single field.
        """

        return SfrProfileView(self.model, package=self.package)

    # -- backing methods for the pre-view spellings ---------------------------
    # ``profile`` replaced two older names when the derived tables became view
    # objects. These bodies preserve the OLD return values exactly (a DataFrame
    # and an unsigned-line figure); the mapping below is the only place the old
    # spellings appear, and they resolve solely through __getattr__ so they stay
    # out of dir()/completion (D12).
    #
    # Deliberately named ``_legacy_*`` rather than echoing the old spelling: a
    # private member called ``_long_profile_frame`` would still surface the
    # retired name in IDE completion, which is exactly what D12 exists to stop.
    def _legacy_profile_frame(self, *, per: int = 0) -> pd.DataFrame:
        """Back the retired frame spelling; returns what it always returned."""

        return self.profile.get(per=per)

    def _legacy_profile_plot(self, **kwargs):
        """Back the retired plot spelling.

        Pins ``signed_exchange=False`` so it keeps drawing the single unsigned
        line it always drew; signed bars are the new ``profile.plot()`` default.
        """

        kwargs.setdefault("signed_exchange", False)
        return self.profile.plot(**kwargs)

    __getattr__ = deprecated_instance_getattr(
        {
            "long_profile": (
                "_legacy_profile_frame",
                "model.sfr.results.profile.get",
                "0.1.0",
            ),
            "plot_long_profile": (
                "_legacy_profile_plot",
                "model.sfr.results.profile.plot",
                "0.1.0",
            ),
        },
        "myflopy.modflow.mf6.package_surface_water.SfrResultsNamespace",
    )


class SurfaceWaterExchangeResultsExplorer(SpatialView):
    """Combined SFR/LAK exchange explorer with one shared physical sign scale."""

    def __init__(self, model: SimulationBase):
        """Bind the combined SFR/LAK exchange explorer to ``model``."""

        self.model = model

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        include: str | Iterable[str] | None = None,
        lak_connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return combined SFR/LAK exchange rows in one unified sign convention.

        One row per period per groundwater cell per package type, each period read
        at its LAST saved step on its budget file's own clock; every SFR and LAK
        package of the model contributes (several of a type add in a cell).
        ``exchange_intensity`` is myflopy's unified sign -- positive is the
        surface-water feature GAINING from groundwater -- but not one unit:
        ``-q / rlen`` on an SFR row (length²/time, per unit of reach length) and
        ``q / flow_area`` on a LAK row (length/time, per unit of bed area).

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period to keep, one row per cell at its last saved
            step. One period (a list raises ``TypeError``); ``None`` (default)
            keeps every period.
        layer : int or iterable of int or None, optional
            Zero-based layer(s) of the groundwater cell to keep. ``None`` (default)
            keeps every layer.
        include : str or iterable of str or None, optional
            Package type(s) contributing rows, ``"sfr"`` and/or ``"lak"``, any
            case; anything else raises ``ValueError``. ``None`` (default) is both.
            A type the model has no package of adds no rows.
        lak_connection_type : str or iterable of str or None, optional
            LAK connection type(s) to keep -- ``"vertical"``, ``"horizontal"``,
            ``"embeddedh"``, ``"embeddedv"``, any case -- applied to the LAK rows
            only; SFR rows are untouched. ``None`` (default) keeps every type, and
            a type no connection has drops every LAK row without a word.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``per``, ``layer``, ``cell``, the summed ``q``
            (each package's raw sign) and ``exchange_intensity`` and ``source``
            (``"sfr"``/``"lak"``), plus ``rlen``/``q_per_length`` on SFR rows and
            ``flow_area``/``q_per_area`` on LAK rows (``NaN`` on the other's).
            The SFR rows sum ``model.sfr.results.q.get()``, which joins records to
            reaches on the CELL, so on a cell holding ``k`` reaches the SFR ``q``
            and ``rlen`` are each ``k`` times their true sums while
            ``exchange_intensity``, their ratio, is right (testing canonical cell
            195 at period 3: ``q`` -1480.68 for -493.56; ledger 440).
        """

        return build_surface_water_exchange_cell_table(
            self.model,
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of combined surface-water exchange rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label="surface_water.results.q",
            value_columns=["exchange_intensity"],
        )

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        per: int = 0,
        layer: int = 0,
        include: str | Iterable[str] | None = None,
        lak_connection_type: str | Iterable[str] | None = None,
        # -- table to cells --------------------------------------------------
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing --------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover -----------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """SFR and LAK exchange on ONE map, one scale, one sign convention.

        A RESULTS noun over both surface-water packages at once. The two write
        their budgets in OPPOSITE frames -- SFR's record is cell-side, LAK's is
        lake-side -- so drawing them on one map means normalizing first. This
        draws ``exchange_intensity``, myflopy's own unified field, in which:

        - **positive** = groundwater gaining INTO the surface-water feature
        - **negative** = the surface-water feature LOSING to groundwater

        and both packages are expressed per unit of contact (reach length for
        SFR, flow area for LAK), so a cell's colour reports intensity rather than
        how much stream or lakebed it holds. The two are NOT one unit: an SFR
        cell is length²/time (per length of reach), a LAK cell length/time (per
        area of bed), so they share a colourbar, not a dimension.

        Because the quantity is SIGNED, the scale is centred on zero:
        ``zmin``/``zmax`` default to ``-max|q|``/``+max|q|`` and the trace gets
        ``zmid=0``, which is what makes gaining and losing different colours
        rather than two shades of one. Pass ``zmin``/``zmax`` yourself to pin a
        scale across several periods -- that also disables the symmetric
        default, so pass both.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into the tooltip on every cell. Passed anyway, they ride
        ``**trace_kwargs`` to the free verb. ``show_layer_elevs`` and
        ``hover_surfaces`` ARE named (ledger 306): the layer frame is the grid's,
        so they act on any map.

        Parameters
        ----------
        per : int, default 0
            Stress period to draw, zero-based, at its LAST saved step on each
            budget file's own clock (ledger 337's rule; the steps were ADDED until
            ledger 348). One period per picture; use :meth:`animate` to flip
            through them.
        layer : int, default 0
            Zero-based layer. Connections in other layers are not drawn.
        include : str or iterable of str, optional
            Which packages contribute -- ``"sfr"``, ``"lak"``, or both. With
            none, every surface-water package the model has. Narrow it when one
            package's exchange is an order of magnitude larger and is flattening
            the other's colours.
        lak_connection_type : str or iterable of str, optional
            Restrict the LAK side to particular connection types
            (``"vertical"``, ``"horizontal"``, ``"embeddedh"``,
            ``"embeddedv"``). Has no effect on the SFR side, which has no such
            distinction.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, most often.
            Note the sign convention here is already myflopy's unified one, so
            there is rarely a reason to pass ``-1.0``. Applied to every cell after
            the per-cell normalization, ``fill_value`` included -- except when the
            selection holds no record at all, which draws the fill unscaled.
        fill_value : float, default 0.0
            Value given to cells with neither a reach nor a lake connection.
            ``float("nan")`` leaves them uncoloured, which usually reads better:
            surface water touches a small fraction of most grids.
        **trace_kwargs
            Every other keyword, forwarded to ``model.plot.map`` after two checks:
            ``values``, ``type``, ``custom_hover``, ``hover_heads`` and ``hover_ks``
            raise ``TypeError`` (this noun supplies its own field and hover), and
            ``hover_spec`` is read as the old spelling of ``hover``, which wins when
            both are given. A free-verb parameter this noun does not name -- ``zoom``,
            ``bgs``, ``kstpkper``, ``per_timestep``, ``hover_layers``,
            ``show_mounding``, ... -- reaches the verb by that name, and ``custom_zs``,
            which this noun computes, raises ``TypeError`` as a second value. The rest
            are properties of the Plotly ``Choroplethmap`` trace -- ``zmid``,
            ``colorbar``, ``reversescale``, ``showscale`` -- checked only when the
            figure is built, so a misspelled one raises ``ValueError`` at ``.fig`` or
            ``.show()``, and is IGNORED under ``backend="mpl"``, which never builds
            that trace. ``zmid`` defaults to ``0.0`` here, the centre of the
            symmetric scale; pass one to move it.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        get : the combined per-cell exchange table behind the picture.
        summary : the same records reduced to one row per package.
        myflopy.modflow.mf6.package_surface_water.SfrBudgetResultsExplorer.map :
            SFR alone, in MODFLOW 6's own cell-side sign convention.
        myflopy.modflow.mf6.package_surface_water.LakBudgetResultsExplorer.map :
            LAK alone, in MODFLOW 6's own lake-side sign convention.

        Examples
        --------
        >>> model.surface_water.results.q.map()
        >>> model.surface_water.results.q.map(per=5)
        >>> model.surface_water.results.q.map(include="sfr")
        >>> model.surface_water.results.q.map(lak_connection_type="vertical")
        >>> model.surface_water.results.q.map(fill_value=float("nan"))
        >>> model.surface_water.results.q.map(backend="mpl").savefig("sw.png")
        """

        refuse_noun_parameters(
            "surface_water.results.q", "exchange_intensity", trace_kwargs
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.get(
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )
        values, cell_hover = build_surface_water_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        if absmax > 0:
            zmin = -absmax if zmin is None else zmin
            zmax = absmax if zmax is None else zmax
        # `zmid` is a Plotly trace property, not a verb parameter, so it rides
        # the open tail the way `colorbar`/`reversescale` do.
        trace_kwargs.setdefault("zmid", 0.0)
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
            hover_spec=surface_water_hover(),
            zmin=zmin,
            zmax=zmax,
            # This map draws ``exchange_intensity``, myflopy's OWN unified field,
            # which is normalized so POSITIVE = the feature gains (see
            # build_surface_water_exchange_cell_table). That matches the
            # "feature" orientation -- blue at the positive end -- regardless of
            # each source package's own raw frame.
            colorscale=colorscale or _exchange_colorscale("feature"),
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class SurfaceWaterResultsNamespace(FieldMappable):
    """Namespace for combined surface-water result explorers (field: ``q``)."""

    _default_field = "q"

    def _field_names(self):
        """The single mappable combined surface-water field: exchange ``q``."""

        return ["q"]

    def __init__(self, model: SimulationBase):
        """Bind the combined surface-water results namespace to ``model``."""

        self.model = model

    @property
    def q(self) -> SurfaceWaterExchangeResultsExplorer:
        """Return one shared SFR/LAK exchange explorer."""

        return SurfaceWaterExchangeResultsExplorer(self.model)


class SurfaceWaterInputFieldExplorer(SpatialView):
    """One cell-mapped LAK or SFR input field."""

    def __init__(self, inputs: SurfaceWaterInputsNamespace, field_name: str):
        """Pin a LAK/SFR inputs namespace to one mappable numeric input field."""

        self.inputs = inputs
        self.model = inputs.model
        self.package_name = inputs.package_name
        self.field_name = str(field_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return rows containing this mapped input field.

        The inputs namespace's table (:meth:`SurfaceWaterInputsNamespace.get`)
        reduced to its identifying columns and this field.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period to keep. Every period carries every reach or
            connection, its period-block settings held until a later block changes
            them. One period (a list raises ``TypeError``); ``None`` (default)
            keeps all ``nper`` copies.
        layer : int or iterable of int or None, optional
            Zero-based layer(s) of the groundwater cell. ``None`` (default) keeps
            every layer.
        cells : int or iterable of int or None, optional
            Zero-based cell number(s) within a layer (the ``cell`` column), ANDed
            with ``layer``. ``None`` (default) keeps every cell.

        Returns
        -------
        pandas.DataFrame
            Whichever of ``model``, ``package``, ``per``, ``lake``, ``reach``,
            ``layer``, ``cell`` the table has, then the field's column.
        """

        frame = self.inputs.get(per=per, layer=layer, cells=cells)
        metadata = [
            column
            for column in ("model", "package", "per", "lake", "reach", "layer", "cell")
            if column in frame.columns
        ]
        return frame.loc[:, [*metadata, self.field_name]].copy()

    def summary(self) -> pd.DataFrame:
        """Return a compact one-field summary of this LAK/SFR input field."""

        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        per: int = 0,
        layer: int = 0,
        # -- table to cells --------------------------------------------------
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str | None = None,
        # -- colour ----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours --------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting ----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing --------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover -----------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """This LAK or SFR input field as a choropleth: one colour per cell.

        A RECORD noun over one declared field of a surface-water package -- a
        reach length, a lakebed leakance, an inflow. The package's table is
        reduced to one value per cell (``agg``), cells the package does not
        touch take ``fill_value``, and the result is drawn like any other
        per-cell array.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into the tooltip on every cell. Passed anyway, they ride
        ``**trace_kwargs`` to the free verb. ``show_layer_elevs`` and
        ``hover_surfaces`` ARE named (ledger 306): the layer frame is the grid's,
        so they act on any map.

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based. Every period carries every reach or
            connection, a period-block setting held until a later block changes
            it, so a field that does not vary with time reads the same in every
            period.
        layer : int, default 0
            Zero-based layer. Records in other layers are not drawn -- for a
            lake spanning several layers this is how you choose.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, most often.
            Applied AFTER ``agg`` (so under ``-1.0``, ``agg="max"`` draws
            minus the largest record) and to every cell, ``fill_value`` included:
            ``fill_value=1.0, multiplier=2.0`` draws 2.0 where there is no record
            -- but 1.0 when the selection holds no record at all.
        fill_value : float, default 0.0
            Value given to cells this package has no record for. Pass
            ``float("nan")`` to leave them uncoloured, which is usually clearer
            here: surface water touches a small fraction of most grids, so a 0.0
            backdrop takes over the colour scale.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, optional
            How several records landing in ONE cell are reduced to one number.
            The default is chosen per field rather than fixed: ``sum`` for the
            extensive ones (``connection_area``, ``rlen``, ``inflow``,
            ``runoff``), ``first`` for everything else, which is what keeps an
            elevation or a leakance from being added to itself. Set it
            explicitly when a cell holds several reaches and you want a
            per-reach statistic rather than the cell total.
        **trace_kwargs
            Every other keyword, forwarded to ``model.plot.map`` after two checks:
            ``values``, ``type``, ``custom_hover``, ``hover_heads`` and ``hover_ks``
            raise ``TypeError`` (this noun supplies its own field and hover), and
            ``hover_spec`` is read as the old spelling of ``hover``, which wins when
            both are given. A free-verb parameter this noun does not name -- ``zoom``,
            ``bgs``, ``kstpkper``, ``per_timestep``, ``hover_layers``,
            ``show_mounding``, ... -- reaches the verb by that name, and ``custom_zs``,
            which this noun computes, raises ``TypeError`` as a second value. The rest
            are properties of the Plotly ``Choroplethmap`` trace -- ``zmid``,
            ``colorbar``, ``reversescale``, ``showscale`` -- checked only when the
            figure is built, so a misspelled one raises ``ValueError`` at ``.fig`` or
            ``.show()``, and is IGNORED under ``backend="mpl"``, which never builds
            that trace.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        get : the records behind the picture, as a DataFrame.
        summary : the same records reduced to one row.

        Examples
        --------
        >>> model.sfr.inputs.rlen.map()
        >>> model.sfr.inputs.rhk.map(agg="mean", logscale=True)
        >>> model.lak.inputs.bedleak.map(fill_value=float("nan"))
        >>> model.sfr.inputs.inflow.map(per=3)
        >>> model.sfr.inputs.rwid.map(select="all_streams")
        >>> model.lak.inputs.strt.map(backend="mpl").savefig("strt.png")
        """

        refuse_noun_parameters(
            f"{self.package_name}.inputs.{self.field_name}",
            self.field_name,
            trace_kwargs,
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.get(per=per, layer=layer)
        if agg is None:
            agg = (
                "sum"
                if self.field_name in {"connection_area", "rlen", "inflow", "runoff"}
                else "first"
            )
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
            hover_spec=cell_input_hover(self.field_name),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale or "earth",
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class SurfaceWaterInputsNamespace(FieldMappable):
    """Consistent input exploration namespace for LAK and SFR.

    Every numeric input column is a first-class field node (``sfr.inputs.rhk``,
    ``lak.inputs.connection_area``) with the unified grammar; the namespace
    verbs take ``field=`` as sugar over them.
    """

    def __init__(self, model: SimulationBase, package_name: str, *, package: str | None = None):
        """Bind a LAK or SFR inputs namespace to ``model`` (``package_name``, the type, lowercased).

        ``package`` is the one package's NAME (``None``: the model's only package of
        the type).
        """

        self.model = model
        self.package_name = str(package_name).lower()
        self.package = package

    @property
    def _default_field(self) -> str:
        """Preferred input field: LAK ``connection_area``, SFR ``rhk``."""

        return "connection_area" if self.package_name == "lak" else "rhk"

    def _field_names(self) -> list[str]:
        """The mappable numeric input field names discovered for this package."""

        return self.fields["field"].tolist()

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return normalized package inputs mapped to groundwater cells.

        One row per reach (SFR) or per lake connection (LAK) per stress period,
        ``0 .. nper-1``: the static ``packagedata`` (and, for LAK, the
        ``connectiondata``) joined with the numeric period-block settings, each
        held from the period that sets it until a later one changes it, as MODFLOW
        6 holds them. A selector nothing matches gives an empty table, not an
        error.

        Parameters
        ----------
        per : int or None, optional
            Zero-based stress period to keep. One period (a list raises
            ``TypeError``); ``None`` (default) keeps all ``nper`` copies.
        layer : int or iterable of int or None, optional
            Zero-based layer(s) of the groundwater cell. ``None`` (default) keeps
            every layer.
        cells : int or iterable of int or None, optional
            Zero-based cell number(s) within a layer (the ``cell`` column), ANDed
            with ``layer``. ``None`` (default) keeps every cell.

        Returns
        -------
        pandas.DataFrame
            ``reach`` (SFR) or ``lake`` with the connection columns (LAK),
            ``layer``, ``cell``, every input field (:attr:`fields` lists the
            numeric ones), ``per``, ``package`` and ``model``.
        """

        builder = (
            build_lak_input_table
            if self.package_name == "lak"
            else build_sfr_input_table
        )
        return _filter_normalized_table(
            builder(self.model, package=self.package), per=per, layer=layer, cells=cells
        )

    @property
    def fields(self) -> pd.DataFrame:
        """Return numeric input fields that can be mapped."""

        frame = self.get()
        excluded = {
            "model",
            "package",
            "per",
            "lake",
            "reach",
            "layer",
            "cell",
            "ifno",
            "iconn",
            "ncon",
            "ndv",
            "nlakeconn",
        }
        fields = [
            column
            for column in frame.columns
            if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
        ]
        return pd.DataFrame({"field": fields})

    def summary(self) -> pd.DataFrame:
        """Return one stacked summary row per mappable input field for this package."""

        frames = [
            getattr(self, field).summary() for field in self.fields["field"].tolist()
        ]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def __getattr__(self, field_name: str) -> SurfaceWaterInputFieldExplorer:
        """Resolve ``inputs.<field>`` to a field-pinned explorer (else ``AttributeError``)."""

        if field_name not in set(self.fields["field"].tolist()):
            raise AttributeError(
                f"{type(self).__name__!s} has no input field {field_name!r}"
            )
        return SurfaceWaterInputFieldExplorer(self, field_name)

    @property
    def default(self) -> SurfaceWaterInputFieldExplorer:
        """Return the preferred package input field."""

        return getattr(self, self._default_field)

    # map/plot/xs/mosaic/animate come from FieldMappable (field= sugar).


class LakPackageExplorer:
    """Top-level LAK package explorer namespace: ONE LAK package's views.

    ``model.lak`` reads the model's only LAK package. A model carrying several --
    a second ``mf.lak`` beside a default-named first one gives ``lak`` and
    ``lak_2`` -- has every view refuse with ``SeveralPackagesError`` naming them,
    and ``model.lak(package="lak_2")`` is the one spelling that picks one: a lake
    is numbered within its package, so no view may choose (ledger 352).
    """

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the top-level LAK explorer (inputs/connections/budget/results) to ``model``.

        ``package`` is the LAK package's NAME; ``None``, the default, is the model's
        only one. Spell it ``model.lak(package=...)``, which checks the name.
        """

        self.model = model
        self.package = package

    def __call__(self, *, package: str) -> LakPackageExplorer:
        """The LAK views of the ONE LAK package named ``package``.

        Parameters
        ----------
        package : str
            A LAK package's name (FloPy's ``package_name``), any case. On a model
            holding ``lak`` and ``lak_2``, ``"lak"`` is the package NAMED ``lak`` --
            never the type. Checked here: no LAK package of that name raises
            ``ValueError`` naming the model's, and a non-string ``TypeError``.

        Returns
        -------
        LakPackageExplorer
            ``inputs``, ``connections``, ``budget`` and ``results``, each reading
            that package's lakes, connections and output files alone.

        Examples
        --------
        >>> model.lak(package="lak_2").results.stage.get()
        >>> model.lak(package="lak_2").results.q.map(per=1)
        >>> model.lak(package="LAK").budget.get(term="STORAGE")     # any case
        """

        return LakPackageExplorer(
            self.model,
            package=checked_package_name(self.model, "lak", package, spelling="model.lak(package=...)"),
        )

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped LAK starting-stage, connection, and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "lak", package=self.package)

    @property
    def connections(self) -> LakConnectionsExplorer:
        """Return LAK connection-geometry exploration helpers."""

        return LakConnectionsExplorer(self.model, package=self.package)

    @property
    def budget(self) -> LakBudgetNamespace:
        """Return LAK package-output budget helpers for all MF6-defined terms."""

        return LakBudgetNamespace(self.model, package=self.package)

    @property
    def results(self) -> LakResultsNamespace:
        """Return the LAK result exploration namespace."""

        return LakResultsNamespace(self.model, package=self.package)


class SfrPackageExplorer:
    """Top-level SFR package explorer namespace: ONE SFR package's views.

    ``model.sfr`` reads the model's only SFR package. On a model carrying several,
    every view refuses with ``SeveralPackagesError`` naming them, and
    ``model.sfr(package="sfr_2")`` is the one spelling that picks one: a reach is
    numbered within its package (ledger 352).
    """

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the top-level SFR explorer (inputs/budget/results) to ``model``.

        ``package`` is the SFR package's NAME; ``None``, the default, is the model's
        only one. Spell it ``model.sfr(package=...)``, which checks the name.
        """

        self.model = model
        self.package = package

    def __call__(self, *, package: str) -> SfrPackageExplorer:
        """The SFR views of the ONE SFR package named ``package``.

        Parameters
        ----------
        package : str
            An SFR package's name (FloPy's ``package_name``), any case -- never the
            type. Checked here: no SFR package of that name raises ``ValueError``
            naming the model's, and a non-string ``TypeError``.

        Returns
        -------
        SfrPackageExplorer
            ``inputs``, ``budget`` and ``results``, each reading that package's
            reaches and output files alone.

        Examples
        --------
        >>> model.sfr(package="sfr_2").results.stage.get()
        >>> model.sfr(package="sfr_2").results.profile.plot(per=2)
        """

        return SfrPackageExplorer(
            self.model,
            package=checked_package_name(self.model, "sfr", package, spelling="model.sfr(package=...)"),
        )

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped SFR reach-hydraulic and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "sfr", package=self.package)

    @property
    def budget(self) -> SfrBudgetNamespace:
        """Return SFR package-output budget helpers for all MF6-defined terms."""

        return SfrBudgetNamespace(self.model, package=self.package)

    @property
    def results(self) -> SfrResultsNamespace:
        """Return the SFR result exploration namespace."""

        return SfrResultsNamespace(self.model, package=self.package)


class SurfaceWaterPackageExplorer:
    """Top-level combined surface-water explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level combined SFR/LAK explorer to ``model``."""

        self.model = model

    @property
    def results(self) -> SurfaceWaterResultsNamespace:
        """Return combined SFR/LAK result explorers."""

        return SurfaceWaterResultsNamespace(self.model)


__all__ = [
    "SfrBudgetResultsExplorer",
    "LakBudgetResultsExplorer",
    "LakStageResultsExplorer",
    "LakStageChangeExplorer",
    "LakConnectionsExplorer",
    "SfrStageResultsExplorer",
    "SfrProfileView",
    "SfrReachProfileView",
    "LakBudgetView",
    "LakResultsNamespace",
    "LakBudgetNamespace",
    "SfrBudgetNamespace",
    "SfrResultsNamespace",
    "SurfaceWaterExchangeResultsExplorer",
    "SurfaceWaterResultsNamespace",
    "SurfaceWaterInputFieldExplorer",
    "SurfaceWaterInputsNamespace",
    "LakPackageExplorer",
    "SfrPackageExplorer",
    "SurfaceWaterPackageExplorer",
]
