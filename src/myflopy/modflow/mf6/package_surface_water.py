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
    get_default_budget_term,
    get_package_result_spec,
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
from myflopy.modflow.utils.datatypes.hover import (
    cell_input_hover,
    lak_hover,
    sfr_hover,
    surface_water_hover,
)

logger = get_logger(__name__)

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


def _join_feature_stage(frame, stage_table, *, per=None):
    """Join per-cell feature stage onto an exchange frame for the hover.

    ``stage_table`` is a normalized stage result table (one row per feature-cell
    per period). Cells touched by more than one feature take the mean stage,
    matching :class:`StageResultsExplorer`'s series aggregate.
    """

    if frame.empty or "stage" in frame.columns or stage_table.empty:
        return frame
    selected = stage_table
    if per is not None and "per" in selected.columns:
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


def join_sfr_stage(model, frame, *, per=None):
    """Add SFR reach stage to an exchange frame; a no-op if stage is unavailable."""

    try:
        stage_table = build_sfr_stage_result_table(model)
    except _STAGE_UNAVAILABLE:
        logger.debug(
            "no SFR stage for %s; the exchange hover will omit it",
            getattr(model, "name", model), exc_info=True,
        )
        return frame
    return _join_feature_stage(frame, stage_table, per=per)


def join_lak_stage(model, frame, *, per=None):
    """Add lake stage to an exchange frame; a no-op if stage is unavailable."""

    try:
        stage_table = build_lak_stage_result_table(model)
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

        model.packages.sfr.results.q.profile.plot(per=3)
        model.packages.sfr.results.q.profile(per=3).plot()
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
        """Return an equivalent view bound to stress period ``per``."""

        return type(self)(
            self.explorer, y_column=self.y_column, y_label=self.y_label, label=self.label, per=per
        )

    def _period(self, per: int | None) -> int:
        """Resolve an explicit ``per`` against the period bound to this view."""

        return self.per if per is None else int(per)

    def get(self, *, per: int | None = None) -> pd.DataFrame:
        """Return the reach-ordered profile table for the selected period."""

        frame = self.explorer.get(per=self._period(per))
        if frame.empty:
            return frame
        return frame.sort_values(["reach", "cell"]).reset_index(drop=True)

    def summary(self, *, per: int | None = None) -> pd.DataFrame:
        """Return a compact digest of the profiled field."""

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
        per
            Zero-based stress period; defaults to the period bound to the view.
        x
            Either ``"distance"`` for cumulative stream distance or ``"reach"``
            for raw reach number.
        plot_fig
            If ``True``, call ``show()`` on the created figure.
        return_fig
            If ``True``, return the created figure.
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
    ):
        """Bind an SFR budget-result explorer (defaults to the ``SFR`` term's ``q_gwf``).

        ``value_name`` is the emitted column (``q_gwf`` -- the SFR exchange is
        aquifer-referenced); the accessor stays ``results.q`` via ``result_name``.
        """

        super().__init__(model, "sfr", budget_text, value_name, result_name="q")

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized SFR budget-result table for selected rows."""

        frame = build_sfr_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
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
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
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
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        the tooltip on every cell.

        Parameters
        ----------
        per : int, default 0
            Stress period to draw, zero-based. One period per picture; use
            :meth:`animate` to flip through them.
        layer : int, default 0
            Zero-based layer. Reaches connected in other layers are not drawn.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, or ``-1.0`` to
            read gaining as positive. That changes the PICTURE only: ``.get()``
            keeps MF6's convention, so the colorbar will then disagree with the
            table, and the caption should say so.
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
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

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
        >>> model.packages.sfr.results.q.map()
        >>> model.packages.sfr.results.q.map(per=5)
        >>> model.packages.sfr.results.q.map(zmin=-0.5, zmax=0.5)
        >>> model.packages.sfr.results.q.map(fill_value=float("nan"))
        >>> model.packages.sfr.results.q.map(select="all_streams", select_style="both")
        >>> model.packages.sfr.results.q.map(backend="mpl").savefig("sfr_q.png")
        """

        del agg
        refuse_noun_parameters(self.label, self.value_name, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = join_sfr_stage(self.model, self.get(per=per, layer=layer), per=per)
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
            hover_spec=sfr_hover(),
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
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
                "model.packages.sfr.results.q.profile.plot",
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
        """Return an equivalent view bound to stress period ``per``."""

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
        """Return the exchange summary grouped by period, lake, and connection type."""

        return self.explorer._build_budget_summary(
            per=self._period(per), connection_type=connection_type
        )

    def summary(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return a compact digest of the budget summary."""

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
        per
            Zero-based stress period; defaults to the period bound to the view.
        value
            Summary field to visualize: ``"q"``, ``"flow_area"``, or ``"q_per_area"``.
        ax
            Optional Matplotlib axes to draw onto.
        return_fig
            If ``True``, return the created figure.
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
    ):
        """Bind a LAK budget-result explorer (defaults to the ``GWF`` exchange ``q_lake``).

        ``value_name`` is the emitted column (``q_lake`` -- LAK's exchange is
        feature-referenced); the accessor stays ``results.q`` via ``result_name``.
        """

        super().__init__(model, "lak", budget_text, value_name, result_name="q")

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK budget-result table for selected rows."""

        frame = build_lak_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
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
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
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
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        the tooltip on every cell.

        Parameters
        ----------
        per : int, default 0
            Stress period to draw, zero-based. One period per picture; use
            :meth:`animate` to flip through them.
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
        fill_value : float, default 0.0
            Value given to cells with no lake connection. Zero is honest for an
            exchange, but ``float("nan")`` leaves them uncoloured, which reads
            better when the lakes cover a small part of the grid.
        agg : str, default 'sum'
            **Accepted and not used.** The normalization is computed explicitly
            from total exchange over total flow area per cell, so no free choice
            of reducer remains. Kept in the signature because removing it would
            break callers that pass it; it is a no-op, not a silent alternative.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

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
        >>> model.packages.lak.results.q.map()
        >>> model.packages.lak.results.q.map(per=5, connection_type="vertical")
        >>> model.packages.lak.results.q.map(layer=1, fill_value=float("nan"))
        >>> model.packages.lak.results.q.map(zmin=-0.2, zmax=0.2)
        >>> model.packages.lak.results.q.map(contours=True, contour_levels=6)
        >>> model.packages.lak.results.q.map(backend="mpl").savefig("lak_q.png")
        """

        del agg
        refuse_noun_parameters(self.label, self.value_name, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = join_lak_stage(
            self.model,
            self.get(per=per, layer=layer, connection_type=connection_type),
            per=per,
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
            hover_spec=lak_hover(),
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
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

        Parameters
        ----------
        per
            Optional zero-based stress period filter. When omitted, all periods
            are summarized.

        Returns
        -------
        pandas.DataFrame
            Summary table with signed volumetric exchange ``q`` and
            area-normalized exchange ``q_per_area`` grouped by period, lake,
            and connection type.
        """

        frame = self.get(per=per, connection_type=connection_type)
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
                "model.packages.lak.results.q.budget.get",
                "0.1.0",
            ),
            "plot_budget": (
                "_legacy_budget_plot",
                "model.packages.lak.results.q.budget.plot",
                "0.1.0",
            ),
        },
        "myflopy.modflow.mf6.package_surface_water.LakBudgetResultsExplorer",
    )


class LakStageResultsExplorer(StageResultsExplorer):
    """LAK stage explorer -- the unified grammar's ``plot()`` draws one line
    per lake by stress period (replaced the old ``plot_timeseries``)."""

    def __init__(self, model: SimulationBase):
        """Bind a LAK stage-result explorer using the lake stage table builder."""

        super().__init__(model, "lak", build_lak_stage_result_table)


class LakStageChangeExplorer:
    """Explorer for lake-stage changes between stress periods."""

    def __init__(self, model: SimulationBase):
        """Bind the lake stage-change explorer to ``model``."""

        self.model = model

    def get(
        self,
        *,
        lake: int | None = None,
        per0: int | None = None,
        per1: int | None = None,
    ) -> pd.DataFrame:
        """Return stage-change rows for one lake and/or one period transition."""

        frame = build_lak_stage_change_table(self.model)
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

    def __init__(self, model: SimulationBase):
        """Bind the LAK connection-geometry explorer to ``model``."""

        self.model = model

    def get(
        self,
        *,
        lake: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK connection table."""

        frame = build_lak_connection_table(self.model)
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
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
        # -- renderer --------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """LAK connection GEOMETRY as a choropleth: one colour per cell.

        An INPUTS noun. It draws how a lake is wired to the aquifer -- the area
        or width of each connection -- not what flows through it. For the flow,
        see ``model.packages.lak.results.q.map()``.

        Connection geometry has no time axis, so there is no ``per``: a lake's
        connections are declared once and do not change with the stress period.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        the tooltip on every cell.

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
            Applied before ``agg``.
        fill_value : float, default 0.0
            Value given to cells with no lake connection. Pass ``float("nan")``
            to leave them uncoloured, which is usually clearer here: a lake
            covers a small part of the grid, and a 0.0 backdrop takes over the
            colour scale.
        per : int, optional
            Stress period to read (0-based). Mutually exclusive with ``kstpkper``;
            with neither, the model's first output time is used.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

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
        >>> model.packages.lak.connections.map()
        >>> model.packages.lak.connections.map(lake=0, layer=1)
        >>> model.packages.lak.connections.map(value_column="connwidth", agg="mean")
        >>> model.packages.lak.connections.map(fill_value=float("nan"))
        >>> model.packages.lak.connections.map(logscale=True)
        >>> model.packages.lak.connections.map(backend="mpl").savefig("conn.png")
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class SfrStageResultsExplorer(StageResultsExplorer):
    """SFR stage explorer with reach-profile helpers."""

    def __init__(self, model: SimulationBase):
        """Bind an SFR stage-result explorer using the reach stage table builder."""

        super().__init__(model, "sfr", build_sfr_stage_result_table)

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
                "model.packages.sfr.results.stage.profile.plot",
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

    def __init__(self, model: SimulationBase):
        """Bind the LAK results namespace to ``model``."""

        self.model = model

    @property
    def stage(self) -> LakStageResultsExplorer:
        """Return the lake stage explorer mapped to connected cells."""

        return LakStageResultsExplorer(self.model)

    @property
    def stage_change(self) -> LakStageChangeExplorer:
        """Return the lake-stage change explorer."""

        return LakStageChangeExplorer(self.model)

    @property
    def q(self) -> LakBudgetResultsExplorer:
        """Return the lake-groundwater exchange result explorer.

        Exchange maps are normalized to lake connection area, so
        ``map(...)`` renders ``sum(q) / sum(flow_area)`` by cell.
        """

        budget_text, value_name = get_default_budget_term("lak") or ("GWF", "q_lake")
        return LakBudgetResultsExplorer(
            self.model, budget_text=budget_text, value_name=value_name
        )


class LakBudgetNamespace:
    """Namespace for all MF6-defined LAK package-output budget terms."""

    def __init__(self, model: SimulationBase):
        """Bind the LAK budget-term namespace to ``model``."""

        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the LAK package-output budget term names available for the model."""

        return [
            str(value).strip().upper() for value in self.model.outputs.lak.bud.types
        ]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of LAK package-output budget terms.

        Parameters
        ----------
        term
            Optional LAK budget term filter such as ``"GWF"`` or
            ``["GWF", "STORAGE"]``.
        per
            Optional zero-based stress period filter.
        lakes
            Optional iterable of zero-based lake ids to keep.
        """

        frame = build_lak_budget_term_table(self.model, term=term)
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
        """Summarize LAK budget terms by selected grouping columns."""

        frame = self.get(term=term, per=per, lakes=lakes)
        group_columns = list(by) if by is not None else ["per", "lake", "term"]
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
        """Pivot the long LAK budget table to one wide table by term."""

        frame = self.get(term=term, per=per, lakes=lakes)
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

    def __init__(self, model: SimulationBase):
        """Bind the SFR budget-term namespace to ``model``."""

        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the SFR package-output budget term names available for the model."""

        return [
            str(value).strip().upper() for value in self.model.outputs.sfr.bud.types
        ]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of SFR package-output budget terms."""

        frame = build_sfr_budget_term_table(self.model, term=term)
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
        """Summarize SFR budget terms by selected grouping columns."""

        frame = self.get(term=term, per=per, reaches=reaches)
        group_columns = list(by) if by is not None else ["per", "reach", "term"]
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
        """Pivot the long SFR budget table to one wide table by term."""

        frame = self.get(term=term, per=per, reaches=reaches)
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

        model.packages.sfr.results.profile.plot(per=3)
        model.packages.sfr.results.profile(per=3).plot()

    The table merges reach geometry, streambed elevations, simulated stage, and
    stream-groundwater exchange, so it feeds custom analysis as readily as it
    feeds ``plot``.
    """

    def __init__(self, model: SimulationBase, *, per: int = 0):
        """Bind the profile view to ``model`` at stress period ``per``."""

        self.model = model
        self.per = int(per)

    def __repr__(self) -> str:
        """Show the bound period, since it is the view's only state."""

        return f"{type(self).__name__}(per={self.per})"

    def __call__(self, *, per: int) -> SfrProfileView:
        """Return an equivalent view bound to stress period ``per``."""

        return type(self)(self.model, per=per)

    def _period(self, per: int | None) -> int:
        """Resolve an explicit ``per`` against the period bound to this view."""

        return self.per if per is None else int(per)

    def get(self, *, per: int | None = None) -> pd.DataFrame:
        """Return the merged reach-ordered profile table."""

        return build_sfr_long_profile_table(self.model, per=self._period(per))

    def summary(self, *, per: int | None = None) -> pd.DataFrame:
        """Return a compact digest of the profile's stage and exchange fields."""

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
        per
            Zero-based stress period; defaults to the period bound to the view.
        x
            Either ``"distance"`` for cumulative stream distance or ``"reach"``
            for raw reach number.
        include_stage
            Whether to show simulated stream stage.
        include_streambed
            Whether to show streambed top and bottom elevations.
        include_exchange
            Whether to show stream-groundwater exchange on a secondary axis.
        signed_exchange
            When ``True`` (the default) draw exchange as per-reach bars colored
            by sign -- blue where the reach gains (NEGATIVE q_gwf: SFR's cell
            record is aquifer-referenced, MF6's raw sign), red where it loses,
            matching the SFR map colorscale. When ``False`` draw one unsigned line.
        plot_fig
            If ``True``, call ``show()`` on the created figure.
        return_fig
            If ``True``, return the created figure.
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

    def __init__(self, model: SimulationBase):
        """Bind the SFR results namespace to ``model``."""

        self.model = model

    @property
    def stage(self) -> StageResultsExplorer:
        """Return the stream stage explorer mapped to reach cells."""

        return SfrStageResultsExplorer(self.model)

    @property
    def q(self) -> SfrBudgetResultsExplorer:
        """Return the stream-groundwater exchange result explorer."""

        budget_text, value_name = get_default_budget_term("sfr") or ("SFR", "q_gwf")
        return SfrBudgetResultsExplorer(
            self.model, budget_text=budget_text, value_name=value_name
        )

    @property
    def profile(self) -> SfrProfileView:
        """Return the long-profile view: ``.get()`` for the table, ``.plot()``.

        Named ``profile`` to match ``results.stage.profile`` one level down --
        inside an SFR namespace a profile can only be longitudinal, so "long"
        carried no information. This namespace-level view merges every field;
        the field-level ones cover a single field.
        """

        return SfrProfileView(self.model)

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
                "model.packages.sfr.results.profile.get",
                "0.1.0",
            ),
            "plot_long_profile": (
                "_legacy_profile_plot",
                "model.packages.sfr.results.profile.plot",
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
        """Return combined SFR/LAK exchange rows in a unified L/T convention."""

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
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
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
        SFR, flow area for LAK) so a stream and a lakebed are comparable numbers
        rather than two different quantities sharing a colourbar.

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
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        the tooltip on every cell.

        Parameters
        ----------
        per : int, default 0
            Stress period to draw, zero-based. One period per picture; use
            :meth:`animate` to flip through them.
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
            there is rarely a reason to pass ``-1.0``.
        fill_value : float, default 0.0
            Value given to cells with neither a reach nor a lake connection.
            ``float("nan")`` leaves them uncoloured, which usually reads better:
            surface water touches a small fraction of most grids.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

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
        >>> model.packages.surface_water.results.q.map()
        >>> model.packages.surface_water.results.q.map(per=5)
        >>> model.packages.surface_water.results.q.map(include="sfr")
        >>> model.packages.surface_water.results.q.map(lak_connection_type="vertical")
        >>> model.packages.surface_water.results.q.map(fill_value=float("nan"))
        >>> model.packages.surface_water.results.q.map(backend="mpl").savefig("sw.png")
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
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
        """Return rows containing this mapped input field."""

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
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover -----------------------------------------------------------
        hover=None,
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
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        the tooltip on every cell.

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based. A field that does not vary with
            time has all its records under period 0.
        layer : int, default 0
            Zero-based layer. Records in other layers are not drawn -- for a
            lake spanning several layers this is how you choose.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, most often.
            Applied before ``agg``.
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
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
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
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
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
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

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
        >>> model.packages.sfr.inputs.rlen.map()
        >>> model.packages.sfr.inputs.rhk.map(agg="mean", logscale=True)
        >>> model.packages.lak.inputs.bedleak.map(fill_value=float("nan"))
        >>> model.packages.sfr.inputs.inflow.map(per=3)
        >>> model.packages.sfr.inputs.rwid.map(select="all_streams")
        >>> model.packages.lak.inputs.strt.map(backend="mpl").savefig("strt.png")
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
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class SurfaceWaterInputsNamespace(FieldMappable):
    """Consistent input exploration namespace for LAK and SFR.

    Every numeric input column is a first-class field node (``sfr.inputs.rhk``,
    ``lak.inputs.connection_area``) with the unified grammar; the namespace
    verbs take ``field=`` as sugar over them.
    """

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a LAK or SFR inputs namespace to ``model`` (``package_name`` lowercased)."""

        self.model = model
        self.package_name = str(package_name).lower()

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
        """Return normalized package inputs mapped to groundwater cells."""

        builder = (
            build_lak_input_table
            if self.package_name == "lak"
            else build_sfr_input_table
        )
        return _filter_normalized_table(
            builder(self.model), per=per, layer=layer, cells=cells
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
    """Top-level LAK package explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level LAK explorer (inputs/connections/budget/results) to ``model``."""

        self.model = model

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped LAK starting-stage, connection, and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "lak")

    @property
    def connections(self) -> LakConnectionsExplorer:
        """Return LAK connection-geometry exploration helpers."""

        return LakConnectionsExplorer(self.model)

    @property
    def budget(self) -> LakBudgetNamespace:
        """Return LAK package-output budget helpers for all MF6-defined terms."""

        return LakBudgetNamespace(self.model)

    @property
    def results(self) -> LakResultsNamespace:
        """Return the LAK result exploration namespace."""

        return LakResultsNamespace(self.model)


class SfrPackageExplorer:
    """Top-level SFR package explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level SFR explorer (inputs/budget/results) to ``model``."""

        self.model = model

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped SFR reach-hydraulic and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "sfr")

    @property
    def budget(self) -> SfrBudgetNamespace:
        """Return SFR package-output budget helpers for all MF6-defined terms."""

        return SfrBudgetNamespace(self.model)

    @property
    def results(self) -> SfrResultsNamespace:
        """Return the SFR result exploration namespace."""

        return SfrResultsNamespace(self.model)


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
