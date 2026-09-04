"""Generic result explorer classes behind model.packages."""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import TYPE_CHECKING

import pandas as pd

from myflopy._logging import get_logger

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_budget import (
    budget_value_units,
    build_budget_result_table,
)
from myflopy.modflow.mf6.package_explorer_utils import (
    _filter_normalized_table,
    _normalize_iterable_filter,
    _normalize_term_filter,
)
from myflopy.modflow.mf6.package_plotting import (
    FieldMappable,
    SpatialView,
    _apply_backend,
    _symmetric_color_limit,
    build_cell_input_map_payload,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_registry import (
    get_default_package_colorscale,
    get_package_explorer_spec,
    get_package_result_spec,
)
from myflopy.modflow.mf6.package_tables import (
    summarize_input_table,
)
from myflopy.modflow.utils.datatypes.hover import result_hover

logger = get_logger(__name__)


class CellBudgetResultsExplorer(SpatialView):
    """Normalized explorer for one cell-based package result term."""

    def __init__(
        self,
        model: SimulationBase,
        package_name: str,
        budget_text: str,
        value_name: str,
        result_name: str | None = None,
        label: str | None = None,
    ):
        """Bind one cell-based result term: its MF6 ``budget_text`` and output column.

        ``result_name`` is the registry/accessor identity (e.g. ``"q"``, the name
        behind ``results.q``); ``value_name`` is the emitted DataFrame COLUMN
        (e.g. ``"q_gwf"``), which names its reference frame. The two differ for
        the signed exchange results and are equal for everything else.
        ``package_name`` is lowercased.

        ``label`` overrides how :meth:`summary` names this view. It defaults to
        the package-grammar spelling, which is right for ``results.<field>`` but
        wrong for a model-budget term reached through ``model.budget.<term>`` --
        that is not a package result and should not claim to be one.
        """

        self.model = model
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)
        self.result_name = str(result_name) if result_name is not None else str(value_name)
        self.label = (
            str(label)
            if label is not None
            else f"{self.package_name}.results.{self.value_name}"
        )

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized result table for this budget term."""

        frame = build_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            package_name=self.package_name,
            value_name=self.value_name,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the available result rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=self.label,
            value_columns=[self.value_name],
        )

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.DataFrame:
        """Pivot this result term to one column per stress period."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        if frame.empty:
            return pd.DataFrame(columns=[*index])
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Wide result index columns were not found: {missing_index}")

        wide = frame.pivot_table(
            index=list(index),
            columns="per",
            values=value_column,
            aggfunc=agg,
        )
        wide.columns = [f"per_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def long(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.Series:
        """Return a long result series indexed by ``kstpkper/layer/cell``."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        index_columns = ["kstpkper", "layer", "cell"]
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays(
                [[] for _ in index_columns],
                names=index_columns,
            )
            return pd.Series([], index=empty_index, dtype=float, name=value_column)
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [
            column for column in index_columns if column not in frame.columns
        ]
        if missing_index:
            raise KeyError(f"Long result index columns were not found: {missing_index}")

        series = (
            frame.groupby(index_columns, dropna=False)[value_column]
            .agg(agg)
            .sort_index()
        )
        series.name = value_column
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`."""

        return self.long(**kwargs)

    # NOTE: the series view of this result is the unified grammar's ``plot()``
    # (SpatialView) -- ``results.q.plot(cells=[...])`` replaced plot_timeseries.

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
        """This budget term as a choropleth: one colour per cell.

        A RESULTS noun over a cell-based budget term. The simulated flows for one
        stress period are reduced to one value per cell (``agg``), cells this
        term never touches take ``fill_value``, and the result is drawn like any
        other per-cell array. That reduction is why this signature carries
        ``multiplier``/``fill_value``/``agg`` on top of the shared drawing
        parameters -- the free verb never does it, because it is handed the
        values directly.

        **The signs are MODFLOW 6's own, and they are not negated here.** A
        positive ``q`` is flow INTO the groundwater cell from the boundary and a
        negative one is flow out of it, exactly as the budget file records it, so
        a number read off the tooltip matches the number in ``get()``. Because
        the quantity is signed, the exchange term ``q`` is drawn on a diverging
        colorscale centred on zero: ``zmin``/``zmax`` default to
        ``-max|q|``/``+max|q|`` and the trace gets ``zmid=0``, which is what
        makes "gaining" and "losing" different colours rather than two shades of
        the same one. Pass ``zmin``/``zmax`` yourself to pin a scale across
        several periods; that also disables the symmetric default, so pass both.

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
            :meth:`animate` to flip through them, or :meth:`wide` to see every
            period side by side as numbers.
        layer : int, default 0
            Zero-based layer. Records in other layers are not drawn, so pass the
            layer the boundary actually sits in -- ``.get()`` shows it in the
            ``layer`` column, and ``.summary()`` counts how many layers this
            term reaches at all.
        multiplier : float, default 1.0
            Scale every value before drawing. Unit conversion (ft³/d to gpm),
            or ``-1.0`` to draw outflow as positive -- which changes the PICTURE
            only: the table keeps MF6's convention, and the colorbar will then
            read backwards from ``.get()``, so say so in the caption if you use
            it. Applied before ``agg``.
        fill_value : float, default 0.0
            Value given to cells this term has no record for. Zero is the honest
            default for a flow -- a cell with no boundary exchanges nothing --
            but pass ``float("nan")`` to leave untouched cells uncoloured, which
            is usually clearer for a sparse package such as WEL.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, default 'sum'
            How several records landing in ONE cell are reduced to one number.
            ``sum`` is right for a flow: two drains in one cell remove the sum of
            what each removes. Change it only when you want a per-feature
            statistic rather than the cell's net exchange.
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
        wide : one column per stress period, for reading the time axis.
        plot : this term as a series by stress period rather than in plan view.

        Examples
        --------
        >>> model.packages.drn.results.q.map()
        >>> model.packages.ghb.results.q.map(per=5, contours=True)
        >>> model.packages.riv.results.q.map(zmin=-25, zmax=25)
        >>> model.packages.wel.results.q.map(fill_value=float("nan"))
        >>> model.packages.uzf.results.gwrch.map(per=3, colorscale="Blues")
        >>> model.budget.sto_ss.map(per=2, backend="mpl").savefig("storage.png")
        """

        refuse_noun_parameters(self.label, self.value_name, trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.get(per=per, layer=layer)
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        if self.result_name == "q":
            # A signed exchange: centre the scale on zero so the sign is the
            # thing the colour reports. MF6's own signs are kept -- only the
            # LIMITS are made symmetric.
            absmax = _symmetric_color_limit(values)
            if absmax > 0:
                zmin = -absmax if zmin is None else zmin
                zmax = absmax if zmax is None else zmax
            # `zmid` is a Plotly trace property, not a verb parameter, so it
            # rides the open tail the way `colorbar`/`reversescale` do.
            trace_kwargs.setdefault("zmid", 0.0)
        result_spec = get_package_result_spec(self.package_name, self.result_name)
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
            hover_spec=result_hover(
                self.value_name,
                title=f"{self.package_name.upper()} {self.value_name}",
                units=(
                    {self.value_name: budget_value_units(self.model)}
                    if self.result_name == "q"
                    else None
                ),
            ),
            zmin=zmin,
            zmax=zmax,
            colorscale=(
                colorscale
                or (result_spec.colorscale if result_spec is not None else None)
                or ("RdBu" if self.value_name == "q" else None)
                or get_default_package_colorscale(self.package_name)
                or "earth"
            ),
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


class StageResultsExplorer(SpatialView):
    """Normalized explorer for cell-mapped stage results such as LAK and SFR."""

    #: value label + series column for the unified grammar
    value_name = "stage"

    def __init__(self, model: SimulationBase, package_name: str, builder):
        """Bind a stage result explorer; ``builder(model)`` yields its normalized table."""

        self.model = model
        self.package_name = str(package_name).lower()
        self._builder = builder

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (stage repeats per connected cell)."""

        return "mean"  # stage repeats per connected cell; summing is meaningless

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized stage table for the selected rows."""

        frame = self._builder(self.model)
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available stage results."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"{self.package_name}.results.stage",
            value_columns=["stage"],
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
        agg: str = "first",
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
        """Simulated stage as a choropleth: one colour per connected cell.

        A RESULTS noun over a surface-water feature's stage. Stage belongs to the
        FEATURE -- a lake or a reach -- not to a cell, so it is painted onto every
        cell that feature connects to; cells with no connection take
        ``fill_value``. That mapping is why this signature carries
        ``multiplier``/``fill_value``/``agg`` on top of the shared drawing
        parameters: the free verb never does it, because it is handed the values
        directly.

        The default ``agg="first"`` is the honest one HERE, where the budget
        nouns default to ``"sum"``: stage is an elevation, so two connections in
        one cell share one water level and adding them would invent a number.

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
            Zero-based layer. Connections in other layers are not drawn -- for a
            lake that occupies several layers this is how you pick the one you
            mean, and ``.get()`` shows which layers exist in its ``layer``
            column.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, or a datum
            shift. Applied before ``agg``.
        fill_value : float, default 0.0
            Value given to cells this package does not connect to. Pass
            ``float("nan")`` to leave them uncoloured, which is almost always
            clearer for stage: zero is a real elevation, and a plain 0.0 backdrop
            will dominate the colour scale.
        agg : {'first', 'mean', 'min', 'max', 'last', 'sum'}, default 'first'
            How several connections landing in ONE cell are reduced to one
            number. ``first`` suits a shared water level; use ``mean`` when two
            different features overlap a cell and you want the average of their
            stages rather than whichever was read first.
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
        get : the per-connection stages behind the picture, as a DataFrame.
        summary : the same records reduced to one row.
        plot : stage as a series by stress period rather than in plan view.

        Examples
        --------
        >>> model.packages.lak.results.stage.map()
        >>> model.packages.lak.results.stage.map(per=5, fill_value=float("nan"))
        >>> model.packages.sfr.results.stage.map(layer=1, colorscale="Blues")
        >>> model.packages.lak.results.stage.map(contours=True, contour_levels=6)
        >>> model.packages.lak.results.stage.map(zmin=95, zmax=105)
        >>> model.packages.sfr.results.stage.map(backend="mpl").savefig("stage.png")
        """

        refuse_noun_parameters(f"{self.package_name}.results.stage", "stage", trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.get(per=per, layer=layer)
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column="stage",
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
            hover_spec=result_hover(
                "stage",
                title=f"{self.package_name.upper()} stage",
                units={"stage": "ft"},
            ),
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


class CellPackageResultsNamespace(FieldMappable):
    """Namespace for cell-based package results represented by one budget term.

    Simple BC packages expose a single field ``q``; ``results.map()`` maps it and
    ``results.map(field="q")`` is explicit.
    """

    _default_field = "q"

    def _field_names(self):
        """Registry-declared result fields, always including ``"q"`` (the default term)."""

        names = self.fields["field"].tolist()
        return names if "q" in names else ["q", *names]  # .q always resolves

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a cell-based results namespace to ``model`` for one package."""

        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported package result fields."""

        spec = get_package_explorer_spec(self.package_name)
        if spec is None:
            return pd.DataFrame(
                columns=[
                    "field",
                    "budget_text",
                    "value_name",
                    "label",
                    "colorscale",
                    "diverging",
                ]
            )
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported package result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(getattr(self, field_name).summary())
        if not frames:
            return pd.DataFrame(
                columns=[
                    "label",
                    "records",
                    "periods",
                    "layers",
                    "cells",
                    "value_columns",
                ]
            )
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(
            lambda values: values[0] if values else None
        )
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def __getattr__(self, result_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed named result explorer."""

        result_spec = get_package_result_spec(self.package_name, result_name)
        if result_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no result {result_name!r}"
            )
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
            result_name=result_name,
        )

    @property
    def q(self) -> CellBudgetResultsExplorer:
        """Return the primary package-exchange result explorer."""

        result_spec = get_package_result_spec(self.package_name, "q")
        if result_spec is None:
            return CellBudgetResultsExplorer(
                self.model, self.package_name, self.package_name.upper(), "q",
                result_name="q",
            )
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
            result_name="q",
        )


class UzfResultsNamespace(FieldMappable):
    """Namespace for UZF result explorers.

    Fields: ``gwrch`` (groundwater recharge, the default) and ``sat``
    (unsaturated-zone saturation). ``results.map(field="sat")`` or
    ``results.sat.map()``.
    """

    _default_field = "gwrch"

    def _field_names(self):
        """Registry-declared UZF result fields, defaulting to ``["gwrch", "sat"]``."""

        names = self.fields["field"].tolist()
        return names or ["gwrch", "sat"]

    def __init__(self, model: SimulationBase):
        """Bind the UZF results namespace to ``model``."""

        self.model = model

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported UZF result fields."""

        spec = get_package_explorer_spec("uzf")
        if spec is None:
            return pd.DataFrame(
                columns=[
                    "field",
                    "budget_text",
                    "value_name",
                    "label",
                    "colorscale",
                    "diverging",
                ]
            )
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported UZF result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(self._field(field_name).summary())
        if not frames:
            return pd.DataFrame(
                columns=[
                    "label",
                    "records",
                    "periods",
                    "layers",
                    "cells",
                    "value_columns",
                ]
            )
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(
            lambda values: values[0] if values else None
        )
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def _field(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return one registry-backed UZF result explorer."""

        result_spec = get_package_result_spec("uzf", field_name)
        if result_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no UZF result field {field_name!r}"
            )
        return CellBudgetResultsExplorer(
            self.model, "uzf", result_spec.budget_text, result_spec.value_name,
            result_name=field_name,
        )

    def __getattr__(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed UZF result explorer."""

        return self._field(field_name)

    @property
    def gwrch(self) -> CellBudgetResultsExplorer:
        """Return groundwater recharge from the UZF package."""

        return self._field("gwrch")

    @property
    def sat(self) -> CellBudgetResultsExplorer:
        """Return normalized unsaturated-zone saturation results."""

        return self._field("sat")


class PackageBudgetTermExplorer:
    """Filtered helper for one package budget term or a small term family."""

    def __init__(
        self,
        namespace,
        *,
        term: str | Iterable[str],
        label: str,
    ):
        """Pin a parent budget namespace to one term (or family) under a display ``label``."""

        self._namespace = namespace
        self.term = term
        self.label = label

    @property
    def types(self) -> list[str]:
        """Return the normalized MF6 LAK term names covered by this helper."""

        return _normalize_term_filter(self.term) or []

    def get(
        self,
        *,
        per: int | Iterable[int] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Return the filtered package budget-term dataframe."""

        return self._namespace.get(term=self.term, per=per, **filters)

    def summary(
        self,
        *,
        per: int | Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Summarize the filtered package budget terms."""

        return self._namespace.summary(term=self.term, per=per, by=by, **filters)

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
        **filters,
    ) -> pd.DataFrame:
        """Pivot the filtered package budget terms to a wide dataframe."""

        return self._namespace.wide(
            term=self.term, per=per, index=index, values=values, **filters
        )


def budget_term_attribute(term: str) -> str:
    """Return the Python attribute name that reaches one MF6 budget record.

    ``"SOURCE-SINK MIX"`` -> ``source_sink_mix``, ``"FLOW-JA-FACE"`` ->
    ``flow_ja_face``, ``"STO-SS"`` -> ``sto_ss``, ``"UZF-GWRCH"`` ->
    ``uzf_gwrch``. MF6 separates words with hyphens and, in one case, a space;
    both collapse to the underscore.

    This is the FIRST normalizer in the codebase to run this direction. The
    package-level ``budget.<term>`` namespaces hand-write both spellings as
    literals, and their only normalizer (``_normalize_term_filter``) upcases
    without converting ``_`` back to ``-`` -- so ``budget.get(term="ext_inflow")``
    silently returns an empty frame there. Deriving the mapping in one place is
    what lets :meth:`ModelBudgetNamespace.__getitem__` accept either spelling.
    """

    return re.sub(r"[^0-9a-z]+", "_", str(term).strip().lower()).strip("_")


class ModelBudgetNamespace:
    """Every term in this model's OWN budget file, as a noun.

    ``model.budget.<term>`` returns a :class:`CellBudgetResultsExplorer`, so each
    term answers the full spatial verb set (``get``/``summary``/``plot``/``map``/
    ``xs``/``mosaic``/``animate``) rather than the reduced set the package-level
    ``<pkg>.budget.<term>`` namespaces offer.

    Terms are **discovered at runtime**, not declared. That is a deliberate
    departure from the hand-written package namespaces, because the model
    budget's term set genuinely varies with kind and configuration: GWT's storage
    term is ``STORAGE-AQUEOUS`` and GWE's is ``STORAGE-CELLBLK``, ``DECAY``
    appears only when MST declares decay or sorption, and a GWF model's terms are
    whichever boundary packages it happens to carry. Declaring them would mean
    hand-maintaining a list that is wrong for most models.

    Available on EVERY model kind. The plumbing underneath
    (``_get_budget_reader``, and the record converter) is kind-neutral, so gating
    this to transport would have been artificial -- and ``STO-SS``, ``STO-SY``,
    ``DATA-SAT`` and ``DATA-SPDIS`` have no package accessor at all, making this
    their only route that is not the legacy ``model.bud(...)`` wrapper.

    Not to be confused with its three neighbours:

    * ``model.bud(pkg)`` -- the legacy compatibility wrapper, raw frames.
    * ``model.budget_cumulative`` / ``model.budget_incremental`` -- whole-model
      totals from the LISTING file, not per-cell.
    * ``model.packages.<pkg>.budget.<term>`` -- a genuinely DIFFERENT file, the
      package-output budget, whose node layout is feature-first. Conflating the
      two is what produced the off-by-one in ledger 92.
    """

    def __init__(self, model: SimulationBase):
        """Bind the model-budget term namespace to ``model``."""

        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the MF6 record names present in this model's budget file."""

        reader = self.model._get_budget_reader()
        return [
            str(name).strip()
            for name in reader.get_unique_record_names(decode=True)
        ]

    def _terms(self) -> dict[str, str]:
        """Map each attribute name to the MF6 record name it reaches."""

        return {budget_term_attribute(name): name for name in self.types}

    def __dir__(self) -> list[str]:
        """List the real terms alongside the normal members, for autocomplete."""

        try:
            terms = self._terms()
        except Exception:  # noqa: BLE001 - dir() must never raise
            # 7.3 left this broad on purpose. `_terms` opens the budget file,
            # and that chain was measured to raise types no builtin tuple
            # covers -- flopy's MFDataException, and NotImplementedError out of
            # the base Grid.shape that CellBudgetFile touches unconditionally
            # when a model has no discretization. An autocomplete that raises
            # is worse than one that comes back short.
            logger.debug("no budget terms for __dir__", exc_info=True)
            terms = {}
        return sorted({*super().__dir__(), *terms})

    def __getattr__(self, name: str):
        """Resolve one budget term to its explorer."""

        # `model` would recurse (it is looked up by _terms below) and dunder /
        # private probes must fail fast rather than open the budget file.
        if name.startswith("_") or name == "model":
            raise AttributeError(name)

        terms = self._terms()
        record = terms.get(name)
        if record is None:
            raise AttributeError(
                f"model {self.model.name!r} has no budget term {name!r}. "
                f"Available terms: {sorted(terms)}"
            )
        return CellBudgetResultsExplorer(
            self.model,
            # The MF6 record name, lowercased -- NOT the attribute spelling. For
            # a term that is also a package ("DRN") this keeps the registry
            # lookup working, so the term inherits that package's colorscale and
            # result spec; for the rest it makes the hover read "SOURCE-SINK MIX"
            # the way MF6 spells it, rather than "SOURCE_SINK_MIX".
            package_name=record.lower(),
            budget_text=record,
            value_name="q",
            result_name="q",
            label=f"budget.{name}",
        )

    def __getitem__(self, term: str):
        """Reach a term by EITHER its attribute name or its MF6 record name.

        ``model.budget["SOURCE-SINK MIX"]`` and ``model.budget["source_sink_mix"]``
        are the same object, so a term name copied straight out of ``types`` (or
        out of an MF6 listing file) always works.
        """

        return getattr(self, budget_term_attribute(term))

    def __repr__(self) -> str:
        """Show which terms this model's budget file actually carries."""

        try:
            terms = ", ".join(sorted(self._terms())) or "no terms"
        except Exception:  # noqa: BLE001 - repr must never raise
            # Broad on purpose, same reasoning as __dir__ above. A __repr__ that
            # raises breaks the debugger and the traceback you were reading when
            # you needed it -- so this one reports the failure in its own text.
            logger.debug("no budget terms for __repr__", exc_info=True)
            terms = "budget file unavailable"
        return f"<{type(self).__name__} {getattr(self.model, 'name', '?')!r}: {terms}>"


__all__ = [
    "CellBudgetResultsExplorer",
    "StageResultsExplorer",
    "CellPackageResultsNamespace",
    "UzfResultsNamespace",
    "PackageBudgetTermExplorer",
    "ModelBudgetNamespace",
    "budget_term_attribute",
]
