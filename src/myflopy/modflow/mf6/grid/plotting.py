"""Plotting helpers for Voronoi grids and grid-derived cross sections."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import shapely as shp
from flopy.discretization.vertexgrid import VertexGrid
from flopy.plot.crosssection import PlotCrossSection

from myflopy import viz as f
from myflopy.modflow.mf6.cross_section_plotting import (
    filled_section,
    mirror_horizontally,
)
from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    as_mpl_figure,
    normalize_backend,
)
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.hover import conc_hover, head_hover, temp_hover, trans_hover
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg
from myflopy.viz import Picture

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from myflopy.modflow.mf6.simulation.base import SimulationBase


#: The sectioned hover each dependent variable gets by default. Without one,
#: `Choro` falls back to a flat `Cell No. / Area / x / y / ...` dump -- readable
#: only if you already know what you are looking at.
_HOVER_FOR_TYPE = {"hds": head_hover, "conc": conc_hover, "temp": temp_hover}


def _default_hover_spec(map_type: str, model=None):
    """The default :class:`HoverSpec` for a map of ``map_type``, or None.

    The three dependent variables and the two transmissivities get one -- heads
    and the transmissivities labelled in the model's own units, read from
    ``model``.
    ``type="custom"`` is the package/group path, which supplies its own
    ``custom_hover``; ``"rch"``/``"ks"`` read fields whose hover is assembled
    inline. Returning None for those leaves their existing behaviour exactly as
    it was.
    """

    key = str(map_type).lower()
    length, time = unit_labels(model) if model is not None else (None, None)
    if key in ("trans", "trans_max"):
        return trans_hover(saturated=key == "trans", length=length, time=time)
    if key == "hds":
        # A head is a length: written in the model's own unit, never assumed
        # (ledger 403).
        return head_hover(unit=length)
    builder = _HOVER_FOR_TYPE.get(key)
    return builder() if builder is not None else None


def _choropleth_factory(
    vor,
    model: SimulationBase = None,
    kstpkper: tuple = None,
    per: int = None,
    layer: int = 0,
    type: str = 'hds',
    custom_hover: dict = None,
    custom_zs: list = None,
    zmin: float | int = None,
    zmax: float | int = None,
    zoom: int = 13,
    show_layer_elevs: bool | None = None,
    show_mounding: bool = False,
    hover_heads: bool = True,
    hover_ks: bool = False,
    features=None,
    colorscale: str | list | tuple = None,
    logscale: bool = False,
    hover_spec=None,
    **choro_kwargs,
) -> Choro:
    """
    Create a Choro wrapper for Voronoi plotting.

    ``show_layer_elevs=None`` (the default) means **decide from the grid**: the
    layer-elevation hover needs ``vor.gdf_topbtm``, and a grid without it makes
    Choro's own ``True`` default raise ``AttributeError`` as soon as the hover is
    built. Pass ``True``/``False`` to force it.

    That resolution used to live at the CALL SITES -- 25 of them repeated
    ``kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(model))``
    around a hardcoded ``False`` here, and the two that forgot silently lost five
    hover rows. It belongs in one place, next to the attribute it depends on. Do
    not "simplify" this back into a bare passthrough.

    ``colorscale``/``logscale``/``hover_spec`` are named because callers reach
    for them constantly; anything else in ``**choro_kwargs`` rides through to the
    ``go.Choroplethmap`` trace (``zmid``, ``colorbar``, ``reversescale``, ...).
    Those are validated LATE, by Plotly at ``plot()`` time, not here -- ``title=``
    in particular is a matplotlib-only argument and raises there.
    """
    head_options = [
        name for name, asked in (("bgs", choro_kwargs.get("bgs")), ("show_mounding", show_mounding))
        if asked
    ]
    if head_options and model is None:
        # Both read a run model's HEADS. With none bound -- a bare grid, an unrun
        # model, a barrier map's cell-number backdrop -- the Choro was built and
        # `.fig` raised `AttributeError: 'NoneType' object has no attribute
        # 'all_values'`, far from the call (ledger 436).
        raise TypeError(
            f"{' and '.join(f'{name}=True' for name in head_options)} "
            f"{'map' if len(head_options) > 1 else 'maps'} a run "
            "model's HEADS (depth below ground, or head above the layer bottom), and "
            "this map has none to read: a bare grid (vor.plot.map), a model not yet "
            "run, or a barrier map's backdrop. Map the run model instead: "
            "model.plot.map(bgs=True) or model.hds.map(show_mounding=True)."
        )
    inputs = model if model is not None else choro_kwargs.get("input_model")
    if hover_spec is None and custom_hover is None and inputs is not None:
        hover_spec = _default_hover_spec(type, inputs)
    if zmin is not None and zmax is not None and zmin >= zmax:
        # Kept from `ModelVisualization.plotly_head_map_animation`, which 8.6b
        # deleted -- and moved DOWN here, so it now guards every map rather than
        # only the animated head map. An inverted range renders an all-one-colour
        # picture with no error, which reads as a broken model.
        raise ValueError(f"zmin must be less than zmax; got zmin={zmin}, zmax={zmax}.")
    if show_layer_elevs is True:
        # An EXPLICIT request also turns on the sectioned hover's surfaces, which
        # decides from the spec's own flag and so ignored this one entirely.
        # Only when explicit: `None` resolves to True merely because the grid HAS
        # a layer frame, which means "we could", not "you asked" -- and defaulting
        # every head map to the full stacked table is exactly what
        # `test_head_map_active_strip_is_default_and_compact` forbids.
        # `setdefault` is not enough: the free verb forwards `hover_surfaces=None`
        # explicitly, so the key is already present. Test the VALUE, which also
        # leaves an explicit `hover_surfaces=False` winning.
        if choro_kwargs.get("hover_surfaces") is None:
            choro_kwargs["hover_surfaces"] = True
    if show_layer_elevs is None:
        show_layer_elevs = getattr(vor, "gdf_topbtm", None) is not None
    return Choro(
        vor=vor,
        model=model,
        kstpkper=kstpkper,
        per=per,
        layer=layer,
        type=type,
        custom_hover=custom_hover,
        custom_zs=custom_zs,
        zmin=zmin,
        zmax=zmax,
        zoom=zoom,
        show_layer_elevs=show_layer_elevs,
        show_mounding=show_mounding,
        hover_heads=hover_heads,
        hover_ks=hover_ks,
        features=features,
        colorscale=colorscale,
        logscale=logscale,
        hover_spec=hover_spec,
        **choro_kwargs,
    )


def _grid_section_factory(vor: VoronoiGridPlus, line: shp.LineString | Path):
    """
    Build a cross-section helper for the Voronoi grid.
    """
    return GridSection(vor=vor, line=line)


def _as_linestring(geometry) -> shp.LineString:
    """Coerce a ``LineString`` or ``MultiLineString`` to a single merged ``LineString`` (raises otherwise)."""

    if isinstance(geometry, shp.LineString):
        return geometry

    if isinstance(geometry, shp.MultiLineString):
        merged = shp.line_merge(geometry)
        if isinstance(merged, shp.LineString):
            return merged

        coords = []
        for line in geometry.geoms:
            line_coords = list(line.coords)
            if coords and coords[-1] == line_coords[0]:
                coords.extend(line_coords[1:])
            else:
                coords.extend(line_coords)
        return shp.LineString(coords)

    raise ValueError(f'line arg must resolve to a LineString, not {type(geometry)}')


class GridSection(Picture):
    """
    Represents a section of a grid and provides tools for creating and plotting
    cross-sections.
    """

    #: Whether `.fig` has assembled its traces. A CLASS attribute, matching the
    #: other pictures, so instances built via `object.__new__` still answer.
    _assembled = False
    _fig = None

    def __init__(self, vor, line: shp.LineString | shp.MultiLineString | Path):
        """Build a cross-section of grid ``vor`` along ``line`` (a geometry or vector file)."""

        self.vor = vor
        props = vor.get_disv_gridprops()
        self.grid = VertexGrid(
            vertices=props['vertices'],
            top=vor.gdf_topbtm[0].values,
            botm=vor.gdf_topbtm.loc[:, 1:].values.T,
            cell2d=props['cell2d'],
            lenuni='feet',
            ncpl=props['ncpl'],
            crs=vor.crs,
            nlay=vor.nlay,
        )

        if isinstance(line, Path):
            geometry = read_shp_gpkg(line).union_all()
            self.coords = _as_linestring(geometry).coords
        elif isinstance(line, (shp.LineString, shp.MultiLineString)):
            self.coords = _as_linestring(line).coords
        else:
            raise ValueError(f'line arg must be a Path or LineString, not {type(line)}')

        self.xy = np.array([xy for xy in self.coords])

    @property
    def polys(self):
        """Return FloPy cross-section polygons for the selected line."""
        polys = PlotCrossSection(
            modelgrid=self.grid,
            line={'line': self.xy},
        ).polygons
        return polys

    @property
    def poly_coords(self):
        """Return raw vertex arrays for each cross-section polygon."""
        poly_coords = []
        for poly in self.polys.values():
            verts = poly[0].get_xy()
            poly_coords.append(verts)
        return poly_coords

    def to_frame(self) -> pd.DataFrame:
        """
        Return section polygon outlines as a long-form DataFrame.
        """
        rows = []
        for i, verts in enumerate(self.poly_coords):
            series_name = f"polygon_{i}"
            for distance, elevation in verts:
                rows.append(
                    {
                        "distance": float(distance),
                        "elevation": float(elevation),
                        "series": series_name,
                        "polygon_id": i,
                    }
                )
        return pd.DataFrame(rows)

    def plot_mpl(self, **kwargs):
        """Plot the cell outlines this section cuts, with Matplotlib.

        One line per cell polygon -- distance along the section against
        elevation, both in the grid's length units -- drawn by the ``figs``
        cross-section helper (:func:`myflopy.viz.plot_cross_section`) from
        :meth:`to_frame`. It always makes a NEW figure: the helper takes no Axes.

        Parameters
        ----------
        **kwargs
            Forwarded to :func:`myflopy.viz.plot_cross_section`. ``show_legend``
            defaults to ``False`` here (a legend would name every polygon);
            ``title``, ``xlabel`` (default ``'Distance'``) and ``ylabel`` (default
            ``'Elevation'``), ``style``, ``x_units_per_inch``/``y_units_per_inch``
            (an exact print scale, applied only when both are set),
            ``line_kwargs`` (to ``Axes.plot``) and ``scale_kwargs`` act; anything
            else goes on to the figure factory (``figsize``, default ``(11, 8.5)``
            inches; ``font_size``; ...), which raises ``TypeError`` on a keyword it
            does not know. ``data``/``x``/``y``/``series_col`` are this method's and raise
            ``TypeError`` as duplicates; ``layer``, ``simplify_tolerance``,
            ``sheet_name`` and ``read_kwargs`` read a FILE source and are ignored,
            since this hands the helper a table.

        Returns
        -------
        tuple of (matplotlib.figure.Figure, matplotlib.axes.Axes)
        """
        from myflopy.viz import plot_cross_section

        data = self.to_frame()
        kwargs.setdefault("show_legend", False)
        return plot_cross_section(
            data=data,
            x="distance",
            y="elevation",
            series_col="series",
            **kwargs,
        )

    @property
    def fig(self):
        """The assembled cross-section figure (built once, then reused).

        Cached because `Picture` requires repeated access to return the SAME
        figure, so `section.fig.update_layout(...)` then `section.show()` acts
        on one figure.
        """
        if self._assembled:
            return self._fig
        fig = f.Fig()
        for verts in self.poly_coords:
            xs, ys = verts[:, 0], verts[:, 1]
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    fill='toself',
                    mode='lines',
                    line=dict(color='black'),
                    name='Polygon',
                )
            )
        self._fig = fig
        self._assembled = True
        return fig



class GridMesh(Picture):
    """The bare mesh -- cell edges, no values, no basemap (plan 8.4).

    The picture ``map()`` cannot give you. A :class:`Choro` colours cells against
    a web basemap and so hard-requires a CRS; this draws the geometry in the
    grid's own coordinates and needs none, which is what makes it the right view
    for a grid you are still refining, before there is a model or a projection to
    speak of.

    Absorbs the old ``vor.plot2d()`` (Plotly) and the inherited FloPy
    ``VoronoiGrid.plot()`` (Matplotlib, still reachable as :meth:`plot_mpl`).
    """

    #: Whether `.fig` has assembled its traces. A CLASS attribute, matching
    #: `Choro`/`XSection`, so instances built via `object.__new__` still answer.
    _assembled = False
    _fig = None

    def __init__(self, vor):
        """Bind a mesh view to grid ``vor``."""

        self.vor = vor

    @property
    def fig(self):
        """Cell edges as one Plotly trace per cell, in the grid's own coordinates.

        Assembled once and cached: `Picture` requires repeated access to return
        the SAME figure, so `mesh.fig.update_layout(...)` then `mesh.show()`
        acts on one figure.
        """

        if self._assembled:
            return self._fig
        fig = f.Fig(layout=getattr(self.vor, "scatt_layout", None))
        for cell in range(len(self.vor.x_coords_by_node)):
            fig.add_scattergl(
                x=self.vor.x_coords_by_node[cell],
                y=self.vor.y_coords_by_node[cell],
                opacity=1,
                mode='lines',
                line_color='black',
                line_width=1,
                showlegend=False,
            )
        self._fig = fig
        self._assembled = True
        return fig

    def plot_mpl(self, ax=None, plot_title: bool = True, **kwargs):
        """Render with Matplotlib, via FloPy's own patch-collection renderer.

        This IS ``VoronoiGrid.plot`` -- called unbound, because ``vor.plot`` is
        now the namespace this object came from. Nothing about the drawing
        changed; only the spelling did. ``vor.plot.grid(backend="mpl")`` is the
        same drawing returned as its ``Figure``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The Axes to draw the cell edges into, in the grid's own coordinates.
            ``None`` (the default) draws into a NEW figure of one equal-aspect
            Axes, as every other static picture does. It used to call FloPy's
            ``plt.subplot(1, 1, 1, aspect="equal")``, which drew onto pyplot's
            CURRENT figure -- over a figure of the caller's own (ledger
            436).
        plot_title : bool, default True
            Title the Axes ``"ncells: <ncpl>; nverts: <nverts>"``. ``False`` leaves
            the title alone.
        **kwargs
            Forwarded through FloPy's ``VoronoiGrid.get_patch_collection`` to
            ``PlotMapView.plot_grid``: ``colors`` (default ``'grey'``; ``color``,
            ``ec`` and ``edgecolor`` are read as it), and everything else to
            Matplotlib's ``LineCollection`` -- ``linewidths``, ``alpha``, ... An
            unknown keyword raises ``AttributeError`` there. ``ax`` cannot be
            passed here; it is the parameter above.

        Returns
        -------
        matplotlib.axes.Axes
            The Axes drawn on.
        """

        from flopy.utils.voronoi import VoronoiGrid

        if ax is None:
            _, ax = f.mpl_axes(subplot_kw={"aspect": "equal"})
        return VoronoiGrid.plot(self.vor, ax=ax, plot_title=plot_title, **kwargs)


class GridPlots:
    """The plotting verbs for a Voronoi grid -- ``vor.plot.map()`` (plan 8.4).

    Three verbs, because a bare grid can answer three questions: what are the
    values over it (:meth:`map`), what does it look like in section
    (:meth:`section`), and what does the mesh itself look like (:meth:`grid`).
    Everything they return is a :class:`~myflopy.viz.Picture`.

    Calling the namespace is shorthand for the mesh: ``vor.plot()`` is
    ``vor.plot.grid()``. That is deliberate -- ``vor.plot`` used to BE FloPy's
    ``VoronoiGrid.plot()``, and "draw the grid" is what people already reach for
    that name to do.

    Replaces eleven aliases. Five were genuinely distinct pictures and survive as
    these verbs or as options on them; the rest were dead, duplicated, or not
    pictures at all -- see the ledger.

    Lives here rather than in :mod:`myflopy.plot` for a layering reason:
    ``voronoi.py`` sits BELOW that module in the import graph, so binding this
    from there would point upward and force a deferred import.
    """

    def __init__(self, vor):
        """Bind the plotting verbs to grid ``vor``."""

        self.vor = vor

    def __repr__(self):
        """Name the verbs, since tab-completion is how this gets discovered."""

        return f"GridPlots({getattr(self.vor, 'ncpl', '?')} cells: map, section, grid)"

    def __call__(self) -> GridMesh:
        """``vor.plot()`` -> the mesh. See :meth:`grid`."""

        return self.grid()

    def map(
        self,
        values=None,
        *,
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        features=None,
        hillshade_path=None,
        zoom: int = 13,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        show_layer_elevs: bool | None = None,
        custom_hover: dict | None = None,
        hover=None,
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """A plan-view map of this grid, on a basemap.

        ``values`` is any per-cell array; with none, cells are keyed by node id
        -- which is what the old ``map_nodes()`` drew.

        Requires a CRS, because the basemap does. Use :meth:`grid` for a grid
        that does not have one yet.

        A DELIBERATE SUBSET of :func:`myflopy.plot.map`. Everything omitted --
        ``per``, ``kstpkper``, ``layer``, ``type``, ``bgs``, ``show_mounding``,
        the sectioned-hover knobs -- reads a model's results, and a bare grid has
        none, so advertising them here would promise something this scope cannot
        do. They still ride the ``**trace_kwargs`` tail if you pass one, because
        narrowing the signature should not narrow what already worked -- bar
        ``bgs`` and ``show_mounding``, which read heads and raise ``TypeError``
        here (they built, and raised ``AttributeError`` at ``.fig``; ledger
        465).

        .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.map`.

        The single map verb. Contours, observation markers, a hillshade and particle
        pathlines are all **options** here rather than verbs of their own, because a
        map is a plan view whatever is drawn on it.

        Every parameter below is spelled out in the signature rather than swept into
        ``**kwargs``, so an editor can complete and type-check them. The defaults are
        not restated here by hand -- each mirrors the default of whichever link in the
        ``map -> _choropleth_factory -> Choro`` chain owns that argument, and
        ``test_plot_vocabulary`` fails if the two ever drift apart.

        .. end generated

        Parameters
        ----------
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight: indices, a boolean mask, a registered region
            name, a vector file, or a geometry to intersect. Replaces
            ``show_selected_cells()`` and ``show_overlapping_geometry()``.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            ``"outline"`` (default) draws the dissolved boundary of the
            selection and leaves the cells at full opacity; ``"dim"`` fades
            everything else instead, which reads well on a bare node-id grid and
            badly on a field; ``"both"`` does each.
        select_color : str, optional
            Highlight colour; defaults to the palette's.
        backend : {'plotly', 'mpl'}, default 'plotly'
            ``"plotly"`` (default) returns the interactive ``Choro`` picture;
            ``"mpl"`` returns a static :class:`matplotlib.figure.Figure`. The
            same switch every other picture verb takes, so the spelling does not
            change with the scope. Note the Matplotlib branch draws in model
            coordinates with no basemap, which makes ``zoom`` inert there.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        values : array-like of float, optional
            ONE layer's field, one value per cell in cell order -- heads, drawdown, K,
            a zone id, a residual. Shape ``(ncpl,)``, ``ncpl`` being the cells the map
            draws (``vor.ncpl``); a single ``(1, ncpl)`` row is read as that row, a
            structured model's ``(nrow, ncol)`` in node order, and a ``pandas.Series``
            by POSITION, not by its index. A masked array's masked cells draw as gaps.
            Overrides whatever ``type`` would have read; ``layer=`` still names the
            layer the values belong to, which the hover and the contour clip read.
            Anything else RAISES at the call, on both backends, naming ``(ncpl,)``: a
            field for every layer -- ``(nlay, ncpl)``, as ``ies.calibrated().get("npf",
            "k")`` returns, or flattened -- names ``values=array[layer], layer=layer``
            (``ValueError``); a dict, which is how a per-period dataset arrives, names
            ``values=data[period]`` (``TypeError``); a ``DataFrame`` or a single number
            is a ``TypeError``, and any other length a ``ValueError`` (ledger
            601). Default ``None``: the field ``type`` names.
        per : int, optional
            Stress period to read (0-based). Mutually exclusive with ``kstpkper``;
            with neither, the model's first output time is used.
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)`` to read, as MODFLOW reports it. Use
            ``model.kstpkper`` to list what is available.
        per_timestep : {'last', 'first'} or int, default 'last'
            Which timestep WITHIN ``per`` to read, when a period has several. Ignored
            when ``kstpkper`` names the timestep outright.
        layer : int, default 0
            Zero-based layer. Layer 0 is the top. It must be one of the layers the
            model's views number -- the grid's, and ONE on a DISU model, whose nodes
            MF6 writes as one record whatever layering FloPy infers -- or the call
            raises ``ValueError`` naming the count (ledger 448); ``-1`` is
            accepted only beside ``show_mounding=True``. A bare grid checks nothing.
        type : {'hds', 'conc', 'temp', 'rch', 'ks', 'trans', 'trans_max', 'custom'}, default 'hds'
            Which field to read when ``values`` is not given. ``'hds'`` heads,
            ``'conc'`` GWT concentration, ``'temp'`` GWE temperature -- these three
            also select the default sectioned hover -- ``'rch'`` the simulated
            recharge rate, and ``'ks'`` the horizontal conductivity ``npf.k``.
            ``'trans'`` is TRANSMISSIVITY, K times the SATURATED thickness at this
            map's ``per``/``kstpkper`` -- ``min(head, top) - bottom`` in a
            convertible layer, the full thickness in a confined one (MF6's own
            ICELLTYPE rule) -- and ``'trans_max'`` is K times the FULL thickness;
            both hover T with the K and thickness it multiplies, in the model's
            units. ``'ks'`` and ``'trans_max'`` are INPUTS, so they draw on a model
            that has not been run; the others are results and raise there. Case does
            not matter; any other value raises, naming this set, rather than drawing
            node numbers. The same two fields are nouns too:
            ``model.npf.transmissivity`` / ``.transmissivity_max``, whose
            ``.mosaic()`` draws one panel per layer.
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
        bgs : bool, default False
            Map depth below ground: each cell's head minus the model top, negative
            below it. Blank where the layer has no head, dry or inactive (ledger 315).
            ``show_mounding`` takes precedence.
        zoom : int, default 13
            Initial map zoom. Ignored when ``fit_bounds`` is True.
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
        hover_layers : {'active', 'active+strip', 'all', 'none'}, optional
            How the hover reads a field with one value per layer (heads,
            concentration, temperature, transmissivity): ``'active'`` the mapped
            layer's value alone; ``'active+strip'`` adds one line with the other
            layers' values; ``'all'`` a table with one row per layer, the mapped one
            accented; ``'none'`` no layer readout -- the value line goes too, leaving
            the title, ``hover_fields`` and the footer (with ``hover_surfaces``, the
            elevations alone). Default ``None`` keeps the map's own default
            (``'active+strip'`` on those fields). Any other value raises
            ``ValueError`` at the call, naming the four. No effect on a field with
            one value per cell, nor on the legacy ``custom_hover``.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        hover_fields : sequence of str, optional
            Extra columns to append to the hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        show_mounding : bool, default False
            Map mounding instead of head: each cell's head above its layer's bottom
            (above the model top with ``layer=-1``), clipped at 0, and add it to the
            hover. A DRY cell reads 0 and an inactive one is blank (ledger 315).
        hover_heads : bool, default True
            Include heads in the legacy flat hover.
        hover_ks : bool, default False
            Include hydraulic conductivity in the legacy flat hover.
        custom_hover : dict, optional
            Legacy flat hover: ``{label: per-cell sequence}``. Supplying it
            suppresses the default sectioned hover.
        rch_scale : float, optional
            Multiplier applied to recharge values when ``type='rch'``.
        animation_kstpkpers : sequence of tuple, optional
            The output times ``.ani`` steps through. Defaults to every time the model
            wrote.
        **trace_kwargs
            Anything else rides through to the ``go.Choroplethmap`` trace --
            ``zmid``, ``colorbar``, ``reversescale``, ``showscale``. These are
            genuinely open-ended and Plotly owns their names, so they are validated
            LATE, at render time, not here.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default), a
            :class:`~myflopy.viz.Picture`: it renders itself in Jupyter, and answers
            ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. It also
            carries ``.plot_mpl()`` for a static rendering and ``.ani`` for the
            animation over periods.

            With ``backend='mpl'``, a bare Matplotlib ``Figure`` -- not a Picture, so
            use ``.savefig(path)`` and ``.axes[0]`` rather than the picture verbs.

        See Also
        --------
        section : the same data as a vertical slice.
        grid : the mesh with no values and no basemap (and no CRS needed).
        myflopy.plot.animate : flip a sequence of these through time.

        Examples
        --------
        >>> model.plot.map(layer=0)                          # this model's heads
        >>> model.plot.map(values=drawdown, layer=0)         # any per-cell array
        >>> model.plot.map(layer=0, contours=True, contour_levels=8)
        >>> model.plot.map(layer=0, features="wells.gpkg", hillshade_path="hs.tif")
        >>> model.plot.map(values=k[2], layer=2, logscale=True,
        ...                colorscale="Viridis")  # one layer of K, over decades
        >>> model.plot.map(layer=0, features={"Wells": "wells.gpkg",
        ...                                   "Ponds": plot.Feature("ponds.gpkg", fill=True)})
        >>> model.plot.map(layer=0, zmin=100, zmax=125).save("heads.png")
        >>> vor.plot.map(values=node_ids)                    # a bare grid
        """

        if "type" in trace_kwargs:
            # It rode the tail to a picture with no model to read it from, and drew
            # node numbers under a K or recharge colorbar (ledger 311).
            raise TypeError(
                f"type={trace_kwargs['type']!r}: a bare grid has no model to read a "
                "field from. Map it from the model (model.plot.map(type=...)), or "
                "pass the numbers as values=."
            )
        picture = _choropleth_factory(
            self.vor,
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale,
            logscale=logscale,
            zoom=zoom,
            features=features,
            custom_hover=custom_hover,
            show_layer_elevs=show_layer_elevs,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            hover=hover,
            select=select,
            select_style=select_style,
            select_color=select_color,
            **({"custom_zs": values} if values is not None else {}),
            **trace_kwargs,
        )
        return _apply_backend(picture, backend)

    def section(
        self,
        line,
        *,
        fill="layer",
        fill_cmap: str = "viridis",
        fill_label: str | None = None,
        reverse: bool = False,
        backend: str = "mpl",
    ):
        """The grid's cross-section along ``line``: its layers, or its cell outlines.

        A bare grid has no results, so this is the geometry alone -- by default
        the cells the line cuts, coloured by layer; with ``fill=None``, their
        outlines. For a section through a model's results use
        ``model.plot.section(...)``.

        .. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.section`.

        ``section`` draws one of two pictures, and ``fill`` chooses which:

        ===============================  =============================================
        ``fill="layer"`` (the default)   The CROSS-SECTION: every cell the line cuts,
                                         drawn from its bottom to its top and coloured
                                         by layer, with the simulated water level on
                                         top.
        ``fill="results"``               The same cells, coloured by the model's own
                                         result -- head on GWF, concentration on GWT,
                                         temperature on GWE -- with a colorbar.
        ``fill=<(nlay, ncpl) array>``    The same cells, coloured by any field of your
                                         own: K, a zone id, a residual.
        ``fill=None``                    The LINE PROFILE: head against distance along
                                         the line, one line per ``layer=``, plus the
                                         model top. No cells.
        ===============================  =============================================

        Both pictures draw on both backends, and each backend draws the same cells.
        ``backend="mpl"`` -- the default -- returns a static Matplotlib ``Figure``;
        ``backend="plotly"`` returns an interactive picture, which for the line
        profile can also step through the saved times (``.ani``).

        Through a bare GRID (``vor.plot.section(line)``) there are no results: the
        cross-section is the layer geometry alone, ``fill="results"`` raises, and
        ``fill=None`` draws the cell outlines.

        **Which arguments go with which picture.** An argument given to the other
        picture RAISES and says which picture it belongs to, rather than being
        quietly ignored:

        - *either picture:* ``line`` or ``cells``, ``per`` or ``kstpkper``,
          ``reverse``, ``section_name``, ``backend``;
        - *the cross-section:* ``layers``, ``head_layers``, ``layer_labels``,
          ``show_grid``, ``legend``, ``fill_cmap``, ``fill_label``;
        - *the line profile (``fill=None``):* ``layer``, ``interpolate``,
          ``interpolator``, ``use_rbf``, ``num_points``, ``clip``, ``x_or_y``,
          ``show_model_top``, ``show_model_btm``, ``surf_type``,
          ``animation_kstpkpers``.

        Mind ``layers=`` against ``layer=``: the cross-section's chooses whose CELLS
        are drawn, the profile's whose head LINES.

        .. end generated

        Parameters
        ----------
        line : LineString, MultiLineString, or Path
            The section line, as a geometry or a vector file to read it from.
        fill : {'layer'}, array-like, or None, default 'layer'
            What colours the cells: their layer (``'layer'``), or a
            ``(nlay, ncpl)`` array. ``None`` draws the cell outlines instead. A
            bare grid has no results, so ``'results'`` raises here --
            ``model.plot.section(fill="results")`` is the picture for that.
        fill_cmap : str, default 'viridis'
            A Matplotlib colormap name for an array ``fill``.
        fill_label : str, optional
            The colorbar title for an array ``fill``.
        backend : {'mpl', 'plotly'}, default 'mpl'
            ``'mpl'`` returns a static :class:`matplotlib.figure.Figure`;
            ``'plotly'`` an interactive picture -- a ``FilledSection``, or with
            ``fill=None`` the :class:`GridSection` outlines. Both backends draw
            the same cells.
        reverse : bool, default False
            Traverse the line the other way -- a section drawn north-south comes
            back south-north, mirrored about its own centre.

            Mirrors the finished DRAWING, not the line, and so does not work the
            way you would do it by hand: FloPy's ``PlotCrossSection`` normalizes
            direction internally, so handing it reversed coordinates gives it
            byte-identical geometry and changes nothing. The view is inverted
            instead and the distance axis relabelled, so it still reads from the
            left edge.

            Both backends mirror the filled section. With ``fill=None`` only
            ``backend="mpl"`` can: the Plotly outline picture is not an axes, and
            raises rather than ignoring the argument.

        Returns
        -------
        matplotlib.figure.Figure, FilledSection or GridSection
            A Matplotlib ``Figure`` under the default ``backend="mpl"``. Under
            ``"plotly"``, a :class:`~myflopy.viz.Picture`: a ``FilledSection``,
            or with ``fill=None`` the ``GridSection`` outlines.

            .. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
        Raises
        ------
        ValueError
            An argument that belongs to the other picture, naming the picture it
            belongs to; an interpolation argument without interpolation, naming
            ``interpolate=True``; ``per`` and ``kstpkper`` together; a period or
            time the output did not save; a model's argument given to a bare grid;
            ``fill='results'`` on a bare grid; a flat fill array; a line that crosses
            no cell.
        TypeError
            A keyword ``section`` does not take at all -- a typo -- on either picture.

        See Also
        --------
        map : the same model in plan view.
        myflopy.plot.grid : the mesh itself.

            .. end generated

        Examples
        --------
        >>> vor.plot.section(line)                        # the layers, by colour
        >>> vor.plot.section(line).savefig("section.png")
        >>> vor.plot.section(line, backend="plotly")      # the same, interactive
        >>> vor.plot.section(line, fill=None)             # the cell outlines
        """

        if fill is not None:
            # Straight to the layer-1 renderer. Routing through
            # `myflopy.plot.section` would be an UPWARD import from this module.
            return filled_section(
                _grid_section_factory(self.vor, line=line).grid, line,
                fill=fill, cmap=fill_cmap, label=fill_label, reverse=reverse,
                backend=normalize_backend(backend),
            )
        if fill_cmap != "viridis" or fill_label is not None:
            raise ValueError(
                "fill_cmap= and fill_label= colour the cells of the "
                "cross-section; fill=None draws only their outlines. Drop "
                "fill=None, or drop them."
            )
        if reverse and normalize_backend(backend) != "mpl":
            raise ValueError(
                "reverse= mirrors the drawn section by inverting its axes, "
                "which is a Matplotlib operation; the Plotly geometry section "
                "is not an axes. Pass backend='mpl'."
            )
        picture = _grid_section_factory(self.vor, line=line)
        if normalize_backend(backend) == "plotly":
            return picture
        figure = as_mpl_figure(picture.plot_mpl())
        if reverse:
            # The geometry section draws cells too, so it mirrors exactly as the
            # filled one does -- `mirror_horizontally` inverts the view rather
            # than the drawn geometry, so it does not care what is on the axes.
            mirror_horizontally(figure.axes[0])
        return figure

    def grid(self, *, backend: str = "plotly"):
        """The bare mesh: cell edges, no values, no basemap, no CRS needed.

        The mesh is fully determined by the grid, so the only argument is which
        renderer draws it. The ``**kwargs`` that used to sit here reached
        ``GridMesh(vor)``, which accepts no other argument, so every one of them
        raised.

        .. generated from myflopy.plot.grid -- edit there, then run scripts/derive_docstrings.py
        Bound form of :func:`myflopy.plot.grid`.

        The one picture :func:`map` cannot give you. A choropleth colours cells
        against a web basemap and so **requires a CRS**; this draws in the grid's own
        coordinates and needs none, which makes it the view for a grid you are still
        refining, before there is a model or a projection.

        Both backends draw the same subject -- this grid -- so ``backend=`` switches
        only the renderer. That is why the 3-D volume is ``grid`` and not
        ``surface``: ``surface`` means a height field.

        .. end generated

        Parameters
        ----------
        backend : {'plotly', 'mpl'}, default 'plotly'
            ``"plotly"`` returns the :class:`GridMesh` picture; ``"mpl"``
            returns FloPy's own patch rendering as a Matplotlib figure -- the
            same thing ``.plot_mpl()`` on the picture gives you, spelled the way
            every other verb spells it. The 3-D volume is
            ``myflopy.plot.grid(stack, backend="vtk")``, which needs a layer
            stack this scope does not have.

        Returns
        -------
        GridMesh or matplotlib.figure.Figure

            .. generated from myflopy.plot.grid -- edit there, then run scripts/derive_docstrings.py
        Raises
        ------
        ValueError
            If ``backend`` is none of ``'plotly'``, ``'mpl'`` or ``'vtk'``, if ``pathlines``
            or any ``vtk only`` argument above is given with the plotly backend, or
            if the VTK backend is asked for without either pathlines or a layer
            stack.

        See Also
        --------
        map : values over the cells, on a basemap.
        surface : a 3-D height field, which is a different shape.
        myflopy.layers.StackPlots.grid : the LAYER-stack 3-D volume, which takes
            ``layers``/``scale``/``color_by``/``cmap`` -- a different builder, and
            not reachable through this function.

            .. end generated

        Examples
        --------
        >>> vor.plot.grid()
        >>> vor.plot.grid(backend="mpl")
        """

        mesh = GridMesh(self.vor)
        if normalize_backend(backend) == "plotly":
            return mesh
        # FloPy's patch renderer hands back an `Axes`, not a figure.
        return as_mpl_figure(mesh.plot_mpl())
