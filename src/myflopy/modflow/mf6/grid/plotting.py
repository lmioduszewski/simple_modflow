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

    The three dependent variables and the two transmissivities get one -- the
    latter labelled in the model's own units, read from ``model``.
    ``type="custom"`` is the package/group path, which supplies its own
    ``custom_hover``; ``"rch"``/``"ks"`` read fields whose hover is assembled
    inline. Returning None for those leaves their existing behaviour exactly as
    it was.
    """

    key = str(map_type).lower()
    if key in ("trans", "trans_max"):
        length, time = unit_labels(model) if model is not None else (None, None)
        return trans_hover(saturated=key == "trans", length=length, time=time)
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
        """
        Plot the grid cross-section with the figs matplotlib cross-section helper.
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
        changed; only the spelling did.
        """

        from flopy.utils.voronoi import VoronoiGrid

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
        bgs: bool = False,
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
        ``per``, ``kstpkper``, ``layer``, ``type``, ``show_mounding``, the
        sectioned-hover knobs -- reads a model's results, and a bare grid has
        none, so advertising them here would promise something this scope cannot
        do. They still ride the ``**trace_kwargs`` tail if you pass one, because
        narrowing the signature should not narrow what already worked.

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
            coordinates with no basemap, which makes ``bgs`` and ``zoom`` inert
            there.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        values : sequence of float, optional
            One value per cell -- heads, drawdown, K, a zone id, a residual, any
            per-cell array. Overrides whatever ``type`` would have read. Length must
            equal ``vor.ncpl``.
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
            Zero-based layer. Layer 0 is the top.
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
            ``model.packages.npf.transmissivity`` / ``.transmissivity_max``, whose
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
            How the per-layer profile renders in the hover.
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
        >>> model.plot.map(values=k, logscale=True, colorscale="Viridis")  # K over decades
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
            bgs=bgs,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            hover=hover,
            select=select,
            select_style=select_style,
            select_color=select_color,
            **({"custom_zs": list(values)} if values is not None else {}),
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

            .. generated from myflopy.plot.section -- edit there, then run scripts/derive_docstrings.py
        cells : int or list of int, optional
            The section path as cell numbers, in order; the line runs through their
            centroids.
        per : int, optional
            Stress period, zero-based: the LAST step saved in it is drawn, so
            ``per=0`` is the end of the first period as the output saved it, on
            both pictures (and on ``map`` and ``surface``). With neither ``per`` nor
            ``kstpkper``, both pictures draw the model's FIRST saved time. Giving
            both raises. A model only.
        kstpkper : tuple of (int, int), optional
            The exact ``(time step, period)`` to draw; one the output did not save
            raises. A model only.
        layers : int or sequence of int, optional
            *(cross-section)* Which layers' CELLS to draw, zero-based; every layer
            by default. The rest are masked out and the vertical extent is cropped
            to what remains. A layer outside the model raises.
        head_layers : int, sequence of int, or None, default 0
            *(cross-section, ``fill='layer'``)* Whose water levels go on top:
            layer 0 by default, a list draws one labelled line each, ``None`` none.
            Ignored by ``'results'`` and an array, which paint the field onto the
            cells instead. A model only.
        layer_labels : sequence of str, optional
            *(cross-section)* Names for the layers, one per layer in model order,
            for the legend and the water levels. Defaults to the names the model was
            built with (``ModelContext(surfaces=stack.build(vor))`` keeps them),
            else ``Layer 1..N``.
        show_grid : bool, default True
            *(cross-section)* Draw the cell edges.
        legend : bool or str, default True
            *(cross-section)* Where the legend goes. ``False`` drops it; ``True``
            puts it at the upper right; ``"auto"`` leaves it to the renderer; a side
            or corner (``"bottom"``, ``"topright"``, or matplotlib's own
            ``"upper right"``) places it inside the plot; ``"outside right"`` or
            ``"outside bottom"`` (also ``left`` and ``top``) places it beside the
            plot -- usually what you want, since inside it sits on the geology it
            describes. An unknown word raises.
        layer : int or list of int, default 0
            *(line profile)* Whose head lines to draw; a list overlays several. An
            interpolated profile draws one layer, so a list of more raises there.
        interpolate : bool, default False
            *(line profile)* Sample an interpolated surface along the line, rather
            than stepping from cell centre to cell centre.
        interpolator : str, optional
            *(line profile)* The interpolation method by name -- ``'griddata'``,
            ``'rbf'``, ``'linearND'`` or ``'cloughTocher2D'``. Giving one turns
            ``interpolate`` on.
        use_rbf : bool, default True
            *(line profile, interpolating)* With no ``interpolator``: radial basis
            functions (``True``) or ``griddata`` (``False``). ``False`` without
            interpolation raises.
        num_points : int, default 100
            *(line profile, interpolating)* How many points to sample along the line.
            Changed without interpolation, it raises: nothing is sampled then.
        clip : Path or geometry, optional
            *(line profile, interpolating)* Clip the interpolated surface to this
            region; raises without interpolation.
        x_or_y : {'x', 'y'}, default 'x'
            *(line profile)* With a single cell in ``cells``, run the line east-west
            (``'x'``) or north-south (``'y'``) through it.
        show_model_top : bool, default True
            *(line profile)* Draw the model top.
        show_model_btm : bool, default False
            *(line profile)* Draw the bottoms of the layers in ``layer``.
        surf_type : {'hds', 'lyr'}, default 'hds'
            *(line profile, interpolating)* Interpolate the head field, or with
            ``'lyr'`` the TOP of ``layer`` -- the model top for ``layer=0``, the bottom
            of layer ``k - 1`` for ``layer=k`` -- from the model's own
            discretization. ``'lyr'`` without interpolation raises.
        animation_kstpkpers : sequence of tuple, optional
            *(line profile, ``backend='plotly'``)* The times ``.ani`` steps through;
            every saved time by default. A Matplotlib figure has no ``.ani``, so
            this raises with ``backend='mpl'``.
        section_name : str, optional
            The cross-section's title, or the line profile's legend name.
        **kwargs
            *(line profile, interpolating)* Passed on to the interpolated surface:
            ``resolution`` (its grid size per axis, default 1000) and ``crs``, or its
            ``xs``/``ys``/``zs``/``vor``/``hds`` inputs. Any of them without
            interpolation raises ``ValueError``; any other keyword raises
            ``TypeError``, on either picture.

            .. end generated

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

            .. generated from myflopy.plot.grid -- edit there, then run scripts/derive_docstrings.py
        pathlines : DataFrame, optional
            Particle track records, drawn as time-coloured tubes over the 3-D mesh.
            Requires ``backend="vtk"``; in plan view the equivalent is
            ``map(pathlines=...)``. Passing it with ``backend="plotly"`` raises.

        vertical_exaggeration : float, default 1.0
            *(vtk only)* Multiplier on z, to make a thin model legible.
        model_style : {'wireframe', 'surface', 'points'}, default 'wireframe'
            *(vtk only)* How the grid itself is drawn beneath the tracks.
        model_opacity : float, default 0.25
            *(vtk only)* Opacity of the grid, so tracks inside it stay visible.
        pathline_cmap : str, default 'viridis'
            *(vtk only)* Colormap for the time-coloured tubes.
        pathline_width : float, default 4.0
            *(vtk only)* Tube width.
        show_edges : bool, default True
            *(vtk only)* Draw cell edges on the mesh.
        off_screen : bool, default True
            *(vtk only)* Render without opening a window -- the right default in a
            notebook or on a headless machine.
        **kwargs
            Forwarded to the backend's builder.

            .. end generated

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
