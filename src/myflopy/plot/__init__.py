"""The plotting front door: one verb per kind of picture (plans 8.3-8.4).

    from myflopy import plot

    plot.map(model, layer=0)                  # plan view
    plot.map(vor, values=drawdown)            # ... of any per-cell array
    plot.section(model, cells=[1653, 651])    # vertical slice through results
    plot.section(vor, line=line)              # ... of the grid itself
    plot.surface(model, layer=0)              # 3-D
    plot.grid(vor)                            # the bare mesh, no basemap
    plot.mosaic([a, b, c])                    # compose any pictures
    plot.animate(frames)                      # flip through pictures

**Geometry chooses the verb.** Not content, and not renderer. A map is a plan
view whatever is drawn on it, so contours, observation locations and a hillshade
are *options* on ``map`` rather than verbs of their own -- which is why there is
no ``plot.contours``. The same rule retires ``plot3d`` (a 3-D view is
``surface``) and ``map_nodes`` (a map whose values are node ids).

``grid`` is the one exception, and it earns it on a constraint rather than a
taste: a choropleth draws over a web basemap and so **requires a CRS**, while the
mesh view needs none. Without it there is no picture at all for a grid you are
still refining, before a projection exists.

**What you pass decides what you get.** Every verb takes a model or a bare grid
as its first argument: a model draws its results, a grid draws itself. Nothing
else changes.

**Everything returned is a Picture** (:class:`myflopy.viz.Picture`): it renders
itself in Jupyter, and answers ``.fig``, ``.show()``, ``.save(path)`` and
``.html(path)``. You never need a trailing ``.plot()``.

Two things deliberately absent:

* **``timeseries``.** Charts belong to a node -- ``model.packages.ghb.results.q
  .plot()`` -- because the series only means something with the model's periods
  attached. There is no useful "chart these bare arrays" that plotly does not
  already do better.
* **Figure construction and theming.** That is :mod:`myflopy.viz` (``Fig``,
  ``Template``, ``PALETTE``, ``subplots``). ``myflopy.plot`` draws models;
  ``myflopy.viz`` builds and styles figures. ``mosaic`` appears in both because
  it is genuinely both -- it is re-exported here, not reimplemented.

**Two of the verbs are COMBINATORS.** ``map``, ``section``, ``surface`` and
``grid`` are picture verbs: their first argument is the subject being drawn.
``mosaic`` and ``animate`` take a collection of finished pictures instead --
they compose rather than draw. They appear on the namespaces too, for
discoverability, but the model is not their subject.

**The same verbs exist on the objects themselves** -- ``model.plot.map()``,
``vor.plot.grid()`` -- and they call these functions, so there is one
implementation behind both spellings. A model answers all six; a bare grid
answers ``map``, ``section`` and ``grid``, which are the three questions it can
answer without results. ``vor.plot()`` is shorthand for ``vor.plot.grid()``.

Binding these replaced ``model.cor()``, ``model.section()``, ``model.srf`` and
eleven ``vor.*`` aliases, and deliberately shadows FloPy's inherited
``VoronoiGrid.plot`` -- whose renderer is still there as
``vor.plot.grid().plot_mpl()``.
"""

from __future__ import annotations

import inspect

import numpy as np

from myflopy._logging import get_logger
from myflopy.modflow.mf6.cross_section_plotting import filled_section
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.plotting import (
    GridMesh,
    GridPlots,
    GridSection,
    _choropleth_factory,
)
from myflopy.modflow.mf6.headsplus import DependentVariableFile
from myflopy.modflow.mf6.interactive_plotting import (
    SliderAnimation,
    build_particle_tracking_scene,
)
from myflopy.modflow.mf6.package_inputs import (
    CellPackageInputFieldExplorer,
    StaticArrayFieldExplorer,
    UzfFieldInputsExplorer,
)
from myflopy.modflow.mf6.package_model import HfbPackageExplorer, HfbResultsExplorer

# The noun tiers live one layer DOWN (package_plotting, L1) so the nouns
# themselves can import them without pointing upward; re-exported here
# because they are a statement about `plot.map`.
from myflopy.modflow.mf6.package_plotting import (
    LAYER_FIELD_MAP_PARAMS,
    MPL_BACKENDS,
    NOUN_INERT_PARAMS,
    NOUN_MAP_PARAMS,
    NOUN_REFUSED_PARAMS,
    _apply_backend,
    as_mpl_figure,
    normalize_backend,
)
from myflopy.modflow.mf6.package_results import (
    CellBudgetResultsExplorer,
    StageResultsExplorer,
)
from myflopy.modflow.mf6.package_surface_water import (
    LakBudgetResultsExplorer,
    LakConnectionsExplorer,
    SfrBudgetResultsExplorer,
    SurfaceWaterExchangeResultsExplorer,
    SurfaceWaterInputFieldExplorer,
)
from myflopy.modflow.mf6.prt_maps import PRTPathlineView
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.xsections import XSection, render_xsections
from myflopy.modflow.utils.inputs import UzfInput
from myflopy.viz import FrameAnimation, Picture, mosaic

logger = get_logger(__name__)

__all__ = [
    "LAYER_FIELD_MAP_PARAMS",
    "NOUN_INERT_PARAMS",
    "NOUN_MAP_PARAMS",
    "NOUN_REFUSED_PARAMS",
    "map",
    "section",
    "surface",
    "grid",
    "mosaic",
    "animate",
    "ModelPlots",
    "FrameAnimation",
    "SliderAnimation",
    "GridPlots",
    "Picture",
    "Choro",
    "XSection",
    "GridSection",
    "GridMesh",
    "InterpolatedSurface",
]


def _grid_of(source):
    """The Voronoi grid for ``source``, and whether ``source`` carries RESULTS.

    Duck-typed rather than isinstance-checked, so loaded runs, live builds and
    group members all work without this module importing three model classes.

    ``.vor`` alone is not enough to mean "model": a :class:`~myflopy.layers
    .LayerStack` and a ``LayerBuildResult`` carry one too, and reading ``.vor``
    as "this is a model" sent them down the results path to die on ``.hds``.
    What actually distinguishes a model is that it HAS results, so that is what
    is asked. A layer stack is geometry -- it answers the grid verbs.
    """

    grid = getattr(source, "vor", None)
    if grid is None:
        # Objects that carry a grid under a different name: an imported
        # `UsgModel` holds `.grid`, and a `BuiltModel` holds it on its context.
        # Without this they were returned AS IF they were grids and died later
        # with `AttributeError: no attribute 'gdf_vorPolys'`, several frames from
        # the call that was actually wrong.
        for owner, attr in ((source, "grid"), (getattr(source, "context", None), "grid")):
            candidate = getattr(owner, attr, None)
            if candidate is not None and hasattr(candidate, "gdf_vorPolys"):
                return candidate, False
    if grid is None:
        # Anything else is taken AS the grid -- this dispatch is duck-typed on
        # purpose (see above), so a stand-in that answers the grid protocol is
        # legitimate and no type check belongs here.
        return source, False
    return grid, _has_results(source)


def _has_results(source) -> bool:
    """Whether ``source`` can serve a results field, without insisting it can.

    This is a CAPABILITY PROBE, so it must not raise. `hasattr` alone is not
    enough: it swallows only `AttributeError`, while `model.hds` on a model that
    has not been run reaches flopy's binary reader and raises `FileNotFoundError`
    -- which escaped this probe and killed every `map()` on an unrun model,
    INCLUDING inputs maps that need no results at all.
    """

    for attr in ("hds", "conc", "temp"):
        try:
            if getattr(source, attr, None) is not None:
                return True
        except (OSError, ValueError, KeyError):  # noqa: BLE001 is not needed --
            # the set IS closed: a missing/short/unreadable output file surfaces
            # as OSError (FileNotFoundError) from flopy's reader, and a present
            # but unparseable one as ValueError/KeyError. Anything else is a real
            # bug and should not be swallowed by a probe.
            logger.debug(
                "results probe: %r is unreadable, treating the source as "
                "input-only", attr,
            )
    return False


def map(      # noqa: A001 - the verb IS `map`
    source,
    /,
    values=None,
    *,
    # -- which numbers ------------------------------------------------------
    per: int | None = None,
    kstpkper: tuple[int, int] | None = None,
    per_timestep: str | int = "last",
    layer: int = 0,
    type: str = "hds",
    # -- colour -------------------------------------------------------------
    zmin: float | None = None,
    zmax: float | None = None,
    colorscale: str | list | tuple | None = None,
    logscale: bool = False,
    # -- contours -----------------------------------------------------------
    contours: bool | str = False,
    contour_values=None,
    contour_levels: int | float | list = 10,
    contour_color: str = "black",
    contour_width: float = 1.5,
    contour_name: str | None = None,
    contour_clip: bool = True,
    contour_resolution: int = 150,
    contour_method: str = "linear",
    # -- selection -----------------------------------------------------------
    select=None,
    select_style: str = "outline",
    select_color: str | None = None,
    # -- overlays and framing ------------------------------------------------
    locs=None,
    hillshade_path=None,
    bgs: bool = False,
    zoom: int = 13,
    fit_bounds: bool = True,
    bounds_padding: float = 0.05,
    # -- hover ---------------------------------------------------------------
    hover=None,
    hover_layers: str | None = None,
    hover_surfaces: bool | None = None,
    hover_fields=None,
    show_layer_elevs: bool | None = None,
    show_mounding: bool = False,
    hover_heads: bool = True,
    hover_ks: bool = False,
    custom_hover: dict | None = None,
    # -- niche ---------------------------------------------------------------
    rch_scale: float | None = None,
    animation_kstpkpers=None,
    # -- renderer ------------------------------------------------------------
    backend: str = "plotly",
    **trace_kwargs,
):
    """Draw a plan-view map of one value per grid cell.

    The single map verb. Contours, observation markers, a hillshade and particle
    pathlines are all **options** here rather than verbs of their own, because a
    map is a plan view whatever is drawn on it.

    Every parameter below is spelled out in the signature rather than swept into
    ``**kwargs``, so an editor can complete and type-check them. The defaults are
    not restated here by hand -- each mirrors the default of whichever link in the
    ``map -> _choropleth_factory -> Choro`` chain owns that argument, and
    ``test_plot_vocabulary`` fails if the two ever drift apart.

    Parameters
    ----------
    source : SimulationBase or VoronoiGridPlus
        A model draws its own results; a bare grid draws the grid. Positional
        only. Passing a grid disables everything that needs model context
        (periods, layer elevations, the sectioned hover).
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
    type : {'hds', 'conc', 'temp', 'rch', 'ks', 'custom'}, default 'hds'
        Which field to read when ``values`` is not given. ``'hds'`` heads,
        ``'conc'`` GWT concentration, ``'temp'`` GWE temperature. The first
        three also select the default sectioned hover.
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
    bgs : bool, default False
        Draw the basemap beneath a semi-transparent cell layer.
    zoom : int, default 13
        Initial map zoom. Ignored when ``fit_bounds`` is True.
    fit_bounds : bool, default True
        Fit the initial view to the grid extent rather than using ``zoom``.
    bounds_padding : float, default 0.05
        Fractional padding around the fitted bounds.
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
        Add head-above-initial (mounding) to the hover.
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
    >>> model.plot.map(layer=0, locs="wells.gpkg", hillshade_path="hs.tif")
    >>> model.plot.map(layer=0, zmin=100, zmax=125).save("heads.png")
    >>> vor.plot.map(values=node_ids)                    # a bare grid
    """

    vor, is_model = _grid_of(source)
    if is_model:
        trace_kwargs.setdefault("model", source)
    if values is not None:
        trace_kwargs["custom_zs"] = list(values)
    # Resolved BEFORE the picture is built, so a misspelled backend raises
    # instead of costing a full render first.
    kind = normalize_backend(backend)
    choro = _choropleth_factory(
        vor,
        per=per,
        kstpkper=kstpkper,
        layer=layer,
        type=type,
        zmin=zmin,
        zmax=zmax,
        colorscale=colorscale,
        logscale=logscale,
        zoom=zoom,
        locs=locs,
        select=select,
        select_style=select_style,
        select_color=select_color,
        custom_hover=custom_hover,
        show_layer_elevs=show_layer_elevs,
        show_mounding=show_mounding,
        hover_heads=hover_heads,
        hover_ks=hover_ks,
        # Not named by the factory; these ride its own **choro_kwargs to `Choro`.
        per_timestep=per_timestep,
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
        hover=hover,
        hover_layers=hover_layers,
        hover_surfaces=hover_surfaces,
        hover_fields=hover_fields,
        rch_scale=rch_scale,
        animation_kstpkpers=animation_kstpkpers,
        **trace_kwargs,
    )
    return _apply_backend(choro, kind)


#: Section arguments that only mean something with results behind them. A bare
#: grid has no periods and no field to interpolate, so `GridSection` accepts none
#: of them -- naming them lets the grid branch say so instead of raising an
#: opaque TypeError from a constructor the caller never named.
def _asked_for(value, default) -> bool:
    """Did the caller actually supply this argument?

    Not simply ``value != default``: ``cells`` and ``layer`` accept sequences,
    and a numpy array compared against ``None`` returns an ELEMENTWISE array,
    which raises "truth value is ambiguous" the moment it is used in a boolean
    context. An array is never a default, so treat that as supplied.
    """

    if value is default:
        return False
    try:
        return bool(value != default)
    except ValueError:
        return True


#: name -> the default that means "not asked for". `GridSection.__init__` takes
#: `(vor, line)` and nothing else, so every one of these is model-only.
_MODEL_ONLY_SECTION_ARGS = {
    "cells": None,
    "per": None,
    "kstpkper": None,
    "layer": 0,
    "x_or_y": "x",
    "spacing": 10,
    "num_points": 100,
    "interpolate": False,
    "use_rbf": True,
    "interpolator": None,
    "extrapolate_beyond_section_ends": False,
    "show_model_top": True,
    "show_model_btm": False,
    "surf_type": "hds",
    "section_name": None,
    "clip": None,
    "animation_kstpkpers": None,
}


#: Section arguments that shape the LINE PROFILE and have no meaning once the
#: cells themselves are drawn. `fill=` renders through FloPy's
#: `PlotCrossSection`, which walks the grid cell by cell -- there is no sampling
#: interval to set, no interpolation to choose, and the model top is an edge of
#: the drawn geometry rather than a separate trace. Naming them lets the filled
#: branch say so instead of accepting a knob it will not turn.
_PROFILE_ONLY_SECTION_ARGS = {
    "x_or_y": "x",
    "spacing": 10,
    "num_points": 100,
    "interpolate": False,
    "use_rbf": True,
    "interpolator": None,
    "extrapolate_beyond_section_ends": False,
    "show_model_top": True,
    "show_model_btm": False,
    "surf_type": "hds",
    "clip": None,
    "animation_kstpkpers": None,
}



def _filled_section(source, vor, is_model, *, line, cells, per, kstpkper,
                    fill, fill_cmap, fill_label, layers, head_layers,
                    layer_labels, title, backend, asked):
    """Draw the section as CELLS -- the grid, its layers, and a field on them.

    The picture `section()` could not draw before 2026-09-02: its line profile
    answers "what is the head along this line", and this answers "what does the
    model look like in cross-section". Same subject, same verb; `fill` chooses
    which. This resolves the verb's arguments; the drawing is
    :func:`~myflopy.modflow.mf6.cross_section_plotting.filled_section`, which the
    grid scope calls directly.
    """

    if backend != "mpl":
        raise ValueError(
            "fill= draws the cells through FloPy's cross-section renderer, "
            "which is Matplotlib; pass backend='mpl'. Without fill=, the "
            "section is a line profile and both backends draw it."
        )
    stray = sorted(
        name for name, default in _PROFILE_ONLY_SECTION_ARGS.items()
        if _asked_for(asked.get(name, default), default)
    )
    if stray:
        raise ValueError(
            f"{', '.join(stray)} shape the line PROFILE; fill= draws the cells "
            f"themselves, so they have nothing to act on. Drop them, or drop "
            f"fill= to get the profile."
        )
    if _asked_for(asked.get("layer", 0), 0):
        raise ValueError(
            "layer= overlays head PROFILES, which the filled section does not "
            "draw. Use layers= to choose which layers' cells to draw, and "
            "head_layers= to choose whose water levels go on top."
        )
    if line is None and cells is None:
        raise ValueError("fill= needs a section line: pass line= or cells=.")
    if line is not None and cells is not None:
        raise ValueError("line= and cells= both give the path; pass one.")

    if cells is not None:
        # A cell path becomes the polyline through those cells' centroids, so
        # `cells=` keeps working here rather than being refused for a reason the
        # caller would find arbitrary.
        chosen = [cells] if isinstance(cells, (int, np.integer)) else list(cells)
        centroids = vor.gdf_vorPolys.geometry.iloc[chosen].centroid
        path = [(point.x, point.y) for point in centroids]
    else:
        path = line

    # A bare grid has no flopy modelgrid of its own; `GridSection` already builds
    # one from the DISV gridprops, so borrow it rather than repeat the
    # construction (top/botm/cell2d/ncpl/crs) a second time.
    grid = source.gwf.modelgrid if is_model else GridSection(vor=vor, line=_as_line(path)).grid
    return filled_section(
        grid, path,
        model=source if is_model else None,
        fill=fill, layers=layers, head_layers=head_layers,
        layer_labels=layer_labels, per=per, kstpkper=kstpkper,
        cmap=fill_cmap, label=fill_label, title=title,
    )


def _as_line(path):
    """A LineString for `GridSection`, which needs a geometry rather than coords."""

    import shapely as shp

    return path if isinstance(path, shp.LineString) else shp.LineString(list(path))


def section(
    source,
    /,
    line=None,
    *,
    cells=None,
    per: int | None = None,
    kstpkper: tuple[int, int] | None = None,
    layer: int | list[int] = 0,
    x_or_y: str = "x",
    spacing: int = 10,
    num_points: int = 100,
    interpolate: bool = False,
    use_rbf: bool = True,
    interpolator: str | None = None,
    extrapolate_beyond_section_ends: bool = False,
    show_model_top: bool = True,
    show_model_btm: bool = False,
    surf_type: str = "hds",
    section_name: str | None = None,
    clip=None,
    animation_kstpkpers=None,
    fill=None,
    fill_cmap: str = "viridis",
    fill_label: str | None = None,
    layers=None,
    head_layers: int | list[int] | None = 0,
    layer_labels=None,
    backend: str = "plotly",
    **kwargs,
):
    """Draw a vertical slice through a model's results, or a grid's geometry.

    Through a MODEL this is the results section -- the field against distance
    along the line. Through a bare GRID it is the geometry section: layers and
    cell edges, no results. They are different classes because they answer
    different questions, and **most arguments below only apply to the model
    branch** -- they are marked. A grid takes ``line`` and nothing else.

    Parameters
    ----------
    source : SimulationBase or VoronoiGridPlus
        Positional only. A model gives results; a grid gives geometry.
    line : LineString or Path, optional
        The section line as a geometry or a vector file. Positional, so
        ``plot.section(vor, line)`` reads naturally. For a bare grid this is the
        only way to specify the path. Mutually exclusive with ``cells``.
    cells : int or list of int, optional
        *(model only)* Cell indices defining the section path, in order.
    per : int, optional
        *(model only)* Stress period (0-based). Mutually exclusive with
        ``kstpkper``.
    kstpkper : tuple of (int, int), optional
        *(model only)* Exact ``(timestep, period)``.
    layer : int or list of int, default 0
        *(model only)* Layer(s) to draw. A list overlays several.
    x_or_y : {'x', 'y'}, default 'x'
        Which coordinate becomes the horizontal axis.
    spacing : int, default 10
        *(model only)* Sample spacing along the line, in model units.
    num_points : int, default 100
        *(model only)* Number of samples when interpolating.
    interpolate : bool, default False
        *(model only)* Interpolate between cell centers rather than stepping
        cell to cell.
    use_rbf : bool, default True
        *(model only)* Use radial-basis interpolation when ``interpolate`` is
        True.
    interpolator : str, optional
        *(model only)* Override the interpolation method by name.
    extrapolate_beyond_section_ends : bool, default False
        *(model only)* Extend the section past the first and last cell centers.
    show_model_top : bool, default True
        *(model only)* Draw the model-top profile.
    show_model_btm : bool, default False
        *(model only)* Draw layer-bottom profiles.
    surf_type : {'hds', 'lyr'}, default 'hds'
        *(model only)* Section the head field, or the layer elevations.
    section_name : str, optional
        *(model only)* Legend name for the traces.
    clip : Path or geometry, optional
        *(model only)* Restrict the section to cells intersecting this region.
    animation_kstpkpers : sequence of tuple, optional
        *(model only)* The periods ``.ani`` steps through; defaults to every
        output time.
    fill : {'layer', 'results'} or array-like, optional
        Draw the CELLS instead of a line profile -- the grid, its layers, and a
        field painted on them. This is the "cross-section of the model" picture:
        without it, ``section`` answers *what is the head along this line*; with
        it, *what does the model look like through here*.

        ``'layer'`` colours each cell by its layer and draws the simulated head
        as a water surface over the geology. ``'results'`` colours the cells by
        the model's own dependent variable instead -- heads on GWF,
        concentration on GWT, temperature on GWE -- with a colorbar. Any
        ``(nlay, ncpl)`` array does the same for a field of your own (K, a zone
        id, a residual); a flat ``(ncpl,)`` array raises rather than colouring
        every layer the same.

        **Requires ``backend="mpl"``**: the renderer is FloPy's
        ``PlotCrossSection`` and there is no Plotly equivalent to defer to.
        Passing ``fill=`` with the default backend raises and says so, rather
        than returning a line profile you did not ask for.

        The line-profile arguments -- ``interpolate``, ``spacing``,
        ``num_points``, ``use_rbf``, ``interpolator``, ``x_or_y``,
        ``show_model_top``, ``show_model_btm``, ``surf_type``, ``clip``,
        ``extrapolate_beyond_section_ends``, ``animation_kstpkpers`` -- have
        nothing to act on here and raise if given. ``line=`` and ``cells=`` both
        work; a cell path becomes the polyline through those cells' centroids.
    fill_cmap : str, default 'viridis'
        Colormap for ``fill='results'`` or an array ``fill``. Ignored for
        ``fill='layer'``, which uses the discrete layer palette.
    fill_label : str, optional
        Colorbar label for ``fill='results'`` or an array ``fill``.
    layers : int or sequence of int, optional
        *(``fill=`` only)* Which layers' CELLS to draw, zero-based. With none,
        every layer. The unselected cells are masked out and the vertical extent
        is cropped to what remains -- leaving the axis at full height would put
        the two layers you asked for in a thin band with empty space above and
        below. A layer outside the model raises rather than drawing nothing.

        This is why ``layer=`` raises here: it means something else. ``layer=``
        overlays head PROFILES on the line-profile section; the filled section
        draws cells, so choosing them needs its own name.
    head_layers : int or sequence of int or None, default 0
        *(``fill='layer'`` only)* Whose water levels to draw over the geology.
        The default, layer 0, is the single water table every filled section drew
        before this parameter existed. A list draws one surface per layer, each
        labelled and coloured through :func:`~myflopy.viz.category_colors` -- so
        a given layer's water level keeps its colour across figures. ``None``
        draws none, for the geology alone.

        Ignored with ``fill='results'`` or an array fill: those paint the field
        onto the cells, so a line of the same quantity on top would say it twice.
    layer_labels : sequence of str, optional
        *(``fill=`` only)* Legend names for the layers -- your unit names rather
        than ``Layer 1..N``. One per layer, in model order.

        Defaults to the names the model's own build context carries
        (``ModelContext(surfaces=stack.build(vor))`` keeps them, so a model
        declared through the spec API already knows its layers are called
        "sand" and "clay"). A model built imperatively, or one whose ``surfaces``
        is a plain frame, has no names to find and falls back to ``Layer N``.
    backend : {'plotly', 'mpl'}, default 'plotly'
        Which renderer draws the section. ``'plotly'`` returns the interactive
        picture; ``'mpl'`` returns a static
        :class:`matplotlib.figure.Figure`. Accepts ``'interactive'`` and
        ``'matplotlib'``/``'static'`` as aliases; anything else raises.

        Applies to both branches: a model section renders through the same
        overlay builder the noun grammar uses
        (``model.hds.section(line=..., backend='mpl')``), and a grid section
        through ``GridSection``'s own cell-outline renderer. Until 2026-09-02
        this parameter did not exist on the verb and a caller who passed it got
        a Plotly figure back with no error at all -- it vanished into the
        ``**kwargs`` tail. Passing it is now the documented spelling; prefer it
        to ``.plot_mpl()`` on the result, which is the same renderer reached a
        second way.
    **kwargs
        Forwarded to the underlying section class.

    Returns
    -------
    XSection or GridSection or matplotlib.figure.Figure
        With ``backend='plotly'`` (the default) a
        :class:`~myflopy.viz.Picture` -- ``XSection`` for a model, which also
        carries ``.ani``, or ``GridSection`` for a bare grid. Both answer
        ``.plot_mpl()``.

        With ``backend='mpl'``, a bare Matplotlib ``Figure``: use
        ``.savefig(path)`` and ``.axes[0]``, not the picture verbs.

    Raises
    ------
    ValueError
        If a model-only argument is given for a bare grid. It names the
        arguments, because the alternative is a ``TypeError`` from a constructor
        the caller never mentioned.

    See Also
    --------
    map : the same data in plan view.
    myflopy.plot.grid : the mesh itself.

    Examples
    --------
    >>> model.plot.section(cells=[1653, 651, 1241])
    >>> model.plot.section(line=line, layer=[0, 1], interpolate=True)
    >>> model.plot.section(cells=cells).save("section.png")
    >>> vor.plot.section(line)                       # geometry, no results
    >>> model.plot.section(line=line, backend="mpl") # a static mpl Figure
    >>> vor.plot.section(line, backend="mpl").savefig("grid_section.png")
    >>> model.plot.section(line=line, fill="layer", backend="mpl")
    >>> model.plot.section(line=line, fill="results", per=15, backend="mpl")
    >>> vor.plot.section(line, fill="layer", backend="mpl")   # geometry only
    >>> model.plot.section(line=line, fill="layer", layers=[0, 1], backend="mpl")
    >>> model.plot.section(line=line, fill="layer", head_layers=[0, 2], backend="mpl")
    >>> model.plot.section(line=line, fill="layer", backend="mpl",
    ...                    layer_labels=["fill", "sand", "clay", "till"])
    """

    vor, is_model = _grid_of(source)
    kind = normalize_backend(backend)
    if fill is not None:
        return _filled_section(
            source, vor, is_model,
            line=line, cells=cells, per=per, kstpkper=kstpkper,
            fill=fill, fill_cmap=fill_cmap, fill_label=fill_label,
            layers=layers, head_layers=head_layers, layer_labels=layer_labels,
            title=section_name, backend=kind,
            asked=locals(),
        )
    if layers is not None or layer_labels is not None or head_layers != 0:
        raise ValueError(
            "layers=, head_layers= and layer_labels= describe the CELLS a "
            "filled section draws; the line profile has none. Pass fill= "
            "(with backend='mpl'), or use layer= to overlay profiles."
        )
    if is_model:
        picture = XSection(
            model=source,
            line=line,
            cells=cells,
            per=per,
            kstpkper=kstpkper,
            layer=layer,
            x_or_y=x_or_y,
            spacing=spacing,
            num_points=num_points,
            interpolate=interpolate,
            use_rbf=use_rbf,
            interpolator=interpolator,
            extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
            show_model_top=show_model_top,
            show_model_btm=show_model_btm,
            surf_type=surf_type,
            section_name=section_name,
            clip=clip,
            animation_kstpkpers=animation_kstpkpers,
            **kwargs,
        )
        if kind == "plotly":
            return picture
        # The SAME renderer the noun grammar uses (`model.hds.section(...)`),
        # not a second matplotlib path -- one overlay builder, so the two
        # spellings cannot drift into drawing different pictures.
        return render_xsections(
            {picture.section_name: picture}, backend="mpl", title=section_name,
        )

    asked = {
        "cells": cells, "per": per, "kstpkper": kstpkper, "layer": layer,
        "x_or_y": x_or_y, "spacing": spacing, "num_points": num_points,
        "interpolate": interpolate, "use_rbf": use_rbf,
        "interpolator": interpolator,
        "extrapolate_beyond_section_ends": extrapolate_beyond_section_ends,
        "show_model_top": show_model_top, "show_model_btm": show_model_btm,
        "surf_type": surf_type, "section_name": section_name, "clip": clip,
        "animation_kstpkpers": animation_kstpkpers,
    }
    given = [
        name for name, value in asked.items()
        if _asked_for(value, _MODEL_ONLY_SECTION_ARGS[name])
    ]
    if given:
        raise ValueError(
            f"{', '.join(given)} need a model's results; a bare grid section "
            f"draws geometry only. Pass a model as the first argument, or drop "
            f"these and give just the line."
        )
    geometry = GridSection(vor=vor, line=line, **kwargs)
    if kind == "plotly":
        return geometry
    return as_mpl_figure(geometry.plot_mpl())


def surface(
    source,
    /,
    *,
    layer: int = 0,
    per: int | None = None,
    kstpkper: tuple[int, int] | None = None,
    surf_type: str = "hds",
    resolution: int = 1000,
    use_rbf: bool = False,
    interpolator: str | None = None,
    clip=None,
    crs: str | None = None,
    xs=None,
    ys=None,
    zs=None,
    **kwargs,
) -> InterpolatedSurface:
    """Draw a 3-D interpolated surface -- a height field ``z(x, y)``.

    Always Plotly, and always a height field. The 3-D grid VOLUME and 3-D
    particle pathlines are different shapes and live on :func:`grid` with
    ``backend="vtk"``: a ``backend=`` switch should change the renderer, not
    what is being drawn.

    Parameters
    ----------
    source : SimulationBase or VoronoiGridPlus
        Positional only.
    layer : int, default 0
        Zero-based layer to interpolate.
    per : int, optional
        Stress period (0-based). Mutually exclusive with ``kstpkper``.
    kstpkper : tuple of (int, int), optional
        Exact ``(timestep, period)``. Defaults to the model's FIRST output time.
    surf_type : {'hds', 'lyr'}, default 'hds'
        Interpolate the head field, or the layer elevation surface.
    resolution : int, default 1000
        Interpolation grid size per axis. Higher is smoother and slower.
    use_rbf : bool, default False
        Radial-basis interpolation instead of linear.
    interpolator : str, optional
        Override the interpolation method by name.
    clip : Path or geometry, optional
        Restrict the surface to cells intersecting this region. See also
        ``.clipped_fig()`` on the result.
    crs : str, optional
        Override the grid's CRS.
    xs, ys, zs : array-like, optional
        Supply the point cloud directly instead of reading it from a model.
    **kwargs
        Forwarded to :class:`~myflopy.modflow.mf6.grid.interpolated_surface
        .InterpolatedSurface`.

    Returns
    -------
    InterpolatedSurface
        A :class:`~myflopy.viz.Picture`; also carries ``.clipped_fig()`` and
        ``.surface_trace()`` for composing into a larger scene.

    See Also
    --------
    grid : the mesh itself, including the 3-D volume via ``backend="vtk"``.

    Examples
    --------
    >>> model.plot.surface(layer=0)
    >>> model.plot.surface(layer=0, surf_type="lyr")     # layer elevations
    >>> model.plot.surface(layer=0, resolution=400).html("surface.html")
    """

    vor, is_model = _grid_of(source)
    common = dict(
        layer=layer,
        per=per,
        kstpkper=kstpkper,
        surf_type=surf_type,
        resolution=resolution,
        use_rbf=use_rbf,
        interpolator=interpolator,
        clip=clip,
        crs=crs,
        xs=xs,
        ys=ys,
        zs=zs,
    )
    if is_model:
        return InterpolatedSurface(model=source, **common, **kwargs)
    return InterpolatedSurface(vor=vor, **common, **kwargs)


#: name -> "not asked for" default, for the 3-D scene builder. Same purpose as
#: `_MODEL_ONLY_SECTION_ARGS`: the plotly branch draws a flat mesh and accepts
#: none of these, so silently dropping them would be worse than saying so.
_VTK_ONLY_GRID_ARGS = {
    "vertical_exaggeration": 1.0,
    "model_style": "wireframe",
    "model_opacity": 0.25,
    "pathline_cmap": "viridis",
    "pathline_width": 4.0,
    "show_edges": True,
    "off_screen": True,
}


def grid(
    source,
    /,
    *,
    backend: str = "plotly",
    pathlines=None,
    vertical_exaggeration: float = 1.0,
    model_style: str = "wireframe",
    model_opacity: float = 0.25,
    pathline_cmap: str = "viridis",
    pathline_width: float = 4.0,
    show_edges: bool = True,
    off_screen: bool = True,
    **kwargs,
):
    """Draw the mesh itself -- flat in 2-D, or the layered volume in 3-D.

    The one picture :func:`map` cannot give you. A choropleth colours cells
    against a web basemap and so **requires a CRS**; this draws in the grid's own
    coordinates and needs none, which makes it the view for a grid you are still
    refining, before there is a model or a projection.

    Both backends draw the same subject -- this grid -- so ``backend=`` switches
    only the renderer. That is why the 3-D volume is ``grid`` and not
    ``surface``: ``surface`` means a height field.

    Parameters
    ----------
    source : SimulationBase or VoronoiGridPlus or LayerBuildResult
        Positional only. Anything that carries a grid.
    backend : {'plotly', 'mpl', 'vtk'}, default 'plotly'
        ``'plotly'`` draws cell edges in 2-D -- fast, CRS-free, no basemap.
        ``'mpl'`` draws the same 2-D mesh through FloPy's own patch renderer and
        returns a static :class:`matplotlib.figure.Figure`; it is the same
        renderer as ``.plot_mpl()`` on the returned picture, offered here so the
        switch is spelled the same way on every verb. ``'vtk'`` renders the cell
        VOLUME in 3-D as an interactive PyVista scene, and needs the ``viz3d``
        extra (``pip install 'myflopy[viz3d]'``).

        All three draw the same subject -- this grid -- which is the rule the
        parameter follows everywhere: ``backend=`` changes the renderer, never
        what is being drawn.
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

    Returns
    -------
    GridMesh or VtkScene
        Both are :class:`~myflopy.viz.Picture`. ``GridMesh`` carries
        ``.plot_mpl()``, which is FloPy's own patch renderer. ``VtkScene``
        exposes ``.scene`` (the PyVista ``Plotter``) instead of ``.fig``.

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

    Examples
    --------
    >>> vor.plot.grid()                         # the 2-D mesh
    >>> vor.plot()                              # the same thing, shorthand
    >>> vor.plot.grid().plot_mpl()              # FloPy's matplotlib renderer
    >>> stack.plot.grid(["sand", "clay"])       # 3-D volume, a subset of layers
    >>> model.plot.grid(pathlines=run.track_records, backend="vtk")
    """

    scene_args = {
        "vertical_exaggeration": vertical_exaggeration,
        "model_style": model_style,
        "model_opacity": model_opacity,
        "pathline_cmap": pathline_cmap,
        "pathline_width": pathline_width,
        "show_edges": show_edges,
        "off_screen": off_screen,
    }

    if str(backend).lower() in MPL_BACKENDS:
        if pathlines is not None:
            raise ValueError(
                "pathlines are only drawn by the 3-D scene; pass backend='vtk' "
                "for tubes over the grid, or use map(pathlines=...) in plan view."
            )
        vor, _ = _grid_of(source)
        # FloPy's patch renderer returns an `Axes`; the contract is a figure.
        return as_mpl_figure(GridMesh(vor, **kwargs).plot_mpl())
    if backend == "plotly":
        if pathlines is not None:
            raise ValueError(
                "pathlines are only drawn by the 3-D scene; pass backend='vtk' "
                "for tubes over the grid, or use map(pathlines=...) in plan view."
            )
        # Named but not accepted by the flat mesh. Forwarding them would raise a
        # TypeError from `GridMesh`; dropping them silently would be worse still.
        stray = [n for n, v in scene_args.items() if _asked_for(v, _VTK_ONLY_GRID_ARGS[n])]
        if stray:
            raise ValueError(
                f"{', '.join(stray)} configure the 3-D scene; the flat mesh has "
                f"no use for them. Pass backend='vtk', or drop them."
            )
        vor, _ = _grid_of(source)
        return GridMesh(vor, **kwargs)
    if backend != "vtk":
        raise ValueError(
            f"backend must be 'plotly', 'mpl' or 'vtk', not {backend!r}."
        )

    if pathlines is None:
        raise ValueError(
            "plot.grid(backend='vtk') needs either pathlines= (a model's particle "
            "tracks) or a layer stack -- use stack.plot.grid() for layer geometry."
        )
    return build_particle_tracking_scene(source, pathlines, **scene_args, **kwargs)


def animate(
    frames,
    *,
    backend: str = "plotly",
    title=None,
    dpi: int = 140,
    interval_ms: int = 700,
    **kwargs,
):
    """Flip through a sequence of pictures.

    A **combinator**, like :func:`mosaic` -- its first argument is the frames,
    not a subject. To animate one model's results over time, prefer the grammar
    (``model.hds.animate(kind="map", over="period")``), which generates the
    frames for you and calls this.

    Parameters
    ----------
    frames : sequence
        The frames, as bare pictures or ``(label, picture)`` pairs -- the same
        shapes :func:`mosaic` accepts. Unlabelled frames are numbered.
    backend : {'plotly', 'png'}, default 'plotly'
        ``'plotly'`` builds one live figure with play/pause and a slider: fast
        and interactive, but every frame must share a trace structure, and a
        choropleth re-embeds its geometry per frame, so the file grows with
        cells x frames. ``'png'`` rasterizes each frame and pages through them
        with a browser slider: frames need share NOTHING, so kinds can be mixed,
        and the size does not grow with cell count. Prefer ``'png'`` on a large
        grid.
    title : str, optional
        Figure title (plotly) or document title (png).
    dpi : int, default 140
        *(png only)* Raster resolution per frame.
    interval_ms : int, default 700
        *(png only)* Milliseconds per frame during playback.
    **kwargs
        Forwarded to the backend's animation class.

    Returns
    -------
    FrameAnimation or SliderAnimation
        Both are :class:`~myflopy.viz.Picture`. ``FrameAnimation.fig`` is the
        plotly figure; ``SliderAnimation`` has no single figure, so its ``.fig``
        raises and names ``.frames`` instead. ``SliderAnimation.export()``
        returns the richer ``StandaloneHtmlSlider`` handle when you want the
        frame manifest or the resume/progress machinery.

    Raises
    ------
    ValueError
        If ``frames`` is empty, if ``backend`` is neither value, or if the
        plotly backend is given frames that do not share a trace structure. That
        last one raises rather than silently falling back to raster -- swapping
        an interactive figure for a static page changes what you get.

    See Also
    --------
    myflopy.viz.mosaic : the same frames side by side instead of in sequence.
    myflopy.export_head_map_slider_html : a FloPy-rendered slider, which is a
        different picture rather than a second spelling of this one.

    Examples
    --------
    >>> frames = [model.plot.map(per=p, layer=0) for p in range(model.nper)]
    >>> plot.animate(frames).show()
    >>> plot.animate(frames, backend="png").html("heads.html")
    >>> plot.animate([("start", first), ("end", last)])
    >>> mixed = [("a map", model.plot.map()), ("a section", model.plot.section(cells=cells))]
    >>> plot.animate(mixed, backend="png")           # only png can mix kinds
    """

    if backend == "plotly":
        # A live plotly figure has no raster step, so these two mean nothing to
        # it. `FrameAnimation` would raise TypeError naming an argument the
        # caller never wrote; say which backend they belong to instead.
        stray = [
            n for n, v, default in
            (("dpi", dpi, 140), ("interval_ms", interval_ms, 700))
            if v != default
        ]
        if stray:
            raise ValueError(
                f"{', '.join(stray)} configure the rasterised slider; the live "
                f"plotly figure has no frames to rasterise. Pass backend='png'."
            )
        return FrameAnimation(frames, title=title, **kwargs)
    if backend == "png":
        return SliderAnimation(
            frames, title=title, dpi=dpi, interval_ms=interval_ms, **kwargs
        )
    raise ValueError(f"backend must be 'plotly' or 'png', not {backend!r}.")


def _bound_args(local_vars: dict, tail: str) -> dict:
    """The named arguments of a bound verb, ready to forward to the free one.

    Called as the FIRST statement of each bound verb, so ``locals()`` holds
    exactly ``self``, the named parameters, and the ``**kwargs`` tail -- nothing
    else has been assigned yet.

    Why not retype the names in the call? Because each bound verb mirrors thirty-
    odd parameters of its free counterpart, and a forwarding list written out by
    hand is a second place to forget one. The signature stays explicit (that is
    the whole point -- an editor reads the ``def``), while the body cannot drift
    from it.
    """

    return {k: v for k, v in local_vars.items() if k not in ("self", tail)}


class ModelPlots:
    """The verbs bound to one model -- ``model.plot.map()`` (plan 8.4).

    Exactly the same functions as the module-level front door, with the model
    already supplied: ``model.plot.map(layer=0)`` is ``plot.map(model, layer=0)``.
    Which spelling you reach for is taste; there is one implementation.

    Bound rather than re-dispatched: the module verbs duck-type on ``.vor`` to
    tell a model from a grid, and a model whose grid has not been resolved yet
    carries ``vor = None`` -- which would silently make it look like a grid. The
    namespace knows what it holds, so it passes the model straight through.

    A fresh instance per access, like ``model.packages``: models are rehydrated
    via ``cls.__new__`` in ``from_built_run``, so a memoized namespace could
    outlive the state it closed over.

    Each verb repeats its free counterpart's parameters in full rather than
    taking ``**kwargs``. That duplication is deliberate: PyCharm and Pylance are
    static, so a ``**kwargs`` forwarder shows the caller nothing, no matter how
    good the docstring is. ``test_plot_vocabulary`` asserts each bound signature
    equals the free one minus ``source``, so the copies cannot drift.
    """

    def __init__(self, model):
        """Bind the plotting verbs to ``model``."""

        self.model = model

    def __repr__(self):
        """Name the subject and the verbs, so tab-completion has a companion."""

        return f"ModelPlots({getattr(self.model, 'name', '?')!r}: map, section, surface, grid, animate, mosaic)"

    # The bare names below resolve to this module's functions, not to these
    # methods -- class scope is not in a method's name-lookup chain.
    def map(      # noqa: A003 - the verb IS `map`
        self,
        values=None,
        *,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        per_timestep: str | int = "last",
        layer: int = 0,
        type: str = "hds",
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
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        locs=None,
        hillshade_path=None,
        bgs: bool = False,
        zoom: int = 13,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        hover=None,
        hover_layers: str | None = None,
        hover_surfaces: bool | None = None,
        hover_fields=None,
        show_layer_elevs: bool | None = None,
        show_mounding: bool = False,
        hover_heads: bool = True,
        hover_ks: bool = False,
        custom_hover: dict | None = None,
        rch_scale: float | None = None,
        animation_kstpkpers=None,
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """This model's plan-view map. See :func:`myflopy.plot.map`.

        The single map verb. Contours, observation markers, a hillshade and particle
        pathlines are all **options** here rather than verbs of their own, because a
        map is a plan view whatever is drawn on it.

        Every parameter below is spelled out in the signature rather than swept into
        ``**kwargs``, so an editor can complete and type-check them. The defaults are
        not restated here by hand -- each mirrors the default of whichever link in the
        ``map -> _choropleth_factory -> Choro`` chain owns that argument, and
        ``test_plot_vocabulary`` fails if the two ever drift apart.

        Parameters
        ----------
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
        type : {'hds', 'conc', 'temp', 'rch', 'ks', 'custom'}, default 'hds'
            Which field to read when ``values`` is not given. ``'hds'`` heads,
            ``'conc'`` GWT concentration, ``'temp'`` GWE temperature. The first
            three also select the default sectioned hover.
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
        bgs : bool, default False
            Draw the basemap beneath a semi-transparent cell layer.
        zoom : int, default 13
            Initial map zoom. Ignored when ``fit_bounds`` is True.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
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
            Add head-above-initial (mounding) to the hover.
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
        >>> model.plot.map(layer=0, locs="wells.gpkg", hillshade_path="hs.tif")
        >>> model.plot.map(layer=0, zmin=100, zmax=125).save("heads.png")
        >>> vor.plot.map(values=node_ids)                    # a bare grid

        Bound form of :func:`myflopy.plot.map`.
        """

        return map(self.model, **_bound_args(locals(), "trace_kwargs"), **trace_kwargs)

    def section(
        self,
        line=None,
        *,
        cells=None,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer: int | list[int] = 0,
        x_or_y: str = "x",
        spacing: int = 10,
        num_points: int = 100,
        interpolate: bool = False,
        use_rbf: bool = True,
        interpolator: str | None = None,
        extrapolate_beyond_section_ends: bool = False,
        show_model_top: bool = True,
        show_model_btm: bool = False,
        surf_type: str = "hds",
        section_name: str | None = None,
        clip=None,
        animation_kstpkpers=None,
        fill=None,
        fill_cmap: str = "viridis",
        fill_label: str | None = None,
        layers=None,
        head_layers: int | list[int] | None = 0,
        layer_labels=None,
        backend: str = "plotly",
        **kwargs,
    ):
        """A vertical slice through this model. See :func:`myflopy.plot.section`.

        Through a MODEL this is the results section -- the field against distance
        along the line. Through a bare GRID it is the geometry section: layers and
        cell edges, no results. They are different classes because they answer
        different questions, and **most arguments below only apply to the model
        branch** -- they are marked. A grid takes ``line`` and nothing else.

        Parameters
        ----------
        line : LineString or Path, optional
            The section line as a geometry or a vector file. Positional, so
            ``plot.section(vor, line)`` reads naturally. For a bare grid this is the
            only way to specify the path. Mutually exclusive with ``cells``.
        cells : int or list of int, optional
            *(model only)* Cell indices defining the section path, in order.
        per : int, optional
            *(model only)* Stress period (0-based). Mutually exclusive with
            ``kstpkper``.
        kstpkper : tuple of (int, int), optional
            *(model only)* Exact ``(timestep, period)``.
        layer : int or list of int, default 0
            *(model only)* Layer(s) to draw. A list overlays several.
        x_or_y : {'x', 'y'}, default 'x'
            Which coordinate becomes the horizontal axis.
        spacing : int, default 10
            *(model only)* Sample spacing along the line, in model units.
        num_points : int, default 100
            *(model only)* Number of samples when interpolating.
        interpolate : bool, default False
            *(model only)* Interpolate between cell centers rather than stepping
            cell to cell.
        use_rbf : bool, default True
            *(model only)* Use radial-basis interpolation when ``interpolate`` is
            True.
        interpolator : str, optional
            *(model only)* Override the interpolation method by name.
        extrapolate_beyond_section_ends : bool, default False
            *(model only)* Extend the section past the first and last cell centers.
        show_model_top : bool, default True
            *(model only)* Draw the model-top profile.
        show_model_btm : bool, default False
            *(model only)* Draw layer-bottom profiles.
        surf_type : {'hds', 'lyr'}, default 'hds'
            *(model only)* Section the head field, or the layer elevations.
        section_name : str, optional
            *(model only)* Legend name for the traces.
        clip : Path or geometry, optional
            *(model only)* Restrict the section to cells intersecting this region.
        animation_kstpkpers : sequence of tuple, optional
            *(model only)* The periods ``.ani`` steps through; defaults to every
            output time.
        fill : {'layer', 'results'} or array-like, optional
            Draw the CELLS instead of a line profile -- the grid, its layers, and a
            field painted on them. This is the "cross-section of the model" picture:
            without it, ``section`` answers *what is the head along this line*; with
            it, *what does the model look like through here*.

            ``'layer'`` colours each cell by its layer and draws the simulated head
            as a water surface over the geology. ``'results'`` colours the cells by
            the model's own dependent variable instead -- heads on GWF,
            concentration on GWT, temperature on GWE -- with a colorbar. Any
            ``(nlay, ncpl)`` array does the same for a field of your own (K, a zone
            id, a residual); a flat ``(ncpl,)`` array raises rather than colouring
            every layer the same.

            **Requires ``backend="mpl"``**: the renderer is FloPy's
            ``PlotCrossSection`` and there is no Plotly equivalent to defer to.
            Passing ``fill=`` with the default backend raises and says so, rather
            than returning a line profile you did not ask for.

            The line-profile arguments -- ``interpolate``, ``spacing``,
            ``num_points``, ``use_rbf``, ``interpolator``, ``x_or_y``,
            ``show_model_top``, ``show_model_btm``, ``surf_type``, ``clip``,
            ``extrapolate_beyond_section_ends``, ``animation_kstpkpers`` -- have
            nothing to act on here and raise if given. ``line=`` and ``cells=`` both
            work; a cell path becomes the polyline through those cells' centroids.
        fill_cmap : str, default 'viridis'
            Colormap for ``fill='results'`` or an array ``fill``. Ignored for
            ``fill='layer'``, which uses the discrete layer palette.
        fill_label : str, optional
            Colorbar label for ``fill='results'`` or an array ``fill``.
        layers : int or sequence of int, optional
            *(``fill=`` only)* Which layers' CELLS to draw, zero-based. With none,
            every layer. The unselected cells are masked out and the vertical extent
            is cropped to what remains -- leaving the axis at full height would put
            the two layers you asked for in a thin band with empty space above and
            below. A layer outside the model raises rather than drawing nothing.

            This is why ``layer=`` raises here: it means something else. ``layer=``
            overlays head PROFILES on the line-profile section; the filled section
            draws cells, so choosing them needs its own name.
        head_layers : int or sequence of int or None, default 0
            *(``fill='layer'`` only)* Whose water levels to draw over the geology.
            The default, layer 0, is the single water table every filled section drew
            before this parameter existed. A list draws one surface per layer, each
            labelled and coloured through :func:`~myflopy.viz.category_colors` -- so
            a given layer's water level keeps its colour across figures. ``None``
            draws none, for the geology alone.

            Ignored with ``fill='results'`` or an array fill: those paint the field
            onto the cells, so a line of the same quantity on top would say it twice.
        layer_labels : sequence of str, optional
            *(``fill=`` only)* Legend names for the layers -- your unit names rather
            than ``Layer 1..N``. One per layer, in model order.

            Defaults to the names the model's own build context carries
            (``ModelContext(surfaces=stack.build(vor))`` keeps them, so a model
            declared through the spec API already knows its layers are called
            "sand" and "clay"). A model built imperatively, or one whose ``surfaces``
            is a plain frame, has no names to find and falls back to ``Layer N``.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the section. ``'plotly'`` returns the interactive
            picture; ``'mpl'`` returns a static
            :class:`matplotlib.figure.Figure`. Accepts ``'interactive'`` and
            ``'matplotlib'``/``'static'`` as aliases; anything else raises.

            Applies to both branches: a model section renders through the same
            overlay builder the noun grammar uses
            (``model.hds.section(line=..., backend='mpl')``), and a grid section
            through ``GridSection``'s own cell-outline renderer. Until 2026-09-02
            this parameter did not exist on the verb and a caller who passed it got
            a Plotly figure back with no error at all -- it vanished into the
            ``**kwargs`` tail. Passing it is now the documented spelling; prefer it
            to ``.plot_mpl()`` on the result, which is the same renderer reached a
            second way.
        **kwargs
            Forwarded to the underlying section class.

        Returns
        -------
        XSection or GridSection or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default) a
            :class:`~myflopy.viz.Picture` -- ``XSection`` for a model, which also
            carries ``.ani``, or ``GridSection`` for a bare grid. Both answer
            ``.plot_mpl()``.

            With ``backend='mpl'``, a bare Matplotlib ``Figure``: use
            ``.savefig(path)`` and ``.axes[0]``, not the picture verbs.

        Raises
        ------
        ValueError
            If a model-only argument is given for a bare grid. It names the
            arguments, because the alternative is a ``TypeError`` from a constructor
            the caller never mentioned.

        See Also
        --------
        map : the same data in plan view.
        myflopy.plot.grid : the mesh itself.

        Examples
        --------
        >>> model.plot.section(cells=[1653, 651, 1241])
        >>> model.plot.section(line=line, layer=[0, 1], interpolate=True)
        >>> model.plot.section(cells=cells).save("section.png")
        >>> vor.plot.section(line)                       # geometry, no results
        >>> model.plot.section(line=line, backend="mpl") # a static mpl Figure
        >>> vor.plot.section(line, backend="mpl").savefig("grid_section.png")
        >>> model.plot.section(line=line, fill="layer", backend="mpl")
        >>> model.plot.section(line=line, fill="results", per=15, backend="mpl")
        >>> vor.plot.section(line, fill="layer", backend="mpl")   # geometry only
        >>> model.plot.section(line=line, fill="layer", layers=[0, 1], backend="mpl")
        >>> model.plot.section(line=line, fill="layer", head_layers=[0, 2], backend="mpl")
        >>> model.plot.section(line=line, fill="layer", backend="mpl",
        ...                    layer_labels=["fill", "sand", "clay", "till"])

        Bound form of :func:`myflopy.plot.section`.
        """

        return section(self.model, **_bound_args(locals(), "kwargs"), **kwargs)

    def surface(
        self,
        *,
        layer: int = 0,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        surf_type: str = "hds",
        resolution: int = 1000,
        use_rbf: bool = False,
        interpolator: str | None = None,
        clip=None,
        crs: str | None = None,
        xs=None,
        ys=None,
        zs=None,
        **kwargs,
    ) -> InterpolatedSurface:
        """A 3-D interpolated surface of this model. See :func:`myflopy.plot.surface`.

        Always Plotly, and always a height field. The 3-D grid VOLUME and 3-D
        particle pathlines are different shapes and live on :func:`grid` with
        ``backend="vtk"``: a ``backend=`` switch should change the renderer, not
        what is being drawn.

        Parameters
        ----------
        layer : int, default 0
            Zero-based layer to interpolate.
        per : int, optional
            Stress period (0-based). Mutually exclusive with ``kstpkper``.
        kstpkper : tuple of (int, int), optional
            Exact ``(timestep, period)``. Defaults to the model's FIRST output time.
        surf_type : {'hds', 'lyr'}, default 'hds'
            Interpolate the head field, or the layer elevation surface.
        resolution : int, default 1000
            Interpolation grid size per axis. Higher is smoother and slower.
        use_rbf : bool, default False
            Radial-basis interpolation instead of linear.
        interpolator : str, optional
            Override the interpolation method by name.
        clip : Path or geometry, optional
            Restrict the surface to cells intersecting this region. See also
            ``.clipped_fig()`` on the result.
        crs : str, optional
            Override the grid's CRS.
        xs, ys, zs : array-like, optional
            Supply the point cloud directly instead of reading it from a model.
        **kwargs
            Forwarded to :class:`~myflopy.modflow.mf6.grid.interpolated_surface
            .InterpolatedSurface`.

        Returns
        -------
        InterpolatedSurface
            A :class:`~myflopy.viz.Picture`; also carries ``.clipped_fig()`` and
            ``.surface_trace()`` for composing into a larger scene.

        See Also
        --------
        grid : the mesh itself, including the 3-D volume via ``backend="vtk"``.

        Examples
        --------
        >>> model.plot.surface(layer=0)
        >>> model.plot.surface(layer=0, surf_type="lyr")     # layer elevations
        >>> model.plot.surface(layer=0, resolution=400).html("surface.html")

        Bound form of :func:`myflopy.plot.surface`.
        """

        return surface(self.model, **_bound_args(locals(), "kwargs"), **kwargs)

    def grid(
        self,
        *,
        backend: str = "plotly",
        pathlines=None,
        vertical_exaggeration: float = 1.0,
        model_style: str = "wireframe",
        model_opacity: float = 0.25,
        pathline_cmap: str = "viridis",
        pathline_width: float = 4.0,
        show_edges: bool = True,
        off_screen: bool = True,
        **kwargs,
    ):
        """This model's mesh -- flat, or the 3-D volume with particle tracks.

        See :func:`myflopy.plot.grid`.

        The one picture :func:`map` cannot give you. A choropleth colours cells
        against a web basemap and so **requires a CRS**; this draws in the grid's own
        coordinates and needs none, which makes it the view for a grid you are still
        refining, before there is a model or a projection.

        Both backends draw the same subject -- this grid -- so ``backend=`` switches
        only the renderer. That is why the 3-D volume is ``grid`` and not
        ``surface``: ``surface`` means a height field.

        Parameters
        ----------
        backend : {'plotly', 'mpl', 'vtk'}, default 'plotly'
            ``'plotly'`` draws cell edges in 2-D -- fast, CRS-free, no basemap.
            ``'mpl'`` draws the same 2-D mesh through FloPy's own patch renderer and
            returns a static :class:`matplotlib.figure.Figure`; it is the same
            renderer as ``.plot_mpl()`` on the returned picture, offered here so the
            switch is spelled the same way on every verb. ``'vtk'`` renders the cell
            VOLUME in 3-D as an interactive PyVista scene, and needs the ``viz3d``
            extra (``pip install 'myflopy[viz3d]'``).

            All three draw the same subject -- this grid -- which is the rule the
            parameter follows everywhere: ``backend=`` changes the renderer, never
            what is being drawn.
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

        Returns
        -------
        GridMesh or VtkScene
            Both are :class:`~myflopy.viz.Picture`. ``GridMesh`` carries
            ``.plot_mpl()``, which is FloPy's own patch renderer. ``VtkScene``
            exposes ``.scene`` (the PyVista ``Plotter``) instead of ``.fig``.

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

        Examples
        --------
        >>> vor.plot.grid()                         # the 2-D mesh
        >>> vor.plot()                              # the same thing, shorthand
        >>> vor.plot.grid().plot_mpl()              # FloPy's matplotlib renderer
        >>> stack.plot.grid(["sand", "clay"])       # 3-D volume, a subset of layers
        >>> model.plot.grid(pathlines=run.track_records, backend="vtk")

        Bound form of :func:`myflopy.plot.grid`.
        """

        return grid(self.model, **_bound_args(locals(), "kwargs"), **kwargs)

    def animate(
        self,
        frames,
        *,
        backend: str = "plotly",
        title=None,
        dpi: int = 140,
        interval_ms: int = 700,
        **kwargs,
    ):
        """Flip through pictures. See :func:`myflopy.plot.animate`.

        Subject-free, like :meth:`mosaic` -- it animates the frames you hand it,
        which need not all come from this model. For this model's results over
        time, the grammar generates the frames for you:
        ``model.hds.animate(kind="map", over="period")``.

        A **combinator**, like :func:`mosaic` -- its first argument is the frames,
        not a subject. To animate one model's results over time, prefer the grammar
        (``model.hds.animate(kind="map", over="period")``), which generates the
        frames for you and calls this.

        Parameters
        ----------
        frames : sequence
            The frames, as bare pictures or ``(label, picture)`` pairs -- the same
            shapes :func:`mosaic` accepts. Unlabelled frames are numbered.
        backend : {'plotly', 'png'}, default 'plotly'
            ``'plotly'`` builds one live figure with play/pause and a slider: fast
            and interactive, but every frame must share a trace structure, and a
            choropleth re-embeds its geometry per frame, so the file grows with
            cells x frames. ``'png'`` rasterizes each frame and pages through them
            with a browser slider: frames need share NOTHING, so kinds can be mixed,
            and the size does not grow with cell count. Prefer ``'png'`` on a large
            grid.
        title : str, optional
            Figure title (plotly) or document title (png).
        dpi : int, default 140
            *(png only)* Raster resolution per frame.
        interval_ms : int, default 700
            *(png only)* Milliseconds per frame during playback.
        **kwargs
            Forwarded to the backend's animation class.

        Returns
        -------
        FrameAnimation or SliderAnimation
            Both are :class:`~myflopy.viz.Picture`. ``FrameAnimation.fig`` is the
            plotly figure; ``SliderAnimation`` has no single figure, so its ``.fig``
            raises and names ``.frames`` instead. ``SliderAnimation.export()``
            returns the richer ``StandaloneHtmlSlider`` handle when you want the
            frame manifest or the resume/progress machinery.

        Raises
        ------
        ValueError
            If ``frames`` is empty, if ``backend`` is neither value, or if the
            plotly backend is given frames that do not share a trace structure. That
            last one raises rather than silently falling back to raster -- swapping
            an interactive figure for a static page changes what you get.

        See Also
        --------
        myflopy.viz.mosaic : the same frames side by side instead of in sequence.
        myflopy.export_head_map_slider_html : a FloPy-rendered slider, which is a
            different picture rather than a second spelling of this one.

        Examples
        --------
        >>> frames = [model.plot.map(per=p, layer=0) for p in range(model.nper)]
        >>> plot.animate(frames).show()
        >>> plot.animate(frames, backend="png").html("heads.html")
        >>> plot.animate([("start", first), ("end", last)])
        >>> mixed = [("a map", model.plot.map()), ("a section", model.plot.section(cells=cells))]
        >>> plot.animate(mixed, backend="png")           # only png can mix kinds

        Bound form of :func:`myflopy.plot.animate`.
        """

        return animate(**_bound_args(locals(), "kwargs"), **kwargs)

    def mosaic(
        self,
        panels,
        *,
        ncols: int = 3,
        title: str | None = None,
        diff: bool = False,
        sync_views: bool = True,
        colorbar=None,
    ):
        """Compose any pictures into one figure. See :func:`myflopy.viz.mosaic`.

        Subject-free -- it takes the panels you hand it, which need not all come
        from this model. It lives here so the five verbs are discoverable in one
        place from a model you already have.

        The free-form composer of the unified view grammar: pass any mix of
        ``Choro`` choropleth maps and Plotly figures (``viz.Fig`` timeseries,
        cross-sections, ...) and get one figure back. The leaf
        ``.mosaic(by=...)`` verbs are sugar over this same composition.

        Parameters
        ----------
        panels
            A list of panels, or of ``(label, panel)`` pairs. A panel is either a
            ``Choro`` (anything exposing ``get_choropleth()``) or a Plotly figure
            whose traces are copied into its grid cell. Map panels contribute their
            cell trace **and** their overlays (contours, location markers,
            pathlines) via ``overlay_traces()``. Unlabeled panels take their figure
            title, else ``Panel <n>``.
        ncols
            Grid width; rows grow as needed.
        title
            Overall figure title.
        diff
            When ``True``, the shared map color scale is centered at zero
            (diverging), as used by the diff surfaces.
        sync_views
            Every map panel always *starts* framed to the same shared extent (the
            union of the panels' grid bounds) so the small multiples line up at the
            site. When ``sync_views`` is ``True`` (default), the panels are also
            wired to pan and zoom **together** live -- dragging or zooming one map
            moves the others (via a ``plotly_relayout`` handler injected at
            ``show()`` / ``write_html()`` time; inline notebook display shows the
            shared start view but not the live linking). Set ``False`` to let each
            map pan/zoom independently after the shared start. Panels without
            geometry (raw Plotly figures) are unaffected.
        colorbar
            Colorbar settings for the shared map color axis -- a dict of Plotly
            ``colorbar`` properties, or a **callable** ``(cmin, cmax) -> dict``.
            Panels are pooled onto one ``coloraxis``, which discards each panel's
            own ``colorbar``; without this a log-scaled mosaic silently reads in
            log10 units. The callable form exists because the useful labels depend
            on the pooled limits, and those are only known here -- e.g.
            ``colorbar=lambda lo, hi: {"tickvals": ..., "ticktext": ...}``.
            Ignored when no panel is a map.

        Returns
        -------
        Fig
            A single Plotly figure holding every panel as a subplot. Not a
            :class:`Picture` -- it is already the assembled figure, so use it
            directly (``.show()``, ``.write_html(...)``).

        See Also
        --------
        myflopy.plot.animate : the same panels in sequence rather than side by side.

        Bound form of :func:`myflopy.plot.mosaic`.

        Examples
        --------
        >>> viz.mosaic([
        ...     group.hds.map("F9b"),                       # a choropleth
        ...     group.packages.lak.results.stage.plot(),     # a timeseries
        ... ], ncols=2)

        Bound form of :func:`myflopy.plot.mosaic`.
        """

        return mosaic(
            panels,
            ncols=ncols,
            title=title,
            diff=diff,
            sync_views=sync_views,
            colorbar=colorbar,
        )


# --- one docstring per verb, shown wherever the verb appears -------------------
_SECTION_HEADS = (
    "Parameters", "Other Parameters", "Returns", "Yields", "Raises", "Warns",
    "See Also", "Notes", "References", "Examples",
)


def _split_sections(doc: str) -> tuple[str, str]:
    """Split a NumPy docstring into (summary + prose, the sectioned remainder).

    The split point is the first section heading -- a known title on its own
    line with a dashed underline beneath it.
    """

    lines = doc.splitlines()
    for i, line in enumerate(lines[:-1]):
        stripped = line.strip()
        underline = lines[i + 1].strip()
        if (
            stripped in _SECTION_HEADS
            and underline
            and set(underline) == {"-"}
            and len(underline) >= len(stripped)
        ):
            return "\n".join(lines[:i]).rstrip(), "\n".join(lines[i:])
    return doc.rstrip(), ""


#: Sections whose entries are parameter names, and so may be filtered by name.
#: Everything else in a NumPy docstring (See Also, Notes, Examples) is prose and
#: must be carried across untouched -- the previous filter walked the WHOLE
#: sectioned remainder, so dropping a parameter called `grid` or `section` also
#: deleted the See Also entry that happened to share its name.
_PARAMETER_SECTIONS = ("Parameters", "Other Parameters")


def _parameter_name(line: str) -> tuple[str, ...]:
    """The names a Parameters entry head declares, or ``()`` if it is not one.

    An entry head sits at column 0 and reads ``name : type``. NumPy allows a
    GROUPED head -- ``zmin, zmax : float, optional`` -- which documents two
    parameters in one entry, so this returns a tuple rather than a name.
    """

    if not line[:1].strip():
        return ()
    if ":" not in line:
        # A VAR-ARG head -- `**kwargs`, `*args` -- is a legal NumPy entry with no
        # type after it, and was invisible to the name-filter, so splicing a
        # docstring that already had one appended a second copy.
        stripped = line.strip()
        return (stripped,) if stripped.startswith("*") else ()
    head = line.split(":", 1)[0]
    return tuple(part.strip() for part in head.split(",") if part.strip())


def _bare_parameter_name(line: str) -> tuple[str, ...]:
    """An entry head with no type, valid only INSIDE a Parameters block.

    NumPy allows `panels` on its own line with the description indented under
    it, and `viz.mosaic` writes every one of its parameters that way. At column
    0 anywhere else the same line is prose or a section heading -- `Returns` is a
    perfectly good identifier -- so this is never applied outside the block.
    """

    stripped = line.strip()
    if not line[:1].strip() or ":" in stripped or stripped in _SECTION_HEADS:
        return ()
    parts = [part.strip() for part in stripped.split(",")]
    return tuple(parts) if parts and all(p.isidentifier() for p in parts) else ()


def _documented_names(sections: str) -> set[str]:
    """Every parameter name a docstring's Parameters blocks declare."""

    names: set[str] = set()
    in_params = False
    previous = ""
    for line in sections.splitlines():
        stripped = line.strip()
        if stripped and set(stripped) == {"-"} and previous.strip() in _SECTION_HEADS:
            in_params = previous.strip() in _PARAMETER_SECTIONS
        if in_params:
            names.update(_parameter_name(line) or _bare_parameter_name(line))
        previous = line
    return names


def _filter_parameters(sections: str, *, keep=None, drop=()) -> str:
    """Filter a docstring's Parameters entries by name, leaving prose alone.

    ``keep`` (when given) is the allowed set; ``drop`` is removed either way. A
    grouped head survives if any of its names does, and is rewritten to just
    those -- dropping ``zmin`` from ``zmin, zmax : float`` has to leave ``zmax``
    documented, where a name-match filter would either keep both or lose both.

    Only :data:`_PARAMETER_SECTIONS` are filtered. This is what makes the
    function safe to point at a whole docstring: a noun that does not accept
    ``values`` still wants the free verb's See Also and Examples.
    """

    kept: list[str] = []
    in_params = False
    skipping = False
    previous = ""
    for line in sections.splitlines():
        stripped = line.strip()
        if stripped and set(stripped) == {"-"} and previous.strip() in _SECTION_HEADS:
            in_params = previous.strip() in _PARAMETER_SECTIONS
            skipping = False
        names = (_parameter_name(line) or _bare_parameter_name(line)) if in_params else ()
        if names:
            wanted = [
                n for n in names
                if n not in drop and (keep is None or n in keep)
            ]
            skipping = not wanted
            if not skipping and len(wanted) != len(names):
                # A grouped head lost some of its names; re-head it with the rest.
                line = f"{', '.join(wanted)} :{line.split(':', 1)[1]}"
        elif skipping and in_params:
            # Continuation lines are indented; column 0 ends the entry.
            if line.startswith(" ") or not stripped:
                previous = line
                continue
            skipping = False
        if not skipping:
            kept.append(line)
        previous = line
    return "\n".join(kept)


def _drop_parameter(sections: str, name: str) -> str:
    """Remove one entry from a Parameters block (see :func:`_filter_parameters`)."""

    return _filter_parameters(sections, drop=(name,))


def inherit_map_docs(method, *, keep, extra: str = "") -> None:
    """Give a noun's ``map()`` the free verb's entries for the names it accepts.

    The noun's own summary and prose stay; the Parameters block is spliced from
    ``myflopy.plot.map`` and trimmed to ``keep``, so one edit to the verb's
    reference reaches every noun. ``extra`` is appended for parameters the verb
    has no notion of (``backend``, ``multiplier``, ``agg``, ...), which are
    genuine noun-local names -- a record noun reduces a table to one value per
    cell, work the verb never does because it is handed ``values`` directly.
    """

    own_doc = getattr(method, "_myflopy_own_doc", None)
    if own_doc is None:
        own_doc = method.__doc__ or ""
        method._myflopy_own_doc = own_doc
    own_prose, own_sections = _split_sections(inspect.cleandoc(own_doc).strip())
    _, sections = _split_sections(inspect.cleandoc(map.__doc__))
    sections = _filter_parameters(sections, keep=set(keep))
    if extra:
        sections = _merge_parameter_sections(
            "Parameters\n----------\n" + inspect.cleandoc(extra), sections
        )
    if own_sections:
        # The noun documents some of its own parameters -- `per`/`multiplier`
        # on `inputs.uzf`, `fill_value`/`agg` on a record field. Merging the
        # inherited block wholesale left the SAME name documented TWICE, the
        # generic copy last and therefore the one a reader ends on: measured on
        # `UzfInput.map`, the noun's "``'all'`` is not accepted here" entry was
        # followed by `plot.map`'s "Mutually exclusive with ``kstpkper``",
        # naming a parameter this noun deliberately does not accept. The
        # method's own entry wins and the inherited block contributes only the
        # names it does not already cover -- the same rule
        # `_inherit_verb_docs` already applies to the bound namespaces.
        own_names = {
            name
            for line in own_sections.splitlines()
            for name in _parameter_name(line)
        }
        sections = _filter_parameters(sections, drop=own_names)
        sections = _merge_parameter_sections(own_sections, sections)
    method.__doc__ = "\n\n".join(p for p in (own_prose, sections) if p)


def _merge_parameter_sections(own: str, inherited: str) -> str:
    """Fold an inherited Parameters block into one the method already has.

    Emits ONE `Parameters` heading: the method's own entries first, then the
    inherited ones, then whatever other sections each side carried. Two headings
    is malformed NumPy -- a reader stops at the first block and never sees the
    rest.
    """

    def _split_param_block(text: str) -> tuple[str, str]:
        lines, start, end = text.splitlines(), None, None
        for i, line in enumerate(lines[:-1]):
            if line.strip() in _PARAMETER_SECTIONS and set(lines[i + 1].strip()) == {"-"}:
                start = i + 2
                break
        if start is None:
            return "", text
        for j in range(start, len(lines) - 1):
            if lines[j].strip() in _SECTION_HEADS and set(lines[j + 1].strip()) == {"-"}:
                end = j
                break
        end = len(lines) if end is None else end
        body = "\n".join(lines[start:end]).rstrip()
        rest = "\n".join(lines[:start - 2] + lines[end:]).strip()
        return body, rest

    own_body, own_rest = _split_param_block(own)
    inh_body, inh_rest = _split_param_block(inherited)
    merged = "\n".join(part for part in (own_body, inh_body) if part.strip())
    parts = ["Parameters\n----------\n" + merged] if merged else []
    # The method's OWN Returns/See Also/Examples win outright; the inherited copy
    # of a section it already has is dropped rather than appended. Concatenating
    # both is the same malformation this function was written to end for
    # `Parameters`, one section along: measured on 7 of 8 spliced nouns, the LAST
    # Examples block a reader saw was the free verb's, so `help()` on a noun
    # ended with `>>> vor.plot.map(values=node_ids)` -- a grid map, and a call
    # the noun itself REFUSES.
    #
    # Emitted in canonical NumPy ORDER rather than own-then-inherited. Appending
    # blindly put an inherited `Returns` after the method's own `Examples`, which
    # is a docstring no renderer lays out the way its author meant.
    sections = dict(_split_named_sections(inh_rest))
    sections.update(dict(_split_named_sections(own_rest)))     # own wins
    parts += [sections[head] for head in _SECTION_HEADS if head in sections]
    return "\n\n".join(parts)


def _split_named_sections(text: str) -> list[tuple[str, str]]:
    """``[(head, block)]`` for a sectioned docstring remainder, in order.

    A block runs from its own heading to the next one, so it carries its
    underline and its body. Text before the first heading is not a section and is
    dropped -- callers pass a remainder that starts at one.
    """

    lines = text.splitlines()
    starts = [
        i for i, line in enumerate(lines[:-1])
        if line.strip() in _SECTION_HEADS
        and lines[i + 1].strip()
        and set(lines[i + 1].strip()) == {"-"}
    ]
    return [
        (lines[start].strip(),
         "\n".join(lines[start:starts[n + 1] if n + 1 < len(starts) else len(lines)]).rstrip())
        for n, start in enumerate(starts)
    ]


def _inherit_verb_docs(namespace) -> None:
    """Give each bound verb its free counterpart's reference, as ONE docstring.

    `model.plot.map(` is what an editor shows you, and a bound method whose
    docstring is "See :func:`myflopy.plot.map`" is a dead end at exactly the
    moment you wanted the parameter list. Duplicating the text onto thirteen
    bound methods would drift within a release, so the prose is written once on
    the free verb and re-headed here.

    Spliced rather than APPENDED. The previous version concatenated the two
    docstrings behind a row of dashes, which is malformed NumPy -- a dashed line
    is a section underline, and one with no title above it makes the whole
    docstring unparseable, so PyCharm drops structured rendering and shows the
    lot as plain text. Here the bound summary replaces the free summary and the
    sections carry over untouched, which is a single well-formed docstring.

    Runtime introspection (`help()`, Jupyter's `?`, most editor hovers) reads
    `__doc__`, so this reaches them. A purely static reader sees only the short
    source docstring, which is why that is written to stand alone -- and why the
    signature, which static readers DO see, is explicit rather than `**kwargs`.
    """

    for verb in ("map", "section", "surface", "grid", "mosaic", "animate"):
        method = getattr(namespace, verb, None)
        free = globals().get(verb)
        if method is None or free is None or not free.__doc__:
            continue
        # Splice from the method's ORIGINAL docstring, cached on first pass.
        # Without this the function is not idempotent -- it re-splices its own
        # output, and a second call took `ModelPlots.map` from 149 lines to 297,
        # with every section duplicated. Nothing calls it twice today; a module
        # reloaded in a notebook does.
        own_doc = getattr(method, "_myflopy_own_doc", None)
        if own_doc is None:
            own_doc = method.__doc__ or ""
            method._myflopy_own_doc = own_doc
        own_prose, own_sections = _split_sections(inspect.cleandoc(own_doc).strip())
        prose, sections = _split_sections(inspect.cleandoc(free.__doc__))
        if not sections:
            continue
        # Keep the free verb's extended prose, drop its summary line: the bound
        # method's own summary says what binding this subject means.
        extended = "\n".join(prose.splitlines()[1:]).strip()
        sections = _drop_parameter(sections, "source")
        if own_sections:
            # The bound method documents some of its own parameters (GridPlots.map
            # does). Appending the inherited block wholesale gave it TWO
            # `Parameters` headings, which is malformed NumPy -- so the method's
            # own entries win and the inherited block contributes only the names
            # it does not already cover.
            own_names = _documented_names(own_sections)
            sections = _filter_parameters(sections, drop=own_names)
            sections = _merge_parameter_sections(own_sections, sections)
        # Both guards make the splice IDEMPOTENT, which it was not: re-running it
        # over its own output grew `ModelPlots.section` from 169 lines to 179,
        # duplicating the extended prose, the `**kwargs` entry and the footer.
        # That matters now the generated text is written back into the source --
        # a module reloaded in a notebook would otherwise compound it.
        footer = f"Bound form of :func:`myflopy.plot.{verb}`."
        if extended and extended in own_prose:
            extended = ""
        parts = [own_prose, extended, sections]
        if footer not in own_prose and not sections.rstrip().endswith(footer):
            parts.append(footer)
        method.__doc__ = "\n\n".join(p for p in parts if p)


def _inherit_noun_docs() -> None:
    """Splice `plot.map`'s reference onto the NOUN verbs, trimmed to what each takes.

    Reaches DOWN from `myflopy.plot` (layer 7) into the explorer modules, exactly
    as `_inherit_verb_docs` already does for the bound namespaces -- the nouns
    cannot import upward to fetch it themselves.
    """

    tier1 = set(NOUN_MAP_PARAMS)
    layer_field = tier1 | set(LAYER_FIELD_MAP_PARAMS) | {"per", "layer"}
    record = tier1 | {"per", "layer"}

    #: Per-layer FIELD nouns: heads, concentration, temperature. They honour the
    #: Tier 2 names as well, because they genuinely have a time axis and a
    #: per-layer profile behind them.
    inherit_map_docs(DependentVariableFile.map, keep=layer_field)

    #: RECORD nouns: a package's period table reduced to one value per cell.
    #: Tier 1 only -- Tier 2 was measured byte-identical on these.
    #:
    #: Only nouns BELOW this module (layer 7) are listed. `GroupLakConnections`
    #: is layer 9, and reaching up to it here pulls in `myflopy.project`
    #: (layer 13) mid-import, which fails on a partially initialised
    #: `simulation.base`. It splices itself instead, at the bottom of its own
    #: module -- the mirror of what `StackPlots` needs, which is below this
    #: module and so cannot fetch the reference itself.
    for method in (
        CellBudgetResultsExplorer.map,
        CellPackageInputFieldExplorer.map,
        HfbPackageExplorer.map,
        HfbResultsExplorer.map,
        LakBudgetResultsExplorer.map,
        LakConnectionsExplorer.map,
        PRTPathlineView.map,
        SfrBudgetResultsExplorer.map,
        StageResultsExplorer.map,
        StaticArrayFieldExplorer.map,
        SurfaceWaterExchangeResultsExplorer.map,
        SurfaceWaterInputFieldExplorer.map,
        UzfFieldInputsExplorer.map,
    ):
        inherit_map_docs(method, keep=record)

    #: `model.inputs.uzf` is a record noun with NO layer axis, so it takes the
    #: `record` set MINUS `layer`. `finf()` reindexes to `ncpl` -- UZF is keyed
    #: by `ifno` and joined back through `packagedata.cellid[1]` -- which makes
    #: infiltration a per-COLUMN quantity, and `map()` has never named `layer`.
    #: Splicing the shared set documented a `layer` the signature does not
    #: accept, and one a reader would then pass: it lands in `**trace_kwargs`,
    #: reaches `plot.map(layer=...)`, and moves the HEAD tooltip while leaving
    #: every colour exactly where it was (measured: `custom_zs` identical,
    #: figure not).
    inherit_map_docs(UzfInput.map, keep=tier1 | {"per"})


_inherit_verb_docs(ModelPlots)
_inherit_verb_docs(GridPlots)

# `StackPlots` lives in `myflopy.layers` (L4) and cannot import this module
# (L5), so the binding happens here, where both are already in scope.
from myflopy.layers import StackPlots as _StackPlots  # noqa: E402

_inherit_verb_docs(_StackPlots)
_inherit_noun_docs()
