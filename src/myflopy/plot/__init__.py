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

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.plotting import (
    GridMesh,
    GridPlots,
    GridSection,
    _choropleth_factory,
)
from myflopy.modflow.mf6.interactive_plotting import (
    SliderAnimation,
    build_particle_tracking_scene,
)
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.viz import FrameAnimation, Picture, mosaic

logger = get_logger(__name__)

__all__ = [
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
    **trace_kwargs,
) -> Choro:
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
        Interpolation grid size used to build the contours. Higher is smoother
        and slower.
    contour_method : {'linear', 'cubic', 'nearest'}, default 'linear'
        Interpolation method for that grid.
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
    **trace_kwargs
        Anything else rides through to the ``go.Choroplethmap`` trace --
        ``zmid``, ``colorbar``, ``reversescale``, ``showscale``. These are
        genuinely open-ended and Plotly owns their names, so they are validated
        LATE, at render time, not here.

    Returns
    -------
    Choro
        A :class:`~myflopy.viz.Picture`: renders itself in Jupyter, and answers
        ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. Also carries
        ``.plot_mpl()`` for a static Matplotlib rendering and ``.ani`` for the
        animation over periods.

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
    return _choropleth_factory(
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
    **kwargs
        Forwarded to the underlying section class.

    Returns
    -------
    XSection or GridSection
        A :class:`~myflopy.viz.Picture`. ``XSection`` (model) also carries
        ``.ani``; both answer ``.plot_mpl()``.

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
    """

    vor, is_model = _grid_of(source)
    if is_model:
        return XSection(
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
    return GridSection(vor=vor, line=line, **kwargs)


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
    backend : {'plotly', 'vtk'}, default 'plotly'
        ``'plotly'`` draws cell edges in 2-D -- fast, CRS-free, no basemap.
        ``'vtk'`` renders the cell VOLUME in 3-D as an interactive PyVista
        scene, and needs the ``viz3d`` extra
        (``pip install 'myflopy[viz3d]'``).
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
        If ``backend`` is neither ``'plotly'`` nor ``'vtk'``, if ``pathlines``
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
        raise ValueError(f"backend must be 'plotly' or 'vtk', not {backend!r}.")

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
        **trace_kwargs,
    ) -> Choro:
        """This model's plan-view map. See :func:`myflopy.plot.map`."""

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
        **kwargs,
    ) -> XSection:
        """A vertical slice through this model. See :func:`myflopy.plot.section`."""

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
        """A 3-D interpolated surface of this model. See :func:`myflopy.plot.surface`."""

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


def _drop_parameter(sections: str, name: str) -> str:
    """Remove one entry from a Parameters block -- its head line and its body.

    A bound verb already supplies ``source``; leaving it documented tells the
    reader to pass an argument the signature does not have.
    """

    lines = sections.splitlines()
    kept, skipping = [], False
    for line in lines:
        head = line.split(":")[0].strip()
        starts_entry = line[:1].strip() != "" and not line.startswith(" ")
        if starts_entry and (head == name or head.startswith(f"{name} ")):
            skipping = True
            continue
        if skipping:
            # Continuation lines of the dropped entry are indented; anything at
            # column 0 begins the next entry (or the next section) and ends it.
            if line.startswith(" ") or not line.strip():
                continue
            skipping = False
        kept.append(line)
    return "\n".join(kept)


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
        own = inspect.cleandoc(method.__doc__ or "").strip()
        prose, sections = _split_sections(inspect.cleandoc(free.__doc__))
        if not sections:
            continue
        # Keep the free verb's extended prose, drop its summary line: the bound
        # method's own summary says what binding this subject means.
        extended = "\n".join(prose.splitlines()[1:]).strip()
        sections = _drop_parameter(sections, "source")
        parts = [own, extended, sections, f"Bound form of :func:`myflopy.plot.{verb}`."]
        method.__doc__ = "\n\n".join(p for p in parts if p)


_inherit_verb_docs(ModelPlots)
_inherit_verb_docs(GridPlots)

# `StackPlots` lives in `myflopy.layers` (L4) and cannot import this module
# (L5), so the binding happens here, where both are already in scope.
from myflopy.layers import StackPlots as _StackPlots  # noqa: E402

_inherit_verb_docs(_StackPlots)
