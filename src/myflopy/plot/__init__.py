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
        return source, False
    has_results = any(hasattr(source, attr) for attr in ("hds", "conc", "temp"))
    return grid, has_results


def map(source, /, values=None, **kwargs) -> Choro:      # noqa: A001 - the verb IS `map`
    """Draw a plan-view map of one value per grid cell.

    The single map verb. Contours, observation markers, a hillshade and particle
    pathlines are all **options** here rather than verbs of their own, because a
    map is a plan view whatever is drawn on it.

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
    contour_width : float, default 1.5
    contour_name : str, optional
        Legend name for the contour trace.
    contour_clip : bool, default True
        Clip contours to the active domain instead of the full grid extent.
    contour_resolution : int, default 150
        Interpolation grid size used to build the contours. Higher is smoother
        and slower.
    contour_method : {'linear', 'cubic', 'nearest'}, default 'linear'
        Interpolation method for that grid.
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
    show_layer_elevs : bool, optional
        Add model-top and per-layer-bottom rows to the hover. Defaults to
        whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
        because forcing it on a grid without them raises.
    show_mounding : bool, default False
        Add head-above-initial (mounding) to the hover.
    hover_heads, hover_ks : bool
        Include heads / hydraulic conductivity in the legacy flat hover.
    hover : HoverSpec, optional
        Replace the sectioned hover outright. See
        :mod:`myflopy.modflow.utils.datatypes.hover`.
    hover_layers : {'active', 'active+strip', 'all', 'none'}, optional
        How the per-layer profile renders in the hover.
    hover_surfaces : bool, optional
        Add the model-top / layer-bottom table to the sectioned hover.
    hover_fields : sequence of str, optional
        Extra columns to append to the hover.
    custom_hover : dict, optional
        Legacy flat hover: ``{label: per-cell sequence}``. Supplying it
        suppresses the default sectioned hover.
    **kwargs
        Anything else rides through to the ``go.Choroplethmap`` trace --
        ``zmid``, ``colorbar``, ``reversescale``, ``showscale``. These are
        validated LATE, by Plotly at render time, not here.

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
        kwargs.setdefault("model", source)
    if values is not None:
        kwargs["custom_zs"] = list(values)
    return _choropleth_factory(vor, **kwargs)


def section(source, /, **kwargs):
    """Draw a vertical slice through a model's results, or a grid's geometry.

    Through a MODEL this is the results section -- the field against distance
    along the line. Through a bare GRID it is the geometry section: layers and
    cell edges, no results. They are different classes because they answer
    different questions.

    Parameters
    ----------
    source : SimulationBase or VoronoiGridPlus
        Positional only. A model gives results; a grid gives geometry.
    cells : int or list of int, optional
        Cell indices defining the section path, in order. Mutually exclusive
        with ``line``.
    line : LineString or Path, optional
        The section line as a geometry or a vector file. For a bare grid this
        is the only way to specify the path, and it is positional.
    per : int, optional
        Stress period (0-based). Mutually exclusive with ``kstpkper``.
    kstpkper : tuple of (int, int), optional
        Exact ``(timestep, period)``.
    layer : int or list of int, default 0
        Layer(s) to draw. A list overlays several.
    x_or_y : {'x', 'y'}, default 'x'
        Which coordinate becomes the horizontal axis.
    spacing : int, default 10
        Sample spacing along the line, in model units.
    num_points : int, default 100
        Number of samples when interpolating.
    interpolate : bool, default False
        Interpolate between cell centers rather than stepping cell to cell.
    use_rbf : bool, default True
        Use radial-basis interpolation when ``interpolate`` is True.
    interpolator : str, optional
        Override the interpolation method by name.
    extrapolate_beyond_section_ends : bool, default False
        Extend the section past the first and last cell centers.
    show_model_top : bool, default True
        Draw the model-top profile.
    show_model_btm : bool, default False
        Draw layer-bottom profiles.
    surf_type : {'hds', 'lyr'}, default 'hds'
        Section the head field, or the layer elevations.
    section_name : str, optional
        Legend name for the traces.
    clip : Path or geometry, optional
        Restrict the section to cells intersecting this region.
    animation_kstpkpers : sequence of tuple, optional
        The periods ``.ani`` steps through; defaults to every output time.

    Returns
    -------
    XSection or GridSection
        A :class:`~myflopy.viz.Picture`. ``XSection`` (model) also carries
        ``.ani``; both answer ``.plot_mpl()``.

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
        return XSection(model=source, **kwargs)
    return GridSection(vor=vor, **kwargs)


def surface(source, /, **kwargs) -> InterpolatedSurface:
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
    if is_model:
        return InterpolatedSurface(model=source, **kwargs)
    return InterpolatedSurface(vor=vor, **kwargs)


def grid(source, /, *, backend: str = "plotly", pathlines=None, **kwargs):
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

    Other Parameters
    ----------------
    layers : str or int or list, optional
        (``backend="vtk"`` on a layer stack.) Which layers to show: a name, an
        index, or a list mixing them. Colours stay keyed to each layer's
        position, so a subset looks the same as it does in the full stack.
    scale : float, default 8
        (``backend="vtk"``.) Vertical exaggeration.
    color_by : str, default 'layer'
        (``backend="vtk"``.) Cell scalar to colour by.
    cmap : str, default 'tab10'
        (``backend="vtk"``.) Colormap for that scalar.

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
        is given with the plotly backend, or if the VTK backend is asked for
        without either pathlines or a layer stack.

    See Also
    --------
    map : values over the cells, on a basemap.
    surface : a 3-D height field, which is a different shape.

    Examples
    --------
    >>> vor.plot.grid()                         # the 2-D mesh
    >>> vor.plot()                              # the same thing, shorthand
    >>> vor.plot.grid().plot_mpl()              # FloPy's matplotlib renderer
    >>> stack.plot.grid(["sand", "clay"])       # 3-D volume, a subset of layers
    >>> model.plot.grid(pathlines=run.track_records, backend="vtk")
    """

    if backend == "plotly":
        if pathlines is not None:
            raise ValueError(
                "pathlines are only drawn by the 3-D scene; pass backend='vtk' "
                "for tubes over the grid, or use map(pathlines=...) in plan view."
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
    return build_particle_tracking_scene(source, pathlines, **kwargs)


def animate(frames, *, backend: str = "plotly", title=None, **kwargs):
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
        (``backend="png"``.) Raster resolution per frame.
    interval_ms : int, default 700
        (``backend="png"``.) Milliseconds per frame during playback.

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
        return FrameAnimation(frames, title=title, **kwargs)
    if backend == "png":
        return SliderAnimation(frames, title=title, **kwargs)
    raise ValueError(f"backend must be 'plotly' or 'png', not {backend!r}.")


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
    """

    def __init__(self, model):
        """Bind the plotting verbs to ``model``."""

        self.model = model

    def __repr__(self):
        """Name the subject and the verbs, so tab-completion has a companion."""

        return f"ModelPlots({getattr(self.model, 'name', '?')!r}: map, section, surface, grid, animate, mosaic)"

    # The bare names below resolve to this module's functions, not to these
    # methods -- class scope is not in a method's name-lookup chain.
    def map(self, values=None, **kwargs) -> Choro:      # noqa: A003 - the verb IS `map`
        """This model's plan-view map. See :func:`myflopy.plot.map`."""

        return map(self.model, values=values, **kwargs)

    def section(self, **kwargs) -> XSection:
        """A vertical slice through this model. See :func:`myflopy.plot.section`."""

        return section(self.model, **kwargs)

    def surface(self, **kwargs) -> InterpolatedSurface:
        """A 3-D interpolated surface of this model. See :func:`myflopy.plot.surface`."""

        return surface(self.model, **kwargs)

    def grid(self, *, backend: str = "plotly", pathlines=None, **kwargs):
        """This model's mesh -- flat, or the 3-D volume with particle tracks.

        See :func:`myflopy.plot.grid`.
        """

        return grid(self.model, backend=backend, pathlines=pathlines, **kwargs)

    def animate(self, frames, **kwargs):
        """Flip through pictures. See :func:`myflopy.plot.animate`.

        Subject-free, like :meth:`mosaic` -- it animates the frames you hand it,
        which need not all come from this model. For this model's results over
        time, the grammar generates the frames for you:
        ``model.hds.animate(kind="map", over="period")``.
        """

        return animate(frames, **kwargs)

    def mosaic(self, panels, **kwargs):
        """Compose any pictures into one figure. See :func:`myflopy.viz.mosaic`.

        Subject-free -- it takes the panels you hand it, which need not all come
        from this model. It lives here so the five verbs are discoverable in one
        place from a model you already have.
        """

        return mosaic(panels, **kwargs)


# --- one docstring per verb, shown wherever the verb appears -------------------
def _inherit_verb_docs(namespace) -> None:
    """Append each free verb's full docstring to its bound counterpart.

    `model.plot.map(` is what an editor shows you, and a bound method whose
    docstring is "See :func:`myflopy.plot.map`" is a dead end at exactly the
    moment you wanted the parameter list. Duplicating the text onto thirteen
    bound methods would drift within a release, so the bound docstring keeps its
    own short note about what binding means and the full reference is appended
    from the single source.

    Runtime introspection -- `help()`, Jupyter's `?`, and most editor hovers --
    reads `__doc__`, so this reaches them. A purely static reader still sees the
    short source docstring, which is why that is written to stand alone.
    """

    for verb in ("map", "section", "surface", "grid", "mosaic", "animate"):
        method = getattr(namespace, verb, None)
        free = globals().get(verb)
        if method is None or free is None or not free.__doc__:
            continue
        own = (method.__doc__ or "").strip()
        method.__doc__ = f"{own}\n\n{'-' * 70}\nFull reference (`myflopy.plot.{verb}`):\n\n{free.__doc__}"


_inherit_verb_docs(ModelPlots)
_inherit_verb_docs(GridPlots)

# `StackPlots` lives in `myflopy.layers` (L4) and cannot import this module
# (L5), so the binding happens here, where both are already in scope.
from myflopy.layers import StackPlots as _StackPlots  # noqa: E402

_inherit_verb_docs(_StackPlots)
