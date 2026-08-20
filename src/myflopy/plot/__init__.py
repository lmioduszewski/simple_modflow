"""The plotting front door: one verb per kind of picture (plans 8.3-8.4).

    from myflopy import plot

    plot.map(model, layer=0)                  # plan view
    plot.map(vor, values=drawdown)            # ... of any per-cell array
    plot.section(model, cells=[1653, 651])    # vertical slice through results
    plot.section(vor, line=line)              # ... of the grid itself
    plot.surface(model, layer=0)              # 3-D
    plot.grid(vor)                            # the bare mesh, no basemap
    plot.mosaic([a, b, c])                    # compose any pictures
    plot.animate(model)                       # frames through time

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
from myflopy.modflow.mf6.interactive_plotting import build_particle_tracking_scene
from myflopy.modflow.utils.animations import Animation
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.viz import Picture, mosaic

__all__ = [
    "map",
    "section",
    "surface",
    "grid",
    "mosaic",
    "animate",
    "ModelPlots",
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
    """A plan-view map: cell values, with optional contours, locations, hillshade.

    ``source`` is a model (draws its results) or a bare grid (draws the grid).
    ``values`` is any per-cell array -- heads, drawdown, K, a zone id, anything
    one-value-per-cell -- which is what makes this the single map verb rather
    than one function per quantity.

    Content is an option, never a separate call::

        plot.map(model, layer=0, contours=True)
        plot.map(model, locs=wells_gpkg, hillshade_path=hillshade)
        plot.map(vor, values=node_ids)        # replaces map_nodes()

    Returns a :class:`~myflopy.modflow.utils.datatypes.choros.Choro`, which is a
    Picture: it renders itself, and answers ``.fig``/``.show()``/``.save()``.
    """

    vor, is_model = _grid_of(source)
    if is_model:
        kwargs.setdefault("model", source)
    if values is not None:
        kwargs["custom_zs"] = list(values)
    return _choropleth_factory(vor, **kwargs)


def section(source, /, **kwargs):
    """A vertical slice.

    Through a MODEL, this is the results section -- a field against distance
    along the line, animatable over periods::

        plot.section(model, cells=[1653, 651, 1241])
        plot.section(model, line=line, layer=[0, 1])

    Through a bare GRID it is the geometry section: layers and cell edges, no
    results::

        plot.section(vor, line=line)

    Both are Pictures. They are different classes because they answer different
    questions, not because the API could not decide.
    """

    vor, is_model = _grid_of(source)
    if is_model:
        return XSection(model=source, **kwargs)
    return GridSection(vor=vor, **kwargs)


def surface(source, /, **kwargs) -> InterpolatedSurface:
    """A 3-D interpolated surface.

    ``plot.surface(model, layer=0)`` interpolates the model's heads;
    ``surf_type="lyr"`` gives the layer elevation instead. Replaces ``plot3d``,
    which drew the same kind of picture for the grid alone.

    Always Plotly, and always a height field ``z(x, y)``. The 3-D VOLUME and
    3-D pathlines are a different shape and live on :func:`grid` with
    ``backend="vtk"`` -- a ``backend=`` switch should change the renderer, not
    what is being drawn.
    """

    vor, is_model = _grid_of(source)
    if is_model:
        return InterpolatedSurface(model=source, **kwargs)
    return InterpolatedSurface(vor=vor, **kwargs)


def grid(source, /, *, backend: str = "plotly", pathlines=None, **kwargs):
    """The mesh itself -- flat in 2-D, or the layered volume in 3-D.

    ``backend="plotly"`` (the default) draws cell edges in the grid's own
    coordinates. It is the one picture ``map`` cannot give you: a choropleth
    colours cells against a web basemap and so requires a CRS, while this needs
    none, which makes it the view for a grid you are still refining. Replaces
    ``vor.plot2d()`` and FloPy's inherited ``VoronoiGrid.plot()``, whose
    matplotlib rendering is still right there as ``.plot_mpl()``.

    ``backend="vtk"`` draws the same subject in 3-D -- the cell VOLUME, as an
    interactive PyVista scene (needs the ``viz3d`` extra). ``pathlines=`` adds
    particle tracks as time-coloured tubes over it, mirroring
    ``map(pathlines=...)`` in plan view::

        model.plot.grid(pathlines=run.track_records, backend="vtk")
        stack.plot.grid(["sand", "clay"])          # vtk is the default there

    Both branches draw this grid, so ``backend=`` switches only the renderer --
    which is why the 3-D volume is ``grid`` and not ``surface``.
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


def animate(model, periods=None, **kwargs) -> Animation:
    """Frames through time for one model.

    Still the model-bound form. Plan 8.6 generalizes it to accept any sequence of
    pictures, so an animation can mix kinds the way :func:`mosaic` already does.
    """

    return Animation(model, periods=periods, **kwargs)


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

    def animate(self, periods=None, **kwargs) -> Animation:
        """Frames through this model's periods. See :func:`myflopy.plot.animate`."""

        return animate(self.model, periods=periods, **kwargs)

    def mosaic(self, panels, **kwargs):
        """Compose any pictures into one figure. See :func:`myflopy.viz.mosaic`.

        Subject-free -- it takes the panels you hand it, which need not all come
        from this model. It lives here so the five verbs are discoverable in one
        place from a model you already have.
        """

        return mosaic(panels, **kwargs)
