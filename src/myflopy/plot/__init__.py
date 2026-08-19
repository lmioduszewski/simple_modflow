"""The plotting front door: one verb per kind of picture (plan 8.3).

    from myflopy import plot

    plot.map(model, layer=0)                  # plan view
    plot.map(vor, values=drawdown)            # ... of any per-cell array
    plot.section(model, cells=[1653, 651])    # vertical slice through results
    plot.section(vor, line=line)              # ... of the grid itself
    plot.surface(model, layer=0)              # 3-D
    plot.mosaic([a, b, c])                    # compose any pictures
    plot.animate(model)                       # frames through time

**Geometry chooses the verb.** Not content, and not renderer. A map is a plan
view whatever is drawn on it, so contours, observation locations and a hillshade
are *options* on ``map`` rather than verbs of their own -- which is why there is
no ``plot.contours``. The same rule retires ``plot3d`` (a 3-D view is
``surface``) and ``map_nodes`` (a map whose values are node ids).

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

The same verbs exist on the objects themselves (``model.plot.map()``,
``vor.plot.map()``) -- plan 8.4. These are the same functions, so there is one
implementation behind both spellings.
"""

from __future__ import annotations

from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.plotting import GridSection, _choropleth_factory
from myflopy.modflow.utils.animations import Animation
from myflopy.modflow.utils.datatypes.choros import Choro
from myflopy.modflow.utils.datatypes.xsections import XSection
from myflopy.viz import Picture, mosaic

__all__ = [
    "map",
    "section",
    "surface",
    "mosaic",
    "animate",
    "ModelPlots",
    "Picture",
    "Choro",
    "XSection",
    "GridSection",
    "InterpolatedSurface",
]


def _grid_of(source):
    """The Voronoi grid for ``source``, and whether ``source`` is a model.

    A model carries its grid as ``.vor``; a grid is its own. Duck-typed rather
    than isinstance-checked so loaded runs, live builds and group members all
    work without this module importing three model classes.
    """

    grid = getattr(source, "vor", None)
    if grid is not None:
        return grid, True
    return source, False


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

    grid, is_model = _grid_of(source)
    if is_model:
        kwargs.setdefault("model", source)
    if values is not None:
        kwargs["custom_zs"] = list(values)
    return _choropleth_factory(grid, **kwargs)


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

    grid, is_model = _grid_of(source)
    if is_model:
        return XSection(model=source, **kwargs)
    return GridSection(vor=grid, **kwargs)


def surface(source, /, **kwargs) -> InterpolatedSurface:
    """A 3-D interpolated surface.

    ``plot.surface(model, layer=0)`` interpolates the model's heads;
    ``surf_type="lyr"`` gives the layer elevation instead. Replaces ``plot3d``,
    which drew the same kind of picture for the grid alone.

    A ``backend="vtk"`` option lands in plan 8.5, which folds the PyVista scenes
    (``vtk_3d``, particle-tracking) in behind this same verb.
    """

    grid, is_model = _grid_of(source)
    if is_model:
        return InterpolatedSurface(model=source, **kwargs)
    return InterpolatedSurface(vor=grid, **kwargs)


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

        return f"ModelPlots({getattr(self.model, 'name', '?')!r}: map, section, surface, animate, mosaic)"

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
