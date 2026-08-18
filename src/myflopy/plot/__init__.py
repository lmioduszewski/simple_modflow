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
from myflopy.modflow.mf6.grid.plotting import GridSection, build_choropleth
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
    return build_choropleth(grid, **kwargs)


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
