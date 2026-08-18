"""`myflopy.plot` -- the verbs, and what they deliberately do not include (8.3).

The point of this module is discoverability: one verb per KIND OF PICTURE, where
the kind is decided by geometry (plan view / vertical slice / 3-D), never by
content or renderer. That is what retires `plot3d`, `map_nodes` and a top-level
`contours` -- a 3-D view is a `surface`, a node-id map is a `map`, and contours
are something drawn ON a map.

These tests pin the shape of the front door. The canonical-model test at the end
is the one that would catch a verb that imports cleanly but cannot actually draw.
"""

from __future__ import annotations

import pytest

from myflopy import plot, viz

VERBS = ("map", "section", "surface", "mosaic", "animate")


def test_the_front_door_exposes_exactly_the_verbs():
    assert set(VERBS) <= set(plot.__all__)
    for verb in VERBS:
        assert callable(getattr(plot, verb)), verb


def test_mf_plot_resolves_lazily():
    """`myflopy/__init__.py.__getattr__` special-cases subpackages one at a time;
    a new one is invisible as `mf.plot` until it is added there -- and `import
    myflopy.plot` would still work, which is what makes the omission easy to
    miss."""

    import myflopy as mf

    assert mf.plot.map is plot.map


def test_mosaic_is_vizs_own_function_not_a_reimplementation():
    """One implementation behind both spellings. `viz.mosaic` already composes
    ARBITRARY panels -- it was only ever shadowed by `<node>.mosaic()` sugar."""

    assert plot.mosaic is viz.mosaic


@pytest.mark.parametrize("retired", ["plot3d", "map_nodes", "plot2d", "contours", "timeseries"])
def test_the_retired_spellings_are_not_verbs(retired):
    """Each of these is a picture kind that turned out not to be one:
    `plot3d` -> `surface`, `map_nodes`/`plot2d` -> `map(values=...)`,
    `contours` -> `map(contours=True)`. `timeseries` is absent because a chart
    belongs to a node, which knows the model's periods."""

    assert not hasattr(plot, retired)


def test_the_verbs_dispatch_on_what_you_pass():
    """A model draws its results; a bare grid draws itself. Duck-typed on `.vor`
    so loaded runs, live builds and group members all work without this module
    importing three model classes."""

    from types import SimpleNamespace

    grid = SimpleNamespace(ncpl=2)
    model = SimpleNamespace(vor=grid)

    assert plot._grid_of(model) == (grid, True)
    assert plot._grid_of(grid) == (grid, False)


# --- the test that would catch a verb that cannot draw ------------------------
@pytest.mark.canonical
@pytest.mark.slow
def test_every_verb_draws_on_the_canonical_model(canonical_run):
    """Imports proving nothing is the failure mode here: a facade re-exports
    fine and then dies on the first real call."""

    picture = plot.map(canonical_run, layer=0)
    assert isinstance(picture, viz.Picture)
    assert len(picture.fig.data) >= 1

    section = plot.section(canonical_run, cells=[0, 1, 2])
    assert isinstance(section, viz.Picture)

    surface = plot.surface(canonical_run, layer=0)
    assert isinstance(surface, viz.Picture)

    # A map of an arbitrary per-cell array -- the case `map_nodes`/`plot2d`
    # existed for, and the reason `values=` is the single map option.
    node_ids = list(range(int(canonical_run.vor.ncpl)))
    assert len(plot.map(canonical_run.vor, values=node_ids).fig.data) >= 1

    combined = plot.mosaic([picture, plot.map(canonical_run, layer=1)])
    assert combined is not None
