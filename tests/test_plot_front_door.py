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


@pytest.mark.canonical
@pytest.mark.slow
def test_the_exact_call_shapes_the_notebooks_use(canonical_run):
    """The notebooks are not executed by the suite, so a migration that reads
    fine and dies on the first run is the failure mode. These are the shapes
    8.3 rewrote `build_choropleth(...)` into, verbatim."""

    import matplotlib

    matplotlib.use("Agg")
    model, vor = canonical_run, canonical_run.vor
    heads = list(model.hds.array(layer=0))

    # canonical_02: `wt_map = plot.map(model, values=list(wt), layer=0)`
    assert len(plot.map(model, values=heads, layer=0).fig.data) >= 1

    # canonical_02 / canonical_fast_tour: grid form, then the mpl backend
    assert plot.map(vor, values=heads, layer=0).plot_mpl() is not None

    # canonical_06 / bearcreek: grid form with extra Choro options
    assert plot.map(vor, values=heads, layer=0, zmin=0, zmax=100) is not None


# --- 8.4: the same verbs, bound to the object ---------------------------------
BOUND_VERBS = ("map", "section", "surface", "animate", "mosaic")


def test_the_model_namespace_delegates_to_the_module_verbs():
    """`model.plot.map()` must BE `plot.map(model)`, not a second code path.

    Checked by call rather than by identity: the methods are bound wrappers, so
    `is` cannot hold. What must hold is that each one routes to this module's
    function with the model as the subject.
    """

    from types import SimpleNamespace

    from myflopy.plot import ModelPlots

    for verb in BOUND_VERBS:
        assert callable(getattr(ModelPlots, verb)), verb

    seen = {}
    model = SimpleNamespace(vor=SimpleNamespace(ncpl=2))
    ns = ModelPlots(model)
    for verb in ("map", "section", "surface", "animate"):
        original = getattr(plot, verb)
        try:
            setattr(plot, verb, lambda source, _v=verb, **kw: seen.setdefault(_v, source))
            getattr(ns, verb)()
        finally:
            setattr(plot, verb, original)
        assert seen[verb] is model, f"{verb} did not pass the bound model through"


def test_the_namespace_does_not_re_enter_grid_dispatch():
    """A model whose grid is unresolved carries `vor = None`, which `_grid_of`
    would read as "this IS a grid". The namespace holds its subject, so it must
    pass the model straight through rather than re-deriving what it is."""

    from types import SimpleNamespace

    from myflopy.plot import ModelPlots

    ungridded = SimpleNamespace(vor=None)
    assert plot._grid_of(ungridded) == (ungridded, False)  # the trap

    captured = []
    original = plot.map
    try:
        plot.map = lambda source, **kw: captured.append(source)
        ModelPlots(ungridded).map()
    finally:
        plot.map = original
    assert captured == [ungridded]


@pytest.mark.parametrize("gone", ["cor", "srf", "section"])
def test_the_replaced_model_accessors_are_gone(gone):
    """8.4 deletes rather than deprecates. `section` is the subtle one -- it was
    a real method on the model, and `model.plot.section()` is now the only
    spelling."""

    from myflopy.modflow.mf6.simulation.base import SimulationBase

    assert not hasattr(SimulationBase, gone)


def test_model_plot_shadows_flopys_delegate_but_leaves_it_reachable():
    """`SimulationBase.plot` used to be a one-line delegate to FloPy's
    `MFSimulation.plot`. The namespace takes the name; the renderer stays
    reachable as `model.sim.plot(...)`, which is where it actually lives."""

    from myflopy.modflow.mf6.simulation.base import SimulationBase

    assert isinstance(SimulationBase.plot, property)


@pytest.mark.canonical
@pytest.mark.slow
def test_a_model_map_keeps_the_layer_elevation_hover(canonical_run):
    """The regression this stage fixes.

    `model.cor()` defaulted `show_layer_elevs=True`; `_choropleth_factory`
    hardcoded `False`, so from 8.3 until now `plot.map(model)` silently dropped
    the 'Top of Model' and 'Layer N Bottom' hover rows. The default now resolves
    from the grid, which is the condition 25 call sites used to spell out.
    """

    picture = plot.map(canonical_run, layer=0)
    assert picture.show_layer_elevs is True
    picture.fig  # assembling the figure is what builds the hover
    hover = picture._hover_dict
    assert "Top of Model" in hover, sorted(hover)
    assert any(k.endswith("Bottom") for k in hover), sorted(hover)


@pytest.mark.canonical
@pytest.mark.slow
def test_a_bare_grid_map_still_opts_out_of_layer_elevations(canonical_run):
    """The other side of the same default: a grid carrying no `gdf_topbtm`
    cannot build that hover, and Choro's own `True` default raises on it. That
    asymmetry is the whole reason the factory resolves rather than passes through.
    """

    from myflopy.modflow.mf6.grid.plotting import _choropleth_factory

    vor = canonical_run.vor
    saved = getattr(vor, "gdf_topbtm", None)
    try:
        vor.gdf_topbtm = None
        assert _choropleth_factory(vor).show_layer_elevs is False
    finally:
        vor.gdf_topbtm = saved
