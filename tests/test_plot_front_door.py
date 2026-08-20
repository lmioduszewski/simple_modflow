"""`myflopy.plot` -- the verbs, and what they deliberately do not include.

Built in 8.3 as module-level functions; 8.4 bound the same functions onto the
objects (`model.plot.*`, `vor.plot.*`) and deleted the accessors they replace.

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

VERBS = ("map", "section", "surface", "grid", "mosaic", "animate")


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
    `plot3d` -> `surface`, `map_nodes` -> `map(values=...)`, `plot2d` -> `grid`
    (8.4 -- it is the mesh, which `map` cannot draw without a CRS),
    `contours` -> `map(contours=True)`. `timeseries` is absent because a chart
    belongs to a node, which knows the model's periods."""

    assert not hasattr(plot, retired)


def test_the_verbs_dispatch_on_what_you_pass():
    """A thing with results draws its results; anything else draws geometry.

    Duck-typed so loaded runs, live builds and group members all work without
    this module importing three model classes. The discriminator is RESULTS, not
    `.vor` -- see the next test for why.
    """

    from types import SimpleNamespace

    grid = SimpleNamespace(ncpl=2)
    model = SimpleNamespace(vor=grid, hds=object())

    assert plot._grid_of(model) == (grid, True)
    assert plot._grid_of(grid) == (grid, False)


def test_a_layer_stack_is_not_mistaken_for_a_model():
    """A `LayerStack`/`LayerBuildResult` carries `.vor` too.

    Dispatching on `.vor` alone read one as a model and sent it down the results
    path, where it died on `.hds`. A layer stack is geometry: it has a grid and
    no results, so it must answer the GRID verbs.
    """

    from types import SimpleNamespace

    grid = SimpleNamespace(ncpl=5)
    stack = SimpleNamespace(vor=grid, top=[], botm=[], names=["sand"])

    assert plot._grid_of(stack) == (grid, False)


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
BOUND_VERBS = ("map", "section", "surface", "grid", "animate", "mosaic")


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
    for verb in ("map", "section", "surface", "grid", "animate"):
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


# --- 8.4b: the grid scope -----------------------------------------------------
RETIRED_GRID_ALIASES = (
    "mapit", "choropleth", "dash_selector", "show", "map_nodes",
    "show_selected_cells", "show_overlapping_geometry",
    "plot3d", "plot2d", "plottri", "cross_section",
)


@pytest.mark.parametrize("alias", RETIRED_GRID_ALIASES)
def test_the_eleven_grid_aliases_are_gone(alias):
    """Deleted, not deprecated. Five were real pictures and live on as
    `vor.plot.map/section/grid` or options on them; the rest were dead
    (`mapit` needs folium, which is not a dependency), duplicated, or not
    pictures at all (`dash_selector` blocked on a Dash server just to be READ).
    """

    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    assert not hasattr(VoronoiGridPlus, alias)


def test_the_grid_namespace_lives_below_myflopy_plot():
    """Not a style point -- a ratchet constraint.

    `voronoi.py` sits below `myflopy.plot` in the import graph, so binding the
    namespace from there would point upward and force a deferred import, which
    `test_import_layering` pins exactly. Building `GridPlots` in
    `grid/plotting.py` keeps the edge downward and the count unchanged.
    """

    import json
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[1]
    layers = json.loads((root / "tests/import_layers.json").read_text())
    assert layers["myflopy.modflow.mf6.grid.plotting"] <= layers["myflopy.modflow.mf6.grid.voronoi"]


def test_calling_the_grid_namespace_draws_the_grid():
    """`vor.plot()` must keep working: it used to be FloPy's `VoronoiGrid.plot`,
    and notebooks call it that way. `__call__` routes it to the mesh."""

    from myflopy.plot import GridMesh, GridPlots

    class FakeVor:
        ncpl = 2
        x_coords_by_node = [[0, 1, 1, 0], [1, 2, 2, 1]]
        y_coords_by_node = [[0, 0, 1, 1], [0, 0, 1, 1]]

    ns = GridPlots(FakeVor())
    assert isinstance(ns(), GridMesh)
    assert len(ns().fig.data) == 2


def test_the_mesh_is_a_picture_and_needs_no_crs():
    """The reason `grid` is a verb rather than a `map` option: Choro requires a
    CRS for its basemap, and a grid under construction may not have one. The
    mesh draws in the grid's own coordinates."""

    from myflopy.plot import GridMesh

    class NoCrsVor:
        crs = None
        x_coords_by_node = [[0, 1, 1, 0]]
        y_coords_by_node = [[0, 0, 1, 1]]

    mesh = GridMesh(NoCrsVor())
    assert isinstance(mesh, viz.Picture)
    assert isinstance(mesh.fig, viz.Fig)
    assert len(mesh.fig.data) == 1


def test_flopys_renderer_survives_the_shadowing():
    """`vor.plot` is now the namespace, so FloPy's method is unreachable by that
    name. It is not gone -- `GridMesh.plot_mpl` calls it unbound, which is what
    makes the shadowing a rename rather than a removal."""

    import inspect

    from flopy.utils.voronoi import VoronoiGrid

    from myflopy.plot import GridMesh

    assert "VoronoiGrid.plot" in inspect.getsource(GridMesh.plot_mpl)
    assert callable(VoronoiGrid.plot)


@pytest.mark.canonical
@pytest.mark.slow
def test_the_grid_verbs_draw_on_the_canonical_grid(canonical_run):
    """Imports prove nothing; these must actually render."""

    import matplotlib

    matplotlib.use("Agg")
    vor = canonical_run.vor

    assert len(vor.plot.grid().fig.data) == int(vor.ncpl)
    assert vor.plot.grid().plot_mpl() is not None
    assert len(vor.plot.map().fig.data) >= 1

    # `select=` absorbs show_selected_cells / show_overlapping_geometry.
    picked = vor.plot.map(select=[0, 1, 2])
    assert tuple(picked.fig.data[0].selectedpoints) == (0, 1, 2)


def test_the_cell_selector_is_still_reachable_from_a_grid():
    """`vor.dash_selector` was deleted in 8.4b, but the workflow it served was
    not: box/lasso-select cells on the map, read back the cell-id list, query
    that batch. That has always lived on `Choro`; the grid attribute was a
    property wrapping a default map.

    Pinned because the reachable path is now two steps
    (`vor.plot.map().dash_selector()`) and nothing else in the suite exercises
    it -- it launches a Dash server, so it cannot be CALLED here.
    """

    from myflopy.modflow.utils.datatypes.choros import Choro
    from myflopy.plot import GridPlots

    assert callable(Choro.dash_selector)
    # A method, not a property: reading it must not start a server. This is the
    # bug the old grid alias had.
    assert not isinstance(Choro.dash_selector, property)
    assert "dash_selector" not in dir(GridPlots)


# --- 8.5b: the 3-D scene -------------------------------------------------------
def test_grid_takes_a_backend_not_a_different_verb():
    """Why the 3-D volume is `grid` and not `surface`.

    `surface` means a height field `z(x, y)`. `vtk_3d` draws a cell VOLUME and
    the particle scene draws POLYLINE tubes, so hanging them off `surface` with
    a `backend=` switch would silently change the SUBJECT, not the renderer --
    the exact thing this module's docstring forbids. `grid` already means "the
    mesh itself", so a 3-D layered mesh is that same subject, redrawn.
    """

    import inspect

    assert "backend" in inspect.signature(plot.grid).parameters
    assert "pathlines" in inspect.signature(plot.grid).parameters
    assert "backend" not in inspect.signature(plot.surface).parameters


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"backend": "opengl"}, "'plotly' or 'vtk'"),
        ({"backend": "plotly", "pathlines": object()}, "map\\(pathlines"),
        ({"backend": "vtk"}, "pathlines="),
    ],
)
def test_grid_says_what_went_wrong(kwargs, match):
    """Three ways to ask for a picture that does not exist, each answered with
    the spelling that does."""

    import pytest as _pytest
    from types import SimpleNamespace

    with _pytest.raises(ValueError, match=match):
        plot.grid(SimpleNamespace(ncpl=2), **kwargs)


def test_a_vtk_scene_is_a_picture_that_is_not_plotly():
    """The third renderer. A PyVista scene answers the same verbs; `.fig` raises
    and names `.scene`, because `fig` is documented as the Plotly figure."""

    pv = pytest.importorskip("pyvista")

    from myflopy.viz import VtkScene

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(pv.Line((0, 0, 0), (1, 1, 1)))
    scene = VtkScene(plotter)
    try:
        assert isinstance(scene, viz.Picture)
        with pytest.raises(TypeError, match=r"\.scene"):
            scene.fig
    finally:
        plotter.close()


def test_the_3d_dependencies_are_named_by_the_optional_helper():
    """`pyvista`/`trame` are the `viz3d` extra. Before 8.5b `_optional.require`
    had never heard of them and one call site imported pyvista with no guard at
    all, so a missing dependency surfaced as a bare ImportError."""

    from myflopy._optional import _EXTRA_FOR_MODULE

    assert _EXTRA_FOR_MODULE["pyvista"] == "viz3d"
    assert _EXTRA_FOR_MODULE["trame"] == "viz3d"


@pytest.mark.parametrize("gone", ["ParticleTrackingScene", "export_particle_tracking_html"])
def test_the_bespoke_particle_exporters_are_gone(gone):
    """Both folded into `VtkScene`: the dataclass became the Picture, and the
    exporter became `.html(path)`. `ModelVisualization`'s two wrappers around
    them had zero callers and went too."""

    from myflopy.modflow.mf6 import interactive_plotting

    assert not hasattr(interactive_plotting, gone)
    assert not hasattr(interactive_plotting.ModelVisualization, "particle_tracking_scene")
    assert not hasattr(interactive_plotting.ModelVisualization, "particle_tracking_html")
