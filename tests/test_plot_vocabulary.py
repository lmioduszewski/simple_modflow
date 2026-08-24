"""The plotting vocabulary, pinned across every scope (plan 8.7).

Phase 8 built one set of verbs and bound it in four places: the module
(`myflopy.plot`), the model, the grid, and the layer stack. The failure this
file exists to catch is DRIFT -- a verb added at one scope and forgotten at
another, or a scope quietly inventing a spelling of its own.

The plan asked for "all three scopes expose the same verb set". They do not, and
should not: a bare grid has no results so it cannot answer `surface`, and a layer
stack has no time so it cannot answer `animate`. What must hold is narrower and
more useful -- every scope draws from ONE vocabulary, and each scope's subset is
declared here on purpose.
"""

from __future__ import annotations

import inspect

import pytest

from myflopy import plot, viz

#: Every verb in the vocabulary. Nothing outside this set may appear on a scope.
VOCABULARY = {"map", "section", "surface", "grid", "mosaic", "animate"}

#: What each scope answers, and why it is not the whole vocabulary.
#: `mosaic`/`animate` are COMBINATORS -- subject-free -- so they live at module
#: level and on the model, which is the common entry point. See the ledger.
SCOPES = {
    "module": VOCABULARY,
    "model": VOCABULARY,
    # A bare grid has no results: no field to interpolate (`surface`) and no
    # periods to flip through (`animate`).
    "grid": {"map", "section", "grid"},
    # A layer stack is geometry with no time, so no `animate`; it DOES have
    # surfaces (layer contacts) and a mesh.
    "stack": {"map", "section", "surface", "grid"},
}


def _verbs(obj) -> set[str]:
    """Public callables on a namespace class, or lowercase names in `__all__`."""

    if inspect.ismodule(obj):
        return {name for name in obj.__all__ if name[0].islower()}
    return {
        name
        for name, _ in inspect.getmembers(obj, inspect.isfunction)
        if not name.startswith("_")
    }


def _namespace(scope):
    from myflopy.layers import StackPlots
    from myflopy.plot import GridPlots, ModelPlots

    return {"module": plot, "model": ModelPlots, "grid": GridPlots, "stack": StackPlots}[scope]


@pytest.mark.parametrize("scope", sorted(SCOPES))
def test_each_scope_answers_exactly_its_declared_verbs(scope):
    """Adding a verb to one scope and forgetting another fails here, naming both."""

    actual = _verbs(_namespace(scope))
    expected = SCOPES[scope]
    assert actual == expected, (
        f"{scope} verbs drifted: missing={sorted(expected - actual)}, "
        f"unexpected={sorted(actual - expected)}"
    )


@pytest.mark.parametrize("scope", sorted(SCOPES))
def test_no_scope_invents_a_verb_outside_the_vocabulary(scope):
    """The stronger half: a scope may answer FEWER verbs, never different ones.

    This is what stops `vor.plot.draw()` or `stack.plot.render()` appearing --
    a second spelling for a picture that already has one.
    """

    assert _verbs(_namespace(scope)) <= VOCABULARY


def test_the_module_is_the_union_of_every_scope():
    """No verb may exist only on an object. If a scope can do it, the free
    function can too -- that is what makes `plot.map(model)` and
    `model.plot.map()` interchangeable rather than two APIs."""

    union = set().union(*SCOPES.values())
    assert _verbs(plot) == union == VOCABULARY


@pytest.mark.parametrize("verb", sorted(VOCABULARY))
def test_every_verb_is_documented_where_it_is_defined(verb):
    """A verb with no docstring is a verb nobody can discover from the REPL."""

    function = getattr(plot, verb)
    assert function.__doc__, f"plot.{verb} has no docstring"
    assert len(function.__doc__.strip()) > 40, f"plot.{verb}'s docstring is a stub"


# --- the contract every verb's OUTPUT answers ---------------------------------
def test_every_picture_class_answers_the_output_verbs():
    """Layer 3: whatever a verb returns, you can show it and write it out.

    Checked against the real subclasses rather than a stub -- the mistake 8.1
    made and ledger 135 records.
    """

    for cls in viz.Picture.__subclasses__():
        if not cls.__module__.startswith("myflopy."):
            continue                                    # test-local stubs
        for method in ("show", "save", "html", "_repr_mimebundle_"):
            assert callable(getattr(cls, method, None)), f"{cls.__name__} lacks {method}"


@pytest.mark.parametrize(
    ("cls_name", "names_instead"),
    [("MplPicture", "axes"), ("VtkScene", "scene"), ("SliderAnimation", "frames")],
)
def test_a_non_plotly_picture_says_what_to_use_instead(cls_name, names_instead):
    """Three renderers now. The two that are not Plotly must not pretend to have
    a `fig` -- and must name the attribute that replaces it, or the error is a
    dead end."""

    import myflopy.modflow.mf6.interactive_plotting as interactive
    from myflopy.viz import MplPicture, VtkScene

    cls = {"MplPicture": MplPicture, "VtkScene": VtkScene,
           "SliderAnimation": interactive.SliderAnimation}[cls_name]
    source = inspect.getsource(cls.fig.fget if isinstance(cls.fig, property) else cls.fig)
    assert "raise" in source, f"{cls_name}.fig should raise"
    assert names_instead in source, f"{cls_name}.fig must name `.{names_instead}`"


# --- what the vocabulary deliberately excludes --------------------------------
@pytest.mark.parametrize(
    "retired",
    ["cor", "xs", "srf", "plot2d", "plot3d", "map_nodes", "plottri", "mapit",
     "choropleth", "vtk_3d", "preview", "views", "thickness_map", "surface_3d"],
)
def test_the_retired_spellings_are_gone_everywhere(retired):
    """One name per picture. Every one of these drew something the vocabulary
    now draws, under a name you had to memorize separately."""

    from myflopy.layers import LayerBuildResult, LayerStack
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
    from myflopy.modflow.mf6.simulation.base import SimulationBase

    for owner in (SimulationBase, VoronoiGridPlus, LayerStack, LayerBuildResult):
        assert not hasattr(owner, retired), f"{owner.__name__}.{retired} survived"


def test_model_visualize_is_gone_but_its_exporters_are_not():
    """8.6b deleted the namespace, not the pictures. The `export_*_slider_html`
    functions render through FloPy's PlotMapView, which `animate(backend="png")`
    does not reproduce -- so they are a distinct picture, not a duplicate."""

    import myflopy as mf
    from myflopy.modflow.mf6 import interactive_plotting
    from myflopy.modflow.mf6.simulation.base import SimulationBase

    assert not hasattr(SimulationBase, "visualize")
    assert not hasattr(interactive_plotting, "ModelVisualization")
    for exporter in (
        "export_matplotlib_slider_html",
        "export_head_map_slider_html",
        "export_head_layer_mosaic_slider_html",
        "export_cross_section_slider_html",
    ):
        assert callable(getattr(mf, exporter)), f"mf.{exporter} went missing"


# --- the allowlist that names tests by string ---------------------------------
def test_the_slow_test_allowlist_names_only_tests_that_exist():
    """`conftest._SLOW_TESTS` marks heavy tests by NAME, so a rename silently
    un-marks one and moves real work into the fast lane.

    Plan 8.7 found two entries orphaned by 8.5b's renames -- both pyvista VTK
    exports. Nothing failed, which is the problem: the tests kept passing, just
    in the wrong lane. This asserts every name still resolves.
    """

    import ast
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parent

    # Read the two sets by PARSING conftest rather than importing it: importing
    # would re-run its module body, and a test that has side effects on the
    # fixture module it is checking is a bad trade for three lines saved.
    def _string_set(name: str) -> set[str]:
        tree = ast.parse((root / "conftest.py").read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == name for target in node.targets
            ):
                return {
                    element.value
                    for element in node.value.elts
                    if isinstance(element, ast.Constant)
                }
        raise AssertionError(f"conftest.py no longer defines {name}")

    slow_tests = _string_set("_SLOW_TESTS")
    retired = _string_set("_RETIRED_SLOW_TESTS")
    defined = set()
    for path in root.glob("test_*.py"):
        defined.update(re.findall(r"^def (test_\w+)", path.read_text(encoding="utf-8"), re.M))

    stale = slow_tests - defined
    assert not stale, (
        f"_SLOW_TESTS names tests that no longer exist: {sorted(stale)} -- "
        "they were renamed or deleted, and are no longer marked slow"
    )

    # And the register of retired entries must stay retired.
    revived = retired & defined
    assert not revived, (
        f"these names are in _RETIRED_SLOW_TESTS but exist again: {sorted(revived)}"
    )


# --- docstrings are the API surface an editor shows ---------------------------
@pytest.mark.parametrize("verb", sorted(VOCABULARY))
def test_every_verb_documents_its_parameters(verb):
    """A verb that takes `**kwargs` and does not list them is undiscoverable.

    These forward to picture classes with 15-39 constructor parameters, so the
    docstring is the only place the caller can learn what is accepted -- there
    is no signature to read.
    """

    import inspect

    doc = inspect.getdoc(getattr(plot, verb)) or ""
    for section in ("Parameters", "Returns", "Examples"):
        assert section in doc, f"plot.{verb} has no {section} section"
    assert ">>>" in doc, f"plot.{verb} has no runnable example"


@pytest.mark.parametrize("scope", ["model", "grid", "stack"])
def test_the_bound_verbs_carry_the_full_reference(scope):
    """`model.plot.map(` is what an editor shows on hover.

    The bound methods are thin forwarders, so their own docstrings are short by
    design. `_inherit_verb_docs` appends the free function's full reference to
    each, from the single source -- duplicating it onto thirteen methods would
    drift within a release.
    """

    import inspect

    namespace = _namespace(scope)
    for verb in SCOPES[scope]:
        doc = inspect.getdoc(getattr(namespace, verb)) or ""
        assert "Parameters" in doc, f"{scope}.plot.{verb} lost its parameter reference"
        assert f"myflopy.plot.{verb}" in doc, f"{scope}.plot.{verb} does not name its source"


# --- the hover a map gets by default ------------------------------------------
@pytest.mark.canonical
@pytest.mark.slow
def test_a_model_map_gets_the_sectioned_hover(canonical_run):
    """The verb and the grammar must agree on what a map's hover looks like.

    `model.hds.map()` built a `HoverSpec` and `model.plot.map()` did not, so the
    same picture reached by the two documented routes carried different hovers --
    the sectioned one, or a flat `Cell No. / Area / x / y` dump. The default now
    resolves from the map's `type`, like `show_layer_elevs` does.
    """

    from myflopy.modflow.utils.datatypes.hover import HoverSpec

    model = canonical_run
    direct = model.plot.map(layer=0)
    grammar = model.hds.map(layer=0)
    for picture in (direct, grammar):
        picture.fig
        assert isinstance(picture._resolved_hover_spec(), HoverSpec)

    assert direct.get_choropleth().hovertemplate == grammar.get_choropleth().hovertemplate


@pytest.mark.canonical
@pytest.mark.slow
def test_a_bare_grid_map_keeps_the_flat_hover(canonical_run):
    """The other side of that default: the sectioned hover needs model context --
    layers, periods, dates -- which a bare grid does not have. `None` is the
    right answer there, not a broken section."""

    vor = canonical_run.vor
    picture = vor.plot.map(values=list(canonical_run.hds.array(layer=0)))
    picture.fig
    assert picture._resolved_hover_spec() is None
    assert "Cell No." in (picture.get_choropleth().hovertemplate or "")


@pytest.mark.canonical
@pytest.mark.slow
def test_every_model_backed_map_gets_the_sectioned_hover(canonical_run):
    """Breadth, not just the two entry points ledger 142 compared.

    A map reached through the verb, the dependent-variable reader, a
    static-array package or a cell-stress package must all hover the same way.
    They arrive by different code paths -- some pass `hover_spec`, some
    `custom_hover`, some neither -- so the only honest check is to build one of
    each and look.
    """

    from myflopy.modflow.utils.datatypes.hover import HoverSpec

    model = canonical_run
    surfaces = {
        "plot.map":        lambda: model.plot.map(layer=0),
        "hds.map":         lambda: model.hds.map(layer=0),
        "npf.k.map":       lambda: model.packages.npf.k.map(),
        "sto.ss.map":      lambda: model.packages.sto.ss.map(),
        "rch.inputs.map":  lambda: model.packages.rch.inputs.map(per=0),
        "ghb.results.q":   lambda: model.packages.ghb.results.q.map(per=0),
        "drn.inputs.elev": lambda: model.packages.drn.inputs.elev.map(per=0),
        "lak.results.q":   lambda: model.packages.lak.results.q.map(per=0),
        "sfr.results.q":   lambda: model.packages.sfr.results.q.map(per=0),
        "uzf.inputs.finf": lambda: model.packages.uzf.inputs.finf.map(per=0),
    }
    flat = []
    for name, build in surfaces.items():
        picture = build()
        picture.fig
        if not isinstance(picture._resolved_hover_spec(), HoverSpec):
            flat.append(name)
    assert not flat, f"these maps fell back to the flat hover: {flat}"


# --- signatures, not just docstrings (plan 8.8) -------------------------------
#
# PyCharm and Pylance are STATIC: they read the `def` line and never execute the
# module, so `__signature__`, `functools.wraps` and a beautifully written
# docstring all reach `help()` and none of them reach the editor. The only thing
# that does is an explicit parameter list. These tests are the price of that --
# they keep the explicit copies honest so the duplication cannot rot.

#: Each verb's forwarding chain, nearest link first. A parameter's default is
#: owned by the FIRST link that names it; that is the value the verb must mirror.
def _chains():
    from myflopy.modflow.mf6.grid import plotting as gp
    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
    from myflopy.modflow.mf6.grid.plotting import Choro
    from myflopy.modflow.mf6.interactive_plotting import (
        SliderAnimation,
        build_particle_tracking_scene,
    )
    from myflopy.modflow.utils.datatypes.xsections import XSection

    return {
        "map": [gp._choropleth_factory, Choro.__init__],
        "section": [XSection.__init__],
        "surface": [InterpolatedSurface.__init__],
        "grid": [build_particle_tracking_scene],
        "animate": [SliderAnimation.__init__],
    }


def _named(func) -> dict:
    """{name: default} for a callable's named parameters, minus var-args."""

    return {
        name: p.default
        for name, p in inspect.signature(func).parameters.items()
        if p.kind not in (p.VAR_KEYWORD, p.VAR_POSITIONAL) and name != "self"
    }


_SECTION_HEADS = {
    "Parameters", "Other Parameters", "Returns", "Yields", "Raises", "Warns",
    "See Also", "Notes", "References", "Examples",
}


def _documented(func) -> set[str]:
    """Parameter names carrying an entry in a NumPy Parameters block.

    Written out rather than regexed because NumPy style groups names on one
    line -- `zmin, zmax : float, optional` documents two parameters, and a
    naive `"\\nzmax "` search calls the second one undocumented.
    """

    lines = (inspect.getdoc(func) or "").splitlines()
    names, inside = set(), False
    for i, line in enumerate(lines):
        head = line.strip()
        underlined = (
            i + 1 < len(lines)
            and lines[i + 1].strip()
            and set(lines[i + 1].strip()) == {"-"}
        )
        if head in _SECTION_HEADS and underlined:
            inside = head in ("Parameters", "Other Parameters")
            continue
        if inside and line and not line[0].isspace():
            for part in line.split(":")[0].split(","):
                part = part.strip().lstrip("*")
                if part.isidentifier():
                    names.add(part)
    return names


@pytest.mark.parametrize("verb", sorted(_chains()))
def test_a_named_default_matches_the_link_that_owns_it(verb):
    """The verb restates its chain's defaults; this proves it restates them right.

    A mirrored default that drifts is worse than no signature at all: the editor
    would confidently show `zoom=13` while the call actually produced something
    else. Nothing here is hand-written -- the expected values are read off the
    forwarding target at runtime.
    """

    free = _named(getattr(plot, verb))
    owner, expected = {}, {}
    for link in _chains()[verb]:
        for name, default in _named(link).items():
            # First link wins: `_choropleth_factory` resolves `show_layer_elevs`
            # from the grid and passes the RESULT to `Choro`, so the factory's
            # `None` is the default a caller sees, not `Choro`'s `True`.
            if name not in expected and default is not inspect._empty:
                owner[name], expected[name] = link.__qualname__, default
    wrong = [
        f"{verb}.{name}: signature says {free[name]!r}, "
        f"{owner[name]} uses {expected[name]!r}"
        for name in sorted(set(free) & set(expected))
        if free[name] != expected[name]
    ]
    assert not wrong, "mirrored defaults drifted:\n  " + "\n  ".join(wrong)


@pytest.mark.parametrize("verb", sorted(VOCABULARY))
def test_the_bound_model_verb_mirrors_the_free_one(verb):
    """`model.plot.map` must accept exactly what `plot.map` does, minus `source`.

    Two hand-written parameter lists is the cost of a facade an editor can read.
    This is what keeps them one list in practice: add a parameter to the free
    verb, forget the bound one, and the diff shows up here by name.
    """

    from myflopy.plot import ModelPlots

    free = _named(getattr(plot, verb) if verb != "mosaic" else viz.mosaic)
    free.pop("source", None)
    bound = _named(getattr(ModelPlots, verb))
    assert bound == free, (
        f"model.plot.{verb} and plot.{verb} disagree.\n"
        f"  only on the free verb:  {sorted(set(free) - set(bound))}\n"
        f"  only on the bound verb: {sorted(set(bound) - set(free))}\n"
        f"  differing defaults:     "
        f"{ {k: (free[k], bound[k]) for k in set(free) & set(bound) if free[k] != bound[k]} }"
    )


@pytest.mark.parametrize("scope", ["grid", "stack"])
def test_a_narrower_scope_stays_a_subset_of_the_free_verb(scope):
    """`vor.plot.map` may offer FEWER parameters than `plot.map`, never other ones.

    A bare grid has no periods, so dropping `per=` is right. Inventing a name
    that the shared factory does not know would not be -- it would typo-check
    clean and then vanish into `**kwargs`.
    """

    #: Genuinely local to the scope: not forwarded to the free verb at all.
    local = {
        "grid": {"select"},
        "stack": {"basemap", "layers", "backend", "color_by", "scale", "cmap",
                  "width", "height", "x", "y", "line", "resolution", "colorscale",
                  "opacity", "show_grid", "legend", "title"},
    }[scope]
    namespace = _namespace(scope)
    strays = {}
    for verb in SCOPES[scope]:
        free = set(_named(getattr(plot, verb))) - {"source"}
        bound = set(_named(getattr(namespace, verb)))
        extra = bound - free - local
        if extra:
            strays[verb] = sorted(extra)
    assert not strays, f"{scope} scope invented parameters the free verb lacks: {strays}"


@pytest.mark.parametrize("verb", sorted(VOCABULARY))
def test_every_signature_parameter_is_documented(verb):
    """A parameter an editor completes but the docstring never mentions is a trap.

    The reverse direction (documented but not accepted) is the one that bit
    hardest: `plot.grid` documented `layers`/`scale`/`color_by`/`cmap`, which
    live on the LAYER-stack builder and were never reachable through it.
    """

    func = getattr(plot, verb) if verb != "mosaic" else viz.mosaic
    # `**trace_kwargs` is documented, and should be -- but it is an open tail,
    # not a parameter, so it belongs in neither set.
    var_args = {
        n for n, p in inspect.signature(func).parameters.items()
        if p.kind in (p.VAR_KEYWORD, p.VAR_POSITIONAL)
    }
    documented = _documented(func) - var_args
    accepted = set(_named(func)) - {"source"}
    assert not accepted - documented, (
        f"plot.{verb} accepts but never documents: {sorted(accepted - documented)}"
    )
    assert not documented - accepted - {"source"}, (
        f"plot.{verb} documents parameters it cannot accept: "
        f"{sorted(documented - accepted - {'source'})}"
    )


def test_the_verbs_a_user_calls_are_not_bare_kwargs_forwarders():
    """The whole point of 8.8, pinned: no `(*args, **kwargs)` on a picture verb.

    `FieldMappable`'s namespace-level `field=` sugar is exempt and stays exempt:
    it dispatches to accessors with genuinely incompatible signatures (a lake's
    `connections` map has no `per=`; `DrnInput.map` takes `per` positionally), so
    one merged signature could only be achieved by lying. Its docstring says so
    and sends you to the leaf. Everything else must be explicit.
    """

    from myflopy.layers import StackPlots
    from myflopy.plot import GridPlots, ModelPlots

    def is_bare(func) -> bool:
        """Var-args and nothing else. Taking NO arguments is not bare -- it is
        the most informative signature there is, and `vor.plot.grid()` earns it:
        a mesh is fully determined by its grid."""

        params = [p for n, p in inspect.signature(func).parameters.items() if n != "self"]
        varargs = [p for p in params if p.kind in (p.VAR_KEYWORD, p.VAR_POSITIONAL)]
        return bool(varargs) and len(varargs) == len(params)

    bare = [
        f"{label}.{verb}"
        for label, namespace in (("model.plot", ModelPlots), ("vor.plot", GridPlots),
                                 ("stack.plot", StackPlots))
        for verb in _verbs(namespace)
        if is_bare(getattr(namespace, verb))
    ]
    bare += [
        f"plot.{verb}" for verb in VOCABULARY
        if is_bare(getattr(plot, verb) if verb != "mosaic" else viz.mosaic)
    ]
    assert not bare, f"these verbs still tell an editor nothing: {bare}"
