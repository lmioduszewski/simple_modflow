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

import functools
import inspect
import re

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
        # `select`/`select_style`/`select_color` were grid-only until ledger 160
        # promoted them to every scope; the grid entry is deliberately empty now.
        "grid": set(),
        "stack": {"basemap", "layers", "backend", "color_by", "scale", "cmap",
                  "width", "height", "x", "y", "line", "resolution", "colorscale",
                  "opacity", "show_grid", "legend", "title",
                  # vtk sheet detail on `surface(backend="vtk")`; the same name
                  # the free `grid` verb uses for its 3-D scene.
                  "show_edges"},
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


# --- the renderer switch, at every scope that has one -------------------------
#: Verbs whose `backend=` chooses between Plotly and Matplotlib, and the scopes
#: that offer them. `surface` is deliberately absent: it draws a 3-D height
#: field, `InterpolatedSurface` has no `plot_mpl`, and matplotlib has no
#: equivalent picture -- offering the parameter there would be a name that
#: cannot act. `animate`'s `backend` is `{'plotly', 'png'}`, a rasteriser over
#: frames rather than a second renderer of one subject, so it is not this switch.
_MPL_BACKEND_VERBS = ("map", "section", "grid")


def _mpl_capable(scope: str, verb: str) -> bool:
    """Whether this scope draws `verb` through a Matplotlib renderer."""

    if verb not in _MPL_BACKEND_VERBS or verb not in SCOPES[scope]:
        return False
    # The layer stack's pictures are Matplotlib NATIVELY (`MplPicture`), so there
    # is no Plotly renderer to switch away from and no `backend=` to offer.
    return scope != "stack"


@pytest.mark.parametrize("scope", sorted(SCOPES))
@pytest.mark.parametrize("verb", _MPL_BACKEND_VERBS)
def test_the_renderer_switch_is_spelled_the_same_at_every_scope(scope, verb):
    """`backend=` must be a NAMED parameter wherever the renderer can switch.

    The bug this pins was silent, which is why a type check is not enough:
    `plot.section` had no `backend` parameter and a `**kwargs` tail, so
    `model.plot.section(line=line, backend="mpl")` stored the argument on the
    section object and returned a Plotly picture with no error at all. The noun
    tier had offered the switch since 8.2; the verb tier never had it.
    """

    if not _mpl_capable(scope, verb):
        pytest.skip(f"{scope}.{verb} has no Plotly/Matplotlib pair")
    function = getattr(_namespace(scope), verb, None) or getattr(plot, verb)
    assert "backend" in inspect.signature(function).parameters, (
        f"{scope}.{verb} can render both ways but does not name `backend`"
    )


@pytest.mark.parametrize("verb", _MPL_BACKEND_VERBS)
def test_an_unknown_backend_raises_rather_than_being_ignored(verb, canonical_run):
    """A misspelled renderer must fail loudly.

    The whole defect was a swallowed argument, so the guard that matters is not
    "does 'mpl' work" but "does anything else refuse to be ignored".
    """

    subject = canonical_run if verb != "grid" else canonical_run.vor
    kwargs = {"cells": [0, 1, 2]} if verb == "section" else {}
    with pytest.raises(ValueError, match="backend must be"):
        getattr(plot, verb)(subject, backend="not-a-renderer", **kwargs)


@pytest.mark.parametrize(
    ("scope", "verb"),
    [("module", "map"), ("module", "section"), ("model", "map"), ("model", "section"),
     ("grid", "map"), ("grid", "section"), ("grid", "grid")],
)
def test_the_mpl_backend_returns_one_bare_figure_everywhere(scope, verb, canonical_run):
    """One return type across the family, so a loop over pictures cannot trip.

    The three Matplotlib renderers behind this switch disagree about what they
    hand back -- `Choro.plot_mpl` a `Figure`, the `figs` cross-section helper a
    `(fig, ax)` tuple, FloPy's patch renderer an `Axes` -- and normalising at
    each call site was how the family ended up with three answers.
    """

    import matplotlib

    matplotlib.use("Agg")
    import shapely as shp
    from matplotlib.figure import Figure

    xmin, ymin, xmax, ymax = canonical_run.vor.gdf_vorPolys.total_bounds
    line = shp.LineString([(xmin + 1, (ymin + ymax) / 2), (xmax - 1, (ymin + ymax) / 2)])

    if scope == "module":
        subject = canonical_run if verb != "grid" else canonical_run.vor
        extra = {"line": line} if verb == "section" else {}
        drawn = getattr(plot, verb)(subject, backend="mpl", **extra)
    else:
        namespace = _namespace(scope)(
            canonical_run if scope == "model" else canonical_run.vor
        )
        if verb == "section":
            drawn = (namespace.section(line, backend="mpl") if scope == "grid"
                     else namespace.section(line=line, backend="mpl"))
        else:
            drawn = getattr(namespace, verb)(backend="mpl")
    assert isinstance(drawn, Figure), (
        f"{scope}.{verb}(backend='mpl') returned {type(drawn).__name__}, not a Figure"
    )


# --- documented VALUES, not just documented names -----------------------------
def _documented_choices(func, parameter: str) -> list[str]:
    """The `{'a', 'b'}` set a NumPy Parameters entry declares for `parameter`."""

    doc = inspect.getdoc(func) or ""
    match = re.search(rf"^{parameter} : \{{(.+?)\}}", doc, re.M)
    assert match, f"{func.__name__} does not document a value set for {parameter}"
    return [v.strip().strip("'\"") for v in match.group(1).split(",")]


#: What each documented value set needs alongside it to be exercised, and why
#: `type` is not here: `'conc'`/`'temp'` need a GWT/GWE model, which the GWF
#: canonical fixture is not, so the set cannot be swept from this test.
_VALUE_SET_EXTRAS = {
    "contour_method": {"contours": True},
    "select_style": {"select": [0, 1]},
    "per_timestep": {},
    "hover_layers": {},
    "backend": {},
}


@pytest.mark.parametrize("parameter", sorted(_VALUE_SET_EXTRAS))
def test_every_documented_value_is_one_the_code_accepts(parameter, canonical_run):
    """`test_every_signature_parameter_is_documented` checks NAMES. This checks
    the values inside them, which is where the same defect hid one level down.

    Measured on 2026-09-02: `contour_method` documented
    `{'linear', 'cubic', 'nearest'}` and `'nearest'` raised
    "contour method must be 'linear' or 'cubic'" -- it has never been
    implemented. The docstring is spliced onto all eighteen noun verbs, so one
    wrong value in one entry became eighteen wrong entries.
    """

    import matplotlib

    matplotlib.use("Agg")

    rejected = []
    for value in _documented_choices(plot.map, parameter):
        kwargs = {parameter: value, **_VALUE_SET_EXTRAS[parameter]}
        try:
            picture = plot.map(canonical_run, layer=0, **kwargs)
            # Choro is lazy; the validation lives in the render.
            if hasattr(picture, "fig"):
                picture.fig
        except (ValueError, TypeError) as error:
            rejected.append(f"{value!r}: {error}")
    assert not rejected, (
        f"plot.map documents {parameter} values it rejects: {rejected}"
    )


def test_contour_resolution_is_documented_as_cubic_only(canonical_run):
    """It is inert under the DEFAULT method, and the docstring has to say so.

    `_linear_contour_segments` triangulates the cell centres and does not take a
    resolution at all; only `_cubic_contour_segments` interpolates onto a
    resolution-square grid. Measured: 338 contour points at 40, 150 and 400
    under `linear`, and a count that actually moves under `cubic`.
    """

    def points(**kwargs):
        picture = plot.map(canonical_run, layer=0, contours=True, **kwargs)
        return sum(
            len(getattr(trace, "lon", []) or [])
            for trace in picture.fig.data
            if trace.type == "scattermap"
        )

    assert points(contour_resolution=40) == points(contour_resolution=400), (
        "contour_resolution now acts under the linear method -- update the "
        "docstring, which says it does not"
    )
    assert "Only acts under" in (inspect.getdoc(plot.map) or ""), (
        "contour_resolution's entry must say which method it applies to"
    )


# --- section(fill=): the cells, not a line profile -----------------------------
def _mid_line(model):
    """A West-East line across the middle of a model's grid."""

    import shapely as shp

    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    middle = (ymin + ymax) / 2
    return shp.LineString([(xmin + 1, middle), (xmax - 1, middle)])


def test_a_filled_section_draws_the_cells_and_the_line_profile_does_not(canonical_run):
    """`fill=` is the picture the verb could not draw before (2026-09-02).

    `section` answered "what is the head along this line" -- two traces, no
    geometry. The picture people mean by "a cross-section of the model" is the
    cells: the grid, its layers, and a field on them. That existed only as
    `myflopy.modflow.mf6.plot_model_cross_section`, reachable by deep import and
    from no scope at all, which is how it was lost twice.

    Counting artists is the falsifiable form: a profile has LINES and no filled
    collections; a filled section has both.
    """

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)

    profile = canonical_run.plot.section(line=line, backend="mpl").axes[0]
    assert not profile.collections, "the line profile drew filled cells"
    assert profile.lines, "the line profile drew no profile"

    filled = canonical_run.plot.section(line=line, fill="layer", backend="mpl").axes[0]
    assert filled.collections, "fill='layer' drew no cells"
    # ... and the RESULTS are still there, as the water surface over the geology.
    assert filled.lines, "fill='layer' dropped the simulated head surface"


def test_a_results_fill_paints_the_cells_and_earns_a_colorbar(canonical_run):
    """`fill='results'` colours the cells by the model's own field.

    Kind-neutral by construction -- it reads through `_field_reader`, so it is
    heads on GWF, concentration on GWT, temperature on GWE. The colorbar is the
    observable difference from `fill='layer'`, whose discrete layer palette gets
    a legend instead.
    """

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)

    layered = canonical_run.plot.section(line=line, fill="layer", backend="mpl")
    painted = canonical_run.plot.section(
        line=line, fill="results", per=1, backend="mpl", fill_label="head (ft)"
    )
    assert len(painted.axes) > len(layered.axes), "fill='results' drew no colorbar"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"fill": "layer"}, "backend='mpl'"),
        ({"fill": "layer", "backend": "mpl", "interpolate": True}, "line PROFILE"),
        ({"fill": "layer", "backend": "mpl", "spacing": 99}, "line PROFILE"),
        ({"fill": "nope", "backend": "mpl"}, "fill must be"),
    ],
)
def test_fill_refuses_what_it_cannot_honour(kwargs, message, canonical_run):
    """Every one of these was silently accepted somewhere in this family before.

    A `fill=` that quietly returned a Plotly line profile, or that took
    `interpolate=True` and ignored it, would be the same defect the whole 8.8
    pass exists to end -- a parameter an editor completes and the code discards.
    """

    with pytest.raises(ValueError, match=message):
        canonical_run.plot.section(line=_mid_line(canonical_run), **kwargs)


def test_a_bare_grid_fills_by_layer_but_has_no_results_to_paint(canonical_run):
    """A grid answers `fill='layer'` and refuses `fill='results'`, naming why."""

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)
    assert canonical_run.vor.plot.section(line, fill="layer", backend="mpl").axes[0].collections

    with pytest.raises(ValueError, match="bare grid has"):
        canonical_run.vor.plot.section(line, fill="results", backend="mpl")


def test_a_flat_array_fill_is_refused_rather_than_broadcast(canonical_run):
    """One value per cell is ambiguous in section -- it would paint every layer
    the same and look like a real answer."""

    import matplotlib

    matplotlib.use("Agg")
    flat = canonical_run.hds.array(layer=0)
    with pytest.raises(ValueError, match="PER LAYER"):
        canonical_run.plot.section(
            line=_mid_line(canonical_run), fill=flat, backend="mpl"
        )


def _drawn_cells(figure):
    """How many section cells actually carry a value (the rest are masked)."""

    import numpy as np

    array = figure.axes[0].collections[-1].get_array()
    return int(np.sum(~np.ma.getmaskarray(array))) if array is not None else -1


def _legend_labels(figure):
    legend = figure.axes[0].get_legend()
    return [text.get_text() for text in legend.get_texts()] if legend else []


def test_layers_selects_which_cells_are_drawn_and_crops_to_them(canonical_run):
    """`layers=` must act. It was accepted and IGNORED on the first cut.

    `layer=` already existed on the profile branch meaning "overlay these
    layers' head profiles", so it rode the filled branch doing nothing -- the
    exact defect (a parameter an editor completes and the code discards) that
    this whole family of work exists to end. The filled branch gets its own
    name and refuses the profile one.

    Cropping is part of the behaviour, not decoration: masking alone leaves the
    two layers you asked for in a thin band of a full-height axis.
    """

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)
    section = functools.partial(
        canonical_run.plot.section, line=line, fill="layer", backend="mpl"
    )

    everything = section()
    subset = section(layers=[0, 1])
    assert _drawn_cells(subset) < _drawn_cells(everything)
    assert len(_legend_labels(subset)) < len(_legend_labels(everything))

    full_height = everything.axes[0].get_ylim()
    cropped = subset.axes[0].get_ylim()
    assert cropped[0] > full_height[0], "the view was not cropped to the drawn layers"

    with pytest.raises(ValueError, match="layers= to choose"):
        section(layer=2)
    with pytest.raises(ValueError, match="outside this model"):
        section(layers=[0, 99])


def test_head_layers_chooses_whose_water_levels_are_drawn(canonical_run):
    """One water table is the default; several are drawn labelled, none is legal."""

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)
    section = functools.partial(
        canonical_run.plot.section, line=line, fill="layer", backend="mpl"
    )

    one = section()
    several = section(head_layers=[0, 2])
    none = section(head_layers=None)

    assert len(none.axes[0].lines) == 0
    assert len(several.axes[0].lines) == 2 * len(one.axes[0].lines)
    # ... and each is named, because several unlabelled lines are unreadable.
    labels = _legend_labels(several)
    assert sum("Head" in text for text in labels) == 2
    assert _legend_labels(one).count("Simulated Head") == 1


def test_layer_labels_name_the_units_and_default_to_the_models_own(canonical_run):
    """The legend should say "sand", not "Layer 2", when anything knows better.

    Explicit names win; otherwise they come from the model's build context
    (`ModelContext(surfaces=stack.build(vor))` carries them through
    `ModelSpec.build`). The canonical model is built imperatively and has no
    context, so it exercises the documented fallback.
    """

    import matplotlib

    matplotlib.use("Agg")
    from myflopy.modflow.mf6.cross_section_plotting import layer_labels_from_model

    line = _mid_line(canonical_run)
    named = canonical_run.plot.section(
        line=line, fill="layer", backend="mpl", layers=[0, 1],
        layer_labels=["sand", "clay", "till", "rock"],
    )
    assert _legend_labels(named)[:2] == ["sand", "clay"]

    assert layer_labels_from_model(canonical_run, canonical_run.gwf.modelgrid.nlay) is None
    fallback = canonical_run.plot.section(line=line, fill="layer", backend="mpl")
    assert _legend_labels(fallback)[0] == "Layer 1"


def test_the_filled_only_arguments_are_refused_on_the_profile(canonical_run):
    """They describe cells; a line profile has none, so they must not be ignored."""

    for kwargs in ({"layers": [0]}, {"layer_labels": ["a"]}, {"head_layers": [0, 1]}):
        with pytest.raises(ValueError, match="filled section draws"):
            canonical_run.plot.section(line=_mid_line(canonical_run), **kwargs)


def test_layers_hides_the_excluded_cells_outlines_too(canonical_run):
    """Excluding a layer must remove its CELLS, not just its fill.

    Reported as "only the layers I list are filled, but cells for all layers
    still draw". FloPy's `plot_grid` outlines every cell and takes no layer
    filter, so the first cut masked the fill and left the outlines behind.

    Cropping the y-axis cannot cover for it, which is why this is a separate
    test from the cropping one: a layer's elevation range overlaps its
    neighbours'. Measured on the canonical model, the bottom layer spans
    -1.3..59.8 while the retained band is 28.5..166.0 -- so it is drawn straight
    through the middle of the picture no matter where the limits sit.
    """

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)
    nlay = canonical_run.gwf.modelgrid.nlay
    section = functools.partial(
        canonical_run.plot.section, line=line, fill="layer", backend="mpl"
    )

    def _patches(figure):
        return sum(len(c.get_paths()) for c in figure.axes[0].collections)

    everything = _patches(section())
    without_bottom = _patches(section(layers=list(range(nlay - 1))))
    # Every cell of the excluded layer is gone -- outline included. The full
    # picture draws each cell twice (grid + fill), the subset once with edges.
    assert without_bottom < everything / 2, (
        f"{without_bottom} patches for {nlay - 1} of {nlay} layers, against "
        f"{everything} for all of them -- the excluded cells are still drawn"
    )

    one_layer = _patches(section(layers=[1]))
    assert one_layer * (nlay - 1) == pytest.approx(without_bottom, rel=0.05)


def test_show_grid_is_named_and_leftover_arguments_are_refused(canonical_run):
    """`show_grid=False` used to vanish into the `**kwargs` tail.

    The filled branch builds no class that takes an open tail -- that is the
    profile branch's `XSection` -- so anything left in it was accepted and
    dropped. Found by a probe passing `show_grid=False` and getting edges
    anyway.
    """

    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    line = _mid_line(canonical_run)
    section = functools.partial(
        canonical_run.plot.section, line=line, fill="layer", backend="mpl",
        layers=[0, 1],
    )

    def _edge_width(figure):
        filled = [c for c in figure.axes[0].collections if c.get_array() is not None][-1]
        return float(np.atleast_1d(filled.get_linewidth())[0])

    # Asserted on the LINE WIDTH the code sets, read off the style rather than
    # written as a literal. The obvious check -- does the collection carry an
    # edge colour -- measures rcParams, not this code: seaborn's theme (applied
    # by `mpl_axes`, so by whichever test ran first) forces one on regardless,
    # and the assertion passed alone and failed in the full suite.
    from myflopy.modflow.mf6.cross_section_plotting import ModelCrossSectionStyle

    styled = ModelCrossSectionStyle().grid_linewidth
    assert _edge_width(section()) == pytest.approx(styled), "cell edges were not drawn"
    assert _edge_width(section(show_grid=False)) != pytest.approx(styled), (
        "show_grid=False still styled the cell edges"
    )

    with pytest.raises(ValueError, match="not arguments of a filled section"):
        section(interpolate_me=True)


def test_a_subset_section_fits_its_axis_to_what_is_drawn(canonical_run):
    """"The axes need to adjust to the new extent" -- and from the DRAWING.

    The first cut cropped using the grid's own `top`/`botm`, which are
    whole-grid statistics: the limits came from cells the section line never
    crosses. Measured on the canonical model, `layers=[0]` gave an axis of
    66.3..164.2 for content spanning 72.0..146.2 -- a quarter of the height
    empty, which reads as "the axis still thinks the other layers are there".

    Water surfaces count as drawn content: they are `Line2D`, not collections,
    and leaving them out cropped `layers=[3]` with the default `head_layers=0`
    to an axis that excluded the very line its legend advertised.
    """

    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    line = _mid_line(canonical_run)
    nlay = canonical_run.gwf.modelgrid.nlay

    def _fit(**kwargs):
        """(fraction of the axis that is empty, everything fits) for one figure."""

        figure = canonical_run.plot.section(
            line=line, fill="layer", backend="mpl", **kwargs
        )
        axes = figure.axes[0]
        drawn = [p.vertices[:, 1] for c in axes.collections for p in c.get_paths()]
        drawn += [np.asarray(ln.get_ydata(), dtype=float) for ln in axes.lines]
        values = np.concatenate(drawn)
        values = values[np.isfinite(values)]
        low, high = axes.get_ylim()
        return 1 - (values.max() - values.min()) / (high - low), (
            low <= values.min() and values.max() <= high
        )

    for kwargs in ({"layers": [0]}, {"layers": list(range(nlay - 1))},
                   {"layers": [nlay - 1]},
                   {"layers": [nlay - 1], "head_layers": None}):
        empty, fits = _fit(**kwargs)
        assert fits, f"{kwargs} clipped part of what it drew"
        # The 5% pad on each side and nothing more.
        assert empty == pytest.approx(0.09, abs=0.02), (
            f"{kwargs} left {empty:.0%} of the axis empty"
        )


def test_legend_placement_actually_moves_the_legend(canonical_run):
    """`legend=` must place it, not merely be accepted.

    Positions are compared against each other rather than to fixed pixels: the
    figure size and the theme both come from elsewhere, so an absolute
    coordinate would pin the wrong thing. What must hold is that "bottom" is
    below "top", "left" is left of "right", and "outside right" leaves the axes
    entirely -- which is the placement a section actually needs, since a layer
    legend anywhere inside sits on top of the geology it describes.
    """

    import matplotlib

    matplotlib.use("Agg")
    line = _mid_line(canonical_run)

    def _centre(legend):
        figure = canonical_run.plot.section(
            line=line, fill="layer", backend="mpl", legend=legend
        )
        figure.canvas.draw()
        drawn = figure.axes[0].get_legend()
        assert drawn is not None, f"legend={legend!r} drew none"
        box = drawn.get_window_extent()
        axes_box = figure.axes[0].get_window_extent()
        return box, axes_box

    top, _ = _centre("top")
    bottom, _ = _centre("bottom")
    assert bottom.y0 < top.y0, "'bottom' was not below 'top'"

    left, _ = _centre("left")
    right, _ = _centre("right")
    assert left.x0 < right.x0, "'left' was not left of 'right'"

    outside, axes_box = _centre("outside right")
    assert outside.x0 >= axes_box.x1 - 1, "'outside right' stayed inside the axes"

    below, axes_box = _centre("outside bottom")
    assert below.y1 <= axes_box.y0 + 1, "'outside bottom' stayed inside the axes"

    assert canonical_run.plot.section(
        line=line, fill="layer", backend="mpl", legend=False
    ).axes[0].get_legend() is None

    for bad in ("nowhere", "outside sideways"):
        with pytest.raises(ValueError):
            canonical_run.plot.section(
                line=line, fill="layer", backend="mpl", legend=bad
            )
