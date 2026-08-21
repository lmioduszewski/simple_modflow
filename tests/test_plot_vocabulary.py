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
