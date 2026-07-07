"""Tests for the grouped map API: deprecated-shortcut hiding, the ``models=``
mosaic option on inputs/results maps, positional model names, and the
``<pkg>.results.map`` namespace passthrough.
"""

from __future__ import annotations

import matplotlib
import pytest

matplotlib.use("Agg")  # headless render for the mosaic smoke test

from myflopy.project.model_group import (
    GroupCellPackageResultsNamespace,
    GroupLakResultsNamespace,
    GroupLakStageResults,
    GroupPackageInputs,
    GroupSfrResultsNamespace,
    GroupSfrStageResults,
    ModelGroup,
    _GroupSpatialView,
)


class _FakeModel:
    def __init__(self, name):
        self.name = name
        self.package_names = []


def _group(names=("a", "b", "c"), reference="a"):
    return ModelGroup({n: _FakeModel(n) for n in names}, reference=reference)


# --- deprecated flat shortcuts hidden but still functional --------------------
def test_deprecated_shortcuts_hidden_from_dir():
    group = _group()
    names = set(dir(group))
    assert not ({"rch", "chd", "drn", "ghb", "wel", "uzf"} & names)  # hidden
    assert "packages" in names  # the preferred path stays visible


def test_deprecated_shortcuts_still_work_with_warning():
    group = _group()
    with pytest.warns(DeprecationWarning):
        accessor = group.rch
    assert accessor is group._rch
    assert isinstance(accessor, GroupPackageInputs)


# --- mosaic(by="model") facets over the group's models -----------------------
def test_inputs_mosaic_facets_over_models():
    accessor = GroupPackageInputs(_group(names=("a", "b", "c")), "rch")
    accessor._spatial_map = lambda **kwargs: f"panel-{kwargs.get('model')}"

    panels = accessor._facet_panels("model", per=0, layer=0, model=None)
    assert [label for label, _ in panels] == ["a", "b", "c"]
    assert [choro for _, choro in panels] == ["panel-a", "panel-b", "panel-c"]


def test_group_inputs_use_unified_grammar_not_old_verbs():
    import inspect

    accessor = GroupPackageInputs(_group(), "rch")
    for verb in ("map", "mosaic", "animate"):
        assert callable(getattr(accessor, verb))
    assert not hasattr(accessor, "subplot_map")  # replaced by mosaic(by="model")
    assert "models" not in inspect.signature(accessor.map).parameters  # replaced by mosaic


# --- <pkg>.results.map namespace passthrough ---------------------------------
class _FakeResultAccessor:
    def map(self, model_name=None, **kwargs):
        return ("map", model_name, kwargs)

    def mosaic(self, **kwargs):
        return ("mosaic", kwargs)

    def animate(self, **kwargs):
        return ("animate", kwargs)


def test_results_namespace_field_dispatch():
    namespace = GroupCellPackageResultsNamespace(_FakeResultAccessor())
    assert namespace.field_names() == ["q"]

    # default field q -> .q.map, with a positional model forwarded
    verb, model_name, kwargs = namespace.map("F9b", per=3)
    assert verb == "map" and model_name == "F9b" and kwargs["per"] == 3

    assert namespace.mosaic(by="model")[0] == "mosaic"  # inherited grammar
    with pytest.raises(ValueError):
        namespace.map(field="not_a_field")


# --- LAK/SFR grouped results expose a spatially mappable ``stage`` field ------
class _StageResultAccessor:
    """Fake result accessor carrying a group (stage props read ``.group``)."""

    def __init__(self, group):
        self.group = group


def test_lak_sfr_grouped_results_expose_stage_field():
    """Grouped LAK/SFR results advertise ``stage`` alongside ``q`` and the stage
    accessor is a full :class:`SpatialView` (map/mosaic/animate), matching the
    single-model ``lak.results`` / ``sfr.results`` grammar."""

    group = _group(names=("a", "b"))
    lak_ns = GroupLakResultsNamespace(_StageResultAccessor(group))
    sfr_ns = GroupSfrResultsNamespace(_StageResultAccessor(group))

    assert lak_ns.field_names() == ["q", "stage"]
    assert sfr_ns.field_names() == ["q", "stage"]

    # field='stage' dispatches to a grouped stage SpatialView (not q).
    assert isinstance(lak_ns._field_accessor("stage"), GroupLakStageResults)
    assert isinstance(sfr_ns._field_accessor("stage"), GroupSfrStageResults)
    with pytest.raises(ValueError):
        lak_ns._field_accessor("not_a_field")

    # The stage accessors inherit the unified grammar and face over the models.
    for stage in (lak_ns.stage, sfr_ns.stage):
        assert isinstance(stage, _GroupSpatialView)
        for verb in ("map", "mosaic", "animate"):
            assert callable(getattr(stage, verb))
        assert not hasattr(stage, "subplot_map")
        stage._spatial_map = lambda **kwargs: f"stage-{kwargs.get('model')}"
        panels = stage._facet_panels("model", per=0, layer=0, model=None)
        assert [label for label, _ in panels] == ["a", "b"]
        assert [choro for _, choro in panels] == ["stage-a", "stage-b"]


# --- real render smoke (matplotlib mosaic on a real group) -------------------
@pytest.mark.slow
def test_group_inputs_results_mosaic_render_on_canonical(canonical_run):
    """A two-model group (run vs reload) renders a shared-scale per-model mosaic
    for both an input and a result field via ``mosaic(by="model")``."""

    import myflopy as mf
    from matplotlib.figure import Figure

    reloaded = mf.load_mf6_run(canonical_run.workspace)
    group = ModelGroup({"a": canonical_run, "b": reloaded}, reference="a")

    inputs_fig = group.packages.rch.inputs.mosaic(by="model", backend="mpl")
    assert isinstance(inputs_fig, Figure)

    results_fig = group.packages.ghb.results.mosaic(by="model", backend="mpl")
    assert isinstance(results_fig, Figure)

    # a single-model input map is still the interactive (Plotly) Choro
    single = group.packages.rch.inputs.map()  # reference default
    assert single is not None and not isinstance(single, Figure)


@pytest.mark.slow
def test_group_lak_sfr_stage_field_mosaic_render_on_canonical(canonical_run):
    """Grouped ``lak.results``/``sfr.results`` map the ``stage`` field spatially:
    ``mosaic(field="stage", by="model")`` renders a per-model mpl mosaic and
    ``map(field="stage")`` returns the interactive Choro -- the same grammar the
    single-model stage explorer exposes (regression for the ``Unknown field
    'stage'`` gap)."""

    import myflopy as mf
    from matplotlib.figure import Figure

    reloaded = mf.load_mf6_run(canonical_run.workspace)
    group = ModelGroup({"a": canonical_run, "b": reloaded}, reference="a")

    for package in ("lak", "sfr"):
        results = getattr(group.packages, package).results
        assert "stage" in results.field_names()

        stage_mosaic = results.mosaic(field="stage", by="model", backend="mpl")
        assert isinstance(stage_mosaic, Figure)

        stage_map = results.map(field="stage")  # reference model, interactive
        assert stage_map is not None and not isinstance(stage_map, Figure)
