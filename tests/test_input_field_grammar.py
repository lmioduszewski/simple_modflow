"""Tests for inputs.<field> nodes + field= sugar across the three surfaces
(Phase 4 of the unified view grammar).

Single-model input namespaces are FieldMappable (fields are the nodes); the
group/diff input accessors are LeafFieldSugar (the accessor is itself a leaf
drawing the default field, with pinned field nodes alongside).
"""

from __future__ import annotations

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from myflopy.modflow.mf6.package_inputs import (
    CellPackageInputFieldExplorer,
    CellPackageInputsExplorer,
    UzfInputsNamespace,
)
from myflopy.modflow.mf6.package_registry import get_package_input_field_names
from myflopy.project.model_diff import PackageDiff
from myflopy.project.model_group import (
    GroupPackageInputField,
    GroupPackageInputs,
    ModelGroup,
)


class _FakeModel:
    def __init__(self, name):
        self.name = name
        self.package_names = []


def _group(names=("a", "b"), reference="a"):
    return ModelGroup({n: _FakeModel(n) for n in names}, reference=reference)


# --- registry field lists -------------------------------------------------------
def test_registry_lists_input_field_names():
    assert get_package_input_field_names("ghb") == ["bhead", "cond"]
    assert get_package_input_field_names("drn") == ["elev", "cond"]
    assert get_package_input_field_names("rch") == ["recharge"]
    assert "finf" in get_package_input_field_names("uzf")
    assert get_package_input_field_names("not_a_package") == []


# --- single model: FieldMappable namespaces --------------------------------------
def test_single_model_inputs_field_dispatch():
    inputs = CellPackageInputsExplorer(model=None, package_name="ghb")
    assert inputs.field_names() == ["bhead", "cond"]
    assert inputs._default_field == "bhead"

    cond = inputs._field_accessor("cond")
    assert isinstance(cond, CellPackageInputFieldExplorer)
    assert cond.field_name == "cond"
    # default field when field=None
    assert inputs._field_accessor(None).field_name == "bhead"
    with pytest.raises(ValueError):
        inputs._field_accessor("head")  # chd field, not ghb
    with pytest.raises(AttributeError):
        _ = inputs.not_a_field


def test_uzf_inputs_namespace_field_dispatch():
    inputs = UzfInputsNamespace(model=None)
    assert inputs._default_field == "finf"
    assert set(inputs.field_names()) >= {"finf", "pet", "extdp"}
    assert inputs._field_accessor("pet").field_name == "pet"
    with pytest.raises(ValueError):
        inputs._field_accessor("stage")


# --- group: LeafFieldSugar on the aggregator leaf --------------------------------
def test_group_inputs_field_nodes_pin_the_value_column():
    accessor = GroupPackageInputs(_group(), "ghb")
    assert accessor.field_names() == ["bhead", "cond"]

    node = accessor.cond  # attribute sugar
    assert isinstance(node, GroupPackageInputField)
    assert node.field_name == "cond" and node.package_name == "ghb"

    captured = {}
    accessor._spatial_map = lambda **kwargs: captured.update(kwargs) or "CHORO"
    fresh = GroupPackageInputField(accessor, "cond")
    assert fresh._spatial_map(per=3, model="b") == "CHORO"
    assert captured["value_column"] == "cond"  # pinned
    assert captured["per"] == 3 and captured["model"] == "b"

    with pytest.raises(ValueError):
        accessor.field("head")
    with pytest.raises(AttributeError):
        _ = accessor.not_a_field  # not RecursionError


def test_group_inputs_field_kwarg_dispatches_to_the_node():
    accessor = GroupPackageInputs(_group(), "ghb")
    seen = {}

    class _Node:
        def map(self, *args, **kwargs):
            seen["args"] = args
            seen["kwargs"] = kwargs
            return "FIG"

    def _node_factory(name):
        seen["field"] = name
        return _Node()

    accessor._field_node = _node_factory
    assert accessor.map("b", field="cond", per=2) == "FIG"
    assert seen["field"] == "cond"
    assert seen["args"] == ("b",) and seen["kwargs"] == {"per": 2}


# --- diff: PackageDiff field pinning ----------------------------------------------
def test_package_diff_field_nodes_pin_value_column_and_diff_series():
    diff = _group().diff()
    pkg = diff.packages.ghb.inputs
    assert pkg.field_names() == ["bhead", "cond"]

    node = pkg.cond
    assert isinstance(node, PackageDiff)
    assert node.field_name == "cond"
    assert node._spatial_value_label() == "cond"

    captured = {}
    node._inputs.compare_map = lambda **kwargs: captured.update(kwargs) or "FIG"
    assert node.map("b") == "FIG"
    assert captured["value_column"] == "cond"
    assert captured["model_name"] == "b"

    # series hook draws the aligned *_diff column
    frame = pd.DataFrame(columns=["model", "per", "cell", "cond", "cond_diff"])
    assert node._series_value_column(frame) == "cond_diff"
    with pytest.raises(KeyError):
        pkg.bhead._series_value_column(frame)  # bhead_diff not present

    with pytest.raises(ValueError):
        pkg.field("head")


# --- field-name parity: one canonical list per package across surfaces -----------
def test_field_name_parity_across_surfaces():
    """Single-model, group, and diff namespaces advertise identical field names."""

    from myflopy.modflow.mf6.package_surface_water import (
        LakResultsNamespace,
        SfrResultsNamespace,
    )
    from myflopy.project.model_group import (
        GroupLakResultsNamespace,
        GroupSfrResultsNamespace,
    )
    from myflopy.project.model_results_diff import (
        LakResultsDiffNamespace,
        SfrResultsDiffNamespace,
        UzfResultsDiffNamespace,
    )

    group = _group()
    diff = group.diff()

    # results fields (q + stage everywhere; the diff-only 'flow' alias is gone)
    assert LakResultsNamespace(None)._field_names() == ["q", "stage"]
    assert GroupLakResultsNamespace(None)._field_names() == ["q", "stage"]
    assert LakResultsDiffNamespace(diff)._field_names() == ["q", "stage"]

    assert SfrResultsNamespace(None)._field_names() == ["q", "stage"]
    assert GroupSfrResultsNamespace(None)._field_names() == ["q", "stage"]
    assert SfrResultsDiffNamespace(diff)._field_names() == ["q", "stage"]

    assert UzfResultsDiffNamespace(diff)._field_names() == ["gwrch", "sat"]

    # input fields come from one registry on all three surfaces
    for pkg in ("ghb", "drn", "chd", "wel", "rch"):
        single = CellPackageInputsExplorer(None, pkg).field_names()
        grouped = GroupPackageInputs(group, pkg).field_names()
        diffed = diff.packages.ghb.inputs.field_names() if pkg == "ghb" else None
        assert single == grouped == get_package_input_field_names(pkg)
        if diffed is not None:
            assert diffed == single


# --- real renders on the canonical model -------------------------------------------
@pytest.mark.slow
def test_input_field_grammar_renders_on_canonical(canonical_run):
    """inputs.<field> nodes and field= sugar render on all three surfaces."""

    import myflopy as mf
    import plotly.graph_objects as go
    from matplotlib.figure import Figure

    model = canonical_run

    # single model: default field, explicit node, field= sugar -- same content
    default_map = model.packages.ghb.inputs.map()
    node_map = model.packages.ghb.inputs.cond.map()
    sugar_map = model.packages.ghb.inputs.map(field="cond")
    for panel in (default_map, node_map, sugar_map):
        assert panel is not None and not isinstance(panel, Figure)
    assert list(node_map.zs) == list(sugar_map.zs)  # same field drawn

    cond_series = model.packages.ghb.inputs.cond.plot()
    assert isinstance(cond_series, go.Figure)

    uzf_map = model.packages.uzf.inputs.map(field="pet")
    assert uzf_map is not None
    sfr_map = model.packages.sfr.inputs.map(field="rhk")
    assert sfr_map is not None

    # group surface
    reloaded = mf.load_mf6_run(model.workspace)
    group = ModelGroup({"a": model, "b": reloaded}, reference="a")
    group_cond = group.packages.ghb.inputs.cond.map("b")
    assert group_cond is not None and not isinstance(group_cond, Figure)
    group_mosaic = group.packages.ghb.inputs.mosaic(
        field="cond", by="model", backend="mpl"
    )
    assert isinstance(group_mosaic, Figure)

    # diff surface (a == b, so flat Δ -- but every path must build)
    diff = group.diff()
    diff_cond = diff.packages.ghb.inputs.cond.map("b")
    assert diff_cond is not None
    diff_series = diff.packages.ghb.inputs.plot(field="cond")
    assert isinstance(diff_series, go.Figure)
