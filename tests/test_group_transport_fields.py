"""``group.conc`` / ``group.temp`` and ``diff().conc`` / ``.temp`` (plan §6.1/6.2 item 4).

The grouped transport twins of ``group.hds``. They are not clones: §6.2 says to
fix the abstraction rather than copy-paste, so ``GroupHeads`` was refactored onto
a shared ``_GroupFieldView`` and all three kinds are now the same code with five
class attributes each.

Fixtures build TWO coupled transport runs that differ in one physical parameter
(porosity), because a grouped comparison over identical models can only ever
prove the plumbing runs, never that it computes the right difference.
"""

from __future__ import annotations

import numpy as np
import pytest

from myflopy.project.group.conc import GroupConc
from myflopy.project.group.spatial import GroupHeads, _GroupFieldView
from myflopy.project.group.temp import GroupTemp
from myflopy.project.model_group import ModelGroup
from tests.test_gwt_gwe_results import _coupled_run


def _transport_group(tmp_path_factory, kind: str) -> ModelGroup:
    """Two coupled transport models differing only in porosity, as a group."""

    root = tmp_path_factory.mktemp(f"group_{kind}")
    reference = _coupled_run(root / "a", kind).model("trans")
    variant = _coupled_run(root / "b", kind, porosity=0.05).model("trans")
    return ModelGroup({"a": reference, "b": variant}, reference="a")


@pytest.fixture(scope="module")
def gwt_group(tmp_path_factory):
    return _transport_group(tmp_path_factory, "gwt")


@pytest.fixture(scope="module")
def gwe_group(tmp_path_factory):
    return _transport_group(tmp_path_factory, "gwe")


def test_all_three_field_views_are_one_class_with_different_knobs():
    """The refactor's whole point: no third copy of the same hundred lines.

    §6.2 says that if this phase is not mostly mechanical, fix the §6.0
    abstraction instead of copy-pasting. So the three grouped field accessors
    share ``_GroupFieldView`` and differ ONLY in their class attributes.
    """

    for view in (GroupHeads, GroupConc, GroupTemp):
        assert issubclass(view, _GroupFieldView)
        # no view re-implements the shared table/compare logic
        for shared in ("get", "compare", "compare_map", "_spatial_map"):
            assert shared not in vars(view), f"{view.__name__} re-implements {shared}"

    knobs = {
        view.__name__: (
            view.reader_attribute,
            view.table_attribute,
            view.value_column,
            view.value_label,
        )
        for view in (GroupHeads, GroupConc, GroupTemp)
    }
    assert knobs == {
        "GroupHeads": ("hds", "all_heads", "elev", "head"),
        "GroupConc": ("conc", "all_conc", "conc", "conc"),
        "GroupTemp": ("temp", "all_temp", "temp", "temp"),
    }


@pytest.mark.slow
def test_the_grouped_transport_tables_carry_the_field_column(gwt_group, gwe_group):
    """``get()`` names its column after the field, not after heads."""

    for group, column in ((gwt_group, "conc"), (gwe_group, "temp")):
        accessor = getattr(group, column)
        frame = accessor.get()
        assert list(frame.columns) == ["model", "kstpkper", "layer", "cell", column]
        assert set(frame["model"]) == {"a", "b"}
        assert frame[column].notna().any()


@pytest.mark.slow
def test_the_grouped_difference_is_model_minus_reference(gwt_group, gwe_group):
    """Orientation, computed independently -- the assertion a sign flip cannot pass.

    Nothing pinned this even for heads: the shared colorscale helper is tested in
    isolation, so reversing the subtraction inverts every difference map while
    ``test_diff_maps_use_rdbu_negative_red_positive_blue`` still passes. Here the
    expected difference is recomputed from ``get()`` and compared element-wise.
    """

    for group, column in ((gwt_group, "conc"), (gwe_group, "temp")):
        accessor = getattr(group, column)
        comparison = accessor.compare().sort_values(["per", "layer", "cell"])
        assert not comparison.empty

        raw = accessor.get()
        per_model = {
            name: sub.sort_values(["layer", "cell"])[column].to_numpy(dtype=float)
            for name, sub in raw.groupby("model")
        }
        expected = per_model["b"] - per_model["a"]

        assert np.allclose(
            comparison["diff"].to_numpy(dtype=float), expected, atol=1e-12
        ), f"{column}: diff is not (model - reference)"
        # the two models really do differ, or the check above proves nothing
        assert np.abs(expected).max() > 1e-9, "porosity knob stopped changing the field"


@pytest.mark.slow
def test_the_grouped_field_answers_the_whole_grammar(gwt_group):
    """map/plot/xs/mosaic come from SpatialView and must work unchanged."""

    from myflopy.viz import Fig

    conc = gwt_group.conc
    assert conc.map() is not None
    assert isinstance(conc.plot(layer=0), Fig)
    # xs was BROKEN on the transport readers until 2026-07-27 (ledger 99):
    # XSection reached for model.hds, which the kind guard refuses on a GWT model
    assert isinstance(conc.section(cells=[0, 5, 10]), Fig)
    assert conc.mosaic(by="model") is not None


@pytest.mark.slow
def test_the_public_diff_verb_reaches_the_transport_fields(gwt_group, gwe_group):
    """``group.diff().conc`` / ``.temp`` -- the ONE public diff route."""

    for group, column, label in (
        (gwt_group, "conc", "Δ conc"),
        (gwe_group, "temp", "Δ temp"),
    ):
        leaf = getattr(group.diff(), column)
        assert type(leaf).__name__ == f"{column.capitalize()}ResultDiff"

        frame = leaf.get()
        assert f"reference_{column}" in frame.columns and "diff" in frame.columns

        choro = leaf.map("b")
        assert label in choro.get_choropleth().hovertemplate

        summary = leaf.summary()
        assert summary.loc[0, "model"] == "b"
        assert summary.loc[0, "max_abs_diff"] > 0  # the two models really differ


@pytest.mark.slow
def test_the_transport_group_refuses_the_wrong_field(gwt_group):
    """A GWT group has no temperature; the kind gate says so rather than guessing."""

    with pytest.raises(AttributeError, match="is a GWT model.*temperature"):
        gwt_group.temp.get()


@pytest.mark.slow
def test_all_conc_and_all_temp_are_kind_gated(gwt_group):
    """Unlike ``all_heads``, the new tables assert their kind.

    ``all_heads`` is ungated for historical reasons; a new transport table should
    not inherit that, or a group holding the wrong kind reads a missing ``.ucn``
    and fails far from the cause.
    """

    model = gwt_group.models["a"]
    assert not model.all_conc.empty

    with pytest.raises(AttributeError, match="is a GWT model.*temperature"):
        model.all_temp
