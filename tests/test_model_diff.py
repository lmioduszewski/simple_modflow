"""Tests for the ModelDiff engine (Phase 1: structural + value tiers).

The fast tests patch ``build_cell_package_input_table`` so the whole diff engine
(structural set difference + value compare + summary + report) runs on small,
controlled tables without building a real MF6 model. One ``slow`` test exercises
the full path end-to-end on real canonical models via the ``model.diff(...)``
door.
"""

from __future__ import annotations

import pandas as pd
import pytest

from myflopy.project import model_diff as md
from myflopy.project.group import inputs as group_inputs
from myflopy.project.model_config import ModelConfig
from myflopy.project.model_group import ModelGroup


def _empty_config() -> ModelConfig:
    return ModelConfig(pd.DataFrame(columns=["section", "setting", "value"]))


class _FakeModel:
    """Minimal stand-in exposing only what the diff engine reads."""

    def __init__(self, name: str, package_names, config: ModelConfig | None = None):
        self.name = name
        self.package_names = list(package_names)
        self.config = config if config is not None else _empty_config()


def _ghb(model_name: str, rows):
    """Build a tidy GHB input table (rows: (per, layer, cell, bhead, cond))."""

    return pd.DataFrame(
        rows, columns=["per", "layer", "cell", "bhead", "cond"]
    ).assign(package="ghb", model=model_name)


@pytest.fixture
def tables(monkeypatch):
    """Patch the input-table builder in both modules; return the table registry."""

    registry: dict[tuple[str, str], pd.DataFrame] = {}

    def stub(model, package_name, *, per=None, layer=None, cells=None):
        key = (model.name, str(package_name).lower())
        if key not in registry:
            # Emulate the real builder raising for a package that isn't attached.
            raise KeyError(f"{package_name} not attached to {model.name}")
        frame = registry[key].copy()
        if per is not None:
            frame = frame[frame["per"] == int(per)]
        if layer is not None:
            wanted = layer if isinstance(layer, (list, tuple, set)) else [layer]
            frame = frame[frame["layer"].isin([int(x) for x in wanted])]
        if cells is not None:
            frame = frame[frame["cell"].isin([int(x) for x in cells])]
        return frame.reset_index(drop=True)

    monkeypatch.setattr(md, "build_cell_package_input_table", stub)
    # Patch where the builder is USED: group/inputs.py (plan 4.2 split).
    monkeypatch.setattr(group_inputs, "build_cell_package_input_table", stub)
    return registry


def _group(tables, model_names, *, reference):
    models = {name: _FakeModel(name, ["GHB"]) for name in model_names}
    return ModelGroup(models, reference=reference)


# --- identical (the faithful-copy happy path) --------------------------------
def test_identical_models_report_no_difference(tables):
    rows = [(0, 0, 10, 5.0, 100.0), (0, 0, 11, 5.0, 100.0)]
    tables[("ref", "ghb")] = _ghb("ref", rows)
    tables[("twin", "ghb")] = _ghb("twin", rows)

    diff = _group(tables, ["ref", "twin"], reference="ref").diff()
    summary = diff.summary()

    assert list(summary["model"]) == ["twin"]
    row = summary.iloc[0]
    assert row["identical"]
    assert row["cells_only_in_reference"] == 0
    assert row["cells_only_in_model"] == 0
    assert row["cells_shared"] == 2
    assert row["value_cells_changed"] == 0
    assert diff.packages.ghb.inputs.cells().empty
    assert "identical to reference" in diff.report()


# --- structural + value tiers together ---------------------------------------
def test_structural_and_value_differences_detected(tables):
    # ref cells 10..13; variant drops 10, adds 14, and changes cond on cell 12.
    tables[("ref", "ghb")] = _ghb(
        "ref",
        [(0, 0, 10, 5.0, 100.0), (0, 0, 11, 5.0, 100.0),
         (0, 0, 12, 5.0, 100.0), (0, 0, 13, 5.0, 100.0)],
    )
    tables[("variant", "ghb")] = _ghb(
        "variant",
        [(0, 0, 11, 5.0, 100.0), (0, 0, 12, 5.0, 999.0),
         (0, 0, 13, 5.0, 100.0), (0, 0, 14, 5.0, 100.0)],
    )

    diff = _group(tables, ["ref", "variant"], reference="ref").diff()

    # structural tier
    cells = diff.packages.ghb.inputs.cells()
    only_ref = set(cells.loc[cells["membership"] == "only_in_reference", "cell"])
    only_model = set(cells.loc[cells["membership"] == "only_in_model", "cell"])
    assert only_ref == {10}
    assert only_model == {14}

    # value tier (aligned x / reference_x / x_diff on shared cells)
    values = diff.packages.ghb.inputs.values()
    changed = values.loc[values["cond_diff"] != 0.0]
    assert set(changed["cell"]) == {12}
    assert float(changed["cond_diff"].iloc[0]) == pytest.approx(899.0)

    # summary rolls both tiers up
    row = diff.summary().iloc[0]
    assert row["cells_only_in_reference"] == 1
    assert row["cells_only_in_model"] == 1
    assert row["cells_shared"] == 3
    assert row["value_cells_changed"] == 1
    assert not row["identical"]


def test_structural_diff_survives_only_in_one_model(tables):
    """The bug the structural tier fixes: a value-only inner join drops the
    cell that exists in only one model; the structural tier surfaces it."""

    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 99, 5.0, 100.0)])

    diff = _group(tables, ["ref", "variant"], reference="ref").diff()
    cells = diff.packages.ghb.inputs.cells()
    assert set(cells.loc[cells["membership"] == "only_in_reference", "cell"]) == {10}
    assert set(cells.loc[cells["membership"] == "only_in_model", "cell"]) == {99}
    # No shared cells -> value compare is empty, but the diff is NOT identical.
    assert not diff.summary().iloc[0]["identical"]


# --- package presence differs ------------------------------------------------
def test_package_absent_in_one_model(tables):
    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0), (0, 0, 11, 5.0, 100.0)])
    # variant has no GHB at all: not in package_names, no table registered.
    models = {"ref": _FakeModel("ref", ["GHB"]), "variant": _FakeModel("variant", [])}
    diff = ModelGroup(models, reference="ref").diff()

    row = diff.summary().iloc[0]
    assert row["present_in_reference"]
    assert not row["present_in_model"]
    assert row["cells_only_in_reference"] == 2
    assert row["cells_only_in_model"] == 0
    assert row["value_cells_changed"] == 0
    assert not row["identical"]


# --- reference-star across 3 models ------------------------------------------
def test_reference_star_three_models(tables):
    base = [(0, 0, 10, 5.0, 100.0), (0, 0, 11, 5.0, 100.0)]
    tables[("ref", "ghb")] = _ghb("ref", base)
    tables[("twin", "ghb")] = _ghb("twin", base)
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 10, 5.0, 100.0)])  # dropped 11

    diff = _group(tables, ["ref", "twin", "variant"], reference="ref").diff()
    assert diff.reference == "ref"
    assert diff.model_names == ["twin", "variant"]

    summary = diff.summary().set_index("model")
    assert summary.loc["twin", "identical"]
    assert not summary.loc["variant", "identical"]
    assert summary.loc["variant", "cells_only_in_reference"] == 1


# --- focus, namespace, filters, errors ---------------------------------------
def test_focus_on_one_model(tables):
    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("twin", "ghb")] = _ghb("twin", [(0, 0, 10, 5.0, 100.0)])
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 99, 5.0, 100.0)])

    diff = _group(tables, ["ref", "twin", "variant"], reference="ref").diff()
    focus = diff.model("variant")
    assert set(focus.summary()["model"]) == {"variant"}
    assert not focus.cells("ghb").empty
    assert "`variant`" in focus.report()


def test_per_and_layer_filters(tables):
    tables[("ref", "ghb")] = _ghb(
        "ref", [(0, 0, 10, 5.0, 100.0), (1, 0, 20, 5.0, 100.0)]
    )
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 10, 5.0, 100.0)])  # missing per=1

    diff = _group(tables, ["ref", "variant"], reference="ref").diff()
    per0 = diff.packages.ghb.inputs.cells(per=0)
    assert per0.empty  # identical at period 0
    per1 = diff.packages.ghb.inputs.cells(per=1)
    assert set(per1.loc[per1["membership"] == "only_in_reference", "cell"]) == {20}


def test_package_names_only_lists_present_packages(tables):
    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 10, 5.0, 100.0)])
    diff = _group(tables, ["ref", "variant"], reference="ref").diff()
    assert diff.package_names == ["ghb"]  # rch/chd/drn/wel not attached


def test_unsupported_package_and_self_diff_errors(tables):
    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 10, 5.0, 100.0)])
    diff = _group(tables, ["ref", "variant"], reference="ref").diff()

    with pytest.raises(AttributeError):
        diff.package("nope")
    with pytest.raises(ValueError):
        diff.packages.ghb.inputs.cells(model_name="ref")  # reference vs itself
    with pytest.raises(KeyError):
        diff.packages.ghb.inputs.cells(model_name="ghost")


# --- the unified tree shape ---------------------------------------------------
def test_diff_tree_mirrors_single_and_group_shape(tables):
    """diff.packages.<pkg>.inputs/.results + diff.hds/diff.bud -- one grammar."""

    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("variant", "ghb")] = _ghb("variant", [(0, 0, 10, 5.0, 100.0)])
    diff = _group(tables, ["ref", "variant"], reference="ref").diff()

    node = diff.packages.ghb
    assert node.inputs.package_name == "ghb"  # PackageDiff (inputs tier)
    assert node.results.field_names() == ["q"]

    # advanced nodes: inputs = geometry diff, results = field namespaces
    assert type(diff.packages.lak.inputs).__name__ == "LakConnectionDiff"
    assert diff.packages.lak.results.field_names() == ["q", "stage"]
    assert type(diff.packages.sfr.inputs).__name__ == "SfrReachDiff"
    assert diff.packages.sfr.results.field_names() == ["q", "stage"]
    assert diff.packages.uzf.results.field_names() == ["gwrch", "sat"]
    assert type(diff.packages.mvr.results).__name__ == "MvrResultDiff"

    # heads/budget leaves at the surface root, like model.hds / group.hds
    assert type(diff.hds).__name__ == "HeadsResultDiff"
    assert type(diff.bud).__name__ == "BudgetResultDiff"

    # the old flat results tree is gone outright
    assert not hasattr(diff, "results")

    # the string door returns the same node as attribute access
    assert diff.package("ghb").package_name == "ghb"
    assert diff.package("ghb").inputs.package_name == "ghb"


# --- map / values accept a positional model name -----------------------------
def test_values_accepts_positional_model_name(tables):
    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("A", "ghb")] = _ghb("A", [(0, 0, 10, 6.0, 100.0)])
    tables[("B", "ghb")] = _ghb("B", [(0, 0, 10, 5.0, 100.0)])
    diff = _group(tables, ["ref", "A", "B"], reference="ref").diff()

    values = diff.packages.ghb.inputs.values("A")  # positional -> only compared model A
    assert set(values["model"]) == {"A"}


def test_map_forwards_positional_model_name(tables):
    tables[("ref", "ghb")] = _ghb("ref", [(0, 0, 10, 5.0, 100.0)])
    tables[("A", "ghb")] = _ghb("A", [(0, 0, 10, 6.0, 100.0)])
    tables[("B", "ghb")] = _ghb("B", [(0, 0, 10, 5.0, 100.0)])
    pkg = _group(tables, ["ref", "A", "B"], reference="ref").diff().packages.ghb.inputs

    captured = {}
    pkg._inputs.compare_map = lambda **kwargs: captured.update(kwargs) or "FIG"
    # positional model name must be accepted (was a TypeError) and forwarded
    assert pkg.map("A") == "FIG"
    assert captured["model_name"] == "A"


def test_map_without_model_name_errors_helpfully(tables):
    for name in ("ref", "A", "B"):
        tables[(name, "ghb")] = _ghb(name, [(0, 0, 10, 5.0, 100.0)])
    diff = _group(tables, ["ref", "A", "B"], reference="ref").diff()

    with pytest.raises(ValueError) as excinfo:
        diff.packages.ghb.inputs.map()  # ambiguous: two non-reference models
    message = str(excinfo.value)
    assert "compare_map" not in message      # no internal-method leak
    assert "A" in message and "B" in message  # lists the choices


# --- end-to-end on a real model via the model.diff() door --------------------
@pytest.mark.slow
def test_model_diff_end_to_end_on_canonical(tmp_path):
    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )

    config = CanonicalModelConfig.validation()
    reference = build_canonical_model(tmp_path / "reference", config=config)
    reference.name = "reference"
    variant = build_canonical_model(tmp_path / "variant", config=config)
    variant.name = "variant"

    # Perturb the variant's GHB in period 0: drop one cell (structural) and bump
    # the conductance on another (value).
    ghb = variant.package("ghb")
    spd = ghb.stress_period_data.get_data()
    period0 = spd[0].copy()
    dropped_cellid = tuple(period0["cellid"][0])
    kept = period0[1:].copy()
    kept["cond"][0] = float(kept["cond"][0]) + 12345.0
    spd[0] = kept
    ghb.stress_period_data.set_data(spd)

    # Perturb the solver configuration too (config tier).
    variant.ims.outer_dvclose = 0.5

    # Door 1: model.diff(other) -- this model is the reference.
    diff = reference.diff(variant)
    summary = diff.summary()
    ghb_row = summary[
        (summary["package"] == "ghb") & (summary["model"] == "variant")
    ].iloc[0]
    assert ghb_row["cells_only_in_reference"] >= 1
    assert ghb_row["value_cells_changed"] >= 1
    assert not ghb_row["identical"]

    # The dropped cell shows up structurally as only-in-reference at period 0.
    only_ref = diff.packages.ghb.inputs.cells(per=0)
    dropped_cell = int(dropped_cellid[-1])
    assert dropped_cell in set(
        only_ref.loc[only_ref["membership"] == "only_in_reference", "cell"]
    )

    # Config tier detects the solver change, and model.config works on a real model.
    config_diffs = diff.config.settings()
    assert (config_diffs["setting"] == "outer_dvclose").any()
    assert not bool(diff.config.summary().iloc[0]["identical"])
    assert "ims" in reference.config.sections
    assert "Configuration differences" in diff.report()

    # Connection tier works on real LAK/SFR networks: unperturbed here, so the
    # lake connections and stream reaches must read as identical.
    assert bool(diff.packages.lak.inputs.summary().iloc[0]["identical"])
    assert diff.packages.lak.inputs.connections().empty
    assert bool(diff.packages.sfr.inputs.summary().iloc[0]["identical"])
    assert diff.packages.sfr.inputs.reaches().empty

    # Door 2: ModelGroup(...).diff() gives the same reference-star result.
    group_summary = (
        ModelGroup({"reference": reference, "variant": variant}, reference="reference")
        .diff()
        .summary()
    )
    group_row = group_summary[
        (group_summary["package"] == "ghb") & (group_summary["model"] == "variant")
    ].iloc[0]
    assert bool(group_row["identical"]) == bool(ghb_row["identical"])
    assert group_row["value_cells_changed"] == ghb_row["value_cells_changed"]

    # The list form of the door works too.
    assert not reference.diff([variant]).summary().empty
