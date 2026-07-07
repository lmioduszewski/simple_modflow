"""Tests for the ModelDiff results tier (5a: heads + overall budget).

Heads and budget readers are faked so the diff/stats/tolerance/report logic runs
on controlled data without a real MF6 run.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from myflopy.project.model_config import ModelConfig
from myflopy.project.model_group import ModelGroup
from myflopy.project.model_results_diff import (
    CellBudgetResultDiff,
    MvrResultDiff,
    StageResultDiff,
)


class _FakeHeads:
    """Stands in for ModelGroup.hds; .compare() returns a fixed frame."""

    def __init__(self, frame: pd.DataFrame):
        self._frame = frame

    def compare(self, **kwargs) -> pd.DataFrame:
        return self._frame.copy()


class _FakeModel:
    def __init__(self, name, budget: pd.DataFrame | None = None):
        self.name = name
        self.package_names = []
        self.config = ModelConfig(pd.DataFrame(columns=["section", "setting", "value"]))
        self._budget = budget

    def budget_incremental(self):
        return self._budget.copy()


def _heads(model_name, rows):
    """rows: (kstpkper, layer, cell, elev, reference_elev)."""
    frame = pd.DataFrame(
        rows, columns=["kstpkper", "layer", "cell", "elev", "reference_elev"]
    )
    frame["model"] = model_name
    frame["reference_model"] = "ref"
    frame["diff"] = frame["elev"].astype(float) - frame["reference_elev"].astype(float)
    return frame[
        ["model", "reference_model", "kstpkper", "layer", "cell",
         "elev", "reference_elev", "diff"]
    ]


def _budget(totims, **terms):
    return pd.DataFrame({"totim": totims, **terms})


def _group(*, heads_frame, budgets, reference="ref"):
    models = {name: _FakeModel(name, budget) for name, budget in budgets.items()}
    group = ModelGroup(models, reference=reference)
    group.hds = _FakeHeads(heads_frame)
    return group


# --- heads -------------------------------------------------------------------
def test_identical_heads_are_within_tolerance():
    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0), ((0, 0), 0, 11, 6.0, 6.0)])
    budgets = {"ref": _budget([1.0], TOTAL_IN=[10.0]), "variant": _budget([1.0], TOTAL_IN=[10.0])}
    diff = _group(heads_frame=frame, budgets=budgets).diff()
    summary = diff.hds.summary().iloc[0]
    assert summary["max_abs_diff"] == 0.0
    assert bool(summary["within_tolerance"])


def test_head_difference_flagged_with_location():
    frame = _heads(
        "variant",
        [((0, 5), 0, 10, 5.0, 5.0), ((0, 5), 1, 42, 5.8, 5.0)],  # 0.8 diff at cell 42
    )
    budgets = {"ref": _budget([1.0], TOTAL_IN=[10.0]), "variant": _budget([1.0], TOTAL_IN=[10.0])}
    diff = _group(heads_frame=frame, budgets=budgets).diff()
    summary = diff.hds.summary().iloc[0]
    assert summary["max_abs_diff"] == pytest.approx(0.8)
    assert summary["argmax_cell"] == 42
    assert summary["argmax_layer"] == 1
    assert summary["argmax_kstpkper"] == (0, 5)
    assert not bool(summary["within_tolerance"])
    assert summary["rmse"] == pytest.approx(np.sqrt((0.0**2 + 0.8**2) / 2))


def test_head_tolerance_is_overridable():
    frame = _heads("variant", [((0, 0), 0, 10, 100.05, 100.0)])  # 0.05 diff, ref 100
    budgets = {"ref": _budget([1.0], TOTAL_IN=[10.0]), "variant": _budget([1.0], TOTAL_IN=[10.0])}
    heads = _group(heads_frame=frame, budgets=budgets).diff().hds
    # default atol+rtol (1e-3): 0.05 <= 1e-3 + 1e-3*100 = 0.101 -> within
    assert bool(heads.summary().iloc[0]["within_tolerance"])
    # tighten: atol=0.01, rtol=0 -> 0.05 > 0.01 -> breached
    assert not bool(heads.summary(atol=0.01, rtol=0.0).iloc[0]["within_tolerance"])


def test_heads_map_forwards_model_name():
    """diff.hds.map('X') -> the group heads compare_map (Δhead map)."""

    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0)])
    budgets = {"ref": _budget([1.0], TOTAL_IN=[10.0]), "variant": _budget([1.0], TOTAL_IN=[10.0])}
    group = _group(heads_frame=frame, budgets=budgets)
    diff = group.diff()

    captured = {}
    group.hds.compare_map = lambda **kwargs: captured.update(kwargs) or "FIG"
    assert diff.hds.map("variant", per=3, layer=1) == "FIG"
    assert captured["model_name"] == "variant"
    assert captured["per"] == 3 and captured["layer"] == 1


# --- budget ------------------------------------------------------------------
def test_identical_budget_terms_within_tolerance():
    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0)])
    budgets = {
        "ref": _budget([1.0, 2.0], RCH_IN=[100.0, 100.0], GHB_OUT=[50.0, 50.0]),
        "variant": _budget([1.0, 2.0], RCH_IN=[100.0, 100.0], GHB_OUT=[50.0, 50.0]),
    }
    summary = _group(heads_frame=frame, budgets=budgets).diff().bud.summary()
    assert bool(summary["within_tolerance"].all())
    assert set(summary["term"]) == {"RCH_IN", "GHB_OUT"}


def test_budget_term_difference_detected():
    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0)])
    budgets = {
        "ref": _budget([1.0, 2.0], RCH_IN=[100.0, 100.0], GHB_OUT=[50.0, 50.0]),
        "variant": _budget([1.0, 2.0], RCH_IN=[100.0, 100.0], GHB_OUT=[60.0, 70.0]),
    }
    summary = (
        _group(heads_frame=frame, budgets=budgets)
        .diff()
        .bud.summary()
        .set_index("term")
    )
    assert bool(summary.loc["RCH_IN", "within_tolerance"])
    assert not bool(summary.loc["GHB_OUT", "within_tolerance"])
    # totals: ref 100, model 130 -> +30, +30%
    assert summary.loc["GHB_OUT", "diff_total"] == pytest.approx(30.0)
    assert summary.loc["GHB_OUT", "pct_change"] == pytest.approx(30.0)


def test_budget_aligns_on_totim_intersection():
    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0)])
    budgets = {
        "ref": _budget([1.0, 2.0, 3.0], RCH_IN=[100.0, 100.0, 100.0]),
        "variant": _budget([2.0, 3.0], RCH_IN=[100.0, 100.0]),  # missing totim 1.0
    }
    data = _group(heads_frame=frame, budgets=budgets).diff().bud.get()
    assert set(data["totim"]) == {2.0, 3.0}  # only the shared timesteps


# --- report + errors ---------------------------------------------------------
def test_results_appear_in_report_only_when_requested():
    frame = _heads("variant", [((0, 5), 1, 42, 5.8, 5.0)])
    budgets = {"ref": _budget([1.0], RCH_IN=[100.0]), "variant": _budget([1.0], RCH_IN=[130.0])}
    diff = _group(heads_frame=frame, budgets=budgets).diff()
    assert "Results differences" not in diff.report()          # opt-in
    report = diff.report(results=True)
    assert "Results differences" in report
    assert "heads" in report and "budget" in report


def test_results_self_reference_errors():
    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0)])
    budgets = {"ref": _budget([1.0], RCH_IN=[100.0]), "variant": _budget([1.0], RCH_IN=[100.0])}
    diff = _group(heads_frame=frame, budgets=budgets).diff()
    with pytest.raises(ValueError):
        diff.bud.summary(model_name="ref")
    with pytest.raises(KeyError):
        diff.hds.summary(model_name="ghost")


# --- per-package cell budget + UZF (5b) --------------------------------------
class _FakeCellAccessor:
    def __init__(self, package_name, value_name, frame):
        self.package_name = package_name
        self.value_name = value_name
        self._frame = frame

    def compare(self, **kwargs):
        return self._frame.copy()


def _cell_frame(rows):
    """rows: (cell, q, reference_q)."""
    frame = pd.DataFrame(rows, columns=["cell", "q", "reference_q"])
    frame["model"] = "variant"
    frame["reference_model"] = "ref"
    frame["per"] = 0
    frame["layer"] = 0
    frame["q_diff"] = frame["q"].astype(float) - frame["reference_q"].astype(float)
    return frame


def _bare_diff(names=("ref", "variant")):
    return ModelGroup({n: _FakeModel(n) for n in names}, reference="ref").diff()


def test_cell_budget_summary_stats_and_tolerance():
    frame = _cell_frame([(10, -5.0, -5.0), (11, -8.0, -5.0)])  # 3.0 diff at cell 11
    cbd = CellBudgetResultDiff(_bare_diff(), _FakeCellAccessor("ghb", "q", frame))
    row = cbd.summary().iloc[0]
    assert row["package"] == "ghb" and row["value"] == "q"
    assert row["max_abs_diff"] == pytest.approx(3.0)
    assert row["argmax_cell"] == 11
    assert not bool(row["within_tolerance"])


def test_cell_budget_identical_within_tolerance():
    frame = _cell_frame([(10, -5.0, -5.0), (11, -8.0, -8.0)])
    cbd = CellBudgetResultDiff(_bare_diff(), _FakeCellAccessor("drn", "q", frame))
    assert bool(cbd.summary().iloc[0]["within_tolerance"])


def test_cell_budget_map_forwards_model_name():
    """diff.packages.sfr.results.q.map('X') / .packages.ghb.map('X') must forward the
    model to the accessor's compare_map (the SFR/LAK results Δ map)."""

    accessor = _FakeCellAccessor("ghb", "q", _cell_frame([(10, -5.0, -5.0)]))
    captured = {}
    accessor.compare_map = lambda **kwargs: captured.update(kwargs) or "FIG"
    cbd = CellBudgetResultDiff(_bare_diff(), accessor)

    assert cbd.map("variant") == "FIG"          # positional model name accepted
    assert captured["model_name"] == "variant"


def test_cell_budget_map_requires_accessor_support():
    accessor = _FakeCellAccessor("ghb", "q", _cell_frame([(10, -5.0, -5.0)]))  # no compare_map
    cbd = CellBudgetResultDiff(_bare_diff(), accessor)
    with pytest.raises(AttributeError):
        cbd.map("variant")


def test_results_packages_namespace_resolves():
    diff = _bare_diff()
    ghb = diff.packages.ghb.results.q
    assert isinstance(ghb, CellBudgetResultDiff)
    assert ghb.package_name == "ghb"
    with pytest.raises(AttributeError):
        _ = diff.packages.definitely_not_a_package


def test_results_uzf_namespace_resolves():
    diff = _bare_diff()
    assert diff.packages.uzf.results.gwrch.value_name == "gwrch"
    assert diff.packages.uzf.results.sat.value_name == "sat"
    assert diff.packages.uzf.results.gwrch.package_name == "uzf"


# --- LAK/SFR stage + MVR (5c) ------------------------------------------------
class _FakeStageAccessor:
    def __init__(self, frame, entity):
        self._frame = frame
        self._entity = entity

    def compare(self, **kwargs):  # ignores per/entity filters for the test
        return self._frame.copy()


def _stage_frame(entity, rows):
    """rows: (per, feature, stage, reference_stage)."""
    frame = pd.DataFrame(rows, columns=["per", entity, "stage", "reference_stage"])
    frame["model"] = "variant"
    frame["reference_model"] = "ref"
    frame["stage_diff"] = frame["stage"].astype(float) - frame["reference_stage"].astype(float)
    return frame[
        ["model", "reference_model", "per", entity, "stage", "reference_stage", "stage_diff"]
    ]


def test_lake_stage_diff_flags_overtop_location():
    # lake 1 rises 0.46 in period 8 (facility-stage style)
    frame = _stage_frame("lake", [(0, 0, 10.0, 10.0), (8, 1, 388.21, 387.75)])
    sd = StageResultDiff(_bare_diff(), _FakeStageAccessor(frame, "lake"), entity="lake")
    row = sd.summary().iloc[0]
    assert row["max_abs_diff"] == pytest.approx(0.46)
    assert row["argmax_lake"] == 1
    assert row["argmax_per"] == 8
    assert not bool(row["within_tolerance"])


def test_sfr_reach_stage_diff_within_tolerance():
    frame = _stage_frame("reach", [(0, 5, 3.0, 3.0), (1, 5, 3.0, 3.0)])
    sd = StageResultDiff(_bare_diff(), _FakeStageAccessor(frame, "reach"), entity="reach")
    summary = sd.summary()
    assert "argmax_reach" in summary.columns
    assert bool(summary.iloc[0]["within_tolerance"])


def test_mvr_summary_empty_and_columns_without_mover_output():
    diff = _bare_diff()
    summary = diff.packages.mvr.results.summary()
    assert summary.empty
    assert list(summary.columns) == MvrResultDiff._COLUMNS


def test_results_lak_sfr_mvr_namespaces_resolve():
    diff = _bare_diff()
    assert isinstance(diff.packages.lak.results.stage, StageResultDiff)
    assert diff.packages.lak.results.stage._entity == "lake"
    assert isinstance(diff.packages.sfr.results.stage, StageResultDiff)
    assert diff.packages.sfr.results.stage._entity == "reach"
    assert isinstance(diff.packages.lak.results.q, CellBudgetResultDiff)
    assert diff.packages.sfr.results.q.package_name == "sfr"
    assert isinstance(diff.packages.mvr.results, MvrResultDiff)


@pytest.mark.slow
def test_results_diff_end_to_end_identical_on_canonical(canonical_run):
    """Real reader path: a run vs a fresh reload of its own outputs must read as
    identical within tolerance (the regression 'faithful-copy' case)."""

    import myflopy as mf

    run = canonical_run
    reloaded = mf.load_mf6_run(run.workspace)
    diff = ModelGroup({"a": run, "b": reloaded}, reference="a").diff()

    heads = diff.hds.summary().iloc[0]
    assert heads["max_abs_diff"] == 0.0
    assert bool(heads["within_tolerance"])

    budget = diff.bud.summary()
    assert not budget.empty
    assert bool(budget["within_tolerance"].all())

    # per-package cell budget + UZF on real outputs (identical -> within tolerance)
    ghb_cells = diff.packages.ghb.results.q.summary()
    assert not ghb_cells.empty
    assert bool(ghb_cells["within_tolerance"].all())
    uzf = diff.packages.uzf.results.gwrch.summary()
    assert bool(uzf["within_tolerance"].all())

    # LAK/SFR stage + SFR flow on real outputs (identical -> within tolerance)
    lak_stage = diff.packages.lak.results.stage.summary()
    assert not lak_stage.empty
    assert bool(lak_stage["within_tolerance"].all())
    sfr_stage = diff.packages.sfr.results.stage.summary()
    assert not sfr_stage.empty
    assert bool(sfr_stage["within_tolerance"].all())
    assert bool(diff.packages.sfr.results.q.summary()["within_tolerance"].all())

    # MVR (canonical routes lake -> stream via the mover); if terms are present,
    # a run vs a reload of itself must be within tolerance.
    mvr = diff.packages.mvr.results.summary()
    if not mvr.empty:
        assert bool(mvr["within_tolerance"].all())

    # spatial diff maps render on real outputs (a==b so flat, but must build)
    assert diff.hds.map("b") is not None
    assert diff.packages.sfr.results.q.map("b") is not None
    assert diff.packages.ghb.results.q.map("b") is not None

    report = diff.report(results=True)
    assert "Results differences" in report
