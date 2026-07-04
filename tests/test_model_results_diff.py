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
    summary = diff.results.heads.summary().iloc[0]
    assert summary["max_abs_diff"] == 0.0
    assert bool(summary["within_tolerance"])


def test_head_difference_flagged_with_location():
    frame = _heads(
        "variant",
        [((0, 5), 0, 10, 5.0, 5.0), ((0, 5), 1, 42, 5.8, 5.0)],  # 0.8 diff at cell 42
    )
    budgets = {"ref": _budget([1.0], TOTAL_IN=[10.0]), "variant": _budget([1.0], TOTAL_IN=[10.0])}
    diff = _group(heads_frame=frame, budgets=budgets).diff()
    summary = diff.results.heads.summary().iloc[0]
    assert summary["max_abs_diff"] == pytest.approx(0.8)
    assert summary["argmax_cell"] == 42
    assert summary["argmax_layer"] == 1
    assert summary["argmax_kstpkper"] == (0, 5)
    assert not bool(summary["within_tolerance"])
    assert summary["rmse"] == pytest.approx(np.sqrt((0.0**2 + 0.8**2) / 2))


def test_head_tolerance_is_overridable():
    frame = _heads("variant", [((0, 0), 0, 10, 100.05, 100.0)])  # 0.05 diff, ref 100
    budgets = {"ref": _budget([1.0], TOTAL_IN=[10.0]), "variant": _budget([1.0], TOTAL_IN=[10.0])}
    heads = _group(heads_frame=frame, budgets=budgets).diff().results.heads
    # default atol+rtol (1e-3): 0.05 <= 1e-3 + 1e-3*100 = 0.101 -> within
    assert bool(heads.summary().iloc[0]["within_tolerance"])
    # tighten: atol=0.01, rtol=0 -> 0.05 > 0.01 -> breached
    assert not bool(heads.summary(atol=0.01, rtol=0.0).iloc[0]["within_tolerance"])


# --- budget ------------------------------------------------------------------
def test_identical_budget_terms_within_tolerance():
    frame = _heads("variant", [((0, 0), 0, 10, 5.0, 5.0)])
    budgets = {
        "ref": _budget([1.0, 2.0], RCH_IN=[100.0, 100.0], GHB_OUT=[50.0, 50.0]),
        "variant": _budget([1.0, 2.0], RCH_IN=[100.0, 100.0], GHB_OUT=[50.0, 50.0]),
    }
    summary = _group(heads_frame=frame, budgets=budgets).diff().results.budget.summary()
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
        .results.budget.summary()
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
    data = _group(heads_frame=frame, budgets=budgets).diff().results.budget.get()
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
        diff.results.budget.summary(model_name="ref")
    with pytest.raises(KeyError):
        diff.results.heads.summary(model_name="ghost")


@pytest.mark.slow
def test_results_diff_end_to_end_identical_on_canonical(canonical_run):
    """Real reader path: a run vs a fresh reload of its own outputs must read as
    identical within tolerance (the regression 'faithful-copy' case)."""

    import myflopy as mf

    run = canonical_run
    reloaded = mf.load_mf6_run(run.workspace)
    diff = ModelGroup({"a": run, "b": reloaded}, reference="a").diff()

    heads = diff.results.heads.summary().iloc[0]
    assert heads["max_abs_diff"] == 0.0
    assert bool(heads["within_tolerance"])

    budget = diff.results.budget.summary()
    assert not budget.empty
    assert bool(budget["within_tolerance"].all())

    report = diff.report(results=True)
    assert "Results differences" in report
