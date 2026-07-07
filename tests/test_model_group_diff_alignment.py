"""Regression tests for ModelGroup diff alignment across models whose runs
differ in time discretization and carry floating-point geometry noise.

These exercise the *real* ``compare()`` merge logic (the results-tier fast tests
in ``test_model_results_diff.py`` fake ``.compare()``, so they never touched the
merge key). Three bugs are covered, all found diffing real F9 model variants:

1. ``LoadedMf6Run._get_budget_reader`` must open the ``.cbc`` as double
   precision -- FloPy's auto-detect misreads MF6's double budget file as single
   and, on Windows, raises ``OSError`` (Errno 22) that escapes its own guard.
2. Cell-budget ``compare()`` must join on stable integer identity only. Float
   geometry (``rlen``/``distance_*``) carries sub-unit per-model noise and the
   ``(kstp, kper)`` tuple varies with a model's time steps; either in the merge
   key silently empties the comparison.
3. Heads ``compare()`` must align on the stress *period*, not the ``(kstp,
   kper)`` tuple, so models with different ``nstp`` still line up.
"""

from __future__ import annotations

import types

import pandas as pd
import pytest

from myflopy.project.model_group import (
    GroupCellPackageResults,
    GroupHeads,
    GroupSfrBudgetResults,
    ModelGroup,
    _stable_compare_keys,
)


class _FakeModel:
    def __init__(self, name):
        self.name = name
        self.package_names = []


def _group(reference="ref", names=("ref", "variant")):
    return ModelGroup({n: _FakeModel(n) for n in names}, reference=reference)


# --- helper ------------------------------------------------------------------
def test_stable_compare_keys_drops_floats_values_and_time_metadata():
    frame = pd.DataFrame(
        {
            "model": ["ref"],
            "package": ["sfr"],
            "kstpkper": [(13, 0)],
            "per": [0],
            "reach": pd.array([475], dtype="Int64"),
            "layer": [0],
            "cell": [105],
            "node": [105],
            "node2": [475],
            "rlen": [16.000493],        # float geometry -> must be excluded
            "distance_mid": [16564.63],  # float geometry -> must be excluded
            "q": [-166.29],             # value -> excluded
            "q_per_length": [-10.4],    # value -> excluded
        }
    )
    keys = _stable_compare_keys(frame, ["q", "q_per_length"])
    assert keys == ["package", "per", "reach", "layer", "cell", "node", "node2"]
    # explicitly: no float, no value, no model, no kstpkper
    for dropped in ("rlen", "distance_mid", "q", "q_per_length", "model", "kstpkper"):
        assert dropped not in keys


# --- bug 2: cell-budget merge ------------------------------------------------
def test_cell_budget_compare_aligns_despite_different_kstp():
    """GHB budget: two models save period 0 at kstp 13 vs kstp 9. The old key
    merged on the ``(kstp, kper)`` tuple and returned nothing; aligning on the
    period recovers the shared cell."""

    acc = GroupCellPackageResults(_group(), "ghb", budget_text="GHB", value_name="q")
    ref = pd.DataFrame(
        {
            "model": ["ref"], "package": ["ghb"], "kstpkper": [(13, 0)],
            "per": [0], "layer": [0], "cell": [100], "node": [100],
            "node2": [0], "q": [-5.0],
        }
    )
    variant = ref.copy()
    variant["model"] = "variant"
    variant["kstpkper"] = [(9, 0)]   # different time-step numbering
    variant["q"] = [-7.0]
    combined = pd.concat([ref, variant], ignore_index=True)
    acc.get = lambda **kwargs: combined.copy()

    comp = acc.compare()
    assert len(comp) == 1
    row = comp.iloc[0]
    assert row["cell"] == 100
    assert row["reference_q"] == pytest.approx(-5.0)
    assert row["q_diff"] == pytest.approx(-2.0)


def test_sfr_flow_compare_survives_float_geometry_noise():
    """SFR exchange: a shared reach whose ``rlen``/``distance_mid`` differ by
    sub-millimeter FP noise between models must NOT be dropped from the diff."""

    acc = GroupSfrBudgetResults(_group())
    ref = pd.DataFrame(
        {
            "model": ["ref"], "package": ["sfr"], "kstpkper": [(13, 0)],
            "per": [0], "reach": pd.array([475], dtype="Int64"), "layer": [0],
            "cell": [105], "rlen": [16.000493], "distance_start": [0.0],
            "distance_mid": [16564.632017], "distance_end": [16580.0],
            "q": [-166.294714], "q_per_length": [-166.294714 / 16.000493],
            "node": [105], "node2": [475],
        }
    )
    variant = ref.copy()
    variant["model"] = "variant"
    variant["kstpkper"] = [(9, 0)]
    variant["rlen"] = [16.000000]            # FP noise vs 16.000493
    variant["distance_mid"] = [16564.632264]  # FP noise vs ...017
    variant["q"] = [-166.331694]
    variant["q_per_length"] = [-166.331694 / 16.0]
    combined = pd.concat([ref, variant], ignore_index=True)
    acc.get = lambda **kwargs: combined.copy()

    comp = acc.compare()
    assert len(comp) == 1  # the shared reach survives the FP-noisy geometry
    row = comp.iloc[0]
    assert row["reach"] == 475
    assert row["q_diff"] == pytest.approx(-166.331694 - (-166.294714))


def test_cell_budget_compare_reduces_multi_step_period_no_cartesian():
    """A model saving several steps per period must not self-join across steps.
    The compare collapses to the period-end flux (regression: the canonical
    model saves multiple timesteps per period, which previously produced a
    2x2 cartesian and spurious non-zero diffs on identical outputs)."""

    acc = GroupCellPackageResults(_group(), "ghb", budget_text="GHB", value_name="q")

    def two_steps(model, q_mid, q_end):
        return pd.DataFrame(
            {
                "model": [model, model], "package": ["ghb", "ghb"],
                "kstpkper": [(1, 0), (5, 0)], "per": [0, 0], "layer": [0, 0],
                "cell": [100, 100], "node": [100, 100], "node2": [0, 0],
                "q": [q_mid, q_end],
            }
        )

    combined = pd.concat(
        [two_steps("ref", -3.0, -5.0), two_steps("variant", -3.0, -6.0)],
        ignore_index=True,
    )
    acc.get = lambda **kwargs: combined.copy()

    comp = acc.compare()
    assert len(comp) == 1                          # period-end only, not a 2x2 join
    assert comp.iloc[0]["q"] == pytest.approx(-6.0)          # variant period-end
    assert comp.iloc[0]["reference_q"] == pytest.approx(-5.0)  # ref period-end
    assert comp.iloc[0]["q_diff"] == pytest.approx(-1.0)


def test_cell_budget_compare_matches_prior_when_discretization_agrees():
    """When both models share time-step numbering the merge is unaffected --
    identical value -> zero diff on the shared cell."""

    acc = GroupCellPackageResults(_group(), "drn", budget_text="DRN", value_name="q")
    ref = pd.DataFrame(
        {
            "model": ["ref"], "package": ["drn"], "kstpkper": [(1, 0)],
            "per": [0], "layer": [0], "cell": [7], "node": [7], "node2": [0],
            "q": [-3.0],
        }
    )
    variant = ref.copy()
    variant["model"] = "variant"
    combined = pd.concat([ref, variant], ignore_index=True)
    acc.get = lambda **kwargs: combined.copy()

    comp = acc.compare()
    assert len(comp) == 1
    assert comp.iloc[0]["q_diff"] == pytest.approx(0.0)


# --- bug 3: heads period alignment -------------------------------------------
def test_heads_compare_aligns_periods_across_different_kstp():
    heads = GroupHeads(_group())
    data = pd.DataFrame(
        {
            "model": ["ref", "ref", "variant", "variant"],
            "kstpkper": [(13, 0), (13, 1), (9, 0), (9, 1)],  # same periods, diff kstp
            "layer": [0, 0, 0, 0],
            "cell": [10, 10, 10, 10],
            "elev": [5.0, 6.0, 5.5, 6.5],
        }
    )
    heads.get = lambda **kwargs: data.copy()

    comp = heads.compare().sort_values("per").reset_index(drop=True)
    assert len(comp) == 2                       # both periods aligned
    assert set(comp["per"]) == {0, 1}
    assert comp.loc[0, "diff"] == pytest.approx(0.5)  # 5.5 - 5.0
    assert comp.loc[1, "diff"] == pytest.approx(0.5)  # 6.5 - 6.0


def test_heads_compare_reduces_to_period_end_when_multiple_steps_saved():
    heads = GroupHeads(_group())
    data = pd.DataFrame(
        {
            "model": ["ref", "ref", "variant", "variant"],
            "kstpkper": [(1, 0), (5, 0), (1, 0), (5, 0)],  # two saves in period 0
            "layer": [0, 0, 0, 0],
            "cell": [10, 10, 10, 10],
            "elev": [1.0, 9.0, 1.0, 9.5],
        }
    )
    heads.get = lambda **kwargs: data.copy()

    comp = heads.compare()
    assert len(comp) == 1                        # collapsed to the period-end step
    assert comp.iloc[0]["kstpkper"] == (5, 0)    # kept the largest-kstp record
    assert comp.iloc[0]["diff"] == pytest.approx(0.5)  # 9.5 - 9.0


# --- bug 1: CBC precision ----------------------------------------------------
def test_loaded_run_opens_cbc_as_double(monkeypatch, tmp_path):
    import flopy

    from myflopy.project.run_model import LoadedMf6Run

    captured = {}

    def fake_cellbudgetfile(path, **kwargs):
        captured["path"] = path
        captured["kwargs"] = kwargs
        return "READER"

    monkeypatch.setattr(flopy.utils, "CellBudgetFile", fake_cellbudgetfile)

    stub = types.SimpleNamespace(_bud=None, workspace=tmp_path, name="mymodel")
    reader = LoadedMf6Run._get_budget_reader(stub)

    assert reader == "READER"
    assert captured["kwargs"].get("precision") == "double"
    assert str(captured["path"]).endswith("mymodel.cbc")
    # cached on the instance
    assert stub._bud == "READER"
