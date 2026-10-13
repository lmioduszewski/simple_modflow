"""What the review says it read, and what a write-back and the forward run read.

The leftovers of the PEST-preparation run lanes:

* ``sensitivity()`` and ``parameters_at_bounds()`` said nothing when PESTPP-IES had
  dropped most realizations (ledger 575);
* a write-back filled templates with a parameter's value whatever its ``scale`` and
  ``offset``;
* the pilot-point forward-run step read its values with pandas' default float parser;
* DRN and RIV flow sets took no set ``weight=``.

The model is the toy of ``tests/test_pest_handed_targets.py``; PESTPP-IES runs for real
where the test is about what it does with a failed realization.
"""

from __future__ import annotations

import io
import json
import shutil
import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest import forward_run
from myflopy.modflow.mf6.pest.calibrated import _complete_values
from myflopy.modflow.mf6.pest.ies import IesResults
from myflopy.modflow.mf6.pest.pilot_points import PILOT_POINT_MANIFEST, idw_factors
from tests.test_pest_handed_targets import DRN_CELLS, run_toy

pyemu = pytest.importorskip("pyemu")

needs_pestpp = pytest.mark.skipif(shutil.which("pestpp-ies") is None,
                                  reason="needs pestpp-ies on PATH")

HEADS = pd.DataFrame({"time": [0, 1, 2], "W2": [91.0, 92.0, 93.0],
                      "W1": [93.0, 93.5, 94.0], "W10": [95.0, 96.0, 97.0]})
LOCATIONS = pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0})


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return run_toy(tmp_path_factory.mktemp("review_leftovers"))


# -- the review says how many realizations its statistics rest on -----------------------


def _review_stub(prior_rows: int, posterior_rows: int, prior_ran: int):
    parameter_data = pd.DataFrame(
        {"parnme": ["p"], "pargp": ["g"], "partrans": ["log"], "parlbnd": [0.5],
         "parubnd": [2.0]}).set_index("parnme")
    rng = np.random.default_rng(0)
    prior = pd.DataFrame({"p": 10 ** rng.normal(0, 0.15, prior_rows)})
    posterior = pd.DataFrame({"p": 10 ** rng.normal(0, 0.05, posterior_rows)})

    class _Stub(IesResults):
        def __init__(self):
            self.pst = SimpleNamespace(parameter_data=parameter_data)

        prior_parameters = property(lambda self: SimpleNamespace(_df=prior))
        posterior_parameters = property(lambda self: SimpleNamespace(_df=posterior))
        prior = property(lambda self: SimpleNamespace(_df=prior.iloc[:prior_ran]))

    return _Stub()


def test_the_review_reports_and_warns_when_most_realizations_were_dropped():
    few = _review_stub(prior_rows=12, posterior_rows=5, prior_ran=9)
    with pytest.warns(UserWarning, match=r"sensitivity\(\): the posterior ensemble holds 5 "
                                         r"of the 12 realizations"):
        table = few.sensitivity()
    assert table.loc["g", ["n_reals", "n_reals_drawn"]].tolist() == [5, 12]
    with pytest.warns(UserWarning, match=r"parameters_at_bounds\(\): the posterior"):
        bounds = few.parameters_at_bounds()
    assert bounds.loc["g", ["n_reals", "n_reals_drawn"]].tolist() == [5, 12]

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        enough = _review_stub(prior_rows=12, posterior_rows=6, prior_ran=12)
        assert enough.sensitivity().loc["g", "n_reals"] == 6


def test_a_prior_table_rests_on_every_realization_drawn():
    """The prior PARAMETER ensemble keeps a realization whose run failed, so a prior
    table averages all of them: ``n_reals`` is every one drawn, and nothing warns
    however few ran -- it read the 3 that ran here and warned the means rested on 3
    (review of ledger 575)."""

    few = _review_stub(prior_rows=12, posterior_rows=5, prior_ran=3)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        table = few.parameters_at_bounds(which="prior")
    assert table.loc["g", ["n_reals", "n_reals_drawn"]].tolist() == [12, 12]
    assert table.loc["g", "n_parameters"] == 1


@needs_pestpp
@pytest.mark.slow
def test_pestpp_drops_failed_realizations_from_the_ensembles_the_review_counts(
        toy, tmp_path):
    """Real PESTPP-IES, its forward run failing every realization whose K multiplier is
    below 1: the prior parameter ensemble keeps every realization drawn, its observation
    ensemble only the runs that finished, and the posterior only those PESTPP-IES kept --
    which is what ``n_reals`` and ``n_reals_drawn`` count."""

    cal = toy.pest("drop", workspace=tmp_path / "drop", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    cal.parameterize("chd", style="constant", bounds=(-1.0, 1.0))
    cal.observe(mf.HeadTargets(locations=LOCATIONS, values=HEADS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("drop.pst", noptmax=0)
    template = cal.template_workspace
    script = template / "forward_run.py"
    script.write_text(
        "import sys as _sys\n"
        "if float(open('mult/k_inst0_constant.csv').read().split()[0]) < 1.0:\n"
        "    _sys.exit(1)\n" + script.read_text())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        review = cal.run_ies(reals=12, iterations=1, ies_lambda_mults=1.0,
                             lambda_scale_fac=1.0)
    drawn = len(pd.read_csv(template / "drop.0.par.csv", index_col=0))
    ran = len(pd.read_csv(template / "drop.0.obs.csv", index_col=0))
    kept = len(pd.read_csv(template / "drop.1.par.csv", index_col=0))
    assert drawn == 12 and ran < drawn and kept <= ran
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = review.sensitivity()
        prior = review.parameters_at_bounds(which="prior")
    assert set(table["n_reals"]) == {kept} and set(table["n_reals_drawn"]) == {drawn}
    # The prior table averages every realization drawn, a failed run's included.
    assert set(prior["n_reals"]) == {drawn}


# -- a write-back refuses a scale or offset ----------------------------------------------


def test_a_write_back_refuses_a_parameter_with_a_scale_or_offset():
    data = pd.DataFrame({"parnme": ["a", "b"], "partrans": ["log", "none"],
                         "parval1": [1.0, 0.0], "scale": [1.0, 1.0],
                         "offset": [0.0, 0.0]}).set_index("parnme", drop=False)
    pst = SimpleNamespace(parameter_data=data)
    values = pd.Series({"a": 1.5, "b": 0.25})
    assert _complete_values(pst, values) == {"a": 1.5, "b": 0.25}
    data.loc["b", "offset"] = 100.0
    with pytest.raises(ValueError, match=r"\['b'\] carry a scale or offset -- 'b': scale 1, "
                                         r"offset 100"):
        _complete_values(pst, values)
    data.loc["b", "offset"] = 0.0
    data.loc["a", "scale"] = 2.0
    with pytest.raises(ValueError, match=r"'a': scale 2, offset 0"):
        _complete_values(pst, values)


# -- the pilot-point step reads its values exactly ---------------------------------------


def test_the_pilot_point_step_reads_each_value_as_written(tmp_path, monkeypatch):
    """``0.12076665792328141`` reads one ulp low through pandas' default parser; on a base
    chosen at a rounding edge of the file's 11 digits that ulp moved the written value."""

    text = "0.12076665792328141"
    exact = float(text)
    default = float(pd.read_csv(io.StringIO(f"value\n{text}"))["value"].iloc[0])
    assert default != exact, "pandas' parser now reads this exactly: pick another value"
    base = None
    for target in np.linspace(1.0, 9.0, 4001):
        edge = float(f"{target:.10E}") + 0.5e-10
        start = edge / exact
        for step in range(-40, 41):
            candidate = start + step * np.spacing(start)
            if f"{candidate * exact:.10E}" != f"{candidate * default:.10E}":
                base = candidate
                break
        if base is not None:
            break
    assert base is not None

    (tmp_path / "s_pp.csv").write_text(f"parnme,value\nsp0000,{text}\nsp0001,{text}\n")
    pd.DataFrame({"parnme": ["sp0000", "sp0001"], "x": [0.0, 1.0], "y": [0.0, 0.0],
                  "layer": 0, "zone": 0}).to_csv(tmp_path / "s_pp.points.csv", index=False)
    index, weight = idw_factors(np.array([0.0]), np.zeros(1), [0.0, 1.0], [0.0, 0.0])
    np.savez(tmp_path / "s_pp.factors.npz", cells=np.arange(1), x=np.zeros(1),
             y=np.zeros(1), index=index, weight=weight)
    np.savetxt(tmp_path / "a.txt.pp_base", np.array([[base]]), fmt="%.18e")
    np.savetxt(tmp_path / "a.txt", np.array([[base]]), fmt="%.18e")
    (tmp_path / PILOT_POINT_MANIFEST).write_text(json.dumps({
        "version": 1,
        "files": [{"array_file": "a.txt", "base_from_file": False,
                   "base_file": "a.txt.pp_base", "sets": ["s"]}],
        "sets": [{"name": "s", "array_files": ["a.txt"], "parameter_csv": "s_pp.csv",
                  "points_csv": "s_pp.points.csv", "factors": "s_pp.factors.npz",
                  "interpolation": "linear", "lower_limit": None, "upper_limit": None}],
    }))
    monkeypatch.chdir(tmp_path)
    forward_run.apply_pilotpoints_to_array(["a.txt"])
    assert (tmp_path / "a.txt").read_text().split()[0] == f"{base * exact:.10E}"


# -- a DRN or RIV set takes a set weight ---------------------------------------------------


def test_a_drn_set_takes_a_set_weight_below_a_zones_own(toy, tmp_path):
    """Row > location > set > 1.0 (ledger 504): the spec's ``weight`` fills a reading that
    has neither its own nor its zone's. RIV sets share the builder."""

    targets = mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["seep", "own"], "cells": [DRN_CELLS[:2], DRN_CELLS[2:]],
                                "weight": [None, 4.0]}),
        values=pd.DataFrame({"time": [0, 1], "seep": [-3.0, -3.0], "own": [-1.0, -1.0]}))
    cal = toy.pest("w", workspace=tmp_path / "w", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    cal.observe(mf.DrnFlowObservationSpec(targets=targets, weight=0.25))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("w.pst", noptmax=0)
    measured = pst.observation_data[pst.observation_data["weight"] > 0]
    # Two zones x two periods: `seep` at the set's 0.25, `own` at its zone's 4.0.
    assert sorted(measured["weight"]) == [0.25, 0.25, 4.0, 4.0]
    assert mf.RivFlowObservationSpec(targets=None, weight=0.5).weight == 0.5
    assert mf.DrnFlowObservationSpec(targets=None).weight is None
