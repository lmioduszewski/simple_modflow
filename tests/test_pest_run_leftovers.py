"""What a PESTPP-IES launch sets, draws and records.

The leftovers of the PEST-preparation run lanes:

* ``noise=True`` set nothing, so an earlier ``noise=False`` stayed on ``cal.pst`` and the
  next run went without noise (ledger 572) -- and a removed option is no noise either:
  PESTPP-IES 5.2.16 given none draws none, so ``noise=True`` writes ``ies_no_noise
  false`` (ledger 590);
* the correlated prior myflopy draws into ``ies_par_en`` came from numpy's global
  generator with no seed, so no run could be drawn again (ledger 573);
* a parallel run given ``master_dir=`` outside the default place was not found by
  ``model.pest_runs`` (ledger 508's "Not done"; ledger 574).

The model is the toy of ``tests/test_pest_handed_targets.py``; PESTPP-IES is stubbed,
since these are about what the launch writes. ``tests/test_pest_review_leftovers.py``
holds the review's and the forward run's.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.runs import METADATA_FILENAME, find_pest_runs
from tests.test_pest_handed_targets import run_toy

pyemu = pytest.importorskip("pyemu")

HEADS = pd.DataFrame({"time": [0, 1, 2], "W2": [91.0, 92.0, 93.0],
                      "W1": [93.0, 93.5, 94.0], "W10": [95.0, 96.0, 97.0]})
LOCATIONS = pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0})


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return run_toy(tmp_path_factory.mktemp("run_leftovers"))


def _calibration(toy, workspace: Path, name: str, *, correlated: bool):
    cal = toy.pest(name, workspace=workspace, start_datetime="2024-01-01")
    if correlated:
        cal.parameterize("k", style="grid", correlation=300.0, bounds=(0.5, 2.0))
    else:
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    cal.observe(mf.HeadTargets(locations=LOCATIONS, values=HEADS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}.pst")
    return cal


@pytest.fixture
def plain(toy, tmp_path):
    cal = _calibration(toy, tmp_path / "pest" / "plain", "plain", correlated=False)
    yield cal
    cal._release_pyemu_log()


@pytest.fixture
def correlated(toy, tmp_path):
    cal = _calibration(toy, tmp_path / "pest" / "corr", "corr", correlated=True)
    yield cal
    cal._release_pyemu_log()


def _stub_pestpp(monkeypatch, cal):
    """PESTPP-IES replaced by a recorder; the review is not read (nothing was run)."""

    import myflopy.modflow.mf6.pest.ies as ies

    commands = []
    monkeypatch.setattr(cal.pyemu.os_utils, "run",
                        lambda command, cwd=None: commands.append((command, cwd)))
    monkeypatch.setattr(ies, "IesResults", lambda workspace, **_: Path(workspace))
    return commands


def _executions(cal) -> list[dict]:
    return json.loads((cal.template_workspace / METADATA_FILENAME).read_text())["executions"]


# -- noise=True clears what noise=False set ----------------------------------------------


def test_noise_true_clears_an_earlier_noise_false(plain, monkeypatch):
    cal = plain
    _stub_pestpp(monkeypatch, cal)
    written = cal.template_workspace / "plain.pst"

    cal.run_ies(reals=3, iterations=1, noise=False)
    assert cal.pst.pestpp_options["ies_no_noise"] is True
    assert "ies_no_noise" in written.read_text()
    cal.run_ies(reals=3, iterations=1)
    # WRITTEN false, not left out: PESTPP-IES 5.2.16 given no noise option draws none
    # (ledger 590).
    assert cal.pst.pestpp_options["ies_no_noise"] is False
    assert _written_option(written, "ies_no_noise") == "false"
    cal.prior(reals=3, noise=False)
    assert _written_option(written, "ies_no_noise") == "true"
    cal.prior(reals=3)
    assert _written_option(written, "ies_no_noise") == "false"


def _written_option(pst_file, name):
    """The value a control file writes for PEST++ option ``name``, lower-cased."""

    lines = [line.split() for line in pst_file.read_text().splitlines()]
    values = [parts[1].lower() for parts in lines if len(parts) == 2 and parts[0] == name]
    assert len(values) == 1, values
    return values[0]


def test_ies_no_noise_by_hand_is_refused_naming_noise(plain, monkeypatch):
    _stub_pestpp(monkeypatch, plain)
    with pytest.raises(ValueError, match=r"run_ies\(\) sets ies_no_noise itself, from noise="):
        plain.run_ies(reals=3, ies_no_noise=True)
    with pytest.raises(ValueError, match=r"prior\(\) sets ies_no_noise itself"):
        plain.prior(reals=3, ies_no_noise=False)


# -- the correlated prior is seeded and the seed recorded -------------------------------


def test_one_seed_draws_one_prior_and_the_launch_records_it(correlated, monkeypatch):
    cal = correlated
    _stub_pestpp(monkeypatch, cal)
    ensemble = cal.template_workspace / "corr.prior_par.csv"
    before = np.random.get_state()[1].copy()

    cal.prior(reals=5, seed=11)
    first = ensemble.read_text()
    assert np.array_equal(np.random.get_state()[1], before), "numpy's generator was moved"
    cal.prior(reals=5, seed=11)
    assert ensemble.read_text() == first
    cal.run_ies(reals=5, iterations=1, seed=12)
    assert ensemble.read_text() != first

    cal.prior(reals=5)
    drawn = ensemble.read_text()
    records = _executions(cal)
    assert [record["prior_seed"] for record in records[:3]] == [11, 11, 12]
    unseeded = records[-1]["prior_seed"]
    assert isinstance(unseeded, int) and records[-1]["prior_ensemble"] == "corr.prior_par.csv"
    cal.prior(reals=5, seed=unseeded)
    assert ensemble.read_text() == drawn, "the recorded seed does not draw the run again"


def test_draw_prior_takes_a_seed_and_none_follows_numpys_generator(correlated):
    cal = correlated
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        one = cal.draw_prior(20, seed=3)._df
        two = cal.draw_prior(20, seed=3)._df
        np.random.seed(5)
        three = cal.draw_prior(20)._df
        np.random.seed(5)
        four = cal.draw_prior(20)._df
    pd.testing.assert_frame_equal(one, two)
    pd.testing.assert_frame_equal(three, four)
    assert not one.equals(three)
    with pytest.raises(TypeError, match=r"seed=True"):
        cal.draw_prior(5, seed=True)
    with pytest.raises(ValueError, match=r"2\*\*32 - 1"):
        cal.draw_prior(5, seed=-1)


def test_a_seed_where_nothing_is_drawn_is_refused_naming_random_seed(
        plain, correlated, monkeypatch):
    _stub_pestpp(monkeypatch, plain)
    with pytest.raises(ValueError, match=r"no grid or array pilotpoints parameter carries "
                                         r"correlation=.*pass random_seed=4"):
        plain.run_ies(reals=3, seed=4)
    _stub_pestpp(monkeypatch, correlated)
    with pytest.raises(ValueError, match=r"ies_par_en is passed"):
        correlated.prior(reals=3, seed=4, ies_par_en="mine.csv")


# -- a master anywhere is found through the launch's record -----------------------------


def test_a_master_given_elsewhere_is_found_and_survives_a_rebuild(plain, monkeypatch, tmp_path):
    cal = plain
    _stub_pestpp(monkeypatch, cal)
    master = tmp_path / "elsewhere" / "big_run"

    def start_workers(worker_dir, exe, pst, **kwargs):
        pd.DataFrame({"iteration": [0, 1], "total_runs": [3, 6], "mean": [9.0, 8.0]}).to_csv(
            Path(kwargs["master_dir"]) / "plain.phi.actual.csv", index=False)

    monkeypatch.setattr(cal.pyemu.os_utils, "start_workers", start_workers)
    cal.run_ies(reals=3, iterations=1, workers=2, master_dir=master)
    record = _executions(cal)[-1]
    assert record["directory"] == str(master.resolve()) and record["workers"] == 2

    runs = find_pest_runs(tmp_path / "pest")
    assert [run.name for run in runs] == ["plain"]
    assert runs[0].executions == {"ies": master.resolve()}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("plain.pst")
    assert find_pest_runs(tmp_path / "pest")[0].executions == {"ies": master.resolve()}
    # Scanning a root that holds the master too lists the run once, not the copy.
    assert [run.name for run in find_pest_runs(tmp_path)] == ["plain"]
