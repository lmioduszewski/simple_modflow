"""A PESTPP-IES run's review reads THAT run, and a long run is launched to survive.

* ledger 507: a run directory holding an earlier run's outputs reviewed that
  run's last iteration as the posterior (a rerun with fewer iterations, a prior after a
  run), and a ``noise=False`` run its noise ensemble;
* ledger 508: ``model.pest_runs`` opened the BUILD's ``.pst`` stem while
  PESTPP-IES writes ``<name>.*``, found masters by a name-prefix glob and classified
  them by a substring, never listed a serial run, fell back silently to another kind,
  and listed an interrupted run's worker directories as runs of their own;
* ledger 509: a caller's ``ies_restart_obs_en`` was kept while a correlated
  prior was drawn into ``ies_par_en``, so PESTPP-IES restarted from observations of
  other parameters;
* ledger 510: PESTPP-IES's ``overdue_giveup_fac`` was left at its own 2.0;
* ledger 511: pyEMU's ``chunk_len`` was left at 50, which spawns a
  20-process pool per realization at 50 or more multiplier files;
* ledger 512: ``ies_save_binary`` ensembles were invisible;
* ledger 513: every ensemble was re-read from disk on every property access;
* ledger 514: ``field()``'s mean/change/reduction were arithmetic on a
  log-normal field, and ``parameters_at_bounds`` took ``log10`` of an arithmetic mean;
* ledger 515: ``best('min_phi')`` returned a realization PESTPP-IES dropped.

Most tests drive ``IesResults`` and run discovery off hand-built directories. Two
launch tests stub PESTPP-IES itself; one runs it for real on the small toy of
``test_pest_handed_targets.py``.
"""

from __future__ import annotations

import json
import os
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest import project as project_module
from myflopy.modflow.mf6.pest.ies import IesResults
from myflopy.modflow.mf6.pest.runs import (
    BUILD_MARKER,
    METADATA_FILENAME,
    execution_kind,
    find_pest_runs,
    pestpp_ies_outputs,
)
from tests.test_pest_handed_targets import run_toy

pyemu = pytest.importorskip("pyemu")

OBS = ["h_a", "h_b"]
PARS = ["p1", "p2"]


def _control_file(directory: Path, case: str, *, noptmax: int = 2, **options) -> Path:
    """A real (if tiny) version-2 control file ``<case>.pst`` in ``directory``."""

    directory.mkdir(parents=True, exist_ok=True)
    pst = pyemu.Pst.from_par_obs_names(PARS, OBS)
    pst.parameter_data.loc[:, "partrans"] = "log"
    pst.parameter_data.loc[:, "parlbnd"] = 0.1
    pst.parameter_data.loc[:, "parubnd"] = 10.0
    pst.parameter_data.loc[:, "parval1"] = 1.0
    pst.control_data.noptmax = noptmax
    for key, value in options.items():
        pst.pestpp_options[key] = value
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst.write(str(directory / f"{case}.pst"), version=2)
    return directory / f"{case}.pst"


def _obs(directory: Path, case: str, iteration: int, value: float, *, binary=False) -> None:
    frame = pd.DataFrame({name: [value, value + 1.0] for name in OBS}, index=["0", "base"])
    _write_ensemble(directory, f"{case}.{iteration}.obs", frame, "obs", binary=binary)


def _par(directory: Path, case: str, iteration: int, value: float, *, binary=False) -> None:
    frame = pd.DataFrame({name: [value, value] for name in PARS}, index=["0", "base"])
    _write_ensemble(directory, f"{case}.{iteration}.par", frame, "par", binary=binary)


def _write_ensemble(directory, stem, frame, kind, *, binary):
    if binary:
        pst = pyemu.Pst.from_par_obs_names(PARS, OBS)
        cls = pyemu.ObservationEnsemble if kind == "obs" else pyemu.ParameterEnsemble
        cls(pst=pst, df=frame).to_binary(str(directory / f"{stem}.jcb"))
    else:
        frame.to_csv(directory / f"{stem}.csv")


def _phi(directory: Path, case: str, iterations, columns=("0", "base"), rows=None) -> None:
    rows = rows or [[float(it), 1.0] for it in iterations]
    pd.DataFrame(
        [[it, 2, 1.0, 0.0, 1.0, 1.0, *row] for it, row in zip(iterations, rows, strict=True)],
        columns=["iteration", "total_runs", "mean", "standard_deviation", "min", "max",
                 *columns],
    ).to_csv(directory / f"{case}.phi.actual.csv", index=False)


# -- the review reads this run's iterations (ledger 507) ----------------------


def test_an_iteration_the_runs_phi_does_not_record_is_not_the_posterior(tmp_path):
    """The rerun wrote iterations 0 and 1; the earlier run's iteration 2 is still there."""

    _control_file(tmp_path, "cal")
    for iteration in (0, 1, 2):
        _obs(tmp_path, "cal", iteration, 10.0 * iteration)
        _par(tmp_path, "cal", iteration, 1.0 + iteration)
    _phi(tmp_path, "cal", [0, 1])

    review = IesResults(tmp_path, case_name="cal")
    assert review.iterations == [0, 1]
    assert review.posterior_iteration == 1
    assert float(review.posterior._df.loc["0", "h_a"]) == 10.0


def test_a_prior_reviews_its_one_iteration_beside_an_earlier_ies_run(tmp_path):
    _control_file(tmp_path, "cal", noptmax=-1)
    for iteration in (0, 1, 2):
        _obs(tmp_path, "cal", iteration, 10.0 * iteration)
    _phi(tmp_path, "cal", [0])

    review = IesResults(tmp_path, case_name="cal")
    assert review.iterations == [0]
    assert review.posterior_iteration == review.prior_iteration == 0


def test_a_run_without_noise_has_none_whatever_an_earlier_run_left(tmp_path):
    _control_file(tmp_path, "cal", ies_no_noise=True)
    _obs(tmp_path, "cal", 0, 1.0)
    _phi(tmp_path, "cal", [0])
    pd.DataFrame({name: [1.0] for name in OBS}, index=["0"]).to_csv(
        tmp_path / "cal.obs+noise.csv"
    )

    assert IesResults(tmp_path, case_name="cal").noise is None


def test_the_outputs_a_launch_clears_are_the_runs_and_nothing_else(tmp_path):
    names = [
        "cal.0.obs.csv", "cal.12.par.jcb", "cal.3.base.par", "cal.1.pcs.csv",
        "cal.phi.actual.csv", "cal.phi.meas.csv", "cal.obs+noise.csv", "cal.base.par",
        "cal.base.rei",
    ]
    kept = [
        "cal.pst", "cal.prior_par.csv", "cal.obs_data.csv", "cal.rec", "calx.0.obs.csv",
        "cal_ss.0.obs.csv", "other.0.obs.csv", "cal.par_data.csv",
    ]
    for name in names + kept:
        (tmp_path / name).write_text("x")
    (tmp_path / "cal.7.dir").mkdir()

    assert [path.name for path in pestpp_ies_outputs(tmp_path, "cal")] == sorted(names)


# -- binary ensembles and the cache (ledgers 512, 513) --------------


def test_a_binary_run_is_reviewed(tmp_path):
    _control_file(tmp_path, "cal", ies_save_binary=True)
    for iteration in (0, 1):
        _obs(tmp_path, "cal", iteration, 5.0 + iteration, binary=True)
        _par(tmp_path, "cal", iteration, 2.0, binary=True)
    _phi(tmp_path, "cal", [0, 1])

    review = IesResults(tmp_path, case_name="cal")
    assert review.iterations == [0, 1]
    assert float(review.posterior._df.loc["0", "h_a"]) == 6.0
    assert float(review.posterior_parameters._df.loc["base", "p1"]) == 2.0


def test_the_binary_form_wins_when_the_run_saved_binary(tmp_path):
    """A CSV an earlier, non-binary run left beside the binary run's own ensemble."""

    _control_file(tmp_path, "cal", ies_save_binary=True)
    _obs(tmp_path, "cal", 0, 1.0)
    _obs(tmp_path, "cal", 0, 7.0, binary=True)
    _phi(tmp_path, "cal", [0])

    assert float(IesResults(tmp_path, case_name="cal").prior._df.loc["0", "h_a"]) == 7.0


def test_an_ensemble_is_read_once_until_its_file_changes(tmp_path, monkeypatch):
    _control_file(tmp_path, "cal")
    _obs(tmp_path, "cal", 0, 1.0)
    _phi(tmp_path, "cal", [0])
    review = IesResults(tmp_path, case_name="cal")

    reads = []
    original = pyemu.ObservationEnsemble.from_csv.__func__

    def counting(cls, *args, **kwargs):
        reads.append(1)
        return original(cls, *args, **kwargs)

    monkeypatch.setattr(pyemu.ObservationEnsemble, "from_csv", classmethod(counting))
    for _ in range(4):
        review.prior
        review.posterior
    assert len(reads) == 1

    _obs(tmp_path, "cal", 0, 3.0)
    path = tmp_path / "cal.0.obs.csv"
    stat = path.stat()
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10**9))
    assert float(review.prior._df.loc["0", "h_a"]) == 3.0
    assert len(reads) == 2


# -- best('min_phi') (ledger 515) --------------------------------------------


def _stub(phi_frame, posterior_labels):
    ensemble = pd.DataFrame({"h_a": 1.0}, index=posterior_labels)

    class _Stub(IesResults):
        def __init__(self):
            pass

        phi = property(lambda self: phi_frame)
        posterior = property(lambda self: SimpleNamespace(_df=ensemble))
        posterior_parameters = property(lambda self: SimpleNamespace(_df=ensemble))

    return _Stub()


def test_min_phi_never_returns_a_realization_pestpp_dropped():
    """An empty phi cell first: `min` over NaN kept it, and its label has no posterior."""

    phi = pd.DataFrame(
        [[1, 9, 1.0, 0.0, 1.0, 1.0, np.nan, 5.0, 3.0, 4.0]],
        columns=["iteration", "total_runs", "mean", "standard_deviation", "min", "max",
                 "0", "1", "2", "base"],
    )
    assert _stub(phi, ["1", "2", "base"]).best(criterion="min_phi") == "2"
    # A realization with a phi but no posterior parameter set is not carried forward.
    assert _stub(phi, ["1", "base"]).best(criterion="min_phi") == "base"
    with pytest.raises(ValueError, match="dropped or failed every one"):
        _stub(phi, ["0"]).best(criterion="min_phi")


# -- field statistics in the parameters' space (ledger 514) ------------------


def _field_stub(transform: str, prior: np.ndarray, posterior: np.ndarray):
    names = [f"oname:kfieldl0_otype:arr_i:{i}_j:0" for i in range(prior.shape[1])]
    parameter_data = pd.DataFrame(
        {"parnme": ["kl0pp0000"], "pargp": ["kl0"], "partrans": [transform],
         "parlbnd": [0.1], "parubnd": [10.0]}
    ).set_index("parnme")

    class _Stub(IesResults):
        def __init__(self):
            self.pst = SimpleNamespace(parameter_data=parameter_data)

        capture_fields = property(lambda self: [
            {"target": "k", "name": "k", "family": "array", "prefix": "kfield",
             "layer_prefixes": {"0": "kfieldl0"}}
        ])
        prior = property(lambda self: SimpleNamespace(_df=pd.DataFrame(prior, columns=names)))
        posterior = property(
            lambda self: SimpleNamespace(_df=pd.DataFrame(posterior, columns=names))
        )

        def _field_cell_index(self, info, layer):
            return {name: index for index, name in enumerate(names)}

    return _Stub()


def _lognormal(sd, *, cells=50, reals=400, seed=0):
    rng = np.random.default_rng(seed)
    return 10.0 ** rng.normal(1.0, sd, size=(reals, cells))


def test_an_unchanged_median_maps_as_no_change():
    """Prior log10 K ~ N(1, 0.5), posterior N(1, 0.1): the data narrowed K, not lowered it."""

    frame = _field_stub("log", _lognormal(0.5), _lognormal(0.1, seed=1)).field("k")

    assert frame.attrs["space"] == "log10"
    # Arithmetic means read this as a ~0.54x drop (the finding's 0.535).
    assert np.median(frame["change"]) == pytest.approx(1.0, abs=0.05)
    assert np.median(frame["reduction"]) == pytest.approx(0.8, abs=0.03)
    assert np.median(frame["prior_mean"]) == pytest.approx(10.0, rel=0.1)
    assert np.median(frame["posterior_std"]) == pytest.approx(0.1, rel=0.1)


def test_the_field_reduction_is_sensitivitys_learned():
    """One parameter per cell: `field`'s reduction and `sensitivity`'s learned agree."""

    prior, posterior = _lognormal(0.5, cells=3), _lognormal(0.2, cells=3, seed=2)
    review = _field_stub("log", prior, posterior)
    parameter_data = pd.DataFrame(
        {"parnme": ["a", "b", "c"], "pargp": ["k"] * 3, "partrans": ["log"] * 3}
    ).set_index("parnme")

    class _Sensitivity(IesResults):
        def __init__(self):
            self.pst = SimpleNamespace(parameter_data=parameter_data)

        prior_parameters = property(
            lambda self: SimpleNamespace(_df=pd.DataFrame(prior, columns=["a", "b", "c"])))
        posterior_parameters = property(
            lambda self: SimpleNamespace(_df=pd.DataFrame(posterior, columns=["a", "b", "c"])))

    learned = _Sensitivity().sensitivity(by_group=False).loc[["a", "b", "c"], "learned"]
    assert review.field("k")["reduction"].to_numpy() == pytest.approx(learned.to_numpy())


def test_a_linear_parameter_keeps_linear_statistics():
    prior, posterior = _lognormal(0.5), _lognormal(0.1, seed=1)
    frame = _field_stub("none", prior, posterior).field("k")

    assert frame.attrs["space"] == "linear"
    assert frame["prior_mean"].to_numpy() == pytest.approx(prior.mean(axis=0))


def test_a_log_parameter_at_its_bound_is_judged_by_its_geometric_mean():
    """Three realizations at the upper bound and one at the lower: the arithmetic mean,
    7.5, sat at 94% of the log range and counted at the bound with tol=0.1."""

    parameter_data = pd.DataFrame(
        {"parnme": ["p"], "pargp": ["g"], "partrans": ["log"], "parlbnd": [0.1],
         "parubnd": [10.0]}
    ).set_index("parnme")
    values = pd.DataFrame({"p": [10.0, 10.0, 10.0, 0.1]})

    class _Stub(IesResults):
        def __init__(self):
            self.pst = SimpleNamespace(parameter_data=parameter_data)

        posterior_parameters = property(lambda self: SimpleNamespace(_df=values))
        prior_parameters = property(lambda self: SimpleNamespace(_df=values))

    assert int(_Stub().parameters_at_bounds(tol=0.1).loc["g", "n_at_upper"]) == 0


# -- run discovery (ledger 508) ----------------------------------------------


def _build_dir(directory: Path, name: str, *, pst_file: str | None = None) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / METADATA_FILENAME).write_text(json.dumps(
        {"project_name": name, "model_name": "m", "pst_file": pst_file or f"{name}.pst"}
    ))
    (directory / BUILD_MARKER).write_text(json.dumps({"directory": str(directory.resolve())}))
    return directory


def _master(template: Path, name: str, suffix: str, *, marker: bool = True) -> Path:
    master = template.parent / f"{name}_{suffix}"
    shutil.copytree(template, master)
    (master / f"{name}.pst").write_text("")
    if not marker:  # a master copied from a template built before ledger 360
        (master / BUILD_MARKER).unlink()
    return master


def test_a_calibration_holds_only_its_own_masters(tmp_path):
    """`tt` was given `tt_ss`'s and `tt_prior_k`'s masters; `tt_prior_k`'s IES master
    was classified 'prior' by the substring in its name."""

    tt = _build_dir(tmp_path / "tt", "tt")
    ss = _build_dir(tmp_path / "tt_ss", "tt_ss")
    pk = _build_dir(tmp_path / "tt_prior_k", "tt_prior_k")
    # Without markers, as a template built before ledger 360 left them: the NAME alone
    # must tell whose master each is.
    _master(ss, "tt_ss", "ies_master", marker=False)
    pk_ies = _master(pk, "tt_prior_k", "ies_master", marker=False)
    pk_prior = _master(pk, "tt_prior_k", "prior_master", marker=False)

    runs = {run.name: run for run in find_pest_runs(tmp_path)}
    assert set(runs) == {"tt", "tt_ss", "tt_prior_k"}
    assert runs["tt"].executions == {}
    assert runs["tt_prior_k"].executions == {"ies": pk_ies, "prior": pk_prior}
    assert runs[tt.name].kinds == []


def test_a_master_copied_from_another_template_is_not_this_runs(tmp_path):
    first = _build_dir(tmp_path / "a" / "cal", "cal")
    second = _build_dir(tmp_path / "b" / "cal", "cal")
    stray = second.parent / "cal_ies_master"
    shutil.copytree(first, stray)
    (stray / "cal.pst").write_text("")

    assert find_pest_runs(tmp_path / "b")[0].executions == {}


def test_review_opens_the_case_pestpp_ran_not_the_builds_file(tmp_path):
    """`model.pest('valley_pest')` + `cal.build('valley.pst')`: PESTPP-IES writes valley_pest.*."""

    template = _build_dir(tmp_path / "valley_pest", "valley_pest", pst_file="valley.pst")
    _control_file(template, "valley")
    _control_file(template, "valley_pest", noptmax=-1)
    _obs(template, "valley_pest", 0, 1.0)
    _phi(template, "valley_pest", [0])

    run = find_pest_runs(tmp_path)[0]
    assert run.case == "valley_pest"
    assert run.executions == {"prior": template}  # a SERIAL run is listed
    assert run.review("prior").iterations == [0]
    with pytest.raises(ValueError, match=r"review\(kind='prior'\)"):
        run.review()  # no silent fallback from 'ies' to the prior


def test_a_serial_run_and_a_master_of_one_kind_open_the_later(tmp_path):
    template = _build_dir(tmp_path / "cal", "cal")
    _control_file(template, "cal", noptmax=3)
    _phi(template, "cal", [0])
    master = _master(template, "cal", "ies_master")
    _phi(master, "cal", [0, 1])
    os.utime(template / "cal.phi.actual.csv", (1, 1))

    assert find_pest_runs(tmp_path)[0].executions == {"ies": master}
    os.utime(master / "cal.phi.actual.csv", (0, 0))
    assert find_pest_runs(tmp_path)[0].executions == {"ies": template}


def test_execution_kind_reads_the_launched_noptmax(tmp_path):
    _control_file(tmp_path, "a", noptmax=-1)
    _control_file(tmp_path, "b", noptmax=3)
    _control_file(tmp_path, "c", noptmax=3)
    for case in ("a", "b"):
        _phi(tmp_path, case, [0])

    assert [execution_kind(tmp_path, case) for case in "abc"] == ["prior", "ies", None]


def test_an_interrupted_runs_workers_are_not_listed_as_runs(tmp_path):
    template = _build_dir(tmp_path / "tt", "tt")
    for index in range(3):
        shutil.copytree(template, tmp_path / f"worker_{index}")
    # A template built before ledger 360 carried no marker, nor did its workers.
    old = tmp_path / "legacy" / "old"
    old.mkdir(parents=True)
    (old / METADATA_FILENAME).write_text(json.dumps({"project_name": "old"}))
    shutil.copytree(old, old.parent / "worker_0")

    assert [(run.name, run.template_dir) for run in find_pest_runs(tmp_path)] == [
        ("old", old), ("tt", template)
    ]



def test_a_moved_calibration_root_keeps_its_runs(tmp_path):
    """A marker keeps the path its build wrote, so after the root moves EVERY template's
    names a directory that is not its own; only a copy names one still there (ledger
    516). Moved: the run, its master, and not its leftover worker."""

    template = _build_dir(tmp_path / "old" / "tt", "tt")
    master = _master(template, "tt", "ies_master")
    shutil.copytree(template, template.parent / "worker_0")
    shutil.move(tmp_path / "old", tmp_path / "new")
    moved = tmp_path / "new" / "tt"

    runs = find_pest_runs(tmp_path / "new")
    assert [(run.name, run.template_dir) for run in runs] == [("tt", moved)]
    assert runs[0].executions == {"ies": moved.parent / master.name}


def test_a_copied_calibration_root_keeps_its_runs_and_a_backup_beside_one_is_not_a_run(
    tmp_path,
):
    template = _build_dir(tmp_path / "project" / "tt", "tt")
    shutil.copytree(tmp_path / "project", tmp_path / "copy")
    # The copy's marker names the ORIGINAL, which is outside the scanned root.
    assert [run.template_dir for run in find_pest_runs(tmp_path / "copy")] == [
        tmp_path / "copy" / "tt"
    ]
    # A backup beside the template names a build this scan sees: a copy, not a run.
    shutil.copytree(template, template.parent / "tt_backup")
    assert [run.template_dir for run in find_pest_runs(tmp_path / "project")] == [template]

# -- the launch (ledgers 507, 509, 510, 511) ------


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return run_toy(tmp_path_factory.mktemp("pest_runs_review"))


HEADS = pd.DataFrame({"time": [0, 1, 2], "W2": [91.0, 92.0, 93.0],
                      "W1": [93.0, 93.5, 94.0], "W10": [95.0, 96.0, 97.0]})


def _calibration(toy, workspace: Path, name: str, *, correlated: bool = False):
    locations = pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0})
    cal = toy.pest(name, workspace=workspace, start_datetime="2024-01-01")
    if correlated:
        cal.parameterize("k", style="grid", correlation=300.0, bounds=(0.5, 2.0),
                         physical=(1e-3, 100.0))
    else:
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(mf.HeadTargets(locations=locations, values=HEADS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}_build.pst", noptmax=0)
    return cal


@pytest.fixture(scope="module")
def correlated(toy, tmp_path_factory):
    return _calibration(toy, tmp_path_factory.mktemp("correlated") / "corr", "corr",
                        correlated=True)


def _stub_pestpp(monkeypatch, cal):
    """PESTPP-IES replaced by a recorder, and the ensembles an earlier test named cleared."""

    for option in ("ies_par_en", "ies_restart_obs_en"):
        cal.pst.pestpp_options.pop(option, None)
    commands = []
    monkeypatch.setattr(
        cal.pyemu.os_utils, "run", lambda command, cwd=None: commands.append((command, cwd))
    )
    return commands


def test_a_restart_ensemble_is_refused_beside_a_drawn_prior(correlated, monkeypatch):
    cal = correlated
    _stub_pestpp(monkeypatch, cal)
    with pytest.raises(ValueError, match=r"ies_restart_obs_en='my\.obs\.csv'.*ies_par_en"):
        cal.prior(reals=3, ies_restart_obs_en="my.obs.csv")
    # It persists on cal.pst from the refused call, so the next is refused too ...
    with pytest.raises(ValueError, match="del cal.pst.pestpp_options"):
        cal.prior(reals=3)
    # ... until it is paired with the parameters that made it, which are used as given.
    cal.prior(reals=3, ies_par_en="my.par.csv", ies_restart_obs_en="my.obs.csv")
    assert cal.pst.pestpp_options["ies_par_en"] == "my.par.csv"
    del cal.pst.pestpp_options["ies_restart_obs_en"]
    cal.prior(reals=3)
    assert cal.pst.pestpp_options["ies_par_en"] == "corr.prior_par.csv"


def test_overdue_giveup_fac_defaults_to_five_and_is_overridable(correlated, monkeypatch):
    cal = correlated
    _stub_pestpp(monkeypatch, cal)
    written = cal.template_workspace / "corr.pst"

    cal.run_ies(reals=3, iterations=1)
    assert cal.pst.pestpp_options["overdue_giveup_fac"] == 5.0
    assert "overdue_giveup_fac" in written.read_text()
    cal.prior(reals=3, overdue_giveup_fac=12)
    assert cal.pst.pestpp_options["overdue_giveup_fac"] == 12.0
    cal.run_ies(reals=3, overdue_giveup_fac=None)
    assert "overdue_giveup_fac" not in cal.pst.pestpp_options
    assert "overdue_giveup_fac" not in written.read_text()


def test_every_multiplier_file_is_applied_in_the_forward_runs_own_process(correlated):
    text = (correlated.template_workspace / "forward_run.py").read_text()
    assert f"chunk_len={project_module.PYEMU_APPLY_CHUNK_LEN}" in text
    assert "chunk_len=50" not in text


def test_a_serial_launch_clears_the_earlier_runs_outputs(correlated, monkeypatch):
    cal = correlated
    _stub_pestpp(monkeypatch, cal)
    template = cal.template_workspace
    stale = ["corr.2.obs.csv", "corr.2.par.csv", "corr.phi.actual.csv", "corr.obs+noise.csv",
             "corr.2.base.par"]
    for name in stale + ["corr.0.par.csv", "corr_build.0.obs.csv"]:
        (template / name).write_text("x")

    cal.run_ies(reals=3, iterations=1, ies_par_en="corr.0.par.csv")

    assert not any((template / name).exists() for name in stale)
    assert (template / "corr.0.par.csv").exists(), "a file an option names is kept"
    assert (template / "corr_build.0.obs.csv").exists(), "another case's is not this run's"


def test_a_parallel_master_never_holds_the_templates_serial_outputs(
    correlated, monkeypatch, tmp_path
):
    """pyEMU copied the whole template into the master, an earlier serial run's outputs
    included, and ledger 507 removed the ones PESTPP-IES did not rewrite after the run.
    The launch makes the master itself now, without them (ledger 532)."""

    cal = correlated
    _stub_pestpp(monkeypatch, cal)
    template = cal.template_workspace
    for name in ("corr.2.obs.csv", "corr.phi.actual.csv", "corr.0.obs.csv"):
        (template / name).write_text("earlier serial run")
    master = tmp_path / "corr_ies_master"

    def start_workers(worker_dir, exe, pst, **kwargs):
        assert kwargs["reuse_master"] is True and Path(kwargs["master_dir"]) == master
        assert not (master / "corr.2.obs.csv").exists()
        assert not (master / "corr.phi.actual.csv").exists()
        (master / "corr.0.obs.csv").write_text("this run")
        raise KeyboardInterrupt  # a run stopped mid-way

    monkeypatch.setattr(cal.pyemu.os_utils, "start_workers", start_workers)
    with pytest.raises(KeyboardInterrupt):
        cal.run_ies(reals=3, workers=2, master_dir=master)

    assert not (master / "corr.2.obs.csv").exists()
    assert (master / "corr.0.obs.csv").read_text() == "this run"
    assert (template / "corr.2.obs.csv").exists(), "the template keeps its serial run"
    for name in ("corr.2.obs.csv", "corr.phi.actual.csv", "corr.0.obs.csv"):
        (template / name).unlink()


def test_reruns_in_one_template_review_their_own_iterations(toy, tmp_path):
    """Real PESTPP-IES, serially in one template, the manual's naming: pest('valley_pest')
    + build('valley.pst')."""

    locations = pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0})
    cal = toy.pest("valley_pest", workspace=tmp_path / "valley_pest",
                   start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(mf.HeadTargets(locations=locations, values=HEADS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("valley.pst", noptmax=0)
    fast = {"ies_lambda_mults": 1.0, "lambda_scale_fac": 1.0}

    first = cal.run_ies(reals=3, iterations=2, **fast)
    assert first.iterations == [0, 1, 2]
    with pytest.warns(UserWarning, match=r"run_ies\(\) deleted \d+ output file"):
        second = cal.run_ies(reals=3, iterations=1, noise=False, ies_save_binary=True, **fast)
    assert second.iterations == [0, 1]
    assert second.posterior_iteration == 1
    assert not (cal.template_workspace / "valley_pest.2.obs.csv").exists()
    assert (cal.template_workspace / "valley_pest.1.obs.jcb").exists()
    assert second.noise is None
    assert len(second.posterior._df) >= 1

    with pytest.warns(UserWarning, match=r"prior\(\) deleted \d+ output file"):
        prior = cal.prior(reals=3)
    assert prior.iterations == [0]

    handle = find_pest_runs(tmp_path)[0]
    assert handle.case == "valley_pest"
    assert handle.kinds == ["prior"]
    assert handle.review("prior").iterations == [0]
