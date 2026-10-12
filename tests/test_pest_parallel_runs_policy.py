"""A PESTPP-IES run keeps its own workers, never deletes finished results unasked, refuses
a live run's folders, cleans up after itself and can resume (ledger 532).

Until then (ledger 361):

* every parallel launch of every calibration of one model made its workers in the
  master's PARENT, ``<model>.pest/worker_<n>`` -- two calibrations, or one calibration's
  ``prior()`` and ``run_ies()`` in two kernels, shared them, and a launch deleted a live
  run's workers (they were copies of the same template, so the guard let them go);
* a relaunch deleted the previous run's master -- its results, a multi-day run's
  included -- with no flag;
* nothing could tell a run still going from a finished one;
* pyEMU deleted the workers whenever its master exited, failed or not, and copied the
  template's earlier serial results into the master;
* an interrupted run could not be continued: PESTPP-IES refuses the ``/r`` switch, and
  nothing set up its ensemble restart;
* a serial run deleted the template's earlier results with a log line nobody saw.

The model: the small toy of ``tests/test_pest_handed_targets.py`` (a 6 x 6 Voronoi grid,
three periods), a constant-K calibration of three head targets. The policy tests launch
pyEMU's own ``start_workers`` with ``true`` as PEST++ -- it exits at once, so a run
"finishes" without PEST -- and the results reader stubbed; the resume tests run
PESTPP-IES for real, with at most 2 workers.
"""

from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import time
import warnings
from pathlib import Path

import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.run_state import (
    RESUME_RECORD,
    RUN_LOCK,
    archive_for_resume,
    pending_resume,
    run_activity,
    splice_resumed_run,
    worker_root,
    write_run_lock,
)
from myflopy.modflow.mf6.pest.runs import find_pest_runs
from tests.test_pest_handed_targets import run_toy

pyemu = pytest.importorskip("pyemu")

_TRUE = shutil.which("true")
needs_true = pytest.mark.skipif(_TRUE is None, reason="needs a `true` executable on PATH")
needs_pestpp = pytest.mark.skipif(
    shutil.which("pestpp-ies") is None, reason="needs pestpp-ies on PATH"
)

HEADS = pd.DataFrame({"time": [0, 1, 2], "W2": [91.0, 92.0, 93.0],
                      "W1": [93.0, 93.5, 94.0], "W10": [95.0, 96.0, 97.0]})
LOCATIONS = pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0})
#: PESTPP-IES options that make an iteration one lambda test, for speed.
FAST = {"ies_lambda_mults": 1.0, "lambda_scale_fac": 1.0}


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return run_toy(tmp_path_factory.mktemp("parallel_runs_policy"))


def _calibration(toy, name: str, workspace: Path | None = None):
    cal = toy.pest(name, workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(mf.HeadTargets(locations=LOCATIONS, values=HEADS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}.pst", noptmax=0)
    return cal


@pytest.fixture
def cal(toy, tmp_path):
    """A calibration built at its own template under ``tmp_path``."""

    built = _calibration(toy, "cal", tmp_path / "pest" / "cal")
    yield built
    built._release_pyemu_log()


@pytest.fixture
def _no_results_reader(monkeypatch):
    """The run's results are not read: `true` wrote none."""

    import myflopy.modflow.mf6.pest.ies as ies

    monkeypatch.setattr(ies, "IesResults", lambda workspace, **_: Path(workspace))


@pytest.fixture
def spy_start_workers(cal, monkeypatch):
    """pyEMU's own ``start_workers``, recording each call's keyword arguments."""

    calls = []
    real = cal.pyemu.os_utils.start_workers

    def start_workers(*args, **kwargs):
        calls.append(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(cal.pyemu.os_utils, "start_workers", start_workers)
    return calls


def _results(directory: Path, case: str, iterations=(0, 1), noptmax: int | None = None):
    """What a PESTPP-IES run of ``case`` leaves: ensembles per iteration and a phi file."""

    directory.mkdir(parents=True, exist_ok=True)
    for n in iterations:
        for kind in ("par", "obs"):
            pd.DataFrame({"x": [float(n)]}, index=["0"]).to_csv(directory / f"{case}.{n}.{kind}.csv")
    pd.DataFrame({"iteration": list(iterations), "total_runs": [3 * (n + 1) for n in iterations],
                  "mean": [10.0 - n for n in iterations]}).to_csv(
        directory / f"{case}.phi.actual.csv", index=False)
    if noptmax is not None:
        text = (directory / f"{case}.pst").read_text(encoding="utf-8")
        lines = [
            f"noptmax {noptmax}" if line.split()[:1] == ["noptmax"] else line
            for line in text.splitlines()
        ]
        (directory / f"{case}.pst").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _copy_of_template(cal, directory: Path) -> Path:
    """What a launch makes of the template: a copy, build marker and all."""

    shutil.copytree(cal.template_workspace, directory)
    return directory


# -- workers per master -------------------------------------------------------------


@needs_true
def test_workers_live_beside_their_own_master_and_go_when_the_run_finishes(
    cal, spy_start_workers, _no_results_reader
):
    master = cal.run_ies(workers=2, exe="true")
    assert master == cal.template_workspace.parent / "cal_ies_master"
    (call,) = spy_start_workers
    assert Path(call["worker_root"]) == worker_root(master) == master.with_name("cal_ies_master.workers")
    assert call["reuse_master"] is True and call["cleanup"] is False
    assert not worker_root(master).exists(), "a clean finish deletes the workers"
    assert not (master / RUN_LOCK).exists(), "the lock goes with the run"
    assert not list(cal.template_workspace.parent.glob("worker_*"))


@needs_true
def test_keep_workers_keeps_them(cal, _no_results_reader):
    master = cal.prior(workers=2, exe="true", keep_workers=True)
    root = worker_root(master)
    assert sorted(path.name for path in root.iterdir() if path.is_dir()) == ["worker_0", "worker_1"]
    assert not (root / RUN_LOCK).exists()
    # Each worker carries the template's build metadata; discovery skips them as copies.
    assert [handle.name for handle in find_pest_runs(cal.template_workspace.parent)] == ["cal"]


@needs_true
def test_a_failed_run_keeps_its_workers(cal, monkeypatch, _no_results_reader):
    false = shutil.which("false")
    if false is None:
        pytest.skip("needs a `false` executable on PATH")
    with pytest.raises(Exception, match="master returned non-zero"):
        cal.run_ies(workers=2, exe="false")
    root = worker_root(cal.template_workspace.parent / "cal_ies_master")
    assert (root / "worker_0").is_dir() and (root / "worker_1").is_dir()
    assert not (root / RUN_LOCK).exists()


@needs_true
def test_two_calibrations_and_a_prior_beside_a_history_match_share_no_worker(
    toy, tmp_path, monkeypatch, _no_results_reader
):
    """Every calibration of a model used ``<model>.pest/worker_<n>``; a second one's run
    deleted the first's workers, live or not."""

    first = _calibration(toy, "first", tmp_path / "pest" / "first")
    second = _calibration(toy, "second", tmp_path / "pest" / "second")
    roots = []
    real = first.pyemu.os_utils.start_workers
    monkeypatch.setattr(first.pyemu.os_utils, "start_workers",
                        lambda *a, **k: roots.append(k["worker_root"]) or real(*a, **k))
    # The first calibration's prior is running: its workers carry a live lock.
    live = worker_root(first.template_workspace.parent / "first_prior_master")
    (live / "worker_0").mkdir(parents=True)
    shutil.copytree(first.template_workspace, live / "worker_0", dirs_exist_ok=True)
    write_run_lock(live, method="prior", calibration="first")
    before = sorted(str(path) for path in live.rglob("*"))

    second.prior(workers=2, exe="true")
    second.run_ies(workers=2, exe="true")
    first.run_ies(workers=2, exe="true")

    assert len(set(roots)) == 3
    assert str(live) not in roots
    assert sorted(str(path) for path in live.rglob("*")) == before
    (live / RUN_LOCK).unlink()
    for each in (first, second):
        each._release_pyemu_log()


# -- finished results are not deleted unasked ---------------------------------------


@needs_true
def test_a_master_holding_a_finished_run_is_refused_unless_overwrite(cal, _no_results_reader):
    master = _copy_of_template(cal, cal.template_workspace.parent / "cal_ies_master")
    _results(master, "cal", iterations=(0, 1, 2), noptmax=2)

    with pytest.raises(FileExistsError, match=r"holds a finished history match \(iterations \[0, 1, 2\]\)") as refused:
        cal.run_ies(workers=2, exe="true")
    message = str(refused.value)
    assert str(master) in message
    assert "overwrite=True" in message and "master_dir=<a new or empty folder>" in message
    assert "resume_ies" not in message
    assert (master / "cal.2.obs.csv").exists()

    assert cal.run_ies(workers=2, exe="true", overwrite=True) == master
    assert not (master / "cal.2.obs.csv").exists()


@needs_true
def test_an_interrupted_masters_refusal_names_resume(cal, _no_results_reader):
    master = _copy_of_template(cal, cal.template_workspace.parent / "cal_ies_master")
    _results(master, "cal", iterations=(0, 1), noptmax=3)
    with pytest.raises(FileExistsError, match="INTERRUPTED history match") as refused:
        cal.run_ies(workers=2, exe="true")
    assert "cal.resume_ies(workers=2)" in str(refused.value)


@needs_true
def test_a_prior_master_holding_results_is_refused_unless_overwrite(cal, _no_results_reader):
    master = _copy_of_template(cal, cal.template_workspace.parent / "cal_prior_master")
    _results(master, "cal", iterations=(0,), noptmax=-1)
    with pytest.raises(FileExistsError, match="prior Monte Carlo's results"):
        cal.prior(workers=2, exe="true")
    cal.prior(workers=2, exe="true", overwrite=True)


@needs_true
def test_another_calibrations_master_is_refused_even_with_overwrite(toy, tmp_path, _no_results_reader):
    first = _calibration(toy, "first", tmp_path / "pest" / "first")
    second = _calibration(toy, "second", tmp_path / "pest" / "second")
    master = _copy_of_template(first, tmp_path / "shared_master")
    with pytest.raises(FileExistsError, match="not a copy of this calibration's template"):
        second.run_ies(workers=2, exe="true", master_dir=master, overwrite=True)
    for each in (first, second):
        each._release_pyemu_log()


# -- a live run's folders are refused --------------------------------------------------


@needs_true
def test_a_master_a_run_is_going_in_is_refused_whatever_overwrite_says(cal, _no_results_reader):
    master = _copy_of_template(cal, cal.template_workspace.parent / "cal_ies_master")
    write_run_lock(master, method="run_ies", calibration="cal")   # this process: alive
    before = sorted(path.name for path in master.iterdir())
    with pytest.raises(RuntimeError, match="may still be going -- a run is going in it: run_ies") as refused:
        cal.run_ies(workers=2, exe="true", overwrite=True)
    assert f"process {os.getpid()}" in str(refused.value)
    assert f"delete {master / RUN_LOCK}" in str(refused.value)
    assert sorted(path.name for path in master.iterdir()) == before


@needs_true
def test_a_workers_folder_a_run_is_going_in_is_refused(cal, _no_results_reader):
    root = worker_root(cal.template_workspace.parent / "cal_ies_master")
    root.mkdir()
    write_run_lock(root, method="run_ies", calibration="cal")
    with pytest.raises(RuntimeError, match="worker directory"):
        cal.run_ies(workers=2, exe="true")
    (root / RUN_LOCK).unlink()


def _lock_of_a_finished_process(directory: Path) -> None:
    """A lock naming a process that has exited: the kernel that wrote it was killed."""

    process = subprocess.Popen(["true"])
    process.wait()
    write_run_lock(directory, method="run_ies", calibration="cal")
    lock = json.loads((directory / RUN_LOCK).read_text(encoding="utf-8"))
    lock["pid"] = process.pid
    (directory / RUN_LOCK).write_text(json.dumps(lock), encoding="utf-8")


@needs_true
def test_a_stale_lock_is_said_so_and_ignored(cal, _no_results_reader):
    master = _copy_of_template(cal, cal.template_workspace.parent / "cal_ies_master")
    _lock_of_a_finished_process(master)
    assert run_activity(master).state == "stale"
    with pytest.warns(UserWarning, match="STALE") as said:
        cal.run_ies(workers=2, exe="true")
    assert [w.filename for w in said if "STALE" in str(w.message)] == [__file__]
    assert not (master / RUN_LOCK).exists()


@pytest.mark.skipif(not Path("/proc/self/stat").exists(), reason="reads /proc")
def test_a_lock_whose_pid_now_names_another_process_is_stale(tmp_path):
    """A pid is reused: the lock records when its process started, and a running process
    with that pid that started at another time is not the run's."""

    write_run_lock(tmp_path, method="run_ies", calibration="cal")
    assert run_activity(tmp_path).state == "live"
    lock = json.loads((tmp_path / RUN_LOCK).read_text(encoding="utf-8"))
    lock["process_start"] = str(int(lock["process_start"]) - 1)
    (tmp_path / RUN_LOCK).write_text(json.dumps(lock), encoding="utf-8")
    assert run_activity(tmp_path).state == "stale"


def test_a_lock_from_another_host_counts_as_live(tmp_path):
    write_run_lock(tmp_path, method="run_ies", calibration="cal")
    lock = json.loads((tmp_path / RUN_LOCK).read_text(encoding="utf-8"))
    lock["host"] = socket.gethostname() + "-elsewhere"
    (tmp_path / RUN_LOCK).write_text(json.dumps(lock), encoding="utf-8")
    activity = run_activity(tmp_path)
    assert activity.state == "unknown" and activity.blocks
    assert "cannot check" in activity.describe()


def test_a_copied_lock_names_another_folder_and_is_not_this_ones(tmp_path):
    """pyEMU copies a template -- a serial run's lock and all -- into each worker."""

    template, copy = tmp_path / "template", tmp_path / "worker_0"
    template.mkdir()
    write_run_lock(template, method="run_ies", calibration="cal")
    shutil.copytree(template, copy)
    assert run_activity(template).state == "live"
    assert run_activity(copy).state == "free"


@pytest.mark.skipif(not Path("/proc/self/cwd").exists(), reason="reads /proc")
def test_a_pestpp_process_working_in_a_folder_makes_it_live_without_a_lock(tmp_path):
    """The kernel that launched a run was killed and PESTPP-IES kept going: no lock to
    read, but a PEST++ process works there. A shell left there is not a run."""

    sleep = shutil.which("sleep")
    if sleep is None:
        pytest.skip("needs `sleep`")
    # A script: the kernel names its process after the file, `pestpp-ies`.
    fake = tmp_path / "bin" / "pestpp-ies"
    fake.parent.mkdir()
    fake.write_text(f"#!/bin/sh\n{sleep} \"$1\"\n", encoding="utf-8")
    fake.chmod(0o755)
    run_dir, shell_dir = tmp_path / "master", tmp_path / "other"
    run_dir.mkdir()
    shell_dir.mkdir()
    pestpp = subprocess.Popen([str(fake), "30"], cwd=run_dir)
    shell = subprocess.Popen([sleep, "30"], cwd=shell_dir)
    try:
        deadline = time.monotonic() + 5
        while run_activity(run_dir).state != "live" and time.monotonic() < deadline:
            time.sleep(0.05)
        activity = run_activity(run_dir)
        assert activity.state == "live" and activity.processes[0][0] == pestpp.pid
        assert "pestpp-ies" in activity.describe()
        assert run_activity(shell_dir).state == "free"
    finally:
        for process in (pestpp, shell):
            process.kill()
            process.wait()


def test_build_refuses_a_template_a_serial_run_is_going_in(cal):
    write_run_lock(cal.template_workspace, method="run_ies", calibration="cal")
    with pytest.raises(RuntimeError, match=r"Refusing to start build\(\): the template"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("cal.pst", noptmax=0)
    assert (cal.template_workspace / "cal.pst").exists()
    (cal.template_workspace / RUN_LOCK).unlink()


@needs_true
def test_no_launch_starts_while_a_serial_run_is_going_in_the_template(cal, _no_results_reader):
    write_run_lock(cal.template_workspace, method="prior", calibration="cal")
    for launch in (lambda: cal.run_ies(exe="true"), lambda: cal.prior(workers=2, exe="true"),
                   lambda: cal.resume_ies()):
        with pytest.raises(RuntimeError, match="the template"):
            launch()
    (cal.template_workspace / RUN_LOCK).unlink()


# -- a serial run says what it deletes; a parallel master inherits nothing ----------------


def test_a_serial_run_deleting_earlier_results_says_so(cal, monkeypatch):
    commands = []
    monkeypatch.setattr(cal.pyemu.os_utils, "run", lambda command, cwd=None: commands.append(cwd))
    _results(cal.template_workspace, "cal", iterations=(0, 1))
    with pytest.warns(UserWarning, match=r"run_ies\(\) deleted 5 output file\(s\) of an earlier PESTPP-IES run"):
        cal.run_ies(reals=3)
    assert commands and not (cal.template_workspace / "cal.1.obs.csv").exists()
    assert not (cal.template_workspace / RUN_LOCK).exists()

    _results(cal.template_workspace, "cal", iterations=(0,))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cal.prior(reals=3, overwrite=True)   # asked for: no warning
    assert not (cal.template_workspace / "cal.0.obs.csv").exists()


@needs_true
def test_a_parallel_master_is_made_without_the_templates_serial_results(cal, monkeypatch, _no_results_reader):
    _results(cal.template_workspace, "cal", iterations=(0, 1, 2))
    seen = {}
    real = cal.pyemu.os_utils.start_workers

    def start_workers(*args, **kwargs):
        seen["master"] = sorted(path.name for path in Path(kwargs["master_dir"]).iterdir())
        return real(*args, **kwargs)

    monkeypatch.setattr(cal.pyemu.os_utils, "start_workers", start_workers)
    master = cal.run_ies(workers=2, exe="true")
    assert "cal.pst" in seen["master"] and RUN_LOCK in seen["master"]
    assert not [name for name in seen["master"] if name.startswith("cal.") and ".obs." in name]
    assert "cal.phi.actual.csv" not in seen["master"]
    assert (cal.template_workspace / "cal.2.obs.csv").exists(), "the template keeps its serial run"
    assert not (master / "cal.2.obs.csv").exists()


# -- resume: the splice --------------------------------------------------------------


def _phi_rows(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def test_a_resumed_runs_iterations_are_joined_after_the_interrupted_ones(tmp_path):
    case = "cal"
    run = tmp_path / "master"
    _results(run, case, iterations=(0, 1, 2))
    (run / f"{case}.obs+noise.csv").write_text("noise\n", encoding="utf-8")
    (run / f"{case}.rec").write_text("interrupted record\n", encoding="utf-8")
    (run / f"{case}.3.obs.csv").write_text("a partial iteration 3\n", encoding="utf-8")

    archive = archive_for_resume(run, case, 2)
    assert archive.name == "cal.resumed_from_2"
    assert not list(run.glob("cal.*.csv"))
    assert pending_resume(run, case) == archive

    # What PESTPP-IES writes, numbering from 0 again: 0 restates iteration 2.
    _results(run, case, iterations=(0, 1, 2))
    (run / f"{case}.obs+noise.csv").write_text("noise\n", encoding="utf-8")
    (run / f"{case}.rec").write_text("resumed record\n", encoding="utf-8")

    assert splice_resumed_run(run, case, archive) == [0, 1, 2, 3, 4]
    phi = _phi_rows(run / f"{case}.phi.actual.csv")
    assert phi["iteration"].tolist() == [0, 1, 2, 3, 4]
    assert phi["mean"].tolist() == [10.0, 9.0, 8.0, 9.0, 8.0]
    assert phi["total_runs"].tolist() == [3, 6, 9, 6 + 9, 9 + 9]
    for n, x in ((0, 0.0), (2, 2.0), (3, 1.0), (4, 2.0)):
        assert pd.read_csv(run / f"{case}.{n}.obs.csv", index_col=0)["x"].iloc[0] == x
    assert (archive / f"{case}.3.obs.csv").read_text() == "a partial iteration 3\n"
    assert (archive / f"{case}.rec").read_text() == "interrupted record\n"
    assert json.loads((archive / RESUME_RECORD).read_text())["spliced"] is True
    assert pending_resume(run, case) is None


def test_a_resume_that_wrote_nothing_puts_everything_back(tmp_path):
    run = tmp_path / "master"
    _results(run, "cal", iterations=(0, 1))
    before = sorted(path.name for path in run.iterdir())
    archive = archive_for_resume(run, "cal", 1)
    splice_resumed_run(run, "cal", archive)
    assert sorted(path.name for path in run.iterdir() if path != archive) == before
    assert _phi_rows(run / "cal.phi.actual.csv")["iteration"].tolist() == [0, 1]


# -- resume: refusals ---------------------------------------------------------------


def test_resume_refuses_what_it_cannot_resume(cal, monkeypatch):
    template = cal.template_workspace
    monkeypatch.setattr(cal.pyemu.os_utils, "run", lambda command, cwd=None: None)
    with pytest.raises(ValueError, match="sets ies_par_en itself"):
        cal.resume_ies(ies_par_en="mine.csv")
    with pytest.raises(ValueError, match="no run of 'cal' to resume") as nothing:
        cal.resume_ies(workers=2)
    assert "cal_ies_master" in str(nothing.value)

    cal.run_ies(reals=3, iterations=3)   # PESTPP-IES stubbed: writes nothing
    with pytest.raises(ValueError, match="completed no iteration"):
        cal.resume_ies()
    _results(template, "cal", iterations=(0, 1))
    with pytest.raises(ValueError, match=r"has completed iteration 1, and iterations=1 asks for no more"):
        cal.resume_ies(iterations=1)

    with pytest.warns(UserWarning, match="deleted 5 output file"):
        cal.prior(reals=3)
    _results(template, "cal", iterations=(0,))
    with pytest.raises(ValueError, match="is a prior Monte Carlo"):
        cal.resume_ies()


def test_resume_refuses_a_run_of_another_build(cal, monkeypatch):
    monkeypatch.setattr(cal.pyemu.os_utils, "run", lambda command, cwd=None: None)
    cal.run_ies(reals=3, iterations=3)
    _results(cal.template_workspace, "cal", iterations=(0, 1))
    cal.pst.parameter_data.loc[cal.pst.par_names[0], "parnme"] = "renamed"
    cal.pst.parameter_data.index = cal.pst.parameter_data["parnme"]
    with pytest.raises(ValueError, match="the template was rebuilt since"):
        cal.resume_ies()


REC = """\
...current lambda: 10
  ---  starting solve for iteration: 1  ---
...current lambda:  10
  ---  not updating lambda (standard deviation reduction criteria not met)  ---
  ---  starting solve for iteration: 2  ---
...current lambda:  10
  ---  updating lambda to  0.75  ---
  ---  EnsembleMethod iteration 2 report  ---
  ---  starting solve for iteration: 3  ---
...current lambda:  0.75
  ---  updating lambda to  0.05625  ---
"""


def test_resume_starts_from_the_lambda_the_run_had_reached(cal, monkeypatch):
    """Without it PESTPP-IES computes its initial lambda again, and the resumed run takes
    another path: on the canonical testing model the resumed iteration 2 read phi 1.651
    against the uninterrupted 1.714, and 1.71477 with it carried
    (``docs/reproducers/ies_resume_interrupted_run.py``)."""

    seen = {}

    def run(command, cwd=None):
        seen["pst"] = pyemu.Pst(str(Path(cwd) / "cal.pst"))

    monkeypatch.setattr(cal.pyemu.os_utils, "run", run)
    cal.run_ies(reals=3, iterations=4)   # stubbed: writes the control file only
    _results(cal.template_workspace, "cal", iterations=(0, 1, 2))
    (cal.template_workspace / "cal.rec").write_text(REC, encoding="utf-8")

    cal.resume_ies()
    options = seen["pst"].pestpp_options
    assert float(options["ies_initial_lambda"]) == 0.75
    assert options["ies_restart_obs_en"] == "cal.resumed_from_2/cal.2.obs.csv"
    assert options["ies_restart_par_en"] == "cal.resumed_from_2/cal.2.par.csv"
    assert options["ies_par_en"] == "cal.resumed_from_2/cal.0.par.csv"
    assert int(seen["pst"].control_data.noptmax) == 2
    assert int(cal.pst.control_data.noptmax) == 4, "cal.pst is not the resumed run's"


def _resume_stub(cal, monkeypatch):
    """PESTPP-IES stubbed: records each launch's NOPTMAX; the SECOND launch (the first
    resume) writes its restated iteration 0 and one iteration, then stops -- interrupted
    again."""

    launched = []

    def run(command, cwd=None):
        launched.append(int(pyemu.Pst(str(Path(cwd) / "cal.pst")).control_data.noptmax))
        if len(launched) == 2:
            _results(Path(cwd), "cal", iterations=(0, 1))

    monkeypatch.setattr(cal.pyemu.os_utils, "run", run)
    return launched


def test_a_resume_interrupted_again_resumes_to_the_runs_own_total(cal, monkeypatch):
    """The restart rewrites the run's control file with NOPTMAX the iterations LEFT; left
    so, a run of 4 cut after 1, resumed and cut again after 2 read as FINISHED -- the
    relaunch refusal said so and offered no resume -- and ``resume_ies()`` asked for 3 in
    all, then refused (review of ledger 533)."""

    launched = _resume_stub(cal, monkeypatch)
    template = cal.template_workspace
    cal.run_ies(reals=3, iterations=4)
    _results(template, "cal", iterations=(0, 1))
    cal.resume_ies()
    run_pst = pyemu.Pst(str(template / "cal.pst"))
    assert int(run_pst.control_data.noptmax) == 4
    for key in ("ies_par_en", "ies_restart_par_en", "ies_restart_obs_en", "ies_obs_en",
                "ies_initial_lambda"):
        assert key not in run_pst.pestpp_options, key

    cal.resume_ies()
    assert launched == [4, 3, 2]


def test_a_resume_killed_before_it_joined_is_joined_with_the_runs_total(cal, monkeypatch):
    """A hard kill skips the ``finally``: the next resume joins the archive AND puts the
    control file back before reading how many iterations the run was launched for."""

    import myflopy.modflow.mf6.pest.project as project

    launched = _resume_stub(cal, monkeypatch)
    template = cal.template_workspace
    cal.run_ies(reals=3, iterations=4)
    _results(template, "cal", iterations=(0, 1))
    with monkeypatch.context() as killed:
        killed.setattr(project, "splice_resumed_run", lambda *args: None)
        killed.setattr(project.PestProject, "_restore_resumed_control_file", lambda *args: None)
        cal.resume_ies()
    assert pending_resume(template, "cal") is not None

    with pytest.warns(UserWarning, match=r"joined the iterations a previous resume left unjoined"):
        cal.resume_ies()
    assert launched == [4, 3, 2]


def test_a_splice_stopped_part_way_is_finished_by_the_next(tmp_path, monkeypatch):
    """Run again after it stopped, a splice took the interrupted run's iterations, moved
    back by then, for the resumed run's: it deleted iteration 0 and renumbered the rest."""

    import myflopy.modflow.mf6.pest.run_state as run_state

    case, run = "cal", tmp_path / "master"
    _results(run, case, iterations=(0, 1, 2))
    archive = archive_for_resume(run, case, 2)
    _results(run, case, iterations=(0, 1, 2))   # the resumed run: 0 restates iteration 2

    def stopped(*args):
        raise OSError("stopped part-way")

    with monkeypatch.context() as broken:
        broken.setattr(run_state, "_join_phi", stopped)
        with pytest.raises(OSError, match="stopped part-way"):
            splice_resumed_run(run, case, archive)
    assert pending_resume(run, case) == archive

    assert splice_resumed_run(run, case, archive) == [0, 1, 2, 3, 4]
    for n, x in ((0, 0.0), (1, 1.0), (2, 2.0), (3, 1.0), (4, 2.0)):
        assert pd.read_csv(run / f"{case}.{n}.par.csv", index_col=0)["x"].iloc[0] == x, n
    assert _phi_rows(run / f"{case}.phi.actual.csv")["mean"].tolist() == [10.0, 9.0, 8.0, 9.0, 8.0]


def test_a_serial_run_deletes_an_earlier_resumes_archive(cal, monkeypatch):
    """A resume killed before it joined leaves the earlier run's iterations in an archive;
    a later serial run must not leave it for the next resume to splice into ITS run."""

    monkeypatch.setattr(cal.pyemu.os_utils, "run", lambda command, cwd=None: None)
    template = cal.template_workspace
    _results(template, "cal", iterations=(0, 1))
    archive = archive_for_resume(template, "cal", 1)
    with pytest.warns(UserWarning, match=r"cal\.resumed_from_1/"):
        cal.run_ies(reals=3)
    assert not archive.exists()
    assert pending_resume(template, "cal") is None


# -- resume: real PESTPP-IES -----------------------------------------------------------


@needs_pestpp
@pytest.mark.slow
def test_a_serial_run_cut_after_iteration_1_resumes_to_the_uninterrupted_runs_phi(toy, tmp_path):
    """A 3-realization run of 2 iterations, then the same run cut after iteration 1 (its
    iteration-2 files gone, as a killed kernel leaves them) and resumed: the resumed
    iteration 2 matches the uninterrupted one."""

    full = _calibration(toy, "cal", tmp_path / "full")
    full.run_ies(reals=3, iterations=2, **FAST)
    uninterrupted = pd.read_csv(full.template_workspace / "cal.phi.actual.csv")

    cut = tmp_path / "cut"
    shutil.copytree(full.template_workspace, cut)
    full._release_pyemu_log()
    for path in cut.glob("cal.2.*"):
        path.unlink()
    for name in ("cal.base.par", "cal.base.rei"):
        (cut / name).unlink(missing_ok=True)
    for phi in cut.glob("cal.phi.*.csv"):
        rows = pd.read_csv(phi)
        rows[rows["iteration"] <= 1].to_csv(phi, index=False)
    resumed_cal = toy.pest("cal", workspace=cut, start_datetime="2024-01-01")
    resumed_cal.pst = pyemu.Pst(str(cut / "cal.pst"))

    review = resumed_cal.resume_ies()
    assert review.iterations == [0, 1, 2]
    resumed = pd.read_csv(cut / "cal.phi.actual.csv")
    assert resumed["iteration"].tolist() == [0, 1, 2]
    assert resumed["mean"].iloc[:2].tolist() == uninterrupted["mean"].iloc[:2].tolist()
    assert resumed["mean"].iloc[2] == pytest.approx(uninterrupted["mean"].iloc[2], rel=1e-3)
    assert (cut / "cal.resumed_from_1" / "cal.rec").exists()
    assert not (cut / RUN_LOCK).exists()


@needs_pestpp
@pytest.mark.slow
def test_a_parallel_run_extended_by_resume_is_one_run_found_by_discovery(cal):
    first = cal.run_ies(reals=3, iterations=1, workers=2, **FAST)
    master = first.workspace
    assert first.iterations == [0, 1]
    assert not worker_root(master).exists()

    review = cal.resume_ies(workers=2, iterations=2)
    assert review.iterations == [0, 1, 2]
    phi = pd.read_csv(master / "cal.phi.actual.csv")
    assert phi["iteration"].tolist() == [0, 1, 2]
    assert phi["total_runs"].is_monotonic_increasing
    assert not worker_root(master).exists() and not (master / RUN_LOCK).exists()
    (handle,) = find_pest_runs(cal.template_workspace.parent)
    assert handle.executions == {"ies": master}
    assert handle.review("ies").iterations == [0, 1, 2]
