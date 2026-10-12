"""An interrupted parallel PESTPP-IES run is refused, not deleted, and resumes to the uninterrupted run's phi (ledger 532).

A history match is launched with ``run_ies(workers=2)`` in a CHILD process -- standing in
for a Jupyter kernel -- and the child's whole process group (Python, the PESTPP-IES master
and both workers) is killed with SIGKILL once iteration 1 is complete: a killed kernel,
a reboot. Then, as a user does in a fresh kernel: the calibration is declared and built
again, a relaunch is attempted, and the run is resumed with ``resume_ies(workers=2)``.
The same calibration run without interruption is the reference.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/ies_resume_interrupted_run.py

Needs ``mf6``, ``pestpp-ies`` and pyEMU; writes only under a temporary directory, which
it removes. Runs at most 2 workers.

THE SETUP. ``build_canonical_calibration_demo`` on the canonical TESTING profile
(``CanonicalModelConfig.testing()``, 21 x 21): the canonical model as truth, 16 head wells
sampled from it, K perturbed x3 as the starting model. The calibration: ``k`` constant per
layer, bounds (0.2, 5); the 16 head targets; ``reals=6``, ``iterations=3``, one lambda
test per iteration (``ies_lambda_mults=1``, ``lambda_scale_fac=1``). PESTPP-IES draws
its own prior and noise from its fixed default seed, so the reference and the
interrupted run start from the same ensembles.

MEASURED (PEST++ 5.2.16, pyEMU 1.4.0, mf6 6.8.0.dev0, 2026-10-12, load average ~1-3)
  reference (uninterrupted)        71 s; phi mean by iteration 35.3446, 3.39681, 1.71448,
                                   1.35244 -- realization 3 FAILED in iteration 3 (its
                                   record: "par:obs realization runs failed ... 3:3"), so 5
                                   realizations there
  killed after iteration 1         the lock stays in the master; run_activity reads it
                                   STALE once the process group is gone; worker_0 and
                                   worker_1 stay (an unclean end keeps them)
  relaunch, run_ies(workers=2)     refused, FileExistsError: "holds an INTERRUPTED history
                                   match", naming cal.resume_ies(workers=2)
  resume_ies(workers=2)            46 s; iterations [0, 1, 2, 3]; phi 35.3446, 3.39681
                                   (the interrupted run's own rows), 1.71477, 1.13853 --
                                   iteration 2 within 1.7e-4 of the reference; iteration 3
                                   differs (16%) because the reference LOST realization 3
                                   there and the resumed run kept all 6; worker
                                   directories gone after the clean finish; the archive
                                   cut.resumed_from_1/ keeps cut.rec and cut.log
  the same, the lambda NOT carried iteration 2 read 1.65128 against 1.71448 (3.7%): the
  (before ies_initial_lambda)      restart computed its own initial lambda, so
                                   resume_ies sets ies_initial_lambda to the lambda the
                                   record shows in force (0.0075 here)
The 6-realization toy of tests/test_pest_handed_targets.py, cut after iteration 1 of 3 by
hand and restarted with the same options, gave 551.955 / 531.533 / 517.839 against the
uninterrupted 551.950 / 531.521 / 517.839, and `pestpp-ies case.pst /r` stops with
"/r option not supported by pestpp-ies".
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import os  # noqa: E402
import signal  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402

import pandas as pd  # noqa: E402

from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo  # noqa: E402
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig  # noqa: E402
from myflopy.modflow.mf6.pest.run_state import RUN_LOCK, run_activity, worker_root  # noqa: E402

REALS, ITERATIONS = 6, 3
FAST = {"ies_lambda_mults": 1.0, "lambda_scale_fac": 1.0}


def _demo(root: pathlib.Path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_canonical_calibration_demo(root / "model", config=CanonicalModelConfig.testing())


def _calibration(demo, name: str):
    cal = demo.model.pest(name, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 1000.0))
    cal.observe(demo.head_targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}.pst", noptmax=0)
    return cal


def _phi(directory: pathlib.Path, case: str) -> pd.DataFrame:
    """Each iteration's phi mean, run count, and how many realizations it holds."""

    phi = pd.read_csv(directory / f"{case}.phi.actual.csv")
    realizations = [c for c in phi.columns if c not in
                    ("iteration", "total_runs", "mean", "standard_deviation", "min", "max")]
    phi["realizations"] = phi[realizations].notna().sum(axis=1)
    return phi[["iteration", "total_runs", "realizations", "mean"]]


def _failures(directory: pathlib.Path, case: str) -> list[str]:
    """The record's lines about failed or dropped realizations."""

    text = (directory / f"{case}.rec").read_text(errors="replace")
    return [line.strip() for line in text.splitlines()
            if any(word in line.lower() for word in ("failed", "dropping", "overdue"))
            and "0 " not in line[:3]][:8]


def child(root: pathlib.Path) -> None:
    """The interrupted run: declared, built and launched as a notebook would."""

    cal = _calibration(_demo(root), "cut")
    cal.run_ies(reals=REALS, iterations=ITERATIONS, workers=2, **FAST)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ies_resume_") as directory:
        root = pathlib.Path(directory)
        demo = _demo(root)

        started = time.monotonic()
        reference = _calibration(demo, "ref").run_ies(
            reals=REALS, iterations=ITERATIONS, workers=2, **FAST)
        print(f"reference run: {time.monotonic() - started:.0f} s")
        print(_phi(pathlib.Path(reference.workspace), "ref").to_string(index=False))
        print("reference record, failed/dropped:", _failures(pathlib.Path(reference.workspace), "ref"))

        master = pathlib.Path(demo.model.workspace).with_name(
            pathlib.Path(demo.model.workspace).name + ".pest") / "cut_ies_master"
        process = subprocess.Popen(
            [sys.executable, __file__, "--child", str(root)],
            start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        deadline = time.monotonic() + 1800
        while time.monotonic() < deadline and process.poll() is None:
            phi = master / "cut.phi.actual.csv"
            try:
                if (master / "cut.1.obs.csv").is_file() and 1 in set(pd.read_csv(phi)["iteration"]):
                    break
            except (OSError, ValueError, KeyError):
                pass
            time.sleep(1)
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
        time.sleep(2)
        print(f"\nkilled the child's process group after iteration 1; master holds "
              f"{sorted(p.name for p in master.glob('cut.*.obs.csv'))}")
        print(f"lock left in the master: {(master / RUN_LOCK).is_file()}; "
              f"run_activity: {run_activity(master).state}; workers left: "
              f"{sorted(p.name for p in worker_root(master).glob('worker_*'))}")

        # A fresh kernel: declare and build again (the template is rebuilt; the master is not).
        cal = _calibration(demo, "cut")
        try:
            cal.run_ies(reals=REALS, iterations=ITERATIONS, workers=2, **FAST)
            print("relaunch: STARTED (the interrupted run's results were deleted)")
        except FileExistsError as refused:
            print("relaunch refused:", str(refused).splitlines()[0][:160], "...")
            print("   ", [line.strip() for line in str(refused).splitlines() if "resume_ies" in line])

        started = time.monotonic()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")   # the stale lock is said so
            resumed = cal.resume_ies(workers=2)
        print(f"\nresume_ies(workers=2): {time.monotonic() - started:.0f} s; "
              f"iterations {resumed.iterations}")
        joined = _phi(master, "cut")
        print(joined.to_string(index=False))
        print("resumed record, failed/dropped:", _failures(master, "cut"))
        reference_phi = _phi(pathlib.Path(reference.workspace), "ref")
        both = joined.merge(reference_phi, on="iteration", suffixes=("_resumed", "_reference"))
        both["relative difference"] = (both["mean_resumed"] / both["mean_reference"] - 1).abs()
        print("\n", both[["iteration", "mean_resumed", "mean_reference", "relative difference"]]
              .to_string(index=False))
        print(f"workers left after the clean finish: {worker_root(master).exists()}; "
              f"archive: {sorted(p.name for p in master.glob('cut.resumed_from_*/*'))}")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--child":
        child(pathlib.Path(sys.argv[2]))
    else:
        main()
