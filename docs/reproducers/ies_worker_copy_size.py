"""A parallel PESTPP-IES launch copies a template's earlier serial run into each worker -- measured on the canonical model (ledger 588).

pyEMU's ``start_workers`` copies the directory it is handed WHOLE into each worker and
takes no ignore list. ``run_ies(workers=N)`` handed it the template, so a calibration
first run serially (as a quick look, or a prior) carried that run's ensembles, phi files
and record into every worker. This measures the per-worker copy before and after.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/ies_worker_copy_size.py

Needs ``mf6``, ``pestpp-ies`` and pyEMU; writes only under a temporary directory, which
it removes.

THE SETUP. ``build_canonical_calibration_demo`` on the canonical TESTING profile
(``CanonicalModelConfig.testing()``, 21 x 21): ``k`` constant per layer and the 16 head
targets, as ``ies_resume_interrupted_run.py``. A SERIAL ``run_ies(reals=50,
iterations=2)`` leaves its outputs in the template. Then a parallel launch with ``true``
as PEST++ (it exits at once, so each worker holds exactly what was copied into it) and
``keep_workers=True``; the size of ``worker_0`` is the per-worker copy.

MEASURED (PEST++ 5.2.16, pyEMU 1.4.0, mf6 6.8.0.dev0, 2026-10-12)
  template as built                         0.72 MB in 157 files
  after the serial run_ies(50, 2)           3.36 MB in 210 files (+2.64 MB)
  per worker, BEFORE (the template whole)   3.36 MB in 210 files -- 53 files, 2.63 MB,
                                            not in the template as built: the run's
                                            ensembles, phi files, record, and the model
                                            outputs of its last forward run
  per worker, AFTER (what the build left)   0.74 MB in 157 files -- 0 such files
On Ten Trails a forward run writes a cell budget file per run, so the left-over model
output alone is a full .cbc per worker.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import shutil  # noqa: E402
import tempfile  # noqa: E402
import warnings  # noqa: E402

import myflopy.modflow.mf6.pest.ies as ies  # noqa: E402
from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo  # noqa: E402
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig  # noqa: E402
from myflopy.modflow.mf6.pest.run_state import worker_root  # noqa: E402

REALS, ITERATIONS = 50, 2
FAST = {"ies_lambda_mults": 1.0, "lambda_scale_fac": 1.0}


def _size(directory: pathlib.Path) -> tuple[int, int]:
    """(bytes, files) of everything under ``directory``."""

    files = [path for path in directory.rglob("*") if path.is_file()]
    return sum(path.stat().st_size for path in files), len(files)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ies_worker_copy_") as directory:
        root = pathlib.Path(directory)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            demo = build_canonical_calibration_demo(
                root / "model", config=CanonicalModelConfig.testing())
            cal = demo.model.pest("cal", start_datetime="2024-01-01")
            cal.parameterize("k", style="constant", bounds=(0.2, 5.0), physical=(1e-3, 1000.0))
            cal.observe(demo.head_targets)
            cal.build("cal.pst", noptmax=0)
        template = pathlib.Path(cal.template_workspace)
        built, built_files = _size(template)
        as_built = {path.relative_to(template) for path in template.rglob("*") if path.is_file()}
        print(f"template as built: {built / 1e6:.2f} MB in {built_files} files")

        cal.run_ies(reals=REALS, iterations=ITERATIONS, **FAST)
        after, after_files = _size(template)
        print(f"template after a serial run_ies(reals={REALS}, iterations={ITERATIONS}): "
              f"{after / 1e6:.2f} MB in {after_files} files "
              f"(+{(after - built) / 1e6:.2f} MB of that run's outputs)")

        # `true` as PEST++ wrote no results for the review to open.
        ies.IesResults = lambda workspace, **_: pathlib.Path(workspace)
        master = cal.run_ies(workers=2, exe=shutil.which("true"), keep_workers=True)
        first = worker_root(master) / "worker_0"
        worker, worker_files = _size(first)
        print(f"per worker, as copied: {worker / 1e6:.2f} MB in {worker_files} files")
        extra = [path for path in first.rglob("*")
                 if path.is_file() and path.relative_to(first) not in as_built]
        print(f"of which not in the template as built: {len(extra)} files, "
              f"{sum(path.stat().st_size for path in extra) / 1e6:.2f} MB: "
              f"{sorted(path.name for path in extra)[:12]}")


if __name__ == "__main__":
    main()
