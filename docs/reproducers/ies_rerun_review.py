"""A PESTPP-IES run's review reads that run, after another ran in the same place (ledgers 507, 508).

A serial ``run_ies``/``prior`` writes ``<name>.*`` into the template, where an earlier
run of the calibration left its own. The review took every ``<name>.<n>.obs.csv`` on
disk and called the highest the posterior; and ``model.pest_runs`` opened the stem of
the file ``cal.build(filename)`` wrote, not ``<name>``. This script drives the manual's
pattern through the public API on whichever tree is on ``PYTHONPATH``:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/ies_rerun_review.py

The "before" column is this script against the parent commit's ``src`` (``git archive
<commit> src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``). Needs ``mf6``,
``pestpp-ies`` and pyEMU; writes only under a temporary directory, which it removes.

THE SETUP. The 6 x 6 irregular Voronoi toy of ``tests/test_pest_handed_targets.py``,
run once; ``model.pest("valley_pest")`` + ``cal.build("valley.pst")``, ``k`` constant
and three head targets; all runs SERIAL in the template, three realizations,
``ies_lambda_mults=1.0, lambda_scale_fac=1.0``:

  1. ``run_ies(iterations=2)``
  2. ``run_ies(iterations=1)``
  3. ``prior()``
  4. ``model.pest_runs`` and its ``review()``

MEASURED (PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-12)   before (claude/pestprep)          after
  2. iterations / posterior / phi rows              [0, 1, 2] / 2 / [0, 1]            [0, 1] / 1 / [0, 1]
  3. iterations / posterior / phi rows              [0, 1, 2] / 2 / [0]               [0] / 0 / [0]
  4. the handle                                     built (not run), case 'valley'    prior, case 'valley_pest'
     review('ies')                                  case 'valley', iterations []      ValueError naming
                                                                                      review(kind='prior')
     review('prior')                                case 'valley', iterations []      iterations [0]
Step 1 reads [0, 1, 2] / 2 / [0, 1, 2] on both.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import tempfile  # noqa: E402
import warnings  # noqa: E402

import pandas as pd  # noqa: E402

import myflopy as mf  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from tests.test_pest_handed_targets import run_toy  # noqa: E402

FAST = {"ies_lambda_mults": 1.0, "lambda_scale_fac": 1.0}


def _show(label, results) -> None:
    print(f"{label}: iterations {results.iterations} posterior "
          f"{results.posterior_iteration} phi rows {list(results.phi['iteration'])}")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ies_rerun_") as directory:
        root = pathlib.Path(directory)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = run_toy(root)
            cal = model.pest("valley_pest", workspace=root / "pest" / "valley_pest",
                             start_datetime="2024-01-01")
            cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
            cal.observe(mf.HeadTargets(
                locations=pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34],
                                        "layer": 0}),
                values=pd.DataFrame({"time": [0, 1, 2], "W2": [91.0, 92.0, 93.0],
                                     "W1": [93.0, 93.5, 94.0], "W10": [95.0, 96.0, 97.0]}),
            ))
            cal.build("valley.pst")
        _show("1. run_ies(iterations=2)", cal.run_ies(reals=3, iterations=2, **FAST))
        _show("2. run_ies(iterations=1)", cal.run_ies(reals=3, iterations=1, **FAST))
        _show("3. prior()", cal.prior(reals=3))
        from myflopy.modflow.mf6.pest.runs import find_pest_runs

        handle = find_pest_runs(root / "pest")[0]
        print(f"4. handle {handle!r}, case {handle.case!r}")
        for kind in ("ies", "prior"):
            try:
                review = handle.review(kind)
                print(f"   review({kind!r}): case {review.case!r}, iterations {review.iterations}")
            except (ValueError, FileNotFoundError) as error:
                print(f"   review({kind!r}): {type(error).__name__}: {str(error)[:120]}")


if __name__ == "__main__":
    main()
