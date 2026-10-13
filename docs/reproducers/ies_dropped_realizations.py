"""Which ensembles keep a realization PESTPP-IES dropped (ledger 575).

``IesResults.sensitivity()`` and ``parameters_at_bounds()`` report how many
realizations their statistics rest on, ``n_reals``, beside how many PESTPP-IES
started with, ``n_reals_drawn``, and warn below half. What each count reads is what
this script measures, on whichever tree is on ``PYTHONPATH``:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/ies_dropped_realizations.py

Needs ``mf6``, ``pestpp-ies`` and pyEMU; writes only under a temporary directory.

THE SETUP. The 6 x 6 Voronoi toy of ``tests/test_pest_handed_targets.py``; ``k`` and
``chd`` constant, three head targets; the template's ``forward_run.py`` made to exit 1
whenever the K multiplier PEST writes is below 1.0; ``run_ies(reals=12, iterations=2)``,
serial, ``ies_lambda_mults=1.0, lambda_scale_fac=1.0``.

MEASURED (PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-12)
  iteration 0  parameter ensemble 12 rows (every realization drawn, ``base`` included)
               observation ensemble 9 rows (the runs that finished)
  iteration 1  parameter and observation ensembles 5 rows each
  iteration 2  5 and 5
So the PRIOR parameter ensemble is what was drawn, the prior observation ensemble what
ran, and a later iteration's parameter ensemble what PESTPP-IES kept: ``n_reals`` 5 of
``n_reals_drawn`` 12 for the posterior, 9 of 12 for the prior.
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

HEADS = pd.DataFrame({"time": [0, 1, 2], "W2": [91.0, 92.0, 93.0],
                      "W1": [93.0, 93.5, 94.0], "W10": [95.0, 96.0, 97.0]})
LOCATIONS = pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0})


def main() -> None:
    with tempfile.TemporaryDirectory() as scratch:
        root = pathlib.Path(scratch)
        toy = run_toy(root / "toy")
        cal = toy.pest("d", workspace=root / "d", start_datetime="2024-01-01")
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
        cal.parameterize("chd", style="constant", bounds=(-1.0, 1.0))
        cal.observe(mf.HeadTargets(locations=LOCATIONS, values=HEADS))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("d.pst")
        template = cal.template_workspace
        script = template / "forward_run.py"
        script.write_text(
            "import sys as _sys\n"
            "if float(open('mult/k_inst0_constant.csv').read().split()[0]) < 1.0:\n"
            "    _sys.exit(1)\n" + script.read_text())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            review = cal.run_ies(reals=12, iterations=2, ies_lambda_mults=1.0,
                                 lambda_scale_fac=1.0)
        for n in review.iterations:
            par = pd.read_csv(template / f"d.{n}.par.csv", index_col=0)
            obs = pd.read_csv(template / f"d.{n}.obs.csv", index_col=0)
            print(f"iteration {n}: parameter ensemble {len(par)} rows, observation "
                  f"ensemble {len(obs)} rows")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            print(review.sensitivity()[["n_reals", "n_reals_drawn"]])
            print(review.parameters_at_bounds(which="prior")[["n_reals", "n_reals_drawn"]])
        cal._release_pyemu_log()


if __name__ == "__main__":
    main()
