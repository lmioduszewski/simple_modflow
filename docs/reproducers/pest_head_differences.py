"""Difference observations score a datum error 0 where absolute heads score it fully.

Ledger 565. ``HeadTargets.differences()`` hands PEST each hydrograph's
shape: each reading minus its label's mean over its measured periods (``kind="mean"``),
or minus the previous measured one (``"previous"``). This script measures it through
PESTPP-IES.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_head_differences.py

Needs ``mf6`` and ``pestpp-ies``.

THE SETUP. ``tests/test_pest_handed_targets.py``'s toy (6 x 6 Voronoi, three 10-day
periods; heads move between periods as the stream's inflow changes). Targets at W1
(cell 1) and W2 (cell 2) are the model's OWN period-end heads plus a constant per well
-- +7 ft and -3 ft, a datum error -- registered twice: absolute (``hds``) and as mean
differences (``hdiff``); one parameter, ``k``; ``cal.prior(reals=3)``. Reported:
PESTPP-IES's ``phi.group.csv`` per realization, and the same with ``kind="previous"``.

MEASURED (mf6 6.8.0.dev0, PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-12)
  kind="mean"      base realization: hds 174.000 (= 3 x 7^2 + 3 x 3^2), hdiff 0.000000;
                   the two K-perturbed realizations: hds 173.949 / 174.288, hdiff
                   0.000119 / 0.003509 -- the shape moves with K, the datum never enters
  kind="previous"  base: hds 174.000, hdiff 0.000000; perturbed: hdiff 0.000076 /
                   0.002239
(the measured numbers are printed by ``main()``.)
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

# PEST++ starts each forward run in the template, where a relative PYTHONPATH would
# import another myflopy: made absolute here, before anything changes directory
# (ledger 425).
make_pythonpath_absolute()

import tempfile  # noqa: E402
import warnings  # noqa: E402

import pandas as pd  # noqa: E402

import myflopy as mf  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from tests.test_pest_handed_targets import run_toy  # noqa: E402


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="pest_head_differences_") as directory:
        root = pathlib.Path(directory)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = run_toy(root / "toy")
            probe = mf.HeadTargets(locations={"W1": 1, "W2": 2},
                                   values={"time": [0], "W1": [1.0], "W2": [1.0]})
            sim = probe.simulated_heads(model).set_index("per")
            targets = mf.HeadTargets(
                locations={"W1": 1, "W2": 2},
                values={"time": [0, 1, 2], "W1": (sim["W1"] + 7.0).tolist(),
                        "W2": (sim["W2"] - 3.0).tolist()},
            )
            for kind in ("mean", "previous"):
                cal = model.pest(f"d_{kind}", workspace=root / kind, start_datetime="2024-01-01")
                cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
                cal.observe(targets)
                cal.observe(targets.differences(kind))
                cal.build(f"d_{kind}.pst")
                ies = cal.prior(reals=3)
                table = pd.read_csv(ies.workspace / f"d_{kind}.phi.group.csv")
                print(f'kind="{kind}":')
                print(table.loc[:, ["par_realization", "hds", "hdiff"]].to_string(index=False))


if __name__ == "__main__":
    main()
