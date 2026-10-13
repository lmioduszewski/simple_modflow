"""The recalibration rehearsal's recovery, swept over ensemble size, iterations and seeds.

Ledger 589. ``canonical_08_recalibration_rehearsal.ipynb`` calibrates a
Ten Trails-shaped start model back to a truth that differs from it by known multipliers
(``myflopy.modflow.mf6.canonical_rehearsal``) and ends in a ``rehearsal_report``. This
script EXECUTES THAT NOTEBOOK -- its code cells, in order, one namespace -- once per
configuration and seed, with the first cell's knobs overridden, and prints each
report's numbers, so every threshold the slow test asserts rests on this table and the
workflow measured is the one the notebook teaches (not a copy of it).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/rehearsal_recovery_budget.py [config ...]

A config is ``REALSxITERATIONS`` (``16x2``, ``32x3``); default: every config below,
seeds 1-3. Needs ``mf6`` and ``pestpp-ies``; minutes per configuration (6 workers).

THE CRITERIA (the user's, Q9). Exact: the truth round trip (the template at the truth's
multipliers against the noise-free data, <= 1e-3), the write-back round trip (the model
rebuilt from ``ies.calibrated().to_library()`` against what PESTPP-IES ran for the base
realization), the equilibrium of the truth's last two cycles (0.05 ft). Stochastic:
measured phi of the last iteration <= 0.2 x the prior's; every scalar truth inside its
posterior 5-95% band; crisp's conductance closing >= 50% of its log gap; the forecast's
truth inside the posterior ensemble's range.

MEASURED (mf6 6.8.0.dev0, PESTPP-IES 5.2.16, pyEMU 1.4.0, 2026-10-12, 6 workers, nice 10,
load average 11-21 from other work; "round trip" and "write-back" are the worst error
over what is allowed it, <= 1 passes; every run's equilibrium passed):

     config seed   min phi ratio crisp gap forecast round trip write-back  outside the 5-95% band
       16x2    1   1.6    0.0138      1.01      yes    1.1e-05    2.1e-01  bedleak_north, bedleak_south, rch_dry, rch_wet
       16x2    2   1.6    0.0148      1.13      yes    1.1e-05    3.0e-01  bedleak_north, crisp, rch_dry, rch_wet, sy
       16x2    3   1.4    0.0201      1.24       NO    1.1e-05    1.6e-01  bedleak_north, bedleak_south, crisp, rch_dry, rch_wet
       32x3    1   3.3    0.0070      1.02       NO    1.1e-05    3.3e-01  bedleak_south
       32x3    2   3.3    0.0095      1.21      yes    1.1e-05    2.3e-01  bedleak_north, crisp, rch_wet
       32x3    3   3.7    0.0097      1.09       NO    1.1e-05    3.2e-01  crisp, rch_dry, sy
       48x3    1   6.7    0.0071      1.08      yes    1.1e-05    3.3e-01  bedleak_south, crisp, sy
       48x3    2   4.6    0.0075      1.03      yes    1.1e-05    3.3e-01  bedleak_north, rch_wet, sy
       48x3    3   6.6    0.0072      1.03      yes    1.1e-05    1.9e-01  bedleak_north, crisp, sy

The misses are far outside, not at a band's edge: 32x3 seed 2 put the north lake's bed
leakance at [1.58, 2.71] for a truth of 0.5 (z -7.7) and crisp at [3.63, 3.86] for 3.0
(z -9.4). The same three 32x3 runs with the notebook's ``cal.observe(heads.differences())``
line removed -- the head noise then counted once -- put 6, 5 and 5 of the 6 truths inside
(|z| <= 1.8 bar one rch_dry at -3.6, bands 3-10x wider), and their forecast ranges missed
the truth by 0.01-0.11 ft on all three seeds (phi ratio 0.008-0.014, crisp gap 0.92-0.98).
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
import tempfile
import time
import warnings

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

ROOT = pathlib.Path(__file__).resolve().parents[2]
NOTEBOOKS = ROOT / "examples" / "mf6" / "notebooks"
NOTEBOOK = NOTEBOOKS / "canonical_08_recalibration_rehearsal.ipynb"
CONFIGS = ("16x2", "24x3", "32x3", "48x3")
SEEDS = (1, 2, 3)


def run(reals: int, iterations: int, seed: int, workspace: pathlib.Path) -> dict:
    """Execute the notebook at these knobs; its ``rehearsal_report``."""

    import IPython.display

    cells = [
        "".join(cell["source"])
        for cell in json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
        if cell["cell_type"] == "code"
    ]
    IPython.display.display = lambda *objects, **options: None
    namespace = {"__name__": "__main__", "display": IPython.display.display}
    overrides = (f"REALS = {reals}\nITERATIONS = {iterations}\nSEED = {seed}\n"
                 f"WORKSPACE = Path({str(workspace)!r})\n")
    for index, source in enumerate(cells):
        exec(compile(source, f"<cell {index + 1}>", "exec"), namespace)
        if index == 0:
            exec(overrides, namespace)
    return namespace["rehearsal_report"]


def main(configs) -> None:
    os.chdir(NOTEBOOKS)
    sys.path.insert(0, str(NOTEBOOKS))
    print(f"{'config':>7} {'seed':>4} {'min':>5} {'phi ratio':>9} {'crisp gap':>9} "
          f"{'forecast':>8} {'round trip':>10} {'write-back':>10} {'equil':>5}  outside the 5-95% band")
    with tempfile.TemporaryDirectory(prefix="rehearsal_budget_") as scratch:
        for config in configs:
            reals, iterations = (int(part) for part in config.split("x"))
            for seed in SEEDS:
                start = time.time()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    report = run(reals, iterations, seed, pathlib.Path(scratch) / f"{config}_{seed}")
                passes = report["passes"]
                low, high = report["forecast_posterior"]
                outside = [f"{name} (z {report['scalars'][name]['z']:+.1f})"
                           for name, inside in report["scalars_inside"].items() if not inside]
                print(f"{config:>7} {seed:>4} {(time.time() - start) / 60:5.1f} "
                      f"{report['phi_ratio']:9.4f} {report['crisp_gap_closed']:9.2f} "
                      f"{'yes' if passes['forecast_covered'] else 'NO':>8} "
                      f"{report['truth_round_trip']:10.1e} {report['write_back']:10.1e} "
                      f"{'yes' if report['equilibrium'] else 'NO':>5}  {', '.join(outside) or '-'}"
                      f"  (forecast {report['forecast_truth']:.3f} in {low:.3f}..{high:.3f})",
                      flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or CONFIGS)
