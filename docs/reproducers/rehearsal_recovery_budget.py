"""The recalibration rehearsal's recovery, swept over ensemble size, iterations and seeds.

Ledger 589. ``canonical_08_recalibration_rehearsal.ipynb`` calibrates a
Ten Trails-shaped start model back to a truth that differs from it by known multipliers
(``myflopy.modflow.mf6.canonical_rehearsal``) and ends in a ``rehearsal_report``. This
script EXECUTES THAT NOTEBOOK -- its code cells, in order, one namespace -- once per
configuration and seed, with the first cell's knobs overridden, and prints each
report's numbers, so every threshold the slow test asserts rests on this table and the
workflow measured is the one the notebook teaches (not a copy of it).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/rehearsal_recovery_budget.py \
        [config ...] [--seeds 1,2,3] [--workers N] [--keep DIR] [--linear] [--reports DIR]

A config is ``REALSxITERATIONS[+flag...]`` (``16x2``, ``64x6+lambda``); default: the
four unflagged configs of the first table, seeds 1-3. ``+lambda`` turns PESTPP-IES's own
lambda testing back on (drops the notebook's ``ies_lambda_mults=1.0,
lambda_scale_fac=1.0``, so its default multipliers and scale factors are tried every
iteration); ``+nodiff`` observes the levels alone at the reading's weight (no
``heads.differences``); ``+double`` restores the notebook's declaration before ledger
593 -- the levels and their differences both at the reading's weight. Each
flag edits the notebook's SOURCE before it runs and refuses to run unless each of its
texts is found exactly once. ``--workers`` overrides the notebook's ``WORKERS``;
``--keep`` keeps every run's workspace and its ``report.json`` (a relative path is
resolved where the script is started); ``--linear`` adds :func:`linear_check` after each
run; ``--reports DIR`` runs nothing and pools the ``report.json`` files a ``--keep`` run
left. After the per-run table it prints, for each configuration, the seeds POOLED:
``canonical_rehearsal.pooled_recovery``'s verdict and how often an exact posterior fails
its thresholds -- by the binomial it uses, and with each run's z correlated as its
posterior is (:func:`correlated_false_fail`). Needs ``mf6`` and ``pestpp-ies``; minutes
per run.

THE CRITERIA (the user's, Q9; the coverage restated 2026-10-14, ledger 593).
Exact: the truth round trip (the template at the truth's multipliers against the
noise-free data, <= 1e-3), the write-back round trip (the model rebuilt from
``ies.calibrated().to_library()`` against what PESTPP-IES ran for the base realization),
the equilibrium of the truth's last two cycles (0.05 ft). Stochastic, per run: measured
phi of the last iteration <= 0.2 x the prior's; crisp's conductance closing >= 50% of its
log gap. Stochastic, POOLED over several converged runs (``pooled_recovery``): the
fraction of scalar truths inside their posterior 5-95% bands, the largest |z| and the
runs whose forecast range holds its truth -- six runs: at least 28 of 36, at most 3.63,
at least 5 of 6. That replaced "every scalar truth inside its band, and the forecast
covered, in every run", which an exact posterior meets only 0.9^6 = 53% of the time.

THE BANDS. Since ledger 593 a band is the normal fitted to the posterior
REALIZATIONS' log10 values (mean +- 1.645 sd; the ``base`` row left out) and ``inside``
is ``|z| <= 1.645``. Every table above the POOLED one below was scored by the earlier
rule -- the ensemble's own 5th/95th percentiles with the ``base`` row in them, and ``z``
from the same moments -- under which the two could disagree (64x6 seed 6: sy z -1.56
reported outside). And every table above it ran the notebook's declaration of then: its
unflagged rows are what ``+double`` is now.

"phi truth" is the phi of the TRUTH's multipliers on the noisy data the run fitted (the
template run at them), "phi base" the posterior base realization's: when the truth fits
better, the posterior has not converged (ledger 590).

MEASURED WITH NO NOISE DRAWN -- before ledger 590 ``run_ies``'s default
left PESTPP-IES's noise off (mf6 6.8.0.dev0, PESTPP-IES 5.2.16, pyEMU 1.4.0, 2026-10-12,
6 workers, nice 10, load average 11-21 from other work; "round trip" and "write-back"
are the worst error over what is allowed it, <= 1 passes; every run's equilibrium
passed; the phi columns were not printed then):

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
line removed put 6, 5 and 5 of the 6 truths inside (|z| <= 1.8 bar one rch_dry at -3.6,
bands 3-10x wider), and their forecast ranges missed the truth by 0.01-0.11 ft on all
three seeds (phi ratio 0.008-0.014, crisp gap 0.92-0.98). It is not the heads' noise
counted twice: 32x3 seed 2 with BOTH head weights divided by sqrt(2) left bedleak_north
(z -6.4), crisp (-8.3) and rch_wet (-4.3) outside, as at full weight.

MEASURED WITH NOISE DRAWN (after ledger 590; 2 workers, load ~4; the
notebook's own knobs bar REALS/ITERATIONS/SEED):

     config seed  phi ratio  crisp gap  forecast  round trip  write-back  phi truth  phi base  outside
       16x2 2026   <= 0.2     >= 0.5     NO        <= 1        <= 1          352.3    1043.6  all 6
       16x2    1   <= 0.2     >= 0.5     yes       <= 1        <= 1          349.7     637.7  5 (crisp in)
       16x2    3   <= 0.2     >= 0.5     NO        <= 1        <= 1          386.3    1048.9  5 (sy in)
       32x3    2   <= 0.2     >= 0.5     NO        <= 1        <= 1          309.6     503.2  bedleak_north, crisp, rch_wet

Every exact gate, phi and crisp's gap pass with noise as without; the posterior base fits
the data 1.6-3.0x worse than the truth does in every run, so the coverage misses are an
ensemble IES has not converged (``ies_lambda_mults=1.0``, ``lambda_scale_fac=1.0``, two or
three iterations), with the realizations collapsed around where it stopped.

MEASURED AT A CONVERGING BUDGET (ledger 591; noise drawn; 12 workers, nice
10, load average 10-17 from other work; 64x6: 736 model runs, 8.5-8.8 min a run; 64x3:
400, 4.8-5.2 min; 32x4: 288, 3.7 min). Every run: truth round trip 1.1e-05, write-back
0.13-0.36, equilibrium yes, phi ratio 0.007-0.014, crisp gap 0.97-1.02. "base/truth" is
phi base over phi truth; below 1, the base fits the noisy data better than the truth
does, as a converged fit must:

     config            seed  phi (measured+noise) by iteration           base/truth  forecast  outside the 5-95% band (z)
     64x6+lambda          1  52488 1731 674 635 629 627 625                    0.92  yes       rch_wet -1.7
     64x6+lambda          2  53106 1785 657 619 613 610 608                    0.95  yes       sy -1.8
     64x6+lambda          3  54191 1988 706 667 661 659 658                    0.89  yes       rch_wet +2.5, sy +4.0
     64x6+lambda          4  48538 2282 658 600 591 588 587                    0.95  yes       rch_dry +2.0
     64x6+lambda          5  47940 1665 655 604 595 591 588                    0.88  yes       -
     64x6+lambda          6  47655 1500 699 663 656 654 652                    0.90  yes       bedleak_north -2.5, rch_dry +2.9, sy -1.6
     64x6+lambda+nodiff   1  50406 1447 401 363 357 355 354                    0.88  yes       rch_wet -2.2
     64x6+lambda+nodiff   2  51073 1433 383 350 343 341 339                    0.92  yes       sy -1.7
     64x6+lambda+nodiff   3  52162 1602 406 366 360 358 357                    0.84  yes       rch_wet +2.0, sy +2.8
     64x6+lambda+nodiff   4  46295 1935 397 342 332 328 326                    0.92  yes       -
     64x6+lambda+nodiff   5  46136 1340 385 335 327 324 322                    0.82  yes       -
     64x6+lambda+nodiff   6  45585 1193 406 372 366 364 362                    0.87  yes       bedleak_north -2.5, rch_dry +2.1
     64x3+lambda          1  52488 1731 674 635                                0.93  yes       rch_wet -1.8
     64x3+lambda          2  53106 1785 657 619                                0.96  yes       -
     64x3+lambda          3  54191 1988 706 667                                0.90  yes       rch_wet +1.4, sy +3.8
     32x4+lambda          1  50659 1419 669 644 640                            0.93  NO        rch_wet -2.4
     32x4+lambda          2  51819 2226 654 614 608                            0.96  yes       bedleak_north -1.6, sy -2.3
     32x4+lambda          3  47307 2171 762 676 663                            0.90  yes       rch_wet +1.6, sy +3.7
     16x2+lambda          1  52671 1474 803                                    1.43  NO        bedleak_north +1.7, bedleak_south +10.6, rch_wet -3.2, sy -2.1

Pooled, the six-seed runs put 28 of 36 scalar truths inside their 5-95% bands with the
differences (RMS z 1.37) and 30 of 36 without (RMS z 1.19); 1 and 2 runs of 6 had all six
inside, where an exact posterior's would about half the time (0.9^6). The forecast's
truth was inside its posterior range in 12 of 12. Phi has nearly stopped falling after the
third iteration (the fourth to sixth lower it 1.4-4.6% in all, the sixth under 0.7%), and
64x3 recovers what 64x6 does; two iterations, lambda testing or not, do not converge.

``--linear`` (``linear_check`` on five 64x6 runs: seeds 1, 3 and 6, seeds 3 and 6
nodiff): the ensemble's posterior sd against the linear posterior's at the base is
0.98-1.16x for sy, crisp and the bed leakances, 0.82-1.0x for rch_wet (0.82-0.83x levels
only) and 0.71-0.91x for rch_dry -- the ensemble is NARROWER than the data say for the
recharge multipliers, in all five runs. Of the four misses the levels alone leave on
seeds 3 and 6, three are the linear posterior's own (|z_lin| > 1.645: seed 3 sy +2.97 and
rch_wet +1.80, seed 6 bedleak_north -2.34), so the noise drawn; the fourth, seed 6
rch_dry (z +2.08 against z_lin +1.40), is that narrowness. The sd of z under the data's
actual noise is 1.15-1.33 for rch_wet, rch_dry and sy with the differences beside the
levels (each reading's noise counted twice in the shape) and 0.82-0.98 for every parameter
with the levels alone; crisp and the bed leakances 0.92-1.02 either way (ledger
592).

REVIEWED (ledger 591; 2026-10-13, load 7-16): ``64x3+lambda --seeds 2``
reran to the digit (phi 53106 1785 657 619, base/truth 0.96, all six inside), and
``128x3+lambda+nodiff --seeds 6 --linear`` (800 model runs, 13.6 min at 8 workers), the
seed whose rch_dry the 64-member ensemble missed:

     config               seed  phi by iteration          base/truth  forecast  outside (z; z_lin)
     128x3+lambda+nodiff     6  52051 1204 405 370              0.90  yes       bedleak_north -1.8 (-2.0), rch_wet -1.6 (-1.5)

rch_dry is inside (z +0.78) and its sd 0.85x the linear's (0.71x at 64x6); rch_wet 0.89x
(0.82x), bedleak_south 1.29x. rch_wet is outside by 0.001 of its band at |z| < 1.645 --
the band is the ensemble's percentiles, z its log mean and sd, and they need not agree.

POOLED (ledger 593; 2026-10-14, noise drawn, 12 workers, nice 10, load average
9-18 from other work; 4.7-4.8 min a run). The scoring of every row here is the one rule
above; "z" in the outside column. Every run: truth round trip 1.1e-05, write-back
0.19-0.30, equilibrium yes, phi ratio 0.012-0.015, crisp gap 0.96-1.01, forecast covered:

     config               seed  base/truth  outside the 5-95% band (z)
     64x3+lambda+double      1        0.93  rch_wet -1.79
     64x3+lambda+double      2        0.96  -
     64x3+lambda+double      3        0.90  sy +3.73
     64x3+lambda+double      4        0.98  -
     64x3+lambda+double      5        0.92  -
     64x3+lambda+double      6        0.92  bedleak_north -2.29, sy -1.69
     64x3+lambda (split)     1        0.92  -
     64x3+lambda (split)     2        0.96  -
     64x3+lambda (split)     3        0.90  sy +2.56
     64x3+lambda (split)     4        0.99  bedleak_north -1.73
     64x3+lambda (split)     5        0.92  -
     64x3+lambda (split)     6        0.93  bedleak_north -2.13

(the ``+double`` rows were measured as the notebook then stood, unflagged, and the split
ones as ``+split`` -- an edit this script carried until the notebook took it.)

     pooled              inside (>= 28)  largest |z| (<= 3.63)  forecast (>= 5)  verdict
     64x3+lambda+double  32 of 36        3.73 (seed 3, sy)      6 of 6           FAILS
     64x3+lambda (split) 33 of 36        2.56 (seed 3, sy)      6 of 6           PASSES

An exact posterior fails these thresholds 3.1% of the time by the binomial and the
independent tails ``pooled_recovery`` computes (coverage 0.77%, |z| 1.0%, forecast 1.35%
at 63 realizations), and about 4% with each run's z correlated as its posterior is --
the rehearsal's carries rch_wet with sy at 0.68-0.71 and with rch_dry at 0.61-0.62,
which takes the coverage threshold's rate from 0.77% to 1.8% (Monte Carlo, both
measurements alike); 29, the binomial's 2.5% point, would have failed an exact
posterior 4.1% of the time on its own. The levels beside their differences at the
reading's weight fail on the largest |z|: the same seed's sy is +3.73 with the shape
counted twice and +2.56 with it counted once (ledger 592: that double count widens a
seasonal parameter's z by 1.15-1.33x). The split -- the levels at ``1 / hypot(0.05,
DATUM=0.05)`` ft, the differences at ``sqrt(0.05**-2 - w_levels**2)`` -- is the
notebook's declaration since. The synthetic data carry no datum error, so it
under-weights the mean levels: a posterior on the safe side. The notebook's CI budget
with it (16x2, seeds 2026, 1, 2): phi ratio 0.021-0.029, crisp gap 1.01-1.16, round trip
1.1e-05, write-back 0.22-0.32 -- its per-run gates hold; the pooled verdict fails (3 of
18 inside, |z| up to 14.6), as an unconverged ensemble must.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
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
#: What a config's ``+flag`` does to the notebook's SOURCE before it runs: one or more
#: (old, new) edits. Each must match its cell exactly once, or the run is refused: a
#: notebook edited since would otherwise be measured without the change the table claims.
FLAGS = {
    # PESTPP-IES's own lambda testing (its default ies_lambda_mults and
    # lambda_scale_fac): drop the notebook's single-lambda overrides.
    "lambda": (("ies_lambda_mults=1.0, lambda_scale_fac=1.0", ""),),
    # The levels alone, at the reading's weight 1 / NOISE["head"]: no shape targets.
    "nodiff": (
        ('    w_levels = 1.0 / np.hypot(NOISE["head"], DATUM)\n',
         '    w_levels = 1.0 / NOISE["head"]\n'),
        ('    w_shape = np.sqrt(NOISE["head"] ** -2 - w_levels**2)\n', ""),
        ("    cal.observe(heads.differences(weight=w_shape))\n", ""),
    ),
    # The notebook's declaration before ledger 593: the levels AND their
    # differences both at the reading's weight, each reading's noise counted twice in
    # the shape (ledger 592).
    "double": (
        ("case.wells.assign(weight=w_levels)", "case.wells"),
        ("    cal.observe(heads.differences(weight=w_shape))\n",
         "    cal.observe(heads.differences())\n"),
    ),
}


def parse_config(config: str) -> tuple[int, int, tuple[str, ...]]:
    """``REALSxITERATIONS[+flag...]`` -> (reals, iterations, flags)."""

    match = re.fullmatch(r"(\d+)x(\d+)((?:\+\w+)*)", config)
    if match is None:
        raise SystemExit(f"config {config!r}: write REALSxITERATIONS[+flag...], e.g. 64x6+lambda")
    flags = tuple(flag for flag in match.group(3).split("+") if flag)
    unknown = sorted(set(flags) - set(FLAGS))
    if unknown:
        raise SystemExit(f"config {config!r}: unknown flag(s) {unknown}; known: {sorted(FLAGS)}")
    return int(match.group(1)), int(match.group(2)), flags


def _cells(flags) -> list[str]:
    cells = [
        "".join(cell["source"])
        for cell in json.loads(NOTEBOOK.read_text(encoding="utf-8"))["cells"]
        if cell["cell_type"] == "code"
    ]
    for flag in flags:
        for old, new in FLAGS[flag]:
            hits = [index for index, source in enumerate(cells) if old in source]
            if len(hits) != 1 or cells[hits[0]].count(old) != 1:
                raise SystemExit(f"+{flag}: {old!r} is not in exactly one place of "
                                 f"{NOTEBOOK.name}; update FLAGS to the notebook's current source")
            cells[hits[0]] = cells[hits[0]].replace(old, new)
    return cells


def run(reals: int, iterations: int, seed: int, workspace: pathlib.Path, *,
        flags=(), workers: int | None = None, linear: bool = False) -> dict:
    """Execute the notebook at these knobs; its ``rehearsal_report``, plus phi detail.

    ``linear=True`` adds ``report["linear"]``: :func:`linear_check` at the posterior base.
    """

    import IPython.display

    IPython.display.display = lambda *objects, **options: None
    namespace = {"__name__": "__main__", "display": IPython.display.display}
    overrides = (f"REALS = {reals}\nITERATIONS = {iterations}\nSEED = {seed}\n"
                 f"WORKSPACE = Path({str(workspace)!r})\n")
    if workers is not None:
        overrides += f"WORKERS = {int(workers)}\n"
    for index, source in enumerate(_cells(flags)):
        exec(compile(source, f"<cell {index + 1}>", "exec"), namespace)
        if index == 0:
            exec(overrides, namespace)
    report = dict(namespace["rehearsal_report"])
    report["phi_truth"], report["phi_base"] = _phi_at_truth_and_base(namespace, workspace)
    phi = namespace["ies"].phi_measured.set_index("iteration")
    report["phi_by_iteration"] = {int(it): {"mean": float(row["mean"]), "base": float(row["base"])
                                            if "base" in row else float("nan"),
                                            "runs": int(row["total_runs"])}
                                  for it, row in phi.iterrows()}
    report["scalar_correlation"] = scalar_correlation(namespace["cal"], namespace["ies"],
                                                      report["scalars"])
    if linear:
        report["linear"] = linear_check(namespace["cal"], namespace["ies"].workspace,
                                        namespace["ies"].case, report["scalars"], workspace / "linear",
                                        workers=workers or 4,
                                        head_noise=namespace["NOISE"]["head"])
    return report


def _weighted_phi(cal, simulated):
    """PEST++'s measured phi of ``simulated`` against ``cal``'s data (a bound holding is 0)."""

    from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal

    table = rehearsal.round_trip(cal, simulated)
    weights = cal.pst.observation_data["weight"].astype(float)
    weights.index = weights.index.str.lower()
    return float(((table["error"] * weights.loc[table.index]) ** 2).sum())


def _phi_at_truth_and_base(namespace, workspace):
    """Phi on the NOISY data at the truth's multipliers, and at the posterior base.

    The test of whether the posterior is a converged answer (ledger 590):
    when the truth fits the data better than the base realization IES ended on, the
    narrow bands are an ensemble that has not converged, not what the data say.
    """

    from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal

    cal = namespace["cal"]
    truth = rehearsal.truth_parameter_values(namespace["case"], cal)
    at_truth = rehearsal.run_template(cal, truth, workspace / "truth_on_noisy_data")
    return _weighted_phi(cal, at_truth), _weighted_phi(cal, rehearsal.base_observations(namespace["ies"]))


def scalar_correlation(cal, ies, scalars) -> dict:
    """The posterior correlation of the scalars' log10 values, ``{group: {group: r}}``.

    From the posterior REALIZATIONS (the ``base`` row left out). Under an exact
    posterior a run's z values are correlated as its parameters are, which the
    binomial behind :func:`~myflopy.modflow.mf6.canonical_rehearsal.pooled_recovery`'s
    coverage threshold ignores; :func:`correlated_false_fail` puts it back.
    """

    import numpy as np

    data = cal.pst.parameter_data
    groups = data.groupby("pargp")["parnme"].apply(list)
    names = {group: str(members[0]).lower() for group, members in groups.items()
             if len(members) == 1 and group in scalars}
    posterior = ies.posterior_parameters._df.rename(columns=str.lower)
    posterior = posterior.drop(index=[i for i in posterior.index if str(i).lower() == "base"])
    logs = np.log10(posterior[list(names.values())].astype(float))
    logs.columns = list(names)
    return logs.corr().to_dict()


def correlated_false_fail(reports, verdict, draws=200_000, seed=0) -> dict:
    """How often an exact posterior fails ``verdict``'s thresholds with each run's z
    correlated as its posterior is (``report["scalar_correlation"]``): Monte Carlo of
    one multivariate normal per run, the forecast at its exact binomial rate."""

    import numpy as np

    from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal

    rng = np.random.default_rng(seed)
    inside = np.zeros(draws, int)
    largest = np.zeros(draws)
    for report in reports:
        order = list(report["scalars"])
        corr = np.array([[report["scalar_correlation"][a][b] for b in order] for a in order])
        z = rng.multivariate_normal(np.zeros(len(order)), corr, size=draws, method="eigh")
        inside += (np.abs(z) <= rehearsal.Z_BAND).sum(axis=1)
        largest = np.maximum(largest, np.abs(z).max(axis=1))
    coverage = float((inside < verdict.inside_min).mean())
    largest_z = float((largest > verdict.z_max).mean())
    return {"coverage": coverage, "largest_z": largest_z,
            "either": float(((inside < verdict.inside_min) | (largest > verdict.z_max)).mean())}


LINEAR_STEP = 0.01  # log10 step of the finite-difference Jacobian (2.3 % of a value)


def linear_check(cal, master, case, scalars, scratch, *, workers=4, head_noise=None) -> dict:
    """The linear (first-order) posterior at the base realization, beside the ensemble's.

    A finite-difference Jacobian of every measured equality observation (the bound
    left out: it holds) against every adjustable parameter in log10, at the posterior
    ``base`` realization -- one template run per parameter (:func:`run_template`) --
    and a prior covariance from the prior ensemble IES drew (``<case>.0.par.csv``;
    the scalars' cross terms set to 0, as drawn). For each scalar parameter group:

    - ``ies_sd``/``lin_sd``: the posterior sd (log10) of the ensemble and of the
      linear posterior ``(J'QJ + Cp^-1)^-1`` -- equal when IES has converged to what
      the data say;
    - ``z_lin``: the truth's distance from the base in linear-posterior sds;
    - ``z_sd`` / ``z_sd_levels``: the sd of that z over repeated draws of the data's
      ACTUAL noise -- independent per (well, month) level, lake stage and flow, and
      each ``hdiff`` difference the SAME draws through its well's mean
      (``synthetic_hydrographs`` bins to truth + draw, and the differences are taken
      from those binned targets) -- with the observations as declared and with the
      levels alone. 1 is a calibrated posterior; above 1 its bands are that many times
      too narrow. A level's noise is ``head_noise`` (ft) when given -- the notebook's
      ``NOISE["head"]``, which the levels' weight no longer is once it carries a datum
      error (ledger 593) -- else 1 / its weight, as every other observation's.
    """

    from concurrent.futures import ThreadPoolExecutor

    import numpy as np
    import pandas as pd

    from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal

    pst = cal.pst
    master = pathlib.Path(master)
    iterations = sorted(int(path.name.split(".")[-3]) for path in master.glob(f"{case}.*.par.csv"))
    posterior = pd.read_csv(master / f"{case}.{iterations[-1]}.par.csv", index_col=0)
    prior = pd.read_csv(master / f"{case}.0.par.csv", index_col=0)
    posterior.columns = posterior.columns.str.lower()
    prior.columns = prior.columns.str.lower()
    data = pst.parameter_data.copy()
    data.index = data.index.str.lower()
    names = [name for name in data.index if data.loc[name, "partrans"] not in ("fixed", "tied")]
    group_of = data.loc[names, "pargp"]
    single = {group: list(members) for group, members in group_of.groupby(group_of).groups.items()
              if len(members) == 1}
    scalar_names = [single[group][0] for group in scalars if group in single]
    base = posterior.loc["base", names].astype(float)
    prior_cov = np.log10(prior.drop(index="base", errors="ignore")[names].astype(float)).cov().to_numpy()
    for name in scalar_names:
        i = names.index(name)
        keep = prior_cov[i, i]
        prior_cov[i, :] = prior_cov[:, i] = 0.0
        prior_cov[i, i] = keep

    observations = pst.observation_data.copy()
    observations.index = observations.index.str.lower()
    weight = observations["weight"].astype(float)
    group = observations["obgnme"].astype(str)
    used = weight.index[(weight > 0) & ~group.str.startswith(("l_", "less", "g_", "greater"))]

    jobs = {"base": base.to_dict()}
    for name in names:
        values = base.to_dict()
        values[name] = base[name] * 10 ** LINEAR_STEP
        jobs[name] = values
    scratch = pathlib.Path(scratch)

    def forward(item):
        key, values = item
        return key, rehearsal.run_template(cal, values, scratch / re.sub(r"\W", "_", key))

    with ThreadPoolExecutor(max(1, int(workers))) as pool:
        simulated = dict(pool.map(forward, jobs.items()))
    jacobian = np.column_stack([
        (simulated[name].loc[used] - simulated["base"].loc[used]).to_numpy() / LINEAR_STEP
        for name in names])
    w = weight.loc[used].to_numpy()
    groups = group.loc[used].to_numpy()

    # The data's actual noise as a map of independent unit draws: noise = A @ draws.
    n = len(used)
    sd = 1.0 / np.where(w > 0, w, 1.0)
    if head_noise is not None:
        sd = np.where(groups == "hds", float(head_noise), sd)
    noise_map = np.diag(np.where(groups == "hdiff", 0.0, sd))
    level_row = {re.sub(r"^oname:hds_", "", name): i for i, name in enumerate(used) if groups[i] == "hds"}
    by_well = {}
    for i, name in enumerate(used):
        if groups[i] == "hdiff":
            source = level_row[re.sub(r"^oname:hdiff_", "", name)]
            by_well.setdefault(re.search(r"usecol:(.+?)_", name).group(1), []).append((i, source))
    for rows in by_well.values():
        mean = sum(noise_map[source] for _, source in rows) / len(rows)
        for i, source in rows:
            noise_map[i] = noise_map[source] - mean
    noise_cov = noise_map @ noise_map.T

    def posterior_and_z_sd(mask):
        J, Q = jacobian[mask], np.diag(w[mask] ** 2)
        cov = np.linalg.inv(J.T @ Q @ J + np.linalg.inv(prior_cov))
        sampling = cov @ J.T @ Q @ noise_cov[np.ix_(mask, mask)] @ Q @ J @ cov
        return cov, np.sqrt(np.diag(sampling) / np.diag(cov))

    declared, z_sd = posterior_and_z_sd(np.ones(n, bool))
    _, z_sd_levels = posterior_and_z_sd(groups != "hdiff")
    log_post = np.log10(posterior.drop(index="base", errors="ignore")[names].astype(float))
    result = {}
    for group_name in scalars:
        if group_name not in single:
            continue
        name = single[group_name][0]
        i = names.index(name)
        truth = np.log10(scalars[group_name]["truth"])
        lin_sd = float(np.sqrt(declared[i, i]))
        result[group_name] = {
            "ies_sd": float(log_post[name].std()), "lin_sd": lin_sd,
            "z_lin": float((truth - np.log10(base[name])) / lin_sd),
            "z_sd": float(z_sd[i]), "z_sd_levels": float(z_sd_levels[i]),
        }
    return result


def _print_run(label: str, seed: int, minutes: float, report: dict) -> None:
    """One run's numbers, in full: what the summary row compresses."""

    phis = report["phi_by_iteration"]
    print(f"\n== {label} seed {seed}: {minutes:.1f} min", flush=True)
    print("   phi (measured+noise) by iteration: "
          + "  ".join(f"{it}: {row['mean']:.1f} (base {row['base']:.1f}, {row['runs']} runs)"
                      for it, row in phis.items()))
    print(f"   phi truth {report['phi_truth']:.1f}  phi base {report['phi_base']:.1f}  "
          f"(base / truth {report['phi_base'] / report['phi_truth']:.2f})  "
          f"phi ratio {report['phi_ratio']:.4f}")
    for name, row in report["scalars"].items():
        inside = report["scalars_inside"][name]
        print(f"   {name:>14}: truth {row['truth']:.3f}  posterior 5-95% "
              f"[{row['posterior_p05']:.3f}, {row['posterior_p95']:.3f}]  z {row['z']:+.2f}"
              f"{'' if inside else '  OUTSIDE'}")
    forecast = report["forecast"]
    print(f"   forecast: truth {forecast['truth']:.3f} in {forecast['posterior_low']:.3f}.."
          f"{forecast['posterior_high']:.3f} ({forecast['members']} realizations): "
          f"{'yes' if forecast['covered'] else 'NO'}")
    print(f"   crisp gap closed {report['crisp_gap_closed']:.2f}  truth round trip "
          f"{report['truth_round_trip']:.1e}  write-back {report['write_back']:.1e}  "
          f"equilibrium {'yes' if report['equilibrium'] else 'NO'}", flush=True)
    for name, row in report.get("linear", {}).items():
        print(f"   linear {name:>14}: sd ensemble {row['ies_sd']:.4f} linear {row['lin_sd']:.4f}  "
              f"z linear {row['z_lin']:+.2f}  sd of z under the data's noise {row['z_sd']:.2f} "
              f"(levels only {row['z_sd_levels']:.2f})", flush=True)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("configs", nargs="*", default=list(CONFIGS),
                        help="REALSxITERATIONS[+lambda][+nodiff] (default: %(default)s)")
    parser.add_argument("--seeds", default=",".join(map(str, SEEDS)),
                        help="comma-separated seeds (default: %(default)s)")
    parser.add_argument("--workers", type=int, default=None,
                        help="PESTPP-IES workers (default: the notebook's own)")
    parser.add_argument("--linear", action="store_true",
                        help="after each run, the linear posterior at its base (linear_check)")
    parser.add_argument("--keep", type=pathlib.Path, default=None,
                        help="keep each run's workspace and report (report.json) under this "
                             "directory instead of a temporary one")
    parser.add_argument("--reports", type=pathlib.Path, default=None,
                        help="run nothing: pool the report.json files a --keep run left under "
                             "this directory, by configuration")
    args = parser.parse_args(argv)
    if args.reports is not None:
        rows = []
        for path in sorted(args.reports.glob("*/report.json")):
            label, _, seed = path.parent.name.rpartition("_")
            rows.append((label, int(seed), float("nan"), json.loads(path.read_text(encoding="utf-8"))))
        _print_pooled(rows)
        return
    seeds = [int(seed) for seed in args.seeds.split(",")]
    configs = [(config, *parse_config(config)) for config in args.configs]
    keep = args.keep.resolve() if args.keep is not None else None  # before the chdir
    os.chdir(NOTEBOOKS)
    sys.path.insert(0, str(NOTEBOOKS))
    rows = []
    with tempfile.TemporaryDirectory(prefix="rehearsal_budget_") as scratch:
        base = keep if keep is not None else pathlib.Path(scratch)
        for label, reals, iterations, flags in configs:
            for seed in seeds:
                workspace = base / f"{label.replace('+', '_')}_{seed}"
                start = time.time()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    report = run(reals, iterations, seed, workspace, flags=flags,
                                 workers=args.workers, linear=args.linear)
                minutes = (time.time() - start) / 60
                _print_run(label, seed, minutes, report)
                if args.keep is not None:
                    (workspace / "report.json").write_text(
                        json.dumps(report, indent=1, default=float), encoding="utf-8")
                rows.append((label, seed, minutes, report))
    print(f"\n{'config':>16} {'seed':>4} {'min':>5} {'phi ratio':>9} {'crisp gap':>9} "
          f"{'forecast':>8} {'round trip':>10} {'write-back':>10} {'equil':>5} {'phi truth':>9} "
          f"{'phi base':>8}  outside the 5-95% band")
    for label, seed, minutes, report in rows:
        outside = [f"{name} (z {report['scalars'][name]['z']:+.1f})"
                   for name, inside in report["scalars_inside"].items() if not inside]
        print(f"{label:>16} {seed:>4} {minutes:5.1f} "
              f"{report['phi_ratio']:9.4f} {report['crisp_gap_closed']:9.2f} "
              f"{'yes' if report['forecast']['covered'] else 'NO':>8} "
              f"{report['truth_round_trip']:10.1e} {report['write_back']:10.1e} "
              f"{'yes' if report['equilibrium'] else 'NO':>5} "
              f"{report['phi_truth']:9.1f} {report['phi_base']:8.1f}  {', '.join(outside) or '-'}")
    _print_pooled(rows)


def _print_pooled(rows) -> None:
    """Each configuration's runs, seeds pooled: ``pooled_recovery``'s verdict, and how
    often an exact posterior fails its thresholds -- by the binomial it uses and with
    each run's z correlated as its posterior is."""

    from myflopy.modflow.mf6 import canonical_rehearsal as rehearsal

    for label in dict.fromkeys(label for label, *_ in rows):
        reports = [report for name, _, _, report in rows if name == label]
        verdict = rehearsal.pooled_recovery(reports)
        correlated = correlated_false_fail(reports, verdict)
        print(f"\npooled {label} ({verdict.runs} runs): "
              f"{'PASSES' if verdict.passed else 'FAILS'}; an exact posterior fails it "
              f"at most {verdict.false_fail['total']:.3f} (binomial; correlated as measured: "
              f"coverage {correlated['coverage']:.3f}, largest |z| {correlated['largest_z']:.3f})")
        print(verdict.summary().to_string())


if __name__ == "__main__":
    main()
