"""A target set can be a one-sided BOUND: what PESTPP-IES does with one, and the review.

Ledger 564. ``inequality="at_most"``/``"at_least"`` on a target set puts
its measured observations in PESTPP-IES's ``l_``/``g_`` groups; the IES review scores
them as PEST++ does. This script measures the PEST++ behaviour the design rests on and
the review before and after.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_inequality_targets.py

Needs ``mf6`` and ``pestpp-ies``.

THE SETUP.
  A. PEST++ alone (no myflopy): one parameter ``x`` (start 8, bounds 0.1-20, no
     transform); a "model" writing ``x`` three times; observations ``o_eq`` (8, weight
     0), ``o_lt`` in group ``l_hds`` (5, weight 2: at most 5) and ``o_gt`` in ``g_hds``
     (3, weight 2: at least 3). PESTPP-IES NOPTMAX=3, 30 reals, seed 7 -- once with no
     noise option, once with a ``standard_deviation`` of 0.5 on both bounds. Reported:
     the noise ensemble's mean and spread of each bound, and the mean phi before and
     after.
  B. ``tests/test_pest_handed_targets.py``'s toy (6 x 6 Voronoi, heads ~100-102 ft):
     ``cap`` -- W1/W2 at most 1000 ft (held by every realization), ``low`` -- W3 at most
     1 ft (violated by every one), and an ordinary W4 at 95 ft; ``cal.prior(reals=3)``.
     Reported: PESTPP-IES's ``phi.group.csv`` base row, the review's
     ``phi_contributions()``, the same scored two-sided (the rule before), and the
     share of ``l_cap`` observations ``conflict()`` flags.
  C. Which weights ``phi.group.csv`` holds (review of ledger 564): A's
     model with ``o_eq`` measured at 6 (weight 1), NOPTMAX=1, 5 reals, once plainly and
     once with an ``ies_phi_factor_file`` (``eq`` 0.9, ``l_hds`` 0.1, ``g_hds`` 0.2).
     Reported: the base row's ``eq`` at iteration 0 against ``(8 - 6)**2 = 4`` at the
     control file's weight. (At NOPTMAX=0 PESTPP-IES applies no phi factor: 4.0 both.)

MEASURED (mf6 6.8.0.dev0, PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-12)
  A  no noise option     PESTPP-IES "resetting to `ies_no_noise` to true"; noise
                         ensemble o_lt mean 5.000 sd 0.000, o_gt mean 3.000 sd 0.000;
                         mean phi 143.6 -> 1.49, one realization at phi 0 (both
                         bounds held: the satisfied band scores 0)
     std 0.5             o_lt mean 4.945 sd 0.513, o_gt mean 3.030 sd 0.421: noise IS
                         drawn onto a bound -- passed through, not suppressed;
                         mean phi 143.6 -> 2.12, six realizations at phi 0
  B  phi.group.csv base  l_cap 0, l_low 10,187.4 (= (sim - 1)^2 by hand, 10,187.4),
                         hds 49.0
     phi_contributions   before (two-sided): l_cap 2,423,185.3 -- the dominant group,
                         for a bound every realization holds; after: l_cap 0,
                         l_low 10,187.4, hds 49.0 (read from the file)
     conflict, l_cap     before 3 of 3 in conflict; after 0 of 3
  C  phi.group.csv eq    plain: 4.0 (the control file's weight); with the phi factor
                         file: 21.6 -- PESTPP-IES's ADJUSTED weights, so
                         ``phi_contributions()`` reports the balance the run fitted while
                         its fallback (no file) recomputes at the control file's weights
(the measured numbers are printed by ``main()``; the table above is one run's.)
"""

from __future__ import annotations

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

# PEST++ starts each forward run in the template, where a relative PYTHONPATH would
# import another myflopy: made absolute here, before anything changes directory
# (ledger 425).
make_pythonpath_absolute()

import os  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402
import warnings  # noqa: E402

import pandas as pd  # noqa: E402
import pyemu  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy.modflow.mf6.pest.ies import IesResults  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from tests.test_pest_handed_targets import run_toy  # noqa: E402

IES = shutil.which("pestpp-ies") or str(pathlib.Path.home() / ".local/bin/pestpp-ies")


def _one_parameter(root: pathlib.Path, *, noise: bool) -> dict:
    root.mkdir(parents=True)
    (root / "model.py").write_text(
        "x = float(open('par.dat').read().split()[0])\n"
        "open('out.dat', 'w').write(f'{x:.8f}\\n{x:.8f}\\n{x:.8f}\\n')\n"
    )
    (root / "par.dat.tpl").write_text("ptf ~\n~   x    ~\n")
    (root / "out.dat.ins").write_text("pif ~\nl1 !o_eq!\nl1 !o_lt!\nl1 !o_gt!\n")
    (root / "par.dat").write_text("8.0\n")
    subprocess.run([sys.executable, "model.py"], cwd=root, check=True)
    here = pathlib.Path.cwd()
    try:
        os.chdir(root)
        pst = pyemu.Pst.from_io_files(["par.dat.tpl"], ["par.dat"], ["out.dat.ins"], ["out.dat"])
    finally:
        os.chdir(here)
    par = pst.parameter_data
    par.loc["x", ["parval1", "parlbnd", "parubnd", "partrans"]] = [8.0, 0.1, 20.0, "none"]
    obs = pst.observation_data
    obs.loc["o_eq", ["obsval", "weight", "obgnme"]] = [8.0, 0.0, "eq"]
    obs.loc["o_lt", ["obsval", "weight", "obgnme"]] = [5.0, 2.0, "l_hds"]
    obs.loc["o_gt", ["obsval", "weight", "obgnme"]] = [3.0, 2.0, "g_hds"]
    if noise:
        obs["standard_deviation"] = [0.5 if name != "o_eq" else 1.0 for name in obs.index]
    # PEST++ drops a command's leading "/", so the model runs on the system python3
    # (it imports nothing).
    pst.model_command = ["python3 model.py"]
    pst.control_data.noptmax = 3
    pst.pestpp_options.update({"ies_num_reals": 30, "random_seed": 7})
    pst.write(str(root / "t.pst"), version=2)
    subprocess.run([IES, "t.pst"], cwd=root, check=True, capture_output=True)
    drawn = pd.read_csv(root / "t.obs+noise.csv", index_col=0)
    phi = pd.read_csv(root / "t.phi.actual.csv")
    return {
        "o_lt": (drawn["o_lt"].mean(), drawn["o_lt"].std()),
        "o_gt": (drawn["o_gt"].mean(), drawn["o_gt"].std()),
        "phi": (phi["mean"].iloc[0], phi["mean"].iloc[-1]),
        "zero_phi_reals": int((phi.iloc[-1, 6:].astype(float) == 0).sum()),
    }


def section_c(root: pathlib.Path) -> None:
    for factors in (None, "eq,0.9\nl_hds,0.1\ng_hds,0.2\n"):
        case = root / f"c_{'factors' if factors else 'plain'}"
        case.mkdir(parents=True)
        (case / "model.py").write_text(
            "x = float(open('par.dat').read().split()[0])\n"
            "open('out.dat', 'w').write(f'{x:.8f}\\n{x:.8f}\\n{x:.8f}\\n')\n"
        )
        (case / "par.dat.tpl").write_text("ptf ~\n~   x    ~\n")
        (case / "out.dat.ins").write_text("pif ~\nl1 !o_eq!\nl1 !o_lt!\nl1 !o_gt!\n")
        (case / "par.dat").write_text("8.0\n")
        subprocess.run([sys.executable, "model.py"], cwd=case, check=True)
        here = pathlib.Path.cwd()
        try:
            os.chdir(case)
            pst = pyemu.Pst.from_io_files(["par.dat.tpl"], ["par.dat"], ["out.dat.ins"], ["out.dat"])
        finally:
            os.chdir(here)
        pst.parameter_data.loc["x", ["parval1", "parlbnd", "parubnd", "partrans"]] = [8.0, 0.1, 20.0, "none"]
        obs = pst.observation_data
        obs.loc["o_eq", ["obsval", "weight", "obgnme"]] = [6.0, 1.0, "eq"]
        obs.loc["o_lt", ["obsval", "weight", "obgnme"]] = [5.0, 2.0, "l_hds"]
        obs.loc["o_gt", ["obsval", "weight", "obgnme"]] = [3.0, 2.0, "g_hds"]
        pst.model_command = ["python3 model.py"]
        pst.control_data.noptmax = 1
        pst.pestpp_options.update({"ies_num_reals": 5, "random_seed": 7})
        if factors:
            (case / "pf.csv").write_text(factors)
            pst.pestpp_options["ies_phi_factor_file"] = "pf.csv"
        pst.write(str(case / "t.pst"), version=2)
        subprocess.run([IES, "t.pst"], cwd=case, check=True, capture_output=True)
        table = pd.read_csv(case / "t.phi.group.csv")
        base = table.loc[(table["par_realization"].astype(str) == "base")
                         & (table["iteration"] == 0)].iloc[0]
        print(f"C  phi factors {'on ' if factors else 'off'}: phi.group.csv base eq "
              f"{float(base['eq']):.1f} (control-file weight: 4.0)")


def _heads(values, *, names, cells, inequality=None):
    locations = pd.DataFrame({"name": list(names), "cell": list(cells), "layer": 0})
    return mf.HeadTargets(locations=locations, values=pd.DataFrame(values), inequality=inequality)


def section_b(root: pathlib.Path) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = run_toy(root / "toy")
        cal = model.pest("bounded", workspace=root / "pest", start_datetime="2024-01-01")
        cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
        cal.observe(_heads({"time": [0, 1], "W1": [1000.0, 1000.0], "W2": [1000.0, None]},
                           names=("W1", "W2"), cells=(1, 2), inequality="at_most"), prefix="cap")
        cal.observe(_heads({"time": [2], "W3": [1.0]}, names=("W3",), cells=(34,),
                           inequality="at_most"), prefix="low")
        cal.observe(_heads({"time": [0], "W4": [95.0]}, names=("W4",), cells=(23,)))
        cal.build("bounded.pst")
        ies = cal.prior(reals=3)
    table = pd.read_csv(ies.workspace / "bounded.phi.group.csv")
    base = table.loc[table["par_realization"].astype(str) == "base"].iloc[0]
    print("B  phi.group.csv base:", {g: float(base[g]) for g in ("l_cap", "l_low", "hds")})
    simulated = pd.read_csv(cal.template_workspace / "low_simulated_heads.csv").set_index("per")
    print(f"   l_low by hand: {(simulated.loc[2, 'W3'] - 1.0) ** 2:.1f}")
    review = IesResults(ies.workspace, case_name="bounded")
    print("   phi_contributions (after):", review.phi_contributions().round(1).to_dict())
    obs = review.pst.observation_data
    obs = obs.loc[obs["weight"] > 0]
    sim = review.posterior._df.loc["base", obs.index].astype(float)
    two_sided = ((obs["weight"] * (sim - obs["obsval"])) ** 2).groupby(obs["obgnme"]).sum()
    print("   phi two-sided (before):", two_sided.round(1).to_dict())
    conflict = review.conflict()
    cap = conflict[conflict["obgnme"] == "l_cap"]
    before = ((cap["measured"] < cap["prior_lo"]) | (cap["measured"] > cap["prior_hi"])).sum()
    print(f"   conflict, l_cap: before {before} of {len(cap)}; after "
          f"{int(cap['in_conflict'].sum())} of {len(cap)}")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="pest_inequality_") as directory:
        root = pathlib.Path(directory)
        for noise in (False, True):
            found = _one_parameter(root / f"a_noise_{noise}", noise=noise)
            print(
                f"A  noise={noise}: o_lt mean {found['o_lt'][0]:.3f} sd {found['o_lt'][1]:.3f}; "
                f"o_gt mean {found['o_gt'][0]:.3f} sd {found['o_gt'][1]:.3f}; mean phi "
                f"{found['phi'][0]:.1f} -> {found['phi'][1]:.2f}; reals at phi 0: "
                f"{found['zero_phi_reals']}"
            )
        section_b(root)
        section_c(root)


if __name__ == "__main__":
    main()
