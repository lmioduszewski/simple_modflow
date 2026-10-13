"""Two PESTPP-IES options myflopy named wrongly: the noise switch and the seed.

A. ``run_ies(noise=True)`` must WRITE ``ies_no_noise false``: PESTPP-IES given no option
draws no noise. B. The seed of PESTPP-IES's own draws is ``random_seed``; it REJECTS
``ies_random_seed``, which ``seed=``'s refusal told the caller to pass.

Ledger 590. Ledger 572 made ``noise=True`` -- the default of
``run_ies``/``prior`` -- REMOVE ``ies_no_noise`` from the control file, on the premise
that PESTPP-IES draws measurement noise unless told not to. It does not: given no
noise-specific option it resets ``ies_no_noise`` to true and writes an "obs+noise"
ensemble that is the measured values in every realization. So every default run fitted
each realization to the same data, and its posterior bands were too narrow.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/ies_noise_and_seed_options.py

Needs ``pestpp-ies``.

THE SETUP (A). PEST++ alone, no myflopy: one parameter ``x`` (start 8, bounds 0.1-20, no
transform), a "model" writing ``x`` to three observations of weight 2 (a standard
deviation of 0.5), PESTPP-IES at NOPTMAX=-1 (the prior only), 30 realizations,
``random_seed`` 7 -- once with ``ies_no_noise`` absent (what ``noise=True`` wrote
before), once ``false`` (what it writes now), once ``true`` (``noise=False``).
Reported: the spread across realizations of each observation in ``t.obs+noise.csv``.

MEASURED (PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-12)
  absent   "resetting to `ies_no_noise` to true"; obs+noise sd 0.000 0.000 0.000
  false    obs+noise sd 0.460 0.508 0.498 (1 / weight = 0.5)
  true     obs+noise sd 0.000 0.000 0.000
On the recalibration rehearsal (``canonical_08``, 1,184 observations, 364 weighted) the
same: absent, every column of ``rehearsal.obs+noise.csv`` had sd 0; ``false``, up to
230 (a baseflow of weight 1 / (0.05 x its mean)).

THE SETUP (B). A's model, ``ies_no_noise`` absent, with ``ies_random_seed 7`` in place
of ``random_seed 7``. MEASURED: PESTPP-IES exits 1 -- "control file parsing error: the
following control data keyword lines were not accepted: ies_random_seed 7 ...
forgive_unknown_args is 'false'"; with ``random_seed 7`` it runs (A).
"""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import pandas as pd  # noqa: E402
import pyemu  # noqa: E402

IES = shutil.which("pestpp-ies") or str(pathlib.Path.home() / ".local/bin/pestpp-ies")


def noise_spread(root: pathlib.Path, option: bool | None, *, seed_key: str = "random_seed"):
    """The obs+noise ensemble's spread per observation, and whether PEST++ said it reset.

    ``None`` when PESTPP-IES refuses the control file (its exit code is not 0).
    """

    root.mkdir(parents=True)
    (root / "model.py").write_text(
        "x = float(open('par.dat').read().split()[0])\n"
        "open('out.dat', 'w').write(f'{x:.8f}\\n{x:.8f}\\n{x:.8f}\\n')\n"
    )
    (root / "par.dat.tpl").write_text("ptf ~\n~   x    ~\n")
    (root / "out.dat.ins").write_text("pif ~\nl1 !o_1!\nl1 !o_2!\nl1 !o_3!\n")
    (root / "par.dat").write_text("8.0\n")
    subprocess.run([sys.executable, "model.py"], cwd=root, check=True)
    here = pathlib.Path.cwd()
    try:
        os.chdir(root)
        pst = pyemu.Pst.from_io_files(["par.dat.tpl"], ["par.dat"], ["out.dat.ins"], ["out.dat"])
    finally:
        os.chdir(here)
    pst.parameter_data.loc["x", ["parval1", "parlbnd", "parubnd", "partrans"]] = [8.0, 0.1, 20.0, "none"]
    pst.observation_data.loc[:, "obsval"] = 8.0
    pst.observation_data.loc[:, "weight"] = 2.0
    # PEST++ drops a command's leading "/", so the model runs on the system python3
    # (it imports nothing).
    pst.model_command = ["python3 model.py"]
    pst.control_data.noptmax = -1
    pst.pestpp_options.update({"ies_num_reals": 30, seed_key: 7})
    if option is not None:
        pst.pestpp_options["ies_no_noise"] = option
    pst.write(str(root / "t.pst"), version=2)
    done = subprocess.run([IES, "t.pst"], cwd=root, check=False, capture_output=True, text=True)
    if done.returncode != 0:
        return None
    drawn = pd.read_csv(root / "t.obs+noise.csv", index_col=0).drop(index="base", errors="ignore")
    return [float(value) for value in drawn.std()], "resetting to `ies_no_noise`" in done.stdout


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="ies_default_noise_") as scratch:
        for label, option in (("absent", None), ("false", False), ("true", True)):
            spread, reset = noise_spread(pathlib.Path(scratch) / label, option)
            print(f"A {label:>7}  {'reset to no noise; ' if reset else ''}obs+noise sd "
                  + " ".join(f"{value:.3f}" for value in spread))
        refused = noise_spread(pathlib.Path(scratch) / "seed", None, seed_key="ies_random_seed")
        print("B ies_random_seed: " + ("PESTPP-IES refused the control file" if refused is None
                                       else "accepted"))


if __name__ == "__main__":
    main()
