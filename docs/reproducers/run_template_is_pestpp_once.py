"""``cal.run_template`` is what PESTPP-IES runs once, and refuses a run of another mesh (ledger 605).

The canonical TESTING model (21x21, 6 periods) calibrated on ``kl0``/``k33l0`` against
its 8 head targets (``hd``), two one-sided bound sets -- ``land``: ``at_most``, weight 2,
three of six cells violated at the start; ``floor``: ``at_least``, weight 3, two of six
violated -- and a zero-weight forecast. Five parts:

A. ``cal.run_template`` (MODFLOW 6 in the copy), its ``outputs_from=`` route (the model's
   own run, no MODFLOW 6) and PESTPP-IES at NOPTMAX 0 in a launch copy of the template,
   at the starting values and at ``kl0 = 0.5``: phi per group against PESTPP-IES's own
   ``t.phi.group.csv``, and the largest simulated difference against its ``t.base.rei``.
B. That NOPTMAX-0 directory as a review: PESTPP-IES wrote a phi file and no ensemble, so
   ``open_ies_run(...).iterations`` is empty and ``phi_contributions`` raises; which is
   why ``run_ies(iterations=0)`` is refused and ``runs.execution_kind`` reads NOPTMAX 0
   as no execution.
C. ``outputs_from=`` a 22x22 run (484 cells, the same model name and clock) against the
   441-cell template: refused, naming both mesh fingerprints. Before ledger 605
   the rehearsal's ``run_template`` returned 48 finite heads off by up to 48.79 ft (median
   20.91 ft) with no error -- the instruction files read cells by position.
D. ``outputs_from=`` an earlier ``run_template`` folder (ledger 607): a calibration
   ``w`` on three wells ``a``/``b``/``c``, run once; the wells MOVED (other cells, same
   names) and ``w`` rebuilt; then ``outputs_from=`` the first folder. At da26644c its
   forward run, target map and instruction file replaced the rebuilt template's in the
   copy and it returned the OLD wells' heads, 58.04 ft from the MF6 route, with no
   error. Now refused: the folder is a PEST working directory.
E. An observation with no finite simulated value (ledger 607): the template's
   forward run made to write ``nan`` for one head. PESTPP-IES at NOPTMAX 0 fails the run
   ("control file parameter value run failed", exit 1); at da26644c ``run_template``
   returned a table with no error -- measured on a calibration of the same 8 wells alone
   (``kl0``, ``hd``), its ``groupby("group")["phi"].sum()`` read hd 10.33, the NaN row's
   phi skipped. Now refused (``RuntimeError``).

``--time-run DIR`` also times the discretization check on a COPY of a real run (the
check ``outputs_from=`` makes: a discretization-only ``MFSimulation.load`` plus
``disv_pin``).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/run_template_is_pestpp_once.py

MEASURED (2026-10-16, pyEMU 1.4.0, PESTPP-IES 5.2.16, MF6 6.8.0.dev0; load average ~16;
the design pass at 55947911 measured the same numbers):

  group     start: run_template  outputs_from  phi.group.csv  kl0=0.5: run_template  phi.group.csv
  hd                  11.192149     11.192149        11.1921              54.474205        54.4742
  l_land              72.000000     72.000000             72               3.078418        3.07842
  g_floor            108.000000    108.000000            108             368.561892        368.562
  fore                 0.000000      0.000000              0               0.000000              0
  largest |simulated - .base.rei|: 5.0e-10 (start), 4.8e-10 (kl0=0.5); outputs_from vs MF6 1.9e-09

  B: iterations [] ; phi_contributions raises FileNotFoundError ; execution_kind None
  C: refused -- model 'viz_prt_master': its mesh is DISV, 4 layer(s), fingerprint
     78695be464bde479 against the template's DISV, 4 layer(s), fingerprint
     49680e8054405558 (another mesh)

  D: refused -- is a PEST working directory, not a model run: it holds
     myflopy_pest_build.json, myflopy_template_run.json, forward_run.py
  E: pestpp-ies exit 1, "run failed" ; run_template refused -- 1 of 126 observation(s)
     read no finite simulated value

  The discretization check on a copy of Ten Trails ``ss`` (5 x 20,323 cells): 1.53 s
  (``--time-run``; the design pass measured 1.15 s at a lower load).
"""

from __future__ import annotations

import argparse
import contextlib
import io
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()


def _demo(root, config):
    from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo

    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_canonical_calibration_demo(root, config=config, n_head_wells=8)


def _calibration(root: Path):
    import numpy as np

    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
    from myflopy.modflow.mf6.observations import HeadTargets

    demo = _demo(root / "model", CanonicalModelConfig.testing())
    model = demo.model
    cells = [int(c) for c in demo.head_targets.locations_gdf["cell"].iloc[:6]]
    probe = HeadTargets.from_cells(cells=cells, values=np.full((model.nper, 6), np.nan),
                                   layers=0, names=[f"b{i}" for i in range(6)])
    start = probe.simulated_heads(model)[[f"b{i}" for i in range(6)]].to_numpy(dtype=float)
    land = HeadTargets.from_cells(cells=cells, values=start + np.array([-1, -1, -1, 1, 1, 1.0]),
                                  layers=0, names=[f"b{i}" for i in range(6)], weights=2.0,
                                  inequality="at_most")
    floor = HeadTargets.from_cells(cells=cells, values=start + np.array([1, 1, -1, -1, -1, -1.0]),
                                   layers=0, names=[f"g{i}" for i in range(6)], weights=3.0,
                                   inequality="at_least")
    cal = model.pest("t", start_datetime="2020-01-01", workspace=root / "pest_t")
    cal.parameterize("k", style="constant", layers=[0], name="kl0")
    cal.parameterize("k33", style="constant", layers=[0], name="k33l0")
    cal.observe(demo.head_targets, prefix="hd")
    cal.observe(land, prefix="land")
    cal.observe(floor, prefix="floor")
    cal.forecast(demo.forecast_targets, prefix="fore")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("t.pst")
    return demo, cal


def _pestpp_once(cal, directory: Path, values=None) -> Path:
    import pyemu

    cal._launch_copy(directory, method="reproducer", announce=False)
    if values:
        pst = pyemu.Pst(str(directory / "t.pst"))
        for name, value in values.items():
            pst.parameter_data.loc[name, "parval1"] = value
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pst.write(str(directory / "t.pst"), version=2)
    done = subprocess.run(["pestpp-ies", "t.pst"], cwd=directory, capture_output=True,
                          text=True, check=False)
    if done.returncode != 0:
        raise SystemExit(f"pestpp-ies failed in {directory}:\n{done.stdout[-2000:]}")
    return directory


def _pestpp_phi(directory: Path):
    import pandas as pd
    import pyemu

    groups = pd.read_csv(directory / "t.phi.group.csv")
    base = groups.loc[groups["par_realization"].astype(str) == "base"].iloc[0]
    columns = [c for c in groups.columns
               if c not in ("iteration", "total_runs", "obs_realization", "par_realization")]
    rei = pyemu.pst_utils.read_resfile(str(directory / "t.base.rei"))
    modelled = rei["modelled"].astype(float)
    modelled.index = modelled.index.str.lower()
    return {column: float(base[column]) for column in columns}, modelled


def part_a_and_b(root: Path, demo, cal) -> None:
    from myflopy.modflow.mf6.pest.ies import open_ies_run
    from myflopy.modflow.mf6.pest.runs import execution_kind

    kl0 = next(n for n in cal.pst.par_names if "kl0" in n and "k33" not in n)
    half = {n: (0.5 if n == kl0 else 1.0) for n in cal.pst.adj_par_names}
    start = cal.run_template(root / "start")
    outputs = cal.run_template(root / "outputs", outputs_from=demo.model)
    moved = cal.run_template(root / "half", values=half)
    phi_start, rei_start = _pestpp_phi(_pestpp_once(cal, root / "pestpp_start"))
    phi_half, rei_half = _pestpp_phi(_pestpp_once(cal, root / "pestpp_half", {kl0: 0.5}))
    mine = {name: table.groupby("group")["phi"].sum()
            for name, table in (("start", start), ("outputs", outputs), ("half", moved))}
    print(f"\nA. {'group':8s} {'start: run_template':>20s} {'outputs_from':>13s} "
          f"{'phi.group.csv':>14s} {'kl0=0.5: run_template':>22s} {'phi.group.csv':>14s}")
    for group in phi_start:
        print(f"   {group:8s} {mine['start'].get(group, 0.0):20.6f} "
              f"{mine['outputs'].get(group, 0.0):13.6f} {phi_start[group]:14.6g} "
              f"{mine['half'].get(group, 0.0):22.6f} {phi_half[group]:14.6g}")
    print("   largest |simulated - .base.rei|: "
          f"{(start['simulated'] - rei_start.reindex(start.index)).abs().max():.1e} (start), "
          f"{(moved['simulated'] - rei_half.reindex(moved.index)).abs().max():.1e} (kl0=0.5); "
          f"outputs_from vs MF6 {(outputs['simulated'] - start['simulated']).abs().max():.1e}")
    once = root / "pestpp_start"
    review = open_ies_run(once, case_name="t")
    try:
        review.phi_contributions()
        raised = "nothing"
    except FileNotFoundError as error:
        raised = f"FileNotFoundError ({str(error)[:60]}...)"
    print(f"\nB. iterations {review.iterations} ; phi_contributions raises {raised} ; "
          f"execution_kind {execution_kind(once, 't')}")


def part_c(root: Path, cal) -> None:
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig

    other = _demo(root / "big", CanonicalModelConfig(
        nrow=22, ncol=22, nlay=4, nper=6, cell_size=100.0, steps_per_period=1)).model
    try:
        cal.run_template(root / "other", outputs_from=other)
        print("\nC. NOT refused")
    except ValueError as error:
        line = next(line for line in str(error).splitlines() if "mesh" in line)
        print(f"\nC. refused -- {line.strip().lstrip('- ')}")


def part_d(root: Path, demo) -> None:
    import numpy as np

    from myflopy.modflow.mf6.observations import HeadTargets

    model = demo.model

    def build(cells):
        wells = HeadTargets.from_cells(cells=cells, values=np.full((model.nper, 3), 100.0),
                                       layers=0, names=["a", "b", "c"])
        cal = model.pest("w", start_datetime="2020-01-01", workspace=root / "pest_w")
        cal.parameterize("k", style="constant", layers=[0], name="kl0")
        cal.observe(wells, prefix="w")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("w.pst")
        return cal

    build([10, 50, 100]).run_template(root / "w_first")
    cal = build([200, 300, 400])  # the wells moved, the same names
    try:
        cal.run_template(root / "w_again", outputs_from=root / "w_first")
        print("\nD. NOT refused")
    except ValueError as error:
        said = str(error).split(" -- a copy")[0].split("' ", 1)[-1]
        print(f"\nD. refused -- {said}")


def part_e(root: Path, cal) -> None:
    forward_run = cal.template_workspace / "forward_run.py"
    original = forward_run.read_text()
    output = str(next(name for name in cal.pst.output_files if "hd" in str(name)))
    hook = (f"\n    import pandas as _pd\n    _f = _pd.read_csv({output!r})\n"
            f"    _f.iloc[0, 1] = float('nan')\n    _f.to_csv({output!r}, index=False, "
            "na_rep='nan')\n")
    cut = original.index("\nif __name__")
    forward_run.write_text(original[:cut].rstrip("\n") + hook + "\n" + original[cut:])
    try:
        directory = root / "pestpp_nan"
        cal._launch_copy(directory, method="reproducer", announce=False)
        done = subprocess.run(["pestpp-ies", "t.pst"], cwd=directory, capture_output=True,
                              text=True, check=False)
        failed = "run failed" in done.stdout or "run failed" in (directory / "t.rec").read_text()
        try:
            cal.run_template(root / "nan")
            verdict = "NOT refused"
        except RuntimeError as error:
            verdict = "refused -- " + str(error).split(" (")[0].split("): ")[-1]
    finally:
        forward_run.write_text(original)
    print(f"\nE. pestpp-ies exit {done.returncode}, \"run failed\" {failed} ; run_template {verdict}")


def time_run(directory: Path) -> None:
    from myflopy.modflow.mf6.pest.project import _discretization_identity, _load_discretization

    with tempfile.TemporaryDirectory() as scratch:
        copy = Path(scratch) / "run"
        shutil.copytree(directory, copy)
        started = time.perf_counter()
        simulation = _load_discretization(copy)
        identities = [_discretization_identity(simulation.get_model(name))
                      for name in simulation.model_names]
        print(f"\nThe discretization check on a copy of {directory}: "
              f"{time.perf_counter() - started:.2f} s ({identities})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--time-run", type=Path, default=None,
                        help="a real MF6 run directory to time the check on (copied first)")
    arguments = parser.parse_args()
    with tempfile.TemporaryDirectory() as scratch:
        root = Path(scratch)
        demo, cal = _calibration(root)
        part_a_and_b(root, demo, cal)
        part_c(root, cal)
        part_d(root, demo)
        part_e(root, cal)
    if arguments.time_run is not None:
        time_run(arguments.time_run)


if __name__ == "__main__":
    main()
