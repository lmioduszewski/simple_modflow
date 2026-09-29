"""The review of ledger 346: what the observation route's period keying cannot see (ledger 347).

Ledger 346 keyed the observation route by stress period, one row per period at its
end, and refuses a target keyed by MF6's simulated time -- the route's key until then
-- naming the period it ends. This script measures the three things the review found
about that, on a real MF6 run:

    /home/lukem/python/envs/gw/bin/python docs/reproducers/obs_route_period_end_review.py

It uses the public API plus ``period_end.period_end_rows`` and the reader's
``_observation_output_table``, so it measures whichever tree is on ``PYTHONPATH`` --
section B's "before" is ``87acb4c``'s ``src`` (``git archive 87acb4c src | tar -x -C
<dir>``, then ``PYTHONPATH=<dir>/src``); A, C and D read the same on both. Needs ``mf6``.

THE SETUP. The DRN part of ``tests/test_obs_route_period_end.py``'s toy: a 6 x 6
irregular Voronoi grid of ~100 ft cells, one confined layer, CHD at 102 ft west and
east, heads filling from 90 ft, and a three-cell DRN zone ``seep`` (cells 20, 21, 27)
at 96 ft declared as DRN observations, so it takes the observation route.

  A. A time key BELOW nper. Twelve 2-day periods, three steps each: the pre-346
     spelling keys the ends of periods 0, 1, 2 as t = 2, 4, 6, which are also the
     period numbers 2, 4, 6. What ``compare()`` reads for them, and for 2, 4, 24.
  B. MF6's time column. The toy with three 10-day periods of three steps, its DRN
     observations declared over one TEXT block (cell 20) and one BINARY block (cells
     21, 27): the time each file holds for one step, and what the zone reads.
  C. That text rounding against ``END_RTOL``: MF6's period ends (the 8-decimal
     lengths FloPy writes, added up) written at 13 significant digits, through
     ``period_end_rows``, for long runs of short periods.
  D. ``PRECISION SINGLE`` in an OBS file's options, as ledger 346 supposed MF6 took;
     and ``DIGITS 4``, which sets the digits of the VALUES.

MEASURED (mf6 6.8.0.dev0, 2026-09-29)
  A  compare() keyed 2, 4, 6        reads periods 2, 4, 6: -630.055, -808.923, -852.023
                                    (the ends of periods 0, 1, 2 are 0.0, -383.248,
                                    -630.055) -- in silence, on both trees
     compare() keyed 2, 4, 24       REFUSED: "24 is the end of period 11 ... key it 11"
  B  one step's time                text 3.333333333333, binary 3.3333333333333335
                                    (13 significant digits against every digit; the
                                    values beside it are written with every digit)
     the zone over both files       87acb4c: REFUSED, "hold different times"; after:
                                    -710.119, -853.68, -864.954, equal to the one-file
                                    declaration to the bit
  C  hourly, 1 year (8,760 ends)    reads all; max |pertim - perlen| 3.33e-09
     hourly, 10 years (87,600)      reads all; 4.33e-09
     1/3 day, 30 years (32,850)     reads all; 4.33e-09
     hourly, 30 years (262,800)     REFUSED ("1 row(s) past the end")
     0.1 day, 50 years (182,500)    REFUSED ("no row in stress period(s) 116909, ...")
  D  PRECISION SINGLE               MF6 stops: "Unrecognized option: PRECISION"
     DIGITS 4                       first row 3.333333333333,-4.576,-5.151,-17.34: the
                                    values at 4 digits, the time still at 13

Why each is so, and what was decided: ledger 347.
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations._shared import _observation_output_table
from myflopy.modflow.mf6.period_end import period_end_rows

NSTP = 3
DRN_CELLS = [20, 21, 27]
MF6 = shutil.which("mf6") or str(pathlib.Path.home() / ".local/bin/mf6")


def _toy(root: pathlib.Path, *, perlen: float, nper: int):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    drain = mf.drn(
        stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]},
        observations=mf.DrnFlowTargets(locations={"seep": DRN_CELLS}).to_flopy_obs(
            None, csv_name="drn_flow_targets.csv"))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0), mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        drain,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(perlen, NSTP, 1.0)] * nper)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _zone(model, values=None):
    return mf.DrnFlowTargets(locations={"seep": DRN_CELLS}, values=values)


def section_a(root: pathlib.Path) -> None:
    print("\nA. a time key below nper: twelve 2-day periods")
    model = _toy(root / "twelve", perlen=2.0, nper=12)
    series = _zone(model).simulated_series(model)
    print(f"   the ends of periods 0, 1, 2: {np.round(series['seep'].to_numpy()[:3], 3).tolist()}")
    old = _zone(model, pd.DataFrame({"time": [2, 4, 6], "seep": -700.0})).compare(model)
    read = old.sort_values("per")
    print(f"   keyed 2, 4, 6 (the old spelling of those ends): read periods "
          f"{read['per'].tolist()}, {np.round(read['sim_flow'].to_numpy(), 3).tolist()}")
    try:
        _zone(model, pd.DataFrame({"time": [2, 4, 24], "seep": -700.0})).compare(model)
    except ValueError as error:
        print(f"   keyed 2, 4, 24: REFUSED -- ...{str(error)[str(error).find('24 is'):][:60]}")


def section_b(root: pathlib.Path) -> None:
    print("\nB. MF6's time column: one zone over a text and a binary file")
    model = _toy(root / "steps", perlen=10.0, nper=3)
    copy = root / "split"
    shutil.copytree(model.workspace, copy)
    records = [f"  seep_c{cell}  DRN  1 {cell + 1}" for cell in DRN_CELLS]
    (copy / "m.drn.obs").write_text(
        "BEGIN options\nEND options\n\n"
        f"BEGIN continuous  FILEOUT  drn_a.csv\n{records[0]}\nEND continuous\n\n"
        "BEGIN continuous  FILEOUT  drn_b.bsv  BINARY\n"
        + "\n".join(records[1:]) + "\nEND continuous\n")
    (copy / "drn_flow_targets.csv").unlink()
    subprocess.run([MF6], cwd=copy, capture_output=True, check=True)
    text = (copy / "drn_a.csv").read_text().splitlines()[1].split(",")[0]
    binary = _observation_output_table(copy / "drn_b.bsv", owner="b")[1][0].to_numpy()
    text_time = _observation_output_table(copy / "drn_a.csv", owner="b")[1][0].to_numpy()
    print(f"   one step's time: text {text}, binary {float(binary[0])!r}")
    print(f"   equal bit for bit (the pre-347 check): {np.array_equal(text_time, binary)}; "
          f"within 1e-12 of the time: {np.allclose(text_time, binary, rtol=1e-12, atol=0.0)}")
    reloaded = mf.load_mf6_run(copy, verbosity_level=0)
    try:
        split = _zone(reloaded).simulated_series(reloaded)
        whole = _zone(model).simulated_series(model)
        print(f"   the zone over both files: {np.round(split['seep'].to_numpy(), 3).tolist()}; "
              f"equal to the one-file declaration: {np.array_equal(split['seep'], whole['seep'])}")
    except ValueError as error:
        print(f"   the zone over both files: REFUSED -- {str(error)[:110]}")


def section_c() -> None:
    print("\nC. MF6's 13-digit time against END_RTOL (synthetic, as MF6 adds and writes it)")
    for label, length, nper in [
        ("hourly, 1 year", 1 / 24, 24 * 365),
        ("hourly, 10 years", 1 / 24, 24 * 3650),
        ("1/3 day, 30 years", 1 / 3, 3 * 365 * 30),
        ("hourly, 30 years", 1 / 24, 24 * 365 * 30),
        ("0.1 day, 50 years", 0.1, 10 * 365 * 50),
    ]:
        ends = np.cumsum(np.full(nper, float(f"{length:15.8f}")))
        times = [float(f"{end:.13G}") for end in ends]
        try:
            period_end_rows(times, [length] * nper)
            worst = np.abs(np.diff(np.concatenate([[0.0], times])) - length).max()
            print(f"   {label:20s} ({nper:>7,} ends): reads all; max |pertim - perlen| {worst:.2e}")
        except ValueError as error:
            print(f"   {label:20s} ({nper:>7,} ends): REFUSED -- {str(error)[:90]}")


def section_d(root: pathlib.Path) -> None:
    print("\nD. PRECISION SINGLE in an OBS file")
    model = _toy(root / "single", perlen=10.0, nper=3)
    copy = root / "single_ws"
    shutil.copytree(model.workspace, copy)
    records = "\n".join(f"  seep_c{cell}  DRN  1 {cell + 1}" for cell in DRN_CELLS)
    (copy / "m.drn.obs").write_text(
        "BEGIN options\n  PRECISION SINGLE\nEND options\n\n"
        f"BEGIN continuous  FILEOUT  drn.bsv  BINARY\n{records}\nEND continuous\n")
    run = subprocess.run([MF6], cwd=copy, capture_output=True, text=True)
    error = [line.strip() for line in run.stdout.splitlines() if "PRECISION" in line]
    print(f"   PRECISION SINGLE: MF6 exit code {run.returncode}: {error[:1]}")
    (copy / "m.drn.obs").write_text(
        "BEGIN options\n  DIGITS 4\nEND options\n\n"
        f"BEGIN continuous  FILEOUT  drn.csv\n{records}\nEND continuous\n")
    subprocess.run([MF6], cwd=copy, capture_output=True, check=True)
    print(f"   DIGITS 4: first row {(copy / 'drn.csv').read_text().splitlines()[1]!r}")


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="obs_route_period_end_review_"))
    try:
        section_a(root)
        section_b(root)
        section_c()
        section_d(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
