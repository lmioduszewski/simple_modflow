"""Array recharge per period: what one forward run's multiplier step costs at Ten Trails scale.

Ledger 547 and 549. ``recharge`` on an array RCH (READASARRAYS,
what ``mf.rch.array`` and Ten Trails write) scales one ``recharge`` file per period
BLOCK; ``periods="month"`` groups them into one parameter set per calendar month; pilot
points on recharge interpolate ONE field per run and multiply it into every period file.
This builds a template the size of the Ten Trails average-year design -- 73 periods (one
steady day, then 6 x 12 months) of 18,107 cells -- with four stacked recharge parameters,
writes one realization into a copy of it, and times the two pre-model steps of its
``forward_run.py`` separately, then MF6 itself for scale, and checks what MF6 was handed.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/array_recharge_period_apply.py

Uses only the public API plus the template's own forward-run lines, so it measures
whichever tree is on ``PYTHONPATH`` (made absolute at start, ledger 425). Needs ``mf6``.
A few minutes; on the parent commit (claude/pestprep 415201b) the build refuses array
recharge outright ("it is the ARRAY form (READASARRAYS), which this target does not
scale").

THE SETUP. DISV, one layer of 19 x 953 = 18,107 square cells (Ten Trails' transient
runs carry 18,107 per layer), ``mf.AverageYear(6, start="2017-10-01", nstp=1)``, an array
RCH with a block every period (``by_month``), one head target. Parameters:
``rch`` constant (x1.5), ``rchm`` ``periods="month"`` (m03 x2, every other month x1),
``rchz`` zone (two zones, x1), ``rchpp`` pilot points every 3,000 length units (x1).

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, 2026-10-12, nice 10, load average 11-13)
  period files scaled                               73
  pyEMU apply_list_and_array_pars (chunk_len 1e6)   2.45 s
  pilot points, one command for all 73 files       1.31 s
  MF6 (73 periods, one step each)                   3.05 s
  a March file written / its base                   3.000000  (1.5 x 2)
  the steady file written / its base                1.500000  (no seasonal multiplier)
  template size                                     76 MB

The design pass measured 6.55 s for the pyEMU apply of 73 x 18,107 with a process pool
per run; ``chunk_len=10**6`` (ledger 511's PstFrom default) keeps it in one process.
The pilot-point step is file I/O, not interpolation: the factor field is one gather and
row sum (a few ms); reading and writing 73 files of 18,107 values is the rest.
"""

from __future__ import annotations

import json
import os
import pathlib
import shutil
import sys
import tempfile
import time
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import flopy  # noqa: E402
import geopandas as gpd  # noqa: E402
import pandas as pd  # noqa: E402
import pyemu  # noqa: E402
from flopy.utils.gridutil import get_disv_kwargs  # noqa: E402

import myflopy as mf  # noqa: E402
from myflopy.modflow.mf6.pest import forward_run  # noqa: E402

NROW, NCOL = 19, 953
NCPL = NROW * NCOL
YEAR = mf.AverageYear(6, start="2017-10-01", nstp=1)


def _run(workspace: pathlib.Path) -> pathlib.Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=YEAR.nper, perioddata=YEAR.perioddata, time_units="days",
                          start_date_time=YEAR.start_date_time)
    flopy.mf6.ModflowIms(sim, complexity="moderate")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(
        1, NROW, NCOL, np.full(NCOL, 50.0), np.full(NROW, 50.0), 100.0, [0.0]))
    flopy.mf6.ModflowGwfnpf(gwf, k=20.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1e-5, iconvert=0,
                            steady_state=YEAR.steady_state, transient=YEAR.transient)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 100.0], [(0, NCPL - 1), 95.0]])
    monthly = {month: np.full(NCPL, 1e-4 * (1 + month / 12)) for month in range(1, 13)}
    flopy.mf6.ModflowGwfrcha(gwf, recharge=YEAR.by_month(monthly))
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "LAST")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


def _timed(label: str, call) -> float:
    started = time.perf_counter()
    call()
    seconds = time.perf_counter() - started
    print(f"  {label:<48}{seconds:6.2f} s")
    return seconds


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="array_recharge_apply_"))
    try:
        run = _run(root / "run")
        model = mf.load_mf6_run(run)
        cal = model.pest("rch", workspace=root / "cal")
        points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                                  geometry=gpd.points_from_xy([500.0], [500.0]))
        cal.observe(mf.HeadTargets(locations=points, values=pd.DataFrame(
            {"time": [1], "W": [99.0]})))
        cal.parameterize("recharge", name="rch", bounds=(0.5, 2.0))
        cal.parameterize("recharge", name="rchm", periods="month", bounds=(0.5, 2.0))
        zones = (np.arange(NCPL) % NCOL >= NCOL // 2).astype(int) + 1
        cal.parameterize("recharge", style="zone", zones=zones, name="rchz", bounds=(0.5, 2.0))
        cal.parameterize("recharge", style="pilotpoints", name="rchpp", pp_spacing=3000.0,
                         bounds=(0.5, 2.0))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pst = cal.build("rch.pst", noptmax=0)

        worker = root / "worker"
        shutil.copytree(cal.template_workspace, worker)
        par = pst.parameter_data.copy()
        par.loc[par["pargp"] == "rch", "parval1"] = 1.5
        par.loc[par["pargp"] == "rchm_m03", "parval1"] = 2.0
        for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
            pyemu.pst_utils.write_to_template(par["parval1"], str(worker / tpl),
                                              str(worker / inp))
        main_text = (worker / "forward_run.py").read_text()
        main_text = main_text[main_text.index("def main():"):]
        apply_line = next(line.strip() for line in main_text.splitlines()
                          if "apply_list_and_array_pars(" in line)
        pilot_lines = [line.strip() for line in main_text.splitlines()
                       if "apply_pilotpoints_to_array(" in line]
        namespace = {"pyemu": pyemu,
                     "apply_pilotpoints_to_array": forward_run.apply_pilotpoints_to_array}
        manifest = json.loads((worker / "pilot_point_sets.json").read_text())
        files = sorted({entry["array_file"] for entry in manifest["files"]})
        print(f"  {'period files scaled':<48}{len(files):6d}")
        cwd = pathlib.Path.cwd()
        os.chdir(worker)
        try:
            _timed("pyEMU apply_list_and_array_pars", lambda: exec(apply_line, namespace))  # noqa: S102
            _timed(f"pilot points, {len(pilot_lines)} command(s)",
                   lambda: [exec(line, namespace) for line in pilot_lines])  # noqa: S102
            _timed("MF6", lambda: flopy.mbase.run_model("mf6", None, model_ws=".", silent=True))
        finally:
            os.chdir(cwd)

        march = YEAR.months_of(0)[3]
        spec = next(spec for spec in cal._native_parameter_specs if spec.name == "rch")
        by_period = {target.period: target.filename for target in spec.resolved_targets}
        for label, period in (("a March file written / its base", march),
                              ("the steady file written / its base", 0)):
            written = np.loadtxt(worker / by_period[period]).reshape(-1)
            base = np.loadtxt(cal.template_workspace / "org" / by_period[period]).reshape(-1)
            print(f"  {label:<48}{float(np.median(written / base)):.6f}")
        size = sum(path.stat().st_size for path in cal.template_workspace.rglob("*")
                   if path.is_file())
        print(f"  {'template size':<48}{size / 1e6:6.0f} MB")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
