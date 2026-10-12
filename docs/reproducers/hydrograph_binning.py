"""Which reading a stress period's target is: nearest its END, against a monthly median (ledger 543).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/hydrograph_binning.py

Needs `mf6` on PATH (measured with 6.8.0.dev0, FloPy 3.10). About 3 s.

1. **Ten Trails' clock, read-only.** Every run under
   ``/home/lukem/models/mf6/tentrails/runs/`` declares ``START_DATE_TIME 2017-09-30``
   with monthly periods from October, so each period ENDS at 00:00 on the LAST day of
   its month (period 0: 2017-10-31). The first three periods are printed from the
   ``baseline`` TDIS when it is there (the file is only read). A reading taken on
   that last day therefore lies AFTER its month's period end: binned by containment,
   ``(start, end]``, it lands in the NEXT period; nearest-end keeps it in its own.

2. **A model's own hydrograph, binned back.** A 5 x 5 unconfined aquifer on that
   clock (12 monthly periods, 6 steps each, seasonal recharge, ``SAVE HEAD ALL``) is
   run, and its simulated head at EVERY time step, dated ``start + totim``, is handed
   to ``HeadTargets.from_hydrographs`` as the readings. Each target is then compared
   with the same run's period-END head (ledger 336). The rule that reads the period
   end returns it; a statistic over the period does not:

       rule          |target - simulated period-end head|, ft (measured 2026-10-12)
                     max        mean over 24 targets
       nearest_end   0.000000   0.000000
       median        0.270502   0.152388
       mean          0.292089   0.164369

   The worst is period 0, where the heads climb from their initial 50 ft all month.
   A monthly statistic of a CHANGING head is not its value at the period's end --
   which is what MF6 reports for the period -- so it biases a calibration by about
   half a month of change. ``tests/test_hydrograph_targets.py`` asserts the
   nearest-end case and that the median misses.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import shutil
import tempfile
import warnings

import flopy
import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.period_calendar import PeriodCalendar

TT_RUN = pathlib.Path("/home/lukem/models/mf6/tentrails/runs/baseline")
LENGTHS = [31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30]
NCPL = 25


def _grid(n=5, d=100.0):
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [[j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
              for j in range(n) for i in range(n)]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def ten_trails_period_ends() -> None:
    """The first periods of the Ten Trails baseline TDIS, read and never written."""

    if not (TT_RUN / "mfsim.nam").exists():
        print(f"(no {TT_RUN}; skipping the Ten Trails calendar)")
        return
    with contextlib.redirect_stdout(io.StringIO()):
        sim = flopy.mf6.MFSimulation.load(sim_ws=str(TT_RUN), load_only=["tdis"], verbosity_level=0)
    tdis = sim.tdis
    calendar = PeriodCalendar(
        tdis.start_date_time.get_data(),
        [row[0] for row in tdis.perioddata.get_data()],
        tdis.time_units.get_data(),
    )
    print("Ten Trails baseline TDIS:")
    print(calendar.table().head(3)[["per", "start", "end", "month"]].to_string(index=False))
    last_day = pd.Timestamp("2017-10-31 12:00")
    print(f"a reading at {last_day}: containment -> period "
          f"{int(np.searchsorted(calendar.ends.asi8, last_day.value, side='left'))}, "
          "nearest end -> period 0\n")


def run_history(root: pathlib.Path):
    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, NCPL), dtype=int))
    rates = {p: 2e-3 * (1.0 + np.cos(2 * np.pi * (p + 9) / 12)) for p in range(12)}
    model = mf.gwf("flow", context=context, packages=[
        mf.disv(nlay=1, top=100.0, botm=[0.0], **vor.get_disv_gridprops()),
        mf.ic(strt=50.0),
        mf.npf(k=5.0, icelltype=1),
        mf.sto(sy=0.05, ss=1e-5, iconvert=1, transient={0: True}),
        mf.ghb(stress_period_data={0: [[(0, c), 50.0, 100.0] for c in range(0, NCPL, 5)]}),
        mf.rch.array(recharge={p: np.full(NCPL, r) for p, r in rates.items()}, irch=0, context=context),
        mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
    ])
    tdis = mf.tdis(nper=12, perioddata=[(float(n), 6, 1.0) for n in LENGTHS],
                   time_units="days", start_date_time="2017-09-30")
    simulation = mf.simulation(model, name="s", packages=[tdis])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").prepare_run("r", simulation)
        success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run.model("flow")


def own_hydrograph(model) -> pd.DataFrame:
    hds = flopy.utils.HeadFile(str(pathlib.Path(model.workspace) / "flow.hds"))
    start = model.calendar.start
    rows = []
    for totim in hds.get_times():
        data = np.asarray(hds.get_data(totim=totim)).reshape(-1)
        rows.append({"date": start + pd.to_timedelta(totim, unit="D"), "W1": data[12], "W2": data[7]})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ten_trails_period_ends()
    root = pathlib.Path(tempfile.mkdtemp(prefix="hydrograph_binning_"))
    try:
        model = run_history(root)
        readings = own_hydrograph(model)
        locations = pd.DataFrame({"name": ["W1", "W2"], "cell": [12, 7], "layer": [0, 0]})
        print(f"{len(readings)} dated readings per well (every time step of 12 monthly periods)")
        for rule in ("nearest_end", "median", "mean"):
            targets = mf.HeadTargets.from_hydrographs(locations, readings, model=model, rule=rule)
            residual = targets.compare(model)["residual"].abs()
            worst = targets.compare(model).loc[residual.idxmax()]
            print(f"{rule:<12} max |target - period-end head| = {residual.max():.6f} ft "
                  f"(worst: {worst['name']} period {int(worst['per'])}), mean {residual.mean():.6f} ft")
    finally:
        shutil.rmtree(root, ignore_errors=True)
