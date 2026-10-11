"""A correlated pilot-point prior spreads as far as its bounds say: measured.

Ledger 424. ``run_ies()``/``prior()`` hand PESTPP-IES a prior ensemble in which each
correlated pilot-point set is drawn from its variogram. That draw took the variogram's
sill (1.0) as the log10 VARIANCE of every pilot value, whatever its bounds, while the
correlated ``grid`` parameters beside it -- drawn by pyEMU's ``PstFrom.draw`` -- are
scaled so the bounds span six standard deviations. PESTPP-IES holds the initial
ensemble to the bounds (``ies_enforce_bounds``), so every pilot value drawn outside
them ran AT a bound. Truncating the draw ourselves would change nothing PESTPP-IES does
not already do; the draw itself was the defect.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/pilot_point_prior_bounds.py

Uses only the public API, so it measures whichever tree is on ``PYTHONPATH`` (made
absolute at start, ledger 425); the "before" column is the parent commit's
``src``. Needs ``mf6`` and ``pestpp-ies``. About 5 minutes (50 realizations, serial).

THE SETUP. A 6 x 6 irregular Voronoi grid of ~100 ft cells, one confined layer, CHD at
102 ft on the west and east edges, heads filling from 90 ft, three 10-day periods.
``k`` by pilot points (``pp_space=1``) and ``k33`` by ``grid`` (one per cell), each
with ``correlation=300``, bounds ``(0.5, 2)``; two head targets. ``np.random.seed(0)``,
``prior(reals=50)``.

MEASURED (mf6 6.8.0.dev0, PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-11; "before" is da33185)
                                            before          after
  pilot points: drawn outside the bounds    76.2%           0.2%
  pilot points: log10 standard deviation    0.964           0.097
  pilot points: run AT a bound by IES       76.1%           0.3%
  grid (k33):   drawn outside the bounds     0.2%           0.2%
  grid (k33):   log10 standard deviation    0.104           0.104
  grid (k33):   run AT a bound by IES        0.2%           0.2%
  (log10(2 / 0.5) / 6 = 0.100)
"""

from __future__ import annotations

import pathlib
import shutil
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

# PEST++ starts each forward run in the template, where a relative PYTHONPATH would
# import another myflopy: made absolute here, before anything changes directory
# (ledger 425).
make_pythonpath_absolute()


def _run(root: pathlib.Path):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, k33=1.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(10.0, 1, 1.0)] * 3)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def main():
    root = pathlib.Path(tempfile.mkdtemp(prefix="pilot_point_prior_bounds_"))
    try:
        model = _run(root)
        cal = model.pest("bounds", workspace=root / "pest", start_datetime="2024-01-01")
        cal.parameterize("k", style="pilotpoints", pp_space=1, correlation=300.0,
                         bounds=(0.5, 2.0), name="kp")
        cal.parameterize("k33", style="grid", correlation=300.0, bounds=(0.5, 2.0), name="kg")
        cal.observe(mf.HeadTargets(
            locations=pd.DataFrame({"name": ["W", "V"], "cell": [14, 21], "layer": 0}),
            values=pd.DataFrame({"time": [0, 1, 2], "W": [98.0] * 3, "V": [97.0] * 3})))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("bounds.pst", noptmax=0)
        np.random.seed(0)
        results = cal.prior(reals=50)
        bounds = cal.pst.parameter_data
        drawn = pd.read_csv(cal.template_workspace / "bounds.prior_par.csv", index_col=0)
        ran = pd.read_csv(pathlib.Path(results.workspace) / "bounds.0.par.csv", index_col=0)
        ran = ran.loc[[label for label in ran.index if str(label) != "base"]]
        for label, prefix in (("pilot points", "kpl"), ("grid (k33)", "pname:kg")):
            columns = [name for name in drawn.columns if name.startswith(prefix)]
            low, high = bounds.loc[columns, "parlbnd"], bounds.loc[columns, "parubnd"]
            values = drawn[columns]
            outside = (values.lt(low, axis=1) | values.gt(high, axis=1)).to_numpy().mean()
            at_bound = (np.isclose(ran[columns], low) | np.isclose(ran[columns], high)).mean()
            print(f"{label:13s} drawn outside {outside:6.1%}   log10 std "
                  f"{np.log10(values.to_numpy()).std():.3f}   run at a bound {at_bound:6.1%}")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
