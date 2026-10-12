"""A variogram's nugget reaches the prior PESTPP-IES receives: measured.

Ledger 489. ``cal.parameterize(..., correlation=, nugget=)`` builds an exponential
variogram of contribution 1.0 and adds the nugget on top. ``run_ies()``/``prior()``/
``draw_prior()`` draw the correlated prior through pyEMU's
``helpers.geostatistical_draws`` (``PstFrom.draw`` for ``style="grid"``, called
directly for pilot points), which builds the structure's covariance -- nugget on the
diagonal -- and then SETS each diagonal entry to the variance the bounds imply,
scaling the off-diagonal by the same standard deviations. The nugget was replaced
with the diagonal, so the draw was the same at every nugget. The project now scales
the structure to a sill of 1 (``geostats.unit_sill``), the rule pyEMU's own
``SpecSim2d`` applies, so the bounds give the total variance and the nugget its share.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/prior_nugget.py

Seconds; needs pyEMU only. It calls ``geostatistical_draws`` exactly as the project
does (``sigma_range=6``, ``scale_offset=True``) on the structure ``_geostruct_for``
built BEFORE (``build_geostruct`` of the declared spec) and AFTER (``unit_sill`` first).

THE SETUP. 40 log-transformed parameters on a line 10 ft apart, bounds ``(0.5, 2)``,
``correlation=30``; 4000 realizations, ``np.random.seed(0)``. Reported: the mean log10
correlation of neighbours (10 ft apart) and the mean log10 standard deviation.

MEASURED (pyEMU 1.4.0, 2026-10-11)
  nugget   neighbours' correlation        log10 sd         exp(-1/3)/(1+n)
           before      after              before  after
  0.0      0.716       0.716              0.100   0.100    0.717
  0.5      0.716       0.477              0.100   0.100    0.478
  1.0      0.716       0.358              0.100   0.100    0.358
  (log10(2 / 0.5) / 6 = 0.100; the same numbers through a full build, pilot points
  and ``style="grid"`` alike: ``tests/test_pilot_point_prior.py``)
"""

from __future__ import annotations

import contextlib
import io
import warnings

import numpy as np
import pandas as pd
import pyemu

from myflopy.modflow.mf6.pest.geostats import build_geostruct, unit_sill
from myflopy.modflow.mf6.pest.specs import ExpGeoStruct

N, SPACING, RANGE, REALS = 40, 10.0, 30.0, 4000


def measure(structure) -> tuple[float, float]:
    names = [f"p{i}" for i in range(N)]
    pst = pyemu.Pst.from_par_obs_names(names, ["o0"])
    par = pst.parameter_data
    par.loc[:, "partrans"] = "log"
    par.loc[:, "parval1"] = 1.0
    par.loc[:, "parlbnd"] = 0.5
    par.loc[:, "parubnd"] = 2.0
    par.loc[:, "pargp"] = "k"
    frame = pd.DataFrame({"parnme": names, "x": np.arange(N) * SPACING, "y": 0.0})
    np.random.seed(0)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        drawn = pyemu.helpers.geostatistical_draws(
            pst, {build_geostruct(structure): [frame]}, num_reals=REALS,
            sigma_range=6, scale_offset=True, verbose=False)
    values = np.log10(drawn._df[names].to_numpy())
    correlation = np.corrcoef(values.T)
    neighbours = float(np.mean([correlation[i, i + 1] for i in range(N - 1)]))
    return neighbours, float(values.std(axis=0).mean())


if __name__ == "__main__":
    print("nugget  before(corr, sd)   after(corr, sd)   exp(-h/a)/(1+n)")
    for nugget in (0.0, 0.5, 1.0):
        declared = ExpGeoStruct(range=RANGE, nugget=nugget, transform="log")
        before = measure(declared)
        after = measure(unit_sill(declared))
        print(f"{nugget:<6}  {before[0]:.3f}, {before[1]:.3f}       "
              f"{after[0]:.3f}, {after[1]:.3f}      {np.exp(-SPACING / RANGE) / (1 + nugget):.3f}")
