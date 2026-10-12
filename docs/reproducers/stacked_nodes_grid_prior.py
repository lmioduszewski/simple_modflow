"""A correlated prior over two parameters at ONE location cannot be drawn, nugget or not.

Ledger 452. A layered DISU model's array file holds every node, so a
``style="grid"`` parameter with ``correlation=`` puts the nodes stacked in a column --
which share their centre -- in one correlated set. Run:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/stacked_nodes_grid_prior.py

THE SETUP. Three log-transformed parameters, two at (0, 0) and one at (5, 0), bounds
(0.5, 2.0), one exponential variogram of range 30 -- drawn by pyEMU's
``geostatistical_draws``, the function ``PstFrom.draw`` (and so ``run_ies``/``prior``)
reaches, once with no nugget and once with a nugget of 1.0.

WHAT IT SHOWS. Both draws fail with ``LinAlgError: Matrix is not positive definite``:
pyEMU builds the covariance (nugget on its diagonal), then sets each parameter's
variance from its bounds, which REPLACES the diagonal and so the nugget, leaving two
rows equal. This calls pyEMU with a structure of sill ``1 + nugget``, as the project
did until ledger 489; the project now scales it to a sill of 1 first
(``geostats.unit_sill``), and the nugget then survives and the draw succeeds -- but
it correlates the stacked nodes at ``1 / (1 + nugget)``, a vertical structure the
map-plane variogram never described, so the build still refuses the set
(``docs/reproducers/prior_nugget.py`` measures the nugget's effect).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyemu


def draw(nugget: float) -> str:
    pst = pyemu.Pst.from_par_obs_names(["p0", "p1", "p2"], ["o0"])
    par = pst.parameter_data
    par.loc[:, "partrans"] = "log"
    par.loc[:, "parval1"] = 1.0
    par.loc[:, "parlbnd"] = 0.5
    par.loc[:, "parubnd"] = 2.0
    par.loc[:, "pargp"] = "k"
    gs = pyemu.geostats.GeoStruct(
        variograms=pyemu.geostats.ExpVario(contribution=1.0, a=30.0), nugget=nugget,
        transform="log")
    frame = pd.DataFrame({"parnme": ["p0", "p1", "p2"], "x": [0.0, 0.0, 5.0],
                          "y": [0.0, 0.0, 0.0]})
    try:
        ensemble = pyemu.helpers.geostatistical_draws(pst, {gs: [frame]}, num_reals=5)
    except np.linalg.LinAlgError as error:
        return f"LinAlgError: {error}"
    return f"drew {ensemble.shape[0]} realizations"


if __name__ == "__main__":
    for nugget in (0.0, 1.0):
        print(f"nugget {nugget}: {draw(nugget)}")
