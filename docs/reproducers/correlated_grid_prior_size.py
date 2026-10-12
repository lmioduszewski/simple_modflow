"""What one correlated ``style="grid"`` block's prior draw costs, by its parameter count.

Ledger 522. ``cal.parameterize(target, style="grid", correlation=...)`` puts
one parameter per cell (per record of a list target) in one correlated block -- a
layer file of an array target -- and ``run_ies()``/``prior()``/``draw_prior()`` draw
its prior through pyEMU's ``helpers.geostatistical_draws``, which builds the block's
DENSE ``n x n`` covariance, scales it by the bounds' variance and decomposes it. On
one Ten Trails layer (20,323 cells) that is 3.3 GB of float64 before any copy, and the
draw raised ``MemoryError`` under a 12 GB address-space limit (finding scale-5). The
build now refuses a correlated block above ``CORRELATED_GRID_LIMIT`` parameters.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/correlated_grid_prior_size.py

Each size is drawn in its OWN process, so the peak resident memory read from
``resource.getrusage`` belongs to that draw alone. About a minute in all; the largest
draw holds ~2 GB, so run it on a machine with that free.

THE SETUP. ``n`` log-transformed parameters at random points in a 30,000 ft square
(State Plane offsets, as Ten Trails), bounds ``(0.1, 10)``, the project's geostruct
(``build_geostruct(unit_sill(ExpGeoStruct(range=3000)))``), 100 realizations, the
project's ``sigma_range=6`` and ``scale_offset=True``.

MEASURED (pyEMU 1.4.0, numpy on 20 cores shared with other work, 2026-10-12)
       n   one n x n   peak RSS     draw
   1,000     0.01 GB     0.40 GB      1.3 s
   2,000     0.03 GB     0.48 GB      1.3 s
   5,000     0.20 GB     0.99 GB      3.7 s
   8,000     0.51 GB     1.92 GB      7.7 s
  (load average ~9 at the start. Memory grows as n**2 -- 1.92 GB at 8,000 puts one
  20,323-cell layer near 12 GB, every layer drawn in turn at every launch -- so the
  limit sits at 5,000, a block's draw near 1 GB and a few seconds.)
"""

from __future__ import annotations

import contextlib
import io
import resource
import subprocess
import sys
import time

SIZES = (1000, 2000, 5000, 8000)


def draw(n: int) -> None:
    import numpy as np
    import pandas as pd
    import pyemu

    from myflopy.modflow.mf6.pest.geostats import build_geostruct, unit_sill
    from myflopy.modflow.mf6.pest.specs import ExpGeoStruct

    rng = np.random.default_rng(0)
    names = [f"p{i}" for i in range(n)]
    pst = pyemu.Pst.from_par_obs_names(names, ["o1"])
    par = pst.parameter_data
    par.loc[:, "partrans"] = "log"
    par.loc[:, "parval1"] = 1.0
    par.loc[:, "parlbnd"] = 0.1
    par.loc[:, "parubnd"] = 10.0
    frame = pd.DataFrame({
        "parnme": names,
        "x": rng.uniform(1.24e6, 1.24e6 + 30000.0, n),
        "y": rng.uniform(7.2e5, 7.2e5 + 30000.0, n),
    })
    structure = build_geostruct(unit_sill(ExpGeoStruct(range=3000.0, transform="log")))
    started = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        pyemu.helpers.geostatistical_draws(
            pst, struct_dict={structure: [frame]}, num_reals=100, sigma_range=6,
            scale_offset=True, subset=pd.Index(names), verbose=False,
        )
    seconds = time.perf_counter() - started
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # kB on Linux
    dense_gb = n * n * 8 / 1e9
    print(f"  {n:>6,}   {dense_gb:6.2f} GB   {peak_gb:6.2f} GB   {seconds:6.1f} s", flush=True)


def main() -> None:
    print("       n   one n x n   peak RSS     draw", flush=True)
    for n in SIZES:
        subprocess.run([sys.executable, __file__, str(n)], check=True)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        draw(int(sys.argv[1]))
    else:
        main()
