"""Pilot-point interpolation: log10 IDW from the 8 nearest points against global linear IDW.

Ledger 526 and 527. The forward run interpolated every cell from
EVERY pilot point of its layer by inverse distance squared, averaging the multipliers
in LINEAR space, and recomputed each cell's distances on every run; the net's spacing
was ``pp_space`` median cell widths. It now interpolates each cell from the 8 nearest
points of its set, in log10 space for a log parameter, through weights computed once
at build (``idw_factors``), and the spacing is ``pp_spacing=`` in model length units.
This measures both rules on one real mesh: the "before" column is the old rule
re-implemented here line for line (``forward_run.py`` at claude/pestprep 3a2d6f9),
the "after" column calls the library.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/pilot_point_interpolation.py [grid.disv.grb]

Reads a binary grid file READ-ONLY; the default is Ten Trails'
(``/home/lukem/models/mf6/tentrails/runs/hm_mounding_final_flow/tentrails.disv.grb``,
18,107 cells, 17,968 x 15,627 ft). About 30 s; the "old loop" timings are machine load.

THE SETUP. Nets laid by ``place_pilot_points(pp_spacing=...)`` on the grid's cell
centres. (1) Points per layer and cost per run of one layer. (2) A 500 ft net; every
point within R of the net's median point at 5, the rest at 1; the cells within 0.7 R.
(3) Every point at ``10 ** N(0, sd)``, five draws: the field's geometric mean over the
points'. sd 0.33 is the prior bounds (0.1, 10) give.

MEASURED (2026-10-12, load average ~12, Ten Trails ``hm_mounding_final_flow`` grid)
  (1) pp_spacing   points/layer   apply per run, after   old per-cell loop
        500 ft         772             1.0 ms              213 ms
       1000 ft         208             1.1 ms              234 ms
       1500 ft         100             1.4 ms              163 ms
       2500 ft          43             1.1 ms              208 ms
      (pp_space=5, the old default, laid 3,984 points per layer on this mesh; 4,501 on
      the 20,323-cell steady run's)
  (2) block            core cells, after (mean / min)   before (mean / min)
       R=1000 ft, 11 points      4.43 / 3.17                3.25 / 2.30
       R=1500 ft, 29 points      4.88 / 4.07                3.57 / 2.89
       R=3000 ft, 114 points     5.00 / 5.00                4.14 / 3.61
  (3) sd     field geometric mean / points', after   before
      0.25          1.010                              1.156
      0.33          0.974                              1.236
      0.50          0.966                              1.705
"""

from __future__ import annotations

import pathlib
import sys
import time
from types import SimpleNamespace

import numpy as np

from myflopy.modflow.mf6.pest.pilot_points import idw_factors, place_pilot_points

DEFAULT = pathlib.Path(
    "/home/lukem/models/mf6/tentrails/runs/hm_mounding_final_flow/tentrails.disv.grb")


def _old_idw(px, py, pv, x, y):
    """The old forward run's rule for one cell: global 1/d^2 over every point, linear."""

    d2 = (px - x) ** 2 + (py - y) ** 2
    zero = np.flatnonzero(d2 == 0)
    if zero.size:
        return float(pv[zero[0]])
    weights = 1.0 / np.maximum(d2, 1.0e-12)
    return float(np.sum(weights * pv) / np.sum(weights))


def main(path: pathlib.Path):
    import flopy

    grid = flopy.mf6.utils.MfGrdFile(str(path), verbose=False).modelgrid
    xc = np.asarray(grid.xcellcenters, dtype=float).reshape(-1)
    yc = np.asarray(grid.ycellcenters, dtype=float).reshape(-1)
    gwf = SimpleNamespace(modelgrid=grid)
    print(f"{xc.size} cells, {np.ptp(xc):,.0f} x {np.ptp(yc):,.0f}")

    for spacing in (500.0, 1000.0, 1500.0, 2500.0):
        pp = place_pilot_points(gwf, pp_spacing=spacing)
        px, py = pp["x"].to_numpy(), pp["y"].to_numpy()
        index, weight = idw_factors(xc, yc, px, py, 8)
        pv = 10 ** np.random.default_rng(0).normal(0, 0.33, len(pp))
        start = time.perf_counter()
        for _ in range(20):
            10.0 ** np.sum(weight * np.log10(pv)[index], axis=1)
        after = (time.perf_counter() - start) / 20
        start = time.perf_counter()
        for i in range(xc.size):
            _old_idw(px, py, pv, xc[i], yc[i])
        before = time.perf_counter() - start
        print(f"(1) pp_spacing={spacing:g}: {len(pp)} points; apply {after * 1e3:.1f} ms, "
              f"old loop {before * 1e3:.0f} ms")

    pp = place_pilot_points(gwf, pp_spacing=500.0)
    px, py = pp["x"].to_numpy(), pp["y"].to_numpy()
    index, weight = idw_factors(xc, yc, px, py, 8)
    centre = np.array([np.median(px), np.median(py)])
    for radius in (1000.0, 1500.0, 3000.0):
        pv = np.where(np.hypot(px - centre[0], py - centre[1]) < radius, 5.0, 1.0)
        core = np.flatnonzero(np.hypot(xc - centre[0], yc - centre[1]) < 0.7 * radius)
        after = (10.0 ** np.sum(weight * np.log10(pv)[index], axis=1))[core]
        before = np.array([_old_idw(px, py, pv, xc[i], yc[i]) for i in core])
        print(f"(2) R={radius:g} ({int((pv == 5).sum())} points): after {after.mean():.2f} / "
              f"{after.min():.2f}, before {before.mean():.2f} / {before.min():.2f}")

    rng = np.random.default_rng(1)
    for sd in (0.25, 0.33, 0.5):
        ratios_after, ratios_before = [], []
        for _ in range(5):
            pv = 10 ** rng.normal(0, sd, len(px))
            mean = 10 ** np.log10(pv).mean()
            sample = rng.choice(xc.size, 3000, replace=False)
            after = (10.0 ** np.sum(weight * np.log10(pv)[index], axis=1))[sample]
            before = np.array([_old_idw(px, py, pv, xc[i], yc[i]) for i in sample])
            ratios_after.append(10 ** np.log10(after).mean() / mean)
            ratios_before.append(10 ** np.log10(before).mean() / mean)
        print(f"(3) sd={sd}: after {np.mean(ratios_after):.3f}, before {np.mean(ratios_before):.3f}")


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT)
