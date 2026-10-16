"""What pyEMU makes of a correlated LIST ``style="grid"`` set whose records share a point (ledger 570).

pyEMU locates each list record at its cell's centre. Where two parameters land on one
point it splits the covariance group by the index column before the cell's -- a DIS
cellid's layer -- and only when the key HAS one: a DISV ``(layer, cell)`` or a DISU
``node`` key has none, so the set stays one group with two parameters at one point.
This script turns myflopy's refusal OFF and shows what pyEMU then does:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/correlated_list_grid_shared_points.py

Needs ``mf6`` and pyEMU (no PEST++); writes only under a temporary directory.

THE SETUP. Two layers of 4 x 5 cells, 10 ft square (bottoms 50 and 0 ft), FloPy-written
and opened with ``mf.load_mf6_run``; CHD records in BOTH layers of the first three cells
("column"), or in layer 0 of five cells ("distinct"); ``chd`` with ``style="grid",
correlation=30.0``, built, then ``np.random.seed(0); draw_prior(400)``.

MEASURED (pyEMU 1.4.0, 2026-10-12)
  grid  case       covariance groups      prior draw
  DISV  distinct   chd                    correlated, off-diagonal 0.30 .. 0.74
  DISV  column     chd                    LinAlgError: Matrix is not positive definite
  DIS   distinct   chd                    correlated, 0.30 .. 0.74
  DIS   column     chd_cov0, chd_cov1     within a layer 0.44 .. 0.74, across -0.05 .. 0.00
  DISU  distinct   chd                    correlated, 0.30 .. 0.74
  DISU  column     chd                    LinAlgError: Matrix is not positive definite
With the refusal ON, both "column" sets on DISV and DISU are refused at build, and DIS
is built (its layers are each correlated).
"""

from __future__ import annotations

import pathlib
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
from flopy.utils.gridutil import get_disu_kwargs, get_disv_kwargs

import myflopy as mf
from myflopy.modflow.mf6.pest import native_parameters

NLAY, NROW, NCOL = 2, 4, 5
NCPL = NROW * NCOL
CASES = {
    "distinct": [(0, number) for number in range(5)],
    "column": [(layer, number) for number in range(3) for layer in (0, 1)],
}


def _run(workspace: pathlib.Path, grid: str, chd) -> pathlib.Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    delr, delc = np.full(NCOL, 10.0), np.full(NROW, 10.0)
    if grid == "disv":
        flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(NLAY, NROW, NCOL, delr, delc, 100.0,
                                                        [50.0, 0.0]))
        cell = lambda layer, number: (layer, number)  # noqa: E731
    elif grid == "dis":
        flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0,
                                top=100.0, botm=[50.0, 0.0])
        cell = lambda layer, number: (layer, number // NCOL, number % NCOL)  # noqa: E731
    else:
        flopy.mf6.ModflowGwfdisu(gwf, **get_disu_kwargs(NLAY, NROW, NCOL, delr, delc, 100.0,
                                                        [50.0, 0.0], return_vertices=True))
        cell = lambda layer, number: (layer * NCPL + number,)  # noqa: E731
    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[
        [cell(layer, number), 100.0 - 0.1 * number] for layer, number in chd])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


def main() -> None:
    refusal = native_parameters._refuse_list_records_pyemu_cannot_correlate
    native_parameters._refuse_list_records_pyemu_cannot_correlate = lambda *a, **k: None
    targets = mf.HeadTargets(
        locations=gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                                   geometry=gpd.points_from_xy([25.0], [25.0])),
        values=pd.DataFrame({"time": [0], "W": [99.0]}))
    try:
        with tempfile.TemporaryDirectory() as scratch:
            root = pathlib.Path(scratch)
            for grid in ("disv", "dis", "disu"):
                for case, cells in CASES.items():
                    run = _run(root / f"{grid}_{case}", grid, cells)
                    cal = mf.load_mf6_run(run).pest("c", workspace=root / f"pest_{grid}_{case}",
                                                    start_datetime="2024-01-01")
                    cal.observe(targets)
                    cal.parameterize("chd", style="grid", correlation=30.0)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        cal.build("c.pst")
                        frame = cal._native_parameter_frames["chd"]
                        try:
                            np.random.seed(0)
                            drawn = cal.draw_prior(400)._df[list(frame["parnme"])]
                        except np.linalg.LinAlgError as error:
                            print(f"{grid:5s} {case:9s} {sorted(set(frame['covgp']))}  "
                                  f"LinAlgError: {error}")
                            continue
                    corr = drawn.corr().to_numpy()
                    group = frame["covgp"].to_numpy()
                    pairs = ~np.eye(len(corr), dtype=bool)
                    same = pairs & (group[:, None] == group[None, :])
                    text = f"within a group {corr[same].min():.2f} .. {corr[same].max():.2f}"
                    if (pairs & ~same).any():
                        across = corr[pairs & ~same]
                        text += f", across groups {across.min():.2f} .. {across.max():.2f}"
                    print(f"{grid:5s} {case:9s} {sorted(set(group))}  {text}")
                    cal._release_pyemu_log()
    finally:
        native_parameters._refuse_list_records_pyemu_cannot_correlate = refusal


if __name__ == "__main__":
    main()
