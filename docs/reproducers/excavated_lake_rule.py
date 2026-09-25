"""The measurements behind `mf.lak`'s excavated rule (ledger 325, the excavation design's C7).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/excavated_lake_rule.py

Needs `mf6` on PATH (measured with 6.8.0.dev0+0418efe4). About 15 s.

PART 1 -- a pit lake below the water table in uniform flow, against the same pit as a
K = 1e4 inclusion. A lake is water: in the saturated zone it should carry the regional
flow like an almost infinitely conductive inclusion, which is a reference MF6 can
compute with no LAK at all. `irregular_voronoi_grid` of 20 ft cells (40 x 40 and
60 x 60); a stack of 100 / 82 / 40 ft; a circular pit of radius 80 ft (65 and 66 lake
cells); CHD 97 ft on the west edge and 93 ft on the east, every active layer; NPF
k = k33 = 10, icelltype 1; NEWTON; steady; IMS COMPLEX, BICGSTAB, outer_dvclose 1e-8.
  R  no lake; K = 1e4 in the pit's layer-1 cells. Through-flow = what enters that set of
     cells from outside it (FLOW-JA-FACE through the model's own IA/JA).
  E  the pit EXCAVATED to 82 ft (`stack.excavate`) and `mf.lak` on it (bed leakance
     3000/d): the lake takes the excavated rule. Through-flow = GWF->lake inflow.
  L  the un-excavated stack and `mf.lak` with the same bed as a cell-indexed Series
     (`"bathy"`, the rule before ledger 325; its rim walls sit on the lake's own column).
Measured: E / R = 1.0054 (40 x 40) and 1.0053 (60 x 60); L / R = 0.5230 and 0.5250 --
0.3920 and 0.3978 until ledger 328, while L's CONNLEN was the full centre-to-centre
spacing (`grid.cl12` then; 1.52-3.16x, median 1.97, the lake cell's own side it is now).
GWF and LAK percent discrepancy 0.00 in every run -- asserted, never Normal
termination alone (the design measured a steady toy at Normal termination with a 53.89%
error). The design's own prototype measured 1.005-1.006 for this rule, 0.924-0.945 with
CONNLEN taken as the full centre-to-centre spacing (what `grid.cl12` held before ledger
328), 0.39-0.42 for L (with that spacing too).

PART 2 -- the canonical lake, which is NOT excavated (the design's F8), under the new
`unexcavated=` guard: 33 of its 33 VERTICAL connections sit in a cell whose top is above
the 96.0 ft floor, by up to 13.45 ft (median 9.669), all 33 tops above its starting
stage of 101. Its connection data is the rule before this change in every column but
CONNLEN (pinned in `tests/data/lak_legacy_golden.json`, written by that code); CONNLEN is
each lake cell's own side of its face since ledger 328.

NOT reproducible from the repository, and recorded in the ledger instead: the legacy rule
was also compared, row for row, on a 612-cell reconstruction of a real pond (1,753
connections, 1,141 horizontal, thinnest face 1.32e-4 ft) built from COPIES of that run's
MF6 files -- identical, before ledger 328 moved CONNLEN.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import footprint_cells


def _run(root, grid, result, k, lake=None):
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    west, east = np.flatnonzero(xy[:, 0] < xmin + 20), np.flatnonzero(xy[:, 0] > xmax - 20)
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    chd = [[(layer, int(c)), head] for cells, head in ((west, 97.0), (east, 93.0))
           for c in cells for layer in range(2) if result.idomain[layer, c] > 0]
    packages = [mf.disv(result), mf.npf(k=k, k33=k, icelltype=1), mf.ic(strt=95.0),
                mf.chd(stress_period_data={0: chd}),
                mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
                      saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])]
    if lake is not None:
        packages.append(lake(context))
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=packages)
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-8, inner_dvclose=1e-10,
               outer_maximum=300, inner_maximum=300, linear_acceleration="BICGSTAB")])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    workspace = pathlib.Path(run.workspace)
    listing = (workspace / "m.lst").read_text()
    discrepancy = [float(line.split()[-1]) for line in listing.splitlines()
                   if "PERCENT DISCREPANCY" in line]
    return workspace, discrepancy


def part_one(root: pathlib.Path) -> None:
    for n in (40, 60):
        grid = irregular_voronoi_grid(nrow=n, ncol=n, cell_size=20.0, refined=False)
        xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
        pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs,
                                geometry=[Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(80.0, 64)])
        stack = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(82.0))
                 .add("lower", bottom=mf.Flat(40.0)))
        plain = stack.build(grid)
        cells = footprint_cells(grid.gdf_vorPolys, pond.geometry)[0]
        k = np.full((2, grid.ncpl), 10.0)
        k_ref = k.copy()
        k_ref[0, cells] = 1.0e4

        workspace, d_ref = _run(root / f"R{n}", grid, plain, k_ref)
        grb = flopy.mf6.utils.MfGrdFile(str(workspace / "m.disv.grb"))
        ia, ja = np.asarray(grb.ia), np.asarray(grb.ja)
        flows = np.asarray(flopy.utils.CellBudgetFile(str(workspace / "m.cbc"), precision="double")
                           .get_data(text="FLOW-JA-FACE")[-1]).ravel()
        inside = set(cells)
        reference = sum(flows[p] for c in cells for p in range(ia[c] + 1, ia[c + 1])
                        if ja[p] not in inside and flows[p] > 0)

        def through(workspace):
            budget = flopy.utils.CellBudgetFile(str(workspace / "m_lak_budget.lak"), precision="double")
            q = np.asarray(budget.get_data(text="GWF")[-1]["q"])
            return float(q[q > 0].sum())

        pit = stack.excavate(bottom=82.0, within=pond, name="pit").build(grid)
        excavated = lambda ctx: mf.lak(context=ctx, nper=1, lakes=pond, lake_id_field="name",  # noqa: E731
                                       starting_stage=95.0, bed_leakance=3000.0, unexcavated="raise")
        workspace, d_e = _run(root / f"E{n}", grid, pit, k, excavated)
        e = through(workspace)
        legacy = lambda ctx: mf.lak(context=ctx, nper=1, lakes=pond, lake_id_field="name",  # noqa: E731
                                    starting_stage=95.0, bed_leakance=3000.0, unexcavated="ignore",
                                    lake_bottom=pd.Series({int(c): 82.0 for c in cells}))
        workspace, d_l = _run(root / f"L{n}", grid, plain, k, legacy)
        el = through(workspace)
        print(f"{n} x {n}, {len(cells)} lake cells: reference {reference:.3f} ft3/d (discrepancy {d_ref})")
        print(f"  excavated rule  E/R = {e / reference:.4f}  (GWF, LAK discrepancy {d_e})")
        print(f"  legacy 'bathy'  L/R = {el / reference:.4f}  (GWF, LAK discrepancy {d_l})")


def part_two(root: pathlib.Path) -> None:
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model

    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        build_canonical_model(root / "canonical", config=CanonicalModelConfig.testing())
    guard = [str(w.message) for w in caught if "NOT excavated" in str(w.message)]
    print("canonical lake, the unexcavated= guard:")
    print("  " + guard[0].splitlines()[1] if guard else "  (no warning)")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        part_one(pathlib.Path(tmp))
        part_two(pathlib.Path(tmp))
