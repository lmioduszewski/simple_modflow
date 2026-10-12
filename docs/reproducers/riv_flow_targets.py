"""RIV exchange-flow targets and the zone budget route: the MF6 facts they rest on.

Ledgers 534 .. -f. ``RivFlowTargets`` and ``DrnFlowTargets`` share one zone
implementation whose budget route reads records by the package that wrote them.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/riv_flow_targets.py [<copy of a Ten Trails run>]

Needs ``mf6``. Public API only, so it measures whichever tree is on ``PYTHONPATH``.

THE SETUP. The 6 x 6 irregular Voronoi grid of ``tests/test_flow_targets_period_end.py``
(``ncpl`` 36), two confined layers (top 110 ft, bottoms 50 and 0 ft), CHD at 102 ft on
the west and east edges of layer 0; three 10-day periods of one step, ``SAVE ... LAST``.
RIV ``riv`` (stage 95 ft) at cells 14, 15, 16 of layer 0, boundname ``gainers``; RIV
``crisp`` (stage 106 ft) at cells 20, 21 of layer 0, boundname ``losers``, absent from
period 2 on. RIV ``quiet`` (cell 9) and DRN ``drn_b`` (cell 22) are declared WITHOUT
``save_flows``; DRN ``drn_a`` (cell 27, layers 0 AND 1) with it. ``riv`` declares RIV
observations by cellid and by boundname; one names a cell with no river.

  A. MF6's RIV sign: aquifer -> river (a gaining reach) negative, river -> aquifer
     positive.
  B. A BOUNDNAME observation equals the sum of its cells' (why boundname zones were
     deferred: they change the API, not the numbers -- ledger 539).
  C. An observation on a cell holding no river: MF6 accepts it and writes 3e30.
  D. A package with no boundary in a period still writes a budget record, with zero
     entries; a package without ``save_flows`` writes none (ledger 536).
  E. ``RivFlowTargets`` budget route: every package, one package, the unsaved one.
  F. A zone drawn over ``drn_a``'s column with no layer: every layer (ledger
     535); a zone over ``crisp``'s cell 21 in period 2: 0.0 (ledger
     537).
  G. Optional: on a copy of a Ten Trails run (READ-ONLY source -- copy it first), the
     south DRN polygon with and without a layer, and ``riv`` / ``crisp`` apart.

MEASURED (2026-10-12, MF6 6.8, this tree):

  A. riv's cells -307.5, -304.2, -322.2 ft3/d at the end of period 0 (negative);
     crisp's +100.0 each (river leaking in).
  B. ``GALL - (G14 + G15 + G16)`` = 1.1e-13.
  C. ``NONE`` reads 3.0e30 in every row (the reader takes it as 0.0).
  D. One RIV and one CRISP record per saved time, CRISP with 0 entries in period 2; no
     QUIET record; DRN_A only, no DRN_B.
  E. package=None raises naming 'quiet' (RIV) and 'drn_b' (DRN) and ``save_flows=True``;
     package='riv' reads -307.49, -323.16, -325.43; package='crisp' 100.0, 100.0, 0.0.
     Before this commit the DRN zone read ``drn_a`` alone with no message (the design's
     toy: -339.52 where both drains removed ~679), and ``package=`` was refused on the
     budget route.
  F. The drawn zone with no layer: -1773.34 (both layers); with layer 0: -243.29 --
     what it read before this commit with no layer. crisp's cell 21: 100.0, 100.0, 0.0
     (it read NaN in period 2 before).
  G. Ten Trails ``ss`` copy: south DRN polygon -2,613,784.4 every layer, -33,300.5 in
     layer 0 (1.3%, before); every river cell -400,119.7 summed, ``riv`` -192,843.0,
     ``crisp`` -207,276.7.
"""

from __future__ import annotations

import sys
import tempfile
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid


def _model(root: Path):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = [int(c) for c in np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))]
    riv_obs = {"riv_obs.csv": [("g14", "RIV", (0, 14)), ("g15", "RIV", (0, 15)),
                               ("g16", "RIV", (0, 16)), ("gall", "RIV", "gainers"),
                               ("none", "RIV", (0, 30))]}
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[50.0, 0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 102.0] for c in edge]}),
        mf.riv.flopy(name="riv", filename="m.riv", save_flows=True, boundnames=True,
                     stress_period_data={0: [[(0, c), 95.0, 50.0, 93.0, "gainers"] for c in (14, 15, 16)]},
                     observations=riv_obs),
        mf.riv.flopy(name="crisp", filename="m.crisp", save_flows=True, boundnames=True,
                     stress_period_data={0: [[(0, c), 106.0, 50.0, 104.0, "losers"] for c in (20, 21)],
                                         2: []}),
        mf.riv.flopy(name="quiet", filename="m.quiet", save_flows=False,
                     stress_period_data={0: [[(0, 9), 95.0, 50.0, 93.0]]}),
        mf.drn.flopy(name="drn_a", filename="m.drn_a", save_flows=True,
                     stress_period_data={0: [[(0, 27), 96.0, 50.0], [(1, 27), 70.0, 50.0]]}),
        mf.drn.flopy(name="drn_b", filename="m.drn_b", save_flows=False,
                     stress_period_data={0: [[(0, 22), 95.0, 50.0]]}),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(10.0, 1, 1.0)] * 3)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def _try(label, call):
    try:
        print(f"  {label}: {call()}")
    except ValueError as error:
        print(f"  {label}: ValueError: {str(error)[:160]}")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        model = _model(Path(tmp))
        workspace = Path(model.workspace)
        observed = pd.read_csv(workspace / "riv_obs.csv")
        print("A/B/C. riv's RIV observations (cellid, boundname, a cell with no river):")
        print(observed.to_string(index=False))
        print("  GALL - (G14 + G15 + G16):",
              (observed["GALL"] - observed[["G14", "G15", "G16"]].sum(axis=1)).abs().max())

        cbc = CellBudgetFile(workspace / "m.cbc", precision="double")
        print("D. budget records (kper, term, package, entries):")
        for index, header in enumerate(cbc.recordarray):
            term = header["text"].decode().strip()
            if term in ("RIV", "DRN"):
                print(f"  {int(header['kper'])} {term} {header['paknam2'].decode().strip()} "
                      f"{len(cbc.get_record(index))}")
        cbc.close()

        print("E. RivFlowTargets / DrnFlowTargets budget route:")
        riv_zone = {"z": [14, 20]}
        _try("riv + crisp, package=None (quiet in scope)",
             lambda: mf.RivFlowTargets(locations=riv_zone).simulated_series(model)["z"].tolist())
        _try("package='riv'",
             lambda: mf.RivFlowTargets(locations=riv_zone, package="riv").simulated_series(model)["z"].tolist())
        _try("package='crisp'",
             lambda: mf.RivFlowTargets(locations=riv_zone, package="crisp").simulated_series(model)["z"].tolist())
        _try("DRN cells 22 + 27, package=None (drn_b unsaved)",
             lambda: mf.DrnFlowTargets(locations={"z": [22, 27]}).simulated_series(model)["z"].tolist())

        print("F. a drawn zone over drn_a's column, and crisp's cell 21 across its empty period:")
        polygon = model.vor.gdf_vorPolys.geometry.loc[27].buffer(-1.0)
        drawn = gpd.GeoDataFrame({"name": ["z"]}, geometry=[polygon], crs=model.vor.crs)
        _try("no layer", lambda: mf.DrnFlowTargets(locations=drawn, package="drn_a")
             .simulated_series(model)["z"].tolist())
        _try("layer 0", lambda: mf.DrnFlowTargets(locations=drawn.assign(layer=0), package="drn_a")
             .simulated_series(model)["z"].tolist())
        _try("crisp cell 21", lambda: mf.RivFlowTargets(locations={"c21": [21]}, package="crisp")
             .simulated_series(model)["c21"].tolist())

    if len(sys.argv) > 1:
        from collections import Counter

        copy = Path(sys.argv[1])
        print(f"G. {copy} (a COPY):")
        tt = mf.load_mf6_run(copy, verbosity_level=0)
        drn = tt.package("drn")
        cellids = drn.stress_period_data.get_data(0)["cellid"]
        print("  DRN records by layer:", sorted(Counter(int(c[0]) for c in cellids).items()))
        columns = sorted({int(c[1]) for c in cellids})
        south = tt.vor.gdf_vorPolys.geometry.loc[columns].union_all().buffer(-1.0)
        for layer in (None, 0):
            frame = {"name": ["south"]} if layer is None else {"name": ["south"], "layer": [layer]}
            zone = mf.DrnFlowTargets(locations=gpd.GeoDataFrame(frame, geometry=[south], crs=tt.vor.crs))
            print(f"  south polygon, layer {layer}: {float(zone.simulated_series(tt)['south'].iloc[-1]):,.1f}")
        rivers = tt.package("riv")
        cells = sorted({int(c[1]) for package in rivers for c in package.stress_period_data.get_data(0)["cellid"]})
        for package in (None, *(p.package_name for p in rivers)):
            series = mf.RivFlowTargets(locations={"all": cells}, package=package).simulated_series(tt)
            print(f"  every river cell, package={package!r}: {float(series['all'].iloc[-1]):,.1f}")


if __name__ == "__main__":
    main()
