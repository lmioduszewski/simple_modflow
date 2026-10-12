"""The measurements behind the build-time activity guard (ledgers 321 and 322, design C5).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/active_cell_guard.py

Needs `mf6` on PATH (measured with 6.8.0.dev0). About 5 s.

The model: a 4 x 4 hand-built Voronoi grid (50 ft cells), a stack of fill 100-95 over
sand 95-65, and a pond cut to 93 ft through the fill at cells 5 and 6 -- so layer 1 is
REMOVED there (idomain 0) and the sand below is the connected cell. NEWTON, K 10, CHD
98 ft on the east edge in the sand, IMS closed tight. myflopy's guard would refuse the
cases below at build, so each is made by adding the FloPy package to an already-built
model directly -- which is what the guard is for.

PART 1 -- what FloPy and MF6 do with a cell the DISV removed. FloPy 3.10 ACCEPTS a CHD
record and a HEAD observation at the dug cell (0, 5) at construction; MF6 refuses both
at read ("Cell is outside active grid domain" / "Error reading data from ID string").
FloPy REFUSES an HFB touching the same cell at construction ("Cellid (0, 5) is outside
of the active model grid") -- while MF6 itself only warns about one and drops it
(read in ``gwf-hfb.f90``, not run: FloPy will not write one).

PART 2 -- what MF6 does with a recharge index naming a removed cell: drops it, silently.
With ``irch`` all layer 1 (what ``irch=None`` means, and what a ``"top_active"``
resolved before the cut writes), MF6 terminates normally with the budget closed and
applies the recharge of every column EXCEPT the two dug ones -- the RCH inflow in its
budget file is the declared total less exactly rate x the dug cells' area.

PART 3 -- the cost of the list walk the guard adds: a 654,480-record list (72 periods x
9,090 cells, the size of the imported Ten Trails ETS) against a 3 x 9,090 idomain.
"""

from __future__ import annotations

import atexit
import pathlib
import shutil
import subprocess
import tempfile
import time
import warnings
from types import SimpleNamespace

import flopy
import numpy as np

import myflopy as mf
from myflopy.grid_identity import _cellids_in, _inactive_cellid_complaints
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.specs import ModelContext

warnings.simplefilter("ignore")
WORK = pathlib.Path(tempfile.mkdtemp(prefix="active_cell_guard_"))
atexit.register(shutil.rmtree, WORK, ignore_errors=True)  # /tmp is a quota'd tmpfs (ledger 460)
RATE = 1e-3


def _base(tag):
    """The excavated model, built by myflopy, returned as a live FloPy simulation."""

    grid = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=50.0)
    pond = mf.Area(grid.gdf_vorPolys.geometry.loc[[5, 6]].union_all())
    result = (mf.LayerStack(top=mf.Flat(100.0)).add("fill", thickness=5.0)
              .add("sand", thickness=30.0)
              .excavate(bottom=93.0, within=pond, name="pond").build(grid, attach=True))
    east = [int(c) for c in np.flatnonzero(grid.gdf_vorPolys.geometry.bounds.maxx >= 200 - 1e-6)]
    model = mf.gwf("m", context=ModelContext(grid=grid, domain=result.idomain),
                   newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, icelltype=1), mf.ic(strt=99.0),
        mf.chd(stress_period_data={0: [[(1, c), 98.0] for c in east]}),
        mf.oc(budget_filerecord="m.cbc", saverecord=[("BUDGET", "ALL")]),
    ])
    spec = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-6, inner_dvclose=1e-8,
               rcloserecord=[1e-4, "STRICT"]),
    ])
    built = spec.build_flopy(WORK / tag)
    return built.simulation, built.simulation.get_model("m"), grid, result


def _run(sim):
    sim.write_simulation(silent=True)
    ws = pathlib.Path(sim.sim_path)
    subprocess.run(["mf6"], cwd=ws, capture_output=True, text=True, check=False)
    listing = (ws / "mfsim.lst").read_text()
    return ws, listing


def part1():
    sim, gwf, _, _ = _base("chd")
    flopy.mf6.ModflowGwfchd(gwf, pname="dug", filename="m.dug.chd",
                            stress_period_data={0: [[(0, 5), 98.0]]})
    print("PART 1  FloPy accepted a CHD at the dug cell (0, 5) at construction")
    _, listing = _run(sim)
    print("        MF6:", next(line.strip() for line in listing.splitlines()
                               if "outside active grid domain" in line))

    sim, gwf, _, _ = _base("obs")
    flopy.mf6.ModflowUtlobs(gwf, filename="m.obs",
                            continuous={"h.csv": [("h5", "HEAD", (0, 5))]})
    print("        FloPy accepted a HEAD observation at (0, 5) at construction")
    _, listing = _run(sim)
    print("        MF6:", next(line.strip() for line in listing.splitlines()
                               if "Error reading data from ID string" in line))

    _, gwf, _, _ = _base("hfb")
    try:
        flopy.mf6.ModflowGwfhfb(gwf, stress_period_data={0: [[(0, 4), (0, 5), 1e-3]]})
        print("        FloPy accepted an HFB touching (0, 5)")
    except flopy.mf6.mfbase.MFDataException as error:
        reason = next(line for line in str(error).splitlines() if "Cellid" in line)
        print("        FloPy refused an HFB touching (0, 5):", reason.strip())


def part2():
    sim, gwf, grid, _ = _base("rch")
    flopy.mf6.ModflowGwfrcha(gwf, recharge=RATE, irch=np.zeros(grid.ncpl, dtype=int),
                             save_flows=True)
    ws, listing = _run(sim)
    m_listing = (ws / "m.lst").read_text()
    discrepancy = float([line for line in m_listing.splitlines()
                         if "PERCENT DISCREPANCY" in line][-1].split()[-1])
    applied = float(np.sum(flopy.utils.CellBudgetFile(str(ws / "m.cbc"))
                           .get_data(text="RCH")[-1]["q"]))
    areas = grid.gdf_vorPolys.area.to_numpy()
    total = RATE * areas.sum()
    dug = RATE * areas[[5, 6]].sum()
    print(f"PART 2  Normal termination: {'Normal termination' in listing}; "
          f"discrepancy {discrepancy:.2f}%")
    print(f"        declared {total:.6f} ft3/d, dug columns {dug:.6f}; MF6 applied "
          f"{applied:.6f} = declared - dug to {abs(applied - (total - dug)) / total:.1e}")


def part3():
    ncpl, nlay, nper = 9090, 3, 72
    records = {p: [[(0, c), 100.0, RATE, 5.0] for c in range(ncpl)] for p in range(nper)}
    spec = mf.PackageSpec("evt", mf.evt.flopy(stress_period_data={0: [[(0, 0), 1.0, RATE, 5.0]]})
                          .builder, {"stress_period_data": records})
    idomain = np.ones((nlay, ncpl), dtype=int)
    idomain[0, :10] = 0
    model = SimpleNamespace(packagelist=[SimpleNamespace(name=["evt"], package_type="evt")])
    start = time.perf_counter()
    count = sum(1 for _ in _cellids_in(spec.options))
    walk = time.perf_counter() - start
    start = time.perf_counter()
    [complaint] = _inactive_cellid_complaints(model, [spec], idomain, {})
    guard = time.perf_counter() - start
    print(f"PART 3  {count:,} records: range walk {walk:.2f} s, activity check {guard:.2f} s")


if __name__ == "__main__":
    part1()
    part2()
    part3()
