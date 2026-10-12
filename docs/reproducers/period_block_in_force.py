"""Which period block MODFLOW 6 runs, and what the input views read, measured.

Ledger 449. Run it::

    python docs/reproducers/period_block_in_force.py [workspace]

It needs ``mf6`` on PATH and writes nothing outside the workspace it is given (a
temporary directory, removed at the end, by default). Under a minute.

WHAT IT MEASURES, each through MF6's own budget file
  1. A LIST package (GHB on a 3 x 3 model, blocks ``{0: two cells, 2: [], 3: one
     cell}``, five periods): MF6 runs the last block at or before each period and an
     EMPTY block runs none -- cells ``[2, 5], [2, 5], [], [6], [6]`` -- and
     ``model.ghb.inputs.get(per=N)`` now reads the same cells.
  2. HFB on the same model (blocks ``{0: three faces, 2: [], 3: one face}``): the faces
     ``model.hfb.results.q.get()`` lists at each step carry no flow, and a period-0
     face the table drops in period 2 does.
  3. UZF (the canonical testing model): a period-2 block listing ten objects at
     ``finf = 0`` leaves every other object infiltrating in period 3 -- a UZF block sets
     only the objects it lists -- and ``uzf.inputs.finf.get(per=3)`` reads that.
  4. Reading every list and UZF input view at every period of the canonical model,
     with GHB written once in period 0, adds no block: written and run, GHB applies 84
     records in every period (it was 84 -> 0 from the period read, before the fix).
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import flopy
import numpy as np

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _grid(n=3, d=10.0):
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _cells_mf6_ran(model, text):
    reader = model._get_budget_reader()
    return {
        int(kper): sorted(int(node) - 1 for record in reader.get_data(text=text, kstpkper=(kstp, kper))
                          for node in record["node"])
        for kstp, kper in reader.get_kstpkper()
    }


def small_model(root: Path):
    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, 9), dtype=int))
    hfb = {0: [[(0, 0), (0, 1), 1e-9], [(0, 3), (0, 4), 1e-9], [(0, 6), (0, 7), 1e-9]],
           2: [], 3: [[(0, 3), (0, 4), 1e-9]]}
    model = mf.gwf("flow", context=context, packages=[
        mf.disv(nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops()),
        mf.ic(strt=5.0), mf.npf(k=1.0),
        mf.chd(stress_period_data={0: [[(0, 0), 6.0], [(0, 8), 4.0]]}),
        mf.ghb(stress_period_data={0: [[(0, 2), 5.0, 1.0], [(0, 5), 5.0, 1.0]], 2: [],
                                   3: [[(0, 6), 5.5, 1.0]]}),
        mf.hfb.flopy(stress_period_data=hfb),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="blocks",
                               packages=[mf.tdis(nper=5, perioddata=[(1.0, 1, 1.0)] * 5)])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root / "small", name="small").prepare_run("r", simulation)
        ok, _ = run.execute(silent=True)
    assert ok
    model = run.model("flow")

    print("1. GHB: the cells MF6 ran, and what inputs.get(per=N) reads")
    print(f"   MF6:   {_cells_mf6_ran(model, 'GHB')}")
    print(f"   table: { {p: sorted(model.ghb.inputs.get(per=p)['cell']) for p in range(5)} }")

    print("2. HFB: the faces listed at each step, and their largest |q|")
    flows = model.hfb.results.q.get()
    print(f"   listed per period: {flows.groupby('per').size().to_dict()}; "
          f"max |q| on a listed face: {flows['q_cell1_to_cell2'].abs().max():.3g}")
    ia, ja = model.hfb.results._connectivity()
    face = next(k for k in range(ia[0], ia[1]) if ja[k] == 1)
    q = np.ravel(model._get_budget_reader().get_data(text="FLOW-JA-FACE", kstpkper=(0, 2))[0])[face]
    print(f"   the period-0 face (0, 1) in period 2, where no block holds it: q = {q:.4g}")


def _canonical(path):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return build_canonical_model(path, config=CanonicalModelConfig.testing())


def _run(model):
    with contextlib.redirect_stdout(io.StringIO()):
        model.sim.write_simulation(silent=True)
        ok, _ = model.sim.run_simulation(silent=True)
    assert ok


def uzf_rule(root: Path):
    print("3. UZF: a period-2 block sets ten objects to finf = 0")
    model = _canonical(root / "uzf")
    perioddata = model.package("uzf").perioddata
    first = perioddata.get_data(key=0).copy()
    for period in range(1, model.nper):
        perioddata.remove_transient_key(period)
    listed = first[:10].copy()
    listed["finf"] = 0.0
    perioddata.set_data({0: first, 2: listed})
    model.package("uzf").budget_filerecord.set_data([("uzf_rule.cbc",)])
    read = model.uzf.inputs.finf.get(per=3)
    _run(model)
    budget = flopy.utils.CellBudgetFile(str(Path(model.sim.sim_path) / "uzf_rule.cbc"))

    def infiltration(period):
        return {int(r["node"]) - 1: float(r["q"])
                for r in budget.get_data(text="INFILTRATION", kstpkper=(0, period))[0]}

    before, after = infiltration(0), infiltration(3)
    was = [i for i in range(10, len(first)) if before[i] > 0]
    still = [i for i in was if after[i] > 0]
    print(f"   MF6 period 3: objects 0-9 infiltrate {sum(after[i] for i in range(10))}; "
          f"{len(still)} of the {len(was)} unlisted objects infiltrating in period 0 still do")
    print(f"   finf.get(per=3): {len(read)} rows, objects 0-9 at "
          f"{sorted(set(read.set_index('ifno').loc[:9, 'finf']))}")


def reads_change_nothing(root: Path):
    print("4. every input view read at every period, then written and run")
    model = _canonical(root / "read")
    for period in range(1, model.nper):
        model.package("ghb").stress_period_data.remove_transient_key(period)
    for pkg in ("chd", "ghb", "rch", "evt", "wel", "drn", "riv"):
        for period in range(model.nper):
            getattr(model, pkg).inputs.get(per=period)
    for period in range(model.nper):
        model.uzf.inputs.finf.get(per=period)
    blocks = sorted(model.package("ghb").stress_period_data.get_active_key_dict())
    _run(model)
    ran = _cells_mf6_ran(model, "GHB")
    print(f"   GHB blocks after the reads: {blocks}; MF6's GHB records per period: "
          f"{[len(ran[p]) for p in range(model.nper)]}")


def main(root: Path) -> None:
    small_model(root)
    uzf_rule(root)
    reads_change_nothing(root)


if __name__ == "__main__":
    given = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    workspace = given or Path(tempfile.mkdtemp(prefix="period_block_in_force_"))
    try:
        main(workspace)
    finally:
        if given is None:
            shutil.rmtree(workspace, ignore_errors=True)
