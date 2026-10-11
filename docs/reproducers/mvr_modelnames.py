"""What MF6 does with ``MODELNAMES`` on a MODEL's mover (ledger 413).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/mvr_modelnames.py

``mvr_spec`` (``mf.mvr``, ``mf.mvr.flopy``) builds the ``MVR`` package of one GWF
model. MF6's ``MODELNAMES`` option makes every record name a model before each
package -- ``(mname, pname)`` in ``PACKAGES``, ``mname1 pname1 id1 mname2 pname2 id2
mvrtype value`` per period -- and is meant for the mover of a GWF-GWF exchange.
Before ledger 413, ``mf.mvr(modelnames=True)`` wrote ``MODELNAMES`` beside
package-only records, and ``mf.mvr.flopy(..., modelnames=True)`` read the
model-name layout as the single-model one (its MODEL name became the package the
spec ``requires=``) and died on a misleading message. Both are refused now.

THE SETUP. The model of ``docs/reproducers/sfr_stream_inputs.py``: a 6 x 6 irregular
Voronoi grid of ~100 ft cells, one confined layer, CHD at 90 ft on the west and east
edges, one steady day, two SFR streams (``source``, 4 reaches, 100 ft3/d entering its
headwater; ``canal``, 2 reaches, unconnected), and a mover sending FACTOR 0.5 of the
source's last reach (reach 3) to the canal's first (reach 4).

WHAT IT MEASURES. The mover file as ``mf.mvr`` writes it, then the same file rewritten
three ways, each run by MF6: (1) as written; (2) ``MODELNAMES`` added, package-only
records (what ``mf.mvr(modelnames=True)`` wrote); (3) ``MODELNAMES`` added, records in
the model-name layout. Measured on MF6 6.8: (1) runs and the canal's first reach
receives 49.998 ft3/d FROM-MVR; (2) and (3) both stop reading ``m.mvr`` with
"MODELNAMES cannot be specified unless the mover package is for an exchange".
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from sfr_stream_inputs import _builder  # noqa: E402

import myflopy as mf  # noqa: E402


def _model(root):
    builder = _builder(mover=True)
    grid = builder.grid
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < 100.0) | (xy[:, 0] > 500.0))
    observations = [(f"from_mvr_{reach}", "FROM-MVR", reach + 1) for reach in range(builder.total_nreaches)]
    sfr = builder.with_updates(options={"observations": {"sfr.obs.csv": observations}}).build()
    mvr = mf.mvr(nper=1, moves=(mf.Move(mf.MoverConnection("sfr", 3), mf.MoverConnection("sfr", 4), value=0.5),))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        sfr, mvr,
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").run("r", simulation)
    return pathlib.Path(run.model().workspace)


def _mf6(workspace):
    (workspace / "sfr.obs.csv").unlink(missing_ok=True)
    result = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
    lines = (workspace / "mfsim.lst").read_text().splitlines()
    start = next((i for i, line in enumerate(lines) if "ERROR REPORT" in line), None)
    if start is not None:
        return "MF6 refused: " + " ".join(line.strip() for line in lines[start + 1 : start + 4] if line.strip())
    if result.returncode != 0:
        return f"MF6 exited {result.returncode} with no error report"
    flows = pd.read_csv(workspace / "sfr.obs.csv").iloc[-1]
    return f"runs: FROM-MVR into reach 4 = {flows['FROM_MVR_4']:.3f} ft3/d"


def main() -> None:
    root = pathlib.Path(tempfile.mkdtemp(prefix="mvr_modelnames_"))
    try:
        workspace = _model(root)
        written = (workspace / "m.mvr").read_text()
        with_option = written.replace("BEGIN options\n", "BEGIN options\n  MODELNAMES\n")
        layouts = {
            "1 as mf.mvr writes it": written,
            "2 MODELNAMES, package-only records": with_option,
            "3 MODELNAMES, model-name layout": with_option.replace("\n  sfr\n", "\n  m  sfr\n").replace(
                "  sfr  4  sfr  5", "  m  sfr  4  m  sfr  5"),
        }
        for label, text in layouts.items():
            (workspace / "m.mvr").write_text(text)
            block = text[text.index("BEGIN packages"):].split("\n")
            print(f"{label}: {' | '.join(line.strip() for line in block if line.strip() and 'BEGIN' not in line and 'END' not in line)}")
            print(f"   {_mf6(workspace)}")
        print("\nmyflopy now:")
        for label, call in {
            "mf.mvr(modelnames=True)": lambda: mf.mvr(nper=1, modelnames=True, moves=(
                mf.Move(mf.MoverConnection("sfr", 3), mf.MoverConnection("sfr", 4)),)),
            "mf.mvr.flopy(model-name layout, modelnames=True)": lambda: mf.mvr.flopy(
                packages=[["m", "sfr"]], perioddata={0: [("m", "sfr", 3, "m", "sfr", 4, "FACTOR", 0.5)]},
                modelnames=True),
        }.items():
            try:
                call()
                print(f"   {label}: built")
            except ValueError as error:
                print(f"   {label}: refused: {str(error)[:110]}...")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()
