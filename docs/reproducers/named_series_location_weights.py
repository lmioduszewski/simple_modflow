"""Lake and SFR location weights, a name given twice, and unchecked location weights.

Ledger 609. ``LakeStageTargets``, ``SfrStageTargets`` and ``SfrFlowTargets``
read a locations table as ``name -> id`` and dropped every other column, so the
``weight`` column heads, concentrations and DRN/RIV zones read was lost, and those
observations reached the control file at 1.0. A name given twice was resolved in
silence, differently per shape. No family checked a location's weight.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/named_series_location_weights.py

It uses only the public API, so it measures whichever tree is on ``PYTHONPATH`` --
the "before" column below is this script run against 55947911's ``src``
(``git archive 55947911 src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``).
Needs ``mf6`` and pyEMU; it builds a control file and never runs PEST++.

THE SETUP. ``tests/test_pest_handed_targets.py``'s model, built here: a 6 x 6
irregular Voronoi grid of ~100 ft cells, one confined layer, CHD at 102 ft west and
east, a four-reach stream, a two-cell lake ``pond`` (lake 0); three 10-day periods,
no MF6 OBS, so every family reads its FILE route.

  A. A lake target ``pond`` (location weight 50, a reading of weight 0.5 at period 0)
     and ``pond2`` (lake 0, no location weight) under
     ``LakeStageObservationSpec(weight=4.0)``; SFR stage ``g`` (location weight 50);
     SFR flow ``gage`` (location weight 0.01). The pst weight of each observation.
  B. A name given twice: each shape, and ``from_records`` naming two lakes.
  C. A location weight that is no number from 0: head -3, DRN zone ``'x'``.
  D. A lake ``from_records`` record carrying ``weight`` 50.

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, 2026-10-16)            before          after
  A  SFR stage g / SFR flow gage, location weight 50 / 0.01  1.0 / 1.0       50 / 0.01
     pond [period 0, 1, 2] (reading, location, location)    [0.5, 4, 4]     [0.5, 50, 50]
     pond2 [period 0, 1] (no location weight: the set's)     [4, 4]          [4, 4]
  B  lake DataFrame / records / Series repeat               kept LAST        REFUSED
     SFR DataFrame repeat                                   kept FIRST       REFUSED
     SFR Series repeat                                      kept LAST        REFUSED
     dict {1: 0, "1": 1}                                    kept lake 1      REFUSED
     from_records: 'a' on lakes 0 and 1                     lake 0           REFUSED
  C  head location weight -3                                pst -3.0         REFUSED
     DRN zone weight 'x'                                    pst 1.0          REFUSED
  D  from_records record weight 50                          dropped (NaN)    50.0

Why each is so, and what was decided: ledger 609.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
import pyemu

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

# The build writes a forward run that imports myflopy; made absolute here, as every
# reproducer that builds a PEST template does (ledger 425).
make_pythonpath_absolute()

REACHES = [13, 14, 15, 16]
LAKE_CELLS = [8, 9]
DRN_CELLS = [20, 21, 27]
PROBE = 123456.789


def build_toy(root: pathlib.Path):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    last = len(REACHES) - 1
    stream = mf.sfr.flopy(
        packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                      1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
        connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                        for i in range(len(REACHES))],
        perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]], 2: [[0, "INFLOW", 20.0]]},
        stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud",
    )
    lake = mf.lak.flopy(
        nlakes=1, packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
        connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                        for i, c in enumerate(LAKE_CELLS)],
        perioddata={0: [[0, "WITHDRAWAL", 10.0]]}, boundnames=True, stage_filerecord="m.lak.stg",
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        stream,
        mf.drn(stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(10.0, 3, 1.0)] * 3)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def build(model, workspace, observe):
    cal = model.pest("weights", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in observe:
        cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("weights.pst")
    return cal, pst


def weight_at(cal, pst, table, period, column):
    """The pst weight of the observation PEST reads from one cell, found through pyEMU."""

    template = cal.template_workspace
    frame = pd.read_csv(template / table)
    probe = frame.astype({column: float})
    probe.loc[frame.iloc[:, 0].astype(float) == period, column] = PROBE
    path = template / f"probe_{table}"
    probe.to_csv(path, index=False)
    read = pyemu.pst_utils.InstructionFile(str(template / f"{table}.ins")).read_output_file(str(path))["obsval"]
    name = read.index[np.isclose(read.to_numpy(dtype=float), PROBE)][0]
    return float(pst.observation_data.loc[name, "weight"])


def kept(make):
    """What a construction keeps, or REFUSED."""

    try:
        return make().get().to_dict("records")
    except ValueError as error:
        return f"REFUSED ({str(error)[:90]}...)"


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        model = build_toy(root)

        print("A. pst weights")
        pond = mf.LakeStageTargets(
            locations=pd.DataFrame({"name": ["pond", "pond2"], "lake": [0, 0], "weight": [50.0, None]}),
            values=pd.DataFrame({"time": [0, 1, 2, 0, 1], "name": ["pond", "pond", "pond", "pond2", "pond2"],
                                 "stage": [99.0, 98.0, 97.0, 99.0, 98.0],
                                 "weight": [0.5, None, None, None, None]}))
        stage = mf.SfrStageTargets(locations=pd.DataFrame({"name": ["g"], "reach": [1], "weight": [50.0]}),
                                   values=pd.DataFrame({"time": [0, 1], "g": [95.0, 95.5]}))
        flow = mf.SfrFlowTargets(locations=pd.DataFrame({"name": ["gage"], "reach": [3], "weight": [0.01]}),
                                 values=pd.DataFrame({"time": [0, 1], "gage": [-40.0, -41.0]}))
        cal, pst = build(model, root / "a", (LakeStageObservationSpec(targets=pond, weight=4.0), stage, flow))
        lake = "lak_stage_simulated_lake_stage.csv"
        print("   pond  ", [weight_at(cal, pst, lake, p, "pond") for p in (0, 1, 2)])
        print("   pond2 ", [weight_at(cal, pst, lake, p, "pond2") for p in (0, 1)])
        print("   g     ", [weight_at(cal, pst, "sfr_stage_simulated_sfr_stage.csv", p, "g") for p in (0, 1)])
        print("   gage  ", [weight_at(cal, pst, "sfr_flow_simulated_sfr_flow.csv", p, "gage") for p in (0, 1)])

        print("B. a name given twice")
        rows = {"lake": [{"name": "a", "lake": 0}, {"name": "a", "lake": 1}],
                "reach": [{"name": "a", "reach": 0}, {"name": "a", "reach": 1}]}
        print("   lake DataFrame ", kept(lambda: mf.LakeStageTargets(locations=pd.DataFrame(rows["lake"]))))
        print("   lake records   ", kept(lambda: mf.LakeStageTargets(locations=rows["lake"])))
        print("   lake Series    ", kept(lambda: mf.LakeStageTargets(locations=pd.Series([0, 1], index=["a", "a"]))))
        print("   lake dict      ", kept(lambda: mf.LakeStageTargets(locations={1: 0, "1": 1})))
        print("   SFR DataFrame  ", kept(lambda: mf.SfrStageTargets(locations=pd.DataFrame(rows["reach"]))))
        print("   SFR Series     ", kept(lambda: mf.SfrStageTargets(locations=pd.Series([0, 1], index=["a", "a"]))))
        print("   from_records   ", kept(lambda: mf.LakeStageTargets.from_records([
            {"name": "a", "lake": 0, "time": 0, "stage": 1.0}, {"name": "a", "lake": 1, "time": 1, "stage": 2.0}])))

        print("C. a location weight that is no number from 0")
        try:
            heads = mf.HeadTargets(locations=pd.DataFrame({"name": ["A"], "cell": [22], "weight": [-3.0]}),
                                   values=pd.DataFrame({"time": [0], "A": [95.0]}))
            cal, pst = build(model, root / "c1", (heads,))
            print("   head -3      pst", weight_at(cal, pst, "hds_simulated_heads.csv", 0, "A"))
        except ValueError as error:
            print("   head -3      REFUSED:", str(error)[:110])
        try:
            zone = mf.DrnFlowTargets(locations=pd.DataFrame({"name": ["z"], "cells": [DRN_CELLS], "weight": ["x"]}),
                                     values=pd.DataFrame({"time": [0], "z": [-3.0]}))
            cal, pst = build(model, root / "c2", (zone,))
            print("   DRN zone 'x' pst", weight_at(cal, pst, "drn_flow_simulated_drn_flow.csv", 0, "z"))
        except ValueError as error:
            print("   DRN zone 'x' REFUSED:", str(error)[:110])

        print("D. from_records keeps a record's weight")
        records = mf.LakeStageTargets.from_records([{"name": "a", "lake": 0, "time": 0, "stage": 1.0, "weight": 50.0}])
        print("   ", records.get().to_dict("records"))


if __name__ == "__main__":
    main()
