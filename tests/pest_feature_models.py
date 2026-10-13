"""A small MF6 model with lakes and streams, for the feature-target PEST tests.

Shared by ``test_pest_lak_bedleak.py`` and ``test_pest_sfr_rhk.py`` (ledger
558). DISV, 2 layers of 4 x 5 cells, two periods, a CHD at each end:

* LAK ``lak``: lake 0 ``Keevie`` -- three connections, TWO of them to cell (0, 7),
  bedleak 0.1, 0.2, 0.25 -- and lake 1 ``marjorie``, two connections, 0.3 and 0.35.
* SFR ``sfr``: reaches 0-2 ``creek_a`` and 3-4 ``Creek_B``; reaches 1 and 2 share cell
  (0, 11). rhk 0.5 .. 0.9, Manning's n 0.030 .. 0.034.

Variants drop the boundnames, put ``NONE`` in lake 1's second bedleak or in reach 4's
cellid -- both of which MF6 runs --, leave reach 4 without a boundname, or add
a second LAK package ``ponds``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

NLAY, NROW, NCOL = 2, 4, 5
#: Lake connections: (lake, connection, cell, claktype, bedleak).
LAKE_CONNECTIONS = [
    (0, 0, 6, "vertical", 0.1),
    (0, 1, 7, "horizontal", 0.2),
    (0, 2, 7, "horizontal", 0.25),
    (1, 0, 13, "vertical", 0.3),
    (1, 1, 14, "vertical", 0.35),
]
#: Reaches: (cell, boundname, rhk, man).
REACHES = [
    (10, "creek_a", 0.5, 0.030),
    (11, "creek_a", 0.6, 0.031),
    (11, "creek_a", 0.7, 0.032),
    (16, "Creek_B", 0.8, 0.033),
    (17, "Creek_B", 0.9, 0.034),
]


def build_model(workspace: Path, *, lake_boundnames=True, sfr_boundnames=True,
                none_bedleak=False, none_reach=False, two_lak=False,
                unnamed_reach=False) -> Path:
    """Write and run the model in ``workspace``; return it."""

    import flopy
    from flopy.utils.gridutil import get_disv_kwargs

    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(
        NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [50.0, 0.0]))
    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 100.0], [(0, 4), 99.0]])
    connections = []
    for lake, number, cell, kind, bedleak in LAKE_CONNECTIONS:
        if none_bedleak and (lake, number) == (1, 1):
            bedleak = "NONE"
        geometry = (0.0, 0.0, 0.0, 0.0) if kind == "vertical" else (90.0, 100.0, 5.0, 10.0)
        connections.append([lake, number, (0, cell), kind, bedleak, *geometry])
    names = ["Keevie", "marjorie"]
    packagedata = [[lake, 101.0, sum(row[0] == lake for row in LAKE_CONNECTIONS)]
                   + ([names[lake]] if lake_boundnames else []) for lake in range(2)]
    flopy.mf6.ModflowGwflak(gwf, pname="lak", filename="m.lak", boundnames=lake_boundnames,
                            nlakes=2, noutlets=0, ntables=0, packagedata=packagedata,
                            connectiondata=connections,
                            perioddata={0: [[0, "status", "active"], [1, "status", "active"]]})
    if two_lak:
        flopy.mf6.ModflowGwflak(gwf, pname="ponds", filename="m.ponds", nlakes=1,
                                noutlets=0, ntables=0, packagedata=[[0, 101.0, 1]],
                                connectiondata=[[0, 0, (0, 2), "vertical", 0.5,
                                                 0.0, 0.0, 0.0, 0.0]],
                                perioddata={0: [[0, "status", "active"]]})
    reaches = []
    downstream = {0: [-1], 1: [0, -2], 2: [1], 3: [-4], 4: [3]}
    for number, (cell, bound, rhk, man) in enumerate(REACHES):
        cellid = "none" if none_reach and number == 4 else (0, cell)
        reaches.append([number, cellid, 10.0, 2.0, 0.001, 95.0 - number * 0.1, 1.0, rhk, man,
                        len(downstream[number]), 1.0, 0]
                       + ([bound] if sfr_boundnames and not (unnamed_reach and number == 4)
                          else []))
    flopy.mf6.ModflowGwfsfr(gwf, pname="sfr", filename="m.sfr", boundnames=sfr_boundnames,
                            nreaches=len(REACHES), packagedata=reaches,
                            connectiondata=[[number, *links]
                                            for number, links in downstream.items()],
                            perioddata={0: [[0, "inflow", 10.0], [3, "inflow", 5.0]]})
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0], (workspace / "mfsim.lst").read_text()[-2000:]
    return workspace


def head_targets():
    """One head target at cell 2, both periods."""

    import myflopy as mf

    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([25.0], [35.0]))
    return mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [0, 1], "W": [99.6, 99.6]}))


def calibration(run: Path, tmp_path: Path, name: str = "feat"):
    """A calibration of the run at ``run`` observing :func:`head_targets`."""

    import myflopy as mf

    model = mf.load_mf6_run(run)
    cal = model.pest(name, workspace=tmp_path / name, start_datetime="2024-01-01")
    cal.observe(head_targets())
    return cal


def build(cal, name: str = "feat"):
    """``cal.build`` at NOPTMAX 0, pyEMU's warnings silenced."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst", noptmax=0)


def forward_run(cal, pst, tmp_path: Path, values: dict[str, float]) -> Path:
    """Run the template's forward run in a copy, as a PEST worker does, each parameter
    GROUP in ``values`` set to its value; assert MF6 terminates normally."""

    import pyemu

    worker = tmp_path / "worker"
    shutil.rmtree(worker, ignore_errors=True)
    shutil.copytree(cal.template_workspace, worker)
    par = pst.parameter_data.copy()
    for group, value in values.items():
        names = par.index[par["pargp"] == group]
        assert len(names), f"no parameter of group {group!r}: {sorted(set(par['pargp']))}"
        par.loc[names, "parval1"] = value
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(worker / tpl), str(worker / inp))
    done = subprocess.run([sys.executable, "forward_run.py"], cwd=worker, capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    assert "Normal termination" in (worker / "mfsim.lst").read_text()
    return worker


def records(path: Path) -> pd.DataFrame:
    """An MF6 list file as text columns."""

    return pd.read_csv(path, sep=r"\s+", header=None, dtype=str, keep_default_na=False)
