"""The measurement behind ledger 469: what a build-side cellid check missed off a DISV.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/build_cellids_dis_disu.py

``ModelSpec.build``'s range and removed-cell checks read a DISV package only: on a DIS
or DISU model nothing was compared, so a cellid past the grid, on a removed cell, or of
another discretization's FORM built clean. Each case below is written by FloPy directly
(what ``ModelSpec.build`` handed it) on a 2 x 3 x 4 grid -- top 10, bottoms 5 and 0,
K 1, two CHD at 9 and 8 at opposite corners of layer 1 -- plus ONE extra CHD record at
9.0, and MF6 is run. Printed per case: whether FloPy wrote it, MF6's verdict, and for a
run that finished, the head MF6 held at the cell it actually read.

MEASURED (2026-10-11):

  DIS  (0, 0, 9)  column past ncol        FloPy writes '1 1 10 9.0'; MF6 stops
  DIS  (2, 0, 0)  layer past nlay         FloPy writes '3 1 1 9.0'; MF6 stops
  DIS  (0, 0, 1)  removed by IDOMAIN      FloPy writes '1 1 2 9.0'; MF6 stops
  DIS  (0, 5)     a DISV pair             FloPy stops: IndexError: tuple index out of range
  DISV (0, 1, 2)  a DIS triple            FloPy writes '1 2 9.0' -- it DROPS the third
                                          part -- and MF6 RUNS, holding (0, 1) at 9.0
  DISU (99,)      node past the grid      FloPy writes '100 9.0'; MF6 stops
  DISU (3, 5)     a DISV pair             FloPy writes '4 9.0' -- it drops the 5 -- and
                                          MF6 RUNS, holding node 3 at 9.0

Every one is now refused by ``ModelSpec.build`` naming the package, the cellid and the
grid's form or extent. The four MF6 stops on are named before it is run; the two that
RAN are why a wrong FORM is refused rather than left to FloPy, which keeps the parts it
expects and drops the rest without a word.
"""

from __future__ import annotations

import pathlib
import shutil
import tempfile

import flopy
import numpy as np
from flopy.utils.gridutil import get_disu_kwargs, get_disv_kwargs

NLAY, NROW, NCOL = 2, 3, 4


def _model(kind: str, extra, idomain=None):
    workspace = pathlib.Path(tempfile.mkdtemp())
    sim = flopy.mf6.MFSimulation(sim_name="s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="COMPLEX")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    if kind == "dis":
        flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0,
                                top=10.0, botm=[5.0, 0.0], idomain=idomain)
        corners = [(0, 0, 0), (1, NROW - 1, NCOL - 1)]
    elif kind == "disv":
        flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(NLAY, NROW, NCOL, 10.0, 10.0, 10.0, [5.0, 0.0]))
        corners = [(0, 0), (1, NROW * NCOL - 1)]
    else:
        kwargs = get_disu_kwargs(NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0),
                                 10.0, [5.0, 0.0], return_vertices=True)
        flopy.mf6.ModflowGwfdisu(gwf, **kwargs)
        corners = [(0,), (NLAY * NROW * NCOL - 1,)]
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=8.5)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[(corners[0], 9.0), (corners[1], 8.0), (extra, 9.0)])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    return sim, workspace


def case(kind: str, extra, read_back=None, idomain=None) -> None:
    label = f"{kind.upper():4s} {str(extra):10s}"
    try:
        sim, workspace = _model(kind, extra, idomain)
        sim.write_simulation(silent=True)
    except (IndexError, ValueError, flopy.mf6.mfbase.MFDataException) as error:
        print(f"{label} FloPy stops: {type(error).__name__}: {str(error).splitlines()[0][:70]}")
        return
    try:
        written = [line.strip() for line in (workspace / "m.chd").read_text().splitlines()
                   if line.strip() and line.strip()[0].isdigit()][-1]
        ok, _ = sim.run_simulation(silent=True)
        if not ok:
            print(f"{label} FloPy writes {written!r}; MF6 STOPS")
            return
        heads = flopy.utils.HeadFile(str(workspace / "m.hds")).get_data().reshape(-1)
        held = f", head at node {read_back}: {heads[read_back]:.4g}" if read_back is not None else ""
        print(f"{label} FloPy writes {written!r}; MF6 RUNS{held}")
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def main() -> None:
    idomain = np.ones((NLAY, NROW, NCOL), dtype=int)
    idomain[0, 0, 1] = 0
    case("dis", (0, 0, 9))
    case("dis", (2, 0, 0))
    case("dis", (0, 0, 1), idomain=idomain)
    case("dis", (0, 5))
    case("disv", (0, 1, 2), read_back=1)   # MF6 node of (layer 0, cell 1)
    case("disu", (99,))
    case("disu", (3, 5), read_back=3)


if __name__ == "__main__":
    main()
