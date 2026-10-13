"""A model-like object whose heads are a REAL MF6 head file, for the target readers (ledger 581).

The head targets read a run model's values through ONE reader -- the binary file
its OC declares, at each period's end (ledger 557). A test that wants chosen heads
used to hand ``compare()`` a ``SimpleNamespace`` carrying an ``all_heads`` table,
which kept a second derivation alive in the library for it alone. Here the heads
are written by MF6 itself: every cell is a CHD at the head asked for, so the file
holds those numbers exactly, and the double carries only what the reader asks of
a run model -- its ``workspace``, ``name``, grid and TDIS lengths.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import flopy
import numpy as np


def head_file_model(workspace: Path, vor, heads, *, name: str = "m", perlen: float = 1.0):
    """Run MF6 so that ``<name>.hds`` holds ``heads``; a model-like object reading it.

    Parameters
    ----------
    workspace : pathlib.Path
        Where the simulation is written and run (created).
    vor : VoronoiGridPlus
        The grid; its DISV properties are written, one cell per polygon.
    heads : sequence of array-like
        One ``(nlay, ncpl)`` array per stress period (a ``(ncpl,)`` array is one
        layer): every cell is a CHD at that head, so MF6 writes it exactly.
    name : str, optional
        The GWF model's name. Default ``"m"``.
    perlen : float, optional
        Each period's length (one step each). Default 1.0.

    Returns
    -------
    types.SimpleNamespace
        ``vor``, ``workspace``, ``name``, ``sim`` (the FloPy simulation, whose TDIS
        the period-end rule reads) and ``gwf`` (the FloPy model, with its DISV).
    """

    workspace = Path(workspace)
    workspace.mkdir(parents=True, exist_ok=True)
    frames = [np.atleast_2d(np.asarray(period, dtype=float)) for period in heads]
    nlay, ncpl = frames[0].shape
    props = vor.get_disv_gridprops()
    top = float(max(frame.max() for frame in frames)) + 10.0
    bottoms = [float(min(frame.min() for frame in frames)) - 10.0 * (k + 1) for k in range(nlay)]
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=len(frames), perioddata=[(perlen, 1, 1.0)] * len(frames))
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdisv(
        gwf, nlay=nlay, ncpl=ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
        cell2d=props["cell2d"], top=top, botm=bottoms, idomain=np.ones((nlay, ncpl), dtype=int),
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=float(frames[0].mean()))
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={
        per: [[(k, c), float(frame[k, c])] for k in range(nlay) for c in range(ncpl)]
        for per, frame in enumerate(frames)
    })
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")])
    sim.write_simulation(silent=True)
    success, report = sim.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])
    return SimpleNamespace(vor=vor, workspace=workspace, name=name, sim=sim, gwf=gwf)
