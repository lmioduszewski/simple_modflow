"""What MF6 does with a NaN layer elevation, by the idomain of the cell it bounds (ledger 587).

A raster covering only part of a grid, used as a stack's top or a unit's bottom, left
NaN in a DISV's ``top``/``botm``. FloPy writes ``NAN``; MF6 reads it. This runs a 2-layer,
3 x 3 DISV with a NaN at cell 8 in four ways and compares every head with the same model
holding a finite value there.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/nan_layer_elevations.py

Needs ``mf6``; writes only under a temporary directory, which it removes.

MEASURED (mf6 6.8.0.dev0, FloPy 3.10, 2026-10-12)
  active      botm[0, 8] NaN, cell 8 active in both layers: normal termination, all 18
              heads NaN (MF6 reports success)
  passthrough botm[0, 8] NaN, cell (0, 8) idomain -1: normal termination, every head
              finite, but (1, 8) reads 90.0 -- its starting head -- against 91.67 with a
              finite botm: the cell below the pass-through silently disconnected
  inactive    top, botm[0, 8], botm[1, 8] NaN, column 8 idomain 0: normal termination,
              every active head equal to the finite model's
  below       botm[1, 8] NaN, (1, 8) idomain 0: as the finite model
So a NaN may sit only where the cell it bounds is INACTIVE: ``surfaces.refuse_unbounded_cells``
refuses one bounding any cell with ``idomain != 0``.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import flopy
import numpy as np
from flopy.utils.gridutil import get_disv_kwargs


def _run(ws: Path, *, nan=(), idomain_zero=(), passthrough=()):
    """Run the model with NaN at the (array, layer, cell) entries in ``nan``."""

    sim = flopy.mf6.MFSimulation("s", sim_ws=str(ws), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    kw = get_disv_kwargs(2, 3, 3, np.full(3, 10.0), np.full(3, 10.0), 100.0, [50.0, 0.0])
    top = np.full(9, 100.0)
    botm = np.array([np.full(9, 50.0), np.full(9, 0.0)])
    idomain = np.ones((2, 9), dtype=int)
    for which, layer, cell in nan:
        if which == "top":
            top[cell] = np.nan
        else:
            botm[layer, cell] = np.nan
    for layer, cell in idomain_zero:
        idomain[layer, cell] = 0
    for layer, cell in passthrough:
        idomain[layer, cell] = -1
    kw.update(top=top, botm=botm, idomain=idomain)
    flopy.mf6.ModflowGwfdisv(gwf, **kw)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=90.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 95.0], [(0, 4), 90.0]])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    heads = gwf.output.head().get_data().reshape(2, 9) if success else None
    return success, heads, idomain


def main() -> None:
    cases = {
        "active": dict(nan=[("botm", 0, 8)]),
        "passthrough": dict(nan=[("botm", 0, 8)], passthrough=[(0, 8)]),
        "inactive": dict(nan=[("top", 0, 8), ("botm", 0, 8), ("botm", 1, 8)],
                         idomain_zero=[(0, 8), (1, 8)]),
        "below": dict(nan=[("botm", 1, 8)], idomain_zero=[(1, 8)]),
    }
    with tempfile.TemporaryDirectory(prefix="nan_elevations_") as directory:
        root = Path(directory)
        for name, case in cases.items():
            finite = {key: value for key, value in case.items() if key != "nan"}
            _, reference, _ = _run(root / f"{name}_finite", **finite)
            success, heads, idomain = _run(root / name, **case)
            active = idomain > 0
            print(
                f"{name:<12} success={success}  NaN heads={int(np.isnan(heads[active]).sum())} of {int(active.sum())}  "
                f"max |head - finite model's| over active cells="
                f"{float(np.nanmax(np.abs(heads[active] - reference[active]))):.4f}  "
                f"head(1, 8)={heads[1, 8]:.4f} vs {reference[1, 8]:.4f}"
            )


if __name__ == "__main__":
    main()
