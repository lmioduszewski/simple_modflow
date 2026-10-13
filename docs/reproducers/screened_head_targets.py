"""Screened head targets on a real NEWTON model, and what one shared period-end reader costs.

Ledgers 556 (a well screened across layers is ONE observation of
its K x saturated-open-length weighted head) and 557 (``compare()``
and the PEST forward run read through one period-end reader of the OC-declared
head file).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/screened_head_targets.py <run copy> [n_targets]

``<run copy>`` is a COPY of a run directory (the user's Ten Trails runs are read
only); only its input files and its head file are read. It prints:

1. how many ACTIVE cells hold a head below their own bottom at the last saved
   time, and how many hold MF6's -1e30 -- under NEWTON a dry cell keeps a finite
   head, so a dry layer has to be weighted out by its saturated thickness;
2. two real columns with a screen from the middle of layer 0 to the middle of
   layer 2 -- one whose water table is inside layer 0, one whose layers 0-2 all
   hold a head below their bottoms -- each layer's head, K, open length and
   saturated open length, and the screened head ``compare()`` reports;
3. the old ``compare()`` route -- the ``HeadsPlus`` value table, period-end rows,
   joined -- against the shared reader, for ``n_targets`` (default 200) random
   layer-0 targets at every period: wall time, and whether every value agrees.

MEASURED 2026-10-12 (load average 1.5-3 on 20 cores):

  Ten Trails ``ss`` copy (1 period, 5 x 20,323 DISV, NEWTON, ICELLTYPE 1):
    1. 101,615 active cells; 16,900 hold a finite head below their bottom
       (8,065 / 7,843 / 992 / 0 / 0 by layer), none -1e30.
    2. cell 13035, screen 455.76-535.98: layer 0 head 528.007, K 374.5, open 8.02,
       saturated 0.05; layer 1 head 476.664, K 0.4427, open 52.37, saturated 1.08;
       layer 2 head 418.508 (below its bottom 435.932), open 19.82, saturated 0 --
       screened head 526.820. Cell 85, all three below their bottoms: 406.612,
       layer 2's head (the deepest open layer).
    3. 200 targets: old route 0.62-1.46 s, shared reader 0.00-0.04 s, identical.
  Ten Trails ``hm_mounding_final_flow`` copy (38 periods, 5 x 18,107):
    1. 90,446 active cells; at the last period 12,185 hold a finite head below
       their bottom (4,946 / 6,303 / 936 / 0 / 0), none -1e30.
    2. cell 15061, screen 456.06-539.75: layer 0 head 531.412 (saturated 0.25 of
       8.60, K 411.6), layer 1 479.492 (4.81 of 56.48, K 0.436), layer 2 416.217
       (dry) -- screened head 530.393. Cell 85, all dry: 404.202, layer 2's head.
    3. 200 targets x 38 periods: old route 3.74-6.82 s, shared reader 0.02-0.05 s,
       every value identical.

  ``compare()`` with 200 layer-0 targets at 38 periods, measured with a separate
  script in two processes (the code at HEAD before the change, then after):
  3.2-4.2 s per call and 1.89 GB peak RSS before (7.0-9.5 s at a load of ~14);
  0.91 s on the first call -- which loads the simulation for its TDIS -- then
  0.06 s, and 0.45 GB, after; the same 7,600 values.
"""

from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.hydraulics import model_flow_properties, saturated_thickness
from myflopy.modflow.mf6.output_files import _declared_output_path
from myflopy.modflow.mf6.period_end import at_period_end, period_ends, period_lengths
from myflopy.modflow.mf6.target_values import read_target_values


def _old_compare_route(model) -> pd.DataFrame:
    """The derivation ``compare()`` used before ledger 557, kept HERE as part 3's baseline.

    It was ``observations._shared._simulated_field_by_period``, which the library
    kept for test doubles alone until ledger 581 deleted it: the
    ``HeadsPlus`` value table, each period's END row on the file's own clock.
    """

    frame = model.all_heads.reset_index().copy()
    ends = period_ends(model, reader=lambda: model.hds, file="the head file", saves="HEAD",
                       fallback=frame["kstpkper"])
    frame = at_period_end(frame, ends).copy()
    frame["per"] = frame["kstpkper"].map(lambda item: int(item[1]))
    frame = frame.sort_values("per", kind="stable").rename(columns={"elev": "sim_head"})
    return frame.loc[:, ["per", "layer", "cell", "sim_head"]].reset_index(drop=True)


def main():
    workspace = Path(sys.argv[1]).resolve()
    count = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    warnings.simplefilter("ignore")
    model = mf.load_mf6_run(workspace, verbosity_level=0)
    flow = model_flow_properties(model)
    nlay = np.asarray(flow.botm).shape[0]
    botm = np.asarray(flow.botm, dtype=float).reshape(nlay, -1)
    top = np.asarray(flow.top, dtype=float).reshape(-1)
    tops = np.vstack([top[None, :], botm[:-1]])
    k = np.asarray(flow.k, dtype=float).reshape(nlay, -1)
    idomain = np.asarray(flow.idomain).reshape(nlay, -1)
    celltype = np.broadcast_to(np.asarray(flow.icelltype).reshape(nlay, -1), botm.shape)
    head_file = _declared_output_path(workspace, model.name, "head")
    import flopy

    hds = flopy.utils.HeadFile(str(head_file))
    last = hds.get_kstpkper()[-1]
    heads = np.asarray(hds.get_data(kstpkper=last), dtype=float).reshape(nlay, -1)
    hds.close()

    active = idomain > 0
    below = active & (heads > -1e29) & (heads < botm)
    print(f"1. {workspace.name}: {int(active.sum())} active cells; at (kstp, kper) "
          f"{tuple(int(v) for v in last)} "
          f"{int(below.sum())} hold a finite head below their bottom "
          f"(per layer {below.sum(axis=1).tolist()}), "
          f"{int((active & (heads <= -1e29)).sum())} hold -1e30.")

    three = active[0] & active[1] & active[2]
    wet0 = three & (heads[0] > botm[0]) & (heads[0] < tops[0])
    dry0 = three & below[0] & below[1] & below[2]
    picks = []
    if wet0.any():  # layer 0 partly saturated: the water table inside it
        picks.append(("water table in layer 0", int(np.flatnonzero(wet0)[
            np.argmax(np.abs(heads[0] - heads[2])[wet0])])))
    if dry0.any():  # every open layer dry: the deepest open layer's head
        picks.append(("layers 0-2 all below their bottoms", int(np.flatnonzero(dry0)[0])))
    for label, cell in picks:
        screen_top = 0.5 * (tops[0, cell] + botm[0, cell])
        screen_bottom = 0.5 * (tops[2, cell] + botm[2, cell])
        print(f"2. cell {cell} ({label}): screen {screen_bottom:.2f} to {screen_top:.2f}")
        for layer in range(3):
            upper = min(screen_top, tops[layer, cell])
            lower = max(screen_bottom, botm[layer, cell])
            saturated = float(saturated_thickness(heads[layer, cell], upper, lower,
                                                  icelltype=celltype[layer, cell]))
            print(f"   layer {layer}: head {heads[layer, cell]:.3f}  bottom "
                  f"{botm[layer, cell]:.3f}  K {k[layer, cell]:.4g}  open {upper - lower:.2f}  "
                  f"saturated {saturated:.2f}  icelltype {celltype[layer, cell]}")
        targets = mf.HeadTargets.from_cells(
            cells=[cell], names=["S"], values=[[np.nan]] * len(period_lengths(model)),
            screen_tops=[screen_top], screen_bottoms=[screen_bottom])
        frame = targets.compare(model)
        print(f"   screened head at the last period: {frame['sim_head'].iloc[-1]:.3f} "
              f"(layers: {', '.join(f'{h:.3f}' for h in heads[:3, cell])})")

    rng = np.random.default_rng(0)
    cells = rng.choice(np.flatnonzero(active[0]), count, replace=False)
    mapping = pd.DataFrame({"obsname": [f"W{i}" for i in range(count)],
                            "name": [f"W{i}" for i in range(count)],
                            "layer": 0, "cell": cells})
    started = time.perf_counter()
    old = _old_compare_route(model)
    old = old.loc[old["layer"] == 0].set_index(["per", "cell"])["sim_head"]
    old_s = time.perf_counter() - started
    started = time.perf_counter()
    new = read_target_values(head_file, mapping, perlen=period_lengths(model),
                             file=head_file.name, saves="HEAD")
    new_s = time.perf_counter() - started
    agree = all(
        np.array_equal(new[f"W{i}"].to_numpy(), old.xs(int(c), level="cell").to_numpy())
        for i, c in enumerate(cells)
    )
    print(f"3. {count} targets x {len(new)} periods: old route {old_s:.2f} s, "
          f"shared reader {new_s:.2f} s; every value identical: {agree}")


if __name__ == "__main__":
    main()
