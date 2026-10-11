"""A cellid read as the cell the model's own grid numbers, on DIS and DISU (ledger 419).

``split_cellid_columns`` -- behind every package INPUT table and the diff's value
tier -- read every FloPy cellid as ``(layer, cell)``, so a DIS ``(layer, row, col)``
named its ROW as its cell. On DISU ``node_to_lni`` followed FloPy's IHC-inferred
layering while the heads table numbers every node in one layer. This script measures
both on real MF6 runs:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/structured_cellids.py

It uses only public spellings present before and after the fix, so it measures
whichever tree is on ``PYTHONPATH`` -- the "before" is ``da33185``'s ``src``
(``git archive da33185 src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``).
Needs ``mf6``.

THE SETUP. The models of ``tests/test_structured_cellids.py``: a 2-layer 3 x 4 DIS
grid of 10 ft cells, CHD 5.0 on the west column and 1.0 on the east one of both
layers; a twin of it; one whose CHD at (1, 2, 3) holds 2.0; and the same layout
written as a DISU (``get_disu_kwargs``) with ``(node,)`` cellids.

MEASURED (mf6 6.8.0.dev0, 2026-10-11), before (da33185) -> after:
  DIS  chd.inputs cells, layer 0       [0, 0, 1, 1, 2, 2] -> [0, 3, 4, 7, 8, 11]
  DIS  self-diff value_cells_changed   12 -> 0
  DIS  one-record change               13 -> 1  (cell 11 of layer 1)
  DISU chd.inputs                      TypeError (int() of a tuple) -> layer 0, nodes
  DISU node 12 in chd.results.q        (1, 0) -> (0, 12); heads table: (0, 12) both

Added in review (the same model with ``hfb=True``: barriers of hydchr 1e-6 across the
column 1 | 2 face of row 1, layer 0, and row 2, layer 1), before -> after:
  DIS  hfb.inputs cells                [(0,1,2), (1,1,2)] -> [(0,5,6), (1,9,10)]
  DIS  hfb.results.q                   -5.557, -14.161 -> -8.54e-05, -2.49e-04
       MF6's own FLOW-JA-FACE across the barriers' faces: -8.54e-05, -2.49e-04
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import myflopy as mf  # noqa: E402
from tests.test_structured_cellids import _dis_run, _disu_run, _face_flow  # noqa: E402


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        base = mf.load_mf6_run(_dis_run(root / "base"))
        twin = mf.load_mf6_run(_dis_run(root / "twin"))
        one = mf.load_mf6_run(_dis_run(root / "one", changed={(1, 2, 3): 2.0}))

        table = base.chd.inputs.get()
        print("DIS  chd.inputs cells, layer 0 ", sorted(table.loc[table.layer == 0, "cell"]))
        for label, other in (("self-diff", twin), ("one-record change", one)):
            summary = mf.ModelGroup({"base": base, "other": other}).diff().summary()
            print(f"DIS  {label:<18} value_cells_changed",
                  int(summary.loc[summary.package == "chd", "value_cells_changed"].iloc[0]))

        disu = mf.load_mf6_run(_disu_run(root / "disu"))
        try:
            inputs = disu.chd.inputs.get()
            print("DISU chd.inputs layers", sorted(set(inputs.layer)), "cells",
                  sorted(inputs.cell))
        except TypeError as error:
            print("DISU chd.inputs TypeError:", error)
        q = disu.chd.results.q.get()
        node_12 = tuple(q.loc[q.node == 12, ["layer", "cell"]].iloc[0])
        print("DISU node 12 in chd.results.q  ", node_12)
        heads = disu.hds.get(per=0)
        print("DISU heads table layers", sorted(set(heads.layer)), "cells 0..", heads.cell.max())

        hfb_ws = _dis_run(root / "hfb", hfb=True)
        barriers = mf.load_mf6_run(hfb_ws)
        table = barriers.hfb.inputs.get()
        print("DIS  hfb.inputs cells ", table[["layer", "cell1", "cell2"]].values.tolist())
        print("DIS  hfb.results.q    ", barriers.hfb.results.q.get()["q_cell1_to_cell2"].tolist())
        ncpl = 12
        print("     MF6 FLOW-JA-FACE ", [
            _face_flow(hfb_ws, layer * ncpl + a, layer * ncpl + b)
            for layer, a, b in ((0, 5, 6), (1, 9, 10))
        ])


if __name__ == "__main__":
    main()
