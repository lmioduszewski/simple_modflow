"""What reconcile does to a unit at or above its own minimum, and to the engine's edge cases.

Behind ledger 377. Run it under two source trees and diff::

    git archive <older-commit> src | tar -x -C /tmp/older
    PYTHONPATH=/tmp/older/src python docs/reproducers/reconcile_present_units.py > before.txt
    PYTHONPATH=src            python docs/reproducers/reconcile_present_units.py > after.txt
    diff before.txt after.txt

Measured 2026-10-11 against 8c99c58:

- a 0.8 ft clay, ``min_thickness=0.5``, ``pinch="inactive"``: 0.1 ft and idomain 0
  before, 0.8 ft and idomain 1 after (the same with ``"passthrough"``, -1 before);
  split in two, slices ``[0.1, 0.7]`` before and ``[0.4, 0.4]`` after; with
  ``trigger_sep=2`` 0.1/0 before, 0.8/1 after; a clay exactly AT its 0.5 minimum
  0.1/0 before, 0.5/1 after; a 0.3 ft clay 0.1 ft and idomain 0 on both trees;
- a 2 ft ``"floor"`` unit at its 2 ft minimum split ``[0.9, 0.1]``: slices
  ``[1.8, 1.0]`` and its base 0.8 ft lower before, ``[1.8, 0.2]`` after; a 1 ft one
  held open ``[1.0, 1.0]`` before and ``[1.8, 0.2]`` (by share) after;
- ``reconcile_surfaces``: index ``[10, 11]`` raised ``IndexError`` before; index
  ``[1, 0]`` moved the clear row (label 1, 5.0 -> 9.9) and left the conflicting one;
  ``which="top"`` on ``[10.0, 9.5, 9.45]`` gave ``[9.6, 9.55, 9.45]`` before and
  ``[9.65, 9.55, 9.45]`` after;
- ``adjust_top_btm_overlaps(min_sep=1.0, buffer=0.5)`` set the lowered bottom 0.1
  below the top before (the final reconcile ran on the defaults) and 1.0 after.

Every LayerStack the test suite builds (129 distinct declarations), the two layer
notebooks' and the reproducers' stacks came out byte-identical on both trees.

Ledger 446, measured 2026-10-11 against 3289af0 (an explicit ``trigger_sep`` no
longer replaces a ``"floor"`` layer's trigger; it is capped at the layer's share of the
unit minimum):

- a 3 ft floor unit, ``min_thickness=2``, ``trigger_sep=5``: 2.0 ft and its base 1 ft
  higher before, 3.0 ft after; split 3, slices ``[0.6667] * 3`` before and
  ``[1.0] * 3`` after; a 2.4 ft unit with min 1 split 3 (``test_trigger_sep_can_be_
  overridden``), ``[0.3333] * 3`` and its base 1.4 ft higher before, ``[0.8] * 3``
  after;
- a 1.5 ft floor unit with min 2 is opened to 2.0 by default and under ``trigger_sep``
  1.6 and 5, and left 1.5 under 1.0, on both trees -- what ``trigger_sep`` still does.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import adjust_top_btm_overlaps, reconcile_surfaces

warnings.simplefilter("ignore")
vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
Flat = mf.Flat


def show(label, stack, **build):
    result = stack.build(vor, **build)
    n = len(result.thickness) - 1
    print(f"{label:46s} thickness {np.round(result.thickness[:n, 0], 4).tolist()} "
          f"idomain {result.idomain[:n, 0].tolist()} base {round(float(result.botm[n - 1, 0]), 4)}")


def clay(bottom, **kw):
    return (mf.LayerStack(top=Flat(10.0)).add("clay", bottom=Flat(bottom), **kw)
            .add("below", bottom=Flat(0.0)))


for pinch in ("inactive", "passthrough"):
    show(f"0.8 ft clay, min 0.5, {pinch}", clay(9.2, min_thickness=0.5, pinch=pinch))
show("0.8 ft clay, min 0.5, inactive, split 2", clay(9.2, min_thickness=0.5, pinch="inactive", split=2))
show("0.8 ft clay, min 0.5, inactive, trigger_sep 2", clay(9.2, min_thickness=0.5, pinch="inactive"),
     trigger_sep=2.0)
show("0.5 ft clay, min 0.5, inactive", clay(9.5, min_thickness=0.5, pinch="inactive"))
show("0.3 ft clay, min 0.5, inactive", clay(9.7, min_thickness=0.5, pinch="inactive"))
show("2 ft floor unit, min 2, split [0.9, 0.1]", clay(8.0, min_thickness=2.0, split=[0.9, 0.1]))
show("1 ft floor unit, min 2, split [0.9, 0.1]", clay(9.0, min_thickness=2.0, split=[0.9, 0.1]))
# ledger 446: an explicit trigger_sep on a "floor" unit
show("3 ft floor unit, min 2, trigger_sep 5", clay(7.0, min_thickness=2.0), trigger_sep=5.0)
show("3 ft floor unit, min 2, split 3, trigger_sep 5", clay(7.0, min_thickness=2.0, split=3),
     trigger_sep=5.0)
show("2.4 ft floor unit, min 1, split 3, trigger_sep 5", clay(7.6, min_thickness=1.0, split=3),
     trigger_sep=5.0)
for trigger in (None, 1.0, 1.6, 5.0):
    show(f"1.5 ft floor unit, min 2, trigger_sep {trigger}", clay(8.5, min_thickness=2.0),
         trigger_sep=trigger)

for index in ([10, 11], [1, 0]):
    frame = pd.DataFrame({0: [10.0, 10.0], 1: [5.0, 9.99]}, index=index)
    try:
        print(f"index {index}:", reconcile_surfaces(None, df=frame)[1].to_dict())
    except IndexError as error:
        print(f"index {index}: IndexError {error}")
top = reconcile_surfaces(None, df=pd.DataFrame({0: [10.0], 1: [9.5], 2: [9.45]}), which="top")
print("which='top' [10.0, 9.5, 9.45]:", np.round(top.iloc[0].to_numpy(), 4).tolist())


class _Grid:
    gdf_topbtm = pd.DataFrame({0: [10.0, 10.0, 10.0], 1: [5.0, 5.0, 5.0]})

    def reconcile_surfaces(self, df=None, min_sep=0.1, trigger_sep=1, which="bottom"):
        return reconcile_surfaces(self, df=df, min_sep=min_sep, trigger_sep=trigger_sep, which=which)

    def get_vor_cells_as_series(self, shp):
        return pd.Series([[0]])

    def find_adjacent_cells(self, cell):
        return [1, 2]


grid = _Grid()
out = adjust_top_btm_overlaps(grid, min_sep=1.0, buffer=0.5)
print("adjust_top_btm_overlaps(min_sep=1.0, buffer=0.5): lowered bottom", round(float(out.loc[0, 1]), 4))
