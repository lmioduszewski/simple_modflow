"""The measurement behind ledger 430: what the every-period boundary export changes.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/usg_boundary_periods.py [OLD_EXPORT]

PART 1 -- the synthetic DRN of ``tests/test_usg_export_periods.py``: drains on cells 0
and 1 in period 0, cell 2 ADDED in period 1, switched off (``ITMP = 0``) in period 2
and reused off in period 3. Before: cells [0, 1], no trace of the switch-off. Now:
cells [0, 1, 2], ``periods`` "0-1"/"0-1"/"1", and periods [2, 3] switched off.

PART 2 -- the Ten Trails MODFLOW-USG model (``flow-tt01_USE_5yr.nam``, 72 periods),
when it is on this machine. Its BC cell sets never change between periods (CHD
defined in all 72, DRN and GHB in period 0 and reused), so the export must not move:
measured DRN 205 / GHB 88 / CHD 69 rows, every value, geometry and cycle column
identical to the export before the change, with ``first_period``/``periods`` added
and no period switched off. Pass the old ``export.py`` (``git show
<before>:src/myflopy/modflow/usg/export.py > old_export.py``) as OLD_EXPORT to compare
column by column; without it the new frames alone are summarised. Reading the model
takes about a minute.
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys
import warnings

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from myflopy.modflow.usg import read_usg  # noqa: E402
from myflopy.modflow.usg.export import _bc_frames, _switched_off  # noqa: E402
from tests.test_usg_export_periods import _drain_model  # noqa: E402

TEN_TRAILS = pathlib.Path(
    "/home/lukem/models/mf6/Ten Trails/SSPA/Model_FINAL_ForLuke/14-Calib40_Iter01_DV"
)


def part_1() -> None:
    model = _drain_model()
    points, _ = _bc_frames(model, "DRN", "EPSG:2927")
    print("PART 1  cells", sorted(points["old_cell"]),
          "periods", points.set_index("old_cell")["periods"].to_dict(),
          "switched off", _switched_off(model, "DRN"))


def part_2(old_export: str | None) -> None:
    nam = TEN_TRAILS / "flow-tt01_USE_5yr.nam"
    if not nam.exists():
        print("PART 2  skipped: the Ten Trails model is not on this machine")
        return
    warnings.filterwarnings("ignore")
    usg = read_usg(nam, gsf=TEN_TRAILS / "flow-tt01.gsf", crs="EPSG:2927")
    old = None
    if old_export:
        spec = importlib.util.spec_from_file_location("old_export", old_export)
        old = importlib.util.module_from_spec(spec)
        sys.modules["old_export"] = old
        spec.loader.exec_module(old)
    for ftype in ("DRN", "GHB", "CHD"):
        points, lines = _bc_frames(usg, ftype, "EPSG:2927")
        line = (f"PART 2  {ftype} rows {len(points)} periods {sorted(set(points['periods']))} "
                f"switched off {_switched_off(usg, ftype)}")
        if old is not None:
            before, before_lines = old._bc_frames(usg, ftype, "EPSG:2927")
            same = all(
                b.drop(columns="geometry").reset_index(drop=True).equals(
                    a[b.columns].drop(columns="geometry").reset_index(drop=True))
                and b.geometry.reset_index(drop=True).geom_equals(
                    a.geometry.reset_index(drop=True)).all()
                for b, a in ((before, points), (before_lines, lines))
            )
            added = [c for c in points.columns if c not in before.columns]
            line += f" | identical to before: {same}, added {added}"
        print(line)


if __name__ == "__main__":
    part_1()
    part_2(sys.argv[1] if len(sys.argv) > 1 else None)
