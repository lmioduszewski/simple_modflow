"""The measurements behind the USG importer fixes of ledgers 367, 368 and 369.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/usg_switch_off_crs_pixels.py

Needs `mf6` on PATH (measured with 6.8.0.dev0). A few seconds. The model is the
synthetic 2-layer, 4-cell MODFLOW-USG model of ``tests/test_usg_import.py`` (10 ft
cells, EPSG:2927 in US survey feet), with its drains switched off in period 2
(``ITMP = 0``).

PART 1 (367) -- MODFLOW 6 on the converted model, period 2 written two ways: LEFT
OUT (what the conversion did) and EMPTY (what it does now). MF6 reads a missing period
as "reuse the previous one", so the left-out form keeps the drains running.
Measured: left out -64.41 ft3/d of DRN in period 2; empty 0.0.

PART 2 (368) -- the CRS labels on one grid after ``to_mf6(crs="EPSG:2926")``.
Before: the grid read 2926 while ``gdf_vorPolys``/``gdf_topbtm`` and the lat/lon frame
stayed on 2927. Now every one reads 2926 and the lat/lon frame is rebuilt.

PART 3 (369) -- the pixel size of the exported rasters when the same feet model
is written to EPSG:2856 (the same projection in metres). Before: 2.5 m pixels for a
3.05 m cell, the default measured in feet and applied in metres. Now 0.762 m.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile

import flopy
import rasterio

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

from myflopy.modflow.usg import read_usg  # noqa: E402
from tests.test_usg_export import workspace as export_workspace  # noqa: E402
from tests.test_usg_import import _SWITCHED_OFF, _write_usg_model  # noqa: E402


def drain_flows(workspace: pathlib.Path) -> dict[int, float]:
    budget = flopy.utils.CellBudgetFile(str(workspace / "tiny.cbc"))
    return {
        int(kper): float(budget.get_data(text="DRN", kstpkper=(kstp, kper))[0]["q"].sum())
        for kstp, kper in budget.get_kstpkper()
    }


def part1(root: pathlib.Path) -> None:
    (root / "usg").mkdir()
    nam = _write_usg_model(root / "usg", packages=_SWITCHED_OFF)
    model = read_usg(nam, gsf=root / "usg" / "tiny.gsf", crs="EPSG:2927")
    for label, keep_empty in (("left out", False), ("empty", True)):
        sim = model.to_mf6("tiny", start_date_time="2020-01-01", include=("chd", "drn"))
        drn = next(p for p in sim.models[0].packages if p.name == "drn")
        data = drn.options["stress_period_data"]
        if not keep_empty:
            data.pop(1, None)
        workspace = root / label.replace(" ", "_")
        ok, _ = sim.build(workspace).execute()
        print(f"PART 1  period 2 {label:>8}: ran={ok}  DRN by period {drain_flows(workspace)}")


def part2(root: pathlib.Path) -> None:
    (root / "crs").mkdir()
    nam = _write_usg_model(root / "crs")
    model = read_usg(nam, gsf=root / "crs" / "tiny.gsf", crs="EPSG:2927")
    grid = model.grid
    _ = grid.gdf_topbtm, grid.gdf_latlon
    model.to_mf6("tiny", start_date_time="2020-01-01", crs="EPSG:2926")
    print(
        f"PART 2  grid {grid.crs}  gdf_vorPolys {grid.gdf_vorPolys.crs}  "
        f"gdf_topbtm {grid.gdf_topbtm.crs}  lat/lon bounds {grid.gdf_latlon.total_bounds.round(4)}"
    )


def part3(root: pathlib.Path) -> None:
    class _Factory:
        def mktemp(self, name):
            path = root / name
            path.mkdir()
            return path

    ws = export_workspace.__wrapped__(_Factory())
    model = read_usg(ws / "tiny.nam", gsf=ws / "tiny.gsf", crs="EPSG:2927")
    out = root / "metres"
    model.export_gis(out, crs="EPSG:2856")
    for name in ("top.tif", "lakebed_bottom.tif"):
        with rasterio.open(out / "arrays" / name) as handle:
            print(f"PART 3  {name:<20} pixel {handle.res[0]:.4f} m  ({handle.width} x {handle.height})")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        part1(root)
        part2(root)
        part3(root)
