"""The measurements behind `mf.Basin`, an assumed lakebed (ledger 554).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/basin_volume.py

Needs `mf6` on PATH for PART 3 (measured with 6.8.0.dev0+0418efe4).

PART 1 -- volume against the profile. A circular lake of radius 300 ft and maximum depth
20 ft (a 512-segment polygon) on `irregular_voronoi_grid` meshes of 200, 100, 50 and
25 ft cells over 1000 x 1000 ft, water surface flat at 100 ft. The exact volume is
pi R^2 D / 2 for the bowl (a paraboloid of revolution) and pi R^2 D / 3 for the cone.
"area" is what `mf.Basin` writes -- each covered cell the AREA MEAN of the depth over
the whole cell; "centroid" is the profile read at each covered cell's centroid instead.
Measured (percent volume error, area / centroid):
  bowl  200 ft: -0.00 / +9.99    100 ft: -0.06 / +2.06    50 ft: +0.02 / +0.03    25 ft: -0.00 / +0.13
  cone  200 ft: -0.00 / +6.56    100 ft: -0.06 / +1.67    50 ft: +0.02 / -0.07    25 ft: +0.00 / +0.06
(the item's design measured -12.0% / -10.3% for the centroid at 200 ft on square cells:
its sign and size depend on where the centroids fall, the area mean's do not). The
deepest cell's mean depth stays under 20 (16.9, 18.3, 19.8, 20.0 for the bowl): an area
mean of a profile that peaks at one point never reaches the peak.

PART 2 -- bowl against cone at one maximum depth: 1.5000 of the cone's volume on the
200, 100 and 50 ft meshes, 1.4999 on the 25 ft one.

PART 3 -- the pairing, run in MF6. 600 x 600 ft of 20 ft Voronoi cells, a stack of
100 / 95 / 80 / 40 ft, a 120 ft-radius pond declared in a GeoPackage with a
`max_depth` column of 12 ft; `bed = mf.Basin(pond, "max_depth", water_surface=100.0)`,
`stack.excavate(bottom=bed, within=pond)`, `mf.lak(lakes=pond)` taking its bed from
the excavation record. CHD 97 ft west and 93 ft east (layers 2-3), k = k33 = 10,
icelltype 1, NEWTON, steady, bed leakance 1/d. Measured: the record's bed equals
`sample_bed(bed, ...)` at every cell (asserted below); 139 lake cells, the deepest
cell's bed 88.24 ft; MF6 normal termination, GWF and LAK percent discrepancy both 0.00,
lake stage 95.05 ft (about halfway between the boundary heads, as a lake in uniform flow
should sit). About 10 s.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import tempfile
import warnings

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import basin_depths
from myflopy.surfaces import sample_bed

RADIUS, DEPTH = 300.0, 20.0


def part_1_and_2():
    print("PART 1/2 -- volume error (%), area mean / centroid; bowl/cone ratio")
    for cell in (200.0, 100.0, 50.0, 25.0):
        n = int(1000 / cell)
        grid = irregular_voronoi_grid(nrow=n, ncol=n, cell_size=cell)
        polys = grid.gdf_vorPolys
        x0, y0, x1, y1 = polys.total_bounds
        disk = Point((x0 + x1) / 2, (y0 + y1) / 2).buffer(RADIUS, 512)
        area = polys.area.to_numpy()
        volumes = {}
        for shape, mean in (("bowl", 1 / 2), ("cone", 1 / 3)):
            bed = mf.Basin(mf.Area(disk), DEPTH, water_surface=100.0, shape=shape).values(grid)
            cells = np.flatnonzero(np.isfinite(bed))
            centres = np.array([g.centroid.coords[0] for g in polys.geometry.loc[cells]])
            centroid = basin_depths(centres[:, 0], centres[:, 1], [(disk, DEPTH)], shape=shape)
            exact = np.pi * RADIUS**2 * DEPTH * mean
            volumes[shape] = float((area[cells] * (100.0 - bed[cells])).sum())
            print(f"  {shape} {cell:5.0f} ft: {100 * (volumes[shape] / exact - 1):+.2f} / "
                  f"{100 * (float((area[cells] * centroid).sum()) / exact - 1):+.2f}   "
                  f"deepest cell {np.nanmax(100.0 - bed):.1f} ft")
        print(f"  bowl / cone at {cell:.0f} ft: {volumes['bowl'] / volumes['cone']:.4f}")


def part_3(root):
    print("PART 3 -- excavate + mf.lak on a basin bed, in MF6")
    grid = irregular_voronoi_grid(nrow=30, ncol=30, cell_size=20.0, refined=False)
    x0, y0, x1, y1 = grid.gdf_vorPolys.total_bounds
    pond = root / "pond.gpkg"
    gpd.GeoDataFrame({"name": ["pond"], "max_depth": [12.0]},
                     geometry=[Point((x0 + x1) / 2, (y0 + y1) / 2).buffer(120.0, 128)],
                     crs=grid.crs).to_file(pond)
    bed = mf.Basin(pond, "max_depth", water_surface=100.0)
    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("fill", bottom=mf.Flat(95.0))
             .add("upper", bottom=mf.Flat(80.0)).add("lower", bottom=mf.Flat(40.0)))
    result = stack.excavate(bottom=bed, within=pond, name="pond").build(grid, attach=True)
    record = result.excavations["pond"]
    sampled = sample_bed(bed, grid, np.full(grid.ncpl, 100.0))[list(record.cells)]
    assert np.array_equal(record.bed, sampled), "the record IS the basin's sampled bed"
    print(f"  {len(record.cells)} lake cells, deepest bed {record.bed.min():.2f} ft")

    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    west, east = np.flatnonzero(xy[:, 0] < x0 + 20), np.flatnonzero(xy[:, 0] > x1 - 20)
    chd = [[(k, int(c)), h] for cells, h in ((west, 97.0), (east, 93.0)) for c in cells
           for k in (1, 2) if result.idomain[k, c] > 0]
    lake = mf.lak(context=context, nper=1, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=1.0, unexcavated="raise")
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: chd}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-8, inner_dvclose=1e-10,
               outer_maximum=300, inner_maximum=300, linear_acceleration="BICGSTAB")])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    workspace = pathlib.Path(run.workspace)
    listing = (workspace / "m.lst").read_text()
    discrepancy = [float(line.split()[-1]) for line in listing.splitlines()
                   if "PERCENT DISCREPANCY" in line]
    stage = flopy.utils.HeadFile(str(workspace / "m_lak_stage.lak"), text="STAGE").get_data().ravel()
    print(f"  success={run.success}  percent discrepancy (GWF, LAK) = {discrepancy}  "
          f"stage = {float(stage[0]):.2f} ft")


if __name__ == "__main__":
    part_1_and_2()
    with tempfile.TemporaryDirectory() as tmp:
        part_3(pathlib.Path(tmp))
