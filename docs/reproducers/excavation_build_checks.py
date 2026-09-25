"""Two measurements behind `LayerStack.excavate` (ledger 320, the excavation design's C4).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/excavation_build_checks.py

Needs `mf6` on PATH for part 2 (measured with 6.8.0.dev0). About 5 s.

PART 1 -- `LayerStack.to_disv` became `build().to_disv()`, and that moved nothing.
It used to sample the stack a SECOND time through `LayerSurfaces.to_disv`, passing the
per-layer separations `build()` used so the two would agree; an excavated stack would
then have written the UN-excavated top beside the excavated idomain. The old route is
still callable -- it is the engine method -- so this script builds both for 12 stacks
(three raster stacks: three units with an inactive pinch; a split unit plus a
pass-through pinch; a floor unit declaring `units="feet"`; each under `method="area"`
and `"centroid"`, and with default and non-default `min_sep`/`default_pinch`/
`default_min_thickness`) and compares the package DIGEST and every option.
Measured: 12 of 12 identical (10 distinct digests).

PART 2 -- a steady toy's budget error is the SOLVER'S, not the cut's. A 92-cell
`GridSpec.voronoi` grid on a 600 ft square; a stack at 100 ft of fill (5), silt (1),
sand (30, split in two) and clay (20); a pond of radius 150 ft cut to a bowl bed
(92 ft at the centre, 96.5 at the edge: 21 columns cut, 7 of them through the whole
fill so that a cell is removed, one silt sliver absorbed).
NEWTON, NPF k=10 k33=1 icelltype=1, CHD 98 ft on the east edge in the clay, array
recharge 0.001 ft/d on the top active cell. Run un-excavated and excavated, under IMS
COMPLEX's default closure and under `outer_dvclose=1e-6, inner_dvclose=1e-8,
rcloserecord=[1e-4, "STRICT"]`.
Measured: -0.21% / -0.41% (un-excavated / excavated) with the default closure; 0.00% /
0.00% closed tight. With CHD at 70 ft (every cell above the clay dry) both read +0.59%
at the default closure and 0.00% tight. Every run reports Normal termination. So
`tests/test_excavate.py` closes the solver tight and asserts the discrepancy for BOTH
models -- asserting Normal termination alone proves nothing (the design measured a
steady LAK toy at Normal termination with a 53.89% error).
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import tempfile
import warnings

import geopandas as gpd
import numpy as np
import rasterio
import shapely as shp
from rasterio.transform import from_origin

import myflopy as mf
from myflopy.layers import _reconcile_args
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.specs import ModelContext, package_digest

warnings.simplefilter("ignore")
WORK = pathlib.Path(tempfile.mkdtemp(prefix="excavation_checks_"))


def _raster(path, grid_bounds, z, res):
    x0, y0, x1, y1 = grid_bounds
    nx, ny = int((x1 - x0) / res) + 3, int((y1 - y0) / res) + 3
    xs = x0 + res * (np.arange(nx) + 0.5)
    ys = y1 + res - res * (np.arange(ny) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    with rasterio.open(path, "w", driver="GTiff", height=ny, width=nx, count=1, dtype="float64",
                       crs="EPSG:2927", transform=from_origin(x0, y1 + res, res, res),
                       nodata=-9999.0) as dst:
        dst.write(np.asarray(z(xx, yy), dtype=float), 1)
    return path


def _second_sampling(stack, vor, **kw):
    """The route `LayerStack.to_disv` took before C4: sample the stack again."""

    min_thk, pinch = stack._per_layer_config(kw["default_min_thickness"], kw["default_pinch"])
    _, units = stack._expanded_layers()
    seps, triggers = stack._reconcile_separations(min_thk, pinch, kw["min_sep"], units)
    trigger = stack._trigger_sep(None)
    triggers = [trigger if value is None else value for value in triggers]
    rec_on, which = _reconcile_args("bottom")
    result = stack.build(vor, **kw)
    return stack._layer_surfaces().to_disv(
        vor, idomain=result.idomain, minimum_thickness=min_thk, pinch=pinch,
        length_units=stack.length_units, name="disv", attach=False, reconcile=rec_on,
        which=which, min_sep=seps, trigger_sep=triggers, method=kw["method"],
    )


def part_one():
    vor = irregular_voronoi_grid(nrow=10, ncol=10, cell_size=40.0)
    bounds = vor.gdf_vorPolys.total_bounds
    x0, y0 = bounds[:2]
    ground = _raster(WORK / "ground.tif", bounds,
                     lambda x, y: 120.0 + 0.03 * (x - x0) - 0.02 * (y - y0) + 2.0 * np.sin(x / 37.0), 7.0)
    sand = _raster(WORK / "sand.tif", bounds,
                   lambda x, y: 95.0 + 0.05 * (x - x0) + 3.0 * np.cos(y / 29.0), 7.0)
    clay = _raster(WORK / "clay.tif", bounds, lambda x, y: 90.0 + 0.02 * (y - y0), 7.0)
    stacks = {
        "three": lambda: (mf.LayerStack(top=mf.Raster(ground)).add("sand", bottom=mf.Raster(sand))
                          .add("clay", bottom=mf.Raster(clay), pinch="inactive", min_thickness=2.0)
                          .add("till", thickness=40.0)),
        "split": lambda: (mf.LayerStack(top=mf.Raster(ground))
                          .add("sand", bottom=mf.Raster(sand), split=[0.3, 0.7])
                          .add("clay", bottom=mf.Raster(clay), pinch="passthrough", min_thickness=3.0)
                          .add("till", thickness=30.0, split=2)),
        "floor": lambda: (mf.LayerStack(top=mf.Raster(ground))
                          .add("sand", bottom=mf.Surface.raster(sand, units="feet"), min_thickness=5.0)
                          .add("rock", thickness=25.0)),
    }
    same, digests = 0, set()
    for name, factory in stacks.items():
        for method in ("area", "centroid"):
            for extra in ({}, {"min_sep": 0.5, "default_pinch": "inactive", "default_min_thickness": 2.0}):
                kw = {"default_min_thickness": 1.0, "default_pinch": "floor", "min_sep": 0.1,
                      "method": method, **extra}
                new = factory().to_disv(vor, attach=False, **kw)
                old = _second_sampling(factory(), vor, **kw)
                equal = package_digest(new) == package_digest(old) and all(
                    np.array_equal(np.asarray(new.options[k]), np.asarray(old.options[k]))
                    for k in ("top", "botm", "idomain")
                )
                same += equal
                digests.add(package_digest(new))
                print(f"  {name:6s} {method:9s} {sorted(extra)!s:55s} identical={equal}")
    print(f"PART 1: {same} of 12 identical ({len(digests)} distinct digests)")


def part_two():
    root = WORK / "toy"
    root.mkdir()
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[shp.box(0, 0, 600, 600)], crs=2927).to_file(boundary)
    grid = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary),
                               max_area=12_000.0).resolve(workspace=root / "grid")
    site = root / "pond.gpkg"
    gpd.GeoDataFrame({"name": ["pond"]}, geometry=[shp.Point(300, 300).buffer(150)],
                     crs=2927).to_file(site)
    bowl = _raster(root / "bowl.tif", (0.0, 0.0, 600.0, 600.0),
                   lambda x, y: 92.0 + ((x - 300) ** 2 + (y - 300) ** 2) / 5000.0, 5.0)

    def stack():
        return (mf.LayerStack(top=mf.Flat(100.0)).add("fill", thickness=5.0)
                .add("silt", thickness=1.0).add("sand", thickness=30.0, split=2)
                .add("clay", thickness=20.0))

    plain = stack().build(grid)
    dug = stack().excavate(mf.Raster(bowl), site, name="pond").build(grid)
    cut = int((dug.top < plain.top).sum())
    removed = int(dug.excavated.any(axis=0).sum())
    print(f"PART 2: {cut} columns cut, {removed} with a removed cell, "
          f"{int((dug.absorbed > 0).sum())} sliver(s) absorbed")
    east = [int(c) for c in np.flatnonzero(grid.gdf_vorPolys.geometry.bounds.maxx >= 600 - 1e-6)]
    for chd in (98.0, 70.0):
        for tight in (False, True):
            row = []
            for tag, result in (("un-excavated", plain), ("excavated", dug)):
                context = ModelContext(grid=grid, domain=result.idomain)
                model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
                    mf.disv(result), mf.npf(k=10.0, k33=1.0, icelltype=1), mf.ic(strt=chd + 1.0),
                    mf.chd(stress_period_data={0: [[(4, c), chd] for c in east]}),
                    mf.rch.array(context=context, recharge=1e-3),
                    mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
                ])
                closure = (dict(outer_dvclose=1e-6, inner_dvclose=1e-8, rcloserecord=[1e-4, "STRICT"])
                           if tight else {})
                simulation = mf.simulation(model, name="s", packages=[
                    mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
                    mf.ims(models=["m"], complexity="COMPLEX", **closure)])
                with contextlib.redirect_stdout(io.StringIO()):
                    run = mf.Project(root / f"p_{chd}_{tight}_{tag}", name="p").run("r", simulation)
                listing = (pathlib.Path(run.workspace) / "m.lst").read_text()
                discrepancy = [ln for ln in listing.splitlines() if "PERCENT DISCREPANCY" in ln][-1].split()[-1]
                row.append(f"{tag} {discrepancy}% (success={run.success})")
            print(f"  CHD {chd:5.1f}  {'tight' if tight else 'default':7s} closure: " + ", ".join(row))


if __name__ == "__main__":
    part_one()
    part_two()
