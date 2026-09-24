"""A head surface with a crater to -1e30 at every dry cell, and nothing above 10,000 ft.

Ledger 315. Found by reading: `InterpolatedSurface.zs` cut MF6's no-value sentinels at
`> 10_000`, where every other reader in `headsplus.py` uses `|value| >= 1e29`.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/head_surface_sentinels.py

THE SETUP. A real triangulated grid -- `GridSpec.voronoi`, `max_area=5000` -- on a
600 x 600 ft square: 207 cells, one layer, `icelltype=1`, NO Newton, steady state, CHD
along the west and east edges. Three MF6 runs on it:

  A. CHD 30 ft west and 50 ft east. The north-east patch (centroid x and y > 400, 25
     cells) has its bottom at 60 ft, above every head, so those cells go DRY.
  B. CHD 9,950 ft west and 10,050 ft east: a model in feet at altitude, whose head
     crosses 10,000 ft mid-domain. Nothing dry.
  C. As A with no raised patch, and an interior block (175 < x, y < 325) at idomain 0.

The OLD rule is applied here by hand to the same raw heads, so this script still shows
the "before" once the library is fixed. Resolution 120.

WHAT IT MEASURES
  1. What MF6 writes: -1e30 for a DRY cell, +1e30 for an INACTIVE one.
  2. The old `> 10_000` cutoff keeps -1e30 as a sample. The cubic interpolators
     spread it across the whole surface, not just the dry patch.
  3. The same cutoff erases every real head above 10,000 ft.
  4. The fixed surface stays inside the real head range under every interpolator.
  5. Two more readers of the same table that are NOT fixed here (ledger 315):
     `hds.map()` and `hds.section()` without `interpolator=` both carry -1e30. A
     section WITH `interpolator=` samples this surface, so it is fixed with it.
"""

from __future__ import annotations

import pathlib
import tempfile
import warnings

import geopandas as gpd
import numpy as np
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

SIDE = 600.0
RESOLUTION = 120
INTERPOLATORS = ["griddata", "rbf", "linearND", "cloughTocher2D"]


def _grid(root: pathlib.Path):
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, SIDE, SIDE)], crs=2927).to_file(boundary)
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=5_000.0
    ).resolve(workspace=root / "grid")


def _run(root, vor, name, *, top, botm, west, east, idomain=None):
    """Build and run one steady non-Newton model; return its model view."""

    props = vor.get_disv_gridprops()
    bounds = vor.gdf_vorPolys.geometry.bounds
    botm = np.broadcast_to(np.asarray(botm, dtype=float), (vor.ncpl,))
    # MF6 refuses a CHD head below its cell's bottom, so a raised cell carries none.
    chd = [[(0, int(c)), west] for c in np.flatnonzero(bounds.minx <= 1e-6) if botm[c] < west]
    chd += [[(0, int(c)), east] for c in np.flatnonzero(bounds.maxx >= SIDE - 1e-6) if botm[c] < east]
    disv = mf.disv(
        nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
        vertices=props["vertices"], cell2d=props["cell2d"], top=top, botm=[botm],
        idomain=None if idomain is None else [idomain],
    )
    model = mf.gwf(name, grid=vor, packages=[
        disv, mf.npf(k=10.0, icelltype=1), mf.ic(strt=(west + east) / 2),
        mf.chd(stress_period_data={0: chd}),
        mf.oc(head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name=name, packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=[name])])
    run = mf.Project(root / name, name=name).run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _raw(model) -> np.ndarray:
    return np.asarray(model.hds.array(masked=False), dtype=float)


def _stats(surface: InterpolatedSurface, lo: float, hi: float) -> str:
    grid = np.asarray(surface.surface, dtype=float)
    finite = grid[np.isfinite(grid)]
    outside = int(((finite < lo - 1.0) | (finite > hi + 1.0)).sum())
    drawn = np.asarray(surface.fig.data[0].z, dtype=float)
    return (f"{finite.size:>6,} finite px, {outside:>6,} outside {lo:,.0f}..{hi:,.0f}; "
            f"drawn z {np.nanmin(drawn):.4g} .. {np.nanmax(drawn):.4g}")


def _compare(model, vor, lo: float, hi: float) -> None:
    raw = _raw(model)
    old = raw.copy()
    old[old > 10_000] = np.nan                                   # the rule before the fix
    for interpolator in INTERPOLATORS:
        before = InterpolatedSurface(zs=old, vor=vor, resolution=RESOLUTION,
                                     interpolator=interpolator)
        after = InterpolatedSurface(model=model, kstpkper=(0, 0), layer=0,
                                    resolution=RESOLUTION, interpolator=interpolator)
        print(f"   {interpolator:>14} before: {_stats(before, lo, hi)}")
        print(f"   {'':>14}  after: {_stats(after, lo, hi)}")


def main() -> None:
    warnings.simplefilter("ignore")
    root = pathlib.Path(tempfile.mkdtemp())
    vor = _grid(root)
    x = np.asarray(vor.centroids_x, dtype=float)
    y = np.asarray(vor.centroids_y, dtype=float)
    patch = (x > 400) & (y > 400)
    block = (x > 175) & (x < 325) & (y > 175) & (y < 325)
    print(f"grid: {vor.ncpl} cells; raised patch {patch.sum()} cells; idomain-0 block {block.sum()}")

    dry = _run(root, vor, "dry", top=100.0, botm=np.where(patch, 60.0, 0.0), west=30.0, east=50.0)
    inactive = _run(root, vor, "inactive", top=100.0, botm=0.0, west=30.0, east=50.0,
                    idomain=np.where(block, 0, 1))
    high = _run(root, vor, "high", top=10_200.0, botm=9_800.0, west=9_950.0, east=10_050.0)

    raw_dry, raw_inactive, raw_high = _raw(dry), _raw(inactive), _raw(high)
    print("\n1. what MF6 wrote")
    print(f"   A: {(raw_dry == -1e30).sum()} cells at -1e30, all of them in the raised patch: "
          f"{bool(np.all(patch[raw_dry == -1e30]))}; {patch.sum()} patch cells")
    print(f"   C: {(raw_inactive == 1e30).sum()} cells at +1e30, the idomain-0 block: "
          f"{bool(np.array_equal(raw_inactive == 1e30, block))}")
    print(f"   B: {(raw_high > 10_000).sum()} of {raw_high.size} real heads above 10,000 ft; "
          f"range {raw_high.min():,.2f} .. {raw_high.max():,.2f}")

    print("\n2 + 4. A, dry cells; real heads 30..50 ft")
    _compare(dry, vor, 30.0, 50.0)
    print("\n3 + 4. B, heads across 10,000 ft")
    _compare(high, vor, 9_950.0, 10_050.0)
    print("\n   C, inactive block (the one case the old rule caught)")
    _compare(inactive, vor, 30.0, 50.0)

    print("\n5. other readers of the A run, not fixed here")
    mapped = np.asarray(dry.hds.map().linear_zs, dtype=float)
    trace = dry.hds.map().fig.data[0]
    print(f"   hds.map(): {(mapped == -1e30).sum()} cells at -1e30; trace zmin/zmax "
          f"{trace.zmin}/{trace.zmax}, so the colour range is {np.nanmin(mapped):.4g} .. "
          f"{np.nanmax(mapped):.4g}")
    mounding = np.asarray(dry.hds.map(show_mounding=True).linear_zs, dtype=float)
    print(f"   hds.map(show_mounding=True): the dry cells read {sorted(set(mounding[patch].tolist()))}")
    line = [(10.0, 450.0), (590.0, 450.0)]
    for kwargs in ({}, {"interpolator": "griddata"}):
        head = np.asarray(dry.hds.section(line=line, **kwargs).data[0].y, dtype=float)
        print(f"   hds.section({', '.join(f'{k}={v!r}' for k, v in kwargs.items())}): "
              f"{(head == -1e30).sum()} of {head.size} head points at -1e30, "
              f"range {np.nanmin(head):.4g} .. {np.nanmax(head):.4g}")


if __name__ == "__main__":
    main()
