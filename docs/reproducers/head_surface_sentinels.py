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
  D. A's raised patch AND C's idomain-0 block, for the readers in section 5.

The OLD rule is applied here by hand to the same raw heads, so this script still shows
the "before" once the library is fixed. Resolution 120.

WHAT IT MEASURES
  1. What MF6 writes: -1e30 for a DRY cell, +1e30 for an INACTIVE one.
  2. The old `> 10_000` cutoff keeps -1e30 as a sample. The cubic interpolators
     spread it across the whole surface, not just the dry patch.
  3. The same cutoff erases every real head above 10,000 ft.
  4. The fixed surface stays inside the real head range under every interpolator.
  5. Every other reader of the same table, on D, fixed in the restatement of ledger
     315. Each "before" applies that reader's own old rule by hand: the map and its
     hover blanked only `== 1e30`, the mounding hover and the depth-below-ground map
     blanked nothing, and neither did `hds.section()` without `interpolator=`. The
     contours are the old values handed to `contour_line_segments`. The user's calls
     are in the "after": on the mounding map a DRY cell reads 0, and on the
     depth-below-ground map it is blank.
"""

from __future__ import annotations

import atexit
import pathlib
import shutil
import tempfile
import warnings

import geopandas as gpd
import numpy as np
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.contour_plotting import contour_line_segments
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.utils.datatypes.xsections import XSection

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
    atexit.register(shutil.rmtree, root, ignore_errors=True)  # /tmp is a quota'd tmpfs (ledger 460)
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

    both = _run(root, vor, "both", top=100.0, botm=np.where(patch, 60.0, 0.0), west=30.0,
                east=50.0, idomain=np.where(block, 0, 1))
    print("\n5. every other reader, on D; before = that reader's old rule, by hand")
    _readers(both, vor, patch, block)


def _tally(values, dry, inactive) -> str:
    """How ``values`` reads at the dry and the inactive cells, and its finite range."""

    def at(cells):
        seen = sorted({"blank" if np.isnan(v) else f"{v:.4g}" for v in values[cells]})
        return "/".join(seen)

    return (f"dry {at(dry)}, inactive {at(inactive)}; "
            f"range {np.nanmin(values):.4g} .. {np.nanmax(values):.4g}")


def _readers(model, vor, patch, block) -> None:
    raw = _raw(model)
    dry, inactive = raw == -1e30, raw == 1e30
    top = np.asarray(model.gwf.modelgrid.top, dtype=float)
    bottom = np.asarray(model.gwf.modelgrid.botm[0], dtype=float)
    print(f"   D: {dry.sum()} dry cells (the patch: {bool(np.array_equal(dry, patch))}), "
          f"{inactive.sum()} inactive (the block: {bool(np.array_equal(inactive, block))})")

    old_map = np.where(raw == 1e30, np.nan, raw)
    plain = model.hds.map()
    print(f"   map + sectioned hover  before: {_tally(old_map, dry, inactive)}")
    print(f"   {'':>22} after: {_tally(np.asarray(plain.linear_zs, float), dry, inactive)}")
    hover = np.asarray(plain._build_hover_context().layer_fields["head"][0], dtype=float)
    print(f"   {'':>16} after, hover: {_tally(hover, dry, inactive)}")

    before = len(contour_line_segments(vor, old_map, levels=10))
    after = len(contour_line_segments(vor, np.asarray(plain.linear_zs, float), levels=10))
    drawn = sum(t.type == "scattermap" for t in model.hds.map(contours=True).fig.data)
    print(f"   contours               before: {before} segments; after: {after} segments, "
          f"drawn as {drawn} traces")

    old_mound_map = np.clip(old_map - bottom, 0.0, None)
    old_mound_hover = np.clip(raw - bottom, 0.0, None)
    mound = model.hds.map(show_mounding=True)
    print(f"   mounding map           before: {_tally(old_mound_map, dry, inactive)}")
    print(f"   mounding hover         before: {_tally(old_mound_hover, dry, inactive)}")
    print(f"   {'':>16} after, both: {_tally(np.asarray(mound.linear_zs, float), dry, inactive)}"
          f" | {_tally(np.asarray(mound._mounding_series(), float), dry, inactive)}")

    print(f"   depth below ground     before: {_tally(raw - top, dry, inactive)}")
    bgs = np.asarray(model.hds.map(bgs=True).linear_zs, dtype=float)
    print(f"   {'':>22} after: {_tally(bgs, dry, inactive)}")

    line = [(10.0, 10.0), (590.0, 590.0)]
    crossed = XSection(model=model, line=line, kstpkper=(0, 0)).xs.index.to_numpy()
    drawn = np.asarray(model.hds.section(line=line).data[0].y, dtype=float)
    print(f"   hds.section()          before: {_tally(raw[crossed], dry[crossed], inactive[crossed])}")
    print(f"   {'':>22} after: {_tally(drawn, dry[crossed], inactive[crossed])}")

if __name__ == "__main__":
    main()
