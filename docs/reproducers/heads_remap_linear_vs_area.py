"""Carrying a SOLVED head field to another Voronoi mesh: linear between centres, or area mean?

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python \
        docs/reproducers/heads_remap_linear_vs_area.py [workspace]

Needs ``mf6`` on PATH (measured with 6.8.0.dev0). Runs MF6 sixteen times, ONE model at a
time, on models of 401 and 2,136 cells: about 11 s in all. Writes nothing outside the
workspace it is given (a temporary directory, removed afterwards, when none is).

THE QUESTION. ``model.hds.remap(new_grid, method=...)`` carries a solved head field onto
another mesh, e.g. to give a model on a new grid its starting heads,
``mf.ic(strt=old.hds.remap(new_vor))``. Two methods are on the table:

  linear  a head from MF6 is a POINT value at its cell's DISV centre (a finite-volume
          node), so interpolate linearly between the OLD cells' DISV centres (on the
          Delaunay triangulation of the old generators) and read it at the NEW cells' DISV
          centres. (The shipped engine adds what this bench cannot exercise: a new centre
          outside the old centres' hull, or in a triangle bridging a gap, takes the
          nearest old cell's value moved along its own gradient, clamped to the layer's
          range -- ``field_remap_rim.py`` measures that part.)
  area    treat each old head as constant over its cell polygon and give each new cell
          the area-weighted MEAN of the old cells under it --
          ``grid.geometry.cell_overlay(old, new).reduce(values, rule="mean")``, the engine
          ``with_grid(remap=True)`` already uses for ``ic.strt``; a new cell no old cell
          reaches takes the nearest old centre's value.

Both run through the shipped ``myflopy.modflow.mf6.grid.remap.remap_field``. The figures
below were first measured with standalone prototypes of the two methods and re-measured
through the engine: identical, because on these meshes no cell needs the rim fallback.

THE BENCH. One steady GWF problem solved on TWO ``GridSpec.voronoi`` meshes (Triangle +
FloPy ``tri2vor``, no point or line constraints) of the SAME 3,000 ft x 3,000 ft square,
EPSG:2927, uniform refinement: ``max_area=40_000`` (coarse, 401 cells, median width 139 ft)
and ``max_area=4_500`` (fine, 2,136 cells, median width 70 ft). Two layers: top 100 ft,
bottoms 50 and 0 ft; K = 20 and 50 ft/d, K33 = 2 and 5 ft/d; ICELLTYPE 1 (convertible) on
layer 1, 0 on layer 2 -- so the solve is non-linear (Picard) and a starting head changes the
iteration count. CHD on every cell touching the west edge (95 ft, both layers) and the east
edge (85 ft, both layers): a regional gradient of 10 ft over 3,000 ft. RCH 0.002 ft/d on
layer 1 (a ~0.7 ft mound). One WEL in layer 2, -25,000 ft3/d, in whichever cell holds the
well point on each mesh -- the local feature, a cone 5-7 ft deep at the well cell (layer 2). IMS
SIMPLE with outer/inner DVCLOSE 1e-6/1e-7. No STO: steady state. No cell goes dry.

Two placements of the well, because MF6 puts a well's sink at its cell's NODE wherever in
the cell the well is:

  off-node  the well at (1800, 1700): 113 ft from the coarse well cell's node and 24 ft
            from the fine one's, so the two models' cones are centred ~100 ft apart -- a
            DISCRETIZATION difference that every remap inherits.
  on-node   the well at the coarse node, within 300 ft of (1800, 1700), that lies closest
            to a fine node: (1694.4, 1967.3), exactly on the coarse well cell's node and
            17.8 ft from the fine one's, so both cones are centred within a quarter of a
            fine cell of each other.

WHAT IT MEASURES, per placement.

  (a) coarse -> fine: the coarse solution remapped onto the fine mesh by each method,
      against the FINE model's own solution, per layer: RMSE, max |error|, mean bias
      (remapped minus model). Split into cells near the well (within 2.5 coarse cell
      widths, 346 ft) and the rest, and -- away from the well -- by where a fine centre
      sits relative to the coarse nodes, ``t = d1 / d2`` (distance to the nearest coarse
      centre over distance to the second nearest: ~0 at a coarse node, ~1 on a coarse face).
  (b) fine -> coarse: the same against the coarse model's own solution.
  (c) the fine model restarted from each remapped field, from a flat 95 ft start and from
      its own solution; the coarse one likewise from the fine field. Outer iterations
      ("CALLS TO NUMERICAL SOLUTION") and total inner iterations from ``mfsim.lst``.
  (d) the same two remaps of an ANALYTIC point field -- gradient + recharge mound, and that
      plus a Thiem cone clipped at r = 10 ft -- evaluated exactly at both meshes' centres.
      This isolates the INTERPOLATION error: (a) and (b) also contain the difference between
      the two discretizations, which no remap can remove (the coarse well cell cannot hold
      the fine model's drawdown).

WHAT IT DOES NOT SHOW. A transient start (only a steady solve, and a mildly non-linear one,
so (c) counts the iterations of one steady solve, not a transient drift from a bad IC). Dry
or inactive cells (none here; the dry -> water-table rule is untested by this bench). A
layer map with a different layer count (both models have the same two layers). Meshes of
different footprints, a non-convex domain or a hole: the hull of a ``tri2vor`` mesh's
generators IS its boundary (boundary generators sit on it), so here the out-of-hull nearest
fill never fires, and nothing tests linear interpolation ACROSS a concave notch, a hole or
an inactive block, which a Delaunay triangulation of the old centres would span. Locally
refined meshes (a ``points=`` constraint at the well made Triangle grade the mesh down to
cells under 0.5 ft2 at the point -- 898 cells at ``max_area=40_000`` against 401 without it
-- so the well is placed on an existing node instead). Coarsening ACROSS a sharp feature -- a K
contrast, an HFB, a drain or river line -- where the head has a kink or a step rather than
smooth curvature.

MEASURED (2026-10-15; layer 2 unless said; RMSE / max |error| in ft):

  on-node, away from the well (> 346 ft), against the other model's own solution:
      coarse -> fine   linear 0.0170 / 0.137    area 0.1094 / 0.510   (6.4x RMSE)
      fine -> coarse   linear 0.0160 / 0.114    area 0.0419 / 0.165   (2.6x)
  on-node, every cell:
      coarse -> fine   linear 0.0532 / 2.068    area 0.1277 / 1.865
      fine -> coarse   linear 0.0617 / 1.176    area 0.0576 / 0.764
  off-node, every cell:
      coarse -> fine   linear 0.1345 / 3.348    area 0.1731 / 1.855
      fine -> coarse   linear 0.0994 / 1.082    area 0.1059 / 1.096
  analytic smooth field (pure interpolation error), every cell:
      coarse -> fine   linear 0.00181 / 0.0069  area 0.1040 / 0.366   (57x)
      fine -> coarse   linear 0.00026 / 0.0008  area 0.0388 / 0.137   (150x)

  * Away from the well the AREA error is its own interpolation error (on-node coarse ->
    fine: 0.109 against the model, 0.109 against the analytic field), while LINEAR's is the
    two discretizations' own difference (0.017 against the model, 0.0075 analytic). The
    area mean is O(gradient x old cell width): a fine cell inside a coarse cell gets the
    coarse NODE's value. Fine -> coarse it returns the head at the new cell's polygon
    CENTROID, while MF6 puts the head at the DISV node -- 11.6 ft from the centroid at the
    median and 47 ft at worst on the coarse mesh; the smooth field differs between the two
    points by RMSE 0.0404 / max 0.1396, which IS area's 0.0388 / 0.137 above.
  * Off-node, both methods' "away" error is dominated by the two models' cones being
    centred ~100 ft apart (linear 0.080, area 0.134) -- no remap can remove that.
  * At the well cell itself both are wrong by 0.8-3.4 ft and NEITHER wins consistently: the
    coarse field does not hold the fine model's drawdown, and the fine field's drawdown is
    deeper than a coarse well cell can be. Linear's max is the larger in 3 of 4 cases.
  * Error by position between old nodes (t), away from the well: linear's is flat
    (0.011-0.012 on-node coarse -> fine at every t); area's peaks between node and face
    (0.056 at the node, 0.100 midway, 0.070 at the face, where a straddling new cell
    averages two old values).
  * Iterations: either remap cuts the inner iterations by 30-43% against a flat 95 ft
    start (fine 204-209 -> 126-142; coarse 95-97 -> 55-65) and outer from 5 to 4; the two
    methods are within 0-7 inner iterations of each other. The problem is only mildly
    non-linear, so this does not separate them.
  * Cost, through the shipped engine (its overlay and frame guard included): 0.07-0.09 s
    per two-layer field on an idle machine, for either method -- the overlay dominates.
    The standalone prototypes' interpolation alone took 0.006-0.054 s.
"""

from __future__ import annotations

import re
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np

SIDE = 3000.0
CRS = 2927
TOP, BOTM = 100.0, [50.0, 0.0]
K, K33, ICELLTYPE = [20.0, 50.0], [2.0, 5.0], [1, 0]
CHD_WEST, CHD_EAST = 95.0, 85.0
RECHARGE = 0.002
WELL_TARGET, WELL_LAYER, WELL_Q = (1800.0, 1700.0), 1, -25_000.0
STRT_FLAT = 95.0
MAX_AREA = {"coarse": 40_000.0, "fine": 4_500.0}
NEAR_WELL_WIDTHS = 2.5
METHODS = ("linear", "area")


# ----------------------------------------------------------------------------- the model


def make_grid(parent: Path, tag: str):
    import geopandas as gpd
    from shapely.geometry import box

    import myflopy as mf

    parent.mkdir(parents=True, exist_ok=True)
    boundary = parent / "boundary.gpkg"
    if not boundary.exists():
        gpd.GeoDataFrame(
            {"n": ["domain"]}, geometry=[box(0.0, 0.0, SIDE, SIDE)], crs=CRS
        ).to_file(boundary)
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=MAX_AREA[tag]
    ).resolve(workspace=parent / f"grid_{tag}")


def centres(vor) -> np.ndarray:
    """The DISV ``cell2d`` centres -- where MF6 places each cell's head."""

    return np.array([[row[1], row[2]] for row in vor.get_disv_gridprops()["cell2d"]], float)


def width(vor) -> float:
    return float(np.sqrt(np.median(vor.gdf_vorPolys.geometry.area)))


def well_cell(vor, xy) -> int:
    from shapely.geometry import Point

    hits = np.flatnonzero(vor.gdf_vorPolys.geometry.intersects(Point(*xy)).to_numpy())
    # A point on a face touches two polygons; take the one whose node is nearest.
    node = centres(vor)[hits]
    return int(hits[np.argmin(np.hypot(node[:, 0] - xy[0], node[:, 1] - xy[1]))])


def on_node_well(coarse, fine):
    """The coarse node within 300 ft of the target that lies closest to a fine node."""

    from scipy.spatial import cKDTree

    c, f = centres(coarse), centres(fine)
    near = np.flatnonzero(np.hypot(c[:, 0] - WELL_TARGET[0], c[:, 1] - WELL_TARGET[1]) < 300)
    gap, _ = cKDTree(f).query(c[near])
    best = near[np.argmin(gap)]
    return (float(c[best, 0]), float(c[best, 1]))


def chd_cells(vor):
    bounds = vor.gdf_vorPolys.geometry.bounds
    west = np.flatnonzero(bounds.minx.to_numpy() <= 1e-6)
    east = np.flatnonzero(bounds.maxx.to_numpy() >= SIDE - 1e-6)
    return west, east


def solve(vor, workspace: Path, *, strt, well_xy) -> dict:
    """Build, write and run one steady model; heads ``(nlay, ncpl)`` and iteration counts."""

    import flopy

    import myflopy as mf

    props = vor.get_disv_gridprops()
    west, east = chd_cells(vor)
    chd = [[(k, int(c)), CHD_WEST] for k in range(2) for c in west]
    chd += [[(k, int(c)), CHD_EAST] for k in range(2) for c in east]
    model = mf.gwf("m", grid=vor, packages=[
        mf.disv(nlay=2, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=TOP, botm=BOTM),
        mf.npf(k=K, k33=K33, icelltype=ICELLTYPE),
        mf.ic(strt=strt),
        mf.chd(stress_period_data={0: chd}),
        mf.rch.array(recharge=RECHARGE, irch=0),
        mf.wel(stress_period_data={0: [[(WELL_LAYER, well_cell(vor, well_xy)), WELL_Q]]},
               boundnames=False),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="SIMPLE", print_option="SUMMARY",
               outer_dvclose=1e-6, inner_dvclose=1e-7, outer_maximum=200, inner_maximum=300),
    ])
    if workspace.exists():
        shutil.rmtree(workspace)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        built = simulation.build(workspace)
        built.simulation.write_simulation(silent=True)
        ok, _ = built.simulation.run_simulation(silent=True)
    listing = (workspace / "mfsim.lst").read_text()
    assert ok and "Normal termination" in listing, listing[-2000:]
    outer = int(re.search(r"(\d+)\s+CALLS TO NUMERICAL SOLUTION", listing).group(1))
    inner = int(re.search(r"(\d+)\s+TOTAL ITERATIONS", listing).group(1))
    heads = flopy.utils.HeadFile(workspace / "m.hds").get_data(kstpkper=(0, 0))
    heads = np.asarray(heads, float).reshape(2, -1)
    assert np.all(np.abs(heads) < 1e29), "this bench has no dry or inactive cell"
    shutil.rmtree(workspace)          # one model's files at a time
    return {"heads": heads, "outer": outer, "inner": inner}


# ----------------------------------------------------------------------------- the remaps


def both_remaps(old_vor, new_vor, field, overlay):
    """``{"linear": (nlay, ncpl), "area": ...}`` plus fill counts and timings.

    Both through the SHIPPED engine, ``myflopy.modflow.mf6.grid.remap.remap_field`` -- the
    one ``model.hds.remap`` calls -- so what this measures is what the library does.
    ``overlay`` is unused (the engine computes its own); it is kept so the callers below
    read as they did when this compared two prototypes.
    """

    from myflopy.modflow.mf6.grid.remap import remap_field

    del overlay
    out, info = {}, {}
    for method in METHODS:
        started = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            out[method], report = remap_field(field, old_vor, new_vor, method=method)
        info[method] = {
            "filled": report["outside"] + sum(report["no_value_under"].values()),
            "seconds": time.perf_counter() - started,
        }
    return out, info


# ----------------------------------------------------------------------------- the measures


def stats(estimate, truth, mask=None):
    error = estimate - truth
    if mask is not None:
        error = error[mask]
    return {
        "n": int(error.size),
        "rmse": float(np.sqrt(np.mean(error ** 2))),
        "max": float(np.max(np.abs(error))),
        "bias": float(np.mean(error)),
    }


def face_position(old_xy, new_xy):
    """``d1 / d2`` per new centre: ~0 at an old node, ~1 on an old cell face."""

    from scipy.spatial import cKDTree

    distance, _ = cKDTree(old_xy).query(new_xy, k=2)
    return distance[:, 0] / np.maximum(distance[:, 1], 1e-12)


def analytic(xy, well_xy, *, cone: bool):
    """Gradient + recharge mound, plus (``cone=True``) a Thiem cone clipped at 10 ft."""

    x, y = xy[:, 0], xy[:, 1]
    field = CHD_WEST - (CHD_WEST - CHD_EAST) * x / SIDE + 0.66 * 4.0 * (x / SIDE) * (1 - x / SIDE)
    if cone:
        transmissivity, radius = 3400.0, 1500.0
        r = np.maximum(np.hypot(x - well_xy[0], y - well_xy[1]), 10.0)
        field = field + np.where(
            r < radius, (-WELL_Q) / (2 * np.pi * transmissivity) * np.log(r / radius), 0.0)
    return field


def report_direction(label, old_vor, new_vor, old_heads, new_heads, overlay, well_xy, near_r):
    print(f"\n  ({label}) {old_vor.ncpl} -> {new_vor.ncpl} cells, against the "
          f"{new_vor.ncpl}-cell model's own solution")
    mapped, info = both_remaps(old_vor, new_vor, old_heads, overlay)
    new_xy = centres(new_vor)
    r_well = np.hypot(new_xy[:, 0] - well_xy[0], new_xy[:, 1] - well_xy[1])
    near = r_well <= near_r
    t = face_position(centres(old_vor), new_xy)
    print("      " + "; ".join(
        f"{m}: {info[m]['filled']} cell-layers filled by nearest, "
        f"{info[m]['seconds']:.3f} s" for m in METHODS)
        + f"; {near.sum()} new cells near the well")
    print(f"      {'method':7} {'lay':3} {'subset':9} {'n':>5} {'RMSE':>7} {'max|e|':>7} {'bias':>8}")
    for method in METHODS:
        for k in range(2):
            for subset, mask in (("all", None), ("near well", near), ("away", ~near)):
                s = stats(mapped[method][k], new_heads[k], mask)
                print(f"      {method:7} {k + 1:3} {subset:9} {s['n']:5d} {s['rmse']:7.4f} "
                      f"{s['max']:7.4f} {s['bias']:+8.4f}")
    print("      away from the well, mean |error| by position between old nodes "
          "(t = d1/d2: ~0 node, ~1 face):")
    bins = [(0.0, 0.33, "t<.33"), (0.33, 0.67, ".33-.67"), (0.67, 1.01, "t>.67")]
    for method in METHODS:
        for k in range(2):
            error = np.abs(mapped[method][k] - new_heads[k])
            cells = [f"{name} {error[(~near) & (t >= lo) & (t < hi)].mean():.4f}"
                     for lo, hi, name in bins]
            print(f"        {method:7} layer {k + 1}: " + "  ".join(cells))
    wc = well_cell(new_vor, well_xy)
    print(f"      new well cell {wc} (node {r_well[wc]:.0f} ft from the well), layer 2: "
          f"model {new_heads[1, wc]:.3f}, linear {mapped['linear'][1, wc]:.3f}, "
          f"area {mapped['area'][1, wc]:.3f}")
    for m in METHODS:
        c = int(np.argmax(np.abs(mapped[m][1] - new_heads[1])))
        print(f"      worst layer-2 cell, {m}: {c}, r = {r_well[c]:.0f} ft, t = {t[c]:.2f}")
    return mapped


def analytic_test(old_vor, new_vor, overlay, well_xy, near_r, label):
    old_xy, new_xy = centres(old_vor), centres(new_vor)
    t = face_position(old_xy, new_xy)
    near = np.hypot(new_xy[:, 0] - well_xy[0], new_xy[:, 1] - well_xy[1]) <= near_r
    for cone in (False, True):
        truth = analytic(new_xy, well_xy, cone=cone)
        mapped, _ = both_remaps(old_vor, new_vor, analytic(old_xy, well_xy, cone=cone)[None],
                                overlay)
        for method in METHODS:
            error = mapped[method][0] - truth
            s = stats(mapped[method][0], truth)
            print(f"      {label:12} {'+cone' if cone else 'smooth':6} {method:7} "
                  f"RMSE {s['rmse']:.5f}  max {s['max']:.4f}  bias {s['bias']:+.5f}  "
                  f"away RMSE {np.sqrt(np.mean(error[~near] ** 2)):.5f}  "
                  f"near RMSE {np.sqrt(np.mean(error[near] ** 2)):.4f}  |e| away at node "
                  f"{np.abs(error[(~near) & (t < 0.33)]).mean():.5f} / at face "
                  f"{np.abs(error[(~near) & (t > 0.67)]).mean():.5f}")


def main(workspace: Path) -> None:
    from myflopy.modflow.mf6.grid.geometry import cell_overlay

    grids = {tag: make_grid(workspace / "grids", tag) for tag in MAX_AREA}
    coarse, fine = grids["coarse"], grids["fine"]
    near_r = NEAR_WELL_WIDTHS * width(coarse)
    print(f"coarse: {coarse.ncpl} cells, median width {width(coarse):.1f} ft; fine: "
          f"{fine.ncpl} cells, median width {width(fine):.1f} ft; CHD cells per layer "
          f"{sum(map(len, chd_cells(coarse)))} / {sum(map(len, chd_cells(fine)))}; "
          f"near-well radius {near_r:.0f} ft")
    started = time.perf_counter()
    overlays = {"a": cell_overlay(coarse.gdf_vorPolys, fine.gdf_vorPolys),
                "b": cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)}
    print(f"two overlays: {time.perf_counter() - started:.3f} s")

    placements = {"off-node": WELL_TARGET, "on-node": on_node_well(coarse, fine)}
    for placement, well_xy in placements.items():
        offsets = {tag: float(np.hypot(*(centres(v)[well_cell(v, well_xy)] - well_xy)))
                   for tag, v in grids.items()}
        print(f"\n=== well {placement} at ({well_xy[0]:.1f}, {well_xy[1]:.1f}): its cell's "
              f"node is {offsets['coarse']:.1f} ft away on the coarse mesh, "
              f"{offsets['fine']:.1f} ft on the fine ===")
        base = {tag: solve(vor, workspace / f"run_{tag}", strt=STRT_FLAT, well_xy=well_xy)
                for tag, vor in grids.items()}
        for tag in grids:
            h = base[tag]["heads"]
            print(f"  {tag} solution: layer 1 {h[0].min():.3f}..{h[0].max():.3f}, layer 2 "
                  f"{h[1].min():.3f}..{h[1].max():.3f}")
        fine_from_coarse = report_direction(
            "a", coarse, fine, base["coarse"]["heads"], base["fine"]["heads"],
            overlays["a"], well_xy, near_r)
        coarse_from_fine = report_direction(
            "b", fine, coarse, base["fine"]["heads"], base["coarse"]["heads"],
            overlays["b"], well_xy, near_r)

        print("\n  (c) iterations (outer = CALLS TO NUMERICAL SOLUTION, inner = TOTAL "
              "ITERATIONS)")
        for tag, vor, mapped in (("fine", fine, fine_from_coarse),
                                 ("coarse", coarse, coarse_from_fine)):
            starts = {"flat 95 ft": STRT_FLAT, "linear remap": mapped["linear"],
                      "area remap": mapped["area"], "own solution": base[tag]["heads"]}
            for name, strt in starts.items():
                result = (base[tag] if name == "flat 95 ft"
                          else solve(vor, workspace / f"run_{tag}", strt=strt, well_xy=well_xy))
                drift = np.abs(result["heads"] - base[tag]["heads"]).max()
                print(f"      {tag:6} from {name:12}: outer {result['outer']:3d}, inner "
                      f"{result['inner']:4d}  (max |head - flat-start head| {drift:.1e})")

        print("\n  (d) an ANALYTIC point field (interpolation error only)")
        analytic_test(coarse, fine, overlays["a"], well_xy, near_r, "coarse->fine")
        analytic_test(fine, coarse, overlays["b"], well_xy, near_r, "fine->coarse")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(Path(sys.argv[1]))
    else:
        with tempfile.TemporaryDirectory(prefix="heads_remap_") as tmp:
            main(Path(tmp))
