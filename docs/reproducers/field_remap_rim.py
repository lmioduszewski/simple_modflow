"""Where ``remap``'s linear method cannot interpolate: the rim, a bridging triangle, the clamp.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/field_remap_rim.py

No MF6; about 3 s. The companion of ``heads_remap_linear_vs_area.py``, which settled the
DEFAULT (linear) on two ``GridSpec.voronoi`` meshes -- whose boundary generators lie ON the
boundary, so every new centre is inside the hull of the old ones and none of what this
script measures ever happens there.

THE SETUP. Two meshes from ``rectangular_voronoi`` (``canonical_example.py``) over the same
1,200 ft square: 6x6 (200 ft, OLD) and 12x12 (100 ft, NEW). Their generators are
JITTERED, so the old boundary cells are not centred on the boundary, and the hull of the old
centres leaves a band of the square uncovered -- 38 of the 144 new centres sit in it (a
structured, DIS-like mesh measured 9% in the scouting). Analytic fields, evaluated exactly at
both meshes' MF6 centres (``grid.points``), so every error is the remap's own.

PART 1 -- the RIM. A new centre outside the old hull cannot be interpolated. Candidates, all
LOCAL (the value comes from the valid old cell covering most of the new cell):

  area        the area mean of the old cells under the new cell (the first version)
  none        that old cell's value moved along its own least-squares gradient (a plane
              through it and its vertex-touching neighbours), unclamped
  local       that, clamped to the range of the cells the plane was fitted to
  layer       that, clamped to the range of the whole old layer (SHIPPED)

on four fields: a PLANE, a MOUND, a sharp FRONT (exp(-x/60), a plume edge) and a RIDGE
(flat at the centre, as at a no-flow divide).

PART 2 -- a BRIDGING triangle. Old field 10 west of x = 600 and 20 east of it, with the cells
in 500 < x < 700 INACTIVE. A Delaunay triangle of the valid old centres spans the strip; no
new cell may read a value between the two, since no old cell holds one.

MEASURED (2026-10-15), PART 1 -- rim cells, max / mean |error|:

  field   area            none            local           layer (shipped)
  plane   2.786 / 1.208   0.000 / 0.000   2.782 / 0.380   2.782 / 0.338
  mound   2.701 / 0.902   0.808 / 0.332   0.937 / 0.348   0.865 / 0.316
  front   0.461 / 0.074   0.406 / 0.065   0.418 / 0.074   0.406 / 0.067
  ridge   4.269 / 1.768   1.281 / 0.467   3.564 / 1.162   1.455 / 0.530

  * Unclamped is the most accurate and it is WRONG in kind: on the front it returns -0.055
    from a field whose least value is 4.7e-9 -- a negative concentration. A plane fitted
    across a sharp front extrapolates past it. Interpolation inside a triangle is a
    weighted average and cannot do that.
  * Clamping to the FITTED cells is as bad as the area mean on a plane (2.782 against
    2.786): at the rim a plane lies past its neighbours by definition. Clamping to the
    LAYER clips only the extreme corners -- mean 0.338 against 0.380 -- and was never the
    worse of the two clamps on any field.

MEASURED, PART 2 (the 12x12 field onto the 6x6 mesh): with the local-simplex test, all 36 new
cells read 10.0 or 20.0 exactly. Without it, 6 read 17.21 to 17.96 -- the cells STRADDLING
the strip's edge, which have valid old cells under one side and a centre inside a triangle
reaching across. (The scouting measured the same failure on a 10 ft offset across a
two-column inactive strip of a structured mesh: up to 6.65 ft wrong.)

PART 3 -- WHICH POINT a value sits at. MF6 solves each head at its cell's centre (the
``xc``/``yc`` the DISV writes, ``grid.points``), which on a Voronoi mesh is not the polygon
centroid. The same remap with both meshes' centres swapped for their centroids.

MEASURED, PART 3 -- mean / max |error|:

  6x6 -> 12x12 (old centre-to-centroid: median 37.8, max 78.7, cells ~191 wide)
    plane   centres 0.0893 / 2.782   centroids 0.5489 / 2.782    (6.1x the mean)
    mound   centres 0.3622 / 2.073   centroids 0.5410 / 3.126    (1.5x)
    ridge   centres 0.4053 / 1.455   centroids 0.6551 / 3.892    (1.6x)
  12x12 -> 6x6 (median 17.5, max 55.5, cells ~95 wide)
    plane   centres 0.0000 / 0.000   centroids 0.5757 / 1.450
    mound   centres 0.1045 / 0.533   centroids 0.5102 / 1.613    (4.9x)
    ridge   centres 0.1460 / 0.501   centroids 0.5312 / 2.572    (3.6x)

PART 4 -- a CHAIN. One band of valid cells (the rest inactive) of a 12x12 mesh onto a 24x24
one: each cell's valid neighbours lie close to ONE line, so the slope across the band is set
by their jitter. The first version fitted it whenever the normal matrix's determinant
cleared 1e-9 of its scale; the shipped rule drops a direction the neighbours span less than
a tenth as widely as the other (``_PLANE_SPREAD = 0.01``).

MEASURED, PART 4 -- reached cells, max / mean |error|:

  band 560-640 (18 cells)   shipped 1.622 / 0.152   first version  1.622 / 0.160   area 2.891 / 0.576
  band 300-380 (12 cells)   shipped 2.926 / 0.305   first version  9.741 / 0.866   area 1.983 / 0.557
  band 820-900 (12 cells)   shipped 1.966 / 0.380   first version 22.841 / 1.917   area 2.488 / 0.627

NOT SHOWN. A real model; heads near a pumping well at the rim; a mesh whose boundary
generators sit on the boundary (no rim at all -- ``heads_remap_linear_vs_area.py``).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.spatial import Delaunay

import myflopy.modflow.mf6.grid.remap as remap
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi

FIELDS = {
    "plane": lambda xy: 50 + 0.01 * xy[:, 0] - 0.02 * xy[:, 1],
    "mound": lambda xy: 100 - 20 * np.exp(-((xy[:, 0] - 600) ** 2 + (xy[:, 1] - 600) ** 2) / 2e5),
    "front": lambda xy: np.exp(-(xy[:, 0] / 60.0)),
    "ridge": lambda xy: 100 - 0.00005 * (xy[:, 0] - 600) ** 2,
}


def _grid(n, cell):
    return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1, cell_size=cell))


def _extend_with(mode):
    """A replacement for ``_Triangulation.extend`` applying one of the rim rules.

    Every rule moves the same anchor along the same gradient (``_Triangulation.moved``);
    they differ only in the clamp.
    """

    shipped = remap._Triangulation.extend

    def extend(self, values_at_ids, anchors, points):
        values = np.asarray(values_at_ids, dtype=float)
        if mode == "layer":
            return shipped(self, values, anchors, points)
        moved = self.moved(values, anchors, points)
        if mode == "none":
            return moved
        here, there = self.pairs[0], self.pairs[1]
        low, high = values.copy(), values.copy()
        np.minimum.at(low, here, values[there])
        np.maximum.at(high, here, values[there])
        return np.clip(moved, low[anchors], high[anchors])

    return extend


def part_1(old, new):
    rim = Delaunay(old.points).find_simplex(new.points) < 0
    print(f"PART 1: {int(rim.sum())} of {new.ncpl} new centres outside the old hull")
    print(f"  {'field':6s} {'area':>15s} {'none':>15s} {'local':>15s} {'layer':>15s}")
    shipped = remap._Triangulation.extend
    for name, field in FIELDS.items():
        data, truth = field(old.points), field(new.points)
        cells = []
        area, _ = remap.remap_field(data[None, :], old, new, method="area")
        cells.append(np.abs(area[0] - truth)[rim])
        lowest = np.inf
        for mode in ("none", "local", "layer"):
            remap._Triangulation.extend = _extend_with(mode)
            try:
                out, _ = remap.remap_field(data[None, :], old, new)
            finally:
                remap._Triangulation.extend = shipped
            cells.append(np.abs(out[0] - truth)[rim])
            if mode == "none":
                lowest = out.min()
        row = "  ".join(f"{e.max():6.3f} / {e.mean():5.3f}" for e in cells)
        print(f"  {name:6s} {row}")
        if name == "front":
            print(f"         unclamped least value {lowest:.3g} from a field whose least is {data.min():.3g}")


def part_2(old_fine, new):
    x = np.asarray(old_fine.points, dtype=float)[:, 0]
    field = np.where(x < 600.0, 10.0, 20.0)
    field[(x > 500.0) & (x < 700.0)] = 1.0e30
    bounds = new.gdf_vorPolys.geometry.bounds
    sides = (bounds.maxx.to_numpy() <= 500.0) | (bounds.minx.to_numpy() >= 700.0)
    print("PART 2: field 10 | inactive strip 500-700 | 20")
    shipped = remap._Triangulation.__init__

    def unfiltered(self, old_xy, ids, new_xy, touch):
        shipped(self, old_xy, ids, new_xy, None)

    for label, init in (("with the local-simplex test", shipped), ("without it", unfiltered)):
        remap._Triangulation.__init__ = init
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                out, _ = remap.remap_field(field[None, :], old_fine, new)
        finally:
            remap._Triangulation.__init__ = shipped
        between = ~np.isin(np.round(out[0], 9), [10.0, 20.0])
        print(f"  {label:28s}: {int(between.sum())} of {new.ncpl} new cells read a value "
              f"between the two sides ({int((between & sides).sum())} of them wholly on "
              "one side)" + (f", from {out[0][between].min():.2f} to "
                              f"{out[0][between].max():.2f}" if between.any() else ""))


class _OnCentroids:
    """``grid`` with its cell centres swapped for its polygon CENTROIDS."""

    def __init__(self, grid):
        self.gdf_vorPolys = grid.gdf_vorPolys
        self.iverts = grid.iverts
        centres = grid.gdf_vorPolys.geometry.centroid
        self.points = np.column_stack((centres.x.to_numpy(), centres.y.to_numpy()))


def part_3(old, new):
    """MF6's cell centres against polygon centroids, as the points a value sits at."""

    print("PART 3: linear on MF6 cell centres vs on polygon centroids (mean / max |error|)")
    offsets = np.hypot(*(np.asarray(old.points) - _OnCentroids(old).points).T)
    print(f"  old centre-to-centroid distance: median {np.median(offsets):.1f}, "
          f"max {offsets.max():.1f} (cells {np.sqrt(old.gdf_vorPolys.area.median()):.0f} wide)")
    for name in ("plane", "mound", "ridge"):
        field = FIELDS[name]
        data, truth = field(old.points), field(new.points)
        rows = []
        for label, source, target in (("centres", old, new),
                                      ("centroids", _OnCentroids(old), _OnCentroids(new))):
            out, _ = remap.remap_field(data[None, :], source, target)
            error = np.abs(out[0] - truth)
            rows.append(f"{label} {error.mean():.4f} / {error.max():.4f}")
        print(f"  {name:6s} " + "   ".join(rows))


def part_4(old, new):
    """A CHAIN of valid cells: a gradient fitted to neighbours along one line."""

    print("PART 4: one band of valid cells, h = 100 - 30 (x / 1200)^2 -- max / mean |error|")
    old_xy, new_xy = np.asarray(old.points), np.asarray(new.points)
    for low, high in ((560.0, 640.0), (300.0, 380.0), (820.0, 900.0)):
        band = (old_xy[:, 1] > low) & (old_xy[:, 1] < high)
        field = np.full(old.ncpl, 1.0e30)
        field[band] = 100.0 - 30.0 * (old_xy[band, 0] / 1200.0) ** 2
        under = old.gdf_vorPolys.geometry[band].union_all()
        reached = new.gdf_vorPolys.geometry.intersection(under).area.to_numpy() > 1e-6
        truth = 100.0 - 30.0 * (new_xy[reached, 0] / 1200.0) ** 2
        rows = []
        for label, spread in (("shipped", remap._PLANE_SPREAD), ("first version", 1e-9)):
            shipped, remap._PLANE_SPREAD = remap._PLANE_SPREAD, spread
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    out, _ = remap.remap_field(field[None, :], old, new)
            finally:
                remap._PLANE_SPREAD = shipped
            error = np.abs(out[0][reached] - truth)
            rows.append(f"{label} {error.max():6.3f} / {error.mean():5.3f}")
        area, _ = remap.remap_field(field[None, :], old, new, method="area")
        error = np.abs(area[0][reached] - truth)
        rows.append(f"area {error.max():6.3f} / {error.mean():5.3f}")
        print(f"  band {low:.0f}-{high:.0f} ({int(band.sum())} cells): " + "   ".join(rows))


if __name__ == "__main__":
    coarse, fine = _grid(6, 200.0), _grid(12, 100.0)
    part_1(coarse, fine)
    part_2(fine, coarse)
    part_3(coarse, fine)
    part_3(fine, coarse)
    part_4(fine, _grid(24, 50.0))
