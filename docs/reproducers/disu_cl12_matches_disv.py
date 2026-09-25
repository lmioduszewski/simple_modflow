"""A DISU written from ``grid.cl12``/``grid.hwva`` reproduces the DISV it came from (ledger 328).

    /home/lukem/python/envs/gw/bin/python docs/reproducers/disu_cl12_matches_disv.py [workspace]

Needs ``mf6`` on PATH (measured with 6.8.0.dev0+0418efe4). About 7 s. Writes nothing
outside the workspace it is given (a temporary directory when none is).

THE CLAIM. A DISV file holds no connection lengths: MF6 computes them itself, per SIDE --
the perpendicular from each cell's ``cell2d`` centre to the line through the two vertices
the cells share (``DisvGeom.f90``: ``shared_edge``, ``distance_normal``). DISU's ``CL12`` is
the same quantity supplied by hand ("the center of cell n [to] the shared face with each
adjacent m cell"). From 2025-07-02 until ledger 328 ``VoronoiGridPlus.cl12`` held the FULL
centre-to-centre spacing (between polygon centroids, projected on the face normal) on both
sides, so ``cl1 + cl2`` was about twice MF6's and every DISU myflopy wrote had about half the
horizontal conductance of the same grid as a DISV.

THE BENCH. Three meshes -- a real ``GridSpec.voronoi`` triangulation of a 600 ft square
(``max_area=12_000``, 92 cells: the ``triangulated_grid_finer`` fixture's construction),
``irregular_voronoi_grid(6, 6, 100 ft)`` (36 cells: the ``rectangular_voronoi`` family) and
8 x 8 perturbed quads whose DISV centres sit OFF their centroids (64 cells, ``asymmetric``).
One layer, top 100, bottom 0; K = 10**U(-1, 1) per cell (seed 7), so the two sides of a face
carry different K and the harmonic weighting reads cl1 and cl2 SEPARATELY, not just their
sum; CHD 95 ft on the cells whose centre is within 60 ft of the west edge, 90 ft on the
east; steady. Each mesh goes through BOTH of myflopy's DISU writers, each against a DISV
built by the same path:

  spec    ``simple_model_spec(SimpleModelConfig(grid_type="disu"))`` vs ``"disv"`` (NPF
          icelltype 1, NEWTON under-relaxation -- the writer's own settings);
  legacy  ``DisuGrid`` vs ``DisvGrid`` on a ``SimulationBase`` (icelltype 0).

and three more times through the spec writer with ``cl12`` replaced: by the PRE-328 value,
rebuilt here from the polygons exactly as ``build_disu_connectivity`` computed it then; by
each side given the OTHER cell's length (``swapped``: right in the sum, wrong per side); and
by half the distance between the two DISV centres on both sides (``half spacing``).

ONLY THE THIRD MESH CAN TELL THE SIDES APART. A Voronoi face is the perpendicular bisector of
its two generators, so on the first two meshes cl1 == cl2 to rounding -- the median
``|cl1 - cl2| / (cl1 + cl2)`` is 1.8e-16 and 3.5e-11 -- and the swapped and half-spacing
tables are the right table. On ``asymmetric`` it is 0.17. The canonical TESTING grid, which
the legacy-lake goldens use, is no better: its two sides differ by at most 5.1e-7 ft (4.3e-9
of their sum), printed last.

Measured (all three meshes, both writers; the connection tables in the DISV's and the
DISU's ``.grb`` files are identical, IA and JA):

  - heads: max |DISU - DISV| 1.9e-10 ft (triangulated), 2.1e-10 ft (irregular) and 8.5e-11
    ft (asymmetric) on heads spanning 90-95 ft; FLOW-JA-FACE max |DISU - DISV| / max |q|
    1.2e-10, 3.4e-11 and 6.8e-11; the CHD through-flow 401.8255 = 401.8255 ft3/d (spec) and
    432.5390 = 432.5390 (legacy) on the first mesh, 562.0607 and 604.0355 on the second,
    440.6019 and 473.9564 on the third -- the two writers' own models differ (icelltype),
    and each DISU matches its DISV to solver precision;
  - with the pre-328 cl12, through the spec writer: the CHD through-flow 0.511
    (triangulated), 0.501 (irregular) and 0.480 (asymmetric) of the DISV's, and per
    connection |q| / |dh| -- the conductance MF6 used, over connections with |dh| > 1e-3
    ft -- 0.515, 0.525 and 0.501 of the DISV's at the median, 0.41-0.64, 0.19-0.74 and
    0.21-0.81 across them (175, 75 and 98 connections); heads off by up to 0.10, 0.36 and
    0.37 ft;
  - the swapped and the half-spacing tables: on the two Voronoi meshes the SAME as the right
    one (max |dh| <= 2.5e-10 ft, FLOW-JA-FACE <= 1.2e-10 of its largest); on
    ``asymmetric``, max |dh| 0.41 ft and 0.24 ft, FLOW-JA-FACE off by 0.16 and 0.10 of its
    largest, through-flow 1.015 and 0.958. So a test on a Voronoi mesh alone cannot catch a
    table that puts each side's length on the wrong side (``tests/test_disu_matches_disv.py``
    runs this mesh too);
  - ``grid.calculate_distance`` -- the pre-2025 rule, polygon centroid to the nearest point
    of the shared boundary, still public -- against ``grid.cl12`` over every connection:
    0.39-1.49, 0.57-4.31 and 0.64-2.12 (medians 0.99, 0.98, 1.00) on the three meshes, and
    0.48-15.89 (median 0.99) on the canonical testing grid. Not MF6's number (ledger 328).
"""

from __future__ import annotations

import contextlib
import io
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simplemodel import SimpleModelConfig, simple_model_spec

TOP, BOTTOM, WEST_HEAD, EAST_HEAD, EDGE = 100.0, 0.0, 95.0, 90.0, 60.0


def triangulated(root: Path):
    """``GridSpec.voronoi`` over a 600 ft square, 92 cells (``triangulated_grid_finer``)."""

    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"name": ["d"]}, geometry=[box(0.0, 0.0, 600.0, 600.0)], crs=2927).to_file(boundary)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=12_000.0
        ).resolve(workspace=root / "grid")


def irregular():
    """``irregular_voronoi_grid(6, 6, 100 ft)``: the ``rectangular_voronoi`` family, 36 cells."""

    return irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0)


def asymmetric():
    """8 x 8 quads of 100 ft whose two sides of a face DIFFER (64 cells).

    Interior vertices moved up to 25 ft each way (seed 3), and each DISV centre a random
    interior point of its quad -- a convex combination of its corners, OFF the centroid.
    On a Voronoi mesh the face is the perpendicular bisector of the two generators, so
    cl1 == cl2 to rounding and nothing can tell which side a table put where; here the
    median ``|cl1 - cl2| / (cl1 + cl2)`` is 0.17.
    """

    rng = np.random.default_rng(3)
    n, size = 8, 100.0
    edge = np.arange(n + 1) * size
    x, y = np.meshgrid(edge, edge)
    interior = lambda a: (a > 0.0) & (a < n * size)  # noqa: E731 -- the boundary stays straight
    x = x + np.where(interior(x), rng.uniform(-25.0, 25.0, x.shape), 0.0)
    y = y + np.where(interior(y), rng.uniform(-25.0, 25.0, y.shape), 0.0)
    verts = np.column_stack([x.ravel(), y.ravel()])
    iverts, centres = [], []
    for i in range(n):  # row, in y
        for j in range(n):  # column, in x
            ring = [i * (n + 1) + j, (i + 1) * (n + 1) + j, (i + 1) * (n + 1) + j + 1, i * (n + 1) + j + 1]
            iverts.append(ring)  # clockwise with y up, as MF6 wants
            weights = rng.dirichlet([1.0, 1.0, 1.0, 1.0]) * 0.6 + 0.1
            centres.append(verts[ring].T @ (weights / weights.sum()))
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centres))


def side_asymmetry(grid) -> float:
    """The median ``|cl1 - cl2| / (cl1 + cl2)`` over the connections."""

    lengths = cl12_by_pair(grid, np.asarray(grid.cl12))
    return float(np.median([abs(a - b) / (a + b) for (n, m), a in lengths.items() if n < m
                            for b in [lengths[(m, n)]]]))


def centroid_to_face(grid) -> np.ndarray:
    """``grid.calculate_distance`` over every connection, as a fraction of ``grid.cl12``.

    The pre-2025 helper, still public: polygon centroid to the nearest point of the shared
    boundary. Not MF6's number, and nothing that writes a model reads it.
    """

    lengths = cl12_by_pair(grid, np.asarray(grid.cl12))
    return np.array([grid.calculate_distance(grid.gdf_vorPolys, n, m) / length
                     for (n, m), length in lengths.items()])


def canonical_sides(root: Path) -> None:
    """How far apart the two sides of a face are on the canonical TESTING grid.

    The grid every lake golden in ``tests/test_excavated_lake.py`` but one is on. Its faces
    are Voronoi bisectors, so a legacy wall that read the NEIGHBOUR's side instead of the
    lake cell's would move by no more than this.
    """

    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model

    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        model = build_canonical_model(root / "canonical", config=CanonicalModelConfig.testing())
    lengths = cl12_by_pair(model.vor, np.asarray(model.vor.cl12))
    gaps = np.array([(abs(a - lengths[(m, n)]), abs(a - lengths[(m, n)]) / (a + lengths[(m, n)]))
                     for (n, m), a in lengths.items() if n < m])
    ratio = centroid_to_face(model.vor)
    print(f"\ncanonical testing grid: {model.vor.ncpl} cells, |cl1 - cl2| at most "
          f"{gaps[:, 0].max():.1e} ft, {gaps[:, 1].max():.1e} of cl1 + cl2; calculate_distance / cl12 "
          f"{ratio.min():.2f}-{ratio.max():.2f} (median {np.median(ratio):.2f})")


def cl12_by_pair(grid, cl12) -> dict:
    """``{(n, m): cl12 at row n's entry for m}``."""

    iac, ja = np.asarray(grid.iac), np.asarray(grid.ja)
    ia = np.concatenate([[0], np.cumsum(iac)])
    return {(n, int(ja[pos])): float(cl12[pos])
            for n in range(iac.size) for pos in range(ia[n] + 1, ia[n + 1])}


def swapped_cl12(grid) -> np.ndarray:
    """Each side given the OTHER cell's length: right in their sum, wrong per side."""

    lengths = cl12_by_pair(grid, np.asarray(grid.cl12))
    iac, ja = np.asarray(grid.iac), np.asarray(grid.ja)
    ia = np.concatenate([[0], np.cumsum(iac)])
    out = np.zeros(ja.size)
    for n in range(iac.size):
        for pos in range(ia[n] + 1, ia[n + 1]):
            out[pos] = lengths[(int(ja[pos]), n)]
    return out


def half_spacing_cl12(grid) -> np.ndarray:
    """Half the distance between the two DISV centres, on both sides."""

    centre = np.array([(row[1], row[2]) for row in grid.get_disv_gridprops()["cell2d"]])
    iac, ja = np.asarray(grid.iac), np.asarray(grid.ja)
    ia = np.concatenate([[0], np.cumsum(iac)])
    out = np.zeros(ja.size)
    for n in range(iac.size):
        for pos in range(ia[n] + 1, ia[n + 1]):
            out[pos] = 0.5 * float(np.hypot(*(centre[int(ja[pos])] - centre[n])))
    return out


def pre_328_cl12(grid) -> np.ndarray:
    """``cl12`` as it was from 2025-07-02 to ledger 328: centroid spacing on the face normal."""

    polygons = grid.gdf_vorPolys.geometry.values
    centroids = np.array([(p.centroid.x, p.centroid.y) for p in polygons])
    iac, ja = np.asarray(grid.iac), np.asarray(grid.ja)
    ia = np.concatenate([[0], np.cumsum(iac)])
    old = np.zeros(ja.size)
    for n in range(iac.size):
        for pos in range(ia[n] + 1, ia[n + 1]):
            m = int(ja[pos])
            face = polygons[n].boundary.intersection(polygons[m].boundary)
            segment = max(getattr(face, "geoms", [face]), key=lambda s: s.length)
            (x1, y1), (x2, y2) = segment.coords[0], segment.coords[-1]
            ex, ey = x2 - x1, y2 - y1
            length = np.hypot(ex, ey)
            d = centroids[m] - centroids[n]
            old[pos] = abs(d[0] * -ey / length + d[1] * ex / length)
    return old


def bench(grid):
    """K, the CHD cells and their heads, from the DISV centres."""

    xy = np.array([(row[1], row[2]) for row in grid.get_disv_gridprops()["cell2d"]])
    k = 10.0 ** np.random.default_rng(7).uniform(-1.0, 1.0, grid.ncpl)
    west = np.flatnonzero(xy[:, 0] < xy[:, 0].min() + EDGE)
    east = np.flatnonzero(xy[:, 0] > xy[:, 0].max() - EDGE)
    return k, [int(c) for c in west] + [int(c) for c in east], \
        [WEST_HEAD] * west.size + [EAST_HEAD] * east.size


def outputs(workspace: Path, name: str, grb: str):
    """Heads, FLOW-JA-FACE, CHD through-flow and the connection table MF6 wrote."""

    heads = flopy.utils.HeadFile(workspace / f"{name}.hds").get_data().ravel()
    budget = flopy.utils.CellBudgetFile(workspace / f"{name}.cbc", precision="double")
    flowja = np.asarray(budget.get_data(text="FLOW-JA-FACE")[0]).ravel()
    chd = budget.get_data(text="CHD")[0]
    inflow = float(np.asarray(chd["q"])[np.asarray(chd["q"]) > 0].sum())
    table = flopy.mf6.utils.MfGrdFile(str(workspace / grb), verbose=False)
    return heads, flowja, inflow, (np.asarray(table.ia), np.asarray(table.ja))


def run_spec(root: Path, grid, grid_type: str, tag: str, cl12=None):
    """``simple_model_spec`` with ``grid_type``; ``cl12`` replaces the grid's (the pre-328 case)."""

    k, cells, heads = bench(grid)
    saved = grid._cl12
    if cl12 is not None:
        grid._cl12 = cl12
    try:
        spec = simple_model_spec(SimpleModelConfig(
            vor=grid, name="m", grid_type=grid_type, top=TOP,
            bottom=BOTTOM if grid_type == "disu" else [BOTTOM], k=list(k),
            initial_heads=[92.5] * grid.ncpl, boundary_mode="chd", boundary_cells=cells,
            boundary_head=heads, save_specific_discharge=False, sto_transient={},
            output_save_record=(("HEAD", "ALL"), ("BUDGET", "ALL")), complexity="COMPLEX",
        ))
        workspace = root / tag
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            simulation = spec.build_flopy(workspace).sim
            simulation.write_simulation(silent=True)
            ok, report = simulation.run_simulation(silent=True)
    finally:
        grid._cl12 = saved
    assert ok, "\n".join(report[-20:])
    return outputs(workspace, "m", f"m.{grid_type}.grb")


def run_legacy(root: Path, grid, grid_type: str, tag: str):
    """``DisuGrid`` / ``DisvGrid`` on a ``SimulationBase``, confined (icelltype 0)."""

    from myflopy.modflow.mf6.simulation import SimulationBase
    from myflopy.modflow.mf6.simulation.discretization import (
        DisuGrid,
        DisvGrid,
        TemporalDiscretization,
    )
    from myflopy.modflow.mf6.simulation.packages import (
        CHD,
        InitialConditions,
        KFlow,
        OutputControl,
        Storage,
    )

    k, cells, heads = bench(grid)
    workspace = root / tag
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        model = SimulationBase(name="m", mf_folder_path=workspace, vor=grid, nper=1)
        model.nlay = 1
        if grid_type == "disu":
            DisuGrid(vor=grid, model=model, top=[TOP] * grid.ncpl, bottom=[BOTTOM] * grid.ncpl)
            cellids = cells
        else:
            DisvGrid(vor=grid, model=model, top=[TOP] * grid.ncpl, bottom=[[BOTTOM] * grid.ncpl], nlay=1)
            cellids = [(0, c) for c in cells]
        TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
        InitialConditions(model=model, vor=grid, nlay=1, strt=[92.5] * grid.ncpl)
        KFlow(model=model, k=list(k), icelltype=0, save_specific_discharge=False)
        Storage(model=model, sto_steady={0: True}, sto_transient={})
        OutputControl(model=model)
        CHD(model=model, stress_period_data={0: [[c, h] for c, h in zip(cellids, heads, strict=True)]})
        ok, report = model.run_simulation()
    assert ok, "\n".join(report[-20:])
    run_dir = workspace / "m"
    budget_name = next(run_dir.glob("*.cbc"))
    heads_name = next(run_dir.glob("*.hds"))
    heads_out = flopy.utils.HeadFile(heads_name).get_data().ravel()
    budget = flopy.utils.CellBudgetFile(budget_name, precision="double")
    flowja = np.asarray(budget.get_data(text="FLOW-JA-FACE")[-1]).ravel()
    chd = budget.get_data(text="CHD")[-1]
    inflow = float(np.asarray(chd["q"])[np.asarray(chd["q"]) > 0].sum())
    table = flopy.mf6.utils.MfGrdFile(str(next(run_dir.glob(f"*.{grid_type}.grb"))), verbose=False)
    return heads_out, flowja, inflow, (np.asarray(table.ia), np.asarray(table.ja))


def compare(disv, disu):
    """Max head difference, max relative FLOW-JA-FACE difference, through-flow ratio, conductance ratios."""

    (h1, q1, i1, (ia1, ja1)), (h2, q2, i2, (ia2, ja2)) = disv, disu
    assert np.array_equal(ia1, ia2) and np.array_equal(ja1, ja2), "connection tables differ"
    row = np.repeat(np.arange(ia1.size - 1), np.diff(ia1))
    off = row < ja1  # each connection once
    dh1 = np.abs(h1[row[off]] - h1[ja1[off]])
    dh2 = np.abs(h2[row[off]] - h2[ja2[off]])
    usable = (dh1 > 1e-3) & (dh2 > 1e-3)
    conductance = (np.abs(q2[off]) / dh2)[usable] / (np.abs(q1[off]) / dh1)[usable]
    return (float(np.abs(h2 - h1).max()), float(np.abs(q2 - q1).max() / np.abs(q1).max()),
            i2 / i1, conductance)


def main(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    print(f"mf6: {shutil.which('mf6')}")
    meshes = (("triangulated", triangulated(root)), ("irregular", irregular()), ("asymmetric", asymmetric()))
    for label, grid in meshes:
        grid.get_disu_connectivity(validate=True)
        heads_disv = run_spec(root, grid, "disv", f"{label}_spec_disv")
        ratio = centroid_to_face(grid)
        print(f"\n{label}: {grid.ncpl} cells, heads {heads_disv[0].min():.2f}-{heads_disv[0].max():.2f} ft, "
              f"median |cl1 - cl2| / (cl1 + cl2) {side_asymmetry(grid):.2g}, calculate_distance / cl12 "
              f"{ratio.min():.2f}-{ratio.max():.2f} (median {np.median(ratio):.2f})")
        for writer, disv, disu in (
            ("spec  ", heads_disv, run_spec(root, grid, "disu", f"{label}_spec_disu")),
            ("legacy", run_legacy(root, grid, "disv", f"{label}_legacy_disv"),
             run_legacy(root, grid, "disu", f"{label}_legacy_disu")),
        ):
            dh, dq, ratio, _ = compare(disv, disu)
            print(f"  {writer} DISU vs DISV: max|dh| {dh:.1e} ft, max|dq|/max|q| {dq:.1e}, "
                  f"through-flow {disu[2]:.4f} / {disv[2]:.4f} ft3/d = {ratio:.9f}")
        old = run_spec(root, grid, "disu", f"{label}_spec_disu_pre328", cl12=pre_328_cl12(grid))
        dh, dq, ratio, conductance = compare(heads_disv, old)
        print(f"  pre-328 cl12:  through-flow {ratio:.3f} of the DISV's, conductance |q|/|dh| "
              f"median {np.median(conductance):.3f} (range {conductance.min():.2f}-"
              f"{conductance.max():.2f}, {conductance.size} connections), max|dh| {dh:.2f} ft")
        for name, table in (("swapped sides", swapped_cl12(grid)), ("half spacing ", half_spacing_cl12(grid))):
            dh, dq, ratio, _ = compare(heads_disv, run_spec(root, grid, "disu", f"{label}_{name[:4]}", cl12=table))
            print(f"  {name}:  max|dh| {dh:.1e} ft, max|dq|/max|q| {dq:.1e}, through-flow {ratio:.6f}")
    canonical_sides(root)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(Path(sys.argv[1]))
    else:
        with tempfile.TemporaryDirectory() as tmp:
            main(Path(tmp))
