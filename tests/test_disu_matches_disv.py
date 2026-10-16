"""A DISU written from ``grid.cl12``/``grid.hwva`` reproduces the DISV it came from (ledger 328).

A DISV file carries no connection lengths: MF6 computes them itself, per SIDE, from each
cell's ``cell2d`` centre to the line through the two vertices the cells share. A DISU's
``CL12`` is that same number supplied by hand. From 2025-07-02 until ledger 328
``VoronoiGridPlus.cl12`` held the full centre-to-centre spacing on both sides, so every DISU
myflopy wrote carried about half the DISV's horizontal conductance
(``docs/reproducers/disu_cl12_matches_disv.py`` measured 0.501 and 0.511 of the through-flow).

What these pin, in MF6 (``slow``): the same grid, K field and CHD written as a DISV and as a
DISU through BOTH of myflopy's DISU writers -- ``simple_model_spec(grid_type="disu")`` and
the legacy ``DisuGrid`` -- give the same connection table, heads and ``FLOW-JA-FACE`` to
solver precision, on a real ``GridSpec.voronoi`` triangulation, an
``irregular_voronoi_grid`` and perturbed quads whose DISV centres sit off their centroids.
K varies by cell (10**U(-1, 1)), so MF6's harmonic weighting reads ``cl1`` and ``cl2``
SEPARATELY -- but only the THIRD mesh gives them different values. A Voronoi face is the
perpendicular bisector of its two generators, so on the first two ``cl1 == cl2`` to rounding
(median ``|cl1 - cl2| / (cl1 + cl2)`` 1.8e-16 and 3.5e-11) and a table with the sides
SWAPPED, or half the centre spacing on both, is the right table there; on the quads (0.17)
those two miss the DISV by 0.41 and 0.24 ft (``docs/reproducers/disu_cl12_matches_disv.py``).
And, without MF6, that the table is MF6's rule: per-side, from the DISV centre (not the
polygon centroid), by hand on every connection of the irregular grid and the quads; and that
``calculate_distance``, the pre-2025 rule still public, says it is not that number.
"""

from __future__ import annotations

import builtins
import contextlib
import io
import warnings

import flopy
import numpy as np
import pytest

from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.connectivity import DisvFaces
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simplemodel import SimpleModelConfig, simple_model_spec

TOP, BOTTOM = 100.0, 0.0


def _bench(grid):
    """K per cell, and CHD 95 ft within 60 ft of the west edge, 90 ft of the east."""

    xy = np.array([(row[1], row[2]) for row in grid.get_disv_gridprops()["cell2d"]])
    k = 10.0 ** np.random.default_rng(7).uniform(-1.0, 1.0, grid.ncpl)
    west = np.flatnonzero(xy[:, 0] < xy[:, 0].min() + 60.0)
    east = np.flatnonzero(xy[:, 0] > xy[:, 0].max() - 60.0)
    cells = [int(c) for c in west] + [int(c) for c in east]
    return k, cells, [95.0] * west.size + [90.0] * east.size


def _results(folder, grid_type):
    heads = flopy.utils.HeadFile(next(folder.glob("*.hds"))).get_data().ravel()
    budget = flopy.utils.CellBudgetFile(next(folder.glob("*.cbc")), precision="double")
    flowja = np.asarray(budget.get_data(text="FLOW-JA-FACE")[-1]).ravel()
    table = flopy.mf6.utils.MfGrdFile(str(next(folder.glob(f"*.{grid_type}.grb"))), verbose=False)
    return heads, flowja, np.asarray(table.ia), np.asarray(table.ja)


def _run_spec(workspace, grid, grid_type):
    k, cells, heads = _bench(grid)
    spec = simple_model_spec(SimpleModelConfig(
        vor=grid, name="m", grid_type=grid_type, top=TOP,
        bottom=BOTTOM if grid_type == "disu" else [BOTTOM], k=list(k),
        initial_heads=[92.5] * grid.ncpl, boundary_mode="chd", boundary_cells=cells,
        boundary_head=heads, save_specific_discharge=False, sto_transient={},
        output_save_record=(("HEAD", "ALL"), ("BUDGET", "ALL")), complexity="COMPLEX",
    ))
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        simulation = spec.build_flopy(workspace).sim
        simulation.write_simulation(silent=True)
        ok, report = simulation.run_simulation(silent=True)
    assert ok, "\n".join(report[-20:])
    return _results(workspace, grid_type)


def _run_legacy(workspace, grid, grid_type):
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

    k, cells, heads = _bench(grid)
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
    return _results(workspace / "m", grid_type)


@pytest.fixture(scope="module")
def irregular():
    return irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0)


@pytest.fixture(scope="module")
def asymmetric():
    """8 x 8 quads of 100 ft whose two sides of a face DIFFER (the reproducer's ``asymmetric``).

    Interior vertices moved up to 25 ft each way, and each DISV centre a random interior
    point of its quad (a convex combination of its corners), off the centroid -- so the
    distance from each centre to a shared face is its own, which no Voronoi mesh gives.
    """

    rng = np.random.default_rng(3)
    n, size = 8, 100.0
    edge = np.arange(n + 1) * size
    x, y = np.meshgrid(edge, edge)
    x = x + np.where((x > 0.0) & (x < n * size), rng.uniform(-25.0, 25.0, x.shape), 0.0)
    y = y + np.where((y > 0.0) & (y < n * size), rng.uniform(-25.0, 25.0, y.shape), 0.0)
    verts = np.column_stack([x.ravel(), y.ravel()])
    iverts, centres = [], []
    for i in range(n):
        for j in range(n):
            ring = [i * (n + 1) + j, (i + 1) * (n + 1) + j, (i + 1) * (n + 1) + j + 1, i * (n + 1) + j + 1]
            iverts.append(ring)
            weights = rng.dirichlet([1.0, 1.0, 1.0, 1.0]) * 0.6 + 0.1
            centres.append(verts[ring].T @ (weights / weights.sum()))
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centres))


_MESHES = {"triangulated": "triangulated_grid_finer", "irregular": "irregular", "asymmetric": "asymmetric"}


@pytest.mark.slow
@pytest.mark.parametrize("writer", ["spec", "legacy"])
@pytest.mark.parametrize("mesh", list(_MESHES))
def test_mf6_a_disu_from_cl12_reproduces_its_disv(mesh, writer, request, tmp_path):
    grid = request.getfixturevalue(_MESHES[mesh])
    grid.get_disu_connectivity(validate=True)
    run = _run_spec if writer == "spec" else _run_legacy
    heads_v, flow_v, ia_v, ja_v = run(tmp_path / "disv", grid, "disv")
    heads_u, flow_u, ia_u, ja_u = run(tmp_path / "disu", grid, "disu")

    assert np.array_equal(ia_u, ia_v) and np.array_equal(ja_u, ja_v), "MF6's own connection table"
    assert heads_v.max() - heads_v.min() > 4.0, "the bench carries flow"
    # Measured 1.9e-10 ft and 1.2e-10 (triangulated), 2.1e-10 and 3.4e-11 (irregular),
    # 8.5e-11 and 6.8e-11 (asymmetric). The pre-328 cl12 missed by 0.10 / 0.36 / 0.37 ft and
    # about half the flow; on the asymmetric quads the sides swapped missed by 0.41 ft and
    # half the centre spacing by 0.24 ft, where on the two Voronoi meshes both pass.
    assert np.abs(heads_u - heads_v).max() < 1e-6
    assert np.abs(flow_u - flow_v).max() < 1e-6 * np.abs(flow_v).max()


@pytest.mark.parametrize("mesh", ["irregular", "asymmetric"])
def test_cl12_is_each_cells_own_side_measured_from_its_disv_centre(mesh, request):
    """Not the polygon centroid (MF6 never sees it), not the spacing between centres, and
    not the other cell's side -- which only the asymmetric quads can tell apart."""

    grid = request.getfixturevalue(mesh)
    iac, ja, cl12, hwva, _ = grid.get_disu_connectivity(validate=True)
    props = grid.get_disv_gridprops()
    centre = np.array([(row[1], row[2]) for row in props["cell2d"]])
    faces = DisvFaces(props)
    ia = np.concatenate([[0], np.cumsum(iac)])
    sides = []
    for n in range(grid.ncpl):
        for pos in range(ia[n] + 1, ia[n + 1]):
            m = int(ja[pos])
            v1, v2 = faces.shared_edge(n, m)
            (x1, y1), (x2, y2) = props["vertices"][v1][1:], props["vertices"][v2][1:]
            # distance from n's DISV centre to the LINE through the shared edge, by hand
            normal = np.array([y2 - y1, x1 - x2]) / np.hypot(x2 - x1, y2 - y1)
            assert cl12[pos] == pytest.approx(abs(np.dot(centre[n] - (x1, y1), normal)), rel=1e-12)
            assert hwva[pos] == pytest.approx(np.hypot(x2 - x1, y2 - y1), rel=1e-12)
            sides.append((cl12[pos], np.hypot(*(centre[m] - centre[n])), n, m))
    lengths = {(int(n), int(m)): length for length, _, n, m in sides}
    asymmetry = np.median([abs(a - lengths[(m, n)]) / (a + lengths[(m, n)])
                           for (n, m), a in lengths.items() if n < m])
    spacing = np.array([(length, gap) for length, gap, _, _ in sides])
    # Each side is SHORTER than the spacing between the two centres.
    assert (spacing[:, 0] < spacing[:, 1]).all()
    if mesh == "irregular":
        # A Voronoi face is the perpendicular bisector: half the spacing, and the SAME on
        # both sides (measured 3.5e-11), so this mesh cannot tell which side is which.
        assert np.median(spacing[:, 0] / spacing[:, 1]) == pytest.approx(0.5, abs=0.02)
        assert asymmetry < 1e-9
    else:
        # Measured 0.17: the two sides differ, so the hand check above reads the side.
        assert asymmetry > 0.1


def test_calculate_distance_is_not_the_connection_length_and_says_so(monkeypatch):
    """``calculate_distance`` is the pre-2025 rule -- a polygon CENTROID to the nearest point
    of the shared boundary -- kept public, so its docstrings must not read as MF6's
    per-side length (ledger 328). Their examples run, and state the difference."""

    import doctest

    from myflopy.modflow.mf6.grid import connectivity, geometry

    verts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float)
    grid = VoronoiGridPlus(verts=verts, iverts=[[0, 3, 4, 1], [1, 4, 5, 2]],
                           xcyc=np.array([[0.25, 0.5], [1.5, 0.5]]))
    _, _, cl12, _, _ = grid.get_disu_connectivity()
    assert grid.calculate_distance(grid.gdf_vorPolys, 0, 1) == 0.5   # from the centroid
    assert cl12[1] == 0.75                                           # MF6's, from (0.25, 0.5)
    for function in (geometry.calculate_distance, VoronoiGridPlus.calculate_distance):
        assert "NOT MODFLOW 6's connection length" in " ".join(function.__doc__.split()), function
    # doctest writes the interpreter-wide ``builtins._``: its display hook stores each
    # example's value there and every run ends with ``builtins._ = None``. Once any test
    # has imported GRASS that is GRASS's translation function, which its own
    # ``setup.finish()`` calls, so leaving it None broke every later contour
    # interpolation in the session (``'NoneType' object is not callable``). Restored at
    # teardown, or removed when it was absent.
    monkeypatch.setattr(builtins, "_", getattr(builtins, "_", None), raising=False)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    finder = doctest.DocTestFinder()
    for subject in (geometry.calculate_distance, VoronoiGridPlus.calculate_distance,
                    connectivity.build_disu_connectivity, connectivity.DisvFaces.length_and_width):
        for test in finder.find(subject, subject.__qualname__, globs={"DisvFaces": connectivity.DisvFaces}):
            runner.run(test, out=lambda _text: None)
    results = runner.summarize(verbose=False)
    assert results.attempted >= 12 and results.failed == 0, results
