"""A USG export chain runs between cells that share a FACE, not a corner.

``usg/export._adjacency`` asked the spatial index for ``predicate="touches"``,
which is true of a shared edge and equally true of a single shared vertex. On a
structured or quadtree mesh -- both of which MODFLOW-USG writes, and a ``.gsf``
describes -- that made every diagonal a chain and every 2x2 block a clique, and
the chain is what decides the exported ``length_ft`` and every ``*_per_ft`` twin:

* the diagonal staircase ``[0, 1, 5, 6, 10, 11]`` of 25 ft cells -- what a diagonal
  drain on a structured grid actually is, every step a real shared face -- measured
  145.71 ft where it occupies 125.0, three of its six segments 30.18 ft long inside
  a 25 ft cell, because the walk took corner shortcuts;
* the pure diagonal ``[0, 5, 10, 15]``, which shares no face at all, came back as
  ONE chain measuring 106.07 ft against the 75.0 ft four such cells can span;
* the 2x2 block ``[0, 1, 4, 5]`` came back as a clique, so the walk stepped 1 -> 4
  straight across the shared corner: 85.36 ft total and ``conductance_per_ft``
  3.31 on two of the four cells, where the ring gives 75.0 ft and 4.00;
* a coarse quadtree cell claimed the fine cell diagonally off its corner as a
  neighbour, inventing a chain step where MODFLOW 6 has no connection at all.

The rule is now ``grid.shared_face``, the one the HFB views resolve a barrier's
own face against, and pair-for-pair the older ``grid.adjacent_cells_idx`` on every
mesh measured. What it costs is that a cell alone on its chain
-- nothing shares a face with it, or the open walk could not reach it through a
branch -- has a degenerate segment and a NaN twin, so the export reports those.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

gpd = pytest.importorskip("geopandas")

from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from myflopy.modflow.usg._io import NameFile  # noqa: E402
from myflopy.modflow.usg.export import (  # noqa: E402
    ExportManifest,
    _adjacency,
    _bc_frames,
    _chains,
    _export_boundaries,
)
from myflopy.modflow.usg.model import UsgModel  # noqa: E402
from myflopy.modflow.usg.packages import BasData, DisuData, LpfData, PeriodRecords  # noqa: E402

LOGGER = "myflopy.modflow.usg.export"
CRS = "EPSG:2927"
WIDTH = 25.0  #: cell side, so a run of n cells spans (n - 1) * WIDTH


def square_grid(n: int = 4, width: float = WIDTH) -> VoronoiGridPlus:
    """``n`` x ``n`` square cells, row-major from the bottom left.

    Cell 0 is 0..25 in both axes, so 0 and 5 meet only at (25, 25) and 0 and 1
    share the face at x = 25.
    """

    verts = np.array([[i * width, j * width] for j in range(n + 1) for i in range(n + 1)], float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    centres = np.array(
        [[(i + 0.5) * width, (j + 0.5) * width] for j in range(n) for i in range(n)], float
    )
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=centres, crs=CRS)


def quadtree_grid() -> VoronoiGridPlus:
    """One coarse 20x20 cell and three 10x10 fine cells, as a refinement step makes.

    Cell 0 spans 0..20 in both axes. Cells 1 and 2 stack against its right face,
    each sharing half of it. Cell 3 sits above cell 2 and meets cell 0 at exactly
    one point, (20, 20) -- the case a corner test gets wrong.
    """

    rings = [
        [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)],
        [(20.0, 0.0), (30.0, 0.0), (30.0, 10.0), (20.0, 10.0)],
        [(20.0, 10.0), (30.0, 10.0), (30.0, 20.0), (20.0, 20.0)],
        [(20.0, 20.0), (30.0, 20.0), (30.0, 30.0), (20.0, 30.0)],
    ]
    verts: list[list[float]] = []
    seen: dict[tuple[float, float], int] = {}
    iverts = []
    for ring in rings:
        indices = []
        for point in ring:
            if point not in seen:
                seen[point] = len(verts)
                verts.append(list(point))
            indices.append(seen[point])
        iverts.append(indices)
    centres = np.array(
        [[np.mean([p[0] for p in r]), np.mean([p[1] for p in r])] for r in rings], float
    )
    return VoronoiGridPlus(verts=np.array(verts, float), iverts=iverts, xcyc=centres, crs=CRS)


def usg_model(grid: VoronoiGridPlus, cells, *, conductance: float = 100.0,
              ftype: str = "DRN") -> UsgModel:
    """A one-layer USG model whose only package is ``ftype`` on ``cells``.

    Everything outside DISU/BAS/LPF's required fields is the smallest value that
    parses, because what is under test is which cells the export chains together.
    CHD is offered because it is the package with no extensive column at all --
    nothing to carry a ``*_per_ft`` twin -- and ledger 169 measures it as the most
    branched of the four.
    """

    ncpl = len(grid.gdf_vorPolys)
    ones = np.ones(ncpl)
    disu = DisuData(
        nodes=ncpl, nlay=1, njag=0, ivsd=-1, nper=1, itmuni=4, lenuni=1,
        nodelay=np.array([ncpl]), top=ones * 10.0, bot=np.zeros(ncpl),
        area=grid.gdf_vorPolys.geometry.area.to_numpy(), iac=np.ones(ncpl, int),
        ja=np.arange(1, ncpl + 1), perlen=np.array([1.0]), nstp=np.array([1]),
        tsmult=np.array([1.0]), steady=np.array([False]),
    )
    bas = BasData(ibound=np.ones(ncpl, int), hnoflo=-999.0, strt=ones * 8.0)
    lpf = LpfData(
        ilpfcb=0, hdry=-1e30, options=(), laytyp=np.zeros(1, int), layavg=np.zeros(1, int),
        chani=np.ones(1), layvka=np.zeros(1, int), laywet=np.zeros(1, int),
        hk=ones * 10.0, vka=ones, ss=ones * 1e-5, sy=ones * 0.2,
    )
    second = 7.5 if ftype == "CHD" else conductance
    block = np.array([[c + 1, 7.5, second] for c in cells], float)
    return UsgModel(
        name_file=NameFile(path=Path("tiny.nam"), workspace=Path(".")),
        disu=disu, bas=bas, lpf=lpf, grid=grid,
        boundaries={ftype: PeriodRecords(ftype=ftype, periods={0: block})},
    )


def neighbours(cells, grid) -> dict[int, list[int]]:
    """``_adjacency`` as plain ints, so a failure reads as cell numbers."""

    return {int(k): sorted(int(c) for c in v) for k, v in _adjacency(cells, grid).items()}


def walks(cells, grid) -> list[list[int]]:
    return [[int(c) for c in walk] for walk in _chains(cells, grid)]


def every_step_is_a_face(chain, grid) -> bool:
    pairs = list(zip(chain, chain[1:], strict=False))
    return all(grid.shared_face(a, b) is not None for a, b in pairs)


# -- adjacency -------------------------------------------------------------


def test_a_diagonal_is_four_cells_and_not_one_chain():
    """[0, 5, 10, 15] meet only at vertices, so nothing is adjacent to anything."""

    grid = square_grid()
    assert neighbours([0, 5, 10, 15], grid) == {0: [], 5: [], 10: [], 15: []}
    assert sorted(walks([0, 5, 10, 15], grid)) == [[0], [5], [10], [15]]


def test_a_2x2_block_is_a_ring_and_not_a_clique():
    """Every cell of the block has two face neighbours; the diagonals are not pairs."""

    grid = square_grid()
    assert neighbours([0, 1, 4, 5], grid) == {0: [1, 4], 1: [0, 5], 4: [0, 5], 5: [1, 4]}


def test_a_block_walk_never_steps_across_a_corner():
    """The walk's order is not pinned, but each step must be a real connection."""

    grid = square_grid()
    chains = walks([0, 1, 4, 5], grid)
    assert len(chains) == 1 and sorted(chains[0]) == [0, 1, 4, 5]
    assert every_step_is_a_face(chains[0], grid), chains


def test_a_straight_run_is_unchanged():
    """The case the export was designed for: one chain, in order."""

    grid = square_grid()
    assert neighbours([0, 1, 2, 3], grid) == {0: [1], 1: [0, 2], 2: [1, 3], 3: [2]}
    assert walks([0, 1, 2, 3], grid) == [[0, 1, 2, 3]]
    assert walks([0, 4, 8, 12], grid) == [[0, 4, 8, 12]]


def test_a_staircase_stays_one_chain_and_takes_no_corner_shortcut():
    """A diagonal drain on a structured grid: every step face-sharing, in order.

    The walk used to be [0, 1, 5, 10, 11, 6] -- two corner jumps -- because the
    corner links gave it extra ways to continue.
    """

    grid = square_grid()
    chains = walks([0, 1, 5, 6, 10, 11], grid)
    assert len(chains) == 1
    assert chains[0] in ([0, 1, 5, 6, 10, 11], [11, 10, 6, 5, 1, 0])
    assert every_step_is_a_face(chains[0], grid), chains


def test_an_l_bend_stays_one_chain():
    """Two face-sharing legs meeting at a cell: still one feature, still ordered."""

    grid = square_grid()
    chains = walks([0, 1, 2, 6, 10], grid)
    assert len(chains) == 1
    assert chains[0] in ([0, 1, 2, 6, 10], [10, 6, 2, 1, 0])


def test_a_branching_feature_keeps_every_cell_on_some_run():
    """A cross has no single open walk, so the leftovers recurse into their own runs.

    Cell 5 is the branch point of four arms. Nothing may be dropped, and no run may
    step across a corner -- which is also what exercises the recursive call.
    """

    grid = square_grid()
    cells = [1, 4, 5, 6, 7, 9, 13]
    chains = walks(cells, grid)
    assert sorted(c for chain in chains for c in chain) == sorted(cells)
    assert len(chains) == 3  # one open walk through the branch point, plus two arms
    for chain in chains:
        assert every_step_is_a_face(chain, grid), chain


def test_a_single_cell_branch_leftover_is_kept_as_its_own_run():
    """The leftover of a branch can be ONE cell, and one cell is not droppable.

    A T: three arms on cell 5. The open walk takes two of them -- [1, 5, 4] -- and
    leaves exactly one cell over, which is the ``else [leftover]`` arm of
    ``_chains``. A two-cell leftover recurses instead and never reaches it, so
    without a T in the suite dropping a lone leftover -- and with it a BC record,
    silently -- passes every other test here.
    """

    grid = square_grid()
    cells = [1, 4, 5, 6]
    chains = walks(cells, grid)
    assert sorted(c for chain in chains for c in chain) == cells
    assert sorted(len(chain) for chain in chains) == [1, 3]


def test_a_plus_cross_leaves_two_arms_over_and_loses_none_of_them():
    """Four single-cell arms on one branch point: three cells walk, two stand alone.

    Both loners share a face with cell 5. They are branch arms, not corner
    contacts, which is what the report has to be careful about.
    """

    grid = square_grid()
    cells = [1, 4, 5, 6, 9]
    chains = walks(cells, grid)
    assert sorted(c for chain in chains for c in chain) == cells
    assert sorted(len(chain) for chain in chains) == [1, 1, 3]
    for chain in chains:
        if len(chain) == 1:
            assert grid.shared_face(chain[0], 5) is not None


def test_a_coarse_cell_meets_two_fine_ones_and_not_the_third():
    """A quadtree step: both halves of the shared face count, the corner does not."""

    grid = quadtree_grid()
    assert neighbours([0, 1, 2, 3], grid) == {0: [1, 2], 1: [0, 2], 2: [0, 1, 3], 3: [2]}
    chains = walks([0, 1, 2, 3], grid)
    assert len(chains) == 1 and sorted(chains[0]) == [0, 1, 2, 3]
    assert every_step_is_a_face(chains[0], grid), chains


def test_over_the_whole_grid_adjacency_equals_the_grids_own_adjacency():
    """Restricted to nothing, ``_adjacency`` IS ``grid.adjacent_cells_idx``."""

    grid = square_grid()
    cells = list(range(len(grid.gdf_vorPolys)))
    mine = neighbours(cells, grid)
    theirs = {c: sorted(int(n) for n in grid.adjacent_cells_idx[c]) for c in cells}
    assert mine == theirs


def test_the_face_test_is_the_grids_method_and_not_a_local_copy(monkeypatch):
    """``grid.shared_face`` is consulted for every candidate pair, not re-derived.

    Pinning the ANSWER is not enough: a hand-rolled
    ``isinstance(a.intersection(b), LineString)`` agrees with the grid on these
    meshes and disagrees on the ones ``shared_face`` was written for -- a face
    interrupted into a MultiLineString (measured: an 11 ft contact the hand-rolled
    test calls no contact at all), its GEOS fallback, ``FACE_TOLERANCE``. So assert
    the delegation itself, not just the answer.
    """

    grid = square_grid()
    real, asked = grid.shared_face, []

    def recording(cell_a, cell_b):
        asked.append((int(cell_a), int(cell_b)))
        return real(cell_a, cell_b)

    monkeypatch.setattr(grid, "shared_face", recording)
    assert neighbours([0, 1, 4, 5], grid) == {0: [1, 4], 1: [0, 5], 4: [0, 5], 5: [1, 4]}
    assert (0, 1) in asked and (0, 5) in asked  # the face pair and the corner pair


def test_two_separate_runs_stay_separate():
    grid = square_grid()
    assert sorted(walks([0, 1, 14, 15], grid)) == [[0, 1], [14, 15]]


# -- the exported columns --------------------------------------------------


def exported_lines(grid, cells, ftype="DRN"):
    _, lines = _bc_frames(usg_model(grid, cells, ftype=ftype), ftype, CRS)
    return lines.set_index("old_cell")


def test_a_straight_run_exports_the_length_its_cells_span():
    """Four 25 ft cells span 75 ft, and the conductance rebuilds from the twin."""

    lines = exported_lines(square_grid(), [0, 1, 2, 3])
    assert [lines.loc[c, "length_ft"] for c in (0, 1, 2, 3)] == pytest.approx([12.5, 25.0, 25.0, 12.5])
    assert float(lines["length_ft"].sum()) == pytest.approx(3 * WIDTH)
    assert [lines.loc[c, "conductance_per_ft"] for c in (0, 1, 2, 3)] == pytest.approx(
        [8.0, 4.0, 4.0, 8.0]
    )
    rebuilt = float((lines["conductance_per_ft"] * lines["length_ft"]).sum())
    assert rebuilt == pytest.approx(float(lines["conductance"].sum()))


def test_a_staircase_exports_the_length_it_occupies_and_a_flat_per_foot_twin():
    """145.71 ft with three 30.18 ft segments becomes 125.0 ft with none.

    This is the shape the fix is FOR: every cell is face-connected, so nothing
    goes NaN and the per-foot recipe is exact before and after -- it was the
    lengths, and therefore where a new mesh puts the conductance, that were wrong.
    """

    lines = exported_lines(square_grid(), [0, 1, 5, 6, 10, 11])
    assert float(lines["length_ft"].sum()) == pytest.approx(5 * WIDTH)
    assert sorted(lines["length_ft"]) == pytest.approx([12.5, 12.5, 25.0, 25.0, 25.0, 25.0])
    assert sorted(lines["conductance_per_ft"]) == pytest.approx([4.0, 4.0, 4.0, 4.0, 8.0, 8.0])
    assert not lines["conductance_per_ft"].isna().any()
    rebuilt = float((lines["conductance_per_ft"] * lines["length_ft"]).sum())
    assert rebuilt == pytest.approx(600.0)


def test_a_2x2_block_no_longer_exports_a_corner_crossing_length():
    """The clique walk measured 85.36 ft and 3.31 per ft; the ring measures 75 and 4."""

    lines = exported_lines(square_grid(), [0, 1, 4, 5])
    assert float(lines["length_ft"].sum()) == pytest.approx(3 * WIDTH)
    assert sorted(lines["length_ft"]) == pytest.approx([12.5, 12.5, 25.0, 25.0])
    assert sorted(lines["conductance_per_ft"]) == pytest.approx([4.0, 4.0, 8.0, 8.0])
    assert float((lines["conductance_per_ft"] * lines["length_ft"]).sum()) == pytest.approx(400.0)


def test_a_diagonal_feature_exports_no_length_rather_than_an_invented_one():
    """Each cell stands alone, so its segment is degenerate and its twin is NaN.

    106.07 ft of chain is the number this replaces. A NaN is not a per-foot
    rebuild -- the point layer is, which is why the manifest says so.
    """

    lines = exported_lines(square_grid(), [0, 5, 10, 15])
    assert list(lines["length_ft"]) == [0.0, 0.0, 0.0, 0.0]
    assert lines["conductance_per_ft"].isna().all()
    assert sorted(lines["chain"]) == [0, 1, 2, 3]


def test_a_mixed_feature_keeps_the_chained_cells_exact():
    """A three-cell run plus one cell off on its own: only the loner degenerates."""

    lines = exported_lines(square_grid(), [0, 1, 2, 15])
    assert [lines.loc[c, "length_ft"] for c in (0, 1, 2, 15)] == pytest.approx([12.5, 25.0, 12.5, 0.0])
    assert [lines.loc[c, "conductance_per_ft"] for c in (0, 1, 2)] == pytest.approx([8.0, 4.0, 8.0])
    assert np.isnan(lines.loc[15, "conductance_per_ft"])


def test_segments_still_tile_each_chain():
    """The invariant ledger 169 rests on, on a mesh with both kinds of contact."""

    lines = exported_lines(square_grid(), [0, 1, 2, 6, 10]).reset_index()
    for _, chain in lines.groupby(["layer", "chain"]):
        assert float(chain["length_ft"].sum()) == pytest.approx(
            chain.geometry.union_all().length, abs=1e-9
        )


# -- what the caller is told -----------------------------------------------


def written(grid, cells, tmp_path, ftype="DRN"):
    manifest = ExportManifest(directory=tmp_path, crs=CRS, nper=1)
    _export_boundaries(usg_model(grid, cells, ftype=ftype), tmp_path, CRS, manifest)
    return manifest


def test_corner_only_cells_are_counted_in_the_manifest_and_warned_about(tmp_path, caplog):
    grid = square_grid()
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        manifest = written(grid, [0, 5, 10, 15], tmp_path)
    entry = manifest.layers["boundaries.gpkg:drn_lines"]
    assert entry["cells_alone_on_a_chain"] == 4
    assert any("alone on their chain" in note for note in manifest.notes)
    assert "4 of 4" in caplog.text
    assert "conductance_per_ft is NaN" in caplog.text


def test_a_mixed_feature_reports_only_the_cell_that_is_alone(tmp_path, caplog):
    """One of four, not four of four -- and the two counts are not interchangeable."""

    grid = square_grid()
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        manifest = written(grid, [0, 1, 2, 15], tmp_path)
    assert manifest.layers["boundaries.gpkg:drn_lines"]["cells_alone_on_a_chain"] == 1
    assert "1 of 4" in caplog.text
    assert "1 of 4" in manifest.notes[0]


def test_a_branch_arm_the_walk_could_not_reach_is_reported_honestly(tmp_path, caplog):
    """Both loners of a plus-cross DO share a face; the report must not deny it.

    The count is of one-cell CHAINS, so it is right here -- but the sentence has to
    offer both causes, because "shares no face with the rest of the feature" is
    false for a branch arm, and CHD (the most branched package ledger 169 measured)
    is the likeliest to hit it.
    """

    grid = square_grid()
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        manifest = written(grid, [1, 4, 5, 6, 9], tmp_path)
    assert manifest.layers["boundaries.gpkg:drn_lines"]["cells_alone_on_a_chain"] == 2
    assert grid.shared_face(9, 5) is not None and grid.shared_face(6, 5) is not None
    note = manifest.notes[0]
    assert "2 of 5" in note and "branch" in note
    assert "branch" in caplog.text


def test_a_chd_report_promises_no_per_foot_column_because_it_has_none(tmp_path, caplog):
    """CHD's values are heads -- nothing extensive, so no twin to go NaN."""

    grid = square_grid()
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        manifest = written(grid, [0, 5, 10, 15], tmp_path, ftype="CHD")
    entry = manifest.layers["boundaries.gpkg:chd_lines"]
    assert entry["cells_alone_on_a_chain"] == 4
    assert not any(c.endswith("_per_ft") for c in entry["columns"])
    note = manifest.notes[0]
    assert "_per_ft" not in note and "NaN" not in note
    assert "_per_ft" not in caplog.text and "NaN" not in caplog.text


def test_a_fully_chained_feature_reports_no_cell_alone(tmp_path, caplog):
    grid = square_grid()
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        manifest = written(grid, [0, 1, 2, 3], tmp_path)
    assert manifest.layers["boundaries.gpkg:drn_lines"]["cells_alone_on_a_chain"] == 0
    assert not manifest.notes
    assert caplog.text == ""
