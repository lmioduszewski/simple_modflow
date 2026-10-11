"""A polygon drawn on cell faces covers the cells inside it, and no others.

``GeoPackageSource._cell_shares`` measured each piece of a feature by the piece's
own extent: area if it had any, otherwise LENGTH. A polygon meets each
face-adjacent cell along the edge they share -- a line, with no area -- so every
neighbour of a polygon drawn on cell faces was admitted, with a pseudo-share of
edge / area (L / L**2). A field not wrapped in ``mf.Spread`` is written whole to
each admitted cell, so a DRN polygon that is exactly one cell became five records
at five times the conductance, and a K zone painted five cells. MF6 terminates
normally either way.

Exact coincidence is common, not contrived: square-cell DISV grids with round
spacing, polygons dissolved or exported from the grid itself, nested refinement.
Its near miss is commoner still -- the same polygon after a CRS round trip, or
read back from a DISV FloPy wrote at 9 significant digits -- so a piece narrower
than ``_SLIVER_WIDTH`` is dropped as well.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import shapely
from shapely.geometry import GeometryCollection, LineString, MultiPoint, MultiPolygon, Point, box

import myflopy as mf
from myflopy.geopackage import GeoPackageSource, Spread
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.specs import ModelContext

LIST_BCS = ("chd", "ghb", "drn", "riv", "wel", "rch", "evt")

# Every registry default column, so one file serves all seven helpers.
_VALUES = {
    "head": 10.0,
    "conductance": 100.0,
    "elevation": 10.0,
    "stage": 10.0,
    "rbot": 5.0,
    "rate": -500.0,
    "recharge": 0.01,
    "surface": 10.0,
    "depth": 2.0,
    "k": 7.0,
}

# Exactly cell 5 of the 4x4 grid below; its face-neighbours are 1, 4, 6 and 9.
CELL_5 = box(25.0, 25.0, 50.0, 50.0)


def _context(n: int = 4, length: float = 100.0, x0: float = 0.0, y0: float = 0.0):
    """A real ``VoronoiGridPlus`` of n x n square cells, row-major from the bottom left."""

    d = length / n
    verts = np.array(
        [[x0 + i * d, y0 + j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float
    )
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array(
        [[x0 + (i + 0.5) * d, y0 + (j + 0.5) * d] for j in range(n) for i in range(n)]
    )
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    return ModelContext(grid=grid, domain=np.ones((1, n * n), dtype=int))


def _gpkg(tmp_path, geometry):
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["feature"], "layer": [1], **{k: [v] for k, v in _VALUES.items()}},
        geometry=[geometry],
        crs="EPSG:2927",
    ).to_file(path, layer="bc", driver="GPKG")
    return path


def _records(spec):
    return spec.options["stress_period_data"][0]


def _cells(spec):
    return sorted(record[0][1] for record in _records(spec))


# --- exact coincidence ------------------------------------------------------ #


@pytest.mark.parametrize("package", LIST_BCS)
def test_a_polygon_that_is_one_cell_resolves_to_that_cell_in_every_list_bc(tmp_path, package):
    spec = getattr(mf, package).gpkg(
        _gpkg(tmp_path, CELL_5), layer="bc", context=_context(), nper=1
    )
    assert _cells(spec) == [5]


@pytest.mark.parametrize(
    "conductance",
    [
        "conductance",
        Spread("conductance"),
        Spread("conductance", mode="retained"),
        Spread("conductance", min_share=0.0),
    ],
    ids=["whole", "clip", "retained", "no-floor"],
)
def test_the_whole_value_lands_in_the_one_cell_however_it_is_split(tmp_path, conductance):
    """The pseudo-share was L / L**2, so whether min_share hid it depended on units."""

    records = _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, CELL_5), layer="bc", context=_context(), nper=1,
            conductance=conductance,
        )
    )
    assert [record[0][1] for record in records] == [5]
    assert records[0][2] == pytest.approx(100.0)


@pytest.mark.parametrize("n", [4, 8])
def test_refining_under_a_face_aligned_polygon_keeps_its_cell_count_and_total(tmp_path, n):
    """On fine faces but across coarse cells: 4 cells either way, never 12."""

    records = _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, box(37.5, 37.5, 62.5, 62.5)), layer="bc", context=_context(n),
            nper=1, conductance=Spread("conductance"),
        )
    )
    assert len(records) == 4
    assert sum(record[2] for record in records) == pytest.approx(100.0)


def test_retained_mode_keeps_the_total_inside_the_polygon(tmp_path):
    """Renormalizing the pseudo-shares conserved the total and put 16% of it outside."""

    records = _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, box(37.5, 37.5, 62.5, 62.5)), layer="bc", context=_context(8),
            nper=1, conductance=Spread("conductance", mode="retained"),
        )
    )
    assert {record[0][1]: record[2] for record in records} == pytest.approx(
        {27: 25.0, 28: 25.0, 35: 25.0, 36: 25.0}
    )


def test_a_k_zone_on_cell_faces_paints_only_its_own_cell(tmp_path):
    source = GeoPackageSource(_gpkg(tmp_path, CELL_5), _context(), 1, layer="bc")
    k = source.k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [5]


# --- near coincidence, and what must NOT be mistaken for it ----------------- #
#
# A perturbation is a WIDTH, not a share of the cell, so the rule is on a piece's
# mean width (2 * area / perimeter). A share rule caught a NAD83 round trip on 500 ft
# cells and missed it on 300 ft ones -- this file only tested 500 ft, and missed it.

# Non-round, and away from the origin, where the NAD83/HARN shift vanishes.
STATE_PLANE = dict(x0=1_300_000.123456, y0=650_000.654321)


def _slivered(moved, polys):
    return [c for c in range(len(polys)) if c != 5 and moved.intersection(polys.loc[c]).area > 0]


@pytest.mark.parametrize("cell_size", [500.0, 100.0, 25.0])
@pytest.mark.parametrize("via", ["EPSG:4326", "EPSG:32610", "EPSG:2285"])
def test_the_same_polygon_after_a_crs_round_trip_still_resolves_to_one_cell(
    tmp_path, via, cell_size
):
    context = _context(length=4 * cell_size, **STATE_PLANE)
    polys = context.grid.gdf_vorPolys.geometry
    moved = gpd.GeoSeries([polys.loc[5]], crs="EPSG:2927").to_crs(via).to_crs("EPSG:2927").iloc[0]
    assert _slivered(moved, polys), "the round trip must leave a sliver, or this tests nothing"

    spec = mf.drn.gpkg(_gpkg(tmp_path, moved), layer="bc", context=context, nper=1)
    assert _cells(spec) == [5]


@pytest.mark.parametrize("cell_size", [100.0, 25.0])
def test_a_polygon_read_back_from_a_written_disv_still_resolves_to_one_cell(tmp_path, cell_size):
    """FloPy writes DISV vertices at 9 significant digits -- 0.01 ft at State Plane --
    so a zone drawn from the WRITTEN cells is off the in-memory mesh by up to 0.005."""

    context = _context(length=4 * cell_size, **STATE_PLANE)
    polys = context.grid.gdf_vorPolys.geometry
    written = shapely.transform(
        polys.loc[5], lambda xy: np.array([[float(f"{v:.8E}") for v in point] for point in xy])
    )
    assert _slivered(written, polys), "rounding must leave a sliver, or this tests nothing"

    source = GeoPackageSource(_gpkg(tmp_path, written), context, 1, layer="bc")
    k = source.k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [5]


def test_a_feature_thinner_than_the_tolerance_is_kept_whole(tmp_path):
    """A pond 0.01 ft across is what was drawn, not a sliver of something else."""

    spec = mf.drn.gpkg(
        _gpkg(tmp_path, box(760.0, 760.0, 760.01, 760.01)), layer="bc",
        context=_context(length=2000.0), nper=1,
    )
    assert _cells(spec) == [5]


def test_a_small_part_of_a_multipolygon_is_kept_like_the_same_part_alone(tmp_path):
    """Dissolving a pond into its zone must not lose it: 0.16 ft2 of a 500 ft cell."""

    zone = MultiPolygon([box(0.0, 0.0, 1000.0, 1000.0), box(1760.0, 1760.0, 1760.4, 1760.4)])
    source = GeoPackageSource(_gpkg(tmp_path, zone), _context(length=2000.0), 1, layer="bc")
    k = source.k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [0, 1, 4, 5, 15]


def test_a_real_overlap_of_a_neighbour_still_reaches_it(tmp_path):
    """A strip 0.015 ft wide into cell 6 is past the tolerance: coverage, not noise."""

    spec = mf.drn.gpkg(
        _gpkg(tmp_path, box(25.0, 25.0, 50.015, 50.0)), layer="bc", context=_context(), nper=1
    )
    assert _cells(spec) == [5, 6]


def test_a_collection_mixing_polygons_and_lines_is_measured_by_its_polygons(tmp_path):
    """It was refused (the first ledger 210); a clip leaves such collections by default.
    The line adds no area, so cells 2 and 3, which only it crosses, are not reached."""

    mixed = GeometryCollection(
        [box(30.0, 30.0, 45.0, 45.0), LineString([(60.0, 10.0), (90.0, 10.0)])]
    )
    spec = mf.chd.gpkg(_gpkg(tmp_path, mixed), layer="bc", context=_context(), nper=1)
    assert _cells(spec) == [5]


# --- lines and points on faces: the tie rule -------------------------------- #
#
# A line lying ON a face is in both cells' pieces, so its shares summed to 2: a
# Spread conductance doubled. It is split evenly between the two cells now. A
# point on a face or vertex touches two to four cells and each took the WHOLE
# value; the lowest-numbered active one takes it now. Decided 2026-09-11.


def _line_records(tmp_path, geometry, context=None, **fields):
    return _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, geometry), layer="bc", context=context or _context(), nper=1,
            conductance=Spread("conductance"), **fields,
        )
    )


def _by_cell(records, column=2):
    return {record[0][1]: record[column] for record in records}


def test_a_line_on_a_face_is_split_evenly_between_its_two_cells(tmp_path):
    records = _line_records(tmp_path, LineString([(0.0, 25.0), (100.0, 25.0)]))
    assert _by_cell(records) == pytest.approx({cell: 12.5 for cell in range(8)})


def test_a_line_on_the_grid_edge_has_one_cell_per_stretch_and_keeps_its_length(tmp_path):
    records = _line_records(tmp_path, LineString([(0.0, 0.0), (100.0, 0.0)]))
    assert _by_cell(records) == pytest.approx({0: 25.0, 1: 25.0, 2: 25.0, 3: 25.0})


def test_a_line_on_a_face_whose_far_side_is_inactive_stays_whole_on_the_near_side(tmp_path):
    near_row_only = ModelContext(
        grid=_context().grid, domain=np.array([[1, 1, 1, 1] + [0] * 12])
    )
    records = _line_records(tmp_path, LineString([(0.0, 25.0), (100.0, 25.0)]), near_row_only)
    assert _by_cell(records) == pytest.approx({0: 25.0, 1: 25.0, 2: 25.0, 3: 25.0})


def test_a_line_partly_on_a_face_conserves_its_total(tmp_path):
    records = _line_records(tmp_path, LineString([(0.0, 25.0), (50.0, 25.0), (100.0, 40.0)]))
    assert sum(record[2] for record in records) == pytest.approx(100.0)
    assert sorted(_by_cell(records)) == [0, 1, 4, 5, 6, 7]


def test_a_diagonal_through_vertices_gives_nothing_to_the_cells_it_only_touches(tmp_path):
    records = _line_records(tmp_path, LineString([(0.0, 0.0), (100.0, 100.0)]))
    assert _by_cell(records) == pytest.approx({0: 25.0, 5: 25.0, 10: 25.0, 15: 25.0})


def test_a_face_shared_with_two_smaller_cells_is_still_split_by_length(tmp_path):
    """A hanging node: one coarse cell against two fine ones across the same face."""

    verts = np.array(
        [[0, 0], [0, 2], [2, 2], [2, 1], [2, 0], [3, 0], [3, 1], [3, 2]], dtype=float
    )
    iverts = [[0, 1, 2, 3, 4], [4, 3, 6, 5], [3, 2, 7, 6]]
    xcyc = np.array([[1.0, 1.0], [2.5, 0.5], [2.5, 1.5]])
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    context = ModelContext(grid=grid, domain=np.ones((1, 3), dtype=int))
    records = _line_records(tmp_path, LineString([(2.0, 0.0), (2.0, 2.0)]), context)
    assert _by_cell(records) == pytest.approx({0: 50.0, 1: 25.0, 2: 25.0})


def test_an_elevation_on_a_face_reaches_both_sides_unchanged(tmp_path):
    """Only a Spread field is divided; an elevation is the same on either side."""

    records = _line_records(tmp_path, LineString([(0.0, 25.0), (100.0, 25.0)]))
    assert {record[1] for record in records} == {10.0}


@pytest.mark.parametrize(
    ("point", "cell"),
    [(Point(12.5, 12.5), 0), (Point(25.0, 12.5), 0), (Point(25.0, 25.0), 0), (Point(75.0, 75.0), 10)],
    ids=["interior", "face", "vertex", "vertex-of-four"],
)
def test_a_point_goes_whole_to_one_cell_the_lowest_numbered(tmp_path, point, cell):
    records = _records(
        mf.wel.gpkg(_gpkg(tmp_path, point), layer="bc", context=_context(), nper=1)
    )
    assert _by_cell(records, column=1) == pytest.approx({cell: -500.0})


def test_a_point_skips_an_inactive_lowest_cell(tmp_path):
    context = ModelContext(grid=_context().grid, domain=np.array([[0] + [1] * 15]))
    records = _records(
        mf.wel.gpkg(_gpkg(tmp_path, Point(25.0, 25.0)), layer="bc", context=context, nper=1)
    )
    assert _by_cell(records, column=1) == pytest.approx({1: -500.0})


def test_each_point_of_a_multipoint_gets_its_own_cell(tmp_path):
    records = _records(
        mf.wel.gpkg(
            _gpkg(tmp_path, MultiPoint([(25.0, 25.0), (75.0, 75.0)])), layer="bc",
            context=_context(), nper=1,
        )
    )
    assert _by_cell(records, column=1) == pytest.approx({0: -500.0, 10: -500.0})


def test_edges_only_picks_the_lowest_numbered_edge_cell(tmp_path):
    """On the face between cells 6 (interior) and 7 (edge): 6 is lower, 7 is the edge."""

    records = _records(
        mf.chd.gpkg(
            _gpkg(tmp_path, Point(75.0, 37.5)), layer="bc", context=_context(), nper=1,
            edges_only=True,
        )
    )
    assert [record[0][1] for record in records] == [7]


def test_the_lowest_numbered_cell_wins_whatever_order_the_index_returns(tmp_path):
    """The spatial index hands hits back in tree order, not cell order (measured:
    481 of 685 queries out of order on a scrambled 20x20), so the choice is sorted."""

    slot_of = [5, 14, 0, 11, 7, 2, 13, 8, 1, 15, 4, 10, 3, 12, 9, 6]  # cell -> position
    d = 25.0
    verts = np.array([[i * d, j * d] for j in range(5) for i in range(5)], dtype=float)
    iverts, xcyc = [], []
    for slot in slot_of:
        i, j = slot % 4, slot // 4
        iverts.append([j * 5 + i, (j + 1) * 5 + i, (j + 1) * 5 + i + 1, j * 5 + i + 1])
        xcyc.append([(i + 0.5) * d, (j + 0.5) * d])
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(xcyc))
    cell_at = {slot: cell for cell, slot in enumerate(slot_of)}
    vertices = [(i, j) for j in range(1, 4) for i in range(1, 4)]
    expected = {
        f"v{i}{j}": min(cell_at[(j - dj) * 4 + (i - di)] for dj in (0, 1) for di in (0, 1))
        for i, j in vertices
    }
    path = tmp_path / "wells.gpkg"
    gpd.GeoDataFrame(
        {"name": list(expected), "layer": [1] * 9, "rate": [-500.0] * 9},
        geometry=[Point(i * d, j * d) for i, j in vertices],
        crs="EPSG:2927",
    ).to_file(path, layer="bc", driver="GPKG")
    context = ModelContext(grid=grid, domain=np.ones((1, 16), dtype=int))
    records = _records(mf.wel.gpkg(path, layer="bc", context=context, nper=1))
    assert {record[-1]: record[0][1] for record in records} == expected


# --- a second review of the width rule (ledger 210, revised) ---------------- #


def test_a_piece_that_also_touches_its_cell_along_a_face_is_measured_by_its_area(tmp_path):
    """Cell 1 exactly, plus a 1.2 ft square inside cell 4: the piece in cell 4 is a
    GeometryCollection holding the shared face as a line. Counted in the perimeter,
    that line made the square a 0.00998 'sliver'."""

    context = _context(n=3, length=900.0)
    zone = MultiPolygon([box(300.0, 0.0, 600.0, 300.0), box(449.4, 449.4, 450.6, 450.6)])
    source = GeoPackageSource(_gpkg(tmp_path, zone), context, 1, layer="bc")
    assert np.flatnonzero(source.k_array(value="k", nlay=1, defaults=1.0)[0] != 1.0).tolist() == [1, 4]


def test_a_lab_scale_grid_keeps_its_coverage(tmp_path):
    """0.01-unit cells: every piece is narrower than the tolerance, and all were lost."""

    context = _context(n=10, length=0.1)
    spec = mf.drn.gpkg(_gpkg(tmp_path, box(0.0, 0.0, 0.1, 0.1)), layer="bc", context=context, nper=1)
    assert len(_cells(spec)) == 100


def test_a_thin_piece_that_is_a_real_share_of_its_cell_is_kept(tmp_path):
    """0.05-unit cells: whole cells are 0.025 wide, and the cut edge pieces 0.007 wide
    -- under the tolerance, but 14% of their cell. That is coverage, not noise."""

    context = _context(n=10, length=0.5)
    zone = box(0.043, 0.043, 0.457, 0.457)
    spec = mf.drn.gpkg(_gpkg(tmp_path, zone), layer="bc", context=context, nper=1)
    assert len(_cells(spec)) == 100


def test_a_strip_centred_on_a_face_reaches_both_rows(tmp_path):
    """0.015 wide on y=25: two halves of 0.0075, both below the tolerance -- all slivers."""

    context = _context()
    spec = mf.drn.gpkg(
        _gpkg(tmp_path, box(0.0, 24.9925, 100.0, 25.0075)), layer="bc", context=context, nper=1
    )
    assert _cells(spec) == list(range(8))


def test_a_collection_with_an_empty_member_is_measured_by_the_rest(tmp_path):
    context = _context()
    spec = mf.drn.gpkg(
        _gpkg(tmp_path, GeometryCollection([CELL_5, LineString()])), layer="bc",
        context=context, nper=1,
    )
    assert _cells(spec) == [5]


def test_a_clipped_zone_keeps_its_polygon_and_drops_its_edge_artifacts(tmp_path):
    """What gpd.clip returns by default where a zone runs along the clip edge: its
    polygon plus zero-area lines. The lines lie on the face y=50 and add no area."""

    context = _context()
    clipped = GeometryCollection([CELL_5, LineString([(50.0, 50.0), (100.0, 50.0)])])
    source = GeoPackageSource(_gpkg(tmp_path, clipped), context, 1, layer="bc")
    assert np.flatnonzero(source.k_array(value="k", nlay=1, defaults=1.0)[0] != 1.0).tolist() == [5]
