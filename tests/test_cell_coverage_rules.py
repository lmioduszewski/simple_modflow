"""The cell-coverage rules have one home, and ``GeoPackageSource`` reads it.

Which cells a geometry reaches -- a polygon covers by AREA, a line on a face is
split evenly, a point goes to the lowest-numbered usable cell, a piece thinner
than ``_SLIVER_WIDTH`` is noise (ledgers 210-212) -- was decided inside
``GeoPackageSource``, where nothing else could reach it, while three other places
answering the same question each kept their own copy. The rules are
``grid/geometry.py`` now (import layer 1, so any of them can import it) and
``GeoPackageSource`` delegates.

Nothing about the answers changed. These tests pin that there is ONE copy, and
exercise it directly on the 4x4 grid of ``test_gpkg_on_cell_faces`` -- 25 ft
square cells numbered row-major from the bottom left, so cell 5 is
``box(25, 25, 50, 50)`` and its face-neighbours are 1, 4, 6 and 9.
"""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import LineString, MultiPoint, Point, box

import myflopy as mf
from myflopy import geopackage
from myflopy.geopackage import GeoPackageSource
from myflopy.modflow.mf6.grid import geometry
from myflopy.modflow.mf6.grid.geometry import cell_shares, domain_active, miss_reason
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from tests.test_gpkg_on_cell_faces import CELL_5, _cells, _context, _gpkg

FAR = box(1_300_000.0, 650_000.0, 1_300_100.0, 650_100.0)
BORDERING = box(-50.0, 0.0, 0.0, 100.0)


@pytest.fixture
def polys():
    """The 4x4 grid's cell frame, which is all ``cell_shares`` needs of a grid."""

    return _context().grid.gdf_vorPolys


# --- one copy, and GeoPackageSource reads it -------------------------------- #


def test_the_rules_are_one_object_the_geopackage_source_imported():
    """A ``from`` import binds the object, so identity is what proves there is one."""

    assert geopackage.cell_shares is geometry.cell_shares
    assert geopackage.domain_active is geometry.domain_active
    assert geopackage.miss_reason is geometry.miss_reason
    assert geopackage._extent is geometry._extent


def test_a_geopackage_source_resolves_its_cells_through_the_shared_rule(tmp_path, monkeypatch):
    """Patch the one copy and the answer moves: a polygon on cell 5 lands on cell 2."""

    monkeypatch.setattr(geopackage, "cell_shares", lambda cells, shape, *, usable=None: [(2, 1.0)])
    spec = mf.drn.gpkg(_gpkg(tmp_path, CELL_5), layer="bc", context=_context(), nper=1)
    assert _cells(spec) == [2]


def test_a_geopackage_source_decides_activity_through_the_shared_rule(tmp_path, monkeypatch):
    """With the shared test saying nothing is active, the source refuses as a full miss."""

    monkeypatch.setattr(geopackage, "domain_active", lambda domain, layer, cell: False)
    with pytest.raises(ValueError, match=r"'feature' \(on no active cell\)"):
        mf.drn.gpkg(_gpkg(tmp_path, CELL_5), layer="bc", context=_context(), nper=1)


def test_a_geopackage_source_words_a_miss_through_the_shared_classifier(tmp_path, monkeypatch):
    monkeypatch.setattr(geopackage, "miss_reason", lambda cells, shape, *, target: "elsewhere")
    with pytest.raises(ValueError, match=r"'feature' \(elsewhere\)"):
        mf.drn.gpkg(_gpkg(tmp_path, FAR), layer="bc", context=_context(), nper=1)


# --- cell_shares, called directly ------------------------------------------- #


def test_a_polygon_covers_the_cell_it_overlaps_and_not_its_face_neighbours(polys):
    assert cell_shares(polys, CELL_5) == [(5, 1.0)]


def test_a_piece_thinner_than_the_tolerance_is_noise_not_coverage(polys):
    """1e-6 ft over the face into cell 6 is a reprojection artefact, not an overlap."""

    shares = cell_shares(polys, box(25.0, 25.0, 50.0 + 1e-6, 50.0))
    assert [cell for cell, _ in shares] == [5]
    assert shares[0][1] == pytest.approx(1.0)


def test_a_real_overlap_of_the_neighbour_is_kept_with_its_share(polys):
    """0.015 ft in is past the tolerance: coverage, and shares of the whole feature."""

    shares = dict(cell_shares(polys, box(25.0, 25.0, 50.015, 50.0)))
    assert shares == pytest.approx({5: 625.0 / 625.375, 6: 0.375 / 625.375})


def test_a_thin_piece_that_is_a_real_share_of_its_cell_is_coverage():
    """The width rule needs its share guard: on 0.05-unit cells a clipped edge piece
    is 0.007 wide -- under the tolerance -- but 14% of its cell, and 36 of the 100
    cells the zone covers are reached by nothing else."""

    small = _context(n=10, length=0.5).grid.gdf_vorPolys
    shares = cell_shares(small, box(0.043, 0.043, 0.457, 0.457))
    assert len(shares) == 100
    assert sum(share for _, share in shares) == pytest.approx(1.0)


def test_a_feature_that_is_all_slivers_is_the_feature_and_is_kept(polys):
    """0.015 ft wide, centred on the face x=50: two halves of 0.0075, each under the
    tolerance and each 0.03% of its cell. Dropping every sliver would lose the whole
    strip, so they go only when something thicker survives them."""

    strip = box(50.0 - 0.0075, 25.0, 50.0 + 0.0075, 50.0)
    assert dict(cell_shares(polys, strip)) == pytest.approx({5: 0.5, 6: 0.5})


def test_the_shared_rule_sorts_whatever_order_a_caller_hands_back(polys):
    """``usable`` is documented lowest-first, but the tie rule does not rest on three
    callers remembering it: a reversed list still gives the vertex to the lowest cell
    and still returns the pairs in cell order."""

    def reversed_order(cells):
        return sorted(cells, reverse=True)

    assert cell_shares(polys, Point(25.0, 25.0), usable=reversed_order) == [(0, 1.0)]
    straddling = cell_shares(polys, box(25.0, 25.0, 50.015, 50.0), usable=reversed_order)
    assert [cell for cell, _ in straddling] == [5, 6]


def test_a_line_on_a_face_is_split_evenly_between_its_two_cells(polys):
    shares = dict(cell_shares(polys, LineString([(0.0, 25.0), (100.0, 25.0)])))
    assert shares == pytest.approx({cell: 0.125 for cell in range(8)})
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_line_along_the_grid_edge_keeps_its_whole_length(polys):
    shares = dict(cell_shares(polys, LineString([(0.0, 0.0), (100.0, 0.0)])))
    assert shares == pytest.approx({0: 0.25, 1: 0.25, 2: 0.25, 3: 0.25})


@pytest.mark.parametrize(
    ("point", "cell"),
    [(Point(12.5, 12.5), 0), (Point(25.0, 12.5), 0), (Point(75.0, 75.0), 10)],
    ids=["interior", "face", "vertex-of-four"],
)
def test_a_point_goes_whole_to_the_lowest_numbered_usable_cell(polys, point, cell):
    assert cell_shares(polys, point) == [(cell, 1.0)]


def test_a_multi_part_point_divides_the_feature_between_its_parts(polys):
    """A point has no extent to divide, so the FEATURE is divided instead.

    Every part used to get a share of 1.0, which contradicts the rule stated three
    paragraphs into `cell_shares`'s own docstring -- that fractions are of the WHOLE
    feature. Measured on the default path: a MultiPoint of 3 wells declared at -1000
    wrote **-3000** under `mf.Spread`, which exists to conserve a feature's total
    (ledger 287). The extent-based branch was always right, which is why only points
    were affected: a MultiLineString with one part off-grid already summed to 0.5.
    """

    three = MultiPoint([(12.5, 12.5), (37.5, 37.5), (62.5, 62.5)])
    shares = dict(cell_shares(polys, three))

    assert len(shares) == 3
    assert sum(shares.values()) == pytest.approx(1.0)
    assert all(share == pytest.approx(1 / 3) for share in shares.values())


def test_a_part_that_reaches_no_cell_leaves_the_feature_short(polys):
    """The second half of the same defect: a truncated feature has to LOOK truncated.

    With every part at 1.0 a MultiPoint whose third part was 50,000 ft away summed to
    2.0 -- indistinguishable from a complete two-part feature -- so nothing downstream
    could tell that a third of the declared pumping had gone. It never entered `missed`
    either, because that records a feature reaching ZERO cells.
    """

    truncated = MultiPoint([(12.5, 12.5), (37.5, 37.5), (50_000.0, 50_000.0)])
    shares = dict(cell_shares(polys, truncated))

    assert len(shares) == 2
    assert sum(shares.values()) == pytest.approx(2 / 3)


def test_two_parts_on_one_cell_owe_it_both_their_shares(polys):
    """Accumulated, not a set. Collapsing them would drop the second part's share and
    the feature would read as truncated when nothing is missing."""

    both_in_cell_zero = MultiPoint([(6.0, 6.0), (18.0, 18.0)])
    shares = dict(cell_shares(polys, both_in_cell_zero))

    assert shares == pytest.approx({0: 1.0})


def test_a_single_point_is_unchanged(polys):
    """The refusing direction: dividing by the part count must not move the ordinary
    case, which is one part and therefore one whole share."""

    assert cell_shares(polys, Point(12.5, 12.5)) == [(0, 1.0)]


def test_usable_chooses_which_cells_are_eligible(polys):
    """The caller's own set -- one layer's active cells, or the grid edge -- goes here."""

    without_the_corner = cell_shares(
        polys, Point(25.0, 25.0), usable=lambda cells: sorted(c for c in cells if c != 0)
    )
    assert without_the_corner == [(1, 1.0)]


def test_the_lowest_numbered_cell_wins_whatever_order_the_index_returns():
    """The default keeps every candidate SORTED, because the index answers in tree
    order: on this grid all nine interior vertices come back out of cell order."""

    slot_of = [5, 14, 0, 11, 7, 2, 13, 8, 1, 15, 4, 10, 3, 12, 9, 6]  # cell -> position
    d = 25.0
    verts = np.array([[i * d, j * d] for j in range(5) for i in range(5)], dtype=float)
    iverts, xcyc = [], []
    for slot in slot_of:
        i, j = slot % 4, slot // 4
        iverts.append([j * 5 + i, (j + 1) * 5 + i, (j + 1) * 5 + i + 1, j * 5 + i + 1])
        xcyc.append([(i + 0.5) * d, (j + 0.5) * d])
    scrambled = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(xcyc)).gdf_vorPolys
    cell_at = {slot: cell for cell, slot in enumerate(slot_of)}

    out_of_order = 0
    for i in range(1, 4):
        for j in range(1, 4):
            point = Point(i * d, j * d)
            hits = scrambled.sindex.query(point, predicate="intersects")
            raw = [int(cell) for cell in scrambled.index[hits]]
            out_of_order += raw != sorted(raw)
            lowest = min(cell_at[(j - dj) * 4 + (i - di)] for dj in (0, 1) for di in (0, 1))
            assert cell_shares(scrambled, point) == [(lowest, 1.0)]
    assert out_of_order, "the index answered in cell order, so this tests nothing"


@pytest.mark.parametrize("feature", [None, LineString()], ids=["none", "empty"])
def test_a_geometry_with_nothing_in_it_reaches_no_cell(polys, feature):
    assert cell_shares(polys, feature) == []


def test_a_row_with_no_geometry_never_walks_the_grid_edge(tmp_path, monkeypatch):
    """``get_grid_edge()`` walks every cell and reads the grid's own ``idomain_path``
    if it has one, so the source resolves it only once a feature reaches something --
    a null-geometry row is answered before that, as it was before the extraction."""

    context = _context()
    source = GeoPackageSource(_gpkg(tmp_path, CELL_5), context, 1, layer="bc")

    def walked(*args, **kwargs):
        raise AssertionError("the grid edge was resolved for a row with no geometry")

    monkeypatch.setattr(context.grid, "get_grid_edge", walked)
    assert source._cell_shares(None, layer=0, edges_only=True) == []


# --- domain_active ---------------------------------------------------------- #


@pytest.mark.parametrize(
    ("domain", "layer", "cell", "active"),
    [
        (None, 0, 0, True),
        (np.array([[1, 0]]), 0, 0, True),
        (np.array([[1, 0]]), 0, 1, False),
        (np.array([[-1, 1]]), 0, 0, False),
        (np.array([1, -1, 0, 1]), 0, 1, False),
        (np.array([1, -1, 0, 1]), 3, 3, True),
        (np.array([[1, 1], [0, 1]]), 1, 0, False),
    ],
    ids=["no domain", "active", "inactive", "passthrough", "1d passthrough",
         "1d applies to every layer", "second layer"],
)
def test_domain_active_reads_an_idomain_as_mf6_does(domain, layer, cell, active):
    """``> 0``, not truthiness: ``bool(-1)`` is True and a passthrough cell holds nothing."""

    assert domain_active(domain, layer, cell) is active


def test_a_layer_beyond_the_domain_raises_rather_than_taking_the_last_one():
    with pytest.raises(ValueError, match="Layer 3 is outside the model domain, which has 1 layers"):
        domain_active(np.ones((1, 16), dtype=int), 3, 0)


# --- miss_reason ------------------------------------------------------------ #


@pytest.mark.parametrize(
    ("feature", "kwargs", "reason"),
    [
        (None, {}, "no geometry"),
        (FAR, {}, "outside the grid"),
        (BORDERING, {}, "only borders cells -- a polygon must overlap a cell to reach it"),
        (CELL_5, {}, "on no active cell"),
        (CELL_5, {"target": "active grid-edge cell"}, "on no active grid-edge cell"),
    ],
    ids=["no geometry", "off the grid", "borders only", "on no active cell", "edges only"],
)
def test_miss_reason_words_each_kind_of_miss(polys, feature, kwargs, reason):
    assert miss_reason(polys, feature, **kwargs) == reason


def test_the_geopackage_refusal_still_carries_the_shared_wording(tmp_path):
    """The wording moved module; the message a user sees did not."""

    with pytest.raises(ValueError) as caught:
        mf.chd.gpkg(
            _gpkg(tmp_path, BORDERING), layer="bc", context=_context(), nper=1, edges_only=True
        )
    assert "only borders cells -- a polygon must overlap a cell to reach it" in str(caught.value)


def test_a_caller_with_no_geopackage_and_no_context_gets_the_same_rules(polys):
    """How the head-target, selector and grid-edge callers reach them: a cell frame,
    and their own eligibility closed over -- here, layer 0 of an idomain."""

    active = np.array([[0] + [1] * 15])
    reached = cell_shares(
        polys,
        Point(25.0, 25.0),
        usable=lambda cells: sorted(c for c in cells if domain_active(active, 0, c)),
    )
    assert reached == [(1, 1.0)]
    assert miss_reason(polys, FAR) == "outside the grid"
