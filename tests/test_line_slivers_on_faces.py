"""A line overshooting a face by floating noise writes one record, not two.

``cell_shares`` has had a noise floor on its AREAL branch since ledger 210: a piece
narrower than ``_SLIVER_WIDTH`` and under ``_SLIVER_SHARE`` of its cell is a
reprojection artefact, not coverage. Its LINE branch never got one, so a centerline
drawn to end ON a face and overshooting it by 1e-6 ft gave the cell beyond a piece
of 1e-6 ft -- and every field not wrapped in ``mf.Spread`` is written whole, so that
was a second DRN record at the same elevation and the same conductance. A
2927 -> 4326 -> 2927 round trip moves an end by ~3e-7 ft and does it.

``SFRBuilder`` has applied the floor to a reach since ledger 214, so the two
disagreed about one geometry: SFR kept one reach, the ``.gpkg`` helpers wrote two
records. The floor is the shared rule's now, worded for a length -- so it reaches
every consumer of that rule, the cell selector included (ledger 235's list of
callers), not only a ``.gpkg`` source.

Ledger 212 is unchanged and pinned here: a stretch of line lying ALONG a face is in
both cells and is split evenly between them.
"""

from __future__ import annotations

import geopandas as gpd
import pytest
from shapely.geometry import LineString, box

import myflopy as mf
from myflopy.geopackage import Spread
from myflopy.modflow.mf6.grid.geometry import cell_shares
from myflopy.modflow.mf6.grid.selection import get_vor_cells_as_series
from tests.test_gpkg_on_cell_faces import LIST_BCS, _context, _gpkg, _records

# The 4x4 grid of 25 ft cells: cell 5 is box(25, 25, 50, 50) and the face at x=25
# is the one it shares with cell 4.
INSIDE_CELL_5 = (37.5, 37.5)


def _line_to(x: float) -> LineString:
    """From the middle of cell 5 west to x, along the middle of the row holding 4 and 5.

    So nothing here lies ON a face: the only thing crossing one is the end of the line.
    """

    return LineString([INSIDE_CELL_5, (x, 37.5)])


def _drn(tmp_path, geometry, context=None, **fields):
    return _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, geometry), layer="bc", context=context or _context(), nper=1, **fields
        )
    )


def _cells_of(records):
    return sorted(record[0][1] for record in records)


@pytest.fixture
def polys():
    return _context().grid.gdf_vorPolys


@pytest.fixture
def pair():
    """Two 25 ft cells meeting at x=0, by hand, because the boundary tests need EXACT.

    A piece measured against the tolerance itself has to BE the tolerance: on the
    grid above ``25.0 - 24.99`` is 0.010000000000001563, while a difference taken
    from zero is the same double as the ``0.01`` literal the rule compares with.
    """

    return gpd.GeoDataFrame(
        geometry=[box(-25.0, 0.0, 0.0, 25.0), box(0.0, 0.0, 25.0, 25.0)], crs="EPSG:2927"
    )


# --- the overshoot, on the rule itself --------------------------------------- #


def test_a_line_ending_exactly_on_the_face_reaches_the_one_cell(polys):
    """The control: no overshoot, no second cell, before or after this change."""

    assert cell_shares(polys, _line_to(25.0)) == [(5, 1.0)]


@pytest.mark.parametrize("overshoot", [1e-6, 1e-3], ids=["1e-6", "1e-3"])
def test_an_overshoot_under_the_floor_is_noise_not_a_second_cell(polys, overshoot):
    """Measured at HEAD: both gave cell 4 a full share of its own and a full record."""

    shares = cell_shares(polys, _line_to(25.0 - overshoot))
    assert [cell for cell, _ in shares] == [5]
    # A share is a fraction of the WHOLE feature, so the dropped tip is simply not
    # applied -- ``Spread``'s "clip", and what the areal branch does with a sliver.
    assert shares[0][1] == pytest.approx(12.5 / (12.5 + overshoot))


def test_an_overshoot_past_the_floor_is_a_real_reach_with_its_share(polys):
    """0.015 ft in is past the tolerance -- and is what SFR keeps as a reach too."""

    shares = dict(cell_shares(polys, _line_to(25.0 - 0.015)))
    assert shares == pytest.approx({4: 0.015 / 12.515, 5: 12.5 / 12.515})


def test_a_crs_round_trip_does_not_add_a_cell(polys):
    """2927 -> 4326 -> 2927 on a State Plane grid moved the end by -2.9e-7 ft."""

    context = _context(x0=1_300_000.0, y0=650_000.0)
    line = LineString([(1_300_037.5, 650_037.5), (1_300_025.0, 650_037.5)])
    moved = gpd.GeoSeries([line], crs="EPSG:2927").to_crs("EPSG:4326").to_crs("EPSG:2927").iloc[0]
    assert moved.coords[-1][0] != 1_300_025.0, "the round trip moved nothing; this tests nothing"
    assert [cell for cell, _ in cell_shares(context.grid.gdf_vorPolys, moved)] == [5]


def test_the_shared_selector_takes_the_floor_too(polys):
    """The floor is the SHARED rule's, so it moves every caller ledger 235 lists --
    regions and ``map(select=...)``, particle release, cross-sections,
    ``adjust_top_btm_overlaps``, the legacy inactive-cell set, and a line head
    target. Measured at HEAD: ``[[4, 5]]``."""

    features = gpd.GeoSeries([_line_to(25.0 - 1e-6)], crs=polys.crs)
    assert get_vor_cells_as_series(polys, features).tolist() == [[5]]


# --- what the floor must NOT take -------------------------------------------- #


def test_a_line_along_a_face_still_reaches_both_cells_in_halves(polys):
    """Ledger 212, unchanged: the stretch is in both cells' pieces and is split."""

    shares = dict(cell_shares(polys, LineString([(0.0, 25.0), (100.0, 25.0)])))
    assert shares == pytest.approx({cell: 0.125 for cell in range(8)})
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_dropped_sliver_does_not_disturb_another_pair_sharing_a_face(polys):
    """20 ft along the face between cells 4 and 5, then a 1.4e-6 ft tail into cell 8.

    Dropping cell 8 must take it out of the FACE bookkeeping as well as out of the
    sizes: the pair 4/5 still has to give back half the stretch it shares. Clearing
    that bookkeeping instead leaves both at the full 20 ft -- 0.99999993 each,
    summing to 2.0, which is exactly the doubling ledger 212 removed.
    """

    shares = dict(cell_shares(polys, LineString([(25.0, 30.0), (25.0, 50.0), (24.999999, 50.000001)])))
    assert sorted(shares) == [4, 5]
    assert shares[4] == pytest.approx(0.5)
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_line_that_is_all_slivers_is_the_feature_and_is_kept(polys):
    """0.017 ft long across the vertex at (25, 25): two pieces of 0.0085, each under
    the floor. Dropping every sliver would lose the whole line, so they go only when
    something longer survives them."""

    shares = dict(cell_shares(polys, LineString([(24.994, 24.994), (25.006, 25.006)])))
    assert shares == pytest.approx({0: 0.5, 5: 0.5})


def test_a_line_shorter_than_the_floor_is_never_chopped():
    """The line analogue of ``thin_feature``: on 1 ft cells the sliver share is 0.005
    ft, so the 0.006 piece would be kept and the 0.003 one dropped -- two thirds of a
    0.009 ft feature that is entirely what was drawn. A feature this short is not a
    sliver of anything."""

    polys = _context(length=4.0).grid.gdf_vorPolys
    shares = dict(cell_shares(polys, LineString([(0.994, 0.5), (1.003, 0.5)])))
    assert shares == pytest.approx({0: 6.0 / 9.0, 1: 3.0 / 9.0})
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_short_piece_that_is_a_real_share_of_its_cell_is_coverage():
    """The floor needs its share guard here as it does for a polygon: on 0.02 ft cells
    a piece of 0.005 ft is half the cell's mean width -- a quarter of its side, real
    coverage of cell 11, and under the absolute tolerance. Without the guard it would
    be dropped, since the line's other pieces are longer and the all-slivers escape
    does not apply."""

    polys = _context(n=10, length=0.2).grid.gdf_vorPolys
    shares = dict(cell_shares(polys, LineString([(0.0, 0.021), (0.05, 0.019)])))
    assert sorted(shares) == [1, 2, 10, 11]
    assert shares[11] == pytest.approx(0.005 / 0.05, rel=1e-3)
    assert sum(shares.values()) == pytest.approx(1.0)


def test_a_sliver_dropped_off_a_face_leaves_the_whole_length_with_the_other_cell():
    """4 thousandths along the face x=25 and then 15 ft into cell 5. Cell 4's piece is
    noise; the stretch it shared is not halved out of cell 5 once it is gone."""

    polys = _context().grid.gdf_vorPolys
    line = LineString([(25.0, 25.0), (25.0, 25.004), (40.0, 30.0)])
    shares = cell_shares(polys, line)
    assert [cell for cell, _ in shares] == [5]
    assert shares[0][1] == pytest.approx(1.0)


# --- which measure the guard divides by, and where its edges are -------------- #


def test_the_guard_is_a_share_of_the_cells_width_not_its_perimeter_or_the_feature():
    """On 0.02 ft cells the guards differ by decades: 1% of the cell's mean width is
    1e-4, of its perimeter 8e-4, of this 0.0203 ft feature 2.03e-4. A spike 1.5e-4 ft
    into cell 10 is 1.5 times the cell's own guard and is coverage; measured against
    either of the others it would be dropped."""

    polys = _context(n=10, length=0.2).grid.gdf_vorPolys
    spike = LineString([(0.001, 0.010), (0.001, 0.02 + 7.5e-5), (0.003, 0.010)])
    shares = dict(cell_shares(polys, spike))
    assert sorted(shares) == [0, 10]
    assert shares[10] == pytest.approx(1.515e-4 / spike.length, rel=1e-3)


def test_a_line_piece_is_a_length_against_a_width_never_against_an_area():
    """The same grid, a spike of 2e-5 ft: under the cell's 1e-4 guard, so noise. 1% of
    the cell's AREA is 4e-6 -- a literal copy of the areal branch's guard would keep
    it, comparing a length with an area."""

    polys = _context(n=10, length=0.2).grid.gdf_vorPolys
    spike = LineString([(0.001, 0.010), (0.001, 0.02 + 1e-5), (0.003, 0.010)])
    assert [cell for cell, _ in cell_shares(polys, spike)] == [0]


def test_a_piece_exactly_the_tolerance_long_is_coverage(pair):
    """SFR drops a reach SHORTER than the tolerance (ledger 214), and the floor here
    is the same comparison: 0.01 ft of line in a cell is the shortest piece that is
    still coverage. The face sits at x=0 so the piece length is exactly the ``0.01``
    the rule holds."""

    line = LineString([(12.5, 12.5), (-0.01, 12.5)])
    assert line.intersection(pair.geometry.loc[0]).length == 0.01
    shares = dict(cell_shares(pair, line))
    assert sorted(shares) == [0, 1]
    assert shares[0] == pytest.approx(0.01 / 12.51)


def test_a_feature_exactly_the_tolerance_long_is_what_was_drawn(pair):
    """``thin_feature`` is inclusive for the same reason: a feature OF the tolerance
    is not a sliver of anything. 0.008 ft along the face and 0.002 into cell 1, so
    cell 1's piece is the whole 0.01 and cell 0's is only the shared stretch -- the
    one arrangement where the escape does not already cover this, since the pieces
    overlap. Exclusive, cell 0's 0.008 becomes a sliver and the record disappears."""

    feature = LineString([(0.0, 0.0), (0.0, 0.008), (0.002, 0.008)])
    assert feature.length == 0.01
    # Halved out of both, as ledger 212 asks: 0.008 - 0.004 and 0.01 - 0.004.
    assert dict(cell_shares(pair, feature)) == pytest.approx({0: 0.4, 1: 0.6})


# --- what a .gpkg source writes ---------------------------------------------- #


@pytest.mark.parametrize("package", LIST_BCS)
def test_an_overshooting_line_writes_one_record_in_every_list_bc(tmp_path, package):
    spec = getattr(mf, package).gpkg(
        _gpkg(tmp_path, _line_to(25.0 - 1e-6)), layer="bc", context=_context(), nper=1
    )
    assert _cells_of(spec.options["stress_period_data"][0]) == [5]


def test_the_second_record_carried_the_whole_conductance(tmp_path):
    """The record that is gone was not a small one: an unsplit field is written
    whole, so cell 4 held the same 100.0 conductance and the same 10.0 elevation."""

    records = _drn(tmp_path, _line_to(25.0 - 1e-6))
    assert [(record[0][1], record[1], record[2]) for record in records] == [(5, 10.0, 100.0)]


def test_a_spread_conductance_is_no_longer_split_with_the_cell_beyond(tmp_path):
    """``Spread``'s own ``min_share`` hid this for a wrapped field, but not with the
    floor turned off, and not at all for the plain column above. Under "retained" the
    whole conductance stays in cell 5 instead of leaking 0.008 of it into cell 4."""

    records = _drn(
        tmp_path, _line_to(25.0 - 1e-3),
        conductance=Spread("conductance", mode="retained", min_share=0.0),
    )
    assert {record[0][1]: record[2] for record in records} == pytest.approx({5: 100.0})


def test_a_diagonal_through_vertices_after_a_round_trip_keeps_its_four_cells(tmp_path):
    """The geometry of ledger 214's own SFR test: a diagonal drawn on grid vertices,
    sent to EPSG:4326 and back on a State Plane grid, clips three corner cells by
    ~9e-7 ft. SFR dropped those three micro-reaches; ``.gpkg`` wrote SEVEN DRN records
    where four are right -- 700.0 of conductance against 400.0, since every record
    carries the whole 100.0."""

    context = _context(x0=1_300_000.0, y0=400_000.0)
    line = LineString([(1_300_000.0, 400_000.0), (1_300_100.0, 400_100.0)])
    moved = gpd.GeoSeries([line], crs="EPSG:2927").to_crs("EPSG:4326").to_crs("EPSG:2927").iloc[0]
    records = _drn(tmp_path, moved, context)
    assert _cells_of(records) == [0, 5, 10, 15]
    assert sum(record[2] for record in records) == pytest.approx(400.0)


def test_the_gpkg_and_sfr_now_agree_about_one_geometry(tmp_path):
    """They did not: SFR kept one reach in cell 5 (ledger 214) where ``.gpkg`` wrote
    a second record for cell 4."""

    overshooting = _line_to(25.0 - 1e-6)
    streams = gpd.GeoDataFrame({"name": ["s0"]}, geometry=[overshooting], crs="EPSG:2927")
    sfr = mf.sfr(
        context=_context(), nper=1, streams=streams, stream_id="name", reach_top=10.0, width=10.0
    )
    reach_cells = sorted(int(record[1][1]) for record in sfr.options["packagedata"])
    assert reach_cells == [5]
    assert _cells_of(_drn(tmp_path, overshooting)) == reach_cells
