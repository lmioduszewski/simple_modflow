"""A stream centerline lying ON a face goes to one cell: the lower-numbered.

``SFRBuilder.reaches`` intersected each cell with the centerline, and a closed cell
polygon includes its boundary, so a centerline on a shared face came back full
length from BOTH cells. A stream along a row line of a 4x4 grid gave 8 reaches and
200 ft of rlen for a 100 ft stream -- streambed conductance doubled -- and the tied
midpoints routed the chain back and forth between the two rows. MF6 accepts every
record. Decided 2026-09-11: the lower-numbered cell takes a face stretch; halving it
between both would make two parallel reaches and break the chain.
"""

from __future__ import annotations

import geopandas as gpd
import pytest
from shapely.geometry import LineString

import myflopy as mf
from tests.test_gpkg_on_cell_faces import _context


def _sfr(context, *coords):
    streams = gpd.GeoDataFrame({"name": ["s0"]}, geometry=[LineString(coords)], crs="EPSG:2927")
    return mf.sfr(
        context=context, nper=1, streams=streams, stream_id="name", reach_top=10.0, width=10.0
    )


def _reaches(spec):
    """(cell, rlen) per reach, in reach order."""

    return [(int(record[1][1]), float(record[2])) for record in spec.options["packagedata"]]


def test_a_stream_through_cell_interiors_is_unchanged():
    reaches = _reaches(_sfr(_context(), (0.0, 37.5), (100.0, 37.5)))
    assert [cell for cell, _ in reaches] == [4, 5, 6, 7]
    assert sum(rlen for _, rlen in reaches) == pytest.approx(100.0)


@pytest.mark.parametrize(
    ("coords", "cells"),
    [
        (((0.0, 25.0), (100.0, 25.0)), [0, 1, 2, 3]),
        (((50.0, 0.0), (50.0, 100.0)), [1, 5, 9, 13]),
    ],
    ids=["along-a-row-line", "along-a-column-line"],
)
def test_a_stream_on_a_face_is_one_reach_per_stretch_in_the_lower_cell(coords, cells):
    reaches = _reaches(_sfr(_context(), *coords))
    assert [cell for cell, _ in reaches] == cells
    assert sum(rlen for _, rlen in reaches) == pytest.approx(100.0)


def test_the_chain_follows_the_stream_instead_of_zigzagging_between_rows():
    spec = _sfr(_context(), (0.0, 25.0), (100.0, 25.0))
    downstream = {
        reach: [-other for other in links if other < 0]
        for reach, *links in (list(row) for row in spec.options["connectiondata"])
    }
    assert downstream == {0: [1], 1: [2], 2: [3], 3: []}


def test_a_stream_partly_on_a_face_counts_that_stretch_once():
    """Up through cell 0, along y=25 for 25 ft, up through cell 5: 75 ft, not 100."""

    reaches = _reaches(
        _sfr(_context(), (12.5, 0.0), (12.5, 25.0), (37.5, 25.0), (37.5, 50.0))
    )
    assert reaches == [(0, pytest.approx(37.5)), (1, pytest.approx(12.5)), (5, pytest.approx(25.0))]


def test_a_crs_round_trip_leaves_no_micro_reaches():
    """A diagonal through vertices, sent to EPSG:4326 and back, clipped three corner
    cells by ~9e-7 ft: three micro-reaches in the chain. MF6 has crashed on less."""

    context = _context(x0=1_300_000.0, y0=400_000.0)
    line = LineString([(1_300_000.0, 400_000.0), (1_300_100.0, 400_100.0)])
    moved = gpd.GeoSeries([line], crs="EPSG:2927").to_crs("EPSG:4326").to_crs("EPSG:2927")
    streams = gpd.GeoDataFrame({"name": ["s0"]}, geometry=moved, crs="EPSG:2927")
    spec = mf.sfr(
        context=context, nper=1, streams=streams, stream_id="name", reach_top=10.0, width=10.0
    )
    reaches = _reaches(spec)
    assert [cell for cell, _ in reaches] == [0, 5, 10, 15]
    assert min(rlen for _, rlen in reaches) > 1.0


def test_a_face_stretch_a_lower_cell_drops_stays_with_the_cell_that_keeps_it():
    """3 ft along the face under cell 4, then 20 ft up into it, with min_reach_length=5:
    cell 0's 3 ft is dropped, so cell 4 keeps all 23 ft rather than losing 3 to nobody."""

    streams = gpd.GeoDataFrame(
        {"name": ["s0"]},
        geometry=[LineString([(0.0, 25.0), (3.0, 25.0), (3.0, 45.0)])],
        crs="EPSG:2927",
    )
    spec = mf.sfr(
        context=_context(), nper=1, streams=streams, stream_id="name", reach_top=10.0,
        width=10.0, min_reach_length=5.0,
    )
    assert _reaches(spec) == [(4, pytest.approx(23.0))]
