"""``mf.sfr(reverse_streams=)`` reads a single string as ONE stream id.

A bare string was iterated character by character, so ``reverse_streams="main"``
reversed streams ``m``, ``a``, ``i`` and ``n`` -- none of them -- and a stream digitized
mouth-first stayed reversed in silence, its reaches numbered outlet-first.
"""

from __future__ import annotations

import geopandas as gpd
from shapely.geometry import LineString

import myflopy as mf
from tests.test_gpkg_on_cell_faces import _context


def _cells(reverse_streams):
    streams = gpd.GeoDataFrame(
        {"name": ["main"]}, geometry=[LineString([(100.0, 37.5), (0.0, 37.5)])], crs="EPSG:2927"
    )
    spec = mf.sfr(
        context=_context(), nper=1, streams=streams, stream_id="name", reach_top=10.0,
        width=10.0, reverse_streams=reverse_streams,
    )
    return [int(record[1][1]) for record in spec.options["packagedata"]]


def test_a_single_string_reverses_that_stream():
    # Drawn east to west, so reversed it runs west to east: cells 4, 5, 6, 7.
    assert _cells("main") == [4, 5, 6, 7]


def test_the_string_reads_as_the_one_item_tuple():
    assert _cells("main") == _cells(("main",))
    assert _cells(()) == [7, 6, 5, 4]
