"""A tributary enters the main stem where it MEETS the main stem.

``receiver_location="nearest"`` -- the default -- projected the tributary's
``representative_point()`` onto the main stem, and for a two-vertex tributary that
point is its upstream START. A tributary coming in from the east and joining
mid-grid was connected to the main stem's LAST reach: its water entered downstream
of the confluence and skipped every reach between. Found by the critic of the
vector-to-cell sweep, 2026-09-11; MF6 accepts the connection as written.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString

from myflopy.modflow.mf6.sfr import SFRBuilder, StreamConnection, StreamDiversion
from myflopy.specs import ModelContext
from tests.test_gpkg_on_cell_faces import _context

MAIN = LineString([(0.0, 37.5), (100.0, 37.5)])  # through cells 4, 5, 6, 7


def _builder(*others, connections=(), **options):
    context = _context()
    grid = context.grid
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: np.full(16, 100.0), 1: np.zeros(16)}, geometry=grid.gdf_vorPolys.geometry, crs=grid.crs
    )
    names = ["main"] + [f"s{i}" for i in range(len(others))]
    streams = gpd.GeoDataFrame({"sid": names}, geometry=[MAIN, *others], crs="EPSG:2927")
    return SFRBuilder(
        context=ModelContext(grid=grid, domain=np.ones((1, 16), dtype=int)), nper=1,
        streams=streams, stream_id="sid", connection_tolerance=1.0, connections=connections,
        **options,
    )


def _receiver_cell(builder):
    (_, receiver), = builder.resolved_connections
    return int(builder.reaches.loc[receiver, "cellid"][1])


def test_a_tributary_joining_mid_cell_enters_the_reach_it_joins():
    """From (95, 100) to (37.5, 37.5): it joins in cell 5, not at the main stem's end."""

    assert _receiver_cell(_builder(LineString([(95.0, 100.0), (37.5, 37.5)]))) == 5


def test_a_tributary_joining_on_a_reach_boundary_enters_the_downstream_reach():
    """It joins at x=50, between the reaches in cells 5 and 6. Water joining there
    enters cell 6's reach at its upstream end; cell 5's would carry it above the
    confluence."""

    assert _receiver_cell(_builder(LineString([(90.0, 100.0), (50.0, 37.5)]))) == 6


def test_the_receiver_is_found_from_where_the_source_discharges():
    """A stream running beside the main stem, discharging at its downstream end,
    enters the main stem there -- not at whatever point of the two lines is nearest."""

    side = LineString([(0.0, 62.5), (100.0, 62.5)])
    builder = _builder(side, connections=(StreamConnection("s0", "main"),))
    assert _receiver_cell(builder) == 7


def test_a_diversion_is_delivered_where_it_leaves_its_source():
    """Diversions resolve their receiver the same way, and had the same bug: water
    leaving the main stem's downstream end was delivered to the canal's START."""

    canal = LineString([(0.0, 62.5), (100.0, 62.5)])  # through cells 8, 9, 10, 11
    builder = _builder(
        canal,
        connection_mode="explicit",
        connections=(StreamConnection("main", None), StreamConnection("s0", None)),
        diversions=(StreamDiversion("main", "s0", amount={0: 0.25}),),
    )
    (_, _, receiver, _), = builder.resolved_diversions
    assert int(builder.reaches.loc[receiver, "cellid"][1]) == 11
