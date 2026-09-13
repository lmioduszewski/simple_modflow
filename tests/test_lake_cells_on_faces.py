"""A lake COVERS its cells; it does not keep every cell it grazes.

``LAKBuilder.lake_cells`` kept any cell whose overlap with the lake polygon had an
area above zero, with no noise floor. On a 4x4 grid of 25 ft cells a lake drawn
exactly on cell 5 gave ``[5]``, but the same lake on a State Plane grid round
tripped 2927 -> 4326 -> 2927 gave ``[4, 5, 8, 9]`` -- the extras overlapping by
7.41e-06, 9.07e-14 and 7.41e-06 square feet -- via EPSG:2285 ``[5, 6, 9, 10]``,
and grown by 1e-06 ft all nine cells around it. Each of those becomes a full LAK
cell with its own vertical connection, so the lake's sidewall topology and its
budget moved with the projection it had been drawn in.

``LAKBuilder._active`` was a second copy of the shared idomain test that counted a
ONE-dimensional domain as active everywhere and raised ``IndexError`` for a layer
beyond the domain.
"""

from __future__ import annotations

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import MultiPolygon, box

from myflopy.modflow.mf6 import lakes as lakes_module
from myflopy.modflow.mf6.grid import geometry as grid_geometry
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.lakes import LAKBuilder
from myflopy.specs import ModelContext

STATE_PLANE = (1.3e6, 6.5e5)  # a real EPSG:2927 location, where the precision goes
KEEP = object()  # "leave the default alone", so a test can ask for domain=None


def _grid(*, x0: float = 0.0, y0: float = 0.0, n: int = 4, length: float = 100.0):
    """n x n square cells, row-major from the bottom left; cell 5 is the (25..50)^2 box."""

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
    ncpl = n * n
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [20.0] * ncpl, 1: [10.0] * ncpl, 2: [0.0] * ncpl},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    grid.get_disu_connectivity(validate=True)
    return grid


def _builder(geometry, *, grid=None, domain=KEEP, **updates) -> LAKBuilder:
    grid = _grid() if grid is None else grid
    ncpl = len(grid.iverts)
    values = {
        "context": ModelContext(
            grid=grid,
            domain=np.ones((2, ncpl), dtype=int) if domain is KEEP else domain,
        ),
        "nper": 1,
        "lakes": gpd.GeoDataFrame(
            {"name": ["pond"]}, geometry=[geometry], crs="EPSG:2927"
        ),
        "lake_id_field": "name",
        "starting_stage": 18.0,
        "lake_bottom": 15.0,
        "bed_leakance": 0.1,
    }
    values.update(updates)
    return LAKBuilder(**values)


def _round_trip(geometry, via: str = "EPSG:4326"):
    """The geometry after a projection round trip, which is where the slivers come from."""

    return gpd.GeoSeries([geometry], crs="EPSG:2927").to_crs(via).to_crs("EPSG:2927").iloc[0]


def _cell_5(x0: float = 0.0, y0: float = 0.0):
    return box(x0 + 25.0, y0 + 25.0, x0 + 50.0, y0 + 50.0)


# --------------------------------------------------------------------------
# The four measured cases.
# --------------------------------------------------------------------------


def test_a_lake_drawn_on_one_cell_is_one_cell():
    """The control: exact coordinates already gave the right answer."""

    assert _builder(_cell_5()).lake_cells == {"pond": [5]}


def test_a_round_trip_through_4326_does_not_add_three_cells():
    """Was [4, 5, 8, 9]: overlaps of 7.41e-06, 9.07e-14 and 7.41e-06 square feet."""

    x0, y0 = STATE_PLANE
    grid = _grid(x0=x0, y0=y0)
    lake = _round_trip(_cell_5(x0, y0))

    polys = grid.gdf_vorPolys.geometry
    grazed = polys.intersection(lake)
    assert sorted(int(c) for c in grazed.index[grazed.area > 0.0]) == [4, 5, 8, 9]
    assert max(float(grazed.loc[cell].area) for cell in (4, 8, 9)) < 1e-5

    assert _builder(lake, grid=grid).lake_cells == {"pond": [5]}


def test_a_round_trip_through_2285_does_not_add_three_cells():
    """Was [5, 6, 9, 10]: overlaps of 1.41e-04, 8.58e-03 and 1.94e-09 square feet.

    The largest of these is 8.58e-03 ft2 -- above ``_SLIVER_WIDTH`` as a number, so
    it is the WIDTH of the piece and its share of the cell that reject it, not its
    area against the tolerance.
    """

    x0, y0 = STATE_PLANE
    grid = _grid(x0=x0, y0=y0)
    lake = _round_trip(_cell_5(x0, y0), "EPSG:2285")

    polys = grid.gdf_vorPolys.geometry
    grazed = polys.intersection(lake)
    assert sorted(int(c) for c in grazed.index[grazed.area > 0.0]) == [5, 6, 9, 10]
    assert max(float(grazed.loc[cell].area) for cell in (6, 9, 10)) > 1e-3

    assert _builder(lake, grid=grid).lake_cells == {"pond": [5]}


def test_a_lake_grown_by_a_millionth_of_a_foot_does_not_reach_nine_cells():
    """Was all nine cells around cell 5, the corner ones overlapping by 1e-12 ft2."""

    grown = box(25.0 - 1e-6, 25.0 - 1e-6, 50.0 + 1e-6, 50.0 + 1e-6)

    polys = _grid().gdf_vorPolys.geometry
    grazed = polys.intersection(grown)
    assert len(grazed.index[grazed.area > 0.0]) == 9

    assert _builder(grown).lake_cells == {"pond": [5]}


# --------------------------------------------------------------------------
# What the spurious cells cost, and what the rule must not take away.
# --------------------------------------------------------------------------


def test_each_spurious_cell_was_a_full_lak_cell_with_its_own_vertical_connection():
    """The consequence: connection count and topology, not just a list of ids."""

    x0, y0 = STATE_PLANE
    grid = _grid(x0=x0, y0=y0)
    builder = _builder(_round_trip(_cell_5(x0, y0)), grid=grid)

    cellids = [connection.cellid for _lake, connection in builder.connections]
    vertical = [
        connection.cellid
        for _lake, connection in builder.connections
        if connection.connection_type == "VERTICAL"
    ]
    assert vertical == [(0, 5)]
    assert {cell for _layer, cell in cellids} == {5}
    assert builder.packagedata[0][2] == len(cellids)


def test_a_lake_reaching_a_tenth_of_a_foot_into_the_next_cell_still_covers_it():
    """The floor is the PIECE's width, not its share of the lake.

    A tenth of a foot into cell 6 is 2.5 ft2, 0.4% of the lake and 0.4% of the
    cell -- but the strip is 0.0996 ft wide, ten times ``_SLIVER_WIDTH``, so it is
    a real clip and the cell is a real lake cell. A rule that went by share of the
    lake would drop it, and the lake would lose a sidewall it has.
    """

    lake = box(25.0, 25.0, 50.0 + 0.1, 50.0)

    assert _builder(lake).lake_cells == {"pond": [5, 6]}


def test_a_lake_that_is_itself_that_thin_keeps_every_piece():
    """The shared rule keeps a feature that IS a sliver: they are the lake, not noise.

    A 0.004 ft strip centred on the face between cells 1 and 5 is two 0.002 ft
    halves -- both under the tolerance, so dropping slivers would empty it.
    """

    strip = box(30.0, 25.0 - 0.002, 45.0, 25.0 + 0.002)

    assert _builder(strip).lake_cells == {"pond": [1, 5]}


def test_a_detached_part_of_a_multipart_lake_is_judged_on_its_own():
    """The one cell a modeller might MEAN that the floor takes, pinned deliberately.

    The measure is per PIECE, and a detached part is its own piece: a 15 x 0.005 ft
    satellite of a one-cell lake (0.075 ft2, mean width 0.005) is under the floor
    on its own even though the lake as a whole is 12.5 ft wide, so its cell goes.
    A 1 x 1 ft pond and a 15 x 0.1 ft strip stay, so the loss is confined to parts
    that are themselves noise -- and this is the shared rule, the same answer a
    ``.gpkg`` zone of that shape gets.
    """

    lake = MultiPolygon([_cell_5(), box(80.0, 80.0, 95.0, 80.005)])
    polys = _grid().gdf_vorPolys.geometry
    grazed = polys.intersection(lake)
    assert sorted(int(c) for c in grazed.index[grazed.area > 0.0]) == [5, 15]

    assert _builder(lake).lake_cells == {"pond": [5]}

    pond = MultiPolygon([_cell_5(), box(80.0, 80.0, 81.0, 81.0)])
    assert _builder(pond).lake_cells == {"pond": [5, 15]}
    strip = MultiPolygon([_cell_5(), box(80.0, 80.0, 95.0, 80.1)])
    assert _builder(strip).lake_cells == {"pond": [5, 15]}


def test_a_multipart_lake_whose_every_part_is_thin_keeps_them_all():
    """Both parts under the floor is the whole-feature escape, so nothing is lost.

    The first part straddles the 1|5 face, so this is three cells' worth of pieces
    and every one of them is a sliver on its own.
    """

    lake = MultiPolygon(
        [box(30.0, 25.0 - 0.002, 45.0, 25.0 + 0.002), box(80.0, 80.0, 95.0, 80.005)]
    )

    assert _builder(lake).lake_cells == {"pond": [1, 5, 15]}


def test_a_lake_that_only_borders_the_grid_reaches_none_and_says_so():
    """A lake beyond the edge, touching it along a face, covers nothing (ledger 210)."""

    outside = box(100.0, 25.0, 125.0, 50.0)  # abuts the grid's east edge at x=100
    builder = _builder(outside)

    assert builder.lake_cells == {"pond": []}
    with pytest.raises(ValueError, match="only borders cells"):
        builder.validate()


def test_a_collapsed_lake_outline_is_refused_rather_than_measured_as_a_line():
    """A polygon with no area would fall to the shared rule's LENGTH branch.

    ``box(25, 25, 25, 50)`` is the face between cells 4 and 5: still a ``Polygon``,
    so the geometry-type check passes, and measured by length it reaches both.
    """

    with pytest.raises(ValueError, match="Lake geometry has no area: pond"):
        _builder(box(25.0, 25.0, 25.0, 50.0)).lake_cells


# --------------------------------------------------------------------------
# A lake that reaches nothing follows ledger 211.
# --------------------------------------------------------------------------


def test_a_lake_off_the_grid_raises_with_both_extents_and_both_crss():
    builder = _builder(box(1_000.0, 1_000.0, 1_050.0, 1_050.0))

    assert builder.lake_cells == {"pond": []}
    with pytest.raises(ValueError) as error:
        builder.validate()
    message = str(error.value)
    assert "1 of 1 lakes cover no grid cell" in message
    assert "'pond' (outside the grid)" in message
    assert "1,000.0" in message and "0.0, 0.0, 100.0, 100.0" in message
    assert "EPSG:2927" in message


def test_a_lake_file_in_degrees_names_the_reprojection_it_went_through():
    """The lakes are reprojected to the grid, so the miss must report the SOURCE CRS."""

    builder = _builder(
        box(-122.5, 47.5, -122.4, 47.6),
        lakes=gpd.GeoDataFrame(
            {"name": ["pond"]},
            geometry=[box(-122.5, 47.5, -122.4, 47.6)],
            crs="EPSG:4326",
        ),
    )

    with pytest.raises(ValueError, match="reprojected from EPSG:4326 to the grid's EPSG:2927"):
        builder.validate()


def test_one_lake_of_two_missing_still_raises_naming_only_that_lake():
    """No partial tier: a lake with no cells carries no connection, which MF6 refuses."""

    builder = _builder(
        _cell_5(),
        lakes=gpd.GeoDataFrame(
            {"name": ["pond", "gone"]},
            geometry=[_cell_5(), box(1_000.0, 1_000.0, 1_050.0, 1_050.0)],
            crs="EPSG:2927",
        ),
    )

    with pytest.raises(ValueError) as error:
        builder.validate()
    assert "1 of 2 lakes cover no grid cell" in str(error.value)
    assert "'gone'" in str(error.value) and "'pond'" not in str(error.value)


# --------------------------------------------------------------------------
# One activity test, not two.
# --------------------------------------------------------------------------


def test_a_one_dimensional_domain_is_honoured_per_cell():
    """It counted as active EVERYWHERE, so a lake connected through switched-off cells."""

    domain = np.ones(16, dtype=int)
    domain[5] = 0
    builder = _builder(_cell_5(), domain=domain)

    assert builder._active(0, 4) is True
    assert builder._active(0, 5) is False
    with pytest.raises(ValueError, match="does not intersect active cell 5"):
        builder.connectiondata


def test_a_passthrough_cell_is_not_active_in_either_domain_shape():
    """MF6 removes an idomain -1 cell; ``bool(-1)`` is True, which is why ``> 0``."""

    flat = np.ones(16, dtype=int)
    flat[5] = -1
    assert _builder(_cell_5(), domain=flat)._active(0, 5) is False

    layered = np.ones((2, 16), dtype=int)
    layered[0, 5] = -1
    builder = _builder(_cell_5(), domain=layered)
    assert builder._active(0, 5) is False
    assert builder._active(1, 5) is True


def test_a_layer_beyond_the_domain_names_the_layer():
    """It raised ``IndexError: index 3 is out of bounds for axis 0 with size 2``."""

    builder = _builder(_cell_5())

    with pytest.raises(ValueError, match="Layer 3 is outside the model domain"):
        builder._active(3, 5)


def test_no_domain_is_a_wholly_active_grid():
    assert _builder(_cell_5(), domain=None)._active(7, 5) is True


# --------------------------------------------------------------------------
# The bathymetry route (ledger 231) keeps working, and asks for less.
# --------------------------------------------------------------------------


def test_a_per_cell_bed_need_not_cover_the_sliver_cells_any_more():
    """Ledger 231 raises for a lake cell a Series misses; the slivers were such cells."""

    x0, y0 = STATE_PLANE
    grid = _grid(x0=x0, y0=y0)
    builder = _builder(
        _round_trip(_cell_5(x0, y0)),
        grid=grid,
        lake_bottom=pd.Series({5: 15.0}),
        connection_modes="bathy",
    )

    assert builder.lake_cells == {"pond": [5]}
    builder.validate()  # the Series covers cell 5 and only cell 5


def test_a_bathy_lake_still_steps_down_toward_its_deeper_neighbour():
    """A real two-cell bathy lake still builds sidewalls, so the rule took nothing."""

    lake = box(25.0, 25.0, 75.0, 50.0)  # cells 5 and 6
    builder = _builder(
        lake,
        lake_bottom=pd.Series({4: 20.0, 5: 15.0, 6: 12.0, 7: 20.0, 1: 20.0, 2: 20.0, 9: 20.0, 10: 20.0}),
        connection_modes="bathy",
    )

    assert builder.lake_cells == {"pond": [5, 6]}
    horizontal = [
        connection
        for _lake, connection in builder.connections
        if connection.connection_type == "HORIZONTAL"
    ]
    assert horizontal, "a stepped bathy lake must still have sidewall connections"
    assert {connection.cellid[1] for connection in horizontal} <= {5, 6}


# --------------------------------------------------------------------------
# One copy of each rule.
# --------------------------------------------------------------------------


def test_lakes_uses_the_shared_rules_rather_than_a_copy():
    assert lakes_module.cell_shares is grid_geometry.cell_shares
    assert lakes_module.domain_active is grid_geometry.domain_active
    assert lakes_module.miss_reason is grid_geometry.miss_reason
    assert lakes_module.where_they_are is grid_geometry.where_they_are


def test_lake_cells_goes_through_the_shared_helper(monkeypatch):
    """Patching the shared name must move the answer -- no second implementation."""

    monkeypatch.setattr(lakes_module, "cell_shares", lambda polys, geometry, **kw: [(11, 1.0)])

    assert _builder(_cell_5()).lake_cells == {"pond": [11]}
