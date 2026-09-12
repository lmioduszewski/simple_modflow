"""``mf.sfr_connection(at=(x, y))`` resolves through the cell holding the point.

Ledger 239. The old rule compared the point to each reach PIECE CENTROID, which
put a point 3 ft inside cell 5's own reach (9.5 ft from that reach's centroid)
into cell 6's reach, gave the water at a reach boundary to the UPSTREAM reach
(the opposite of ledger 215), and answered with the nearest reach for a point
nowhere near the grid.
"""

from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
import pytest
import yaml
from shapely.geometry import LineString

from myflopy import ModelContext, PackageSpec
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.sfr import SFRBuilder
from myflopy.package_api import _sfr_mover_index, sfr_connection

# A 4x4 grid of 25 ft square cells, cells numbered row by row from the origin:
# row 1 (y 25..50) is cells 4, 5, 6, 7, so a stream along y = 37.5 runs through
# their middles and the faces between them are x = 25, 50, 75.
CELLS = 4
WIDTH = 25.0


def _grid() -> VoronoiGridPlus:
    verts = np.array(
        [[i * WIDTH, j * WIDTH] for j in range(CELLS + 1) for i in range(CELLS + 1)],
        dtype=float,
    )
    iverts = [
        [j * (CELLS + 1) + i, (j + 1) * (CELLS + 1) + i,
         (j + 1) * (CELLS + 1) + i + 1, j * (CELLS + 1) + i + 1]
        for j in range(CELLS) for i in range(CELLS)
    ]
    centers = np.array(
        [[(i + 0.5) * WIDTH, (j + 0.5) * WIDTH] for j in range(CELLS) for i in range(CELLS)]
    )
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=centers, crs="EPSG:2927")


def _builder(grid: VoronoiGridPlus, *lines: LineString) -> SFRBuilder:
    streams = gpd.GeoDataFrame(
        {"name": [f"s{number}" for number in range(len(lines))]},
        geometry=list(lines),
        crs=grid.crs,
    )
    return SFRBuilder(
        context=ModelContext(grid=grid, domain=np.ones((1, CELLS * CELLS), dtype=int)),
        nper=1,
        streams=streams,
        stream_id="name",
        reach_top=10.0,
        width=10.0,
        gradient=0.001,
        roughness=0.03,
        streambed_k=0.1,
        streambed_thickness=1.0,
    )


def _spec(builder: SFRBuilder) -> PackageSpec:
    return PackageSpec(
        "sfr", flopy.mf6.ModflowGwfsfr, metadata={"sfr_index": _sfr_mover_index(builder)}
    )


def _row_spec(*lines: LineString) -> PackageSpec:
    return _spec(_builder(_grid(), *lines))


FULL_ROW = LineString([(0.0, 37.5), (100.0, 37.5)])      # reaches 0..3 in cells 4..7
SHORT_ROW = LineString([(0.0, 37.5), (60.0, 37.5)])      # reaches 0..2, the last 10 ft long
# Drawn mouth-first, so reach order runs AGAINST cell order -- the only way to
# tell "the downstream reach" from "the highest-numbered cell".
REVERSED_ROW = LineString([(100.0, 37.5), (0.0, 37.5)])  # reaches 0..3 in cells 7..4


def test_reach_table_is_the_one_these_cases_were_measured_on():
    """Guard the fixture: every expectation below is stated in reach numbers."""

    reaches = _builder(_grid(), FULL_ROW).reaches
    assert [int(cellid[1]) for cellid in reaches["cellid"]] == [4, 5, 6, 7]
    assert list(reaches["rno"]) == [0, 1, 2, 3]

    short = _builder(_grid(), SHORT_ROW).reaches
    assert [int(cellid[1]) for cellid in short["cellid"]] == [4, 5, 6]
    assert [float(value) for value in short["rlen"]] == [25.0, 25.0, 10.0]

    backwards = _builder(_grid(), REVERSED_ROW).reaches
    assert [int(cellid[1]) for cellid in backwards["cellid"]] == [7, 6, 5, 4]
    assert list(backwards["rno"]) == [0, 1, 2, 3]


def test_a_point_inside_a_reach_is_that_reach_not_the_nearest_centroid():
    """Measured: (47, 37.5) picked reach 2, whose cell-6 piece is centred at x = 55."""

    spec = _row_spec(SHORT_ROW)
    assert sfr_connection(spec, "s0", at=(47.0, 37.5)).index == 1


def test_every_reach_answers_for_a_point_inside_it():
    spec = _row_spec(FULL_ROW)
    for rno, x in enumerate((12.5, 37.5, 62.5, 87.5)):
        assert sfr_connection(spec, "s0", at=(x, 37.5)).index == rno
        # ...anywhere in the cell, not only on the stream line.
        assert sfr_connection(spec, "s0", at=(x, 30.0)).index == rno


def test_a_point_on_a_reach_boundary_goes_downstream():
    """Ledger 215: water joining at a reach boundary enters the DOWNSTREAM reach."""

    spec = _row_spec(FULL_ROW)
    assert sfr_connection(spec, "s0", at=(25.0, 37.5)).index == 1   # face 4|5, was 0
    assert sfr_connection(spec, "s0", at=(50.0, 37.5)).index == 2   # face 5|6, was 1
    assert sfr_connection(spec, "s0", at=(75.0, 37.5)).index == 3   # face 6|7, was 2


def test_downstream_is_down_the_stream_not_the_higher_cell_number():
    """A stream drawn mouth-first: its reach order runs against its cell order.

    On every left-to-right fixture the downstream reach also sits in the
    highest-numbered cell, so "downstream-most" and "largest cell id" agree and
    neither pins the other. Here they disagree: cell 6 holds reach 1 and cell 5
    holds reach 2, and ledger 215 wants reach 2.
    """

    spec = _row_spec(REVERSED_ROW)
    assert sfr_connection(spec, "s0", at=(50.0, 37.5)).index == 2   # face 6|5
    assert sfr_connection(spec, "s0", at=(25.0, 37.5)).index == 3   # face 5|4
    # A fractional coordinate is not truncated onto that face on the way in:
    # x = 25.4 is INSIDE cell 5, which is reach 2, not the 5|4 boundary's reach 3.
    assert sfr_connection(spec, "s0", at=(25.4, 37.5)).index == 2
    # ...and a point inside each of the four reaches, in reversed order.
    for rno, x in zip((3, 2, 1, 0), (12.5, 37.5, 62.5, 87.5), strict=True):
        assert sfr_connection(spec, "s0", at=(x, 37.5)).index == rno


def test_a_point_on_a_cell_face_away_from_the_stream_takes_the_reach_cell():
    """Only cells carrying a reach of this stream are eligible, as ``usable`` says."""

    spec = _row_spec(FULL_ROW)
    # x = 37.5 is inside cell 1 (below) and cell 5 (the reach); y = 25 is their face.
    assert sfr_connection(spec, "s0", at=(37.5, 25.0)).index == 1
    assert sfr_connection(spec, "s0", at=(37.5, 50.0)).index == 1
    # A vertex touches cells 5, 6, 9 and 10; two carry reaches, and downstream wins.
    assert sfr_connection(spec, "s0", at=(50.0, 50.0)).index == 2


def test_a_point_off_the_grid_raises_naming_the_extent_and_the_crs():
    """Measured: (500, 500) answered reach 3, the nearest centroid, in silence."""

    spec = _row_spec(FULL_ROW)
    with pytest.raises(ValueError, match="outside the grid"):
        sfr_connection(spec, "s0", at=(500.0, 500.0))
    with pytest.raises(ValueError, match=r"\(0.0, 0.0\) to \(100.0, 100.0\) in EPSG:2927"):
        sfr_connection(spec, "s0", at=(-5.0, 37.5))
    # A bare (x, y) carries no CRS of its own, so the message names the one it
    # has to be in -- a coordinate in the wrong frame is the commonest way here.
    with pytest.raises(ValueError, match="coordinate in the grid's own CRS"):
        sfr_connection(spec, "s0", at=(-5.0, 37.5))


def test_a_point_on_a_cell_with_no_reach_raises_naming_the_reach_cells():
    """Measured: (37.5, 12) -- cell 1, a row below the stream -- answered reach 1."""

    spec = _row_spec(FULL_ROW)
    with pytest.raises(ValueError, match="on no cell carrying one of its reaches"):
        sfr_connection(spec, "s0", at=(37.5, 12.0))
    with pytest.raises(ValueError, match=r"cells 4, 5, 6, 7"):
        sfr_connection(spec, "s0", at=(37.5, 12.0))
    # The nearest reach cell is named: it is what the old rule returned in silence.
    with pytest.raises(ValueError, match=r"cell 5, holding reach 1, centred at \(37.5, 37.5\)"):
        sfr_connection(spec, "s0", at=(37.5, 12.0))
    # The grid's own corner is INSIDE the grid -- it is the reach cells it is not
    # on -- so the extent test has to include the boundary.
    with pytest.raises(ValueError, match="on no cell carrying one of its reaches"):
        sfr_connection(spec, "s0", at=(0.0, 0.0))


def test_one_stream_does_not_borrow_another_streams_reach():
    """Two streams, one grid: a point in the other stream's cells is a miss, not a hit."""

    spec = _row_spec(FULL_ROW, LineString([(37.5, 0.0), (37.5, 20.0)]))
    assert sfr_connection(spec, "s1", at=(37.5, 12.0)).index == 4
    with pytest.raises(ValueError, match="on no cell carrying one of its reaches"):
        sfr_connection(spec, "s0", at=(37.5, 12.0))
    with pytest.raises(ValueError, match="on no cell carrying one of its reaches"):
        sfr_connection(spec, "s1", at=(62.5, 37.5))


def test_two_reaches_in_one_cell_go_downstream():
    """The tie ``cell_shares`` cannot see: both reaches are in the SAME cell."""

    spec = PackageSpec(
        "sfr", flopy.mf6.ModflowGwfsfr,
        metadata={"sfr_index": {
            "package": "sfr",
            "outlets": {"main": 4}, "heads": {"main": 3},
            "by_stream": {"main": [3, 4]},
            "cells": {"3": 0, "4": 0},
            "cell_rings": {"0": [[0.0, 0.0], [5.0, 0.0], [5.0, 5.0], [0.0, 5.0]]},
            "bounds": [0.0, 0.0, 5.0, 5.0],
        }},
    )
    assert sfr_connection(spec, "main", at=(2.5, 2.5)).index == 4
    # No "crs" in this hand-built index: the miss names the extent and stops
    # there rather than naming a frame it does not know.
    with pytest.raises(ValueError, match=r"to \(5.0, 5.0\) -- at= is a coordinate"):
        sfr_connection(spec, "main", at=(50.0, 50.0))


def test_the_index_carries_the_reach_cells_not_the_whole_grid():
    """A stream's index is proportional to the STREAM: the grid has 16 cells, not 4."""

    index = _sfr_mover_index(_builder(_grid(), FULL_ROW))
    assert sorted(int(cell) for cell in index["cell_rings"]) == [4, 5, 6, 7]
    assert index["cells"] == {"0": 4, "1": 5, "2": 6, "3": 7}
    assert index["bounds"] == [0.0, 0.0, 100.0, 100.0]
    assert index["crs"] == "EPSG:2927"


def test_the_index_survives_the_spec_round_trip():
    """``PackageSpec.to_dict`` writes JSON, where an integer key comes back a string.

    Ledger 179 promotes those back to ``int`` only for ``options``/``inputs``,
    deliberately never for ``metadata``, so the old int-keyed ``centroids`` lookup
    raised ``KeyError: 0`` on any spec that had been saved and reloaded -- through
    JSON or YAML alike, since ``_json_value`` stringifies the keys on write.
    """

    saved = _row_spec(FULL_ROW).to_dict()
    spec = PackageSpec.from_dict(saved)
    assert sfr_connection(spec, "s0", at=(47.0, 37.5)).index == 1
    assert sfr_connection(spec, "s0", at=(50.0, 37.5)).index == 2
    assert sfr_connection(spec, "s0", at=(62.5, 37.5)).index == 2

    through_yaml = PackageSpec.from_dict(yaml.safe_load(yaml.safe_dump(saved)))
    assert sfr_connection(through_yaml, "s0", at=(47.0, 37.5)).index == 1
    assert sfr_connection(through_yaml, "s0", at=(50.0, 37.5)).index == 2


def test_an_index_written_before_this_change_says_so():
    spec = PackageSpec(
        "sfr", flopy.mf6.ModflowGwfsfr,
        metadata={"sfr_index": {
            "package": "sfr", "outlets": {"main": 1}, "heads": {"main": 0},
            "by_stream": {"main": [0, 1]},
            "centroids": {0: [0.0, 0.0], 1: [5.0, 0.0]},
        }},
    )
    assert sfr_connection(spec, "main").index == 1
    with pytest.raises(ValueError, match="predates locating a point by cell"):
        sfr_connection(spec, "main", at=(1.0, 0.0))


def test_named_locations_and_unknown_streams_are_unchanged():
    spec = _row_spec(FULL_ROW)
    assert sfr_connection(spec, "s0").index == 3
    assert sfr_connection(spec, "s0", at="upstream").index == 0
    with pytest.raises(KeyError):
        sfr_connection(spec, "nope", at=(37.5, 37.5))
    with pytest.raises(ValueError, match="Unknown reach location"):
        sfr_connection(spec, "s0", at="sideways")
