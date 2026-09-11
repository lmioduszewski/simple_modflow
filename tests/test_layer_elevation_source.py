"""Three surface-aware builders, one elevation table.

SFR read ``grid.gdf_topbtm`` directly and never looked at ``context.surfaces``.
LAK and ``CellSurfaceOffset`` read ``context.surfaces`` first. Because
``ModelContext(surfaces=vor.gdf_topbtm)`` captures the frame OBJECT, re-running a
layer-build cell in a notebook -- which re-publishes ``vor.gdf_topbtm`` -- left
two packages in one model on different elevations: measured SFR reach top 70.0
against ``CellSurfaceOffset`` 60.0, truth 30.0, no exception, normal termination.
That was the only route found to wrong numbers in an MF6 input file.

They now share ``specs._layer_elevations``. ``context.surfaces`` still wins when
given -- a deliberate pin and a stale snapshot look identical -- but every
builder agrees on it, and a disagreement with the grid's live table is logged.
"""

from __future__ import annotations

import logging

import geopandas as gpd
import numpy as np
import pytest

import myflopy as mf
from myflopy.specs import ModelContext, _layer_elevations
from tests.test_sfr_builder import _builder, _grid


def _shifted_table(grid, shift):
    """A table shaped like ``grid.gdf_topbtm`` with every elevation moved."""

    live = grid.gdf_topbtm
    return gpd.GeoDataFrame(
        {0: np.asarray(live[0], dtype=float) + shift, 1: np.asarray(live[1], dtype=float)},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )


# --- the resolver ----------------------------------------------------------- #


def test_with_no_context_surfaces_the_grids_live_table_is_used():
    grid = _grid()
    assert _layer_elevations(ModelContext(grid=grid)) is grid.gdf_topbtm


def test_context_surfaces_wins_when_given():
    grid = _grid()
    pinned = _shifted_table(grid, 100.0)
    assert _layer_elevations(ModelContext(grid=grid, surfaces=pinned)) is pinned


def test_a_disagreement_with_the_live_table_is_announced(caplog):
    grid = _grid()
    context = ModelContext(grid=grid, surfaces=_shifted_table(grid, 100.0))
    with caplog.at_level(logging.WARNING, logger="myflopy.specs"):
        _layer_elevations(context)
    assert "disagree" in caplog.text


def test_the_same_numbers_in_a_different_object_are_not_a_disagreement(caplog):
    """A copy is not a stale snapshot; warning on identity would cry wolf."""

    grid = _grid()
    context = ModelContext(grid=grid, surfaces=grid.gdf_topbtm.copy())
    with caplog.at_level(logging.WARNING, logger="myflopy.specs"):
        _layer_elevations(context)
    assert caplog.text == ""


def test_a_layer_build_result_is_turned_into_the_table():
    """LAK used to die on this: 'LayerBuildResult' object has no attribute 'drop'."""

    grid = _grid()
    layers = (
        mf.LayerStack(grid, top=mf.Surface.flat(50.0))
        .add("sand", thickness=10.0)
        .add("clay", thickness=10.0)
        .build(grid)
    )
    table = _layer_elevations(ModelContext(grid=grid, surfaces=layers))
    assert list(table.columns) == [0, 1, 2]
    assert table[0].tolist() == pytest.approx([50.0] * grid.ncpl)


def test_resolving_a_layer_build_result_does_not_publish_it_onto_the_grid():
    grid = _grid()
    before = grid.gdf_topbtm
    layers = (
        mf.LayerStack(grid, top=mf.Surface.flat(50.0)).add("sand", thickness=10.0).build(grid)
    )
    _layer_elevations(ModelContext(grid=grid, surfaces=layers))
    assert grid.gdf_topbtm is before


# --- the bug, through SFR's public output ----------------------------------- #


def test_sfr_reach_tops_come_from_context_surfaces_like_every_other_builder():
    """SFR used to read the grid here while LAK read the context."""

    grid = _grid()
    pinned = _shifted_table(grid, 100.0)
    builder = _builder(
        context=ModelContext(grid=grid, surfaces=pinned, domain=np.ones((1, 6), dtype=int))
    )
    live_tops = np.asarray(grid.gdf_topbtm[0], dtype=float)
    tops = builder._default_reach_tops()

    assert min(tops) >= 100.0            # read from the pinned table
    assert max(tops) > live_tops.max()   # and not from the live one


def test_sfr_still_reads_the_grid_when_no_surfaces_are_given():
    grid = _grid()
    builder = _builder(context=ModelContext(grid=grid, domain=np.ones((1, 6), dtype=int)))
    tops = builder._default_reach_tops()
    assert max(tops) <= float(np.asarray(grid.gdf_topbtm[0], dtype=float).max())


def test_sfr_reads_a_unit_labelled_table_by_position():
    """Columns are [top, botm0, ...] by position; labels must not matter."""

    grid = _grid()
    live = grid.gdf_topbtm
    labelled = gpd.GeoDataFrame(
        {"top": np.asarray(live[0], dtype=float), "sand": np.asarray(live[1], dtype=float)},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    by_label = _builder(
        context=ModelContext(grid=grid, surfaces=labelled, domain=np.ones((1, 6), dtype=int))
    )._default_reach_tops()
    by_integer = _builder(
        context=ModelContext(grid=grid, domain=np.ones((1, 6), dtype=int))
    )._default_reach_tops()
    assert by_label == pytest.approx(by_integer)
