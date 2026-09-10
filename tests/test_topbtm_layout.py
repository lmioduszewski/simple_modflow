"""``gdf_topbtm`` is read positionally, so the geometry column's PLACE is a contract.

Several consumers decode the layer-elevation frame by position -- ``iloc[:, 1]``
for the model top, ``iloc[:, 2:]`` for the bottoms -- and nothing enforced where
the geometry column sat. ``LayerSurfaces.sample(reconcile=True)`` (the default)
dropped geometry, reconciled, and re-wrapped with ``gpd.GeoDataFrame(df,
geometry=...)``, which appends geometry at the END and undid the explicit
geometry-first ordering three lines above it.

Measured on a three-layer frame before the fix: the layer-elevation hover read
the first bottom as the model top (5.0 where the truth was 10.0) and handed a
Polygon back as the last bottom, raising nothing. Integer column labels do not
protect against it, which is why unifying the column NAMES would not have helped.

Normalization lives on the GETTER, not the setter, because two routes bypass the
setter: ``__setstate__`` writes ``__dict__`` directly, and an in-place edit
mutates the frame the getter handed back. The getter stores what it normalized
and returns that same object -- returning a normalized copy would make an
in-place edit write into a temporary and vanish.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.grid.geometry import geometry_first
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.surfaces import LayerSurfaces, Surface


def _geometry_last(n: int = 3) -> gpd.GeoDataFrame:
    """A three-surface frame in the shape ``reconcile=True`` used to produce."""

    frame = gpd.GeoDataFrame(
        {
            0: [10.0] * n,
            1: [5.0] * n,
            2: [0.0] * n,
            "geometry": [Point(float(i), 0.0) for i in range(n)],
        },
        geometry="geometry",
        crs="EPSG:2927",
    )
    assert list(frame.columns)[-1] == "geometry"  # the shape under test
    return frame


def _small_grid() -> VoronoiGridPlus:
    """A 6-cell grid built without triangle, cheap enough for a unit test."""

    verts = np.array(
        [
            [0, 0], [1, 0], [2, 0], [3, 0],
            [0, 1], [1, 1], [2, 1], [3, 1],
            [0, 2], [1, 2], [2, 2], [3, 2],
        ],
        dtype=float,
    )
    iverts = [
        [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6],
        [4, 5, 9, 8], [5, 6, 10, 9], [6, 7, 11, 10],
    ]
    centers = np.array(
        [[0.5, 0.5], [1.5, 0.5], [2.5, 0.5], [0.5, 1.5], [1.5, 1.5], [2.5, 1.5]]
    )
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=centers)


# --- the helper ------------------------------------------------------------ #


def test_geometry_first_moves_geometry_to_the_front():
    assert list(geometry_first(_geometry_last()).columns) == ["geometry", 0, 1, 2]


def test_geometry_first_is_idempotent():
    once = geometry_first(_geometry_last())
    assert list(geometry_first(once).columns) == list(once.columns)


def test_an_already_ordered_frame_is_returned_unchanged():
    """Not merely equal -- the SAME object, so a read costs no copy."""

    ordered = geometry_first(_geometry_last())
    assert geometry_first(ordered) is ordered


def test_geometry_is_found_by_dtype_not_by_the_name_geometry():
    """A frame whose geometry column is named `geom` was swept in as a bottom."""

    frame = gpd.GeoDataFrame(
        {0: [10.0], 1: [5.0], "geom": [Point(0.0, 0.0)]},
        geometry="geom",
        crs="EPSG:2927",
    )
    assert list(geometry_first(frame).columns) == ["geom", 0, 1]


def test_geometry_first_passes_through_none_and_non_frames():
    assert geometry_first(None) is None
    payload = {"top": 1.0}
    assert geometry_first(payload) is payload


# --- what the bug actually did --------------------------------------------- #


def test_a_positional_reader_gets_the_right_elevations():
    """The exact decode `choros.py` performs, against the known truth."""

    frame = geometry_first(_geometry_last())
    top = frame.iloc[:, 1].to_numpy()
    botm = frame.iloc[:, 2:].to_numpy()

    assert top.tolist() == [10.0, 10.0, 10.0]
    assert botm[0].tolist() == [5.0, 0.0]
    # Before the fix this read top=5.0 and a Polygon as the last bottom.
    assert all(isinstance(value, float) for value in botm[0])


# --- the writer sites ------------------------------------------------------ #


@pytest.mark.parametrize("reconcile", [True, False])
def test_sample_puts_geometry_first_under_either_reconcile(reconcile):
    """`reconcile=` decides physics; it must not decide the column order."""

    from tests.test_surfaces import _FakeVor

    gdf = LayerSurfaces([Surface.flat(50), Surface.flat(40)]).sample(
        _FakeVor(3), reconcile=reconcile
    )
    assert list(gdf.columns)[0] == "geometry"


def test_layer_stack_to_disv_publishes_integer_columns():
    """`to_disv` is the path the cheatsheet teaches; it published string labels.

    Those made `mf.sfr` and `vor.plot.section` raise `KeyError: 0` on a grid the
    user had just built, so the frame must carry this build's integer columns.
    """

    grid = _small_grid()
    stack = (
        mf.LayerStack(grid, top=Surface.flat(50.0))
        .add("sand", thickness=10.0)
        .add("clay", thickness=10.0)
    )
    stack.to_disv(grid)

    columns = list(grid.gdf_topbtm.columns)
    assert columns[0] == "geometry"
    assert columns[1:] == [0, 1, 2]
    assert grid.nlay == 2


# --- the getter, which is where the enforcement lives ---------------------- #


def test_the_getter_normalizes_a_frame_written_straight_into_the_instance():
    """`__setstate__` does `__dict__.update`, so a reloaded grid skips the setter."""

    grid = _small_grid()
    grid.__dict__["_gdf_topbtm"] = _geometry_last(grid.ncpl)
    assert list(grid.gdf_topbtm.columns)[0] == "geometry"


def test_an_in_place_edit_survives_the_normalizing_getter():
    """The trap in the obvious fix.

    A getter that returned a normalized COPY would hand this assignment a
    temporary and lose it silently -- measured, a bottom set to 99.0 read back as
    5.0. That would be a worse bug than the one being fixed, and the canonical
    notebook edits the frame exactly this way.
    """

    grid = _small_grid()
    grid.__dict__["_gdf_topbtm"] = _geometry_last(grid.ncpl)

    grid.gdf_topbtm[1] = [99.0] * grid.ncpl

    assert grid.gdf_topbtm[1].tolist() == [99.0] * grid.ncpl


def test_the_getter_returns_the_stored_object_not_a_copy():
    grid = _small_grid()
    grid.__dict__["_gdf_topbtm"] = _geometry_last(grid.ncpl)
    assert grid.gdf_topbtm is grid.gdf_topbtm


def test_the_setter_stores_an_already_normalized_frame():
    """Redundant with the getter by behaviour, deliberate, and pinned here.

    The getter is the enforcement point -- it has to be, since two routes bypass
    the setter. But normalizing on write too means anything reading the private
    ``_gdf_topbtm`` sees an ordered frame without a property access, so this
    asserts against the storage rather than through the getter. Mutation-tested:
    without it every behavioural test still passes.
    """

    grid = _small_grid()
    grid.gdf_topbtm = _geometry_last(grid.ncpl)
    assert list(grid.__dict__["_gdf_topbtm"].columns)[0] == "geometry"


def test_reattaching_a_frame_drops_the_cached_layer_count():
    """`_nlay` cached on first read and nothing invalidated it."""

    grid = _small_grid()
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [12.0] * 6, 1: [8.0] * 6, 2: [0.0] * 6},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    assert grid.nlay == 2

    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [12.0] * 6, 1: [8.0] * 6},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    assert grid.nlay == 1
