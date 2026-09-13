"""A PRT release point on a cell face goes in the lowest ACTIVE cell, not the lowest cell.

``PRTReleasePoints.from_points`` located each point with FloPy's
``modelgrid.intersect``, which answers with the lowest CELL2D and knows nothing
about idomain. With cells 0 and 4 inactive in layer 0, a point on the face
between cells 4 and 5 was released into cell 4 and a point on the vertex shared
by 0, 1, 4 and 5 into cell 0 -- both inactive, both in silence. A point off the
grid entirely got whatever cell ``intersect`` returned for it.

The coordinates WRITTEN are a separate question, settled by ledger 218: MF6 reads
a PRP record in the model's own frame, so an offset grid's release points are
converted with ``get_local_coords``. That conversion is unchanged here and pinned
below on an offset grid of its own.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pytest
from flopy.discretization import VertexGrid
from shapely.geometry import Point

from myflopy.modflow.mf6 import prt as prt_module
from myflopy.modflow.mf6.grid import domain as grid_domain
from myflopy.modflow.mf6.grid import geometry as grid_geometry
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.prt import PRTReleasePoints

LOGGER = "myflopy.modflow.mf6.prt"
N, LENGTH = 4, 100.0
D = LENGTH / N


def _mesh(x0: float = 0.0, y0: float = 0.0):
    """4x4 square cells, row-major from the bottom left (cell 5 is 25..50, 25..50)."""

    verts = np.array(
        [[x0 + i * D, y0 + j * D] for j in range(N + 1) for i in range(N + 1)], dtype=float
    )
    iverts = [
        [j * (N + 1) + i, (j + 1) * (N + 1) + i, (j + 1) * (N + 1) + i + 1, j * (N + 1) + i + 1]
        for j in range(N)
        for i in range(N)
    ]
    xcyc = np.array(
        [[x0 + (i + 0.5) * D, y0 + (j + 0.5) * D] for j in range(N) for i in range(N)],
        dtype=float,
    )
    return verts, iverts, xcyc


def _model(idomain=None, *, xoff: float = 0.0, yoff: float = 0.0):
    """A model-like object: a world-frame ``vor``, and a DISV whose mesh is local.

    ``xoff``/``yoff`` reproduce a USG import (``to_mf6(local_origin=True)``): the
    ``vor`` carries the world coordinates the .gsf held, the DISV the shifted mesh
    plus the origin it was shifted by.
    """

    verts, iverts, xcyc = _mesh(xoff, yoff)
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    local_verts, local_cells, _ = _mesh()
    nlay = 1 if idomain is None else np.asarray(idomain).shape[0]
    modelgrid = VertexGrid(
        vertices=[[i, float(x), float(y)] for i, (x, y) in enumerate(local_verts)],
        cell2d=[
            [cell, float(x), float(y), len(iv), *iv]
            for cell, ((x, y), iv) in enumerate(zip(_mesh()[2], local_cells, strict=True))
        ],
        nlay=nlay,
        top=np.full(N * N, 10.0),
        botm=np.zeros((nlay, N * N)),
        idomain=None if idomain is None else np.asarray(idomain, dtype=int),
        xoff=xoff,
        yoff=yoff,
    )
    disv = None if idomain is None else SimpleNamespace(
        idomain=SimpleNamespace(array=np.asarray(idomain, dtype=int))
    )
    return SimpleNamespace(vor=vor, gwf=SimpleNamespace(disv=disv, modelgrid=modelgrid))


def _domain(inactive: dict[int, list[int]], nlay: int = 1):
    """An idomain of ones with ``{layer: [cells]}`` switched off."""

    values = np.ones((nlay, N * N), dtype=int)
    for layer, cells in inactive.items():
        for cell in cells:
            values[layer, cell] = 0
    return values


# -- the two measured cases ------------------------------------------------


def test_a_point_on_a_face_takes_the_lowest_ACTIVE_cell():
    """(25, 37.5) is on the face between cells 4 and 5; cell 4 is inactive."""

    model = _model(_domain({0: [0, 4]}))

    points = PRTReleasePoints.from_points(model, [(25.0, 37.5)])

    assert [row[1] for row in points.packagedata] == [(0, 5)]


def test_a_point_on_a_vertex_takes_the_lowest_active_cell():
    """(25, 25) is the vertex of cells 0, 1, 4 and 5; 0 and 4 are inactive."""

    model = _model(_domain({0: [0, 4]}))

    points = PRTReleasePoints.from_points(model, [(25.0, 25.0)])

    assert [row[1] for row in points.packagedata] == [(0, 1)]


# -- the rest of the rule --------------------------------------------------


def test_a_point_inside_a_cell_is_that_cell_and_keeps_its_coordinates():
    """``local_z`` is the caller's, not the default -- the record is decided here now."""

    model = _model(_domain({0: [0, 4]}))

    (row,) = PRTReleasePoints.from_points(
        model, [Point(37.5, 37.5)], local_z=0.75
    ).packagedata

    assert row[1] == (0, 5)
    assert (row[2], row[3]) == pytest.approx((37.5, 37.5))
    assert row[4] == 0.75
    assert PRTReleasePoints.from_points(model, [Point(37.5, 37.5)]).packagedata[0][4] == 0.5


def test_the_lowest_cell_wins_the_tie_when_every_cell_is_active():
    """With no domain the answer is still ONE cell, and the lowest of the four."""

    assert PRTReleasePoints.from_points(_model(), [(25.0, 25.0)]).packagedata[0][1] == (0, 0)


def test_the_release_layer_decides_which_activity_is_read():
    """Cell 5 is inactive in layer 0 and active in layer 1, so the layer decides."""

    model = _model(_domain({0: [5]}, nlay=2))

    assert PRTReleasePoints.from_points(model, [(37.5, 37.5)], layer=1).packagedata[0][1] == (1, 5)
    with pytest.raises(ValueError, match="on no active cell in layer 0"):
        PRTReleasePoints.from_points(model, [(37.5, 37.5)], layer=0)


def test_a_release_layer_beyond_the_domain_raises():
    """A domain of one layer cannot answer for layer 3; it used to write ``(3, 5)``."""

    with pytest.raises(ValueError, match="outside the model domain"):
        PRTReleasePoints.from_points(_model(_domain({})), [(37.5, 37.5)], layer=3)


def test_a_passthrough_cell_is_not_an_active_cell():
    """idomain -1 is a vertical passthrough, and ``bool(-1)`` is True."""

    domain = np.ones((1, N * N), dtype=int)
    domain[0, 4] = -1

    points = PRTReleasePoints.from_points(_model(domain), [(25.0, 37.5)])

    assert [row[1] for row in points.packagedata] == [(0, 5)]


def test_every_point_off_the_grid_raises_naming_both_extents():
    """It used to take whatever FloPy's intersect answered with (ledger 211)."""

    model = _model(_domain({0: [0, 4]}))

    with pytest.raises(ValueError, match="No release point lands on an active cell"):
        PRTReleasePoints.from_points(model, [(1_300_000.0, 650_000.0)])
    with pytest.raises(ValueError, match="The release points span"):
        PRTReleasePoints.from_points(model, [(1_300_000.0, 650_000.0)])


def test_one_point_off_the_grid_warns_by_name_and_renumbers_the_rest(caplog):
    """MF6 wants release points 1..NRELEASEPTS, so the survivors are renumbered."""

    model = _model(_domain({0: [0, 4]}))

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        points = PRTReleasePoints.from_points(
            model, [(37.5, 37.5), (500.0, 500.0), (62.5, 62.5)]
        )

    assert [row[0] for row in points.packagedata] == [0, 1]
    assert [row[1] for row in points.packagedata] == [(0, 5), (0, 10)]
    assert any("#1 (500.0, 500.0)" in record.getMessage() for record in caplog.records)


def test_a_dropped_point_does_not_shift_the_group_labels(caplog):
    """The labels are indexed by the point GIVEN, not by the record written."""

    model = _model(_domain({0: [0, 4]}))

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        points = PRTReleasePoints.from_points(
            model,
            [(37.5, 37.5), (500.0, 500.0), (62.5, 62.5)],
            group=["west", "nowhere", "east"],
        )

    assert points.groups == ("west", "east")


def test_points_are_reprojected_to_the_grid_before_they_are_located():
    """Raw degrees over a State Plane grid reach no cell at all."""

    model = _model(_domain({0: [0, 4]}))
    lonlat = gpd.GeoSeries([Point(37.5, 37.5)], crs=model.vor.crs).to_crs("EPSG:4326")

    (row,) = PRTReleasePoints.from_points(model, lonlat).packagedata

    assert row[1] == (0, 5)
    assert (row[2], row[3]) == pytest.approx((37.5, 37.5), abs=1e-3)


def test_a_geodataframe_of_points_is_read_by_its_geometry():
    """The realistic input is a well file, and it arrives as a GeoDataFrame.

    Read as a plain iterable a GeoDataFrame yields its COLUMN NAMES, so this is
    what holds the ``geometry`` unwrap in place -- and with it the only branch
    that can carry a CRS to reproject.
    """

    model = _model(_domain({0: [0, 4]}))
    wells = gpd.GeoDataFrame(
        {"name": ["mw-1"]}, geometry=[Point(25.0, 37.5)], crs=model.vor.crs
    )

    (row,) = PRTReleasePoints.from_points(model, wells).packagedata

    assert row[1] == (0, 5)


def test_a_miss_is_diagnosed_where_the_point_LANDED_not_where_it_was_given():
    """The reason is the grid's frame's; the coordinates named are the caller's.

    A point reprojected onto an inactive cell has missed for a reason the user can
    act on -- "on no active cell" -- and measuring it in the source frame instead
    would report "outside the grid" of a point that is inside it, which is the one
    distinction ledger 211's wording exists to make.
    """

    model = _model(_domain({0: [0, 4]}))
    lonlat = gpd.GeoSeries([Point(12.5, 12.5)], crs=model.vor.crs).to_crs("EPSG:4326")

    with pytest.raises(ValueError) as raised:
        PRTReleasePoints.from_points(model, lonlat)

    assert "(on no active cell in layer 0)" in str(raised.value)
    assert "outside the grid" not in str(raised.value)


def test_more_than_ten_misses_are_listed_to_ten_and_then_counted(caplog):
    """Twelve bad rows in a well file are not twelve lines of coordinates."""

    model = _model(_domain({0: [0, 4]}))
    points = [(37.5, 37.5)] + [(500.0 + 10.0 * i, 500.0) for i in range(12)]

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        released = PRTReleasePoints.from_points(model, points)

    (message,) = [record.getMessage() for record in caplog.records]
    assert len(released.packagedata) == 1
    assert "#10 (590.0, 500.0)" in message
    assert "and 2 more" in message
    assert "#12" not in message


# -- the frame the record is written in (ledger 218) ------------------------


def test_release_points_stay_in_the_models_own_frame_on_an_offset_grid():
    """A point inside cell 5, so the only thing this can fail on is the frame."""

    x0, y0 = 1_340_000.0, 650_000.0
    model = _model(_domain({0: [0, 4]}), xoff=x0, yoff=y0)

    (row,) = PRTReleasePoints.from_points(model, [(x0 + 37.5, y0 + 37.5)]).packagedata

    assert row[1] == (0, 5)
    assert (row[2], row[3]) == pytest.approx((37.5, 37.5))


def test_a_face_on_an_offset_grid_is_resolved_in_the_world_frame():
    """The cell is picked from the world point; cells 0 and 4 are inactive."""

    x0, y0 = 1_340_000.0, 650_000.0
    model = _model(_domain({0: [0, 4]}), xoff=x0, yoff=y0)

    (row,) = PRTReleasePoints.from_points(model, [(x0 + 25.0, y0 + 37.5)]).packagedata

    assert row[1] == (0, 5)
    assert (row[2], row[3]) == pytest.approx((25.0, 37.5))


def test_an_offset_grids_own_coordinates_reach_no_cell():
    """The frame is not guessed: (25, 37.5) on a grid at State Plane is off it."""

    model = _model(_domain({0: [0, 4]}), xoff=1_340_000.0, yoff=650_000.0)

    with pytest.raises(ValueError, match="outside the grid"):
        PRTReleasePoints.from_points(model, [(25.0, 37.5)])


# -- one copy of the rule --------------------------------------------------


def test_the_rule_is_the_shared_one_not_a_copy():
    assert prt_module.cell_shares is grid_geometry.cell_shares
    assert prt_module.domain_active is grid_geometry.domain_active
    assert prt_module.miss_reason is grid_geometry.miss_reason
    assert prt_module.where_they_are is grid_geometry.where_they_are
    assert prt_module.model_domain is grid_domain.model_domain


def _structured_model():
    """The same 4x4 mesh as a DIS grid, with a ``vor`` the shared rule could have used."""

    from flopy.discretization import StructuredGrid

    verts, iverts, xcyc = _mesh()
    grid = StructuredGrid(
        delr=np.full(N, D), delc=np.full(N, D), top=np.full((N, N), 10.0),
        botm=np.zeros((1, N, N)), nlay=1,
    )
    return SimpleNamespace(
        vor=VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc),
        gwf=SimpleNamespace(dis=None, modelgrid=grid),
    )


def test_a_structured_grid_keeps_flopys_own_intersect():
    """The shared rule is indexed by cell id; a DIS cell is ``(row, col)``.

    The same 4x4 mesh, so the cell is the same PLACE either way and the only
    difference is the cellid's shape: FloPy's ``(0, 2, 1)`` against the DISV
    rule's ``(0, 5)``. The grid KIND is what decides, which is why the ``vor``
    here carries a cell frame the rule could otherwise have used.
    """

    (row,) = PRTReleasePoints.from_points(_structured_model(), [(37.5, 37.5)]).packagedata

    assert row[1] == (0, 2, 1)


def test_a_structured_grid_keeps_flopys_refusal_as_well_as_its_answer():
    """Ledger 211 is DISV-only here: on DIS, FloPy's own raise stands, unimproved.

    It names neither extent nor either CRS, and it ends the whole call over one
    bad row -- deliberately left, because the activity half of the rule cannot
    reach a ``(row, col)`` cell either. Pinned so nobody changes half of it.
    """

    with pytest.raises(ValueError, match="outside of the model area"):
        PRTReleasePoints.from_points(_structured_model(), [(37.5, 37.5), (500.0, 500.0)])
