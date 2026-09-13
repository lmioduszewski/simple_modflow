"""``edges_only=`` on a GeoPackage source is the edge of the ACTIVE domain.

``_edge_set`` called ``get_grid_edge()`` with no argument, which left it the
grid's own ``idomain_path`` and nothing else -- so a ``ModelContext`` domain, the
documented package-first route, never reached the edge test. Ledger 236 made a
cell ringing an inactive hole part of the edge; this caller could not see the
hole. Measured before the fix on a 5x5 grid with one inactive cell: 16 cells
written where the active edge is 20, and with the outer ring inactive the source
was refused outright rather than using the ring inside it.

Every expectation is row/column arithmetic on a regular grid, not a reading of
the code under test.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

from myflopy.geopackage import GeoPackageSource
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.specs import ModelContext

N = 5                 # 5 x 5 cells, row-major from the bottom left
D = 25.0              # each 25 ft square, so the grid spans 0..125
CRS = 2927

# Worked out by hand: row r holds cells 5r..5r+4.
OUTER = [0, 1, 2, 3, 4, 5, 9, 10, 14, 15, 19, 20, 21, 22, 23, 24]
CENTRE = 12
FACE_NEIGHBOURS_OF_CENTRE = [7, 11, 13, 17]     # not the corners 6, 8, 16, 18
INTERIOR = [6, 7, 8, 11, 12, 13, 16, 17, 18]


def _grid(idomain_path=None) -> VoronoiGridPlus:
    """N x N square cells, row-major from the bottom left."""

    verts = np.array([[i * D, j * D] for j in range(N + 1) for i in range(N + 1)], dtype=float)
    iverts = [
        [j * (N + 1) + i, (j + 1) * (N + 1) + i, (j + 1) * (N + 1) + i + 1, j * (N + 1) + i + 1]
        for j in range(N)
        for i in range(N)
    ]
    xcyc = [[(i + 0.5) * D, (j + 0.5) * D] for j in range(N) for i in range(N)]
    kwargs = {} if idomain_path is None else {"idomain_path": idomain_path}
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(xcyc), crs=CRS, **kwargs)


def _covering_source(tmp_path, vor, domain, *, layers=(1,), name="whole") -> GeoPackageSource:
    """A source whose feature(s) cover the WHOLE grid, one per entry in ``layers``."""

    whole = box(0, 0, N * D, N * D)
    gdf = gpd.GeoDataFrame(
        {
            "name": [f"{name}{i}" for i in range(len(layers))],
            "layer": list(layers),
            "head": [9.0] * len(layers),
        },
        geometry=[whole] * len(layers),
        crs=CRS,
    )
    path = tmp_path / "band.gpkg"
    gdf.to_file(path, layer="chd", driver="GPKG")
    return GeoPackageSource(
        path=path, context=ModelContext(grid=vor, domain=domain), nper=1, layer="chd"
    )


def _cells_written(source, **kwargs) -> list[int]:
    """The cell ids in the built package's first stress period."""

    spec = source.chd(head="head", edges_only=True, **kwargs)
    return sorted({cell for (_layer, cell), *_rest in spec.options["stress_period_data"][0]})


def _domain(inactive: list[int], nlay: int = 1) -> np.ndarray:
    """An MF6 idomain with ``inactive`` switched off in every layer."""

    domain = np.ones((nlay, N * N), dtype=int)
    for cell in inactive:
        domain[:, cell] = 0
    return domain


def test_a_cell_ringing_an_inactive_hole_is_on_the_edge(tmp_path):
    """The four face-neighbours of a one-cell hole join the outer ring."""

    vor = _grid()
    source = _covering_source(tmp_path, vor, _domain([CENTRE]))
    assert _cells_written(source) == sorted(OUTER + FACE_NEIGHBOURS_OF_CENTRE)


def test_the_corners_of_a_hole_are_not_on_the_edge(tmp_path):
    """A cell sharing only a VERTEX with the hole is still interior."""

    vor = _grid()
    source = _covering_source(tmp_path, vor, _domain([CENTRE]))
    written = _cells_written(source)
    assert [c for c in (6, 8, 16, 18) if c in written] == []


def test_an_inactive_outer_ring_leaves_the_ring_inside_it(tmp_path):
    """With the outer ring off, the edge is the interior's own rim -- not a refusal."""

    vor = _grid()
    source = _covering_source(tmp_path, vor, _domain(OUTER))
    # Eight of the nine interior cells border the inactive ring. The CENTRE does
    # not -- its four face-neighbours are all active -- so it is enclosed, and not
    # an edge cell. That asymmetry is the whole point of the topological rule.
    assert _cells_written(source) == [c for c in INTERIOR if c != CENTRE]


def test_the_edge_is_resolved_per_layer(tmp_path):
    """Two features, two layers, two different active domains -> two different edges."""

    vor = _grid()
    domain = np.ones((2, N * N), dtype=int)
    domain[0, CENTRE] = 0                      # layer 0 has the hole
    source = _covering_source(tmp_path, vor, domain, layers=(1, 2))
    spec = source.chd(head="head", edges_only=True)
    by_layer: dict[int, set[int]] = {}
    for (layer, cell), *_rest in spec.options["stress_period_data"][0]:
        by_layer.setdefault(layer, set()).add(cell)
    assert sorted(by_layer[0]) == sorted(OUTER + FACE_NEIGHBOURS_OF_CENTRE)
    assert sorted(by_layer[1]) == OUTER        # fully active: the outer ring only


def test_a_one_dimensional_domain_is_the_same_activity_in_every_layer(tmp_path):
    """A flat (ncpl,) domain reaches the edge test too."""

    vor = _grid()
    flat = np.ones(N * N, dtype=int)
    flat[CENTRE] = 0
    source = _covering_source(tmp_path, vor, flat)
    assert _cells_written(source) == sorted(OUTER + FACE_NEIGHBOURS_OF_CENTRE)


def test_a_passthrough_cell_is_a_hole_in_the_edge_too(tmp_path):
    """MF6's idomain is three-valued: a -1 cell holds no boundary, so it is a hole."""

    vor = _grid()
    domain = np.ones((1, N * N), dtype=int)
    domain[0, CENTRE] = -1
    source = _covering_source(tmp_path, vor, domain)
    assert _cells_written(source) == sorted(OUTER + FACE_NEIGHBOURS_OF_CENTRE)


def test_without_a_context_domain_the_grids_own_idomain_file_still_decides(tmp_path):
    """The pre-existing route -- a grid built with idomain_path -- keeps working."""

    idomain = tmp_path / "idomain.gpkg"
    x, y = (CENTRE % N) * D, (CENTRE // N) * D
    gpd.GeoDataFrame(
        {"name": ["hole"]}, geometry=[box(x, y, x + D, y + D)], crs=CRS
    ).to_file(idomain, driver="GPKG")
    vor = _grid(idomain_path=idomain)
    source = _covering_source(tmp_path, vor, None)
    assert _cells_written(source) == sorted(OUTER + FACE_NEIGHBOURS_OF_CENTRE)


def test_the_context_domain_wins_over_the_grids_idomain_file(tmp_path):
    """The model's declared domain is the domain; the grid's file does not override it."""

    idomain = tmp_path / "idomain.gpkg"
    gpd.GeoDataFrame(                          # the FILE says cell 6 is the hole
        {"name": ["hole"]}, geometry=[box(D, D, 2 * D, 2 * D)], crs=CRS
    ).to_file(idomain, driver="GPKG")
    vor = _grid(idomain_path=idomain)
    source = _covering_source(tmp_path, vor, _domain([CENTRE]))
    written = _cells_written(source)
    assert written == sorted(OUTER + FACE_NEIGHBOURS_OF_CENTRE)
    assert 11 in written and 7 in written       # the CONTEXT's hole ring
    assert 6 not in written                     # the file's hole is not consulted


def test_the_edge_is_computed_once_per_layer(tmp_path):
    """``get_grid_edge`` walks every cell, so a ten-feature source must not pay ten times."""

    vor = _grid()
    calls: list[object] = []
    real = vor.get_grid_edge

    def counted(*args, **kwargs):
        calls.append(kwargs.get("idomain"))
        return real(*args, **kwargs)

    vor.get_grid_edge = counted
    # four features, two of them sharing a layer
    source = _covering_source(tmp_path, vor, np.ones((2, N * N), dtype=int), layers=(1, 1, 2, 2))
    source.chd(head="head", edges_only=True)
    assert len(calls) == 2


def test_a_layer_past_the_domain_still_names_itself_under_edges_only(tmp_path):
    """The three-valued domain rule's own error is not swallowed by the edge lookup."""

    vor = _grid()
    source = _covering_source(tmp_path, vor, _domain([CENTRE], nlay=1), layers=(2,))
    with pytest.raises(ValueError, match="Layer 1 is outside the model domain"):
        source.chd(head="head", edges_only=True)


def test_edges_only_off_is_unaffected_by_the_edge_set(tmp_path):
    """Without ``edges_only`` every active cell the feature covers is written."""

    vor = _grid()
    source = _covering_source(tmp_path, vor, _domain([CENTRE]))
    spec = source.chd(head="head", edges_only=False)
    written = sorted({cell for (_l, cell), *_r in spec.options["stress_period_data"][0]})
    assert written == [c for c in range(N * N) if c != CENTRE]
