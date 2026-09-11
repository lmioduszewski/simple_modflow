"""A polygon drawn on cell faces covers the cells inside it, and no others.

``GeoPackageSource._cell_shares`` measured each piece of a feature by the piece's
own extent: area if it had any, otherwise LENGTH. A polygon meets each
face-adjacent cell along the edge they share -- a line, with no area -- so every
neighbour of a polygon drawn on cell faces was admitted, with a pseudo-share of
edge / area (L / L**2). A field not wrapped in ``mf.Spread`` is written whole to
each admitted cell, so a DRN polygon that is exactly one cell became five records
at five times the conductance, and a K zone painted five cells. MF6 terminates
normally either way.

Exact coincidence is common, not contrived: square-cell DISV grids with round
spacing, polygons dissolved or exported from the grid itself, nested refinement.
Its near miss is commoner still -- the same polygon after a CRS round trip, or
read back from a DISV FloPy wrote at 9 significant digits -- so a piece narrower
than ``_SLIVER_WIDTH`` is dropped as well.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import shapely
from shapely.geometry import GeometryCollection, LineString, MultiPolygon, box

import myflopy as mf
from myflopy.geopackage import GeoPackageSource, Spread
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.specs import ModelContext

LIST_BCS = ("chd", "ghb", "drn", "riv", "wel", "rch", "evt")

# Every registry default column, so one file serves all seven helpers.
_VALUES = {
    "head": 10.0,
    "conductance": 100.0,
    "elevation": 10.0,
    "stage": 10.0,
    "rbot": 5.0,
    "rate": -500.0,
    "recharge": 0.01,
    "surface": 10.0,
    "depth": 2.0,
    "k": 7.0,
}

# Exactly cell 5 of the 4x4 grid below; its face-neighbours are 1, 4, 6 and 9.
CELL_5 = box(25.0, 25.0, 50.0, 50.0)


def _context(n: int = 4, length: float = 100.0, x0: float = 0.0, y0: float = 0.0):
    """A real ``VoronoiGridPlus`` of n x n square cells, row-major from the bottom left."""

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
    return ModelContext(grid=grid, domain=np.ones((1, n * n), dtype=int))


def _gpkg(tmp_path, geometry):
    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {"name": ["feature"], "layer": [1], **{k: [v] for k, v in _VALUES.items()}},
        geometry=[geometry],
        crs="EPSG:2927",
    ).to_file(path, layer="bc", driver="GPKG")
    return path


def _records(spec):
    return spec.options["stress_period_data"][0]


def _cells(spec):
    return sorted(record[0][1] for record in _records(spec))


# --- exact coincidence ------------------------------------------------------ #


@pytest.mark.parametrize("package", LIST_BCS)
def test_a_polygon_that_is_one_cell_resolves_to_that_cell_in_every_list_bc(tmp_path, package):
    spec = getattr(mf, package).gpkg(
        _gpkg(tmp_path, CELL_5), layer="bc", context=_context(), nper=1
    )
    assert _cells(spec) == [5]


@pytest.mark.parametrize(
    "conductance",
    [
        "conductance",
        Spread("conductance"),
        Spread("conductance", mode="retained"),
        Spread("conductance", min_share=0.0),
    ],
    ids=["whole", "clip", "retained", "no-floor"],
)
def test_the_whole_value_lands_in_the_one_cell_however_it_is_split(tmp_path, conductance):
    """The pseudo-share was L / L**2, so whether min_share hid it depended on units."""

    records = _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, CELL_5), layer="bc", context=_context(), nper=1,
            conductance=conductance,
        )
    )
    assert [record[0][1] for record in records] == [5]
    assert records[0][2] == pytest.approx(100.0)


@pytest.mark.parametrize("n", [4, 8])
def test_refining_under_a_face_aligned_polygon_keeps_its_cell_count_and_total(tmp_path, n):
    """On fine faces but across coarse cells: 4 cells either way, never 12."""

    records = _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, box(37.5, 37.5, 62.5, 62.5)), layer="bc", context=_context(n),
            nper=1, conductance=Spread("conductance"),
        )
    )
    assert len(records) == 4
    assert sum(record[2] for record in records) == pytest.approx(100.0)


def test_retained_mode_keeps_the_total_inside_the_polygon(tmp_path):
    """Renormalizing the pseudo-shares conserved the total and put 16% of it outside."""

    records = _records(
        mf.drn.gpkg(
            _gpkg(tmp_path, box(37.5, 37.5, 62.5, 62.5)), layer="bc", context=_context(8),
            nper=1, conductance=Spread("conductance", mode="retained"),
        )
    )
    assert {record[0][1]: record[2] for record in records} == pytest.approx(
        {27: 25.0, 28: 25.0, 35: 25.0, 36: 25.0}
    )


def test_a_k_zone_on_cell_faces_paints_only_its_own_cell(tmp_path):
    source = GeoPackageSource(_gpkg(tmp_path, CELL_5), _context(), 1, layer="bc")
    k = source.k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [5]


# --- near coincidence, and what must NOT be mistaken for it ----------------- #
#
# A perturbation is a WIDTH, not a share of the cell, so the rule is on a piece's
# mean width (2 * area / perimeter). A share rule caught a NAD83 round trip on 500 ft
# cells and missed it on 300 ft ones -- this file only tested 500 ft, and missed it.

# Non-round, and away from the origin, where the NAD83/HARN shift vanishes.
STATE_PLANE = dict(x0=1_300_000.123456, y0=650_000.654321)


def _slivered(moved, polys):
    return [c for c in range(len(polys)) if c != 5 and moved.intersection(polys.loc[c]).area > 0]


@pytest.mark.parametrize("cell_size", [500.0, 100.0, 25.0])
@pytest.mark.parametrize("via", ["EPSG:4326", "EPSG:32610", "EPSG:2285"])
def test_the_same_polygon_after_a_crs_round_trip_still_resolves_to_one_cell(
    tmp_path, via, cell_size
):
    context = _context(length=4 * cell_size, **STATE_PLANE)
    polys = context.grid.gdf_vorPolys.geometry
    moved = gpd.GeoSeries([polys.loc[5]], crs="EPSG:2927").to_crs(via).to_crs("EPSG:2927").iloc[0]
    assert _slivered(moved, polys), "the round trip must leave a sliver, or this tests nothing"

    spec = mf.drn.gpkg(_gpkg(tmp_path, moved), layer="bc", context=context, nper=1)
    assert _cells(spec) == [5]


@pytest.mark.parametrize("cell_size", [100.0, 25.0])
def test_a_polygon_read_back_from_a_written_disv_still_resolves_to_one_cell(tmp_path, cell_size):
    """FloPy writes DISV vertices at 9 significant digits -- 0.01 ft at State Plane --
    so a zone drawn from the WRITTEN cells is off the in-memory mesh by up to 0.005."""

    context = _context(length=4 * cell_size, **STATE_PLANE)
    polys = context.grid.gdf_vorPolys.geometry
    written = shapely.transform(
        polys.loc[5], lambda xy: np.array([[float(f"{v:.8E}") for v in point] for point in xy])
    )
    assert _slivered(written, polys), "rounding must leave a sliver, or this tests nothing"

    source = GeoPackageSource(_gpkg(tmp_path, written), context, 1, layer="bc")
    k = source.k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [5]


def test_a_feature_thinner_than_the_tolerance_is_kept_whole(tmp_path):
    """A pond 0.01 ft across is what was drawn, not a sliver of something else."""

    spec = mf.drn.gpkg(
        _gpkg(tmp_path, box(760.0, 760.0, 760.01, 760.01)), layer="bc",
        context=_context(length=2000.0), nper=1,
    )
    assert _cells(spec) == [5]


def test_a_small_part_of_a_multipolygon_is_kept_like_the_same_part_alone(tmp_path):
    """Dissolving a pond into its zone must not lose it: 0.16 ft2 of a 500 ft cell."""

    zone = MultiPolygon([box(0.0, 0.0, 1000.0, 1000.0), box(1760.0, 1760.0, 1760.4, 1760.4)])
    source = GeoPackageSource(_gpkg(tmp_path, zone), _context(length=2000.0), 1, layer="bc")
    k = source.k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [0, 1, 4, 5, 15]


def test_a_real_overlap_of_a_neighbour_still_reaches_it(tmp_path):
    """A strip 0.015 ft wide into cell 6 is past the tolerance: coverage, not noise."""

    spec = mf.drn.gpkg(
        _gpkg(tmp_path, box(25.0, 25.0, 50.015, 50.0)), layer="bc", context=_context(), nper=1
    )
    assert _cells(spec) == [5, 6]


def test_a_collection_mixing_polygons_and_lines_is_refused(tmp_path):
    """It has no single measure. Measured as a polygon, its line was dropped silently."""

    mixed = GeometryCollection(
        [box(30.0, 30.0, 45.0, 45.0), LineString([(60.0, 10.0), (90.0, 10.0)])]
    )
    with pytest.raises(ValueError, match="one feature per geometry type"):
        mf.chd.gpkg(_gpkg(tmp_path, mixed), layer="bc", context=_context(), nper=1)
