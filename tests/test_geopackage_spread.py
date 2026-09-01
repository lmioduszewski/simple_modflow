"""``mf.Spread``: splitting an extensive GeoPackage field across the cells it covers.

The behaviour being pinned is that a conductance belongs to the *feature*. Without
``Spread`` a feature's value is written to every cell it intersects, which is right
for an elevation and multiplies a conductance by the cell count -- on real data, by
2.02 on the very grid the values came from.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point, Polygon

import myflopy as mf
from myflopy.geopackage import Spread, _extent


@pytest.fixture
def grid_and_context(tmp_path):
    """A 4x1 strip of unit square cells, and a context over it."""

    from myflopy.specs import ModelContext

    polys = gpd.GeoDataFrame(
        {"geometry": [Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)]) for i in range(4)]},
        crs="EPSG:2927",
    )

    class _Grid:
        gdf_vorPolys = polys
        ncpl = 4

        def get_grid_edge(self):
            return [0, 3]

    grid = _Grid()
    return grid, ModelContext(grid=grid, domain=np.ones((1, 4), dtype=int))


def _write(tmp_path, geometry, **columns):
    """Write a one-row GeoPackage and return its path."""

    path = tmp_path / "bc.gpkg"
    gpd.GeoDataFrame(
        {**{k: [v] for k, v in columns.items()}, "geometry": [geometry]}, crs="EPSG:2927"
    ).to_file(path, layer="bc", driver="GPKG")
    return path


def _records(spec):
    """The period-0 records from a built PackageSpec."""

    return spec.options["stress_period_data"][0]


def test_extent_prefers_area_then_length():
    """A polygon measures by area, a line by length, a point by neither."""

    assert _extent(Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])) == pytest.approx(4.0)
    assert _extent(LineString([(0, 0), (3, 0)])) == pytest.approx(3.0)
    assert _extent(Point(1, 1)) == 0.0


def test_without_spread_the_value_is_broadcast(tmp_path, grid_and_context):
    """The documented prior behaviour: every intersected cell gets the whole value."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1,
                       elevation="elevation", conductance="cond", boundnames=False)
    records = _records(spec)
    assert len(records) == 4
    assert sum(r[-1] for r in records) == pytest.approx(400.0)


def test_spread_conserves_the_total_across_cells(tmp_path, grid_and_context):
    """A line crossing four cells hands each its own share, summing to the original."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    records = _records(spec)
    assert sum(r[-1] for r in records) == pytest.approx(100.0)
    # 0.5 + 1 + 1 + 0.5 of a 3-unit line
    assert sorted(round(r[-1], 6) for r in records) == [
        pytest.approx(100 / 6), pytest.approx(100 / 6),
        pytest.approx(100 / 3), pytest.approx(100 / 3),
    ]


def test_spread_leaves_intensive_fields_alone(tmp_path, grid_and_context):
    """Only the wrapped field is split; an elevation stays whole on every cell."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    assert {r[1] for r in _records(spec)} == {10.0}


def test_a_point_keeps_its_whole_value(tmp_path, grid_and_context):
    """A point has no extent to divide, so Spread is a no-op rather than a zero."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, Point(1.5, 0.5), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    records = _records(spec)
    assert len(records) == 1
    assert records[0][-1] == pytest.approx(100.0)


def test_clip_loses_the_part_outside_the_grid(tmp_path, grid_and_context):
    """Half the line hangs past the mesh, so half the conductance is not applied."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(2.0, 0.5), (6.0, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond", mode="clip"), boundnames=False)
    assert sum(r[-1] for r in _records(spec)) == pytest.approx(50.0)


def test_retained_keeps_the_total_on_what_is_left(tmp_path, grid_and_context):
    """The same line under mode='retained' concentrates its whole total inside."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(2.0, 0.5), (6.0, 0.5)]), layer=1, elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond", mode="retained"), boundnames=False)
    assert sum(r[-1] for r in _records(spec)) == pytest.approx(100.0)


def test_min_share_drops_a_corner_clip(tmp_path, grid_and_context):
    """A cell the feature barely grazes is dropped rather than given a sliver."""

    grid, ctx = grid_and_context
    # 99.5% in cell 1, a 0.5% nick of cell 2
    path = _write(tmp_path, LineString([(1.01, 0.5), (2.005, 0.5)]), layer=1,
                  elevation=10.0, cond=100.0)
    loose = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                        conductance=Spread("cond", min_share=0.0), boundnames=False)
    tight = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                        conductance=Spread("cond", min_share=0.01), boundnames=False)
    assert len(_records(loose)) == 2
    assert len(_records(tight)) == 1


def test_spread_rejects_an_unknown_mode():
    """A typo in mode is refused at construction, not silently treated as clip."""

    with pytest.raises(ValueError, match="clip.*retained"):
        Spread("cond", mode="area")


def test_spread_rejects_an_out_of_range_min_share():
    """min_share is a fraction; 1.0 would drop every cell."""

    with pytest.raises(ValueError, match="min_share"):
        Spread("cond", min_share=1.0)


def test_ghb_conductance_also_spreads(tmp_path, grid_and_context):
    """Spread is a value spec, not a DRN special case -- it rides any package's field."""

    grid, ctx = grid_and_context
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1, head=12.0, cond=90.0)
    spec = mf.ghb.gpkg(path, layer="bc", context=ctx, nper=1, head="head",
                       conductance=Spread("cond"), boundnames=False)
    records = _records(spec)
    assert sum(r[-1] for r in records) == pytest.approx(90.0)
    assert {r[1] for r in records} == {12.0}


def test_every_row_value_dataclass_survives_a_run_manifest():
    """Each dataclass in the RowValue union must convert to plain JSON.

    A run manifest is written with ``json.dumps`` and no ``default=``, so a spec
    that reaches it unconverted raises ``TypeError`` at ``prepare_run`` -- long
    after the package itself built cleanly. This fails the moment a new variant
    is added to the union without a ``_metadata_value`` branch.
    """

    import dataclasses
    import json
    import typing

    from myflopy.geopackage import CellSurfaceOffset, RowValue, _metadata_value

    dataclass_variants = [
        arg for arg in typing.get_args(RowValue)
        if isinstance(arg, type) and dataclasses.is_dataclass(arg)
    ]
    assert dataclass_variants, "RowValue should carry at least one dataclass variant"

    samples = {
        CellSurfaceOffset: CellSurfaceOffset("cell_top", offset=-1.0),
        Spread: Spread("conductance"),
    }
    missing = [cls.__name__ for cls in dataclass_variants if cls not in samples]
    assert not missing, f"add a sample for {missing} so this test still covers the union"

    for cls in dataclass_variants:
        converted = _metadata_value(samples[cls])
        assert isinstance(converted, dict), f"{cls.__name__} did not convert"
        assert converted["type"] == cls.__name__
        json.dumps(converted)  # the actual failure mode


def test_spread_metadata_records_its_settings():
    """The manifest keeps enough to reproduce the split, not just the field name."""

    from myflopy.geopackage import _metadata_value

    assert _metadata_value(Spread("cond", mode="retained", min_share=0.05)) == {
        "type": "Spread",
        "field": "cond",
        "mode": "retained",
        "min_share": 0.05,
    }


def test_passthrough_cells_do_not_receive_boundaries(tmp_path):
    """idomain < 0 is vertical passthrough, not active -- MF6 refuses a boundary there.

    ``bool(-1)`` is ``True``, so a truthiness test lets these through and the run
    fails at read time with "Cell is outside active grid domain". A ``LayerStack``
    with ``pinch="passthrough"`` produces them in quantity.
    """

    from myflopy.specs import ModelContext

    polys = gpd.GeoDataFrame(
        {"geometry": [Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)]) for i in range(4)]},
        crs="EPSG:2927",
    )

    class _Grid:
        gdf_vorPolys = polys
        ncpl = 4

        def get_grid_edge(self):
            return [0, 3]

    # cell 1 active, cell 2 passthrough, cell 3 inactive, cell 4 active
    domain = np.array([[1, -1, 0, 1]])
    ctx = ModelContext(grid=_Grid(), domain=domain)

    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1,
                  elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance=Spread("cond"), boundnames=False)
    cells = sorted(r[0][1] for r in _records(spec))
    assert cells == [0, 3], "only the two truly active cells may carry a record"


def test_one_dimensional_domain_also_excludes_passthrough(tmp_path):
    """The single-layer branch of _active needs the same test, not truthiness."""

    from myflopy.specs import ModelContext

    polys = gpd.GeoDataFrame(
        {"geometry": [Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)]) for i in range(4)]},
        crs="EPSG:2927",
    )

    class _Grid:
        gdf_vorPolys = polys
        ncpl = 4

        def get_grid_edge(self):
            return [0, 3]

    ctx = ModelContext(grid=_Grid(), domain=np.array([1, -1, 0, 1]))
    path = _write(tmp_path, LineString([(0.5, 0.5), (3.5, 0.5)]), layer=1,
                  elevation=10.0, cond=100.0)
    spec = mf.drn.gpkg(path, layer="bc", context=ctx, nper=1, elevation="elevation",
                       conductance="cond", boundnames=False)
    assert sorted(r[0][1] for r in _records(spec)) == [0, 3]
