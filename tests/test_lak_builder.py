from __future__ import annotations

from unittest import mock

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Polygon

from myflopy import (
    LAKBuilder,
    LakeConnection,
    LakeOutlet,
    LakeTableBuilder,
    ModelContext,
    ModelSpec,
    PackageSpec,
    SimulationSpec,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _grid() -> VoronoiGridPlus:
    vertices = np.array(
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
    grid = VoronoiGridPlus(verts=vertices, iverts=iverts, xcyc=centers)
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [12.0] * 6, 1: [8.0] * 6, 2: [0.0] * 6},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    grid.get_disu_connectivity(validate=True)
    return grid


def _lakes() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "name": ["natural", "trench"],
            "stage": [11.0, 10.0],
            "bottom": [9.0, 6.0],
        },
        geometry=[
            Polygon([(0.05, 1.05), (0.95, 1.05), (0.95, 1.95), (0.05, 1.95)]),
            Polygon([(1.05, 0.05), (1.95, 0.05), (1.95, 0.95), (1.05, 0.95)]),
        ],
        crs=_grid().crs,
    )


def _builder(**updates) -> LAKBuilder:
    values = {
        "context": ModelContext(grid=_grid(), domain=np.ones((2, 6), dtype=int)),
        "nper": 2,
        "lakes": _lakes(),
        "lake_id_field": "name",
        "starting_stage": "stage",
        "lake_bottom": "bottom",
        "lake_top": {"trench": 10.0},  # rectangular facilities need an explicit flat top
        "bed_leakance": 0.1,
        "connection_modes": {"natural": "bathy", "trench": "rectangular"},
        "status": {"natural": ["ACTIVE", "ACTIVE"], "trench": ["ACTIVE", "INACTIVE"]},
    }
    values.update(updates)
    return LAKBuilder(**values)


def test_lak_builder_uses_stable_ids_and_mixed_generated_modes():
    builder = _builder()

    assert builder.lake_ids == ("natural", "trench")
    assert builder.lake_numbers == {"natural": 0, "trench": 1}
    assert builder.lake_cells == {"natural": [3], "trench": [1]}
    assert any(lake_id == "trench" and item.connection_type == "HORIZONTAL" for lake_id, item in builder.connections)
    assert builder.packagedata[0][-1] == "natural"
    assert builder.perioddata[1] == [[0, "STATUS", "ACTIVE"], [1, "STATUS", "INACTIVE"]]


def test_lak_builder_accepts_explicit_connections_inside_connection_modes():
    connection = LakeConnection(
        cellid=(1, 1),
        connection_type="HORIZONTAL",
        bottom_elevation=4.0,
        top_elevation=7.0,
        connection_length=2.0,
        connection_width=3.0,
    )
    builder = _builder(
        lakes=_lakes().iloc[[1]].copy(),
        connection_modes={"trench": (connection,)},
        starting_stage={"trench": 10.0},
        lake_bottom=None,
        status=None,
    )

    assert builder.connections == (("trench", connection),)
    assert builder.connectiondata[0][2:] == [(1, 1), "HORIZONTAL", 0.1, 4.0, 7.0, 2.0, 3.0]


def test_lak_builder_supports_connection_specific_leakance_and_outlets():
    builder = _builder(
        bed_leakance={
            "natural": 0.1,
            "trench": 0.2,
            ("trench", "horizontal"): 0.05,
        },
        outlets=(LakeOutlet("natural", "trench", rate=[2.0, 3.0]),),
    )

    trench_horizontal = next(
        row for row in builder.connectiondata if row[0] == 1 and row[3] == "HORIZONTAL"
    )
    assert trench_horizontal[4] == 0.05
    assert builder.outletdata[0][1:4] == [0, 1, "SPECIFIED"]
    assert builder.perioddata[1][-1] == [0, "RATE", 3.0]


def test_lake_table_builder_and_zero_argument_build():
    builder = _builder(
        tables={
            "natural": LakeTableBuilder(area=100.0, bottom=9.0, top=11.0, stage_step=1.0),
        }
    )

    table = builder.prepared_tables["natural"]
    spec = builder.build()

    assert table.rows == ((9.0, 0.0, 100.0), (10.0, 100.0, 100.0), (11.0, 200.0, 100.0))
    assert spec.metadata["builder"] == "LAKBuilder"
    assert spec.options["ntables"] == 1
    with pytest.raises(TypeError):
        builder.build(name="other")


def test_lake_table_builder_can_derive_rows_from_dem(tmp_path):
    dem = tmp_path / "lake_dem.tif"
    with rasterio.open(
        dem,
        "w",
        driver="GTiff",
        height=2,
        width=2,
        count=1,
        dtype="float32",
        crs="EPSG:2927",
        transform=from_origin(0.0, 2.0, 1.0, 1.0),
    ) as destination:
        destination.write(np.array([[8.0, 9.0], [9.0, 10.0]], dtype="float32"), 1)

    table = LakeTableBuilder(
        dem=dem,
        footprint=Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),
        stages=(9.0, 10.0),
    ).build()

    assert table.rows == ((9.0, 1.0, 1.0), (10.0, 4.0, 3.0))


def test_lak_builder_writes_real_flopy_310_package(tmp_path):
    builder = _builder(
        connection_modes="automatic",
        tables={"natural": LakeTableBuilder(rows=[(9.0, 0.0, 10.0), (11.0, 20.0, 10.0)])},
    )
    grid = builder.grid
    gridprops = grid.get_gridprops_vertexgrid()
    flow = ModelSpec(
        "flow",
        "gwf",
        packages=(
            PackageSpec(
                "disv",
                flopy.mf6.ModflowGwfdisv,
                {
                    "nlay": 2,
                    "ncpl": grid.ncpl,
                    "nvert": len(gridprops["vertices"]),
                    "vertices": gridprops["vertices"],
                    "cell2d": gridprops["cell2d"],
                    "top": grid.gdf_topbtm[0].tolist(),
                    "botm": [grid.gdf_topbtm[1].tolist(), grid.gdf_topbtm[2].tolist()],
                },
            ),
            builder.build(),
        ),
    )
    simulation = SimulationSpec(
        "lak",
        models=(flow,),
        packages=(
            PackageSpec("tdis", flopy.mf6.ModflowTdis, {"nper": 2, "perioddata": [(1.0, 1, 1.0)] * 2}),
        ),
    )

    built = simulation.build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert (tmp_path / "flow.lak").exists()
    assert (tmp_path / "lak_natural.lak.tab").exists()


def _wide_lake_builder() -> LAKBuilder:
    """A single lake covering >1 cell -- exercises per-cell sidewall building."""
    grid = _grid()
    # Spans cells 3 (centre 0.5,1.5) and 4 (centre 1.5,1.5) on the top row.
    lakes = gpd.GeoDataFrame(
        {"name": ["wide"]},
        geometry=[Polygon([(0.05, 1.05), (1.95, 1.05), (1.95, 1.95), (0.05, 1.95)])],
        crs=grid.crs,
    )
    return LAKBuilder(
        context=ModelContext(grid=grid, domain=np.ones((2, 6), dtype=int)),
        nper=1,
        lakes=lakes,
        lake_id_field="name",
        starting_stage={"wide": 11.0},
        lake_bottom={"wide": 8.0},
        # A flat-bottomed lake is 'rectangular' (edges-only sidewalls), which is
        # what actually exercises per-cell sidewall building -- flat 'bathy' has
        # no exposed steps and would build zero horizontal connections.
        lake_top={"wide": 11.0},
        connection_modes="rectangular",
        bed_leakance=0.1,
    )


def test_lake_cells_is_cached_not_recomputed_per_cell():
    """lake_cells does a full-grid intersection per lake; building connections must
    not recompute it once per lake cell (that O(cells x lakes x grid) blowup hung
    automatic builds over large lakes)."""
    builder = _wide_lake_builder()
    assert len(builder.lake_cells["wide"]) >= 2  # multi-cell, so per-cell recompute would show

    # Cached: repeated access returns the same object (no recompute).
    assert builder.lake_cells is builder.lake_cells

    original = gpd.GeoSeries.intersection
    calls = {"n": 0}

    def counting_intersection(self, *args, **kwargs):
        calls["n"] += 1
        return original(self, *args, **kwargs)

    object.__setattr__(builder, "_lake_cells", None)  # force one recompute
    with mock.patch.object(gpd.GeoSeries, "intersection", counting_intersection):
        _ = builder.connections  # calls _sidewall_connections for every lake cell

    # One intersection for the single lake -- NOT one per lake cell.
    assert calls["n"] == 1


def test_lake_cells_cache_is_not_shared_after_with_updates():
    builder = _wide_lake_builder()
    _ = builder.lake_cells  # populate the cache
    assert builder._lake_cells is not None

    updated = builder.with_updates(nper=2)
    assert updated._lake_cells is None  # a fresh builder recomputes from its own inputs
    assert updated.lake_cells == builder.lake_cells


# --- connection geometry: bathy vs rectangular ------------------------------
# Regression coverage for the connection-generation rewrite. The bug being pinned:
# horizontal faces were clipped to the (transient) starting stage instead of the
# facility top / neighbor bottom, which choked lake-aquifer leakage and let lakes
# mound tens of feet above their rim. Grid _grid() is 3x2, surfaces [12, 8, 0]
# -> layer 0 = [8, 12], layer 1 = [0, 8]. Cells: bottom row 0,1,2 / top row 3,4,5.


def _horizontals(builder, lake_id):
    return [c for lid, c in builder.connections
            if lid == lake_id and c.connection_type == "HORIZONTAL"]


def _verticals(builder, lake_id):
    return [c for lid, c in builder.connections
            if lid == lake_id and c.connection_type == "VERTICAL"]


def _rect_builder(**updates) -> LAKBuilder:
    grid = _grid()
    lakes = gpd.GeoDataFrame(
        {"name": ["basin"]},
        geometry=[Polygon([(0.05, 1.05), (1.95, 1.05), (1.95, 1.95), (0.05, 1.95)])],
        crs=grid.crs,
    )  # covers cells 3 and 4 (top row, left + middle)
    values = dict(
        context=ModelContext(grid=grid, domain=np.ones((2, 6), dtype=int)),
        nper=1,
        lakes=lakes,
        lake_id_field="name",
        starting_stage={"basin": 7.0},   # deliberately BELOW lake_top
        lake_bottom={"basin": 6.0},
        lake_top={"basin": 11.0},
        bed_leakance=0.1,
        connection_modes="rectangular",
    )
    values.update(updates)
    return LAKBuilder(**values)


def test_rectangular_telev_is_lake_top_edges_only_and_multilayer():
    builder = _rect_builder()
    horizontals = _horizontals(builder, "basin")

    # telev comes from lake_top (11) / layer boundary (8) -- NEVER the start stage (7).
    tops = {c.top_elevation for c in horizontals}
    assert 7.0 not in tops
    assert tops == {11.0, 8.0}

    # Edges only: the shared 3<->4 interior face is skipped. Perimeter faces are
    # 3->0, 4->1, 4->5 (3 faces); the box [6, 11] spans both layers -> 6 horizontals.
    assert len(horizontals) == 6
    assert {c.cellid[0] for c in horizontals} == {0, 1}  # both layers

    # belev is clipped to each layer: layer 0 -> max(8,6)=8, layer 1 -> max(0,6)=6.
    assert {c.bottom_elevation for c in horizontals} == {8.0, 6.0}
    assert len(_verticals(builder, "basin")) == 2  # one per lake cell


def test_rectangular_vault_has_no_horizontals():
    builder = _rect_builder(only_vertical=True)
    assert _horizontals(builder, "basin") == []
    assert len(_verticals(builder, "basin")) == 2


def test_bathy_flat_bottom_multi_cell_is_guarded():
    """A multi-cell 'bathy' lake with a flat bottom has no exposed steps, so it
    emits zero horizontal connections -- a lake with no lateral aquifer exchange.
    The geometry is internally consistent, so it WARNS (rather than raising) and
    names the three real choices."""

    builder = _rect_builder(connection_modes="bathy")  # flat bottom 6.0, cells 3 & 4
    with pytest.warns(UserWarning, match="no horizontal") as record:
        connections = builder.connections
    assert _horizontals(builder, "basin") == []  # still builds, just no sidewalls
    assert connections  # vertical connections remain
    message = str(record[0].message)
    assert "rectangular" in message and "only_vertical" in message and "bathymetry" in message


def test_bathy_single_cell_flat_lake_is_not_guarded():
    """The guard is only for the multi-cell degeneracy; a single-cell bathy lake
    legitimately has no interior faces and must still build without warning."""

    grid = _grid()
    lakes = gpd.GeoDataFrame(
        {"name": ["pond"]},
        geometry=[Polygon([(0.05, 1.05), (0.95, 1.05), (0.95, 1.95), (0.05, 1.95)])],
        crs=grid.crs,
    )  # one cell (cell 3)
    builder = LAKBuilder(
        context=ModelContext(grid=grid, domain=np.ones((2, 6), dtype=int)),
        nper=1,
        lakes=lakes,
        lake_id_field="name",
        starting_stage={"pond": 7.0},
        lake_bottom={"pond": 6.0},
        lake_top={"pond": 11.0},
        bed_leakance=0.1,
        connection_modes="bathy",
    )
    import warnings

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        verticals = _verticals(builder, "pond")
    assert len(verticals) == 1  # builds
    assert not [w for w in caught if "no horizontal" in str(w.message)]  # no guard warning


def test_rectangular_requires_lake_top():
    builder = _rect_builder(lake_top=None)
    with pytest.raises(ValueError, match="lake_top is required"):
        builder.validate()


def test_rectangular_interior_connects_all_faces():
    # A permeable-fill (gravel) basin connects through interior faces too, not just the
    # perimeter -- reproduces the legacy all-faces connectivity and stabilizes a small
    # basin that drains near-empty. The lake covers cells 3 & 4 (see _rect_builder).
    edges = _horizontals(_rect_builder(), "basin")                          # default: edges only
    filled = _horizontals(_rect_builder(rectangular_interior=True), "basin")  # all faces

    # Edges-only skips the shared 3<->4 face; interior adds it (2 faces x 2 layers).
    assert len(edges) == 6
    assert len(filled) == 10
    assert len(filled) > len(edges)

    # A per-lake mapping works too, and the face geometry is unchanged (telev == lake_top
    # / layer boundary, never the starting stage).
    mapped = _horizontals(_rect_builder(rectangular_interior={"basin": True}), "basin")
    assert len(mapped) == 10
    assert {c.top_elevation for c in filled} == {11.0, 8.0}


def test_bathy_connects_only_up_exposed_steps_and_spans_layers():
    grid = _grid()
    lakes = gpd.GeoDataFrame(
        {"name": ["nat"]},
        geometry=[Polygon([(0.05, 1.05), (1.95, 1.05), (1.95, 1.95), (0.05, 1.95)])],
        crs=grid.crs,
    )  # cells 3 and 4
    # Per-cell lake bottom: cell 3 is deep (2), cell 4 shallow (9); neighbors high (10).
    lake_bottom = pd.Series({0: 10.0, 1: 10.0, 2: 10.0, 3: 2.0, 4: 9.0, 5: 10.0})
    builder = LAKBuilder(
        context=ModelContext(grid=grid, domain=np.ones((2, 6), dtype=int)),
        nper=1,
        lakes=lakes,
        lake_id_field="name",
        starting_stage={"nat": 10.0},
        lake_bottom={"nat": lake_bottom},
        bed_leakance=0.1,
        connection_modes="bathy",
    )
    by_cell: dict[int, list] = {}
    for connection in _horizontals(builder, "nat"):
        by_cell.setdefault(connection.cellid[1], []).append(connection)

    # Cell 3 (deepest) has an exposed step toward BOTH neighbors 0 and 4; the face
    # [2, neighbor_bottom] crosses both model layers -> 2 neighbors x 2 layers = 4.
    assert len(by_cell[3]) == 4
    assert {c.cellid[0] for c in by_cell[3]} == {0, 1}
    # telev tracks the neighbor's lake bottom (9 from cell 4), not the stage.
    assert any(c.top_elevation == 9.0 for c in by_cell[3])

    # Cell 4 (lb 9) connects only UP to the higher neighbors 1 and 5 -- never toward
    # the deeper lake cell 3 -- and only in the single overlapping layer.
    assert len(by_cell[4]) == 2
    assert {c.cellid[0] for c in by_cell[4]} == {0}


def test_bathy_scalar_bottom_makes_no_exposed_steps():
    # A flat (scalar) bottom in bathy mode has no steps -> vertical connections
    # only. Internally consistent but degenerate, so it warns (see
    # test_bathy_flat_bottom_multi_cell_is_guarded).
    builder = _rect_builder(connection_modes="bathy", lake_bottom={"basin": 6.0},
                            lake_top=None)
    with pytest.warns(UserWarning, match="no horizontal"):
        assert _horizontals(builder, "basin") == []
    assert len(_verticals(builder, "basin")) == 2
