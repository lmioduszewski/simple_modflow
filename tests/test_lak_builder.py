from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
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
        "bed_leakance": 0.1,
        "connection_modes": {"natural": "automatic", "trench": "rectangular"},
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
