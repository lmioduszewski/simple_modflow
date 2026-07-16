from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Polygon

from myflopy import ModelContext, ModelSpec, PackageSpec, SimulationSpec
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.lakes import LAKBuilder
from myflopy.modflow.mf6.mvr import Move, MoverConnection, MVRBuilder
from myflopy.modflow.mf6.sfr import SFRBuilder


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


def _lake_builder(grid: VoronoiGridPlus) -> LAKBuilder:
    lakes = gpd.GeoDataFrame(
        {"name": ["pond"], "stage": [11.0], "bottom": [9.0]},
        geometry=[Polygon([(0.05, 1.05), (0.95, 1.05), (0.95, 1.95), (0.05, 1.95)])],
        crs=grid.crs,
    )
    return LAKBuilder(
        context=ModelContext(grid=grid, domain=np.ones((2, 6), dtype=int)),
        nper=2,
        lakes=lakes,
        lake_id_field="name",
        starting_stage="stage",
        lake_bottom="bottom",
        bed_leakance=0.1,
        status="ACTIVE",
        mover=True,
    )


def _sfr_builder(grid: VoronoiGridPlus) -> SFRBuilder:
    streams = gpd.GeoDataFrame(
        {"name": ["main_stem"]},
        geometry=[LineString([(0.05, 0.5), (2.95, 0.5)])],
        crs=grid.crs,
    )
    return SFRBuilder(
        context=ModelContext(grid=grid, domain=np.ones((2, 6), dtype=int)),
        nper=2,
        streams=streams,
        stream_id="name",
        inflow={0: [(0, 1.0)], 1: [(0, 1.0)]},
        width=1.0,
        gradient=0.001,
        roughness=0.03,
        streambed_k=0.1,
        streambed_thickness=1.0,
        mover=True,
    )


def test_mvr_builder_normalizes_moves_and_declares_packages():
    move = Move(
        source=MoverConnection("SFR", 2),
        receiver=MoverConnection("LAK", 0),
        value=0.25,
    )
    builder = MVRBuilder(nper=3, moves=(move,))

    assert builder.perioddata == {
        0: [["sfr", 2, "lak", 0, "FACTOR", 0.25]],
        1: [["sfr", 2, "lak", 0, "FACTOR", 0.25]],
        2: [["sfr", 2, "lak", 0, "FACTOR", 0.25]],
    }
    assert builder.packages == [["sfr"], ["lak"]]
    assert builder.build().metadata["builder"] == "MVRBuilder"
    with pytest.raises(TypeError):
        builder.build(moves=())


def test_mvr_endpoint_helpers_resolve_stable_ids_and_locations():
    grid = _grid()
    lak = _lake_builder(grid)
    sfr = _sfr_builder(grid)

    assert lak.connection("pond") == MoverConnection("lak", 0)
    assert sfr.connection("main_stem", "upstream") == MoverConnection("sfr", sfr.stream_reaches["main_stem"][0])
    assert sfr.connection("main_stem") == MoverConnection("sfr", sfr.stream_reaches["main_stem"][-1])


def test_mvr_builder_writes_real_flopy_310_package(tmp_path):
    grid = _grid()
    lak = _lake_builder(grid)
    sfr = _sfr_builder(grid)
    mvr = MVRBuilder(
        nper=2,
        moves={
            period: [Move(sfr.connection("main_stem"), lak.connection("pond"), value=0.5)]
            for period in range(2)
        },
    )
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
            lak.build(),
            sfr.build(),
            mvr.build(),
        ),
    )
    simulation = SimulationSpec(
        "mvr",
        models=(flow,),
        packages=(
            PackageSpec("tdis", flopy.mf6.ModflowTdis, {"nper": 2, "perioddata": [(1.0, 1, 1.0)] * 2}),
        ),
    )

    built = simulation.build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert (tmp_path / "flow.mvr").exists()
