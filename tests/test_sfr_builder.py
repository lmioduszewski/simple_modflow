from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point

from myflopy import (
    ModelContext,
    ModelSpec,
    PackageSpec,
    SFRBuilder,
    SimulationSpec,
    StreamConnection,
    StreamDiversion,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _grid() -> VoronoiGridPlus:
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
    xcyc = np.array(
        [[0.5, 0.5], [1.5, 0.5], [2.5, 0.5], [0.5, 1.5], [1.5, 1.5], [2.5, 1.5]]
    )
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10, 9, 8, 11, 10, 9], 1: [0] * 6},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    return grid


def _streams() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "name": ["tributary", "main"],
            "from_node": ["headwater", "main_start"],
            "to_node": ["junction", "outlet"],
        },
        geometry=[
            LineString([(0.2, 1.5), (1.5, 0.5)]),
            LineString([(0.1, 0.5), (2.9, 0.5)]),
        ],
        crs=_grid().crs,
    )


def _builder(**updates) -> SFRBuilder:
    values = {
        "context": ModelContext(grid=_grid(), domain=np.ones((1, 6), dtype=int)),
        "nper": 2,
        "streams": _streams(),
        "stream_id": "name",
        "connection_tolerance": 0.05,
        "width": {"tributary": 4.0, "main": 8.0},
        "gradient": 0.001,
        "roughness": 0.03,
        "streambed_k": 1.0,
        "streambed_thickness": 1.0,
    }
    values.update(updates)
    return SFRBuilder(**values)


def test_sfr_builder_automatically_resolves_confluence_before_grid_mapping():
    builder = _builder()

    assert builder.network.connections == (StreamConnection("tributary", "main"),)
    assert builder.network.unresolved == ("main",)
    assert builder.stream_ids == ("tributary", "main")

    source, receiver = builder.resolved_connections[0]
    assert source == builder.outlet_reach("tributary")
    assert receiver in builder.stream_reaches["main"]
    assert -receiver in builder.connectiondata[source]
    assert source in builder.connectiondata[receiver]


def test_sfr_builder_supports_explicit_locations_and_node_connections():
    explicit = _builder(
        connection_mode="explicit",
        connections=(
            StreamConnection(
                "tributary",
                "main",
                source_location="downstream",
                receiver_location=Point(1.5, 0.5),
            ),
            StreamConnection("main", None),
        ),
    )
    nodes = _builder(
        connection_mode="nodes",
        from_node="from_node",
        to_node="to_node",
        connections=(StreamConnection("tributary", "main"),),
    )

    assert explicit.resolved_connections[0][1] in explicit.stream_reaches["main"]
    assert nodes.network.connections[0] == StreamConnection("tributary", "main")


def test_sfr_builder_broadcasts_scalar_rates_but_keeps_inflow_a_point_source():
    """A scalar RATE (rainfall/evaporation) applies to EVERY reach; a scalar
    INFLOW is a volumetric point source at the single headwater reach.

    Regression: a scalar rate used to land on the first reach only, so e.g.
    ``rainfall=0.002`` silently wetted one reach of the whole network.
    """

    builder = _builder(rainfall=0.002, inflow=500.0)
    all_reaches = sorted(int(r) for r in builder.reaches.index)
    period0 = builder.perioddata[0]

    rain = [row for row in period0 if row[1] == "RAINFALL"]
    assert sorted(row[0] for row in rain) == all_reaches  # one row per reach
    assert all(row[2] == 0.002 for row in rain)

    inflow = [row for row in period0 if row[1] == "INFLOW"]
    assert [row[0] for row in inflow] == [int(builder.reaches.index[0])]
    assert inflow[0][2] == 500.0


def test_sfr_builder_supports_multiple_reaches_in_one_cell_and_diversions():
    streams = gpd.GeoDataFrame(
        {"name": ["source", "canal"]},
        geometry=[
            LineString([(0.1, 0.3), (2.9, 0.3)]),
            LineString([(1.1, 0.6), (2.9, 0.6)]),
        ],
        crs=_grid().crs,
    )
    builder = _builder(
        streams=streams,
        connection_mode="explicit",
        connections=(StreamConnection("source", None), StreamConnection("canal", None)),
        diversions=(StreamDiversion("source", "canal", amount={0: 0.25, 1: 0.5}),),
    )

    assert builder.reaches["cellid"].duplicated().any()
    assert builder.diversiondata[0][3] == "FRACTION"
    assert builder.perioddata[1][-1][-1] == 0.5


def test_sfr_builder_build_is_immutable_and_accepts_no_configuration():
    base = _builder()
    wide = base.with_updates(width=20.0)

    assert base.packagedata[0][3] == 4.0
    assert wide.packagedata[0][3] == 20.0
    assert wide.build().metadata["builder"] == "SFRBuilder"
    with pytest.raises(TypeError):
        wide.build(name="other")


def test_sfr_builder_preserves_explicit_unit_conversions():
    builder = _builder(length_conversion=3.28081, time_conversion=86_400.0)

    spec = builder.build()

    assert spec.options["length_conversion"] == 3.28081
    assert spec.options["time_conversion"] == 86_400.0


def test_sfr_builder_defaults_to_feet_days_from_context_units():
    # The default ModelContext is feet/days, so Manning's runs in the model's units
    # without the caller setting the conversions. (The old SI 1.0/1.0 default silently
    # ran Manning's in meters/seconds -> bogus stream stage on a feet/days model.)
    spec = _builder().build()

    assert spec.options["length_conversion"] == pytest.approx(3.28081)
    assert spec.options["time_conversion"] == pytest.approx(86_400.0)


def test_sfr_builder_derives_conversions_from_context_units():
    ctx = ModelContext(
        grid=_grid(), domain=np.ones((1, 6), dtype=int),
        length_units="meters", time_units="seconds",
    )
    spec = _builder(context=ctx).build()

    assert spec.options["length_conversion"] == pytest.approx(1.0)
    assert spec.options["time_conversion"] == pytest.approx(1.0)


def test_mf6_unit_conversion_helpers():
    from myflopy.specs import mf6_length_conversion, mf6_time_conversion

    assert mf6_length_conversion("feet") == pytest.approx(3.28081)
    assert mf6_length_conversion("FEET") == pytest.approx(3.28081)  # case-insensitive
    assert mf6_length_conversion("meters") == 1.0
    assert mf6_length_conversion("unknown") == 1.0
    assert mf6_time_conversion("days") == 86_400.0
    assert mf6_time_conversion("hours") == 3600.0
    assert mf6_time_conversion("seconds") == 1.0
    with pytest.raises(ValueError):
        mf6_length_conversion("furlongs")
    with pytest.raises(ValueError):
        mf6_time_conversion("fortnights")


def test_sfr_builder_writes_real_flopy_310_package(tmp_path):
    builder = _builder(connection_mode="explicit", connections=(StreamConnection("tributary", "main"), StreamConnection("main", None)))
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
                    "nlay": 1,
                    "ncpl": grid.ncpl,
                    "nvert": len(gridprops["vertices"]),
                    "vertices": gridprops["vertices"],
                    "cell2d": gridprops["cell2d"],
                    "top": grid.gdf_topbtm[0].tolist(),
                    "botm": [grid.gdf_topbtm[1].tolist()],
                },
            ),
            builder.build(),
        ),
    )
    simulation = SimulationSpec(
        "sfr",
        models=(flow,),
        packages=(
            PackageSpec("tdis", flopy.mf6.ModflowTdis, {"nper": 2, "perioddata": [(1.0, 1, 1.0)] * 2}),
        ),
    )

    built = simulation.build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert (tmp_path / "flow.sfr").exists()
