from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Polygon

from myflopy import (
    CellSurfaceOffset,
    GeoPackageSource,
    ModelContext,
    Project,
    SimpleModelConfig,
    simple_model_spec,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _two_cell_grid() -> VoronoiGridPlus:
    verts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ],
        dtype=float,
    )
    iverts = [[0, 3, 2, 1], [1, 2, 5, 4]]
    xcyc = np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float)
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0, 9.0], 1: [5.0, 6.0]},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    return grid


def _write_inputs(path):
    gdf = gpd.GeoDataFrame(
        {
            "name": ["left", "right"],
            "layer": [1, 1],
            "head_0": [10.0, 9.0],
            "head_1": [10.5, 8.5],
            "conductance": [2.0, 3.0],
            "elevation": [8.0, 7.0],
            "stage": [9.5, 8.5],
            "rbot": [7.5, 6.5],
            "surface": [10.0, 9.0],
            "et_rate": [0.002, 0.003],
            "depth": [2.5, 2.0],
            "rate": [-1.0, -2.0],
            "recharge": [0.001, 0.002],
            "k": [5.0, 10.0],
            "top_offset": [-1.0, -0.5],
            "height_above_bottom": [1.0, 2.0],
            "min_elev": [8.5, 7.5],
        },
        geometry=[
            Polygon([(0.0, 0.0), (0.99, 0.0), (0.99, 1.0), (0.0, 1.0)]),
            Polygon([(1.01, 0.0), (2.0, 0.0), (2.0, 1.0), (1.01, 1.0)]),
        ],
        crs="EPSG:2927",
    )
    gdf.to_file(path, driver="GPKG")


def test_geopackage_source_builds_specs_and_arrays(tmp_path):
    vor = _two_cell_grid()
    path = tmp_path / "boundaries.gpkg"
    _write_inputs(path)
    source = GeoPackageSource(
        path,
        ModelContext(grid=vor, domain=np.array([[1, 0]])),
        nper=2,
    )

    chd = source.chd(head=["head_0", "head_1"])
    ghb = source.ghb(head="head_0")
    drn = source.drn()
    riv = source.riv()
    evt = source.evt(rate="et_rate")
    wel = source.wel()
    rch = source.rch()
    k = source.k_array(value="k", nlay=1, defaults=1.0)

    assert chd.options["stress_period_data"][0] == [[(0, 0), 10.0, "left"]]
    assert chd.options["stress_period_data"][1] == [[(0, 0), 10.5, "left"]]
    assert ghb.options["stress_period_data"][0] == [[(0, 0), 10.0, 2.0, "left"]]
    assert drn.options["stress_period_data"][0] == [[(0, 0), 8.0, 2.0, "left"]]
    assert riv.options["stress_period_data"][0] == [[(0, 0), 9.5, 2.0, 7.5, "left"]]
    assert riv.metadata["fields"] == {
        "stage": "stage", "conductance": "conductance", "rbot": "rbot",
    }
    assert evt.options["stress_period_data"][0] == [[(0, 0), 10.0, 0.002, 2.5, "left"]]
    assert evt.options["nseg"] == 1
    assert evt.metadata["fields"] == {
        "surface": "surface", "rate": "et_rate", "depth": "depth",
    }
    assert wel.options["stress_period_data"][0] == [[(0, 0), -1.0, "left"]]
    assert rch.options["stress_period_data"][0] == [[(0, 0), 0.001, "left"]]
    assert k.tolist() == [[5.0, 1.0]]
    assert chd.metadata["source_type"] == "geopackage"
    assert chd.metadata["fields"] == {"head": ["head_0", "head_1"]}


def test_geopackage_specs_run_through_project(tmp_path):
    vor = _two_cell_grid()
    path = tmp_path / "boundaries.gpkg"
    _write_inputs(path)

    baseline = simple_model_spec(
        SimpleModelConfig(
            vor=vor,
            name="gpkg_flow",
            top=[10.0, 9.0],
            bottom=[[0.0, 0.0]],
            initial_heads=[10.0, 9.0],
            k=[1.0, 1.0],
            save_specific_discharge=False,
            sto_transient={},
        )
    )
    source = GeoPackageSource(path, baseline.model("gpkg_flow").context, nper=1)
    flow = baseline.model("gpkg_flow")
    flow = flow.with_package(source.chd(head="head_0"))
    flow = flow.with_package(source.ghb(head="head_0"))
    flow = flow.with_package(source.riv())
    flow = flow.with_package(source.evt(rate="et_rate"))
    simulation = baseline.with_model(flow)

    run = Project(tmp_path / "project").run("baseline", simulation)

    assert run.success is True
    assert isinstance(run.built.models["gpkg_flow"].packages["ghb"], flopy.mf6.ModflowGwfghb)
    assert isinstance(run.built.models["gpkg_flow"].packages["riv"], flopy.mf6.ModflowGwfriv)
    assert isinstance(run.built.models["gpkg_flow"].packages["evt"], flopy.mf6.ModflowGwfevt)

    # the riv/evt explorer surfaces work end-to-end on the completed run:
    # registry fields in the input table, earth input map, RdBu signed-q results
    from myflopy import load_mf6_run

    view = load_mf6_run(run.workspace)
    inputs = view.packages.riv.inputs.get()
    assert {"stage", "cond", "rbot"}.issubset(inputs.columns)
    assert float(inputs.loc[inputs["cell"] == 0, "stage"].iloc[0]) == 9.5
    stage_map = view.packages.riv.inputs.map(per=0)
    assert stage_map.colorscale == "earth"
    q = view.packages.riv.results.q.get()
    assert not q.empty  # RIV budget term recorded by the run
    q_map = view.packages.riv.results.q.map(per=0)
    assert q_map.colorscale == "RdBu"

    evt_inputs = view.packages.evt.inputs.get()
    assert {"surface", "rate", "depth"}.issubset(evt_inputs.columns)
    assert view.packages.evt.inputs.map(per=0).colorscale == "earth"
    evt_q = view.packages.evt.results.q.get()
    assert not evt_q.empty  # EVT budget term recorded by the run
    assert view.packages.evt.results.q.map(per=0).colorscale == "RdBu"


def test_geopackage_source_supports_long_period_format_and_clear_field_errors(tmp_path):
    vor = _two_cell_grid()
    path = tmp_path / "long_wells.gpkg"
    gpd.GeoDataFrame(
        {
            "name": ["well", "well"],
            "layer": [1, 1],
            "period": [1, 2],
            "rate": [-1.0, -2.0],
        },
        geometry=[
            Polygon([(0.0, 0.0), (0.99, 0.0), (0.99, 1.0), (0.0, 1.0)]),
            Polygon([(0.0, 0.0), (0.99, 0.0), (0.99, 1.0), (0.0, 1.0)]),
        ],
        crs="EPSG:2927",
    ).to_file(path, driver="GPKG")
    source = GeoPackageSource(
        path,
        ModelContext(grid=vor),
        nper=2,
        period_field="period",
        period_base=1,
    )

    well = source.wel()

    assert well.options["stress_period_data"] == {
        0: [[(0, 0), -1.0, "well"]],
        1: [[(0, 0), -2.0, "well"]],
    }
    with pytest.raises(ValueError, match="missing fields: conductance, head"):
        source.ghb()


def test_geopackage_source_supports_cell_surface_offsets(tmp_path):
    vor = _two_cell_grid()
    path = tmp_path / "relative_boundaries.gpkg"
    _write_inputs(path)
    source = GeoPackageSource(
        path,
        ModelContext(grid=vor, domain=np.array([[1, 1]])),
        nper=1,
    )

    drn = source.drn(
        elevation=CellSurfaceOffset(reference="cell_top", offset=-2.0),
        conductance="conductance",
    )
    ghb = source.ghb(
        head=CellSurfaceOffset(
            reference="cell_bottom",
            offset="height_above_bottom",
            minimum="min_elev",
        ),
        conductance="conductance",
    )
    chd = source.chd(
        head=CellSurfaceOffset(
            reference="cell_top",
            offset="top_offset",
            minimum=8.25,
        )
    )

    assert drn.options["stress_period_data"][0] == [
        [(0, 0), 8.0, 2.0, "left"],
        [(0, 1), 7.0, 3.0, "right"],
    ]
    assert ghb.options["stress_period_data"][0] == [
        [(0, 0), 8.5, 2.0, "left"],
        [(0, 1), 8.0, 3.0, "right"],
    ]
    assert chd.options["stress_period_data"][0] == [
        [(0, 0), 9.0, "left"],
        [(0, 1), 8.5, "right"],
    ]
    assert ghb.metadata["fields"]["head"] == {
        "type": "CellSurfaceOffset",
        "reference": "cell_bottom",
        "offset": "height_above_bottom",
        "minimum": "min_elev",
    }
