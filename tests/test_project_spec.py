from pathlib import Path

import geopandas as gpd
import myflopy as mf
from shapely.geometry import LineString, Point, Polygon

from myflopy.modflow.mf6.grid.triangle import TriangleGrid


def _basic_simulation() -> mf.SimulationSpec:
    grid = mf.GridSpec.structured(
        name="small_grid",
        nlay=1,
        nrow=1,
        ncol=1,
        delr=100.0,
        delc=100.0,
        top=10.0,
        botm=0.0,
    )
    gwf = (
        mf.gwf("gwf")
        .with_grid(grid)
        .with_package(mf.ic(strt=9.0))
        .with_package(mf.npf(k=1.0))
    )
    return mf.SimulationSpec(
        "baseline",
        models=(gwf,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )


def test_project_layout_defaults_to_home_mf6(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    layout = mf.ProjectLayout.for_project("Elk Creek")

    assert layout.root == tmp_path / "mf6" / "Elk_Creek"
    assert layout.project_spec_path == layout.root / "specs" / "project_spec.json"
    assert (
        layout.simulation_workspace("baseline")
        == layout.root / "simulations" / "baseline"
    )


def test_simulation_spec_with_models_and_derivation_replace_package():
    baseline = _basic_simulation()

    high_k = baseline.derive("high_k").replace_package("gwf", mf.npf(k=35.0))

    assert high_k.name == "high_k"
    assert high_k.derived_from == "baseline"
    assert high_k.model("gwf").package("npf").options["k"] == 35.0
    assert baseline.model("gwf").package("npf").options["k"] == 1.0
    assert high_k.lineage[-1] == {
        "operation": "replace_package",
        "model": "gwf",
        "package": "npf",
    }


def test_grid_and_source_specs_round_trip():
    grid = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec("inputs/domain.gpkg", layer="boundary"),
        refinement=mf.GeoPackageSourceSpec(
            "inputs/refinement.gpkg",
            layer="zones",
            fields={"area": "max_area"},
        ),
        breaklines=[
            mf.ShapeSource("inputs/streams.shp", crs="EPSG:26915"),
        ],
        points=[
            mf.TableSource("inputs/wells.csv"),
        ],
        crs="EPSG:26915",
        snap_tolerance=5.0,
        min_angle=30,
    )

    loaded = mf.GridSpec.from_dict(grid.to_dict())

    assert loaded.method == "voronoi"
    assert loaded.engine == "triangle_voronoi_plus"
    assert loaded.boundary.layer == "boundary"
    assert loaded.refinement.fields == {"area": "max_area"}
    assert loaded.breaklines[0].path == Path("inputs/streams.shp")
    assert loaded.options["snap_tolerance"] == 5.0


def test_package_spec_inputs_round_trip_separately_from_options():
    package = mf.PackageSpec(
        "rch",
        dict,
        options={"fixed_cell": True},
        inputs={
            "recharge": mf.RasterSource(
                "inputs/recharge.tif",
                band=1,
                map_to="grid",
                method="mean",
            )
        },
    )

    loaded = mf.PackageSpec.from_dict(package.to_dict())

    assert loaded.options == {"fixed_cell": True}
    assert loaded.inputs["recharge"].path == Path("inputs/recharge.tif")
    assert loaded.inputs["recharge"].map_to == "grid"


def test_voronoi_grid_spec_resolves_to_triangle_setup(tmp_path):
    inputs = tmp_path / "inputs.gpkg"
    gpd.GeoDataFrame(
        {"name": ["domain"]},
        geometry=[Polygon([(0, 0), (40, 0), (40, 30), (0, 30)])],
        crs="EPSG:2927",
    ).to_file(inputs, layer="boundary", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["inner"], "max_area": [25.0], "priority": [2]},
        geometry=[Polygon([(8, 8), (28, 8), (28, 22), (8, 22)])],
        crs="EPSG:2927",
    ).to_file(inputs, layer="refinement", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["stream"], "max_area": [8.0]},
        geometry=[LineString([(3, 26), (35, 4)])],
        crs="EPSG:2927",
    ).to_file(inputs, layer="streams", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["well"]},
        geometry=[Point(12, 12)],
        crs="EPSG:2927",
    ).to_file(inputs, layer="points", driver="GPKG")

    grid = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec("inputs.gpkg", layer="boundary"),
        refinement=mf.GeoPackageSourceSpec(
            "inputs.gpkg",
            layer="refinement",
            fields={"area": "max_area", "label": "name", "priority": "priority"},
        ),
        breaklines=[
            mf.GeoPackageSourceSpec(
                "inputs.gpkg",
                layer="streams",
                fields={"area": "max_area", "label": "name"},
            )
        ],
        points=[mf.GeoPackageSourceSpec("inputs.gpkg", layer="points")],
        boundary_max_area=200.0,
        breakline_buffer=2.0,
        profile="clean",
    )

    tri = grid.resolve(
        project_root=tmp_path, workspace=tmp_path / "triangle", build=False
    )
    preview = tri.preview_regions()

    assert isinstance(tri, TriangleGrid)
    assert tri.domain_geometry.area == 1200.0
    assert set(preview["label"]) == {"domain_size", "inner", "stream"}
    assert tri._nodes is not None
    assert tri._nodes.shape == (1, 2)


def test_project_spec_save_load_round_trip(tmp_path):
    baseline = _basic_simulation()
    high_k = baseline.derive("high_k").replace_package("gwf", mf.npf(k=35.0))
    project = mf.ProjectSpec(
        "elk_creek",
        root=tmp_path / "elk_creek",
    ).with_simulations(baseline, high_k)

    saved_path = project.save()
    loaded = mf.ProjectSpec.load(project.layout.root)

    assert saved_path == project.layout.project_spec_path
    assert (project.layout.root / "project.json").exists()
    assert (project.layout.root / "specs" / "simulations" / "baseline.json").exists()
    assert (project.layout.root / "specs" / "simulations" / "high_k.json").exists()
    assert not (project.layout.root / "snapshots").exists()
    assert loaded.simulation("high_k").derived_from == "baseline"
    assert loaded.simulation("high_k").model("gwf").package("npf").options["k"] == 35.0


def test_project_spec_build_uses_simulations_workspace(tmp_path):
    project = mf.ProjectSpec(
        "elk_creek",
        root=tmp_path / "elk_creek",
    ).with_simulation(_basic_simulation())

    run = project.build("baseline")

    assert run.workspace == project.layout.root / "simulations" / "baseline"
    assert run.simulation is not None
    assert run.flopy_model("gwf") is not None
    assert run.spec.model("gwf").options["model_rel_path"] == "gwf"
