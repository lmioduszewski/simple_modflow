from pathlib import Path

import geopandas as gpd
import myflopy as mf
from shapely.geometry import LineString, Point, Polygon

from myflopy.modflow.mf6.grid.triangle import TriangleGrid


def _basic_simulation() -> mf.SimulationSpec:
    gwf = mf.gwf("gwf").with_package(mf.ic(strt=9.0)).with_package(mf.npf(k=1.0))
    return mf.SimulationSpec(
        "baseline",
        models=(gwf,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )


def _write_voronoi_inputs(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame(
        {"name": ["domain"]},
        geometry=[Polygon([(0, 0), (40, 0), (40, 30), (0, 30)])],
        crs="EPSG:2927",
    ).to_file(path, layer="boundary", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["inner"], "max_area": [25.0], "priority": [2]},
        geometry=[Polygon([(8, 8), (28, 8), (28, 22), (8, 22)])],
        crs="EPSG:2927",
    ).to_file(path, layer="refinement", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["stream"], "max_area": [8.0]},
        geometry=[LineString([(3, 26), (35, 4)])],
        crs="EPSG:2927",
    ).to_file(path, layer="streams", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["well"]},
        geometry=[Point(12, 12)],
        crs="EPSG:2927",
    ).to_file(path, layer="points", driver="GPKG")
    return path


def _voronoi_grid_spec(
    path: str = "inputs.gpkg",
    *,
    include_breaklines: bool = True,
) -> mf.GridSpec:
    breaklines = (
        [
            mf.GeoPackageSourceSpec(
                path,
                layer="streams",
                fields={"area": "max_area", "label": "name"},
            )
        ]
        if include_breaklines
        else []
    )
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path, layer="boundary"),
        refinement=mf.GeoPackageSourceSpec(
            path,
            layer="refinement",
            fields={"area": "max_area", "label": "name", "priority": "priority"},
        ),
        breaklines=breaklines,
        points=[mf.GeoPackageSourceSpec(path, layer="points")],
        boundary_max_area=200.0,
        breakline_buffer=2.0,
        profile="clean",
    )


def _write_python_grid_builder(project_root: Path) -> Path:
    script = project_root / "grid" / "build_grid.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(
        "\n".join(
            [
                "from pathlib import Path",
                "from types import SimpleNamespace",
                "",
                "def build_grid(project_root, workspace, spec):",
                "    project_root = Path(project_root)",
                "    workspace = Path(workspace)",
                "    workspace.mkdir(parents=True, exist_ok=True)",
                "    marker = workspace / 'built_grid.txt'",
                "    marker.write_text(f'{spec.name}|{project_root.name}')",
                "    return SimpleNamespace(",
                "        ncpl=1,",
                "        name=spec.name,",
                "        project_root=project_root,",
                "        workspace=workspace,",
                "    )",
                "",
            ]
        )
    )
    return script


def test_project_layout_defaults_to_home_mf6(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)

    layout = mf.ProjectLayout.for_project("Elk Creek")

    assert layout.root == tmp_path / "mf6" / "Elk_Creek"
    assert layout.project_spec_path == layout.root / "specs" / "project_spec.json"
    assert (
        layout.package_spec_path("npf/base")
        == layout.root / "specs" / "packages" / "npf" / "base.json"
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
    assert loaded.mesh_options["snap_tolerance"] == 5.0
    assert loaded.triangle_options["angle"] == 30


def test_python_grid_spec_round_trip():
    grid = mf.GridSpec.python(
        "grid/build_grid.py",
        function="build_grid",
        inputs=[
            "inputs/domain.gpkg",
            mf.RasterSource("inputs/dem.tif", band=1),
        ],
        options={"profile": "project"},
        metadata={"purpose": "advanced grid"},
    )

    loaded = mf.GridSpec.from_dict(grid.to_dict())

    assert loaded.method == "python"
    assert loaded.script == "grid/build_grid.py"
    assert loaded.function == "build_grid"
    assert loaded.inputs[0] == "inputs/domain.gpkg"
    assert loaded.inputs[1].path == Path("inputs/dem.tif")
    assert loaded.options == {"profile": "project"}
    assert loaded.metadata == {"purpose": "advanced grid"}


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


def test_package_ref_round_trip_and_model_package_lookup():
    model = mf.gwf("gwf").with_package("npf/base")
    loaded = mf.ModelSpec.from_dict(model.to_dict())

    assert loaded.package_names == ("npf/base",)
    assert loaded.package("npf").key == "npf/base"
    assert loaded.package("npf/base").key == "npf/base"


def test_voronoi_grid_pickles_leanly(tmp_path):
    import pickle

    _write_voronoi_inputs(tmp_path / "inputs.gpkg")
    grid = _voronoi_grid_spec(include_breaklines=False).resolve(
        project_root=tmp_path, workspace=tmp_path / "g", build=True
    )

    # A freshly built grid carries its Triangle builder.
    assert grid.tri is not None

    restored = pickle.loads(pickle.dumps(grid))

    # The builder is dropped, but the essential mesh geometry survives.
    assert restored.tri is None
    assert restored.ncpl == grid.ncpl
    assert len(restored.iverts) == len(grid.iverts)
    assert restored.name == grid.name


def test_voronoi_grid_spec_resolves_to_triangle_setup(tmp_path):
    _write_voronoi_inputs(tmp_path / "inputs.gpkg")
    grid = _voronoi_grid_spec()

    tri = grid.resolve(
        project_root=tmp_path, workspace=tmp_path / "triangle", build=False
    )
    preview = tri.preview_regions()

    assert isinstance(tri, TriangleGrid)
    assert tri.domain_geometry.area == 1200.0
    assert set(preview["label"]) == {"domain_size", "inner", "stream"}
    assert tri._nodes is not None
    assert tri._nodes.shape == (1, 2)


def test_python_grid_spec_resolves_builder_script(tmp_path):
    _write_python_grid_builder(tmp_path)
    (tmp_path / "inputs").mkdir()
    (tmp_path / "inputs" / "domain.gpkg").write_text("placeholder")
    grid = mf.GridSpec.python(
        "grid/build_grid.py",
        function="build_grid",
        inputs=["inputs/domain.gpkg"],
    )

    resolved = grid.resolve(project_root=tmp_path, workspace=tmp_path / "grid_ws")

    assert resolved.name == "grid"
    assert resolved.project_root == tmp_path
    assert resolved.workspace == tmp_path / "grid_ws"
    marker = tmp_path / "grid_ws" / "built_grid.txt"
    assert marker.read_text() == f"grid|{tmp_path.name}"
    assert resolved.myflopy_grid_spec is grid


def test_model_spec_with_grid_resolves_into_model_context(tmp_path):
    _write_voronoi_inputs(tmp_path / "inputs.gpkg")
    model = mf.gwf("gwf").with_grid(_voronoi_grid_spec())
    simulation = mf.SimulationSpec(
        "baseline",
        models=(model,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )

    built = simulation.build_flopy(
        tmp_path / "sim",
        build_context=mf.SpecBuildContext(
            project_root=tmp_path,
            simulation_workspace=tmp_path / "sim",
            grid_workspace=tmp_path / "sim" / "_grid",
            build_grids=False,
        ),
    )

    resolved_grid = built.models["gwf"].context.grid
    assert isinstance(resolved_grid, TriangleGrid)
    assert built.model("gwf").myflopy_context.grid is resolved_grid
    assert Path(resolved_grid.model_ws) == tmp_path / "sim" / "_grid" / "gwf"


def test_manual_model_context_grid_is_preserved_without_grid_spec(tmp_path):
    manual_grid = {"kind": "manual-grid"}
    model = mf.gwf("gwf", context=mf.ModelContext(grid=manual_grid))
    simulation = mf.SimulationSpec(
        "baseline",
        models=(model,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )

    built = simulation.build_flopy(tmp_path / "sim")

    assert built.models["gwf"].context.grid is manual_grid
    assert built.model("gwf").myflopy_context.grid is manual_grid


def test_project_save_load_round_trip(tmp_path):
    baseline = _basic_simulation()
    high_k = baseline.derive("high_k").replace_package("gwf", mf.npf(k=35.0))
    project = mf.Project(tmp_path / "elk_creek", name="elk_creek")
    project.add_simulation(baseline)
    project.add_simulation(high_k)

    saved_path = project.save()
    loaded = mf.Project.load(project.root)

    assert saved_path == project.layout.project_spec_path
    assert (project.root / "project.json").exists()
    assert (project.root / "specs" / "simulations" / "baseline.json").exists()
    assert (project.root / "specs" / "simulations" / "high_k.json").exists()
    assert loaded.simulations["high_k"].derived_from == "baseline"
    assert (
        loaded.simulations["high_k"].model("gwf").package("npf").options["k"] == 35.0
    )


def test_project_saves_and_loads_package_library(tmp_path):
    project = mf.Project(tmp_path / "elk_creek", name="elk_creek")
    project.add_package("npf/base", mf.npf(k=25.0))
    project.add_package(
        "rch/big_pond", mf.PackageSpec("rch", dict, {"recharge": 0.001})
    )

    project.save()
    loaded = mf.Project.load(project.root)

    assert (project.root / "specs" / "packages" / "npf" / "base.json").exists()
    assert (project.root / "specs" / "packages" / "rch" / "big_pond.json").exists()
    assert loaded.packages["npf/base"].options["k"] == 25.0
    assert loaded.packages["rch/big_pond"].options["recharge"] == 0.001


def test_project_build_resolves_package_refs(tmp_path):
    gwf = mf.gwf("gwf").with_package("npf/base")
    simulation = mf.SimulationSpec(
        "baseline",
        models=(gwf,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )
    project = mf.Project(tmp_path / "elk_creek", name="elk_creek")
    project.add_package("npf/base", mf.npf(k=25.0))
    project.add_simulation(simulation)

    run = project.prepare_run("baseline", "baseline").build()

    assert "npf" in run.built.models["gwf"].packages
    assert run.spec.model("gwf").package("npf/base").key == "npf/base"


def test_project_build_resolves_project_relative_grid_sources(tmp_path):
    project = mf.Project(tmp_path / "elk_creek", name="elk_creek")
    _write_voronoi_inputs(project.root / "inputs" / "grid.gpkg")
    model = mf.gwf("gwf").with_grid(
        _voronoi_grid_spec("inputs/grid.gpkg", include_breaklines=False)
    )
    simulation = mf.SimulationSpec(
        "baseline",
        models=(model,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )
    project.add_simulation(simulation)

    run = project.prepare_run("baseline", "baseline").build()
    resolved_grid = run.built.models["gwf"].context.grid

    assert run.workspace == project.runs_dir / "baseline"
    assert (run.workspace / "_grid" / "gwf").exists()
    assert Path(resolved_grid.tri.model_ws) == run.workspace / "_grid" / "gwf"
    assert run.flopy_model("gwf").myflopy_context.grid is resolved_grid


def test_project_build_resolves_python_grid_builder(tmp_path):
    project = mf.Project(tmp_path / "elk_creek", name="elk_creek")
    _write_python_grid_builder(project.root)
    grid = mf.GridSpec.python("grid/build_grid.py", function="build_grid")
    model = mf.gwf("gwf").with_grid(grid)
    simulation = mf.SimulationSpec(
        "baseline",
        models=(model,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )
    project.add_simulation(simulation)

    run = project.prepare_run("baseline", "baseline").build()
    resolved_grid = run.built.models["gwf"].context.grid

    assert resolved_grid.workspace == run.workspace / "_grid" / "gwf"
    marker = run.workspace / "_grid" / "gwf" / "built_grid.txt"
    assert marker.read_text() == "grid|elk_creek"
    assert run.flopy_model("gwf").myflopy_context.grid is resolved_grid
