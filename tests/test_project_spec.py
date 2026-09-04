from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import LineString, Point, Polygon

import myflopy as mf
from myflopy.builders import build_gwf_gwt_exchange
from myflopy.modflow.mf6.grid.triangle import TriangleGrid
from myflopy.specs import ExchangeSpec, PostBuildHook


def _example_post_build_hook(model, packages, context):
    """Module-level hook callback (importable, so it serializes by reference)."""

    return {"ran": True}


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


def test_list_bc_partial_builders_round_trip(tmp_path):
    """The 7 list BCs build via a ``functools.partial(_build_named, cls)`` -- a
    picklable-but-not-importable builder. ``_callable_ref`` now serializes that
    partial (func ref + a ``$callable`` class ref), so CHD/GHB/DRN/RIV/WEL/RCH/EVT
    survive ``to_dict`` -> JSON -> ``from_dict`` and still build (plan §5.6A)."""

    import json

    bcs = {
        "chd": mf.chd(stress_period_data={0: [[(0, 0), 10.0]]}),
        "ghb": mf.ghb(stress_period_data={0: [[(0, 1), 9.0, 1.0]]}),
        "drn": mf.drn(stress_period_data={0: [[(0, 1), 9.0, 1.0]]}),
        "riv": mf.riv(stress_period_data={0: [[(0, 0), 10.0, 1.0, 9.0]]}),
        "wel": mf.wel(stress_period_data={0: [[(0, 1), -5.0]]}),
        "rch": mf.rch(stress_period_data={0: [[(0, 0), 0.001]]}),
        "evt": mf.evt(stress_period_data={0: [[(0, 0), 10.0, 0.002, 1.0]]}),
    }
    # The serialized builder is the tagged-partial shape, and it is JSON-safe.
    chd_payload = bcs["chd"].to_dict()
    assert chd_payload["builder"] == {
        "partial": "myflopy.advanced:_build_named",
        "args": [{"$callable": "flopy.mf6.modflow.mfgwfchd:ModflowGwfchd"}],
        "keywords": {},
    }
    for name, spec in bcs.items():
        loaded = mf.PackageSpec.from_dict(json.loads(json.dumps(spec.to_dict())))
        assert loaded.to_dict() == spec.to_dict(), name

    # And a full model carrying every list BC round-trips and builds real packages.
    model = mf.gwf("gwf", packages=[
        mf.disv(nlay=1, ncpl=2, nvert=6,
                vertices=[[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 1.0, 1.0], [3, 0.0, 1.0], [4, 2.0, 0.0], [5, 2.0, 1.0]],
                cell2d=[[0, 0.5, 0.5, 4, 0, 1, 2, 3], [1, 1.5, 0.5, 4, 1, 4, 5, 2]], top=10.0, botm=0.0),
        mf.ic(strt=5.0), mf.npf(k=1.0), *bcs.values(),
        mf.oc(budget_filerecord="gwf.cbc", saverecord=[("BUDGET", "ALL")]),
    ])
    sim = mf.SimulationSpec("s", models=[model],
                            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["gwf"])])
    loaded_sim = mf.SimulationSpec.from_dict(json.loads(json.dumps(sim.to_dict())))
    built = loaded_sim.build_flopy(tmp_path)
    flopy_model = built.simulation.get_model("gwf")
    for name in bcs:
        assert flopy_model.get_package(name) is not None, name
        # ... and it still HAS its rows. `get_package` is not enough: a package
        # whose period keys came back as strings exists but is empty, which is
        # how this round-tripped stably while losing every boundary (ledger 179).
        period_data = flopy_model.get_package(name).stress_period_data.get_data()
        assert list(period_data) == [0], name
        assert len(period_data[0]) == 1, name

    # The written MF6 files carry a data row, not just a MAXBOUND that claims one.
    built.simulation.write_simulation(silent=True)
    for name in bcs:
        body = (tmp_path / f"gwf.{name}").read_text()
        rows = [
            line.strip()
            for line in body.splitlines()
            if line.strip()
            and not line.strip().upper().startswith(("BEGIN", "END", "#"))
        ]
        assert any(row[0].isdigit() for row in rows), f"{name} wrote an empty period"


def test_package_ref_round_trip_and_model_package_lookup():
    model = mf.gwf("gwf").with_package("npf/base")
    loaded = mf.ModelSpec.from_dict(model.to_dict())

    assert loaded.package_names == ("npf/base",)
    assert loaded.package("npf").key == "npf/base"
    assert loaded.package("npf/base").key == "npf/base"


def test_model_spec_post_build_hooks_round_trip():
    model = mf.ModelSpec(
        "gwf",
        "gwf",
        hooks=(PostBuildHook("record", _example_post_build_hook),),
    )

    payload = model.to_dict()
    assert payload["hooks"] == [
        {
            "name": "record",
            "callback": "tests.test_project_spec:_example_post_build_hook",
        }
    ]

    loaded = mf.ModelSpec.from_dict(payload)
    assert len(loaded.hooks) == 1
    assert loaded.hooks[0].name == "record"
    # Resolved back to the very same importable function object.
    assert loaded.hooks[0].callback is _example_post_build_hook


def test_simulation_spec_exchanges_round_trip():
    gwf = mf.gwf("flow").with_package(mf.ic(strt=1.0)).with_package(mf.npf(k=1.0))
    gwt = mf.gwt("transport")
    exchange = ExchangeSpec(
        "flow-transport",
        build_gwf_gwt_exchange,
        models=("flow", "transport"),
        options={"print_flows": True},
    )
    simulation = mf.SimulationSpec(
        "coupled",
        models=(gwf, gwt),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
        exchanges=(exchange,),
    )

    payload = simulation.to_dict()
    assert payload["exchanges"] == [
        {
            "name": "flow-transport",
            "builder": "myflopy.builders:build_gwf_gwt_exchange",
            "models": ["flow", "transport"],
            "options": {"print_flows": True},
        }
    ]

    loaded = mf.SimulationSpec.from_dict(payload)
    assert len(loaded.exchanges) == 1
    restored = loaded.exchanges[0]
    assert restored.name == "flow-transport"
    assert restored.builder is build_gwf_gwt_exchange
    assert restored.models == ("flow", "transport")
    assert restored.options == {"print_flows": True}


def test_exchange_with_non_importable_builder_fails_loud():
    """A lambda/closure builder cannot serialize; to_dict must say so clearly,
    replacing the old blanket 'exchanges cannot be serialized' guard."""

    exchange = ExchangeSpec(
        "bad",
        lambda simulation, models, **options: None,
        models=("flow", "transport"),
    )
    simulation = mf.SimulationSpec(
        "coupled",
        models=(mf.gwf("flow"), mf.gwt("transport")),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
        exchanges=(exchange,),
    )

    with pytest.raises(ValueError, match="cannot be serialized"):
        simulation.to_dict()


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

    report = project.save()
    loaded = mf.Project.load(project.root)

    # save() returns a SaveReport, not a Path. It is os.fspath-able, so it can
    # still stand in for one -- but it is deliberately NOT equal to a Path.
    assert Path(report) == project.layout.project_spec_path
    assert report.simulations == 2
    assert report.skipped == ()
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

    run = project.prepare_run("baseline", "baseline")

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

    run = project.prepare_run("baseline", "baseline")
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

    run = project.prepare_run("baseline", "baseline")
    resolved_grid = run.built.models["gwf"].context.grid

    assert resolved_grid.workspace == run.workspace / "_grid" / "gwf"
    marker = run.workspace / "_grid" / "gwf" / "built_grid.txt"
    assert marker.read_text() == "grid|elk_creek"
    assert run.flopy_model("gwf").myflopy_context.grid is resolved_grid


# --- Integer stress-period keys survive serialization (ledger 179) ----------
#
# ``_json_value`` must stringify dict keys to reach JSON, and nothing reversed
# it. FloPy accepts a string period key and writes an EMPTY period block under a
# ``MAXBOUND`` that still claims a boundary, so a saved list BC came back
# silently empty and MF6 terminated normally without it.


def test_int_period_keys_restored_only_where_they_belong():
    """``_spec_options`` promotes period keys; free-form fields keep their strings."""

    from myflopy.specs import _is_int_key, _spec_options, _spec_value

    assert [_is_int_key(k) for k in ("0", "-1", "12")] == [True, True, True]
    # Non-canonical spellings are refused, so the promotion cannot merge two
    # distinct keys onto one (``"007"`` and ``"7"`` would both become ``7``).
    assert not any(
        _is_int_key(k) for k in ("007", "1_0", "1.0", " 1", "", "head.csv", 0)
    )

    restored = _spec_options(
        {
            "stress_period_data": {"0": ["a"], "2": ["b"]},
            "nested": {"0": {"1": "deep"}},
            "negative": {"0": "a", "-1": "b"},
            "padded": {"007": "a", "8": "b"},
            "continuous": {"head.csv": ["a"]},
            "empty": {},
            # Promotion is recursive through LISTS as well as dicts.
            "in_list": [{"0": "a"}, {"1": "b"}],
        }
    )
    assert sorted(restored["stress_period_data"]) == [0, 2]
    assert restored["nested"] == {0: {1: "deep"}}
    assert sorted(restored["negative"]) == [-1, 0]
    # All-or-nothing per dict: one non-canonical key leaves the whole dict alone.
    assert sorted(restored["padded"]) == ["007", "8"]
    assert list(restored["continuous"]) == ["head.csv"]
    assert restored["empty"] == {}
    assert restored["in_list"] == [{0: "a"}, {1: "b"}]
    # The kwargs mapping's own keys are FloPy parameter names -- never promoted.
    assert all(isinstance(key, str) for key in restored)
    # And the default path (metadata, lineage) is untouched.
    assert list(_spec_value({"0": "a", "1": "b"})) == ["0", "1"]


def test_package_metadata_keeps_string_keys_through_round_trip():
    """Provenance is free-form user data, so ``"0"`` must not become ``0``."""

    package = mf.npf(k=1.0).with_metadata(**{"0": "zero", "source": "gis"})
    loaded = mf.PackageSpec.from_dict(package.to_dict())

    assert loaded.metadata == {"0": "zero", "source": "gis"}


def test_project_save_load_preserves_stress_period_rows(tmp_path):
    """The reported failure end to end: save -> load -> the BC still has its rows."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package(
        "chd/base",
        mf.chd(stress_period_data={0: [[(0, 0), 10.0]], 3: [[(0, 1), 9.0]]}),
    )
    project.save()

    reloaded = mf.Project.load(tmp_path / "proj")
    period_data = reloaded.packages["chd/base"].options["stress_period_data"]

    assert sorted(period_data) == [0, 3]
    assert all(isinstance(key, int) for key in period_data)


def test_int_period_keys_restored_through_inputs_and_bound_builder_kwargs():
    """The other two channels carrying builder kwargs, both easy to leave behind.

    ``PackageSpec.build`` splats ``**inputs, **options`` into one call, and a
    ``functools.partial`` builder can carry bound keywords (the §5.6A list-BC
    machinery). All three are the same semantic class, so all three restore.
    """

    import functools

    import flopy

    from myflopy.advanced import _build_named

    period = {0: [[(0, 0), 10.0]]}

    via_inputs = mf.PackageSpec(
        "chd", flopy.mf6.ModflowGwfchd, inputs={"stress_period_data": period}
    )
    assert list(
        mf.PackageSpec.from_dict(via_inputs.to_dict()).inputs["stress_period_data"]
    ) == [0]

    via_partial = mf.PackageSpec(
        "chd",
        functools.partial(
            _build_named, flopy.mf6.ModflowGwfchd, stress_period_data=period
        ),
    )
    revived = mf.PackageSpec.from_dict(via_partial.to_dict())
    assert list(revived.builder.keywords["stress_period_data"]) == [0]


def test_from_dict_still_rejects_a_null_options_block():
    """An absent key defaults; an explicit ``null`` must stay a loud ``TypeError``.

    A hand-edited YAML with an emptied ``options:`` would otherwise build a
    package with no period data at all -- the same "exists but is empty" outcome
    the int-key restoration exists to prevent.
    """

    base = {
        "kind": "PackageSpec",
        "name": "chd",
        "builder": "flopy.mf6.modflow.mfgwfchd:ModflowGwfchd",
    }
    assert mf.PackageSpec.from_dict(base).options == {}
    for field in ("options", "inputs", "metadata", "requires"):
        with pytest.raises(TypeError):
            mf.PackageSpec.from_dict({**base, field: None})


def test_source_payloads_are_rebuilt_by_registry_and_keep_int_keys():
    """Two defects in the same guard, both silent.

    ``_spec_value`` decided "is this a source?" with ``kind.endswith("Source")``.
    ``DataSourceSpec``'s own kind ends in "Spec", so the public base class came
    back as a plain ``dict``. And the guard returned BEFORE the int-key
    restoration, so a period-keyed ``LiteralSource`` lost its keys through JSON
    while YAML kept them -- the two documented formats disagreeing about one spec.
    """

    import json

    import yaml

    from myflopy.sources import DataSourceSpec, LiteralSource, RasterSource
    from myflopy.specs import _json_value, _spec_options, _spec_value

    for source in (DataSourceSpec(path="x.tif"), RasterSource(path="r.tif")):
        revived = mf.PackageSpec.from_dict(
            json.loads(json.dumps(mf.npf(k=1.0).with_input("k", source).to_dict()))
        ).inputs["k"]
        assert type(revived) is type(source)

    # An unknown "...Source" kind must still fail loudly rather than degrade.
    with pytest.raises(ValueError, match="Unknown source kind"):
        _spec_value({"kind": "MadeUpSource"})

    payload = _json_value({"steady_state": LiteralSource({0: True, 1: False})})
    from_json = _spec_options(json.loads(json.dumps(payload)))["steady_state"]
    from_yaml = _spec_options(yaml.safe_load(yaml.safe_dump(payload)))["steady_state"]
    assert from_json.value == from_yaml.value == {0: True, 1: False}


def test_exchange_options_restore_int_period_keys():
    """``ModflowGwfgwf.perioddata`` is an exchange-level MVR, so it IS period-keyed.

    Latent today -- FloPy 3.10 cannot construct the exchange child package at all
    -- but the round trip lost the keys, so a fixed FloPy would have shipped a
    silently no-op inter-model mover.
    """

    import json

    from myflopy.builders import build_gwf_gwf_exchange

    exchange = ExchangeSpec(
        name="gwfgwf",
        builder=build_gwf_gwf_exchange,
        models=("a", "b"),
        options={"nexg": 1, "perioddata": {0: [["a_sfr", 0, "b_lak", 0, "FACTOR", 1.0]]}},
    )
    revived = ExchangeSpec.from_dict(json.loads(json.dumps(exchange.to_dict())))

    assert list(revived.options["perioddata"]) == [0]


def test_tuples_survive_serialization_and_observations_write_real_cellids(tmp_path):
    """Tuples are tagged, not flattened -- the difference reaches the MF6 files.

    FloPy re-tuples cellids into a ``stress_period_data`` recarray, which hid
    this; ``ModflowUtlobs`` instead writes each record field with ``repr``, so an
    untagged cellid was emitted as ``west chd [0, 0, 0]``. MF6 parses that,
    terminates normally, and returns DNODATA for every observation (ledger 179).
    """

    import json

    def build() -> mf.SimulationSpec:
        model = mf.gwf(
            "gwf",
            packages=[
                mf.dis(nlay=1, nrow=1, ncol=3, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                mf.ic(strt=5.0),
                mf.npf(k=1.0),
                mf.chd(
                    stress_period_data={0: [[(0, 0, 0), 10.0], [(0, 0, 2), 8.0]]},
                    observations={
                        "chd_obs.csv": [("west", "chd", (0, 0, 0)), ("east", "chd", (0, 0, 2))]
                    },
                ),
            ],
        )
        return mf.SimulationSpec(
            "s",
            models=[model],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["gwf"])],
        )

    original = build()
    # The round trip is now a true fixed point, tuples included.
    assert mf.SimulationSpec.from_dict(
        json.loads(json.dumps(original.to_dict()))
    ).to_dict() == original.to_dict()

    revivals = {
        "direct": build(),
        "json": mf.SimulationSpec.from_dict(json.loads(json.dumps(build().to_dict()))),
        "yaml": mf.SimulationSpec.from_yaml(build().to_yaml()),
    }
    for label, spec in revivals.items():
        workspace = tmp_path / label
        built = spec.build_flopy(workspace)
        built.simulation.write_simulation(silent=True)
        obs_file = next(workspace.rglob("*.obs"))
        rows = [
            line.split()
            for line in obs_file.read_text().splitlines()
            if line.strip().startswith(("west", "east"))
        ]
        assert rows == [["west", "chd", "1", "1", "1"], ["east", "chd", "1", "1", "3"]], label


def test_grid_spec_tuple_field_is_not_double_tagged():
    """``GridSpec.inputs`` is a tuple FIELD whose serialized shape is a plain list.

    Tagging the field itself would make ``from_dict`` walk ``{"$tuple": [...]}``
    as a mapping and rebuild the tuple as ``("$tuple",)``.
    """

    grid = mf.GridSpec.python(
        "grid/build_grid.py",
        function="build_grid",
        inputs=["inputs/domain.gpkg", mf.RasterSource("inputs/dem.tif", band=1)],
    )
    payload = grid.to_dict()

    assert isinstance(payload["inputs"], list)
    assert mf.GridSpec.from_dict(payload).inputs[0] == "inputs/domain.gpkg"
