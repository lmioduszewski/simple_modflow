from __future__ import annotations

import json
import shutil

import flopy
import pytest

import myflopy as mf
from myflopy import ModelSpec, PackageSpec, Project, Run, SimulationSpec, build_ims
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.workspace import ModelView


def _tiny_flow_spec() -> SimulationSpec:
    """Return a complete one-cell steady-state GWF simulation."""

    flow = ModelSpec(
        "flow",
        "gwf",
        packages=(
            PackageSpec(
                "dis",
                flopy.mf6.ModflowGwfdis,
                {
                    "nlay": 1,
                    "nrow": 1,
                    "ncol": 1,
                    "delr": 1.0,
                    "delc": 1.0,
                    "top": 10.0,
                    "botm": 0.0,
                },
            ),
            PackageSpec("ic", flopy.mf6.ModflowGwfic, {"strt": 9.0}),
            PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0}),
            PackageSpec(
                "chd",
                flopy.mf6.ModflowGwfchd,
                {"stress_period_data": {0: [[(0, 0, 0), 9.0]]}},
            ),
            PackageSpec(
                "oc",
                flopy.mf6.ModflowGwfoc,
                {
                    "head_filerecord": "flow.hds",
                    "budget_filerecord": "flow.cbc",
                    "saverecord": [("HEAD", "ALL"), ("BUDGET", "ALL")],
                },
            ),
        ),
    )
    return SimulationSpec(
        "tiny_flow",
        models=(flow,),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
            PackageSpec(
                "ims",
                build_ims,
                {"models": ("flow",), "complexity": "SIMPLE"},
            ),
        ),
    )


def _tiny_flow_with_npf_ref() -> SimulationSpec:
    """Return the tiny flow simulation with NPF supplied as a project ref."""

    flow = ModelSpec(
        "flow",
        "gwf",
        packages=(
            PackageSpec(
                "dis",
                flopy.mf6.ModflowGwfdis,
                {
                    "nlay": 1,
                    "nrow": 1,
                    "ncol": 1,
                    "delr": 1.0,
                    "delc": 1.0,
                    "top": 10.0,
                    "botm": 0.0,
                },
            ),
            PackageSpec("ic", flopy.mf6.ModflowGwfic, {"strt": 9.0}),
            mf.ref("npf/base"),
            PackageSpec(
                "chd",
                flopy.mf6.ModflowGwfchd,
                {"stress_period_data": {0: [[(0, 0, 0), 9.0]]}},
            ),
            PackageSpec(
                "oc",
                flopy.mf6.ModflowGwfoc,
                {"head_filerecord": "flow.hds", "saverecord": [("HEAD", "ALL")]},
            ),
        ),
    )
    return SimulationSpec(
        "tiny_flow",
        models=(flow,),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
            PackageSpec(
                "ims", build_ims, {"models": ("flow",), "complexity": "SIMPLE"}
            ),
        ),
    )


def test_project_package_library_resolves_and_swaps_references(tmp_path):
    import numpy as np

    project = Project(tmp_path / "demo", name="demo")
    project.add_package(
        "npf/base", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0})
    )
    project.add_package(
        "npf/high_k", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 100.0})
    )

    baseline = project.add_simulation(_tiny_flow_with_npf_ref())
    high_k = baseline.derive("high_k").replace_package("flow", mf.ref("npf/high_k"))
    project.add_simulation(high_k)

    base_run = project.prepare_run("baseline", "tiny_flow").build()
    hk_run = project.prepare_run("high_k", "high_k").build()

    def npf_k(run: Run) -> float:
        npf = run.built.built_model("flow").package("npf")
        return float(np.asarray(npf.k.array).reshape(-1)[0])

    # Same reference machinery, different library entry per simulation.
    assert npf_k(base_run) == 1.0
    assert npf_k(hk_run) == 100.0


def test_project_unresolved_package_reference_raises(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    # Library is empty: the "npf/base" reference cannot be resolved.
    project.add_simulation(_tiny_flow_with_npf_ref())

    run = project.prepare_run("baseline", "tiny_flow")
    with pytest.raises((KeyError, ValueError)):
        run.build()


def test_project_save_rejects_unresolved_references(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    # The simulation references "npf/base", which is not in the library.
    project.add_simulation(_tiny_flow_with_npf_ref())

    issues = project.validate()
    assert any("npf/base" in issue for issue in issues)
    with pytest.raises(ValueError, match="npf/base"):
        project.save()

    # Once defined, the project is valid and saves.
    project.add_package(
        "npf/base", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0})
    )
    assert project.validate() == []
    project.save()


def test_provenance_keeps_package_reference_variant(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    project.add_package(
        "npf/base", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0})
    )
    project.add_package(
        "npf/high_k", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 50.0})
    )

    high_k = (
        _tiny_flow_with_npf_ref()
        .derive("high_k")
        .replace_package("flow", mf.ref("npf/high_k"))
    )

    # Lineage records the library key, not just the package name.
    assert high_k.lineage[-1] == {
        "operation": "replace_package",
        "model": "flow",
        "package": "npf/high_k",
    }

    run = project.prepare_run("high_k", high_k)
    manifest = json.loads(run.manifest_path.read_text(encoding="utf-8"))
    flow_packages = manifest["spec"]["models"][0]["packages"]
    assert "npf/high_k" in flow_packages
    assert "npf" not in flow_packages


def test_project_owns_reusable_specs_and_prepares_runs(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    simulation = project.add_simulation(_tiny_flow_spec())
    package = project.add_package("npf/base", simulation.models[0].packages[2])

    run = project.prepare_run("baseline", simulation, metadata={"scenario": "baseline"})

    assert project.simulations["tiny_flow"] is simulation
    assert project.packages["npf/base"] is package
    assert run.workspace == project.runs_dir / "baseline"
    assert run.status == "created"
    assert run.manifest_path.exists()

    # Persistence is opt-in: nothing on disk until save().
    assert not project.layout.project_spec_path.exists()
    project.save()
    reloaded = Project.load(project.root)
    assert reloaded.name == "demo"
    assert set(reloaded.simulations) == {"tiny_flow"}
    assert set(reloaded.packages) == {"npf/base"}


def test_run_builds_writes_executes_and_reopens_with_flopy_310(tmp_path):
    if shutil.which("mf6") is None:
        pytest.skip("mf6 executable is required for the run lifecycle integration test")

    project = Project(tmp_path / "demo", name="demo")
    project.add_simulation(_tiny_flow_spec())

    run = project.run("baseline", "tiny_flow")

    assert run.success is True
    assert run.status == "completed"
    assert (run.workspace / "mfsim.nam").exists()
    assert (run.workspace / "flow.hds").exists()

    manifest = json.loads(run.manifest_path.read_text(encoding="utf-8"))
    assert manifest["run"]["status"] == "completed"
    assert manifest["spec"]["models"][0]["packages"] == ["dis", "ic", "npf", "chd", "oc"]
    assert manifest["spec"]["models"][0]["package_metadata"] == {}
    assert manifest["spec"]["models"][0]["context_metadata"] == {}
    assert manifest["spec"]["models"][0]["post_build_hooks"] == []

    reopened = project.reopen_run("baseline")

    assert reopened.spec is None
    assert reopened.status == "completed"
    assert reopened.success is True
    assert reopened.simulation.get_model("flow").name == "flow"
    assert [item.name for item in project.discover_runs()] == ["baseline"]
    assert [item.name for item in project.discover_native_runs()] == ["baseline"]


def test_simulation_spec_build_returns_run_with_preferred_model_view(tmp_path):
    simulation = _tiny_flow_spec().with_workspace(tmp_path / "direct")

    run = simulation.build()
    model = run.model("flow")
    default_model = run.model()

    assert isinstance(run, Run)
    assert run.status == "built"
    assert run.workspace == tmp_path / "direct"
    assert run.model_names == ("flow",)
    assert isinstance(model, ModelView)
    assert default_model is model
    assert model.name == "flow"
    assert model.gwf is run.built.gwf("flow")
    assert run.flopy_model() is model.gwf
    assert model.sim is run.simulation
    assert model.sim.plot.__func__.__name__ == "_simulation_plot_compat"
    assert model.package("npf").parent.name == "flow"
    assert model.regions.model is model

    base_model = SimulationBase.from_built_run(run, "flow")

    assert isinstance(base_model, SimulationBase)
    assert base_model.name == "flow"
    assert base_model.gwf is model.gwf
    assert base_model.sim is run.simulation
    assert base_model.regions.model is base_model


def test_native_workspace_load_returns_file_backed_model_view(tmp_path):
    simulation = _tiny_flow_spec().with_workspace(tmp_path / "native")
    run = simulation.build()
    run.write()

    loaded = mf.load_run(run.workspace)
    model = loaded.model("flow")
    default_model = loaded.model()

    assert loaded.status == "written"
    assert loaded.model_names == ("flow",)
    assert default_model is model
    assert model.name == "flow"
    assert model.gwf.name == "flow"
    assert loaded.flopy_model("flow").name == "flow"
    assert model.package("npf").parent.name == "flow"


def test_run_model_method_is_typed_for_ide_completion():
    assert Run.model.__annotations__["return"] == "SimulationBase"


def test_reopened_run_without_spec_cannot_be_rebuilt(tmp_path):
    project = Project(tmp_path / "demo")
    run = project.prepare_run("baseline", _tiny_flow_spec())
    reopened = project.reopen_run("baseline", load=False)

    with pytest.raises(ValueError, match="without a spec cannot be rebuilt"):
        reopened.build()


def _one_cell_gwf(name: str) -> ModelSpec:
    """Return a minimal one-cell GWF model spec."""

    return ModelSpec(
        name,
        "gwf",
        packages=(
            PackageSpec(
                "dis",
                flopy.mf6.ModflowGwfdis,
                {
                    "nlay": 1,
                    "nrow": 1,
                    "ncol": 1,
                    "delr": 1.0,
                    "delc": 1.0,
                    "top": 10.0,
                    "botm": 0.0,
                },
            ),
            PackageSpec("ic", flopy.mf6.ModflowGwfic, {"strt": 9.0}),
            PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0}),
        ),
    )


def _two_model_sim() -> SimulationSpec:
    """Return a two-model simulation (separate solvers, no exchange)."""

    return SimulationSpec(
        "coupled",
        models=(_one_cell_gwf("north"), _one_cell_gwf("south")),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
            PackageSpec("ims_north", build_ims, {"models": ("north",), "complexity": "SIMPLE"}),
            PackageSpec("ims_south", build_ims, {"models": ("south",), "complexity": "SIMPLE"}),
        ),
    )


def test_project_multi_model_run_uses_per_model_subdirs(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    project.add_simulation(_two_model_sim())

    run = project.prepare_run("high_k", "coupled")
    run.write()

    # mfsim.nam stays at the run root; each model writes into its own subdir.
    assert (run.workspace / "mfsim.nam").exists()
    assert (run.workspace / "north" / "north.dis").exists()
    assert (run.workspace / "south" / "south.dis").exists()
    assert run.spec.model("north").options["model_rel_path"] == "north"
    assert run.spec.model("south").options["model_rel_path"] == "south"


def test_project_single_model_run_stays_flat(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    project.add_simulation(_tiny_flow_spec())

    run = project.prepare_run("baseline", "tiny_flow")
    run.write()

    # One model: no subdirectory, files at the run root.
    assert (run.workspace / "flow.dis").exists()
    assert not (run.workspace / "flow").is_dir()
    assert run.spec.model("flow").options.get("model_rel_path", ".") == "."


def test_two_solver_simulation_writes_separate_solutions(tmp_path):
    # Two models with two solvers must produce two solutions in mfsim.nam.
    # Regression: IMS used to build before models and collapse into one
    # solution, which MF6 rejects for coupled simulations.
    simulation = SimulationSpec(
        "multi",
        models=(_one_cell_gwf("north"), _one_cell_gwf("south")),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
            mf.ims(name="ims_north", models=["north"], complexity="SIMPLE"),
            mf.ims(name="ims_south", models=["south"], complexity="SIMPLE"),
        ),
    ).with_workspace(tmp_path / "multi")

    run = simulation.build()
    run.write()

    nam = (run.workspace / "mfsim.nam").read_text()
    solutions = nam[nam.find("BEGIN solutiongroup") :]
    assert solutions.count("ims6") == 2
    assert "ims_north.ims" in solutions and "ims_south.ims" in solutions
    # The first model's solver is listed first.
    assert solutions.index("ims_north.ims") < solutions.index("ims_south.ims")


def test_project_pickles_array_bearing_packages_and_reuses_them(tmp_path):
    import numpy as np

    project = Project(tmp_path / "demo", name="demo")
    # Scalar package: JSON-serializable.
    project.add_package(
        "npf/scalar", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 10.0})
    )
    # Computed-array package: pickled as an artifact.
    karr = np.full((1, 1, 1), 7.0)
    project.add_package(
        "npf/array", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": karr})
    )

    project.save()
    pkgs = project.layout.package_specs_dir

    # Scalar -> JSON only; array -> JSON pointer + pickle + version sidecar.
    assert (pkgs / "npf" / "scalar.json").exists()
    assert not (pkgs / "npf" / "scalar.pkl").exists()
    assert (pkgs / "npf" / "array.json").exists()
    assert (pkgs / "npf" / "array.pkl").exists()
    assert (pkgs / "npf" / "array.versions.json").exists()

    # Reload in a fresh project: both packages come back, array intact.
    loaded = Project.load(project.root)
    assert loaded.packages["npf/scalar"].options["k"] == 10.0
    reused = np.asarray(loaded.packages["npf/array"].options["k"])
    assert float(reused.reshape(-1)[0]) == 7.0
