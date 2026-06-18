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

    reopened_project = Project.reopen(project.root)
    assert reopened_project.name == "demo"
    assert reopened_project.simulations == {}


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
