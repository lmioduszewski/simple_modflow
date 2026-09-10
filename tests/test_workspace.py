from __future__ import annotations

import json
import shutil

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy import (
    ModelSpec,
    PackageRef,
    PackageSpec,
    Project,
    Run,
    SimulationSpec,
    build_ims,
)
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
    high_k = baseline.derive("high_k").with_packages(mf.ref("npf/high_k"), on="flow")
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

    # prepare_run builds by default, so the unresolved reference raises there.
    with pytest.raises((KeyError, ValueError)):
        project.prepare_run("baseline", "tiny_flow")


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


def test_project_custom_layout_round_trips(tmp_path):
    layout = mf.ProjectLayout(
        tmp_path / "demo", specs_dir_name="recipes", packages_dir_name="pkgs"
    )
    project = Project(tmp_path / "demo", name="demo", layout=layout)
    project.add_package(
        "npf/base", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0})
    )
    project.save()

    assert (project.root / "recipes" / "pkgs" / "npf" / "base.json").exists()

    # Load without passing a layout: the saved custom dir names are restored.
    loaded = Project.load(project.root)
    assert loaded.layout.specs_dir_name == "recipes"
    assert loaded.layout.packages_dir_name == "pkgs"
    assert loaded.packages["npf/base"].options["k"] == 1.0


def test_project_load_fails_loud_on_missing_file(tmp_path):
    project = Project(tmp_path / "demo", name="demo")
    project.add_package(
        "npf/base", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0})
    )
    project.save()
    project.layout.package_spec_path("npf/base").unlink()

    with pytest.raises(FileNotFoundError, match="npf/base"):
        Project.load(project.root)


def test_add_package_accepts_semantic_keys(tmp_path):
    # A key is a free-form library address; the package keeps its own name and is
    # addressed by the full key, so a semantic key whose prefix is not the MF6
    # package type is valid (no warning, no coupling).
    project = Project(tmp_path / "demo", name="demo")

    spec = project.add_package(
        "k/calibrated", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0})
    )
    assert spec.name == "npf"
    assert project.packages["k/calibrated"] is spec


def test_replace_package_targets_slot_for_semantic_key_ref(tmp_path):
    # A semantically-keyed ref ("k/high" -> builds npf) replaces the npf slot
    # when the slot is named explicitly.
    project = Project(tmp_path / "demo", name="demo")
    project.add_package("k/high", mf.npf(k=99.0))

    variant = (
        _tiny_flow_spec()
        .derive("v")
        .with_packages(mf.ref("k/high"), on="flow", slot="npf")
    )
    assert variant.model("flow").package("k/high").key == "k/high"
    # Targeting a non-existent slot is a clear error.
    with pytest.raises(KeyError, match="no 'does_not_exist' package"):
        _tiny_flow_spec().with_packages(mf.npf(k=1.0), on="flow", slot="does_not_exist")


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
        .with_packages(mf.ref("npf/high_k"), on="flow")
    )

    # Lineage records the library key, not just the package name.
    assert high_k.lineage[-1] == {
        "operation": "with_packages",
        "model": "flow",
        "packages": ["npf/high_k"],
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
    assert run.status == "built"  # prepare_run builds in memory by default
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

    # The manifest stores the full recipe, so a reopened run gets a real spec back.
    assert reopened.spec is not None
    # the SPEC's name is the simulation's; "baseline" is the RUN's
    assert reopened.spec.name == manifest["spec"]["name"]
    assert reopened.spec.model_names == ("flow",)
    # ...and the readable summary block is retained too, so saving cannot erase it.
    assert reopened.spec_block == manifest["spec"]
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


def test_a_run_whose_manifest_has_no_recipe_cannot_be_rebuilt(tmp_path):
    """Manifests written before recipes existed -- and runs whose spec embedded a
    raw array -- carry a readable summary only. Those stay summary-only: openable
    and inspectable, but not rebuildable, which is the pre-recipe status quo."""

    project = Project(tmp_path / "demo")
    run = project.prepare_run("baseline", _tiny_flow_spec())
    manifest = json.loads(run.manifest_path.read_text())
    manifest.pop("spec_recipe")
    run.manifest_path.write_text(json.dumps(manifest))

    reopened = project.reopen_run("baseline", load=False)

    assert reopened.spec is None
    assert reopened.spec_block is not None       # it still knows what it WAS
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


def test_reopening_a_run_does_not_erase_its_provenance(tmp_path):
    """A reopened run has no ``SimulationSpec``, and used to write that ignorance
    to disk: ``to_dict`` emitted ``"spec": null`` over the record of what the run
    contained. ``_record_status`` saves on every state change, so merely
    re-executing a reopened run destroyed its own history.
    """

    project = mf.Project(tmp_path / "proj", name="p")
    simulation = _tiny_flow_spec()
    run = project.prepare_run("baseline", simulation)
    original = json.loads(run.manifest_path.read_text())["spec"]
    assert original["models"][0]["packages"]

    reopened = project.reopen_run("baseline", load=False)
    reopened.save_manifest()

    assert json.loads(run.manifest_path.read_text())["spec"] == original

    # Every status transition re-writes the whole manifest, so the block has to
    # survive a sequence, not just one save.
    for status in ("built", "written", "completed"):
        reopened._record_status(status)
    persisted = json.loads(run.manifest_path.read_text())
    assert persisted["spec"] == original
    assert persisted["run"]["status"] == "completed"


def test_a_live_spec_still_wins_over_a_carried_forward_block(tmp_path):
    """The carry-forward must not shadow a real spec: a run that HAS one
    re-summarizes it, so an edited simulation is recorded as edited."""

    project = mf.Project(tmp_path / "proj", name="p")
    run = project.prepare_run("baseline", _tiny_flow_spec())
    stale = json.loads(run.manifest_path.read_text())["spec"]

    # The edit adds a package to the MODEL. It used to hand an `oc` to
    # `with_package` at the SIMULATION tier -- which FloPy rejects at build
    # ("Invalid package type"), invisible here because this test only writes a
    # manifest, and hidden further by an assertion that accepted the package in
    # either tier. The simulation tier now refuses it at declaration.
    # `sto`, not `oc`: the model already carries an `oc`, and the manifest
    # summary records package NAMES, so replacing one changes nothing to assert.
    edited = _tiny_flow_spec()
    edited = edited.with_model(
        edited.model("flow").with_packages(mf.sto(steady_state={0: True}))
    )
    reopened = project.reopen_run("baseline", load=False)
    reopened.spec = edited
    reopened.save_manifest()

    written = json.loads(run.manifest_path.read_text())["spec"]
    assert written != stale
    assert any("sto" in model["packages"] for model in written["models"])


def test_a_native_workspace_carries_no_spec_block(tmp_path):
    """``Run.load`` on a manifest-less MF6 workspace has nothing to carry."""

    project = mf.Project(tmp_path / "proj", name="p")
    run = project.prepare_run("baseline", _tiny_flow_spec())
    run.write()
    run.manifest_path.unlink()

    native = mf.Run.load(run.workspace, load=False)

    assert native.spec_block is None
    assert native.to_dict()["spec"] is None


def test_one_unserializable_simulation_no_longer_costs_the_whole_library(tmp_path):
    """``save()`` used to refuse everything when any simulation embedded an array.

    A project with a full package library and one such simulation therefore wrote
    NOTHING -- measured: 0 package files -- and the library did not survive the
    session even though every package would have persisted on its own.
    """

    project = mf.Project(tmp_path / "proj", name="p")
    array_package = mf.npf(k=np.arange(9.0))
    project.add_package("npf/base", array_package)
    project.add_package("ic/steady", mf.ic(strt=1.0))
    project.add_simulation(
        mf.SimulationSpec(
            "good",
            models=[mf.gwf("flow", packages=[mf.ref("npf/base"), mf.ref("ic/steady")])],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )
    # embeds the array directly rather than referencing the library
    project.add_simulation(
        mf.SimulationSpec(
            "inline",
            models=[mf.gwf("flow", packages=[array_package])],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )

    report = project.save()

    assert report.packages == 2
    assert report.pickled == 1              # the array package took the pickle path
    assert report.simulations == 1
    assert [name for name, _ in report.skipped] == ["inline"]
    assert bool(report) is True

    loaded = mf.Project.load(project.root)
    assert sorted(loaded.packages) == ["ic/steady", "npf/base"]
    assert np.array_equal(loaded.packages["npf/base"].options["k"], np.arange(9.0))
    # The skipped simulation is left OUT of the listing, not listed and absent --
    # load() raises FileNotFoundError for a listed entry whose file is missing.
    assert sorted(loaded.simulations) == ["good"]


def test_strict_save_still_refuses_everything(tmp_path):
    """``strict=True`` is the old behaviour, and must raise BEFORE writing."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_simulation(
        mf.SimulationSpec(
            "inline",
            models=[mf.gwf("flow", packages=[mf.npf(k=np.arange(9.0))])],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )

    with pytest.raises(ValueError, match="cannot be saved"):
        project.save(strict=True)

    assert not project.layout.project_spec_path.exists()


def test_a_dangling_reference_still_raises_but_after_the_write(tmp_path):
    """An embedded array is nothing to fix; a dangling ref is a mistake with a fix.

    So it still raises -- but only after the library has been persisted, so the
    error costs a message rather than the whole project.
    """

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("ic/steady", mf.ic(strt=1.0))
    project.add_simulation(
        mf.SimulationSpec(
            "dangling",
            models=[mf.gwf("flow", packages=[mf.ref("npf/base")])],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )

    with pytest.raises(ValueError, match="npf/base"):
        project.save()

    # ...and the library is on disk anyway.
    assert sorted(mf.Project.load(project.root).packages) == ["ic/steady"]


def test_a_reopened_run_can_be_rebuilt_and_forked_in_a_later_session(tmp_path):
    """The manifest stores the rebuildable recipe, not just a readable summary.

    Before this, ``Run.reopen`` had no ``SimulationSpec`` at all, so a run could
    be read but never re-run or varied without its original build script. The
    summary block records package NAMES; the recipe records what builds them.
    """

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    simulation = mf.SimulationSpec(
        "baseline",
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    mf.ic(strt=1.0),
                    mf.ref("npf/base"),
                    mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )
    project.add_simulation(simulation)
    project.save()
    project.prepare_run("baseline", simulation)

    # A genuinely fresh session: nothing but the directory.
    later = mf.Project.load(tmp_path / "proj")
    reopened = later.reopen_run("baseline", load=False)

    assert reopened.spec is not None
    # The library reference survives as a reference -- not flattened into a copy.
    assert [p.key for p in reopened.spec.models[0].packages if isinstance(p, PackageRef)] == [
        "npf/base"
    ]

    reopened.build()
    assert reopened.built.simulation.get_model("flow").npf.k.array[0][0][0] == 1.0

    # ...and it can be varied, which is what the recipe really buys.
    variant = reopened.spec.derive("high_k").with_packages(mf.npf(k=99.0), on="flow")
    forked = later.prepare_run("high_k", variant)
    assert forked.built.simulation.get_model("flow").npf.k.array[0][0][0] == 99.0
    assert variant.derived_from == "baseline"


def test_a_spec_that_cannot_be_json_encoded_stays_summary_only(tmp_path):
    """An embedded array is not an error here -- the run just stays as
    rebuildable as it was before, with its readable summary intact."""

    project = mf.Project(tmp_path / "proj", name="p")
    simulation = mf.SimulationSpec(
        "inline",
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    mf.npf(k=np.arange(4.0)),
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )
    run = project.prepare_run("inline", simulation, build=False)
    manifest = json.loads(run.manifest_path.read_text())

    assert manifest["spec_recipe"] is None
    assert manifest["spec"]["models"][0]["packages"] == ["dis", "npf"]

    reopened = project.reopen_run("inline", load=False)
    assert reopened.spec is None
    assert reopened.spec_block == manifest["spec"]


def test_an_unrevivable_recipe_degrades_instead_of_failing_to_open(tmp_path):
    """A recipe naming a builder this version cannot import must not make the
    run unopenable -- it degrades to summary-only and warns."""

    project = mf.Project(tmp_path / "proj", name="p")
    run = project.prepare_run("baseline", _tiny_flow_spec())
    manifest = json.loads(run.manifest_path.read_text())
    manifest["spec_recipe"]["packages"][0]["builder"] = "no.such.module:Gone"
    run.manifest_path.write_text(json.dumps(manifest))

    reopened = project.reopen_run("baseline", load=False)

    assert reopened.spec is None
    assert reopened.status == manifest["run"]["status"]
    assert reopened.spec_block is not None
