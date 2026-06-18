from __future__ import annotations

import json
import pickle
from types import SimpleNamespace

import pytest

import myflopy as mf
from myflopy.specs import GridRef


def _sim_with_grid(name: str, grid) -> mf.SimulationSpec:
    """Return a one-model simulation whose grid comes from ``grid``."""

    gwf = mf.gwf("gwf", grid=grid, packages=[mf.ic(strt=1.0)])
    return mf.SimulationSpec(
        name,
        models=(gwf,),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),),
    )


def test_grid_ref_round_trip():
    reference = mf.grid_ref("base")

    assert reference.key == "base"
    assert GridRef.from_dict(reference.to_dict()) == reference


def test_grid_from_object_resolves_to_the_object():
    obj = SimpleNamespace(tag="A")
    spec = mf.GridSpec.from_object(obj)

    assert spec.method == "object"
    assert spec.resolve() is obj


def test_grid_from_object_is_not_directly_serializable():
    spec = mf.GridSpec.from_object(SimpleNamespace())

    with pytest.raises(ValueError, match="cannot be serialized"):
        spec.to_dict()


def test_model_grid_ref_serializes_and_round_trips():
    model = mf.gwf("gwf", grid=mf.grid_ref("base"))

    loaded = mf.ModelSpec.from_dict(model.to_dict())

    assert isinstance(loaded.grid, GridRef)
    assert loaded.grid.key == "base"


def test_bare_grid_object_is_wrapped_for_a_model():
    obj = SimpleNamespace(tag="bare")
    model = mf.gwf("gwf", grid=obj)

    assert isinstance(model.grid, mf.GridSpec)
    assert model.grid.method == "object"
    assert model.grid.resolve() is obj


def test_project_grid_library_resolves_and_swaps(tmp_path):
    project = mf.Project(tmp_path / "demo", name="demo")
    project.add_grid("a", mf.GridSpec.from_object(SimpleNamespace(tag="A")))
    project.add_grid("b", mf.GridSpec.from_object(SimpleNamespace(tag="B")))

    baseline = _sim_with_grid("baseline", mf.grid_ref("a"))
    project.add_simulation(baseline)
    swapped = baseline.derive("swapped").replace_grid("gwf", mf.grid_ref("b"))
    project.add_simulation(swapped)

    run_a = project.prepare_run("a_run", "baseline").build()
    run_b = project.prepare_run("b_run", "swapped").build()

    # Same reference machinery, different grid resolved per simulation.
    assert run_a.built.models["gwf"].context.grid.tag == "A"
    assert run_b.built.models["gwf"].context.grid.tag == "B"
    assert swapped.lineage[-1] == {"operation": "replace_grid", "model": "gwf"}


def test_project_unresolved_grid_reference_raises(tmp_path):
    project = mf.Project(tmp_path / "demo", name="demo")
    # Library has no grid named "missing".
    project.add_simulation(_sim_with_grid("baseline", mf.grid_ref("missing")))

    run = project.prepare_run("baseline", "baseline")
    with pytest.raises(KeyError):
        run.build()


def test_grid_from_pickle_resolves(tmp_path):
    obj = SimpleNamespace(tag="pickled")
    path = tmp_path / "g.pkl"
    with path.open("wb") as handle:
        pickle.dump(obj, handle)

    spec = mf.GridSpec.from_pickle(path)

    assert spec.method == "pickle"
    assert spec.resolve().tag == "pickled"


def test_grid_pickle_spec_round_trips():
    spec = mf.GridSpec.from_pickle("specs/grids/base.pkl", name="base")

    loaded = mf.GridSpec.from_dict(spec.to_dict())

    assert loaded.method == "pickle"
    assert loaded.pickle_path == "specs/grids/base.pkl"
    assert loaded.name == "base"


def test_project_saves_and_loads_object_grid(tmp_path):
    project = mf.Project(tmp_path / "demo", name="demo")
    project.add_grid("base", mf.GridSpec.from_object(SimpleNamespace(tag="built")))
    project.add_simulation(_sim_with_grid("baseline", mf.grid_ref("base")))

    project.save()

    assert (project.root / "specs" / "grids" / "base.pkl").exists()
    assert (project.root / "specs" / "grids" / "base.json").exists()
    assert (project.root / "specs" / "grids" / "base.versions.json").exists()

    loaded = mf.Project.load(project.root)
    grid = loaded.grids["base"]

    assert isinstance(grid, mf.GridSpec)
    assert grid.method == "pickle"
    assert grid.resolve(project_root=loaded.root).tag == "built"

    # The reloaded project still builds, resolving the unpickled grid.
    run = loaded.prepare_run("run", "baseline").build()
    assert run.built.models["gwf"].context.grid.tag == "built"


def test_pickled_grid_version_mismatch_warns(tmp_path):
    project = mf.Project(tmp_path / "demo", name="demo")
    project.add_grid("base", mf.GridSpec.from_object(SimpleNamespace(tag="x")))
    project.save()

    sidecar = project.root / "specs" / "grids" / "base.versions.json"
    data = json.loads(sidecar.read_text())
    data["flopy"] = "0.0.0-not-real"
    sidecar.write_text(json.dumps(data))

    with pytest.warns(UserWarning, match="different library"):
        mf.Project.load(project.root)
