"""YAML (de)serialization for SimulationSpec (plan §5.6).

``specs_io`` is a thin file-format wrapper over the already-complete
``to_dict``/``from_dict`` round-trip (proved package-by-package in
``test_project_spec.py``), so these tests focus on the YAML/file/Project seam:
that YAML preserves the dict faithfully across every spec kind, that a
hand-written minimal YAML builds a real model, and that the file + Project
entry points work.
"""

from __future__ import annotations

import pytest

import myflopy as mf
from myflopy.builders import build_gwf_gwt_exchange
from myflopy.specs import ExchangeSpec, PostBuildHook


def _record_hook(model, packages, context):
    """Module-level hook callback (importable, so it serializes by reference)."""

    return {"ran": True}


def _disv_kwargs() -> dict:
    return dict(
        nlay=1,
        ncpl=2,
        nvert=6,
        vertices=[[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 1.0, 1.0], [3, 0.0, 1.0], [4, 2.0, 0.0], [5, 2.0, 1.0]],
        cell2d=[[0, 0.5, 0.5, 4, 0, 1, 2, 3], [1, 1.5, 0.5, 4, 1, 4, 5, 2]],
        top=10.0,
        botm=0.0,
    )


def _rich_simulation() -> mf.SimulationSpec:
    """A simulation exercising list BCs, a package ref, exchanges, and a hook."""

    flow = mf.gwf(
        "flow",
        packages=[
            mf.disv(**_disv_kwargs()),
            mf.ic(strt=5.0),
            "npf/base",  # a PackageRef, resolved from a project library
            mf.chd(stress_period_data={0: [[(0, 0), 10.0]]}),
            mf.wel(stress_period_data={0: [[(0, 1), -5.0]]}),
            mf.rch(stress_period_data={0: [[(0, 0), 0.001]]}),
            mf.oc(budget_filerecord="flow.cbc", saverecord=[("BUDGET", "ALL")]),
        ],
    ).with_hook(PostBuildHook("record", _record_hook))
    transport = mf.gwt("transport", packages=[mf.disv(**_disv_kwargs()), mf.ic(strt=0.0)])
    return mf.SimulationSpec(
        "coupled",
        models=(flow, transport),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow", "transport"])),
        exchanges=(ExchangeSpec("flow-transport", build_gwf_gwt_exchange, models=("flow", "transport")),),
    )


def test_yaml_string_round_trip_covers_bcs_refs_exchanges_hooks():
    """to_yaml -> from_yaml is identical to the dict for every spec kind (§5.6)."""

    sim = _rich_simulation()
    text = sim.to_yaml()
    assert text.startswith("kind: SimulationSpec\n")  # readable, field order preserved

    loaded = mf.SimulationSpec.from_yaml(text)
    assert loaded.to_dict() == sim.to_dict()
    # the pieces the plan calls out explicitly survived the YAML trip
    assert loaded.exchanges[0].builder is build_gwf_gwt_exchange
    assert loaded.model("flow").hooks[0].callback is _record_hook
    assert loaded.model("flow").package("npf").key == "npf/base"  # PackageRef


def test_yaml_file_round_trip_accepts_path_and_filename(tmp_path):
    """to_yaml(path) writes it; from_yaml reads a Path or a filename string (§5.6)."""

    sim = _rich_simulation()
    dest = tmp_path / "coupled.yaml"
    returned = sim.to_yaml(dest)
    assert dest.read_text() == returned  # to_yaml returns the same text it wrote

    from_path = mf.SimulationSpec.from_yaml(dest)
    from_name = mf.SimulationSpec.from_yaml(str(dest))
    assert from_path.to_dict() == from_name.to_dict() == sim.to_dict()


def test_hand_written_minimal_yaml_builds(tmp_path):
    """A YAML file typed by hand (no Python spec objects) builds a real MF6 model."""

    yaml_text = """
kind: SimulationSpec
name: minimal
models:
- kind: ModelSpec
  name: gwf
  model_type: gwf
  packages:
  - kind: PackageSpec
    name: dis
    builder: myflopy.builders:build_dis
    options: {nlay: 1, nrow: 2, ncol: 2, delr: 100.0, delc: 100.0, top: 10.0, botm: 0.0}
  - kind: PackageSpec
    name: ic
    builder: myflopy.builders:build_ic
    options: {strt: 5.0}
  - kind: PackageSpec
    name: npf
    builder: flopy.mf6.modflow.mfgwfnpf:ModflowGwfnpf
    options: {k: 1.0}
  - kind: PackageSpec
    name: chd
    builder:
      partial: myflopy.advanced:_build_named
      args:
      - {$callable: flopy.mf6.modflow.mfgwfchd:ModflowGwfchd}
      keywords: {}
    options:
      pname: chd
      filename: "{model_name}.chd"
      stress_period_data: {0: [[[0, 0, 0], 10.0], [[0, 1, 1], 9.0]]}
  - kind: PackageSpec
    name: oc
    builder: myflopy.builders:build_oc
    options: {budget_filerecord: gwf.cbc, saverecord: [[BUDGET, ALL]]}
packages:
- kind: PackageSpec
  name: tdis
  builder: flopy.mf6.modflow.mftdis:ModflowTdis
  options: {nper: 1, perioddata: [[1.0, 1, 1.0]]}
- kind: PackageSpec
  name: ims
  builder: myflopy.builders:build_ims
  options: {models: [gwf], pname: ims}
"""
    sim = mf.SimulationSpec.from_yaml(yaml_text)
    assert sim.name == "minimal"
    built = sim.build_flopy(tmp_path)
    flopy_model = built.simulation.get_model("gwf")
    assert flopy_model.get_package("chd") is not None
    assert flopy_model.get_package("dis") is not None


def test_project_add_simulation_from_yaml(tmp_path):
    """Project.add_simulation_from_yaml reads a file and registers the spec (§5.6)."""

    sim = _rich_simulation()
    dest = tmp_path / "coupled.yaml"
    sim.to_yaml(dest)

    project = mf.Project(tmp_path / "proj", name="p")
    registered = project.add_simulation_from_yaml(dest)
    assert registered.name == "coupled"
    assert project.simulation("coupled").to_dict() == sim.to_dict()


def test_from_yaml_rejects_non_mapping_and_bad_type():
    """A YAML scalar/list is not a SimulationSpec; a non-str/Path source is a TypeError."""

    with pytest.raises(ValueError, match="expected a mapping"):
        mf.SimulationSpec.from_yaml("- just\n- a\n- list\n")
    with pytest.raises(TypeError, match="str or Path"):
        mf.SimulationSpec.from_yaml(123)  # type: ignore[arg-type]
