"""Load a SimulationSpec from a YAML file, build it, and run MODFLOW 6 (plan §5.6).

    python run.py

`SimulationSpec.from_yaml` is the inverse of `.to_yaml()`: the YAML is just the
JSON-safe `to_dict()` payload written to disk, so no Python spec-building code is
needed to describe the model -- see `model.yaml`. `Project.add_simulation_from_yaml`
does the same and registers the spec in a project's simulation library.
"""

from __future__ import annotations

from pathlib import Path

import myflopy as mf

HERE = Path(__file__).resolve().parent

# 1. Read the model straight from YAML (a Path, a filename str, or YAML text).
sim = mf.SimulationSpec.from_yaml(HERE / "model.yaml")
print("loaded simulation:", sim.name, "with models", [m.name for m in sim.models])

# 2. Build + write + run MODFLOW 6.
built = sim.build_flopy(HERE / "runs" / "yaml_demo")
built.simulation.write_simulation()
success, _ = built.simulation.run_simulation(silent=True)
print("run success:", success)

# 3. Round-trip the other way: to_yaml() reproduces the same dict.
assert mf.SimulationSpec.from_yaml(sim.to_yaml()).to_dict() == sim.to_dict()
print("to_yaml/from_yaml round-trips exactly")

# 4. The Project entry point registers a YAML spec in a reusable library.
project = mf.Project(HERE / "runs" / "project", name="yaml_demo")
project.add_simulation_from_yaml(HERE / "model.yaml")
print("registered in project:", list(project.simulations))
