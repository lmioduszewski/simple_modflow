"""Everything that serializes a grid, tested on a grid flopy actually triangulated.

The suite's other grids come from `rectangular_voronoi`, which hand-builds its
`iverts` from python `range` objects. A real `GridSpec.voronoi` grid does not: its
`cell2d` vertex indices are a MIX of `int` and `numpy.int64`, and its cells have
differing vertex counts.

`_json_value` passed `np.float64` (a `float` subclass) and refused `np.int64` (not an
`int` subclass), so on a real grid `to_dict` raised, `spec_recipe` was None,
`freeze_simulation` died, `Run._freeze` swallowed it as a warning, and the run could
not be rebuilt -- while the run itself executed with Normal termination. 3,315 tests
passed throughout.

This module is the one that would have failed.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import myflopy as mf
from myflopy.project.freeze import thaw_simulation
from myflopy.specs import ModelContext, _json_value, disv_pin, grid_pin
from myflopy.workspace import Project


def _simulation(grid, name="real"):
    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    return mf.simulation(
        model, name=name, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )


def test_the_fixture_really_is_the_awkward_shape(triangulated_grid):
    """If this ever stops holding, the module below stops testing anything."""

    cell2d = triangulated_grid.get_disv_gridprops()["cell2d"]
    kinds = {type(value).__name__ for row in cell2d for value in row[3:]}
    assert "int64" in kinds, "the fixture no longer produces numpy ints"
    widths = {int(row[3]) for row in cell2d}
    assert len(widths) > 1, "the fixture no longer has cells of differing widths"


def test_a_numpy_scalar_serializes():
    """It is a number, and whether one got through was an accident of subclassing."""

    assert _json_value(np.int64(7)) == 7
    assert isinstance(_json_value(np.int64(7)), int)
    assert _json_value(np.float32(1.5)) == pytest.approx(1.5)
    assert _json_value(np.bool_(True)) is True
    assert _json_value(np.array(3)) == 3            # a 0-d array is a scalar too


def test_a_real_grid_model_serializes(triangulated_grid):
    """`to_dict` is the gate every other one of these stands behind."""

    assert _simulation(triangulated_grid).to_dict()["kind"] == "SimulationSpec"


def test_a_real_grid_run_is_frozen_and_rebuildable(tmp_path, triangulated_grid):
    """The claim ledger 251 makes, on the grid a user actually has."""

    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_simulation(triangulated_grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = project.prepare_run("r1", "real")
        run.built.simulation.write_simulation(silent=True)

    assert run.frozen is not None
    assert run.frozen["readable"] is True
    thawed = thaw_simulation(run.workspace, run_name="r1")
    assert thawed is not None
    assert thawed.models[0].context.grid.ncpl == triangulated_grid.ncpl


def test_a_real_grid_run_reopens_and_rebuilds(tmp_path, triangulated_grid):
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_simulation(triangulated_grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = project.prepare_run("r1", "real")
        run.built.simulation.write_simulation(silent=True)

    reopened = Project(tmp_path / "proj", name="p").reopen_run("r1")
    reopened.build()
    assert reopened.spec is not None


def test_a_real_grid_fingerprints_consistently(tmp_path, triangulated_grid):
    """`grid_pin` and `disv_pin` must agree on a mesh with ragged cell widths."""

    built = _simulation(triangulated_grid).build_flopy(tmp_path / "ws")
    disv = built.built_model("m").model.get_package("disv")
    assert disv_pin(disv) == grid_pin(triangulated_grid)


def test_a_real_grid_project_saves_and_reloads(tmp_path, triangulated_grid):
    """The whole chain: save, reload in a fresh Project, rebuild."""

    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_simulation(triangulated_grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = project.save()
    assert report.skipped == ()

    reloaded = Project.load(tmp_path / "proj")
    assert reloaded.prepare_run("r", "real") is not None
