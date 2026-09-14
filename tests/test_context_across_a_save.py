"""A model's grid guard must survive being written to a project file and read back.

`ModelSpec.to_dict` cannot carry a grid, a domain or a surfaces table. Inlining them
was considered and refused (ledger 193): it floods `Project.validate`, breaks
`to_yaml`, and measures ~3.3 MB per model at 9,405 cells. But they were dropped in
SILENCE, and every branch of the build-time guard is gated on the context holding
them -- so a model that raised live built clean after a save and a load, and a
package pin survived the round trip and was then never compared.

Sixteen characters survive instead: the grid's own fingerprint, which is enough to
tell one mesh from another.
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import (
    GRID_PIN,
    ModelContext,
    disv_pin,
    dropped_context_fields,
    grid_pin,
)
from myflopy.workspace import Project


def _grid(n=5):
    return rectangular_voronoi(CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1))


def _sim(disv_grid, context_grid, *, name="m"):
    """A simulation whose DISV and whose context need not agree."""

    model = mf.gwf(
        "flow",
        context=ModelContext(
            grid=context_grid,
            domain=np.ones((1, context_grid.ncpl), dtype=int),
        ),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **disv_grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    return mf.simulation(
        model, name=name, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )


def _saved_and_reloaded(tmp_path, simulation):
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(simulation)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        project.save()
    return Project.load(tmp_path / "proj")


# ---------------------------------------------------------------- the fingerprint


def test_a_built_disv_fingerprints_as_its_own_grid():
    """The whole idea: the package can answer for the mesh once the grid is gone."""

    for size in (3, 5, 8):
        grid = _grid(size)
        built = _sim(grid, grid).build_flopy(f"/tmp/_pin_{size}").built_model("flow").model
        assert disv_pin(built.get_package("disv")) == grid_pin(grid)


def test_two_meshes_fingerprint_differently():
    grid, other = _grid(4), _grid(5)
    built = _sim(grid, grid).build_flopy("/tmp/_pin_diff").built_model("flow").model
    assert disv_pin(built.get_package("disv")) != grid_pin(other)


def test_a_disv_that_cannot_be_read_is_unpinned_rather_than_an_error():
    from types import SimpleNamespace

    hostile = SimpleNamespace(
        cell2d=SimpleNamespace(get_data=lambda: (_ for _ in ()).throw(RuntimeError("no"))),
        ncpl=SimpleNamespace(get_data=lambda: 4),
        nvert=SimpleNamespace(get_data=lambda: 9),
    )
    assert disv_pin(hostile) is None
    assert disv_pin(None) is None


# ---------------------------------------------------------------- across the save


def test_the_pin_is_recorded_when_a_model_is_written(tmp_path):
    grid = _grid()
    reloaded = _saved_and_reloaded(tmp_path, _sim(grid, grid))
    context = reloaded.simulations["m"].models[0].context
    assert context.grid is None                      # still not carried
    assert (context.metadata or {})[GRID_PIN] == grid_pin(grid)


def test_a_mismatched_model_is_still_refused_after_a_reload(tmp_path):
    """The regression this exists to close: it raised live and built clean reloaded."""

    reloaded = _saved_and_reloaded(tmp_path, _sim(_grid(4), _grid(5)))
    with pytest.raises(ValueError, match=r"is not the mesh this model was declared against"):
        reloaded.prepare_run("r", "m")


def test_a_matching_model_still_builds_after_a_reload(tmp_path):
    """The control -- the fingerprint must be invisible to a model that agrees."""

    grid = _grid()
    reloaded = _saved_and_reloaded(tmp_path, _sim(grid, grid))
    assert reloaded.prepare_run("r", "m") is not None


def test_the_refusal_names_both_fingerprints_and_the_way_out(tmp_path):
    reloaded = _saved_and_reloaded(tmp_path, _sim(_grid(4), _grid(5)))
    with pytest.raises(ValueError) as caught:
        reloaded.prepare_run("r", "m")
    message = str(caught.value)
    assert grid_pin(_grid(5)) in message             # what the model recorded
    assert "mf.grid_ref(key)" in message             # how to get the rest back


def test_a_model_saved_before_pins_existed_still_builds(tmp_path):
    """No recorded pin means no comparison, not a refusal."""

    import json

    grid = _grid()
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_sim(grid, grid))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        project.save()
    spec_path = project.layout.simulation_spec_path("m")
    payload = json.loads(spec_path.read_text())
    payload["models"][0]["context"]["metadata"].pop(GRID_PIN)
    spec_path.write_text(json.dumps(payload))

    reloaded = Project.load(tmp_path / "proj")
    assert reloaded.prepare_run("r", "m") is not None


# ---------------------------------------------------------------- saying so


def test_saving_names_the_context_fields_it_cannot_write(tmp_path, caplog):
    """`skipped=0` used to read as "nothing was lost" while the guard went quiet."""

    grid = _grid()
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_sim(grid, grid))
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        report = project.save()

    assert ("m.flow", "grid") in report.context_dropped
    assert ("m.flow", "domain") in report.context_dropped
    assert "reload as None" in caplog.text
    assert "mf.grid_ref(key)" in caplog.text


def test_a_project_with_nothing_to_drop_is_quiet(tmp_path, caplog):
    project = Project(tmp_path / "proj", name="p")
    model = mf.gwf(
        "flow",
        packages=[
            mf.dis(nlay=1, nrow=2, ncol=2, delr=10.0, delc=10.0, top=10.0, botm=[0.0]),
            mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds"),
        ],
    )
    project.add_simulation(
        mf.simulation(model, name="plain",
                      packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    )
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        report = project.save()
    assert report.context_dropped == ()
    assert "reload as None" not in caplog.text


def test_dates_are_not_reported_because_nothing_reads_them():
    """`dropped_context_fields` names what is LOST, not merely what is absent.

    ``context.dates`` is not serialized either, and no module in the library reads
    it, so warning about it would be noise. Ledger 251 records the same omission on
    the frozen snapshot for the same reason.
    """

    import datetime

    context = ModelContext(dates=[datetime.datetime(2020, 1, 1)])
    assert dropped_context_fields(context) == ()


def test_the_report_pairs_a_field_with_its_model(tmp_path):
    """Two models, one bad: the report has to say which."""

    grid = _grid()
    good = mf.gwf("plain", packages=[
        mf.dis(nlay=1, nrow=2, ncol=2, delr=10.0, delc=10.0, top=10.0, botm=[0.0]),
        mf.npf(k=1.0), mf.ic(strt=5.0)])
    with_grid = mf.gwf("gridded", context=ModelContext(grid=grid), packages=[
        mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
        mf.npf(k=1.0), mf.ic(strt=5.0)])
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(mf.simulation(
        good, with_grid, name="two",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims()],
    ))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        report = project.save()
    models = {model for model, _field in report.context_dropped}
    assert models == {"two.gridded"}


# ------------------------------------------------- what the review found


def test_a_present_grid_overwrites_a_carried_pin(tmp_path):
    """Otherwise the first pin ever written freezes and outlives its own mesh.

    Load, re-grid the model coherently, save again: the file recorded the ORIGINAL
    fingerprint while the project no longer contained that mesh, and the model was
    then refused forever -- naming a mesh nobody could produce.
    """

    import dataclasses

    first, second = _grid(4), _grid(5)
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_sim(first, first))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        project.save()

    reloaded = Project.load(tmp_path / "proj")
    simulation = reloaded.simulations["m"]
    model = simulation.models[0]
    regridded = dataclasses.replace(
        model,
        context=dataclasses.replace(
            model.context, grid=second, domain=np.ones((1, second.ncpl), dtype=int)
        ),
        packages=tuple(
            mf.disv(nlay=1, top=10.0, botm=[0.0], **second.get_disv_gridprops())
            if package.name == "disv"
            else package
            for package in model.packages
        ),
    )
    reloaded.simulations["m"] = dataclasses.replace(simulation, models=(regridded,))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reloaded.save()

    again = Project.load(tmp_path / "proj")
    assert (again.simulations["m"].models[0].context.metadata or {})[GRID_PIN] == grid_pin(second)
    assert again.prepare_run("r", "m") is not None


def test_a_disv_that_leaves_nvert_to_flopy_pins_as_its_grid(tmp_path):
    """FloPy computes NVERT at write time, so ``nvert.get_data()`` is None.

    Hashing it made a correct model fingerprint differently from its own grid and be
    refused after a reload. The vertex indices already carry the connectivity.
    """

    import flopy

    from myflopy.specs import PackageSpec

    grid = _grid(4)
    props = {k: v for k, v in grid.get_disv_gridprops().items() if k != "nvert"}
    model = mf.gwf("flow", context=ModelContext(grid=grid), packages=[
        PackageSpec("disv", flopy.mf6.ModflowGwfdisv,
                    {"nlay": 1, "top": 10.0, "botm": [0.0], **props}),
        mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds")])
    sim = mf.simulation(model, name="nv",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    built = sim.build_flopy(tmp_path / "ws").built_model("flow").model
    assert disv_pin(built.get_package("disv")) == grid_pin(grid)


def test_one_mesh_written_two_ways_pins_the_same():
    """A ring rotated, or wound the other way, is the same polygon.

    Hashing the raw list made the fingerprint mean "this serialization of this mesh",
    so two identical meshes disagreed and the reloaded guard REFUSED something the
    live guard accepts -- the opposite of a weaker check.
    """

    from myflopy.specs import _canonical_ring

    assert _canonical_ring([0, 1, 2, 3]) == _canonical_ring([1, 2, 3, 0])
    assert _canonical_ring([0, 1, 2, 3]) == _canonical_ring([0, 3, 2, 1])
    assert _canonical_ring([0, 1, 2, 3]) != _canonical_ring([0, 1, 3, 2])
    assert _canonical_ring([]) == ()

    grid = _grid(4)
    props = dict(grid.get_disv_gridprops())
    rows = [list(row) for row in props["cell2d"]]
    rows[0] = rows[0][:4] + rows[0][5:] + [rows[0][4]]        # rotate one ring

    class Rotated:
        crs = grid.crs
        ncpl = grid.ncpl

        def get_disv_gridprops(self):
            return {**props, "cell2d": rows}

    assert grid_pin(Rotated()) == grid_pin(grid)


def test_the_fingerprint_cannot_see_a_coordinate_only_change():
    """A known, stated blind spot -- pinned so it cannot be forgotten.

    ``grid_pin`` hashes CONNECTIVITY, which is what makes it exact and frame-free
    (ledger 250). The cost is that a rescale, a mirror or a unit conversion leaves it
    unchanged. The live guard catches those through cell centres; the reloaded guard
    cannot, and both the refusal text and the save warning say so.
    """

    grid = _grid(4)
    props = grid.get_disv_gridprops()
    metres = [[i, x * 0.3048, y * 0.3048] for i, x, y in props["vertices"]]

    class Rescaled:
        crs = grid.crs
        ncpl = grid.ncpl

        def get_disv_gridprops(self):
            return {**props, "vertices": metres}

    assert grid_pin(Rescaled()) == grid_pin(grid)


def test_the_warning_counts_models_not_field_pairs(tmp_path, caplog):
    """One model carrying grid and domain is one context, not two."""

    grid = _grid()
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_sim(grid, grid))
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        project.save()
    assert "1 model context(s)" in caplog.text
    assert "domain and grid" in caplog.text        # sorted, so the order is stable


def test_the_warning_states_what_the_fingerprint_cannot_catch(tmp_path, caplog):
    grid = _grid()
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_sim(grid, grid))
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        project.save()
    assert "rescale" in caplog.text and "unit conversion" in caplog.text
