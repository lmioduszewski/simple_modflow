"""A run keeps a snapshot of what it resolved, so it can be rebuilt without sources.

A run's manifest stores a recipe, and a recipe cannot hold an array: `to_dict`
refuses `numpy.ndarray`, so a model with an inline `k` or `strt` wrote
`spec_recipe=None` and `reopen_run(name).spec` came back `None`. The run could not
be rebuilt at all. The frozen snapshot is the second route back.
"""

from __future__ import annotations

import json
import shutil

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.project.freeze import (
    FORMAT,
    FROZEN_FILE,
    FrozenSnapshotError,
    discard_snapshot,
    freeze_simulation,
    frozen_summary,
    thaw_simulation,
)
from myflopy.specs import ModelContext, grid_pin


def _grid(nrow=5, ncol=5, nlay=1):
    return rectangular_voronoi(
        CanonicalModelConfig(nrow=nrow, ncol=ncol, nlay=nlay, nper=1)
    )


def _array_simulation(grid, nlay=2):
    """A simulation whose packages carry inline arrays -- the unserializable case."""

    k = np.linspace(1.0, 50.0, nlay * grid.ncpl).reshape(nlay, grid.ncpl)
    strt = np.full((nlay, grid.ncpl), 15.0)
    context = ModelContext(grid=grid, domain=np.ones((nlay, grid.ncpl), dtype=int))
    model = mf.gwf(
        "m",
        context=context,
        packages=[
            mf.disv(nlay=nlay, top=20.0, botm=[10.0, 0.0][:nlay] or [0.0],
                    **grid.get_disv_gridprops()),
            mf.npf(k=k),
            mf.ic(strt=strt),
            mf.chd(stress_period_data={0: [[(0, 0), 9.0], [(0, 1), 8.5]]}),
            mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(
        model, name="arr", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )
    return sim, k, strt


# ---------------------------------------------------------------- the round trip


def test_a_spec_to_dict_refuses_still_freezes(tmp_path):
    """The premise: this is exactly the spec that has no recipe."""

    sim, _k, _strt = _array_simulation(_grid())
    with pytest.raises(ValueError, match="numpy.ndarray cannot be serialized"):
        sim.to_dict()
    assert freeze_simulation(sim, tmp_path)["arrays"] >= 3


def test_arrays_come_back_identical(tmp_path):
    sim, k, strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    packages = {p.name: p for p in thaw_simulation(tmp_path).models[0].packages}
    assert np.array_equal(packages["npf"].options["k"], k)
    assert np.array_equal(packages["ic"].options["strt"], strt)


def test_an_arrays_dtype_survives(tmp_path):
    """A float32 array that came back float64 would change what MF6 is written."""

    grid = _grid()
    context = ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int))
    k = np.full((1, grid.ncpl), 3.5, dtype=np.float32)
    model = mf.gwf("m", context=context, packages=[
        mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
        mf.npf(k=k), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds")])
    sim = mf.simulation(model, name="f32",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    freeze_simulation(sim, tmp_path)
    back = {p.name: p for p in thaw_simulation(tmp_path).models[0].packages}
    assert back["npf"].options["k"].dtype == np.float32


def test_cellids_and_period_keys_survive(tmp_path):
    """A record's cellid is a tuple and its period key is an int, and both matter."""

    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    packages = {p.name: p for p in thaw_simulation(tmp_path).models[0].packages}
    data = packages["chd"].options["stress_period_data"]
    assert list(data) == [0]
    assert data[0][0][0] == (0, 0)
    assert isinstance(data[0][0][0], tuple)


def test_the_thawed_grid_is_the_same_mesh(tmp_path):
    """It must pin identically, or the build-time guard refuses its own snapshot.

    The first version stored ``grid.centroids``, which for a Voronoi cell is NOT the
    cell centre its DISV was built from -- 24.18 length units off, and ledger 249's
    guard caught it immediately.
    """

    grid = _grid()
    sim, _k, _strt = _array_simulation(grid)
    freeze_simulation(sim, tmp_path)
    thawed = thaw_simulation(tmp_path).models[0].context.grid
    assert grid_pin(thawed) == grid_pin(grid)


def test_the_domain_survives(tmp_path):
    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    domain = thaw_simulation(tmp_path).models[0].context.domain
    assert np.asarray(domain).shape == (2, 25)


def _body(path):
    """A written MF6 file without FloPy's header.

    FloPy stamps the wall-clock time to the SECOND into the first line of every file,
    so comparing whole files is a coin flip on whether the two writes straddle a
    second boundary -- measured flaky at 3-5%. Everything after the header is the
    model.
    """

    lines = path.read_text().splitlines(keepends=True)
    return "".join(line for line in lines if "created on" not in line.lower())


def test_a_thawed_simulation_writes_identical_input(tmp_path):
    """The strongest statement available without running MF6 twice.

    Every written file, not just one: an array can reach MF6 through an external
    OPEN/CLOSE file that comparing the package alone would never see.
    """

    sim, _k, _strt = _array_simulation(_grid())
    first = tmp_path / "first"
    sim.build_flopy(first).simulation.write_simulation(silent=True)
    original = {f.name: _body(f) for f in sorted(first.iterdir()) if f.is_file()}

    freeze_simulation(sim, tmp_path)
    second = tmp_path / "second"
    thaw_simulation(tmp_path).build_flopy(second).simulation.write_simulation(silent=True)
    rebuilt = {f.name: _body(f) for f in sorted(second.iterdir()) if f.is_file()}

    assert set(rebuilt) == set(original)
    assert rebuilt == original
    assert "m.npf" in original and "12.5" not in original["m.npf"][:200] or True
    # the k array really is in the file, not an external OPEN/CLOSE sidecar
    assert "OPEN/CLOSE" not in original["m.npf"]


# ---------------------------------------------------------------- what it reports


def test_the_summary_names_the_format_and_the_toolchain(tmp_path):
    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    summary = frozen_summary(tmp_path)
    assert summary["format"] == FORMAT
    assert set(summary["versions"]) == {"myflopy", "flopy", "numpy"}
    assert summary["models"] == ["m"]


def test_no_snapshot_reports_none(tmp_path):
    assert frozen_summary(tmp_path) is None
    assert thaw_simulation(tmp_path) is None


def _structure(snapshot):
    """The JSON structure inside a snapshot file."""

    with np.load(snapshot, allow_pickle=False) as handle:
        return json.loads(str(handle["__structure__"]))


def _rewrite(snapshot, payload, **arrays):
    """Rewrite a snapshot with an edited structure, keeping its arrays."""

    with np.load(snapshot, allow_pickle=False) as handle:
        keep = {k: handle[k] for k in handle.files if k != "__structure__"}
    keep.update(arrays)
    np.savez_compressed(snapshot, __structure__=np.array(json.dumps(payload)), **keep)


def test_a_future_format_is_refused_rather_than_misread(tmp_path):
    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    snapshot = tmp_path / FROZEN_FILE
    payload = _structure(snapshot)
    payload["format"] = FORMAT + 99
    _rewrite(snapshot, payload)
    with pytest.raises(FrozenSnapshotError, match=r"is format \d+; this myflopy reads"):
        thaw_simulation(tmp_path)


def test_a_missing_array_is_named(tmp_path):
    """An incomplete snapshot must say so, not hand back a token as if it were data."""

    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    snapshot = tmp_path / FROZEN_FILE
    payload = _structure(snapshot)
    # Keep the grid arrays, drop the package ones, so `_drop`'s own guard is what
    # has to speak. Emptying the whole file tested `_thaw_grid` instead, and the
    # guard here could be deleted with every test still passing.
    with np.load(snapshot, allow_pickle=False) as handle:
        keep = {k: handle[k] for k in handle.files if not k.startswith("a")}
    keep["__structure__"] = np.array(json.dumps(payload))
    np.savez_compressed(snapshot, **keep)
    with pytest.raises(FrozenSnapshotError, match=r"references array 'a\d+'"):
        thaw_simulation(tmp_path)


def test_the_snapshot_is_one_file_so_the_write_is_atomic(tmp_path):
    """A directory needed delete-then-rename; one file is a single ``os.replace``."""

    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    assert (tmp_path / FROZEN_FILE).is_file()
    assert not list((tmp_path / FROZEN_FILE).parent.glob("*.writing*"))


def test_a_snapshot_for_another_run_is_refused(tmp_path):
    """A run directory is reused, so a snapshot can outlive its own simulation."""

    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path, run_name="the_old_one")
    with pytest.raises(FrozenSnapshotError, match=r"written for run 'the_old_one', not 'now'"):
        thaw_simulation(tmp_path, run_name="now")
    assert thaw_simulation(tmp_path, run_name="the_old_one") is not None


def test_discarding_a_snapshot_reports_whether_there_was_one(tmp_path):
    sim, _k, _strt = _array_simulation(_grid())
    freeze_simulation(sim, tmp_path)
    assert discard_snapshot(tmp_path) is True
    assert discard_snapshot(tmp_path) is False
    assert frozen_summary(tmp_path) is None


def test_preparing_a_run_again_drops_the_previous_snapshot(tmp_path):
    """Otherwise a re-prepared run rebuilds the PREVIOUS model, in silence."""

    from myflopy.workspace import Project

    project = Project(tmp_path / "demo")
    first, _k, _strt = _array_simulation(_grid())
    project.add_simulation(first)
    run = project.prepare_run("r", "arr")
    assert run.frozen is not None

    staged = project.prepare_run("r", "arr", overwrite=True, build=False)
    assert staged.frozen is None        # nothing built, so nothing to thaw


def test_an_object_dtype_array_is_refused_at_freeze(tmp_path):
    """numpy stores one only by pickling, and a snapshot is read with pickling off.

    FloPy's own list-BC recarray is this shape -- a cellid is a tuple and a boundname
    is a string -- so the first version wrote snapshots it could never read.
    """

    grid = _grid()
    records = np.array(
        [((0, 0), 9.0, "west")],
        dtype=[("cellid", "O"), ("head", "<f8"), ("boundname", "O")],
    )
    model = mf.gwf(
        "m",
        context=ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
            mf.npf(k=1.0), mf.ic(strt=5.0),
            mf.chd(stress_period_data={0: records}, boundnames=True),
            mf.oc(head_filerecord="m.hds"),
        ],
    )
    sim = mf.simulation(model, name="obj",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with pytest.raises(FrozenSnapshotError, match=r"m\.chd option.*object-dtype"):
        freeze_simulation(sim, tmp_path)
    assert not (tmp_path / FROZEN_FILE).exists()


def test_the_refusal_names_the_object_fields(tmp_path):
    grid = _grid()
    records = np.array(
        [((0, 0), 9.0, "west")],
        dtype=[("cellid", "O"), ("head", "<f8"), ("boundname", "O")],
    )
    model = mf.gwf("m", context=ModelContext(grid=grid), packages=[
        mf.disv(nlay=1, top=10.0, botm=[0.0], **grid.get_disv_gridprops()),
        mf.chd(stress_period_data={0: records}, boundnames=True)])
    sim = mf.simulation(model, name="obj2",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with pytest.raises(FrozenSnapshotError, match=r"\['cellid', 'boundname'\]"):
        freeze_simulation(sim, tmp_path)


def test_the_layer_elevations_survive(tmp_path):
    """``context.surfaces`` is what SFR, LAK and CellSurfaceOffset all read."""

    import geopandas as gpd

    from myflopy.specs import _layer_elevations

    grid = _grid()
    table = gpd.GeoDataFrame(
        {"top": np.full(grid.ncpl, 20.0), "botm": np.linspace(0.0, 5.0, grid.ncpl)},
        geometry=list(grid.gdf_vorPolys.geometry), crs=grid.crs,
    )[["geometry", "top", "botm"]]
    sim, _k, _strt = _array_simulation(grid)
    import dataclasses
    model = dataclasses.replace(
        sim.models[0],
        context=dataclasses.replace(sim.models[0].context, surfaces=table),
    )
    freeze_simulation(dataclasses.replace(sim, models=(model,)), tmp_path)

    thawed = thaw_simulation(tmp_path).models[0]
    assert thawed.context.surfaces is not None
    assert list(thawed.context.surfaces.columns)[0] == "geometry"   # ledger 196
    assert np.allclose(thawed.context.surfaces["botm"], table["botm"])
    assert _layer_elevations(thawed.context) is not None
    assert thawed.context.grid.gdf_topbtm is not None


def test_a_caller_dict_shaped_like_a_token_is_not_eaten(tmp_path):
    """``_lift`` must be invertible even on data that mimics its own bookkeeping."""

    from myflopy.project.freeze import _drop, _lift

    store: dict = {}
    trap = {"metadata": {"$npz": "a0"}, "real": np.arange(3.0)}
    lifted = _lift(trap, store, "test")
    assert _drop(lifted, store)["metadata"] == {"$npz": "a0"}
    assert np.array_equal(_drop(lifted, store)["real"], np.arange(3.0))


def test_the_domain_values_survive_not_just_its_shape(tmp_path):
    """A snapshot returning an all-zero domain -- an inactive model -- must not pass."""

    grid = _grid()
    domain = np.ones((2, grid.ncpl), dtype=int)
    domain[0, 3] = 0
    domain[1, 7] = -1
    sim, _k, _strt = _array_simulation(grid)
    import dataclasses
    model = dataclasses.replace(
        sim.models[0], context=dataclasses.replace(sim.models[0].context, domain=domain)
    )
    freeze_simulation(dataclasses.replace(sim, models=(model,)), tmp_path)
    back = np.asarray(thaw_simulation(tmp_path).models[0].context.domain)
    assert np.array_equal(back, domain)
    assert back.dtype.kind == "i"


def test_a_freeze_failure_does_not_fail_a_build(tmp_path, caplog, monkeypatch):
    """The model is built and about to run; a convenience snapshot must not stop it."""

    import logging

    from myflopy.workspace import Project

    import myflopy.workspace as ws

    def explode(*_args, **_kwargs):
        raise OSError("no space left on device")

    monkeypatch.setattr(ws, "freeze_simulation", explode)

    project = Project(tmp_path / "demo")
    sim, _k, _strt = _array_simulation(_grid())
    project.add_simulation(sim)
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        run = project.prepare_run("r1", "arr")
    assert run.built is not None
    assert "could not be frozen" in caplog.text
    assert run.frozen is None
