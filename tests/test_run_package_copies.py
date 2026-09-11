"""A run rebuilds what it was prepared with, not the library as edited since.

Ledger 230. ``prepare_run`` copies each library package the simulation
references into ``<run>/_library/``, and a rebuild resolves ``mf.ref`` from
those copies first. Before this, reopening a run in a later session and
rebuilding it silently used the library as it stood THEN -- the new K, under the
fingerprint of the old one.
"""

from __future__ import annotations

import json
import logging
import threading

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy.specs import PackageSpec


def _simulation(*packages) -> mf.SimulationSpec:
    return mf.SimulationSpec(
        "baseline",
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    *packages,
                    mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )


def _prepared(tmp_path, npf=None) -> mf.Project:
    """A saved project whose run 'baseline' references npf/base (k=1 unless given)."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", npf if npf is not None else mf.npf(k=1.0))
    project.prepare_run("baseline", _simulation(mf.ic(strt=1.0), mf.ref("npf/base")))
    project.save()
    return project


def _edited_later(tmp_path, **options) -> mf.Project:
    """A later session edits npf/base and saves; returns the session after that."""

    session = mf.Project.load(tmp_path / "proj")
    session.packages.edit("npf/base", **options)
    session.save()
    return mf.Project.load(tmp_path / "proj")


def _k(run) -> list[float]:
    return run.built.simulation.get_model("flow").npf.k.array.ravel().tolist()


def _manifest(run_dir) -> dict:
    return json.loads((run_dir / "run.json").read_text())


def _copy_of(run_dir, key):
    return run_dir / _manifest(run_dir)["package_snapshot"][key]


def _as_if_prepared_before_copies(project) -> None:
    path = project.runs_dir / "baseline" / "run.json"
    manifest = json.loads(path.read_text())
    manifest.pop("package_snapshot")
    path.write_text(json.dumps(manifest))


def _warnings(caplog) -> list[str]:
    return [
        record.getMessage()
        for record in caplog.records
        if record.name == "myflopy.workspace" and record.levelno >= logging.WARNING
    ]


def test_a_reopened_run_builds_what_it_was_prepared_with(tmp_path):
    """The defect: this came back 50.0 -- the edited library -- while the run's
    manifest still carried the fingerprint of 1.0."""

    _prepared(tmp_path)
    project = _edited_later(tmp_path, k=50.0)

    run = project.reopen_run("baseline", load=False).build()

    assert _k(run) == [1.0] * 4
    # drift still says the LIBRARY moved on; it no longer means the run did
    assert project.runs.drift().set_index("key").loc["npf/base", "state"] == "edited"


def test_an_array_package_is_copied_as_a_pickle(tmp_path):
    _prepared(tmp_path, npf=mf.npf(k=np.arange(1.0, 5.0).reshape(1, 2, 2)))
    project = _edited_later(tmp_path, k=np.full((1, 2, 2), 9.0))

    run = project.reopen_run("baseline", load=False).build()

    assert _k(run) == [1.0, 2.0, 3.0, 4.0]
    assert json.loads(_copy_of(run.workspace, "npf/base").read_text())["method"] == "pickle"


def test_a_run_opened_without_its_project_rebuilds_from_its_copies(tmp_path):
    """mf.load_run has no library in scope, so a ref could not resolve at all."""

    project = _prepared(tmp_path)

    run = mf.load_run(project.runs_dir / "baseline", load=False).build()

    assert _k(run) == [1.0] * 4


def test_the_copies_survive_the_manifest_being_rewritten(tmp_path):
    """Every status change rewrites run.json; a reopened run must carry the index
    forward, or its second rebuild falls back to the library."""

    _prepared(tmp_path)
    mf.Project.load(tmp_path / "proj").reopen_run("baseline", load=False).build()
    project = _edited_later(tmp_path, k=50.0)

    assert _k(project.reopen_run("baseline", load=False).build()) == [1.0] * 4


def test_a_run_prepared_in_this_session_reads_nothing_back(tmp_path):
    """It builds from the objects it copied; reading them back would unpickle
    every array package a second time for nothing."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    run = project.prepare_run(
        "baseline", _simulation(mf.ic(strt=1.0), mf.ref("npf/base")), build=False
    )
    for path in (run.workspace / "_library").rglob("*.json"):
        path.unlink()

    assert _k(run.build()) == [1.0] * 4


def test_a_run_prepared_before_copies_warns_when_the_library_moved(tmp_path, caplog):
    """It has nothing else to build from, so it still uses the library -- but it
    says so, naming the key, instead of silently building something else."""

    project = _prepared(tmp_path)
    _as_if_prepared_before_copies(project)
    project = _edited_later(tmp_path, k=50.0)
    run = project.reopen_run("baseline", load=False)

    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        run.build()

    assert _k(run) == [50.0] * 4
    assert any("'baseline'" in m and "'npf/base'" in m for m in _warnings(caplog))


def test_no_warning_when_the_library_has_not_moved(tmp_path, caplog):
    project = _prepared(tmp_path)
    _as_if_prepared_before_copies(project)
    run = mf.Project.load(tmp_path / "proj").reopen_run("baseline", load=False)

    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        run.build()

    assert _warnings(caplog) == []


def test_rewriting_such_a_run_records_what_it_now_holds(tmp_path):
    """Its inputs on disk now come from the library's version, so its digests
    say so -- and it keeps a copy, so from here on it is frozen like a new run."""

    project = _prepared(tmp_path)
    _as_if_prepared_before_copies(project)
    project = _edited_later(tmp_path, k=50.0)
    run = project.reopen_run("baseline", load=False).build()

    run.write(silent=True)

    manifest = _manifest(run.workspace)
    assert manifest["package_digests"]["npf/base"] == project.packages.digests()["npf/base"]
    assert "npf/base" in manifest["package_snapshot"]
    drift = project.refresh_runs().drift().set_index("key")
    assert drift.loc["npf/base", "state"] == "unchanged"

    project = _edited_later(tmp_path, k=7.0)
    assert _k(project.reopen_run("baseline", load=False).build()) == [50.0] * 4


def test_a_missing_copy_is_refused_not_replaced_from_the_library(tmp_path):
    project = _prepared(tmp_path)
    _copy_of(project.runs_dir / "baseline", "npf/base").unlink()
    run = mf.Project.load(tmp_path / "proj").reopen_run("baseline", load=False)

    with pytest.raises(FileNotFoundError, match="own copy"):
        run.build()


def test_an_unreadable_copy_is_refused_too(tmp_path):
    project = _prepared(tmp_path, npf=mf.npf(k=np.arange(1.0, 5.0).reshape(1, 2, 2)))
    run_dir = project.runs_dir / "baseline"
    doc = json.loads(_copy_of(run_dir, "npf/base").read_text())
    (run_dir / doc["pickle_path"]).write_bytes(b"not a pickle")
    run = mf.Project.load(tmp_path / "proj").reopen_run("baseline", load=False)

    with pytest.raises(ValueError, match="could not be read back"):
        run.build()


def test_keys_that_share_a_file_name_keep_separate_copies(tmp_path):
    """_slug is lossy: 'base/flow k' and 'base/flow_k' name one file."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("base/flow k", mf.npf(k=1.0))
    project.add_package("base/flow_k", mf.ic(strt=3.0))
    project.prepare_run(
        "baseline", _simulation(mf.ref("base/flow k"), mf.ref("base/flow_k"))
    )

    index = _manifest(project.runs_dir / "baseline")["package_snapshot"]
    assert len(set(index.values())) == 2
    run = mf.load_run(project.runs_dir / "baseline", load=False).build()
    model = run.built.simulation.get_model("flow")
    assert model.npf.k.array.ravel().tolist() == [1.0] * 4
    assert model.ic.strt.array.ravel().tolist() == [3.0] * 4


def test_a_package_that_cannot_be_copied_does_not_stop_the_run(tmp_path, caplog):
    """Neither JSON-encodable nor picklable: left out, with a warning."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package(
        "npf/odd", PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": threading.Lock()})
    )

    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        run = project.prepare_run(
            "baseline", _simulation(mf.ic(strt=1.0), mf.ref("npf/odd")), build=False
        )

    assert "npf/odd" not in _manifest(run.workspace)["package_snapshot"]
    assert any("'npf/odd'" in m and "could not be copied" in m for m in _warnings(caplog))
