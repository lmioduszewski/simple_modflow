"""Reading MF6's own verdict back out of ``mfsim.lst`` (design step 1).

A run manifest's ``status`` is myflopy's verb -- what the library did. Nothing
recorded what MF6 did afterwards, so a run written by one session and executed by
another sat at ``status="written", success=None`` forever. Measured on a real
project: four runs, all reading ``written``/``None``, three of which had reached
``Normal termination``.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

import myflopy as mf
from myflopy.run_record import ListingOutcome, read_listing_outcome

_NORMAL_TAIL = """\
 Run end date and time (yyyy/mm/dd hh:mm:ss): 2026/09/03 12:50:17
 Elapsed run time:  3 Hours, 53 Minutes,  1 Seconds

 Normal termination of simulation.
"""

_PREMATURE_TAIL = """\
 ERROR REPORT:

  1. Failure to converge.
 Premature termination of simulation.
"""

# What a killed run actually looks like: it just stops mid-iteration. Taken from
# the shape of a real interrupted run rather than invented.
_INTERRUPTED_TAIL = """\
 Backtracking                      36                     1          1       241446.42
 Model                             36         40                                148.16560
"""


def _listing(workspace: Path, tail: str) -> Path:
    """Write a workspace holding one mfsim.lst with the given ending."""

    workspace.mkdir(parents=True, exist_ok=True)
    listing = workspace / "mfsim.lst"
    listing.write_text(tail, encoding="utf-8")
    return listing


def test_the_four_outcomes(tmp_path):
    assert read_listing_outcome(_listing(tmp_path / "ok", _NORMAL_TAIL).parent).outcome == "normal"
    assert read_listing_outcome(_listing(tmp_path / "bad", _PREMATURE_TAIL).parent).outcome == "failed"
    assert (
        read_listing_outcome(_listing(tmp_path / "cut", _INTERRUPTED_TAIL).parent).outcome
        == "incomplete"
    )
    assert read_listing_outcome(tmp_path / "nothing_here").outcome == "absent"


def test_an_unfinished_run_has_no_verdict_rather_than_a_bad_one(tmp_path):
    """``succeeded`` is None, not False. Reporting an interrupted run as failure
    is the same class of lie as reporting ``status='written'`` as success."""

    assert read_listing_outcome(_listing(tmp_path / "ok", _NORMAL_TAIL).parent).succeeded is True
    assert read_listing_outcome(_listing(tmp_path / "bad", _PREMATURE_TAIL).parent).succeeded is False
    assert read_listing_outcome(_listing(tmp_path / "cut", _INTERRUPTED_TAIL).parent).succeeded is None
    assert read_listing_outcome(tmp_path / "gone").succeeded is None


def test_mf6s_own_timing_is_reported_verbatim(tmp_path):
    outcome = read_listing_outcome(_listing(tmp_path / "ok", _NORMAL_TAIL).parent)

    assert outcome.finished_at == "2026/09/03 12:50:17"
    assert outcome.elapsed == "3 Hours, 53 Minutes,  1 Seconds"


def test_only_the_tail_is_read(tmp_path):
    """A long transient run writes one block per time step; the verdict is in the
    last few lines, so a huge listing must not be read whole."""

    workspace = tmp_path / "big"
    workspace.mkdir()
    padding = "x" * 5_000_000
    (workspace / "mfsim.lst").write_text(padding + "\n" + _NORMAL_TAIL, encoding="utf-8")

    started = time.perf_counter()
    outcome = read_listing_outcome(workspace)
    elapsed = time.perf_counter() - started

    assert outcome.outcome == "normal"
    assert elapsed < 0.5, f"read took {elapsed:.2f}s -- the whole file was read"


def test_inputs_newer_than_the_listing_mark_it_stale(tmp_path):
    """A listing that predates its inputs describes a superseded run: the verdict
    is real, but no longer about what is on disk now. Observed on a real
    project, where baseline's inputs were rewritten after it ran."""

    workspace = tmp_path / "run"
    listing = _listing(workspace, _NORMAL_TAIL)
    assert read_listing_outcome(workspace).stale is False

    package = workspace / "flow.rch"
    package.write_text("BEGIN period 1\nEND period 1\n", encoding="utf-8")
    os.utime(package, (time.time() + 60, time.time() + 60))

    assert read_listing_outcome(workspace).stale is True


def test_outputs_written_during_the_run_do_not_mark_it_stale(tmp_path):
    """Results are legitimately newer than the listing -- counting them would
    make every completed run report stale."""

    workspace = tmp_path / "run"
    _listing(workspace, _NORMAL_TAIL)
    for name in ("flow.hds", "flow.cbc", "run.json"):
        output = workspace / name
        output.write_text("x", encoding="utf-8")
        os.utime(output, (time.time() + 60, time.time() + 60))

    assert read_listing_outcome(workspace).stale is False


def test_an_unreadable_listing_is_classified_not_raised(tmp_path):
    """A run you cannot classify must still be listable beside the ones you can."""

    workspace = tmp_path / "run"
    workspace.mkdir()
    (workspace / "mfsim.lst").write_bytes(b"\xff\xfe not utf-8 at all")

    assert read_listing_outcome(workspace).outcome == "incomplete"


# --- the Run integration -----------------------------------------------------


def _executed_run(tmp_path: Path, tail: str) -> mf.Run:
    """A prepared run whose workspace holds a listing with the given ending."""

    project = mf.Project(tmp_path / "proj", name="p")
    simulation = mf.SimulationSpec(
        "baseline",
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    mf.ic(strt=1.0),
                    mf.npf(k=1.0),
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )
    run = project.prepare_run("baseline", simulation)
    _listing(run.workspace, tail)
    return run


def test_status_and_outcome_are_two_axes(tmp_path):
    """The measured defect: a run reads ``written``/``None`` while its listing
    says the simulation terminated normally."""

    run = _executed_run(tmp_path, _NORMAL_TAIL)

    assert run.status == "built"          # what myflopy did
    assert run.success is None
    assert run.outcome.outcome == "normal"  # what MF6 did
    assert run.outcome.succeeded is True


def test_reconcile_adopts_the_verdict_and_persists_it(tmp_path):
    run = _executed_run(tmp_path, _NORMAL_TAIL)
    original_spec = json.loads(run.manifest_path.read_text())["spec"]

    outcome = run.reconcile()

    assert outcome.succeeded is True
    assert run.status == "completed"
    assert run.success is True
    manifest = json.loads(run.manifest_path.read_text())
    assert manifest["run"]["status"] == "completed"
    assert manifest["run"]["success"] is True
    # ...and reconciling must not erase what the run contained (ledger 179/181).
    assert manifest["spec"] == original_spec
    assert manifest["spec_recipe"] is not None


def test_reconcile_leaves_an_unfinished_run_alone(tmp_path):
    """No verdict means no change -- not a recorded failure."""

    run = _executed_run(tmp_path, _INTERRUPTED_TAIL)
    before = run.status

    outcome = run.reconcile()

    assert outcome.outcome == "incomplete"
    assert run.status == before
    assert run.success is None


def test_reconcile_records_a_real_failure(tmp_path):
    run = _executed_run(tmp_path, _PREMATURE_TAIL)

    run.reconcile()

    assert run.status == "failed"
    assert run.success is False


def test_the_outcome_is_never_cached_into_the_manifest(tmp_path):
    """The listing is the authority and is right there; a stored verdict would be
    a second source of truth that can disagree with the file it came from."""

    run = _executed_run(tmp_path, _NORMAL_TAIL)
    run.reconcile()

    manifest = json.loads(run.manifest_path.read_text())
    assert "outcome" not in manifest
    assert "elapsed" not in json.dumps(manifest)

    # Rewriting the listing changes the answer immediately, with no invalidation.
    _listing(run.workspace, _PREMATURE_TAIL)
    assert run.outcome.outcome == "failed"


# --- the runs library --------------------------------------------------------


def _run_with(project: mf.Project, name: str, tail: str | None) -> Path:
    """Prepare a run and optionally give it a listing."""

    simulation = mf.SimulationSpec(
        name,
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    mf.ic(strt=1.0),
                    mf.npf(k=1.0),
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )
    run = project.prepare_run(name, simulation)
    if tail is not None:
        _listing(run.workspace, tail)
    return run.workspace


def test_the_runs_library_reports_both_axes(tmp_path):
    project = mf.Project(tmp_path / "proj", name="p")
    _run_with(project, "good", _NORMAL_TAIL)
    _run_with(project, "bad", _PREMATURE_TAIL)
    _run_with(project, "killed", _INTERRUPTED_TAIL)
    _run_with(project, "never_ran", None)

    frame = project.runs.summary().set_index("key")

    assert set(frame.index) == {"good", "bad", "killed", "never_ran"}
    assert frame.loc["good", "outcome"] == "normal"
    assert frame.loc["bad", "outcome"] == "failed"
    assert frame.loc["killed", "outcome"] == "incomplete"
    assert frame.loc["never_ran", "outcome"] == "absent"
    # status is myflopy's axis and is the same for all four
    assert set(frame["status"]) == {"built"}


def test_a_run_without_a_manifest_is_listed_not_hidden(tmp_path):
    """Measured under examples/: 39 MF6 workspaces against 14 run.json files.
    A ledger that omitted two thirds of real runs would be worse than none."""

    project = mf.Project(tmp_path / "proj", name="p")
    workspace = _run_with(project, "owned", _NORMAL_TAIL)
    (workspace.parent / "bare").mkdir()
    (workspace.parent / "bare" / "mfsim.nam").write_text("# hand-run\n", encoding="utf-8")
    _listing(workspace.parent / "bare", _NORMAL_TAIL)

    frame = project.runs.summary().set_index("key")

    assert set(frame.index) == {"owned", "bare"}
    assert frame.loc["bare", "status"] == "discovered"
    # ...keyed the same way whether or not it has a manifest
    assert "runs/" not in " ".join(frame.index)


def test_workspaces_outside_runs_keep_their_path_so_names_cannot_collide(tmp_path):
    project = mf.Project(tmp_path / "proj", name="p")
    _run_with(project, "from_usg", _NORMAL_TAIL)
    stray = project.root / "scratch" / "from_usg"
    stray.mkdir(parents=True)
    (stray / "mfsim.nam").write_text("# elsewhere\n", encoding="utf-8")

    keys = set(project.runs.summary()["key"])

    assert keys == {"from_usg", "scratch/from_usg"}


def test_one_unreadable_manifest_does_not_hide_the_others(tmp_path):
    """`_write_json` is a non-atomic write_text called several times per run, so
    a truncated manifest is a real state. It must be a ROW, not an exception."""

    project = mf.Project(tmp_path / "proj", name="p")
    _run_with(project, "good", _NORMAL_TAIL)
    broken = _run_with(project, "broken", _NORMAL_TAIL)
    (broken / "run.json").write_text('{"run": ', encoding="utf-8")

    frame = project.runs.summary().set_index("key")

    assert frame.loc["broken", "status"] == "unreadable"
    assert frame.loc["good", "outcome"] == "normal"     # still listed, still right


def test_the_library_is_cached_with_an_explicit_refresh(tmp_path):
    """`project.runs["x"]` must return the same object twice -- a property that
    rescanned would make `run.metadata[...] = ...` silently vanish."""

    project = mf.Project(tmp_path / "proj", name="p")
    _run_with(project, "good", _NORMAL_TAIL)

    assert project.runs["good"] is project.runs["good"]

    _run_with(project, "later", _NORMAL_TAIL)
    assert "later" not in project.runs          # cached, deliberately

    assert "later" in project.refresh_runs()


def test_a_missing_run_names_the_near_miss_like_every_other_library(tmp_path):
    project = mf.Project(tmp_path / "proj", name="p")
    _run_with(project, "baseline", _NORMAL_TAIL)

    with pytest.raises(KeyError) as caught:
        project.runs["baselien"]

    assert "Project 'p' has no run 'baselien'." in caught.value.args[0]
    assert "Did you mean: 'baseline'?" in caught.value.args[0]


# --- drift: what an edit invalidated (design step 4 / §6a) ------------------


def _library_project(tmp_path: Path) -> mf.Project:
    """A project whose simulation references library packages, then runs it."""

    import numpy as np

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_package("rch/series", mf.rch.array(recharge={0: np.zeros(4)}, irch=None))
    project.add_package("rch_pond", mf.rch.array(recharge={0: np.ones(4)}, irch=None))
    project.add_simulation(
        mf.SimulationSpec(
            "baseline",
            models=[
                mf.gwf(
                    "flow",
                    packages=[
                        mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                        mf.ic(strt=1.0),
                        mf.ref("npf/base"),
                        mf.ref("rch/series"),
                    ],
                )
            ],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
        )
    )
    project.prepare_run("baseline", project.simulations["baseline"])
    return project


def test_only_library_references_are_fingerprinted(tmp_path):
    """An inline package is baked into the run's own recipe and cannot drift; a
    reference names a library entry that can be edited afterwards."""

    project = _library_project(tmp_path)
    manifest = json.loads((project.runs_dir / "baseline" / "run.json").read_text())

    assert set(manifest["package_digests"]) == {"npf/base", "rch/series"}
    assert "dis" not in manifest["package_digests"]      # inline, cannot drift


def test_editing_a_package_shows_up_as_drift(tmp_path):
    """The defect this exists for: the manifest records package KEYS, so editing
    npf/base leaves every earlier run still naming it -- true of the key, false
    of the model. Measured before this: two models differing 50x in K produced
    byte-identical manifests."""

    project = _library_project(tmp_path)
    assert set(project.runs.drift()["state"]) == {"unchanged"}

    project.packages.edit("npf/base", k=50.0, note="post-calibration")
    project.packages.pop("rch/series")

    drift = project.refresh_runs().drift().set_index("key")

    assert drift.loc["npf/base", "state"] == "edited"
    assert drift.loc["npf/base", "ran_with"] != drift.loc["npf/base", "library_now"]
    assert drift.loc["rch/series", "state"] == "removed"
    assert drift.loc["rch/series", "library_now"] == "—"


def test_the_summary_names_the_drifted_keys(tmp_path):
    project = _library_project(tmp_path)
    assert project.runs.summary().set_index("key").loc["baseline", "drift"] == "—"

    project.packages.edit("npf/base", k=50.0)

    row = project.refresh_runs().summary().set_index("key").loc["baseline"]
    assert row["drift"] == "npf/base"


def test_a_run_recorded_before_digests_reports_unknown(tmp_path):
    """A missing fingerprint is not evidence of sameness. Reporting it as
    ``unchanged`` would be the same class of lie as ``status='written'`` reading
    as success -- which is the defect two commits back."""

    project = _library_project(tmp_path)
    manifest_path = project.runs_dir / "baseline" / "run.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.pop("package_digests")
    manifest_path.write_text(json.dumps(manifest))

    refreshed = project.refresh_runs()

    assert refreshed.drift().iloc[0]["state"] == "unknown"
    assert refreshed.summary().set_index("key").loc["baseline", "drift"] == "unknown"


def test_drift_needs_a_project_bound_library(tmp_path):
    """A library detached by dict.copy() has no packages to compare against, and
    says so rather than reporting everything unchanged."""

    from myflopy.workspace import RunLibrary

    project = _library_project(tmp_path)
    assert project.runs.drift() is not None

    with pytest.raises(RuntimeError, match="not bound to a project"):
        RunLibrary(project.runs).drift()


def test_editing_does_not_disturb_a_runs_own_recipe(tmp_path):
    """Drift reports; it must not rewrite history. The run's stored digests and
    spec blocks are what it RAN with and stay put."""

    project = _library_project(tmp_path)
    manifest_path = project.runs_dir / "baseline" / "run.json"
    before = json.loads(manifest_path.read_text())

    project.packages.edit("npf/base", k=50.0)
    project.refresh_runs().drift()

    assert json.loads(manifest_path.read_text()) == before
