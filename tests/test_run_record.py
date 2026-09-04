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
