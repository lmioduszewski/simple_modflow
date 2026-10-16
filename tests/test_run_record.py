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

import numpy as np
import pytest

import myflopy as mf
from myflopy.run_record import read_listing_outcome
from myflopy.specs import ModelContext
from myflopy.workspace import Project, Run

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


def _nam_workspace(workspace: Path) -> None:
    """A workspace whose name files declare their inputs, as MF6 writes them."""

    workspace.mkdir(parents=True, exist_ok=True)
    (workspace / "mfsim.nam").write_text(
        "BEGIN timing\n  TDIS6  s.tdis\nEND timing\n"
        "BEGIN models\n  gwf6  flow.nam  flow\nEND models\n",
        encoding="utf-8",
    )
    (workspace / "flow.nam").write_text(
        "BEGIN packages\n  DIS6 flow.dis\n  NPF6 flow.npf\nEND packages\n",
        encoding="utf-8",
    )
    for name in ("s.tdis", "flow.dis", "flow.npf"):
        (workspace / name).write_text("x", encoding="utf-8")


def test_stale_asks_the_name_files_not_the_file_extension(tmp_path):
    """An earlier version deny-listed seven output suffixes and treated
    everything else as an input, so a figure or a note saved into a run
    directory marked the run stale -- and it MISSED MF6's own ``.cbb`` budget
    file, an output the list did not name. Measured: adding ``heads.png``,
    ``notes.txt`` or ``flow.cbb`` to a clean finished run all flipped it.

    The name files declare the inputs, so they are what gets asked.
    """

    from myflopy.run_record import _scan_simulation_inputs

    workspace = tmp_path / "run"
    _nam_workspace(workspace)
    past = time.time() - 600
    for item in workspace.iterdir():
        os.utime(item, (past, past))
    _listing(workspace, _NORMAL_TAIL)          # newest

    assert {f.name for f in _scan_simulation_inputs(workspace)[0]} == {
        "mfsim.nam", "flow.nam", "s.tdis", "flow.dis", "flow.npf",
    }
    assert read_listing_outcome(workspace).stale is False

    for name in ("heads.png", "notes.txt", "flow.cbb", "flow.hds"):
        stray = workspace / name
        stray.write_text("x", encoding="utf-8")
        os.utime(stray, (time.time() + 60, time.time() + 60))
        assert read_listing_outcome(workspace).stale is False, name
        stray.unlink()

    # ...but editing a declared INPUT still marks it, which is the whole point.
    os.utime(workspace / "flow.npf", (time.time() + 60, time.time() + 60))
    assert read_listing_outcome(workspace).stale is True


def test_a_workspace_without_a_name_file_cannot_tell(tmp_path):
    """No mfsim.nam means no declared inputs, which reads as "cannot tell" --
    not as "nothing changed"."""

    from myflopy.run_record import _scan_simulation_inputs

    workspace = tmp_path / "run"
    _listing(workspace, _NORMAL_TAIL)

    assert _scan_simulation_inputs(workspace)[0] == []
    assert read_listing_outcome(workspace).stale is False


def test_inputs_newer_than_the_listing_mark_it_stale(tmp_path):
    """A listing that predates its inputs describes a superseded run: the verdict
    is real, but no longer about what is on disk now. Observed on a real
    project, where baseline's inputs were rewritten after it ran."""

    workspace = tmp_path / "run"
    _nam_workspace(workspace)
    _listing(workspace, _NORMAL_TAIL)
    assert read_listing_outcome(workspace).stale is False

    package = workspace / "flow.npf"          # declared by flow.nam
    package.write_text("BEGIN options\nEND options\n", encoding="utf-8")
    os.utime(package, (time.time() + 60, time.time() + 60))

    assert read_listing_outcome(workspace).stale is True


def test_stale_reads_the_inputs_the_pest_build_copies(tmp_path):
    """One rule for "the files MF6 reads" (ledger 611): ``stale`` read the name
    files' entries alone, so an external array (OPEN/CLOSE) or a FILEIN file edited after
    the run left the run reading as current."""

    from myflopy import run_record
    from myflopy.modflow.mf6.pest import project as pest_project

    assert run_record._scan_simulation_inputs is pest_project._scan_simulation_inputs
    workspace = tmp_path / "run"
    _nam_workspace(workspace)
    (workspace / "flow.npf").write_text(
        "BEGIN griddata\n  k\n    OPEN/CLOSE  'k.txt'\nEND griddata\n", encoding="utf-8"
    )
    (workspace / "k.txt").write_text("1.0\n", encoding="utf-8")
    past = time.time() - 600
    for item in workspace.iterdir():
        os.utime(item, (past, past))
    _listing(workspace, _NORMAL_TAIL)
    assert read_listing_outcome(workspace).stale is False

    os.utime(workspace / "k.txt", (time.time() + 60, time.time() + 60))
    assert read_listing_outcome(workspace).stale is True


def test_outputs_written_during_the_run_do_not_mark_it_stale(tmp_path):
    """Results are legitimately newer than the listing -- counting them would
    make every completed run report stale."""

    workspace = tmp_path / "run"
    _nam_workspace(workspace)
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


@pytest.mark.skipif(__import__("shutil").which("mf6") is None, reason="needs mf6")
def test_the_summary_has_an_inputs_column(tmp_path):
    """Whether each run's files are what its record writes (``run.inputs``, ledger
    611), beside MF6's ``outcome`` and ``stale``."""

    project = mf.Project(tmp_path / "proj", name="p")
    _run_with(project, "never_ran", None)
    for name in ("clean", "edited", "old"):
        workspace = _run_with(project, name, None)
        Run.reopen(workspace, load=False).write(silent=True)
    ic = next(project.runs_dir.glob("edited/*.ic"))
    ic.write_text(ic.read_text().replace("CONSTANT", "constant  "), encoding="utf-8")
    manifest = project.runs_dir / "old" / "run.json"
    data = json.loads(manifest.read_text())
    data["metadata"].pop("inputs")
    manifest.write_text(json.dumps(data))

    frame = project.refresh_runs().summary().set_index("key")

    assert list(frame.columns).index("inputs") == list(frame.columns).index("stale") + 1
    assert frame.loc["clean", "inputs"] == "ok"
    assert frame.loc["edited", "inputs"] == f"edited ({ic.name})"
    assert frame.loc["old", "inputs"] == "unfingerprinted"
    assert frame.loc["never_ran", "inputs"] == "—"


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


# --- defects found by adversarial verification of drift ----------------------


def test_the_digest_is_insensitive_to_representation_not_content(tmp_path):
    """``Project.save()``/``load()`` round-trips options through JSON and hands
    them back in a different order, so an order-sensitive digest reported EVERY
    reloaded package as edited -- which would have made drift useless on exactly
    the workflow it exists for."""

    import flopy

    from myflopy.specs import PackageSpec, package_digest

    original = mf.npf(k=1.0, icelltype=1)
    reordered = PackageSpec(
        "npf",
        original.builder,
        {key: original.options[key] for key in reversed(list(original.options))},
    )
    assert package_digest(original) == package_digest(reordered)

    # requires order is not observable by MF6 either
    left = PackageSpec("x", flopy.mf6.ModflowGwfnpf, {}, requires=("a", "b"))
    right = PackageSpec("x", flopy.mf6.ModflowGwfnpf, {}, requires=("b", "a"))
    assert package_digest(left) == package_digest(right)

    # ...but a real change still moves it
    assert package_digest(original) != package_digest(mf.npf(k=50.0, icelltype=1))


def test_the_digest_survives_a_real_save_and_reload(tmp_path):
    project = _library_project(tmp_path)
    project.save()

    reloaded = mf.Project.load(project.root)
    reloaded._runs = None

    assert set(reloaded.runs.drift()["state"]) == {"unchanged"}


def test_a_package_that_cannot_be_fingerprinted_does_not_break_prepare_run(tmp_path):
    """A bookkeeping column must never gate the primary build verb. A closure
    builder was buildable and runnable before drift existed, and still is."""

    import flopy

    from myflopy.specs import PackageSpec

    def make_builder():
        def build_npf(model, **options):
            return flopy.mf6.ModflowGwfnpf(model, **options)

        return build_npf

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/custom", PackageSpec("npf", make_builder(), {"k": 1.0}))
    project.add_simulation(
        mf.SimulationSpec(
            "s",
            models=[
                mf.gwf(
                    "flow",
                    packages=[
                        mf.dis(nlay=1, nrow=1, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                        mf.ref("npf/custom"),
                    ],
                )
            ],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )
    assert project.validate() == []

    run = project.prepare_run("r", "s", build=False)      # used to raise

    assert json.loads(run.manifest_path.read_text())["package_digests"] == {
        "npf/custom": None
    }
    # ...and drift says it cannot tell, rather than guessing at "unchanged"
    assert project.runs.drift().iloc[0]["state"] == "unfingerprinted"


def test_a_detached_library_does_not_report_a_clean_bill_of_health(tmp_path):
    """``drift()`` refuses when unbound; ``summary()`` used to print "—" for the
    same library, which reads as "nothing changed"."""

    from myflopy.workspace import RunLibrary

    project = _library_project(tmp_path)
    detached = RunLibrary(project.runs)
    detached._packages = None

    assert detached.summary().set_index("key").loc["baseline", "drift"] == "unbound"

    # ...and a copy() keeps its binding rather than becoming one of these
    assert project.runs.copy()._packages is project.packages


def test_the_drift_column_fingerprints_the_library_once_per_call(tmp_path):
    """Digesting per (run, package) made summary() O(runs x packages) and
    re-hashed the same arrays for every run -- measured ~1 ms to ~3 s."""

    import time

    import numpy as np

    project = mf.Project(tmp_path / "proj", name="p")
    for index in range(6):
        project.add_package(f"npf/v{index}", mf.npf(k=np.random.rand(9405)))
    refs = [mf.ref(f"npf/v{index}") for index in range(6)]
    project.add_simulation(
        mf.SimulationSpec(
            "s",
            models=[
                mf.gwf(
                    "flow",
                    packages=[
                        mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                        *refs,
                    ],
                )
            ],
            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
        )
    )
    for index in range(20):
        project.prepare_run(f"r{index}", "s", build=False)

    started = time.perf_counter()
    frame = project.runs.summary()
    elapsed = time.perf_counter() - started

    assert len(frame) == 20
    assert elapsed < 1.0, f"summary() took {elapsed:.2f}s -- the memo is not working"


def test_metadata_is_excluded_from_the_fingerprint(tmp_path):
    """Provenance says where a spec came from, not what it builds -- stamping a
    note must not read as a change to the model."""

    from myflopy.specs import package_digest

    project = _library_project(tmp_path)
    before = package_digest(project.package("npf/base"))

    project.packages.edit("npf/base", note="just a note")

    assert package_digest(project.package("npf/base")) == before
    assert project.package("npf/base").metadata["edited"]["note"] == "just a note"
    assert project.runs.summary().set_index("key").loc["baseline", "drift"] == "—"


def test_a_re_prepared_run_does_not_inherit_the_previous_verdict(tmp_path):
    """``prepare_run(overwrite=True)`` reuses the directory and leaves the old
    listing and results in place (design D5, still open). Step 1's outcome reader
    then attributed the PREVIOUS execution's verdict to a run that had not been
    executed -- measured: ``outcome='normal', elapsed='3 Hours, 53 Minutes'`` for
    a model that had just been re-prepared. That is the class of lie this whole
    line of work exists to remove, reintroduced one axis over.
    """

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
    _listing(run.workspace, _NORMAL_TAIL)          # written AFTER preparation
    run.reconcile()
    assert run.status == "completed"

    # Backdate the listing so the re-preparation is unambiguously later, then
    # re-prepare -- exactly what overwrite=True does to a directory that already
    # holds a finished run.
    stale_time = time.time() - 600
    os.utime(run.workspace / "mfsim.lst", (stale_time, stale_time))
    again = project.prepare_run("baseline", simulation, overwrite=True)

    assert again.outcome.outcome == "superseded"
    assert again.outcome.succeeded is None
    # ...so reconcile() refuses to re-adopt a verdict that is not this run's
    again.reconcile()
    assert again.success is None
    assert project.refresh_runs().summary().iloc[0]["outcome"] == "superseded"


def test_a_discovered_run_reports_its_files_timestamps_not_the_scans(tmp_path):
    """``Run.created_at`` defaults to construction time, so a workspace scanned
    without a manifest reported the SCAN's wall clock in the ``updated`` column.
    Measured on the real project: two discovered runs 55 microseconds apart --
    that is the loop, not the runs."""

    project = mf.Project(tmp_path / "proj", name="p")
    bare = project.runs_dir / "hand_run"
    bare.mkdir(parents=True)
    (bare / "mfsim.nam").write_text("# hand-run\n", encoding="utf-8")
    _listing(bare, _NORMAL_TAIL)
    old = time.time() - 86_400
    for item in bare.iterdir():
        os.utime(item, (old, old))

    row = project.runs.summary().set_index("key").loc["hand_run"]

    assert row["status"] == "discovered"
    assert row["updated"].startswith(
        __import__("datetime").datetime.fromtimestamp(
            old, __import__("datetime").timezone.utc
        ).strftime("%Y-%m-%d")
    )
    # ...and its listing is NOT called superseded: a scanned run has no known
    # preparation time to be older than.
    assert row["outcome"] == "normal"


# --- PEST calibrations are not runs ------------------------------------------


def _calibration_tree(runs_dir: Path) -> Path:
    """A calibration beside a run: template, both masters, three workers.

    Each holds a copy of the model's MF6 inputs (the whole model workspace until
    ledger 383), so each carries an ``mfsim.nam`` -- which is why they flooded
    the ledger.
    """

    pest = runs_dir / "baseline.pest"
    for sub in (
        "cal", "cal_ies_master", "cal_prior_master",
        "worker_0", "worker_1", "worker_2",
    ):
        workspace = pest / sub
        workspace.mkdir(parents=True)
        (workspace / "mfsim.nam").write_text("x", encoding="utf-8")
    (pest / "cal" / "myflopy_pest_metadata.json").write_text("{}", encoding="utf-8")
    return pest


def test_a_calibration_is_reported_not_listed_as_runs(tmp_path):
    """One calibration would otherwise add six rows to a ledger whose job is
    telling you what you ran. It is excluded, counted, and named in the repr --
    absent from the table, never silently absent."""

    project = mf.Project(tmp_path / "proj", name="p")
    simulation = mf.SimulationSpec(
        "s",
        models=[
            mf.gwf(
                "flow",
                packages=[
                    mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
                    mf.ic(strt=1.0),
                ],
            )
        ],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )
    project.add_simulation(simulation)
    project.prepare_run("baseline", simulation).write()
    _calibration_tree(project.runs_dir)

    runs = project.refresh_runs()

    assert list(runs) == ["baseline"]
    calibrations = runs.calibrations
    assert list(calibrations["calibration"]) == ["runs/baseline.pest"]
    assert list(calibrations["workspaces"]) == [6]     # ONE row, not six
    assert "1 calibration excluded" in repr(runs)
    assert "model.pest_runs" in runs._repr_html_()


def test_the_outermost_marker_wins_so_one_calibration_is_one_row(tmp_path):
    """The template carries the metadata marker and sits INSIDE the .pest root,
    so returning the first match up the tree split one calibration into two --
    the tree, and its own template counted separately."""

    from myflopy.run_record import pest_tree_root

    root = tmp_path
    template = root / "runs" / "baseline.pest" / "cal"
    template.mkdir(parents=True)
    (template / "myflopy_pest_metadata.json").write_text("{}", encoding="utf-8")

    assert pest_tree_root(template, root) == root / "runs" / "baseline.pest"
    assert pest_tree_root(root / "runs" / "baseline.pest" / "worker_0", root) == (
        root / "runs" / "baseline.pest"
    )
    assert pest_tree_root(root / "runs" / "baseline", root) is None


def test_the_legacy_pest_dirname_and_a_custom_workspace_are_both_caught(tmp_path):
    from myflopy.run_record import pest_tree_root

    legacy = tmp_path / "runs" / "m" / "pest" / "cal"
    legacy.mkdir(parents=True)
    assert pest_tree_root(legacy, tmp_path) == tmp_path / "runs" / "m" / "pest"

    custom = tmp_path / "scratch" / "elsewhere"
    custom.mkdir(parents=True)
    (custom / "myflopy_pest_metadata.json").write_text("{}", encoding="utf-8")
    assert pest_tree_root(custom, tmp_path) == custom


def test_the_pest_markers_match_the_pest_module():
    """They are DUPLICATED rather than imported -- see the note in run_record --
    so this is what keeps them honest."""

    from myflopy import run_record
    from myflopy.modflow.mf6.pest import runs as pest_runs

    assert run_record.PEST_ROOT_SUFFIX == pest_runs.PEST_ROOT_SUFFIX
    assert run_record.LEGACY_PEST_DIRNAME == pest_runs.LEGACY_PEST_DIRNAME
    assert run_record.PEST_METADATA == pest_runs.METADATA_FILENAME


def test_importing_myflopy_does_not_drag_in_the_pest_package():
    """The reason those constants are duplicated. Importing
    myflopy.modflow.mf6.pest.runs executes the pest package __init__, measured at
    1.43 s against an `import myflopy` that costs 0.02 s -- a 70x regression to
    reach three strings. A subprocess, because another test may already have
    imported pest into this one."""

    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-c",
         "import sys, myflopy; "
         "print('myflopy.modflow.mf6.pest' in sys.modules)"],
        capture_output=True, text=True, check=True,
    )

    assert result.stdout.strip() == "False", (
        "import myflopy now loads the pest package; the run_record constants "
        "exist to avoid exactly that"
    )


@pytest.mark.parametrize(
    "package",
    [
        lambda: mf.oc(saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
        lambda: mf.wel(
            stress_period_data={0: [[(0, 1), -1.0, "north"], [(0, 2), -2.0, "north"]]},
            boundnames=True,
        ),
        lambda: mf.wel(stress_period_data={0: [[(0, 1), np.float64(-1.5)]]}),
    ],
    ids=["shared-strings", "shared-boundname", "numpy-scalar"],
)
def test_a_json_round_trip_does_not_change_the_digest(package):
    """What Project.save / load does to a package, and why drift() cried 'edited'."""

    import json

    from myflopy.specs import PackageSpec, package_digest

    original = package()
    back = PackageSpec.from_dict(json.loads(json.dumps(original.to_dict())))
    assert package_digest(back) == package_digest(original)


def test_an_edited_value_still_changes_the_digest():
    """The normalization must not fold real differences away."""

    from myflopy.specs import package_digest

    assert package_digest(mf.wel(stress_period_data={0: [[(0, 1), np.float64(-1.5)]]})) != (
        package_digest(mf.wel(stress_period_data={0: [[(0, 1), -2.5]]}))
    )


def _gpkg_library_project(tmp_path):
    """A project whose CHD is resolved from a GeoPackage on disk."""

    import geopandas as gpd
    from shapely.geometry import box

    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        rectangular_voronoi,
    )

    grid = rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1))
    path = tmp_path / "bc.gpkg"
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds

    def write(head):
        gpd.GeoDataFrame(
            {"name": ["west"], "layer": [1], "head": [head]},
            geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)],
            crs=grid.crs,
        ).to_file(path, layer="chd", driver="GPKG")

    write(9.0)
    context = ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int))
    project = Project(tmp_path / "proj", name="p")
    project.add_package(
        "chd/west", mf.chd.gpkg(path, layer="chd", nper=1, head="head", context=context)
    )
    project.add_package("npf/base", mf.npf(k=1.0))
    model = mf.gwf(
        "m", context=context, packages=[mf.ref("chd/west"), mf.ref("npf/base")]
    )
    project.add_simulation(mf.simulation(model, name="s", packages=[mf.tdis(nper=1)]))
    project.prepare_run("r", "s", build=False)
    return project, path, write


def test_drift_sees_the_source_file_change_under_a_package(tmp_path):
    """`state` answers about the LIBRARY; `source` answers about the world outside it.

    A package resolved from a GeoPackage records `source_path` and nothing else, so
    editing that file left every run reporting `unchanged` -- and `drift()`'s own
    docstring promises `unchanged` means "reproducible from the library as it stands".
    It still is, which is why this is a SECOND column rather than a fifth state: the
    key means what it meant, and the file it was built from does not.
    """

    project, path, write = _gpkg_library_project(tmp_path)
    before = project.runs.drift().set_index("key")

    assert before.loc["chd/west", "source"] == "ok"
    assert before.loc["npf/base", "source"] == "—", "a package with no source"

    write(42.0)
    after = project.refresh_runs().drift().set_index("key")

    assert after.loc["chd/west", "source"] == "edited"
    assert after.loc["chd/west", "state"] == "unchanged", (
        "the LIBRARY key is genuinely unchanged -- collapsing the two would hide one"
    )

    path.unlink()
    gone = project.refresh_runs().drift().set_index("key")

    assert gone.loc["chd/west", "source"] == "missing"


def test_a_LAZY_package_reports_its_source_state_too(tmp_path):
    """The column was INVERTED, and the wrong half of the pair was blind (ledger 293).

    `_library_source_memo` keys off `metadata['source_path']`, which only the EAGER
    declaration stamped. So the package whose records are already frozen -- the one
    that rebuilds fine with its GeoPackage deleted -- reported `missing` when the file
    moved, while the LAZY one, which cannot build at all without that file, reported
    `—`: "names no source". Exactly backwards, and backwards in the direction that
    stays quiet.
    """

    import geopandas as gpd
    from shapely.geometry import box

    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        rectangular_voronoi,
    )

    grid = rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1))
    path = tmp_path / "bc.gpkg"
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds

    def write(head):
        gpd.GeoDataFrame(
            {"name": ["west"], "layer": [1], "head": [head]},
            geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)],
            crs=grid.crs,
        ).to_file(path, layer="chd", driver="GPKG")

    write(9.0)
    project = Project(tmp_path / "proj", name="p")
    # NO `context=` -- this is the lazy form, and the recipe it holds IS the path.
    project.add_package(
        "chd/west", mf.chd.gpkg(path, layer="chd", nper=1, head="head")
    )

    assert project.runs._library_source_memo()["chd/west"] == "ok"

    write(42.0)
    assert project.runs._library_source_memo()["chd/west"] == "edited"

    path.unlink()
    assert project.runs._library_source_memo()["chd/west"] == "missing"


def _bundled_project(tmp_path, head=9.0):
    """A project ARCHIVED through the public verb, with its author's tree deleted.

    Built by `save(bundle_sources=True)` rather than by hand. The first version of
    this fixture fabricated the state with `dataclasses.replace(spec, metadata=...)`
    -- and an adversarial review measured that the library could not actually PRODUCE
    that state, because bundling repointed the recipe and left the metadata behind.
    So the guard passed while the thing it guarded was unreachable: it certified.
    """

    import shutil

    import geopandas as gpd
    from shapely.geometry import box

    gis = tmp_path / "gis"
    gis.mkdir()
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [head]},
        geometry=[box(0, 0, 10, 10)], crs=2927,
    ).to_file(gis / "bc.gpkg", layer="chd", driver="GPKG")

    project = Project(tmp_path / "proj", name="p")
    project.add_package(
        "chd/west", mf.chd.gpkg(gis / "bc.gpkg", layer="chd", nper=1, head="head")
    )
    project.save(bundle_sources=True)
    shutil.rmtree(gis)                      # the archive is shipped; the author's tree is gone
    return project


def test_a_bundled_relative_source_is_read_from_the_project_not_the_cwd(tmp_path, monkeypatch):
    """`missing` must be a claim about the world, not about where the reader stands."""

    from myflopy.specs import SOURCE_SHA

    project = _bundled_project(tmp_path)

    assert project.package("chd/west").metadata["source_path"] == "inputs/bc.gpkg", (
        "bundling must repoint the metadata, or this test has nothing to anchor"
    )
    assert project.package("chd/west").metadata[SOURCE_SHA] is not None

    monkeypatch.chdir(tmp_path)  # anywhere but the project root

    assert project.runs._library_source_memo()["chd/west"] == "ok", (
        "the relative path must anchor to the project, not to os.getcwd()"
    )


def test_bundling_repoints_the_metadata_so_drift_watches_the_file_the_BUILD_reads(tmp_path):
    """The false clean bill, which is the worse half (found in review).

    Bundling moved the RECIPE into `inputs/` and left `metadata['source_path']` naming
    the original. `_library_source_memo` keys off the metadata -- so editing the
    bundled copy, the only file the build ever opens, read `source=ok`. That is the
    `unchanged`-is-the-dangerous-verdict failure occurring inside the column that was
    added to remove it.
    """

    import geopandas as gpd
    from shapely.geometry import box

    project = _bundled_project(tmp_path)

    assert project.runs._library_source_memo()["chd/west"] == "ok"

    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [42.0]},
        geometry=[box(0, 0, 10, 10)], crs=2927,
    ).to_file(tmp_path / "proj" / "inputs" / "bc.gpkg", layer="chd", driver="GPKG")

    assert project.runs._library_source_memo()["chd/west"] == "edited", (
        "the column must watch the bundled copy, which is what the build opens"
    )


def test_a_package_declared_before_fingerprints_reports_unfingerprinted(tmp_path):
    """A missing fingerprint is not evidence of sameness -- the module's own rule."""

    import dataclasses

    from myflopy.specs import SOURCE_SHA

    project, _, _ = _gpkg_library_project(tmp_path)
    old = project.package("chd/west")
    project.add_package(
        "chd/west",
        dataclasses.replace(
            old, metadata={k: v for k, v in old.metadata.items() if k != SOURCE_SHA}
        ),
        overwrite=True,
    )

    row = project.refresh_runs().drift().set_index("key").loc["chd/west"]

    assert row["source"] == "unfingerprinted"


def test_recording_a_source_fingerprint_moves_no_digest(tmp_path):
    """Load-bearing: `metadata` is excluded from `package_digest` (ledgers 186, 250).

    If it were not, adding this would report every stored `.gpkg` package as drifted
    on the day it landed -- which is precisely why ledger 250 kept the grid pin out.
    """

    import dataclasses

    from myflopy.specs import SOURCE_SHA, package_digest

    project, _, _ = _gpkg_library_project(tmp_path)
    spec = project.package("chd/west")

    assert spec.metadata[SOURCE_SHA] is not None
    without = dataclasses.replace(
        spec, metadata={k: v for k, v in spec.metadata.items() if k != SOURCE_SHA}
    )

    assert package_digest(spec) == package_digest(without)


def test_an_unreadable_source_does_not_stop_a_package_being_declared(tmp_path):
    """A bookkeeping column must never decide whether a model can be built."""

    from myflopy.specs import source_fingerprint

    assert source_fingerprint(tmp_path / "nope.gpkg") is None
    assert source_fingerprint(tmp_path) is None, "a directory is not readable as bytes"


# ------------------------------------ what the RUN read, not what the library says


def _sourced_project(tmp_path, head=9.0, *, inline=False):
    """A project whose CHD resolves LAZILY from a GeoPackage, by ref or inline."""

    import geopandas as gpd
    from shapely.geometry import box

    gis = tmp_path / "gis"
    gis.mkdir(exist_ok=True)

    def write(value):
        gpd.GeoDataFrame(
            {"name": ["west"], "layer": [1], "head": [value]},
            geometry=[box(0, 0, 600, 600)], crs=2927,
        ).to_file(gis / "bc.gpkg", layer="chd", driver="GPKG")

    write(head)
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(
        gis / "b.gpkg"
    )

    project = Project(tmp_path / "proj", name="p")
    project.add_grid(
        "g",
        mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=gis / "b.gpkg"), max_area=40000.0
        ),
    )
    chd = mf.chd.gpkg(gis / "bc.gpkg", layer="chd", nper=1, head="head")
    if not inline:
        project.add_package("chd/west", chd)
    stack = mf.LayerStack(top=mf.Surface.flat(100.0)).add("sand", thickness=50)
    model = mf.gwf(
        "m",
        grid=mf.grid_ref("g"),
        packages=[
            mf.disv(stack), mf.npf(k=1.0), mf.ic(strt=95.0),
            mf.oc(head_filerecord="m.hds"),
            chd if inline else mf.ref("chd/west"),
        ],
    )
    project.add_simulation(
        mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    )
    return project, write


def _head_of(run):
    """The one distinct head the CHD wrote (the polygon covers every cell)."""

    heads = run.model("m").gwf.get_package("chd").stress_period_data.get_data(0)["head"]
    return sorted(set(heads.tolist()))


def test_two_runs_either_side_of_an_edit_are_DISTINGUISHABLE(tmp_path):
    """The correctness bug (ledger 294): the library's digest is the same for both.

    `package_digests` records what the LIBRARY declared, and both runs reference one
    declaration -- so two runs that built materially different models recorded the
    identical fingerprint and `drift()` gave them the same verdict, INCLUDING the run
    that agrees with the file. `built_from` records what each run actually READ.
    """

    project, write = _sourced_project(tmp_path)

    before = project.prepare_run("before", "s")
    assert _head_of(before) == pytest.approx([9.0])

    write(42.0)
    after = project.prepare_run("after", "s")
    assert _head_of(after) == pytest.approx([42.0])

    assert before.package_digests == after.package_digests, (
        "the premise: one library declaration, so the digests CANNOT differ"
    )
    assert before.built_from != after.built_from

    drift = project.refresh_runs().drift().set_index("run")

    assert drift.loc["after", "source"] == "ok", "this run agrees with the file"
    assert drift.loc["before", "source"] == "edited", "and this one does not"


def test_an_INLINE_lazy_package_gets_a_row_and_a_verdict(tmp_path):
    """`_package_digests` records only `mf.ref` entries, so inline had NO row.

    Its docstring's premise -- "a package written inline is baked into the run's own
    recipe and cannot drift" -- was true until lazy packages, which hold a PATH and
    open it at build. Measured before the fix: edit the file and drift reported
    `unchanged` with source `—` ("names no source") while the rebuild gave 42.0.
    """

    project, write = _sourced_project(tmp_path, inline=True)
    run = project.prepare_run("r", "s")

    assert run.package_digests == {}, "the premise: nothing is recorded for inline"
    assert "chd" in run.built_from

    clean = project.refresh_runs().drift().set_index("key")
    assert clean.loc["chd", "source"] == "ok"
    assert clean.loc["chd", "state"] == "—", (
        "there is no library key, so `state` has nothing to answer about"
    )

    write(42.0)
    moved = project.refresh_runs().drift().set_index("key")

    assert moved.loc["chd", "source"] == "edited"


def test_an_EAGER_package_is_not_re_hashed_at_build(tmp_path):
    """Only a lazy package reads a file at build, so only a lazy one is re-hashed.

    Re-reading an eager package's GeoPackage would describe a read that never
    happened -- and it is what keeps the cost honest: hashing is 0.15% of the work it
    describes for an ordinary GeoPackage, but a full read of a multi-gigabyte DEM is
    ~1.9 s, and the build only samples such a raster.
    """

    import geopandas as gpd
    from shapely.geometry import box

    from myflopy.specs import SOURCE_SHA, ModelContext

    project, _write = _sourced_project(tmp_path)
    grid = project.grid("g").resolve(workspace=tmp_path / "gridws")
    context = ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int))
    eager = mf.chd.gpkg(
        tmp_path / "gis" / "bc.gpkg", layer="chd", nper=1, head="head", context=context
    )
    project.add_package("chd/west", eager, overwrite=True)

    run = project.prepare_run("r", "s")
    recorded = run.built_from["chd/west"]

    assert list(recorded.values()) == [eager.metadata[SOURCE_SHA]], (
        "an eager package carries its DECLARATION hash -- it read nothing at build"
    )

    # And it still answers, because that hash is what this run was built from.
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [42.0]},
        geometry=[box(0, 0, 600, 600)], crs=2927,
    ).to_file(tmp_path / "gis" / "bc.gpkg", layer="chd", driver="GPKG")

    assert project.refresh_runs().drift().set_index("key").loc["chd/west", "source"] == (
        "edited"
    )


def test_a_package_naming_SEVERAL_files_reports_the_worst_with_a_count(tmp_path):
    """A lazy discretization names one file per surface plus the top.

    One cell per package, so it must not read clean because two of three files held
    still. Worst-first, because `missing` outranks `edited`.
    """

    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import box

    gis = tmp_path / "gis"
    gis.mkdir()

    def raster(name, value):
        path = gis / name
        with rasterio.open(
            path, "w", driver="GTiff", height=20, width=20, count=1, dtype="float32",
            crs="EPSG:2927", transform=from_origin(0, 600, 30, 30),
        ) as handle:
            handle.write(np.full((20, 20), value, dtype="float32"), 1)
        return path

    top, mid = raster("top.tif", 100.0), raster("mid.tif", 60.0)
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(
        gis / "b.gpkg"
    )

    project = Project(tmp_path / "proj", name="p")
    project.add_grid(
        "g",
        mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=gis / "b.gpkg"), max_area=40000.0
        ),
    )
    stack = mf.LayerStack(top=mf.Surface.raster(top)).add(
        "sand", bottom=mf.Surface.raster(mid)
    )
    model = mf.gwf(
        "m", grid=mf.grid_ref("g"),
        packages=[mf.disv(stack), mf.npf(k=1.0), mf.ic(strt=95.0),
                  mf.oc(head_filerecord="m.hds")],
    )
    project.add_simulation(
        mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    )
    run = project.prepare_run("r", "s")

    assert len(run.built_from["disv"]) == 2, "two rasters, one cell"
    assert project.refresh_runs().drift().set_index("key").loc["disv", "source"] == "ok"

    raster("mid.tif", 55.0)
    assert project.refresh_runs().drift().set_index("key").loc["disv", "source"] == (
        "edited (1 of 2)"
    )

    top.unlink()
    assert project.refresh_runs().drift().set_index("key").loc["disv", "source"] == (
        "missing (1 of 2)"
    ), "worst-first: `missing` outranks `edited`"


def test_built_from_survives_the_manifest_round_trip(tmp_path):
    """It rides in `metadata`, which round-trips wholesale.

    A top-level manifest key would be dropped by any rewrite from code that does not
    know it -- and `_record_status` rewrites on every status change, so the record
    would vanish the first time a run was written or executed.
    """

    project, _write = _sourced_project(tmp_path)
    run = project.prepare_run("r", "s")
    recorded = run.built_from

    assert recorded
    run._record_status("written")          # the rewrite that drops a top-level key

    written = json.loads((run.workspace / "run.json").read_text())

    assert written["metadata"][Run.BUILT_FROM] == recorded
    assert Run.BUILT_FROM not in written, "it must NOT also sit at the top level"
