"""A run's RECORD is not its FILES, and neither replaces the other in silence.

Ledger 611. A run folder holds the record -- ``run.spec``, its frozen snapshot,
its library copies and digests -- which a rebuild, ``drift()`` and the write-backs read;
and the MF6 input files, which are what MF6 runs. Before this:

* ``execute()`` on a completed run that held a simulation edited through FloPy ran the
  OLD files (S1 of the design's reproducer);
* ``write()`` of such an edit left the record saying the old value with nothing marked;
* a run reopened without loading -- ``project.runs[key]``, ``load_run(load=False)`` --
  REBUILT its record over its files on ``execute()``, reverting a written edit and a
  hand edit alike (S4-S6).

Every model here is a 2x2 steady grid with no boundary, so its heads equal ``ic.strt``
and the head file says which ``strt`` MF6 ran.
"""

from __future__ import annotations

import json
import os
import shutil
import time
import warnings

import flopy
import pytest

import myflopy as mf
from myflopy import workspace as workspace_module
from myflopy.run_record import (
    InputsVerdict,
    _same_token,
    _scan_simulation_inputs,
    written_differences,
)

needs_mf6 = pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")


def _simulation(strt: float = 1.0) -> mf.SimulationSpec:
    return mf.SimulationSpec(
        "baseline",
        models=[mf.gwf("flow", packages=[
            mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
            mf.ic(strt=strt),
            mf.ref("npf/base"),
            mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
        ])],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )


def _executed(tmp_path):
    """A project with one run, ``baseline``, written and executed: ``(project, run)``."""

    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    run = project.prepare_run("baseline", _simulation())
    ok, report = run.execute()
    assert ok, "\n".join(report[-20:])
    return project, run


def _heads(run) -> float:
    heads = flopy.utils.HeadFile(run.workspace / "flow.hds")
    try:
        return float(heads.get_data()[0, 0, 0])
    finally:
        heads.close()  # an unclosed reader warns inside another test's gc.collect()


def _set_strt(run, value: float) -> None:
    run.simulation.get_model("flow").ic.strt.set_data(value)


def _hand_edit(run, old: str, new: str, name: str = "flow.ic") -> None:
    path = run.workspace / name
    text = path.read_text()
    assert old in text
    path.write_text(text.replace(old, new))


def _write_edited(run, value: float = 50.0) -> None:
    _set_strt(run, value)
    with pytest.warns(UserWarning, match="edited through FloPy"):
        run.write(silent=True)


# -- what execute() runs ---------------------------------------------------------------


@needs_mf6
def test_execute_runs_an_in_memory_edit_of_a_completed_run(tmp_path):
    """S1: the held simulation was edited after the run completed; MF6 ran the old files."""

    _, run = _executed(tmp_path)
    _set_strt(run, 50.0)
    with pytest.warns(UserWarning, match="edited through FloPy after it was built"):
        run.execute()
    assert _heads(run) == pytest.approx(50.0)


@needs_mf6
def test_writing_a_flopy_edit_marks_the_files_edited_and_warns(tmp_path):
    """S2: written, the edit is in the files and NOT in the record -- and now it says so."""

    _, run = _executed(tmp_path)
    _write_edited(run)
    verdict = run.inputs
    assert verdict.verdict == "edited"
    assert verdict.departs == {"flow.ic": "edited through FloPy after run.build()"}
    assert verdict.matches is False
    # Persisted: a later session reads the same mark, on either reopen.
    assert mf.load_run(run.workspace, load=False).inputs.departs == verdict.departs
    assert mf.load_run(run.workspace).inputs.verdict == "edited"
    # Writing it again is not a NEW edit, so it does not warn twice.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        run.write(silent=True)


@needs_mf6
@pytest.mark.parametrize("door", ["load_run", "runs", "reopen_run"])
def test_a_run_reopened_without_loading_runs_its_files_not_its_record(tmp_path, door):
    """S4/S5: a reopened run's execute() rebuilt its RECORD over its files."""

    project, run = _executed(tmp_path)
    _write_edited(run)
    frozen = run.workspace / ".myflopy" / "frozen.npz"
    stamp = frozen.stat().st_mtime_ns
    if door == "load_run":
        reopened = mf.load_run(run.workspace, load=False)
    elif door == "runs":
        reopened = mf.Project(tmp_path / "proj", name="p").runs["baseline"]
    else:
        reopened = project.reopen_run("baseline", load=False)
    assert reopened.simulation is None
    reopened.execute()
    assert _heads(reopened) == pytest.approx(50.0)
    assert "50.00000000" in (run.workspace / "flow.ic").read_text()
    assert frozen.stat().st_mtime_ns == stamp  # nothing was rebuilt or re-frozen
    assert reopened.inputs.verdict == "edited"


@needs_mf6
def test_a_hand_edited_input_reads_edited_and_is_never_reverted(tmp_path):
    """S6 -- ss's route: a file changed outside myflopy after the write."""

    _, run = _executed(tmp_path)
    _hand_edit(run, "1.00000000", "77.00000000")
    verdict = mf.load_run(run.workspace, load=False).inputs
    assert verdict.verdict == "edited"
    assert verdict.departs == {"flow.ic": "changed outside myflopy after its last write"}
    reopened = mf.load_run(run.workspace, load=False)
    reopened.execute()
    assert _heads(reopened) == pytest.approx(77.0)
    assert "77.00000000" in (run.workspace / "flow.ic").read_text()


@needs_mf6
def test_write_with_nothing_held_refuses_over_files_that_depart(tmp_path):
    """Building the record over files not known to match it is a choice, so it is asked for."""

    _, run = _executed(tmp_path)
    _hand_edit(run, "1.00000000", "77.00000000")
    reopened = mf.load_run(run.workspace, load=False)
    with pytest.raises(ValueError, match=r"run\.open\(\).*run\.build\(\) then run\.write\(\)"):
        reopened.write(silent=True)
    assert "77.00000000" in (run.workspace / "flow.ic").read_text()
    # The one spelling that puts the record over the files.
    reopened.build()
    reopened.write(silent=True)
    reopened.execute()
    assert reopened.inputs.verdict == "ok"
    assert _heads(reopened) == pytest.approx(1.0)


@needs_mf6
def test_an_unedited_run_reads_ok_and_reexecutes_without_rewriting(tmp_path):
    _, run = _executed(tmp_path)
    assert run.inputs == InputsVerdict("ok", {}, run.inputs.written_at)
    assert run.inputs.written_at is not None
    names = _scan_simulation_inputs(run.workspace)[0]
    before = {name: (run.workspace / name).stat().st_mtime_ns for name in names}
    time.sleep(0.01)
    run.execute()
    mf.load_run(run.workspace, load=False).execute()
    assert {name: (run.workspace / name).stat().st_mtime_ns for name in names} == before


@needs_mf6
def test_a_run_written_before_fingerprints_is_unfingerprinted_not_ok(tmp_path):
    """No record of what was written is not evidence of agreement -- the full check asks."""

    _, run = _executed(tmp_path)
    manifest = run.workspace / "run.json"
    data = json.loads(manifest.read_text())
    data["metadata"].pop("inputs")
    manifest.write_text(json.dumps(data))
    reopened = mf.load_run(run.workspace, load=False)
    assert reopened.inputs.verdict == "unfingerprinted"
    assert reopened.inputs.matches is None
    # The full check builds and writes the record beside the files: they agree.
    assert reopened._record_departures() == {}
    _hand_edit(run, "1.00000000", "77.00000000")
    departs = mf.load_run(run.workspace, load=False)._record_departures()
    assert list(departs) == ["flow.ic"]
    assert "77.00000000 (files) vs 1.00000000 (record)" in departs["flow.ic"]
    with pytest.raises(ValueError, match="unfingerprinted"):
        mf.load_run(run.workspace, load=False).write(silent=True)


@needs_mf6
def test_files_changed_under_a_held_simulation_are_refused(tmp_path):
    _, run = _executed(tmp_path)
    _hand_edit(run, "1.00000000", "77.00000000")
    with pytest.raises(ValueError, match=r"flow\.ic.*run\.open\(\) then run\.execute\(\)"):
        run.execute()
    run.open()
    run.execute()
    assert _heads(run) == pytest.approx(77.0)


@needs_mf6
def test_an_external_array_edit_is_seen(tmp_path):
    """An OPEN/CLOSE file is an input too: the scan follows it from its package file."""

    _, run = _executed(tmp_path)
    run.simulation.set_all_data_external()
    run.write(silent=True)  # same values, stored elsewhere: not an edit
    assert run.inputs.verdict == "ok"
    external = [
        name.as_posix() for name in _scan_simulation_inputs(run.workspace)[0]
        if "strt" in name.as_posix()
    ]
    assert external, _scan_simulation_inputs(run.workspace)[0]
    path = run.workspace / external[0]
    path.write_text(path.read_text().replace("1.00000000", "9.00000000"))
    verdict = run.inputs
    assert verdict.verdict == "edited"
    assert verdict.departs == {external[0]: "changed outside myflopy after its last write"}
    owners = _scan_simulation_inputs(run.workspace)[4]
    assert owners[path.relative_to(run.workspace)].as_posix() == "flow.ic"


@needs_mf6
def test_an_opened_runs_edit_is_marked_after_open_and_carries_what_it_held(tmp_path):
    """The files lineage: departures found at open() ride into the next write's mark."""

    _, run = _executed(tmp_path)
    _hand_edit(run, "1.00000000", "77.00000000")
    reopened = mf.load_run(run.workspace)
    reopened.simulation.get_model("flow").ic.strt.set_data(60.0)
    with pytest.warns(UserWarning, match="after it was opened"):
        reopened.write(silent=True)
    # Both of what happened to flow.ic, in order.
    assert reopened.inputs.departs == {"flow.ic": (
        "changed outside myflopy after its last write, "
        "then edited through FloPy after run.open()")}
    reopened.simulation.get_model("flow").npf.k.set_data(3.0)
    with pytest.warns(UserWarning, match="flow.npf"):
        reopened.write(silent=True)
    assert reopened.inputs.departs["flow.npf"] == "edited through FloPy after run.open()"


# -- the fingerprint -------------------------------------------------------------------


@needs_mf6
def test_the_fingerprint_moves_only_the_edited_file(tmp_path):
    _, run = _executed(tmp_path)
    before = workspace_module._simulation_fingerprint(run.simulation)
    assert before == workspace_module._simulation_fingerprint(run.simulation)
    _set_strt(run, 2.0)
    after_ic = workspace_module._simulation_fingerprint(run.simulation)
    assert [name for name in before if before[name] != after_ic[name]] == ["flow.ic"]
    run.simulation.tdis.perioddata.set_data([(2.0, 1, 1.0)])
    after_tdis = workspace_module._simulation_fingerprint(run.simulation)
    assert [name for name in before if after_ic[name] != after_tdis[name]] == ["baseline.tdis"]


def test_the_fingerprint_covers_every_child_package_and_is_stable(canonical_model):
    """OBS children sit in a model's packagelist; every file the model writes is keyed."""

    sim = canonical_model.sim
    first = workspace_module._simulation_fingerprint(sim)
    assert first == workspace_module._simulation_fingerprint(sim)
    assert any(name.endswith(".obs") for name in first), sorted(first)
    assert "mfsim.nam" in first


@needs_mf6
def test_a_list_package_on_a_real_grid_is_not_edited_by_its_own_first_write(
        tmp_path, triangulated_grid):
    """FloPy fills MAXBOUND at write; the fingerprint must not call that an edit.

    Measured on a Ten Trails run: five list packages' fingerprints moved across the first
    write until the fingerprint applied FloPy's size rule itself. On a REAL grid (ledgers
    263/265), so the record also freezes, thaws and rebuilds a triangulated mesh.
    """

    grid = triangulated_grid
    props = grid.get_disv_gridprops()
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    model = mf.gwf("flow", context=mf.ModelContext(grid=grid), packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10.0, botm=[0.0]),
        mf.ic(strt=5.0),
        mf.ref("npf/base"),
        mf.chd(stress_period_data={0: [[(0, 0), 5.0], [(0, grid.ncpl - 1), 4.0]]}),
        mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    run = project.prepare_run("r", simulation)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        ok, report = run.execute()
        assert ok, "\n".join(report[-20:])
        # The write filled MAXBOUND; a second execute must neither call that an edit
        # nor rewrite the files.
        stamp = (run.workspace / "flow.chd").stat().st_mtime_ns
        time.sleep(0.01)
        run.execute()
        assert (run.workspace / "flow.chd").stat().st_mtime_ns == stamp
        run.write(silent=True)
    assert run.inputs.verdict == "ok"
    assert run.inputs.departs == {}
    reopened = mf.load_run(run.workspace, load=False)
    reopened.metadata.pop("inputs")
    assert reopened._record_departures(full=True) == {}


# -- the comparison --------------------------------------------------------------------


def test_numbers_compare_by_value_to_the_coarser_precision_and_words_by_case():
    assert _same_token("1", "1.00000000")
    assert _same_token("1.0D+00", "1.00000000")
    assert _same_token("1.10000002", "1.1")
    assert not _same_token("1", "31.00000000")
    assert not _same_token("1.00000000", "1.00000002")
    assert _same_token("BEGIN", "begin")
    assert not _same_token("days", "seconds")


def test_written_differences_names_the_file_block_and_values(tmp_path):
    files, record = tmp_path / "files", tmp_path / "record"
    for root, perioddata, extra in (
        (files, "  1  1  1.10000002\n", True),
        (record, "  31.00000000  10  1.10000002\n", False),
    ):
        root.mkdir()
        (root / "mfsim.nam").write_text(
            "# File generated by Flopy at whatever time\n"
            "BEGIN timing\n  TDIS6  s.tdis\nEND timing\n"
            + ("BEGIN models\n  gwf6  extra.nam  extra\nEND models\n" if extra else "")
        )
        (root / "s.tdis").write_text(
            "# a header that differs\nBEGIN dimensions\n  NPER 1\nEND dimensions\n"
            f"BEGIN perioddata\n{perioddata}END perioddata\n"
        )
    (files / "extra.nam").write_text("BEGIN packages\nEND packages\n")
    rows = written_differences(files, record)
    assert ("s.tdis", "PERIODDATA",
            "2 of 3 values differ: 1 (files) vs 31.00000000 (record); "
            "1 (files) vs 10 (record)") in rows
    assert ("extra.nam", "", "only in the files") in rows
    assert any(row[0] == "mfsim.nam" and "written differently" in row[2] for row in rows)
    assert written_differences(files, files) == []


def test_the_scan_names_each_inputs_owner(tmp_path):
    (tmp_path / "mfsim.nam").write_text(
        "BEGIN timing\n  TDIS6  s.tdis\nEND timing\n"
        "BEGIN models\n  gwf6  m.nam  m\nEND models\n"
    )
    (tmp_path / "s.tdis").write_text("BEGIN perioddata\n 1 1 1\nEND perioddata\n")
    (tmp_path / "m.nam").write_text(
        "BEGIN packages\n  NPF6  m.npf  npf\n  CHD6  m.chd  chd\nEND packages\n"
    )
    (tmp_path / "m.npf").write_text("BEGIN griddata\n k\n  OPEN/CLOSE  'k.txt'\nEND griddata\n")
    (tmp_path / "m.chd").write_text("BEGIN options\n TS6 FILEIN m.ts\nEND options\n")
    (tmp_path / "m.ts").write_text("BEGIN timeseries\n OPEN/CLOSE ts.dat\nEND timeseries\n")
    for name in ("k.txt", "ts.dat"):
        (tmp_path / name).write_text("1.0\n")
    owners = {k.as_posix(): v.as_posix() for k, v in _scan_simulation_inputs(tmp_path)[4].items()}
    assert owners == {
        "mfsim.nam": "mfsim.nam", "s.tdis": "s.tdis", "m.nam": "m.nam", "m.npf": "m.npf",
        "m.chd": "m.chd", "k.txt": "m.npf", "m.ts": "m.chd", "ts.dat": "m.chd",
    }


def test_a_bare_workspace_has_no_record(tmp_path):
    run = workspace_module.Run(name="bare", workspace=tmp_path, status="loaded")
    assert run.inputs.verdict == "no record"
    assert run._record_departures() == {}


@needs_mf6
def test_a_prepared_run_is_unwritten_until_written(tmp_path):
    project = mf.Project(tmp_path / "proj", name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    run = project.prepare_run("baseline", _simulation())
    assert run.inputs.verdict == "unwritten"
    run.write(silent=True)
    assert run.inputs.verdict == "ok"
    os.utime(run.workspace / "flow.ic")  # a touch is not an edit: content is hashed
    assert run.inputs.verdict == "ok"
