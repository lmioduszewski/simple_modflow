"""A PEST build deletes its ``workspace=`` only when a build wrote it there (ledger 360).

``model.pest(name, workspace=W)`` makes ``W`` the TEMPLATE directory, and ``build()``
deleted it whatever it held -- ``_clear_stale_template`` ran ``shutil.rmtree(W)`` and
pyEMU's ``PstFrom(remove_existing=True)`` removed it again. So a ``W`` that was an
existing folder went with everything in it: measured, ``W`` = the Project root holding
the model's own run was deleted, and the build then failed with "No MF6 simulation in
.../runs/r". The staging directory beside it (``W.staging``, ledger 356) was deleted
the same way.

A build now deletes each of the two only when it is absent, empty, or one a previous
build wrote there in that role -- marked by the ``.myflopy_pest_build.json`` the build
writes into it, naming the directory -- and refuses ``W`` (or its staging directory)
that is, contains or sits inside the run directory whatever it holds. Everything is
checked before anything is deleted, and ``PstFrom`` deletes nothing itself.

The model: a 4 x 4 irregular Voronoi grid of 10 ft cells, one layer 10 - 0 ft, K 10
ft/d, CHD at 9 ft west and 8 ft east, one steady period, run through
``mf.Project(root / "p").run("r", ...)`` -- the run directory is ``root/p/runs/r``.
Every test runs its own copy under its own ``tmp_path``.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.project import (
    BUILD_MARKER,
    METADATA_FILENAME,
    PestProject,
    _write_build_marker,
)

pyemu = pytest.importorskip("pyemu")


def _simulation():
    grid = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=10.0, refined=False)
    props = grid.get_disv_gridprops()
    x = np.array([cell[1] for cell in props["cell2d"]])
    west, east = np.flatnonzero(x < x.min() + 5.0), np.flatnonzero(x > x.max() - 5.0)
    chd = [[(0, int(c)), 9.0] for c in west] + [[(0, int(c)), 8.0] for c in east]
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=8.5),
        mf.chd(stress_period_data={0: chd}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _model(root: Path):
    """The model in the module docstring, run under ``root``; its run model."""

    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", _simulation())
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _calibration(model, workspace: Path | None, *, name: str = "cal") -> PestProject:
    cal = model.pest(name, workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "cell": [5], "layer": [0]}),
        values=pd.DataFrame({"time": [0], "W": [8.6]}),
    ))
    return cal


def _build(cal: PestProject):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{cal.name}.pst", noptmax=0)


def _tree(directory) -> dict[str, str]:
    """Every file under ``directory``, by relative path, to the SHA-256 of its bytes."""

    directory = Path(directory)
    return {str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(directory.rglob("*")) if path.is_file()}


def _assert_built_at(cal: PestProject, workspace: Path):
    """A whole template at ``workspace``, marked as the build's, and no staging copy left."""

    workspace = Path(workspace)
    assert (workspace / f"{cal.name}.pst").is_file() and (workspace / "mfsim.nam").is_file()
    marker = json.loads((workspace / BUILD_MARKER).read_text(encoding="utf-8"))
    assert marker["directory"] == str(workspace.resolve()) and marker["role"] == "template"
    assert not cal._staging_workspace().exists()


# -- refused: a folder no build wrote, and anything over the run ---------------------------


def test_a_project_root_as_the_workspace_is_refused_and_left_byte_identical(tmp_path):
    """The reviewer's case: ``W`` = the Project root holding the model's own run.
    It was deleted -- the run with it -- and the build then died on the missing run.
    Refused when the project is made, naming the default, and again by ``build()`` if
    the attribute is pointed there afterwards; not one byte of the project moves."""

    model = _model(tmp_path)
    root = tmp_path / "p"
    before = _tree(root)

    with pytest.raises(ValueError, match=r"CONTAINS the model's run directory") as refused:
        model.pest("cal", workspace=root, start_datetime="2024-01-01")
    message = str(refused.value)
    assert "model.pest('cal', ...)" in message and str(model.workspace.parent / "r.pest" / "cal") in message
    assert "workspace=<a new or empty folder>" in message

    cal = _calibration(model, tmp_path / "elsewhere")
    cal.template_workspace = root
    with pytest.raises(ValueError, match=r"CONTAINS the model's run directory"):
        _build(cal)

    assert _tree(root) == before
    assert not (tmp_path / "p.staging").exists() and not (tmp_path / "elsewhere").exists()


def test_a_non_empty_folder_no_build_wrote_is_refused_and_left_as_it_was(tmp_path):
    """A data folder beside the project: its files were deleted and the template
    built over them. Refused, naming the folder, what it holds and the spellings."""

    model = _model(tmp_path)
    data = tmp_path / "data"
    data.mkdir()
    (data / "wells.csv").write_text("name,x,y\nW,1,2\n", encoding="utf-8")
    (data / "notes").mkdir()
    (data / "notes" / "field.txt").write_text("field notes\n", encoding="utf-8")
    # A metadata file alone is not this calibration's template: it names no project.
    (data / METADATA_FILENAME).write_text("{}", encoding="utf-8")
    before, run_before = _tree(data), _tree(model.workspace)

    cal = _calibration(model, data)
    with pytest.raises(FileExistsError, match=r"already holds 3 entries") as refused:
        _build(cal)
    message = str(refused.value)
    assert f"workspace={str(data)!r}" in message and "notes" in message and "wells.csv" in message
    assert "model.pest('cal', ...)" in message and "delete it yourself" in message

    assert _tree(data) == before
    assert _tree(model.workspace) == run_before
    assert not cal._staging_workspace().exists()


def test_a_folder_holding_the_run_is_refused_even_marked_as_a_template(tmp_path):
    """``W`` that contains the run, or sits inside it, is refused whatever it holds:
    a build's marker -- forged here, naming the folder -- does not make it safe."""

    model = _model(tmp_path)
    run = Path(model.workspace)
    for workspace, relation in ((tmp_path, "CONTAINS"), (run / "pest" / "cal", "sits INSIDE"),
                                (run, "IS")):
        workspace.mkdir(parents=True, exist_ok=True)
        _write_build_marker(workspace, "template")
        before = _tree(tmp_path)
        with pytest.raises(ValueError, match=f"{relation} the model's run directory"):
            model.pest("cal", workspace=workspace, start_datetime="2024-01-01")
        assert _tree(tmp_path) == before, workspace


def test_a_workspace_whose_staging_directory_is_the_run_is_refused(tmp_path):
    """The staging copy goes beside ``W`` as ``W.staging``, and the build deleted
    whatever was there: a run directory named ``m.staging``, built at ``W = m``, was
    deleted before it was copied. Refused, and the run is left as it was."""

    run = tmp_path / "m.staging"
    shutil.copytree(Path(_model(tmp_path / "source").workspace), run)
    model = mf.load_mf6_run(run, verbosity_level=0)
    before = _tree(run)

    with pytest.raises(ValueError, match=r"its staging directory .*IS the model's run directory"):
        model.pest("cal", workspace=tmp_path / "m", start_datetime="2024-01-01")
    assert _tree(run) == before


def test_a_user_folder_where_the_staging_copy_goes_is_refused(tmp_path):
    """``W`` absent, ``W.staging`` a folder of the user's: it was deleted, used and
    removed. Refused before ``W`` is made."""

    model = _model(tmp_path)
    staging = tmp_path / "tmpl.staging"
    staging.mkdir()
    (staging / "mine.txt").write_text("the user's\n", encoding="utf-8")

    cal = _calibration(model, tmp_path / "tmpl")
    with pytest.raises(FileExistsError, match=r"the staging directory beside it"):
        _build(cal)
    assert (staging / "mine.txt").read_text(encoding="utf-8") == "the user's\n"
    assert not (tmp_path / "tmpl").exists()


def test_a_folder_that_appears_where_the_staging_copy_goes_is_not_deleted(tmp_path, monkeypatch):
    """Between the check and the copy a folder of the user's appears at ``W.staging``:
    the build makes the staging directory with ``mkdir``, which refuses, and its
    ``finally`` removes only a staging copy carrying the build's own marker -- so the
    folder survives the failed build. Both used ``shutil.rmtree``."""

    model = _model(tmp_path)
    staging = tmp_path / "tmpl.staging"
    clear = PestProject._clear_stale_template

    def _then_a_folder_appears(self):
        clear(self)
        staging.mkdir()
        (staging / "mine.txt").write_text("the user's\n", encoding="utf-8")

    monkeypatch.setattr(PestProject, "_clear_stale_template", _then_a_folder_appears)
    with pytest.raises(FileExistsError):
        _build(_calibration(model, tmp_path / "tmpl"))
    assert (staging / "mine.txt").read_text(encoding="utf-8") == "the user's\n"


def test_a_copy_of_a_template_is_not_taken_for_one(tmp_path):
    """The marker names the directory it was written for, so a COPY of a template --
    what pyEMU's ``start_workers`` makes of it as the IES master, or a backup -- is
    refused, not cleared."""

    model = _model(tmp_path)
    _build(_calibration(model, tmp_path / "cal"))
    master = tmp_path / "cal_ies_master"
    shutil.copytree(tmp_path / "cal", master)
    (master / "cal.phi.actual.csv").write_text("iteration,mean\n0,1.0\n", encoding="utf-8")
    before = _tree(master)

    with pytest.raises(FileExistsError, match=r"no myflopy PEST build wrote it"):
        _build(_calibration(model, master))
    assert _tree(master) == before


def test_pstfrom_itself_deletes_nothing(tmp_path, monkeypatch):
    """The second line: with the guard bypassed, ``PstFrom`` -- made with
    ``remove_existing=False`` -- refuses an occupied ``W`` rather than deleting it."""

    model = _model(tmp_path)
    data = tmp_path / "data"
    data.mkdir()
    (data / "keep.txt").write_text("keep\n", encoding="utf-8")
    monkeypatch.setattr(PestProject, "_clear_stale_template", lambda self: None)

    with pytest.raises(Exception, match="already exists"):
        _build(_calibration(model, data))
    assert (data / "keep.txt").read_text(encoding="utf-8") == "keep\n"


def test_a_workspace_that_is_a_file_is_refused(tmp_path):
    model = _model(tmp_path)
    path = tmp_path / "cal"
    path.write_text("not a folder\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="is a FILE"):
        _build(_calibration(model, path))
    assert path.read_text(encoding="utf-8") == "not a folder\n"


# -- used: absent, empty, and a template a build wrote ------------------------------------


def test_a_template_a_previous_build_wrote_is_rebuilt(tmp_path):
    """Re-running the calibration -- a NEW ``model.pest`` on the same ``W``, as a
    re-run notebook cell makes -- clears the old template, a serial IES run's output
    in it included, and a staging copy a killed build left beside it."""

    model = _model(tmp_path)
    workspace = tmp_path / "tmpl"
    first = _calibration(model, workspace)
    _build(first)
    _assert_built_at(first, workspace)
    (workspace / "cal.phi.actual.csv").write_text("iteration,mean\n0,1.0\n", encoding="utf-8")
    staging = tmp_path / "tmpl.staging"
    staging.mkdir()
    _write_build_marker(staging, "staging")
    (staging / "mfsim.nam").write_text("half a copy\n", encoding="utf-8")

    again = _calibration(model, workspace)
    _build(again)

    _assert_built_at(again, workspace)
    assert not (workspace / "cal.phi.actual.csv").exists()
    assert not staging.exists()


def test_the_default_workspace_is_rebuilt(tmp_path):
    """The default, ``<model workspace>.pest/<name>``, built twice."""

    model = _model(tmp_path)
    for _ in range(2):
        cal = _calibration(model, None)
        _build(cal)
        _assert_built_at(cal, Path(model.workspace).parent / "r.pest" / "cal")


def test_an_empty_folder_is_used(tmp_path):
    model = _model(tmp_path)
    empty = tmp_path / "empty"
    empty.mkdir()
    cal = _calibration(model, empty)
    _build(cal)
    _assert_built_at(cal, empty)


def test_a_template_built_before_the_marker_is_rebuilt_only_as_this_calibrations(tmp_path):
    """A template built before ledger 360 has no marker. It is cleared when its
    ``myflopy_pest_metadata.json`` names THIS calibration and this run and it is named
    as ``model.pest(name)`` names one -- the template a re-run notebook cell meets
    after the upgrade. A copy of it under another calibration's name (an IES master
    is ``<name>_ies_master``) is not."""

    model = _model(tmp_path)
    default = Path(model.workspace).parent / "r.pest" / "cal"
    _build(_calibration(model, None))
    (default / BUILD_MARKER).unlink()
    master = default.parent / "cal_ies_master"
    shutil.copytree(default, master)
    before = _tree(master)

    # Another calibration's: its metadata names "cal", not "cal_ies_master".
    with pytest.raises(FileExistsError, match="no myflopy PEST build wrote it"):
        _build(_calibration(model, None, name="cal_ies_master"))
    # This calibration's metadata, but not where `model.pest("cal")` puts a template.
    with pytest.raises(FileExistsError, match="no myflopy PEST build wrote it"):
        _build(_calibration(model, master))
    assert _tree(master) == before

    cal = _calibration(model, None)
    _build(cal)
    _assert_built_at(cal, default)


def test_the_staging_copy_is_marked_as_the_builds_while_the_build_runs(tmp_path, monkeypatch):
    """A build killed while it runs leaves ``W.staging`` behind; the next build clears
    it only because it carries the build's marker naming it, written BEFORE the run is
    copied in -- and a marker the run itself holds (a run made in a copy of a
    template) is not copied over it."""

    model = _model(tmp_path)
    run = Path(model.workspace)
    (run / BUILD_MARKER).write_text(json.dumps(
        {"written_by": "myflopy PEST build", "role": "template", "directory": "/elsewhere"}),
        encoding="utf-8")
    seen = {}
    build_pstfrom = PestProject._build_pstfrom

    def _look(self, staging):
        seen.update(json.loads((Path(staging) / BUILD_MARKER).read_text(encoding="utf-8")))
        return build_pstfrom(self, staging)

    monkeypatch.setattr(PestProject, "_build_pstfrom", _look)
    cal = _calibration(model, tmp_path / "tmpl")
    _build(cal)

    assert seen == {"written_by": "myflopy PEST build", "role": "staging",
                    "directory": str((tmp_path / "tmpl.staging").resolve())}
    _assert_built_at(cal, tmp_path / "tmpl")
