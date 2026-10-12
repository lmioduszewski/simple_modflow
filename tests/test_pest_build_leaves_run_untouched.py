"""A PEST build leaves the source run and the model object untouched (ledger 356).

``PestProject.build()`` re-stored the CALLER's model -- each declared array target
made layered (``relayer_array_target``, ledger 110), every array and list moved into
an external file (``set_all_data_external``) -- and wrote the simulation back into
the RUN directory, which ``PstFrom`` then copied into the template. Measured on the
model below: every input file rewritten (7 with their arrays moved out, ``m.disv``
from 502 lines to 30) and 22 files added, 18 storage entries of the model object
changed; with a lake whose periods are empty, the build died inside FloPy having
already written 21 files into the run directory; on a run with inline arrays opened
with ``mf.load_mf6_run`` (loaded lazily) its array targets were refused; and pyEMU's
``PstFrom.log`` landed in the working directory.

It now copies the run's input files (the whole directory until ledger 383) to a
staging directory beside the template, loads the simulation from that copy -- keeping
LAK ``IMPLICIT``, which FloPy's loader drops (ledger 326) -- relayers and externalizes
THAT, and hands it to ``PstFrom``.

The model: a 12 x 12 irregular Voronoi grid of 20 ft cells, two layers (100 - 82 -
40 ft) with a pond excavated to 82 ft holding an ``IMPLICIT`` lake, under NEWTON;
scalar K 10 and K33 2 ft/d, so ``k`` and ``k33`` must be relayered (ledger 110);
CHD at 94 ft on the west and east edges of layer 2; recharge on every seventh active
cell; one steady and one two-step transient period, heads and budgets saved at every
step. Each test runs its own copy, so a build that writes into one cannot hide behind
another test's.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.lak_implicit import has_implicit, keep_implicit

pyemu = pytest.importorskip("pyemu")

#: Every file MF6 writes for this model, removed before a template is run again.
OUTPUTS = ("m.hds", "m.cbc", "m_lak_stage.lak", "m_lak_budget.lak", "m_lak_budget.csv",
           "m.lst", "mfsim.lst")
#: The output files compared byte for byte: MF6 writes no time stamp into them.
ANSWERS = ("m.hds", "m.cbc", "m_lak_stage.lak", "m_lak_budget.lak", "m_lak_budget.csv")
#: Cells the explicit pilot points sit on.
PP_CELLS = [3, 70, 140]


def _simulation(*, inflow):
    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs, geometry=[
        Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 32)])
    result = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(82.0))
              .add("lower", bottom=mf.Flat(40.0))
              .excavate(bottom=82.0, within=pond, name="pit").build(grid))
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    lake = mf.lak(context=context, nper=2, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=3000.0, implicit=True,
                  unexcavated="raise", inflow=inflow)
    recharge = [[(0, int(c)), 1e-3] for c in range(0, grid.ncpl, 7) if result.idomain[0][c] > 0]
    model = mf.gwf("m", context=context, packages=[
        mf.disv(result), mf.npf(k=10.0, k33=2.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 94.0] for c in edge]}),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, steady_state={0: True}, transient={1: True}),
        mf.rch(stress_period_data={0: recharge}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ], newtonoptions="NEWTON")
    return mf.simulation(model, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(1.0, 1, 1.0), (10.0, 2, 1.0)])])


def _run(root: Path, *, inflow=5.0):
    """The model in the module docstring, built and run; returns the ``Run``."""

    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and this run
        # is never reopened.
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", _simulation(inflow=inflow))
    assert run.success, "\n".join(run.report[-20:])
    return run


def _calibrate(model, workspace: Path | None, *, name="s1"):
    """A calibration of k (scalar: relayered), CHD (a list target) and k33 pilot points, built.

    ``workspace=None`` is ``model.pest``'s default, ``<model workspace>.pest/<name>``."""

    cal = model.pest(name, workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), capture=True)
    cal.parameterize("chd", style="constant", bounds=(-1.0, 1.0))
    grid = model.gwf.modelgrid
    points = np.c_[np.asarray(grid.xcellcenters)[PP_CELLS], np.asarray(grid.ycellcenters)[PP_CELLS]]
    cal.parameterize("k33", style="pilotpoints", pp_points=points, bounds=(0.5, 2.0))
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "cell": [30], "layer": [1]}),
        values=pd.DataFrame({"time": [0, 1], "W": [94.5, 94.6]}),
    ))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{name}.pst", noptmax=0)
    return cal


def _tree(directory) -> dict[str, str]:
    """Every file under ``directory``, by relative path, to the SHA-256 of its bytes."""

    directory = Path(directory)
    return {str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(directory.rglob("*")) if path.is_file()}


def _assert_same_tree(before: dict[str, str], after: dict[str, str]):
    changed = sorted(name for name in before.keys() & after.keys() if before[name] != after[name])
    added, removed = sorted(after.keys() - before.keys()), sorted(before.keys() - after.keys())
    assert not (changed or added or removed), (
        f"the build wrote into the run directory: changed {changed}, added {added}, "
        f"removed {removed}"
    )


def _storage(simulation) -> dict:
    """How every dataset of every package is STORED: layered or not, and each layer's form."""

    found = {}
    for model_name in simulation.model_names:
        for package in simulation.get_model(model_name).packagelist:
            for block in package.blocks.values():
                for name, data in block.datasets.items():
                    held = getattr(data, "_data_storage", None)
                    for key, storage in (held.items() if isinstance(held, dict) else [(None, held)]):
                        if storage is None:
                            continue
                        if hasattr(storage, "layer_storage"):
                            form = (bool(getattr(storage, "layered", False)), tuple(
                                str(layer.data_storage_type).rsplit(".", 1)[-1]
                                for layer in storage.layer_storage.elements()))
                        else:
                            form = (str(storage.data_storage_type).rsplit(".", 1)[-1],)
                        found[(model_name, package.filename, name, key)] = form
    return found


def _mf6() -> str:
    exe = shutil.which("mf6")
    if exe is None:
        pytest.skip("mf6 not on PATH")
    return exe


def _run_again(template: Path, destination: Path) -> Path:
    """MF6 run directly in a copy of a template -- the model PEST is handed, no multiplier applied."""

    shutil.copytree(template, destination)
    for name in OUTPUTS:
        (destination / name).unlink(missing_ok=True)
    result = subprocess.run([_mf6()], cwd=destination, capture_output=True, text=True)
    assert "Normal termination" in result.stdout, result.stdout[-2000:]
    return destination


def _assert_same_answers(source: Path, rerun: Path):
    """Heads, the cell budget, the lake's stage and its budget: the same bytes."""

    for name in ANSWERS:
        assert (rerun / name).read_bytes() == (Path(source) / name).read_bytes(), (
            f"{name} differs between the source run and its template's model"
        )


# -- (a) the run directory, (b) the model object ------------------------------------------


def test_a_build_leaves_every_file_of_the_run_directory_as_it_was(tmp_path):
    """Hashed before and after: nothing changed, added or removed. The build used to
    rewrite 7 of this run's files with every array external (``m.disv`` 502 lines to
    30) and add 22."""

    run = _run(tmp_path)
    before = _tree(run.workspace)
    _calibrate(run.model(), tmp_path / "pest")
    _assert_same_tree(before, _tree(run.workspace))


def test_a_build_leaves_the_model_object_as_it_was(tmp_path):
    """Every dataset keeps its storage -- internal stays internal, a scalar stays a
    scalar, nothing is made layered -- and the simulation keeps its path. The build
    used to relayer ``k`` and ``k33`` and externalize all 18 stored entries."""

    run = _run(tmp_path)
    model = run.model()
    storage, path = _storage(model.sim), Path(model.sim.sim_path)
    k = ("m", "m.npf", "k", None)
    assert storage[k] == (False, ("internal_constant",)), "the fixture must hold a scalar k"

    _calibrate(model, tmp_path / "pest")

    after = _storage(model.sim)
    assert after[k] == (False, ("internal_constant",))
    assert after == storage, {key: (storage.get(key), after.get(key))
                              for key in storage.keys() | after.keys()
                              if storage.get(key) != after.get(key)}
    assert Path(model.sim.sim_path) == path


# -- (c)/(d) the staged copy is the run, and answers as it did ------------------------------


@pytest.mark.parametrize("inflow", [5.0, None], ids=["inflow", "empty_lake_periods"])
def test_the_template_gives_the_source_runs_answer_and_keeps_implicit(tmp_path, inflow):
    """MF6 run in the template gives the source run's heads, cell budget, lake stage
    and lake budget, byte for byte, with ``IMPLICIT`` in its LAK file and the scalar
    ``k`` stored one file per layer.

    ``empty_lake_periods``: an ``mf.lak`` with no per-period setting writes empty
    period blocks, and the old build died on the in-memory model inside FloPy's
    ``set_all_data_external`` (``MFDataException`` on LAK ``perioddata``) after it
    had written external files into the run directory; the copy LOADED from the
    files externalizes."""

    run = _run(tmp_path, inflow=inflow)
    cal = _calibrate(run.model(), tmp_path / "pest")
    template = cal.template_workspace

    assert has_implicit(template / "m.lak"), "the staged load dropped LAK IMPLICIT"
    assert sorted(path.name for path in template.glob("m.npf_k_layer*.txt")) == [
        "m.npf_k_layer1.txt", "m.npf_k_layer2.txt"]
    _assert_same_answers(run.workspace, _run_again(template, tmp_path / "rerun"))


def test_a_run_opened_with_load_mf6_run_builds_its_array_targets(tmp_path):
    """``mf.load_mf6_run`` loads a run LAZILY, so on the caller's object the old build
    found no NPF package to relayer or externalize, and refused an array target of a
    run whose arrays are inline: ``FileNotFoundError: No external input file found for
    target 'k'`` -- on this model and on a real 411 MB run. The staged copy is loaded
    whole."""

    run = _run(tmp_path)
    model = mf.load_mf6_run(Path(run.workspace), verbosity_level=0)
    before = _tree(run.workspace)
    cal = _calibrate(model, tmp_path / "pest")

    _assert_same_tree(before, _tree(run.workspace))
    template = cal.template_workspace
    assert has_implicit(template / "m.lak")
    assert sorted(path.name for path in template.glob("m.npf_k_layer*.txt")) == [
        "m.npf_k_layer1.txt", "m.npf_k_layer2.txt"]
    _assert_same_answers(run.workspace, _run_again(template, tmp_path / "rerun"))


def test_a_loaded_run_whose_arrays_are_already_external(tmp_path):
    """A run whose arrays are in external files (whole-grid ``k``, not layered), opened
    with ``mf.load_mf6_run``: its directory is left as it was, and its template answers
    as it did."""

    run = _run(tmp_path)
    external = tmp_path / "external"
    shutil.copytree(run.workspace, external)
    simulation = flopy.mf6.MFSimulation.load(sim_ws=str(external), verbosity_level=0)
    keep_implicit(simulation)
    simulation.set_all_data_external()
    simulation.write_simulation(silent=True)
    for name in OUTPUTS:
        (external / name).unlink(missing_ok=True)
    result = subprocess.run([_mf6()], cwd=external, capture_output=True, text=True)
    assert "Normal termination" in result.stdout, result.stdout[-2000:]
    assert (external / "m.npf_k.txt").is_file() and not list(external.glob("m.npf_k_layer*"))

    model = mf.load_mf6_run(external, verbosity_level=0)
    before = _tree(external)
    cal = _calibrate(model, tmp_path / "pest")

    _assert_same_tree(before, _tree(external))
    template = cal.template_workspace
    assert has_implicit(template / "m.lak")
    assert sorted(path.name for path in template.glob("m.npf_k_layer*.txt")) == [
        "m.npf_k_layer1.txt", "m.npf_k_layer2.txt"]
    _assert_same_answers(external, _run_again(template, tmp_path / "rerun"))


def test_the_template_is_the_runs_files_not_an_edit_made_in_memory_after_the_run(tmp_path):
    """The run's FILES win, because its output -- the base observations -- was
    computed from them: ``k`` and ``k33`` edited on the model object after the run
    reach neither the template's ``k`` files nor the pilot points' base array, and the
    edits stay on the object, untouched."""

    run = _run(tmp_path)
    model = run.model()
    model.gwf.npf.k.set_data(99.0)
    model.gwf.npf.k33.set_data(7.0)

    cal = _calibrate(model, tmp_path / "pest")
    template = cal.template_workspace

    for layer in (1, 2):
        assert np.all(np.loadtxt(template / f"m.npf_k_layer{layer}.txt") == 10.0)
    for layer in (0, 1):
        cells = pd.read_csv(template / f"k33l{layer}_pp.cells.csv")
        assert np.all(cells["base_value"] == 2.0), "pilot points read the object's k33"
    assert np.all(np.asarray(model.gwf.npf.k.get_data()) == 99.0)
    assert np.all(np.asarray(model.gwf.npf.k33.get_data()) == 7.0)


# -- (e) what the build reads is the staged copy, and parameters still apply ---------------


def test_pilot_points_and_the_spatial_reference_read_the_staged_grid(tmp_path):
    """The cell centres pilot points interpolate between, and pyEMU's spatial
    reference, are the template's own DISV as written -- FloPy writes 8 decimals --
    not the object's full-precision grid (which they read until the object was left
    alone; before that, the externalization had rebuilt the object's grid from the
    files)."""

    run = _run(tmp_path)
    model = run.model()
    cal = _calibrate(model, tmp_path / "pest")
    written = flopy.mf6.MFSimulation.load(
        sim_ws=str(cal.template_workspace), verbosity_level=0, load_only=["disv"],
    ).get_model("m").modelgrid
    cells = pd.read_csv(cal.template_workspace / "k33l0_pp.cells.csv")

    assert np.array_equal(cells["x"].to_numpy(), np.asarray(written.xcellcenters, dtype=float))
    assert np.array_equal(cells["y"].to_numpy(), np.asarray(written.ycellcenters, dtype=float))
    assert not np.array_equal(cells["x"].to_numpy(),
                              np.asarray(model.gwf.modelgrid.xcellcenters, dtype=float)), (
        "the fixture's in-memory grid must differ from the written one, or this checks nothing")
    assert np.array_equal(np.asarray(cal.pf._spatial_reference.xcellcenters, dtype=float),
                          np.asarray(written.xcellcenters, dtype=float))


def test_a_multiplier_still_reaches_the_model_in_a_forward_run(tmp_path):
    """A realization with ``k`` doubled: pyEMU's forward run writes 2x the base ``k``
    into both layer files and the heads move."""

    run = _run(tmp_path)
    cal = _calibrate(run.model(), tmp_path / "pest")
    template = cal.template_workspace

    def forward():
        result = subprocess.run([sys.executable, "forward_run.py"], cwd=template,
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]
        return pd.read_csv(template / "hds_simulated_heads.csv")

    base_heads = forward()
    base_k = [np.loadtxt(template / f"m.npf_k_layer{layer}.txt") for layer in (1, 2)]
    parameters = cal.pst.parameter_data
    parameters.loc[parameters["pargp"] == "k", "parval1"] = 2.0
    for tpl, written in zip(cal.pst.template_files, cal.pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(parameters["parval1"], str(template / tpl),
                                          str(template / written))
    moved_heads = forward()

    for layer, base in zip((1, 2), base_k, strict=True):
        assert np.allclose(np.loadtxt(template / f"m.npf_k_layer{layer}.txt"), 2.0 * base)
    assert not np.allclose(moved_heads["W"], base_heads["W"]), "k doubled moved no head"


# -- (f) nothing in the working directory; the staging copy goes ---------------------------


def test_nothing_is_written_into_the_working_directory(tmp_path, monkeypatch):
    """pyEMU's ``PstFrom`` opens ``PstFrom.log`` in the working directory -- a test run
    left one in the repository root. It is in the template, and still being written
    after the build moved it there."""

    run = _run(tmp_path)
    model = run.model()
    working = tmp_path / "cwd"
    working.mkdir()
    monkeypatch.chdir(working)

    cal = _calibrate(model, tmp_path / "pest")

    assert sorted(path.name for path in working.iterdir()) == []
    log = (cal.template_workspace / "PstFrom.log").read_text()
    assert "starting PstFrom process" in log, "the log's start was lost in the move"
    assert "finished: setting up dirs" in log, "the log's copy step was lost in the move"
    assert "style parameters for file(s)" in log, "nothing logged after the move"
    assert not cal._staging_workspace().exists()


def test_a_relative_model_workspace_builds_its_template_where_it_names(tmp_path, monkeypatch):
    """Every canonical PEST notebook runs its model under ``Path('../artifacts/...')``,
    so the model's workspace -- and the default template beside it -- are RELATIVE.
    ``PstFrom`` is constructed with the staging directory as the working directory, and
    the template path was resolved INSIDE that: it was made in the staging copy and
    deleted with it, and the build died moving the log (``FileNotFoundError``) -- ledger
    357, the review of 356. The template is where its relative path names, from the
    caller's working directory, which holds nothing else."""

    monkeypatch.chdir(tmp_path)
    run = _run(Path("rel"))
    model = run.model()
    assert not Path(model.workspace).is_absolute(), "the fixture must give a relative workspace"
    before = _tree(run.workspace)

    cal = _calibrate(model, None)
    template = cal.template_workspace

    assert not template.is_absolute()
    assert template.resolve() == (tmp_path / "rel" / "p" / "runs" / "r.pest" / "s1").resolve()
    assert (template / "s1.pst").is_file() and has_implicit(template / "m.lak")
    assert "style parameters for file(s)" in (template / "PstFrom.log").read_text()
    assert not cal._staging_workspace().exists()
    assert sorted(path.name for path in tmp_path.iterdir()) == ["rel"]
    _assert_same_tree(before, _tree(run.workspace))
    _assert_same_answers(run.workspace, _run_again(template, tmp_path / "rerun"))


def test_a_second_build_closes_the_first_builds_log_before_deleting_its_template(tmp_path):
    """The log has lived IN the template since ledger 356, and a rebuild of one
    project deletes the template first, while the first ``PstFrom`` still held its log
    open -- which Windows refuses to delete (reasoned from Windows' file semantics, not
    run here). The rebuild closes it first (ledger 357), and logs to its own."""

    run = _run(tmp_path)
    cal = _calibrate(run.model(), tmp_path / "pest")
    first = cal.pf.logger.f
    assert not first.closed

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("s1.pst", noptmax=0)

    assert first.closed, "the first build's PstFrom.log is still open in the template it deletes"
    assert not cal.pf.logger.f.closed
    assert "style parameters for file(s)" in (cal.template_workspace / "PstFrom.log").read_text()


def test_a_model_never_written_is_refused_and_its_directory_left_alone(tmp_path):
    """A run built in memory has no MF6 files: the build refuses, naming the spelling,
    where it used to write the simulation into the run directory."""

    built = _simulation(inflow=5.0).build(tmp_path / "unwritten")
    before = _tree(built.workspace)
    cal = built.model().pest("s1", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant")
    with pytest.raises(FileNotFoundError, match=r"no mfsim\.nam.*Write and run the model first"):
        cal.build()
    _assert_same_tree(before, _tree(built.workspace))


def test_a_failed_build_leaves_no_staging_copy_and_the_run_as_it_was(tmp_path):
    """A target at a period the run does not have is refused after the copy is staged;
    the staging directory goes with the failure, and the run is untouched."""

    run = _run(tmp_path)
    before = _tree(run.workspace)
    cal = run.model().pest("s1", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant")
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "cell": [30], "layer": [1]}),
        values=pd.DataFrame({"time": [9], "W": [94.5]}),
    ))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="W"):
            cal.build()
    assert not cal._staging_workspace().exists()
    _assert_same_tree(before, _tree(run.workspace))
