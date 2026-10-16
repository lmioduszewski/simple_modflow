"""A run's own packages, simulation and mesh go into a project's library (ledger 613).

``run.to_library(project)`` registers the run's RECORD -- its ``_library`` copies under the
keys it recorded, its simulation from the frozen snapshot with every library slot back as
``mf.ref(key)`` and every inline package as declared, and each model's mesh in
``project.grids`` behind ``mf.grid_ref(key)`` -- after building and writing that import in
a scratch directory and comparing it, file by file, with the files the run wrote.

The model is a one-layer DISV on the session's TRIANGULATED grid (ledgers 263/265: the mesh
is pinned, frozen, pickled and saved here), with a lazy discretization and a lazy
``mf.chd.gpkg`` among its library packages. Nothing here runs MF6: every check is on the
written input files.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import flopy
import geopandas as gpd
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.grid_identity import grid_pin
from myflopy.modflow.mf6.pest.calibrated import CalibratedInputs
from myflopy.run_record import written_differences
from myflopy.specs import GridRef, PackageRef, PackageSpec, package_digest
from myflopy.workspace import LibraryReport, Run

KEYS = ("dis/base", "npf/base", "ic/base", "chd/west", "oc/base")


def _write_source(path: Path, grid, head: float) -> None:
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [head]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)], crs=grid.crs,
    ).to_file(path, layer="chd", driver="GPKG")


def _simulation(grid, *, packages=None, tdis=None) -> mf.SimulationSpec:
    return mf.simulation(
        mf.gwf("m", context=mf.ModelContext(grid=grid),
               packages=packages or [mf.ref(key) for key in KEYS]),
        name="s",
        packages=[tdis or mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    )


def _library(project, grid, source: Path) -> None:
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    project.add_package("dis/base", mf.disv(stack))
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_package("ic/base", mf.ic(strt=95.0))
    project.add_package("chd/west", mf.chd.gpkg(source, layer="chd", nper=1, head="head"))
    project.add_package("oc/base", mf.oc(head_filerecord="m.hds",
                                         saverecord=[("HEAD", "ALL")]))


@pytest.fixture
def made(tmp_path, triangulated_grid):
    """A project whose run ``r`` of simulation ``s`` is written: a namespace."""

    grid = triangulated_grid
    source = tmp_path / "bc.gpkg"
    _write_source(source, grid, 99.0)
    root = tmp_path / "proj"
    project = mf.Project(root, name="p")
    _library(project, grid, source)
    simulation = project.add_simulation(_simulation(grid))
    run = project.prepare_run("r", simulation)
    run.write(silent=True)
    return SimpleNamespace(project=project, run=run, root=root, grid=grid, source=source)


def _fresh(made, name: str = "fresh") -> mf.Project:
    """The same runs, an empty library."""

    return mf.Project(made.root, name=name)


def _labels(packages) -> list[str]:
    return [entry.key if isinstance(entry, PackageRef) else entry.name for entry in packages]


def _edit(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    assert old in text, (path, old)
    path.write_text(text.replace(old, new, 1))


# -- what goes in --------------------------------------------------------------------


def test_a_run_comes_back_into_an_empty_library_under_the_keys_it_recorded(made):
    fresh = _fresh(made)
    run = fresh.runs["r"]
    report = run.to_library(fresh)

    assert isinstance(report, LibraryReport) and mf.LibraryReport is LibraryReport
    assert report.library_keys == {key: key for key in sorted(KEYS)}
    assert {key: package_digest(fresh.packages[key]) for key in KEYS} == run.package_digests
    assert report.simulation == "s" and fresh.simulations.names() == ["s"]
    model = fresh.simulations["s"].models[0]
    assert _labels(model.packages) == list(run.package_keys["m"])
    assert report.inputs.verdict == "ok" and report.differences == ()
    assert report.unsaveable is None
    # The mesh: registered once, unprefixed, under the model and its pin, and the model
    # declared on it rather than carrying it in its context.
    pin = grid_pin(made.grid)
    assert report.grids == {"m": f"m/{pin}"}
    assert model.grid == GridRef(f"m/{pin}") and model.context.grid is None
    assert grid_pin(fresh.grids[f"m/{pin}"].obj) == pin
    assert sorted(report.added) == sorted(
        [("package", key) for key in KEYS] + [("grid", f"m/{pin}"), ("simulation", "s")])
    assert report.identical == report.replaced == ()
    assert "inputs='ok'" in repr(report) and "<table>" in report._repr_html_()


def test_what_it_imported_writes_the_files_the_run_ran(made):
    fresh = _fresh(made)
    fresh.runs["r"].to_library(fresh)
    again = fresh.prepare_run("again", "s")
    again.write(silent=True)
    assert written_differences(made.run.workspace, again.workspace) == []


@pytest.mark.parametrize("door", ["load_run", "reopen_run", "runs"])
def test_every_door_hands_the_same_verb(made, door):
    fresh = _fresh(made)
    run = {
        "load_run": lambda: mf.load_run(made.run.workspace, load=False),
        "reopen_run": lambda: made.project.reopen_run("r", load=False),
        "runs": lambda: fresh.runs["r"],
    }[door]()
    assert isinstance(run, Run)
    report = run.to_library(fresh)
    assert report.library_keys == {key: key for key in sorted(KEYS)}
    assert report.inputs.verdict == "ok"


def test_importing_twice_changes_nothing(made):
    fresh = _fresh(made)
    fresh.runs["r"].to_library(fresh)
    held = fresh.simulations["s"]
    packages = {key: fresh.packages[key] for key in KEYS}
    report = fresh.runs["r"].to_library(fresh)
    assert report.added == () and report.replaced == ()
    assert len(report.identical) == len(KEYS) + 2
    assert fresh.simulations["s"] is held
    assert all(fresh.packages[key] is packages[key] for key in KEYS)


def test_the_project_that_prepared_it_already_holds_it(made):
    """Its own library: every package and the simulation are the run's -- identical."""

    report = made.project.runs["r"].to_library(made.project, overwrite=False)
    assert {key for kind, key in report.identical if kind == "package"} == set(KEYS)
    assert ("simulation", "s") in report.identical
    assert report.added == (("grid", report.grids["m"]),)


# -- collisions, prefix, overwrite -----------------------------------------------------


def test_a_key_holding_something_else_is_refused_and_nothing_is_registered(made):
    fresh = _fresh(made)
    fresh.add_package("npf/base", mf.npf(k=5.0))
    with pytest.raises(KeyError, match=r"package 'npf/base'.*prefix=.*overwrite=True.*"
                                       r"WITHOUT a prefix"):
        fresh.runs["r"].to_library(fresh)
    assert fresh.packages.names() == ["npf/base"]
    assert fresh.simulations.names() == [] and fresh.grids.names() == []
    # A simulation of the same name building something else names how it differs.
    fresh.add_simulation(_simulation(made.grid, tdis=mf.tdis(nper=1, perioddata=[(9.0, 1, 1.0)])))
    with pytest.raises(KeyError, match=r"simulation 's' \(simulation package 'tdis' changed"):
        fresh.runs["r"].to_library(fresh)


def test_prefix_keeps_both_and_renames_the_refs(made):
    fresh = _fresh(made)
    fresh.add_package("npf/base", mf.npf(k=5.0))
    report = fresh.runs["r"].to_library(fresh, prefix="old/")
    assert report.library_keys["npf/base"] == "old/npf/base"
    assert fresh.package("npf/base").options["k"] == 5.0  # untouched
    imported = fresh.simulations["old/s"]
    assert report.simulation == "old/s"
    assert _labels(imported.models[0].packages) == [f"old/{key}" for key in KEYS]
    assert imported.derived_from is None
    assert imported.lineage[-1]["operation"] == "to_library"
    assert imported.lineage[-1]["prefix"] == "old/"
    # The mesh is shared, never prefixed.
    assert report.grids["m"] == f"m/{grid_pin(made.grid)}"
    again = fresh.prepare_run("again", "old/s")
    again.write(silent=True)
    assert written_differences(made.run.workspace, again.workspace) == []


def test_overwrite_replaces_what_differs_and_reports_it(made):
    fresh = _fresh(made)
    fresh.add_package("npf/base", mf.npf(k=5.0))
    report = fresh.runs["r"].to_library(fresh, overwrite=True)
    assert report.replaced == (("package", "npf/base"),)
    assert package_digest(fresh.packages["npf/base"]) == made.run.package_digests["npf/base"]


def test_a_grid_key_holding_another_mesh_is_only_replaced(made, triangulated_grid_finer):
    fresh = _fresh(made)
    key = f"m/{grid_pin(made.grid)}"
    fresh.add_grid(key, triangulated_grid_finer)
    with pytest.raises(KeyError, match=r"grid 'm/.*never prefixed.*overwrite=True"):
        fresh.runs["r"].to_library(fresh, prefix="x/")
    assert fresh.simulations.names() == []
    report = fresh.runs["r"].to_library(fresh, overwrite=True)
    assert ("grid", key) in report.replaced
    assert grid_pin(fresh.grids[key].obj) == grid_pin(made.grid)


# -- the files against the record ------------------------------------------------------


def _tdis_changed_outside_myflopy(run):
    _edit(run.workspace / "s.tdis", "1.00000000  1       1.00000000",
          "2.00000000  3       1.00000000")


def _tdis_changed_through_flopy(run):
    run.simulation.tdis.perioddata.set_data([(2.0, 3, 1.0)])
    with pytest.warns(UserWarning, match="edited through FloPy"):
        run.write(silent=True)


@pytest.mark.parametrize("change", [_tdis_changed_outside_myflopy, _tdis_changed_through_flopy])
def test_a_file_changed_after_the_run_wrote_it_is_refused_naming_keep_record(made, change):
    change(made.run)
    fresh = _fresh(made)
    with pytest.raises(ValueError, match=r"(?s)s\.tdis PERIODDATA: 2 of 3 values differ.*"
                                         r'files_differ="keep_record"'):
        fresh.runs["r"].to_library(fresh)
    assert fresh.packages.names() == [] and fresh.simulations.names() == []


def test_keep_record_imports_the_record_warns_and_lists_the_differences(made):
    _tdis_changed_outside_myflopy(made.run)
    fresh = _fresh(made)
    with pytest.warns(UserWarning, match=r"imported its RECORD.*s\.tdis"):
        report = fresh.runs["r"].to_library(fresh, files_differ="keep_record")
    assert report.inputs.verdict == "edited" and list(report.inputs.departs) == ["s.tdis"]
    assert [(file, block) for file, block, _ in report.differences] == [
        ("s.tdis", "PERIODDATA")]
    tdis = fresh.simulations["s"].packages[0]
    assert [tuple(row) for row in tdis.options["perioddata"]] == [(1.0, 1, 1.0)]
    # The correction the docs show, right after: the library then says what MF6 ran.
    fresh.add_simulation(fresh.simulations["s"].with_packages(
        mf.tdis(nper=1, perioddata=[(2.0, 3, 1.0)])), overwrite=True)
    again = fresh.prepare_run("again", "s")
    again.write(silent=True)
    assert written_differences(made.run.workspace, again.workspace) == []


def test_an_array_edited_through_flopy_is_seen(made):
    """The case ``ModelDiff`` calls identical: it compares no array package."""

    made.run.simulation.get_model("m").ic.strt.set_data(50.0)
    with pytest.warns(UserWarning, match="edited through FloPy"):
        made.run.write(silent=True)
    fresh = _fresh(made)
    with pytest.raises(ValueError, match=r"m\.ic GRIDDATA"):
        fresh.runs["r"].to_library(fresh)


def test_a_lazy_library_package_comes_back_as_its_recipe_in_its_own_slot(made):
    fresh = _fresh(made)
    fresh.runs["r"].to_library(fresh)
    slot = fresh.simulations["s"].models[0].packages[KEYS.index("chd/west")]
    assert slot == PackageRef("chd/west")
    chd = fresh.packages["chd/west"]
    assert package_digest(chd) == package_digest(made.project.packages["chd/west"])
    # Its source edited since the run: the library would resolve it again, so the
    # record no longer writes the run's files -- here, unlike the write-backs' check
    # against what the run froze (ledger 612).
    _write_source(made.source, made.grid, 42.0)
    other = _fresh(made, "other")
    with pytest.raises(ValueError, match=r"m\.chd PERIOD 1"):
        other.runs["r"].to_library(other)


def test_an_inline_lazy_package_comes_back_as_its_declaration(made):
    """The snapshot holds ``mf.disv(stack)`` RESOLVED -- arrays no project file holds."""

    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    declared = mf.disv(stack)
    made.project.add_simulation(mf.simulation(
        mf.gwf("m", context=mf.ModelContext(grid=made.grid),
               packages=[declared, *map(mf.ref, KEYS[1:])]),
        name="s2", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    ))
    run = made.project.prepare_run("r2", "s2")
    run.write(silent=True)
    fresh = _fresh(made)
    report = fresh.runs["r2"].to_library(fresh)
    disv = fresh.simulations["s2"].models[0].packages[0]
    assert isinstance(disv, PackageSpec)
    assert package_digest(disv) == package_digest(declared)
    assert report.unsaveable is None and report.inputs.verdict == "ok"
    assert fresh.save().simulations == 1


# -- refusals --------------------------------------------------------------------------


def test_a_bare_mf6_workspace_is_refused_naming_what_is_not_built(tmp_path):
    workspace = tmp_path / "bare"
    sim = flopy.mf6.MFSimulation(sim_ws=str(workspace))
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="flow")
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=2, ncol=2, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    sim.write_simulation(silent=True)
    project = mf.Project(tmp_path / "proj", name="p")
    with pytest.raises(ValueError, match=r"no run\.json.*Phase-2 decision D2.*not built"):
        mf.load_run(workspace, load=False).to_library(project)


def test_a_run_prepared_before_runs_kept_copies_is_refused(made):
    manifest = made.run.manifest_path
    data = json.loads(manifest.read_text())
    data.pop("package_snapshot")
    manifest.write_text(json.dumps(data))
    fresh = _fresh(made)
    with pytest.raises(ValueError, match=r"holds no copy of them.*ledger 230"):
        fresh.runs["r"].to_library(fresh)


def test_a_copy_edited_on_disk_is_refused(made):
    copy = made.run.workspace / made.run.package_snapshot["npf/base"]
    data = json.loads(copy.read_text())
    data["options"]["k"] = 7.0
    copy.write_text(json.dumps(data))
    fresh = _fresh(made)
    with pytest.raises(ValueError, match=r"copy of library package 'npf/base'.*no longer "
                                         r"digests"):
        fresh.runs["r"].to_library(fresh)


def test_a_run_never_written_imports_and_says_nothing_was_compared(made):
    run = made.project.prepare_run("r_unwritten", "s")
    fresh = _fresh(made)
    report = fresh.runs[run.name].to_library(fresh)
    assert report.inputs.verdict == "unwritten" and report.differences == ()
    assert fresh.simulations.names() == ["s"]


def test_files_differ_and_prefix_are_checked_at_the_call(made):
    fresh = _fresh(made)
    run = fresh.runs["r"]
    with pytest.raises(ValueError, match=r'files_differ= is "raise".*"keep_record"'):
        run.to_library(fresh, files_differ="keep_files")
    with pytest.raises(ValueError, match="non-empty str"):
        run.to_library(fresh, prefix="")
    with pytest.raises(TypeError, match="non-empty str"):
        run.to_library(fresh, prefix=3)
    assert fresh.packages.names() == []


def test_a_relative_source_is_refused_across_project_roots(made, tmp_path):
    made.project.save(bundle_sources=True)
    run = made.project.prepare_run("r_bundled", "s")
    run.write(silent=True)
    elsewhere = mf.Project(tmp_path / "elsewhere", name="e")
    with pytest.raises(ValueError, match=r"relative to its project root.*name other files"):
        mf.load_run(run.workspace, load=False).to_library(elsewhere)
    assert elsewhere.packages.names() == []
    report = _fresh(made).runs["r_bundled"].to_library(_fresh(made, "same_root"))
    assert report.inputs.verdict == "ok"


# -- the mesh through a save -----------------------------------------------------------


def test_the_runs_mesh_survives_a_save(made):
    fresh = _fresh(made)
    report = fresh.runs["r"].to_library(fresh)
    saved = fresh.save()
    assert saved.simulations == 1 and saved.grids == 1 and saved.context_dropped == ()
    later = mf.Project.load(made.root)
    key = report.grids["m"]
    assert later.simulations["s"].models[0].grid == GridRef(key)
    assert grid_pin(later.grids[key].resolve(project_root=made.root, build=False)) == (
        grid_pin(made.grid))
    again = later.prepare_run("again", "s")
    again.write(silent=True)
    assert written_differences(made.run.workspace, again.workspace) == []
    # Importing again in the later session finds every entry already there.
    assert later.runs["r"].to_library(later).added == ()


def _imported_with_k_array(made, *, disv=None) -> mf.Project:
    """A run whose K is one value per cell, imported into a fresh library, saved, reloaded."""

    import numpy as np

    made.project.add_package("npf/base", mf.npf(k=np.linspace(1.0, 9.0, made.grid.ncpl)),
                             overwrite=True)
    if disv is not None:
        made.project.add_package("dis/base", disv, overwrite=True)
    made.project.prepare_run("rk", "s").write(silent=True)
    fresh = _fresh(made)
    fresh.runs["rk"].to_library(fresh)
    fresh.save()
    return mf.Project.load(made.root)


def test_a_saved_import_still_has_the_mesh_its_arrays_are_on(made, triangulated_grid_finer):
    """Through ``grid_ref`` the reloaded model still knows its mesh, so ``remap=True``
    carries the per-cell K; a context grid would have reloaded as None ("holds no grid")."""

    later = _imported_with_k_array(made)
    with pytest.warns(UserWarning, match="remapped k onto the new mesh"):
        moved = later.simulations["s"].with_grid(
            triangulated_grid_finer, remap=True, library=later)
    npf = moved.models[0].packages[KEYS.index("npf/base")]
    assert isinstance(npf, PackageSpec)
    assert len(npf.options["k"]) == triangulated_grid_finer.ncpl


@pytest.mark.xfail(strict=True, reason="needs the REMAP lane")
def test_a_saved_import_with_a_built_discretization_moves(made, triangulated_grid_finer):
    """An eager ``mf.disv(result)`` holds top/botm/idomain per cell; ``grid_free`` says it
    is re-declared, which the remap at this commit does not read (RemapRefused)."""

    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    later = _imported_with_k_array(made, disv=mf.disv(stack.build(made.grid)))
    with pytest.warns(UserWarning, match="remapped k onto the new mesh"):
        later.simulations["s"].with_grid(
            triangulated_grid_finer, remap=True, library=later, grid_free=("dis/base",))


# -- the write-back's comparison is the same one ---------------------------------------


def _refusal(project, run, source="s", swap=None):
    """``CalibratedInputs._refuse_a_moved_simulation`` on ``run``, without a calibration."""

    stub = SimpleNamespace(_run=lambda: Run.reopen(run.workspace, load=False))
    CalibratedInputs._refuse_a_moved_simulation(stub, project, run, source, swap or {})


def test_the_write_back_compares_a_lazy_library_package_by_its_key(made):
    """The snapshot holds a lazy ref's slot RESOLVED: every such run was refused as
    "inline package 'chd' changed" against the library that prepared it."""

    _refusal(made.project, made.run)
    moved = made.project.simulations["s"].with_packages(mf.ic(strt=1.0), slot="ic/base", on="m")
    made.project.simulations["s"] = moved
    with pytest.raises(ValueError, match=r"model 'm' carries .* where the run carried"):
        _refusal(made.project, made.run)


def test_the_write_back_finds_an_imported_run_identical(made):
    fresh = _fresh(made)
    fresh.runs["r"].to_library(fresh)
    _refusal(fresh, fresh.runs["r"])
    fresh.add_simulation(fresh.simulations["s"].with_packages(
        mf.tdis(nper=1, perioddata=[(5.0, 1, 1.0)])), overwrite=True)
    with pytest.raises(ValueError, match=r"simulation package 'tdis' changed"):
        _refusal(fresh, fresh.runs["r"])


def test_to_library_never_writes_into_the_run(made):
    """It reads the run; the check is built and written in a scratch directory."""

    before = {path: path.stat().st_mtime_ns for path in made.run.workspace.rglob("*")}
    fresh = _fresh(made)
    fresh.runs["r"].to_library(fresh)
    assert {path: path.stat().st_mtime_ns for path in made.run.workspace.rglob("*")} == before


def test_a_run_on_a_deferred_grid_ref_is_not_refused_by_its_mesh(tmp_path,
                                                                  triangulated_grid_finer):
    """A model declared on ``mf.grid_ref(key)`` to a RECIPE: the library names no mesh
    until it is built, the snapshot names the one the build resolved -- not a difference
    (ledger 614). It read ``mesh ref:vor where the run had pin:...``."""

    root = tmp_path / "proj"
    root.mkdir()
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"name": ["d"]}, geometry=[box(0.0, 0.0, 600.0, 600.0)],
                     crs=2927).to_file(boundary)
    project = mf.Project(root, name="p")
    project.add_grid("vor", mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=40_000.0))
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    project.add_package("dis/base", mf.disv(stack))
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_package("ic/base", mf.ic(strt=95.0))
    project.add_simulation(mf.simulation(
        mf.gwf("m", grid=mf.grid_ref("vor"),
               packages=[mf.ref(key) for key in ("dis/base", "npf/base", "ic/base")]),
        name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])],
    ))
    run = project.prepare_run("r", "s")
    run.write(silent=True)
    _refusal(project, run)
    # A fresh library takes the mesh under the run's own grid_ref key.
    fresh = _fresh(SimpleNamespace(root=root))
    report = fresh.runs["r"].to_library(fresh)
    assert report.grids == {"m": "vor"} and report.inputs.verdict == "ok"
    _refusal(fresh, fresh.runs["r"])
    # A ref to another BUILT mesh still differs, by pin.
    fresh.add_grid("finer", triangulated_grid_finer)
    model = fresh.simulations["s"].models[0]
    fresh.simulations["s"] = replace(fresh.simulations["s"], models=(
        replace(model, grid=GridRef("finer")),))
    with pytest.raises(ValueError, match=r"mesh pin:\w+ where the run had pin:"):
        _refusal(fresh, fresh.runs["r"])
