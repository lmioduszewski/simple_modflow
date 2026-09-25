"""``mf.lak(implicit=True)``: MF6's LAK ``IMPLICIT`` option, which FloPy 3.10 cannot write (ledger 326).

The stage of each lake becomes one more unknown of the flow matrix, solved WITH the
heads instead of by turns. FloPy's ``ModflowGwflak`` has no keyword for it, so a second
writer inserts it into the written file -- ``lak_implicit.ModflowGwflakImplicit``, the
package's own ``write``, which every write path reaches. What these tests pin:

- the flag rides in the spec's OPTIONS (serialized, frozen, digested; ``False`` adds no
  key, so no existing digest moves) and FloPy never sees it;
- the keyword lands inside ``BEGIN/END OPTIONS``, exactly once, from every write path:
  ``Run.write``/``execute``, a direct ``write_simulation`` on a built simulation, a
  run reopened from its frozen snapshot, a run LOADED by ``Run.open`` or
  ``mf.load_mf6_run`` -- and NOT from a bare ``flopy.mf6.MFSimulation.load``, which
  drops it silently (measured; pinned as the documented limit);
- MF6's requirements: a CG solution is refused at build (MF6 stops), the automatic IMS
  picks BICGSTAB, a model without NEWTON warns, a lake with no connection is refused;
- the routes that REBUILD a LAK package rather than write it carry the class: FloPy's
  model splitter (``ParallelModelWorkflow.split_model``), a package artifact, and
  ``model.config``, which shows it as the lake's ``implicit`` setting so a config diff
  tells two runs apart by it (each dropped it silently, the review of ledger 326);
- a library lake reached through ``mf.ref``: the IMS ``mf.simulation`` generated is given
  BICGSTAB once the build resolves it, and a declared CG solver is still refused;
- the executable probe: run once per executable file, cached, and a build that does not
  accept the option WARNS naming it and its version, and the run goes ahead;
- and MF6 (slow): an implicit lake runs to normal termination at the default
  formulation's stage, with its budget closed.
"""

from __future__ import annotations

import contextlib
import io
import os
import pathlib
import stat
import warnings

import flopy
import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point, box

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.lak_implicit import (
    ModflowGwflakImplicit,
    has_implicit,
    implicit_support,
    insert_implicit,
    warn_if_implicit_unsupported,
)
from myflopy.specs import package_digest


def _count(path) -> int:
    """How many lines of the file's OPTIONS block are IMPLICIT (and none outside it)."""

    text = pathlib.Path(path).read_text().upper()
    options = text.split("BEGIN OPTIONS", 1)[1].split("END OPTIONS", 1)[0]
    assert text.count("IMPLICIT") == options.count("IMPLICIT"), "only inside the block"
    return sum(1 for line in options.splitlines() if line.split()[:1] == ["IMPLICIT"])


def _pit(stack_top=100.0):
    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs,
                            geometry=[Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 32)])
    result = (mf.LayerStack(top=mf.Flat(stack_top)).add("upper", bottom=mf.Flat(82.0))
              .add("lower", bottom=mf.Flat(40.0))
              .excavate(bottom=82.0, within=pond, name="pit").build(grid))
    return grid, pond, result


def _simulation(grid, pond, result, *, implicit=True, newton=True, ims=None, inflow=None, nper=1):
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    lake = mf.lak(context=context, nper=nper, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=3000.0, implicit=implicit,
                  unexcavated="raise", inflow=inflow)
    model = mf.gwf("m", context=context, packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 94.0] for c in edge]}),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, transient={0: True}) if nper > 1 else mf.sto(steady_state={0: True}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ], **({"newtonoptions": "NEWTON"} if newton else {}))
    packages = [mf.tdis(nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)]
    if ims is not None:
        packages.append(ims)
    return mf.simulation(model, name="s", packages=packages)


@pytest.fixture(scope="module")
def pit():
    return _pit()


# -- the declaration ------------------------------------------------------------------------


def test_implicit_rides_in_the_options_and_moves_the_digest(pit):
    grid, pond, result = pit
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    kw = dict(context=context, nper=1, lakes=pond, lake_id_field="name", starting_stage=95.0)
    on, off = mf.lak(**kw, implicit=True), mf.lak(**kw)
    assert type(on.options["implicit"]) is bool and on.options["implicit"] is True
    assert "implicit" not in off.options, "False writes no key: no existing digest moves"
    assert package_digest(on) != package_digest(off)
    assert package_digest(off) == package_digest(mf.lak(**kw, implicit=False))
    assert package_digest(mf.lak(**kw, implicit=np.bool_(True))) == package_digest(on)
    # MF6 does not write the package-convergence CSV under IMPLICIT (the stage is in the
    # IMS check), so its file record is not asked for.
    assert "package_convergence_filerecord" in off.options
    assert "package_convergence_filerecord" not in on.options
    for back in (mf.PackageSpec.from_dict(on.to_dict()),):
        assert back.options["implicit"] is True and package_digest(back) == package_digest(on)
    with pytest.raises(TypeError, match="implicit must be True or False"):
        mf.lak(**kw, implicit="yes")


def test_insert_implicit_writes_it_once_inside_the_options_block(tmp_path):
    path = tmp_path / "m.lak"
    path.write_text("# header\nBEGIN options\n  BOUNDNAMES\nEND options\n\nBEGIN dimensions\n"
                    "  NLAKES 1\nEND dimensions\n")
    assert insert_implicit(path) and not insert_implicit(path)
    assert _count(path) == 1
    lines = path.read_text().splitlines()
    assert lines.index("  IMPLICIT") == lines.index("END options") - 1
    # A lower-case keyword a user wrote counts: nothing is added.
    path.write_text("BEGIN OPTIONS\n  implicit\nEND OPTIONS\n")
    assert not insert_implicit(path)
    # No options block at all: one is made, above the first block.
    path.write_text("BEGIN dimensions\n  NLAKES 1\nEND dimensions\n")
    assert insert_implicit(path) and _count(path) == 1
    assert path.read_text().startswith("BEGIN options\n  IMPLICIT\nEND options\n")
    assert has_implicit(path) and not has_implicit(tmp_path / "missing.lak")


# -- every write path ---------------------------------------------------------------------


def _lak_file(workspace) -> pathlib.Path:
    return pathlib.Path(workspace) / "m.lak"


def test_every_write_path_writes_it_exactly_once(pit, tmp_path):
    grid, pond, result = pit
    simulation = _simulation(grid, pond, result)
    project = mf.Project(tmp_path / "p", name="p")
    with contextlib.redirect_stdout(io.StringIO()):
        run = project.prepare_run("r", simulation)
        run.write(silent=True)
    lak = _lak_file(run.workspace)
    assert isinstance(run.simulation.get_model("m").get_package("lak"), ModflowGwflakImplicit)
    assert _count(lak) == 1, "Run.write"
    run.simulation.write_simulation(silent=True)
    assert _count(lak) == 1, "a second, direct write_simulation"
    run.write(silent=True)
    assert _count(lak) == 1, "Run.write again"
    # A BuiltSimulation, written by hand.
    built = simulation.build_flopy(tmp_path / "flopy")
    built.simulation.write_simulation(silent=True)
    assert _count(_lak_file(tmp_path / "flopy")) == 1, "build_flopy + write_simulation"
    # A run reopened from its frozen snapshot rebuilds the package through the builder.
    reopened = mf.load_run(run.workspace)
    with contextlib.redirect_stdout(io.StringIO()):
        reopened.write(silent=True)
    assert isinstance(reopened.simulation.get_model("m").get_package("lak"), ModflowGwflakImplicit)
    assert _count(lak) == 1, "a reopened run"
    # The declared False writes none.
    off = _simulation(grid, pond, result, implicit=False).build_flopy(tmp_path / "off")
    off.simulation.write_simulation(silent=True)
    assert "IMPLICIT" not in _lak_file(tmp_path / "off").read_text().upper()


def test_a_loaded_run_keeps_it_through_myflopy_and_loses_it_through_bare_flopy(pit, tmp_path):
    grid, pond, result = pit
    with contextlib.redirect_stdout(io.StringIO()):
        run = mf.Project(tmp_path / "p", name="p").prepare_run("r", _simulation(grid, pond, result))
        run.write(silent=True)
    workspace = pathlib.Path(run.workspace)
    # Run.open: FloPy's loader, then keep_implicit.
    opened = mf.load_run(workspace)
    simulation = opened.open()
    assert isinstance(simulation.get_model("m").get_package("lak"), ModflowGwflakImplicit)
    simulation.set_sim_path(str(tmp_path / "opened"))
    simulation.write_simulation(silent=True)
    assert _count(tmp_path / "opened" / "m.lak") == 1
    # mf.load_mf6_run: the same, on its lazy load.
    loaded = mf.load_mf6_run(workspace)
    loaded.load_all()
    assert isinstance(loaded.sim.get_model("m").get_package("lak"), ModflowGwflakImplicit)
    # A bare FloPy load keeps nothing -- the documented limit (ledger 326).
    bare = flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0)
    assert type(bare.get_model("m").get_package("lak")) is flopy.mf6.ModflowGwflak
    bare.set_sim_path(str(tmp_path / "bare"))
    bare.write_simulation(silent=True)
    assert "IMPLICIT" not in (tmp_path / "bare" / "m.lak").read_text().upper()


def test_defining_the_subclass_leaves_flopys_loader_alone(pit, tmp_path):
    """FloPy registers every ``MFPackage`` subclass as the class its LOADER builds, so
    defining ``ModflowGwflakImplicit`` made every LAK file loaded afterwards one --
    and a rewrite then inserted IMPLICIT into files that never had it (measured). The
    registry is restored; a file WITHOUT the keyword loads as FloPy's own class, through
    myflopy's load paths too, and writes none."""

    grid, pond, result = pit
    assert flopy.mf6.mfbase.PackageContainer.packages_by_abbr["gwflak"] is flopy.mf6.ModflowGwflak
    with contextlib.redirect_stdout(io.StringIO()):
        run = mf.Project(tmp_path / "p", name="p").prepare_run(
            "r", _simulation(grid, pond, result, implicit=False))
        run.write(silent=True)
    workspace = pathlib.Path(run.workspace)
    for simulation in (flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0),
                       mf.load_run(workspace).open()):
        lak = simulation.get_model("m").get_package("lak")
        assert type(lak) is flopy.mf6.ModflowGwflak
        out = tmp_path / f"re{id(simulation)}"
        simulation.set_sim_path(str(out))
        simulation.write_simulation(silent=True)
        assert "IMPLICIT" not in (out / "m.lak").read_text().upper()


# -- the routes that rebuild a LAK package ---------------------------------------------------


def _written(pit, root, *, implicit=True, inflow=None):
    grid, pond, result = pit
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").prepare_run(
            "r", _simulation(grid, pond, result, implicit=implicit, inflow=inflow))
        run.write(silent=True)
    return run


def test_a_split_model_keeps_it(pit, tmp_path):
    """FloPy's splitter rebuilds every package through ``package_factory``, which answers
    FloPy's own ``ModflowGwflak``: a split lake was solved by the default formulation,
    silently (the review of ledger 326, measured). The split package takes its source's
    class; a source without IMPLICIT splits without it."""

    from myflopy.modflow.mf6.parallel import ParallelModelWorkflow

    grid = pit[0]
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    mask = (xy[:, 0] > np.median(xy[:, 0])).astype(int)
    for implicit in (True, False):
        # Inflow, because FloPy's splitter refuses a LAK with no period data.
        run = _written(pit, tmp_path / f"p{implicit}", implicit=implicit, inflow=[10.0])
        loaded = mf.load_mf6_run(run.workspace)
        loaded.load_all()
        out = tmp_path / f"split{implicit}"
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            split = ParallelModelWorkflow(loaded).split_model(
                workspace=out, mask=mask, write=True, save_mapping=False)
        lakes = sorted(out.glob("*.lak"))
        assert lakes, "the split wrote a lake"
        for path in lakes:
            if implicit:
                assert _count(path) == 1, path.name
            else:
                assert "IMPLICIT" not in path.read_text().upper(), path.name
        classes = {type(split.simulation.get_model(name).get_package("lak"))
                   for name in split.simulation.model_names
                   if split.simulation.get_model(name).get_package("lak") is not None}
        assert classes == {ModflowGwflakImplicit if implicit else flopy.mf6.ModflowGwflak}


def test_a_package_artifact_keeps_it(pit, tmp_path):
    """``build_package_artifact`` read every option through FloPy's datasets, and IMPLICIT
    is none: a restore built a plain ``ModflowGwflak``. It also carried no BOUNDNAMES, so
    restoring any ``mf.lak`` lake (its id is a boundname) died in FloPy's recarray."""

    from myflopy.project.components import apply_package_artifact, build_package_artifact

    for implicit in (True, False):
        run = _written(pit, tmp_path / f"p{implicit}", implicit=implicit)
        model = run.model("m")
        artifact = build_package_artifact(model, artifact_id="lake", package_name="lak")
        assert artifact.package_data.get("implicit") is (True if implicit else None)
        before = model.gwf.get_package("lak").connectiondata.get_data()
        model.gwf.remove_package("lak")
        apply_package_artifact(model, artifact, validate=True)
        restored = model.gwf.get_package("lak")
        assert type(restored) is (ModflowGwflakImplicit if implicit else flopy.mf6.ModflowGwflak)
        assert len(restored.connectiondata.get_data()) == len(before)
        out = tmp_path / f"restored{implicit}"
        model.sim.set_sim_path(str(out))
        model.sim.write_simulation(silent=True)
        text = (out / "m.lak").read_text().upper()
        if implicit:
            assert _count(out / "m.lak") == 1
            assert "PACKAGE_CONVERGENCE" not in text, "MF6 writes none under IMPLICIT"
        else:
            assert "IMPLICIT" not in text and "PACKAGE_CONVERGENCE" in text


def test_the_config_shows_it_and_a_diff_tells_two_runs_apart(pit, tmp_path):
    """``model.config`` reads each package's options from FloPy's datasets, so two runs
    differing ONLY in ``implicit=True`` compared EQUAL -- ledger 286's NEWTON trap,
    reopened. It reads the class now: a built model, a loaded one, and the diff."""

    from myflopy.project.model_group import ModelGroup

    runs = {implicit: _written(pit, tmp_path / f"p{implicit}", implicit=implicit)
            for implicit in (False, True)}
    assert runs[True].model("m").config.section("lak")["implicit"] is True
    assert "implicit" not in runs[False].model("m").config.section("lak")
    loaded = {}
    for implicit, run in runs.items():
        loaded[implicit] = mf.load_mf6_run(run.workspace)
        loaded[implicit].load_all()
    assert loaded[True].config.section("lak")["implicit"] is True
    difference = ModelGroup({"off": loaded[False], "on": loaded[True]}).diff().config.settings()
    rows = difference[difference["section"] == "lak"]
    assert rows[["setting", "model_value"]].values.tolist() == [["implicit", True]]
    assert len(difference) == 1, "nothing else differs between the two runs"


# -- MF6's requirements ---------------------------------------------------------------------


def test_a_cg_solution_is_refused_and_the_automatic_one_is_bicgstab(pit, tmp_path):
    grid, pond, result = pit
    for ims in (mf.ims(models=["m"], complexity="SIMPLE"),
                mf.ims(models=["m"]),
                mf.ims(models=["m"], complexity="COMPLEX", linear_acceleration="CG")):
        with pytest.raises(ValueError, match="CONJUGATE GRADIENT"):
            _simulation(grid, pond, result, ims=ims).build_flopy(tmp_path / "cg")
    for ims in (mf.ims(models=["m"], complexity="SIMPLE", linear_acceleration="BICGSTAB"),
                mf.ims(models=["m"], complexity="MODERATE")):
        _simulation(grid, pond, result, ims=ims).build_flopy(tmp_path / "ok")
    # No solver declared, and no NEWTON: the automatic IMS still has to be BICGSTAB.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        built = _simulation(grid, pond, result, newton=False).build_flopy(tmp_path / "auto")
    ims = next(p for p in built.simulation.sim_package_list if p.package_type == "ims")
    assert str(ims.linear_acceleration.get_data()).upper() == "BICGSTAB"
    # ...and CG without an implicit lake is untouched.
    _simulation(grid, pond, result, implicit=False,
                ims=mf.ims(models=["m"], complexity="SIMPLE")).build_flopy(tmp_path / "plain")


def test_a_library_lake_gets_bicgstab_from_the_generated_solver(pit, tmp_path):
    """``mf.simulation`` cannot see through an ``mf.ref``, so a library lake got the SIMPLE
    solver it generates, CG by default, and the build refused a solver the caller never
    declared (the review of ledger 326, measured). The generated solver is marked, and the
    build gives it BICGSTAB once the ref resolves to an implicit lake -- on the reopened run
    and through a save too. A declared CG solver is still refused, and a library lake
    WITHOUT IMPLICIT leaves the generated solver alone."""

    from myflopy.grid_identity import GENERATED_SOLVER

    grid, pond, result = pit
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    project = mf.Project(tmp_path / "p", name="p")
    for implicit in (True, False):
        project.add_package(f"lak/{implicit}", mf.lak(
            context=context, nper=1, lakes=pond, lake_id_field="name", starting_stage=95.0,
            bed_leakance=3000.0, implicit=implicit, unexcavated="raise"))

    def model(name, *lake):
        # No NEWTON: with it `mf.simulation` picks BICGSTAB whatever the lake is, which
        # is not the case under test.
        return mf.gwf(name, context=context, packages=[
            mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1), mf.ic(strt=95.0),
            mf.chd(stress_period_data={0: [[(1, int(c)), 94.0] for c in edge]}),
            mf.sto(steady_state={0: True}), *lake,
            mf.oc(head_filerecord=f"{name}.hds", budget_filerecord=f"{name}.cbc",
                  saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])])

    def simulation(key, *packages):
        return mf.simulation(model("m", mf.ref(key)), name="s", packages=[mf.tdis(), *packages])

    def acceleration(simulation, solver="m_ims"):
        ims = next(p for p in simulation.sim_package_list if p.filename.startswith(solver))
        value = ims.linear_acceleration.get_data()
        return None if value is None else str(value).upper()

    def prepare(name, spec):
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*without NEWTON")
            return project.prepare_run(name, spec)

    declared = simulation("lak/True")
    generated = next(p for p in declared.packages if p.name == "m_ims")
    assert generated.metadata[GENERATED_SOLVER] is True
    assert "linear_acceleration" not in generated.options, "mf.simulation cannot see the ref"
    run = prepare("r", declared)
    assert acceleration(run.simulation) == "BICGSTAB"
    with contextlib.redirect_stdout(io.StringIO()):
        run.write(silent=True)
    ims_file = next(pathlib.Path(run.workspace).glob("*.ims"))
    assert "BICGSTAB" in ims_file.read_text().upper()
    # ...on the run reopened from its frozen snapshot (the marker is metadata, which the
    # snapshot's to_dict carries), and on the solver serialized by itself.
    reopened = mf.load_run(run.workspace)
    ims_file.unlink()
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*without NEWTON")
        reopened.write(silent=True)
    assert "BICGSTAB" in ims_file.read_text().upper()
    assert mf.PackageSpec.from_dict(generated.to_dict()).metadata[GENERATED_SOLVER] is True
    # A solver the caller declared is the caller's: CG is still refused, by name -- and so
    # is the generated one once the caller has chosen its acceleration.
    with pytest.raises(ValueError, match="CONJUGATE GRADIENT"):
        prepare("r_cg", simulation("lak/True", mf.ims(models=["m"], complexity="SIMPLE")))
    edited = declared.with_packages(generated.with_options(linear_acceleration="CG"))
    assert next(p for p in edited.packages if p.name == "m_ims").metadata[GENERATED_SOLVER]
    with pytest.raises(ValueError, match="CONJUGATE GRADIENT.*as declared"):
        prepare("r_edited", edited)
    # And a library lake without IMPLICIT leaves the generated solver as it was -- as does
    # a second model's own generated solver beside an implicit lake's.
    assert acceleration(prepare("r_plain", simulation("lak/False")).simulation) is None
    both = prepare("r_two", mf.simulation(model("m", mf.ref("lak/True")), model("n"), name="s",
                                          packages=[mf.tdis()])).simulation
    assert (acceleration(both), acceleration(both, "n_ims")) == ("BICGSTAB", None)


def test_no_newton_warns_and_a_lake_with_no_connection_is_refused(pit, tmp_path):
    grid, pond, result = pit
    with pytest.warns(UserWarning, match="IMPLICIT.*without NEWTON"):
        _simulation(grid, pond, result, newton=False).build_flopy(tmp_path / "n")
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    two = gpd.GeoDataFrame({"name": ["pit", "tank"]}, crs=grid.crs,
                           geometry=[pond.geometry.iloc[0], box(0, 0, 20, 20)])
    with pytest.raises(ValueError, match="lake 'tank' has no connection"):
        mf.lak(context=context, nper=1, lakes=two, lake_id_field="name", starting_stage=95.0,
               lake_bottom={"tank": 82.0}, implicit=True, unexcavated="ignore",
               connection_modes={"pit": "bathy", "tank": []})


# -- the executable probe -------------------------------------------------------------------


def _fake_mf6(root: pathlib.Path, *, version="mf6: 6.6.2") -> tuple[pathlib.Path, pathlib.Path]:
    """An executable that answers --version and refuses IMPLICIT the way an old MF6 would."""

    counter = root / "probes.txt"
    exe = root / "old_mf6"
    exe.write_text(
        "#!/bin/sh\n"
        f'if [ "$1" = "--version" ]; then echo "{version}"; exit 0; fi\n'
        f'echo run >> "{counter}"\n'
        'printf "ERROR REPORT:\\n  1. UNKNOWN LAK OPTION: IMPLICIT\\n" > mfsim.lst\n'
        'echo "ERROR REPORT:"; echo "  1. UNKNOWN LAK OPTION: IMPLICIT"\n'
        "exit 1\n"
    )
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return exe, counter


def test_the_probe_runs_once_per_executable_file(tmp_path):
    exe, counter = _fake_mf6(tmp_path)
    first = implicit_support(str(exe))
    assert first.supported is False and first.version == "mf6: 6.6.2"
    assert "UNKNOWN LAK OPTION: IMPLICIT" in first.detail
    assert implicit_support(str(exe)) is first and counter.read_text().count("run") == 1
    os.utime(exe, ns=(exe.stat().st_atime_ns, exe.stat().st_mtime_ns + 10**9))   # "rebuilt"
    implicit_support(str(exe))
    assert counter.read_text().count("run") == 2
    missing = implicit_support(str(tmp_path / "no_such_mf6"))
    assert missing.supported is None and missing.path is None


def test_an_executable_without_the_option_is_warned_about_and_the_run_goes_ahead(pit, tmp_path):
    grid, pond, result = pit
    exe, _ = _fake_mf6(tmp_path)
    run = mf.Project(tmp_path / "p", name="p").prepare_run(
        "r", _simulation(grid, pond, result), executable=str(exe))
    with pytest.warns(UserWarning, match=r"does not accept it") as caught, \
            contextlib.redirect_stdout(io.StringIO()):
        success, report = run.execute()
    text = str(caught[0].message)
    assert str(exe) in text and "mf6: 6.6.2" in text and "UNKNOWN LAK OPTION: IMPLICIT" in text
    assert not success, "it ran, and the executable stopped with its own message"
    # ...and a run loaded with mf.load_mf6_run, rewritten and run again, is probed too.
    loaded = mf.load_mf6_run(run.workspace)
    loaded.load_all()
    loaded.sim.exe_name = str(exe)
    with pytest.warns(UserWarning, match=r"does not accept it"), \
            contextlib.redirect_stdout(io.StringIO()):
        loaded.run_simulation()
    built = _simulation(grid, pond, result).build_flopy(tmp_path / "flopy")
    with pytest.warns(UserWarning, match="could not be checked"):
        warn_if_implicit_unsupported(built.simulation, str(tmp_path / "no_such_mf6"))
    off = _simulation(grid, pond, result, implicit=False).build_flopy(tmp_path / "off")
    assert warn_if_implicit_unsupported(off.simulation, str(exe)) is None


# -- MF6 ------------------------------------------------------------------------------------


def _stage_and_discrepancy(workspace):
    workspace = pathlib.Path(workspace)
    stage = flopy.utils.HeadFile(str(workspace / "m_lak_stage.lak"), text="STAGE",
                                 precision="double").get_alldata().ravel()
    listing = (workspace / "m.lst").read_text()
    discrepancy = [abs(float(line.split()[-1])) for line in listing.splitlines()
                   if "PERCENT DISCREPANCY" in line]
    return stage, discrepancy


@pytest.mark.slow
def test_mf6_an_implicit_lake_matches_the_default_formulation(pit, tmp_path):
    """Transient, the pit filled at 2,000 ft3/d for three days: MF6 solves the stage in the
    matrix (the listing says so) and lands on the default formulation's stages, with the
    GWF and LAK budgets closed. Normal termination alone is not the assertion."""

    grid, pond, result = pit
    ims = mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-8, inner_dvclose=1e-10,
                 linear_acceleration="BICGSTAB")
    stages = {}
    for implicit in (False, True):
        simulation = _simulation(grid, pond, result, implicit=implicit, ims=ims,
                                 inflow=[2000.0, 2000.0, 2000.0], nper=3)
        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            run = mf.Project(tmp_path / f"p{implicit}", name="p").run("r", simulation)
        assert run.success, "\n".join(run.report[-20:])
        stage, discrepancy = _stage_and_discrepancy(run.workspace)
        assert len(discrepancy) == 6 and max(discrepancy) < 0.1, "GWF and LAK, each period"
        listing = (pathlib.Path(run.workspace) / "m.lst").read_text()
        assert ("IMPLICIT FORMULATION" in listing) is implicit
        stages[implicit] = stage
    assert _count(pathlib.Path(run.workspace) / "m.lak") == 1
    np.testing.assert_allclose(stages[True], stages[False], atol=1e-6)   # 1.0e-8 measured
    assert stages[True][-1] > stages[True][0] + 0.01, "the pit filled: the comparison means something"
