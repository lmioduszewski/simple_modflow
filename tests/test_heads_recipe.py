"""``mf.ic(strt=mf.Heads(old_run))``: a finished run's heads as a RECIPE (ledger 596).

``model.hds.remap`` (ledger 595) hands back an array, and an array makes a simulation
unsaveable and pins it to its mesh. The recipe is the same read and the same engine,
resolved at BUILD against whatever grid the model has -- so the model saves, moves
(``rebuilt``), freezes what it carried, bundles the two files it reads and drifts by their
content. What this file pins:

* **It IS the eager remap**, on the file's grid, on any mesh -- one helper, not a copy.
* **The run is resolved at the DECLARATION**: head file, grid file, periods, CRS, hashes;
  a bad ``per``, an unconverged step, a DISU or a missing file is refused there.
* **The layer map is checked against the BUILT model** -- which ``remap`` cannot do.
* **It fills ``mf.ic(strt=)`` of a GWF model and nothing else.**
* **Save, YAML, digest-by-content, freeze, bundle, drift, move**, as every recipe.
"""

from __future__ import annotations

import shutil
import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.grid import remap as remap_module
from myflopy.sources import FIELD_METHODS, FIELD_UNCOVERED, Heads
from myflopy.specs import ModelContext, SimulationSpec, _json_value, _spec_value, package_digest

pytestmark = pytest.mark.slow                       # every test here reads an MF6 run


def _boundary(root, side=600.0):
    path = root / f"boundary_{int(side)}.gpkg"
    if not path.exists():
        gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, side, side)], crs=2927).to_file(path)
    return path


def _mesh(root, max_area, tag, side=600.0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=_boundary(root, side)), max_area=max_area
        ).resolve(workspace=root / tag)


def _flow_model(vor, strt, *, nlay=2, idomain=None, name="m", context=False,
                disv_options=None, gwf_options=None):
    """CHD 9,990 west / 10,010 east in the bottom layer; layer 0's north-east corner DRY."""

    props = vor.get_disv_gridprops()
    x = np.asarray(vor.centroids_x, dtype=float)
    y = np.asarray(vor.centroids_y, dtype=float)
    raised = (x > 400) & (y > 400)
    edges = vor.gdf_vorPolys.geometry.bounds
    east_edge = edges.maxx.max()
    keep = np.ones(vor.ncpl, dtype=bool) if idomain is None else np.asarray(idomain)[nlay - 1] > 0
    chd = [[(nlay - 1, int(c)), 9_990.0] for c in np.flatnonzero((edges.minx <= 1e-6) & keep)]
    chd += [[(nlay - 1, int(c)), 10_010.0]
            for c in np.flatnonzero((edges.maxx >= east_edge - 1e-6) & keep)]
    botm = [np.where(raised, 10_050.0, 9_950.0)] + [
        np.full(vor.ncpl, 9_800.0 - 50.0 * k) for k in range(nlay - 1)
    ]
    disv = dict(nlay=nlay, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10_100.0, botm=botm)
    if idomain is not None:
        disv["idomain"] = idomain
    disv.update(disv_options or {})
    where = dict(context=ModelContext(grid=vor)) if context else dict(grid=vor)
    return mf.gwf(name, **where, **(gwf_options or {}), packages=[
        mf.disv(**disv), mf.npf(k=10.0, icelltype=1), mf.ic(strt=strt),
        mf.chd(stress_period_data={0: chd}),
        mf.oc(head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]),
    ])


def _simulation(model, name="s"):
    return mf.simulation(model, name=name, packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["m"])])


def _run(project, name, simulation):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Run '.*' built but its frozen snapshot", UserWarning)
        run = project.run(name, simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """An OLD run on one mesh, and two other meshes of the same square."""

    root = tmp_path_factory.mktemp("heads_recipe")
    old_mesh = _mesh(root, 12_000.0, "old")
    project = mf.Project(root / "old_project", name="old")
    old_run = _run(project, "old", _simulation(_flow_model(old_mesh, 10_000.0)))
    return {
        "root": root,
        "old_mesh": old_mesh,
        "old_run": old_run,
        "new": _mesh(root, 5_000.0, "new"),
        "third": _mesh(root, 20_000.0, "third"),
    }


def _strt(run, name="m"):
    return np.asarray(run.model(name).fl.ic.strt.array, dtype=float)


def _eager(world, grid, **kwargs):
    """The eager remap on the grid MF6 WROTE (what the recipe reads), warnings off."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return mf.load_mf6_run(world["old_run"].workspace).hds.remap(grid, **kwargs)


# -- the point: it IS the eager remap, on any mesh -------------------------------------------


@pytest.mark.parametrize("target", ["new", "third"])
def test_one_declaration_builds_the_eager_remap_on_any_mesh(world, target, tmp_path):
    grid = world[target]
    project = mf.Project(tmp_path / "p", name="p")
    run = _run(project, "r", _simulation(_flow_model(grid, mf.Heads(world["old_run"]))))
    np.testing.assert_allclose(_strt(run), _eager(world, grid), atol=1e-9)


@pytest.mark.parametrize("given", ["run", "model", "directory", "text"])
def test_a_run_a_model_or_a_directory_declares_the_same_recipe(world, given):
    run = world["old_run"]
    value = {"run": run, "model": run.model(), "directory": run.workspace,
             "text": str(run.workspace)}[given]
    stamped = mf.ic(strt=mf.Heads(value)).inputs["strt"]
    assert stamped.resolved and stamped.model == "m" and stamped.crs == "EPSG:2927"
    assert stamped.path == (run.workspace / "m.hds").resolve()
    assert stamped.grid_path == (run.workspace / "m.disv.grb").resolve()
    assert stamped.perlen == (1.0,) and stamped.source_sha and stamped.grid_sha


def test_the_vocabularies_are_the_engine_s():
    assert FIELD_METHODS == remap_module.FIELD_METHODS
    assert FIELD_UNCOVERED == remap_module.UNCOVERED_MODES


# -- refused at the declaration ----------------------------------------------------------


def test_what_is_not_a_finished_run_is_refused_at_the_declaration(world, tmp_path):
    run = world["old_run"]
    with pytest.raises(FileNotFoundError, match="no run at"):
        mf.ic(strt=mf.Heads(tmp_path / "absent"))
    with pytest.raises(ValueError, match="not a file"):
        mf.ic(strt=mf.Heads(run.workspace / "m.hds"))
    (tmp_path / "empty").mkdir()
    with pytest.raises(FileNotFoundError, match="no mfsim.nam"):
        mf.ic(strt=mf.Heads(tmp_path / "empty"))
    with pytest.raises(ValueError, match="no model of that name"):
        mf.ic(strt=mf.Heads(run, model="nope"))
    with pytest.raises(ValueError, match="is not a stress period"):
        mf.ic(strt=mf.Heads(run, per=4))
    with pytest.raises(TypeError, match="give the finished run"):
        mf.Heads(object())


def test_a_step_mf6_did_not_converge_on_is_refused_at_the_declaration(world, tmp_path):
    copy = tmp_path / "copy"
    shutil.copytree(world["old_run"].workspace, copy)
    with (copy / "mfsim.lst").open("a") as listing:
        listing.write("\n Solution Group 1 did not converge for stress period 1 and time step 1\n")
    with pytest.raises(ValueError, match="did NOT converge on"):
        mf.ic(strt=mf.Heads(copy))


@pytest.mark.parametrize(
    "declare, says",
    [
        (lambda h: mf.npf(k=h), r"mf.npf\(k=mf.Heads"),
        (lambda h: mf.sto(sy=h), r"mf.sto\(sy=mf.Heads"),
        (lambda h: mf.ic(strt=[h, 9_990.0]), "fills every layer itself"),
        (lambda h: mf.rch.array(recharge=h), "STARTING heads"),
        (lambda h: h * 2.0, "cannot be multiplied"),
        (lambda h: 2.0 * h, "cannot be multiplied"),
        (lambda h: mf.ArrayProduct.of(2.0, h), "cannot be a factor"),
    ],
)
def test_it_fills_mf_ic_strt_and_nothing_else(world, declare, says):
    with pytest.raises(TypeError, match=says):
        declare(mf.Heads(world["old_run"]))


def test_a_recipe_in_options_is_refused(world):
    with pytest.raises(TypeError, match=r"holds mf.Heads\(\.\.\.\).*mf.ic\(strt=mf.Heads"):
        mf.PackageSpec("ic", lambda model, **kw: None, {"strt": mf.Heads(world["old_run"])})


# -- refused at the build ------------------------------------------------------------------


def test_the_layer_map_is_checked_against_the_built_model(world, tmp_path):
    old = world["old_run"]
    grid = world["new"]

    def build(recipe, nlay):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"Run '.*' built but its frozen snapshot", UserWarning)
            return _simulation(_flow_model(grid, recipe, nlay=nlay)).build(tmp_path / f"b{nlay}")

    with pytest.raises(ValueError, match=r"2 layer\(s\) and this model has 3.*layers=\[0, 1, 1\]"):
        build(mf.Heads(old), 3)
    with pytest.raises(ValueError, match=r"layers=\[0, 0, 0\]"):
        build(mf.Heads(old, layers=0), 3)
    with pytest.raises(ValueError, match="one old layer per layer of this model"):
        build(mf.Heads(old, layers=[0, 1]), 3)
    built = build(mf.Heads(old, layers=[0, 0, 1]), 3)
    strt = np.asarray(built.simulation.get_model("m").ic.strt.array)
    expected = _eager(world, grid, layers=[0, 0, 1])
    np.testing.assert_allclose(strt, expected, atol=1e-9)
    np.testing.assert_array_equal(strt[0], strt[1])


def test_a_transport_model_is_refused(world, tmp_path):
    """``mf.ic`` on a GWT model is a starting CONCENTRATION; a head is not one."""

    grid = world["new"]
    props = grid.get_disv_gridprops()
    transport = mf.gwt("t", context=ModelContext(grid=grid), packages=[
        mf.disv(nlay=2, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10_100.0,
                botm=[9_950.0, 9_800.0]),
        mf.ic(strt=mf.Heads(world["old_run"])),
    ])
    simulation = SimulationSpec("st", models=[transport], packages=[
        mf.tdis(nper=1), mf.ims(name="ims_t", models=["t"])])
    with pytest.raises(TypeError, match="a GWF model's mf.ic"):
        simulation.build(tmp_path / "gwt")


def test_a_cell_the_new_model_keeps_inactive_is_not_counted(world, tmp_path):
    """The new idomain reaches the engine: cells off the old model and inactive say nothing."""

    big = _mesh(world["root"], 20_000.0, "big", side=800.0)
    footprint = world["old_mesh"].gdf_vorPolys.geometry.union_all()
    outside = big.gdf_vorPolys.geometry.intersection(footprint).area.to_numpy() < 1e-6
    assert outside.any()
    for idomain, warned in ((None, True), (np.vstack([~outside, ~outside]).astype(int), False)):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _simulation(_flow_model(big, mf.Heads(world["old_run"]), idomain=idomain)).build(
                tmp_path / f"active_{warned}"
            )
        assert any("lie outside the old model" in str(w.message) for w in caught) is warned


# -- it saves, round-trips, digests by content, moves, freezes, bundles and drifts ---------


def _registered(root, grid, recipe):
    project = mf.Project(root, name="p")
    project.add_grid("mesh", mf.GridSpec.from_object(grid))
    simulation = _simulation(mf.gwf("m", grid=mf.grid_ref("mesh"), packages=[
        mf.disv(mf.LayerStack(top=mf.Surface.flat(10_100.0))
                .add("upper", thickness=150.0).add("lower", thickness=150.0)),
        mf.npf(k=10.0, icelltype=1), mf.ic(strt=recipe),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ]))
    project.add_simulation(simulation)
    return project, simulation


def test_a_project_holding_the_recipe_saves_loads_and_builds_the_same_heads(world, tmp_path):
    project, _ = _registered(tmp_path / "p", world["new"], mf.Heads(world["old_run"]))
    assert project.save().simulations == 1
    loaded = mf.Project.load(tmp_path / "p")
    simulation = next(iter(loaded.simulations.values()))
    recipe = next(p for p in simulation.models[0].packages if p.name == "ic").inputs["strt"]
    assert isinstance(recipe, Heads) and recipe.resolved
    built = loaded.prepare_run("r", simulation).build()
    strt = np.asarray(built.simulation.get_model("m").ic.strt.array)
    np.testing.assert_allclose(strt, _eager(world, world["new"]), atol=1e-9)


def test_the_recipe_round_trips_through_json_and_yaml(world, tmp_path):
    recipe = mf.ic(strt=mf.Heads(world["old_run"], per=0, layers=[1, 0], method="area")).inputs["strt"]
    assert _spec_value(_json_value(recipe)) == recipe
    _, simulation = _registered(tmp_path / "p", world["new"], mf.Heads(world["old_run"], layers=[0, 1]))
    assert SimulationSpec.from_yaml(simulation.to_yaml()).to_dict() == simulation.to_dict()


def test_the_digest_reads_the_files_not_where_they_are(world, tmp_path):
    copy = tmp_path / "copy"
    shutil.copytree(world["old_run"].workspace, copy)
    here = package_digest(mf.ic(strt=mf.Heads(world["old_run"])))
    assert package_digest(mf.ic(strt=mf.Heads(copy))) == here
    assert package_digest(mf.ic(strt=mf.Heads(copy, per=0, layers=[1, 1]))) != here


def test_a_model_starting_from_the_recipe_moves_rebuilt(world, tmp_path):
    simulation = _simulation(_flow_model(world["new"], mf.Heads(world["old_run"])))
    report = {row["package"]: row["how"] for row in simulation.move_report()}
    assert report["ic"] == "rebuilt"


def test_a_frozen_run_rebuilds_with_the_old_run_deleted(world, tmp_path):
    """The run freezes the CARRIED heads, not the recipe (ledgers 251, 270)."""

    from myflopy.project.freeze import freeze_simulation, thaw_simulation

    copy = tmp_path / "old_copy"
    shutil.copytree(world["old_run"].workspace, copy)
    simulation = _simulation(_flow_model(world["new"], mf.Heads(copy), context=True))
    built = simulation.build(tmp_path / "built").built
    carried = built.models["m"].resolved["ic"].options["strt"]
    assert isinstance(carried, np.ndarray)
    freeze_simulation(
        simulation, tmp_path / "snap", run_name="s",
        contexts={name: model.context for name, model in built.models.items()},
        resolved={name: dict(model.resolved or {}) for name, model in built.models.items()},
    )
    shutil.rmtree(copy)
    packages = {p.name: p for p in thaw_simulation(tmp_path / "snap").models[0].packages}
    np.testing.assert_array_equal(np.asarray(packages["ic"].options["strt"]), carried)


def test_bundling_copies_the_head_file_and_the_grid_file(world, tmp_path):
    copy = tmp_path / "outside" / "old"
    shutil.copytree(world["old_run"].workspace, copy)
    project, _ = _registered(tmp_path / "p", world["new"], mf.Heads(copy))
    report = project.save(bundle_sources=True)
    names = {str(destination).rsplit("/", 1)[-1] for _, destination in report.bundled}
    assert {"m.hds", "m.disv.grb"} <= names
    shutil.rmtree(copy)
    loaded = mf.Project.load(project.root)
    built = loaded.prepare_run("r", next(iter(loaded.simulations.values()))).build()
    strt = np.asarray(built.simulation.get_model("m").ic.strt.array)
    np.testing.assert_allclose(strt, _eager(world, world["new"]), atol=1e-9)


def test_drift_reads_both_files_and_notices_an_edit(world, tmp_path):
    copy = tmp_path / "old_copy"
    shutil.copytree(world["old_run"].workspace, copy)
    project, simulation = _registered(tmp_path / "p", world["new"], mf.Heads(copy))
    run = project.prepare_run("r1", simulation)
    run.build()
    assert set(project.runs.drift()["source"]) == {"ok"}
    read = run.metadata["built_from"]["ic"]
    assert {path.rsplit("/", 1)[-1] for path in read} == {"m.hds", "m.disv.grb"}
    with (copy / "m.hds").open("ab") as head_file:
        head_file.write(b"\0")
    assert set(project.runs.drift()["source"]) == {"edited (1 of 2)"}


def test_a_step_that_stopped_converging_after_the_declaration_is_refused_at_build(world, tmp_path):
    """The old run re-run between declaring and building: the build reads its listing too."""

    copy = tmp_path / "copy"
    shutil.copytree(world["old_run"].workspace, copy)
    simulation = _simulation(_flow_model(world["new"], mf.Heads(copy)))
    with (copy / "mfsim.lst").open("a") as listing:
        listing.write("\n Solution Group 1 did not converge for stress period 1 and time step 1\n")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Run '.*' built but its frozen snapshot", UserWarning)
        with pytest.raises(ValueError, match="did NOT converge on"):
            simulation.build(tmp_path / "built")



# -- the review of ledger 596: layouts, grid files, cut records, CRS, digests, options ------


def test_a_coupled_run_s_flow_model_is_found_and_its_heads_carried(world, tmp_path):
    """A coupled run keeps each model in a subdirectory, and its name files name every
    package relative to the SIMULATION directory (``flow/flow.disv``): joined onto the
    model's folder the paths doubled and every coupled run was refused."""

    from myflopy.modflow.mf6.output_files import _name_file_oc
    from tests.test_gwt_gwe_results import _coupled_run

    run = _coupled_run(tmp_path, "gwt")
    oc_file = _name_file_oc(run.workspace, "flow")[1]
    assert oc_file is not None and oc_file.is_file(), "the coupled run's OC was not found"
    recipe = mf.ic(strt=mf.Heads(run)).inputs["strt"]
    assert recipe.model == "flow" and recipe.grid_path.is_file()
    with pytest.raises(ValueError, match="that model is GWT6"):
        mf.ic(strt=mf.Heads(run, model="trans"))
    big = _mesh(tmp_path, 20_000.0, "big", side=800.0)
    built = _simulation(_flow_model(big, mf.Heads(run), nlay=1)).build(tmp_path / "b")
    expected = mf.load_mf6_run(run.workspace, model_name="flow").hds.remap(big, layers=[0])
    np.testing.assert_allclose(np.asarray(built.simulation.get_model("m").ic.strt.array),
                               expected, atol=1e-9)


@pytest.mark.parametrize(
    "variant, grid_file",
    [
        (dict(disv_options={"grb_filerecord": "mesh.grb"}), "mesh.grb"),
        (dict(gwf_options={"model_rel_path": "flow"}), "m.disv.grb"),
    ],
)
def test_a_named_grid_file_and_a_model_folder_are_found(world, tmp_path, variant, grid_file):
    project = mf.Project(tmp_path / "p", name="p")
    run = _run(project, "old", _simulation(_flow_model(world["old_mesh"], 10_000.0, **variant)))
    recipe = mf.ic(strt=mf.Heads(run)).inputs["strt"]
    assert recipe.grid_path.name == grid_file and recipe.grid_path.is_file()
    built = _simulation(_flow_model(world["new"], mf.Heads(run))).build(tmp_path / "b")
    expected = mf.load_mf6_run(run.workspace).hds.remap(world["new"])
    np.testing.assert_allclose(np.asarray(built.simulation.get_model("m").ic.strt.array),
                               expected, atol=1e-9)


def test_a_run_written_with_nogrb_is_refused_by_name(world, tmp_path):
    project = mf.Project(tmp_path / "p", name="p")
    run = _run(project, "old", _simulation(
        _flow_model(world["old_mesh"], 10_000.0, disv_options={"nogrb": True})))
    with pytest.raises(ValueError, match="says NOGRB"):
        mf.ic(strt=mf.Heads(run))


def test_missing_files_and_a_disu_are_refused_at_the_declaration(world, tmp_path):
    def copy(tag):
        target = tmp_path / tag
        shutil.copytree(world["old_run"].workspace, target)
        return target

    no_grid = copy("no_grid")
    (no_grid / "m.disv.grb").unlink()
    with pytest.raises(FileNotFoundError, match="binary grid file .* does not exist"):
        mf.ic(strt=mf.Heads(no_grid))
    no_heads = copy("no_heads")
    (no_heads / "m.hds").unlink()
    with pytest.raises(FileNotFoundError, match="head file .* does not exist"):
        mf.ic(strt=mf.Heads(no_heads))
    disu = copy("disu")
    name_file = disu / "m.nam"
    name_file.write_text(name_file.read_text().replace("DISV6", "DISU6"))
    with pytest.raises(ValueError, match="is DISU"):
        mf.ic(strt=mf.Heads(disu))
    several = copy("several")
    simulation_file = several / "mfsim.nam"
    simulation_file.write_text(
        simulation_file.read_text().replace("  gwf6  m.nam  m", "  gwf6  m.nam  m\n  gwf6  m2.nam  m2")
    )
    with pytest.raises(ValueError, match=r"2 GWF models \(m, m2\).*model=\.\.\."):
        mf.ic(strt=mf.Heads(several))


def test_a_record_the_file_was_cut_inside_is_refused(world, tmp_path):
    """A run killed -- or still running -- mid-write leaves a header with no numbers."""

    cut = tmp_path / "cut"
    shutil.copytree(world["old_run"].workspace, cut)
    head_file = cut / "m.hds"
    head_file.write_bytes(head_file.read_bytes()[:-40])
    with pytest.raises(ValueError, match="stops inside the record"):
        mf.ic(strt=mf.Heads(cut))
    with pytest.raises(ValueError, match="stops inside the record"):
        mf.load_mf6_run(cut).hds.remap(world["new"])


def test_an_old_layer_the_run_lacks_is_refused_at_the_declaration(world):
    with pytest.raises(ValueError, match="outside the old model's 2 layers"):
        mf.ic(strt=mf.Heads(world["old_run"], layers=[0, 2]))


def test_the_stamped_fields_are_refused_beside_a_run_and_kept_on_a_reload(world, tmp_path):
    run = world["old_run"]
    for given in (dict(grid_path=run.workspace / "m.disv.grb"), dict(perlen=[1.0]),
                  dict(source_sha="0" * 16)):
        with pytest.raises(TypeError, match="filled in by mf.ic"):
            mf.Heads(run, **given)
    stamped = mf.ic(strt=mf.Heads(run, uncovered="raise", crs="EPSG:2927")).inputs["strt"]
    reloaded = _spec_value(_json_value(stamped))
    assert reloaded == stamped and reloaded.resolved and reloaded.uncovered == "raise"


def test_a_run_re_run_after_the_declaration_is_counted_on_its_own_periods(world, tmp_path):
    """The build counts ``per`` on the run's TDIS as it NOW is: here it says 2 periods
    where the head file holds 1, so ``per=-1`` is the second -- which the file lacks."""

    copy = tmp_path / "copy"
    shutil.copytree(world["old_run"].workspace, copy)
    simulation = _simulation(_flow_model(world["new"], mf.Heads(copy)))
    tdis = next(copy.glob("*.tdis"))
    text = tdis.read_text()
    text = text.replace("NPER  1", "NPER  2").replace("nper  1", "nper  2")
    lines = text.splitlines()
    at = next(i for i, line in enumerate(lines) if line.strip().upper().startswith("END PERIODDATA"))
    lines.insert(at, lines[at - 1])
    tdis.write_text("\n".join(lines) + "\n")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Run '.*' built but its frozen snapshot", UserWarning)
        with pytest.raises(ValueError, match="saved no record in stress period 1"):
            simulation.build(tmp_path / "built")


def test_a_domain_that_does_not_fit_is_refused_by_name(world):
    from myflopy.modflow.mf6.heads_recipe import resolve_heads

    recipe = mf.ic(strt=mf.Heads(world["old_run"])).inputs["strt"]
    with pytest.raises(ValueError, match="the model's domain has shape"):
        resolve_heads(recipe, world["new"], nlay=2, domain=np.ones((3, world["new"].ncpl)))


def test_the_crs_comes_from_the_call_the_model_the_manifest_or_the_default(world, tmp_path):
    run = world["old_run"]
    assert mf.ic(strt=mf.Heads(run, crs="EPSG:2926")).inputs["strt"].crs == "EPSG:2926"
    bare = tmp_path / "bare"
    shutil.copytree(run.workspace, bare)
    (bare / "run.json").unlink(missing_ok=True)
    viewed = mf.load_mf6_run(bare, crs="EPSG:2926")
    assert mf.ic(strt=mf.Heads(viewed)).inputs["strt"].crs == "EPSG:2926"
    assert mf.ic(strt=mf.Heads(bare)).inputs["strt"].crs == "EPSG:2927"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Run '.*' built but its frozen snapshot", UserWarning)
        with pytest.raises(ValueError, match="different coordinate systems"):
            _simulation(_flow_model(world["new"], mf.Heads(run, crs="EPSG:2926"))).build(
                tmp_path / "frame")


def test_the_digest_reads_both_files_and_not_the_tdis(world, tmp_path):
    def digest(path):
        return package_digest(mf.ic(strt=mf.Heads(path)))

    base = digest(world["old_run"].workspace)
    grid = tmp_path / "grid"
    shutil.copytree(world["old_run"].workspace, grid)
    with (grid / "m.disv.grb").open("ab") as handle:
        handle.write(b"\0")
    assert digest(grid) != base, "a different grid file must move the digest"
    stamp = tmp_path / "stamp"
    shutil.copytree(world["old_run"].workspace, stamp)
    tdis = next(stamp.glob("*.tdis"))
    tdis.write_text("# rewritten by FloPy at another time\n" + tdis.read_text())
    assert digest(stamp) == base, "FloPy's timestamp in the TDIS must not move it"


def test_every_option_reaches_the_build(world, tmp_path):
    grid = world["new"]

    def strt(recipe, nlay=2, tag="b"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            built = _simulation(_flow_model(grid, recipe, nlay=nlay)).build(tmp_path / tag)
        return np.asarray(built.simulation.get_model("m").ic.strt.array)

    np.testing.assert_allclose(strt(mf.Heads(world["old_run"], method="area"), tag="area"),
                               _eager(world, grid, method="area"), atol=1e-9)
    np.testing.assert_allclose(strt(mf.Heads(world["old_run"], layers=1), nlay=1, tag="one"),
                               _eager(world, grid, layers=[1]), atol=1e-9)
    big = _mesh(world["root"], 20_000.0, "big_options", side=800.0)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"Run '.*' built but its frozen snapshot", UserWarning)
        with pytest.raises(ValueError, match="uncovered='nearest'"):
            _simulation(_flow_model(big, mf.Heads(world["old_run"], uncovered="raise"))).build(
                tmp_path / "raise")


def test_a_gwe_model_is_refused_too(world, tmp_path):
    grid = world["new"]
    props = grid.get_disv_gridprops()
    energy = mf.gwe("e", context=ModelContext(grid=grid), packages=[
        mf.disv(nlay=2, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10_100.0,
                botm=[9_950.0, 9_800.0]),
        mf.ic(strt=mf.Heads(world["old_run"])),
    ])
    simulation = SimulationSpec("se", models=[energy], packages=[
        mf.tdis(nper=1), mf.ims(name="ims_e", models=["e"])])
    with pytest.raises(TypeError, match="a GWF model's mf.ic"):
        simulation.build(tmp_path / "gwe")


def test_a_run_that_started_from_the_recipe_reopens_with_the_old_run_deleted(world, tmp_path):
    """A real run: frozen by Run.build, reopened from its folder, rebuilt -- no old run."""

    copy = tmp_path / "old_copy"
    shutil.copytree(world["old_run"].workspace, copy)
    project = mf.Project(tmp_path / "p", name="p")
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)        # the snapshot must be written
        run = project.run("r", _simulation(_flow_model(world["new"], mf.Heads(copy), context=True)))
    assert run.success
    first = _strt(run)
    shutil.rmtree(copy)
    reopened = mf.load_run(run.workspace)
    rebuilt = reopened.build().simulation.get_model("m").ic.strt.array
    np.testing.assert_array_equal(np.asarray(rebuilt), first)


def test_a_relative_run_path_is_stored_absolute(world, tmp_path, monkeypatch):
    monkeypatch.chdir(world["old_run"].workspace.parent)
    recipe = mf.ic(strt=mf.Heads(world["old_run"].workspace.name)).inputs["strt"]
    assert recipe.path.is_absolute() and recipe.grid_path.is_absolute()
