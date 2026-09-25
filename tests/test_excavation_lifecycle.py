"""An excavated layer stack you can write down, save, archive, move and reopen (ledger 323).

The excavation design's C6 (section 6). Until it, an excavated stack REFUSED ``to_dict``
(ledger 320's interim), so ``mf.disv(<gridless excavated stack>)`` raised at the call and
such a model could not be saved. Its excavations now ride in the ``$layer_stack``
payload -- the footprint FILE under the ``"zone"`` key ``Surface.where`` already uses, so
every walker that finds a stack's files (the bundler, the build-time anchor,
``runs.drift()``'s ``built_from``, the hash stamping) finds it with no wiring of its own
-- and each site of the lifecycle is pinned here:

- the payload round-trips every footprint form (path + layer, a column floor, an
  ``mf.Area``, a GeoDataFrame, a bare geometry) through JSON and YAML on a real
  ``GridSpec.voronoi`` grid, in plain python types (``type(x) is``, ledger 265);
- a lazy ``mf.disv`` hashes the footprint and the bed, so moving them keeps the digest
  and editing either moves it -- and ``runs.drift()`` reports an edited footprint;
- ``save(bundle_sources=True)`` carries both, and the archive rebuilds with every
  source deleted; a run's frozen snapshot rebuilds without them too;
- one declaration builds on two meshes, and ``with_grid`` 28 -> 92 re-covers the new
  mesh's cells -- the same cells ``mf.lak`` would take there;
- a lazy DISV publishes the elevations and records of what it BUILT, not the grid's
  table, whatever ``attach=`` says.
"""

from __future__ import annotations

import json
import shutil
import warnings

import geopandas as gpd
import numpy as np
import pytest
import shapely as shp
import yaml

import myflopy as mf
from myflopy import LAKBuilder, ModelContext
from myflopy.grid_identity import excavation_cells
from myflopy.layers import LayerStack
from myflopy.project.freeze import thaw_simulation
from myflopy.specs import _json_value, _spec_value, package_digest, source_paths
from myflopy.workspace import Project
from tests.test_excavate import _raster

SIDE = 600.0


def _gis(root):
    """A boundary, a two-layer footprint file (with a per-feature floor column) and a
    bowl bed -- 92 ft at the centre, 96.5 at the pond's 150 ft edge."""

    root.mkdir(parents=True, exist_ok=True)
    boundary = root / "dom.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[shp.box(0, 0, SIDE, SIDE)], crs=2927).to_file(boundary)
    site = root / "site.gpkg"
    _write_pond(site)
    gpd.GeoDataFrame({"name": ["basin"], "floor": [96.0]}, geometry=[shp.box(470, 40, 590, 160)],
                     crs=2927).to_file(site, layer="basins")
    bowl = _raster(root / "bowl.tif", lambda x, y: 92.0 + ((x - 300) ** 2 + (y - 300) ** 2) / 5000.0)
    return type("Gis", (), {"root": root, "boundary": boundary, "site": site, "bowl": bowl})


def _write_pond(path, centre=(300.0, 300.0)):
    gpd.GeoDataFrame({"name": ["pond"], "floor": [93.0]},
                     geometry=[shp.Point(*centre).buffer(150)], crs=2927).to_file(path, layer="ponds")


@pytest.fixture
def gis(tmp_path):
    return _gis(tmp_path / "gis")


def _stack():
    return (mf.LayerStack(top=mf.Flat(100.0))
            .add("fill", thickness=5.0)
            .add("silt", thickness=1.0)
            .add("sand", thickness=30.0, split=2)
            .add("clay", thickness=20.0))


def _pond_stack(gis, **kw):
    kw.setdefault("bottom", mf.Raster(gis.bowl))
    return _stack().excavate(within=gis.site, layer="ponds", name="pond", **kw)


def _grid_spec(boundary, max_area):
    return mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=max_area)


def _simulation(gis, max_area=40_000.0, *, stack=None, grid=None, **disv_kwargs):
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(stack if stack is not None else _pond_stack(gis), **disv_kwargs),
        mf.npf(k=1.0), mf.ic(strt=99.0), mf.oc(head_filerecord="m.hds"),
    ])
    if grid is None:
        model = model.with_grid(_grid_spec(gis.boundary, max_area))
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _plain(value, where="payload"):
    """Every leaf a plain python type -- `type(x) is`, never isinstance (ledger 265)."""

    if type(value) is dict:
        for key, item in value.items():
            assert type(key) is str, f"{where}: key {key!r}"
            _plain(item, f"{where}.{key}")
    elif type(value) is list:
        for i, item in enumerate(value):
            _plain(item, f"{where}[{i}]")
    else:
        assert type(value) in (str, float, int, bool, type(None)), f"{where}: {type(value)}"


def _lake_cells(grid, gis):
    lake = LAKBuilder(context=ModelContext(grid=grid), nper=1, lakes=gis.site, layer="ponds",
                      lake_id_field="name", starting_stage=100.0, lake_bottom=90.0)
    return tuple(sorted({c for cells in lake.lake_cells.values() for c in cells}))


def _same_build(a, b):
    np.testing.assert_array_equal(a.top, b.top)
    np.testing.assert_array_equal(a.botm, b.botm)
    np.testing.assert_array_equal(a.idomain, b.idomain)
    assert list(a.excavations) == list(b.excavations)
    for name in a.excavations:
        assert a.excavations[name].cells == b.excavations[name].cells
        np.testing.assert_array_equal(a.excavations[name].bed, b.excavations[name].bed)


# -- the payload ---------------------------------------------------------------------


def _forms(gis):
    frame = gpd.read_file(gis.site, layer="ponds")
    bare = frame.geometry.iloc[0]
    return {
        "path": dict(within=gis.site, layer=np.str_("ponds")),
        "column": dict(within=gis.site, layer="ponds", bottom="floor"),
        "area": dict(within=mf.Area.coerce(frame, note="surveyed 2026"),
                     bottom=np.float64(93.0)),
        "frame": dict(within=frame),
        "geometry": dict(within=bare, min_remnant=np.float64(0.25)),
    }


@pytest.mark.parametrize("form", ["path", "column", "area", "frame", "geometry"])
def test_every_footprint_form_round_trips_through_json_and_yaml(triangulated_grid, gis, form):
    """On a real ``GridSpec.voronoi`` grid (int64 cell2d), with numpy scalars in."""

    kw = _forms(gis)[form]
    kw.setdefault("bottom", mf.Raster(gis.bowl))
    stack = _stack().excavate(name=np.str_("pond"), **kw)
    payload = stack.to_dict()
    _plain(payload)

    written = _json_value(stack)
    for route, revived in (
        ("json", _spec_value(json.loads(json.dumps(written)))),
        ("yaml", _spec_value(yaml.safe_load(yaml.safe_dump(written)))),
    ):
        assert isinstance(revived, LayerStack), route
        assert revived.to_dict() == payload, route
        _same_build(stack.build(triangulated_grid), revived.build(triangulated_grid))


def test_a_payload_names_the_file_by_path_and_an_area_inline(gis):
    by_path = _pond_stack(gis).to_dict()["excavations"][0]
    assert by_path["zone"] == gis.site.resolve().as_posix() and by_path["layer"] == "ponds"
    assert "area" not in by_path
    inline = _stack().excavate(93.0, _forms(gis)["area"]["within"], name="p").to_dict()
    entry = inline["excavations"][0]
    assert "zone" not in entry and entry["area"]["crs"] == "EPSG:2927"
    assert shp.from_wkt(entry["area"]["wkt"]).equals(
        gpd.read_file(gis.site, layer="ponds").geometry.iloc[0])


def test_from_dict_builds_directly_and_reads_no_file(tmp_path, gis):
    """A bundled payload names a project-relative path that only the build anchors; the
    load must not look for it against the working directory (ledger 293)."""

    payload = _pond_stack(gis).to_dict()
    payload["excavations"][0]["zone"] = "inputs/site.gpkg"          # does not exist here
    revived = LayerStack.from_dict(payload)
    assert revived._excavations[0].zone.as_posix() == "inputs/site.gpkg"


@pytest.mark.parametrize("edit, match", [
    (lambda p: p.update(ponds=[]), r"keys this version does not know: \['ponds'\]"),
    (lambda p: p["excavations"][0].update(depth=3), r"excavation entry carries keys.*'depth'"),
    (lambda p: p["excavations"][0].update(area={"kind": "Area", "wkt": "POINT (0 0)"}),
     "not both"),
    (lambda p: p["excavations"][0].pop("zone"), "not both, and not neither"),
])
def test_from_dict_refuses_what_it_does_not_know(gis, edit, match):
    payload = _pond_stack(gis).to_dict()
    edit(payload)
    with pytest.raises(ValueError, match=match):
        LayerStack.from_dict(payload)


def _inline_payload(gis):
    """A stack with one inline pond (``area``) and one pond from the file (``zone``)."""

    stack = _stack().excavate(93.0, _forms(gis)["area"]["within"], name="inline")
    return stack.excavate(within=gis.site, layer="ponds", bottom="floor", name="filed").to_dict()


def _second(entry, **changes):
    return [entry, {**entry, **changes}]


@pytest.mark.parametrize("edit, error, match", [
    # the same name twice: one record would silently replace the other at build while
    # both cuts stayed in the arrays -- and the overlap refusal, keyed by name, missed it
    (lambda p: p.update(excavations=_second(p["excavations"][0])), ValueError,
     "an excavation named 'inline' already exists"),
    (lambda p: p.update(excavations=_second(p["excavations"][0], bottom={"kind": "flat", "value": 85.0},
                                            area=mf.Area(shp.box(0, 0, 60, 60)).to_dict())),
     ValueError, "an excavation named 'inline' already exists"),
    (lambda p: p["excavations"][0].update(bottom=mf.Surface.offset_below(9.0).to_dict()),
     ValueError, "measured from the surface ABOVE"),
    (lambda p: p["excavations"][0].update(
        bottom=(mf.Surface.raster("x.tif", fill="propagate") - 2.0).to_dict()),
     ValueError, "built on a fill='propagate' source"),
    (lambda p: p["excavations"][0].update(layer="ponds"), ValueError, "selects a layer of a FILE"),
    (lambda p: p["excavations"][0].update(bottom={"column": "floor"}), ValueError,
     "a column floor needs within= to be a FILE"),
    (lambda p: p["excavations"][1].update(layer=0), TypeError, "VECTOR layer of the within file"),
    (lambda p: p["excavations"][0].update(name="sand"), ValueError, "names a unit or layer"),
    (lambda p: p["excavations"][0].update(name="sand_2"), ValueError, "names a unit or layer"),
    (lambda p: p["excavations"][0].update(name=""), ValueError, "may not be empty"),
    (lambda p: p["excavations"][0].update(min_remnant=-0.5), ValueError, "finite number >= 0"),
    (lambda p: p["excavations"][0].update(min_remnant="0.1"), TypeError, "must be a number"),
], ids=["same-name", "same-name-elsewhere", "relative-bed", "propagate-algebra",
        "layer-without-file", "column-without-file", "integer-layer", "unit-name",
        "split-layer-name", "empty-name", "negative-remnant", "string-remnant"])
def test_from_dict_refuses_what_excavate_refuses(gis, edit, error, match):
    """A hand-written payload (``SimulationSpec.from_yaml``) reaches ``from_dict``, which
    builds entries without reading their files -- but not without ``excavate``'s other
    refusals: the verb and the payload share them."""

    payload = _inline_payload(gis)
    assert LayerStack.from_dict(json.loads(json.dumps(payload))).to_dict() == payload
    edit(payload)
    with pytest.raises(error, match=match):
        LayerStack.from_dict(payload)


def test_the_verb_refuses_an_empty_name_as_the_payload_does(gis):
    with pytest.raises(ValueError, match="may not be empty"):
        _stack().excavate(93.0, _forms(gis)["area"]["within"], name="")


# -- hashes and digests --------------------------------------------------------------


def test_a_lazy_disv_hashes_the_footprint_and_the_bed(gis):
    held = mf.disv(_pond_stack(gis)).inputs["layer_stack"]
    [excavation] = held._excavations
    assert len(excavation.source_shas["zone"]) == 16
    assert len(excavation.bottom.source_shas["path"]) == 16
    entry = held.digest_payload()["excavations"][0]
    assert "zone" not in entry and "path" not in entry["bottom"], "the hash replaces the path"


def _digest(gis, site=None, bowl=None):
    stack = _stack().excavate(bottom=mf.Raster(bowl or gis.bowl), within=site or gis.site,
                              layer="ponds", name="pond")
    return package_digest(mf.disv(stack))


def test_moving_the_files_keeps_the_digest_and_editing_either_moves_it(tmp_path, gis):
    moved = tmp_path / "elsewhere"
    moved.mkdir()
    site, bowl = shutil.copy2(gis.site, moved / "site.gpkg"), shutil.copy2(gis.bowl, moved / "bowl.tif")
    before = _digest(gis)
    assert _digest(gis, site, bowl) == before, "moving is not editing"

    _write_pond(site, centre=(310.0, 300.0))
    assert _digest(gis, site, bowl) != before, "an edited footprint moves the digest"
    _raster(bowl, lambda x, y: 91.0 + 0 * x)
    assert _digest(gis, gis.site, bowl) != before, "an edited bed moves the digest"


def test_an_areas_note_and_source_are_not_content(gis):
    frame = gpd.read_file(gis.site, layer="ponds")

    def digest(area):
        return package_digest(mf.disv(_stack().excavate(93.0, area, name="pond")))

    plain = digest(mf.Area.coerce(frame))
    assert digest(mf.Area.coerce(gis.site, layer="ponds", note="drawn by hand")) == plain
    assert digest(mf.Area.coerce(frame.translate(xoff=5.0))) != plain


def test_an_excavated_lazy_disv_digests_the_same_across_json_and_a_save(tmp_path, gis):
    spec = mf.disv(_pond_stack(gis))
    reloaded = mf.PackageSpec.from_dict(json.loads(json.dumps(spec.to_dict())))
    assert package_digest(spec) == package_digest(reloaded)

    project = Project(tmp_path / "p", name="p")
    project.add_package("disv/pond", spec)
    project.save()
    assert package_digest(Project.load(tmp_path / "p").package("disv/pond")) == package_digest(spec)


def test_an_unexcavated_stack_writes_no_excavations_key():
    """So its payload -- and every stored lazy-DISV digest -- is what it was."""

    assert set(_stack().to_dict()) == {"top", "layers", "length_units", "time_units"}


def test_every_file_an_excavation_reads_is_found_by_the_walkers(tmp_path, gis):
    """The bundler, the anchor, ``built_from`` and the stamping share ONE walker
    (``specs.repoint_source_paths``); it finds the footprint under ``zone`` and a
    contour bed's files with no excavation-specific wiring. No GRASS run: only hashed."""

    contours = tmp_path / "bed_contours.gpkg"
    gpd.GeoDataFrame({"z": [92.0]}, geometry=[shp.LineString([(250, 300), (350, 300)])],
                     crs=2927).to_file(contours)
    bed = mf.Contours(contours, z="z", region_vector=gis.site, out=tmp_path / "bed.tif")
    spec = mf.disv(_stack().excavate(bed, gis.site, layer="ponds", name="pond"))
    found = {path.split("/")[-1] for path in source_paths(spec.inputs)}
    assert {"site.gpkg", "bed_contours.gpkg"} <= found and "bed.tif" not in found


# -- the archive, the snapshot, drift ------------------------------------------------


def _project(tmp_path, gis):
    project = Project(tmp_path / "proj", name="p")
    project.add_grid("g", _grid_spec(gis.boundary, 40_000.0))
    project.add_simulation(_simulation(gis, grid=mf.grid_ref("g")))
    return project


def _disv(run):
    disv = run.built.built_model("m").model.get_package("disv")
    return {name: np.asarray(getattr(disv, name).get_data()) for name in ("top", "botm", "idomain")}


def test_a_bundled_excavated_project_rebuilds_with_every_source_deleted(tmp_path, gis):
    project = _project(tmp_path, gis)
    before = _disv(project.prepare_run("before", "s"))
    assert (before["idomain"][0] == 0).any(), "the premise: the pond is dug"

    report = project.save(bundle_sources=True)
    carried = {str(source).split("/")[-1] for source, _ in report.bundled}
    assert {"site.gpkg", "bowl.tif", "dom.gpkg"} <= carried, sorted(carried)

    shutil.rmtree(gis.root)
    reopened = Project.load(tmp_path / "proj")
    after = _disv(reopened.prepare_run("after", "s"))
    for name in before:
        np.testing.assert_array_equal(after[name], before[name])


def test_a_run_rebuilds_from_its_snapshot_with_every_source_deleted(tmp_path, gis):
    project = _project(tmp_path, gis)
    run = project.prepare_run("r", "s")
    before = _disv(run)
    live = run.built.built_model("m").context
    shutil.rmtree(gis.root)

    thawed = thaw_simulation(run.workspace)
    rebuilt = thawed.build_flopy(workspace=tmp_path / "again")
    disv = rebuilt.built_model("m").model.get_package("disv")
    for name, values in before.items():
        np.testing.assert_array_equal(np.asarray(getattr(disv, name).get_data()), values)

    # The RECORDS come back with the elevations, where the build and `mf.lak` read them:
    # the snapshot kept the cut and dropped the names of what it cut.
    records = live.surfaces.attrs["excavations"]
    context = thawed.models[0].context
    assert context.surfaces.attrs["excavations"] == records
    assert context.grid.excavations == records
    assert excavation_cells(rebuilt.built_model("m").context) == excavation_cells(live)
    assert len(excavation_cells(live)) > 0


def test_a_package_on_a_dug_cell_of_a_reopened_run_is_refused_naming_the_excavation(tmp_path, gis):
    project = _project(tmp_path, gis)
    run = project.prepare_run("r", "s")
    idomain = _disv(run)["idomain"]
    record = run.built.built_model("m").context.surfaces.attrs["excavations"]["pond"]
    cell = next(c for c in record.cells if idomain[0, c] == 0)

    thawed = thaw_simulation(run.workspace)
    chd = mf.chd(stress_period_data={0: [[(0, cell), 98.0]]})
    with pytest.raises(ValueError, match=r"1 inside excavation 'pond'"):
        thawed.with_packages(chd, on="m").build_flopy(workspace=tmp_path / "again")


def test_drift_reports_an_edited_footprint(tmp_path, gis):
    project = _project(tmp_path, gis)
    project.prepare_run("r1", "s")
    rows = project.runs.drift().set_index("key")
    assert rows.loc["disv", "source"] == "ok"

    _write_pond(gis.site, centre=(320.0, 300.0))
    rows = project.refresh_runs().drift().set_index("key")
    assert rows.loc["disv", "source"].startswith("edited")


# -- meshes ----------------------------------------------------------------------------


@pytest.mark.parametrize("max_area, ncpl", [(40_000.0, 28), (12_000.0, 92)])
def test_one_excavated_declaration_builds_on_two_meshes(tmp_path, gis, max_area, ncpl):
    """Serialized, the stack re-reads the footprint and re-covers the cells on whatever
    mesh the build supplies -- the cells ``mf.lak`` takes on THAT mesh."""

    built = _simulation(gis, max_area).build_flopy(tmp_path / f"ws{ncpl}").built_model("m")
    grid = built.context.grid
    assert grid.ncpl == ncpl
    record = built.context.surfaces.attrs["excavations"]["pond"]
    assert record.cells == _lake_cells(grid, gis)
    disv = built.model.get_package("disv")
    idomain = np.asarray(disv.idomain.get_data())
    top = np.asarray(disv.top.get_data())
    botm = np.asarray(disv.botm.get_data())
    tops = np.vstack([top[None, :], botm[:-1]])
    for cell, bed in zip(record.cells, record.bed, strict=True):
        if bed < 100.0:
            k = int(np.flatnonzero(idomain[:, cell] > 0)[0])
            assert tops[k, cell] == bed and (idomain[:k, cell] == 0).all()


def test_with_grid_moves_an_excavated_model_28_to_92(tmp_path, gis):
    spec = _simulation(gis)
    moved = spec.with_grid(_grid_spec(gis.boundary, 12_000.0), name="s_fine")
    assert {row["package"]: row["how"] for row in moved.move_report()}["disv"] == "rebuilt"

    built = moved.build_flopy(tmp_path / "fine").built_model("m")
    assert built.context.grid.ncpl == 92
    record = built.context.surfaces.attrs["excavations"]["pond"]
    assert record.cells == _lake_cells(built.context.grid, gis)


# -- what the lazy DISV publishes ------------------------------------------------------


@pytest.mark.parametrize("attach", [True, False])
def test_the_lazy_disv_publishes_what_it_built_whatever_attach_says(tmp_path, gis, attach):
    grid = _grid_spec(gis.boundary, 40_000.0).resolve(workspace=tmp_path / "g")
    mf.LayerStack(top=mf.Flat(500.0)).add("x", thickness=10.0).build(grid, attach=True)  # stale
    drains = tmp_path / "drn.gpkg"
    gpd.GeoDataFrame({"cond": [10.0], "layer": [1], "name": ["d"]},
                     geometry=[shp.LineString([(20, 20), (20, 580)])], crs=2927).to_file(drains)
    built = _simulation(gis, stack=_pond_stack(gis), attach=attach,
                        ).with_grid(mf.GridSpec.from_object(grid))
    # A package that READS the elevations, lazily, so it reads the published table.
    built = built.with_packages(mf.drn.gpkg(
        drains, conductance="cond", elevation=mf.CellSurfaceOffset("model_top", -2.0), nper=1,
    ), on="m")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = built.build_flopy(tmp_path / "ws").built_model("m")
    # The published table IS this model's DISV, so the grid's older one (kept by
    # attach=False) disagreeing with it is not a stale CONTEXT -- the warning that said
    # so advised rebuilding the context, which would put 500.0 back.
    assert not [w for w in caught if "disagree" in str(w.message)]
    drn = model.model.get_package("drn").stress_period_data.get_data()[0]
    top = np.asarray(model.model.get_package("disv").top.get_data())
    np.testing.assert_allclose(drn["elev"], top[[cellid[1] for cellid in drn["cellid"]]] - 2.0)

    table = model.context.surfaces
    np.testing.assert_array_equal(np.asarray(table[0], dtype=float),
                                  np.asarray(model.model.get_package("disv").top.get_data()))
    assert list(table.columns)[0] == "geometry", "geometry first (ledger 196)"
    assert excavation_cells(model.context) == {c: "pond" for c in table.attrs["excavations"]["pond"].cells}
    stale_top = float(np.asarray(model.context.grid.gdf_topbtm[0])[0])
    assert (stale_top == 500.0) is (not attach), "attach decides only the GRID"


def test_a_published_record_is_shared_not_copied_by_every_derived_frame(tmp_path, gis):
    """pandas deep-copies ``attrs`` into every frame derived from the table; a record is
    immutable (frozen, tuple cells, read-only bed), so it copies as itself."""

    import copy

    built = _simulation(gis).build_flopy(tmp_path / "ws").built_model("m")
    table = built.context.surfaces
    record = table.attrs["excavations"]["pond"]
    assert copy.deepcopy(record) is record and copy.copy(record) is record
    assert not record.bed.flags.writeable
    assert table[[0, 1]].attrs["excavations"]["pond"] is record


def test_a_stale_package_on_a_lazy_model_is_refused_naming_the_excavation(tmp_path, gis):
    """The C5 guard reads the record the lazy DISV published on the context's table."""

    grid = _grid_spec(gis.boundary, 40_000.0).resolve(workspace=tmp_path / "g")
    dug = _pond_stack(gis).build(grid)
    cell = next(c for c in dug.excavations["pond"].cells if dug.idomain[0, c] == 0)
    model = mf.gwf("m", context=ModelContext(grid=grid), packages=[
        mf.disv(_pond_stack(gis), attach=False), mf.npf(k=1.0), mf.ic(strt=99.0),
        mf.chd(stress_period_data={0: [[(0, cell), 98.0]]}),
    ])
    spec = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    grid.excavations = {}
    with pytest.raises(ValueError, match=r"1 inside excavation 'pond'"):
        spec.build_flopy(tmp_path / "ws")


def test_an_excavation_bed_does_not_set_the_draft_grid_extent(gis):
    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
             .excavate(bottom=mf.Raster(gis.bowl), within=gis.site, layer="ponds"))
    assert stack._georeferenced_surfaces() == []
    with pytest.raises(ValueError, match="needs an extent"):
        stack.draft_grid()
