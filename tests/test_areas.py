"""Areas of interest: `model.aoi` / `project.aoi` and their names (ledger 310).

The user's design, each pinned here: areas are DEFINED ON THE MODEL and propagate
to the project when the model belongs to one (`model.aoi` IS the live
`project.aoi`); a model with no live project keeps them for the session only --
only `project.save()` writes a project's files; and a name then works wherever a
place does: `extent=`, `select=`, `features=`.
"""

from __future__ import annotations

import gc
import json
import shutil

import geopandas as gpd
import numpy as np
import pytest
import shapely as shp

import myflopy as mf
from myflopy import areas as areas_module
from myflopy.areas import Area
from myflopy.library import AoiLibrary

# --- a tiny DISV project that runs -----------------------------------------------


def _simulation(vor, name):
    from myflopy.specs import ModelContext

    props = dict(vor.get_disv_gridprops())
    model = mf.gwf(
        "flow",
        context=ModelContext(grid=vor, domain=np.ones((1, vor.ncpl), dtype=int)),
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **props),
            mf.npf(k=1.0),
            mf.ic(strt=5.0),
            mf.chd(stress_period_data=[((0, 0), 6.0), ((0, vor.ncpl - 1), 4.0)]),
            mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
        ],
    )
    return mf.simulation(model, name=name, packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _build(root):
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi

    vor = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=3, nlay=1, nper=1))
    project = mf.Project(root, name="p")
    run = project.prepare_run("base", _simulation(vor, "base"))
    assert run.execute()[0]
    return project


@pytest.fixture
def project(tmp_path):
    return _build(tmp_path / "proj")


def _saved_and_gone(root, **areas):
    """Build, add ``areas``, save -- and leave NO live Project behind.

    A fixture would keep the Project alive for the whole test, and a model would
    then find it in the live registry: the "no live project" tests passed for the
    wrong reason until this helper replaced the fixture there.
    """

    project = _build(root)
    for name, place in areas.items():
        project.add_aoi(name, place)
    project.save()
    workspace = project.runs["base"].workspace
    del project
    gc.collect()
    assert areas_module.live_areas(root) is None, "the project is still alive"
    return workspace


SITE = (0.0, 0.0, 150.0, 150.0)       # the south-west quarter of a 0..300 grid
NORTH = (0.0, 200.0, 300.0, 300.0)    # the top row of cells


# --- the value ---------------------------------------------------------------------


@pytest.mark.parametrize(
    ("spelling", "bounds"),
    [
        (SITE, SITE),
        ([(0, 0), (150, 150), (40, 90)], SITE),
        (shp.LineString([(0, 150), (150, 0)]), SITE),
    ],
)
def test_an_area_is_built_from_any_place_spelling(spelling, bounds):
    area = Area.coerce(spelling)
    assert area.geometry.bounds == pytest.approx(bounds)
    assert area.crs is None  # the grid's frame -- never a default filled in


def test_a_frame_keeps_its_own_crs_and_refuses_a_relabel():
    frame = gpd.GeoDataFrame(geometry=[shp.box(*SITE)], crs="EPSG:2927")
    assert Area.coerce(frame).crs == "EPSG:2927"
    with pytest.raises(ValueError, match="relabelling"):
        Area.coerce(frame, crs="EPSG:4326")


def test_an_area_reads_a_vector_file_once_and_keeps_where_it_came_from(tmp_path):
    path = tmp_path / "site.gpkg"
    gpd.GeoDataFrame(geometry=[shp.box(*SITE)], crs="EPSG:2927").to_file(path, layer="aoi")
    area = Area.coerce(path, layer="aoi")
    assert area.geometry.bounds == pytest.approx(SITE)
    assert area.source == str(path)


@pytest.mark.parametrize(
    ("value", "error"),
    [((150, 0, 0, 150), ValueError), (object(), TypeError)],
)
def test_what_is_not_a_place_raises(value, error):
    with pytest.raises(error):
        Area.coerce(value)


def test_the_saved_document_is_plain_json_at_full_precision():
    """WKT defaults to 6 decimals in shapely; a lon/lat vertex would move ~0.1 m."""

    area = Area.coerce((-122.123456789012, 47.0, -122.0, 47.5), crs="EPSG:4326", note="n")
    doc = area.to_dict()
    assert all(type(value) is str for value in doc.values())
    assert json.loads(json.dumps(doc)) == doc
    back = Area.from_dict(doc)
    assert back.geometry.equals_exact(area.geometry, tolerance=0.0)
    assert (back.crs, back.note) == ("EPSG:4326", "n")


# --- the library -------------------------------------------------------------------


def test_assignment_coerces_so_the_dict_and_the_verb_agree():
    library = AoiLibrary()
    library["site"] = SITE
    library.update({"north": NORTH})
    library.setdefault("line", shp.LineString([(0, 0), (10, 5)]))
    assert all(isinstance(area, Area) for area in library.values())
    assert library["site"].name == "site"


def test_registering_an_existing_name_refuses_unless_asked():
    library = AoiLibrary()
    library._define("site", SITE)
    with pytest.raises(KeyError, match="overwrite=True"):
        library._define("site", NORTH)
    library._define("site", NORTH, overwrite=True)
    assert library["site"].geometry.bounds == pytest.approx(NORTH)


def test_a_missed_name_suggests_the_near_one():
    library = AoiLibrary()
    library["site"] = SITE
    with pytest.raises(KeyError, match="Did you mean: 'site'"):
        library["sitte"]


# --- a project keeps them ----------------------------------------------------------


def test_a_project_saves_its_areas_inline_and_loads_them_back(project):
    project.add_aoi("site", SITE, note="SW quarter")
    report = project.save()
    assert report.aoi == 1 and bool(report)
    document = json.loads(project.layout.project_spec_path.read_text())
    assert document["aoi"]["site"]["kind"] == "Area"
    # Inline: no per-area file, so an archive needs nothing extra.
    assert not list(project.layout.specs_dir.glob("**/*site*"))

    loaded = mf.Project.load(project.root)
    assert loaded.aoi["site"].geometry.bounds == pytest.approx(SITE)
    assert loaded.aoi["site"].note == "SW quarter"
    assert "aoi=1" in repr(loaded)
    assert "area" in set(loaded.summary()["library"])


def test_a_project_with_only_areas_saves_something(tmp_path):
    empty = mf.Project(tmp_path / "areas_only")
    empty.add_aoi("site", SITE)
    assert bool(empty.save())


def test_a_project_saved_before_areas_existed_still_loads(project):
    project.save()
    spec = project.layout.project_spec_path
    document = json.loads(spec.read_text())
    document.pop("aoi")
    spec.write_text(json.dumps(document))
    assert len(mf.Project.load(project.root).aoi) == 0


# --- a model shares them -----------------------------------------------------------


def test_a_models_areas_are_its_live_projects(project):
    project.add_aoi("site", SITE)
    model = project.runs["base"].model("flow")  # a SEPARATE Run, rebuilt from disk
    assert model.aoi is project.aoi
    model.add_aoi("north", NORTH)
    assert sorted(project.aoi) == ["north", "site"]


def test_a_prepared_runs_model_shares_them_too(project, tmp_path):
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi

    vor = rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=3, nlay=1, nper=1))
    run = project.prepare_run("second", _simulation(vor, "second"))
    assert run.model("flow").aoi is project.aoi


def test_without_a_live_project_the_saved_areas_are_read_and_additions_stay_in_session(tmp_path):
    root = tmp_path / "proj"
    workspace = _saved_and_gone(root, site=SITE)

    model = mf.load_run(workspace).model("flow")
    assert sorted(model.aoi) == ["site"]
    model.add_aoi("scratch", NORTH)
    assert "scratch" in model.aoi
    assert sorted(mf.Project.load(root).aoi) == ["site"]  # nothing was written


def test_a_moved_project_is_still_found_by_its_run_layout(tmp_path):
    """A copied project's runs still record the ORIGINAL root."""

    root = tmp_path / "proj"
    _saved_and_gone(root, site=SITE)
    copy = tmp_path / "copy"
    shutil.copytree(root, copy)
    shutil.rmtree(root)
    model = mf.load_run(copy / "runs" / "base").model("flow")
    assert model._project_root is not None and not model._project_root.exists()
    assert sorted(model.aoi) == ["site"]


def test_a_model_outside_any_project_keeps_areas_for_the_session(canonical_run):
    canonical_run.add_aoi("canonical_site", (500, 600, 1000, 1100))
    assert "canonical_site" in canonical_run.aoi
    del canonical_run.aoi["canonical_site"]  # the session model is shared


# --- names work wherever a place does ----------------------------------------------


def test_a_name_frames_a_map_selects_cells_and_draws_as_a_feature(project):
    project.add_aoi("site", SITE)
    project.add_aoi("north", NORTH)
    model = project.runs["base"].model("flow")

    static = model.plot.map(layer=0, extent="site", bounds_padding=0.0, backend="mpl")
    axis = static.axes[0]
    assert (*axis.get_xlim(), *axis.get_ylim()) == pytest.approx((0.0, 150.0, 0.0, 150.0))

    picture = model.plot.map(layer=0, select="north")
    cells = sorted(int(c) for c in model.vor.get_vor_cells_as_series(shp.box(*NORTH)).explode())
    assert picture._select_cells == sorted(set(cells))

    assert len(model.plot.map(layer=0, features=model.aoi).overlay_traces()) == 2

    figure = model.hds.mosaic(by="layer", extent="north")
    assert figure.layout.map.zoom > model.plot.map(layer=0).fig.layout.map.zoom


def test_a_name_that_is_both_a_region_and_an_area_is_refused(project):
    project.add_aoi("site", SITE)
    model = project.runs["base"].model("flow")
    model.add_region("site", cellids=[0])
    with pytest.raises(ValueError, match="BOTH a region"):
        model.plot.map(layer=0, select="site")


def test_an_unknown_name_lists_what_is_defined(project):
    project.add_aoi("site", SITE)
    model = project.runs["base"].model("flow")
    with pytest.raises(ValueError, match=r"areas: \['site'\]"):
        model.plot.map(layer=0, select="sitte")
    with pytest.raises(KeyError, match="site"):
        model.plot.map(layer=0, extent="sitte")


def test_a_map_without_a_named_extent_never_reads_the_saved_areas(tmp_path, monkeypatch):
    """The lookup can open the project's spec; only a NAME may pay for it."""

    workspace = _saved_and_gone(tmp_path / "proj")
    model = mf.load_run(workspace).model("flow")

    def explode(root):
        raise AssertionError("read the saved areas for a map that named none")

    monkeypatch.setattr(areas_module, "saved_areas", explode)
    import myflopy.modflow.mf6.simulation.base as base

    monkeypatch.setattr(base, "saved_areas", explode)
    model.plot.map(layer=0, extent=SITE)
    model.plot.map(layer=0, select=[0])


# --- review fixes (ledger 312) -----------------------------------------------------


def test_every_dict_spelling_coerces(project):
    """`|=`, a reflected `|` and the constructor inserted raw tuples, which broke
    `save()` half-way and made `select=name` read the box as cell ids."""

    project.aoi |= {"sw": SITE}
    merged = {"n": NORTH} | project.aoi
    built = AoiLibrary({"x": SITE})
    for library in (project.aoi, merged, built):
        assert all(isinstance(area, Area) for area in library.values())
    assert bool(project.save())


def test_a_save_that_cannot_build_its_document_writes_nothing(project, monkeypatch):
    project.add_aoi("site", SITE)

    def broken(self):
        raise RuntimeError("cannot serialize")

    monkeypatch.setattr(Area, "to_dict", broken)
    with pytest.raises(RuntimeError):
        project.save()
    assert not project.manifest_path.exists()  # no manifest naming a missing spec


def test_a_copied_projects_models_use_the_copys_areas(tmp_path):
    """A copy's runs still record the ORIGINAL root; trusting that gave the copy's
    models the original's areas, and add_aoi wrote into the original."""

    original = tmp_path / "orig"
    _saved_and_gone(original, site=SITE)
    copy = tmp_path / "copy"
    shutil.copytree(original, copy)
    variant = mf.Project.load(copy)
    variant.add_aoi("variant_site", NORTH)
    variant.save()
    del variant
    gc.collect()

    model = mf.load_run(copy / "runs" / "base").model("flow")  # no live project
    assert sorted(model.aoi) == ["site", "variant_site"]

    live_original = mf.Project.load(original)
    live_copy = mf.Project.load(copy)
    model = live_copy.runs["base"].model("flow")
    assert model.aoi is live_copy.aoi
    model.add_aoi("only_the_copy", SITE)
    assert "only_the_copy" not in live_original.aoi


def test_a_short_lived_second_project_does_not_unregister_the_first(project):
    project.add_aoi("site", SITE)
    project.save()
    mf.Project.load(project.root).runs.summary()  # a notebook cell, then collected
    gc.collect()
    model = project.runs["base"].model("flow")
    assert model.aoi is project.aoi


def test_rebinding_the_library_keeps_the_models_on_it(project):
    project.add_aoi("site", SITE)
    project.aoi = project.aoi | {"north": NORTH}
    gc.collect()
    model = project.runs["base"].model("flow")
    assert model.aoi is project.aoi and sorted(model.aoi) == ["north", "site"]


def test_a_run_found_anywhere_under_the_project_shares_its_areas(project):
    """Discovered runs record no root, and the old guess looked exactly two levels
    up -- `runs/scenarios/x` or `imports/usg/converted` never found the project."""

    project.add_aoi("site", SITE)  # never saved: only the LIVE project knows it
    workspace = project.runs["base"].workspace
    assert mf.load_mf6_run(workspace).aoi is project.aoi
    nested = project.root / "runs" / "scenarios" / "copy"
    shutil.copytree(workspace, nested)
    assert mf.load_mf6_run(nested).aoi is project.aoi


def test_an_area_can_be_styled_as_a_feature(project):
    project.add_aoi("site", SITE)
    model = project.runs["base"].model("flow")
    styled = mf.plot.Feature(model.aoi["site"], color="red", fill=True)
    picture = model.plot.map(layer=0, features=styled)
    assert picture.overlay_traces()
    assert model.plot.map(layer=0, features=styled, backend="mpl") is not None


def test_the_same_crs_spelled_differently_is_not_a_relabel():
    area = Area.coerce(SITE, crs="EPSG:26910")
    assert Area.coerce(area, crs=26910).crs == "EPSG:26910"
    with pytest.raises(ValueError, match="relabel"):
        Area.coerce(area, crs=4326)
