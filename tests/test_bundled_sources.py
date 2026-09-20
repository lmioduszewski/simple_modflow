"""A saved project can carry the files its grids are built from.

`ModelSpec.build` re-resolves a `GridSpec` on EVERY build, so a project whose
boundary GeoPackage moved is a project that no longer builds -- and the failure is a
raw `pyogrio.errors.DataSourceError` with no myflopy frame in it.
`save(bundle_sources=True)` copies those files into `inputs/` and repoints the
recipes, which also makes relative paths the normal case: hence the anchoring half.

LAZY package sources are bundled too (ledger 293) -- a `.gpkg` declared without a
`context=`, and the layer stack behind a lazy `mf.disv`. An EAGER package still names
nothing: it resolves at declaration and holds finished records, measured as zero file
reads inside `PackageSpec.build` across the suite, so copying its file would archive
something nothing opens.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.project.bundle import bundle_grid_sources, grid_sources
from myflopy.project.freeze import discard_snapshot
from myflopy.workspace import Project, Run


def _gpkg(path, size=600.0):
    """A one-polygon GeoPackage usable as a grid boundary."""

    path.parent.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame(
        {"name": ["d"]}, geometry=[box(0, 0, size, size)], crs=2927
    ).to_file(path)
    return path


def _grid_spec(boundary, **kwargs):
    return mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=str(boundary)), max_area=40000.0, **kwargs
    )


def _deferred_sim(boundary):
    model = mf.gwf(
        "m", packages=[mf.npf(k=1.0), mf.ic(strt=5.0), mf.oc(head_filerecord="m.hds")]
    ).with_grid(_grid_spec(boundary))
    return mf.simulation(
        model, name="rel", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]
    )


# ---------------------------------------------------------------- what it copies


def test_a_bundled_project_resolves_after_its_source_is_deleted(tmp_path):
    """The whole point: the project stops depending on where it was authored."""

    outside = _gpkg(tmp_path / "elsewhere" / "domain.gpkg")
    root = tmp_path / "proj"
    project = Project(root, name="p")
    project.add_grid("coarse", _grid_spec(outside))

    report = project.save(bundle_sources=True)
    assert [to for _from, to in report.bundled] == ["inputs/domain.gpkg"]

    shutil.rmtree(tmp_path / "elsewhere")
    loaded = Project.load(root)
    grid = loaded.grids["coarse"].resolve(project_root=root, workspace=root / "_g")
    assert grid.ncpl > 0


def test_the_stored_path_becomes_project_relative(tmp_path):
    """An absolute path pins a project to one machine; a relative one travels."""

    outside = _gpkg(tmp_path / "elsewhere" / "domain.gpkg")
    project = Project(tmp_path / "proj", name="p")
    project.add_grid("coarse", _grid_spec(outside))
    project.save(bundle_sources=True)

    boundary = project.grids["coarse"].boundary
    assert str(boundary.path) == "inputs/domain.gpkg"
    assert boundary.external is False        # it lives inside the project now


def test_saving_without_the_flag_copies_nothing(tmp_path):
    outside = _gpkg(tmp_path / "elsewhere" / "domain.gpkg")
    project = Project(tmp_path / "proj", name="p")
    project.add_grid("coarse", _grid_spec(outside))

    report = project.save()
    assert report.bundled == ()
    assert not list((tmp_path / "proj" / "inputs").glob("*.gpkg"))


def test_an_inline_grid_on_a_model_is_bundled_too(tmp_path):
    """Both tiers: the grid library, and a GridSpec written straight onto a model."""

    outside = _gpkg(tmp_path / "elsewhere" / "domain.gpkg")
    project = Project(tmp_path / "proj", name="p")
    project.add_simulation(_deferred_sim(outside))
    project.save(bundle_sources=True)

    stored = project.simulations["rel"].models[0].grid.boundary.path
    assert str(stored) == "inputs/domain.gpkg"


def test_the_registry_is_repointed_so_a_later_save_cannot_undo_it(tmp_path):
    """A plain save() after a bundle must not write the original paths back."""

    outside = _gpkg(tmp_path / "elsewhere" / "domain.gpkg")
    root = tmp_path / "proj"
    project = Project(root, name="p")
    project.add_grid("coarse", _grid_spec(outside))
    project.save(bundle_sources=True)
    project.save()                                  # no flag this time

    shutil.rmtree(tmp_path / "elsewhere")
    assert str(Project.load(root).grids["coarse"].boundary.path) == "inputs/domain.gpkg"


def test_two_sources_with_one_basename_do_not_collide(tmp_path):
    """Otherwise the second silently overwrites the first, and both read the survivor."""

    first = _gpkg(tmp_path / "a" / "domain.gpkg", size=600.0)
    second = _gpkg(tmp_path / "b" / "domain.gpkg", size=900.0)
    project = Project(tmp_path / "proj", name="p")
    project.add_grid("one", _grid_spec(first))
    project.add_grid("two", _grid_spec(second))

    project.save(bundle_sources=True)
    stored = {
        str(project.grids["one"].boundary.path),
        str(project.grids["two"].boundary.path),
    }
    assert len(stored) == 2
    assert len(list((tmp_path / "proj" / "inputs").glob("*.gpkg"))) == 2


def test_a_source_already_inside_inputs_is_not_copied_again(tmp_path):
    root = tmp_path / "proj"
    inside = _gpkg(root / "inputs" / "domain.gpkg")
    project = Project(root, name="p")
    project.add_grid("coarse", _grid_spec(inside))

    report = project.save(bundle_sources=True)
    assert report.bundled == ()
    assert inside.exists()


def test_a_missing_source_is_reported_not_fatal(tmp_path, caplog):
    """One moved file must not cost the project everything else it would have saved."""

    import logging

    gone = tmp_path / "elsewhere" / "domain.gpkg"
    gone.parent.mkdir(parents=True)
    project = Project(tmp_path / "proj", name="p")
    project.add_grid("coarse", _grid_spec(gone))

    with caplog.at_level(logging.WARNING, logger="myflopy.project.bundle"):
        report = project.save(bundle_sources=True)
    assert report.bundled == ()
    assert "could not be found" in caplog.text


def test_a_directory_source_is_copied_whole(tmp_path):
    """A .gdb is a DIRECTORY; copying only its "file" bundles nothing."""

    from myflopy.specs import GridSpec

    gdb = tmp_path / "elsewhere" / "data.gdb"
    gdb.mkdir(parents=True)
    (gdb / "a00000001.gdbtable").write_bytes(b"x")
    spec = GridSpec(
        name="g", grid_type="disv", method="voronoi",
        boundary=mf.GeoPackageSourceSpec(path=str(gdb)),
    )
    bundled, report = bundle_grid_sources(
        spec, root=tmp_path / "proj", inputs_dir=tmp_path / "proj" / "inputs"
    )
    assert str(bundled.boundary.path) == "inputs/data.gdb"
    assert (tmp_path / "proj" / "inputs" / "data.gdb" / "a00000001.gdbtable").exists()
    assert len(report["copied"]) == 1


def test_every_named_path_is_reported(tmp_path):
    """``grid_sources`` is the list bundling acts on, so it must be complete."""

    from myflopy.specs import GridSpec

    raster = tmp_path / "top.tif"
    idomain = tmp_path / "idom.gpkg"
    spec = GridSpec(
        name="g", grid_type="disv", method="voronoi",
        boundary=mf.GeoPackageSourceSpec(path=str(tmp_path / "b.gpkg")),
        refinement=mf.GeoPackageSourceSpec(path=str(tmp_path / "r.gpkg")),
        breaklines=(mf.GeoPackageSourceSpec(path=str(tmp_path / "bl.gpkg")),),
        points=(mf.GeoPackageSourceSpec(path=str(tmp_path / "pt.gpkg")),),
        voronoi_options={"rasters": [str(raster)], "idomain_path": str(idomain)},
    )
    named = {path.name for path in grid_sources(spec)}
    assert named == {"b.gpkg", "r.gpkg", "bl.gpkg", "pt.gpkg", "top.tif", "idom.gpkg"}


def test_the_lazy_option_paths_are_bundled(tmp_path):
    """``rasters`` and ``idomain_path`` are read long after the build, not during it.

    Leaving them out means a bundled project builds and then dies on its first
    cross-section or edge query, which is not archived.
    """

    from myflopy.specs import GridSpec

    raster = tmp_path / "elsewhere" / "top.tif"
    raster.parent.mkdir(parents=True)
    raster.write_bytes(b"not really a tif")
    spec = GridSpec(
        name="g", grid_type="disv", method="voronoi",
        boundary=mf.GeoPackageSourceSpec(path=str(_gpkg(tmp_path / "elsewhere" / "b.gpkg"))),
        voronoi_options={"rasters": [str(raster)]},
    )
    bundled, _report = bundle_grid_sources(
        spec, root=tmp_path / "proj", inputs_dir=tmp_path / "proj" / "inputs"
    )
    assert [str(p) for p in bundled.voronoi_options["rasters"]] == ["inputs/top.tif"]


# ---------------------------------------------------------------- the anchor


def test_a_relative_source_anchors_to_the_project_not_the_directory_you_are_in(tmp_path):
    """Bundling writes relative paths, so this is what makes a bundle portable.

    ``Project.reopen_run`` sets a build context with the project root;
    ``Run.reopen`` and ``load_run`` had none, so the resolver fell back to the
    current working directory and the same run failed from a different folder.
    """

    root = tmp_path / "proj"
    _gpkg(root / "inputs" / "domain.gpkg")
    project = Project(root, name="p")
    project.add_simulation(_deferred_sim("inputs/domain.gpkg"))
    run = project.prepare_run("r1", "rel")
    run.built.simulation.write_simulation(silent=True)
    discard_snapshot(run.workspace)          # force the RECIPE path

    assert run.project_root == root

    here = os.getcwd()
    os.chdir(tmp_path)                       # anywhere but the project
    try:
        reopened = Run.reopen(root / "runs" / "r1")
        assert reopened.project_root == root
        reopened.build()
    finally:
        os.chdir(here)


def test_a_manifest_without_a_project_root_still_opens(tmp_path):
    """Runs prepared before this was recorded must keep working."""

    import json

    root = tmp_path / "proj"
    _gpkg(root / "inputs" / "domain.gpkg")
    project = Project(root, name="p")
    project.add_simulation(_deferred_sim("inputs/domain.gpkg"))
    run = project.prepare_run("r1", "rel")
    run.built.simulation.write_simulation(silent=True)

    manifest = json.loads(run.manifest_path.read_text())
    manifest.pop("project_root")
    run.manifest_path.write_text(json.dumps(manifest))

    reopened = Run.reopen(root / "runs" / "r1")
    assert reopened.project_root is None
    reopened.build()                         # the snapshot still carries it through


# --- three defects in the shipped bundler (ledger 292) ------------------------


def _shapefile(path):
    """A one-polygon SHAPEFILE, which is five files sharing a stem."""

    path.parent.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame(
        {"name": ["d"]}, geometry=[box(0, 0, 600.0, 600.0)], crs=2927
    ).to_file(path)
    return path


def test_a_multi_file_source_is_bundled_WHOLE(tmp_path):
    """Copying only the file the recipe names produced a project that reported success
    and could no longer be OPENED.

    Measured before this, on a tree where NOTHING was deleted: a shapefile written with
    five sidecars (.cpg .dbf .prj .shp .shx) was bundled as the `.shp` alone, the grid
    that resolved before the save then failed `DataSourceError`, and `SaveReport` said it
    had bundled. That is strictly worse than not bundling -- it takes a project that
    builds and makes it not build (ledger 292). A `.gdb` directory was already carried
    whole by `copytree`; this is the other shape of the same mistake.
    """

    source = _shapefile(tmp_path / "gis" / "domain.shp")
    written = {path.suffix for path in source.parent.iterdir()}
    assert {".shp", ".shx", ".dbf"} <= written, "the fixture must really be multi-file"

    root = tmp_path / "proj"
    project = Project(root, name="p")
    project.add_grid("g", _grid_spec(source))
    project.save(bundle_sources=True)

    bundled = {path.suffix for path in (root / "inputs").iterdir()}
    assert written <= bundled, f"sidecars left behind: {sorted(written - bundled)}"
    # The real assertion: GDAL can open the copy. A .shp without its .shx cannot.
    assert len(gpd.read_file(root / "inputs" / "domain.shp")) == 1


def test_a_single_file_source_does_not_drag_in_its_neighbours(tmp_path):
    """The refusing direction. A `.gpkg` is self-contained, so a same-stem neighbour is
    somebody else's file and must stay where it is."""

    source = _gpkg(tmp_path / "gis" / "domain.gpkg")
    (source.parent / "domain.notes.txt").write_text("not mine")

    root = tmp_path / "proj"
    project = Project(root, name="p")
    project.add_grid("g", _grid_spec(source))
    project.save(bundle_sources=True)

    assert sorted(p.name for p in (root / "inputs").iterdir()) == ["domain.gpkg"]


def test_a_source_the_project_ALREADY_holds_is_not_copied_into_itself(tmp_path):
    """`already_inside` tested `inputs/` only, never the project root.

    A `method="object"` grid becomes a pickle on save, and the bundler then copied the
    project's own `specs/grids/g.pkl` into `inputs/` and repointed the recipe at the
    copy, orphaning the original. On a real grid that pickle is the largest thing a
    saved project holds, so it doubled it (ledger 292). It is also the shape
    `usg.export_gis(dir)` encourages -- writing sources INTO the project.
    """

    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi

    root = tmp_path / "proj"
    project = Project(root, name="p")
    project.add_grid("g", rectangular_voronoi(CanonicalModelConfig(nrow=3, ncol=3, nlay=1, nper=1)))

    # SAVE, LOAD, then bundle -- the live grid is still an OBJECT and names no pickle, so
    # a single save cannot reach this. Only the RELOADED recipe carries
    # `pickle_path="specs/grids/g.pkl"`, which is the path the bundler then copied into
    # `inputs/`. Writing this as one save passed with the fix reverted.
    project.save()
    reloaded = Project.load(root)
    assert str(reloaded.grid("g").pickle_path) == "specs/grids/g.pkl", "the premise"
    reloaded.save(bundle_sources=True)

    pickles = sorted(path.relative_to(root).as_posix() for path in root.rglob("*.pkl"))
    assert pickles == ["specs/grids/g.pkl"], f"duplicated: {pickles}"


def test_a_REFUSED_save_copies_nothing_and_moves_nothing(tmp_path):
    """Bundling was the first statement of `save`, ahead of `validate()`.

    So a `strict=True` refusal left `inputs/` populated and the LIVE registry repointed
    into it, with no spec file written -- a state nothing on disk records, and one a
    later plain `save()` would then persist without the user ever asking for an archive.
    """

    source = _gpkg(tmp_path / "gis" / "domain.gpkg")
    root = tmp_path / "proj"
    project = Project(root, name="p")
    project.add_grid("g", _grid_spec(source))
    project.add_simulation(
        mf.simulation(
            mf.gwf("m", packages=[mf.ref("nope")]), name="s", packages=[mf.tdis(nper=1)]
        )
    )

    with pytest.raises(ValueError, match="cannot be saved"):
        project.save(strict=True, bundle_sources=True)

    assert not (root / "inputs").exists()
    assert str(project.grid("g").boundary.path) == str(source)


# ----------------------------------------------- lazy packages: the second carrier


def _bc_gpkg(path, head=9.0):
    """A one-point GeoPackage usable as a boundary condition source."""

    from shapely.geometry import Point

    path.parent.mkdir(parents=True, exist_ok=True)
    gpd.GeoDataFrame(
        {"name": ["w"], "head": [head], "cond": [1.0], "layer": [1]},
        geometry=[Point(300.0, 300.0)], crs=2927,
    ).to_file(path, driver="GPKG")
    return path


def _raster(path, value=100.0):
    """A small GeoTIFF covering the `_gpkg` boundary, usable as a stack top."""

    import rasterio
    from rasterio.transform import from_origin

    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path, "w", driver="GTiff", height=20, width=20, count=1,
        dtype="float32", crs="EPSG:2927", transform=from_origin(0, 600, 30, 30),
    ) as handle:
        handle.write(np.full((20, 20), value, dtype="float32"), 1)
    return path


def _lazy_project(tmp_path):
    """A project carrying BOTH lazy carriers: a stack raster and a `.gpkg` recipe."""

    gis = tmp_path / "gis"
    boundary, bcs, top = _gpkg(gis / "b.gpkg"), _bc_gpkg(gis / "bc.gpkg"), _raster(gis / "top.tif")

    project = Project(tmp_path / "proj", name="p")
    project.add_grid("g", _grid_spec(boundary))
    stack = mf.LayerStack(top=mf.Surface.raster(top)).add("sand", thickness=50)
    model = mf.gwf(
        "m",
        grid=mf.grid_ref("g"),
        packages=[
            mf.disv(stack),
            mf.npf(k=1.0),
            mf.ic(strt=95.0),
            mf.oc(head_filerecord="m.hds"),
            mf.ghb.gpkg(bcs, nper=1, head="head", conductance="cond"),
        ],
    )
    project.add_simulation(
        mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    )
    return project, gis


def test_a_lazy_model_rebuilds_after_EVERY_source_it_reads_is_deleted(tmp_path):
    """The headline, and it needs both carriers to be true.

    A lazy `mf.disv` holds a gridless `LayerStack`, which per ledger 270 is the only
    shape a layered model saves in -- and its raster is opened inside the
    discretization, BEFORE any boundary condition is built. So this same assertion
    passes vacuously on a bundler that carries only the `.gpkg` recipes, right up
    until a model has a layer stack.
    """

    project, gis = _lazy_project(tmp_path)
    report = project.save(bundle_sources=True)
    carried = {os.path.basename(source) for source, _ in report.bundled}

    assert {"top.tif", "bc.gpkg"} <= carried, f"a carrier was missed: {sorted(carried)}"

    shutil.rmtree(gis)
    reopened = Project.load(tmp_path / "proj")
    gwf = reopened.prepare_run("r", reopened.simulations["s"]).model("m").gwf

    assert float(np.mean(gwf.get_package("disv").top.array)) == pytest.approx(100.0)
    assert len(gwf.get_package("ghb").stress_period_data.get_data(0)) == 1


def test_the_RASTER_is_what_fails_first_when_a_carrier_is_missed(tmp_path):
    """Pins the ORDER, which is the whole argument for doing both at once.

    Measured by deleting only the stack's raster and leaving the GeoPackage: the
    failure is `RasterioIOError` from inside the discretization, never a GeoPackage
    error -- so a bundler that carried only the boundary recipes would report a
    complete archive and produce one that cannot open.
    """

    project, gis = _lazy_project(tmp_path)
    project.save(bundle_sources=True)

    # Un-bundle the raster only: put the recipe back on a file that will not be there.
    (tmp_path / "proj" / "inputs" / "top.tif").unlink()
    shutil.rmtree(gis)
    reopened = Project.load(tmp_path / "proj")

    with pytest.raises(Exception) as caught:
        reopened.prepare_run("r", reopened.simulations["s"])

    assert "top.tif" in str(caught.value), (
        f"the raster must be what names itself first, got {caught.value!r}"
    )


def test_a_bundled_package_source_is_read_from_the_project_not_the_cwd(tmp_path, monkeypatch):
    """The anchoring half, on the PACKAGE tier.

    Bundling stores project-relative paths; nothing resolved them for a package, so a
    bundled project rebuilt from any directory but its own root opened
    `inputs/bc.gpkg` against `os.getcwd()`. Grids never had the bug --
    `grid_spec_resolver._source_path` has anchored them all along -- which is exactly
    why the package tier went unnoticed.
    """

    project, gis = _lazy_project(tmp_path)
    project.save(bundle_sources=True)
    shutil.rmtree(gis)

    monkeypatch.chdir(tmp_path.parent)

    reopened = Project.load(tmp_path / "proj")
    gwf = reopened.prepare_run("r", reopened.simulations["s"]).model("m").gwf

    assert len(gwf.get_package("ghb").stress_period_data.get_data(0)) == 1


def test_an_EAGER_package_names_no_file_to_bundle(tmp_path):
    """The negative, and it is what keeps the walk cheap and safe to run on everything.

    An eager `.gpkg` resolved at declaration and holds records, so there is nothing to
    copy and nothing to repoint. If this ever fails, the bundler has started archiving
    files nothing opens -- and repointing a spec whose path is now only provenance.
    """

    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        rectangular_voronoi,
    )
    from myflopy.project.bundle import bundle_package_sources
    from myflopy.specs import ModelContext

    grid = rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=1, nper=1))
    context = ModelContext(grid=grid, domain=np.ones((1, grid.ncpl), dtype=int))
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "gis" / "eager.gpkg"
    path.parent.mkdir(parents=True)
    gpd.GeoDataFrame(
        {"name": ["w"], "layer": [1], "head": [9.0], "cond": [1.0]},
        geometry=[box(xmin, ymin, xmax, ymax)], crs=grid.crs,
    ).to_file(path, driver="GPKG")

    eager = mf.ghb.gpkg(path, nper=1, head="head", conductance="cond", context=context)

    assert eager.inputs == {}, "an eager package holds records, not a recipe"

    # The project root is a SIBLING of `gis/`, not its parent. With the source inside
    # the root, `already_inside` (ledger 292) declines to copy it and this test would
    # pass on a bundler that walks eager packages happily -- certifying rather than
    # checking.
    root = tmp_path / "proj"
    inputs = root / "inputs"
    inputs.mkdir(parents=True)
    bundled, report = bundle_package_sources(
        eager, root=root, inputs_dir=inputs, project_root=root
    )

    assert report["copied"] == []
    assert bundled is eager
    assert list(inputs.iterdir()) == []


def test_a_stack_holding_ARRAYS_is_reported_rather_than_raised_over(tmp_path):
    """A bundler that raised would cost the caller every other simulation.

    `Surface.from_array` holds per-cell data, which ledger 193 refused to inline, so
    its stack cannot be serialized at all. Such a model is skipped by the save anyway
    -- the honest answer is a line in `not_bundled`, not an exception that takes the
    whole project down with it.
    """

    from myflopy.project.bundle import bundle_package_sources
    from myflopy.specs import PackageSpec

    stack = mf.LayerStack(top=mf.Surface.from_array(np.full(9, 100.0))).add(
        "sand", thickness=50
    )
    spec = PackageSpec("disv", lambda model, **kw: None, inputs={"layer_stack": stack})
    inputs = tmp_path / "inputs"
    inputs.mkdir()

    bundled, report = bundle_package_sources(
        spec, root=tmp_path, inputs_dir=inputs, project_root=tmp_path
    )

    assert bundled is spec, "nothing was repointed"
    assert report["notes"], "and the caller is TOLD, rather than reading an empty list"
    assert "LayerStack" in report["notes"][0][0]


def test_the_GRASS_cache_is_carried_when_it_exists_and_is_silent_when_it_does_not(tmp_path):
    """`Surface.out` is a WRITE target, so it takes the opposite rule to an input.

    Carried when present: the cache IS the interpolated raster, and an archive without
    it needs GRASS installed to rebuild. Not reported when absent: it has simply not
    been produced yet, and listing it would put a fault in the report for a project in
    a perfectly ordinary state.
    """

    from myflopy.project.bundle import bundle_package_sources
    from myflopy.specs import PackageSpec

    gis = tmp_path / "gis"
    contours = _gpkg(gis / "contours.gpkg")
    cache = _raster(gis / "interp.tif")
    root = tmp_path / "proj"          # a SIBLING of gis/, so a copy really has to happen
    inputs = root / "inputs"
    inputs.mkdir(parents=True)

    def bundle(out):
        surface = mf.Surface(kind="contours", contours=contours, out=out)
        spec = PackageSpec("disv", lambda model, **kw: None, inputs={"s": surface})
        return bundle_package_sources(
            spec, root=root, inputs_dir=inputs, project_root=root
        )

    bundled, report = bundle(cache)
    carried = {os.path.basename(item["from"]) for item in report["copied"]}

    assert carried == {"contours.gpkg", "interp.tif"}
    assert bundled.inputs["s"].out == Path("inputs/interp.tif"), "the recipe now names it"

    _bundled, report = bundle(gis / "never_made.tif")

    assert report["missing"] == [], (
        "a cache that was never produced is not a missing input"
    )
    assert {os.path.basename(item["from"]) for item in report["copied"]} == {
        "contours.gpkg"
    }, "the INPUT is still carried when the cache beside it is absent"
