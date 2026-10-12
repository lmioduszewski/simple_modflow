"""``mf.Basin``: an assumed lakebed where no bathymetry exists (ledger 554).

``Surface.basin(within, depth, *, water_surface, shape="bowl")`` is the water surface
less a depth that is 0 at the shoreline and ``depth`` at each footprint part's pole of
inaccessibility -- a paraboloid (``"bowl"``, mean depth D/2 on a circle) or a cone (mean
D/3). Pinned here:

- the VOLUME it holds is the profile's on two meshes (the area-mean rule; a centroid
  sample missed it by 10% at 200 ft cells), and a bowl holds 1.5x a cone;
- off the footprint's cells it is NaN; a per-feature depth column gives each lake its own;
- what it refuses, at the call and at sampling, naming the feature;
- it serializes (file under ``zone``, an Area inline), digests by content, moves no
  existing surface's payload, and a bundled/drifted project finds its footprint file;
- end to end: ``stack.excavate(bottom=basin)`` + ``mf.lak`` runs in MF6 with a stage,
  and the excavation record's bed IS ``sample_bed(basin)``.
"""

from __future__ import annotations

import contextlib
import io
import json
import pathlib
import shutil
import warnings
from dataclasses import fields

import flopy
import geopandas as gpd
import numpy as np
import pytest
import rasterio
import shapely as shp
import shapely.affinity
from rasterio.transform import from_origin
from shapely.geometry import Point

import myflopy as mf
from myflopy import ModelContext
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.geometry import basin_depths, basin_radius, footprint_cells
from myflopy.specs import _json_value, _spec_value, package_digest, source_paths
from myflopy.surfaces import Surface, sample_bed
from myflopy.workspace import Project

RADIUS, DEPTH = 300.0, 20.0


def _grid(cell):
    n = int(1000 / cell)
    return irregular_voronoi_grid(nrow=n, ncol=n, cell_size=cell)


def _disk(grid, radius=RADIUS):
    x0, y0, x1, y1 = grid.gdf_vorPolys.total_bounds
    return Point((x0 + x1) / 2, (y0 + y1) / 2).buffer(radius, 512)


def _volume(grid, bed, water=100.0):
    """Sum of cell area times depth below ``water`` over the cells the bed covers."""

    finite = np.isfinite(bed)
    return float((grid.gdf_vorPolys.area.to_numpy()[finite] * (water - bed[finite])).sum())


# -- the volume ------------------------------------------------------------------------


@pytest.mark.parametrize("cell", [100.0, 50.0])
@pytest.mark.parametrize("shape, mean", [("bowl", 1 / 2), ("cone", 1 / 3)])
def test_the_basin_holds_the_profiles_volume_on_two_meshes(cell, shape, mean):
    grid = _grid(cell)
    bed = mf.Basin(mf.Area(_disk(grid)), DEPTH, water_surface=100.0, shape=shape).values(grid)
    expected = np.pi * RADIUS**2 * DEPTH * mean
    assert _volume(grid, bed) == pytest.approx(expected, rel=0.005)


def test_a_bowl_holds_one_and_a_half_times_a_cone():
    grid = _grid(50.0)
    lake = mf.Area(_disk(grid))
    bowl = _volume(grid, mf.Basin(lake, DEPTH, water_surface=100.0).values(grid))
    cone = _volume(grid, mf.Basin(lake, DEPTH, water_surface=100.0, shape="cone").values(grid))
    assert bowl / cone == pytest.approx(1.5, rel=0.005)


def test_a_cell_reads_the_area_mean_not_the_centroid():
    """The area mean is what conserves volume; the centroid misses it on coarse cells."""

    grid = _grid(200.0)
    disk = _disk(grid)
    bed = mf.Basin(mf.Area(disk), DEPTH, water_surface=100.0).values(grid)
    cells = np.flatnonzero(np.isfinite(bed))
    centres = np.array([g.centroid.coords[0] for g in grid.gdf_vorPolys.geometry.loc[cells]])
    centroid_depth = basin_depths(centres[:, 0], centres[:, 1], [(disk, DEPTH)])
    area = grid.gdf_vorPolys.area.to_numpy()[cells]
    expected = np.pi * RADIUS**2 * DEPTH / 2
    assert _volume(grid, bed) == pytest.approx(expected, rel=0.005)
    assert abs(float((area * centroid_depth).sum()) / expected - 1) > 0.05, (
        "the premise: a centroid sample misses the volume at 200 ft cells"
    )
    assert np.nanmax(100.0 - bed) < DEPTH, "no cell reaches the declared maximum"


def test_the_deepest_cell_is_at_the_pole_and_cells_off_the_lake_are_nan():
    grid = _grid(50.0)
    disk = _disk(grid)
    bed = mf.Basin(mf.Area(disk), DEPTH, water_surface=100.0).values(grid)
    polys = grid.gdf_vorPolys
    covered = footprint_cells(polys, [disk])[0]
    assert sorted(np.flatnonzero(np.isfinite(bed)).tolist()) == covered
    deepest = int(np.nanargmin(bed))
    assert polys.geometry.loc[deepest].centroid.distance(disk.centroid) < 50.0
    assert basin_radius(disk) == pytest.approx(RADIUS, rel=1e-3)


def test_a_depth_column_gives_each_lake_its_own_depth(tmp_path):
    grid = _grid(50.0)
    x0, y0, _, _ = grid.gdf_vorPolys.total_bounds
    lakes = tmp_path / "lakes.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west", "east"], "max_depth": [10.0, 30.0], "stage": [100.0, 90.0]},
        geometry=[Point(x0 + 250, y0 + 500).buffer(150, 256), Point(x0 + 750, y0 + 500).buffer(150, 256)],
        crs=grid.crs,
    ).to_file(lakes)
    bed = mf.Basin(lakes, "max_depth", water_surface="stage").values(grid)
    polys = grid.gdf_vorPolys
    west, east = footprint_cells(polys, gpd.read_file(lakes).geometry)
    area = polys.area.to_numpy()
    for cells, depth, stage in ((west, 10.0, 100.0), (east, 30.0, 90.0)):
        volume = float((area[cells] * (stage - bed[cells])).sum())
        assert volume == pytest.approx(np.pi * 150**2 * depth / 2, rel=0.01)


def test_a_raster_water_surface_is_sampled_with_its_units(tmp_path):
    grid = _grid(100.0)
    x0, y0, x1, y1 = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "ws_m.tif"
    with rasterio.open(path, "w", driver="GTiff", height=50, width=50, count=1, dtype="float64",
                       crs=grid.crs, transform=from_origin(x0, y1, (x1 - x0) / 50, (y1 - y0) / 50)) as out:
        out.write(np.full((50, 50), 30.0), 1)
    lake = mf.Area(_disk(grid))
    in_metres = mf.Basin(lake, DEPTH, water_surface=mf.Raster(path, units="meters"))
    in_feet = mf.Basin(lake, DEPTH, water_surface=30.0 / 0.3048)
    np.testing.assert_allclose(
        sample_bed(in_metres, grid, np.full(grid.ncpl, 200.0), length_units="feet"),
        sample_bed(in_feet, grid, np.full(grid.ncpl, 200.0), length_units="feet"),
    )


# -- refusals --------------------------------------------------------------------------


@pytest.fixture
def lakes_file(tmp_path):
    grid = _grid(100.0)
    path = tmp_path / "lakes.gpkg"
    gpd.GeoDataFrame({"name": ["a", "b"], "max_depth": [10.0, None], "stage": [100.0, 99.0]},
                     geometry=[_disk(grid, 200.0), shp.affinity.translate(_disk(grid, 50.0), 400.0)],
                     crs=grid.crs).to_file(path)
    return grid, path


def test_a_null_depth_is_refused_naming_the_feature(lakes_file):
    _grid_, path = lakes_file
    with pytest.raises(ValueError, match=r"feature 'row 1' has max_depth=.*NULL depth"):
        mf.Basin(path, "max_depth", water_surface=100.0)


@pytest.mark.parametrize("water, match", [
    (mf.Surface.offset_below(2.0), r"measured from the surface ABOVE"),
    (mf.Flat(100.0).thickness(), r"measured from the surface ABOVE"),
    (mf.Surface.raster("x.tif", fill="propagate"), r"fill='propagate'"),
    (mf.Surface.raster("x.tif", fill="propagate") - 2.0, r"fill='propagate'"),
])
def test_a_relative_or_propagating_water_surface_is_refused(water, match):
    with pytest.raises(ValueError, match=match):
        mf.Basin(mf.Area(shp.box(0, 0, 10, 10)), 2.0, water_surface=water)


@pytest.mark.parametrize("kwargs, error, match", [
    (dict(depth=0.0), ValueError, r"positive, finite length"),
    (dict(depth=float("nan")), ValueError, r"positive, finite length"),
    (dict(depth=True), TypeError, r"not a bool"),
    (dict(depth="max_depth"), ValueError, r"needs within= to be a FILE"),
    (dict(water_surface="stage"), ValueError, r"needs within= to be a FILE"),
    (dict(water_surface=True), TypeError, r"not a bool"),
    (dict(water_surface=[1.0]), TypeError, r"got a list"),
    (dict(shape="dish"), ValueError, r"\['bowl', 'cone'\]"),
    (dict(layer="lakes"), ValueError, r"layer of the footprint FILE"),
    (dict(layer=0), TypeError, r"VECTOR layer"),
    (dict(within=[1, 2, 3]), TypeError, r"not cell numbers"),
    (dict(within=shp.LineString([(0, 0), (5, 5)])), ValueError, r"polygonal with area"),
    (dict(within="no/such.gpkg"), FileNotFoundError, r"does not exist"),
])
def test_what_the_call_refuses(kwargs, error, match):
    arguments = dict(within=mf.Area(shp.box(0, 0, 10, 10)), depth=2.0, water_surface=100.0)
    arguments.update(kwargs)
    within = arguments.pop("within")
    depth = arguments.pop("depth")
    with pytest.raises(error, match=match):
        mf.Basin(within, depth, **arguments)


def test_a_column_or_raster_that_is_neither_is_refused_listing_the_columns(lakes_file):
    _grid_, path = lakes_file
    with pytest.raises(ValueError, match=r"not a column of lakes.gpkg \(name, max_depth, stage\)"):
        mf.Basin(path, 5.0, water_surface="stages")


def test_units_on_the_basin_itself_is_refused():
    basin = mf.Basin(mf.Area(shp.box(0, 0, 10, 10)), 2.0, water_surface=100.0)
    with pytest.raises(ValueError, match=r"declare units= on the water_surface"):
        Surface.from_dict({**basin.to_dict(), "units": "meters"})


def test_overlapping_features_and_a_footprint_off_the_grid_are_refused_at_sampling(tmp_path):
    grid = _grid(100.0)
    overlap = tmp_path / "overlap.gpkg"
    disk = _disk(grid, 200.0)
    gpd.GeoDataFrame({"n": ["a", "b"]}, geometry=[disk, shp.affinity.translate(disk.buffer(-50), 30.0)],
                     crs=grid.crs).to_file(overlap)
    with pytest.raises(ValueError, match=r"overlap \('row 0' and 'row 1'\)"):
        mf.Basin(overlap, 5.0, water_surface=100.0).values(grid)
    far = mf.Basin(mf.Area(shp.box(1e6, 1e6, 1e6 + 50, 1e6 + 50)), 5.0, water_surface=100.0)
    with pytest.raises(ValueError, match=r"cover no grid cell, so they carry no bed"):
        far.values(grid)


# -- serialization and digests ---------------------------------------------------------


def test_a_basin_round_trips_from_a_file_and_from_an_area(lakes_file):
    _grid_, path = lakes_file
    for basin in (
        mf.Basin(path, 5.0, water_surface="stage", shape="cone"),
        mf.Basin(mf.Area(shp.box(0, 0, 10, 10), crs="EPSG:2927"), 2.0,
                 water_surface=mf.Raster("ws.tif", units="meters")),
    ):
        payload = json.loads(json.dumps(_json_value(basin)))
        revived = _spec_value(payload)
        assert revived == basin
        assert type(revived.zone) is type(basin.zone)


def test_from_dict_refuses_what_the_call_refuses():
    payload = mf.Basin(mf.Area(shp.box(0, 0, 10, 10)), 2.0, water_surface=100.0).to_dict()
    with pytest.raises(ValueError, match=r"\['bowl', 'cone'\]"):
        Surface.from_dict({**payload, "shape": "dish"})
    with pytest.raises(ValueError, match=r"positive, finite length"):
        Surface.from_dict({**payload, "value": -1.0})
    with pytest.raises(ValueError, match=r"names a column of"):
        Surface.from_dict({**payload, "depth_column": "d", "value": None})


def test_the_new_fields_write_nothing_for_every_other_surface():
    """So no stored payload or digest of an existing surface moves (ledger 293)."""

    new = {"depth_column", "surface_column", "shape", "zone_layer"}
    assert new <= {item.name for item in fields(Surface)}
    for surface in (
        mf.Flat(1.0), mf.Raster("a.tif", units="meters"), mf.Surface.offset_below(2.0),
        mf.Contours("c.gpkg", z="z", region_vector="r.gpkg"), mf.Min(1.0, "a.tif"), mf.Where("z.gpkg", 1.0, 2.0),
        mf.Isopach(3.0), mf.Toward(1.0, 0.5),
    ):
        assert not new & set(surface.to_dict()), surface.kind


def test_a_lazy_disv_digests_the_basins_footprint_by_content(tmp_path, lakes_file):
    _grid_, path = lakes_file

    def digest(footprint):
        bed = mf.Basin(footprint, 5.0, water_surface="stage")
        stack = mf.LayerStack(top=mf.Flat(100.0)).add("a", bottom=mf.Flat(50.0))
        return package_digest(mf.disv(stack.excavate(bed, mf.Area(shp.box(0, 0, 1, 1)), name="p")))

    moved = tmp_path / "moved"
    moved.mkdir()
    copy = pathlib.Path(shutil.copy2(path, moved / "lakes.gpkg"))
    assert digest(path) == digest(copy), "moving the file is not editing it"
    frame = gpd.read_file(copy)
    frame.loc[0, "stage"] = 101.0
    frame.to_file(copy)
    assert digest(path) != digest(copy), "editing the footprint moves the digest"


# -- the lake, end to end --------------------------------------------------------------


@pytest.fixture
def pond(tmp_path):
    grid = irregular_voronoi_grid(nrow=30, ncol=30, cell_size=20.0, refined=False)
    x0, y0, x1, y1 = grid.gdf_vorPolys.total_bounds
    path = tmp_path / "gis" / "pond.gpkg"
    path.parent.mkdir()
    gpd.GeoDataFrame({"name": ["pond"], "max_depth": [12.0]},
                     geometry=[Point((x0 + x1) / 2, (y0 + y1) / 2).buffer(120.0, 128)],
                     crs=grid.crs).to_file(path)
    return grid, path


def _stack():
    """100 / 95 / 80 / 40: a 5 ft top layer, so a bed deeper than 95 removes its cell."""

    return (mf.LayerStack(top=mf.Flat(100.0)).add("fill", bottom=mf.Flat(95.0))
            .add("upper", bottom=mf.Flat(80.0)).add("lower", bottom=mf.Flat(40.0)))


def test_the_excavation_record_is_the_basins_sampled_bed(pond):
    grid, path = pond
    bed = mf.Basin(path, "max_depth", water_surface=100.0)
    result = _stack().excavate(bottom=bed, within=path, name="pond").build(grid)
    record = result.excavations["pond"]
    sampled = sample_bed(bed, grid, np.full(grid.ncpl, 100.0))
    np.testing.assert_array_equal(record.bed, sampled[list(record.cells)])
    assert np.isfinite(record.bed).all()
    assert (result.idomain[0, list(record.cells)] == 0).any(), "the deep cells are dug"


@pytest.mark.slow
def test_mf6_runs_an_excavated_basin_lake_with_a_stage(pond, tmp_path):
    grid, path = pond
    bed = mf.Basin(path, "max_depth", water_surface=100.0)
    result = _stack().excavate(bottom=bed, within=path, name="pond").build(grid, attach=True)
    context = ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    west, east = np.flatnonzero(xy[:, 0] < xmin + 20), np.flatnonzero(xy[:, 0] > xmax - 20)
    chd = [[(k, int(c)), h] for cells, h in ((west, 97.0), (east, 93.0)) for c in cells
           for k in (1, 2) if result.idomain[k, c] > 0]
    lake = mf.lak(context=context, nper=1, lakes=path, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=1.0, unexcavated="raise")
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: chd}), lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-8, inner_dvclose=1e-10,
               outer_maximum=300, inner_maximum=300, linear_acceleration="BICGSTAB")])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = Project(tmp_path / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    workspace = pathlib.Path(run.workspace)
    listing = (workspace / "m.lst").read_text()
    discrepancy = [abs(float(line.split()[-1])) for line in listing.splitlines()
                   if "PERCENT DISCREPANCY" in line]
    assert len(discrepancy) == 2 and max(discrepancy) < 0.1, "GWF and LAK budgets both close"
    stage = flopy.utils.HeadFile(str(workspace / "m_lak_stage.lak"), text="STAGE").get_data().ravel()
    deepest = float(np.nanmin(record_bed := result.excavations["pond"].bed))
    assert deepest < float(stage[0]) < 100.0, (deepest, float(stage[0]), record_bed.size)
    assert 93.0 < float(stage[0]) < 97.0, "the lake sits between the boundary heads"


# -- the archive and drift -------------------------------------------------------------


def _project(tmp_path, pond_path, boundary):
    bed = mf.Basin(pond_path, "max_depth", water_surface=100.0)
    # The excavation's own footprint is INLINE, so the basin is the only carrier of the
    # file: whatever the walkers find of it, they find through the surface.
    area = mf.Area.coerce(gpd.read_file(pond_path))
    model = mf.gwf("m", grid=mf.grid_ref("g"), packages=[
        mf.disv(_stack().excavate(bottom=bed, within=area, name="pond")),
        mf.npf(k=1.0), mf.ic(strt=99.0), mf.oc(head_filerecord="m.hds")])
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    project = Project(tmp_path / "proj", name="p")
    project.add_grid("g", mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=boundary),
                                              max_area=2_000.0))
    project.add_simulation(simulation)
    return project


@pytest.fixture
def archive(tmp_path):
    gis = tmp_path / "gis"
    gis.mkdir()
    boundary = gis / "dom.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[shp.box(0, 0, 400, 400)], crs=2927).to_file(boundary)
    pond = gis / "pond.gpkg"
    gpd.GeoDataFrame({"name": ["pond"], "max_depth": [12.0]},
                     geometry=[Point(200, 200).buffer(100.0, 64)], crs=2927).to_file(pond)
    return gis, boundary, pond


def test_the_walkers_find_the_basins_footprint(archive):
    _gis, _boundary, pond = archive
    bed = mf.Basin(pond, "max_depth", water_surface=100.0)
    stack = _stack().excavate(bottom=bed, within=mf.Area(shp.box(0, 0, 1, 1)), name="p")
    assert {p.split("/")[-1] for p in source_paths(mf.disv(stack).inputs)} == {"pond.gpkg"}


def test_a_bundled_basin_project_rebuilds_with_its_sources_deleted(tmp_path, archive):
    gis, boundary, pond = archive
    project = _project(tmp_path, pond, boundary)
    before = project.prepare_run("before", "s").built.built_model("m").model.get_package("disv")
    report = project.save(bundle_sources=True)
    assert "pond.gpkg" in {str(source).split("/")[-1] for source, _ in report.bundled}
    shutil.rmtree(gis)
    after = Project.load(tmp_path / "proj").prepare_run("after", "s")
    disv = after.built.built_model("m").model.get_package("disv")
    for name in ("top", "botm", "idomain"):
        np.testing.assert_array_equal(np.asarray(getattr(disv, name).get_data()),
                                      np.asarray(getattr(before, name).get_data()))
    assert (np.asarray(disv.idomain.get_data())[0] == 0).any(), "the pond is dug"


def test_drift_reports_an_edited_basin_footprint(tmp_path, archive):
    _gis, boundary, pond = archive
    project = _project(tmp_path, pond, boundary)
    project.prepare_run("r1", "s")
    assert project.runs.drift().set_index("key").loc["disv", "source"] == "ok"
    frame = gpd.read_file(pond)
    frame.loc[0, "max_depth"] = 15.0
    frame.to_file(pond)
    assert project.refresh_runs().drift().set_index("key").loc["disv", "source"].startswith("edited")


# -- a basin is not a layer surface on its own (ledger 555) -------------------


def _small_lake():
    grid = irregular_voronoi_grid(nrow=20, ncol=20, cell_size=20.0, refined=False)
    x0, y0, x1, y1 = grid.gdf_vorPolys.total_bounds
    lake = mf.Area(Point((x0 + x1) / 2, (y0 + y1) / 2).buffer(80.0, 64))
    return grid, mf.Basin(lake, 12.0, water_surface=100.0)


@pytest.mark.parametrize("where", ["unit bottom", "top", "shifted unit bottom"])
def test_a_basin_as_a_layer_surface_is_refused_where_it_has_no_bed(where):
    """Off its lake a basin is NaN, and a stack wrote that NaN into the DISV: MF6 ran it
    to normal termination with all 800 heads NaN (measured, review of R3-basin)."""

    grid, bed = _small_lake()
    stack = {
        "unit bottom": lambda: mf.LayerStack(top=mf.Flat(100.0)).add("a", bottom=bed),
        "top": lambda: mf.LayerStack(top=bed).add("a", bottom=mf.Flat(50.0)),
        "shifted unit bottom": lambda: mf.LayerStack(top=mf.Flat(100.0)).add("a", bottom=bed - 1.0),
    }[where]().add("b", bottom=mf.Flat(40.0))
    with pytest.raises(ValueError, match=r"reads an mf.Basin, which has no bed off its lake.*stack.excavate"):
        stack.build(grid)


def test_a_basin_composed_onto_a_surface_with_a_value_everywhere_is_a_layer_surface():
    grid, bed = _small_lake()
    result = (mf.LayerStack(top=mf.Flat(100.0)).add("a", bottom=mf.Min(mf.Flat(95.0), bed))
              .add("b", bottom=mf.Flat(40.0)).build(grid))
    lake = np.isfinite(bed.values(grid))
    assert np.isfinite(result.botm).all()
    np.testing.assert_allclose(result.botm[0, ~lake], 95.0)
    assert (result.botm[0, lake] < 95.0).any(), "the bed is baked in under the lake"
