"""``HeadsPlus.get_obs_cells`` / ``get_obs_heads``: what a points file and ``crs=`` mean.

Ledger 476. ``crs=`` changed no cell, a polygon over several cells got the
first of them in silence, ``get_obs_heads(<str>)`` read the constructor's file in
place of the one given, and a second path on one reader raised.
"""

from __future__ import annotations

import warnings

import geopandas as gpd
import pytest
import shapely

from myflopy.modflow.mf6.headsplus import HeadsPlus


@pytest.fixture(scope="module")
def points(canonical_run, tmp_path_factory):
    model = canonical_run
    centroids = model.vor.gdf_vorPolys.geometry.centroid
    root = tmp_path_factory.mktemp("head_obs_locations")
    crs = model.vor.crs
    paths = {}
    paths["two"] = root / "points.gpkg"
    gpd.GeoDataFrame({"ExploName": ["b", "a"]}, geometry=[centroids.iloc[10], centroids.iloc[200]],
                     crs=crs).to_file(paths["two"])
    paths["other"] = root / "other.gpkg"
    gpd.GeoDataFrame({"ExploName": ["z"]}, geometry=[centroids.iloc[300]], crs=crs).to_file(paths["other"])
    in_wgs = gpd.GeoSeries([centroids.iloc[42]], crs=crs).to_crs(4326).iloc[0]
    paths["naive_wgs"] = root / "naive_wgs.gpkg"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpd.GeoDataFrame({"ExploName": ["w"]}, geometry=[in_wgs]).to_file(paths["naive_wgs"])
    paths["polygon"] = root / "polygon.gpkg"
    over = model.vor.gdf_vorPolys.loc[[0, 1, 22]].union_all().convex_hull.buffer(-1)
    gpd.GeoDataFrame({"ExploName": ["poly"]}, geometry=[over], crs=crs).to_file(paths["polygon"])
    assert isinstance(over, shapely.Polygon)
    return model, paths


def test_crs_labels_a_file_that_carries_none(points):
    model, paths = points
    assert model.hds.get_obs_cells(paths["naive_wgs"], crs="EPSG:4326") == {"w": 42}


def test_crs_naming_another_frame_than_the_files_is_refused(points):
    model, paths = points
    with pytest.raises(ValueError, match="carries its own"):
        model.hds.get_obs_cells(paths["two"], crs="EPSG:4326")


def test_a_location_over_several_cells_is_refused(points):
    model, paths = points
    with pytest.raises(ValueError, match="'poly' \\(\\d+ cells\\) reach several"):
        model.hds.get_obs_cells(paths["polygon"])


def test_get_obs_heads_reads_the_path_it_is_given(points):
    model, paths = points
    reader = HeadsPlus(model=model, obs_path=paths["two"])
    assert list(reader.get_obs_heads(str(paths["other"])).columns) == ["z"]
    assert list(reader.get_obs_heads(paths["two"]).columns) == ["a", "b"]
    assert list(reader.get_obs_heads().columns) == ["a", "b"]  # the locations resolved last


def test_get_obs_heads_columns_are_each_points_head(points):
    model, paths = points
    heads = HeadsPlus(model=model).get_obs_heads(paths["two"])
    table = model.hds.all_heads
    for name, cell in (("a", 200), ("b", 10)):
        expected = table.xs(cell, level="cell")["elev"]
        assert heads[name].sort_index().tolist() == expected.reorder_levels(
            ["layer", "kstpkper"]).sort_index().tolist()
