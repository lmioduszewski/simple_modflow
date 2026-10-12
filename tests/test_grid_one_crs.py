"""One grid keeps one CRS, whatever writes to it (ledger 431).

Ledger 368 made ``grid.crs = ...`` relabel the grid and every frame it caches. Four
routes still split one grid between two CRSs:

* ``grid.gdf_vorPolys = frame`` stored a frame in any CRS as the cell cache;
* ``grid.gdf_topbtm = frame`` likewise for the layer elevations;
* ``grid.get_gdf_vorPolys(crs=other)`` built the cell frame under ``other``'s label
  and installed it as the cache;
* ``HeadsPlus`` and the legacy ``Boundaries`` copied ``vor.crs`` when built, so a
  grid relabelled afterwards left them on the old label.

The first three are refused naming ``grid.crs = ...`` (relabel) and ``to_crs``
(reproject); the last two read the grid's CRS live.
"""

from __future__ import annotations

import copy

import geopandas as gpd
import pandas as pd
import pytest

from myflopy.modflow.mf6.boundaries import Boundaries
from myflopy.modflow.mf6.headsplus import HeadsPlus
from tests.test_usg_export_adjacency import square_grid


def test_a_cell_frame_in_another_crs_is_refused_and_the_cache_kept():
    grid = square_grid()
    cells = grid.gdf_vorPolys
    with pytest.raises(ValueError, match=r"grid.gdf_vorPolys = frame.*grid.crs = 'EPSG:4326'"):
        grid.gdf_vorPolys = cells.to_crs("EPSG:4326")
    assert grid.gdf_vorPolys is cells and grid.gdf_vorPolys.crs == "EPSG:2927"


def test_a_cell_frame_with_no_crs_takes_the_grids_label_and_stays_the_same_object():
    grid = square_grid()
    frame = gpd.GeoDataFrame(geometry=grid.gdf_vorPolys.geometry.to_numpy())
    grid.gdf_vorPolys = frame
    assert grid.gdf_vorPolys is frame and frame.crs == "EPSG:2927"


def test_the_grids_own_crs_under_another_spelling_is_accepted():
    grid = square_grid()
    grid.gdf_vorPolys = grid.gdf_vorPolys.set_crs(2927, allow_override=True)
    assert grid.get_gdf_vorPolys(crs=2927).crs == "EPSG:2927"


def test_a_layer_table_in_another_crs_is_refused():
    grid = square_grid()
    table = gpd.GeoDataFrame(
        {"geometry": grid.gdf_vorPolys.geometry.to_numpy(), 0: 10.0, 1: 0.0}, crs="EPSG:2926"
    )
    with pytest.raises(ValueError, match=r"grid.gdf_topbtm = frame.*frame.to_crs\('EPSG:2927'\)"):
        grid.gdf_topbtm = table


def test_building_the_cell_frame_under_another_label_is_refused():
    grid = square_grid()
    with pytest.raises(ValueError, match=r"get_gdf_vorPolys\(crs='EPSG:2926'\).*grid.crs = 'EPSG:2926'"):
        grid.get_gdf_vorPolys(crs="EPSG:2926")
    assert grid.gdf_vorPolys.crs == "EPSG:2927"


def test_a_legacy_boundary_follows_a_relabel_of_its_grid(tmp_path):
    """Its source is reprojected to the NEW label, not kept on the old one."""

    grid = square_grid()
    source = tmp_path / "drain.gpkg"
    gpd.GeoDataFrame(
        {"uid": ["d"]}, geometry=[grid.gdf_vorPolys.geometry.iloc[5]], crs="EPSG:2927"
    ).to_file(source)
    boundary = Boundaries(vor=grid, shp_gpkg=source, uid="uid")
    assert boundary.gdf.crs == "EPSG:2927"
    assert boundary.intersections.loc["d"]

    grid.crs = "EPSG:2926"
    assert boundary.crs == "EPSG:2926"
    assert boundary.gdf.crs == "EPSG:2926"
    # The source moved ~1.1 km on being reprojected to the new label, off this grid.
    assert boundary.intersections.loc["d"] == []


def test_a_heads_reader_follows_a_relabel_of_its_grid(canonical_run):
    grid = copy.deepcopy(canonical_run.vor)
    heads = HeadsPlus(canonical_run.hds.hds_path, vor=grid)
    assert heads.crs == grid.crs
    grid.crs = "EPSG:2926"
    assert heads.crs == "EPSG:2926"
    with pytest.raises(AttributeError, match=r"grid.crs = 'EPSG:2855'"):
        heads.crs = "EPSG:2855"


@pytest.mark.parametrize("builder", ["recharge", "ghb"])
def test_a_legacy_builders_own_cell_sets_follow_a_relabel_of_its_grid(tmp_path, builder):
    """A subclass's cached cells followed the OLD label beside a refreshed source.

    ``RechargeFromShp.cell_ids`` and ``GHBFromVector.cell_ids`` cache the cells the
    source covers; the base class dropped only its own caches on a relabel, so the
    object answered ``intersections`` on the new label and ``cell_ids`` on the old
    one; and ``GHBFromVector.line_ghb`` was built with a COPY of the grid's CRS
    (Review of U6, ledger 431).
    """

    from myflopy.modflow.mf6.ghb import GHBFromVector
    from myflopy.modflow.mf6.recharge import RechargeFromShp

    grid = square_grid()
    source = tmp_path / "zone.gpkg"
    gpd.GeoDataFrame(
        {"uid": ["z"], "rate": [1.0]}, geometry=[grid.gdf_vorPolys.geometry.iloc[5]],
        crs="EPSG:2927",
    ).to_file(source)
    if builder == "recharge":
        boundary = RechargeFromShp(vor=grid, shp_gpkg=source, uid="uid", rch_fields=["rate"])
    else:
        boundary = GHBFromVector(
            vor=grid, shp_gpkg=source, line_shp=source, uid="uid", crs=None, fields=["rate"]
        )
    assert 5 in boundary.cell_ids["z"]
    line = boundary.line_ghb if builder == "ghb" else None

    grid.crs = "EPSG:2926"
    if line is not None:  # the line source it builds follows the grid too
        assert line.crs == "EPSG:2926"
    assert boundary.intersections.loc["z"] == []
    assert boundary.cell_ids == {"z": []}


def test_a_legacy_boundary_refuses_a_crs_other_than_its_grids(tmp_path):
    """``Boundaries(crs=X)`` reprojected the source to X and intersected it against
    cells in the grid's CRS -- two frames, so a source over cell 5 found no cell
    (ledger 466). Another spelling of the grid's own CRS is accepted."""

    grid = square_grid()
    source = tmp_path / "drain.gpkg"
    gpd.GeoDataFrame(
        {"uid": ["d"]}, geometry=[grid.gdf_vorPolys.geometry.iloc[5]], crs="EPSG:2927"
    ).to_file(source)
    with pytest.raises(ValueError, match=r"Boundaries\(crs='EPSG:2926'\).*leave crs= out"):
        Boundaries(vor=grid, shp_gpkg=source, uid="uid", crs="EPSG:2926")
    assert 5 in Boundaries(vor=grid, shp_gpkg=source, uid="uid", crs=2927).intersections.loc["d"]

    # A relabel that makes the two differ is refused on the next read, not followed.
    boundary = Boundaries(vor=grid, shp_gpkg=source, uid="uid", crs=2927)
    grid.crs = "EPSG:2926"
    with pytest.raises(ValueError, match=r"grid.crs = 2927"):
        boundary.intersections


def test_budget_polygons_in_another_crs_find_their_cells(canonical_run, tmp_path):
    """``budget_obs_df`` passed the polygons' OWN CRS as ``Boundaries(crs=)``, so
    polygons stored in another CRS than the grid's were intersected in their own
    frame against cells in the grid's, and summed no flow (ledger 466)."""

    from myflopy.modflow.mf6.budget_tables import budget_obs_df

    model = canonical_run
    budget = model.bud("chd")
    cells = model.vor.gdf_vorPolys
    zone = cells.geometry.loc[sorted(set(budget.df.index.get_level_values(0)))[:3]].union_all()
    frames = {}
    for crs in (cells.crs, "EPSG:4326"):
        path = tmp_path / f"zone_{str(crs).replace(':', '_')}.gpkg"
        gpd.GeoDataFrame({"name": ["zone"]}, geometry=[zone], crs=cells.crs).to_crs(crs).to_file(path)
        frames[crs] = budget_obs_df(budget, model=model, shp_gpkg=path)
    grid_frame, other = frames[cells.crs], frames["EPSG:4326"]
    assert grid_frame["zone"].abs().sum() > 0
    pd.testing.assert_frame_equal(grid_frame, other)


def test_a_geopackage_source_follows_a_relabel_of_its_grid(tmp_path):
    """``GeoPackageSource.gdf`` cached the layer reprojected to the grid's CRS at its
    first read, so a grid relabelled afterwards kept resolving the layer in the OLD
    frame (ledger 466). It is reprojected from the file's own frame now."""

    from myflopy.geopackage import GeoPackageSource
    from myflopy.specs import ModelContext

    grid = square_grid()
    path = tmp_path / "drain.gpkg"
    drawn = gpd.GeoDataFrame(
        {"name": ["d"]}, geometry=[grid.gdf_vorPolys.geometry.iloc[5].centroid], crs="EPSG:2927"
    )
    drawn.to_file(path)
    source = GeoPackageSource(path, ModelContext(grid=grid), nper=1, layer_field=None)
    first = source.gdf
    assert first.crs == "EPSG:2927" and source.gdf is first

    grid.crs = "EPSG:2926"
    moved = source.gdf
    assert moved.crs == "EPSG:2926"
    expected = drawn.to_crs("EPSG:2926").geometry.iloc[0]
    assert moved.geometry.iloc[0].distance(expected) < 1e-6
    # Back to the first label: the file's own frame again, not a round trip of the second.
    grid.crs = "EPSG:2927"
    assert source.gdf.geometry.iloc[0].equals_exact(drawn.geometry.iloc[0], 1e-9)
