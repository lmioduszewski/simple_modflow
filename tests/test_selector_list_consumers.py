"""Three callers read the shared selector's answer as a flat list of cells.

``get_vor_cells_as_series`` answers one entry PER FEATURE, each a LIST of cells
(the shape ledger 192 kept, measured by ledger 235). ``build_idomain`` via
``assign_voronoi_cells_to_layers``, ``InterpolatedSurface(clip=)`` and
``get_overlapping_area(shp_gpkg=)`` each treated that as a flat list of cells, so
all three raised on EVERY input: ``TypeError: cannot use 'list' as a set element``
for the first, and ``KeyError: None of [Index([(5,)], ...)]`` for the other two.

That made ``SimulationBase(idomain_path=...)`` unusable: reading ``model.idomain``,
``model.idomain_gdf`` or ``model.inactive_cells`` raised.

Three more callers had written the flattening out for themselves and worked --
``adjust_top_btm_overlaps``, ``_inactive_from_file`` (ledger 236) and
``RegionRegistry.add_from_geometry`` -- so the fix is one ``selected_cells``, and
the tests below pin that each of the six goes through it rather than growing its
own copy back. ``_inactive_from_file`` had a second copy too: the extents-and-CRS
sentence of ledger 211, which was already two copies before it.
"""

from __future__ import annotations

import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, box

from myflopy.modflow.mf6.grid import geometry as geometry_module
from myflopy.modflow.mf6.grid import selection as selection_module
from myflopy.modflow.mf6.grid.geometry import selected_cells, where_they_are
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.selection import get_grid_edge_cells
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simulation import indexing as indexing_module
from myflopy.modflow.mf6.simulation import regions as regions_module
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.indexing import build_idomain
from myflopy.modflow.utils.datatypes.readers import assign_voronoi_cells_to_layers

SELECTION_LOGGER = "myflopy.modflow.mf6.grid.selection"

#: Exactly cell 5 of the 4x4 grid below: the shape that made every one of these
#: callers pick nine cells before ledger 235, and raise after it.
CELL_5 = box(25.0, 25.0, 50.0, 50.0)
#: Well inside cell 5, covering 36% of it -- the area rule keeps it, and only it.
PART_OF_CELL_5 = box(30.0, 30.0, 45.0, 45.0)
AWAY = box(500.0, 500.0, 600.0, 600.0)


def _grid(n: int = 4, length: float = 100.0) -> VoronoiGridPlus:
    """n x n square cells, row-major from the bottom left; cell 5 is (25..50, 25..50)."""

    d = length / n
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    # Two layers, so `build_idomain` has something to shape its answer to, and so
    # a `surf_type='lyr'` surface has an elevation to read.
    polys = vor.gdf_vorPolys
    vor.gdf_topbtm = gpd.GeoDataFrame(
        {0: np.full(len(polys), 100.0), 1: np.full(len(polys), 50.0), 2: np.zeros(len(polys))},
        geometry=polys.geometry,
        crs=polys.crs,
    )
    return vor


@pytest.fixture
def grid():
    return _grid()


def _idomain_file(path: Path, features, crs) -> Path:
    """A layer-tagged polygon file, the shape `idomain_path=` takes."""

    gpd.GeoDataFrame(
        {"layer": [layer for layer, _ in features]},
        geometry=[geometry for _, geometry in features],
        crs=crs,
    ).to_file(path, driver="GPKG")
    return path


# --------------------------------------------------------------------------- #
# 1. build_idomain
# --------------------------------------------------------------------------- #


def test_build_idomain_switches_off_the_cell_a_polygon_covers(grid, tmp_path):
    """Was `TypeError: cannot use 'list' as a set element` for every file."""

    path = _idomain_file(tmp_path / "idomain.gpkg", [("0 1", CELL_5)], grid.crs)

    idomain, _ = build_idomain(grid, path)

    assert np.asarray(idomain).shape == (2, 16)
    assert [cell for cell, value in enumerate(idomain[0]) if value == 0] == [5]
    assert [cell for cell, value in enumerate(idomain[1]) if value == 0] == [5]


def test_build_idomain_covers_by_area_so_a_partial_polygon_still_takes_one_cell(grid, tmp_path):
    """The area rule decides: a polygon inside one cell reaches that cell, alone."""

    path = _idomain_file(tmp_path / "idomain.gpkg", [("0", PART_OF_CELL_5)], grid.crs)

    idomain, _ = build_idomain(grid, path)

    assert [cell for cell, value in enumerate(idomain[0]) if value == 0] == [5]
    assert sum(idomain[1]) == 16


def test_build_idomain_switches_off_only_the_layers_a_feature_names(grid, tmp_path):
    path = _idomain_file(tmp_path / "idomain.gpkg", [("1", CELL_5)], grid.crs)

    idomain, _ = build_idomain(grid, path)

    assert sum(idomain[0]) == 16
    assert [cell for cell, value in enumerate(idomain[1]) if value == 0] == [5]


def test_build_idomain_asks_the_selector_once_for_the_whole_file(grid, tmp_path):
    """Once per ROW, a row landing off the grid was a whole-input miss, so the
    ledger-211 policy fired on a single row."""

    path = _idomain_file(
        tmp_path / "idomain.gpkg",
        [("0 1", CELL_5), ("0 1", box(50.0, 25.0, 75.0, 50.0)), ("1", box(0.0, 0.0, 25.0, 25.0))],
        grid.crs,
    )
    calls = []
    original = grid.get_vor_cells_as_series
    grid.get_vor_cells_as_series = lambda *a, **k: (calls.append(1), original(*a, **k))[1]

    idomain, _ = build_idomain(grid, path)

    assert len(calls) == 1
    assert [cell for cell, value in enumerate(idomain[0]) if value == 0] == [5, 6]
    assert [cell for cell, value in enumerate(idomain[1]) if value == 0] == [0, 5, 6]


def test_build_idomain_keeps_the_features_that_land_and_names_the_one_that_does_not(
    grid, tmp_path, caplog
):
    path = _idomain_file(tmp_path / "idomain.gpkg", [("0", CELL_5), ("0", AWAY)], grid.crs)

    with caplog.at_level(logging.WARNING, logger=SELECTION_LOGGER):
        idomain, _ = build_idomain(grid, path)

    assert [cell for cell, value in enumerate(idomain[0]) if value == 0] == [5]
    assert any("reach no grid cell" in record.getMessage() for record in caplog.records)


def test_build_idomain_keeps_a_layer_that_lands_when_another_layers_feature_misses(
    grid, tmp_path, caplog
):
    """A file is refused only when it lands NOWHERE. One layer's feature missing
    is a partial miss: the guard reads `not any(...)`, and `not all(...)` would
    throw away a file whose other layers are perfectly good."""

    path = _idomain_file(tmp_path / "idomain.gpkg", [("0", CELL_5), ("1", AWAY)], grid.crs)

    with caplog.at_level(logging.WARNING, logger=SELECTION_LOGGER):
        idomain, _ = build_idomain(grid, path)

    assert [cell for cell, value in enumerate(idomain[0]) if value == 0] == [5]
    assert sum(idomain[1]) == 16
    assert any("reach no grid cell" in record.getMessage() for record in caplog.records)


def test_build_idomain_refuses_a_file_that_lands_nowhere(grid, tmp_path):
    """Ledger 211. It used to answer "every cell active", which runs to normal
    termination with the idomain doing nothing."""

    path = _idomain_file(tmp_path / "idomain.gpkg", [("0 1", AWAY)], grid.crs)

    with pytest.raises(ValueError) as raised:
        build_idomain(grid, path)

    message = str(raised.value)
    assert "lands on a grid cell" in message
    assert "outside the grid" in message
    # Both extents and both CRSs, which is what a location mismatch is found by.
    assert "500.0, 500.0, 600.0, 600.0" in message
    assert "0.0, 0.0, 100.0, 100.0" in message
    assert "EPSG:2927" in message


def test_build_idomain_says_so_when_the_file_carries_no_crs(grid, tmp_path):
    path = _idomain_file(tmp_path / "idomain.gpkg", [("0", AWAY)], None)

    with pytest.raises(ValueError, match="the file has NO CRS"):
        build_idomain(grid, path)


def test_build_idomain_names_the_reprojection_when_the_file_is_in_another_crs(grid, tmp_path):
    """The reason ledger 211 prints two CRSs: a file in degrees over a projected
    grid is the usual cause of a total miss, and "both are in EPSG:2927" of a file
    that was in fact reprojected sends the reader looking at the wrong thing."""

    degrees = gpd.GeoSeries([AWAY], crs=grid.crs).to_crs("EPSG:4326")
    path = _idomain_file(tmp_path / "degrees.gpkg", [("0", degrees.iloc[0])], "EPSG:4326")

    with pytest.raises(ValueError) as raised:
        build_idomain(grid, path)

    assert "reprojected from EPSG:4326 to the grid's EPSG:2927" in str(raised.value)


def test_a_reprojected_file_that_does_land_still_switches_its_cell_off(grid, tmp_path):
    """The other half of that: the selector reprojects, so a 4326 file WORKS."""

    degrees = gpd.GeoSeries([CELL_5], crs=grid.crs).to_crs("EPSG:4326")
    path = _idomain_file(tmp_path / "degrees.gpkg", [("0", degrees.iloc[0])], "EPSG:4326")

    idomain, _ = build_idomain(grid, path)

    assert [cell for cell, value in enumerate(idomain[0]) if value == 0] == [5]


def test_a_model_given_an_idomain_path_now_has_a_domain(grid, tmp_path):
    """End to end. `SimulationBase(idomain_path=)` raised `TypeError` on EVERY
    read before this -- `model.idomain`, `model.idomain_gdf` and
    `model.inactive_cells` alike -- so the parameter could not be used at all."""

    path = _idomain_file(tmp_path / "idomain.gpkg", [("0 1", CELL_5)], grid.crs)
    model = SimulationBase(
        name="idomain_end_to_end", mf_folder_path=tmp_path, vor=grid, idomain_path=path
    )

    assert np.asarray(model.idomain).shape == (2, 16)
    # Flattened over both layers: cell 5 of layer 0, and cell 5 of layer 1.
    assert model.inactive_cells == [5, 21]
    assert list(model.idomain_gdf["0"]) == [[5]]


def test_assign_voronoi_cells_to_layers_puts_a_flat_list_of_cells_in_each_column(grid):
    """The value is the cells THAT ROW reaches -- `[5]`, never `[[5]]`."""

    gdf = gpd.GeoDataFrame(
        {"layer": ["0 1", "1"]}, geometry=[CELL_5, box(0.0, 0.0, 25.0, 25.0)], crs=grid.crs
    )

    out, columns = assign_voronoi_cells_to_layers(gdf, grid, return_col_names=True)

    assert columns == ["0", "1"]
    assert list(out["0"]) == [[5], None]
    assert list(out["1"]) == [[5], [0]]


def test_assign_voronoi_cells_to_layers_maps_back_by_position_not_by_label(grid):
    """A repeated index label would otherwise give both rows one row's cells."""

    gdf = gpd.GeoDataFrame(
        {"layer": ["0", "0"]},
        geometry=[CELL_5, box(0.0, 0.0, 25.0, 25.0)],
        crs=grid.crs,
        index=[7, 7],
    )

    out = assign_voronoi_cells_to_layers(gdf, grid)

    assert list(out["0"]) == [[5], [0]]


def test_assigning_cells_leaves_the_frames_crs_alone(grid):
    """`build_idomain`'s refusal reads the CRS off the frame AFTER this ran, and
    says "the file has NO CRS" when there is none. That only stays true while the
    selector does not write the grid's CRS back into the frame it was handed."""

    gdf = gpd.GeoDataFrame({"layer": ["0"]}, geometry=[CELL_5], crs=None)

    out = assign_voronoi_cells_to_layers(gdf, grid)

    assert gdf.crs is None and out.crs is None


# --------------------------------------------------------------------------- #
# 2. InterpolatedSurface(clip=)
# --------------------------------------------------------------------------- #


def test_a_clipped_surface_takes_the_cells_the_clip_covers(grid):
    """Was `KeyError: None of [Index([(5,)], ...)] are in the [index]`."""

    surface = InterpolatedSurface(vor=grid, surf_type="lyr", clip=CELL_5)

    assert surface.clipped_cells == [5]
    assert surface.xs.tolist() == [37.5]
    assert surface.ys.tolist() == [37.5]


def test_a_clip_covering_part_of_a_cell_still_takes_that_cell(grid):
    surface = InterpolatedSurface(vor=grid, surf_type="lyr", clip=PART_OF_CELL_5)

    assert surface.clipped_cells == [5]


def test_a_clip_read_from_a_file_reaches_the_same_cell(grid, tmp_path):
    path = tmp_path / "clip.gpkg"
    gpd.GeoDataFrame(geometry=[CELL_5], crs=grid.crs).to_file(path, driver="GPKG")

    surface = InterpolatedSurface(vor=grid, surf_type="lyr", clip=path)

    assert surface.clipped_cells == [5]


def test_zs_reads_the_clipped_cells_even_when_asked_for_first(grid):
    """`zs` read `_clipped_cells`, which only `xs`/`ys` ever filled, so asking for
    it first took all 16 cells against `xs`'s one."""

    surface = InterpolatedSurface(vor=grid, surf_type="lyr", clip=CELL_5)

    assert len(surface.zs) == 1
    assert len(surface.zs) == len(surface.xs)


def test_no_clip_still_takes_every_cell(grid):
    surface = InterpolatedSurface(vor=grid, surf_type="lyr")

    assert surface.clipped_cells is None
    assert len(surface.xs) == 16
    assert len(surface.zs) == 16


def test_a_clip_that_covers_no_cell_is_refused_with_both_extents(grid):
    surface = InterpolatedSurface(vor=grid, surf_type="lyr", clip=AWAY)

    with pytest.raises(ValueError) as raised:
        _ = surface.xs

    message = str(raised.value)
    assert "covers no grid cell" in message
    assert "outside the grid" in message
    assert "500.0, 500.0, 600.0, 600.0" in message
    assert "0.0, 0.0, 100.0, 100.0" in message


# --------------------------------------------------------------------------- #
# 3. get_overlapping_area
# --------------------------------------------------------------------------- #


def test_overlapping_area_is_the_area_of_the_cells_a_polygon_covers(grid):
    """Was `KeyError: None of [Index([(5,)], ...)] are in the [index]`."""

    assert grid.get_overlapping_area(shp_gpkg=CELL_5) == pytest.approx(625.0)


def test_overlapping_area_counts_a_whole_cell_a_polygon_only_partly_covers(grid):
    """It is the CELLS' area, not the polygon's: 625, not 225."""

    assert grid.get_overlapping_area(shp_gpkg=PART_OF_CELL_5) == pytest.approx(625.0)


def test_overlapping_area_reads_a_file(grid, tmp_path):
    path = _idomain_file(tmp_path / "zone.gpkg", [("0", box(25.0, 25.0, 75.0, 50.0))], grid.crs)

    assert grid.get_overlapping_area(shp_gpkg=path) == pytest.approx(1250.0)


def test_overlapping_area_of_a_polygon_that_reaches_nothing_is_zero(grid, caplog):
    with caplog.at_level(logging.WARNING, logger=SELECTION_LOGGER):
        area = grid.get_overlapping_area(shp_gpkg=AWAY)

    assert area == 0.0
    assert any("reach no grid cell" in record.getMessage() for record in caplog.records)


def test_explicit_cell_ids_still_work(grid):
    assert grid.get_overlapping_area(cell_list=[5, 6]) == pytest.approx(1250.0)


# --------------------------------------------------------------------------- #
# 4. One copy of the flattening, and one copy of the "where they are" wording
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "selection, expected",
    [
        (pd.Series([[5], [6, 5]], name="cells"), [5, 6]),
        (pd.Series([], name="cells", dtype=object), []),
        ([[5], [0]], [0, 5]),
        ([5, 0], [0, 5]),
        ([np.array([3, 1])], [1, 3]),
        ([None, [2]], [2]),
    ],
)
def test_selected_cells_flattens_sorts_and_deduplicates(selection, expected):
    assert selected_cells(selection) == expected


def test_every_list_of_lists_consumer_imports_the_one_flattener():
    """Identity, so that a caller quietly growing its own copy fails here.

    Six callers want the flat answer. Two live in geometry.py beside it
    (`get_overlapping_area`, `adjust_top_btm_overlaps`); these four import it --
    the clipped surface, `build_idomain`, `_inactive_from_file` (selection.py) and
    `add_from_geometry` (regions.py)."""

    from myflopy.modflow.mf6.grid import interpolated_surface as surface_module

    for module in (surface_module, selection_module, indexing_module, regions_module):
        assert module.selected_cells is geometry_module.selected_cells


def test_get_overlapping_area_goes_through_the_one_flattener(grid, monkeypatch):
    seen = []
    original = geometry_module.selected_cells
    monkeypatch.setattr(
        geometry_module, "selected_cells", lambda s: (seen.append(s), original(s))[1]
    )

    grid.get_overlapping_area(shp_gpkg=CELL_5)

    assert len(seen) == 1


def test_adjust_top_btm_overlaps_goes_through_the_one_flattener(grid, tmp_path, monkeypatch):
    """It carried its own nine-line copy of the flattening."""

    path = tmp_path / "adjust.gpkg"
    gpd.GeoDataFrame(geometry=[CELL_5], crs=grid.crs).to_file(path, driver="GPKG")
    seen = []
    original = geometry_module.selected_cells
    monkeypatch.setattr(
        geometry_module, "selected_cells", lambda s: (seen.append(s), original(s))[1]
    )

    geometry_module.adjust_top_btm_overlaps(grid, shp_path=path, layer_bottom_name=1)

    assert len(seen) == 1


def test_the_clipped_surface_goes_through_the_one_flattener(grid, monkeypatch):
    """And through it ONCE per surface. `clipped_cells` memoises, so reading
    `xs`, `ys` and `zs` does not set the spatial join up three more times --
    dropping the memo leaves every test above green and the selector running per
    read."""

    from myflopy.modflow.mf6.grid import interpolated_surface as surface_module

    seen = []
    original = surface_module.selected_cells
    monkeypatch.setattr(
        surface_module, "selected_cells", lambda s: (seen.append(s), original(s))[1]
    )

    surface = InterpolatedSurface(vor=grid, surf_type="lyr", clip=CELL_5)
    assert surface.clipped_cells == [5]
    assert len(surface.xs) == 1
    assert len(surface.ys) == 1
    assert len(surface.zs) == 1

    assert len(seen) == 1


def test_the_grid_edge_reads_an_idomain_file_through_the_one_flattener(grid, tmp_path, monkeypatch):
    """`_inactive_from_file` (ledger 236) is the third copy of the flattening."""

    path = _idomain_file(tmp_path / "edge.gpkg", [("0", CELL_5)], grid.crs)
    seen = []
    original = selection_module.selected_cells
    monkeypatch.setattr(
        selection_module, "selected_cells", lambda s: (seen.append(s), original(s))[1]
    )

    edge = get_grid_edge_cells(grid, idomain_path=path)

    assert len(seen) == 1
    assert 5 not in edge


def test_build_idomain_reads_its_columns_through_the_one_flattener(grid, tmp_path, monkeypatch):
    """It carried its own `sorted(set(cell for sublist in ...))`, five lines above
    the refusal this change added."""

    path = _idomain_file(tmp_path / "idomain.gpkg", [("0", CELL_5)], grid.crs)
    seen = []
    original = indexing_module.selected_cells
    monkeypatch.setattr(
        indexing_module, "selected_cells", lambda s: (seen.append(s), original(s))[1]
    )

    build_idomain(grid, path)

    # One call per layer column, and the file names one.
    assert len(seen) == 1


def test_a_region_from_geometry_goes_through_the_one_flattener(grid, tmp_path, monkeypatch):
    """`RegionRegistry.add_from_geometry` is the fifth copy."""

    model = SimulationBase(name="region_flattener", mf_folder_path=tmp_path, vor=grid)
    seen = []
    original = regions_module.selected_cells
    monkeypatch.setattr(
        regions_module, "selected_cells", lambda s: (seen.append(s), original(s))[1]
    )

    region = model.add_region_from_geometry("cell5", CELL_5)

    assert len(seen) == 1
    assert region.cells == [5]


def test_where_they_are_names_both_extents_and_both_crss(grid):
    features = gpd.GeoDataFrame(geometry=[CELL_5], crs=grid.crs)

    said = where_they_are(features, grid.gdf_vorPolys, source_crs=features.crs)

    assert said == (
        "The features span (25.0, 25.0, 50.0, 50.0) and the grid "
        "(0.0, 0.0, 100.0, 100.0); both are in EPSG:2927."
    )


def test_where_they_are_takes_a_bare_geometry_too(grid):
    assert "The clip bounds span (25.0, 25.0, 50.0, 50.0)" in where_they_are(
        CELL_5, grid.gdf_vorPolys, noun="clip bounds"
    )


def test_where_they_are_reports_a_reprojection_and_a_missing_crs(grid):
    features = gpd.GeoDataFrame(geometry=[Point(0.0, 0.0)], crs="EPSG:4326")

    assert "reprojected from EPSG:4326 to the grid's EPSG:2927" in where_they_are(
        features, grid.gdf_vorPolys, source_crs="EPSG:4326", grid_crs="EPSG:2927"
    )
    assert "the features carry NO CRS, so they were read as the grid's" in where_they_are(
        features, grid.gdf_vorPolys, grid_crs="EPSG:2927"
    )
    assert "the file has NO CRS, so its coordinates were read" in where_they_are(
        features,
        grid.gdf_vorPolys,
        grid_crs="EPSG:2927",
        no_crs_subject="the file has NO CRS, so its coordinates",
    )


def test_the_grid_crs_defaults_to_the_polygons_own(grid):
    """A caller that reads it off a grid OBJECT gets `None` when that object has
    no `.crs`, and the sentence then claimed "both are in ..." of a file that had
    in fact been reprojected. The polygons always carry it."""

    features = gpd.GeoDataFrame(geometry=[Point(0.0, 0.0)], crs="EPSG:4326")

    assert "reprojected from EPSG:4326 to the grid's EPSG:2927" in where_they_are(
        features, grid.gdf_vorPolys, source_crs="EPSG:4326"
    )


def test_every_existing_copy_of_that_wording_now_imports_it():
    """`GeoPackageSource._where` (ledger 211), the observation targets' copy (233,
    moved into `observations/_shared.py` by 238) and `_inactive_from_file` (236) were
    THREE copies before this had a fourth caller. Their own suites --
    `test_gpkg_misses.py`, `test_head_targets_on_faces.py`, `test_drn_zones_on_faces.py`,
    `test_grid_edge.py` -- pin the words."""

    from myflopy import geopackage as geopackage_module
    from myflopy.modflow.mf6.observations import _shared as shared_module

    for module in (geopackage_module, shared_module, selection_module, indexing_module):
        assert module.where_they_are is geometry_module.where_they_are


def test_the_grid_edge_refusal_delegates_to_the_shared_wording(grid, tmp_path, monkeypatch):
    """Ledger 236's copy said "the file is in X and the grid in Y" -- which reads
    as a mismatch even when both are the same -- and had its own `_bounds` helper
    next to geometry.py's. Now it says which of the three cases it is."""

    path = _idomain_file(tmp_path / "away.gpkg", [("0", AWAY)], grid.crs)
    seen = []
    original = selection_module.where_they_are
    monkeypatch.setattr(
        selection_module,
        "where_they_are",
        lambda *a, **k: (seen.append(a), original(*a, **k))[1],
    )

    with pytest.raises(ValueError) as raised:
        get_grid_edge_cells(grid, idomain_path=path)

    message = str(raised.value)
    assert len(seen) == 1
    assert "nothing would be inactive" in message
    assert "The features span (500.0, 500.0, 600.0, 600.0)" in message
    assert "the grid (0.0, 0.0, 100.0, 100.0)" in message
    assert "both are in EPSG:2927" in message


def test_selection_no_longer_carries_its_own_extent_helper():
    """Two sibling modules, one importing the other, each with a private
    `_bounds` -- one returning a formatted string and one a tuple -- is the state
    this was in for one commit."""

    assert not hasattr(selection_module, "_bounds")


def test_the_geopackage_refusal_delegates_to_the_shared_wording(tmp_path, monkeypatch):
    """`tests/test_gpkg_misses.py` pins the WORDS; this pins that they come from
    the one copy, so a caller cannot quietly grow its own again."""

    import myflopy as mf
    from myflopy import geopackage as geopackage_module
    from tests.test_gpkg_misses import FAR, _write
    from tests.test_gpkg_on_cell_faces import _context

    seen = []
    original = geopackage_module.where_they_are
    monkeypatch.setattr(
        geopackage_module,
        "where_they_are",
        lambda *a, **k: (seen.append(a), original(*a, **k))[1],
    )

    with pytest.raises(ValueError, match="The features span"):
        mf.drn.gpkg(_write(tmp_path, FAR), layer="bc", context=_context(), nper=1)

    assert len(seen) == 1


def test_the_head_target_refusal_delegates_to_the_shared_wording(monkeypatch):
    """Same, for ledger 233's copy."""

    from myflopy.modflow.mf6.observations import _shared as shared_module
    from tests.test_head_targets_on_faces import _grid as _faces_grid
    from tests.test_head_targets_on_faces import _model, _targets

    seen = []
    original = shared_module.where_they_are
    monkeypatch.setattr(
        shared_module,
        "where_they_are",
        lambda *a, **k: (seen.append(a), original(*a, **k))[1],
    )
    targets = _targets([("W1", Point(500.0, 500.0))])

    with pytest.raises(ValueError, match="The targets span"):
        targets.match_to_model(_model(_faces_grid()))

    assert len(seen) == 1
