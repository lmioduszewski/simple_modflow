"""``layer=`` on the grid's feature-to-cell selectors.

A GeoPackage may hold several layers. ``geopandas.read_file`` reads the FIRST one
and merely emits a ``UserWarning``, so before this parameter existed there was no
way to ask for any other -- and picking up the wrong table was silent in any code
that does not surface warnings.

``layer=None`` is GDAL's own default, so every call written before the parameter
behaves exactly as it did; that equivalence is what ``test_default_is_unchanged``
pins.
"""

from __future__ import annotations

import geopandas as gpd
import pytest
import shapely as shp

from myflopy.modflow.mf6.canonical import irregular_voronoi_grid


@pytest.fixture(scope="module")
def vor():
    """A small deterministic Voronoi grid."""

    return irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0)


@pytest.fixture(scope="module")
def multi_layer(tmp_path_factory, vor):
    """A GeoPackage with two layers whose footprints do not overlap."""

    x0, y0, x1, y1 = vor.gdf_vorPolys.total_bounds
    third = (x1 - x0) / 3.0
    path = tmp_path_factory.mktemp("vectors") / "multi.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"]},
        geometry=[shp.box(x0, y0, x0 + third, y1)],
        crs=vor.crs,
    ).to_file(path, layer="ponds", driver="GPKG")
    gpd.GeoDataFrame(
        {"name": ["east"]},
        geometry=[shp.box(x1 - third, y0, x1, y1)],
        crs=vor.crs,
    ).to_file(path, layer="basins", driver="GPKG")
    return path


@pytest.fixture(scope="module")
def single_layer(tmp_path_factory, vor):
    """A one-layer GeoPackage, the shape most callers pass."""

    x0, y0, x1, y1 = vor.gdf_vorPolys.total_bounds
    path = tmp_path_factory.mktemp("vectors") / "one.gpkg"
    gpd.GeoDataFrame(
        {"name": ["all"]},
        geometry=[shp.box(x0, y0, x1, y1)],
        crs=vor.crs,
    ).to_file(path, layer="only", driver="GPKG")
    return path


def test_default_is_unchanged(vor, multi_layer):
    """Omitting ``layer`` still reads the first layer -- old calls are untouched."""

    assert vor.get_vor_cells_as_dict(multi_layer) == vor.get_vor_cells_as_dict(
        multi_layer, layer="ponds"
    )


def test_layer_selects_a_later_layer(vor, multi_layer):
    """The second layer is reachable, and covers different cells from the first."""

    first = vor.get_vor_cells_as_dict(multi_layer, layer="ponds")
    second = vor.get_vor_cells_as_dict(multi_layer, layer="basins")
    assert first and second
    # values are nested one level: {name: [[cell, cell, ...]]}, one inner list
    # per matched geometry. That shape predates this change and two callers rely
    # on it, so it is flattened here rather than altered.
    west = {cell for group in first.values() for cells in group for cell in cells}
    east = {cell for group in second.values() for cells in group for cell in cells}
    assert west and east
    assert not west & east, "the two layers should select disjoint cells"


def test_layer_accepts_an_index(vor, multi_layer):
    """An integer indexes the layers in file order."""

    assert vor.get_vor_cells_as_dict(multi_layer, layer=1) == vor.get_vor_cells_as_dict(
        multi_layer, layer="basins"
    )


def test_unknown_layer_names_the_available_ones(vor, multi_layer):
    """GDAL says only 'could not be opened'; say what the file actually holds."""

    with pytest.raises(ValueError, match="no layer 'lakes'") as caught:
        vor.get_vor_cells_as_dict(multi_layer, layer="lakes")
    assert "'ponds'" in str(caught.value) and "'basins'" in str(caught.value)


def test_out_of_range_index_is_refused(vor, multi_layer):
    """An index past the end reports the count rather than failing inside GDAL."""

    with pytest.raises(ValueError, match="index 5 is out of range"):
        vor.get_vor_cells_as_dict(multi_layer, layer=5)


def test_single_layer_file_still_works_by_name(vor, single_layer):
    """Naming the only layer is legal, not a special case."""

    assert vor.get_vor_cells_as_dict(single_layer, layer="only")


def test_series_default_is_unchanged(vor, multi_layer):
    """The series form has the same default and the same equivalence."""

    plain = vor.get_vor_cells_as_series(multi_layer, name_field="name")
    named = vor.get_vor_cells_as_series(multi_layer, name_field="name", layer="ponds")
    assert list(plain.index) == list(named.index) == ["west"]
    assert plain.to_list() == named.to_list()


def test_series_layer_selects_a_later_layer(vor, multi_layer):
    """``layer=`` reaches the second layer through the series form too."""

    series = vor.get_vor_cells_as_series(multi_layer, name_field="name", layer="basins")
    assert list(series.index) == ["east"]


def test_series_refuses_layer_for_an_in_memory_geometry(vor):
    """A geometry has no layers to choose between, so the parameter cannot act."""

    with pytest.raises(ValueError, match="layer= selects a layer inside a vector FILE"):
        vor.get_vor_cells_as_series(shp.box(0, 0, 10, 10), layer="ponds")


def test_series_unknown_layer_names_the_available_ones(vor, multi_layer):
    """Same message from the series form -- one helper does the reading."""

    with pytest.raises(ValueError, match="no layer 'lakes'"):
        vor.get_vor_cells_as_series(multi_layer, layer="lakes")
