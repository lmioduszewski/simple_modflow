"""``VoronoiGridPlus.to_gpkg`` -- the grid's cells as a GeoPackage QGIS can profile.

Asked for 2026-09-14, after ``to_xugrid(surfaces=True)`` turned out not to solve
the actual problem. A QGIS elevation profile reads a MESH layer's elevation from
its vertex Z, a UGRID file carries none, and no setting on the Elevation tab can
point it at a dataset group -- so a mesh profiles as a flat line at zero however
it is configured (ledger 261/262, both measured against QGIS 4.2.2).

A VECTOR layer profiles from its own geometry Z with one setting, ``Clamping:
Absolute``. Hence this: one GeoPackage layer per surface, each carrying that
surface's elevation as Z.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest

from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)


@pytest.fixture(scope="module")
def bare_grid():
    """A grid with NO layer surfaces attached."""

    return rectangular_voronoi(CanonicalModelConfig(nrow=5, ncol=5, nlay=2, nper=1))


@pytest.fixture(scope="module")
def grid(bare_grid):
    """The same grid with a two-layer frame attached (20 ft over 30 ft)."""

    import myflopy as mf
    from myflopy.layers import Array

    vor = rectangular_voronoi(CanonicalModelConfig(nrow=5, ncol=5, nlay=2, nper=1))
    ncpl = int(vor.ncpl)
    (
        mf.LayerStack(vor, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("upper", thickness=20.0)
        .add("lower", thickness=30.0)
        .build(attach=True)
    )
    return vor


def _layers(path) -> list[str]:
    return list(gpd.list_layers(path)["name"])


def test_the_cells_layer_is_written_with_no_options(grid, tmp_path):
    """The plain case: the grid's polygons and a cell index, nothing else."""

    path = grid.to_gpkg(tmp_path / "grid.gpkg")
    assert path.exists()
    assert _layers(path) == ["cells"]
    cells = gpd.read_file(path, layer="cells")
    assert len(cells) == int(grid.ncpl)
    assert list(cells["cell"]) == list(range(int(grid.ncpl)))
    assert cells.crs == grid.crs
    assert not cells.geometry.iloc[0].has_z


def test_per_cell_data_becomes_attribute_columns(grid, tmp_path):
    ncpl = int(grid.ncpl)
    values = np.linspace(1.0, 2.0, ncpl)
    path = grid.to_gpkg(tmp_path / "grid.gpkg", data={"k": values})
    cells = gpd.read_file(path, layer="cells")
    assert np.allclose(cells["k"].to_numpy(), values)


def test_a_wrong_length_array_is_refused_by_name(grid, tmp_path):
    with pytest.raises(ValueError, match=r"'k' has 3 value\(s\)"):
        grid.to_gpkg(tmp_path / "grid.gpkg", data={"k": [1.0, 2.0, 3.0]})


def test_surfaces_write_one_z_valued_layer_each(grid, tmp_path):
    """The whole point: Z on the geometry, so QGIS needs only Clamping=Absolute."""

    path = grid.to_gpkg(tmp_path / "grid.gpkg", surfaces=True)
    assert _layers(path) == ["cells", "top", "botm_1", "botm_2"]

    for name in ("top", "botm_1", "botm_2"):
        frame = gpd.read_file(path, layer=name)
        assert frame.geometry.iloc[0].has_z, f"{name} carries no Z"
        # every vertex of every cell sits at that cell's own elevation
        first = np.array([geom.exterior.coords[0][2] for geom in frame.geometry])
        assert np.allclose(first, frame[name].to_numpy())
        assert list(frame.columns) == ["cell", name, "geometry"]

    # the surfaces are ordered, so the Z is not scrambled between layers
    top = gpd.read_file(path, layer="top")["top"].to_numpy()
    botm_1 = gpd.read_file(path, layer="botm_1")["botm_1"].to_numpy()
    botm_2 = gpd.read_file(path, layer="botm_2")["botm_2"].to_numpy()
    assert np.allclose(top - botm_1, 20.0)
    assert np.allclose(botm_1 - botm_2, 30.0)


def test_the_cells_layer_carries_the_surface_values_as_columns(grid, tmp_path):
    """Flat, for labelling and queries -- the Z-valued layers are for the profile."""

    path = grid.to_gpkg(tmp_path / "grid.gpkg", surfaces=True, data={"k": np.ones(int(grid.ncpl))})
    cells = gpd.read_file(path, layer="cells")
    assert {"cell", "k", "top", "botm_1", "botm_2"} <= set(cells.columns)
    assert not cells.geometry.iloc[0].has_z


def test_a_grid_without_surfaces_writes_nothing_at_all(bare_grid, tmp_path):
    """Resolved BEFORE the first write, so a refusal does not leave a half-built
    GeoPackage behind -- the file would otherwise hold a `cells` layer and none
    of the surfaces the caller asked for."""

    assert bare_grid.gdf_topbtm is None
    path = tmp_path / "grid.gpkg"
    with pytest.raises(ValueError, match="no layer surfaces"):
        bare_grid.to_gpkg(path, surfaces=True)
    assert not path.exists()
    # ...while the plain export still works on that same grid
    assert bare_grid.to_gpkg(path).exists()


def test_the_layer_name_is_configurable(grid, tmp_path):
    path = grid.to_gpkg(tmp_path / "grid.gpkg", layer="model_cells", surfaces=True)
    assert _layers(path)[0] == "model_cells"


def test_a_missing_parent_directory_is_created(grid, tmp_path):
    path = grid.to_gpkg(tmp_path / "nested" / "deeper" / "grid.gpkg")
    assert path.exists()
