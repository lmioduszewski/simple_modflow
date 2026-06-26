"""Tests for VoronoiGridPlus.to_xugrid() (optional xugrid/xarray export)."""

from __future__ import annotations

import numpy as np
import pytest

# xugrid + xarray are optional dependencies; skip the whole module if absent.
xu = pytest.importorskip("xugrid")

from myflopy.modflow.mf6.canonical_example import (  # noqa: E402
    CanonicalModelConfig,
    rectangular_voronoi,
)


@pytest.fixture(scope="module")
def vor():
    """A small real DISV Voronoi grid (built once for the module)."""
    return rectangular_voronoi(CanonicalModelConfig(nrow=5, ncol=5, nlay=2, nper=1))


def test_topology_matches_grid(vor):
    """The exported mesh has one face per cell, all nodes, and the grid CRS."""
    gp = vor.get_disv_gridprops()
    uds = vor.to_xugrid()
    grid = uds.ugrid.grids[0]
    assert grid.n_face == int(gp["ncpl"])
    assert grid.n_node == int(gp["nvert"])
    assert str(grid.crs).upper().endswith(str(vor.crs).split(":")[-1])


def test_face_areas_match_geopandas(vor):
    """Face geometry is correct: xugrid areas equal the polygon areas."""
    uds = vor.to_xugrid()
    grid = uds.ugrid.grids[0]
    xu_area = np.asarray(grid.area)
    gpd_area = np.asarray(vor.area_list)
    assert np.allclose(xu_area, gpd_area, rtol=1e-9)


def test_single_field_returns_ugriddataarray(vor):
    """A 1-D array becomes a UgridDataArray on the face dimension."""
    ncpl = int(vor.ncpl)
    values = np.linspace(90.0, 100.0, ncpl)
    uda = vor.to_xugrid(values, name="head")
    assert isinstance(uda, xu.UgridDataArray)
    assert uda.name == "head"
    assert uda.dims == (uda.ugrid.grid.face_dimension,)
    assert np.allclose(uda.to_numpy(), values)


def test_dict_and_layered_returns_ugriddataset(vor):
    """A dict with a (nlay, ncpl) array becomes a layered UgridDataset."""
    ncpl = int(vor.ncpl)
    heads = np.linspace(90.0, 100.0, ncpl)
    k = np.vstack([np.full(ncpl, 10.0), np.full(ncpl, 1.0)])
    uds = vor.to_xugrid({"head": heads, "k": k}, layer_dim="layer")
    assert isinstance(uds, xu.UgridDataset)
    assert set(uds.data_vars) == {"head", "k"}
    assert uds["k"].dims == ("layer", uds.ugrid.grids[0].face_dimension)
    assert uds["k"].shape == (2, ncpl)


def test_netcdf_roundtrip(vor, tmp_path):
    """The exported dataset writes to UGRID-NetCDF and reopens intact."""
    ncpl = int(vor.ncpl)
    uds = vor.to_xugrid({"head": np.arange(ncpl, dtype=float)})
    out = tmp_path / "grid.nc"
    uds.ugrid.to_netcdf(out)
    assert out.exists()
    reopened = xu.open_dataset(out)
    assert reopened.ugrid.grids[0].n_face == ncpl
    assert "head" in reopened.data_vars


def test_wrong_length_raises(vor):
    """A cell axis that does not match ncpl is rejected with a clear error."""
    with pytest.raises(ValueError, match="ncpl"):
        vor.to_xugrid(np.zeros(int(vor.ncpl) + 1))
