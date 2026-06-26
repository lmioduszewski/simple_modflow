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


def test_ugrid2d_topology(vor):
    """The shared ugrid2d() helper returns a valid mesh with the grid CRS."""
    gp = vor.get_disv_gridprops()
    grid = vor.ugrid2d()
    assert grid.n_face == int(gp["ncpl"])
    assert grid.n_node == int(gp["nvert"])
    assert str(grid.crs).upper().endswith(str(vor.crs).split(":")[-1])


def test_plot_requires_single_field(vor):
    """The documented reduce-then-plot rule: isel a (time, layer, cell) array
    down to one field to plot; plotting the full array raises."""
    import matplotlib
    matplotlib.use("Agg")
    import xarray as xr

    ncpl = int(vor.ncpl)
    grid = vor.ugrid2d()
    face_dim = grid.face_dimension
    uda = xu.UgridDataArray(
        xr.DataArray(
            np.random.rand(1, 2, ncpl),
            dims=("time", "layer", face_dim),
            coords={"time": [0], "layer": [0, 1]},
            name="head",
        ),
        grid,
    )
    # Reduced to one (cell,) field -> plots.
    art = uda.isel(time=-1, layer=0).ugrid.plot()
    assert art is not None
    # The full (time, layer, cell) array cannot be plotted directly.
    with pytest.raises(ValueError, match="non-topology dimensions"):
        uda.ugrid.plot()


# --------------------------------------------------------------------------- #
# model.to_xugrid() -- needs a real MF6 run, so these are slow.
# --------------------------------------------------------------------------- #

import myflopy as mf  # noqa: E402
from myflopy.layers import Array  # noqa: E402
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig  # noqa: E402


def _build_two_layer_run(root):
    """Build + run a tiny 2-layer steady model with a CHD gradient."""
    grid = rectangular_voronoi(CanonicalModelConfig(nrow=6, ncol=6, nlay=2, nper=1))
    gp = grid.get_disv_gridprops()
    ncpl = int(grid.ncpl)
    layers = (
        mf.LayerStack(grid, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("upper", thickness=20.0)
        .add("lower", thickness=30.0)
        .build(attach=True)
    )
    ctx = mf.ModelContext(grid=grid, domain=layers.idomain, surfaces=grid.gdf_topbtm)
    chd = {0: [[(0, 0), 95.0], [(0, ncpl - 1), 90.0]]}
    flow = mf.gwf("hello", context=ctx, packages=[
        mf.disv(nlay=2, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                vertices=gp["vertices"], cell2d=gp["cell2d"],
                top=layers.top, botm=layers.botm, idomain=layers.idomain,
                length_units="FEET"),
        mf.ic(strt=92.0), mf.npf(k=10.0), mf.sto(steady_state={0: True}),
        mf.chd(stress_period_data=chd),
        mf.oc(head_filerecord="hello.hds", budget_filerecord="hello.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    sim = mf.SimulationSpec("baseline", models=(flow,), packages=(
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=("hello",)),
    ))
    project = mf.Project(root)
    project.add_simulation(sim)
    run = project.run("baseline", "baseline")
    assert run.success, "tiny model failed to converge"
    return run.model("hello"), ncpl


@pytest.fixture(scope="module")
def ran_model(tmp_path_factory):
    root = tmp_path_factory.mktemp("model_xugrid")
    model, ncpl = _build_two_layer_run(root)
    return model, ncpl


@pytest.mark.slow
def test_model_to_xugrid_stacks_time_layer_cell(ran_model):
    """model.to_xugrid() returns a (time, layer, cell) UgridDataArray."""
    model, ncpl = ran_model
    uda = model.to_xugrid()
    assert isinstance(uda, xu.UgridDataArray)
    face_dim = uda.ugrid.grid.face_dimension
    assert uda.dims == ("time", "layer", face_dim)
    assert uda.shape == (1, 2, ncpl)
    assert {"time", "kstp", "kper", "layer"}.issubset(set(uda.coords))
    assert list(uda["layer"].values) == [0, 1]


@pytest.mark.slow
def test_model_to_xugrid_matches_hds_array(ran_model):
    """Each layer of the export equals model.hds.array for that layer."""
    model, _ = ran_model
    uda = model.to_xugrid()
    for layer in (0, 1):
        exported = uda.isel(time=-1, layer=layer).to_numpy()
        reference = model.hds.array(layer=layer)
        assert np.allclose(exported, reference, equal_nan=True)


@pytest.mark.slow
def test_model_to_xugrid_layer_subset(ran_model):
    """layers=[...] restricts the export to the requested layers."""
    model, ncpl = ran_model
    uda = model.to_xugrid(layers=[1])
    assert uda.shape == (1, 1, ncpl)
    assert list(uda["layer"].values) == [1]


@pytest.mark.slow
def test_model_to_xugrid_netcdf_roundtrip(ran_model, tmp_path):
    """The exported heads write to UGRID-NetCDF and reopen intact."""
    model, ncpl = ran_model
    out = tmp_path / "heads.nc"
    model.to_xugrid().isel(time=-1).ugrid.to_netcdf(out)
    assert out.exists()
    reopened = xu.open_dataset(out)
    assert reopened.ugrid.grids[0].n_face == ncpl


@pytest.mark.slow
def test_model_to_xugrid_bad_time_raises(ran_model):
    """An unavailable kstpkper key is rejected with a clear error."""
    model, _ = ran_model
    with pytest.raises(ValueError, match="unavailable"):
        model.to_xugrid(times=[(99, 99)])
