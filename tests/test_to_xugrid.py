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


# ---------------------------------------------------------------------------
# `surfaces=True` -- the layer frame on the mesh (ledger 261)
#
# Asked for 2026-09-14: "can I export to xugrid and have all layer elevations
# attached to show in QGIS?" Nothing attached them; `model.to_xugrid()` carried
# the head field alone. The QGIS half is what dictates the SHAPE: a mesh dataset
# there is one value per face, so a `(nlay, ncpl)` array is dropped silently --
# measured against the QGIS 4.2.2 on the development machine, which listed no
# elevation datasets at all for a layered array and every one of them once
# flattened. Hence one variable per layer, never a layer dimension.
# ---------------------------------------------------------------------------


def _surface_names(nlay: int) -> list[str]:
    return ["top", *[f"botm_{n + 1}" for n in range(nlay)]]


@pytest.fixture(scope="module")
def vor_surfaces():
    """A grid with a two-layer frame attached -- no model build or run needed.

    The module's `vor` fixture deliberately has NONE, which is what the
    refusal below is asserted against.
    """

    grid = rectangular_voronoi(CanonicalModelConfig(nrow=5, ncol=5, nlay=2, nper=1))
    ncpl = int(grid.ncpl)
    (
        mf.LayerStack(grid, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("upper", thickness=20.0)
        .add("lower", thickness=30.0)
        .build(attach=True)
    )
    return grid


def test_layer_surface_fields_are_flat_and_one_per_layer(vor_surfaces):
    """`top` plus one `botm_<n>` per layer, each a bare per-cell array."""

    grid = vor_surfaces
    fields = grid.layer_surface_fields()
    nlay = int(grid.gdf_topbtm.shape[1]) - 2
    assert list(fields) == _surface_names(nlay) == ["top", "botm_1", "botm_2"]
    assert all(values.shape == (int(grid.ncpl),) for values in fields.values())
    # read positionally from the geometry-first frame (ledger 196)
    assert np.allclose(fields["top"], grid.gdf_topbtm.iloc[:, 1].to_numpy(dtype=float))
    assert np.allclose(fields["botm_1"], grid.gdf_topbtm.iloc[:, 2].to_numpy(dtype=float))
    # the declared thicknesses, so the columns are not transposed or reversed
    assert np.allclose(fields["top"] - fields["botm_1"], 20.0)
    assert np.allclose(fields["botm_1"] - fields["botm_2"], 30.0)


def test_a_grid_without_surfaces_says_so(vor):
    """A caller asking for surfaces that are not there has a problem to fix,
    not a quietly emptier file. The module's `vor` genuinely has none."""

    assert vor.gdf_topbtm is None
    with pytest.raises(ValueError, match="no layer surfaces"):
        vor.layer_surface_fields()
    with pytest.raises(ValueError, match="no layer surfaces"):
        vor.to_xugrid(surfaces=True)


def test_grid_surfaces_export_is_a_dataset_of_flat_variables(vor_surfaces):
    """`surfaces=True` forces a Dataset: there is more than one variable, and
    none of them carries a layer dimension."""

    nlay = int(vor_surfaces.gdf_topbtm.shape[1]) - 2
    uds = vor_surfaces.to_xugrid(surfaces=True)
    assert isinstance(uds, xu.UgridDataset)
    assert sorted(uds.data_vars) == sorted(_surface_names(nlay))
    face_dim = uds.ugrid.grids[0].face_dimension
    assert all(uds[name].dims == (face_dim,) for name in uds.data_vars)


def test_grid_surfaces_ride_alongside_a_field(vor_surfaces):
    """A single array plus the frame is still one Dataset, not a DataArray."""

    ncpl = int(vor_surfaces.ncpl)
    uds = vor_surfaces.to_xugrid(np.linspace(90.0, 100.0, ncpl), name="head", surfaces=True)
    assert isinstance(uds, xu.UgridDataset)
    assert "head" in uds.data_vars and "top" in uds.data_vars
    # ...and without surfaces the single-array return type is unchanged
    assert isinstance(vor_surfaces.to_xugrid(np.zeros(ncpl), name="head"), xu.UgridDataArray)


def test_an_explicit_value_beats_the_derived_surface(vor_surfaces):
    """What the caller handed over wins over what the grid derived."""

    own = np.full(int(vor_surfaces.ncpl), -999.0)
    uds = vor_surfaces.to_xugrid({"top": own}, surfaces=True)
    assert np.allclose(np.asarray(uds["top"]), own)


@pytest.mark.slow
def test_model_surfaces_export_carries_the_field_and_the_frame(ran_model):
    """`model.to_xugrid(surfaces=True)` is the head field plus flat elevations."""

    model, ncpl = ran_model
    uds = model.to_xugrid(surfaces=True)
    assert isinstance(uds, xu.UgridDataset)
    assert "head" in uds.data_vars
    assert sorted(n for n in uds.data_vars if n != "head") == sorted(_surface_names(2))
    face_dim = uds.ugrid.grids[0].face_dimension
    assert uds["head"].dims == ("time", "layer", face_dim)
    assert uds["top"].dims == (face_dim,)
    # unchanged without it
    assert isinstance(model.to_xugrid(), xu.UgridDataArray)


@pytest.mark.slow
def test_one_isel_flattens_the_field_and_leaves_the_frame_alone(ran_model, tmp_path):
    """The documented QGIS recipe. `isel` cannot touch the elevations -- they
    carry neither dimension -- so one call yields a file whose every variable is
    one value per face."""

    model, ncpl = ran_model
    flat = model.to_xugrid(surfaces=True).isel(time=-1, layer=0)
    face_dim = flat.ugrid.grids[0].face_dimension
    assert all(flat[name].dims == (face_dim,) for name in flat.data_vars)

    path = tmp_path / "model_flat.nc"
    flat.ugrid.to_netcdf(path)
    import xarray as xr

    with xr.open_dataset(path) as ds:
        for name in ["head", *_surface_names(2)]:
            assert name in ds.variables, f"{name} did not survive the round trip"
            assert ds[name].shape == (ncpl,)
