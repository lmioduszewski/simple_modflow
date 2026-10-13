"""An edit's override REPLACES a recipe that fills the same option (ledger 561).

``fork_package(mf.npf(k=mf.Polygons(...)), k=array)`` -- and ``packages.fork``/``edit``
and ``with_options`` -- kept ``k`` in both ``options`` and ``inputs``, which
``PackageSpec.build`` refused as "supplied twice"; over a product or a per-layer list of
recipes the edit itself was refused. The override now replaces the recipe, named in the
stamp's ``replaced_inputs``. A recipe the builder REQUIRES (a lazy ``.gpkg``
declaration's records) is not a source for an option and is still refused.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, box

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.specs import fork_package


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    return grid, grid.get_disv_gridprops()


def test_a_fork_override_replaces_an_input_recipe(tmp_path):
    """``fork_package(npf(k=mf.Polygons(...)), k=array)`` built 'supplied twice'."""

    grid, props = _grid()
    path = tmp_path / "units.gpkg"
    gpd.GeoDataFrame({"k": [3.0]}, geometry=[box(*grid.gdf_vorPolys.total_bounds)],
                     crs=grid.crs).to_file(path, layer="u")
    npf = mf.npf(k=mf.Polygons(path, column="k", layer="u", default=1.0), icelltype=0)
    forked = fork_package(npf, k=np.full(grid.ncpl, 7.0))
    assert "k" not in forked.inputs
    assert forked.metadata["forked_from"]["replaced_inputs"] == ["k"]
    simulation = mf.simulation(mf.gwf("m", context=mf.ModelContext(grid=grid), packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10.0, botm=[0.0]),
        forked, mf.ic(strt=5.0)]), name="s",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    built = simulation.build_flopy(tmp_path / "built")
    assert np.allclose(built.models["m"].model.npf.k.array, 7.0)
    # The same rule on with_options; a per-layer list of recipes is replaced too.
    layered = mf.npf(k=[mf.Polygons(path, column="k", layer="u", default=1.0), 1.0])
    assert "k" not in layered.with_options(k=[2.0, 1.0]).inputs
    assert "k" not in fork_package(layered, k=[2.0, 1.0]).inputs


def test_a_recipe_the_builder_requires_is_still_refused(tmp_path):
    """A lazy ``.gpkg`` declaration's records are its builder's input, not an option's source."""

    grid, _props = _grid()
    path = tmp_path / "b.gpkg"
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame({"head": [1.0]}, geometry=[LineString([(xmin + 5, ymin), (xmin + 5, ymax)])],
                     crs=grid.crs).to_file(path, layer="b")
    lazy = mf.chd.gpkg(path, layer="b", nper=1, head="head")
    with pytest.raises(TypeError, match="part of this package's recipe"):
        fork_package(lazy, stress_period_data={0: [[(0, 0), 1.0]]})
