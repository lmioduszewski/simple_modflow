"""Grid inputs are read by one rule each and refused by name (ledgers 457, 458, 459).

- an ``idomain_path`` given as a ``str`` failed an ``assert`` on the grid and on the
  legacy model;
- ``generate_dissipating_point_cloud`` checked its polygon with a bare ``assert``;
- ``GridSpec.voronoi(min_angle=, triangle_options={"angle": ...})`` and a spec
  holding ``options["min_angle"]`` beside ``triangle_options["angle"]`` meshed at
  whichever one won, in silence;
- three copies of "a finite number above zero" disagreed about ``"2.5"`` and ``True``;
- the legacy ``DisuGrid`` read ``gdf_topbtm["top"]``, which a ``LayerStack``-attached
  table (integer labels) does not have, and the ``KeyError`` escaped its handler.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
import shapely as shp

import myflopy as mf
from myflopy._numbers import positive_number
from myflopy.modflow.mf6.grid.geometry import (
    generate_grid_around_point,
    generate_grid_polygons,
)
from myflopy.modflow.mf6.grid.triangle import TriangleGrid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _two_cells(**kwargs) -> VoronoiGridPlus:
    verts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]])
    return VoronoiGridPlus(
        verts=verts, iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]]), **kwargs,
    )


# --- idomain_path ------------------------------------------------------------ #


@pytest.mark.parametrize("given", ["inactive.gpkg", Path("inactive.gpkg")])
def test_a_grid_takes_an_idomain_path_as_str_or_path(given):
    grid = _two_cells(idomain_path=given)
    assert grid.idomain_path == Path("inactive.gpkg")
    assert type(grid.idomain_path) is not str


def test_a_grid_refuses_an_idomain_path_that_is_not_a_path():
    with pytest.raises(TypeError, match=r"idomain_path must be a path .*got 3 \(int\)"):
        _two_cells(idomain_path=3)


def test_the_legacy_model_takes_a_str_idomain_path_and_clears_on_none(tmp_path):
    from myflopy.modflow.mf6.simulation.base import SimulationBase

    model = SimulationBase(name="m", mf_folder_path=tmp_path, vor=_two_cells(), nper=1,
                           idomain_path=str(tmp_path / "off.gpkg"))
    assert model.idomain_path == tmp_path / "off.gpkg"
    model.idomain_path = None
    assert model.idomain_path is None
    with pytest.raises(TypeError, match="idomain_path must be a path"):
        model.idomain_path = 1.5


# --- the dissipating point cloud --------------------------------------------- #


def test_the_dissipating_cloud_refuses_a_multipolygon_by_name(tmp_path):
    parts = shp.MultiPolygon([shp.box(0, 0, 1, 1), shp.box(5, 5, 6, 6)])
    with pytest.raises(TypeError, match=r"polygon= must be ONE shapely Polygon.*MultiPolygon of 2 parts"):
        TriangleGrid.generate_dissipating_point_cloud(parts)
    path = tmp_path / "two.gpkg"
    gpd.GeoDataFrame(geometry=list(parts.geoms), crs="EPSG:2927").to_file(path)
    with pytest.raises(TypeError, match=r"MultiPolygon of 2 parts \(the union of '.*two.gpkg'\)"):
        TriangleGrid.generate_dissipating_point_cloud(str(path))
    with pytest.raises(TypeError, match=r"got a NoneType"):
        TriangleGrid.generate_dissipating_point_cloud()
    rings = TriangleGrid.generate_dissipating_point_cloud(
        shp.box(0, 0, 1000, 1000), buffer_dist=500, num_buffers=2, min_spacing=100,
        max_spacing=200,
    )
    assert len(rings) == 3


# --- two spellings of one grid option in two places ------------------------- #


@pytest.fixture
def boundary():
    return mf.GeoPackageSourceSpec("boundary.gpkg")


@pytest.mark.parametrize(("kwargs", "named"), [
    ({"min_angle": 30, "triangle_options": {"angle": 20}},
     r"min_angle= and triangle_options\['angle'\]"),
    ({"options": {"min_angle": 30}, "triangle_options": {"angle": 20}},
     r"options\['min_angle'\] and triangle_options\['angle'\]"),
    ({"angle": 30, "options": {"min_angle": 20}}, r"angle= and options\['min_angle'\]"),
    ({"min_angle": 30, "options": {"angle": 20}}, r"min_angle= and options\['angle'\]"),
    ({"options": {"model_ws": "a"}, "triangle_options": {"workspace": "b"}},
     r"options\['model_ws'\] and triangle_options\['workspace'\]"),
    ({"triangle_options": {"model_ws": "a", "workspace": "b"}},
     r"triangle_options\['model_ws'\] and triangle_options\['workspace'\]"),
])
def test_gridspec_voronoi_refuses_one_option_spelt_twice_in_two_places(boundary, kwargs, named):
    with pytest.raises(TypeError, match=named + r" are spellings of ONE option"):
        mf.GridSpec.voronoi(boundary=boundary, **kwargs)


def test_one_key_in_two_places_keeps_its_documented_precedence(boundary):
    """A keyword wins over the dict it lands in; that is documented, not refused."""

    spec = mf.GridSpec.voronoi(boundary=boundary, angle=30, triangle_options={"angle": 20})
    assert spec.triangle_options["angle"] == 30


def test_a_saved_spec_holding_both_spellings_is_refused_at_resolve(boundary, tmp_path):
    import dataclasses

    spec = mf.GridSpec.voronoi(boundary=boundary, min_angle=30)
    saved = dataclasses.replace(spec, options={"min_angle": 25},
                                triangle_options={"angle": 20})
    with pytest.raises(TypeError, match=r"options\['min_angle'\] and triangle_options\['angle'\]"):
        saved.resolve(workspace=tmp_path)


# --- one positive-number check ---------------------------------------------- #


@pytest.mark.parametrize("bad", ["2.5", True, np.bool_(True), 0, -1.0, np.nan, np.inf, None])
def test_every_grid_step_refuses_what_is_not_a_number_above_zero(bad):
    with pytest.raises(ValueError, match="must be a finite number"):
        positive_number(bad, "spacing")
    with pytest.raises(ValueError, match=r"spacing must be a finite number above zero"):
        generate_grid_around_point(shp.Point(0, 0), bad, 2, "EPSG:2927")


@pytest.mark.parametrize("bad", ["2", True, np.bool_(True), -1, np.nan, None])
def test_a_lattice_size_is_read_by_the_same_rule(bad):
    """``float(size)`` took ``"2"`` as 2 and ``True`` as 1 (review of ledger 457)."""

    with pytest.raises(ValueError, match=r"size must be a finite number of steps, zero or more"):
        generate_grid_around_point(shp.Point(0, 0), 1.0, bad, "EPSG:2927")
    with pytest.raises(ValueError, match=r"size must be a finite number of steps"):
        generate_grid_polygons(shp.Point(0, 0), 1.0, bad, 0.0)


@pytest.mark.parametrize("bad", ["0.5", True, np.bool_(False), None])
def test_a_box_gap_is_read_by_the_same_rule(bad):
    with pytest.raises(ValueError, match=r"gap must be a finite number, zero or more"):
        generate_grid_polygons(shp.Point(0, 0), 2.5, 4, bad)


def test_zero_is_a_size_and_a_gap_but_not_a_step():
    assert len(generate_grid_around_point(shp.Point(0, 0), 1.0, 0, "EPSG:2927")) == 1
    assert len(generate_grid_polygons(shp.Point(0, 0), 2.5, np.int64(2), 0).geoms) == 4
    assert positive_number(0, "gap", allow_zero=True) == 0.0
    with pytest.raises(ValueError, match="gap must be a finite number, zero or more, got -1"):
        positive_number(-1, "gap", allow_zero=True)
    with pytest.raises(ValueError, match="spacing must be a finite number above zero, got 0"):
        positive_number(0, "spacing")


def test_a_lattice_takes_a_numpy_number():
    points = generate_grid_around_point(shp.Point(0, 0), np.float32(2.5), 2, "EPSG:2927")
    assert len(points) == 9
    assert positive_number(np.int64(3), "size") == 3.0


# --- the legacy DisuGrid reads gdf_topbtm by position ------------------------ #


def _disu(tmp_path, table):
    from myflopy.modflow.mf6.simulation.base import SimulationBase
    from myflopy.modflow.mf6.simulation.discretization import DisuGrid

    vor = _two_cells()
    vor.get_disu_connectivity(validate=True)
    vor.gdf_topbtm = table(vor)
    model = SimulationBase(name="disu", mf_folder_path=tmp_path, vor=vor, nper=1)
    return DisuGrid(vor=vor, model=model).disu


def _table(columns):
    def make(vor):
        frame = gpd.GeoDataFrame(geometry=vor.gdf_vorPolys.geometry.values, crs=vor.crs)
        for label, values in columns.items():
            frame[label] = values
        return frame

    return make


@pytest.mark.parametrize("labels", [(0, 1), ("top", "bottom"), ("ground", "base")])
def test_disugrid_reads_top_and_bottom_whatever_the_columns_are_called(tmp_path, labels):
    disu = _disu(tmp_path, _table({labels[0]: [10.0, 11.0], labels[1]: [0.0, 1.0]}))
    assert disu.top.get_data().tolist() == [10.0, 11.0]
    assert disu.bot.get_data().tolist() == [0.0, 1.0]


def test_disugrid_refuses_to_pick_one_of_several_bottoms(tmp_path):
    with pytest.raises(ValueError, match=r"DisuGrid is ONE layer.*2 layer bottoms.*bottom="):
        _disu(tmp_path, _table({0: [10.0, 11.0], 1: [5.0, 6.0], 2: [0.0, 1.0]}))
