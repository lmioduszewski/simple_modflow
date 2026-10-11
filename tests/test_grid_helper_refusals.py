"""Grid helpers that crashed, dropped an argument, or contradicted their annotation.

From the 2026-09-27 docstring pass, which documented each as current behaviour:

* ``set_k_vor()`` with no mapping raised ``AttributeError`` (``k_dict.keys()`` on
  ``None``) instead of giving every cell ``k_default``; and its ``"Kh"`` column was
  seeded with the int ``100``, so a float K took pandas' incompatible-dtype upcast,
  which warns today and is deprecated to raise.
* ``get_vor_idx_from_geometry_idx()`` with no query raised ``AttributeError``
  (``None.iloc``) before the query check its sibling raises could run.
* ``generate_grid_around_point`` / ``generate_grid_polygons`` were annotated
  ``spacing: float`` and drove ``range()`` with it, so ``2.5`` raised ``TypeError``;
  and ``int()`` of each end truncated toward zero, so a lattice around a
  non-integer centre was shifted and lost a row.
* ``VoronoiGridPlus(crs=None)`` was refused only after the triangulation's cells
  were built, by the ``grid.crs =`` setter's message.
* ``VoronoiGridPlus(idomain=...)`` -- and ``vor_from_disu``/``from_gsf``'s, and
  ``GridSpec.voronoi(idomain=...)`` routed to it -- was accepted and dropped in
  silence (ledger 388); the array form also dropped ``qhull_options`` and any
  other keyword.

Every expected lattice here is written out by hand, not computed by the code under
test.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy.grid_spec_resolver import _voronoi_options
from myflopy.modflow.mf6.grid import geometry
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

VERTS = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
IVERTS = [[0, 3, 4, 1], [1, 4, 5, 2]]
XCYC = np.array([[0.5, 0.5], [1.5, 0.5]])


def _grid(**kwargs) -> VoronoiGridPlus:
    return VoronoiGridPlus(verts=VERTS, iverts=IVERTS, xcyc=XCYC, **kwargs)


# --- set_k_vor ---------------------------------------------------------------


def test_set_k_vor_with_no_mapping_gives_every_cell_the_default():
    grid = _grid()
    assert grid.set_k_vor() == [100.0, 100.0]
    assert grid.set_k_vor(k_default=3.5) == [3.5, 3.5]
    assert geometry.set_k_vor(grid, None, k_default=7) == [7.0, 7.0]


def test_set_k_vor_writes_a_float_k_over_the_int_default_without_an_upcast():
    grid = _grid()
    with warnings.catch_warnings():
        warnings.simplefilter("error")      # pandas' FutureWarning would raise here
        k = grid.set_k_vor({5.5: [0]}, k_default=100)
    assert k == [5.5, 100.0]
    assert all(type(value) is float for value in k)
    assert grid.gdf_vorPolys["Kh"].dtype == np.float64


# --- get_vor_idx_from_geometry_idx ------------------------------------------


@pytest.mark.parametrize("call", ["free", "method"])
def test_get_vor_idx_from_geometry_idx_with_no_query_names_the_missing_frame(call):
    grid = _grid()
    with pytest.raises(ValueError, match="no query was given: pass the features as gdf_to_query"):
        if call == "free":
            geometry.get_vor_idx_from_geometry_idx(grid)
        else:
            grid.get_vor_idx_from_geometry_idx()


# --- the point / box lattice -------------------------------------------------


def _axis(series, attr):
    return sorted({round(getattr(point, attr), 9) for point in series})


def test_a_fractional_spacing_builds_the_lattice_it_describes():
    points = geometry.generate_grid_around_point(Point(0, 0), 2.5, 4, "EPSG:2927")
    assert len(points) == 25
    assert _axis(points, "x") == [-5.0, -2.5, 0.0, 2.5, 5.0]
    assert _axis(points, "y") == [-5.0, -2.5, 0.0, 2.5, 5.0]
    assert str(points.crs) == "EPSG:2927"

    boxes = geometry.generate_grid_polygons(Point(0, 0), 2.5, 4, 0.1)
    assert len(boxes.geoms) == 16
    assert boxes.bounds == pytest.approx((-4.95, -4.95, 4.95, 4.95))
    assert {round(b.area, 9) for b in boxes.geoms} == {round(2.4 ** 2, 9)}

    # The VoronoiGridPlus spellings are the same lattice.
    assert VoronoiGridPlus.generate_grid_around_point(Point(0, 0), 2.5, 4, 2927).equals(points)
    assert VoronoiGridPlus.generate_grid_polygons(Point(0, 0), 2.5, 4, 0.1).equals(boxes)


def test_an_integer_lattice_on_an_integer_centre_is_unchanged():
    # What range(int(minx), int(maxx) + 1, spacing) gave, which callers rely on.
    points = geometry.generate_grid_around_point(Point(10, 20), 2, 4, "EPSG:2927")
    assert _axis(points, "x") == [6.0, 8.0, 10.0, 12.0, 14.0]
    assert _axis(points, "y") == [16.0, 18.0, 20.0, 22.0, 24.0]


def test_the_lattice_is_centred_on_a_non_integer_centre():
    # int() truncated -3.3 to -3 and 4.7 to 4: x = -3, -1, 1, 3 -- shifted and a row short.
    points = geometry.generate_grid_around_point(Point(0.7, 0.3), 2, 4, "EPSG:2927")
    assert _axis(points, "x") == [-3.3, -1.3, 0.7, 2.7, 4.7]
    assert _axis(points, "y") == [-3.7, -1.7, 0.3, 2.3, 4.3]


def test_an_odd_size_keeps_both_ends():
    # range(int(-1.5), int(1.5) + 1, 1) gave -1, 0, 1: three points for three steps.
    points = geometry.generate_grid_around_point(Point(0, 0), 1, 3, "EPSG:2927")
    assert _axis(points, "x") == [-1.5, -0.5, 0.5, 1.5]


@pytest.mark.parametrize(
    "spacing, size, match",
    [(0, 4, "spacing"), (-1.0, 4, "spacing"), (float("nan"), 4, "spacing"), (1.0, -2, "size")],
)
def test_a_lattice_that_cannot_exist_is_refused(spacing, size, match):
    with pytest.raises(ValueError, match=match):
        geometry.generate_grid_around_point(Point(0, 0), spacing, size, "EPSG:2927")
    with pytest.raises(ValueError, match=match):
        geometry.generate_grid_polygons(Point(0, 0), spacing, size, 0.0)


@pytest.mark.parametrize("gap", [-0.1, 2.5, 3.0, float("nan"), float("inf")])
def test_a_gap_that_leaves_no_box_or_overlaps_is_refused(gap):
    # Ledger 415: shp.box accepts min > max, so a gap of spacing or more
    # gave inverted boxes and a negative one overlapping boxes, in silence.
    with pytest.raises(ValueError, match=r"gap must be a finite number, zero or more and below spacing"):
        geometry.generate_grid_polygons(Point(0, 0), 2.5, 4, gap)
    with pytest.raises(ValueError, match=r"gap must be"):
        VoronoiGridPlus.generate_grid_polygons(Point(0, 0), 2.5, 4, gap)


@pytest.mark.parametrize("size", [0, 0.5, 0.999])
def test_a_size_under_one_step_holds_no_box_and_is_refused(size):
    # It returned an empty MultiPolygon.
    with pytest.raises(ValueError, match=r"is under one step, so the lattice holds no box"):
        geometry.generate_grid_polygons(Point(0, 0), 2.5, size, 0.0)
    assert len(geometry.generate_grid_polygons(Point(0, 0), 2.5, 1, 0.0).geoms) == 1


# --- generate_grid_coordinates --------------------------------------------------


def _extent_grid(x0, x1, y0, y1) -> VoronoiGridPlus:
    verts = np.array([[x0, y0], [x1, y0], [x0, y1], [x1, y1]], dtype=float)
    return VoronoiGridPlus(
        verts=verts, iverts=[[0, 2, 3, 1]], xcyc=np.array([[(x0 + x1) / 2, (y0 + y1) / 2]])
    )


def test_grid_coordinates_take_a_fractional_spacing():
    # Ledger 415: range() refused 2.5 with TypeError.
    xs, ys = geometry.generate_grid_coordinates(_extent_grid(0.0, 5.0, 0.0, 2.5), 2.5)
    assert sorted(set(xs)) == [0.0, 2.5, 5.0]
    assert sorted(set(ys)) == [0.0, 2.5]
    assert len(xs) == len(ys) == 6


def test_grid_coordinates_start_at_a_non_integer_extent_not_its_truncation():
    # int() truncated 0.7 to 0 and 4.6 to 4: x = 0, 2, 4, the last short of the edge.
    grid = _extent_grid(0.7, 4.6, -1.3, 0.2)
    xs, ys = grid.generate_grid_coordinates(2)
    assert sorted({round(x, 9) for x in xs}) == [0.7, 2.7, 4.7]
    assert sorted({round(y, 9) for y in ys}) == [-1.3, 0.7]


def test_grid_coordinates_on_an_integer_extent_are_unchanged():
    # What range(int(x_min), int(x_max) + s, s) gave.
    assert _grid().generate_grid_coordinates(1) == ([0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1])
    xs, ys = _extent_grid(0.0, 10.0, 0.0, 10.0).generate_grid_coordinates(4)
    assert sorted(set(xs)) == [0.0, 4.0, 8.0, 12.0]


@pytest.mark.parametrize("spacing", [0, -1.0, float("nan"), True, "2"])
def test_grid_coordinates_refuse_a_spacing_that_is_not_a_positive_number(spacing):
    with pytest.raises(ValueError, match=r"grid_spacing must be a finite number above zero"):
        _grid().generate_grid_coordinates(spacing)


# --- set_k_vor refuses an unknown cell -------------------------------------------


def test_set_k_vor_refuses_a_cell_id_the_grid_does_not_have_before_writing():
    # Ledger 415: index.isin matched only the ids the grid has, so a typo'd id
    # left its cell at k_default in silence.
    grid = _grid()
    grid.set_k_vor({5.0: [1]})
    with pytest.raises(ValueError, match=r"k_dict names cell ids the grid does not have \(K 9.0: \[7\]\)"):
        grid.set_k_vor({1.0: [0], 9.0: [7]})
    assert grid.gdf_vorPolys["Kh"].tolist() == [100.0, 5.0]


# --- VoronoiGridPlus(crs=None) -------------------------------------------------


def test_crs_none_is_refused_at_construction_naming_crs_before_any_cell_is_built():
    with pytest.raises(ValueError, match=r"VoronoiGridPlus\(crs=None\): a grid needs"):
        _grid(crs=None)
    # Before the triangulation is read: a tri FloPy would refuse never gets there.
    with pytest.raises(ValueError, match=r"VoronoiGridPlus\(crs=None\): a grid needs"):
        VoronoiGridPlus(tri=object(), crs=None)
    with pytest.raises(ValueError, match=r"VoronoiGridPlus\(crs='not a crs'\): pyproj"):
        _grid(crs="not a crs")


def test_the_default_crs_is_never_refused():
    assert _grid().crs == "EPSG:2927"
    assert _grid(crs=2926).crs == "EPSG:2926"


# --- idomain= is refused everywhere it was dropped ----------------------------


def test_voronoigridplus_refuses_idomain_naming_the_model_context():
    with pytest.raises(TypeError, match=r"ModelContext\(grid=\.\.\., domain=\.\.\.\)"):
        _grid(idomain=[0])
    assert not hasattr(_grid(), "idomain")


@pytest.mark.parametrize("builder", ["vor_from_disu", "from_gsf"])
def test_the_array_classmethods_take_no_idomain_or_qhull_options(builder):
    import inspect

    params = inspect.signature(getattr(VoronoiGridPlus, builder)).parameters
    assert "idomain" not in params
    assert "qhull_options" not in params
    assert "idomain_path" in params


def test_the_array_form_refuses_options_it_cannot_use():
    with pytest.raises(TypeError, match="qhull_options with the array form"):
        _grid(qhull_options="Qbb")
    with pytest.raises(TypeError, match="furthest_site with the array form"):
        _grid(furthest_site=True)


@pytest.mark.parametrize(
    "kwargs, spelling",
    [
        ({"idomain": [0]}, "idomain="),
        ({"options": {"idomain": [0]}}, "options={'idomain': ...}"),
        ({"voronoi_options": {"idomain": [0]}}, "voronoi_options={'idomain': ...}"),
    ],
)
def test_gridspec_voronoi_refuses_idomain_naming_the_model_context(tmp_path, kwargs, spelling):
    with pytest.raises(TypeError, match=r"takes no .*ModelContext\(domain=\.\.\.\)") as raised:
        mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=tmp_path / "b.gpkg"), boundary_max_area=10.0, **kwargs)
    assert spelling in str(raised.value)


def test_a_spec_carrying_idomain_reaches_the_grids_refusal_at_resolve(tmp_path):
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=tmp_path / "b.gpkg"), boundary_max_area=10.0)
    payload = spec.to_dict()
    payload["voronoi_options"] = {"idomain": [0]}
    restored = mf.GridSpec.from_dict(payload)
    options = _voronoi_options(restored, dict(restored.options), project_root=None)
    assert options["idomain"] == [0]
    with pytest.raises(TypeError, match="takes no idomain="):
        VoronoiGridPlus(**options)
