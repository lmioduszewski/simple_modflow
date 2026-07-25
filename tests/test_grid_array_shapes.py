import os

os.environ["MPLBACKEND"] = "Agg"

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, box

from myflopy.modflow.mf6.contour_plotting import _resolve_contour_levels, contour_line_segments
from myflopy.modflow.mf6.headsplus import _as_layer_cell_heads
from myflopy.modflow.utils.datatypes.choros import _as_cell_vector


class _DummyVor:
    def __init__(self, gdf):
        self.gdf_vorPolys = gdf
        self.ncpl = len(gdf)
        centers = gdf.geometry.centroid
        self.centroids_x = centers.x.to_list()
        self.centroids_y = centers.y.to_list()


def test_heads_shape_helper_preserves_disv_layer_cell_arrays():
    values = np.array(
        [
            [10.0, 9.5, 9.0],
            [8.0, 7.5, 7.0],
        ]
    )

    shaped = _as_layer_cell_heads(values, nlay=2, ncpl=3)

    assert shaped.shape == (2, 3)
    assert np.allclose(shaped, values)


def test_heads_shape_helper_flattens_structured_layer_row_col_arrays():
    values = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)

    shaped = _as_layer_cell_heads(values, nlay=2, ncpl=12)

    assert shaped.shape == (2, 12)
    assert np.allclose(shaped[0], values[0].reshape(-1))
    assert np.allclose(shaped[1], values[1].reshape(-1))


def test_cell_vector_helper_flattens_structured_top_or_bottom_arrays():
    values = np.arange(3 * 4, dtype=float).reshape(3, 4)

    shaped = _as_cell_vector(values, ncpl=12, label="top")

    assert shaped.shape == (12,)
    assert np.allclose(shaped, values.reshape(-1))


def test_contour_segments_clip_to_active_domain_and_support_cubic_method():
    cells = []
    values = []
    for row in range(3):
        for col in range(3):
            cells.append(box(col, row, col + 1, row + 1))
            values.append(row + 0.5)
    gdf = gpd.GeoDataFrame({"geometry": cells}, crs="EPSG:3857")
    vor = _DummyVor(gdf)
    active_domain = gdf.drop(index=4).union_all()
    full_domain = gdf.union_all()

    segments = contour_line_segments(
        vor,
        values,
        levels=[1.5],
        clip_geometry=active_domain,
    )
    smooth_segments = contour_line_segments(
        vor,
        values,
        levels=[1.5],
        clip_geometry=full_domain,
        method="cubic",
        resolution=40,
    )

    assert segments
    assert smooth_segments
    assert {segment["level"] for segment in smooth_segments} == {segment["level"] for segment in segments}
    assert max(len(segment["x"]) for segment in smooth_segments) > max(len(segment["x"]) for segment in segments)
    for segment in segments:
        line = LineString(zip(segment["x"], segment["y"], strict=False))
        assert active_domain.covers(line)
    for segment in smooth_segments:
        line = LineString(zip(segment["x"], segment["y"], strict=False))
        assert full_domain.covers(line)


def test_scalar_contour_levels_are_intervals():
    levels = _resolve_contour_levels(10, np.asarray([93.2, 117.8]))
    decimal_levels = _resolve_contour_levels(2.5, np.asarray([1.1, 8.2]))
    explicit = _resolve_contour_levels([101.0, 107.0], np.asarray([90.0, 120.0]))

    assert np.allclose(levels, [90.0, 100.0, 110.0, 120.0])
    assert np.allclose(decimal_levels, [0.0, 2.5, 5.0, 7.5, 10.0])
    assert np.allclose(explicit, [101.0, 107.0])


def test_points_to_latlon_reprojects_model_coordinates_onto_the_map_frame():
    """Track records and observation points arrive in model x/y; maps live in WGS84.

    The grid's own polygons are reprojected the same way, so a point placed by
    this helper lands inside the cell footprints the choropleth draws.
    """

    import pytest

    from myflopy.modflow.mf6.grid.geometry import points_to_latlon

    # An asymmetric footprint far from (0, 0), so a lon/lat swap or a collapse to
    # one point cannot pass: at the origin lon and lat are numerically alike.
    gdf = gpd.GeoDataFrame(
        {"geometry": [box(1_100_000.0, 500_000.0, 1_120_000.0, 560_000.0)]},
        crs="EPSG:2927",
    )
    vor = _DummyVor(gdf)
    vor.crs = "EPSG:2927"
    vor.crs_latlon = "EPSG:4326"

    xs = [1_105_000.0, 1_115_000.0]
    ys = [510_000.0, 550_000.0]
    lon, lat = points_to_latlon(vor, xs, ys)

    west, south, east, north = gdf.to_crs("EPSG:4326").total_bounds
    assert all(west <= value <= east for value in lon)
    assert all(south <= value <= north for value in lat)
    # distinct points stay distinct, and the transform round-trips to model coords
    assert lon[0] != lon[1] and lat[0] != lat[1]
    back = gpd.GeoSeries(gpd.points_from_xy(lon, lat), crs="EPSG:4326").to_crs("EPSG:2927")
    assert np.allclose(back.x.to_numpy(), xs, atol=0.01)
    assert np.allclose(back.y.to_numpy(), ys, atol=0.01)

    # a grid with no CRS cannot be mapped -- say so instead of raising from pyproj
    vor.crs = None
    with pytest.raises(ValueError, match="no CRS"):
        points_to_latlon(vor, [500.0], [500.0])
