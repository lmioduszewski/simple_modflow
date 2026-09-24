"""`map(extent=...)` -- the window a map opens on, at every scope and backend
(ledger 309).

The user's rules, each pinned here: an extent sets the VIEW only (no crop, no
colour re-range); one extent frames every panel of a mosaic and a LIST frames each
panel on its own; it is a box, points, a geometry (its bounds) or a name. Matplotlib
is where it earns its keep -- a static map cannot be zoomed -- so every static
renderer is held to the same padded box.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import shapely as shp

from myflopy import plot, viz
from myflopy.modflow.utils.datatypes.extent import is_extent_list, resolve_extent

BOX = (500.0, 600.0, 1000.0, 1100.0)  # inside the canonical 0..2100 ft grid
PADDED = (475.0, 575.0, 1025.0, 1125.0)  # bounds_padding=0.05 on each side


def _limits(axis):
    """``(xmin, ymin, xmax, ymax)`` -- shapely's order, like the extents."""

    (x0, x1), (y0, y1) = axis.get_xlim(), axis.get_ylim()
    return x0, y0, x1, y1


def _panel_axes(figure):
    return [axis for axis in figure.axes if axis.collections and axis.get_label() != "<colorbar>"]


def _mercator_y(latitude):
    radians = np.radians(latitude)
    return (1.0 - np.log(np.tan(radians) + 1.0 / np.cos(radians)) / np.pi) / 2.0


def _fits(view, canvas, wgs84, *, padding=0.05):
    """The view shows the whole box, fitted on its tighter axis (see test_map_framing)."""

    west, south, east, north = wgs84
    zoom = view.zoom
    x = 360.0 * canvas[0] / (512.0 * 2.0**zoom) / (east - west)
    y = canvas[1] / (512.0 * 2.0**zoom) / abs(_mercator_y(north) - _mercator_y(south))
    center_ok = view.center.lon == pytest.approx((west + east) / 2) and view.center.lat == pytest.approx(
        (south + north) / 2
    )
    return center_ok and min(x, y) == pytest.approx(1.0 + 2.0 * padding) and x >= 1 and y >= 1


# --- what counts as an extent ------------------------------------------------------


@pytest.fixture(scope="module")
def grid(canonical_run):
    vor = canonical_run.vor
    return {"crs": vor.crs, "bounds": vor.gdf_vorPolys.total_bounds}


def test_four_numbers_are_shapely_order(grid):
    assert resolve_extent(BOX, **grid).native == BOX


@pytest.mark.parametrize(
    "spelling",
    [
        [(500, 600), (1000, 1100), (700, 800)],        # a run of points -> their box
        np.array([[500, 600], [1000, 1100]]),
        shp.box(*BOX),                                  # a polygon -> its bounds
        shp.LineString([(500, 1100), (1000, 600)]),     # a transect -> its bounds
    ],
)
def test_points_and_geometries_frame_their_bounds(grid, spelling):
    assert resolve_extent(spelling, **grid).native == pytest.approx(BOX)


def test_a_frame_in_another_crs_is_reprojected(grid):
    """A GeoDataFrame carries its own CRS; bare numbers are the grid's."""

    wgs84 = gpd.GeoDataFrame(geometry=[shp.box(*BOX)], crs=grid["crs"]).to_crs(4326)
    assert resolve_extent(wgs84, **grid).native == pytest.approx(BOX, abs=1e-3)


def test_the_wgs84_box_contains_the_native_one(grid):
    extent = resolve_extent(BOX, **grid)
    corners = gpd.GeoSeries(shp.points([(BOX[0], BOX[1]), (BOX[2], BOX[3])]), crs=grid["crs"]).to_crs(4326)
    west, south, east, north = extent.wgs84
    assert west <= corners.x.min() and east >= corners.x.max()
    assert south <= corners.y.min() and north >= corners.y.max()


@pytest.mark.parametrize(
    ("value", "error", "match"),
    [
        ((1000, 600, 500, 1100), ValueError, "NOT FloPy's"),          # xmin > xmax
        ((5000, 5000, 6000, 6000), ValueError, "both read in"),       # off the grid
        ((-122.5, 47.0, -122.4, 47.1), ValueError, "off the grid"),   # lon/lat as feet
        (shp.Point(1, 1), ValueError, "has no area"),
        ([BOX, BOX], TypeError, "only a mosaic"),
        ("nowhere", KeyError, "names no area"),
        (object(), TypeError, "is not an extent"),
    ],
)
def test_what_cannot_frame_a_map_raises_by_name(grid, value, error, match):
    with pytest.raises(error, match=match):
        resolve_extent(value, **grid)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (BOX, False),                                   # four numbers: ONE box
        ([(0, 0), (1, 1)], False),                      # points: ONE box
        ([(0, 0), (1, 1), (2, 2), (3, 3)], False),      # four POINTS, still one
        ([BOX, BOX], True),                             # boxes: one per panel
        (["site", "north"], True),                      # names: one per panel
        ([BOX, None, shp.box(*BOX)], True),
        ("site", False),
        (None, False),
    ],
)
def test_one_extent_is_told_from_a_list_by_its_items(value, expected):
    assert is_extent_list(value) is expected


# --- a single map, both backends ---------------------------------------------------


def test_a_plotly_map_opens_on_the_extent(canonical_run):
    picture = canonical_run.plot.map(layer=0, extent=BOX)
    view = picture.fig.layout.map
    assert _fits(view, viz.map_canvas(picture.fig.layout), picture.extent.wgs84)
    # ...and the view is the only thing that moved.
    whole = canonical_run.plot.map(layer=0)
    assert np.array_equal(
        np.asarray(picture.linear_zs, dtype=float), np.asarray(whole.linear_zs, dtype=float),
        equal_nan=True,
    )
    assert len(picture.fig.data[0].z) == canonical_run.vor.ncpl


@pytest.mark.parametrize(
    "draw",
    [
        lambda m: m.plot.map(layer=0, extent=BOX, backend="mpl"),
        lambda m: plot.map(m, layer=0, extent=BOX, backend="mpl"),
        lambda m: m.vor.plot.map(extent=BOX, backend="mpl"),
        lambda m: m.packages.npf.k.map(extent=BOX, backend="mpl"),
        lambda m: m.hds.map(extent=BOX, backend="mpl"),
        lambda m: m.plot.map(type="trans_max", extent=BOX, backend="mpl"),
        lambda m: m.plot.map(layer=0, extent=BOX).plot_mpl(),
    ],
    ids=["bound", "free", "grid", "noun", "heads", "trans", "plot_mpl"],
)
def test_a_static_map_is_clipped_to_the_padded_extent(canonical_run, draw):
    assert _limits(draw(canonical_run).axes[0]) == pytest.approx(PADDED)


def test_a_static_map_without_an_extent_still_autoscales(canonical_run):
    """Only an extent sets limits; a default picture must not move."""

    axis = canonical_run.plot.map(layer=0, backend="mpl").axes[0]
    assert axis.get_autoscalex_on() and axis.get_autoscaley_on()


def test_features_far_off_the_grid_cannot_widen_a_framed_static_map(canonical_run):
    far = shp.Point(50_000.0, 50_000.0)
    axis = canonical_run.plot.map(layer=0, extent=BOX, features=far, backend="mpl").axes[0]
    assert _limits(axis) == pytest.approx(PADDED)


def test_bounds_padding_pads_the_extent(canonical_run):
    axis = canonical_run.plot.map(layer=0, extent=BOX, bounds_padding=0.0, backend="mpl").axes[0]
    assert _limits(axis) == pytest.approx(BOX)


def test_an_extent_with_no_fit_is_refused(canonical_run):
    with pytest.raises(ValueError, match="fit_bounds=False"):
        canonical_run.plot.map(layer=0, extent=BOX, fit_bounds=False)


def test_a_png_export_is_refitted_on_the_extent(canonical_run):
    picture = canonical_run.plot.map(layer=0, extent=BOX)
    view = picture.export_map_view((700, 500))
    center = view["center"]
    west, south, east, north = picture.extent.wgs84
    assert (center["lon"], center["lat"]) == pytest.approx(((west + east) / 2, (south + north) / 2))


# --- composers ---------------------------------------------------------------------

SUBPLOTS = ("map", "map2", "map3", "map4")


def test_one_extent_frames_every_mosaic_panel_and_keeps_them_linked(canonical_run):
    figure = canonical_run.hds.mosaic(by="layer", extent=BOX)
    views = {(round(figure.layout[s].zoom, 9), round(figure.layout[s].center.lon, 9)) for s in SUBPLOTS}
    assert len(views) == 1
    assert figure._post_scripts  # the live pan/zoom link


def test_a_list_of_extents_frames_each_panel_on_its_own(canonical_run):
    windows = [BOX, None, (0, 0, 2100, 2100), (1500, 1500, 2000, 2000)]
    figure = canonical_run.hds.mosaic(by="layer", extent=windows)
    zooms = [figure.layout[s].zoom for s in SUBPLOTS]
    assert zooms[0] == pytest.approx(zooms[3])      # two boxes of one size
    assert zooms[1] == pytest.approx(zooms[2])      # the grid, twice
    assert zooms[0] > zooms[1]
    assert not getattr(figure, "_post_scripts", None)  # a pan would snap them together


def test_a_static_mosaic_frames_each_panel_on_its_own(canonical_run):
    windows = [BOX, None, (0, 0, 2100, 2100), (1500, 1500, 2000, 2000)]
    figure = canonical_run.hds.mosaic(by="layer", extent=windows, backend="mpl")
    limits = [_limits(axis) for axis in _panel_axes(figure)]
    assert limits[0] == pytest.approx(PADDED)
    assert limits[3] == pytest.approx((1475.0, 1475.0, 2025.0, 2025.0))


def test_syncing_panels_on_different_windows_is_refused(canonical_run):
    with pytest.raises(ValueError, match="DIFFERENT extents"):
        canonical_run.hds.mosaic(by="layer", extent=[BOX, None, None, None], sync_views=True)


def test_a_list_of_the_wrong_length_names_the_panels(canonical_run):
    with pytest.raises(ValueError, match="2 windows for 4 panels"):
        canonical_run.hds.mosaic(by="layer", extent=[BOX, BOX])


def test_the_free_mosaic_takes_a_list_with_none_for_a_panels_own(canonical_run):
    own = (1500.0, 1500.0, 2000.0, 2000.0)
    figure = viz.mosaic(
        [canonical_run.plot.map(layer=0), canonical_run.plot.map(layer=1, extent=own)],
        extent=[BOX, None],
    )
    first = canonical_run.plot.map(layer=0, extent=BOX).extent.wgs84
    second = canonical_run.plot.map(layer=1, extent=own).extent.wgs84
    size = viz.map_canvas(figure.layout, "map")
    assert _fits(figure.layout.map, size, first)
    assert _fits(figure.layout.map2, size, second)


def test_a_free_mosaic_refuses_a_window_for_a_panel_that_is_not_a_map(canonical_run):
    series = viz.Fig()
    with pytest.raises(TypeError, match="not a map"):
        viz.mosaic([canonical_run.plot.map(layer=0), series], extent=[BOX, BOX])


def test_an_animation_opens_on_its_extent_on_both_backends(canonical_run):
    figure = canonical_run.hds.animate(extent=BOX)
    extent = canonical_run.plot.map(layer=0, extent=BOX).extent
    assert _fits(figure.layout.map, viz.map_canvas(figure.layout), extent.wgs84)
    static = canonical_run.hds.animate(extent=BOX, backend="mpl")
    assert _limits(static._fig.axes[0]) == pytest.approx(PADDED)
    with pytest.raises(TypeError, match="ONE view"):
        canonical_run.hds.animate(extent=[BOX, BOX])


def test_a_png_animation_frame_keeps_its_extent(canonical_run):
    """`animate(backend="png")` renders EVERY frame through `_mpl_figure_of`, which
    calls the picture's `plot_mpl` -- so a framed picture is a framed frame."""

    from myflopy.modflow.mf6.interactive_plotting import _mpl_figure_of

    frames = [canonical_run.plot.map(layer=lay, extent=BOX) for lay in (0, 1)]
    assert type(plot.animate(frames, backend="png")).__name__ == "SliderAnimation"
    for picture in frames:
        assert _limits(_mpl_figure_of(picture).axes[0]) == pytest.approx(PADDED)
