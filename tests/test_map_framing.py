"""A plotly map opens showing the grid it draws (ledger 307).

The fit assumed a 1000x700 canvas that nothing draws. Measured in a browser on the
canonical model: a single map at 1084x584 showed 88.6% of the grid's height, every
panel of a 4-panel mosaic ~10% of the grid (each subplot got the single-map zoom),
and `.save("x.png")` at kaleido's 700x500 78.6% of its height. A hillshade then
reset every view to zoom 12 centred on the RASTER.

The check here is the one the critic validated against MapLibre to three digits:
at zoom z a canvas `w` pixels wide spans `360 * w / (512 * 2**z)` degrees of
longitude, and `h` pixels span `h / (512 * 2**z)` of the Web Mercator height.
"""

from __future__ import annotations

import numpy as np
import pytest

from myflopy import viz
from myflopy.modflow.utils.datatypes.choros import ChoroplethHillshadeBackground


def _mercator_y(latitude):
    radians = np.radians(latitude)
    return (1.0 - np.log(np.tan(radians) + 1.0 / np.cos(radians)) / np.pi) / 2.0


def _visible_share(zoom, canvas, bounds):
    """``(x, y)``: how many times the box fits across a ``canvas`` at ``zoom``."""

    west, south, east, north = bounds
    width, height = canvas
    lon_span = 360.0 * width / (512.0 * 2.0**zoom)
    y_span = height / (512.0 * 2.0**zoom)
    return lon_span / (east - west), y_span / abs(_mercator_y(north) - _mercator_y(south))


@pytest.fixture(scope="module")
def grid_bounds(canonical_run):
    return tuple(float(b) for b in canonical_run.vor.gdf_latlon.total_bounds)


def _contains(zoom, canvas, bounds, *, padding=0.05):
    x, y = _visible_share(zoom, canvas, bounds)
    # Fitted exactly on the tighter axis, with the padding on each side.
    return min(x, y) == pytest.approx(1.0 + 2.0 * padding) and x >= 1.0 and y >= 1.0


def test_a_single_map_fits_plotlys_default_canvas(canonical_run, grid_bounds):
    figure = canonical_run.plot.map(layer=0).fig
    canvas = viz.map_canvas(figure.layout)
    assert canvas == (700.0, 430.0)  # 700x450 less the map's 20 px title margin
    assert _contains(figure.layout.map.zoom, canvas, grid_bounds)


def test_a_bigger_real_canvas_only_adds_margin(canonical_run, grid_bounds):
    """The browser measurement that cropped: 995x564 of map at 1084x584."""

    zoom = canonical_run.plot.map(layer=0).fig.layout.map.zoom
    x, y = _visible_share(zoom, (995.0, 564.0), grid_bounds)
    assert x > 1.0 and y > 1.0


def test_a_single_map_is_centred_on_its_box_midpoint(canonical_run, grid_bounds):
    """The standalone map centred on the area centroid, the subplots on the
    midpoint -- two framings of one irregular domain."""

    west, south, east, north = grid_bounds
    center = canonical_run.plot.map(layer=0).fig.layout.map.center
    assert (center.lon, center.lat) == pytest.approx(((west + east) / 2, (south + north) / 2))


def test_every_mosaic_panel_fits_its_own_subplot(canonical_run, grid_bounds):
    figure = canonical_run.hds.mosaic(by="layer")  # 4 panels, 3 columns, 2 rows
    assert figure.layout.height == 2 * viz.MOSAIC_ROW_HEIGHT + 100 + 80
    for subplot in ("map", "map2", "map3", "map4"):
        canvas = viz.map_canvas(figure.layout, subplot)
        assert _contains(figure.layout[subplot].zoom, canvas, grid_bounds)


def test_an_animation_fits_its_own_canvas(canonical_run, grid_bounds):
    figure = canonical_run.hds.animate()
    assert _contains(figure.layout.map.zoom, viz.map_canvas(figure.layout), grid_bounds)


def test_a_png_is_refitted_for_its_export_size_and_the_screen_view_restored(
    canonical_run, grid_bounds, monkeypatch, tmp_path
):
    picture = canonical_run.plot.map(layer=0)
    figure = picture.fig
    screen = (figure.layout.map.zoom, figure.layout.map.center.lat)
    written = {}

    def fake_write_image(path, **kwargs):
        written["zoom"] = figure.layout.map.zoom
        written["size"] = (kwargs.get("width"), kwargs.get("height"))

    monkeypatch.setattr(figure, "write_image", fake_write_image)
    picture.save(tmp_path / "heads.png")
    # Kaleido's default 700x500, less the map's own 20 px margin.
    assert _contains(written["zoom"], (700.0, 480.0), grid_bounds)
    assert written["zoom"] != screen[0]
    assert (figure.layout.map.zoom, figure.layout.map.center.lat) == screen

    picture.save(tmp_path / "wide.png", width=1600, height=500)
    assert _contains(written["zoom"], (1600.0, 480.0), grid_bounds)


def test_a_fixed_zoom_is_not_refitted_for_export(canonical_run):
    picture = canonical_run.plot.map(layer=0, fit_bounds=False, zoom=12)
    assert picture.export_map_view((700, 500)) is None


def test_a_hillshade_adds_a_layer_and_leaves_the_view_alone(monkeypatch):
    """It set zoom 12 centred on the raster after the fit, overriding
    `fit_bounds` and an explicit `zoom` alike."""

    hillshade = ChoroplethHillshadeBackground.__new__(ChoroplethHillshadeBackground)
    hillshade._bounds = (-127.0, 45.0, -126.0, 46.0)
    monkeypatch.setattr(hillshade, "file_to_data_uri", lambda: "data:image/png;base64,")
    figure = viz.Fig()
    figure.update_layout(map={"center": {"lat": 1.0, "lon": 2.0}, "zoom": 17.0})

    hillshade.update_fig_layout(figure)
    assert figure.layout.map.zoom == 17.0
    assert (figure.layout.map.center.lat, figure.layout.map.center.lon) == (1.0, 2.0)
    assert figure.layout.map.layers[0].name == "hillshade"


def test_map_canvas_reads_the_figure_its_template_and_the_domain():
    figure = viz.Fig()
    figure.update_layout(width=1000, height=600, margin={"l": 10, "r": 10, "t": 30, "b": 20})
    assert viz.map_canvas(figure.layout) == (980.0, 550.0)
    # Unset margins fall back to the TEMPLATE's (the house template names them).
    bare = viz.Fig()
    template = bare.layout.template.layout.margin
    assert viz.map_canvas(bare.layout) == (
        700.0 - template.l - template.r,
        450.0 - template.t - template.b,
    )
    # A subplot is its domain's share of the plot area.
    grid = viz.subplots(1, 2, specs=[[{"type": "map"}, {"type": "map"}]])
    width, height = viz.map_canvas(grid.layout, "map2")
    domain = grid.layout.map2.domain
    assert width == pytest.approx((700.0 - 160.0) * (domain.x[1] - domain.x[0]))
    assert height == pytest.approx(450.0 - 180.0)
