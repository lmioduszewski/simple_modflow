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


#: The canvases plotly actually DREW, measured in a browser at 700x450 (ledger 311,
#: docs/reproducers/map_canvas_margins.py). Asserting against `map_canvas` alone
#: certifies the rule by itself -- that is how the colorbar's margin went unseen.
DRAWN_SINGLE = (619.0, 430.0)
DRAWN_ANIMATION = (540.0, 310.0)


def test_a_single_map_fits_plotlys_default_canvas(canonical_run, grid_bounds):
    figure = canonical_run.plot.map(layer=0).fig
    canvas = viz.map_canvas(figure.layout, colorbar=True)
    # 700x450 less the 20 px title margin and the colorbar's grown right margin.
    assert canvas == (700.0 - viz.COLORBAR_MARGIN, 430.0)
    assert _contains(figure.layout.map.zoom, canvas, grid_bounds)
    x, y = _visible_share(figure.layout.map.zoom, DRAWN_SINGLE, grid_bounds)
    assert x >= 1.0 and y >= 1.0


#: A WIDE window, so WIDTH is the tighter axis -- on the square canonical grid the
#: height always is, and a test of the grid alone cannot see a width error.
WIDE = (0.0, 900.0, 2100.0, 1200.0)


def test_a_wide_window_fits_the_canvases_plotly_really_drew(canonical_run):
    """Measured in a browser (docs/reproducers/map_canvas_margins.py): before the
    fix a standalone map at 700x450 showed 97.3% of this window's width and an
    animation 91.4% -- plotly grew the margins for the colorbar and controls after
    the fit. Literal measured canvases, not `map_canvas`'s own numbers."""

    # No padding: the window must fit EXACTLY, or the padding hides a margin error.
    single = canonical_run.plot.map(layer=0, extent=WIDE, bounds_padding=0.0)
    x, y = _visible_share(single.fig.layout.map.zoom, DRAWN_SINGLE, single.extent.wgs84)
    assert x >= 1.0 and y >= 1.0
    animation = canonical_run.hds.animate(extent=WIDE, bounds_padding=0.0)
    x, y = _visible_share(animation.layout.map.zoom, DRAWN_ANIMATION, single.extent.wgs84)
    assert x >= 1.0 and y >= 1.0


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
        canvas = viz.map_canvas(figure.layout, subplot, colorbar=True)
        assert _contains(figure.layout[subplot].zoom, canvas, grid_bounds)


def test_an_animation_fits_the_canvas_its_controls_leave(canonical_run, grid_bounds):
    """The slider and play buttons grow both side margins to 80 px at draw time;
    fitting the declared 50/0 cut a wide window to 91.4% (ledger 311)."""

    figure = canonical_run.hds.animate()
    canvas = viz.map_canvas(figure.layout, colorbar=True, controls=True)
    assert _contains(figure.layout.map.zoom, canvas, grid_bounds)
    x, y = _visible_share(figure.layout.map.zoom, DRAWN_ANIMATION, grid_bounds)
    assert x >= 1.0 and y >= 1.0


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
    # Kaleido's default 700x500, less the 20 px title margin and the colorbar's.
    assert _contains(written["zoom"], (700.0 - viz.COLORBAR_MARGIN, 480.0), grid_bounds)
    assert written["zoom"] != screen[0]
    assert (figure.layout.map.zoom, figure.layout.map.center.lat) == screen

    picture.save(tmp_path / "wide.png", width=1600, height=500)
    assert _contains(written["zoom"], (1600.0 - viz.COLORBAR_MARGIN, 480.0), grid_bounds)


def test_an_explicit_export_size_wins_over_the_figures_own(canonical_run, grid_bounds, monkeypatch, tmp_path):
    """plotly's `to_image(width=)` beats `layout.width`; the refit fitted the
    figure's 1400x900 and the 500x400 PNG showed 47.5% of the grid (ledger 311)."""

    picture = canonical_run.plot.map(layer=0)
    figure = picture.fig
    figure.update_layout(width=1400, height=900)
    written = {}
    monkeypatch.setattr(figure, "write_image", lambda path, **kw: written.update(zoom=figure.layout.map.zoom))
    # A TALL, narrow export: width is the tighter axis here and not at 1400x900,
    # so fitting the wrong one gives a different zoom.
    picture.save(tmp_path / "small.png", width=300, height=800)
    assert _contains(written["zoom"], (300.0 - viz.COLORBAR_MARGIN, 780.0), grid_bounds)


def test_an_export_restores_the_screen_centre_and_keeps_a_view_the_caller_set(
    canonical_run, monkeypatch, tmp_path
):
    """The saved 'screen' centre was the live object the refit mutated, so the
    export centre stuck; and a view set on `.fig` was overwritten (ledger 311)."""

    picture = canonical_run.plot.map(layer=0)
    figure = picture.fig
    monkeypatch.setattr(figure, "write_image", lambda path, **kw: None)
    before = figure.layout.map.center.to_plotly_json()
    picture.save(tmp_path / "a.png", width=1600, height=300)
    assert figure.layout.map.center.to_plotly_json() == before

    figure.update_layout(map_center={"lat": 45.0, "lon": -126.0}, map_zoom=9.0)
    assert picture.export_map_view((700, 500)) is None  # the caller's view stands
    picture.save(tmp_path / "b.png")
    assert (figure.layout.map.center.lat, figure.layout.map.zoom) == (45.0, 9.0)


def test_save_restores_a_centre_it_moved_for_the_export(tmp_path):
    """`layout.map.center` is a live object `update_layout` mutates in place; the
    restore held it, and wrote the EXPORT centre back (ledger 311). Pinned on a
    stub whose export centre differs, since a map's own never does."""

    class Stub(viz.Picture):
        def __init__(self):
            self._figure = viz.Fig()
            self._figure.update_layout(map={"center": {"lat": 1.0, "lon": 2.0}, "zoom": 5.0})

        @property
        def fig(self):
            return self._figure

        def export_map_view(self, size):
            return {"center": {"lat": 9.0, "lon": 9.0}, "zoom": 7.0}

    stub = Stub()
    seen = {}
    stub._figure.write_image = lambda path, **kw: seen.update(
        center=stub._figure.layout.map.center.to_plotly_json(), zoom=stub._figure.layout.map.zoom
    )
    stub.save(tmp_path / "x.png")
    assert seen == {"center": {"lat": 9.0, "lon": 9.0}, "zoom": 7.0}
    assert stub.fig.layout.map.center.to_plotly_json() == {"lat": 1.0, "lon": 2.0}
    assert stub.fig.layout.map.zoom == 5.0


def test_a_fixed_zoom_is_not_refitted_for_export(canonical_run):
    picture = canonical_run.plot.map(layer=0, fit_bounds=False, zoom=12)
    assert picture.export_map_view((700, 500)) is None


def test_a_hillshade_adds_a_layer_and_leaves_the_view_alone(monkeypatch):
    """It set zoom 12 centred on the raster after the fit, overriding
    `fit_bounds` and an explicit `zoom` alike."""

    hillshade = ChoroplethHillshadeBackground.__new__(ChoroplethHillshadeBackground)
    hillshade._bounds = (-127.0, 45.0, -126.0, 46.0)
    monkeypatch.setattr(hillshade, "data_uri", lambda: "data:image/png;base64,")
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
