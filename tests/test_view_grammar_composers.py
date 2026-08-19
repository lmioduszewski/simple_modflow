"""Tests for the unified grammar's composers -- ``mosaic``/``animate`` with
``kind="map"|"plot"|"section"`` -- and the free-form ``viz.mosaic`` panel composer
(Phase 2 of the panel-verbs/composer redesign).
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import pytest
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure

matplotlib.use("Agg")

from myflopy import viz
from myflopy.modflow.mf6.package_plotting import FieldMappable, SpatialView


# --- fakes ----------------------------------------------------------------------
class _FakeChoro:
    """Minimal Choro stand-in: exposes get_choropleth() like the real one."""

    def __init__(self, z):
        self._z = list(z)

    def get_choropleth(self):
        return go.Choroplethmap(z=self._z, colorscale="Viridis")


class _FakeGeoChoro(_FakeChoro):
    """A map panel that also carries geometry, so mosaic can fit its view.

    ``map_view`` echoes whatever extent it is framed with (its own bounds, or a
    shared/union extent passed by the composer) as a bbox-centered view, which
    lets the tests assert exactly how each subplot was framed.
    """

    def __init__(self, z, *, bounds):
        super().__init__(z)
        self._bounds = tuple(float(b) for b in bounds)

    @property
    def latlon_bounds(self):
        return self._bounds

    def map_view(self, *, bounds=None):
        west, south, east, north = (
            self._bounds if bounds is None else tuple(float(b) for b in bounds)
        )
        return {
            "style": "carto-voyager",
            "center": {"lat": (south + north) / 2.0, "lon": (west + east) / 2.0},
            "zoom": float(max(east - west, north - south)),
        }


class _FakeSeriesHost(SpatialView):
    value_name = "q"

    def __init__(self, frame, models=None):
        self._frame = frame
        self._models = models

    def get(self):
        return self._frame

    def _spatial_models(self):
        return self._models


def _series_frame():
    return pd.DataFrame(
        {
            "model": ["a"] * 4 + ["b"] * 4,
            "per": [0, 0, 1, 1] * 2,
            "layer": [0] * 8,
            "cell": [1, 2, 1, 2] * 2,
            "lake": [0, 1, 0, 1] * 2,
            "q": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        }
    )


# --- viz.mosaic: the free-form composer -----------------------------------------
def test_viz_mosaic_composes_maps_and_figures():
    timeseries = viz.Fig()
    timeseries.add_scatter(x=[0, 1], y=[1.0, 2.0], name="stage")
    panels = [
        ("map a", _FakeChoro([1.0, 2.0, 3.0])),
        ("series", timeseries),
        _FakeChoro([4.0, 5.0, np.nan]),  # unlabeled -> "Panel 3"
    ]
    fig = viz.mosaic(panels, ncols=2, title="mixed")
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 3  # two map traces + one scatter
    titles = [annotation.text for annotation in fig.layout.annotations]
    assert titles[:3] == ["map a", "series", "Panel 3"]
    # maps share one coloraxis with global finite limits
    assert fig.layout.coloraxis.cmin == 1.0 and fig.layout.coloraxis.cmax == 5.0


def test_viz_mosaic_diff_centers_the_shared_scale():
    fig = viz.mosaic([_FakeChoro([-1.0, 4.0])], diff=True)
    assert fig.layout.coloraxis.cmin == -4.0
    assert fig.layout.coloraxis.cmax == 4.0
    assert fig.layout.coloraxis.cmid == 0.0


def test_viz_mosaic_rejects_bad_input():
    with pytest.raises(ValueError):
        viz.mosaic([])
    with pytest.raises(TypeError):
        viz.mosaic(["not a panel"])


def test_viz_mosaic_shares_start_view_and_wires_live_sync_by_default():
    # two maps over disjoint extents -> both framed to the union (0,0,4,4)
    panels = [
        _FakeGeoChoro([1.0, 2.0], bounds=(0.0, 0.0, 2.0, 2.0)),
        _FakeGeoChoro([3.0, 4.0], bounds=(2.0, 2.0, 4.0, 4.0)),
    ]
    fig = viz.mosaic(panels, ncols=2)
    assert fig.layout.map.center.lat == 2.0 and fig.layout.map.center.lon == 2.0
    assert fig.layout.map2.center.lat == 2.0 and fig.layout.map2.center.lon == 2.0
    assert fig.layout.map.zoom == fig.layout.map2.zoom  # one shared start view
    # live pan/zoom sync handler attached, referencing both map subplots
    scripts = getattr(fig, "_post_scripts", [])
    assert scripts and "plotly_relayout" in scripts[0]
    assert '"map"' in scripts[0] and '"map2"' in scripts[0]


def test_viz_mosaic_sync_off_still_shares_start_view_without_live_js():
    # disjoint extents -> both STILL start at the shared union view (site-aligned)
    panels = [
        _FakeGeoChoro([1.0, 2.0], bounds=(0.0, 0.0, 2.0, 2.0)),
        _FakeGeoChoro([3.0, 4.0], bounds=(10.0, 10.0, 12.0, 12.0)),
    ]
    fig = viz.mosaic(panels, ncols=2, sync_views=False)
    assert fig.layout.map.center.lon == fig.layout.map2.center.lon == 6.0
    assert not getattr(fig, "_post_scripts", [])  # no live linking


def test_viz_mosaic_without_geometry_leaves_maps_unfitted():
    # a plain Choro (no geometry) can't be fitted -> map subplot stays default
    fig = viz.mosaic([_FakeChoro([1.0, 2.0])])
    assert fig.layout.map.center.lat is None
    assert fig.layout.map.zoom is None
    assert not getattr(fig, "_post_scripts", [])


def test_synced_mosaic_emits_pan_zoom_js_into_html(tmp_path):
    panels = [
        _FakeGeoChoro([1.0, 2.0], bounds=(0.0, 0.0, 2.0, 2.0)),
        _FakeGeoChoro([3.0, 4.0], bounds=(2.0, 2.0, 4.0, 4.0)),
    ]
    out = tmp_path / "mosaic.html"
    viz.mosaic(panels, ncols=2).write_html(str(out))
    html = out.read_text(encoding="utf-8")
    # the sync handler survives the full render pipeline into the HTML output
    assert "plotly_relayout" in html and "Plotly.relayout" in html


def test_leaf_map_mosaic_threads_sync_views_through_to_viz():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    boxes = {
        "a": _FakeGeoChoro([1.0, 2.0], bounds=(0.0, 0.0, 2.0, 2.0)),
        "b": _FakeGeoChoro([3.0, 4.0], bounds=(10.0, 10.0, 12.0, 12.0)),
    }
    host._spatial_map = lambda *, per, layer, model, **kwargs: boxes[model]

    synced = host.mosaic(by="model")  # sync_views defaults to True
    assert synced.layout.map.center.lon == synced.layout.map2.center.lon == 6.0
    assert getattr(synced, "_post_scripts", [])  # live sync wired

    independent = host.mosaic(by="model", sync_views=False)
    assert independent.layout.map.center.lon == 6.0  # still shares the start view
    assert not getattr(independent, "_post_scripts", [])  # but not linked live


def test_map_animation_frames_are_fitted_to_data():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    boxes = {
        "a": _FakeGeoChoro([1.0, 2.0], bounds=(0.0, 0.0, 2.0, 2.0)),
        "b": _FakeGeoChoro([3.0, 4.0], bounds=(4.0, 4.0, 6.0, 6.0)),
    }
    host._spatial_map = lambda *, per, layer, model, **kwargs: boxes[model]
    fig = host.animate(kind="map", over="model")
    assert isinstance(fig, viz.Fig)  # house figure, not a bare go.Figure
    assert len(fig.frames) == 2
    # the single animated map is framed to data (union 0,0,6,6), not the world
    assert fig.layout.map.center is not None and fig.layout.map.zoom is not None
    assert fig.layout.map.center.lon == 3.0


def test_choro_map_view_fits_and_shares_bounds():
    """The real Choro.map_view frames to data, and honors a shared extent."""

    from myflopy.modflow.utils.datatypes.choros import Choro

    class _FakeGDF:
        total_bounds = np.array([-122.0, 47.0, -121.0, 48.0])

    class _FakeVor:
        gdf_latlon = _FakeGDF()

    choro = Choro.__new__(Choro)  # bypass the heavy __init__
    choro.vor = _FakeVor()
    choro.fit_bounds = True
    choro.zoom = 13
    choro.bounds_padding = 0.05

    assert choro.latlon_bounds == (-122.0, 47.0, -121.0, 48.0)

    own = choro.map_view()
    assert own["center"] == {"lat": 47.5, "lon": -121.5}
    assert 0.0 < own["zoom"] <= 20.0

    wider = choro.map_view(bounds=(-124.0, 46.0, -120.0, 49.0))
    assert wider["center"] == {"lat": 47.5, "lon": -122.0}
    assert wider["zoom"] < own["zoom"]  # a bigger extent zooms out further

    choro.fit_bounds = False
    assert choro.map_view()["zoom"] == 13  # fixed zoom when not fitting

    choro.vor = None
    assert choro.latlon_bounds is None
    assert choro.map_view() == {"style": "carto-voyager"}


def test_leaf_plotly_map_mosaic_is_sugar_over_viz_mosaic():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    host._spatial_map = lambda **kwargs: _FakeChoro([1.0, float(2 + len(str(kwargs)))])
    fig = host.mosaic(by="model")
    assert isinstance(fig, go.Figure) and len(fig.data) == 2
    assert fig.layout.coloraxis.cauto is False  # shared scale applied


# --- mosaic(kind="plot") ---------------------------------------------------------
def test_plot_mosaic_facets_by_model_default_on_groups():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    fig = host.mosaic(kind="plot")
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 4  # 2 panels x 2 lake lines
    titles = [annotation.text for annotation in fig.layout.annotations]
    assert titles[:2] == ["a", "b"]
    # legend deduplicated: each lake shown once
    assert sum(1 for trace in fig.data if trace.showlegend) == 2


def test_plot_mosaic_facets_by_entity_and_layer():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    by_lake = host.mosaic(kind="plot", by="lake")
    titles = [annotation.text for annotation in by_lake.layout.annotations]
    assert titles[:2] == ["Lake 0", "Lake 1"]

    single = _FakeSeriesHost(_series_frame().drop(columns=["model"]))
    by_default = single.mosaic(kind="plot")  # single model -> entity facet
    titles = [annotation.text for annotation in by_default.layout.annotations]
    assert titles[:2] == ["Lake 0", "Lake 1"]


def test_plot_mosaic_mpl_backend_and_errors():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    fig = host.mosaic(kind="plot", backend="mpl")
    assert isinstance(fig, Figure)
    with pytest.raises(ValueError):
        host.mosaic(kind="plot", by="period")
    with pytest.raises(ValueError):
        host.mosaic(kind="plot", by="reach")  # column not present
    with pytest.raises(ValueError):
        host.mosaic(kind="hologram")


# --- animate(kind="plot") --------------------------------------------------------
def test_plot_animation_flips_over_models():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    fig = host.animate(kind="plot")
    # a house viz.Fig: pan-drag + scroll-zoom defaults, not a bare go.Figure
    assert isinstance(fig, viz.Fig)
    assert fig.layout.template.layout.dragmode == "pan"
    assert fig._config.get("scrollZoom") is True
    assert [frame.name for frame in fig.frames] == ["a", "b"]
    assert len(fig.frames[0].data) == 2  # lake lines per frame

    anim = host.animate(kind="plot", backend="mpl")
    assert isinstance(anim, FuncAnimation)

    with pytest.raises(ValueError):
        host.animate(kind="plot", over="period")
    with pytest.raises(ValueError):
        host.animate(kind="hologram")


def test_plot_animation_requires_model_axis():
    single = _FakeSeriesHost(_series_frame().drop(columns=["model"]))
    with pytest.raises(ValueError):
        single.animate(kind="plot")


# --- kind="section" plumbing -----------------------------------------------------------
def test_xs_kind_raises_cleanly_off_heads_leaves():
    host = _FakeSeriesHost(_series_frame(), models=["a", "b"])
    with pytest.raises(NotImplementedError):
        host.section()
    with pytest.raises(NotImplementedError):
        host.mosaic(kind="section")
    with pytest.raises(NotImplementedError):
        host.animate(kind="section")


def test_field_mappable_xs_dispatch():
    sentinel = object()

    class _Leaf:
        def section(self, *args, **kwargs):
            return sentinel

    class _Namespace(FieldMappable):
        _default_field = "q"

        def _field_names(self):
            return ["q"]

        @property
        def q(self):
            return _Leaf()

    assert _Namespace().section(line=None) is sentinel


# --- real renders on the canonical model ------------------------------------------
@pytest.mark.slow
def test_composers_render_on_canonical(canonical_run):
    """kind= composers render real figures/animations on all three surfaces."""

    from shapely.geometry import LineString

    import myflopy as mf
    from myflopy.project.model_group import ModelGroup

    model = canonical_run
    bounds = model.vor.gdf_vorPolys.total_bounds
    line = LineString(
        [
            (bounds[0] + 0.25 * (bounds[2] - bounds[0]), bounds[1] + 0.25 * (bounds[3] - bounds[1])),
            (bounds[0] + 0.75 * (bounds[2] - bounds[0]), bounds[1] + 0.75 * (bounds[3] - bounds[1])),
        ]
    )

    # -- single model: animated + mosaicked cross sections -------------------
    xs_animation = model.hds.animate(kind="section", line=line)  # over periods
    # the house viz.Fig carries the interaction defaults (pan, scroll-zoom)
    assert isinstance(xs_animation, viz.Fig) and len(xs_animation.frames) >= 2
    assert xs_animation.layout.template.layout.dragmode == "pan"
    assert isinstance(
        model.hds.animate(kind="section", line=line, backend="mpl"), FuncAnimation
    )
    xs_mosaic = model.hds.mosaic(kind="section", line=line, per=[0, 1])
    assert isinstance(xs_mosaic, go.Figure)
    assert isinstance(
        model.hds.mosaic(kind="section", line=line, per=[0, 1], backend="mpl"), Figure
    )

    # -- single model: entity-faceted series mosaic --------------------------
    stage_mosaic = model.packages.lak.results.mosaic(kind="plot", field="stage")
    assert isinstance(stage_mosaic, go.Figure)  # one panel per lake

    # -- group surface --------------------------------------------------------
    reloaded = mf.load_mf6_run(model.workspace)
    group = ModelGroup({"a": model, "b": reloaded}, reference="a")

    group_xs_mosaic = group.hds.mosaic(kind="section", line=line)  # panel per member
    titles = [annotation.text for annotation in group_xs_mosaic.layout.annotations]
    assert titles[:2] == ["a", "b"]

    # group map mosaic: real geometry -> subplots zoom to data, synced by default
    group_map_mosaic = group.hds.mosaic(by="model")
    assert group_map_mosaic.layout.map.zoom and group_map_mosaic.layout.map.zoom > 0
    assert group_map_mosaic.layout.map.center.lat == group_map_mosaic.layout.map2.center.lat
    assert group_map_mosaic.layout.map.zoom == group_map_mosaic.layout.map2.zoom
    # opting out fits each panel independently (still zoomed to data, not the world)
    independent = group.hds.mosaic(by="model", sync_views=False)
    assert independent.layout.map.zoom and independent.layout.map.zoom > 0

    group_plot_anim = group.packages.ghb.results.animate(kind="plot")
    assert [frame.name for frame in group_plot_anim.frames] == ["a", "b"]

    # -- diff surface ----------------------------------------------------------
    diff = group.diff()
    diff_xs_mosaic = diff.hds.mosaic(kind="section", line=line)
    assert isinstance(diff_xs_mosaic, go.Figure)

    # -- free-form composer: a map and a timeseries in one grid ----------------
    mixed = viz.mosaic(
        [
            ("heads", model.hds.map()),
            ("lake stage", model.packages.lak.results.stage.plot()),
        ],
        ncols=2,
    )
    assert isinstance(mixed, go.Figure) and len(mixed.data) >= 2


def test_viz_mosaic_carries_map_overlays_into_each_panel():
    """A panel's overlays (contours, location markers, pathlines) ride into the mosaic.

    Copying only the cell trace used to drop every overlay silently, so a mosaic
    of contoured maps -- or of pathline maps -- showed half the figure.
    """

    class _OverlaidChoro(_FakeChoro):
        def overlay_traces(self):
            return [go.Scattermap(mode="lines", lon=[0.0, 1.0], lat=[0.0, 1.0], name="path")]

    fig = viz.mosaic([_OverlaidChoro([1.0, 2.0]), _OverlaidChoro([3.0, 4.0])], ncols=2)

    assert [trace.type for trace in fig.data] == [
        "choroplethmap",
        "scattermap",
        "choroplethmap",
        "scattermap",
    ]
    # only the cell traces join the shared color scale; the overlays keep their own
    assert fig.layout.coloraxis.cmin == 1.0 and fig.layout.coloraxis.cmax == 4.0
    assert [trace.subplot for trace in fig.data] == ["map", "map", "map2", "map2"]


def test_real_choro_overlays_survive_composition(canonical_run):
    """The bug ledger 62 records: a mosaic of contoured maps lost its contours.

    Exercised on a real ``Choro`` rather than a fake, so it pins what
    ``overlay_traces()`` actually collects -- contour polylines and location
    markers, not just whatever a test double chooses to return.
    """

    import geopandas as gpd
    from shapely.geometry import LineString, Point

    vor = canonical_run.vor
    centroid = vor.gdf_vorPolys.geometry.iloc[0].centroid
    locs = gpd.GeoDataFrame(
        {"ExploName": ["well A", "a transect"]},
        geometry=[
            Point(centroid.x, centroid.y),
            # a geometry type the marker overlay does not draw: skipped, not fatal
            LineString([(centroid.x, centroid.y), (centroid.x + 10.0, centroid.y)]),
        ],
        crs=vor.crs,
    )
    choro = canonical_run.plot.map(layer=0, contours=True, contour_levels=3, locs=locs)

    overlays = choro.overlay_traces()
    names = [trace.name for trace in overlays]
    assert "well A" in names, names
    assert "a transect" not in names  # LineString skipped
    contours = [trace for trace in overlays if trace.mode == "lines"]
    assert contours, names
    assert all(trace.type == "scattermap" for trace in overlays)

    # ...and every one of them reaches the composed panel
    fig = viz.mosaic([("heads", choro)])
    assert len(fig.data) == 1 + len(overlays)
    assert sum(trace.type == "choroplethmap" for trace in fig.data) == 1
