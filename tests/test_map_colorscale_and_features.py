"""Log color scales that read in real units, every Plotly colorscale name, and
custom geometry over a map (``features=``) -- at every scope (2026-09-13).

Three things were measured before this work:

* ``logscale=True`` existed on every map but labelled the colorbar in log10
  units (``-3 … 2`` for a K field running 0.001 to 100) while the hover showed
  real values (ledger 74), took ``zmin``/``zmax`` in log10 units, and drew the
  static backend as a linear bar of exponents.
* ``colorscale=`` accepted eighteen names out of Plotly's ninety-odd, warning and
  drawing ``earth`` for the rest, and the static backend mapped names through a
  hand table that rendered ``RdBu`` mirrored (ledger 69/70).
* ``locs=`` drew a file as unlabelled black markers with no legend, no hover, no
  styling, skipped lines, and was ignored by the static backend.
"""

from __future__ import annotations

import inspect
import os
import warnings

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, MultiPoint, MultiPolygon, Point, Polygon

from myflopy import plot, viz
from myflopy.modflow.mf6.grid.plotting import GridPlots, _choropleth_factory
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.utils.datatypes.choros import mpl_colormap_for, plotly_colorscale_stops
from myflopy.modflow.utils.datatypes.features import (
    Feature,
    feature_traces,
    mpl_color,
    resolve_features,
)
from myflopy.plot import NOUN_MAP_PARAMS, ModelPlots


def _two_cell_vor():
    """The smallest grid a Choro will build on (mirrors test_colorscale_policy.py)."""

    verts = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        dtype=float,
    )
    return VoronoiGridPlus(
        verts=verts,
        iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float),
    )


@pytest.fixture
def vor():
    return _two_cell_vor()


def _wells(crs):
    return gpd.GeoDataFrame(
        {"name": ["W1", "W2"], "q": [-100.0, -250.5]},
        geometry=[Point(0.5, 0.5), Point(1.5, 0.5)],
        crs=crs,
    )


POND = Polygon([(0.2, 0.2), (0.8, 0.2), (0.8, 0.8), (0.2, 0.8)])
TRANSECT = LineString([(0.0, 0.5), (2.0, 0.5)])


def _log_map(vor, **kwargs):
    return _choropleth_factory(vor, custom_zs=[0.001, 100.0], logscale=True, **kwargs)


# ---------------------------------------------------------------------------
# logscale
# ---------------------------------------------------------------------------


def test_a_log_map_carries_log10_values_but_labels_its_colorbar_in_real_units(vor):
    """Plotly has no log color axis, so the trace carries log10 -- and the bar
    must then say what each decade IS, or it reads ``-3 … 2``."""

    trace = _log_map(vor).get_choropleth()
    assert list(trace.z) == pytest.approx([-3.0, 2.0])
    assert trace.colorbar.tickvals == (-3.0, -2.0, -1.0, 0.0, 1.0, 2.0)
    assert trace.colorbar.ticktext == ("0.001", "0.01", "0.1", "1", "10", "100")


def test_log_limits_are_given_in_data_units(vor):
    """``zmin=1e-4, zmax=1e3`` means K from 0.0001 to 1000 -- the same units the
    bar and the hover use -- and reaches the trace as log10."""

    trace = _log_map(vor, zmin=1e-4, zmax=1e3).get_choropleth()
    assert (trace.zmin, trace.zmax) == (-4.0, 3.0)
    assert trace.colorbar.ticktext[0] == "0.0001"
    assert trace.colorbar.ticktext[-1] == "1000"  # not "1e+03"


def test_a_log_map_transforms_zmid_too(vor):
    trace = _log_map(vor, zmid=1.0).get_choropleth()
    assert trace.zmid == 0.0


@pytest.mark.parametrize("kwargs", [{"zmin": 0}, {"zmax": -1.0}, {"zmid": 0.0}])
def test_a_non_positive_log_limit_is_refused_at_the_call(vor, kwargs):
    with pytest.raises(ValueError, match="data units and must be finite and positive"):
        _log_map(vor, **kwargs)


def test_a_caller_colorbar_key_wins_over_the_decade_labels(vor):
    """A title rides along; the caller's own ticks replace the default ones."""

    trace = _log_map(vor, colorbar={"title": "K (m/d)"}).get_choropleth()
    assert trace.colorbar.title.text == "K (m/d)"
    assert trace.colorbar.ticktext[0] == "0.001"

    own = _log_map(vor, colorbar={"tickvals": [-3, 2], "ticktext": ["lo", "hi"]}).get_choropleth()
    assert own.colorbar.ticktext == ("lo", "hi")


def test_a_log_map_with_no_finite_value_still_builds(vor):
    trace = _choropleth_factory(
        vor, custom_zs=[float("nan"), -1.0], logscale=True
    ).get_choropleth()
    assert all(not np.isfinite(z) for z in trace.z)
    assert trace.colorbar.tickvals is None


def test_the_static_backend_draws_a_log_map_through_a_log_norm(vor):
    """A real log axis in data units, not a linear bar of exponents."""

    figure = _log_map(vor).plot_mpl()
    bar = figure.axes[-1]
    assert bar.get_yscale() == "log"
    lo, hi = bar.get_ylim()
    assert (lo, hi) == pytest.approx((0.001, 100.0))

    limited = _log_map(vor, zmin=1e-4, zmax=1e3).plot_mpl()
    assert limited.axes[-1].get_ylim() == pytest.approx((1e-4, 1e3))


def test_the_static_backend_applies_caller_ticks_in_data_units(vor):
    """Caller ticks are the trace's log10 positions; on the data-unit log axis
    they land at ``10 ** tick`` and carry the caller's text (the PEST change
    map's ``0.1x / 1x / 10x``)."""

    figure = _log_map(
        vor, colorbar={"tickvals": [-2, 0, 2], "ticktext": ["0.01x", "1x", "100x"]}
    ).plot_mpl()
    bar = figure.axes[-1]
    assert list(bar.get_yticks()) == pytest.approx([0.01, 1.0, 100.0])
    assert [t.get_text() for t in bar.get_yticklabels()] == ["0.01x", "1x", "100x"]


def test_a_log_map_contours_the_data_units_not_the_exponents(canonical_run):
    """``contours=True`` contours the mapped field; on a log map that field is
    the real one, so a level reads ``0.01`` and never ``-2``."""

    ncpl = canonical_run.vor.ncpl
    values = np.logspace(-3, 2, ncpl)
    choro = canonical_run.plot.map(values=values, logscale=True, contours=True, contour_levels=4)
    choro.overlay_traces()
    levels = [segment["level"] for segment in choro._contour_segments]
    assert levels, "no contours were built"
    assert min(levels) > 0.0


def test_a_mosaic_of_log_maps_relabels_its_pooled_colorbar_by_itself(vor):
    figure = viz.mosaic(
        [("a", _log_map(vor)), ("b", _choropleth_factory(vor, custom_zs=[0.01, 10.0], logscale=True))],
        ncols=2,
    )
    coloraxis = figure.layout.coloraxis
    assert (coloraxis.cmin, coloraxis.cmax) == pytest.approx((-3.0, 2.0))
    assert coloraxis.colorbar.ticktext == ("0.001", "0.01", "0.1", "1", "10", "100")


def test_a_mosaic_refuses_to_pool_log_and_linear_maps(vor):
    with pytest.raises(ValueError, match="same logscale"):
        viz.mosaic([("a", _log_map(vor)), ("b", _choropleth_factory(vor, custom_zs=[1.0, 2.0]))])


def test_an_explicit_mosaic_colorbar_still_wins(vor):
    figure = viz.mosaic(
        [("a", _log_map(vor)), ("b", _log_map(vor))], ncols=2,
        colorbar={"tickvals": [-3, 2], "ticktext": ["lo", "hi"]},
    )
    assert figure.layout.coloraxis.colorbar.ticktext == ("lo", "hi")


def test_an_animation_of_a_log_map_labels_the_shared_axis_in_real_units(canonical_run):
    """``.ani`` pools every frame onto one color axis; its bar is rebuilt over
    the SHARED range, not copied from the first frame."""

    animation = canonical_run.plot.map(layer=0, logscale=True).ani
    coloraxis = animation.layout.coloraxis
    assert np.isfinite([coloraxis.cmin, coloraxis.cmax]).all()
    labels = coloraxis.colorbar.ticktext
    assert labels and all(float(text) > 0 for text in labels)


def test_the_hover_keeps_the_real_value_on_a_log_map(vor):
    """What the reader hovers is the field, never its logarithm."""

    trace = _choropleth_factory(
        vor, custom_zs=[0.001, 100.0], custom_hover={"K": [0.001, 100.0]}, logscale=True
    ).get_choropleth()
    assert any("0.001" in str(row) for row in trace.customdata)


def test_a_log_map_reads_the_same_at_the_free_bound_and_noun_scopes(canonical_run):
    for picture in (
        plot.map(canonical_run, layer=0, logscale=True),
        canonical_run.plot.map(layer=0, logscale=True),
        canonical_run.hds.map(logscale=True),
    ):
        colorbar = picture.get_choropleth().colorbar
        assert colorbar.ticktext, "log map without real-unit labels"
        assert all(float(text) > 0 for text in colorbar.ticktext)


# ---------------------------------------------------------------------------
# colorscale names
# ---------------------------------------------------------------------------


def test_every_plotly_named_colorscale_is_accepted(vor):
    import plotly.colors as pc

    for name in pc.named_colorscales():
        assert _choropleth_factory(vor, custom_zs=[1.0, 2.0], colorscale=name).colorscale == name
    for spelled in ("Viridis", "VIRIDIS", "viridis_r", "Turbo", "RdBu_r"):
        choro = _choropleth_factory(vor, custom_zs=[1.0, 2.0], colorscale=spelled)
        assert choro.colorscale == spelled
        # and Plotly itself takes it, which is what the picture depends on
        assert choro.get_choropleth().colorscale


def test_an_unknown_colorscale_name_warns_and_falls_back_to_earth(vor, caplog):
    with caplog.at_level("WARNING"):
        choro = _choropleth_factory(vor, custom_zs=[1.0, 2.0], colorscale="nonsense")
    assert choro.colorscale == "earth"
    assert "nonsense" in caplog.text and "viridis" in caplog.text


def test_explicit_stops_still_pass_straight_through(vor):
    stops = [[0.0, "#d62728"], [1.0, "#1f77b4"]]
    assert _choropleth_factory(vor, custom_zs=[1.0, 2.0], colorscale=stops).colorscale == stops


def test_the_static_colormap_is_built_from_the_plotly_stops():
    """A name renders the SAME way on both backends -- which the hand table did
    not manage: it sent ``RdBu`` to matplotlib's reversed map (ledger 69/70)."""

    stops = plotly_colorscale_stops("RdBu")
    assert stops is not None and stops[0][0] == 0.0
    colormap = mpl_colormap_for("RdBu")
    low = np.asarray(colormap(0.0)[:3])
    assert low == pytest.approx(np.asarray(mpl_color(stops[0][1])[:3]), abs=1e-6)
    assert low[0] > low[2]  # red end low, on matplotlib too
    reversed_map = mpl_colormap_for("RdBu_r")
    assert reversed_map(0.0)[2] > reversed_map(0.0)[0]  # blue end low once reversed

    turbo = mpl_colormap_for("Turbo")
    assert turbo(0.0)[:3] != mpl_colormap_for("earth")(0.0)[:3]  # not gist_earth any more
    assert plotly_colorscale_stops("nonsense") is None
    assert mpl_colormap_for("nonsense").name == "plotly_earth"


def test_mpl_color_accepts_plotly_rgb_strings():
    assert mpl_color("rgb(214,39,40)") == pytest.approx((214 / 255, 39 / 255, 40 / 255, 1.0))
    assert mpl_color("rgba(0, 0, 255, 0.5)") == pytest.approx((0.0, 0.0, 1.0, 0.5))
    assert mpl_color("#0072B2")[:3] == pytest.approx((0.0, 0x72 / 255, 0xB2 / 255))


# ---------------------------------------------------------------------------
# features
# ---------------------------------------------------------------------------


def _by_name(traces, name):
    return [trace for trace in traces if trace.name == name]


def test_features_take_a_geodataframe_a_geometry_and_a_file(vor, tmp_path):
    path = tmp_path / "sites.gpkg"
    _wells(vor.crs).to_file(path, layer="wells", driver="GPKG")
    gpd.GeoDataFrame({"id": [1]}, geometry=[POND], crs=vor.crs).to_file(path, layer="ponds", driver="GPKG")

    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features=[_wells(vor.crs), POND, Feature(path, layer="ponds"), str(path)],
    )
    names = [trace.name for trace in choro.overlay_traces()]
    assert names == ["features", "features 2", "sites:ponds", "sites"]
    # the bare path read the FIRST layer (the wells), as GDAL does
    assert _by_name(choro.overlay_traces(), "sites")[0].mode == "markers"


def test_a_dict_names_the_layers_and_a_feature_carries_its_own_name(vor):
    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features={"Wells": _wells(vor.crs), "Pond": Feature(POND, name="ignored"), "Transect": TRANSECT},
    )
    assert [trace.name for trace in choro.overlay_traces()] == ["Wells", "Pond", "Transect"]
    named = _choropleth_factory(vor, custom_zs=[1.0, 2.0], features=Feature(TRANSECT, name="A-A'"))
    assert [trace.name for trace in named.overlay_traces()] == ["A-A'"]


def test_each_layer_is_one_legend_entry_with_a_house_color(vor):
    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features={"Wells": _wells(vor.crs), "Sites": gpd.GeoDataFrame(
            geometry=[Point(0.5, 0.5), POND, TRANSECT], crs=vor.crs
        )},
    )
    traces = choro.overlay_traces()
    wells, sites = _by_name(traces, "Wells"), _by_name(traces, "Sites")
    assert [trace.showlegend for trace in wells] == [True]
    # three geometry kinds, one legend entry, one legendgroup
    assert len(sites) == 3 and [trace.showlegend for trace in sites] == [True, False, False]
    assert {trace.legendgroup for trace in sites} == {"Sites"}
    colors = {wells[0].marker.color, sites[0].line.color}
    assert len(colors) == 2 and colors <= set(viz.PALETTE.categorical)
    # a layer keeps its colour on the next figure
    again = _choropleth_factory(vor, custom_zs=[1.0, 2.0], features={"Wells": _wells(vor.crs)})
    assert again.overlay_traces()[0].marker.color == wells[0].marker.color


def test_points_lines_and_polygons_become_the_right_traces(vor):
    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features={"Wells": _wells(vor.crs), "Pond": POND, "Transect": TRANSECT},
    )
    traces = {trace.name: trace for trace in choro.overlay_traces()}
    assert traces["Wells"].mode == "markers" and len(traces["Wells"].lon) == 2
    assert traces["Transect"].mode == "lines" and traces["Transect"].fill is None
    assert traces["Pond"].mode == "lines" and traces["Pond"].lon[0] == traces["Pond"].lon[-1]


def test_multi_part_geometry_draws_every_part_in_one_trace(vor):
    two_ponds = MultiPolygon([POND, Polygon([(1.2, 0.2), (1.8, 0.2), (1.8, 0.8), (1.2, 0.8)])])
    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features={"Ponds": two_ponds, "Pts": MultiPoint([(0.5, 0.5), (1.5, 0.5)])},
    )
    traces = {trace.name: trace for trace in choro.overlay_traces()}
    assert len(traces) == 2
    assert traces["Ponds"].lon.count(None) == 1  # two rings, one separator
    assert len(traces["Pts"].lon) == 2


def test_a_filled_polygon_adds_a_fill_beneath_its_outline(vor):
    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features=Feature(POND, name="Pond", fill=True, opacity=0.5, color="#ff0000"),
    )
    fill, outline = choro.overlay_traces()
    assert fill.fill == "toself" and fill.fillcolor == "rgba(255,0,0,0.500)"
    assert fill.showlegend is False and fill.hoverinfo == "skip"
    assert outline.fill is None and outline.line.color == "#ff0000" and outline.showlegend is True


def test_the_hover_lists_the_attribute_columns(vor):
    wells = _wells(vor.crs)
    default = feature_traces(resolve_features({"Wells": wells}, crs=vor.crs))[0]
    assert default.text[0] == "<b>W1</b><br>q: −100"
    assert default.hovertemplate == "%{text}<extra></extra>"

    only_name = feature_traces(resolve_features(Feature(wells, hover=False), crs=vor.crs))[0]
    assert only_name.text[1] == "<b>W2</b>"

    chosen = feature_traces(resolve_features(Feature(wells, label="q", hover=["name"]), crs=vor.crs))[0]
    assert chosen.text[1] == "<b>-250.5</b><br>name: W2"

    with pytest.raises(KeyError, match="available"):
        resolve_features(Feature(wells, hover=["depth"]), crs=vor.crs)
    with pytest.raises(KeyError, match="label column"):
        resolve_features(Feature(wells, label="depth"), crs=vor.crs)


def test_a_geometry_with_no_attributes_hovers_its_layer_name(vor):
    trace = feature_traces(resolve_features({"Pond": POND}, crs=vor.crs))[0]
    assert set(trace.text) == {"<b>Pond</b>"}


def test_a_source_without_a_crs_takes_the_grids(vor):
    bare = _wells(None)
    placed = resolve_features(bare, crs=vor.crs)[0].frame
    declared = resolve_features(_wells(vor.crs), crs=vor.crs)[0].frame
    assert placed.crs.to_epsg() == 4326
    assert list(placed.geometry.x) == pytest.approx(list(declared.geometry.x))

    with pytest.raises(ValueError, match="no CRS"):
        resolve_features(bare, crs=None)
    # Feature(crs=) declares it for a source that has none
    explicit = resolve_features(Feature(bare, crs=vor.crs), crs=None)[0].frame
    assert list(explicit.geometry.x) == pytest.approx(list(declared.geometry.x))


def test_a_source_with_a_crs_is_reprojected_and_an_override_is_refused(vor):
    declared = resolve_features(_wells(vor.crs), crs=vor.crs)[0].frame
    already_latlon = _wells(vor.crs).to_crs(epsg=4326)
    same = resolve_features(already_latlon, crs=vor.crs)[0].frame
    assert list(same.geometry.x) == pytest.approx(list(declared.geometry.x))
    with pytest.raises(ValueError, match="already declares"):
        resolve_features(Feature(_wells(vor.crs), crs="EPSG:4326"), crs=vor.crs)


def test_a_missing_file_and_an_unknown_layer_raise_at_the_call(vor, tmp_path):
    with pytest.raises(FileNotFoundError):
        _choropleth_factory(vor, custom_zs=[1.0, 2.0], features=tmp_path / "nope.gpkg")
    path = tmp_path / "sites.gpkg"
    _wells(vor.crs).to_file(path, layer="wells", driver="GPKG")
    with pytest.raises(ValueError, match="wells"):
        _choropleth_factory(vor, custom_zs=[1.0, 2.0], features=Feature(path, layer="ponds"))


def test_an_empty_source_only_warns(vor, caplog):
    with caplog.at_level("WARNING"):
        choro = _choropleth_factory(
            vor, custom_zs=[1.0, 2.0], features=gpd.GeoDataFrame(geometry=[], crs=vor.crs)
        )
    assert choro.overlay_traces() == []
    assert "no geometry" in caplog.text


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"opacity": 1.5}, "opacity"),
        ({"width": -1.0}, "width"),
        ({"hover": 3}, "hover"),
        ({"layer": "x"}, "layer= selects"),
    ],
)
def test_bad_feature_arguments_are_refused(kwargs, message):
    with pytest.raises((ValueError, TypeError), match=message):
        Feature(POND, **kwargs)


def test_an_unknown_features_spelling_is_refused(vor):
    with pytest.raises(TypeError, match="features= accepts"):
        _choropleth_factory(vor, custom_zs=[1.0, 2.0], features=42)


def test_features_ride_overlay_traces_into_mosaic_and_animate(vor):
    a = _choropleth_factory(vor, custom_zs=[1.0, 2.0], features={"Wells": _wells(vor.crs)})
    b = _choropleth_factory(vor, custom_zs=[2.0, 3.0], features={"Wells": _wells(vor.crs)})
    mosaic = viz.mosaic([("a", a), ("b", b)], ncols=2)
    assert sum(trace.type == "scattermap" for trace in mosaic.data) == 2
    animation = plot.animate([("a", a), ("b", b)])
    assert [trace.type for trace in animation.fig.data] == ["choroplethmap", "scattermap"]
    assert all(len(frame.data) == 2 for frame in animation.fig.frames)


def test_the_static_backend_draws_the_features(vor):
    choro = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0],
        features={"Wells": _wells(vor.crs), "Pond": Feature(POND, fill=True), "Transect": TRANSECT},
    )
    bare = _choropleth_factory(vor, custom_zs=[1.0, 2.0]).plot_mpl()
    figure = choro.plot_mpl()
    assert len(figure.axes[0].collections) > len(bare.axes[0].collections)
    assert [t.get_text() for t in figure.axes[0].get_legend().get_texts()] == ["Wells", "Pond", "Transect"]
    hidden = _choropleth_factory(
        vor, custom_zs=[1.0, 2.0], features=Feature(POND, show_legend=False)
    ).plot_mpl()
    assert hidden.axes[0].get_legend() is None


def test_features_reach_every_scope(canonical_run):
    vor = canonical_run.vor
    centroid = vor.gdf_vorPolys.geometry.iloc[0].centroid
    wells = gpd.GeoDataFrame({"name": ["W1"]}, geometry=[Point(centroid.x, centroid.y)], crs=vor.crs)
    for picture in (
        plot.map(canonical_run, layer=0, features=wells),
        canonical_run.plot.map(layer=0, features=wells),
        vor.plot.map(values=list(range(vor.ncpl)), features=wells),
        canonical_run.hds.map(features=wells),
    ):
        assert [trace.type for trace in picture.fig.data] == ["choroplethmap", "scattermap"]
    assert "features" in NOUN_MAP_PARAMS and "locs" not in NOUN_MAP_PARAMS
    static = canonical_run.plot.map(layer=0, features=wells, backend="mpl")
    assert static.axes[0].get_legend() is not None


# ---------------------------------------------------------------------------
# the retired spelling
# ---------------------------------------------------------------------------


def test_locs_is_deprecated_but_still_drawn(vor):
    with pytest.warns(DeprecationWarning, match=r"myflopy.plot.map\(locs=\).*features="):
        choro = _choropleth_factory(vor, custom_zs=[1.0, 2.0], locs=_wells(vor.crs))
    assert [trace.name for trace in choro.overlay_traces()] == ["locations"]

    with pytest.warns(DeprecationWarning):
        both = _choropleth_factory(
            vor, custom_zs=[1.0, 2.0], locs=_wells(vor.crs), features={"Pond": POND}
        )
    assert [trace.name for trace in both.overlay_traces()] == ["Pond", "locations"]


def test_locs_reaches_the_warning_from_the_verbs_too(vor):
    with pytest.warns(DeprecationWarning, match="features="):
        plot.map(vor, values=[1, 2], locs=_wells(vor.crs))
    with pytest.warns(DeprecationWarning, match="features="):
        vor.plot.map(values=[1, 2], locs=_wells(vor.crs))


def test_locs_is_no_longer_named_by_any_map_signature():
    """A deprecated spelling must not echo anywhere an editor reads."""

    for function in (plot.map, ModelPlots.map, GridPlots.map, _choropleth_factory):
        params = inspect.signature(function).parameters
        assert "locs" not in params, function.__qualname__
        assert "features" in params, function.__qualname__
    from myflopy.modflow.utils.datatypes.choros import Choro

    assert "locs" not in inspect.signature(Choro.__init__).parameters
    assert not any(name.startswith("locs") or "_locs" in name for name in vars(Choro))


def test_a_map_can_be_built_without_features_or_a_warning(vor):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        choro = _choropleth_factory(vor, custom_zs=[1.0, 2.0])
    assert choro.features == [] and choro.overlay_traces() == []


# ---------------------------------------------------------------------------
# the review pass (ledger 252/254 addenda, 255)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "module",
    ["myflopy.modflow.utils.datatypes.choros", "myflopy.modflow.utils.datatypes.features"],
)
def test_the_choropleth_and_feature_modules_import_cold(module):
    """A fresh interpreter, nothing primed. The suite never sees this: conftest
    imports the canonical model first, which imports the grid package, which is
    exactly the edge the cycle went through (grid/__init__ -> plotting -> choros
    -> features -> grid.selection -> grid/__init__)."""

    import subprocess
    import sys

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    completed = subprocess.run(
        [sys.executable, "-c", f"import {module}"],
        capture_output=True, text=True,
        env={**os.environ, "PYTHONPATH": os.path.join(root, "src")},
    )
    assert completed.returncode == 0, completed.stderr[-800:]


def test_a_bare_grid_log_map_hovers_the_real_value(vor):
    """The LEGACY hover (no hover_spec: a grid map) listed the log10 exponent."""

    trace = plot.map(vor, values=[0.001, 100.0], logscale=True).get_choropleth()
    flat = " ".join(str(row) for row in trace.customdata)
    assert "0.001" in flat and "100" in flat
    assert "-3" not in flat


def test_features_draw_the_active_geometry_column(vor):
    """After ``set_geometry`` a stale ``geometry`` column remains; the active one
    is what is drawn, on both backends, and neither column is a hover attribute."""

    wells = _wells(vor.crs)
    wells["buffer"] = wells.geometry.buffer(0.1)
    moved = wells.set_geometry("buffer")
    layer = resolve_features({"Buffers": moved}, crs=vor.crs)[0]
    assert layer.hover_columns == ("q",)
    trace = feature_traces([layer])[0]
    assert trace.mode == "lines"                      # polygons, not the stale points
    assert min(abs(x) for x in trace.lon if x is not None) > 100  # reprojected, not raw feet


def test_ani_carries_the_features_into_every_frame(canonical_run):
    vor = canonical_run.vor
    centroid = vor.gdf_vorPolys.geometry.iloc[0].centroid
    wells = gpd.GeoDataFrame({"name": ["W1"]}, geometry=[Point(centroid.x, centroid.y)], crs=vor.crs)
    animation = canonical_run.plot.map(layer=0, features=wells).ani
    assert [trace.type for trace in animation.data] == ["choroplethmap", "scattermap"]
    assert all([trace.type for trace in frame.data] == ["choroplethmap", "scattermap"]
               for frame in animation.frames)


def test_animate_pools_log_frames_onto_one_labelled_axis(vor):
    a = _log_map(vor)
    b = _choropleth_factory(vor, custom_zs=[0.01, 10.0], logscale=True)
    figure = plot.animate([("a", a), ("b", b)]).fig
    coloraxis = figure.layout.coloraxis
    assert (coloraxis.cmin, coloraxis.cmax) == pytest.approx((-3.0, 2.0))
    assert coloraxis.colorbar.ticktext == ("0.001", "0.01", "0.1", "1", "10", "100")
    assert all(trace.coloraxis == "coloraxis" for trace in [figure.data[0], *(f.data[0] for f in figure.frames)])
    with pytest.raises(ValueError, match="same logscale"):
        plot.animate([("a", a), ("b", _choropleth_factory(vor, custom_zs=[1.0, 2.0]))]).fig


def test_the_noun_tier_static_composers_read_real_units(vor):
    """``SpatialView._mpl_mosaic`` / ``_mpl_animation`` drew ``choro.zs`` (log10)
    through a hard-coded linear viridis, so a static log mosaic read -3 ... 2."""

    from matplotlib.colors import LogNorm

    from myflopy.modflow.mf6.package_plotting import SpatialView, _mpl_shared_scale

    a = _log_map(vor, colorscale="Viridis")
    b = _choropleth_factory(vor, custom_zs=[0.01, 10.0], logscale=True, colorscale="Viridis")
    arrays, cmap, norm = _mpl_shared_scale([a, b], diff=False)
    assert isinstance(norm, LogNorm) and (norm.vmin, norm.vmax) == pytest.approx((0.001, 100.0))
    assert list(arrays[0]) == pytest.approx([0.001, 100.0])
    assert cmap.name == "plotly_Viridis"
    with pytest.raises(ValueError, match="same logscale"):
        _mpl_shared_scale([a, _choropleth_factory(vor, custom_zs=[1.0, 2.0])], diff=False)

    class _Stub:
        def _spatial_is_diff(self):
            return False

        def _spatial_value_label(self):
            return "K (m/d)"

    figure = SpatialView._mpl_mosaic(_Stub(), [("a", a), ("b", b)], ncols=2, title="t")
    assert figure.axes[-1].get_yscale() == "log"
    assert figure.axes[0].collections[0].get_clim() == pytest.approx((0.001, 100.0))
    animation = SpatialView._mpl_animation(_Stub(), [("a", a), ("b", b)], title="t")
    assert animation._fig.axes[0].collections[0].get_clim() == pytest.approx((0.001, 100.0))


@pytest.mark.parametrize("kwargs", [{"zmin": float("inf")}, {"zmax": float("nan")}])
def test_a_non_finite_log_limit_is_refused_at_the_call(vor, kwargs):
    with pytest.raises(ValueError, match="finite and positive"):
        _log_map(vor, **kwargs)


def test_the_static_backend_refuses_a_non_positive_override_on_a_log_map(vor):
    with pytest.raises(ValueError, match="vmin=0"):
        _log_map(vor).plot_mpl(vmin=0)


def test_decade_labels_are_written_out_in_full():
    from myflopy.viz import _decade_label, log_decade_colorbar

    assert _decade_label(6, "") == "1000000"
    assert _decade_label(-5, "") == "0.00001"
    assert _decade_label(0.5, "x") == "3.16x"
    assert log_decade_colorbar(float("-inf"), 2.0) == {}


def test_layer_names_never_collide_and_the_retired_layer_never_clobbers_a_key(vor):
    from myflopy.modflow.utils.datatypes.features import _unique, combine_features

    names = _unique(["a", "a", "a 2"])
    assert names[:2] == ["a", "a 2"] and len(set(names)) == 3   # the explicit "a 2" is re-suffixed
    combined = combine_features({"locations": POND}, Feature(_wells(vor.crs), name="locations"))
    assert list(combined) == ["locations", "locations 2"]
    with pytest.warns(DeprecationWarning):
        choro = _choropleth_factory(
            vor, custom_zs=[1.0, 2.0], features={"locations": POND}, locs=_wells(vor.crs)
        )
    assert [trace.name for trace in choro.overlay_traces()] == ["locations", "locations 2"]


def test_hover_values_are_escaped_and_missing_values_are_blank(vor):
    frame = gpd.GeoDataFrame(
        {"name": ["D<b>x</b>", None], "flag": [True, False], "n": [pd.NA, 3]},
        geometry=[Point(0.5, 0.5), Point(1.5, 0.5)], crs=vor.crs,
    ).astype({"n": "Int64"})
    trace = feature_traces(resolve_features({"Pts": frame}, crs=vor.crs))[0]
    assert trace.text[0] == "<b>D&lt;b&gt;x&lt;/b&gt;</b><br>flag: True<br>n: "
    assert trace.text[1] == "<b>Pts</b><br>flag: False<br>n: 3"


def test_an_empty_geometry_row_is_dropped_without_a_warning(vor):
    frame = gpd.GeoDataFrame(geometry=[Point(0.5, 0.5), Polygon(), None], crs=vor.crs)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        layer = resolve_features(frame, crs=vor.crs)[0]
    assert len(layer.frame) == 1


def test_a_color_only_one_backend_parses_is_refused_at_construction():
    with pytest.raises(ValueError, match="both backends"):
        Feature(POND, color="hsl(120,50%,50%)")
    assert Feature(POND, color="rgb(1,2,3)").color == "rgb(1,2,3)"


def test_the_static_backend_explodes_nested_collections(vor):
    from shapely.geometry import GeometryCollection

    nested = GeometryCollection([Point(0.5, 0.5), GeometryCollection([Point(1.5, 0.5), POND])])
    figure = _choropleth_factory(vor, custom_zs=[1.0, 2.0], features=nested).plot_mpl()
    bare = _choropleth_factory(vor, custom_zs=[1.0, 2.0]).plot_mpl()
    # two marker collections' worth of points plus a polygon patch beyond the cells
    assert len(figure.axes[0].collections) >= len(bare.axes[0].collections) + 2


def test_the_deprecation_warning_is_attributed_to_the_caller(vor):
    with pytest.warns(DeprecationWarning) as record:
        vor.plot.map(values=[1, 2], locs=_wells(vor.crs))
    warned = [w for w in record if issubclass(w.category, DeprecationWarning)]
    assert warned and warned[0].filename == __file__
