from __future__ import annotations

from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point, box

from myflopy.layers import (
    Array,
    Clamp,
    Flat,
    Isopach,
    LayerBuildResult,
    LayerStack,
    Max,
    Min,
)
from myflopy.surfaces import LayerSurfaces, Surface

_PLANE = dict(xs=[-1, 3, -1, 3], ys=[-1, -1, 1, 1], zs=[15, -5, 15, -5])  # z = 10 - 5x


def _fake_vor(n=3):
    """Point-backed grid stand-in (enough for flat/thickness/points surfaces)."""
    pts = [Point(float(i), 0.0) for i in range(n)]
    gdf = gpd.GeoDataFrame({"geometry": pts}, crs="EPSG:2927")
    ns = SimpleNamespace(
        gdf_vorPolys=gdf,
        centroids=([p.x for p in pts], [p.y for p in pts]),
        crs="EPSG:2927",
        ncpl=n,
    )
    ns.get_disv_gridprops = lambda: {
        "ncpl": n,
        "nvert": n + 1,
        "vertices": [[i, float(i), 0.0] for i in range(n + 1)],
        "cell2d": [[i, float(i), 0.0, 1, i] for i in range(n)],
    }
    return ns


def _poly_vor(polys, crs="EPSG:2927"):
    gdf = gpd.GeoDataFrame({"geometry": list(polys)}, crs=crs)
    c = gdf.geometry.centroid
    return SimpleNamespace(
        gdf_vorPolys=gdf, centroids=(c.x.to_numpy(), c.y.to_numpy()), crs=crs, ncpl=len(gdf)
    )


def _write_raster(path, data, *, nodata=None, west=0.0, north=4.0, px=1.0, crs="EPSG:2927"):
    import rasterio
    from rasterio.transform import from_origin

    data = np.asarray(data, dtype="float64")
    h, w = data.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=h, width=w, count=1, dtype="float64",
        crs=crs, transform=from_origin(west, north, px, px), nodata=nodata,
    ) as dst:
        dst.write(data, 1)


def test_layerstack_basic_build():
    vor = _fake_vor(3)
    result = (
        LayerStack(vor, top=Flat(100))
        .add("a", thickness=10)
        .add("b", thickness=20)
        .build(reconcile=False)
    )
    assert isinstance(result, LayerBuildResult)
    assert result.nlay == 2
    assert result.names == ["a", "b"]
    assert list(result.top) == [100, 100, 100]
    assert list(result.botm[0]) == [90, 90, 90]
    assert list(result.botm[1]) == [70, 70, 70]
    assert (result.idomain == 1).all()
    assert result.length_units == "feet" and result.time_units == "days"


def test_facade_matches_layersurfaces_engine():
    vor = _fake_vor(3)
    facade = (
        LayerStack(vor, top=Flat(100)).add("a", thickness=10).add("b", bottom=Flat(50))
        .build(reconcile=False)
    )
    top, botm = LayerSurfaces(
        [Surface.flat(100), Surface.constant_thickness(10), Surface.flat(50)]
    ).top_botm(vor, reconcile=False)
    assert np.allclose(facade.top, top)
    assert np.allclose(facade.botm, botm)


def test_layerstack_coerces_path_and_number():
    vor = _fake_vor(2)
    result = LayerStack(vor, top=200).add("x", bottom=100).build(reconcile=False)
    assert list(result.top) == [200, 200]
    assert list(result.botm[0]) == [100, 100]


def test_coerce_rejects_bad_type():
    with pytest.raises(TypeError):
        LayerStack(_fake_vor(2), top=object())


def test_add_requires_exactly_one_of_bottom_or_thickness():
    vor = _fake_vor(2)
    with pytest.raises(ValueError, match="exactly one"):
        LayerStack(vor, top=Flat(10)).add("x")
    with pytest.raises(ValueError, match="exactly one"):
        LayerStack(vor, top=Flat(10)).add("x", bottom=Flat(5), thickness=5)


def test_duplicate_layer_name_raises():
    s = LayerStack(_fake_vor(2), top=Flat(10)).add("x", thickness=1)
    with pytest.raises(ValueError, match="already exists"):
        s.add("x", thickness=1)


def test_per_layer_pinch_policy_override():
    vor = _fake_vor(3)
    base = LayerStack(vor, top=Flat(10)).add("a", bottom=Surface.from_points(**_PLANE))
    assert base.build(default_min_thickness=1.0, reconcile=False).idomain[0].tolist() == [-1, 1, 1]

    override = LayerStack(vor, top=Flat(10)).add(
        "a", bottom=Surface.from_points(**_PLANE), pinch="inactive"
    )
    assert override.build(default_min_thickness=1.0, reconcile=False).idomain[0].tolist() == [0, 1, 1]


def test_per_layer_min_thickness_override():
    vor = _fake_vor(3)  # thickness [0, 5, 10]
    # default threshold 1.0 -> only cell0 thin; raise threshold for this layer to 6
    stack = LayerStack(vor, top=Flat(10)).add(
        "a", bottom=Surface.from_points(**_PLANE), min_thickness=6.0
    )
    idom = stack.build(reconcile=False).idomain[0]
    assert idom.tolist() == [-1, -1, 1]  # 0 and 5 are < 6, 10 is not


def test_replace_insert_remove_by_name():
    vor = _fake_vor(2)
    s = LayerStack(vor, top=Flat(100)).add("a", thickness=10).add("b", thickness=10)
    s.replace("a", thickness=40)
    assert list(s.build(reconcile=False).botm[0]) == [60, 60]
    s.insert_below("a", "a2", thickness=5)
    assert s.names == ["a", "a2", "b"]
    s.remove("a2")
    assert s.names == ["a", "b"]


def test_remove_unknown_layer_raises():
    s = LayerStack(_fake_vor(2), top=Flat(10)).add("a", thickness=1)
    with pytest.raises(KeyError, match="no layer named"):
        s.remove("missing")


def test_units_conversion_meters_to_feet(tmp_path):
    path = tmp_path / "dem_m.tif"
    _write_raster(path, np.full((4, 4), 100.0))  # 100 metres everywhere
    vor = _poly_vor([box(0, 1, 3, 4)])
    result = (
        LayerStack(vor, top=Surface.raster(path, units="meters"), length_units="feet")
        .add("x", bottom=Flat(0))
        .build(reconcile=False)
    )
    assert result.top[0] == pytest.approx(100 * 3.280839895, rel=1e-6)


def test_units_on_relative_surface_raises(tmp_path):
    vor = _fake_vor(2)
    bad = Surface.constant_thickness(10)
    bad = type(bad)(**{**bad.__dict__, "units": "meters"})  # force units on a relative kind
    with pytest.raises(ValueError, match="absolute surfaces"):
        LayerSurfaces([Surface.flat(10), bad]).sample(vor, reconcile=False, length_units="feet")


def test_to_disv_sets_length_units_and_idomain():
    vor = _fake_vor(3)
    stack = LayerStack(vor, top=Flat(10), length_units="feet").add(
        "a", bottom=Surface.from_points(**_PLANE)
    )
    spec = stack.to_disv(reconcile=False, default_min_thickness=1.0)
    assert spec.options["length_units"] == "FEET"
    assert spec.options["idomain"][0].tolist() == [-1, 1, 1]


def test_report_is_a_string():
    vor = _fake_vor(3)
    result = LayerStack(vor, top=Flat(10)).add("clay", thickness=5).build(reconcile=False)
    report = result.report()
    assert "LayerStack" in report and "clay" in report


def test_cache_status_reports_only_derived_surfaces(tmp_path):
    vor = _fake_vor(2)
    contour = Surface.from_contours(tmp_path / "c.gpkg", z="e", out=tmp_path / "c.interp.tif")
    stack = LayerStack(vor, top=Flat(10)).add("bed", bottom=contour).add("x", thickness=5)
    # Only the contour layer is derived; no GRASS is invoked by status().
    assert stack.cache_status() == {"bed": "missing"}


def test_facade_isopach_layer():
    vor = _fake_vor(3)
    result = (
        LayerStack(vor, top=Flat(100))
        .add("a", thickness=Isopach(Flat(20)))  # 20-thick layer -> bottom 80
        .build(reconcile=False)
    )
    assert list(result.botm[0]) == [80, 80, 80]


def test_facade_surface_algebra_bottom():
    vor = _fake_vor(3)
    # bottom = max(flat 40, flat 60) = 60 everywhere
    result = (
        LayerStack(vor, top=Flat(100))
        .add("a", bottom=Max(Flat(40), Flat(60)))
        .add("b", bottom=Clamp(Flat(10), lower=20, upper=80))  # -> 20
        .build(reconcile=False)
    )
    assert list(result.botm[0]) == [60, 60, 60]
    assert list(result.botm[1]) == [20, 20, 20]


def test_thickness_map_is_a_picture_over_matplotlib_axes():
    """8.5a: `preview()`/`thickness_map()` became `plot.map()`.

    It stays Matplotlib and basemap-free -- a stack under construction is often
    on synthetic coordinates, where a web basemap lands in the ocean -- but it is
    now a Picture, so `.show()`/`.save()`/`.html()` work like everywhere else.
    """
    from myflopy.viz import MplPicture

    vor = _poly_vor([box(0, 0, 1, 1), box(1, 0, 2, 1)])
    picture = LayerStack(vor, top=Flat(100)).add("a", thickness=10).build(reconcile=False).plot.map()
    assert isinstance(picture, MplPicture)
    assert hasattr(picture.axes, "set_title")   # the Matplotlib axes, as before
    assert picture.axes is picture.axes         # idempotent, per the contract


def test_algebra_aliases_build_the_right_kinds():
    assert Min(Flat(1), Flat(2)).kind == "min"
    assert Max(Flat(1), Flat(2)).kind == "max"
    assert Clamp(Flat(1), lower=0).kind == "clamp"
    assert Isopach(Flat(5)).kind == "isopach"


# --- views (need a real Voronoi grid; built once for the module) -------------
@pytest.fixture(scope="module")
def real_vor():
    import tempfile

    import myflopy as mf

    ws = tempfile.mkdtemp(prefix="vorgrid_")
    tri = mf.TriangleGrid(model_ws=ws, angle=30)
    tri.set_domain_rectangle(x_dist=400, y_dist=300, origin=(0, 0))
    tri.build()
    return mf.VoronoiGridPlus(tri)


def test_result_vertex_grid(real_vor):
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    vg = res.vertex_grid()
    assert vg.nlay == 1 and vg.ncpl == real_vor.ncpl


def test_section_is_a_picture_over_matplotlib_axes(real_vor):
    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(30)).add("b", bottom=Flat(10))
        .build()
    )
    assert hasattr(res.plot.section(y=150).axes, "set_title")            # y= shorthand
    assert hasattr(res.plot.section(x=200, legend=False).axes, "set_title")  # x= shorthand
    assert hasattr(res.plot.section(line=[(0, 150), (400, 150)]).axes, "set_title")


def test_section_defaults_to_center_line(real_vor):
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    assert hasattr(res.plot.section().axes, "set_title")  # no line/x/y -> W-E centre


def test_thickness_map_totals_or_one_layer(real_vor):
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    assert hasattr(res.plot.map().axes, "set_title")
    assert hasattr(res.plot.map(layer="a").axes, "set_title")


def test_thickness_map_can_opt_into_the_georeferenced_choropleth(real_vor):
    """`basemap=True` routes through the shared choropleth instead.

    The default is basemap-free on purpose (synthetic coordinates), but a grid
    that really is georeferenced should get the same map everything else draws.
    """
    from myflopy.viz import Fig

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    picture = res.plot.map(basemap=True)
    assert isinstance(picture.fig, Fig)
    assert len(picture.fig.data) >= 1


def test_views_is_gone_compose_with_mosaic_instead(real_vor):
    """`views()` displayed four pictures and returned None -- a Layer-3 concern
    fused into Layer 1. Composing is `myflopy.plot.mosaic`, which takes any
    pictures you like rather than a fixed four."""
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    assert not hasattr(res, "views")
    assert not hasattr(res, "thickness_map")
    assert not hasattr(res, "cross_section")
    assert not hasattr(res, "surface_3d")


def test_surface_3d_returns_plotly_figure(real_vor):
    import plotly.graph_objects as go

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    fig_top = res.plot.surface("top", resolution=30).fig
    assert isinstance(fig_top, go.Figure)
    assert len(fig_top.data) == 1                       # single surface
    assert isinstance(res.plot.surface("a", resolution=30).fig, go.Figure)


def test_surface_names_lists_top_and_each_bottom(real_vor):
    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(30)).add("b", bottom=Flat(10))
        .build()
    )
    assert res.surface_names == ["top", "a", "b"]


def test_surface_3d_accepts_a_list_of_layers(real_vor):
    import plotly.graph_objects as go

    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(30)).add("b", bottom=Flat(10))
        .build()
    )
    fig = res.plot.surface(["top", "b"], resolution=30).fig
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 2
    assert {tr.name for tr in fig.data} == {"top", "b"}
    assert fig.layout.showlegend is True               # multi -> legend on


def test_surface_3d_all_draws_every_surface(real_vor):
    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(30)).add("b", bottom=Flat(10))
        .build()
    )
    fig = res.plot.surface("all", resolution=30).fig
    assert len(fig.data) == len(res.surface_names)      # top + a + b == 3


def test_surface_3d_color_by_elevation_shares_one_colorbar(real_vor):
    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(30)).add("b", bottom=Flat(10))
        .build()
    )
    fig = res.plot.surface("all", color_by="elevation", resolution=30).fig
    shown = [tr for tr in fig.data if tr.showscale]
    assert len(shown) == 1                              # exactly one colorbar


def test_surface_3d_unknown_name_raises(real_vor):
    import pytest

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    with pytest.raises(KeyError):
        res.plot.surface("nope", resolution=30).fig  # lazy: raises on assembly


def test_surface_3d_flat_single_layer_stays_visible(real_vor):
    # A perfectly flat go.Surface has zero vertical extent and vanishes in some
    # WebGL viewers. A flat layer is given faint relief, a non-degenerate colour
    # range, and a real z-axis so the sheet always draws (solid-filled).
    import numpy as np

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    fig = res.plot.surface("a", resolution=30).fig          # flat -> solid fill + relief
    tr = fig.data[0]
    z = np.asarray(tr.z, dtype=float)
    assert len(tr.colorscale) == 2                        # solid 2-stop scale, not elevation
    assert tr.cmin < tr.cmax                              # non-degenerate colour range
    assert float(np.nanmax(z) - np.nanmin(z)) > 0         # real vertical extent -> renders
    zr = fig.layout.scene.zaxis.range
    assert zr is not None and zr[1] > zr[0]               # explicit, widened z-axis


def test_surface_3d_height_defaults_to_fill_container(real_vor):
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    assert res.plot.surface("top", resolution=20).fig.layout.height is None   # fills page
    assert res.plot.surface("top", resolution=20, height=700).fig.layout.height == 700


def test_surface_writes_standalone_html_through_the_picture(real_vor, tmp_path):
    """`html_path=`/`browser=` were Layer-3 concerns baked into the builder.
    Writing the file is `.html(path)`, the same call on every picture."""
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    out = tmp_path / "s.html"
    picture = res.plot.surface("a", resolution=20)
    assert picture.html(out, include_plotlyjs=True) == out
    assert out.exists() and out.stat().st_size > 0       # standalone HTML written
    assert hasattr(picture.fig, "data")                   # and the figure is still there


def test_resolve_layer_indices_by_name_and_index(real_vor):
    import pytest

    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(40))
        .add("b", bottom=Flat(20)).add("c", bottom=Flat(5)).build()
    )
    assert res._resolve_layer_indices(None) == [0, 1, 2]
    assert res._resolve_layer_indices("all") == [0, 1, 2]
    assert res._resolve_layer_indices("b") == [1]
    assert res._resolve_layer_indices(["c", "a"]) == [0, 2]   # sorted + unique
    assert res._resolve_layer_indices([0, 2, 2]) == [0, 2]
    with pytest.raises(KeyError):
        res._resolve_layer_indices("nope")
    with pytest.raises(IndexError):
        res._resolve_layer_indices(5)


def test_cross_section_uses_shared_layered_renderer(real_vor):
    # The pre-model facade and the model-aware plotter share one renderer; this
    # exercises that core directly on a VertexGrid (no built model needed).
    from myflopy.modflow.mf6.cross_section_plotting import plot_layered_cross_section

    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(30)).add("b", bottom=Flat(10))
        .build()
    )
    fig, ax = plot_layered_cross_section(
        res.vertex_grid(), line=[(0, 150), (400, 150)], layer_labels=res.names,
    )
    assert hasattr(ax, "set_title")
    assert len(ax.get_legend().get_texts()) == res.nlay  # one legend entry per layer


def test_surface_trace_builds_one_go_surface(real_vor):
    import numpy as np
    import plotly.graph_objects as go

    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    isurf = InterpolatedSurface(
        xs=np.asarray(real_vor.centroids[0]), ys=np.asarray(real_vor.centroids[1]),
        zs=np.asarray(res.top), surf_type="lyr", resolution=30,
    )
    assert isinstance(isurf.surface_trace(), go.Surface)
    assert isinstance(isurf.surface_trace(showlegend=True, name="x"), go.Surface)


def test_the_3d_grid_is_a_scene_that_writes_nothing_by_itself(real_vor, tmp_path, monkeypatch):
    """8.5b: `vtk_3d()` became `plot.grid(backend="vtk")`.

    It used to write a standalone HTML file on EVERY call -- into `Path.cwd()`
    when given no path -- and hand back an `IFrame` pointing at it. That litter
    was real enough that `.gitignore` carried a line naming the file. Building a
    picture is now Layer 1 only; writing it is `.html(path)`.
    """
    from myflopy.viz import VtkScene

    monkeypatch.chdir(tmp_path)
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    scene = res.plot.grid()

    assert isinstance(scene, VtkScene)
    assert list(tmp_path.glob("*.html")) == []      # nothing written just by building
    try:
        out = scene.html(tmp_path / "v.html")
        assert out.exists() and out.stat().st_size > 0
    finally:
        scene.scene.close()


def test_the_3d_grid_shows_a_subset_of_layers(real_vor, tmp_path):
    from myflopy.viz import VtkScene

    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(40))
        .add("b", bottom=Flat(20)).add("c", bottom=Flat(5)).build()
    )
    one = res.plot.grid("b")            # single, by name
    some = res.plot.grid(["a", "c"])    # subset, by name
    try:
        assert isinstance(one, VtkScene) and isinstance(some, VtkScene)
        assert one.html(tmp_path / "one.html").stat().st_size > 0
        assert some.html(tmp_path / "some.html").stat().st_size > 0
    finally:
        one.scene.close()
        some.scene.close()


def test_the_3d_grid_rejects_an_unknown_backend(real_vor):
    """`backend=` is a renderer switch with exactly two settings; a typo should
    say so rather than silently drawing the wrong thing."""
    import pytest

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    with pytest.raises(ValueError, match="'vtk' or 'plotly'"):
        res.plot.grid(backend="opengl")


def test_stack_level_view_oneliners(real_vor):
    """`stack.plot` builds with defaults, then hands over the result's namespace.

    The one-liner shape survives 8.5a -- you still never call `.build()` for a
    quick look -- it is just spelled through the same verbs as everything else.
    """
    import plotly.graph_objects as go

    stack = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20))
    assert hasattr(stack.plot.section(y=150).axes, "set_title")  # builds, then views
    assert isinstance(stack.plot.surface().fig, go.Figure)
    assert hasattr(stack.plot.map().axes, "set_title")


# --- QC / validation ---------------------------------------------------------
def test_qc_clean_stack_is_ok(real_vor):
    res = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).add("b", bottom=Flat(5))
        .build()
    )
    r = res.qc()
    assert r.ok
    assert r.nan_active_cells == 0 and r.isolated_active == []
    assert r.n_active_components == 1
    assert "OK" in str(r)


def test_qc_detects_and_prunes_isolated_active_cells(real_vor):
    import dataclasses

    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    idom = res.idomain.copy()
    victim = 0
    for j in real_vor.adjacent_cells_idx[victim]:   # cut the cell off from its neighbours
        idom[0, int(j)] = 0
    res = dataclasses.replace(res, idomain=idom)
    r = res.qc()
    assert not r.ok
    assert (0, victim) in r.isolated_active
    pruned = res.prune_isolated()
    assert pruned.idomain[0, victim] == 0           # isolated cell deactivated
    assert pruned.qc().isolated_active == []        # and none remain


def test_qc_flags_nan_bounded_active_cells(real_vor):
    arr = np.full(real_vor.ncpl, 20.0)
    arr[1] = np.nan                                  # one cell has no bottom coverage
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Array(arr)).build(reconcile=False)
    r = res.qc()
    assert not r.ok
    assert r.nan_active_cells >= 1 and r.nan_botm[0] >= 1


def test_validate_raises_on_problems(real_vor):
    arr = np.full(real_vor.ncpl, 20.0)
    arr[0] = np.nan
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Array(arr)).build(reconcile=False)
    with pytest.raises(ValueError, match="failed QC"):
        res.validate()


def test_stack_qc_adds_reconcile_diagnostics(real_vor):
    # bottom above top forces reconcile to push it down.
    r = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(60)).qc()
    assert r.reconcile_adjusted is not None
    assert r.reconcile_adjusted[0] > 0 and r.reconcile_max_shift[0] > 0


# --- reading existing MODFLOW arrays (from_array / from_modflow) --------------
def test_from_array_round_trips_on_same_grid(real_vor):
    arr = np.linspace(10, 30, real_vor.ncpl)
    res = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Array(arr)).build(reconcile=False)
    assert np.allclose(res.botm[0], arr)


def test_from_array_wrong_length_raises(real_vor):
    stack = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Array(np.zeros(real_vor.ncpl + 1)))
    with pytest.raises(ValueError, match="from_array surface has"):
        stack.build(reconcile=False)


def test_from_modflow_round_trips_existing_arrays_verbatim(real_vor):
    base = (
        LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).add("b", bottom=Flat(5))
        .build()
    )
    vg = base.vertex_grid()                          # stand-in for an existing MF6 model grid
    stack = LayerStack.from_modflow(real_vor, vg, names=["a", "b"], resample=False)
    res = stack.build()
    assert np.allclose(res.top, base.top, equal_nan=True)
    assert np.allclose(res.botm, base.botm, equal_nan=True)


def test_from_modflow_resamples_and_defaults_names(real_vor):
    base = LayerStack(real_vor, top=Flat(50)).add("a", bottom=Flat(20)).build()
    stack = LayerStack.from_modflow(real_vor, base.vertex_grid())   # resample=True default
    assert stack.names == ["layer1"]
    res = stack.build()
    assert np.allclose(res.top, base.top, atol=1e-6, equal_nan=True)
