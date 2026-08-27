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
    Raster,
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
    contour = Surface.from_contours(
        tmp_path / "c.gpkg", z="e", out=tmp_path / "c.interp.tif",
        region_vector=tmp_path / "domain.gpkg",
    )
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


# --- the grid can arrive after the layering (2026-08-27) --------------------- #
#
# `LayerSurfaces` was always declarable before the grid; `LayerStack` was not,
# because it took `vor` in its constructor. That forced a choice between the
# ordering (declare layering first, project-first) and the per-layer control
# (`thickness=`, per-layer `min_thickness`/`pinch`) that only the facade offers.
# `vor` is now optional, and `build`/`qc`/`to_disv` accept one -- mirroring the
# override `to_disv` already had.

def test_a_stack_can_be_declared_with_no_grid():
    """The layering is lazy: `add` resolves nothing, so no grid is needed yet."""

    stack = (
        LayerStack(top=Flat(200))
        .add("sand", thickness=40.0, min_thickness=2.0, pinch="inactive")
        .add("clay", thickness=25.0, pinch="passthrough")
    )
    assert stack.vor is None
    assert stack.names == ["sand", "clay"]


@pytest.mark.parametrize("consume", ["build", "to_disv"])
def test_the_grid_can_be_supplied_at_use_time(consume):
    vor = _fake_vor(3)
    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    assert getattr(stack, consume)(vor) is not None


def test_qc_takes_the_grid_too(real_vor):
    """Separate from the others: `qc`'s isolated-cell check needs a real grid
    (`adjacent_cells_idx`), which the lightweight `_fake_vor` does not carry."""

    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    assert "1 layers" in str(stack.qc(real_vor))


def test_for_grid_binds_a_copy_and_leaves_the_original_deferred():
    """One declaration, several grids -- the layering is the expensive part."""

    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    coarse, fine = _fake_vor(2), _fake_vor(5)
    a, b = stack.for_grid(coarse), stack.for_grid(fine)

    assert stack.vor is None, "for_grid must not mutate the declaration"
    assert a.vor is coarse and b.vor is fine
    assert a.build().top.size == 2
    assert b.build().top.size == 5


def test_a_bound_stack_still_takes_the_old_positional_form():
    """`LayerStack(vor, top)` is in every notebook and doc; it must keep working."""

    vor = _fake_vor(3)
    assert LayerStack(vor, Flat(200)).add("a", thickness=10.0).build().nlay == 1
    assert LayerStack(vor, top=Flat(200)).add("a", thickness=10.0).build().nlay == 1


@pytest.mark.parametrize(
    ("consume", "pattern"),
    [("build", r"stack\.build\(vor\)"), ("qc", r"stack\.qc\(vor\)"),
     ("to_disv", r"stack\.to_disv\(vor\)")],
)
def test_using_a_deferred_stack_without_a_grid_names_the_fix(consume, pattern):
    """Otherwise this surfaces as `NoneType has no attribute ncpl`, three frames down."""

    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    with pytest.raises(ValueError, match=pattern):
        getattr(stack, consume)()


def test_plot_without_a_grid_falls_back_rather_than_refusing():
    """Superseded the "name both spellings" error one commit later: `.plot` now
    draws on a draft grid instead of raising. It can still only fail for want of
    an EXTENT -- these surfaces are `Flat`, so there is nothing to infer -- and
    that message points at the extent, not at binding a grid."""

    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    with pytest.raises(ValueError, match="no surface in this stack carries"):
        stack.plot


def test_passing_a_surface_as_the_grid_is_caught_at_construction():
    """`LayerStack(ground)` reads as the deferred form and is not; without this
    guard `ground` binds to `vor` and fails much later, somewhere else."""

    with pytest.raises(TypeError, match="passed the top surface positionally"):
        LayerStack(Flat(200))
    with pytest.raises(TypeError, match="first argument is the grid"):
        LayerStack(Flat(200), Flat(100))
    with pytest.raises(TypeError, match="needs a `top` surface"):
        LayerStack()


# --- the documented reconcile options are the accepted ones ----------------- #

@pytest.mark.parametrize("value", ["bottom", "top", True, False, None])
def test_every_documented_reconcile_option_is_accepted(value):
    """The docstring on `build` lists these five; `_reconcile_args` must take them."""

    vor = _fake_vor(3)
    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    assert stack.build(vor, reconcile=value) is not None


def test_an_undocumented_reconcile_value_is_rejected():
    vor = _fake_vor(3)
    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    with pytest.raises(ValueError, match="reconcile must be"):
        stack.build(vor, reconcile="sideways")


def test_the_docstring_lists_exactly_the_accepted_options():
    """Catches the drift this test file exists to prevent: a new option added to
    `_reconcile_args` and never documented, or a documented one that never
    worked. Compared against the parser's own source, not a hand-copied list."""

    import inspect

    from myflopy.layers import _reconcile_args

    accepted = set()
    for token in ("'bottom'", "'top'", "True", "False", "None"):
        assert token.strip("'") in inspect.getsource(_reconcile_args), token
        accepted.add(token)

    doc = inspect.getdoc(LayerStack.build) or ""
    header = next(ln for ln in doc.splitlines() if ln.startswith("reconcile :"))
    for token in accepted:
        bare = token.strip("'")
        assert bare in header, f"{bare!r} is accepted but missing from the signature line"


def test_bottom_holds_the_model_top_and_top_can_move_it():
    """The behavioural claim the docstring makes, pinned with real numbers."""

    import pandas as pd

    from myflopy.modflow.mf6.grid.geometry import reconcile_surfaces

    crossed = pd.DataFrame({"top": [100.0], "a": [105.0], "b": [90.0]})
    lowered = reconcile_surfaces(None, df=crossed.copy(), which="bottom")
    raised = reconcile_surfaces(None, df=crossed.copy(), which="top")

    assert lowered["top"][0] == 100.0, "'bottom' must leave the model top alone"
    assert lowered["a"][0] < lowered["top"][0]
    assert raised["top"][0] > 100.0, "'top' moves the model top -- documented as such"


# --- erosion: a surface cutting down through the layers below ---------------- #

def _channel_stack(pinch, cap=False):
    """Ground at 200 with a channel incised to 115, cutting sand AND clay."""

    ground_vals = np.array([200.0, 115.0, 200.0])          # middle cell is the channel
    ground = Surface.from_array(ground_vals)
    sand_base = Flat(150).capped_at(ground - 2.0) if cap else Flat(150)
    return (
        LayerStack(top=ground, length_units="feet")
        .add("sand", bottom=sand_base, min_thickness=2.0, pinch=pinch)
        .add("clay", bottom=Flat(120), min_thickness=2.0, pinch=pinch)
        .add("till", bottom=Flat(60), min_thickness=2.0, pinch=pinch)
    )


@pytest.mark.parametrize(
    ("pinch", "cut_idomain"), [("inactive", 0), ("passthrough", -1)]
)
def test_an_incised_channel_cuts_every_layer_it_passes_through(pinch, cut_idomain):
    """The documented erosion mechanism: reconcile pushes the cut contacts down
    and `pinch` decides what the emptied cells become. No separate "cut" call."""

    layers = _channel_stack(pinch).build(_fake_vor(3))
    ch = 1                                                  # the channel cell

    assert layers.idomain[0][ch] == cut_idomain, "sand is cut out"
    assert layers.idomain[1][ch] == cut_idomain, "clay is cut out too -- it cascades"
    assert layers.idomain[2][ch] == 1, "till survives and floors the channel"
    assert layers.thickness[2][ch] > 50, "till reaches up to the channel floor"


def test_capping_the_contact_keeps_a_veneer_instead_of_cutting_it_out():
    """The explicit alternative. Same channel, different geology: `capped_at`
    thins the unit to a remnant that stays ACTIVE, where reconcile removes it."""

    vor = _fake_vor(3)
    ch = 1
    cut = _channel_stack("inactive").build(vor)
    capped = _channel_stack("inactive", cap=True).build(vor)

    assert cut.idomain[0][ch] == 0 and cut.thickness[0][ch] < 1
    assert capped.idomain[0][ch] == 1
    assert capped.thickness[0][ch] == pytest.approx(2.0)


# --- looking at a stack before its grid exists (2026-08-27) ------------------ #
#
# `.plot` on a gridless stack draws on a coarse `draft_grid` over the surfaces'
# own extent. Pictures only: an approximate picture is useful, invented model
# geometry is not, so build/qc/to_disv still demand a real grid.

@pytest.fixture
def dem(tmp_path):
    """A small georeferenced GeoTIFF: 3000 x 2000 at EPSG:2927."""

    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path = tmp_path / "ground.tif"
    values = np.array(
        [[200.0 - 0.02 * (x * 25) for x in range(120)] for _ in range(80)],
        dtype="float32",
    )
    with rasterio.open(
        path, "w", driver="GTiff", height=80, width=120, count=1, dtype="float32",
        crs="EPSG:2927", transform=from_origin(0, 2000, 25, 25),
    ) as handle:
        handle.write(values, 1)
    return path


def _gridless(dem):
    return (
        LayerStack(top=Raster(dem), length_units="feet")
        .add("sand", thickness=40.0, pinch="inactive")
        .add("clay", bottom=Flat(120), min_thickness=2.0, pinch="passthrough")
    )


@pytest.mark.parametrize("verb", ["map", "section", "surface"])
def test_a_gridless_stack_can_still_be_looked_at(dem, verb):
    stack = _gridless(dem)
    assert stack.vor is None
    assert getattr(stack.plot, verb)() is not None


def test_the_draft_grid_takes_its_extent_from_the_surfaces(dem):
    grid = _gridless(dem).draft_grid(cells=200)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    assert (round(xmin), round(ymin), round(xmax), round(ymax)) == (0, 0, 3000, 2000)
    assert "2927" in str(grid.crs)


def test_the_draft_grid_is_built_once_and_reused(dem, monkeypatch):
    """Otherwise every picture pays for a fresh mesh."""

    stack = _gridless(dem)
    calls = []
    original = LayerStack.draft_grid
    monkeypatch.setattr(
        LayerStack, "draft_grid",
        lambda self, **kw: (calls.append(1), original(self, **kw))[1],
    )
    stack.plot.map()
    stack.plot.section(y=1000)
    assert len(calls) == 1


@pytest.mark.parametrize("consume", ["build", "qc", "to_disv"])
def test_the_draft_grid_is_never_used_for_model_geometry(dem, consume):
    """The line that makes the fallback safe: pictures may be approximate,
    a DISV may not be invented."""

    with pytest.raises(ValueError, match="needs a grid"):
        getattr(_gridless(dem), consume)()


def test_an_explicit_extent_overrides_the_surfaces():
    """Also the escape hatch for a stack whose surfaces carry no extent."""

    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    grid = stack.draft_grid(cells=100, extent=(0, 0, 500, 400), crs="EPSG:2927")
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    assert (round(xmax - xmin), round(ymax - ymin)) == (500, 400)


def test_a_stack_with_no_extent_anywhere_says_so():
    """`Flat`/`Array` describe thickness with no notion of where -- nothing to infer."""

    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    with pytest.raises(ValueError, match="no surface in this stack carries"):
        stack.draft_grid()


def test_an_empty_extent_is_rejected():
    stack = LayerStack(top=Flat(200)).add("sand", thickness=40.0)
    with pytest.raises(ValueError, match="empty extent"):
        stack.draft_grid(extent=(10, 10, 10, 50))


def test_the_extent_walk_descends_into_composed_surfaces(dem):
    """A stack whose only georeferenced surface is buried inside algebra --
    `Flat(150).capped_at(Raster(...) - 2)` -- must still find the raster."""

    ground = Raster(dem)
    stack = LayerStack(top=Flat(300)).add(
        "sand", bottom=Flat(150).capped_at(ground - 2.0)
    )
    assert stack._georeferenced_surfaces(), "the walk missed a nested raster"
    assert round(stack.draft_grid(cells=100).gdf_vorPolys.total_bounds[2]) == 3000


# --- contact surfaces as individual VTK sheets (2026-08-27) ----------------- #

@pytest.mark.slow
def test_surface_renders_each_contact_as_its_own_vtk_sheet(real_vor):
    """`grid(backend="vtk")` fuses the layers into one VOLUME, so you cannot look
    at a contact on its own. This renders one actor per surface, which is what
    lets a viewer hide the sheets above and see underneath."""

    pytest.importorskip("pyvista")
    from myflopy.viz import VtkScene

    stack = (LayerStack(real_vor, top=Flat(100))
             .add("sand", thickness=20).add("clay", thickness=30))
    scene = stack.plot.surface("all", backend="vtk", resolution=40)

    assert isinstance(scene, VtkScene)
    meshes = [a for a in scene.scene.renderer.actors.values()
              if a.GetMapper() is not None and a.GetMapper().GetInput() is not None
              and a.GetMapper().GetInput().GetNumberOfCells() > 1]
    assert len(meshes) == 3, "one sheet per surface: top + two bottoms"


@pytest.mark.slow
def test_a_single_named_contact_can_be_drawn_alone(real_vor):
    pytest.importorskip("pyvista")

    stack = (LayerStack(real_vor, top=Flat(100))
             .add("sand", thickness=20).add("clay", thickness=30))
    scene = stack.plot.surface("clay", backend="vtk", resolution=30)
    meshes = [a for a in scene.scene.renderer.actors.values()
              if a.GetMapper() is not None and a.GetMapper().GetInput() is not None
              and a.GetMapper().GetInput().GetNumberOfCells() > 1]
    assert len(meshes) == 1


def test_surface_rejects_an_unknown_backend(real_vor):
    stack = LayerStack(real_vor, top=Flat(100)).add("sand", thickness=20)
    with pytest.raises(ValueError, match="backend must be 'plotly' or 'vtk'"):
        stack.plot.surface(backend="opengl")


@pytest.mark.slow
def test_no_data_cells_do_not_become_a_sheet_at_zero(tmp_path):
    """A NaN cell rendered as z=0 reads as a real contact at sea level, which is
    exactly the failure `nodata=` was added to stop. Those cells are dropped.

    Builds its own denser grid: `real_vor` has four cells, too few to interpolate
    from once any are missing.
    """

    pytest.importorskip("pyvista")
    import myflopy as mf

    tri = mf.TriangleGrid(model_ws=str(tmp_path), angle=30)
    tri.set_domain_rectangle(x_dist=400, y_dist=300, origin=(0, 0), max_area=800)
    tri.build()
    vor = mf.VoronoiGridPlus(tri, crs="EPSG:2927")

    zs = np.full(vor.ncpl, 100.0)
    zs[: vor.ncpl // 3] = np.nan            # a third of the domain has no data
    stack = LayerStack(vor, top=Surface.from_array(zs)).add("a", thickness=10)
    scene = stack.plot.surface("top", backend="vtk", resolution=40)

    sheet = next(a.GetMapper().GetInput() for a in scene.scene.renderer.actors.values()
                 if a.GetMapper() is not None and a.GetMapper().GetInput() is not None
                 and a.GetMapper().GetInput().GetNumberOfCells() > 1)
    bounds = sheet.GetBounds()          # (xmin, xmax, ymin, ymax, zmin, zmax)
    assert bounds[4] > 0, "a sheet reaching z=0 means the NaN cells were drawn"
