"""A drawn surface is masked to its domain, not to the convex hull of its points.

Ledger 314. `griddata` fills the convex hull of the points it is handed and is NaN
outside it, so a surface left to it LOOKS clipped while spanning every concavity of a
domain that is not convex. `.show()` drew exactly that -- `.fig` was
`clipped_fig(clip=False)` -- and only `plot_heatmap` masked, through a raster whose
transform sat half a pixel off, filling the clipped-away pixels with 0 and then
blanking every real 0 with them. The heatmap itself was drawn upside down.

Every test here uses an L-shaped grid of square cells, its north-east block removed.
On a CONVEX domain the hull of the cell centres sits inside the footprint, so a surface
with no mask at all is indistinguishable from a masked one, and a test on the usual
rectangular fixture would pass with the whole fix reverted.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import plotly.graph_objects as go
import pytest
import rasterio
import rasterio.transform
import shapely
from shapely.geometry import MultiPolygon, Polygon, box

from myflopy import viz
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

#: The grid below, drawn independently of it: 100 x 100 with the NE 60 x 60 cut out.
L_SHAPE = Polygon([(0, 0), (100, 0), (100, 40), (40, 40), (40, 100), (0, 100)])
#: The removed block's interior, well clear of its edges.
NOTCH = box(50.0, 50.0, 95.0, 95.0)
#: Pixel centres run 5..95 (the cell centres' extent). At 60 per axis none lands
#: exactly on x = 40 or y = 40, where "inside" would turn on a tie-break rule.
RESOLUTION = 60


def _l_grid(n: int = 10, cell: float = 10.0, notch: int = 4) -> VoronoiGridPlus:
    """``n`` x ``n`` square cells with the block ``i, j >= notch`` removed -- an L."""

    keep = [(i, j) for j in range(n) for i in range(n) if not (i >= notch and j >= notch)]
    index: dict[tuple[int, int], int] = {}
    verts: list[tuple[float, float]] = []

    def vertex(i: int, j: int) -> int:
        if (i, j) not in index:
            index[(i, j)] = len(verts)
            verts.append((i * cell, j * cell))
        return index[(i, j)]

    iverts = [[vertex(i, j), vertex(i, j + 1), vertex(i + 1, j + 1), vertex(i + 1, j)] for i, j in keep]
    xcyc = np.array([((i + 0.5) * cell, (j + 0.5) * cell) for i, j in keep])
    grid = VoronoiGridPlus(verts=np.array(verts, dtype=float), iverts=iverts, xcyc=xcyc)
    polys = grid.gdf_vorPolys
    # A top that rises to the north, so orientation is readable off the values.
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: 100.0 + 0.1 * xcyc[:, 1], 1: np.full(len(keep), 50.0)},
        geometry=polys.geometry,
        crs=polys.crs,
    )
    return grid


@pytest.fixture(scope="module")
def grid():
    return _l_grid()


def _points(grid):
    xs = np.asarray(grid.centroids_x, dtype=float)
    ys = np.asarray(grid.centroids_y, dtype=float)
    return xs, ys, 100.0 + 0.1 * xs + 0.2 * ys


def _surface(grid, **kwargs) -> InterpolatedSurface:
    xs, ys, zs = _points(grid)
    return InterpolatedSurface(xs, ys, zs, grid, resolution=RESOLUTION, **kwargs)


def _in(region, surface: InterpolatedSurface) -> np.ndarray:
    grid_x, grid_y = surface.xy_meshgrid
    return shapely.contains_xy(region, grid_x, grid_y)


def _drawn(figure) -> np.ndarray:
    return np.asarray(figure.data[0].z, dtype=float)


def test_the_fixture_is_one_a_missing_mask_would_show_on(grid):
    """The precondition: unmasked, the interpolation reaches into the notch."""

    surface = _surface(grid)
    bare = np.asarray(surface.surface, dtype=float)

    assert (np.isfinite(bare) & _in(NOTCH, surface)).sum() > 100


@pytest.mark.parametrize("how", ["points supplied with the grid", "points read from the grid"])
def test_show_draws_the_domain_and_nothing_outside_it(grid, how):
    """The reported call: `InterpolatedSurface(xs, ys, zs, vor).show()`.

    Asserted in BOTH directions -- nothing outside the L, and every value the
    interpolator had inside it kept -- so masking too much fails as well.
    """

    if how == "points supplied with the grid":
        surface = _surface(grid)
    else:
        surface = InterpolatedSurface(vor=grid, surf_type="lyr", resolution=RESOLUTION)
    drawn = _drawn(surface.fig)
    bare = np.asarray(surface.surface, dtype=float)
    inside = _in(L_SHAPE, surface)

    assert not (np.isfinite(drawn) & ~inside).any()
    np.testing.assert_array_equal(np.isfinite(drawn), np.isfinite(bare) & inside)
    np.testing.assert_array_equal(drawn[inside], bare[inside])


def test_masking_does_not_write_into_the_cached_interpolation(grid):
    surface = _surface(grid)
    before = np.array(surface.surface, dtype=float)
    _ = surface.fig

    np.testing.assert_array_equal(np.asarray(surface.surface, dtype=float), before)
    np.testing.assert_array_equal(_drawn(surface.clipped_fig(clip=False)), before)


def test_surface_trace_draws_what_fig_draws(grid):
    surface = _surface(grid)

    np.testing.assert_array_equal(
        np.asarray(surface.surface_trace().z, dtype=float), _drawn(surface.fig)
    )


def test_clip_masks_the_drawn_surface_even_for_points_supplied_directly(grid):
    """`clip=` did nothing to `.show()`: supplied points are never subset, and the
    figure was never masked, so the only thing that read the clip was the heatmap."""

    region = box(0.0, 0.0, 70.0, 30.0)
    surface = _surface(grid, clip=region)
    drawn = _drawn(surface.fig)

    assert np.isfinite(drawn).any()
    assert not (np.isfinite(drawn) & ~_in(region, surface)).any()


@pytest.mark.parametrize("spelling", ["str", "MultiPolygon"])
def test_a_clip_the_setter_used_to_drop_is_used(grid, tmp_path, spelling):
    """A `str` path and a `MultiPolygon` fell through a setter with no branch for
    them, leaving `clip` None -- an unclipped surface, and nothing said."""

    parts = MultiPolygon([box(0.0, 0.0, 30.0, 30.0), box(60.0, 0.0, 90.0, 30.0)])
    if spelling == "str":
        path = tmp_path / "clip.gpkg"
        gpd.GeoDataFrame(geometry=[box(0.0, 0.0, 30.0, 30.0)], crs=grid.crs).to_file(path)
        clip, region = str(path), box(0.0, 0.0, 30.0, 30.0)
    else:
        clip, region = parts, parts
    surface = _surface(grid, clip=clip)
    drawn = _drawn(surface.fig)

    assert surface.clip is not None
    assert np.isfinite(drawn).any()
    assert not (np.isfinite(drawn) & ~_in(region, surface)).any()


def test_a_clip_of_another_type_is_refused(grid):
    with pytest.raises(TypeError, match="clip must be a Polygon"):
        _surface(grid, clip=gpd.GeoSeries([box(0.0, 0.0, 30.0, 30.0)]))


def test_supplied_values_do_not_ask_which_cells_the_clip_covers(grid):
    """`zs` resolved the clip's cells before looking at the values it was given, so a
    clip away from the GRID refused a point cloud that never came from it."""

    xs, ys, zs = _points(grid)
    surface = InterpolatedSurface(xs, ys, zs, grid, clip=box(500.0, 500.0, 600.0, 600.0))

    np.testing.assert_array_equal(surface.zs, zs)


def test_a_mask_that_leaves_nothing_to_draw_says_where_both_are(grid):
    surface = _surface(grid, clip=box(500.0, 500.0, 600.0, 600.0))

    with pytest.raises(ValueError, match="leaves nothing to draw") as raised:
        _ = surface.fig
    assert "(500.0, 500.0, 600.0, 600.0)" in str(raised.value)


def test_a_bare_point_cloud_is_drawn_as_interpolated_and_refuses_a_forced_clip(grid):
    xs, ys, zs = _points(grid)
    surface = InterpolatedSurface(xs, ys, zs, resolution=RESOLUTION)

    assert surface.region is None
    np.testing.assert_array_equal(_drawn(surface.fig), np.asarray(surface.surface, dtype=float))
    with pytest.raises(ValueError, match="nothing to clip to"):
        surface.clipped_fig(clip=True)
    with pytest.raises(ValueError, match="nothing to clip to"):
        surface.clip_raster_with_polygon()


def test_clipped_fig_masks_a_surface_it_is_handed(grid):
    """`clipped_fig(surface=..., clip=True)` threw the caller's surface away."""

    surface = _surface(grid)
    other = np.full((RESOLUTION, RESOLUTION), 7.0)
    drawn = _drawn(surface.clipped_fig(surface=other, clip=True))
    inside = _in(L_SHAPE, surface)

    assert (drawn[inside] == 7.0).all()
    assert np.isnan(drawn[~inside]).all()


def test_the_raster_is_registered_to_the_pixels_the_surface_draws(grid):
    """`from_origin` takes a CORNER; the transform handed it the first CENTRE."""

    surface = _surface(grid)
    grid_x, grid_y = surface.xy_meshgrid
    last = RESOLUTION - 1
    for row, col in [(0, 0), (0, last), (last, 0), (last, last)]:
        x, y = rasterio.transform.xy(surface.transform, row, col)
        assert (x, y) == pytest.approx((grid_x[row, col], grid_y[row, col]), abs=1e-9)


def test_the_raster_clip_keeps_exactly_the_pixels_inside_the_domain(grid):
    """Half a pixel off, the clip kept a row of the notch and dropped a row of the L."""

    surface = _surface(grid)
    clipped = surface.clip_raster_with_polygon()[0][0]
    bare = np.asarray(surface.surface, dtype=float)

    np.testing.assert_array_equal(np.isfinite(clipped), np.isfinite(bare) & _in(L_SHAPE, surface))
    np.testing.assert_array_equal(np.isfinite(clipped), np.isfinite(_drawn(surface.fig)))


def test_a_real_zero_survives_every_clip(grid, tmp_path):
    """The clip filled with 0 and callers then turned every 0 into NaN -- a sea-level
    datum or an undrawn-down cell included. Linear interpolation reproduces an
    exactly-zero region exactly, which is what makes the loss visible."""

    xs, ys, _ = _points(grid)
    zs = np.where(xs < 50.0, 0.0, xs - 50.0)
    surface = InterpolatedSurface(
        xs, ys, zs, grid, resolution=RESOLUTION, interpolator="linearND"
    )
    inside = _in(L_SHAPE, surface)
    zeros = (np.asarray(surface.surface, dtype=float) == 0.0) & inside
    assert zeros.sum() > 100

    output = tmp_path / "clipped.tif"
    surface.save_raster(output_tif=output)
    with rasterio.open(output) as saved:
        written = saved.read(1)
        assert np.isnan(saved.nodata)

    for label, array in [
        (".fig", _drawn(surface.fig)),
        ("clipped_fig(clip=True)", _drawn(surface.clipped_fig(clip=True))),
        ("clip_raster_with_polygon", surface.clip_raster_with_polygon()[0][0]),
        ("save_raster", written),
    ]:
        assert (array[zeros] == 0.0).all(), label
        assert np.isnan(array[~inside]).all(), label


def test_save_raster_takes_the_pieces_clip_raster_with_polygon_returns(grid, tmp_path):
    """`any([ndarray, ...])` asked an array for its truth value, so passing them raised."""

    surface = _surface(grid)
    image, transform, meta = surface.clip_raster_with_polygon()
    output = tmp_path / "given.tif"

    surface.save_raster(output_tif=output, clipped_image=image, clipped_transform=transform, meta=meta)
    with rasterio.open(output) as saved:
        np.testing.assert_array_equal(saved.read(1), image[0])
    with pytest.raises(ValueError, match="all three"):
        surface.save_raster(output_tif=tmp_path / "part.tif", clipped_image=image)


def test_plot_heatmap_is_north_up_on_the_model_coordinates(grid, monkeypatch):
    """Row 0 of the raster is its NORTH edge, and a heatmap given no `y` draws row 0
    at the BOTTOM: the L's north-east notch appeared in the south-east."""

    shown = []
    monkeypatch.setattr(viz.Fig, "show", lambda self, *args, **kwargs: shown.append(self))
    surface = InterpolatedSurface(vor=grid, surf_type="lyr", resolution=RESOLUTION)

    assert surface.plot_heatmap() is None
    (figure,) = shown
    trace = figure.data[0]
    x = np.asarray(trace.x, dtype=float)
    y = np.asarray(trace.y, dtype=float)
    z = np.asarray(trace.z, dtype=float)

    assert (np.diff(y) > 0).all(), "rows must be handed to plotly south to north"
    assert np.nanmean(z[-1]) > np.nanmean(z[0]), "the top rises to the north"
    north_east = (y[:, None] > 50.0) & (x[None, :] > 50.0)
    south_east = (y[:, None] < 35.0) & (x[None, :] > 50.0)
    assert np.isnan(z[north_east]).all()
    assert np.isfinite(z[south_east]).all()
    assert figure.layout.yaxis.scaleanchor == "x"


def test_the_stack_surface_is_masked_to_the_grid(grid):
    """`stack.plot.surface` built its surfaces with no grid, so nothing could mask them."""

    from myflopy.layers import Flat, LayerStack

    stack = LayerStack(grid, top=Flat(100.0)).add("sand", thickness=20.0).build()
    figure = stack.plot.surface("all", resolution=RESOLUTION).fig
    probe = InterpolatedSurface(*_points(grid), resolution=RESOLUTION)
    inside = _in(L_SHAPE, probe)

    traces = [trace for trace in figure.data if isinstance(trace, go.Surface)]
    assert len(traces) == 2
    for trace in traces:
        drawn = np.asarray(trace.z, dtype=float)
        assert np.isfinite(drawn).any()
        assert not (np.isfinite(drawn) & ~inside).any(), trace.name


@pytest.mark.slow
def test_the_stack_vtk_sheets_are_masked_like_the_plotly_ones(grid):
    """The VTK sheets are documented as the same numbers as their Plotly twins."""

    pytest.importorskip("pyvista")
    from myflopy.layers import Flat, LayerStack

    stack = LayerStack(grid, top=Flat(100.0)).add("sand", thickness=20.0).build()
    scene = stack.plot.surface("all", backend="vtk", resolution=RESOLUTION)
    # A kept cell may have corners up to one pixel outside the edge it straddles.
    pixel = 90.0 / (RESOLUTION - 1)
    near = L_SHAPE.buffer(pixel * 1.5)

    assert len(scene.meshes) == 2
    for sheet in scene.meshes:
        xy = np.asarray(sheet.points)[:, :2]
        assert len(xy) > 0
        assert shapely.contains_xy(near, xy[:, 0], xy[:, 1]).all()
