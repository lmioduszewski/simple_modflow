from __future__ import annotations


def test_choro_plot_returns_figure_without_calling_show():
    import plotly.graph_objects as go

    from myflopy.modflow.utils.datatypes.choros import Choro

    choro = object.__new__(Choro)
    figure = go.Figure()
    choro._fig = figure
    choro.add_choropleth = lambda: None
    choro.add_contours = lambda: None
    choro._locs = None
    choro._overlays = []
    choro.hillshade_path = None
    # `.fig` also draws the `select=` highlight (ledger 160). Stubbed rather than
    # defended against with getattr: a half-built Choro should fail loudly here,
    # not silently draw a map with its highlight missing.
    choro._select_cells = []
    choro.select_style = "outline"

    assert choro.fig is figure

import os
import shutil
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from myflopy.modflow.mf6.grid.plotting import _choropleth_factory  # noqa: E402
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from myflopy.modflow.mf6.interactive_plotting import (  # noqa: E402
    FrameExportProgress,
    ModelMapStyle,
    StandaloneHtmlSlider,
    _select_plotly_frames,
    build_particle_tracking_scene,
    export_cross_section_slider_html,
    export_head_layer_mosaic_slider_html,
    export_head_map_slider_html,
    export_matplotlib_slider_html,
    plot_model_head_map,
)
from myflopy.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from myflopy.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.packages import KFlow  # noqa: E402
from myflopy.plot import ModelPlots  # noqa: E402
from myflopy.viz import (
    VtkScene,  # noqa: E402
    _write_plotly_choropleth_restyle_html,  # noqa: E402
)


def _workspace(name: str) -> Path:
    path = ROOT / ".pytest-work" / f"{name}_{time.time_ns()}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _two_cell_model(name: str, workspace: Path, *, nlay: int = 1, nper: int = 2):
    verts = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        dtype=float,
    )
    vor = VoronoiGridPlus(
        verts=verts,
        iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float),
    )
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=nper)
    bottom = [[2.0, 2.0] for _ in range(nlay)]
    DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=bottom, nlay=nlay, idomain=[[1, 1]] * nlay)
    TemporalDiscretization(model=model, per_len=1, num_steps=1, multiplier=1.0)
    KFlow(model=model, k=np.ones((nlay, 2)), k33_vert=np.full((nlay, 2), 0.25))
    return model


def _large_multilayer_model(name: str, workspace: Path, *, nlay: int = 8, nper: int = 12):
    nrow, ncol = 25, 40
    vertices = [(float(col), float(row)) for row in range(nrow + 1) for col in range(ncol + 1)]
    iverts = []
    centers = []
    for row in range(nrow):
        for col in range(ncol):
            lower_left = row * (ncol + 1) + col
            lower_right = lower_left + 1
            upper_left = (row + 1) * (ncol + 1) + col
            upper_right = upper_left + 1
            iverts.append([lower_left, upper_left, upper_right, lower_right])
            centers.append((col + 0.5, row + 0.5))
    vor = VoronoiGridPlus(
        verts=np.asarray(vertices, dtype=float),
        iverts=iverts,
        xcyc=np.asarray(centers, dtype=float),
    )
    ncpl = nrow * ncol
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=nper)
    DisvGrid(
        vor=vor,
        model=model,
        top=np.full(ncpl, 100.0),
        bottom=[np.full(ncpl, 100.0 - 10.0 * (layer + 1)) for layer in range(nlay)],
        nlay=nlay,
        idomain=np.ones((nlay, ncpl), dtype=int),
    )
    TemporalDiscretization(model=model, per_len=30.0, num_steps=2, multiplier=1.0)
    KFlow(model=model, k=np.ones((nlay, ncpl)), k33_vert=np.full((nlay, ncpl), 0.25))
    return model


def _large_head_frames(*, nlay: int = 8, ncpl: int = 1_000, nper: int = 12):
    cell_gradient = np.linspace(0.0, 12.0, ncpl)
    frames = []
    for period in range(nper):
        frame = np.vstack(
            [95.0 - 7.0 * layer - 0.5 * period - cell_gradient for layer in range(nlay)]
        )
        frames.append(frame)
    frames[2][1, 10] = 1.0e30
    frames[4][3, 75] = np.nan
    return frames


def _assert_standalone_slider(path: Path, *, expected_frames: int):
    text = path.read_text(encoding="utf-8")
    assert "<!doctype html>" in text
    assert 'type="range"' in text
    assert "data:image/png;base64," in text
    assert "setInterval" in text
    assert text.count("data:image/png;base64,") == expected_frames
    assert "ipywidgets" not in text


def test_generic_matplotlib_slider_is_self_contained(tmp_path):
    output = tmp_path / "generic_slider.html"

    def render(value, _index):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [value, value + 1])
        ax.set_title(f"Frame {value}")
        return fig

    result = export_matplotlib_slider_html(
        render,
        [0, 1, 2],
        output,
        labels=["start", "middle", "end"],
        title="Training-style frame slider",
    )

    assert isinstance(result, StandaloneHtmlSlider)
    assert result.frame_count == 3
    assert result.labels == ("start", "middle", "end")
    _assert_standalone_slider(output, expected_frames=3)


def test_external_slider_supports_selection_progress_and_resume(tmp_path):
    output = tmp_path / "resumable.html"
    rendered_values = []
    events: list[FrameExportProgress] = []

    def render(value, _index):
        rendered_values.append(value)
        fig, ax = plt.subplots()
        ax.plot([0, 1], [value, value + 1])
        return fig

    first = export_matplotlib_slider_html(
        render,
        list(range(10)),
        output,
        embed_frames=False,
        frame_stride=2,
        max_frames=3,
        progress=events.append,
    )

    assert first.labels == ("0", "2", "4")
    assert first.rendered_frames == 3
    assert first.reused_frames == 0
    assert rendered_values == [0, 2, 4]
    assert [event.status for event in events] == ["rendered"] * 3
    assert (first.frame_directory / "frames.json").exists()

    (first.frame_directory / "frame_00001.png").unlink()
    rendered_values.clear()
    events.clear()
    resumed = export_matplotlib_slider_html(
        render,
        list(range(10)),
        output,
        embed_frames=False,
        frame_stride=2,
        max_frames=3,
        resume=True,
        progress=events.append,
    )

    assert rendered_values == [2]
    assert resumed.rendered_frames == 1
    assert resumed.reused_frames == 2
    assert [event.status for event in events] == ["reused", "rendered", "reused"]


def test_external_slider_rejects_manifest_mismatch_and_concurrent_writer(tmp_path):
    output = tmp_path / "locked.html"

    def render(value, _index):
        fig, ax = plt.subplots()
        ax.plot([0, 1], [value, value])
        return fig

    export_matplotlib_slider_html(render, [0, 1], output, embed_frames=False)
    with pytest.raises(ValueError, match="manifest"):
        export_matplotlib_slider_html(
            render,
            [0, 1, 2],
            output,
            embed_frames=False,
            resume=True,
        )

    lock_path = output.with_suffix(".html.lock")
    lock_path.write_text("another process", encoding="utf-8")
    try:
        with pytest.raises(RuntimeError, match="Another export"):
            export_matplotlib_slider_html(render, [0], output)
    finally:
        lock_path.unlink()


def test_slider_selection_validation(tmp_path):
    with pytest.raises(ValueError, match="frame_stride"):
        export_matplotlib_slider_html(lambda *_: plt.figure(), [1], tmp_path / "stride.html", frame_stride=0)
    with pytest.raises(IndexError, match="frame_indices"):
        export_matplotlib_slider_html(
            lambda *_: plt.figure(),
            [1],
            tmp_path / "indices.html",
            frame_indices=[1],
        )
    with pytest.raises(ValueError, match="embed_frames=False"):
        export_matplotlib_slider_html(
            lambda *_: plt.figure(),
            [1],
            tmp_path / "resume.html",
            resume=True,
        )


def test_head_map_renderer_uses_flopy_plot_map_view_style():
    workspace = _workspace("head_map_renderer")
    try:
        model = _two_cell_model("head_map_renderer", workspace)
        fig, ax = plot_model_head_map(
            model,
            np.array([[8.0, 9.0]]),
            style=ModelMapStyle(show_contours=False),
            title="Heads",
        )
        assert ax.get_title() == "Heads"
        assert len(ax.collections) >= 1
        plt.close(fig)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_head_map_renderer_normalizes_real_multilayer_disv_shape():
    workspace = _workspace("head_map_disv_shape")
    try:
        model = _two_cell_model("head_map_disv_shape", workspace, nlay=2)
        real_disv_shape = np.array([[[8.0, 9.0]], [[6.0, 7.0]]])
        fig, ax = plot_model_head_map(
            model,
            real_disv_shape,
            layer=1,
            style=ModelMapStyle(show_contours=False),
        )
        assert len(ax.collections) >= 1
        plt.close(fig)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_head_map_and_layer_mosaic_sliders_match_training_notebook_patterns(tmp_path):
    workspace = _workspace("head_map_sliders")
    try:
        model = _two_cell_model("head_map_sliders", workspace, nlay=2)
        frames = [
            np.array([[8.0, 9.0], [6.0, 7.0]]),
            np.array([[8.5, 9.5], [6.5, 7.5]]),
        ]
        style = ModelMapStyle(show_contours=False, dpi=70)

        map_result = export_head_map_slider_html(
            model,
            tmp_path / "head_map.html",
            head_frames=frames,
            labels=["period 1", "period 2"],
            layer=0,
            style=style,
        )
        mosaic_result = export_head_layer_mosaic_slider_html(
            model,
            tmp_path / "mosaic.html",
            head_frames=frames,
            labels=["period 1", "period 2"],
            layers=[0, 1],
            ncols=2,
            style=style,
        )

        assert map_result.frame_count == 2
        assert mosaic_result.frame_count == 2
        _assert_standalone_slider(map_result.path, expected_frames=2)
        _assert_standalone_slider(mosaic_result.path, expected_frames=2)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_head_layer_mosaic_has_high_resolution_default_and_direct_controls(monkeypatch, tmp_path):
    workspace = _workspace("head_mosaic_resolution")
    try:
        model = _two_cell_model("head_mosaic_resolution", workspace, nlay=2)
        captured = {}
        original = export_matplotlib_slider_html

        def capture_export(render_frame, values, output_path, **kwargs):
            captured["dpi"] = kwargs["dpi"]
            fig = render_frame(values[0], 0)
            captured["figsize"] = tuple(fig.get_size_inches())
            plt.close(fig)
            return original(render_frame, values[:1], output_path, **kwargs)

        def draw_without_contours(_model, _data, *, ax, title, **_kwargs):
            ax.set_title(title)
            ax.scatter([0, 1], [0, 1])
            return ax.figure, ax

        monkeypatch.setattr(
            "myflopy.modflow.mf6.interactive_plotting.export_matplotlib_slider_html",
            capture_export,
        )
        monkeypatch.setattr(
            "myflopy.modflow.mf6.interactive_plotting.plot_model_head_map",
            draw_without_contours,
        )
        export_head_layer_mosaic_slider_html(
            model,
            tmp_path / "default_mosaic.html",
            head_frames=[np.array([[8.0, 9.0], [6.0, 7.0]])],
            layers=[0, 1],
            ncols=2,
        )
        assert captured["dpi"] == 200

        export_head_layer_mosaic_slider_html(
            model,
            tmp_path / "controlled_mosaic.html",
            head_frames=[np.array([[8.0, 9.0], [6.0, 7.0]])],
            layers=[0, 1],
            ncols=2,
            dpi=275,
            panel_figsize=(6, 4),
        )
        assert captured["dpi"] == 275
        assert captured["figsize"] == pytest.approx((12, 4))

        with pytest.raises(ValueError, match="dpi"):
            export_head_layer_mosaic_slider_html(
                model,
                tmp_path / "bad_dpi.html",
                head_frames=[np.array([[8.0, 9.0], [6.0, 7.0]])],
                dpi=0,
            )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_large_multilayer_results_use_external_frames_and_mask_dry_values(tmp_path):
    workspace = _workspace("large_multilayer_slider")
    try:
        # 4 layers x 6 periods keeps every exercised behavior (external
        # frames, the 1e30 dry mask in frame 2, the NaN in frame 4 on layer 3)
        # while rendering a fraction of the PNG panels.
        nlay, nper = 4, 6
        model = _large_multilayer_model(
            "large_multilayer_slider", workspace, nlay=nlay, nper=nper
        )
        frames = _large_head_frames(nlay=nlay, nper=nper)
        style = ModelMapStyle(show_contours=True, dpi=45, figsize=(5, 4))

        map_result = export_head_map_slider_html(
            model,
            tmp_path / "large_map.html",
            head_frames=frames,
            layer=1,
            style=style,
            embed_frames=False,
        )
        mosaic_result = export_head_layer_mosaic_slider_html(
            model,
            tmp_path / "large_mosaic.html",
            head_frames=frames,
            layers=list(range(nlay)),
            ncols=2,
            style=style,
            embed_frames=False,
        )
        cross_section_result = export_cross_section_slider_html(
            model,
            LineString([(0.5, 0.5), (39.5, 24.5)]),
            tmp_path / "large_cross_section.html",
            head_frames=frames,
            labels=[f"period {period + 1}" for period in range(nper)],
            dpi=45,
            embed_frames=False,
            show_legend=False,
        )

        for result in (map_result, mosaic_result, cross_section_result):
            assert result.frame_count == nper
            assert not result.embedded_frames
            assert result.frame_directory is not None
            assert len(list(result.frame_directory.glob("*.png"))) == nper
            html = result.path.read_text(encoding="utf-8")
            assert "data:image/png;base64," not in html
            assert f"frame_{nper - 1:05d}.png" in html  # last frame referenced
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_large_multilayer_shared_color_scale_ignores_dry_and_nan_values(monkeypatch, tmp_path):
    workspace = _workspace("large_multilayer_scale")
    try:
        model = _large_multilayer_model("large_multilayer_scale", workspace)
        frames = _large_head_frames()
        limits = []
        original = plot_model_head_map

        def capture_limits(*args, **kwargs):
            limits.append((kwargs["vmin"], kwargs["vmax"]))
            return original(*args, **kwargs)

        monkeypatch.setattr(
            "myflopy.modflow.mf6.interactive_plotting.plot_model_head_map",
            capture_limits,
        )
        export_head_map_slider_html(
            model,
            tmp_path / "shared_scale.html",
            head_frames=frames,
            layer=1,
            style=ModelMapStyle(show_contours=False, show_colorbar=False, dpi=35, figsize=(4, 3)),
        )

        expected_min = min(np.nanmin(np.where(np.abs(frame[1]) < 1.0e29, frame[1], np.nan)) for frame in frames)
        expected_max = max(np.nanmax(np.where(np.abs(frame[1]) < 1.0e29, frame[1], np.nan)) for frame in frames)
        assert limits == [(expected_min, expected_max)] * 12
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_cross_section_slider_reuses_canonical_matplotlib_cross_section(tmp_path):
    workspace = _workspace("cross_section_slider")
    try:
        model = _two_cell_model("cross_section_slider", workspace)
        result = export_cross_section_slider_html(
            model,
            LineString([(0.5, -0.5), (0.5, 1.5)]),
            tmp_path / "cross_section.html",
            head_frames=[np.array([[8.0, 9.0]]), np.array([[8.5, 9.5]])],
            labels=["initial", "later"],
            show_legend=False,
            dpi=70,
        )
        assert result.labels == ("initial", "later")
        _assert_standalone_slider(result.path, expected_frames=2)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_the_model_namespace_is_gone_but_its_exporters_are_not(tmp_path):
    """8.6b deletes `ModelVisualization` -- a thin namespace over module-level
    functions that were ALREADY public. The functions stay: they render through
    flopy's `PlotMapView` (grid lines, contour overlays, a `ModelMapStyle`),
    which `animate(backend="png")` does NOT reproduce, so they are a distinct
    picture rather than a duplicate spelling.
    """

    from myflopy.modflow.mf6 import interactive_plotting
    from myflopy.modflow.mf6.simulation.base import SimulationBase

    assert not hasattr(interactive_plotting, "ModelVisualization")
    assert not hasattr(SimulationBase, "visualize")

    workspace = _workspace("model_visualize")
    try:
        model = _two_cell_model("model_visualize", workspace)
        result = export_head_map_slider_html(
            model,
            tmp_path / "bound_map.html",
            head_frames=[np.array([[8.0, 9.0]])],
            style=ModelMapStyle(show_contours=False, dpi=70),
        )
        assert result.frame_count == 1
        _assert_standalone_slider(result.path, expected_frames=1)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_a_non_choropleth_animation_falls_through_to_plotly(tmp_path):
    """`FrameAnimation.html` writes the restyle page only for choropleths.

    The restyle controller reads `frame.data[0].z`, so it only means anything for
    a map. A section animation has Scatter frames and must take plotly's own
    writer instead -- the branch this exercises.
    """
    import plotly.graph_objects as go

    from myflopy.viz import FrameAnimation

    class _ScatterPicture:
        def __init__(self, ys):
            self.fig = go.Figure(go.Scatter(x=[0, 1], y=ys))

    animation = FrameAnimation(
        [("first", _ScatterPicture([1, 2])), ("second", _ScatterPicture([2, 3]))]
    )
    output = animation.html(tmp_path / "section_animation.html")

    assert output.exists()
    text = output.read_text(encoding="utf-8")
    assert "plotly" in text.lower()
    assert "smApplyFrame" not in text          # the restyle controller is map-only
    assert len(animation.fig.frames) == 2


def test_plotly_animation_export_can_select_frames(tmp_path):
    """Frame subselection survives `ModelVisualization`'s deletion.

    It was never that class's logic -- `_select_plotly_frames` did the work and
    still does. What went is the method that called it.
    """
    import plotly.graph_objects as go

    frames = [
        go.Frame(name=str(index), data=[go.Scatter(x=[0, 1], y=[index, index + 1])])
        for index in range(6)
    ]
    figure = go.Figure(data=frames[0].data, frames=frames)
    figure.update_layout(
        sliders=[
            {
                "steps": [
                    {"label": str(index), "method": "animate", "args": [[str(index)]]}
                    for index in range(6)
                ]
            }
        ]
    )

    result = _select_plotly_frames(figure, frame_stride=2, max_frames=2)

    assert [frame.name for frame in result.frames] == ["0", "2"]
    assert [step.label for step in result.layout.sliders[0].steps] == ["0", "2"]


def test_plotly_cross_section_frames_preserve_current_axis_zoom():
    from types import SimpleNamespace

    from myflopy.modflow.utils.datatypes.xsections import XSection

    periods = [(0, 0), (0, 1), (0, 2)]

    class DummyXSection(XSection):
        def __init__(self):
            self._model = SimpleNamespace(name="axis-lock", kstpkper=periods)
            self._kstpkper = periods[0]
            self.animation_kstpkpers = periods
            self.interpolate = False
            self.interpolator = None
            self.section_name = "axis-lock"
            self._layer = [0]
            self.show_model_top = False
            self.show_model_btm = False

        @property
        def xsect(self):
            offset = float(self.kstpkper[1])
            return [0.0, 10.0], [[100.0 + offset, 95.0 + offset]]

    figure = DummyXSection().ani

    assert figure.layout.uirevision == "lock"
    assert figure.layout.xaxis.uirevision == "lock"
    assert figure.layout.yaxis.uirevision == "lock"
    assert figure.layout.yaxis.range is not None
    assert all(frame.layout.xaxis.range is None for frame in figure.frames)
    assert all(frame.layout.yaxis.range is None for frame in figure.frames)
    assert all(frame.layout.xaxis.autorange is None for frame in figure.frames)
    assert all(frame.layout.yaxis.autorange is None for frame in figure.frames)


def test_plotly_head_map_animation_keeps_all_frames_and_slider(monkeypatch, tmp_path):
    from types import SimpleNamespace

    import plotly.graph_objects as go

    from myflopy.modflow.utils.datatypes.choros import Choro

    periods = [(0, index) for index in range(8)]

    class DummyHeads:
        kstpkper = periods

    class DummyModel:
        name = "plotly-map"
        kstpkper = periods
        hds = DummyHeads()

    choro = object.__new__(Choro)
    choro._model = DummyModel()
    choro._vor = SimpleNamespace(
        grid_centroid=SimpleNamespace(x=-122.1, y=47.2),
        gdf_latlon=SimpleNamespace(total_bounds=np.array([-122.2, 47.1, -122.0, 47.3])),
    )
    choro._kstpkper = periods[0]
    choro.fit_bounds = True
    choro.bounds_padding = 0.05
    choro.zoom = 13
    choro.get_choropleth = lambda: go.Choropleth(
        z=np.array([float(choro.kstpkper[1]), float(choro.kstpkper[1] + 1)]),
        locations=["a", "b"],
        customdata=np.array(
            [
                [choro.kstpkper[0], choro.kstpkper[1], 0],
                [choro.kstpkper[0], choro.kstpkper[1], 1],
            ]
        ),
        hovertemplate="Time Step: %{customdata[0]}<br>Stress Period: %{customdata[1]}",
        geojson={"type": "FeatureCollection", "features": []},
        colorscale="Earth",
    )

    model = object.__new__(SimulationBase)
    model._visualize = None
    received_choro_kwargs = {}

    class DummyMap:
        @property
        def ani(self):
            return choro.ani

    def build_dummy_map(**kwargs):
        received_choro_kwargs.update(kwargs)
        choro._zmin = kwargs.get("zmin")
        choro._zmax = kwargs.get("zmax")
        return DummyMap()

    monkeypatch.setattr(ModelPlots, "map", lambda _self, **kwargs: build_dummy_map(**kwargs))
    output = tmp_path / "plotly_map.html"
    # 8.6b: `plotly_head_map_animation` is gone. The animation figure and its
    # export are the picture's own -- `.ani` builds it, `.html()` writes it.
    result = model.plot.map(layer=0, zmin=-5, zmax=20).ani
    _write_plotly_choropleth_restyle_html(result, output)

    assert len(result.frames) == len(periods)
    assert len(result.layout.sliders[0].steps) == len(periods)
    assert result.data[0].geojson is not None
    assert all(frame.data[0].geojson is not None for frame in result.frames)
    assert all(frame.traces is None for frame in result.frames)
    assert all(frame.baseframe == str(periods[0]) for frame in result.frames)
    assert result.frames[0].data[0].z[0] != result.frames[-1].data[0].z[0]
    assert result.frames[0].data[0].customdata[0][1] != result.frames[-1].data[0].customdata[0][1]
    assert "Stress Period" in result.data[0].hovertemplate
    assert received_choro_kwargs["zmin"] == -5
    assert received_choro_kwargs["zmax"] == 20
    assert result.data[0].coloraxis == "coloraxis"
    assert result.data[0].zmin is None
    assert result.data[0].zmax is None
    assert result.data[0].colorscale is None
    assert result.layout.coloraxis.cmin == -5
    assert result.layout.coloraxis.cmax == 20
    assert result.layout.coloraxis.colorscale is not None
    assert result.layout.coloraxis.cauto is False
    for frame in result.frames:
        frame_json = frame.data[0].to_plotly_json()
        assert {"type", "z", "customdata", "locations", "hovertemplate", "coloraxis", "geojson"}.issubset(
            frame_json
        )
        assert frame_json["coloraxis"] == "coloraxis"
        assert "colorscale" not in frame_json
        assert "zmin" not in frame_json
        assert "zmax" not in frame_json
        assert "colorbar" not in frame_json
    assert len(result.layout.updatemenus[0].buttons[0].args) == 1
    assert result.layout.updatemenus[0].buttons[0].args[0] is None
    assert result.layout.sliders[0].steps[0].args[1]["frame"]["redraw"] is False
    assert result.layout.uirevision == "lock"
    assert result.layout.map.center.lon == pytest.approx(-122.1)
    assert result.layout.map.zoom > 0
    assert result.layout.map.bounds.west is None
    assert result.layout.map.bounds.east is None
    assert output.exists()
    output_html = output.read_text(encoding="utf-8")
    assert "Plotly.restyle" in output_html
    assert "smApplyFrame" in output_html
    assert "smPlayToken" in output_html
    assert "smApplyFrame(index).then" in output_html
    assert "setTimeout(() => smPlayFrame(index + 1, token)" in output_html
    assert "const z = frame.z.slice()" in output_html
    assert "row.slice()" in output_html
    assert "plotly_buttonclicked" in output_html
    assert '"scrollZoom": true' in output_html
    assert '"displaylogo": false' in output_html

    # `_fig`, not `.fig`: this block tests `update_layout()`'s effect on the
    # ACCUMULATOR, and under the 8.1 Picture contract `.fig` assembles the
    # traces -- which this stubbed choro has no real grid to do.
    fitted_zoom = choro._fig.layout.map.zoom
    choro.bounds_padding = 0.50
    choro.update_layout()
    assert choro._fig.layout.map.zoom < fitted_zoom
    assert choro._fig.layout.map.bounds.west is None

    # 8.6b moved this guard DOWN into `_choropleth_factory`, so it now covers
    # every map rather than only the animated head map.
    with pytest.raises(ValueError, match="zmin must be less than zmax"):
        _choropleth_factory(choro.vor, zmin=20, zmax=20)

    choro.fit_bounds = False
    choro._fig = go.Figure()
    choro.update_layout()
    assert choro._fig.layout.map.zoom == 13
    assert choro._fig.layout.map.bounds.west is None


def test_particle_tracking_scene_uses_flopy_vtk_and_pyvista(monkeypatch):
    pv = pytest.importorskip("pyvista")
    from flopy.export import vtk as vtk_module

    calls = {}
    model_mesh = pv.Cube()
    path_mesh = pv.Line((0, 0, 0), (1, 1, 1), resolution=3)
    path_mesh["time"] = np.linspace(0.0, 1.0, path_mesh.n_points)

    class DummyVtk:
        def __init__(self, **kwargs):
            calls["init"] = kwargs

        def add_model(self, gwf):
            calls["model"] = gwf

        def add_pathline_points(self, pathlines):
            calls["pathlines"] = pathlines

        def to_pyvista(self):
            return [model_mesh, path_mesh]

    monkeypatch.setattr(vtk_module, "Vtk", DummyVtk)

    class DummyModel:
        gwf = object()

    pathlines = pd.DataFrame(
        {"particleid": [0, 0], "time": [0.0, 1.0], "k": [0, 0], "x": [0.0, 1.0], "y": [0.0, 1.0], "z": [0.0, 1.0]}
    )
    scene = build_particle_tracking_scene(DummyModel(), pathlines)
    try:
        assert isinstance(scene, VtkScene)
        assert len(scene.meshes) == 2
        assert calls["model"] is DummyModel.gwf
        assert calls["pathlines"].equals(pathlines)
    finally:
        scene.scene.close()


def test_particle_tracking_scene_builds_with_real_flopy_vtk():
    pytest.importorskip("pyvista")
    workspace = _workspace("real_vtk_scene")
    try:
        model = _two_cell_model("real_vtk_scene", workspace)
        pathlines = pd.DataFrame(
            {
                "particleid": [0, 0, 0],
                "time": [0.0, 0.5, 1.0],
                "k": [0, 0, 0],
                "x": [0.25, 0.75, 1.25],
                "y": [0.5, 0.5, 0.5],
                "z": [8.0, 7.5, 7.0],
            }
        )
        scene = build_particle_tracking_scene(model, pathlines)
        try:
            assert len(scene.meshes) >= 2
            assert any("time" in mesh.array_names for mesh in scene.meshes)
        finally:
            scene.scene.close()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_particle_tracking_scene_exports_real_standalone_html(tmp_path):
    """8.5b: `ParticleTrackingScene.export_html` became `VtkScene.html` -- the
    same verb every other picture answers."""
    pv = pytest.importorskip("pyvista")

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(pv.Line((0, 0, 0), (1, 1, 1)), line_width=4)
    scene = VtkScene(plotter)
    output = scene.html(tmp_path / "particle_scene.html")
    try:
        assert output.exists()
        assert output.stat().st_size > 100_000
        text = output.read_text(encoding="utf-8")
        assert "<html" in text.lower()
        assert "vtk" in text.lower()
    finally:
        plotter.close()


def test_slider_validation_rejects_empty_frames_and_label_mismatch(tmp_path):
    with pytest.raises(ValueError, match="at least one"):
        export_matplotlib_slider_html(lambda *_: plt.figure(), [], tmp_path / "empty.html")
    with pytest.raises(ValueError, match="same length"):
        export_matplotlib_slider_html(
            lambda *_: plt.figure(),
            [1, 2],
            tmp_path / "bad_labels.html",
            labels=["one"],
        )


def test_the_pathline_scene_is_reachable_through_the_grid_verb():
    """The 8.5b entry point, end to end.

    The error paths are pinned in `test_plot_front_door`, but only this proves
    `plot.grid(model, pathlines=..., backend="vtk")` and its bound twin actually
    build a scene -- the failure mode being a verb that dispatches fine and dies
    on the first real call.
    """

    pytest.importorskip("pyvista")

    import myflopy.plot as plot
    from myflopy.viz import VtkScene

    workspace = _workspace("grid_verb_scene")
    try:
        model = _two_cell_model("grid_verb_scene", workspace)
        pathlines = pd.DataFrame(
            {
                "particleid": [0, 0, 0],
                "time": [0.0, 0.5, 1.0],
                "k": [0, 0, 0],
                "x": [0.25, 0.75, 1.25],
                "y": [0.5, 0.5, 0.5],
                "z": [8.0, 7.5, 7.0],
            }
        )
        free = plot.grid(model, pathlines=pathlines, backend="vtk")
        bound = model.plot.grid(pathlines=pathlines, backend="vtk")
        try:
            assert isinstance(free, VtkScene) and isinstance(bound, VtkScene)
            assert len(free.meshes) >= 2
        finally:
            free.scene.close()
            bound.scene.close()
    finally:
        shutil.rmtree(workspace, ignore_errors=True)
