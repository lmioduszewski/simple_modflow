"""What a picture writes is what it says it writes (ledger 481).

``Picture.html`` did not create its directory where every other picture's ``html``
did; ``FrameAnimation.html`` sent any trace carrying ``z`` -- a surface, a heatmap --
to the choropleth restyle page, whose script broke Play; ``VtkScene.save("shot")``
returned a path PyVista did not write (``shot.png``) and refused ``.PNG``, and its
``html`` ignored every keyword; ``Choro.add_overlay`` after ``.fig`` never reached the
figure; ``Choro.add_hillshade`` reused one fixed ``hillshade.png`` whatever raster made
it and raised rasterio's error for its default ``None``; ``save_raster`` ignored
``polygon_clip``/``buffer`` beside the clipped image.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import pytest
import rasterio
from rasterio.transform import from_origin

from myflopy import viz
from myflopy.modflow.utils.datatypes.choros import ChoroplethHillshadeBackground


class _Stub(viz.Picture):
    def __init__(self, trace):
        self._fig = viz.Fig()
        self._fig.add_trace(trace)

    @property
    def fig(self):
        return self._fig


def test_picture_html_creates_its_directory(tmp_path):
    target = tmp_path / "new" / "dir" / "page.html"
    assert _Stub(go.Scatter(x=[0, 1], y=[0, 1])).html(target) == target
    assert target.exists()


def test_a_surface_animation_is_not_written_as_a_choropleth_page(tmp_path, monkeypatch):
    def _refuse(*args, **kwargs):
        raise AssertionError("a surface reached the choropleth restyle page")

    monkeypatch.setattr(viz, "_write_plotly_choropleth_restyle_html", _refuse)
    frames = [(f"f{i}", _Stub(go.Surface(z=np.full((2, 2), float(i))))) for i in range(2)]
    page = viz.FrameAnimation(frames).html(tmp_path / "surface.html", include_plotlyjs="cdn")
    assert page.exists() and "surface" in page.read_text()


class _FakePlotter:
    """PyVista's ``screenshot`` reads its suffix case-sensitively, lower case only."""

    def __init__(self):
        self.calls = []

    def screenshot(self, filename, **kwargs):
        suffix = filename.rsplit(".", 1)[-1] if "." in filename.rsplit("/", 1)[-1] else ""
        if suffix not in ("png", "jpg", "jpeg", "bmp", "tif", "tiff"):
            raise ValueError(f"unsupported {suffix!r}")
        self.calls.append((filename, kwargs))
        with open(filename, "wb") as handle:
            handle.write(b"image")


def test_a_scene_screenshot_names_the_file_it_wrote(tmp_path):
    scene = viz.VtkScene(_FakePlotter())
    upper = scene.save(tmp_path / "shot.PNG", scale=2)
    assert upper.exists() and upper.name == "shot.PNG"
    with pytest.raises(ValueError, match="names no image format"):
        scene.save(tmp_path / "shot")
    with pytest.raises(TypeError, match="takes no options"):
        scene.save(tmp_path / "scene.html", scale=2)
    with pytest.raises(TypeError):
        scene.html(tmp_path / "scene.html", scale=2)


def test_an_overlay_added_after_the_figure_reaches_it(canonical_run):
    choro = canonical_run.plot.map(layer=0)
    before = len(choro.fig.data)
    choro.add_overlay(go.Scattermap(lon=[0.0], lat=[0.0], name="late"))
    assert len(choro.fig.data) == before + 1
    assert choro.fig.data[-1].name == "late"


def _tif(path, value):
    with rasterio.open(
        path, "w", driver="GTiff", height=2, width=2, count=1, dtype="float32",
        crs="EPSG:2927", transform=from_origin(0, 2, 1, 1),
    ) as dst:
        dst.write(np.array([[[0, value], [value, 0]]], dtype="float32"))
    return path


def test_a_hillshade_is_drawn_from_its_own_raster_and_writes_nothing(
    canonical_run, tmp_path, monkeypatch,
):
    monkeypatch.chdir(tmp_path)
    first = ChoroplethHillshadeBackground(_tif(tmp_path / "a.tif", 1.0)).data_uri()
    second = ChoroplethHillshadeBackground(_tif(tmp_path / "b.tif", -5.0)).data_uri()
    assert first.startswith("data:image/png;base64,") and first != second
    output = canonical_run.model_output_folder_path
    before = sorted(path.name for path in output.iterdir())
    choro = canonical_run.plot.map(layer=0)
    choro.fig  # assembled
    choro.add_hillshade(tmp_path / "a.tif")
    assert choro.fig.layout.map.layers[0].source == first
    assert sorted(path.name for path in tmp_path.iterdir()) == ["a.tif", "b.tif"]
    assert sorted(path.name for path in output.iterdir()) == before
    with pytest.raises(TypeError, match="pass tif_path="):
        choro.add_hillshade()
