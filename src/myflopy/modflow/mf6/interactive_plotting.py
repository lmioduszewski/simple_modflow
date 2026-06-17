"""Standalone interactive visualization helpers for MF6 results.

The matplotlib exporters intentionally pre-render figures and use a small
browser-side slider. The resulting HTML files do not require Jupyter, a Python
kernel, or a web server.
"""

from __future__ import annotations

import base64
import copy
import html
import io
import json
import os
from contextlib import contextmanager
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
from flopy.plot import PlotMapView

from myflopy.modflow.mf6.cross_section_plotting import (
    ModelCrossSectionStyle,
    plot_model_cross_section,
)

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from myflopy.modflow.mf6.simulation.base import SimulationBase


_DEFAULT_PLOTLY_CONFIG = {
    "scrollZoom": True,
    "displaylogo": False,
}
_DEFAULT_MOSAIC_DPI = 200


def _plotly_config(fig, config: dict[str, Any] | None = None) -> dict[str, Any]:
    """Merge the standard figs Plotly config with optional export overrides."""

    merged = dict(_DEFAULT_PLOTLY_CONFIG)
    figure_config = getattr(fig, "_config", None)
    if figure_config:
        merged.update(figure_config)
    if config:
        merged.update(config)
    return merged


@dataclass(frozen=True)
class ModelMapStyle:
    """Styling controls shared by standalone matplotlib head-map sliders."""

    figsize: tuple[float, float] = (10, 8)
    cmap: str = "viridis"
    grid_color: str = "#3c4652"
    grid_linewidth: float = 0.2
    contour_color: str = "white"
    contour_linewidth: float = 0.7
    contour_levels: int | Sequence[float] = 10
    show_grid: bool = True
    show_contours: bool = True
    show_colorbar: bool = True
    dpi: int = 140


@dataclass(frozen=True)
class StandaloneHtmlSlider:
    """Description of a generated standalone frame-slider HTML document."""

    path: Path
    labels: tuple[str, ...]
    frame_count: int
    title: str
    embedded_frames: bool = True
    frame_directory: Path | None = None
    rendered_frames: int = 0
    reused_frames: int = 0


@dataclass(frozen=True)
class FrameExportProgress:
    """Progress event emitted while rendering a standalone slider."""

    index: int
    total: int
    label: str
    status: str
    path: Path | None = None


@dataclass
class ParticleTrackingScene:
    """FloPy VTK / PyVista particle-tracking scene."""

    plotter: Any
    meshes: tuple[Any, ...]

    def export_html(self, path: str | Path) -> Path:
        """Export the interactive PyVista scene as standalone HTML."""

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            self.plotter.export_html(path)
        except ImportError as error:
            raise ImportError(
                "PyVista standalone HTML export requires its trame dependencies. "
                "Install them with `python -m pip install trame trame-vtk trame-vuetify`."
            ) from error
        return path


def _figure_to_data_uri(fig: "Figure", *, dpi: int) -> str:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def _write_figure_png(fig: "Figure", path: Path, *, dpi: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    fig.savefig(temporary, format="png", dpi=dpi, bbox_inches="tight")
    temporary.replace(path)
    return path


def _write_text_atomic(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


@contextmanager
def _export_lock(output_path: Path):
    lock_path = output_path.with_suffix(f"{output_path.suffix}.lock")
    try:
        descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as error:
        raise RuntimeError(
            f"Another export appears to be using {output_path}. "
            f"Remove the stale lock only after confirming no export is running: {lock_path}"
        ) from error
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as lock:
            lock.write(str(os.getpid()))
        yield
    finally:
        lock_path.unlink(missing_ok=True)


def _emit_progress(progress, event: FrameExportProgress) -> None:
    if callable(progress):
        progress(event)
    elif progress:
        print(f"[{event.index + 1}/{event.total}] {event.status}: {event.label}", flush=True)


def _select_sequence(
    values: Sequence[Any],
    labels: Sequence[str],
    *,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
) -> tuple[list[Any], list[str]]:
    if frame_stride < 1:
        raise ValueError("frame_stride must be at least 1.")
    if max_frames is not None and max_frames < 1:
        raise ValueError("max_frames must be at least 1.")
    indices = list(range(len(values))) if frame_indices is None else [int(index) for index in frame_indices]
    if any(index < 0 or index >= len(values) for index in indices):
        raise IndexError("frame_indices contains an index outside the available frame range.")
    indices = indices[::frame_stride]
    if max_frames is not None:
        indices = indices[:max_frames]
    if not indices:
        raise ValueError("Frame selection must contain at least one frame.")
    return [values[index] for index in indices], [str(labels[index]) for index in indices]


def _slider_html(*, frames: Sequence[str], labels: Sequence[str], title: str, interval_ms: int) -> str:
    frame_json = json.dumps(list(frames))
    label_json = json.dumps([str(label) for label in labels])
    safe_title = html.escape(title)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{safe_title}</title>
<style>
:root {{ color-scheme: dark; --bg:#10151c; --card:#18212b; --line:#344352; --text:#e9f0f5; --accent:#55c2a3; }}
* {{ box-sizing: border-box; }}
body {{ margin:0; background:var(--bg); color:var(--text); font:15px/1.4 Segoe UI, sans-serif; }}
main {{ width:min(1400px, 100%); margin:auto; padding:18px; }}
.card {{ background:var(--card); border:1px solid var(--line); border-radius:12px; padding:14px; box-shadow:0 12px 34px #0006; }}
h1 {{ margin:0 0 12px; font-size:20px; }}
img {{ display:block; width:100%; height:auto; background:white; border-radius:7px; }}
.controls {{ display:grid; grid-template-columns:auto 1fr auto; gap:12px; align-items:center; margin-top:12px; }}
button {{ border:1px solid var(--line); border-radius:6px; background:#22303c; color:var(--text); padding:7px 13px; cursor:pointer; }}
input[type=range] {{ width:100%; accent-color:var(--accent); }}
#frame-label {{ min-width:120px; text-align:right; font-variant-numeric:tabular-nums; }}
</style>
</head>
<body><main><section class="card">
<h1>{safe_title}</h1>
<img id="frame" alt="{safe_title}">
<div class="controls">
<button id="play" type="button">Play</button>
<input id="slider" type="range" min="0" max="{len(frames) - 1}" value="0" step="1">
<span id="frame-label"></span>
</div>
</section></main>
<script>
const frames = {frame_json};
const labels = {label_json};
const image = document.getElementById("frame");
const slider = document.getElementById("slider");
const label = document.getElementById("frame-label");
const play = document.getElementById("play");
let timer = null;
function show(index) {{
  const i = Number(index);
  image.src = frames[i];
  label.textContent = labels[i];
  slider.value = i;
}}
slider.addEventListener("input", event => show(event.target.value));
play.addEventListener("click", () => {{
  if (timer !== null) {{ clearInterval(timer); timer = null; play.textContent = "Play"; return; }}
  play.textContent = "Pause";
  timer = setInterval(() => show((Number(slider.value) + 1) % frames.length), {int(interval_ms)});
}});
show(0);
</script></body></html>
"""


def export_matplotlib_slider_html(
    render_frame: Callable[[Any, int], "Figure | tuple[Figure, Any]"],
    frame_values: Sequence[Any],
    output_path: str | Path,
    *,
    labels: Sequence[str] | None = None,
    title: str = "Model results through time",
    dpi: int = 140,
    interval_ms: int = 700,
    close_figures: bool = True,
    embed_frames: bool = True,
    frame_directory: str | Path | None = None,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
    resume: bool = False,
    progress: bool | Callable[[FrameExportProgress], None] = False,
) -> StandaloneHtmlSlider:
    """Render matplotlib frames into a browser slider.

    Embedded frames produce one self-contained HTML file. For long simulations
    or large figures, ``embed_frames=False`` writes PNG assets beside the HTML
    and avoids holding all base64-encoded frames in memory.
    """

    values = list(frame_values)
    if not values:
        raise ValueError("frame_values must contain at least one frame.")
    resolved_labels = [str(value) for value in values] if labels is None else [str(value) for value in labels]
    if len(resolved_labels) != len(values):
        raise ValueError("labels must have the same length as frame_values.")
    values, resolved_labels = _select_sequence(
        values,
        resolved_labels,
        frame_indices=frame_indices,
        frame_stride=frame_stride,
        max_frames=max_frames,
    )
    if resume and embed_frames:
        raise ValueError("resume=True requires embed_frames=False.")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    resolved_frame_directory = None
    if not embed_frames:
        resolved_frame_directory = (
            output_path.parent / f"{output_path.stem}_frames"
            if frame_directory is None
            else Path(frame_directory)
        )
        resolved_frame_directory.mkdir(parents=True, exist_ok=True)

    rendered_count = 0
    reused_count = 0
    frames: list[str] = []
    with _export_lock(output_path):
        if resolved_frame_directory is not None:
            manifest_path = resolved_frame_directory / "frames.json"
            manifest = {"labels": resolved_labels, "dpi": int(dpi), "frame_count": len(values)}
            if resume and manifest_path.exists():
                existing = json.loads(manifest_path.read_text(encoding="utf-8"))
                if existing != manifest:
                    raise ValueError(
                        "Existing frame manifest does not match this export. "
                        "Use a different frame_directory or restart without resume=True."
                    )
            _write_text_atomic(manifest_path, json.dumps(manifest, indent=2))

        for index, value in enumerate(values):
            frame_path = (
                None
                if resolved_frame_directory is None
                else resolved_frame_directory / f"frame_{index:05d}.png"
            )
            if resume and frame_path is not None and frame_path.exists() and frame_path.stat().st_size > 0:
                reused_count += 1
                frames.append(Path(os.path.relpath(frame_path, output_path.parent)).as_posix())
                _emit_progress(
                    progress,
                    FrameExportProgress(index, len(values), resolved_labels[index], "reused", frame_path),
                )
                continue

            rendered = render_frame(value, index)
            fig = rendered[0] if isinstance(rendered, tuple) else rendered
            try:
                if embed_frames:
                    frames.append(_figure_to_data_uri(fig, dpi=dpi))
                else:
                    _write_figure_png(fig, frame_path, dpi=dpi)
                    frames.append(Path(os.path.relpath(frame_path, output_path.parent)).as_posix())
                rendered_count += 1
                _emit_progress(
                    progress,
                    FrameExportProgress(index, len(values), resolved_labels[index], "rendered", frame_path),
                )
            finally:
                if close_figures:
                    plt.close(fig)

        _write_text_atomic(
            output_path,
            _slider_html(frames=frames, labels=resolved_labels, title=title, interval_ms=interval_ms),
        )
    return StandaloneHtmlSlider(
        output_path,
        tuple(resolved_labels),
        len(frames),
        title,
        embedded_frames=embed_frames,
        frame_directory=resolved_frame_directory,
        rendered_frames=rendered_count,
        reused_frames=reused_count,
    )


def _resolve_frames(
    model: "SimulationBase",
    *,
    kstpkpers: Sequence[tuple[int, int]] | None,
    head_frames: Sequence[Any] | None,
) -> tuple[list[Any], list[tuple[int, int]] | None]:
    if head_frames is not None:
        values = [_normalize_head_frame(frame) for frame in head_frames]
        if not values:
            raise ValueError("head_frames must contain at least one array.")
        return values, None
    if kstpkpers is None:
        kstpkpers = list(model.gwf.output.head().get_kstpkper())
    values = [tuple(value) for value in kstpkpers]
    if not values:
        raise ValueError("No head-output times were found.")
    return values, values


def _frame_labels(
    model: "SimulationBase",
    values: Sequence[Any],
    *,
    kstpkpers: Sequence[tuple[int, int]] | None,
    labels: Sequence[str] | None,
) -> list[str]:
    if labels is not None:
        return [str(label) for label in labels]
    if kstpkpers is None:
        return [f"Frame {index + 1}" for index in range(len(values))]
    return [f"kstpkper={tuple(value)}" for value in kstpkpers]


def _normalize_head_frame(values):
    array = np.asarray(values)
    if array.ndim == 3 and array.shape[1] == 1:
        return array[:, 0, :]
    return array


def _head_frame(model: "SimulationBase", value, resolved_kstpkpers):
    data = value if resolved_kstpkpers is None else model.gwf.output.head().get_data(kstpkper=value)
    return _normalize_head_frame(data)


def _clean_head_values(values):
    array = np.asarray(values, dtype=float)
    return np.ma.masked_where(~np.isfinite(array) | (np.abs(array) >= 1.0e29), array)


def _shared_head_limits(
    model: "SimulationBase",
    values: Sequence[Any],
    resolved_kstpkpers,
    layers: Sequence[int],
) -> tuple[float | None, float | None]:
    minimum = np.inf
    maximum = -np.inf
    for value in values:
        frame = np.asarray(_head_frame(model, value, resolved_kstpkpers), dtype=float)
        for layer in layers:
            layer_values = frame if frame.ndim == 1 else frame[layer]
            valid = layer_values[np.isfinite(layer_values) & (np.abs(layer_values) < 1.0e29)]
            if valid.size:
                minimum = min(minimum, float(valid.min()))
                maximum = max(maximum, float(valid.max()))
    if not np.isfinite(minimum) or not np.isfinite(maximum):
        return None, None
    return minimum, maximum


def export_cross_section_slider_html(
    model: "SimulationBase",
    line,
    output_path: str | Path,
    *,
    kstpkpers: Sequence[tuple[int, int]] | None = None,
    head_frames: Sequence[Any] | None = None,
    labels: Sequence[str] | None = None,
    style: ModelCrossSectionStyle | None = None,
    title: str = "Head cross section through time",
    dpi: int = 140,
    interval_ms: int = 700,
    embed_frames: bool = True,
    frame_directory: str | Path | None = None,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
    resume: bool = False,
    progress: bool | Callable[[FrameExportProgress], None] = False,
    **cross_section_kwargs,
) -> StandaloneHtmlSlider:
    """Export existing ``plot_model_cross_section`` styling as a standalone slider."""

    values, resolved_kstpkpers = _resolve_frames(model, kstpkpers=kstpkpers, head_frames=head_frames)
    frame_labels = _frame_labels(model, values, kstpkpers=resolved_kstpkpers, labels=labels)
    values, frame_labels = _select_sequence(
        values,
        frame_labels,
        frame_indices=frame_indices,
        frame_stride=frame_stride,
        max_frames=max_frames,
    )
    if resolved_kstpkpers is not None:
        resolved_kstpkpers = values

    def render(value, index):
        kwargs = dict(cross_section_kwargs)
        kwargs["title"] = frame_labels[index]
        if resolved_kstpkpers is None:
            kwargs["head_data"] = value
        else:
            kwargs["kstpkper"] = value
        return plot_model_cross_section(model, line, style=style, **kwargs)

    return export_matplotlib_slider_html(
        render,
        values,
        output_path,
        labels=frame_labels,
        title=title,
        dpi=dpi,
        interval_ms=interval_ms,
        embed_frames=embed_frames,
        frame_directory=frame_directory,
        resume=resume,
        progress=progress,
    )


def plot_model_head_map(
    model: "SimulationBase",
    head_data,
    *,
    layer: int = 0,
    ax=None,
    style: ModelMapStyle | None = None,
    title: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
):
    """Render one FloPy-style head map using the canonical matplotlib style."""

    style = ModelMapStyle() if style is None else style
    if ax is None:
        fig, ax = plt.subplots(figsize=style.figsize)
    else:
        fig = ax.figure
    array = _normalize_head_frame(head_data)
    layer_data = _clean_head_values(array if array.ndim == 1 else array[layer])
    view = PlotMapView(model=model.gwf, modelgrid=model.gwf.modelgrid, layer=layer, ax=ax)
    image = view.plot_array(layer_data, cmap=style.cmap, vmin=vmin, vmax=vmax)
    if style.show_grid:
        view.plot_grid(color=style.grid_color, linewidth=style.grid_linewidth)
    if style.show_contours and layer_data.count() > 1 and float(layer_data.max()) > float(layer_data.min()):
        contours = view.contour_array(
            layer_data,
            levels=style.contour_levels,
            colors=style.contour_color,
            linewidths=style.contour_linewidth,
        )
        ax.clabel(contours, inline=True, fontsize=7)
    if style.show_colorbar:
        fig.colorbar(image, ax=ax, shrink=0.8, label="Head")
    if title:
        ax.set_title(title)
    ax.set_aspect("equal")
    return fig, ax


def plot_particle_pathlines(
    model: "SimulationBase",
    pathlines,
    *,
    layer: int | str = "all",
    ax=None,
    head_data=None,
    head_layer: int = 0,
    pathline_color: str | None = None,
    pathline_alpha: float = 0.7,
    pathline_linewidth: float = 1.0,
    show_grid: bool = True,
    grid_alpha: float = 0.25,
    title: str = "Particle pathlines",
):
    """Plot PRT, MODPATH, or compatible particle pathlines on a FloPy map."""

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig = ax.figure
    view = PlotMapView(model=model.gwf, modelgrid=model.gwf.modelgrid, layer=head_layer, ax=ax)
    if head_data is not None:
        array = np.asarray(head_data)
        view.plot_array(array if array.ndim == 1 else array[head_layer], alpha=0.35)
    if show_grid:
        view.plot_grid(alpha=grid_alpha)
    kwargs = {"layer": layer, "alpha": pathline_alpha, "linewidth": pathline_linewidth}
    if pathline_color is not None:
        kwargs["colors"] = pathline_color
    view.plot_pathline(pathlines, **kwargs)
    ax.set_title(title)
    ax.set_aspect("equal")
    return fig, ax


def export_head_map_slider_html(
    model: "SimulationBase",
    output_path: str | Path,
    *,
    layer: int = 0,
    kstpkpers: Sequence[tuple[int, int]] | None = None,
    head_frames: Sequence[Any] | None = None,
    labels: Sequence[str] | None = None,
    style: ModelMapStyle | None = None,
    title: str = "Head map through time",
    vmin: float | None = None,
    vmax: float | None = None,
    shared_color_scale: bool = True,
    embed_frames: bool = True,
    frame_directory: str | Path | None = None,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
    resume: bool = False,
    progress: bool | Callable[[FrameExportProgress], None] = False,
) -> StandaloneHtmlSlider:
    """Export a FloPy-style head map with a standalone browser time slider."""

    style = ModelMapStyle() if style is None else style
    values, resolved_kstpkpers = _resolve_frames(model, kstpkpers=kstpkpers, head_frames=head_frames)
    frame_labels = _frame_labels(model, values, kstpkpers=resolved_kstpkpers, labels=labels)
    values, frame_labels = _select_sequence(
        values,
        frame_labels,
        frame_indices=frame_indices,
        frame_stride=frame_stride,
        max_frames=max_frames,
    )
    if resolved_kstpkpers is not None:
        resolved_kstpkpers = values
    if shared_color_scale and (vmin is None or vmax is None):
        shared_min, shared_max = _shared_head_limits(model, values, resolved_kstpkpers, [layer])
        vmin = shared_min if vmin is None else vmin
        vmax = shared_max if vmax is None else vmax

    def render(value, index):
        data = _head_frame(model, value, resolved_kstpkpers)
        return plot_model_head_map(
            model,
            data,
            layer=layer,
            style=style,
            title=frame_labels[index],
            vmin=vmin,
            vmax=vmax,
        )

    return export_matplotlib_slider_html(
        render,
        values,
        output_path,
        labels=frame_labels,
        title=title,
        dpi=style.dpi,
        embed_frames=embed_frames,
        frame_directory=frame_directory,
        resume=resume,
        progress=progress,
    )


def export_head_layer_mosaic_slider_html(
    model: "SimulationBase",
    output_path: str | Path,
    *,
    layers: Sequence[int] | None = None,
    kstpkpers: Sequence[tuple[int, int]] | None = None,
    head_frames: Sequence[Any] | None = None,
    labels: Sequence[str] | None = None,
    style: ModelMapStyle | None = None,
    title: str = "Head layer mosaic through time",
    ncols: int = 3,
    dpi: int | None = None,
    panel_figsize: tuple[float, float] | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    shared_color_scale: bool = True,
    embed_frames: bool = True,
    frame_directory: str | Path | None = None,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
    resume: bool = False,
    progress: bool | Callable[[FrameExportProgress], None] = False,
) -> StandaloneHtmlSlider:
    """Export the layer-mosaic slider pattern used in the FloPy training notebooks.

    ``dpi`` controls the PNG resolution directly. ``panel_figsize`` controls the
    physical size of each layer panel before rasterization.
    """

    if style is None:
        style = ModelMapStyle(dpi=_DEFAULT_MOSAIC_DPI)
    resolved_dpi = style.dpi if dpi is None else int(dpi)
    resolved_panel_figsize = style.figsize if panel_figsize is None else panel_figsize
    if resolved_dpi <= 0:
        raise ValueError("dpi must be greater than zero.")
    if len(resolved_panel_figsize) != 2 or any(float(value) <= 0 for value in resolved_panel_figsize):
        raise ValueError("panel_figsize must contain two positive values.")
    layers = list(range(model.gwf.modelgrid.nlay)) if layers is None else [int(layer) for layer in layers]
    if not layers:
        raise ValueError("layers must contain at least one layer.")
    values, resolved_kstpkpers = _resolve_frames(model, kstpkpers=kstpkpers, head_frames=head_frames)
    frame_labels = _frame_labels(model, values, kstpkpers=resolved_kstpkpers, labels=labels)
    values, frame_labels = _select_sequence(
        values,
        frame_labels,
        frame_indices=frame_indices,
        frame_stride=frame_stride,
        max_frames=max_frames,
    )
    if resolved_kstpkpers is not None:
        resolved_kstpkpers = values
    if shared_color_scale and (vmin is None or vmax is None):
        shared_min, shared_max = _shared_head_limits(model, values, resolved_kstpkpers, layers)
        vmin = shared_min if vmin is None else vmin
        vmax = shared_max if vmax is None else vmax

    def render(value, index):
        data = _head_frame(model, value, resolved_kstpkpers)
        nrows = int(np.ceil(len(layers) / ncols))
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(
                float(resolved_panel_figsize[0]) * ncols,
                float(resolved_panel_figsize[1]) * nrows,
            ),
        )
        axes = np.atleast_1d(axes).reshape(-1)
        mosaic_style = replace(style, show_colorbar=False)
        images = []
        for ax, layer in zip(axes, layers):
            plot_model_head_map(
                model,
                data,
                layer=layer,
                ax=ax,
                style=mosaic_style,
                title=f"Layer {layer + 1}",
                vmin=vmin,
                vmax=vmax,
            )
            images.append(ax.collections[0])
        for ax in axes[len(layers) :]:
            ax.set_visible(False)
        if style.show_colorbar and images:
            fig.colorbar(images[0], ax=axes[: len(layers)].tolist(), shrink=0.75, label="Head")
        fig.suptitle(frame_labels[index])
        return fig

    return export_matplotlib_slider_html(
        render,
        values,
        output_path,
        labels=frame_labels,
        title=title,
        dpi=resolved_dpi,
        embed_frames=embed_frames,
        frame_directory=frame_directory,
        resume=resume,
        progress=progress,
    )


def build_particle_tracking_scene(
    model: "SimulationBase",
    pathlines,
    *,
    vertical_exaggeration: float = 1.0,
    model_style: str = "wireframe",
    model_opacity: float = 0.25,
    pathline_cmap: str = "viridis",
    pathline_width: float = 4.0,
    show_edges: bool = True,
    off_screen: bool = True,
) -> ParticleTrackingScene:
    """Build the FloPy VTK / PyVista particle scene used by the training notebooks."""

    try:
        import pyvista as pv
    except ImportError as error:
        raise ImportError("Particle-tracking 3D scenes require pyvista.") from error
    from flopy.export.vtk import Vtk

    vtk = Vtk(
        model=model.gwf,
        binary=False,
        vertical_exageration=vertical_exaggeration,
        smooth=False,
    )
    vtk.add_model(model.gwf)
    vtk.add_pathline_points(copy.deepcopy(pathlines))
    raw_meshes = vtk.to_pyvista()
    meshes = tuple(raw_meshes if isinstance(raw_meshes, list) else [raw_meshes])
    plotter = pv.Plotter(off_screen=off_screen)
    for index, mesh in enumerate(meshes):
        is_pathline = index == len(meshes) - 1 and len(meshes) > 1
        if is_pathline:
            scalars = next((name for name in ("time", "t") if name in mesh.array_names), None)
            plotter.add_mesh(
                mesh,
                scalars=scalars,
                cmap=pathline_cmap,
                line_width=pathline_width,
                render_lines_as_tubes=True,
            )
        else:
            plotter.add_mesh(
                mesh,
                style=model_style,
                opacity=model_opacity,
                show_edges=show_edges,
                color="#c9d6df",
            )
    plotter.add_axes()
    plotter.reset_camera()
    return ParticleTrackingScene(plotter=plotter, meshes=meshes)


def export_particle_tracking_html(
    model: "SimulationBase",
    pathlines,
    output_path: str | Path,
    **scene_kwargs,
) -> Path:
    """Build and export a standalone FloPy VTK / PyVista particle scene."""

    scene = build_particle_tracking_scene(model, pathlines, **scene_kwargs)
    try:
        return scene.export_html(output_path)
    finally:
        scene.plotter.close()


def _select_plotly_frames(
    fig,
    *,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
):
    frames = list(fig.frames)
    if not frames:
        return fig
    if frame_indices is None and frame_stride == 1 and max_frames is None:
        return fig
    selected, _ = _select_sequence(
        frames,
        [str(frame.name) for frame in frames],
        frame_indices=frame_indices,
        frame_stride=frame_stride,
        max_frames=max_frames,
    )
    import plotly.graph_objects as go

    figure_type = type(fig)
    try:
        fig = figure_type(data=selected[0].data, frames=selected, layout=fig.layout)
    except (TypeError, ValueError):
        fig = go.Figure(data=selected[0].data, frames=selected, layout=fig.layout)
    selected_names = {str(frame.name) for frame in selected}
    if fig.layout.sliders:
        for slider in fig.layout.sliders:
            slider.steps = tuple(
                step
                for step in slider.steps
                if str(step.label) in selected_names
            )
    return fig


def _write_plotly_choropleth_restyle_html(
    fig,
    output_path: str | Path,
    *,
    include_plotlyjs: bool | str = True,
    interval_ms: int = 250,
    config: dict[str, Any] | None = None,
) -> Path:
    """Export a choropleth animation using repeatable trace-only restyles."""

    import plotly.graph_objects as go
    import plotly.io as pio
    from plotly.utils import PlotlyJSONEncoder

    if not fig.frames:
        raise ValueError("A choropleth restyle export requires at least one frame.")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame_payload = [
        {
            "name": str(frame.name),
            "z": frame.data[0].z,
            "customdata": frame.data[0].customdata,
        }
        for frame in fig.frames
    ]
    payload_json = json.dumps(frame_payload, cls=PlotlyJSONEncoder)

    export_fig = go.Figure(data=fig.data, layout=fig.layout)
    export_fig.frames = ()
    export_fig.update_layout(
        updatemenus=[
            {
                "type": "buttons",
                "buttons": [
                    {"label": "Play", "method": "skip", "args": []},
                    {"label": "Pause", "method": "skip", "args": []},
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "currentvalue": {"prefix": "Frame:"},
                "steps": [
                    {"label": frame["name"], "method": "skip", "args": []}
                    for frame in frame_payload
                ],
            }
        ],
    )
    post_script = f"""
const smGraph = document.getElementById("{{plot_id}}");
const smFrames = {payload_json};
let smFrameIndex = 0;
let smTimer = null;
let smPlayToken = 0;

function smStop() {{
  smPlayToken += 1;
  if (smTimer !== null) {{
    clearTimeout(smTimer);
    smTimer = null;
  }}
}}

function smApplyFrame(index) {{
  smFrameIndex = index;
  const frame = smFrames[index];
  const z = frame.z.slice();
  const customdata = frame.customdata.map(
    (row) => Array.isArray(row) ? row.slice() : row
  );
  return Plotly.restyle(
    smGraph,
    {{z: [z], customdata: [customdata]}},
    [0]
  ).then(() => Plotly.relayout(smGraph, {{"sliders[0].active": index}}));
}}

function smPlayFrame(index, token) {{
  if (token !== smPlayToken) return;
  smApplyFrame(index).then(() => {{
    if (token !== smPlayToken || index >= smFrames.length - 1) return;
    smTimer = setTimeout(() => smPlayFrame(index + 1, token), {int(interval_ms)});
  }});
}}

smGraph.on("plotly_sliderchange", (event) => {{
  smStop();
  const index = smFrames.findIndex((frame) => frame.name === String(event.step.label));
  if (index >= 0) smApplyFrame(index);
}});

smGraph.on("plotly_buttonclicked", (event) => {{
  if (event.button.label === "Pause") {{
    smStop();
    return;
  }}
  if (event.button.label !== "Play") return;
  smStop();
  const token = smPlayToken;
  smPlayFrame(0, token);
}});
smGraph.dataset.simpleModflowRestyleReady = "true";
"""
    pio.write_html(
        export_fig,
        file=output_path,
        include_plotlyjs=include_plotlyjs,
        post_script=post_script,
        config=_plotly_config(fig, config),
        auto_play=False,
        auto_open=False,
    )
    return output_path


class ModelVisualization:
    """Model-bound entry point for canonical standalone visualization exports."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def cross_section_slider_html(self, line, output_path, **kwargs) -> StandaloneHtmlSlider:
        return export_cross_section_slider_html(self.model, line, output_path, **kwargs)

    def head_map_slider_html(self, output_path, **kwargs) -> StandaloneHtmlSlider:
        return export_head_map_slider_html(self.model, output_path, **kwargs)

    def head_layer_mosaic_slider_html(self, output_path, **kwargs) -> StandaloneHtmlSlider:
        return export_head_layer_mosaic_slider_html(self.model, output_path, **kwargs)

    def particle_tracking_scene(self, pathlines, **kwargs) -> ParticleTrackingScene:
        return build_particle_tracking_scene(self.model, pathlines, **kwargs)

    def particle_tracking_html(self, pathlines, output_path, **kwargs) -> Path:
        return export_particle_tracking_html(self.model, pathlines, output_path, **kwargs)

    def plotly_cross_section_animation(
        self,
        *,
        output_path: str | Path | None = None,
        frame_indices: Sequence[int] | None = None,
        frame_stride: int = 1,
        max_frames: int | None = None,
        include_plotlyjs: bool | str = True,
        config: dict[str, Any] | None = None,
        **xsection_kwargs,
    ):
        """Build a Plotly cross-section animation and optionally export HTML.

        HTML exports use the standard figs Plotly config by default. Supply
        ``config`` to override or extend options such as ``scrollZoom``.
        """

        selection_built_in = False
        if frame_indices is not None or frame_stride != 1 or max_frames is not None:
            try:
                periods = list(self.model.kstpkper)
                selected, _ = _select_sequence(
                    periods,
                    [str(period) for period in periods],
                    frame_indices=frame_indices,
                    frame_stride=frame_stride,
                    max_frames=max_frames,
                )
                xsection_kwargs["animation_kstpkpers"] = selected
                selection_built_in = True
            except (AttributeError, TypeError):
                pass
        fig = self.model.xs(**xsection_kwargs).ani
        if not selection_built_in:
            fig = _select_plotly_frames(
                fig,
                frame_indices=frame_indices,
                frame_stride=frame_stride,
                max_frames=max_frames,
            )
        if output_path is not None:
            import plotly.io as pio

            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            pio.write_html(
                fig,
                file=output_path,
                include_plotlyjs=include_plotlyjs,
                config=_plotly_config(fig, config),
                auto_open=False,
            )
        return fig

    def plotly_head_map_animation(
        self,
        *,
        output_path: str | Path | None = None,
        frame_indices: Sequence[int] | None = None,
        frame_stride: int = 1,
        max_frames: int | None = None,
        zmin: float | int | None = None,
        zmax: float | int | None = None,
        include_plotlyjs: bool | str = True,
        config: dict[str, Any] | None = None,
        **choro_kwargs,
    ):
        """Build a Plotly head-map animation and optionally export HTML.

        ``zmin`` and ``zmax`` define a fixed global colorscale range for every
        frame. When omitted, the range is calculated across selected frames.
        HTML exports use the standard figs Plotly config by default.
        """

        if zmin is not None and zmax is not None and zmin >= zmax:
            raise ValueError("zmin must be less than zmax.")
        choro_kwargs["zmin"] = zmin
        choro_kwargs["zmax"] = zmax

        selection_built_in = False
        if frame_indices is not None or frame_stride != 1 or max_frames is not None:
            try:
                periods = list(self.model.kstpkper)
                selected, _ = _select_sequence(
                    periods,
                    [str(period) for period in periods],
                    frame_indices=frame_indices,
                    frame_stride=frame_stride,
                    max_frames=max_frames,
                )
                choro_kwargs["animation_kstpkpers"] = selected
                selection_built_in = True
            except (AttributeError, TypeError):
                pass
        fig = self.model.cor(**choro_kwargs).ani
        if not selection_built_in:
            fig = _select_plotly_frames(
                fig,
                frame_indices=frame_indices,
                frame_stride=frame_stride,
                max_frames=max_frames,
            )
        if output_path is not None:
            _write_plotly_choropleth_restyle_html(
                fig,
                output_path,
                include_plotlyjs=include_plotlyjs,
                config=config,
            )
        return fig
