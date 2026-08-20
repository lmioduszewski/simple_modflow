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
from collections.abc import Callable, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
from flopy.plot import PlotMapView

from myflopy._optional import require
from myflopy.modflow.mf6.cross_section_plotting import (
    ModelCrossSectionStyle,
    plot_model_cross_section,
)
from myflopy.viz import Picture, VtkScene, _normalize_frames, mpl_axes

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from myflopy.modflow.mf6.simulation.base import SimulationBase


_DEFAULT_MOSAIC_DPI = 200


@dataclass(frozen=True)
class ModelMapStyle:
    """Appearance settings for the matplotlib head-map slider exporters.

    A single bundle of styling knobs passed to the standalone HTML head-map
    exporters (e.g. :func:`export_head_map_slider_html`) so every rendered frame
    looks consistent: figure size/DPI, colormap, grid and contour styling, and
    which overlays (grid, contours, colorbar) to draw. All fields have sensible
    defaults, so override only what you need.

    Attributes
    ----------
    figsize, dpi
        Matplotlib figure size and resolution for each frame.
    cmap
        Colormap for the head field.
    grid_color, grid_linewidth, show_grid
        Cell-edge overlay styling and toggle.
    contour_color, contour_linewidth, contour_levels, show_contours
        Head-contour overlay styling, level count/values, and toggle.
    show_colorbar
        Whether to draw the colorbar.
    """

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
    """A handle to a generated standalone frame-slider HTML file and its stats.

    The return value of the ``export_*_slider_html`` functions: a self-contained
    HTML document with a time/layer slider that pages through pre-rendered frames
    (heads, cross-sections, ...), shareable without a Python kernel. This object
    records where it was written (``path``), the per-frame ``labels`` and
    ``frame_count``, whether frames are inlined as data URIs (``embedded_frames``)
    or stored alongside in ``frame_directory``, and how many frames were freshly
    rendered vs reused from cache.

    Attributes
    ----------
    path
        The written HTML document.
    labels, frame_count
        Per-frame slider labels and the number of frames.
    title
        Document title.
    embedded_frames, frame_directory
        Whether frames are inlined, else the directory holding the frame images.
    rendered_frames, reused_frames
        Counts of newly rendered vs cache-reused frames.
    """

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
    """One progress update emitted while a slider's frames are being rendered.

    Passed to the optional ``progress`` callback of the ``export_*_slider_html``
    functions once per frame, so callers can show a progress bar or log. Reports
    which frame (``index`` of ``total``), its ``label``, the ``status`` (e.g.
    rendered vs reused-from-cache), and the frame's output ``path`` when written.

    Attributes
    ----------
    index, total
        Zero-based frame index and total frame count.
    label
        The frame's slider label.
    status
        Short status string for this frame.
    path
        Output path for the frame, when applicable.
    """

    index: int
    total: int
    label: str
    status: str
    path: Path | None = None


def _figure_to_data_uri(fig: Figure, *, dpi: int) -> str:
    """Render a matplotlib figure to a base64 PNG ``data:`` URI for inline embedding."""

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def _write_figure_png(fig: Figure, path: Path, *, dpi: int) -> Path:
    """Save a figure to ``path`` as PNG via an atomic temp-file replace."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    fig.savefig(temporary, format="png", dpi=dpi, bbox_inches="tight")
    temporary.replace(path)
    return path


def _write_text_atomic(path: Path, text: str) -> None:
    """Write ``text`` to ``path`` atomically via a temp-file replace."""

    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


@contextmanager
def _export_lock(output_path: Path):
    """An exclusive-file-lock context so two exports cannot write the same output at once."""

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
    """Report a frame-export progress event via a callback, or a printed line if ``progress`` is truthy."""

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
    """Subselect ``(values, labels)`` by explicit indices, stride, and a max-frame cap (validated)."""

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
    """Build a standalone HTML page embedding image ``frames`` with a browser-side play/slider."""

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
    render_frame: Callable[[Any, int], Figure | tuple[Figure, Any]],
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
    model: SimulationBase,
    *,
    kstpkpers: Sequence[tuple[int, int]] | None,
    head_frames: Sequence[Any] | None,
) -> tuple[list[Any], list[tuple[int, int]] | None]:
    """Resolve the animation frames: explicit ``head_frames`` arrays, else the model's saved kstpkper.

    Returns ``(values, resolved_kstpkpers)`` where the second is ``None`` when
    frames are supplied directly.
    """

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
    model: SimulationBase,
    values: Sequence[Any],
    *,
    kstpkpers: Sequence[tuple[int, int]] | None,
    labels: Sequence[str] | None,
) -> list[str]:
    """Per-frame display labels: explicit ``labels``, else the kstpkper, else ``Frame N``."""

    if labels is not None:
        return [str(label) for label in labels]
    if kstpkpers is None:
        return [f"Frame {index + 1}" for index in range(len(values))]
    return [f"kstpkper={tuple(value)}" for value in kstpkpers]


def _normalize_head_frame(values):
    """Squeeze a single-row middle axis so a head array is ``(nlay, ncpl)``."""

    array = np.asarray(values)
    if array.ndim == 3 and array.shape[1] == 1:
        return array[:, 0, :]
    return array


def _head_frame(model: SimulationBase, value, resolved_kstpkpers):
    """One normalized head array: ``value`` itself, or read from the model at that kstpkper."""

    data = value if resolved_kstpkpers is None else model.gwf.output.head().get_data(kstpkper=value)
    return _normalize_head_frame(data)


def _clean_head_values(values):
    """Mask non-finite and MODFLOW dry/no-data sentinels (``|value| >= 1e29``) in a head array."""

    array = np.asarray(values, dtype=float)
    return np.ma.masked_where(~np.isfinite(array) | (np.abs(array) >= 1.0e29), array)


def _shared_head_limits(
    model: SimulationBase,
    values: Sequence[Any],
    resolved_kstpkpers,
    layers: Sequence[int],
) -> tuple[float | None, float | None]:
    """The shared (min, max) head across all frames and ``layers``, for a fixed color scale."""

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
    model: SimulationBase,
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
    """Write a self-contained HTML slider of head cross-sections through time.

    Renders the head field along ``line`` for each selected step (reusing the
    ``plot_model_cross_section`` styling) and packs the frames into one standalone
    HTML document with a time slider -- shareable without a Python kernel. Frame
    selection (``frame_indices``/``frame_stride``/``max_frames``), caching
    (``resume``/``frame_directory``), and ``progress`` callbacks are supported.

    Parameters
    ----------
    model
        The model whose heads are sectioned.
    line
        The cross-section line (a 2-point LineString or coordinate pair).
    output_path
        Destination HTML file.
    kstpkpers, head_frames, labels
        Which time steps (or explicit head frames) to render and their labels.
    style, title, dpi, interval_ms
        Cross-section styling and animation/output options.
    embed_frames, frame_directory, resume
        Whether to inline frames; else where to store them and whether to reuse.
    progress
        ``True`` to print progress, or a :class:`FrameExportProgress` callback.
    **cross_section_kwargs
        Extra keyword arguments forwarded to ``plot_model_cross_section``.

    Returns
    -------
    StandaloneHtmlSlider
        A handle to the written document and its render statistics.
    """

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
        """Render one cross-section frame figure for the slider export."""

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
    model: SimulationBase,
    head_data,
    *,
    layer: int = 0,
    ax=None,
    style: ModelMapStyle | None = None,
    title: str | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
):
    """Draw one head map (plan view) for a layer onto a matplotlib axis.

    Renders a single head field as a filled cell map with optional grid and
    contour overlays, using a :class:`ModelMapStyle`. The static, single-frame
    building block beneath :func:`export_head_map_slider_html`; call it directly
    to drop a head map onto your own figure (pass ``ax`` to compose).

    Parameters
    ----------
    model
        The model providing the grid (``model.gwf.modelgrid``).
    head_data
        A head array -- per-layer ``(nlay, ncpl)`` or a single layer ``(ncpl,)``.
    layer
        Layer index to draw when ``head_data`` is 2-D.
    ax
        Existing axis to draw on; a new figure is created when ``None``.
    style
        Appearance settings; defaults to :class:`ModelMapStyle`.
    title, vmin, vmax
        Title and fixed color limits.

    Returns
    -------
    tuple
        ``(fig, ax)`` for further composition.
    """

    style = ModelMapStyle() if style is None else style
    if ax is None:
        fig, ax = mpl_axes(figsize=style.figsize)
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
    model: SimulationBase,
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
    """Draw particle pathlines (plan view) on a matplotlib map, over an optional head field.

    Plots 2-D pathlines from a PRT/MODPATH-style result on top of the model grid,
    optionally shaded by a head field. Accepts pathlines from any
    FloPy-compatible source. The static map counterpart to the 3-D
    :func:`build_particle_tracking_scene`.

    Parameters
    ----------
    model
        The model providing the grid.
    pathlines
        Particle pathline records (PRT/MODPATH-compatible).
    layer
        Which layer's pathlines to draw, or ``"all"``.
    ax
        Existing axis to draw on; a new figure is created when ``None``.
    head_data, head_layer
        Optional head field to shade beneath the pathlines, and its layer.
    pathline_color, pathline_alpha, pathline_linewidth
        Pathline styling.
    show_grid, grid_alpha, title
        Grid overlay toggle/opacity and plot title.

    Returns
    -------
    tuple
        ``(fig, ax)`` for further composition.
    """

    if ax is None:
        fig, ax = mpl_axes(figsize=(10, 8))
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
    model: SimulationBase,
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
    """Write a self-contained HTML slider of plan-view head maps through time.

    Renders one :func:`plot_model_head_map` per selected step for a given layer and
    packs them into a standalone HTML document with a time slider -- shareable
    without a Python kernel. By default a shared color scale is computed across all
    frames so the colormap is comparable step to step. Supports frame
    selection/caching and ``progress`` callbacks like the other slider exporters.

    Parameters
    ----------
    model
        The model whose heads are mapped.
    output_path
        Destination HTML file.
    layer
        Layer index to map.
    kstpkpers, head_frames, labels
        Which time steps (or explicit head frames) to render and their labels.
    style, title
        Map appearance (:class:`ModelMapStyle`) and document title.
    vmin, vmax, shared_color_scale
        Fixed color limits, or auto-share one scale across frames.
    embed_frames, frame_directory, frame_indices, frame_stride, max_frames, resume
        Frame inlining, caching, and selection controls.
    progress
        ``True`` to print progress, or a :class:`FrameExportProgress` callback.

    Returns
    -------
    StandaloneHtmlSlider
        A handle to the written document and its render statistics.
    """

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
        """Render one head-map frame figure for the slider export."""

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
    model: SimulationBase,
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
        """Render one multi-layer head-mosaic frame figure for the slider export."""

        data = _head_frame(model, value, resolved_kstpkpers)
        nrows = int(np.ceil(len(layers) / ncols))
        fig, axes = mpl_axes(
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
        for ax, layer in zip(axes, layers, strict=False):
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
    model: SimulationBase,
    pathlines,
    *,
    vertical_exaggeration: float = 1.0,
    model_style: str = "wireframe",
    model_opacity: float = 0.25,
    pathline_cmap: str = "viridis",
    pathline_width: float = 4.0,
    show_edges: bool = True,
    off_screen: bool = True,
) -> VtkScene:
    """Build an interactive 3-D PyVista scene of particle pathlines over the grid.

    Exports the model grid and ``pathlines`` to VTK, renders the grid as a
    translucent wireframe with the pathlines drawn as time-colored tubes, and
    returns a :class:`~myflopy.viz.VtkScene` -- display it, ``.show()`` it, or
    ``.html(path)`` it like any other picture. Requires ``pyvista``.

    Parameters
    ----------
    model
        The flow model providing the grid.
    pathlines
        Particle pathline records to render.
    vertical_exaggeration
        Vertical scale factor for the 3-D view.
    model_style, model_opacity, show_edges
        Grid rendering style/opacity and cell-edge toggle.
    pathline_cmap, pathline_width
        Colormap and tube width for the time-colored pathlines.
    off_screen
        Render without opening a window (required for headless export).

    Returns
    -------
    VtkScene
        The configured plotter + meshes, as a Picture.
    """

    pv = require("pyvista", feature="interactive 3-D particle-tracking scenes")
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
    return VtkScene(plotter, meshes=meshes, title="particle pathlines")


def _select_plotly_frames(
    fig,
    *,
    frame_indices: Sequence[int] | None = None,
    frame_stride: int = 1,
    max_frames: int | None = None,
):
    """Subselect a Plotly figure's animation frames by indices/stride/max (no-op if unfiltered)."""

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


def _mpl_figure_of(picture):
    """A Matplotlib ``Figure`` for ``picture``, whatever renderer it is.

    The bridge that makes a PNG animation work across renderers. Deliberately
    NOT via Plotly's static export: ``kaleido`` needs a Chrome install that this
    package does not depend on, and every plotly picture here already carries a
    ``plot_mpl()``. Those return a Figure, an Axes, or a ``(fig, ax)`` pair
    depending on the class, so the shapes are normalized here rather than at
    four call sites.
    """

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    if isinstance(picture, Figure):
        return picture
    axes = getattr(picture, "axes", None)
    if isinstance(axes, Axes):                      # viz.MplPicture
        return axes.figure
    if hasattr(picture, "plot_mpl"):
        drawn = picture.plot_mpl()
        if isinstance(drawn, Figure):
            return drawn
        if isinstance(drawn, Axes):
            return drawn.figure
        if isinstance(drawn, (tuple, list)) and drawn:
            first = drawn[0]
            return first if isinstance(first, Figure) else first.figure
    raise TypeError(
        f"{type(picture).__name__} cannot be rendered to a frame image: it has "
        "no `plot_mpl()` and is not a Matplotlib picture. Use "
        'backend="plotly" for a live plotly animation.'
    )


class SliderAnimation(Picture):
    """Frames pre-rendered to images, paged by a browser slider (plan 8.6a).

    The general form: each frame is rasterized independently, so the frames need
    share nothing at all -- a map, a cross-section and a 3-D scene can sit in one
    animation. It is also the form whose size does not grow with cell count,
    which is what makes it the practical choice for a large grid.

    There is no single figure behind it, so ``.fig`` raises; the picture IS the
    exported page. ``.html(path)`` writes it and returns the path, and
    ``.export(path, ...)`` returns the richer
    :class:`StandaloneHtmlSlider` handle when you want the frame manifest or the
    resume/progress machinery.
    """

    def __init__(self, frames, *, title=None, dpi: int = 140, interval_ms: int = 700):
        """Hold ``[(label, picture), ...]`` to be rasterized on export."""

        self.frames = _normalize_frames(frames)
        self.title = title or "Model results through time"
        self.dpi = dpi
        self.interval_ms = interval_ms

    def __repr__(self):
        """Name the frame count, which is what you check when one is missing."""

        return f"SliderAnimation({len(self.frames)} frames)"

    @property
    def labels(self) -> list[str]:
        """The per-frame slider labels, in order."""

        return [label for label, _ in self.frames]

    @property
    def fig(self):
        """Not available: the frames are images, not one figure."""

        raise TypeError(
            "A SliderAnimation has no single figure -- its frames are rendered "
            "images. Use `.html(path)` to write the page, `.frames` for the "
            '(label, picture) pairs, or backend="plotly" for a live figure.'
        )

    def export(self, path, **kwargs) -> StandaloneHtmlSlider:
        """Write the slider page and return the full :class:`StandaloneHtmlSlider`.

        ``.html(path)`` returns a ``Path`` like every other picture; this is the
        same write with the richer handle -- frame counts, the frame directory,
        and how many frames were rendered vs reused. ``resume=``, ``progress=``,
        ``embed_frames=`` and the rest ride through to
        :func:`export_matplotlib_slider_html`.
        """

        pictures = [picture for _, picture in self.frames]
        kwargs.setdefault("dpi", self.dpi)
        kwargs.setdefault("interval_ms", self.interval_ms)
        return export_matplotlib_slider_html(
            lambda picture, _index: _mpl_figure_of(picture),
            pictures,
            path,
            labels=self.labels,
            title=self.title,
            **kwargs,
        )

    def html(self, path, **kwargs):
        """Write the standalone slider page and return its path."""

        return self.export(path, **kwargs).path

    def save(self, path, **kwargs):
        """Write the animation out; only ``.html`` can hold every frame."""

        from pathlib import Path as _Path

        path = _Path(path)
        if path.suffix.lower() in {".html", ".htm"}:
            return self.html(path, **kwargs)
        raise ValueError(
            f"An animation cannot be saved as {path.suffix!r} -- a still image "
            "holds one frame. Save to .html, or save a single frame's picture."
        )

    def show(self, *args, **kwargs):
        """Display the animation inline (Jupyter) by embedding the slider page."""

        import tempfile

        from IPython.display import HTML, display

        with tempfile.TemporaryDirectory() as tmp:
            written = self.html(Path(tmp) / "animation.html")
            return display(HTML(written.read_text(encoding="utf-8")))

    def _repr_mimebundle_(self, *args, **kwargs):
        """Render inline in Jupyter as the self-contained slider page."""

        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            written = self.html(Path(tmp) / "animation.html")
            return {"text/html": written.read_text(encoding="utf-8")}
