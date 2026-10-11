"""Standalone interactive visualization helpers for MF6 results.

The matplotlib exporters intentionally pre-render figures and use a small
browser-side slider. The resulting HTML files do not require Jupyter, a Python
kernel, or a web server.
"""

from __future__ import annotations

import base64
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
import pandas as pd
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
    defaults, so override only what you need. Frozen: derive a variant with
    :func:`dataclasses.replace`.

    Parameters
    ----------
    figsize : tuple of (float, float), optional
        Matplotlib figure size in inches, ``(width, height)``. Default ``(10, 8)``.
        :func:`plot_model_head_map` uses it only when it creates the figure (no
        ``ax=``); :func:`export_head_layer_mosaic_slider_html` uses it as the size
        of EACH layer panel when its ``panel_figsize`` is not given.
    cmap : str or matplotlib.colors.Colormap, optional
        Colormap for the head cells, handed to FloPy's ``plot_array``. Default
        ``"viridis"``.
    grid_color : str, optional
        Matplotlib color of the cell-edge overlay. Default ``"#3c4652"``.
    grid_linewidth : float, optional
        Line width of the cell-edge overlay, in points. Default ``0.2``.
    contour_color : str, optional
        Matplotlib color of the head contours. Default ``"white"``.
    contour_linewidth : float, optional
        Line width of the head contours, in points. Default ``0.7``.
    contour_levels : int or sequence of float, optional
        Handed to Matplotlib's contouring as ``levels``: an int is a level count
        (Matplotlib picks up to that many plus one round values), a sequence
        gives the head values to contour, in model length units. Default ``10``.
    show_grid : bool, optional
        Draw the cell edges over the head field. Default ``True``.
    show_contours : bool, optional
        Draw labelled head contours. Default ``True``. Skipped on a frame whose
        layer has fewer than two valid cells or a flat field, since there is
        nothing to contour.
    show_colorbar : bool, optional
        Draw a colorbar labelled ``"Head"``. Default ``True``. The layer mosaic
        draws ONE shared colorbar for the whole frame instead of one per panel.
    dpi : int, optional
        PNG resolution, in dots per inch, of each frame that
        :func:`export_head_map_slider_html` writes, and the mosaic's resolution
        when its ``dpi`` is not given. Default ``140``;
        :func:`export_head_layer_mosaic_slider_html` with no ``style`` uses 200.
        :func:`plot_model_head_map` ignores it -- it returns a figure, not a PNG.
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
    path : pathlib.Path
        The written HTML document.
    labels : tuple of str
        The slider label of each frame, in slider order -- after frame selection
        (``frame_indices``/``frame_stride``/``max_frames``).
    frame_count : int
        Number of frames in the page, which is ``len(labels)``.
    title : str
        The page's ``<title>`` and heading.
    embedded_frames : bool, optional
        ``True`` when every frame is inlined in the HTML as a base64 PNG ``data:``
        URI; ``False`` when the page references PNG files in ``frame_directory``
        and must be moved together with that directory. Default ``True``.
    frame_directory : pathlib.Path or None, optional
        The directory holding the frame PNGs and their ``frames.json`` manifest;
        ``None`` when the frames are embedded. Default ``None``.
    rendered_frames : int, optional
        Frames drawn during this export. Default ``0``.
    reused_frames : int, optional
        Frames taken from an earlier export's PNGs under ``resume=True`` instead of
        being drawn again. Default ``0``.
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
    index : int
        Zero-based position of this frame among the SELECTED frames (after
        ``frame_indices``/``frame_stride``/``max_frames``), not its index in the
        original sequence.
    total : int
        Number of selected frames in the export.
    label : str
        The frame's slider label.
    status : str
        ``"rendered"`` when the frame was drawn in this export, or ``"reused"``
        when ``resume=True`` took an existing PNG instead.
    path : pathlib.Path or None, optional
        The frame's PNG file when frames are written to a directory
        (``embed_frames=False``); ``None`` when the frame is embedded in the HTML.
        Default ``None``.
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


def _checked_labels(labels: Sequence[Any], count: int, *, frames: str) -> list[str]:
    """``labels`` as strings, refused unless there is exactly one per frame.

    The ONE length check every slider exporter shares: a list one too long used
    to be cut short in silence, so every label after a missing one titled the
    wrong frame, and one too short died in an ``IndexError`` naming neither
    count. ``frames`` says where the ``count`` frames came from. A ``str`` is
    refused: it iterates by character, so ``labels="ab"`` passed the count check
    on two frames and labelled them ``"a"`` and ``"b"``.
    """

    if isinstance(labels, (str, bytes)):
        raise TypeError(
            f"labels={labels!r} is one string, which would label the frames one "
            "character each; give one label per frame as a list, e.g. "
            f"labels=[{labels!r}, ...]."
        )
    resolved = [str(label) for label in labels]
    if len(resolved) != int(count):
        raise ValueError(
            f"labels has {len(resolved)} entries but {frames} has {int(count)} "
            f"frames; give one label per frame, before frame_indices/"
            f"frame_stride/max_frames select from them."
        )
    return resolved


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

    Parameters
    ----------
    render_frame : callable
        ``render_frame(value, index) -> Figure`` draws one frame. ``value`` is an
        item of ``frame_values`` and ``index`` is its zero-based position among
        the SELECTED frames (after ``frame_indices``/``frame_stride``/
        ``max_frames``), not its index in ``frame_values``. It may return a
        ``(fig, ...)`` tuple, e.g. ``(fig, ax)``; the first item is the figure
        rasterized. Not called for a frame reused under ``resume=True``.
    frame_values : iterable
        One value per frame, handed to ``render_frame`` in order. Consumed with
        ``list()``, so any iterable works. Must hold at least one item, else
        ``ValueError``.
    output_path : str or pathlib.Path
        The HTML file to write. Its parent directory is created. A
        ``<output_path>.lock`` file is held for the duration of the export.
    labels : sequence of str, optional
        One slider label per item of ``frame_values`` (before selection); each is
        passed through ``str()``. Default ``None`` labels a frame ``str(value)``.
        A length different from ``frame_values`` raises ``ValueError`` naming
        both counts; a single ``str`` raises ``TypeError``.
    title : str, optional
        The page's ``<title>`` and heading (HTML-escaped). Default
        ``"Model results through time"``.
    dpi : int, optional
        Resolution, in dots per inch, each figure is saved at as a PNG. Default
        ``140``.
    interval_ms : int, optional
        Milliseconds between frames while the page's Play button runs. Default
        ``700``.
    close_figures : bool, optional
        Close each figure with ``plt.close`` once it is rasterized, even if
        rasterizing fails. Default ``True``; ``False`` leaves every frame's figure
        open, which holds its memory.
    embed_frames : bool, optional
        ``True`` (default) inlines each frame as a base64 PNG ``data:`` URI, one
        self-contained file. ``False`` writes ``frame_00000.png``, ... plus a
        ``frames.json`` manifest to ``frame_directory`` and references them by
        path relative to the HTML.
    frame_directory : str or pathlib.Path, optional
        Where the PNGs go when ``embed_frames=False``; created if missing.
        Default ``None`` is ``<output stem>_frames`` beside the HTML. Ignored
        when ``embed_frames=True``.
    frame_indices : sequence of int, optional
        Zero-based positions in ``frame_values`` to keep, in the order given
        (so it may reorder or repeat frames). Applied before ``frame_stride``
        and ``max_frames``. An index out of range raises ``IndexError``, a
        negative one included (no counting from the end). Default ``None`` keeps
        every frame.
    frame_stride : int, optional
        Keep every ``frame_stride``-th frame of the selection. Default ``1``
        (all); below 1 raises ``ValueError``.
    max_frames : int, optional
        Keep at most this many frames, after ``frame_indices`` and
        ``frame_stride``. Default ``None`` (no cap); below 1 raises
        ``ValueError``.
    resume : bool, optional
        Reuse a frame whose PNG already exists (non-empty) in ``frame_directory``
        instead of calling ``render_frame`` for it. Requires
        ``embed_frames=False`` (``ValueError`` otherwise), and the existing
        ``frames.json`` must match this export's labels, ``dpi`` and frame count
        (``ValueError`` otherwise). Default ``False`` redraws every frame.
    progress : bool or callable, optional
        ``False`` (default) reports nothing; ``True`` prints one
        ``[i/n] status: label`` line per frame; a callable is called once per
        frame with a :class:`FrameExportProgress`.

    Returns
    -------
    StandaloneHtmlSlider
        A handle to the written document and its render statistics.

    Raises
    ------
    ValueError
        On an empty ``frame_values`` or selection, a ``labels`` length mismatch,
        ``frame_stride``/``max_frames`` below 1, ``resume=True`` with embedded
        frames, or a ``frames.json`` that does not match a resumed export.
    TypeError
        When ``labels`` is a single ``str`` rather than one label per frame.
    IndexError
        When ``frame_indices`` names a frame that does not exist.
    RuntimeError
        When ``<output_path>.lock`` already exists, i.e. another export appears
        to be writing the same file.
    """

    values = list(frame_values)
    if not values:
        raise ValueError("frame_values must contain at least one frame.")
    resolved_labels = (
        [str(value) for value in values]
        if labels is None
        else _checked_labels(labels, len(values), frames="frame_values")
    )
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
    """Per-frame display labels: explicit ``labels``, else the kstpkper, else ``Frame N``.

    Explicit ``labels`` must hold one entry per frame (:func:`_checked_labels`).
    """

    if labels is not None:
        if kstpkpers is None:
            source = "head_frames"
        else:
            source = "the saved head times (kstpkpers)"
        return _checked_labels(labels, len(values), frames=source)
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
    model : SimulationBase
        The myflopy model whose heads are sectioned -- what ``run.model()``
        returns, or a legacy ``SimulationBase``. Its ``gwf`` supplies the grid and,
        unless ``head_frames`` is given, the head file.
    line : LineString, dict, sequence of (x, y), str or pathlib.Path
        The section line, in the grid's coordinates: a shapely ``LineString`` (a
        ``MultiLineString`` is merged, keeping the longest piece if it does not
        join), FloPy's ``{"line": [...]}`` dict, a sequence of ``(x, y)`` vertices,
        or a path to a vector file whose geometries are unioned into one line.
    output_path : str or pathlib.Path
        The HTML file to write.
    kstpkpers : sequence of tuple of (int, int), optional
        The saved output times to render, as FloPy's zero-based
        ``(kstp, kper)`` pairs, read from the model's head file. Default ``None``
        renders every time the head file holds. Ignored when ``head_frames`` is
        given.
    head_frames : sequence of array_like, optional
        Head arrays to render instead of reading the head file, one per frame,
        each ``(nlay, ncpl)`` (a ``(nlay, 1, ncpl)`` array is squeezed) or
        ``(ncpl,)``. Takes precedence over ``kstpkpers``. Default ``None``.
    labels : sequence of str, optional
        One slider label per frame BEFORE selection; each also titles its frame.
        Default ``None`` labels a frame ``"kstpkper=(kstp, kper)"``, or
        ``"Frame N"`` (one-based) for ``head_frames``. A list whose length
        differs from the frame count raises ``ValueError`` naming both counts;
        a single ``str`` raises ``TypeError``.
    style : ModelCrossSectionStyle, optional
        The section's appearance, forwarded to ``plot_model_cross_section``.
        Default ``None`` uses its default style.
    title : str, optional
        The page's ``<title>`` and heading. Default
        ``"Head cross section through time"``. Each frame's own title is its
        label.
    dpi : int, optional
        Resolution, in dots per inch, of each frame's PNG. Default ``140``.
    interval_ms : int, optional
        Milliseconds between frames while the page's Play button runs. Default
        ``700``.
    embed_frames : bool, optional
        Inline every frame in the HTML (default ``True``), or write PNGs to
        ``frame_directory``. See :func:`export_matplotlib_slider_html`.
    frame_directory : str or pathlib.Path, optional
        Where the PNGs go when ``embed_frames=False``. Default ``None`` is
        ``<output stem>_frames`` beside the HTML. Ignored when frames are
        embedded.
    frame_indices : sequence of int, optional
        Zero-based positions among the resolved frames to keep, in the order
        given. Default ``None`` keeps all. See
        :func:`export_matplotlib_slider_html`.
    frame_stride : int, optional
        Keep every ``frame_stride``-th selected frame. Default ``1``.
    max_frames : int, optional
        Keep at most this many frames after the index and stride selection.
        Default ``None`` (no cap).
    resume : bool, optional
        Reuse PNGs an earlier export left in ``frame_directory``; requires
        ``embed_frames=False``. Default ``False``.
    progress : bool or callable, optional
        ``True`` prints one line per frame; a callable receives a
        :class:`FrameExportProgress` per frame. Default ``False``.
    **cross_section_kwargs
        Forwarded to every frame's ``plot_model_cross_section`` call, e.g.
        ``head_layer=``, ``ylim=``, ``xlim=``, ``layer_colors=``,
        ``layer_labels=``, ``show_grid=``, ``show_legend=``. This function sets
        that call's ``title`` (the frame label) and its ``head_data`` (with
        ``head_frames``) or ``kstpkper`` (otherwise) per frame, overriding the
        same name given here. A name ``plot_model_cross_section`` does not accept
        raises ``TypeError`` when the first frame is drawn.

    Returns
    -------
    StandaloneHtmlSlider
        A handle to the written document and its render statistics.

    Raises
    ------
    ValueError
        When ``labels`` does not hold one entry per resolved frame (the
        ``head_frames``, or the saved times), naming both counts; and as
        :func:`export_matplotlib_slider_html` raises for the frame selection.
    TypeError
        When ``labels`` is a single ``str`` rather than one label per frame.
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
    model : SimulationBase
        The myflopy model providing the grid -- what ``run.model()`` returns, or a
        legacy ``SimulationBase``. Its ``gwf.modelgrid`` is drawn.
    head_data : array_like
        The head field, in model length units: ``(nlay, ncpl)``, ``(nlay, 1,
        ncpl)`` (squeezed to ``(nlay, ncpl)``), or one layer's ``(ncpl,)``. A
        non-finite value or a no-value marker (``|value| >= 1e29``, MF6's dry and
        inactive cells) is left blank.
    layer : int, optional
        Zero-based layer to draw from a multi-layer ``head_data``; with a
        ``(ncpl,)`` array the array itself is drawn. Also the layer FloPy's
        ``PlotMapView`` is built for. Default ``0``.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. Default ``None`` creates a new figure of
        ``style.figsize``.
    style : ModelMapStyle, optional
        Colormap, grid and contour overlays and colorbar. Default ``None`` is
        ``ModelMapStyle()``. Its ``dpi`` is not used here.
    title : str, optional
        Axes title. Default ``None`` (or ``""``) sets none.
    vmin : float, optional
        Lower color limit, in head units. Default ``None`` lets Matplotlib take it
        from the drawn layer's valid values.
    vmax : float, optional
        Upper color limit, in head units. Default ``None`` lets Matplotlib take it
        from the drawn layer's valid values.

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
    model : SimulationBase
        The myflopy model providing the grid -- what ``run.model()`` returns, or a
        legacy ``SimulationBase``. Its ``gwf.modelgrid`` is drawn.
    pathlines : pandas.DataFrame, numpy.recarray, or list of them
        Particle tracks, handed to FloPy's ``PlotMapView.plot_pathline``: MF6 PRT
        track data or MODPATH 6/7 pathline data (e.g. FloPy's
        ``PathlineFile.get_alldata()``), all particles in one table or one table
        per particle. Their ``x``/``y`` are in the model's own frame; FloPy
        places them with the grid's offset and rotation.
    layer : int or str, optional
        Which track points to draw. ``"all"`` (default) draws every point. An int
        of 0 or more keeps the points whose layer column equals it, compared as
        given (``k`` for MODPATH data, PRT's ``ilay`` otherwise -- no base
        conversion); a negative int draws every point, as ``"all"`` does. Any
        other string makes FloPy use ``head_layer`` instead.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. Default ``None`` creates a new 10 x 8 inch figure.
    head_data : array_like, optional
        A head field shaded (alpha 0.35, FloPy's default ``viridis``) under the
        tracks: ``(nlay, ncpl)``, of which layer ``head_layer`` is drawn, or one
        layer's ``(ncpl,)``. FloPy masks NaN and MF6's ``+/-1e30`` no-value
        markers. Default ``None`` draws no head field.
    head_layer : int, optional
        Zero-based layer of ``head_data`` to shade, and the layer FloPy's
        ``PlotMapView`` is built for. Default ``0``.
    pathline_color : str, optional
        Matplotlib color of every track. Default ``None`` leaves FloPy's grey
        (``"0.5"``).
    pathline_alpha : float, optional
        Track opacity, 0 to 1. Default ``0.7``.
    pathline_linewidth : float, optional
        Track line width, in points. Default ``1.0``.
    show_grid : bool, optional
        Draw the cell edges. Default ``True``.
    grid_alpha : float, optional
        Opacity of the cell edges, 0 to 1. Default ``0.25``.
    title : str, optional
        Axes title. Default ``"Particle pathlines"``.

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
    model : SimulationBase
        The myflopy model whose heads are mapped -- what ``run.model()`` returns,
        or a legacy ``SimulationBase``. Its ``gwf`` supplies the grid and, unless
        ``head_frames`` is given, the head file.
    output_path : str or pathlib.Path
        The HTML file to write.
    layer : int, optional
        Zero-based layer to map. Default ``0``. With ``(ncpl,)`` head frames the
        frame itself is drawn.
    kstpkpers : sequence of tuple of (int, int), optional
        The saved output times to render, as FloPy's zero-based
        ``(kstp, kper)`` pairs, read from the model's head file. Default ``None``
        renders every time the head file holds. Ignored when ``head_frames`` is
        given.
    head_frames : sequence of array_like, optional
        Head arrays to render instead of reading the head file, one per frame,
        each ``(nlay, ncpl)`` (a ``(nlay, 1, ncpl)`` array is squeezed) or
        ``(ncpl,)``. Takes precedence over ``kstpkpers``. Default ``None``.
    labels : sequence of str, optional
        One slider label per frame BEFORE selection; each also titles its frame.
        Default ``None`` labels a frame ``"kstpkper=(kstp, kper)"``, or
        ``"Frame N"`` (one-based) for ``head_frames``. A list whose length
        differs from the frame count raises ``ValueError`` naming both counts;
        a single ``str`` raises ``TypeError``.
    style : ModelMapStyle, optional
        Map appearance, and the PNG resolution (``style.dpi``) -- this function
        has no ``dpi`` of its own. Default ``None`` is ``ModelMapStyle()``.
    title : str, optional
        The page's ``<title>`` and heading. Default ``"Head map through time"``.
    vmin : float, optional
        Lower color limit for every frame, in head units. Default ``None``: taken
        from the frames when ``shared_color_scale``, else per frame.
    vmax : float, optional
        Upper color limit for every frame, in head units. Default ``None``: taken
        from the frames when ``shared_color_scale``, else per frame.
    shared_color_scale : bool, optional
        When ``True`` (default), a ``vmin``/``vmax`` left as ``None`` is filled
        with the minimum/maximum head of ``layer`` over ALL selected frames
        (ignoring non-finite values and ``|value| >= 1e29``), so colors compare
        across frames; from the head file that costs one extra read per frame.
        ``False`` lets each frame scale to its own values. An explicit
        ``vmin``/``vmax`` always wins.
    embed_frames : bool, optional
        Inline every frame in the HTML (default ``True``), or write PNGs to
        ``frame_directory``. See :func:`export_matplotlib_slider_html`.
    frame_directory : str or pathlib.Path, optional
        Where the PNGs go when ``embed_frames=False``. Default ``None`` is
        ``<output stem>_frames`` beside the HTML. Ignored when frames are
        embedded.
    frame_indices : sequence of int, optional
        Zero-based positions among the resolved frames to keep, in the order
        given. Default ``None`` keeps all. See
        :func:`export_matplotlib_slider_html`.
    frame_stride : int, optional
        Keep every ``frame_stride``-th selected frame. Default ``1``.
    max_frames : int, optional
        Keep at most this many frames after the index and stride selection.
        Default ``None`` (no cap).
    resume : bool, optional
        Reuse PNGs an earlier export left in ``frame_directory``; requires
        ``embed_frames=False``. Default ``False``.
    progress : bool or callable, optional
        ``True`` prints one line per frame; a callable receives a
        :class:`FrameExportProgress` per frame. Default ``False``.

    Returns
    -------
    StandaloneHtmlSlider
        A handle to the written document and its render statistics.

    Raises
    ------
    ValueError
        When ``labels`` does not hold one entry per resolved frame (the
        ``head_frames``, or the saved times), naming both counts; and as
        :func:`export_matplotlib_slider_html` raises for the frame selection.
    TypeError
        When ``labels`` is a single ``str`` rather than one label per frame.
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

    Each frame is one figure holding a :func:`plot_model_head_map` panel per
    layer (titled ``Layer 1``, ``Layer 2``, ... one-based), laid out ``ncols``
    wide, with at most one colorbar for the whole figure (``style.show_colorbar``),
    and the frames are packed into a standalone HTML slider like
    :func:`export_head_map_slider_html`. ``dpi`` controls the PNG resolution
    directly. ``panel_figsize`` controls the physical size of each layer panel
    before rasterization.

    Parameters
    ----------
    model : SimulationBase
        The myflopy model whose heads are mapped -- what ``run.model()`` returns,
        or a legacy ``SimulationBase``. Its ``gwf`` supplies the grid, the layer
        count and, unless ``head_frames`` is given, the head file.
    output_path : str or pathlib.Path
        The HTML file to write.
    layers : sequence of int, optional
        Zero-based layers to draw, one panel each, in the order given. Default
        ``None`` draws every layer of the grid. An empty sequence raises
        ``ValueError``.
    kstpkpers : sequence of tuple of (int, int), optional
        The saved output times to render, as FloPy's zero-based
        ``(kstp, kper)`` pairs, read from the model's head file. Default ``None``
        renders every time the head file holds. Ignored when ``head_frames`` is
        given.
    head_frames : sequence of array_like, optional
        Head arrays to render instead of reading the head file, one per frame,
        each ``(nlay, ncpl)`` (a ``(nlay, 1, ncpl)`` array is squeezed). Takes
        precedence over ``kstpkpers``. Default ``None``.
    labels : sequence of str, optional
        One slider label per frame BEFORE selection; each also becomes that
        frame's figure title. Default ``None`` labels a frame
        ``"kstpkper=(kstp, kper)"``, or ``"Frame N"`` (one-based) for
        ``head_frames``. A list whose length differs from the frame count
        raises ``ValueError`` naming both counts; a single ``str`` raises
        ``TypeError``.
    style : ModelMapStyle, optional
        Each panel's appearance. Its colorbar is drawn once for the whole frame
        (``style.show_colorbar``), never per panel. Default ``None`` is
        ``ModelMapStyle(dpi=200)``.
    title : str, optional
        The page's ``<title>`` and heading. Default
        ``"Head layer mosaic through time"``.
    ncols : int, optional
        Panels per row; the row count is ``ceil(len(layers) / ncols)`` and
        unused grid cells are hidden. Default ``3``.
    dpi : int, optional
        PNG resolution, in dots per inch. Default ``None`` uses ``style.dpi``
        (200 with the default style). Zero or less raises ``ValueError``.
    panel_figsize : tuple of (float, float), optional
        Size of ONE panel in inches, ``(width, height)``; the figure is
        ``ncols`` panels wide and one panel high per row. Default ``None`` uses
        ``style.figsize``. A length other than two, or a value of zero or less,
        raises ``ValueError``.
    vmin : float, optional
        Lower color limit for every panel of every frame, in head units. Default
        ``None``: taken from the frames when ``shared_color_scale``, else per
        panel.
    vmax : float, optional
        Upper color limit for every panel of every frame, in head units. Default
        ``None``: taken from the frames when ``shared_color_scale``, else per
        panel.
    shared_color_scale : bool, optional
        When ``True`` (default), a ``vmin``/``vmax`` left as ``None`` is filled
        with the minimum/maximum head over ``layers`` and ALL selected frames
        (ignoring non-finite values and ``|value| >= 1e29``), so every panel
        shares one scale; from the head file that costs one extra read per
        frame. ``False`` lets each panel scale to its own values, which makes the
        single colorbar describe only the first panel. An explicit
        ``vmin``/``vmax`` always wins.
    embed_frames : bool, optional
        Inline every frame in the HTML (default ``True``), or write PNGs to
        ``frame_directory``. See :func:`export_matplotlib_slider_html`.
    frame_directory : str or pathlib.Path, optional
        Where the PNGs go when ``embed_frames=False``. Default ``None`` is
        ``<output stem>_frames`` beside the HTML. Ignored when frames are
        embedded.
    frame_indices : sequence of int, optional
        Zero-based positions among the resolved frames to keep, in the order
        given. Default ``None`` keeps all. See
        :func:`export_matplotlib_slider_html`.
    frame_stride : int, optional
        Keep every ``frame_stride``-th selected frame. Default ``1``.
    max_frames : int, optional
        Keep at most this many frames after the index and stride selection.
        Default ``None`` (no cap).
    resume : bool, optional
        Reuse PNGs an earlier export left in ``frame_directory``; requires
        ``embed_frames=False``. Default ``False``.
    progress : bool or callable, optional
        ``True`` prints one line per frame; a callable receives a
        :class:`FrameExportProgress` per frame. Default ``False``.

    Returns
    -------
    StandaloneHtmlSlider
        A handle to the written document and its render statistics.

    Raises
    ------
    ValueError
        When ``labels`` does not hold one entry per resolved frame (the
        ``head_frames``, or the saved times), naming both counts; and as
        :func:`export_matplotlib_slider_html` raises for the frame selection.
    TypeError
        When ``labels`` is a single ``str`` rather than one label per frame.
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


def _tracks_in_world_frame(pathlines, modelgrid):
    """A copy of PRT ``pathlines`` in the world frame the exported grid mesh is in.

    MF6 writes tracks in the model's OWN frame. On an offset grid the grid mesh sat
    at x 1,340,000 and the pathline mesh at x 5..20. Accepts what
    ``Vtk.add_pathline_points`` does: a DataFrame, a record array, or a list of them.
    """

    if isinstance(pathlines, list):
        return [_tracks_in_world_frame(part, modelgrid) for part in pathlines]
    if isinstance(pathlines, pd.DataFrame):
        frame = pathlines.copy()
        frame["x"], frame["y"] = modelgrid.get_coords(frame["x"].to_numpy(), frame["y"].to_numpy())
        return frame
    records = np.array(pathlines, copy=True)
    records["x"], records["y"] = modelgrid.get_coords(records["x"], records["y"])
    return records


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
    model : SimulationBase
        The myflopy model providing the grid -- what ``run.model()`` returns, or a
        legacy ``SimulationBase``. Its ``gwf`` is exported to VTK with FloPy's
        ``Vtk.add_model``.
    pathlines : pandas.DataFrame, numpy.recarray, or list of them
        Particle tracks, handed to FloPy's ``Vtk.add_pathline_points``: MF6 PRT
        track data or MODPATH pathline data, all particles in one table or one
        table per particle; each needs ``x`` and ``y`` fields. Their ``x``/``y``
        are read in the model's own frame and moved to the grid's world frame
        (its offset and rotation) first, so tracks and grid line up on an
        offset grid. The input is copied, not modified.
    vertical_exaggeration : float, optional
        Vertical scale factor for the exported points (FloPy's
        ``vertical_exageration``). Default ``1.0``, no exaggeration.
    model_style : str, optional
        PyVista ``style`` for the grid mesh: ``"wireframe"`` (default; PyVista
        draws the wireframe of the OUTER geometry only), ``"surface"``,
        ``"points"`` or ``"points_gaussian"``.
    model_opacity : float, optional
        Grid mesh opacity, 0 to 1. Default ``0.25``.
    pathline_cmap : str, optional
        Colormap for the track mesh, coloured by its ``time`` (or ``t``) array
        when it has one. Default ``"viridis"``.
    pathline_width : float, optional
        PyVista ``line_width`` of the tracks, which are drawn as tubes. Default
        ``4.0``.
    show_edges : bool, optional
        Draw the grid mesh's cell edges. Default ``True``.
    off_screen : bool, optional
        Build the PyVista ``Plotter`` off screen, opening no window -- what a
        headless export needs. Default ``True``.

    Returns
    -------
    VtkScene
        The configured plotter + meshes, as a Picture.

    Raises
    ------
    ModuleNotFoundError
        When ``pyvista`` is not installed.
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
    vtk.add_pathline_points(_tracks_in_world_frame(pathlines, model.gwf.modelgrid))
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
