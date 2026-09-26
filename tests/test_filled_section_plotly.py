"""`section(fill=...)` on Plotly draws what Matplotlib draws (2026-09-26).

The filled section -- the cells a line cuts, coloured by layer or by a field, with
the water levels on top -- was Matplotlib-only from ledger 177 on, and `fill=` with
the default backend raised. It is now one picture, `FilledSection`, with two
renderers. The Plotly one does not re-derive the geometry: both read one FloPy
``PlotCrossSection`` and one set of resolved arrays, colours and labels.

"Equivalent" is asserted, not assumed. Each test draws the SAME call on both
backends, reads back what was actually drawn -- every filled cell's polygon and
colour, the separate cell outlines, every water-level segment and its colour, the
legend entries in order, and both axis ranges -- and compares the two. Matplotlib
leaves a reversed section's data alone and inverts the VIEW, so its x values are
read as distance along the view (``lo + hi - x``) before comparing.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402
import shapely as shp  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402

from myflopy.modflow.mf6.cross_section_plotting import FilledSection  # noqa: E402
from myflopy.viz import mpl_color  # noqa: E402

DECIMALS = 6


def _mid_line(model):
    """A West-East line across the middle of the model's grid."""

    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    middle = (ymin + ymax) / 2
    return shp.LineString([(xmin + 1, middle), (xmax - 1, middle)])


def _polygon(xs, ys):
    """A polygon as a hashable key: closing vertex dropped, rounded."""

    points = np.column_stack([np.asarray(xs, float), np.asarray(ys, float)])
    if len(points) > 1 and np.allclose(points[0], points[-1]):
        points = points[:-1]
    return tuple(map(tuple, np.round(points, DECIMALS)))


def _colour(value):
    """Any colour spelling as ``(r, g, b, a)``: 0-255 ints and a rounded alpha."""

    red, green, blue, alpha = mpl_color(value) if isinstance(value, str) else value
    return (round(red * 255), round(green * 255), round(blue * 255), round(alpha, 3))


def _split(xs, ys):
    """Pieces of a ``None``-separated Plotly trace."""

    pieces, current_x, current_y = [], [], []
    for x, y in zip(list(xs or ()), list(ys or ()), strict=True):
        if x is None or y is None:
            if current_x:
                pieces.append((current_x, current_y))
            current_x, current_y = [], []
        else:
            current_x.append(float(x))
            current_y.append(float(y))
    if current_x:
        pieces.append((current_x, current_y))
    return pieces


def _mpl_drawing(figure):
    """What the Matplotlib renderer drew, in along-the-view coordinates."""

    axes = figure.axes[0]
    figure.canvas.draw()                    # a collection maps its colours at draw time
    low, high = sorted(axes.get_xlim())
    inverted = axes.xaxis_inverted()

    def along(x):
        x = np.asarray(x, dtype=float)
        return low + high - x if inverted else x

    cells, outlines = [], []
    for collection in axes.collections:
        paths = collection.get_paths()
        if collection.get_array() is not None:                  # the fill
            faces = collection.get_facecolors()
            for index, path in enumerate(paths):
                face = faces[index] if len(faces) > 1 else faces[0]
                cells.append((_polygon(along(path.vertices[:, 0]), path.vertices[:, 1]),
                              _colour(face)))
        else:                                                    # plot_grid's outlines
            outlines += [_polygon(along(p.vertices[:, 0]), p.vertices[:, 1]) for p in paths]
    segments = {}
    for line in axes.lines:
        x = along(line.get_xdata())
        y = np.asarray(line.get_ydata(), dtype=float)
        segments.setdefault(_colour(to_rgba(line.get_color())), []).append(
            (round(float(x.min()), DECIMALS), round(float(x.max()), DECIMALS),
             round(float(y[0]), DECIMALS))
        )
    legend = axes.get_legend()
    return {
        "cells": sorted(cells),
        "outlines": sorted(outlines),
        "segments": {colour: sorted(v) for colour, v in segments.items()},
        "legend": [text.get_text() for text in legend.get_texts()] if legend else [],
        "x": (low, high),
        "y": tuple(sorted(axes.get_ylim())),
        "colorbar": len(figure.axes) > 1,
    }


def _plotly_drawing(picture):
    """What the Plotly renderer drew, from the traces' own tags."""

    figure = picture.fig
    cells, outlines, segments, legend = [], [], {}, []
    for trace in figure.data:
        if trace.meta == "cells":
            cells += [(_polygon(x, y), _colour(trace.fillcolor)) for x, y in _split(trace.x, trace.y)]
        elif trace.meta == "grid":
            outlines += [_polygon(x, y) for x, y in _split(trace.x, trace.y)]
        elif trace.meta == "surface":
            segments.setdefault(_colour(trace.line.color), []).extend(
                (round(min(x), DECIMALS), round(max(x), DECIMALS), round(y[0], DECIMALS))
                for x, y in _split(trace.x, trace.y)
            )
        if trace.meta in ("cells", "surface") and trace.showlegend:
            legend.append(trace.name)
    return {
        "cells": sorted(cells),
        "outlines": sorted(outlines),
        "segments": {colour: sorted(v) for colour, v in segments.items()},
        "legend": legend if figure.layout.showlegend else [],
        "x": tuple(figure.layout.xaxis.range),
        "y": tuple(figure.layout.yaxis.range),
        "colorbar": any(t.meta == "colorbar" for t in figure.data),
    }


def _assert_same_drawing(mpl, plotly):
    assert mpl["cells"], "the Matplotlib renderer drew no cells -- nothing is being compared"
    assert plotly["cells"] == mpl["cells"], "the two renderers filled different cells or colours"
    assert plotly["outlines"] == mpl["outlines"], "the two renderers outlined different cells"
    assert plotly["segments"] == mpl["segments"], "the water levels differ"
    assert plotly["legend"] == mpl["legend"]
    assert plotly["colorbar"] == mpl["colorbar"]
    np.testing.assert_allclose(plotly["x"], mpl["x"], rtol=0, atol=1e-6)
    np.testing.assert_allclose(plotly["y"], mpl["y"], rtol=0, atol=1e-6)


#: One per option that changes WHAT is drawn. `fill=<array>` is added in the
#: test itself, because the array needs the model's shape.
SCENARIOS = {
    "layer": {"fill": "layer"},
    "layer-labels-and-two-water-levels": {
        "fill": "layer", "head_layers": [0, 2], "layer_labels": ["a", "b", "c", "d"],
    },
    "layer-subset-crops": {"fill": "layer", "layers": [1, 2]},
    "layer-reversed": {"fill": "layer", "reverse": True},
    "layer-subset-reversed": {"fill": "layer", "layers": [0], "reverse": True},
    "layer-no-grid-no-water": {"fill": "layer", "show_grid": False, "head_layers": None},
    "layer-no-legend": {"fill": "layer", "legend": False},
    "results": {"fill": "results", "per": 1, "fill_label": "head (ft)"},
    "results-subset": {"fill": "results", "layers": [0]},
    "results-other-colormap": {"fill": "results", "fill_cmap": "magma"},
}


@pytest.mark.parametrize("options", SCENARIOS.values(), ids=SCENARIOS.keys())
def test_both_renderers_draw_the_same_section(options, canonical_run):
    """Cell for cell: polygons, colours, outlines, water levels, legend, ranges."""

    import matplotlib.pyplot as plt

    line = _mid_line(canonical_run)
    options = dict(options)
    if "layer_labels" in options:
        options["layer_labels"] = options["layer_labels"][: canonical_run.gwf.modelgrid.nlay]
    figure = canonical_run.plot.section(line=line, backend="mpl", **options)
    picture = canonical_run.plot.section(line=line, backend="plotly", **options)
    try:
        _assert_same_drawing(_mpl_drawing(figure), _plotly_drawing(picture))
    finally:
        plt.close(figure)


def test_an_array_fill_is_the_same_on_both(canonical_run):
    """A field of your own -- here a ramp over every cell -- colours the same cells alike."""

    import matplotlib.pyplot as plt

    grid = canonical_run.gwf.modelgrid
    ramp = np.arange(grid.nlay * grid.ncpl, dtype=float).reshape(grid.nlay, grid.ncpl)
    line = _mid_line(canonical_run)
    figure = canonical_run.plot.section(line=line, fill=ramp, backend="mpl", fill_label="ramp")
    picture = canonical_run.plot.section(line=line, fill=ramp, backend="plotly", fill_label="ramp")
    try:
        _assert_same_drawing(_mpl_drawing(figure), _plotly_drawing(picture))
    finally:
        plt.close(figure)


def test_a_bare_grid_draws_the_same_on_both(canonical_run):
    """The grid scope reaches the same picture: geometry by layer, no water levels."""

    import matplotlib.pyplot as plt

    line = _mid_line(canonical_run)
    figure = canonical_run.vor.plot.section(line, fill="layer", backend="mpl")
    picture = canonical_run.vor.plot.section(line, fill="layer", backend="plotly")
    try:
        drawn = _plotly_drawing(picture)
        _assert_same_drawing(_mpl_drawing(figure), drawn)
        assert not drawn["segments"], "a bare grid has no results, so no water level"
    finally:
        plt.close(figure)


def test_the_plotly_section_is_a_picture_with_hover(canonical_run):
    """It answers the picture contract, and every drawn cell answers a hover."""

    from myflopy.viz import Fig, Picture

    picture = canonical_run.plot.section(
        line=_mid_line(canonical_run), fill="layer", backend="plotly"
    )
    assert isinstance(picture, FilledSection) and isinstance(picture, Picture)
    assert isinstance(picture.fig, Fig), "a house figure, not a bare go.Figure"
    hover = [t for t in picture.fig.data if t.meta == "hover"]
    assert len(hover) == 1
    cells = sum(len(_split(t.x, t.y)) for t in picture.fig.data if t.meta == "cells")
    assert len(hover[0].x) == cells, "a drawn cell with no hover"
    names = {row[0] for row in hover[0].customdata}
    nlay = canonical_run.gwf.modelgrid.nlay
    assert names <= {f"Layer {k + 1}" for k in range(nlay)}


@pytest.mark.parametrize(
    ("legend", "expected"),
    [
        ("outside right", {"xanchor": "left"}),
        ("outside bottom", {"yanchor": "top", "orientation": "h"}),
        ("topleft", {"x": 0.0, "y": 1.0, "xanchor": "left", "yanchor": "top"}),
        (True, {"x": 1.0, "y": 1.0, "xanchor": "right", "yanchor": "top"}),
    ],
)
def test_one_legend_word_places_the_legend_on_both(legend, expected, canonical_run):
    """The Plotly legend reads the same resolved placement Matplotlib does."""

    layout = canonical_run.plot.section(
        line=_mid_line(canonical_run), fill="layer", backend="plotly", legend=legend
    ).fig.layout.legend
    for key, value in expected.items():
        assert layout[key] == value, f"legend={legend!r}: {key} is {layout[key]!r}"
    if str(legend).startswith("outside"):
        outside = layout.x > 1 or layout.x < 0 or layout.y > 1 or layout.y < 0
        assert outside, f"legend={legend!r} was placed inside the plot"


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"legend": "sideways"}, "not a placement"),
        ({"fill_cmap": "not-a-colormap"}, "not a Matplotlib colormap"),
        ({"layers": [99]}, "outside this model"),
        ({"fill": "nope"}, "fill must be"),
    ],
)
def test_a_bad_argument_raises_at_the_call_not_at_first_draw(options, message, canonical_run):
    """`.fig` is built lazily, so the Plotly picture validates in its constructor."""

    with pytest.raises(ValueError, match=message):
        canonical_run.plot.section(
            line=_mid_line(canonical_run), backend="plotly", **{"fill": "layer", **options}
        )


def test_a_line_that_misses_the_grid_raises_at_the_call(canonical_run):
    xmin, ymin, xmax, ymax = canonical_run.vor.gdf_vorPolys.total_bounds
    far = shp.LineString([(xmax + 1e5, ymin), (xmax + 2e5, ymax)])
    with pytest.raises(ValueError, match="crosses no cells"):
        canonical_run.plot.section(line=far, fill="layer", backend="plotly")


@pytest.mark.parametrize("backend", ["mpl", "plotly"])
def test_layer_labels_name_the_water_levels_too(backend, canonical_run):
    """`layer_labels=` reached the layer entries and stopped there: with two water
    levels the legend read "sand, clay, ... Head, Layer 1, Head, Layer 2". Both
    renderers take the levels' names from the same list now."""

    import matplotlib.pyplot as plt

    nlay = canonical_run.gwf.modelgrid.nlay
    names = [f"unit {k}" for k in range(nlay)]
    drawn = canonical_run.plot.section(
        line=_mid_line(canonical_run), fill="layer", backend=backend,
        layer_labels=names, head_layers=[0, 1],
    )
    legend = _mpl_drawing(drawn)["legend"] if backend == "mpl" else _plotly_drawing(drawn)["legend"]
    if backend == "mpl":
        plt.close(drawn)
    assert legend == [*names, "Head, unit 0", "Head, unit 1"]
