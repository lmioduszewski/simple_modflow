"""`section()` draws the model's cross-section by default (2026-09-26).

The user's call: "I want model.plot.section(line=line) to default to fill='layer'
and backend='mpl'". So `fill="layer"` and `backend="mpl"` are the defaults at every
scope that names them -- the free verb, `model.plot`, `vor.plot` -- and `fill=None`
is the one spelling of the line profile. A profile argument given without it
raises and names `fill=None`, rather than the picture depending on which other
arguments were passed.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402
import shapely as shp  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from myflopy import plot  # noqa: E402


def _mid_line(model):
    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    middle = (ymin + ymax) / 2
    return shp.LineString([(xmin + 1, middle), (xmax - 1, middle)])


def _png(figure) -> bytes:
    import io

    import matplotlib.pyplot as plt

    buffer = io.BytesIO()
    figure.savefig(buffer, format="png")
    plt.close(figure)
    return buffer.getvalue()


@pytest.mark.parametrize("scope", ["free", "model", "grid"])
def test_the_default_is_the_matplotlib_cross_section_by_layer(scope, canonical_run):
    """A bare call is byte-for-byte `fill="layer", backend="mpl"`."""

    line = _mid_line(canonical_run)
    call = {
        "free": lambda **kw: plot.section(canonical_run, line, **kw),
        "model": lambda **kw: canonical_run.plot.section(line=line, **kw),
        "grid": lambda **kw: canonical_run.vor.plot.section(line, **kw),
    }[scope]
    default = call()
    assert isinstance(default, Figure), f"{scope}: the default is not a Matplotlib figure"
    assert any(c.get_array() is not None for c in default.axes[0].collections), (
        f"{scope}: the default drew no cells"
    )
    assert _png(default) == _png(call(fill="layer", backend="mpl"))


def test_fill_none_is_the_line_profile(canonical_run):
    """No cells: the head lines along the line, and the model top."""

    figure = canonical_run.plot.section(line=_mid_line(canonical_run), fill=None)
    axes = figure.axes[0]
    assert not axes.collections, "fill=None drew cells"
    assert axes.lines, "fill=None drew no profile"


def test_the_matplotlib_profile_honours_its_top_and_bottom_flags(canonical_run):
    """Both were forced by the shared frame builder -- the top always drawn, the
    bottoms never -- so on Matplotlib they were accepted and ignored while the
    Plotly picture honoured them."""

    line = _mid_line(canonical_run)

    def labels(**kw):
        figure = canonical_run.plot.section(line=line, fill=None, **kw)
        legend = figure.axes[0].get_legend()
        return [t.get_text() for t in legend.get_texts()] if legend else []

    assert "model top" in labels()
    assert "model top" not in labels(show_model_top=False)
    assert "Lyr 0 Btm" in labels(show_model_btm=True)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        # a profile argument on the default picture names the profile's spelling
        ({"interpolate": True}, "Pass fill=None"),
        ({"num_points": 50}, "Pass fill=None"),
        ({"layer": [0, 1]}, r"line profile draws \(fill=None\)"),
        # the cross-section's arguments on the profile name the cross-section
        ({"fill": None, "layers": [0]}, "cells of the CROSS-SECTION"),
        ({"fill": None, "fill_cmap": "magma"}, "cells of the CROSS-SECTION"),
        # `.ani` belongs to the Plotly profile alone
        ({"fill": None, "animation_kstpkpers": [(0, 0)]}, "backend='plotly'"),
    ],
)
def test_an_argument_for_the_other_picture_raises_and_names_it(kwargs, message, canonical_run):
    with pytest.raises(ValueError, match=message):
        canonical_run.plot.section(line=_mid_line(canonical_run), **kwargs)


@pytest.mark.parametrize(
    ("kwargs", "name"),
    [({"per": 1}, "per"), ({"kstpkper": (0, 0)}, "kstpkper"), ({"head_layers": [0, 1]}, "head_layers")],
)
def test_a_bare_grids_cross_section_refuses_what_needs_results(kwargs, name, canonical_run):
    """With the cross-section the default, `plot.section(vor, line, per=3)` would
    draw the same geometry as without `per` -- taken and dropped."""

    with pytest.raises(ValueError, match=f"{name} need a model's results"):
        plot.section(canonical_run.vor, _mid_line(canonical_run), **kwargs)


def test_the_grid_scope_outlines_refuse_the_cross_sections_colours(canonical_run):
    with pytest.raises(ValueError, match="fill=None draws only their outlines"):
        canonical_run.vor.plot.section(_mid_line(canonical_run), fill=None, fill_cmap="magma")


def test_the_profile_still_animates_on_plotly(canonical_run):
    """`.ani` survives the move: it is the Plotly line profile's, as before."""

    picture = canonical_run.plot.section(
        line=_mid_line(canonical_run), fill=None, backend="plotly"
    )
    assert len(picture.ani.frames) == len(canonical_run.kstpkper)
    assert np.isfinite(len(picture.fig.data))
