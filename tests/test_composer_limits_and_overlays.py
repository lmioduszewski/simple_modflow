"""The composers keep what each panel was asked for (ledger 306).

`mosaic`/`animate` NAME `zmin`/`zmax` and forward them to every panel (ledger
258), and the panel honours them -- then the composer pooled the colour axis over
the data regardless, so `hds.mosaic(zmin=100, zmax=120)` drew 79.04..153.76 on
BOTH backends while `animate` honoured the same arguments. And the static
mosaic/animation re-plotted the cells themselves instead of calling the panel's
own `plot_mpl`, which dropped every `features=` layer and the `select=` outline.
"""

from __future__ import annotations

import numpy as np
import pytest
from shapely.geometry import Point

from myflopy import viz

WELL = Point(1050.0, 1050.0)


def _cell_norm(axis):
    """The colour norm of the cell collection a static panel drew."""

    return axis.collections[0].norm


# --- the colour limits reach the shared axis ------------------------------------


def test_a_plotly_mosaic_honours_zmin_and_zmax(canonical_run):
    fig = canonical_run.hds.mosaic(by="layer", zmin=100.0, zmax=120.0)
    assert (fig.layout.coloraxis.cmin, fig.layout.coloraxis.cmax) == (100.0, 120.0)


def test_a_static_mosaic_honours_zmin_and_zmax(canonical_run):
    fig = canonical_run.hds.mosaic(by="layer", zmin=100.0, zmax=120.0, backend="mpl")
    panels = [axis for axis in fig.axes if axis.get_label() != "<colorbar>" and axis.collections]
    assert len(panels) == 4
    for axis in panels:
        norm = _cell_norm(axis)
        assert (norm.vmin, norm.vmax) == (100.0, 120.0)


def test_one_limit_leaves_the_other_pooled(canonical_run):
    """`zmin` alone must not freeze the top as well -- or drop the bottom."""

    pooled = canonical_run.hds.mosaic(by="layer")
    fig = canonical_run.hds.mosaic(by="layer", zmin=100.0)
    assert fig.layout.coloraxis.cmin == 100.0
    assert fig.layout.coloraxis.cmax == pytest.approx(pooled.layout.coloraxis.cmax)


def test_a_log_mosaic_takes_its_limits_in_data_units(canonical_run):
    fig = canonical_run.packages.npf.k.mosaic(logscale=True, zmin=1.0, zmax=100.0)
    # The trace axis is log10; the caller wrote data units.
    assert (fig.layout.coloraxis.cmin, fig.layout.coloraxis.cmax) == (0.0, 2.0)
    static = canonical_run.packages.npf.k.mosaic(
        logscale=True, zmin=1.0, zmax=100.0, backend="mpl"
    )
    norm = _cell_norm(static.axes[0])
    assert (norm.vmin, norm.vmax) == (1.0, 100.0)


def test_the_free_mosaic_takes_the_widest_panel_limits(canonical_run):
    """Panels a caller built separately may disagree; clipping one would misreport it."""

    fig = viz.mosaic([
        canonical_run.plot.map(layer=0, zmin=100.0, zmax=120.0),
        canonical_run.plot.map(layer=1, zmin=90.0, zmax=110.0),
    ])
    assert (fig.layout.coloraxis.cmin, fig.layout.coloraxis.cmax) == (90.0, 120.0)


# --- the static composers draw what the static map draws ------------------------


def test_every_static_mosaic_panel_draws_its_features_and_selection(canonical_run):
    single = canonical_run.packages.npf.k.map(backend="mpl", features=WELL, select=[0, 1, 2])
    expected = len(single.axes[0].collections)
    assert expected >= 3  # cells, the selection outline, the well

    fig = canonical_run.packages.npf.k.mosaic(backend="mpl", features=WELL, select=[0, 1, 2])
    panels = [axis for axis in fig.axes if axis.get_label() != "<colorbar>" and axis.collections]
    assert len(panels) == 4
    assert [len(axis.collections) for axis in panels] == [expected] * 4


def test_a_static_animation_frame_draws_its_features_under_the_shared_scale(canonical_run):
    animation = canonical_run.hds.animate(backend="mpl", features=WELL, zmin=100.0, zmax=120.0)
    axis = animation._fig.axes[0]
    assert len(axis.collections) >= 2  # cells and the well
    norm = _cell_norm(axis)
    assert (norm.vmin, norm.vmax) == (100.0, 120.0)


def test_plot_mpl_takes_a_shared_norm(canonical_run):
    from matplotlib.colors import Normalize

    norm = Normalize(vmin=0.0, vmax=1.0)
    fig = canonical_run.plot.map(layer=0).plot_mpl(norm=norm, colorbar=False)
    assert _cell_norm(fig.axes[0]) is norm
    assert np.isclose(_cell_norm(fig.axes[0]).vmax, 1.0)


def test_symmetric_panel_limits_stay_centred_in_a_mosaic(canonical_run):
    """Ledger 71's centering half: a `change` map is log-scaled with SYMMETRIC
    limits about 1x, and pooling over the data put the white point at 0.231 of
    the diverging scale, so increases rendered red. Honouring each panel's limits
    (widest of symmetric pairs is symmetric) keeps it at 0.5 without `diff=True`."""

    ncpl = canonical_run.vor.ncpl
    rng = np.random.default_rng(0)
    panels = [
        canonical_run.plot.map(
            values=10.0 ** rng.uniform(-0.3, 1.0, ncpl), logscale=True,
            zmin=10.0 ** -half, zmax=10.0 ** half, colorscale="RdBu",
        )
        for half in (1.0, 0.5)
    ]
    axis = viz.mosaic(panels).layout.coloraxis
    assert (axis.cmin, axis.cmax) == (-1.0, 1.0)
    assert (0.0 - axis.cmin) / (axis.cmax - axis.cmin) == 0.5
