"""A correlated prior for pilot points (plan §5.8 item 3, swapped).

The backlog item was "Tikhonov/preferred-value regularization helpers". Measured
against the runner this repo actually ships, that item is not merely low-value —
it is harmful: PESTPP-IES ignores prior-information equations, and a version-2
control file in ``regularization`` mode fails to PARSE, so a ``cal.regularize()``
would have broken ``run_ies`` outright.

The real gap underneath it: ``correlation=`` was documented for
``style="pilotpoints"`` (``pilot_points.py:14``) and read nowhere in that module,
so the variogram a caller asked for had no effect on the prior at all.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig


def _prior_correlation_by_distance(cal, reals=60):
    """Mean prior correlation of near vs far pilot pairs, in log space."""

    cal._inject_geostatistical_prior(reals)
    path = cal.template_workspace / f"{cal.name}.prior_par.csv"
    if not path.exists():
        return None
    ensemble = pd.read_csv(path, index_col=0)
    frame = cal._pilot_point_frames["k"][0][0].set_index("parnme")
    names = [name for name in frame.index if name in ensemble.columns]
    if len(names) < 4:
        return None

    x = frame.loc[names, "x"].to_numpy(dtype=float)
    y = frame.loc[names, "y"].to_numpy(dtype=float)
    correlations = np.corrcoef(np.log10(ensemble[names].to_numpy()).T)
    distances = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])

    upper = np.triu_indices(len(names), 1)
    pair_distance, pair_correlation = distances[upper], correlations[upper]
    near = pair_distance < np.percentile(pair_distance, 25)
    far = pair_distance > np.percentile(pair_distance, 75)
    return float(pair_correlation[near].mean()), float(pair_correlation[far].mean())


@pytest.mark.slow
@pytest.mark.parametrize(
    "correlation,expect_correlated",
    [(150.0, False), (3000.0, True)],
    ids=["range-below-point-spacing", "range-spanning-the-domain"],
)
def test_the_variogram_range_shapes_the_pilot_point_prior(
    tmp_path, correlation, expect_correlated
):
    """The behavioural proof that ``correlation=`` is read.

    Nearest pilot pair on this fixture is ~372 m. A 150 m range is shorter than
    that, so neighbouring points should be effectively independent; a 3000 m
    range spans the domain, so they should move together and decay with
    distance. Before this, both produced the same uncorrelated draw.
    """

    pytest.importorskip("pyemu")
    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    cal = demo.model.pest("pp", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_space=5, correlation=correlation,
                     bounds=(0.1, 10.0), physical=(1e-3, 300.0))
    cal.observe(demo.head_targets)
    cal.build("pp.pst", noptmax=0)

    measured = _prior_correlation_by_distance(cal)
    assert measured is not None, "no prior ensemble was written for pilot points"
    near, far = measured

    if expect_correlated:
        assert near > 0.5, f"neighbouring pilot points are not correlated ({near:.3f})"
        assert near > far, "correlation should decay with distance"
    else:
        assert abs(near) < 0.25, (
            f"a variogram range below the point spacing should leave points "
            f"nearly independent, got {near:.3f}"
        )


@pytest.mark.slow
def test_no_correlation_means_no_injected_prior(tmp_path):
    """`correlation=` is what opts a run into a drawn prior. Without it,
    PESTPP-IES draws from bounds itself and there is nothing to inject —
    writing an ensemble anyway would silently change a run nobody asked to
    change."""

    pytest.importorskip("pyemu")
    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    cal = demo.model.pest("nocorr", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_space=5, bounds=(0.1, 10.0))
    cal.observe(demo.head_targets)
    cal.build("nocorr.pst", noptmax=0)

    cal._inject_geostatistical_prior(20)
    assert not (cal.template_workspace / "nocorr.prior_par.csv").exists()
    assert "ies_par_en" not in cal.pst.pestpp_options


def test_a_regularization_mode_control_file_is_refused_before_launch():
    """PESTPP-IES ignores prior-information equations AND, in the version-2
    format myflopy writes, fails to parse a `* regularization` block at all —
    measured: `pestpp-ies reg.pst` exits 1 with `control file parsing error`
    before any model run. Catching it here names the cause; letting it through
    names a keyword the caller never typed."""

    from types import SimpleNamespace

    from myflopy.modflow.mf6.pest.project import PestProject

    project = object.__new__(PestProject)
    project.pst = SimpleNamespace(
        control_data=SimpleNamespace(pestmode="regularization")
    )
    with pytest.raises(ValueError, match="ies_reg_factor"):
        project._refuse_regularization_mode()

    # The ordinary mode passes straight through.
    project.pst.control_data.pestmode = "estimation"
    assert project._refuse_regularization_mode() is None
