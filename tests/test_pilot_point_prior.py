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
    frame = cal._pilot_point_sets[0]["frame"].set_index("parnme")
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
    cal.parameterize("k", style="pilotpoints", pp_spacing=500.0, correlation=correlation,
                     bounds=(0.1, 10.0), physical=(1e-3, 300.0))
    cal.observe(demo.head_targets)
    cal.build("pp.pst")

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
    cal.parameterize("k", style="pilotpoints", pp_spacing=500.0, bounds=(0.1, 10.0))
    cal.observe(demo.head_targets)
    cal.build("nocorr.pst")

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


def test_the_pilot_point_prior_spreads_as_far_as_its_bounds_say(tmp_path):
    """The pilot draw took the variogram's sill of 1 as its log10 variance, whatever the
    bounds: with bounds (0.5, 2) 77% of the values handed to PESTPP-IES lay outside
    them and it ran each AT a bound. It is pyEMU's geostatistical draw now, the one a
    correlated ``grid`` parameter of the same bounds gets: log10 standard deviation
    ``log10(ub / lb) / 6`` (ledger 424)."""

    pytest.importorskip("pyemu")
    import myflopy as mf
    from tests.test_pest_structured_grid import _build, _dis_run, _heads

    model = mf.load_mf6_run(_dis_run(tmp_path / "run"))
    cal = model.pest("bounds", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_spacing=10.0, correlation=30.0,
                     bounds=(0.5, 2.0), name="kp")
    cal.parameterize("k33", style="grid", correlation=30.0, bounds=(0.5, 2.0), name="kg")
    cal.observe(_heads())
    _build(cal, "bounds")
    np.random.seed(11)
    cal._inject_geostatistical_prior(200)
    ensemble = pd.read_csv(cal.template_workspace / "bounds.prior_par.csv", index_col=0)
    expected = np.log10(4.0) / 6
    for prefix in ("kpl", "pname:kg"):
        columns = [name for name in ensemble.columns if name.startswith(prefix)]
        values = ensemble[columns].to_numpy()
        assert len(columns) == 60
        assert np.log10(values).std() == pytest.approx(expected, rel=0.15), prefix
        assert ((values < 0.5) | (values > 2.0)).mean() < 0.01, prefix


def _adjacent_correlation(ensemble, names, x, y, groups, spacing):
    """Mean log10 correlation of parameters one ``spacing`` apart in one group."""

    values = np.log10(ensemble[names].to_numpy())
    correlations = np.corrcoef(values.T)
    distances = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    adjacent = (groups[:, None] == groups[None, :]) & np.isclose(distances, spacing)
    return float(correlations[adjacent].mean()), float(values.std(axis=0).mean())


@pytest.mark.parametrize("nugget", [0.0, 0.5, 1.0])
def test_the_nugget_reaches_the_prior_pestpp_ies_receives(tmp_path, nugget):
    """pyEMU's draw replaced the covariance's diagonal with the bound-derived variance,
    discarding a nugget added on top of the sill-1 variogram: neighbours 10 ft apart
    under a 30 ft range correlated 0.718 at nuggets 0, 0.5 and 1 alike, for pilot points
    and a ``grid`` parameter. The structure is scaled to a sill of 1 now, so the
    correlation is ``exp(-h / range) / (1 + nugget)`` and the variance stays the one the
    bounds give (ledger 489)."""

    pytest.importorskip("pyemu")
    import myflopy as mf
    from tests.test_pest_structured_grid import _build, _dis_run, _heads

    model = mf.load_mf6_run(_dis_run(tmp_path / "run"))
    cal = model.pest("nug", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_spacing=10.0, correlation=30.0,
                     bounds=(0.5, 2.0), name="kp", nugget=nugget)
    cal.parameterize("k33", style="grid", correlation=30.0, bounds=(0.5, 2.0), name="kg",
                     nugget=nugget)
    cal.observe(_heads())
    _build(cal, "nug")
    np.random.seed(7)
    cal._inject_geostatistical_prior(1000)
    ensemble = pd.read_csv(cal.template_workspace / "nug.prior_par.csv", index_col=0)

    points = pd.concat(item["frame"] for item in cal._pilot_point_sets
                       if item["spec"].name == "kp")
    grid = cal.pst.parameter_data.loc[
        [name for name in ensemble.columns if name.startswith("pname:kg")]]
    expected = np.exp(-10.0 / 30.0) / (1.0 + nugget)
    for what, names, x, y, groups in (
        ("pilot points", list(points["parnme"]), points["x"], points["y"],
         points["parnme"].str[:4]),
        ("grid", list(grid.index), grid["x"], grid["y"], grid["pargp"]),
    ):
        correlation, spread = _adjacent_correlation(
            ensemble, names, np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(groups), 10.0)
        assert correlation == pytest.approx(expected, abs=0.04), what
        assert spread == pytest.approx(np.log10(4.0) / 6, rel=0.05), what


@pytest.mark.parametrize("nugget", [-0.1, float("nan"), float("inf")])
def test_a_nugget_that_is_not_a_variance_is_refused(tmp_path, nugget):
    """A negative nugget would make two distinct points more correlated than one point
    with itself once the structure is scaled to a sill of 1 (ledger 489)."""

    pytest.importorskip("pyemu")
    from myflopy.modflow.mf6.pest.project import PestProject

    project = object.__new__(PestProject)
    project._native_parameter_specs = []
    with pytest.raises(ValueError, match=r"nugget=.*finite number of at least 0\.0"):
        project.parameterize("k", style="grid", correlation=30.0, nugget=nugget)
    assert project._native_parameter_specs == []


@pytest.mark.slow
def test_zoned_pilot_points_are_drawn_correlated_within_a_zone_and_not_across(tmp_path):
    """One covariance block per pilot-point SET -- a (layer, zone) of one call -- so a
    long variogram ties neighbours of one zone together and leaves neighbours across a
    zone edge independent (ledger 527)."""

    pytest.importorskip("pyemu")
    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    xc = np.asarray(demo.model.gwf.modelgrid.xcellcenters, dtype=float).reshape(-1)
    zones = np.where(xc < np.median(xc), 1, 2)
    cal = demo.model.pest("ppzone", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_spacing=300.0, correlation=3000.0,
                     zones=zones, layers=[0], bounds=(0.1, 10.0))
    cal.observe(demo.head_targets)
    cal.build("ppzone.pst")
    np.random.seed(5)
    cal._inject_geostatistical_prior(500)
    ensemble = pd.read_csv(cal.template_workspace / "ppzone.prior_par.csv", index_col=0)

    points = pd.concat(item["frame"] for item in cal._pilot_point_sets)
    names = list(points["parnme"])
    x, y = points["x"].to_numpy(dtype=float), points["y"].to_numpy(dtype=float)
    zone = points["zone"].to_numpy()
    correlations = np.corrcoef(np.log10(ensemble[names].to_numpy()).T)
    distances = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    near = (distances > 0) & (distances < 450.0)
    within = near & (zone[:, None] == zone[None, :])
    across = near & (zone[:, None] != zone[None, :])
    assert within.any() and across.any()
    assert correlations[within].mean() > 0.5
    assert abs(correlations[across].mean()) < 0.15
