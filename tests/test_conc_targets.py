"""``ConcTargets`` -- concentration observations for transport calibration.

Plan §6.1/6.2 item 5 (the observation half). ``ConcTargets`` is a configuration
of ``HeadTargets`` rather than a copy, because concentration is sampled exactly
the way head is: a point value at a ``(layer, cell)`` per stress period.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.canonical_transport import (
    build_canonical_transport_calibration_demo,
    monitoring_well_cells,
)
from myflopy.modflow.mf6.observations import ConcTargets, HeadTargets


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    """Truth-sampled transport calibration demo (built + run once)."""

    return build_canonical_transport_calibration_demo(
        tmp_path_factory.mktemp("conc_cal") / "model",
        config=CanonicalModelConfig.testing(),
    )


def test_conc_targets_reads_concentration_not_heads():
    """The one thing that must differ from HeadTargets, asserted directly."""

    assert issubclass(ConcTargets, HeadTargets)
    assert (ConcTargets._TABLE_ATTRIBUTE, ConcTargets._STORE_COLUMN) == (
        "all_conc",
        "conc",
    )
    assert (HeadTargets._TABLE_ATTRIBUTE, HeadTargets._STORE_COLUMN) == (
        "all_heads",
        "elev",
    )
    # inherited wholesale -- a copy would show these in ConcTargets' own body
    for shared in ("match_to_model", "compare", "residuals", "stats", "to_long"):
        assert shared not in vars(ConcTargets), f"ConcTargets re-implements {shared}"


@pytest.mark.slow
def test_monitoring_wells_are_placed_where_concentration_can_move(demo):
    """Wells on the SOURCE cells would be pinned at 1.0 and constrain nothing.

    CNC holds the source cells at the source concentration for the whole run, so
    a well there reports the same number no matter what K is. The selector keeps
    intermediate-concentration cells downgradient of the source instead.
    """

    source = set(demo.transport.source_cells)
    wells = set(demo.well_cells)

    assert wells, "no monitoring wells were selected"
    assert not (wells & source), "a monitoring well sits on a source cell"

    long = demo.conc_targets.to_long()
    values = pd.to_numeric(long["head_target"], errors="coerce").dropna()
    assert values.max() < demo.transport.source_concentration, (
        "a target sits at the source concentration -- that well cannot move"
    )
    assert values.max() > 0.01, "every target is ~zero: the plume never reached the wells"


@pytest.mark.slow
def test_the_spoiled_model_actually_misfits_the_truth_targets(demo):
    """The demo is only a calibration if the starting model is WRONG.

    Truth-derived targets sampled from the model that produced them would fit
    perfectly; the point of the K perturbation is to open a gap for history
    matching to close. A demo whose residuals are all zero would run PEST to
    completion and prove nothing.
    """

    comparison = demo.conc_targets.compare(demo.transport.transport_view())
    residual = comparison["residual"].abs()

    assert len(comparison) > 0
    assert int((residual > 1e-6).sum()) == len(comparison), (
        "some residuals are zero -- the perturbed model still matches truth there"
    )
    assert residual.max() > 0.01, f"largest misfit is only {residual.max()}"
    assert demo.start_k_factor != 1.0


@pytest.mark.slow
def test_truth_sampled_targets_round_trip_to_zero_residual(demo):
    """Sampling a model and comparing against it must give exactly zero.

    This is the control for the test above: it proves the misfit there comes
    from the K perturbation, not from a units or indexing error in the sampler.
    """

    view = demo.transport.transport_view()
    cells = monitoring_well_cells(demo.transport, count=4)
    names = [f"chk_{i:02d}" for i in range(len(cells))]
    locations = [
        {"name": n, "layer": 0, "cell": int(c)}
        for n, c in zip(names, cells, strict=False)
    ]
    nper = int(demo.model.nper)
    placeholder = pd.DataFrame(
        {"per": list(range(nper)), **{n: [np.nan] * nper for n in names}}
    )

    sampler = ConcTargets(locations=locations, values=placeholder, time_column="per")
    simulated = sampler.simulated_heads(view)
    assert list(simulated.columns) == ["per", *names]

    round_trip = ConcTargets(locations=locations, values=simulated, time_column="per")
    assert round_trip.compare(view)["residual"].abs().max() == pytest.approx(0.0, abs=1e-9)


@pytest.mark.slow
def test_the_forecast_is_a_separate_downgradient_prediction(demo):
    """The forecast cell must not also be an observation, or it is not a forecast."""

    assert demo.forecast_cell not in demo.well_cells
    assert not demo.forecast_targets.to_long().empty
