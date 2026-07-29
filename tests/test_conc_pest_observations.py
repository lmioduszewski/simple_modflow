"""Concentration observations through the PEST build (plan §6.1/6.2 item 5).

The end-to-end proof: build a real ``.pst`` whose history-matching targets are
CONCENTRATIONS, run the injected forward run once, and check it regenerates the
simulated file pyEMU's instruction file expects.
"""

from __future__ import annotations

import json
import subprocess
import sys

import geopandas as gpd
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.canonical_transport import (
    build_canonical_transport_calibration_demo,
)
from myflopy.modflow.mf6.observations import ConcTargets, HeadTargets
from myflopy.modflow.mf6.pest.project import (
    _OBSERVATION_SPEC_TYPES,
    _TARGET_SPEC_TYPES,
    METADATA_FILENAME,
)
from myflopy.modflow.mf6.pest.specs import ConcObservationSpec


def test_conc_targets_are_matched_before_head_targets():
    """Ladder ORDER is load-bearing, because ConcTargets subclasses HeadTargets.

    ``_coerce_observation_spec`` walks ``_TARGET_SPEC_TYPES`` with ``isinstance``.
    With ``HeadTargets`` first, every ``ConcTargets`` would match it and be
    wrapped as a HEAD spec -- silently reading ``.hds`` and calibrating against
    the wrong physics.
    """

    order = [target.__name__ for target, _ in _TARGET_SPEC_TYPES]
    assert order.index("ConcTargets") < order.index("HeadTargets")
    assert issubclass(ConcTargets, HeadTargets)
    assert ConcObservationSpec in _OBSERVATION_SPEC_TYPES


def test_observe_wraps_bare_conc_targets(tmp_path):
    """``cal.observe(conc_targets)`` must produce a ConcObservationSpec."""

    locations = [{"name": "cw_00", "layer": 0, "cell": 3}]
    values = pd.DataFrame({"per": [0, 1], "cw_00": [0.1, 0.2]})
    targets = ConcTargets(locations=locations, values=values, time_column="per")

    from myflopy.modflow.mf6.pest.project import PestProject

    spec = PestProject._coerce_observation_spec(
        object.__new__(PestProject), targets, prefix=None
    )
    assert isinstance(spec, ConcObservationSpec)
    assert spec.prefix == "conc"


@pytest.mark.slow
def test_a_concentration_calibration_builds_and_its_forward_run_works(tmp_path):
    """The whole point: a control file history-matched against concentrations.

    Asserts three things a broken wiring would each fail differently:

    * ``obsval`` is never NaN. The PEST value assignment reads
      ``getattr(row, "value", getattr(row, "head_target", np.nan))``, so a target
      whose value column was renamed without updating that fallback yields a
      control file full of NaN observations -- which BUILDS and RUNS and
      calibrates to nothing (ledger 102).
    * the injected forward run exits 0 and regenerates the simulated CSV, reading
      the TRANSPORT model's ``.ucn`` rather than the flow model's ``.hds``.
    * the regenerated values are not all zero, which is what reading the wrong
      binary (or the wrong layer) would produce.
    * the location snapshot and its metadata are written, which is what puts
      concentration on ``plot_obs_residuals`` (ledger 105). ``match_to_model``
      already computed the coordinates; before this they were discarded and a
      concentration-only calibration reviewed as a blank map.
    """

    demo = build_canonical_transport_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )
    cal = demo.model.pest("transport_pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.05, 2.0), physical=(0.01, 300.0))
    cal.observe(demo.conc_targets)
    cal.forecast(demo.forecast_targets)
    pst = cal.build("transport_pest.pst", noptmax=0)

    assert pst.npar_adj == 1
    assert pst.nnz_obs == len(demo.conc_targets.to_long())

    observations = pst.observation_data
    assert int(observations["obsval"].isna().sum()) == 0, (
        "NaN obsval: the target value column never reached the PEST builder"
    )
    assert any("oname:conc" in name for name in observations.index)

    template = cal.template_workspace
    forward_run = (template / "forward_run.py").read_text()
    assert "_write_conc_target_csv(" in forward_run
    assert demo.transport.transport_name in forward_run, (
        "the forward run does not name the GWT model, so it cannot find the .ucn"
    )

    result = subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    regenerated = template / "conc_simulated_conc.csv"
    assert regenerated.exists()
    frame = pd.read_csv(regenerated)
    assert list(frame.columns)[0] == "per"
    assert (frame.drop(columns=["per"]).abs().to_numpy() > 1e-9).any(), (
        "every regenerated concentration is zero -- the wrong binary was read"
    )

    # --- the residual map's inputs (ledger 105) ---
    metadata = json.loads((template / METADATA_FILENAME).read_text())
    conc_set = next(
        entry for entry in metadata["observation_sets"]
        if entry["kind"] == "conc_targets"
    )
    assert conc_set["geometry"] == "points"

    locations = gpd.read_file(template / conc_set["locations_file"])
    assert not locations.empty
    assert locations.geometry.x.notna().all(), (
        "no coordinates snapshotted -- every concentration target would be "
        "dropped from the residual map"
    )
    mapping = pd.read_csv(template / conc_set["mapping_file"])
    assert set(mapping["name"]) == set(locations["name"])


@pytest.mark.slow
def test_a_flow_only_simulation_says_why_concentration_cannot_be_observed(tmp_path):
    """A clear error beats reading a .ucn that does not exist."""

    from myflopy.modflow.mf6.canonical_example import build_canonical_model
    from myflopy.modflow.mf6.pest.observations import resolve_transport_model_name

    model = build_canonical_model(
        tmp_path / "flow_only", config=CanonicalModelConfig.testing()
    )

    class _Stub:
        pass

    stub = _Stub()
    stub.model = model
    with pytest.raises(ValueError, match="No GWT model found"):
        resolve_transport_model_name(stub)


def test_conc_targets_are_on_the_public_surface():
    """Both names ship from the top level, like their head counterparts."""

    assert mf.ConcTargets is ConcTargets
    assert mf.ConcObservationSpec is ConcObservationSpec
    assert {"ConcTargets", "ConcObservationSpec"} <= set(mf.__preferred__)
