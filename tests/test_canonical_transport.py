"""The canonical valley model with a solute-transport sibling (plan §6.1/6.2 item 6).

Ledger 56 deferred this as "large -- the canonical runs on the single-model
``SimulationBase``, so a coupled GWT needs the multi-model spec path". Measured
2026-07-28: that was wrong twice over. ``SimulationBase`` already owns a real
``MFSimulation`` that accepts a GWT sibling, and the spec path would have made
calibration WORSE, not possible -- it puts each model's arrays in a subdirectory
where PEST's root-globbing parameter resolution cannot find them.
"""

from __future__ import annotations

import pytest

from myflopy.modflow.mf6.canonical import CANONICAL_MODEL_CONTRACT
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.canonical_transport import (
    MAX_MF6_MODEL_NAME,
    build_canonical_transport_model,
    transport_model_name,
)


@pytest.fixture(scope="module")
def transport_model(tmp_path_factory):
    """Build + run the canonical model with transport attached (once per module)."""

    return build_canonical_transport_model(
        tmp_path_factory.mktemp("canon_transport") / "model",
        config=CanonicalModelConfig.testing(),
    )


def test_the_sibling_name_fits_mf6s_limit():
    """MF6 rejects a MODELNAME over 16 characters, and the canonical name is 14.

    A plain ``f"{flow}_trans"`` produced ``viz_prt_master_trans`` (20) and MF6
    refused the whole simulation before solving anything.
    """

    assert transport_model_name("viz_prt_master") == "viz_prt_master_t"
    assert len(transport_model_name("viz_prt_master")) <= MAX_MF6_MODEL_NAME
    long_name = "a" * 40
    assert len(transport_model_name(long_name)) == MAX_MF6_MODEL_NAME


@pytest.mark.slow
def test_the_canonical_contract_still_holds_with_transport_attached(transport_model):
    """Attaching a sibling must not disturb the model every other test shares.

    ``CANONICAL_MODEL_CONTRACT.validate()`` is entirely ``model.gwf``-scoped, so
    a GWT sibling in the same simulation is invisible to it. Asserting that is
    what makes it safe to attach transport to the canonical family at all.
    """

    CANONICAL_MODEL_CONTRACT.validate(transport_model.model)
    assert transport_model.model.gwf.model_type == "gwf6"


@pytest.mark.slow
def test_the_simulation_stays_flat_so_pest_can_find_its_files(transport_model):
    """The layout decision this whole design turns on.

    ``Project.prepare_run`` gives each model its own subdirectory, so external
    arrays land at ``<ws>/flow/flow.npf_k.txt`` -- but PEST's parameter-file
    resolution globs the workspace ROOT, so calibration would break. A sibling on
    the existing ``SimulationBase`` simulation stays flat. Both binaries also sit
    at the root, which is what lets a forward-run post-processor read the ``.ucn``
    the same way it reads the ``.hds``.
    """

    workspace = transport_model.workspace
    flow_name = transport_model.model.name
    transport_name = transport_model.transport_name

    assert (workspace / f"{flow_name}.hds").exists()
    assert (workspace / f"{transport_name}.ucn").exists()

    # ...and BOTH models' external arrays land beside them, not one level down.
    # (K splits into per-layer files, which is what `parameterize('k')` globs.)
    transport_model.model.sim.set_all_data_external()
    transport_model.model.sim.write_simulation(silent=True)

    assert not [p for p in workspace.iterdir() if p.is_dir()], (
        "a per-model subdirectory appeared -- PEST's resolver globs the root only"
    )
    assert sorted(p.name for p in workspace.glob(f"{flow_name}.npf_k_layer*.txt")), (
        "the flow model's K arrays are not at the workspace root"
    )
    assert (workspace / f"{transport_name}.mst_porosity.txt").exists(), (
        "the transport model's porosity array is not at the workspace root"
    )


@pytest.mark.slow
def test_the_transport_model_carries_a_real_plume(transport_model):
    """A source that produces no plume would make every downstream test vacuous."""

    view = transport_model.transport_view()
    assert view.model_type == "gwt6"

    frame = view.conc.get()
    assert not frame.empty
    assert frame["conc"].max() == pytest.approx(
        transport_model.source_concentration, rel=1e-6
    )
    # the plume has actually spread beyond the handful of source cells
    spread = int((frame["conc"] > 0.01).sum())
    assert spread > 10 * len(transport_model.source_cells), (
        f"only {spread} cell-periods above 1% -- the plume did not develop"
    )


@pytest.mark.slow
def test_every_advanced_flow_package_got_its_transport_counterpart(transport_model):
    """MF6 refuses to run without them, and each one shows up in the budget.

    The canonical model carries LAK, SFR, UZF and MVR; a GWT sibling without
    LKT/SFT/UZT/MVT fails with "GWF water mover is active but the GWT MVT package
    has not been specified". Asserting on the BUDGET rather than on the package
    list proves they are actually participating, not merely present.
    """

    view = transport_model.transport_view()
    terms = {str(term).strip().upper() for term in view.budget.types}

    assert {"LKT", "SFT", "UZT"} <= terms, f"advanced transport terms missing: {terms}"
    assert "CNC" in terms, "the contaminant source is not in the budget"
    assert "SOURCE-SINK MIX" in terms
