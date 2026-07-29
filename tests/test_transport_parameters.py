"""Transport properties as calibration targets (ledger 103).

``cal.parameterize("porosity")`` adjusts a property that lives on the GWT
sibling of the model the calibration hangs off. The interesting part is not the
plumbing -- pyEMU's ``apply_list_and_array_pars`` is filename-driven and does
not care which model wrote a file -- but the three ways it can go SILENTLY
wrong, each of which is pinned below.
"""

from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from pyemu.pst.pst_utils import write_to_template

from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.canonical_transport import (
    build_canonical_transport_calibration_demo,
)
from myflopy.modflow.mf6.pest.native_parameters import (
    NativeParameterSpec,
    _select_layer_files,
    resolve_target,
)
from myflopy.modflow.mf6.pest.project import METADATA_FILENAME


def test_porosity_is_declared_as_a_transport_target():
    """The recipe carries WHICH model owns the file; nothing else has to know."""

    recipe = resolve_target("porosity")
    assert recipe.model == "transport"
    assert (recipe.package, recipe.variable) == ("mst", "porosity")
    assert resolve_target("mst.porosity") is recipe
    assert resolve_target("n") is recipe
    # every flow target keeps the default
    assert resolve_target("k").model == "flow"


def test_pilot_points_on_a_transport_target_are_refused():
    """``add_pilot_point_parameter`` reads ``project.model.gwf.npf.k`` as the
    array to interpolate against, whatever the target is. Pointed at porosity it
    would build, run and calibrate -- against the K field. Refuse at declaration
    time, when the user is looking at the call, not at build time."""

    with pytest.raises(NotImplementedError, match="NPF K"):
        NativeParameterSpec(target="porosity", style="pilotpoints")

    # the same style on the target it WAS written for stays available
    assert NativeParameterSpec(target="k", style="pilotpoints").style == "pilotpoints"


def test_layers_on_a_whole_grid_array_is_refused_not_ignored():
    """The old code filtered per-layer files and fell back to ALL files when
    nothing matched (``... or files`` / ``if selected:``). On a target MF6 writes
    as one whole-grid array there is nothing to match, so ``layers=[0]``
    parameterized every layer while reporting success -- the worst kind of
    wrong, because the control file looks exactly as intended."""

    recipe = resolve_target("porosity")
    with pytest.raises(ValueError, match="single whole-grid array"):
        _select_layer_files(["m_t.mst_porosity.txt"], [0], recipe)


def test_layers_that_match_no_file_names_the_ones_that_exist():
    """Same fallthrough, per-layer case: ``layers=[9]`` on a 4-layer model used
    to silently parameterize all four."""

    files = [f"m.npf_k_layer{i}.txt" for i in (1, 2, 3, 4)]
    recipe = resolve_target("k")

    assert _select_layer_files(files, [0, 2], recipe) == [
        "m.npf_k_layer1.txt", "m.npf_k_layer3.txt"
    ]
    with pytest.raises(ValueError, match=r"Available layers: \[0, 1, 2, 3\]"):
        _select_layer_files(files, [9], recipe)


def test_a_whole_grid_capture_file_does_not_claim_to_be_layer_zero():
    """``layer_prefixes`` names which file holds which layer, and
    ``IesResults._cell_map`` switches on it: populated means "the array index IS
    the cell id", empty means "flat = layer * ncpl + cell". A single whole-grid
    file recorded as ``{0: prefix}`` takes the first branch and hands
    ``plot_field`` nlay*ncpl values for ncpl cells, every one of them labelled
    layer 0.

    Driven directly because the build path relayers first, so this state is not
    reachable from ``build()`` -- it is the guard for an array MF6 cannot store
    layered.
    """

    from types import SimpleNamespace

    from myflopy.modflow.mf6.pest.project import PestProject

    def _capture(filename):
        project = SimpleNamespace(
            pf=SimpleNamespace(add_observations=lambda *a, **k: None),
            _capture_field_specs={},
        )
        spec = NativeParameterSpec(target="porosity", capture=True)
        spec.resolved_files = [filename]
        PestProject._add_capture_field_observations(project, spec)
        return project._capture_field_specs["porosity"]["layer_prefixes"]

    assert _capture("m_t.mst_porosity.txt") == {}
    assert _capture("m_t.mst_porosity_layer3.txt") == {2: "porosityfieldl2"}


@pytest.mark.slow
def test_a_porosity_calibration_builds_and_the_multiplier_reaches_the_model(tmp_path):
    """The end-to-end proof, and the one that would catch a no-op parameter.

    A multiplier that never reaches the model input builds a control file that
    runs, reports phi and calibrates NOTHING (ledger 102). So this asserts the
    transport input file actually changed, and that concentration moved with it
    -- porosity is only worth estimating because it does.
    """

    demo = build_canonical_transport_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )
    cal = demo.model.pest("poros", start_datetime="2024-01-01")
    # `physical=` deliberately sits ABOVE any value this model reaches: it
    # clamps the FINAL value, so a bound at the model's own starting porosity
    # would swallow the multiplier and make this test pass or fail for the
    # wrong reason.
    spec = cal.parameterize(
        "porosity", style="constant", bounds=(0.5, 2.0), physical=(0.01, 0.9)
    )
    cal.observe(demo.conc_targets)
    pst = cal.build("poros.pst", noptmax=0)

    transport_name = demo.transport.transport_name
    assert all(name.startswith(transport_name) for name in spec.resolved_files), (
        f"porosity resolved to {spec.resolved_files}, which is not the GWT "
        "model's files -- the flow model's name was used"
    )
    assert pst.npar_adj == 1

    template = cal.template_workspace
    porosity_file = template / f"{transport_name}.mst_porosity_layer1.txt"

    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template, check=True, capture_output=True, text=True,
    )
    baseline_porosity = np.loadtxt(porosity_file).copy()
    baseline_conc = pd.read_csv(template / "conc_simulated_conc.csv")

    parameters = pst.parameter_data
    parameters.loc[parameters.index[parameters["pargp"] == "porosity"], "parval1"] = 1.6
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=False):
        write_to_template(parameters["parval1"], str(template / tpl), str(template / inp))
    subprocess.run(
        [sys.executable, "forward_run.py"],
        cwd=template, check=True, capture_output=True, text=True,
    )

    assert np.allclose(np.loadtxt(porosity_file), baseline_porosity * 1.6), (
        "the porosity multiplier never reached the transport model input"
    )
    moved = pd.read_csv(template / "conc_simulated_conc.csv")
    difference = (
        moved.drop(columns=["per"]) - baseline_conc.drop(columns=["per"])
    ).abs().to_numpy().max()
    assert difference > 1e-3, (
        "concentration did not respond to porosity, so the parameter is "
        f"unidentifiable from this data (max |dC| = {difference})"
    )


@pytest.mark.slow
def test_a_scalar_griddata_array_is_re_stored_one_file_per_layer(tmp_path):
    """MST porosity is supplied as a scalar, so MF6 externalizes ONE file of
    ``nlay * ncpl`` values. pyEMU looks up ``get_xy([i, j])`` for every row of an
    array against a spatial reference holding ``ncpl`` entries, so that file
    raises ``IndexError: index 441 is out of bounds`` from inside
    ``add_parameters`` -- naming neither the target nor the cause. Re-storing it
    layered gives it the shape ``npf_k`` already has, with the same numbers.
    """

    demo = build_canonical_transport_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )
    transport = demo.model.sim.get_model(demo.transport.transport_name)
    before = np.asarray(transport.get_package("mst").porosity.get_data(), dtype=float)

    cal = demo.model.pest("relayer", start_datetime="2024-01-01")
    spec = cal.parameterize("porosity", style="constant", bounds=(0.5, 2.0))
    cal.observe(demo.conc_targets)
    cal.build("relayer.pst", noptmax=0)

    assert len(spec.resolved_files) == int(transport.modelgrid.nlay)
    after = np.loadtxt(cal.template_workspace / spec.resolved_files[0])
    assert np.allclose(after, before[0]), "relayering changed the values"


@pytest.mark.slow
def test_capture_on_a_whole_grid_array_records_every_layer(tmp_path):
    """``capture=True`` names each captured file's layer so ``plot_field`` can
    ask for one. Recording a whole-grid array as ``{0: prefix}`` would label
    every layer's values "layer 0" and hand plot_field nlay*ncpl values for ncpl
    cells; after relayering there is a real file per layer."""

    demo = build_canonical_transport_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )
    cal = demo.model.pest("captured", start_datetime="2024-01-01")
    cal.parameterize("porosity", style="constant", bounds=(0.5, 2.0), capture=True)
    cal.observe(demo.conc_targets)
    cal.build("captured.pst", noptmax=0)

    metadata = json.loads(
        (cal.template_workspace / METADATA_FILENAME).read_text()
    )
    captured = next(
        entry for entry in metadata["capture_fields"] if entry["name"] == "porosity"
    )
    nlay = int(demo.model.gwf.modelgrid.nlay)
    assert sorted(int(k) for k in captured["layer_prefixes"]) == list(range(nlay))


@pytest.mark.slow
def test_the_transport_demo_carries_heads_as_well_as_concentration(tmp_path):
    """Porosity is absent from the flow equation, and transport velocity is
    ``v = Ki/n``, so K and porosity are near-collinear from concentration alone.
    Heads respond to K and not at all to porosity, which is what makes the pair
    identifiable. A demo that shipped concentration only would teach an
    ill-posed calibration."""

    demo = build_canonical_transport_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )

    assert demo.start_porosity_factor != 1.0, "porosity is not perturbed"
    assert not demo.head_targets.to_long().empty
    head_names = set(demo.head_targets.match_to_model(demo.model)["cell"])
    conc_names = set(demo.conc_targets.match_to_model(demo.model)["cell"])
    assert head_names == conc_names, "heads and concentration sample the same wells"

    transport = demo.model.sim.get_model(demo.transport.transport_name)
    porosity = np.asarray(transport.get_package("mst").porosity.get_data(), dtype=float)
    assert np.allclose(porosity, 0.25 * demo.start_porosity_factor), (
        "the model handed to PEST still has the truth porosity, so calibrating "
        "it would start at the answer"
    )
