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


@pytest.mark.slow
def test_pilot_points_on_porosity_interpolate_the_porosity_array(tmp_path):
    """This target was REFUSED for one day (ledger 110) on the grounds that the
    pilot-point interpolation was hardwired to the flow model's NPF K.

    That was a true description of a bug, not a real limitation: the same
    hardwiring silently overwrote K33 with horizontal K for a FLOW target too
    (ledger 115). With the base array resolved from the recipe, transport pilot
    points are simply correct, so this asserts the capability rather than the
    refusal -- specifically that unit multipliers reproduce the MST porosity
    array and not the K field, which is what a regression would break.
    """

    demo = build_canonical_transport_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing()
    )
    cal = demo.model.pest("pp_poros", start_datetime="2024-01-01")
    cal.parameterize("porosity", style="pilotpoints", pp_space=6,
                     bounds=(0.5, 2.0), physical=(0.01, 0.9))
    cal.observe(demo.conc_targets)
    cal.build("pp_poros.pst", noptmax=0)

    result = subprocess.run(
        [sys.executable, "forward_run.py"], cwd=cal.template_workspace,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]

    transport = demo.model.sim.get_model(demo.transport.transport_name)
    truth = np.asarray(transport.get_package("mst").porosity.get_data(), dtype=float)
    written = np.loadtxt(
        cal.template_workspace
        / f"{demo.transport.transport_name}.mst_porosity_layer1.txt"
    )
    assert np.allclose(written, truth[0], rtol=1e-6), (
        "pilot points did not interpolate the porosity array"
    )
    k = np.asarray(demo.model.gwf.npf.k.get_data(), dtype=float)
    assert not np.allclose(written, k[0], rtol=1e-3), (
        "porosity was overwritten with the hydraulic conductivity field"
    )


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


def _relayer_stub(*, layered, nlay, values):
    """A project whose one array reports the given storage state."""

    from types import SimpleNamespace

    calls = []

    array = SimpleNamespace(
        supports_layered=lambda: True,
        _get_storage_obj=lambda: SimpleNamespace(layered=layered),
        get_data=lambda: values,
        store_internal=lambda: calls.append("store_internal"),
        make_layered=lambda: calls.append("make_layered"),
        set_data=lambda data: calls.append(("set_data", len(data))),
    )
    transport = SimpleNamespace(
        get_package=lambda name: SimpleNamespace(porosity=array),
        modelgrid=SimpleNamespace(nlay=nlay),
        model_type="gwt6",
    )
    simulation = SimpleNamespace(
        model_names=["flow", "flow_t"],
        get_model=lambda name: transport if name == "flow_t" else SimpleNamespace(
            model_type="gwf6"
        ),
    )
    return SimpleNamespace(model=SimpleNamespace(sim=simulation)), calls


def test_relayering_is_skipped_when_it_would_be_wrong_or_pointless():
    """Found by the two-calibration test, not by inspection: on the SECOND
    build the arrays are already external, and FloPy refuses to make external
    data layered (``Converting external file data into layered data currently
    not support``). Relayering unconditionally therefore broke every
    second-build-on-one-model, which is exactly the capability ledger 107 was
    restoring.

    Also skipped for a single-layer model: a whole-grid file there already holds
    exactly ``ncpl`` values, so pyEMU indexes it fine.
    """

    from myflopy.modflow.mf6.pest.native_parameters import relayer_array_target

    recipe = resolve_target("porosity")
    values = np.full((4, 9), 0.25)

    project, calls = _relayer_stub(layered=True, nlay=4, values=values)
    assert relayer_array_target(project, recipe) is False
    assert calls == [], "an already-layered array must be left alone"

    project, calls = _relayer_stub(layered=False, nlay=1, values=np.full((1, 9), 0.25))
    assert relayer_array_target(project, recipe) is False
    assert calls == []

    project, calls = _relayer_stub(layered=False, nlay=4, values=values)
    assert relayer_array_target(project, recipe) is True
    # store_internal FIRST: FloPy refuses make_layered on external data.
    assert calls == ["store_internal", "make_layered", ("set_data", 4)]


def test_a_stale_whole_grid_file_never_wins_over_the_per_layer_files(tmp_path):
    """Relayering an array that was ALREADY external leaves the old whole-grid
    file on disk -- unreferenced by the package, but still there. Resolving the
    exact name first would point every parameter at values MODFLOW no longer
    reads, and the calibration would run happily against them."""

    from myflopy.modflow.mf6.pest.native_parameters import _resolve_files

    recipe = resolve_target("porosity")
    (tmp_path / "m_t.mst_porosity.txt").write_text("0.25\n", encoding="utf-8")
    assert _resolve_files(tmp_path, "m_t", recipe) == ["m_t.mst_porosity.txt"]

    for layer in (1, 2):
        (tmp_path / f"m_t.mst_porosity_layer{layer}.txt").write_text("0.25\n", encoding="utf-8")
    assert _resolve_files(tmp_path, "m_t", recipe) == [
        "m_t.mst_porosity_layer1.txt", "m_t.mst_porosity_layer2.txt"
    ]


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
