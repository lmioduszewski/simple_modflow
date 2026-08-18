from __future__ import annotations

import json
import pickle
import shutil
from pathlib import Path

import numpy as np
import pytest

from examples.mf6.visualization_prt_master_support import (
    MasterExampleConfig,
    build_transient_model,
)
from myflopy.modflow.mf6.interactive_plotting import ModelMapStyle

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.canonical


def test_master_example_defaults_to_large_complex_model():
    config = MasterExampleConfig()

    assert config.ncpl >= 10_000
    assert config.nlay >= 4
    assert config.nper >= 6
    assert config.ncpl * config.nlay >= 40_000


def test_canonical_fixture_enforces_complete_model_contract(canonical_model):
    canonical_model.canonical_contract.validate(canonical_model)
    attached = {name.lower() for name in canonical_model.package_names}
    assert attached >= canonical_model.canonical_contract.required_packages
    icelltype = np.asarray(canonical_model.gwf.npf.icelltype.get_data())
    assert tuple(int(icelltype[layer].reshape(-1)[0]) for layer in range(4)) == (1, 1, 0, 0)
    conductivity = np.asarray(canonical_model.gwf.npf.k.get_data())
    assert np.median(conductivity[2]) < np.median(conductivity[1]) * 0.01
    assert canonical_model.vor.gdf_vorPolys.area.std() > 0
    areas = canonical_model.vor.gdf_vorPolys.area
    assert areas.loc[canonical_model.get_region_cells("all_lakes")].median() < areas.median() * 0.90
    assert areas.loc[canonical_model.get_region_cells("all_streams")].median() < areas.median() * 0.90


def test_canonical_evt_and_uzf_footprints_stay_disjoint(canonical_model):
    """EVT (groundwater ET) must never overlap UZF (unsaturated-zone ET).

    The canonical UZF runs with ``simulate_et`` auto-enabled (it is passed
    ``pet``/``extdp``) but WITHOUT ``linear_gwet``/``square_gwet``, so it removes
    ET from the vadose zone only, while EVT removes it from groundwater. Applying
    both to the same cells would double-count a single PET demand. UZF covers the
    valley floor, EVT the walls -- if a future change grows either footprint into
    the other, this fails instead of silently double-counting.
    """

    evt_cells = set(canonical_model.get_region_cells("upland_et"))
    uzf_cells = set(canonical_model.get_region_cells("uzf_active"))

    assert evt_cells, "canonical model exposes no EVT region"
    assert uzf_cells, "canonical model exposes no UZF region"
    assert not (evt_cells & uzf_cells), (
        "EVT and UZF share cells -> ET is double-counted: "
        f"{sorted(evt_cells & uzf_cells)[:10]}"
    )
    # the un-routed outlet river is likewise carved out of the UZF footprint,
    # exactly as the lake and stream cells are
    riv_cells = set(canonical_model.get_region_cells("outlet_river"))
    assert riv_cells and not (riv_cells & uzf_cells)


def test_canonical_riv_and_evt_carry_real_flux(canonical_run):
    """riv/evt are contract packages, so they must do physical work, not just exist."""

    riv = canonical_run.packages.riv.results.q.get()
    assert not riv.empty
    # riv/evt are aquifer-referenced .cbc records, so the exchange column names
    # its frame: q_gwf (negative = out of the aquifer). An outlet river in
    # equilibrium with the water table both gains and loses, which also keeps the
    # signed-q RdBu colorscale exercised in both directions.
    assert riv["q_gwf"].min() < 0.0 < riv["q_gwf"].max()

    evt = canonical_run.packages.evt.results.q.get()
    assert not evt.empty
    assert (evt["q_gwf"] <= 0.0).all(), "ET may only remove water"
    period0 = evt[evt["per"] == evt["per"].min()]
    active = int((period0["q_gwf"].abs() > 1e-9).sum())
    # a strict subset transpires: cells shallower than the extinction depth draw
    # water, deeper ones yield exactly zero -- that partition IS the feature
    assert 0 < active < len(period0)

    expected_fields = {
        "riv": {"stage", "cond", "rbot"},
        "evt": {"surface", "rate", "depth"},
    }
    for package, fields in expected_fields.items():
        columns = set(getattr(canonical_run.packages, package).inputs.get().columns)
        assert fields <= columns, f"{package} inputs missing {sorted(fields - columns)}"


def test_master_example_full_profile_boundary_heads_stay_above_cell_bottoms():
    workspace = ROOT / ".pytest-work" / "master_example_full_boundary_check"
    shutil.rmtree(workspace, ignore_errors=True)
    try:
        model = build_transient_model(workspace, config=MasterExampleConfig(), name="master_bounds")
        bottoms = model.gwf.disv.botm.get_data()
        for records in model.gwf.chd.stress_period_data.get_data().values():
            for record in records:
                layer, cell = record["cellid"]
                assert float(record["head"]) > float(bottoms[layer, cell])
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_master_example_validation_profile_runs_complex_package_topology():
    workspace = ROOT / ".pytest-work" / "master_example_validation"
    shutil.rmtree(workspace, ignore_errors=True)
    try:
        config = MasterExampleConfig.validation()
        model = build_transient_model(workspace, config=config, name="master_valid")

        required = {
            "disv", "npf", "sto", "chd", "ghb", "drn", "rch", "wel",
            "lak", "sfr", "mvr", "uzf", "oc",
        }
        assert required <= {name.lower() for name in model.package_names}
        assert model.gwf.modelgrid.nlay == config.nlay
        assert model.gwf.modelgrid.ncpl == config.ncpl
        assert model.get_region_cells("all_lakes")
        assert model.get_region_cells("all_streams")
        assert model.get_region_cells("uzf_active")

        success, report = model.run_simulation()
        assert success, "\n".join(report[-20:])
        with (model.model_output_folder_path / f"{model.name}.model").open("rb") as file:
            reopened = pickle.load(file)
        assert set(reopened.targets.keys()) == {"drn_flow", "heads", "lake_stage", "sfr_flow", "sfr_stage"}
        heads = model.gwf.output.head()
        assert len(heads.get_kstpkper()) == config.nper * config.steps_per_period
        assert heads.get_data().shape == (config.nlay, 1, config.ncpl)
        assert (model.model_output_folder_path / "master_heads.csv").exists()
        assert (model.model_output_folder_path / "master_sfr_stage.csv").exists()
        assert (model.model_output_folder_path / "master_sfr_flow.csv").exists()
        assert (model.model_output_folder_path / "master_drn.csv").exists()

        simulated = {
            "heads": model.targets.heads.targets.simulated_heads(model),
            "lake_stage": model.targets.lake_stage.targets.simulated_series(model),
            "sfr_stage": model.targets.sfr_stage.targets.simulated_series(model),
            "sfr_flow": model.targets.sfr_flow.targets.simulated_series(model),
            "drn_flow": model.targets.drn_flow.targets.simulated_series(model),
        }
        assert all(not frame.empty for frame in simulated.values())
        assert not model.packages.lak.results.stage.get().empty
        assert not model.packages.sfr.results.q.get().empty
        assert not model.outputs.uzf.ifno_to_cellid.empty
        from myflopy.modflow.mf6.canonical import canonical_head_signals

        head_signals = canonical_head_signals(model)
        assert min(head_signals["spatial_range_by_layer"]) > 10.0
        assert max(head_signals["temporal_range_by_layer"]) > 3.0
        assert head_signals["maximum_drawdown_by_layer"][3] > 1.0
        from myflopy.modflow.mf6.canonical import canonical_feature_signals

        feature_signals = canonical_feature_signals(model)
        assert feature_signals["head_change"]["pond_mound"] > 0.5
        assert feature_signals["head_change"]["shallow_pumping"] > 1.0
        assert feature_signals["head_change"]["deep_pumping"] > 1.0
        assert feature_signals["seepage_peak"]["north_springs"] > 0.0
        assert feature_signals["seepage_peak"]["south_springs"] > 0.0
        from myflopy.modflow.mf6.canonical import canonical_sfr_signals

        sfr_signals = canonical_sfr_signals(model)
        assert sfr_signals["reach_count"] >= 40
        assert sfr_signals["wet_reach_count"] == sfr_signals["reach_count"]
        assert sfr_signals["minimum_depth"] > 0.02
        assert sfr_signals["gaining_reach_count"] >= 3
        assert sfr_signals["losing_reach_count"] >= 3
        assert sfr_signals["minimum_routed_flow"] > 5_000.0

        slider = model.visualize.head_map_slider_html(
            workspace / "validation_head_slider.html",
            kstpkpers=heads.get_kstpkper()[:2],
            layer=0,
            style=ModelMapStyle(show_contours=False, dpi=40, figsize=(4, 3)),
            embed_frames=False,
        )
        assert slider.frame_count == 2
        assert len(list(slider.frame_directory.glob("*.png"))) == 2
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_master_notebook_documents_full_integration_surface():
    path = ROOT / "examples" / "mf6" / "notebooks" / "master_large_model_visualization_prt_workflow.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    source = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
    )

    for required in (
        "10,000-cell",
        "head_map_slider_html",
        "head_layer_mosaic_slider_html",
        "cross_section_slider_html",
        "plotly_head_map_animation",
        "plotly_cross_section_animation",
        "PRTReleasePoints",
        "open_prt_run",
        "export_3d_html",
    ):
        assert required in source


def test_canonical_master_notebook_set_uses_one_builder():
    notebook_root = ROOT / "examples" / "mf6" / "notebooks"
    # The numbered model run set (canonical_00..canonical_03) builds the shared
    # mf.build_canonical_model() model directly. The PEST notebooks
    # (canonical_04/05/06) calibrate the same model via
    # build_canonical_calibration_demo, and canonical_model_template.ipynb is an
    # editable scaffold -- all are intentionally excluded here.
    paths = sorted(
        path
        for path in notebook_root.glob("canonical_0[0-3]*.ipynb")
        if not path.name.endswith(".executed.ipynb")
    )
    assert len(paths) == 4
    combined = ""
    for path in paths:
        notebook = json.loads(path.read_text(encoding="utf-8"))
        source = "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])
        assert "build_canonical_model" in source
        assert "notebook_header" in source
        assert sum(cell["cell_type"] == "markdown" for cell in notebook["cells"]) >= 2
        assert sum(cell["cell_type"] == "code" for cell in notebook["cells"]) >= 2
        combined += source
    for required in (
        "infiltration",
        "seepage",
        "observations",
        "visual",
        "PRT",
        "parallel",
    ):
        assert required.lower() in combined.lower()


def test_pest_notebooks_calibrate_the_canonical_model():
    notebook_root = ROOT / "examples" / "mf6" / "notebooks"
    pest_paths = sorted(
        path
        for path in notebook_root.glob("canonical_0[4-6]*.ipynb")
        if not path.name.endswith(".executed.ipynb")
    )
    assert len(pest_paths) == 3
    sources = {}
    for path in pest_paths:
        notebook = json.loads(path.read_text(encoding="utf-8"))
        source = "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])
        # Every PEST notebook calibrates the canonical model itself; the old
        # standalone demo models are gone.
        assert "build_canonical_calibration_demo" in source
        assert "gold_standard" not in source
        assert "build_calibration_demo(" not in source
        # PEST projects are constructed via the model.pest(...) front door, not
        # by importing PestProject directly.
        assert "model.pest(" in source
        assert not any(
            output.get("output_type") == "error"
            for cell in notebook["cells"]
            for output in cell.get("outputs", [])
        )
        sources[path.name] = source

    combined = "\n".join(sources.values())
    for required in ("parameterize", "observe", "forecast", "settings", "run_ies"):
        assert required in combined

    # The uncertainty notebook exercises the prior Monte Carlo + IES diagnostics
    # and the spatially-varying-K "property pattern" maps.
    ies_source = sources["canonical_06_pest_ies_uncertainty.ipynb"]
    for required in (
        "cal.prior",
        "plot_prior_vs_obs",
        "conflict",
        "plot_phi_distribution",
        "plot_phi_contributions",
        "parameters_at_bounds",
        "style='pilotpoints'",
        "capture=True",
        "plot_field",
        "plot.map(",
        # 6.4B: the "did the data inform this?" and "where is it biased?" figures
        "stat='reduction'",
        "plot_field_mosaic",
        "plot_obs_residuals",
        "obs_residuals()",
    ):
        assert required in ies_source


def test_transport_notebook_teaches_the_transport_calibration_surface():
    """``canonical_07`` is checked here because it falls outside both other globs.

    ``test_canonical_master_notebook_set_uses_one_builder`` globs
    ``canonical_0[0-3]`` and asserts exactly 4; the PEST test globs
    ``canonical_0[4-6]`` and asserts exactly 3. A ``canonical_07`` matches neither,
    so without this it would ship with NO static checking at all -- and nothing in
    this repo executes notebooks, so static checking is all there is.
    """

    import json
    from pathlib import Path

    notebooks = sorted(
        (Path(__file__).resolve().parents[1] / "examples/mf6/notebooks").glob(
            "canonical_07*.ipynb"
        )
    )
    assert len(notebooks) == 1, f"expected one transport notebook, found {notebooks}"
    notebook = json.loads(notebooks[0].read_text(encoding="utf-8"))
    source = "\n".join(
        "".join(cell["source"]) for cell in notebook["cells"]
    )

    # the house scaffolding every canonical notebook shares
    assert "notebook_header" in source
    assert sum(cell["cell_type"] == "markdown" for cell in notebook["cells"]) >= 2
    assert sum(cell["cell_type"] == "code" for cell in notebook["cells"]) >= 2

    # the surface this notebook exists to teach -- each name is a feature that
    # would otherwise have no worked example anywhere
    for required in (
        "build_canonical_transport_calibration_demo",  # the transport fixture
        "transport_view",                              # reading the GWT sibling
        "view.conc.map",                               # the concentration reader
        "view.budget",                                 # model.budget.<term>
        "particle_tracking.prt",                       # PRT on the coupled model
        "pathlines.map",
        "travel_time.map",
        "cal.observe(demo.conc_targets)",              # calibration ON CONCENTRATION
        "cal.parameterize('k'",
        "forward_run.py",
    ):
        assert required in source, f"canonical_07 no longer teaches {required!r}"

    # The heavy cell is GATED -- assert the flag exists, not its current value.
    # Asserting `RUN_IES = False` would fail whenever someone flips it to True to
    # actually run the calibration, which is the whole point of the gate.
    assert "RUN_IES" in source

    # No ERROR outputs -- matching the rule the other notebook tests use. Asserting
    # "no outputs at all" would fail the moment someone RUNS the notebook locally,
    # which is the normal way to work on one; the pre-commit nbstripout hook is what
    # keeps outputs out of the committed file.
    assert not any(
        output.get("output_type") == "error"
        for cell in notebook["cells"]
        for output in cell.get("outputs", [])
    ), "canonical_07 has a cell that raised; fix it before committing"
